"""
Final peak-calling script for a calling card pipeline.

This script takes the per-group raw insertion count bedGraphs, expands them
for peak calling only, calls peaks for each group, then uses pybedtools intersect
to find overlapping insertions from each group's bedgraph, and finally sums the
insertion counts (column 4) for each peak.
"""

import sys
import argparse
import logging
import os
import shutil
from pathlib import Path
from contextlib import redirect_stdout, redirect_stderr
from concurrent.futures import ProcessPoolExecutor
from typing import Optional, Dict, List, Tuple

import polars as pl
import pandas as pd
import pybedtools
from tqdm import tqdm
import pycallingcards as cc

VERSION = "1.0.0" 

PEAK_CALLING_PARAMS_DEFAULTS = {
    'method': 'CCcaller',
    'reference': 'hg38',
    'pvalue_cutoff': 0.01,
    'pvalue_adj_cutoff': 0.01,
    'min_insertions': 1,
    'minlen': 0,
    'extend': 0,
    'minnum': 0,
    'maxbetween': 250,
    'lam_win_size': 1000000
}

def setup_logging(output_dir: Path) -> None:
    """Initializes logging to a file and stdout."""
    log_file = output_dir / "consensus_peaks.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] - %(message)s",
        handlers=[logging.FileHandler(log_file, mode="w"), logging.StreamHandler(sys.stdout)],
    )
    logging.info(f"Logging initialized. Log file: {log_file}")


def process_group_peaks(
    group_name: str,
    bedgraph_path: Path,
    peak_calling_params: Dict,
    all_group_bedgraphs: List[Tuple[str, Path]]
) -> Optional[pl.DataFrame]:
    """
    Worker: (1) Calls peaks for one group, then in a loop:
    (2A) Intersects peaks with each group's raw bedgraph insertions.
    (2B) Identifies which insertion rows intersect which peak.
    (2C) Sums the insertion counts (column 4) for each peak.
    """
    try:
        # Load the raw bedgraph data for the primary group
        raw = pd.read_csv(
            bedgraph_path, sep="\t", header=None,
            names=["chrom", "start", "end", "unique_count"],
            dtype={"chrom": str, "start": int, "end": int, "unique_count": int}
        )
        if raw.empty:
            logging.warning(f"No insertions for '{group_name}'. Skipping.")
            return None

        # (1) Call peaks for the group
        # Expand insertions in memory only for pycallingcards, which requires one row per insertion
        expanded = raw.loc[raw.index.repeat(raw.unique_count)][["chrom", "start", "end"]]
        pycc_in = expanded.rename(columns={"chrom": "Chr", "start": "Start", "end": "End"})

        with open(os.devnull, 'w') as f, redirect_stdout(f), redirect_stderr(f):
            peaks_pd = cc.pp.call_peaks(expdata=pycc_in, **peak_calling_params)

        if peaks_pd is None or peaks_pd.empty:
            logging.warning(f"No peaks for '{group_name}'.")
            return None

        # Convert to Polars for efficient manipulation
        peaks_pl = pl.from_pandas(peaks_pd).rename({
            "Chr": "chrom", "Start": "start", "End": "end",
            "Experiment Insertions": "group_insertions",
            "Reference Insertions": "reference_insertions",
            "Expected Insertions": "expected_insertions",
            "pvalue": "pval", "pvalue_adj": "pval_adj"
        })

        # Add metadata columns and a unique peak ID for grouping
        peaks_pl = peaks_pl.with_columns([
            pl.lit(group_name).alias("group_name"),
            (pl.col("end") - pl.col("start")).alias("width"),
            # Create a unique ID based on original coordinates for reliable joins later
            (pl.col("chrom") + ":" + pl.col("start").cast(str) + "-" + pl.col("end").cast(str)).alias("peak_id")
        ])

        # Standardize peak widths to a minimum of 200bp
        peaks_pl = peaks_pl.with_columns([
            pl.when(pl.col("width") < 200)
              .then((pl.col("start") - ((200 - pl.col("width") + 1) // 2)).cast(pl.Int64))
              .otherwise(pl.col("start")).alias("start"),
            pl.when(pl.col("width") < 200)
              .then((pl.col("end") + ((200 - pl.col("width")) // 2)).cast(pl.Int64))
              .otherwise(pl.col("end")).alias("end")
        ])

        # Create a BedTool object from the peaks, including the unique ID
        peaks_bt = pybedtools.BedTool.from_dataframe(
            peaks_pl.select(["chrom", "start", "end", "peak_id"]).to_pandas()
        ).sort()

        # (2) Loop to intersect with all groups' insertions
        for other_name, other_bed in all_group_bedgraphs:
            col_name = f"{other_name}_insertion_count"
            other_insertions_bt = pybedtools.BedTool(str(other_bed))

            # (2A) Use pybedtools intersect
            # `wa=True` writes the original peak. `wb=True` writes the original intersecting insertion.
            intersect_result = peaks_bt.intersect(other_insertions_bt, wa=True, wb=True)

            if intersect_result.count() > 0:
                # Convert to pandas DataFrame for easy aggregation
                intersect_df = intersect_result.to_dataframe(
                    names=['chrom_peak', 'start_peak', 'end_peak', 'peak_id',
                           'chrom_ins', 'start_ins', 'end_ins', 'insertion_value']
                )
                
                # (2B & 2C) Identify intersecting rows and sum column 4
                # Group by the unique peak identifier and sum the insertion counts
                summed_counts = intersect_df.groupby('peak_id')['insertion_value'].sum().reset_index()
                summed_counts.rename(columns={'insertion_value': col_name}, inplace=True)
                
                # Join the summed counts back to the main Polars DataFrame
                peaks_pl = peaks_pl.join(pl.from_pandas(summed_counts), on='peak_id', how='left')

            else:
                # If no intersections were found for this group, add a column of zeros
                peaks_pl = peaks_pl.with_columns(pl.lit(0).cast(pl.Int64).alias(col_name))

        # Ensure all newly added count columns have nulls filled with 0
        all_count_cols = [f"{g}_insertion_count" for g, _ in all_group_bedgraphs]
        peaks_pl = peaks_pl.with_columns(
            [pl.col(c).fill_null(0) for c in all_count_cols if c in peaks_pl.columns]
        )

        # Update peak_id and width to reflect any resizing
        peaks_pl = peaks_pl.with_columns([
            (pl.col("chrom") + ":" + pl.col("start").cast(str) + "-" + pl.col("end").cast(str)).alias("peak_id"),
            (pl.col("end") - pl.col("start")).alias("peak_width")
        ])

        logging.info(f"Processed '{group_name}' with {peaks_pl.height} peaks.")
        return peaks_pl

    except Exception as e:
        logging.error(f"Error in '{group_name}': {e}", exc_info=True)
        return None

# ✨ NEW: Top-level helper function for multiprocessing
def run_process_group_peaks(task_args):
    """Helper function to unpack arguments for ProcessPoolExecutor.map."""
    return process_group_peaks(*task_args)


def main():
    """Main function to orchestrate the peak calling and counting workflow."""
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input_dir", required=True, help="Directory with per-group bedGraph files.")
    p.add_argument("--output_dir", required=True, help="Directory to save output files.")
    p.add_argument("--annotation_file", required=True, help="Annotation file mapping samples to groups.")
    p.add_argument("--workers", type=int, default=12, help="Number of parallel worker processes.")
    p.add_argument('--merge_with_distance', action='store_true', help="If set, merge peaks within the specified distance.")
    
    # Add peak calling parameters from defaults
    for k, v in PEAK_CALLING_PARAMS_DEFAULTS.items():
        if k != "maxbetween":
            p.add_argument(f"--{k}", type=type(v), default=v)
    p.add_argument("--maxbetween", type=int, default=250, help="Maximum distance between peaks to merge.")
    args = p.parse_args()

    # --- Setup ---
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    setup_logging(output_dir)
    pybedtools.helpers.set_tempdir(str(output_dir / "pybedtools_tmp"))

    # --- Prepare Inputs ---
    ann = pl.read_csv(
        args.annotation_file,
        separator="\t" if args.annotation_file.endswith(".tsv") else ","
    )
    group_names = ann["group_name"].unique().to_list()
    all_bedgraphs = [
        (g, Path(args.input_dir) / f"{g}.bedgraph")
        for g in group_names
        if (Path(args.input_dir) / f"{g}.bedgraph").is_file()
    ]

    # --- Create Tasks for Parallel Processing ---
    peak_params = {k: getattr(args, k) for k in PEAK_CALLING_PARAMS_DEFAULTS}
    tasks = [
        (g, bg, peak_params, all_bedgraphs)
        for g, bg in all_bedgraphs
    ]

    # --- Execute in Parallel ---
    logging.info(f"Calling peaks & counting insertions using {args.workers} workers...")
    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as exe:
    
        #  Use the picklable helper function
        future_results = exe.map(run_process_group_peaks, tasks)
        for res in tqdm(future_results, total=len(tasks)):
            if res is not None:
                results.append(res)

    if not results:
        logging.error("No peaks were generated across all groups. Exiting.")
        sys.exit(1)

    # --- Consolidate and Save Results ---
    per_group = pl.concat(results)
    per_group_dir = output_dir / "per_group_peaks_FINAL"
    per_group_dir.mkdir(exist_ok=True)
    for g in group_names:
        df_g = per_group.filter(pl.col("group_name") == g)
        if df_g.height > 0:
            df_g.write_csv(per_group_dir / f"{g}_peaks_matrix.tsv", separator="\t")

    # --- Generate Consensus Peaks ---
    peaks_bt = pybedtools.BedTool.from_dataframe(
        per_group.select(["chrom", "start", "end"]).to_pandas()
    ).sort()
    
    # Merge with  max between allowable distance if merge_with_distance is used.
    if args.merge_with_distance:
        merged = peaks_bt.merge(d=args.maxbetween)
    else:
        merged = peaks_bt.merge()
        
    cons_df = pl.from_pandas(
        merged.to_dataframe(names=["chrom", "start", "end"])
    ).with_columns([
        (pl.col("chrom") + ":" + pl.col("start").cast(str) + "-" + pl.col("end").cast(str)).alias("peak_id"),
        (pl.col("end") - pl.col("start")).alias("peak_width")
    ])

    cons_df.write_csv(output_dir / "consensus_peaks_only.tsv", separator="\t")
    logging.info(f"✅ Saved consensus peaks: {cons_df.height} intervals.")


if __name__ == '__main__':
    main()
