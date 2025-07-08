# generate_final_binding_regions.py 
"""
Final peak-processing script for a calling card pipeline.

This script takes pre-called peaks for each group, resizes them, and generates
per-group intersection matrices. It then creates a final consensus peak set and
uses it to generate a final consensus count matrix with insertion counts from all groups.
"""

import sys
import argparse
import logging
import os
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from typing import Optional, List, Tuple

import polars as pl
import pandas as pd
import pybedtools
from tqdm import tqdm

VERSION = "3.3.1"

def setup_logging(output_dir: Path) -> None:
    """Initializes logging to a file and stdout."""
    log_file = output_dir / "consensus_peaks.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] - %(message)s",
        handlers=[logging.FileHandler(log_file, mode="w"), logging.StreamHandler(sys.stdout)],
    )
    logging.info(f"Logging initialized. Log file: {log_file}")


def process_group_intersections(
    group_name: str,
    peak_file_path: Path,
    all_group_bedgraphs: List[Tuple[str, Path]]
) -> Optional[pl.DataFrame]:
    """
    Worker: (1) Loads and resizes peaks for one group, then in a loop:
    (2) Intersects those peaks with each group's raw bedgraph insertions
        and sums the insertion counts (column 4) for each peak.
    """
    try:
        peaks_pl = pl.read_csv(
            peak_file_path, separator="\t", has_header=False,
            new_columns=["chrom", "start", "end"]
        ).select(["chrom", "start", "end"])

        if peaks_pl.is_empty():
            logging.warning(f"No peaks in file for '{group_name}'. Skipping.")
            return None

        peaks_pl = peaks_pl.with_columns(
            (pl.col("end") - pl.col("start")).alias("width")
        ).with_columns(
            (pl.col("start") + (pl.col("width") / 2)).cast(pl.Int64).alias("midpoint")
        )

        peaks_pl = peaks_pl.with_columns(
            pl.when(pl.col("width") < 200)
              .then(pl.col("midpoint") - 100)
              .otherwise(pl.col("start"))
              .alias("start"),
            pl.when(pl.col("width") < 200)
              .then(pl.col("midpoint") + 100)
              .otherwise(pl.col("end"))
              .alias("end")
        ).drop("width", "midpoint")

        peaks_pl = peaks_pl.with_columns(
            pl.when(pl.col("start") < 0).then(0).otherwise(pl.col("start")).alias("start")
        )

        peaks_pl = peaks_pl.with_columns([
            pl.lit(group_name).alias("group_name"),
            (pl.col("chrom") + ":" + pl.col("start").cast(str) + "-" + pl.col("end").cast(str)).alias("peak_id")
        ])

        peaks_bt = pybedtools.BedTool.from_dataframe(
            peaks_pl.select(["chrom", "start", "end", "peak_id"]).to_pandas()
        ).sort()

        for other_name, other_bed in all_group_bedgraphs:
            col_name = f"{other_name}_insertion_count"
            other_insertions_bt = pybedtools.BedTool(str(other_bed))
            intersect_result = peaks_bt.intersect(other_insertions_bt, wa=True, wb=True)

            if intersect_result.count() > 0:
                intersect_df = intersect_result.to_dataframe(
                    names=['chrom_peak', 'start_peak', 'end_peak', 'peak_id',
                           'chrom_ins', 'start_ins', 'end_ins', 'insertion_value']
                )
                summed_counts = intersect_df.groupby('peak_id')['insertion_value'].sum().reset_index()
                summed_counts.rename(columns={'insertion_value': col_name}, inplace=True)
                peaks_pl = peaks_pl.join(pl.from_pandas(summed_counts), on='peak_id', how='left')
            else:
                peaks_pl = peaks_pl.with_columns(pl.lit(0).cast(pl.Int64).alias(col_name))

        all_count_cols = [f"{g}_insertion_count" for g, _ in all_group_bedgraphs]
        peaks_pl = peaks_pl.with_columns(
            [pl.col(c).fill_null(0).cast(pl.Int64) for c in all_count_cols if c in peaks_pl.columns]
        )
        peaks_pl = peaks_pl.with_columns((pl.col("end") - pl.col("start")).alias("peak_width"))

        logging.info(f"Processed intersections for '{group_name}' using {peaks_pl.height} resized peaks.")
        return peaks_pl

    except Exception as e:
        logging.error(f"Error in worker for '{group_name}': {e}", exc_info=True)
        return None

def count_consensus_worker(
    group_name: str, bedgraph_path: Path, consensus_peaks_bt: pybedtools.BedTool
) -> Tuple[str, pl.Series]:
    """Worker to count insertions from one group in the final consensus peaks."""
    col_name = f"{group_name}_insertion_count"
    bedgraph_bt = pybedtools.BedTool(str(bedgraph_path))
    
    mapped = consensus_peaks_bt.map(b=bedgraph_bt, c=4, o='sum').to_dataframe(
        names=["chrom", "start", "end", "peak_id", "count"]
    )
    counts = pl.from_pandas(mapped)["count"].str.replace(".", "0").cast(pl.Int64)
    counts = counts.rename(col_name)
    return group_name, counts


def run_worker(task_args):
    """Helper function to unpack arguments for the first ProcessPoolExecutor."""
    return process_group_intersections(*task_args)

def run_consensus_worker(task_args):
    """Helper function to unpack arguments for the consensus counting worker."""
    return count_consensus_worker(*task_args)


def main():
    """Main function to orchestrate the peak counting and consensus workflow."""
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--peak_dir", required=True, help="Directory containing pre-called SPAN peak files.")
    p.add_argument("--peak_suffix", required=True, help="Suffix of the peak files to use (e.g., '_b50.span.peak').")
    p.add_argument("--bedgraph_dir", required=True, help="Directory with per-group insertion count bedGraph files.")
    p.add_argument("--output_dir", required=True, help="Directory to save output files.")
    p.add_argument("--annotation_file", required=True, help="Annotation file mapping samples to groups.")
    p.add_argument("--workers", type=int, default=12, help="Number of parallel worker processes.")
    p.add_argument('--max_dist_merge', type=int, default=100, help="If > 0, merge final consensus peaks within this distance.")
    args = p.parse_args()

    # --- Setup ---
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    setup_logging(output_dir)
    pybedtools_tmp_dir = output_dir / "pybedtools_tmp"
    pybedtools_tmp_dir.mkdir(exist_ok=True)
    pybedtools.helpers.set_tempdir(str(pybedtools_tmp_dir))

    # --- Prepare Inputs ---
    ann = pl.read_csv(args.annotation_file, separator="\t" if args.annotation_file.endswith(".tsv") else ",")
    group_names = ann["group_name"].unique().to_list()
    
    all_bedgraphs = [(g, Path(args.bedgraph_dir) / f"{g}.bedgraph") for g in group_names if (Path(args.bedgraph_dir) / f"{g}.bedgraph").is_file()]

    tasks = []
    for group_name, _ in all_bedgraphs:
        peak_file = Path(args.peak_dir) / f"{group_name}{args.peak_suffix}"
        if peak_file.is_file():
            tasks.append((group_name, peak_file, all_bedgraphs))
        else:
            logging.warning(f"Could not find peak file for group '{group_name}' at: {peak_file}. Skipping group.")

    # --- Execute Per-Group Processing in Parallel ---
    logging.info(f"Processing intersections for {len(tasks)} groups using {args.workers} workers...")
    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as exe:
        future_results = exe.map(run_worker, tasks)
        for res in tqdm(future_results, total=len(tasks)):
            if res is not None:
                results.append(res)

    if not results:
        logging.error("No peak intersections could be processed. Exiting.")
        sys.exit(1)

    # --- Consolidate and Save Per-Group Results ---
    per_group = pl.concat(results)
    per_group_dir = output_dir / "per_group_peak_matrices"
    per_group_dir.mkdir(exist_ok=True)
    for g in per_group["group_name"].unique().to_list():
        df_g = per_group.filter(pl.col("group_name") == g)
        if df_g.height > 0:
            df_g.write_csv(per_group_dir / f"{g}_peak_matrix.tsv", separator="\t")

    # --- Generate Final Consensus Peaks ---
    logging.info("Generating final consensus peak set...")
    consensus_bt = pybedtools.BedTool.from_dataframe(per_group.select(["chrom", "start", "end"]).to_pandas()).sort()
    
    if args.max_dist_merge > 0:
        merged = consensus_bt.merge(d=args.max_dist_merge)
    else:
        merged = consensus_bt.merge()
    
    cons_df = pl.from_pandas(merged.to_dataframe(names=["chrom", "start", "end"]))
    cons_df = cons_df.with_columns([
        (pl.col("chrom") + ":" + pl.col("start").cast(str) + "-" + pl.col("end").cast(str)).alias("peak_id"),
        (pl.col("end") - pl.col("start")).alias("peak_width")
    ])

    output_path = output_dir / "consensus_peaks_only.tsv"
    cons_df.write_csv(output_path, separator="\t")
    logging.info(f"✅ Saved {cons_df.height} consensus peaks to: {output_path}")

    # --- Generate Final Consensus Count Matrix ---
    logging.info("Generating final consensus count matrix...")
    final_consensus_bt = pybedtools.BedTool.from_dataframe(
        cons_df.select(["chrom", "start", "end", "peak_id"]).to_pandas()
    ).sort()
    
    consensus_tasks = [(g, bg, final_consensus_bt) for g, bg in all_bedgraphs]
    
    with ProcessPoolExecutor(max_workers=args.workers) as exe:
        # FIX: Use the named helper function instead of a lambda
        future_results = exe.map(run_consensus_worker, consensus_tasks)
        for res in tqdm(future_results, total=len(consensus_tasks)):
            if res is not None:
                group_name, counts_series = res
                cons_df = cons_df.with_columns(counts_series)

    matrix_path = output_dir / "consensus_peak_counts.tsv"
    cons_df.write_csv(matrix_path, separator="\t")
    logging.info(f"✅ Successfully generated consensus count matrix at: {matrix_path}")


if __name__ == '__main__':
    main()
