# generate_final_binding_regions.py
"""
Final peak-calling script for a calling card pipeline.

This script takes the per-group raw insertion count bedGraphs, expands them
to represent each individual insertion, calls peaks for each group, creates
a comprehensive matrix of all peaks vs. all insertions, and finally merges
them to generate a final consensus peak set with aggregated counts.
"""
# --- Standard Library Imports ---
import sys
import argparse
import logging
import os
from pathlib import Path
from contextlib import redirect_stdout, redirect_stderr
from concurrent.futures import ProcessPoolExecutor
from typing import Optional, Tuple, Dict

# --- Third-party Imports ---
import polars as pl
import pandas as pd
import pybedtools
from tqdm import tqdm
import pycallingcards as cc

# --- Script Version ---
VERSION = "1.0.0"

# --- Peak Calling Configuration ---
PEAK_CALLING_PARAMS_DEFAULTS = {
    'method': 'CCcaller',
    'reference': 'hg38',
    'pvalue_cutoff': 0.01,
    'pvalue_adj_cutoff': 0.01,
    'min_insertions': 5,
    'minlen': 0,
    'extend': 0,
    'maxbetween': 250,
    'minnum': 0,
    'lam_win_size': 1000000
}

# --- Function Definitions ---

def setup_logging(output_dir: Path) -> None:
    log_file = output_dir / "consensus_peaks.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] - %(message)s",
        handlers=[logging.FileHandler(log_file, mode="w"), logging.StreamHandler(sys.stdout)],
    )
    logging.info(f"Logging initialized. Log file: {log_file}")


def process_group_bedgraph(task: tuple) -> Optional[Tuple[pl.DataFrame, pl.DataFrame]]:
    """
    Simplified worker function. Reads data for one group and calls peaks.
    All intersections are handled in the main process after parallelization.
    """
    group_name, bedgraph_path, peak_calling_params = task
    try:
        raw_pd = pd.read_csv(
            bedgraph_path, sep="\t", header=None,
            names=["chrom", "start", "end", "unique_count"],
            dtype={"chrom": str, "start": int, "end": int, "unique_count": int},
        )
        if raw_pd.empty:
            logging.warning(f"Bedgraph for group '{group_name}' is empty. Skipping.")
            return None

        expanded_pd = raw_pd.loc[raw_pd.index.repeat(raw_pd["unique_count"])]
        if len(expanded_pd) != raw_pd['unique_count'].sum():
            logging.error(f"✖ Integrity check FAILED for group '{group_name}'.")
            return None
        logging.info(f"Group '{group_name}': Expanded to {len(expanded_pd):,} insertion events. ✔")

        pycc_input_df = expanded_pd[["chrom", "start", "end"]].rename(columns={"chrom": "Chr", "start": "Start", "end": "End"})
        
        with open(os.devnull, 'w') as f, redirect_stdout(f), redirect_stderr(f):
            group_peaks_pd = cc.pp.call_peaks(expdata=pycc_input_df, **peak_calling_params)

        expanded_pl = pl.from_pandas(expanded_pd).with_columns(pl.lit(group_name).alias("group_name"))
        
        if group_peaks_pd is None or group_peaks_pd.empty:
            logging.warning(f"No peaks were called for group '{group_name}'.")
            return (pl.DataFrame(), expanded_pl)

        group_peaks_pl = pl.from_pandas(group_peaks_pd)
        rename_map = {
            "Chr": "chrom", "Start": "start", "End": "end",
            "Experiment Insertions": "group_insertions",
            "Reference Insertions": "reference_insertions",
            "Expected Insertions": "expected_insertions",
            "pvalue": "pval", "pvalue_adj": "pval_adj",
        }
        group_peaks_pl = group_peaks_pl.rename({k: v for k, v in rename_map.items() if k in group_peaks_pl.columns})
        
        group_peaks_pl = group_peaks_pl.with_columns(pl.lit(group_name).alias("group_name"))
        
        # Expand peaks to at least 200 bp if under that
        group_peaks_pl = group_peaks_pl.with_columns([
            (pl.col("end") - pl.col("start")).alias("width")
        ])

        # Compute padded coordinates
        group_peaks_pl = group_peaks_pl.with_columns([
            pl.when(pl.col("width") < 200)
              .then((pl.col("start") - ((200 - pl.col("width") + 1) // 2)).cast(pl.Int64))
              .otherwise(pl.col("start")).alias("start"),
            pl.when(pl.col("width") < 200)
              .then((pl.col("end") + ((200 - pl.col("width")) // 2)).cast(pl.Int64))
              .otherwise(pl.col("end")).alias("end")
        ])


        logging.info(f"Group '{group_name}': Successfully called {len(group_peaks_pl)} peaks.")
        return group_peaks_pl, expanded_pl

    except Exception as e:
        logging.error(f"A failure occurred in group {group_name}: {e}", exc_info=True)
        return None

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("--input_dir", required=True, help="Directory of raw group bedgraphs.")
    p.add_argument("--output_dir", required=True, help="Directory to save outputs.")
    p.add_argument("--annotation_file", required=True, help="CSV/TSV with group_name column.")
    p.add_argument("--workers", type=int, default=12, help="Number of parallel processes.")
    p.add_argument('--merge_with_distance', action='store_true', help='If set, merge consensus peaks within the --maxbetween distance.')
    for k, v in PEAK_CALLING_PARAMS_DEFAULTS.items():
        p.add_argument(f"--{k}", type=type(v), default=v)
    args = p.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    setup_logging(output_dir)
    logging.info(f"Final Consensus Peak Caller v{VERSION}")
    
    logging.info("--- Peak Calling Parameters ---")
    peak_params_for_logging = {k: getattr(args, k) for k in PEAK_CALLING_PARAMS_DEFAULTS}
    for key, value in peak_params_for_logging.items():
        logging.info(f"  {key}: {value}")
    logging.info("-----------------------------")
    
    input_dir = Path(args.input_dir)
    per_group_dir = output_dir / "per_group_peaks_FINAL"; per_group_dir.mkdir(exist_ok=True)
    pybedtools.helpers.set_tempdir(str(output_dir / "pybedtools_tmp"))

    ann = pl.read_csv(args.annotation_file, separator="\t" if args.annotation_file.endswith(".tsv") else ",")
    all_group_names = ann["group_name"].unique().to_list()
    peak_params = {k: getattr(args, k) for k in PEAK_CALLING_PARAMS_DEFAULTS}
    
    tasks = [(g, input_dir / f"{g}.bedgraph", peak_params) for g in all_group_names if (input_dir / f"{g}.bedgraph").is_file()]
    
    # --- Phase 1: Call peaks for each group in parallel ---
    logging.info(f"Starting per-group peak calling with {args.workers} workers…")
    all_group_peaks, all_group_insertions = [], []
    with ProcessPoolExecutor(max_workers=args.workers) as exe:
        results = list(tqdm(exe.map(process_group_bedgraph, tasks), total=len(tasks), desc="Calling Peaks Per Group"))

    for res in results:
        if res:
            peaks_df, ins_df = res
            if peaks_df is not None and not peaks_df.is_empty():
                all_group_peaks.append(peaks_df)
            if ins_df is not None and not ins_df.is_empty():
                all_group_insertions.append(ins_df)

    if not all_group_peaks:
        logging.error("No peaks were generated for any group. Exiting.")
        sys.exit(1)

    master_peaks_df = pl.concat(all_group_peaks)

    # --- Phase 2: Create All-vs-All Insertion Matrix ---
    logging.info("Building all-vs-all insertion count matrix...")
    all_insertions_dict: Dict[str, pd.DataFrame] = {}
    insertion_summary_data = []
    for group_name in all_group_names:
        bedgraph_path = input_dir / f"{group_name}.bedgraph"
        if bedgraph_path.is_file():
            raw_pd = pd.read_csv(bedgraph_path, sep="\t", header=None, names=["chrom", "start", "end", "unique_count"])
            if not raw_pd.empty:
                all_insertions_dict[group_name] = raw_pd.loc[raw_pd.index.repeat(raw_pd["unique_count"])][["chrom", "start", "end"]]
                insertion_summary_data.append({"group_name": group_name, "total_insertions": raw_pd['unique_count'].sum()})

    peaks_bed = pybedtools.BedTool.from_dataframe(master_peaks_df.select(["chrom", "start", "end"]).to_pandas()).sort()
    
    for group_name, insertions_df in all_insertions_dict.items():
        insertions_bed = pybedtools.BedTool.from_dataframe(insertions_df)
        new_col_name = f"{group_name}_insertions"
        intersected = peaks_bed.intersect(insertions_bed, c=True, wa=True).to_dataframe(names=["chrom", "start", "end", new_col_name])
        counts_series = pl.Series(new_col_name, intersected[new_col_name].values)
        master_peaks_df = master_peaks_df.with_columns(counts_series)

    master_peaks_df = master_peaks_df.with_columns([
        (pl.col("chrom") + ":" + pl.col("start").cast(str) + "-" + pl.col("end").cast(str)).alias("peak_id"),
        (pl.col("end") - pl.col("start")).alias("peak_width")
    ])
    
    logging.info(f"Saving final per-group peak files with all-vs-all counts to {per_group_dir}")
    for g in all_group_names:
        if g in master_peaks_df['group_name'].unique():
            dfg = master_peaks_df.filter(pl.col("group_name") == g)
            dfg.write_csv(per_group_dir / f"{g}_peaks_matrix.tsv", separator="\t")

    # --- Phase 3: Consensus and Final Aggregation ---
    logging.info("— Generating Final Consensus Peak Set —")
    tool_to_merge = pybedtools.BedTool.from_dataframe(master_peaks_df.select(["chrom", "start", "end"]).to_pandas()).sort()
    
    if args.merge_with_distance:
        logging.info(f"Merging consensus peaks within {args.maxbetween} bp of each other...")
        merged = tool_to_merge.merge(d=args.maxbetween)
    else:
        logging.info("Merging only overlapping consensus peaks...")
        merged = tool_to_merge.merge()

    consensus_df = pl.from_pandas(
        merged.to_dataframe(names=["chrom", "start", "end"])
    ).with_columns([
        (pl.col("chrom") + ":" + pl.col("start").cast(str) + "-" + pl.col("end").cast(str)).alias("peak_id"),
        (pl.col("end") - pl.col("start")).alias("peak_width")
    ])
    logging.info(f"Created {consensus_df.height:,} consensus intervals.")

    all_ins_df = pl.concat(all_group_insertions)
    insert_cols = ["chrom", "start", "end", "group_name"]
    consensus_cols = ["chrom", "start", "end", "peak_id"]
    ins_bt = pybedtools.BedTool.from_dataframe(all_ins_df.select(insert_cols).to_pandas())
    cons_bt = pybedtools.BedTool.from_dataframe(consensus_df.select(consensus_cols).to_pandas())
    intersection_cols = insert_cols + [f"{c}_c" for c in consensus_cols]
    mapped = pl.from_pandas(ins_bt.intersect(cons_bt, wa=True, wb=True).to_dataframe(names=intersection_cols))

    final_df = mapped.group_by(["peak_id_c", "group_name"]).agg(
        pl.len().alias("num_insertions_in_consensus_peak")
    ).sort(["peak_id_c", "group_name"])

    final_output = consensus_df.join(
        final_df, left_on="peak_id", right_on="peak_id_c"
    ).select([
        "chrom","start","end","peak_id","peak_width",
        "group_name","num_insertions_in_consensus_peak"
    ]).sort(["chrom","start","group_name"])

    final_output_path = output_dir / "final_results_with_insertions_per_group.tsv"
    final_output.write_csv(final_output_path, separator="\t")
    
    consensus_peaks_only_path = output_dir / "consensus_peaks_only.tsv"
    consensus_df.select(["chrom", "start", "end", "peak_id", "peak_width"]).sort(["chrom","start"])\
        .write_csv(consensus_peaks_only_path, separator="\t")

    summary_df = pl.DataFrame(insertion_summary_data).sort("group_name")
    summary_path = output_dir / "total_insertions_per_group.tsv"
    summary_df.write_csv(summary_path, separator="\t")

    logging.info(f"Summary of insertions per group saved to: {summary_path}")
    logging.info(f"Aggregated consensus results saved to: {final_output_path}")
    logging.info(f"Consensus peak definitions saved to: {consensus_peaks_only_path}")
    logging.info(f"Per-group peak matrices saved in: {per_group_dir}")
    
    # Cleanup temporary files and directories
    try:
        shutil.rmtree(tmp_dir)
    except Exception:
        pass
    pybedtools.helpers.cleanup()


if __name__ == "__main__":
    main()
