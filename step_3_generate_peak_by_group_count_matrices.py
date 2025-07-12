# step_3_generate_peak_by_group_count_matrices.py
"""
This script takes the peaks from the SPAN peak caller and the raw unique insertion count
bedgraph files for each group to generate per‑group intersection matrices. These matrices
have peak_id as rows and the insertion count of each group as columns.

--merge_distance (default 200) argument lets you merge any resized peaks that are ≤ N bp apart.
   Internally this is done with pybedtools' BedTool.merge(d=N). A value of 0 (default)
   disables merging.

"""

# ==============================================================================
# LIBRARY IMPORTS
# ==============================================================================
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

# ==============================================================================
# SCRIPT VERSION
# ==============================================================================
VERSION = "1.0.0"

# ==============================================================================
# HELPER & UTILITY FUNCTIONS
# ==============================================================================

def setup_logging(output_dir: Path) -> None:
    """Initializes logging to a file and stdout."""
    log_file = output_dir / "consensus_peaks.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] - %(message)s",
        handlers=[logging.FileHandler(log_file, mode="w"), logging.StreamHandler(sys.stdout)],
    )
    logging.info(f"Logging initialized. Log file: {log_file}")


# ---------------------------------------------------------------------------
# Core worker function -------------------------------------------------------
# ---------------------------------------------------------------------------

def process_group_intersections(
    group_name: str,
    peak_file_path: Path,
    all_group_bedgraphs: List[Tuple[str, Path]],
    merge_distance: int,
) -> Optional[pl.DataFrame]:
    """Resize, optionally merge, and intersect peaks for one group.

    1.  Peaks <200 bp are sym‑metrically padded to 200 bp around the midpoint.
    2.  If *merge_distance* > 0, adjacent/overlapping peaks within that distance are merged
        via ``bedtools merge -d N`` (pybedtools wrapper).
    3.  The resized/merged peak set is intersected against every group's insertion‑map
        bedgraph to build an all‑vs‑all matrix (rows = this group's peaks;
        columns = insertion counts per group).
    """
    try:
        # ------------------------------------------------------------------
        # 1. Load & resize peaks -------------------------------------------
        # ------------------------------------------------------------------
        peaks_pl = pl.read_csv(
            peak_file_path, separator="\t", has_header=False,
            new_columns=["chrom", "start", "end"],
        ).select(["chrom", "start", "end"])

        if peaks_pl.is_empty():
            logging.warning(f"No peaks in file for '{group_name}'. Skipping.")
            return None

        # midpoint & width
        peaks_pl = peaks_pl.with_columns(
            (pl.col("end") - pl.col("start")).alias("width"),
        ).with_columns(
            (pl.col("start") + (pl.col("width") / 2)).cast(pl.Int64).alias("midpoint")
        )
        # pad to 200 bp
        peaks_pl = peaks_pl.with_columns(
            pl.when(pl.col("width") < 200)
              .then(pl.col("midpoint") - 100)
              .otherwise(pl.col("start")).alias("start"),
            pl.when(pl.col("width") < 200)
              .then(pl.col("midpoint") + 100)
              .otherwise(pl.col("end")).alias("end"),
        ).drop(["width", "midpoint"])
        # ensure non‑negative starts
        peaks_pl = peaks_pl.with_columns(
            pl.when(pl.col("start") < 0).then(0).otherwise(pl.col("start")).alias("start")
        )

        # ------------------------------------------------------------------
        # 2. Optional merging ----------------------------------------------
        # ------------------------------------------------------------------
        peaks_bt = pybedtools.BedTool.from_dataframe(
            peaks_pl.to_pandas()
        ).sort()

        if merge_distance > 0:
            merged_bt = peaks_bt.merge(d=merge_distance)  # bedtools merge ‑d N
            peaks_bt = merged_bt.sort()
            logging.info(
                f"Merged peaks for '{group_name}' with d={merge_distance}. "
                f"Count now: {peaks_bt.count()}"
            )

        # Re‑load into Polars and reassign peak IDs
        merged_df = peaks_bt.to_dataframe(names=["chrom", "start", "end"])
        merged_df["peak_id"] = (
            merged_df["chrom"] + ":" + merged_df["start"].astype(str) + "-" + merged_df["end"].astype(str)
        )
        merged_df["group_name"] = group_name
        peaks_pl = pl.from_pandas(merged_df)

        # ------------------------------------------------------------------
        # 3. Intersections --------------------------------------------------
        # ------------------------------------------------------------------
        peaks_bt = pybedtools.BedTool.from_dataframe(merged_df).sort()
        for other_name, other_bed in all_group_bedgraphs:
            col_name = f"{other_name}_insertion_count"
            other_bt = pybedtools.BedTool(str(other_bed))
            intersect_result = peaks_bt.intersect(other_bt, wa=True, wb=True)

            if intersect_result.count() > 0:
                intersect_df = intersect_result.to_dataframe(
                    names=[
                        "chrom_peak", "start_peak", "end_peak", "peak_id",
                        "chrom_ins", "start_ins", "end_ins", "insertion_value",
                    ]
                )
                summed = (
                    intersect_df.groupby("peak_id")["insertion_value"].sum().reset_index()
                )
                summed.rename(columns={"insertion_value": col_name}, inplace=True)
                peaks_pl = peaks_pl.join(pl.from_pandas(summed), on="peak_id", how="left")
            else:
                peaks_pl = peaks_pl.with_columns(pl.lit(0).cast(pl.Int64).alias(col_name))

        # fill nulls -> 0
        count_cols = [f"{g}_insertion_count" for g, _ in all_group_bedgraphs]
        peaks_pl = peaks_pl.with_columns(
            [pl.col(c).fill_null(0).cast(pl.Int64) for c in count_cols if c in peaks_pl.columns]
        )
        peaks_pl = peaks_pl.with_columns((pl.col("end") - pl.col("start")).alias("peak_width"))

        logging.info(
            f"Processed intersections for '{group_name}' using {peaks_pl.height} peaks (merge d={merge_distance})."
        )
        return peaks_pl

    except Exception as e:
        logging.error(f"Error in worker for '{group_name}': {e}", exc_info=True)
        return None


# Helper for ProcessPoolExecutor ------------------------------------------------

def run_worker(task_args):
    return process_group_intersections(*task_args)

# ==============================================================================
# MAIN -------------------------------------------------------------------------
# ==============================================================================

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--peak_dir", required=True, help="Directory containing pre‑called SPAN peak files.")
    p.add_argument("--peak_suffix", required=True, help="Suffix of the peak files to use (e.g., '_b50.span.peak').")
    p.add_argument("--bedgraph_dir", required=True, help="Directory with per‑group insertion count bedGraph files.")
    p.add_argument("--output_dir", required=True, help="Directory to save output files.")
    p.add_argument("--annotation_file", required=True, help="Annotation file mapping samples to groups.")
    p.add_argument("--workers", type=int, default=12, help="Number of parallel worker processes.")
    p.add_argument("--merge_distance", type=int, default=200,
                   help="Distance (bp) within which to merge peaks after resizing; 0 disables.")

    args = p.parse_args()

    # --- Setup --------------------------------------------------------------
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    setup_logging(output_dir)

    pybedtools_tmp = output_dir / "pybedtools_tmp"
    pybedtools_tmp.mkdir(exist_ok=True)
    pybedtools.helpers.set_tempdir(str(pybedtools_tmp))

    # --- Prepare inputs -----------------------------------------------------
    ann = pl.read_csv(args.annotation_file, separator="\t" if args.annotation_file.endswith(".tsv") else ",")
    group_names = ann["group_name"].unique().to_list()

    all_bedgraphs = [
        (g, Path(args.bedgraph_dir) / f"{g}.bedgraph")
        for g in group_names if (Path(args.bedgraph_dir) / f"{g}.bedgraph").is_file()
    ]

    tasks = []
    for group_name, _ in all_bedgraphs:
        peak_file = Path(args.peak_dir) / f"{group_name}{args.peak_suffix}"
        if peak_file.is_file():
            tasks.append((group_name, peak_file, all_bedgraphs, args.merge_distance))
        else:
            logging.warning(f"Missing peak file for group '{group_name}' at {peak_file}. Skipping.")

    # --- Parallel processing -----------------------------------------------
    logging.info(f"Processing {len(tasks)} groups with {args.workers} workers (merge d={args.merge_distance}).")
    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as exe:
        for res in tqdm(exe.map(run_worker, tasks), total=len(tasks)):
            if res is not None:
                results.append(res)

    if not results:
        logging.error("No peak intersections processed. Exiting.")
        sys.exit(1)

    # --- Consolidate & save --------------------------------------------------
    per_group = pl.concat(results)
    per_group_dir = output_dir / "per_group_peak_matrices"
    per_group_dir.mkdir(exist_ok=True)

    for g in per_group["group_name"].unique().to_list():
        df_g = per_group.filter(pl.col("group_name") == g)
        if df_g.height > 0:
            df_g.write_csv(per_group_dir / f"{g}_peak_matrix.tsv", separator="\t")


if __name__ == "__main__":
    main()
