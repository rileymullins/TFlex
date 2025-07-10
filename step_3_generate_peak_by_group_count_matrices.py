# step_3_generate_peak_by_group_count_matrices.py
"""
Final peak-processing script for a calling card pipeline.

This script takes pre-called peaks for each group, resizes them, and generates
per-group intersection matrices. It then creates a final consensus peak set and
uses it to generate a final consensus count matrix with insertion counts from all groups.
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


def process_group_intersections(
    group_name: str,
    peak_file_path: Path,
    all_group_bedgraphs: List[Tuple[str, Path]]
) -> Optional[pl.DataFrame]:
    """
    Worker: (1) Loads and makes the min size of the peaks 200bp per group, and then in a loop:
    (2) Intersects those peaks with all groups' raw bedgraph insertions
        and sums the insertion counts (column 4) for each peak.
        
    This will result in a matrix of each group's peaks as rows and the unique insertion count of each groups by column.
    """
    try:
        # ----------------------------------------------------------------------
        # 1. Load and Resize the Group's Own Peaks
        # ----------------------------------------------------------------------
        
        # Load the 3-column peak file for the current group
        peaks_pl = pl.read_csv(
            peak_file_path, separator="\t", has_header=False,
            new_columns=["chrom", "start", "end"]
        ).select(["chrom", "start", "end"])

        if peaks_pl.is_empty():
            logging.warning(f"No peaks in file for '{group_name}'. Skipping.")
            return None

        # Calculate the midpoint of each peak
        peaks_pl = peaks_pl.with_columns(
            (pl.col("end") - pl.col("start")).alias("width")
        ).with_columns(
            (pl.col("start") + (pl.col("width") / 2)).cast(pl.Int64).alias("midpoint")
        )

        # Resize any peak smaller than 200bp to be 200bp centered at its midpoint. 
        # This adds a window to idenitfy other group's insertions in the peak if it is very small.
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

        # Ensure no peak coordinates are less than 0 after resizing.
        peaks_pl = peaks_pl.with_columns(
            pl.when(pl.col("start") < 0).then(0).otherwise(pl.col("start")).alias("start")
        )

        # Add metadata columns: the group this peak set belongs to and a unique peak_id.
        peaks_pl = peaks_pl.with_columns([
            pl.lit(group_name).alias("group_name"),
            (pl.col("chrom") + ":" + pl.col("start").cast(str) + "-" + pl.col("end").cast(str)).alias("peak_id")
        ])
        
        # Convert the Polars DataFrame to a BedTool object for efficient intersection
        peaks_bt = pybedtools.BedTool.from_dataframe(
            peaks_pl.select(["chrom", "start", "end", "peak_id"]).to_pandas()
        ).sort()

        # ----------------------------------------------------------------------
        # 2. Intersect These Peaks Against All Groups' Insertions
        # ----------------------------------------------------------------------
        # This loop creates an "all-vs-all" matrix from the perspective of one group's peaks.
        for other_name, other_bed in all_group_bedgraphs:
            col_name = f"{other_name}_insertion_count"
            other_insertions_bt = pybedtools.BedTool(str(other_bed))
            
            # Intersect this group's peaks with another group's insertions.
            # wa=True and wb=True to store the original A and B features for each overlap.
            intersect_result = peaks_bt.intersect(other_insertions_bt, wa=True, wb=True)

            # Check if any intersections were found
            if intersect_result.count() > 0:
                # If intersections were found, convert results to a pandas DataFrame for aggregation
                intersect_df = intersect_result.to_dataframe(
                    names=['chrom_peak', 'start_peak', 'end_peak', 'peak_id',
                           'chrom_ins', 'start_ins', 'end_ins', 'insertion_value']
                )
                # For each peak, sum the counts of all insertions that fell within it
                summed_counts = intersect_df.groupby('peak_id')['insertion_value'].sum().reset_index()
                summed_counts.rename(columns={'insertion_value': col_name}, inplace=True)
                
                # Join the summed counts back to the main peak table
                peaks_pl = peaks_pl.join(pl.from_pandas(summed_counts), on='peak_id', how='left')
            else:
                # If no intersections, create a column of zeros for this group comparison
                peaks_pl = peaks_pl.with_columns(pl.lit(0).cast(pl.Int64).alias(col_name))

        # ----------------------------------------------------------------------
        # 3. Final Cleanup and Return
        # ----------------------------------------------------------------------
        # After joining, peaks may with no overlaps with have null values. In that case, fill with 0.
        all_count_cols = [f"{g}_insertion_count" for g, _ in all_group_bedgraphs]
        peaks_pl = peaks_pl.with_columns(
            [pl.col(c).fill_null(0).cast(pl.Int64) for c in all_count_cols if c in peaks_pl.columns]
        )
        
        # Add a peak width column
        peaks_pl = peaks_pl.with_columns((pl.col("end") - pl.col("start")).alias("peak_width"))

        logging.info(f"Processed intersections for '{group_name}' using {peaks_pl.height} resized peaks.")
        return peaks_pl

    except Exception as e:
        logging.error(f"Error in worker for '{group_name}': {e}", exc_info=True)
        return None


def run_worker(task_args):
    """Helper function to unpack arguments for the first ProcessPoolExecutor."""
    return process_group_intersections(*task_args)



def main():
    """Main function to orchestrate the peak counting and consensus workflow."""
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--peak_dir", required=True, help="Directory containing pre-called SPAN peak files.")
    p.add_argument("--peak_suffix", required=True, help="Suffix of the peak files to use (e.g., '_b50.span.peak').")
    p.add_argument("--bedgraph_dir", required=True, help="Directory with per-group insertion count bedGraph files.")
    p.add_argument("--output_dir", required=True, help="Directory to save output files.")
    p.add_argument("--annotation_file", required=True, help="Annotation file mapping samples to groups.")
    p.add_argument("--workers", type=int, default=12, help="Number of parallel worker processes.")
    args = p.parse_args()

    # --- Setup ---
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    setup_logging(output_dir)
    pybedtools_tmp_dir = output_dir / "pybedtools_tmp"
    pybedtools_tmp_dir.mkdir(exist_ok=True)
    pybedtools.helpers.set_tempdir(str(pybedtools_tmp_dir))

    # --- Prepare Inputs ---
    # Get all unique group names from the annotation file
    ann = pl.read_csv(args.annotation_file, separator="\t" if args.annotation_file.endswith(".tsv") else ",")
    group_names = ann["group_name"].unique().to_list()
    
    # Find all existing bedgraph files for the groups
    all_bedgraphs = [(g, Path(args.bedgraph_dir) / f"{g}.bedgraph") for g in group_names if (Path(args.bedgraph_dir) / f"{g}.bedgraph").is_file()]

    # Prepare tasks for the first parallel step: generating per-group matrices
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
        # exe.map efficiently applies the run_worker function to each task in parallel
        future_results = exe.map(run_worker, tasks)
        for res in tqdm(future_results, total=len(tasks)):
            if res is not None:
                results.append(res)

    if not results:
        logging.error("No peak intersections could be processed. Exiting.")
        sys.exit(1)

    # --- Consolidate and Save Per-Group Results ---
    # Combine the list of DataFrames from each worker into a single large DataFrame
    per_group = pl.concat(results)
    
    # Save a separate count matrix for each group's peak set
    per_group_dir = output_dir / "per_group_peak_matrices"
    per_group_dir.mkdir(exist_ok=True)
    for g in per_group["group_name"].unique().to_list():
        df_g = per_group.filter(pl.col("group_name") == g)
        if df_g.height > 0:
            df_g.write_csv(per_group_dir / f"{g}_peak_matrix.tsv", separator="\t")

if __name__ == '__main__':
    main()
