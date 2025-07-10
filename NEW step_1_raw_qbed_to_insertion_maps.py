# step_1_raw_qbed_to_insertion_maps.py
"""
An annotation-driven calling card data processing pipeline.

Major steps:
1.  Input loading, barcode correction, annotation, and per-sample partitioning.
2.  For each sample:
    a. Fragment-based peak calling.
    b. Generation of 1-bp unique insertion site maps.
3.  Aggregation of insertion maps by experimental group.
4.  Generation of bedGraph, BigWig, and qbed files for visualization and downstream analysis.

This script combines the functionalities of peak calling and insertion mapping into a single, streamlined workflow.
"""

# ==============================================================================
# LIBRARY IMPORTS
# ==============================================================================

import os
import sys
import argparse
import logging
import io
import glob
import shutil
import subprocess
import math
from pathlib import Path
from contextlib import redirect_stdout, redirect_stderr
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from typing import List, Sequence, Optional, Dict
import polars as pl
import pandas as pd
import numpy as np
import pybedtools
from tqdm import tqdm
from umi_tools import UMIClusterer
import pycallingcards as cc
import pyBigWig


# ==============================================================================
# SCRIPT CONFIGURATION
# ==============================================================================

VERSION = "3.0.0" # Updated version for fix

# Define the set of valid chromosome names to retain for downstream processing
VALID_CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

# Default parameters for pycallingcards peak calling
PEAK_CALLING_PARAMS_DEFAULTS = {
    'method': 'CCcaller',
    'reference': 'hg38',
    'pvalue_cutoff': 0.01,
    'pvalue_adj_cutoff': 0.01,
    'min_insertions': 5,
    'minlen': 0,
    'extend': 200,
    'maxbetween': 150,
    'minnum': 0,
    'lam_win_size': 1000000
}

# Schema expected for the raw input qbed/bed files
INPUT_QBED_SCHEMA = {
    'chrom': pl.Utf8,
    'start': pl.Int64,
    'end': pl.Int64,
    'reads': pl.Int64,
    'strand': pl.Utf8,
    'name': pl.Utf8
}

# ==============================================================================
# HELPER & UTILITY FUNCTIONS
# ==============================================================================

def setup_logging(output_dir: Path) -> None:
    """
    Initializes logging to both file and stdout for tracking progress and errors.
    """
    log_file_path = output_dir / "pipeline.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] - %(message)s",
        handlers=[
            logging.FileHandler(log_file_path, mode="w"),
            logging.StreamHandler(sys.stdout)
        ]
    )
    logging.info(f"Logging initialized. Log file: {log_file_path}")

def print_version_info() -> None:
    """
    Logs versions of the script and key dependencies.
    """
    import umi_tools
    logging.info(
        f"Pipeline v{VERSION} | Python v{sys.version.split()[0]} | Polars v{pl.__version__} | "
        f"PyCallingCards v{cc.__version__} | UMI-tools v{umi_tools.__version__} | pybedtools v{pybedtools.__version__}"
    )

def validate_args(args: argparse.Namespace) -> None:
    """
    Validates input arguments and ensures all input/output files and directories exist.
    Exits with an error if validation fails.
    """
    if not (1 <= args.workers <= os.cpu_count()):
        raise ValueError(f"Number of workers must be between 1 and {os.cpu_count()}.")
    if not Path(args.input_dir).is_dir():
        raise FileNotFoundError(f"Input directory not found: {args.input_dir}")
    if args.annotation_file and not Path(args.annotation_file).is_file():
        raise FileNotFoundError(f"Annotation file not found: {args.annotation_file}")
    if args.chrom_sizes and not Path(args.chrom_sizes).is_file():
        raise FileNotFoundError(f"Chromosome sizes file not found: {args.chrom_sizes}")
    if not shutil.which("bedGraphToBigWig"):
        logging.error("`bedGraphToBigWig` utility not found in your system's PATH. Please install it (part of UCSC utilities).")
        sys.exit(1)
    logging.info("All arguments are valid.")

def perform_sample_barcode_correction(
    raw_reads_df: pl.DataFrame,
    whitelist_barcodes: Sequence[str],
    correction_threshold: int
) -> pl.DataFrame:
    """
    Corrects raw sample barcodes to those in the annotation whitelist using Hamming distance.
    """
    unique_raw_barcodes = raw_reads_df["sample_barcode_raw"].unique().to_list()
    whitelist_set = set(whitelist_barcodes)
    barcode_correction_map = {}

    for raw_bc in unique_raw_barcodes:
        if raw_bc is None: continue
        if raw_bc in whitelist_set:
            barcode_correction_map[raw_bc] = raw_bc
            continue
        min_dist = float('inf')
        best_match = None
        for whitelist_bc in whitelist_barcodes:
            if len(raw_bc) != len(whitelist_bc): continue
            dist = sum(c1 != c2 for c1, c2 in zip(raw_bc, whitelist_bc))
            if dist < min_dist:
                min_dist = dist
                best_match = whitelist_bc
        if min_dist <= correction_threshold:
            barcode_correction_map[raw_bc] = best_match

    barcode_map_df = pl.DataFrame([
        {"sample_barcode_raw": k, "sample_barcode": v}
        for k, v in barcode_correction_map.items() if v is not None
    ])

    if barcode_map_df.is_empty():
        return raw_reads_df.clear()

    return raw_reads_df.join(barcode_map_df, on="sample_barcode_raw", how="inner")

def calculate_total_count_size_factors(group_maps: Dict[str, pl.DataFrame]) -> Dict[str, float]:
    """
    Calculates size factors for normalization and prints the total unique
    insertion count for each experimental group.
    """
    logging.info("Calculating Total Count size factors and printing total unique insertions...")
    
    # Calculate the sum of 'unique_count' for each group
    library_sizes = {name: df["unique_count"].sum() for name, df in group_maps.items()}
    
    # Create a DataFrame to display the total counts and log it
    total_counts_df = pl.DataFrame({
        "group": list(library_sizes.keys()),
        "total_unique_insertions": list(library_sizes.values())
    })
    logging.info("Total Unique Insertions per Group:\n%s", total_counts_df)
    
    # Continue with the original size factor calculation
    log_totals = np.log([count + 1 for count in library_sizes.values()])
    geo_mean_total = np.exp(np.mean(log_totals))
    
    if geo_mean_total <= 1:
        return {name: 1.0 for name in group_maps.keys()}
        
    factors = {name: total / geo_mean_total for name, total in library_sizes.items()}
    
    size_factor_df = pl.DataFrame({
        "group": list(factors.keys()),
        "size_factor": list(factors.values())
    })
    logging.info("Calculated Size Factors:\n%s", size_factor_df)
    
    return factors

def generate_binned_sum_track(group_name: str, normalized_bw_path: Path, window_size: int, output_dir_bw: Path) -> None:
    logging.info(f"--- Generating binned sum track for group '{group_name}' (window: {window_size} bp) ---")
    binned_bw_path = output_dir_bw / f"{group_name}.sum_{window_size}bp.bw"
    if binned_bw_path.exists():
        logging.info(f"Skipping existing binned track: {binned_bw_path.name}")
        return
    try:
        bw_in = pyBigWig.open(str(normalized_bw_path))
        if not bw_in: return
        bw_out = pyBigWig.open(str(binned_bw_path), "w")
        if not bw_out:
            bw_in.close()
            return
        bw_out.addHeader(list(bw_in.chroms().items()))
        for chrom, chrom_len in bw_in.chroms().items():
            if chrom_len == 0: continue
            arr = np.nan_to_num(bw_in.values(chrom, 0, chrom_len, numpy=True))
            num_bins = math.ceil(chrom_len / window_size)
            all_starts, all_ends, all_values = [], [], []
            for i in range(num_bins):
                start = i * window_size
                end = min(start + window_size, chrom_len)
                bin_sum = arr[start:end].sum()
                if bin_sum > 0:
                    all_starts.append(start); all_ends.append(end); all_values.append(float(bin_sum))
            if all_starts:
                bw_out.addEntries([chrom] * len(all_starts), all_starts, ends=all_ends, values=all_values)
        bw_in.close(); bw_out.close()
        logging.info(f"Successfully created binned BigWig: {binned_bw_path}")
    except Exception as e:
        logging.error(f"Binned sum generation for '{group_name}' failed: {e}", exc_info=True)

def convert_bedgraphs_to_qbed(bedgraph_dir: Path, qbed_dir: Path):
    logging.info(f"--- Converting bedGraphs in '{bedgraph_dir.name}' to qbeds in '{qbed_dir.name}' ---")
    qbed_dir.mkdir(parents=True, exist_ok=True)
    bedgraph_files = sorted(bedgraph_dir.glob("*.bedgraph"))
    if not bedgraph_files:
        logging.warning(f"No .bedgraph files found in {bedgraph_dir} to convert.")
        return
    for bedgraph_file in tqdm(bedgraph_files, desc="Converting to qbed"):
        qbed_path = qbed_dir / (bedgraph_file.stem + ".qbed")
        if qbed_path.exists():
            continue
        df = pd.read_csv(bedgraph_file, sep="\t", header=None, names=["chrom", "start", "end", "count"], dtype={"chrom": str, "start": int, "end": int, "count": int})
        if df.empty: continue
        df_expanded = df.loc[df.index.repeat(df["count"])]
        df_expanded[["chrom", "start", "end"]].to_csv(qbed_path, sep="\t", index=False, header=False)
    logging.info("✅ Finished qbed conversion.")


# ==============================================================================
# CORE WORKER FUNCTION
# ==============================================================================

def process_pooled_sample_data(task: tuple, annotation_df: pl.DataFrame, output_dir: Path) -> Optional[pl.DataFrame]:
    """
    Processes pooled fragments for a single sample. This is a hybrid function.
    - Prepares data exactly like the original step_1 script.
    - Calls fragment-based peaks using pycallingcards.
    - Determines precise 1-bp insertion sites using the logic from the original step_2 script.
    - Returns a DataFrame of unique insertion sites for the sample (format: bedGraph).
    """
    sample_parquet_path_str, sample_name, library_name_concat, peak_calling_params, srt_bc_clustering_threshold = task
    sample_parquet_path = Path(sample_parquet_path_str)

    try:
        # ====== 1. LOAD AND PREPARE SAMPLE DATA ======
        sample_fragments_df = pl.read_parquet(sample_parquet_path)
        sample_fragments_df = sample_fragments_df.with_columns(pl.lit(library_name_concat).alias("library_name"))
        
        row = annotation_df.filter(pl.col("sample_name") == sample_name)
        group_name = row["group_name"][0] if row.height > 0 else None
        
        logging.info(f"[{sample_name}] Performing fragment-level deduplication...")
        fragment_deduplication_grouping_keys = [
            col for col in sample_fragments_df.columns if col not in ["sample_barcode_raw", "reads"]
        ]
        initial_fragment_count = sample_fragments_df.height
        sample_fragments_df = sample_fragments_df.group_by(fragment_deduplication_grouping_keys).agg(
            pl.sum("reads").alias("reads"),
            pl.col("sample_barcode_raw").first().alias("sample_barcode_raw")
        )
        logging.info(f"[{sample_name}] Fragments before deduplication: {initial_fragment_count}, after: {sample_fragments_df.height}")

        # ====== 2. FRAGMENT-BASED PEAK CALLING ======
        pycc_fragments_df = (
            sample_fragments_df
            .select(["chrom", "start", "end"])
            .to_pandas()
            .rename(columns={"chrom": "Chr", "start": "Start", "end": "End"})
        )
        with redirect_stderr(io.StringIO()), redirect_stdout(io.StringIO()):
            initial_peaks_pd = cc.pp.call_peaks(expdata=pycc_fragments_df, **peak_calling_params)

        if initial_peaks_pd is None or initial_peaks_pd.empty:
            logging.warning(f"[{sample_name}] No fragment peaks were called.")
            return None

        initial_peaks_df = pl.from_pandas(initial_peaks_pd.iloc[:, :3])
        initial_peaks_df.columns = ["chrom", "fragment_peak_start_for_sample", "fragment_peak_end_for_sample"]
        initial_peaks_df = initial_peaks_df.with_columns(
            (pl.lit(f"{sample_name}_peak_") + pl.arange(0, initial_peaks_df.height).cast(pl.Utf8)).alias("sample_peak_id")
        )

        fragment_peaks_dir = output_dir / "fragment_peaks_per_sample"
        fragment_peaks_dir.mkdir(exist_ok=True, parents=True)
        fragment_peaks_path = fragment_peaks_dir / f"{sample_name}_fragment_peaks.parquet"
        initial_peaks_df.write_parquet(fragment_peaks_path)
        logging.info(f"[{sample_name}] Saved {initial_peaks_df.height} fragment peaks to {fragment_peaks_path.name}")


        # ====== 3. MAP FRAGMENTS TO PEAKS ======
        bed_fragments_df = sample_fragments_df.select(
            ["chrom", "start", "end", "strand", "srt_bc", "reads"]
        ).to_pandas()
        bed_fragments = pybedtools.BedTool.from_dataframe(bed_fragments_df)
        bed_peaks = pybedtools.BedTool.from_dataframe(initial_peaks_df.to_pandas())

        fragment_columns = bed_fragments_df.columns.tolist()
        peak_columns = [f"{col}_peak" for col in initial_peaks_df.columns]
        intersection_columns = fragment_columns + peak_columns
        
        intersected_df = pl.from_pandas(
            bed_fragments.intersect(bed_peaks, wa=True, wb=True).to_dataframe(names=intersection_columns)
        )
        if intersected_df.is_empty():
            logging.warning(f"[{sample_name}] No fragments intersected with peaks.")
            return None




        # ============================================================ #
        # CRITICAL LOGIC OF DEFINING THE 1-BP INSERTION SITE           #
        # ============================================================ #
        sites = []
        for _, peak_group_df in intersected_df.group_by("sample_peak_id_peak"):
            srt_bc_counts = peak_group_df.group_by("srt_bc").agg(pl.sum("reads"))
            srt_bc_read_count_dict = {
                row['srt_bc']: row['reads']
                for row in srt_bc_counts.to_dicts() if row['srt_bc']
            }
            if not srt_bc_read_count_dict:
                continue

            # Cluster SRT barcodes by Hamming distance with UMI tools library
            clusterer = UMIClusterer(cluster_method="directional")
            clustered_srt_bcs = clusterer(
                {k.encode('utf-8'): v for k, v in srt_bc_read_count_dict.items()},
                threshold=srt_bc_clustering_threshold
            )
            if not clustered_srt_bcs:
                continue
            
            # Map each SRT to a cluster representative (deduplication)
            srt_bc_rep_map = {
                srt_bc.decode('utf-8'): cluster[0].decode('utf-8')
                for cluster in clustered_srt_bcs for srt_bc in cluster
            }
            
            # Collapse fragments by coordinates, strand, and UMI cluster (unique molecules)
            deduped_frags_df = (
                peak_group_df
                .with_columns(pl.col("srt_bc").replace(srt_bc_rep_map).alias("srt_bc_rep"))
                .filter(pl.col("srt_bc_rep").is_not_null())
                .group_by("chrom", "start", "end", "strand", "srt_bc_rep")
                .agg(pl.sum("reads").alias("reads_in_molecule"))
            )

            if deduped_frags_df.height < 5:
                # Enforce minimum unique fragments in peak (removes noisy peaks)
                continue

            summary = (
                deduped_frags_df.group_by("srt_bc_rep")
                .agg([
                    pl.col("start").max().alias("max_start"),
                    pl.col("end").min().alias("min_end"),
                    pl.first("strand").alias("strand"),
                    pl.first("chrom").alias("chrom")
                ])
                .with_columns(
                    # For each SRT barcode, identify the single coordinate most proximal to the actual insertion event.
                    
                    pl.when(pl.col("strand") == "+")
                    .then(pl.col("min_end"))
                    # For '+' strand fragments, the reads move 5' <- 3' (opposite of convention for + strand).
                    # In the input qbed, the 'end' coordinate is the true end of the read from the sequencer for '+' strand fragments.
                    # Therefore, if the SRT barcode is on '+' strand fragments, the smallest 'end' coordinate is selected.
                    
                    # For '-' strand fragments ('.otherwise' means it is '-' strand since there are only two strands, + and -), the reads move 5' -> 3' (opposite of convention for - strand).
                    # In the input qbed, the 'start' coordinate is the true end of the read from the sequencer for '-' strand fragments.
                    # Therefore, if the SRT barcode is on '-' strand fragments, the largest 'start' coordinate is selected.
                    .otherwise(pl.col("max_start"))
                    .alias("site")
                )
                .select(["chrom", "site"])
            )
            sites.append(summary)

        if not sites:
            logging.warning(f"[{sample_name}] No unique insertion sites were identified after filtering.")
            return None
        
        final_sites_df = pl.concat(sites)
        
        return (
            final_sites_df.group_by(["chrom", "site"])
            .agg(pl.count().alias("unique_count"))
            .rename({"site": "start"})
            .with_columns((pl.col("start") + 1).alias("end"))
            .select(["chrom", "start", "end", "unique_count"])
            .sort(["chrom", "start"])
        )

    except Exception as exc:
        logging.error(f"Error processing {sample_parquet_path.name}: {exc}", exc_info=True)
        return None
        
        
        
# ==============================================================================
# MAIN PIPELINE
# ==============================================================================

def main():
    """ Main function to orchestrate the entire workflow. """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--input_dir", required=True, help="Directory containing input qbed/bed files.")
    parser.add_argument("--output_dir", required=True, help="Directory to store all output files and subdirectories.")
    parser.add_argument("--annotation_file", required=True, type=str, help="Path to a 4-column TSV or CSV file with header: 'library_name,sample_barcode,sample_name,group_name'.")
    parser.add_argument("--chrom_sizes", required=True, help="Path to a chromosome sizes file (e.g., hg38.chrom.sizes).")
    parser.add_argument("--workers", type=int, default=10, help="Number of parallel workers.")
    parser.add_argument("--srt_bc_dist_threshold", type=int, default=1, help="Max Hamming distance for SRT-BC (UMI) clustering.")
    parser.add_argument("--sample_barcode_dist_threshold", type=int, default=2, help="Max Hamming distance for sample barcode correction.")
    parser.add_argument("--min_rows_threshold", type=int, default=50000, help="Min fragments for a sample to be processed.")
    parser.add_argument("--sum_window_size", type=int, default=50, help="Window size in bp for binned summary BigWig tracks.")

    for param_key, param_default in PEAK_CALLING_PARAMS_DEFAULTS.items():
        parser.add_argument(f"--{param_key}", type=type(param_default), default=param_default)
    args = parser.parse_args()

    # ==== Setup ====
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    setup_logging(output_dir)
    print_version_info()
    validate_args(args)
    
    pybedtools_temp_dir = output_dir / "pybedtools_tmp"
    pybedtools_temp_dir.mkdir(exist_ok=True, parents=True)
    pybedtools.helpers.set_tempdir(pybedtools_temp_dir)

    # --- Phase 1: Load, Correct, Annotate, and Partition Reads ---
    logging.info("--- Phase 1: Loading, Correcting, Annotating, and Partitioning Reads ---")
    
    # Using original step_1 logic for partitioning: No extra cleaning.
    annotation_file_separator = '\t' if args.annotation_file.endswith(('.tsv', '.txt')) else ','
    annotation_df = pl.read_csv(args.annotation_file, separator=annotation_file_separator, has_header=True)
    
    whitelist_sample_barcodes = annotation_df["sample_barcode"].unique().to_list()
    logging.info(f"Derived a whitelist of {len(whitelist_sample_barcodes)} barcodes from the annotation file.")

    input_files = [p for ext in ("*.qbed", "*.bed", "*.qbed.gz", "*.bed.gz") for p in Path(args.input_dir).glob(ext)]
    if not input_files:
        logging.error(f"No input files found directly in {args.input_dir}.")
        sys.exit(1)

    partitioned_temp_dir = output_dir / "partitioned_temp"
    partitioned_temp_dir.mkdir(exist_ok=True, parents=True)
    
    observed_sample_library_pairs = set()

    logging.info(f"Processing {len(input_files)} input files...")
    for file_path in tqdm(input_files, desc="Reading Input Files and Partitioning"):
        try:
            fragment_df = pl.read_csv(
                file_path, has_header=False, separator="\t",
                new_columns=list(INPUT_QBED_SCHEMA.keys()),
                schema_overrides=INPUT_QBED_SCHEMA, ignore_errors=True
            )
            
            fragment_df = (
                fragment_df
                .filter(pl.col("chrom").is_in(VALID_CHROMOSOMES) & pl.col("name").is_not_null())
                .with_columns([
                    pl.col("name").str.split("/").list.get(0).alias("library_name"),
                    pl.col("name").str.split("/").list.get(1).alias("sample_barcode_raw"),
                    pl.col("name").str.split("/").list.get(2).alias("srt_bc")
                ])
                .drop("name")
                .filter(pl.col("sample_barcode_raw").is_not_null())
            )
            if fragment_df.is_empty(): continue

            corrected_reads_df = perform_sample_barcode_correction(
                fragment_df,
                whitelist_barcodes=whitelist_sample_barcodes,
                correction_threshold=args.sample_barcode_dist_threshold
            )
            if corrected_reads_df.is_empty(): continue

            # Using original step_1 join logic
            annotated_reads_df = corrected_reads_df.join(
                annotation_df, on=["library_name", "sample_barcode"], how="inner"
            )
            if annotated_reads_df.is_empty(): continue

            pairs_in_file = annotated_reads_df.select(["sample_name", "library_name"]).unique()
            for row in pairs_in_file.iter_rows():
                observed_sample_library_pairs.add(row)
            
            for sample_name_tuple, group_df in annotated_reads_df.group_by("sample_name"):
                sample_name = sample_name_tuple[0]
                temp_out_path = partitioned_temp_dir / f"{sample_name}_{file_path.stem}.parquet"
                group_df.write_parquet(temp_out_path)
        except Exception as exc:
            logging.error(f"Failed to process file {file_path}: {exc}", exc_info=True)

    # --- Phase 2: Consolidate Partitions and Process Samples ---
    logging.info("--- Phase 2: Consolidating and Processing Samples ---")
    collapsed_dir = output_dir / "collapsed_fragments_per_sample"
    collapsed_dir.mkdir(exist_ok=True, parents=True)

    sample_to_libs_map = {}
    for sample_name, library_name in observed_sample_library_pairs:
        if sample_name not in sample_to_libs_map:
            sample_to_libs_map[sample_name] = set()
        sample_to_libs_map[sample_name].add(library_name)

    sample_name_to_librarystr = {
        sample: "__".join(sorted(list(libs)))
        for sample, libs in sample_to_libs_map.items()
    }
    
    unique_sample_names = sorted(list(sample_to_libs_map.keys()))
    
    for sample_name in tqdm(unique_sample_names, desc="Consolidating Samples"):
        sample_files = list(partitioned_temp_dir.glob(f"{sample_name}_*.parquet"))
        if not sample_files: continue
        
        sample_df = pl.read_parquet(sample_files)
        
        out_path = collapsed_dir / f"{sample_name}.parquet"
        if sample_df.height >= args.min_rows_threshold:
            sample_df.write_parquet(out_path)
        else:
            logging.warning(f"Skipping sample '{sample_name}': has {sample_df.height} reads, less than threshold of {args.min_rows_threshold}.")
    
    shutil.rmtree(partitioned_temp_dir)
    collapsed_files = sorted(list(collapsed_dir.glob("*.parquet")))
    if not collapsed_files:
        logging.error(f"No samples passed the min_rows_threshold. Exiting.")
        sys.exit(1)

    # --- Run per-sample peak calling and insertion mapping ---
    peak_params = {k: getattr(args, k) for k in PEAK_CALLING_PARAMS_DEFAULTS}
    tasks = [
        (str(f), f.stem, sample_name_to_librarystr[f.stem], peak_params, args.srt_bc_dist_threshold)
        for f in collapsed_files if f.stem in sample_name_to_librarystr
    ]
    
    sample_maps: Dict[str, pl.DataFrame] = {}
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(process_pooled_sample_data, task, annotation_df, output_dir): task for task in tasks}
        for f in tqdm(as_completed(futures), total=len(futures), desc="Calling Peaks & Mapping Insertions"):
            result_df = f.result()
            task_info = futures[f]
            sample_name = task_info[1]
            if result_df is not None and not result_df.is_empty():
                sample_maps[sample_name] = result_df

    if not sample_maps:
        logging.error("No valid insertion maps were generated from any sample. Exiting.")
        sys.exit(1)

    # --- Phase 3: Generate Insertion Map Outputs ---
    logging.info("--- Phase 3: Generating Final Insertion Map Outputs ---")
    
    dirs = {
        "sample_bg": output_dir / "raw_unique_insertions_per_sample_bedgraph",
        "sample_bw": output_dir / "raw_unique_insertions_per_sample_bigwig",
        "group_bg": output_dir / "raw_unique_insertion_count_per_group_bedgraph",
        "group_bw": output_dir / "raw_unique_insertion_count_per_group_bigwig",
        "norm_bg": output_dir / "size_normalized_unique_insertions_per_group_bedgraph",
        "norm_bw": output_dir / "size_normalized_unique_insertions_per_group_bigwig",
        "binned_sum_bw": output_dir / f"binned_normalized_count_sum_bigwig_window_{args.sum_window_size}bp" if args.sum_window_size else None,
        "qbed": output_dir / "raw_unique_insertion_count_per_group_qbed"
    }
    for d in dirs.values():
        if d: d.mkdir(exist_ok=True, parents=True)
        
    # Using original step_2 logic for grouping: Read and clean annotation file.
    logging.info("Loading and cleaning annotation for group aggregation...")
    cleaned_ann_df = pl.read_csv(
        args.annotation_file,
        separator='\t' if args.annotation_file.endswith(('.tsv', '.txt')) else ',',
        has_header=True
    )
    # This logic is from step_2
    cleaned_ann_df = cleaned_ann_df.rename(
        {col: col.strip() for col in cleaned_ann_df.columns}
    ).with_columns(pl.col("sample_name").str.strip_chars())

    group_to_samples: List[tuple[str, List[str]]] = cleaned_ann_df.group_by("group_name").agg(pl.col("sample_name")).rows()

    # --- STAGE 3a: Write Sample BedGraphs & BigWigs ---
    logging.info("--- Generating Per-Sample bedGraphs and BigWigs ---")
    chrom_sizes_path = Path(args.chrom_sizes)
    for sample_name, df in tqdm(sample_maps.items(), desc="Saving Sample Maps"):
        bg_path = dirs["sample_bg"] / f"{sample_name}.bedgraph"
        bw_path = dirs["sample_bw"] / f"{sample_name}.bw"
        df.write_csv(bg_path, separator="\t", include_header=False)
        if not bw_path.exists():
            try:
                subprocess.run(["bedGraphToBigWig", str(bg_path), str(chrom_sizes_path), str(bw_path)], check=True, capture_output=True, text=True)
            except subprocess.CalledProcessError as e:
                logging.error(f"Failed to create BigWig for sample {sample_name}: {e.stderr.strip()}")

    # --- STAGE 3b: Aggregate to Group BedGraphs & BigWigs ---
    logging.info("--- Aggregating Samples into Groups ---")
    group_maps: Dict[str, pl.DataFrame] = {}
    for group_name, sample_name_list in tqdm(group_to_samples, desc="Aggregating Groups"):
        samp_list = [s for s in sample_name_list if s in sample_maps]
        if not samp_list: continue
        
        concat = pl.concat([sample_maps[s] for s in samp_list])
        agg = concat.group_by(["chrom", "start", "end"]).agg(pl.sum("unique_count")).sort(["chrom", "start"])
        group_maps[group_name] = agg

        bg_path = dirs["group_bg"] / f"{group_name}.bedgraph"
        bw_path = dirs["group_bw"] / f"{group_name}.bw"
        agg.write_csv(bg_path, separator="\t", include_header=False)
        if not bw_path.exists():
            try:
                subprocess.run(["bedGraphToBigWig", str(bg_path), str(chrom_sizes_path), str(bw_path)], check=True, capture_output=True, text=True)
            except subprocess.CalledProcessError as e:
                logging.error(f"Failed to create BigWig for group {group_name}: {e.stderr.strip()}")

    # --- STAGE 3c: Normalize Group Maps ---
    logging.info("--- Normalizing Group Maps ---")
    if len(group_maps) > 1:
        factors = calculate_total_count_size_factors(group_maps)
        for g, df in tqdm(group_maps.items(), desc="Normalizing Group Maps"):
            fct = factors.get(g, 1.0)
            if fct > 0:
                norm_df = df.with_columns((pl.col("unique_count") / fct).alias("normalized_count")).select(["chrom", "start", "end", "normalized_count"])
                bg_path = dirs["norm_bg"] / f"{g}.bedgraph"
                bw_path = dirs["norm_bw"] / f"{g}.bw"
                norm_df.write_csv(bg_path, separator="\t", include_header=False)
                if not bw_path.exists():
                    try:
                        subprocess.run(["bedGraphToBigWig", str(bg_path), str(chrom_sizes_path), str(bw_path)], check=True, capture_output=True, text=True)
                    except subprocess.CalledProcessError as e:
                        logging.error(f"Failed to create normalized BigWig for {g}: {e.stderr.strip()}")
    else:
        logging.warning("Skipping normalization as only one group was found.")

    # --- STAGE 3d: Convert Group BedGraphs to Qbeds ---
    convert_bedgraphs_to_qbed(bedgraph_dir=dirs["group_bg"], qbed_dir=dirs["qbed"])

    # --- STAGE 3e: Binned Sum Track Generation ---
    if args.sum_window_size and args.sum_window_size > 0 and dirs["binned_sum_bw"]:
        logging.info("--- Generating Binned Sum Tracks ---")
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            for group_name, _ in group_to_samples:
                bw_path = dirs["norm_bw"] / f"{group_name}.bw"
                if bw_path.is_file():
                    executor.submit(generate_binned_sum_track, group_name, bw_path, args.sum_window_size, dirs["binned_sum_bw"])

    logging.info("Pipeline complete")
    pybedtools.helpers.cleanup()

if __name__ == "__main__":
    main()
