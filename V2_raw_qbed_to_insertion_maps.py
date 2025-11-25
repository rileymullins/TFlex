# V2_raw_qbed_to_insertion_maps.py
"""
Major steps:
1.  Input loading, barcode correction, annotation.
2.  Filtering: Remove reads overlapping Blacklist (+/- 150bp).
3.  Partitioning: Split data by Sample and Chromosome.
4.  Processing per sample and chromosome: 
    - Bin fragments based on a max gap size to next fragment (controlled by --variable_bin_gap).
    - Correct SRT barcode with UMI correction tools.
    - Deduplicate fragments of the same insertion based on matching SRT barcode and strand with the same variable width bin.
5.  Outputs:
    - Per-sample variable width bin BEDs with fragment counts.
    - QBED where each row is a unique insertion.
    - QBED with stats where unique insertions at the same site and strand are collapsed stats.
        - Columns are chrom, start, end, unique_count, total_reads, strand.
    - Per-sample and per-group raw and counts per million (CPM) insertion count versions of the following:
        - Bedgraph (unique insertions at the same site regardless of strand are collapsed to one row, and the fourth column is insertion count).
        - Binned BigWig (bin sized controlled by --sum_window_size, default 100bp).
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
import pyBigWig


# ==============================================================================
# SCRIPT CONFIGURATION
# ==============================================================================

VERSION = "2.0.0"

# Define the set of valid chromosome names.
VALID_CHROMOSOMES = set([f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"])

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
    import umi_tools
    logging.info(
        f"Pipeline v{VERSION} | Python v{sys.version.split()[0]} | Polars v{pl.__version__} | "
        f"UMI-tools v{umi_tools.__version__} | pybedtools v{pybedtools.__version__}"
    )

def validate_args(args: argparse.Namespace) -> None:
    if not (1 <= args.workers <= os.cpu_count()):
        raise ValueError(f"Number of workers must be between 1 and {os.cpu_count()}.")
    if not Path(args.input_dir).is_dir():
        raise FileNotFoundError(f"Input directory not found: {args.input_dir}")
    if args.annotation_file and not Path(args.annotation_file).is_file():
        raise FileNotFoundError(f"Annotation file not found: {args.annotation_file}")
    if args.chrom_sizes and not Path(args.chrom_sizes).is_file():
        raise FileNotFoundError(f"Chromosome sizes file not found: {args.chrom_sizes}")
    if args.blacklist_file and not Path(args.blacklist_file).is_file():
        raise FileNotFoundError(f"Blacklist file not found: {args.blacklist_file}")
    if not shutil.which("bedGraphToBigWig"):
        logging.error("`bedGraphToBigWig` utility not found in your system's PATH.")
        sys.exit(1)
    logging.info("All arguments are valid.")

def load_chrom_sizes(path: Path) -> Dict[str, int]:
    """Loads chromosome sizes into a dictionary for clipping."""
    sizes = {}
    with open(path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                sizes[parts[0]] = int(parts[1])
    return sizes

def safe_write_1bp_bedgraph(df: pl.DataFrame, bg_path: Path, score_col: str = "unique_count"):
    """Writes a 1-bp resolution BedGraph."""
    clean_df = (
        df.group_by(["chrom", "start", "end"])
        .agg(pl.sum(score_col))
        .sort(["chrom", "start"])
    )
    if clean_df.height == 0: return
    clean_df.select(["chrom", "start", "end", score_col]).write_csv(bg_path, separator="\t", include_header=False)

def bin_and_write_bigwig(df: pl.DataFrame, bw_path: Path, chrom_sizes_path: Path, chrom_sizes_map: Dict[str, int], score_col: str = "unique_count", bin_size: int = 50):
    """Bins data into fixed windows, clips coordinates, and converts to BigWig."""
    binned_df = (
        df.with_columns((pl.col("start") // bin_size * bin_size).alias("bin_start"))
        .group_by(["chrom", "bin_start"])
        .agg(pl.sum(score_col).alias("score"))
        .with_columns((pl.col("bin_start") + bin_size).alias("bin_end"))
        .rename({"bin_start": "start", "bin_end": "end"})
        .select(["chrom", "start", "end", "score"])
        .sort(["chrom", "start"])
    )
    
    size_df = pl.DataFrame({"chrom": list(chrom_sizes_map.keys()), "max_len": list(chrom_sizes_map.values())})
    binned_df = (
        binned_df.join(size_df, on="chrom", how="inner")
        .with_columns(pl.min_horizontal("end", "max_len").alias("end"))
        .filter(pl.col("start") < pl.col("end"))
        .drop("max_len")
    )

    if binned_df.height == 0: return

    temp_bg = bw_path.with_suffix(".temp.bedgraph")
    binned_df.write_csv(temp_bg, separator="\t", include_header=False)
    try:
        subprocess.run(["bedGraphToBigWig", str(temp_bg), str(chrom_sizes_path), str(bw_path)], check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"bedGraphToBigWig FAILED for {bw_path.name}")
        logging.error(f"STDERR: {e.stderr.strip()}")
    finally:
        if temp_bg.exists(): temp_bg.unlink()

def clip_coordinates(df: pl.DataFrame, chrom_sizes_map: Dict[str, int]) -> pl.DataFrame:
    """ Clips 'end' coordinates to ensure they do not exceed chromosome size. """
    size_df = pl.DataFrame({"chrom": list(chrom_sizes_map.keys()), "max_len": list(chrom_sizes_map.values())})
    return (
        df.join(size_df, on="chrom", how="inner")
        .with_columns(pl.min_horizontal("end", "max_len").alias("end"))
        .filter(pl.col("start") < pl.col("end"))
        .drop("max_len")
    )

def normalize_and_filter_chromosomes(df: pl.DataFrame, valid_chroms: set) -> pl.DataFrame:
    return (
        df
        .with_columns(pl.col("chrom").str.replace("GRCh38_", ""))
        .with_columns(
            pl.when(pl.col("chrom").str.starts_with("chr").not_())
            .then(pl.lit("chr") + pl.col("chrom"))
            .otherwise(pl.col("chrom"))
            .alias("chrom")
        )
        .filter(pl.col("chrom").is_in(valid_chroms))
    )

def filter_dataframe_by_blacklist(df: pl.DataFrame, blacklist_bt: pybedtools.BedTool) -> pl.DataFrame:
    """
    Filters rows in Polars DataFrame that overlap with regions in the blacklist BedTool.
    Uses 'intersect -v' logic (keep only rows that do NOT overlap).
    """
    if df.is_empty():
        return df

    # Convert Polars DataFrame to CSV string for BedTool creation
    # This avoids writing to disk for every chunk
    csv_string = df.write_csv(separator="\t", include_header=False)
    
    try:
        input_bt = pybedtools.BedTool(csv_string, from_string=True)
        
        # Intersect with -v to get only entries NOT overlapping the blacklist
        filtered_bt = input_bt.intersect(blacklist_bt, v=True)
        
        if len(filtered_bt) == 0:
            return df.clear()
            
        # Read back into Polars
        filtered_df = pl.read_csv(
            io.StringIO(str(filtered_bt)),
            separator="\t",
            has_header=False,
            new_columns=df.columns,
            schema_overrides=df.schema
        )
        return filtered_df
    except Exception as e:
        logging.error(f"Error during blacklist filtering: {e}")
        raise e

def perform_sample_barcode_correction(raw_reads_df: pl.DataFrame, whitelist_barcodes: Sequence[str], correction_threshold: int) -> pl.DataFrame:
    unique_raw_barcodes = raw_reads_df["sample_barcode_raw"].unique().to_list()
    whitelist_set = set(whitelist_barcodes)
    barcode_correction_map = {}
    
    exact_matches = 0
    corrected_matches = 0
    rejected_matches = 0

    for raw_bc in unique_raw_barcodes:
        if raw_bc is None: continue
        if raw_bc in whitelist_set:
            barcode_correction_map[raw_bc] = raw_bc
            exact_matches += 1
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
            corrected_matches += 1
        else:
            rejected_matches += 1

    logging.info(f"Barcode Correction Stats: {exact_matches} exact, {corrected_matches} corrected (dist<={correction_threshold}), {rejected_matches} rejected.")
    barcode_map_df = pl.DataFrame([{"sample_barcode_raw": k, "sample_barcode": v} for k, v in barcode_correction_map.items() if v is not None])
    if barcode_map_df.is_empty(): return raw_reads_df.clear()
    return raw_reads_df.join(barcode_map_df, on="sample_barcode_raw", how="inner")

def fill_genomic_gaps(data_bins_df: pl.DataFrame, chrom_sizes_map: Dict[str, int], count_col: str) -> pl.DataFrame:
    present_chroms = data_bins_df["chrom"].unique().to_list()
    gap_dfs = []
    for chrom in present_chroms:
        cdf = data_bins_df.filter(pl.col("chrom") == chrom).sort("start")
        if cdf.height == 0: continue
        chrom_len = chrom_sizes_map.get(chrom)
        if not chrom_len: continue
        
        starts = cdf["start"].to_list()
        ends = cdf["end"].to_list()
        gap_starts = [0] + ends
        gap_ends = starts + [chrom_len]
        
        gaps = pl.DataFrame({
            "chrom": chrom,
            "start": gap_starts,
            "end": gap_ends,
            count_col: 0.0
        })
        valid_gaps = gaps.filter(pl.col("end") > pl.col("start"))
        if valid_gaps.height > 0: gap_dfs.append(valid_gaps)
            
    if not gap_dfs: return data_bins_df.sort(["chrom", "start"])
    
    all_gaps = pl.concat(gap_dfs)
    if data_bins_df.schema[count_col] != all_gaps.schema[count_col]:
        all_gaps = all_gaps.with_columns(pl.col(count_col).cast(data_bins_df.schema[count_col]))

    return pl.concat([
        data_bins_df.select(["chrom", "start", "end", count_col]),
        all_gaps.select(["chrom", "start", "end", count_col])
    ]).sort(["chrom", "start"])

def generate_variable_width_bins(df: pl.DataFrame, max_gap: int, count_col: str, chrom_sizes_map: Dict[str, int]) -> pl.DataFrame:
    """
    Merges adjacent insertion sites, tiles the genome with 0s, returns full coverage BED.
    """
    islands_df = (
        df.sort(["chrom", "start"])
        .with_columns([
            (pl.col("start") - pl.col("end").shift(1)).fill_null(9999999).alias("gap"),
            (pl.col("chrom") != pl.col("chrom").shift(1)).fill_null(True).alias("new_chrom")
        ])
        .with_columns(((pl.col("gap") > max_gap) | pl.col("new_chrom")).alias("is_break"))
        .with_columns(pl.col("is_break").cum_sum().alias("cluster_id"))
        .group_by("cluster_id")
        .agg([
            pl.col("chrom").first(),
            pl.col("start").min(),
            pl.col("end").max(),
            pl.col(count_col).sum()
        ])
        .select(["chrom", "start", "end", count_col])
    )
    tiled_df = fill_genomic_gaps(islands_df, chrom_sizes_map, count_col)
    return tiled_df.with_columns((pl.col("end") - pl.col("start")).alias("width"))

def convert_bedgraphs_to_qbed(bedgraph_dir: Path, qbed_dir: Path):
    qbed_dir.mkdir(parents=True, exist_ok=True)
    bedgraph_files = sorted(bedgraph_dir.glob("*.bedgraph"))
    if not bedgraph_files: return
    for bedgraph_file in tqdm(bedgraph_files, desc="Converting to qbed"):
        qbed_path = qbed_dir / (bedgraph_file.stem + ".qbed")
        if qbed_path.exists(): continue
        df = pd.read_csv(bedgraph_file, sep="\t", header=None, names=["chrom", "start", "end", "count"], dtype={"chrom": str, "start": int, "end": int, "count": int})
        if df.empty: continue
        df_expanded = df.loc[df.index.repeat(df["count"])]
        df_expanded[["chrom", "start", "end"]].to_csv(qbed_path, sep="\t", index=False, header=False)

# ==============================================================================
# CORE WORKER FUNCTION
# ==============================================================================

def process_chromosome_chunk(task: tuple, output_dir: Path) -> Optional[pl.DataFrame]:
    parquet_path_str, sample_name, library_name_concat, srt_bc_clustering_threshold, variable_bin_gap = task
    parquet_path = Path(parquet_path_str)

    try:
        chunk_df = pl.read_parquet(parquet_path)
        chunk_df = chunk_df.with_columns(pl.lit(library_name_concat).alias("library_name"))
        
        fragment_deduplication_grouping_keys = [col for col in chunk_df.columns if col not in ["sample_barcode_raw", "reads"]]
        chunk_df = chunk_df.group_by(fragment_deduplication_grouping_keys).agg(
            pl.sum("reads").alias("reads"),
            pl.col("sample_barcode_raw").first().alias("sample_barcode_raw")
        )

        # Proximity Clustering
        clustered_df = (
            chunk_df
            .sort(["start"])
            .with_columns([
                (pl.col("start") - pl.col("end").shift(1)).fill_null(9999999).alias("gap")
            ])
            .with_columns((pl.col("gap") > variable_bin_gap).alias("is_break"))
            .with_columns(pl.col("is_break").cum_sum().alias("variable_bin_id"))
        )

        # UMI Correction
        sites = []
        for _, bin_group_df in clustered_df.group_by("variable_bin_id"):
            srt_bc_counts = bin_group_df.group_by("srt_bc").agg(pl.sum("reads"))
            srt_bc_read_count_dict = {row['srt_bc']: row['reads'] for row in srt_bc_counts.to_dicts() if row['srt_bc']}
            if not srt_bc_read_count_dict: continue

            clusterer = UMIClusterer(cluster_method="directional")
            clustered_srt_bcs = clusterer({k.encode('utf-8'): v for k, v in srt_bc_read_count_dict.items()}, threshold=srt_bc_clustering_threshold)
            
            if not clustered_srt_bcs:
                if len(srt_bc_read_count_dict) > 0:
                     srt_bc_rep_map = {k: k for k in srt_bc_read_count_dict}
                else:
                     continue
            else:
                srt_bc_rep_map = {srt_bc.decode('utf-8'): cluster[0].decode('utf-8') for cluster in clustered_srt_bcs for srt_bc in cluster}
            
            deduped_frags_df = (
                bin_group_df
                .with_columns(pl.col("srt_bc").replace(srt_bc_rep_map).alias("srt_bc_rep"))
                .filter(pl.col("srt_bc_rep").is_not_null())
                .group_by("chrom", "start", "end", "strand", "srt_bc_rep")
                .agg(pl.sum("reads").alias("reads_in_molecule"))
            )

            summary = (
                deduped_frags_df.group_by("srt_bc_rep")
                .agg([
                    pl.col("start").max().alias("max_start"),
                    pl.col("end").min().alias("min_end"),
                    pl.first("strand").alias("strand"),
                    pl.first("chrom").alias("chrom"),
                    pl.sum("reads_in_molecule").alias("reads_for_fragment")
                ])
                .with_columns(
                    pl.when(pl.col("strand") == "+").then(pl.col("min_end"))
                    .otherwise(pl.col("max_start")).alias("site")
                )
                .select(["chrom", "site", "strand", "reads_for_fragment"])
            )
            sites.append(summary)

        if not sites: return None
        final_sites_df = pl.concat(sites)
        
        return (
            final_sites_df.group_by(["chrom", "site", "strand"])
            .agg([pl.count().alias("unique_count"), pl.sum("reads_for_fragment").alias("total_reads")])
            .rename({"site": "start"})
            .with_columns((pl.col("start") + 1).alias("end"))
            .select(["chrom", "start", "end", "unique_count", "total_reads", "strand"])
            .sort(["start"])
        )

    except Exception as exc:
        logging.error(f"Error processing {parquet_path.name}: {exc}", exc_info=True)
        return None

# ==============================================================================
# MAIN PIPELINE
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--input_dir", required=True, help="Directory containing input qbed/bed files.")
    parser.add_argument("--output_dir", required=True, help="Directory to store all output files.")
    parser.add_argument("--annotation_file", required=True, type=str, help="Path to annotation file.")
    parser.add_argument("--chrom_sizes", required=True, help="Path to a chromosome sizes file.")
    parser.add_argument("--blacklist_file", type=str, help="Path to a blacklist BED file. Regions (+/- 150bp) will be excluded.")
    parser.add_argument("--workers", type=int, default=12, help="Number of parallel workers.")
    parser.add_argument("--srt_bc_dist_threshold", type=int, default=2, help="Max Hamming distance for SRT-BC (UMI) clustering.")
    parser.add_argument("--sample_barcode_dist_threshold", type=int, default=3, help="Max Hamming distance for sample barcode correction.")
    parser.add_argument("--min_rows_threshold", type=int, default=0, help="Min fragments for a sample to be processed.")
    parser.add_argument("--sum_window_size", type=int, default=100, help="Window size in bp for binned summary BigWig tracks (Default 100).")
    parser.add_argument("--variable_bin_gap", type=int, default=1000, help="Max gap to merge adjacent insertions for variable-width binning.")

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    setup_logging(output_dir)
    print_version_info()
    validate_args(args)
    
    pybedtools_temp_dir = output_dir / "pybedtools_tmp"
    pybedtools_temp_dir.mkdir(exist_ok=True, parents=True)
    pybedtools.helpers.set_tempdir(pybedtools_temp_dir)

    chrom_sizes_map = load_chrom_sizes(Path(args.chrom_sizes))

    # --- Phase 1: Load, Filter (Blacklist), and Partition ---
    logging.info("--- Phase 1: Loading, Correcting, Annotating, and Partitioning Reads ---")
    annotation_df = pl.read_csv(args.annotation_file, separator='\t' if args.annotation_file.endswith(('.tsv', '.txt')) else ',', has_header=True)
    whitelist_sample_barcodes = annotation_df["sample_barcode"].unique().to_list()
    input_files = [p for ext in ("*.qbed", "*.bed", "*.qbed.gz", "*.bed.gz") for p in Path(args.input_dir).glob(ext)]
    
    # Prepare Blacklist if provided
    expanded_blacklist_bt = None
    if args.blacklist_file:
        logging.info(f"Loading blacklist file: {args.blacklist_file}")
        try:
            # Load blacklist, slop by 150bp using chromosome sizes
            bl_bt = pybedtools.BedTool(args.blacklist_file)
            expanded_blacklist_bt = bl_bt.slop(b=150, g=str(args.chrom_sizes)).merge()
            logging.info("Blacklist loaded and expanded by +/- 150bp.")
        except Exception as e:
            logging.error(f"Failed to process blacklist file: {e}")
            sys.exit(1)

    partitioned_temp_dir = output_dir / "partitioned_temp"
    partitioned_temp_dir.mkdir(exist_ok=True, parents=True)
    observed_sample_library_pairs = set()

    for file_path in tqdm(input_files, desc="Reading Input Files"):
        try:
            fragment_df = pl.read_csv(file_path, has_header=False, separator="\t", new_columns=list(INPUT_QBED_SCHEMA.keys()), schema_overrides=INPUT_QBED_SCHEMA, ignore_errors=True)
            
            initial_count = fragment_df.height
            fragment_df = normalize_and_filter_chromosomes(fragment_df, VALID_CHROMOSOMES)
            
            # --- APPLY BLACKLIST FILTERING HERE ---
            if expanded_blacklist_bt is not None and not fragment_df.is_empty():
                count_before_blacklist = fragment_df.height
                fragment_df = filter_dataframe_by_blacklist(fragment_df, expanded_blacklist_bt)
                count_after_blacklist = fragment_df.height
                dropped_reads = count_before_blacklist - count_after_blacklist
                if dropped_reads > 0:
                     # Optional: Debug log, reduced to avoid spamming console
                     pass

            fragment_df = (
                fragment_df.with_columns([
                    pl.col("name").str.split("/").list.get(0).alias("library_name"),
                    pl.col("name").str.split("/").list.get(1).alias("sample_barcode_raw"),
                    pl.col("name").str.split("/").list.get(2).alias("srt_bc")
                ])
                .drop("name")
                .filter(pl.col("sample_barcode_raw").is_not_null())
            )
            
            if fragment_df.is_empty():
                logging.warning(f"File {file_path.name}: 0 reads survived filtering (chrom + blacklist).")
                continue

            corrected_reads_df = perform_sample_barcode_correction(fragment_df, whitelist_sample_barcodes, args.sample_barcode_dist_threshold)
            if corrected_reads_df.is_empty():
                logging.warning(f"File {file_path.name}: No barcodes matched whitelist.")
                continue

            annotated_reads_df = corrected_reads_df.join(annotation_df, on=["library_name", "sample_barcode"], how="inner")
            if annotated_reads_df.is_empty():
                logging.warning(f"File {file_path.name}: Annotation join failed.")
                continue

            for row in annotated_reads_df.select(["sample_name", "library_name"]).unique().iter_rows():
                observed_sample_library_pairs.add(row)
            
            for (sample_name,), group_df in annotated_reads_df.group_by("sample_name"):
                group_df.write_parquet(partitioned_temp_dir / f"{sample_name}_{file_path.stem}.parquet")
        except Exception as exc:
            logging.error(f"Failed to process {file_path}: {exc}")

    # --- Phase 2: Consolidate by Chromosome ---
    logging.info("--- Phase 2: Consolidating Samples by Chromosome ---")
    collapsed_dir = output_dir / "collapsed_fragments_per_sample_chrom"
    collapsed_dir.mkdir(exist_ok=True)
    
    sample_to_libs_map = {}
    for s, l in observed_sample_library_pairs: sample_to_libs_map.setdefault(s, set()).add(l)
    
    processing_tasks = []
    unique_samples = sorted(sample_to_libs_map.keys())
    for sample_name in tqdm(unique_samples, desc="Consolidating & Splitting"):
        sample_files = list(partitioned_temp_dir.glob(f"{sample_name}_*.parquet"))
        if not sample_files: continue
        full_sample_df = pl.read_parquet(sample_files)
        if full_sample_df.height < args.min_rows_threshold:
            logging.warning(f"Skipping sample '{sample_name}': {full_sample_df.height} reads.")
            continue
        lib_str = "__".join(sorted(sample_to_libs_map[sample_name]))
        for (chrom,), chrom_df in full_sample_df.group_by("chrom"):
            chunk_path = collapsed_dir / f"{sample_name}__{chrom}.parquet"
            chrom_df.write_parquet(chunk_path)
            processing_tasks.append((str(chunk_path), sample_name, lib_str, args.srt_bc_dist_threshold, args.variable_bin_gap))

    shutil.rmtree(partitioned_temp_dir)
    if not processing_tasks:
        logging.error("No samples passed criteria. Exiting.")
        sys.exit(1)

    # --- Phase 3: Run Parallel Processing ---
    logging.info(f"Processing {len(processing_tasks)} chromosome chunks with {args.workers} workers...")
    sample_results_accumulator = {}
    
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(process_chromosome_chunk, t, output_dir): t[1] for t in processing_tasks}
        for f in tqdm(as_completed(futures), total=len(futures), desc="Mapping Insertions"):
            res = f.result()
            sample_name = futures[f]
            if res is not None and not res.is_empty():
                sample_results_accumulator.setdefault(sample_name, []).append(res)

    # --- Phase 4: Aggregation & Final Outputs ---
    logging.info("--- Phase 4: Generating Final Outputs ---")
    dirs = {
        "sample_bg": output_dir / "raw_per_sample_bedgraph",
        "sample_bw": output_dir / "raw_per_sample_bigwig",
        "sample_cpm_bg": output_dir / "per_sample_bedgraph_cpm",
        "sample_cpm_bw": output_dir / "per_sample_bigwig_cpm",
        "sample_var_bw": output_dir / f"per_sample_variable_width_bins_{args.variable_bin_gap}bp",
        "sample_stats": output_dir / "per_sample_qbed_with_stats",
        "sample_qbed": output_dir / "raw_per_sample_qbed",
        "group_bg": output_dir / "raw_per_group_bedgraph",
        "group_bw": output_dir / "raw_per_group_bigwig",
        "group_cpm_bg": output_dir / "per_group_bedgraph_cpm",
        "group_cpm_bw": output_dir / "per_group_bigwig_cpm",
        "group_stats": output_dir / "per_group_qbed_with_stats",
        "qbed": output_dir / "raw_per_group_qbed",
    }
    for d in dirs.values(): d.mkdir(exist_ok=True, parents=True)

    chrom_sizes_path = Path(args.chrom_sizes)
    cleaned_ann_df = pl.read_csv(args.annotation_file, separator='\t' if args.annotation_file.endswith(('.tsv', '.txt')) else ',', has_header=True).rename({c: c.strip() for c in annotation_df.columns}).with_columns(pl.col("sample_name").str.strip_chars())
    group_to_samples = cleaned_ann_df.group_by("group_name").agg(pl.col("sample_name")).rows()

    # 4a. Combine Chromosomes back into Samples
    sample_maps = {}
    for sample_name, chrom_dfs in tqdm(sample_results_accumulator.items(), desc="Merging Chromosomes"):
        if not chrom_dfs: continue
        full_sample_df = pl.concat(chrom_dfs).sort(["chrom", "start"])
        # Clip Coordinates
        full_sample_df = clip_coordinates(full_sample_df, chrom_sizes_map)
        sample_maps[sample_name] = full_sample_df
        
        # 1. Raw: 1-bp BedGraph, Fixed BigWig
        bg_path, bw_path = dirs["sample_bg"] / f"{sample_name}.bedgraph", dirs["sample_bw"] / f"{sample_name}.bw"
        safe_write_1bp_bedgraph(full_sample_df, bg_path, score_col="unique_count")
        
        # Pass sum_window_size (default 100) here
        bin_and_write_bigwig(full_sample_df, bw_path, chrom_sizes_path, chrom_sizes_map, score_col="unique_count", bin_size=args.sum_window_size)
        
        # 2. CPM: 1-bp BedGraph, Fixed BigWig
        total = full_sample_df["unique_count"].sum()
        if total > 0:
            cpm_df = full_sample_df.with_columns((pl.col("unique_count") * 1e6 / total).alias("cpm")).select(["chrom", "start", "end", "cpm"])
            cpm_bg_path, cpm_bw_path = dirs["sample_cpm_bg"] / f"{sample_name}.bedgraph", dirs["sample_cpm_bw"] / f"{sample_name}.bw"
            safe_write_1bp_bedgraph(cpm_df, cpm_bg_path, score_col="cpm")
            
            # Pass sum_window_size (default 100) here
            bin_and_write_bigwig(cpm_df, cpm_bw_path, chrom_sizes_path, chrom_sizes_map, score_col="cpm", bin_size=args.sum_window_size)

            # 3. Variable Width Bins (BED ONLY, Tiled 0-ChromSize)
            var_bin_df = generate_variable_width_bins(full_sample_df, args.variable_bin_gap, "unique_count", chrom_sizes_map)
            var_bed_path = dirs["sample_var_bw"] / f"{sample_name}.bed"
            var_bin_df.write_csv(var_bed_path, separator="\t", include_header=False)

        # 4. Stats
        full_sample_df.select(["chrom", "start", "end", "unique_count", "total_reads", "strand"]).write_csv(dirs["sample_stats"] / f"{sample_name}.qbed", separator="\t", include_header=False)

    # 4b. Groups
    for group_name, sample_list in tqdm(group_to_samples, desc="Aggregating Groups"):
        valid_samples = [sample_maps[s] for s in sample_list if s in sample_maps]
        if not valid_samples: continue
        agg_df = pl.concat(valid_samples).group_by(["chrom", "start", "end", "strand"]).agg([pl.sum("unique_count"), pl.sum("total_reads")]).sort(["chrom", "start"])
        agg_df = clip_coordinates(agg_df, chrom_sizes_map)

        # 1. Raw: 1-bp BedGraph, Fixed BigWig
        bg_path, bw_path = dirs["group_bg"] / f"{group_name}.bedgraph", dirs["group_bw"] / f"{group_name}.bw"
        safe_write_1bp_bedgraph(agg_df, bg_path, score_col="unique_count")
        
        # Pass sum_window_size (default 100) here
        bin_and_write_bigwig(agg_df, bw_path, chrom_sizes_path, chrom_sizes_map, score_col="unique_count", bin_size=args.sum_window_size)

        # 2. CPM: 1-bp BedGraph, Fixed BigWig
        total = agg_df["unique_count"].sum()
        if total > 0:
            cpm_df = agg_df.with_columns((pl.col("unique_count") * 1e6 / total).alias("cpm")).select(["chrom", "start", "end", "cpm"])
            cpm_bg_path, cpm_bw_path = dirs["group_cpm_bg"] / f"{group_name}.bedgraph", dirs["group_cpm_bw"] / f"{group_name}.bw"
            safe_write_1bp_bedgraph(cpm_df, cpm_bg_path, score_col="cpm")
            
            # Pass sum_window_size (default 100) here
            bin_and_write_bigwig(cpm_df, cpm_bw_path, chrom_sizes_path, chrom_sizes_map, score_col="cpm", bin_size=args.sum_window_size)

        # 3. Stats
        agg_df.select(["chrom", "start", "end", "unique_count", "total_reads", "strand"]).write_csv(dirs["group_stats"] / f"{group_name}.qbed", separator="\t", include_header=False)

    convert_bedgraphs_to_qbed(dirs["group_bg"], dirs["qbed"])
    convert_bedgraphs_to_qbed(dirs["sample_bg"], dirs["sample_qbed"])
    
    logging.info("Pipeline complete")
    pybedtools.helpers.cleanup()

if __name__ == "__main__":
    main()
