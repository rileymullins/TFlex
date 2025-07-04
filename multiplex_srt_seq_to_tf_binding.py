import os
import sys
import argparse
import logging
import io
import glob
import shutil
from pathlib import Path
from contextlib import redirect_stdout, redirect_stderr
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Sequence, Optional
import polars as pl
import pandas as pd
import pybedtools
from tqdm import tqdm
from umi_tools import UMIClusterer
import pycallingcards as cc

# ==============================================================================
# SCRIPT CONFIGURATION
# ==============================================================================

VERSION = "1.0.0"

# Define the set of valid chromosome names to retain for downstream processing
VALID_CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

# Default parameters for pycallingcards peak calling, controlling stringency and peak detection thresholds
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
    Logs versions of the script and dependencies.
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
    logging.info("All arguments are valid.")

def write_column_definitions(output_dir: Path) -> None:
    """
    Writes a table explaining all output columns in the final results file.
    Pipeline step: Output aggregation (column_definitions.tsv).
    """
    column_definitions = """Variable    Description
consensus_peak_id    A unique identifier for the final merged peak, formatted as chrom:start-end.
chrom    The chromosome of the final merged peak.
consensus_peak_start    The start coordinate of the final merged, pan-dataset consensus peak.
consensus_peak_end    The end coordinate of the final merged, pan-dataset consensus peak.
consensus_peak_width    The width (end - start) of the final merged, pan-dataset consensus peak.
num_samples_in_consensus_peak    The number of unique sample_name values that contributed reads to this consensus peak.
num_groups_in_consensus_peak    The number of unique group_name values that contributed reads to this consensus peak.
sample_name    The sample identifier, as defined in the annotation file.
group_name    The experimental group identifier, as defined in the annotation file.
library_name    Double-underscore concatenation of all library identifiers for this sample_name.
fragment_peak_start_for_sample    The start coordinate of the initial raw peak called on the data for this sample_name.
fragment_peak_end_for_sample    The end coordinate of the initial raw peak for this sample_name.
fragment_peak_width_for_sample    The width of the initial raw peak for this sample_name.
srt_bc_peak_start_for_sample    The start coordinate of the peak for this sample_name after SRT-BC refinement.
srt_bc_peak_end_for_sample    The end coordinate of the SRT-BC-refined peak for this sample_name.
srt_bc_peak_width_for_sample    The width of the SRT-BC-refined peak for this sample_name.
sample_total_reads_in_consensus_peak    The total number of reads from this sample_name that fall within the final consensus peak.
sample_total_fragments_in_consensus_peak    The number of unique molecular fragments (qbed rows after SRT barcode correction and coordinate deduplication) from this sample_name within the final consensus peak.
sample_total_insertions_in_consensus_peak    The number of unique SRT barcodes from this sample_name within the final consensus peak.
group_total_reads_in_consensus_peak    The sum of sample_total_reads_in_consensus_peak for all sample_names belonging to this group_name within the consensus peak.
group_total_fragments_in_consensus_peak    The sum of sample_total_fragments_in_consensus_peak for all sample_names belonging to this group_name within the consensus peak.
group_total_insertions_in_consensus_peak    The sum of sample_total_insertions_in_consensus_peak for all sample_names belonging to this group_name within the consensus peak.
"""
    definitions_path = output_dir / "column_definitions.tsv"
    with open(definitions_path, "w") as f:
        f.write(column_definitions.strip())
    logging.info(f"Wrote column definitions to {definitions_path}")

def perform_sample_barcode_correction(
    raw_reads_df: pl.DataFrame,
    whitelist_barcodes: Sequence[str],
    correction_threshold: int
) -> pl.DataFrame:
    """
    Corrects raw sample barcodes to those in the annotation whitelist (all unique sample_barcode values act as the whitelist) using Hamming distance.
    Only barcodes within the allowed threshold are retained (barcode correction step).
    Returns a DataFrame with corrected sample barcodes.
    """
    # ====== MAJOR STEP: BARCODE CORRECTION LOGIC ======
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
        # Calculate Hamming distance to each whitelist barcode
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

# ==============================================================================
# CORE WORKER FUNCTION
# ==============================================================================

def process_pooled_sample_data(task: tuple, annotation_df: Optional[pl.DataFrame]=None, output_dir: Optional[Path]=None) -> Optional[pl.DataFrame]:
    """
    Processes pooled reads for a single sample_name (across all libraries for that sample).
    - Calls fragment-based peaks using pycallingcards (removes noise by requiring at least 5 unique fragments per peak).
    - Performs SRT barcode-based peak refinement (strand-aware, ±100bp), assigns peak stats, and saves results as parquet.
    - Counts number of unique SRT barcodes in each peak (proxy for unique transposon insertions).
    - Returns a DataFrame of detailed per-sample peak statistics.
    Pipeline steps: Fragment-based peak calling, SRT barcode-based peak refinement.
    """
    sample_parquet_path_str, sample_name, library_name_concat, peak_calling_params, srt_bc_clustering_threshold = task
    sample_parquet_path = Path(sample_parquet_path_str)

    try:
        # ====== MAJOR STEP: 5a. LOAD AND PREPARE SAMPLE DATA ======
        # Load the pooled fragment data for this sample_name from the .parquet file
        sample_fragments_df = pl.read_parquet(sample_parquet_path)
        
        # Set library_name to double-underscore-concatenated string for this sample_name
        sample_fragments_df = sample_fragments_df.with_columns(pl.lit(library_name_concat).alias("library_name"))

        # Get group_name for this sample_name from the annotation
        group_name = None
        if annotation_df is not None:
            row = annotation_df.filter(pl.col("sample_name") == sample_name)
            group_name = row["group_name"][0] if row.height > 0 else None

        # ====== MAJOR STEP: 5b. FRAGMENT-LEVEL DEDUPLICATION ======
        # Deduplicate rows/fragments that are unique only due to a sequencing error in the same barcode (i.e., only unique in 'sample_barcode_raw' value).
        # This deduplicated dataframe is used for intersecting fragments with fragment-based peaks too.
        logging.info(f"[{sample_name}] Performing fragment-level deduplication: merging rows identical in all fields except only 'sample_barcode_raw' and summing 'reads'.")

        # Define the columns to group by. These are all columns present in sample_fragments_df except 'sample_barcode_raw' AND 'reads'.
        fragment_deduplication_grouping_keys = [
            col for col in sample_fragments_df.columns if col not in ["sample_barcode_raw", "reads"]
        ]

        # Perform the grouping and sum the 'reads' of these identical rows.
        # The first 'sample_barcode_raw' encountered in the merge is retained. However, it's not used downstream.
        initial_fragment_count = sample_fragments_df.height
        sample_fragments_df = sample_fragments_df.group_by(fragment_deduplication_grouping_keys).agg(
            pl.sum("reads").alias("reads"),
            pl.col("sample_barcode_raw").first().alias("sample_barcode_raw")
        )
        logging.info(f"[{sample_name}] Fragments before this deduplication: {initial_fragment_count}")
        logging.info(f"[{sample_name}] Fragments after this deduplication: {sample_fragments_df.height}")



        # ====== MAJOR STEP: 6. FRAGMENT-BASED PEAK CALLING ======
        # Run peak caller on all fragments/qbed rows (ignoring SRT barcodes at this stage)
        pycc_fragments_df = (
            sample_fragments_df
            .select(["chrom", "start", "end"])
            .to_pandas()
            .rename(columns={"chrom": "Chr", "start": "Start", "end": "End"})
        )
        # Silence stdout/stderr from pycallingcards for clean logs
        with redirect_stderr(io.StringIO()), redirect_stdout(io.StringIO()):
            initial_peaks_pd = cc.pp.call_peaks(expdata=pycc_fragments_df, **peak_calling_params)

        if initial_peaks_pd is None or initial_peaks_pd.empty:
            return None

        # Format peak coordinates and assign unique sample_peak_id
        initial_peaks_df = pl.from_pandas(initial_peaks_pd.iloc[:, :3])
        initial_peaks_df.columns = ["chrom", "fragment_peak_start_for_sample", "fragment_peak_end_for_sample"]
        initial_peaks_df = initial_peaks_df.with_columns(
            (pl.lit(f"{sample_name}_peak_") + pl.arange(0, initial_peaks_df.height).cast(pl.Utf8)).alias("sample_peak_id")
        )

        # Save these initial fragment peaks.
        if output_dir is not None:
            fragment_peaks_dir = output_dir / "fragment_peaks"
            fragment_peaks_dir.mkdir(exist_ok=True, parents=True)
            fragment_peaks_path = fragment_peaks_dir / f"{sample_name}_fragment_peaks.parquet"
            initial_peaks_df.write_parquet(fragment_peaks_path)


        # ====== MAJOR STEP: 7. MAP FRAGMENTS TO PEAKS ======
        # Determines which fragments are in each called peak (to be refined by SRT barcodes)
        # Intersect fragments with fragment-based peaks to map fragments to fragment-based peaks
        bed_fragments_df = sample_fragments_df.select(
            ["chrom", "start", "end", "strand", "srt_bc", "reads"]
        ).to_pandas()
        bed_fragments = pybedtools.BedTool.from_dataframe(bed_fragments_df)
        bed_peaks = pybedtools.BedTool.from_dataframe(initial_peaks_df.to_pandas())

        fragment_columns = bed_fragments_df.columns.tolist()
        peak_columns = [f"{col}_peak" for col in initial_peaks_df.columns]
        intersection_columns = fragment_columns + peak_columns

        # Each row now represents a fragment assigned to a specific fragment-based peak
        intersected_df = pl.from_pandas(
            bed_fragments.intersect(bed_peaks, wa=True, wb=True).to_dataframe(names=intersection_columns)
        )
        if intersected_df.is_empty():
            return None



        # ====== MAJOR STEP: 8a. REFINE PEAKS WITH SRT BARCODES ======
        # SRT barcode-based peak refinement and per-sample stats
        # Loop over all fragment-based peaks, refining boundaries by SRT barcode logic
        refined_peak_boundaries = []
        for peak_id, peak_group_df in intersected_df.group_by("sample_peak_id_peak"):
            # Collapse reads by SRT barcode (deduplication per SRT barcode)
            srt_bc_counts = peak_group_df.group_by("srt_bc").agg(pl.sum("reads"))
            srt_bc_read_count_dict = {
                row['srt_bc']: row['reads']
                for row in srt_bc_counts.to_dicts()
                if row['srt_bc'] is not None
            }
            if not srt_bc_read_count_dict:
                continue

            # Cluster SRT barcodes by Hamming distance (UMI correction)
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





            # ============================================================ #
            # CRITICAL LOGIC OF DEFINING SRT BARCODE-BASED PEAK BOUNDARIES #
            # ============================================================ #

            
            # Make a dataframe of '+' strand fragments that has the min 'end' value across all fragments of each SRT barcode listed.
                # The 'end' value for '+' strand fragments is the true end of the read as it came off the sequencer.
            # The reads of '+' strand fragments move 5' <- 3' (opposite of convention).
                # The most 5' (upstream, smallest) value is the the most proximal to the true transposon/genomic DNA junction.
                
            plus_strand_frags = deduped_frags_df.filter(pl.col("strand") == "+")
            plus_strand_coords = (
                plus_strand_frags.group_by("srt_bc_rep")
                .agg(pl.min("end"))["end"]  # The "end" value for '+' strand coords is the true end of the read as it came off the sequencer. This selects the smallest (most upstream)'end' value per SRT barcode.
                .to_list()
            )
            
            
            
            # Make a dataframe of '-' strand fragments that has the max start value of all fragments of each SRT barcode listed.
                # The 'start' value for '-' strand coords is the true end of the read as it came off the sequencer.
            # The reads of '-' strand fragments move 5' -> 3' (opposite of convention).
                # The most 3' (downstream, largest) value is the the most proximal to the true transposon/genomic DNA junction.
                
            minus_strand_frags = deduped_frags_df.filter(pl.col("strand") == "-")
            minus_strand_coords = (
                minus_strand_frags.group_by("srt_bc_rep")
                .agg(pl.max("start"))["start"]# The "start" value for '-' strand coords is the true end of the read as it came off the sequencer. This selects the largest (most downstream) 'start' value per SRT barcode.
                .to_list()
            )



            # Combine all coordinates into one to make the peak boundary definitions simpler
            all_coords = plus_strand_coords + minus_strand_coords
            
            
            # Define SRT barcode-based peak boundaries
            if all_coords:
            
                # min of all_coords can be used because the strand-aware selection of the most proximal coordinate per SRT barcode has already been selected. Therefore, the min is the most upstream (smallst) position of the fragment that is the closest to the transpsoson per SRT barcode; it is the best approximation of where the most 5' unique insertion occurred. Thus, it's the best approximation of the SRT barcode-bound start of the peak.
                min_coord = min(all_coords)
                
                
                # Same rationale as for the min(all_coords), but now max is defined as the most downstream (largest) position; it is the best approximation of where the most 3' unique insertion occurred. Thus, it's the best approximation of the SRT barcode-bound end of the peak.
                max_coord = max(all_coords)
                
                
                # The SRT barcode-based peak has a hard coded 100 bp extension to the min_coord (-100) and max_coord (+100). This is to build in a margin of error when defining the SRT barcode-based peak.
                # Generation of the consensus peaks through merging overlapping peaks is a second margin of error layer used.
            
                srt_bc_peak_start = max(0, min_coord - 100) # This prevents the the start coordinate being below 0.
                srt_bc_peak_end = max_coord + 100
                
            else:
                # This case should not be reached if deduped_frags_df is not empty, but included anyways.
                continue
            
            # ============================================================ #
            # END SECTION ON DEFINING SRT BARCODE-BASED PEAK BOUNDARIES    #
            # ============================================================ #




            # ====== MAJOR STEP: 8b. COMPUTE PER-SAMPLE PEAK STATISTICS ======
            # Compute all per-sample stats for this peak
            fragment_peak_start = peak_group_df["fragment_peak_start_for_sample_peak"][0]
            fragment_peak_end = peak_group_df["fragment_peak_end_for_sample_peak"][0]
            fragment_peak_width = fragment_peak_end - fragment_peak_start
            srt_bc_peak_width = srt_bc_peak_end - srt_bc_peak_start

            sample_total_reads = int(peak_group_df["reads"].sum())
            sample_total_fragments = deduped_frags_df.height
            sample_total_insertions = deduped_frags_df["srt_bc_rep"].n_unique()

            refined_peak_boundaries.append({
                "chrom": peak_group_df["chrom"][0],
                "sample_name": sample_name,
                "group_name": group_name,
                "library_name": library_name_concat,
                "fragment_peak_start_for_sample": fragment_peak_start,
                "fragment_peak_end_for_sample": fragment_peak_end,
                "fragment_peak_width_for_sample": fragment_peak_width,
                "srt_bc_peak_start_for_sample": srt_bc_peak_start,
                "srt_bc_peak_end_for_sample": srt_bc_peak_end,
                "srt_bc_peak_width_for_sample": srt_bc_peak_width,
                "sample_total_reads_in_consensus_peak": sample_total_reads,
                "sample_total_fragments_in_consensus_peak": sample_total_fragments,
                "sample_total_insertions_in_consensus_peak": sample_total_insertions,
                "sample_peak_id": str(peak_id) if not isinstance(peak_id, tuple) else str(peak_id[0])
            })

        if not refined_peak_boundaries:
            return None

        sample_peak_df = pl.DataFrame(refined_peak_boundaries)

        # Save per-sample peak stats: ~/sample_peaks/sample_name_peaks.parquet
        if output_dir is not None:
            out_path = output_dir / "sample_peaks" / f"{sample_name}_peaks.parquet"
            out_path.parent.mkdir(parents=True, exist_ok=True)
            sample_peak_df.write_parquet(out_path)

        return sample_peak_df

    except Exception as exc:
        logging.error(f"Error processing {sample_parquet_path.name}: {exc}", exc_info=True)
        return None

# ==============================================================================
# MAIN PIPELINE
# ==============================================================================

def main():
    """
    
    Major steps:
    1. Input loading
    2. Barcode correction
    3. Annotation
    4. Per-sample partitioning
    5. Fragment-based peak calling
    6. Fragments-to-fragment-based peak mapping (intersection)
    7. SRT barcode-based peak refinement
    8. Consensus peak generation
    9. Sample peak set-to-consensus peak mapping (intersection)
    10. Output aggregation
    
    """
    parser = argparse.ArgumentParser(description="An annotation-driven calling card data processing pipeline.")
    parser.add_argument("--input_dir", required=True, help="Directory containing input qbed/bed files.")
    parser.add_argument("--output_dir", required=True, help="Directory to store output files.")
    parser.add_argument("--workers", type=int, default=10, help="Number of parallel workers.")
    parser.add_argument("--srt_bc_dist_threshold", type=int, default=1, help="Max Hamming distance for SRT-BC (UMI) clustering.")
    parser.add_argument("--sample_barcode_dist_threshold", type=int, default=2, help="Max Hamming distance for sample barcode correction.")
    parser.add_argument("--min_rows_threshold", type=int, default=50000, help="Min fragments for a sample to be processed.")
    parser.add_argument("--annotation_file", required=True, type=str, help="Path to a 4-column TSV or CSV file with header: 'library_name,sample_barcode,sample_name,group_name'.")
    for param_key, param_default in PEAK_CALLING_PARAMS_DEFAULTS.items():
        parser.add_argument(f"--{param_key}", type=type(param_default), default=param_default)
    args = parser.parse_args()




    # ==== Set up logging and output directories ====
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    setup_logging(output_dir)
    print_version_info()
    validate_args(args)
    write_column_definitions(output_dir)




    # ==== Set a dedicated temporary directory for pybedtools ====
    pybedtools_temp_dir = output_dir / "pybedtools_tmp"
    pybedtools_temp_dir.mkdir(exist_ok=True)
    pybedtools.helpers.set_tempdir(pybedtools_temp_dir)
    
    
    
    logging.info("--- Phase 1: Loading, Correcting, Annotating, and Partitioning Reads (Memory-Efficiently) ---")
    
    # ==== Set a dedicated temporary directory for pybedtools ====
    # Load annotation table and create whitelist
    annotation_file_separator = '\t' if args.annotation_file.endswith(('.tsv', '.txt')) else ','
    annotation_df = pl.read_csv(args.annotation_file, separator=annotation_file_separator, has_header=True)
    whitelist_sample_barcodes = annotation_df["sample_barcode"].unique().to_list()
    logging.info(f"Derived a whitelist of {len(whitelist_sample_barcodes)} barcodes from the annotation file.")

    # Find input files
    input_files = [
        file_path for ext in ("*.qbed", "*.bed", "*.qbed.gz", "*.bed.gz")
        for file_path in Path(args.input_dir).glob(ext)
    ]
    if not input_files:
        logging.error(f"No input files found directly in {args.input_dir}.")
        sys.exit(1)

    # Setup a temporary directory for intermediate partitioned files
    partitioned_temp_dir = output_dir / "partitioned_temp"
    partitioned_temp_dir.mkdir(exist_ok=True, parents=True)
    
    # This will store all (sample_name, library_name) pairs found in the data
    observed_sample_library_pairs = set()

    # Process each input file individually to create intermediate partitioned files
    logging.info(f"Processing {len(input_files)} input files and writing intermediate data...")
    for file_path in tqdm(input_files, desc="Reading Input Files and Partitioning by Sample"):
        try:
            # ====== STEP 1: INPUT LOADING ======
            # Load one file
            fragment_df = pl.read_csv(
                file_path, has_header=False, separator="\t",
                new_columns=list(INPUT_QBED_SCHEMA.keys()),
                schema_overrides=INPUT_QBED_SCHEMA, ignore_errors=True
            )
            
            # Parse name field and filter
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

            if fragment_df.is_empty():
                continue
                
            # ====== STEP 2: BARCODE CORRECTION ======
            # Barcode correction (on the single file's data)
            corrected_reads_df = perform_sample_barcode_correction(
                fragment_df,
                whitelist_barcodes=whitelist_sample_barcodes,
                correction_threshold=args.sample_barcode_dist_threshold
            )
            
            if corrected_reads_df.is_empty():
                continue

            # ====== STEP 3: ANNOTATION ======
            # Annotation (on the single file's data)
            annotated_reads_df = corrected_reads_df.join(
                annotation_df, on=["library_name", "sample_barcode"], how="inner"
            )
            
            if annotated_reads_df.is_empty():
                continue

            # Collect the observed (sample_name, library_name) pairs from this file's data
            pairs_in_file = annotated_reads_df.select(["sample_name", "library_name"]).unique()
            for row in pairs_in_file.iter_rows():
                observed_sample_library_pairs.add(row)

            # ====== STEP 4: PER-SAMPLE PARTITIONING ======
            # Partition this file's data by sample_name and write to temp directory
            for sample_name_tuple, group_df in annotated_reads_df.group_by("sample_name"):
                sample_name = sample_name_tuple[0]
                temp_out_path = partitioned_temp_dir / f"{sample_name}_{file_path.stem}.parquet"
                group_df.write_parquet(temp_out_path)

        except Exception as exc:
            logging.error(f"Failed to process file {file_path}: {exc}", exc_info=True)

    # Consolidate the intermediate files for each sample
    logging.info("--- Phase 2: Consolidating Sample Data and Calling Peaks ---")
    collapsed_dir = output_dir / "collapsed_per_sample"
    collapsed_dir.mkdir(exist_ok=True, parents=True)

    intermediate_files = list(partitioned_temp_dir.glob("*.parquet"))
    if not intermediate_files:
        logging.error("No valid reads could be processed from input files.")
        sys.exit(1)
        
    # Build the library string map from the pairs observed in the actual data
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
    logging.info(f"Found data for {len(unique_sample_names)} samples. Consolidating files.")

    # Consolidate partitioned files into one file per sample
    for sample_name in tqdm(unique_sample_names, desc="Consolidating Samples"):
        sample_files = list(partitioned_temp_dir.glob(f"{sample_name}_*.parquet"))
        if not sample_files:
            continue
            
        sample_df = pl.read_parquet(sample_files)
        
        out_path = collapsed_dir / f"{sample_name}.parquet"
        if sample_df.height >= args.min_rows_threshold:
            sample_df.write_parquet(out_path)
        else:
            logging.warning(f"Skipping sample '{sample_name}': has {sample_df.height} reads, which is less than the threshold of {args.min_rows_threshold}.")

    # Clean up the temporary directory
    shutil.rmtree(partitioned_temp_dir)

    collapsed_files = list(collapsed_dir.glob("*.parquet"))
    if not collapsed_files:
        logging.error(f"No samples passed the min_rows_threshold. Exiting.")
        sys.exit(1)





    # ====== 5. Fragment-based peak calling, 6. Fragments-to-fragment-based peak mapping (intersection), and 7. SRT barcode-based peak refinement ======
    # See def process_pooled_sample_data for how 5, 6, and 7 are coded.
    
    # For each sample, call peaks and perform SRT barcode-based refinement in parallel
    peak_params = {k: getattr(args, k) for k in PEAK_CALLING_PARAMS_DEFAULTS}
    tasks = [
        (str(f), f.stem, sample_name_to_librarystr[f.stem], peak_params, args.srt_bc_dist_threshold)
        for f in collapsed_files
    ]

    sample_peaks_list = []
    sample_peaks_dir = output_dir / "sample_peaks"
    sample_peaks_dir.mkdir(exist_ok=True, parents=True)

    # Check if sample peak files already exist (to allow skipping this phase)
    expected_sample_peak_files = [
        sample_peaks_dir / f"{Path(f).stem}_peaks.parquet"
        for f in collapsed_files
    ]
    required_columns = {
        "sample_name", "group_name", "library_name",
        "fragment_peak_start_for_sample", "fragment_peak_end_for_sample", "fragment_peak_width_for_sample",
        "srt_bc_peak_start_for_sample", "srt_bc_peak_end_for_sample", "srt_bc_peak_width_for_sample",
        "sample_total_reads_in_consensus_peak", "sample_total_fragments_in_consensus_peak", "sample_total_insertions_in_consensus_peak",
        "sample_peak_id"
    }
    sample_peaks_list = []
    can_skip_phase2 = True
    for f in expected_sample_peak_files:
        if not f.exists():
            can_skip_phase2 = False
            break
        df = pl.read_parquet(f)
        if not required_columns.issubset(set(df.columns)):
            can_skip_phase2 = False
            break
        sample_peaks_list.append(df)

    if can_skip_phase2:
        logging.info("Skipping Peak Calling Phase: All sample peaks already called and contain required columns.")

    else:
        # Run per-sample worker for each sample (in parallel, peak calling and SRT refinement)
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = {
                executor.submit(process_pooled_sample_data, task, annotation_df, output_dir): task
                for task in tasks
            }
            for f in tqdm(as_completed(futures), total=len(futures), desc="Processing Samples (fragment peaks -> SRT-BC correction -> SRT-BC peaks)"):
                result_df = f.result()
                if result_df is not None and not result_df.is_empty():
                    # Ensure sample_peak_id is always a string (never a tuple)
                    if "sample_peak_id" in result_df.columns:
                        if isinstance(result_df["sample_peak_id"][0], tuple):
                            result_df = result_df.with_columns(
                                pl.col("sample_peak_id").apply(lambda x: x[0] if isinstance(x, tuple) else x)
                            )
                    sample_peaks_list.append(result_df)

    if not sample_peaks_list:
        logging.error("No valid peaks generated from any sample.")
        sys.exit(1)

    all_sample_peaks_df = pl.concat(sample_peaks_list, rechunk=True)
    logging.info(f"Defined {all_sample_peaks_df.height} sample-specific peaks.")

    # Save all sample peaks for downstream use
    all_sample_peaks_df.write_parquet(sample_peaks_dir / "all_sample_peaks.parquet")




    # ====== 8. Consensus peak generation ======
    logging.info("--- Phase 3: Generating Consensus Peaks and Final Output ---")

    # Merge all sample-specific SRT-BC peaks to define consensus peaks across all samples/groups
    consensus_peaks_to_merge = all_sample_peaks_df.select("chrom", "srt_bc_peak_start_for_sample", "srt_bc_peak_end_for_sample").to_pandas()
    consensus_bed = pybedtools.BedTool.from_dataframe(consensus_peaks_to_merge).sort().merge()
    consensus_peaks_df = pl.from_pandas(
        consensus_bed.to_dataframe(names=["chrom", "consensus_peak_start", "consensus_peak_end"])
    ).with_columns(
        (pl.col("chrom") + pl.lit(":") + pl.col("consensus_peak_start").cast(pl.String) + pl.lit("-") + pl.col("consensus_peak_end").cast(pl.String)).alias("consensus_peak_id")
    )
    logging.info(f"Created {consensus_peaks_df.height} final consensus peaks.")




    # ====== 9. Sample peak-to-consensus peak mapping ======
    # Intersect each sample's SRT-BC peaks with consensus peaks to propagate per-sample stats to consensus peaks
    # Each (sample_peak, consensus_peak) overlap results in a row with all per-sample stats mapped to the consensus peak
    sample_peaks_for_bed = all_sample_peaks_df.select(
        "chrom",
        pl.col("srt_bc_peak_start_for_sample").alias("start"),
        pl.col("srt_bc_peak_end_for_sample").alias("end"),
        *[c for c in all_sample_peaks_df.columns if c not in {"chrom", "srt_bc_peak_start_for_sample", "srt_bc_peak_end_for_sample", "start", "end"}]
    )
    sample_peaks_bedtool = pybedtools.BedTool.from_dataframe(sample_peaks_for_bed.to_pandas())

    consensus_peaks_for_bed = consensus_peaks_df.select(
        "chrom",
        pl.col("consensus_peak_start").alias("start"),
        pl.col("consensus_peak_end").alias("end"),
        "consensus_peak_id"
    )
    consensus_peaks_bedtool = pybedtools.BedTool.from_dataframe(consensus_peaks_for_bed.to_pandas())

    sample_peak_columns = sample_peaks_for_bed.columns
    consensus_peak_columns = [f"{c}_c" for c in consensus_peaks_for_bed.columns]

    # Each row in mapped_peaks_df contains a sample peak, its mapped consensus peak, and all per-sample stats
    mapped_peaks_df = pl.from_pandas(
        sample_peaks_bedtool.intersect(consensus_peaks_bedtool, wa=True, wb=True).to_dataframe(
            names=sample_peak_columns + consensus_peak_columns
        )
    ).select(
        pl.all().exclude(["chrom_c", "start_c", "end_c"]),
        pl.col("consensus_peak_id_c").alias("consensus_peak_id")
    )

    mapped_peaks_df = mapped_peaks_df.rename({
        "start": "srt_bc_peak_start_for_sample",
        "end": "srt_bc_peak_end_for_sample"
    })




    # ====== 10. Output aggregation ======
    # For each consensus peak and group, aggregate supporting sample stats to compute group-level statistics from all samples of the same group.
    stats_per_group_df = (
        mapped_peaks_df
        .group_by(["consensus_peak_id", "group_name"])
        .agg([
            pl.sum("sample_total_reads_in_consensus_peak").alias("group_total_reads_in_consensus_peak"),
            pl.sum("sample_total_fragments_in_consensus_peak").alias("group_total_fragments_in_consensus_peak"),
            pl.sum("sample_total_insertions_in_consensus_peak").alias("group_total_insertions_in_consensus_peak"),
        ])
    )

    # For each consensus peak, count how many unique samples/groups support it
    stats_per_peak_df = (
        mapped_peaks_df
        .group_by("consensus_peak_id")
        .agg([
            pl.n_unique("sample_name").alias("num_samples_in_consensus_peak"),
            pl.n_unique("group_name").alias("num_groups_in_consensus_peak")
        ])
    )

    # Join consensus peak coordinates and widths to mapped peaks
    mapped_peaks_df = mapped_peaks_df.join(consensus_peaks_df, on="consensus_peak_id", how="left")

    # Add widths for consensus and sample peaks
    mapped_peaks_df = mapped_peaks_df.with_columns([
        (pl.col("consensus_peak_end") - pl.col("consensus_peak_start")).alias("consensus_peak_width"),
        (pl.col("srt_bc_peak_end_for_sample") - pl.col("srt_bc_peak_start_for_sample")).alias("srt_bc_peak_width_for_sample"),
        (pl.col("fragment_peak_end_for_sample") - pl.col("fragment_peak_start_for_sample")).alias("fragment_peak_width_for_sample"),
    ])

    # Join per-group and per-peak summary stats to the final mapped table
    final_results_df = (
        mapped_peaks_df
        .join(stats_per_group_df, on=["consensus_peak_id", "group_name"], how="left")
        .join(stats_per_peak_df, on="consensus_peak_id", how="left")
    )

    # Output final TSV (final_results.tsv) with all required columns and sorted by genomic coordinates
    final_output_columns = [
        "consensus_peak_id",
        "chrom",
        "consensus_peak_start",
        "consensus_peak_end",
        "consensus_peak_width",
        "num_samples_in_consensus_peak",
        "num_groups_in_consensus_peak",
        "sample_name",
        "group_name",
        "library_name",
        "fragment_peak_start_for_sample",
        "fragment_peak_end_for_sample",
        "fragment_peak_width_for_sample",
        "srt_bc_peak_start_for_sample",
        "srt_bc_peak_end_for_sample",
        "srt_bc_peak_width_for_sample",
        "sample_total_reads_in_consensus_peak",
        "sample_total_fragments_in_consensus_peak",
        "sample_total_insertions_in_consensus_peak",
        "group_total_reads_in_consensus_peak",
        "group_total_fragments_in_consensus_peak",
        "group_total_insertions_in_consensus_peak",
    ]

    final_output_file_path = output_dir / "final_results.tsv"
    logging.info(f"Writing final output to {final_output_file_path}")

    final_results_df.select(final_output_columns).sort(["chrom", "consensus_peak_start"]).write_csv(final_output_file_path, separator="\t", null_value="NA")

    logging.info("🧬 Pipeline complete! Final results and column definitions are written.")
    pybedtools.helpers.cleanup()

if __name__ == "__main__":
    main()
