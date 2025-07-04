# generate_insertion_maps.py
"""
Post-alignment script to generate files for genome browser/gene track 
visualizations of unique transposon insertions.

Outputs bedGraph and bigWig files of raw and normalized insertion events
per-sample and per-group.

This script uses a count normalization based on the geometric mean of 
unique SRT barcode counts of all groups, calculated directly from the 
raw counts generated in this script.
"""
# Import necessary standard libraries.
import sys  # For system-specific parameters and functions, like exiting the script.
import argparse  # For parsing command-line arguments.
import logging  # For logging events and errors.
import shutil  # For high-level file operations (e.g., checking for executables).
import subprocess  # For running external commands (like bedGraphToBigWig).
from pathlib import Path  # For object-oriented filesystem paths.
from concurrent.futures import ProcessPoolExecutor, as_completed  # For running tasks in parallel.
from typing import Dict, Optional  # For type hinting.

# Import third-party libraries for data manipulation and bioinformatics.
import polars as pl  # A fast DataFrame library for manipulating structured data.
import pybedtools  # A library for working with genomic interval files (like BED).
from tqdm import tqdm  # A library for creating progress bars.
from umi_tools import UMIClusterer  # For clustering Unique Molecular Identifiers (UMIs).
import numpy as np  # For numerical operations, especially for normalization.

# Define the version of the script.
VERSION = "1.0.0" 


def setup_logging(output_dir: Path) -> None:
    """
    Initializes logging to both a file and standard output (stdout).

    This function sets up a global logger that will write informational
    messages and errors to a log file within the specified output directory,
    and also print them to the console.

    Args:
        output_dir: The directory where the log file 'insertion_maps.log' will be created.
    """
    # Define the full path for the log file.
    log_file = output_dir / "insertion_maps.log"
    
    # Configure the basic logging settings.
    logging.basicConfig(
        level=logging.INFO,  # Set the minimum level of messages to log.
        format="%(asctime)s [%(levelname)s] - %(message)s",  # Define the format for log messages.
        handlers=[
            # A handler to write logs to the specified file, overwriting it if it exists ('w').
            logging.FileHandler(log_file, mode="w"),
            # A handler to stream logs to the console.
            logging.StreamHandler(sys.stdout),
        ],
    )
    logging.info(f"Logging initialized. Log file: {log_file}")


# --- Normalization function ---
def calculate_total_count_size_factors(
    group_maps: Dict[str, pl.DataFrame]
) -> Dict[str, float]:
    """
    Calculates size factors for normalization based on the total number of insertions per group.

    This method computes size factors to normalize for differences in sequencing
    depth or library size between experimental groups. It uses the geometric mean
    of total counts as a stable reference against which each group's total count is compared.

    Args:
        group_maps: A dictionary where keys are group names (str) and values are
                    polars DataFrames containing the raw insertion counts for that group.

    Returns:
        A dictionary where keys are group names and values are their calculated size factors (float).
    """
    logging.info("Calculating Total Count size factors from generated raw counts...")
    
    # 1. Calculate the total library size for each group.
    # This is done by summing the 'unique_count' column for each group's DataFrame.
    library_sizes = {
        name: df["unique_count"].sum() for name, df in group_maps.items()
    }

    # Check if all library sizes are zero. If so, normalization is not possible.
    if not any(library_sizes.values()):
        logging.warning("All groups have zero counts. Cannot calculate size factors.")
        # Return a default factor of 1.0 for each group.
        return {name: 1.0 for name in group_maps.keys()}

    # 2. Use the geometric mean of library sizes as a stable reference point.
    # The geometric mean is less sensitive to extreme outliers than the arithmetic mean.
    # A pseudocount of 1 is added to handle cases where a library size is zero.
    log_totals = np.log([count + 1 for count in library_sizes.values()])
    geo_mean_total = np.exp(np.mean(log_totals))

    # If the geometric mean is close to zero, it means most counts are zero.
    if geo_mean_total <= 1:
        logging.warning("Geometric mean of library sizes is near zero. Defaulting factors to 1.0.")
        return {name: 1.0 for name in group_maps.keys()}

    # 3. The size factor for each group is its total count divided by the geometric mean.
    factors = {
        name: total / geo_mean_total for name, total in library_sizes.items()
    }

    # Create and print a DataFrame for easy inspection of the calculated size factors.
    size_factor_df = pl.DataFrame({
        "group": list(factors.keys()),
        "size_factor": list(factors.values())
    })
    logging.info("Calculated Size Factors:")
    print(size_factor_df)
    
    return factors

# --- Validation function ---
def perform_validation_check(
    group_bgs: Dict[str, pl.DataFrame],
    final_results: Path,
):
    """
    Performs a validation check by comparing this script's raw group counts
    against the main pipeline's final summary file.

    This ensures data consistency and integrity between different stages of the analysis.

    Args:
        group_bgs: A dictionary of the raw group counts calculated by this script.
        final_results: The path to the 'final_results.tsv' file from the main pipeline.
    """
    logging.info("=== Validation that Total Unique Insertions Matches Final Output of multiplex_srt_seq_to_tf_binding.py ===")
    
    # Check if the validation file exists.
    if not final_results.is_file():
        logging.error(f"Validation skipped: final_results.tsv not found at {final_results}")
        return
        
    # Read the final results file into a polars DataFrame.
    fr = pl.read_csv(final_results, separator="\t", null_values="NA")
    
    # Iterate through each group calculated in this script.
    for g, df in group_bgs.items():
        # Sum the total unique insertions for the group from this script's data.
        script_sum = int(df["unique_count"].sum()) if not df.is_empty() else 0
        
        # Calculate the sum from the main pipeline's output for the same group.
        # This requires filtering by group, grouping by peak ID to get the unique count per peak,
        # and then summing those counts.
        pipeline_sum = int(
            fr.filter(pl.col("group_name") == g)
            .group_by("consensus_peak_id")
            .agg(pl.first("group_total_insertions_in_consensus_peak"))[
                "group_total_insertions_in_consensus_peak"
            ]
            .sum()
        )
        
        # Compare the two sums and log the result.
        if script_sum == pipeline_sum:
            logging.info(f"✔ Group '{g}': ✅ MATCH! (Raw Count: {script_sum})")
        else:
            logging.warning(f"✖ Group '{g}': MISMATCH! script_raw_count={script_sum} pipeline_raw_count={pipeline_sum}")


# --- Core function to convert each sample's fragment file ---
def process_sample_insertions(task: tuple) -> Optional[pl.DataFrame]:
    """
    Worker function to generate a bedGraph DataFrame of unique insertions for a single sample.

    This function reads fragment data, identifies insertion sites, deduplicates them
    based on SRT barcodes using UMI clustering, and returns a DataFrame of unique
    insertion locations and their counts.

    Args:
        task: A tuple containing (sample_name, parquet_path, fragment_peaks_path, dist_thresh).

    Returns:
        A polars DataFrame with columns ['chrom', 'start', 'end', 'unique_count']
        representing the unique insertions, or None if an error occurs or no insertions are found.
    """
    # Unpack the task tuple for clarity.
    sample_name, parquet_path, fragment_peaks_path, dist_thresh = task
    try:
        # Read the sample's fragment data from a Parquet file.
        frags = pl.read_parquet(parquet_path)
        # Deduplicate fragment records that are identical except for read count and raw barcode.
        dedup_keys = [c for c in frags.columns if c not in ["sample_barcode_raw", "reads"]]
        frags = frags.group_by(dedup_keys).agg(
            pl.sum("reads").alias("reads"),
            pl.col("sample_barcode_raw").first().alias("sample_barcode_raw")
        )

        # If the fragment peaks file doesn't exist or is empty, we can't proceed.
        if not Path(fragment_peaks_path).is_file(): return None
        fp = pl.read_parquet(fragment_peaks_path)
        if fp.is_empty(): return None

        # Convert polars DataFrames to pandas for compatibility with pybedtools.
        df_frag = frags.select(["chrom", "start", "end", "strand", "srt_bc", "reads"]).to_pandas()
        df_peak = fp.select(["chrom", "fragment_peak_start_for_sample", "fragment_peak_end_for_sample", "sample_peak_id"]).to_pandas()

        # Create BedTool objects from the pandas DataFrames.
        bed_f = pybedtools.BedTool.from_dataframe(df_frag)
        bed_p = pybedtools.BedTool.from_dataframe(df_peak)
        
        # Define column names for the upcoming intersection result.
        cols = list(df_frag.columns) + [c + "_peak" for c in df_peak.columns]
        
        # Intersect fragments with their corresponding peaks to associate them.
        # `wa=True` writes the original entry in A. `wb=True` writes the original entry in B.
        inter = pl.from_pandas(bed_f.intersect(bed_p, wa=True, wb=True).to_dataframe(names=cols))

        # If there's no intersection, there are no insertions within peaks.
        if inter.is_empty(): return None

        # A list to hold the final insertion site summaries for this sample.
        sites = []
        # Group fragments by the peak they fall into.
        for _, grp in inter.group_by("sample_peak_id_peak"):
            # Within each peak, count reads per SRT barcode.
            cnts = grp.group_by("srt_bc").agg(pl.sum("reads"))
            # Create a dictionary of {srt_bc: read_count}.
            rd = {r["srt_bc"]: r["reads"] for r in cnts.to_dicts() if r["srt_bc"]}
            if not rd: continue

            # Initialize a UMI clusterer. 'directional' allows for mismatches in barcodes.
            cl = UMIClusterer(cluster_method="directional")
            # Cluster SRT barcodes based on sequence similarity (threshold) to identify unique molecules.
            # This corrects for sequencing errors in the barcodes.
            clusters = cl({k.encode(): v for k, v in rd.items()}, threshold=dist_thresh)
            
            # Create a map from each original barcode to its representative (the first barcode in its cluster).
            rep = {bc.decode(): cluster[0].decode() for cluster in clusters for bc in cluster}

            # Replace original SRT barcodes with their representative barcode.
            # This effectively deduplicates the fragments based on the clustered barcodes.
            deduped_frags_df = grp.with_columns(pl.col("srt_bc").replace(rep).alias("srt_bc_rep")) \
                .filter(pl.col("srt_bc_rep").is_not_null()) \
                .group_by("chrom", "start", "end", "strand", "srt_bc_rep") \
                .agg(pl.sum("reads").alias("reads_in_molecule"))

            # Skip peaks with very few uniquely fragmented molecules as they might be noise.
            if deduped_frags_df.height < 5: continue

            # For each unique molecule (representative barcode), determine the best approximation of the true insertion site.
                # The 'end' value for '+' strand fragments is the true end of the read as it came off the sequencer.
                    # The reads of '+' strand fragments move 5' <- 3' (opposite of convention).
                    # The most 5' (upstream, smallest) value is the the most proximal to the true transposon/genomic DNA junction.

                # The 'start' value for '-' strand coords is the true end of the read as it came off the sequencer.
                    # The reads of '-' strand fragments move 5' -> 3' (opposite of convention).
                    # The most 3' (downstream, largest) value is the the most proximal to the true transposon/genomic DNA junction.
            
            summary = deduped_frags_df.group_by("srt_bc_rep").agg([
                pl.col("start").max().alias("max_start"), # Get the maximum start coordinate.
                pl.col("end").min().alias("min_end"),     # Get the minimum end coordinate.
                pl.first("strand").alias("strand"),       # Get the strand.
                pl.first("chrom").alias("chrom"),         # Get the chromosome.
            ]).with_columns(

                # For the '+' strand, this is the start of the fragment.
                # For the '-' strand, this is the end of the fragment.
                pl.when(pl.col("strand") == "+").then(pl.col("min_end")).otherwise(pl.col("max_start")).alias("site")
            ).select(["chrom", "site"]) # Keep only chromosome and the calculated site.
            sites.append(summary)

        
        # If no sites were identified across all peaks, return None.
        if not sites: return None
        
        # Concatenate the results from all peaks into a single DataFrame.
        final = pl.concat(sites)
        
        # Group by the exact insertion site to count how many unique molecules map to each site.
        return final.group_by(["chrom", "site"]).agg(pl.count().alias("unique_count")) \
            .rename({"site": "start"}) \
            .with_columns((pl.col("start") + 1).alias("end")) \
            .select(["chrom", "start", "end", "unique_count"]).sort(["chrom", "start"])
            
    except Exception as e:
        # Log any errors that occur during processing for a specific sample.
        logging.error(f"Error in sample {sample_name}: {e}", exc_info=True)
        return None


def main():
    """The main entry point for the script."""
    # --- Argument Parsing ---
    # Set up the argument parser with a description.
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("--output_dir_of_multiplex_srt_seq_to_tf_binding", required=True, help="Path to the main pipeline's output directory.")
    p.add_argument("--ins_map_output_dir", required=True, help="Path to the directory where insertion maps will be saved.")
    p.add_argument("--annotation_file", required=True, help="Path to the annotation file mapping samples to groups.")
    p.add_argument("--chrom_sizes", required=True, help="Path to a two-column chromosome sizes file (e.g., 'chr1\t248956422').")
    p.add_argument("--final_results_tsv", default=None, help="Path to final_results.tsv for validation. If not provided, it's inferred from the main output directory.")
    p.add_argument("--srt_bc_dist_threshold", type=int, default=1, help="Hamming distance threshold for clustering SRT barcodes.")
    p.add_argument("--workers", type=int, default=8, help="Number of parallel processes to use.")
    args = p.parse_args()

    # --- Setup ---
    # Create the main output directory.
    out = Path(args.ins_map_output_dir)
    out.mkdir(parents=True, exist_ok=True)
    
    # Initialize logging.
    setup_logging(out)
    logging.info(f"Insertion Mapper v{VERSION}")

    # Check for the required `bedGraphToBigWig` tool in the system's PATH.
    if not shutil.which("bedGraphToBigWig"):
        logging.error("`bedGraphToBigWig` not found in PATH. This tool is required to generate bigWig files.")
        sys.exit(1)
    logging.info("Found `bedGraphToBigWig` executable.")

    # --- Path and Directory Management ---
    # Define paths based on the main pipeline's output structure.
    base = Path(args.output_dir_of_multiplex_srt_seq_to_tf_binding)
    collapsed = base / "collapsed_per_sample"
    fragment_peaks_dir = base / "fragment_peaks"
    sizes = Path(args.chrom_sizes)

    # Infer the path to the final results TSV if not explicitly provided.
    if args.final_results_tsv:
        final_tsv = Path(args.final_results_tsv)
    else:
        final_tsv = base / "final_results.tsv"
        logging.info(f"--final_results_tsv not provided, defaulting to: {final_tsv}")

    # Ensure the necessary input directories from the main pipeline exist.
    if not (collapsed.is_dir() and fragment_peaks_dir.is_dir()):
        logging.error("Missing required 'collapsed_per_sample/' or 'fragment_peaks/' directories.")
        sys.exit(1)
    
    # Define and create all necessary output subdirectories.
    dirs = {
        "sample_bg": out / "raw_unique_insertions_per_sample_bedgraph",
        "sample_bw": out / "raw_unique_insertions_per_sample_bigwig",
        "group_bg":  out / "raw_unique_insertion_count_per_group_bedgraph",
        "group_bw":  out / "raw_unique_insertion_count_per_group_bigwig",
        "norm_bg":   out / "size_normalized_unique_insertions_per_group_bedgraph",
        "norm_bw":   out / "size_normalized_unique_insertions_per_group_bigwig"
    }
    for d in dirs.values():
        d.mkdir(exist_ok=True, parents=True)

    # Set a temporary directory for pybedtools to avoid cluttering the system's /tmp.
    pybedtools_tmp_dir = out / "pybedtools_tmp"
    pybedtools_tmp_dir.mkdir(exist_ok=True, parents=True)
    pybedtools.helpers.set_tempdir(pybedtools_tmp_dir)

    # --- Parallel Processing of Samples ---
    # Create a list of tasks, where each task is a tuple of arguments for the worker function.
    tasks = [(f.stem, str(f), str(fragment_peaks_dir / f"{f.stem}_fragment_peaks.parquet"), args.srt_bc_dist_threshold) for f in collapsed.glob("*.parquet")]
    
    # A dictionary to store the resulting DataFrames from each sample.
    sample_maps: Dict[str, pl.DataFrame] = {}

    # Use a ProcessPoolExecutor to run `process_sample_insertions` in parallel.
    with ProcessPoolExecutor(args.workers) as exe:
        # Submit all tasks to the pool.
        futures = {exe.submit(process_sample_insertions, t): t[0] for t in tasks}
        
        # Process results as they are completed, with a progress bar.
        for fut in tqdm(as_completed(futures), total=len(futures), desc="Mapping Insertions Per Sample"):
            s = futures[fut]  # Get the sample name associated with the future.
            df = fut.result() # Get the result (DataFrame or None).
            
            # If the processing was successful and yielded a non-empty DataFrame...
            if df is not None and not df.is_empty():
                sample_maps[s] = df # Store the result.
                bg = dirs["sample_bg"] / f"{s}.bedgraph"
                bw = dirs["sample_bw"] / f"{s}.bw"
                # Write the raw bedGraph file for the sample.
                df.write_csv(bg, separator="\t", include_header=False)
                
                # Convert the bedGraph to a BigWig file using the external tool.
                try:
                    subprocess.run(["bedGraphToBigWig", str(bg), str(sizes), str(bw)], check=True, capture_output=True, text=True)
                except subprocess.CalledProcessError as e:
                    logging.error(f"Failed to create BigWig for {s}: {e.stderr.strip()}")

    # If no samples were successfully processed, exit the script.
    if not sample_maps:
        logging.error("No sample bedGraphs could be created. Exiting.")
        sys.exit(1)

    # --- Grouping and Normalization ---
    # Read the annotation file to map samples to their respective groups.
    with open(args.annotation_file, 'r') as f:
        # Auto-detect the separator (tab or comma).
        sep_char = '\t' if '\t' in f.readline() else ','
    ann = pl.read_csv(args.annotation_file, separator=sep_char, has_header=True)
    # Clean up column names and sample names.
    ann = ann.rename({col: col.strip() for col in ann.columns}).with_columns(pl.col("sample_name").str.strip_chars())
    # Create a mapping from group names to a list of their associated sample names.
    group_to_samples = ann.group_by("group_name").agg(pl.col("sample_name"))
    
    # A dictionary to store the aggregated DataFrames for each group.
    group_maps: Dict[str, pl.DataFrame] = {}
    
    # Iterate over each group to aggregate sample data.
    for group_name, sample_name_list in group_to_samples.iter_rows():
        # Get a list of samples in this group that were successfully processed.
        samp_list = [s for s in sample_name_list if s in sample_maps]
        if not samp_list: continue # Skip if no samples from this group are available.
        
        # Concatenate DataFrames from all samples in the group.
        concat = pl.concat([sample_maps[s] for s in samp_list])
        
        # Aggregate the counts for identical insertion sites across the group.
        agg = concat.group_by(["chrom", "start", "end"]).agg(pl.sum("unique_count")).sort(["chrom", "start"])
        group_maps[group_name] = agg
        
        # Write the raw bedGraph and BigWig files for the group.
        bg = dirs["group_bg"] / f"{group_name}.bedgraph"
        bw = dirs["group_bw"] / f"{group_name}.bw"
        agg.write_csv(bg, separator="\t", include_header=False)
        try:
            subprocess.run(["bedGraphToBigWig", str(bg), str(sizes), str(bw)], check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"Failed to create BigWig for {group_name}: {e.stderr.strip()}")

    # Perform normalization only if there is more than one group to compare.
    if len(group_maps) > 1:
        # Calculate size factors based on the raw group counts generated in this script.
        factors = calculate_total_count_size_factors(group_maps)
        
        # Apply normalization to each group.
        for g, df in group_maps.items():
            fct = factors.get(g, 1.0) # Get the size factor for the group.
            if fct > 0:
                # Normalize the counts by dividing by the size factor.
                norm_df = df.with_columns((pl.col("unique_count") / fct).alias("normalized_count"))
                norm_output_df = norm_df.select(["chrom", "start", "end", "normalized_count"])
                
                # Write the normalized bedGraph and BigWig files.
                bg = dirs["norm_bg"] / f"{g}.bedgraph"
                bw = dirs["norm_bw"] / f"{g}.bw"
                norm_output_df.write_csv(bg, separator="\t", include_header=False)
                try:
                    subprocess.run(["bedGraphToBigWig", str(bg), str(sizes), str(bw)], check=True, capture_output=True, text=True)
                except subprocess.CalledProcessError as e:
                    logging.error(f"Failed to create normalized BigWig for {g}: {e.stderr.strip()}")
            else:
                logging.warning(f"Cannot normalize {g}, invalid factor {fct}")
    else:
        logging.warning(f"Skipping normalization: only {len(group_maps)} group(s) found.")

    # --- Final Validation ---
    # Perform the validation check at the end of the script.
    perform_validation_check(group_maps, final_tsv)

    logging.info("✅ All insertion maps generated.")


if __name__ == "__main__":
    main()
