# generate_insertion_maps.py
"""
Post-processing script for a calling card pipeline.
Generates per-sample and per-group insertion maps (bedGraph/BigWig).
It also generates binned and summed signal tracks from normalized data within a specified window size.

This script uses a total count normalization method based on the geometric mean of total unique insertions across samples.

Because the exact insertion site is approximate, it helps to bin for visualization. 
This is helpful because we do not have the exact position of where the representative read of each
SRT barcode would hit the junction. Those that do hit the junction tend to "stack up" already.


"""
# --- Standard Library Imports ---
import sys
import argparse
import logging
import shutil
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from typing import Dict, Optional, List
import math
import polars as pl
import pybedtools
from tqdm import tqdm
from umi_tools import UMIClusterer
import numpy as np
import pyBigWig

# --- Script Version ---
VERSION = "1.0.0"


# --- Function Definitions ---

def setup_logging(output_dir: Path) -> None:
    """
    Initializes logging to write messages to both a file and the console (stdout).
    
    Args:
        output_dir: The directory where the log file 'insertion_maps.log' will be saved.
    """
    log_file = output_dir / "insertion_maps.log"
    logging.basicConfig(
        level=logging.INFO,  # Set the minimum level of messages to log.
        format="%(asctime)s [%(levelname)s] - %(message)s",  # Define the log message format.
        handlers=[
            logging.FileHandler(log_file, mode="w"),  # Handler to write to a file, overwriting old logs.
            logging.StreamHandler(sys.stdout),      # Handler to print logs to the console.
        ],
    )
    logging.info(f"Logging initialized. Log file: {log_file}")


def calculate_total_count_size_factors(
    group_maps: Dict[str, pl.DataFrame]
) -> Dict[str, float]:
    """
    Calculates size factors to normalize for library depth between experimental groups.

    This method finds a normalization factor by calculating the geometric mean of the total unique
    insertion counts across all groups. 
    Each group's total unique count is then divided by this normalization factor.

    Args:
        group_maps: A dictionary where keys are group names and values are DataFrames
                    containing the raw, aggregated insertion counts for that group.

    Returns:
        A dictionary mapping each group name to its calculated size factor.
    """
    logging.info("Calculating Total Count size factors from generated raw counts...")
    
    # 1. Calculate the total number of unique insertions for each group.
    library_sizes = {
        name: df["unique_count"].sum() for name, df in group_maps.items()
    }


    # 2. Use the geometric mean of library sizes as a stable reference.
    # The geometric mean is less sensitive to extreme outliers than the arithmetic mean,
    # providing a more robust reference for normalization.
    # A pseudocount of +1 is added to handle groups with zero counts.
    log_totals = np.log([count + 1 for count in library_sizes.values()])
    geo_mean_total = np.exp(np.mean(log_totals))

    # Handle case where the geometric mean is near zero.
    if geo_mean_total <= 1:
        logging.warning("Geometric mean of library sizes is near zero. Defaulting factors to 1.0.")
        return {name: 1.0 for name in group_maps.keys()}

    # 3. The size factor is each library's total size divided by the reference.
    factors = {
        name: total / geo_mean_total for name, total in library_sizes.items()
    }

    # Log the results.
    size_factor_df = pl.DataFrame({
        "group": list(factors.keys()),
        "size_factor": list(factors.values())
    })
    logging.info("Calculated Size Factors:")
    print(size_factor_df)
    
    return factors


def perform_validation_check(
    group_bgs: Dict[str, pl.DataFrame],
    final_results: Path,
):
    """
    Performs an integrity check by comparing the script's raw group counts
    against the main pipeline's final summary file ('final_results.tsv').
    
    They should be equal because the processing to count unique insertions is identical.
    This script adds the extra feature of "pushing" all fragments of the same SRT barcode to the single 
    position that is the closest to the transpsoson, which may be at the transpsoson.

    Args:
        group_bgs: A dictionary of the raw group counts calculated by this script.
        final_results: Path to the 'final_results.tsv' file from the main pipeline.
    """
    logging.info("=== Validation Check ===")
    if not final_results.is_file():
        logging.error(f"Validation skipped: final_results.tsv not found at {final_results}")
        return
        
    fr = pl.read_csv(final_results, separator="\t", null_values="NA")
    
    # For each group, compare the sum calculated here with the sum from the main pipeline.
    for g, df in group_bgs.items():
        script_sum = int(df["unique_count"].sum()) if not df.is_empty() else 0
        
        # Calculate the sum from the pipeline output for the same group.
        # This ensures that the total number of insertions matches between pipeline stages.
        pipeline_sum = int(
            fr.filter(pl.col("group_name") == g)
            .group_by("consensus_peak_id")
            .agg(pl.first("group_total_insertions_in_consensus_peak"))[
                "group_total_insertions_in_consensus_peak"
            ]
            .sum()
        )
        
        if script_sum == pipeline_sum:
            logging.info(f"✔ Group '{g}': Total counts MATCH! (Raw Count: {script_sum})")
        else:
            logging.warning(f"✖ Group '{g}': Total counts MISMATCH! script_raw_count={script_sum} pipeline_raw_count={pipeline_sum}")


def process_sample_insertions(task: tuple) -> Optional[pl.DataFrame]:
    """
    Worker function that processes a single sample's fragment data to identify
    unique transposon insertion sites. This is the main computational task that
    is run in parallel.
    This is very similar to the function of the other script to count unique insertions
    per sample in its sample peak.

    Args:
        task: A tuple containing (sample_name, parquet_path, fragment_peaks_path, dist_thresh).

    Returns:
        A Polars DataFrame with columns ['chrom', 'start', 'end', 'unique_count']
        representing unique insertions, or None if an error occurs.
    """
    sample_name, parquet_path, fragment_peaks_path, dist_thresh = task
    try:
        # Read and pre-process fragment data for the sample.
        frags = pl.read_parquet(parquet_path)
        dedup_keys = [c for c in frags.columns if c not in ["sample_barcode_raw", "reads"]]
        frags = frags.group_by(dedup_keys).agg(
            pl.sum("reads").alias("reads"),
            pl.col("sample_barcode_raw").first().alias("sample_barcode_raw")
        )

        # Skip if the associated fragment peak file doesn't exist or is empty.
        if not Path(fragment_peaks_path).is_file(): return None
        fp = pl.read_parquet(fragment_peaks_path)
        if fp.is_empty(): return None

        # Convert to pandas for pybedtools compatibility.
        df_frag = frags.select(["chrom", "start", "end", "strand", "srt_bc", "reads"]).to_pandas()
        df_peak = fp.select(["chrom", "fragment_peak_start_for_sample", "fragment_peak_end_for_sample", "sample_peak_id"]).to_pandas()

        # Create BedTool objects and find fragments that fall within defined peak regions.
        bed_f = pybedtools.BedTool.from_dataframe(df_frag)
        bed_p = pybedtools.BedTool.from_dataframe(df_peak)
        cols = list(df_frag.columns) + [c + "_peak" for c in df_peak.columns]
        inter = pl.from_pandas(bed_f.intersect(bed_p, wa=True, wb=True).to_dataframe(names=cols))

        if inter.is_empty(): return None

        sites = []
        # Process fragments on a peak-by-peak basis.
        for _, grp in inter.group_by("sample_peak_id_peak"):
            # Count reads associated with each SRT barcode within this peak.
            cnts = grp.group_by("srt_bc").agg(pl.sum("reads"))
            rd = {r["srt_bc"]: r["reads"] for r in cnts.to_dicts() if r["srt_bc"]}
            if not rd: continue

            # Cluster SRT barcodes to correct for sequencing errors. This groups similar
            # barcodes together, treating them as originating from the same molecule.
            cl = UMIClusterer(cluster_method="directional")
            clusters = cl({k.encode(): v for k, v in rd.items()}, threshold=dist_thresh)
            rep = {bc.decode(): cluster[0].decode() for cluster in clusters for bc in cluster}

            # Replace original SRT barcodes with their representative from the cluster.
            deduped_frags_df = grp.with_columns(pl.col("srt_bc").replace(rep).alias("srt_bc_rep")) \
                .filter(pl.col("srt_bc_rep").is_not_null()) \
                .group_by("chrom", "start", "end", "strand", "srt_bc_rep") \
                .agg(pl.sum("reads").alias("reads_in_molecule"))

            # Skip peaks with very few unique molecules to reduce noise.
            if deduped_frags_df.height < 5: continue

            # For each unique molecule (srt_bc_rep), determine its single insertion site that is most proximal to the transpsoson.
            # For '+' strand fragments, that is the most upstream position. Since the true end of the read for '+' strand reads is
            # the 'end' coordinate, we take the smallest (most upstream) end coordiante of that SRT barcode as the best approximation.
            
            # We do the same for the '-' strand, just with the opposit logic (max start value)
            
            summary = deduped_frags_df.group_by("srt_bc_rep").agg([
                pl.col("start").max().alias("max_start"),
                pl.col("end").min().alias("min_end"),
                pl.first("strand").alias("strand"),
                pl.first("chrom").alias("chrom"),
            ]).with_columns(
                pl.when(pl.col("strand") == "+").then(pl.col("min_end")).otherwise(pl.col("max_start")).alias("site")
            ).select(["chrom", "site"])
            sites.append(summary)

        if not sites: return None
        
        # Aggregate all identified sites across all peaks for the sample.
        final = pl.concat(sites)
        
        # Count how many unique barcodes mapped to each exact coordinate.
        # This is the basis of generating the insertion maps for genome browser views.
        
        return final.group_by(["chrom", "site"]).agg(pl.count().alias("unique_count")) \
            .rename({"site": "start"}).with_columns((pl.col("start") + 1).alias("end")) \
            .select(["chrom", "start", "end", "unique_count"]).sort(["chrom", "start"])
            
    except Exception as e:
        logging.error(f"Error in sample {sample_name}: {e}", exc_info=True)
        return None
        
        
def generate_binned_sum_track(
    group_name: str,
    normalized_bw_path: Path,
    window_size: int,
    output_dir_bw: Path
) -> None:
    """
    Generates a BigWig track of summed signals in fixed windows from an input BigWig file.
    Uses true half-open intervals [start, end) so that adjacent bins tile the genome
    without overlap or gaps, and render correctly in genome browsers.

    Args:
        group_name: Name of the experimental group.
        normalized_bw_path: Path to the input normalized BigWig file.
        window_size: Size of each bin in base pairs.
        output_dir_bw: Directory for the output binned BigWig.
    """
    logging.info(f"--- Generating binned sum track for group '{group_name}' (window: {window_size} bp) ---")
    binned_bw_path = output_dir_bw / f"{group_name}.sum_{window_size}bp.bw"

    try:
        bw_in = pyBigWig.open(str(normalized_bw_path))
        if not bw_in:
            logging.error(f"Could not open input BigWig: {normalized_bw_path}")
            return

        bw_out = pyBigWig.open(str(binned_bw_path), "w")
        if not bw_out:
            logging.error(f"Could not open output BigWig for writing: {binned_bw_path}")
            bw_in.close()
            return

        # Transfer header (chrom sizes)
        bw_out.addHeader(list(bw_in.chroms().items()))

        for chrom, chrom_len in bw_in.chroms().items():
            if chrom_len == 0:
                continue

            # Load per-base signal, replace NaN with 0
            arr = np.nan_to_num(bw_in.values(chrom, 0, chrom_len, numpy=True))
            num_bins = math.ceil(chrom_len / window_size)

            all_starts, all_ends, all_values = [], [], []
            for i in range(num_bins):
                start = i * window_size
                # Use half-open intervals [start, end)
                end = min(start + window_size, chrom_len)

                # Sum exactly window_size bases (or fewer for last bin)
                bin_sum = arr[start:end].sum()
                if bin_sum > 0:
                    all_starts.append(start)
                    all_ends.append(end)
                    all_values.append(float(bin_sum))

            if all_starts:
                bw_out.addEntries(
                    [chrom] * len(all_starts),
                    all_starts,
                    ends=all_ends,
                    values=all_values
                )

        bw_in.close()
        bw_out.close()
        logging.info(f"✅ Successfully created half-open binned BigWig: {binned_bw_path}")
    except Exception as e:
        logging.error(f"Binned sum generation for '{group_name}' failed: {e}", exc_info=True)



def main():
    """Main function to orchestrate the entire workflow."""
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("--output_dir_of_multiplex_srt_seq_to_tf_binding", required=True, help="Path to the main pipeline's output directory.")
    p.add_argument("--ins_map_output_dir", required=True, help="Directory to save the output insertion maps.")
    p.add_argument("--annotation_file", required=True, help="Annotation file mapping samples to groups (CSV or TSV).")
    p.add_argument("--chrom_sizes", required=True, help="Path to a two-column chromosome sizes file (e.g., 'chr1\\t248956422').")
    p.add_argument("--final_results_tsv", default=None, help="Path to final_results.tsv for validation. Inferred if not provided.")
    p.add_argument("--srt_bc_dist_threshold", type=int, default=1, help="Hamming distance for clustering SRT barcodes. This must be the same as used for multiplex_srt_seq_to_tf_binding.py.")
    p.add_argument("--workers", type=int, default=8, help="Number of parallel processes to use.")
    p.add_argument("--sum_window_size", type=int, default=50, help="If set, generates additional tracks with counts summed in N bp windows.")
    args = p.parse_args()

    # --- Initial Setup ---
    out = Path(args.ins_map_output_dir)
    out.mkdir(parents=True, exist_ok=True)
    setup_logging(out)
    logging.info(f"Insertion Mapper v{VERSION}")

    # Check for required external command-line tools.
    if not shutil.which("bedGraphToBigWig"):
        logging.error("`bedGraphToBigWig` not found in PATH. This tool is required.")
        sys.exit(1)
    logging.info("Found `bedGraphToBigWig` executable.")

    base = Path(args.output_dir_of_multiplex_srt_seq_to_tf_binding)
    sizes = Path(args.chrom_sizes)

    # Define all output directories.
    dirs = {
        "sample_bg": out / "raw_unique_insertions_per_sample_bedgraph",
        "sample_bw": out / "raw_unique_insertions_per_sample_bigwig",
        "group_bg": out / "raw_unique_insertion_count_per_group_bedgraph",
        "group_bw": out / "raw_unique_insertion_count_per_group_bigwig",
        "norm_bg": out / "size_normalized_unique_insertions_per_group_bedgraph",
        "norm_bw": out / "size_normalized_unique_insertions_per_group_bigwig",
        "binned_sum_bw": out / f"binned_sum_bigwig_{args.sum_window_size}bp" if args.sum_window_size else None,
    }
    for d in dirs.values():
        if d: d.mkdir(exist_ok=True, parents=True)

    # Load annotation file to map samples to groups.
    with open(args.annotation_file, 'r') as f:
        sep_char = '\t' if '\t' in f.readline() else ','
    ann = pl.read_csv(args.annotation_file, separator=sep_char, has_header=True)
    ann = ann.rename({col: col.strip() for col in ann.columns}).with_columns(pl.col("sample_name").str.strip_chars())
    group_to_samples: List[tuple[str, List[str]]] = ann.group_by("group_name").agg(pl.col("sample_name")).rows()

    # --- Logic to skip directly to the binning bigwig function if the other files exist ---
    # Check if the necessary input files already exist. If so, we can skip the entire main processing pipeline to save time.
    skip_main_processing = False
    if args.sum_window_size:
        norm_files_exist = []
        for group_name, _ in group_to_samples:
            norm_bw_path = dirs["norm_bw"] / f"{group_name}.bw"
            norm_files_exist.append(norm_bw_path.is_file())
        
        if all(norm_files_exist):
            logging.info("All required normalized BigWig files found. Skipping main processing and proceeding to binned sum generation.")
            skip_main_processing = True

    # --- Main Processing Workflow (per-sample -> per-group -> normalized) ---
    if not skip_main_processing:
        # Define where to find the per-sample parquet files and their corresponding peak files.
        collapsed = base / "collapsed_per_sample"
        fragment_peaks_dir = base / "fragment_peaks"
        
        # Decide which final_results.tsv to use for validation of total counts - either user-supplied or default.
        if args.final_results_tsv:
            final_tsv = Path(args.final_results_tsv)
        else:
            final_tsv = base / "final_results.tsv"
            logging.info(f"--final_results_tsv not provided, defaulting to: {final_tsv}")
            
        # Ensure the two key input directories from multiplex_srt_seq_to_tf_binding.py actually exist before proceeding.
        if not (collapsed.is_dir() and fragment_peaks_dir.is_dir()):
            logging.error("Missing required 'collapsed_per_sample/' or 'fragment_peaks/' directories.")
            sys.exit(1)
        
        # Set a temporary directory for pybedtools to avoid cluttering system /tmp.
        pybedtools_tmp_dir = out / "pybedtools_tmp"
        pybedtools_tmp_dir.mkdir(exist_ok=True, parents=True)
        pybedtools.helpers.set_tempdir(pybedtools_tmp_dir)

        # 1. Process individual SAMPLES in parallel.
        #   N workers are made, each reading one parquet -> intersecting peaks -> deduplicating UMIs -> writing a bedGraph + BigWig.
        logging.info(f"Starting per-sample insertion mapping with {args.workers} workers...")
        tasks = [
            (
                f.stem,  # sample_name
                str(f),  # path to collapsed_per_sample/{sample}.parquet
                str(fragment_peaks_dir / f"{f.stem}_fragment_peaks.parquet"),
                args.srt_bc_dist_threshold
            )
            for f in collapsed.glob("*.parquet")
        ]
        sample_maps: Dict[str, pl.DataFrame] = {}

        with ProcessPoolExecutor(args.workers) as exe:
            # Submit each sample’s processing job; map Future -> sample name for progress tracking.

            futures = {exe.submit(process_sample_insertions, t): t[0] for t in tasks}
            for fut in tqdm(as_completed(futures), total=len(futures), desc="Mapping Insertions Per Sample"):
                s = futures[fut]
                df = fut.result()
                
                # Only keep samples that have insertion sites.
                if df is not None and not df.is_empty():
                    sample_maps[s] = df
                    
                    # Write out the raw-unique-insertions bedGraph, then convert to BigWig.
                    bg = dirs["sample_bg"] / f"{s}.bedgraph"
                    bw = dirs["sample_bw"] / f"{s}.bw"
                    df.write_csv(bg, separator="\t", include_header=False)
                    try:
                        subprocess.run(["bedGraphToBigWig", str(bg), str(sizes), str(bw)], check=True, capture_output=True, text=True)
                    except subprocess.CalledProcessError as e:
                        logging.error(f"Failed to create BigWig for sample {s}: {e.stderr.strip()}")

        # If no sample produced data, stop.
        if not sample_maps:
            logging.error("No sample bedGraphs could be created. Exiting.")
            sys.exit(1)

        # 2. Aggregate sample maps into GROUP maps.
        #    For each defined group, concatenate its samples’ per-site counts and sum them.
        logging.info("Aggregating samples into groups...")
        group_maps: Dict[str, pl.DataFrame] = {}
        for group_name, sample_name_list in group_to_samples:
        
            # Pick only samples that are actually in sample_name_list
            samp_list = [s for s in sample_name_list if s in sample_maps]
            if not samp_list: continue
            
            # Concatenate all the sample DataFrames of the group, then sum unique_counts per (chrom,start,end).
            concat = pl.concat([sample_maps[s] for s in samp_list])
            agg = concat.group_by(["chrom", "start", "end"]).agg(pl.sum("unique_count")).sort(["chrom", "start"])
            group_maps[group_name] = agg
            
            # Write out raw group bedGraph + BigWig files for the group
            bg = dirs["group_bg"] / f"{group_name}.bedgraph"
            bw = dirs["group_bw"] / f"{group_name}.bw"
            agg.write_csv(bg, separator="\t", include_header=False)
            try:
                subprocess.run(["bedGraphToBigWig", str(bg), str(sizes), str(bw)], check=True, capture_output=True, text=True)
            except subprocess.CalledProcessError as e:
                logging.error(f"Failed to create BigWig for group {group_name}: {e.stderr.strip()}")

        # 3. NORMALIZE group maps.
        #    If there’s more than one group, compute size factors and divide each group’s counts.
        if len(group_maps) > 1:
            logging.info("Normalizing group-level insertion maps...")
            factors = calculate_total_count_size_factors(group_maps)
            
            for g, df in group_maps.items():
                fct = factors.get(g, 1.0)
                if fct > 0:
                    norm_df = df.with_columns((pl.col("unique_count") / fct).alias("normalized_count"))
                    norm_output_df = norm_df.select(["chrom", "start", "end", "normalized_count"])
                    
                    # Write normalized bedGraph + BigWig.
                    bg = dirs["norm_bg"] / f"{g}.bedgraph"
                    bw = dirs["norm_bw"] / f"{g}.bw"
                    norm_output_df.write_csv(bg, separator="\t", include_header=False)
                    try:
                        subprocess.run(["bedGraphToBigWig", str(bg), str(sizes), str(bw)], check=True, capture_output=True, text=True)
                    except subprocess.CalledProcessError as e:
                        logging.error(f"Failed to create normalized BigWig for {g}: {e.stderr.strip()}")
                else:
                    logging.warning(f"Cannot normalize group {g}, invalid factor {fct}")
        else:
            logging.warning(f"Skipping normalization: only {len(group_maps)} group found.")
        
        # 4. Perform final VALIDATION check.
        #    Compare the raw total unique insertion counts in this script agains the final_results.tsv from multiplex_srt_seq_to_tf_binding.py
        perform_validation_check(group_maps, final_tsv)

    # --- Binned Sum Track Generation ---
    # This step runs if --sum_window_size is provided.
    if args.sum_window_size and args.sum_window_size > 0:
    
        logging.info(f"--- Starting Parallel Binned Sum Track Generation with threads (Workers: {args.workers}) ---")

        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = {}
            for group_name, _ in group_to_samples:
                bw_path = dirs["norm_bw"] / f"{group_name}.bw"
                if not bw_path.is_file():
                    logging.warning(f"Skipping binned sum for '{group_name}': normalized BigWig not found.")
                    continue

                # Submit each group to be processed by the generate_binned_sum_track function in a separate thread.
                fut = executor.submit(
                    generate_binned_sum_track,
                    group_name,
                    bw_path,
                    args.sum_window_size,
                    dirs["binned_sum_bw"],
                )
                futures[fut] = group_name

            # Wait for all submitted tasks to complete and check for any exceptions.
            for fut in as_completed(futures):
                g = futures[fut]
                try:
                    fut.result()
                except Exception as e:
                    logging.error(f"Binned sum for '{g}' failed: {e}")

    logging.info("✅ All insertion maps generated.")


if __name__ == "__main__":
    main()
