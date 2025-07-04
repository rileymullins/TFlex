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
import sys
import argparse
import logging
import shutil
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, Optional
import polars as pl
import pybedtools
from tqdm import tqdm
from umi_tools import UMIClusterer
import numpy as np

VERSION = "1.0.0" 


def setup_logging(output_dir: Path) -> None:
    """Initializes logging to both file and stdout."""
    log_file = output_dir / "insertion_maps.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] - %(message)s",
        handlers=[
            logging.FileHandler(log_file, mode="w"),
            logging.StreamHandler(sys.stdout),
        ],
    )
    logging.info(f"Logging initialized. Log file: {log_file}")


# --- START: UPDATED NORMALIZATION FUNCTION ---
def calculate_total_count_size_factors(
    group_maps: Dict[str, pl.DataFrame]
) -> Dict[str, float]:
    """
    Calculates size factors based on the total number of insertions per group
    using the provided group_maps dictionary.
    """
    logging.info("Calculating Total Count size factors from generated raw counts...")
    
    # 1. Calculate the total library size for each group from the input dictionary.
    library_sizes = {
        name: df["unique_count"].sum() for name, df in group_maps.items()
    }

    if not any(library_sizes.values()):
        logging.warning("All groups have zero counts. Cannot calculate size factors.")
        return {name: 1.0 for name in group_maps.keys()}

    # 2. Use the geometric mean of library sizes as a stable reference point.
    log_totals = np.log([count + 1 for count in library_sizes.values()])
    geo_mean_total = np.exp(np.mean(log_totals))

    if geo_mean_total <= 1: # geo_mean_total can be 1 if all counts are 0
        logging.warning("Geometric mean of library sizes is near zero. Defaulting factors to 1.0.")
        return {name: 1.0 for name in group_maps.keys()}

    # 3. The size factor is each library's total size divided by the reference.
    factors = {
        name: total / geo_mean_total for name, total in library_sizes.items()
    }

    size_factor_df = pl.DataFrame({
        "group": list(factors.keys()),
        "size_factor": list(factors.values())
    })
    logging.info("Calculated Size Factors:")
    print(size_factor_df)
    
    return factors
# --- END: UPDATED NORMALIZATION FUNCTION ---


def perform_validation_check(
    group_bgs: Dict[str, pl.DataFrame],
    final_results: Path,
):
    """
    Performs a validation check by comparing the script's raw group counts
    against the main pipeline's final_results.tsv.
    """
    logging.info("=== Validation that Total Unique Insertions Matches Final Output of multiplex_srt_seq_to_tf_binding.py ===")
    if not final_results.is_file():
        logging.error(f"Validation skipped: final_results.tsv not found at {final_results}")
        return
        
    fr = pl.read_csv(final_results, separator="\t", null_values="NA")
    
    for g, df in group_bgs.items():
        script_sum = int(df["unique_count"].sum()) if not df.is_empty() else 0
        
        # Calculate the sum from the pipeline output for the same group
        pipeline_sum = int(
            fr.filter(pl.col("group_name") == g)
            .group_by("consensus_peak_id")
            .agg(pl.first("group_total_insertions_in_consensus_peak"))[
                "group_total_insertions_in_consensus_peak"
            ]
            .sum()
        )
        
        if script_sum == pipeline_sum:
            logging.info(f"✔ Group '{g}': MATCH! (Raw Count: {script_sum})")
        else:
            logging.warning(f"✖ Group '{g}': MISMATCH! script_raw_count={script_sum} pipeline_raw_count={pipeline_sum}")


def process_sample_insertions(task: tuple) -> Optional[pl.DataFrame]:
    """Worker: generate per-sample bedGraph DataFrame of unique insertions."""
    sample_name, parquet_path, fragment_peaks_path, dist_thresh = task
    try:
        frags = pl.read_parquet(parquet_path)
        dedup_keys = [c for c in frags.columns if c not in ["sample_barcode_raw", "reads"]]
        frags = frags.group_by(dedup_keys).agg(
            pl.sum("reads").alias("reads"),
            pl.col("sample_barcode_raw").first().alias("sample_barcode_raw")
        )

        if not Path(fragment_peaks_path).is_file(): return None
        fp = pl.read_parquet(fragment_peaks_path)
        if fp.is_empty(): return None

        df_frag = frags.select(["chrom", "start", "end", "strand", "srt_bc", "reads"]).to_pandas()
        df_peak = fp.select(["chrom", "fragment_peak_start_for_sample", "fragment_peak_end_for_sample", "sample_peak_id"]).to_pandas()

        bed_f = pybedtools.BedTool.from_dataframe(df_frag)
        bed_p = pybedtools.BedTool.from_dataframe(df_peak)
        cols = list(df_frag.columns) + [c + "_peak" for c in df_peak.columns]
        inter = pl.from_pandas(bed_f.intersect(bed_p, wa=True, wb=True).to_dataframe(names=cols))

        if inter.is_empty(): return None

        sites = []
        for _, grp in inter.group_by("sample_peak_id_peak"):
            cnts = grp.group_by("srt_bc").agg(pl.sum("reads"))
            rd = {r["srt_bc"]: r["reads"] for r in cnts.to_dicts() if r["srt_bc"]}
            if not rd: continue

            cl = UMIClusterer(cluster_method="directional")
            clusters = cl({k.encode(): v for k, v in rd.items()}, threshold=dist_thresh)
            rep = {bc.decode(): cluster[0].decode() for cluster in clusters for bc in cluster}

            deduped_frags_df = grp.with_columns(pl.col("srt_bc").replace(rep).alias("srt_bc_rep")) \
                .filter(pl.col("srt_bc_rep").is_not_null()) \
                .group_by("chrom", "start", "end", "strand", "srt_bc_rep") \
                .agg(pl.sum("reads").alias("reads_in_molecule"))

            if deduped_frags_df.height < 5: continue

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
        final = pl.concat(sites)
        return final.group_by(["chrom", "site"]).agg(pl.count().alias("unique_count")) \
            .rename({"site": "start"}).with_columns((pl.col("start") + 1).alias("end")) \
            .select(["chrom", "start", "end", "unique_count"]).sort(["chrom", "start"])
    except Exception as e:
        logging.error(f"Error in sample {sample_name}: {e}", exc_info=True)
        return None


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("--output_dir_of_multiplex_srt_seq_to_tf_binding", required=True)
    p.add_argument("--ins_map_output_dir", required=True)
    p.add_argument("--annotation_file", required=True)
    p.add_argument("--chrom_sizes", required=True, help="Path to two-column chrom.sizes (chr<tab>length).")
    p.add_argument("--final_results_tsv", default=None, help="Path to final_results.tsv for validation. If not provided, it's inferred from the main output directory.")
    p.add_argument("--srt_bc_dist_threshold", type=int, default=1)
    p.add_argument("--workers", type=int, default=8)
    args = p.parse_args()

    out = Path(args.ins_map_output_dir)
    out.mkdir(parents=True, exist_ok=True)
    setup_logging(out)
    logging.info(f"Insertion Mapper v{VERSION}")

    if not shutil.which("bedGraphToBigWig"):
        logging.error("`bedGraphToBigWig` not found in PATH.")
        sys.exit(1)
    logging.info("Found `bedGraphToBigWig` executable.")

    base = Path(args.output_dir_of_multiplex_srt_seq_to_tf_binding)
    collapsed = base / "collapsed_per_sample"
    fragment_peaks_dir = base / "fragment_peaks"
    sizes = Path(args.chrom_sizes)

    if args.final_results_tsv:
        final_tsv = Path(args.final_results_tsv)
    else:
        final_tsv = base / "final_results.tsv"
        logging.info(f"--final_results_tsv not provided, defaulting to: {final_tsv}")

    if not (collapsed.is_dir() and fragment_peaks_dir.is_dir()):
        logging.error("Missing required 'collapsed_per_sample/' or 'fragment_peaks/' directories.")
        sys.exit(1)
    
    dirs = { "sample_bg": out / "raw_unique_insertions_per_sample_bedgraph", "sample_bw": out / "raw_unique_insertions_per_sample_bigwig", "group_bg":  out / "raw_unique_insertion_count_per_group_bedgraph", "group_bw":  out / "raw_unique_insertion_count_per_group_bigwig", "norm_bg":   out / "size_normalized_unique_insertions_per_group_bedgraph", "norm_bw":   out / "size_normalized_unique_insertions_per_group_bigwig" }
    for d in dirs.values():
        d.mkdir(exist_ok=True, parents=True)

    pybedtools_tmp_dir = out / "pybedtools_tmp"
    pybedtools_tmp_dir.mkdir(exist_ok=True, parents=True)
    pybedtools.helpers.set_tempdir(pybedtools_tmp_dir)

    tasks = [(f.stem, str(f), str(fragment_peaks_dir / f"{f.stem}_fragment_peaks.parquet"), args.srt_bc_dist_threshold) for f in collapsed.glob("*.parquet")]
    sample_maps: Dict[str, pl.DataFrame] = {}

    with ProcessPoolExecutor(args.workers) as exe:
        futures = {exe.submit(process_sample_insertions, t): t[0] for t in tasks}
        for fut in tqdm(as_completed(futures), total=len(futures), desc="Mapping Insertions Per Sample"):
            s = futures[fut]
            df = fut.result()
            if df is not None and not df.is_empty():
                sample_maps[s] = df
                bg = dirs["sample_bg"] / f"{s}.bedgraph"
                bw = dirs["sample_bw"] / f"{s}.bw"
                df.write_csv(bg, separator="\t", include_header=False)
                try: subprocess.run(["bedGraphToBigWig", str(bg), str(sizes), str(bw)], check=True, capture_output=True, text=True)
                except subprocess.CalledProcessError as e: logging.error(f"Failed to create BigWig for {s}: {e.stderr.strip()}")

    if not sample_maps:
        logging.error("No sample bedGraphs could be created. Exiting.")
        sys.exit(1)

    with open(args.annotation_file, 'r') as f:
        sep_char = '\t' if '\t' in f.readline() else ','
    ann = pl.read_csv(args.annotation_file, separator=sep_char, has_header=True)
    ann = ann.rename({col: col.strip() for col in ann.columns}).with_columns(pl.col("sample_name").str.strip_chars())
    group_to_samples = ann.group_by("group_name").agg(pl.col("sample_name"))
    
    group_maps: Dict[str, pl.DataFrame] = {}
    for group_name, sample_name_list in group_to_samples.iter_rows():
        samp_list = [s for s in sample_name_list if s in sample_maps]
        if not samp_list: continue
        concat = pl.concat([sample_maps[s] for s in samp_list])
        agg = concat.group_by(["chrom", "start", "end"]).agg(pl.sum("unique_count")).sort(["chrom", "start"])
        group_maps[group_name] = agg
        bg = dirs["group_bg"] / f"{group_name}.bedgraph"
        bw = dirs["group_bw"] / f"{group_name}.bw"
        agg.write_csv(bg, separator="\t", include_header=False)
        try: subprocess.run(["bedGraphToBigWig", str(bg), str(sizes), str(bw)], check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e: logging.error(f"Failed to create BigWig for {group_name}: {e.stderr.strip()}")

    if len(group_maps) > 1:
        # --- Call the normalization function using the script's own generated counts ---
        factors = calculate_total_count_size_factors(group_maps)
        
        for g, df in group_maps.items():
            fct = factors.get(g, 1.0)
            if fct > 0:
                norm_df = df.with_columns((pl.col("unique_count") / fct).alias("normalized_count"))
                norm_output_df = norm_df.select(["chrom", "start", "end", "normalized_count"])
                bg = dirs["norm_bg"] / f"{g}.bedgraph"
                bw = dirs["norm_bw"] / f"{g}.bw"
                norm_output_df.write_csv(bg, separator="\t", include_header=False)
                try: subprocess.run(["bedGraphToBigWig", str(bg), str(sizes), str(bw)], check=True, capture_output=True, text=True)
                except subprocess.CalledProcessError as e: logging.error(f"Failed to create normalized BigWig for {g}: {e.stderr.strip()}")
            else:
                logging.warning(f"Cannot normalize {g}, invalid factor {fct}")
    else:
        logging.warning(f"Skipping normalization: only {len(group_maps)} group(s) found.")

    # Validation is performed at the end using the specified final_tsv path
    perform_validation_check(group_maps, final_tsv)

    logging.info("✅ All insertion maps generated.")


if __name__ == "__main__":
    main()
