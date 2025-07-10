# step_2_fragments_to_insertion_maps.py

"""
Post-processing script for a calling card pipeline.
Generates per-sample and per-group insertion maps (bedGraph/BigWig).
As a final step, it converts the per-group bedGraphs into expanded qbed files.
Each stage checks for existing output to avoid re-running completed steps.
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
import pandas as pd
import pybedtools
from tqdm import tqdm
from umi_tools import UMIClusterer
import numpy as np
import pyBigWig

# --- Script Version ---
VERSION = "1.3.0"


# --- Function Definitions ---
# Note: All helper functions (setup_logging, calculate_total_count_size_factors, etc.)
# remain the same as the previous version. They are included here for completeness.

def setup_logging(output_dir: Path) -> None:
    log_file = output_dir / "insertion_maps.log"
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] - %(message)s", handlers=[logging.FileHandler(log_file, mode="w"), logging.StreamHandler(sys.stdout)])
    logging.info(f"Logging initialized. Log file: {log_file}")

def calculate_total_count_size_factors(group_maps: Dict[str, pl.DataFrame]) -> Dict[str, float]:
    logging.info("Calculating Total Count size factors from generated raw counts...")
    library_sizes = {name: df["unique_count"].sum() for name, df in group_maps.items()}
    log_totals = np.log([count + 1 for count in library_sizes.values()])
    geo_mean_total = np.exp(np.mean(log_totals))
    if geo_mean_total <= 1:
        return {name: 1.0 for name in group_maps.keys()}
    factors = {name: total / geo_mean_total for name, total in library_sizes.items()}
    size_factor_df = pl.DataFrame({"group": list(factors.keys()), "size_factor": list(factors.values())})
    logging.info("Calculated Size Factors:\n%s", size_factor_df)
    return factors

def perform_validation_check(group_bgs: Dict[str, pl.DataFrame], final_results: Path):
    logging.info("=== Validation Check ===")
    if not final_results.is_file(): return
    fr = pl.read_csv(final_results, separator="\t", null_values="NA")
    for g, df in group_bgs.items():
        script_sum = int(df["unique_count"].sum()) if not df.is_empty() else 0
        pipeline_sum = int(fr.filter(pl.col("group_name") == g).group_by("consensus_peak_id").agg(pl.first("group_total_insertions_in_consensus_peak"))["group_total_insertions_in_consensus_peak"].sum())
        if script_sum == pipeline_sum:
            logging.info(f"✔ Group '{g}': Total counts MATCH! (Raw Count: {script_sum})")
        else:
            logging.warning(f"✖ Group '{g}': Total counts MISMATCH! script_raw_count={script_sum} pipeline_raw_count={pipeline_sum}")

def process_sample_insertions(task: tuple) -> Optional[pl.DataFrame]:
    sample_name, parquet_path, fragment_peaks_path, dist_thresh = task
    try:
        frags = pl.read_parquet(parquet_path)
        dedup_keys = [c for c in frags.columns if c not in ["sample_barcode_raw", "reads"]]
        frags = frags.group_by(dedup_keys).agg(pl.sum("reads").alias("reads"), pl.col("sample_barcode_raw").first().alias("sample_barcode_raw"))
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
            deduped_frags_df = grp.with_columns(pl.col("srt_bc").replace(rep).alias("srt_bc_rep")).filter(pl.col("srt_bc_rep").is_not_null()).group_by("chrom", "start", "end", "strand", "srt_bc_rep").agg(pl.sum("reads").alias("reads_in_molecule"))
            if deduped_frags_df.height < 5: continue
            summary = deduped_frags_df.group_by("srt_bc_rep").agg([pl.col("start").max().alias("max_start"), pl.col("end").min().alias("min_end"), pl.first("strand").alias("strand"), pl.first("chrom").alias("chrom")]).with_columns(pl.when(pl.col("strand") == "+").then(pl.col("min_end")).otherwise(pl.col("max_start")).alias("site")).select(["chrom", "site"])
            sites.append(summary)
        if not sites: return None
        final = pl.concat(sites)
        return final.group_by(["chrom", "site"]).agg(pl.count().alias("unique_count")).rename({"site": "start"}).with_columns((pl.col("start") + 1).alias("end")).select(["chrom", "start", "end", "unique_count"]).sort(["chrom", "start"])
    except Exception as e:
        logging.error(f"Error in sample {sample_name}: {e}", exc_info=True)
        return None

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
        logging.info(f"✅ Successfully created binned BigWig: {binned_bw_path}")
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

def main():
    """Main function to orchestrate the entire workflow."""
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    # --- Args ---
    p.add_argument("--output_dir_of_step_1", required=True)
    p.add_argument("--ins_map_output_dir", required=True)
    p.add_argument("--annotation_file", required=True)
    p.add_argument("--chrom_sizes", required=True)
    p.add_argument("--final_results_tsv", default=None)
    p.add_argument("--srt_bc_dist_threshold", type=int, default=1)
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--sum_window_size", type=int, default=50)
    args = p.parse_args()

    # --- Setup ---
    out = Path(args.ins_map_output_dir)
    out.mkdir(parents=True, exist_ok=True)
    setup_logging(out)
    logging.info(f"Insertion Mapper v{VERSION}")
    if not shutil.which("bedGraphToBigWig"):
        logging.error("`bedGraphToBigWig` not found in PATH."); sys.exit(1)
    
    base = Path(args.output_dir_of_step_1)
    sizes = Path(args.chrom_sizes)
    pybedtools_tmp_dir = out / "pybedtools_tmp"
    pybedtools_tmp_dir.mkdir(exist_ok=True, parents=True)
    pybedtools.helpers.set_tempdir(pybedtools_tmp_dir)

    dirs = {
        "sample_bg": out / "raw_unique_insertions_per_sample_bedgraph",
        "sample_bw": out / "raw_unique_insertions_per_sample_bigwig",
        "group_bg": out / "raw_unique_insertion_count_per_group_bedgraph",
        "group_bw": out / "raw_unique_insertion_count_per_group_bigwig",
        "norm_bg": out / "size_normalized_unique_insertions_per_group_bedgraph",
        "norm_bw": out / "size_normalized_unique_insertions_per_group_bigwig",
        "binned_sum_bw": out / f"binned_normalized_count_sum_bigwig_window_{args.sum_window_size}bp" if args.sum_window_size else None,
        "qbed": out / "raw_unique_insertion_count_per_group_qbed"
    }
    for d in dirs.values():
        if d: d.mkdir(exist_ok=True, parents=True)

    ann = pl.read_csv(args.annotation_file, separator='\t' if '\t' in open(args.annotation_file).readline() else ',', has_header=True)
    ann = ann.rename({col: col.strip() for col in ann.columns}).with_columns(pl.col("sample_name").str.strip_chars())
    group_to_samples: List[tuple[str, List[str]]] = ann.group_by("group_name").agg(pl.col("sample_name")).rows()
    
    # --- STAGE 1: Per-Sample Unique Insertion Mapping ---
    logging.info("--- STAGE 1: Per-Sample Insertion Mapping ---")
    sample_maps: Dict[str, pl.DataFrame] = {}
    collapsed = base / "collapsed_per_sample"
    fragment_peaks_dir = base / "fragment_peaks"
    
    all_samples = [f.stem for f in collapsed.glob("*.parquet")]
    expected_sample_bgs = {f"{s}.bedgraph" for s in all_samples}
    existing_sample_bgs = {f.name for f in dirs["sample_bg"].glob("*.bedgraph")}

    if not expected_sample_bgs.issubset(existing_sample_bgs):
        if not (collapsed.is_dir() and fragment_peaks_dir.is_dir()):
            logging.error("Missing required 'collapsed_per_sample/' or 'fragment_peaks/' directories."); sys.exit(1)
        
        logging.info(f"Starting per-sample insertion mapping with {args.workers} workers...")
        tasks = [(s, str(collapsed / f"{s}.parquet"), str(fragment_peaks_dir / f"{s}_fragment_peaks.parquet"), args.srt_bc_dist_threshold) for s in all_samples]
        with ProcessPoolExecutor(args.workers) as exe:
            futures = {exe.submit(process_sample_insertions, t): t[0] for t in tasks}
            for fut in tqdm(as_completed(futures), total=len(futures), desc="Mapping Insertions Per Sample"):
                s, df = futures[fut], fut.result()
                if df is not None and not df.is_empty():
                    sample_maps[s] = df
                    df.write_csv(dirs["sample_bg"] / f"{s}.bedgraph", separator="\t", include_header=False)
    else:
        logging.info("All per-sample bedgraphs found, loading from disk...")
        for s in tqdm(all_samples, desc="Loading sample maps"):
            df = pl.read_csv(dirs["sample_bg"] / f"{s}.bedgraph", separator="\t", has_header=False, new_columns=["chrom", "start", "end", "unique_count"])
            if not df.is_empty():
                sample_maps[s] = df
    
    # Generate BigWigs for samples (can be skipped if they exist)
    for s in all_samples:
        bg = dirs["sample_bg"] / f"{s}.bedgraph"
        bw = dirs["sample_bw"] / f"{s}.bw"
        if bg.exists() and not bw.exists():
            try:
                subprocess.run(["bedGraphToBigWig", str(bg), str(sizes), str(bw)], check=True, capture_output=True, text=True)
            except subprocess.CalledProcessError as e:
                logging.error(f"Failed to create BigWig for sample {s}: {e.stderr.strip()}")

    # --- STAGE 2: Aggregate Samples into Groups ---
    logging.info("--- STAGE 2: Aggregating Samples into Groups ---")
    group_maps: Dict[str, pl.DataFrame] = {}
    expected_group_bgs = {f"{g}.bedgraph" for g, _ in group_to_samples}
    existing_group_bgs = {f.name for f in dirs["group_bg"].glob("*.bedgraph")}

    if not expected_group_bgs.issubset(existing_group_bgs):
        if not sample_maps: logging.error("No sample data available for group aggregation."); sys.exit(1)
        logging.info("Generating group-level bedgraphs...")
        for group_name, sample_name_list in tqdm(group_to_samples, desc="Aggregating groups"):
            samp_list = [s for s in sample_name_list if s in sample_maps]
            if not samp_list: continue
            concat = pl.concat([sample_maps[s] for s in samp_list])
            agg = concat.group_by(["chrom", "start", "end"]).agg(pl.sum("unique_count")).sort(["chrom", "start"])
            group_maps[group_name] = agg
            agg.write_csv(dirs["group_bg"] / f"{group_name}.bedgraph", separator="\t", include_header=False)
    else:
        logging.info("All group-level bedgraphs found, loading from disk...")
        for group_name, _ in tqdm(group_to_samples, desc="Loading group maps"):
            df = pl.read_csv(dirs["group_bg"] / f"{group_name}.bedgraph", separator="\t", has_header=False, new_columns=["chrom", "start", "end", "unique_count"])
            group_maps[group_name] = df
    
    # Generate BigWigs for groups (can be skipped if they exist)
    for group_name, _ in group_to_samples:
        bg = dirs["group_bg"] / f"{group_name}.bedgraph"
        bw = dirs["group_bw"] / f"{group_name}.bw"
        if bg.exists() and not bw.exists():
             try:
                subprocess.run(["bedGraphToBigWig", str(bg), str(sizes), str(bw)], check=True, capture_output=True, text=True)
             except subprocess.CalledProcessError as e:
                logging.error(f"Failed to create BigWig for group {group_name}: {e.stderr.strip()}")

    # --- STAGE 3: Normalize Group Maps ---
    logging.info("--- STAGE 3: Normalizing Group Maps ---")
    if len(group_maps) > 1:
        expected_norm_bws = {f"{g}.bw" for g, _ in group_to_samples}
        existing_norm_bws = {f.name for f in dirs["norm_bw"].glob("*.bw")}
        if not expected_norm_bws.issubset(existing_norm_bws):
            logging.info("Generating normalized group-level tracks...")
            factors = calculate_total_count_size_factors(group_maps)
            for g, df in group_maps.items():
                fct = factors.get(g, 1.0)
                if fct > 0:
                    norm_df = df.with_columns((pl.col("unique_count") / fct).alias("normalized_count")).select(["chrom", "start", "end", "normalized_count"])
                    bg = dirs["norm_bg"] / f"{g}.bedgraph"
                    bw = dirs["norm_bw"] / f"{g}.bw"
                    norm_df.write_csv(bg, separator="\t", include_header=False)
                    try:
                        subprocess.run(["bedGraphToBigWig", str(bg), str(sizes), str(bw)], check=True, capture_output=True, text=True)
                    except subprocess.CalledProcessError as e:
                        logging.error(f"Failed to create normalized BigWig for {g}: {e.stderr.strip()}")
        else:
            logging.info("All normalized BigWig files found, skipping normalization.")
    
    final_tsv = Path(args.final_results_tsv) if args.final_results_tsv else base / "final_results.tsv"
    perform_validation_check(group_maps, final_tsv)

    # --- STAGE 4: Convert Group BedGraphs to Expanded Qbeds ---
    convert_bedgraphs_to_qbed(bedgraph_dir=dirs["group_bg"], qbed_dir=dirs["qbed"])

    # --- STAGE 5: Binned Sum Track Generation ---
    if args.sum_window_size and args.sum_window_size > 0:
        logging.info("--- STAGE 5: Binned Sum Track Generation ---")
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            for group_name, _ in group_to_samples:
                bw_path = dirs["norm_bw"] / f"{group_name}.bw"
                if bw_path.is_file():
                    executor.submit(generate_binned_sum_track, group_name, bw_path, args.sum_window_size, dirs["binned_sum_bw"])

    logging.info("✅ All processing complete.")


if __name__ == "__main__":
    main()
