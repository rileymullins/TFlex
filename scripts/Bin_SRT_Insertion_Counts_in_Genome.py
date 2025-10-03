#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
from multiprocessing import Pool, cpu_count
from functools import partial
import logging

# --- Configuration ---
BIN_SIZE = 200
N_CORES = 10

# --- Path Setup ---
qbed_dir = Path("<path to per replicate standard CCs processed data in qbed file format>")
bedgraph_dir = Path("<path to per replicate TFlex data in bedgraph file format>")
output_dir = Path("<output_dir>")
output_dir.mkdir(parents=True, exist_ok=True)

# --- Logging ---
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def qbed_to_bedgraph(qbed_path: Path):
    """
    FIXED: Reads the qbed file as a standard CSV with a header.
    """
    try:
        # The file is a standard CSV, so pandas can handle it with default settings.
        df = pd.read_csv(qbed_path)
        if df.empty:
            return None
        
        # Rename columns to be consistent (lowercase 'chr', 'start', 'end').
        df = df.rename(columns={'Chr': 'chr', 'Start': 'start', 'End': 'end'})
        
        # Count unique insertions per coordinate to create a 'count' column.
        bg = df.groupby(['chr', 'start', 'end']).size().reset_index(name='count')
        return bg
    except Exception as e:
        logging.error(f"Failed to convert qbed to bedgraph: {qbed_path}, error: {e}")
        return None

def bin_bedgraph_df(df: pd.DataFrame, sample_name: str, bin_size: int):
    """
    Efficiently bins a bedgraph-like DataFrame using a vectorized operation.
    """
    total_input = df['count'].sum()
    df['bin_start'] = (df['start'] // bin_size) * bin_size
    binned = df.groupby(['chr', 'bin_start'])['count'].sum().reset_index()
    binned['bin_end'] = binned['bin_start'] + bin_size
    binned['sample'] = sample_name
    binned = binned[['chr', 'bin_start', 'bin_end', 'count', 'sample']]

    total_binned = binned['count'].sum()
    if total_input == total_binned:
        logging.info(f"{sample_name}: ✅ Counts match ({total_input})")
    else:
        logging.warning(f"{sample_name}: ⚠️ Count mismatch! Input={total_input}, Binned={total_binned}")
    return binned

def process_file(file_path: Path, bin_size: int):
    """Main worker function to process a single file."""
    sample_name = file_path.stem
    logging.info(f"Processing {sample_name}...")
    
    df = None
    if file_path.suffix.lower() == '.qbed':
        df = qbed_to_bedgraph(file_path)
    else:  # Assumes .bedgraph
        # FIX: Use a raw string r'\t' to avoid the SyntaxWarning.
        df = pd.read_csv(file_path, sep=r'\t', header=None, names=['chr', 'start', 'end', 'count'], engine='python')

    if df is None or df.empty:
        logging.warning(f"No valid data to process for {sample_name}")
        return None

    binned_df = bin_bedgraph_df(df, sample_name, bin_size)

    if binned_df is None or binned_df.empty:
        logging.warning(f"No non-zero bins for {sample_name}")
        return None

    out_file = output_dir / f"{sample_name}_binned{bin_size}bp.bedgraph"
    binned_df[['chr', 'bin_start', 'bin_end', 'count']].to_csv(out_file, sep='\t', index=False, header=False)
    
    return binned_df

if __name__ == "__main__":
    all_files = list(qbed_dir.glob("*.qbed")) + list(bedgraph_dir.glob("*.bedgraph"))
    if not all_files:
        logging.warning("No .qbed or .bedgraph files found. Exiting.")
    else:
        n_cores_to_use = max(1, cpu_count() - 2)
        logging.info(f"Starting parallel processing with {n_cores_to_use} cores...")
        
        worker_func = partial(process_file, bin_size=BIN_SIZE)
        with Pool(n_cores_to_use) as pool:
            results = pool.map(worker_func, all_files)

        successful_results = [df for df in results if df is not None]
        if successful_results:
            logging.info("Combining all samples into a single file...")
            combined_df = pd.concat(successful_results, ignore_index=True)
            combined_file = output_dir / f"combined_binned_{BIN_SIZE}bp_all_samples.tsv"
            combined_df.to_csv(combined_file, sep='\t', index=False)
            logging.info(f"✅ Successfully written combined file: {combined_file}")
        else:
            logging.warning("No files were successfully processed.")
