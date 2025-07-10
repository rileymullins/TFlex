# Multiplexed TF Mapping with SRTs - A Bioinformatics Pipeline

This repository contains a series of scripts to perform an end-to-end analysis for multiplexed transcription factor (TF) mapping using self-reporting transposons (SRTs). The pipeline processes raw sequencing data (in qbed/bed format), calls peaks, performs differential analysis, and annotates the results to identify TF-specific binding sites and enriched motifs.

## Table of Contents

- [Pipeline Overview](#pipeline-overview)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
  - [Step 1: Raw Data to Insertion Maps](#step-1-raw-data-to-insertion-maps)
  - [Step 2: Peak Calling with SPAN](#step-2-peak-calling-with-span)
  - [Step 3: Generate Peak-by-Group Count Matrices](#step-3-generate-peak-by-group-count-matrices)
  - [Step 4: Differential Analysis, Annotation, and Motif Discovery](#step-4-differential-analysis-annotation-and-motif-discovery)
- [Output Directory Structure](#output-directory-structure)
- [Dependencies](#dependencies)

## Pipeline Overview

The workflow is divided into four main stages, each handled by a dedicated script:

1. **`step_1_raw_qbed_to_insertion_maps.py`**: This script is the entry point of the pipeline. It takes raw qbed files, performs barcode correction, annotates reads based on a sample sheet, and generates 1-bp unique insertion site maps. It outputs raw and normalized signal tracks (bedGraph/BigWig) for each sample and experimental group.

2. **`step_2_run_span_peak_caller.sh`**: Using the group-level qbed files generated in Step 1, this script runs the SPAN peak caller to identify regions of significant insertion enrichment for each experimental group.

3. **`step_3_generate_peak_by_group_count_matrices.py`**: This script takes the SPAN peaks from Step 2 and the insertion data from Step 1. It resizes peaks to a minimum width and generates count matrices, quantifying the number of insertions from every group that fall within each group's peak set.

4. **`step_4_DESeq2_Diff_Peaks_HOMER-Annotations_and_Motifs.R`**: The final step uses the count matrices from Step 3 to perform differential binding analysis with DESeq2. It identifies group-specific peaks, annotates them using HOMER, performs motif discovery, and generates publication-quality plots and heatmaps.

## Prerequisites

### Software

- Python 3.8+
- R 4.0+
- **SPAN Peak Caller**: The JAR file must be downloaded and accessible.
- **HOMER**: Must be installed and its bin directory added to your system's PATH.
- **UCSC Kent Utilities**: Specifically, `bedGraphToBigWig` must be in your system's PATH.

### Python Libraries

- polars
- pandas
- numpy
- pybedtools
- tqdm
- umi_tools
- pycallingcards
- pyBigWig

### R Libraries

- DESeq2
- ComplexHeatmap
- tidyverse (and its components like dplyr, purrr, ggplot2)
- data.table
- RColorBrewer
- circlize
- viridis
- janitor

## Installation

1. **Clone the repository:**
   ```bash
   git clone <your-repo-url>
   cd <your-repo-directory>
   ```

2. **Install Python dependencies:**
   ```bash
   pip install -r requirements.txt
   ```
   *Note: You may need to create a requirements.txt file based on the library versions listed in the scripts.*

3. **Install R dependencies:**
   Open an R session and run:
   ```r
   # Install from CRAN
   install.packages(c("tidyverse", "data.table", "RColorBrewer", "circlize", "viridis", "janitor", "parallel"))

   # Install from Bioconductor
   if (!requireNamespace("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   BiocManager::install(c("DESeq2", "ComplexHeatmap", "GenomicRanges"))
   ```

4. **Setup External Tools:**
   - Ensure the `span-*.jar` file is downloaded.
   - Install HOMER.
   - Install UCSC Kent Utilities.

## Usage

### Step 1: Raw Data to Insertion Maps

This script processes raw qbed files to generate unique insertion maps and signal tracks.

**Script:** `step_1_raw_qbed_to_insertion_maps.py`

**Arguments:**

- `--input_dir` (required): Directory containing input qbed/bed files.
- `--output_dir` (required): Directory to store all output files.
- `--annotation_file` (required): Path to a TSV/CSV file with columns: `library_name`, `sample_barcode`, `sample_name`, `group_name`.
- `--chrom_sizes` (required): Path to a chromosome sizes file (e.g., `hg38.chrom.sizes`).
- `--workers` (optional, default: 10): Number of parallel workers.
- `--srt_bc_dist_threshold` (optional, default: 1): Max Hamming distance for SRT-BC (UMI) clustering.
- `--sample_barcode_dist_threshold` (optional, default: 2): Max Hamming distance for sample barcode correction.
- `--min_rows_threshold` (optional, default: 50000): Minimum number of fragments for a sample to be processed.
- `--sum_window_size` (optional, default: 50): Window size (bp) for binned summary BigWig tracks.
- **Peak Calling Parameters:** Several arguments like `--pvalue_cutoff`, `--extend`, etc., are available to configure the pycallingcards peak caller.

**Example Command:**
```bash
python step_1_raw_qbed_to_insertion_maps.py \
  --input_dir /path/to/raw_qbeds \
  --output_dir /path/to/step1_output \
  --annotation_file /path/to/samples.tsv \
  --chrom_sizes /path/to/hg38.chrom.sizes \
  --workers 16
```

### Step 2: Peak Calling with SPAN

This script runs the SPAN peak caller on the group-level qbed files generated in Step 1.

**Script:** `step_2_run_span_peak_caller.sh`

**Configuration:**
You must edit the variables at the top of the script:

- `SPAN_JAR`: Path to the `span-*.jar` file.
- `INPUT_DIR`: Path to the `raw_unique_insertion_count_per_group_qbed` directory from Step 1.
- `OUTPUT_DIR`: Directory to save the output peak files.
- `CHROM_SIZES`: Path to your chromosome sizes file.
- `JAVA_MEM` (optional): Memory to allocate to Java (e.g., `-Xmx16G`).
- `BIN_SIZE` (optional, default: 50): Bin size for SPAN analysis.
- `FDR` (optional, default: 0.01): False Discovery Rate cutoff.
- `MAX_JOBS` (optional, default: 2): Maximum number of parallel SPAN jobs.

**Example Command:**
```bash
# After configuring the script internally
bash step_2_run_span_peak_caller.sh
```

### Step 3: Generate Peak-by-Group Count Matrices

This script intersects the SPAN peaks with insertion data to create count matrices.

**Script:** `step_3_generate_peak_by_group_count_matrices.py`

**Arguments:**

- `--peak_dir` (required): Directory containing pre-called SPAN peak files.
- `--peak_suffix` (required): Suffix of the peak files to use (e.g., `_b50.span.peak`).
- `--bedgraph_dir` (required): Directory with per-group insertion count bedGraph files.
- `--output_dir` (required): Directory to save output files.
- `--annotation_file` (required): Annotation file mapping samples to groups.
- `--workers` (optional, default: 12): Number of parallel worker processes.

**Example Command:**
```bash
python step_3_generate_peak_by_group_count_matrices.py \
  --peak_dir /path/to/step2_output \
  --peak_suffix _b50.span.peak \
  --bedgraph_dir /path/to/step1_output/bedgraph \
  --output_dir /path/to/step3_output \
  --annotation_file /path/to/samples.tsv \
  --workers 16
```

### Step 4: Differential Analysis, Annotation, and Motif Discovery

*Note: This section appears to be incomplete in the original document. The script `step_4_DESeq2_Diff_Peaks_HOMER-Annotations_and_Motifs.R` is mentioned but detailed usage instructions are not provided.*

## Output Directory Structure

*To be documented based on the actual output structure of each step.*

## Dependencies

The complete list of dependencies is provided in the [Prerequisites](#prerequisites) section above. Ensure all software and libraries are properly installed before running the pipeline.
