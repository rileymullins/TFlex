# multiplexed_srt_tf_mapping

Identify peaks from multiplexed self-reporting transposon data for identifying transcription factor binding sites.

---

## Overview of method
This method is capable of antibody-independent semi-multiplxed TF binding site experiments. That is, after initial transfection of a TF-tranposase fusion and a pool of self-reporting transposons containing many SRT barcodes and a unique sample barcode, all samples can be combined. 

This method leverages the self-reporting transposon technology. It is "self-reporting" because it contains a promoter that generates transcripts containing the junction of the terminal repeat and downstream genomic DNA. Thus, the sample barcode and SRT barcode are encoded into the genomic DNA in the transposon, which is then transcribed. Thus, the approximate TF binding site can be identifed through RNA sequencing.

After pooling all samples, the workflow is RNA extraction -> cDNA generation -> amplification of the self-reporting transpsons -> Illumina Nextera XT kit -> sequencing.
![Picture1](https://github.com/user-attachments/assets/e02c5bdd-cfae-4988-9c55-deff1e93ea60)

---

## Description of data processing script
This pipeline processes multiplexed self-reporting transposon data to identify transcription factor (TF) binding sites. It supports high-throughput experiments with many barcoded samples, corrects technical errors in barcode assignment, collapses redundant reads, calls peaks per sample, and generates a consensus set of reproducible peaks across all experimental conditions.
![Picture1](https://github.com/user-attachments/assets/dd45c0cc-8132-4d65-af3a-1b24befb073c)

## Pipeline steps

1. **Input loading and barcode parsing:**  
   Reads all qbed/bed files and parses the `name` field into `library_name`, `sample_barcode_raw`, and `srt_bc`.
2. **Barcode correction:**  
   Corrects sample barcodes against the annotation whitelist, using Hamming distance up to the `--sample_barcode_dist_threshold`.
3. **Annotation:**  
   Joins corrected reads with the annotation file, assigning `sample_name` and `group_name`.
4. **Per-sample partitioning:**  
   Collapses any duplicate rows for each sample, summing the reads, and saved ass a `.parquet` file for efficient downstream processing.
5. **Peak calling:**  
   For each sample, calls peaks with pycallingcards, then refines peak boundaries using strand-aware logic of the position per SRT barcode that is the most proximal to the junction of the transposon and genomic DNA.
6. **Consensus peak generation:**  
   Merges sample-specific peaks to produce consensus peaks across all samples.
7. **Read-to-peak mapping:**  
  Intersects all sample peak sets with the consensus peaks to map which sample peaks were merged in the consensus peak.
8. **Output aggregation:**  
  Each row is a unique consensus peak and sample_name pair populated with the per-sample peak statistics (e.g., fragment-based peak coordinates, SRT barcode-based peak coordinates, total reads, total fragments, and unique insertions within each of those peaks). The attributes of all samples in the same group are summed per consensus peak to produces the per-group statistics (`final_results.tsv`).
   A `column_definitions.tsv` file describing all output columns is also saved in the output directory.

---
### Supported input file formats (a .qbed file is generated from the data alignment pipeline but a .bed could be used)

- `.qbed`, `.bed`, `.qbed.gz`, `.bed.gz`
- Each file must contain columns in the order:  
  `chrom`, `start`, `end`, `reads`, `strand`, `name`
- The `name` field should be formatted as:  
  `library_name/sample_barcode/srt_barcode`
  
### Example input file row

```
chr    start    end    reads    strand        name
chr1   12345	  12359	 100	    +	        [library_name]/[sample_barcode]/[srt_barcode]
```

## Interpretation of 'strand', 'start', 'end', and read orientation in relation to the tranpsoson
- The R2 read end is reported in the input file.
- The strand in the qbed refers the strand that the **transposon** was inserted into, not the strand of the read.
- By convention of the alignemnt pipeline, **'end' coordinate** for a **'+'** strand row is the true end of the read, and the **'start' coordinate** for a **'-'** strand row is the true end of the read.
- Reads that hit the transposon are truncated to the last base prior to beginning of the transposon.

- Summary:
  - For a **'+'** strand row in the qbed:
    - **Transposon** is inserted on the **'+'** strand.
    - R2 read is moving  **5' <- 3'**.
    - **'end' coordinate** is the true end of the read.
    - **Smallest 'end' coordinate** (most upstream) is the closest to the transposon.

  - For a **'-'** strand row in the qbed:
    - **Transposon** is inserted on the **'-'** strand.
    - R2 read is moving  **5' -> 3'**.
    - **'start' coordinate** is the true end of the read.
    - **Largest start coordinate** (most downstream) is the closest to the transposon.
      
## Key logic for converting the fragment-based peaks to SRT barcode-based peak boundaries
- **Fragment-based peaks:** Peaks called from the initial input qbed file where every row represents one uniquely fragmented molecule for which the read ended at a different position.
- **SRT barcode-based peaks:** Refined fragment-based peaks where the start and end coordinates are defined based on the unique SRT barcode positions within the peak.
  
  - Peak start = the **smallest** value between the minimum 'end' coordinate of '+' strand fragments and maximum 'start' coordinate of '-' strand fragments.
    - This represents the most **upstream** value that is the most proximal to the transposon between the '+' and '-' strand fragments.
  - Peak end = the **largest** value between the minimum 'end' coordinate of '+' strand fragments and maximum 'start' coordinate of '-' strand fragments.
    - This represents the most **downstream** value that is the most proximal to the transposon between the '+' and '-' strand fragments.
  
- Note that in both instances, the minimum 'end' coordinate of '+' strand fragments and maximum 'start' coordinate of '-' strand fragments are used.
  - This is because these value represent the most proximal position of all fragments of a given SRT barcode to the transposon.
  - It is the best approximation of the true transposon insertion site if the read does not make contact with the transposon.


### Features

- **Flexible input:** Accepts multiple `.qbed`, `.bed`, or gzipped equivalents from a directory.
- **Barcode correction:** Corrects sample barcodes against a whitelist (from the annotation file) using defined Hamming distance cutoff.
- **Customizable:** Uses a user-supplied annotation file to map barcodes to sample and group names.
- **Fragment deduplication:** Collapses reads within each peak that are identical in coordinates and SRT barcodes using UMI-tools clustering-based correction.
- **Parallel peak calling:** Efficiently calls peaks per sample using pycallingcards, leveraging multiprocessing.
- **Consensus peak construction:** Merges sample-specific peaks to produce a unified, non-redundant set of consensus peaks across the dataset.
- **Comprehensive reporting:** Outputs per-sample and per-group statistics, final peak sets, and detailed column definitions.
- **Robust logging:** Logs all steps and errors for traceability.
- **Memory efficient:** Capable of running on most laptops, but will need to monitor memory usage and lower the number of workers as needed.

---

## Usage
**Default parameters command line**
```bash
python ccaller_to_assign_reproducible_peaks_with_full_pan_barcode_consensus_peaks_and_umi_based_peak_boundaries.py \
  --input_dir /path/to/qbed_files \
  --output_dir /path/to/results \
  --annotation_file /path/to/annotation_file.tsv \
  --workers 8 \
  --sample_barcode_dist_threshold 2 \
  --srt_bc_dist_threshold 1 \
  --min_rows_threshold 50000 \
  --method CCcaller \
  --reference hg38 \
  --pvalue_cutoff 0.01 \
  --pvalue_adj_cutoff 0.01 \
  --min_insertions 5 \
  --minlen 0 \
  --extend 200 \
  --maxbetween 150 \
  --minnum 0 \
  --lam_win_size 1000000
```
min_insertions ,minlen ,extend, maxbetween, minnum, lam_win_size all refer to the pycallingcards peak caller settings. Full description available [here](https://pycallingcards.readthedocs.io/en/latest/api/reference/pycallingcards.preprocessing.call_peaks.html)

---

### Required arguments

- `--input_dir`  
  Directory containing input `.qbed`, `.bed`, `.qbed.gz`, or `.bed.gz` files. Each file should be in standard qbed/bed format.
- `--output_dir`  
  Directory where output files will be written.
- `--annotation_file`  
  Path to the annotation file (see format below).

### Optional arguments

- `--workers`  
  Number of parallel worker processes for peak calling (default: 10). One worker will be assigned a sample_name and perform all functions on that sample_name.
- `--sample_barcode_dist_threshold`  
  Maximum Hamming distance allowed for correcting sample barcodes (default: 2). A value of 2 means that only sample barcodes that have ≤ 2 mismatches from the sample barcode are assigned to the whitelisted sample barcode.
- `--srt_bc_dist_threshold`  
  Maximum Hamming distance for SRT-barcode clustering (default: 1). A value of 1 means that only SRT-barcodes with ≤ 1 mismatches from each other will be grouped.
- Additional parameters for fragment-based peak calling (see script source for defaults and descriptions):
  - `--pvalue_cutoff`, `--pvalue_adj_cutoff`, `--min_insertions`, `--minlen`, `--extend`, `--maxbetween`, `--minnum`, `--lam_win_size`
  - It is recommended to leave these parameters at the default values.
  - The purpose of these parameters is to define all regions of concentrated fragments for subsequent SRT-barcode-based peak boundary definitions and donwstream analysis.
  - Fragment peaks with <5 fragments following SRT barcode correction are considered noise and discarded from the analysis.

---

## Annotation file format

This file should be a `.csv` or `.tsv` (tab-separated values) format. The columns must appear **in this order**:  
**library_name**, **sample_barcode**, **sample_name**, **group_name**.

**Example annotation file (tab-separated):**

```
library_name	sample_barcode	sample_name	group_name
Library_A	AAGGCAGACG	Rep1_HyPBase	HyPBase
Library_A	TACCGCTGAC	Rep1_TF_A	TF_A
Library_A	AAGATTAGAC	Rep1_TF_B	TF_B
Library_A	AGTATGACCG	Rep1_TF_C	TF_C
Library_A	TACTAGGTTG	Rep1_TF_D	TF_D
Library_A	CTCTGGATTA	Rep1_TF_E	TF_E
Library_A	AGTTAACGAT	Rep1_TF_F	TF_F
Library_A	GAGACACGTC	Rep1_TF_G	TF_G
Library_A	GTAATGCTTC	Rep1_TF_H	TF_H
Library_A	TTAACGATCG	Rep1_TF_I	TF_I
Library_B	AAGGCAGACG	Rep2_HyPBase	HyPBase
Library_B	TACCGCTGAC	Rep2_TF_A	TF_A
Library_B	AAGATTAGAC	Rep2_TF_B	TF_B
Library_B	AGTATGACCG	Rep2_TF_C	TF_C
Library_B	TACTAGGTTG	Rep2_TF_D	TF_D
Library_B	CTCTGGATTA	Rep2_TF_E	TF_E
Library_B	AGTTAACGAT	Rep2_TF_F	TF_F
Library_B	GAGACACGTC	Rep2_TF_G	TF_G
Library_B	GTAATGCTTC	Rep2_TF_H	TF_H
Library_B	TTAACGATCG	Rep2_TF_I	TF_I
```

- **All columns must be present, and tab-separated.**
- **All `sample_name` values must be unique.**

---


## Output files

- `final_results.tsv`  
  Tab-separated table containing all peak and sample statistics. See variable definitions below.
- `column_definitions.tsv`  
  Tab-separated table describing all columns in `final_results.tsv`.
- `pipeline.log`  
  Log file of all steps, warnings, and errors.
- `collapsed_per_sample/*.parquet`  
  Per-sample data files (intermediate).
- `pybedtools_tmp/`  
  Temporary files for bedtools operations (auto-cleaned at completion).

---

## Variable names

| Variable | Description |
|---|---|
| consensus_peak_id | chrom:consensus_peak_start-consensus_peak_end |
| chrom | chromosome |
| consensus_peak_start | start coordinate of final merged, pan-dataset consensus peak |
| consensus_peak_end | end coordinate of final merged, pan-dataset consensus peak |
| consensus_peak_width | width (end – start) of final merged, pan-dataset consensus peak |
| num_samples_in_consensus_peak | number of unique sample_name values in this consensus peak |
| num_groups_in_consensus_peak | number of unique group_name values in this consensus peak |
| sample_name | sample identifier (from annotation file). **Must be unique.** This is the peak calling unit; could be a replicate (e.g., replicate-1_TF-A), a unique time point (e.g., timepoint1_replicate1_TF-A), or a unique experimental condition (e.g., drugX_replicate1_TF-A). |
| group_name | group identifier (from annotation file). Used to aggregate stats across samples belonging to the same broader group. Examples: "TF_A", "timepoint1_TF-A", or "drugX_TF-A". |
| fragment_peak_start_for_sample | fragment-based peak start coordinate for this sample_name |
| fragment_peak_end_for_sample | fragment-based peak end coordinate for this sample_name |
| fragment_peak_width_for_sample | fragment-based peak width (end – start) for this sample_name |
| srt_bc_peak_start_for_sample | SRT-barcode-based peak start coordinate for this sample_name |
| srt_bc_peak_end_for_sample | SRT-barcode-based peak end coordinate for this sample_name |
| srt_bc_peak_width_for_sample | SRT-barcode-based peak width (end – start) for this sample_name |
| sample_total_reads_in_consensus_peak | total read count in consensus peak for this sample_name |
| sample_total_fragments_in_consensus_peak | total fragment count (merged qbed rows after SRT-barcode correction and deduplication) in consensus peak for this sample_name |
| sample_total_insertions_in_consensus_peak | total unique insertions (unique SRT barcodes) in consensus peak for this sample_name |
| group_total_reads_in_consensus_peak | sum of sample_total_reads_in_consensus_peak across all sample_name values in this group_name |
| group_total_fragments_in_consensus_peak | sum of sample_total_fragments_in_consensus_peak across all sample_name values in this group_name |
| group_total_insertions_in_consensus_peak | sum of sample_total_insertions_in_consensus_peak across all sample_name values in this group_name |

---

## Dependencies

The following packages and versions were used to create and test this script using **Python 3.12.10**:

| Package         | Version    | Purpose                        |
|-----------------|------------|--------------------------------|
| [polars](https://docs.pola.rs/user-guide/installation/)          | 1.31.0     | Fast DataFrames manipulation                |
| [pandas](https://pandas.pydata.org/docs/getting_started/install.html)          | 2.2.3      | DataFrames, pybedtools input   |
| [pybedtools](https://daler.github.io/pybedtools/main.html)      | 0.12.0     | Genomic intervals intersections             |
| [pycallingcards](https://pycallingcards.readthedocs.io/en/latest/#)  | 1.0.0      | Fragment-based peak calling                   |
| [UMI-tools](https://umi-tools.readthedocs.io/en/latest/INSTALL.html)       | 1.1.6      | SRT barcode clustering to deduplicate within peak             |
| [tqdm](https://pypi.org/project/tqdm/)            | 4.67.1     | Progress bars                  |
| [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html)        | 2.31.1     | (external, required for pybedtools) |

### Installation

Install all required Python packages with pip:
```bash
pip install polars==1.31.0 pandas==2.2.3 pybedtools==0.12.0 pycallingcards==1.0.0 umi_tools==1.1.6 tqdm==4.67.1
```
> **Note:**  
> pybedtools requires bedtools be installed on your system and available in your `$PATH`.  
> Install with conda:
> ```bash
> conda install -c bioconda bedtools
> ```
> Or with apt (on Ubuntu):
> ```bash
> sudo apt-get install bedtools
> ```

---

## Checking your installed versions

To check package versions, run the following in the terminal:

```bash
python -c "import polars, pandas, pybedtools, tqdm, umi_tools, pycallingcards as cc; print('polars:', polars.__version__); print('pandas:', pandas.__version__); print('pybedtools:', pybedtools.__version__); print('tqdm:', tqdm.__version__); print('umi_tools:', umi_tools.__version__); print('pycallingcards:', cc.__version__)"
bedtools --version
```

---

## Options

- All pycallingcards peak calling parameters can be overridden via the command line.
- The pipeline can be parallelized across as many CPUs as available using `--workers`.
- Barcode and SRT clustering stringency can be tuned via `--sample_barcode_dist_threshold` and `--srt_bc_dist_threshold`.

---


## Contact

For questions, bugs, or feature requests, please open an issue or contact the repository maintainer.
