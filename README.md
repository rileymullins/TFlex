# Processing Multiplex SRT Sequencing Data to TF Binding Sites and Visualization Files

This repository contains a bioinformatics pipeline for processing multiplexed self-reporting transposon (SRT) sequencing data. The primary goal is to identify putative transcription factor (TF) binding sites (peaks) from raw sequencing reads and generate standard visualization files (bedGraph, BigWig) for genomic analysis.

---

## Table of Contents

- [Overview of the Method](#overview-of-the-method)
- [Script 1: `multiplex_srt_seq_to_tf_binding.py`](#script-1-multiplex_srt_seq_to_tf_bindingpy)
  - [Pipeline Steps](#pipeline-steps)
  - [Usage](#usage-of-multiplex_srt_seq_to_tf_bindingpy)
  - [Input Formats](#input-formats)
  - [Output Files](#output-files-of-multiplex_srt_seq_to_tf_bindingpy)
  - [Output Column Definitions](#output-column-definitions)
- [Script 2: `generate_insertion_maps.py`](#script-2-generate_insertion_mapspy)
  - [Workflow](#workflow)
  - [Usage](#usage-of-generate_insertion_mapspy)
  - [Output Files](#output-files-of-generate_insertion_mapspy)
- [Key Data Processing Details](#key-data-processing-details)
  - [Read Orientation](#strand-start-end-and-read-orientation)
  - [Peak Refinement Logic](#key-logic-for-converting-fragment-based-peaks-to-srt-barcode-based-peaks)

---

## Overview of the Method

This method leverages self-reporting transposon (SRT) technology. An SRT is "self-reporting" because it contains a promoter that drives the transcription of RNA containing the junction of the transposon's terminal repeat and the downstream genomic DNA. Sequencing this RNA precisely maps the transposon insertion site, which serves as a proxy for TF binding.

The method is antibody-independent and semi-multiplexed, as both the sample barcode and the SRT barcode are encoded in the transposon's DNA.
 
 ![Picture1](https://github.com/user-attachments/assets/46acac4e-7843-49a7-9e7e-2a5b616febc6)

**Workflow Summary:**

1.  **Transfection:** Samples are transfected with a large, diverse pool of plasmids. Each pool contains a unique sample barcode but many different SRT barcodes.
2.  **Pooling:** All samples, each with its unique sample barcode, are pooled after transfection.
3.  **Library Prep & Sequencing:** The workflow proceeds through RNA extraction, cDNA synthesis, amplification, and finally, library preparation with the Illumina Nextera XT kit for sequencing.
4.  **Demultiplexing:** Post-sequencing, reads are assigned to their original samples based on the sample barcode.
5.  **Mapping:** The unique transposon insertion sites are identified by their genomic coordinates and the SRT barcode, representing the location and signal of a TF binding event.

---

## Script 1: `multiplex_srt_seq_to_tf_binding.py`

This script is the core of the pipeline. It converts aligned SRT sequencing data into a set of putative TF binding sites. The primary output is a count matrix suitable for differential binding analysis with tools like DESeq2, allowing for the identification of binding sites enriched over a control (e.g., unfused transposase) and differential peaks between TFs.

 ![Picture1](https://github.com/user-attachments/assets/dd45c0cc-8132-4d65-af3a-1b24befb073c)

### Pipeline Steps

#### Phase 1: Loading, Correcting, and Annotating All Reads

1.  **Input Loading & Barcode Parsing**: Reads all input `qbed`/`bed` files. The `name` field is parsed to extract `library_name`, `sample_barcode_raw`, and `srt_bc`. To manage memory, files are processed sequentially, and intermediate results are partitioned by `sample_name` into temporary `.parquet` files.
2.  **Barcode Correction**: Corrects raw sample barcodes against a whitelist using Hamming distance (configurable with `--sample_barcode_dist_threshold`).
3.  **Annotation**: Joins the corrected data with the annotation file to assign `sample_name` and `group_name`.
4.  **Per-Sample Partitioning**: Saves the processed data for each `sample_name` into a dedicated `.parquet` file (`<output_dir>/collapsed_per_sample/`) for efficient downstream processing.

#### Phase 2: Calling Peaks on Each Sample's Data

5.  **Fragment-based Peak Calling**: For each sample, merges identical fragments and calls peaks using `pycallingcards`. This identifies genomic regions with significant transposon insertion activity. Per-sample peaks are saved in `<output_dir>/fragment_peaks/`.
6.  **Fragment-to-Peak Mapping**: Intersects the deduplicated fragments with their corresponding peak set to map which fragments belong to which peak.
7.  **SRT Barcode-based Peak Refinement**: Refines the initial peak boundaries using strand-aware logic to pinpoint the most likely insertion site for each unique SRT barcode. This step also counts the unique insertions (unique SRT barcodes) within the peak, which is the primary signal for TF binding. Refined peaks are saved in `<output_dir>/sample_peaks/`.

#### Phase 3: Generating Consensus Peaks and Final Output

8.  **Consensus Peak Generation**: Merges all sample-specific peaks to create a unified set of consensus peaks across all samples and groups.
9.  **Sample-to-Consensus Peak Mapping**: Intersects each sample's peak set with the consensus peak set.
10. **Output Aggregation**: Generates the final `final_results.tsv` file. Each row represents a consensus peak/sample pair, populated with comprehensive statistics (read counts, fragment counts, insertion counts). Stats are also aggregated at the group level. A `column_definitions.tsv` file is also generated.

### Usage of `multiplex_srt_seq_to_tf_binding.py`

#### Dependencies

This script was tested with **Python 3.12.10**.

| Package                                                                  | Version | Purpose                               |
| ------------------------------------------------------------------------ | ------- | ------------------------------------- |
| [polars](https://docs.pola.rs/user-guide/installation/)                  | 1.31.0  | Fast DataFrame manipulation           |
| [pandas](https://pandas.pydata.org/docs/getting_started/install.html)    | 2.2.3   | DataFrames, pybedtools input          |
| [pybedtools](https://daler.github.io/pybedtools/main.html)               | 0.12.0  | Genomic interval intersections        |
| [pycallingcards](https://pycallingcards.readthedocs.io/en/latest/)       | 1.0.0   | Fragment-based peak calling         |
| [UMI-tools](https://umi-tools.readthedocs.io/en/latest/INSTALL.html)     | 1.1.6   | SRT barcode clustering              |
| [tqdm](https://pypi.org/project/tqdm/)                                   | 4.67.1  | Progress bars                         |
| [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html) | 2.31.1  | (External, required by pybedtools)    |

*Other package versions may be compatible.*

#### Installation

1.  **Install Python packages:**
    ```bash
    pip install polars==1.31.0 pandas==2.2.3 pybedtools==0.12.0 pycallingcards==1.0.0 umi_tools==1.1.6 tqdm==4.67.1
    ```

2.  **Install `bedtools`:** `pybedtools` is a wrapper and requires `bedtools` to be installed and in your system's `$PATH`.
    * **With conda (recommended):**
        ```bash
        conda install -c bioconda bedtools
        ```
    * **With apt (Ubuntu/Debian):**
        ```bash
        sudo apt-get install bedtools
        ```

3.  **Verify installation:**
    ```bash
    # Check Python packages
    python -c "import polars, pandas, pybedtools, tqdm, umi_tools, pycallingcards as cc; print(f'polars: {polars.__version__}\npandas: {pandas.__version__}\npybedtools: {pybedtools.__version__}\ntqdm: {tqdm.__version__}\numi_tools: {umi_tools.__version__}\npycallingcards: {cc.__version__}')"

    # Check bedtools
    bedtools --version
    ```

#### Command Line Example

```bash
python multiplex_srt_seq_to_tf_binding.py \
    --input_dir /path/to/qbed_files \
    --output_dir /path/to/results \
    --annotation_file /path/to/annotation_file.tsv \
    --workers 10 \
    --sample_barcode_dist_threshold 2 \
    --srt_bc_dist_threshold 1 \
    --min_rows_threshold 50000 \
    --reference hg38
```
*Note: Parameters for the `pycallingcards` peak caller (`min_insertions`, `extend`, `maxbetween`, etc.) can also be set. See the [pycallingcards documentation](https://pycallingcards.readthedocs.io/en/latest/api/reference/pycallingcards.preprocessing.call_peaks.html) for details.*

#### Arguments

**Required:**
* `--input_dir`: Directory with input `.qbed`, `.bed`, `.qbed.gz`, or `.bed.gz` files.
* `--output_dir`: Directory where output files will be written.
* `--annotation_file`: Path to the annotation file.

**Optional:**
* `--workers`: Number of parallel worker processes for peak calling (default: `10`).
* `--sample_barcode_dist_threshold`: Max Hamming distance for correcting sample barcodes (default: `2`).
* `--srt_bc_dist_threshold`: Max Hamming distance for SRT barcode clustering (default: `1`).
* `--min_rows_threshold`: Minimum fragments required for a sample to be processed (default: `50000`).
* `--reference`: Reference genome (e.g., `hg38`, `mm10`). Required by `pycallingcards`.

### Input Formats

#### Qbed/Bed File Format
-   **No header** row.
-   Supported extensions: `.qbed`, `.bed`, `.qbed.gz`, `.bed.gz`.
-   Columns must be: `chrom`, `start`, `end`, `reads`, `strand`, `name`.
-   The `name` field must be formatted as: `library_name/sample_barcode/srt_barcode`.

*Example Row:*
```
chr1	12345	12359	100	+	library_A/AAGGCAGACG/GATTACA
```

#### Annotation File Format
-   A `.csv` or `.tsv` file with a **header**.
-   Columns must be: `library_name`, `sample_barcode`, `sample_name`, `group_name`.
-   **All `sample_name` values must be unique.**

*Example (tab-separated):*
```tsv
library_name	sample_barcode	sample_name	group_name
Library_A	AAGGCAGACG	Rep1_HyPBase	HyPBase
Library_A	TACCGCTGAC	Rep1_TF_A	TF_A
Library_B	AAGGCAGACG	Rep2_HyPBase	HyPBase
Library_B	TACCGCTGAC	Rep2_TF_A	TF_A
```

### Output Files of `multiplex_srt_seq_to_tf_binding.py`

-   `<output_dir>/final_results.tsv`: The main output table with all peak and sample statistics.
-   `<output_dir>/column_definitions.tsv`: Describes each column in `final_results.tsv`.
-   `<output_dir>/collapsed_per_sample/`: Directory with per-sample fragment data (`.parquet`).
-   `<output_dir>/sample_peaks/`: Directory with SRT barcode-refined peaks for each sample (`.parquet`).
-   `<output_dir>/pipeline.log`: A log file of all steps, warnings, and errors.

### Output Column Definitions

| Variable                                  | Description                                                                                             |
| ----------------------------------------- | ------------------------------------------------------------------------------------------------------- |
| `consensus_peak_id`                       | Unique identifier for the final merged peak (`chrom:start-end`).                                        |
| `chrom`                                   | Chromosome of the consensus peak.                                                                       |
| `consensus_peak_start`                    | Start coordinate of the consensus peak.                                                                 |
| `consensus_peak_end`                      | End coordinate of the consensus peak.                                                                   |
| `sample_name`                             | The sample identifier from the annotation file.                                                         |
| `group_name`                              | The experimental group identifier from the annotation file.                                             |
| `srt_bc_peak_start_for_sample`            | Start coordinate of the peak for this sample after SRT barcode refinement.                              |
| `srt_bc_peak_end_for_sample`              | End coordinate of the SRT barcode-refined peak for this sample.                                         |
| `sample_total_reads_in_consensus_peak`    | Total reads from this sample within the consensus peak.                                                 |
| `sample_total_fragments_in_consensus_peak`| Number of unique molecular fragments from this sample within the consensus peak.                        |
| `sample_total_insertions_in_consensus_peak`| Number of unique SRT barcodes from this sample within the consensus peak. **(Primary signal)** |
| `group_total_insertions_in_consensus_peak` | Sum of `sample_total_insertions_in_consensus_peak` for all samples in this group within the consensus peak. |
| ...                                       | *(Other columns provide additional detail on peak widths, counts, etc.)* |

---

## Script 2: `generate_insertion_maps.py`

This post-processing script generates `bedGraph` and `BigWig` files for visualizing insertion data in a genome browser. It identifies the single most likely insertion site for each unique SRT barcode and can normalize signal between experimental groups.

### Workflow

1.  **Per-Sample Insertion Mapping**: For each sample, processes the collapsed fragment data to determine the single most probable insertion coordinate for each unique SRT barcode using strand-aware logic.
2.  **Group Aggregation**: Aggregates the per-sample insertion maps based on `group_name`.
3.  **Normalization**: If multiple groups exist, calculates size factors based on the geometric mean of total unique insertions to normalize signal across groups, enabling direct comparison.
4.  **Validation**: Compares its calculated total insertion counts against the results from Script 1 as an integrity check.
5.  **Binned Sum Track Generation (Optional)**: If `--sum_window_size` is provided, generates additional BigWig tracks where the normalized signal is summed across fixed-width genomic windows.

### Usage of `generate_insertion_maps.py`

#### Dependencies
Requires all dependencies from Script 1, plus:

| Package                                                                  | Version | Purpose                             |
| ------------------------------------------------------------------------ | ------- | ----------------------------------- |
| [numpy](https://numpy.org/install/)                                      | 1.26.4+ | Numerical operations for binning    |
| [pyBigWig](https://pypi.org/project/pyBigWig/)                           | 0.3.22+ | Reading/writing BigWig files        |
| [bedGraphToBigWig](http://hgdownload.soe.ucsc.edu/admin/exe/)             | (any)   | (UCSC tool, required externally)    |

#### Installation

1.  Install additional Python packages:
    ```bash
    pip install numpy pyBigWig
    ```
2.  Install `bedGraphToBigWig` from the [UCSC binary utilities directory](http://hgdownload.soe.ucsc.edu/admin/exe/) or with conda:
    ```bash
    conda install -c bioconda ucsc-bedgraphtobigwig
    ```

#### Command Line Example

```bash
python generate_insertion_maps.py \
    --output_dir_of_multiplex_srt_seq_to_tf_binding /path/to/main/pipeline/output \
    --ins_map_output_dir /path/to/insertion_map_results \
    --annotation_file /path/to/annotation_file.tsv \
    --chrom_sizes /path/to/hg38.chrom.sizes \
    --srt_bc_dist_threshold 1 \
    --workers 8 \
    --sum_window_size 50
```

#### Arguments
* `--output_dir_of_multiplex_srt_seq_to_tf_binding`: **(Required)** Path to the output directory from Script 1.
* `--ins_map_output_dir`: **(Required)** Directory to save the output insertion maps.
* `--annotation_file`: **(Required)** The same annotation file used for Script 1.
* `--chrom_sizes`: **(Required)** Path to a two-column, tab-separated chromosome sizes file.
* `--srt_bc_dist_threshold`: **(Required)** Must be the same value used for Script 1.
* `--sum_window_size`: (Optional) If set, generates binned summary tracks (e.g., `50`).

### Output Files of `generate_insertion_maps.py`

All files are saved within the `--ins_map_output_dir`.
-   `raw_unique_insertions_per_sample_bedgraph/` & `..._bigwig/`
-   `raw_unique_insertion_count_per_group_bedgraph/` & `..._bigwig/`
-   `size_normalized_unique_insertions_per_group_bedgraph/` & `..._bigwig/`
-   `binned_sum_bigwig_{window_size}bp/` (Optional)

---

## Key Data Processing Details

*(Applies to both scripts)*

### 'strand', 'start', 'end', and Read Orientation

The coordinates in the input `qbed` file relate to the R2 read and the strand of the transposon insertion.

-   **For a `+` strand row:**
    -   The transposon inserted on the **`+`** strand of the reference.
    -   The R2 read is oriented **5' ← 3'** relative to the reference.
    -   The **`end` coordinate** is the true 3' end of the read.
    -   The **smallest `end` coordinate** for a given SRT barcode is the position closest to the insertion site.

-   **For a `-` strand row:**
    -   The transposon inserted on the **`-`** strand of the reference.
    -   The R2 read is oriented **5' → 3'** relative to the reference.
    -   The **`start` coordinate** is the true 3' end of the read.
    -   The **largest `start` coordinate** for a given SRT barcode is the position closest to the insertion site.

### Key Logic for Converting Fragment-based Peaks to SRT Barcode-based Peaks

This logic is central to `multiplex_srt_seq_to_tf_binding.py`.

1.  Initial **fragment-based peaks** are broad regions called from all unique fragment coordinates.
2.  To create refined **SRT barcode-based peaks**, all fragments within a fragment-based peak are grouped by their SRT barcode.
3.  For each SRT barcode, the single coordinate most proximal to the transposon is found using the strand-aware logic described above.
4.  The **SRT barcode-based peak start** is the *minimum* of all these proximal coordinates.
5.  The **SRT barcode-based peak end** is the *maximum* of all these proximal coordinates.
6.  A 100 bp extension is added to the final boundaries to provide a margin of error.

This process narrows the broad peak to a more precise region defined by the unique insertion events.



# Processing Multiplex SRT Sequencing Data to TF Binding Sites (Peaks) and Visualization Files (bedGraph, bigWig)
1. ***multiplex_srt_seq_to_tf_binding.py*** - Convert aligned sequencing data of multiplexed self-reporting transposon sequencing data into putative transcription factor binding sites.
   - This output is used to generate a count matrix DESeq2 to define TF binding sites that are enriched above the unfused transposase (HyPBase) control and differential peak among transcription factors.
   - The peaks enriched above above the unfused transposase (HyPBase) control are considered the true binding sites.

2. ***generate_insertion_maps.py*** - Generate bedGraph and bigwig files of unique insertions.
   - This output is used to for Visualization of unique insertions on gene tracks.

---

## Overview of method
- This method leverages the self-reporting transposon technology.
   - It is "self-reporting" because it contains a promoter that generates RNA containing the junction of the terminal repeat and downstream genomic DNA.
   - This RNA can be sequenced to define the approximate transcritpoin factor binding position.
     
- Because the sample barcode and SRT barcode are encoded into the genomic DNA in the transposon, this method is capable of antibody-independent, semi-multiplexed TF binding site experiments. 
   - Samples are transfected/electroporated with a large pool of plasmids. 
   - Each plasmid pool consists of a large, diverse set of SRT barcodes all with the same sample barcode.
   - All samples transfected with a unique sample barcode pool are combined after intial transfection/electroporation.
     
- After pooling all samples, the workflow is RNA extraction -> cDNA generation -> amplification of the self-reporting transpsons -> library preparation from amplicons with the Illumina Nextera XT kit -> sequencing.
   - Following sequencing, the reads originating from each sample are demultiplexed on the basis of the sample barcode.
   - The unique transposon insertions and their position, representing location and signal of transcription factor, of each sample are defined by the SRT barcode and genomic coordinate of the transposon.
 
 ![Picture1](https://github.com/user-attachments/assets/46acac4e-7843-49a7-9e7e-2a5b616febc6)

---

## Description of multiplex_srt_seq_to_tf_binding.py script
- This pipeline processes multiplexed self-reporting transposon data to identify transcription factor (TF) binding sites.
- It supports high-throughput experiments with many barcoded samples, corrects technical errors in barcode assignment, collapses redundant reads, calls peaks per sample, and generates a consensus set of reproducible peaks across all experimental conditions.

 ![Picture1](https://github.com/user-attachments/assets/dd45c0cc-8132-4d65-af3a-1b24befb073c)

## Pipeline steps 
### Phase 1: Loading, Correcting, and Annotating All Reads
1. **Input loading and barcode parsing:**  
   - Reads all qbed/bed files and parses the `name` field into `library_name`, `sample_barcode_raw`, and `srt_bc`.
   - To avoid memory issues with a lot of data, each input file is looped through one at a time. For each file, barcode correction and annotation is performed. Then the data is partitioned by sample_name to temporary .parquet files in a partitioned_temp/ directory.
   - After all input files have been processed, the intermediate files are grouped sample_name and concatenated to create the final collapsed_per_sample/{sample_name}.parquet files.
   - Temporary directory with the intermediate files is removed after completion.
2. **Barcode correction:**  
   - Corrects sample barcodes against the annotation whitelist, using Hamming distance up to the `--sample_barcode_dist_threshold`.
3. **Annotation:**  
   - Joins corrected reads with the annotation file, assigning `sample_name` and `group_name`.
4. **Per-sample partitioning:**  
   - Saves all rows of each `sample_name` as a `.parquet` file for efficient downstream processing.
   - The .parquet of each sample's fragments/qbed partition is saved in the `<output_dir>/collapsed_per_sample directory` within the specified output directory.
     
### Phase 2: Calling Peaks on Each Sample's Data
5. **Fragment-based peak calling:**  
   - For each sample, loads each partitioned `.parquet` file and merges rows identical in all fields except `sample_barcode_raw` and sums their reads.
   - Then peaks are called with pycallingcards on these fragments, which may or may not be associated with more than one SRT barcode.
   - The goal here is to define regions in the genome of at least on transpsoson insertion that is supported by at least 5 differentially fragmented molecules to remove noise.
   - Per sample fragment-based peaks are saved in `<output_dir>/fragment_peaks` within the specified output directory.
6. **Fragments-to-fragment-based peak mapping (intersection)**
   - The deduplicated fragments of each sample are intersected to their own peak set to map the fragments to the peaks.
   - This is required for step 7.
7. **SRT barcode-based peak refinement:**
   - The fragment-based peaks refines peak boundaries using strand-aware logic of the position per SRT barcode that is the most proximal to the junction of the transposon and genomic DNA.
   - This step will also count the number of unique transposon insertions (equal to unique SRT barcodes) in the fragment-based peak. The unique transposon insertions acts as the signal of TF binding.
   - Per sampel SRT barcode-based peaks are saved in `<output_dir>/sample_peaks directory` within the specified output directory as a .parquet file.

### Phase 3: Generating Consensus Peaks and Final Output
8. **Consensus peak generation:**  
   - Merges sample-specific peaks to produce consensus peaks across all samples and groups (i.e., everything specified in the annotation file).
9. **Sample peak set-to-consensus peak mapping (intersection):**  
   - Intersects all sample peak sets with the consensus peaks to map which sample peaks were merged in the consensus peak.
10. **Output aggregation:**  
   - Each row is a unique consensus peak and sample_name pair populated with the per-sample peak statistics (e.g., fragment-based peak coordinates, SRT barcode-based peak coordinates, total reads, total fragments, and unique insertions within each of those peaks). The attributes of all samples in the same group are summed per consensus peak to produces the per-group statistics (`final_results.tsv`).
   - A `column_definitions.tsv` file describing all output columns is also saved in the output directory.


## Usage of multiplex_srt_seq_to_tf_binding.py
### Dependencies

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

*Other versions may be compatible*

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

**Checking installed versions**
To check package versions, run the following in the terminal:
```bash
python -c "import polars, pandas, pybedtools, tqdm, umi_tools, pycallingcards as cc; print('polars:', polars.__version__); print('pandas:', pandas.__version__); print('pybedtools:', pybedtools.__version__); print('tqdm:', tqdm.__version__); print('umi_tools:', umi_tools.__version__); print('pycallingcards:', cc.__version__)"
bedtools --version
```

---

### Default parameters command line
```bash
python multiplex_srt_seq_to_tf_binding.py \
  --input_dir /path/to/qbed_files \
  --output_dir /path/to/results \
  --annotation_file /path/to/annotation_file.tsv \
  --workers 10 \
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
min_insertions, minlen, extend, maxbetween, minnum, lam_win_size all refer to the pycallingcards peak caller settings. Full description available [here](https://pycallingcards.readthedocs.io/en/latest/api/reference/pycallingcards.preprocessing.call_peaks.html).


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
  - Maximum Hamming distance allowed for correcting sample barcodes (default: 2).
  - The value set here must be less than the minimal Hamming distance across all sample barcodes.
  - A value of 2 means that sample barcodes with ≤ 2 mismatches from each sample barcode in the annotation file will be corrected to the sample barcode in the annotation file.
  - A higher value here means a less stringent filter to correct sample barcode sequencing errors into the sample barcode in the annotation file.
  
- `--srt_bc_dist_threshold`  
  - Maximum Hamming distance for SRT barcode clustering and correction (default: 1).
  - A value of 1 means that only SRT barcodes with ≤ 1 mismatches from each other will be considered the same SRT barcode.
  - A higher value means more aggressive SRT barcode clustering and correction.

- Additional parameters for fragment-based peak calling (see script source for defaults and descriptions):
  - `--pvalue_cutoff`, `--pvalue_adj_cutoff`, `--min_insertions`, `--minlen`, `--extend`, `--maxbetween`, `--minnum`, `--lam_win_size`
  - It is recommended to leave these parameters at the default values.
  - The purpose of these parameters is to define all regions of concentrated fragments for subsequent SRT-barcode-based peak boundary definitions and donwstream analysis.
  - Fragment peaks with <5 fragments following SRT barcode correction are considered noise and discarded from the analysis.

### Other considerations
- All pycallingcards peak calling parameters can be overridden via the command line.
- The pipeline can be parallelized across as many CPUs as available using `--workers`.
- Sample barcode and SRT barcode correction stringency can be adjusted by `--sample_barcode_dist_threshold` and `--srt_bc_dist_threshold`.



### Supported input file formats 
- The file should ***not*** have a header.
- `.qbed`, `.bed`, `.qbed.gz`, `.bed.gz`
- Each file must contain columns in the order:  
  `chrom`, `start`, `end`, `reads`, `strand`, `name`
- The `name` field should be formatted as:  
  `library_name/sample_barcode/srt_barcode`

**Example input file row (no header)**
```
chr1   12345	  12359	 100	    +	        [library_name]/[sample_barcode]/[srt_barcode]
```

### Annotation file format

This file should be a `.csv` or `.tsv` (tab-separated values) format. The columns must appear **in this order**:  
**library_name**, **sample_barcode**, **sample_name**, **group_name**.

**Example annotation file (tab-separated):**
```
library_name	sample_barcode	sample_name	group_name
Library_A	AAGGCAGACG	Rep1_HyPBase	HyPBase
Library_A	TACCGCTGAC	Rep1_TF_A	TF_A
Library_A	AAGATTAGAC	Rep1_TF_B	TF_B
Library_A	AGTATGACCG	Rep1_TF_C	TF_C
Library_A	AGTATGACCG	Rep2_TF_A	TF_A
Library_B	AAGGCAGACG	Rep2_HyPBase	HyPBase
Library_B	TACCGCTGAC	Rep2_TF_A	TF_A
Library_B	AAGATTAGAC	Rep2_TF_B	TF_B
Library_B	AGTATGACCG	Rep2_TF_C	TF_C
```
- **All `sample_name` values must be unique.**
- **The same `sample_name` may be used across different libraries.**
   - In the example above, the library name will be Library_A__Library_B for Rep2_TF_A.
   - If you want the separate libraries to be treated as independent, you must use unique labels, such as Rep2_TF_A_LibA and Rep2_TF_A_LibB.

---

## Output files of multiplex_srt_seq_to_tf_binding.py

- `<output_dir>/final_results.tsv`  
  Tab-separated table containing all peak and sample statistics. See variable definitions below.
- `<output_dir>/column_definitions.tsv`  
  Tab-separated table describing all columns in `final_results.tsv`.
- `<output_dir>/collapsed_per_sample/sample.parquet`
  The subset of the qbed rows corresponding to the sample_name in the annotation file.
- `<output_dir>/sample_peaks/sample_peaks.parquet`
  The peaks called after SRT barcode deduplication and defining the SRT barcode peak boundaries for each sample_name in the annotation file.
- `<output_dir>/pipeline.log`  
  Log file of all steps, warnings, and errors.
- `<output_dir>/pybedtools_tmp/`  
  Temporary files for bedtools operations (auto-cleaned at completion).
  
---
## Key details of data processing (applies to both multiplex_srt_seq_to_tf_binding.py and generate_insertion_maps.py)

### 'strand', 'start', 'end', and read orientation in relation to the tranpsoson
- The R2 read end is reported in the input file.
- The strand in the qbed refers the strand that the **transposon** was inserted into, not the strand of the read.
- By convention of the alignment pipeline,
   - **'end' coordinate** for a **'+'** strand row is the true end of the read
   - **'start' coordinate** for a **'-'** strand row is the true end of the read.
     
- Reads that hit the transposon are truncated to the last base prior to beginning of the transposon.

### Summary:
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
      
### Key logic for converting the fragment-based peaks to SRT barcode-based peaks
- **Fragment-based peaks:** Peaks called from the initial input qbed file where every row represents one uniquely fragmented molecule for which the read ended at a different position.
- **SRT barcode-based peaks:** Refined fragment-based peaks where the start and end coordinates are defined based on the unique SRT barcode positions within the peak.
  
  - Peak start = the **smallest** value between the minimum 'end' coordinate of '+' strand fragments and maximum 'start' coordinate of '-' strand fragments.
    - This represents the most **upstream** value that is the most proximal to the transposon between the '+' and '-' strand fragments.
  - Peak end = the **largest** value between the minimum 'end' coordinate of '+' strand fragments and maximum 'start' coordinate of '-' strand fragments.
    - This represents the most **downstream** value that is the most proximal to the transposon between the '+' and '-' strand fragments.
  - There is a 100bp extension applied to both ends of the SRT barcode-based peak.
    - If the SRT barcode-based peak has only fragments of one SRT barcode, then the width would be 0. With 100bp extension to either side, this becomes 200.

       
- Note that in both instances, the minimum 'end' coordinate of '+' strand fragments and maximum 'start' coordinate of '-' strand fragments are used.
  - This is because these value represent the most proximal position of all fragments of a given SRT barcode to the transposon.
  - It is the best approximation of the true transposon insertion site if the read does not make contact with the transposon.

- The SRT barcode-based peak width will nearly always be smaller than the fragment-based peak width.
  - However on rare occasion, the fragment-based peak width will be larger. The following logic explains how this can happen:
    - The peak caller will always apply an extension to both side of the fragment-based peak (set by extend parameter, default 200bp to both sides)
      - A 100bp extension is applied to both sides of the SRT barcode-based peak.
      - In some cases, a fragment will occur in the extension window of the fragment-based peak *and* will be <100bp from the peak start or end.
      - If that fragment is from the strand-aware most proximal position of the SRT barcode to the transposon, then the 100bp extension applied after defining the SRT barcode peak will generate a SRT barcode-based peak width > fragment-based peak width.


---


## Output column names with descriptions of multiplex_srt_seq_to_tf_binding.py

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
