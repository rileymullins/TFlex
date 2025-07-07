# Processing Multiplex SRT Sequencing Data to TF Binding Sites and Visualization 
1. ***multiplex_srt_seq_to_tf_binding.py***
   - Convert aligned sequencing data of multiplexed self-reporting transposon sequencing data into putative transcription factor binding sites.
   - This output is used to define TF binding sites that are enriched above the unfused transposase (HyPBase) control and differential peak among transcription factors.
   - The peaks enriched above above the unfused transposase (HyPBase) control are considered the true binding sites.

2. ***generate_insertion_maps.py***
   - Generate bedGraph and bigwig files of unique transposon insertions.
   - This output is used to for visualization of unique insertions on gene tracks.

---

## Table of Contents

- [Overview of the Method](#overview-of-the-method)
- [Script 1: `multiplex_srt_seq_to_tf_binding.py`](#script-1-multiplex_srt_seq_to_tf_bindingpy)
    - [Pipeline steps](#pipeline-steps)
    - [Usage of multiplex_srt_seq_to_tf_binding.py](#usage-of-multiplex_srt_seq_to_tf_bindingpy)
    - [Output files of multiplex_srt_seq_to_tf_binding.py](#output-files-of-multiplex_srt_seq_to_tf_bindingpy)
    - [Output column names with descriptions of multiplex_srt_seq_to_tf_binding.py](#output-column-names-with-descriptions-of-multiplex_srt_seq_to_tf_bindingpy)
- [Script 2: `generate_insertion_maps.py`](#script-2-generate_insertion_mapspy)
    - [Workflow](#workflow)
    - [Usage of `generate_insertion_maps.py`](#usage-of-generate_insertion_mapspy)
    - [Output Files of `generate_insertion_maps.py`](#output-files-of-generate_insertion_mapspy)
- [Key details of data processing](#key-details-of-data-processing-applies-to-both-multiplex_srt_seq_to_tf_bindingpy-and-generate_insertion_mapspy)
    - [Read orientation in relation to the tranpsoson](#strand-start-end-and-read-orientation-in-relation-to-the-tranpsoson)
    - [Logic for converting fragment-based peaks to SRT barcode-based peaks](#key-logic-for-converting-the-fragment-based-peaks-to-srt-barcode-based-peaks)
    - [Logic for visualization](#key-logic-for-visualization)

---

## Overview of the Method

- This method leverages self-reporting transposon (SRT) technology.
- An SRT is "self-reporting" because it contains a promoter that drives the transcription of RNA containing the junction of the transposon's terminal repeat and the downstream genomic DNA.
- Sequencing this RNA maps the transposon insertion site, which represents the approximate TF binding site.

The method is antibody-independent and semi-multiplexed. Because both the sample barcode and the SRT barcode are inserted into genomic DNA, cells from different conditions can be combined following transfection/electroporation if they have distinct sample barcodes.

![Picture1](https://github.com/user-attachments/assets/46acac4e-7843-49a7-9e7e-2a5b616febc6)

**Workflow Summary:**

1.  **Transfection:** Samples are transfected with a large, diverse pool of plasmids. Each pool contains a unique sample barcode but many SRT barcodes.
2.  **Pooling:** All samples, each with its unique sample barcode, are pooled after transfection.
3.  **Library Prep & Sequencing:** The workflow is RNA extraction, cDNA synthesis, SRT amplification, library preparation, and 2x150bp sequencing.
4.  **Demultiplexing:** Post-sequencing, reads are assigned to their original samples based on the sample barcode.
5.  **Mapping:** The unique transposon insertion sites are identified by their genomic coordinates and the SRT barcode, representing the location and signal of a TF binding event, respectively.

---

## Script 1: `multiplex_srt_seq_to_tf_binding.py`
- This pipeline processes multiplexed self-reporting transposon data to identify transcription factor (TF) binding sites.
- It supports high-throughput experiments with many barcoded samples, corrects technical errors in barcode assignment, collapses redundant reads, calls peaks per sample, and generates a consensus set of reproducible peaks across all experimental conditions.

![Picture1](https://github.com/user-attachments/assets/dd45c0cc-8132-4d65-af3a-1b24befb073c)

### Pipeline steps
#### Phase 1: Loading, Correcting, and Annotating All Reads
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

#### Phase 2: Calling Peaks on Each Sample's Data
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
   - Per sample SRT barcode-based peaks are saved in `<output_dir>/sample_peaks directory` within the specified output directory as a .parquet file.

#### Phase 3: Generating Consensus Peaks and Final Output
8. **Consensus peak generation:**
   - Merges sample-specific peaks to produce SRT barcode consensus peaks across all samples and groups (i.e., everything specified in the annotation file).
9. **Sample peak set-to-consensus peak mapping (intersection):**
   - Intersects all sample peak sets with the SRT barcode consensus peaks to map which sample peaks were merged in the SRT barcodeconsensus peak.
10. **Output aggregation:**
    - Each row is a unique consensus peak and sample_name pair populated with the per-sample peak statistics (e.g., fragment-based peak coordinates, SRT barcode-based peak coordinates, total reads, total fragments, and unique insertions within each of those peaks). The attributes of all samples in the same group are summed per consensus peak to produces the per-group statistics (`final_results.tsv`).
    - A `column_definitions.tsv` file describing all output columns is also saved in the output directory.


### Usage of multiplex_srt_seq_to_tf_binding.py
#### Dependencies

The following packages and versions were used to create and test this script using **Python 3.12.10**:

| Package | Version | Purpose |
|---|---|---|
| [polars](https://docs.pola.rs/user-guide/installation/) | 1.31.0 | Fast DataFrames manipulation |
| [pandas](https://pandas.pydata.org/docs/getting_started/install.html) | 2.2.3 | DataFrames, pybedtools input |
| [pybedtools](https://daler.github.io/pybedtools/main.html) | 0.12.0 | Genomic intervals intersections |
| [pycallingcards](https://pycallingcards.readthedocs.io/en/latest/#) | 1.0.0 | Fragment-based peak calling |
| [UMI-tools](https://umi-tools.readthedocs.io/en/latest/INSTALL.html) | 1.1.6 | SRT barcode clustering to deduplicate within peak |
| [tqdm](https://pypi.org/project/tqdm/) | 4.67.1 | Progress bars |
| [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html) | 2.31.1 | (external, required for pybedtools) |

*Other versions may be compatible*

#### Installation

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
#### Default parameters command line
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


#### Arguments
- `--input_dir`
  Directory containing input `.qbed`, `.bed`, `.qbed.gz`, or `.bed.gz` files. Each file should be in standard qbed/bed format.
- `--output_dir`
  Directory where output files will be written.
- `--annotation_file`
  Path to the annotation file (see format below).

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

- Additional parameters for fragment-based peak calling (see pycallingcards website for more details):
  - `--pvalue_cutoff`, `--pvalue_adj_cutoff`, `--min_insertions`, `--minlen`, `--extend`, `--maxbetween`, `--minnum`, `--lam_win_size`
  - It is recommended to leave these parameters at the default values.
  - The purpose of these parameters is to define all regions of concentrated fragments for subsequent SRT-barcode-based peak boundary definitions and donwstream analysis.
  - Fragment peaks with <5 fragments following SRT barcode correction are considered noise and discarded from the analysis.


#### Supported input file formats
- The file should ***not*** have a header.
- `.qbed`, `.bed`, `.qbed.gz`, `.bed.gz`
- Each file must contain columns in the order:
  `chrom`, `start`, `end`, `reads`, `strand`, `name`
- The `name` field should be formatted as:
  `library_name/sample_barcode/srt_barcode`

**Example input file row (no header)**
```
chr1    12345    12359   100     +         [library_name]/[sample_barcode]/[srt_barcode]
```

#### Annotation file format

This file should be a `.csv` or `.tsv` (tab-separated values) format. The columns must appear **in this order**:
**library_name**, **sample_barcode**, **sample_name**, **group_name**.

**Example annotation file (tab-separated):**
```
library_name	sample_barcode	sample_name	group_name
Library_A	AAGGCAGACG	Rep1_HyPBase	HyPBase
Library_A	TACCGCTGAC	Rep1_TF_A	   TF_A
Library_A	AAGATTAGAC	Rep1_TF_B	   TF_B
Library_A	AGTATGACCG	Rep1_TF_C	   TF_C
Library_A	AGTATGACCG	Rep2_TF_A	   TF_A
Library_B	AAGGCAGACG	Rep2_HyPBase	HyPBase
Library_B	TACCGCTGAC	Rep2_TF_A	   TF_A
Library_B	AAGATTAGAC	Rep2_TF_B	   TF_B
Library_B	AGTATGACCG	Rep2_TF_C	   TF_C
```
- **The same `sample_name` may be used across different libraries.**
  - In the example above, the library name will be Library_A__Library_B for Rep2_TF_A.
  - If you want the separate libraries to be treated as independent, you must use unique labels, such as Rep2_TF_A_LibA and Rep2_TF_A_LibB.

---

### Output files of multiplex_srt_seq_to_tf_binding.py

- `<output_dir>/merged_srt_barcode_based_peaks.tsv`
  Tab-separated table containing all merged (bedtools merge) peaks with per group and per sample and sample statistics. See variable definitions below.
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

### Output column names with descriptions of multiplex_srt_seq_to_tf_binding.py

| Variable | Description |
|---|---|
| consensus_peak_id | chrom:consensus_peak_start-consensus_peak_end |
| chrom | chromosome |
| consensus_peak_start | start coordinate of final merged, pan-dataset consensus peak |
| consensus_peak_end | end coordinate of final merged, pan-dataset consensus peak |
| consensus_peak_width | width (end – start) of final merged, pan-dataset consensus peak |
| num_samples_in_consensus_peak | number of unique sample_name values in this consensus peak |
| num_groups_in_consensus_peak | number of unique group_name values in this consensus peak |
| sample_name | sample identifier (from annotation file). This is the peak calling unit; could be a replicate (e.g., replicate-1_TF-A), a unique time point (e.g., timepoint1_replicate1_TF-A), or a unique experimental condition (e.g., drugX_replicate1_TF-A). |
| group_name | group identifier (from annotation file). Used to aggregate stats across samples belonging to the same broader group. Examples: "TF_A", "timepoint1_TF-A", or "drugX_TF-A". |
| fragment_peak_start_for_sample | fragment-based peak start coordinate for this sample_name |
| fragment_peak_end_for_sample | fragment-based peak end coordinate for this sample_name |
| fragment_peak_width_for_sample | fragment-based peak width (end – start) for this sample_name |
| srt_bc_peak_start_for_sample | SRT-barcode-based peak start coordinate for this sample_name |
| srt_bc_peak_end_for_sample | SRT-barcode-based peak end coordinate for this sample_name |
| srt_bc_peak_width_for_sample | SRT-barcode-based peak width (end – start) for this sample_name |
| sample_total_reads_in_consensus_peak | total read count in consensus peak for this sample_name |
| sample_total_fragments_in_consensus_peak | total fragment count (merged qbed rows after SRT-barcode correction and deduplication) in consensus peak for this sample_name |
| sample_total_insertions_in_consensus_peak | total unique insertions (unique SRT barcodes) in consensus peak for this sample_name |
| group_total_reads_in_consensus_peak | sum of sample_total_reads_in_consensus_peak across all sample_name values in this group_name |
| group_total_fragments_in_consensus_peak | sum of sample_total_fragments_in_consensus_peak across all sample_name values in this group_name |
| group_total_insertions_in_consensus_peak | sum of sample_total_insertions_in_consensus_peak across all sample_name values in this group_name |

---

## Script 2: `generate_insertion_maps.py`

This post-processing script generates `bedGraph` and `BigWig` files for visualizing insertion data in a genome browser. It identifies the single most likely insertion site for each unique SRT barcode and can normalize signal between experimental groups.

### Workflow

1.  **Per-Sample Insertion Mapping**: For each sample, processes the collapsed fragment data to determine the single most probable insertion coordinate for each unique SRT barcode using strand-aware logic.
2.  **Group Aggregation**: Aggregates the per-sample insertion maps based on `group_name`.
3.  **Normalization**: If multiple groups exist, calculates size factors based on the geometric mean of total unique insertions to normalize signal across groups, enabling direct comparison.
4.  **Validation**: Compares its calculated total insertion counts against the results from Script 1 to confirm total unique counts for all conditions is the same.
5.  **Binned Sum Track Generation**: If `--sum_window_size` is provided (default 50bp), generates additional BigWig tracks where the normalized signal is summed across fixed-width genomic windows.
    - Once the script runs fully the first time, the output files will be detected, and the script will skip to the binning step. This allows for testing of multiple bin sizes.
    - This is recommended because the insertions site positions used are the best approximation of where the true transposon/DNA junction lies.
    - Summing signal in bins only affects visualization and is simply a way to integrate the unique insertions that likely end in the same position if we were to know exactly where the read hits the junction.

### Usage of `generate_insertion_maps.py`

#### Dependencies
Requires all dependencies from Script 1, plus:

| Package | Version | Purpose |
|---|---|---|
| [numpy](https://numpy.org/install/) | 1.26.4+ | Numerical operations for binning |
| [pyBigWig](https://pypi.org/project/pyBigWig/) | 0.3.22+ | Reading/writing BigWig files |
| [bedGraphToBigWig](http://hgdownload.soe.ucsc.edu/admin/exe/) | (any) | (UCSC tool for bedGraph to bigWig conversion) |

- Also, chrom.sizes files are needed bigWig file generation. Those for hg38 are provided in this repo.


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
- `raw_unique_insertions_per_sample_bedgraph/` & `..._bigwig/`
- `raw_unique_insertion_count_per_group_bedgraph/` & `..._bigwig/`
- `size_normalized_unique_insertions_per_group_bedgraph/` & `..._bigwig/`
- `binned_sum_bigwig_{window_size}bp/` (Optional)

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

  - **Peak start** = the **smallest** value between the minimum 'end' coordinate of '+' strand fragments and maximum 'start' coordinate of '-' strand fragments.
    - This represents the most **upstream** value that is the most proximal to the transposon between the '+' and '-' strand fragments.
  - **Peak end** = the **largest** value between the minimum 'end' coordinate of '+' strand fragments and maximum 'start' coordinate of '-' strand fragments.
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

### Key logic for visualization
- **Frequently, multiple fragments of the same SRT barcode will be present in the same peak.**
  - When generating SRT barcode-based peaks, the count of unique SRT barcodes suffices.
  - However, for visualization a single representative genomic coordinate must be chosen for that SRT barcode.
  - Using the same strand-aware logic explained directly above,
    - The most representative coordinate for a **'+' strand SRT barcode** is the *minimum* (most upstream) 'end' coordinate among all fragments of that SRT barcode in the SRT barcode-based peak.
    - The most representative coordinate for a **'-' strand SRT barcode** is the *maximum* (most downstream) 'start' coordinate among all fragments of that SRT barcode in the SRT barcode-based peak.
