# **Annotation-Driven Calling Card Analysis Pipeline**

This repository contains a suite of scripts to perform a complete annotation-driven calling card data processing and analysis workflow. The pipeline starts with raw sequencing data in qBED format, identifies precise transposon insertion sites, calls peaks of enrichment, performs differential analysis, and finishes with functional annotation and motif discovery.

The workflow is designed to be modular, with four main scripts that pass their outputs to the next stage.

## **Outline**

* [Workflow Diagram](https://www.google.com/search?q=%23workflow-diagram)  
* [Installation Requirements](https://www.google.com/search?q=%23installation-requirements)  
  * [Conda Environment](https://www.google.com/search?q=%231-conda-environment-recommended)  
  * [Python Dependencies](https://www.google.com/search?q=%232-python-dependencies)  
  * [R Dependencies](https://www.google.com/search?q=%233-r-dependencies)  
  * [External Software](https://www.google.com/search?q=%234-external-software--utilities)  
* [Scripts Overview](https://www.google.com/search?q=%23scripts-overview)  
  * [step\_1\_raw\_qbed\_to\_insertion\_maps.py](https://www.google.com/search?q=%23step_1_raw_qbed_to_insertion_mapspy)  
  * [step\_2\_run\_span\_peak\_caller.sh](https://www.google.com/search?q=%23step_2_run_span_peak_callersh)  
  * [step\_3\_generate\_peak\_by\_group\_count\_matrices.py](https://www.google.com/search?q=%23step_3_generate_peak_by_group_count_matricespy)  
  * [step\_4\_DESeq2\_Diff\_Peaks\_HOMER-Annotations\_and\_Motifs.R](https://www.google.com/search?q=%23step_4_deseq2_diff_peaks_homer-annotations_and_motifsr)

## **Workflow Diagram**

graph TD  
    A\[Raw qBED Files\] \--\> B(step\_1\_raw\_qbed\_to\_insertion\_maps.py);  
    B \--\> C{Group-level\<br\>Insertion Maps\<br\>(.bedgraph/.qbed)};  
    C \--\> D(step\_2\_run\_span\_peak\_caller.sh);  
    D \--\> E{SPAN Peaks\<br\>(.peak)};  
    subgraph " "  
        C \--\> F;  
        E \--\> F;  
    end  
    F(step\_3\_generate\_peak\_by\_group\_count\_matrices.py) \--\> G{Per-Group\<br\>Count Matrices\<br\>(.tsv)};  
    G \--\> H(step\_4\_DESeq2\_Diff\_Peaks\_HOMER-Annotations\_and\_Motifs.R);  
    H \--\> I{Differential Peaks\<br\>Annotations\<br\>Motifs & Plots};

## **Installation Requirements**

This pipeline relies on Python, R, and several external command-line tools. The following instructions outline how to set up the necessary environment. We recommend using a conda environment to manage Python and R dependencies.

#### **1\. Conda Environment (Recommended)**

First, install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) if you don't have it. Then, create and activate a new environment.

\# Create a new conda environment with Python and R  
conda create \-n callingcard\_env python=3.9 r-base=4.2 \-y

\# Activate the environment  
conda activate callingcard\_env

#### **2\. Python Dependencies**

With the conda environment active, install the required Python packages using pip.

pip install polars pandas numpy pybedtools tqdm umi-tools pycallingcards pybigwig

#### **3\. R Dependencies**

Install the required R packages from CRAN and Bioconductor. Run the R interpreter from your activated conda environment by simply typing R in your terminal, then run the following commands.

\# Install from CRAN  
install.packages(c("purrr", "stringr", "tibble", "tidyverse", "parallel", "RColorBrewer", "circlize", "viridis", "janitor", "BiocManager"))

\# Install from Bioconductor  
BiocManager::install(c("DESeq2", "ComplexHeatmap", "GenomicRanges"))

#### **4\. External Software & Utilities**

This pipeline requires four key pieces of external software that must be installed and accessible from your system's PATH.

* UCSC Kent Utilities: The bedGraphToBigWig utility is essential for Step 1\.  
  * Instructions: Download the pre-compiled binary for your system from the [UCSC downloads page](http://hgdownload.soe.ucsc.edu/admin/exe/). Place the binary in a directory that is included in your system's PATH (e.g., /usr/local/bin).  
* Java Runtime Environment (JRE): Required to run the SPAN peak caller in Step 2\.  
  * Instructions: Java can be installed via your system's package manager or directly from [Oracle](https://www.java.com/en/download/). You can also install it via conda:  
    conda install \-c conda-forge openjdk

* SPAN Peak Caller: This is the Java-based peak caller used in Step 2\.  
  * Instructions: Download the span.jar file from the [SPAN GitHub repository's release page](https://www.google.com/search?q=https://github.com/AndersenLab/span/releases). You will provide the path to this .jar file when running the script.  
* HOMER (Hypergeometric Optimization of Motif EnRichment): Required for the annotation and motif analysis in Step 4\.  
  * Instructions: Follow the detailed installation guide on the [HOMER website](http://homer.ucsd.edu/homer/introduction/install.html). The installation typically involves running a configureHomer.pl script, which will also download the necessary genome data. Ensure the HOMER bin directory is added to your system's PATH.

## **Scripts Overview**

### step\_1\_raw\_qbed\_to\_insertion\_maps.py

This is the foundational script of the pipeline. It takes raw qBED/BED files, which contain fragment coordinates and associated barcode information, and processes them into high-resolution, 1-bp unique insertion site maps. It combines initial data cleaning, annotation, per-sample peak calling, and precise insertion site determination into a single, parallelized workflow.

#### **Key Features**

* Barcode Correction: Corrects sequencing errors in sample barcodes by comparing them against a user-provided whitelist (from the annotation file) using Hamming distance.  
* Data Annotation & Partitioning: Annotates reads with sample and group information and partitions the raw data into separate files for each sample.  
* Fragment-Based Peak Calling: For each individual sample, it performs a preliminary fragment-based peak call using pycallingcards. This step helps to focus the subsequent analysis on regions with at least some signal.  
* Precise 1-bp Insertion Site Logic: This is the core of the script. Within each preliminary peak, it performs UMI/SRT barcode clustering and deduplication. Based on the strand of the fragment, it identifies the most probable 1-bp insertion coordinate (the 3' end of the transposon insertion event).  
* Aggregation & Normalization: Aggregates the unique 1-bp insertion maps for all samples within an experimental group. It then calculates size factors and generates both raw and normalized signal tracks.  
* Parallel Processing: Uses concurrent.futures.ProcessPoolExecutor to process multiple samples simultaneously, significantly speeding up the workflow.

#### **Inputs**

* \--input\_dir: A directory containing raw .qbed or .bed files. The name column (column 6\) must be in the format: {library\_name}/{sample\_barcode\_raw}/{srt\_bc}.  
* \--output\_dir: The main directory where all results will be stored.  
* \--annotation\_file: A CSV or TSV file with four columns: library\_name, sample\_barcode, sample\_name, group\_name.  
* \--chrom\_sizes: A tab-separated file containing chromosome names and their lengths (e.g., hg38.chrom.sizes).  
* \--workers: Number of CPU cores to use for parallel processing.

#### **Outputs**

This script generates a structured output directory with the following subdirectories:

* raw\_unique\_insertions\_per\_sample\_bedgraph/: Raw 1-bp insertion counts for each individual sample.  
* raw\_unique\_insertions\_per\_sample\_bigwig/: BigWig versions of the per-sample counts.  
* raw\_unique\_insertion\_count\_per\_group\_bedgraph/: Raw 1-bp insertion counts aggregated by experimental group.  
* raw\_unique\_insertion\_count\_per\_group\_bigwig/: BigWig versions of the raw group counts.  
* raw\_unique\_insertion\_count\_per\_group\_qbed/: The group-level data converted back to qBED format, suitable for input into the SPAN peak caller (Step 2).  
* size\_normalized\_unique\_insertions\_per\_group\_bedgraph/: Group counts normalized by library size (size factors).  
* size\_normalized\_unique\_insertions\_per\_group\_bigwig/: BigWig versions of the normalized group counts.  
* binned\_normalized\_count\_sum\_bigwig\_window\_.../: Normalized signal summed across genomic bins of a specified size for easier visualization.  
* fragment\_peaks\_per\_sample/: Intermediate files containing the fragment-based peaks called for each sample.  
* pipeline.log: A log file detailing the script's execution.

#### **Usage**

python step\_1\_raw\_qbed\_to\_insertion\_maps.py \\  
    \--input\_dir /path/to/raw\_qbeds \\  
    \--output\_dir /path/to/step1\_output \\  
    \--annotation\_file /path/to/annotation.csv \\  
    \--chrom\_sizes /path/to/hg38.chrom.sizes \\  
    \--workers 16 \\  
    \--srt\_bc\_dist\_threshold 1 \\  
    \--sample\_barcode\_dist\_threshold 2 \\  
    \--pvalue\_cutoff 0.01

### step\_2\_run\_span\_peak\_caller.sh

This script is a Bash wrapper designed to run the SPAN (Signal-based Peak ANalyzer) Java program. It iterates through the group-level qBED files generated by Step 1 and calls peaks, which are regions of statistically significant insertion enrichment.

#### **Key Features**

* Automated Iteration: Automatically finds all .qbed files in the input directory.  
* Simple Configuration: Key SPAN parameters like bin size and FDR are easily configurable.  
* Wrapper: Simplifies the execution of the SPAN command-line tool.

#### **Inputs**

* \-j: Path to the span.jar file.  
* \-i: Input directory containing the .qbed files from Step 1 (raw\_unique\_insertion\_count\_per\_group\_qbed/).  
* \-o: Directory where the output peak files will be saved.  
* \-c: Path to the chromosome sizes file.  
* \-b: (Optional) Bin size for SPAN analysis. Default: 50.  
* \-f: (Optional) False Discovery Rate (FDR) cutoff. Default: 0.01.

#### **Outputs**

* For each input qBED file (e.g., GroupName.qbed), it produces a corresponding peak file (e.g., GroupName\_b50.span.peak). This is a 3-column BED file (chrom, start, end) defining the peak regions.

#### **Usage**

bash step\_2\_run\_span\_peak\_caller.sh \\  
    \-j /path/to/span.jar \\  
    \-i /path/to/step1\_output/raw\_unique\_insertion\_count\_per\_group\_qbed \\  
    \-o /path/to/step2\_output/span\_peaks \\  
    \-c /path/to/hg38.chrom.sizes \\  
    \-b 50 \\  
    \-f 0.01

### step\_3\_generate\_peak\_by\_group\_count\_matrices.py

This script integrates the results from the first two steps to produce the final count matrices required for differential analysis. It takes the SPAN-called peaks (Step 2\) and quantifies the number of raw insertions (from Step 1\) from every experimental group that fall within those peaks.

#### **Key Features**

* Peak Resizing: Ensures that all peaks have a minimum width (e.g., 200bp) by symmetrically extending smaller peaks from their midpoint. This creates a consistent window for counting insertions across all peaks.  
* Intersection-Based Counting: Uses the efficient pybedtools library to intersect each group's resized peak set against the raw insertion bedGraph files of *all* other groups.  
* Matrix Generation: For each group, it generates a comprehensive matrix where rows are the peaks called for that group, and columns contain the total number of insertions from every other group within those peaks.  
* Parallel Processing: Processes each group's peak file in parallel to speed up computation.

#### **Inputs**

* \--peak\_dir: Directory containing the SPAN peak files from Step 2\.  
* \--peak\_suffix: The unique suffix of the peak files to be processed (e.g., \_b50.span.peak).  
* \--bedgraph\_dir: Directory with per-group insertion count bedGraph files from Step 1 (raw\_unique\_insertion\_count\_per\_group\_bedgraph/).  
* \--output\_dir: Directory to save the output files.  
* \--annotation\_file: The same annotation file used in Step 1\.  
* \--workers: Number of parallel worker processes.

#### **Outputs**

* per\_group\_peak\_matrices/: A directory containing a separate count matrix for each group (e.g., GroupName\_peak\_matrix.tsv). Each file contains the peak coordinates and columns of insertion counts from all groups.  
* consensus\_peaks.log: A log file for the run.

#### **Usage**

python step\_3\_generate\_peak\_by\_group\_count\_matrices.py \\  
    \--peak\_dir /path/to/step2\_output/span\_peaks \\  
    \--peak\_suffix \_b50.span.peak \\  
    \--bedgraph\_dir /path/to/step1\_output/raw\_unique\_insertion\_count\_per\_group\_bedgraph \\  
    \--output\_dir /path/to/step3\_output \\  
    \--annotation\_file /path/to/annotation.csv \\  
    \--workers 16

### step\_4\_DESeq2\_Diff\_Peaks\_HOMER-Annotations\_and\_Motifs.R

This is the final, comprehensive R script for downstream analysis and interpretation. It performs differential peak analysis, identifies group-specific peaks, annotates them to genomic features, performs motif discovery, and generates a suite of publication-ready plots and tables.

#### **Workflow**

1. Configuration: User must define input paths and experimental group names at the top of the script.  
2. Normalization: Calculates size factors based on total raw insertions per group and normalizes the count matrices from Step 3\.  
3. Differential Analysis: Uses DESeq2 to perform pairwise differential accessibility analysis. It runs two types of comparisons:  
   * Each experimental group vs. a designated control\_group.  
   * Each experimental group vs. all other groups to find group-specific peaks.  
4. Peak Filtering & Merging: Filters for significantly differential peaks based on log2FoldChange and normalized count cutoffs. It then merges overlapping significant peaks (vs. control) to create a consensus peak set.  
5. Heatmap Visualization: Generates Z-score scaled heatmaps of normalized counts for both the "group-specific" and the "vs. Control" peak sets using the ComplexHeatmap package.  
6. Summary Plots: Creates bar plots summarizing total insertions and the number of significant peaks per group.  
7. HOMER Integration (Manual Step):  
   * The script generates BED files for significant peaks and prints the annotatePeaks.pl (for annotation) and findMotifsGenome.pl (for motif discovery) commands to the console.  
   * The user must copy these commands and run them in a terminal where HOMER is installed.  
8. Post-HOMER Analysis: Once the HOMER jobs are complete, the R script can be re-run or continued to:  
   * Parse the HOMER annotation files to plot the genomic distribution of peaks (e.g., promoter, intron, intergenic).  
   * Extract lists of promoter-bound genes.  
   * Parse the HOMER motif analysis results to plot the top enriched transcription factor motifs for each group.

#### **Inputs**

* This script does not take command-line arguments. The user must manually edit the paths in the PART 0: CONFIGURATION section of the script.  
* It requires the output directories from Step 1 (raw\_unique\_insertion\_count\_per\_group\_bedgraph/) and Step 3 (per\_group\_peak\_matrices/).

#### **Outputs**

A master output directory (STEP\_4\_output\_deseq2\_peak\_analysis/) containing:

* DESeq2\_results\_by\_group/: Detailed DESeq2 output for every pairwise comparison.  
* COMBINED\_... .csv: Tables consolidating all differential analysis results.  
* merged\_significant\_peaks\_vs\_HyPBase\_matrix.csv: The final count matrix for consensus peaks significantly enriched over the control group.  
* group\_specific\_peaks\_summary.csv: The final count matrix for peaks found to be specific to each experimental group.  
* heatmap\_... .png/pdf: Heatmap visualizations.  
* Various bar plots (.png/.pdf) summarizing the data.  
* HOMER/: A subdirectory containing:  
  * Per\_Group\_BEDs/: BED files used as input for HOMER.  
  * Per\_Group\_Annotations/: The raw output from annotatePeaks.pl.  
  * HOMER\_motif\_discovery/: The raw output from findMotifsGenome.pl.  
  * Plots summarizing genomic distribution and top enriched motifs.  
  * Lists of promoter-bound genes.

#### **Usage**

1. Configure the Script: Open step\_4\_... .R and modify the paths and group definitions in PART 0.  
2. Run in R: Execute the script in an R environment.  
3. Run HOMER: The script will pause after generating HOMER commands. Copy these commands into your terminal, run them, and wait for them to complete.  
4. Complete R Script: Once HOMER is finished, the R script can continue to process the results and generate the final plots.

