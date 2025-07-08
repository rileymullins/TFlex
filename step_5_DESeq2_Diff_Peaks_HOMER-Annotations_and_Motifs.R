# step_5_DESeq2_Diff_Peaks_HOMER-Annotations_and_Motifs.R


# ==============================================================================
# PART 0: CONFIGURATION
# ==============================================================================

# ==== 1. Load Libraries ====
suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(annotatr)
  library(biomaRt)
  library(circlize)
  library(ComplexHeatmap)
  library(data.table)
  library(DESeq2)
  library(dplyr)
  library(forcats)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(ggplot2)
  library(ggrepel)
  library(grid)
  library(parallel)
  library(purrr)
  library(readr)
  library(RColorBrewer)
  library(stringr)
  library(tibble)
  library(tidyr)
})

# ==== 2. Define Inputs & Experimental Groups ====
# Point to the final output file from the Python script.
input_file <- "/path/to/your/output_dir/final_results.tsv"
annotation_file <- "/path/to/your/annotation_file"
output_dir <- "/path/to/your/output_dir/analysis_plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Manually define the groups for the analysis.
# The names here (e.g., "HyPBase") MUST EXACTLY MATCH the values in the 'group_name' column of your final_results and annotation_file

# Print the possible groups to choose from.
annotation_dataframe <- read.csv(annotation_file)
cat("Groups to choose from:\n", paste(unique(annotation_dataframe$group_name), collapse = "\n"))

# This group_map should be all groups that you want to compare against background control (HyPBase), compare all pairwise against each other,
 # and compare each against all others including the background control (HyPBase) to find group-specific peaks.
group_map <- c(
 "group_1",
 "group_2", 
 "etc."
)

names(group_map) <- group_map # Name elements by themselves for compatibility.

# Define experiment and control groups from your manual list.
control_group <- "HyPBase" # This is the background control group. It will always be HyPBase. Could be HyPBase, HyPBase_timepointA, HyPBase_treatmentX, etc.
experimental_groups <- setdiff(group_map, control_group)


###############################################################################################
#### SKIP TO 'PART 5: LOAD EXISTING RESULTS FROM PRIOR STEPS' IF 1-4 HAVE ALREADY BEEN RUN ####
###############################################################################################




# ==== 3. Data Loading and Matrix Creation ====
# Read the `final_results.tsv` file.
full_data <- fread(input_file)

# Filter the data to ONLY include the groups specified in your manual group_map.
full_data <- full_data %>% filter(group_name %in% group_map)

# CREATE COUNT DATA: Use the group_total_insertions_in_consensus_peak column.
count_data <- full_data %>%
  select(consensus_peak_id, group_name, group_total_insertions_in_consensus_peak) %>%
  distinct(consensus_peak_id, group_name, .keep_all = TRUE)

# CREATE METADATA TO JOIN: Select relevant group-level statistics.
metadata_to_join <- full_data %>%
  select(consensus_peak_id, group_name, num_samples_in_consensus_peak,
         group_total_fragments_in_consensus_peak, group_total_reads_in_consensus_peak) %>%
  distinct(consensus_peak_id, group_name, .keep_all = TRUE)


# CREATE THE COUNT MATRIX
count_matrix <- dcast(count_data, consensus_peak_id ~ group_name,
                      value.var = "group_total_insertions_in_consensus_peak", fill = 0)

count_matrix_df <- as.data.frame(count_matrix)
rownames(count_matrix_df) <- count_matrix_df$consensus_peak_id
count_matrix_df$consensus_peak_id <- NULL

# Ensure all groups from the manual map are present as columns and enforce order.
missing_groups <- setdiff(group_map, colnames(count_matrix_df))
for (group in missing_groups) {
    count_matrix_df[[group]] <- 0
}
count_matrix_df <- count_matrix_df[, group_map]

count_matrix_df[] <- lapply(count_matrix_df, as.integer)
cat("Count matrix created with", nrow(count_matrix_df), "peaks and", ncol(count_matrix_df), "groups.\n\n")

# ==== 4. Global Normalization ====
# The column names of count_matrix_df represent the sample or group names.
# We create a data.frame where:
#   row.names = sample names (colnames of the count matrix)
#   condition = condition is set to be the same as the sample/group name (this is arbitrary here, since differential analysis is not being done here)
# This allows DESeq2 to process the samples and assign size factors for normalization.
col_data_all <- data.frame(
  row.names = colnames(count_matrix_df),      # Set sample names as row names
  condition = colnames(count_matrix_df)       # A dummy condition (each sample its own "condition")
)

# countData = count_matrix_df: the raw count matrix (rows = consensus peaks, columns = samples)
# colData = col_data_all: sample metadata (created above)
# design = ~ 1:
#   No groups or conditions, just an intercept.
#   This means DESeq2 does not model any condition effect and only will normalize counts.
dds_all <- DESeqDataSetFromMatrix(
  countData = count_matrix_df,
  colData = col_data_all,
  design = ~ 1
)

# DESeq2 computes size factors to account for differences in total unique insertions
# These factors are used to scale counts so groups can be fairly compared.
# Size factors are stored inside dds_all.
dds_all <- estimateSizeFactors(dds_all)
sizeFactors(dds_all) # Print size factors to confirm they look right.

# Generate a normalized counts matrix to save
normalized_counts <- counts(dds_all, normalized = TRUE)

# Save file
norm_counts_file <- file.path(output_dir, "normalized_srt_barcode_counts_per_group_in_consensus_peaks.csv")
fwrite(as.data.table(normalized_counts, keep.rownames = "consensus_peak_id"), file = norm_counts_file)
cat("Normalized SRT barcode counts for groups in consensus peaks saved to:", norm_counts_file, "\n")

# Calculate normalized counts per kb
peak_widths <- full_data %>%
  select(consensus_peak_id, consensus_peak_width) %>%
  distinct() %>%
  mutate(peak_width_kb = consensus_peak_width / 1000)

peak_widths_ordered <- peak_widths[match(rownames(normalized_counts), peak_widths$consensus_peak_id), ]
normalized_counts_per_kb <- normalized_counts / peak_widths_ordered$peak_width_kb

norm_counts_pkb_file <- file.path(output_dir, "normalized_srt_barcode_counts_per_kb_per_group_in_consensus_peaks.csv")
fwrite(as.data.table(normalized_counts_per_kb, keep.rownames = "consensus_peak_id"), file = norm_counts_pkb_file)
cat("Normalized SRT barcode counts per kb for groups in consensus peaks saved to:", norm_counts_pkb_file, "\n\n")





# ==============================================================================
# PART 1: GROUP vs. CONTROL ANALYSIS (ONE TAIL)
# ==============================================================================
all_control_results_list <- list()

for (group in experimental_groups) {
  cat("------------------------------------------------------------\n")
  cat("Comparing:", group, "vs.", control_group, "\n")
  cat("------------------------------------------------------------\n")

  pairwise_matrix <- count_matrix_df[, c(control_group, group)] # Matrix of HyPBase and the group being processed in the loop.

  # Identify rows to keep: at least one count > 0 in either the control OR the treatment group.
  # This omits peaks that have no signal in any of the samples for this specific comparison.
  keep_rows <- (rowSums(pairwise_matrix[, control_group, drop=FALSE], na.rm = TRUE) > 0) |  
               (rowSums(pairwise_matrix[, group, drop=FALSE], na.rm = TRUE) > 0)

  # Filter the matrix to keep only the identified rows.
  pairwise_matrix <- pairwise_matrix[keep_rows, ]

  # Generate col_data object for DESeq2
  col_data <- data.frame(condition = factor(c("control", "treatment")), row.names = c(control_group, group))

  dds <- DESeqDataSetFromMatrix(countData = pairwise_matrix, colData = col_data, design = ~ condition) # Make DESeq dataset from raw counts matrix.
  dds <- estimateSizeFactors(dds) # Determine size factors for the comparison.

  # Dispersion calculation
   # Make a copy of the main object.
   dds_for_disp <- dds
   # Calculate dispersion. The ~ 1 means to treat all samples as if they belong to a single group.
   # This effectively does the following:
    # (1) Finds the formula that models dispersion in both datasets in the comparison and fits a curve to model that. The plot is dispersion on Y axis, total counts on X per peak.
    # (2) Raise data points below the line up to the fitted line, increasing their dispersion.
    # (3) Leave data points above the line where they are so that the dispersion is not decreased.

   # DESeq2 will give lower log2FC values in situations when the group has 10 insertions in the peak and HyPBase has 0, so 10 v 0 will have a lower log2FC than 100 v 0.
   design(dds_for_disp) <- formula(~ 1)
   dds_for_disp <- estimateDispersions(dds_for_disp, fitType = "local")
   dispersions(dds) <- dispersions(dds_for_disp)

  # Run test to find peaks that are greater in the 'treatment' (group of experimental_group set of loop) relative to 'control' (HyPBase).
  # One tailed ('greater') because the question is is the group higher than HyPBase.
  dds <- nbinomWaldTest(dds)
  res <- results(dds, contrast = c("condition", "treatment", "control"), altHypothesis = "greater")

  # Make a plot of the dispersion for reference of what DESeq2 is doing for dispersion.
  disp_file <- file.path(output_dir, paste0(group, "_vs_", control_group, "_diagnostic_plot_dispersion_estimates.png"))
  png(disp_file, width = 6, height = 6, units = "in", res = 300)
  plotDispEsts(dds_for_disp)
  dev.off()

  # Make a dataframe out of the results
  # Note that the p value is meaningless because the comparison does not use replicates.
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("consensus_peak_id") %>%
    mutate(comparison = paste0(group, "_vs_", control_group))

  # Generate objects to add the the metadata in the final_output.tsv file to the results dataframe
  group_metadata <- metadata_to_join %>% filter(group_name == group) %>% select(-group_name) %>% rename_with(~ paste0(., "_", group), -consensus_peak_id)
  control_metadata <- metadata_to_join %>% filter(group_name == control_group) %>% select(-group_name) %>% rename_with(~ paste0(., "_", control_group), -consensus_peak_id)

  # Generate object to add the normalized counts in the consensus peak to results dataframe
  norm_counts_to_add <- as.data.frame(normalized_counts) %>%
    tibble::rownames_to_column("consensus_peak_id") %>%
    select(consensus_peak_id, all_of(c(group, control_group))) %>%
    rename_with(~ paste0("norm_count_", .), -consensus_peak_id)

  # Generate object to add the normalized counts per kb of consensus peak to results dataframe
  norm_counts_pkb_to_add <- as.data.frame(normalized_counts_per_kb) %>%
    tibble::rownames_to_column("consensus_peak_id") %>%
    select(consensus_peak_id, all_of(c(group, control_group))) %>%
    rename_with(~ paste0("norm_count_pkb_", .), -consensus_peak_id)

  # Generate the final results dataframe that has the metadata and counts for each comparison.
  res_df_with_meta <- res_df %>%
    left_join(group_metadata, by = "consensus_peak_id") %>%
    left_join(control_metadata, by = "consensus_peak_id") %>%
    left_join(norm_counts_to_add, by = "consensus_peak_id") %>%
    left_join(norm_counts_pkb_to_add, by = "consensus_peak_id") %>%
    arrange(desc(log2FoldChange), pvalue)

  # Store in a list
  all_control_results_list[[group]] <- res_df_with_meta

  # Output the results dataframe.
  output_file_full <- file.path(output_dir, paste0(group, "_vs_", control_group, "_full_results_one_tailed.csv"))
  fwrite(res_df_with_meta, file = output_file_full, na = "NA")
  cat("Full one-tailed results saved to:", output_file_full, "\n")
}
# Note this will generate a warning: 'In qf(0.99, p, m - p) : NaNs produced'
# Reason is that m = total samples = 1 control + 1 treatment = 2 and p = model parameters = 1 intercept + 1 condition effect = 2.
# qf() = function to compute quantiles of the F-distribution
# Therefore, qf(0.99, p, m - p) -> qf(0.99, p = 2, (m = 2) - (p = 2)) -> qf(0.99, 2, 0). 
# This cannot be computed, so NAs are producted








# ==============================================================================
# PART 2: ALL-vs-ALL PAIRWISE GROUP ANALYSIS (TWO TAIL)
# SEE PART 1 FOR FULL DETAILS ON THE RATIONALE OF THE DESEQ2 STRATEGY.
# ==============================================================================

cat("PART 2: Starting All-vs-All Group pairwise analysis...\n")
all_pairwise_results_list <- list()
group_vs_group_pairs <- combn(experimental_groups, 2, simplify = FALSE)

for (pair in group_vs_group_pairs) {
  group1 <- pair[1]
  group2 <- pair[2]

  cat("------------------------------------------------------------\n")
  cat("Comparing:", group2, "vs.", group1, "\n")
  cat("------------------------------------------------------------\n")

  pairwise_matrix <- count_matrix_df[, c(group1, group2)]
  pairwise_matrix <- pairwise_matrix[rowSums(pairwise_matrix) > 0, ] # Keep peaks where there is at least 1 count for a group. Since we are doing TF v all, this will be all rows.
  col_data <- data.frame(condition = factor(c("group1", "group2")), row.names = c(group1, group2))

  dds <- DESeqDataSetFromMatrix(countData = pairwise_matrix, colData = col_data, design = ~ condition)
  dds <- estimateSizeFactors(dds)
  dds_for_disp <- dds
  design(dds_for_disp) <- formula(~ 1)
  dds_for_disp <- estimateDispersions(dds_for_disp, fitType = "local")
  dispersions(dds) <- dispersions(dds_for_disp)
  dds <- nbinomWaldTest(dds)
  res <- results(dds, contrast = c("condition", "group2", "group1"))

  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("consensus_peak_id") %>%
    mutate(comparison = paste0(group2, "_vs_", group1))

  group1_metadata <- metadata_to_join %>% filter(group_name == group1) %>% select(-group_name) %>% rename_with(~ paste0(., "_", group1), -consensus_peak_id)
  group2_metadata <- metadata_to_join %>% filter(group_name == group2) %>% select(-group_name) %>% rename_with(~ paste0(., "_", group2), -consensus_peak_id)

  norm_counts_to_add <- as.data.frame(normalized_counts) %>%
    tibble::rownames_to_column("consensus_peak_id") %>%
    select(consensus_peak_id, all_of(c(group1, group2))) %>%
    rename_with(~ paste0("norm_count_", .), -consensus_peak_id)

  norm_counts_pkb_to_add <- as.data.frame(normalized_counts_per_kb) %>%
    tibble::rownames_to_column("consensus_peak_id") %>%
    select(consensus_peak_id, all_of(c(group1, group2))) %>%
    rename_with(~ paste0("norm_count_pkb_", .), -consensus_peak_id)

  res_df_with_meta <- res_df %>%
    left_join(group1_metadata, by = "consensus_peak_id") %>%
    left_join(group2_metadata, by = "consensus_peak_id") %>%
    left_join(norm_counts_to_add, by = "consensus_peak_id") %>%
    left_join(norm_counts_pkb_to_add, by = "consensus_peak_id") %>%
    arrange(desc(abs(log2FoldChange)))

  all_pairwise_results_list[[paste0(group2, "_vs_", group1)]] <- res_df_with_meta

  output_file_pair <- file.path(output_dir, paste0(group2, "_vs_", group1, "_pairwise_results.csv"))
  fwrite(res_df_with_meta, file = output_file_pair, na = "NA")
  cat("Pairwise results saved to:", output_file_pair, "\n")
}







# ==============================================================================
# PART 3: GROUP vs. ALL OTHER GROUPS ANALYSIS (ONE TAIL)
# SEE PART 1 FOR FULL DETAILS ON THE RATIONALE OF THE DESEQ2 STRATEGY.
# ==============================================================================

all_vs_all_others_results_list <- list()
experimental_groups2 <- c(experimental_groups, control_group) # This is set so that the comparison is the group vs. all other groups including HyPBase.

for (group1 in experimental_groups2) {
  other_groups <- setdiff(experimental_groups2, group1)
  comparison_name <- paste0(group1, "_vs_ALL_OTHERS")
  cat("------------------------------------------------------------\n")
  cat("Comparing:", comparison_name, "\n")
  cat("------------------------------------------------------------\n")

  group1_counts <- count_matrix_df[, group1]
  other_groups_counts <- if (length(other_groups) > 1) rowSums(count_matrix_df[, other_groups]) else count_matrix_df[, other_groups] # sum the counts of the other groups

  vs_all_matrix <- data.frame(group1_counts, other_groups_counts)
  colnames(vs_all_matrix) <- c(group1, "ALL_OTHERS")
  vs_all_matrix <- vs_all_matrix[rowSums(vs_all_matrix) > 0, ] # Keep peaks where there is at least 1 count for a group. Since we are doing TF v all, this will be all rows.
  col_data <- data.frame(condition = factor(c("group1", "all_others")), row.names = c(group1, "ALL_OTHERS"))

  dds <- DESeqDataSetFromMatrix(countData = vs_all_matrix, colData = col_data, design = ~ condition)
  dds <- estimateSizeFactors(dds)
  dds_for_disp <- dds
  design(dds_for_disp) <- formula(~ 1)
  dds_for_disp <- estimateDispersions(dds_for_disp, fitType = "local")
  dispersions(dds) <- dispersions(dds_for_disp)
  dds <- nbinomWaldTest(dds)
  res <- results(dds, contrast = c("condition", "group1", "all_others"), altHypothesis = "greater")

  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("consensus_peak_id") %>%
    mutate(comparison = comparison_name)

  group1_metadata <- metadata_to_join %>% filter(group_name == group1) %>% select(-group_name) %>% rename_with(~ paste0(., "_", group1), -consensus_peak_id)

  norm_counts_to_add <- as.data.frame(normalized_counts) %>%
    tibble::rownames_to_column("consensus_peak_id") %>%
    select(consensus_peak_id, all_of(group1)) %>%
    rename_with(~ paste0("norm_count_", .), -consensus_peak_id)

  norm_counts_pkb_to_add <- as.data.frame(normalized_counts_per_kb) %>%
    tibble::rownames_to_column("consensus_peak_id") %>%
    select(consensus_peak_id, all_of(group1)) %>%
    rename_with(~ paste0("norm_count_pkb_", .), -consensus_peak_id)

  res_df_with_meta <- res_df %>%
    left_join(group1_metadata, by = "consensus_peak_id") %>%
    left_join(norm_counts_to_add, by = "consensus_peak_id") %>%
    left_join(norm_counts_pkb_to_add, by = "consensus_peak_id") %>%
    arrange(desc(log2FoldChange))

  all_vs_all_others_results_list[[comparison_name]] <- res_df_with_meta

  output_file_vs_all <- file.path(output_dir, paste0(comparison_name, "_results.csv"))
  fwrite(res_df_with_meta, file = output_file_vs_all, na = "NA")
  cat("Group vs ALL OTHERS results saved to:", output_file_vs_all, "\n")
}







# ==============================================================================
# PART 4: SAVE COMBINED OUTPUTS
# ==============================================================================
combined_control_results <- bind_rows(all_control_results_list)
fwrite(combined_control_results, file = file.path(output_dir, "COMBINED_all_Group_vs_Control_results.csv"), na = "NA")
cat("Combined Group vs. Control results table saved.\n")

combined_pairwise_results <- bind_rows(all_pairwise_results_list)
fwrite(combined_pairwise_results, file = file.path(output_dir, "COMBINED_all_Group_vs_Group_results.csv"), na = "NA")
cat("Combined All-vs-All Group results table saved.\n")

combined_vs_all_others_results <- bind_rows(all_vs_all_others_results_list)
fwrite(combined_vs_all_others_results, file = file.path(output_dir, "COMBINED_Group_vs_ALL_OTHERS_results.csv"), na = "NA")
cat("Combined Group vs. ALL OTHERS results table saved.\n\n")







# ==============================================================================
# PART 5: LOAD EXISTING RESULTS FROM PRIOR STEPS
# ==============================================================================

combined_control_file <- file.path(output_dir, "COMBINED_all_Group_vs_Control_results.csv")
combined_pairwise_file <- file.path(output_dir, "COMBINED_all_Group_vs_Group_results.csv")
combined_vs_all_file <- file.path(output_dir, "COMBINED_Group_vs_ALL_OTHERS_results.csv")
norm_counts_pkb_file <- file.path(output_dir, "normalized_srt_barcode_counts_per_kb_per_group_in_consensus_peaks.csv")
norm_counts_file <- file.path(output_dir, "normalized_srt_barcode_counts_per_group_in_consensus_peaks.csv")

if (file.exists(combined_control_file) && file.exists(combined_pairwise_file) && file.exists(combined_vs_all_file) && file.exists(norm_counts_pkb_file) && file.exists(norm_counts_file)) {
  combined_control_results <- fread(combined_control_file)
  combined_pairwise_results <- fread(combined_pairwise_file)
  combined_vs_all_others_results <- fread(combined_vs_all_file)
  normalized_counts_per_kb <- fread(norm_counts_pkb_file)
  normalized_counts <- fread(norm_counts_file)
  cat("Successfully loaded all combined result files.\n\n")
} else {
  stop("Result files not found. Run the analysis parts first.")
}








# ==============================================================================
# PART 6: VISUALIZATIONS/PLOTS
# ==============================================================================
# The norm_count_pkb and log2 fold change cutoffs are used to filter to high-confidence sites where the peak is supported by a lot of insertions.
# It's best to look at the DESeq2 results and understand the relation of the normalized counts and Log2FC to determine the best
      # combination of log2FC and counts cutoffs to get high confidence sites.
# Can also use total counts.

# Define cutoffs. Log2FC cutoff 1 and norm count per kb cutoff 20 is usually optimal.
norm_count_pkb_cutoff = 20
log2FoldChange_cutoff = 1

# Set color_scale  
color_scale  <- colorRamp2(
    breaks = seq(from = q_low_specific, to = q_high_specific, length.out = 100),
    colors = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4")))(100)
  )


# -- Prepare Data. This only needs to be ran once. --
# vs Control Group
top_control_hits_list <- combined_control_results %>%
  filter(log2FoldChange >= log2FoldChange_cutoff) %>%
  rowwise() %>%
  mutate(group = str_remove(comparison, "_vs_.*")) %>%
  filter(get(paste0("norm_count_pkb_", group)) >= norm_count_pkb_cutoff) %>%
  ungroup() %>%
  arrange(comparison, desc(log2FoldChange))

unique_top_control_peaks <- unique(top_control_hits_list$consensus_peak_id)
cat(length(unique_top_control_peaks), "unique top peaks from Group vs. Control collected for heatmap.\n")

# Group v Group pairwise (see above for comment on norm_count_pkb cutoff)
top_pairwise_hits_list <- combined_pairwise_results %>%
  filter(abs(log2FoldChange) >= log2FoldChange_cutoff) %>%
  rowwise() %>%
  mutate(
    group = str_remove(comparison, "_vs_.*"),
    control = str_remove(comparison, ".*_vs_")
  ) %>%
  filter(
    get(paste0("norm_count_pkb_", group)) >= norm_count_pkb_cutoff |
    get(paste0("norm_count_pkb_", control)) >= norm_count_pkb_cutoff
  ) %>%
  ungroup() %>%
  arrange(comparison, desc(abs(log2FoldChange)))

unique_top_pairwise_peaks <- unique(top_pairwise_hits_list$consensus_peak_id)
cat(length(unique_top_pairwise_peaks), "unique top peaks from All-vs-All Group collected for heatmap.\n")

# Group vs all others (see above for comment on norm_count_pkb cutoff)
top_group_v_all_hits_list <- combined_vs_all_others_results %>%
  filter(log2FoldChange >= log2FoldChange_cutoff) %>%
  rowwise() %>%
  mutate(group = str_remove(comparison, "_vs_.*")) %>%
  filter(!group %in% control_group) %>%
  filter(get(paste0("norm_count_pkb_", group)) >= norm_count_pkb_cutoff) %>%
  ungroup() %>%
  arrange(comparison, desc(abs(log2FoldChange)))

union_top_group_v_all_hits <- unique(top_group_v_all_hits_list$consensus_peak_id)
cat(length(union_top_group_v_all_hits), "unique top peaks from Group vs. All collected for heatmap.\n")


# Read in data to generate heatmaps of norm_counts_pkb.
# Plotting the per kb normalization helps to control for large differences in peak size with a higher chance of harboring insertions.
# The values are z scored across groups per peak, so it's not so important here.
# If interested in relative absolute signals, norm_counts_pkb is the best choice.

normalized_counts_per_kb <- as.data.frame(fread(norm_counts_pkb_file))
rownames(normalized_counts_per_kb) <- normalized_counts_per_kb$consensus_peak_id
normalized_counts_per_kb <- normalized_counts_per_kb[, -1]



# Make heatmaps of log2 scaled normalized counts per kb among the peaks.
# -- Heatmap 1: From Group vs. Control Results --
if (length(unique_top_control_peaks) > 1) {
  
  # Generate log2 scaled matrix of normalized counts per kb
  heatmap_matrix_control <- normalized_counts_per_kb[unique_top_control_peaks, ]
  heatmap_matrix_log_scaled_control <- t(scale(t(log2(heatmap_matrix_control + 1))))
  
  # Define 5th and 95th quantile for color scale
  q_low_control <- quantile(heatmap_matrix_log_scaled_control, 0.05, na.rm = TRUE)
  q_high_control <- quantile(heatmap_matrix_log_scaled_control, 0.95, na.rm = TRUE)
  
  # Generate heatmap
  ht_control <- Heatmap(
    heatmap_matrix_log_scaled_control,
    name = "Z-score",
    col = color_scale,
    row_km = 8, # Adjust Km to find what captures the patterns the best
    cluster_columns = TRUE,
    cluster_rows = FALSE,
    show_row_names = FALSE,
    row_title = NULL,
    row_dend_reorder = T,
    column_dend_reorder = T,
    border ="black",
    column_title = paste("Peaks Above", control_group,
                         "\n(ins/kb>",norm_count_pkb_cutoff, "log2fc>",log2FoldChange_cutoff,")"),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),

    heatmap_legend_param = list(title = "Z-score",
                                title_gp = gpar(fontsize = 8),
                                labels_gp = gpar(fontsize = 8)),
    row_names_gp = gpar(fontsize = 9),
    column_names_gp = gpar(fontsize = 9)
  )

  # Save heatmap
  heatmap_file_control <- file.path(output_dir, "heatmap_top_vs_Control_peaks.png")
  png(heatmap_file_control, width = 3, height = 8, units = "in", res = 300)
  draw(ht_control)
  dev.off()
  cat("Heatmap for Group vs. Control results saved to:", heatmap_file_control, "\n")
} else {
  cat("Not enough unique differential peaks found in Group vs Control comparisons to generate a heatmap.\n\n")
}

# -- Heatmap 2: From All-vs-All Group Pairwise Results --
if (length(unique_top_pairwise_peaks) > 1) {
  
  # Generate log2 scaled matrix of normalized counts per kb
  heatmap_matrix_pairwise <- normalized_counts_per_kb[unique_top_pairwise_peaks, experimental_groups]
  heatmap_matrix_log_scaled_pairwise <- t(scale(t(log2(heatmap_matrix_pairwise + 1))))
  
  # Define 5th and 95th quantile for color scale
  q_low_pairwise <- quantile(heatmap_matrix_log_scaled_pairwise, 0.05, na.rm = TRUE)
  q_high_pairwise <- quantile(heatmap_matrix_log_scaled_pairwise, 0.95, na.rm = TRUE)

  # Generate heatmap
  ht_pairwise <- Heatmap(
    heatmap_matrix_log_scaled_pairwise,
    name = "Z-score",
    col = color_scale,
    row_km = 10, # Adjust Km to find what captures the patterns the best.
    border ="black", # Border around the Km blocks of peaks
    cluster_columns = TRUE,
    show_row_names = FALSE,
    cluster_rows = F, # Row Km clustering only
    row_title = NULL,
    row_dend_reorder = T,
    column_dend_reorder = T,
    column_title = paste("Pairwise Group vs Group Peaks",
                         "\n(ins/kb>",norm_count_pkb_cutoff, "log2fc>",log2FoldChange_cutoff,")"),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    heatmap_legend_param = list(title = "Z-score",
                                title_gp = gpar(fontsize = 8),
                                labels_gp = gpar(fontsize = 8)),
    row_names_gp = gpar(fontsize = 9),
    column_names_gp = gpar(fontsize = 9)
  )
  
  # Save heatmap
  heatmap_file_pairwise <- file.path(output_dir, "heatmap_top_Group_vs_Group_peaks.png")
  png(heatmap_file_pairwise, width = 3, height = 8, units = "in", res = 300)
  draw(ht_pairwise)
  dev.off()
  cat("Heatmap for All-vs-All Group results saved to:", heatmap_file_pairwise, "\n\n")
} else {
  cat("Not enough unique differential peaks found in pairwise comparisons to generate a heatmap.\n\n")
}

# -- Heatmap 3: From Group vs All Others Results --
if (length(union_top_group_v_all_hits) > 1) {
  
  # Generate log2 scaled matrix of normalized counts per kb
  heatmap_matrix_specific <- normalized_counts_per_kb[union_top_group_v_all_hits, experimental_groups]
  heatmap_matrix_log_scaled_specific <- t(scale(t(log2(heatmap_matrix_specific + 1))))
  
  # Define 5th and 95th quantile for color scale
  q_low_specific <- quantile(heatmap_matrix_log_scaled_specific, 0.05, na.rm = TRUE)
  q_high_specific <- quantile(heatmap_matrix_log_scaled_specific, 0.95, na.rm = TRUE)
  
  # Generate heatmap
  ht_specific <- Heatmap(
    heatmap_matrix_log_scaled_specific,
    name = "Z-score",
    col = color_scale,
    row_km = 9, # Adjust Km to find what captures the patterns the best.
    border ="black", # Border around the Km blocks of peaks
    cluster_columns = TRUE,
    show_row_names = FALSE,
    row_dend_reorder = T,
    column_dend_reorder = T,
    cluster_rows = F,
    row_title = NULL,
    column_title = paste("Group-Specific Peaks",
                         "\n(ins/kb>",norm_count_pkb_cutoff, "log2fc>",log2FoldChange_cutoff,")"),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    heatmap_legend_param = list(title = "Z-score",
                                title_gp = gpar(fontsize = 8),
                                labels_gp = gpar(fontsize = 8)),
    row_names_gp = gpar(fontsize = 9),
    column_names_gp = gpar(fontsize = 9)
  )
  
  # Save heatmap
  heatmap_file_specific <- file.path(output_dir, "heatmap_top_Group_specific_peaks.png")
  png(heatmap_file_specific, width = 3, height = 8, units = "in", res = 300)
  draw(ht_specific)
  dev.off()
  cat("Heatmap for Group-specific results saved to:", heatmap_file_specific, "\n\n")
} else {
  cat("Not enough unique differential peaks found in group-specific comparisons to generate a heatmap.\n\n")
}

# ==============================================================================
# PART 7: SUMMARY PLOTS
# ==============================================================================
# ==== Define the order of groups for plots and their colors ====
group_order <- c(
  "HyPBase",
  "JUN",
  "BATF",
  "FOSL1",
  "JUN.FOS.BATF",
  "MAFF.ATF3.BATF2",
  "BATF3.JDP2.BATF2",
  "BATF3.BATF2.MAFF",
  "FOSL1.BATF.MAFA"
)

# Create a palette matching the number of groups, named by group_order
palette_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(group_order))
names(palette_colors) <- group_order

# ==== Bar Plot of Total Unique Insertions ====
bar_data_insertions <- full_data %>%
  group_by(group_name) %>%
  distinct(consensus_peak_id, .keep_all = T) %>%
  summarise(unique_insertions = sum(group_total_insertions_in_consensus_peak, na.rm = TRUE)) %>%
  ungroup() %>%
  # ensure only groups in our defined order are plotted
  filter(group_name %in% group_order) %>%
  mutate(group_name = factor(group_name, levels = group_order))

p_insertions <- ggplot(bar_data_insertions,
                       aes(x = group_name, y = unique_insertions, fill = group_name)) +
  geom_col(color = "black", linewidth = 0.1, width = 0.75) +
  geom_text(aes(label = unique_insertions),
            vjust = -0.4, size = 2, color = "black") +
  scale_fill_manual(values = palette_colors, guide = FALSE) +
  labs(x = 'Group', y = 'Total Unique Insertions') +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 8),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.2)
  )
p_insertions
ggsave(file.path(output_dir, 'unique_insertions_per_group_barplot.png'),
       p_insertions, width = 5, height = 5, dpi = 300, bg = "white")
cat("Total insertions bar plot saved.\n")

# ==== Bar Plot of Total Significant Peaks (No Bins) ====
bar_data_simple <- top_control_hits_list %>%
  count(group, name = "num_peaks") %>%
  mutate(group = factor(group, levels = group_order))

p_simple <- ggplot(bar_data_simple,
                   aes(x = group, y = num_peaks, fill = group)) +
  geom_col(color = "black", linewidth = 0.1, width = 0.75) +
  geom_text(aes(label = num_peaks),
            vjust = -0.3, size = 3, color = "black") +
  scale_fill_manual(values = palette_colors, guide = FALSE) +
  labs(x = "Group",
       y = "Number of Significant Peaks\n(log2FC ≥ cutoff, Norm Ins/kb ≥ cutoff)") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.2)
  )
p_simple
ggsave(file.path(output_dir, "peaks_per_group_barplot_simple.png"),
       p_simple, width = 5, height = 5, dpi = 300, bg = "white")
cat("Simple bar plot saved to:\n", file.path(output_dir, "peaks_per_group_barplot_simple.png"), "\n")


# ==== Bar Plot with Log2FC Bins ====
bar_data_binned <- top_control_hits_list %>%
  mutate(
    log2fc_bin = case_when(
      log2FoldChange >= 1 & log2FoldChange < 2 ~ "1-2",
      log2FoldChange >= 2 & log2FoldChange < 3 ~ "2-3",
      log2FoldChange >= 3 & log2FoldChange < 4 ~ "3-4",
      log2FoldChange >= 4 & log2FoldChange < 5 ~ "4-5",
      log2FoldChange >= 5 & log2FoldChange < 6 ~ "5-6",
      log2FoldChange >= 6 ~ "6+",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(log2fc_bin)) %>%
  count(group, log2fc_bin, name = "num_peaks") %>%
  mutate(group = factor(group, levels = group_order))

bar_data_binned$log2fc_bin <- factor(bar_data_binned$log2fc_bin, levels = rev(c("1-2", "2-3", "3-4", "4-5", "5-6", "6+")))

bin_colors <- c("1-2" = "#8dd3c7", "2-3" = "#ffffb3", "3-4" = "#bebada", "4-5" = "#fb8072", "5-6" = "#80b1d3", "6+" = "#fdb462")

p_binned <- ggplot(bar_data_binned, aes(x = group, y = num_peaks, fill = log2fc_bin)) +
  geom_col(color = "black", linewidth = 0.1, width = 0.75) +
  geom_text(aes(label = num_peaks), position = position_stack(vjust = 0.5), size = 2.5, fontface = "plain", color = "black") +
  scale_fill_manual(values = bin_colors, name = "log2FC Bin", drop = FALSE) +
  labs(x = 'Group', y = 'Number of Significant Peaks\n(Binned log2FC, Norm Ins/kb > 20)') +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "plain", color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10), panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.2)
  )
bar_data_binned
bar_file_binned <- file.path(output_dir, 'peaks_per_group_barplot_log2fc_bins.png')
ggsave(bar_file_binned, p_binned, width = 7, height = 5, dpi = 300, bg = "white")
cat("Binned bar plot saved.\n")



# ==========================================================================================================
# Define Significant Peaks (used in both PART 8 and PART 9)

# IMPORTANT - Make sure HOMER is installed and HOMER's bin directory is in your PATH environment variable.
# ==========================================================================================================
# This object has every peak that is significantly above background for at least one group.
significant_peak_events <- combined_control_results %>%
  as_tibble() %>%
  filter(log2FoldChange >= log2FoldChange_cutoff) %>%
  pivot_longer(starts_with("norm_count_pkb_"), names_to = "norm_col", values_to = "norm_val") %>%
  mutate(
    group_name = sub("^norm_count_pkb_", "", norm_col),
    group_from_comparison = sub(paste0("_vs_", control_group, "$"), "", comparison)
  ) %>%
  filter(group_name == group_from_comparison, !is.na(norm_val), norm_val >= norm_count_pkb_cutoff) %>%
  select(consensus_peak_id, group_name) %>%
  distinct()

cat("Identified", nrow(significant_peak_events), "total binding events above background across all groups.\n")

# This object is the union of all significant peaks, for annotation and for the motif background.
all_significant_peaks_df <- significant_peak_events %>%
  distinct(consensus_peak_id) %>%
  separate(consensus_peak_id, into = c("chr", "start", "end"), sep = "[:-]", remove = FALSE) %>%
  mutate(name = consensus_peak_id, score = 0, strand = ".") %>%
  select(chr, start, end, name, score, strand)

# This object will be used to join binding site info to the annotation table from HOMER
significant_hits <- significant_peak_events %>%
  mutate(is_bind_site = TRUE) %>%
  pivot_wider(names_from = group_name, values_from = is_bind_site, names_glue = "{group_name}_bind_site")





# ==============================================================================
# PART 8: PEAK ANNOTATION WITH HOMER
# ==============================================================================

output_dir_homer <- file.path(output_dir, "HOMER")
dir.create(output_dir_homer, showWarnings = FALSE, recursive = TRUE)
annotation_output_dir <- file.path(output_dir_homer, "PeakAnnotation")
dir.create(annotation_output_dir, showWarnings = FALSE, recursive = TRUE)

# --- 8a. Export BED file of all significant peaks ---
all_peaks_bed_path <- file.path(annotation_output_dir, "all_significant_peaks_for_annotation.bed")
write.table(all_significant_peaks_df, all_peaks_bed_path, sep = "\t", quote = F, row.names = F, col.names = F)
cat("Exported", nrow(all_significant_peaks_df), "peaks for annotation to:", all_peaks_bed_path, "\n")

# --- 8b. Generate annotatePeaks.pl Command ---
annotation_output_file <- file.path(annotation_output_dir, "all_significant_peaks_annotated.txt")
cmd_annotate <- paste(
  "annotatePeaks.pl", shQuote(all_peaks_bed_path), "hg38",
  ">", shQuote(annotation_output_file)
)
# Run this in the the termanl
cat(cmd_annotate) 



# --- 8c. Process and Plot HOMER Annotation Results ---
homer_annotation_file <- annotation_output_file # Use dynamic path

    homer_annots <- fread(homer_annotation_file)
    names(homer_annots)[1] <- "PeakID"
    colnames(homer_annots) <- tolower(gsub(" ", "_", colnames(homer_annots)))
    homer_annots <- homer_annots %>%
        rename(consensus_peak_id = peakid) %>%
        select(-any_of(c("strand", "peak_score", "focus_ratio/region_size", "nearest_unigene")))

    # Join with binding site info
    homer_annots <- homer_annots %>%
        left_join(significant_hits, by = "consensus_peak_id") %>%
        mutate(across(contains("_bind_site"), ~ if_else(is.na(.), FALSE, .)))

    # Simplify annotations
    homer_annots <- homer_annots %>%
        mutate(annot_type_simple = case_when(
            distance_to_tss <= 1000 & distance_to_tss >= -5000 ~ "-5kb to +1kb of TSS",
            grepl("3' UTR", annotation) ~ "3' UTR",
            grepl("5' UTR", annotation) ~ "5' UTR",
            grepl("exon", annotation) ~ "Exon",
            grepl("intron", annotation) ~ "Intron",
            grepl("Intergenic", annotation) ~ "Intergenic/non-coding",
            grepl("non-coding", annotation) ~ "Intergenic/non-coding",
            grepl("TTS", annotation) ~ "TTS",
            TRUE ~ "Other"
        ))
    
   
    # --- Plot Genomic Distribution ---
    plot_data_homer_dist <- homer_annots %>%
      select(consensus_peak_id, annot_type_simple, ends_with("_bind_site")) %>%
      pivot_longer(cols = ends_with("_bind_site"), names_to = "group", values_to = "is_bound") %>%
      filter(is_bound) %>%
      mutate(group = str_remove(group, "_bind_site")) %>%
      distinct(consensus_peak_id, group, .keep_all = TRUE) %>%
      count(group, annot_type_simple, name = "count") %>%
      group_by(group) %>%
      mutate(proportion = count / sum(count)) %>%
      ungroup()

    # Set order for plot
    label_order <- c("-5kb to +1kb of TSS", "5' UTR", "Exon", "Intron", "3' UTR", "TTS", "Intergenic/non-coding")
    plot_data_homer_dist$annot_type_simple <- factor(plot_data_homer_dist$annot_type_simple, levels = rev(label_order))
    plot_data_homer_dist$group <- factor(plot_data_homer_dist$group, levels = rev(intersect(experimental_groups, unique(plot_data_homer_dist$group))))

    # Set colors
    plot_colors <- RColorBrewer::brewer.pal(n = length(unique(plot_data_homer_dist$annot_type_simple)), name = "Accent")
    
    # Make plot
    genomic_distribution_plot <- ggplot(plot_data_homer_dist, aes(x = proportion, y = group, fill = annot_type_simple)) +
      geom_bar(stat = "identity", position = "stack", color = "black", width = 0.7) +
      scale_fill_manual(values = plot_colors, name = "Genomic Region", drop = FALSE) +
      labs(x = "Proportion of Peaks", y = NULL, title = "Genomic Distribution of Significant Peaks") +
      theme_minimal(base_size = 10) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = "black")
      ) +
      guides(fill = guide_legend(reverse = TRUE))
    genomic_distribution_plot
    ggsave(file.path(output_dir, "HOMER_genomic_annotation_proportions_plot.png"), genomic_distribution_plot, width = 8, height = 5, dpi = 300, bg = "white")
    cat("HOMER genomic distribution plot saved.\n")

    # --- Save Promoter-Bound Gene Lists from HOMER ---
    if ("gene_name" %in% names(homer_annots)) {
      promoter_bound_genes_homer <- homer_annots %>%
        filter(annot_type_simple == "-5kb to +1kb of TSS", !is.na(gene_name), trimws(gene_name) != "") %>%
        select(gene_name, ends_with("_bind_site")) %>%
        pivot_longer(cols = ends_with("_bind_site"), names_to = "group", values_to = "is_bound") %>%
        filter(is_bound) %>%
        mutate(group = str_remove(group, "_bind_site")) %>%
        distinct(group, gene_name) %>%
        group_by(group) %>%
        summarise(genes = list(sort(unique(gene_name))), .groups = "drop") %>%
        { setNames(.$genes, .$group) }

      if (length(promoter_bound_genes_homer) > 0) {
        max_len <- max(lengths(promoter_bound_genes_homer))
        promoter_genes_df_homer <- as.data.frame(lapply(promoter_bound_genes_homer, `length<-`, max_len))
        output_file_promoter_genes_homer <- file.path(output_dir, "HOMER_promoter_bound_genes_by_Group.csv")
        fwrite(promoter_genes_df_homer, output_file_promoter_genes_homer, row.names = FALSE, na = "")
        cat("HOMER promoter-bound gene lists saved to:", output_file_promoter_genes_homer, "\n")
      }
    }

    
    # Save the home annotation file
    write.csv(homer_annots,file.path(output_dir,
                                     paste0("HOMER_Annotations_Peak_Set_norm_count_pkb_",norm_count_pkb_cutoff,
                                            "_Log2FC_",log2FoldChange_cutoff,".csv")
                                     )
    )




# ==============================================================================
# PART 9: MOTIF ANALYSIS WITH HOMER. This takes a long time to run (1+ hrs)
# ==============================================================================
motif_analysis_dir <- file.path(output_dir_homer, "MotifAnalysis")
dir.create(motif_analysis_dir, showWarnings = FALSE, recursive = TRUE)


# Loop through each group to export its specific peaks and generate commands
groups_to_export <- unique(significant_peak_events$group_name)
for (current_group in groups_to_export) {
  # --- Filter to the peaks of that group ---
  group_peaks_df <- significant_peak_events %>%
    filter(group_name == current_group) %>%
    distinct(consensus_peak_id) %>%
    separate(consensus_peak_id, into = c("chr", "start", "end"), sep = "[:-]", remove = FALSE) %>%
    mutate(name = consensus_peak_id, score = 0, strand = ".") %>%
    select(chr, start, end, name, score, strand)

  # --- Save those peaks as a .bed in the motif_analysis_dir --- #
  bed_path <- file.path(motif_analysis_dir, paste0(current_group, "_significant_peaks.bed"))
  write.table(group_peaks_df, bed_path, sep = "\t", quote = F, row.names = F, col.names = F)
  cat("Exported", nrow(group_peaks_df), "peaks for Group '", current_group, "' for motif analysis.\n")
}

# --- Run the run_homer_motifs_script.sh script ---
# General usage is:       bash run_homer_motifs.sh <BED_DIR> <GENOME> <NUM WORKERS>
# - Recommended to use one worker per TF if possible.
# Each worker takes one BED file and runs the entire findMotifsGenome.pl command on it from start to finish.

# Specific usage:         bash <path_where_script_is>/run_homer_motifs.sh <motif_analysis_dir> hg38 <one worker per tf>

# bash /Users/rileymullins/Documents/test_scripts_for_multiplex_cc_analysis/TESTING_06302025/testing_07032025/run_homer_motif_script.sh /Users/rileymullins/Documents/test_scripts_for_multiplex_cc_analysis/TESTING_06302025/testing_07032025/full_dataset/analysis_files/HOMER/MotifAnalysis hg38


# --- 9c. Process and Plot HOMER Motif Results ---
motif_result_files <- list.files(motif_analysis_dir, pattern = "knownResults.txt", recursive = TRUE, full.names = TRUE)

if (length(motif_result_files) > 0) {
  motif_results <- map_dfr(motif_result_files, function(file) {
    group_name <- basename(dirname(dirname(file))) # Assumes structure is .../MotifAnalysis/<group_name>/MotifOutput/
    df <- tryCatch(read_tsv(file, comment = "#", show_col_types = FALSE), error = function(e) NULL)
    if (!is.null(df)) {
      df$group_name <- group_name
    }
    df
  })

  if (nrow(motif_results) > 0) {
    motif_results_clean <- motif_results %>%
      rename(Motif_Name = `Motif Name`, P_value = `P-value`, Log_P_value = `Log P-value`) %>%
      mutate(neglogP = Log_P_value * -1) %>%
      mutate(
        Motif_Name_Clean = Motif_Name %>%
          str_remove_all(regex("ChIP-Seq|\\(ChIP-Seq\\)", ignore_case = TRUE)) %>%
          str_remove_all(regex("\\bHOMER\\b", ignore_case = TRUE)) %>%
          str_replace_all("[-/]+", " ") %>%
          str_squish()
      )

    top_n <- 10 # Top n motifs to plot
    motif_plot_data <- motif_results_clean %>%
      group_by(group_name) %>%
      arrange(-neglogP) %>%
      slice_head(n = top_n) %>%
      ungroup() %>%
      mutate(Motif_Name_Clean = fct_reorder(Motif_Name_Clean, neglogP))

    motif_plot <- ggplot(motif_plot_data, aes(x = neglogP, y = Motif_Name_Clean, fill = group_name)) +
      geom_col(show.legend = FALSE, width = 0.7, color = "black", linewidth = 0.25) +
      facet_wrap(~group_name, scales = "free_y", ncol = 2) +
      labs(x = expression(-log[10](P) ~ "enrichment"), y = NULL, title = "Top Enriched Motifs per Group") +
      theme_bw(base_size = 10) +
      theme(
        strip.text = element_text(face = "bold", size = 10, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12, color = "black"),
        panel.spacing = unit(1, "lines"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.3)
      ) +
      scale_fill_brewer(palette = "Set1")

    motif_plot_file <- file.path(output_dir, "homer_top_motifs_plot.png")
    ggsave(motif_plot_file, motif_plot, width = 8, height = 10, dpi = 300, bg = "white")
    cat("HOMER motif plot saved to:", motif_plot_file, "\n")
  }
} else {
  cat("No 'knownResults.txt' files found. Skipping motif plotting.\n")
}
