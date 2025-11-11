# step_4_DESeq2_Diff_Peaks_HOMER-Annotations_and_Motifs.R


# ==============================================================================
# PART 0: CONFIGURATION
# ==============================================================================

# ---- 0a. Load Libraries ----
suppressPackageStartupMessages({
  library(data.table)
  library(purrr)
  library(stringr)
  library(tibble)
  library(DESeq2)
  library(tidyverse)
  library(parallel)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(circlize)
  library(GenomicRanges)
  library(viridis)
  library(janitor)
})

# ---- 0b. Define Input & Output Directories ----
python_output_dir          <- <dir of steps 1-3 outputs>
raw_insertion_bedgraph_dir <- file.path(python_output_dir, "<STEP 1 output dir>/raw_unique_insertion_count_per_group_bedgraph")
per_group_peaks_dir        <- file.path(python_output_dir, ,"<STEP 3 output dir>/per_group_peak_matrices")
output_dir                 <- file.path(python_output_dir, "STEP_4_output_deseq2_peak_analysis")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 0c. Calculate Total Insertions per Group ----
bedgraph_files <- list.files(
  raw_insertion_bedgraph_dir,
  pattern = "\\.bedgraph$",
  full.names = TRUE
)

# Initialize result named numeric vector
total_counts <- numeric()

# Loop over files
for (file in bedgraph_files) {
  # Read only the 4th column (V4) using fread
  dt <- fread(file, select = 4, col.names = "V4")
  
  # Sum the V4 values
  total <- sum(dt$V4)
  
  # Strip .bedgraph extension from filename
  name <- sub("\\.bedgraph$", "", basename(file))
  
  # Store the result
  total_counts[name] <- total
}

cat("Total insertion counts per group:\n")
print(total_counts)

# ---- 0d. Define Experimental Groups ----
deseq_output_dir <- file.path(output_dir, "DESeq2_results_by_group")
dir.create(deseq_output_dir, recursive = TRUE, showWarnings = FALSE)

all_groups          <- names(total_counts)
control_group       <- "HyPBase"
experimental_groups <- setdiff(all_groups, control_group)

cat("Control group:", control_group, "\n")
cat("Experimental groups:", paste(experimental_groups, collapse = ", "), "\n\n")

# ==============================================================================
# PART 1: PRE-PROCESSING & NORMALIZATION
# ==============================================================================

# ---- 1a. Manual Size Factor Calculation ----
log_geo_means       <- log(total_counts[total_counts > 0])
geo_mean            <- exp(mean(log_geo_means))
manual_size_factors <- total_counts / geo_mean
names(manual_size_factors) <- names(total_counts)

print(manual_size_factors)

# ---- 1b. Generate Normalized Peak Matrices For Each Group ----
normalized_peaks_dir <- file.path(output_dir, "normalized_peak_matrices")
dir.create(normalized_peaks_dir, recursive = TRUE, showWarnings = FALSE)

# Loop over each of the peak matrices
for (group in experimental_groups) {
  # Define the path for the input raw count matrix
  matrix_file <- file.path(per_group_peaks_dir, paste0(group, "_peak_matrix.tsv"))
  if (!file.exists(matrix_file)) {
    warning(paste("Matrix file not found for", group, "- skipping normalization."))
    next
  }
  
  # Read the raw count matrix
  group_peak_matrix <- fread(matrix_file)
  
  # Make a copy to store the normalized counts
  normalized_matrix <- copy(group_peak_matrix)
  
  # Find all columns that contain raw insertion counts
  count_cols <- grep("_insertion_count$", colnames(normalized_matrix), value = TRUE)
  
  # Loop through each count column to normalize it
  for (col_name in count_cols) {
    # 1. Get the group name from the column name
    group_name <- sub("_insertion_count$", "", col_name)
    
    # 2. Find the matching size factor for that group
    size_factor <- manual_size_factors[group_name]
    
    # 3. Divide the entire column by its size factor to get normalized counts.
    normalized_matrix[, (col_name) := .SD[[col_name]] / size_factor]
  }
  
  # Rename the columns to indicate they are now normalized
  new_col_names <- sub("_insertion_count$", "_normalized_count", count_cols)
  setnames(normalized_matrix, old = count_cols, new = new_col_names)
  
  # Define the output file path and save the normalized matrix
  output_file <- file.path(normalized_peaks_dir, paste0(group, "_normalized_peak_matrix.csv"))
  fwrite(normalized_matrix, output_file)
  
  cat("Saved normalized matrix for:", group, "\n")
}
cat("--- Normalization complete. Files saved in:", normalized_peaks_dir, "---\n\n")

# ---- 1c. Load and Combine All Normalized Matrices ----
norm_matrix_files <- list.files(
  normalized_peaks_dir,
  pattern = "_normalized_peak_matrix.csv$",
  full.names = TRUE
)

# Read all files into a list and combine
all_norm_matrices_list <- lapply(norm_matrix_files, fread)
combined_norm_matrix   <- rbindlist(all_norm_matrices_list)

# Get names of normalized count columns
norm_count_cols <- grep("_normalized_count$", colnames(combined_norm_matrix), value = TRUE)

# Create a final, tidy, long-format table of all normalized counts
final_norm_matrix_long <- combined_norm_matrix %>%
  select(peak_id, source_group = group_name, all_of(norm_count_cols)) %>%
  pivot_longer(
    cols = -c(peak_id, source_group),
    names_to = "group",
    values_to = "normalized_count"
  ) %>%
  mutate(group = sub("_normalized_count$", "", group))

# Save the combined long-format matrix
output_file <- file.path(normalized_peaks_dir, "combined_long_format_normalized_peak_matrix.csv")
fwrite(final_norm_matrix_long, output_file)

# ==============================================================================
# PART 2: WRAPPER FUNCTION FOR DESEQ2 ANALYSIS
# ==============================================================================

run_all_comparisons_for_group <- function(current_group, all_groups, control_group, per_group_peaks_dir, output_dir, manual_size_factors) {
  
  # Create a dedicated subdirectory for the current group's results
  group_output_dir <- file.path(output_dir, "DESeq2_results_by_group", current_group)
  dir.create(group_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  matrix_file <- file.path(per_group_peaks_dir, paste0(current_group, "_peak_matrix.tsv"))
  if (!file.exists(matrix_file)) {
    warning(paste("Matrix file not found for", current_group, "- skipping."))
    return(NULL)
  }
  
  group_peak_matrix <- fread(matrix_file)
  
  disp_dir <- file.path(output_dir, "dispersion_plots")
  dir.create(disp_dir, showWarnings = FALSE, recursive = TRUE)
  
  # --- Pairwise: group vs. every other group ---
  for (other_group in setdiff(all_groups, current_group)) {
    comparison_name <- paste0(current_group, "_vs_", other_group)
    
    group_col_name      <- paste0(current_group, "_insertion_count")
    comparator_col_name <- paste0(other_group, "_insertion_count")
    
    count_matrix <- group_peak_matrix %>%
      select(all_of(c(group_col_name, comparator_col_name))) %>%
      as.data.frame()
    rownames(count_matrix) <- group_peak_matrix$peak_id
    colnames(count_matrix) <- c(current_group, other_group)
    
    keep_rows       <- (rowSums(count_matrix[, current_group, drop = FALSE], na.rm = TRUE) > 0) |
                       (rowSums(count_matrix[, other_group, drop = FALSE], na.rm = TRUE) > 0)
    pairwise_matrix <- count_matrix[keep_rows, , drop = FALSE]
    
    if (nrow(pairwise_matrix) == 0) next
    
    col_data <- data.frame(
      condition = factor(c("treatment", "control")),
      row.names = c(current_group, other_group)
    )
    
    dds <- DESeqDataSetFromMatrix(
      countData = pairwise_matrix,
      colData = col_data,
      design = ~condition
    )
    
    # Apply pre-calculated manual size factors for 1-vs-1 comparisons
    sizeFactors(dds) <- manual_size_factors[c(current_group, other_group)]
    
    # Estimate dispersions
    dds_for_disp           <- dds
    design(dds_for_disp)   <- formula(~1)
    dds_for_disp           <- estimateDispersions(dds_for_disp, fitType = "local")
    dispersions(dds)       <- dispersions(dds_for_disp)
    
    dds <- nbinomWaldTest(dds)
    
    res <- results(dds, contrast = c("condition", "treatment", "control"), altHypothesis = "greater")
    
    res_df <- as.data.frame(res) %>%
      rownames_to_column("peak_id") %>%
      mutate(comparison = comparison_name)
    
    # Get normalized counts for this specific comparison
    norm_counts_df <- as.data.frame(counts(dds, normalized = TRUE)) %>%
      rownames_to_column("peak_id")
    
    # Rename normalized count columns for clarity
    old_grp_col <- current_group
    old_cmp_col <- other_group
    colnames(norm_counts_df)[colnames(norm_counts_df) == old_grp_col] <- "norm_group_counts"
    colnames(norm_counts_df)[colnames(norm_counts_df) == old_cmp_col] <- "norm_comparator_counts"
    
    # Prepare metadata to join
    metadata_to_join <- group_peak_matrix %>%
      select(peak_id, peak_width, all_of(c(group_col_name, comparator_col_name)))
    
    # Join all data for a comprehensive results table
    joined <- res_df %>%
      left_join(metadata_to_join, by = "peak_id") %>%
      left_join(norm_counts_df,   by = "peak_id")
    
    # Rename raw-count columns for clarity
    old_grp_col <- group_col_name
    old_cmp_col <- comparator_col_name
    colnames(joined)[colnames(joined) == old_grp_col] <- "group_counts"
    colnames(joined)[colnames(joined) == old_cmp_col] <- "comparator_counts"
    
    # Write results to file
    fwrite(joined, file.path(group_output_dir, paste0(comparison_name, "_DESeq2_results.csv")))
    
    # Save dispersion plot for comparisons against the main control
    if (other_group == control_group) {
      png(
        file.path(disp_dir, paste0(comparison_name, "_dispersion_plot.png")),
        width = 6, height = 6, units = "in", res = 300
      )
      plotDispEsts(dds_for_disp, main = comparison_name)
      dev.off()
    }
  }
  
  return(TRUE)
}

# ==============================================================================
# PART 3: EXECUTE PARALLEL ANALYSIS
# ==============================================================================
cat("Starting differential analyses in parallel\n")
num_cores <- max(1, parallel::detectCores() - 1)

# 1. Create a socket cluster
cl <- parallel::makeCluster(num_cores)

# 2. Export global variables and load libraries on each worker
parallel::clusterExport(cl, varlist = c(
  "run_all_comparisons_for_group",
  "all_groups", "control_group",
  "per_group_peaks_dir", "output_dir", "manual_size_factors",
  "total_counts"
))
parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(DESeq2)
    library(tidyverse)
  })
})

# 3. Run the analysis in parallel
parallel::parLapply(
  cl,
  experimental_groups,
  run_all_comparisons_for_group,
  all_groups = all_groups,
  control_group = control_group,
  per_group_peaks_dir = per_group_peaks_dir,
  output_dir = output_dir,
  manual_size_factors = manual_size_factors
)

# 4. Shut down the cluster
parallel::stopCluster(cl)

# ==============================================================================
# PART 4: COMBINE AND SAVE ALL DESEQ2 RESULTS
# ==============================================================================

# Find all result files
all_files <- list.files(
  path = deseq_output_dir,
  pattern = "_DESeq2_results.csv$",
  recursive = TRUE,
  full.names = TRUE
)

# Read all files into a single data frame
all_results <- map_dfr(all_files, fread)

# Separate the results into 'vs. Control' and 'vs. Group'
combined_control_results <- all_results %>%
  filter(grepl(paste0("_vs_", control_group, "$"), comparison))

combined_pairwise_results <- all_results %>%
  filter(!grepl(paste0("_vs_", control_group, "$"), comparison))

# Add a 'group' column for easier filtering
combined_control_results <- combined_control_results %>%
  mutate(group = str_remove(comparison, "_vs_.*"))
combined_pairwise_results <- combined_pairwise_results %>%
  mutate(group = str_remove(comparison, "_vs_.*"))

# Print summary and save combined files
cat("Combined", nrow(combined_control_results), "Group vs. Control results.\n")
cat("Combined", nrow(combined_pairwise_results), "Group vs. Group results.\n")
cat("Combined", nrow(all_results), "total results.\n")

fwrite(combined_control_results, file = file.path(output_dir, "COMBINED_all_Group_vs_Control_results.csv"))
fwrite(combined_pairwise_results, file = file.path(output_dir, "COMBINED_all_Group_vs_Group_results.csv"))
fwrite(all_results, file = file.path(output_dir, "COMBINED_all_Group_vs_Group_and_Group_v_Control_results.csv"))

# ==============================================================================
# PART 5: ANALYSIS OF PEAKS SIGNIFICANTLY HIGHER THAN CONTROL
# ==============================================================================

# ---- 5a. Identify and Save Significant Peaks (vs. Control) ----
norm_count_cutoff     <- 10
log2FoldChange_cutoff <- 1

group_output_dir <- file.path(output_dir, "DESeq2_results_peaks_above_control_by_group")
dir.create(group_output_dir, recursive = TRUE, showWarnings = FALSE)

all_sig <- bind_rows(lapply(experimental_groups, function(group_name) {
  df <- combined_control_results %>%
    filter(
      group == group_name,
      log2FoldChange >= log2FoldChange_cutoff,
      norm_group_counts >= norm_count_cutoff
    )
  
  # Save filtered hits for this group
  df_save <- df %>%
    separate(peak_id, into = c("chr", "start", "end"), sep = ":|-", convert = TRUE)
  
  group_output_file <- file.path(
    output_dir,
    "DESeq2_results_peaks_above_control_by_group",
    paste0(
      group_name, "_vs_control_DESeq2_log2fc_cutoff_",
      log2FoldChange_cutoff, "_norm_count_cutoff_",
      norm_count_cutoff, ".csv"
    )
  )
  
  fwrite(df_save, group_output_file)
  
  return(df)
}))
cat("Found", nrow(all_sig), "significant hits (vs control) to merge.\n")

# ---- 5b. Merge Significant Peaks into Consensus Regions ----
gr_orig <- all_sig %>%
  separate(peak_id, into = c("chr", "start", "end"), sep = ":|-", convert = TRUE) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

merged_gr <- reduce(gr_orig, with.revmap = TRUE)
merged_gr$merged_peak_id <- paste0("merged_peak_", seq_along(merged_gr))
cat("Merged into", length(merged_gr), "regions.\n")

# ---- 5c. Aggregate Normalized Counts for Merged Peaks ----
# 1. Map each merged peak back to its original constituent peaks
map_merged_to_original <- data.frame(
  merged_peak_id = rep(merged_gr$merged_peak_id, lengths(mcols(merged_gr)$revmap)),
  orig_idx       = unlist(mcols(merged_gr)$revmap)
) %>%
  left_join(
    all_sig %>%
      mutate(orig_idx = row_number()) %>%
      select(orig_idx, peak_id, source_group = group),
    by = "orig_idx"
  )

# 2. Join map with comprehensive normalized counts and aggregate
norm_matrix <- map_merged_to_original %>%
  left_join(final_norm_matrix_long, by = c("peak_id", "source_group")) %>%
  group_by(merged_peak_id, group) %>%
  summarise(total_norm = sum(normalized_count, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = group, values_from = total_norm, values_fill = 0)

# 3. Assemble and save the final merged matrix
final_output_df_control <- as.data.frame(merged_gr)[, c("merged_peak_id", "seqnames", "start", "end")] %>%
  dplyr::rename(chrom = seqnames) %>%
  left_join(norm_matrix, by = "merged_peak_id") %>%
  replace(is.na(.), 0)

fwrite(
  final_output_df_control,
  file.path(output_dir, "merged_significant_peaks_vs_HyPBase_matrix.csv")
)
cat(
  "Saved final matrix for peaks above control:",
  file.path(output_dir, "merged_significant_peaks_vs_HyPBase_matrix.csv"),
  "\n\n"
)

# ===================================================================================
# PART 6: FIND GROUP-SPECIFIC PEAKS
# ===================================================================================

# ---- 6a. Find Peaks More Abundant in One Group vs. ALL Others ----
specific_peak_list <- lapply(experimental_groups, function(current_group) {
  other_groups <- setdiff(all_groups, current_group)
  
  # Get a list of significant peak vectors, one for each comparison
  list_of_peak_id_vectors <- lapply(other_groups, function(other_group) {
    comp_name_1 <- paste(current_group, other_group, sep = "_vs_")
    all_results %>%
      filter(
        comparison == comp_name_1,
        log2FoldChange >= log2FoldChange_cutoff,
        norm_group_counts >= norm_count_cutoff
      ) %>%
      pull(peak_id)
  })
  
  # Find the intersection of all peak sets
  specific_peaks <- Reduce(intersect, list_of_peak_id_vectors)
  
  return(data.frame(peak_id = specific_peaks, specific_for_group = current_group))
})

all_specific_peaks <- bind_rows(specific_peak_list)
cat("\nFound a total of", nrow(all_specific_peaks), "group-specific peaks across all groups.\n")

# ---- 6b. Assemble Final Table for Group-Specific Peaks ----
# 1. Join specific peaks with the long-format normalized count table
specific_peaks_with_counts <- all_specific_peaks %>%
  left_join(
    final_norm_matrix_long,
    by = c("peak_id", "specific_for_group" = "source_group")
  )


# 2. Pivot data to create a wide matrix
final_summary_df <- specific_peaks_with_counts %>%
  pivot_wider(
    id_cols = c(peak_id, specific_for_group),
    names_from = group,
    values_from = normalized_count,
    values_fill = 0
  )

# 3. Add coordinates and reorder columns
final_output_df_specific <- final_summary_df %>%
  separate(peak_id, into = c("chrom", "start", "end"), sep = "[:-]", remove = FALSE, convert = TRUE) %>%
  select(chrom, start, end, peak_id, specific_for_group, everything())

# 4. Save the final summary matrix
fwrite(
  final_output_df_specific,
  file.path(output_dir, "group_specific_peaks_summary.csv")
)
cat(
  "Saved final matrix for group-specific peaks:",
  file.path(output_dir, "group_specific_peaks_summary.csv"),
  "\n\n"
)

# ---- 6c. Save Individual Group-Specific Peak Files ----
group_output_dir <- file.path(output_dir, "DESeq2_results_peaks_group_specific")
dir.create(group_output_dir, recursive = TRUE, showWarnings = FALSE)

unique_groups <- unique(final_output_df_specific$specific_for_group)

for (group in unique_groups) {
  group_df <- final_output_df_specific %>%
    filter(specific_for_group == group) %>%
    select(chrom, start, end, everything()) # Ensures chrom/start/end are first
  
  group_output_file <- file.path(
    output_dir,
    "DESeq2_results_peaks_group_specific",
    paste0(
      group, "_DESeq2_log2fc_cutoff_",
      log2FoldChange_cutoff, "_norm_count_cutoff_",
      norm_count_cutoff, ".csv"
    )
  )
  fwrite(group_df, group_output_file)
}


# ==============================================================================
# PART 7: SUMMARY PLOTS
# ==============================================================================

# ---- 7a. Define Group Order & Palette ----
group_order_plots <- c() # Define order
palette_colors <- c() # Define color palette
names(palette_colors) <- group_order_plots

  # This will be used in part 10c too.

# ---- 7b. Bar Plot of Total Raw Insertions ----
bar_data_insertions <- data.frame(
  group_name = names(total_counts),
  total_insertions = total_counts
) %>%
  filter(group_name %in% group_order_plots) %>%
  mutate(group_name = factor(group_name, levels = group_order_plots))

p_insertions <- ggplot(bar_data_insertions, aes(x = group_name, y = total_insertions, fill = group_name)) +
  geom_col(color = "black", linewidth = 0.1, width = 0.75) +
  geom_text(aes(label = total_insertions), vjust = -0.4, size = 3, color = "black") +
  scale_fill_manual(values = palette_colors, guide = "none") +
  labs(x = NULL, y = "Total Raw Insertions") +
  theme_minimal(base_size = 10) +
  theme(
    axis.title.x       = element_blank(),
    axis.title.y       = element_text(color = "black", size = 10),
    axis.text.x        = element_text(angle = 45, hjust = 1, color = "black", size = 10),
    axis.text.y        = element_text(color = "black", size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.line          = element_line(color = "black", linewidth = 0.2)
  )

ggsave(
  file.path(output_dir, 'total_raw_insertions_per_group_barplot.png'),
  p_insertions, width = 5, height = 4.5, dpi = 300, bg = "white"
)
ggsave(
  file.path(output_dir, 'total_raw_insertions_per_group_barplot.pdf'),
  p_insertions, width = 5, height = 4.5, bg = "white"
)

# ---- 7c. Bar Plot of Total Significant Peaks (vs Control) ----
bar_data_simple <- all_sig %>%
  dplyr::count(group, name = "num_peaks") %>%
  mutate(group = factor(group, levels = group_order_plots))

p_simple <- ggplot(bar_data_simple, aes(x = group, y = num_peaks, fill = group)) +
  geom_col(color = "black", linewidth = 0.1, width = 0.75) +
  geom_text(aes(label = num_peaks), vjust = -0.3, size = 3, color = "black") +
  scale_fill_manual(values = palette_colors, guide = "none") +
  labs(x = NULL, y = "Number of Peaks\nabove HyPBase") +
  theme_minimal(base_size = 10) +
  theme(
    axis.title.x       = element_blank(),
    axis.title.y       = element_text(color = "black", size = 10),
    axis.text.x        = element_text(angle = 45, hjust = 1, color = "black", size = 10),
    axis.text.y        = element_text(color = "black", size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.line          = element_line(color = "black", linewidth = 0.2)
  )

ggsave(
  file.path(output_dir, "significant_peaks_per_group_barplot.png"),
  p_simple, width = 5, height = 4.5, dpi = 300, bg = "white"
)
ggsave(
  file.path(output_dir, "significant_peaks_per_group_barplot.pdf"),
  p_simple, width = 5, height = 4.5, bg = "white"
)

# ---- 7d. Bar Plot with Log2FC Bins ----
bar_data_binned <- all_sig %>%
  mutate(
    log2fc_bin = case_when(
      log2FoldChange >= 1 & log2FoldChange < 2 ~ "1-2",
      log2FoldChange >= 2 & log2FoldChange < 3 ~ "2-3",
      log2FoldChange >= 3 & log2FoldChange < 4 ~ "3-4",
      log2FoldChange >= 4 & log2FoldChange < 5 ~ "4-5",
      log2FoldChange >= 5 & log2FoldChange < 6 ~ "5-6",
      log2FoldChange >= 6                    ~ "6+",
      TRUE                                   ~ NA_character_
    )
  ) %>%
  filter(!is.na(log2fc_bin)) %>%
  dplyr::count(group, log2fc_bin, name = "num_peaks") %>%
  mutate(group = factor(group, levels = group_order_plots))

bar_data_binned$log2fc_bin <- factor(
  bar_data_binned$log2fc_bin,
  levels = rev(c("1-2", "2-3", "3-4", "4-5", "5-6", "6+"))
)

# Define bin labels and colors
bin_labels <- c("1-2", "2-3", "3-4", "4-5", "5-6", "6+")
bin_colors <- viridis(6)
names(bin_colors) <- bin_labels

p_binned <- ggplot(bar_data_binned, aes(x = group, y = num_peaks, fill = log2fc_bin)) +
  geom_col(color = "black", linewidth = 0.1, width = 0.75) +
  scale_fill_manual(values = bin_colors, name = "log2FC bin", drop = FALSE) +
  labs(x = NULL, y = "Number of Peaks\nabove HyPBase") +
  theme_minimal(base_size = 10) +
  theme(
    axis.title.x     = element_blank(),
    axis.title.y     = element_text(color = "black", size = 10),
    axis.text.x      = element_text(angle = 45, hjust = 1, color = "black", size = 10),
    axis.text.y      = element_text(color = "black", size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(color = "black", linewidth = 0.2)
  )

ggsave(
  file.path(output_dir, 'peaks_per_group_barplot_log2fc_bins.png'),
  p_binned, width = 5, height = 4.5, dpi = 300, bg = "white"
)
ggsave(
  file.path(output_dir, 'peaks_per_group_barplot_log2fc_bins.pdf'),
  p_binned, width = 5, height = 4.5, bg = "white"
)

# ==============================================================================
# PART 8: PER-GROUP PEAK ANNOTATION & MOTIF ANALYSIS WITH HOMER
# ==============================================================================

# ---- 8a. Per-Group Peak Annotation ----
output_dir_homer      <- file.path(output_dir, "HOMER")
bed_output_dir        <- file.path(output_dir_homer, "Per_Group_BEDs")
annotation_output_dir <- file.path(output_dir_homer, "Per_Group_Annotations")
dir.create(output_dir_homer, showWarnings = FALSE, recursive = TRUE)
dir.create(bed_output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(annotation_output_dir, showWarnings = FALSE, recursive = TRUE)

# Loop through each group to generate BED files and print HOMER commands
# Copy and paste HOMER commands into terminal.
# Make sure HOMER is installed and the HOMER bin directory is in your PATH.
for (current_group in experimental_groups) {
  # 1. Filter for significant peaks
  group_sig_peaks <- all_sig %>%
    filter(group == current_group)
  
  if (nrow(group_sig_peaks) == 0) {
    cat(paste0("No significant peaks for ", current_group, ", skipping.\n"))
    next
  }
  
  # 2. Create a BED-formatted data frame
  group_bed_df <- group_sig_peaks %>%
    distinct(peak_id) %>%
    separate(peak_id, into = c("chr", "start", "end"), sep = "[:-]", remove = FALSE) %>%
    mutate(name = peak_id, score = 0, strand = ".") %>%
    select(chr, start, end, name, score, strand)
  
  # 3. Define paths
  group_bed_path        <- file.path(bed_output_dir, paste0(current_group, "_significant_peaks.bed"))
  group_annotation_file <- file.path(annotation_output_dir, paste0(current_group, "_annotated.txt"))
  
  # 4. Write BED file
  write.table(group_bed_df, group_bed_path, sep = "\t", quote = F, row.names = F, col.names = F)
  
  # 5. Generate and print the HOMER command for annotation
  cmd_annotate <- paste(
    "annotatePeaks.pl", shQuote(group_bed_path), "hg38",
    ">", shQuote(group_annotation_file)
  )
  cat(cmd_annotate, "\n\n")
}




# COPY AND PASTE THE COMMANDS IN THE TERMINAL, RUN, AND THEN CONTINUE ON TO STEP 9B.





# ---- 8b. Process HOMER Annotations and Plot Genomic Distribution ----
# 1. Read and process each group's HOMER annotation file
all_annotation_files <- list.files(
  annotation_output_dir,
  pattern = "_annotated.txt$",
  full.names = TRUE
)
processed_annotations_per_group <- list()

for (file in all_annotation_files) {
  group_name <- sub("_annotated.txt$", "", basename(file))
  annots     <- fread(file) %>%
    mutate(group = group_name)
  
  # Standardize column names and select
  colnames(annots) <- tolower(gsub(" ", "_", colnames(annots)))
  annots <- annots %>%
    select(-any_of(c("strand", "peak_score", "focus_ratio/region_size", "nearest_unigene")))
  
  # Simplify annotation categories
  annots <- annots %>%
    mutate(
      annot_type_simple = case_when(
        distance_to_tss <= 1000 & distance_to_tss >= -5000 ~ "-5kb to +1kb of TSS",
        grepl("3' UTR", annotation)                         ~ "3' UTR",
        grepl("5' UTR", annotation)                         ~ "5' UTR",
        grepl("exon", annotation)                           ~ "Exon",
        grepl("intron", annotation)                         ~ "Intron",
        grepl("Intergenic", annotation)                     ~ "Intergenic/non-coding",
        grepl("non-coding", annotation)                     ~ "Intergenic/non-coding",
        grepl("TTS", annotation)                            ~ "TTS",
        TRUE                                                ~ "Other"
      )
    )
  
  # Clean up peak_id column name
  colnames(annots)[1] <- "peak_id"
  
  processed_annotations_per_group[[group_name]] <- annots
}

# 2. Combine processed annotations and compute proportions for plotting
plot_data_list <- list()
for (group_name in names(processed_annotations_per_group)) {
  annots <- processed_annotations_per_group[[group_name]]
  df     <- annots %>%
    dplyr::count(group, annot_type_simple, name = "count") %>%
    group_by(group) %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup()
  
  plot_data_list[[group_name]] <- df
}
plot_data_homer_dist <- bind_rows(plot_data_list)

# 3. Factor ordering for consistent plotting
label_order <- c(
  "-5kb to +1kb of TSS", "5' UTR", "Exon",
  "Intron", "3' UTR", "TTS", "Intergenic/non-coding"
)
plot_data_homer_dist <- plot_data_homer_dist %>%
  mutate(
    annot_type_simple = factor(annot_type_simple, levels = rev(label_order)),
    group             = factor(group, levels = rev(intersect(experimental_groups, unique(group))))
  )

# 4. Generate palette and create plot
plot_colors <- RColorBrewer::brewer.pal(
  n    = length(levels(plot_data_homer_dist$annot_type_simple)),
  name = "Accent"
)

genomic_distribution_plot <- ggplot(
  plot_data_homer_dist,
  aes(x = group, y = proportion, fill = annot_type_simple)
) +
  geom_col(color = "black", linewidth = 0.2, width = 0.75) +
  scale_fill_manual(
    values = plot_colors,
    name   = "Genomic Region",
    drop   = FALSE,
    guide  = guide_legend(reverse = TRUE)
  ) +
  coord_flip() +
  labs(x = NULL, y = "Proportion of Peaks") +
  theme_minimal(base_size = 10) +
  theme(
    axis.title.x     = element_blank(),
    axis.title.y     = element_text(color = "black", size = 10),
    axis.text.x      = element_text(color = "black", size = 10),
    axis.text.y      = element_text(color = "black", size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line        = element_line(color = "black", size = 0.5),
    plot.title       = element_text(hjust = 0.5, face = "bold")
  )

# 5. Save the plot
ggsave(
  file.path(output_dir, "HOMER_genomic_annotation_proportions_plot_accent.png"),
  genomic_distribution_plot, width = 6.5, height = 4, dpi = 300, bg = "white"
)
ggsave(
  file.path(output_dir, "HOMER_genomic_annotation_proportions_plot_accent.pdf"),
  genomic_distribution_plot, width = 6.5, height = 4, bg = "white"
)

# ---- 8c. Save Promoter-Bound Gene Lists from HOMER ----
promoter_bound_genes_homer <- lapply(processed_annotations_per_group, function(annots) {
  annots %>%
    filter(
      annot_type_simple == "-5kb to +1kb of TSS",
      !is.na(gene_name),
      trimws(gene_name) != "",
      gene_type %in% "protein-coding"
    ) %>%
    pull(gene_name) %>%
    unique() %>%
    sort()
})

# Drop groups with zero genes and pad to a rectangle for CSV export
promoter_bound_genes_homer <- promoter_bound_genes_homer[lengths(promoter_bound_genes_homer) > 0]

if (length(promoter_bound_genes_homer) > 0) {
  max_len <- max(lengths(promoter_bound_genes_homer))
  promoter_genes_df_homer <- as.data.frame(
    lapply(promoter_bound_genes_homer, `length<-`, max_len),
    stringsAsFactors = FALSE
  )
  output_file_promoter_genes_homer <- file.path(
    output_dir,
    "HOMER_promoter_bound_genes_protein_coding_genes_by_Group.csv"
  )
  fwrite(
    promoter_genes_df_homer,
    output_file_promoter_genes_homer,
    row.names = FALSE,
    na = ""
  )
  cat(
    "HOMER promoter-bound gene lists saved to:",
    output_file_promoter_genes_homer, "\n"
  )
}

# ---- 8d. Join Annotations Back to DESeq2 Results ----
# Combine all processed annotation files into one dataframe
combined_annotations <- bind_rows(processed_annotations_per_group)

# (Confirmation checks removed for neatness)

# Join annotations to the significant 'vs. Control' results
all_sig_annots <- left_join(all_sig, combined_annotations, by = c("peak_id", "group")) %>%
  select(-group_counts, -comparator_counts)

# Save the final annotated results table
fwrite(
  all_sig_annots,
  file = file.path(output_dir, "HOMER_annotations_COMBINED_all_Group_vs_Control_results.csv")
)

# ==============================================================================
# PART 9: MOTIF ANALYSIS WITH HOMER
# ==============================================================================

# ---- 9a. Create Directory and Export BED Files for Motif Discovery ----
motif_analysis_dir <- file.path(output_dir_homer, "HOMER_motif_discovery")
dir.create(motif_analysis_dir, showWarnings = FALSE, recursive = TRUE)

groups_to_export <- unique(all_sig_annots$group)

for (current_group in groups_to_export) {
  group_peaks_df <- all_sig_annots %>%
    filter(group == current_group) %>%
    distinct(peak_id, .keep_all = TRUE) %>%
    transmute(chr, start, end, name = peak_id, score = 0, strand = ".")
  
  bed_path <- file.path(motif_analysis_dir, paste0(current_group, "_significant_peaks.bed"))
  write.table(group_peaks_df, bed_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  cat("Exported", nrow(group_peaks_df), "peaks for Group '", current_group, "' to:", bed_path, "\n")
}

# --- 9b. Run the run_homer_motifs_script.sh script ---
# General usage is:       bash run_homer_motifs.sh <BED_DIR> <GENOME>

# Specific usage:         bash <path to script>/run_homer_motifs.sh <motif_analysis_dir> hg38 <num_workers>

# Example: bash run_homer_motifs.sh /path/to/HOMER_motif_analysis hg38 12


# ---- 10c. Process and Plot HOMER Motif Results ----
# Read and process all 'knownResults.txt' files
motif_result_files <- list.files(motif_analysis_dir, pattern = "knownResults.txt", recursive = TRUE, full.names = TRUE)

motif_results <- purrr::map_dfr(motif_result_files, function(file) {
  readr::read_tsv(file, comment = "#", show_col_types = FALSE) %>%
    select(!starts_with("...")) %>%
    mutate(group_name = basename(dirname(file)) %>% str_remove("_motif_output$"))
}) %>%
  janitor::clean_names()

# Prepare top 10 motifs per group for plotting
motif_plot_data <- motif_results %>%
  mutate(
    neglogP          = log_p_value * -1,
    motif_name_clean = motif_name %>%
      str_remove_all(regex("ChIP-Seq|\\(ChIP-Seq\\)", ignore_case = TRUE)) %>%
      str_remove_all(regex("\\bHOMER\\b", ignore_case = TRUE)) %>%
      str_replace_all("[-/]+", " ") %>%
      str_squish()
  ) %>%
  mutate(group_name = factor(group_name, levels = group_order_plots)) %>%
  group_by(group_name) %>%
  arrange(-neglogP) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(motif_name_clean = reorder(motif_name_clean, neglogP))

# Build and save the faceted bar plot
motif_plot <- ggplot(motif_plot_data, aes(x = neglogP, y = motif_name_clean, fill = group_name)) +
  geom_col(show.legend = FALSE, width = 0.7, color = "black", linewidth = 0.25) +
  facet_wrap(~ group_name, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = palette_colors) +
  labs(x = expression(-log[10](P) ~ "enrichment"), y = NULL) +
  theme_bw(base_size = 10) +
  theme(
    strip.text       = element_text(face = "bold", size = 10),
    axis.text        = element_text(color = "black", size = 8),
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 12),
    panel.spacing    = unit(1, "lines"),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.3)
  )

ggsave(
  file.path(output_dir, "homer_top_motifs_plot.png"),
  motif_plot, width = 10, height = 10, dpi = 300, bg = "white"
)
ggsave(
  file.path(output_dir, "homer_top_motifs_plot.pdf"),
  motif_plot, width = 10, height = 10, bg = "white"
)
