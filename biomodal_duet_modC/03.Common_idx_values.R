# ==============================================================================
# Extract Common Sites Across All Samples
# ==============================================================================
# Unite all 10 samples and filter to sites present in ALL samples
# ==============================================================================

library(methylKit)
library(dplyr)
library(pheatmap)

# Load methylKit objects
# load("methylkit_objects.RData")

cat("=== Extracting Sites Across All Samples ===\n\n")

# ------------------------------------------------------------------------------
# 1. Check sample information
# ------------------------------------------------------------------------------

cat("5mC samples:\n")
print(mc_methylkit)

cat("\n5hmC samples:\n")
print(hmc_methylkit)

# ------------------------------------------------------------------------------
# 2. Unite all samples - 5mC
# ------------------------------------------------------------------------------

cat("\n[1] Uniting all 5mC samples\n")

# First unite with default settings (min 1 per group)
mc_all_united_raw <- unite(
  mc_methylkit,
  destrand = FALSE,
  min.per.group = NULL
)

cat(sprintf("  Total samples: %d\n", length(mc_methylkit)))
cat(sprintf("  CpGs (≥1 sample per group): %d\n", nrow(mc_all_united_raw)))

# Filter to sites present in ALL samples (no NAs in coverage)
mc_data_raw <- getData(mc_all_united_raw)
keep_rows_mc <- complete.cases(mc_data_raw[, mc_all_united_raw@coverage.index])
mc_all_united <- mc_all_united_raw[keep_rows_mc, ]

cat(sprintf("  CpGs (present in ALL samples): %d\n", nrow(mc_all_united)))

# ------------------------------------------------------------------------------
# 3. Unite all samples - 5hmC
# ------------------------------------------------------------------------------

cat("\n[2] Uniting all 5hmC samples\n")

hmc_all_united_raw <- unite(
  hmc_methylkit,
  destrand = FALSE,
  min.per.group = NULL
)

cat(sprintf("  Total samples: %d\n", length(hmc_methylkit)))
cat(sprintf("  CpGs (≥1 sample per group): %d\n", nrow(hmc_all_united_raw)))

# Filter to sites present in ALL samples
hmc_data_raw <- getData(hmc_all_united_raw)
keep_rows_hmc <- complete.cases(hmc_data_raw[, hmc_all_united_raw@coverage.index])
hmc_all_united <- hmc_all_united_raw[keep_rows_hmc, ]

cat(sprintf("  CpGs (present in ALL samples): %d\n", nrow(hmc_all_united)))

# ------------------------------------------------------------------------------
# 4. Extract methylation values for each sample
# ------------------------------------------------------------------------------

cat("\n[3] Extracting methylation values per sample\n")

# 5mC - Get data
mc_data <- getData(mc_all_united)

# Extract numCs and coverage columns for each sample
n_samples_mc <- length(mc_all_united@sample.ids)

# Calculate beta values for each sample
mc_beta_list <- list()
mc_coverage_list <- list()

for (i in 1:n_samples_mc) {
  sample_name <- mc_all_united@sample.ids[i]
  
  # numCs column index for this sample
  numCs_col <- mc_all_united@numCs.index[i]
  # coverage column index for this sample
  coverage_col <- mc_all_united@coverage.index[i]
  
  # Calculate beta value (%)
  mc_beta_list[[sample_name]] <- (mc_data[, numCs_col] / mc_data[, coverage_col]) * 100
  
  # Get coverage
  mc_coverage_list[[sample_name]] <- mc_data[, coverage_col]
}

# Create beta value matrix
mc_matrix <- data.frame(
  chr = mc_data$chr,
  start = mc_data$start,
  end = mc_data$end,
  strand = mc_data$strand,
  do.call(cbind, mc_beta_list)
)

cat(sprintf("  5mC matrix: %d sites × %d samples\n", 
            nrow(mc_matrix), n_samples_mc))

# Create coverage matrix
mc_coverage_matrix <- data.frame(
  chr = mc_data$chr,
  start = mc_data$start,
  end = mc_data$end,
  strand = mc_data$strand,
  do.call(cbind, mc_coverage_list)
)
colnames(mc_coverage_matrix)[5:ncol(mc_coverage_matrix)] <- 
  paste0(mc_all_united@sample.ids, "_coverage")

# 5hmC - Get data
hmc_data <- getData(hmc_all_united)

n_samples_hmc <- length(hmc_all_united@sample.ids)

hmc_beta_list <- list()
hmc_coverage_list <- list()

for (i in 1:n_samples_hmc) {
  sample_name <- hmc_all_united@sample.ids[i]
  
  numCs_col <- hmc_all_united@numCs.index[i]
  coverage_col <- hmc_all_united@coverage.index[i]
  
  hmc_beta_list[[sample_name]] <- (hmc_data[, numCs_col] / hmc_data[, coverage_col]) * 100
  hmc_coverage_list[[sample_name]] <- hmc_data[, coverage_col]
}

hmc_matrix <- data.frame(
  chr = hmc_data$chr,
  start = hmc_data$start,
  end = hmc_data$end,
  strand = hmc_data$strand,
  do.call(cbind, hmc_beta_list)
)

cat(sprintf("  5hmC matrix: %d sites × %d samples\n", 
            nrow(hmc_matrix), n_samples_hmc))

hmc_coverage_matrix <- data.frame(
  chr = hmc_data$chr,
  start = hmc_data$start,
  end = hmc_data$end,
  strand = hmc_data$strand,
  do.call(cbind, hmc_coverage_list)
)
colnames(hmc_coverage_matrix)[5:ncol(hmc_coverage_matrix)] <- 
  paste0(hmc_all_united@sample.ids, "_coverage")

# ------------------------------------------------------------------------------
# 5. Save matrices
# ------------------------------------------------------------------------------

cat("\n[4] Saving matrices\n")

output_dir <- "Common_Sites_Analysis"
dir.create(output_dir, showWarnings = FALSE)

# Save as CSV
write.csv(mc_matrix, 
          file.path(output_dir, "5mC_common_sites_beta_values.csv"),
          row.names = FALSE)

write.csv(hmc_matrix,
          file.path(output_dir, "5hmC_common_sites_beta_values.csv"),
          row.names = FALSE)

write.csv(mc_coverage_matrix,
          file.path(output_dir, "5mC_common_sites_coverage.csv"),
          row.names = FALSE)

write.csv(hmc_coverage_matrix,
          file.path(output_dir, "5hmC_common_sites_coverage.csv"),
          row.names = FALSE)

# Save as RData
save(mc_all_united, hmc_all_united, mc_matrix, hmc_matrix,
     mc_coverage_matrix, hmc_coverage_matrix,
     file = file.path(output_dir, "Common_Sites_All_Samples.RData"))

cat(sprintf("  Files saved in: %s/\n", output_dir))

# ------------------------------------------------------------------------------
# 6. Summary statistics per sample
# ------------------------------------------------------------------------------

cat("\n[5] Calculating summary statistics\n")

# 5mC sample statistics
mc_beta_df <- mc_matrix[, 5:ncol(mc_matrix)]
mc_coverage_df <- mc_coverage_matrix[, 5:ncol(mc_coverage_matrix)]

mc_sample_stats <- data.frame(
  Sample = mc_all_united@sample.ids,
  Treatment = mc_all_united@treatment,
  Mean_Methylation = colMeans(mc_beta_df, na.rm = TRUE),
  Median_Methylation = apply(mc_beta_df, 2, median, na.rm = TRUE),
  SD_Methylation = apply(mc_beta_df, 2, sd, na.rm = TRUE),
  Mean_Coverage = colMeans(mc_coverage_df, na.rm = TRUE),
  Median_Coverage = apply(mc_coverage_df, 2, median, na.rm = TRUE)
)

# 5hmC sample statistics
hmc_beta_df <- hmc_matrix[, 5:ncol(hmc_matrix)]
hmc_coverage_df <- hmc_coverage_matrix[, 5:ncol(hmc_coverage_matrix)]

hmc_sample_stats <- data.frame(
  Sample = hmc_all_united@sample.ids,
  Treatment = hmc_all_united@treatment,
  Mean_Methylation = colMeans(hmc_beta_df, na.rm = TRUE),
  Median_Methylation = apply(hmc_beta_df, 2, median, na.rm = TRUE),
  SD_Methylation = apply(hmc_beta_df, 2, sd, na.rm = TRUE),
  Mean_Coverage = colMeans(hmc_coverage_df, na.rm = TRUE),
  Median_Coverage = apply(hmc_coverage_df, 2, median, na.rm = TRUE)
)

write.csv(mc_sample_stats,
          file.path(output_dir, "5mC_sample_statistics.csv"),
          row.names = FALSE)

write.csv(hmc_sample_stats,
          file.path(output_dir, "5hmC_sample_statistics.csv"),
          row.names = FALSE)

cat("\n5mC Sample Statistics:\n")
print(mc_sample_stats)

cat("\n5hmC Sample Statistics:\n")
print(hmc_sample_stats)

# ------------------------------------------------------------------------------
# 7. Optional: More flexible unite with min.per.group
# ------------------------------------------------------------------------------

cat("\n[6] Alternative: Unite with relaxed criteria\n")

mc_relaxed_united <- unite(
  mc_methylkit,
  destrand = FALSE,
  min.per.group = 2L
)

hmc_relaxed_united <- unite(
  hmc_methylkit,
  destrand = FALSE,
  min.per.group = 2L
)

cat(sprintf("  5mC (min 2 per group, may have NAs): %d CpGs\n", nrow(mc_relaxed_united)))
cat(sprintf("  5hmC (min 2 per group, may have NAs): %d CpGs\n", nrow(hmc_relaxed_united)))

# ------------------------------------------------------------------------------
# 8. Create sample correlation heatmap
# ------------------------------------------------------------------------------

cat("\n[7] Creating sample correlation plots\n")

# 5mC correlation
mc_cor <- cor(mc_beta_df, use = "pairwise.complete.obs")

pdf(file.path(output_dir, "5mC_sample_correlation_heatmap.pdf"), 
    width = 10, height = 10)
pheatmap(mc_cor,
         main = "5mC Sample Correlation (Common Sites)",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         display_numbers = TRUE,
         number_format = "%.2f",
         fontsize_number = 8)
dev.off()

# 5hmC correlation
hmc_cor <- cor(hmc_beta_df, use = "pairwise.complete.obs")

pdf(file.path(output_dir, "5hmC_sample_correlation_heatmap.pdf"), 
    width = 10, height = 10)
pheatmap(hmc_cor,
         main = "5hmC Sample Correlation (Common Sites)",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         display_numbers = TRUE,
         number_format = "%.2f",
         fontsize_number = 8)
dev.off()

cat("  Correlation heatmaps saved\n")

# ------------------------------------------------------------------------------
# 9. Summary
# ------------------------------------------------------------------------------

cat("\n", rep("=", 80), "\n", sep = "")
cat("SUMMARY\n")
cat(rep("=", 80), "\n", sep = "")

summary_table <- data.frame(
  Metric = c(
    "Total samples",
    "CpGs present in ALL samples",
    "CpGs (min 2 per group, may have NAs)",
    "Output files per modification"
  ),
  `5mC` = c(
    length(mc_methylkit),
    nrow(mc_all_united),
    nrow(mc_relaxed_united),
    "4 files"
  ),
  `5hmC` = c(
    length(hmc_methylkit),
    nrow(hmc_all_united),
    nrow(hmc_relaxed_united),
    "4 files"
  ),
  check.names = FALSE
)

print(summary_table)

cat("\nOutput files (9 total):\n")
cat("  1. 5mC_common_sites_beta_values.csv         - Methylation % for all samples\n")
cat("  2. 5mC_common_sites_coverage.csv            - Coverage for all samples\n")
cat("  3. 5mC_sample_statistics.csv                - Per-sample summary statistics\n")
cat("  4. 5mC_sample_correlation_heatmap.pdf       - Sample correlation heatmap\n")
cat("  5. 5hmC_common_sites_beta_values.csv        - Hydroxymethylation % for all samples\n")
cat("  6. 5hmC_common_sites_coverage.csv           - Coverage for all samples\n")
cat("  7. 5hmC_sample_statistics.csv               - Per-sample summary statistics\n")
cat("  8. 5hmC_sample_correlation_heatmap.pdf      - Sample correlation heatmap\n")
cat("  9. Common_Sites_All_Samples.RData           - R objects for further analysis\n")

cat("\nNote:\n")
cat("  - 'CpGs present in ALL samples' = sites with no missing data across all 10 samples\n")
cat("  - 'min 2 per group' = sites present in ≥2 samples per treatment group (may have NAs)\n")

cat("\n", rep("=", 80), "\n", sep = "")
cat("COMPLETE\n")
cat(rep("=", 80), "\n\n", sep = "")
