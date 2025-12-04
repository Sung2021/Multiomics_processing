# ==============================================================================
# DMR Analysis Pipeline - Treatment2 vs Control
# ==============================================================================
# Comparing: Treatment2 (treatment=1) vs Control (treatment=0)
# Complete analysis: QC, filtering, normalization, DM, DMR, visualization
# ==============================================================================

library(methylKit)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyr)

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

# Load methylKit objects (from previous script)
# load("methylkit_objects.RData")

# Output directory
output_dir <- "DMR_Analysis_Results"
dir.create(output_dir, showWarnings = FALSE)

comp_dir <- file.path(output_dir, "Treatment2_vs_Control")
dir.create(comp_dir, showWarnings = FALSE)

# Analysis parameters
params <- list(
  min_coverage = 10,
  max_percentile = 99.9,
  diff_threshold = 25,
  qvalue_threshold = 0.01,
  dmr_window_size = 1000,
  dmr_step_size = 1000,
  dmr_min_cpg = 10,
  n_cores = 4
)

# Comparison info
groups <- c(2, 0)  # Treatment2=2, Control=0
labels <- c("Treatment2", "Control")

cat(rep("=", 80), "\n", sep = "")
cat("ANALYSIS: Treatment2 vs Control\n")
cat(rep("=", 80), "\n", sep = "")

# ==============================================================================
# 5mC ANALYSIS
# ==============================================================================

cat("\n### 5mC Analysis ###\n")

# ------------------------------------------------------------------------------
# STEP 1: Subset samples
# ------------------------------------------------------------------------------
cat("\n[STEP 1] Subsetting 5mC samples\n")

selected_idx_mc <- which(mc_methylkit@treatment %in% groups)
mc_subset <- mc_methylkit[selected_idx_mc]

# Reassign treatment (Treatment2=1, Control=0)
new_treatment_mc <- ifelse(mc_methylkit@treatment[selected_idx_mc] == groups[1], 1, 0)
mc_subset@treatment <- new_treatment_mc

cat(sprintf("  Selected %d samples\n", length(mc_subset)))
cat(sprintf("  Treatment2: %d samples\n", sum(new_treatment_mc == 1)))
cat(sprintf("  Control: %d samples\n", sum(new_treatment_mc == 0)))

# ------------------------------------------------------------------------------
# STEP 2: QC
# ------------------------------------------------------------------------------
cat("\n[STEP 2] Quality Control - 5mC\n")

qc_dir_mc <- file.path(comp_dir, "5mC_01_QC")
dir.create(qc_dir_mc, showWarnings = FALSE)

# Methylation stats
pdf(file.path(qc_dir_mc, "methylation_stats.pdf"), width = 15, height = 10)
par(mfrow = c(ceiling(length(mc_subset)/3), 3))
for (i in 1:length(mc_subset)) {
  getMethylationStats(mc_subset[[i]], plot = TRUE, both.strands = FALSE)
  title(main = paste("5mC -", mc_subset[[i]]@sample.id))
}
dev.off()

# Coverage stats
pdf(file.path(qc_dir_mc, "coverage_stats.pdf"), width = 15, height = 10)
par(mfrow = c(ceiling(length(mc_subset)/3), 3))
for (i in 1:length(mc_subset)) {
  getCoverageStats(mc_subset[[i]], plot = TRUE, both.strands = FALSE)
  title(main = paste("5mC -", mc_subset[[i]]@sample.id))
}
dev.off()

# Coverage summary
cov_summary_mc <- data.frame(
  Sample = sapply(mc_subset, function(x) x@sample.id),
  Treatment = sapply(mc_subset, function(x) x@treatment),
  NumCpGs = sapply(mc_subset, function(x) nrow(x)),
  MeanCoverage = sapply(mc_subset, function(x) mean(x$coverage)),
  MedianCoverage = sapply(mc_subset, function(x) median(x$coverage))
)
write.csv(cov_summary_mc, file.path(qc_dir_mc, "coverage_summary.csv"), row.names = FALSE)

cat("  QC complete\n")

# ------------------------------------------------------------------------------
# STEP 3: Filtering
# ------------------------------------------------------------------------------
cat("\n[STEP 3] Filtering 5mC by coverage\n")

mc_filtered <- filterByCoverage(
  mc_subset,
  lo.count = params$min_coverage,
  lo.perc = NULL,
  hi.count = NULL,
  hi.perc = params$max_percentile
)

filter_stats_mc <- data.frame(
  Sample = sapply(mc_filtered, function(x) x@sample.id),
  Before = sapply(mc_subset, function(x) nrow(x)),
  After = sapply(mc_filtered, function(x) nrow(x)),
  Retained_pct = round(sapply(mc_filtered, function(x) nrow(x)) / 
                       sapply(mc_subset, function(x) nrow(x)) * 100, 2)
)
write.csv(filter_stats_mc, file.path(qc_dir_mc, "filtering_stats.csv"), row.names = FALSE)

cat(sprintf("  Average retention: %.1f%%\n", mean(filter_stats_mc$Retained_pct)))

# ------------------------------------------------------------------------------
# STEP 4: Normalization
# ------------------------------------------------------------------------------
cat("\n[STEP 4] Normalizing 5mC\n")

mc_norm <- normalizeCoverage(mc_filtered, method = "median")

# Normalization effect plot
pdf(file.path(qc_dir_mc, "normalization_effect.pdf"), width = 12, height = 6)
par(mfrow = c(1, 2))

boxplot(lapply(mc_filtered, function(x) log10(x$coverage + 1)),
        main = "Before Normalization",
        ylab = "log10(Coverage + 1)",
        las = 2,
        names = sapply(mc_filtered, function(x) x@sample.id))

boxplot(lapply(mc_norm, function(x) log10(x$coverage + 1)),
        main = "After Normalization",
        ylab = "log10(Coverage + 1)",
        las = 2,
        names = sapply(mc_norm, function(x) x@sample.id))
dev.off()

cat("  Normalization complete\n")

# ------------------------------------------------------------------------------
# STEP 5: Unite
# ------------------------------------------------------------------------------
cat("\n[STEP 5] Merging 5mC samples\n")

mc_united <- unite(mc_norm, destrand = FALSE, min.per.group = NULL)

cat(sprintf("  Common CpGs: %d\n", nrow(mc_united)))

save(mc_united, file = file.path(comp_dir, "5mC_united.RData"))

# ------------------------------------------------------------------------------
# STEP 6: Correlation
# ------------------------------------------------------------------------------
cat("\n[STEP 6] Sample correlation - 5mC\n")

corr_dir_mc <- file.path(comp_dir, "5mC_02_Correlation")
dir.create(corr_dir_mc, showWarnings = FALSE)

pdf(file.path(corr_dir_mc, "correlation.pdf"), width = 8, height = 8)
getCorrelation(mc_united, plot = TRUE, method = "pearson")
dev.off()

pdf(file.path(corr_dir_mc, "clustering.pdf"), width = 10, height = 6)
clusterSamples(mc_united, dist = "correlation", method = "ward.D2", plot = TRUE)
dev.off()

cat("  Correlation complete\n")

# ------------------------------------------------------------------------------
# STEP 7: PCA
# ------------------------------------------------------------------------------
cat("\n[STEP 7] PCA - 5mC\n")

pca_dir_mc <- file.path(comp_dir, "5mC_03_PCA")
dir.create(pca_dir_mc, showWarnings = FALSE)

pdf(file.path(pca_dir_mc, "PCA.pdf"), width = 12, height = 6)
par(mfrow = c(1, 2))
PCASamples(mc_united, screeplot = TRUE)
PCASamples(mc_united, screeplot = FALSE)
dev.off()

cat("  PCA complete\n")

# ------------------------------------------------------------------------------
# STEP 8: Differential Methylation
# ------------------------------------------------------------------------------
cat("\n[STEP 8] Differential methylation - 5mC\n")

dm_dir_mc <- file.path(comp_dir, "5mC_04_Differential_Methylation")
dir.create(dm_dir_mc, showWarnings = FALSE)

mc_diff <- calculateDiffMeth(
  mc_united,
  overdispersion = "MN",
  test = "Chisq",
  mc.cores = params$n_cores
)

cat(sprintf("  Tested %d CpGs\n", nrow(mc_diff)))

write.table(mc_diff, file.path(dm_dir_mc, "all_results.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Extract significant sites
mc_hyper <- getMethylDiff(mc_diff, difference = params$diff_threshold,
                          qvalue = params$qvalue_threshold, type = "hyper")
mc_hypo <- getMethylDiff(mc_diff, difference = params$diff_threshold,
                         qvalue = params$qvalue_threshold, type = "hypo")
mc_all <- getMethylDiff(mc_diff, difference = params$diff_threshold,
                        qvalue = params$qvalue_threshold, type = "all")

cat(sprintf("  Hyper: %d, Hypo: %d, Total: %d\n", 
            nrow(mc_hyper), nrow(mc_hypo), nrow(mc_all)))

write.table(mc_hyper, file.path(dm_dir_mc, "hyper_methylated.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mc_hypo, file.path(dm_dir_mc, "hypo_methylated.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mc_all, file.path(dm_dir_mc, "significant_sites.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------------------------------------------------------------
# STEP 9: DMR
# ------------------------------------------------------------------------------
cat("\n[STEP 9] DMR analysis - 5mC\n")

dmr_dir_mc <- file.path(comp_dir, "5mC_05_DMR")
dir.create(dmr_dir_mc, showWarnings = FALSE)

mc_tiles <- tileMethylCounts(
  mc_united,
  win.size = params$dmr_window_size,
  step.size = params$dmr_step_size,
  cov.bases = params$dmr_min_cpg
)

cat(sprintf("  Created %d tiles\n", nrow(mc_tiles)))

mc_dmr_diff <- calculateDiffMeth(mc_tiles, overdispersion = "MN",
                                  test = "Chisq", mc.cores = params$n_cores)

mc_dmr_sig <- getMethylDiff(mc_dmr_diff, difference = params$diff_threshold,
                            qvalue = params$qvalue_threshold, type = "all")

cat(sprintf("  Significant DMRs: %d\n", nrow(mc_dmr_sig)))

if (nrow(mc_dmr_sig) > 0) {
  write.table(mc_dmr_sig, file.path(dmr_dir_mc, "DMR.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  mc_dmr_bed <- data.frame(
    chr = mc_dmr_sig$chr,
    start = mc_dmr_sig$start,
    end = mc_dmr_sig$end,
    name = paste0("DMR_", 1:nrow(mc_dmr_sig)),
    score = round(mc_dmr_sig$meth.diff, 2),
    strand = "."
  )
  
  write.table(mc_dmr_bed, file.path(dmr_dir_mc, "DMR.bed"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# ------------------------------------------------------------------------------
# STEP 10: Visualization
# ------------------------------------------------------------------------------
cat("\n[STEP 10] Creating visualizations - 5mC\n")

viz_dir_mc <- file.path(comp_dir, "5mC_06_Visualization")
dir.create(viz_dir_mc, showWarnings = FALSE)

# Methylation difference distribution
pdf(file.path(viz_dir_mc, "meth_diff_distribution.pdf"), width = 10, height = 6)
hist(mc_diff$meth.diff, breaks = 100,
     main = "5mC Methylation Difference Distribution\nTreatment2 vs Control",
     xlab = "Methylation Difference (%)",
     col = "lightblue", border = "white")
abline(v = c(-params$diff_threshold, params$diff_threshold), 
       col = "red", lty = 2, lwd = 2)
legend("topright", legend = sprintf("Threshold = ±%d%%", params$diff_threshold),
       col = "red", lty = 2, lwd = 2)
dev.off()

# Volcano plot
pdf(file.path(viz_dir_mc, "volcano_plot.pdf"), width = 10, height = 8)

volcano_data_mc <- data.frame(
  meth_diff = mc_diff$meth.diff,
  qvalue = mc_diff$qvalue,
  sig = ifelse(abs(mc_diff$meth.diff) >= params$diff_threshold & 
               mc_diff$qvalue < params$qvalue_threshold,
               ifelse(mc_diff$meth.diff > 0, "Hyper", "Hypo"), "NS")
)

p_volcano_mc <- ggplot(volcano_data_mc, aes(x = meth_diff, y = -log10(qvalue), color = sig)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c("Hyper" = "red", "Hypo" = "blue", "NS" = "gray")) +
  geom_vline(xintercept = c(-params$diff_threshold, params$diff_threshold),
             linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(params$qvalue_threshold),
             linetype = "dashed", color = "black") +
  labs(title = "5mC Volcano Plot\nTreatment2 vs Control",
       x = "Methylation Difference (%)", y = "-log10(q-value)",
       color = "Significance") +
  theme_bw() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

print(p_volcano_mc)
dev.off()

# MA plot
pdf(file.path(viz_dir_mc, "MA_plot.pdf"), width = 10, height = 8)

mean_meth_mc <- rowMeans(getData(mc_united)[, mc_united@numCs.index] / 
                         getData(mc_united)[, mc_united@coverage.index] * 100)

ma_data_mc <- data.frame(
  mean_meth = mean_meth_mc,
  meth_diff = mc_diff$meth.diff,
  sig = volcano_data_mc$sig
)

p_ma_mc <- ggplot(ma_data_mc, aes(x = mean_meth, y = meth_diff, color = sig)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c("Hyper" = "red", "Hypo" = "blue", "NS" = "gray")) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_hline(yintercept = c(-params$diff_threshold, params$diff_threshold),
             linetype = "dashed", color = "black") +
  labs(title = "5mC MA Plot\nTreatment2 vs Control",
       x = "Mean Methylation (%)", y = "Methylation Difference (%)",
       color = "Significance") +
  theme_bw() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

print(p_ma_mc)
dev.off()

# Chromosome distribution
if (nrow(mc_all) > 0) {
  pdf(file.path(viz_dir_mc, "chr_distribution.pdf"), width = 12, height = 6)
  
  chr_counts_mc <- data.frame(table(mc_all$chr))
  colnames(chr_counts_mc) <- c("Chromosome", "Count")
  
  p_chr_mc <- ggplot(chr_counts_mc, aes(x = Chromosome, y = Count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = "5mC: Significant Sites by Chromosome\nTreatment2 vs Control",
         x = "Chromosome", y = "Number of Significant Sites") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  print(p_chr_mc)
  dev.off()
}

cat("  Visualizations complete\n")

# ------------------------------------------------------------------------------
# STEP 11: Summary
# ------------------------------------------------------------------------------
cat("\n[STEP 11] Generating summary - 5mC\n")

summary_mc <- data.frame(
  Metric = c(
    "Samples (Treatment/Control)",
    "Common CpGs after merging",
    "Hyper-methylated sites",
    "Hypo-methylated sites",
    "Total significant sites",
    "Significant DMRs",
    "Mean coverage (Treatment)",
    "Mean coverage (Control)"
  ),
  Value = c(
    sprintf("%d / %d", sum(new_treatment_mc == 1), sum(new_treatment_mc == 0)),
    nrow(mc_united),
    nrow(mc_hyper),
    nrow(mc_hypo),
    nrow(mc_all),
    nrow(mc_dmr_sig),
    sprintf("%.1f", mean(cov_summary_mc$MeanCoverage[cov_summary_mc$Treatment == 1])),
    sprintf("%.1f", mean(cov_summary_mc$MeanCoverage[cov_summary_mc$Treatment == 0]))
  )
)

write.csv(summary_mc, file.path(comp_dir, "5mC_summary.csv"), row.names = FALSE)

print(summary_mc)

# ==============================================================================
# 5hmC ANALYSIS
# ==============================================================================

cat("\n\n### 5hmC Analysis ###\n")

# ------------------------------------------------------------------------------
# STEP 1: Subset samples
# ------------------------------------------------------------------------------
cat("\n[STEP 1] Subsetting 5hmC samples\n")

selected_idx_hmc <- which(hmc_methylkit@treatment %in% groups)
hmc_subset <- hmc_methylkit[selected_idx_hmc]

new_treatment_hmc <- ifelse(hmc_methylkit@treatment[selected_idx_hmc] == groups[1], 1, 0)
hmc_subset@treatment <- new_treatment_hmc

cat(sprintf("  Selected %d samples\n", length(hmc_subset)))
cat(sprintf("  Treatment2: %d samples\n", sum(new_treatment_hmc == 1)))
cat(sprintf("  Control: %d samples\n", sum(new_treatment_hmc == 0)))

# ------------------------------------------------------------------------------
# STEP 2: QC
# ------------------------------------------------------------------------------
cat("\n[STEP 2] Quality Control - 5hmC\n")

qc_dir_hmc <- file.path(comp_dir, "5hmC_01_QC")
dir.create(qc_dir_hmc, showWarnings = FALSE)

pdf(file.path(qc_dir_hmc, "methylation_stats.pdf"), width = 15, height = 10)
par(mfrow = c(ceiling(length(hmc_subset)/3), 3))
for (i in 1:length(hmc_subset)) {
  getMethylationStats(hmc_subset[[i]], plot = TRUE, both.strands = FALSE)
  title(main = paste("5hmC -", hmc_subset[[i]]@sample.id))
}
dev.off()

pdf(file.path(qc_dir_hmc, "coverage_stats.pdf"), width = 15, height = 10)
par(mfrow = c(ceiling(length(hmc_subset)/3), 3))
for (i in 1:length(hmc_subset)) {
  getCoverageStats(hmc_subset[[i]], plot = TRUE, both.strands = FALSE)
  title(main = paste("5hmC -", hmc_subset[[i]]@sample.id))
}
dev.off()

cov_summary_hmc <- data.frame(
  Sample = sapply(hmc_subset, function(x) x@sample.id),
  Treatment = sapply(hmc_subset, function(x) x@treatment),
  NumCpGs = sapply(hmc_subset, function(x) nrow(x)),
  MeanCoverage = sapply(hmc_subset, function(x) mean(x$coverage)),
  MedianCoverage = sapply(hmc_subset, function(x) median(x$coverage))
)
write.csv(cov_summary_hmc, file.path(qc_dir_hmc, "coverage_summary.csv"), row.names = FALSE)

cat("  QC complete\n")

# ------------------------------------------------------------------------------
# STEP 3: Filtering
# ------------------------------------------------------------------------------
cat("\n[STEP 3] Filtering 5hmC by coverage\n")

hmc_filtered <- filterByCoverage(
  hmc_subset,
  lo.count = params$min_coverage,
  lo.perc = NULL,
  hi.count = NULL,
  hi.perc = params$max_percentile
)

filter_stats_hmc <- data.frame(
  Sample = sapply(hmc_filtered, function(x) x@sample.id),
  Before = sapply(hmc_subset, function(x) nrow(x)),
  After = sapply(hmc_filtered, function(x) nrow(x)),
  Retained_pct = round(sapply(hmc_filtered, function(x) nrow(x)) / 
                       sapply(hmc_subset, function(x) nrow(x)) * 100, 2)
)
write.csv(filter_stats_hmc, file.path(qc_dir_hmc, "filtering_stats.csv"), row.names = FALSE)

cat(sprintf("  Average retention: %.1f%%\n", mean(filter_stats_hmc$Retained_pct)))

# ------------------------------------------------------------------------------
# STEP 4: Normalization
# ------------------------------------------------------------------------------
cat("\n[STEP 4] Normalizing 5hmC\n")

hmc_norm <- normalizeCoverage(hmc_filtered, method = "median")

pdf(file.path(qc_dir_hmc, "normalization_effect.pdf"), width = 12, height = 6)
par(mfrow = c(1, 2))

boxplot(lapply(hmc_filtered, function(x) log10(x$coverage + 1)),
        main = "Before Normalization",
        ylab = "log10(Coverage + 1)",
        las = 2,
        names = sapply(hmc_filtered, function(x) x@sample.id))

boxplot(lapply(hmc_norm, function(x) log10(x$coverage + 1)),
        main = "After Normalization",
        ylab = "log10(Coverage + 1)",
        las = 2,
        names = sapply(hmc_norm, function(x) x@sample.id))
dev.off()

cat("  Normalization complete\n")

# ------------------------------------------------------------------------------
# STEP 5: Unite
# ------------------------------------------------------------------------------
cat("\n[STEP 5] Merging 5hmC samples\n")

hmc_united <- unite(hmc_norm, destrand = FALSE, min.per.group = NULL)

cat(sprintf("  Common CpGs: %d\n", nrow(hmc_united)))

save(hmc_united, file = file.path(comp_dir, "5hmC_united.RData"))

# ------------------------------------------------------------------------------
# STEP 6: Correlation
# ------------------------------------------------------------------------------
cat("\n[STEP 6] Sample correlation - 5hmC\n")

corr_dir_hmc <- file.path(comp_dir, "5hmC_02_Correlation")
dir.create(corr_dir_hmc, showWarnings = FALSE)

pdf(file.path(corr_dir_hmc, "correlation.pdf"), width = 8, height = 8)
getCorrelation(hmc_united, plot = TRUE, method = "pearson")
dev.off()

pdf(file.path(corr_dir_hmc, "clustering.pdf"), width = 10, height = 6)
clusterSamples(hmc_united, dist = "correlation", method = "ward.D2", plot = TRUE)
dev.off()

cat("  Correlation complete\n")

# ------------------------------------------------------------------------------
# STEP 7: PCA
# ------------------------------------------------------------------------------
cat("\n[STEP 7] PCA - 5hmC\n")

pca_dir_hmc <- file.path(comp_dir, "5hmC_03_PCA")
dir.create(pca_dir_hmc, showWarnings = FALSE)

pdf(file.path(pca_dir_hmc, "PCA.pdf"), width = 12, height = 6)
par(mfrow = c(1, 2))
PCASamples(hmc_united, screeplot = TRUE)
PCASamples(hmc_united, screeplot = FALSE)
dev.off()

cat("  PCA complete\n")

# ------------------------------------------------------------------------------
# STEP 8: Differential Methylation
# ------------------------------------------------------------------------------
cat("\n[STEP 8] Differential methylation - 5hmC\n")

dm_dir_hmc <- file.path(comp_dir, "5hmC_04_Differential_Methylation")
dir.create(dm_dir_hmc, showWarnings = FALSE)

hmc_diff <- calculateDiffMeth(
  hmc_united,
  overdispersion = "MN",
  test = "Chisq",
  mc.cores = params$n_cores
)

cat(sprintf("  Tested %d CpGs\n", nrow(hmc_diff)))

write.table(hmc_diff, file.path(dm_dir_hmc, "all_results.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

hmc_hyper <- getMethylDiff(hmc_diff, difference = params$diff_threshold,
                           qvalue = params$qvalue_threshold, type = "hyper")
hmc_hypo <- getMethylDiff(hmc_diff, difference = params$diff_threshold,
                          qvalue = params$qvalue_threshold, type = "hypo")
hmc_all <- getMethylDiff(hmc_diff, difference = params$diff_threshold,
                         qvalue = params$qvalue_threshold, type = "all")

cat(sprintf("  Hyper: %d, Hypo: %d, Total: %d\n", 
            nrow(hmc_hyper), nrow(hmc_hypo), nrow(hmc_all)))

write.table(hmc_hyper, file.path(dm_dir_hmc, "hyper_methylated.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(hmc_hypo, file.path(dm_dir_hmc, "hypo_methylated.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(hmc_all, file.path(dm_dir_hmc, "significant_sites.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------------------------------------------------------------
# STEP 9: DMR
# ------------------------------------------------------------------------------
cat("\n[STEP 9] DMR analysis - 5hmC\n")

dmr_dir_hmc <- file.path(comp_dir, "5hmC_05_DMR")
dir.create(dmr_dir_hmc, showWarnings = FALSE)

hmc_tiles <- tileMethylCounts(
  hmc_united,
  win.size = params$dmr_window_size,
  step.size = params$dmr_step_size,
  cov.bases = params$dmr_min_cpg
)

cat(sprintf("  Created %d tiles\n", nrow(hmc_tiles)))

hmc_dmr_diff <- calculateDiffMeth(hmc_tiles, overdispersion = "MN",
                                   test = "Chisq", mc.cores = params$n_cores)

hmc_dmr_sig <- getMethylDiff(hmc_dmr_diff, difference = params$diff_threshold,
                             qvalue = params$qvalue_threshold, type = "all")

cat(sprintf("  Significant DMRs: %d\n", nrow(hmc_dmr_sig)))

if (nrow(hmc_dmr_sig) > 0) {
  write.table(hmc_dmr_sig, file.path(dmr_dir_hmc, "DMR.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  hmc_dmr_bed <- data.frame(
    chr = hmc_dmr_sig$chr,
    start = hmc_dmr_sig$start,
    end = hmc_dmr_sig$end,
    name = paste0("DMR_", 1:nrow(hmc_dmr_sig)),
    score = round(hmc_dmr_sig$meth.diff, 2),
    strand = "."
  )
  
  write.table(hmc_dmr_bed, file.path(dmr_dir_hmc, "DMR.bed"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# ------------------------------------------------------------------------------
# STEP 10: Visualization
# ------------------------------------------------------------------------------
cat("\n[STEP 10] Creating visualizations - 5hmC\n")

viz_dir_hmc <- file.path(comp_dir, "5hmC_06_Visualization")
dir.create(viz_dir_hmc, showWarnings = FALSE)

# Methylation difference distribution
pdf(file.path(viz_dir_hmc, "meth_diff_distribution.pdf"), width = 10, height = 6)
hist(hmc_diff$meth.diff, breaks = 100,
     main = "5hmC Methylation Difference Distribution\nTreatment2 vs Control",
     xlab = "Methylation Difference (%)",
     col = "lightblue", border = "white")
abline(v = c(-params$diff_threshold, params$diff_threshold),
       col = "red", lty = 2, lwd = 2)
legend("topright", legend = sprintf("Threshold = ±%d%%", params$diff_threshold),
       col = "red", lty = 2, lwd = 2)
dev.off()

# Volcano plot
pdf(file.path(viz_dir_hmc, "volcano_plot.pdf"), width = 10, height = 8)

volcano_data_hmc <- data.frame(
  meth_diff = hmc_diff$meth.diff,
  qvalue    = hmc_diff$qvalue,
  sig       = ifelse(
    abs(hmc_diff$meth.diff) >= params$diff_threshold &
      hmc_diff$qvalue < params$qvalue_threshold,
    ifelse(hmc_diff$meth.diff > 0, "Hyper", "Hypo"),
    "NS"
  )
)

p_volcano_hmc <- ggplot(volcano_data_hmc,
                        aes(x = meth_diff, y = -log10(qvalue), color = sig)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c("Hyper" = "red",
                                "Hypo"  = "blue",
                                "NS"    = "gray")) +
  geom_vline(xintercept = c(-params$diff_threshold, params$diff_threshold),
             linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(params$qvalue_threshold),
             linetype = "dashed", color = "black") +
  labs(title = "5hmC Volcano Plot\nTreatment2 vs Control",
       x = "Methylation Difference (%)", y = "-log10(q-value)",
       color = "Significance") +
  theme_bw() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

print(p_volcano_hmc)
dev.off()

# MA plot
pdf(file.path(viz_dir_hmc, "MA_plot.pdf"), width = 10, height = 8)

mean_meth_hmc <- rowMeans(
  getData(hmc_united)[, hmc_united@numCs.index] /
    getData(hmc_united)[, hmc_united@coverage.index] * 100
)

ma_data_hmc <- data.frame(
  mean_meth = mean_meth_hmc,
  meth_diff = hmc_diff$meth.diff,
  sig       = volcano_data_hmc$sig
)

p_ma_hmc <- ggplot(ma_data_hmc,
                   aes(x = mean_meth, y = meth_diff, color = sig)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c("Hyper" = "red",
                                "Hypo"  = "blue",
                                "NS"    = "gray")) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_hline(yintercept = c(-params$diff_threshold, params$diff_threshold),
             linetype = "dashed", color = "black") +
  labs(title = "5hmC MA Plot\nTreatment2 vs Control",
       x = "Mean Methylation (%)", y = "Methylation Difference (%)",
       color = "Significance") +
  theme_bw() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

print(p_ma_hmc)
dev.off()

# Chromosome distribution
if (nrow(hmc_all) > 0) {
  pdf(file.path(viz_dir_hmc, "chr_distribution.pdf"), width = 12, height = 6)

  chr_counts_hmc <- data.frame(table(hmc_all$chr))
  colnames(chr_counts_hmc) <- c("Chromosome", "Count")

  p_chr_hmc <- ggplot(chr_counts_hmc, aes(x = Chromosome, y = Count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = "5hmC: Significant Sites by Chromosome\nTreatment2 vs Control",
         x = "Chromosome", y = "Number of Significant Sites") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

  print(p_chr_hmc)
  dev.off()
}

cat("  Visualizations complete\n")

# ------------------------------------------------------------------------------
# STEP 11: Summary
# ------------------------------------------------------------------------------
cat("\n[STEP 11] Generating summary - 5hmC\n")

summary_hmc <- data.frame(
  Metric = c(
    "Samples (Treatment/Control)",
    "Common CpGs after merging",
    "Hyper-hydroxymethylated sites",
    "Hypo-hydroxymethylated sites",
    "Total significant sites",
    "Significant DMRs",
    "Mean coverage (Treatment)",
    "Mean coverage (Control)"
  ),
  Value = c(
    sprintf("%d / %d", sum(new_treatment_hmc == 1), sum(new_treatment_hmc == 0)),
    nrow(hmc_united),
    nrow(hmc_hyper),
    nrow(hmc_hypo),
    nrow(hmc_all),
    nrow(hmc_dmr_sig),
    sprintf("%.1f",
            mean(cov_summary_hmc$MeanCoverage[cov_summary_hmc$Treatment == 1])),
    sprintf("%.1f",
            mean(cov_summary_hmc$MeanCoverage[cov_summary_hmc$Treatment == 0]))
  )
)

write.csv(summary_hmc, file.path(comp_dir, "5hmC_summary.csv"), row.names = FALSE)
print(summary_hmc)

# ==============================================================================
# Combined Summary
# ==============================================================================
cat("\n", rep("=", 80), "\n", sep = "")
cat("COMBINED SUMMARY: Treatment2 vs Control\n")
cat(rep("=", 80), "\n", sep = "")

combined_summary <- rbind(
  data.frame(Modification = "5mC", summary_mc),
  data.frame(Modification = "5hmC", summary_hmc)
)

print(combined_summary)

write.csv(
  combined_summary,
  file.path(comp_dir, "Combined_Summary.csv"),
  row.names = FALSE
)

# Combined visualization
pdf(file.path(comp_dir, "Combined_Volcano_Plot.pdf"), width = 14, height = 6)
gridExtra::grid.arrange(p_volcano_mc, p_volcano_hmc, ncol = 2)
dev.off()

pdf(file.path(comp_dir, "Combined_MA_Plot.pdf"), width = 14, height = 6)
gridExtra::grid.arrange(p_ma_mc, p_ma_hmc, ncol = 2)
dev.off()

# Save all results
save(mc_united, mc_diff, mc_hyper, mc_hypo, mc_all, mc_dmr_sig,
     hmc_united, hmc_diff, hmc_hyper, hmc_hypo, hmc_all, hmc_dmr_sig,
     file = file.path(comp_dir, "All_Results.RData"))

cat("\n", rep("=", 80), "\n", sep = "")
cat("ANALYSIS COMPLETE: Treatment2 vs Control\n")
cat(rep("=", 80), "\n", sep = "")
cat(sprintf("\nResults saved in: %s\n", comp_dir))
cat("\nDirectory structure:\n")
cat(sprintf("  %s/\n", comp_dir))
cat("  ├── 5mC_01_QC/\n")
cat("  ├── 5mC_02_Correlation/\n")
cat("  ├── 5mC_03_PCA/\n")
cat("  ├── 5mC_04_Differential_Methylation/\n")
cat("  ├── 5mC_05_DMR/\n")
cat("  ├── 5mC_06_Visualization/\n")
cat("  ├── 5hmC_01_QC/\n")
cat("  ├── 5hmC_02_Correlation/\n")
cat("  ├── 5hmC_03_PCA/\n")
cat("  ├── 5hmC_04_Differential_Methylation/\n")
cat("  ├── 5hmC_05_DMR/\n")
cat("  ├── 5hmC_06_Visualization/\n")
cat("  ├── 5mC_summary.csv\n")
cat("  ├── 5hmC_summary.csv\n")
cat("  ├── Combined_Summary.csv\n")
cat("  ├── Combined_Volcano_Plot.pdf\n")
cat("  ├── Combined_MA_Plot.pdf\n")
cat("  └── All_Results.RData\n")
cat("\n")
