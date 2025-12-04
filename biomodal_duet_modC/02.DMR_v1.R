# ==============================================================================
# DMR Analysis Pipeline with Complete Analysis for Each Pairwise Comparison
# ==============================================================================
# Treatment groups: 0=Control, 1=Treatment1, 2=Treatment2
# Comparisons: T2 vs C, T1 vs C, T2 vs T1
# Each comparison gets full QC, filtering, normalization, DM analysis
# ==============================================================================

library(methylKit)
library(dplyr)
library(ggplot2)
library(gridExtra)

# ------------------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------------------

# 출력 디렉토리 생성
output_dir <- "DMR_Analysis_Results"
dir.create(output_dir, showWarnings = FALSE)

# 비교 조합 정의
comparisons <- list(
  list(name = "Treatment2_vs_Control", 
       groups = c(2, 0),
       labels = c("Treatment2", "Control")),
  list(name = "Treatment1_vs_Control", 
       groups = c(1, 0),
       labels = c("Treatment1", "Control")),
  list(name = "Treatment2_vs_Treatment1", 
       groups = c(2, 1),
       labels = c("Treatment2", "Treatment1"))
)

# 분석 파라미터
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

# ------------------------------------------------------------------------------
# 1. Load methylKit objects
# ------------------------------------------------------------------------------

cat("=== Loading methylKit objects ===\n")
# load("methylkit_objects.RData")  # 필요시

cat("\n5mC methylRawList:\n")
print(mc_methylkit)

cat("\n5hmC methylRawList:\n")
print(hmc_methylkit)

cat("\nTreatment groups:\n")
print(table(mc_methylkit@treatment))

# ------------------------------------------------------------------------------
# 2. Function: Complete analysis for one comparison
# ------------------------------------------------------------------------------

perform_complete_analysis <- function(meth_obj, comp_info, mod_type, params, output_dir) {
  
  comp_name <- comp_info$name
  groups <- comp_info$groups
  labels <- comp_info$labels
  
  cat("\n", rep("=", 80), "\n", sep = "")
  cat(sprintf("ANALYSIS: %s - %s\n", mod_type, comp_name))
  cat(sprintf("Comparing: %s (treatment=1) vs %s (treatment=0)\n", labels[1], labels[2]))
  cat(rep("=", 80), "\n", sep = "")
  
  # 출력 디렉토리 생성
  comp_dir <- file.path(output_dir, comp_name)
  dir.create(comp_dir, showWarnings = FALSE)
  
  # ============================================================================
  # STEP 1: Subset samples for this comparison
  # ============================================================================
  cat("\n[STEP 1] Subsetting samples\n")
  
  selected_idx <- which(meth_obj@treatment %in% groups)
  selected_samples <- meth_obj[selected_idx]
  
  # Treatment 재할당 (첫번째 그룹=1, 두번째 그룹=0)
  new_treatment <- ifelse(meth_obj@treatment[selected_idx] == groups[1], 1, 0)
  selected_samples@treatment <- new_treatment
  
  cat(sprintf("  Selected %d samples\n", length(selected_samples)))
  cat(sprintf("  Treatment group (%s): %d samples\n", 
              labels[1], sum(new_treatment == 1)))
  cat(sprintf("  Control group (%s): %d samples\n", 
              labels[2], sum(new_treatment == 0)))
  
  # ============================================================================
  # STEP 2: QC - Individual samples
  # ============================================================================
  cat("\n[STEP 2] Quality Control\n")
  
  qc_dir <- file.path(comp_dir, "01_QC")
  dir.create(qc_dir, showWarnings = FALSE)
  
  # QC for each sample
  pdf(file.path(qc_dir, sprintf("%s_methylation_stats.pdf", mod_type)), 
      width = 15, height = 10)
  par(mfrow = c(ceiling(length(selected_samples)/3), 3))
  for (i in 1:length(selected_samples)) {
    getMethylationStats(selected_samples[[i]], plot = TRUE, both.strands = FALSE)
    title(main = paste(mod_type, "-", selected_samples[[i]]@sample.id))
  }
  dev.off()
  
  pdf(file.path(qc_dir, sprintf("%s_coverage_stats.pdf", mod_type)), 
      width = 15, height = 10)
  par(mfrow = c(ceiling(length(selected_samples)/3), 3))
  for (i in 1:length(selected_samples)) {
    getCoverageStats(selected_samples[[i]], plot = TRUE, both.strands = FALSE)
    title(main = paste(mod_type, "-", selected_samples[[i]]@sample.id))
  }
  dev.off()
  
  # Coverage summary
  cov_summary <- data.frame(
    Sample = sapply(selected_samples, function(x) x@sample.id),
    Treatment = sapply(selected_samples, function(x) x@treatment),
    NumCpGs = sapply(selected_samples, function(x) nrow(x)),
    MeanCoverage = sapply(selected_samples, function(x) mean(x$coverage)),
    MedianCoverage = sapply(selected_samples, function(x) median(x$coverage))
  )
  write.csv(cov_summary, 
            file.path(qc_dir, sprintf("%s_coverage_summary.csv", mod_type)),
            row.names = FALSE)
  
  cat("  QC plots and summary saved\n")
  
  # ============================================================================
  # STEP 3: Filtering by coverage
  # ============================================================================
  cat("\n[STEP 3] Filtering by coverage\n")
  cat(sprintf("  Min coverage: %d\n", params$min_coverage))
  cat(sprintf("  Max percentile: %.1f\n", params$max_percentile))
  
  filtered <- filterByCoverage(
    selected_samples,
    lo.count = params$min_coverage,
    lo.perc = NULL,
    hi.count = NULL,
    hi.perc = params$max_percentile
  )
  
  # Filtering stats
  filter_stats <- data.frame(
    Sample = sapply(filtered, function(x) x@sample.id),
    Before = sapply(selected_samples, function(x) nrow(x)),
    After = sapply(filtered, function(x) nrow(x)),
    Retained_pct = round(sapply(filtered, function(x) nrow(x)) / 
                         sapply(selected_samples, function(x) nrow(x)) * 100, 2)
  )
  write.csv(filter_stats,
            file.path(qc_dir, sprintf("%s_filtering_stats.csv", mod_type)),
            row.names = FALSE)
  
  cat(sprintf("  Average retention: %.1f%%\n", mean(filter_stats$Retained_pct)))
  
  # ============================================================================
  # STEP 4: Normalization
  # ============================================================================
  cat("\n[STEP 4] Normalization (median)\n")
  
  normalized <- normalizeCoverage(filtered, method = "median")
  
  # Normalization effect plot
  pdf(file.path(qc_dir, sprintf("%s_normalization_effect.pdf", mod_type)),
      width = 12, height = 6)
  par(mfrow = c(1, 2))
  
  # Before normalization
  boxplot(lapply(filtered, function(x) log10(x$coverage + 1)),
          main = "Before Normalization",
          ylab = "log10(Coverage + 1)",
          las = 2,
          names = sapply(filtered, function(x) x@sample.id))
  
  # After normalization
  boxplot(lapply(normalized, function(x) log10(x$coverage + 1)),
          main = "After Normalization",
          ylab = "log10(Coverage + 1)",
          las = 2,
          names = sapply(normalized, function(x) x@sample.id))
  dev.off()
  
  cat("  Normalization complete\n")
  
  # ============================================================================
  # STEP 5: Merge samples (unite)
  # ============================================================================
  cat("\n[STEP 5] Merging samples\n")
  
  united <- unite(
    normalized,
    destrand = FALSE,
    min.per.group = NULL
  )
  
  cat(sprintf("  Common CpGs: %d\n", nrow(united)))
  
  # Save united object
  save(united, file = file.path(comp_dir, sprintf("%s_united.RData", mod_type)))
  
  # ============================================================================
  # STEP 6: Sample correlation and clustering
  # ============================================================================
  cat("\n[STEP 6] Sample correlation and clustering\n")
  
  corr_dir <- file.path(comp_dir, "02_Correlation")
  dir.create(corr_dir, showWarnings = FALSE)
  
  # Correlation heatmap
  pdf(file.path(corr_dir, sprintf("%s_correlation.pdf", mod_type)),
      width = 8, height = 8)
  getCorrelation(united, plot = TRUE, method = "pearson")
  dev.off()
  
  # Hierarchical clustering
  pdf(file.path(corr_dir, sprintf("%s_clustering.pdf", mod_type)),
      width = 10, height = 6)
  clusterSamples(united, dist = "correlation", method = "ward.D2", plot = TRUE)
  dev.off()
  
  cat("  Correlation analysis complete\n")
  
  # ============================================================================
  # STEP 7: PCA
  # ============================================================================
  cat("\n[STEP 7] Principal Component Analysis\n")
  
  pca_dir <- file.path(comp_dir, "03_PCA")
  dir.create(pca_dir, showWarnings = FALSE)
  
  pdf(file.path(pca_dir, sprintf("%s_PCA.pdf", mod_type)),
      width = 12, height = 6)
  par(mfrow = c(1, 2))
  PCASamples(united, screeplot = TRUE)
  PCASamples(united, screeplot = FALSE)
  dev.off()
  
  cat("  PCA complete\n")
  
  # ============================================================================
  # STEP 8: Differential Methylation - Base level
  # ============================================================================
  cat("\n[STEP 8] Differential methylation analysis (base level)\n")
  
  dm_dir <- file.path(comp_dir, "04_Differential_Methylation")
  dir.create(dm_dir, showWarnings = FALSE)
  
  diff_meth <- calculateDiffMeth(
    united,
    overdispersion = "MN",
    test = "Chisq",
    mc.cores = params$n_cores
  )
  
  cat(sprintf("  Tested %d CpGs\n", nrow(diff_meth)))
  
  # Save all results
  write.table(diff_meth,
              file.path(dm_dir, sprintf("%s_all_results.txt", mod_type)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Extract significant sites
  diff_hyper <- getMethylDiff(
    diff_meth,
    difference = params$diff_threshold,
    qvalue = params$qvalue_threshold,
    type = "hyper"
  )
  
  diff_hypo <- getMethylDiff(
    diff_meth,
    difference = params$diff_threshold,
    qvalue = params$qvalue_threshold,
    type = "hypo"
  )
  
  diff_all <- getMethylDiff(
    diff_meth,
    difference = params$diff_threshold,
    qvalue = params$qvalue_threshold,
    type = "all"
  )
  
  cat(sprintf("  Significant sites:\n"))
  cat(sprintf("    Hyper: %d\n", nrow(diff_hyper)))
  cat(sprintf("    Hypo: %d\n", nrow(diff_hypo)))
  cat(sprintf("    Total: %d\n", nrow(diff_all)))
  
  # Save significant sites
  write.table(diff_hyper,
              file.path(dm_dir, sprintf("%s_hyper_methylated.txt", mod_type)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  write.table(diff_hypo,
              file.path(dm_dir, sprintf("%s_hypo_methylated.txt", mod_type)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  write.table(diff_all,
              file.path(dm_dir, sprintf("%s_significant_sites.txt", mod_type)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # ============================================================================
  # STEP 9: Differential Methylation - Region level (DMR)
  # ============================================================================
  cat("\n[STEP 9] Differential methylation regions (DMR)\n")
  
  dmr_dir <- file.path(comp_dir, "05_DMR")
  dir.create(dmr_dir, showWarnings = FALSE)
  
  # Tile windows
  tiles <- tileMethylCounts(
    united,
    win.size = params$dmr_window_size,
    step.size = params$dmr_step_size,
    cov.bases = params$dmr_min_cpg
  )
  
  cat(sprintf("  Created %d tiles\n", nrow(tiles)))
  
  # Calculate differential methylation for regions
  dmr_diff <- calculateDiffMeth(
    tiles,
    overdispersion = "MN",
    test = "Chisq",
    mc.cores = params$n_cores
  )
  
  # Extract significant DMRs
  dmr_sig <- getMethylDiff(
    dmr_diff,
    difference = params$diff_threshold,
    qvalue = params$qvalue_threshold,
    type = "all"
  )
  
  cat(sprintf("  Significant DMRs: %d\n", nrow(dmr_sig)))
  
  if (nrow(dmr_sig) > 0) {
    # Save DMR results
    write.table(dmr_sig,
                file.path(dmr_dir, sprintf("%s_DMR.txt", mod_type)),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Save as BED format
    dmr_bed <- data.frame(
      chr = dmr_sig$chr,
      start = dmr_sig$start,
      end = dmr_sig$end,
      name = paste0("DMR_", 1:nrow(dmr_sig)),
      score = round(dmr_sig$meth.diff, 2),
      strand = "."
    )
    
    write.table(dmr_bed,
                file.path(dmr_dir, sprintf("%s_DMR.bed", mod_type)),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  # ============================================================================
  # STEP 10: Visualization
  # ============================================================================
  cat("\n[STEP 10] Creating visualizations\n")
  
  viz_dir <- file.path(comp_dir, "06_Visualization")
  dir.create(viz_dir, showWarnings = FALSE)
  
  # 10.1 Methylation difference distribution
  pdf(file.path(viz_dir, sprintf("%s_meth_diff_distribution.pdf", mod_type)),
      width = 10, height = 6)
  hist(diff_meth$meth.diff,
       breaks = 100,
       main = sprintf("%s Methylation Difference Distribution\n%s", 
                     mod_type, comp_name),
       xlab = "Methylation Difference (%)",
       col = "lightblue",
       border = "white")
  abline(v = c(-params$diff_threshold, params$diff_threshold), 
         col = "red", lty = 2, lwd = 2)
  legend("topright", 
         legend = sprintf("Threshold = ±%d%%", params$diff_threshold),
         col = "red", lty = 2, lwd = 2)
  dev.off()
  
  # 10.2 Volcano plot
  pdf(file.path(viz_dir, sprintf("%s_volcano_plot.pdf", mod_type)),
      width = 10, height = 8)
  
  volcano_data <- data.frame(
    meth_diff = diff_meth$meth.diff,
    qvalue = diff_meth$qvalue,
    sig = ifelse(abs(diff_meth$meth.diff) >= params$diff_threshold & 
                 diff_meth$qvalue < params$qvalue_threshold,
                 ifelse(diff_meth$meth.diff > 0, "Hyper", "Hypo"),
                 "NS")
  )
  
  p_volcano <- ggplot(volcano_data, aes(x = meth_diff, y = -log10(qvalue), color = sig)) +
    geom_point(alpha = 0.5, size = 1.5) +
    scale_color_manual(values = c("Hyper" = "red", "Hypo" = "blue", "NS" = "gray")) +
    geom_vline(xintercept = c(-params$diff_threshold, params$diff_threshold),
               linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(params$qvalue_threshold),
               linetype = "dashed", color = "black") +
    labs(title = sprintf("%s Volcano Plot\n%s", mod_type, comp_name),
         x = "Methylation Difference (%)",
         y = "-log10(q-value)",
         color = "Significance") +
    theme_bw() +
    theme(legend.position = "top",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  print(p_volcano)
  dev.off()
  
  # 10.3 MA plot
  pdf(file.path(viz_dir, sprintf("%s_MA_plot.pdf", mod_type)),
      width = 10, height = 8)
  
  # Calculate mean methylation
  mean_meth <- rowMeans(getData(united)[, united@numCs.index] / 
                        getData(united)[, united@coverage.index] * 100)
  
  ma_data <- data.frame(
    mean_meth = mean_meth,
    meth_diff = diff_meth$meth.diff,
    sig = volcano_data$sig
  )
  
  p_ma <- ggplot(ma_data, aes(x = mean_meth, y = meth_diff, color = sig)) +
    geom_point(alpha = 0.5, size = 1.5) +
    scale_color_manual(values = c("Hyper" = "red", "Hypo" = "blue", "NS" = "gray")) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black") +
    geom_hline(yintercept = c(-params$diff_threshold, params$diff_threshold),
               linetype = "dashed", color = "black") +
    labs(title = sprintf("%s MA Plot\n%s", mod_type, comp_name),
         x = "Mean Methylation (%)",
         y = "Methylation Difference (%)",
         color = "Significance") +
    theme_bw() +
    theme(legend.position = "top",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  print(p_ma)
  dev.off()
  
  # 10.4 Chromosome distribution
  if (nrow(diff_all) > 0) {
    pdf(file.path(viz_dir, sprintf("%s_chr_distribution.pdf", mod_type)),
        width = 12, height = 6)
    
    chr_counts <- data.frame(table(diff_all$chr))
    colnames(chr_counts) <- c("Chromosome", "Count")
    
    p_chr <- ggplot(chr_counts, aes(x = Chromosome, y = Count)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      labs(title = sprintf("%s: Significant Sites by Chromosome\n%s", 
                          mod_type, comp_name),
           x = "Chromosome",
           y = "Number of Significant Sites") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    
    print(p_chr)
    dev.off()
  }
  
  cat("  Visualizations complete\n")
  
  # ============================================================================
  # STEP 11: Summary report
  # ============================================================================
  cat("\n[STEP 11] Generating summary report\n")
  
  summary_stats <- data.frame(
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
      sprintf("%d / %d", sum(new_treatment == 1), sum(new_treatment == 0)),
      nrow(united),
      nrow(diff_hyper),
      nrow(diff_hypo),
      nrow(diff_all),
      nrow(dmr_sig),
      sprintf("%.1f", mean(cov_summary$MeanCoverage[cov_summary$Treatment == 1])),
      sprintf("%.1f", mean(cov_summary$MeanCoverage[cov_summary$Treatment == 0]))
    )
  )
  
  write.csv(summary_stats,
            file.path(comp_dir, sprintf("%s_summary.csv", mod_type)),
            row.names = FALSE)
  
  print(summary_stats)
  
  # ============================================================================
  # Return results
  # ============================================================================
  
  return(list(
    united = united,
    diff = diff_meth,
    hyper = diff_hyper,
    hypo = diff_hypo,
    all = diff_all,
    dmr = dmr_sig,
    summary = summary_stats
  ))
}

# ------------------------------------------------------------------------------
# 3. Run complete analysis for all comparisons
# ------------------------------------------------------------------------------

cat("\n", rep("=", 80), "\n", sep = "")
cat("STARTING COMPLETE ANALYSIS FOR ALL COMPARISONS\n")
cat(rep("=", 80), "\n", sep = "")

# 결과 저장용
all_results <- list(
  mc = list(),
  hmc = list()
)

for (comp in comparisons) {
  
  # 5mC analysis
  all_results$mc[[comp$name]] <- perform_complete_analysis(
    meth_obj = mc_methylkit,
    comp_info = comp,
    mod_type = "5mC",
    params = params,
    output_dir = output_dir
  )
  
  # 5hmC analysis
  all_results$hmc[[comp$name]] <- perform_complete_analysis(
    meth_obj = hmc_methylkit,
    comp_info = comp,
    mod_type = "5hmC",
    params = params,
    output_dir = output_dir
  )
}

# ------------------------------------------------------------------------------
# 4. Cross-comparison summary
# ------------------------------------------------------------------------------

cat("\n", rep("=", 80), "\n", sep = "")
cat("GENERATING CROSS-COMPARISON SUMMARY\n")
cat(rep("=", 80), "\n", sep = "")

# Compile summary across all comparisons
summary_list <- list()

for (comp_name in names(all_results$mc)) {
  summary_list[[paste0(comp_name, "_5mC")]] <- data.frame(
    Comparison = comp_name,
    Modification = "5mC",
    Common_CpGs = nrow(all_results$mc[[comp_name]]$united),
    Hyper = nrow(all_results$mc[[comp_name]]$hyper),
    Hypo = nrow(all_results$mc[[comp_name]]$hypo),
    Total_Sig = nrow(all_results$mc[[comp_name]]$all),
    DMR_count = nrow(all_results$mc[[comp_name]]$dmr)
  )
  
  summary_list[[paste0(comp_name, "_5hmC")]] <- data.frame(
    Comparison = comp_name,
    Modification = "5hmC",
    Common_CpGs = nrow(all_results$hmc[[comp_name]]$united),
    Hyper = nrow(all_results$hmc[[comp_name]]$hyper),
    Hypo = nrow(all_results$hmc[[comp_name]]$hypo),
    Total_Sig = nrow(all_results$hmc[[comp_name]]$all),
    DMR_count = nrow(all_results$hmc[[comp_name]]$dmr)
  )
}

summary_all <- do.call(rbind, summary_list)
rownames(summary_all) <- NULL

print(summary_all)

write.csv(summary_all,
          file.path(output_dir, "Summary_All_Comparisons.csv"),
          row.names = FALSE)

# Summary barplot
pdf(file.path(output_dir, "Summary_Barplot.pdf"), width = 14, height = 8)

library(tidyr)
summary_plot <- summary_all %>%
  select(Comparison, Modification, Hyper, Hypo) %>%
  pivot_longer(cols = c(Hyper, Hypo), names_to = "Direction", values_to = "Count")

p_summary <- ggplot(summary_plot, aes(x = Comparison, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Modification, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) +
  labs(title = "Differential Methylation Summary Across All Comparisons",
       y = "Number of Sites",
       fill = "Direction") +
  scale_fill_manual(values = c("Hyper" = "red", "Hypo" = "blue"))

print(p_summary)
dev.off()

# ------------------------------------------------------------------------------
# 5. Save complete results
# ------------------------------------------------------------------------------

save(all_results, comparisons, params,
     file = file.path(output_dir, "All_Results.RData"))

cat("\n", rep("=", 80), "\n", sep = "")
cat("ANALYSIS COMPLETE!\n")
cat(rep("=", 80), "\n", sep = "")
cat(sprintf("\nAll results saved in: %s/\n", output_dir))
cat("\nDirectory structure:\n")
cat("  DMR_Analysis_Results/\n")
cat("  ├── Summary_All_Comparisons.csv\n")
cat("  ├── Summary_Barplot.pdf\n")
cat("  ├── All_Results.RData\n")
for (comp in comparisons) {
  cat(sprintf("  ├── %s/\n", comp$name))
  cat("  │   ├── 01_QC/\n")
  cat("  │   ├── 02_Correlation/\n")
  cat("  │   ├── 03_PCA/\n")
  cat("  │   ├── 04_Differential_Methylation/\n")
  cat("  │   ├── 05_DMR/\n")
  cat("  │   └── 06_Visualization/\n")
}

cat("\n")



QC (개별 샘플별)
Filtering 및 통계
Normalization 효과 확인
Sample correlation/clustering
PCA
Base-level DM
DMR 분석
완전한 시각화 (volcano, MA, chromosome distribution)
요약 리포트
