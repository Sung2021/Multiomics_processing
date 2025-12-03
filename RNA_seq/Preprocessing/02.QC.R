# 필요한 패키지
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(reshape2)

# ============================================================================
# 데이터 준비
# ============================================================================

# Count matrix 준비 (이전 코드로 생성된 merged_counts.txt 사용)
count_data <- read.delim("merged_counts.txt", stringsAsFactors = FALSE)

# Count matrix 추출 (gene annotation 제외)
count_matrix <- count_data %>%
  select(-Geneid, -Chr, -Start, -End, -Strand, -Length) %>%
  as.matrix()

rownames(count_matrix) <- count_data$Geneid

# 샘플명
sample_names <- colnames(count_matrix)

# ============================================================================
# 1. Basic Data Quality Checks
# ============================================================================

# 1.1 Library Size Analysis
# ----------------------------------------------------------------------------
library_sizes <- colSums(count_matrix)

# Barplot
pdf("1.1_library_size.pdf", width = 10, height = 6)
par(mar = c(8, 5, 4, 2))
barplot(library_sizes / 1e6, 
        names.arg = sample_names,
        las = 2,
        ylab = "Library Size (millions)",
        main = "Library Size per Sample",
        col = "steelblue")
abline(h = median(library_sizes) / 1e6, col = "red", lty = 2)
legend("topright", legend = "Median", col = "red", lty = 2)
dev.off()

# Log scale barplot
pdf("1.1_library_size_log.pdf", width = 10, height = 6)
par(mar = c(8, 5, 4, 2))
barplot(log10(library_sizes), 
        names.arg = sample_names,
        las = 2,
        ylab = "log10(Library Size)",
        main = "Library Size per Sample (log scale)",
        col = "steelblue")
dev.off()

# Summary statistics
library_size_summary <- data.frame(
  Sample = sample_names,
  Total_Counts = library_sizes,
  Million_Reads = library_sizes / 1e6
)
write.table(library_size_summary, "1.1_library_size_summary.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


# 1.2 Gene Detection Rate
# ----------------------------------------------------------------------------
detected_genes_0 <- colSums(count_matrix > 0)
detected_genes_10 <- colSums(count_matrix > 10)

# Barplot
pdf("1.2_gene_detection_rate.pdf", width = 12, height = 6)
par(mfrow = c(1, 2), mar = c(8, 5, 4, 2))

barplot(detected_genes_0,
        names.arg = sample_names,
        las = 2,
        ylab = "Number of Genes",
        main = "Detected Genes (count > 0)",
        col = "darkgreen")
abline(h = median(detected_genes_0), col = "red", lty = 2)

barplot(detected_genes_10,
        names.arg = sample_names,
        las = 2,
        ylab = "Number of Genes",
        main = "Detected Genes (count > 10)",
        col = "darkgreen")
abline(h = median(detected_genes_10), col = "red", lty = 2)
dev.off()

# Summary
detection_summary <- data.frame(
  Sample = sample_names,
  Detected_count_gt_0 = detected_genes_0,
  Detected_count_gt_10 = detected_genes_10,
  Detection_Rate_0 = detected_genes_0 / nrow(count_matrix) * 100,
  Detection_Rate_10 = detected_genes_10 / nrow(count_matrix) * 100
)
write.table(detection_summary, "1.2_gene_detection_summary.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


# 1.3 Expression Distribution Analysis
# ----------------------------------------------------------------------------
log_counts <- log2(count_matrix + 1)

# Boxplot
pdf("1.3_expression_distribution_boxplot.pdf", width = 12, height = 6)
par(mar = c(8, 5, 4, 2))
boxplot(log_counts,
        names = sample_names,
        las = 2,
        ylab = "log2(counts + 1)",
        main = "Expression Distribution per Sample",
        col = "lightblue",
        outline = FALSE)
dev.off()

# Density plot
pdf("1.3_expression_distribution_density.pdf", width = 10, height = 6)
plot(density(log_counts[, 1]), 
     xlim = c(0, max(log_counts)),
     ylim = c(0, max(sapply(1:ncol(log_counts), function(i) max(density(log_counts[, i])$y)))),
     main = "Expression Distribution Density",
     xlab = "log2(counts + 1)",
     ylab = "Density",
     col = 1,
     lwd = 2)
for (i in 2:ncol(log_counts)) {
  lines(density(log_counts[, i]), col = i, lwd = 2)
}
legend("topright", legend = sample_names, col = 1:ncol(log_counts), lwd = 2, cex = 0.7)
dev.off()

# Violin plot (ggplot2)
log_counts_long <- melt(log_counts)
colnames(log_counts_long) <- c("Gene", "Sample", "Expression")

p <- ggplot(log_counts_long, aes(x = Sample, y = Expression, fill = Sample)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "Expression Distribution (Violin Plot)",
       y = "log2(counts + 1)")
ggsave("1.3_expression_distribution_violin.pdf", p, width = 12, height = 6)


# ============================================================================
# 2. Further Quality Control Metrics
# ============================================================================

# 2.1 Outlier Detection by Distance
# ----------------------------------------------------------------------------

# Sample correlation matrix
cor_matrix <- cor(count_matrix, method = "pearson")

# Correlation heatmap
pdf("2.1_sample_correlation_heatmap.pdf", width = 10, height = 9)
pheatmap(cor_matrix,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         display_numbers = TRUE,
         number_format = "%.2f",
         fontsize_number = 8,
         main = "Sample Correlation (Pearson)")
dev.off()

# Hierarchical clustering dendrogram
sample_dist <- dist(t(count_matrix))
sample_hclust <- hclust(sample_dist, method = "complete")

pdf("2.1_hierarchical_clustering.pdf", width = 10, height = 6)
plot(sample_hclust,
     main = "Hierarchical Clustering of Samples",
     xlab = "",
     sub = "")
dev.off()

# Distance matrix heatmap
dist_matrix <- as.matrix(sample_dist)
pdf("2.1_sample_distance_heatmap.pdf", width = 10, height = 9)
pheatmap(dist_matrix,
         clustering_distance_rows = sample_dist,
         clustering_distance_cols = sample_dist,
         main = "Sample Distance Matrix")
dev.off()


# 2.2 Batch Effect Assessment (if batch info exists)
# ----------------------------------------------------------------------------
# batch_info가 있다면 사용
# batch_info <- read.delim("batch_info.txt")
# 
# # RLE plot
# library(EDASeq)
# rle_data <- betweenLaneNormalization(count_matrix, which = "median")
# pdf("2.2_RLE_plot.pdf", width = 12, height = 6)
# plotRLE(rle_data, col = batch_info$color, outline = FALSE)
# dev.off()


# 2.3 Gene Expression Filtering Analysis
# ----------------------------------------------------------------------------

# 필터링 기준별 유전자 수
filter_thresholds <- c(0, 1, 5, 10, 20, 50, 100)
filter_summary <- data.frame(
  Threshold = filter_thresholds,
  Genes_Remaining = sapply(filter_thresholds, function(x) {
    sum(rowSums(count_matrix > x) >= ncol(count_matrix) * 0.1)  # 10% 샘플에서 발현
  })
)

pdf("2.3_filtering_threshold.pdf", width = 8, height = 6)
plot(filter_summary$Threshold, filter_summary$Genes_Remaining,
     type = "b",
     pch = 19,
     col = "darkblue",
     xlab = "Count Threshold",
     ylab = "Number of Genes Remaining",
     main = "Gene Filtering Impact\n(expressed in ≥10% samples)")
grid()
dev.off()

write.table(filter_summary, "2.3_filtering_summary.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# 일반적 필터링: count > 10 in at least 10% samples
keep <- rowSums(count_matrix > 10) >= ncol(count_matrix) * 0.1
filtered_counts <- count_matrix[keep, ]

cat("Original genes:", nrow(count_matrix), "\n")
cat("Filtered genes:", nrow(filtered_counts), "\n")
cat("Removed genes:", nrow(count_matrix) - nrow(filtered_counts), "\n")


# 2.4 Pearson Correlation between Samples (after filtering)
# ----------------------------------------------------------------------------
cor_matrix_filtered <- cor(filtered_counts, method = "pearson")

pdf("2.4_correlation_heatmap_filtered.pdf", width = 10, height = 9)
pheatmap(cor_matrix_filtered,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         display_numbers = TRUE,
         number_format = "%.2f",
         fontsize_number = 8,
         main = "Sample Correlation (after gene filtering)")
dev.off()

# Correlation distribution
cor_values <- cor_matrix_filtered[lower.tri(cor_matrix_filtered)]
pdf("2.4_correlation_distribution.pdf", width = 8, height = 6)
hist(cor_values,
     breaks = 30,
     col = "lightblue",
     main = "Distribution of Sample Correlations",
     xlab = "Pearson Correlation",
     ylab = "Frequency")
abline(v = median(cor_values), col = "red", lty = 2, lwd = 2)
legend("topleft", legend = paste("Median =", round(median(cor_values), 3)),
       col = "red", lty = 2)
dev.off()


# 2.5 PCA Analysis
# ----------------------------------------------------------------------------

# PCA on log-transformed filtered data
log_filtered <- log2(filtered_counts + 1)
pca_result <- prcomp(t(log_filtered), scale. = TRUE)

# Variance explained
variance_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100

# Scree plot
pdf("2.5_scree_plot.pdf", width = 8, height = 6)
plot(1:length(variance_explained), variance_explained,
     type = "b",
     pch = 19,
     col = "darkblue",
     xlab = "Principal Component",
     ylab = "Variance Explained (%)",
     main = "PCA Scree Plot")
grid()
dev.off()

# PC1 vs PC2 plot
pca_data <- data.frame(
  Sample = sample_names,
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2]
)

p <- ggplot(pca_data, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point(size = 4, color = "steelblue") +
  geom_text(vjust = -1, hjust = 0.5, size = 3) +
  theme_bw() +
  labs(title = "PCA Plot (PC1 vs PC2)",
       x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 2), "%)"))
ggsave("2.5_PCA_plot.pdf", p, width = 10, height = 8)

# Outlier detection based on PCA scores
# Mahalanobis distance
pca_scores <- pca_result$x[, 1:2]
center <- colMeans(pca_scores)
cov_matrix <- cov(pca_scores)
mahal_dist <- mahalanobis(pca_scores, center, cov_matrix)

# Chi-square cutoff (95% confidence)
cutoff <- qchisq(0.95, df = 2)
outliers <- mahal_dist > cutoff

outlier_summary <- data.frame(
  Sample = sample_names,
  Mahalanobis_Distance = mahal_dist,
  Is_Outlier = outliers
)
write.table(outlier_summary, "2.5_PCA_outlier_detection.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# PCA plot with outliers highlighted
pca_data$Outlier <- outliers
p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Outlier, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, hjust = 0.5, size = 3) +
  scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "red")) +
  theme_bw() +
  labs(title = "PCA Plot with Outlier Detection",
       x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 2), "%)"))
ggsave("2.5_PCA_outliers.pdf", p, width = 10, height = 8)

# PC1 vs PC3, PC2 vs PC3
pca_data$PC3 <- pca_result$x[, 3]

p1 <- ggplot(pca_data, aes(x = PC1, y = PC3, label = Sample)) +
  geom_point(size = 4, color = "steelblue") +
  geom_text(vjust = -1, hjust = 0.5, size = 3) +
  theme_bw() +
  labs(title = "PCA Plot (PC1 vs PC3)",
       x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
       y = paste0("PC3 (", round(variance_explained[3], 2), "%)"))

p2 <- ggplot(pca_data, aes(x = PC2, y = PC3, label = Sample)) +
  geom_point(size = 4, color = "steelblue") +
  geom_text(vjust = -1, hjust = 0.5, size = 3) +
  theme_bw() +
  labs(title = "PCA Plot (PC2 vs PC3)",
       x = paste0("PC2 (", round(variance_explained[2], 2), "%)"),
       y = paste0("PC3 (", round(variance_explained[3], 2), "%)"))

library(gridExtra)
pdf("2.5_PCA_additional_plots.pdf", width = 14, height = 6)
grid.arrange(p1, p2, ncol = 2)
dev.off()


# ============================================================================
# Deliverables Summary
# ============================================================================

cat("\n=== QC Analysis Complete ===\n")
cat("Generated files:\n")
cat("1.1: Library size plots and summary\n")
cat("1.2: Gene detection rate plots and summary\n")
cat("1.3: Expression distribution plots (boxplot, density, violin)\n")
cat("2.1: Sample correlation and clustering plots\n")
cat("2.2: (Batch effect - if applicable)\n")
cat("2.3: Gene filtering analysis\n")
cat("2.4: Sample correlation after filtering\n")
cat("2.5: PCA plots and outlier detection\n")
