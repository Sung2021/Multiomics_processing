# 필요한 패키지
library(dplyr)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(gridExtra)
library(scales)
library(dendextend)

# ============================================================================
# 데이터 준비
# ============================================================================

# Count matrix 준비
count_data <- read.delim("merged_counts.txt", stringsAsFactors = FALSE)

# Count matrix 추출
count_matrix <- count_data %>%
  select(-Geneid, -Chr, -Start, -End, -Strand, -Length) %>%
  as.matrix()

rownames(count_matrix) <- count_data$Geneid
sample_names <- colnames(count_matrix)

# ============================================================================
# 1. Basic Data Quality Checks
# ============================================================================

# 1.1 Library Size Analysis
# ----------------------------------------------------------------------------
library_sizes <- colSums(count_matrix)

lib_df <- data.frame(
  Sample = factor(sample_names, levels = sample_names),
  Library_Size = library_sizes,
  Million_Reads = library_sizes / 1e6
)

# Barplot
p1 <- ggplot(lib_df, aes(x = Sample, y = Million_Reads)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = median(lib_df$Million_Reads), 
             color = "red", linetype = "dashed", linewidth = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Library Size per Sample",
       x = "Sample",
       y = "Library Size (millions)") +
  annotate("text", x = length(sample_names) * 0.9, 
           y = median(lib_df$Million_Reads) * 1.1,
           label = "Median", color = "red")

ggsave("1.1_library_size.pdf", p1, width = 10, height = 6)

# Log scale barplot
p2 <- ggplot(lib_df, aes(x = Sample, y = Library_Size)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  scale_y_log10(labels = comma) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Library Size per Sample (log scale)",
       x = "Sample",
       y = "Library Size (log10)")

ggsave("1.1_library_size_log.pdf", p2, width = 10, height = 6)

# Summary
write.table(lib_df, "1.1_library_size_summary.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


# 1.2 Gene Detection Rate
# ----------------------------------------------------------------------------
detected_genes_0 <- colSums(count_matrix > 0)
detected_genes_10 <- colSums(count_matrix > 10)

detect_df <- data.frame(
  Sample = rep(factor(sample_names, levels = sample_names), 2),
  Threshold = rep(c("count > 0", "count > 10"), each = length(sample_names)),
  Detected = c(detected_genes_0, detected_genes_10)
)

p3 <- ggplot(detect_df, aes(x = Sample, y = Detected, fill = Threshold)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(data = detect_df %>% group_by(Threshold) %>% 
               summarise(median = median(Detected)),
             aes(yintercept = median, color = Threshold),
             linetype = "dashed", linewidth = 1) +
  scale_fill_manual(values = c("darkgreen", "forestgreen")) +
  scale_color_manual(values = c("darkgreen", "forestgreen"), guide = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Detection Rate",
       x = "Sample",
       y = "Number of Genes Detected")

ggsave("1.2_gene_detection_rate.pdf", p3, width = 12, height = 6)

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
# Log transformation
log_counts <- log2(count_matrix + 1)

# Boxplot
log_df <- melt(log_counts)
colnames(log_df) <- c("Gene", "Sample", "Expression")

p4 <- ggplot(log_df, aes(x = Sample, y = Expression)) +
  geom_boxplot(fill = "lightblue", outlier.shape = NA) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Expression Distribution per Sample",
       x = "Sample",
       y = "log2(counts + 1)")

ggsave("1.3_expression_distribution_boxplot.pdf", p4, width = 12, height = 6)

# Density plot
p5 <- ggplot(log_df, aes(x = Expression, color = Sample)) +
  geom_density(linewidth = 1) +
  theme_bw() +
  theme(legend.position = "right") +
  labs(title = "Expression Distribution Density",
       x = "log2(counts + 1)",
       y = "Density")

ggsave("1.3_expression_distribution_density.pdf", p5, width = 10, height = 6)

# Violin plot
p6 <- ggplot(log_df, aes(x = Sample, y = Expression, fill = Sample)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "Expression Distribution (Violin Plot)",
       x = "Sample",
       y = "log2(counts + 1)")

ggsave("1.3_expression_distribution_violin.pdf", p6, width = 12, height = 6)


# ============================================================================
# 2. Further Quality Control Metrics
# ============================================================================

# 2.1 Outlier Detection by Distance
# ----------------------------------------------------------------------------

# Sample correlation matrix (log-transformed data 사용)
cor_matrix <- cor(log_counts, method = "pearson")

# Correlation heatmap
pdf("2.1_sample_correlation_heatmap.pdf", width = 10, height = 9)
pheatmap(cor_matrix,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         display_numbers = TRUE,
         number_format = "%.2f",
         fontsize_number = 8,
         main = "Sample Correlation (Pearson, log2 counts)")
dev.off()

# Correlation heatmap (ggplot version)
cor_df <- melt(cor_matrix)
colnames(cor_df) <- c("Sample1", "Sample2", "Correlation")

p7 <- ggplot(cor_df, aes(x = Sample1, y = Sample2, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", Correlation)), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0.5, limits = c(0, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Sample Correlation (Pearson, log2 counts)",
       x = "", y = "") +
  coord_fixed()

ggsave("2.1_sample_correlation_heatmap_ggplot.pdf", p7, width = 10, height = 9)

# Hierarchical clustering dendrogram (log-transformed data 사용)
sample_dist <- dist(t(log_counts))
sample_hclust <- hclust(sample_dist, method = "complete")

pdf("2.1_hierarchical_clustering.pdf", width = 10, height = 6)
plot(sample_hclust,
     main = "Hierarchical Clustering of Samples (log2 counts)",
     xlab = "", sub = "")
dev.off()

# Distance matrix heatmap
dist_matrix <- as.matrix(sample_dist)
dist_df <- melt(dist_matrix)
colnames(dist_df) <- c("Sample1", "Sample2", "Distance")

p8 <- ggplot(dist_df, aes(x = Sample1, y = Sample2, fill = Distance)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Sample Distance Matrix (log2 counts)",
       x = "", y = "") +
  coord_fixed()

ggsave("2.1_sample_distance_heatmap.pdf", p8, width = 10, height = 9)


# Distance-based outlier detection using hclust
# ----------------------------------------------------------------------------

# 각 샘플의 평균 거리 계산
mean_distances <- rowMeans(dist_matrix)

# 이상치 기준: 평균 + 1.96 * SD (95% threshold)
threshold_95 <- mean(mean_distances) + 1.96 * sd(mean_distances)

# 이상치 판별
outliers_dist <- mean_distances > threshold_95

distance_outlier_summary <- data.frame(
  Sample = sample_names,
  Mean_Distance = mean_distances,
  Z_Score = (mean_distances - mean(mean_distances)) / sd(mean_distances),
  Is_Outlier_95 = outliers_dist,
  Threshold_95 = threshold_95
)

# 정렬해서 저장
distance_outlier_summary <- distance_outlier_summary %>%
  arrange(desc(Mean_Distance))

write.table(distance_outlier_summary, "2.1_distance_outlier_detection.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# 이상치 샘플만 추출
outlier_samples_dist <- distance_outlier_summary %>%
  filter(Is_Outlier_95 == TRUE)

if (nrow(outlier_samples_dist) > 0) {
  cat("\n=== Distance-based Outliers (95% threshold) ===\n")
  print(outlier_samples_dist)
  
  write.table(outlier_samples_dist, "2.1_distance_outliers_list.txt",
              sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  cat("\n=== No distance-based outliers detected ===\n")
}

# 시각화: Mean distance barplot with threshold
dist_plot_df <- distance_outlier_summary
dist_plot_df$Sample <- factor(dist_plot_df$Sample, 
                               levels = dist_plot_df$Sample)

p_dist1 <- ggplot(dist_plot_df, aes(x = Sample, y = Mean_Distance, 
                                     fill = Is_Outlier_95)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = threshold_95, 
             color = "red", linetype = "dashed", linewidth = 1) +
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "red"),
                    name = "Outlier") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Mean Distance from Other Samples",
       subtitle = paste("Threshold (95%):", round(threshold_95, 2)),
       x = "Sample",
       y = "Mean Distance")

ggsave("2.1_distance_outlier_barplot.pdf", p_dist1, width = 12, height = 6)

# Z-score plot
p_dist2 <- ggplot(dist_plot_df, aes(x = Sample, y = Z_Score, 
                                     fill = Is_Outlier_95)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 1.96, 
             color = "red", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = -1.96, 
             color = "red", linetype = "dashed", linewidth = 1) +
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "red"),
                    name = "Outlier") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Distance Z-Score per Sample",
       subtitle = "Threshold: ±1.96 (95% confidence)",
       x = "Sample",
       y = "Z-Score")

ggsave("2.1_distance_zscore_plot.pdf", p_dist2, width = 12, height = 6)

# Dendrogram with outliers highlighted
pdf("2.1_hierarchical_clustering_outliers.pdf", width = 12, height = 6)
dend <- as.dendrogram(sample_hclust)
labels_colors(dend) <- ifelse(outliers_dist[order.dendrogram(dend)], "red", "black")
plot(dend,
     main = "Hierarchical Clustering with Distance-based Outliers (Red = Outlier)")
dev.off()

# Combined plot: distance heatmap with dendrogram and annotation
pdf("2.1_distance_heatmap_with_dendrogram.pdf", width = 12, height = 10)
pheatmap(dist_matrix,
         clustering_distance_rows = sample_dist,
         clustering_distance_cols = sample_dist,
         annotation_row = data.frame(
           Outlier = ifelse(outliers_dist, "Yes", "No"),
           row.names = sample_names
         ),
         annotation_colors = list(Outlier = c("Yes" = "red", "No" = "steelblue")),
         main = "Sample Distance Matrix with Outliers (log2 counts)")
dev.off()


# 2.3 Gene Expression Filtering Analysis
# ----------------------------------------------------------------------------

# 필터링 기준별 유전자 수 (raw counts 사용)
filter_thresholds <- c(0, 1, 5, 10, 20, 50, 100)
filter_summary <- data.frame(
  Threshold = filter_thresholds,
  Genes_Remaining = sapply(filter_thresholds, function(x) {
    sum(rowSums(count_matrix > x) >= ncol(count_matrix) * 0.1)
  })
)

p9 <- ggplot(filter_summary, aes(x = Threshold, y = Genes_Remaining)) +
  geom_line(color = "darkblue", linewidth = 1.5) +
  geom_point(color = "darkblue", size = 3) +
  theme_bw() +
  labs(title = "Gene Filtering Impact\n(expressed in ≥10% samples)",
       x = "Count Threshold",
       y = "Number of Genes Remaining")

ggsave("2.3_filtering_threshold.pdf", p9, width = 8, height = 6)

write.table(filter_summary, "2.3_filtering_summary.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# 일반적 필터링
keep <- rowSums(count_matrix > 10) >= ncol(count_matrix) * 0.1
filtered_counts <- count_matrix[keep, ]

cat("\nOriginal genes:", nrow(count_matrix), "\n")
cat("Filtered genes:", nrow(filtered_counts), "\n")
cat("Removed genes:", nrow(count_matrix) - nrow(filtered_counts), "\n")

# Log transformation of filtered data
log_filtered <- log2(filtered_counts + 1)


# 2.4 Pearson Correlation between Samples (after filtering)
# ----------------------------------------------------------------------------
cor_matrix_filtered <- cor(log_filtered, method = "pearson")

# Correlation heatmap
cor_filt_df <- melt(cor_matrix_filtered)
colnames(cor_filt_df) <- c("Sample1", "Sample2", "Correlation")

p10 <- ggplot(cor_filt_df, aes(x = Sample1, y = Sample2, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", Correlation)), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0.5, limits = c(0, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Sample Correlation (after gene filtering, log2 counts)",
       x = "", y = "") +
  coord_fixed()

ggsave("2.4_correlation_heatmap_filtered.pdf", p10, width = 10, height = 9)

# Correlation distribution
cor_values <- cor_matrix_filtered[lower.tri(cor_matrix_filtered)]
cor_dist_df <- data.frame(Correlation = cor_values)

p11 <- ggplot(cor_dist_df, aes(x = Correlation)) +
  geom_histogram(bins = 30, fill = "lightblue", color = "black") +
  geom_vline(xintercept = median(cor_values), 
             color = "red", linetype = "dashed", linewidth = 1) +
  theme_bw() +
  labs(title = "Distribution of Sample Correlations",
       x = "Pearson Correlation",
       y = "Frequency") +
  annotate("text", x = median(cor_values), y = Inf, 
           label = paste("Median =", round(median(cor_values), 3)),
           vjust = 2, color = "red")

ggsave("2.4_correlation_distribution.pdf", p11, width = 8, height = 6)


# 2.5 PCA Analysis
# ----------------------------------------------------------------------------

# PCA on log-transformed filtered data
pca_result <- prcomp(t(log_filtered), scale. = TRUE)

# Variance explained
variance_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100

# Scree plot
scree_df <- data.frame(
  PC = 1:length(variance_explained),
  Variance = variance_explained
)

p12 <- ggplot(scree_df, aes(x = PC, y = Variance)) +
  geom_line(color = "darkblue", linewidth = 1.5) +
  geom_point(color = "darkblue", size = 3) +
  theme_bw() +
  labs(title = "PCA Scree Plot",
       x = "Principal Component",
       y = "Variance Explained (%)")

ggsave("2.5_scree_plot.pdf", p12, width = 8, height = 6)

# PC1 vs PC2 plot
pca_data <- data.frame(
  Sample = sample_names,
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  PC3 = pca_result$x[, 3]
)

p13 <- ggplot(pca_data, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point(size = 4, color = "steelblue") +
  geom_text(vjust = -1, hjust = 0.5, size = 3) +
  theme_bw() +
  labs(title = "PCA Plot (PC1 vs PC2)",
       x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 2), "%)"))

ggsave("2.5_PCA_plot.pdf", p13, width = 10, height = 8)

# Outlier detection based on PCA scores
pca_scores <- pca_result$x[, 1:2]
center <- colMeans(pca_scores)
cov_matrix <- cov(pca_scores)
mahal_dist <- mahalanobis(pca_scores, center, cov_matrix)

# Chi-square cutoff (95% confidence)
cutoff <- qchisq(0.95, df = 2)
outliers_pca <- mahal_dist > cutoff

outlier_summary <- data.frame(
  Sample = sample_names,
  Mahalanobis_Distance = mahal_dist,
  Is_Outlier = outliers_pca
)
write.table(outlier_summary, "2.5_PCA_outlier_detection.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# 이상치 샘플만 추출
outlier_samples_pca <- outlier_summary %>%
  filter(Is_Outlier == TRUE)

if (nrow(outlier_samples_pca) > 0) {
  cat("\n=== PCA-based Outliers (95% threshold) ===\n")
  print(outlier_samples_pca)
  
  write.table(outlier_samples_pca, "2.5_PCA_outliers_list.txt",
              sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  cat("\n=== No PCA-based outliers detected ===\n")
}

# PCA plot with outliers highlighted
pca_data$Outlier <- outliers_pca

p14 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Outlier, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, hjust = 0.5, size = 3) +
  scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "red")) +
  theme_bw() +
  labs(title = "PCA Plot with Outlier Detection",
       x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 2), "%)"))

ggsave("2.5_PCA_outliers.pdf", p14, width = 10, height = 8)

# PC1 vs PC3
p15 <- ggplot(pca_data, aes(x = PC1, y = PC3, label = Sample)) +
  geom_point(size = 4, color = "steelblue") +
  geom_text(vjust = -1, hjust = 0.5, size = 3) +
  theme_bw() +
  labs(title = "PCA Plot (PC1 vs PC3)",
       x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
       y = paste0("PC3 (", round(variance_explained[3], 2), "%)"))

# PC2 vs PC3
p16 <- ggplot(pca_data, aes(x = PC2, y = PC3, label = Sample)) +
  geom_point(size = 4, color = "steelblue") +
  geom_text(vjust = -1, hjust = 0.5, size = 3) +
  theme_bw() +
  labs(title = "PCA Plot (PC2 vs PC3)",
       x = paste0("PC2 (", round(variance_explained[2], 2), "%)"),
       y = paste0("PC3 (", round(variance_explained[3], 2), "%)"))

pdf("2.5_PCA_additional_plots.pdf", width = 14, height = 6)
grid.arrange(p15, p16, ncol = 2)
dev.off()


# ============================================================================
# Deliverables Summary
# ============================================================================

cat("\n=== QC Analysis Complete ===\n")
cat("Generated files:\n")
cat("1.1: Library size plots and summary (raw counts)\n")
cat("1.2: Gene detection rate plots and summary (raw counts)\n")
cat("1.3: Expression distribution plots (log2 counts)\n")
cat("2.1: Sample correlation, clustering, and distance-based outliers (log2 counts)\n")
cat("2.3: Gene filtering analysis (raw counts)\n")
cat("2.4: Sample correlation after filtering (log2 counts)\n")
cat("2.5: PCA plots and outlier detection (log2 counts)\n")
