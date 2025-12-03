# 필요한 패키지
library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)

# ============================================================================
# 데이터 준비
# ============================================================================

# Count matrix 불러오기
count_data <- read.delim("merged_counts.txt", stringsAsFactors = FALSE)

count_matrix <- count_data %>%
  select(-Geneid, -Chr, -Start, -End, -Strand, -Length) %>%
  as.matrix()

rownames(count_matrix) <- count_data$Geneid

# Gene filtering
keep <- rowSums(count_matrix > 10) >= ncol(count_matrix) * 0.1
filtered_counts <- count_matrix[keep, ]

sample_names <- colnames(filtered_counts)

# 메타데이터 추출
extract_condition <- function(sample_name) {
  condition <- gsub("_[0-9]+$", "", sample_name)
  return(condition)
}

metadata <- data.frame(
  Sample = sample_names,
  Condition = sapply(sample_names, extract_condition),
  row.names = sample_names
)

metadata$Condition <- factor(metadata$Condition)
if ("Control" %in% levels(metadata$Condition)) {
  metadata$Condition <- relevel(metadata$Condition, ref = "Control")
}

cat("\n=== Extracted Metadata ===\n")
print(metadata)

write.table(metadata, "metadata.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# Grid parameters (실제 Fold Change 기준)
param_grid <- expand.grid(
  padj_cutoff = c(0.01, 0.05, 0.1),
  fc_cutoff = c(1.2, 1.5, 2.0, 3.0),  # Fold Change
  stringsAsFactors = FALSE
)

# log2FC cutoff 계산
param_grid$log2fc_cutoff <- log2(param_grid$fc_cutoff)

cat("\n=== Parameter Grid ===\n")
print(param_grid)

dir.create("DEG_results", showWarnings = FALSE)

# ============================================================================
# 1. DESeq2 Setting
# ============================================================================

cat("\n=== 1. DESeq2 Setting ===\n")

# DESeq2 object 생성
dds <- DESeqDataSetFromMatrix(
  countData = filtered_counts,
  colData = metadata,
  design = ~ Condition
)

# DESeq2 실행
dds <- DESeq(dds)

cat("\nAvailable comparisons:\n")
print(resultsNames(dds))

# 비교 조합
conditions <- levels(metadata$Condition)
comparisons <- conditions[conditions != "Control"]

cat("\nComparisons to perform:\n")
for (comp in comparisons) {
  cat("  Control vs", comp, "\n")
}

# ============================================================================
# 2. Differential Expression Analysis
# ============================================================================

cat("\n=== 2. Differential Expression Analysis ===\n")

# 전체 결과 요약
grid_summary <- data.frame()

# 2.1 Grid Search across parameters
cat("\n--- 2.1 Grid Search across parameters ---\n")

for (i in 1:nrow(param_grid)) {
  padj_cut <- param_grid$padj_cutoff[i]
  fc_cut <- param_grid$fc_cutoff[i]
  lfc_cut <- param_grid$log2fc_cutoff[i]
  
  cat(sprintf("\nGrid %d/%d: padj < %.2f, FC > %.1f (log2FC > %.2f)\n", 
              i, nrow(param_grid), padj_cut, fc_cut, lfc_cut))
  
  for (comp in comparisons) {
    cat(sprintf("  Analyzing: Control vs %s\n", comp))
    
    contrast_name <- paste0("Condition_", comp, "_vs_Control")
    
    res <- results(dds, name = contrast_name, alpha = padj_cut)
    
    res_df <- as.data.frame(res)
    res_df$Gene <- rownames(res_df)
    res_df <- res_df[, c("Gene", "baseMean", "log2FoldChange", 
                         "lfcSE", "stat", "pvalue", "padj")]
    
    # 필터되지 않은 전체 결과 저장
    file_prefix_unfiltered <- sprintf("2.1_padj%.2f_fc%.1f_%s_vs_Control", 
                                     padj_cut, fc_cut, comp)
    file_prefix_unfiltered <- gsub("\\.", "p", file_prefix_unfiltered)
    
    res_all <- res_df %>% arrange(pvalue)
    write.table(res_all,
                file.path("DEG_results", paste0(file_prefix_unfiltered, "_unfiltered.txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    # NA 제거 후 분석
    res_df <- res_df[!is.na(res_df$padj), ]
    
    # DEG 분류
    res_df$Regulation <- "Not significant"
    res_df$Regulation[res_df$padj < padj_cut & 
                      res_df$log2FoldChange > lfc_cut] <- "Up"
    res_df$Regulation[res_df$padj < padj_cut & 
                      res_df$log2FoldChange < -lfc_cut] <- "Down"
    
    n_up <- sum(res_df$Regulation == "Up")
    n_down <- sum(res_df$Regulation == "Down")
    n_total <- n_up + n_down
    
    cat(sprintf("    DEGs: %d (Up: %d, Down: %d)\n", n_total, n_up, n_down))
    
    # 파일명
    file_prefix <- sprintf("2.1_padj%.2f_fc%.1f_%s_vs_Control", 
                          padj_cut, fc_cut, comp)
    file_prefix <- gsub("\\.", "p", file_prefix)
    
    # 필터된 전체 결과 저장
    res_sorted <- res_df %>% arrange(padj, desc(abs(log2FoldChange)))
    write.table(res_sorted, 
                file.path("DEG_results", paste0(file_prefix, "_all.txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    # DEG만 저장
    deg_only <- res_sorted %>% filter(Regulation != "Not significant")
    if (nrow(deg_only) > 0) {
      write.table(deg_only,
                  file.path("DEG_results", paste0(file_prefix, "_DEG_only.txt")),
                  sep = "\t", quote = FALSE, row.names = FALSE)
    }
    
    # 요약 정보 저장
    grid_summary <- rbind(grid_summary, data.frame(
      Comparison = paste0("Control_vs_", comp),
      padj_cutoff = padj_cut,
      FC_cutoff = fc_cut,
      log2FC_cutoff = lfc_cut,
      Total_DEGs = n_total,
      Upregulated = n_up,
      Downregulated = n_down,
      Total_Tested = nrow(res_df)
    ))
  }
}

# Grid search summary 저장
write.table(grid_summary, 
            "DEG_results/2.1_grid_search_summary.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# 2.2 Volcano Plot
cat("\n--- 2.2 Volcano Plot ---\n")

for (i in 1:nrow(param_grid)) {
  padj_cut <- param_grid$padj_cutoff[i]
  fc_cut <- param_grid$fc_cutoff[i]
  lfc_cut <- param_grid$log2fc_cutoff[i]
  
  for (comp in comparisons) {
    contrast_name <- paste0("Condition_", comp, "_vs_Control")
    res <- results(dds, name = contrast_name, alpha = padj_cut)
    
    res_df <- as.data.frame(res)
    res_df$Gene <- rownames(res_df)
    res_df <- res_df[!is.na(res_df$padj), ]
    
    res_df$Regulation <- "Not significant"
    res_df$Regulation[res_df$padj < padj_cut & 
                      res_df$log2FoldChange > lfc_cut] <- "Up"
    res_df$Regulation[res_df$padj < padj_cut & 
                      res_df$log2FoldChange < -lfc_cut] <- "Down"
    
    n_up <- sum(res_df$Regulation == "Up")
    n_down <- sum(res_df$Regulation == "Down")
    n_total <- n_up + n_down
    
    file_prefix <- sprintf("2.2_padj%.2f_fc%.1f_%s_vs_Control", 
                          padj_cut, fc_cut, comp)
    file_prefix <- gsub("\\.", "p", file_prefix)
    
    p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), 
                                     color = Regulation)) +
      geom_point(alpha = 0.6, size = 1.5) +
      scale_color_manual(values = c("Up" = "red", 
                                    "Down" = "blue", 
                                    "Not significant" = "grey")) +
      geom_vline(xintercept = c(-lfc_cut, lfc_cut), 
                 linetype = "dashed", color = "black") +
      geom_hline(yintercept = -log10(padj_cut), 
                 linetype = "dashed", color = "black") +
      theme_bw() +
      labs(title = sprintf("Volcano Plot: Control vs %s", comp),
           subtitle = sprintf("padj < %.2f, FC > %.1f (log2FC > %.2f) | DEGs: %d (↑%d, ↓%d)",
                            padj_cut, fc_cut, lfc_cut, n_total, n_up, n_down),
           x = "log2 Fold Change",
           y = "-log10(adjusted p-value)") +
      theme(legend.position = "right")
    
    ggsave(file.path("DEG_results", paste0(file_prefix, "_volcano.pdf")),
           p_volcano, width = 10, height = 8)
  }
}

# 2.3 MA Plot
cat("\n--- 2.3 MA Plot ---\n")

for (i in 1:nrow(param_grid)) {
  padj_cut <- param_grid$padj_cutoff[i]
  fc_cut <- param_grid$fc_cutoff[i]
  lfc_cut <- param_grid$log2fc_cutoff[i]
  
  for (comp in comparisons) {
    contrast_name <- paste0("Condition_", comp, "_vs_Control")
    res <- results(dds, name = contrast_name, alpha = padj_cut)
    
    res_df <- as.data.frame(res)
    res_df$Gene <- rownames(res_df)
    res_df <- res_df[!is.na(res_df$padj), ]
    
    res_df$Regulation <- "Not significant"
    res_df$Regulation[res_df$padj < padj_cut & 
                      res_df$log2FoldChange > lfc_cut] <- "Up"
    res_df$Regulation[res_df$padj < padj_cut & 
                      res_df$log2FoldChange < -lfc_cut] <- "Down"
    
    n_up <- sum(res_df$Regulation == "Up")
    n_down <- sum(res_df$Regulation == "Down")
    n_total <- n_up + n_down
    
    file_prefix <- sprintf("2.3_padj%.2f_fc%.1f_%s_vs_Control", 
                          padj_cut, fc_cut, comp)
    file_prefix <- gsub("\\.", "p", file_prefix)
    
    p_ma <- ggplot(res_df, aes(x = log10(baseMean + 1), y = log2FoldChange,
                                color = Regulation)) +
      geom_point(alpha = 0.6, size = 1.5) +
      scale_color_manual(values = c("Up" = "red",
                                    "Down" = "blue",
                                    "Not significant" = "grey")) +
      geom_hline(yintercept = c(-lfc_cut, lfc_cut),
                 linetype = "dashed", color = "black") +
      geom_hline(yintercept = 0, color = "black") +
      theme_bw() +
      labs(title = sprintf("MA Plot: Control vs %s", comp),
           subtitle = sprintf("padj < %.2f, FC > %.1f (log2FC > %.2f) | DEGs: %d (↑%d, ↓%d)",
                            padj_cut, fc_cut, lfc_cut, n_total, n_up, n_down),
           x = "log10(Base Mean)",
           y = "log2 Fold Change") +
      theme(legend.position = "right")
    
    ggsave(file.path("DEG_results", paste0(file_prefix, "_MA.pdf")),
           p_ma, width = 10, height = 8)
  }
}

# 2.4 Summary Visualization
cat("\n--- 2.4 Summary Visualization ---\n")

# DEG count comparison per comparison
for (comp in comparisons) {
  comp_data <- grid_summary %>%
    filter(grepl(comp, Comparison))
  
  if (nrow(comp_data) > 0) {
    p_summary <- ggplot(comp_data, 
                        aes(x = factor(FC_cutoff), 
                            y = Total_DEGs,
                            fill = factor(padj_cutoff))) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_brewer(palette = "Set2", name = "padj cutoff") +
      theme_bw() +
      labs(title = sprintf("DEG Count Comparison: Control vs %s", comp),
           x = "Fold Change cutoff",
           y = "Number of DEGs") +
      theme(axis.text.x = element_text(angle = 0))
    
    ggsave(file.path("DEG_results", 
                    sprintf("2.4_DEG_count_comparison_%s.pdf", comp)),
           p_summary, width = 10, height = 6)
  }
}

# Heatmap of DEG counts
p_heatmap <- ggplot(grid_summary, 
                    aes(x = factor(FC_cutoff),
                        y = factor(padj_cutoff),
                        fill = Total_DEGs)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Total_DEGs), color = "black", size = 4) +
  scale_fill_gradient(low = "white", high = "darkred", name = "DEG Count") +
  facet_wrap(~Comparison, nrow = 1) +
  theme_bw() +
  labs(title = "DEG Count Heatmap across Parameters",
       x = "Fold Change cutoff",
       y = "padj cutoff") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("DEG_results/2.4_DEG_count_heatmap.pdf",
       p_heatmap, width = 12, height = 6)

# ============================================================================
# 3. Top DEGs Analysis
# ============================================================================

cat("\n=== 3. Top DEGs Analysis ===\n")

# 3.1 Top 50 DEGs (strict criteria: padj < 0.01, FC > 2)
cat("\n--- 3.1 Top 50 DEGs Extraction ---\n")

for (comp in comparisons) {
  contrast_name <- paste0("Condition_", comp, "_vs_Control")
  res <- results(dds, name = contrast_name, alpha = 0.01)
  res_df <- as.data.frame(res)
  res_df$Gene <- rownames(res_df)
  res_df <- res_df[!is.na(res_df$padj), ]
  
  top_deg <- res_df %>%
    filter(padj < 0.01 & abs(log2FoldChange) > 1) %>%  # FC > 2
    arrange(padj) %>%
    head(50)
  
  if (nrow(top_deg) > 0) {
    write.table(top_deg,
                file.path("DEG_results",
                         sprintf("3.1_top50_DEGs_Control_vs_%s.txt", comp)),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    cat(sprintf("  Control vs %s: %d top DEGs\n", comp, nrow(top_deg)))
  }
}

# 3.2 Heatmap of Top DEGs
cat("\n--- 3.2 Heatmap of Top DEGs ---\n")

for (comp in comparisons) {
  contrast_name <- paste0("Condition_", comp, "_vs_Control")
  res <- results(dds, name = contrast_name, alpha = 0.01)
  res_df <- as.data.frame(res)
  res_df$Gene <- rownames(res_df)
  res_df <- res_df[!is.na(res_df$padj), ]
  
  top_deg <- res_df %>%
    filter(padj < 0.01 & abs(log2FoldChange) > 1) %>%
    arrange(padj) %>%
    head(50)
  
  if (nrow(top_deg) >= 10) {
    norm_counts <- counts(dds, normalized = TRUE)
    top_genes <- top_deg$Gene
    top_mat <- norm_counts[top_genes, ]
    
    top_mat_log <- log2(top_mat + 1)
    top_mat_scaled <- t(scale(t(top_mat_log)))
    
    anno_col <- data.frame(
      Condition = metadata$Condition,
      row.names = rownames(metadata)
    )
    
    pdf(file.path("DEG_results",
                 sprintf("3.2_heatmap_top50_Control_vs_%s.pdf", comp)),
        width = 10, height = 12)
    pheatmap(top_mat_scaled,
             annotation_col = anno_col,
             show_rownames = TRUE,
             show_colnames = TRUE,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             fontsize_row = 6,
             main = sprintf("Top 50 DEGs: Control vs %s", comp))
    dev.off()
  }
}

# ============================================================================
# 4. PCA on DEGs
# ============================================================================

cat("\n=== 4. PCA on DEGs ===\n")

# VST transformation
vsd <- vst(dds, blind = FALSE)

for (comp in comparisons) {
  contrast_name <- paste0("Condition_", comp, "_vs_Control")
  res <- results(dds, name = contrast_name, alpha = 0.05)
  res_df <- as.data.frame(res)
  
  deg_genes <- rownames(res_df)[!is.na(res_df$padj) & 
                                 res_df$padj < 0.05 & 
                                 abs(res_df$log2FoldChange) > log2(2)]  # FC > 2
  
  if (length(deg_genes) > 10) {
    vsd_deg <- assay(vsd)[deg_genes, ]
    pca_deg <- prcomp(t(vsd_deg), scale. = TRUE)
    
    variance_explained <- (pca_deg$sdev^2 / sum(pca_deg$sdev^2)) * 100
    
    pca_data <- data.frame(
      Sample = colnames(vsd_deg),
      PC1 = pca_deg$x[, 1],
      PC2 = pca_deg$x[, 2],
      Condition = metadata$Condition
    )
    
    p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, 
                                  color = Condition, label = Sample)) +
      geom_point(size = 4) +
      geom_text_repel(size = 3) +
      theme_bw() +
      labs(title = sprintf("PCA on DEGs: Control vs %s", comp),
           subtitle = sprintf("Using %d DEGs (padj<0.05, FC>2)", 
                            length(deg_genes)),
           x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
           y = paste0("PC2 (", round(variance_explained[2], 2), "%)"))
    
    ggsave(file.path("DEG_results",
                    sprintf("4_PCA_DEGs_Control_vs_%s.pdf", comp)),
           p_pca, width = 10, height = 8)
  }
}

# ============================================================================
# Deliverables
# ============================================================================

cat("\n=== DESeq2 Analysis Complete ===\n")
cat("\nGenerated files:\n")
cat("1. DESeq2 Setting\n")
cat("2. Differential Expression Analysis\n")
cat("  2.1: Grid search results (unfiltered, all, DEG_only) and summary\n")
cat("  2.2: Volcano plots\n")
cat("  2.3: MA plots\n")
cat("  2.4: Summary visualizations\n")
cat("3. Top DEGs Analysis\n")
cat("  3.1: Top 50 DEG lists\n")
cat("  3.2: Heatmaps\n")
cat("4. PCA on DEGs\n")
