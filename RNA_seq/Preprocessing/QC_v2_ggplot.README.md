# RNA-seq QC Pipeline — `02.QC_v2_ggplot.R`

This script performs **comprehensive QC for bulk RNA-seq count matrices** using base R, `dplyr`, `ggplot2`, `pheatmap`, and clustering/distance metrics.
It expects a **raw merged count matrix (`merged_counts.txt`)** containing feature metadata columns followed by sample count columns.

---

## 1. Input

**File:** `merged_counts.txt` (tab-delimited)
Required columns:

* `Geneid`, `Chr`, `Start`, `End`, `Strand`, `Length`
* One or more sample columns with raw counts.

The script constructs:

* `count_matrix`: numeric sample-by-gene matrix
* `log_counts`: log2(counts + 1)
* `filtered_counts`: post-filtering gene matrix

---

## 2. QC Modules

### **1. Library Size Analysis**

Outputs:

* Barplot of library sizes (raw + log10)
* Summary table with total reads and reads in millions
  Checks for:
* Library size imbalance
* Median-based deviation

Files:

* `1.1_library_size.pdf`
* `1.1_library_size_log.pdf`
* `1.1_library_size_summary.txt`

---

### **2. Gene Detection Rate**

Counts number of expressed genes per sample under two thresholds:

* count > 0
* count > 10

Outputs:

* Barplot with median reference lines
* Detection summary table

Files:

* `1.2_gene_detection_rate.pdf`
* `1.2_gene_detection_summary.txt`

---

### **3. Expression Distribution Analysis**

Using log2(count+1):

* Boxplots
* Density curves
* Violin plots

Files:

* `1.3_expression_distribution_boxplot.pdf`
* `1.3_expression_distribution_density.pdf`
* `1.3_expression_distribution_violin.pdf`

---

## 3. Additional QC Metrics

### **4. Sample Correlation, Clustering, Distance-based Outliers**

Calculates:

* Pearson correlation matrix
* Hierarchical clustering
* Distance matrix (Euclidean on log-counts)
* Outlier detection via:

  * Mean distance > mean + 1.96*SD
  * Z-score threshold ±1.96

Outputs:

* Correlation heatmaps (pheatmap + ggplot)
* Clustering dendrogram
* Distance heatmaps
* Outlier summary + plots

Files:

* `2.1_sample_correlation_heatmap.pdf`
* `2.1_sample_correlation_heatmap_ggplot.pdf`
* `2.1_hierarchical_clustering.pdf`
* `2.1_sample_distance_heatmap.pdf`
* `2.1_distance_outlier_detection.txt`
* `2.1_distance_outliers_list.txt` (if any)
* `2.1_distance_outlier_barplot.pdf`
* `2.1_distance_zscore_plot.pdf`
* `2.1_hierarchical_clustering_outliers.pdf`
* `2.1_distance_heatmap_with_dendrogram.pdf`

---

## 4. Gene Filtering Analysis

### Filtering rule:

A gene is kept if:

* expressed above threshold **in ≥10% of samples**

Thresholds tested: 0,1,5,10,20,50,100.

Outputs:

* Line plot showing remaining genes vs threshold
* Summary table
* Filtered matrix summary printed to console

Files:

* `2.3_filtering_threshold.pdf`
* `2.3_filtering_summary.txt`

---

## 5. Correlation After Filtering

Recomputes Pearson correlation using filtered genes only.

Outputs:

* Correlation heatmap
* Correlation distribution histogram

Files:

* `2.4_correlation_heatmap_filtered.pdf`
* `2.4_correlation_distribution.pdf`

---

## 6. PCA + Outlier Detection

Performs PCA on log-transformed filtered counts.
Outlier detection uses Mahalanobis distance vs χ² cutoff (95%).

Outputs:

* Scree plot
* PCA plots (PC1–PC2, PC1–PC3, PC2–PC3)
* PCA outlier summary + highlight plot

Files:

* `2.5_scree_plot.pdf`
* `2.5_PCA_plot.pdf`
* `2.5_PCA_outlier_detection.txt`
* `2.5_PCA_outliers_list.txt`
* `2.5_PCA_outliers.pdf`
* `2.5_PCA_additional_plots.pdf`

---

## 7. Summary of All Deliverables

The script generates QC figures and tables covering:

1. Library size
2. Gene detection
3. Expression distribution
4. Correlation & clustering
5. Distance-based outliers
6. Gene filtering impact
7. PCA + multivariate outliers

These outputs provide a complete EDA/QC package for bulk RNA-seq before downstream differential expression or normalization workflows.

---
