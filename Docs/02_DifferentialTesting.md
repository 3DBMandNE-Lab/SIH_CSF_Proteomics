# 02_DifferentialTesting.R 
### Proteomics Differential Abundance Analysis Pipeline

## Overview

This script implements an advanced differential abundance analysis pipeline for proteomics data using the limma framework. It performs empirical Bayes moderated t-tests to identify proteins with significant abundance changes between experimental groups, with comprehensive statistical analysis, multiple testing correction, and publication-ready visualizations.

## Purpose

The differential abundance analysis pipeline addresses critical challenges in proteomics data analysis:
- **Statistical rigor**: Implements empirical Bayes moderated t-tests for robust statistical inference
- **Multiple testing correction**: Applies Benjamini-Hochberg FDR correction
- **Visualization**: Generates publication-ready volcano plots and heatmaps
- **Quality assessment**: Provides comprehensive statistical summaries and quality metrics
- **Flexible design**: Handles various experimental designs and group comparisons

## Technical Approach

### Empirical Bayes Moderated t-tests (limma)
- **Method**: Uses `lmFit()` and `eBayes()` functions from limma package
- **Rationale**: Borrows information across proteins to improve variance estimation
- **Benefits**: More powerful than standard t-tests, especially for small sample sizes

### Multiple Testing Correction
- **Method**: Benjamini-Hochberg (BH) false discovery rate correction
- **Rationale**: Controls false positives in high-dimensional data
- **Threshold**: Configurable FDR threshold (default: 0.05)

### Statistical Thresholds
- **FDR threshold**: Controls false discovery rate (default: 0.05)
- **Log fold-change threshold**: Minimum effect size (default: 1.0)
- **Combined criteria**: Proteins must meet both significance and effect size thresholds

## Dependencies

```r
# Required R packages
library(limma)           # Linear models for microarray data (empirical Bayes)
library(EnhancedVolcano) # Enhanced volcano plots
library(ggplot2)         # Advanced plotting
library(ggrepel)         # Label positioning for plots
library(dplyr)           # Data manipulation
library(pheatmap)        # Heatmap generation
```

## Input Requirements

### Prerequisites
- **Required input**: `mat_scaled` (from `01_DataNormalization.R`)
- **Required input**: `sample_anno` (sample annotations with 'Group' column)
- **Data format**: Normalized and scaled numeric matrix with proteins as rows, samples as columns

### Data Quality Requirements
- **No missing values**: Input matrix should be complete (no NA values)
- **Numeric data**: All values should be normalized intensity measurements
- **Group information**: Sample annotations must contain 'Group' column with at least 2 groups
- **Valid dimensions**: Matrix should have >1 protein and >1 sample per group

## Configuration Parameters

```r
diff_config <- list(
  # Analysis settings
  reference_group = "NoLeak",     # Reference group for comparisons
  test_group = "Leak",           # Test group for comparisons
  comparison_name = "Leak_vs_NoLeak",  # Name for output files
  
  # Statistical thresholds
  pvalue_threshold = 0.05,        # P-value threshold
  fdr_threshold = 0.05,          # FDR threshold
  logfc_threshold = 1.0,         # Log2 fold-change threshold
  
  # Visualization settings
  volcano_point_alpha = 0.7,     # Point transparency
  volcano_label_size = 3,        # Label font size
  volcano_max_labels = 20,       # Maximum labels to show
  heatmap_top_n = 20,           # Number of top proteins for heatmap
  heatmap_cell_size = 10,       # Cell size for heatmap
  
  # Output settings
  verbose = TRUE,
  save_results = TRUE,
  save_plots = FALSE,
  plot_format = "pdf"
)
```

## Pipeline Steps

 - Data Validation
 - Experimental Design
 - Differential Abundance Analysis
 - Results Summary

## Output Objects

### Primary Outputs
- `de_res`: Differential abundance results (limma topTable output)
- `fit2`: Fitted limma model with empirical Bayes results

### Visualization Objects
- `volcano_plot`: Enhanced volcano plot with background shading
- `heatmap_plot`: Heatmap of top differential proteins

### Configuration and Metadata
- `diff_config`: Configuration parameters used

## Quality Metrics

### Statistical Summary
- **Total proteins**: Number of proteins analyzed
- **Significant proteins**: Proteins meeting FDR and logFC thresholds
- **Upregulated proteins**: Proteins with positive logFC and significant FDR
- **Downregulated proteins**: Proteins with negative logFC and significant FDR

### Effect Size Distribution
- **Log fold-change range**: Distribution of effect sizes
- **Significance distribution**: Distribution of adjusted p-values
- **Volcano plot quadrants**: Visualization of significance vs effect size

### Quality Assessment
- **FDR control**: False discovery rate control assessment
- **Effect size distribution**: Assessment of meaningful biological changes
- **Sample size adequacy**: Power analysis for detected effects

## Next Steps

After running this differential analysis pipeline, the results (`de_res`) are ready for:

1. **Pathway enrichment analysis**: Functional annotation of significant proteins
2. **Gene ontology analysis**: Biological process and molecular function analysis
3. **Protein-protein interaction networks**: Network analysis of significant proteins
4. **Biomarker validation**: Validation of candidate biomarkers
5. **Functional annotation**: Detailed functional analysis of significant proteins
6. **Meta-analysis**: Integration with other datasets or studies
