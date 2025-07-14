# 01_DataNormalization.R

## Overview

This script implements an advanced normalization and quality assessment pipeline for proteomics data. It applies Variance Stabilizing Normalization (VSN) to handle intensity-dependent variance, performs Z-score scaling for downstream analysis, and implements PCA-based outlier detection with comprehensive quality assessment visualizations.

## Purpose

The normalization pipeline addresses critical challenges in proteomics data analysis:
- **Variance stabilization**: Handles heteroscedasticity in proteomics intensity data
- **Data scaling**: Prepares data for statistical analysis and machine learning
- **Outlier detection**: Identifies and removes problematic samples
- **Quality assessment**: Comprehensive evaluation of data quality and structure

## Technical Approach

### Variance Stabilizing Normalization (VSN)
- **Method**: Uses `vsn2()` function for intensity-dependent variance stabilization
- **Rationale**: Proteomics data often exhibits heteroscedasticity (variance depends on intensity)
- **Benefits**: Stabilizes variance across intensity range, improves statistical power

### Z-score Scaling and Centering
- **Method**: Applies `scale()` function with centering and scaling
- **Rationale**: Standardizes data for downstream statistical analysis
- **Benefits**: Ensures all features contribute equally to analysis

### PCA-based Outlier Detection
- **Method**: Uses Mahalanobis distance in principal component space
- **Threshold**: Configurable threshold (default: 9, equivalent to >3 SD)
- **Rationale**: Identifies samples that deviate significantly from the majority

## Dependencies

```r
# Required R packages
library(vsn)         # Variance stabilizing normalization
library(matrixStats) # Matrix statistics for variance calculation
library(factoextra)  # PCA visualization
library(pheatmap)    # Heatmap generation
library(tidyverse)   # Data manipulation and visualization
```

## Input Requirements

### Prerequisites
- **Required input**: `mat_clean` (from `00_PreprocessData.R`)
- **Required input**: `sample_anno` (sample annotations)
- **Data format**: Numeric matrix with proteins as rows, samples as columns

### Data Quality Requirements
- **No missing values**: Input matrix should be complete (no NA values)
- **Numeric data**: All values should be numeric intensity measurements
- **Valid dimensions**: Matrix should have >1 protein and >1 sample

## Configuration Parameters

```r
norm_config <- list(
  # Normalization settings
  vsn_method = "vsn2",           # VSN method to use
  
  # Scaling parameters
  center_data = TRUE,             # Center the data
  scale_data = TRUE,              # Scale the data
  
  # Outlier detection
  outlier_threshold = 9,          # Mahalanobis distance threshold (>3 SD)
  pca_components = 2,             # Number of PCs for outlier detection
  
  # Visualization settings
  top_variable_proteins = 50,     # Number of proteins for heatmap
  heatmap_cell_size = 8,         # Cell size for heatmap
  pca_ellipses = TRUE,           # Add confidence ellipses to PCA
  
  # Output settings
  verbose = TRUE,
  save_plots = FALSE,
  plot_format = "pdf"
)
```

## Pipeline Steps
 - Data Validation
 - Variance Stabilizing Normalization (VSN)
 - Scaling and Centering
 - Outlier Detection
 - Quality Assessment Visualizations

## Output Objects

### Primary Outputs
- `mat_scaled`: Normalized and scaled intensity matrix (proteins × samples)
- `mat_norm`: VSN-normalized matrix (before scaling)
- `vsn_fit`: VSN fit object for future predictions

### Quality Assessment Objects
- `pca_res`: PCA results for quality assessment
- `sample_anno`: Updated sample annotations (outliers removed)

### Configuration and Metadata
- `norm_config`: Configuration parameters used

## Quality Metrics

### Normalization Statistics
- **Before/after comparisons**: Mean, standard deviation, and variance for each step
- **Variance stabilization**: Reduction in heteroscedasticity
- **Scaling effectiveness**: Achievement of zero mean and unit variance

### PCA Quality Assessment
- **Variance explained**: Percentage of variance explained by each principal component
- **Sample clustering**: Group separation in PCA space
- **Outlier detection**: Samples with high Mahalanobis distance

### Data Quality Metrics
- **Mean intensity**: Should be approximately 0 after scaling
- **Standard deviation**: Should be approximately 1 after scaling
- **Coefficient of variation**: Measure of relative variability
