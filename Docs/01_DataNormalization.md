# 01_DataNormalization.R - Proteomics Data Normalization and Quality Assessment

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

### 1. Data Validation
```r
# Check if input data exists
if (!exists("mat_clean")) {
  stop("Input data 'mat_clean' not found. Please run 00_PreprocessData.R first.")
}

# Validate input data structure
validate_normalization_input(mat_clean)
```

### 2. Variance Stabilizing Normalization (VSN)
```r
# Apply VSN for intensity-dependent variance stabilization
vsn_fit <- vsn2(mat_clean)
mat_norm <- predict(vsn_fit, mat_clean)

# Analyze VSN effect
analyze_normalization_effect(mat_clean, mat_norm, "VSN")
```

### 3. Scaling and Centering
```r
# Apply Z-score scaling and centering
mat_scaled <- t(scale(t(mat_norm), 
                      center = norm_config$center_data, 
                      scale = norm_config$scale_data))

# Analyze scaling effect
analyze_normalization_effect(mat_norm, mat_scaled, "Z-score scaling")
```

### 4. Outlier Detection
```r
# Perform PCA for outlier detection
pca_tmp <- prcomp(t(mat_scaled), center = FALSE, scale. = FALSE)
pc_scores <- pca_tmp$x[, 1:norm_config$pca_components, drop = FALSE]

# Calculate Mahalanobis distance
dist_z <- rowSums(scale(pc_scores)^2)
outliers <- which(dist_z > norm_config$outlier_threshold)

# Remove outliers if detected
if (length(outliers) > 0) {
  outlier_samples <- rownames(pc_scores)[outliers]
  mat_scaled <- mat_scaled[, -outliers, drop = FALSE]
  sample_anno <- sample_anno[-match(outlier_samples, rownames(sample_anno)), , drop = FALSE]
}
```

### 5. Quality Assessment Visualizations
```r
# PCA for quality assessment
pca_res <- prcomp(t(mat_scaled), center = TRUE, scale. = FALSE)
pca_plot <- fviz_pca_ind(pca_res,
                          geom.ind = "point",
                          col.ind = sample_anno$Group,
                          addEllipses = norm_config$pca_ellipses,
                          legend.title = "Group",
                          title = "PCA: Sample Quality Assessment")

# Heatmap of top variable proteins
top_proteins <- names(sort(rowVars(mat_scaled), decreasing = TRUE))[1:norm_config$top_variable_proteins]
heatmap_plot <- pheatmap(mat_scaled[top_proteins, ],
                         annotation_col = sample_anno["Group", drop = FALSE],
                         scale = "row",
                         show_rownames = TRUE,
                         border_color = "white",
                         cellwidth = norm_config$heatmap_cell_size,
                         cellheight = norm_config$heatmap_cell_size,
                         show_colnames = FALSE,
                         main = sprintf("Top %d Variable Proteins", norm_config$top_variable_proteins))
```

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

## Example Output

```
[14:31:00] NORM_INFO: Starting proteomics data normalization pipeline
[14:31:00] NORM_INFO: Input validation passed: 1950 proteins × 45 samples
[14:31:00] NORM_INFO: === VARIANCE STABILIZING NORMALIZATION ===
[14:31:01] NORM_INFO: Applying VSN for intensity-dependent variance stabilization
[14:31:01] NORM_INFO: VSN - Before: mean=12.456, sd=2.345, var=5.498
[14:31:01] NORM_INFO: VSN - After:  mean=12.123, sd=1.234, var=1.523
[14:31:01] NORM_INFO: VSN completed: 1950 proteins × 45 samples
[14:31:01] NORM_INFO: === SCALING AND CENTERING ===
[14:31:01] NORM_INFO: Applying Z-score scaling and centering
[14:31:01] NORM_INFO: Z-score scaling - Before: mean=12.123, sd=1.234, var=1.523
[14:31:01] NORM_INFO: Z-score scaling - After:  mean=0.000, sd=1.000, var=1.000
[14:31:01] NORM_INFO: Scaling completed: 1950 proteins × 45 samples
[14:31:01] NORM_INFO: === OUTLIER DETECTION ===
[14:31:01] NORM_INFO: Performing PCA-based outlier detection
[14:31:01] NORM_INFO: Outlier detection: 45 samples analyzed, threshold = 9.0
[14:31:01] NORM_INFO: No outliers detected; all samples retained
[14:31:02] NORM_INFO: === QUALITY ASSESSMENT VISUALIZATIONS ===
[14:31:02] NORM_INFO: Generating PCA plot for quality assessment
[14:31:02] NORM_INFO: Generating heatmap of top variable proteins
[14:31:02] NORM_INFO: === NORMALIZATION SUMMARY ===
[14:31:02] NORM_INFO: Final normalized dataset: 1950 proteins × 45 samples
[14:31:02] NORM_INFO: Total variance explained by all PCs: 100.0%
[14:31:02] NORM_INFO: PC1 explains 35.2% of variance
[14:31:02] NORM_INFO: PC2 explains 18.7% of variance
[14:31:02] NORM_INFO: Data quality metrics:
[14:31:02] NORM_INFO: - Mean intensity: 0.000
[14:31:02] NORM_INFO: - Standard deviation: 1.000
[14:31:02] NORM_INFO: - Coefficient of variation: 0.000
[14:31:02] NORM_INFO: Normalization pipeline completed successfully
```

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

## Error Handling

### Input Validation
- Checks for required input objects (`mat_clean`, `sample_anno`)
- Validates data matrix structure and dimensions
- Ensures no missing values in input data

### VSN Convergence
- Handles VSN convergence issues
- Provides alternative normalization methods if needed
- Validates VSN fit quality

### Outlier Detection
- Configurable outlier threshold
- Handles edge cases with few samples
- Provides detailed outlier analysis

### Error Messages
- `"Input data 'mat_clean' not found"`: Run preprocessing pipeline first
- `"VSN failed to converge"`: Check data quality or try alternative normalization
- `"Outlier threshold too strict/lenient"`: Adjust `outlier_threshold` parameter

## Usage

### Basic Usage
```r
# Run the normalization pipeline
source("01_DataNormalization.R")

# Access normalized data
dim(mat_scaled)  # Check dimensions
head(mat_scaled) # View first few rows
```

### Custom Configuration
```r
# Modify configuration before running
norm_config$outlier_threshold <- 12  # More lenient outlier detection
norm_config$top_variable_proteins <- 100  # More proteins in heatmap
norm_config$save_plots <- TRUE  # Save plots to files

# Run pipeline with custom settings
source("01_DataNormalization.R")
```

### Advanced Usage
```r
# Access VSN fit for new data
new_data_norm <- predict(vsn_fit, new_data_matrix)

# Use PCA results for additional analysis
pca_scores <- pca_res$x
pca_loadings <- pca_res$rotation
```

## Next Steps

After running this normalization pipeline, the scaled data (`mat_scaled`) is ready for:

1. **Differential expression analysis**: Statistical testing for protein abundance changes
2. **Pathway enrichment analysis**: Functional annotation and pathway analysis
3. **Clustering analysis**: Unsupervised learning for pattern discovery
4. **Machine learning applications**: Classification, regression, or feature selection
5. **Network analysis**: Protein-protein interaction networks
6. **Biomarker discovery**: Identification of diagnostic or prognostic markers

## Technical Notes

### VSN Implementation
- **Statistical foundation**: Based on generalized log transformation
- **Intensity-dependent**: Adapts to local variance patterns
- **Robust estimation**: Handles outliers and extreme values

### Outlier Detection Methodology
- **PCA-based**: Uses principal components for dimensionality reduction
- **Mahalanobis distance**: Multivariate distance measure
- **Configurable threshold**: Adjustable sensitivity for outlier detection

### Performance Considerations
- **Memory efficient**: Optimized for large proteomics datasets
- **Fast execution**: Efficient matrix operations
- **Scalable**: Handles datasets with thousands of proteins and samples

### Reproducibility
- **VSN fit object**: Can be applied to new data
- **Detailed logging**: Timestamped progress tracking
- **Configuration tracking**: All parameters documented and accessible

## Troubleshooting

### Common Issues
1. **VSN convergence**: If VSN fails, check data quality or try log transformation
2. **Outlier detection**: Adjust threshold if too many/few samples flagged
3. **Memory issues**: For large datasets, increase R memory limit
4. **Plot generation**: Ensure sufficient memory for visualization

### Performance Tips
1. **Large datasets**: Consider processing in batches
2. **Memory optimization**: Close unnecessary R objects before running
3. **Plot saving**: Use `save_plots = TRUE` for publication-quality figures

### Alternative Approaches
1. **Log transformation**: Use `log2()` if VSN fails
2. **Quantile normalization**: Alternative to VSN for some datasets
3. **Robust scaling**: Use median-based scaling for outlier-resistant normalization

## Citation

If you use this normalization pipeline in your research, please cite:

```
Joseph, K. (2024). Proteomics Data Normalization Pipeline. 
GitHub repository. https://github.com/yourusername/proteomics-preprocessing
```

## Contact

For questions or issues with this normalization pipeline, please open an issue on GitHub or contact the maintainer. 