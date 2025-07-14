# 02_DifferentialTesting.R - Proteomics Differential Abundance Analysis Pipeline

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

### 1. Data Validation
```r
# Check if input data exists
if (!exists("mat_scaled")) {
  stop("Input data 'mat_scaled' not found. Please run 01_DataNormalization.R first.")
}

# Validate input data structure
validate_differential_input(mat_scaled, sample_anno)
```

### 2. Experimental Design
```r
# Create design matrix
design_result <- create_design_matrix(sample_anno, diff_config$reference_group)
design <- design_result$design
sample_anno <- design_result$sample_annotations
```

### 3. Differential Abundance Analysis
```r
# Fit linear model
fit <- lmFit(mat_scaled, design)

# Apply empirical Bayes shrinkage
fit2 <- eBayes(fit)

# Extract results
de_res <- topTable(
  fit2,
  coef = paste0("Group", diff_config$test_group),
  number = Inf,
  adjust.method = "BH"
)
```

### 4. Results Summary
```r
# Calculate summary statistics
total_proteins <- nrow(de_res)
significant_proteins <- sum(de_res$adj.P.Val < diff_config$fdr_threshold)
upregulated <- sum(de_res$adj.P.Val < diff_config$fdr_threshold & de_res$logFC > diff_config$logfc_threshold)
downregulated <- sum(de_res$adj.P.Val < diff_config$fdr_threshold & de_res$logFC < -diff_config$logfc_threshold)
```

### 5. Visualization
```r
# Generate volcano plot
volcano_plot <- create_volcano_plot(de_res, diff_config)

# Generate heatmap of top differential proteins
top_proteins <- rownames(de_res)[1:diff_config$heatmap_top_n]
heatmap_plot <- pheatmap(mat_scaled[top_proteins, ],
                         annotation_col = sample_anno["Group", drop = FALSE],
                         scale = "row",
                         show_rownames = TRUE)
```

## Output Objects

### Primary Outputs
- `de_res`: Differential abundance results (limma topTable output)
- `fit2`: Fitted limma model with empirical Bayes results

### Visualization Objects
- `volcano_plot`: Enhanced volcano plot with background shading
- `heatmap_plot`: Heatmap of top differential proteins

### Configuration and Metadata
- `diff_config`: Configuration parameters used

## Example Output

```
[14:32:00] DIFF_INFO: Starting proteomics differential abundance analysis
[14:32:00] DIFF_INFO: Input validation passed: 1950 proteins × 45 samples
[14:32:00] DIFF_INFO: === EXPERIMENTAL DESIGN ===
[14:32:00] DIFF_INFO: Design matrix created with 2 groups
[14:32:00] DIFF_INFO: Reference group: NoLeak
[14:32:00] DIFF_INFO: Design matrix columns: (Intercept), GroupLeak
[14:32:00] DIFF_INFO: === DIFFERENTIAL ABUNDANCE ANALYSIS ===
[14:32:00] DIFF_INFO: Fitting linear model with limma
[14:32:00] DIFF_INFO: Applying empirical Bayes shrinkage
[14:32:00] DIFF_INFO: Extracting differential abundance results
[14:32:00] DIFF_INFO: === RESULTS SUMMARY ===
[14:32:00] DIFF_INFO: Total proteins analyzed: 1950
[14:32:00] DIFF_INFO: Significant proteins (FDR < 0.050): 156
[14:32:00] DIFF_INFO: Upregulated proteins: 89
[14:32:00] DIFF_INFO: Downregulated proteins: 67
[14:32:00] DIFF_INFO: Top significant proteins: PROTEIN1, PROTEIN2, PROTEIN3, ...
[14:32:00] DIFF_INFO: === VISUALIZATION ===
[14:32:00] DIFF_INFO: Generating volcano plot
[14:32:00] DIFF_INFO: Generating heatmap of top differential proteins
[14:32:00] DIFF_INFO: === RESULTS EXPORT ===
[14:32:00] DIFF_INFO: Results saved to: Differential_Proteins_Leak_vs_NoLeak.csv
[14:32:00] DIFF_INFO: === DIFFERENTIAL ANALYSIS SUMMARY ===
[14:32:00] DIFF_INFO: Analysis completed for Leak vs NoLeak comparison
[14:32:00] DIFF_INFO: Statistical thresholds: FDR < 0.050, |logFC| > 1.0
[14:32:00] DIFF_INFO: Quality metrics:
[14:32:00] DIFF_INFO: - Total proteins: 1950
[14:32:00] DIFF_INFO: - Significant proteins: 156 (8.0%)
[14:32:00] DIFF_INFO: - Upregulated: 89
[14:32:00] DIFF_INFO: - Downregulated: 67
[14:32:00] DIFF_INFO: Differential analysis pipeline completed successfully
```

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

## Error Handling

### Input Validation
- Checks for required input objects (`mat_scaled`, `sample_anno`)
- Validates data matrix structure and dimensions
- Ensures group information is available and valid

### Statistical Validation
- Validates experimental design matrix
- Checks for sufficient sample sizes per group
- Ensures statistical assumptions are met

### Error Messages
- `"Input data 'mat_scaled' not found"`: Run normalization pipeline first
- `"Sample annotations must contain 'Group' column"`: Check sample annotation format
- `"At least 2 groups required"`: Ensure proper experimental design

## Usage

### Basic Usage
```r
# Run the differential analysis pipeline
source("02_DifferentialTesting.R")

# Access results
head(de_res)  # View top results
dim(de_res)   # Check dimensions
```

### Custom Configuration
```r
# Modify configuration before running
diff_config$fdr_threshold <- 0.01  # More stringent FDR control
diff_config$logfc_threshold <- 0.5  # Lower effect size threshold
diff_config$reference_group <- "Control"  # Change reference group
diff_config$test_group <- "Treatment"  # Change test group

# Run pipeline with custom settings
source("02_DifferentialTesting.R")
```

### Advanced Usage
```r
# Access fitted model for additional analysis
coef_names <- colnames(fit2$coefficients)
fit2$coefficients[1:5, ]  # View coefficients for top proteins

# Extract specific contrasts
specific_results <- topTable(fit2, coef = "GroupLeak", number = 50)

# Access model diagnostics
plotSA(fit2)  # Mean-variance trend
```

## Next Steps

After running this differential analysis pipeline, the results (`de_res`) are ready for:

1. **Pathway enrichment analysis**: Functional annotation of significant proteins
2. **Gene ontology analysis**: Biological process and molecular function analysis
3. **Protein-protein interaction networks**: Network analysis of significant proteins
4. **Biomarker validation**: Validation of candidate biomarkers
5. **Functional annotation**: Detailed functional analysis of significant proteins
6. **Meta-analysis**: Integration with other datasets or studies

## Technical Notes

### Limma Implementation
- **Linear modeling**: Uses linear models for differential expression
- **Empirical Bayes**: Borrows information across features for improved variance estimation
- **Moderated statistics**: Provides more robust statistical inference

### Multiple Testing Correction
- **Benjamini-Hochberg**: Controls false discovery rate
- **Adaptive thresholds**: Can be adjusted based on study requirements
- **Interpretation**: FDR < 0.05 means <5% false positives among significant results

### Performance Considerations
- **Memory efficient**: Optimized for large proteomics datasets
- **Fast execution**: Efficient matrix operations and statistical computations
- **Scalable**: Handles datasets with thousands of proteins and samples

### Reproducibility
- **Detailed logging**: Timestamped progress tracking
- **Configuration tracking**: All parameters documented and accessible
- **Export capabilities**: Results and plots can be saved for publication

## Troubleshooting

### Common Issues
1. **No significant results**: Check statistical thresholds or sample sizes
2. **Convergence issues**: Ensure sufficient sample sizes per group
3. **Memory issues**: For large datasets, increase R memory limit
4. **Design matrix errors**: Verify group assignments and experimental design

### Performance Tips
1. **Large datasets**: Consider processing in batches
2. **Memory optimization**: Close unnecessary R objects before running
3. **Plot saving**: Use `save_plots = TRUE` for publication-quality figures

### Alternative Approaches
1. **Different statistical tests**: Consider t-tests or non-parametric methods
2. **Alternative corrections**: Use Bonferroni or other multiple testing corrections
3. **Effect size measures**: Use different effect size thresholds based on biological relevance

## Citation

If you use this differential analysis pipeline in your research, please cite:

```
Joseph, K. (2024). Proteomics Differential Abundance Analysis Pipeline. 
GitHub repository. https://github.com/yourusername/proteomics-preprocessing
```

## Contact

For questions or issues with this differential analysis pipeline, please open an issue on GitHub or contact the maintainer. 