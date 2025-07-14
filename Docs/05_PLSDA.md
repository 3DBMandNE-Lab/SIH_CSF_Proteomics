# PLS-DA Analysis Pipeline for Proteomics Data Classification

## Overview

The `05_PLSDA.R` script performs Partial Least Squares Discriminant Analysis (PLS-DA) on preprocessed proteomics data to classify samples and identify key discriminatory proteins. This pipeline is designed for supervised classification analysis that combines dimensionality reduction with discriminant analysis to maximize separation between predefined groups.

## Purpose

PLS-DA analysis serves as a powerful supervised classification tool that:
- Classifies samples into predefined groups with high accuracy
- Identifies key discriminatory proteins through VIP scores
- Provides feature selection through sparse PLS-DA
- Offers biological interpretation through GO enrichment analysis
- Evaluates model performance through cross-validation
- Creates comprehensive visualizations for sample classification

## Prerequisites

### Input Data Requirements

The script expects the following data structures from previous pipeline steps:

1. **`mat_scaled`**: A normalized and scaled numeric matrix (proteins × samples) with:
   - No missing values
   - Properly scaled and centered data
   - Row names as protein/gene identifiers
   - Column names as sample identifiers

2. **`sample_anno`**: A data frame with:
   - Row names matching sample names in `mat_scaled`
   - A `Group` column indicating sample conditions
   - At least 3 samples per group for reliable classification

### Dependencies

The script requires the following R packages:
- **Core PLS-DA**: `mixOmics` for PLS-DA analysis
- **Visualization**: `ggplot2` for plotting
- **Data manipulation**: `dplyr` for data processing
- **Biological analysis**: `clusterProfiler`, `org.Hs.eg.db`
- **Additional plotting**: `viridis` for color palettes

## Configuration

### Key Parameters

The script includes a comprehensive configuration system:

```r
config <- list(
  # Input/Output paths
  input_dir = "data/processed",
  output_dir = "results/plsda",
  
  # PLS-DA parameters
  n_components = 2,              # Number of PLS components
  validation_folds = 5,          # Number of folds for cross-validation
  n_repeats = 10,               # Number of repeats for performance evaluation
  test_keepX = c(5, 10, 15, 20), # Test values for sparse PLS-DA
  
  # Feature selection parameters
  top_vip_count = 20,           # Number of top VIP proteins to analyze
  top_loading_count = 20,       # Number of top loading proteins to analyze
  
  # Visualization parameters
  plot_width = 12,
  plot_height = 10,
  dpi = 300,
  
  # Statistical parameters
  pvalue_cutoff = 0.05,
  qvalue_cutoff = 0.05,
  
  # Random seed for reproducibility
  seed = 42
)
```

### Parameter Descriptions

- **`n_components`**: Number of PLS components to extract (typically 2-3)
- **`validation_folds`**: Number of folds for cross-validation
- **`n_repeats`**: Number of repeats for performance evaluation
- **`test_keepX`**: Test values for sparse PLS-DA feature selection
- **`top_vip_count`**: Number of top VIP proteins to analyze
- **`top_loading_count`**: Number of top loading proteins to analyze

## Pipeline Steps

### 1. Data Preparation

The script begins by preparing data for PLS-DA analysis:

- Transposes the protein matrix to samples × features format
- Creates factor variable for group classification
- Validates sample group sizes and data quality
- Ensures sufficient samples per group for reliable classification

### 2. PLS-DA Model Fitting

Performs the core PLS-DA analysis:

- Fits PLS-DA model with specified number of components
- Extracts sample scores and loadings
- Calculates variance explained by each component
- Provides detailed logging of model fitting process

### 3. Sample Score Visualization

Creates comprehensive visualizations of sample classification:

- **Score plots**: Shows sample separation in PLS component space
- **Confidence ellipses**: 95% confidence regions for each group
- **Variance explained**: Displays percentage of variance explained by each component
- **Group separation**: Visual assessment of classification quality

### 4. VIP Score Analysis

Analyzes Variable Importance in Projection (VIP) scores:

- **VIP calculation**: Computes VIP scores for all proteins
- **Top VIP proteins**: Identifies most discriminatory proteins
- **VIP visualizations**: Creates barplots and lollipop plots
- **Feature ranking**: Ranks proteins by discriminatory power

### 5. Loading Analysis

Examines protein loadings on PLS components:

- **Loading extraction**: Extracts loadings for each component
- **Top loadings**: Identifies proteins with highest loadings
- **Loading visualization**: Creates barplots of top loadings
- **Direction analysis**: Distinguishes positive and negative loadings

### 6. Performance Evaluation

Evaluates model performance through cross-validation:

- **Cross-validation**: Performs k-fold cross-validation
- **Performance metrics**: Calculates classification error rates
- **Component selection**: Determines optimal number of components
- **Performance plots**: Visualizes error rates across components

### 7. Sparse PLS-DA for Feature Selection

Implements sparse PLS-DA for feature selection:

- **Parameter tuning**: Optimizes sparsity parameters
- **Feature selection**: Selects optimal number of features per component
- **Sparse model**: Fits sparse PLS-DA with selected features
- **Sparse visualization**: Creates plots of sparse PLS-DA results

### 8. GO Enrichment Analysis

Performs biological interpretation of discriminatory proteins:

- **Protein selection**: Combines top VIP and loading proteins
- **GO enrichment**: Performs GO Biological Process enrichment
- **Network analysis**: Creates gene-concept networks
- **Enrichment visualization**: Generates dotplots and lollipop plots

## Output Files

### Main Results
- **`plsda_results.RData`**: Complete PLS-DA results including models and scores
- **`top_discriminatory_proteins.csv`**: Top discriminatory proteins with scores
- **`plsda_config.rds`**: Configuration parameters used
- **`plsda_summary.txt`**: Comprehensive analysis summary

### Visualizations
- **`plsda_sample_scores.pdf`**: Sample classification scores
- **`vip_scores_barplot.pdf`**: VIP scores barplot
- **`vip_scores_lollipop.pdf`**: VIP scores lollipop plot
- **`loadings_barplot.pdf`**: Protein loadings barplot
- **`plsda_performance.pdf`**: Cross-validation performance
- **`sparse_plsda_plot.pdf`**: Sparse PLS-DA results
- **`go_enrichment_dotplot.pdf`**: GO enrichment dotplot
- **`go_enrichment_network.pdf`**: GO enrichment network
- **`go_enrichment_lollipop.pdf`**: GO enrichment lollipop plot

## Biological Interpretation

### Classification Performance

PLS-DA provides several performance metrics:
- **Classification accuracy**: Percentage of correctly classified samples
- **Cross-validation error**: Error rate from cross-validation
- **Component importance**: Variance explained by each component
- **Group separation**: Visual and statistical assessment of group separation

### Discriminatory Proteins

The pipeline identifies key discriminatory proteins through:
- **VIP scores**: Variable Importance in Projection scores
- **Loading weights**: Direct contribution to component construction
- **Feature selection**: Sparse PLS-DA for optimal feature selection
- **Biological relevance**: GO enrichment of discriminatory proteins

### Model Validation

Comprehensive validation includes:
- **Cross-validation**: K-fold cross-validation for performance assessment
- **Sparse modeling**: Feature selection to reduce overfitting
- **Performance curves**: Error rates across different numbers of components
- **Biological validation**: GO enrichment for biological relevance

## Usage Examples

### Basic Usage
```r
# Run with default parameters
source("05_PLSDA.R")
```

### Custom Configuration
```r
# Modify configuration before running
config$n_components <- 3
config$validation_folds <- 10
source("05_PLSDA.R")
```

### Integration with Previous Steps
```r
# Ensure preprocessing and normalization are complete
source("01_PreprocessData.R")
source("02_NormalizeData.R")
source("05_PLSDA.R")
```

## Troubleshooting

### Common Issues

1. **Insufficient samples per group**:
   - Ensure at least 3 samples per group
   - Consider combining similar groups if necessary

2. **Poor classification performance**:
   - Check data quality and preprocessing
   - Adjust number of components
   - Verify group assignments

3. **No significant GO enrichment**:
   - Verify gene symbols match annotation database
   - Check if discriminatory proteins are sufficient for enrichment

4. **Memory issues with large datasets**:
   - Reduce number of components
   - Use sparse PLS-DA for feature selection
   - Limit cross-validation folds

### Performance Optimization

For large datasets:
- Use sparse PLS-DA for feature selection
- Reduce number of cross-validation folds
- Limit number of components
- Use parallel processing when available

## Advanced Features

### Custom Component Selection
The script automatically uses 2 components, but you can adjust:
```r
# Manual component selection
config$n_components <- 3
```

### Custom Cross-Validation
Modify cross-validation parameters:
```r
config$validation_folds <- 10
config$n_repeats <- 20
```

### Integration with Differential Expression
Results can be integrated with differential expression from `02_DifferentialTesting.R`:
```r
# Load differential results
load("results/differential/differential_results.RData")
# Compare with PLS-DA discriminatory proteins
```

## Technical Details

### Algorithm Specifications

- **PLS-DA Algorithm**: Uses mixOmics implementation
- **Cross-Validation**: K-fold cross-validation with multiple repeats
- **Feature Selection**: Sparse PLS-DA with optimal parameter tuning
- **VIP Calculation**: Standard VIP score computation

### Statistical Methods

- **Classification**: PLS-DA for supervised classification
- **Performance Evaluation**: Cross-validation with error rate calculation
- **Feature Selection**: Sparse PLS-DA with parameter optimization
- **Multiple Testing**: Benjamini-Hochberg correction for GO enrichment

### Reproducibility

- Fixed random seeds for reproducible results
- Comprehensive logging of all parameters
- Saved configuration for exact reproduction

## Integration with Workflow

This script is designed to integrate seamlessly with the complete proteomics analysis pipeline:

1. **Preprocessing** (`00_PreprocessData.R`) → Clean data
2. **Normalization** (`01_DataNormalization.R`) → Normalize data
3. **Differential Testing** (`02_DifferentialTesting.R`) → Identify DE proteins
4. **Pathway Enrichment** (`03_PathwayEnrichment.R`) → Biological interpretation
5. **NMF Analysis** (`04_NMF.R`) → Module discovery
6. **PLS-DA Analysis** (`05_PLSDA.R`) → Sample classification ← **Current Step**

The PLS-DA results can be used to:
- Validate sample classification and group separation
- Identify key discriminatory proteins for biomarker discovery
- Compare with differential expression results
- Guide biological interpretation through GO enrichment
- Provide quality control for experimental design

## Citation

When using this PLS-DA pipeline, please cite:
- mixOmics: Rohart, F., et al. (2017). mixOmics: An R package for 'omics feature selection and multiple data integration. PLoS Computational Biology, 13(11), e1005752.
- clusterProfiler: Yu, G., et al. (2012). clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS, 16(5), 284-287.

## Support

For questions or issues with this pipeline:
1. Check the troubleshooting section above
2. Review the generated log files in the output directory
3. Verify input data format and requirements
4. Consult the main README for workflow integration 