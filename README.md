# Proteomics Data Analysis Pipeline

A comprehensive R pipeline for preprocessing proteomics intensity data with robust filtering, validation, and detailed logging.

## Overview

This pipeline implements a systematic approach to proteomics data preprocessing with the following key features:

- **Aggressive protein filtering**: Removes proteins with excessive missingness (>20% by default)
- **Minimal sample removal**: Only removes samples with extreme missingness (>50% by default)
- **Comprehensive validation**: File structure and data integrity checks
- **Detailed logging**: Step-by-step progress tracking with timestamps
- **Configurable parameters**: Easy modification of thresholds and settings
- **Reproducible analysis**: Session information capture for reproducibility

## Requirements

### R Packages
```r
install.packages(c("readxl", "limma", "vsn", "tidyverse", "pheatmap", "factoextra"))
```

### Data Format
The pipeline expects an Excel file with two sheets:
- **Raw_Intensity**: Contains protein intensity data with columns for gene names and sample intensities
- **SampleInfo**: Contains sample metadata with a "SampleID" column

## Quick Start

1. **Prepare your data**: Place your Excel file in the `data/` directory as `Raw_Proteomics.xlsx`

2. **Run the preprocessing pipeline**:
```r
source("00_PreprocessData.R")
```

3. **Run the normalization pipeline**:
```r
source("01_DataNormalization.R")
```

4. **Run the differential analysis pipeline**:
```r
source("02_DifferentialTesting.R")
```

5. **Access results**: 
   - Preprocessed data: `mat_clean` (cleaned intensity matrix)
   - Normalized data: `mat_scaled` (normalized and scaled matrix)
   - Differential results: `de_res` (differential abundance results)
   - Sample annotations: `sample_anno` (updated after outlier removal)

## Configuration

Modify the `config` list in the script to adjust parameters:

```r
config <- list(
  # File paths
  input_file = "data/Raw_Proteomics.xlsx",
  
  # Filtering thresholds
  gene_missing_threshold = 0.20,    # 20% missing threshold for genes
  sample_missing_threshold = 0.50,  # 50% missing threshold for samples
  
  # Output settings
  verbose = TRUE,                   # Detailed logging
  save_intermediate = FALSE         # Save session info
)
```

## Pipeline Steps

### Preprocessing Pipeline (00_PreprocessData.R)

#### 1. Data Import and Validation
- Validates file existence and Excel structure
- Checks for required sheets
- Imports raw intensity data and sample annotations

#### 2. Data Matrix Construction
- Converts raw data to numeric matrix format
- Handles duplicate gene names
- Creates protein × sample intensity matrix

#### 3. Missingness Analysis and Filtering
- **Gene-level filtering**: Removes proteins with >20% missing values
- **Sample-level analysis**: Identifies samples with excessive missingness
- **Final cleanup**: Removes any remaining genes with NA values

#### 4. Quality Assessment
- Calculates retention statistics for genes and samples
- Reports data completeness metrics
- Provides detailed summary of filtering results

### Normalization Pipeline (01_DataNormalization.R)

#### 1. Variance Stabilizing Normalization (VSN)
- Applies VSN for intensity-dependent variance stabilization
- Handles heteroscedasticity in proteomics data
- Provides statistical analysis of normalization effects

#### 2. Scaling and Centering
- Performs Z-score scaling and centering
- Prepares data for downstream statistical analysis
- Analyzes scaling effects on data distribution

#### 3. Outlier Detection
- Implements PCA-based outlier detection
- Uses Mahalanobis distance with configurable thresholds
- Automatically removes outlier samples

#### 4. Quality Assessment Visualizations
- Generates PCA plots for sample quality assessment
- Creates heatmaps of top variable proteins
- Provides comprehensive quality metrics

### Differential Analysis Pipeline (02_DifferentialTesting.R)

#### 1. Experimental Design
- Creates design matrix for group comparisons
- Handles various experimental designs
- Validates group assignments and sample sizes

#### 2. Statistical Analysis
- Implements empirical Bayes moderated t-tests
- Applies multiple testing correction (FDR)
- Provides robust statistical inference

#### 3. Results Analysis
- Identifies differentially abundant proteins
- Calculates effect sizes and significance
- Generates comprehensive statistical summaries

#### 4. Visualization
- Creates publication-ready volcano plots
- Generates heatmaps of top differential proteins
- Provides quality assessment visualizations

## Output Objects

After running the preprocessing pipeline, the following objects are available:

- `mat_clean`: Cleaned intensity matrix (proteins × samples)
- `sample_anno`: Sample annotation data frame
- `config`: Configuration parameters used
- `gene_missing_rates`: Original gene missingness rates
- `sample_missing_rates`: Sample missingness rates after gene filtering

After running the normalization pipeline, additional objects are available:

- `mat_scaled`: Normalized and scaled intensity matrix (proteins × samples)
- `mat_norm`: VSN-normalized matrix (before scaling)
- `vsn_fit`: VSN fit object for future predictions
- `pca_res`: PCA results for quality assessment
- `norm_config`: Normalization configuration parameters used

After running the differential analysis pipeline, additional objects are available:

- `de_res`: Differential abundance results (limma topTable output)
- `fit2`: Fitted limma model with empirical Bayes results
- `volcano_plot`: Enhanced volcano plot with background shading
- `heatmap_plot`: Heatmap of top differential proteins
- `diff_config`: Differential analysis configuration parameters used

## License

This project is licensed under the MIT License - see the LICENSE file for details.
