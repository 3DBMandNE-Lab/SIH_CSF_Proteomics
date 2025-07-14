# Proteomics Data Preprocessing Pipeline

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

## Configuration

Modify the `config` list in the script to adjust parameters:

```r
config <- list(input_file = "data/Raw_Proteomics.xlsx", # File paths
               # Filtering thresholds
               gene_missing_threshold = 0.20,    # 20% missing threshold for genes
               sample_missing_threshold = 0.50,  # 50% missing threshold for samples
               # Output settings
               verbose = TRUE,                   # Detailed logging
               save_intermediate = FALSE         # Save session info
               )
```

## Pipeline Steps

### 1. Data Import and Validation
- Validates file existence and Excel structure
- Checks for required sheets
- Imports raw intensity data and sample annotations

### 2. Data Matrix Construction
- Converts raw data to numeric matrix format
- Handles duplicate gene names
- Creates protein × sample intensity matrix

### 3. Missingness Analysis and Filtering
- **Gene-level filtering**: Removes proteins with >20% missing values
- **Sample-level analysis**: Identifies samples with excessive missingness
- **Final cleanup**: Removes any remaining genes with NA values

### 4. Quality Assessment
- Calculates retention statistics for genes and samples
- Reports data completeness metrics
- Provides detailed summary of filtering results

## Output Objects

After running the pipeline, the following objects are available:

- `mat_clean`: Cleaned intensity matrix (proteins × samples)
- `sample_anno`: Sample annotation data frame
- `config`: Configuration parameters used
- `gene_missing_rates`: Original gene missingness rates
- `sample_missing_rates`: Sample missingness rates after gene filtering

## Next Steps

The cleaned data is ready for downstream analysis including:

1. **Normalization**: Variance stabilizing normalization (VSN) or quantile normalization
2. **Log transformation**: Log2 transformation for statistical analysis
3. **Outlier detection**: Sample and protein outlier identification
4. **Principal Component Analysis (PCA)**: Dimensionality reduction and quality assessment
5. **Differential expression analysis**: Statistical testing for protein abundance changes
