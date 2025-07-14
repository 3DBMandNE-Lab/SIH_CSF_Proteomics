# 00_PreprocessData.R

## Overview

This script implements a comprehensive preprocessing pipeline for proteomics intensity data with aggressive protein filtering and minimal sample removal strategies. The pipeline is designed to handle high-throughput proteomics data while preserving maximum sample information.

## Purpose

The preprocessing pipeline addresses common challenges in proteomics data analysis:
- **Missing value handling**: Systematic approach to missing data in proteomics
- **Data quality assessment**: Comprehensive evaluation of data completeness
- **Feature selection**: Intelligent protein filtering based on missingness patterns
- **Sample preservation**: Conservative sample removal strategy

## Technical Approach

### Aggressive Protein Filtering
- **Primary filter**: Removes proteins with >20% missing values across samples
- **Secondary filter**: Eliminates proteins with any remaining NA values after sample filtering
- **Rationale**: Proteins with excessive missingness provide limited statistical power

### Minimal Sample Removal
- **Threshold**: Only removes samples with >50% missing values
- **Strategy**: Preserves maximum sample information for downstream analysis
- **Rationale**: Sample removal should be conservative to maintain statistical power

## Dependencies

```r
# Required R packages
library(readxl)      # Excel file reading
library(limma)       # Linear models for microarray data
library(vsn)         # Variance stabilizing normalization
library(tidyverse)   # Data manipulation and visualization
library(pheatmap)    # Heatmap generation
library(factoextra)  # PCA visualization
```

## Input Requirements

### Data Format
- **File type**: Excel (.xlsx) file
- **Location**: `data/Raw_Proteomics.xlsx` (configurable)
- **Required sheets**: 
  - `Raw_Intensity`: Protein intensity data
  - `SampleInfo`: Sample metadata

## Configuration Parameters

```r
config <- list(
  # File paths
  input_file = "data/Raw_Proteomics.xlsx",
  raw_sheet = "Raw_Intensity",
  info_sheet = "SampleInfo",
  
  # Filtering thresholds
  gene_missing_threshold = 0.20,    # 20% missing threshold for genes
  sample_missing_threshold = 0.50,  # 50% missing threshold for samples
  
  # Output settings
  verbose = TRUE,
  save_intermediate = FALSE
)
```

## Pipeline Steps

 - Data Import and Validation
 - Data Matrix Construction
 - Missingness Analysis and Filtering
 - Quality Assessment


## Output Objects

### Primary Outputs
- `mat_clean`: Cleaned intensity matrix (proteins × samples)
- `sample_anno`: Sample annotation data frame

### Configuration and Metadata
- `config`: Configuration parameters used
- `gene_missing_rates`: Original gene missingness rates
- `sample_missing_rates`: Sample missingness rates after gene filtering

## Usage

### Basic Usage
```r
# Run the preprocessing pipeline
source("00_PreprocessData.R")

# Access cleaned data
dim(mat_clean)  # Check dimensions
head(mat_clean) # View first few rows
```
