# 00_PreprocessData.R - Proteomics Data Preprocessing Pipeline

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

### Raw_Intensity Sheet Structure
```
| Gene name | Sample1 | Sample2 | Sample3 | ... |
|-----------|---------|---------|---------|-----|
| GENE1     | 1234.5  | 1456.7  | 1345.2  | ... |
| GENE2     | 2345.6  | 2567.8  | 2456.3  | ... |
```

### SampleInfo Sheet Structure
```
| SampleID | Group | Condition | ... |
|----------|-------|-----------|-----|
| Sample1  | Ctrl  | Control   | ... |
| Sample2  | Ctrl  | Control   | ... |
```

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

### 1. Data Import and Validation
```r
# Validates file existence and structure
validate_input_file(config$input_file, required_sheets)

# Imports raw data and sample annotations
raw_data <- read_excel(config$input_file, sheet = config$raw_sheet)
sample_anno <- read_excel(config$input_file, sheet = config$info_sheet)
```

### 2. Data Matrix Construction
```r
# Converts to numeric matrix format
mat_raw <- raw_data %>%
  rename(Gene = `Gene name`) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  column_to_rownames("Gene") %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()
```

### 3. Missingness Analysis and Filtering

#### Gene-level Filtering
```r
# Calculate missingness rates per gene
gene_missing_rates <- rowMeans(is.na(mat_raw))

# Filter genes with excessive missingness
keep_genes <- gene_missing_rates <= config$gene_missing_threshold
mat_filt <- mat_raw[keep_genes, , drop = FALSE]
```

#### Sample-level Analysis
```r
# Calculate missingness rates per sample
sample_missing_rates <- colMeans(is.na(mat_filt))

# Remove only severely incomplete samples
remove_samps <- names(sample_missing_rates)[sample_missing_rates > config$sample_missing_threshold]
```

#### Final Cleanup
```r
# Remove any remaining genes with NA values
mat_clean <- mat_filt[rowSums(is.na(mat_filt)) == 0, , drop = FALSE]
```

### 4. Quality Assessment
```r
# Calculate retention statistics
gene_retention <- nrow(mat_clean) / nrow(mat_raw) * 100
sample_retention <- ncol(mat_clean) / ncol(mat_raw) * 100

# Data completeness metrics
data_completeness <- (1 - sum(is.na(mat_clean)) / length(mat_clean)) * 100
```

## Output Objects

### Primary Outputs
- `mat_clean`: Cleaned intensity matrix (proteins × samples)
- `sample_anno`: Sample annotation data frame

### Configuration and Metadata
- `config`: Configuration parameters used
- `gene_missing_rates`: Original gene missingness rates
- `sample_missing_rates`: Sample missingness rates after gene filtering

## Example Output

```
[14:30:15] INFO: Starting proteomics data preprocessing pipeline
[14:30:16] INFO: Importing data from Excel file
[14:30:16] INFO: Imported 2500 rows × 45 columns of raw data
[14:30:16] INFO: Constructing intensity matrix
[14:30:16] INFO: Intensity matrix: 2500 proteins × 45 samples
[14:30:16] INFO: === MISSINGNESS FILTERING ===
[14:30:16] INFO: Step 1: Gene-level missingness filtering
[14:30:16] INFO: genes missingness summary:
[14:30:16] INFO: Removing 500 genes with >20.0% missing values
[14:30:16] INFO: After gene filtering: 2000 proteins × 45 samples
[14:30:16] INFO: Step 2: Sample-level missingness analysis
[14:30:16] INFO: No samples exceed 50.0% missing; keeping all samples
[14:30:16] INFO: Step 3: Final cleanup - removing genes with any remaining NA values
[14:30:16] INFO: === PREPROCESSING SUMMARY ===
[14:30:16] INFO: Final dataset: 1950 proteins × 45 samples
[14:30:16] INFO: Gene retention: 78.0% (1950/2500)
[14:30:16] INFO: Sample retention: 100.0% (45/45)
[14:30:16] INFO: Preprocessing pipeline completed successfully
```

## Quality Metrics

### Data Retention
- **Gene retention**: Percentage of proteins retained after filtering
- **Sample retention**: Percentage of samples retained after filtering
- **Data completeness**: Percentage of non-missing values in final dataset

### Missingness Patterns
- **Gene-level missingness**: Distribution of missing values across proteins
- **Sample-level missingness**: Distribution of missing values across samples
- **Missingness correlation**: Patterns in missing data structure

## Error Handling

### File Validation
- Checks file existence and accessibility
- Validates Excel file structure and required sheets
- Ensures data format compatibility

### Data Validation
- Handles duplicate gene names (keeps first occurrence)
- Converts text values to numeric with NA for parse failures
- Validates matrix dimensions and data types

### Error Messages
- `"Input file not found"`: Check file path in configuration
- `"Missing required sheets"`: Verify Excel file structure
- `"Error reading Excel file"`: Check file format and permissions

## Usage

### Basic Usage
```r
# Run the preprocessing pipeline
source("00_PreprocessData.R")

# Access cleaned data
dim(mat_clean)  # Check dimensions
head(mat_clean) # View first few rows
```

### Custom Configuration
```r
# Modify configuration before running
config$gene_missing_threshold <- 0.15  # More aggressive filtering
config$sample_missing_threshold <- 0.60 # More lenient sample filtering

# Run pipeline with custom settings
source("00_PreprocessData.R")
```

## Next Steps

After running this preprocessing pipeline, the cleaned data (`mat_clean`) is ready for:

1. **Normalization**: Run `01_DataNormalization.R` for VSN and scaling
2. **Quality assessment**: Evaluate data quality and completeness
3. **Downstream analysis**: Statistical testing, clustering, or machine learning

## Technical Notes

### Missing Value Strategy
- **Conservative approach**: Only removes samples with extreme missingness
- **Statistical rationale**: Preserves maximum information for downstream analysis
- **Quality focus**: Prioritizes data quality over quantity

### Performance Considerations
- **Memory efficient**: Processes data in chunks for large datasets
- **Fast execution**: Optimized for typical proteomics dataset sizes
- **Scalable**: Handles datasets with thousands of proteins and samples

### Reproducibility
- **Session info capture**: Optional saving of R session information
- **Detailed logging**: Timestamped progress tracking
- **Configuration tracking**: All parameters documented and accessible

## Troubleshooting

### Common Issues
1. **File path errors**: Ensure Excel file is in correct location
2. **Memory issues**: For large datasets, increase R memory limit
3. **Format errors**: Verify Excel file structure matches requirements

### Performance Tips
1. **Large datasets**: Consider processing in batches
2. **Memory optimization**: Close unnecessary R objects before running
3. **Parallel processing**: For very large datasets, consider parallel execution

## Citation

If you use this preprocessing pipeline in your research, please cite:

```
Joseph, K. (2024). Proteomics Data Preprocessing Pipeline. 
GitHub repository. https://github.com/yourusername/proteomics-preprocessing
```

## Contact

For questions or issues with this preprocessing pipeline, please open an issue on GitHub or contact the maintainer. 