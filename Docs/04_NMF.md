# NMF-Based Module Generation Pipeline

## Overview

The `04_NMF.R` script performs Non-negative Matrix Factorization (NMF) analysis on preprocessed proteomics data to identify co-expression modules and their biological significance. This pipeline is designed to discover latent patterns in protein expression data and provide comprehensive biological interpretation through module-based analysis.

## Purpose

NMF analysis serves as a powerful dimensionality reduction and pattern discovery tool that:
- Identifies co-expression modules in proteomics data
- Provides biological interpretation through GO enrichment analysis
- Enables integration with differential expression results
- Offers quality control metrics for module assessment

## Prerequisites

### Input Data Requirements

The script expects the following data structures from previous pipeline steps:

1. **`mat_clean`**: A numeric matrix (proteins × samples) with:
   - No missing values (NAs)
   - Non-negative values (NMF requirement)
   - Row names as protein/gene identifiers
   - Column names as sample identifiers

2. **`sample_anno`**: A data frame with:
   - Row names matching sample names in `mat_clean`
   - A `Group` column indicating sample conditions
   - Additional annotation columns as needed

3. **`de_res`** (optional): Differential expression results for integration analysis

### Dependencies

The script requires the following R packages:
- **Core NMF**: `NMF` for matrix factorization
- **Visualization**: `pheatmap`, `ggplot2`, `patchwork`
- **Data manipulation**: `dplyr`, `tibble`
- **Biological analysis**: `clusterProfiler`, `org.Hs.eg.db`
- **Statistical analysis**: `cluster`, `ggpubr`
- **Additional plotting**: `ggridges`, `ggalluvial`, `tidytext`

## Configuration

### Key Parameters

The script includes a comprehensive configuration system:

```r
config <- list(
  # Input/Output paths
  input_dir = "data/processed",
  output_dir = "results/nmf",
  
  # NMF parameters
  rank_range = 2:6,           # Range of ranks to test
  n_runs_consensus = 20,      # Number of runs for consensus clustering
  n_runs_final = 50,          # Number of runs for final NMF
  top_genes_per_module = 15,  # Top genes to visualize per module
  max_genes_enrichment = 200, # Maximum genes per module for enrichment
  
  # Visualization parameters
  plot_width = 12,
  plot_height = 16,
  dpi = 300,
  
  # Statistical parameters
  pvalue_cutoff = 0.05,
  qvalue_cutoff = 0.05,
  
  # Random seed for reproducibility
  seed = 123
)
```

### Parameter Descriptions

- **`rank_range`**: Range of ranks to test for optimal module number selection
- **`n_runs_consensus`**: Number of NMF runs for consensus clustering stability
- **`n_runs_final`**: Number of runs for final NMF at optimal rank
- **`top_genes_per_module`**: Number of top genes to visualize per module
- **`max_genes_enrichment`**: Maximum genes per module for GO enrichment analysis

## Pipeline Steps

### 1. Optimal Rank Estimation

The script begins by estimating the optimal number of modules using consensus clustering:

- Runs NMF across multiple ranks (default: 2-6)
- Computes cophenetic coefficients for rank selection
- Generates consensus clustering plots
- Automatically selects optimal rank based on stability metrics

### 2. NMF Factorization

Performs NMF at the optimal rank to decompose the data:

- Factorizes the data matrix into W (genes × modules) and H (modules × samples) matrices
- Assigns meaningful names to modules
- Provides detailed logging of factorization process

### 3. Gene-to-Module Assignment

Assigns each gene to its most representative module:

- Uses maximum loading criterion for assignment
- Computes loading scores for each gene
- Generates summary statistics for each module
- Saves assignments to CSV file

### 4. Module Metagene Visualization

Creates comprehensive visualizations of module patterns:

- **Metagene heatmap**: Shows module activity across samples
- **Gene loadings heatmap**: Displays top genes per module
- **Module scores barplot**: Individual sample scores by module

### 5. GO Enrichment Analysis

Performs biological interpretation of modules:

- Uses top genes from each module for enrichment
- Performs GO Biological Process enrichment
- Creates network plots for each module
- Handles cases with no significant enrichment

### 6. Quality Control and Statistical Analysis

Provides comprehensive quality assessment:

- **Statistical comparisons**: Module scores between groups
- **Silhouette analysis**: Module quality assessment
- **Integration analysis**: Links with differential expression results
- **Alluvial plots**: Flow of DE genes into modules

## Output Files

### Main Results
- **`nmf_results.RData`**: Complete NMF results including W, H matrices
- **`gene_module_assignments.csv`**: Gene-to-module assignments with loading scores
- **`nmf_config.rds`**: Configuration parameters used
- **`nmf_summary.txt`**: Comprehensive analysis summary

### Visualizations
- **`consensus_clustering.pdf`**: Consensus clustering results
- **`module_metagenes_heatmap.pdf`**: Module activity across samples
- **`top_gene_loadings_heatmap.pdf`**: Top genes per module
- **`module_scores_barplot.pdf`**: Sample scores by module
- **`module_GO_enrichment.pdf`**: GO enrichment networks
- **`module_comparison_boxplot.pdf`**: Statistical comparisons
- **`module_silhouette_analysis.pdf`**: Module quality assessment
- **`DE_genes_alluvial.pdf`**: Integration with DE results (if available)

## Biological Interpretation

### Module Characteristics

Each identified module represents:
- **Co-expressed proteins**: Proteins with similar expression patterns
- **Biological processes**: Enriched GO terms indicating functional coherence
- **Sample patterns**: Distinct activity patterns across experimental conditions

### Quality Metrics

The pipeline provides several quality metrics:
- **Silhouette width**: Measures how well genes cluster within modules
- **Statistical significance**: P-values for group comparisons
- **Enrichment significance**: Adjusted p-values for GO terms

### Integration with Differential Expression

When DE results are available, the pipeline:
- Maps differentially expressed genes to modules
- Creates alluvial plots showing gene flow
- Provides biological context for DE findings

## Usage Examples

### Basic Usage
```r
# Run with default parameters
source("04_NMF.R")
```

### Custom Configuration
```r
# Modify configuration before running
config$rank_range <- 3:8
config$n_runs_final <- 100
source("04_NMF.R")
```

### Integration with Previous Steps
```r
# Ensure preprocessing and normalization are complete
source("01_PreprocessData.R")
source("02_NormalizeData.R")
source("03_DifferentialTesting.R")
source("04_NMF.R")
```

## Troubleshooting

### Common Issues

1. **Negative values in data**:
   - Ensure preprocessing step removed negative values
   - Check data transformation in previous steps

2. **Missing GO enrichment**:
   - Verify gene symbols match annotation database
   - Check if modules have sufficient genes for enrichment

3. **Memory issues with large datasets**:
   - Reduce `n_runs_consensus` and `n_runs_final`
   - Limit `max_genes_enrichment`

### Performance Optimization

For large datasets:
- Reduce the number of NMF runs
- Limit the number of genes per module for enrichment
- Use parallel processing when available

## Advanced Features

### Custom Rank Selection
The script automatically selects optimal rank, but you can override:
```r
# Manual rank selection based on biological knowledge
opt_rank <- 4  # Set before NMF factorization
```

### Custom Group Comparisons
Modify group comparisons for statistical tests:
```r
my_comp <- list(c("Control", "Treatment"))
```

### Integration with Pathway Analysis
Results can be integrated with pathway enrichment from `03_PathwayEnrichment.R`:
```r
# Load pathway results
load("results/pathway/pathway_results.RData")
# Integrate with NMF modules
```

## Technical Details

### Algorithm Specifications

- **NMF Algorithm**: Uses the default NMF algorithm from the NMF package
- **Consensus Clustering**: Multiple runs with different random seeds
- **Rank Selection**: Based on cophenetic coefficient stability
- **Gene Assignment**: Maximum loading criterion

### Statistical Methods

- **Module Comparison**: Wilcoxon rank-sum tests
- **Multiple Testing**: Benjamini-Hochberg correction
- **Quality Assessment**: Silhouette analysis

### Reproducibility

- Fixed random seeds for reproducible results
- Comprehensive logging of all parameters
- Saved configuration for exact reproduction

## Integration with Workflow

This script is designed to integrate seamlessly with the complete proteomics analysis pipeline:

1. **Preprocessing** (`01_PreprocessData.R`) → Clean data
2. **Normalization** (`02_NormalizeData.R`) → Normalize data
3. **Differential Testing** (`03_DifferentialTesting.R`) → Identify DE proteins
4. **NMF Analysis** (`04_NMF.R`) → Module discovery ← **Current Step**
5. **Pathway Enrichment** (`03_PathwayEnrichment.R`) → Biological interpretation

The NMF results can be used to:
- Guide pathway analysis by focusing on module-specific genes
- Provide biological context for differential expression results
- Identify novel protein interaction networks
- Generate hypotheses for follow-up experiments

## Citation

When using this NMF pipeline, please cite:
- NMF package: Gaujoux, R., & Seoighe, C. (2010). A flexible R package for nonnegative matrix factorization. BMC Bioinformatics, 11(1), 367.
- clusterProfiler: Yu, G., et al. (2012). clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS, 16(5), 284-287.

## Support

For questions or issues with this pipeline:
1. Check the troubleshooting section above
2. Review the generated log files in the output directory
3. Verify input data format and requirements
4. Consult the main README for workflow integration 