# 03_PathwayEnrichment.R - Proteomics Pathway Enrichment Analysis Pipeline

## Overview

This script implements a comprehensive pathway enrichment analysis pipeline for proteomics data using advanced bioinformatics tools. It performs both over-representation analysis (ORA) and gene set enrichment analysis (GSEA) to identify biological pathways and processes associated with differentially abundant proteins, with publication-ready visualizations and integrated dashboards.

## Purpose

The pathway enrichment analysis pipeline addresses critical challenges in proteomics data interpretation:
- **Functional annotation**: Identifies biological processes associated with differential proteins
- **Pathway discovery**: Discovers enriched pathways using multiple statistical approaches
- **Visualization**: Generates comprehensive publication-ready visualizations
- **Integration**: Combines ORA and GSEA for comprehensive pathway analysis
- **Interpretation**: Provides biological context for differential abundance results

## Technical Approach

### Over-Representation Analysis (ORA)
- **Method**: Uses `enrichGO()` from clusterProfiler package
- **Rationale**: Tests whether differentially abundant proteins are over-represented in specific GO terms
- **Benefits**: Identifies specific biological processes enriched in differential proteins

### Gene Set Enrichment Analysis (GSEA)
- **Method**: Uses `gseGO()` for genome-wide pathway analysis
- **Rationale**: Tests whether gene sets show coordinated changes in expression
- **Benefits**: Captures subtle but coordinated changes across gene sets

### Multiple Visualization Types
- **Enhanced volcano plots**: Publication-ready differential abundance visualization
- **Network plots**: Gene-concept networks showing relationships
- **Enrichment maps**: Network visualization of enriched terms
- **Heatmaps**: Core enrichment visualization
- **UpSet plots**: Gene overlap analysis across terms

## Dependencies

```r
# Required R packages
library(clusterProfiler)  # Gene set enrichment analysis
library(enrichplot)       # Enhanced enrichment visualizations
library(GOplot)          # Circular network visualizations
library(ggplot2)         # Advanced plotting
library(patchwork)       # Multi-panel plot assembly
library(viridis)         # Color palettes
library(ggrepel)         # Label positioning
library(pheatmap)        # Heatmap generation
library(dplyr)           # Data manipulation
```

## Input Requirements

### Prerequisites
- **Required input**: `de_res` (from `02_DifferentialTesting.R`)
- **Required input**: `mat_scaled` (normalized data matrix)
- **Required input**: `sample_anno` (sample annotations)
- **Data format**: Differential analysis results with logFC and adj.P.Val columns

### Data Quality Requirements
- **Differential results**: Must contain logFC and adj.P.Val columns
- **Gene identifiers**: Must be compatible with organism annotation database
- **Statistical thresholds**: Configurable FDR and logFC thresholds
- **Valid dimensions**: Sufficient genes for meaningful enrichment analysis

## Configuration Parameters

```r
enrich_config <- list(
  # Analysis settings
  organism = "org.Hs.eg.db",     # Organism annotation database
  key_type = "SYMBOL",           # Gene identifier type
  ontology = "ALL",              # GO ontology (BP, MF, CC, ALL)
  
  # Statistical thresholds
  pvalue_cutoff = 0.05,         # P-value threshold for enrichment
  fdr_cutoff = 0.05,            # FDR threshold for enrichment
  logfc_threshold = 1.0,        # Log2 fold-change threshold
  
  # GSEA settings
  min_gssize = 10,              # Minimum gene set size
  max_gssize = 500,             # Maximum gene set size
  gsea_pvalue_cutoff = 0.05,    # GSEA p-value threshold
  
  # Visualization settings
  show_category = 20,            # Number of categories to show
  dotplot_show = 20,            # Categories for dotplot
  network_show = 10,            # Categories for network plots
  heatmap_show = 15,            # Categories for heatmap
  upset_show = 10,              # Categories for upset plot
  
  # Plot settings
  volcano_point_size = 2,        # Point size for volcano plot
  volcano_alpha = 0.8,          # Point transparency
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
if (!exists("de_res")) {
  stop("Differential analysis results 'de_res' not found. Please run 02_DifferentialTesting.R first.")
}

# Validate input data structure
validate_enrichment_input(de_res, mat_scaled, sample_anno)
```

### 2. Gene Set Preparation
```r
# Prepare gene lists for enrichment analysis
gene_sets <- prepare_gene_lists(de_res, enrich_config)

# Create named gene list for GSEA
gene_list <- differential_results$logFC
names(gene_list) <- rownames(differential_results)
gene_list <- sort(gene_list, decreasing = TRUE)
```

### 3. Over-Representation Analysis (ORA)
```r
# Load organism database
org_db <- get(enrich_config$organism)

# Perform GO enrichment for up-regulated genes
ego_up <- enrichGO(
  gene = gene_sets$up_genes,
  OrgDb = org_db,
  ont = enrich_config$ontology,
  pAdjustMethod = "BH",
  keyType = enrich_config$key_type,
  pvalueCutoff = enrich_config$pvalue_cutoff
)

# Perform GO enrichment for down-regulated genes
ego_down <- enrichGO(
  gene = gene_sets$down_genes,
  OrgDb = org_db,
  ont = enrich_config$ontology,
  pAdjustMethod = "BH",
  keyType = enrich_config$key_type,
  pvalueCutoff = enrich_config$pvalue_cutoff
)
```

### 4. Gene Set Enrichment Analysis (GSEA)
```r
# Perform GSEA on GO Biological Process
gsea_res <- gseGO(
  geneList = gene_sets$gene_list,
  keyType = enrich_config$key_type,
  OrgDb = org_db,
  ont = "BP",
  minGSSize = enrich_config$min_gssize,
  maxGSSize = enrich_config$max_gssize,
  pvalueCutoff = enrich_config$gsea_pvalue_cutoff,
  verbose = enrich_config$verbose
)
```

### 5. Visualization Generation
```r
# Create enhanced volcano plot
volcano_plot <- create_enhanced_volcano(de_res, enrich_config)

# Generate GSEA dotplot
gsea_dotplot <- dotplot(gsea_res,
                        showCategory = enrich_config$dotplot_show,
                        title = "GSEA of GO Biological Processes")

# Create enrichment map
sim <- pairwise_termsim(gsea_res)
enrich_map <- emapplot(sim, showCategory = enrich_config$network_show)

# Generate gene-concept network
cnet_plot <- cnetplot(gsea_res, 
                      foldChange = gene_sets$gene_list, 
                      showCategory = enrich_config$network_show)
```

## Output Objects

### Primary Analysis Results
- `ego_up`: GO enrichment results for up-regulated genes
- `ego_down`: GO enrichment results for down-regulated genes
- `gsea_res`: GSEA results for all genes
- `gene_sets`: Prepared gene lists for analysis

### Visualization Objects
- `volcano_plot`: Enhanced volcano plot with background shading
- `heatmap_plot`: Heatmap of top differential proteins
- `gsea_dotplot`: GSEA enrichment dotplot
- `enrich_map`: Enrichment map network
- `cnet_plot`: Gene-concept network
- `heat_plot`: Core enrichment heatmap
- `upset_plot`: Gene overlap upset plot
- `dashboard`: Integrated multi-panel dashboard

### Configuration and Metadata
- `enrich_config`: Configuration parameters used

## Example Output

```
[14:33:00] ENRICH_INFO: Starting proteomics pathway enrichment analysis
[14:33:00] ENRICH_INFO: Input validation passed: 1950 proteins in differential results
[14:33:00] ENRICH_INFO: === GENE SET PREPARATION ===
[14:33:00] ENRICH_INFO: Gene sets prepared: 89 up-regulated, 67 down-regulated
[14:33:00] ENRICH_INFO: === OVER-REPRESENTATION ANALYSIS ===
[14:33:00] ENRICH_INFO: Loading organism annotation database
[14:33:00] ENRICH_INFO: Performing GO enrichment for up-regulated genes
[14:33:00] ENRICH_INFO: Found 15 enriched terms for up-regulated genes
[14:33:00] ENRICH_INFO: Performing GO enrichment for down-regulated genes
[14:33:00] ENRICH_INFO: Found 12 enriched terms for down-regulated genes
[14:33:00] ENRICH_INFO: === GENE SET ENRICHMENT ANALYSIS ===
[14:33:00] ENRICH_INFO: Performing GSEA on GO Biological Process
[14:33:00] ENRICH_INFO: GSEA completed: 25 enriched gene sets found
[14:33:00] ENRICH_INFO: === VISUALIZATION ===
[14:33:00] ENRICH_INFO: Generating enhanced volcano plot
[14:33:00] ENRICH_INFO: Generating heatmap of top differential proteins
[14:33:00] ENRICH_INFO: Generating GSEA dotplot
[14:33:00] ENRICH_INFO: Generating enrichment map
[14:33:00] ENRICH_INFO: Generating gene-concept network
[14:33:00] ENRICH_INFO: Generating core enrichment heatmap
[14:33:00] ENRICH_INFO: Generating upset plot
[14:33:00] ENRICH_INFO: === INTEGRATED DASHBOARD ===
[14:33:00] ENRICH_INFO: Creating integrated enrichment dashboard
[14:33:00] ENRICH_INFO: === RESULTS EXPORT ===
[14:33:00] ENRICH_INFO: Up-regulated GO enrichment results saved
[14:33:00] ENRICH_INFO: Down-regulated GO enrichment results saved
[14:33:00] ENRICH_INFO: GSEA results saved
[14:33:00] ENRICH_INFO: === ENRICHMENT ANALYSIS SUMMARY ===
[14:33:00] ENRICH_INFO: Analysis completed for 1950 total genes
[14:33:00] ENRICH_INFO: Up-regulated genes: 89
[14:33:00] ENRICH_INFO: Down-regulated genes: 67
[14:33:00] ENRICH_INFO: Up-regulated enriched terms: 15
[14:33:00] ENRICH_INFO: Down-regulated enriched terms: 12
[14:33:00] ENRICH_INFO: GSEA enriched gene sets: 25
[14:33:00] ENRICH_INFO: Pathway enrichment analysis pipeline completed successfully
```

## Quality Metrics

### Enrichment Statistics
- **Total genes**: Number of genes analyzed
- **Up-regulated genes**: Genes with positive logFC and significant FDR
- **Down-regulated genes**: Genes with negative logFC and significant FDR
- **Enriched terms**: Number of significantly enriched GO terms

### GSEA Metrics
- **Enriched gene sets**: Number of gene sets with significant enrichment
- **Normalized enrichment scores**: Distribution of NES values
- **Gene set sizes**: Range of gene set sizes analyzed

### Visualization Quality
- **Publication-ready plots**: High-resolution visualizations
- **Multi-panel dashboards**: Integrated analysis summaries
- **Interactive elements**: Configurable plot parameters

## Error Handling

### Input Validation
- Checks for required input objects (`de_res`, `mat_scaled`, `sample_anno`)
- Validates differential results structure and required columns
- Ensures gene identifiers are compatible with annotation database

### Analysis Validation
- Validates gene set sizes for meaningful enrichment
- Checks for sufficient genes in each direction
- Ensures organism database is available

### Error Messages
- `"Differential analysis results 'de_res' not found"`: Run differential analysis first
- `"Missing required columns in differential results"`: Check differential results format
- `"No enriched terms found"`: Adjust thresholds or check gene sets

## Usage

### Basic Usage
```r
# Run the pathway enrichment pipeline
source("03_PathwayEnrichment.R")

# Access results
head(ego_up)      # View up-regulated enrichment results
head(ego_down)    # View down-regulated enrichment results
head(gsea_res)    # View GSEA results
```

### Custom Configuration
```r
# Modify configuration before running
enrich_config$fdr_cutoff <- 0.01  # More stringent FDR control
enrich_config$logfc_threshold <- 0.5  # Lower effect size threshold
enrich_config$organism <- "org.Mm.eg.db"  # Change to mouse database
enrich_config$save_plots <- TRUE  # Save plots to files

# Run pipeline with custom settings
source("03_PathwayEnrichment.R")
```

### Advanced Usage
```r
# Access specific enrichment results
up_terms <- ego_up@result
down_terms <- ego_down@result
gsea_terms <- gsea_res@result

# Create custom visualizations
custom_dotplot <- dotplot(ego_up, showCategory = 10)
custom_network <- cnetplot(ego_up, foldChange = gene_sets$gene_list)

# Export specific results
write.csv(up_terms, "custom_up_enrichment.csv")
```

## Next Steps

After running this pathway enrichment pipeline, the results are ready for:

1. **KEGG pathway analysis**: Additional pathway database analysis
2. **Reactome pathway analysis**: Curated pathway database analysis
3. **Protein-protein interaction networks**: Network analysis of enriched proteins
4. **Functional annotation**: Detailed functional analysis of enriched terms
5. **Biomarker prioritization**: Selection of candidate biomarkers
6. **Meta-analysis**: Integration with other datasets or studies

## Technical Notes

### ORA Implementation
- **Statistical foundation**: Hypergeometric test for over-representation
- **Multiple testing correction**: Benjamini-Hochberg FDR correction
- **Gene set filtering**: Configurable size limits for meaningful analysis

### GSEA Implementation
- **Rank-based analysis**: Uses gene ranking by differential abundance
- **Permutation testing**: Statistical significance via permutation
- **Normalized enrichment scores**: Standardized effect size measures

### Performance Considerations
- **Memory efficient**: Optimized for large gene sets
- **Fast execution**: Efficient statistical computations
- **Scalable**: Handles datasets with thousands of genes

### Reproducibility
- **Detailed logging**: Timestamped progress tracking
- **Configuration tracking**: All parameters documented and accessible
- **Export capabilities**: Results and plots can be saved for publication

## Troubleshooting

### Common Issues
1. **No enriched terms**: Check statistical thresholds or gene set sizes
2. **Database errors**: Ensure organism annotation database is installed
3. **Memory issues**: For large datasets, increase R memory limit
4. **Gene identifier mismatches**: Verify gene ID format compatibility

### Performance Tips
1. **Large datasets**: Consider processing in batches
2. **Memory optimization**: Close unnecessary R objects before running
3. **Plot saving**: Use `save_plots = TRUE` for publication-quality figures

### Alternative Approaches
1. **Different databases**: Use KEGG, Reactome, or other pathway databases
2. **Alternative ontologies**: Use Molecular Function or Cellular Component
3. **Custom gene sets**: Use user-defined gene sets for enrichment

## Citation

If you use this pathway enrichment pipeline in your research, please cite:

```
Joseph, K. (2024). Proteomics Pathway Enrichment Analysis Pipeline. 
GitHub repository. https://github.com/yourusername/proteomics-preprocessing
```

## Contact

For questions or issues with this pathway enrichment pipeline, please open an issue on GitHub or contact the maintainer. 