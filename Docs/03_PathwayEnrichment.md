# 03_PathwayEnrichment.R

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

 - Data Validation
 - Gene Set Preparation
 - Over-Representation Analysis (ORA)
 - Gene Set Enrichment Analysis (GSEA)
 - Visualization Generation

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

## Next Steps

After running this pathway enrichment pipeline, the results are ready for:

1. **KEGG pathway analysis**: Additional pathway database analysis
2. **Reactome pathway analysis**: Curated pathway database analysis
3. **Protein-protein interaction networks**: Network analysis of enriched proteins
4. **Functional annotation**: Detailed functional analysis of enriched terms
5. **Biomarker prioritization**: Selection of candidate biomarkers
6. **Meta-analysis**: Integration with other datasets or studies
