# =============================================================================
# Proteomics Pathway Enrichment Analysis Pipeline
# =============================================================================
# 
# This script performs comprehensive pathway enrichment analysis of differential
# proteomics results using Gene Ontology and Gene Set Enrichment Analysis.
#
# Input Requirements:
#   - de_res: Differential abundance results from differential analysis pipeline
#   - mat_scaled: Normalized intensity matrix for background analysis
#   - sample_anno: Sample annotation data frame
#
# Outputs:
#   - ego_up: GO enrichment results for up-regulated genes
#   - ego_down: GO enrichment results for down-regulated genes
#   - gsea_res: GSEA results for all genes
#   - gene_sets: Prepared gene lists for analysis
#   - Comprehensive visualizations (volcano plots, network plots, heatmaps)
#   - Integrated enrichment dashboard
#   - Comprehensive logging of all analysis steps
#
# Dependencies: clusterProfiler, enrichplot, GOplot, ggplot2, patchwork, viridis, org.Hs.eg.db
#
# Author: Dr.-Ing. Kevin Joseph, Neurosurgery, Medical Center - University of Freiburg
# Date: 2024
# Version: 1.0
# =============================================================================

# =============================================================================
# 1. SETUP AND CONFIGURATION
# =============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(GOplot)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(viridis)
  library(ggrepel)
  library(pheatmap)
})

# Configuration parameters
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

# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================

#' Print formatted status message for enrichment analysis
#' @param message Message to print
#' @param level Message level (info, warning, error)
print_enrich_status <- function(message, level = "info") {
  if (!enrich_config$verbose) return()
  
  timestamp <- format(Sys.time(), "%H:%M:%S")
  prefix <- switch(level,
                   info = "ENRICH_INFO",
                   warning = "ENRICH_WARNING", 
                   error = "ENRICH_ERROR")
  
  cat(sprintf("[%s] %s: %s\n", timestamp, prefix, message))
}

#' Validate input data for enrichment analysis
#' @param differential_results Differential analysis results
#' @param data_matrix Input data matrix
#' @param sample_annotations Sample annotation data frame
#' @return Logical indicating if data is valid
validate_enrichment_input <- function(differential_results, data_matrix, sample_annotations) {
  if (is.null(differential_results)) {
    stop("Differential analysis results are NULL")
  }
  
  if (is.null(data_matrix)) {
    stop("Input data matrix is NULL")
  }
  
  if (is.null(sample_annotations)) {
    stop("Sample annotations are NULL")
  }
  
  if (!is.data.frame(differential_results) && !is.matrix(differential_results)) {
    stop("Differential results must be a data frame or matrix")
  }
  
  if (nrow(differential_results) == 0) {
    stop("Differential results are empty")
  }
  
  required_columns <- c("logFC", "adj.P.Val")
  missing_columns <- setdiff(required_columns, colnames(differential_results))
  if (length(missing_columns) > 0) {
    stop("Missing required columns in differential results: ", paste(missing_columns, collapse = ", "))
  }
  
  print_enrich_status(sprintf("Input validation passed: %d proteins in differential results", 
                             nrow(differential_results)))
  return(TRUE)
}

#' Prepare gene lists for enrichment analysis
#' @param differential_results Differential analysis results
#' @param config Configuration parameters
#' @return List containing up-regulated and down-regulated gene sets
prepare_gene_lists <- function(differential_results, config) {
  # Create named gene list for GSEA
  gene_list <- differential_results$logFC
  names(gene_list) <- rownames(differential_results)
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Prepare up and down gene sets for ORA
  up_genes <- rownames(differential_results)[
    differential_results$logFC > config$logfc_threshold & 
    differential_results$adj.P.Val < config$fdr_cutoff
  ]
  
  down_genes <- rownames(differential_results)[
    differential_results$logFC < -config$logfc_threshold & 
    differential_results$adj.P.Val < config$fdr_cutoff
  ]
  
  print_enrich_status(sprintf("Gene sets prepared: %d up-regulated, %d down-regulated", 
                             length(up_genes), length(down_genes)))
  
  return(list(
    gene_list = gene_list,
    up_genes = up_genes,
    down_genes = down_genes
  ))
}

#' Create enhanced volcano plot
#' @param differential_results Differential analysis results
#' @param config Configuration parameters
#' @return ggplot object
create_enhanced_volcano <- function(differential_results, config) {
  # Prepare data for plotting
  plot_data <- differential_results %>%
    mutate(
      Gene = rownames(.),
      SignifUp = adj.P.Val < config$fdr_cutoff & logFC > config$logfc_threshold,
      SignifDn = adj.P.Val < config$fdr_cutoff & logFC < -config$logfc_threshold
    )
  
  # Create enhanced volcano plot
  volcano_plot <- ggplot(plot_data, aes(x = logFC, y = -log10(P.Value))) +
    
    # Background shading for enriched regions
    annotate("rect", xmin = -Inf, xmax = -config$logfc_threshold, ymin = 0, ymax = Inf,
             fill = "#B3E2CD", alpha = 0.3) +  # soft green
    annotate("rect", xmin = config$logfc_threshold, xmax = Inf, ymin = 0, ymax = Inf,
             fill = "#FDCDAC", alpha = 0.3) +  # soft red
    
    # Point layer
    geom_point(aes(color = adj.P.Val), size = config$volcano_point_size, alpha = config$volcano_alpha) +
    scale_color_viridis_c(option = "plasma", direction = -1, limits = c(0, config$fdr_cutoff),
                          oob = scales::squish, name = "adj. P.Val") +
    
    # Labels for significant genes
    geom_text_repel(
      data = filter(plot_data, abs(logFC) > config$logfc_threshold),
      aes(label = Gene),
      size = 3,
      segment.size = 0.3,
      segment.alpha = 0.6,
      box.padding = 0.4,
      point.padding = 0.3
    ) +
    
    # Cutoff lines
    geom_vline(xintercept = c(-config$logfc_threshold, config$logfc_threshold), 
               linetype = "dashed", color = "grey60") +
    geom_hline(yintercept = -log10(config$pvalue_cutoff), 
               linetype = "dashed", color = "grey60") +
    
    # Axis and title
    scale_x_continuous(breaks = seq(-2, 2, 0.5), 
                      limits = c(min(plot_data$logFC)-0.2, max(plot_data$logFC)+0.2)) +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.1))) +
    labs(
      title = "Volcano Plot: Differential Protein Abundance",
      subtitle = "Shaded: left = down-regulated, right = up-regulated",
      x = expression(log[2]~Fold~Change),
      y = expression(-log[10]~p-value)
    ) +
    
    # Clean theme
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 10)),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90"),
      legend.position = "right"
    ) +
    coord_fixed()
  
  return(volcano_plot)
}

# =============================================================================
# 3. DATA VALIDATION
# =============================================================================

print_enrich_status("Starting proteomics pathway enrichment analysis")

# Check if input data exists
if (!exists("de_res")) {
  stop("Differential analysis results 'de_res' not found. Please run 02_DifferentialTesting.R first.")
}

if (!exists("mat_scaled")) {
  stop("Normalized data 'mat_scaled' not found. Please run 01_DataNormalization.R first.")
}

if (!exists("sample_anno")) {
  stop("Sample annotations 'sample_anno' not found. Please run 01_DataNormalization.R first.")
}

# Validate input data
validate_enrichment_input(de_res, mat_scaled, sample_anno)

# =============================================================================
# 4. GENE SET PREPARATION
# =============================================================================

print_enrich_status("=== GENE SET PREPARATION ===")

# Prepare gene lists for enrichment analysis
gene_sets <- prepare_gene_lists(de_res, enrich_config)

# =============================================================================
# 5. OVER-REPRESENTATION ANALYSIS (ORA)
# =============================================================================

print_enrich_status("=== OVER-REPRESENTATION ANALYSIS ===")

# Load organism database
print_enrich_status("Loading organism annotation database")
org_db <- get(enrich_config$organism)

# Perform GO enrichment for up-regulated genes
if (length(gene_sets$up_genes) > 0) {
  print_enrich_status("Performing GO enrichment for up-regulated genes")
  ego_up <- enrichGO(
    gene = gene_sets$up_genes,
    OrgDb = org_db,
    ont = enrich_config$ontology,
    pAdjustMethod = "BH",
    keyType = enrich_config$key_type,
    pvalueCutoff = enrich_config$pvalue_cutoff
  )
  
  if (nrow(ego_up) > 0) {
    print_enrich_status(sprintf("Found %d enriched terms for up-regulated genes", nrow(ego_up)))
  } else {
    print_enrich_status("No enriched terms found for up-regulated genes")
  }
} else {
  print_enrich_status("No up-regulated genes found for enrichment analysis")
  ego_up <- NULL
}

# Perform GO enrichment for down-regulated genes
if (length(gene_sets$down_genes) > 0) {
  print_enrich_status("Performing GO enrichment for down-regulated genes")
  ego_down <- enrichGO(
    gene = gene_sets$down_genes,
    OrgDb = org_db,
    ont = enrich_config$ontology,
    pAdjustMethod = "BH",
    keyType = enrich_config$key_type,
    pvalueCutoff = enrich_config$pvalue_cutoff
  )
  
  if (nrow(ego_down) > 0) {
    print_enrich_status(sprintf("Found %d enriched terms for down-regulated genes", nrow(ego_down)))
  } else {
    print_enrich_status("No enriched terms found for down-regulated genes")
  }
} else {
  print_enrich_status("No down-regulated genes found for enrichment analysis")
  ego_down <- NULL
}

# =============================================================================
# 6. GENE SET ENRICHMENT ANALYSIS (GSEA)
# =============================================================================

print_enrich_status("=== GENE SET ENRICHMENT ANALYSIS ===")

print_enrich_status("Performing GSEA on GO Biological Process")
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

if (nrow(gsea_res) > 0) {
  print_enrich_status(sprintf("GSEA completed: %d enriched gene sets found", nrow(gsea_res)))
} else {
  print_enrich_status("No enriched gene sets found in GSEA")
}

# =============================================================================
# 7. VISUALIZATION
# =============================================================================

print_enrich_status("=== VISUALIZATION ===")

# Create enhanced volcano plot
print_enrich_status("Generating enhanced volcano plot")
volcano_plot <- create_enhanced_volcano(de_res, enrich_config)
print(volcano_plot)

# Create heatmap of top differential proteins
print_enrich_status("Generating heatmap of top differential proteins")
top_n <- min(20, nrow(de_res))
top_proteins <- rownames(de_res)[1:top_n]

heatmap_plot <- pheatmap(
  mat_scaled[top_proteins, ],
  annotation_col = sample_anno["Group", drop = FALSE],
  border_color = "white",
  cellwidth = enrich_config$heatmap_cell_size,
  cellheight = enrich_config$heatmap_cell_size,
  show_colnames = FALSE,
  main = sprintf("Top %d Differential Proteins", top_n),
  fontsize_row = 8,
  fontsize_col = 8
)

print(heatmap_plot)

# Create enrichment visualizations
if (!is.null(gsea_res) && nrow(gsea_res) > 0) {
  print_enrich_status("Generating GSEA dotplot")
  gsea_dotplot <- dotplot(gsea_res,
                          showCategory = enrich_config$dotplot_show,
                          title = "GSEA of GO Biological Processes",
                          font.size = 12,
                          label_format = 30) +
    scale_color_gradient2(
      low = "steelblue",
      mid = "grey80",
      high = "firebrick",
      midpoint = 0,
      name = "NES"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.text.y = element_text(size = 12),
      legend.position = "right"
    )
  
  print(gsea_dotplot)
  
  # Create enrichment map
  print_enrich_status("Generating enrichment map")
  sim <- pairwise_termsim(gsea_res)
  enrich_map <- emapplot(sim, showCategory = enrich_config$network_show, layout = "fr") +
    labs(title = "Enrichment Map of GO BP") +
    theme_bw()
  
  print(enrich_map)
  
  # Create gene-concept network
  print_enrich_status("Generating gene-concept network")
  cnet_plot <- cnetplot(gsea_res, 
                        foldChange = gene_sets$gene_list, 
                        showCategory = enrich_config$network_show,
                        circular = FALSE, 
                        colorEdge = TRUE) +
    labs(title = "Gene-Concept Network") +
    theme_bw(base_size = 12) +
    coord_fixed()
  
  print(cnet_plot)
  
  # Create heatplot
  print_enrich_status("Generating core enrichment heatmap")
  heat_plot <- heatplot(gsea_res, 
                        foldChange = gene_sets$gene_list, 
                        showCategory = enrich_config$heatmap_show) +
    labs(title = "Core Enrichment Heatmap") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  print(heat_plot)
  
  # Create upset plot
  print_enrich_status("Generating upset plot")
  upset_plot <- upsetplot(gsea_res, showCategory = enrich_config$upset_show) +
    labs(title = "Gene Overlap Across Top Terms") +
    theme_bw(base_size = 12)
  
  print(upset_plot)
}

# =============================================================================
# 8. INTEGRATED DASHBOARD
# =============================================================================

print_enrich_status("=== INTEGRATED DASHBOARD ===")

if (!is.null(gsea_res) && nrow(gsea_res) > 0) {
  print_enrich_status("Creating integrated enrichment dashboard")
  
  # Assemble multi-panel dashboard
  dashboard <- (gsea_dotplot + enrich_map) /
    (cnet_plot + heat_plot) /
    (upset_plot + volcano_plot) +
    plot_annotation(
      title = "Integrated Pathway Enrichment Dashboard",
      subtitle = "GSEA + ORA + Visualization Summary",
      theme = theme(plot.title = element_text(face = "bold", size = 16))
    )
  
  print(dashboard)
}

# =============================================================================
# 9. RESULTS EXPORT
# =============================================================================

print_enrich_status("=== RESULTS EXPORT ===")

if (enrich_config$save_results) {
  # Save ORA results
  if (!is.null(ego_up) && nrow(ego_up) > 0) {
    write.csv(as.data.frame(ego_up), 
              file = "GO_Enrichment_Up_Regulated.csv", 
              row.names = FALSE)
    print_enrich_status("Up-regulated GO enrichment results saved")
  }
  
  if (!is.null(ego_down) && nrow(ego_down) > 0) {
    write.csv(as.data.frame(ego_down), 
              file = "GO_Enrichment_Down_Regulated.csv", 
              row.names = FALSE)
    print_enrich_status("Down-regulated GO enrichment results saved")
  }
  
  # Save GSEA results
  if (!is.null(gsea_res) && nrow(gsea_res) > 0) {
    write.csv(as.data.frame(gsea_res), 
              file = "GSEA_Results.csv", 
              row.names = FALSE)
    print_enrich_status("GSEA results saved")
  }
  
  # Save plots if requested
  if (enrich_config$save_plots) {
    plot_dir <- "plots"
    if (!dir.exists(plot_dir)) dir.create(plot_dir)
    
    ggsave(file.path(plot_dir, "volcano_plot_enhanced.pdf"), 
            volcano_plot, width = 10, height = 8, dpi = 300)
    
    if (!is.null(gsea_res) && nrow(gsea_res) > 0) {
      ggsave(file.path(plot_dir, "gsea_dotplot.pdf"), 
              gsea_dotplot, width = 12, height = 10, dpi = 300)
      ggsave(file.path(plot_dir, "enrichment_dashboard.pdf"), 
              dashboard, width = 16, height = 18, dpi = 300)
    }
    
    print_enrich_status(sprintf("Plots saved to %s/ directory", plot_dir))
  }
}

# =============================================================================
# 10. FINAL SUMMARY
# =============================================================================

print_enrich_status("=== ENRICHMENT ANALYSIS SUMMARY ===")

# Summary statistics
total_genes <- length(gene_sets$gene_list)
up_genes_count <- length(gene_sets$up_genes)
down_genes_count <- length(gene_sets$down_genes)

print_enrich_status(sprintf("Analysis completed for %d total genes", total_genes))
print_enrich_status(sprintf("Up-regulated genes: %d", up_genes_count))
print_enrich_status(sprintf("Down-regulated genes: %d", down_genes_count))

if (!is.null(ego_up) && nrow(ego_up) > 0) {
  print_enrich_status(sprintf("Up-regulated enriched terms: %d", nrow(ego_up)))
}

if (!is.null(ego_down) && nrow(ego_down) > 0) {
  print_enrich_status(sprintf("Down-regulated enriched terms: %d", nrow(ego_down)))
}

if (!is.null(gsea_res) && nrow(gsea_res) > 0) {
  print_enrich_status(sprintf("GSEA enriched gene sets: %d", nrow(gsea_res)))
}

print_enrich_status("Pathway enrichment analysis pipeline completed successfully")

# =============================================================================
# 11. OUTPUT OBJECTS
# =============================================================================
# 
# The following objects are now available for downstream analysis:
# - ego_up: GO enrichment results for up-regulated genes
# - ego_down: GO enrichment results for down-regulated genes
# - gsea_res: GSEA results for all genes
# - gene_sets: Prepared gene lists for analysis
# - volcano_plot: Enhanced volcano plot
# - heatmap_plot: Heatmap of top differential proteins
# - enrich_config: Configuration parameters used
#
# Next steps typically include:
# - KEGG pathway analysis
# - Reactome pathway analysis
# - Protein-protein interaction networks
# - Functional annotation
# - Biomarker prioritization
# =============================================================================