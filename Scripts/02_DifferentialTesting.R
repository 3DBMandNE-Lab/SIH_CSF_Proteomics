# =============================================================================
# Proteomics Differential Abundance Analysis Pipeline
# =============================================================================
# 
# This script performs comprehensive differential abundance analysis of normalized
# proteomics data using empirical Bayes moderated t-tests and multiple testing correction.
#
# Input Requirements:
#   - mat_scaled: Normalized and scaled intensity matrix from normalization pipeline
#   - sample_anno: Sample annotation data frame with Group column
#
# Outputs:
#   - de_res: Differential abundance results (limma topTable output)
#   - fit2: Fitted limma model with empirical Bayes results
#   - volcano_plot: Enhanced volcano plot with background shading
#   - heatmap_plot: Heatmap of top differential proteins
#   - Statistical summaries and quality assessment metrics
#   - Comprehensive logging of all analysis steps
#
# Dependencies: limma, EnhancedVolcano, ggplot2, ggrepel, dplyr, pheatmap
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
  library(limma)
  library(EnhancedVolcano)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(pheatmap)
})

# Configuration parameters
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

# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================

#' Print formatted status message for differential analysis
#' @param message Message to print
#' @param level Message level (info, warning, error)
print_diff_status <- function(message, level = "info") {
  if (!diff_config$verbose) return()
  
  timestamp <- format(Sys.time(), "%H:%M:%S")
  prefix <- switch(level,
                   info = "DIFF_INFO",
                   warning = "DIFF_WARNING", 
                   error = "DIFF_ERROR")
  
  cat(sprintf("[%s] %s: %s\n", timestamp, prefix, message))
}

#' Validate input data for differential analysis
#' @param data_matrix Input data matrix
#' @param sample_annotations Sample annotation data frame
#' @return Logical indicating if data is valid
validate_differential_input <- function(data_matrix, sample_annotations) {
  if (is.null(data_matrix)) {
    stop("Input data matrix is NULL")
  }
  
  if (is.null(sample_annotations)) {
    stop("Sample annotations are NULL")
  }
  
  if (!is.matrix(data_matrix) && !is.data.frame(data_matrix)) {
    stop("Input must be a matrix or data frame")
  }
  
  if (nrow(data_matrix) == 0 || ncol(data_matrix) == 0) {
    stop("Input data matrix is empty")
  }
  
  if (!"Group" %in% colnames(sample_annotations)) {
    stop("Sample annotations must contain 'Group' column")
  }
  
  if (length(unique(sample_annotations$Group)) < 2) {
    stop("At least 2 groups required for differential analysis")
  }
  
  print_diff_status(sprintf("Input validation passed: %d proteins × %d samples", 
                           nrow(data_matrix), ncol(data_matrix)))
  return(TRUE)
}

#' Create design matrix for differential analysis
#' @param sample_annotations Sample annotation data frame
#' @param reference_group Reference group name
#' @return Design matrix for limma
create_design_matrix <- function(sample_annotations, reference_group) {
  # Ensure reference group is first level
  sample_annotations$Group <- relevel(
    factor(sample_annotations$Group, 
           levels = c(reference_group, 
                     unique(sample_annotations$Group[sample_annotations$Group != reference_group]))), 
    ref = reference_group
  )
  
  # Create design matrix
  design <- model.matrix(~ Group, data = sample_annotations)
  
  print_diff_status(sprintf("Design matrix created with %d groups", ncol(design)))
  print_diff_status(sprintf("Reference group: %s", reference_group))
  
  return(list(design = design, sample_annotations = sample_annotations))
}

#' Generate enhanced volcano plot
#' @param results Differential analysis results
#' @param config Configuration parameters
#' @return ggplot object
create_volcano_plot <- function(results, config) {
  # Prepare results for plotting
  plot_data <- results %>%
    mutate(
      Gene = rownames(.),
      SignifUp = adj.P.Val < config$fdr_threshold & logFC > config$logfc_threshold,
      SignifDn = adj.P.Val < config$fdr_threshold & logFC < -config$logfc_threshold
    )
  
  # Create volcano plot with background shading
  volcano_plot <- ggplot(plot_data, aes(x = logFC, y = -log10(P.Value), color = adj.P.Val)) +
    
    # Background shading for enriched regions
    annotate("rect", xmin = -Inf, xmax = -config$logfc_threshold,
             ymin = -Inf, ymax = Inf,
             fill = "lightgreen", alpha = 0.2) +
    annotate("rect", xmin = config$logfc_threshold, xmax = Inf,
             ymin = -Inf, ymax = Inf,
             fill = "lightcoral", alpha = 0.2) +
    
    # Points
    geom_point(alpha = config$volcano_point_alpha) +
    scale_color_gradient(low = "red", high = "blue",
                         limits = c(0, config$fdr_threshold), 
                         oob = scales::squish) +
    
    # Labels for significant proteins
    geom_text_repel(
      data = filter(plot_data, abs(logFC) > config$logfc_threshold),
      aes(label = Gene),
      size = config$volcano_label_size,
      box.padding = 0.3,
      point.padding = 0.2,
      max.overlaps = config$volcano_max_labels
    ) +
    
    # Cutoff lines
    geom_vline(xintercept = c(-config$logfc_threshold, config$logfc_threshold), 
               linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(config$pvalue_threshold), 
               linetype = "dashed", color = "grey50") +
    
    # Labels and theme
    labs(
      x = sprintf("log2 Fold Change (%s vs %s)", 
                  config$test_group, config$reference_group),
      y = "-log10(p-value)",
      color = "adj. P.Val",
      title = "Volcano Plot of Differential Abundance",
      subtitle = sprintf("Shaded regions: left = %s-enriched, right = %s-enriched", 
                        config$reference_group, config$test_group)
    ) +
    theme_bw() +
    coord_fixed()
  
  return(volcano_plot)
}

# =============================================================================
# 3. DATA VALIDATION
# =============================================================================

print_diff_status("Starting proteomics differential abundance analysis")

# Check if input data exists
if (!exists("mat_scaled")) {
  stop("Input data 'mat_scaled' not found. Please run 01_DataNormalization.R first.")
}

if (!exists("sample_anno")) {
  stop("Sample annotations 'sample_anno' not found. Please run 01_DataNormalization.R first.")
}

# Validate input data
validate_differential_input(mat_scaled, sample_anno)

# =============================================================================
# 4. EXPERIMENTAL DESIGN
# =============================================================================

print_diff_status("=== EXPERIMENTAL DESIGN ===")

# Create design matrix
design_result <- create_design_matrix(sample_anno, diff_config$reference_group)
design <- design_result$design
sample_anno <- design_result$sample_annotations

print_diff_status(sprintf("Design matrix columns: %s", paste(colnames(design), collapse = ", ")))

# =============================================================================
# 5. DIFFERENTIAL ABUNDANCE ANALYSIS
# =============================================================================

print_diff_status("=== DIFFERENTIAL ABUNDANCE ANALYSIS ===")

# Step 1: Fit linear model
print_diff_status("Fitting linear model with limma")
fit <- lmFit(mat_scaled, design)

# Step 2: Apply empirical Bayes shrinkage
print_diff_status("Applying empirical Bayes shrinkage")
fit2 <- eBayes(fit)

# Step 3: Extract results
print_diff_status("Extracting differential abundance results")
de_res <- topTable(
  fit2,
  coef = paste0("Group", diff_config$test_group),
  number = Inf,
  adjust.method = "BH"
)

# =============================================================================
# 6. RESULTS SUMMARY
# =============================================================================

print_diff_status("=== RESULTS SUMMARY ===")

# Calculate summary statistics
total_proteins <- nrow(de_res)
significant_proteins <- sum(de_res$adj.P.Val < diff_config$fdr_threshold)
upregulated <- sum(de_res$adj.P.Val < diff_config$fdr_threshold & de_res$logFC > diff_config$logfc_threshold)
downregulated <- sum(de_res$adj.P.Val < diff_config$fdr_threshold & de_res$logFC < -diff_config$logfc_threshold)

print_diff_status(sprintf("Total proteins analyzed: %d", total_proteins))
print_diff_status(sprintf("Significant proteins (FDR < %.3f): %d", diff_config$fdr_threshold, significant_proteins))
print_diff_status(sprintf("Upregulated proteins: %d", upregulated))
print_diff_status(sprintf("Downregulated proteins: %d", downregulated))

# Show top proteins
if (significant_proteins > 0) {
  top_proteins <- head(rownames(de_res), 10)
  print_diff_status(sprintf("Top significant proteins: %s", paste(top_proteins, collapse = ", ")))
} else {
  print_diff_status("No significant proteins found at current thresholds")
}

# =============================================================================
# 7. VISUALIZATION
# =============================================================================

print_diff_status("=== VISUALIZATION ===")

# Generate volcano plot
print_diff_status("Generating volcano plot")
volcano_plot <- create_volcano_plot(de_res, diff_config)
print(volcano_plot)

# Generate heatmap of top differential proteins
print_diff_status("Generating heatmap of top differential proteins")
top_n <- min(diff_config$heatmap_top_n, nrow(de_res))
top_proteins <- rownames(de_res)[1:top_n]

heatmap_plot <- pheatmap(
  mat_scaled[top_proteins, ],
  annotation_col = sample_anno["Group", drop = FALSE],
  show_rownames = TRUE,
  border_color = "white",
  cellwidth = diff_config$heatmap_cell_size,
  cellheight = diff_config$heatmap_cell_size,
  show_colnames = FALSE,
  main = sprintf("Top %d Differential Proteins", top_n),
  fontsize_row = 8,
  fontsize_col = 8
)

print(heatmap_plot)

# =============================================================================
# 8. RESULTS EXPORT
# =============================================================================

print_diff_status("=== RESULTS EXPORT ===")

if (diff_config$save_results) {
  # Prepare results for export
  export_results <- de_res %>%
    mutate(
      Gene = rownames(.),
      Comparison = diff_config$comparison_name,
      Significant = adj.P.Val < diff_config$fdr_threshold,
      Upregulated = adj.P.Val < diff_config$fdr_threshold & logFC > diff_config$logfc_threshold,
      Downregulated = adj.P.Val < diff_config$fdr_threshold & logFC < -diff_config$logfc_threshold
    )
  
  # Save results
  output_file <- sprintf("Differential_Proteins_%s.csv", diff_config$comparison_name)
  write.csv(export_results, file = output_file, row.names = FALSE)
  
  print_diff_status(sprintf("Results saved to: %s", output_file))
  
  # Save plots if requested
  if (diff_config$save_plots) {
    plot_dir <- "plots"
    if (!dir.exists(plot_dir)) dir.create(plot_dir)
    
    ggsave(file.path(plot_dir, sprintf("volcano_plot_%s.pdf", diff_config$comparison_name)), 
            volcano_plot, width = 10, height = 8, dpi = 300)
    
    pdf(file.path(plot_dir, sprintf("heatmap_top_proteins_%s.pdf", diff_config$comparison_name)), 
        width = 12, height = 10)
    print(heatmap_plot)
    dev.off()
    
    print_diff_status(sprintf("Plots saved to %s/ directory", plot_dir))
  }
}

# =============================================================================
# 9. FINAL SUMMARY
# =============================================================================

print_diff_status("=== DIFFERENTIAL ANALYSIS SUMMARY ===")

print_diff_status(sprintf("Analysis completed for %s vs %s comparison", 
                         diff_config$test_group, diff_config$reference_group))
print_diff_status(sprintf("Statistical thresholds: FDR < %.3f, |logFC| > %.1f", 
                         diff_config$fdr_threshold, diff_config$logfc_threshold))
print_diff_status(sprintf("Quality metrics:"))
print_diff_status(sprintf("- Total proteins: %d", total_proteins))
print_diff_status(sprintf("- Significant proteins: %d (%.1f%%)", 
                         significant_proteins, significant_proteins/total_proteins*100))
print_diff_status(sprintf("- Upregulated: %d", upregulated))
print_diff_status(sprintf("- Downregulated: %d", downregulated))

print_diff_status("Differential analysis pipeline completed successfully")

# =============================================================================
# 10. OUTPUT OBJECTS
# =============================================================================
# 
# The following objects are now available for downstream analysis:
# - de_res: Differential abundance results (limma topTable output)
# - fit2: Fitted limma model with empirical Bayes results
# - volcano_plot: Enhanced volcano plot
# - heatmap_plot: Heatmap of top differential proteins
# - diff_config: Configuration parameters used
#
# Next steps typically include:
# - Pathway enrichment analysis
# - Gene ontology analysis
# - Protein-protein interaction networks
# - Biomarker validation
# - Functional annotation
# =============================================================================