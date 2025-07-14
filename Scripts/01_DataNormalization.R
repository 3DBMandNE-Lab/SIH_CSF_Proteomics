# =============================================================================
# Proteomics Data Normalization and Quality Assessment
# =============================================================================
# 
# Description: Advanced normalization and quality assessment pipeline for proteomics data
# - Applies Variance Stabilizing Normalization (VSN) for intensity-dependent variance
# - Performs Z-score scaling and centering for downstream analysis
# - Implements PCA-based outlier detection with configurable thresholds
# - Generates comprehensive quality assessment visualizations
# - Provides detailed logging and progress tracking
#
# Author: Kevin Joseph
# Date: 07/2024
# 
# Dependencies:
# - vsn: Variance stabilizing normalization
# - matrixStats: Matrix statistics for variance calculation
# - factoextra: PCA visualization
# - pheatmap: Heatmap generation
# - tidyverse: Data manipulation and visualization
#
# Input: mat_clean (from 00_PreprocessData.R)
# Output: Normalized data matrix and quality assessment plots
# =============================================================================

# =============================================================================
# 1. SETUP AND CONFIGURATION
# =============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(vsn)
  library(matrixStats)
  library(factoextra)
  library(pheatmap)
  library(tidyverse)
})

# Configuration parameters
norm_config <- list(
  # Normalization settings
  vsn_method = "vsn2",           # VSN method to use
  
  # Scaling parameters
  center_data = TRUE,             # Center the data
  scale_data = TRUE,              # Scale the data
  
  # Outlier detection
  outlier_threshold = 9,          # Mahalanobis distance threshold (>3 SD)
  pca_components = 2,             # Number of PCs for outlier detection
  
  # Visualization settings
  top_variable_proteins = 50,     # Number of proteins for heatmap
  heatmap_cell_size = 8,         # Cell size for heatmap
  pca_ellipses = TRUE,           # Add confidence ellipses to PCA
  
  # Output settings
  verbose = TRUE,
  save_plots = FALSE,
  plot_format = "pdf"
)

# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================

#' Print formatted status message for normalization pipeline
#' @param message Message to print
#' @param level Message level (info, warning, error)
print_norm_status <- function(message, level = "info") {
  if (!norm_config$verbose) return()
  
  timestamp <- format(Sys.time(), "%H:%M:%S")
  prefix <- switch(level,
                   info = "NORM_INFO",
                   warning = "NORM_WARNING", 
                   error = "NORM_ERROR")
  
  cat(sprintf("[%s] %s: %s\n", timestamp, prefix, message))
}

#' Validate input data for normalization
#' @param data_matrix Input data matrix
#' @return Logical indicating if data is valid
validate_normalization_input <- function(data_matrix) {
  if (is.null(data_matrix)) {
    stop("Input data matrix is NULL")
  }
  
  if (!is.matrix(data_matrix) && !is.data.frame(data_matrix)) {
    stop("Input must be a matrix or data frame")
  }
  
  if (nrow(data_matrix) == 0 || ncol(data_matrix) == 0) {
    stop("Input data matrix is empty")
  }
  
  if (any(is.infinite(data_matrix))) {
    warning("Infinite values detected in data matrix")
  }
  
  print_norm_status(sprintf("Input validation passed: %d proteins × %d samples", 
                           nrow(data_matrix), ncol(data_matrix)))
  return(TRUE)
}

#' Calculate and display normalization statistics
#' @param before_data Data before normalization
#' @param after_data Data after normalization
#' @param step_name Name of the normalization step
analyze_normalization_effect <- function(before_data, after_data, step_name) {
  before_stats <- c(
    mean = mean(before_data, na.rm = TRUE),
    sd = sd(before_data, na.rm = TRUE),
    var = var(as.vector(before_data), na.rm = TRUE)
  )
  
  after_stats <- c(
    mean = mean(after_data, na.rm = TRUE),
    sd = sd(after_data, na.rm = TRUE),
    var = var(as.vector(after_data), na.rm = TRUE)
  )
  
  print_norm_status(sprintf("%s - Before: mean=%.3f, sd=%.3f, var=%.3f", 
                           step_name, before_stats["mean"], before_stats["sd"], before_stats["var"]))
  print_norm_status(sprintf("%s - After:  mean=%.3f, sd=%.3f, var=%.3f", 
                           step_name, after_stats["mean"], after_stats["sd"], after_stats["var"]))
}

# =============================================================================
# 3. DATA VALIDATION
# =============================================================================

print_norm_status("Starting proteomics data normalization pipeline")

# Check if input data exists
if (!exists("mat_clean")) {
  stop("Input data 'mat_clean' not found. Please run 00_PreprocessData.R first.")
}

if (!exists("sample_anno")) {
  stop("Sample annotations 'sample_anno' not found. Please run 00_PreprocessData.R first.")
}

# Validate input data
validate_normalization_input(mat_clean)

# =============================================================================
# 4. VARIANCE STABILIZING NORMALIZATION (VSN)
# =============================================================================

print_norm_status("=== VARIANCE STABILIZING NORMALIZATION ===")

print_norm_status("Applying VSN for intensity-dependent variance stabilization")
vsn_fit <- vsn2(mat_clean)
mat_norm <- predict(vsn_fit, mat_clean)

# Analyze VSN effect
analyze_normalization_effect(mat_clean, mat_norm, "VSN")

print_norm_status(sprintf("VSN completed: %d proteins × %d samples", 
                         nrow(mat_norm), ncol(mat_norm)))

# =============================================================================
# 5. SCALING AND CENTERING
# =============================================================================

print_norm_status("=== SCALING AND CENTERING ===")

print_norm_status("Applying Z-score scaling and centering")
mat_scaled <- t(scale(t(mat_norm), 
                      center = norm_config$center_data, 
                      scale = norm_config$scale_data))

# Analyze scaling effect
analyze_normalization_effect(mat_norm, mat_scaled, "Z-score scaling")

print_norm_status(sprintf("Scaling completed: %d proteins × %d samples", 
                         nrow(mat_scaled), ncol(mat_scaled)))

# =============================================================================
# 6. OUTLIER DETECTION
# =============================================================================

print_norm_status("=== OUTLIER DETECTION ===")

print_norm_status("Performing PCA-based outlier detection")

# Perform PCA for outlier detection
pca_tmp <- prcomp(t(mat_scaled), center = FALSE, scale. = FALSE)
pc_scores <- pca_tmp$x[, 1:norm_config$pca_components, drop = FALSE]

# Calculate Mahalanobis distance
dist_z <- rowSums(scale(pc_scores)^2)
outliers <- which(dist_z > norm_config$outlier_threshold)

print_norm_status(sprintf("Outlier detection: %d samples analyzed, threshold = %.1f", 
                         length(dist_z), norm_config$outlier_threshold))

if (length(outliers) > 0) {
  outlier_samples <- rownames(pc_scores)[outliers]
  print_norm_status(sprintf("Outlier samples detected (%d): %s", 
                           length(outliers), paste(outlier_samples, collapse = ", ")), "warning")
  
  # Remove outliers
  mat_scaled <- mat_scaled[, -outliers, drop = FALSE]
  sample_anno <- sample_anno[-match(outlier_samples, rownames(sample_anno)), , drop = FALSE]
  
  print_norm_status(sprintf("Outliers removed: %d proteins × %d samples remaining", 
                           nrow(mat_scaled), ncol(mat_scaled)))
} else {
  print_norm_status("No outliers detected; all samples retained")
}

# =============================================================================
# 7. QUALITY ASSESSMENT VISUALIZATIONS
# =============================================================================

print_norm_status("=== QUALITY ASSESSMENT VISUALIZATIONS ===")

# PCA for quality assessment
print_norm_status("Generating PCA plot for quality assessment")
pca_res <- prcomp(t(mat_scaled), center = TRUE, scale. = FALSE)

pca_plot <- fviz_pca_ind(pca_res,
                          geom.ind = "point",
                          col.ind = sample_anno$Group,
                          addEllipses = norm_config$pca_ellipses,
                          legend.title = "Group",
                          title = "PCA: Sample Quality Assessment") +
  coord_fixed() + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

print(pca_plot)

# Heatmap of top variable proteins
print_norm_status("Generating heatmap of top variable proteins")
top_proteins <- names(sort(rowVars(mat_scaled), decreasing = TRUE))[1:norm_config$top_variable_proteins]

heatmap_plot <- pheatmap(mat_scaled[top_proteins, ],
                         annotation_col = sample_anno["Group", drop = FALSE],
                         scale = "row",
                         show_rownames = TRUE,
                         border_color = "white",
                         cellwidth = norm_config$heatmap_cell_size,
                         cellheight = norm_config$heatmap_cell_size,
                         show_colnames = FALSE,
                         main = sprintf("Top %d Variable Proteins", norm_config$top_variable_proteins),
                         fontsize_row = 8,
                         fontsize_col = 8)

print(heatmap_plot)

# =============================================================================
# 8. FINAL SUMMARY AND OUTPUT
# =============================================================================

print_norm_status("=== NORMALIZATION SUMMARY ===")

# Calculate quality metrics
total_variance_explained <- sum(pca_res$sdev^2 / sum(pca_res$sdev^2) * 100)
pc1_variance <- pca_res$sdev[1]^2 / sum(pca_res$sdev^2) * 100
pc2_variance <- pca_res$sdev[2]^2 / sum(pca_res$sdev^2) * 100

print_norm_status(sprintf("Final normalized dataset: %d proteins × %d samples", 
                         nrow(mat_scaled), ncol(mat_scaled)))
print_norm_status(sprintf("Total variance explained by all PCs: %.1f%%", total_variance_explained))
print_norm_status(sprintf("PC1 explains %.1f%% of variance", pc1_variance))
print_norm_status(sprintf("PC2 explains %.1f%% of variance", pc2_variance))

# Data quality metrics
print_norm_status("Data quality metrics:")
print_norm_status(sprintf("- Mean intensity: %.3f", mean(mat_scaled, na.rm = TRUE)))
print_norm_status(sprintf("- Standard deviation: %.3f", sd(mat_scaled, na.rm = TRUE)))
print_norm_status(sprintf("- Coefficient of variation: %.3f", 
                         sd(mat_scaled, na.rm = TRUE) / mean(mat_scaled, na.rm = TRUE)))

# Save plots if requested
if (norm_config$save_plots) {
  plot_dir <- "plots"
  if (!dir.exists(plot_dir)) dir.create(plot_dir)
  
  ggsave(file.path(plot_dir, "pca_quality_assessment.pdf"), pca_plot, 
          width = 10, height = 8, dpi = 300)
  
  pdf(file.path(plot_dir, "heatmap_top_proteins.pdf"), width = 12, height = 10)
  print(heatmap_plot)
  dev.off()
  
  print_norm_status(sprintf("Plots saved to %s/ directory", plot_dir))
}

print_norm_status("Normalization pipeline completed successfully")

# =============================================================================
# 9. OUTPUT OBJECTS
# =============================================================================
# 
# The following objects are now available for downstream analysis:
# - mat_scaled: Normalized and scaled intensity matrix (proteins × samples)
# - mat_norm: VSN-normalized matrix (before scaling)
# - vsn_fit: VSN fit object for future predictions
# - pca_res: PCA results for quality assessment
# - sample_anno: Updated sample annotations (outliers removed)
# - norm_config: Configuration parameters used
#
# Next steps typically include:
# - Differential expression analysis
# - Pathway enrichment analysis
# - Clustering analysis
# - Machine learning applications
# =============================================================================

