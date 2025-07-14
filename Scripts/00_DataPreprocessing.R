
# =============================================================================
# Proteomics Data Preprocessing Pipeline
# =============================================================================
# 
# Description: Comprehensive preprocessing pipeline for proteomics data analysis
# - Reads raw proteomics intensity data from Excel files
# - Performs aggressive protein filtering based on missingness thresholds
# - Implements minimal sample removal strategy (only extreme cases)
# - Handles data validation and error checking
# - Provides detailed logging of preprocessing steps
#
# Author: Kevin Joseph
# Date: 07/2024
# 
# Dependencies:
# - readxl: Excel file reading
# - limma: Linear models for microarray data
# - vsn: Variance stabilizing normalization
# - tidyverse: Data manipulation and visualization
# - pheatmap: Heatmap generation
# - factoextra: PCA visualization
#
# Input: Excel file with "Raw_Intensity" and "SampleInfo" sheets
# Output: Cleaned intensity matrix ready for downstream analysis
# =============================================================================

# =============================================================================
# 1. SETUP AND CONFIGURATION
# =============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(readxl)
  library(limma)
  library(vsn)
  library(tidyverse)
  library(pheatmap)
  library(factoextra)
})

# Configuration parameters
config <- list(
  # File paths (modify these for your setup)
  input_file = "data/Raw_Proteomics.xlsx",  # Relative path for GitHub
  raw_sheet = "Raw_Intensity",
  info_sheet = "SampleInfo",
  
  # Filtering thresholds
  gene_missing_threshold = 0.20,  # Remove genes with >20% missing values
  sample_missing_threshold = 0.50, # Remove samples with >50% missing values
  
  # Output settings
  verbose = TRUE,
  save_intermediate = FALSE
)

# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================

#' Validate file existence and structure
#' @param file_path Path to the Excel file
#' @param required_sheets Vector of required sheet names
#' @return Logical indicating if file is valid
validate_input_file <- function(file_path, required_sheets) {
  if (!file.exists(file_path)) {
    stop("Input file not found: ", file_path)
  }
  
  # Check if file is readable Excel
  tryCatch({
    sheets <- excel_sheets(file_path)
    missing_sheets <- setdiff(required_sheets, sheets)
    if (length(missing_sheets) > 0) {
      stop("Missing required sheets: ", paste(missing_sheets, collapse = ", "))
    }
    TRUE
  }, error = function(e) {
    stop("Error reading Excel file: ", e$message)
  })
}

#' Print formatted status message
#' @param message Message to print
#' @param level Message level (info, warning, error)
print_status <- function(message, level = "info") {
  if (!config$verbose) return()
  
  timestamp <- format(Sys.time(), "%H:%M:%S")
  prefix <- switch(level,
                   info = "INFO",
                   warning = "WARNING", 
                   error = "ERROR")
  
  cat(sprintf("[%s] %s: %s\n", timestamp, prefix, message))
}

#' Calculate and display missingness statistics
#' @param data_matrix Input data matrix
#' @param data_type Type of data ("genes" or "samples")
analyze_missingness <- function(data_matrix, data_type) {
  if (data_type == "genes") {
    missing_rates <- rowMeans(is.na(data_matrix))
    summary_stats <- summary(missing_rates)
  } else {
    missing_rates <- colMeans(is.na(data_matrix))
    summary_stats <- summary(missing_rates)
  }
  
  print_status(sprintf("%s missingness summary:", data_type), "info")
  print(summary_stats)
  
  return(missing_rates)
}

# =============================================================================
# 3. DATA IMPORT AND VALIDATION
# =============================================================================

print_status("Starting proteomics data preprocessing pipeline")

# Validate input file
required_sheets <- c(config$raw_sheet, config$info_sheet)
validate_input_file(config$input_file, required_sheets)

# Import data
print_status("Importing data from Excel file")
raw_data <- read_excel(config$input_file, sheet = config$raw_sheet)
sample_anno <- read_excel(config$input_file, sheet = config$info_sheet) %>%
  column_to_rownames("SampleID")

print_status(sprintf("Imported %d rows × %d columns of raw data", 
                    nrow(raw_data), ncol(raw_data)))

# Construct intensity matrix
print_status("Constructing intensity matrix")
mat_raw <- raw_data %>%
  rename(Gene = `Gene name`) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  column_to_rownames("Gene") %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

print_status(sprintf("Intensity matrix: %d proteins × %d samples", 
                    nrow(mat_raw), ncol(mat_raw)))

# =============================================================================
# 4. MISSINGNESS FILTERING
# =============================================================================

print_status("=== MISSINGNESS FILTERING ===")

# Step 1: Analyze and filter genes based on missingness
print_status("Step 1: Gene-level missingness filtering")
gene_missing_rates <- analyze_missingness(mat_raw, "genes")

gene_thresh <- config$gene_missing_threshold
keep_genes <- gene_missing_rates <= gene_thresh
removed_genes <- sum(!keep_genes)

print_status(sprintf("Removing %d genes with >%.1f%% missing values", 
                    removed_genes, gene_thresh * 100))

mat_filt <- mat_raw[keep_genes, , drop = FALSE]
print_status(sprintf("After gene filtering: %d proteins × %d samples", 
                    nrow(mat_filt), ncol(mat_filt)))

# Step 2: Analyze and optionally filter samples
print_status("Step 2: Sample-level missingness analysis")
sample_missing_rates <- analyze_missingness(mat_filt, "samples")

samp_thresh <- config$sample_missing_threshold
remove_samps <- names(sample_missing_rates)[sample_missing_rates > samp_thresh]

if (length(remove_samps) > 0) {
  print_status(sprintf("Removing %d samples with >%.1f%% missing values: %s", 
                      length(remove_samps), samp_thresh * 100,
                      paste(remove_samps, collapse = ", ")))
  
  mat_filt <- mat_filt[, !(colnames(mat_filt) %in% remove_samps), drop = FALSE]
  sample_anno <- sample_anno[!rownames(sample_anno) %in% remove_samps, , drop = FALSE]
} else {
  print_status(sprintf("No samples exceed %.1f%% missing; keeping all samples", 
                      samp_thresh * 100))
}

# Step 3: Final cleanup - remove any remaining genes with NA values
print_status("Step 3: Final cleanup - removing genes with any remaining NA values")
genes_with_na <- rowSums(is.na(mat_filt)) > 0
if (any(genes_with_na)) {
  print_status(sprintf("Removing %d genes with remaining NA values", sum(genes_with_na)))
}

mat_clean <- mat_filt[rowSums(is.na(mat_filt)) == 0, , drop = FALSE]

# =============================================================================
# 5. FINAL SUMMARY AND OUTPUT
# =============================================================================

print_status("=== PREPROCESSING SUMMARY ===")
print_status(sprintf("Final dataset: %d proteins × %d samples", 
                    nrow(mat_clean), ncol(mat_clean)))

# Calculate data retention statistics
gene_retention <- nrow(mat_clean) / nrow(mat_raw) * 100
sample_retention <- ncol(mat_clean) / ncol(mat_raw) * 100

print_status(sprintf("Gene retention: %.1f%% (%d/%d)", 
                    gene_retention, nrow(mat_clean), nrow(mat_raw)))
print_status(sprintf("Sample retention: %.1f%% (%d/%d)", 
                    sample_retention, ncol(mat_clean), ncol(mat_raw)))

# Data quality metrics
print_status("Data quality metrics:")
print_status(sprintf("- Total missing values: %d", sum(is.na(mat_clean))))
print_status(sprintf("- Data completeness: %.2f%%", 
                    (1 - sum(is.na(mat_clean)) / length(mat_clean)) * 100))

# Save session info for reproducibility
if (config$save_intermediate) {
  session_info_file <- "session_info_preprocessing.txt"
  sink(session_info_file)
  print(sessionInfo())
  sink()
  print_status(sprintf("Session info saved to: %s", session_info_file))
}

print_status("Preprocessing pipeline completed successfully")

# =============================================================================
# 6. OUTPUT OBJECTS
# =============================================================================
# 
# The following objects are now available for downstream analysis:
# - mat_clean: Cleaned intensity matrix (proteins × samples)
# - sample_anno: Sample annotation data frame
# - config: Configuration parameters used
# - gene_missing_rates: Original gene missingness rates
# - sample_missing_rates: Sample missingness rates after gene filtering
#
# Next steps typically include:
# - Normalization (e.g., vsn, quantile normalization)
# - Log transformation
# - Outlier detection
# - Principal component analysis (PCA)
# =============================================================================

