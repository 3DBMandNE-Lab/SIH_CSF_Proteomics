# =============================================================================
# PLS-DA Analysis for Proteomics Data Classification
# =============================================================================
# 
# This script performs Partial Least Squares Discriminant Analysis (PLS-DA) 
# on preprocessed proteomics data to classify samples and identify key 
# discriminatory proteins.
#
# Input Requirements:
#   - mat_scaled: normalized and scaled protein matrix (proteins × samples)
#   - sample_anno: data.frame with rownames = sample names, column Group
#
# Outputs:
#   - PLS-DA model with optimal number of components
#   - Sample classification scores and visualizations
#   - VIP (Variable Importance in Projection) scores for feature selection
#   - Sparse PLS-DA for feature selection
#   - GO enrichment analysis of discriminatory proteins
#   - Performance evaluation and cross-validation results
#
# Dependencies: mixOmics, ggplot2, dplyr, clusterProfiler, org.Hs.eg.db
#
# Author: Kevin Joseph
# Date: 07/2025
# Version: 1.0
# =============================================================================

# =============================================================================
# Configuration and Setup
# =============================================================================

# Configuration parameters
config <- list(
  # Input/Output paths
  input_dir = "data/processed",
  output_dir = "results/plsda",
  
  # PLS-DA parameters
  n_components = 2,              # Number of PLS components
  validation_folds = 5,          # Number of folds for cross-validation
  n_repeats = 10,               # Number of repeats for performance evaluation
  test_keepX = c(5, 10, 15, 20), # Test values for sparse PLS-DA
  
  # Feature selection parameters
  top_vip_count = 20,           # Number of top VIP proteins to analyze
  top_loading_count = 20,       # Number of top loading proteins to analyze
  
  # Visualization parameters
  plot_width = 12,
  plot_height = 10,
  dpi = 300,
  
  # Statistical parameters
  pvalue_cutoff = 0.05,
  qvalue_cutoff = 0.05,
  
  # Random seed for reproducibility
  seed = 42
)

# Create output directory
if (!dir.exists(config$output_dir)) {
  dir.create(config$output_dir, recursive = TRUE)
}

# =============================================================================
# Load Required Libraries
# =============================================================================

required_packages <- c(
  "mixOmics", "ggplot2", "dplyr", "clusterProfiler", "org.Hs.eg.db"
)

# Install missing packages
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages)
}

# Load all packages
invisible(lapply(required_packages, library, character.only = TRUE))

# =============================================================================
# Data Loading and Validation
# =============================================================================

cat("Loading and validating input data...\n")

# Load normalized data
tryCatch({
  load(file.path(config$input_dir, "normalized_data.RData"))
  cat("✓ Normalized data loaded successfully\n")
}, error = function(e) {
  stop("Failed to load normalized data: ", e$message)
})

# Validate input data
if (!exists("mat_scaled") || !is.matrix(mat_scaled)) {
  stop("mat_scaled matrix not found or invalid")
}

if (!exists("sample_anno") || !is.data.frame(sample_anno)) {
  stop("sample_anno data.frame not found or invalid")
}

if (!"Group" %in% colnames(sample_anno)) {
  stop("sample_anno must contain 'Group' column")
}

# Check for sufficient samples per group
group_counts <- table(sample_anno$Group)
if (any(group_counts < 3)) {
  stop("Each group must have at least 3 samples for PLS-DA")
}

cat("✓ Input data validation passed\n")
cat("  - Matrix dimensions:", nrow(mat_scaled), "×", ncol(mat_scaled), "\n")
cat("  - Sample groups:", paste(names(group_counts), "(", group_counts, ")", collapse = ", "), "\n")

# =============================================================================
# 1. Data Preparation
# =============================================================================

cat("\n=== Step 1: Preparing data for PLS-DA ===\n")

# Prepare predictor matrix (samples × proteins)
X <- t(mat_scaled)  # Transpose to samples × features
Y <- factor(sample_anno$Group, levels = unique(sample_anno$Group))

cat("✓ Data preparation completed\n")
cat("  - X matrix dimensions:", nrow(X), "samples ×", ncol(X), "proteins\n")
cat("  - Y factor levels:", paste(levels(Y), collapse = ", "), "\n")

# =============================================================================
# 2. PLS-DA Model Fitting
# =============================================================================

cat("\n=== Step 2: Fitting PLS-DA model ===\n")

set.seed(config$seed)

# Run PLS-DA with specified number of components
cat("Running PLS-DA with", config$n_components, "components...\n")

plsda_res <- plsda(X, Y, ncomp = config$n_components)

# Extract variance explained
var_explained <- round(plsda_res$explained_variance$X * 100, 1)
cat("✓ PLS-DA model fitted successfully\n")
cat("  - Variance explained by components:", paste(var_explained, "%", collapse = ", "), "\n")

# =============================================================================
# 3. Sample Score Visualization
# =============================================================================

cat("\n=== Step 3: Visualizing sample scores ===\n")

# Extract sample scores
scores <- as.data.frame(plsda_res$variates$X)
colnames(scores) <- paste0("PLS", 1:config$n_components)
scores$Group <- Y

# Create score plot with confidence ellipses
score_plot <- ggplot(scores, aes(x = PLS1, y = PLS2, color = Group, shape = Group)) +
  # Filled confidence ellipses
  stat_ellipse(aes(fill = Group), geom = "polygon",
               level = 0.95, alpha = 0.2, color = NA) +
  # Points on top
  geom_point(size = 4, stroke = 1.2) +
  # Custom scales
  scale_color_manual(values = c(NoLeak = "#1f78b4", Leak = "#e31a1c")) +
  scale_fill_manual(values = c(NoLeak = "#a6cee3", Leak = "#fb9a99")) +
  scale_shape_manual(values = c(NoLeak = 16, Leak = 17)) +
  # Axis labels with % variance
  labs(
    title = "PLS-DA: Sample Classification",
    subtitle = paste0("PLS1 (", var_explained[1], "%) vs PLS2 (", var_explained[2], "%)"),
    x = paste0("PLS1 (", var_explained[1], "% variance)"),
    y = paste0("PLS2 (", var_explained[2], "% variance)"),
    color = NULL, shape = NULL, fill = NULL
  ) +
  theme_bw() +
  coord_fixed() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.major = element_line(color = "grey90", size = 0.4),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.margin = margin(b = -10)
  )

ggsave(
  file.path(config$output_dir, "plsda_sample_scores.pdf"),
  score_plot,
  width = config$plot_width, height = config$plot_height, dpi = config$dpi
)

cat("✓ Sample score visualization completed\n")

# =============================================================================
# 4. VIP Score Analysis
# =============================================================================

cat("\n=== Step 4: Analyzing VIP scores ===\n")

# Calculate VIP scores for all components
vip_scores <- vip(plsda_res)

# Extract VIP for component 1 and rank
vip1 <- vip_scores[, 1]
top_vip <- sort(vip1, decreasing = TRUE)[1:config$top_vip_count]
top_vip_genes <- names(top_vip)

# Create VIP score data frame
df_vip <- data.frame(
  Gene = factor(top_vip_genes, levels = top_vip_genes),
  VIP = top_vip
)

# Create VIP score barplot
vip_plot <- ggplot(df_vip, aes(x = Gene, y = VIP)) +
  geom_col(fill = "#2c7bb6") +
  coord_flip() +
  labs(
    title = paste0("Top ", config$top_vip_count, " Proteins by VIP (PLS1)"),
    y = "VIP score",
    x = ""
  ) +
  theme_bw(base_size = 14)

ggsave(
  file.path(config$output_dir, "vip_scores_barplot.pdf"),
  vip_plot,
  width = config$plot_width, height = config$plot_height, dpi = config$dpi
)

# Create VIP score lollipop plot
df_vip_ordered <- df_vip %>%
  arrange(VIP) %>%
  mutate(Gene = factor(Gene, levels = Gene))

vip_lollipop <- ggplot(df_vip_ordered, aes(x = Gene, y = VIP)) +
  geom_segment(aes(x = Gene, xend = Gene, y = 0, yend = VIP),
               color = "gray70", size = 0.5) +
  geom_point(color = "#2c7bb6", size = 4) +
  coord_flip() +
  labs(
    title = paste0("Top ", config$top_vip_count, " VIP Scores"),
    x = NULL,
    y = "VIP score"
  ) +
  theme_bw()

ggsave(
  file.path(config$output_dir, "vip_scores_lollipop.pdf"),
  vip_lollipop,
  width = config$plot_width, height = config$plot_height, dpi = config$dpi
)

cat("✓ VIP score analysis completed\n")
cat("  - Top VIP protein:", names(top_vip)[1], "(VIP =", round(top_vip[1], 3), ")\n")

# =============================================================================
# 5. Loading Analysis
# =============================================================================

cat("\n=== Step 5: Analyzing loadings ===\n")

# Extract loadings for PLS1
load1 <- plsda_res$loadings$X[, 1]
df_load <- data.frame(
  Gene = names(load1),
  Loading = load1
) %>%
  arrange(desc(abs(Loading))) %>%
  dplyr::slice(1:config$top_loading_count)

# Create loading plot
loading_plot <- ggplot(df_load, aes(x = reorder(Gene, Loading), y = Loading)) +
  geom_col(fill = ifelse(df_load$Loading > 0, "#d73027", "#4575b4")) +
  coord_flip() +
  labs(
    title = paste0("Top ", config$top_loading_count, " Loadings on PLS1"),
    y = "Loading weight",
    x = ""
  ) +
  theme_bw(base_size = 14)

ggsave(
  file.path(config$output_dir, "loadings_barplot.pdf"),
  loading_plot,
  width = config$plot_width, height = config$plot_height, dpi = config$dpi
)

cat("✓ Loading analysis completed\n")
cat("  - Top positive loading:", df_load$Gene[1], "(loading =", round(df_load$Loading[1], 3), ")\n")

# =============================================================================
# 6. Performance Evaluation
# =============================================================================

cat("\n=== Step 6: Evaluating model performance ===\n")

set.seed(config$seed)

# Perform cross-validation
cat("Running cross-validation with", config$validation_folds, "folds and", config$n_repeats, "repeats...\n")

perf_res <- perf(plsda_res, 
                 validation = "Mfold", 
                 folds = config$validation_folds, 
                 nrepeat = config$n_repeats,
                 progressBar = FALSE)

# Create performance plot
perf_plot <- plot(perf_res, col = c("#1b9e77", "#d95f02", "firebrick"))

ggsave(
  file.path(config$output_dir, "plsda_performance.pdf"),
  perf_plot,
  width = config$plot_width, height = config$plot_height, dpi = config$dpi
)

cat("✓ Performance evaluation completed\n")

# =============================================================================
# 7. Sparse PLS-DA for Feature Selection
# =============================================================================

cat("\n=== Step 7: Sparse PLS-DA for feature selection ===\n")

# Tune sparsity parameters
cat("Tuning sparse PLS-DA parameters...\n")

tune_res <- tune.splsda(
  X, Y, ncomp = config$n_components,
  validation = 'Mfold', 
  folds = config$validation_folds,
  test.keepX = config$test_keepX,
  nrepeat = 5, 
  progressBar = FALSE
)

# Extract optimal keepX
keepX <- tune_res$choice.keepX
cat("  - Optimal keepX values:", paste(keepX, collapse = ", "), "\n")

# Run sparse PLS-DA
splsda_res <- splsda(X, Y, ncomp = config$n_components, keepX = keepX)

# Create sparse PLS-DA plot
splsda_plot <- plotIndiv(splsda_res, comp = c(1, 2), group = Y,
                         ind.names = FALSE, ellipse = TRUE,
                         legend = TRUE, title = "Sparse PLS-DA")

ggsave(
  file.path(config$output_dir, "sparse_plsda_plot.pdf"),
  splsda_plot,
  width = config$plot_width, height = config$plot_height, dpi = config$dpi
)

cat("✓ Sparse PLS-DA completed\n")

# =============================================================================
# 8. GO Enrichment Analysis
# =============================================================================

cat("\n=== Step 8: GO enrichment analysis of discriminatory proteins ===\n")

# Combine top VIP and Loading hits
top_drivers <- unique(c(top_vip_genes, df_load$Gene))

cat("  - Total discriminatory proteins:", length(top_drivers), "\n")

# Perform GO enrichment
ego_drivers <- enrichGO(
  gene = top_drivers,
  keyType = "SYMBOL",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = config$pvalue_cutoff
)

if (nrow(ego_drivers@result) > 0) {
  # Create dotplot
  dot_plot <- dotplot(ego_drivers, showCategory = 10) +
    labs(title = "GO BP Enrichment of PLS-DA Drivers")
  
  ggsave(
    file.path(config$output_dir, "go_enrichment_dotplot.pdf"),
    dot_plot,
    width = config$plot_width, height = config$plot_height, dpi = config$dpi
  )
  
  # Create network plot
  network_plot <- cnetplot(ego_drivers) +
    theme_bw() +
    coord_fixed()
  
  ggsave(
    file.path(config$output_dir, "go_enrichment_network.pdf"),
    network_plot,
    width = config$plot_width, height = config$plot_height, dpi = config$dpi
  )
  
  # Create lollipop plot of top enriched terms
  df_topFE <- ego_drivers@result %>%
    arrange(desc(FoldEnrichment)) %>%
    dplyr::slice(1:20) %>%
    mutate(
      Term = factor(Description, levels = rev(Description)),
      GeneRatioNum = as.numeric(sub("/.*", "", GeneRatio)) /
        as.numeric(sub(".*/", "", GeneRatio))
    )
  
  lollipop_plot <- ggplot(df_topFE, aes(x = FoldEnrichment, y = Term)) +
    geom_segment(aes(x = 0, xend = FoldEnrichment, y = Term, yend = Term),
                 color = "gray80") +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_viridis_c(option = "plasma", direction = -1, 
                          name = "adj. P", limits = c(0, 0.05)) +
    scale_size(range = c(3, 7), name = "Count") +
    labs(
      title = "Top GO BP Terms: PLS-DA Drivers",
      x = "Fold Enrichment",
      y = NULL
    ) +
    theme_bw()
  
  ggsave(
    file.path(config$output_dir, "go_enrichment_lollipop.pdf"),
    lollipop_plot,
    width = config$plot_width, height = config$plot_height, dpi = config$dpi
  )
  
  cat("✓ GO enrichment analysis completed\n")
  cat("  - Enriched terms found:", nrow(ego_drivers@result), "\n")
} else {
  cat("Warning: No significant GO enrichment found\n")
}

# =============================================================================
# 9. Save Results and Summary
# =============================================================================

cat("\n=== Step 9: Saving results and summary ===\n")

# Save PLS-DA results
save(plsda_res, splsda_res, perf_res, tune_res, vip_scores, top_drivers,
     file = file.path(config$output_dir, "plsda_results.RData"))

# Save configuration
saveRDS(config, file.path(config$output_dir, "plsda_config.rds"))

# Save top discriminatory proteins
top_proteins_df <- data.frame(
  Gene = top_drivers,
  VIP_Score = vip1[top_drivers],
  Loading_Score = load1[top_drivers],
  Rank_VIP = match(top_drivers, names(sort(vip1, decreasing = TRUE))),
  Rank_Loading = match(top_drivers, names(sort(abs(load1), decreasing = TRUE)))
)

write.csv(top_proteins_df, 
          file.path(config$output_dir, "top_discriminatory_proteins.csv"), 
          row.names = FALSE)

# Generate summary report
sink(file.path(config$output_dir, "plsda_summary.txt"))
cat("PLS-DA Analysis Summary\n")
cat("=====================\n\n")
cat("Input Data:\n")
cat("- Matrix dimensions:", nrow(mat_scaled), "×", ncol(mat_scaled), "\n")
cat("- Sample groups:", paste(names(group_counts), "(", group_counts, ")", collapse = ", "), "\n\n")
cat("PLS-DA Results:\n")
cat("- Number of components:", config$n_components, "\n")
cat("- Variance explained:", paste(var_explained, "%", collapse = ", "), "\n")
cat("- Top discriminatory proteins:", length(top_drivers), "\n\n")
cat("Performance Metrics:\n")
cat("- Cross-validation folds:", config$validation_folds, "\n")
cat("- Number of repeats:", config$n_repeats, "\n\n")
if (nrow(ego_drivers@result) > 0) {
  cat("GO Enrichment Results:\n")
  cat("- Enriched terms:", nrow(ego_drivers@result), "\n")
  cat("- Top enriched term:", ego_drivers@result$Description[1], "\n")
}
cat("\nAnalysis completed successfully!\n")
sink()

cat("✓ PLS-DA analysis completed successfully!\n")
cat("Results saved to:", config$output_dir, "\n")
cat("Files generated:\n")
cat("- plsda_results.RData: Main results\n")
cat("- top_discriminatory_proteins.csv: Top discriminatory proteins\n")
cat("- Various visualization PDFs\n")
cat("- plsda_summary.txt: Analysis summary\n")

# =============================================================================
# End of Script
# =============================================================================
