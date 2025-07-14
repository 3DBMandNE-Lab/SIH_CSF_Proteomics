# =============================================================================
# NMF-Based Module Generation from Cleaned Proteomics Matrix
# =============================================================================
# 
# This script performs Non-negative Matrix Factorization (NMF) analysis on 
# preprocessed proteomics data to identify co-expression modules and their 
# biological significance.
#
# Input Requirements:
#   - mat_clean: proteins × samples numeric matrix (no NAs, non-negative)
#   - sample_anno: data.frame with rownames = sample names, column Group
#   - de_res: differential expression results (optional, for integration)
#
# Outputs:
#   - Optimal rank estimation via consensus clustering
#   - NMF factorization at chosen rank (W and H matrices)
#   - Gene-to-module assignments with loading scores
#   - Comprehensive visualizations of module metagenes and gene loadings
#   - GO enrichment analysis for each module
#   - Quality control metrics and statistical comparisons
#
# Dependencies: NMF, pheatmap, dplyr, ggplot2, clusterProfiler, org.Hs.eg.db,
#              patchwork, ggridges, ggalluvial, cluster, ggpubr, tidytext
#
# Author: [Your Name]
# Date: [Current Date]
# Version: 1.0
# =============================================================================

# =============================================================================
# Configuration and Setup
# =============================================================================

# Configuration parameters
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

# Create output directory
if (!dir.exists(config$output_dir)) {
  dir.create(config$output_dir, recursive = TRUE)
}

# =============================================================================
# Load Required Libraries
# =============================================================================

required_packages <- c(
  "NMF", "pheatmap", "dplyr", "ggplot2", "patchwork",
  "clusterProfiler", "org.Hs.eg.db", "ggridges", 
  "ggalluvial", "cluster", "ggpubr", "tidytext"
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

# Load preprocessed data
tryCatch({
  load(file.path(config$input_dir, "preprocessed_data.RData"))
  cat("✓ Preprocessed data loaded successfully\n")
}, error = function(e) {
  stop("Failed to load preprocessed data: ", e$message)
})

# Validate input data
if (!exists("mat_clean") || !is.matrix(mat_clean)) {
  stop("mat_clean matrix not found or invalid")
}

if (!exists("sample_anno") || !is.data.frame(sample_anno)) {
  stop("sample_anno data.frame not found or invalid")
}

if (any(is.na(mat_clean))) {
  stop("mat_clean contains NA values - please run preprocessing first")
}

if (any(mat_clean < 0)) {
  stop("mat_clean contains negative values - NMF requires non-negative data")
}

cat("✓ Input data validation passed\n")
cat("  - Matrix dimensions:", nrow(mat_clean), "×", ncol(mat_clean), "\n")
cat("  - Sample groups:", paste(unique(sample_anno$Group), collapse = ", "), "\n")

# =============================================================================
# 1. Optimal Rank Estimation via Consensus Clustering
# =============================================================================

cat("\n=== Step 1: Estimating optimal number of modules ===\n")

set.seed(config$seed)

# Run consensus clustering across different ranks
cat("Running consensus clustering across ranks", min(config$rank_range), "to", max(config$rank_range), "...\n")

consensus_est <- nmf(
  mat_clean, 
  rank = config$rank_range, 
  nrun = config$n_runs_consensus, 
  seed = config$seed, 
  .options = 'p'
)

# Plot consensus clustering results
consensus_plot <- plot(consensus_est)
ggsave(
  file.path(config$output_dir, "consensus_clustering.pdf"),
  consensus_plot,
  width = 10, height = 8, dpi = config$dpi
)

# Extract cophenetic coefficients for rank selection
cophenetic_coef <- sapply(consensus_est, cophenetic)
cat("Cophenetic coefficients by rank:\n")
for (i in seq_along(cophenetic_coef)) {
  cat("  Rank", config$rank_range[i], ":", round(cophenetic_coef[i], 4), "\n")
}

# Select optimal rank (you may need to adjust this based on your data)
# Rule of thumb: choose rank where cophenetic coefficient starts to drop
opt_rank <- config$rank_range[which.max(cophenetic_coef)]
cat("Selected optimal rank:", opt_rank, "\n")

# =============================================================================
# 2. NMF Factorization at Optimal Rank
# =============================================================================

cat("\n=== Step 2: Running NMF at optimal rank ===\n")

set.seed(config$seed + 1)

# Run NMF at optimal rank
cat("Running NMF with", config$n_runs_final, "runs at rank", opt_rank, "...\n")

nmf_res <- nmf(
  mat_clean, 
  rank = opt_rank, 
  nrun = config$n_runs_final, 
  seed = config$seed + 1, 
  .options = 'vP'
)

# Extract basis (W) and coefficient (H) matrices
W <- basis(nmf_res)    # genes × opt_rank
H <- coef(nmf_res)     # opt_rank × samples

# Add meaningful column/row names
colnames(W) <- paste0("Module", 1:opt_rank)
rownames(H) <- paste0("Module", 1:opt_rank)

cat("✓ NMF factorization completed\n")
cat("  - W matrix dimensions:", nrow(W), "×", ncol(W), "\n")
cat("  - H matrix dimensions:", nrow(H), "×", ncol(H), "\n")

# =============================================================================
# 3. Gene-to-Module Assignment
# =============================================================================

cat("\n=== Step 3: Assigning genes to modules ===\n")

# Assign each gene to its dominant module
gene_module <- apply(W, 1, function(x) which.max(x))
module_assign <- tibble(
  Gene   = rownames(W),
  Module = paste0("Module", gene_module),
  Loading = pmax(0, apply(W, 1, max))    # max loading per gene
)

# Summary statistics
module_summary <- module_assign %>%
  group_by(Module) %>%
  summarise(
    n_genes = n(),
    mean_loading = mean(Loading),
    max_loading = max(Loading),
    .groups = 'drop'
  )

cat("Module assignment summary:\n")
print(module_summary)

# Save gene-module assignments
write.csv(module_assign, 
          file.path(config$output_dir, "gene_module_assignments.csv"), 
          row.names = FALSE)
cat("✓ Gene-module assignments saved\n")

# =============================================================================
# 4. Module Metagene Visualization
# =============================================================================

cat("\n=== Step 4: Visualizing module metagenes ===\n")

# Create heatmap of module metagenes across samples (H matrix)
metagene_heatmap <- pheatmap(
  H,
  annotation_col = sample_anno["Group"],
  cluster_rows = TRUE,
  clustering_method = "average",
  border_color = "white",
  scale = "row",
  cellwidth = 10,
  cellheight = 10,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  main = "Module Metagenes Across Samples",
  fontsize = 10
)

ggsave(
  file.path(config$output_dir, "module_metagenes_heatmap.pdf"),
  metagene_heatmap,
  width = 10, height = 6, dpi = config$dpi
)

# =============================================================================
# 5. Top Gene Loadings Visualization
# =============================================================================

cat("\n=== Step 5: Visualizing top gene loadings ===\n")

# Select top genes by loading in each module
top_genes <- module_assign %>%
  group_by(Module) %>%
  slice_max(order_by = Loading, n = config$top_genes_per_module) %>%
  pull(Gene)

# Create annotation for heatmap
anno_row <- module_assign %>%
  filter(Gene %in% top_genes) %>%
  column_to_rownames("Gene") %>%
  dplyr::select(Module)

# Create heatmap of top gene loadings
gene_loadings_heatmap <- pheatmap(
  W[top_genes, ],
  annotation_row = anno_row,
  cluster_rows = TRUE,
  border_color = "white",
  cellwidth = 10,
  cluster_cols = TRUE,
  main = "Top Gene Loadings per Module"
)

ggsave(
  file.path(config$output_dir, "top_gene_loadings_heatmap.pdf"),
  gene_loadings_heatmap,
  width = 10, height = 8, dpi = config$dpi
)

# =============================================================================
# 6. Module Scores Visualization
# =============================================================================

cat("\n=== Step 6: Visualizing module scores ===\n")

# Prepare data for plotting
df_H <- as.data.frame(t(H)) %>%
  rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "Module", values_to = "Score") %>%
  left_join(sample_anno, by = c("Sample" = "rowname"))

# Create barplot of module scores per sample
module_scores_plot <- ggplot(df_H, aes(x = Sample, y = Score, fill = Module)) +
  geom_col(position = "dodge") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(
    title = "Module Scores per Sample", 
    y = "Metagene Score",
    x = "Sample"
  )

ggsave(
  file.path(config$output_dir, "module_scores_barplot.pdf"),
  module_scores_plot,
  width = 12, height = 6, dpi = config$dpi
)

# =============================================================================
# 7. GO Enrichment Analysis
# =============================================================================

cat("\n=== Step 7: Performing GO enrichment analysis ===\n")

# Build list of genes per module (top genes by loading)
module_genes <- module_assign %>%
  group_by(Module) %>%
  slice_max(Loading, n = config$max_genes_enrichment) %>%
  summarize(Genes = list(Gene)) %>%
  deframe()

# Perform enrichment for each module
ego_list <- lapply(module_genes, function(genes) {
  tryCatch({
    enrichGO(
      gene = genes,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = config$pvalue_cutoff,
      qvalueCutoff = config$qvalue_cutoff
    )
  }, error = function(e) {
    cat("Warning: GO enrichment failed for module with", length(genes), "genes\n")
    return(NULL)
  })
})

# Remove NULL results
ego_list <- ego_list[!sapply(ego_list, is.null)]

if (length(ego_list) > 0) {
  # Create network plots for each module
  plots <- imap(ego_list, function(ego, mod) {
    if (nrow(ego@result) > 0) {
      cnetplot(ego, showCategory = 8, node_label = "category") +
        scale_color_viridis_c(option = "magma", direction = -1) +
        labs(title = paste0(mod, " – Top BP")) +
        theme_bw() +
        coord_fixed()
    } else {
      ggplot() + 
        annotate("text", x = 0.5, y = 0.5, label = paste("No significant terms in", mod)) +
        theme_void()
    }
  })
  
  # Combine plots
  fig_modules <- wrap_plots(plots, ncol = 1)
  
  ggsave(
    file.path(config$output_dir, "module_GO_enrichment.pdf"),
    fig_modules,
    width = 6, height = 12, dpi = config$dpi
  )
  
  cat("✓ GO enrichment analysis completed\n")
} else {
  cat("Warning: No significant GO enrichment found\n")
}

# =============================================================================
# 8. Quality Control and Statistical Analysis
# =============================================================================

cat("\n=== Step 8: Quality control and statistical analysis ===\n")

# Module-wise sample scores boxplot with statistical tests
if (exists("de_res")) {
  # Integration with differential expression results
  de_genes <- rownames(de_res)[abs(de_res$logFC) > 1 & de_res$adj.P.Val < 0.05]
  
  if (length(de_genes) > 0) {
    df_alluv <- data.frame(
      Gene = de_genes,
      Direction = ifelse(de_res[de_genes, "logFC"] > 0, "Up", "Down"),
      Module = gene_module$Module[match(de_genes, gene_module$Gene)]
    )
    
    # Alluvial plot
    alluvial_plot <- ggplot(df_alluv, aes(axis1 = Direction, axis2 = Module, y = 1)) +
      geom_alluvium(aes(fill = Module), width = 1/12) +
      geom_stratum(width = 1/12, fill = "grey80") +
      geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
      labs(title = "Flow of Top DE Genes into Modules") +
      theme_bw()
    
    ggsave(
      file.path(config$output_dir, "DE_genes_alluvial.pdf"),
      alluvial_plot,
      width = 10, height = 6, dpi = config$dpi
    )
  }
}

# Statistical comparison of module scores between groups
my_comp <- list(c("NoLeak", "Leak"))  # Adjust based on your group names

module_comparison_plot <- ggplot(df_H, aes(x = Module, y = Score, fill = Group)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.6, outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              size = 1.5, alpha = 0.7) +
  stat_compare_means(
    comparisons = my_comp,
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = TRUE,
    bracket.size = 0.3,
    tip.length = 0.02,
    position = position_dodge(width = 0.75)
  ) +
  scale_fill_manual(values = c(NoLeak = "#0072B2", Leak = "#D55E00")) +
  labs(
    title = "Module Metagene Scores: Group Comparison",
    x = "NMF Module",
    y = "Metagene Score"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

ggsave(
  file.path(config$output_dir, "module_comparison_boxplot.pdf"),
  module_comparison_plot,
  width = 10, height = 6, dpi = config$dpi
)

# Silhouette analysis for module quality
distW <- dist(W)
cl <- as.numeric(factor(gene_module$Module))
sil <- silhouette(cl, distW)
df_sil <- data.frame(sil_width = sil[, 3], Module = gene_module$Module)

silhouette_plot <- ggplot(df_sil, aes(x = Module, y = sil_width, fill = Module)) +
  geom_boxplot(outlier.shape = NA) +
  xlab("") + 
  ylab("Silhouette Width") +
  labs(title = "Silhouette Widths of Gene Clustering") + 
  theme_bw()

ggsave(
  file.path(config$output_dir, "module_silhouette_analysis.pdf"),
  silhouette_plot,
  width = 8, height = 6, dpi = config$dpi
)

# =============================================================================
# 9. Save Results and Summary
# =============================================================================

cat("\n=== Step 9: Saving results and summary ===\n")

# Save NMF results
save(nmf_res, W, H, module_assign, consensus_est, opt_rank,
     file = file.path(config$output_dir, "nmf_results.RData"))

# Save configuration
saveRDS(config, file.path(config$output_dir, "nmf_config.rds"))

# Generate summary report
sink(file.path(config$output_dir, "nmf_summary.txt"))
cat("NMF Analysis Summary\n")
cat("==================\n\n")
cat("Input Data:\n")
cat("- Matrix dimensions:", nrow(mat_clean), "×", ncol(mat_clean), "\n")
cat("- Sample groups:", paste(unique(sample_anno$Group), collapse = ", "), "\n\n")
cat("NMF Results:\n")
cat("- Optimal rank:", opt_rank, "\n")
cat("- Number of modules:", opt_rank, "\n")
cat("- Total genes assigned:", nrow(module_assign), "\n\n")
cat("Module Summary:\n")
print(module_summary)
cat("\nAnalysis completed successfully!\n")
sink()

cat("✓ NMF analysis completed successfully!\n")
cat("Results saved to:", config$output_dir, "\n")
cat("Files generated:\n")
cat("- nmf_results.RData: Main results\n")
cat("- gene_module_assignments.csv: Gene-module assignments\n")
cat("- Various visualization PDFs\n")
cat("- nmf_summary.txt: Analysis summary\n")

# =============================================================================
# End of Script
# =============================================================================