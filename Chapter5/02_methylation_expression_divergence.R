# --- thesis_scripts/Chapter5/02_methylation_expression_divergence.R ---
# Purpose: Analyzes the divergence (e.g., standard deviation, variance, range)
#          of methylation and expression levels across different conditions
#          (e.g., tissues, developmental stages, possibly individuals).
#          Relates divergence metrics to gene features or functional annotations.
# Chapter Relevance: 5.3 DNA Methylation and Expression Divergence

# Load required libraries
library(dplyr)
library(tidyr) # For pivot_longer/gather
library(ggplot2)
library(matrixStats) # For rowSds, rowVars etc.
library(stats)
library(ggpubr) # For stat_cor

# --- Configuration & Input --- 

# Define input file paths (Replace with actual paths or loaded objects)
# Assumes input matrices/dataframes `methylation_matrix` and `expression_matrix` are loaded.
# - Rows should be genes (gene_id as rownames or first column).
# - Columns should be samples/conditions.
# Also assumes `gene_annotations` dataframe is loaded (needs gene_id and other features).

output_dir <- "results/chapter5/"
divergence_stats_file <- file.path(output_dir, "methylation_expression_divergence_stats.csv")
div_corr_plot_file <- file.path(output_dir, "methylation_vs_expression_divergence_correlation.pdf")
div_feature_plot_prefix <- file.path(output_dir, "divergence_vs_feature") # Prefix for multiple plots

# --- Data Preparation & Divergence Calculation --- 

# Check if input data exists (replace with your actual data objects)
if (!exists("methylation_matrix") || !exists("expression_matrix") || !exists("gene_annotations")) {
    stop("Input data objects (methylation_matrix, expression_matrix, gene_annotations) not found.")
}

# Ensure matrices have gene IDs accessible (e.g., as rownames)
# Convert to matrix if they are dataframes, preserving rownames
if (!is.matrix(methylation_matrix)) {
    if("gene_id" %in% colnames(methylation_matrix)){
        rownames(methylation_matrix) <- methylation_matrix$gene_id
        methylation_matrix$gene_id <- NULL
        methylation_matrix <- as.matrix(methylation_matrix)
    } else {
        stop("methylation_matrix needs gene IDs as rownames or a 'gene_id' column.")
    }
}
if (!is.matrix(expression_matrix)) {
     if("gene_id" %in% colnames(expression_matrix)){
        rownames(expression_matrix) <- expression_matrix$gene_id
        expression_matrix$gene_id <- NULL
        expression_matrix <- as.matrix(expression_matrix)
    } else {
        stop("expression_matrix needs gene IDs as rownames or a 'gene_id' column.")
    }
}

# Ensure annotations have gene_id
if (!"gene_id" %in% colnames(gene_annotations)) {
    stop("gene_annotations needs a 'gene_id' column.")
}

# Align genes between matrices and annotations
common_genes <- intersect(rownames(methylation_matrix), rownames(expression_matrix))
common_genes <- intersect(common_genes, gene_annotations$gene_id)

cat("Found", length(common_genes), "common genes across all inputs.\n")

methylation_matrix_aligned <- methylation_matrix[common_genes, , drop = FALSE]
expression_matrix_aligned <- expression_matrix[common_genes, , drop = FALSE]
gene_annotations_aligned <- gene_annotations %>% filter(gene_id %in% common_genes)

# Calculate divergence metrics (Standard Deviation chosen here, could use var, IQR, range etc.)
cat("Calculating divergence metrics (Standard Deviation)...\n")
meth_divergence <- rowSds(methylation_matrix_aligned, na.rm = TRUE)
expr_divergence <- rowSds(expression_matrix_aligned, na.rm = TRUE)

# Combine divergence stats with annotations
divergence_stats <- data.frame(
    gene_id = common_genes,
    methylation_sd = meth_divergence,
    expression_sd = expr_divergence
) %>% 
    left_join(gene_annotations_aligned, by = "gene_id")

# Handle potential NaN or Inf values resulting from SD calculation (e.g., constant values)
divergence_stats <- divergence_stats %>% 
    filter(is.finite(methylation_sd) & is.finite(expression_sd))

cat("Divergence calculated for", nrow(divergence_stats), "genes.\n")

# Save divergence stats
write.csv(divergence_stats, divergence_stats_file, row.names = FALSE)
cat("Saved divergence statistics to:", divergence_stats_file, "\n")

# --- Analysis & Visualization --- 

# 1. Correlation between Methylation Divergence and Expression Divergence
cat("\nAnalyzing correlation between methylation and expression divergence...\n")

cor_test_divergence <- cor.test(~ methylation_sd + expression_sd,
                                data = divergence_stats,
                                method = "spearman") # Use Spearman for potentially non-linear relationship
print(cor_test_divergence)

p_div_corr <- ggplot(divergence_stats, aes(x = methylation_sd, y = expression_sd)) +
    geom_point(alpha = 0.2, size = 0.8) +
    geom_smooth(method = "lm", se = FALSE, color = "red") + # Add linear fit line (optional)
    stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
    labs(
        title = "Methylation Divergence vs. Expression Divergence",
        x = "Methylation Standard Deviation (across samples)",
        y = "Expression Standard Deviation (across samples)"
    ) +
    # Optional: Log scale if data is highly skewed
    # scale_x_log10() + 
    # scale_y_log10() +
    theme_minimal()

print(p_div_corr)
ggsave(div_corr_plot_file, p_div_corr, width = 6, height = 5)
cat("Saved divergence correlation plot to:", div_corr_plot_file, "\n")


# 2. Relationship between Divergence and Gene Features
cat("\nAnalyzing divergence vs. gene features...\n")

# Example: Divergence vs. Mean Expression/Methylation Level (requires merging mean levels)
# Need to calculate mean levels first if not already in annotations
if (!all(c("mean_methylation", "mean_expression") %in% colnames(divergence_stats))) {
    cat("Calculating mean levels for plotting vs divergence...\n")
    mean_meth <- rowMeans(methylation_matrix_aligned, na.rm = TRUE)
    mean_expr <- rowMeans(expression_matrix_aligned, na.rm = TRUE)
    mean_levels <- data.frame(gene_id = common_genes, 
                              mean_methylation = mean_meth, 
                              mean_expression = mean_expr)
    divergence_stats <- left_join(divergence_stats, mean_levels, by="gene_id")
}

features_to_plot_vs_div <- c("mean_methylation", "mean_expression") # Add others like length, GC content if available

plot_list_div_feat <- list()
for (feature in features_to_plot_vs_div) {
    if (feature %in% colnames(divergence_stats)) {
        # Methylation divergence vs feature
        title_meth <- paste("Methylation Divergence vs.", feature)
        p_meth <- ggplot(divergence_stats, aes(x = !!sym(feature), y = methylation_sd)) +
                    geom_point(alpha = 0.1, size = 0.5) +
                    geom_smooth(method="loess", color="blue", se=FALSE) +
                    labs(title = title_meth, x = feature, y = "Methylation SD") + theme_minimal()
        plot_list_div_feat[[paste0("meth_",feature)]] <- p_meth
        ggsave(paste0(div_feature_plot_prefix, "_meth_vs_", feature, ".pdf"), p_meth, width=5, height=4)

        # Expression divergence vs feature
        title_expr <- paste("Expression Divergence vs.", feature)
        p_expr <- ggplot(divergence_stats, aes(x = !!sym(feature), y = expression_sd)) +
                    geom_point(alpha = 0.1, size = 0.5) +
                    geom_smooth(method="loess", color="red", se=FALSE) + 
                    labs(title = title_expr, x = feature, y = "Expression SD") + theme_minimal()
        plot_list_div_feat[[paste0("expr_",feature)]] <- p_expr
        ggsave(paste0(div_feature_plot_prefix, "_expr_vs_", feature, ".pdf"), p_expr, width=5, height=4)
    } else {
        warning("Feature '", feature, "' not found in divergence_stats.")
    }
}
cat("Saved divergence vs feature plots with prefix:", div_feature_plot_prefix, "\n")

# Display some plots (optional)
# if (length(plot_list_div_feat) > 0) { print(plot_list_div_feat[[1]]) }

cat("--- 02_methylation_expression_divergence.R finished ---\n") 