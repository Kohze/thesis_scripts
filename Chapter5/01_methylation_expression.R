# --- thesis_scripts/Chapter5/01_methylation_expression.R ---
# Purpose: Analyzes the relationship between DNA methylation levels and gene expression.
#          Calculates correlations, performs statistical tests (e.g., t-tests)
#          between methylation groups (e.g., high/low) and expression levels,
#          and generates visualizations like scatter plots and box plots.
# Chapter Relevance: 5.2 DNA Methylation and Gene Expression

# Load required libraries
library(dplyr)
library(ggplot2)
library(ggpubr) # For stat_compare_means, stat_cor
library(stats) # For cor.test, t.test
library(tidyr) # For pivot_longer or gather
library(broom) # For tidying test results

# --- Configuration & Input --- 

# Define input file paths (Replace with actual paths or loaded objects)
# Assumes input objects `methylation_data` and `expression_data` are loaded.
# - methylation_data: Needs columns like gene_id, methylation_level (e.g., avg beta value)
# - expression_data: Needs columns like gene_id, expression_value (e.g., log2(TPM+1))
# - gene_annotations: Needs gene_id and potentially gene_biotype or other categories.

output_dir <- "results/chapter5/"
correlation_plot_file <- file.path(output_dir, "methylation_expression_correlation.pdf")
boxplot_file <- file.path(output_dir, "methylation_group_expression_boxplot.pdf")
stats_output_csv <- file.path(output_dir, "methylation_expression_stats.csv")

# Threshold for classifying methylation groups (e.g., based on beta values)
low_meth_threshold <- 0.3
high_meth_threshold <- 0.7

# --- Data Preparation --- 

# Check if input data exists (replace with your actual data objects)
if (!exists("methylation_data") || !exists("expression_data") || !exists("gene_annotations")) {
    stop("Input data objects (methylation_data, expression_data, gene_annotations) not found.")
}

# Merge methylation and expression data by gene_id
cat("Merging methylation and expression data...\n")
merged_data <- inner_join(methylation_data, expression_data, by = "gene_id")

# Merge with annotations (optional, if needed for stratification)
merged_data <- left_join(merged_data, gene_annotations, by = "gene_id")

# Remove NA values that would interfere with analysis
merged_data <- merged_data %>% 
    filter(!is.na(methylation_level) & !is.na(expression_value))

# Classify genes into methylation groups
cat("Classifying genes into methylation groups...\n")
merged_data <- merged_data %>%
    mutate(methylation_group = case_when(
        methylation_level <= low_meth_threshold ~ "Low",
        methylation_level >= high_meth_threshold ~ "High",
        TRUE ~ "Intermediate" # Genes between thresholds
    )) %>%
    # Convert to factor with desired levels
    mutate(methylation_group = factor(methylation_group, levels = c("Low", "Intermediate", "High")))

cat("Data prepared. Genes processed:", nrow(merged_data), "\n")
print(table(merged_data$methylation_group))

# --- Analysis --- 

# 1. Correlation Analysis (Continuous Methylation vs. Expression)
cat("\nPerforming correlation analysis...\n")
cor_test_result <- cor.test(~ methylation_level + expression_value,
                            data = merged_data,
                            method = "pearson") # Or "spearman"

print(cor_test_result)

# Generate scatter plot with correlation
p_correlation <- ggplot(merged_data, aes(x = methylation_level, y = expression_value)) +
    geom_point(alpha = 0.2, size = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color="blue") +
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
    labs(
        title = "Methylation vs. Gene Expression",
        x = "DNA Methylation Level (e.g., Beta Value)",
        y = "Gene Expression (e.g., log2(TPM+1))"
    ) +
    theme_minimal()

print(p_correlation)
ggsave(correlation_plot_file, p_correlation, width = 6, height = 5)
cat("Saved correlation plot to:", correlation_plot_file, "\n")

# 2. Group Comparison (Expression across Methylation Groups)
cat("\nComparing expression across methylation groups...\n")

# Filter for groups to compare (e.g., Low vs. High)
comparison_data <- merged_data %>%
    filter(methylation_group %in% c("Low", "High"))

# Perform t-test (or Wilcoxon test)
comp_test_result <- t.test(expression_value ~ methylation_group, data = comparison_data)
# Or: comp_test_result <- wilcox.test(expression_value ~ methylation_group, data = comparison_data)

print(comp_test_result)

# Generate boxplot with comparison
p_boxplot <- ggplot(comparison_data, aes(x = methylation_group, y = expression_value, fill = methylation_group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width=0.1, alpha=0.1, size=0.5) +
    stat_compare_means(method = "t.test", comparisons = list(c("Low", "High"))) +
    labs(
        title = "Gene Expression by Methylation Group",
        x = "Methylation Group",
        y = "Gene Expression (e.g., log2(TPM+1))"
    ) +
    theme_minimal() +
    theme(legend.position = "none")

print(p_boxplot)
ggsave(boxplot_file, p_boxplot, width = 5, height = 5)
cat("Saved boxplot to:", boxplot_file, "\n")

# --- Optional: Stratified Analysis (e.g., by gene biotype) --- 
# if ("gene_biotype" %in% colnames(merged_data)) {
#     cat("\nPerforming stratified analysis by gene biotype...\n")
#     # Example: Correlation for protein coding genes only
#     coding_data <- filter(merged_data, gene_biotype == "protein_coding")
#     if(nrow(coding_data) > 2) {
#         cor_test_coding <- cor.test(~ methylation_level + expression_value, data = coding_data)
#         print("Correlation for Protein Coding Genes:")
#         print(cor_test_coding)
#         # Add similar plots and tests for other biotypes if needed
#     } else {
#         cat("Not enough protein coding genes for stratified analysis.\n")
#     }
# }

# --- Save Statistics --- 
cat("\nSaving summary statistics...\n")
stats_summary <- tidy(cor_test_result) %>%
    mutate(analysis = "Overall Correlation") %>%
    bind_rows(tidy(comp_test_result) %>% mutate(analysis = "Low vs High Group Comparison")) %>%
    select(analysis, estimate, statistic, p.value, parameter, conf.low, conf.high, method, alternative)

write.csv(stats_summary, stats_output_csv, row.names = FALSE)
cat("Saved statistics to:", stats_output_csv, "\n")

cat("--- 01_methylation_expression.R finished ---\n") 