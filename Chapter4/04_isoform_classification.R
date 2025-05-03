# --- thesis_scripts/Chapter4/03_isoform_classification.R ---
# Purpose: Classifies genes based on the number of isoforms they produce
#          (e.g., single-isoform vs. multi-isoform genes).
#          Performs statistical tests (e.g., t-tests, Wilcoxon tests)
#          to compare features between these classes.
#          Generates box plots visualizing feature distributions.
# Chapter Relevance: 4.5 Splicing Effects on Protein Properties (comparing feature distributions)

# Load required libraries
library(ggplot2)
library(dplyr)
library(ggpubr) # For stat_compare_means
library(tidyr) # For pivot_longer
library(broom) # For tidying statistical test results

# --- Function Definitions ---

# Function to perform statistical comparison between groups and create boxplot
# Inputs: data = data frame
#         feature_col = string name of the feature column to compare
#         group_col = string name of the grouping column (e.g., 'IsoformClass')
#         ylab = string label for the y-axis (feature name)
#         title = string plot title
#         test_method = string specifying the statistical test (e.g., "t.test", "wilcox.test")
# Output: ggplot object
create_comparison_boxplot <- function(data, feature_col, group_col, ylab, title, test_method = "wilcox.test") {
  # Check if columns exist
  if (!all(c(feature_col, group_col) %in% colnames(data))) {
    warning(paste("Columns", feature_col, "or", group_col, "not found. Skipping plot:", title))
    return(NULL)
  }
  
  # Ensure group column is a factor
  data[[group_col]] <- as.factor(data[[group_col]])
  
  # Check if there are at least two groups
  if (length(levels(data[[group_col]])) < 2) {
    warning(paste("Need at least two groups in column", group_col, "for comparison. Skipping plot:", title))
    return(NULL)
  }
  
  # Check if feature column is numeric
  if (!is.numeric(data[[feature_col]])) {
      warning(paste("Feature column", feature_col, "is not numeric. Skipping plot:", title))
      return(NULL)
  }

  # Create boxplot
  p <- ggplot(data, aes(x = !!sym(group_col), y = !!sym(feature_col), fill = !!sym(group_col))) +
    geom_boxplot(outlier.shape = NA) + # Hide outliers for clarity, maybe add jitter points later
    geom_jitter(width = 0.2, alpha = 0.3, size=0.5) + # Add jittered points
    stat_compare_means(method = test_method, label.y.npc = 0.9) + # Add comparison p-value
    labs(title = title, x = "Isoform Class", y = ylab) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
          legend.position = "none") # Remove legend if fill is redundant with x-axis
          
  return(p)
}

# Function to classify genes based on isoform count
# Input: data = data frame with a column for isoform count (e.g., 'isoforms_per_gene')
#        isoform_count_col = string name of the isoform count column
#        threshold = integer, the threshold to distinguish single vs multi (e.g., 1 means >1 is multi)
# Output: data frame with a new 'IsoformClass' column
classify_genes_by_isoforms <- function(data, isoform_count_col = "isoforms_per_gene", threshold = 1) {
  if (!isoform_count_col %in% colnames(data)) {
    stop(paste("Isoform count column", isoform_count_col, "not found."))
  }
  
  data <- data %>%
    mutate(IsoformClass = ifelse(!!sym(isoform_count_col) > threshold,
                                 paste0("Multi-isoform (>", threshold, ")"),
                                 paste0("Single-isoform (<", threshold+1, ")"))) %>%
    mutate(IsoformClass = factor(IsoformClass, levels = c(paste0("Single-isoform (<", threshold+1, ")"), 
                                                          paste0("Multi-isoform (>", threshold, ")")))) # Set factor levels
  return(data)
}

# --- Execution ---

# Define Input/Output paths
# Assumes input data object `gene_features` (dataframe) is loaded.
# This object should contain gene-level features including 'isoforms_per_gene'.
output_dir <- "."
classification_plot_pdf <- file.path(output_dir, "gene_feature_by_isoform_class_boxplots.pdf")
classification_stats_csv <- file.path(output_dir, "gene_feature_by_isoform_class_stats.csv")

# Check if input data exists (replace with your actual data object)
if (!exists("gene_features")) {
  stop("Input data object 'gene_features' not found. Please load the required data.")
}

# 1. Classify Genes
cat("Classifying genes based on isoform count...\n")
gene_features_classified <- classify_genes_by_isoforms(gene_features, isoform_count_col = "isoforms_per_gene", threshold = 1)

# Print summary of classification
print(table(gene_features_classified$IsoformClass))

# 2. Define features to compare between classes
features_to_compare <- c(
  "gene_length_mean", 
  "exon_number_mean", 
  "prot_length_mean", 
  "hydropathy_mean", 
  "pi_mean", 
  "disorder_mean"
  # Add other relevant numeric features from gene_features
)

# Provide user-friendly names for plotting
feature_labels <- c(
  "Mean Gene Length (bp)",
  "Mean Exon Count",
  "Mean Protein Length (aa)",
  "Mean Hydropathy (GRAVY)",
  "Mean Isoelectric Point (pI)",
  "Mean Disorder Content (%)"
)
names(feature_labels) <- features_to_compare

# 3. Generate Comparison Boxplots
plot_list <- list()
stats_list <- list()

for (feature in features_to_compare) {
  if (feature %in% colnames(gene_features_classified) && is.numeric(gene_features_classified[[feature]])) {
    cat("Generating comparison plot for:", feature, "\n")
    
    # Create plot
    plot_obj <- create_comparison_boxplot(
      data = gene_features_classified,
      feature_col = feature,
      group_col = "IsoformClass",
      ylab = feature_labels[feature],
      title = paste(feature_labels[feature], "by Isoform Class"),
      test_method = "wilcox.test" # Use Wilcoxon test as default
    )
    
    if (!is.null(plot_obj)) {
      plot_list[[feature]] <- plot_obj
      
      # Perform the test and store results
      test_result <- tryCatch({
         wilcox.test(as.formula(paste(feature, "~ IsoformClass")), data = gene_features_classified)
      }, error = function(e) {warning(paste("Wilcoxon test failed for", feature, ":", e$message)); NULL})
      
      if (!is.null(test_result)) {
         stats_list[[feature]] <- tidy(test_result) %>% 
                                   mutate(feature = feature, 
                                          comparison = paste(levels(gene_features_classified$IsoformClass), collapse = " vs "))
      }
      
    } else {
       cat("Skipping plot generation for", feature, "due to earlier warnings.\n")
    }
  } else {
    warning(paste("Feature '", feature, "' not found or not numeric in the classified data. Skipping."))
  }
}

# 4. Save Plots
if (length(plot_list) > 0) {
  n_plots <- length(plot_list)
  n_cols <- min(3, n_plots)
  n_rows <- ceiling(n_plots / n_cols)

  pdf(classification_plot_pdf, width = n_cols * 4, height = n_rows * 4)
  gridExtra::grid.arrange(grobs = plot_list, ncol = n_cols)
  dev.off()
  
  cat("Saved comparison boxplots to:", classification_plot_pdf, "\n")
} else {
  cat("No comparison boxplots were generated.", "\n")
}

# 5. Save Statistics
if(length(stats_list) > 0) {
  combined_stats <- bind_rows(stats_list) %>%
                    select(feature, comparison, statistic, p.value, method, alternative)
  write.csv(combined_stats, classification_stats_csv, row.names = FALSE)
  cat("Saved comparison statistics to:", classification_stats_csv, "\n")
} else {
   cat("No statistical tests were performed or results saved.", "\n")
}

cat("--- 03_isoform_classification.R finished ---\n") 