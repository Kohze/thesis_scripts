# --- thesis_scripts/Chapter4/02_pearson_corr.R ---
# Purpose: Calculates Pearson correlations between various protein features
#          (lengths, hydropathy, PI, disorder, etc.) and isoform diversity metrics
#          (number of isoforms per gene, proportion of protein-coding isoforms).
#          Generates scatter plots with correlation coefficients and p-values.
# Chapter Relevance: 4.5 Splicing Effects on Protein Properties

# Load required libraries
library(ggplot2)
library(dplyr)
library(stats) # For cor.test
library(ggpubr) # For stat_cor
library(gridExtra) # For arranging plots

# --- Function Definitions ---

# Function to create a scatter plot with correlation test results
# Inputs: data = data frame containing the two variables
#         x_var = string name of the x-axis variable column
#         y_var = string name of the y-axis variable column
#         xlab = string label for the x-axis
#         ylab = string label for the y-axis
#         title = string plot title
# Output: ggplot object
create_correlation_plot <- function(data, x_var, y_var, xlab, ylab, title) {
  # Check if columns exist
  if (!all(c(x_var, y_var) %in% colnames(data))) {
    warning(paste("Columns", x_var, "or", y_var, "not found in the data. Skipping plot:", title))
    return(NULL)
  }
  
  # Remove rows with NA in either variable
  plot_data <- data %>% dplyr::select(!!sym(x_var), !!sym(y_var)) %>% na.omit()
  
  if (nrow(plot_data) < 3) {
      warning(paste("Insufficient non-NA data points (<3) for correlation. Skipping plot:", title))
      return(NULL)
  }
  
  # Perform correlation test
  cor_test_result <- tryCatch({
    cor.test(plot_data[[x_var]], plot_data[[y_var]], method = "pearson")
  }, error = function(e) {
    warning(paste("Correlation test failed for", title, ":", e$message))
    return(NULL)
  })
  
  if (is.null(cor_test_result)) {
    return(NULL)
  }
  
  # Build the plot
  p <- ggplot(plot_data, aes(x = !!sym(x_var), y = !!sym(y_var))) +
    geom_point(alpha = 0.5, size=1) + # Make points slightly transparent
    geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") + # Add linear model fit line
    stat_cor(method = "pearson", # Add correlation coefficient and p-value
             label.x.npc = "left", # Position label
             label.y.npc = "top",
             aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")))
    labs(x = xlab, y = ylab, title = title) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size=11, face="bold"))
    
  return(p)
}

# --- Execution ---

# Define Input/Output paths
# Assumes input data object `gene_features` (dataframe) is loaded.
# This object should contain gene-level features and isoform diversity metrics.
output_dir <- "."
corr_plot_pdf <- file.path(output_dir, "gene_feature_isoform_diversity_correlations.pdf")

# Check if input data exists (replace with your actual data object)
if (!exists("gene_features")) {
  stop("Input data object 'gene_features' not found. Please load the required data.")
}

# Define pairs of variables for correlation plots
# Structure: list(list(x_var, y_var, xlab, ylab, title), ...)
plot_definitions <- list(
  # --- Correlations with Number of Isoforms --- 
  list(x_var = "gene_length_mean", y_var = "isoforms_per_gene",
       xlab = "Mean Gene Length (bp)", ylab = "Isoforms per Gene",
       title = "Gene Length vs. Isoform Count"),
  list(x_var = "exon_number_mean", y_var = "isoforms_per_gene",
       xlab = "Mean Exon Count", ylab = "Isoforms per Gene",
       title = "Exon Count vs. Isoform Count"),
  list(x_var = "prot_length_mean", y_var = "isoforms_per_gene",
       xlab = "Mean Protein Length (aa)", ylab = "Isoforms per Gene",
       title = "Protein Length vs. Isoform Count"),
  list(x_var = "hydropathy_mean", y_var = "isoforms_per_gene",
       xlab = "Mean Hydropathy (GRAVY)", ylab = "Isoforms per Gene",
       title = "Hydropathy vs. Isoform Count"),
  list(x_var = "pi_mean", y_var = "isoforms_per_gene",
       xlab = "Mean Isoelectric Point (pI)", ylab = "Isoforms per Gene",
       title = "pI vs. Isoform Count"),
  list(x_var = "disorder_mean", y_var = "isoforms_per_gene",
       xlab = "Mean Disorder Content (%)", ylab = "Isoforms per Gene",
       title = "Disorder Content vs. Isoform Count"),

  # --- Correlations with Proportion of Protein-Coding Isoforms --- 
  # (Assuming a column like 'prop_coding_isoforms' exists in gene_features)
  list(x_var = "gene_length_mean", y_var = "prop_coding_isoforms",
       xlab = "Mean Gene Length (bp)", ylab = "% Coding Isoforms",
       title = "Gene Length vs. % Coding Isoforms"),
  list(x_var = "exon_number_mean", y_var = "prop_coding_isoforms",
       xlab = "Mean Exon Count", ylab = "% Coding Isoforms",
       title = "Exon Count vs. % Coding Isoforms"),
  list(x_var = "prot_length_mean", y_var = "prop_coding_isoforms",
       xlab = "Mean Protein Length (aa)", ylab = "% Coding Isoforms",
       title = "Protein Length vs. % Coding Isoforms"),
  list(x_var = "hydropathy_mean", y_var = "prop_coding_isoforms",
       xlab = "Mean Hydropathy (GRAVY)", ylab = "% Coding Isoforms",
       title = "Hydropathy vs. % Coding Isoforms"),
  list(x_var = "pi_mean", y_var = "prop_coding_isoforms",
       xlab = "Mean Isoelectric Point (pI)", ylab = "% Coding Isoforms",
       title = "pI vs. % Coding Isoforms"),
  list(x_var = "disorder_mean", y_var = "prop_coding_isoforms",
       xlab = "Mean Disorder Content (%)", ylab = "% Coding Isoforms",
       title = "Disorder Content vs. % Coding Isoforms")
)

# Generate all plots
plot_list <- list()
for (pdef in plot_definitions) {
  cat("Generating plot:", pdef$title, "\n")
  plot_obj <- create_correlation_plot(gene_features,
                                      pdef$x_var, pdef$y_var,
                                      pdef$xlab, pdef$ylab, pdef$title)
  if (!is.null(plot_obj)) {
    plot_list[[length(plot_list) + 1]] <- plot_obj
  }
}

# Arrange plots into a grid and save to PDF
if (length(plot_list) > 0) {
  # Determine grid layout (e.g., 3 columns)
  n_plots <- length(plot_list)
  n_cols <- min(3, n_plots) # Max 3 columns
  n_rows <- ceiling(n_plots / n_cols)

  pdf(corr_plot_pdf, width = n_cols * 4, height = n_rows * 3.5)
  grid.arrange(grobs = plot_list, ncol = n_cols)
  dev.off()
  
  cat("Saved correlation plots to:", corr_plot_pdf, "\n")
} else {
  cat("No plots were generated.", "\n")
}

cat("--- 02_pearson_corr.R finished ---\n") 