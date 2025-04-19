# --- thesis_scripts/Chapter5/05_shape_analysis_post_analysis.R ---
# Purpose: Performs post-analysis on DNA shape predictions, potentially integrating
#          shape features with motif occurrences, chromatin states, or other genomic data.
#          May involve summarizing shape parameters around motifs or within specific regions.
# Chapter Relevance: 5.5 Motif Analysis and DNA Shape

# Load required libraries
library(GenomicRanges)
library(dplyr)
library(ggplot2)
# Potentially needs libraries for reading shape predictions (e.g., DNAshapeR results)
# or specific plotting packages depending on the analysis.

# --- Configuration & Input --- 

# Define input file paths (Replace with actual paths or loaded objects)
# Assumes the following might be loaded:
# - shape_predictions: Data structure containing DNA shape predictions 
#                      (e.g., output from DNAshapeR getShape() function, perhaps 
#                       a list of matrices or a combined data frame/GRanges).
#                       Needs to be associated with genomic coordinates.
# - regions_of_interest: GRanges object defining regions for shape analysis 
#                        (e.g., motif sites from monaLisa/motifmatchr, HMM states,
#                         promoters). Needs a column identifying the feature/motif.
# - annotations: Optional data frame or GRanges with additional info about the regions.

output_dir <- "results/chapter5/shape_analysis/"
shape_summary_file <- file.path(output_dir, "dna_shape_summary.csv")
shape_plot_prefix <- file.path(output_dir, "dna_shape_plot") # Prefix for plots

# Parameters
window_size <- 10 # Example: +/- window around a central point (e.g., motif center)
shape_parameters <- c("MGW", "ProT", "Roll", "HelT") # Parameters predicted by DNAshapeR

# --- Data Preparation --- 

# Check if input data exists
if (!exists("shape_predictions")) {
    stop("Input data 'shape_predictions' not found.")
}
if (!exists("regions_of_interest") || !is(regions_of_interest, "GRanges")) {
    stop("Input GRanges object 'regions_of_interest' not found or is not a GRanges object.")
}
# Add checks for specific structure/content of shape_predictions if known.
# Example: Check if it's a list and contains expected shape parameters.
# if (!is.list(shape_predictions) || !all(shape_parameters %in% names(shape_predictions))) {
#     stop("'shape_predictions' should be a list containing matrices for: ", 
#          paste(shape_parameters, collapse=", "))
# }

# Ensure output directory exists
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# --- Shape Analysis Functions --- 

# Example Function: Summarize shape parameters around region centers
# This is a placeholder and needs significant adaptation based on the actual
# structure of `shape_predictions` and `regions_of_interest`.

summarize_shape_around_centers <- function(shape_data, regions_gr, shape_params, window) {
    # This function needs to:
    # 1. Identify the center of each region in regions_gr.
    # 2. Extract shape parameters from shape_data for a window around each center.
    # 3. Align and average shape profiles for regions grouped by type (e.g., motif name).
    # 4. Return a summary data frame.
    
    # Placeholder - requires actual implementation based on data format
    warning("Function 'summarize_shape_around_centers' is a placeholder and needs implementation.")
    
    # Example structure of output (needs calculation)
    summary_df <- data.frame(
        feature_id = character(), # e.g., motif name or region type
        position = integer(),    # Position relative to center (-window to +window)
        shape_param = character(),# MGW, ProT, etc.
        mean_value = numeric(),
        sd_value = numeric(),
        n_regions = integer()
    )
    
    return(summary_df)
}

# Example Function: Plot average shape profiles
plot_average_shape <- function(shape_summary_df, shape_param_to_plot) {
    
    plot_data <- filter(shape_summary_df, shape_param == shape_param_to_plot)
    
    if (nrow(plot_data) == 0) {
        warning("No data found for shape parameter: ", shape_param_to_plot)
        return(NULL)
    }
    
    p <- ggplot(plot_data, aes(x = position, y = mean_value, color = feature_id)) +
        geom_line() +
        # Optional: Add ribbon for standard deviation/error
        # geom_ribbon(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value, fill = feature_id), alpha = 0.2) +
        facet_wrap(~ feature_id) + # Plot each feature/motif separately
        labs(
            title = paste("Average DNA Shape Profile:", shape_param_to_plot),
            x = "Position relative to region center (bp)",
            y = paste("Mean", shape_param_to_plot)
        ) +
        theme_minimal() +
        theme(legend.position = "none")
        
    return(p)
}

# --- Execution --- 

cat("Starting DNA shape post-analysis...\n")

# 1. Summarize Shape Features
# NOTE: This requires the actual implementation of summarize_shape_around_centers
#       based on how DNAshapeR results (`shape_predictions`) are stored and how
#       regions (`regions_of_interest`, e.g., motif locations) are defined.
#       The code below assumes the placeholder function is implemented.

shape_summary <- summarize_shape_around_centers(
    shape_data = shape_predictions,
    regions_gr = regions_of_interest,
    shape_params = shape_parameters,
    window = window_size
)

if (nrow(shape_summary) > 0) {
    # Save summary table
    write.csv(shape_summary, shape_summary_file, row.names = FALSE)
    cat("Saved DNA shape summary statistics to:", shape_summary_file, "\n")

    # 2. Generate Plots
    cat("Generating average shape profile plots...\n")
    plot_list_shape <- list()
    for (param in shape_parameters) {
        p_shape <- plot_average_shape(shape_summary, param)
        if (!is.null(p_shape)) {
            plot_list_shape[[param]] <- p_shape
            ggsave(paste0(shape_plot_prefix, "_", param, ".pdf"), p_shape, width=7, height=5)
        }
    }
    cat("Saved DNA shape plots with prefix:", shape_plot_prefix, "\n")
    
    # Display first plot (optional)
    # if(length(plot_list_shape) > 0) { print(plot_list_shape[[1]]) }
    
} else {
    cat("Shape summary data frame is empty. Cannot save or plot results. Check 'summarize_shape_around_centers' implementation.\n")
}

cat("--- 05_shape_analysis_post_analysis.R finished ---\n") 