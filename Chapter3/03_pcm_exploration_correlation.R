# --- thesis_scripts/Chapter3/03_pcm_exploration_correlation.R ---
# Purpose: Load calculated PCM features and perform exploratory analysis.
#          Calculates and visualizes feature correlations (corrplot).
#          Compares feature distributions between nORFs, ORFs, and random sequences.
# Chapter Relevance: 3.3 Proteochemometric Analysis, 3.4 High-Dimensional Statistical Analysis (partially)

# Load necessary packages
library(data.table)
library(ggplot2)
library(corrplot)
library(dplyr)
library(tidyr) # For gather, needed by ggbetweenstats
library(ggpubr) # For ggboxplot, stat_compare_means
library(gridExtra) # For grid.arrange
library(ggplot2)
library(ggExtra)
library(scales) # For scaleFUN
library(RColorBrewer) # For corrplot colors

# Define Input/Output files (Ensure these match output from 02_...R)
norf_features_file <- "nORFs_pcm_features.csv"
orf_features_file <- "ORFs_pcm_features.csv"
random_features_file <- "Random_pcm_features.csv"

# --- Load PCM Features --- 
cat("Loading PCM feature files...\n")
if (!file.exists(norf_features_file)) stop("nORF features file not found: ", norf_features_file)
if (!file.exists(orf_features_file)) stop("ORF features file not found: ", orf_features_file)
if (!file.exists(random_features_file)) stop("Random features file not found: ", random_features_file)

nORF <- fread(norf_features_file)
ORF <- fread(orf_features_file)
Random <- fread(random_features_file)
cat("Features loaded.\n")

# Clean up potential index columns from csv writing/reading
nORF[, V1 := NULL] # Remove default data.table index column if present
ORF[, V11 := NULL]
Random[, V1 := NULL]
# Original script had X.1 and X, remove if they exist from potential earlier saves
if("X" %in% names(nORF)) nORF[, X := NULL]
if("X.1" %in% names(nORF)) nORF[, X.1 := NULL]
if("X" %in% names(ORF)) ORF[, X := NULL]
if("X.1" %in% names(ORF)) ORF[, X.1 := NULL]
if("X" %in% names(Random)) Random[, X := NULL]
if("X.1" %in% names(Random)) Random[, X.1 := NULL]

# Assign consistent ID for random sequences
Random[, id := "random"] # Original script did this, ensures consistency

# --- Combine Data & Prepare for Correlation --- 
cat("Combining nORF and ORF features for correlation analysis...\n")
corrFrame <- rbindlist(list(nORF, ORF), fill = TRUE) # Use fill=TRUE in case of minor column mismatches

# Optional: Check for and handle duplicated column names before correlation
dup_cols <- duplicated(names(corrFrame))
if (any(dup_cols)) {
    warning("Duplicate column names found: ", paste(names(corrFrame)[dup_cols], collapse=", "), ". Removing duplicates.")
    corrFrame <- corrFrame[, !dup_cols, with = FALSE]
}


# Prepare numeric frame for correlation
# Exclude non-numeric columns (id, type) and handle potential NAs/Infs
non_numeric_cols <- c("id", "type")
feature_cols <- setdiff(names(corrFrame), non_numeric_cols)
numeric_corrFrame <- corrFrame[, ..feature_cols] # Select feature columns using data.table syntax

# Convert all columns to numeric, coercing errors
for (j in seq_along(numeric_corrFrame)) {
    set(numeric_corrFrame, i = NULL, j = j, value = suppressWarnings(as.numeric(numeric_corrFrame[[j]])))
}

# Remove columns that are all NA after coercion
all_na_cols <- names(which(colSums(is.na(numeric_corrFrame)) == nrow(numeric_corrFrame)))
if (length(all_na_cols) > 0) {
    cat("Removing columns with all NAs:", paste(all_na_cols, collapse=", "), "\n")
    numeric_corrFrame <- numeric_corrFrame[, .SD, .SDcols = !all_na_cols]
}

# Impute remaining NAs (e.g., with column mean) or remove rows with NAs
# Option 1: Remove rows with any NAs (simpler)
# numeric_corrFrame <- na.omit(numeric_corrFrame)
# Option 2: Impute NAs (e.g., with mean)
for (j in names(numeric_corrFrame)) {
    set(numeric_corrFrame, i = which(is.na(numeric_corrFrame[[j]])), j = j, value = mean(numeric_corrFrame[[j]], na.rm = TRUE))
}

# Remove zero-variance columns
variances <- apply(numeric_corrFrame, 2, var, na.rm = TRUE)
zero_var_cols <- names(which(variances < 1e-10))
if (length(zero_var_cols) > 0) {
    cat("Removing zero/low variance columns:", paste(zero_var_cols, collapse=", "), "\n")
    numeric_corrFrame <- numeric_corrFrame[, .SD, .SDcols = !names(numeric_corrFrame) %in% zero_var_cols]
}

# --- Correlation Plot --- 
cat("Calculating and plotting correlation matrix...\n")
if (ncol(numeric_corrFrame) > 1) {
    corrTable <- cor(numeric_corrFrame, method = "pearson", use = "pairwise.complete.obs")

    # Plotting the correlation matrix
    png("pcm_feature_correlation_plot.png", width=1200, height=1200)
    corrplot(corrTable,
             method = "color",
             order = "hclust", # Hierarchical clustering order
             type = "upper",      # Show upper triangle
             tl.cex = 0.5,      # Label size
             tl.col = "black",
             tl.srt = 45,       # Label rotation
             addCoef.col = "grey30", # Add correlation coefficient
             number.cex = 0.4,     # Size of coefficients
             cl.cex = 0.6,       # Legend scale size
             #addrect = 4,       # Add rectangles around clusters (adjust number)
             diag = FALSE,        # Don't plot diagonal
             col=brewer.pal(n=8, name="RdYlBu") # Color palette
             )
    dev.off()
    cat("Saved correlation plot to pcm_feature_correlation_plot.png\n")
} else {
    cat("Skipping correlation plot - insufficient valid numeric columns.\n")
}

# --- Feature Distribution Comparison Plots --- 
cat("Preparing data for feature distribution plots...\n")
# Combine all three groups
jointFrame <- rbindlist(list(nORF, ORF, Random), fill = TRUE)

# Rename 'type' column to 'group' for clarity in plots
setnames(jointFrame, "type", "group")

# Define the plotting function from original script
plotTripletPlot <- function(dataset, dimension_name, output_dir = ".") {
    if (!dimension_name %in% names(dataset)) {
        warning("Dimension ", dimension_name, " not found in dataset.")
        return(NULL)
    }

    # Ensure dimension is numeric
    dataset[, (dimension_name) := suppressWarnings(as.numeric(get(dimension_name)))]
    plot_data <- dataset[!is.na(get(dimension_name)), .SD, .SDcols = c("group", dimension_name)]

    if (nrow(plot_data) < 10 || length(unique(plot_data$group)) < 2) {
        warning("Insufficient data or groups for plotting dimension: ", dimension_name)
        return(NULL)
    }

    cat(" Plotting distribution for feature:", dimension_name, "\n")

    # Calculate plot limits and comparison positions dynamically
    plot_min <- min(plot_data[[dimension_name]], na.rm = TRUE)
    plot_max <- max(plot_data[[dimension_name]], na.rm = TRUE)
    plot_range <- plot_max - plot_min

    # Position significance labels above the max data point
    y_pos_scale <- 0.05 # Adjust this multiplier to control spacing
    sig_pos_1 <- plot_max + y_pos_scale * plot_range
    sig_pos_2 <- plot_max + 2 * y_pos_scale * plot_range
    sig_pos_3 <- plot_max + 3 * y_pos_scale * plot_range

    # Define comparisons (ensure groups exist)
    existing_groups <- unique(plot_data$group)
    my_comparisons <- list()
    if (all(c("nORF", "ORF") %in% existing_groups)) my_comparisons <- c(my_comparisons, list(c("nORF", "ORF")))
    if (all(c("nORF", "Random") %in% existing_groups)) my_comparisons <- c(my_comparisons, list(c("nORF", "Random")))
    if (all(c("ORF", "Random") %in% existing_groups)) my_comparisons <- c(my_comparisons, list(c("ORF", "Random")))

    # Assign y positions for labels (adjust logic if fewer than 3 comparisons)
    label_y_positions <- c(sig_pos_1, sig_pos_2, sig_pos_3)[1:length(my_comparisons)]

    # Rounding function for scales
    scaleFUN <- function(x) sprintf("%.2f", x)

    # Adjust y-axis limits for boxplot to accommodate labels
    boxplot_ymax <- max(plot_max, label_y_positions, na.rm = TRUE) * 1.1
    boxplot_ymin <- plot_min - 0.1 * plot_range

    # Boxplot with significance tests
    plot_box <- ggboxplot(plot_data, x = "group", y = dimension_name, xlab = " ", size = 0.1,
                          color = "group") +
        scale_color_manual(values = c("nORF" = "#00AFBB", "ORF" = "#E7B800", "Random" = "#D22B2B")) +
        scale_fill_manual(values = c("nORF" = "#00AFBB", "ORF" = "#E7B800", "Random" = "#D22B2B")) +
        theme(axis.text.y = element_text(size = 10)) +
        scale_y_continuous(limits = c(boxplot_ymin, boxplot_ymax)) +
        labs(y = dimension_name)

    # Add comparisons if any exist
    if (length(my_comparisons) > 0) {
        plot_box <- plot_box + stat_compare_means(comparisons = my_comparisons, label = "p.signif",
                                                    label.y = label_y_positions, method = "t.test")
    }

    # Density plot (absolute)
    plot_density_abs <- ggplot(plot_data, aes_string(x = dimension_name)) +
        geom_density(aes(color = group), position = "identity", alpha = 0.8, size = 1) +
        scale_y_continuous(labels = scaleFUN) +
        scale_color_manual(values = c("nORF" = "#00AFBB", "ORF" = "#E7B800", "Random" = "#D22B2B")) +
        scale_fill_manual(values = c("nORF" = "#00AFBB", "ORF" = "#E7B800", "Random" = "#D22B2B")) +
        theme(legend.position = "none", text = element_text(size = 15), axis.text.x = element_blank(),
              axis.title.x = element_blank(), axis.title.y = element_text(size=12))

    # Density plot (relative/fill)
    plot_density_rel <- ggplot(plot_data, aes_string(x = dimension_name)) +
        geom_density(aes(color = group, fill = group), position = "fill", alpha = 0.5) +
        scale_color_manual(values = c("nORF" = "#00AFBB", "ORF" = "#E7B800", "Random" = "#D22B2B")) +
        scale_fill_manual(values = c("nORF" = "#00AFBB", "ORF" = "#E7B800", "Random" = "#D22B2B")) +
        theme(legend.position = "none", text = element_text(size = 15), axis.title.y = element_text(size=12)) +
        labs(y = "Relative Density", x = dimension_name)

    # Arrange plots
    figure <- grid.arrange(arrangeGrob(plot_density_abs, plot_density_rel), plot_box, ncol = 2, widths = c(2, 1))

    # Save the plot
    output_filename <- file.path(output_dir, paste0("feature_distribution_", dimension_name, ".png"))
    ggsave(filename = output_filename, plot = figure, width = 10, height = 5, dpi = 150)

    return(figure)
}

# --- Select Top Features and Plot Distributions ---
# Use features identified in the original script's XGBoost importance or select manually
# Example: Using features mentioned in the plotTripletPlot calls
top_features_to_plot <- c("prop3.G1.residue100", "VS577", "DAYM780201.lag30",
                          "BHAR880101.lag29", "CIDH920105.lag10", "BIGC670101.lag30",
                          "V1093", "V532", "V240", "V534", "V1122", "V1119",
                          "V1007", "V1011", "V944", "V1027", "V702") # Replace Vxxx with actual feature names if possible

# Check which features actually exist in the combined frame
feature_cols_available <- names(jointFrame)
top_features_to_plot <- intersect(top_features_to_plot, feature_cols_available)

cat("\nPlotting distributions for selected top features:\n")
if (length(top_features_to_plot) > 0) {
    plots_list <- lapply(top_features_to_plot, function(feature) {
        plotTripletPlot(jointFrame, feature)
    })
} else {
    cat("Warning: None of the selected top features were found in the data.\n")
}


# --- Optional: 2D Density Plot (Example from original script) ---
# Needs actual feature names 'length' and 'v240' (or equivalent)
# if ("length" %in% names(jointFrame) && "V240" %in% names(jointFrame)) {
#     cat("\nCreating 2D density plot...\n")
#     p_2d <- ggplot(jointFrame, aes(x=length, y=V240) ) +
#         geom_density_2d() +
#         labs(title="2D Density: Sequence Length vs Feature V240")
#     print(p_2d)
#     ggsave("feature_vs_length_2d_density.png", plot=p_2d)
# } else {
#     cat("Skipping 2D density plot - required columns not found.\n")
# }

cat("--- 03_pcm_exploration_correlation.R finished ---\n") 