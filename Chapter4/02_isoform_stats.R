# --- thesis_scripts/Chapter4/01_isoform_stats.R ---
# Purpose: Calculates and compares correlations between different alternative splicing
#          event types (e.g., ES, IR, A3, A5) at both the isoform and gene levels.
#          Helps understand co-occurrence or mutual exclusivity of splicing events.
#          Generates correlation heatmaps and summary statistics tables.
# Chapter Relevance: 4.4 Isoform Expression and Distribution / 4.6 Classification and Combinatorial Analysis

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2) # For melt function
library(patchwork) # For combining plots
library(xtable) # For creating LaTeX tables
library(kableExtra) # For formatting tables
library(stats) # For cor function

# --- Function Definitions ---

# Function to calculate p-value for correlation
# Inputs: r = correlation coefficient, n = sample size
# Output: p-value
cor_to_p <- function(r, n) {
  if (is.na(r) || abs(r) == 1 || n <= 2) {
    return(NA_real_) # Return NA if correlation is undefined or sample size is too small
  }
  t_stat <- (r * sqrt(n - 2)) / sqrt(1 - r^2)
  p_val <- 2 * pt(abs(t_stat), df = n - 2, lower.tail = FALSE)
  return(p_val)
}

# Function to process splicing data for a brain region
# Inputs: splicing_data = data frame with isoform_id and binary columns for splice events
#         gene_data = data frame mapping isoform_id to gene_id
# Output: list containing isoform correlation matrix, gene correlation matrix,
#         available splice types, number of genes, and number of isoforms
process_brain_region <- function(splicing_data, gene_data) {
  # Define standard splicing types
  splicing_types <- c("ES", "MEE", "MES", "IR", "A5", "A3", "ATSS", "ATTS")

  # Identify splicing types actually present in the data
  available_types <- intersect(splicing_types, colnames(splicing_data))
  if(length(available_types) == 0) {
    stop("None of the expected splicing type columns are present in the input splicing_data.")
  }
  print(paste("Available splicing types:", paste(available_types, collapse=", ")))

  # --- Isoform-wise analysis ---
  # Ensure data is numeric
  isoform_matrix <- sapply(splicing_data[, available_types], as.numeric)
  rownames(isoform_matrix) <- splicing_data$isoform_id # Keep track of isoforms

  # Handle potential NAs introduced by as.numeric coercion (though unlikely for 0/1 data)
  isoform_matrix[is.na(isoform_matrix)] <- 0

  # Calculate isoform correlation matrix
  print("Calculating isoform correlation matrix...")
  isoform_cor <- cor(isoform_matrix, use = "pairwise.complete.obs")
  n_isoforms <- nrow(isoform_matrix)
  print(paste("Isoform correlation matrix calculated for", n_isoforms, "isoforms."))

  # --- Gene-wise analysis ---
  print("Performing gene-wise analysis...")
  # Ensure gene_data has the necessary columns
  if (!all(c("isoform_id", "gene_id") %in% colnames(gene_data))) {
      stop("gene_data must contain 'isoform_id' and 'gene_id' columns.")
  }
  # Ensure splicing_data has 'isoform_id'
  if (!"isoform_id" %in% colnames(splicing_data)) {
      stop("splicing_data must contain 'isoform_id' column.")
  }

  # Merge splicing data with gene mapping
  gene_splicing <- merge(splicing_data[, c("isoform_id", available_types)],
                         gene_data[, c("isoform_id", "gene_id")],
                         by = "isoform_id")

  # Aggregate splicing events at the gene level (presence/absence: 1 if any isoform has it)
  gene_matrix_df <- gene_splicing %>%
    group_by(gene_id) %>%
    summarise(across(all_of(available_types), ~as.numeric(any(. > 0, na.rm = TRUE))), .groups = "drop")

  # Convert to matrix for correlation
  gene_matrix <- as.matrix(gene_matrix_df[, -1]) # Exclude gene_id column
  rownames(gene_matrix) <- gene_matrix_df$gene_id

  # Calculate gene correlation matrix
  print("Calculating gene correlation matrix...")
  gene_cor <- cor(gene_matrix, use = "pairwise.complete.obs")
  n_genes <- nrow(gene_matrix)
  print(paste("Gene correlation matrix calculated for", n_genes, "genes."))

  return(list(isoform_cor = isoform_cor, gene_cor = gene_cor,
              available_types = available_types,
              n_genes = n_genes, n_isoforms = n_isoforms))
}


# Function to calculate statistics (Mean, SD, CV, p-value, FDR) for correlations
# Inputs: forebrain, hindbrain, midbrain = results from process_brain_region
#         cor_type = "gene" or "isoform"
# Output: data frame with summary statistics for each pair of splice types
calculate_stats <- function(forebrain, hindbrain, midbrain, cor_type) {
  # Select the appropriate correlation matrix and sample size
  if(cor_type == "gene") {
    fb_cor <- forebrain$gene_cor
    hb_cor <- hindbrain$gene_cor
    mb_cor <- midbrain$gene_cor
    # Use the minimum sample size across tissues for conservative p-value calculation
    n <- min(forebrain$n_genes, hindbrain$n_genes, midbrain$n_genes)
  } else {
    fb_cor <- forebrain$isoform_cor
    hb_cor <- hindbrain$isoform_cor
    mb_cor <- midbrain$isoform_cor
    n <- min(forebrain$n_isoforms, hindbrain$n_isoforms, midbrain$n_isoforms)
  }
  available_types <- forebrain$available_types # Assume types are the same across regions

  # Get all unique pairs of splice types
  all_pairs <- expand.grid(row = available_types, col = available_types, stringsAsFactors = FALSE)
  # Keep only pairs where row < col alphabetically to avoid duplicates and self-correlations
  all_pairs <- all_pairs[as.character(all_pairs$row) < as.character(all_pairs$col), ]

  # Calculate statistics for each pair
  stats <- all_pairs %>%
    rowwise() %>%
    mutate(
      Pair = paste(row, col, sep = "-"),
      Forebrain = fb_cor[row, col],
      Hindbrain = hb_cor[row, col],
      Midbrain = mb_cor[row, col],
      Mean = mean(c(Forebrain, Hindbrain, Midbrain), na.rm = TRUE),
      SD = sd(c(Forebrain, Hindbrain, Midbrain), na.rm = TRUE),
      CV = ifelse(abs(Mean) > 1e-6, (SD / abs(Mean)) * 100, NA), # Avoid division by zero
      # Calculate p-value based on the mean correlation and minimum sample size
      p_value = cor_to_p(Mean, n)
    ) %>%
    ungroup() # End rowwise processing

  # Select and order columns
  stats <- stats[, c("Pair", "Mean", "SD", "CV", "p_value")]
  stats <- stats[order(stats$p_value, -abs(stats$Mean)), ] # Order by p-value, then by magnitude of correlation

  # Calculate FDR adjusted p-value
  stats$FDR <- p.adjust(stats$p_value, method = "BH")

  # Assign significance levels based on FDR
  stats$Significance <- cut(stats$FDR,
                            breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), # Use -Inf for lower bound
                            labels = c("***", "**", "*", "ns"),
                            include.lowest = TRUE, right = FALSE) # Intervals are [low, high)

  return(stats)
}


# Function to create heatmap with significance levels
# Inputs: stats_data = output from calculate_stats
# Output: ggplot object
create_heatmap_with_significance <- function(stats_data, title) {
  # Reshape the data for the heatmap
  # Separate Pair column into Type1 and Type2
  heatmap_data <- stats_data %>%
    separate(Pair, into = c("Type1", "Type2"), sep = "-", remove = FALSE) # Keep Pair column

  # Create a symmetric matrix for plotting
  all_types <- sort(unique(c(heatmap_data$Type1, heatmap_data$Type2)))
  heatmap_matrix <- matrix(NA, nrow = length(all_types), ncol = length(all_types),
                           dimnames = list(all_types, all_types))

  # Fill the matrix with Mean correlation values
  for (i in 1:nrow(heatmap_data)) {
    heatmap_matrix[heatmap_data$Type1[i], heatmap_data$Type2[i]] <- heatmap_data$Mean[i]
    heatmap_matrix[heatmap_data$Type2[i], heatmap_data$Type1[i]] <- heatmap_data$Mean[i]
  }
  diag(heatmap_matrix) <- 1 # Set diagonal to 1 (self-correlation)

  # Melt the matrix for ggplot
  melted_matrix <- melt(heatmap_matrix, na.rm = FALSE) # Keep NAs for diagonal
  colnames(melted_matrix) <- c("Type1", "Type2", "Correlation")

  # Add significance levels back to the melted data
  # Need to match pairs in both directions (e.g., A-B and B-A)
  stats_data_sym <- stats_data %>%
    select(Pair, Significance) %>%
    separate(Pair, into = c("T1", "T2"), sep="-") %>%
    bind_rows(., rename(., T1 = T2, T2 = T1)) # Add symmetric pairs

  melted_matrix <- melted_matrix %>%
    left_join(stats_data_sym, by = c("Type1" = "T1", "Type2" = "T2")) %>%
    # Set significance for diagonal to empty string
    mutate(Significance = ifelse(Type1 == Type2, "", as.character(Significance)))

  # Create the heatmap plot
  ggplot(melted_matrix, aes(x = Type1, y = Type2, fill = Correlation)) +
    geom_tile(color = "grey80") + # Add tile borders
    # Add text for correlation value and significance level
    geom_text(aes(label = ifelse(Type1 != Type2 & !is.na(Correlation),
                                 paste(sprintf("%.2f", Correlation), Significance, sep = "\n"),
                                 "1.00")), # Show 1.00 on diagonal
              size = 3, color = "black") +
    # Use a diverging color scale
    scale_fill_gradient2(low = "#313695", high = "#A50026", mid = "#FFFFBF", # Blue-White-Red
                         midpoint = 0, limit = c(-1, 1), space = "Lab",
                         name="Mean Correlation") +
    theme_minimal(base_size = 10) + # Adjust base font size
    # Improve axis labels and theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=9),
          axis.text.y = element_text(size=9),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5, face="bold"),
          plot.subtitle = element_text(hjust = 0.5, size=9),
          panel.grid.major = element_blank(), # Remove grid lines
          panel.border = element_blank(),
          legend.position = "right") +
    coord_fixed() + # Ensure square tiles
    labs(title = title,
         subtitle = "*** FDR < 0.001, ** FDR < 0.01, * FDR < 0.05, ns: not significant")
}

# --- Execution ---

# Define Input/Output paths
# Assumes input data objects like combinedForebrain, combinedHindBrain, combinedMidBrain
# are loaded into the R environment beforehand.
# These objects should contain $AlternativeSplicingAnalysis (dataframe) and $exons (dataframe or GRanges)
output_dir <- "."
isof_stats_csv <- file.path(output_dir, "isoform_wise_splicing_stats.csv")
gene_stats_csv <- file.path(output_dir, "gene_wise_splicing_stats.csv")
gene_heatmap_pdf <- file.path(output_dir, "gene_wise_correlation_heatmap_significance.pdf")
isoform_heatmap_pdf <- file.path(output_dir, "isoform_wise_correlation_heatmap_significance.pdf")

# Check if input data exists (replace with your actual data objects)
if (!exists("combinedForebrain") || !exists("combinedHindBrain") || !exists("combinedMidBrain")) {
    stop("Input data objects (combinedForebrain, combinedHindBrain, combinedMidBrain) not found.")
}

# Process each brain region
result_forebrain <- process_brain_region(combinedForebrain$AlternativeSplicingAnalysis, combinedForebrain$exons)
result_hindbrain <- process_brain_region(combinedHindBrain$AlternativeSplicingAnalysis, combinedHindBrain$exons)
result_midbrain <- process_brain_region(combinedMidBrain$AlternativeSplicingAnalysis, combinedMidBrain$exons)

# Calculate gene-wise and isoform-wise statistics
gene_stats <- calculate_stats(result_forebrain, result_hindbrain, result_midbrain, "gene")
isoform_stats <- calculate_stats(result_forebrain, result_hindbrain, result_midbrain, "isoform")

# Print the stats tables
print("Gene-wise Statistics:")
print(gene_stats)
print("Isoform-wise Statistics:")
print(isoform_stats)

# Save stats tables to CSV
write.csv(isoform_stats, isof_stats_csv, row.names = FALSE)
write.csv(gene_stats, gene_stats_csv, row.names = FALSE)
cat("Saved statistics tables to:", isof_stats_csv, "and", gene_stats_csv, "\n")

# Create the heatmaps
gene_heatmap <- create_heatmap_with_significance(gene_stats, "Gene-wise Splicing Type Correlations")
isoform_heatmap <- create_heatmap_with_significance(isoform_stats, "Isoform-wise Splicing Type Correlations")

# Save the heatmaps
ggsave(gene_heatmap_pdf, gene_heatmap, width = 8, height = 7)
ggsave(isoform_heatmap_pdf, isoform_heatmap, width = 8, height = 7)
cat("Saved heatmap plots to:", gene_heatmap_pdf, "and", isoform_heatmap_pdf, "\n")

# Display the plots (optional)
# print(gene_heatmap)
# print(isoform_heatmap)

# --- Create formatted tables for output --- 
# Function to format tables using kableExtra
format_kable_table <- function(stats_data, caption) {
  stats_data %>%
    mutate(
      Mean = cell_spec(sprintf("%.3f", Mean),
                       color = case_when(
                           Mean > 0.5 ~ "darkgreen",
                           Mean < -0.5 ~ "darkred",
                           Mean > 0.2 ~ "green",
                           Mean < -0.2 ~ "red",
                           TRUE ~ "black" # Default color
                           ),
                       bold = ifelse(Significance != "ns", TRUE, FALSE)),
      SD = sprintf("%.3f", SD),
      CV = sprintf("%.3f", CV),
      p_value = sprintf("%.2e", p_value),
      FDR = sprintf("%.2e", FDR),
      Significance = cell_spec(Significance,
                               color = case_when(
                                   Significance == "***" ~ "#A50026", # Dark Red
                                   Significance == "**" ~ "#F46D43", # Orange-Red
                                   Significance == "*" ~ "#FDAE61", # Light Orange
                                   TRUE ~ "#66BD63" # Green for ns (or choose grey)
                                   ),
                               bold = TRUE)
    ) %>%
    kbl(escape = FALSE, caption = caption, booktabs = TRUE) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                  full_width = FALSE, font_size = 10) %>%
    column_spec(1, bold = TRUE) # Make Pair column bold
}

# Create formatted tables
kable_gene_stats <- format_kable_table(gene_stats, "Gene Aggregated Splicing Statistics")
kable_isoform_stats <- format_kable_table(isoform_stats, "Isoform Splicing Statistics")

# Print or save the kable tables (saving requires webshot/webshot2 package)
print(kable_gene_stats)
# save_kable(kable_gene_stats, "gene_stats_table.pdf") # Requires webshot2
print(kable_isoform_stats)
# save_kable(kable_isoform_stats, "isoform_stats_table.pdf") # Requires webshot2

cat("--- 01_isoform_stats.R finished ---\n") 