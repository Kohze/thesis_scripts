# --- thesis_scripts/Chapter5/03_hmm_enrichment.R ---
# Purpose: Performs enrichment analysis of genomic features (e.g., genes, enhancers,
#          CpG islands, repeats) within Hidden Markov Model (HMM) chromatin states.
#          Uses Fisher's exact test or hypergeometric test to assess significance.
# Chapter Relevance: 5.4 Chromatin State Analysis

# Load required libraries
library(GenomicRanges)
library(dplyr)
library(stats) # For fisher.test, phyper
library(parallel) # Optional: For parallelizing tests

# --- Configuration & Input --- 

# Define input file paths (Replace with actual paths or loaded objects)
# Assumes the following GRanges objects are loaded:
# - hmm_states: GRanges object with HMM state assignments (e.g., column 'state')
# - features_list: A named list of GRanges objects, where each element represents
#                  a feature type to test for enrichment (e.g., list(Genes = gene_gr, CpG = cpg_gr))
# - background_regions: Optional GRanges object representing the genomic background
#                     (e.g., the whole genome or specific accessible regions).
#                     If NULL, the extent of hmm_states might be used, but defining
#                     an appropriate background is crucial for correct interpretation.

output_dir <- "results/chapter5/"
enrichment_results_csv <- file.path(output_dir, "hmm_state_feature_enrichment.csv")

# Parameters
num_cores <- max(1, detectCores() - 1) # Number of cores for parallel processing

# --- Data Preparation --- 

# Check if input data exists
if (!exists("hmm_states") || !is(hmm_states, "GRanges")) {
    stop("Input GRanges object 'hmm_states' not found or is not a GRanges object.")
}
if (!exists("features_list") || !is.list(features_list) || length(features_list) == 0 || 
    !all(sapply(features_list, is, "GRanges"))) {
    stop("Input 'features_list' must be a non-empty named list of GRanges objects.")
}
if (exists("background_regions") && !is.null(background_regions) && !is(background_regions, "GRanges")) {
    stop("If provided, 'background_regions' must be a GRanges object.")
}
if (!"state" %in% names(mcols(hmm_states))) {
     stop("'hmm_states' GRanges object must have a metadata column named 'state'.")
}

# Define the background set
# Using the full span of HMM states if no specific background is given.
# WARNING: Using HMM states span might inflate enrichment if states are sparse.
# It's generally better to provide a defined background (e.g., mappable genome).
if (!exists("background_regions") || is.null(background_regions)) {
    warning("No specific 'background_regions' provided. Using the union of HMM state regions as background. Enrichment results may be less reliable.")
    background_regions <- reduce(hmm_states)
    # Calculate total background size carefully, sum of widths of reduced regions
    total_background_bases <- sum(width(background_regions))
} else {
    cat("Using provided background regions.\n")
    total_background_bases <- sum(width(background_regions))
}

if (total_background_bases == 0) {
    stop("Total background size is zero. Cannot perform enrichment analysis.")
}

# --- Enrichment Analysis Function --- 

# Function to perform enrichment test for one feature type in one HMM state
calculate_enrichment <- function(state_regions, feature_regions, background_regions, total_background_bases) {
    
    # 1. Overlap between state and feature
    overlap_state_feature <- sum(width(intersect(state_regions, feature_regions)))
    
    # 2. Feature bases within the state but not overlapping the specific feature regions (used in 2x2 table)
    # This calculation might be incorrect for the standard Fisher test setup. Let's rethink.

    # Standard Fisher's Exact Test approach using overlaps:
    # We need counts of regions, not bases, for the standard contingency table approach.
    # Let's switch to region counts. Define 'universe' as background regions.
    
    # Alternative: Base-pair overlap approach (requires careful definition of non-feature/non-state)
    
    # Let's stick to the base-pair overlap Fisher test approach:
    # a = bases in state AND feature (overlap_state_feature)
    # b = bases in state BUT NOT feature
    # c = bases NOT in state BUT IN feature
    # d = bases NOT in state AND NOT in feature

    # Calculate total bases in the current state
    total_state_bases <- sum(width(state_regions))
    
    # Calculate total bases for the current feature within the background
    feature_in_background <- intersect(feature_regions, background_regions)
    total_feature_bases_in_background <- sum(width(feature_in_background))
    
    a <- overlap_state_feature
    b <- total_state_bases - a
    c <- total_feature_bases_in_background - a
    d <- total_background_bases - a - b - c

    # Ensure non-negative values (can happen with slight overlap inconsistencies)
    a <- max(0, a); b <- max(0, b); c <- max(0, c); d <- max(0, d)

    # Create the contingency table
    contingency_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    
    # Expected overlap
    expected_overlap <- (total_state_bases / total_background_bases) * total_feature_bases_in_background
    fold_enrichment <- ifelse(expected_overlap > 0, a / expected_overlap, NA)
    
    # Perform Fisher's exact test
    fisher_result <- tryCatch({
        fisher.test(contingency_table, alternative = "greater") # Test for enrichment
    }, error = function(e) {
        warning(paste("Fisher test failed:", e$message))
        list(p.value = NA, estimate = NA) # Return NAs on error
    })
    
    return(list(
        observed_overlap_bp = a,
        state_bases = total_state_bases,
        feature_bases_in_background = total_feature_bases_in_background,
        background_bases = total_background_bases,
        expected_overlap_bp = expected_overlap,
        fold_enrichment = fold_enrichment,
        p_value = fisher_result$p.value,
        odds_ratio = fisher_result$estimate
    ))
}

# --- Run Enrichment Analysis --- 

all_states <- unique(hmm_states$state)
all_feature_names <- names(features_list)

cat("Starting enrichment analysis for", length(all_states), "states and", length(all_feature_names), "features...\n")

# Create list of tasks for parallel processing
tasks <- list()
for (state_name in all_states) {
    state_regions_subset <- hmm_states[hmm_states$state == state_name]
    # Reduce state regions to avoid double counting bases within the same state
    state_regions_reduced <- reduce(state_regions_subset)
    
    for (feature_name in all_feature_names) {
        feature_regions_subset <- features_list[[feature_name]]
        tasks[[length(tasks) + 1]] <- list(
            state_name = state_name,
            feature_name = feature_name,
            state_regions = state_regions_reduced,
            feature_regions = feature_regions_subset
        )
    }
}

# Run tasks in parallel
# Consider using BiocParallel::bplapply for GRanges compatibility if available
cluster <- makeCluster(num_cores)
clusterExport(cluster, c("calculate_enrichment", "background_regions", "total_background_bases", "intersect", "width", "reduce", "fisher.test")) # Export necessary functions/objects

results_list <- parLapply(cluster, tasks, function(task) {
    enrichment_result <- calculate_enrichment(
        state_regions = task$state_regions,
        feature_regions = task$feature_regions,
        background_regions = background_regions,
        total_background_bases = total_background_bases
    )
    # Add identifiers
    enrichment_result$state <- task$state_name
    enrichment_result$feature <- task$feature_name
    return(enrichment_result)
})

stopCluster(cluster)

# Combine results into a data frame
enrichment_df <- bind_rows(results_list)

# Adjust P-values for multiple testing (across all state-feature combinations)
enrichment_df$fdr <- p.adjust(enrichment_df$p_value, method = "BH")

# Order results
enrichment_df <- enrichment_df %>% 
    arrange(fdr, p_value)

# --- Output Results --- 

cat("Enrichment analysis complete.\n")
print(head(enrichment_df))

# Save results to CSV
write.csv(enrichment_df, enrichment_results_csv, row.names = FALSE)
cat("Saved enrichment results to:", enrichment_results_csv, "\n")


# --- Optional: Visualization (e.g., Heatmap of Fold Enrichment) ---
# library(ggplot2)
# library(reshape2)
# 
# # Example: Heatmap of log2 fold enrichment
# heatmap_data <- enrichment_df %>% 
#                 mutate(log2_fold_enrichment = log2(fold_enrichment)) %>% 
#                 # Handle Inf/-Inf/NA from log2(0) or NA fold_enrichment
#                 mutate(log2_fold_enrichment = ifelse(is.finite(log2_fold_enrichment), log2_fold_enrichment, 0)) %>% 
#                 # Cap extreme values for better visualization
#                 mutate(log2_fold_enrichment = pmax(-4, pmin(4, log2_fold_enrichment))) %>% 
#                 select(state, feature, log2_fold_enrichment)
# 
# # Add significance stars based on FDR
# heatmap_data <- heatmap_data %>% 
#                 left_join(enrichment_df %>% select(state, feature, fdr), by = c("state", "feature")) %>%
#                 mutate(significance = case_when(
#                                         fdr < 0.001 ~ "***",
#                                         fdr < 0.01  ~ "**",
#                                         fdr < 0.05  ~ "*",
#                                         TRUE ~ ""
#                                     ))
# 
# p_heatmap <- ggplot(heatmap_data, aes(x = feature, y = state, fill = log2_fold_enrichment)) + 
#     geom_tile(color = "grey90") + 
#     geom_text(aes(label = significance), size = 3, color = "black") + # Add significance stars
#     scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Log2(Fold Enrichment)") + 
#     labs(title = "HMM State Enrichment Analysis", x = "Genomic Feature", y = "HMM State") + 
#     theme_minimal() + 
#     theme(axis.text.x = element_text(angle = 45, hjust = 1), 
#           plot.title = element_text(hjust = 0.5))
# 
# print(p_heatmap)
# ggsave(file.path(output_dir, "hmm_enrichment_heatmap.pdf"), p_heatmap, width = 8, height = max(4, length(all_states)*0.5))
# cat("Saved enrichment heatmap (optional visualization).\n")

cat("--- 03_hmm_enrichment.R finished ---\n") 