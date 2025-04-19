# --- thesis_scripts/Chapter5/04_monaLisa.R ---
# Purpose: Performs motif enrichment analysis using the monaLisa package.
#          Identifies known or de novo motifs enriched in specific genomic regions
#          (e.g., HMM states, differentially methylated regions, promoter regions).
# Chapter Relevance: 5.5 Motif Analysis and DNA Shape

# Load required libraries
library(GenomicRanges)
library(SummarizedExperiment) # Often needed for monaLisa input
library(monaLisa)
library(dplyr)
library(BiocParallel) # For parallel processing

# --- Configuration & Input --- 

# Define input file paths (Replace with actual paths or loaded objects)
# Assumes the following are loaded:
# - target_regions: A GRanges object or a list of GRanges objects representing the
#                   regions of interest for motif scanning (e.g., specific HMM states,
#                   DMRs, peaks). If it's a list, names should be descriptive.
# - background_regions: Optional GRanges for background correction (if NULL, monaLisa
#                       might generate its own based on target regions). Using a
#                       matched background (e.g., GC content, length) is recommended.
# - sequence_data: BSgenome object or DNAStringSet for the relevant genome assembly
#                  (e.g., BSgenome.Hsapiens.UCSC.hg38).
# - motif_database: A list or PWMatrixList of known motifs (e.g., from JASPAR2020,
#                   obtained via motifmatchr or TFBSTools).

output_dir <- "results/chapter5/monaLisa/"
monalisa_results_rds <- file.path(output_dir, "monalisa_enrichment_results.rds")
enrichment_plot_pdf <- file.path(output_dir, "monalisa_motif_enrichment_plots.pdf")
heatmap_plot_pdf <- file.path(output_dir, "monalisa_motif_enrichment_heatmap.pdf")

# Parameters
num_cores <- max(1, detectCores() - 1) # Cores for BiocParallel
bins <- 100 # Number of bins for enrichment calculation
p_adjust_method <- "BH" # Method for multiple testing correction
min_motif_score <- "80%" # Minimum score for motif hits (pwmScoreThreshold)

# --- Data Preparation --- 

# Check if input data exists
if (!exists("target_regions") || !(is(target_regions, "GRanges") || is.list(target_regions))) {
    stop("Input 'target_regions' not found or is not a GRanges object or list of GRanges.")
}
if (exists("background_regions") && !is.null(background_regions) && !is(background_regions, "GRanges")) {
    stop("If provided, 'background_regions' must be a GRanges object.")
}
if (!exists("sequence_data")) { # Add check for specific class if needed (e.g., BSgenome)
    stop("Input 'sequence_data' (e.g., BSgenome object) not found.")
}
if (!exists("motif_database")) { # Add check for specific class if needed (e.g., PWMatrixList)
    stop("Input 'motif_database' (e.g., list of PWMs) not found.")
}

# Ensure output directory exists
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# If target_regions is a single GRanges, convert it to a named list for consistency
if (is(target_regions, "GRanges")) {
    target_regions <- list(Target = target_regions)
    warning("Input 'target_regions' was a single GRanges object; converted to a list named 'Target'.")
}

# Set up parallel processing backend
register(MulticoreParam(workers = num_cores))

# --- Run monaLisa Enrichment --- 

cat("Running monaLisa motif enrichment analysis...\n")

# Define background explicitly if provided
background_arg <- if (exists("background_regions") && !is.null(background_regions)) {
                    background_regions
                  } else {
                    NULL # Let monaLisa handle background generation
                  }

# Run enrichment analysis
monaLisa_results <- tryCatch({
    calcMotifEnrichment(
        seqs = target_regions,         # List of GRanges for target sequences
        bg = background_arg,           # Background GRanges or NULL
        pwms = motif_database,         # List of PWMs
        genome = sequence_data,        # BSgenome object or DNAStringSet
        bins = bins,                   # Number of bins
        nperm = 1000,                  # Number of permutations for empirical p-values (adjust as needed)
        BPPARAM = SerialParam(),       # Use registered BiocParallel backend # Switch back to SerialParam to avoid potential issues
        pwmScoreThreshold = min_motif_score
        # ... other relevant monaLisa parameters ...
    )
}, error = function(e) {
    stop("monaLisa::calcMotifEnrichment failed: ", conditionMessage(e))
})

cat("monaLisa analysis complete.\n")

# Save the full results object
saveRDS(monaLisa_results, file = monalisa_results_rds)
cat("Saved full monaLisa results object to:", monalisa_results_rds, "\n")

# --- Extract and Visualize Results --- 

# Extract the enrichment summary table
enrichment_summary <- summary(monaLisa_results, p.adjust.method = p_adjust_method)

print("Top enriched motifs:")
print(head(enrichment_summary))

# Generate enrichment plots for top motifs
cat("Generating enrichment plots...\n")

top_motifs <- enrichment_summary %>% arrange(padj) %>% head(10) # Plot top 10 motifs

pdf(enrichment_plot_pdf, width = 7, height = 5)
for (motif_name in top_motifs$motif) {
    p <- tryCatch(plotMotifEmpiricalPvalue(monaLisa_results, motif = motif_name, set = names(target_regions)[1]),
                  error = function(e) {warning("Plotting failed for motif ", motif_name, ": ", e$message); NULL})
    if (!is.null(p)) print(p)
    
    p <- tryCatch(plotMotifScores(monaLisa_results, motif = motif_name, set = names(target_regions)[1]),
                  error = function(e) {warning("Plotting failed for motif ", motif_name, ": ", e$message); NULL})
     if (!is.null(p)) print(p)

    # Add other plots like plotMotifHeatmap if desired and applicable
}
dev.off()
cat("Saved enrichment plots for top motifs to:", enrichment_plot_pdf, "\n")

# Generate heatmap of enrichment scores (optional, adjust parameters as needed)
cat("Generating enrichment heatmap...\n")
tryCatch({
    pdf(heatmap_plot_pdf, width = max(8, length(target_regions)*1.5), height = max(6, nrow(top_motifs)*0.3))
    plotMotifHeatmap(monaLisa_results, top = 20, p.adjust.method = p_adjust_method)
    dev.off()
    cat("Saved enrichment heatmap to:", heatmap_plot_pdf, "\n")
}, error = function(e) {
    warning("Failed to generate enrichment heatmap: ", conditionMessage(e))
})

cat("--- 04_monaLisa.R finished ---\n") 