# --- thesis_scripts/Chapter5/06_gviz_plotting.R ---
# Purpose: Generates genomic visualizations using the Gviz package.
#          Combines various data tracks (e.g., gene models, methylation, HMM states,
#          motif sites, experimental data) for specific genomic regions of interest.
# Chapter Relevance: Integrates and visualizes findings from previous Chapter 5 analyses.

# Load required libraries
library(Gviz)
library(GenomicRanges)
library(rtracklayer) # Potentially needed for importing BED/GFF/BigWig tracks
library(dplyr)

# --- Configuration & Input --- 

# Define input file paths and parameters (Replace with actual paths or loaded objects)
# Assumes various GRanges or data frame objects representing tracks are loaded.
# Examples:
# - gene_models_gr: GRanges with gene/transcript structures (e.g., from TxDb object)
# - methylation_gr: GRanges with methylation scores (e.g., beta values in a metadata col)
# - hmm_states_gr: GRanges with HMM state assignments (metadata col 'state')
# - motif_sites_gr: GRanges with locations of specific motif instances
# - regions_to_plot: A data frame or GRanges specifying the genomic coordinates
#                    (chr, start, end) of regions to visualize.
# - genome_assembly: String identifier for the genome assembly (e.g., "hg38", "hg19")

output_dir <- "results/chapter5/gviz_plots/"
plot_file_prefix <- file.path(output_dir, "genomic_region_plot") # Prefix for output files

# Plotting parameters
extend_plot_window <- 5000 # Extend region boundaries by this amount (bp)

# --- Data Preparation & Track Creation --- 

# Check if input data exists
if (!exists("gene_models_gr")) stop("Input 'gene_models_gr' not found.")
if (!exists("regions_to_plot")) stop("Input 'regions_to_plot' specifying plot regions not found.")
if (!exists("genome_assembly")) stop("Variable 'genome_assembly' specifying the genome is not defined.")
# Add checks for other optional track data if they are intended to be plotted

# Ensure output directory exists
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Create Gviz tracks (examples - customize based on available data)

# 1. Genome Axis Track
gtrack <- GenomeAxisTrack()

# 2. Gene Region Track
# Ensure gene_models_gr has appropriate columns like 'gene', 'transcript', 'symbol'
# These might need to be fetched from an annotation database if not present.
geneTrack <- GeneRegionTrack(gene_models_gr, 
                             genome = genome_assembly, 
                             name = "Genes", 
                             transcriptAnnotation = "symbol", # Show gene symbols
                             background.title = "#440154FF", # Viridis color
                             col.title = "white",
                             showId = TRUE, geneSymbol = TRUE)
                             

# 3. Methylation Data Track (Example using DataTrack)
# Assumes methylation_gr has a 'score' column with methylation levels
if (exists("methylation_gr") && "score" %in% names(mcols(methylation_gr))) {
    methTrack <- DataTrack(range = methylation_gr, 
                           genome = genome_assembly,
                           name = "Methylation", 
                           type = "p", # Plot points
                           cex = 0.5,
                           pch = 16,
                           ylim = c(0, 1), # Assuming beta values 0-1
                           background.title = "#31688EFF") # Viridis color
} else {
    methTrack <- NULL
    warning("Methylation data track not created. 'methylation_gr' with 'score' column not found.")
}

# 4. HMM States Track (Example using AnnotationTrack)
# Assumes hmm_states_gr has a 'state' column
if (exists("hmm_states_gr") && "state" %in% names(mcols(hmm_states_gr))) {
    # Define colors for states (customize based on your states)
    state_colors <- setNames(RColorBrewer::brewer.pal(length(unique(hmm_states_gr$state)), "Set3"), 
                             unique(hmm_states_gr$state))
    hmmTrack <- AnnotationTrack(range = hmm_states_gr, 
                                genome = genome_assembly, 
                                name = "HMM States", 
                                id = hmm_states_gr$state, # Use state name as ID
                                # shape = "box", # Default shape
                                fill = state_colors[hmm_states_gr$state],
                                stacking = "dense",
                                background.title = "#35B779FF") # Viridis color
} else {
    hmmTrack <- NULL
    warning("HMM states track not created. 'hmm_states_gr' with 'state' column not found.")
}

# 5. Motif Sites Track (Example using AnnotationTrack)
# Assumes motif_sites_gr has a 'name' or 'motif_id' column
if (exists("motif_sites_gr")) {
    motifTrack <- AnnotationTrack(range = motif_sites_gr, 
                                  genome = genome_assembly, 
                                  name = "Motif Sites",
                                  # id = motif_sites_gr$name, # Optional: show motif names
                                  shape = "arrow", # Use arrows to show strand if available
                                  fill = "orange",
                                  stacking = "squish",
                                  background.title = "#FDE725FF", # Viridis color
                                  col.title="black")
} else {
    motifTrack <- NULL
    warning("Motif sites track not created. 'motif_sites_gr' not found.")
}

# --- Plotting Function --- 

# Function to plot a specific genomic region with selected tracks
plot_genomic_region <- function(chromosome, start_pos, end_pos, 
                                all_tracks, filename, 
                                title = paste(chromosome, start_pos, end_pos, sep = ":")) {
    
    cat("Plotting region:", title, "to file:", filename, "\n")
    
    # Filter out NULL tracks
    tracks_to_plot <- Filter(Negate(is.null), all_tracks)
    
    if (length(tracks_to_plot) < 2) { # Need at least axis + one data track
        warning("Not enough valid tracks to plot region: ", title)
        return()
    }
    
    # Define plot dimensions (adjust as needed)
    plot_heights <- c(0.5, 2, rep(1, length(tracks_to_plot) - 2)) # Relative heights: axis, genes, others
    plot_width_cm <- 20 # Adjust width
    plot_height_cm <- sum(plot_heights) * 1.5 # Adjust height multiplier
    
    # Save plot to PDF
    pdf(filename, width = plot_width_cm / 2.54, height = plot_height_cm / 2.54) # Convert cm to inches
    tryCatch({
        plotTracks(tracks_to_plot, 
                   chromosome = chromosome, 
                   from = start_pos, 
                   to = end_pos, 
                   main = title,
                   sizes = plot_heights, # Control relative track heights
                   background.panel = "#F7F7F7",
                   fontcolor.title = "white",
                   col.axis="black", fontcolor.axis="black"
                   # Add more Gviz styling options here
                   )
    }, error = function(e) {
        warning("Failed to plot region ", title, ": ", conditionMessage(e))
        plot.new() # Create blank plot on error to avoid corrupt PDF
        title(main = paste("Error plotting region:", title), col.main = "red")
    })
    dev.off()
}

# --- Generate Plots for Regions of Interest --- 

cat("\nGenerating Gviz plots for specified regions...\n")

# Combine all tracks into a list (order matters for plotting)
all_available_tracks <- list(gtrack, geneTrack, methTrack, hmmTrack, motifTrack)

# Iterate through the regions defined in `regions_to_plot`
# Assumes regions_to_plot has columns: chr, start, end, name (optional)
if (nrow(regions_to_plot) > 0) {
    for (i in 1:nrow(regions_to_plot)) {
        region_chr <- regions_to_plot$chr[i]
        region_start <- regions_to_plot$start[i] - extend_plot_window
        region_end <- regions_to_plot$end[i] + extend_plot_window
        region_name <- ifelse("name" %in% colnames(regions_to_plot), 
                              gsub("[^[:alnum:]_.-]", "_", regions_to_plot$name[i]), # Sanitize name
                              paste0("region_", i))
        
        plot_filename <- paste0(plot_file_prefix, "_", region_name, ".pdf")
        plot_title <- paste("Region:", region_name, "(", region_chr, ":", regions_to_plot$start[i], "-", regions_to_plot$end[i], ")")
        
        # Ensure start is not negative
        region_start <- max(1, region_start)
        
        plot_genomic_region(chromosome = region_chr, 
                            start_pos = region_start, 
                            end_pos = region_end, 
                            all_tracks = all_available_tracks,
                            filename = plot_filename,
                            title = plot_title)
    }
} else {
    cat("No regions specified in 'regions_to_plot'. No plots generated.\n")
}

cat("--- 06_gviz_plotting.R finished ---\n") 