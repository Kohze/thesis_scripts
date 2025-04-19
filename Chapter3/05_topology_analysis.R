# --- thesis_scripts/Chapter3/05_topology_analysis.R ---
# Purpose: Analyze nORF structural topology based on predicted PDB structures
#          and DSSP secondary structure assignments. Calculates solvent accessible area,
#          Ramachandran angles, and secondary structure element proportions.
#          Includes Normal Mode Analysis (NMA) functions (execution part commented out).
# Chapter Relevance: 3.6 Topology Modelling

# Load necessary packages
library(bio3d)      # For reading PDB, DSSP, NMA, Ramachandran plots
library(dplyr)
library(tidyr)      # For gather
library(ggplot2)
library(ggstatsplot) # For ggbetweenstats

# Define Input/Output Directories
pdb_input_dir <- "./pdb"       # Directory containing PDB files (e.g., from ESMFold)
dssp_input_dir <- "./dssp"      # Directory containing DSSP output files for the PDBs
output_dir <- "."

# --- Normal Mode Analysis (NMA) --- 
# Functions for NMA - Note: Running NMA can be computationally intensive.

# Example function from original script to calculate flexibility index
# (Note: This function definition seems unusual for IFlex, typically IFlex is calculated differently in bio3d context)
# A more standard approach uses fluctuations directly.
IFlex_original <- function(weights) {
    N <- length(weights)
    if (N == 0) return(NA)
    numerator <- sum(weights^2, na.rm = TRUE)
    denominator <- N
    IFlex <- sqrt(numerator / denominator)
    return(IFlex)
}

# --- Optional: Run NMA on PDB files --- 
# WARNING: This loop can take a very long time and may fail for certain structures.
# It requires significant computational resources.
run_nma_analysis <- FALSE # Set to TRUE to attempt NMA calculation

if (run_nma_analysis) {
    if (!dir.exists(pdb_input_dir)) {
        warning("PDB input directory not found: ", pdb_input_dir, ". Skipping NMA.")
    } else {
        cat("Starting NMA analysis (this may take a long time)...\n")
        all_pdb_files <- list.files(pdb_input_dir, pattern = "\\.pdb$", full.names = TRUE)
        all_nma_results <- list()

        # Consider using parallel processing (e.g., future.apply) for speedup
        for (pdb_file in all_pdb_files) {
            pdb_id <- tools::file_path_sans_ext(basename(pdb_file))
            cat(" Processing NMA for:", pdb_id, "\n")
            pdb_data <- tryCatch(read.pdb(pdb_file), error = function(e) {
                warning("Error reading PDB ", pdb_id, ": ", conditionMessage(e))
                NULL
            })

            if (!is.null(pdb_data)) {
                nma_result <- tryCatch(nma(pdb_data), error = function(e) {
                    warning("Error running NMA for ", pdb_id, ": ", conditionMessage(e))
                    NULL
                })

                if (!is.null(nma_result)) {
                    all_nma_results[[pdb_id]] <- nma_result
                    # Example: Calculate and store fluctuations
                    # fluctuations <- nma_result$fluctuations
                    # flexibility <- sd(fluctuations) # Example flexibility metric
                    # cat("  NMA Frequencies:", paste(round(nma_result$frequencies, 2), collapse=", "), "\n")
                } else {
                    cat("  NMA failed for", pdb_id, "\n")
                }
            } else {
                 cat("  Skipping NMA due to PDB read error for", pdb_id, "\n")
            }
            # Add progress reporting if running on many files
        }
        cat("NMA analysis finished.\n")
        # Optional: Save NMA results
        # save(all_nma_results, file = file.path(output_dir, "nma_results.RData"))
    }
} else {
    cat("Skipping NMA analysis as run_nma_analysis is FALSE.\n")
}


# --- DSSP Analysis (Secondary Structure) --- 
cat("\nStarting DSSP analysis from directory:", dssp_input_dir, "\n")

if (!dir.exists(dssp_input_dir)) {
    stop("DSSP input directory not found: ", dssp_input_dir)
}

dssp_files <- list.files(dssp_input_dir, pattern = "\\.dssp$")
if (length(dssp_files) == 0) {
    stop("No DSSP files found in ", dssp_input_dir)
}

# Create an empty dataframe to store the summary information
dssp_summary_df <- data.frame(
    id = character(),
    sol_area_mean = numeric(),
    sol_area_sd = numeric(),
    phi_mean = numeric(),
    psi_mean = numeric(),
    prop_turn = numeric(),    # Proportion Turn (T, S)
    prop_helix = numeric(),   # Proportion Helix (H, G, I) - Grouping G/I with H based on common practice
    prop_strand = numeric(),  # Proportion Beta Strand (E, B)
    prop_coil = numeric(),    # Proportion Coil (represented by ' ' in dssp$sse)
    stringsAsFactors = FALSE
)

# Dataframe to store all phi/psi angles for combined Ramachandran plot
all_phi_psi_df <- data.frame(phi = numeric(), psi = numeric())

# Loop through all DSSP files in the folder
processed_count <- 0
error_count <- 0
for (file in dssp_files) {
    pdb_name <- tools::file_path_sans_ext(file)
    dssp_result <- tryCatch({
        dssp_data <- read.dssp(file.path(dssp_input_dir, file))
        
        # Basic checks
        if(is.null(dssp_data$acc) || is.null(dssp_data$phi) || is.null(dssp_data$psi) || is.null(dssp_data$sse)){
            warning("Missing expected columns in DSSP object for: ", file)
            return(NULL)
        }
        
        # Extract solvent accessible area stats
        sol_area_mean <- mean(dssp_data$acc, na.rm = TRUE)
        sol_area_sd <- sd(dssp_data$acc, na.rm = TRUE)

        # Extract Ramachandran angles stats
        phi_mean <- mean(dssp_data$phi, na.rm = TRUE)
        psi_mean <- mean(dssp_data$psi, na.rm = TRUE)

        # Store all valid phi/psi pairs for combined plot
        valid_angles <- !is.na(dssp_data$phi) & !is.na(dssp_data$psi)
        if(sum(valid_angles) > 0) {
             all_phi_psi_df <- rbind(all_phi_psi_df, data.frame(phi = dssp_data$phi[valid_angles], psi = dssp_data$psi[valid_angles]))
        }
       
        # Calculate secondary structure proportions
        sse <- dssp_data$sse
        total_count <- length(sse)
        if(total_count == 0) return(NULL) # Skip if no SSE info

        prop_turn <- sum(sse %in% c("T", "S"), na.rm=TRUE) / total_count
        prop_helix <- sum(sse %in% c("H", "G", "I"), na.rm=TRUE) / total_count
        prop_strand <- sum(sse %in% c("E", "B"), na.rm=TRUE) / total_count
        prop_coil <- sum(sse == " ", na.rm=TRUE) / total_count # Coil often represented by space

        # Append the information to the summary dataframe
        summary_row <- data.frame(
            id = pdb_name,
            sol_area_mean = sol_area_mean,
            sol_area_sd = sol_area_sd,
            phi_mean = phi_mean,
            psi_mean = psi_mean,
            prop_turn = prop_turn,
            prop_helix = prop_helix,
            prop_strand = prop_strand,
            prop_coil = prop_coil
        )
        return(summary_row)

    }, error = function(e) {
        warning("Error processing DSSP file ", file, ": ", conditionMessage(e))
        return(NULL)
    })

    if (!is.null(dssp_result)) {
        dssp_summary_df <- rbind(dssp_summary_df, dssp_result)
        processed_count <- processed_count + 1
    } else {
        error_count <- error_count + 1
    }
     if ((processed_count + error_count) %% 100 == 0) {
         cat(" Processed", processed_count, "DSSP files (", error_count, "errors )...\n")
     }
}

cat("Finished DSSP processing. Successfully processed:", processed_count, "Errored:", error_count, "\n")

# Save the summary data
summary_output_file <- file.path(output_dir, "dssp_summary_stats.csv")
cat("Saving DSSP summary stats to:", summary_output_file, "\n")
fwrite(as.data.table(dssp_summary_df), summary_output_file)

# --- Ramachandran Plots --- 
cat("Generating Ramachandran plots...\n")

# Plot 1: Averaged Phi/Psi (less informative but present in original script)
if (nrow(dssp_summary_df) > 0) {
    avg_phi_psi_matrix <- as.matrix(dssp_summary_df[, c("phi_mean", "psi_mean")])
    png(file.path(output_dir, "ramachandran_average.png"), width = 800, height = 800)
    tryCatch(ramachandran(avg_phi_psi_matrix, xBins = 150, yBins = 150,
                          main = "Averaged nORF Ramachandran Plot (Mean Phi/Psi per structure)"),
             error=function(e){ plot.new(); title("Error plotting average Ramachandran") })
    dev.off()
} else {
    cat ("Skipping average Ramachandran plot - no summary data.\n")
}

# Plot 2: Combined Ramachandran plot from all residues (more standard)
if (nrow(all_phi_psi_df) > 0) {
    phi_psi_matrix <- as.matrix(all_phi_psi_df)
    # Subsample if too many points for performance
    if (nrow(phi_psi_matrix) > 500000) {
         cat("Subsampling Ramachandran points for plotting (original points:", nrow(phi_psi_matrix), ")...\n")
         set.seed(456)
         sample_indices <- sample(1:nrow(phi_psi_matrix), 500000)
         phi_psi_matrix <- phi_psi_matrix[sample_indices, ]
    }
    png(file.path(output_dir, "ramachandran_combined.png"), width = 800, height = 800)
    tryCatch(ramachandran(phi_psi_matrix, xBins = 150, yBins = 150,
                          main = "Combined nORF Ramachandran Plot (All Residues)"),
             error=function(e){ plot.new(); title("Error plotting combined Ramachandran") })
    dev.off()
} else {
    cat("Skipping combined Ramachandran plot - no angle data collected.\n")
}

cat("Ramachandran plots saved.\n")

# --- Secondary Structure Distribution Plot --- 
cat("Generating secondary structure distribution plot...\n")
if (nrow(dssp_summary_df) > 0 && ncol(dssp_summary_df) > 5) { # Check if SSE columns exist
    # Prepare data for ggbetweenstats (long format)
    sse_cols <- c("prop_turn", "prop_helix", "prop_strand", "prop_coil")
    if (all(sse_cols %in% names(dssp_summary_df))) {
        sse_data_long <- dssp_summary_df %>% 
            select(id, all_of(sse_cols)) %>% 
            pivot_longer(cols = all_of(sse_cols), names_to = "sse_type", values_to = "proportion") %>% 
            mutate(sse_type = gsub("prop_", "", sse_type)) # Clean names for plot

        # Create plot
        p_sse_dist <- ggbetweenstats(
          data  = sse_data_long,
          x     = sse_type,
          y     = proportion,
          title = "nORF Secondary Structure Element Proportions (per Structure)",
          xlab = "Secondary Structure Element",
          ylab = "Proportion"
          # Add pairwise comparisons if desired
          # pairwise.comparisons = TRUE, 
          # pairwise.display = "significant", 
          # p.adjust.method = "BH"
        )
        
        print(p_sse_dist)
        ggsave(file.path(output_dir, "secondary_structure_distribution.png"), plot = p_sse_dist, width = 8, height = 6)
        cat("Saved secondary structure distribution plot.\n")
        
    } else {
         cat("Skipping SSE distribution plot - missing required proportion columns.\n")
    }
} else {
    cat("Skipping SSE distribution plot - no summary data or SSE columns.\n")
}

cat("--- 05_topology_analysis.R finished ---\n") 