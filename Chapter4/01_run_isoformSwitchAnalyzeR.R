# --- run_isoformSwitchAnalyzeR.R ---
# Purpose: Performs Differential Transcript Usage (DTU) analysis using the
#          isoformSwitchAnalyzeR package for brain development time course data.
#          This script processes data for Forebrain, Midbrain, and Hindbrain separately.
#          Placeholders need to be customized by the user.
# Chapter Relevance: Chapter 4 (Generates core DTU results)

# --- 1. Load Libraries ---
print("Loading necessary libraries...")
suppressPackageStartupMessages({
  library(isoformSwitchAnalyzeR)
  library(dplyr)
  library(Biostrings)
  library(ggplot2)
  # Add other libraries if needed, e.g., for specific external tool integration
})

# --- 2. Define Paths and Parameters (USER MUST CUSTOMIZE) ---
print("Defining paths and parameters...")

# Input Data Paths
kallisto_path_template <- "/path/to/kallisto_quant/{sample_id}" # Template for Kallisto output dirs
annotation_path        <- "/path/to/your/annotation.gtf"       # Path to GTF/GFF3 annotation file
fasta_path             <- "/path/to/your/genome.fasta"         # Path to genome FASTA file

# External Tool Paths (if using sequence analysis features)
cpat_path      <- "/path/to/CPAT/Installation/cpat.py"        # Path to CPAT script
pfam_path      <- "/path/to/pfam_scan.pl"                     # Path to pfam_scan.pl script
signalp_path   <- "/path/to/signalp-5.0/bin/signalp"          # Path to SignalP executable (adjust version)
netsurfp2_path <- "/path/to/netsurfP-2/run_netsurfP-2.sh"     # Path to NetSurfP-2 script

# Output Directory
output_dir <- "./thesis_scripts/chapter4/isoformSwitchAnalyzeR_output" # Directory to save results
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  print(paste("Created output directory:", output_dir))
}

# Analysis Parameters
analysis_threads   <- 4          # Number of cores for parallel processing
alpha_value        <- 0.05       # Significance threshold for switch testing
dIF_cutoff         <- 0.1        # Minimum delta Isoform Fraction cutoff
expression_cutoff  <- 1          # Minimum TPM expression cutoff for filtering
replicate_count_cutoff <- 2      # Minimum number of replicates a transcript needs expression in

# Define which sequence analysis features to run (Set to FALSE if external tools are not installed/configured)
run_orf_analysis    <- TRUE
run_cpat_analysis   <- FALSE # Requires CPAT installation
run_signalp_analysis<- FALSE # Requires SignalP installation
run_pfam_analysis   <- FALSE # Requires PfamScan and Pfam DB installation
run_netsurfp2_analysis <- FALSE # Requires NetSurfP-2 installation

# --- 3. Define Experimental Design (USER MUST CUSTOMIZE) ---
# Create a data frame mapping samples to conditions and Kallisto paths.
# Example structure - adapt to your actual sample IDs and timepoints/conditions.
# Ensure 'sampleID' column matches the directory names in kallisto_path_template.
# 'condition' column defines the groups being compared (e.g., timepoints, tissues).
# Add a 'tissue' column to loop through different brain regions.

print("Defining experimental design...")
sample_data <- data.frame(
  sampleID = c(
    # Forebrain
    "FB_E10.5_rep1", "FB_E10.5_rep2", "FB_E11.5_rep1", "FB_E11.5_rep2",
    "FB_E12.5_rep1", "FB_E12.5_rep2", "FB_E13.5_rep1", "FB_E13.5_rep2",
    "FB_E14.5_rep1", "FB_E14.5_rep2", "FB_E15.5_rep1", "FB_E15.5_rep2",
    "FB_E16.5_rep1", "FB_E16.5_rep2", "FB_P0_rep1", "FB_P0_rep2",
    # Midbrain
    "MB_E10.5_rep1", "MB_E10.5_rep2", "MB_E11.5_rep1", "MB_E11.5_rep2",
    "MB_E12.5_rep1", "MB_E12.5_rep2", "MB_E13.5_rep1", "MB_E13.5_rep2",
    "MB_E14.5_rep1", "MB_E14.5_rep2", "MB_E15.5_rep1", "MB_E15.5_rep2",
    "MB_E16.5_rep1", "MB_E16.5_rep2", "MB_P0_rep1", "MB_P0_rep2",
    # Hindbrain
    "HB_E10.5_rep1", "HB_E10.5_rep2", "HB_E11.5_rep1", "HB_E11.5_rep2",
    "HB_E12.5_rep1", "HB_E12.5_rep2", "HB_E13.5_rep1", "HB_E13.5_rep2",
    "HB_E14.5_rep1", "HB_E14.5_rep2", "HB_E15.5_rep1", "HB_E15.5_rep2",
    "HB_E16.5_rep1", "HB_E16.5_rep2", "HB_P0_rep1", "HB_P0_rep2"
  ),
  stringsAsFactors = FALSE
)

# Add tissue and condition (timepoint) columns based on sampleID
sample_data <- sample_data %>%
  mutate(
    tissue = case_when(
      grepl("^FB", sampleID) ~ "Forebrain",
      grepl("^MB", sampleID) ~ "Midbrain",
      grepl("^HB", sampleID) ~ "Hindbrain",
      TRUE ~ NA_character_
    ),
    condition = sub("^(FB|MB|HB)_|_rep\d$", "", sampleID) # Extract timepoint as condition
  )

# Add full Kallisto paths
sample_data$path <- file.path(dirname(kallisto_path_template), sample_data$sampleID)

# Verify the sample_data structure
print("Sample data preview:")
print(head(sample_data))
print("Tissue summary:")
print(table(sample_data$tissue))

# Define the tissues to process
tissues_to_process <- unique(sample_data$tissue)

# --- 4. Loop Through Tissues and Run Analysis ---
all_results_list <- list() # Store results for each tissue

for (current_tissue in tissues_to_process) {
  print(paste(">>> Starting analysis for:", current_tissue))

  # Subset sample data for the current tissue
  tissue_sample_data <- sample_data %>% filter(tissue == current_tissue)

  # Verify subset
  if (nrow(tissue_sample_data) == 0) {
    warning(paste("No samples found for tissue:", current_tissue, "- Skipping."))
    next
  }
  print(paste("Processing", nrow(tissue_sample_data), "samples for", current_tissue))

  # --- 4.1 Import Kallisto Data ---
  print(paste("Importing Kallisto data for", current_tissue, "..."))
  # Check if paths exist before importing
  existing_paths <- tissue_sample_data$path[sapply(tissue_sample_data$path, dir.exists)]
  if (length(existing_paths) < nrow(tissue_sample_data)) {
      warning(paste("Missing Kallisto directories for:", paste(setdiff(tissue_sample_data$path, existing_paths), collapse=", ")))
  }
  if (length(existing_paths) == 0) {
      warning(paste("No valid Kallisto directories found for", current_tissue, "- Skipping import."))
      next
  }

  # Subset data for existing paths only
  tissue_sample_data_filtered <- tissue_sample_data %>% filter(path %in% existing_paths)

  kallisto_abundance <- tryCatch({
      importIsoformExpression(
          sampleVector = tissue_sample_data_filtered$path,
          addIsofomIdAsColumn = TRUE # Recommended for joining later
      )
  }, error = function(e) {
      message(paste("Error importing Kallisto data for", current_tissue, ":", e$message))
      NULL
  })

  if (is.null(kallisto_abundance)) next # Skip if import failed

  # --- 4.2 Create switchAnalyzeRlist ---
  print(paste("Creating switchAnalyzeRlist for", current_tissue, "..."))
  switchList <- tryCatch({
    importRdata(
      isoformCountMatrix   = kallisto_abundance$counts,
      isoformRepExpression = kallisto_abundance$abundance,
      designMatrix         = tissue_sample_data_filtered[, c("sampleID", "condition")],
      isoformExonAnnoation = annotation_path,
      isoformNtFasta       = fasta_path,
      comparisonsToMake = NULL, # Define comparisons later if needed, or let DEXSeq handle all pairs
      showProgress = TRUE
    )
  }, error = function(e) {
    message(paste("Error creating switchAnalyzeRlist for", current_tissue, ":", e$message))
    NULL
  })

  if (is.null(switchList)) next # Skip if creation failed

  print("Initial switchAnalyzeRlist summary:")
  print(summary(switchList))

  # --- 4.3 Pre-filtering ---
  print(paste("Pre-filtering data for", current_tissue, "..."))
  switchList_filtered <- tryCatch({
    preFilter(
      switchAnalyzeRlist = switchList,
      geneExpressionCutoff = expression_cutoff,     # Min gene TPM
      isoformExpressionCutoff = expression_cutoff,  # Min isoform TPM
      IFcutoff = 0.01,             # Min isoform fraction (prevents issues with very low IF isoforms)
      removeSingleIsoformGenes = TRUE,
      reduceToSwitchingGenes = FALSE, # Keep all genes initially
      filterIsoformsWithoutReplicates = TRUE,
      minNrIsoforms = 2, # Genes must have at least 2 isoforms to be considered for switching
      nCutoff = replicate_count_cutoff # Isoforms must be expressed in at least this many replicates per condition
    )
  }, error = function(e) {
    message(paste("Error during pre-filtering for", current_tissue, ":", e$message))
    NULL
  })

   if (is.null(switchList_filtered)) {
     warning(paste("Pre-filtering failed or removed all data for", current_tissue, ". Saving unfiltered list."))
     # Save the unfiltered list if filtering failed drastically
     saveRDS(switchList, file = file.path(output_dir, paste0(current_tissue, "_switchList_unfiltered.rds")))
     next
   }

  print("SwitchAnalyzeRlist summary after filtering:")
  print(summary(switchList_filtered))
  if (nrow(switchList_filtered$isoformFeatures) == 0) {
     warning(paste("No data remaining after filtering for", current_tissue, "- Skipping further analysis."))
     next
  }

  # --- 4.4 Statistical Analysis (DEXSeq) ---
  # Compares IF values across all conditions defined in the design matrix
  print(paste("Running DEXSeq test for", current_tissue, "..."))
  switchList_tested <- tryCatch({
    isoformSwitchTestDEXSeq(
      switchAnalyzeRlist = switchList_filtered,
      alpha = alpha_value,   # Significance level for screening step
      dIFcutoff = dIF_cutoff, # Significance level for stage-wise analysis
      reduceFurtherToGenesWithConsequences = FALSE, # Keep all significant switches initially
      nCores = analysis_threads,
      showProgress = TRUE
    )
  }, error = function(e) {
    message(paste("Error running isoformSwitchTestDEXSeq for", current_tissue, ":", e$message))
    NULL
  })

  if (is.null(switchList_tested)) next # Skip if testing failed

  print("SwitchAnalyzeRlist summary after testing:")
  print(summary(switchList_tested))
  if (nrow(switchList_tested$isoformFeatures) == 0 || is.null(switchList_tested$switchAnalyzeRlist)) {
     warning(paste("No significant switches found or error during testing for", current_tissue, "- Skipping sequence analysis."))
     # Still save the tested list even if no switches are significant
     saveRDS(switchList_tested, file = file.path(output_dir, paste0(current_tissue, "_switchList_tested_noSig.rds")))
     next
  }

  # --- 4.5 Sequence Analysis (Part 1: ORF Prediction) ---
  if (run_orf_analysis) {
      print(paste("Running ORF analysis for", current_tissue, "..."))
      switchList_analyzed <- tryCatch({
          analyzeORF(
              switchAnalyzeRlist = switchList_tested,
              genomeObject = NULL, # Already provided via fasta_path in importRdata
              orfMethod = "longest", # or 'longest.AnnotatedWhenPossible' or 'mostUpstream'
              showProgress = TRUE,
              quiet = FALSE
          )
      }, error = function(e) {
          message(paste("Error during ORF analysis for", current_tissue, ":", e$message))
          switchList_tested # Return the state before the failed analysis
      })
  } else {
      switchList_analyzed <- switchList_tested # Skip ORF analysis
  }


  # --- 4.6 Sequence Analysis (Part 2: External Tools - Conditional) ---
  # Only run if the corresponding external tools are configured and flags are TRUE

  analysis_args <- list(
      switchAnalyzeRlist = switchList_analyzed, # Start with potentially ORF-analyzed list
      pathToOutput = file.path(output_dir, paste0(current_tissue, "_external_analsysis")), # Subdir for tool outputs
      nCores = analysis_threads,
      quiet = FALSE
  )

  if (!dir.exists(analysis_args$pathToOutput)) {
      dir.create(analysis_args$pathToOutput, recursive = TRUE)
  }

  if (run_cpat_analysis && file.exists(cpat_path)) {
      print(paste("Running CPAT analysis for", current_tissue, "..."))
      analysis_args$pathToCPAT <- cpat_path
      switchList_analyzed <- tryCatch(analyzeCPAT(analysis_args), error = function(e) {
          message(paste("CPAT analysis failed:", e$message)); switchList_analyzed })
  }
  if (run_signalp_analysis && file.exists(signalp_path)) {
      print(paste("Running SignalP analysis for", current_tissue, "..."))
      analysis_args$pathToSignalP <- signalp_path
      analysis_args$signalPversion <- 5 # Adjust if using a different version
      switchList_analyzed <- tryCatch(analyzeSignalP(analysis_args), error = function(e) {
          message(paste("SignalP analysis failed:", e$message)); switchList_analyzed })
  }
  if (run_pfam_analysis && file.exists(pfam_path)) {
      print(paste("Running Pfam analysis for", current_tissue, "..."))
      analysis_args$pathToPfam <- pfam_path
      analysis_args$pfamDb <- "path/to/pfam/db" # USER MUST PROVIDE PATH TO PFAM DB DIRECTORY
      switchList_analyzed <- tryCatch(analyzePFAM(analysis_args), error = function(e) {
          message(paste("Pfam analysis failed:", e$message)); switchList_analyzed })
  }
    if (run_netsurfp2_analysis && file.exists(netsurfp2_path)) {
      print(paste("Running NetSurfP-2 analysis for", current_tissue, "..."))
      analysis_args$pathToNetSurfP2 <- netsurfp2_path
      switchList_analyzed <- tryCatch(analyzeNetSurfP2(analysis_args), error = function(e) {
          message(paste("NetSurfP-2 analysis failed:", e$message)); switchList_analyzed })
  }

  # Final list object after potential sequence analyses
  switchList_final <- switchList_analyzed


  # --- 4.7 Consequence Analysis ---
  print(paste("Analyzing switch consequences for", current_tissue, "..."))
  # Define features for consequence analysis (can be customized)
  features_to_analyze = c(
      'isoform_length',
      'exon_number',
      'intron_structure',
      'domains_identified', # Requires Pfam analysis
      'ORF_genomic',
      'ORF_length',
      'ORF_contains_start_codon',
      'ORF_contains_stop_codon',
      'NMD_status',          # Requires ORF analysis
      'tss',
      'last_exon_length',
      'signal_peptide_identified', # Requires SignalP analysis
      'coding_potential'     # Requires CPAT analysis
  )

  # Filter features based on which analyses were actually run
  if (!run_pfam_analysis) features_to_analyze <- setdiff(features_to_analyze, 'domains_identified')
  if (!run_orf_analysis) features_to_analyze <- setdiff(features_to_analyze, c('NMD_status', 'ORF_genomic', 'ORF_length', 'ORF_contains_start_codon', 'ORF_contains_stop_codon'))
  if (!run_signalp_analysis) features_to_analyze <- setdiff(features_to_analyze, 'signal_peptide_identified')
  if (!run_cpat_analysis) features_to_analyze <- setdiff(features_to_analyze, 'coding_potential')

  switchList_consequences <- tryCatch({
    analyzeSwitchConsequences(
      switchAnalyzeRlist = switchList_final,
      consequencesToAnalyze = features_to_analyze,
      alpha = alpha_value,
      dIFcutoff = dIF_cutoff,
      showProgress = TRUE,
      quiet = FALSE
    )
  }, error = function(e) {
    message(paste("Error analyzing switch consequences for", current_tissue, ":", e$message))
    switchList_final # Return the list state before the failed consequence analysis
  })


  # --- 4.8 Save Results ---
  print(paste("Saving results for", current_tissue, "..."))

  # Save the final switchAnalyzeRlist object
  final_rds_path <- file.path(output_dir, paste0(current_tissue, "_switchAnalyzeRlist_final.rds"))
  saveRDS(switchList_consequences, file = final_rds_path)
  print(paste("Saved final switchAnalyzeRlist to:", final_rds_path))

  # Extract and save key results tables (optional)
  if (!is.null(switchList_consequences$switchAnalysis)) {
    switch_results_table <- switchList_consequences$switchAnalysis
    write.csv(switch_results_table, file = file.path(output_dir, paste0(current_tissue, "_significant_switches.csv")), row.names = FALSE)
    print(paste("Saved significant switches table for", current_tissue))
  } else {
     print("No switch analysis results found to save.")
  }

  if (!is.null(switchList_consequences$switchConsequence)) {
    consequence_table <- switchList_consequences$switchConsequence
    write.csv(consequence_table, file = file.path(output_dir, paste0(current_tissue, "_switch_consequences.csv")), row.names = FALSE)
    print(paste("Saved switch consequences table for", current_tissue))
  } else {
     print("No switch consequence results found to save.")
  }


  # Store the final object in the overall list (useful if loaded later)
  all_results_list[[current_tissue]] <- switchList_consequences

  print(paste(">>> Finished analysis for:", current_tissue))

} # End of tissue loop

# --- 5. Combine Results (Example - if needed later) ---
# The script `fullscript.txt` might load these RDS files individually.
# Or you could combine them here if needed, e.g., for cross-tissue comparisons.
# Example:
# combinedForebrain <- all_results_list[["Forebrain"]]
# combinedMidbrain  <- all_results_list[["Midbrain"]]
# combinedHindbrain <- all_results_list[["Hindbrain"]]
# print("Combined objects created (example).")

print("--- isoformSwitchAnalyzeR Analysis Complete ---") 