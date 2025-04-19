# --- thesis_scripts/Chapter3/01_setup_load_data.R ---
# Purpose: Load required packages and initial data for Chapter 3 analyses.
#          Reads nORF GTF, canonical ORF FASTA, and generates/loads random sequences.
#          Extracts sequence attributes.
# Chapter Relevance: 3.3, 3.4, 3.5, 3.6 (Data Setup)

## Package Installation (Commented out - Use main README instructions / renv)
# # Install and Load the BiocManager package
# # install.packages("BiocManager")
# # library(BiocManager)
#
# # install.packages("protr")
# # install.packages('svglite')
# # BiocManager::install("kebabs")
# # install.packages("stringi")
# # BiocManager::install("AnnotationHub")
# # BiocManager::install("protr")
# # BiocManager::install("ballgown") # ballgown seems related to RNA-seq, maybe not needed here?
# # install.packages("ggplot2")
# # install.packages("corrplot")
# # install.packages("ggExtra")
# # install.packages("plyr")
# # install.packages("ggpubr")
# # install.packages("caret")
# # install.packages("gridExtra")
# # install.packages("pROC")
# # install.packages("xgboost")
# # install.packages("data.table")
# # install.packages("drat", repos="https://cran.rstudio.com")
# # drat:::addRepo("dmlc")
# # install.packages("xgboost") # Redundant install
# # install.packages("phylotools") # Used later for PDB reading?
# # BiocManager::install("phastCons100way.UCSC.hg38") # Used later?

## Package initialisation
library(kebabs)
library(stringi)
library(protr)
# library(ballgown) # Confirm if needed for Ch3
library(ggplot2)
library(corrplot)
# library(ggExtra) # Needed later for plotting
# library(plyr) # Needed later?
# library(ggpubr) # Needed later for plotting
# library(caret) # Needed later for ML
# library(gridExtra) # Needed later for plotting
# library(pROC) # Needed later for ML
# library(xgboost) # Needed later for ML
library(data.table)
library(Biostrings)
# library(phylotools) # Needed later?
# library(phastCons100way.UCSC.hg38) # Needed later?


## Data Loading

# Define input file paths (Modify these paths as needed)
norf_gtf_file <- "nORFsDB.1.1.gtf"
canonical_fasta_file <- "UP000005640_9606_20k.fasta"
random_seq_file <- "1000RandomSequences.csv" # Or set to NULL to generate random seqs
fasta_output_dir <- "."
missing_seq_output_file <- "nORFs_missing.fasta"
pdb_dir <- './pdb'

# nORFs from nORFsDB 1.1 GTF
cat("Loading nORFs from:", norf_gtf_file, "\n")
if (!file.exists(norf_gtf_file)) stop("nORF GTF file not found: ", norf_gtf_file)
norfsDB <- tryCatch(ballgown::gffRead(norf_gtf_file), error = function(e) {
    stop("Error reading GTF file: ", conditionMessage(e))
})

cat("Extracting nORF attributes...\n")
norfsAttributesAAseq <- ballgown::getAttributeField(norfsDB$attributes, "AA_seq", attrsep = "; ")
norfsAttributesGeneId <- ballgown::getAttributeField(norfsDB$attributes, "gene_id", attrsep = "; ")
norfsAttributesStart <- ballgown::getAttributeField(norfsDB$attributes, "start_codon", attrsep = "; ")

# Clean extracted attributes
nORFAAseq <- gsub("[^[:alnum:]=\\.]", "", norfsAttributesAAseq)
nORFGeneIds <- gsub("[^[:alnum:]=\\.]", "", norfsAttributesGeneId)
nORFStartCodons <- gsub("[^[:alnum:]=\\.]", "", norfsAttributesStart)

cat("Loaded", length(nORFGeneIds), "nORFs.\n")

# --- Optional: Check PDB file existence and write missing sequences ---
# This section assumes PDB files named by gene_id exist in pdb_dir
if (dir.exists(pdb_dir)) {
    cat("Checking for missing PDB files in:", pdb_dir, "\n")
    pdb_files_present <- list.files(path = pdb_dir, pattern = "\\.pdb$", all.files = TRUE, recursive = TRUE)
    pdb_ids_present <- sub(".pdb", "", basename(pdb_files_present))
    pdb_ids_missing <- setdiff(nORFGeneIds, pdb_ids_present)
    cat("Found", length(pdb_ids_missing), "nORFs with missing PDB files.\n")

    if (length(pdb_ids_missing) > 0) {
        missing_indices <- which(nORFGeneIds %in% pdb_ids_missing)
        nORFGeneIds_missing <- nORFGeneIds[missing_indices]
        nORFAAseq_missing <- nORFAAseq[missing_indices]

        cat("Writing missing sequences to:", missing_seq_output_file, "\n")
        write.table(paste0(">", nORFGeneIds_missing, "\n", nORFAAseq_missing),
                    file = file.path(fasta_output_dir, missing_seq_output_file),
                    sep = "", row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
} else {
    cat("PDB directory not found:", pdb_dir, "- Skipping PDB check.\n")
}

# --- Create .fasta file for all nORFs ---
norf_fasta_output_file <- "nORFs_all.fasta"
cat("Writing all nORF sequences to:", norf_fasta_output_file, "\n")
write.table(paste0(">", nORFGeneIds, "\n", nORFAAseq),
            file = file.path(fasta_output_dir, norf_fasta_output_file),
            sep = "", row.names = FALSE, col.names = FALSE, quote = FALSE)


# --- Load/Generate Random Sequences ---
if (!is.null(random_seq_file) && file.exists(random_seq_file)) {
    cat("Loading random sequences from:", random_seq_file, "\n")
    randomSeqDF <- read.csv(random_seq_file, header = TRUE)
    randomSequencesAA <- randomSeqDF$aa # Assuming column name is 'aa'
    cat("Loaded", length(randomSequencesAA), "random sequences.\n")
} else {
    cat("Generating 1000 random sequences (50-1000 AA length)...\n")
    set.seed(11) # for reproducibility
    seqLength <- runif(1000, min = 50, max = 1000)
    randomSequences <- kebabs::genRandBioSeqs(seqType = "AA", numSeqs = 1000, seqLength = seqLength,
                                             biostring = TRUE, seed = 11)
    randomSequencesAA <- as.character(randomSequences)
    # Optional: Save generated sequences
    # write.csv(data.frame(aa = randomSequencesAA), "generated_random_sequences.csv", row.names = FALSE)
    cat("Generated", length(randomSequencesAA), "random sequences.\n")
}

# --- Load Canonical ORFs ---
cat("Loading canonical ORFs from:", canonical_fasta_file, "\n")
if (!file.exists(canonical_fasta_file)) stop("Canonical ORF FASTA file not found: ", canonical_fasta_file)
cORFfastaFile <- tryCatch(Biostrings::readDNAStringSet(canonical_fasta_file), error = function(e) {
    # Try reading as AAStringSet if DNAStringSet fails
    tryCatch(Biostrings::readAAStringSet(canonical_fasta_file), error = function(e2) {
        stop("Error reading canonical FASTA (tried DNA and AA): ", conditionMessage(e2))
    })
})

cORFSeqNames <- names(cORFfastaFile)
cORFSequencesAA <- as.character(cORFfastaFile)
cat("Loaded", length(cORFSeqNames), "canonical ORFs.\n")


# --- Plot Start Codon Usage (Optional Visualization) ---
cat("Plotting nORF start codon usage...\n")
codon_counts <- table(nORFStartCodons)
codon_data <- data.frame(codon = names(codon_counts), count = as.numeric(codon_counts))

p_codon_usage <- ggplot(codon_data, aes(x = reorder(codon, -count), y = count)) +
    geom_col() +
    labs(title = "nORF Start Codon Usage", x = "Start Codon", y = "Count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

print(p_codon_usage)
ggsave("nORF_start_codon_usage.png", plot = p_codon_usage, width = 6, height = 4)

# --- Save Loaded Data (Optional) ---
# It might be useful to save the key loaded objects for faster reloading in subsequent scripts
save_data_file <- "chapter3_loaded_data.RData"
cat("Saving loaded sequences and IDs to:", save_data_file, "\n")
save(nORFAAseq, nORFGeneIds, nORFStartCodons,
     cORFSequencesAA, cORFSeqNames,
     randomSequencesAA,
     file = save_data_file)

cat("--- 01_setup_load_data.R finished ---\n") 