"""
AlphaGenome nORF Functional Analysis Script

This script performs functional importance mapping of non-coding ORFs (nORFs) using
the AlphaGenome deep learning model. It processes GTF annotation files, filters
nORFs by length, and generates functional importance maps showing how genomic
variations affect RNA-seq expression.

The script implements a chunked, resumable analysis approach to handle large
datasets efficiently and provides checkpoint functionality to resume interrupted
analyses.

Author: Robin Gounder
Date: 2025-07-01
Version: 1.0

Dependencies:
    - alphagenome: Deep learning model for genomic analysis
    - pandas: Data manipulation and analysis
    - numpy: Numerical computing
    - matplotlib: Plotting and visualization
    - tqdm: Progress bars
    - IPython: Jupyter notebook integration

Usage:
    This script is designed to run in a Jupyter notebook environment with
    the required GTF file uploaded to the working directory.
"""

# @title Step 1: Install Libraries and Import
from IPython.display import clear_output
!pip install alphagenome -q
!pip install tqdm -q
clear_output()

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import os
import time
import json
from tqdm.auto import tqdm

# AlphaGenome specific imports
from alphagenome.data import gene_annotation, genome
from alphagenome.models import dna_client
from alphagenome.models import variant_scorers
from alphagenome.interpretation import ism

# @title Step 2: Load AlphaGenome Model
try:
    api_key = 'PLACEHOLDER API KEY' # IMPORTANT: Get AlphaGenome API key from https://alphagenome.com/
    dna_model = dna_client.create(api_key)
    print("✅ AlphaGenome model loaded successfully.")
except Exception as e:
    print(f"❌ Error loading AlphaGenome model: {e}")
    raise

# @title Step 3: Load nORFs from GTF File and Filter by Length
GTF_FILENAME = 'nORFsDB.1.1.gtf' # <-- CHANGE THIS IF YOU HAVE A DIFFERENT GTF FILE
MAX_NORF_LENGTH = 1000 # <-- SET a maximum length for analysis

def parse_gtf_attributes(attribute_string):
    attributes = {}
    for key, value in re.findall(r'(\w+)\s"([^"]*)"', attribute_string):
        attributes[key] = value
    return attributes
try:
    print(f"Loading and parsing {GTF_FILENAME}...")
    gtf_cols = ['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
    df_norfs = pd.read_csv(GTF_FILENAME, sep='\t', comment='#', header=None, names=gtf_cols, low_memory=False)
    all_norf_data = []
    for index, row in df_norfs.iterrows():
        attrs = parse_gtf_attributes(row['attributes'])
        chrom = str(row['chromosome'])
        if not chrom.startswith('chr'): chrom = 'chr' + chrom
        if pd.isna(row['start']) or pd.isna(row['end']): continue
        norf_entry = { 'norf_id': attrs.get('gene_id'), 'chromosome': chrom, 'start': int(row['start']) - 1, 'end': int(row['end']) }
        all_norf_data.append(norf_entry)
    print(f"Successfully loaded {len(all_norf_data)} total nORFs.")
    norf_data_filtered = [ n for n in all_norf_data if (n['end'] - n['start']) <= MAX_NORF_LENGTH ]
    print(f"✅ Kept {len(norf_data_filtered)} nORFs that are shorter than {MAX_NORF_LENGTH} bp for analysis.")
except FileNotFoundError:
    print(f"❌ ERROR: File not found! Please make sure '{GTF_FILENAME}' is uploaded.")
    norf_data_filtered = []

# @title Step 4: Run Resumable, Chunked Analysis for All Filtered nORFs

def chunk_interval(interval: genome.Interval, chunk_size: int):
    """Breaks a large interval into a list of smaller ones."""
    chunks = []
    for start in range(interval.start, interval.end, chunk_size):
        end = min(start + chunk_size, interval.end)
        if start < end: # Ensure we don't create empty intervals
            chunks.append(genome.Interval(interval.chromosome, start, end))
    return chunks

if not norf_data_filtered:
    print("Skipping analysis because no nORFs remained after filtering.")
else:
    OUTPUT_DIR = 'alphagenome_images'
    CHECKPOINT_DIR = 'alphagenome_checkpoints'
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    print(f"\nOutputs will be saved to '{OUTPUT_DIR}' and '{CHECKPOINT_DIR}'.")

    ONTOLOGY_TERM = 'EFO:0002067'
    CHUNK_SIZE = 100  # Process 100bp at a time
    THROTTLE_SECONDS = 1 # Pause 1 second between chunks
    SUPPORTED_LENGTHS = sorted(dna_client.SUPPORTED_SEQUENCE_LENGTHS.values())

    for norf in tqdm(norf_data_filtered, desc="Total nORF Progress"):
        NORF_ID, NORF_CHROMOSOME, NORF_START, NORF_END = (
            norf['norf_id'], norf['chromosome'], norf['start'], norf['end']
        )
        safe_norf_id = re.sub(r'[/\\?%*:|"<>]', '_', NORF_ID)
        output_path = os.path.join(OUTPUT_DIR, f'{safe_norf_id}_functional_map.svg')
        checkpoint_path = os.path.join(CHECKPOINT_DIR, f'{safe_norf_id}_chunks_checkpoint.json')

        if os.path.exists(output_path):
            continue

        try:
            print(f"\n--- Analyzing {NORF_ID} ---")
            
            # 1. Define the full region and break it into chunks
            full_ism_interval = genome.Interval(NORF_CHROMOSOME, NORF_START - 512, NORF_END + 512)
            interval_chunks = chunk_interval(full_ism_interval, CHUNK_SIZE)

            # 2. Load checkpoint for completed chunks
            if os.path.exists(checkpoint_path):
                with open(checkpoint_path, 'r') as f:
                    try:
                        completed_chunks_data = json.load(f)
                    except json.JSONDecodeError: # Handle corrupted file
                        completed_chunks_data = {}
            else:
                completed_chunks_data = {}

            # 3. Process each chunk
            for chunk in tqdm(interval_chunks, desc=f"Scoring Chunks for {NORF_ID}", leave=False, ncols=100):
                chunk_key = str(chunk)
                if chunk_key in completed_chunks_data:
                    continue

                # Smartly size the API interval for this small chunk
                api_length = next((l for l in SUPPORTED_LENGTHS if l >= chunk.width), None)
                if api_length is None: continue
                api_interval = chunk.resize(api_length)

                rnaseq_scorer = variant_scorers.CenterMaskScorer(
                    requested_output=dna_client.OutputType.RNA_SEQ, width=501,
                    aggregation_type=variant_scorers.AggregationType.DIFF_MEAN,
                )
                
                # We can use the simple retry logic now because each job is small
                for attempt in range(3): # Simple 3-retry loop
                    try:
                        variant_scores = dna_model.score_ism_variants(
                            interval=api_interval, ism_interval=chunk, variant_scorers=[rnaseq_scorer]
                        )
                        # Process and store the successful result for this chunk
                        def extract_score(adata):
                            values = adata.X[:, adata.var['ontology_curie'] == ONTOLOGY_TERM]
                            return values.mean() if values.size > 0 else 0
                        
                        chunk_results = []
                        for score_obj in variant_scores:
                            variant = score_obj[0].uns['variant']
                            score = float(extract_score(score_obj[0]))
                            chunk_results.append({'pos': variant.position, 'score': score})
                        
                        completed_chunks_data[chunk_key] = chunk_results
                        with open(checkpoint_path, 'w') as f:
                            json.dump(completed_chunks_data, f)
                        
                        time.sleep(THROTTLE_SECONDS) # Pause after a successful chunk
                        break # Success
                    except Exception as e:
                        if attempt < 2 and 'quota' in str(e).lower():
                            time.sleep(60) # Longer pause for quota error
                        elif attempt == 2:
                             print(f"Warning: Failed to process chunk {chunk_key} after 3 attempts.")
                        else:
                            raise e # A different kind of error

            # 4. All chunks are processed, now assemble the full results and plot
            print(f"✅ Assembling and plotting results for {NORF_ID}...")
            position_scores = {}
            for chunk_key, results in completed_chunks_data.items():
                for item in results:
                    pos, score = item['pos'], item['score']
                    if pos not in position_scores: position_scores[pos] = []
                    position_scores[pos].append(abs(score))
            
            if not position_scores:
                print(f"--- No data to plot for {NORF_ID}. Skipping. ---")
                continue

            final_scores = {pos: np.mean(vals) for pos, vals in position_scores.items()}
            positions, scores = sorted(final_scores.keys()), [final_scores[p] for p in sorted(final_scores.keys())]

            fig, ax = plt.subplots(figsize=(18, 6))
            ax.plot(positions, scores, color='navy', alpha=0.8)
            ax.axvspan(NORF_START, NORF_END, color='palegreen', alpha=0.5, label=f'{NORF_ID} Body')
            ax.set_title(f'Functional Importance Map for {NORF_ID}', fontsize=16)
            ax.set_xlabel(f'Position on {NORF_CHROMOSOME}', fontsize=12)
            ax.set_ylabel('Functional Importance Score\n(Mean absolute effect on RNA-SEQ)', fontsize=12)
            ax.grid(True, linestyle='--', alpha=0.6)
            ax.legend()
            
            plt.savefig(output_path, format='svg', bbox_inches='tight')
            print(f"✅ Plot saved to '{output_path}'")
            plt.show()
            plt.close(fig)

        except Exception as e:
            print(f"❌ An unrecoverable error occurred during setup/analysis for {NORF_ID}: {e}")
            continue
    
    print("\n\n--- All analyses complete! ---")