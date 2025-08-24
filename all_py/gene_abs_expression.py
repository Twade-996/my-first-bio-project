#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import re
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, wilcoxon
from functools import reduce

# --- Helper Functions (Self-contained) ---

def ReadCompartmentScores(filename):
    """Reads a two-column file (bin_index, score) and returns a list of scores."""
    scores = []
    try:
        with open(filename, 'r') as score_file:
            for line in score_file:
                parts = line.strip().split()
                if len(parts) >= 2:
                    try: scores.append(float(parts[1]))
                    except ValueError: pass
    except FileNotFoundError: print(f"Error: Score file '{filename}' not found."); sys.exit(1)
    return scores

def create_id_map_from_gtf(gtf_file_path):
    """Parses a GTF file to create a map from ENSEMBL ID (version stripped) to gene name."""
    print(f"Building ENSEMBL ID to Gene Name map from '{gtf_file_path}'...")
    id_map = {}
    try:
        with open(gtf_file_path, 'r') as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')
                if len(parts) < 9 or parts[2] != 'gene': continue
                attributes_str = parts[8]
                gene_id_match = re.search(r'gene_id "([^"]+)"', attributes_str)
                gene_name_match = re.search(r'gene_name "([^"]+)"', attributes_str)
                if gene_id_match and gene_name_match:
                    ensembl_id_no_version = gene_id_match.group(1).split('.')[0]
                    id_map[ensembl_id_no_version] = gene_name_match.group(1)
    except FileNotFoundError: print(f"Error: GTF file '{gtf_file_path}' not found."); sys.exit(1)
    print(f"Map created successfully with {len(id_map)} entries.")
    return id_map

def process_rsem_to_tpm_map(rsem_files, id_to_name_map):
    """Reads RSEM files, converts IDs, merges, and returns a gene_name -> mean_TPM dictionary."""
    # ... This function remains correct and unchanged ...
    all_dataframes = []
    print("\n--- Step 1b: Processing RSEM quantification files ---")
    for i, file_path in enumerate(rsem_files):
        print(f"  - Reading file: {file_path}")
        try:
            df = pd.read_csv(file_path, sep='\t')
            df['gene_id_no_version'] = df['gene_id'].str.split('.').str[0]
            df['gene_name'] = df['gene_id_no_version'].map(id_to_name_map)
            df_processed = df[['gene_name', 'TPM']].copy().rename(columns={'TPM': f'TPM_rep{i+1}'})
            df_processed.dropna(subset=['gene_name'], inplace=True)
            all_dataframes.append(df_processed)
        except Exception as e:
            print(f"  Warning: Could not process file '{file_path}'. Error: {e}. Skipping.")
    if not all_dataframes: print("Error: No valid RSEM files processed."); sys.exit(1)
    print("  - Merging replicates and calculating mean TPM...")
    df_merged = reduce(lambda left, right: pd.merge(left, right, on='gene_name', how='outer'), all_dataframes)
    df_merged.set_index('gene_name', inplace=True)
    mean_tpm_series = df_merged.fillna(0).mean(axis=1)
    print("  - Mean TPM calculation complete.")
    return mean_tpm_series.to_dict()

def parse_bin_gene_map(map_file):
    """Parses bedtools output to create a bin -> [genes] map. Keys are 0-based integers."""
    bin_to_gene = {}
    print(f"\n--- Step 2: Parsing bin-to-gene map from '{map_file}' ---")
    try:
        with open(map_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 13 or parts[6] != 'gene': continue
                # CORRECTED LOGIC: Convert bin_ID string to 0-based integer index
                bin_id_str = parts[3] # e.g., 'bin_1'
                try:
                    bin_index = int(bin_id_str.split('_')[1]) - 1 # 'bin_1' -> 0
                except (IndexError, ValueError):
                    continue # Skip if the format is not as expected
                
                attributes_str = parts[12]
                match = re.search(r'gene_name "([^"]+)"', attributes_str)
                if match:
                    gene_name = match.group(1)
                    if bin_index not in bin_to_gene: bin_to_gene[bin_index] = set()
                    bin_to_gene[bin_index].add(gene_name)
    except FileNotFoundError: print(f"Error: Bin-gene map file '{map_file}' not found."); sys.exit(1)
    for bin_index in bin_to_gene:
        bin_to_gene[bin_index] = list(bin_to_gene[bin_index])
    print(f"Done. Found genes in {len(bin_to_gene)} bins.")
    return bin_to_gene

def main():
    parser = argparse.ArgumentParser(description="A comprehensive script to analyze absolute expression levels vs. compartment switches.")
    # ... (all arguments remain the same) ...
    parser.add_argument('--rsem-files', nargs='+', required=True, help="List of RSEM .tsv files for the target cell line.")
    parser.add_argument('--gtf', required=True, help="Path to the GENCODE annotation GTF file for ID conversion.")
    parser.add_argument('--bin-gene-map', required=True, help="Path to the bin-to-gene map file from bedtools.")
    parser.add_argument('--scores1', required=True, help="Path to c-score file for the starting condition.")
    parser.add_argument('--scores2', required=True, help="Path to c-score file for the target condition.")
    parser.add_argument('--target-name', default="TargetCell", help="Name of the target cell line.")
    parser.add_argument('--skip-range', nargs=2, type=int, metavar=('S', 'E'), help="1-based bin range to skip.")
    parser.add_argument('-t', '--threshold', type=float, default=0.1, help="Absolute c-score threshold.")
    parser.add_argument('-o', '--output', default="absolute_expression_levels.png", help="Output plot file name.")
    args = parser.parse_args()

    # --- Step 1 & 2: Load all raw data ---
    id_map = create_id_map_from_gtf(args.gtf)
    tpm_map = process_rsem_to_tpm_map(args.rsem_files, id_map)
    bin_to_gene = parse_bin_gene_map(args.bin_gene_map)
    scores1 = ReadCompartmentScores(args.scores1)
    scores2 = ReadCompartmentScores(args.scores2)

    # --- Step 3: Core logic with corrected indexing and skipping ---
    print("\n--- Step 3: Associating expression with compartments ---")
    results = {'B -> A': [], 'Stable A': []}
    
    # The total number of bins is determined by the length of the shortest score file
    num_total_bins = min(len(scores1), len(scores2))
    print(f"Analyzing {num_total_bins} total bins based on score files.")
    
    # Create a set of bins to skip for fast lookup
    skip_indices = set()
    if args.skip_range:
        # Convert 1-based skip range to 0-based indices
        start_idx, end_idx = args.skip_range[0] - 1, args.skip_range[1] - 1
        skip_indices = set(range(start_idx, end_idx + 1))
        print(f"Will skip {len(skip_indices)} bins in range {args.skip_range}.")

    # Iterate through every bin using a universal 0-based index
    for i in range(num_total_bins):
        # CORRECTED LOGIC: Check if the current bin should be skipped
        if i in skip_indices:
            continue
            
        # Get scores for the current valid bin
        s1 = scores1[i]
        s2 = scores2[i]
        
        # Categorize the bin
        state1 = 'A' if s1 > args.threshold else ('B' if s1 < -args.threshold else 'T')
        state2 = 'A' if s2 > args.threshold else ('B' if s2 < -args.threshold else 'T')
        category = None
        if state1 == 'B' and state2 == 'A': category = 'B -> A'
        elif state1 == 'A' and state2 == 'A': category = 'Stable A'
        
        # If the category is one we care about, and this bin contains genes...
        if category and i in bin_to_gene:
            for gene in bin_to_gene[i]:
                if gene in tpm_map:
                    results[category].append(np.log2(tpm_map[gene] + 1))

    # --- Step 4: Plotting and Statistics (unchanged) ---
    print("\n--- Step 4: Generating final plot and statistics ---")
    plot_data = [{'Category': cat, 'log2(TPM+1)': val} for cat, val_list in results.items() for val in val_list]
    plot_df = pd.DataFrame(plot_data)
    
    if plot_df.empty:
        print("Error: No data to plot. This might happen if no genes were found in the specified categories."); sys.exit(1)

    plt.figure(figsize=(6, 8))
    sns.boxplot(data=plot_df, x='Category', y='log2(TPM+1)', order=['Stable A', 'B -> A'],
                palette={'Stable A': 'lightcoral', 'B -> A': 'red'})
    plt.title(f'Absolute Expression Levels in {args.target_name}', fontsize=16)
    plt.ylabel(f'Expression Level (log2(TPM+1)) in {args.target_name}', fontsize=12)
    plt.xlabel('Compartment Category', fontsize=12)
    
    try:
        u_stat, p_value = mannwhitneyu(results['B -> A'], results['Stable A'], alternative='two-sided')
        stats_text = f"Mann-Whitney U Test:\np-value = {p_value:.2e}"
        print(f"\nStatistical Test Results:\n{stats_text}")
        plt.text(0.5, 0.95, stats_text, ha='center', va='center', transform=plt.gca().transAxes, fontsize=12,
                 bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
    except ValueError as e:
        print(f"\nCould not perform statistical test: {e}")

    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    print(f"\nAnalysis complete! Plot saved to '{args.output}'")
    
if __name__ == "__main__":
    main()