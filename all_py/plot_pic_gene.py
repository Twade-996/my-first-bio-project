#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import re
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D, proj3d

# --- Helper Functions (self-contained, same as before) ---
def ReadPDB(filename):
    x,y,z=[],[],[]; f=open(filename,'r'); [ (x.append(float(l[30:38])),y.append(float(l[38:46])),z.append(float(l[46:54]))) for l in f if l.startswith("ATOM")]; f.close(); return x,y,z,None
def ReadCompartmentScores(filename):
    s=[]; f=open(filename,'r'); [s.append(float(l.strip().split()[1])) for l in f if len(l.strip().split())>=2]; f.close(); return s
def parse_bin_gene_map(map_file):
    bin_to_gene, gene_to_bin = {}, {}
    print(f"Parsing bin-to-gene map from '{map_file}'...")
    try:
        with open(map_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 13 or parts[6] != 'gene': continue
                bin_id_str, attributes_str = parts[3], parts[12]
                try: bin_index = int(bin_id_str.split('_')[1]) - 1
                except (IndexError, ValueError): continue
                match = re.search(r'gene_name "([^"]+)"', attributes_str)
                if match:
                    gene_name = match.group(1)
                    if bin_index not in bin_to_gene: bin_to_gene[bin_index] = set()
                    bin_to_gene[bin_index].add(gene_name)
                    if gene_name not in gene_to_bin: gene_to_bin[gene_name] = bin_index
    except FileNotFoundError: print(f"Error: Bin-gene map file '{map_file}' not found."); sys.exit(1)
    for bin_index in bin_to_gene: bin_to_gene[bin_index] = list(bin_to_gene[bin_index])
    print(f"Done. Found {len(gene_to_bin)} unique genes in {len(bin_to_gene)} bins.")
    return bin_to_gene, gene_to_bin
def center_coords(coords): return coords - coords.mean(axis=0)

def main():
    parser = argparse.ArgumentParser(description="Generate a publication-quality 3D figure with automated gene expression hotspot highlighting.")
    parser.add_argument("pdb"); parser.add_argument("scores1"); parser.add_argument("scores2"); parser.add_argument('--bin-gene-map', required=True)
    parser.add_argument('--de-results', required=True, help="Path to the FULL, GENOME-WIDE differential expression results (.csv).")
    parser.add_argument("-o", '--output', required=True); parser.add_argument("-lw", "--linewidth", type=float, default=1.5)
    parser.add_argument("--skip-range", nargs=2, type=int, metavar=('S', 'E')); parser.add_argument("--focus-range", nargs=2, type=int, metavar=('S', 'E'))
    parser.add_argument("-t", "--threshold", type=float, default=0.1); parser.add_argument('--show-bins', action='store_true'); parser.add_argument('--bin-size', type=float, default=15.0)
    parser.add_argument('--top-genes', type=int, metavar='N', help="Highlight the bins containing the top N most up-regulated genes on this chromosome.")
    parser.add_argument('--view', nargs=2, type=float, default=[30, -60], metavar=('ELEV', 'AZIM'))
    args = parser.parse_args()

    # --- Data Loading and Processing ---
    print("--- Step 1 & 2: Loading and Processing Data ---")
    bin_to_gene, gene_to_bin = parse_bin_gene_map(args.bin_gene_map)
    coords = np.array(ReadPDB(args.pdb)[:3]).T; scores1_raw = ReadCompartmentScores(args.scores1); scores2_raw = ReadCompartmentScores(args.scores2)
    if args.skip_range: start, end = args.skip_range[0]-1, args.skip_range[1]-1; scores1, scores2 = scores1_raw[:start]+scores1_raw[end+1:], scores2_raw[:start]+scores2_raw[end+1:]
    else: scores1, scores2 = scores1_raw, scores2_raw
    min_len = min(len(coords), len(scores1), len(scores2)); coords, scores1, scores2 = coords[:min_len], np.array(scores1[:min_len]), np.array(scores2[:min_len])
    coords_centered = center_coords(coords)

    # --- Step 3: Categorize Bins and Prepare Colors ---
    print("--- Step 3: Determining Colors ---")
    bin_categories = []
    for s1, s2 in zip(scores1, scores2):
        state1='A' if s1>args.threshold else ('B' if s1<-args.threshold else 'T'); state2='A' if s2>args.threshold else ('B' if s2<-args.threshold else 'T')
        if state1=='B' and state2=='A': bin_categories.append('B -> A')
        elif state1=='A' and state2=='B': bin_categories.append('A -> B')
        else: bin_categories.append('Stable/Other')
    
    background_color = 'plum'
    # CORRECTED: Create the numpy array with dtype=object to prevent string truncation.
    line_colors = np.array(['red' if c == 'B -> A' else 'blue' if c == 'A -> B' else background_color for c in bin_categories], dtype=object)
    bin_colors = np.copy(line_colors)

    # --- Step 4: Find Top Genes and Override Bin Colors ---
    if args.top_genes:
        print("--- Step 4: Finding Top Genes to Highlight ---")
        de_df = pd.read_csv(args.de_results)
        genes_on_this_chromosome = list(gene_to_bin.keys())
        de_df_chr = de_df[de_df['gene_name'].isin(genes_on_this_chromosome)]
        if not de_df_chr.empty:
            top_genes_df = de_df_chr.sort_values(by='log2FoldChange', ascending=False).head(args.top_genes)
            print(f"  - Highlighting bins for the top {len(top_genes_df)} up-regulated genes on this chromosome:")
            for _, row in top_genes_df.iterrows():
                gene, log2fc = row['gene_name'], row['log2FoldChange']
                if gene in gene_to_bin:
                    bin_idx = gene_to_bin[gene]
                    if bin_idx < len(bin_colors):
                        bin_colors[bin_idx] = 'yellow'
                        print(f"    - Bin {bin_idx+1} (Category: {bin_categories[bin_idx]}) -> Gene: {gene} (log2FC: {log2fc:.2f}). Will be colored YELLOW.")
                    else:
                        print(f"    - (Skipping) Top Gene {gene} in Bin {bin_idx+1} is outside the modeled structure range.")
        else:
            print("  - Warning: No matching genes found to highlight.")
            
    # --- Step 5, 6, 7: Plotting and Saving ---
    focus_start_idx, focus_end_idx = 0, len(coords)
    if args.focus_range: focus_start_idx, focus_end_idx = args.focus_range[0]-1, args.focus_range[1]
    coords_to_plot, line_colors_to_plot, bin_colors_to_plot = coords_centered[focus_start_idx:focus_end_idx], line_colors[focus_start_idx:focus_end_idx], bin_colors[focus_start_idx:focus_end_idx]
    print("\n--- Step 5: Generating Plot ---")
    fig = plt.figure(figsize=(12, 12)); ax = fig.add_subplot(111, projection='3d')
    if len(coords_to_plot) > 1:
        start_idx, current_color = 0, line_colors_to_plot[0]
        for i in range(1, len(coords_to_plot)):
            if line_colors_to_plot[i] != current_color or i == len(coords_to_plot) - 1:
                segment = coords_to_plot[start_idx:i+1]; ax.plot(segment[:,0],segment[:,1],segment[:,2],color=current_color,linewidth=args.linewidth, zorder=1)
                start_idx, current_color = i, line_colors_to_plot[i]
    if args.show_bins:
        ax.scatter(coords_to_plot[:,0],coords_to_plot[:,1],coords_to_plot[:,2], c=bin_colors_to_plot, s=args.bin_size, depthshade=True, alpha=1.0, zorder=2)
    print("--- Step 6: Finalizing and Saving ---")
    ax.set_box_aspect([1,1,1]); ax.view_init(elev=args.view[0], azim=args.view[1]); ax.axis('off')
    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    print(f"\nAnalysis complete! Plot saved to '{args.output}'")

if __name__ == "__main__":
    main()