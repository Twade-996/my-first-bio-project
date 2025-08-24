#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import re
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider
from functools import reduce

# --- Helper Functions (self-contained) ---
def ReadPDB(filename):
    x, y, z = [], [], []
    try:
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith("ATOM"):
                    x.append(float(line[30:38])); y.append(float(line[38:46])); z.append(float(line[46:54]))
    except FileNotFoundError: print(f"Error: PDB file '{filename}' not found."); sys.exit(1)
    except (ValueError, IndexError): print(f"Error parsing coordinates from '{filename}'."); sys.exit(1)
    return x, y, z, None

def ReadCompartmentScores(filename):
    scores = []
    try:
        with open(filename, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    try: scores.append(float(parts[1]))
                    except ValueError: pass
    except FileNotFoundError: print(f"Error: Score file '{filename}' not found."); sys.exit(1)
    return scores

class PanHandler:
    # ... [PanHandler class remains unchanged] ...
    def __init__(self, fig, ax): self.fig, self.ax, self.press, self.panning = fig, ax, None, False; self.connect()
    def connect(self):
        self.cid_press=self.fig.canvas.mpl_connect('button_press_event',self.on_press)
        self.cid_release=self.fig.canvas.mpl_connect('button_release_event',self.on_release)
        self.cid_motion=self.fig.canvas.mpl_connect('motion_notify_event',self.on_motion)
    def on_press(self, e):
        if e.inaxes != self.ax or e.button != 2: return
        self.panning, self.press = True, (e.x, e.y)
    def on_motion(self, e):
        if not self.panning or e.inaxes != self.ax: return
        dx, dy = (e.x-self.press[0])*-1, (e.y-self.press[1])*-1; xlim, ylim = self.ax.get_xlim(),self.ax.get_ylim()
        self.ax.set_xlim(xlim[0]+dx*(xlim[1]-xlim[0])*0.003, xlim[1]+dx*(xlim[1]-xlim[0])*0.003)
        self.ax.set_ylim(ylim[0]+dy*(ylim[1]-ylim[0])*0.003, ylim[1]+dy*(ylim[1]-ylim[0])*0.003)
        self.press = (e.x, e.y); self.fig.canvas.draw()
    def on_release(self, e):
        if e.button == 2: self.panning, self.press = False, None

def CalRotateM(S1, S2):
    A=np.dot(S1.T,S2); u,_,v=np.linalg.svd(A)
    if np.linalg.det(u)*np.linalg.det(v)<0: u[:,-1]*=-1
    return np.dot(u,v)
def center_coords(coords): return coords - coords.mean(axis=0)

# --- The Main Application Class ---
class Interactive3DViewer:
    # ... [init and most methods are the same] ...
    def __init__(self, args):
        self.args = args
        self.load_and_process_data()
        self.setup_plot()

    @staticmethod
    def create_id_map_from_gtf(gtf_file_path):
        # ... [Unchanged] ...
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

    @staticmethod
    def process_rsem_to_tpm_map(rsem_files, id_to_name_map):
        # ... [Unchanged] ...
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
            except Exception as e: print(f"  Warning: Could not process file '{file_path}'. Error: {e}. Skipping.")
        if not all_dataframes: print("Error: No valid RSEM files processed."); sys.exit(1)
        print("  - Merging replicates and calculating mean TPM...")
        df_merged = reduce(lambda left, right: pd.merge(left, right, on='gene_name', how='outer'), all_dataframes)
        df_merged.set_index('gene_name', inplace=True)
        mean_tpm_series = df_merged.fillna(0).mean(axis=1)
        print("  - Mean TPM calculation complete.")
        return mean_tpm_series.to_dict()

    @staticmethod
    def parse_bin_gene_map(map_file):
        # ... [Unchanged] ...
        bin_to_gene = {}
        print(f"\n--- Step 2: Parsing bin-to-gene map from '{map_file}' ---")
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
        except FileNotFoundError: print(f"Error: Bin-gene map file '{map_file}' not found."); sys.exit(1)
        for bin_index in bin_to_gene: bin_to_gene[bin_index] = list(bin_to_gene[bin_index])
        print(f"Done. Found genes in {len(bin_to_gene)} bins.")
        return bin_to_gene

    def load_and_process_data(self):
        # ... [Unchanged] ...
        print("--- Step 1: Loading and Processing Data ---")
        id_map = self.create_id_map_from_gtf(self.args.gtf)
        self.tpm_map = self.process_rsem_to_tpm_map(self.args.rsem_files, id_map) if self.args.rsem_files else {}
        self.de_map = self.load_de_results(self.args.rna_seq) if self.args.rna_seq else {}
        self.bin_to_gene = self.parse_bin_gene_map(self.args.bin_gene_map)
        scores1_raw, scores2_raw = ReadCompartmentScores(self.args.scores1), ReadCompartmentScores(self.args.scores2)
        if self.args.skip_range:
            start, end = self.args.skip_range[0]-1, self.args.skip_range[1]-1
            scores1, scores2 = scores1_raw[:start]+scores1_raw[end+1:], scores2_raw[:start]+scores2_raw[end+1:]
        else: scores1, scores2 = scores1_raw, scores2_raw
        coords1, coords2 = np.array(ReadPDB(self.args.pdb1)[:3]).T, np.array(ReadPDB(self.args.pdb2)[:3]).T
        min_len = min(len(coords1), len(scores1), len(coords2), len(scores2))
        self.scores1, self.scores2 = np.array(scores1[:min_len]), np.array(scores2[:min_len])
        self.coords1_centered = center_coords(coords1[:min_len])
        coords2_centered = center_coords(coords2[:min_len])
        self.coords2_aligned = np.dot(coords2_centered, CalRotateM(coords2_centered, self.coords1_centered))
        self.focus_start_idx, self.focus_end_idx = 0, min_len
        if self.args.focus_range: self.focus_start_idx, self.focus_end_idx = self.args.focus_range[0]-1, self.args.focus_range[1]
        self.bin_categories = self.categorize_bins()

    def setup_plot(self):
        # ... [Unchanged] ...
        self.fig = plt.figure(figsize=(10, 10))
        self.ax = self.fig.add_axes([0, 0.05, 1, 0.95], projection='3d')
        self.fig.canvas.mpl_connect('pick_event', self.on_pick)
        ax_slider = self.fig.add_axes([0.2, 0.02, 0.65, 0.03])
        self.morph_slider = Slider(ax=ax_slider, label='Morph', valmin=0.0, valmax=1.0, valinit=0.0)
        self.morph_slider.on_changed(self.update)
        self.pan_handler = PanHandler(self.fig, self.ax)
        self.update(0.0)
        plt.show()

    def update(self, val):
        self.ax.cla()
        t = val
        coords_interp = (1 - t) * self.coords1_centered + t * self.coords2_aligned
        coords_to_plot = coords_interp[self.focus_start_idx:self.focus_end_idx]
        categories_to_plot = self.bin_categories[self.focus_start_idx:self.focus_end_idx]
        
        # Determine colors (this logic is unchanged)
        if self.args.highlight_switches:
            colors = np.array(['red' if c=='B -> A' else 'blue' if c=='A -> B' else 'lightgray' for c in categories_to_plot])
        else:
            scores_to_plot = self.scores1[self.focus_start_idx:self.focus_end_idx] if t < 0.5 else self.scores2[self.focus_start_idx:self.focus_end_idx]
            colors = np.array(['red' if s > self.args.threshold else 'blue' if s < -self.args.threshold else 'gray' for s in scores_to_plot])

        # Plot lines (unchanged)
        if len(coords_to_plot) > 1:
            start_idx, current_color = 0, colors[0]
            for i in range(1, len(coords_to_plot)):
                if colors[i] != current_color or i == len(coords_to_plot) - 1:
                    segment = coords_to_plot[start_idx : i + 1]
                    self.ax.plot(segment[:,0], segment[:,1], segment[:,2], color=current_color, linewidth=self.args.linewidth)
                    start_idx, current_color = i, colors[i]
        
        # MODIFIED: Add visible spheres if requested
        if self.args.show_bins and len(coords_to_plot) > 0:
            self.ax.scatter(coords_to_plot[:,0], coords_to_plot[:,1], coords_to_plot[:,2], 
                            c=colors, s=self.args.bin_size, depthshade=True, alpha=0.8)

        # Plot invisible points for picking (unchanged)
        self.ax.scatter(coords_to_plot[:,0], coords_to_plot[:,1], coords_to_plot[:,2], s=20, alpha=0.0, picker=5)
        
        # Set view limits and title (unchanged)
        all_coords = np.vstack([self.coords1_centered, self.coords2_aligned])[self.focus_start_idx:self.focus_end_idx]
        if all_coords.size > 0:
            max_range = (all_coords.max(axis=0)-all_coords.min(axis=0)).max()*0.55; mid = all_coords.mean(axis=0)
            self.ax.set_xlim(mid[0]-max_range, mid[0]+max_range); self.ax.set_ylim(mid[1]-max_range, mid[1]+max_range); self.ax.set_zlim(mid[2]-max_range, mid[2]+max_range)
        self.ax.axis('off'); self.ax.set_title(f"{(1-t)*100:.0f}% {self.args.name1} | {t*100:.0f}% {self.args.name2}")
        self.fig.canvas.draw_idle()

    # ... [on_pick, categorize_bins, load_de_results methods are all unchanged] ...
    def on_pick(self, event):
        if not event.ind: return
        clicked_local_index = event.ind[0]; clicked_global_index = self.focus_start_idx + clicked_local_index
        bin_number = clicked_global_index + 1; category = self.bin_categories[clicked_global_index]
        genes = self.bin_to_gene.get(clicked_global_index, ["No genes found"])
        print("\n--- Bin Information ---")
        print(f" Bin Number: {bin_number} (Index: {clicked_global_index})"); print(f" Compartment Switch: {category}"); print(" Genes in Bin:")
        for gene in genes:
            log2fc = f"{self.de_map.get(gene, 'N/A'):.2f}" if isinstance(self.de_map.get(gene), float) else "N/A"
            tpm = f"{self.tpm_map.get(gene, 'N/A'):.2f}" if isinstance(self.tpm_map.get(gene), float) else "N/A"
            print(f"  - {gene:<15} | log2FC: {log2fc:<7} | TPM in {self.args.name2}: {tpm}")
    def categorize_bins(self):
        categories = []
        for s1, s2 in zip(self.scores1, self.scores2):
            state1='A' if s1>self.args.threshold else ('B' if s1<-self.args.threshold else 'T'); state2='A' if s2>self.args.threshold else ('B' if s2<-self.args.threshold else 'T')
            if state1=='B' and state2=='A': categories.append('B -> A')
            elif state1=='A' and state2=='B': categories.append('A -> B')
            elif state1=='A' and state2=='A': categories.append('Stable A')
            elif state1=='B' and state2=='B': categories.append('Stable B')
            else: categories.append('Other/Transition')
        return categories
    def load_de_results(self, filepath):
        try:
            df=pd.read_csv(filepath); return pd.Series(df.log2FoldChange.values, index=df.gene_name).to_dict()
        except Exception as e: print(f"Warning: Could not load DE results: {e}"); return {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Interactively morph structures with advanced coloring and querying.")
    # ... [Argument parsing section with new arguments added] ...
    parser.add_argument("pdb1"); parser.add_argument("scores1"); parser.add_argument("pdb2"); parser.add_argument("scores2")
    parser.add_argument("--name1", default="Cond1"); parser.add_argument("--name2", default="Cond2")
    parser.add_argument('--gtf', required=True); parser.add_argument('--bin-gene-map', required=True)
    parser.add_argument('--rna-seq'); parser.add_argument('--rsem-files', nargs='+')
    parser.add_argument("-lw", "--linewidth", type=float, default=2.0); parser.add_argument("--skip-range", nargs=2, type=int)
    parser.add_argument("--focus-range", nargs=2, type=int); parser.add_argument("-t", "--threshold", type=float, default=0.1)
    parser.add_argument('--highlight-switches', action='store_true', help="Highlight A->B and B->A switches.")
    # MODIFIED: Add new arguments for displaying bins
    parser.add_argument('--show-bins', action='store_true', help="Display bins as small spheres on the structure.")
    parser.add_argument('--bin-size', type=float, default=5.0, help="Size of the spheres if --show-bins is used. Default: 5.0")
    
    args = parser.parse_args()
    # Encapsulate the entire logic in the class
    viewer = Interactive3DViewer(args)