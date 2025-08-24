#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse, re, sys, numpy as np, pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider
from functools import reduce

# --- Helper Functions (unchanged) ---
def ReadPDB(filename):
    x, y, z = [], [], []; f=open(filename, 'r'); [ (x.append(float(l[30:38])), y.append(float(l[38:46])), z.append(float(l[46:54]))) for l in f if l.startswith("ATOM") ]; f.close(); return x, y, z, None
def ReadCompartmentScores(filename):
    s=[]; f=open(filename, 'r'); [ s.append(float(l.strip().split()[1])) for l in f if len(l.strip().split()) >= 2 ]; f.close(); return s
class PanHandler:
    def __init__(self, fig, ax): self.fig, self.ax, self.press, self.panning = fig, ax, None, False; self.connect()
    def connect(self): self.cid_press=self.fig.canvas.mpl_connect('button_press_event',self.on_press); self.cid_release=self.fig.canvas.mpl_connect('button_release_event',self.on_release); self.cid_motion=self.fig.canvas.mpl_connect('motion_notify_event',self.on_motion)
    def on_press(self, e):
        if e.inaxes != self.ax or e.button != 2: return
        self.panning, self.press = True, (e.x, e.y)
    def on_motion(self, e):
        if not self.panning or e.inaxes != self.ax: return
        dx, dy = (e.x - self.press[0])*-1, (e.y - self.press[1])*-1; xlim, ylim = self.ax.get_xlim(), self.ax.get_ylim()
        self.ax.set_xlim(xlim[0] + dx*(xlim[1]-xlim[0])*0.003, xlim[1] + dx*(xlim[1]-xlim[0])*0.003); self.ax.set_ylim(ylim[0] + dy*(ylim[1]-ylim[0])*0.003, ylim[1] + dy*(ylim[1]-ylim[0])*0.003)
        self.press = (e.x, e.y); self.fig.canvas.draw()
    def on_release(self, e):
        if e.button == 2: self.panning, self.press = False, None
def center_coords(coords): return coords - coords.mean(axis=0)
def CalRotateM(S1, S2):
    A=np.dot(S1.T,S2); u,_,v=np.linalg.svd(A)
    if np.linalg.det(u)*np.linalg.det(v)<0: u[:,-1]*=-1
    return np.dot(u,v)

# --- The Main Application Class ---
class ShowcaseViewer:
    def __init__(self, args):
        self.args = args
        self.load_and_process_data()
        self.identify_highlight_bins() # Renamed function
        self.setup_plot()

    @staticmethod
    def parse_bin_gene_map(map_file): # ... (Unchanged)
        bin_to_gene = {}; print(f"\n--- Parsing bin-to-gene map from '{map_file}' ---")
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

    def load_and_process_data(self): # ... (Unchanged)
        print("--- Step 1: Loading all data ---")
        self.bin_to_gene = self.parse_bin_gene_map(self.args.bin_gene_map)
        self.de_map = self.load_de_results(self.args.rna_seq)
        scores1_raw, scores2_raw = ReadCompartmentScores(self.args.scores1), ReadCompartmentScores(self.args.scores2)
        if self.args.skip_range:
            start, end = self.args.skip_range[0]-1, self.args.skip_range[1]-1
            scores1, scores2 = scores1_raw[:start]+scores1_raw[end+1:], scores2_raw[:start]+scores2_raw[end+1:]
        else: scores1, scores2 = scores1_raw, scores2_raw
        coords1 = np.array(ReadPDB(self.args.pdb1)[:3]).T
        coords2 = np.array(ReadPDB(self.args.pdb2)[:3]).T
        min_len = min(len(coords1), len(scores1), len(coords2), len(scores2))
        self.scores1, self.scores2 = np.array(scores1[:min_len]), np.array(scores2[:min_len])
        self.coords1_centered = center_coords(coords1[:min_len])
        coords2_centered = center_coords(coords2[:min_len])
        self.coords2_aligned = np.dot(coords2_centered, CalRotateM(coords2_centered, self.coords1_centered))
        self.bin_categories = self.categorize_bins()

    # MODIFIED: This function now finds ALL relevant bins, not just the top N
    def identify_highlight_bins(self):
        """Identifies all bins that match the specified highlight criteria."""
        print(f"\n--- Step 2: Identifying all '{self.args.highlight_type}' bins with gene data ---")
        self.highlight_bins = []
        target_category = self.args.highlight_type

        for i, category in enumerate(self.bin_categories):
            if category == target_category and i in self.bin_to_gene:
                # Check if at least one gene in the bin has RNA-seq data
                genes_in_bin = self.bin_to_gene[i]
                found_gene_with_data = any(gene in self.de_map for gene in genes_in_bin)
                
                if found_gene_with_data:
                    self.highlight_bins.append({'bin_index': i})
        
        if not self.highlight_bins:
            print(f"Warning: Could not find any '{target_category}' bins containing genes with expression data.")
        else:
            print(f"Found {len(self.highlight_bins)} bins to highlight.")
            # Optionally print the first few found bins
            for b in self.highlight_bins[:5]:
                print(f"  - Bin {b['bin_index']+1}")
            if len(self.highlight_bins) > 5: print("  - ... and more.")

    def setup_plot(self): # ... (Unchanged)
        self.fig = plt.figure(figsize=(12, 9)); self.ax_3d = self.fig.add_axes([0.0, 0.05, 0.75, 0.95], projection='3d')
        self.ax_legend = self.fig.add_axes([0.75, 0.1, 0.2, 0.8]); self.draw_legend()
        self.fig.canvas.mpl_connect('button_press_event', self.on_legend_click)
        ax_slider = self.fig.add_axes([0.15, 0.02, 0.5, 0.03])
        self.morph_slider = Slider(ax=ax_slider, label='Morph', valmin=0.0, valmax=1.0, valinit=0.0)
        self.morph_slider.on_changed(self.update); self.pan_handler = PanHandler(self.fig, self.ax_3d)
        self.update(0.0); plt.show()
    
    def draw_legend(self):
        self.ax_legend.set_title(f"Highlighted '{self.args.highlight_type}' Bins", fontsize=12)
        self.ax_legend.set_axis_off()
        
        # MODIFIED: We only use one color for all highlighted bins
        highlight_color = 'red' if self.args.highlight_type == 'B -> A' else 'blue'
        
        self.legend_patches = []
        # MODIFIED: The legend now might need to be scrollable if too many bins are found.
        # For simplicity, we'll just list the first ~15-20.
        max_legend_items = 20 
        for i, bin_info in enumerate(self.highlight_bins[:max_legend_items]):
            y_pos = 0.95 - i * 0.045 # Denser packing for more items
            # Create a simple dot instead of a big rectangle
            self.ax_legend.plot(0.05, y_pos, 'o', markersize=8, color=highlight_color)
            label = f"Bin {bin_info['bin_index']+1}"
            text_obj = self.ax_legend.text(0.15, y_pos, label, va='center', fontsize=9)
            # Store the text object to check for clicks later
            self.legend_patches.append({'text_obj': text_obj, 'info': bin_info})
        if len(self.highlight_bins) > max_legend_items:
            self.ax_legend.text(0.1, 0.95 - max_legend_items * 0.045, f"...and {len(self.highlight_bins)-max_legend_items} more.", fontsize=9)


    def on_legend_click(self, event):
        if event.inaxes != self.ax_legend: return
        for item in self.legend_patches:
            # Check if the click was inside this specific text's bounding box
            if item['text_obj'].get_window_extent().contains(event.x, event.y):
                bin_info = item['info']; bin_index = bin_info['bin_index']
                genes = self.bin_to_gene.get(bin_index, ["No genes found"])
                print("\n" + "="*40 + f"\nDetails for Highlighted Bin {bin_index+1}\n" + "="*40)
                print(" All Genes in Bin:")
                for gene in sorted(genes):
                    log2fc = self.de_map.get(gene, "N/A")
                    log2fc_str = f"{log2fc:.2f}" if isinstance(log2fc, float) else "N/A"
                    print(f"  - {gene:<15} | log2FC: {log2fc_str}")
                print("="*40)
                return

    def update(self, val):
        self.ax_3d.cla(); t = val
        coords_interp = (1 - t) * self.coords1_centered + t * self.coords2_aligned
        colors = np.array(['lightgray'] * len(self.bin_categories), dtype=object)
        
        # MODIFIED: Color all identified highlight bins with one color
        highlight_color = 'red' if self.args.highlight_type == 'B -> A' else 'blue'
        highlight_indices = [b['bin_index'] for b in self.highlight_bins]
        for idx in highlight_indices: colors[idx] = highlight_color
        
        # ... (Plotting logic is the same) ...
        if len(coords_interp) > 0:
            for i in range(len(coords_interp) - 1):
                p1, p2 = coords_interp[i], coords_interp[i+1]; self.ax_3d.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], color=colors[i], linewidth=self.args.linewidth)
        max_range=(np.ptp(self.coords1_centered,axis=0).max())*0.55; mid=self.coords1_centered.mean(axis=0)
        self.ax_3d.set_xlim(mid[0]-max_range, mid[0]+max_range); self.ax_3d.set_ylim(mid[1]-max_range, mid[1]+max_range); self.ax_3d.set_zlim(mid[2]-max_range, mid[2]+max_range)
        self.ax_3d.axis('off'); self.ax_3d.set_title(f"{(1-t)*100:.0f}% {self.args.name1} | {t*100:.0f}% {self.args.name2}")
        self.fig.canvas.draw_idle()
    
    def categorize_bins(self): # ... (Unchanged)
        categories = [];
        for s1, s2 in zip(self.scores1, self.scores2):
            state1 = 'A' if s1>self.args.threshold else ('B' if s1<-self.args.threshold else 'T')
            state2 = 'A' if s2>self.args.threshold else ('B' if s2<-self.args.threshold else 'T')
            if state1=='B' and state2=='A': categories.append('B -> A')
            elif state1=='A' and state2=='B': categories.append('A -> B')
            else: categories.append('Stable/Other')
        return categories

    def load_de_results(self, filepath): # ... (Unchanged)
        try:
            df = pd.read_csv(filepath); return pd.Series(df.log2FoldChange.values, index=df.gene_name).to_dict()
        except Exception as e: print(f"Warning: Could not load DE results: {e}"); return {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Showcase gene switches in an interactive 3D morph.")
    parser.add_argument("pdb1"); parser.add_argument("scores1"); parser.add_argument("pdb2"); parser.add_argument("scores2")
    parser.add_argument("--name1", default="Cond1"); parser.add_argument("--name2", default="Cond2")
    parser.add_argument('--gtf', required=True); parser.add_argument('--bin-gene-map', required=True)
    parser.add_argument('--rna-seq', required=True)
    parser.add_argument("-lw", "--linewidth", type=float, default=2.0)
    parser.add_argument("--skip-range", nargs=2, type=int)
    parser.add_argument("-t", "--threshold", type=float, default=0.1)
    
    # MODIFIED: Replaced --top-n with a more descriptive choice
    parser.add_argument("--highlight-type", choices=['B -> A', 'A -> B'], default='B -> A', 
                        help="Specify which type of compartment switch to highlight. Default: B -> A")
    
    args = parser.parse_args()
    viewer = ShowcaseViewer(args)