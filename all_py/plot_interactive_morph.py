#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider

# --- Helper Functions (self-contained) ---
# ... [The helper functions like ReadPDB, CalRotateM, PanHandler, etc., remain exactly the same as the previous version] ...
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
def ReadPDB(filename):
    x, y, z = [], [], []
    try:
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith("ATOM"):
                    x.append(float(line[30:38])); y.append(float(line[38:46])); z.append(float(line[46:54]))
    except FileNotFoundError: print(f"Error: PDB file '{filename}' not found."); sys.exit(1)
    except (ValueError, IndexError): print(f"Error parsing coordinates from '{filename}'."); sys.exit(1)
    return x, y, z, None # Return None for connectivity as we handle it sequentially
class PanHandler:
    def __init__(self, fig, ax): self.fig, self.ax, self.press, self.panning = fig, ax, None, False; self.connect()
    def connect(self):
        self.cid_press = self.fig.canvas.mpl_connect('button_press_event', self.on_press)
        self.cid_release = self.fig.canvas.mpl_connect('button_release_event', self.on_release)
        self.cid_motion = self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)
    def on_press(self, event):
        if event.inaxes != self.ax or event.button != 2: return
        self.panning, self.press = True, (event.x, event.y)
    def on_motion(self, event):
        if not self.panning or event.inaxes != self.ax: return
        dx, dy = (event.x - self.press[0])*-1, (event.y - self.press[1])*-1
        xlim, ylim = self.ax.get_xlim(), self.ax.get_ylim()
        self.ax.set_xlim(xlim[0] + dx*(xlim[1]-xlim[0])*0.003, xlim[1] + dx*(xlim[1]-xlim[0])*0.003)
        self.ax.set_ylim(ylim[0] + dy*(ylim[1]-ylim[0])*0.003, ylim[1] + dy*(ylim[1]-ylim[0])*0.003)
        self.press = (event.x, event.y); self.fig.canvas.draw()
    def on_release(self, event):
        if event.button == 2: self.panning, self.press = False, None
def CalRotateM(S1, S2):
    A = np.dot(S1.T, S2); u, _, v = np.linalg.svd(A)
    if np.linalg.det(u)*np.linalg.det(v) < 0: u[:, -1] *= -1
    return np.dot(u, v)
def center_coords(coords): return coords - coords.mean(axis=0)

# --- Main Analysis and Visualization Script ---
def main():
    parser = argparse.ArgumentParser(description="Interactively morph between two PDB structures, with advanced coloring options.")
    # --- Input Files ---
    parser.add_argument("pdb1", help="The starting PDB file (e.g., GM12878).")
    parser.add_argument("scores1", help="Compartment score file for the starting PDB.")
    parser.add_argument("pdb2", help="The ending PDB file (e.g., K562).")
    parser.add_argument("scores2", help="Compartment score file for the ending PDB.")
    # --- Parameters & Output ---
    parser.add_argument("-lw", "--linewidth", type=float, default=2.0, help="Line width.")
    parser.add_argument("--skip-range", nargs=2, type=int, metavar=('S', 'E'), help="1-based bin range to skip.")
    parser.add_argument("--focus-range", nargs=2, type=int, metavar=('S', 'E'), help="1-based bin range to focus on (after skipping).")
    parser.add_argument("-t", "--threshold", type=float, default=0.1, help="Absolute c-score threshold for A/B definition.")
    # MODIFIED: Add a new mode for highlighting
    parser.add_argument('--highlight-switches', action='store_true', help="Enable highlight mode: colors only A->B (blue) and B->A (red) switches, with a gray background.")
    
    args = parser.parse_args()

    # --- Step 1: Load and Process Data ---
    print("--- Step 1: Loading and Processing Data ---")
    coords1 = np.array(ReadPDB(args.pdb1)[:3]).T; scores1_raw = ReadCompartmentScores(args.scores1)
    coords2 = np.array(ReadPDB(args.pdb2)[:3]).T; scores2_raw = ReadCompartmentScores(args.scores2)
    
    if args.skip_range:
        start_idx, end_idx = args.skip_range[0] - 1, args.skip_range[1] - 1
        scores1 = scores1_raw[:start_idx] + scores1_raw[end_idx + 1:]
        scores2 = scores2_raw[:start_idx] + scores2_raw[end_idx + 1:]
    else:
        scores1, scores2 = scores1_raw, scores2_raw
    
    min_len = min(len(coords1), len(scores1), len(coords2), len(scores2))
    coords1, scores1 = coords1[:min_len], np.array(scores1[:min_len])
    coords2, scores2 = coords2[:min_len], np.array(scores2[:min_len])

    # --- Step 2: Align Structures ---
    print("--- Step 2: Aligning Structures ---")
    coords1_centered = center_coords(coords1); coords2_centered = center_coords(coords2)
    rotation_matrix = CalRotateM(coords2_centered, coords1_centered)
    coords2_aligned = np.dot(coords2_centered, rotation_matrix)
    
    if args.focus_range:
        start_bin, end_bin = args.focus_range
        start_idx, end_idx = start_bin - 1, end_bin
        if start_idx < 0 or end_idx > len(coords1):
            print(f"Error: Focus range [{start_bin}-{end_bin}] is out of bounds."); sys.exit(1)
        print(f"--- Step 2.5: Focusing on bin range {start_bin} to {end_bin} ---")
        coords1_centered, coords2_aligned = coords1_centered[start_idx:end_idx], coords2_aligned[start_idx:end_idx]
        scores1, scores2 = scores1[start_idx:end_idx], scores2[start_idx:end_idx]
        
    # --- Step 3: Setup Plot and Colors ---
    print("--- Step 3: Setting up Interactive Plot ---")
    fig = plt.figure(figsize=(9, 9))
    ax = fig.add_axes([0, 0.05, 1, 0.95], projection='3d')
    
    # MODIFIED: Generate colors based on the chosen mode
    bin_categories = []
    for s1, s2 in zip(scores1, scores2):
        state1 = 'A' if s1 > args.threshold else ('B' if s1 < -args.threshold else 'T')
        state2 = 'A' if s2 > args.threshold else ('B' if s2 < -args.threshold else 'T')
        if state1 == 'B' and state2 == 'A': bin_categories.append('B -> A')
        elif state1 == 'A' and state2 == 'B': bin_categories.append('A -> B')
        else: bin_categories.append('Stable/Other')

    if args.highlight_switches:
        colors = np.array(['red' if c == 'B -> A' else 'blue' if c == 'A -> B' else 'lightgray' for c in bin_categories])
        print(f"Highlight mode enabled. Highlighting {np.sum(colors=='red')} 'B->A' bins and {np.sum(colors=='blue')} 'A->B' bins.")
    else:
        colors_t0 = np.array(['red' if s > args.threshold else 'blue' if s < -args.threshold else 'gray' for s in scores1])
        colors_t1 = np.array(['red' if s > args.threshold else 'blue' if s < -args.threshold else 'gray' for s in scores2])

    # --- Step 4: Core Update Function ---
    def update(val):
        while ax.lines: ax.lines[0].remove()
        t = val
        coords_interp = (1 - t) * coords1_centered + t * coords2_aligned
        
        # MODIFIED: Select color array based on mode
        if args.highlight_switches:
            current_colors = colors
        else:
            current_colors = colors_t0 if t < 0.5 else colors_t1

        if len(coords_interp) > 0:
            start_idx, current_color = 0, current_colors[0]
            for i in range(1, len(coords_interp)):
                if current_colors[i] != current_color or i == len(coords_interp) - 1:
                    segment = coords_interp[start_idx : i + 1]
                    ax.plot(segment[:, 0], segment[:, 1], segment[:, 2], color=current_color, linewidth=args.linewidth)
                    start_idx, current_color = i, current_colors[i]
        
        ax.set_title(f"{(1-t)*100:.0f}% GM12878 | {t*100:.0f}% K562")
        fig.canvas.draw_idle()

    # --- Setup Plot View, Slider, and Pan Handler ---
    all_coords = np.vstack([coords1_centered, coords2_aligned])
    if all_coords.size > 0:
        max_range = (all_coords.max(axis=0) - all_coords.min(axis=0)).max() * 0.55
        mid = all_coords.mean(axis=0); ax.set_xlim(mid[0]-max_range, mid[0]+max_range)
        ax.set_ylim(mid[1]-max_range, mid[1]+max_range); ax.set_zlim(mid[2]-max_range, mid[2]+max_range)
    ax.axis('off')
    ax_slider = fig.add_axes([0.2, 0.02, 0.65, 0.03])
    morph_slider = Slider(ax=ax_slider, label='Morph', valmin=0.0, valmax=1.0, valinit=0.0)
    morph_slider.on_changed(update)
    pan_handler = PanHandler(fig, ax)

    print("--- Step 5: Launching Interactive Window ---")
    update(0.0)
    plt.show()

if __name__ == "__main__":
    main()