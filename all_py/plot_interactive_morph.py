#!/usr/bin/env python
#-*- coding: utf-8 -*-
'''
Interactively explore the morph between two chromosome structures using a slider.
This version includes full interactive controls:
- Left Mouse: Rotate
- Right Mouse: Zoom
- Middle Mouse: Pan (Translate)
- Slider: Morph between structures
'''
import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider

from IO import ReadPDB, ReadCompartmentScores

A_COMP_COLOR = "red"
B_COMP_COLOR = "blue"

# MODIFIED: Added a dedicated class to handle panning
class PanHandler:
    """A class to handle panning the 3D plot with the middle mouse button."""
    def __init__(self, fig, ax):
        self.fig = fig
        self.ax = ax
        self.press = None
        self.panning = False
        
        self.connect()

    def connect(self):
        """Connects the mouse events to their callback functions."""
        self.cid_press = self.fig.canvas.mpl_connect('button_press_event', self.on_press)
        self.cid_release = self.fig.canvas.mpl_connect('button_release_event', self.on_release)
        self.cid_motion = self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def on_press(self, event):
        """Callback for when a mouse button is pressed."""
        if event.inaxes != self.ax: return
        # Check for middle mouse button (button == 2)
        if event.button == 2:
            self.panning = True
            self.press = (event.x, event.y)

    def on_motion(self, event):
        """Callback for when the mouse is moved."""
        if not self.panning or event.inaxes != self.ax: return
        
        # Get the change in mouse position
        dx = event.x - self.press[0]
        dy = event.y - self.press[1]
        
        # Invert the direction for a natural feel
        dx *= -1
        dy *= -1
        
        # Get current limits
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()
        
        # Estimate a scaling factor based on the view range
        x_range = xlim[1] - xlim[0]
        y_range = ylim[1] - ylim[0]
        
        # Adjust the pan speed based on the size of what's on screen
        scale_factor = 0.003
        
        # Update the limits
        self.ax.set_xlim(xlim[0] + dx * x_range * scale_factor, xlim[1] + dx * x_range * scale_factor)
        self.ax.set_ylim(ylim[0] + dy * y_range * scale_factor, ylim[1] + dy * y_range * scale_factor)
        
        # Store current mouse position for the next motion event
        self.press = (event.x, event.y)
        
        self.fig.canvas.draw()

    def on_release(self, event):
        """Callback for when a mouse button is released."""
        if event.button == 2:
            self.panning = False
            self.press = None

def CalRotateM(S1, S2):
    A = np.dot(S1.T, S2)
    u, _, v = np.linalg.svd(A)
    d = np.linalg.det(u) * np.linalg.det(v)
    if d < 0: u[:, -1] *= -1
    return np.dot(u, v)

def center_coords(coords):
    return coords - coords.mean(axis=0)

def main():
    parser = argparse.ArgumentParser(description="Interactively morph between two PDB structures using a slider.")
    parser.add_argument("pdb1", help="The starting PDB file (e.g., GM12878).")
    parser.add_argument("scores1", help="The compartment score file for the starting PDB.")
    parser.add_argument("pdb2", help="The ending PDB file (e.g., K562).")
    parser.add_argument("scores2", help="The compartment score file for the ending PDB.")
    parser.add_argument("-lw", "--linewidth", type=float, default=2.0, help="Line width.")
    parser.add_argument("--skip-range", nargs=2, type=int, metavar=('START', 'END'),
                        help="1-based bin range to skip in score files (e.g., for centromeres).")
    parser.add_argument("--focus-range", nargs=2, type=int, metavar=('START_BIN', 'END_BIN'),
                        help="1-based bin range to focus on. These are bin numbers *after* skipping the centromere.")

    args = parser.parse_args()

    # --- Data Loading and Processing ---
    print("--- Step 1: Loading and Processing Data ---")
    x1, y1, z1, _ = ReadPDB(args.pdb1); coords1 = np.array([x1, y1, z1]).T
    scores1_raw = ReadCompartmentScores(args.scores1)
    x2, y2, z2, _ = ReadPDB(args.pdb2); coords2 = np.array([x2, y2, z2]).T
    scores2_raw = ReadCompartmentScores(args.scores2)

    if args.skip_range:
        start_idx, end_idx = args.skip_range[0] - 1, args.skip_range[1] - 1
        scores1 = scores1_raw[:start_idx] + scores1_raw[end_idx + 1:]
        scores2 = scores2_raw[:start_idx] + scores2_raw[end_idx + 1:]
    else:
        scores1, scores2 = scores1_raw, scores2_raw

    min_len = min(len(coords1), len(scores1), len(coords2), len(scores2))
    print(f"Common length is {min_len} bins. Using this length for all data.")
    coords1, scores1 = coords1[:min_len], np.array(scores1[:min_len])
    coords2, scores2 = coords2[:min_len], np.array(scores2[:min_len])

    # --- Structure Alignment ---
    print("--- Step 2: Aligning Structures ---")
    coords1_centered = center_coords(coords1); coords2_centered = center_coords(coords2)
    rotation_matrix = CalRotateM(coords2_centered, coords1_centered)
    coords2_aligned = np.dot(coords2_centered, rotation_matrix)
    
    if args.focus_range:
        start_bin, end_bin = args.focus_range
        start_idx, end_idx = start_bin - 1, end_bin
        if start_idx < 0 or end_idx > len(coords1):
            print(f"Error: Focus range [{start_bin}-{end_bin}] is out of bounds for data of length {len(coords1)}.")
            sys.exit(1)
        print(f"--- Step 2.5: Focusing on bin range {start_bin} to {end_bin} ---")
        coords1_centered, coords2_aligned = coords1_centered[start_idx:end_idx], coords2_aligned[start_idx:end_idx]
        scores1, scores2 = scores1[start_idx:end_idx], scores2[start_idx:end_idx]
        
    # --- Interactive Plotting Setup ---
    print("--- Step 3: Setting up Interactive Plot ---")
    fig = plt.figure(figsize=(9, 9))
    ax = fig.add_axes([0, 0.05, 1, 0.95], projection='3d')
    
    colors_t0 = np.array([A_COMP_COLOR if s > 0 else B_COMP_COLOR for s in scores1])
    colors_t1 = np.array([A_COMP_COLOR if s > 0 else B_COMP_COLOR for s in scores2])

    def update(val):
        while ax.lines: ax.lines[0].remove()
        t = val
        coords_interp = (1 - t) * coords1_centered + t * coords2_aligned
        current_colors = colors_t0 if t < 0.5 else colors_t1
        if len(coords_interp) > 0:
            start_idx, current_color = 0, current_colors[0]
            for i in range(1, len(coords_interp)):
                if current_colors[i] != current_color or i == len(coords_interp) - 1:
                    segment = coords_interp[start_idx : i + 1]
                    ax.plot(segment[:, 0], segment[:, 1], segment[:, 2], color=current_color, linewidth=args.linewidth)
                    start_idx, current_color = i, current_colors[i]
        
        cell_line1, cell_line2 = args.pdb1.split('/')[-1].split('_')[0], args.pdb2.split('/')[-1].split('_')[0]
        ax.set_title(f"{(1-t)*100:.0f}% {cell_line1} | {t*100:.0f}% {cell_line2}")
        fig.canvas.draw_idle()

    all_coords = np.vstack([coords1_centered, coords2_aligned])
    if all_coords.size > 0:
        max_range = (all_coords.max(axis=0) - all_coords.min(axis=0)).max() * 0.55
        mid = all_coords.mean(axis=0)
        ax.set_xlim(mid[0] - max_range, mid[0] + max_range)
        ax.set_ylim(mid[1] - max_range, mid[1] + max_range)
        ax.set_zlim(mid[2] - max_range, mid[2] + max_range)
    ax.axis('off')

    ax_slider = fig.add_axes([0.2, 0.02, 0.65, 0.03])
    morph_slider = Slider(ax=ax_slider, label='Morph', valmin=0.0, valmax=1.0, valinit=0.0)
    morph_slider.on_changed(update)
    
    # MODIFIED: Activate the PanHandler
    pan_handler = PanHandler(fig, ax)

    print("--- Step 4: Launching Interactive Window ---")
    update(0.0)
    plt.show()

if __name__ == "__main__":
    main()
