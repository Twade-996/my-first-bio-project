#!/usr/bin/env python
#-*- coding: utf-8 -*-
'''
Interactively morph between structures, with A/B/Transitional compartments.
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
TRANSITION_COLOR = "DarkSeaGreen"

class PanHandler: 
    def __init__(self, fig, ax):
        self.fig, self.ax, self.press, self.panning = fig, ax, None, False
        self.connect()
    def connect(self):
        self.cid_press = self.fig.canvas.mpl_connect('button_press_event', self.on_press)
        self.cid_release = self.fig.canvas.mpl_connect('button_release_event', self.on_release)
        self.cid_motion = self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)
    def on_press(self, event):
        if event.inaxes != self.ax or event.button != 2: return
        self.panning, self.press = True, (event.x, event.y)
    def on_motion(self, event):
        if not self.panning or event.inaxes != self.ax: return
        dx, dy = (event.x - self.press[0]) * -1, (event.y - self.press[1]) * -1
        xlim, ylim = self.ax.get_xlim(), self.ax.get_ylim()
        x_range, y_range = xlim[1] - xlim[0], ylim[1] - ylim[0]
        scale = 0.003
        self.ax.set_xlim(xlim[0] + dx * x_range * scale, xlim[1] + dx * x_range * scale)
        self.ax.set_ylim(ylim[0] + dy * y_range * scale, ylim[1] + dy * y_range * scale)
        self.press = (event.x, event.y)
        self.fig.canvas.draw()
    def on_release(self, event):
        if event.button == 2: self.panning, self.press = False, None

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
    parser.add_argument("pdb1", help="Starting PDB file.")
    parser.add_argument("scores1", help="Compartment score file for PDB1.")
    parser.add_argument("pdb2", help="Ending PDB file.")
    parser.add_argument("scores2", help="Compartment score file for PDB2.")
    parser.add_argument("-lw", "--linewidth", type=float, default=2.0, help="Line width.")
    parser.add_argument("--skip-range", nargs=2, type=int, metavar=('S', 'E'), help="1-based bin range to skip in score files.")
    parser.add_argument("--focus-range", nargs=2, type=int, metavar=('S', 'E'), help="1-based bin range to focus on (after skipping).")
    # MODIFIED: Add threshold argument
    parser.add_argument("-t", "--threshold", type=float, default=0.0, help="Absolute c-score value to define A/B compartments. Default: 0.0")
    
    args = parser.parse_args()

    print("--- Step 1: Loading and Processing Data ---")
    coords1 = np.array(ReadPDB(args.pdb1)[:3]).T
    scores1_raw = ReadCompartmentScores(args.scores1)
    coords2 = np.array(ReadPDB(args.pdb2)[:3]).T
    scores2_raw = ReadCompartmentScores(args.scores2)

    if args.skip_range:
        start_idx, end_idx = args.skip_range[0] - 1, args.skip_range[1] - 1
        scores1 = scores1_raw[:start_idx] + scores1_raw[end_idx + 1:]
        scores2 = scores2_raw[:start_idx] + scores2_raw[end_idx + 1:]
    else:
        scores1, scores2 = scores1_raw, scores2_raw

    min_len = min(len(coords1), len(scores1), len(coords2), len(scores2))
    coords1, scores1 = coords1[:min_len], np.array(scores1[:min_len])
    coords2, scores2 = coords2[:min_len], np.array(scores2[:min_len])

    print("--- Step 2: Aligning Structures ---")
    coords1_centered = center_coords(coords1); coords2_centered = center_coords(coords2)
    rotation_matrix = CalRotateM(coords2_centered, coords1_centered)
    coords2_aligned = np.dot(coords2_centered, rotation_matrix)
    
    if args.focus_range:
        start_idx, end_idx = args.focus_range[0] - 1, args.focus_range[1]
        coords1_centered, coords2_aligned = coords1_centered[start_idx:end_idx], coords2_aligned[start_idx:end_idx]
        scores1, scores2 = scores1[start_idx:end_idx], scores2[start_idx:end_idx]
        
    print("--- Step 3: Setting up Interactive Plot ---")
    fig = plt.figure(figsize=(9, 9))
    ax = fig.add_axes([0, 0.05, 1, 0.95], projection='3d')
    
    # MODIFIED: Color generation now includes the threshold
    def get_colors_with_threshold(scores, threshold):
        colors = []
        for s in scores:
            if s > threshold: colors.append(A_COMP_COLOR)
            elif s < -threshold: colors.append(B_COMP_COLOR)
            else: colors.append(TRANSITION_COLOR)
        return np.array(colors)

    colors_t0 = get_colors_with_threshold(scores1, args.threshold)
    colors_t1 = get_colors_with_threshold(scores2, args.threshold)

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
    
    pan_handler = PanHandler(fig, ax)

    print("--- Step 4: Launching Interactive Window ---")
    update(0.0)
    plt.show()

if __name__ == "__main__":
    main()