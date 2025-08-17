#!/usr/bin/env python
#-*- coding: utf-8 -*-
'''
Creates a dynamic animation (morph) between two chromosome structures,
coloring them by A/B compartment status and highlighting their changes.
This version automatically handles PDBs with slightly different bin counts by truncating to the shorter length.
'''
import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

from IO import ReadPDB, ReadCompartmentScores

A_COMP_COLOR = "red"
B_COMP_COLOR = "blue"

def CalRotateM(S1, S2):
    A = np.dot(S1.T, S2)
    u, _, v = np.linalg.svd(A)
    d = np.linalg.det(u) * np.linalg.det(v)
    if d < 0: u[:, -1] *= -1
    return np.dot(u, v)

def center_coords(coords):
    return coords - coords.mean(axis=0)

def main():
    parser = argparse.ArgumentParser(description="Create a morphing animation between two PDB structures, colored by compartment scores.")
    parser.add_argument("pdb1", help="The starting PDB file (e.g., GM12878).")
    parser.add_argument("scores1", help="The compartment score file for the starting PDB.")
    parser.add_argument("pdb2", help="The ending PDB file (e.g., K562).")
    parser.add_argument("scores2", help="The compartment score file for the ending PDB.")
    parser.add_argument("-o", help="Output video file name (e.g., morph.mp4).", required=True)
    parser.add_argument("-f", "--frames", type=int, default=100, help="Number of frames in the animation.")
    parser.add_argument("-lw", "--linewidth", type=float, default=2.0, help="Line width.")
    parser.add_argument("--skip-range", nargs=2, type=int, metavar=('START', 'END'),
                        help="1-based bin range to skip in score files (e.g., for centromeres).")
    
    args = parser.parse_args()

    print("--- Step 1: Loading Data ---")
    x1, y1, z1, _ = ReadPDB(args.pdb1)
    coords1 = np.array([x1, y1, z1]).T
    scores1_raw = ReadCompartmentScores(args.scores1)

    x2, y2, z2, _ = ReadPDB(args.pdb2)
    coords2 = np.array([x2, y2, z2]).T
    scores2_raw = ReadCompartmentScores(args.scores2)

    print(f"Initial counts: PDB1={len(coords1)}, Scores1={len(scores1_raw)}, PDB2={len(coords2)}, Scores2={len(scores2_raw)}")

    print("\n--- Step 2: Filtering Scores (if needed) ---")
    if args.skip_range:
        start_idx, end_idx = args.skip_range[0] - 1, args.skip_range[1] - 1
        print(f"Skipping bins {args.skip_range[0]}-{args.skip_range[1]} from score files.")
        scores1 = scores1_raw[:start_idx] + scores1_raw[end_idx + 1:]
        scores2 = scores2_raw[:start_idx] + scores2_raw[end_idx + 1:]
        print(f"Filtered score counts: Scores1={len(scores1)}, Scores2={len(scores2)}")
    else:
        scores1, scores2 = scores1_raw, scores2_raw

    print("\n--- Step 3: Validating and Truncating Data ---")
    # MODIFIED: Find the minimum length across all datasets and truncate
    min_len = min(len(coords1), len(scores1), len(coords2), len(scores2))
    print(f"Common minimum length found: {min_len} bins. Truncating all datasets to this length.")

    coords1 = coords1[:min_len]
    scores1 = np.array(scores1[:min_len])
    coords2 = coords2[:min_len]
    scores2 = np.array(scores2[:min_len])

    if min_len == 0:
        print("Error: Common length is zero after filtering. Cannot proceed.")
        sys.exit(1)

    print("\n--- Step 4: Aligning Structures ---")
    coords1_centered = center_coords(coords1)
    coords2_centered = center_coords(coords2)
    
    rotation_matrix = CalRotateM(coords2_centered, coords1_centered)
    coords2_aligned = np.dot(coords2_centered, rotation_matrix)

    print("\n--- Step 5: Preparing Animation ---")
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    colors_t0 = np.array([A_COMP_COLOR if s > 0 else B_COMP_COLOR for s in scores1])
    colors_t1 = np.array([A_COMP_COLOR if s > 0 else B_COMP_COLOR for s in scores2])

    def update(frame):
        ax.cla()
        t = frame / (args.frames - 1)
        coords_interp = (1 - t) * coords1_centered + t * coords2_aligned
        current_colors = colors_t0 if t < 0.5 else colors_t1

        start_idx = 0
        current_color = current_colors[0]
        for i in range(1, len(coords_interp)):
            if current_colors[i] != current_color or i == len(coords_interp) - 1:
                segment = coords_interp[start_idx : i + 1]
                ax.plot(segment[:, 0], segment[:, 1], segment[:, 2], color=current_color, linewidth=args.linewidth)
                start_idx = i
                current_color = current_colors[i]
        
        all_coords = np.vstack([coords1_centered, coords2_aligned])
        max_range = (all_coords.max(axis=0) - all_coords.min(axis=0)).max() * 0.55
        mid = all_coords.mean(axis=0)
        ax.set_xlim(mid[0] - max_range, mid[0] + max_range)
        ax.set_ylim(mid[1] - max_range, mid[1] + max_range)
        ax.set_zlim(mid[2] - max_range, mid[2] + max_range)
        ax.axis('off')
        
        cell_line1 = args.pdb1.split('_')[0]
        cell_line2 = args.pdb2.split('_')[0]
        title_text = cell_line1 if t < 0.5 else cell_line2
        ax.set_title(f"Morphing to: {title_text} (Frame {frame+1}/{args.frames})")

    print(f"\n--- Step 6: Creating and Saving Animation ({args.frames} frames) ---")
    anim = FuncAnimation(fig, update, frames=args.frames, interval=50)
    
    try:
        anim.save(args.o, writer='ffmpeg', dpi=150)
        print(f"\nAnimation successfully saved to {args.o}")
    except Exception as e:
        print("\n--- ERROR ---")
        print("Failed to save animation. This usually means 'ffmpeg' is not installed.")
        print("Please install it. If you use Conda, run: 'conda install -c conda-forge ffmpeg'")
        print(f"Original error: {e}")

if __name__ == "__main__":
    main()
