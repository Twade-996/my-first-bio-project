#!/usr/bin/env python
#-*- coding: utf-8 -*-
'''
This is a simple script for chromosome structure visualization using Matplotlib and structure alignment with Kabsch algorithm.
It can also color the structure based on A/B compartment scores, and skip regions like centromeres.

An example with compartment coloring: python plot.py a.pdb -cf a.scores --show

An example skipping a centromere region in the score file:
python plot.py a.pdb -cf a.scores --skip-range 1215 1425 --show
'''
import argparse
import random

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

from IO import ReadPDB, ReadCompartmentScores

COLORS = ["orange", "purple", "cyan", "magenta", "lime"]
A_COMP_COLOR = "red"
B_COMP_COLOR = "blue"

def CalRotateM(S1, S2):
    A = np.dot(S1.coordinate.T, S2.coordinate)
    u, s, v = np.linalg.svd(A)
    d = np.linalg.det(u) * np.linalg.det(v)
    if d < 0: u[:, -1] *= -1
    return np.dot(u, v)

def rmsd(S1, S2):
    S1_p, S2_p = S1.coordinate, S2.coordinate
    return np.sqrt(np.sum((S1_p - S2_p)**2.0) / S2_p.shape[0])

class Structure:
    is_setline = is_setpoint = False

    def __init__(self, filename=None, ax=None):
        self.scores = None
        if filename: self.AddStructureFromFile(filename)
        if ax: self.ax = ax

    def AddStructureFromFile(self, filename):
        x, y, z, c = ReadPDB(filename)
        self.coordinate = np.array([[x[i], y[i], z[i]] for i in range(len(x))])
        self.Center2O()

    def AddScores(self, scores):
        num_coords = len(self.coordinate)
        num_scores = len(scores)

        if num_scores != num_coords:
            print(f"Warning: Mismatch between coordinates ({num_coords}) and scores ({num_scores}).")
            print("         If using --skip-range, please ensure the range is correct.")
            print("         The plot will be truncated to the smaller of the two lists.")
            min_len = min(num_coords, num_scores)
            self.coordinate = self.coordinate[:min_len]
            self.scores = np.array(scores[:min_len])
        else:
            print(f"Successfully matched {num_coords} coordinates with {num_scores} scores.")
            self.scores = np.array(scores)

    def Center2O(self):
        self.Translation(self.GetCenter())

    def GetCenter(self):
        return self.coordinate.mean(axis=0)

    def Translation(self, pos):
        self.coordinate -= pos

    def Rotate(self, m):
        self.coordinate = np.dot(self.coordinate, m)

    # MODIFIED: Optimized drawing function
    def DrawLine(self):
        """
        Draws the structure. Uses compartment colors if scores are available.
        This version is optimized to draw contiguous segments of the same color in a single plot call.
        """
        if not self.is_setline: self.SetLineArgs()

        # If scores are not available or empty, draw with a single color
        if self.scores is None or len(self.scores) == 0:
            self.ax.plot(self.coordinate[:, 0], self.coordinate[:, 1], self.coordinate[:, 2], color=self.line_color, linewidth=self.line_size)
            return

        # Optimized drawing logic for compartments
        coords = self.coordinate
        start_index = 0
        current_color = A_COMP_COLOR if self.scores[0] > 0 else B_COMP_COLOR
        
        for i in range(1, len(coords)):
            score = self.scores[i]
            new_color = A_COMP_COLOR if score > 0 else B_COMP_COLOR

            # If the color changes, or we are at the end, plot the previous segment
            if new_color != current_color or i == len(coords) - 1:
                # The segment to plot is from start_index to i (inclusive)
                segment = coords[start_index : i + 1]
                self.ax.plot(segment[:, 0], segment[:, 1], segment[:, 2], color=current_color, linewidth=self.line_size)
                
                # Start a new segment
                start_index = i
                current_color = new_color


    def SetLineArgs(self, color=None, size=None):
        self.line_color = color if color else "gray"
        self.line_size = float(size) if size else 1.0
        self.is_setline = True

    def SetDrawArgs(self, color=None, l_size=None):
        self.SetLineArgs(color, l_size)

    def Draw(self, line=True, point=False):
        if line: self.DrawLine()
        
# --- Main execution block remains the same as the previous version ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize PDB structures and optionally color by compartment score.")
    parser.add_argument("i", help="PDB file(s), split with ':'.")
    parser.add_argument("-o", help="Output picture file name.")
    parser.add_argument("-lw", "--linewidth", help="Line width(s). Default: 2.0. Split with ':'.", default="2.0")
    parser.add_argument("-c", "--color", help="Default color(s) if not using compartment file. Split with ':'.", default="gray")
    parser.add_argument("-cf", "--compartment_file", help="A/B compartment score file(s). Split with ':'.")
    
    parser.add_argument("--skip-range", nargs=2, type=int, metavar=('START_BIN', 'END_BIN'),
                        help="The 1-based start and end bin range to skip in the score file (e.g., for a centromere). Example: --skip-range 1215 1425")

    parser.add_argument("--withsuperposition", action='store_true', help="Superpose structures.")
    parser.add_argument("--withRMSD", action='store_true', help="Print RMSD values.")
    parser.add_argument("--show", action='store_true', help="Show plot on screen.")
    
    args = parser.parse_args()

    pdbs = args.i.split(':')
    line_widths = args.linewidth.split(':')
    colors = args.color.split(':')
    compartment_files = args.compartment_file.split(':') if args.compartment_file else []

    if len(line_widths) < len(pdbs): line_widths.extend([line_widths[-1]] * (len(pdbs) - len(line_widths)))
    if len(colors) < len(pdbs): colors.extend([colors[-1]] * (len(pdbs) - len(colors)))

    figure = plt.figure(figsize=(8, 8))
    ax = Axes3D(figure)

    structures = []
    for i, pdb in enumerate(pdbs):
        s = Structure(pdb, ax)
        if i < len(compartment_files):
            print(f"Loading scores from '{compartment_files[i]}' for '{pdb}'...")
            scores = ReadCompartmentScores(compartment_files[i])
            
            if args.skip_range:
                start_bin, end_bin = args.skip_range
                start_index = start_bin - 1
                end_index = end_bin - 1
                
                print(f"Original score count: {len(scores)}")
                print(f"Skipping score range for bins {start_bin}-{end_bin} (indices {start_index}-{end_index})")
                
                filtered_scores = scores[:start_index] + scores[end_index + 1:]
                scores = filtered_scores
                print(f"Filtered score count: {len(scores)}")
            
            s.AddScores(scores)
        structures.append(s)

    for i, s in enumerate(structures):
        s.SetDrawArgs(color=colors[i], l_size=line_widths[i])

    if args.withsuperposition and len(structures) > 1:
        s1 = structures[0]
        for i in range(1, len(structures)):
            r_m = CalRotateM(structures[i], s1)
            structures[i].Rotate(r_m)

    if args.withRMSD and len(structures) > 1:
        print("\nStructure1\tStructure2\tRMSD value")
        for i in range(len(structures)):
            for j in range(i + 1, len(structures)):
                rmsd_value = rmsd(structures[i], structures[j])
                print(f"{pdbs[i]}\t{pdbs[j]}\t{rmsd_value:.6f}")

    print("\nDrawing structure(s)...")
    for s in structures:
        s.Draw()
    
    plt.axis('off')

    if args.o:
        print(f"Saving figure to {args.o}")
        plt.savefig(args.o, dpi=300, bbox_inches='tight')

    if args.show:
        plt.show()
