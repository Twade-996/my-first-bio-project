#!/usr/bin/env python
#-*- coding: utf-8 -*-
'''
Visualizes chromosome structures with A/B/Transitional compartment coloring
based on a c-score threshold.
'''
import argparse
import random

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

from IO import ReadPDB, ReadCompartmentScores

# Define global colors for clarity
A_COMP_COLOR = "red"
B_COMP_COLOR = "blue"
TRANSITION_COLOR = "DarkSeaGreen"

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
    def __init__(self, filename=None, ax=None):
        self.scores = None
        if filename: self.AddStructureFromFile(filename)
        if ax: self.ax = ax

    def AddStructureFromFile(self, filename):
        x, y, z, _ = ReadPDB(filename)
        self.coordinate = np.array([x, y, z]).T
        self.Center2O()
    
    def AddScores(self, scores):
        num_coords = len(self.coordinate)
        num_scores = len(scores)

        if num_scores != num_coords:
            print(f"Warning: Mismatch between coordinates ({num_coords}) and scores ({num_scores}). Truncating to the shorter length.")
            min_len = min(num_coords, num_scores)
            self.coordinate = self.coordinate[:min_len]
            self.scores = np.array(scores[:min_len])
        else:
            self.scores = np.array(scores)

    def Center2O(self):
        self.Translation(self.coordinate.mean(axis=0))

    def Translation(self, pos):
        self.coordinate -= pos

    def Rotate(self, m):
        self.coordinate = np.dot(self.coordinate, m)

    def DrawLine(self, threshold=0.0):
        # Optimized drawing for compartments with a threshold
        if self.scores is None or len(self.scores) == 0:
            self.ax.plot(self.coordinate[:, 0], self.coordinate[:, 1], self.coordinate[:, 2], color=TRANSITION_COLOR, linewidth=self.line_size)
            return

        coords = self.coordinate
        start_index = 0
        
        def get_color(score):
            if score > threshold: return A_COMP_COLOR
            if score < -threshold: return B_COMP_COLOR
            return TRANSITION_COLOR
            
        current_color = get_color(self.scores[0])
        
        for i in range(1, len(coords)):
            new_color = get_color(self.scores[i])
            if new_color != current_color or i == len(coords) - 1:
                segment = coords[start_index : i + 1]
                self.ax.plot(segment[:, 0], segment[:, 1], segment[:, 2], color=current_color, linewidth=self.line_size)
                start_index = i
                current_color = new_color

    def SetDrawArgs(self, l_size=None):
        self.line_size = float(l_size) if l_size else 2.0

    def Draw(self, threshold=0.0):
        self.DrawLine(threshold=threshold)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize PDB structures with compartment coloring based on a threshold.")
    parser.add_argument("i", help="The PDB file(s). Split with ':'.")
    parser.add_argument("-o", help="Output picture file name.")
    parser.add_argument("-lw", "--linewidth", help="Line width(s). Default: 2.0. Split with ':'.", default="2.0")
    parser.add_argument("-cf", "--compartment_file", help="A/B compartment score file(s). Split with ':'.")
    parser.add_argument("--skip-range", nargs=2, type=int, metavar=('START', 'END'), help="1-based bin range to skip in the score file.")
    
    # MODIFIED: Add threshold argument
    parser.add_argument("-t", "--threshold", type=float, default=0.0, help="Absolute c-score value to define A/B compartments. Bins with scores between [-t, t] are colored gray. Default: 0.0")

    parser.add_argument("-ws", "--withsuperposition", action='store_true', help="Superpose structures.")
    parser.add_argument("--withRMSD", action='store_true', help="Print RMSD values.")
    parser.add_argument("--show", action='store_true', help="Show plot on screen.")
    
    args = parser.parse_args()

    pdbs = args.i.split(':')
    line_widths = args.linewidth.split(':')
    compartment_files = args.compartment_file.split(':') if args.compartment_file else []

    if len(line_widths) < len(pdbs): line_widths.extend([line_widths[-1]] * (len(pdbs) - len(line_widths)))
    
    figure = plt.figure(figsize=(8, 8))
    ax = Axes3D(figure)

    structures = []
    for i, pdb in enumerate(pdbs):
        s = Structure(pdb, ax)
        if i < len(compartment_files):
            scores = ReadCompartmentScores(compartment_files[i])
            if args.skip_range:
                start_idx, end_idx = args.skip_range[0] - 1, args.skip_range[1] - 1
                scores = scores[:start_idx] + scores[end_idx + 1:]
            s.AddScores(scores)
        structures.append(s)

    for i, s in enumerate(structures):
        s.SetDrawArgs(l_size=line_widths[i])

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

    print(f"\nDrawing structure(s) with c-score threshold: {args.threshold}...")
    for s in structures:
        s.Draw(threshold=args.threshold)
    
    plt.axis('off')

    if args.o:
        plt.savefig(args.o, dpi=300, bbox_inches='tight')

    if args.show:
        plt.show()
