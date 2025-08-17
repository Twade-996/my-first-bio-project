# -*- coding: utf-8 -*-
"""
This module is for file I/O operations, for reading PDB and score files.
"""

def ReadPDB(filename):
    """
    Parses a PDB file to extract coordinates and establish sequential connectivity.
    It reads lines starting with "ATOM".
    """
    x, y, z = [], [], []
    try:
        with open(filename, 'r') as pdb_file:
            for line in pdb_file:
                if line.startswith("ATOM"):
                    x.append(float(line[30:38].strip()))
                    y.append(float(line[38:46].strip()))
                    z.append(float(line[46:54].strip()))
    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found.")
        exit(1)
    except (ValueError, IndexError):
        print(f"Error: Could not parse coordinates from '{filename}'. Check PDB format.")
        exit(1)

    c = {i: i + 1 for i in range(len(x) - 1)}
    return x, y, z, c

def ReadCompartmentScores(filename):
    """
    Reads a two-column file (bin_index, score) and returns a list of scores.
    It assumes the scores are in the second column.
    """
    scores = []
    try:
        with open(filename, 'r') as score_file:
            for line in score_file:
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        scores.append(float(parts[1]))
                    except ValueError:
                        print(f"Warning: Could not parse score from line: '{line.strip()}' in {filename}")
    except FileNotFoundError:
        print(f"Error: The compartment file '{filename}' was not found.")
        exit(1)
    return scores
