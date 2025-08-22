#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""
A command-line tool to analyze the statistical distribution of a c-score file
and suggest appropriate thresholds for A/B compartment definition.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt

# We assume IO.py is in the same directory and contains ReadCompartmentScores
try:
    from IO import ReadCompartmentScores
except ImportError:
    print("Error: Could not import ReadCompartmentScores from IO.py.")
    print("Please ensure IO.py is in the same directory as this script.")
    exit(1)

def main():
    # --- 1. Set up Argument Parser ---
    parser = argparse.ArgumentParser(
        description="Analyze a c-score file to help determine A/B compartment thresholds.",
        formatter_class=argparse.RawTextHelpFormatter # For better help text formatting
    )
    parser.add_argument(
        "cscore_file", 
        help="Path to the c-score file to be analyzed."
    )
    parser.add_argument(
        "--skip-range", 
        nargs=2, 
        type=int, 
        metavar=('START', 'END'),
        help="Optional: 1-based start and end bin range to skip (e.g., for a centromere)."
    )
    parser.add_argument(
        '--no-plot',
        action='store_true',
        help="Optional: Suppress the histogram plot and only show text statistics."
    )
    
    args = parser.parse_args()

    # --- 2. Load and Filter Data ---
    print(f"\nLoading scores from: {args.cscore_file}")
    scores_raw = ReadCompartmentScores(args.cscore_file)
    
    if not scores_raw:
        print("Error: No scores were loaded. The file might be empty or in the wrong format.")
        exit(1)

    if args.skip_range:
        start_idx = args.skip_range[0] - 1
        end_idx = args.skip_range[1] - 1 # This is inclusive in the file, exclusive in slicing
        print(f"Skipping bin range {args.skip_range[0]}-{args.skip_range[1]} (indices {start_idx}-{end_idx})")
        scores = scores_raw[:start_idx] + scores_raw[end_idx + 1:]
    else:
        scores = scores_raw
        
    scores = np.array(scores)

    if scores.size == 0:
        print("Error: All scores were removed after filtering. Check your --skip-range.")
        exit(1)

    # --- 3. Calculate and Print Statistics ---
    mean_val = np.mean(scores)
    std_val = np.std(scores)
    percentile_25 = np.percentile(scores, 25)
    percentile_50 = np.percentile(scores, 50) # Median
    percentile_75 = np.percentile(scores, 75)
    
    print("\n" + "="*40)
    print("          STATISTICAL ANALYSIS          ")
    print("="*40)
    print(f"Total Bins Analyzed (after filtering): {len(scores)}")
    print(f"Mean:                {mean_val:.4f}")
    print(f"Median (50th percentile): {percentile_50:.4f}")
    print(f"Standard Deviation:  {std_val:.4f}")
    print(f"25th Percentile (Q1): {percentile_25:.4f}")
    print(f"75th Percentile (Q3): {percentile_75:.4f}")
    
    print("\n" + "-"*40)
    print("        SUGGESTED THRESHOLDS        ")
    print("-"*40)
    print("These values can be used with the '-t' flag in the visualization scripts.")
    print(f"Based on 0.5 * Std Dev:   {0.5 * std_val:.4f}")
    print(f"Based on 1.0 * Std Dev:   {std_val:.4f}")
    print(f"Based on Quartiles (Q1/Q3): Consider a value between {abs(percentile_25):.4f} and {abs(percentile_75):.4f}")
    print("="*40 + "\n")

    # --- 4. Plot Histogram (unless suppressed) ---
    if not args.no_plot:
        plt.figure(figsize=(10, 6))
        plt.hist(scores, bins=100, color='skyblue', edgecolor='black', alpha=0.7)
        
        # Add vertical lines for statistical markers
        plt.axvline(0.5 * std_val, color='red', linestyle='--', label=f'0.5 Std Dev ({0.5*std_val:.2f})')
        plt.axvline(-0.5 * std_val, color='red', linestyle='--')
        
        plt.axvline(percentile_75, color='green', linestyle=':', label=f'75th Percentile ({percentile_75:.2f})')
        plt.axvline(percentile_25, color='green', linestyle=':', label=f'25th Percentile ({percentile_25:.2f})')
        
        plt.title(f'C-Score Distribution for\n{args.cscore_file.split("/")[-1]}', fontsize=14)
        plt.xlabel('C-Score (PC1 Value)', fontsize=12)
        plt.ylabel('Number of Bins', fontsize=12)
        plt.legend()
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)
        
        print("Displaying histogram plot. Close the plot window to exit the script.")
        plt.show()

if __name__ == "__main__":
    main()
