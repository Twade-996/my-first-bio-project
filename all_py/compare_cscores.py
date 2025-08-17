#!/usr/bin/env python3
import os
import argparse
from typing import List, Tuple

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# ===================================================================
# 1. Utility Functions
# ===================================================================

def load_cscore(file_path: str) -> pd.DataFrame:
    """Loads C-score data from a bedgraph file."""
    try:
        df = pd.read_csv(
            file_path, 
            sep="\t", 
            header=None, 
            names=["chr", "start", "end", "cscore"],
            comment='t' # Ignore track lines starting with 't'
        )
        df["bin_index"] = range(len(df))
        return df
    except FileNotFoundError:
        print(f"Error: Input file not found -> {file_path}")
        exit(1)
    except pd.errors.EmptyDataError:
        print(f"Error: Input file is empty -> {file_path}")
        exit(1)


def call_compartment_status(cscores: pd.Series, threshold: float = 0) -> np.ndarray:
    """Determines the compartment status (A/B) based on a C-score threshold."""
    return np.where(cscores >= threshold, "A", "B")


def find_boundary_indices(status: np.ndarray) -> List[int]:
    """
    Identifies the indices of boundaries where compartment status changes.
    A boundary is found at index i if status[i] is different from status[i+1].
    """
    status_series = pd.Series(status)
    is_boundary = status_series.shift(-1) != status_series
    # Get indices of all boundaries and remove the last invalid one created by shift()
    boundary_indices = is_boundary[is_boundary].index[:-1].tolist()
    return boundary_indices

# ===================================================================
# 2. Plotting Functions
# ===================================================================

def plot_heatmap(df: pd.DataFrame, name1: str, name2: str, outdir: str):
    """Plots a heatmap of the C-score correlation."""
    plt.figure(figsize=(4, 3))
    corr_matrix = df[[f"cscore_{name1}", f"cscore_{name2}"]].corr()
    sns.heatmap(corr_matrix, annot=True, cmap="vlag", vmin=-1, vmax=1, fmt=".3f")
    plt.title("C-score Correlation")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "cscore_correlation_heatmap.png"), dpi=300)
    plt.close()


def plot_scatter(df: pd.DataFrame, name1: str, name2: str, outdir: str):
    """Plots a scatter plot of C-scores."""
    plt.figure(figsize=(6, 6))
    sns.scatterplot(data=df, x=f"cscore_{name1}", y=f"cscore_{name2}", s=10, alpha=0.4, edgecolor=None)
    plt.xlabel(f"{name1} C-score")
    plt.ylabel(f"{name2} C-score")
    plt.title("C-score Scatter Plot")
    # Add reference lines
    plt.axline((0, 0), slope=1, color='gray', linestyle='--', linewidth=1)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "cscore_scatter.png"), dpi=300)
    plt.close()


def plot_distribution(df: pd.DataFrame, name1: str, name2: str, outdir: str):
    """Plots a kernel density estimate of C-score distributions."""
    plt.figure(figsize=(8, 4))
    sns.kdeplot(df[f"cscore_{name1}"], label=name1, fill=True, alpha=0.5, linewidth=1.5)
    sns.kdeplot(df[f"cscore_{name2}"], label=name2, fill=True, alpha=0.5, linewidth=1.5)
    plt.legend()
    plt.title("C-score Distribution")
    plt.xlabel("C-score")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "cscore_distribution.png"), dpi=300)
    plt.close()


def plot_delta_bar(df: pd.DataFrame, name1: str, name2: str, outdir: str):
    """Plots a bar chart of the ΔC-scores."""
    plt.figure(figsize=(12, 4))
    delta_col = f"delta_cscore_{name1}_vs_{name2}"
    colors = np.where(df[delta_col] > 0, 'red', 'blue')
    plt.bar(df["bin_index"], df[delta_col], color=colors, width=1.0)
    plt.xlabel("Bin index")
    plt.ylabel(f"ΔC-score ({name1} - {name2})")
    plt.title("ΔC-score per Bin")
    plt.axhline(0, color='black', linewidth=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "cscore_delta_bar.png"), dpi=300)
    plt.close()


def plot_switch_rate(switch_rate: float, outdir: str):
    """Plots the compartment switch rate."""
    no_switch_rate = 1 - switch_rate
    labels = ["Switch", "No Switch"]
    proportions = [switch_rate, no_switch_rate]
    
    plt.figure(figsize=(4, 4))
    bars = plt.bar(labels, proportions, color=["red", "grey"])
    plt.ylabel("Proportion")
    plt.title(f"Compartment switch rate: {switch_rate:.2%}")
    plt.ylim(0, 1)
    # Add percentage labels on top of bars
    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2.0, yval, f'{yval:.2%}', va='bottom', ha='center')
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "compartment_switch_rate.png"), dpi=300)
    plt.close()


def plot_boundary_drift(drift: float, outdir: str):
    """Plots the mean boundary drift."""
    plt.figure(figsize=(4, 4))
    plt.bar(["Mean boundary drift"], [drift], color="green")
    plt.ylabel("Mean absolute drift (bins)")
    plt.title("Mean Boundary Drift")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "boundary_drift.png"), dpi=300)
    plt.close()

# ===================================================================
# 3. Main Workflow
# ===================================================================
def main():
    parser = argparse.ArgumentParser(
        description="Compare C-scores between two samples, calculate quantitative metrics, and generate visualizations.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("cscore1", help="C-score bedgraph file for sample 1")
    parser.add_argument("cscore2", help="C-score bedgraph file for sample 2")
    parser.add_argument("--outdir", required=True, help="Output directory to store all results")
    # --- Optional arguments ---
    parser.add_argument("--name1", default="File1", help="Name for sample 1 (used in plots and column headers)")
    parser.add_argument("--name2", default="File2", help="Name for sample 2 (used in plots and column headers)")
    
    args = parser.parse_args()

    print(f"[*] Creating output directory: {args.outdir}")
    os.makedirs(args.outdir, exist_ok=True)

    # --- Load and merge data ---
    print("[*] Loading and merging data...")
    df1 = load_cscore(args.cscore1)
    df2 = load_cscore(args.cscore2)
    # Use more descriptive suffixes
    df = pd.merge(
        df1, df2, 
        on=["chr", "start", "end", "bin_index"], 
        suffixes=(f"_{args.name1}", f"_{args.name2}")
    )

    # --- Calculate core metrics ---
    print("[*] Calculating core metrics...")
    cscore1_col = f"cscore_{args.name1}"
    cscore2_col = f"cscore_{args.name2}"
    status1_col = f"status_{args.name1}"
    status2_col = f"status_{args.name2}"
    delta_col = f"delta_cscore_{args.name1}_vs_{args.name2}"

    # 1. Calculate ΔC-score (defined as name1 - name2)
    df[delta_col] = df[cscore1_col] - df[cscore2_col]

    # 2. Call compartment status
    df[status1_col] = call_compartment_status(df[cscore1_col])
    df[status2_col] = call_compartment_status(df[cscore2_col])

    # 3. Calculate compartment switch rate
    switch_bins = (df[status1_col] != df[status2_col]).sum()
    switch_rate = switch_bins / len(df)

    # 4. Calculate boundary drift
    boundaries1 = find_boundary_indices(df[status1_col].values)
    boundaries2 = find_boundary_indices(df[status2_col].values)
    
    # Only calculate drift if both samples have boundaries
    if len(boundaries1) > 0 and len(boundaries2) > 0:
        # For each boundary in sample 1, find the nearest neighbor in sample 2, 
        # calculate the distance, and then find the mean of these distances.
        drift_distances = [min(abs(b1 - b2) for b2 in boundaries2) for b1 in boundaries1]
        mean_drift = np.mean(drift_distances)
    else:
        mean_drift = np.nan # Drift cannot be calculated if one sample has no boundaries

    print(f"    - Compartment switch rate: {switch_rate:.2%}")
    print(f"    - Mean boundary drift: {mean_drift:.2f} bins")

    # --- Save result tables ---
    print("[*] Saving result tables...")
    # Save detailed bin-by-bin comparison results
    df.to_csv(os.path.join(args.outdir, "cscore_comparison_detailed.tsv"), sep="\t", index=False)
    
    # Save summary statistics
    summary = pd.DataFrame({
        "switch_rate": [switch_rate],
        "mean_boundary_drift_bins": [mean_drift],
        "mean_abs_delta_cscore": [df[delta_col].abs().mean()],
        "pearson_correlation": [df[cscore1_col].corr(df[cscore2_col], method='pearson')]
    })
    summary.to_csv(os.path.join(args.outdir, "cscore_comparison_summary.csv"), index=False)

    # --- Generate all plots ---
    print("[*] Generating all plots...")
    plot_heatmap(df, args.name1, args.name2, args.outdir)
    plot_scatter(df, args.name1, args.name2, args.outdir)
    plot_distribution(df, args.name1, args.name2, args.outdir)
    plot_delta_bar(df, args.name1, args.name2, args.outdir)
    plot_switch_rate(switch_rate, args.outoutdir)
    plot_boundary_drift(mean_drift, args.outdir)

    print(f"\n[✓] Analysis complete! All results have been saved to: {args.outdir}")

if __name__ == "__main__":
    main()
