#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr

def read_bedgraph(file):
    df = pd.read_csv(file, sep='\t', comment='t', header=None)
    df.columns = ['chr','start','end','score']
    df['bin'] = df['chr'] + ':' + df['start'].astype(str) + '-' + df['end'].astype(str)
    return df[['bin','score']]

def merge_bedgraph(file1, file2):
    df1 = read_bedgraph(file1)
    df2 = read_bedgraph(file2)
    merged = pd.merge(df1, df2, on='bin', how='outer', suffixes=('1','2'))
    return merged.fillna(np.nan)

def calc_correlation(df, method='pearson'):
    valid = df.dropna()
    if method=='pearson':
        r, p = pearsonr(valid['score1'], valid['score2'])
    elif method=='spearman':
        r, p = spearmanr(valid['score1'], valid['score2'])
    else:
        raise ValueError("Method must be 'pearson' or 'spearman'")
    return r, p

def plot_scatter(df, outdir):
    plt.figure(figsize=(6,6))
    plt.scatter(df['score1'], df['score2'], alpha=0.5, s=10)
    plt.xlabel('File1 C-score')
    plt.ylabel('File2 C-score')
    plt.title('C-score Scatter Plot')
    plt.tight_layout()
    plt.savefig(os.path.join(outdir,'cscore_scatter.png'))
    plt.close()

def plot_distribution(df, outdir):
    plt.figure(figsize=(6,4))
    sns.kdeplot(df['score1'].dropna(), label='File1', fill=True)
    sns.kdeplot(df['score2'].dropna(), label='File2', fill=True)
    plt.xlabel('C-score')
    plt.ylabel('Density')
    plt.title('C-score Distribution')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir,'cscore_distribution.png'))
    plt.close()

def plot_heatmap(df, outdir):
    plt.figure(figsize=(8,6))
    sns.heatmap(df[['score1','score2']].transpose(), cmap='coolwarm', cbar_kws={'label':'C-score'})
    plt.xlabel('Bin index')
    plt.ylabel('File')
    plt.title('C-score Heatmap')
    plt.tight_layout()
    plt.savefig(os.path.join(outdir,'cscore_heatmap.png'))
    plt.close()

def plot_delta_cscore(df, outdir):
    df['delta'] = df['score1'] - df['score2']
    plt.figure(figsize=(12,4))
    colors = ['red' if v>=0 else 'blue' for v in df['delta']]
    plt.bar(range(len(df)), df['delta'], color=colors, width=1.0)
    plt.xlabel('Bin index')
    plt.ylabel('ΔC-score (File1 - File2)')
    plt.title('ΔC-score per bin')
    plt.tight_layout()
    plt.savefig(os.path.join(outdir,'cscore_delta_bar.png'))
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Compare two C-score bedGraph files')
    parser.add_argument('file1')
    parser.add_argument('file2')
    parser.add_argument('--outdir', default='cscore_compare_output')
    parser.add_argument('--analysis', nargs='+', choices=['pearson','spearman','scatter','distribution','heatmap','delta'], 
                        help='Which analyses to output')
    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    merged_df = merge_bedgraph(args.file1, args.file2)

    # 默认输出全部
    analysis_to_run = args.analysis if args.analysis else ['pearson','spearman','scatter','distribution','heatmap','delta']

    # 相关性
    if 'pearson' in analysis_to_run:
        r, p = calc_correlation(merged_df, method='pearson')
        print(f"Pearson r={r:.4f}, p={p:.4e}")
        with open(os.path.join(args.outdir,'pearson.txt'),'w') as f:
            f.write(f"Pearson r={r:.4f}, p={p:.4e}\n")
    if 'spearman' in analysis_to_run:
        r, p = calc_correlation(merged_df, method='spearman')
        print(f"Spearman rho={r:.4f}, p={p:.4e}")
        with open(os.path.join(args.outdir,'spearman.txt'),'w') as f:
            f.write(f"Spearman rho={r:.4f}, p={p:.4e}\n")

    # 可视化
    if 'scatter' in analysis_to_run:
        plot_scatter(merged_df, args.outdir)
    if 'distribution' in analysis_to_run:
        plot_distribution(merged_df, args.outdir)
    if 'heatmap' in analysis_to_run:
        plot_heatmap(merged_df, args.outdir)
    if 'delta' in analysis_to_run:
        plot_delta_cscore(merged_df, args.outdir)

    # 保存合并表
    merged_df.to_csv(os.path.join(args.outdir,'cscore_compare.txt'), sep='\t', index=False)

if __name__ == '__main__':
    main()
