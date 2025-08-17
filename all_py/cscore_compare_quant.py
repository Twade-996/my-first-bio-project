#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
C-score comparison automation: delta C-score, compartment switch rate, boundary drift.
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def read_cscore_bedgraph(file):
    """读取 C-score bedgraph 文件，返回 DataFrame: ['chr', 'start', 'end', 'cscore']"""
    df = pd.read_csv(file, sep='\t', comment='t', header=None, names=['chr','start','end','cscore'])
    return df

def classify_A_B(df, threshold=0.0):
    """根据 C-score 分类 A/B 区室"""
    return df['cscore'].apply(lambda x: 'A' if x > threshold else 'B')

def delta_cscore(df1, df2):
    """计算 ΔC-score"""
    return df1['cscore'] - df2['cscore']

def compartment_switch_rate(status1, status2):
    """计算区室交换率"""
    switch = status1 != status2
    return switch.sum() / len(switch)

def boundary_positions(status):
    """获取 bin 边界位置（A→B 或 B→A），返回 bin 的位置索引"""
    arr = status.to_numpy()  # 转成 numpy array，保证位置对齐
    boundaries = np.where(arr[:-1] != arr[1:])[0].tolist()
    return boundaries

def boundary_drift(bound1, bound2):
    """计算平均绝对边界漂移，按最小长度对齐"""
    n = min(len(bound1), len(bound2))
    if n == 0:
        return np.nan
    return np.mean(np.abs(np.array(bound1[:n]) - np.array(bound2[:n])))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("cscore1", help="C-score file for state 1")
    parser.add_argument("cscore2", help="C-score file for state 2")
    parser.add_argument("--outdir", help="Output folder", default="cscore_compare_out")
    parser.add_argument("--threshold", type=float, default=0.0, help="Threshold for A/B compartment classification")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df1 = read_cscore_bedgraph(args.cscore1)
    df2 = read_cscore_bedgraph(args.cscore2)

    if len(df1) != len(df2):
        # 合并对齐, 缺失补 NA
        df = pd.merge(df1, df2, on=['chr','start','end'], how='outer', suffixes=('_1','_2'))
    else:
        df = pd.DataFrame({
            'chr': df1['chr'],
            'start': df1['start'],
            'end': df1['end'],
            'cscore_1': df1['cscore'],
            'cscore_2': df2['cscore']
        })

    # 分类 A/B
    status1 = df['cscore_1'].apply(lambda x: 'A' if x > args.threshold else 'B')
    status2 = df['cscore_2'].apply(lambda x: 'A' if x > args.threshold else 'B')

    # ΔC-score
    df['delta_cscore'] = df['cscore_1'] - df['cscore_2']

    # 区室交换率
    switch_rate = compartment_switch_rate(status1, status2)

    # 边界漂移
    boundaries1 = boundary_positions(status1)
    boundaries2 = boundary_positions(status2)
    drift = boundary_drift(boundaries1, boundaries2)

    # 输出 CSV
    df_outfile = os.path.join(args.outdir, "cscore_comparison_summary.csv")
    df.to_csv(df_outfile, sep='\t', index=False)

    # 输出区室交换率
    plt.figure(figsize=(4,4))
    plt.bar(['Switch', 'No Switch'], [switch_rate, 1-switch_rate], color=['red','grey'])
    plt.ylabel("Proportion")
    plt.title("Compartment switch rate: {:.2%}".format(switch_rate))
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir,"compartment_switch_rate.png"))
    plt.close()

    # 输出边界漂移
    plt.figure(figsize=(4,4))
    plt.bar(['Mean boundary drift'], [drift], color='green')
    plt.ylabel("Bins")
    plt.title("Mean boundary drift")
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir,"boundary_drift.png"))
    plt.close()

    # 打印汇总
    print(f"Summary saved to {df_outfile}")
    print(f"Compartment switch rate: {switch_rate:.2%}")
    print(f"Mean boundary drift: {drift:.2f} bins")

if __name__ == "__main__":
    main()
