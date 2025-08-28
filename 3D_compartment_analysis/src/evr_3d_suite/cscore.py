# 文件名: src/evr_3d_suite/cscore.py
# 描述: C-score分析、比较和2D可视化

import argparse
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

# 从同一个包中导入IO工具
from . import io_utils

# --- 核心功能函数 ---

def analyze_thresholds_logic(scores_raw, skip_range=None):
    """
    分析C-score分布并返回统计数据。
    """
    if not scores_raw:
        print("警告: 未加载任何分数。文件可能为空或格式错误。", file=sys.stderr)
        return {}, np.array([])

    if skip_range:
        # 修正: 正确地从元组中获取索引
        start_idx = skip_range[0] - 1
        end_idx = skip_range[1] - 1
        scores = scores_raw[:start_idx] + scores_raw[end_idx + 1:]
    else:
        scores = scores_raw
    
    scores = np.array(scores)
    if scores.size == 0:
        print("警告: 过滤后所有分数都被移除。请检查您的 --skip-range。", file=sys.stderr)
        return {}, np.array([])

    stats = {
        "count": len(scores),
        "mean": np.mean(scores),
        "median": np.median(scores), # 使用np.median更直接
        "std": np.std(scores),
        "q1": np.percentile(scores, 25),
        "q3": np.percentile(scores, 75)
    }
    return stats, scores

def compare_scores_logic(df1, df2, name1, name2, threshold=0.0):
    """
    比较两个C-score数据集并计算核心指标。
    """
    df = pd.merge(
        df1, df2, 
        on=["chr", "start", "end"], 
        suffixes=(f"_{name1}", f"_{name2}")
    )
    if df.empty:
        print("警告: 两个C-score文件没有共同的bins，无法比较。", file=sys.stderr)
        return pd.DataFrame(), {}

    cscore1_col, cscore2_col = f"cscore_{name1}", f"cscore_{name2}"
    status1_col, status2_col = f"status_{name1}", f"status_{name2}"
    delta_col = f"delta_cscore"

    df[delta_col] = df[cscore2_col] - df[cscore1_col]
    
    # 改进: 使用可配置的阈值
    df[status1_col] = np.select([df[cscore1_col] > threshold, df[cscore1_col] < -threshold], ['A', 'B'], default='T')
    df[status2_col] = np.select([df[cscore2_col] > threshold, df[cscore2_col] < -threshold], ['A', 'B'], default='T')

    switch_bins = (df[status1_col]!= df[status2_col]).sum()
    switch_rate = switch_bins / len(df)

    # 边界漂移计算
    boundaries1 = df.index[df[status1_col] != df[status1_col].shift()].tolist()[1:]
    boundaries2 = df.index[df[status2_col] != df[status2_col].shift()].tolist()[1:]
    
    mean_drift = np.nan
    if boundaries1 and boundaries2: # 健壮性检查
        drift_distances = [min(abs(b1 - b2) for b2 in boundaries2) for b1 in boundaries1]
        if drift_distances: # 健壮性检查
            mean_drift = np.mean(drift_distances)

    summary = {
        "total_bins_compared": len(df),
        "switch_rate": switch_rate,
        "mean_boundary_drift_bins": mean_drift,
        "mean_abs_delta_cscore": df[delta_col].abs().mean(),
        "pearson_correlation": df[cscore1_col].corr(df[cscore2_col], method='pearson')
    }
    
    return df, summary

# --- 命令行入口点函数 ---

def main_threshold():
    """命令行工具: evr-cscore-threshold"""
    parser = argparse.ArgumentParser(description="分析C-score文件以帮助确定A/B区室阈值。")
    parser.add_argument("cscore_file", help="要分析的C-score文件路径。")
    parser.add_argument("--skip-range", nargs=2, type=int, metavar=('START', 'END'), help="可选：要跳过的基于1的bin范围。")
    parser.add_argument('--no-plot', action='store_true', help="可选：不显示直方图，只显示文本统计信息。")
    args = parser.parse_args()

    scores_raw = io_utils.read_compartment_scores(args.cscore_file)
    stats, scores = analyze_thresholds_logic(scores_raw, args.skip_range)

    if stats:
        # 改进: 格式化输出
        print("\n" + "="*40)
        print("          C-Score 统计分析          ")
        print("="*40)
        print(f"{'分析的总bins数':<20}: {stats['count']}")
        print(f"{'均值':<20}: {stats['mean']:.4f}")
        print(f"{'中位数':<20}: {stats['median']:.4f}")
        print(f"{'标准差':<20}: {stats['std']:.4f}")
        print(f"{'25百分位数 (Q1)':<20}: {stats['q1']:.4f}")
        print(f"{'75百分位数 (Q3)':<20}: {stats['q3']:.4f}")
        print("\n" + "-"*40)
        print("        建议的阈值        ")
        print("-"*40)
        print(f"{'基于 0.5 * 标准差':<20}: {0.5 * stats['std']:.4f}")
        print(f"{'基于 1.0 * 标准差':<20}: {stats['std']:.4f}")
        print(f"{'基于四分位数 (Q1/Q3)':<20}: 考虑 {abs(stats['q1']):.4f} - {abs(stats['q3']):.4f}")
        print("="*40 + "\n")

    if not args.no_plot and scores.size > 0:
        plt.figure(figsize=(10, 6))
        sns.histplot(scores, bins=100, color='skyblue', kde=True)
        plt.title(f'C-Score分布\n{os.path.basename(args.cscore_file)}', fontsize=14)
        plt.xlabel('C-Score', fontsize=12)
        plt.ylabel('Bin数量', fontsize=12)
        plt.grid(True, linestyle='--', linewidth=0.5)
        plt.show()

def main_compare():
    """命令行工具: evr-cscore-compare"""
    parser = argparse.ArgumentParser(description="比较两个样本的C-score，计算量化指标并生成可视化图表。")
    parser.add_argument("cscore1", help="样本1的C-score bedgraph文件")
    parser.add_argument("cscore2", help="样本2的C-score bedgraph文件")
    parser.add_argument("--outdir", required=True, help="存储所有结果的输出目录")
    parser.add_argument("--name1", default="Sample1", help="样本1的名称")
    parser.add_argument("--name2", default="Sample2", help="样本2的名称")
    parser.add_argument("-t", "--threshold", type=float, default=0.0, help="定义A/B区室的绝对C-score阈值。")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    
    try:
        df1 = pd.read_csv(args.cscore1, sep="\t", header=None, names=["chr", "start", "end", "cscore"], comment='t')
        df2 = pd.read_csv(args.cscore2, sep="\t", header=None, names=["chr", "start", "end", "cscore"], comment='t')
    except Exception as e:
        print(f"错误: 无法读取输入文件。请确保它们是有效的bedgraph格式。错误: {e}", file=sys.stderr)
        sys.exit(1)

    df_results, summary_stats = compare_scores_logic(df1, df2, args.name1, args.name2, args.threshold)

    if not df_results.empty:
        df_results.to_csv(os.path.join(args.outdir, "cscore_comparison_detailed.tsv"), sep="\t", index=False)
        pd.DataFrame([summary_stats]).to_csv(os.path.join(args.outdir, "cscore_comparison_summary.csv"), index=False)
        
        # 改进: 打印摘要到终端
        print("\n" + "="*40)
        print(f"     比较摘要: {args.name1} vs {args.name2}     ")
        print("="*40)
        for key, value in summary_stats.items():
            print(f"{key.replace('_', ' ').title():<25}: {value:.4f}")
        print("="*40 + "\n")

        print("分析完成，正在生成图表...")
        # 示例：生成一个散点图
        plt.figure(figsize=(8, 8))
        sns.scatterplot(data=df_results, x=f"cscore_{args.name1}", y=f"cscore_{args.name2}", alpha=0.5)
        plt.title(f'C-Score Correlation ({args.name1} vs {args.name2})', fontsize=14)
        plt.xlabel(f'C-Score {args.name1}', fontsize=12)
        plt.ylabel(f'C-Score {args.name2}', fontsize=12)
        plt.grid(True, linestyle='--', linewidth=0.5)
        plt.axhline(0, color='black', lw=0.5); plt.axvline(0, color='black', lw=0.5)
        plt.savefig(os.path.join(args.outdir, "cscore_correlation_scatter.png"), dpi=300)
        plt.close() # 关闭图像，防止在循环中使用时弹出
        
        print(f"所有结果已保存至: {args.outdir}")

def main_plot_profile():
    """命令行工具: evr-cscore-plot"""
    parser = argparse.ArgumentParser(description="从bedgraph文件绘制C-score曲线图。")
    parser.add_argument("input_file", help="输入的bedgraph文件路径。")
    parser.add_argument("-o", "--output", help="保存输出图像文件的路径。")
    parser.add_argument("-t", "--title", help="图表的自定义标题。")
    parser.add_argument("--show", action="store_true", help="交互式显示图表而不是仅保存。")
    args = parser.parse_args()
    
    try:
        df = pd.read_csv(args.input_file, sep=r"\s+", header=None, comment="t", engine="python", names=["chr", "start", "end", "cscore"])
    except Exception as e:
        print(f"错误: 无法读取输入文件 '{args.input_file}'。错误: {e}", file=sys.stderr)
        sys.exit(1)
        
    plt.figure(figsize=(15, 5))
    plt.plot(df["start"], df["cscore"], color="#3477eb", linewidth=0.9)
    plt.axhline(0, ls="--", color="gray", zorder=0)
    plt.grid(True, linestyle=':', linewidth=0.6, alpha=0.7)
    plt.title(args.title or f"C-score profile for {os.path.basename(args.input_file)}")
    plt.xlabel("基因组位置 (bp)")
    plt.ylabel("C-score")
    plt.tight_layout()

    if args.output:
        plt.savefig(args.output, dpi=300)
    if args.show:
        plt.show()
    plt.close()