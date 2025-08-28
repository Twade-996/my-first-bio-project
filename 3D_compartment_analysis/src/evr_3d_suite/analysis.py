# 文件名: src/evr_3d_suite/analysis.py
# 描述: 整合Hi-C区室状态与RNA-seq基因表达数据进行功能分析

import argparse
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, wilcoxon

# 从同一个包中导入IO工具
from . import io_utils

# ===================================================================
# 1. 核心逻辑与命令行入口点函数
# ===================================================================

def main_analyze_expression_change():
    """
    命令行工具: evr-analyze-change
    关联区室转换与基因表达变化 (log2FoldChange)。
    """
    parser = argparse.ArgumentParser(
        description="关联区室转换与基因表达变化 (log2FoldChange)。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('--bin-gene-map', required=True, help="bedtools生成的bin-基因映射文件路径。")
    parser.add_argument('--rna-seq', required=True, help="差异RNA-seq结果文件路径 (.csv或.tsv)。")
    parser.add_argument('--scores1', required=True, help="条件1 (基线) 的c-score文件路径。")
    parser.add_argument('--scores2', required=True, help="条件2 (比较对象) 的c-score文件路径。")
    parser.add_argument('--name1', default="Cond1", help="条件1的名称 (例如, GM12878)。")
    parser.add_argument('--name2', default="Cond2", help="条件2的名称 (例如, K562)。")
    parser.add_argument('-t', '--threshold', type=float, default=0.1, help="定义A/B区室的c-score绝对值阈值。")
    parser.add_argument('-o', '--output', default="compartment_expression_change.png", help="输出图表文件名。")
    args = parser.parse_args()

    print("[*] 正在加载数据...")
    bin_to_gene, _ = io_utils.parse_bin_gene_map(args.bin_gene_map)
    de_map = io_utils.load_de_results(args.rna_seq)
    scores1 = io_utils.read_compartment_scores(args.scores1)
    scores2 = io_utils.read_compartment_scores(args.scores2)

    if not all((bin_to_gene, de_map, scores1, scores2)):
        print("错误: 一个或多个输入文件加载失败或为空。正在中止。", file=sys.stderr)
        sys.exit(1)

    print("[*] 正在分类bins并关联基因表达...")
    # 补全: 初始化字典
    categories = ['Stable A', 'Stable B', 'A -> B', 'B -> A', 'Other/Transition']
    results = {cat: [] for cat in categories}
    min_len = min(len(scores1), len(scores2))

    for i in range(min_len):
        s1, s2 = scores1[i], scores2[i]
        state1 = 'A' if s1 > args.threshold else ('B' if s1 < -args.threshold else 'T')
        state2 = 'A' if s2 > args.threshold else ('B' if s2 < -args.threshold else 'T')
        
        category = 'Other/Transition'
        if state1 == 'A' and state2 == 'A': category = 'Stable A'
        elif state1 == 'B' and state2 == 'B': category = 'Stable B'
        elif state1 == 'A' and state2 == 'B': category = 'A -> B'
        elif state1 == 'B' and state2 == 'A': category = 'B -> A'
        
        if i in bin_to_gene:
            for gene in bin_to_gene[i]:
                if gene in de_map:
                    results[category].append(de_map[gene])

    # 补全: 过滤掉没有数据的类别，并保持一个固定的绘图顺序
    results = {k: v for k, v in results.items() if v}
    order = [cat for cat in ['Stable B', 'A -> B', 'Other/Transition', 'Stable A', 'B -> A'] if cat in results]

    print("[*] 正在生成图表和统计数据...")
    plot_data = [{'Category': cat, 'log2FoldChange': fc} for cat, fc_list in results.items() for fc in fc_list]
    if not plot_data:
        print("错误: 没有基因表达数据可以匹配。分析无法继续。", file=sys.stderr)
        sys.exit(1)
    plot_df = pd.DataFrame(plot_data)

    plt.figure(figsize=(12, 8))
    sns.boxplot(data=plot_df, x='Category', y='log2FoldChange', order=order,
                palette={'Stable B':'lightblue', 'A -> B':'blue', 'Stable A':'lightcoral', 
                         'B -> A':'red', 'Other/Transition':'lightgray'})
    plt.axhline(0, color='black', linestyle='--')
    title = f'Gene Expression Changes vs. Compartment Switches ({args.name1} to {args.name2})'
    plt.title(title, fontsize=16); plt.ylabel(f'log2 Fold Change ({args.name2} / {args.name1})', fontsize=12); plt.xlabel('Compartment Category', fontsize=12)
    
    # --- 统计分析 ---
    print("\n--- 统计检验结果 ---")
    stats_left, stats_right = [], []
    text_bg = dict(facecolor='white', alpha=0.7, edgecolor='none', pad=2)

    print("组内Wilcoxon符号秩检验 (中位数 vs 0):")
    # 补全: 迭代变化的类别
    for cat in ['A -> B', 'B -> A']:
        if cat in results and len(results[cat]) > 1:
            # 修正: 移除0值以进行符号秩检验
            non_zero_vals = [x for x in results[cat] if x != 0]
            if not non_zero_vals: continue
            _, p_val = wilcoxon(non_zero_vals, alternative='two-sided')
            stats_left.append(f"{cat} vs 0: p={p_val:.2e}")
            print(f"  {cat} vs 0: p-value = {p_val:.2e}")

    print("组间Mann-Whitney U检验 (vs 稳定状态):")
    if 'A -> B' in results and 'Stable B' in results:
        # 修正: 比较正确的组
        _, p_val_ab = mannwhitneyu(results['A -> B'], results['Stable B'], alternative='two-sided')
        stats_right.append(f"A->B vs Stable B: p={p_val_ab:.2e}")
        print(f"  A->B vs Stable B: p-value = {p_val_ab:.2e}")
    if 'B -> A' in results and 'Stable A' in results:
        # 修正: 比较正确的组
        _, p_val_ba = mannwhitneyu(results['B -> A'], results['Stable A'], alternative='two-sided')
        stats_right.append(f"B->A vs Stable A: p={p_val_ba:.2e}")
        print(f"  B->A vs Stable A: p-value = {p_val_ba:.2e}")

    plt.figtext(0.13, 0.01, "One-sample (Wilcoxon Signed-Rank):\n" + "\n".join(stats_left), ha="left", fontsize=9, va="bottom", bbox=text_bg)
    plt.figtext(0.87, 0.01, "Two-sample (Mann-Whitney U):\n" + "\n".join(stats_right), ha="right", fontsize=9, va="bottom", bbox=text_bg)

    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    print(f"\n[✓] 分析完成! 图表已保存至 '{args.output}'")


def main_analyze_absolute_expression():
    """
    命令行工具: evr-analyze-absolute
    分析绝对表达水平与区室转换的关系。
    """
    parser = argparse.ArgumentParser(description="分析绝对表达水平与区室转换的关系。")
    parser.add_argument('--rsem-files', nargs='+', required=True, help="目标细胞系的RSEM.tsv文件列表。")
    parser.add_argument('--gtf', required=True, help="用于ID转换的GENCODE注释GTF文件路径。")
    parser.add_argument('--bin-gene-map', required=True, help="bedtools生成的bin-基因映射文件路径。")
    parser.add_argument('--scores1', required=True, help="起始条件的c-score文件路径。")
    parser.add_argument('--scores2', required=True, help="目标条件的c-score文件路径。")
    parser.add_argument('--name2', default="TargetCond", help="目标条件的名称。")
    parser.add_argument('-t', '--threshold', type=float, default=0.1, help="定义A/B区室的c-score绝对值阈值。")
    parser.add_argument('-o', '--output', default="absolute_expression_levels.png", help="输出图表文件名。")
    args = parser.parse_args()

    print("[*] 正在加载数据...")
    id_map = io_utils.create_id_map_from_gtf(args.gtf)
    tpm_map = io_utils.process_rsem_to_tpm_map(args.rsem_files, id_map)
    bin_to_gene, _ = io_utils.parse_bin_gene_map(args.bin_gene_map)
    scores1 = io_utils.read_compartment_scores(args.scores1)
    scores2 = io_utils.read_compartment_scores(args.scores2)

    if not all((id_map, tpm_map, bin_to_gene, scores1, scores2)):
        print("错误: 一个或多个输入文件加载失败或为空。正在中止。", file=sys.stderr); sys.exit(1)

    print("[*] 正在关联表达量与区室状态...")
    # 补全: 初始化字典
    results = {'B -> A': [], 'Stable A': []}
    min_len = min(len(scores1), len(scores2))

    for i in range(min_len):
        s1, s2 = scores1[i], scores2[i]
        state1 = 'A' if s1 > args.threshold else ('B' if s1 < -args.threshold else 'T')
        state2 = 'A' if s2 > args.threshold else ('B' if s2 < -args.threshold else 'T')
        
        category = None
        if state1 == 'B' and state2 == 'A': category = 'B -> A'
        elif state1 == 'A' and state2 == 'A': category = 'Stable A'
        
        if category and i in bin_to_gene:
            for gene in bin_to_gene[i]:
                if gene in tpm_map:
                    results[category].append(np.log2(tpm_map[gene] + 1))

    print("[*] 正在生成图表和统计数据...")
    # 补全: 创建DataFrame
    plot_data = [{'Category': cat, 'log2(TPM+1)': val} for cat, val_list in results.items() for val in val_list]
    if not plot_data:
        print("错误: 没有数据可供绘图。", file=sys.stderr); sys.exit(1)
    plot_df = pd.DataFrame(plot_data)
    
    plt.figure(figsize=(6, 8))
    # 补全: 指定绘图顺序
    sns.boxplot(data=plot_df, x='Category', y='log2(TPM+1)', order=['Stable A', 'B -> A'],
                palette={'Stable A': 'lightcoral', 'B -> A': 'red'})
    plt.title(f'Absolute Expression Levels in {args.name2}', fontsize=16)
    plt.ylabel(f'Expression Level (log2(TPM+1)) in {args.name2}', fontsize=12)
    plt.xlabel('Compartment Category', fontsize=12)
    
    try:
        # 修正: 比较正确的组
        u_stat, p_value = mannwhitneyu(results['B -> A'], results['Stable A'], alternative='two-sided')
        stats_text = f"Mann-Whitney U Test:\np-value = {p_value:.2e}"
        print(f"\n统计检验结果:\n{stats_text}")
        plt.text(0.5, 0.95, stats_text, ha='center', va='center', transform=plt.gca().transAxes, fontsize=12,
                 bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
    except ValueError as e:
        print(f"\n无法执行统计检验: {e}", file=sys.stderr)

    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    print(f"\n[✓] 分析完成! 图表已保存至 '{args.output}'")