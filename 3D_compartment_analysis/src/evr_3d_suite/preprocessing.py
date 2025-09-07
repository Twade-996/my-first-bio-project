# 文件名: src/evr_3d_suite/preprocessing.py
# 描述: Hi-C数据预处理和格式转换工具集

import argparse
import sys
import subprocess
import os
import tempfile
import numpy as np
import pandas as pd

# ===================================================================
# 1. 核心逻辑函数 - 这些函数可以被命令行工具调用
# ===================================================================

def _run_juicer_dump(hic_file, juicer_jar, chrom, binsize, norm, output_file):
    """使用juicer_tools从.hic文件导出稀疏矩阵。这是一个核心辅助函数。"""
    print(f"[*] 正在调用Juicer: {os.path.basename(hic_file)} -> {os.path.basename(output_file)}")
    print(f"    (染色体: {chrom}, 分辨率: {binsize} bp, 归一化: {norm})")
    
    command = ['java', '-jar', juicer_jar, 'dump', 'observed', norm, hic_file, chrom, chrom, 'BP', str(binsize), output_file]
    
    try:
        process = subprocess.run(command, check=True, capture_output=True, text=True)
        # 打印Juicer的输出，对于调试很有用
        if process.stdout: print(f"    Juicer STDOUT:\n---\n{process.stdout.strip()}\n---")
        if process.stderr: print(f"    Juicer STDERR:\n---\n{process.stderr.strip()}\n---")
        return True
    except FileNotFoundError:
        print("\n错误: 'java' 命令未找到。请确保Java已安装并配置在您的系统PATH中。", file=sys.stderr)
        return False
    except subprocess.CalledProcessError as e:
        print("\n错误: Juicer命令执行失败。", file=sys.stderr)
        print(f"  返回码: {e.returncode}\n  标准错误:\n{e.stderr}", file=sys.stderr)
        print("  请检查您的参数，特别是染色体名称是否与.hic文件中的一致 (例如 '1' vs 'chr1')。", file=sys.stderr)
        return False
    return False

def _convert_sparse_to_dense(input_file, output_file):
    """将三列稀疏矩阵转换为EVR所需的密集矩阵。"""
    print(f"[*] 正在转换稀疏矩阵为密集矩阵: {os.path.basename(input_file)} -> {os.path.basename(output_file)}")
    try:
        df = pd.read_csv(input_file, sep=r"\s+", header=None, names=["start1", "start2", "value"], dtype={'start1': int, 'start2': int, 'value': float})
        if df.empty:
            print("错误: 稀疏矩阵文件为空。", file=sys.stderr); return False
        
        bins = np.unique(np.concatenate([df["start1"].values, df["start2"].values]))
        bin2idx = {b: i for i, b in enumerate(bins)}
        n = len(bins)
        mat = np.zeros((n, n), dtype=float)

        for start1, start2, value in df.values:
            i, j = bin2idx[start1], bin2idx[start2]
            mat[i, j] = mat[j, i] = value
        
        np.savetxt(output_file, mat, fmt="%.6f")
        print(f"[✓] 成功! {n}x{n} 密集矩阵已保存。")
        return True
    except Exception as e:
        print(f"错误: 转换到密集矩阵失败: {e}", file=sys.stderr); return False

def _convert_sparse_to_cscore_format(input_file, output_file, chrom, binsize, output_bed):
    """将三列稀疏矩阵转换为CscoreTool所需的七列格式。"""
    print(f"[*] 正在转换为CscoreTool格式: {os.path.basename(input_file)} -> {os.path.basename(output_file)}")
    max_coord, valid_line_count = 0, 0
    try:
        with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
            for line in fin:
                fields = line.strip().split()
                if len(fields)!= 3: continue
                try: start1, start2, count = int(fields[0]), int(fields[1]), float(fields[2])
                except ValueError: continue
                if count <= 0 or start1 == start2: continue
                max_coord = max(max_coord, start1, start2)
                fout.write(f"read{valid_line_count}\t{chrom}\t{start1}\t+\t{chrom}\t{start2}\t+\n")
                valid_line_count += 1
        
        if valid_line_count == 0: print("[!] 警告: 未找到任何有效的非对角线互作数据。")
        else: print(f"[✓] 成功! {valid_line_count} 行有效数据已保存。")

        if output_bed:
            if max_coord > 0:
                print(f"[*] 正在生成BED文件: {output_bed}")
                num_bins = (max_coord // binsize) + 1
                with open(output_bed, 'w') as f_bed:
                    for b in range(num_bins):
                        f_bed.write(f"{chrom}\t{b * binsize}\t{(b + 1) * binsize}\n")
                print(f"[✓] 成功! {num_bins} 个bins的坐标已保存。")
            else: print("[!] 警告: 未能确定最大坐标，无法生成BED文件。")
        return True
    except Exception as e:
        print(f"错误: 转换为CscoreTool格式失败: {e}", file=sys.stderr); return False

# ===================================================================
# 2. 最终的、合并后的命令行入口点
# ===================================================================

def main_hic_to_evr():
    """命令行工具: evr-prep-hic-to-evr"""
    parser = argparse.ArgumentParser(description="一键式工具：从.hic文件直接生成EVR所需的密集矩阵。")
    parser.add_argument("hic_file", help="输入的.hic文件路径。")
    parser.add_argument("juicer_jar", help="juicer_tools.jar文件的路径。")
    parser.add_argument("output_file", help="EVR所需的密集矩阵输出文件路径。")
    parser.add_argument("--chrom", required=True, help="要导出的染色体名称 (例如, 'chr1')。")
    parser.add_argument("--binsize", required=True, type=int, help="分辨率，单位为bp (例如, 100000)。")
    parser.add_argument("--norm", default="KR", help="要使用的矩阵归一化方法 (例如, KR, VC)。默认为KR。")
    args = parser.parse_args()

    # 使用临时文件来存储中间的稀疏矩阵
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".txt") as tmp_file:
        temp_sparse_path = tmp_file.name

    try:
        # 步骤 1: .hic -> 稀疏矩阵
        if not _run_juicer_dump(args.hic_file, args.juicer_jar, args.chrom, args.binsize, args.norm, temp_sparse_path):
            sys.exit(1) # 如果Juicer失败，则中止
        
        # 步骤 2: 稀疏矩阵 -> 密集矩阵
        if not _convert_sparse_to_dense(temp_sparse_path, args.output_file):
            sys.exit(1) # 如果转换失败，则中止
            
    finally:
        # 关键步骤: 确保无论成功与否，都删除临时文件
        if os.path.exists(temp_sparse_path):
            os.remove(temp_sparse_path)
            print(f"[*] 临时文件 '{temp_sparse_path}' 已被清理。")

def main_hic_to_cscore():
    """命令行工具: evr-prep-hic-to-cscore"""
    parser = argparse.ArgumentParser(description="一键式工具：从.hic文件直接生成CscoreTool所需的七列格式文件。")
    parser.add_argument("hic_file", help="输入的.hic文件路径。")
    parser.add_argument("juicer_jar", help="juicer_tools.jar文件的路径。")
    parser.add_argument("output_file", help="CscoreTool所需的七列格式输出文件。")
    parser.add_argument("--chrom", required=True, help="要导出的染色体名称 (例如, 'chr1')。")
    parser.add_argument("--binsize", required=True, type=int, help="分辨率，单位为bp (例如, 100000)。")
    parser.add_argument("--output-bed", help="可选：保存对应bin坐标的BED文件路径。")
    args = parser.parse_args()

    # CscoreTool需要原始计数，所以norm总是'NONE'
    norm_for_cscore = "NONE"

    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".txt") as tmp_file:
        temp_sparse_path = tmp_file.name
    
    try:
        # 步骤 1: .hic -> 稀疏矩阵 (原始计数)
        if not _run_juicer_dump(args.hic_file, args.juicer_jar, args.chrom, args.binsize, norm_for_cscore, temp_sparse_path):
            sys.exit(1)
        
        # 步骤 2: 稀疏矩阵 -> CscoreTool格式
        if not _convert_sparse_to_cscore_format(temp_sparse_path, args.output_file, args.chrom, args.binsize, args.output_bed):
            sys.exit(1)
            
    finally:
        if os.path.exists(temp_sparse_path):
            os.remove(temp_sparse_path)
            print(f"[*] 临时文件 '{temp_sparse_path}' 已被清理。")