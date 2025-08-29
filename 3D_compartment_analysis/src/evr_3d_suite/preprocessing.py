# 文件名: src/evr_3d_suite/preprocessing.py
# 描述: Hi-C数据预处理和格式转换工具集

import argparse
import sys
import numpy as np
import pandas as pd

def main_juicer_to_evr():
    """
    命令行工具: evr-juicer-to-evr
    将Juicer dump的三列稀疏矩阵转换为EVR所需的密集矩阵。
    """
    parser = argparse.ArgumentParser(description="将Juicer的三列稀疏矩阵转换为EVR所需的密集矩阵。")
    parser.add_argument("input_file", help="Juicer dump的输入文件路径 (start1, start2, value)。")
    parser.add_argument("output_file", help="EVR所需的密集矩阵输出文件路径。")
    args = parser.parse_args()

    print(f"[*] 正在从 '{args.input_file}' 读取稀疏矩阵...")
    try:
        df = pd.read_csv(args.input_file, sep=r"\s+", header=None, names=["start1", "start2", "value"])
        df = df.dropna()
        if df.empty:
            print("错误: 输入文件为空或格式不正确，读取数据失败。", file=sys.stderr)
            sys.exit(1)

        # 获取所有唯一的bin坐标并创建索引映射
        bins = np.unique(np.concatenate([df["start1"].values, df["start2"].values])).astype(int)
        bin2idx = {b: i for i, b in enumerate(bins)}

        n = len(bins)
        mat = np.zeros((n, n), dtype=float)

        print(f"[*] 发现 {n} 个唯一的bins，正在构建 {n}x{n} 密集矩阵...")
        # 填充矩阵
        for _, row in df.iterrows():
            i = bin2idx[int(row["start1"])]
            j = bin2idx[int(row["start2"])]
            mat[i, j] = row["value"]
            mat[j, i] = row["value"] # 确保矩阵对称

        # 保存矩阵
        np.savetxt(args.output_file, mat, fmt="%.6f")
        print(f"[✓] 成功! 密集矩阵已保存至 '{args.output_file}'")

    except FileNotFoundError:
        print(f"错误: 找不到输入文件 '{args.input_file}'", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"处理文件时发生错误: {e}", file=sys.stderr)
        sys.exit(1)

def main_juicer_to_cscore():
    """
    命令行工具: evr-juicer-to-cscore
    将Juicer的原始三列矩阵转换为CscoreTool所需的七列格式，并可选生成bin的BED文件。
    """
    parser = argparse.ArgumentParser(
        description="将Juicer的三列原始矩阵转换为CscoreTool所需的七列格式。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("input_file", help="Juicer dump的输入原始矩阵文件 (start1, start2, count)。")
    parser.add_argument("output_file", help="CscoreTool所需的七列格式输出文件。")
    parser.add_argument("--chrom", required=True, help="染色体名称 (例如, 'chr1')。")
    parser.add_argument("--binsize", required=True, type=int, help="分辨率，单位为bp (例如, 100000)。")
    parser.add_argument("--output-bed", help="可选：保存对应bin坐标的BED文件路径。")
    args = parser.parse_args()

    print(f"[*] 正在处理文件: {args.input_file}")
    print(f"[*] 使用染色体: {args.chrom}，分辨率: {args.binsize} bp")

    max_coord = 0
    line_count = 0
    valid_line_count = 0
    try:
        with open(args.input_file, 'r') as fin, open(args.output_file, 'w') as fout:
            for line_count, line in enumerate(fin, 1):
                fields = line.strip().split()
                if len(fields) != 3: continue
                try:
                    start1, start2, count = int(fields[0]), int(fields[1]), float(fields[2])
                except ValueError:
                    continue
                if count <= 0:
                    continue
                
                max_coord = max(max_coord, start1, start2)
                read_name = f"read{valid_line_count}"
                fout.write(f"{read_name}\t{args.chrom}\t{start1}\t+\t{args.chrom}\t{start2}\t+\n")
                valid_line_count += 1
        
        if valid_line_count == 0:
            print(f"[!] 警告: 输入文件 '{args.input_file}' 中未找到任何有效的互作数据行。")
        else:
            print(f"[✓] 成功! 已处理 {valid_line_count} 行有效数据，已保存至: {args.output_file}")

    except FileNotFoundError:
        print(f"错误: 找不到输入文件 '{args.input_file}'", file=sys.stderr); sys.exit(1)
    except Exception as e:
        print(f"处理文件时发生错误: {e}", file=sys.stderr); sys.exit(1)

    if args.output_bed:
        if max_coord == 0:
            print(f"[!] 警告: 未能从输入文件中确定最大坐标。将不会生成BED文件。")
        else:
            print(f"[*] 正在生成BED文件: {args.output_bed}")
            try:
                num_bins = (max_coord // args.binsize) + 1
                with open(args.output_bed, 'w') as f_bed:
                    for b in range(num_bins):
                        start = b * args.binsize
                        end = start + args.binsize
                        f_bed.write(f"{args.chrom}\t{start}\t{end}\n")
                print(f"[✓] 成功! {num_bins} 个bins的坐标已保存至: {args.output_bed}")
            except Exception as e:
                print(f"写入BED文件时发生错误: {e}", file=sys.stderr)