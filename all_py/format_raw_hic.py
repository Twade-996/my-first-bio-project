#!/usr/bin/env python3
import argparse
import sys

def main():
    """
    一个整合脚本，用于将三列的原始Hi-C矩阵文件转换为特定的七列格式，
    并可以选择性地生成一个与之匹配的基因组坐标BED文件。
    """
    # --- 使用 argparse 设置更强大的命令行界面 ---
    parser = argparse.ArgumentParser(
        description="Converts a 3-column raw interaction matrix (start1, start2, count) into a 7-column formatted file and optionally generates a corresponding BED file.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    # --- 位置参数 (必需) ---
    parser.add_argument("input_file", help="Path to the input raw matrix file (3 columns).")
    parser.add_argument("output_file", help="Path for the final formatted output file.")
    
    # --- 选项参数 ---
    parser.add_argument("--chrom", required=True, help="Chromosome name (e.g., 'chr1', 'chrX').")
    parser.add_argument("--binsize", required=True, type=int, help="Resolution in base pairs (e.g., 100000).")
    # --- 新增功能：可选的BED文件输出 ---
    parser.add_argument("--output_bed", help="(Optional) Path to save the corresponding BED file of genomic bins.")
    
    args = parser.parse_args()

    print(f"[*] Processing file: {args.input_file}")
    print(f"[*] Using chromosome: {args.chrom} with binsize: {args.binsize} bp")

    # 初始化一个变量，用于在遍历文件时追踪遇到的最大坐标
    max_coord = 0
    # 记录处理的行数，用于判断文件是否为空
    line_count = -1

    try:
        with open(args.input_file, 'r') as fin, open(args.output_file, 'w') as fout:
            for line_count, line in enumerate(fin):
                fields = line.strip().split()
                if len(fields) != 3:
                    continue

                start1_str, start2_str, count_str = fields
                
                try:
                    start1 = int(start1_str)
                    start2 = int(start2_str)
                    count = float(count_str)
                except ValueError:
                    continue
                
                if count <= 0:
                    continue

                # --- 核心整合点：在处理数据的同时，更新最大坐标 ---
                current_max = max(start1, start2)
                if current_max > max_coord:
                    max_coord = current_max

                # 执行原有的格式化任务
                read_name = f"read{line_count}"
                fout.write(f"{read_name}\t{args.chrom}\t{start1}\t+\t{args.chrom}\t{start2}\t+\n")

        print(f"[✓] Success! Formatted data saved to: {args.output_file}")

    except FileNotFoundError:
        print(f"错误：找不到输入文件 '{args.input_file}'", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"处理文件时发生错误: {e}", file=sys.stderr)
        sys.exit(1)

    # --- 新增功能：在文件处理完毕后，生成BED文件 ---
    # 检查用户是否请求了 BED 文件输出
    if args.output_bed:
        # 只有在文件包含有效数据时才生成
        if line_count == -1 or max_coord == 0:
            print(f"[!] Warning: Input file '{args.input_file}' was empty or contained no valid data. BED file will not be generated.")
        else:
            print(f"[*] Generating BED file: {args.output_bed}")
            try:
                # 根据找到的最大坐标和分辨率，计算总的bin数量
                num_bins = (max_coord // args.binsize) + 1
                with open(args.output_bed, 'w') as f_bed:
                    for b in range(num_bins):
                        start = b * args.binsize
                        end = start + args.binsize
                        f_bed.write(f"{args.chrom}\t{start}\t{end}\n")
                print(f"[✓] Success! BED file saved to: {args.output_bed}")
            except Exception as e:
                print(f"写入BED文件时发生错误: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()
