#!/usr/bin/env python3
import sys
import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt

def main():
    """
    主函数，用于解析命令行参数并从 bedgraph 文件绘制 C-score 曲线图。
    """
    # --- 1. 使用 argparse 替换 sys.argv，提供更强大的命令行接口 ---
    parser = argparse.ArgumentParser(
        description="A script to plot C-score profiles from a bedgraph file.",
        formatter_class=argparse.RawTextHelpFormatter # 保持帮助信息格式
    )
    
    # 定义位置参数 (必需)
    parser.add_argument(
        "input_file", 
        help="Path to the input bedgraph file."
    )
    
    # 定义选项参数 (可选)
    parser.add_argument(
        "-o", "--output",
        help="Path to save the output image file.\n(Default: automatically derived from input file, e.g., 'input.png')"
    )
    parser.add_argument(
        "-t", "--title",
        help="Custom title for the plot.\n(Default: based on the input filename)"
    )
    # --- 2. 新增 --show 选项 ---
    # action="store_true" 表示这是一个开关，如果出现，其值为True
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display the plot interactively instead of just saving it."
    )
    
    args = parser.parse_args()

    # --- 3. 增加健壮的错误处理 ---
    try:
        print(f"[*] Reading data from: {args.input_file}")
        # comment='t' 会跳过以 't' 开头的行 (例如 bedgraph 的 track line)
        df = pd.read_csv(
            args.input_file, 
            sep=r"\s+", 
            header=None, 
            comment="t", 
            engine="python"
        )
        df.columns = ["chr", "start", "end", "cscore"]
        
        if df.empty:
            print("[!] Warning: The input file is empty or contains no valid data rows.")
            return

    except FileNotFoundError:
        print(f"[-] Error: Input file not found at '{args.input_file}'", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"[-] An error occurred while reading the file: {e}", file=sys.stderr)
        sys.exit(1)

    # --- 4. 智能处理标题和输出文件名 ---
    # 如果用户未提供标题，则根据输入文件名自动生成
    if args.title:
        plot_title = args.title
    else:
        # e.g., "my_data.bedgraph" -> "C-score profile for my_data"
        base_name = os.path.basename(args.input_file)
        clean_name = os.path.splitext(base_name)[0]
        plot_title = f"C-score profile for {clean_name}"
        
    # 如果用户未提供输出文件名，则根据输入文件名自动生成
    if args.output:
        output_path = args.output
    else:
        # e.g., "/path/to/my_data.bedgraph" -> "/path/to/my_data.png"
        base_path = os.path.splitext(args.input_file)[0]
        output_path = base_path + ".png"

    # --- 5. 绘图逻辑 (加入美化) ---
    print("[*] Generating plot...")
    plt.figure(figsize=(15, 5))
    plt.plot(df["start"], df["cscore"], color="#3477eb", linewidth=0.9)
    
    # 添加参考线和网格
    plt.axhline(0, ls="--", color="gray", zorder=0)
    plt.grid(True, linestyle=':', linewidth=0.6, alpha=0.7)
    
    # 设置标题和标签
    plt.title(plot_title, fontsize=16, weight='bold')
    plt.xlabel("Genomic Position (bp)", fontsize=12)
    plt.ylabel("C-score", fontsize=12)
    
    # 优化坐标轴显示，使用科学记数法显示大数字
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    
    plt.tight_layout()

    # --- 6. 实现 --show 和保存逻辑 ---
    # 总是先保存文件
    try:
        print(f"[*] Saving plot to: {output_path}")
        plt.savefig(output_path, dpi=300)
    except Exception as e:
        print(f"[-] Error saving file: {e}", file=sys.stderr)

    # 如果用户指定了 --show，则显示图像
    if args.show:
        print("[*] Displaying plot. Close the plot window to exit.")
        plt.show()

    # 释放内存
    plt.close()
    
    print(f"[✓] Done.")

if __name__ == "__main__":
    main()