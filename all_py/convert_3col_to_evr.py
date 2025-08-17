import sys
import numpy as np
import pandas as pd

def main():
    """
    主函数，用于从命令行读取文件路径并生成矩阵。
    """
    # 检查命令行参数数量
    if len(sys.argv) != 3:
        print(f"用法: python {sys.argv[0]} <输入文件名> <输出文件名>")
        sys.exit(1) # 退出程序

    input_path = sys.argv[1]
    output_path = sys.argv[2]

    try:
        # 读文件
        df = pd.read_csv(input_path, sep="\t", header=None, names=["start1", "start2", "value"])

        # 过滤NaN
        df = df.dropna()

        # 获取bins并创建索引
        bins = np.unique(np.concatenate([df["start1"].values, df["start2"].values]))
        bin2idx = {b: i for i, b in enumerate(bins)}

        n = len(bins)
        mat = np.zeros((n, n), dtype=float)

        # 填充矩阵
        for _, row in df.iterrows():
            i = bin2idx[row["start1"]]
            j = bin2idx[row["start2"]]
            mat[i, j] = row["value"]
            mat[j, i] = row["value"]

        # 保存矩阵
        np.savetxt(output_path, mat, fmt="%.6f")

        print(f"矩阵已成功从 '{input_path}' 创建并保存至 '{output_path}'")

    except FileNotFoundError:
        print(f"错误：找不到输入文件 '{input_path}'")
    except Exception as e:
        print(f"处理文件时发生错误: {e}")

if __name__ == "__main__":
    main()
