#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
auto_plot_images.py

自动读取指定目录下所有 PNG 图片，并在一个窗口中排列显示。
"""

import os
import glob
import math
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import argparse

def main():
    parser = argparse.ArgumentParser(description="自动读取指定目录所有 PNG 图片并排列显示")
    parser.add_argument("img_dir", help="图片目录路径")
    parser.add_argument("--save", help="是否保存整合后的图片，例如 output.png", default=None)
    args = parser.parse_args()

    img_dir = args.img_dir
    img_files = sorted(glob.glob(os.path.join(img_dir, "*.png")))

    if not img_files:
        print(f"目录 {img_dir} 下没有找到 PNG 图片！")
        return

    n = len(img_files)
    # 自动计算行列数：尽量接近平方形
    cols = math.ceil(math.sqrt(n))
    rows = math.ceil(n / cols)

    fig, axes = plt.subplots(rows, cols, figsize=(4*cols, 3*rows))
    axes = axes.flatten()  # 方便循环

    for ax, img_file in zip(axes, img_files):
        img = mpimg.imread(img_file)
        ax.imshow(img)
        ax.axis("off")
        ax.set_title(os.path.basename(img_file), fontsize=10)

    # 隐藏多余子图
    for i in range(len(img_files), len(axes)):
        axes[i].axis("off")

    plt.tight_layout()
    plt.show()

    if args.save:
        fig.savefig(args.save, dpi=300)
        print(f"整合图片已保存为 {args.save}")

if __name__ == "__main__":
    main()
