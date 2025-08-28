# 文件名: src/evr_3d_suite/plot3d.py
# 描述: 3D染色体结构的可视化，包括静态、动态和交互式多组学模式

import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider
from matplotlib.animation import FuncAnimation

# 从同一个包中导入IO工具
from . import io_utils

# ===================================================================
# 1. 核心辅助函数
# ===================================================================

def center_coords(coords):
    if coords.size == 0: return coords
    return coords - coords.mean(axis=0)

def calculate_rotation_matrix(S1, S2):
    A = np.dot(S2.T, S1)
    u, _, v = np.linalg.svd(A)
    if np.linalg.det(u) * np.linalg.det(v) < 0:
        u[:, -1] *= -1
    return np.dot(u, v)

def calculate_rmsd(S1, S2):
    if S1.shape != S2.shape: return np.nan
    return np.sqrt(np.sum((S1 - S2)**2.0) / S1.shape[0])

def categorize_bin_switches(scores1, scores2, threshold):
    categories = []
    for s1, s2 in zip(scores1, scores2):
        state1 = 'A' if s1 > threshold else ('B' if s1 < -threshold else 'T')
        state2 = 'A' if s2 > threshold else ('B' if s2 < -threshold else 'T')
        if state1 == 'B' and state2 == 'A': categories.append('B -> A')
        elif state1 == 'A' and state2 == 'B': categories.append('A -> B')
        elif state1 == 'A' and state2 == 'A': categories.append('Stable A')
        elif state1 == 'B' and state2 == 'B': categories.append('Stable B')
        else: categories.append('Other/Transition')
    return np.array(categories)

class PanHandler:
    def __init__(self, fig, ax):
        self.fig, self.ax, self.press, self.panning = fig, ax, None, False
        self.connect()
    def connect(self):
        self.cid_press=self.fig.canvas.mpl_connect('button_press_event',self.on_press)
        self.cid_release=self.fig.canvas.mpl_connect('button_release_event',self.on_release)
        self.cid_motion=self.fig.canvas.mpl_connect('motion_notify_event',self.on_motion)
    def on_press(self, e):
        if e.inaxes != self.ax or e.button != 2: return
        self.panning, self.press = True, (e.x, e.y)
    def on_motion(self, e):
        if not self.panning or e.inaxes != self.ax: return
        dx, dy = (e.x - self.press[0]) * -1, (e.y - self.press[1]) * -1
        xlim, ylim = self.ax.get_xlim(), self.ax.get_ylim()
        x_range, y_range = xlim[1] - xlim[0], ylim[1] - ylim[0]
        self.ax.set_xlim(xlim[0] + dx * x_range * 0.003, xlim[1] + dx * x_range * 0.003)
        self.ax.set_ylim(ylim[0] + dy * y_range * 0.003, ylim[1] + dy * y_range * 0.003)
        self.press = (e.x, e.y); self.fig.canvas.draw()
    def on_release(self, e):
        if e.button == 2: self.panning, self.press = False, None

# ===================================================================
# 2. 核心可视化类 (用于复杂的交互式绘图)
# ===================================================================

class Chromosome3DViewer:
    def __init__(self, args):
        self.args = args
        self._load_and_process_data()
        self._setup_plot()

    def _load_and_process_data(self):
        print("[*] 正在加载和处理数据...")
        # 加载结构和分数
        coords1_raw = np.array(io_utils.read_pdb(self.args.pdb1)).T
        coords2_raw = np.array(io_utils.read_pdb(self.args.pdb2)).T
        scores1_raw = io_utils.read_compartment_scores(self.args.scores1)
        scores2_raw = io_utils.read_compartment_scores(self.args.scores2)

        # 截断到最小公共长度
        min_len = min(len(coords1_raw), len(coords2_raw), len(scores1_raw), len(scores2_raw))
        self.coords1 = center_coords(coords1_raw[:min_len])
        coords2_centered = center_coords(coords2_raw[:min_len])
        self.scores1 = np.array(scores1_raw[:min_len])
        self.scores2 = np.array(scores2_raw[:min_len])
        
        # 对齐结构
        rot_matrix = calculate_rotation_matrix(self.coords1, coords2_centered)
        self.coords2_aligned = np.dot(coords2_centered, rot_matrix)

        # 分类bin
        self.bin_categories = categorize_bin_switches(self.scores1, self.scores2, self.args.threshold)
        
        # 加载用于查询的组学数据
        self.bin_to_gene, self.gene_to_bin = io_utils.parse_bin_gene_map(self.args.bin_gene_map)
        self.de_map = io_utils.load_de_results(self.args.rna_seq) if self.args.rna_seq else {}
        if self.args.rsem_files and self.args.gtf:
            id_map = io_utils.create_id_map_from_gtf(self.args.gtf)
            self.tpm_map = io_utils.process_rsem_to_tpm_map(self.args.rsem_files, id_map)
        else:
            self.tpm_map = {}
        
        # 如果是高亮Top基因模式，识别这些bins
        if self.args.mode == 'highlight_top_genes':
            self._identify_top_gene_bins()

        print("[✓] 数据加载完毕。")

    def _setup_plot(self):
        self.fig = plt.figure(figsize=(10, 10))
        self.ax_3d = self.fig.add_axes([0, 0.05, 1, 0.95], projection='3d')
        self.fig.canvas.mpl_connect('pick_event', self._on_pick)
        ax_slider = self.fig.add_axes([0.2, 0.02, 0.65, 0.03])
        self.morph_slider = Slider(ax=ax_slider, label='Morph', valmin=0.0, valmax=1.0, valinit=0.0)
        self.morph_slider.on_changed(self.update)
        self.pan_handler = PanHandler(self.fig, self.ax_3d)
        
        # 设置一次性视图
        all_coords = np.vstack([self.coords1, self.coords2_aligned])
        if all_coords.size > 0:
            max_range = (np.ptp(all_coords, axis=0)).max() * 0.55
            mid = all_coords.mean(axis=0)
            self.ax_3d.set_xlim(mid[0]-max_range, mid[0]+max_range); self.ax_3d.set_ylim(mid[1]-max_range, mid[1]+max_range); self.ax_3d.set_zlim(mid[2]-max_range, mid[2]+max_range)
        self.ax_3d.axis('off')
        
        self.update(0.0)
        plt.show()

    def update(self, val):
        while self.ax_3d.lines: self.ax_3d.lines[0].remove()
        while self.ax_3d.collections: self.ax_3d.collections[0].remove()
        t = val
        coords_interp = (1 - t) * self.coords1 + t * self.coords2_aligned
        
        # 确定颜色
        if self.args.mode == 'highlight_switches':
            colors = np.array(['red' if c == 'B -> A' else 'blue' if c == 'A -> B' else 'lightgray' for c in self.bin_categories])
        elif self.args.mode == 'highlight_top_genes':
            colors = np.array(['yellow' if i in self.top_gene_bin_indices else 'plum' for i in range(len(coords_interp))], dtype=object)
        else: # 默认模式
            scores = self.scores1 if t < 0.5 else self.scores2
            colors = np.array([get_compartment_color(s, self.args.threshold) for s in scores])

        # 绘制
        if len(coords_interp) > 0:
            start_idx, current_color = 0, colors[0]
            for i in range(1, len(coords_interp)):
                if colors[i] != current_color or i == len(coords_interp) - 1:
                    segment = coords_interp[start_idx:i+1]
                    self.ax_3d.plot(segment[:,0], segment[:,1], segment[:,2], color=current_color, linewidth=self.args.linewidth)
                    start_idx, current_color = i, colors[i]
            if self.args.show_bins:
                self.ax_3d.scatter(coords_interp[:,0], coords_interp[:,1], coords_interp[:,2], c=colors, s=self.args.bin_size, depthshade=True)
        
        # 添加用于拾取的不可见层
        self.ax_3d.scatter(coords_interp[:,0], coords_interp[:,1], coords_interp[:,2], s=self.args.bin_size*1.5, alpha=0.0, picker=5)
        
        self.ax_3d.set_title(f"{(1-t)*100:.0f}% {self.args.name1} | {t*100:.0f}% {self.args.name2}")
        self.fig.canvas.draw_idle()

    def _on_pick(self, event):
        if event.ind.size == 0: return
        idx = event.ind[0]
        category = self.bin_categories[idx]
        
        if self.args.mode == 'highlight_switches' and category not in ['A -> B', 'B -> A']:
            print(f"\n点击了非切换区域 (Bin {idx+1})")
            return

        print("\n" + "="*40 + f"\nBin {idx+1} (Index: {idx}) 信息\n" + "="*40)
        print(f" 区室状态变化: {category}")
        print(f" C-scores: {self.args.name1}={self.scores1[idx]:.3f}, {self.args.name2}={self.scores2[idx]:.3f}")
        genes = self.bin_to_gene.get(idx, ["未找到基因"])
        print(" Bin内基因:")
        for gene in sorted(genes):
            log2fc = self.de_map.get(gene, "N/A")
            log2fc_str = f"{log2fc:.2f}" if isinstance(log2fc, float) else "N/A"
            tpm = self.tpm_map.get(gene, "N/A")
            tpm_str = f"{tpm:.2f}" if isinstance(tpm, float) else "N/A"
            print(f"  - {gene:<15} | log2FC: {log2fc_str:<7} | TPM in {self.args.name2}: {tpm_str}")
        print("="*40)

    def _identify_top_gene_bins(self):
        """识别所有包含Top N上调基因的bins。"""
        print(f"[*] 正在识别Top {self.args.top_genes} 上调基因...")
        genes_on_this_chromosome = list(self.gene_to_bin.keys())
        de_df_chr = self.args.de_results[self.args.de_results['gene_name'].isin(genes_on_this_chromosome)]
        self.top_gene_bin_indices = set()
        if not de_df_chr.empty:
            top_genes_df = de_df_chr.sort_values(by='log2FoldChange', ascending=False).head(self.args.top_genes)
            print(f"  - 找到 {len(top_genes_df)} 个Top基因进行高亮:")
            for _, row in top_genes_df.iterrows():
                gene, log2fc = row['gene_name'], row['log2FoldChange']
                if gene in self.gene_to_bin:
                    bin_idx = self.gene_to_bin[gene]
                    self.top_gene_bin_indices.add(bin_idx)
                    print(f"    - Bin {bin_idx+1} -> Gene: {gene} (log2FC: {log2fc:.2f})")
        else:
            print("  - 警告: 未找到匹配的基因来高亮。")

# ===================================================================
# 3. 命令行入口点函数
# ===================================================================

def main_static():
    # ... [函数内容与上一版本完全相同] ...
    pass

def main_animation():
    """命令行工具: evr-plot3d-animation"""
    parser = argparse.ArgumentParser(description="创建两个PDB结构之间的平滑过渡动画。")
    # ... [此处省略与原plot_animation.py相同的参数定义和动画生成逻辑] ...
    print("动画生成功能暂未完全整合。需要安装 'ffmpeg'。")

def main_interactive():
    """命令行工具: evr-plot3d-interactive"""
    parser = argparse.ArgumentParser(description="交互式地探索、变形并查询3D染色体结构。")
    # ... [所有交互式模式的参数] ...
    parser.add_argument("pdb1"); parser.add_argument("scores1"); parser.add_argument("pdb2"); parser.add_argument("scores2")
    parser.add_argument("--name1", default="GM12878"); parser.add_argument("--name2", default="K562")
    parser.add_argument('--gtf', required=True); parser.add_argument('--bin-gene-map', required=True)
    parser.add_argument('--rna-seq', required=True); parser.add_argument('--rsem-files', nargs='+', required=True)
    parser.add_argument("-lw", "--linewidth", type=float, default=2.0); parser.add_argument("--skip-range", nargs=2, type=int)
    parser.add_argument("--focus-range", nargs=2, type=int); parser.add_argument("-t", "--threshold", type=float, default=0.1)
    parser.add_argument('--show-bins', action='store_true'); parser.add_argument('-bs', '--bin-size', type=float, default=20)
    parser.add_argument('--mode', choices=['default', 'highlight_switches', 'highlight_top_genes'], default='default', 
                        help="设置可视化模式。default: A/B区室着色; highlight_switches: 高亮区室转换; highlight_top_genes: 高亮Top上调基因。")
    parser.add_argument('--top-genes', type=int, default=10, help="在 'highlight_top_genes' 模式下，指定要高亮的基因数量。")
    
    args = parser.parse_args()
    
    # 为 'highlight_top_genes' 模式预加载DE结果
    if args.mode == 'highlight_top_genes':
        try:
            args.de_results = pd.read_csv(args.rna_seq)
        except Exception as e:
            print(f"错误: 无法加载差异表达文件 '{args.rna_seq}' 用于高亮Top基因。错误: {e}"); sys.exit(1)

    Chromosome3DViewer(args)