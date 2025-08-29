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
    """命令行工具: evr-plot3d-static"""
    parser = argparse.ArgumentParser(description="可视化PDB结构，并基于阈值进行区室着色。")
    # ... (参数定义与您的版本一致) ...
    parser.add_argument("pdbs", help="要可视化的PDB文件列表，用':'分隔。")
    parser.add_argument("-cf", "--compartment_files", help="对应的A/B区室得分文件列表，用':'分隔。")
    parser.add_argument("-o", "--output", help="输出图片文件名。")
    parser.add_argument("-lw", "--linewidth", default="2.0", help="线宽列表，用':'分隔。")
    parser.add_argument("-t", "--threshold", type=float, default=0.0, help="定义A/B区室的c-score绝对值阈值。")
    parser.add_argument("--superpose", action='store_true', help="将所有结构叠加到第一个结构上。")
    parser.add_argument("--rmsd", action='store_true', help="计算并打印RMSD值。")
    parser.add_argument("--show", action='store_true', help="在屏幕上显示图像。")
    args = parser.parse_args()

    pdb_files = args.pdbs.split(':')
    score_files = args.compartment_files.split(':') if args.compartment_files else []
    linewidths = args.linewidth.split(':')
    
    fig = plt.figure(figsize=(8, 8)); ax = fig.add_subplot(111, projection='3d')
    structures = []
    for i, pdb_file in enumerate(pdb_files):
        coords_raw = np.array(io_utils.read_pdb(pdb_file)).T
        scores_raw = io_utils.read_compartment_scores(score_files[i]) if i < len(score_files) else None
        min_len = min(len(coords_raw), len(scores_raw)) if scores_raw else len(coords_raw)
        structures.append({'coords': center_coords(coords_raw[:min_len]), 'scores': np.array(scores_raw[:min_len]) if scores_raw else None, 'lw': float(linewidths[i % len(linewidths)])})

    if args.superpose and len(structures) > 1:
        for i in range(1, len(structures)):
            # 修正: 正确的Kabsch算法调用
            rot_m = calculate_rotation_matrix(structures[0]['coords'], structures[i]['coords'])
            structures[i]['coords'] = np.dot(structures[i]['coords'], rot_m)

    for s in structures:
        coords, scores, lw = s['coords'], s['scores'], s['lw']
        colors = [get_compartment_color(score, args.threshold) for score in scores] if scores is not None else ['gray'] * len(coords)
        for i in range(len(coords) - 1):
            ax.plot(coords[i:i+2, 0], coords[i:i+2, 1], coords[i:i+2, 2], color=colors[i], linewidth=lw)

    if args.rmsd and len(structures) > 1:
        print("\n结构1\t结构2\tRMSD值")
        for i in range(len(structures)):
            for j in range(i + 1, len(structures)):
                print(f"{pdb_files[i]:<25}\t{pdb_files[j]:<25}\t{calculate_rmsd(structures[i]['coords'], structures[j]['coords']):.4f}")

    plt.axis('off')
    if args.output: plt.savefig(args.output, dpi=300, bbox_inches='tight')
    if args.show or not args.output: plt.show()

def main_animation():
    """命令行工具: evr-plot3d-animation"""
    parser = argparse.ArgumentParser(description="创建两个PDB结构之间的平滑过渡动画。")
    parser.add_argument("pdb1", help="起始PDB文件 (例如, GM12878)。")
    parser.add_argument("scores1", help="起始PDB的区室分数文件。")
    parser.add_argument("pdb2", help="结束PDB文件 (例如, K562)。")
    parser.add_argument("scores2", help="结束PDB的区室分数文件。")
    parser.add_argument("--name1", default="Cond1", help="条件1的名称。")
    parser.add_argument("--name2", default="Cond2", help="条件2的名称。")
    parser.add_argument("-o", "--output", required=True, help="输出视频文件名 (例如, morph.mp4)。")
    parser.add_argument("-f", "--frames", type=int, default=120, help="动画的总帧数。")
    parser.add_argument("-lw", "--linewidth", type=float, default=2.0, help="线条宽度。")
    parser.add_argument("-t", "--threshold", type=float, default=0.1, help="定义A/B区室的c-score绝对值阈值。")
    parser.add_argument("--skip-range", nargs=2, type=int, metavar=('S', 'E'), help="要跳过的1-based bin范围。")
    
    args = parser.parse_args()

    print("[*] 正在加载数据...")
    coords1_raw = np.array(io_utils.read_pdb(args.pdb1)).T
    scores1_raw = io_utils.read_compartment_scores(args.scores1)
    coords2_raw = np.array(io_utils.read_pdb(args.pdb2)).T
    scores2_raw = io_utils.read_compartment_scores(args.scores2)

    if args.skip_range:
        start_idx, end_idx = args.skip_range[0] - 1, args.skip_range[1] - 1
        print(f"正在从分数文件中跳过bins {args.skip_range[0]}-{args.skip_range[1]}...")
        scores1 = scores1_raw[:start_idx] + scores1_raw[end_idx + 1:]
        scores2 = scores2_raw[:start_idx] + scores2_raw[end_idx + 1:]
    else:
        scores1, scores2 = scores1_raw, scores2_raw

    min_len = min(len(coords1_raw), len(scores1), len(coords2_raw), len(scores2))
    print(f"检测到公共bin数量为 {min_len}。将所有数据截断至此长度。")
    
    coords1 = center_coords(coords1_raw[:min_len])
    coords2_centered = center_coords(coords2_raw[:min_len])
    scores1 = np.array(scores1[:min_len])
    scores2 = np.array(scores2[:min_len])

    print("[*] 正在对齐结构...")
    rot_matrix = calculate_rotation_matrix(coords1, coords2_centered)
    coords2_aligned = np.dot(coords2_centered, rot_matrix)

    print("[*] 正在准备动画...")
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    colors_t0 = np.array([get_compartment_color(s, args.threshold) for s in scores1])
    colors_t1 = np.array([get_compartment_color(s, args.threshold) for s in scores2])

    def update(frame):
        ax.cla()
        t = frame / (args.frames - 1)
        coords_interp = (1 - t) * coords1 + t * coords2_aligned
        current_colors = colors_t0 if t < 0.5 else colors_t1
        
        # 优化绘图
        start_idx, current_color = 0, current_colors[0]
        for i in range(1, len(coords_interp)):
            if current_colors[i] != current_color or i == len(coords_interp) - 1:
                segment = coords_interp[start_idx : i + 1]
                ax.plot(segment[:, 0], segment[:, 1], segment[:, 2], color=current_color, linewidth=args.linewidth)
                start_idx, current_color = i, current_colors[i]
        
        # 设置固定的视图范围以防止抖动
        all_coords = np.vstack([coords1, coords2_aligned])
        max_range = (np.ptp(all_coords, axis=0)).max() * 0.55
        mid = all_coords.mean(axis=0)
        ax.set_xlim(mid[0] - max_range, mid[0] + max_range)
        ax.set_ylim(mid[1] - max_range, mid[1] + max_range)
        ax.set_zlim(mid[2] - max_range, mid[2] + max_range)
        ax.axis('off')
        
        title_text = args.name1 if t < 0.5 else args.name2
        ax.set_title(f"Morphing to: {title_text} (Frame {frame+1}/{args.frames})")
        # 实时更新进度到终端
        sys.stdout.write(f"\r  > 正在渲染帧: {frame+1}/{args.frames}")
        sys.stdout.flush()

    print(f"[*] 正在创建并保存动画 ({args.frames} 帧)...")
    anim = FuncAnimation(fig, update, frames=args.frames, interval=50)
    
    try:
        anim.save(args.output, writer='ffmpeg', dpi=150, progress_callback=lambda i, n: sys.stdout.write(f"\r  > 正在编码帧: {i+1}/{n}"))
        print(f"\n[✓] 动画成功保存至 '{args.output}'")
    except Exception as e:
        print("\n\n--- 错误 ---", file=sys.stderr)
        print("动画保存失败。这通常意味着 'ffmpeg' 未被安装或找不到。", file=sys.stderr)
        print("如果您使用Conda, 请运行: 'conda install -c conda-forge ffmpeg'", file=sys.stderr)
        print(f"原始错误: {e}", file=sys.stderr)

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

def main_publication_plot():
    """
    命令行工具: evr-plot3d-publication
    生成一张静态的、带有预设视角和Top N基因高亮的出版级3D图像。
    """
    parser = argparse.ArgumentParser(description="生成带有基因表达热点高亮的出版级3D图像。")
    parser.add_argument("pdb", help="要绘制的PDB文件 (通常是第二个条件，如K562)。")
    parser.add_argument("scores1", help="条件1的C-score文件。")
    parser.add_argument("scores2", help="条件2的C-score文件。")
    parser.add_argument('--bin-gene-map', required=True, help="Bin到基因的映射文件。")
    parser.add_argument('--de-results', required=True, help="全基因组差异表达结果 (.csv)。")
    parser.add_argument("-o", '--output', required=True, help="输出图片文件名。")
    parser.add_argument("-lw", "--linewidth", type=float, default=1.5)
    parser.add_argument("--skip-range", nargs=2, type=int, metavar=('S', 'E'))
    parser.add_argument("--focus-range", nargs=2, type=int, metavar=('S', 'E'))
    parser.add_argument("-t", "--threshold", type=float, default=0.1)
    parser.add_argument('--show-bins', action='store_true', help="在高亮区域显示小球。")
    parser.add_argument('--bin-size', type=float, default=30.0, help="高亮bin的小球大小。")
    parser.add_argument('--top-genes', type=int, metavar='N', required=True, help="高亮显示包含Top N个上调基因的bins。")
    parser.add_argument('--view', nargs=2, type=float, default=[30, -60], metavar=('ELEV', 'AZIM'), help="设置相机视角 (仰角, 方位角)。")
    args = parser.parse_args()

    # --- 步骤 1: 加载和处理所有数据 ---
    print("[*] 正在加载和处理数据...")
    bin_to_gene, gene_to_bin = io_utils.parse_bin_gene_map(args.bin_gene_map)
    coords_raw = np.array(io_utils.read_pdb(args.pdb)).T
    scores1_raw = io_utils.read_compartment_scores(args.scores1)
    scores2_raw = io_utils.read_compartment_scores(args.scores2)
    try:
        de_df = pd.read_csv(args.de_results)
    except Exception as e:
        print(f"错误: 无法加载差异表达文件 '{args.de_results}'。错误: {e}"); sys.exit(1)

    # 过滤和截断
    if args.skip_range:
        start, end = args.skip_range[0]-1, args.skip_range[1]-1
        scores1, scores2 = scores1_raw[:start]+scores1_raw[end+1:], scores2_raw[:start]+scores2_raw[end+1:]
    else:
        scores1, scores2 = scores1_raw, scores2_raw
    min_len = min(len(coords_raw), len(scores1), len(scores2))
    coords, scores1, scores2 = coords_raw[:min_len], np.array(scores1[:min_len]), np.array(scores2[:min_len])
    coords_centered = center_coords(coords)
    
    # --- 步骤 2: 识别Top基因并准备颜色 ---
    print(f"[*] 正在识别Top {args.top_genes} 上调基因...")
    bin_categories = categorize_bin_switches(scores1, scores2, args.threshold)
    background_color = 'lightgray' # 使用一个更中性的背景色
    line_colors = np.array([background_color] * len(coords), dtype=object)
    
    genes_on_this_chromosome = list(gene_to_bin.keys())
    de_df_chr = de_df[de_df['gene_name'].isin(genes_on_this_chromosome)]
    highlight_bins = set()

    if not de_df_chr.empty:
        top_genes_df = de_df_chr.sort_values(by='log2FoldChange', ascending=False).head(args.top_genes)
        print(f"  - 找到 {len(top_genes_df)} 个Top基因进行高亮:")
        for _, row in top_genes_df.iterrows():
            gene, log2fc = row['gene_name'], row['log2FoldChange']
            if gene in gene_to_bin:
                bin_idx = gene_to_bin[gene]
                if bin_idx < len(line_colors):
                    category = bin_categories[bin_idx]
                    color = 'red' if 'B -> A' in category else ('blue' if 'A -> B' in category else 'yellow')
                    line_colors[bin_idx] = color
                    highlight_bins.add(bin_idx)
                    print(f"    - Bin {bin_idx+1} (状态: {category}) -> 基因: {gene} (log2FC: {log2fc:.2f}). 将高亮为 {color.upper()}.")
    else:
        print("  - 警告: 在差异表达结果中未找到任何与当前染色体匹配的基因。")

    # --- 步骤 3: 绘图 ---
    print("[*] 正在生成图像...")
    focus_start, focus_end = (args.focus_range[0]-1, args.focus_range[1]) if args.focus_range else (0, len(coords))
    coords_to_plot, colors_to_plot = coords_centered[focus_start:focus_end], line_colors[focus_start:focus_end]

    fig = plt.figure(figsize=(12, 12)); ax = fig.add_subplot(111, projection='3d')

    if len(coords_to_plot) > 1:
        # 绘制背景线条
        ax.plot(coords_to_plot[:,0], coords_to_plot[:,1], coords_to_plot[:,2], color=background_color, linewidth=args.linewidth, zorder=1)
        # 单独绘制高亮线段，使其更突出
        for i in range(len(coords_to_plot) - 1):
            if colors_to_plot[i] != background_color or colors_to_plot[i+1] != background_color:
                ax.plot(coords_to_plot[i:i+2, 0], coords_to_plot[i:i+2, 1], coords_to_plot[i:i+2, 2], 
                        color=colors_to_plot[i], linewidth=args.linewidth * 1.5, zorder=2)

    if args.show_bins:
        highlight_coords = coords_to_plot[[idx - focus_start for idx in highlight_bins if focus_start <= idx < focus_end]]
        highlight_colors = line_colors[[idx for idx in highlight_bins if focus_start <= idx < focus_end]]
        if highlight_coords.size > 0:
            ax.scatter(highlight_coords[:,0], highlight_coords[:,1], highlight_coords[:,2], 
                       c=highlight_colors, s=args.bin_size, depthshade=True, alpha=1.0, zorder=3, edgecolor='black')

    # --- 步骤 4: 设置视图并保存 ---
    print("[*] 正在设置最终视图并保存...")
    ax.set_box_aspect([1,1,1]); ax.view_init(elev=args.view[0], azim=args.view[1]); ax.axis('off')
    plt.savefig(args.output, dpi=300, bbox_inches='tight', transparent=True)
    print(f"\n[✓] 分析完成! 出版级图像已保存至 '{args.output}'")