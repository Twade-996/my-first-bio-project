# EVR 3D Suite: A Chromosome 3D Analysis & Visualization Toolkit

![EVR 3D Suite Banner](https://i.imgur.com/your_banner_image_url.png) 
<!-- 这是一个占位符, 您可以自己创建一个简单的Logo或Banner图上传到imgur等图床替换链接 -->

**EVR 3D Suite** 是一个全面的Python软件包，用于分析和可视化Hi-C数据衍生的3D染色质结构。它提供了一套集成的命令行工具，能够处理从原始`.hic`文件到最终出版级3D交互式模型的全部分析流程。

本工具套件整合了3D基因组学分析的几个关键阶段：

-   **数据预处理**: 直接从`.hic`文件开始，一键生成下游分析所需的各种矩阵格式。
-   **区室分析**: 对A/B区室分数进行统计分析、比较和可视化。
-   **多组学整合**: 将染色质结构动态与RNA-seq基因表达数据（差异变化与绝对水平）进行深度关联。
-   **高级3D可视化**: 创建可交互的、可在两种细胞状态间平滑“变形”的3D模型，并高亮显示结构与功能上的关键变化区域。

---

## 目录

- [安装](#安装)
- [分析工作流程](#分析工作流程)
  - [步骤 0: 数据预处理 (从 .hic 文件开始)](#步骤-0-数据预处理)
  - [步骤 1: 分析区室分数](#步骤-1-分析区室分数)
  - [步骤 2: 关联结构与基因表达](#步骤-2-关联结构与基因表达)
  - [步骤 3: 3D交互式可视化](#步骤-3-3d交互式可视化)
- [命令行工具参考](#命令行工具参考)
- [致谢与引用](#致谢与引用)

---

## 安装

我们强烈建议在一个专用的Conda环境中安装本软件包。

1.  **创建并激活新的Conda环境:**
    ```bash
    conda create -n evr_env python=3.9
    conda activate evr_env
    ```

2.  **安装核心依赖:**
    本软件包依赖于一些常见的科学计算库，以及用于动画制作的`ffmpeg`。
    ```bash
    # 安装Python依赖
    pip install numpy pandas seaborn matplotlib
    # 安装ffmpeg
    conda install -c conda-forge ffmpeg
    ```

3.  **克隆并安装本软件包:**
    克隆此代码仓库，并以“可编辑”模式安装`evr_3d_suite`。这使您对源代码的任何修改都能立即生效。
    ```bash
    git clone https://github.com/your_username/your_repository_name.git
    cd your_repository_name # 进入项目根目录
    pip install -e .
    ```
    > 安装完成后，所有命令行工具 (例如, `evr-hic2evr`, `evr-plot3d-interactive`) 将在您的终端中直接可用。

---

## 分析工作流程

本节将引导您完成一个从原始数据到最终可视化的典型分析流程。

### 步骤 0: 数据预处理 (从 .hic 文件开始)

这是整个分析的起点。您需要一个`.hic`文件和一个`juicer_tools.jar`文件。

**A. 为 EVR 3D 建模准备密集矩阵:**
此命令直接从`.hic`文件中提取KR归一化的互作矩阵，并转换为EVR所需的密集矩阵格式。
```bash
evr-hic2evr path/to/your_data.hic path/to/juicer_tools.jar \
    results/evr_matrix.txt \
    --chrom 1 \
    --binsize 50000 \
    --norm KR
```
> **提示**: `.hic`文件中的染色体命名可能是 `1`, `2`... 或 `chr1`, `chr2`...。请使用 `juicer_tools dump headers` 命令检查您的`.hic`文件以确认。

**B. 为 CscoreTool 分析准备输入文件:**
此命令从`.hic`文件中提取原始互作计数，转换为CscoreTool所需的7列格式，并同时生成一个定义bin坐标的BED文件，这对于后续所有分析至关重要。
```bash
evr-hic2cscore path/to/your_data.hic path/to/juicer_tools.jar \
    results/cscore_input.summary \
    --chrom 1 \
    --binsize 50000 \
    --output-bed results/chr1_50kb_bins.bed
```

**在完成这些预处理后，您就可以使用EVR和CscoreTool来生成后续分析所需的`.pdb`和`.cscore.txt`文件了。**

### 步骤 1: 分析区室分数

在进行复杂分析前，先了解您的区室数据。

**A. 探索C-score分布并确定阈值:**
```bash
evr-cscore-threshold results/your_cscore.txt --skip-range 2430 2578
```
> 此工具会输出统计摘要并绘制直方图，帮助您为定义A/B区室选择一个合理的阈值 (例如 `-t 0.1`)。

**B. 比较两种条件下的区室分数:**
```bash
evr-cscore-compare \
    path/to/gm12878_cscore.bedgraph \
    path/to/k562_cscore.bedgraph \
    --name1 GM12878 \
    --name2 K562 \
    --outdir results/cscore_comparison
```

### 步骤 2: 关联结构与基因表达

这是将结构与功能联系起来的核心分析步骤。

**A. 关联区室*转换*与基因表达*变化***:
```bash
evr-analyze-change \
    --bin-gene-map path/to/bin_gene_map.txt \
    --rna-seq path/to/rna_seq_de_results.csv \
    --scores1 path/to/gm12878_cscore.txt \
    --scores2 path/to/k562_cscore.txt \
    --name1 GM12878 --name2 K562 \
    -o results/expression_vs_switches.png
```

**B. 关联区室*状态*与*绝对*表达水平**:
```bash
evr-analyze-absolute \
    --rsem-files path/to/k562_rep1.tsv path/to/k562_rep2.tsv \
    --gtf path/to/gencode.v19.annotation.gtf \
    --bin-gene-map path/to/bin_gene_map.txt \
    --scores1 path/to/gm12878_cscore.txt \
    --scores2 path/to/k562_cscore.txt \
    --target-name K562 \
    -o results/absolute_expression_levels.png
```

### 步骤 3: 3D交互式可视化

这是所有分析的集大成者，在一个统一的、可交互的视图中呈现您的发现。

```bash
evr-plot3d-interactive \
    path/to/gm12878.pdb path/to/gm12878.cscore \
    path/to/k562.pdb path/to/k562.cscore \
    --name1 GM12878 --name2 K562 \
    --gtf path/to/gencode.v19.annotation.gtf \
    --bin-gene-map path/to/bin_gene_map.txt \
    --rna-seq path/to/rna_seq_de_results.csv \
    --rsem-files path/to/k562_rep1.tsv path/to/k562_rep2.tsv \
    --highlight-switches \
    --show-bins
```
**交互控制:**
- **鼠标左键**: 旋转
- **鼠标右键**: 缩放
- **鼠标中键**: 平移
- **底部滑块**: 在两种结构间平滑变形
- **在模型上左键单击**: 在终端中打印该bin的详细基因和表达信息

---

## 命令行工具参考

所有工具在安装后均可直接从命令行调用。使用 `-h` 或 `--help` 标志查看每个命令的详细选项。

**数据预处理 (`evr-prep-`)**
- `evr-hic2evr`: [一键式] 从.hic文件生成EVR密集矩阵。
- `evr-hic2cscore`: [一键式] 从.hic文件生成CscoreTool输入和bin定义文件。

**区室分数分析 (`evr-cscore-`)**
- `evr-cscore-threshold`: 分析C-score分布以辅助选择阈值。
- `evr-cscore-compare`: 量化比较两个样本的C-score。
- `evr-cscore-plot`: 绘制沿染色体的C-score曲线图。

**功能基因组学分析 (`evr-analyze-`)**
- `evr-analyze-change`: 绘制区室转换 vs 基因表达变化的箱线图。
- `evr-analyze-absolute`: 绘制区室状态 vs 基因绝对表达水平的箱线图。

**3D可视化 (`evr-plot3d-`)**
- `evr-plot3d-static`: 生成静态的3D结构图。
- `evr-plot3d-animation`: 生成结构变形的MP4动画。
- `evr-plot3d-interactive`: 启动功能完备的交互式3D查看器。
- `evr-plot3d-publication`: 生成带有Top N基因高亮的出版级静态图。

---

## 致谢与引用

本研究的实现得益于开源社区和公共数据集的巨大贡献。我们对所有相关的开发者、贡献者和数据提供方表示诚挚的感谢。

### 如何引用

如果您在您的研究中使用了本软件包，请考虑引用以下核心工具的原始文献。

### 数据来源
- **Hi-C Data**: Rao, S. S. P., et al. (2014). *Cell*, [GEO: GSE63525](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525).
- **RNA-seq Data**: The ENCODE Project Consortium. (2012). *Nature*. 数据从 [ENCODE Portal](https://www.encodeproject.org/) 获取，细胞系为 K562 ([ENCSR594NJP](https://www.encodeproject.org/experiments/ENCSR594NJP/)) 和 GM12878 ([ENCSR843RJV](https://www.encodeproject.org/experiments/ENCSR843RJV/))。
- **Gene Annotations**: [GENCODE](https://www.gencodegenes.org/) Release 19 (GRCh37/hg19).

### 核心分析软件
- **Juicer**: Durand, N. C., et al. (2016). *Cell Systems*. [GitHub: aidenlab/juicer](https://github.com/aidenlab/juicer).
- **CscoreTool**: Z. B. Xu, et al. (2021). *Bioinformatics*. [GitHub: scoutzxb/CscoreTool](https://github.com/scoutzxb/CscoreTool).
- **EVR**: G. Liverpool, et al. (2021). *Biophysical Journal*. [GitHub: mbglab/EVR](https://github.com/mbglab/EVR).