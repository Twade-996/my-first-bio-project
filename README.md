# My First Bio Project
## HiC_chromosome_evr_3D
From hic data to chromosome 3D dynamic models

----------------------------------------------------------------------------
# 使用EVR计算bin坐标

.hic

java -jar juicer_tools.jar dump observed KR GSE63525_K562_combined_30.hic chr1 chr1 BP 100000 chr1_100kb.KRnorm.txt

chr1_100kb.KRnorm.txt

convert_3col_to_evr.py chr1_100kb.KRnorm.txt chr1_100kb_evr_matrix.txt

chr1_100kb_evr_matrix.txt

python evr.py -i evr_matrix.txt -o evr_structure.pdb --se 6 -t 6 --seed 1 --ChrType 1

Chr1_100kb_evr_structure.pdb

----------------------------------------------------------------------------
# 使用CscoreTool划分区室

.hic

java -Xmx16G -jar juicer_tools.jar dump observed NONE GSE63525_K562_combined_30.hic chr1 chr1 BP 100000 chr1_100kb_raw.txt

python format_raw_hic.py chr1_100kb_raw.txt chr1_100kb.summary --chrom chr1 --binsize 100000 --output_bed chr1_100kb_bins.bed

chr1_100kb.summary & chr1_100kb_bins.bed

CscoreTool1.1 < windows.bed> < input.summary> < outputPrefix> < session> < minDis> [chrName]

Cscore(txt/bedgraph/…)

----------------------------------------------------------------------------
# 使用PCA验证区室划分

.hic

hicConvertFormat --matrices GSE63525_GM12878_insitu_DpnII_combined_30.hic --outFileName GM12878_50kb.cool --inputFormat hic --outputFormat cool --resolutions 50000

hicConvertFormat --matrices GM12878_50kb_50000.cool --outFileName GM12878_50kb.h5 --inputFormat cool --outputFormat h5 --chromosome 1 --resolutions 50000

hicPCA --matrix GM12878_50kb.h5 --outputFileName gm12878_chr1_50kb_pca.bedgraph --numberOfEigenvectors 1 --format bedgraph --chromosomes 1

---------------------------------------------------------------------------
# 判断A/B区室的cscore阈值

cscore.txt

python cscores_threshold.py cscore.txt

suggested thresholds

---------------------------------------------------------------------------
# 对比分析cscore,并汇总图像

sample_1.bedgraph & sample_2.bedgraph

compare_cscores.py [-h] --outdir OUTDIR [--name1 NAME1] [--name2 NAME2] cscore1 cscore2

boundary_drift.png
compartment_switch_rate.png
cscore_compare.txt
cscore_comparison_summary.csv
cscore_delta_bar.png
cscore_distribution.png
cscore_heatmap.png
cscore_scatter.png
pearson.txt
spearman.txt

python auto_plot_images.py /your/output/path

.png

----------------------------------------------------------------------------
## 将IO.py放置到与绘图脚本同一目录中
# 标注区室的染色体3D模型

sample_1.pdb & .bedgraph

python plot_with_compartment.py chr1_100kb_evr_structure.pdb -cf chr1_100kb__cscore.txt --show

----------------------------------------------------------------------------
# 可交互、可平滑变化的一对3D染色体模型

sample_1.pdb/cscore.txt & sample_2.pdb/cscore.txt

python plot_interactive_morph.py sample_1.pdb cscore.txt sample_2.pdb cscore.txt [...]

----------------------------------------------------------------------------
# 平滑变化的一对3D染色体模型动画

sample_1.pdb/cscore.txt & sample_2.pdb/cscore.txt

python plot_animation.py sample_1.pdb cscore.txt sample_2.pdb cscore.txt -o .mp4

.mp4

-------------------------------------------------------------------------------
# 创建bin与基因的映射关系，并检验区室变化和基因表达之间的关系

awk '{print $0 "\tbin_" FNR}' chr1_100kb_bins.bed > chr1_100kb_bins_named.bed

bedtools intersect -a chr1_100kb_bins_named.bed -b gencode.v19.annotation.gtf -wa -wb > chr1_100kb_bin_gene_map.txt

python gene_and_compartment_analysis.py \
    --bin-gene-map bin_gene_map.txt \
    --rna-seq rna_seq_de_results.csv \
    --scores1 /path/to/your/GM_chr1_100kb__cscore.txt \
    --scores2 /path/to/your/k562_chr1_100kb__cscore.txt \
    --name1 GM12878 \
    --name2 K562 \
    -t 0.1 \
    -o gm12878_vs_k562_compartment_expression.png

# 检验基因绝对表达量的差异

python gene_abs_expression.py \
    --rsem-files /path/to/K562_rep1.tsv /path/to/K562_rep2.tsv \
    --gtf gencode.v19.annotation.gtf \
    --bin-gene-map bin_gene_map_50kb.txt \
    --scores1 /path/to/GM_chr1_50kb__cscore.txt \
    --scores2 /path/to/k562_chr1_50kb__cscore.txt \
    --target-name K562 \
    --skip-range 2430 2578 \
    -t 0.1 \
    -o k562_absolute_expression_levels_50kb_final.png

--------------------------------------------------------------------------------
# 构建区室B->A且含基因的模型

# 标注所有区室变化的bin，点击可输出高亮bin中所含基因及log2FC、TPM
python new_plot.py evr_structure.pdb cscore.txt evr_structure.pdb cscore.txt --gtf gencode.v19.annotation.gtf --bin-gene-map chr1_100kb_bin_gene_map.txt --rna-seq rna_seq_de_results.csv --rsem-files k562_rep1.tsv k562_rep2.tsv -t 0.2 --highlight-switches

# 标注B->A且含有基因的bin
python plot_showcase.py \
    cscore.txt \
    cscore.txt \
    --name1 GM12878 --name2 K562 \
    --gtf gencode.v19.annotation.gtf \
    --bin-gene-map bin_gene_map_50kb.txt \
    --rna-seq rna_seq_de_results.csv \
    --skip-range 2430 2578 \
    -t 0.1 \
    --highlight-type "B -> A"

# 绘制带有区室变化的3D模型，并标出含基因表达变化量最大的前n个bin

python all_py/plot_pic_gene.py \
evr_structure.pdb cscore1.txt cscore2.txt \
--bin-gene-map chr1_100kb_bin_gene_map.txt \
-t 0.2 \
--show-bins \
--de-results rna_seq_de_results.csv \
--top-genes 20 \
--bin-size 5
-o gene_top_20.png

--------------------------------------------------------------------------------
## 致谢 (Acknowledgments)

本研究的实现离不开以下优秀的开源软件和公开数据集的支持。我们对所有开发者、贡献者及数据提供方表示最诚挚的感谢。

### 数据来源 (Data Sources)

#### Hi-C 数据

*   **NCBI GEO: GSE63525**: 本项目所使用的 Hi-C 数据集来源于 Rao 等人于 2014 年发表在 *Cell* 上的研究。此数据集是本项目分析工作的基石。
    *   **访问链接**: [GEO Accession GSE63525](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525)
    *   **主要出版物**: Rao, S. S. P., Huntley, M. H., Durand, N. C., et al. (2014). A 3D map of the human genome at kilobase resolution reveals principles of chromatin looping. *Cell*, *159*(7), 1665-1680.

#### RNA-seq 数据

*   **ENCODE Project**: 本项目使用了来自 ENCODE (The Encyclopedia of DNA Elements) 项目的 RNA-seq 数据。我们感谢 ENCODE 联盟以及数据生成实验室 (Thomas Gingeras, CSHL) 的贡献。
    *   **项目引用**: The ENCODE Project Consortium. (2012). An integrated encyclopedia of DNA elements in the human genome. *Nature*, *489*(7414), 57–74.
    *   **使用的数据集**:
        *   **ENCSR594NJP**: K562 细胞系的 RNA-seq 数据。 [访问链接](https://www.encodeproject.org/experiments/ENCSR594NJP/)
        *   **ENCSR843RJV**: GM12878 细胞系的 RNA-seq 数据。 [访问链接](https://www.encodeproject.org/experiments/ENCSR843RJV/)

### 核心分析软件 (Core Analysis Software)

*   **Juicer**: 由 Aiden Lab 开发的一站式 Hi-C 数据处理流程。
    *   **许可证**: MIT License (商业用途需额外许可)
    *   **GitHub**: [aidenlab/juicer](https://github.com/aidenlab/juicer)
    *   **引用**: Durand, N. C., et al. (2016). Juicer provides a one-click system for analyzing loop-resolution Hi-C experiments. *Cell systems*, *3*(1), 95-98.

*   **HiCExplorer**: 一套功能强大的 Hi-C 数据分析和可视化工具集。
    *   **许可证**: GPL-3.0 License
    *   **GitHub**: [deeptools/HiCExplorer](https://github.com/deeptools/HiCExplorer)
    *   **引用**: Ramirez, F., et al. (2018). Galaxy HiCExplorer: a web server for reproducible Hi-C data analysis, quality control and visualization. *Nucleic acids research*, *46*(W1), W11-W16.

*   **CscoreTool**: 用于 Hi-C 数据 A/B 区室分析的工具。
    *   **许可证**: 未指定 (No license specified)
    *   **GitHub**: [scoutzxb/CscoreTool](https://github.com/scoutzxb/CscoreTool)
    *   **引用**: Zheng, X., & Zheng, Y. (2018). CscoreTool: fast Hi-C compartment analysis at high resolution. *Bioinformatics*, *34*(9), 1568–1570.

*   **EVR (Error-Vector Resultant)**: 用于从互作数据中进行染色质三维结构重建的工具。
    *   **许可证**: MIT License
    *   **GitHub**: [mbglab/EVR](https://github.com/mbglab/EVR)
    *   **引用**: Hua, K.-J., & Ma, B.-G. (2019). EVR: reconstruction of bacterial chromosome 3D structure models using error-vector resultant algorithm. *BMC Genomics*, *20*(1), 738.