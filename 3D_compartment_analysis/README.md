# EVR 3D Suite: A Chromosome 3D Analysis & Visualization Toolkit

**EVR 3D Suite** is a comprehensive Python package for analyzing 3D chromosome structures from Hi-C data. It provides a suite of command-line tools to process compartment scores, integrate gene expression data, and generate publication-quality static, animated, and fully interactive 3D visualizations of chromosome conformations.

This project integrates several key stages of 3D genomics analysis:
- **Compartment Analysis**: Analyze, compare, and visualize A/B compartment scores.
- **Functional Genomics Integration**: Correlate compartment dynamics with RNA-seq data (differential expression and absolute levels).
- **Advanced 3D Visualization**: Create dynamic and interactive 3D models that morph between cellular states, highlighting structural and functional changes.

---

## Table of Contents

- [Installation](#installation)
- [Analysis Workflow](#analysis-workflow)
  - [Prerequisites: Getting 3D Structures and Scores](#prerequisites)
  - [Step 1: Analyze Compartment Scores](#step-1-analyze-compartment-scores)
  - [Step 2: Correlate Structure with Gene Expression](#step-2-correlate-structure-with-gene-expression)
  - [Step 3: Visualize in 3D](#step-3-visualize-in-3d)
- [Command-Line Interface (CLI) Reference](#cli-reference)
- [Acknowledgments](#acknowledgments)

---

## Installation

It is highly recommended to install the package in a dedicated Conda environment.

1.  **Create and activate a new Conda environment:**
    ```bash
    conda create -n evr_env python=3.9
    conda activate evr_env
    ```

2.  **Install FFmpeg:**
    This package requires the external tool FFmpeg for creating animations. You can install it via Conda:
    ```bash
    conda install -c conda-forge ffmpeg
    ```

3.  **Clone and install the package:**
    Clone this repository and install the `evr_3d_suite` package in editable mode. This allows you to modify the source code and have the changes immediately available.
    ```bash
    git clone https://your-repo-url/3D_compartment_analysis.git
    cd 3D_compartment_analysis
    pip install -e .
    ```
    After installation, all command-line tools (e.g., `evr-cscore-compare`, `evr-plot3d-interactive`) will be available in your terminal.

---

## Analysis Workflow

This section outlines a typical analysis workflow using the tools provided in this package.

### Step 0: Preprocessing Juicer Output

The first step is to convert the raw text output from `juicer_tools dump` into formats compatible with downstream tools.

**A. For EVR 3D modeling:**
Convert the 3-column KR-normalized matrix into a dense matrix required by EVR.
```bash
evr-juicer-to-evr path/to/juicer_kr.txt path/to/evr_matrix.txt
```

**B. For CscoreTool analysis:**
Convert the 3-column raw interaction matrix into the 7-column summary format required by CscoreTool. This command also generates a corresponding BED file defining the genomic coordinates for each bin, which is essential for downstream analysis.
```bash
evr-juicer-to-cscore path/to/juicer_raw.txt path/to/cscore_input.summary \
    --chrom chr1 \
    --binsize 100000 \
    --output-bed path/to/chr1_100kb_bins.bed
```

### Prerequisites: Getting 3D Structures and Scores

The `evr_3d_suite` package focuses on the analysis and visualization *after* you have generated the initial 3D structures and compartment scores. You can generate these input files using standard tools:

-   **3D Structure (`.pdb`)**: Use **EVR** to reconstruct the 3D coordinates from a Hi-C contact matrix.
-   **Compartment Scores (`.cscore.txt` or `.bedgraph`)**: Use tools like **CscoreTool** or **HiCExplorer's `hicPCA`** to calculate the first principal component (PC1), which represents the A/B compartment score.

### Step 1: Analyze Compartment Scores

Before diving into 3D, it's crucial to understand your compartment data.

**A. Find a suitable C-score threshold:**
```bash
evr-cscore-threshold path/to/your/cscore.txt
```
This will print statistical summaries and show a histogram to help you decide on a threshold for defining A/B compartments (e.g., `-t 0.1`).

**B. Compare scores between two conditions:**
```bash
evr-cscore-compare \
    path/to/gm12878_cscore.bedgraph \
    path/to/k562_cscore.bedgraph \
    --name1 GM12878 \
    --name2 K562 \
    --outdir results/cscore_comparison
```
This generates a comprehensive set of plots and statistics comparing the two compartment profiles.

### Step 2: Correlate Structure with Gene Expression

This is the core functional analysis. You'll need your RNA-seq results and a bin-to-gene map.

**A. Correlate compartment *switches* with expression *changes***:
```bash
evr-analyze-change \
    --bin-gene-map path/to/bin_gene_map.txt \
    --rna-seq path/to/rna_seq_de_results.csv \
    --scores1 path/to/gm12878_cscore.txt \
    --scores2 path/to/k562_cscore.txt \
    --name1 GM12878 --name2 K562 \
    -o results/expression_vs_switches.png
```

**B. Correlate compartment *status* with *absolute* expression levels**:
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

### Step 3: Visualize in 3D

This is the final step where all data comes together in a powerful, interactive visualization.

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
**Interactive Controls:**
- **Left Mouse**: Rotate
- **Right Mouse**: Zoom
- **Middle Mouse**: Pan
- **Slider**: Morph between structures
- **Left-Click on Bin**: Print detailed gene and expression info to the console.

---

## Command-Line Interface (CLI) Reference

All tools are available directly from the command line after installation. Use the `-h` or `--help` flag for detailed options for each command.

-   **`evr-cscore-threshold`**: Analyzes C-score distribution.
-   **`evr-cscore-compare`**: Compares C-scores between two samples.
-   **`evr-cscore-plot`**: Plots a C-score profile along a chromosome.
-   **`evr-analyze-change`**: Plots compartment switches vs. log2FoldChange.
-   **`evr-analyze-absolute`**: Plots compartment status vs. absolute expression (log2(TPM+1)).
-   **`evr-plot3d-static`**: Generates a static 3D image.
-   **`evr-plot3d-interactive`**: Launches the main interactive 3D viewer.

---

## Acknowledgments

This research was made possible by the excellent open-source software and public datasets provided by the scientific community. We express our sincere gratitude to all developers, contributors, and data providers.

### Data Sources
- **Hi-C Data**: Rao, S. S. P., et al. (2014). *Cell*, [GEO: GSE63525](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525).
- **RNA-seq Data**: The ENCODE Project Consortium. (2012). *Nature*. Data was sourced from the [ENCODE Portal](https://www.encodeproject.org/) for K562 ([ENCSR594NJP](https://www.encodeproject.org/experiments/ENCSR594NJP/)) and GM12878 ([ENCSR843RJV](https://www.encodeproject.org/experiments/ENCSR843RJV/)) cell lines.
- **Gene Annotations**: [GENCODE](https://www.gencodegenes.org/) Release 19 (GRCh37/hg19).

### Core Analysis Software
- **Juicer**: [aidenlab/juicer](https://github.com/aidenlab/juicer)
- **HiCExplorer**: [deeptools/HiCExplorer](https://github.com/deeptools/HiCExplorer)
- **CscoreTool**: [scoutzxb/CscoreTool](https://github.com/scoutzxb/CscoreTool)
- **EVR**: [mbglab/EVR](https://github.com/mbglab/EVR)