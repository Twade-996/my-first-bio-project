# src/evr_3d_suite/__init__.py

"""
EVR 3D Suite: A comprehensive package for 3D chromosome structure analysis, 
visualization, and multi-omics integration.
"""

# 1. 定义包的版本号
# 当您更新包时，只需要修改这一个地方
__version__ = "0.1.0"

# 2. 提升核心功能函数到包的顶层命名空间
# 这让用户可以更方便地导入和使用您包中的关键功能。

# --- 从 io_utils.py 提升 ---
# 这些是几乎所有分析都会用到的基础工具，非常适合提升。
from .io_utils import (
    read_pdb,
    read_compartment_scores,
    parse_bin_gene_map,
    create_id_map_from_gtf,
    process_rsem_to_tpm_map,
    load_de_results
)

# --- 从 cscore.py 提升 ---
# 提升核心的计算逻辑函数，而不是命令行入口函数。
from .cscore import (
    analyze_thresholds_logic,
    compare_scores_logic
)

# --- 从 plot3d.py 提升 ---
# 提升核心的辅助函数。
from .plot3d import (
    center_coords,
    calculate_rotation_matrix,
    calculate_rmsd
)

# 3. 定义 __all__ (可选但推荐)
# 这明确地告诉Python，当其他代码使用 `from evr_3d_suite import *` 时，
# 应该导入哪些名字。这是一种良好的编程实践。
__all__ = [
    # 版本号
    "__version__",
    
    # io_utils
    "read_pdb",
    "read_compartment_scores",
    "parse_bin_gene_map",
    "create_id_map_from_gtf",
    "process_rsem_to_tpm_map",
    "load_de_results",
    
    # cscore
    "analyze_thresholds_logic",
    "compare_scores_logic",

    # plot3d
    "center_coords",
    "calculate_rotation_matrix",
    "calculate_rmsd"
]