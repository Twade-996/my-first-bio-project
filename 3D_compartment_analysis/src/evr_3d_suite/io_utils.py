# 文件名: src/evr_3d_suite/io_utils.py
# 描述: 整合所有数据读取和解析功能，提供健壮的错误处理

import sys
import re
import pandas as pd
from functools import reduce

def read_pdb(filename):
    """
    解析PDB文件以提取每个原子的x, y, z坐标。
    """
    x, y, z = [], [], []
    try:
        with open(filename, 'r') as pdb_file:
            for line in pdb_file:
                if line.startswith("ATOM"):
                    x.append(float(line[30:38].strip()))
                    y.append(float(line[38:46].strip()))
                    z.append(float(line[46:54].strip()))
    except FileNotFoundError:
        print(f"错误: PDB文件 '{filename}' 未找到。", file=sys.stderr)
        sys.exit(1)
    except (ValueError, IndexError):
        print(f"错误: 无法从 '{filename}' 解析坐标。请检查PDB文件格式。", file=sys.stderr)
        sys.exit(1)
    if not x:
        print(f"警告: 在 '{filename}' 中未找到任何 'ATOM' 行。", file=sys.stderr)
    return x, y, z

def read_compartment_scores(filename):
    """
    读取一个分数文件，返回一个分数列表。
    智能处理两列格式、bedgraph格式，并跳过track行。
    """
    scores = []
    try:
        with open(filename, 'r') as score_file:
            for line in score_file:
                if line.lower().startswith('track'):
                    continue
                parts = line.strip().split()
                # 优先使用第4列（bedgraph），如果不存在则使用第2列
                score_val = None
                if len(parts) >= 4:
                    score_val = parts[3]
                elif len(parts) >= 2:
                    score_val = parts[1]
                
                if score_val is not None:
                    try:
                        scores.append(float(score_val))
                    except ValueError:
                        pass # 静默跳过无法转换的行
    except FileNotFoundError:
        print(f"错误: 得分文件 '{filename}' 未找到。", file=sys.stderr)
        sys.exit(1)
    return scores

def parse_bin_gene_map(map_file):
    """
    解析bedtools的输出，创建bin->[genes]映射和gene->bin映射。
    """
    bin_to_gene = {}
    gene_to_bin = {}
    print(f"正在从 '{map_file}' 解析bin-基因映射...")
    try:
        with open(map_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 13: continue
                
                # 修正: 正确的列索引
                feature_type = parts[6]
                if feature_type != 'gene':
                    continue
                
                bin_id_str = parts[3]
                try:
                    bin_index = int(bin_id_str.split('_')[1]) - 1
                except (IndexError, ValueError):
                    continue
                
                attributes_str = parts[12]
                match = re.search(r'gene_name "([^"]+)"', attributes_str)
                if match:
                    gene_name = match.group(1)
                    if bin_index not in bin_to_gene:
                        bin_to_gene[bin_index] = set()
                    bin_to_gene[bin_index].add(gene_name)
                    if gene_name not in gene_to_bin:
                        gene_to_bin[gene_name] = bin_index
                    # 改进: 处理一个基因跨越多个bin的情况
                    elif gene_to_bin[gene_name] != bin_index:
                        # print(f"警告: 基因 '{gene_name}' 出现在多个bins中。反向映射只保留第一次出现的位置。", file=sys.stderr)
                        pass

    except FileNotFoundError:
        print(f"错误: bin-基因映射文件 '{map_file}' 未找到。", file=sys.stderr)
        sys.exit(1)
    
    for bin_index in bin_to_gene:
        bin_to_gene[bin_index] = sorted(list(bin_to_gene[bin_index]))
        
    print(f"完成。在 {len(bin_to_gene)} 个bins中找到了 {len(gene_to_bin)} 个唯一基因。")
    return bin_to_gene, gene_to_bin

def create_id_map_from_gtf(gtf_file_path):
    """
    解析GTF文件，创建一个从ENSEMBL ID (去除版本号) 到基因名的映射。
    """
    print(f"正在从 '{gtf_file_path}' 构建ENSEMBL ID到基因名的映射...")
    id_map = {}
    try:
        with open(gtf_file_path, 'r') as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')
                
                # 修正: 正确的列索引
                if len(parts) < 9 or parts[2] != 'gene':
                    continue
                
                attributes_str = parts[8]
                gene_id_match = re.search(r'gene_id "([^"]+)"', attributes_str)
                gene_name_match = re.search(r'gene_name "([^"]+)"', attributes_str)
                if gene_id_match and gene_name_match:
                    # 修正: 正确地从字符串中提取ID
                    ensembl_id_no_version = gene_id_match.group(1).split('.')[0]
                    id_map[ensembl_id_no_version] = gene_name_match.group(1)
    except FileNotFoundError:
        print(f"错误: GTF文件 '{gtf_file_path}' 未找到。", file=sys.stderr)
        sys.exit(1)
    print(f"成功创建映射，包含 {len(id_map)} 个条目。")
    return id_map

def process_rsem_to_tpm_map(rsem_files, id_to_name_map):
    """
    读取RSEM文件列表，转换ID，合并，并返回一个 gene_name -> mean_TPM 的字典。
    """
    all_dataframes = []
    print("\n--- 正在处理RSEM定量文件 ---")
    for i, file_path in enumerate(rsem_files):
        print(f"  - 读取文件: {file_path}")
        try:
            df = pd.read_csv(file_path, sep='\t')
            if 'gene_id' not in df.columns:
                print(f"  警告: 文件 '{file_path}' 缺少 'gene_id' 列。已跳过。", file=sys.stderr)
                continue
            
            # 修正: 正确的语法
            df['gene_id_no_version'] = df['gene_id'].str.split('.').str[0]
            df['gene_name'] = df['gene_id_no_version'].map(id_to_name_map)
            df_processed = df[['gene_name', 'TPM']].copy().rename(columns={'TPM': f'TPM_rep{i+1}'})
            df_processed.dropna(subset=['gene_name'], inplace=True)
            all_dataframes.append(df_processed)
        except Exception as e:
            print(f"  警告: 无法处理文件 '{file_path}'。错误: {e}。已跳过。", file=sys.stderr)
    
    if not all_dataframes:
        print("错误: 没有处理任何有效的RSEM文件。将返回一个空的TPM映射。", file=sys.stderr)
        return {}
        
    print("  - 正在合并重复样本并计算平均TPM...")
    df_merged = reduce(lambda left, right: pd.merge(left, right, on='gene_name', how='outer'), all_dataframes)
    df_merged.set_index('gene_name', inplace=True)
    # 改进: 填充NaN为0
    mean_tpm_series = df_merged.fillna(0).mean(axis=1)
    print("  - 平均TPM计算完成。")
    return mean_tpm_series.to_dict()

def load_de_results(filepath):
    """
    从CSV或TSV文件加载差异表达结果。
    """
    try:
        df = pd.read_csv(filepath, sep=None, engine='python')
        if 'gene_name' in df.columns and 'log2FoldChange' in df.columns:
            # 改进: 去除重复的基因名
            df.drop_duplicates(subset=['gene_name'], keep='first', inplace=True)
            return pd.Series(df.log2FoldChange.values, index=df.gene_name).to_dict()
        else:
            print(f"警告: '{filepath}' 缺少 'gene_name' 或 'log2FoldChange' 列。", file=sys.stderr)
            return {}
    except Exception as e:
        print(f"警告: 无法加载差异表达结果文件 '{filepath}'。错误: {e}", file=sys.stderr)
        return {}