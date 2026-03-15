#!/usr/bin/env python3
"""
脚本 03: 将 BLAST 命中映射回基因存在/缺失矩阵

功能:
- 读取 LDH 和 MCM 的 BLAST 结果
- 从 pan_genome_reference.fa 中提取基因 ID 到物种的映射
- 生成物种 × 代谢酶的 0/1 矩阵
"""

import pandas as pd
import re
from pathlib import Path


def extract_gene_to_species_mapping(pan_genome_fa: Path) -> dict:
    """
    从 pan_genome_reference.fa 中提取基因 ID 到物种的映射
    
    假设基因 ID 格式：gene_id|species_strain
    返回：{gene_id: species_strain}
    """
    mapping = {}
    
    if not pan_genome_fa.exists():
        print(f"警告：未找到 pan_genome_reference.fa: {pan_genome_fa}")
        return mapping
    
    with open(pan_genome_fa, 'r', encoding='utf-8') as f:
        current_gene = None
        for line in f:
            if line.startswith('>'):
                # 提取基因 ID
                header = line[1:].strip()
                # 假设格式：>gene_id|species_strain 或 >species_strain|gene_id
                current_gene = header.split('|')[0].strip()
                
                # 尝试从 header 中提取物种信息
                if '|' in header:
                    parts = header.split('|')
                    if len(parts) >= 2:
                        # 假设第二部分包含物种信息
                        species_info = parts[1].strip()
                        mapping[current_gene] = species_info
    return mapping


def parse_blast_results(blast_file: Path) -> set:
    """
    解析 BLAST 结果，返回命中的基因 ID 集合
    """
    if not blast_file.exists():
        print(f"警告：未找到 BLAST 结果文件：{blast_file}")
        return set()
    
    try:
        # BLAST outfmt 6: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs
        df = pd.read_csv(
            blast_file,
            sep='\t',
            header=None,
            names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                   'gapopen', 'qstart', 'qend', 'sstart', 'send', 
                   'evalue', 'bitscore', 'qcovs']
        )
        
        # 返回命中的目标基因 ID（sseqid）
        return set(df['sseqid'].unique())
    
    except Exception as e:
        print(f"错误：解析 BLAST 文件失败 {blast_file}: {e}")
        return set()


def map_to_species(ldh_genes: set, mcm_genes: set, gene_mapping: dict) -> pd.DataFrame:
    """
    将基因命中映射到物种水平
    
    返回：DataFrame (species, LDH, MCM)
    """
    # 收集所有物种
    all_species = set()
    
    # 从基因映射中获取所有物种
    for gene_id, species in gene_mapping.items():
        # 提取物种名（去除菌株信息）
        # 假设格式：Veillonella_parvula_strain_X
        parts = species.split('_')
        if len(parts) >= 3:
            species_name = '_'.join(parts[:3])  # Veillonella_parvula
        else:
            species_name = species
        all_species.add(species_name)
    
    # 初始化结果 DataFrame
    result = pd.DataFrame({
        'species': sorted(all_species),
        'LDH': 0,
        'MCM': 0
    })
    result = result.set_index('species')
    
    # 标记 LDH 阳性物种
    for gene_id in ldh_genes:
        if gene_id in gene_mapping:
            species_info = gene_mapping[gene_id]
            parts = species_info.split('_')
            if len(parts) >= 3:
                species_name = '_'.join(parts[:3])
            else:
                species_name = species_info
            
            if species_name in result.index:
                result.loc[species_name, 'LDH'] = 1
    
    # 标记 MCM 阳性物种
    for gene_id in mcm_genes:
        if gene_id in gene_mapping:
            species_info = gene_mapping[gene_id]
            parts = species_info.split('_')
            if len(parts) >= 3:
                species_name = '_'.join(parts[:3])
            else:
                species_name = species_info
            
            if species_name in result.index:
                result.loc[species_name, 'MCM'] = 1
    
    return result.reset_index()


def main():
    print("=== 生成代谢基因存在/缺失矩阵 ===\n")
    
    # 输入文件
    ldh_blast = Path('results/metabolism/ldh_blast.tsv')
    mcm_blast = Path('results/metabolism/mcm_blast.tsv')
    pan_genome_fa = Path('results/pangenome/pan_genome_reference.fa')
    output_file = Path('results/metabolism/metabolism_matrix.csv')
    
    # 1. 提取基因到物种的映射
    print("1. 从 pan_genome_reference.fa 提取基因映射...")
    gene_mapping = extract_gene_to_species_mapping(pan_genome_fa)
    print(f"   找到 {len(gene_mapping)} 个基因的映射关系")
    
    # 2. 解析 BLAST 结果
    print("\n2. 解析 BLAST 结果...")
    ldh_genes = parse_blast_results(ldh_blast)
    mcm_genes = parse_blast_results(mcm_blast)
    print(f"   LDH 命中 {len(ldh_genes)} 个基因")
    print(f"   MCM 命中 {len(mcm_genes)} 个基因")
    
    # 3. 映射到物种水平
    print("\n3. 映射到物种水平...")
    result_df = map_to_species(ldh_genes, mcm_genes, gene_mapping)
    
    # 4. 保存结果
    output_file.parent.mkdir(parents=True, exist_ok=True)
    result_df.to_csv(output_file, index=False)
    print(f"\n✓ 保存代谢矩阵到 {output_file}")
    
    # 打印摘要
    print("\n=== 代谢基因分布摘要 ===")
    print(f"总物种数：{len(result_df)}")
    print(f"LDH 阳性：{result_df['LDH'].sum()}")
    print(f"MCM 阳性：{result_df['MCM'].sum()}")
    print(f"LDH+MCM 双阳性：{len(result_df[(result_df['LDH']==1) & (result_df['MCM']==1)])}")
    
    print("\n=== 完成 ===")


if __name__ == '__main__':
    main()
