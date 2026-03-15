#!/usr/bin/env python3
"""
脚本 01: 从 NCBI datasets 解压产物中提取并重命名基因组

功能:
- 遍历临时目录中的所有 .fna 文件
- 从 assembly_data_report.json 读取物种名和菌株名
- 重命名为统一格式：Veillonella_species_strain.fna
- 如果是外群，命名为 outgroup.fna
"""

import json
import shutil
from pathlib import Path
import sys


def parse_assembly_report(report_path: Path) -> dict:
    """解析 assembly_data_report.json，返回 {accession: {species, strain}}"""
    with open(report_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    result = {}
    for assembly in data.get('assemblies', []):
        # 提取关键信息
        info = assembly.get('assembly_info', {})
        accession = info.get('accession', '')
        
        # 物种名称（属 + 种）
        species_name = info.get('species_name', '')
        
        # 菌株名
        strain = info.get('strain', '')
        
        # RefSeq 类别
        refseq_category = info.get('refseq_category', '')
        
        if accession and species_name:
            result[accession] = {
                'species': species_name.replace(' ', '_'),
                'strain': strain.replace(' ', '_').replace(',', '') if strain else 'unknown',
                'category': refseq_category
            }
    
    return result


def rename_genomes(input_dir: Path, output_dir: Path):
    """重命名基因组文件"""
    report_path = input_dir / 'assembly_data_report.json'
    
    if not report_path.exists():
        print(f"错误：未找到 assembly_data_report.json: {report_path}")
        sys.exit(1)
    
    # 解析报告
    assembly_info = parse_assembly_report(report_path)
    print(f"找到 {len(assembly_info)} 个基因组")
    
    # 确保输出目录存在
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 遍历所有 .fna 文件
    renamed_count = 0
    for fna_file in input_dir.rglob('*.fna'):
        # 从路径中提取 accession
        # 典型路径：ncbi_dataset/data/GCF_000001/GCF_000001.fna
        parts = fna_file.parent.parts
        accession = None
        for part in parts:
            if part.startswith('GCF_') or part.startswith('GCA_'):
                accession = part
                break
        
        if not accession or accession not in assembly_info:
            print(f"跳过未知基因组：{fna_file}")
            continue
        
        info = assembly_info[accession]
        
        # 检查是否为外群（Megasphaera）
        if 'Megasphaera' in info['species']:
            new_name = 'outgroup.fna'
        else:
            # 韦荣氏菌：Veillonella_种名_菌株名.fna
            species = info['species'].replace('Veillonella_', '')
            strain = info['strain'] if info['strain'] != 'unknown' else accession
            new_name = f"Veillonella_{species}_{strain}.fna"
        
        # 复制并重命名
        output_path = output_dir / new_name
        shutil.copy2(fna_file, output_path)
        print(f"✓ {accession} → {new_name}")
        renamed_count += 1
    
    print(f"\n完成：重命名 {renamed_count} 个基因组文件")


def main():
    if len(sys.argv) != 3:
        print("用法：python 01_download_genomes.py <input_dir> <output_dir>")
        print("  input_dir: NCBI datasets 解压后的目录")
        print("  output_dir: 输出目录（存放重命名后的 .fna 文件）")
        sys.exit(1)
    
    input_dir = Path(sys.argv[1])
    output_dir = Path(sys.argv[2])
    
    if not input_dir.exists():
        print(f"错误：输入目录不存在：{input_dir}")
        sys.exit(1)
    
    rename_genomes(input_dir, output_dir)


if __name__ == '__main__':
    main()
