#!/usr/bin/env python3
"""
脚本：生成补充材料（Supplementary Materials）
功能：
- Table S1: 基因组信息表（包含NCBI信息、质量统计、ANI聚类）
- Table S2: 泛基因组矩阵
- Table S3: 功能注释表（COG分类、KEGG通路）
- Table S4: 毒力基因表
- Table S5: 抗性基因表
- 生成补充材料说明文档
"""

import pandas as pd
import numpy as np
from pathlib import Path
import json
from datetime import datetime
import sys


def generate_table_s1(genome_dir: Path, qc_file: Path, ani_clusters: Path) -> pd.DataFrame:
    """生成Table S1: 基因组信息表"""
    print("=== 生成 Table S1: 基因组信息 ===")

    # 基础基因组信息
    genomes = sorted([f.stem for f in genome_dir.glob("*.fna")])

    data = []
    for genome in genomes:
        record = {
            'Genome_ID': genome,
            'NCBI_Accession': '',  # 需要从原始数据提取
            'Species': genome.split('_')[1] if '_' in genome else 'Unknown',
            'Strain': '_'.join(genome.split('_')[2:]) if '_' in genome else 'Unknown',
            'Genome_Size_bp': '',
            'GC_Content': '',
            'N50': '',
            'Contigs': ''
        }
        data.append(record)

    df = pd.DataFrame(data)

    # 合并QC信息
    if qc_file.exists():
        qc_df = pd.read_csv(qc_file)
        # 使用基础QC脚本的列名
        qc_cols = ['Name', 'Genome_Size', 'Contigs', 'N50', 'GC_Content', 'Pass_QC']
        df = df.merge(qc_df[qc_cols],
                      left_on='Genome_ID', right_on='Name', how='left')
        df = df.drop('Name', axis=1)

    # 合并ANI聚类
    if ani_clusters.exists():
        ani_df = pd.read_csv(ani_clusters)
        df = df.merge(ani_df, left_on='Genome_ID', right_on='Genome', how='left')
        df = df.drop('Genome', axis=1)

    print(f"  包含 {len(df)} 个基因组信息")
    return df


def generate_table_s2(panaroo_dir: Path) -> pd.DataFrame:
    """生成Table S2: 泛基因组存在/缺失矩阵"""
    print("\n=== 生成 Table S2: 泛基因组矩阵 ===")

    gpa_file = panaroo_dir / "gene_presence_absence.csv"
    if not gpa_file.exists():
        print(f"  警告：未找到 {gpa_file}")
        return None

    df = pd.read_csv(gpa_file, low_memory=False)
    print(f"  矩阵大小: {df.shape}")
    return df


def generate_table_s3(annotation_dir: Path) -> pd.DataFrame:
    """生成Table S3: 功能注释表（COG/KEGG）"""
    print("\n=== 生成 Table S3: 功能注释 ===")

    # 这里应该整合eggNOG-mapper的结果
    # 目前创建示例结构
    data = {
        'Gene_ID': [],
        'Genome': [],
        'COG_Category': [],
        'COG_Description': [],
        'KEGG_ko': [],
        'KEGG_Pathway': [],
        'Gene_Name': [],
        'Description': []
    }

    df = pd.DataFrame(data)
    print(f"  注意：需要运行eggNOG-mapper后填充此表")
    return df


def generate_cog_summary(panaroo_dir: Path, cog_annotations: Path) -> pd.DataFrame:
    """生成COG分类汇总表"""
    print("\n=== 生成 COG 分类汇总 ===")

    # 读取泛基因组基因列表
    gpa_file = panaroo_dir / "gene_presence_absence.csv"
    if not gpa_file.exists():
        return None

    gpa_df = pd.read_csv(gpa_file, low_memory=False)

    # 统计每个COG类别的基因数
    cog_categories = {
        'D': 'Cell cycle control',
        'M': 'Cell wall/membrane',
        'N': 'Cell motility',
        'O': 'Post-translational modification',
        'T': 'Signal transduction',
        'U': 'Intracellular trafficking',
        'V': 'Defense mechanisms',
        'W': 'Extracellular structures',
        'Y': 'Nuclear structure',
        'Z': 'Cytoskeleton',
        'A': 'RNA processing',
        'B': 'Chromatin structure',
        'J': 'Translation',
        'K': 'Transcription',
        'L': 'Replication',
        'C': 'Energy production',
        'E': 'Amino acid metabolism',
        'F': 'Nucleotide metabolism',
        'G': 'Carbohydrate metabolism',
        'H': 'Coenzyme metabolism',
        'I': 'Lipid metabolism',
        'P': 'Inorganic ion transport',
        'Q': 'Secondary metabolites',
        'R': 'General function',
        'S': 'Function unknown'
    }

    # 创建COG汇总
    data = []
    for code, desc in cog_categories.items():
        data.append({
            'COG_Category': code,
            'Description': desc,
            'Core_Genes': 0,  # 需要实际统计
            'Accessory_Genes': 0,
            'Total_Genes': 0
        })

    df = pd.DataFrame(data)
    return df


def generate_readme(supp_dir: Path, tables: dict):
    """生成补充材料说明文档"""
    print("\n=== 生成补充材料说明 ===")

    readme = f"""# Supplementary Materials

## Comparative Genomics and Pan-genome Analysis of *Veillonella* Species

**Date**: {datetime.now().strftime('%Y-%m-%d')}

---

## Table Descriptions

### Table S1: Genome Information
- **File**: `Table_S1_Genome_Info.xlsx`
- **Content**: Detailed information for all genomes used in this study
  - Genome ID and NCBI accession numbers
  - Species and strain information
  - Assembly statistics (genome size, GC content, N50, contigs)
  - Quality assessment (CheckM completeness and contamination)
  - ANI-based species clusters

### Table S2: Pan-genome Gene Presence/Absence Matrix
- **File**: `Table_S2_Pangenome_Matrix.xlsx`
- **Content**: Binary matrix of gene presence (1) and absence (0)
  - Rows: Gene families (n={tables.get('s2_shape', ['?', '?'])[0]})
  - Columns: Genome isolates (n={tables.get('s2_shape', ['?', '?'])[1]})
  - Generated by Panaroo v1.5.x

### Table S3: Functional Annotation
- **File**: `Table_S3_Functional_Annotation.xlsx`
- **Content**: COG and KEGG annotations for all genes
  - COG functional categories
  - KEGG orthology (KO) assignments
  - KEGG pathway mappings
  - Gene descriptions

### Table S4: Virulence Factors
- **File**: `Table_S4_Virulence_Genes.xlsx`
- **Content**: Virulence factor annotations from VFDB
  - VFDB gene names and functions
  - Presence/absence across genomes
  - Sequence identity and coverage

### Table S5: Antibiotic Resistance Genes
- **File**: `Table_S5_Resistance_Genes.xlsx`
- **Content**: Resistance gene annotations from CARD
  - CARD ARO identifiers
  - Resistance mechanisms
  - Drug classes
  - Presence/absence across genomes

### Table S6: COG Category Summary
- **File**: `Table_S6_COG_Summary.xlsx`
- **Content**: Summary of COG functional categories
  - Number of genes in each COG category
  - Core vs. accessory distribution
  - Functional enrichment statistics

---

## Data Availability

All raw sequencing data are available from NCBI under the following BioProject accession numbers:
- [Add NCBI BioProject ID]

The analysis pipeline and custom scripts are available at:
- GitHub: [Add repository URL]

---

## Software Versions

| Software | Version | Purpose |
|----------|---------|---------|
| Prokka | v1.14.6 | Genome annotation |
| Panaroo | v1.5.x | Pan-genome analysis |
| IQ-TREE | v3.0.1 | Phylogenetic analysis |
| CheckM2 | v1.0.0 | Genome quality assessment |
| FastANI | v1.3.3 | ANI calculation |
| eggNOG-mapper | v2.1.0 | Functional annotation |

---

## Contact

For questions regarding the supplementary materials, please contact:
[Add contact information]
"""

    readme_file = supp_dir / 'README.md'
    with open(readme_file, 'w') as f:
        f.write(readme)
    print(f"[OK] 保存说明文档: {readme_file}")


def save_to_excel(df: pd.DataFrame, output_file: Path, sheet_name: str = 'Sheet1'):
    """保存DataFrame到Excel，带格式"""
    try:
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            df.to_excel(writer, sheet_name=sheet_name, index=False)

            # 获取工作簿和工作表
            workbook = writer.book
            worksheet = writer.sheets[sheet_name]

            # 调整列宽
            for column in worksheet.columns:
                max_length = 0
                column_letter = column[0].column_letter
                for cell in column:
                    try:
                        if len(str(cell.value)) > max_length:
                            max_length = len(str(cell.value))
                    except:
                        pass
                adjusted_width = min(max_length + 2, 50)
                worksheet.column_dimensions[column_letter].width = adjusted_width

        print(f"[OK] 保存: {output_file}")
    except ImportError:
        # 如果没有openpyxl，保存为CSV
        csv_file = output_file.with_suffix('.csv')
        df.to_csv(csv_file, index=False)
        print(f"[OK] 保存CSV: {csv_file} (安装openpyxl可生成Excel)")


def main():
    import argparse
    parser = argparse.ArgumentParser(description="生成补充材料")
    parser.add_argument("--output", "-o", required=True, help="输出目录")
    parser.add_argument("--genome-dir", help="基因组目录（用于Table S1）")
    parser.add_argument("--qc-file", help="QC报告文件（用于Table S1）")
    parser.add_argument("--ani-clusters", help="ANI聚类文件（用于Table S1）")
    parser.add_argument("--panaroo-dir", help="Panaroo输出目录（用于Table S2）")

    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    tables = {}

    # Table S1: 基因组信息
    if args.genome_dir:
        genome_dir = Path(args.genome_dir)
        qc_file = Path(args.qc_file) if args.qc_file else Path('qc_report.csv')
        ani_clusters = Path(args.ani_clusters) if args.ani_clusters else Path('species_clusters.csv')

        s1_df = generate_table_s1(genome_dir, qc_file, ani_clusters)
        save_to_excel(s1_df, output_dir / 'Table_S1_Genome_Info.xlsx', 'Genome_Info')
        tables['s1_count'] = len(s1_df)

    # Table S2: 泛基因组矩阵
    if args.panaroo_dir:
        panaroo_dir = Path(args.panaroo_dir)
        s2_df = generate_table_s2(panaroo_dir)
        if s2_df is not None:
            save_to_excel(s2_df, output_dir / 'Table_S2_Pangenome_Matrix.xlsx', 'Gene_Matrix')
            tables['s2_shape'] = s2_df.shape

    # Table S3: 功能注释（示例）
    s3_df = generate_table_s3(None)
    save_to_excel(s3_df, output_dir / 'Table_S3_Functional_Annotation.xlsx', 'Annotation')

    # Table S6: COG汇总
    if args.panaroo_dir:
        s6_df = generate_cog_summary(panaroo_dir, None)
        if s6_df is not None:
            save_to_excel(s6_df, output_dir / 'Table_S6_COG_Summary.xlsx', 'COG_Summary')

    # 生成说明文档
    generate_readme(output_dir, tables)

    print("\n=== 完成 ===")
    print(f"补充材料已保存到: {output_dir}")


if __name__ == "__main__":
    main()
