#!/usr/bin/env python3
"""
脚本：基础基因组质量评估
功能：
- 计算基础基因组统计（大小、GC含量、N50、contig数）
- 基于统计信息进行质量过滤
- 生成QC报告

注意：这是CheckM的简化版本，用于快速演示。正式发表建议使用CheckM。
"""

import pandas as pd
from pathlib import Path
import sys
from Bio import SeqIO
import numpy as np


def calculate_n50(contig_lengths):
    """计算N50"""
    sorted_lengths = sorted(contig_lengths, reverse=True)
    total = sum(sorted_lengths)
    cumulative = 0
    for length in sorted_lengths:
        cumulative += length
        if cumulative >= total / 2:
            return length
    return sorted_lengths[-1] if sorted_lengths else 0


def calculate_stats(genome_file: Path) -> dict:
    """计算单个基因组的基础统计"""
    contig_lengths = []
    gc_count = 0
    total_bases = 0

    for record in SeqIO.parse(genome_file, "fasta"):
        seq = str(record.seq)
        contig_lengths.append(len(seq))
        gc_count += seq.count('G') + seq.count('C') + seq.count('g') + seq.count('c')
        total_bases += len(seq)

    return {
        'Genome_Size': total_bases,
        'Contigs': len(contig_lengths),
        'N50': calculate_n50(contig_lengths),
        'GC_Content': (gc_count / total_bases * 100) if total_bases > 0 else 0,
        'Longest_Contig': max(contig_lengths) if contig_lengths else 0
    }


def run_qc(input_dir: Path, output_dir: Path):
    """运行QC分析"""
    print("=== 运行基础基因组QC ===")

    genomes = sorted(input_dir.glob("*.fna"))
    print(f"找到 {len(genomes)} 个基因组")

    results = []
    for genome in genomes:
        print(f"  分析: {genome.stem}", flush=True)
        stats = calculate_stats(genome)
        stats['Name'] = genome.stem
        results.append(stats)

    df = pd.DataFrame(results)

    # 质量判断（基于经验阈值）
    # 由于韦荣氏菌基因组约2Mb，我们可以设置合理阈值
    df['Completeness_Est'] = ((df['Genome_Size'] >= 1800000) & (df['Genome_Size'] <= 2800000)) * 100
    # 污染估算基于contig数量（简单估算）
    df['Contamination_Est'] = np.where(df['Contigs'] <= 100, 0,
                                       np.where(df['Contigs'] <= 500, 2, 5))

    # QC判断
    df['Pass_QC'] = (df['Genome_Size'] >= 1500000) & \
                    (df['Genome_Size'] <= 3500000) & \
                    (df['Contigs'] <= 500)

    return df


def generate_report(df: pd.DataFrame, output_dir: Path):
    """生成QC报告"""
    print("\n=== 生成QC报告 ===")

    output_dir.mkdir(parents=True, exist_ok=True)

    # 保存完整报告
    report_file = output_dir / "qc_report.csv"
    df.to_csv(report_file, index=False)
    print(f"[OK] 保存报告: {report_file}")

    # 保存通过QC的基因组列表
    passed_list = output_dir / "passed_qc_genomes.txt"
    passed_genomes = df[df['Pass_QC']]['Name'].tolist()
    with open(passed_list, 'w') as f:
        for genome in passed_genomes:
            f.write(f"{genome}\n")
    print(f"[OK] 通过QC的基因组: {len(passed_genomes)}/{len(df)}")

    # 统计信息
    print(f"\n质量统计:")
    print(f"  基因组大小范围: {df['Genome_Size'].min():,} - {df['Genome_Size'].max():,} bp")
    print(f"  平均基因组大小: {df['Genome_Size'].mean():,.0f} bp")
    print(f"  GC含量范围: {df['GC_Content'].min():.1f}% - {df['GC_Content'].max():.1f}%")
    print(f"  Contig数范围: {df['Contigs'].min()} - {df['Contigs'].max()}")
    print(f"  平均N50: {df['N50'].mean():,.0f} bp")

    # 生成R绘图脚本
    r_script = output_dir / "plot_qc.R"
    with open(r_script, 'w') as f:
        f.write("""
library(ggplot2)
library(dplyr)

df <- read.csv("qc_report.csv")

# 基因组大小分布
p1 <- ggplot(df, aes(x=Genome_Size/1e6, fill=Pass_QC)) +
  geom_histogram(binwidth=0.1, alpha=0.7) +
  scale_fill_manual(values=c("TRUE"="#2E7D32", "FALSE"="#C62828"),
                    labels=c("Pass", "Fail")) +
  labs(x="Genome Size (Mb)", y="Count",
       title="Genome Size Distribution") +
  theme_minimal() +
  theme(legend.position="bottom")

ggsave("qc_genome_size.pdf", p1, width=6, height=4)

# GC含量分布
p2 <- ggplot(df, aes(x=GC_Content, fill=Pass_QC)) +
  geom_histogram(binwidth=1, alpha=0.7) +
  scale_fill_manual(values=c("TRUE"="#2E7D32", "FALSE"="#C62828")) +
  labs(x="GC Content (%)", y="Count",
       title="GC Content Distribution") +
  theme_minimal() +
  theme(legend.position="bottom")

ggsave("qc_gc_content.pdf", p2, width=6, height=4)

# N50 vs Contigs
p3 <- ggplot(df, aes(x=Contigs, y=N50/1e6, color=Pass_QC)) +
  geom_point(size=3, alpha=0.6) +
  scale_color_manual(values=c("TRUE"="#2E7D32", "FALSE"="#C62828"),
                     labels=c("Pass", "Fail")) +
  labs(x="Number of Contigs", y="N50 (Mb)",
       title="Assembly Quality") +
  theme_minimal() +
  theme(legend.position="bottom")

ggsave("qc_assembly.pdf", p3, width=6, height=4)

print("QC plots saved!")
""")
    print(f"[OK] 生成R绘图脚本: {r_script}")


def main():
    import argparse
    parser = argparse.ArgumentParser(description="基础基因组QC")
    parser.add_argument("--input", "-i", required=True, help="输入基因组目录")
    parser.add_argument("--output", "-o", required=True, help="输出目录")

    args = parser.parse_args()

    input_dir = Path(args.input)
    output_dir = Path(args.output)

    if not input_dir.exists():
        print(f"错误：输入目录不存在: {input_dir}")
        sys.exit(1)

    # 运行QC
    df = run_qc(input_dir, output_dir)

    # 生成报告
    generate_report(df, output_dir)

    print("\n=== 完成 ===")


if __name__ == "__main__":
    main()
