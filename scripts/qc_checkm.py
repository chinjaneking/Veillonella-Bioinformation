#!/usr/bin/env python3
"""
脚本：基因组质量评估 (CheckM/CheckM2)
功能：
- 评估所有基因组的完整性和污染
- 过滤低质量基因组
- 生成质量报告表
"""

import pandas as pd
import subprocess
import json
from pathlib import Path
import sys


def run_checkm(input_dir: Path, output_dir: Path, threads: int = 8):
    """运行CheckM v1 质量评估"""
    print("=== 运行 CheckM v1 质量评估 ===")

    output_dir.mkdir(parents=True, exist_ok=True)

    # 准备输入文件列表
    genomes = list(input_dir.glob("*.fna"))
    print(f"找到 {len(genomes)} 个基因组")

    # 首先创建checkm的输入文件列表（bin文件夹格式）
    checkm_input = output_dir / "input_bins"
    checkm_input.mkdir(exist_ok=True)

    # 创建软链接而不是复制（节省空间）
    for genome in genomes:
        link_path = checkm_input / f"{genome.stem}.fna"
        if not link_path.exists():
            link_path.symlink_to(genome.absolute())

    print(f"创建了 {len(genomes)} 个软链接到输入目录")

    # 运行CheckM lineage_wf
    cmd = [
        "checkm", "lineage_wf",
        "-x", "fna",
        "-t", str(threads),
        "--tab_table",
        "-f", str(output_dir / "checkm_output.tsv"),
        str(checkm_input),
        str(output_dir / "checkm_results")
    ]

    print(f"运行命令: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"警告：CheckM 返回非零状态: {result.returncode}")
            print(f"stderr: {result.stderr}")
        print("✓ CheckM 运行完成")
    except subprocess.CalledProcessError as e:
        print(f"错误：CheckM 运行失败: {e}")
        print(f"stderr: {e.stderr}")
        sys.exit(1)
    except FileNotFoundError:
        print("错误：未找到 checkm 命令。请确保 CheckM 已安装。")
        sys.exit(1)


def filter_genomes(checkm_output: Path, completeness_threshold: float,
                   contamination_threshold: float) -> pd.DataFrame:
    """根据质量阈值过滤基因组"""
    print("\n=== 过滤基因组 ===")

    quality_file = checkm_output / "checkm_output.tsv"
    if not quality_file.exists():
        print(f"错误：未找到质量报告: {quality_file}")
        print("请确保CheckM已成功运行")
        sys.exit(1)

    # 读取质量报告（CheckM v1输出格式）
    # 列：Bin Id, Marker lineage, # genomes, # markers, # marker sets, 0, 1, 2, 3, 4, 5+, Completeness, Contamination, Strain heterogeneity
    df = pd.read_csv(quality_file, sep='\t')

    # 重命名列以便处理
    df = df.rename(columns={
        'Bin Id': 'Name',
        'Completeness': 'Completeness',
        'Contamination': 'Contamination'
    })

    # 确保数值类型正确
    df['Completeness'] = pd.to_numeric(df['Completeness'], errors='coerce')
    df['Contamination'] = pd.to_numeric(df['Contamination'], errors='coerce')

    # 应用过滤
    df['Pass_QC'] = (df['Completeness'] >= completeness_threshold) & \
                    (df['Contamination'] <= contamination_threshold)

    # 统计
    passed = df['Pass_QC'].sum()
    failed = len(df) - passed

    print(f"通过 QC: {passed}/{len(df)} ({passed/len(df)*100:.1f}%)")
    print(f"未通过: {failed}/{len(df)} ({failed/len(df)*100:.1f}%)")

    if len(df) > 0:
        print(f"\n质量统计:")
        print(f"  完整性均值: {df['Completeness'].mean():.2f}%")
        print(f"  完整性范围: {df['Completeness'].min():.2f}% - {df['Completeness'].max():.2f}%")
        print(f"  污染均值: {df['Contamination'].mean():.2f}%")
        print(f"  污染范围: {df['Contamination'].min():.2f}% - {df['Contamination'].max():.2f}%")

    return df


def generate_qc_report(qc_df: pd.DataFrame, output_file: Path):
    """生成QC报告"""
    print(f"\n=== 生成 QC 报告 ===")

    # 保存完整报告
    qc_df.to_csv(output_file, index=False)
    print(f"✓ 保存报告: {output_file}")

    # 保存通过QC的基因组列表
    passed_list = output_file.parent / "passed_qc_genomes.txt"
    passed_genomes = qc_df[qc_df['Pass_QC']]['Name'].tolist()
    with open(passed_list, 'w') as f:
        for genome in passed_genomes:
            f.write(f"{genome}\n")
    print(f"✓ 保存通过列表: {passed_list} ({len(passed_genomes)} 个)")

    # 生成质量分布图（如R可用）
    r_script = output_file.parent / "plot_qc_stats.R"
    with open(r_script, 'w') as f:
        f.write("""
library(ggplot2)
library(dplyr)

# 读取数据
df <- read.csv("qc_report.csv")

# 质量分布图
p1 <- ggplot(df, aes(x=Completeness, y=Contamination, color=Pass_QC)) +
  geom_point(alpha=0.6, size=3) +
  geom_vline(xintercept=90, linetype="dashed", color="red") +
  geom_hline(yintercept=5, linetype="dashed", color="red") +
  scale_color_manual(values=c("TRUE"="#2E7D32", "FALSE"="#C62828"),
                     labels=c("TRUE"="Pass", "FALSE"="Fail")) +
  labs(x="Completeness (%)", y="Contamination (%)",
       title="Genome Quality Assessment",
       color="QC Status") +
  theme_minimal() +
  theme(legend.position="bottom")

ggsave("qc_scatter.pdf", p1, width=8, height=6)

# 直方图
p2 <- ggplot(df, aes(x=Completeness, fill=Pass_QC)) +
  geom_histogram(binwidth=2, alpha=0.7) +
  scale_fill_manual(values=c("TRUE"="#2E7D32", "FALSE"="#C62828")) +
  labs(x="Completeness (%)", y="Count",
       title="Completeness Distribution") +
  theme_minimal()

ggsave("qc_completeness_hist.pdf", p2, width=6, height=4)

print("QC plots saved!")
""")
    print(f"✓ 生成R绘图脚本: {r_script}")


def main():
    import argparse
    parser = argparse.ArgumentParser(description="基因组质量评估")
    parser.add_argument("--input", "-i", required=True, help="输入基因组目录")
    parser.add_argument("--output", "-o", required=True, help="输出目录")
    parser.add_argument("--completeness", "-c", type=float, default=90,
                        help="最小完整性阈值 (默认: 90)")
    parser.add_argument("--contamination", "-x", type=float, default=5,
                        help="最大污染阈值 (默认: 5)")
    parser.add_argument("--threads", "-t", type=int, default=8,
                        help="线程数 (默认: 8)")
    parser.add_argument("--skip-checkm", action="store_true",
                        help="跳过CheckM运行，只进行过滤")

    args = parser.parse_args()

    input_dir = Path(args.input)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not args.skip_checkm:
        run_checkm(input_dir, output_dir / "checkm_output", args.threads)

    # 过滤基因组
    qc_df = filter_genomes(
        output_dir / "checkm_output",
        args.completeness,
        args.contamination
    )

    # 生成报告
    generate_qc_report(qc_df, output_dir / "qc_report.csv")

    print("\n=== 完成 ===")


if __name__ == "__main__":
    main()
