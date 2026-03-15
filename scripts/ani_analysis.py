#!/usr/bin/env python3
"""
脚本：ANI分析 (fastANI)
功能：
- 计算所有基因组两两之间的ANI
- 生成ANI矩阵和热图
- 基于ANI阈值进行物种聚类
- 生成ANI分布统计
"""

import pandas as pd
import numpy as np
import subprocess
from pathlib import Path
import sys
from itertools import combinations
import json


def get_genome_list(genome_dir: Path) -> list:
    """获取基因组列表"""
    genomes = sorted(genome_dir.glob("*.fna"))
    if not genomes:
        print(f"错误：在 {genome_dir} 中未找到 .fna 文件")
        sys.exit(1)
    return genomes


def run_fastani_pairwise(genomes: list, output_dir: Path, threads: int = 8) -> Path:
    """运行fastANI两两比对"""
    print("=== 运行 fastANI 两两比对 ===")
    print(f"基因组数量: {len(genomes)}")
    print(f"比对组合数: {len(genomes) * (len(genomes) - 1) // 2}")

    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "ani_results.tsv"

    # 创建查询文件列表
    query_list = output_dir / "query_list.txt"
    with open(query_list, 'w') as f:
        for genome in genomes:
            f.write(f"{genome}\n")

    # 运行fastANI
    cmd = [
        "fastANI",
        "--ql", str(query_list),
        "--rl", str(query_list),
        "-o", str(output_file),
        "-t", str(threads)
    ]

    print(f"运行命令: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print("✓ fastANI 运行完成")
    except subprocess.CalledProcessError as e:
        print(f"错误：fastANI 运行失败: {e}")
        print(f"stderr: {e.stderr}")
        sys.exit(1)
    except FileNotFoundError:
        print("错误：未找到 fastANI 命令。请确保已安装 fastANI。")
        sys.exit(1)

    return output_file


def parse_ani_results(ani_file: Path, genomes: list) -> pd.DataFrame:
    """解析ANI结果并构建矩阵"""
    print("\n=== 解析ANI结果 ===")

    # 读取ANI结果
    columns = ['Query', 'Reference', 'ANI', 'Bidirectional_fragments',
               'Total_fragments_Query', 'Total_fragments_Reference']
    df = pd.read_csv(ani_file, sep='\t', header=None, names=columns)

    # 获取所有基因组名称（不含路径和扩展名）
    genome_names = sorted(set([g.stem for g in genomes]))
    n = len(genome_names)

    print(f"解析到 {len(df)} 条ANI结果")
    print(f"基因组数量: {n}")

    # 构建矩阵
    matrix = pd.DataFrame(index=genome_names, columns=genome_names, dtype=float)

    # 填充矩阵（对称矩阵）
    for _, row in df.iterrows():
        query = Path(row['Query']).stem
        ref = Path(row['Reference']).stem
        ani = row['ANI']

        matrix.loc[query, ref] = ani
        matrix.loc[ref, query] = ani

    # 对角线设为100
    for name in genome_names:
        matrix.loc[name, name] = 100.0

    # 转换为数值类型
    matrix = matrix.astype(float)

    return matrix


def classify_species_by_ani(matrix: pd.DataFrame, species_threshold: float = 95.0) -> dict:
    """基于ANI进行物种聚类"""
    print(f"\n=== 物种聚类 (ANI > {species_threshold}%) ===")

    genome_names = matrix.index.tolist()
    n = len(genome_names)

    # 构建邻接表（同种关系）
    same_species = {name: set() for name in genome_names}

    for i, name1 in enumerate(genome_names):
        for j, name2 in enumerate(genome_names):
            if i < j:  # 避免重复
                ani = matrix.loc[name1, name2]
                if ani >= species_threshold:
                    same_species[name1].add(name2)
                    same_species[name2].add(name1)

    # 使用并查集进行聚类
    parent = {name: name for name in genome_names}

    def find(x):
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]

    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    # 合并同种基因组
    for name1 in genome_names:
        for name2 in same_species[name1]:
            union(name1, name2)

    # 收集聚类结果
    clusters = {}
    for name in genome_names:
        root = find(name)
        if root not in clusters:
            clusters[root] = []
        clusters[root].append(name)

    # 按大小排序
    sorted_clusters = sorted(clusters.values(), key=len, reverse=True)

    print(f"识别到 {len(sorted_clusters)} 个物种组:")
    for i, cluster in enumerate(sorted_clusters, 1):
        print(f"  物种组 {i}: {len(cluster)} 个基因组")
        if len(cluster) <= 5:
            for name in cluster:
                print(f"    - {name}")
        else:
            for name in cluster[:3]:
                print(f"    - {name}")
            print(f"    ... 等共 {len(cluster)} 个")

    # 返回聚类标签
    labels = {}
    for i, cluster in enumerate(sorted_clusters, 1):
        for name in cluster:
            labels[name] = f"Species_{i}"

    return labels


def generate_ani_stats(matrix: pd.DataFrame, output_dir: Path):
    """生成ANI统计信息"""
    print("\n=== 生成ANI统计 ===")

    genome_names = matrix.index.tolist()
    n = len(genome_names)

    # 收集所有两两ANI值（不包括对角线）
    ani_values = []
    for i in range(n):
        for j in range(i+1, n):
            ani = matrix.iloc[i, j]
            if not np.isnan(ani):
                ani_values.append(ani)

    ani_series = pd.Series(ani_values)

    stats = {
        'total_comparisons': len(ani_values),
        'mean_ani': ani_series.mean(),
        'median_ani': ani_series.median(),
        'std_ani': ani_series.std(),
        'min_ani': ani_series.min(),
        'max_ani': ani_series.max(),
        'q25': ani_series.quantile(0.25),
        'q75': ani_series.quantile(0.75),
        'ani_90_95': len([x for x in ani_values if 90 <= x < 95]),
        'ani_95_100': len([x for x in ani_values if 95 <= x <= 100]),
        'ani_below_90': len([x for x in ani_values if x < 90])
    }

    print(f"\nANI统计:")
    print(f"  总比较数: {stats['total_comparisons']}")
    print(f"  均值: {stats['mean_ani']:.2f}%")
    print(f"  中位数: {stats['median_ani']:.2f}%")
    print(f"  标准差: {stats['std_ani']:.2f}%")
    print(f"  范围: {stats['min_ani']:.2f}% - {stats['max_ani']:.2f}%")
    print(f"\n分布:")
    print(f"  ANI < 90%: {stats['ani_below_90']} ({stats['ani_below_90']/len(ani_values)*100:.1f}%)")
    print(f"  90% ≤ ANI < 95%: {stats['ani_90_95']} ({stats['ani_90_95']/len(ani_values)*100:.1f}%)")
    print(f"  ANI ≥ 95%: {stats['ani_95_100']} ({stats['ani_95_100']/len(ani_values)*100:.1f}%)")

    # 保存统计
    stats_file = output_dir / "ani_statistics.json"
    with open(stats_file, 'w') as f:
        json.dump(stats, f, indent=2)
    print(f"\n✓ 保存统计: {stats_file}")

    return stats


def generate_r_heatmap_script(matrix: pd.DataFrame, output_dir: Path):
    """生成R热图绘制脚本"""
    print("\n=== 生成热图绘制脚本 ===")

    r_script = output_dir / "plot_ani_heatmap.R"

    with open(r_script, 'w') as f:
        f.write("""library(ggplot2)
library(reshape2)
library(dplyr)

# 读取ANI矩阵
ani_matrix <- read.csv("ani_matrix.csv", row.names=1, check.names=FALSE)

# 转换为长格式
ani_long <- melt(as.matrix(ani_matrix))
colnames(ani_long) <- c("Query", "Reference", "ANI")

# 热图
p1 <- ggplot(ani_long, aes(x=Query, y=Reference, fill=ANI)) +
  geom_tile() +
  scale_fill_gradient2(low="#2166AC", mid="#F7F7F7", high="#B2182B",
                       midpoint=92.5, limits=c(80, 100),
                       name="ANI (%)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=6),
    axis.text.y = element_text(size=6),
    panel.grid = element_blank()
  ) +
  labs(title="ANI Heatmap", x="", y="")

ggsave("ani_heatmap.pdf", p1, width=12, height=10)

# 带聚类的热图 (pheatmap)
library(pheatmap)

# 处理NA值
ani_mat <- as.matrix(ani_matrix)
ani_mat[is.na(ani_mat)] <- 80  # 用80填充NA

# 绘制热图
pdf("ani_heatmap_clustered.pdf", width=12, height=10)
pheatmap(ani_mat,
         clustering_method="average",
         color=colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
         breaks=seq(80, 100, length.out=101),
         main="ANI Heatmap (Hierarchical Clustering)",
         fontsize_row=6,
         fontsize_col=6)
dev.off()

# ANI分布直方图
p2 <- ggplot(ani_long %>% filter(Query != Reference), aes(x=ANI)) +
  geom_histogram(binwidth=1, fill="steelblue", alpha=0.7) +
  geom_vline(xintercept=95, color="red", linetype="dashed") +
  annotate("text", x=95, y=Inf, vjust=2, label="Species boundary (95%)",
           color="red", size=3) +
  labs(x="ANI (%)", y="Count", title="ANI Distribution") +
  theme_minimal()

ggsave("ani_distribution.pdf", p2, width=6, height=4)

print("ANI plots saved!")
""")

    print(f"✓ 保存R脚本: {r_script}")
    return r_script


def main():
    import argparse
    parser = argparse.ArgumentParser(description="ANI分析")
    parser.add_argument("--input", "-i", required=True, help="输入基因组目录")
    parser.add_argument("--output", "-o", required=True, help="输出目录")
    parser.add_argument("--threads", "-t", type=int, default=8, help="线程数")
    parser.add_argument("--species-threshold", "-s", type=float, default=95.0,
                        help="物种边界阈值 (默认: 95%)")
    parser.add_argument("--skip-fastani", action="store_true",
                        help="跳过fastANI运行（使用已有结果）")

    args = parser.parse_args()

    input_dir = Path(args.input)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 获取基因组列表
    genomes = get_genome_list(input_dir)

    # 运行fastANI
    if not args.skip_fastani:
        ani_file = run_fastani_pairwise(genomes, output_dir, args.threads)
    else:
        ani_file = output_dir / "ani_results.tsv"
        if not ani_file.exists():
            print(f"错误：未找到ANI结果文件: {ani_file}")
            sys.exit(1)

    # 解析结果
    matrix = parse_ani_results(ani_file, genomes)

    # 保存矩阵
    matrix_file = output_dir / "ani_matrix.csv"
    matrix.to_csv(matrix_file)
    print(f"✓ 保存ANI矩阵: {matrix_file}")

    # 物种聚类
    species_labels = classify_species_by_ani(matrix, args.species_threshold)
    labels_df = pd.DataFrame(list(species_labels.items()),
                             columns=['Genome', 'Species_Group'])
    labels_df.to_csv(output_dir / "species_clusters.csv", index=False)
    print(f"✓ 保存物种聚类: {output_dir / 'species_clusters.csv'}")

    # 生成统计
    generate_ani_stats(matrix, output_dir)

    # 生成R脚本
    generate_r_heatmap_script(matrix, output_dir)

    print("\n=== 完成 ===")
    print(f"输出目录: {output_dir}")


if __name__ == "__main__":
    main()
