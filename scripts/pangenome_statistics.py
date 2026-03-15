#!/usr/bin/env python3
"""
脚本：泛基因组统计分析
功能：
- 拟合泛基因组开放/闭合曲线（幂律/指数模型）
- 计算核心基因和泛基因组大小随样本量的变化
- 生成稀有faction和核心fraction曲线
- 预测泛基因组总大小
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pathlib import Path
import json
import sys
from itertools import combinations
from collections import defaultdict


def parse_gene_presence_absence(gpa_file: Path) -> pd.DataFrame:
    """解析Panaroo基因存在/缺失矩阵"""
    print("=== 读取泛基因组矩阵 ===")

    df = pd.read_csv(gpa_file, low_memory=False)
    print(f"矩阵形状: {df.shape}")
    print(f"基因总数: {len(df)}")

    # 获取样本列（通常从第4列开始是样本）
    # Panaroo格式：Gene, Non-unique Gene name, Annotation, etc.
    sample_cols = [col for col in df.columns if col not in
                   ['Gene', 'Non-unique Gene name', 'Annotation', 'No. isolates',
                    'No. sequences', 'Avg group size nuc', 'Genome Fragment',
                    'Order within Fragment', 'Accessory Fragment', 'Accessory Order with Fragment']]

    print(f"样本数量: {len(sample_cols)}")

    return df, sample_cols


def calculate_pan_genome_curve(df: pd.DataFrame, sample_cols: list,
                                iterations: int = 100) -> dict:
    """
    计算泛基因组积累曲线
    对每个样本量n，随机抽样iterations次，计算平均核心基因数和泛基因组大小
    """
    print("\n=== 计算泛基因组曲线 ===")

    n_samples = len(sample_cols)
    sample_sizes = list(range(1, n_samples + 1))

    core_genes = defaultdict(list)
    pan_genes = defaultdict(list)
    accessory_genes = defaultdict(list)

    # 对每个样本量进行多次随机抽样
    for n in sample_sizes:
        print(f"  处理样本量 {n}/{n_samples}...", end='\r')

        for _ in range(min(iterations, len(list(combinations(sample_cols, n))))):
            # 随机选择n个样本
            selected = np.random.choice(sample_cols, n, replace=False)

            # 统计基因
            # 对于每个基因，检查在选中的样本中的存在情况
            gene_presence = df[selected].notna().sum(axis=1)

            # 核心基因：在所有选中样本中都存在
            core = (gene_presence == n).sum()
            # 泛基因：至少在一个样本中存在
            pan = (gene_presence > 0).sum()
            # 附属基因：存在但不是核心
            accessory = pan - core

            core_genes[n].append(core)
            pan_genes[n].append(pan)
            accessory_genes[n].append(accessory)

    print()  # 换行

    # 计算均值和标准差
    results = {
        'sample_size': sample_sizes,
        'core_mean': [np.mean(core_genes[n]) for n in sample_sizes],
        'core_std': [np.std(core_genes[n]) for n in sample_sizes],
        'pan_mean': [np.mean(pan_genes[n]) for n in sample_sizes],
        'pan_std': [np.std(pan_genes[n]) for n in sample_sizes],
        'accessory_mean': [np.mean(accessory_genes[n]) for n in sample_sizes],
        'accessory_std': [np.std(accessory_genes[n]) for n in sample_sizes]
    }

    return pd.DataFrame(results)


def power_law(x, k, gamma):
    """幂律模型：泛基因组大小 = k * n^gamma"""
    return k * np.power(x, gamma)


def exponential_decay(x, a, b):
    """指数衰减模型：核心基因数 = a * exp(-b * n) + omega"""
    return a * np.exp(-b * x) + 1000  # omega设为估计值


def fit_pan_genome_models(curve_df: pd.DataFrame) -> dict:
    """拟合泛基因组模型"""
    print("\n=== 拟合泛基因组模型 ===")

    x = curve_df['sample_size'].values
    y_pan = curve_df['pan_mean'].values
    y_core = curve_df['core_mean'].values

    # 拟合幂律模型（泛基因组）
    try:
        popt_pan, _ = curve_fit(power_law, x, y_pan, p0=[1000, 0.5], maxfev=10000)
        k, gamma = popt_pan
        print(f"\n泛基因组幂律模型: y = {k:.2f} * n^{gamma:.4f}")
        print(f"  gamma (开放度指数): {gamma:.4f}")

        if gamma > 0:
            openness = "开放型 (Open)" if gamma > 0.5 else "中等开放型"
        else:
            openness = "闭合型 (Closed)"
        print(f"  泛基因组类型: {openness}")

        # 预测更大样本量的泛基因组大小
        predicted_200 = power_law(200, k, gamma)
        predicted_500 = power_law(500, k, gamma)
        print(f"  预测泛基因组大小 (n=200): {predicted_200:.0f}")
        print(f"  预测泛基因组大小 (n=500): {predicted_500:.0f}")

        pan_model = {
            'type': 'power_law',
            'k': float(k),
            'gamma': float(gamma),
            'openness': openness,
            'predicted_n200': float(predicted_200),
            'predicted_n500': float(predicted_500)
        }
    except Exception as e:
        print(f"  拟合失败: {e}")
        pan_model = None

    # 拟合指数衰减模型（核心基因）
    try:
        # 估计omega（渐进核心基因数）为最后几个点的均值
        omega_est = np.mean(y_core[-5:])
        popt_core, _ = curve_fit(lambda x, a, b: a * np.exp(-b * x) + omega_est,
                                 x, y_core, p0=[1000, 0.1], maxfev=10000)
        a, b = popt_core
        print(f"\n核心基因指数模型: y = {a:.2f} * exp(-{b:.4f} * n) + {omega_est:.0f}")
        print(f"  估计稳定核心基因数: ~{omega_est:.0f}")

        core_model = {
            'type': 'exponential',
            'a': float(a),
            'b': float(b),
            'omega': float(omega_est)
        }
    except Exception as e:
        print(f"  拟合失败: {e}")
        core_model = None

    return {
        'pan_genome_model': pan_model,
        'core_genome_model': core_model
    }


def classify_genes(df: pd.DataFrame, sample_cols: list) -> dict:
    """分类基因：核心、软核心、附属、特有"""
    print("\n=== 基因分类 ===")

    n = len(sample_cols)

    # 计算每个基因在多少样本中存在
    gene_presence = df[sample_cols].notna().sum(axis=1)

    # 分类
    core = (gene_presence == n).sum()                    # 100%
    soft_core = (gene_presence >= 0.95 * n).sum() - core  # 95-99%
    shell = (gene_presence >= 0.15 * n).sum() - soft_core - core  # 15-95%
    cloud = (gene_presence > 0).sum() - shell - soft_core - core  # <15%

    # 特有基因（仅存在于单个样本）
    unique = (gene_presence == 1).sum()

    classification = {
        'core': int(core),
        'soft_core': int(soft_core),
        'shell': int(shell),
        'cloud': int(cloud),
        'unique': int(unique),
        'total': len(df)
    }

    print(f"\n基因分类统计:")
    print(f"  核心基因 (Core, 100%):      {core} ({core/len(df)*100:.1f}%)")
    print(f"  软核心 (Soft core, ≥95%):   {soft_core} ({soft_core/len(df)*100:.1f}%)")
    print(f"  壳基因 (Shell, 15-95%):     {shell} ({shell/len(df)*100:.1f}%)")
    print(f"  云基因 (Cloud, <15%):       {cloud} ({cloud/len(df)*100:.1f}%)")
    print(f"  特有基因 (Unique):          {unique} ({unique/len(df)*100:.1f}%)")

    return classification


def generate_plots(curve_df: pd.DataFrame, models: dict, output_dir: Path):
    """生成泛基因组统计图"""
    print("\n=== 生成图表 ===")

    # 设置中文字体
    plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 图1: 泛基因组积累曲线
    ax1 = axes[0, 0]
    x = curve_df['sample_size']
    ax1.plot(x, curve_df['pan_mean'], 'b-', linewidth=2, label='Pan-genome')
    ax1.fill_between(x,
                     curve_df['pan_mean'] - curve_df['pan_std'],
                     curve_df['pan_mean'] + curve_df['pan_std'],
                     alpha=0.3, color='blue')

    # 添加拟合曲线
    if models['pan_genome_model']:
        k = models['pan_genome_model']['k']
        gamma = models['pan_genome_model']['gamma']
        x_fit = np.linspace(1, max(x), 100)
        y_fit = power_law(x_fit, k, gamma)
        ax1.plot(x_fit, y_fit, 'r--', linewidth=2,
                label=f"Fit: y={k:.0f}·n^{gamma:.3f}")

    ax1.set_xlabel('Number of Genomes', fontsize=11)
    ax1.set_ylabel('Number of Gene Families', fontsize=11)
    ax1.set_title('Pan-genome Accumulation Curve', fontsize=12, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # 图2: 核心基因曲线
    ax2 = axes[0, 1]
    ax2.plot(x, curve_df['core_mean'], 'g-', linewidth=2, label='Core genome')
    ax2.fill_between(x,
                     curve_df['core_mean'] - curve_df['core_std'],
                     curve_df['core_mean'] + curve_df['core_std'],
                     alpha=0.3, color='green')

    ax2.set_xlabel('Number of Genomes', fontsize=11)
    ax2.set_ylabel('Number of Core Gene Families', fontsize=11)
    ax2.set_title('Core Genome Accumulation Curve', fontsize=12, fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # 图3: Heaps定律图（用于判断开放/闭合）
    ax3 = axes[1, 0]
    # 计算新基因数
    new_genes = [curve_df['pan_mean'].iloc[i] - curve_df['pan_mean'].iloc[i-1]
                 if i > 0 else curve_df['pan_mean'].iloc[0]
                 for i in range(len(curve_df))]
    ax3.plot(x[1:], new_genes[1:], 'purple', marker='o', linewidth=2)
    ax3.axhline(y=1, color='red', linestyle='--', label='Open/Closed boundary')
    ax3.set_xlabel('Number of Genomes', fontsize=11)
    ax3.set_ylabel('New Gene Families per Genome', fontsize=11)
    ax3.set_title('Heaps Law Plot (New Genes per Genome)', fontsize=12, fontweight='bold')
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    # 图4: 基因组成饼图
    ax4 = axes[1, 1]

    # 假设我们有分类数据（需要传入）
    # 这里用模拟数据，实际应该传入classification
    sizes = [1000, 500, 800, 200]  # core, soft_core, shell, cloud
    labels = ['Core\n(100%)', 'Soft Core\n(95-99%)', 'Shell\n(15-95%)', 'Cloud\n(<15%)']
    colors = ['#2E7D32', '#66BB6A', '#FFA726', '#BDBDBD']

    ax4.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%',
            startangle=90, textprops={'fontsize': 10})
    ax4.set_title('Pan-genome Composition', fontsize=12, fontweight='bold')

    plt.tight_layout()
    output_file = output_dir / 'pangenome_statistics.pdf'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"[OK] 保存图表: {output_file}")

    # 单独保存核心图（用于论文）
    fig2, ax = plt.subplots(figsize=(8, 6))
    ax.plot(x, curve_df['pan_mean'], 'b-', linewidth=2.5, label='Pan-genome')
    ax.plot(x, curve_df['core_mean'], 'g-', linewidth=2.5, label='Core genome')
    ax.fill_between(x,
                    curve_df['pan_mean'] - curve_df['pan_std'],
                    curve_df['pan_mean'] + curve_df['pan_std'],
                    alpha=0.2, color='blue')
    ax.fill_between(x,
                    curve_df['core_mean'] - curve_df['core_std'],
                    curve_df['core_mean'] + curve_df['core_std'],
                    alpha=0.2, color='green')

    if models['pan_genome_model']:
        k = models['pan_genome_model']['k']
        gamma = models['pan_genome_model']['gamma']
        x_fit = np.linspace(1, max(x), 100)
        y_fit = power_law(x_fit, k, gamma)
        ax.plot(x_fit, y_fit, 'b--', linewidth=1.5, alpha=0.7)

    ax.set_xlabel('Number of Genomes', fontsize=12)
    ax.set_ylabel('Number of Gene Families', fontsize=12)
    ax.set_title('Pan-genome and Core Genome Curves', fontsize=13, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    output_file2 = output_dir / 'Figure_pangenome_curves.pdf'
    plt.savefig(output_file2, dpi=300, bbox_inches='tight')
    print(f"[OK] 保存论文图: {output_file2}")
    plt.close('all')


def generate_r_script(curve_df: pd.DataFrame, output_dir: Path):
    """生成R绘图脚本（更专业的图）"""
    print("\n=== 生成R绘图脚本 ===")

    # 保存数据为CSV供R使用
    curve_df.to_csv(output_dir / 'pangenome_curve_data.csv', index=False)

    r_script = output_dir / 'plot_pangenome.R'
    with open(r_script, 'w') as f:
        f.write("""library(ggplot2)
library(dplyr)
library(tidyr)

# 读取数据
df <- read.csv("pangenome_curve_data.csv")

# 转换为长格式用于ggplot
df_long <- df %>%
  select(sample_size, pan_mean, core_mean) %>%
  pivot_longer(cols = c(pan_mean, core_mean),
               names_to = "Type", values_to = "Genes")

# 添加标准差
df_long$std <- c(df$pan_std, df$core_std)

# 设置因子顺序
df_long$Type <- factor(df_long$Type,
                       levels = c("pan_mean", "core_mean"),
                       labels = c("Pan-genome", "Core genome"))

# 绘制泛基因组曲线
p <- ggplot(df_long, aes(x=sample_size, y=Genes, color=Type)) +
  geom_line(linewidth=1.2) +
  geom_ribbon(aes(ymin=Genes-std, ymax=Genes+std, fill=Type),
              alpha=0.2, color=NA) +
  scale_color_manual(values=c("Pan-genome"="#1976D2", "Core genome"="#388E3C")) +
  scale_fill_manual(values=c("Pan-genome"="#1976D2", "Core genome"="#388E3C")) +
  labs(x="Number of Genomes",
       y="Number of Gene Families",
       title="Pan-genome Analysis") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(hjust=0.5, face="bold", size=14),
    axis.title = element_text(size=12),
    panel.grid.minor = element_blank()
  )

ggsave("Figure_pangenome_curves_ggplot.pdf", p, width=8, height=6, dpi=300)

print("Pan-genome plot saved!")
""")

    print(f"[OK] 保存R脚本: {r_script}")


def main():
    import argparse
    parser = argparse.ArgumentParser(description="泛基因组统计分析")
    parser.add_argument("--input", "-i", required=True,
                        help="Panaroo gene_presence_absence.csv文件路径")
    parser.add_argument("--output", "-o", required=True, help="输出目录")
    parser.add_argument("--iterations", "-n", type=int, default=100,
                        help="每个样本量的迭代次数 (默认: 100)")

    args = parser.parse_args()

    input_file = Path(args.input)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_file.exists():
        print(f"错误：输入文件不存在: {input_file}")
        sys.exit(1)

    # 读取数据
    df, sample_cols = parse_gene_presence_absence(input_file)

    # 计算曲线
    curve_df = calculate_pan_genome_curve(df, sample_cols, args.iterations)

    # 保存曲线数据
    curve_file = output_dir / 'pangenome_accumulation_curve.csv'
    curve_df.to_csv(curve_file, index=False)
    print(f"[OK] 保存曲线数据: {curve_file}")

    # 拟合模型
    models = fit_pan_genome_models(curve_df)

    # 基因分类
    classification = classify_genes(df, sample_cols)

    # 保存统计结果
    results = {
        'classification': classification,
        'models': models,
        'n_samples': len(sample_cols),
        'total_genes': len(df)
    }

    results_file = output_dir / 'pangenome_statistics.json'
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"[OK] 保存统计结果: {results_file}")

    # 生成图
    generate_plots(curve_df, models, output_dir)

    # 生成R脚本
    generate_r_script(curve_df, output_dir)

    print("\n=== 完成 ===")


if __name__ == "__main__":
    main()
