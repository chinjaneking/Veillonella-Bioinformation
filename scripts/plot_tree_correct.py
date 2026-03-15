#!/usr/bin/env python3
"""
正确绘制系统发育树 - 使用Bio.Phylo的draw_graphviz方式
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import Phylo
import pandas as pd
from pathlib import Path
import sys
import numpy as np


def load_species_groups(clusters_file):
    """加载物种分组"""
    df = pd.read_csv(clusters_file)
    return dict(zip(df['Genome'], df['Species_Group']))


def get_color_map(species_groups):
    """为物种组分配颜色"""
    unique_species = sorted(set(species_groups.values()))
    # 使用更好的颜色方案
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
              '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    return dict(zip(unique_species, colors[:len(unique_species)]))


def plot_tree_correct(tree_file, clusters_file, output_prefix):
    """正确绘制系统发育树"""
    print("=== 绘制系统发育树（修正版）===")

    # 加载树
    tree = Phylo.read(tree_file, "newick")
    print(f"树包含 {len(tree.get_terminals())} 个末端节点")

    # 加载物种分组
    species_groups = load_species_groups(clusters_file)
    color_map = get_color_map(species_groups)

    # 创建更大的图形
    fig, ax = plt.subplots(figsize=(14, 10))

    # 使用Phylo.draw绘制，但定制样式
    # 先获取树的布局信息
    tree.ladderize()  # 梯形化树

    # 手动绘制树
    def get_coordinates(tree):
        """计算树的坐标"""
        # 获取所有终端节点
        terminals = tree.get_terminals()
        n_terminals = len(terminals)

        # 为每个终端节点分配y坐标
        y_coords = {}
        for i, clade in enumerate(terminals):
            y_coords[clade] = n_terminals - i - 1

        # 计算内部节点的y坐标（子节点的中点）
        for clade in tree.find_clades(order='postorder'):
            if not clade.is_terminal():
                children = clade.clades
                y_coords[clade] = sum([y_coords[c] for c in children]) / len(children)

        # 计算x坐标（分支长度累加）
        x_coords = {}
        x_coords[tree.root] = 0

        for clade in tree.find_clades(order='level'):
            if clade != tree.root:
                parent = next(tree.find_clades(lambda c: clade in c.clades))
                branch_length = clade.branch_length if clade.branch_length else 0
                x_coords[clade] = x_coords[parent] + branch_length

        return x_coords, y_coords, terminals

    x_coords, y_coords, terminals = get_coordinates(tree)

    # 绘制分支
    for clade in tree.find_clades():
        if clade != tree.root:
            parent = next(tree.find_clades(lambda c: clade in c.clades))
            x_parent = x_coords[parent]
            x_clade = x_coords[clade]
            y_parent = y_coords[parent]
            y_clade = y_coords[clade]

            # 水平线
            ax.plot([x_parent, x_clade], [y_clade, y_clade],
                   'k-', linewidth=1.5, alpha=0.7)
            # 垂直线
            if not clade.is_terminal():
                children = clade.clades
                y_min = min([y_coords[c] for c in children])
                y_max = max([y_coords[c] for c in children])
                ax.plot([x_clade, x_clade], [y_min, y_max],
                       'k-', linewidth=1.5, alpha=0.7)

    # 添加末端标签
    for i, clade in enumerate(terminals):
        name = clade.name
        x = x_coords[clade]
        y = y_coords[clade]

        # 获取颜色
        species = species_groups.get(name, 'Unknown')
        color = color_map.get(species, 'black')

        # 添加标签
        ax.text(x + 0.001, y, name, fontsize=9,
               color=color, fontweight='bold', va='center')

        # 添加点标记
        ax.plot(x, y, 'o', color=color, markersize=6)

    # 设置坐标轴
    ax.set_xlim(-0.005, max(x_coords.values()) + 0.05)
    ax.set_ylim(-1, len(terminals))

    # 移除坐标轴
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    # 添加比例尺
    scale_x = max(x_coords.values()) * 0.05
    scale_y = -0.5
    scale_length = 0.05  # 5% divergence
    ax.plot([scale_x, scale_x + scale_length], [scale_y, scale_y],
           'k-', linewidth=2)
    ax.text(scale_x + scale_length/2, scale_y - 0.5,
           f'{scale_length:.3f} substitutions/site',
           ha='center', fontsize=9)

    # 添加图例
    legend_patches = [mpatches.Patch(color=color, label=species.replace('_', ' '))
                      for species, color in sorted(color_map.items())]
    ax.legend(handles=legend_patches, loc='lower right',
             title='Species', frameon=True, fontsize=10)

    # 添加标题
    ax.set_title('Phylogenetic Tree of Veillonella Species\n(Based on 1,125 Core Genes)',
                fontsize=14, fontweight='bold', pad=20)

    # 保存
    output_pdf = f"{output_prefix}.pdf"
    output_png = f"{output_prefix}.png"

    plt.tight_layout()
    plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
    plt.savefig(output_png, dpi=300, bbox_inches='tight')

    print(f"[OK] 保存PDF: {output_pdf}")
    print(f"[OK] 保存PNG: {output_png}")

    plt.close()
    print("\n=== 完成 ===")


def main():
    import argparse
    parser = argparse.ArgumentParser(description="正确绘制系统发育树")
    parser.add_argument("--tree", "-t", required=True, help="树文件(Newick)")
    parser.add_argument("--clusters", "-c", required=True, help="物种聚类CSV")
    parser.add_argument("--output", "-o", required=True, help="输出前缀")

    args = parser.parse_args()

    plot_tree_correct(args.tree, args.clusters, args.output)


if __name__ == "__main__":
    main()
