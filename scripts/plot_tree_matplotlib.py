#!/usr/bin/env python3
"""
使用matplotlib绘制系统发育树
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import Phylo
import pandas as pd
from pathlib import Path
import sys


def load_species_groups(clusters_file):
    """加载物种分组"""
    df = pd.read_csv(clusters_file)
    return dict(zip(df['Genome'], df['Species_Group']))


def get_color_map(species_groups):
    """为物种组分配颜色"""
    unique_species = sorted(set(species_groups.values()))
    colors = plt.cm.Set2(range(len(unique_species)))
    return dict(zip(unique_species, colors))


def plot_tree(tree_file, clusters_file, output_prefix):
    """绘制系统发育树"""
    print("=== 绘制系统发育树 ===")

    # 加载树
    tree = Phylo.read(tree_file, "newick")
    print(f"树包含 {len(tree.get_terminals())} 个末端节点")

    # 加载物种分组
    species_groups = load_species_groups(clusters_file)
    color_map = get_color_map(species_groups)

    # 创建图形
    fig, ax = plt.subplots(figsize=(12, 8))

    # 绘制树
    # 使用直角布局
    Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False)

    # 获取所有叶子节点
    leaves = [cl.name for cl in tree.get_terminals()]

    # 为叶子标签着色
    for i, text in enumerate(ax.get_yticklabels()):
        label = text.get_text()
        if label in species_groups:
            species = species_groups[label]
            color = color_map[species]
            text.set_color(color)
            text.set_fontweight('bold')

    # 添加图例
    legend_patches = [mpatches.Patch(color=color, label=species)
                      for species, color in color_map.items()]
    ax.legend(handles=legend_patches, loc='upper left', title='Species')

    # 设置标题
    ax.set_title('Phylogenetic Tree of Veillonella Species', fontsize=14, fontweight='bold')

    # 保存
    output_pdf = f"{output_prefix}.pdf"
    output_png = f"{output_prefix}.png"

    plt.tight_layout()
    plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
    plt.savefig(output_png, dpi=300, bbox_inches='tight')

    print(f"[OK] 保存PDF: {output_pdf}")
    print(f"[OK] 保存PNG: {output_png}")

    plt.close()

    # 创建环形树（使用简单的圆形布局）
    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection='polar'))

    # 简化版：只显示叶子节点和标签
    n_leaves = len(leaves)
    angles = [i * 2 * 3.14159 / n_leaves for i in range(n_leaves)]

    for i, (leaf, angle) in enumerate(zip(leaves, angles)):
        species = species_groups.get(leaf, 'Unknown')
        color = color_map.get(species, 'gray')
        ax.plot([0, angle], [0, 1], color=color, alpha=0.6, linewidth=2)
        ax.text(angle, 1.1, leaf, rotation=angle*180/3.14159-90,
                fontsize=6, ha='center', va='center', color=color)

    ax.set_ylim(0, 1.5)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_title('Circular Phylogenetic Tree of Veillonella', fontsize=14, fontweight='bold', pad=20)

    output_circ = f"{output_prefix}_circular.pdf"
    plt.savefig(output_circ, dpi=300, bbox_inches='tight')
    print(f"[OK] 保存环形树: {output_circ}")

    print("\n=== 完成 ===")


def main():
    import argparse
    parser = argparse.ArgumentParser(description="绘制系统发育树")
    parser.add_argument("--tree", "-t", required=True, help="树文件(Newick)")
    parser.add_argument("--clusters", "-c", required=True, help="物种聚类CSV")
    parser.add_argument("--output", "-o", required=True, help="输出前缀")

    args = parser.parse_args()

    plot_tree(args.tree, args.clusters, args.output)


if __name__ == "__main__":
    main()
