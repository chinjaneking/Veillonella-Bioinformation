#!/usr/bin/env python3
"""
最终版系统发育树绘制 - 显示GCA编号和物种名
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import Phylo
import pandas as pd
from pathlib import Path


def load_species_info(clusters_file):
    """加载物种信息"""
    df = pd.read_csv(clusters_file)
    # 返回 {GCA: (full_species, group)}
    return dict(zip(df['Genome'], zip(df['Species'], df['Species_Group'])))


def get_color_map(species_groups_list):
    """为物种组分配颜色"""
    unique_species = sorted(set(species_groups_list))
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    return dict(zip(unique_species, colors[:len(unique_species)]))


def plot_tree_final(tree_file, clusters_file, output_prefix):
    """绘制最终的系统发育树"""
    print("=== 绘制系统发育树（显示物种名）===")

    # 加载树
    tree = Phylo.read(tree_file, "newick")
    print(f"树包含 {len(tree.get_terminals())} 个末端节点")

    # 加载物种信息
    species_info = load_species_info(clusters_file)
    color_map = get_color_map([g for _, g in species_info.values()])

    # 准备标签颜色和标签名称
    label_colors = {}
    label_names = {}
    for t in tree.get_terminals():
        name = t.name
        if name in species_info:
            full_species, group = species_info[name]
            label_colors[name] = color_map.get(group, 'black')
            # 标签格式：GCA编号 + 物种名缩写
            short_species = group.replace('V_', 'V. ')  # V_parvula -> V. parvula
            label_names[name] = f"{name}\n({short_species})"
        else:
            label_colors[name] = 'black'
            label_names[name] = name

    # 自定义标签函数
    def label_func(clade):
        if clade.is_terminal() and clade.name in label_names:
            return label_names[clade.name]
        return ""

    # 创建图形
    fig, ax = plt.subplots(figsize=(18, 12))

    # 使用Bio.Phylo.draw绘制
    Phylo.draw(tree, axes=ax, do_show=False,
               label_colors=label_colors,
               label_func=label_func,
               branch_labels=lambda c: "")

    # 获取所有文本对象并修改样式
    for text in ax.texts:
        txt = text.get_text()
        # 检查是否是我们的标签（包含GCA）
        for gca, color in label_colors.items():
            if gca in txt:
                text.set_color(color)
                text.set_fontweight('bold')
                text.set_fontsize(8)
                break

    # 设置标题
    ax.set_title('Phylogenetic Tree of Veillonella Species\n(Based on 1,125 Core Genes, IQ-TREE with UFBoot)',
                fontsize=14, fontweight='bold', pad=20)

    # 添加图例（放在右下方避免遮挡数据）
    legend_patches = [mpatches.Patch(color=color,
                                     label=species.replace('V_', 'V. '))
                      for species, color in sorted(color_map.items())]
    ax.legend(handles=legend_patches, loc='lower right',
             title='Species', frameon=True, fontsize=10,
             fancybox=True, shadow=True)

    # 调整布局
    plt.tight_layout()

    # 保存
    output_pdf = f"{output_prefix}.pdf"
    output_png = f"{output_prefix}.png"

    plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
    plt.savefig(output_png, dpi=300, bbox_inches='tight')

    print(f"[OK] 保存PDF: {output_pdf}")
    print(f"[OK] 保存PNG: {output_png}")

    plt.close()
    print("\n=== 完成 ===")


def main():
    import argparse
    parser = argparse.ArgumentParser(description="绘制最终版系统发育树")
    parser.add_argument("--tree", "-t", required=True, help="树文件(Newick)")
    parser.add_argument("--clusters", "-c", required=True, help="物种聚类CSV")
    parser.add_argument("--output", "-o", required=True, help="输出前缀")

    args = parser.parse_args()

    plot_tree_final(args.tree, args.clusters, args.output)


if __name__ == "__main__":
    main()