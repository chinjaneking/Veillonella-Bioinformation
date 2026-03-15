#!/usr/bin/env Rscript
"""
脚本：绘制注释的系统发育树
功能：
- 读取Newick格式的树文件
- 根据ANI聚类结果着色分支
- 添加Bootstrap支持值
- 生成期刊级质量的图

用法：
Rscript plot_annotated_tree.R <treefile> <clusters.csv> <output_prefix>
"""

library(ggtree)
library(ggplot2)
library(ape)
library(dplyr)
library(RColorBrewer)

# 获取命令行参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("用法: Rscript plot_annotated_tree.R <treefile> <clusters.csv> <output_prefix>\n")
  quit(status = 1)
}

tree_file <- args[1]
clusters_file <- args[2]
output_prefix <- args[3]

cat("=== 绘制注释系统发育树 ===\n")
cat(sprintf("树文件: %s\n", tree_file))
cat(sprintf("聚类文件: %s\n", clusters_file))

# 读取树
tree <- read.tree(tree_file)
cat(sprintf("树包含 %d 个taxa\n", length(tree$tip.label)))

# 读取聚类信息
clusters <- read.csv(clusters_file)
cat(sprintf("读取 %d 个聚类标签\n", nrow(clusters)))

# 创建物种分组映射
species_groups <- setNames(clusters$Species_Group, clusters$Genome)

# 为每个tip添加物种信息
tree$tip.species <- species_groups[tree$tip.label]

# 定义颜色方案（使用ColorBrewer的Set2调色板）
unique_species <- unique(clusters$Species_Group)
n_species <- length(unique_species)
if (n_species <= 8) {
  colors <- brewer.pal(n_species, "Set2")
} else {
  colors <- rainbow(n_species)
}
names(colors) <- unique_species

# 基础树图
p <- ggtree(tree, layout = "rectangular", size = 0.8)

# 添加Bootstrap值（从节点数据中读取）
# 注意：需要树的节点有bootstrap信息
# 这里我们添加圆形标记表示高支持度

# 添加tip标签，按物种着色
p <- p + geom_tiplab(
  aes(color = tree$tip.species[match(label, tree$tip.label)]),
  size = 3,
  hjust = -0.1,
  align = TRUE,
  linesize = 0.2
) +
  scale_color_manual(
    name = "Species Group",
    values = colors,
    na.value = "gray50"
  )

# 添加比例尺
p <- p + theme_tree2() +
  geom_treescale(x = 0, y = -2, fontsize = 3)

# 添加标题和主题
p <- p + labs(
  title = "Phylogenetic Tree of Veillonella Species",
  subtitle = sprintf("Based on %d core genes", length(tree$tip.label))
) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

# 保存PDF
output_pdf <- paste0(output_prefix, ".pdf")
ggsave(output_pdf, p, width = 12, height = max(8, length(tree$tip.label) * 0.3), dpi = 300)
cat(sprintf("✓ 保存PDF: %s\n", output_pdf))

# 保存PNG
output_png <- paste0(output_prefix, ".png")
ggsave(output_png, p, width = 12, height = max(8, length(tree$tip.label) * 0.3), dpi = 300)
cat(sprintf("✓ 保存PNG: %s\n", output_png))

# 创建带Bootstrap值的版本（如果树文件中有）
# 读取IQ-TREE的.contree文件（包含bootstrap值）
contree_file <- sub("\\.treefile$", ".contree", tree_file)
if (file.exists(contree_file)) {
  cat("\n=== 创建带Bootstrap值的树 ===\n")

  tree_bs <- read.tree(contree_file)

  # 解析bootstrap值（存储在节点标签中）
  bootstrap_values <- as.numeric(tree_bs$node.label)
  bootstrap_values <- bootstrap_values[!is.na(bootstrap_values)]

  cat(sprintf("Bootstrap范围: %.0f - %.0f\n", min(bootstrap_values), max(bootstrap_values)))

  # 绘制带bootstrap的树
  p_bs <- ggtree(tree_bs, layout = "rectangular", size = 0.8)

  # 添加bootstrap标签（只显示>=70的）
  p_bs <- p_bs + geom_nodelab(
    aes(subset = as.numeric(label) >= 70),
    size = 2,
    color = "darkred",
    hjust = 1.2,
    vjust = -0.5
  )

  # 添加tip标签
  p_bs <- p_bs + geom_tiplab(
    aes(color = species_groups[match(label, names(species_groups))]),
    size = 3,
    hjust = -0.1,
    align = TRUE
  ) +
    scale_color_manual(
      name = "Species Group",
      values = colors,
      na.value = "gray50"
    )

  p_bs <- p_bs + theme_tree2() +
    labs(
      title = "Phylogenetic Tree with Bootstrap Support",
      subtitle = "Ultrafast Bootstrap (UFBoot) values shown for nodes ≥70%"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "right"
    )

  output_bs_pdf <- paste0(output_prefix, "_with_bootstrap.pdf")
  ggsave(output_bs_pdf, p_bs, width = 12, height = max(8, length(tree_bs$tip.label) * 0.3), dpi = 300)
  cat(sprintf("✓ 保存带Bootstrap的树: %s\n", output_bs_pdf))
}

# 创建环形树（circular layout）
cat("\n=== 创建环形树布局 ===\n")
p_circ <- ggtree(tree, layout = "circular", size = 0.6) +
  geom_tiplab(
    aes(color = species_groups[match(label, names(species_groups))]),
    size = 2.5,
    hjust = -0.1
  ) +
  scale_color_manual(
    name = "Species Group",
    values = colors,
    na.value = "gray50"
  ) +
  labs(title = "Circular Phylogenetic Tree of Veillonella") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right"
  )

output_circ <- paste0(output_prefix, "_circular.pdf")
ggsave(output_circ, p_circ, width = 10, height = 10, dpi = 300)
cat(sprintf("✓ 保存环形树: %s\n", output_circ))

cat("\n=== 完成 ===\n")
