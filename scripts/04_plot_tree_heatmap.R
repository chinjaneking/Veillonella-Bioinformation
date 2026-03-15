#!/usr/bin/env Rscript
# 脚本 04: 系统发育树 + 代谢基因热图联画
#
# 功能:
# - 读取系统发育树 (Newick 格式)
# - 读取代谢基因存在/缺失矩阵
# - 使用 ggtree 绘制左侧树 + 右侧热图
# - 输出 PDF 和 PNG 格式

# 加载依赖
library(ggtree)
library(tidyverse)
library(pheatmap)
library(ape)
library(readr)

# ========== 配置 ==========
TREE_FILE <- "results/phylogeny/veillonella.treefile"
MATRIX_FILE <- "results/metabolism/metabolism_matrix.csv"
OUTPUT_PDF <- "results/phylogeny/tree_heatmap.pdf"
OUTPUT_PNG <- "results/phylogeny/tree_heatmap.png"

# ========== 读取数据 ==========
cat("=== 读取输入文件 ===\n")

# 读取系统发育树
if (!file.exists(TREE_FILE)) {
  stop(paste("错误：未找到树文件:", TREE_FILE))
}
tree <- read.tree(TREE_FILE)
cat(sprintf("✓ 读取树文件：%d 个 tip\n", length(tree$tip.label)))

# 读取代谢矩阵
if (!file.exists(MATRIX_FILE)) {
  stop(paste("错误：未找到代谢矩阵文件:", MATRIX_FILE))
}
matrix_df <- read_csv(MATRIX_FILE, show_col_types = FALSE)
cat(sprintf("✓ 读取代谢矩阵：%d 个物种\n", nrow(matrix_df)))

# 打印矩阵内容
cat("\n代谢矩阵内容:\n")
print(matrix_df)

# ========== 数据预处理 ==========
cat("\n=== 数据预处理 ===\n")

# 设置行名
rownames(matrix_df) <- matrix_df$species
matrix_data <- matrix_df[, -1]  # 移除 species 列

# 确保矩阵中的物种名与树中的 tip 标签匹配
# 尝试多种匹配策略
match_species <- function(tree_tips, matrix_species) {
  # 策略 1: 完全匹配
  matched <- match(tree_tips, matrix_species)
  if (all(!is.na(matched))) {
    return(tree_tips)
  }
  
  # 策略 2: 部分匹配（矩阵物种名包含树的 tip 名）
  result <- sapply(tree_tips, function(tip) {
    hits <- grep(tip, matrix_species, fixed = TRUE)
    if (length(hits) > 0) {
      return(matrix_species[hits[1]])
    }
    return(NA)
  })
  
  return(result)
}

# 尝试匹配
matched_species <- match_species(tree$tip.label, rownames(matrix_data))

# 重命名矩阵行以匹配树
if (!all(is.na(matched_species))) {
  # 创建映射
  species_map <- setNames(names(matched_species), matched_species)
  
  # 重命名矩阵行
  new_rownames <- sapply(rownames(matrix_data), function(x) {
    if (x %in% names(species_map)) {
      return(species_map[x])
    }
    return(x)
  })
  rownames(matrix_data) <- new_rownames
}

# 排序矩阵以匹配树的 tip 顺序
tree_tip_order <- tree$tip.label
matrix_ordered <- matrix_data[tree_tip_order, , drop = FALSE]

# 移除 NA 行
matrix_ordered <- matrix_ordered[!is.na(rownames(matrix_ordered)), ]

cat(sprintf("✓ 匹配后矩阵：%d 行 %d 列\n", nrow(matrix_ordered), ncol(matrix_ordered)))

# ========== 绘图 ==========
cat("\n=== 绘制树 + 热图 ===\n")

# 基础树
p <- ggtree(tree, layout = "rectangular") +
  geom_tiplab(size = 3, hjust = -0.1, align = TRUE) +
  theme_tree2() +
  theme(
    legend.position = "right",
    plot.margin = margin(10, 10, 10, 10)
  )

# 添加热图
# 注意：gheatmap 要求矩阵行名与树的 tip 标签完全匹配
tryCatch({
  p2 <- gheatmap(
    p,
    matrix_ordered,
    offset = 0.1,
    width = 0.6,
    colnames = TRUE,
    colnames_angle = 0,
    colnames_offset_y = 0,
    hjust = 0,
    font.size = 3
  ) +
    scale_fill_manual(
      name = "基因存在",
      values = c("0" = "#CCCCCC", "1" = "#E41A1C"),
      labels = c("缺失", "存在")
    )
  
  # 保存 PDF
  ggsave(
    OUTPUT_PDF,
    p2,
    width = 10,
    height = 8,
    units = "in"
  )
  cat(sprintf("✓ 保存 PDF: %s\n", OUTPUT_PDF))
  
  # 保存 PNG
  ggsave(
    OUTPUT_PNG,
    p2,
    width = 10,
    height = 8,
    units = "in",
    dpi = 300
  )
  cat(sprintf("✓ 保存 PNG: %s\n", OUTPUT_PNG))
  
}, error = function(e) {
  cat(sprintf("警告：gheatmap 失败：%s\n", e$message))
  cat("尝试替代方案...\n")
  
  # 替代方案：分别绘制树和热图
  p_tree <- ggtree(tree) +
    geom_tiplab(size = 3) +
    theme_tree2()
  
  ggsave(
    gsub(".pdf", "_tree_only.pdf", OUTPUT_PDF),
    p_tree,
    width = 8,
    height = 6
  )
  
  # 绘制热图
  pdf(gsub(".pdf", "_heatmap_only.pdf", OUTPUT_PDF), width = 6, height = 8)
  pheatmap(
    matrix_ordered,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = TRUE,
    color = c("0" = "#CCCCCC", "1" = "#E41A1C")
  )
  dev.off()
  
  cat("已保存分开的树和热图文件\n")
})

cat("\n=== 完成 ===\n")
