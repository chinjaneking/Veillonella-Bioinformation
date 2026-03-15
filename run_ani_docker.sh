#!/bin/bash
# ANI分析脚本（Docker版本）

cd /app

echo "=== 运行 fastANI 分析 ==="

# 获取基因组列表
GENOMES=$(ls data/genomes/*.fna)
N=$(echo $GENOMES | wc -w)
echo "基因组数量: $N"

# 创建查询列表
mkdir -p results/ani
for g in $GENOMES; do
    echo $g
done > results/ani/query_list.txt

# 运行fastANI
echo "运行 fastANI..."
fastANI \
    --ql results/ani/query_list.txt \
    --rl results/ani/query_list.txt \
    -o results/ani/ani_results.tsv \
    -t 4

echo "=== fastANI 完成 ==="
