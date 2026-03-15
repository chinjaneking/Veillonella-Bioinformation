"""
Snakefile - Veillonella 系统发育与代谢功能分析流程

极轻量化设计（磁盘占用<10GB）
Docker 模式：工具已全局安装，无需 conda 环境
"""

import os
import pandas as pd
from pathlib import Path

# 读取配置
configfile: "config.yaml"

# 动态获取样本列表（下载完成后才能确定）
def get_samples():
    """从 data/genomes 目录动态获取样本列表"""
    genomes_dir = Path("data/genomes")
    if not genomes_dir.exists():
        return []
    return [f.stem for f in genomes_dir.glob("*.fna")]

# 默认规则
rule all:
    input:
        "data/genomes/assembly_data_report.json",  # 强制先下载基因组
        "results/phylogeny/veillonella.treefile",
        "results/metabolism/metabolism_matrix.csv",
        "results/phylogeny/tree_heatmap.pdf"

# ========================================
# 阶段 1: 数据下载
# ========================================
rule download_genomes:
    output:
        zip="data/genomes/raw_datasets.zip",
        report="data/genomes/assembly_data_report.json"
    shell:
        """
        echo "=== 下载韦荣氏菌属基因组 ==="
        # 使用 ncbi-datasets 下载 Complete Genome 级别的 reference/representative 基因组
        datasets download genome taxon Veillonella \
            --assembly-level chromosome,complete \
            --include genome \
            --filename {output.zip} \
            --no-progressbar
        
        echo "=== 解压并提取 .fna 文件 ==="
        # 使用 Python 解压（更可靠）
        python -c "import zipfile; z = zipfile.ZipFile('data/genomes/raw_datasets.zip'); z.extractall('data/genomes/temp')"
        
        # 运行 Python 脚本处理
        python scripts/01_download_genomes.py data/genomes/temp data/genomes
        
        echo "=== 清理临时文件 ==="
        rm -rf data/genomes/temp
        rm -f {output.zip}
        """

# ========================================
# 阶段 2: Prokka 注释
# ========================================
rule prokka_annotate:
    input:
        "data/genomes/{sample}.fna"
    output:
        gff="results/annotation/{sample}/annotation.gff",
        faa="results/annotation/{sample}/annotation.faa",
        fna="results/annotation/{sample}/annotation.fna"
    threads: 4
    log:
        "results/annotation/{sample}/prokka.log"
    shell:
        """
        export PERL5LIB=/opt/miniforge3/lib/perl5/site_perl/5.22.0:/opt/miniforge3/lib/perl5/5.22.0
        export PATH=/opt/miniforge3/bin:$PATH
        echo "=== 注释 {wildcards.sample} ==="
        prokka {input} \
            --outdir results/annotation/{wildcards.sample} \
            --prefix annotation \
            --kingdom {config[prokka][kingdom]} \
            --cpus {threads} \
            --force \
            --quiet \
            2> {log}
        
        echo "=== 修复 GFF 序列名 ==="
        bash scripts/fix_gff.sh {wildcards.sample} results/annotation/{wildcards.sample}/annotation.gff {input}
        
        echo "=== 清理冗余文件 ==="
        cd results/annotation/{wildcards.sample}/
        rm -f *.ffn *.frn *.sqn *.err *.log *.tbl
        """

# ========================================
# 阶段 3: 泛基因组分析
# ========================================
rule panaroo:
    input:
        "results/annotation/GCA_000024945.1/annotation.gff",
        "results/annotation/GCA_002005185.1/annotation.gff",
        "results/annotation/GCA_002082765.1/annotation.gff",
        "results/annotation/GCA_013393365.1/annotation.gff",
        "results/annotation/GCA_016127175.1/annotation.gff",
        "results/annotation/GCA_018406505.1/annotation.gff",
        "results/annotation/GCA_026057735.1/annotation.gff",
        "results/annotation/GCA_028743475.1/annotation.gff",
        "results/annotation/GCA_031190775.1/annotation.gff",
        "results/annotation/GCA_032695835.1/annotation.gff",
        "results/annotation/GCA_036456085.1/annotation.gff",
        "results/annotation/GCA_042850485.1/annotation.gff",
        "results/annotation/GCA_900186885.1/annotation.gff",
        "results/annotation/GCA_900187285.1/annotation.gff",
        "results/annotation/GCA_900637515.1/annotation.gff",
        "results/annotation/GCA_902810435.1/annotation.gff"
    output:
        gpa="results/pangenome/gene_presence_absence.csv",
        core="results/pangenome/core_gene_alignment.aln",
        ref="results/pangenome/pan_genome_reference.fa"
    threads: 8
    log:
        "results/pangenome/panaroo.log"
    shell:
        """
        echo "=== 运行 Panaroo 泛基因组分析（16 个基因组）==="
        # Create file list for Panaroo (one file per line)
        for f in {input}; do echo "$f"; done > /tmp/panaroo_gff_list.txt
        
        panaroo \
            --clean-mode sensitive \
            -a core \
            --refind-mode off \
            -i /tmp/panaroo_gff_list.txt \
            -o results/pangenome \
            -t {threads} \
            2> {log}
        """
# ========================================
# 阶段 4: 系统发育树构建
# ========================================
rule trim_alignment:
    input:
        "results/pangenome/core_gene_alignment.aln"
    output:
        "results/phylogeny/core_alignment_trimmed.aln"
    shell:
        """
        echo "=== 修剪核心基因比对 ==="
        trimal \
            -in {input} \
            -{config[phylogeny][trim_method]} \
            -out {output}
        """

rule build_tree:
    input:
        "results/phylogeny/core_alignment_trimmed.aln"
    output:
        tree="results/phylogeny/veillonella.treefile",
        log="results/phylogeny/veillonella.iqtree"
    threads: 8
    shell:
        """
        echo "=== 构建最大似然树 ==="
        iqtree \
            -s {input} \
            -m {config[phylogeny][model]} \
            -B {config[phylogeny][bootstrap]} \
            -T {threads} \
            -pre results/phylogeny/veillonella \
            -nt {threads}
        
        # 仅保留 treefile，清理其他文件
        mv results/phylogeny/veillonella.treefile {output.tree}
        mv results/phylogeny/veillonella.iqtree {output.log}
        rm -f results/phylogeny/veillonella.*
        """

# ========================================
# 阶段 5: 靶向代谢功能挖掘
# ========================================
rule fetch_reference:
    output:
        ldh="data/reference/ldh_reference.faa",
        mcm="data/reference/mcm_reference.faa"
    script:
        "scripts/02_fetch_reference.py"

rule make_blast_db:
    input:
        "results/pangenome/pan_genome_reference.fa"
    output:
        "results/pangenome/pan_genome_db.pin"
    shell:
        """
        echo "=== 构建 BLAST 数据库 ==="
        makeblastdb -in {input} -dbtype prot -out results/pangenome/pan_genome_db
        """

rule blast_ldh:
    input:
        db="results/pangenome/pan_genome_db.pin",
        ref="data/reference/ldh_reference.faa"
    output:
        "results/metabolism/ldh_blast.tsv"
    shell:
        """
        echo "=== BLAST 搜索 LDH ==="
        blastp \
            -query {input.ref} \
            -db results/pangenome/pan_genome_db \
            -out {output} \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
            -evalue {config[blast][evalue]} \
            -qcov_hsp_perc {config[blast][qcov_hsp]}
        """

rule blast_mcm:
    input:
        db="results/pangenome/pan_genome_db.pin",
        ref="data/reference/mcm_reference.faa"
    output:
        "results/metabolism/mcm_blast.tsv"
    shell:
        """
        echo "=== BLAST 搜索 MCM ==="
        blastp \
            -query {input.ref} \
            -db results/pangenome/pan_genome_db \
            -out {output} \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
            -evalue {config[blast][evalue]} \
            -qcov_hsp_perc {config[blast][qcov_hsp]}
        """

rule metabolism_matrix:
    input:
        ldh="results/metabolism/ldh_blast.tsv",
        mcm="results/metabolism/mcm_blast.tsv",
        gpa="results/pangenome/gene_presence_absence.csv"
    output:
        "results/metabolism/metabolism_matrix.csv"
    script:
        "scripts/03_blast_metabolism.py"

# ========================================
# 阶段 6: 可视化
# ========================================
rule plot_tree_heatmap:
    input:
        tree="results/phylogeny/veillonella.treefile",
        matrix="results/metabolism/metabolism_matrix.csv"
    output:
        "results/phylogeny/tree_heatmap.pdf"
    script:
        "scripts/04_plot_tree_heatmap.R"

# ========================================
# 清理规则（可选手动运行）
# ========================================
rule cleanup:
    input:
        "results/phylogeny/tree_heatmap.pdf"
    shell:
        """
        echo "=== 开始清理临时文件和缓存 ==="
        
        # 清理 BLAST 数据库索引
        rm -rf results/pangenome/pan_genome_db*
        
        # 清理 IQ-TREE 临时文件
        rm -f results/phylogeny/veillonella.*
        
        # 清理 Conda 缓存
        echo "清理 Conda 缓存..."
        conda clean --all -y
        
        echo "=== 清理完成 ==="
        """