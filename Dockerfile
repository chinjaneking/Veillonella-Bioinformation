# 韦荣氏菌分析流程 - Docker 镜像（最终版）
# 分步安装 R 包以避免依赖问题

FROM ubuntu:22.04

LABEL maintainer="Veillonella Project"
LABEL description="Veillonella phylogeny and metabolism analysis"
LABEL version="6.0"

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Asia/Shanghai
ENV PATH="/opt/miniforge3/bin:$PATH"

WORKDIR /app

# 1. 安装系统依赖
RUN apt-get update && apt-get install -y \
    wget curl bzip2 ca-certificates git tzdata unzip \
    perl cpanminus \
    libncurses5-dev libbz2-dev zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# 2. 安装 Miniforge3
RUN wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
    && bash Miniforge3-Linux-x86_64.sh -b -p /opt/miniforge3 \
    && rm Miniforge3-Linux-x86_64.sh
RUN wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
    && bash Miniforge3-Linux-x86_64.sh -b -p /opt/miniforge3 \
    && rm Miniforge3-Linux-x86_64.sh

# 3. 配置 conda
RUN conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda config --set channel_priority flexible

# 4. 安装 mamba
RUN conda install -y mamba

# 5. 安装 Python 和生信工具
RUN mamba install -y \
    python=3.10 \
    prokka \
    panaroo \
    mafft \
    trimal \
    iqtree \
    blast \
    hmmer \
    ncbi-datasets-cli \
    pandas \
    requests \
    pyyaml \
    snakemake \
    && conda clean --all -y

# 5.1 修复 Prokka 版本比较 bug
# Bug: Perl 字符串比较 "2.17" < "2.2" 为 true (字符级比较)
# Fix: 设置 MINVER 为 0 绕过版本检查
RUN sed -i "s/MINVER  => \"2.2\",/MINVER  => \"0.1\",  # patched: bypass version check/g" /opt/miniforge3/bin/prokka \
    && sed -i "s/MINVER  => \"2.2\"/MINVER  => \"0.1\"  # patched/g" /opt/miniforge3/bin/prokka

# 5.2 安装 Perl XML 模块（通过 conda）
RUN mamba install -y perl-xml-simple perl-xml-parser && conda clean --all -y

# 6. 安装 R
RUN mamba install -y r-base r-tidyverse r-ape r-readr r-remotes && conda clean --all -y

# 7. 安装 R 包（ggtree 从 Bioconductor，pheatmap 从 CRAN）
RUN R -e "remotes::install_bioc('ggtree', ask=FALSE)"
RUN R -e "install.packages('pheatmap', repos='https://cloud.r-project.org/', quiet=TRUE)"

# 8. 复制项目文件
COPY Snakefile config.yaml ./
COPY scripts/ scripts/
COPY envs/ envs/

# 9. 创建目录
RUN mkdir -p data/genomes data/reference results/annotation results/pangenome results/phylogeny results/metabolism

ENTRYPOINT ["snakemake"]
CMD ["--cores", "4", "--use-conda"]