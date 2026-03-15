#!/usr/bin/env python3
"""
Snakemake 流程验证脚本

功能:
- 检查所有必需文件是否存在
- 验证配置文件格式
- 检查 Conda 环境
- 预览流程（dry-run）
"""

import sys
import os
from pathlib import Path
import yaml


def check_file(filepath: str, required: bool = True) -> bool:
    """检查文件是否存在"""
    exists = Path(filepath).exists()
    if exists:
        print(f"  [OK] {filepath}")
        return True
    else:
        if required:
            print(f"  [MISSING] {filepath} (必需)")
            return False
        else:
            print(f"  [-] {filepath} (可选，运行后生成)")
            return True


def validate_config() -> bool:
    """验证配置文件"""
    try:
        with open('config.yaml', 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)
        print("  [OK] config.yaml 格式正确")
        return True
    except Exception as e:
        print(f"  [ERROR] config.yaml 解析失败：{e}")
        return False


def check_conda_env() -> bool:
    """检查 Conda 环境"""
    import subprocess
    try:
        result = subprocess.run(
            ['conda', 'env', 'list'],
            capture_output=True,
            text=True,
            timeout=30
        )
        if 'veillonella' in result.stdout:
            print("  [OK] Conda 环境 'veillonella' 已安装")
            return True
        else:
            print("  [WARN] Conda 环境 'veillonella' 未安装")
            print("    运行：conda env create -f envs/veillonella.yaml")
            return False
    except Exception as e:
        print(f"  [ERROR] 无法检查 Conda 环境：{e}")
        return False


def main():
    print("=" * 60)
    print("韦荣氏菌分析流程 - 验证检查")
    print("=" * 60)
    print()
    
    all_passed = True
    
    # 1. 检查必需文件
    print("1. 检查必需文件:")
    required_files = [
        'Snakefile',
        'config.yaml',
        'envs/veillonella.yaml',
        'scripts/01_download_genomes.py',
        'scripts/02_fetch_reference.py',
        'scripts/03_blast_metabolism.py',
        'scripts/04_plot_tree_heatmap.R'
    ]
    
    for f in required_files:
        if not check_file(f):
            all_passed = False
    
    print()
    
    # 2. 检查可选文件（运行后生成）
    print("2. 检查输出目录:")
    optional_dirs = [
        'data/genomes',
        'data/reference',
        'results/annotation',
        'results/pangenome',
        'results/phylogeny',
        'results/metabolism'
    ]
    
    for d in optional_dirs:
        check_file(d, required=False)
    
    print()
    
    # 3. 验证配置
    print("3. 验证配置文件:")
    if not validate_config():
        all_passed = False
    
    print()
    
    # 4. 检查 Conda 环境
    print("4. 检查 Conda 环境:")
    if not check_conda_env():
        all_passed = False
    
    print()
    print("=" * 60)
    
    if all_passed:
        print("[OK] 所有检查通过！可以运行流程")
        print()
        print("运行命令:")
        print("  snakemake --cores 8 --use-conda")
    else:
        print("[ERROR] 部分检查未通过，请先修复问题")
        sys.exit(1)


if __name__ == '__main__':
    main()
