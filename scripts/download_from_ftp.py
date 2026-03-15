#!/usr/bin/env python3
"""
从 NCBI RefSeq FTP 下载韦荣氏菌基因组
"""
import requests
from pathlib import Path
import re

# RefSeq FTP 基础 URL
FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/11257/101/"

# 已知的韦荣氏菌基因组列表（RefSeq representative）
GENOMES = [
    {
        "name": "Veillonella_parvula_DSM_2008",
        "accession": "GCF_000172375.1",
        "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/172/375/GCF_000172375.1_ASM17237v1/GCF_000172375.1_ASM17237v1_genomic.fna.gz"
    },
    {
        "name": "Veillonella_atypica_DSM_20736",
        "accession": "GCF_000437155.1",
        "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/437/155/GCF_000437155.1_ASM43715v2/GCF_000437155.1_ASM43715v2_genomic.fna.gz"
    },
    {
        "name": "Veillonella_dispar_ATCC_17748",
        "accession": "GCF_000189315.1",
        "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/189/315/GCF_000189315.1_ASM18931v1/GCF_000189315.1_ASM18931v1_genomic.fna.gz"
    },
]

def download_genome(genome_info, output_dir):
    """下载单个基因组"""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    output_file = output_dir / f"{genome_info['name']}.fna.gz"
    
    print(f"Downloading {genome_info['name']}...")
    response = requests.get(genome_info['url'], stream=True)
    
    if response.status_code == 200:
        with open(output_file, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print(f"  ✓ Downloaded: {output_file}")
        return True
    else:
        print(f"  ✗ Failed: {response.status_code}")
        return False

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        output_dir = "data/genomes"
    else:
        output_dir = sys.argv[1]
    
    print(f"=== Downloading Veillonella genomes to {output_dir} ===\n")
    
    success_count = 0
    for genome in GENOMES:
        if download_genome(genome, output_dir):
            success_count += 1
    
    print(f"\n=== Download complete: {success_count}/{len(GENOMES)} genomes ===")
