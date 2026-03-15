# Veillonella Comparative Genomics

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.10+-green.svg)](https://www.python.org/)
[![Snakemake](https://img.shields.io/badge/Snakemake-8.x-orange.svg)](https://snakemake.readthedocs.io/)

Comparative genomics and pan-genome analysis of *Veillonella* species.

## Overview

*Veillonella* species are anaerobic, Gram-negative cocci that are ubiquitous inhabitants of the human oral cavity and gastrointestinal tract. This project provides a comprehensive comparative genomic analysis of 15 *Veillonella* strains representing 6 species, revealing a moderately open pan-genome structure with substantial genetic diversity.

## Key Findings

| Metric | Value |
|--------|-------|
| Pan-genome size | 4,456 gene families |
| Core genes | 1,125 (25.2%) |
| Shell genes | 1,049 (23.5%) |
| Cloud genes | 2,282 (51.2%) |
| Pan-genome openness (γ) | 0.3386 (moderately open) |
| Species analyzed | 6 |
| Genomes analyzed | 15 |

## Project Structure

```
Veillonella-Bioinformation/
├── config_paper.yaml          # Analysis configuration
├── Snakefile_paper            # Snakemake workflow
├── Dockerfile_paper           # Docker environment
│
├── scripts/                   # Analysis scripts
│   ├── qc_basic.py           # Quality control
│   ├── pangenome_statistics.py
│   ├── plot_tree_final.py    # Phylogenetic tree visualization
│   └── generate_supplementary.py
│
├── results/                   # Analysis results
│   ├── qc/                   # QC reports
│   ├── pangenome/            # Pan-genome analysis
│   │   ├── gene_presence_absence.csv
│   │   └── pangenome_statistics.json
│   └── phylogeny/            # Phylogenetic analysis
│       ├── veillonella.treefile
│       └── species_clusters.csv
│
├── figures/                   # Publication figures
│   ├── Figure1_Phylogeny_v4.pdf
│   └── Figure2_Pangenome.pdf
│
├── supplementary/             # Supplementary tables
│   ├── Table_S1_Genome_Info.xlsx
│   ├── Table_S2_Pangenome_Matrix.xlsx
│   └── Table_S6_COG_Summary.xlsx
│
├── manuscript/                # Publication manuscript
│   ├── Veillonella_Pangenome_Manuscript.md
│   ├── Supplementary_Materials.md
│   ├── Cover_Letter.md
│   └── figures/
│
└── data/                      # Genome data (not tracked)
    └── genomes/
```

## Species Analyzed

| Species | Strains | Representative Genome |
|---------|---------|----------------------|
| *V. parvula* | 10 | GCA_000024945.1 (DSM 2008T) |
| *V. rogosae* | 1 | GCA_002082765.1 |
| *V. denticariosi* | 1 | GCA_028743475.1 |
| *V. montpellierensis* | 1 | GCA_031190775.1 |
| *V. dispar* | 1 | GCA_032695835.1 |
| *V. atypica* | 1 | GCA_902810435.1 |

## Methods

### Workflow

1. **Data Collection**: Genomes downloaded from NCBI RefSeq
2. **Quality Control**: Assembly statistics assessment
3. **Annotation**: Prokka v1.14.6
4. **Pan-genome Analysis**: Panaroo v1.5.x (sensitive mode)
5. **Phylogenetic Analysis**: IQ-TREE v3.0.1 (MFP, 1000 UFBoot)

### Pan-genome Modeling

Power-law model: `F(n) = 1755.90 × n^0.3386`

The exponent γ = 0.3386 indicates a **moderately open pan-genome**, suggesting ongoing horizontal gene transfer and environmental adaptation.

## Requirements

### Software

- Python 3.10+
- Snakemake 8.x
- Prokka 1.14.6
- Panaroo 1.5.x
- IQ-TREE 3.0.1
- MAFFT 7.475

### Python Packages

```
biopython
pandas
numpy
scipy
matplotlib
seaborn
```

## Quick Start

### Using Docker (Recommended)

```bash
# Build Docker image
docker build -f Dockerfile_paper -t veillonella-paper:latest .

# Run analysis
docker run --rm -v $(pwd):/data veillonella-paper:latest snakemake --cores 8
```

### Using Conda

```bash
# Create environment
conda env create -f environment.yml
conda activate veillonella

# Run workflow
snakemake --cores 8 --use-conda
```

### Manual Execution

```bash
# 1. Quality control
python scripts/qc_basic.py --input data/genomes/ --output results/qc/

# 2. Pan-genome analysis (requires Panaroo)
panaroo -i annotations/*.gff -o results/pangenome/ --clean-mode sensitive -t 8

# 3. Pan-genome statistics
python scripts/pangenome_statistics.py --input results/pangenome/ --output results/pangenome/

# 4. Phylogenetic tree
iqtree2 -s core_alignment.faa -m MFP -bb 1000 -nt AUTO

# 5. Tree visualization
python scripts/plot_tree_final.py -t results/phylogeny/veillonella.treefile -c results/phylogeny/species_clusters.csv -o figures/Figure1_Phylogeny
```

## Results Summary

### Pan-genome Composition

```
Core genes (100%):     1,125 (25.2%)  ████████████████████████
Shell genes (15-95%):  1,049 (23.5%)  ██████████████████████
Cloud genes (<15%):    2,282 (51.2%)  ████████████████████████████████████████████████████
```

### Genome Statistics

| Metric | Range | Mean |
|--------|-------|------|
| Genome size (Mb) | 2.04 - 2.21 | 2.12 |
| GC content (%) | 38.6 - 43.0 | 39.1 |
| Completeness (%) | 100 | 100 |
| Contamination (%) | 0 | 0 |

## Citation

If you use this code or data, please cite:

```
chinjaneking (2026). Comparative Genomics and Pan-genome Analysis of Veillonella Species.
GitHub repository: https://github.com/chinjaneking/Veillonella-Bioinformation
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Author

**chinjaneking**
- Affiliation: Neuma Biotech (Suzhou) Ltd.
- Email: chinjaneking@gmail.com

## Acknowledgments

- NCBI RefSeq for genome data
- The developers of Prokka, Panaroo, IQ-TREE, and other open-source tools

## Data Availability

All genome sequences are publicly available from NCBI RefSeq. Accession numbers are listed in the supplementary materials.