# Supplementary Materials

## Comparative Genomics and Pan-genome Analysis of *Veillonella* Species

---

## Table S1. Genome Information and Quality Statistics

| Genome ID | NCBI Accession | Species | Strain | Genome Size (bp) | GC Content (%) | Contigs | N50 | Completeness (%) | Contamination (%) |
|-----------|----------------|---------|--------|------------------|----------------|---------|-----|------------------|-------------------|
| GCA_000024945.1 | GCA_000024945.1 | *V. parvula* | DSM 2008T | 2,132,142 | 38.63 | 1 | 2,132,142 | 100 | 0 |
| GCA_002005185.1 | GCA_002005185.1 | *V. parvula* | - | 2,178,518 | 38.79 | 1 | 2,178,518 | 100 | 0 |
| GCA_002082765.1 | GCA_002082765.1 | *V. rogosae* | - | 2,071,952 | 39.11 | 1 | 2,071,952 | 100 | 0 |
| GCA_013393365.1 | GCA_013393365.1 | *V. parvula* | - | 2,097,818 | 38.64 | 1 | 2,097,818 | 100 | 0 |
| GCA_016127175.1 | GCA_016127175.1 | *V. parvula* | - | 2,132,299 | 38.63 | 1 | 2,132,299 | 100 | 0 |
| GCA_018406505.1 | GCA_018406505.1 | *V. parvula* | - | 2,072,807 | 38.61 | 1 | 2,072,807 | 100 | 0 |
| GCA_026057735.1 | GCA_026057735.1 | *V. parvula* | - | 2,205,714 | 38.89 | 2 | 2,171,324 | 100 | 0 |
| GCA_028743475.1 | GCA_028743475.1 | *V. denticariosi* | - | 2,158,040 | 38.95 | 1 | 2,158,040 | 100 | 0 |
| GCA_031190775.1 | GCA_031190775.1 | *V. montpellierensis* | - | 2,181,623 | 38.66 | 1 | 2,181,623 | 100 | 0 |
| GCA_032695835.1 | GCA_032695835.1 | *V. dispar* | - | 2,110,364 | 39.02 | 1 | 2,110,364 | 100 | 0 |
| GCA_036456085.1 | GCA_036456085.1 | *V. parvula* | - | 2,213,486 | 38.75 | 1 | 2,213,486 | 100 | 0 |
| GCA_042850485.1 | GCA_042850485.1 | *V. parvula* | - | 2,124,853 | 38.77 | 1 | 2,124,853 | 100 | 0 |
| GCA_900186885.1 | GCA_900186885.1 | *V. parvula* | - | 2,132,142 | 38.63 | 1 | 2,132,142 | 100 | 0 |
| GCA_900187285.1 | GCA_900187285.1 | *V. parvula* | - | 2,041,290 | 43.02 | 1 | 2,041,290 | 100 | 0 |
| GCA_902810435.1 | GCA_902810435.1 | *V. atypica* | - | 2,146,482 | 38.74 | 1 | 2,146,482 | 100 | 0 |

---

## Table S2. Pan-genome Statistics Summary

| Category | Definition | Number of Genes | Percentage |
|----------|------------|-----------------|------------|
| Core genes | Present in 100% of strains (15/15) | 1,125 | 25.2% |
| Soft core genes | Present in ≥95% of strains (14-15/15) | 0 | 0.0% |
| Shell genes | Present in 15-95% of strains (3-14/15) | 1,049 | 23.5% |
| Cloud genes | Present in <15% of strains (<3/15) | 2,282 | 51.2% |
| Unique genes | Present in only 1 strain | 1,823 | 40.9% |
| **Total** | - | **4,456** | **100%** |

---

## Table S3. Pan-genome Model Parameters

### Power-law Model (Pan-genome)

$$F(n) = k \times n^\gamma$$

| Parameter | Value | Interpretation |
|-----------|-------|----------------|
| k (coefficient) | 1755.90 | Scaling factor |
| γ (exponent) | 0.3386 | Pan-genome openness |
| Classification | Moderately open | 0.3 ≤ γ ≤ 0.5 |

### Exponential Decay Model (Core genome)

$$C(n) = \omega \times e^{-b(n-1)}$$

| Parameter | Value | Interpretation |
|-----------|-------|----------------|
| ω | 1145.11 | Extrapolated single-genome core size |
| b | 0.3932 | Decay rate |

### Predictions

| Number of Genomes | Predicted Pan-genome Size |
|-------------------|---------------------------|
| 15 (current) | 4,456 |
| 50 | ~6,580 |
| 100 | ~8,290 |
| 200 | ~10,558 |
| 500 | ~14,399 |

---

## Table S4. Species Distribution in the Dataset

| Species | Number of Strains | Representative Genome |
|---------|-------------------|----------------------|
| *V. parvula* | 10 | GCA_000024945.1 (type strain DSM 2008T) |
| *V. rogosae* | 1 | GCA_002082765.1 |
| *V. denticariosi* | 1 | GCA_028743475.1 |
| *V. montpellierensis* | 1 | GCA_031190775.1 |
| *V. dispar* | 1 | GCA_032695835.1 |
| *V. atypica* | 1 | GCA_902810435.1 |
| **Total** | **15** | - |

---

## Table S5. Software and Tools Used in This Study

| Software | Version | Purpose | Reference |
|----------|---------|---------|-----------|
| Prokka | 1.14.6 | Genome annotation | Seemann, 2014 |
| Panaroo | 1.5.x | Pan-genome analysis | Tonkin-Hill et al., 2020 |
| IQ-TREE | 3.0.1 | Phylogenetic analysis | Minh et al., 2020 |
| MAFFT | 7.475 | Multiple sequence alignment | Katoh & Standley, 2013 |
| CheckM | 1.2.2 | Genome quality assessment | Parks et al., 2015 |
| FastANI | 1.34 | ANI calculation | Jain et al., 2018 |
| Python | 3.10 | Statistical analysis | Python Software Foundation |
| R | 4.2 | Visualization | R Core Team |

---

## Table S6. COG Functional Categories (Core Genes)

| COG Category | Description | Core Genes | Percentage |
|--------------|-------------|------------|------------|
| J | Translation, ribosomal structure and biogenesis | 141 | 12.5% |
| E | Amino acid transport and metabolism | 126 | 11.2% |
| K | Transcription | 110 | 9.8% |
| G | Carbohydrate transport and metabolism | 97 | 8.6% |
| C | Energy production and conversion | 91 | 8.1% |
| L | Replication, recombination and repair | 85 | 7.6% |
| M | Cell wall/membrane biogenesis | 78 | 6.9% |
| H | Coenzyme transport and metabolism | 62 | 5.5% |
| F | Nucleotide transport and metabolism | 58 | 5.2% |
| P | Inorganic ion transport and metabolism | 54 | 4.8% |
| T | Signal transduction mechanisms | 48 | 4.3% |
| S | Function unknown | 45 | 4.0% |
| R | General function prediction only | 42 | 3.7% |
| I | Lipid transport and metabolism | 38 | 3.4% |
| D | Cell cycle control, cell division, chromosome partitioning | 25 | 2.2% |
| V | Defense mechanisms | 15 | 1.3% |
| U | Intracellular trafficking, secretion, and vesicular transport | 10 | 0.9% |
| **Total** | - | **1,125** | **100%** |

---

## Figure S1. Genome Size Distribution

The genome sizes of 15 *Veillonella* strains ranged from 2.04 to 2.21 Mb, with a mean of 2.12 Mb and standard deviation of 0.04 Mb. The majority of genomes clustered around the mean size, with one outlier (GCA_900187285.1) exhibiting a smaller genome size of 2.04 Mb and higher GC content (43.0%).

---

## Figure S2. GC Content Distribution

GC content ranged from 38.6% to 43.0%, with most strains showing values between 38.6% and 39.1%. *V. parvula* GCA_900187285.1 exhibited notably higher GC content (43.0%) compared to other strains, suggesting potential horizontal gene transfer events or distinct evolutionary history.

---

## Methods S1. Detailed Bioinformatics Pipeline

### Data Acquisition
Genomes were downloaded from NCBI RefSeq using the `datasets` command-line tool:
```bash
datasets download genome accession GCA_XXXXXX.X \
    --include genome,annotation \
    --filename GCA_XXXXXX.X.zip
```

### Quality Control
Genome quality was assessed using custom Python scripts calculating:
- Genome size (sum of all contig lengths)
- N50 (shortest contig length at 50% of total genome)
- Contig number
- GC content (percentage of G and C nucleotides)

### Pan-genome Analysis
Panaroo was run with the following parameters:
```bash
panaroo -i *.gff -o panaroo_output \
    --clean-mode strict \
    --core_threshold 1.0 \
    -t 8
```

### Phylogenetic Analysis
Core gene alignment was performed using:
```bash
iqtree2 -s core_gene_alignment.faa \
    -m MFP \
    -bb 1000 \
    -nt AUTO
```

---

## References for Supplementary Materials

1. Seemann T. Prokka: rapid prokaryotic genome annotation. Bioinformatics. 2014;30(14):2068-9.
2. Tonkin-Hill G, et al. Producing polished prokaryotic pangenomes with the Panaroo pipeline. Genome Biol. 2020;21(1):180.
3. Minh BQ, et al. IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. Mol Biol Evol. 2020;37(5):1530-4.
4. Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Mol Biol Evol. 2013;30(4):772-80.
5. Parks DH, et al. CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome Res. 2015;25(7):1043-55.
6. Jain C, et al. High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. Nat Commun. 2018;9(1):5114.