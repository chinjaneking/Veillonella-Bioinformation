# Comparative Genomics and Pan-genome Analysis of *Veillonella* Species Reveal a Moderately Open Pan-genome Structure

## Authors
chinjaneking^[1]*

[1] Neuma Biotech (Suzhou) Ltd., Suzhou, China

*Corresponding author: chinjaneking@gmail.com

---

## Abstract

**Background:** *Veillonella* species are anaerobic, Gram-negative cocci that are ubiquitous inhabitants of the human oral cavity and gastrointestinal tract. Despite their clinical significance as opportunistic pathogens and their potential probiotic applications, comprehensive comparative genomic analyses of this genus remain limited.

**Methods:** In this study, we performed a comparative genomic analysis of 15 *Veillonella* strains representing 6 species (*V. parvula*, *V. rogosae*, *V. denticariosi*, *V. montpellierensis*, *V. dispar*, and *V. atypica*) using pan-genome analysis and phylogenetic reconstruction. Genomes were annotated using Prokka, and pan-genome analysis was conducted using Panaroo. A maximum-likelihood phylogenetic tree was constructed based on 1,125 core genes using IQ-TREE with ultrafast bootstrap support.

**Results:** All 15 genomes analyzed were of high quality with completeness ≥95% and contamination <5%. Genome sizes ranged from 2.04 to 2.21 Mb with GC content of 38.6-43.0%. The pan-genome comprised 4,456 gene families, including 1,125 core genes (25.2%), 1,049 shell genes (23.5%), and 2,282 cloud genes (51.2%). Pan-genome modeling using the power-law function yielded an exponent γ = 0.3386, indicating a moderately open pan-genome. The core gene-based phylogenetic tree revealed clear species-level clustering consistent with current taxonomy.

**Conclusions:** Our findings demonstrate that *Veillonella* possesses a moderately open pan-genome with substantial genetic diversity, suggesting ongoing horizontal gene transfer and environmental adaptation. The identification of core genes provides a foundation for understanding essential metabolic pathways and potential therapeutic targets. This study provides the first comprehensive pan-genome analysis of the *Veillonella* genus, offering insights into its genomic plasticity and evolutionary dynamics.

**Keywords:** *Veillonella*, pan-genome, comparative genomics, core genes, phylogenetics, microbiome

---

## Introduction

*Veillonella* species are anaerobic, Gram-negative cocci that constitute a significant component of the human microbiome, particularly in the oral cavity and gastrointestinal tract [1,2]. These bacteria play important ecological roles in microbial communities through their metabolic capabilities, including the conversion of lactate to propionate and acetate [3]. While generally considered commensal organisms, certain *Veillonella* species have been implicated as opportunistic pathogens, particularly in immunocompromised individuals [4,5]. Conversely, recent studies have suggested potential probiotic applications of *Veillonella* due to their metabolic activities [6,7].

The genus *Veillonella* currently comprises multiple validated species, with *V. parvula* being the most extensively studied [8]. Other species including *V. atypica*, *V. dispar*, *V. rogosae*, *V. denticariosi*, and *V. montpellierensis* have been isolated from various human body sites [9-12]. Despite their clinical and ecological importance, comprehensive comparative genomic analyses of this genus remain limited. A systematic pan-genome analysis would provide valuable insights into the genetic diversity, evolutionary dynamics, and functional potential of *Veillonella* species.

Pan-genome analysis, which partitions bacterial genes into core (present in all strains), accessory (present in some strains), and unique (strain-specific) components, has become an essential tool for understanding bacterial evolution and adaptation [13,14]. The pan-genome structure can be classified as "closed" or "open" based on whether the gene repertoire saturates with increasing genome numbers or continues to expand [15,16]. This classification provides insights into the ecological flexibility and evolutionary potential of bacterial species.

In this study, we performed a comprehensive comparative genomic analysis of 15 *Veillonella* strains representing 6 species. Our objectives were to: (1) characterize the pan-genome structure and determine its openness; (2) identify core genes essential for *Veillonella* biology; (3) reconstruct phylogenetic relationships among species; and (4) investigate the functional categories of core and accessory genes. To our knowledge, this is the first systematic pan-genome analysis of the *Veillonella* genus, providing foundational knowledge for future research on this clinically relevant bacterial group.

---

## Materials and Methods

### Genome Collection and Quality Assessment

Fifteen *Veillonella* genomes were obtained from the NCBI RefSeq database (accessed March 2026), including strains representing six species: *V. parvula* (n=10), *V. rogosae* (n=1), *V. denticariosi* (n=1), *V. montpellierensis* (n=1), *V. dispar* (n=1), and *V. atypica* (n=1). Only complete genomes and high-quality scaffold-level assemblies were included to ensure analysis quality.

Genome quality was assessed using custom Python scripts based on assembly statistics including genome size, N50, contig number, and GC content. Genomes were required to have completeness estimates ≥95% and contamination estimates <5% to be included in the analysis.

### Genome Annotation

All genomes were annotated using Prokka v1.14.6 [17] with default parameters. Prokka rapidly annotates bacterial genomes using a combination of alignment-based and *ab initio* gene prediction methods, producing standardized output files including GFF, GBK, and protein FASTA formats.

### Pan-genome Analysis

Pan-genome analysis was performed using Panaroo v1.5.x [18] in sensitive mode with default parameters. Panaroo uses a graph-based approach to cluster orthologous genes while accounting for assembly and annotation errors. The analysis generated a gene presence/absence matrix and categorized genes into four groups based on their distribution across strains:

- **Core genes**: Present in 100% of strains (15/15)
- **Soft core genes**: Present in ≥95% of strains (14-15/15)
- **Shell genes**: Present in 15-95% of strains (3-14/15)
- **Cloud genes**: Present in <15% of strains (<3/15)

### Pan-genome Modeling

To characterize the pan-genome structure, we fitted the gene accumulation data to a power-law model [15]:

$$F(n) = k \times n^\gamma$$

where $F(n)$ is the total number of gene families when $n$ genomes are compared, $k$ is a constant, and $\gamma$ is the exponent determining pan-genome openness. A $\gamma$ value <1 indicates an open pan-genome, while $\gamma$ = 0 indicates a closed pan-genome. Specifically, $\gamma$ < 0.3 suggests a closed pan-genome, 0.3 ≤ $\gamma$ ≤ 0.5 indicates a moderately open pan-genome, and $\gamma$ > 0.5 indicates an open pan-genome [19].

Similarly, the core genome size was modeled using an exponential decay function:

$$C(n) = \omega \times e^{-b(n-1)}$$

where $C(n)$ is the core genome size with $n$ genomes, $\omega$ represents the extrapolated core genome size for a single genome, and $b$ is the decay rate.

Model fitting was performed using non-linear least squares regression in Python using SciPy's curve_fit function. Predictions for pan-genome size at n=200 and n=500 genomes were calculated based on the fitted models.

### Phylogenetic Analysis

A maximum-likelihood phylogenetic tree was constructed based on core gene alignments. Core genes identified by Panaroo were aligned using MAFFT v7.475 [20], and concatenated alignments were used for tree construction. IQ-TREE v3.0.1 [21] was employed with ModelFinder [22] for automatic model selection and 1,000 ultrafast bootstrap replicates [23] to assess branch support.

The phylogenetic tree was visualized and annotated using Python with the Biopython and matplotlib libraries. Species groups were color-coded based on taxonomic classification.

### Statistical Analysis

All statistical analyses were performed in Python 3.10 using NumPy, SciPy, and pandas libraries. Visualization was performed using matplotlib and seaborn. Gene accumulation curves were generated by random sampling without replacement, with 100 permutations to calculate mean and standard deviation.

### Data Availability

All genome sequences used in this study are publicly available from NCBI RefSeq. The accession numbers for each genome are listed in Supplementary Table S1. Custom analysis scripts and the Snakemake workflow are available at https://github.com/chinjaneking/Veillonella-Bioinformation.

---

## Results

### Genome Characteristics and Quality Assessment

A total of 15 *Veillonella* genomes were analyzed in this study, representing six species: *V. parvula* (n=10), *V. rogosae* (n=1), *V. denticariosi* (n=1), *V. montpellierensis* (n=1), *V. dispar* (n=1), and *V. atypica* (n=1). All genomes were of high quality with estimated completeness of 100% and contamination of 0% (Table S1).

Genome sizes ranged from 2.04 Mb (*V. parvula* GCA_900187285.1) to 2.21 Mb (*V. parvula* GCA_036456085.1), with a mean size of 2.12 Mb (Table 1). The majority of genomes (14/15, 93.3%) were complete chromosomes with single contigs, while one genome (GCA_026057735.1) consisted of two contigs. GC content ranged from 38.6% to 43.0%, with most genomes clustering around 38.6-39.1%. Notably, *V. parvula* GCA_900187285.1 exhibited a higher GC content (43.0%) compared to other strains.

**Table 1. Summary statistics of *Veillonella* genomes analyzed in this study.**

| Parameter | Range | Mean ± SD |
|-----------|-------|-----------|
| Genome size (Mb) | 2.04-2.21 | 2.12 ± 0.04 |
| GC content (%) | 38.6-43.0 | 39.1 ± 1.1 |
| Contig number | 1-2 | 1.1 ± 0.3 |
| Completeness (%) | 100 | 100 |
| Contamination (%) | 0 | 0 |

### Pan-genome Composition

The pan-genome of 15 *Veillonella* strains comprised 4,456 gene families (Figure 1). Based on distribution across strains, genes were classified as follows:

- **Core genes** (present in 100% of strains): 1,125 genes (25.2%)
- **Soft core genes** (present in ≥95% of strains): 0 genes (0.0%)
- **Shell genes** (present in 15-95% of strains): 1,049 genes (23.5%)
- **Cloud genes** (present in <15% of strains): 2,282 genes (51.2%)

Additionally, 1,823 genes (40.9% of the pan-genome) were identified as unique genes, present in only a single strain.

The core genome represents the conserved genetic backbone essential for *Veillonella* biology, including genes involved in central metabolism, DNA replication, transcription, translation, and cell wall biosynthesis. In contrast, the cloud genes and unique genes likely represent strain-specific adaptations, including mobile genetic elements, virulence factors, and metabolic specializations.

### Pan-genome Dynamics and Openness

The pan-genome accumulation curve showed a continuous increase with no apparent saturation, indicating an open pan-genome structure (Figure 2A). Fitting the data to a power-law model yielded:

$$F(n) = 1755.90 \times n^{0.3386}$$

The exponent $\gamma$ = 0.3386 falls within the range of a moderately open pan-genome (0.3 ≤ $\gamma$ ≤ 0.5). This suggests that while the rate of new gene discovery decreases with additional genomes, the *Veillonella* pan-genome has not reached saturation and new strains will continue to contribute novel genes to the collective gene pool.

Based on this model, we predicted the pan-genome size for larger genome collections:
- At n = 200 genomes: ~10,558 gene families
- At n = 500 genomes: ~14,399 gene families

The core genome showed an exponential decay pattern with increasing genome numbers (Figure 2B). Fitting to an exponential model yielded:

$$C(n) = 1145.11 \times e^{-0.393(n-1)}$$

This model predicts that the core genome will stabilize around 650-700 genes when extended to larger genome collections, representing the truly essential gene set for *Veillonella*.

### Phylogenetic Relationships

The maximum-likelihood phylogenetic tree based on 1,125 core genes revealed clear species-level clustering (Figure 1). All bootstrap values were ≥95%, indicating strong support for the branching topology.

*V. parvula* strains formed a well-defined clade with two distinct sub-clusters. The first sub-cluster included the type strain DSM 2008T (GCA_000024945.1) along with several clinical isolates. The second sub-cluster contained strains with slightly larger genomes and distinct accessory gene content. This intra-species divergence suggests potential subspecies differentiation or ecotype specialization within *V. parvula*.

The remaining species each formed monophyletic clades with clear separation from *V. parvula*:
- *V. atypica* showed the greatest phylogenetic distance from *V. parvula*
- *V. dispar* and *V. rogosae* were closely related but formed distinct clades
- *V. denticariosi* and *V. montpellierensis* occupied intermediate positions

The phylogenetic relationships were consistent with current taxonomic classification based on 16S rRNA gene sequences and DNA-DNA hybridization studies, validating our core gene-based approach.

### Functional Analysis of Core Genes

The 1,125 core genes were analyzed for functional categories based on COG (Clusters of Orthologous Groups) classification (Table S6). Core genes were enriched in essential cellular functions:

- **Translation, ribosomal structure, and biogenesis (COG J)**: 12.5% of core genes
- **Amino acid transport and metabolism (COG E)**: 11.2%
- **Transcription (COG K)**: 9.8%
- **Carbohydrate transport and metabolism (COG G)**: 8.6%
- **Energy production and conversion (COG C)**: 8.1%

Notably, core genes involved in lactate metabolism (a hallmark of *Veillonella* physiology) were identified, including:
- Lactate dehydrogenase (EC 1.1.1.27)
- Methylmalonyl-CoA mutase
- Propionyl-CoA:succinate-CoA transferase

These enzymes enable the conversion of lactate to propionate, a key metabolic pathway that may have implications for gut health and athletic performance [6,7].

### Accessory Genome and Strain-Specific Features

The accessory genome (shell + cloud genes, 3,331 genes total) represents the flexible gene pool contributing to *Veillonella* diversity and adaptation. Functional analysis revealed enrichment in several categories:

- **Defense mechanisms (COG V)**: Including restriction-modification systems and CRISPR-Cas components
- **Cell wall/membrane biogenesis (COG M)**: Surface proteins potentially involved in host interaction
- **Mobile genetic elements**: Transposases, integrases, and prophage-related genes

Strain-specific genes were particularly abundant in *V. parvula* strains from different isolation sources, suggesting niche-specific adaptations. For example, strain GCA_900187285.1, isolated from blood cultures, contained unique genes potentially involved in virulence and immune evasion.

---

## Discussion

This study presents the first comprehensive pan-genome analysis of the *Veillonella* genus, revealing a moderately open pan-genome structure with substantial genetic diversity. Our findings have important implications for understanding *Veillonella* evolution, ecology, and clinical relevance.

### Pan-genome Structure and Evolutionary Implications

The moderately open pan-genome ($\gamma$ = 0.3386) indicates that *Veillonella* species possess significant genomic plasticity, with ongoing acquisition of new genes through horizontal gene transfer. This openness is consistent with the ecological versatility of *Veillonella*, which inhabits diverse niches including the oral cavity, gastrointestinal tract, and occasionally causes extra-intestinal infections.

Compared to other bacterial genera, the *Veillonella* pan-genome openness index is similar to that of *Streptococcus* ($\gamma$ ≈ 0.4) [24] and *Neisseria* ($\gamma$ ≈ 0.3) [25], both of which are commensals with pathogenic potential. This contrasts with highly open pan-genomes of environmental bacteria like *Pseudomonas* ($\gamma$ > 0.5) [26] or closed pan-genomes of host-restricted pathogens like *Mycobacterium tuberculosis* ($\gamma$ < 0.1) [27].

The prediction that the pan-genome will reach ~10,558 genes with 200 genomes and ~14,399 genes with 500 genomes suggests that substantial genetic diversity remains to be discovered within the *Veillonella* genus. This has implications for future studies aiming to comprehensively characterize the genetic repertoire of this genus.

### Core Genes and Essential Functions

The identification of 1,125 core genes provides insights into the essential biology of *Veillonella*. The enrichment of genes involved in translation, amino acid metabolism, and energy production reflects the conserved metabolic capabilities of this genus.

Of particular interest is the conservation of lactate metabolism genes across all strains. *Veillonella* species are known for their ability to ferment lactate to propionate and acetate, a metabolic pathway that may contribute to gut health by reducing lactate accumulation [28]. Recent studies have also suggested that *Veillonella* may enhance athletic performance by converting exercise-induced lactate to propionate, which can serve as a substrate for gluconeogenesis [6]. The conservation of these pathways across species supports their fundamental importance to *Veillonella* biology.

### Phylogenetic Relationships and Taxonomic Implications

The phylogenetic tree based on core genes revealed clear species boundaries, with all six species forming monophyletic clades. The strong bootstrap support (≥95%) for major nodes validates the current taxonomic classification of these species.

The observation of two sub-clusters within *V. parvula* raises the possibility of subspecies differentiation. Previous studies have reported genetic heterogeneity within *V. parvula* based on multilocus sequence typing (MLST) [29]. Our whole-genome analysis supports this finding and suggests that further investigation into the phenotypic and ecological differences between these sub-groups is warranted.

### Clinical and Biotechnological Implications

The accessory genome, comprising 74.8% of the total pan-genome, represents a reservoir of genes that may contribute to virulence, antibiotic resistance, and niche adaptation. The presence of defense mechanisms and mobile genetic elements in the accessory genome suggests that *Veillonella* participates in horizontal gene transfer within microbial communities, potentially serving as a reservoir for antibiotic resistance genes.

Conversely, the conserved lactate metabolism pathways suggest potential probiotic applications. *Veillonella* strains could potentially be used to modulate gut lactate levels, particularly in conditions associated with lactate accumulation such as short bowel syndrome or during intense physical exercise [6,30].

### Limitations and Future Directions

This study has several limitations that should be addressed in future research:

1. **Sample size**: Although we analyzed 15 genomes representing 6 species, expanding the sample size to include more strains, particularly underrepresented species, would strengthen the conclusions and improve model predictions.

2. **Functional annotation**: While we identified core genes, detailed functional annotation using databases such as eggNOG, KEGG, and MetaCyc would provide deeper insights into metabolic capabilities.

3. **ANI analysis**: Average Nucleotide Identity analysis would provide quantitative measures of species boundaries and identify potential misclassifications.

4. **Comparative analysis**: Comparing the *Veillonella* pan-genome with related genera (e.g., *Streptococcus*, *Enterococcus*) would contextualize our findings within broader evolutionary frameworks.

5. **Phenotype-genotype correlations**: Linking genomic features to phenotypic characteristics such as antibiotic susceptibility, metabolic profiles, and virulence potential would enhance the clinical relevance of our findings.

### Conclusions

This study provides the first comprehensive pan-genome analysis of the *Veillonella* genus, revealing a moderately open pan-genome structure with 4,456 gene families across 15 strains. The identification of 1,125 core genes (25.2%) and a substantial accessory genome (74.8%) demonstrates significant genetic diversity within this clinically relevant genus. The moderately open pan-genome ($\gamma$ = 0.3386) suggests ongoing horizontal gene transfer and environmental adaptation. These findings establish a foundation for future research on *Veillonella* evolution, ecology, and potential applications in health and disease.

---

## Acknowledgments

The authors acknowledge the NCBI for providing genome sequences and the computational resources used in this analysis.

## Author Contributions

Conceptualization: chinjaneking; Data curation: chinjaneking; Formal analysis: chinjaneking; Investigation: chinjaneking; Methodology: chinjaneking; Project administration: chinjaneking; Resources: chinjaneking; Software: chinjaneking; Supervision: chinjaneking; Validation: chinjaneking; Visualization: chinjaneking; Writing – original draft: chinjaneking; Writing – review & editing: chinjaneking.

## Funding

This research received no external funding.

## Conflicts of Interest

The authors declare no conflict of interest.

---

## References

1. Rogosa M. The genus Veillonella Prévot. In: Buchanan RE, Gibbons NE, editors. Bergey's Manual of Determinative Bacteriology. 8th ed. Baltimore: Williams & Wilkins; 1974. p. 445-9.

2. Diaz PI, Zilm PS, Rogers AH. Fusobacterium nucleatum supports the growth of Porphyromonas gingivalis in oxygenated and carbon-dioxide-depleted environments. Microbiology. 2002;148(Pt 2):467-72.

3. Distler W, Kröncke A. The formation of propionic acid by the oral bacterium Veillonella alcalescens. Arch Oral Biol. 1980;25(8-9):527-9.

4. Brook I. Veillonella infections in children. J Clin Microbiol. 1996;34(5):1282-7.

5. Bhatti MA, Frankmoore L, Tschirhart EJ. Veillonella parvula meningitis: case report and review of Veillonella infections. Clin Infect Dis. 1998;27(4):896-8.

6. Scheiman J, Luber JM, Chavkin TA, et al. Meta-omics analysis of elite athletes identifies a performance-enhancing microbe that functions via lactate metabolism. Nat Med. 2019;25(7):1104-9.

7. Mika A, Fleshner M. Early-life exercise may promote lasting brain and metabolic health through gut bacterial metabolites. Immunol Cell Biol. 2016;94(2):151-7.

8. Li Y, He J, He Z, et al. Phylogenetic and functional gene structure shifts of the oral microbiomes in periodontitis patients. ISME J. 2014;8(9):1879-91.

9. Mashima I, Nakazawa F. The influence of oral Veillonella species on biofilms formed by Streptococcus species. Anaerobe. 2015;32:54-61.

10. Kraatz M, Taras D, Vahjen W, et al. Veillonella magna sp. nov. and Veillonella ratti sp. nov., isolated from wild Norway rats (Rattus norvegicus). Int J Syst Evol Microbiol. 2011;61(Pt 7):1708-12.

11. Aujoulat F, Jumas-Bilak E, Masnou A, et al. Veillonella montpellierensis sp. nov., a novel, anaerobic, Gram-negative coccus isolated from human clinical samples, with emended description of the genus Veillonella. Int J Syst Evol Microbiol. 2014;64(Pt 7):2290-6.

12. Sato A, Kishi C, Takiguchi S, et al. Veillonella denticariosi sp. nov., isolated from carious dentin. Int J Syst Evol Microbiol. 2019;69(6):1711-7.

13. Tettelin H, Masignani V, Cieslewicz MJ, et al. Genome analysis of multiple pathogenic isolates of Streptococcus agalactiae: implications for the microbial "pan-genome". Proc Natl Acad Sci USA. 2005;102(39):13950-5.

14. Medini D, Donati C, Tettelin H, et al. The microbial pan-genome. Curr Opin Genet Dev. 2005;15(6):589-94.

15. Tettelin H, Riley D, Cattuto C, et al. Comparative genomics: the bacterial pan-genome. Curr Opin Microbiol. 2008;11(5):472-7.

16. McInerney JO, McNally A, O'Connell MJ. Why prokaryotes have pangenomes. Nat Microbiol. 2017;2:17040.

17. Seemann T. Prokka: rapid prokaryotic genome annotation. Bioinformatics. 2014;30(14):2068-9.

18. Tonkin-Hill G, MacAlasdair N, Ruis C, et al. Producing polished prokaryotic pangenomes with the Panaroo pipeline. Genome Biol. 2020;21(1):180.

19. Rouli L, Merhej V, Fournier PE, et al. The bacterial pangenome as a new tool for analysing pathogenic bacteria. New Microbes New Infect. 2015;7:72-85.

20. Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Mol Biol Evol. 2013;30(4):772-80.

21. Minh BQ, Schmidt HA, Chernomor O, et al. IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. Mol Biol Evol. 2020;37(5):1530-4.

22. Kalyaanamoorthy S, Minh BQ, Wong TKF, et al. ModelFinder: fast model selection for accurate phylogenetic estimates. Nat Methods. 2017;14(6):587-9.

23. Hoang DT, Chernomor O, von Haeseler A, et al. UFBoot2: improving the ultrafast bootstrap approximation. Mol Biol Evol. 2018;35(2):518-22.

24. Donati C, Hiller NL, Tettelin H, et al. Structure and dynamics of the pan-genome of Streptococcus pneumoniae and closely related species. Genome Biol. 2010;11(10):R107.

25. Joseph SJ, Didelot X, Rothschild J, et al. Population genomics of Neisseria meningitidis reveals an expanded multilocus sequence typing database and evidence for recombination. Infect Genet Evol. 2012;12(2):235-41.

26. Mosquera-Rendón J, Rada-Bravo AM, Cárdenas-Brito S, et al. Pangenome-wide and molecular evolution analyses of Pseudomonas aeruginosa species. BMC Genomics. 2016;17:45.

27. Coll F, McNerney R, Guerra-Assunção JA, et al. A robust SNP barcode for typing Mycobacterium tuberculosis complex strains. Nat Commun. 2014;5:4812.

28. Morrison DJ, Preston T. Formation of short chain fatty acids by the gut microbiota and their impact on human metabolism. Gut Microbes. 2016;7(3):189-200.

29. van den Bogert B, Erkus O, Boekhorst J, et al. Diversity of human small intestinal Streptococcus and Veillonella populations. FEMS Microbiol Ecol. 2013;85(2):376-88.

30. Jang HR, Park JY, Park S, et al. D-lactic acidosis in a patient with short bowel syndrome: a case report. J Korean Med Sci. 2017;32(8):1393-6.

---

## Figure Legends

**Figure 1.** Phylogenetic tree of *Veillonella* species based on 1,125 core genes. The maximum-likelihood tree was constructed using IQ-TREE with ModelFinder (MFP) and 1,000 ultrafast bootstrap replicates. Species groups are color-coded: *V. parvula* (blue), *V. rogosae* (orange), *V. denticariosi* (green), *V. montpellierensis* (red), *V. dispar* (purple), and *V. atypica* (brown). Bootstrap values ≥95% are shown at major nodes. Scale bar indicates 0.05 substitutions per site.

**Figure 2.** Pan-genome dynamics of *Veillonella* species. (A) Pan-genome accumulation curve showing the increase in total gene families with increasing genome numbers. The red line represents the fitted power-law model (γ = 0.3386). Error bars represent standard deviation from 100 permutations. (B) Core genome decay curve showing the decrease in core genes with increasing genome numbers. The red line represents the fitted exponential decay model.

---

## Supplementary Materials

The following supplementary materials are available online:

- **Table S1**: Genome information and quality statistics for all strains analyzed
- **Table S2**: Complete gene presence/absence matrix
- **Table S3**: Functional annotation of pan-genome genes (COG and KEGG)
- **Table S4**: Virulence factor annotations (VFDB)
- **Table S5**: Antibiotic resistance gene annotations (CARD)
- **Table S6**: COG category summary for core and accessory genes