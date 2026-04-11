# Methods

## Bacterial strain and growth conditions

*Streptomyces coelicolor* A3(2) strain M145 was cultured under standard conditions. Three developmental stages were sampled: exponential growth (T1), transition phase (T2), and stationary phase (T3), each with three biological replicates (nine samples total).

## Nanopore sequencing and DNA methylation detection

Genomic DNA was sequenced using Oxford Nanopore Technologies (ONT) long-read sequencing. Reads were aligned to the *S. coelicolor* A3(2) M145 reference genome (NCBI RefSeq GCF_000203835.1) using minimap2 with the `-y` flag to preserve base modification tags (MM, ML, MN). Sequencing coverage ranged from 14 to 52x per sample (eight samples at 27–52x; one sample at 14.3x).

DNA methylation was called using modkit (ONT) in pileup mode, generating per-site methylation calls from aligned BAM files. Two modification types were detected with high confidence: N4-methylcytosine (4mC; 3,294 total sites) and N6-methyladenine (6mA; 2,734 total sites). 5-methylcytosine (5mC) was not detected at high confidence. High-confidence methylation sites were defined as sites passing modkit's default quality thresholds. Motif assignment was performed by extracting the sequence context of each methylation site and identifying enriched motifs: GCCGGC for 4mC and AAGCCCG for 6mA.

For samples with initially low coverage (samples 3-2 and 3-4), additional ONT sequencing data were generated and merged with existing alignments using a standardized pipeline to improve methylation calling confidence.

## Methylation site classification by timepoint

Methylation sites were classified by the timepoints at which they were detected. "T1 sites" were those called as methylated at T1 (with or without methylation at later timepoints); "T2 sites" and "T3 sites" were defined analogously. "Lost" sites were those methylated at T1 but not at T2; "Never" sites were those not methylated at either timepoint. Transition categories (e.g., concordant_loss, discordant_gain_up) were assigned based on the combination of methylation status change and expression change direction.

## RNA-seq library preparation and sequencing

Total RNA was extracted from each of the nine samples (three timepoints x three biological replicates). Paired-end Illumina sequencing (150 bp reads) was performed. Raw reads were quality-controlled and adapter-trimmed using fastp.

## RNA-seq alignment and quantification

Trimmed reads were aligned to the *S. coelicolor* A3(2) M145 reference genome (NCBI RefSeq GCF_000203835.1) using HISAT2 v2.2.1 in prokaryotic mode (no spliced alignment). Overall alignment rate was 98.54%, with 96.55% uniquely mapped reads. Gene-level read counts were quantified using featureCounts.

## Differential expression analysis

Differential expression analysis was performed using DESeq2 with the design formula `~ condition`, where condition corresponds to the three developmental stages (T1, T2, T3). Size factor normalization was applied using the default DESeq2 method. Log2 fold change (LFC) estimates were shrunk using the apeglm method. Three pairwise contrasts were computed: T2 vs T1, T3 vs T1, and T3 vs T2. Genes with fewer than 10 counts in fewer than 3 samples were filtered prior to analysis, yielding 7,646 genes tested. Differentially expressed genes (DEGs) were defined as those with Benjamini-Hochberg adjusted *p* < 0.05 and |LFC| > 1. Mean normalized expression (baseMean) was computed as the average of DESeq2-normalized counts across all nine samples.

## Regulatory gene annotation

Regulatory genes were identified through a two-stage process. First, 873 genes were automatically classified as regulators based on product annotation pattern matching against known regulatory protein families (TFs, TCS components, sigma factors). Second, nine additional genes were manually curated based on established literature: *actII-orf4* (SCO5082), *redD* (SCO5877), *redZ* (SCO5881), *absA1* (SCO3225), *absA2* (SCO3226), *kasO* (SCO6282), *scbR* (SCO6265), *scbR2* (SCO6286), and *atrA* (SCO4677). Regulatory genes were further classified into families: sigma factors (76 genes), TCS sensor kinases (83), TCS response regulators (81), SARP-family activators (7), and other TFs (786). The final set comprised 1,055 regulatory genes with assigned TSS positions.

## Transcription start site assignment

TSS positions were assigned by integrating two data sources. The primary source was the gene start coordinate from the NCBI RefSeq GFF annotation. Where available, experimentally determined TSS positions from dRNA-seq data (Jeong et al. 2016) were used instead (262 genes with experimental TSS data). For genes on the minus strand, distance calculations were adjusted so that negative values indicate upstream positions and positive values indicate positions within the gene body. The final comprehensive TSS table covered 7,996 protein-coding genes.

## TSS-centered methylation density profiles

Methylation density profiles were computed in 200 bp bins spanning ±5 kb from each TSS. For each bin, the number of methylation sites per kb per gene was calculated separately for regulatory genes, non-regulatory genes, and all genes. Gaussian smoothing (sigma = 2 bins, effective 400 bp) was applied to reduce noise. Bootstrap confidence intervals (1,000 resamples) were computed for each bin. The protection zone was defined as the contiguous region where regulatory gene methylation density was depleted relative to flanking regions, with the depletion ratio calculated as the density within the zone divided by the mean density of flanking regions (±3–5 kb).

## Shielded/Exposed classification

For each of the 1,055 regulatory genes, the nearest methylation site distance was computed as the minimum distance between any high-confidence methylation site (across all timepoints and types) and the gene's TSS. ROC analysis was performed using the nearest methylation distance as a predictor of a binary classification variable. The optimal threshold was determined using Youden's *J* statistic (maximizing sensitivity + specificity - 1). AUC with 95% CI was computed using bootstrap resampling. Genes with nearest methylation distance below the optimal threshold (293 bp) were classified as "exposed"; all others were classified as "shielded." The classification was validated by comparison with baseMean (mean expression level) as a negative control predictor. Expression quintile analysis used the Jonckheere-Terpstra (JT) test for trend.

## Sequence determinants of protection

R-M recognition motif density was quantified in the DNA sequence (independent of methylation status) within ±300 bp of each regulatory gene TSS. The AAGCCCG and TGGCCGGC (GCCGGC palindrome) motifs were scanned using exact string matching. Enrichment at exposed versus shielded TF promoters was tested using Fisher's exact test. A logistic regression model incorporating all sequence features (individual motif counts and combined features) was trained to predict exposed/shielded status, with performance evaluated by 5-fold cross-validation AUC. The sequence contribution to the Shielded/Exposed dichotomy was estimated as the ratio of the sequence model's cross-validated AUC to the nearest methylation distance AUC.

## TF binding site methylation analysis

Known TF binding motifs were obtained from CollecTF and CIS-BP databases and scanned across the M145 genome using FIMO (Find Individual Motif Occurrences). Methylation density profiles were computed centered on predicted TF BSs (±2 kb), analogous to the TSS-centered profiles. The enrichment fold at TF BSs was calculated as the observed methylation density divided by the expected density from randomly positioned control sites, matched for GC content and chromosomal position.

## Geographic stratification and Simpson's paradox analysis

The *S. coelicolor* linear chromosome was divided into core (1.5–7.17 Mb) and arm regions (0–1.5 Mb, left arm; 7.17–8.67 Mb, right arm) based on established genomic architecture. For each methylation-expression association test, analyses were performed three ways: unstratified (all genes), core-only, and arm-only. The Cochran-Mantel-Haenszel (CMH) test was used for geographic stratification of categorical methylation-expression associations. For continuous outcomes (LFC comparisons between Lost and Never gene sets), Mann-Whitney *U* tests were performed separately within each geographic stratum. Simpson's paradox was diagnosed when a significant unstratified association became non-significant in both geographic strata.

## Co-expression analysis and module detection

A 62 x 62 Spearman correlation matrix was computed for the 62 exposed regulators using rlog-transformed expression values across all nine samples. Hierarchical clustering was performed using Ward's linkage method with Euclidean distance on the correlation matrix. Modules were identified by cutting the dendrogram at a level that maximized within-module coherence while maintaining biologically interpretable group sizes. Blocs were assigned by merging modules with concordant temporal trajectories: modules with net positive LFC (T3 vs T1) were assigned to the activation bloc, and modules with net negative LFC were assigned to the repression bloc.

Within-type and between-type co-expression medians were compared using Wilcoxon rank-sum tests. Genomic clustering was assessed by computing the mean pairwise chromosomal distance between exposed regulators and comparing it to a null distribution from 10,000 random gene sets of equal size. Operonic pair status was determined from operon predictions in the reference annotation.

## Temporal dynamics analysis

For each of the 62 exposed regulators, expression values were *z*-score normalized across the three timepoints. The phase ratio was defined as the proportion of total absolute expression change (|LFC_T3vsT1|) attributable to the T1-to-T2 transition: phase ratio = |LFC_T2vsT1| / (|LFC_T2vsT1| + |LFC_T3vsT2|). Genes with phase ratio > 0.6 were classified as "early responders." The Mann-Whitney *U* test was used to compare phase ratio distributions between the activation and repression blocs. Module-internal temporal coherence was assessed by computing the median Spearman correlation of temporal trajectories (T1, T2, T3 *z*-scores) among all gene pairs within each module. Methylation-expression timing correlation was assessed by Spearman correlation between the magnitude of methylation change and expression change at each transition (T1-T2 and T2-T3).

## TF family enrichment analysis

TF family annotations were assigned based on protein domain architecture and product description. For each TF family represented among the 62 exposed regulators, Fisher's exact test was used to assess enrichment in the activation versus repression bloc. The TetR family was the only family achieving statistical significance after considering multiple comparisons.

For TCS pair analysis, cognate SK-RR pairs were identified from genomic adjacency and domain complementarity. The asymmetry of exposed/shielded status within pairs was tested using a binomial test (null: 50% probability of both members sharing the same status). The SK/RR composition of exposed partners was tested for bias using a binomial test (null: 50% SK).

## Neighborhood expression analysis

For each exposed TF, the mean absolute LFC (T3 vs T1) of flanking genes (±10 genes) was computed. A null distribution was generated by computing the same statistic for 10,000 randomly sampled gene sets of size 62 from the full genome. The permutation *p*-value was the fraction of random sets with mean neighbor |LFC| greater than or equal to the observed value. Distance decay was assessed by Spearman correlation between chromosomal distance (in bp) from the exposed TF and the neighbor's |LFC|. Transcriptional concordance rate was defined as the fraction of neighbors with the same sign of LFC as the adjacent exposed TF.

## Conservation analysis

Evolutionary conservation of exposed versus shielded regulators was assessed using multiple proxy measures available from the genome annotation: (1) SCO locus tag assignment rate (presence in the original *S. coelicolor* annotation), (2) GC content at third codon positions (GC3), (3) effective number of codons (Nc), and (4) rare codon frequency. Each metric was compared between exposed and shielded groups using Mann-Whitney *U* tests. A composite conservation score was computed as the mean of z-scored individual metrics, and ROC analysis was performed to assess the discriminative power of conservation for predicting exposed status.

## Statistical analysis and software

All analyses were performed in Python 3 using NumPy, pandas, SciPy (v1.11+), statsmodels, and scikit-learn. Differential expression analysis used DESeq2 (v1.38+) in R (v4.3+). Visualization was performed using matplotlib and seaborn. Multiple testing correction used the Benjamini-Hochberg procedure where applicable. *p*-values are reported in scientific notation for *p* < 0.001 and in decimal notation for *p* >= 0.001. Effect sizes (AUC, *rho*, odds ratio, fold change) are reported alongside *p*-values throughout. All code is available at [repository URL].
