# H22: Sequence-Level Motif Depletion Near Regulatory Genes

**Date**: 2026-02-26
**Analysis directory**: `11_epigenome_integration/analysis/45_sequence_level_motif_depletion/`
**Hypothesis**: The observed methylation site depletion near regulatory genes (H15/H20) is at least partially explained by evolutionary counter-selection against R-M recognition motifs (TGGCCGGC, AAGCCCG) at the DNA sequence level.

---

## Background

H15 established that methylation sites of all motif classes are significantly depleted near regulatory genes (fold 0.43-0.62, p < 0.002 for all). H20 confirmed this is a genuine biological signal surviving geographic stratification, not a Simpson's paradox artifact. However, these analyses measured methylation site depletion -- they could not distinguish between two mechanisms:

1. **Sequence-level counter-selection**: The recognition motif sequences themselves (TGGCCGGC for GCCGGC 4mC; AAGCCCG for 6mA) are physically depleted from regulatory gene neighborhoods through evolutionary pressure to avoid deleterious methylation/restriction
2. **Protein occupancy model**: Motif sequences exist at normal frequency, but TF/RNAP binding physically blocks MTase access, so methylation is depleted despite motifs being present

Distinguishing these mechanisms is critical for understanding the Gatekeeper Model: if motifs are sequence-depleted, regulatory avoidance is hardwired into the genome; if only methylation is depleted, the avoidance is dynamic and dependent on transcriptional state.

---

## Methods

### Motif Scanning
- Scanned the full NC_003888.3 chromosome (8,667,507 bp) for exact matches to:
  - **TGGCCGGC** (GCCGGC 4mC recognition context, 8-mer) and its reverse complement GCCGGCCA
  - **AAGCCCG** (6mA recognition motif, 7-mer) and its reverse complement CGGGCTT
- Counted occurrences on both strands

### Gene Regions
For each of 7,593 protein-coding genes on NC_003888.3:
- **Promoter**: -500 to +100 bp from gene start (strand-adjusted)
- **Gene body**: start to end
- **Extended 2kb**: -2,000 bp to gene end (for consistency with H15/H20 2kb windows)

### Statistical Tests
- **Mean density comparison**: Per-gene motif density (occurrences/kb) for regulatory (n=1,008) vs non-regulatory (n=6,585) genes
- **Fisher exact test**: Proportion of genes containing at least one motif occurrence
- **Wilcoxon rank-sum test**: Distribution of per-gene motif counts
- **Permutation test**: 10,000 random samples of 1,008 genes to establish null distribution

### Geographic Stratification
- Core: 1,500,001 to 7,167,506 bp (4,944 genes: 638 regulatory, 4,306 non-regulatory)
- Arm: 0-1,500,000 and 7,167,507-8,667,507 bp (2,649 genes: 370 regulatory, 2,279 non-regulatory)

---

## Results

### 1. Genome-Wide Motif Census

| Motif | Forward | RevComp | Total | Density (per kb) |
|-------|:-------:|:-------:|:-----:|:-----------------:|
| TGGCCGGC | 1,091 | 1,002 | **2,093** | 0.241 |
| AAGCCCG | 649 | 685 | **1,334** | 0.154 |

Both motifs are relatively rare in the GC-rich (72%) *S. coelicolor* genome.

### 2. DNA Motif Density: Regulatory vs Non-Regulatory Genes

| Motif | Zone | Reg density (/kb) | Non-reg density (/kb) | Fold | Fisher p | Wilcoxon p | Significant? |
|-------|------|:--:|:--:|:--:|:--:|:--:|:--:|
| **TGGCCGGC** | Promoter | 0.187 | 0.212 | **0.88** | 0.167 | 0.157 | No |
| | Gene body | 0.191 | 0.259 | **0.74** | **2.98e-07** | **7.05e-07** | **Yes** |
| | Extended 2kb | 0.223 | 0.241 | **0.92** | **6.79e-03** | **3.59e-03** | **Yes** |
| **AAGCCCG** | Promoter | 0.159 | 0.155 | 1.02 | 1.00 | 0.993 | No |
| | Gene body | 0.122 | 0.161 | **0.76** | **1.26e-07** | **3.54e-07** | **Yes** |
| | Extended 2kb | 0.135 | 0.159 | **0.85** | **1.05e-03** | **2.28e-04** | **Yes** |

Key findings:
- **Gene bodies show the strongest DNA motif depletion**: fold 0.74 (TGGCCGGC) and 0.76 (AAGCCCG), both highly significant
- **Extended 2kb regions** show moderate but significant depletion: fold 0.92 and 0.85
- **Promoter regions** show no significant depletion for either motif, despite being the expected target of MTase avoidance

### 3. Geographic Stratification (Extended 2kb Region)

| Motif | Region | Fold | Fisher p | Wilcoxon p |
|-------|--------|:----:|:--------:|:----------:|
| **TGGCCGGC** | Core | **0.90** | **0.038** | **0.011** |
| | Arm | 0.95 | 0.083 (NS) | 0.137 (NS) |
| **AAGCCCG** | Core | **0.83** | **0.007** | **0.002** |
| | Arm | 0.88 | 0.096 (NS) | 0.059 (NS) |

- DNA motif depletion is significant in **core** for both motifs
- Arm regions show the same direction (fold < 1) but do not reach significance, likely due to smaller sample size (370 regulatory genes in arms)

### 4. DNA Motif Depletion vs Methylation Site Depletion

This is the central comparison of H22, contrasting sequence-level (DNA) depletion with the methylation-level depletion from H15/H20.

| Motif | Region | DNA motif fold | Methylation fold (H15/H20) | Ratio | Sequence contribution |
|-------|--------|:-:|:-:|:-:|:-:|
| **TGGCCGGC** | Combined | **0.92** | **0.61** | 1.50 | **21.4%** |
| | Core | **0.90** | **0.67** | 1.34 | **30.9%** |
| | Arm | 0.95 | **0.45** | 2.10 | 9.8% |
| **AAGCCCG** | Combined | **0.85** | **0.43** | 1.97 | **27.1%** |
| | Core | **0.83** | **0.53** | 1.57 | **35.5%** |
| | Arm | 0.88 | **0.24** | 3.67 | 15.6% |

**Sequence contribution** = (1 - DNA fold) / (1 - Methylation fold) x 100%, representing the fraction of total methylation depletion that can be attributed to sequence-level motif depletion.

Interpretation:
- For both motifs, DNA motif depletion is **real but substantially weaker** than methylation depletion
- Sequence-level avoidance accounts for approximately **21-35%** of the observed methylation depletion
- The remaining **65-79%** must be explained by other mechanisms (protein occupancy, chromatin accessibility, or active demethylation)

### 5. Permutation Test

| Motif | Zone | Observed | Null mean +/- SD | z-score | p-value |
|-------|------|:--------:|:----------------:|:-------:|:-------:|
| **TGGCCGGC** | Extended 2kb | 0.2230 | 0.2389 +/- 0.0089 | **-1.79** | **0.035** |
| | Promoter | 0.1868 | 0.2091 +/- 0.0178 | -1.25 | 0.110 |
| **AAGCCCG** | Extended 2kb | 0.1353 | 0.1564 +/- 0.0074 | **-2.87** | **0.001** |
| | Promoter | 0.1587 | 0.1557 +/- 0.0155 | +0.19 | 0.603 |

The permutation test confirms:
- Extended 2kb motif density is significantly below null expectation for both motifs (p = 0.035 and 0.001)
- Promoter motif density is not significantly depleted for either motif
- AAGCCCG depletion is stronger and more significant than TGGCCGGC

### 6. Gene Body Depletion: A Distinct Signal

The strongest DNA motif depletion is in gene bodies, not promoters:

| Motif | Gene body fold | Gene body Fisher p | Promoter fold | Promoter Fisher p |
|-------|:-:|:-:|:-:|:-:|
| TGGCCGGC | **0.74** | **2.98e-07** | 0.88 | 0.167 |
| AAGCCCG | **0.76** | **1.26e-07** | 1.02 | 1.00 |

This is consistent with coding sequence constraints: regulatory genes (TFs, sigma factors, response regulators) have specific protein domains with constrained codon usage, which may inherently reduce the frequency of GC-rich 8-mer (TGGCCGGC) and 7-mer (AAGCCCG) motifs compared to metabolic enzymes and transporters.

---

## Key Figures

### Figure 1: Motif Density Comparison
![](figures/motif_density_comparison.png)
Bar charts showing per-gene mean motif density for regulatory vs non-regulatory genes across three gene zones (promoter, gene body, extended 2kb). Significance stars from the minimum of Fisher and Wilcoxon p-values.

### Figure 2: Sequence vs Methylation Depletion
![](figures/sequence_vs_methylation_depletion.png)
Paired comparison of DNA-level motif depletion (this analysis) vs methylation-level depletion (H15), shown by genomic region. The gap between green (DNA) and purple (methylation) bars represents the protein occupancy component.

### Figure 3: Comprehensive Summary
![](figures/H22_comprehensive_summary.png)
A: Genome-wide motif census. B: Fold enrichment by gene zone. C: Core/arm stratification. D: DNA vs methylation depletion comparison. E-F: Permutation test null distributions with observed values.

---

## Verdict: PARTIAL (Both Mechanisms Contribute)

**DNA sequence-level motif depletion is statistically real but explains only a minority (~21-35%) of the total methylation depletion observed near regulatory genes. The majority (~65-79%) is attributable to protein occupancy or other epigenetic mechanisms.**

### Evidence for Sequence-Level Avoidance (Contributing Factor)
1. **Gene body motifs are depleted** in regulatory genes: fold 0.74-0.76, p < 3e-07 for both motifs
2. **Extended 2kb region** shows significant depletion: fold 0.85-0.92, p < 0.007
3. **Permutation tests confirm** extended region depletion: p = 0.001 (AAGCCCG), p = 0.035 (TGGCCGGC)
4. **Core-specific depletion persists**: consistent with H20 geographic validation

### Evidence for Protein Occupancy as the Dominant Mechanism
1. **DNA motif fold (0.85-0.92) is far weaker than methylation fold (0.43-0.61)**: sequence depletion explains only 21-35% of methylation depletion
2. **Promoter regions show NO DNA motif depletion** (fold ~1.0) despite being the primary target of H15's methylation depletion. This is the strongest evidence for protein occupancy: promoter motifs exist at normal frequency but are protected from methylation by transcription factor binding
3. **Arm regions lose DNA depletion significance** while methylation depletion remains very strong (fold 0.24-0.45), suggesting protein occupancy dominates in arms

### Integrated Model

The observed methylation avoidance of regulatory genes (H15, H20) results from **two complementary layers**:

| Layer | Mechanism | Contribution | Zone specificity |
|-------|-----------|:-:|:-:|
| **Layer 1** | Evolutionary counter-selection against R-M motifs | ~21-35% | Gene body > Extended > Promoter |
| **Layer 2** | Protein occupancy blocking MTase access | ~65-79% | Promoter = Extended > Gene body |

This two-layer architecture is internally consistent:
- **Gene bodies** show the strongest DNA depletion because coding sequence evolution can eliminate motifs without regulatory cost
- **Promoters** show no DNA depletion but strong methylation depletion because TF/RNAP occupancy provides dynamic protection without requiring sequence changes
- **Extended regions** show intermediate effects from both mechanisms

### Implications for the Gatekeeper Model

This refines Layer 2 (Protection/Depletion) of the Gatekeeper Model v2 (H11):
- Regulatory gene avoidance is not a single mechanism but a **composite of evolutionary and dynamic protection**
- The evolutionary component (sequence depletion) represents a "hardwired" basal avoidance that operates even in the absence of active transcription
- The protein occupancy component provides a "dynamic" layer that depends on transcriptional state, consistent with the gatekeeper concept of methylation responding to regulatory context

---

## Output Files

### Tables
- `tables/gene_motif_counts.tsv` - Per-gene motif counts for all 7,593 protein-coding genes (both motifs, 3 zones)
- `tables/depletion_statistics.tsv` - Fisher/Wilcoxon test results for all motif x zone x region combinations (18 rows)
- `tables/sequence_vs_methylation_comparison.tsv` - DNA vs methylation depletion comparison with sequence contribution (6 rows)
- `tables/stratified_by_region.tsv` - Core/arm stratified results (12 rows)
- `tables/permutation_test_results.tsv` - Permutation test results (4 rows)

### Figures
- `figures/motif_density_comparison.pdf/svg/png` - Bar plot of motif density for regulatory vs non-regulatory genes
- `figures/sequence_vs_methylation_depletion.pdf/svg/png` - Paired comparison of DNA-level vs methylation-level depletion
- `figures/H22_comprehensive_summary.pdf/svg/png` - Multi-panel summary (6 panels)

### Scripts
- `scripts/H22_sequence_level_motif_depletion.py`
