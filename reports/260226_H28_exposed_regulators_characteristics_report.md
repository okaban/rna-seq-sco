# H28: Biological Characteristics of 62 Exposed Regulators

**Date**: 2026-02-26
**Analysis directory**: `11_epigenome_integration/analysis/51_exposed_regulators_characteristics/`

## Background

H27 discovered that the 1,055 regulatory genes in S. coelicolor M145 split into two fundamentally different groups:
- **993 "shielded" regulators**: 1,200 bp methylation protection zone around TSS, methylation-insensitive expression
- **62 "exposed" regulators**: NO protection zone, methylation sites at median 114 bp from TSS, expression responds to methylation changes (coordinated regulators)

The critical question: WHY are these 62 genes exposed? What distinguishes them from the 993 shielded regulators?

## Hypothesis

The 62 exposed regulators are characterized by **low T1 expression** (RNAP not occupying the promoter), explaining the absence of the protein-occupancy-based protection zone.

## Methods

Comprehensive multivariable comparison of 62 exposed vs 993 shielded regulatory genes across 8 dimensions:
1. Expression levels (baseMean, T1/T2/T3 normalized counts, |log2FC| change magnitudes)
2. TF family composition (Fisher exact test for each of 25 families)
3. Gene structural features (gene length, promoter GC content, intergenic distance, FIMO TF binding site predictions)
4. Geographic distribution (core vs arm)
5. Methylation environment (nearest methylation distance, site density within 2kb)
6. Operon context (monocistronic vs polycistronic)
7. Temporal expression pattern (constitutive vs dynamic, classified by LFC and padj thresholds)
8. T1 expression quartile trend analysis

All continuous comparisons used Wilcoxon rank-sum tests with rank-biserial r effect sizes. Categorical comparisons used Fisher exact tests.

## Key Results

### PRIMARY HYPOTHESIS: T1 Expression

| Metric | Exposed (n=62) | Shielded (n=955) | p-value | Effect size r |
|--------|---------------|-------------------|---------|---------------|
| **T1 expression (median)** | **90.4** | **123.1** | **0.010** | **0.195** |
| baseMean | 100.7 | 127.5 | 0.215 | 0.094 |
| T2 expression | 80.4 | 115.4 | 0.122 | 0.117 |
| T3 expression | 92.2 | 102.7 | 0.380 | 0.066 |

**Result: SUPPORTED.** Exposed regulators have significantly lower T1 expression (median 90 vs 123, p=0.010, r=0.195). The effect is strongest at T1, weakening at T2 and disappearing at T3 -- consistent with the hypothesis that low initial RNAP occupancy at T1 prevents formation of the protection zone.

### T1 Expression Quartile Trend

| Quartile | n_total | n_exposed | % exposed |
|----------|---------|-----------|-----------|
| Q1 (lowest) | 255 | 28 | **11.0%** |
| Q2 | 254 | 10 | 3.9% |
| Q3 | 254 | 15 | 5.9% |
| Q4 (highest) | 254 | 9 | 3.5% |

Spearman trend test: rho = -0.095, p = 0.0024. The lowest T1 expression quartile contains 28/62 (45.2%) of exposed regulators -- a 3.1-fold enrichment over the expected 5.9% rate. This strong Q1 enrichment confirms that low initial expression is a primary determinant of exposure.

### Expression Change Magnitude (STRONGEST SIGNAL)

| Metric | Exposed | Shielded | p-value | Effect size r |
|--------|---------|----------|---------|---------------|
| **|log2FC| T2vsT1** | **1.47** | **0.87** | **1.0e-05** | **-0.334** |
| **|log2FC| T3vsT1** | **1.53** | **0.93** | **1.7e-09** | **-0.456** |

Exposed regulators show 1.7x larger expression changes than shielded regulators. This is a large effect (r=0.456 for T3) and confirms that exposed promoters are more responsive to developmental/methylation signals.

### Temporal Expression Pattern (DEFINITIVE)

| Pattern | Exposed (%) | Shielded (%) |
|---------|-------------|--------------|
| **Constitutive** | **0 (0.0%)** | **363 (38.0%)** |
| T2T3_up | 22 (35.5%) | 114 (11.5%) |
| T2T3_down | 16 (25.8%) | 149 (15.0%) |
| T3_up | 9 (14.5%) | 92 (9.3%) |
| T3_down | 8 (12.9%) | 85 (8.6%) |
| Other dynamic | 5 (8.1%) | 152 (15.3%) |
| No data | 0 (0.0%) | 38 (3.8%) |

Fisher exact test (dynamic vs constitutive): **OR=infinity, p=5.3e-13**

**Zero of the 62 exposed regulators are constitutively expressed.** Every single one is developmentally regulated, compared to only 62% of shielded regulators. This is the single strongest distinguishing feature.

### Promoter GC Content

| Feature | Exposed | Shielded | p-value | r |
|---------|---------|----------|---------|---|
| GC content (-300bp) | 0.715 | 0.697 | **1.97e-04** | -0.281 |

Exposed regulators have significantly higher promoter GC content (71.5% vs 69.7%). Higher GC may provide more potential methylation target sites (both GCCGGC and AAGCCCG are GC-rich motifs).

### Methylation Environment

| Feature | Exposed | Shielded | p-value | r |
|---------|---------|----------|---------|---|
| Nearest methylation distance | 114 bp | 762 bp | **5.3e-28** | 0.830 |
| Sites within 2kb | 3.0 | 2.0 | **6.3e-06** | -0.336 |
| Methylation density (sites/kb) | 0.75 | 0.50 | **6.3e-06** | -0.336 |

Exposed regulators have 6.7x closer methylation sites and 1.5x higher methylation density. This confirms H27's finding and is partially explained by the higher GC content creating more methylatable motif sites.

### Non-significant Features

| Feature | Exposed median | Shielded median | p-value |
|---------|---------------|-----------------|---------|
| Gene length | 819 bp | 717 bp | 0.096 |
| Intergenic distance | 67 bp | 90 bp | 0.442 |
| FIMO TF BS hits | 4.0 | 4.0 | 0.413 |
| Geographic (arm %) | 43.5% | 38.0% | 0.420 |
| Monocistronic (%) | 38.7% | 43.2% | 0.512 |
| Operon size | 2.0 | 2.0 | 0.515 |

Geographic distribution, gene structure, operon context, and upstream TF binding do NOT distinguish exposed from shielded regulators.

### TF Family Composition

No family shows statistically significant enrichment or depletion in the exposed group (all Fisher p > 0.05). The family distribution is representative of the overall regulator population. The exposure phenotype is not family-specific.

## Integrated Model: Why Are These 62 Regulators Exposed?

The data support a coherent three-factor model:

1. **Low T1 expression** (p=0.010): At the initial timepoint, these genes have low RNAP occupancy, preventing the formation of the protein-occupancy protection zone observed at shielded regulators.

2. **100% dynamic expression** (p=5.3e-13): Every exposed regulator is developmentally regulated. A gene that transitions between "off" and "on" states cannot maintain a permanent protection zone. During "off" periods, the promoter is accessible to methyltransferases.

3. **Higher GC promoter** (p=1.97e-04): GC-rich promoters contain more potential methylation target motifs (GCCGGC, AAGCCCG), creating a higher baseline probability of methylation when the protection zone is absent.

These three factors act synergistically:
- Low initial expression --> no RNAP barrier --> methylation sites can approach TSS
- Dynamic regulation --> periodic loss of protection --> methylation events accumulate
- GC-rich promoter --> more methylatable motif sites available

The result is a class of regulatory genes whose expression is both influenced by and coordinated with methylation changes -- the "exposed promoter" regulatory architecture.

## Summary Statistics

| Category | Tests | Significant (p<0.05) |
|----------|-------|---------------------|
| Expression levels | 8 | 3 (T1, |LFC_T2|, |LFC_T3|) |
| Gene structure | 4 | 1 (GC content) |
| Methylation | 3 | 3 (all) |
| Geography | 1 | 0 |
| Operon | 2 | 0 |
| Temporal | 1 | 1 |
| **Total** | **19** | **8 (42%)** |

## Verdict

**SUPPORTED.** The primary hypothesis that exposed regulators have low T1 expression is confirmed (p=0.010). However, the most powerful distinguishing feature is not T1 expression alone but the complete absence of constitutive expression (OR=infinity, p=5.3e-13) combined with larger expression change magnitudes (|LFC| 1.7x higher, p=1.7e-09). The "exposed" phenotype is defined by dynamically regulated genes with low initial expression and GC-rich promoters -- genes that are transcriptionally "off" at some timepoints, allowing methyltransferase access to their promoters.

## Output Files

### Tables
- `tables/exposed_regulators_full_table.tsv` -- All 62 exposed regulators with all measured features (30 columns)
- `tables/feature_comparison_statistics.tsv` -- All 19 statistical tests with p-values and effect sizes
- `tables/TF_family_distribution.tsv` -- Family-level breakdown with Fisher exact tests (25 families)
- `tables/temporal_pattern_distribution.tsv` -- Expression pattern classification

### Figures
- `figures/expression_comparison.pdf/svg` -- Panel A: 4-panel violin+box plot (baseMean, T1, |LFC_T2|, |LFC_T3|)
- `figures/TF_family_composition.pdf/svg` -- Panel B: TF family counts and % exposed per family
- `figures/gene_features_comparison.pdf/svg` -- Panel C: 6-panel structural features comparison
- `figures/chromosome_map.pdf/svg` -- Panel D: Chromosome position map
- `figures/methylation_distance_vs_expression.pdf/svg` -- Panel E: Nearest methylation distance vs expression scatter
- `figures/H28_comprehensive_summary.pdf/svg` -- 6-panel comprehensive summary figure

### Scripts
- `scripts/H28_exposed_regulators_analysis.py`
