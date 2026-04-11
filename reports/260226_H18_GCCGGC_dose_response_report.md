# H18: GCCGGC Methylation Density-Expression Dose-Response Relationship

**Date:** 2026-02-26
**Analysis directory:** `11_epigenome_integration/analysis/41_GCCGGC_dose_response/`
**Hypothesis:** Genes with more GCCGGC 4mC sites within 2 kb show progressively stronger expression suppression, demonstrating a dose-response relationship consistent with methylation acting as a direct transcriptional brake.

## Verdict: PARTIAL

The dose-response is statistically significant at the population level (Jonckheere-Terpstra p = 4.6e-14 for T2vsT1, p = 6.8e-14 for T3vsT1), but the effect is driven almost entirely by the binary contrast of 0 vs >=1 sites. Among genes that already have methylation, additional sites provide negligible further suppression (Spearman rho = -0.037, p = 0.026 for T2vsT1; rho = -0.020, p = 0.215 for T3vsT1). The dose-response curve is not monotonically decreasing: the 4+ group rebounds toward baseline. This is consistent with methylation as a threshold switch rather than a graded rheostat.

---

## 1. Background

H17 discovered that GCCGGC 4mC proximity is associated with significant expression suppression (p = 4.1e-10, r = 0.09). That analysis was binary (proximal vs not proximal). A true dose-response relationship -- where MORE methylation sites near a gene correlate with STRONGER suppression -- would provide much stronger evidence for direct causation and distinguish between a threshold switch vs. a proportional brake model.

## 2. Methods

### 2.1 Site-to-gene mapping
- 1,289 T1 GCCGGC 4mC sites mapped to all protein-coding genes within 2 kb using GFF coordinates
- Sites classified by location: promoter (<=500 bp upstream of gene start on coding strand), gene body, upstream, downstream
- Per-gene site counts computed (0, 1, 2, 3, 4+ groups)

### 2.2 Expression data
- DESeq2 results: T2vsT1 and T3vsT1 (v1, 2026-01-28)
- 7,496 genes with non-NA log2FoldChange values after merge

### 2.3 Statistical tests
- **Spearman correlation**: site_count vs LFC (all genes; genes with >=1 site separately)
- **Jonckheere-Terpstra trend test**: ordered monotonic dose-response across 0/1/2/3/4+ groups
- **Kruskal-Wallis test**: any-group differences in LFC distribution
- **Partial correlation**: controlling for gene length
- **Stratified analysis**: core vs arm genomic regions separately

### 2.4 Comparisons
- AAGCCCG 6mA at T1 analyzed with identical methods as a control motif
- Promoter vs gene body location-specific dose-response

## 3. Results

### 3.1 GCCGGC site distribution

| Sites per gene | N genes | % of genome |
|:-:|:-:|:-:|
| 0 | 4,272 | 52.8% |
| 1 | 2,392 | 29.6% |
| 2 | 1,003 | 12.4% |
| 3 | 299 | 3.7% |
| 4+ | 117 | 1.4% |

3,811 genes (47.2%) have at least one GCCGGC 4mC site within 2 kb. The 5,791 total site-gene mappings break down as: downstream 2,314 (40.0%), upstream 1,752 (30.3%), gene body 1,225 (21.2%), promoter 500 (8.6%).

### 3.2 Dose-response: T2 vs T1

| Dose group | N genes | Median LFC | Mean LFC | % DEG | % Up | % Down |
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| 0 | 3,819 | **+0.120** | +0.432 | 67.2% | 36.0% | 31.3% |
| 1 | 2,294 | **-0.108** | +0.255 | 69.3% | 31.9% | 37.4% |
| 2 | 977 | **-0.331** | +0.107 | 74.2% | 29.6% | 44.6% |
| 3 | 291 | **-0.318** | +0.169 | 71.5% | 28.5% | 43.0% |
| 4+ | 115 | **-0.089** | +0.355 | 61.7% | 32.2% | 29.6% |

**Key pattern**: Clear suppression from 0 -> 1 -> 2 sites (median LFC drops from +0.120 to -0.331). However, the trend **plateaus at 3 sites** and **reverses at 4+**, where median LFC rebounds to -0.089. The fraction of downregulated DEGs also peaks at dose 2 (44.6%) and drops back at 4+ (29.6%).

### 3.3 Dose-response: T3 vs T1

| Dose group | N genes | Median LFC | Mean LFC | % DEG | % Up | % Down |
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| 0 | 3,819 | **+0.335** | +0.675 | 79.8% | 45.8% | 34.0% |
| 1 | 2,294 | **-0.027** | +0.389 | 81.3% | 40.4% | 40.9% |
| 2 | 977 | **-0.294** | +0.182 | 83.0% | 36.3% | 46.7% |
| 3 | 291 | **+0.018** | +0.531 | 84.2% | 42.3% | 41.9% |
| 4+ | 115 | **+0.111** | +0.696 | 89.6% | 46.1% | 43.5% |

The T3 pattern is weaker. The 0 -> 1 -> 2 suppression is present but the 3 and 4+ groups revert to near-zero or positive median LFC. This indicates the dose-response is not sustained during sporulation.

### 3.4 Statistical tests

| Test | T2 vs T1 | T3 vs T1 |
|---|---|---|
| **Spearman (all genes)** | rho = **-0.087**, p = **4.6e-14** | rho = **-0.086**, p = **7.1e-14** |
| **Spearman (>=1 site only)** | rho = **-0.037**, p = **0.026** | rho = **-0.020**, p = **0.215** |
| **Jonckheere-Terpstra** | p = **4.6e-14** | p = **6.8e-14** |
| **Kruskal-Wallis** | H = 61.0, p = **1.8e-12** | H = 65.9, p = **1.7e-13** |

**Critical observation**: The all-genes Spearman is highly significant (rho ~ -0.087) because the 0-site vs >=1-site contrast is large. But restricting to genes with >=1 site, the dose-response within the methylated group is very weak: rho = -0.037 (T2, barely significant) and rho = -0.020 (T3, non-significant). This demonstrates that the statistical significance is overwhelmingly driven by the **presence/absence** binary, not by a graded dose-response.

### 3.5 Location-specific analysis

| Location | Comparison | 0 vs >=1 median diff | MWU p-value |
|---|---|---|---|
| **Promoter** | T2vsT1 | -0.142 | 0.252 |
| **Gene body** | T2vsT1 | -0.313 | **9.2e-08** |
| **Promoter** | T3vsT1 | -0.062 | 0.934 |
| **Gene body** | T3vsT1 | -0.384 | **7.4e-06** |

Contrary to typical eukaryotic models where promoter methylation is most suppressive, here **gene body methylation** shows the significant dose-dependent suppression effect. Promoter-proximal methylation (<=500 bp upstream) shows no significant dose-response (Spearman rho = +0.052 for T2, +0.065 for T3, both non-significant). This may reflect the fact that GCCGGC sites are predominantly gene-body or intergenic, with only 500 mappings (8.6%) falling in the narrow promoter window.

### 3.6 Confound controls

**Gene length**: Partial Spearman correlation controlling for gene length is virtually identical to the raw correlation (rho_partial = -0.087 for T2, -0.093 for T3), confirming gene length is not a confound. The site_count vs gene_length correlation is weak (rho = 0.091).

**Arm vs Core stratification**: Within-region analysis reveals that the dose-response among genes with >=1 site **disappears entirely** when stratifying:

| Region | Comparison | Spearman rho | p-value | n |
|---|---|---|---|---|
| Core | T2vsT1 | +0.010 | 0.58 | 2,924 |
| Arm | T2vsT1 | -0.039 | 0.29 | 753 |
| Core | T3vsT1 | **+0.050** | **0.007** | 2,924 |
| Arm | T3vsT1 | -0.061 | 0.09 | 753 |

Strikingly, in core genes at T3, the correlation is **positive** (more sites = higher expression), opposite to the hypothesized direction. This suggests the overall dose-response signal may partly reflect the arm/core compositional difference (core genes have more T1 GCCGGC sites AND different baseline expression patterns) rather than a direct causal effect of methylation density.

### 3.7 AAGCCCG comparison

| Motif | Comparison | Spearman (all) | p-value | JT p-value |
|---|---|---|---|---|
| **GCCGGC** | T2vsT1 | -0.087 | 4.6e-14 | 4.6e-14 |
| **GCCGGC** | T3vsT1 | -0.086 | 7.1e-14 | 6.8e-14 |
| **AAGCCCG** | T2vsT1 | -0.003 | 0.800 | 0.800 |
| **AAGCCCG** | T3vsT1 | -0.018 | 0.124 | 0.125 |

AAGCCCG shows **no population-level dose-response** (JT and Spearman p > 0.1 for all genes). However, among genes with >=1 AAGCCCG site, there IS a significant negative correlation (rho = -0.120, p = 2.8e-04 for T2; rho = -0.146, p = 1.0e-05 for T3). This is an intriguing reversal of the GCCGGC pattern: GCCGGC has a strong binary signal but weak within-methylated dose-response, while AAGCCCG has no binary signal but a relatively stronger within-methylated dose-response. The AAGCCCG dose-response is particularly striking at high doses (3 sites: median LFC = -0.86 for T2, -1.11 for T3; 4+ sites: -1.05, -1.86) but the sample sizes are very small (n=34 and n=7).

## 4. Interpretation

### 4.1 Threshold switch, not a graded rheostat

The data is most consistent with GCCGGC 4mC acting as a **threshold switch**: the mere presence of methylation near a gene shifts expression downward, but additional sites provide diminishing returns. The dose-response curve shape -- steep drop from 0 to 1-2 sites, plateau at 3, rebound at 4+ -- argues against a simple linear "methylation brake" model.

### 4.2 The 4+ group anomaly

The rebound at 4+ sites (115 genes) may reflect:
1. **Genomic context bias**: Genes in GCCGGC-dense regions may have distinct regulatory properties
2. **Gene function bias**: Highly methylated genes may include essential housekeeping genes that resist suppression
3. **Statistical noise**: n=115 is relatively small

### 4.3 Confound concern: arm/core composition

The strongest caveat is the arm/core stratification result. The within-region dose-response disappears, and core T3 shows a positive correlation. This suggests the overall signal may partially reflect the fact that GCCGGC sites are concentrated in the core genome, where genes tend to be more stably expressed, rather than a direct causal dose-response of methylation on transcription.

### 4.4 AAGCCCG shows a different pattern

The AAGCCCG comparison reveals motif-specific biology. While GCCGGC has a strong binary effect and weak graded response, AAGCCCG shows the opposite: no binary effect but a stronger graded response among methylated genes. This suggests fundamentally different mechanisms: GCCGGC may primarily mark genomic regions (core vs arm), while AAGCCCG may have more gene-specific regulatory roles at high densities.

## 5. Relationship to Gatekeeper Model

The dose-response findings refine the Gatekeeper Model (H11):
- **Layer 1 (Landscape Remodeling)**: The binary 0-vs->=1 GCCGGC effect and its arm/core confound are consistent with regional landscape-level effects rather than gene-level regulation
- **Layer 2 (Protection)**: The gene body methylation effect (not promoter) suggests a mechanism distinct from sigma factor competition
- The threshold switch model is more compatible with methylation as a **regional chromatin-like mark** that establishes zones of altered transcriptional activity, rather than a gene-by-gene regulatory mechanism

## 6. Output files

### Figures
- `figures/dose_response_boxplot.pdf/svg` -- Violin+box plots of LFC by dose group (T2 and T3)
- `figures/scatter_regression.pdf/svg` -- Scatter plots with trend lines for genes with >=1 site
- `figures/location_specific.pdf/svg` -- Promoter vs gene body location-specific dose-response
- `figures/H18_comprehensive_summary.pdf/svg` -- 4-panel summary figure

### Tables
- `tables/GCCGGC_gene_site_counts.tsv` -- Per-gene site counts (8,083 genes)
- `tables/dose_response_summary.tsv` -- Dose-group summary statistics for both motifs
- `tables/statistical_tests.tsv` -- All statistical test results
- `tables/location_analysis.tsv` -- Promoter vs gene body analysis
- `tables/confound_arm_core_analysis.tsv` -- Arm/core stratified results

### Scripts
- `scripts/H18_dose_response_analysis.py`

## 7. Key numbers for reference

| Metric | Value |
|---|---|
| T1 GCCGGC sites | 1,289 |
| Genes with >=1 site within 2 kb | 3,811 (47.2%) |
| Total site-gene mappings | 5,791 |
| Spearman rho (all genes, T2) | -0.087 (p = 4.6e-14) |
| Spearman rho (>=1 site, T2) | -0.037 (p = 0.026) |
| Spearman rho (>=1 site, T3) | -0.020 (p = 0.215, NS) |
| JT trend p-value (T2) | 4.6e-14 |
| JT trend p-value (T3) | 6.8e-14 |
| Dose-response shape | 0 -> 1 -> 2: monotonic decrease; 3, 4+: plateau/rebound |
| Arm/Core stratified rho | All NS or reversed |
| AAGCCCG within-methylated rho (T3) | -0.146 (p = 1.0e-05) |

---

*Verdict: **PARTIAL** -- Statistically significant overall trend (JT p < 10^-13, Spearman all-genes rho = -0.087, p < 10^-13) but the dose-response is a threshold switch (0 vs >=1) rather than a graded proportional brake. Within-methylated genes, additional sites provide no significant further suppression for T3 (rho = -0.020, NS) and minimal for T2 (rho = -0.037, p = 0.026). Non-monotonic curve shape (4+ group rebounds) and loss of signal under arm/core stratification weaken the case for direct causation.*
