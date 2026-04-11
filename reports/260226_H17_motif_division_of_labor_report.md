# H17: Motif-Specific "Division of Labor" in Functional Gene Targeting

**Date:** 2026-02-26
**Analysis directory:** `11_epigenome_integration/analysis/40_motif_division_of_labor/`
**Status:** REJECTED

---

## Background

H15 (SUPPORTED) demonstrated that all methylation systems in *S. coelicolor* M145 universally avoid regulatory genes, but positive enrichment patterns differ by motif:
- **GCCGGC 4mC**: enriched near Transport (fold=1.45, p=0.0002) and Secondary metabolism (fold=1.58, p=0.038)
- **AAGCCCG 6mA**: trend toward Stress/Defense (fold=1.96, n.s. after Bonferroni)
- **CCGG other 4mC**: enriched near DNA/RNA metabolism (fold=2.82, p=0.005)

This suggested a "division of labor" where different R-M systems target different functional gene categories, potentially leading to distinct transcriptomic consequences.

## Hypothesis

GCCGGC-proximal genes and AAGCCCG-proximal genes show distinct expression dynamics: GCCGGC targets (Transport/SecMet enriched) undergo T2/T3 developmental upregulation, while AAGCCCG targets (Stress/Defense enriched) maintain constitutive or T1-biased expression.

## Methods

### Gene set construction
- **AAGCCCG-proximal (T1):** Extracted from pre-computed site-gene mapping (260 unique sites, within 2 kb) -> 974 unique genes
- **GCCGGC-proximal (T1):** Mapped 1,289 T1 GCCGGC sites to protein-coding genes within 2 kb using GFF annotations -> 3,695 unique genes
- **Overlap:** 451 genes proximal to both motifs ("dual-targeted")
- **Exclusive sets:** AAGCCCG-only = 523 genes, GCCGGC-only = 3,244 genes
- **Background:** 3,778 protein-coding genes not proximal to either motif

### Expression data
- DESeq2 results for T2vsT1 and T3vsT1 (7,646 genes)
- DEG threshold: |LFC| > 1 and padj < 0.05
- Functional categories assigned using same classifier as H13/H15

### Statistical methods
- Wilcoxon rank-sum tests for LFC distribution comparisons
- Chi-squared tests for temporal class distributions
- Bonferroni correction applied to all pairwise tests (n=12)

## Results

### 1. Overall expression dynamics -- opposite to prediction

The results **contradict** the hypothesis. AAGCCCG-proximal genes show *higher* upregulation (not constitutive), while GCCGGC-proximal genes show *downward* shifts relative to background:

| Gene Set | n | Median LFC T2vsT1 | Median LFC T3vsT1 | DEG Rate T2 | DEG Rate T3 | Up:Down T3 |
|----------|---|--------------------|--------------------|-------------|-------------|------------|
| **AAGCCCG-only** | 459 | **+0.073** | **+0.568** | 50.8% | 64.9% | 187:111 (63:37) |
| **GCCGGC-only** | 3,197 | **-0.171** | **-0.028** | 51.2% | 64.6% | 1098:968 (53:47) |
| **Dual-targeted** | 448 | -0.174 | **-0.373** | 51.3% | 64.5% | 125:164 (43:57) |
| Background | 3,542 | +0.115 | +0.324 | 49.4% | 61.8% | 1356:832 (62:38) |

Key findings:
- **AAGCCCG-only** genes have the **highest** median LFC at T3 (+0.568), not the predicted constitutive pattern
- **GCCGGC-only** genes have **negative** median LFC at T2 (-0.171) and near-zero at T3 (-0.028), opposite to the predicted developmental upregulation
- **Dual-targeted** genes show the **most negative** T3 median (-0.373), with a strong downward bias (43% up : 57% down)

### 2. Statistically significant but opposite direction

| Comparison | T2vsT1 p-value | T3vsT1 p-value | Direction |
|------------|---------------|----------------|-----------|
| AAGCCCG-only vs GCCGGC-only | **0.003** (r=-0.086) | **1.3e-04** (r=-0.110) | AAGCCCG higher |
| GCCGGC-only vs Background | **1.4e-10** (r=+0.090) | **4.1e-10** (r=+0.088) | GCCGGC lower |
| AAGCCCG-only vs Background | 0.88 (n.s.) | 0.35 (n.s.) | No difference |
| Dual vs Background | **0.003** (r=+0.086) | **2.7e-09** (r=+0.172) | Dual lower |

The AAGCCCG vs GCCGGC comparison is significant at both timepoints after Bonferroni correction (T3vsT1 p=0.0016), but the direction is opposite to the prediction: **AAGCCCG-proximal genes are MORE upregulated, not less**.

Notably, AAGCCCG-only genes are indistinguishable from background, while GCCGGC-only genes are significantly *suppressed* relative to background.

### 3. Functional category-specific expression

Within-category tests of GCCGGC-proximal vs background:

| Category | n (GCCGGC) | n (Bkgd) | Med LFC T3 (GCCGGC) | Med LFC T3 (Bkgd) | p-value |
|----------|-----------|----------|---------------------|-------------------|---------|
| Transport | 296 | 282 | +0.793 | +0.568 | 0.633 n.s. |
| Secondary metabolism | 119 | 65 | -0.967 | +0.078 | 0.157 n.s. |
| Primary metabolism | 752 | 797 | +0.096 | +0.472 | **0.017** |
| DNA/RNA metabolism | 154 | 156 | -0.907 | +0.073 | **0.004** |
| Stress/Defense (AAGCCCG) | 5 | 119 | +2.019 | +0.109 | 0.188 n.s. |

- GCCGGC-proximal **DNA/RNA metabolism** genes are significantly *downregulated* at T3 (median LFC -0.907 vs +0.073, p=0.004)
- GCCGGC-proximal **Primary metabolism** genes are significantly lower than background (median LFC +0.096 vs +0.472, p=0.017)
- GCCGGC-proximal Transport genes show no significant difference
- Stress/Defense genes near AAGCCCG show an interesting trend toward upregulation (median LFC +2.02) but n=5 is too small for significance

### 4. Temporal expression pattern distribution

Chi-squared tests reveal significantly different temporal class distributions:

| Comparison | Chi-squared | p-value | Significance |
|------------|------------|---------|--------------|
| AAGCCCG-only vs GCCGGC-only | 25.47 | **0.0045** | * |
| GCCGGC-only vs Background | 60.73 | **2.6e-09** | *** |
| AAGCCCG-only vs Background | 10.46 | 0.40 | n.s. |

Key temporal differences:
- **Monotonic increasing:** AAGCCCG-only 10.9% vs GCCGGC-only 7.4% (AAGCCCG has more progressive upregulation)
- **T3-up (late):** AAGCCCG-only **19.2%** vs GCCGGC-only 14.7% (AAGCCCG has more late upregulation)
- **Monotonic decreasing:** GCCGGC-only **12.4%** vs AAGCCCG-only 11.1% (GCCGGC has more progressive downregulation)
- **T3-down (late):** GCCGGC-only 11.4% vs AAGCCCG-only 8.5%
- **Constitutive:** Similar rates (21.6% vs 22.8%)

### 5. Dual-targeted gene analysis

- 451 genes are proximal to both AAGCCCG and GCCGGC sites (Jaccard index = 0.107)
- Dual-targeted genes show the **most negative** T3 LFC (median -0.373) and the highest proportion of T3-down genes (15.4%)
- Dual-targeted genes significantly differ from background at T3 (p=2.7e-09, r=0.172)
- The dual-targeted expression pattern resembles GCCGGC-only more than AAGCCCG-only, consistent with the much higher GCCGGC site density (1,289 vs 260 sites at T1)

## Interpretation

### Why the hypothesis is rejected

The prediction was based on a syllogistic argument: GCCGGC sites are enriched near Transport/SecMet genes -> Transport/SecMet genes are upregulated during development -> therefore GCCGGC-proximal genes should be upregulated. However:

1. **Compositional confound:** GCCGGC has 1,289 T1 sites covering 3,695 genes (46% of protein-coding genes). This massive footprint means GCCGGC-proximal genes largely reflect the genome average, diluting any category-specific effects.

2. **Suppressive association:** GCCGGC-proximal genes are systematically shifted toward lower expression relative to background (p<1e-10 at both timepoints). This is consistent with methylation acting as a **repressive** mark, not a permissive one.

3. **AAGCCCG-proximal genes behave like background:** With only 260 T1 sites and 974 proximal genes, the AAGCCCG footprint is much smaller. These genes are not distinguishable from background, meaning AAGCCCG proximity does not measurably alter expression dynamics.

4. **Functional enrichment =/= expression bias:** The H15 enrichment patterns (Transport near GCCGGC, Stress/Defense trending near AAGCCCG) reflect *genomic positioning*, not expression regulation. The "division of labor" in site placement does not translate to a "division of labor" in expression control.

### Unexpected finding: GCCGGC as a repressive modifier

The most robust finding is that GCCGGC 4mC proximity is associated with **suppressed expression** relative to the genome background:
- Median LFC at T2: -0.171 (GCCGGC) vs +0.115 (background), p=1.4e-10
- Median LFC at T3: -0.028 (GCCGGC) vs +0.324 (background), p=4.1e-10
- Effect size is small (rank-biserial r ~ 0.09) but highly reproducible across both timepoints

This is consistent with the Gatekeeper Model v2 (H11), where methylation acts as a constraint on gene expression rather than an activator.

## Verdict

**REJECTED** -- The hypothesis that GCCGGC targets undergo developmental upregulation while AAGCCCG targets remain constitutive is directly contradicted by the data. The actual pattern is the reverse: AAGCCCG-proximal genes show higher upregulation at T3, while GCCGGC-proximal genes are systematically suppressed. The functional enrichment "division of labor" from H15 reflects genomic positioning, not expression regulation. However, the statistically significant suppressive effect of GCCGGC proximity (p<1e-10) reinforces the Gatekeeper Model's framing of methylation as a repressive/constraining modification.

## Output Files

### Figures
| File | Description |
|------|-------------|
| `figures/violin_LFC_by_motif.pdf/svg` | LFC distributions by motif proximity (T2vsT1 and T3vsT1) |
| `figures/heatmap_category_motif.pdf/svg` | Functional category x motif proximity median LFC heatmap |
| `figures/temporal_class_distribution.pdf/svg` | Temporal expression pattern distribution by motif group |
| `figures/H17_comprehensive_summary.pdf/svg` | 4-panel summary figure |

### Tables
| File | Description |
|------|-------------|
| `tables/expression_summary_stats.tsv` | Summary statistics per motif group |
| `tables/gene_set_comparison.tsv` | Full gene-level data with group assignments (7,646 genes) |
| `tables/functional_category_by_motif.tsv` | Expression metrics by category x motif |
| `tables/dual_targeted_genes.tsv` | 448 dual-targeted genes with expression data |
| `tables/statistical_tests.tsv` | All pairwise Wilcoxon tests with Bonferroni correction |
| `tables/specific_category_tests.tsv` | Within-category expression tests |

### Script
- `scripts/H17_motif_division_of_labor.py`

---

*Analysis performed: 2026-02-26*
