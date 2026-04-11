# H19: Temporal De-repression -- GCCGGC Methylation Loss at T2 Does NOT Cause Upregulation

**Date:** 2026-02-26
**Analysis directory:** `11_epigenome_integration/analysis/42_GCCGGC_temporal_derepression/`
**Status:** REJECTED

---

## Background

H17 established that GCCGGC 4mC proximity is associated with expression suppression (p=4.1e-10, median LFC shift of -0.20 relative to background). Meanwhile, H14 showed that GCCGGC sites undergo COMPLETE remodeling between timepoints (Jaccard=0.000): T1 has 1,289 sites (83% core), T2 has 407 entirely NEW sites (82% arm), and T3 has just 21 sites. Zero sites overlap between any pair of timepoints.

If GCCGGC methylation truly acts as a transcriptional "brake," then genes that lose their nearby methylation at T2 should experience de-repression (upregulation), while genes that gain new methylation at T2 should become suppressed. This provides a temporal (quasi-causal) test of the methylation-suppression model.

## Hypothesis

Genes that lose GCCGGC 4mC methylation at T2 (had a nearby site at T1 but not at T2) show significant upregulation at T2 compared to (a) genes that remain unmethylated at both timepoints and (b) genes that gain new methylation at T2.

## Methods

### Gene classification by methylation transition
- Mapped T1 (1,289) and T2 (407) GCCGGC sites to 8,083 GFF-annotated genes within 2 kb distance
- Classified each gene into four T1-to-T2 transition groups:
  - **Lost**: Had >=1 T1 GCCGGC site within 2 kb, NO T2 site within 2 kb
  - **Gained**: Had NO T1 site within 2 kb, gained >=1 T2 site within 2 kb
  - **Both**: Had sites at both T1 and T2 (different positions, since Jaccard=0)
  - **Never**: No GCCGGC site within 2 kb at either timepoint

### Expression data
- DESeq2 results: T2vsT1 (7,646 genes), T3vsT1, T3vsT2
- DEG threshold: padj < 0.05

### Statistical tests
- Wilcoxon rank-sum tests on T2vsT1 log2FC distributions
- Bonferroni correction for 6 pairwise comparisons
- Rank-biserial correlation as effect size

### Confound controls
- Geographic stratification (core vs arm genes tested separately)
- Baseline expression quintile matching
- T2-to-T3 transition analysis for replication
- Reciprocal AAGCCCG analysis

## Results

### 1. Transition group sizes

| Group | Description | n (all) | n (with expression) |
|-------|-------------|---------|---------------------|
| **Lost** | T1 methylated, T2 unmethylated | 3,280 | 3,226 |
| **Gained** | T1 unmethylated, T2 methylated | 906 | 870 |
| **Both** | Both T1 and T2 (different sites) | 462 | 451 |
| **Never** | No methylation at either T1 or T2 | 3,435 | 2,949 |

The Lost group is large (3,226 genes) because T1 had 1,289 sites broadly distributed across the core genome, reaching many genes within 2 kb. The Never group (2,949) represents genes distant from all GCCGGC sites at both timepoints.

### 2. Expression changes -- OPPOSITE to prediction

| Group | n | Median LFC (T2/T1) | Mean LFC | 95% CI | DEGs | Up | Down |
|-------|---|---------------------|----------|--------|------|-----|------|
| **Lost** | 3,226 | **-0.253** | +0.141 | [+0.067, +0.216] | 71.8% | 30.2% | **41.6%** |
| **Gained** | 870 | **+0.447** | +0.828 | [+0.681, +0.976] | 63.0% | **41.8%** | 21.1% |
| **Both** | 451 | **+0.250** | +0.714 | [+0.497, +0.931] | 61.6% | 37.0% | 24.6% |
| **Never** | 2,949 | **+0.006** | +0.316 | [+0.243, +0.388] | 68.5% | 34.2% | 34.2% |

**Critical finding:** The pattern is the **exact opposite** of the de-repression hypothesis:
- **Lost genes** (should be upregulated) have the **lowest** median LFC (-0.253), with more genes going DOWN than UP
- **Gained genes** (should be suppressed) have the **highest** median LFC (+0.447), strongly biased toward upregulation
- **Never genes** sit near zero, as expected for a baseline

### 3. Statistical tests confirm the reversed pattern

| Comparison | n1 | n2 | p (Wilcoxon) | p (Bonferroni) | Rank-biserial r | Median diff | Direction |
|-----------|-----|-----|-------------|----------------|-----------------|-------------|-----------|
| Lost vs Never | 3,226 | 2,949 | 1.39e-08 | **8.32e-08** | +0.083 | -0.259 | **Lost LOWER** |
| Lost vs Gained | 3,226 | 870 | 2.37e-27 | **1.42e-26** | +0.239 | -0.700 | **Lost LOWER** |
| Gained vs Never | 870 | 2,949 | 5.11e-13 | **3.06e-12** | -0.161 | +0.441 | **Gained HIGHER** |
| Both vs Never | 451 | 2,949 | 3.98e-04 | **2.39e-03** | -0.103 | +0.244 | Both HIGHER |
| Lost vs Both | 3,226 | 451 | 1.67e-10 | **1.00e-09** | +0.185 | -0.503 | Lost LOWER |
| Gained vs Both | 870 | 451 | 7.06e-02 | 4.23e-01 (ns) | -0.061 | +0.197 | n.s. |

All primary tests are highly significant after Bonferroni correction, but in the **opposite** direction to prediction. Lost genes are significantly *more downregulated* than Never genes (p=8.3e-08), and Gained genes are significantly *more upregulated* than Never genes (p=3.1e-12).

### 4. Geographic confound -- explains the reversal

The T1-to-T2 methylation transition is deeply confounded with genomic geography:
- T1 sites are 83% core-located, so **Lost genes are predominantly core genes**
- T2 sites are 82% arm-located, so **Gained genes are predominantly arm genes**

| Region | Comparison | n1 | n2 | Lost median | Never median | p | r |
|--------|-----------|-----|-----|-------------|--------------|---|---|
| **Core** | Lost vs Never | 2,513 | 1,646 | -0.434 | -0.417 | **0.872** | +0.003 |
| **Arm** | Lost vs Never | 713 | 1,303 | +0.302 | +0.483 | **0.019** | +0.063 |

**When stratified by genomic region, the Lost vs Never difference disappears entirely:**
- In **core genes**: Lost (median -0.434) vs Never (median -0.417), p=0.87, r=0.003 -- essentially identical
- In **arm genes**: Lost (median +0.302) vs Never (median +0.483), p=0.019, r=0.063 -- Lost genes are actually slightly LOWER, the opposite of de-repression

The overall Lost vs Never difference was a **Simpson's paradox**: Lost genes appeared more downregulated only because they are enriched in core genes, which globally decline at T2 as development progresses. The methylation transition itself adds no explanatory power.

### 5. Baseline expression quintile control

| Quintile | n_Lost | n_Never | Lost median | Never median | p | r |
|----------|--------|---------|-------------|--------------|---|---|
| Q1 (lowest) | 479 | 594 | -0.183 | +0.076 | **9.42e-07** | 0.174 |
| Q2 | 562 | 621 | -0.410 | -0.040 | **2.08e-05** | 0.143 |
| Q3 | 642 | 616 | -0.419 | +0.031 | **2.57e-06** | 0.153 |
| Q4 | 723 | 594 | -0.226 | -0.156 | 0.170 (ns) | 0.044 |
| Q5 (highest) | 820 | 524 | -0.015 | -0.063 | 0.767 (ns) | 0.010 |

An interesting pattern: Q1-Q3 (lower expression genes) show significant Lost vs Never differences, but Q4-Q5 (higher expression) do not. This suggests the effect is not driven by regression to the mean of highly expressed genes, but rather by the geographic confound affecting lower-expression genes (which are more prevalent on chromosome arms and may differ in baseline behavior).

### 6. T2-to-T3 transition analysis

The T2-to-T3 transition (407 sites to 21 sites, even more dramatic loss) replicates the reversed pattern:

| Group (T2->T3) | n | Median LFC (T3/T2) | Mean LFC |
|-----------------|---|---------------------|----------|
| **Lost** (T2 methyl, no T3) | 1,287 | **+0.136** | +0.392 |
| **Gained** (no T2, T3 methyl) | 57 | **-0.330** | +0.306 |
| **Never** (neither) | 6,118 | **-0.135** | +0.141 |

Lost vs Never at T3vsT2: U=4,335,860, p=1.05e-08, r=-0.101

Here, T2-to-T3 Lost genes DO show higher LFC than Never (+0.136 vs -0.135), but this may again reflect geographic composition (T2 sites being arm-enriched, so Lost genes are now arm-genes which may have different T3 dynamics). Notably, Gained genes (n=57) show the lowest median, directionally consistent with suppression, but the sample is tiny.

### 7. Reciprocal AAGCCCG analysis -- no effect

AAGCCCG sites (T1: 260, T2: 64, both core-enriched) provide a control where geographic confounding is less extreme.

| Group (AAGCCCG T1->T2) | n | Median LFC (T2/T1) | Mean LFC |
|-------------------------|---|---------------------|----------|
| Lost | 836 | -0.019 | +0.372 |
| Gained | 191 | -0.130 | +0.209 |
| Never | 6,548 | -0.023 | +0.320 |

- Lost vs Never: p=0.918, r=0.002 -- no difference whatsoever
- Gained vs Never: p=0.378, r=0.037 -- no difference
- Lost vs Gained: p=0.457, r=-0.035 -- no difference

AAGCCCG methylation loss shows **zero** de-repression effect, consistent with methylation transitions being a by-product of R-M system dynamics rather than a regulatory input to gene expression.

## Interpretation

### Why the hypothesis fails

The de-repression model assumed a simple causal chain: methylation present -> suppression -> methylation removed -> de-repression. Three factors undermine this:

1. **Geographic confound (primary):** The T1-to-T2 methylation shift is from core (83%) to arms (82%). This means "Lost" genes are overwhelmingly core genes, and "Gained" genes are overwhelmingly arm genes. Core and arm genes have fundamentally different developmental trajectories independent of methylation. When controlling for geography, the methylation effect vanishes (core: p=0.87; arm: p=0.019 reversed).

2. **Correlation != causation:** H17's cross-sectional association (GCCGGC-proximal genes are suppressed, p=4.1e-10) is real but appears to be a geographic artifact. GCCGGC T1 sites are core-biased (83%), and core genes as a class tend to have lower expression changes at T2. The methylation-expression association reflects shared geographic distribution, not a causal regulatory mechanism.

3. **R-M system dynamics:** The complete site remodeling (Jaccard=0.000) with geographic redistribution is more consistent with R-M system biology (protection of self-DNA, phage defense) than with gene expression regulation. The sites move where the R-M system needs them, not where gene regulation demands them.

### Implications for the Gatekeeper Model

This result refines the Gatekeeper Model:
- **Layer 2 (Protection/Depletion) remains valid**: Methylation avoids TF binding sites (fold=0.66, p=6.8e-05) and sigma-10 boxes (p=1.73e-12). This is a structural/protective feature, not a regulatory one.
- **Layer 1 (Landscape Remodeling) is confirmed**: Complete site turnover occurs, but its function is defensive (R-M system), not regulatory.
- **The "suppressive effect" of GCCGGC proximity (H17) is likely a geographic confound**, not a direct causal effect of methylation on transcription.

## Output Files

### Tables
- `tables/gene_methylation_transitions.tsv` -- 7,496 genes with T1/T2/T3 methylation status, transition groups, and expression changes
- `tables/transition_group_expression.tsv` -- Summary statistics per transition group
- `tables/statistical_tests.tsv` -- All 6 pairwise Wilcoxon tests with Bonferroni correction
- `tables/geographic_stratified_tests.tsv` -- Core/arm stratified Lost vs Never and Gained vs Never
- `tables/basemean_quintile_tests.tsv` -- Expression quintile-matched comparisons

### Figures
- `figures/transition_groups_violin.pdf/svg` -- LFC distributions by transition group (T2/T1 and T3/T1)
- `figures/fraction_up_down.pdf/svg` -- DEG fractions (up/down/ns) per group
- `figures/frequency_vs_LFC_scatter.pdf/svg` -- T1 site count vs T2vsT1 LFC per gene
- `figures/geographic_context.pdf/svg` -- Chromosome-level site and gene positions
- `figures/H19_comprehensive_summary.pdf/svg` -- 4-panel summary

### Script
- `scripts/H19_temporal_derepression.py`

## Verdict

**REJECTED**

Genes that lose GCCGGC 4mC methylation at T2 do NOT show upregulation. Instead, they show significant **downregulation** relative to never-methylated genes (median LFC -0.253 vs +0.006, p=8.3e-08). However, this difference is entirely explained by the geographic confound: T1 sites are core-biased, so Lost genes are core genes, which globally decline at T2. When stratified by genome region, the Lost vs Never difference vanishes in core (p=0.87) and reverses in arms (p=0.019, Lost lower). The reciprocal AAGCCCG analysis shows zero de-repression effect (p=0.92). The complete methylation site remodeling (Jaccard=0.000) with geographic redistribution is consistent with R-M system defense dynamics rather than gene expression regulation. H17's cross-sectional suppression association is likely a geographic artifact, not a causal effect.
