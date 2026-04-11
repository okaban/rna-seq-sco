# H33: Exposed TF Genomic Neighborhood Transcriptional Impact Analysis

**Date:** 2026-02-27
**Hypothesis:** H33 - Genes in the genomic neighborhood of exposed TFs show expression changes correlated with the TF's methylation coordination type
**Status:** REJECTED (primary hypothesis), PARTIAL (activation/repression bloc directionality)

---

## Background

62 "exposed" regulatory genes lack methylation protection zones and show methylation-expression coordination (H27). H31 demonstrated that their downstream regulons cannot be directly mapped via FIMO binding motifs. In *Streptomyces*, regulatory genes often co-localize with their target genes (pathway-specific regulators near BGC clusters, local regulators near operons). This analysis tests whether genes in the genomic neighborhood of exposed TFs show expression changes correlated with the TF's methylation coordination type, providing indirect evidence for local regulatory activity.

## Methods

### Data sources
- Gene annotation: 8,275 genes from `gene_annotation_basic.tsv`
- 62 exposed TFs from `exposed_regulators_full_table.tsv`
- 1,017 regulatory genes (955 shielded + 62 exposed) from `all_genes_features.tsv`
- DESeq2 differential expression: T2vsT1, T3vsT1, T3vsT2 (7,646 genes with expression data)
- Coordination types from `coordinated_regulatory_genes.tsv`
- Co-expression modules from `coexpression_modules.tsv`

### Analysis pipeline
1. Define neighborhoods at three window sizes: +/-10kb, +/-20kb, +/-50kb
2. Exclude the TF itself and other exposed TFs from neighborhoods
3. Compare neighborhood expression magnitude (|LFC|), DEG rate, and directional concordance
4. Controls: 955 shielded TF neighborhoods + genome-wide background (7,258 non-TF genes)
5. 1,000 permutation tests (random 62-gene sampling)

---

## Key Results

### 1. Neighborhood Expression Magnitude: No Enrichment

| Window | Metric | Exposed | Shielded | p-value (MWU) |
|--------|--------|---------|----------|---------------|
| +/-10kb | Mean |LFC| T2vsT1 | 1.449 | 1.478 | 0.590 |
| +/-10kb | Mean |LFC| T3vsT1 | 1.752 | 1.778 | 0.329 |
| +/-20kb | Mean |LFC| T2vsT1 | 1.445 | 1.531 | 0.726 |
| +/-20kb | Mean |LFC| T3vsT1 | 1.727 | 1.834 | 0.942 |
| +/-50kb | Mean |LFC| T2vsT1 | 1.415 | 1.540 | 0.050 |
| +/-50kb | Mean |LFC| T3vsT1 | 1.749 | 1.854 | 0.407 |

**Result:** Exposed TF neighborhoods do NOT show higher expression change magnitude than shielded TF neighborhoods at any window size. If anything, the trend is slightly lower |LFC| for exposed neighborhoods. The +/-50kb T2vsT1 comparison approaches significance (p=0.050) but in the direction of *lower* |LFC| for exposed neighborhoods.

### 2. DEG Rate: No Difference

| Window | Exposed DEG rate T3 | Shielded DEG rate T3 | p-value |
|--------|---------------------|---------------------|---------|
| +/-10kb | 0.779 | 0.767 | 0.529 |
| +/-20kb | 0.778 | 0.770 | 0.507 |
| +/-50kb | 0.776 | 0.775 | 0.848 |

DEG rates are remarkably similar between exposed and shielded TF neighborhoods across all window sizes.

### 3. Directional Concordance: Marginally Above Random

- **Overall mean concordance:** 0.599 (exposed) vs 0.556 (shielded)
- **Mann-Whitney p = 0.191** (not significant)
- **TFs with significant concordance (binomial p<0.05):** 28/62 (45.2%)
- **Random expectation:** ~50% concordance

**By coordination type (T3):**

| Coordination type | n TFs | Mean concordance | n significant |
|-------------------|-------|-----------------|---------------|
| discordant_gain_up | 15 | 0.669 | 9 |
| ambiguous | 15 | 0.650 | 8 |
| methyl_change_no_expr_change | 5 | 0.564 | 2 |
| concordant_derepression | 9 | 0.548 | 3 |
| concordant_repression | 8 | 0.545 | 3 |
| discordant_loss_down | 10 | 0.527 | 3 |

The `discordant_gain_up` type shows the highest concordance (66.9%), indicating that TFs which gain methylation and go up tend to have neighbors that also go up. However, the overall exposed vs shielded difference is not statistically significant.

### 4. Activation vs Repression Bloc: Strong Directional Divergence

**Module assignment:** 35 activation bloc (M1-3), 26 repression bloc (M4), 1 unassigned

| Metric | Activation (M1-3) | Repression (M4) | p-value |
|--------|-------------------|-----------------|---------|
| Mean |LFC| T3vsT1 | 1.792 | 1.631 | 0.181 |
| DEG rate T3 | 0.721 | 0.815 | 0.104 |
| Neighbor direction (% up) | **63.7%** | **44.3%** | **9.2e-19** |
| Median baseMean | 90.5 | 205.5 | **2.2e-27** |

**Critical finding:** The neighbor direction divergence is extremely significant (chi-squared p=9.2e-19). Activation bloc neighbors are enriched for upregulated genes (63.7% up), while repression bloc neighbors show a slight downward bias (44.3% up = 55.7% down). This is consistent with the TFs having local regulatory influence matching their expression programs.

**BaseMean difference:** Repression bloc neighbors have 2.3-fold higher median expression than activation bloc neighbors (p=2.2e-27), consistent with the repression bloc being located in more highly expressed genomic regions.

**Product keyword differences (p<0.05):**

| Keyword | Activation | Repression | OR | p |
|---------|-----------|------------|-----|---|
| signal | 0.0% | 0.7% | 0.00 | 0.004 |
| hydrolase | 5.9% | 3.7% | 1.65 | 0.018 |
| regulator | 10.3% | 7.7% | 1.38 | 0.038 |
| ribosom | 0.4% | 1.6% | 0.23 | 0.004 |

Activation neighbors are enriched for hydrolases and regulators; repression neighbors are enriched for signaling and ribosomal genes.

### 5. Distance Decay: No Distance-Dependent Effect

| Metric | Exposed | Shielded |
|--------|---------|----------|
| Mean Spearman rho (distance vs |LFC|) | +0.026 | +0.024 |
| t-test vs 0 | p=0.197 | p=1.5e-04 |
| Exposed vs Shielded rho | p=0.997 | - |

No distance-dependent decay of expression change magnitude is observed for exposed TFs. The slight positive rho (+0.026) suggests that more distant neighbors show marginally *higher* |LFC|, opposite to the expected cis-regulatory pattern. The shielded TF rho is nearly identical (+0.024). Neither exposed nor shielded TFs exert a detectable distance-dependent transcriptional effect on their genomic neighborhoods.

### 6. BGC Proximity

- **No exposed TFs** have BGC-related product annotations (none are pathway-specific BGC regulators)
- 22/62 exposed TFs have >=2 BGC-related neighbor genes within +/-20kb, but these matches are based on broad keyword matching (synthase, synthetase, etc.) and largely reflect general metabolism rather than dedicated BGC proximity
- Notable: SC_RS03250 (SCO0266) near a lanthipeptide cluster, SC_RS36820 (SCO6924) near lanthionine synthetase genes, SC_RS39260 (SCO7411) near siderophore transporter genes

### 7. Permutation Test: No Enrichment

| Metric | Observed | Permutation mean | Permutation SD | Z-score | Empirical p |
|--------|----------|-----------------|----------------|---------|-------------|
| Mean |LFC| T3vsT1 | 1.727 | 1.826 | 0.096 | -1.03 | 0.857 |
| Concordance rate | 0.599 | 0.602 | 0.028 | -0.08 | 0.549 |
| DEG rate T3 | 0.778 | 0.747 | 0.022 | +1.43 | 0.062 |

The 1,000-permutation test confirms that exposed TF neighborhoods are NOT enriched for expression change relative to random 62-gene samples. The DEG rate shows a marginal trend (empirical p=0.062) but is not significant after considering the multiple metrics tested.

---

## Interpretation

### Primary hypothesis: REJECTED

Exposed TF neighborhoods do not show elevated expression change magnitude, higher DEG rates, or strong directional concordance compared to shielded TF neighborhoods or genome-wide background. The permutation test definitively confirms that the observed neighborhood statistics are within the null distribution.

This result has several important implications:

1. **Exposed TFs are not acting as classical local cis-regulators** that co-evolve with nearby target genes in the manner of Streptomyces pathway-specific regulators.

2. **The methylation-expression coordination of exposed TFs is gene-autonomous**, not neighborhood-spreading. The methylation response at each exposed TF affects that gene's expression individually, without propagating to nearby genes.

3. **This is consistent with H31** (no FIMO binding motifs for exposed TFs) and H32 (no genomic clustering of exposed TFs) -- the exposed TFs constitute a distributed regulatory layer rather than clustered operonic units.

### Partial support: Activation/Repression bloc directionality

The one strong positive finding is the directional divergence between activation and repression bloc neighborhoods (chi-squared p=9.2e-19). While this could reflect true local regulatory influence, a more parsimonious explanation is **genomic geography**: activation bloc TFs are enriched in chromosomal arms (which tend to be upregulated at T3), while repression bloc TFs are more centrally located (where core genes tend to be downregulated). This mirrors the Simpson's paradox identified in H19 and the general arm/core expression dynamics.

### Integration with previous findings

| Finding | Status | Reference |
|---------|--------|-----------|
| No FIMO binding motifs | Confirmed | H31 |
| No genomic clustering | Confirmed | H32 |
| No neighborhood |LFC| enrichment | **New (H33)** |
| Activation/repression directional divergence | **New (H33)** |
| Gene-autonomous methylation response | Supported | H27, H28, H33 |

The exposed TF regulatory mechanism operates through **trans-acting effects on dispersed targets** rather than cis-acting effects on nearby genes. The binding sites and target genes of these 62 regulators remain unknown, likely requiring experimental approaches (ChIP-seq, genetic knockouts) rather than computational prediction.

---

## Output Files

### Tables
| File | Description |
|------|-------------|
| `tables/per_TF_neighborhood.tsv` | Per-TF statistics: 62 TFs with neighbor counts, mean |LFC|, concordance rates, distance correlations |
| `tables/neighborhood_genes.tsv` | All 2,312 TF-neighbor pairs at +/-20kb with expression data |
| `tables/bloc_comparison.tsv` | 23 product keyword enrichments: activation vs repression bloc neighbors |
| `tables/statistical_tests.tsv` | All 27 statistical tests with test statistics, p-values, effect sizes |
| `tables/permutation_results.tsv` | 1,000-permutation test results for 3 metrics |

### Figures
| File | Description |
|------|-------------|
| `figures/neighborhood_LFC_comparison.pdf/svg` | Violin plots: exposed vs shielded vs background |LFC| at 3 window sizes |
| `figures/directional_concordance.pdf/svg` | Bar charts: concordance rates by coordination type and exposed vs shielded |
| `figures/activation_vs_repression_neighborhoods.pdf/svg` | Bloc comparison: |LFC|, DEG rate, neighbor direction |
| `figures/distance_decay.pdf/svg` | Distance from TF vs mean |LFC| decay curves |
| `figures/chromosome_neighborhood_map.pdf/svg` | Chromosome ideogram with 62 TF positions, module colors, expression direction |
| `figures/H33_comprehensive_summary.pdf/svg` | Multi-panel summary (9 panels) |

### Script
| File | Description |
|------|-------------|
| `scripts/H33_neighborhood_impact.py` | Complete analysis pipeline (Steps 1-10) |

---

## Analysis Directory

`/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/56_exposed_TF_neighborhood/`

---

*Analysis performed: 2026-02-27*
