# H35: AAGCCCG Exposed TF Causal Pathway Analysis

**Date**: 2026-02-27
**Analysis directory**: `11_epigenome_integration/analysis/58_AAGCCCG_exposed_TF_causal/`
**Script**: `scripts/H35_AAGCCCG_causal_pathway.py`

## Background and Rationale

This hypothesis tests the causal pathway from the AAGCCCG methyltransferase (SC_RS17645/SCO3104) to the 62 "exposed" transcription factors, integrating findings from multiple prior analyses:

- **H5**: SC_RS17645 expression drops from T1 to T2 (LFC = -2.19), coinciding with AAGCCCG 6mA site loss (260 to 64 sites, 75% loss)
- **H21**: Genome-wide AAGCCCG temporal de-repression test was completely null (Lost vs Never: p=0.91, r=0.003)
- **H30**: AAGCCCG **DNA sequence motif** is 5.0x enriched at exposed TF TSS +/-300bp (p=1.1e-08), the strongest sequence discriminator between exposed and shielded regulators
- **H34**: 63% of exposed TFs are early responders (T1 to T2 transition), precisely when SC_RS17645 expression drops

**The key insight motivating H35**: H21 tested ALL ~800 genes with AAGCCCG loss. If only the 62 exposed TFs (which lack protection zones) are susceptible to AAGCCCG-mediated regulation, the signal would be diluted below detectability. This hypothesis restricts the test to exposed TFs only.

## Critical Methodological Discovery

A fundamental distinction emerged immediately in the analysis:

| Concept | H30 Finding | H35 Finding |
|---------|-------------|-------------|
| What was measured | AAGCCCG **DNA sequence** occurrence | AAGCCCG **methylated sites** |
| Window | +/-300bp of TSS | +/-500bp of TSS |
| Exposed TFs with hits | ~16/62 (25.8%) | **2/62 (3.2%)** |
| Enrichment | 5.0x (p=1.1e-08) | N/A (too few for enrichment test) |

**The AAGCCCG sequence motif is enriched near exposed TFs, but the vast majority of these motif instances are NOT methylated.** Of the ~1,334 AAGCCCG motif occurrences genome-wide, only 260 are methylated at T1 (19.5%). This means the sequence enrichment (H30) does not translate into a proportional methylation enrichment, severely limiting the statistical power of the causal test.

## Results

### 1. AAGCCCG Methylated Site Mapping (62 Exposed TFs)

| AAGCCCG Status | n TFs | % | Window |
|---------------|-------|---|--------|
| Lost (T1 only) | 2 | 3.2% | +/-500bp |
| Never | 60 | 96.8% | +/-500bp |
| Gained | 0 | 0.0% | +/-500bp |
| Both | 0 | 0.0% | +/-500bp |

At the broader +/-2000bp window:

| Status | n TFs | % |
|--------|-------|---|
| Lost | 6 | 9.7% |
| Never | 56 | 90.3% |

The 2 exposed TFs with methylated AAGCCCG within +/-500bp of TSS:

| Locus Tag | SCO ID | Product | AAGCCCG Sites | LFC T2vsT1 | LFC T3vsT1 |
|-----------|--------|---------|---------------|------------|------------|
| SC_RS22135 | SCO4005 | SigE family sigma factor | 1 (206bp upstream) | +0.265 | +3.232 |
| SC_RS27300 | SCO5027 | Winged helix DNA-binding protein | 1 (108bp downstream) | -2.254 | -1.101 |

These two TFs show opposite expression changes: SCO4005 is mildly upregulated (consistent with de-repression), while SCO5027 is strongly downregulated (inconsistent). The two cases are too few and too contradictory for any causal inference.

### 2. KEY TEST: Exposed-Restricted De-repression (+/-500bp)

| Metric | AAGCCCG-Lost (n=2) | AAGCCCG-Never (n=60) | Test |
|--------|--------------------|-----------------------|------|
| Median LFC T2vsT1 | -0.995 | +0.149 | |
| Mean LFC T2vsT1 | -0.995 | +0.507 | |
| Mann-Whitney U | 37.0 | | p=0.402 |
| Rank-biserial r | 0.383 | | |
| Cohen's d | -0.603 | | |
| One-sided (Lost > Never) | | | p=0.809 |

**Direction is opposite to prediction**: Lost TFs have LOWER expression (median -0.995) than Never TFs (+0.149). The de-repression hypothesis predicts Lost should be HIGHER. The test is non-significant (p=0.402) and severely underpowered (n=2 vs n=60).

Bootstrap 95% CI for median difference (Lost - Never): **[-3.118, +0.678]** (includes zero).

### 3. GCCGGC 4mC Comparison

| GCCGGC Status | n | Median LFC T2vsT1 |
|--------------|---|-------------------|
| Lost | 8 | -0.030 |
| Never | 47 | +0.116 |
| Gained | 7 | (various) |

GCCGGC Lost vs Never: p=0.916, r=0.027. No de-repression effect for GCCGGC either.

### 4. Dose-Response

| Test | rho | p |
|------|-----|---|
| n_AAGCCCG_lost vs LFC T2vsT1 (all 62) | -0.117 | 0.364 |
| n_AAGCCCG_lost vs LFC T3vsT1 (all 62) | +0.071 | 0.581 |

The dose-response is non-significant and in the wrong direction for T2vsT1. The distribution is heavily skewed: 60 TFs have 0 lost sites, 2 have 1 lost site each, making dose-response analysis uninformative.

### 5. SC_RS17645 Expression Correlation

| Gene Set | n | Median Spearman rho | Fraction Negative |
|----------|---|--------------------|--------------------|
| Exposed TFs | 62 | -0.142 | 51.6% |
| Shielded TFs | 955 | +0.233 | 44.3% |
| **MWU two-sided** | | | **p=0.346** |

Exposed TFs tend to have more negative correlation with SC_RS17645 (median rho = -0.142 vs +0.233), consistent with the prediction that when MTase drops, exposed TFs should rise. However, this difference is not statistically significant (p=0.346).

Within exposed TFs, split by AAGCCCG status:
- AAGCCCG-Lost (n=2): median rho = +0.308 (opposite to prediction)
- AAGCCCG-Never (n=60): median rho = -0.142

### 6. Temporal Specificity

| Transition | Median Diff (Lost - Never) | p-value | r |
|-----------|---------------------------|---------|---|
| T1 to T2 (MTase drops) | -1.144 | 0.402 | 0.383 |
| T2 to T3 (post-drop) | +1.845 | 0.068 | -0.767 |
| T1 to T3 (cumulative) | -0.034 | 0.609 | -0.233 |

The T2 to T3 comparison shows the most suggestive signal: Lost TFs show median +2.06 LFC T3vsT2, while Never TFs show +0.22. Bootstrap 95% CI for T3vsT2 median difference: [+0.701, +3.074] (excludes zero). However, with only n=2 in the Lost group, this finding is not robust.

### 7. H21 Comparison: Signal by Gene Set Restriction

| Gene Set | n_Lost | n_Never | Median Diff | p-value | r |
|----------|--------|---------|-------------|---------|---|
| All genes (genome-wide) | 162 | 7,431 | -0.004 | 0.610 | 0.023 |
| All regulatory (1,017) | 15 | 996 | +0.420 | 0.327 | -0.147 |
| Shielded only (955) | 13 | 936 | +0.464 | 0.152 | -0.231 |
| **Exposed only (62)** | **2** | **60** | **-1.144** | **0.402** | **0.383** |

The progressive restriction from genome-wide to regulatory-only shows a suggestive trend: shielded regulatory genes show a weak non-significant upregulation of Lost TFs (median diff +0.464, p=0.152, r=-0.231). This is interesting because the 13 Lost shielded regulators show the predicted de-repression direction, in contrast to the exposed result where the Lost direction is reversed. However, the exposed group's result is dominated by n=2 and is not interpretable.

### 8. Bootstrap Confidence Intervals

| Comparison | Bootstrap Median | 95% CI | Includes Zero? |
|-----------|-----------------|--------|----------------|
| LFC T2vsT1 median diff | -1.144 | [-3.118, +0.678] | Yes |
| LFC T2vsT1 mean diff | -1.487 | [-3.176, +0.161] | Yes |
| abs(LFC) T2vsT1 median diff | -0.219 | [-1.438, +0.965] | Yes |
| LFC T3vsT2 median diff | +1.845 | [+0.701, +3.074] | **No** |

## Interpretation

### Why the hypothesis was untestable

The fundamental problem is a gap between two levels of analysis:

1. **H30 (DNA sequence level)**: 25.8% of exposed TFs have the AAGCCCG 7-mer in their +/-300bp TSS region. This is a property of the DNA sequence that is evolutionarily conserved.

2. **H35 (methylation level)**: Only 3.2% of exposed TFs have a **methylated** AAGCCCG site within +/-500bp of TSS. This is because the SC_RS17645 methyltransferase only methylates a small fraction of its recognition sites (260 out of ~1,334 genome-wide, ~19.5%).

The hypothesis assumed that sequence enrichment would translate to methylation enrichment. It does not. The MTase appears to have site selectivity beyond simple sequence recognition, or many sites are protected from methylation by other factors (e.g., protein binding, chromatin context).

### What the data do show

1. **No evidence for AAGCCCG methylation as a direct causal mechanism for exposed TF regulation** -- with only 2/62 exposed TFs having methylated AAGCCCG near their TSS, methylation cannot be the proximate cause of expression changes in the majority of exposed TFs.

2. **The AAGCCCG sequence motif enrichment (H30) is a structural/evolutionary property** -- it indicates that exposed TF promoters have sequence features compatible with methyltransferase access, but this does not mean they are actively methylated.

3. **Suggestive shielded regulatory trend** -- among shielded regulators (n_Lost=13), Lost genes show higher LFC than Never (median diff +0.464, p=0.152), weakly consistent with de-repression. This is the opposite of the exposed result and suggests that if any regulatory genes experience AAGCCCG-mediated de-repression, it may be the shielded ones (which have more methylated sites due to higher absolute numbers).

4. **SC_RS17645 correlation difference** -- exposed TFs tend toward negative correlation with MTase expression (median rho=-0.142 vs shielded +0.233), but this is not significant (p=0.346) and may reflect the general property that exposed TFs are dynamic/activated while shielded TFs are more stable.

### Relationship to prior hypotheses

| Prior Hypothesis | Relationship to H35 |
|-----------------|---------------------|
| **H5** (SC_RS17645 drop) | Confirmed that MTase expression drops. But methylation loss does not translate to exposed TF de-repression |
| **H21** (Genome-wide null) | H35 cannot improve on H21 because the exposed TF subset has even fewer Lost genes (2 vs 162) |
| **H30** (Sequence enrichment) | Key distinction: sequence motif enrichment != methylation enrichment. The 5x sequence enrichment does not create a 5x methylation enrichment |
| **H34** (Temporal dynamics) | The simultaneous switch at T1 to T2 is not explained by AAGCCCG methylation loss |

## Verdict

**NOT SUPPORTED (UNTESTABLE)**

The causal pathway hypothesis from SC_RS17645 to exposed TFs via AAGCCCG methylation is effectively untestable with the available data because:

1. Only 2/62 exposed TFs have methylated AAGCCCG sites near their TSS (3.2%), making any statistical test severely underpowered
2. The direction of the 2-gene result is opposite to prediction (Lost TFs show repression, not de-repression)
3. The H30 sequence enrichment finding does not translate to methylation enrichment
4. Bootstrap confidence intervals include zero for all T2vsT1 comparisons

**This negative result is itself informative**: it establishes that AAGCCCG 6mA methylation is NOT the proximate mechanism by which the 62 exposed TFs are regulated, despite the strong sequence motif enrichment. The exposed TF response must be mediated by a different mechanism -- potentially the AAGCCCG sequence motif itself (as a TF binding site feature rather than a methylation substrate), or a broader chromatin/protein occupancy change that correlates with but is not caused by methylation.

## Output Files

### Tables
| File | Description |
|------|-------------|
| `tables/exposed_TF_methylation_status.tsv` | Per-TF AAGCCCG/GCCGGC methylation status, 62 rows |
| `tables/derepression_test.tsv` | Group-level LFC statistics for Lost/Never/Gained |
| `tables/SC_RS17645_correlations.tsv` | Per-exposed-TF Spearman correlation with MTase |
| `tables/statistical_tests.tsv` | All 13 statistical tests with effect sizes |
| `tables/H21_comparison.tsv` | Genome-wide vs restricted comparison |

### Figures
| File | Description |
|------|-------------|
| `figures/exposed_vs_genomewide_derepression.pdf/svg` | 4-panel: genome-wide to exposed-restricted signal comparison |
| `figures/AAGCCCG_lost_vs_never_LFC.pdf/svg` | Violin/box: LFC by AAGCCCG status (exposed TFs) |
| `figures/dose_response.pdf/svg` | Scatter: n_sites_lost vs LFC |
| `figures/SC_RS17645_correlation.pdf/svg` | MTase expression correlation analysis |
| `figures/temporal_specificity.pdf/svg` | T2vsT1 vs T3vsT2 effect comparison |
| `figures/H35_comprehensive_summary.pdf/svg` | 9-panel comprehensive summary |

## Key Statistics

| Metric | Value |
|--------|-------|
| Exposed TFs with methylated AAGCCCG (500bp) | 2/62 (3.2%) |
| Exposed TFs with AAGCCCG sequence motif (300bp, H30) | ~16/62 (25.8%) |
| Lost vs Never LFC T2vsT1 p-value | 0.402 |
| Lost vs Never rank-biserial r | 0.383 (wrong direction) |
| Lost vs Never Cohen's d | -0.603 (wrong direction) |
| Bootstrap 95% CI (median diff) | [-3.118, +0.678] |
| H21 genome-wide p-value (benchmark) | 0.610 |
| Dose-response rho | -0.117 (NS) |
| SC_RS17645 correlation difference p | 0.346 |

---

*Analysis performed: 2026-02-27*
