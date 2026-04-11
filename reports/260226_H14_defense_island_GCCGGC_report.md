# H14: Defense Island--GCCGGC T3 Coupling

**日付**: 2026-02-26
**解析ディレクトリ**: `11_epigenome_integration/analysis/37_defense_island_GCCGGC/`
**判定**: **Partial** -- T3 arm-enrichment is moderate (61.9%) and one T3 site maps directly within the defense island, but 21 sites are too few for strong coupling and the dramatic geographic shift is predominantly a T1-core to T2-arm transition rather than a T3 defense-island phenomenon.

---

## 1. Background and Hypothesis

SC_RS36410 encodes a DNA cytosine methyltransferase located at ~7.6 Mb within a defense island that includes an Argonaute, PD-(D/E)XK nuclease, DUF5655, SPDY, and DnaB-like helicase (SC_RS36385--SC_RS36410). All 6 genes show coordinated T3 induction (LFC +1.13 to +5.34). SC_RS19770, the BLAST #1 hit for M.Svi27968I (E=0.007), resides at ~3.9 Mb adjacent to a pseudogene MTase (SC_RS19765, only 102 bp -- 10.5% the size of the functional copy).

**Hypothesis**: SC_RS36410 defense island activation at T3 produces arm-enriched GCCGGC methylation. The geographic shift from core (T1) to arms (T2/T3) reflects different enzymes active at different timepoints.

## 2. Results

### 2.1 GCCGGC 4mC Site Extraction

All GCCGGC-containing 4mC sites carry the TGGCCGGC motif (no other GCCGGC-containing motif was detected).

| Timepoint | GCCGGC Sites | % of Total |
|-----------|-------------|------------|
| T1 | 1,289 | 75.1% |
| T2 | 407 | 23.7% |
| T3 | 21 | 1.2% |
| **Total** | **1,717** | **100%** |

### 2.2 Geographic Distribution by Timepoint

| Timepoint | Core Sites | Arm Sites | % Core | % Arm |
|-----------|-----------|-----------|--------|-------|
| T1 | 1,070 | 219 | 83.0% | 17.0% |
| T2 | 74 | 333 | 18.2% | 81.8% |
| T3 | 8 | 13 | 38.1% | 61.9% |

**Chi-squared test** (3x2 contingency): chi2 = 596.93, p = 2.39e-130, dof = 2.

**Pairwise Fisher exact tests**:

| Comparison | Odds Ratio | p-value |
|-----------|-----------|---------|
| T1 vs T2 | 21.99 | 1.24e-127 |
| T1 vs T3 | 7.94 | 6.62e-06 |
| T2 vs T3 | 0.36 | 4.07e-02 |

Key findings:
- **T1 is overwhelmingly core** (83.0%), consistent with a housekeeping/constitutive MTase active in exponential phase
- **T2 is overwhelmingly arm** (81.8%), a dramatic geographic reversal (OR=22.0, p~10^-127)
- **T3 is moderately arm-enriched** (61.9%), significantly less than T2 (OR=0.36, p=0.041) and significantly more than T1 (OR=7.94, p=6.6e-6)
- T3 represents a partial regression toward the core pattern, not an extreme arm pattern

### 2.3 T3 Site Proximity to Defense Island

Of 21 T3 GCCGGC sites:
- **1 site falls directly within the defense island**: position 7,618,158 (minus strand, within SC_RS36410 coding region, frequency 76.7%)
- Nearest additional site: 132,982 bp away (position 7,476,470)
- Only 1 of 21 T3 sites (4.8%) is within 100 kb of the defense island

For comparison, the defense island region also hosts:
- 2 T2 sites within the DI (pos 7,617,176 at 88.4% frequency; pos 7,621,049 at 71.2%)
- 1 T1 site at 750 bp distance (pos 7,608,702 at 100% frequency)

The T3 site at 7,618,158 maps to the SC_RS36410 coding region itself (7,617,112--7,618,383), potentially representing self-methylation by the defense island MTase.

### 2.4 SC_RS19770 Locus Analysis

**SC_RS19765 pseudogene**:
- Length: 102 bp (10.5% of SC_RS19770's 975 bp)
- Gap to SC_RS19770: only 52 bp
- This is clearly a degraded duplicate, retaining only a fragment of the MTase gene

**GCCGGC sites near SC_RS19770 (within 50 kb)**:
- **21 T1 sites** exclusively -- no T2 or T3 sites
- Distances range from 1,043 bp to 46,062 bp
- All T1 sites have high frequencies (58--98%)

This is notable: the SC_RS19770 locus (core genome, ~3.9 Mb) is surrounded exclusively by T1 GCCGGC sites. SC_RS19770 expression at T1 is low (11.0 normalized counts) but detectable, with massive T3 induction (62.1, LFC=+2.32, padj=9.5e-6).

### 2.5 Defense Island Gene Expression Correlation

Mean expression (normalized counts) across timepoints:

| Gene | Product | T1 | T2 | T3 | T3/T1 LFC |
|------|---------|-----|-----|-----|-----------|
| SC_RS36385 | PD-(D/E)XK nuclease | 33.4 | 27.4 | 77.5 | +1.22 |
| SC_RS36390 | Argonaute | 17.3 | 20.3 | 114.7 | +2.73 |
| SC_RS36395 | DUF5655 | 17.2 | 46.7 | 73.4 | +2.09 |
| SC_RS36400 | SPDY | 8.3 | 6.2 | 57.8 | +2.80 |
| SC_RS36405 | DnaB-like helicase | 4.7 | 6.1 | 59.3 | +3.65 |
| SC_RS36410 | MTase (Dcm) | 1.3 | 4.2 | 64.9 | +5.59 |

**Spearman correlation matrix** (9 samples):

| | Nuclease | Argonaute | DUF5655 | SPDY | Helicase | MTase |
|---|---------|-----------|---------|------|----------|-------|
| Nuclease | 1.00 | 0.87 | 0.40 | 0.80 | 0.45 | 0.50 |
| Argonaute | 0.87 | 1.00 | 0.67 | 0.75 | 0.60 | 0.78 |
| DUF5655 | 0.40 | 0.67 | 1.00 | 0.58 | 0.95 | 0.93 |
| SPDY | 0.80 | 0.75 | 0.58 | 1.00 | 0.68 | 0.58 |
| Helicase | 0.45 | 0.60 | 0.95 | 0.68 | 1.00 | 0.87 |
| MTase | 0.50 | 0.78 | 0.93 | 0.58 | 0.87 | 1.00 |

- **Mean pairwise rho = 0.694** (moderate-to-strong co-expression)
- Two sub-clusters emerge:
  - SC_RS36395 (DUF5655) + SC_RS36405 (Helicase) + SC_RS36410 (MTase): rho = 0.87--0.95
  - SC_RS36385 (Nuclease) + SC_RS36390 (Argonaute) + SC_RS36400 (SPDY): rho = 0.75--0.87
- All 6 genes show unambiguous T3 induction, consistent with coordinated defense island activation

### 2.6 Temporal Methylation Model

**T1 (1,289 sites, 83% core)**:
- Responsible enzyme unclear. SC_RS19770 expression is very low at T1 (11.0 counts). SC_RS36410 is nearly silent (1.3 counts).
- Among highly expressed core MTases at T1, none have obvious GCCGGC specificity annotation.
- The 21 T1 sites clustering near SC_RS19770 suggest it could contribute locally, but 1,289 genome-wide sites require a more abundant enzyme.
- Candidate: a constitutive MTase with GCCGGC recognition not yet annotated in REBASE for M145.

**T2 (407 sites, 82% arm)**:
- The dramatic core-to-arm geographic shift (OR=22.0) is the most striking feature of the GCCGGC methylation system.
- T2 sites are arm-enriched in both left and right arms.
- The transition from 1,289 to 407 sites suggests either reduced overall methylation activity or increased replication-mediated dilution.

**T3 (21 sites, 62% arm)**:
- SC_RS36410 is maximally expressed at T3 (64.9 counts, LFC=+5.59), yet only 21 GCCGGC sites are detected.
- One T3 site maps within SC_RS36410 itself (self-methylation?).
- The low site count despite high expression suggests either: (a) the defense island MTase has limited substrate access in vivo, (b) the enzyme targets different modifications not captured here, or (c) detection sensitivity at T3 is limiting.

### 2.7 Distance from Defense Island

| Timepoint | Median Distance (kb) |
|-----------|---------------------|
| T1 | 3,436 |
| T2 | 988 |
| T3 | 2,603 |

T2 sites are significantly closer to the defense island than T1 sites (median 988 vs 3,436 kb), reflecting the arm-enrichment since the DI is at 7.6 Mb (right arm). T3 sites show intermediate distance (2,603 kb), not consistent with a simple "defense island proximity" model.

## 3. Figures and Tables

### Figures
| File | Description |
|------|-------------|
| `H14_defense_island_GCCGGC_analysis.pdf` | Multi-panel figure (A: geographic distribution stacked bar, B: defense island expression heatmap, C: chromosome map, D: distance histogram) |
| `defense_island_correlation_heatmap.pdf` | Full Spearman correlation heatmap of 6 defense island genes |
| `GCCGGC_position_frequency_by_timepoint.pdf` | Position vs frequency scatter for each timepoint |

### Tables
| File | Description |
|------|-------------|
| `GCCGGC_sites_by_timepoint.tsv` | All 1,717 GCCGGC sites with region and distance annotations |
| `geographic_distribution_summary.tsv` | Arm/core counts and percentages by timepoint |
| `defense_island_expression.tsv` | Normalized counts for 6 defense island genes (9 samples) |
| `defense_island_spearman_correlation.tsv` | Pairwise Spearman correlation matrix |
| `T3_GCCGGC_sites_detail.tsv` | Detailed annotation of all 21 T3 sites |

## 4. Hypothesis Assessment

### Evidence Supporting Hypothesis

1. **T3 arm enrichment exists** (61.9% arm, OR=7.94 vs T1, p=6.6e-6)
2. **One T3 site maps within the defense island** (pos 7,618,158 in SC_RS36410)
3. **All 6 defense island genes are coordinately T3-induced** (mean pairwise rho=0.694)
4. **T1-to-T2-to-T3 geographic shift** is highly significant (chi2=597, p~10^-130)
5. **SC_RS19770 locus is exclusively surrounded by T1 sites** (21 T1 sites within 50 kb)

### Evidence Against or Complicating Hypothesis

1. **T3 has only 21 sites** -- too few for a defense-island-driven methylation program
2. **T3 arm enrichment is weaker than T2** (61.9% vs 81.8%, OR=0.36, p=0.041)
3. **Only 1/21 T3 sites is near the defense island** -- no local enrichment pattern
4. **T3 median distance to DI (2,603 kb) is intermediate**, not short
5. **SC_RS36410 expression is maximal at T3** (64.9 counts) but GCCGGC sites are minimal (21), suggesting the enzyme may not be the primary GCCGGC methyltransferase
6. **The T1-enzyme paradox remains**: 1,289 T1 core sites require a highly active, constitutive MTase that is not SC_RS19770 (11 counts) or SC_RS36410 (1.3 counts)

### Revised Model

The GCCGGC methylation landscape likely involves **at least three phases**:

1. **T1 phase**: An unidentified constitutive MTase methylates ~1,300 core-genome GCCGGC sites. This enzyme is likely a highly expressed, abundant MTase whose GCCGGC specificity has not been annotated.

2. **T2 phase**: The core MTase activity declines (or is diluted), while arm-specific methylation appears (possibly replication-timing or chromatin-accessibility mediated). The dramatic geographic shift (OR=22) is the hallmark of the T1-to-T2 transition.

3. **T3 phase**: The defense island activates (all 6 genes LFC +1.1 to +5.3), and SC_RS36410 may contribute locally (self-methylation at pos 7,618,158), but the overall GCCGGC site count collapses to 21. Defense island activation correlates temporally with the late stationary phase but does not drive widespread GCCGGC methylation.

## 5. Key Insights

1. The **GCCGGC geographic flip** (83% core at T1 to 82% arm at T2, OR=22) is one of the most dramatic epigenomic transitions in this dataset, far exceeding what was reported for 6mA.
2. The defense island's **coordinated T3 induction** (6 genes, rho=0.694) is real but produces minimal detectable methylation (21 sites).
3. SC_RS19765 (pseudogene, 102 bp = 10.5% of functional gene) confirms the **MTase duplication-degradation** pattern at the SC_RS19770 locus.
4. The **identity of the T1 core GCCGGC methyltransferase** remains the central unresolved question -- 1,289 sites require an enzyme far more active than SC_RS19770 or SC_RS36410 at T1.

## 6. Script and Reproducibility

**Analysis script**: `11_epigenome_integration/analysis/37_defense_island_GCCGGC/scripts/defense_island_GCCGGC_analysis.py`

```bash
python3 11_epigenome_integration/analysis/37_defense_island_GCCGGC/scripts/defense_island_GCCGGC_analysis.py
```

---

*Generated: 2026-02-26*
