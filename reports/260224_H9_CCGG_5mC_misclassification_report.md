# H9: Systematic verification of whether CCGG "4mC" calls are actually 5mC misclassified by Nanopore basecalling

**Date:** 2026-02-24
**Analysis directory:** `11_epigenome_integration/analysis/32_CCGG_5mC_misclassification/`

---

## 1. Background and Hypothesis

In Loop 3, H7 analysis revealed that CCGG-context 4mC sites exhibit highly anomalous behavior:

- **Zero positional overlap** between timepoints (Jaccard = 0.000 for T1-T2)
- **No identified cytosine MTase** with adequate T1 expression to explain 1,516 T1 sites
- **Dramatic core-to-arm geographic shift** from T1 (77% core) to T2 (89% arms)
- This pattern is biologically implausible for enzymatic 4mC

**Hypothesis H9:** The CCGG "4mC" detections are actually 5mC produced by Dcm-like enzymes, misclassified by the Nanopore basecaller (modkit). Evidence should include: (a) modification at the C2 position (internal cytosine of CCGG, the Dcm target), (b) systematic differences in detection scores between CCGG-context 4mC and other true 4mC sites, and (c) zero 5mC signal at these same positions despite being putative 5mC.

---

## 2. Key Findings

### 2.1 C Position within CCGG: Strong C2 enrichment (84.5%)

| C Position | Count | Percentage | Interpretation |
|------------|-------|------------|----------------|
| **C2 (internal C = Dcm target)** | **1,717** | **84.5%** | Consistent with Dcm-type modification |
| C1 (first C) | 162 | 8.0% | CCGG and GGCCGG motif sites |
| not_in_CCGG | 153 | 7.5% | CCGC motif sites (different context) |

The dominant TGGCCGGC motif (n=1,717) shows **100% methylation at the C2 position** (position 4 in the octamer = internal C of the GCCGGC palindrome). This is precisely the position where a Dcm-like or Type II cytosine methyltransferase would act.

**By motif:**
- TGGCCGGC (n=1,717): 100% at C2 -- the internal C of the embedded GCCGGC palindrome
- CCGG (n=146): 100% at C1
- GGCCGG (n=16): 100% at C1
- CCGC (n=153): 100% not_in_CCGG (separate context)

### 2.2 5mC Signal at CCGG 4mC Positions: Effectively Zero

This is the single most informative test. If modkit were misclassifying 5mC as 4mC, we would expect at least partial 5mC signal at the same positions.

| Metric | Value |
|--------|-------|
| Total replicate-level records examined | 18,282 |
| Mean 5mC frequency at CCGG 4mC positions | **0.01%** |
| Records with any 5mC signal (mod_freq > 0) | 26 / 18,282 (0.14%) |
| Records with 5mC > 5% | 16 / 18,282 (0.09%) |
| Maximum 5mC frequency observed | 20% |
| Mean 4mC frequency at same positions | **79.9%** |
| Correlation between 4mC and 5mC frequency | r = -0.015 (no relationship) |

**Interpretation:** Modkit (ONT basecaller) decisively classifies these positions as 4mC, not 5mC. The complete absence of 5mC signal is **inconsistent with simple 5mC misclassification**, because even erroneous classification should produce some proportion of reads called as 5mC rather than 4mC. The model confidently assigns 4mC at ~80% frequency, with no ambiguity toward 5mC.

### 2.3 Modification Frequency Comparison

| Group | T1 n | T1 Mean Freq | T1 Median | T2 n | T2 Mean Freq |
|-------|------|-------------|-----------|------|-------------|
| CCGG 4mC | 1,516 | 82.7% | 85.5% | 486 | 77.5% |
| Non-CCGG 4mC | 471 | 74.9% | 75.0% | 1,960 | 79.8% |
| 6mA | 1,934 | 59.2% | 56.8% | 2,120 | 58.6% |

CCGG 4mC sites have **higher** modification frequency than both non-CCGG 4mC and 6mA sites (T1: CCGG vs non-CCGG p < 1e-30, CCGG vs 6mA p < 1e-300). This is the **opposite** of what misclassification would predict -- artifacts typically show lower/borderline scores, not the highest frequencies observed.

**Per-replicate analysis** (from raw pileup files, 9 samples):
- CCGG 4mC: mean 81.5% (T1), 79.9% (T2), 79.2% (T3)
- Non-CCGG 4mC: mean 71.4% (T1), 70.7% (T2), 67.5% (T3)
- 6mA: mean 52.4% (T1), 53.2% (T2), 54.2% (T3)

The CCGG 4mC signal is consistently the strongest of all modification types across all timepoints and all replicates.

### 2.4 Replicate Consistency: No Difference

All high-confidence sites (all groups) show 100% 3-replicate detection. This is because the current dataset uses MIN_REPS=3 filtering, so by definition all HC sites have 3/3 replicate support. There is no difference in replicate consistency between CCGG 4mC, non-CCGG 4mC, and 6mA.

### 2.5 CCWGG (Dcm Substrate) Context: Zero Match

| Context | Genome Count | Methylated CCGG Sites |
|---------|-------------|----------------------|
| CCGG (4bp) | 151,175 | 1,879 (92.5% of CCGG-related) |
| CCWGG (CCAGG+CCTGG, 5bp) | 51,504 | **0 (0.0%)** |
| GCCGGC (6bp palindrome) | 15,726 | 1,717 (91.4%) |
| TGGCCGGC (8bp) | 2,093 | 1,717 (82.0%) |

**Critical finding:** NONE of the CCGG 4mC sites are in CCWGG context. E. coli Dcm methylates CC**W**GG (W = A or T), but these sites are exclusively in CC**G**G context (specifically GCCGGC). This **rules out Dcm-type activity** as the source, since Dcm substrates and CCGG sites are non-overlapping motifs.

### 2.6 Flanking Context Analysis

The dominant 6-bp context around methylated CCGG sites:

| 6-bp Context | Count | Percentage |
|-------------|-------|------------|
| **GCCGGC** | 1,717 | **91.4%** |
| GCCGGG | 76 | 4.0% |
| CCCGGG | 41 | 2.2% |
| ACCGGG | 26 | 1.4% |
| TCCGGG | 19 | 1.0% |

The overwhelmingly dominant context is **GCCGGC**, which is a perfect palindromic sequence. This is consistent with recognition by a Type II restriction-modification system targeting GCCGGC, not a Dcm-type enzyme.

---

## 3. Evidence Summary

| Evidence Line | Observation | Supports H9? | Weight |
|--------------|-------------|-------------|--------|
| C2 position (Dcm target) | 84.5% at C2 | YES | Moderate |
| 5mC signal at same sites | 0.01% mean (effectively zero) | **NO** | **Strong** |
| Modification frequency | 82.7% mean (highest of all types) | **NO** | **Strong** |
| Replicate consistency | 100% 3-rep for all groups | Uninformative | -- |
| CCWGG context | 0% at Dcm substrate | **NO** | **Strong** |
| GCCGGC context | 91.4% at palindromic site | **NO** (suggests specific enzyme) | **Strong** |
| Temporal overlap (H7) | Jaccard = 0.000 | YES | Moderate |
| MTase expression (H7) | No identified C-MTase | YES | Moderate |
| Core-arm shift (H7) | 77% core (T1) to 89% arms (T2) | YES | Moderate |

**FOR H9 (3 lines):** C2 position, temporal instability, no known MTase
**AGAINST H9 (4 lines):** Zero 5mC signal, high 4mC frequency, not CCWGG context, GCCGGC palindrome

---

## 4. Verdict: H9 is **NOT SUPPORTED**

Despite the suggestive positional evidence (C2 methylation), the hypothesis that CCGG "4mC" calls are misclassified 5mC is **not supported** by the signal-level analysis:

1. **Zero 5mC co-detection** -- If modkit were confusing 5mC with 4mC, some fraction of reads should still be called as 5mC. Instead, 5mC frequency is 0.01% while 4mC is 79.9% at the same sites. The model is decisively and consistently calling 4mC.

2. **Highest modification frequency of any type** -- At 82.7% mean weighted frequency (T1), CCGG 4mC sites show the strongest modification signal in the entire dataset, exceeding both non-CCGG 4mC (74.9%) and 6mA (59.2%). Misclassification artifacts are expected to have weak, borderline signals.

3. **Not CCWGG context** -- The sites are in CCGG/GCCGGC context, not CCWGG. E. coli Dcm targets CCWGG (5bp), which is a completely different motif from CCGG (4bp). This eliminates Dcm as the candidate enzyme.

4. **GCCGGC palindrome recognition** -- 91.4% of sites are in the palindromic GCCGGC context, which is a canonical Type II R-M recognition site. This suggests a specific enzyme activity, not random misclassification.

---

## 5. Revised Interpretation

The CCGG 4mC phenomenon at GCCGGC sites is **genuine N4-methylcytosine** (or at minimum, a modification that the ONT model confidently and reproducibly identifies as 4mC), but the responsible enzyme and the biological explanation for the anomalous temporal behavior remain unknown. The key paradox from H7 persists:

- **Strong, reproducible signal** (80%+ frequency, 3/3 replicates) argues for genuine modification
- **Zero temporal overlap** and **geographic instability** argue against maintained enzymatic methylation
- **No identified MTase** with matching specificity and adequate expression

Possible explanations that remain viable:

1. **Uncharacterized MTase** with N4-cytosine specificity for GCCGGC, expressed at low but sufficient levels from an unannotated gene or a known gene with undetermined specificity
2. **Population-level stochastic methylation** where the enzyme methylates a subset of GCCGGC sites in each cell, with different subsets in different growth phases, appearing as complete non-overlap in bulk sequencing
3. **Coverage/threshold effects** where true modification is present at all sites but only crosses the detection threshold at different positions in different timepoints due to varying coverage depth (T1: 59x combined, T2: 41x, T3: 34x)

---

## 6. Output Files

### Figures

| File | Description |
|------|-------------|
| `fig1_C_position_in_CCGG.pdf/svg` | C position distribution within CCGG motif (overall, by timepoint, by motif) |
| `fig2_score_frequency_comparison.pdf/svg` | Violin plots comparing modification frequency across groups and timepoints |
| `fig3_CCGG_flanking_context.pdf/svg` | Flanking base analysis around methylated CCGG sites |
| `fig4_replicate_consistency.pdf/svg` | Replicate consistency (n_reps) comparison |
| `fig5_comprehensive_summary.pdf/svg` | Multi-panel summary figure |
| `fig6_5mC_signal_check.pdf/svg` | 4mC vs 5mC signal at CCGG positions; TGGCCGGC motif position |

### Tables

| File | Description |
|------|-------------|
| `CCGG_C_position_analysis.tsv` | Per-site C position within CCGG for all 2,032 sites |
| `replicate_level_scores.tsv` | Per-replicate modification data from raw pileup files (65,949 records) |
| `score_comparison_statistics.tsv` | Mann-Whitney U and KS test results for frequency comparisons |
| `weighted_freq_comparison_statistics.tsv` | High-confidence level frequency comparison statistics |
| `replicate_consistency_statistics.tsv` | Replicate consistency (n_reps) statistical comparisons |
| `CCGG_flanking_context.tsv` | 6-bp flanking context for each CCGG 4mC site |
| `dual_4mC_5mC_signal.tsv` | Paired 4mC/5mC measurements at CCGG positions |
| `TGGCCGGC_position_analysis.tsv` | Position-within-motif analysis for TGGCCGGC sites |
| `H9_evidence_summary.tsv` | Evidence summary table |
| `H9_context_analysis_summary.tsv` | Comprehensive context analysis summary |

### Scripts

| File | Description |
|------|-------------|
| `scripts/H9_CCGG_5mC_misclassification.py` | Main analysis script |

---

## 7. Methods

### Data sources
- High-confidence methylation sites: `11_epigenome_integration/analysis/01_integration/high_confidence_sites_weighted.csv` (12,429 sites, MIN_REPS=3, weighted frequency >= 50%)
- 4mC motif assignments: `11_epigenome_integration/analysis/23_expanded_motif_search/4mC_motif_assignment.csv` (2,693 sites)
- Raw pileup files: 9 samples from `methyl/260102_M145/analysis/cursor_results/20260108/pileup/` (modkit bedMethyl format)
- Reference genome: *S. coelicolor* A3(2) M145 (NC_003888.3, GCF_000203835.1)

### C position analysis
For each CCGG-related 4mC site, the genomic position was mapped to the nearest CCGG tetranucleotide. The C position within CCGG was determined: C1 (first C, position 0) or C2 (second C, position 1 = internal cytosine). For TGGCCGGC sites, the position was verified by direct genomic coordinate comparison against TGGCCGGC/GCCGGCCA motif occurrences.

### 5mC co-detection analysis
For all 2,032 CCGG-context 4mC positions, the corresponding 5mC signal was extracted from the same pileup files. Modkit reports both 4mC (code "21839") and 5mC (code "m") independently for each cytosine, allowing direct comparison of competing modification calls at the same position.

### Statistical tests
Mann-Whitney U tests (two-sided) and Kolmogorov-Smirnov tests were used for frequency distribution comparisons. All tests were performed using scipy.stats. P-values < 0.05 were considered significant.
