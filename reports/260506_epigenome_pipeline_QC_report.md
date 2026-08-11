# Epigenome Pipeline QC Report — *S. coelicolor* M145 Nanopore Methylation

**Date:** 2026-05-06  
**Study:** *S. coelicolor* A3(2) M145 — Nanopore m4C/6mA methylation, 3 timepoints (T1/T2/T3) × 3 biological replicates  
**Samples:** 1-1, 1-2, 1-3 (T1) / 2-1, 2-3, 2-4 (T2) / 3-2, 3-3, 3-4 (T3)  
**Key data file:** `11_epigenome_integration/analysis/01_integration/high_confidence_sites_weighted.csv` (12,429 site-timepoint entries, all n_reps=3)

---

## Check 1: Coverage and Depth — ⚠️ WARN (1 sub-issue)

**Post-resequencing coverage per sample** (from `reports/260224_T3_resequencing_reanalysis_report.md`):

| Sample | TP | Coverage | Status |
|--------|----|----------|--------|
| 1-1 | T1 | 87.7x | ✓ |
| 1-2 | T1 | ~82x | ✓ |
| 1-3 | T1 | ~125x | ✓ |
| 2-1 | T2 | 99.2x | ✓ |
| 2-3 | T2 | 89.9x | ✓ |
| 2-4 | T2 | 86.1x | ✓ |
| **3-2** | **T3** | **57.6x** | ⚠️ lower |
| 3-3 | T3 | 79.6x | ✓ |
| 3-4 | T3 | 75.4x | ✓ |

All samples pass the minimum coverage threshold applied in analysis (MIN_COVERAGE=5 per sample, MIN_TOTAL_COVERAGE=30 combined, producing a de facto mean of 60–80x per timepoint). All 12,429 high-confidence sites have total_coverage ≥30x (range: 30–263x per site, from `high_confidence_sites_weighted.csv`).

**⚠️ WARN — T3 systematically lower coverage.** Mean coverage per site: T1=80.9x, T2=71.5x, T3=60.7x (calculated from `high_confidence_sites_weighted.csv` column 7). Sample 3-2 (57.6x) remains the lowest in the study after resequencing, even though it improved dramatically from the original 14.3x. The T3 4mC site count (1,647) is ~17% lower than T1 (1,987) and ~33% lower than T2 (2,446), which may partly reflect this residual coverage gap rather than purely biological decrease. The manuscript should acknowledge that T3 4mC counts are a conservative lower-bound estimate.

**⚠️ WARN — Coverage discrepancy for 1-1.** The early report `260202_epigenome_transcriptome_integration_report.md` states 1-1 = 37.1x, while `260224_T3_resequencing_reanalysis_report.md` shows 87.7x. T1 was *not* resequenced. This discrepancy (a ~2.4× difference) is unexplained and likely reflects different measurement methods (e.g., mean genome depth from `samtools coverage` vs. `mosdepth median`), but should be reconciled with a consistent metric before manuscript submission.

---

## Check 2: Alignment Quality — ✅ PASS

**RNA-seq (HISAT2 v2.7.11b, reference GCF_000203835.1):** All 9 samples 98.41–98.70% overall alignment rate (from `02_alignment/analysis/02_alignment_260127_v1/logs/*.hisat2.summary.txt`). Concordant unique alignment rate: 96.42–96.66% for all samples. No outliers.

**Nanopore (minimap2, reference NC_003888.3):** Mapping rates 85.1–92.1% (from `260224_T3_resequencing_reanalysis_report.md`). Acceptable for long-read ONT data. Sample 3-2 has the lowest rate (85.1%), consistent with its resequencing run composition.

**⚠️ WARN — Reference genome accession inconsistency.** The study context specifies alignment to **AL645882.2** (EMBL accession), but all analysis scripts, data files, and reports consistently use **NC_003888.3** (NCBI RefSeq). These are the same sequence (8,667,507 bp, *S. coelicolor* A3(2) complete genome, confirmed by `11_epigenome_integration/data/NC_003888.3.fna`), but the two accessions have different coordinate systems for gene annotations. If AL645882.2 coordinates were used at any point in the pipeline (e.g., gene body annotation from EMBL GFF), there could be positional offsets affecting TSS-distance calculations. The manuscript methods section must state a single, unambiguous accession.

**Note:** No BAM flagstat files are present in the rna-seq directory tree for the Nanopore alignments; they reside outside at `/Users/okaban/bioinfo/methyl/260102_M145/`. The alignment stats cited here come from the T3 resequencing report.

---

## Check 3: Replicate Consistency — ✅ PASS (with caveats)

All 12,429 entries in `high_confidence_sites_weighted.csv` have **n_reps=3** (verified: `awk -F',' 'NR>1{print $6}' ... | sort | uniq -c` → 12429 × "3"). The current pipeline enforces strict 3/3 replicate agreement (MIN_REPS=3 as of post-resequencing reanalysis documented in `260224_T3_resequencing_reanalysis_report.md`). This is the most stringent reproducibility criterion possible.

Per-timepoint n_reps=3 counts: T1=3,921, T2=4,566, T3=3,942 — all three timepoints have robust 3/3 coverage.

**No formal quantitative replicate correlation matrix** (Pearson r or Spearman ρ between all replicate pairs) was found in any report. The 3/3 concordance metric is a strict binary filter, not a graded similarity score. For a manuscript, a supplementary figure showing pairwise replicate scatter plots (or heatmap of Pearson r values between replicates) would strengthen credibility, particularly for T3 where one replicate was resequenced.

---

## Check 4: Modification Model and Probability Thresholds — ⚠️ WARN

**Modkit pileup thresholds** (from `260204_summary_manuscript.md`, lines 975–978):
- 4mC threshold: percentile 0.9375 (modkit default)
- 6mA threshold: percentile ~0.64 (modkit default)
- Post-filter: MIN_MOD_FREQ=50%, MIN_TOTAL_COVERAGE=30x, MIN_REPS=3

Column parsing in `weighted_methylation_analysis.py` (and `unweighted_methylation_analysis.py`): parts[3]=mod_code, parts[5]=strand, parts[9]=coverage, parts[10]=mod_freq_percent. This is **consistent** with the modkit bedMethyl format (N_valid_cov=col10, fraction_modified_as_percent=col11 in 0-indexed). ✅

**⚠️ WARN — Basecalling model documentation inconsistency.** The comment in `11_epigenome_integration/scripts/5mC_vs_4mC_at_CCGG.py` (line 6–7) explicitly states: *"Dorado model: dna_r10.4.1_e8.2_400bps_sup@v5.2.0_4mC_5mC@v1 → Both 4mC and 5mC are called."* The `4mC_5mC@v1` modification model detects 4mC (CHEBI code 21839) and 5mC (code "m"), **not** 6mA (code "a"). Yet 6mA sites are present throughout the dataset (e.g., 1,934 6mA sites at T1). This implies a second modification model was applied during basecalling for 6mA detection, but this is not documented in any script, CLAUDE.md, or report. The `260207_TF_binding_site_methylation_report.md` mentions only "R10.4.1, SUP basecalling v5.2.0" with no 6mA model specification. **Before manuscript submission, the exact Dorado/MinKNOW command used, including all `--modified-bases` flags, must be documented.** The probable model is `dna_r10.4.1_e8.2_400bps_sup@v5.2.0_4mC@v1` + `_6mA@v2` applied simultaneously (or a combined `4mC_6mA` model), but this cannot be confirmed from the rna-seq project directory alone.

**⚠️ WARN — Pileup directory path inconsistency across scripts.** Two different input paths are used:
- `cursor_results/20260108/pileup` — used by `weighted_methylation_analysis.py`, `unweighted_methylation_analysis.py`, `reanalyze_T3_2reps.py` (the three primary analysis scripts)
- `analysis/pileup` — used by `5mC_vs_4mC_at_CCGG.py`

The `260224_H9_CCGG_5mC_misclassification_report.md` (line 193) confirms the raw pileup files reside in `cursor_results/20260108/pileup/`. It is therefore likely that `5mC_vs_4mC_at_CCGG.py` reads from a **different (older) version** of the pileup files — potentially pre-resequencing pileups that do not include the additional T3 data from February 2026. Any claims from that script (e.g., the 5mC-at-CCGG analysis) should be verified as using the correct, post-resequencing pileup files.

---

## Check 5: Site Counts and Motif Calling — ✅ PASS (with one open biological issue)

**High-confidence site counts** (from `high_confidence_sites_weighted.csv`, n_reps=3, ≥50% freq, ≥30x total coverage):

| Timepoint | 4mC | 6mA | Total |
|-----------|-----|-----|-------|
| T1 | 1,987 | 1,934 | 3,921 |
| T2 | 2,446 | 2,120 | 4,566 |
| T3 | 1,647 | 2,295 | 3,942 |
| **Total** | **6,080** | **6,349** | **12,429** |

**Census of unique positions** (from `260224_T3_resequencing_reanalysis_report.md`): 2,693 unique 4mC positions, 3,248 unique 6mA positions.

**Motif attribution for 4mC** (2,693 unique positions):
- TGGCCGGC (GCCGGC context): 1,717 (63.8%)
- AAGCCCG: 814 (30.2%)
- CCGG other: 146 (5.4%)

**⚠️ WARN — GCCGGC 4mC collapse and non-overlap (H7 paradox, open issue).** The 60× genome-wide collapse of GCCGGC-4mC from T1 (4.10% of sites methylated) to T3 (0.067%) is documented in `reports/260504_v_defense_MTase_expression_dynamics_report.md`. More critically, `260224_H7_CCGG_MTase_paradox_report.md` establishes **zero overlap** between T1 and T2 GCCGGC 4mC sites (Jaccard=0; T1=1,516 sites, T2=486 sites, all mutually exclusive). The paradox report (H7, hypothesis (d)) acknowledges a detection artifact is "partially supported" but not ruled out. The interpretation of these sites as true biological methylation vs. model-specific signal noise must be resolved before publication, as it directly affects claims about GCCGGC-m4C dynamics.

---

## Check 6: Known Outlier — T3 Sample 3-2 ✅ RESOLVED; T1 Rep1 Claim ⚠️ UNVERIFIABLE

**T3 sample 3-2 (the documented outlier):** Fully resolved. Original coverage 14.3x → post-resequencing 57.6x (+303%). Impact: 4mC T3 sites increased by +570 (+52.9%), from 1,077 to 1,647. Documented in detail in `reports/260224_T3_resequencing_reanalysis_report.md` and `reports/260129_T3_2rep_comparison_report.md`. Current pipeline uses the resequenced data with MIN_REPS=3. ✅

**T1 replicate 1 (sample 1-1) — "1,249 vs ~3,700 sites" claim:** This specific per-sample site count discrepancy **cannot be confirmed** from any file within the rna-seq project directory. Sample 1-1 has adequate coverage (37.1x in early report, 87.7x in T3 resequencing report). No report in `reports/` flags sample 1-1 as an outlier or documents it having ~1,249 sites vs ~3,700 for the other T1 replicates. The raw modkit pileup BEDs (at `/Users/okaban/bioinfo/methyl/260102_M145/analysis/cursor_results/20260108/pileup/1-1_pileup.bed`) are required to check this per-sample site count and are not accessible within the rna-seq directory. **If this is a known issue, it must be explicitly documented in a report.**

---

## Check 7: Integration with RNA-seq — ✅ PASS

**Sample set consistency (critical):** The RNA-seq alignment BAMs in `02_alignment/analysis/02_alignment_260127_v1/bam/` use **exactly the same 9 samples** as the methylation analysis: M145_1_1/1_2/1_3, M145_2_1/2_3/2_4, M145_3_2/3_3/3_4. The early `260204_summary_manuscript.md` table (line ~963) lists a different set (including M145_2-2 and M145_3-1), but this is an **outdated section** — the actual HISAT2 logs confirm the current correct set. Both datatypes are paired at the timepoint level. ✅

**TSS file:** Jeong et al. 2016 dRNA-seq data loaded as `11_epigenome_integration/analysis/18_tss_analyses/jeong2016_all_tss.csv` (3,570 entries). Primary TSSs (category=P) mapped to SC_RS gene IDs (`tss_analyses.py` lines 118–131). ✅

**Regulatory gene list (1,055 genes):** Extracted from GFF annotation using 25 TF family keywords (confirmed in `reports/260224_H6_genomewide_TF_screen_report.md`). This is the in silico annotation-based list, not a manually curated list. ✅

**⚠️ WARN — Stale sample names in 260204_summary_manuscript.md.** Lines 959–966 of that report still list M145_2-2 (40.5x) and M145_3-1 (33.9x) as Nanopore samples, and the RNA-seq table at line ~2730 lists M145_2-2 and M145_3-1. These samples are NOT in the final analysis. This section of the manuscript must be updated to reflect the correct sample set (2-4 and 3-4 replace 2-2 and 3-1 respectively). Reviewers will notice the mismatch.

---

## Check 8: Integrity Script — ❌ NOT FOUND

`analysis/check_integrity.py` does **not exist** at any location within the rna-seq directory tree. No integrity check script was found (`find` returned no results). The QC/integrity checking for this pipeline has been performed through the reports system (multiple markdown reports in `reports/`), not through an automated script. No automated file integrity checks (e.g., checking file sizes, checksums, or expected site count ranges) are implemented.

---

## Check 9: Errors and Warnings in Logs/Outputs — ✅ PASS (no hard failures)

No `ERROR` or `FAIL` flags were found in `11_epigenome_integration/analysis/01_integration/`. RNA-seq HISAT2 logs contain no error messages. All `fastp` QC passed (Q30: 94.8–95.1%, GC ~68%, retention ≥93.4% for all samples, from `260204_summary_manuscript.md`).

**Documented biological/analytical warnings** (from reports, not pipeline errors):
- `260224_H7_CCGG_MTase_paradox_report.md`: zero inter-timepoint GCCGGC 4mC site overlap (see Check 5)
- `260226_H18/H19_reports`: GCCGGC-expression association is partially a geographic (core/arm) artifact
- `260224_H9_report`: 5mC misclassification concern for CCGG sites — assessed and resolved (CCGG 4mC mean frequency 82.7% > expected for artifacts)
- `260129_T3_2rep_comparison_report.md`: "low-coverage samples can dramatically bias methylation analysis" (resolved by resequencing)
- `reports/260504_methylation_correlation_sensitivity_allangles_report.md`: sensitivity analysis shows absolute T1 4mC × LFC correlation is not null (Spearman ρ borderline at ±200–500 bp TSS window) — the "no correlation" headline result is not fully robust across all analysis angles

---

## Summary Table

| # | Area | Status | Key Finding |
|---|------|--------|-------------|
| 1 | Coverage / depth | ⚠️ WARN | T3 mean coverage 60.7x < T1 80.9x; 1-1 coverage value inconsistent across reports (37x vs 88x) |
| 2 | Alignment quality | ⚠️ WARN | RNA-seq 98.4–98.7% ✓; Nanopore 85–92% ✓; but AL645882.2 vs NC_003888.3 accession mismatch must be reconciled |
| 3 | Replicate consistency | ✅ PASS | All 12,429 sites require 3/3 replicates; no quantitative replicate r-matrix documented |
| 4 | Mod model / thresholds | ⚠️ WARN | 6mA model identity undocumented (4mC_5mC@v1 cannot call 6mA); two conflicting pileup directory paths across scripts |
| 5 | Site counts / motif | ⚠️ WARN | Site counts self-consistent; GCCGGC 4mC non-overlap paradox (H7) unresolved — potential artifact must be addressed |
| 6 | T1 rep1 outlier | ⚠️ WARN | T3 3-2 outlier resolved ✓; T1 rep1 "1,249 sites" claim unverifiable from rna-seq directory |
| 7 | RNA-seq integration | ✅ PASS | Sample sets matched; old manuscript has stale sample IDs (2-2, 3-1) that must be corrected |
| 8 | Integrity script | ❌ FAIL | `analysis/check_integrity.py` does not exist |
| 9 | Errors / warnings | ✅ PASS | No pipeline errors; biological/analytical caveats documented in reports |

---

## Priority Action Items Before Manuscript Submission

1. **[HIGH] Document the exact 6mA basecalling model** — record the precise Dorado/MinKNOW command used, including all `--modified-bases` flags, and add to Methods.
2. **[HIGH] Reconcile reference genome accession** — choose AL645882.2 or NC_003888.3 and use it consistently throughout manuscript, scripts, and supplementary tables.
3. **[HIGH] Address GCCGGC 4mC non-overlap paradox (H7)** — either provide a mechanistic explanation or explicitly acknowledge it as an open question with potential technical contributions.
4. **[MEDIUM] Update 260204_summary_manuscript.md sample tables** — replace M145_2-2/M145_3-1 with M145_2-4/M145_3-4 throughout.
5. **[MEDIUM] Verify 5mC_vs_4mC_at_CCGG.py pileup path** — confirm this script uses the post-resequencing pileup files (`cursor_results/20260108/pileup/`), not the old `analysis/pileup/` directory.
6. **[MEDIUM] Document T1 rep1 per-sample site counts** — verify whether sample 1-1 has an anomalously low per-sample site count and document in QC.
7. **[LOW] Add quantitative replicate concordance figures** — pairwise Pearson r scatter plots or heatmap for all 9 × 9 replicate pairs for both 4mC and 6mA.
8. **[LOW] Reconcile T1 sample 1-1 coverage value** — standardize measurement method and update all relevant tables.

---

*Generated: 2026-05-06 | Based on: `11_epigenome_integration/`, `02_alignment/`, `reports/` directories*
