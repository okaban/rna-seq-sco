# Numerical Fixes — Manuscript Verification Summary

**Date:** 2026-05-24  
**Data source:** `/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/`  
**Files modified:** `manuscript/01_results.md`, `manuscript/full_manuscript.md`

---

## Fix #1 — Shielded TF count

| | Value |
|---|---|
| Manuscript | n = 998 |
| Figure label | n = 994 |
| **Data (all_genes_features_unified_n57.tsv)** | Total TFs = 1055, exposed = 57, shielded = **998** |

**Verdict: Manuscript (998) is CORRECT. No change needed.**

The figure displays 994 because 4 shielded genes have `NaN` in `nearest_methyl_distance` and are excluded from the Mann-Whitney test. The true shielded count is 998.

---

## Fix #2 — Mann-Whitney p-value ✅ FIXED

| | Value |
|---|---|
| Manuscript (old) | *p* = 5.3 × 10⁻²⁸ |
| Figure | *p* = 1.7 × 10⁻²⁵ |
| **Recomputed from data** | *p* = **1.69 × 10⁻²⁵** |

**Data:** `scipy.stats.mannwhitneyu(shielded_dist, exposed_dist, alternative='greater')` on `nearest_methyl_distance`, n_shielded = 994 (non-NaN), n_exposed = 57.

**Action:** Changed `5.3 x 10^-28` → `1.7 x 10^-25` in both manuscript files.

**Locations changed:**
- `01_results.md` line 39
- `full_manuscript.md` line 61

---

## Fix #3 — Distance AUC (pooled)

| | Value |
|---|---|
| Manuscript | 0.917 |
| Fig 3b (bar chart, ROC_analysis.tsv T1 subset, n=364) | 0.923 |
| Fig 4 (dynamic ROC on full data, n=1017) | 0.910 |
| **cross_validation.tsv: 5-fold CV mean AUC** | **0.9166 → 0.917** |

**Verdict: Manuscript (0.917) is CORRECT — it reports the 5-fold cross-validated mean AUC from `cross_validation.tsv`. No change needed.**

Note: The figure discrepancy exists because:
- Fig 3b bar chart reads from `ROC_analysis.tsv` (T1 timepoint subset, n_exposed=11, n_shielded=353) → 0.923
- Fig 4 ROC curve computes dynamically on all valid data (1017 rows) → 0.910
- The manuscript's 0.917 (CV mean) is the statistically appropriate value for model performance.

---

## Fix #4 — Early response gene count ✅ FIXED

| | Value |
|---|---|
| Manuscript (old) | 39/57 (63%) |
| Figure (Fig 5) | 36/57 |
| **Data (temporal_classification.tsv)** | `temporal_class == 'early'` = **36** |
| Phase ratio > 0.6 (recomputed) | **36** |

Also corrected adjacent count:

| | Manuscript (old) | Data |
|---|---|---|
| Gradual responders | 8 | **7** |
| Late responders | 12 | 12 ✓ |
| Non-responders | 2 | 2 ✓ |

Note: The percentage "63%" was already correct for 36/57 (= 63.2%); only the raw count was wrong.

**Locations changed:**
- `01_results.md` line 61
- `full_manuscript.md` line 83

---

## Fix #5 — Bloc sizes ✅ FIXED

| Bloc | Manuscript (old) | Figure (Fig 5b) | **Data** |
|---|---|---|---|
| Activation | ~35 genes, modules 1–3 | 34 | **34 genes** (module 1: 33 + 1 unassigned) |
| Repression | ~26 genes, module 4 | 23 | **23 genes** (module 2) |
| Total | 61 (inconsistent!) | 57 | **57** ✓ |

Data: `temporal_classification.tsv`, `direction` column: up=34, down=23.

Also corrected: `nine TCS components (25% of bloc members)` → `nine TCS components (26% of bloc members)` (9/34 = 26.5%).

**Locations changed:**
- `01_results.md` lines 51, 53
- `full_manuscript.md` lines 73, 75, 109

---

## Fix #6 — Expression AUC (baseMean)

| | Value |
|---|---|
| Manuscript | 0.547 |
| Fig 4 (script hardcode) | 0.543 |
| Fig 3 (script hardcode) | 0.547 |
| **cross_validation.tsv: 5-fold CV mean AUC** | **0.5468 → 0.547** |
| In-sample pooled AUC (−bm direction) | 0.5429 → 0.543 |
| "Fig 3: 0.444" (user-reported) | **Cannot reproduce** — not found in any data file or script |

**Verdict: Manuscript (0.547) is CORRECT — it reports the 5-fold CV mean AUC. No change needed.**

The 0.543 in Fig 4's hardcoded label is the in-sample pooled AUC and is a minor figure-text discrepancy. The 0.444 value reported by user cannot be reproduced from current data.

---

## Summary of all changes

| # | Item | Old value | New value | Files changed |
|---|---|---|---|---|
| 2 | Mann-Whitney p | 5.3 × 10⁻²⁸ | **1.7 × 10⁻²⁵** | 01_results.md, full_manuscript.md |
| 4 | Early responders | 39/57 (63%) | **36/57 (63%)** | 01_results.md, full_manuscript.md |
| 4 | Gradual responders | 8 genes | **7 genes** | 01_results.md, full_manuscript.md |
| 5 | Activation bloc size | ~35 genes, modules 1–3 | **34 genes, module 1** | 01_results.md, full_manuscript.md |
| 5 | Activation bloc TCS % | 25% | **26%** | 01_results.md, full_manuscript.md |
| 5 | Repression bloc size | ~26 genes, module 4 | **23 genes, module 2** | 01_results.md, full_manuscript.md |
| 5 | Activation/repression in Discussion | ~35 / ~26 | **34 / 23** | full_manuscript.md |

**No changes:** Shielded count (998 ✓), Distance AUC (0.917 ✓), Expression AUC (0.547 ✓)
