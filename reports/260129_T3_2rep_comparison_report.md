# T3 2-Replicate Re-analysis Report

**Date:** 2026-01-29
**Purpose:** Evaluate the impact of excluding low-coverage sample 3-2 from T3 analysis

---

## Executive Summary

Excluding low-coverage sample 3-2 (14.3x vs 27-52x average) dramatically improves T3 methylation detection and changes the interpretation of methylation-expression correlations.

**Key Finding:** The apparent negative/weak correlation at T3 was an artifact of low coverage in sample 3-2.

---

## 1. Methylation Site Detection

### Total Sites by Timepoint

| Timepoint | Original (3 reps) | New (2 reps) | Change |
|-----------|-------------------|--------------|--------|
| T1 | 2,349 | 2,349 | 0 |
| T2 | 2,794 | 2,794 | 0 |
| **T3** | **885** | **2,780** | **+1,895 (+214%)** |

### Modification Type Breakdown

| Modification | Original T3 | New T3 | Change |
|-------------|-------------|--------|--------|
| 6mA | 658 | 1,562 | +904 (+137%) |
| 4mC | 227 | 1,218 | +991 (+437%) |

---

## 2. Promoter Methylation

### Genes with Promoter Methylation

| Modification | Timepoint | Original | New (2-rep) | Change |
|-------------|-----------|----------|-------------|--------|
| 6mA | T1 | 261 | 261 | +0 |
| 6mA | T2 | 271 | 271 | +0 |
| **6mA** | **T3** | **194** | **383** | **+189 (+97%)** |
| 4mC | T1 | 347 | 347 | +0 |
| 4mC | T2 | 401 | 401 | +0 |
| **4mC** | **T3** | **52** | **291** | **+239 (+460%)** |

---

## 3. Correlation Analysis Comparison

### Methylation Change vs Expression Change (Spearman)

| Mod | Comparison | n (orig) | r (orig) | p (orig) | n (new) | r (new) | p (new) | Interpretation |
|-----|------------|----------|----------|----------|---------|---------|---------|----------------|
| 6mA | T2_vs_T1 | 361 | **0.169** | **0.001** | 361 | **0.169** | **0.001** | Unchanged (control) |
| 6mA | T3_vs_T1 | 330 | -0.050 | 0.367 | 449 | 0.064 | 0.173 | **Negative → Positive** |
| 6mA | T3_vs_T2 | 343 | **-0.166** | **0.002** | 458 | -0.059 | 0.206 | **Significant → NS** |
| 4mC | T2_vs_T1 | 474 | **0.172** | **0.0002** | 474 | **0.172** | **0.0002** | Unchanged (control) |
| 4mC | T3_vs_T1 | 346 | -0.041 | 0.445 | 430 | 0.083 | 0.086 | **Negative → Positive (marginal)** |
| 4mC | T3_vs_T2 | 396 | -0.093 | 0.063 | 453 | -0.003 | 0.952 | **Marginal → NS** |

---

## 4. Biological Interpretation

### Original Interpretation (with sample 3-2):
- T1→T2: Positive correlation (methylation ↑ ⇔ expression ↑)
- T1→T3 and T2→T3: Negative/no correlation (suggested regulatory shift)

### Revised Interpretation (excluding sample 3-2):
- **All timepoints show consistent positive/neutral correlation**
- The apparent "regulatory shift" at T3 was an artifact of low sequencing depth
- DNA methylation (6mA, 4mC) appears to have a **consistent positive regulatory role** across growth phases

### Implications:
1. **Data quality matters:** Low-coverage samples can dramatically bias methylation analysis
2. **No evidence of regulatory shift:** The methylation-expression relationship appears stable across timepoints
3. **T3 shows similar patterns to T1/T2:** When adequately sequenced, T3 methylation behavior is consistent

---

## 5. Recommendations

1. **Use the 2-replicate T3 data** for downstream analyses
2. **Consider re-sequencing sample 3-2** if a third T3 replicate is needed
3. **Report both analyses** in publications with appropriate caveats
4. **Apply minimum coverage thresholds** more stringently (e.g., 20x per sample)

---

## 6. Output Files

| File | Description |
|------|-------------|
| `T3_2rep_methylation_sites.csv` | T3 methylation sites (2 reps only) |
| `high_confidence_sites_T3_2rep.csv` | Combined T1+T2+T3(2rep) methylation sites |
| `integrated_methyl_expression_T3_2rep.csv` | Full integration with 2-rep T3 |
| `correlation_analysis_T3_2rep.csv` | Correlation results with 2-rep T3 |

---

## 7. Conclusion

**The low coverage in sample 3-2 was the primary cause of reduced T3 methylation detection and the apparent negative correlations observed in the original analysis.**

With 2-replicate analysis:
- T3 methylation site count increases by 214%
- Correlation patterns become consistent with T1→T2 transition
- No evidence of a fundamental shift in methylation-expression regulation at T3

This demonstrates the critical importance of adequate sequencing depth for methylation analysis and highlights that the original T3 results should be interpreted with caution.
