# Spurious Correlation Validation Report

**Date:** 2026-02-03
**Project:** *Streptomyces coelicolor* A3(2) M145

## Executive Summary

### Key Finding: 4mC T2vsT1 Correlation

| Test | Result | Interpretation |
|------|--------|----------------|
| Observed r | 0.1248 | Positive correlation |
| Permutation p | 0.0020 | Significant |
| Partial r retained | 111.2% | Robust |
| 95% Bootstrap CI | [0.046, 0.194] | Excludes 0 |
| Z vs random | 3.19 | Strong |

## Detailed Results

### 1. Permutation Test

| Comparison | Mod | N | Observed r | Perm p |
|------------|-----|---|------------|--------|
| T2vsT1 | 4mC | 626 | 0.1248 | 0.0020 |
| T2vsT1 | 6mA | 621 | -0.0378 | 0.3487 |
| T3vsT1 | 4mC | 542 | -0.0034 | 0.9369 |
| T3vsT1 | 6mA | 662 | -0.0932 | 0.0150 |

### 2. Partial Correlation

| Comparison | Mod | Simple r | Partial r | Retained |
|------------|-----|----------|-----------|----------|
| T2vsT1 | 4mC | 0.1248 | 0.1388 | 111.2% |
| T2vsT1 | 6mA | -0.0378 | -0.0367 | 97.2% |
| T3vsT1 | 4mC | -0.0034 | 0.0044 | 130.6% |
| T3vsT1 | 6mA | -0.0932 | -0.0858 | 92.0% |

### 3. Bootstrap CI

| Comparison | Mod | r | 95% CI | Zero in CI |
|------------|-----|---|--------|------------|
| T2vsT1 | 4mC | 0.1248 | [0.046, 0.194] | **No** |
| T2vsT1 | 6mA | -0.0378 | [-0.114, 0.047] | Yes |
| T3vsT1 | 4mC | -0.0034 | [-0.088, 0.089] | Yes |
| T3vsT1 | 6mA | -0.0932 | [-0.168, -0.012] | **No** |

## Conclusion

**2/4 comparisons show robust correlation**

The **4mC-expression correlation at T2vsT1** passes validation:
- Permutation test significant
- Correlation robust after controlling confounders
- The correlation is **unlikely to be spurious**

## Output Files

| File | Description |
|------|-------------|
| `permutation_test_results.png` | Null distribution vs observed |
| `partial_correlation_comparison.png` | Simple vs partial correlation |
| `negative_controls.png` | Observed vs random |
| `region_specificity.png` | Promoter vs gene body |
