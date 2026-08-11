# A-2: Protection Zone Boundary Sensitivity + T2/T3 Independent Classifier Models

**Date**: 2026-04-20  
**Analysis dir**: `11_epigenome_integration/analysis/69_boundary_sensitivity/`

## Results Summary

### Analysis 1: Boundary Sensitivity (100–500 bp sweep)

Feature: methylation site count within upstream window (all-TP combined), 5-fold stratified CV.

| Window (bp) | CV5 AUC | Exposed captured |
|------------|---------|----------------|
| 100 | 0.388 | 1/11 |
| 150 | 0.536 | 2/11 |
| 200 | 0.599 | 5/11 |
| 250 | 0.763 | 8/11 |
| **293** | **0.808** | **9/11** |
| 350 | 0.792 | 9/11 |
| 400 | 0.794 | 9/11 |
| **450** | **0.834** | **10/11** |
| 500 | 0.826 | 10/11 |

293 bp is near-optimal; 450 bp gives marginal improvement (+0.026) by capturing one more exposed gene.

### Analysis 2: Timepoint-Independent Models

Feature: `-nearest_methyl_distance` per timepoint, 5-fold stratified CV.

| Model | CV5 AUC | Cross-TP AUC |
|-------|---------|-------------|
| All-TP (H29 orig) | 0.911 | — |
| T1 independent | 0.670 | — |
| T2 independent | 0.774 | 0.644 (T1→T2) |
| T3 independent | 0.829 | 0.533 (T1→T3) |

T3 independent outperforms T1→T3 cross-TP by +0.296.

## Key Conclusions

1. 293 bp boundary is a valid choice but not uniquely optimal; the 250–450 bp range all give good performance
2. T3-only methylation (0.829) is the best single-timepoint predictor; all-TP combination (0.911) is best
3. Per-timepoint independent models substantially outperform cross-timepoint prediction at T3

## Output Files

- `tables/A2_boundary_sensitivity_AUC.tsv`
- `tables/A2_timepoint_models_AUC.tsv`
- `figures/A2_boundary_sensitivity.png`
- `figures/A2_timepoint_ROC.png`
