# C6: Exposed(62) vs Shielded(989) comprehensive metric scan (2026-06-29)

Question (reviewer C6): did we comprehensively check which metrics differ, or only |log2FC|?
Source: 52_.../SuppTable_regulatory_gene_classification_n1051.tsv (canonical 62/989). Seed 42, BH-FDR.

## Significant after FDR (BH<0.05)
- position core: OR=3.87 (p_BH=3e-4)
- distance to chromosome centre: Cliff δ=−0.316 (Exposed closer; p_BH=3e-4)
- TF family MerR: OR=4.71 (p_BH=0.022); LysR: OR=3.47 (p_BH=0.035)
- signed LFC_T2vsT1: Cliff δ=−0.193 (p_BH=0.038) — weak repression-direction bias, == the reported geography-controlled r=−0.09 modulatory effect

## NOT significant (FDR)
- baseline expression (baseMean); |LFC| magnitude/variability T2 & T3; DEG fraction T2 & T3; LacI/TetR/Sigma/sensor-kinase families

## Interpretation (permissive, supported & honest)
Exposed regulators differ from Shielded in **chromosomal position (core/centre) and TF family**, and show a **weak signed-direction bias at T2** (the same r=−0.09 effect) — but NOT in expression magnitude, variability, DEG fraction, or baseline. The partition does not select a transcriptionally distinct cohort by magnitude; it selects a positionally/family-distinct cohort with only a weak directional bias = permissive.
