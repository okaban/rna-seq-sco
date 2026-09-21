# 94 — main-figure restructure (D8–D10), 2026-09-21

Author-approved figure reorganisation after the handling-editor review (paper-narrative).
All panels are computed from the CANONICAL per-timepoint file
`01_integration/high_confidence_sites_weighted.csv` and from the modkit pileups
`/Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup/*_pileup.bed`
(0-based `start`; GCCGGC modified C at motif offset 2 on +, reference G at offset 3 on −).

## What each new panel shows, and where the number comes from

| panel | quantity | table |
|---|---|---|
| Fig 1a | per-site 4mC and 5mC frequency at the 1,289 T1 GCCGGC sites, three replicates pooled, coverage >= 10 (4mC median 86.7%, min 44.1; 5mC median 0, max 1.96) | `tables/canonicalT1_GCCGGC_4mC_5mC_pooled.tsv`, `tables/reassignment_summary_T1.tsv` |
| Fig 1b | core fraction 0.830/0.660/0.687 vs motif-matched null 0.631 (perm p 2e-4 / 0.016 / 2e-4) | `90_per_timepoint_census_audit/tables/core_arm_fraction_permutation.tsv` |
| Fig 1c / Fig 4c | occupancy (+-2 kb, T1) vs log2FC(T2/T1), n = 1,019: raw Spearman -0.127, region-controlled -0.090, core-only -0.123, arm-only -0.015 (2,000-replicate bootstrap CIs) | `tables/occupancy_vs_lfc.tsv`, `tables/occ_lfc_stats.json`, `tables/fig4_correlation_forest.tsv` |
| Fig 2b | presence-pattern counts: 858 all three, 330 T1+T2, 232 T2-only (arm-skewed), 1,717 distinct | `tables/gccggc_persistence_patterns.tsv` |
| Fig 2c | Exposed-promoter retention 62/56/45 | `90_.../tables/exposed62_promoter_methylation_by_timepoint.tsv` |
| Fig 4a/4b | 27/62 vs 435/989 with |log2FC| >= 1; Cliff's delta -0.010, TOST rejects |delta| > 0.33 | `52_shielded_exposed_boundary/tables/all_genes_features_unified_n57.tsv` (recomputed here) |

Note on Fig 4a denominators: the main text uses all 62 / all 989 as denominators (27/62 = 43.5%,
435/989 = 44.0%); among genes WITH a fold-change the fractions are 27/59 = 46% and 435/960 = 45%.
The figure uses the main-text convention and the legend states it.

Old main figures are archived at `Writing/fig_images/archive/*_pre260921_restructure.png`;
old Figure 1 is now Supplementary Figure 12 (the previously unused number).

Deferred by the author: D11 (merge Fig 5/6), D12 (dip-test provenance of the 293 bp cutoff;
R-M-presence-stratified O/E), D13 (renumber supplementary figures). Dropping the T3 tracks was
declined — story A's stability claim spans T3.
