# 91_comod_stratified_and_hemi — mock-referee recomputations (EPI_ANALYST, 2026-09-20)

Facts only; no manuscript edits. env epi-rna.

## REFA-02 AAGCCCG 4mC(C4)/6mA(A0/A1) coupling, T1
- `scan_perread_C4_keyed.py` → `tables/perread_pairs_T1_C4.tsv.gz` (200,312 instance×read pairs, 1,333 instances, 119,999 reads; identical unit to 79_/comod_full_denominator_C4.py, plus replicate/instance/read keys).
- `stratified_stats_C4.py` → `tables/REFA02b_perread_OR_stratified_T1_C4.tsv`, `tables/REFA02b_per_instance_OR_margins_ge5.tsv`.
  - pooled OR reproduces 79_: 4.84 / 3.17 / 2.40 at p≥0.5 / 0.75 / 0.9.
  - per-replicate 1-1/1-2/1-3: 4.82/4.84/4.85 (0.5); 3.11/3.21/3.19 (0.75); 2.28/2.51/2.42 (0.9).
  - **Mantel–Haenszel OR stratified by instance: 1.07 [1.04–1.10] (0.5); 1.02 [0.99–1.05] (0.75); 0.97 [0.92–1.01] (0.9).**
  - instances with every 2×2 margin ≥5: 647/642/571; individual OR>1 in 53.6% / 48.4% / 43.3%.
  - within-read cross-instance null (4mC at instance i vs 6mA at a different instance j on the same read; 44,857 reads covering ≥2 instances): OR 4.69 / 3.17 / 2.49 — indistinguishable from the same-instance OR on the same reads (4.75 / 3.06 / 2.30).
  - Conclusion (data): the per-read OR is a between-instance/between-read effect (instances and reads with high 4mC also have high 6mA); within an instance, whether a read carries 4mC at C4 is independent of whether it carries 6mA at A0/A1. Per-read coupling does NOT survive stratification.
- `tables/REFA02a_site_level_within_motif_fisher_T1.tsv`: instance level, N=1,334: 698 4mC(C4) × 366 6mA(A0|A1) instances, 194 both, expected 191.5; Fisher OR 1.04 [0.82–1.32], p=0.81. (418 in 79_ summary is the A0+A1 SITE count, not instance count.)

## REFA-06 GCCGGC strand symmetry (`gccggc_strand_symmetry.py`, modkit pileups, 4mC code 21839, cov≥10 both strands; 4mC at offset 2 (+) / 3 (−) of GCCGGC)
`tables/REFA06_GCCGGC_strand_symmetry_pileup.tsv` (8 pileups). Among methylated instances (≥50% on ≥1 strand): hemi 98.3–100%, full 0–1.7%; per-strand frequency at methylated strands median 85–91%, mean 82–88% (fraction ≥99%: 9–35%). Pearson r(+,−) among methylated −0.92 to −0.97; across all instances ≈ −0.02. Replicate 1-3 reproduces every L146 number (1,493 methylated; 98.46% hemi / 1.54% full = 23/1,493; r = −0.924; mean +/− 4.49/4.00%) → Supp Fig 13 source = 78_reviewer_robustness/scripts/G3G5_background_and_symmetry.py on 1-3_pileup.bed. No code in the tree produces "99.8–100% per-strand frequency" or "full-methylation 0%" (0% only in 3-2).

## GATE-03 (`tables/GATE03_coordinated_distance_contrast_n1051.tsv`)
50_ used 29_genomewide_TF_screen all_regulatory_genes (n=1,055 at the time) = 1,051 canonical + 4 (SC_RS02470, SC_RS06460, SC_RS34150, SC_RS41025; all non-coordinated). Its `nearest_methyl_dist` is the distance to the nearest high-confidence site of ANY mod type/motif pooled over T1–T3 (01_integration/high_confidence_sites_weighted.csv, no filter) — not "nearest GCCGGC 4mC site". On n=1,051: any-site metric 113.5 vs 762 bp, p=5.2e-28 (62 vs 989); canonical nearest_GCCGGC_T1 (52_ SuppTable n1051): 2,798 vs 2,677 bp, p=0.65. The 62 "co-varying/coordinated" genes overlap the 62 Exposed in only 6 genes. Coordinated selection (29_) uses gene-assigned site counts changing between timepoints + DESeq2 |LFC|≥1, padj<0.05; no TSS-distance term, but the metric it is contrasted on (any-site distance, T1–T3 pooled) is the same site set used to select them.

## FIG-03 Supp Fig 18c
Script 78_reviewer_robustness/scripts/P1_4_5_equivalence_continuous.py → tables/P1_4_equivalence.tsv; figure P1_figures.py. Two tests coexist: Welch TOST on |LFC| means, margin ±0.5 log2 (T2vsT1: diff +0.025, 90% CI [−0.305, +0.356], p=0.0097; T3vsT1: p=0.0155); Cliff's δ = −0.010, bootstrap 95% CI [−0.164, +0.147], judged against ±0.147 band → NOT within band (delta_equiv=False). The plotted band is ±0.147 (Cliff's δ scale); the TOST p belongs to the ±0.5 log2 margin. Universe n = 59 vs 960 (genes with defined LFC).
