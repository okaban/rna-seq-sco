# 90_per_timepoint_census_audit — per-timepoint GCCGGC/AAGCCCG site census: canonical file vs first-appearance tables

Generated 2026-09-21 by the EPI_ANALYST session (BLOCKER-0 follow-up). Facts only; no manuscript edits.

## Inputs
- Canonical per-timepoint sites: `01_integration/high_confidence_sites_weighted.csv` (12,429 rows; n_reps = 3 for all rows; weighted_mod_freq ≥ 50; coverage ≥ 30). `position` is **0-based** (ref[position] is the modified base for all rows: 6,080/6,080 4mC = C; 6,343/6,349 6mA = A).
- Reference: `/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa` (8,667,507 bp).
- First-appearance chain (do NOT use for per-timepoint counts): `07_motif_analysis/methylation_site_sequences.csv` (5,941 rows = union of unique positions) → `23_expanded_motif_search/01_expanded_motif_discovery.py` L139 `drop_duplicates` → `4mC_final_census.csv` → `37_defense_island_GCCGGC/tables/GCCGGC_sites_by_timepoint.tsv` (1,289/407/21).
- Motif mapping: '+' GCCGGC internal C = ref[p-2:p+4]=='GCCGGC'; '-' = ref[p-3:p+3]. AAGCCCG: 6mA at A0/A1, 4mC at C4 (reverse-complement CGGGCTT for '-' strand). All 4mC sites map to GCCGGC or AAGCCCG (0 'other').
- Core = 1.5–7.17 Mb; Compartment A = 2.3–6.2 Mb.

## Tables
| file | content |
|---|---|
| `tables/per_timepoint_census.tsv` | per mod type × site class (all/GCCGGC/AAGCCCG/other): n, core n/frac, median freq per T1/T2/T3; pairwise shared counts + Jaccard; union / all-three |
| `tables/first_appearance_decomposition.tsv` | new-at-T2, new-at-T3, persisting, lost, only-in-one, all-three sets with core fractions — shows how 407 (core 0.182) and 21 (core 0.381) arise |
| `tables/core_arm_fraction_permutation.tsv` | core fraction per timepoint (canonical and dedup) vs 10,000-draw null over all 32,118 genomic GCCGGC C positions (genomic core fraction 0.631); Compartment-A fraction |
| `tables/protection_zone_ratio_by_timepoint.tsv` | GCCGGC 4mC density in TSS −300..+50 (strand-aware) / far-upstream −5000..−2000, per timepoint, for Jeong2016 TSS (3,570), regulatory n=1,051, Exposed 62, Shielded 989; plus the 64_ definition (−500..0 / −2000..−500) applied to GCCGGC only |
| `tables/exposed62_promoter_methylation_by_timepoint.tsv` | Exposed-62 (and Shielded-989) nearest GCCGGC 4mC distance per timepoint: n within 293 bp / 1 kb / 2 kb, median; same nd() definition as `52_/scripts/Exposed_dynamicA_define_and_verify.py` |
| `tables/suppfig8c_dmeth_vs_dfire_recompute.tsv` | ΔT1→T2 GCCGGC count per 5-kb window vs ΔFIRE (M→L), Spearman + circular-shift p, canonical vs dedup |
| `tables/aagcccg_pooled_core_fraction_canonical.tsv` | AAGCCCG 4mC+6mA pooled core fraction per timepoint (Supp Fig 11 c/d replacement) |
| `tables/manuscript_affected_numbers.tsv` | 44 EN-manuscript locations whose numbers depend on per-timepoint site sets; status dedup_artifact (33) / canonical_ok (8) / not_recomputed (2) / untraceable (1) |

## Headline numbers (canonical)
GCCGGC 4mC: 1,289 / 1,595 / 1,073; core 0.830 / 0.660 / 0.687 (null 0.631; all three core-enriched, perm p 0.0002 / 0.016 / 0.0002); Jaccard T1–T2 0.700 (1,188 shared), T2–T3 0.632, T1–T3 0.591; all-three 858.
First-appearance decomposition: new at T2 = 407 (core 0.182), new at T3 = 21 (core 0.381), lost after T1 = 101 (core 0.911).
Exposed-62 within 293 bp: 62 / 56 / 45; median nearest 172 / 181 / 184 bp (dedup: 62 / 0 / 1; 172 / 57,884 / 196,230 bp).
Compartment A (2.3–6.2 Mb): 61.5% / 46.8% / 51.1% (dedup 61.5% / 4.2% / 19.0%).
BGC (29 antiSMASH regions) GCCGGC 4mC: 140 / 214 / 118 (dedup 140 / 80 / 4 — reproduces the manuscript exactly).
Supp Fig 8c: ρ = 0.397, circular-shift p = 0.030 (dedup 0.428 / 0.027).
Dynamic sites (≥10 pp, HC at all three timepoints): GCCGGC 244/858 (A5b: 243/855); AAGCCCG 6mA 97/309 (A5b: 47/147 — A5 treats `position` as 1-based; the 147 set = 119 A0 + 27 A1 + 1 non-motif).

## Figure regeneration 2026-09-21 (EPI_ANALYST) — see `tables/figure_regeneration_260921.tsv`
Six manuscript PNGs were regenerated from the canonical file through the shared loader
`canonical_sites.py` (motif assignment identical to `per_timepoint_census.py`); each
producing script now imports it instead of the 37_/23_/07_ first-appearance tables.
Previous PNGs are in `Writing/fig_images/archive/<name>_pre260921.png`.
- `Figure2_merged.png` a/b: 1,289/1,595/1,073, 83/66/69 % core; panel c schematic redrawn
  data-driven (T1 83 % / T2 66 % core, Jaccard 0.70) — message changed, needs author sign-off.
- `Figure5.png` a: 62/56/45 within 293 bp (medians 172/181/184 bp).
- `SuppFigure20.png` c: GCCGGC-Exposed 100/90.3/72.6 % only. The former 'AAGCCCG-prom
  9.1/4.5/0 %' series was `near_frac(AAGp)` on the **GCCGGC** distance columns (fraction of
  the 22 AAGCCCG-promoter regulators with a GCCGGC site ≤293 bp of the TSS), not AAGCCCG
  retention; it is not drawn. Candidate replacements: `52_/tables/S20c_aagcccg_series_candidates.tsv`
  (canonical TSS-based AAGCCCG-6mA class n = 21: 100/85.7/81.0 %). The 22-gene class itself
  (panels a/b) comes from the 23_ 6mA census, whose T1 AAGCCCG set is 260 of the 418
  canonical sites, with distance to gene start/end (gene body = 0) rather than to the TSS.
- `SuppFigure8.png` c: ρ = 0.397 (circular-shift p 0.030 from `tables/suppfig8c_*`).
- `SuppFigure11.png`: GCCGGC 83/66/69 %; AAGCCCG pooled 82/70/73 % (1,116/1,302/1,017).
- `SuppFigure14.png`: n 1,289/1,595/1,073; 698/851/574; 418/451/443 (`79_/tables/C4_region_composition.tsv`).
- `SuppFigure17.png`: not regenerated (T1-only); legend premise only.

## Figure PNGs embedding first-appearance numbers (regenerated 2026-09-21; see above)
- `Writing/fig_images/Figure2_merged.png` (main Fig 2a,b) ← `15_paper_figures/scripts/31_merged_figure2_geography.py` → `02b_figure2_RM_redistribution.py` ← 37_ table
- `Writing/fig_images/Figure5.png` (main Fig 4a, 62/0/1) ← `52_shielded_exposed_boundary/scripts/F4_F5_reframe_figures.py` ← 37_ table
- `Writing/fig_images/SuppFigure20.png` (S20c 100/0/1.6%; 9.1/4.5/0%) ← same script ← 37_ table + `36_AAGCCCG_distribution/tables/AAGCCCG_site_gene_mapping.tsv` (23_ 6mA census)
- `Writing/fig_images/SuppFigure8.png` (S8c ρ 0.43) ← `76_FIRE_methylation_crossref/scripts/make_suppfig8_assembled.py` / `fire_methyl_crossref.py` ← 37_ table (m_n_T2)
- `Writing/fig_images/SuppFigure11.png` (83/18/38; 82/38/61) ← `77_reviewer_figures/fig_modpos_and_redistribution.py` (Figure2_redistribution_4panel) ← 37_ table + 07_ file
- `Writing/fig_images/SuppFigure14.png` (1,289→407→21; T3 n 15–38) ← `79_C4_region_by_timepoint/scripts/c4_region_by_timepoint.py` ← 23_ 4mC/6mA_final_census.csv
- `Writing/fig_images/SuppFigure17.png` (coverage control; T2 arm-localisation premise) ← `78_reviewer_robustness/scripts/C1_copynumber_figure.py` (T1 sites from 37_ = canonical T1; coverage from pileup) — figure data unaffected, legend premise affected
- NOT affected: `Figure1.png` panel b (`01_figure1_landscape.py`, true per-timepoint counts 1,987/2,446/1,647); `SuppFigure5.png` (`70_/A5b_dynamics_10pp_SuppFig5.py`, canonical file; see 1-bp caveat)
- Also dedup-derived downstream (not manuscript figures): `57_temporal_dynamics_exposed_TF/scripts/v_defense_MTase_expression_dynamics.py` (L255 MTase concordance ranking), `47_BGC_methylation_geographic_test` (T1 only → unaffected), `52_/tables/SuppTable1_Exposed62_identity.tsv` columns nearest_GCCGGC_T2/T3.

## Caveats found in passing (not dedup-related)
1. `64_timepoint_TSS` protection ratio (0.839/0.844/1.046) is computed on ALL 4mC (GCCGGC + AAGCCCG-4mC), window −500..0 vs −2000..−500 in 200-bp bins; the Methods text (L160) says GCCGGC with a −5 kb..−2 kb background.
2. `70_methylation_dynamics/scripts/A5_methylation_dynamics.py` L139 `pos0 = position − 1` assumes 1-based coordinates; the canonical file is 0-based. GCCGGC counts are robust (855 ⊂ 858); AAGCCCG 6mA denominator is 147 instead of 309.
3. Denominator wording: 'sufficient coverage at all timepoints' (L47) actually means 'called high-confidence (≥50%) at all three timepoints'.
4. Canonical AAGCCCG 6mA HC site counts are flat across development (418/451/443) while the canon occupancy series (88_) declines 27.4→6.8→3.7%; the two are different metrics (instance-level occupancy vs HC site count) — not reconciled here.
