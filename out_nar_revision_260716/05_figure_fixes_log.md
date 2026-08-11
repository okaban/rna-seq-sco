# Figure fixes log — 修正記録（証拠付き）
# All edits to figure SCRIPTS (git-reversible). Manuscript image slots restored to
# original after each proofing render; the author regenerates + commits when approved.

## Figure 1 (script: 01_figure1_landscape.py) — C4, C5, C6 ✅
- **C6** logo titles: `color=COL_4mC/COL_6mA` → `color='black'` (lines ~477, 491).
- **C4** oriC: green `#2E7D32` → grey `#333333` (line + label); reduces palette so colour is reserved for 4mC/6mA data.
- **C4** core boundary: added on-panel annotation "1.50–7.17 Mb" under the (left-shifted) Core label at ~2.9 Mb, clear of the oriC line.
- **C5** Fig1C asterisks: significance markers re-anchored to a per-region common ceiling `max(4mC,6mA,Genome)+2.5` instead of each bar's own top, so the Intergenic 4mC (5.8) marker no longer collides with the 6mA (10.7) value label. Value-label offset tightened +0.8→+0.4.
- Verified visually: figs_review/Figure1_{ORIGINAL,REVISED}.png.

### C4 core-boundary annotation — placement refined (post-QA)
The on-panel core annotation now reads **"core 1.50–7.17 Mb"** (value = ARM_LEFT/1e6 = 1.50,
matching locked canon), placed above the x-axis line under the "Core" label with a white
background box so the axis line no longer strikes through the text.

## Figure 2 (script: 02b_figure2_RM_redistribution.py) — C14 ✅ (substantive)
### C14 — protection-zone band was mis-placed. FIXED with data evidence.
**Evidence (analysis/48_TSS_methylation_gradient/spatial_profile_data.tsv, GCCGGC_4mC, T1):**
Baseline = median density at |dist|>2500 bp = 0.159 (all genes).
| Region | mean density | % of baseline |
|---|---|---|
| −2200 to 0 bp (OLD shaded band) | 0.157 | **99%** (NOT depleted) |
| −2200 to −300 bp (upstream only) | 0.163 | **103%** (above baseline) |
| ±300 bp (TSS-proximal) | 0.122 | **77%** (depleted) |
Per-bin: trough at −100 bp = 64% baseline, +100 bp = 69%; back to baseline (≥90%) by ±300 bp.
Same pattern in regulatory (±300 = 83%) and non-regulatory (76%) subsets.
**Conclusion:** the depletion is a narrow TSS-proximal dip (~±200–300 bp), NOT a broad ~2,200 bp
upstream zone. The old red band covered a region at baseline — the reviewer (C14) is correct.
**Fix:** band redrawn symmetric at ±300 bp, coinciding with the density trough and bracketed by the
293 bp classification cut-off lines. Label changed to "TSS-proximal depletion (±300 bp)".
Verified: figs_review/Figure2_{ORIGINAL,REVISED}.png.

### ⚠️ AUTHOR FLAG — "~2,200 bp protection zone" text claim (Key-Claims R2-2 / H25)
The 2,200 bp figure comes from the SEPARATE regulatory-gene *symmetric* protection-zone analysis
(Claim R2-2: "2,200 bp symmetric protection zone around regulatory TSS", depletion center +300 bp,
reg/non-reg ratio 0.541). That is a DIFFERENT quantity from the all-genes GCCGGC metagene shown in
Fig 2B. The manuscript currently lets Fig 2B's caption/legend carry the "~2,200 bp" number, which
this panel's own data does not support. **Author decision needed:**
  (a) keep Fig 2B as the all-genes metagene with the corrected ±300 bp band (done), AND move the
      "~2,200 bp symmetric zone" claim to the figure/panel that actually shows the regulatory-gene
      symmetric zone (Fig 3 protection-zone panel / SuppFig), or
  (b) replace Fig 2B with the regulatory-gene symmetric-zone metagene if the 2,200 bp claim is meant
      to be the headline here.
My default (a): the ±300 bp all-genes depletion + 293 bp cut-off is the internally consistent story
for Fig 2B; the 2,200 bp regulatory symmetric zone belongs with the regulatory-gene analysis.

### C9 / C10 / C12 / C13 — see 06_figure2_restructure_proposal.md (panel composition = author decision)

---
## Discrepancy flags surfaced during QA (AUTHOR must reconcile) ⚠️

### SuppFig10 AAGCCCG 6mA site count — RESOLVED ✅ (was 447; corrected to 441)
The figure hardcoded **447**; the manuscript text says **441**; the high-confidence census
gives 260. Traced to the canonical publication-figure data table
`02_publication_figures/timepoint_motif_prevalence.csv`, which lists T1 6mA AAGCCCG = **441**
(441/1934 = **22.80%**, exact match to the manuscript's stated occupancy) and 4mC AAGCCCG = 699.
**441 is canonical; the figure's 447 was a transcription error.** Fixed the generator
(fig_modpos_and_redistribution.py: 447->441) and regenerated; the bar now reads 441. The 260
figure is the stricter high-confidence-census subset (different threshold), not the number the
manuscript commits to — no text change needed.

### SuppFig15 undecylprodigiosin (Red) fold-change — RESOLVED from figure data ✅
Manuscript states "~30-fold from T1 to T3". Author confirms only the figure data exist (no source
spreadsheet), so the fold is judged from the reconstructed figure values. Red T1 sits at the
detection floor (points 0.002/0.002/0.004/0.009 mg/L; mean 0.00425), so any T1-denominated ratio
is unstable — T3/T1 (mean) ≈ 1214×, using max(T1) ≈ 573×; neither is ~30×. The only stable,
figure-supported fold is the detectable-range change **T3/T2 ≈ 17×**.
**Recommendation (text-wording change for author):** replace "~30-fold from T1 to T3" with a
floor-aware statement, e.g. "rose from near the detection limit at T1 to 5.2 ± 0.7 mg/L at T3
(T2->T3 ≈ 17-fold; a T1->T3 ratio is not quoted because T1 is at the detection floor)." DCW/Act
panels unaffected (they reproduce the manuscript ANOVA).

### make_growth_merged.py — post-hoc test method
Corrected to match manuscript Methods: DCW/Act now use Tukey HSD (statsmodels pairwise_tukeyhsd);
Red uses Welch pairwise + Holm as a Dunnett-T3 surrogate (statsmodels lacks Dunnett T3; Red's
near-zero T1/T2 variance makes all pairwise verdicts identical either way). Categorical brackets
unchanged (DCW ****/****/ns; Act ns/****/****; Red ****/**/***).
