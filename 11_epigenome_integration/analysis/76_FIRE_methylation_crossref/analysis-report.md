# FIRE × methylation cross-reference — analysis report (2026-06-15)

**Question.** Can the vegetative methylome serve as a Hi-C-free proxy for the 3D
local-interaction signal (FIRE) that Deng et al. (2023, PNAS, PMID 36877856) showed
predicts transcription and integration-site performance in *S. coelicolor*?

**Data.** Deng Dataset S4 (per-5-kb FIRE value + compartment PC, M-phase and L-phase;
1,724 windows) and Table S1 (20 HCR integration loci with FIRE). Our GCCGGC m4C sites
(T1 n=1,289; T2 n=407). Timepoint alignment: Deng M-phase ~14 h ≈ our T1 (12 h);
L-phase ~22 h ≈ our T2 (24 h). Same chromosome coordinates (NC_003888.3 / AL645882.2).

## Key findings

1. **Test A (genome-wide).** Vegetative (T1) GCCGGC density **positively tracks FIRE_M**:
   Spearman ρ = 0.319 (p = 5.6×10⁻⁴²; n = 1,724 bins), surviving core/arm control
   (partial r = 0.179, p = 7.6×10⁻¹⁴). Methylated bins have higher FIRE than unmethylated
   (median 1.13 vs 0.69; MWU p = 2.1×10⁻³⁵).

2. **Sign flip at the switch.** T2 methylation **anti-correlates** with L-phase FIRE
   (ρ = −0.420, p = 1.3×10⁻⁷⁴): the core→arm relocation moves methylation *out* of the
   high-FIRE active compartment. Directly mirrors the reframe's relocation thesis.

3. **Test B (Exposed).** The 62 methylation-marked Exposed regulators sit in significantly
   **higher-FIRE** promoters than Shielded (median 1.27 vs 0.95; MWU p = 8.6×10⁻⁶).

4. **Test C (integration loci — the decisive test).** At the 10 HCR-M loci Deng *selected and
   validated* for FIRE-guided expression, our T1 methylation density **recovers their FIRE
   ranking**: Spearman ρ = 0.661, p = 0.038 (n = 10). Because Deng established FIRE→expression
   at these loci, methylation density nominates the same high-expression integration
   neighbourhoods — from one Nanopore run, without Hi-C.

5. **Compartment cross-check.** 61.6% of T1 GCCGGC sites fall in the 2.3–6.2 Mb core
   (reproduces the manuscript's "62% in Compartment A"); core bins are enriched for
   methylation (OR = 2.73, p = 2.7×10⁻²⁴).

## Interpretation / decision change

The Discussion's "testable prediction" (methylome as Hi-C-free integration-site proxy) is
now an **empirical result**: methylation density tracks FIRE genome-wide and, critically,
ranks Deng's own validated integration loci (ρ = 0.66). Magnitude is moderate, not
deterministic (genome-wide controlled r ≈ 0.18; the HCR fit has outliers, e.g. M7) — so the
honest claim is **pre-screening / search-narrowing**, not a replacement for FIRE mapping.
→ Upgrade Discussion engineering paragraph from prediction to result; add as Supplementary figure.

## Caveats / blockers
- Direct methylation→*gus*/RK-682 expression untestable: Deng's reporter values are in
  bar-chart figures, not tabulated. We use FIRE (their published predictor) as the proxy
  ground truth at the integration loci; the n=10 rank test is significant but small.
- PC eigenvector (sd04) is non-discriminating here (PC>0 for 92% of bins; ρ≈0.03 NS) — we
  rely on FIRE and the explicit 2.3–6.2 Mb core, not PC sign.
- FIRE windows are 55 kb sliding 5 kb (overlapping); binning uses the 5 kb anchor.

## Outputs
- `tables/fire_methyl_bins.tsv`, `tables/regulatory_FIRE_at_TSS.tsv`
- `figures/fig1_genomewide_FIRE_vs_methylation.{pdf,png}`
- `figures/fig2_exposed_vs_shielded_FIRE.{pdf,png}`
- `figures/fig3_HCR_integration_FIRE_vs_methylation.{pdf,png}`
- `scripts/fire_methyl_crossref.py`, `scripts/fire_methyl_figures.py`
