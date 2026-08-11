# Paper Tone-Down Proposal: From "Protection-Zone Collapse → Derepression" to "Timing Correlation"

**Date:** 2026-05-04
**Context:** Reviewer concern that the manuscript's strongest claim — that protection-zone collapse drives derepression of the 57 Exposed TFs — lacks a direct, gene-level scatter showing methylation change versus log2FC. New analysis in `11_epigenome_integration/analysis/57_temporal_dynamics_exposed_TF/` provides this scatter and reveals that the relationship is (a) weakly positive (Spearman *rho* approx. +0.13 to +0.30), (b) opposite in sign to the "protection collapse → derepression" prediction, (c) sensitive to TSS window choice, and (d) not significant when restricted to the canonical GCCGGC + AAGCCCG motifs.
This proposal rewrites the directional/causal language in the Abstract, Results R3/R4, and Discussion as a "timing correlation" framing, and adds an explicit limitation paragraph.

---

## 0. New empirical evidence summary

| Window (bp) | Sites used | rho (Exposed, n=57) | p | rho (Shielded, n=355) | p |
|---|---|---|---|---|---|
| ±200 | GCCGGC + AAGCCCG | +0.243 | 0.069 | -0.014 | 0.79 |
| ±200 | All 4mC + 6mA | +0.210 | 0.12 | -0.011 | 0.83 |
| ±293 | GCCGGC + AAGCCCG | +0.164 | 0.22 | -0.015 | 0.78 |
| ±293 | All 4mC + 6mA | +0.152 | 0.26 | +0.016 | 0.76 |
| ±500 | GCCGGC + AAGCCCG | +0.168 | 0.21 | +0.033 | 0.54 |
| ±500 | All 4mC + 6mA | +0.257 | 0.053 | +0.037 | 0.49 |
| ±1000 | GCCGGC + AAGCCCG | +0.247 | 0.064 | -0.016 | 0.76 |
| ±1000 | All 4mC + 6mA | +0.305 | 0.021 | +0.082 | 0.12 |

Per-motif decomposition at ±500 bp:
- **GCCGGC (4mC):** Exposed *rho* = +0.133, *p* = 0.33. 44/57 Exposed TFs have zero GCCGGC change near TSS.
- **AAGCCCG (6mA):** Exposed *rho* = +0.131, *p* = 0.33. 55/57 Exposed TFs have zero AAGCCCG change near TSS.

Quadrant counts among Exposed TFs (combined sites, ±500 bp):
- Methyl gain + LFC up: 20 (35%) — dominant
- Methyl loss + LFC down: 14 (25%)
- Methyl gain + LFC down: 9 (16%)
- Methyl loss + LFC up: 6 (11%) — predicted quadrant under "protection collapse → derepression"

Three points follow:

1. The correlation, when present, is in the *opposite* direction to the protection-collapse hypothesis. Exposed TFs that gain TSS-proximal methylation between T1 and T2 tend to go up in expression, not down.
2. The correlation is not driven by the canonical R-M motifs. When sites are restricted to GCCGGC + AAGCCCG, *rho* drops to 0.13–0.25 and never reaches *p* < 0.05.
3. The signal is window-sensitive: at the H29-derived 293 bp boundary it is at its weakest (*rho* = 0.15–0.16, *p* > 0.2). It only crosses *p* < 0.05 at ±1000 bp with all 4mC + 6mA sites pooled.

These observations support a "timing co-occurrence, but no directional causal link" framing rather than the current "structural protection channels methylation effects" framing.

---

## 1. Abstract rewrite

**Current passage** (`05_abstract.md` lines 5–7):

> Instead, we discover a four-layer "Gatekeeper" architecture. The restriction-modification systems undergo complete geographic redistribution (T1-T2 Jaccard = 0.000) without transcriptional consequences. A ~2,200 bp methylation-free protection zone at regulatory gene promoters, created by evolutionary counter-selection (~33%) and collective protein occupancy (~67%), shields 94.6% of regulatory genes (998/1,055). The remaining 57 exposed regulators (5.4%), defined by a 293 bp distance boundary (AUC = 0.917), constitute a distributed vegetative-to-developmental switch. These 57 genes are organized into two antagonistic blocs — an activation bloc (TCS, sigma factors, RamR) and a repression bloc (TetR-enriched, *p* = 0.028) — that engage simultaneously at the exponential-to-transition boundary.
>
> DNA methylation does not control transcription directly in *S. coelicolor*. Rather, a structural dichotomy at regulatory gene promoters channels methylation effects through a specific minority of regulators encoding a coordinated developmental transition, revealing a previously unrecognized organizational principle of bacterial regulatory genomes.

**Proposed rewrite:**

> Instead, we describe a four-layer "Gatekeeper" architecture. The restriction-modification systems undergo complete geographic redistribution (T1-T2 Jaccard = 0.000) without detectable transcriptional consequences. A ~2,200 bp methylation-free zone at regulatory gene promoters — produced by evolutionary counter-selection of recognition motifs (~33%) and collective protein occupancy (~67%) — covers 94.6% of regulatory genes (998/1,055). The remaining 57 regulatory genes (5.4%), defined by a 293 bp distance boundary (AUC = 0.917), lack this zone and *co-occur* with a coordinated vegetative-to-developmental transcriptional shift: their expression changes are organized into two antagonistic groups (activation: TCS, sigma factors, RamR; repression: TetR-enriched, *p* = 0.028) that engage simultaneously at the exponential-to-transition boundary.
>
> DNA methylation does not control transcription directly in *S. coelicolor*. Direct gene-level analysis confirms this: TSS-proximal methylation change between T1 and T2 is at most weakly correlated with the corresponding log2 fold change (Spearman *rho* = 0.13–0.30 across windows of ±200–1000 bp; not significant at the H29-derived 293 bp boundary). We therefore propose a structural co-occurrence model in which the same 57 regulators that lack the protection zone also participate in the developmental transition, without claiming a direct causal link from methylation change to expression change. This dichotomy at regulatory gene promoters reveals a previously unrecognized organizational principle of bacterial regulatory genomes.

Key edits:
- "without transcriptional consequences" → "without detectable transcriptional consequences" (retains evidence boundary).
- "constitute a distributed vegetative-to-developmental switch" → "co-occur with a coordinated vegetative-to-developmental transcriptional shift" (removes implied agency of the 57 genes as a switch driven by methylation).
- "channels methylation effects through" → removed; replaced with explicit *rho* range and the explicit statement that no direct causal link is claimed.
- Adds one sentence quantifying the new gene-level evidence.

---

## 2. Results section rewrites

### 2.1 Section "Temporal remodeling of the protection zone" (`01_results.md` line 33)

**Current passage:**

> However, this average trend concealed a striking gene-category-specific divergence: among regulatory genes, the 4mC protection ratio collapsed to 1.046 at T3 — a value exceeding 1.0, indicating that by late development the region immediately upstream of regulatory gene TSSs is, if anything, *enriched* for 4mC relative to the flanking regions (Fig. S16c, d). Non-regulatory genes showed no comparable collapse (ratio = 0.821 at T1 to 0.845 at T3), demonstrating that the protection zone breakdown is specific to the regulatory gene class.

**Proposed rewrite:**

> The category-specific trend was distinct: among regulatory genes, the 4mC protection ratio rose to 1.046 at T3 — a value above 1.0, indicating that by late development the region immediately upstream of regulatory gene TSSs is, on average, no longer depleted (and is marginally *enriched*) for 4mC relative to the flanking regions (Fig. S16c, d). Non-regulatory genes showed no comparable shift (ratio = 0.821 at T1 to 0.845 at T3). The change is therefore specific to regulatory genes, but we describe it as a *reorganization* rather than a collapse: the average shift (0.823 → 1.046) is a population-level trend, and a gene-level scatter of TSS-proximal methylation change against expression log2FC for the 57 Exposed regulators yields only a weak, window-sensitive correlation (Spearman *rho* = 0.13–0.30 across ±200–1000 bp; *rho* = 0.16, *p* = 0.22 at the H29-derived 293 bp boundary, *rho* = 0.17, *p* = 0.21 when restricted to GCCGGC + AAGCCCG sites at ±500 bp) (Fig. SX). The data therefore document a category-specific temporal redistribution of 4mC near regulatory gene TSSs, but do not establish that this redistribution causes the corresponding expression changes.

Key edits:
- "protection ratio collapsed" → "protection ratio rose ... no longer depleted ... marginally enriched" (descriptive, no implied mechanism).
- "the protection zone breakdown is specific" → "The change is therefore specific ... we describe it as a *reorganization* rather than a collapse" (removes "breakdown", replaces with neutral term, and explicitly flags scope).
- Adds one sentence pointing to the new gene-level scatter and its quantitative limits.

### 2.2 Section "Negative results define the scope of the model" — strengthen (`01_results.md` line 77, last paragraph)

The current paragraph already includes the temporal-decoupling observation (*rho* = 0.136, *p* = 0.299), which is the strongest tone-down anchor in the manuscript. We propose extending it with the new gene-level evidence:

**Insertion at end of paragraph (after "...will require experimental approaches..."):**

> A direct gene-level test reinforces this conclusion. Plotting TSS-proximal methylation change (T2 minus T1, sites within TSS ±500 bp) against log2 fold change for the 57 Exposed regulators yields a Spearman *rho* of +0.26 (*p* = 0.053; *n* = 57 with valid expression data). The sign is positive — Exposed TFs that *gain* methylation tend to go up in expression, the opposite of a "protection-zone collapse → derepression" model — and the value is sensitive to choices of TSS window (*rho* = 0.21 at ±200 bp; *rho* = 0.15 at the H29 boundary of ±293 bp; *rho* = 0.31 at ±1000 bp) and to the inclusion of non-canonical methylation sites (when restricted to GCCGGC + AAGCCCG sites the correlation drops to *rho* = 0.17, *p* = 0.21 at ±500 bp). Among the 57 Exposed TFs with computable methylation change, the dominant joint outcome is "methylation gain + LFC up" (*n* = 20) rather than the "methylation loss + LFC up" pattern (*n* = 6) predicted by a simple protection-collapse model. Together with the temporal-decoupling result, these gene-level patterns rule out a direct, directional causal model and are consistent with a co-occurrence relationship in which TSS-proximal methylation reorganization and the developmental transcriptional shift are concurrent features of the T1-to-T2 transition rather than cause and effect.

### 2.3 Section "57 exposed regulators form a distributed developmental switch" — soften the "switch" language (`01_results.md` line 47)

Two minor edits to remove agency:

- "If the 57 exposed regulators are the sole conduit through which methylation interacts with the transcriptional program" → "If the 57 exposed regulators *spatially co-localize* with the methylation events that interact with the transcriptional program"
- "constituting a distributed developmental switch" → "constituting a regulatory subset that co-varies with the developmental transition"

(Both retain the same factual content but remove the directional "channel/switch" framing.)

---

## 3. Discussion: explicit limitation paragraph

Insert as the second paragraph under "Limitations and future directions" (`03_discussion.md` line 49), after the existing "First, all conclusions are based on correlational evidence..." sentence.

**Proposed new paragraph:**

> An additional limitation concerns the directionality of the protection-zone–expression relationship. A direct gene-level scatter of TSS-proximal methylation change (T2 minus T1) versus expression log2 fold change for the 57 Exposed regulators yields only a weak, window-sensitive correlation: Spearman *rho* ranges from +0.15 to +0.31 across ±200–1000 bp windows and is not significant at the H29-derived 293 bp boundary (*rho* = 0.16, *p* = 0.22) or when restricted to canonical GCCGGC and AAGCCCG sites (*rho* = 0.17, *p* = 0.21 at ±500 bp). The sign of the correlation is also incompatible with a naive "protection-zone collapse → derepression" mechanism: Exposed regulators that gain TSS-proximal methylation between T1 and T2 tend to *increase* in expression rather than decrease (20 of 57 Exposed TFs are in the "methylation gain + LFC up" quadrant versus 6 in the "methylation loss + LFC up" quadrant). We therefore frame the relationship between protection-zone reorganization and the developmental transition as a *temporal co-occurrence* rather than a causal chain. The gene-level evidence does not support a model in which loss of TSS-proximal methylation directly triggers transcriptional activation of the 57 Exposed regulators. Establishing the actual directionality, if any, will require methyltransferase knockout experiments (predicted to alter expression of Exposed but not Shielded regulators if the model is correct) and methylation-targeted CRISPR-dCas9 perturbations at individual Exposed TF promoters.

This paragraph explicitly:
- States the new quantitative limit (*rho* range, *p* values, window dependence).
- Notes the sign mismatch with the simple causal model.
- Describes the relationship as "temporal co-occurrence" rather than causation.
- Specifies which experiments would resolve the question.

---

## 4. Optional: Discussion "Reframing bacterial methylation-transcription interactions" — single-sentence edit

`03_discussion.md` line 11 currently reads:

> Methylation's regulatory influence is thus channeled through a narrow bottleneck of 5.4% of regulatory genes, rather than operating genome-wide as previously assumed for many bacterial methylation systems.

**Proposed rewrite:**

> The structural relationship between methylation and the transcriptional program is therefore concentrated at a narrow subset (5.4%) of regulatory genes — the 57 Exposed regulators — rather than operating genome-wide; gene-level analysis within this subset shows only weak correlations between methylation change and expression change (Spearman *rho* = 0.13–0.30 across windows; not significant at the 293 bp boundary or when restricted to canonical motifs), so the term "regulatory influence" should be read as "spatial co-occurrence" rather than as a demonstrated causal mechanism.

---

## 5. Suggested figure / supplement updates

- Add the new scatter `exposed_TF_methyl_change_vs_LFC` as a main-text panel (Fig. 4d or 4e) or as a Supplementary Figure cited from Section 2.1 above. Source files:
  - `11_epigenome_integration/analysis/figures/exposed_TF_methyl_change_vs_LFC.{pdf,svg}`
  - `11_epigenome_integration/analysis/57_temporal_dynamics_exposed_TF/figures/exposed_TF_methyl_change_vs_LFC.{pdf,svg}`
- Add the motif decomposition `exposed_TF_motif_decomposition.{pdf,svg}` as a supplementary panel.
- Add the window sensitivity plot `exposed_TF_window_sensitivity.{pdf,svg}` as a supplementary panel; cite from the new Discussion limitation paragraph.
- Source data table for these figures: `11_epigenome_integration/analysis/57_temporal_dynamics_exposed_TF/tables/exposed_TF_methyl_change_vs_LFC.tsv` and `exposed_TF_window_sensitivity.tsv`.

---

## 6. Net impact of the rewrite

- The Gatekeeper model is preserved as a *structural* model: the protection zone exists, the 57 Exposed regulators exist, and they participate in the developmental transition.
- What is removed is the implicit causal arrow from "loss of protection" to "increase in expression". The data do not support this arrow at the gene level.
- The model becomes harder to attack on the "show the scatter" axis because the manuscript now reports the scatter, its window sensitivity, and its incompatibility with the naive causal model.
- Reviewers asking "is this just spatial co-occurrence?" now have an explicit answer: yes, that is the framing.

End of proposal.
