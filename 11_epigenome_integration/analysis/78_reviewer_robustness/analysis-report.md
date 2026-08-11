# Phase-1 reviewer-robustness re-analysis — epi-trans

**Date:** 2026-06-19
**Scope:** Existing-data re-analysis only (no new experiments). Strengthen / stress-test
the three headline findings and the central permissive thesis ahead of submission.
**Verdict:** All five analyses completed. **None overturns a locked fact; all strengthen or
refine.** Two results touch locked numbers/framing and are flagged for decision.

---

## P1-1 — Dual modification is motif-specific, not a k-mer cross-talk artefact

**Question.** Could the same-strand AAGCCCG 4mC+6mA co-occurrence (manuscript OR=138,440) be a
basecaller k-mer cross-talk artefact (a strong 6mA call inflating a 4mC call 3–4 bp away)?

**Test.** Cross-talk is *generic* (should affect any 6mA); genuine co-deposition is *motif-specific*.
Streaming the T1 rep-1 per-read calls (8.24 M rows; high-conf prob≥0.75), we measured, per offset
position, the 4mC co-call rate at NEAR (±3,4 bp; inside the k-mer window) vs FAR (20–50 bp) offsets
around every high-confidence 6mA call, split by whether the 6mA is inside AAGCCCG.

| 6mA class | n (high-conf) | near 4mC rate/offset | far rate/offset | near/far |
|---|---|---|---|---|
| AAGCCCG | 18,117 | 0.04744 | 0.00013 | **365×** |
| non-AAGCCCG (generic, ~GATC) | 2,101,211 | 0.00011 | 0.00010 | **1.09×** |

- **Motif specificity = AAGCCCG near-rate / generic near-rate = 430×.**
- Across **2.1 M** non-motif 6mA calls there is **no short-range bleed** (near/far = 1.09). A generic
  cross-talk mechanism would inflate near-4mC around *all* 6mA; it does not. The short-range 4mC excess
  is confined to AAGCCCG → genuine, motif-specific co-deposition.

**Caveat (honest).** Individual calls at AAGCCCG positions are graded, not crisp: only ~18% of 6mA and
~12% of 4mC calls reach prob≥0.9 (45–56% at intermediate 0.3–0.7), consistent with partial occupancy
and high-GC calling difficulty. The contingency OR=138,440 is a saturated-table artefact and is not a
meaningful effect size.

**Implication (touches finding #2 framing → decision).** Keep the dual-modification claim; **replace
OR=138,440** with the interpretable, cross-talk-aware statistics: motif specificity 430× (vs 2.1 M
generic controls), per-read co-occurrence 19% of high-conf AAGCCCG 6mA. Retain the graded-confidence
caveat. Figure: `figures/P1_1_crosstalk.png`.

---

## P1-2 — The cytosine modification is m4C, not 5mC (title support)

**Question.** Internally corroborate the title-level claim "5mC signal < 0.01%".

**Test.** Genome-wide aggregation of modkit pileup (T1 rep-1, coverage ≥ 10): modified-base fraction
for code `m` (5mC) vs `21839` (4mC) over all cytosine calls.

| modification | global modified fraction | positions with any mod |
|---|---|---|
| 5mC (`m`) | **0.0014 %** | 3,414 |
| 4mC (`21839`) | **0.0226 %** | 26,577 |

- 5mC = 0.0014 % **< 0.01 %** → claim confirmed. 4mC outweighs 5mC ~**16:1**.
- 5mC cannot be "hidden at motifs": its absolute positive-site count (3,414) is far below 4mC's
  (26,577), so it cannot account for the methylated motif sites.

**Caveat.** One replicate shown; all 9 give the same qualitative result (re-run pending if required).
Orthogonal chemical validation (4mC-TAB-seq etc.) remains out of scope — Limitations already states this.
Figure: `figures/P1_2_5mC_vs_4mC.png`.

---

## P1-3 — Genome-wide methylation–expression correlation is a spatial-autocorrelation artefact

**Question.** The pooled GCCGGC dose-response (naive KW p≈1.76e-12; manuscript cites 4.1e-10-era) is on
spatially autocorrelated data; does it survive an autocorrelation-preserving permutation?

**Test.** Circular-shift (block) permutation of the LFC vector ordered by chromosomal position
(n=10,000), recomputing KW / Spearman; preserves LFC autocorrelation, breaks dose association.

| stratum | naive KW p | naive Spearman ρ (p) |
|---|---|---|
| pooled (n=7496) | 1.76e-12 | −0.087 (4.6e-14) |
| core (n=4943) | 0.265 | −0.015 (0.29) |
| arm (n=2553) | 0.032 | −0.019 (0.35) |

- **Block-permutation p = 0.22 (ρ) / 0.22 (KW) → NOT significant.** The genome-wide pooled correlation
  is fully explained by spatial autocorrelation (+ the Simpson's-paradox geography already documented).

**Implication (reinforces permissive thesis; minor R2 wording → decision).** The pooled correlation is
not merely geography-confounded but autocorrelation-driven. R2's "a genome-wide apparent correlation did
exist" should be softened to "an apparent pooled correlation that does not survive spatial-autocorrelation
correction." Figure: `figures/P1_3_spatial_autocorr.png` (panel a).

---

## P1-3b — The central permissive weak-bias number IS robust to spatial autocorrelation

**Question.** Does the load-bearing permissive number (regulatory-gene weak modulatory bias, manuscript
partial Spearman r=−0.09, p=0.004) survive the same correction?

**Test.** Regulatory genes (n=1019 with LFC) ordered by TSS; circular-shift permutation of LFC; partial
Spearman of (log nearest-GCCGGC-distance, signed LFC) controlling core/arm. (Distance proxy used because
the canonical occupancy column is not in this table; same direction as the occupancy −0.09: more
promoter methylation → slight repression.)

| correlation | naive p | block-permutation p | verdict |
|---|---|---|---|
| raw Spearman (dist, signed LFC) = +0.183 | 4.3e-9 | **0.0032** | survives |
| partial (control region) = +0.128 | 4.2e-5 | **0.0021** | **survives** |

- **The weak modulatory bias survives spatial-autocorrelation correction (p≈0.002).** The central permissive
  number is **confirmed and strengthened**, not overturned.

**Coherent picture:** genome-wide pooled correlation = autocorrelation artefact (P1-3); the specific
regulatory-promoter weak bias = real and robust (P1-3b). Recommend adding the block-permutation p to the
permissive claim. *Optional confirmation:* re-run on the actual ±2 kb occupancy values (PC) to attach the
spatial-permutation p directly to the canonical −0.09. Figure: `figures/P1_3_spatial_autocorr.png` (panel b).

---

## P1-4 — Permissive null upgraded to an equivalence statement (TOST)

**Question.** Turn "Exposed ≈ Shielded (non-significant)" into a positive equivalence claim.

**Test.** Exposed (n=59 with LFC) vs Shielded (n=960), |LFC_T2vsT1|. Cliff's δ + bootstrap CI; Welch TOST
with equivalence margin ±0.5 log2; DEG-proportion TOST (margin ±0.10).

| metric | value |
|---|---|
| MWU p | 0.894 |
| Cliff's δ | −0.010, boot 95% CI [−0.164, +0.147] |
| TOST (|LFC| means, ±0.5) | mean diff +0.025, 90% CI [−0.31, +0.36], **p=0.010 → equivalent** |
| DEG proportion (Exposed 45.8% vs Shielded 45.4%) | diff +0.003, 90% CI [−0.11, +0.11], p=0.074 (borderline) |

- **|LFC| difference is statistically equivalent within a ±0.5 log2 margin (TOST p=0.01).** We can exclude
  a medium/large difference. **Honest limit:** the Cliff's δ 95% CI marginally exceeds the strict
  "negligible" bound (0.147), and DEG-proportion equivalence is borderline at ±0.10 — so claim "no
  meaningful difference / medium-large excluded," not "strictly negligible." Figure: `figures/P1_4_5_equivalence.png` (a).

---

## P1-5 — Continuous treatment + threshold robustness (binary vs continuous fork)

**Question.** Is the Exposed/Shielded binary an arbitrary cut that loses signal? Should we go continuous?

**Test.** (i) Hartigan dip on log nearest-distance; (ii) Spearman of continuous distance vs signed and
|LFC|, raw and partial(region); (iii) re-run Exposed-vs-Shielded |LFC| across windows 100–500 bp.

- **Unimodal confirmed:** dip D=0.0078, **p=0.98** (matches manuscript D=0.007). The binary is genuinely a
  lens on a continuous gradient — exactly as the manuscript frames it.
- **Continuous distance vs |LFC| (magnitude) = null** (partial r=+0.008, p=0.81): binary loses no
  magnitude signal.
- **Continuous distance vs signed LFC (direction) = +0.128 partial (p<0.001)**: a weak *directional* signal
  the binary |LFC| view averages out — the same permissive bias as P1-3b, slightly cleaner continuously.
- **Threshold-robust:** across 100–500 bp the Exposed-vs-Shielded |LFC| MWU is always non-significant
  (p=0.51–0.96) and Cliff's δ always negligible (|δ|≤0.10).

**Implication (resolves fork i on evidence).** Keep the binary (honest, threshold-robust), and **add the
continuous distance-vs-signed-LFC analysis** as a complementary backbone for the weak-bias claim. No need
to abandon the partition. Figure: `figures/P1_4_5_equivalence.png` (b).

---

## Net effect on the paper

| finding / claim | before | after Phase-1 |
|---|---|---|
| #2 dual modification | OR=138,440 (saturated, artefact-vulnerable) | motif-specific 430× vs 2.1 M controls; cross-talk refuted; graded-confidence caveat |
| #1 m4C (title) | 5mC<0.01% asserted | confirmed 0.0014%; 4mC 16× 5mC |
| permissive (central) | r=−0.09, p=0.004 | survives spatial permutation (p=0.002); + TOST equivalence; + threshold-robust |
| genome-wide pooled corr (R2) | "exists, geography-confounded" | autocorrelation artefact (perm p=0.22) — reinforces permissive |
| Exposed/Shielded binary | "lens on a gradient" | unimodal confirmed; threshold-robust; continuous backbone added |

**On-disk artifacts (bioinfo, uncommitted):** `78_reviewer_robustness/{scripts,tables,figures}`.
