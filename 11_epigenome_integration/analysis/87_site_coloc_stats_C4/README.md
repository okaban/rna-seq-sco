# 87 — site-level AAGCCCG co-localisation statistics, corrected offsets (ledger B01b, 2026-09-13)

`site_coloc_stats_C4.py` -> `tables/site_coloc_2x2_C4.tsv`, `tables/AAGCCCG_occupancy_T1_C4.tsv`

Replaces 68_or_permutation (OR 138,440; 244 in / 1,090 not; universe 1,238,215 hard-coded). Same 2x2 /
Fisher / Haldane / log-OR CI / hypergeometric-permutation design, but with an explicit universe: every
same-strand (A, C) position pair with the A 3 or 4 bp 5' of the C (the A1–C4 / A0–C4 geometry), enumerated
from the reference on both strands: 1,696,280 candidate pairs, 2,668 inside AAGCCCG (2 per instance).

| scope | in: co-mod / not | out: co-mod / not | OR (Haldane) [95% CI] | expected in-motif co-mod under H0 | fold | perm p (10,000, seed 42) |
|---|---|---|---|---|---|---|
| pooled T1–T3 (dedup) | 406 / 2,262 | 0 / 1,693,612 | 6.1 × 10⁵ [3.8 × 10⁴, 9.7 × 10⁶] | 0.64 | 636 | < 1 × 10⁻⁴ (null max 5) |
| T1 only | 214 / 2,454 | 0 / 1,693,612 | 3.0 × 10⁵ [1.8 × 10⁴, 4.7 × 10⁶] | 0.34 | 636 | < 1 × 10⁻⁴ (null max 4) |

Pairs by spacing: pooled A0–C4 (4 bp) 317, A1–C4 (3 bp) 89; T1 181 / 33. Instances with ≥ 1 pair: 345 (pooled), 194 (T1).
Notes: (i) out-of-motif co-modified pairs = 0, so the OR is saturated (Haldane-corrected) and uninformative as an
effect size — as the manuscript already states for 1.4 × 10⁵. (ii) With c = 0 the fold over expectation equals
n_candidate / n_in-motif = 636 regardless of scope; it measures how concentrated the universe is, not an effect
size. (iii) The '910-fold' at EN L150 is a per-read read-identity permutation from the earlier (offset-3/5) per-read
analysis; the corrected per-read full-denominator scan (79_, C4) gives 25,125 co-modified instance × read pairs vs
12,771 expected under independence at ≥ 0.5 (fold 1.97; OR 4.84 [4.73, 4.95]).

6mA occupancy at AAGCCCG, T1 (corrected): 418 canonical 6mA sites over 2,668 A0/A1 positions (15.7%); 366 of 1,334
instances (27.4%) carry 6mA at A0 or A1; 698 of 1,334 C4 positions (52.3%) are canonical 4mC. The manuscript's
'447 6mA sites at T1, 31.6% of genomic AAGCCCG positions' has no traceable denominator (447/0.316 = 1,415).

## Proposed replacement for EN L35 (first two sentences + OR parenthetical) — for the editor, not applied
The AAGCCCG motif also carried N6-methyladenine (6mA) at its two adenines (A₀ and A₁): at T1, 418 of the 2,668
AAGCCCG adenine positions (15.7%) were canonical 6mA sites and 366 of the 1,334 motif instances (27.4%) carried 6mA
at A₀ or A₁, whereas the 4mC call lay at the central cytosine (C₄; 698 of 1,334 instances, 52.3%). At the site
level, canonical 4mC and 6mA calls co-occurred on the same strand of the same motif instance far more often than
expected by chance: pooling canonical sites across timepoints, 406 same-instance pairs carried both marks (A₀–C₄,
4 bp spacing, 317 pairs; A₁–C₄, 3 bp spacing, 89 pairs; 345 of 1,334 instances; 214 pairs at T1 alone), while none
of the 1,693,612 same-strand A/C pairs at 3–4 bp spacing outside the motif carried both (permutation p < 10⁻⁴).
[...] (The contingency odds ratio for this comparison is saturated—no co-modified pair occurs outside the
motif—so the Haldane-corrected value, OR ≈ 6.1 × 10⁵ (95% CI 3.8 × 10⁴–9.7 × 10⁶), is uninformative as an effect
size; it is reported only in Supplementary Table 2.)
Per-read sentence (replaces the 'Single-molecule per-read analysis ... 244 read–position pairs' framing, which
describes a site-level census): Per-read analysis with the full read denominator (200,312 instance × read pairs
covering A₀/A₁ and C₄ at T1) showed 6mA on 47% of reads carrying 4mC at C₄ versus 24% of all reads (OR 4.84, 95% CI
4.73–4.95, ≥ 0.5 probability cutoff; OR 2.40 at ≥ 0.9; Supplementary Table 9).
