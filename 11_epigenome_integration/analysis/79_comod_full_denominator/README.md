# 79 — AAGCCCG per-read co-modification, full read denominator (2026-09-13)

Why: `73_comod_threshold_ROC/threshold_sensitivity.tsv` (source of Supp Table 9)
built its 2x2 only over reads with >=1 modification-call row, so unmodified
reads were structurally excluded ("neither" = 222 at t=0.5). That Berkson-type
selection drives the odds ratio below 1 at every threshold (0.0067 ... 0.35),
which an external reviewer read as negative association.

Here every T1 read covering both an A (A0/A1) and a C (C3/C5) position of the
same motif instance enters the table, from MM/ML tags via pysam
(464,450 reads scanned; 200,830 instance x read pairs; 1,334 instances,
649 + / 685 -).

| t | both | A only | C only | neither | OR | fold vs independence | p |
|---|---|---|---|---|---|---|---|
| 0.50 | 1,408 | 46,866 | 1,401 | 151,155 | 3.24 | 2.09 | 5.5e-199 |
| 0.75 | 298 | 25,977 | 787 | 173,768 | 2.53 | 2.10 | 2.6e-36 |
| 0.90 | 40 | 11,849 | 384 | 188,557 | 1.66 | 1.59 | 3.8e-03 |

Conclusion: per-read co-modification is a POSITIVE association at every
threshold, ~2-fold over independence; P(6mA | 4mC on the same read) ~ 0.50 vs
P(6mA) ~ 0.24 at t=0.5. Supp Table 9 must be rebuilt from this table and the
body text must not quote the 73_ ORs. Effect size is modest (OR 2-3) and
distinct from the site-level OR (138,440) — say so.

Run: `python comod_full_denominator.py /Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa`

## 2026-09-13 additions (ledger B01)

- `build_suppTable9.py` -> `SuppTable9_comod_full_denominator_T1.{tsv,md}`: Supp Table 9 replacement from
  `comod_full_denominator_T1.tsv` with log-OR normal 95% CI, P(6mA|4mC), P(6mA).
- **Offset caveat (critical).** The modkit T1 pileup (`methyl/260102_M145/analysis/pileup/1-1_pileup.bed`, cov ≥ 10)
  places the AAGCCCG 4mC call at motif offset 4 (AAGCC*CG: mean 69%, ≥50% at 88% (+) / 81% (−) of sites) and ~0%
  at offsets 3 and 5 (mean 0.03–0.22%). `comod_full_denominator.py` scans C_OFF = (3,5), so its "C-positive" reads
  (P(4mC) = 1.4% at t = 0.5) are not the AAGCCCG 4mC mark. `comod_full_denominator_C4.py` repeats the scan at offset 4
  -> `comod_full_denominator_T1_C4.tsv`.
- `site_census_AAGCCCG_pairs.py` -> `site_census_AAGCCCG_pairs{,_summary}.tsv`, `site_census_AAGCCCG_instances_summary.tsv`:
  site-level co-localisation census with correct 0-based motif-offset mapping. The `position` column of
  `07_motif_analysis/methylation_site_sequences.csv` is 0-based (ref[position] is the modified base for all four
  mod/strand classes) while its 31-bp `sequence` window is centred one base upstream; the earlier censuses (244/155/85 in
  23_expanded_motif_search/02 & 68_or_permutation; 254 in 65_per_read_comod) classified sites through that window and
  therefore dropped sites and mislabelled offsets. Canonical AAGCCCG sites: 4mC only at C4 (482 +, 494 −); 6mA at A0
  (208 +, 212 −) and A1 (73 +, 74 −). Same-strand, same-instance 4mC–6mA pairs: T1 214 (A0–C4 181, A1–C4 33; 194 of
  1,334 instances), T2 26 (19/7), T3 0; pooled over timepoints 406 (A0–C4 317, A1–C4 89; 345 instances). With the
  65_ rule (same strand, ≤ 10 bp) pooled = 415 (9 cross-instance). The pairs therefore have 4 bp (A0–C4) and 3 bp
  (A1–C4) spacing — same spacings as the manuscript's "A₁↔C₅ / A₀↔C₃", but the offset labels differ.
