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
