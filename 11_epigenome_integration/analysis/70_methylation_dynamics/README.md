
## 2026-09-21 — coordinate-offset fix (Supp Fig 5 / Supp Table 4 regenerated)

`high_confidence_sites_weighted.csv` `position` is 0-based (modkit bedMethyl start). `assign_motifs_fast` subtracted 1
(treating it as 1-based) and matched motif offsets without regard to strand, so it recovered only 147 of the AAGCCCG 6mA
sites present at all three timepoints (309 with the corrected, strand-aware matcher identical to 90_'s canonical census).
GCCGGC is unaffected (855 both ways: the palindromic site is symmetric under the shift). Result at the ≥ 10 pp criterion:
AAGCCCG 97 of 309 dynamic (31.4%; was 47 of 147, 32.0%) — same fraction, doubled denominator; GCCGGC 243 of 855 unchanged.
Old outputs: `tables/superseded_260921_1based/`, `figures/superseded_260921_1based/`.
