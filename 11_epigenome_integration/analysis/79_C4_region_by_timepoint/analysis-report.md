# C4 — Genomic-region composition of methylation by system × timepoint

**Question.** Does the developmental redistribution of methylation (core→arm,
and the T1→T3 decline in site number) change *which genomic regions* the marks
occupy, or only *where on the chromosome* they sit?

**Data / unit.** Canonical census sites (depth ≥ 10, freq ≥ 50%), per timepoint,
for three systems: GCCGGC m4C (`final_motif=TGGCCGGC`), AAGCCCG m4C, AAGCCCG 6mA.
Counts match the canonical numbers used throughout the paper
(GCCGGC 1289/407/21; AAGCCCG-4mC 587/212/15; AAGCCCG-6mA 260/64/38).

**Method.** Each site classified into promoter (TSS −500..−1), 5′UTR (in-gene
0..+100), CDS-internal, or intergenic. Classification is **in-gene first** then
promoter — order-independent, and validated to reproduce the SuppFig 9 classifier
(`62_GO_KEGG_enrichment`) at **99.9% agreement** on its own positions. (A naive
order-dependent first-match rule over-assigns promoter by ~4× when a site lies
inside one gene and within −500 bp of a neighbour's TSS; that rule was rejected.)

**Result.** All three systems are overwhelmingly CDS-internal at every timepoint
(GCCGGC m4C 86/84/86%; AAGCCCG m4C 80/84/53%; AAGCCCG 6mA 86/86/68%), with a
small, broadly stable promoter fraction (~4–9%). The region composition is
essentially constant across the developmental time course even as the absolute
site number collapses (e.g. GCCGGC 1289→21) and the marks relocate from core to
arms. T3 fractions for the small-n strata (GCCGGC n=21, AAGCCCG-4mC n=15,
AAGCCCG-6mA n=38) are noisy and shown only for completeness.

**Interpretation / decision.** The redistribution is a **spatial** relocation,
not a redirection of methylation toward or away from particular functional
regions. This is consistent with the permissive / non-targeted model: the
systems mark active chromosomal *territory*, and their functional-region profile
does not change with development.

**Caveat.** T3 strata are small (n=15–38); per-region percentages there are
unstable and not interpreted beyond "still predominantly CDS-internal."

Figure: `figures/C4_region_by_timepoint.png` · Table:
`tables/C4_region_composition.tsv` · Script: `scripts/c4_region_by_timepoint.py`.
Manuscript: Supplementary Figure 14.
