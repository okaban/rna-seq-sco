# 86 — gene-set counts for ledger B10 (2026-09-13)

`geneset_counts_B10.py` -> `tables/B10_exposed_shielded_split.tsv`, `tables/B10_geneset_nesting.tsv`

- Regulatory genes: 1,051 TSS-assignable = 62 Exposed / 989 Shielded; with defined LFC_T2vsT1: 1,019 = 59 Exposed / 960 Shielded
  (identical for LFC_T3vsT1).
- Motif-proximal gene sets (`40_motif_division_of_labor/tables/gene_set_comparison.tsv`, column `exclusive_group`) are
  MUTUALLY EXCLUSIVE by construction: GCCGGC-only / AAGCCCG-only / Dual-targeted (pairwise overlaps = 0). The GO/KEGG
  enrichment (62_) was run on these exclusive sets. "Dual ⊂ GCCGGC-proximal" holds only for the inclusive definition
  (GCCGGC-only ∪ Dual).
- Current table (2026-04-19 21:59; n = 7,374 genes with LFC in both contrasts): 3,106 / 438 / 443 (inclusive 3,549 / 881).
  Superseded `_s0_archive` (21:57; n = 7,646): 3,197 / 459 / 448 (inclusive 3,645 / 907) — these are the numbers in EN L94.
  Group labels are identical for all shared genes; the difference is the expression universe (278 genes dropped, 6 added).
  The KEGG statistics quoted at EN L96 (sco02020: 64/115, enrichment 1.32, OR 1.74, FDR 0.129) match the CURRENT
  `E1_KEGG_enrichment_GCCGGC-proximal_4mC.tsv`, not the archive (64/116, 1.32, 1.73, 0.139).
