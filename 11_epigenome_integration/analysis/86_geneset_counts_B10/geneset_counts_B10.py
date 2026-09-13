#!/usr/bin/env python3
"""86 — ledger B10: (a) Exposed/Shielded split among regulatory genes with a defined
LFC_T2vsT1; (b) nesting/overlap of the three motif-proximal gene sets used for the
GO/KEGG enrichment (EN L94: GCCGGC-proximal 3,197 / AAGCCCG-proximal 459 / dual-targeted 448).

Sources:
 (a) 52_shielded_exposed_boundary/tables/SuppTable_regulatory_gene_classification_n1051.tsv
     (canonical identity table; column `class`, `LFC_T2vsT1`).
 (b) 40_motif_division_of_labor/tables/gene_set_comparison.tsv (current, written 2026-04-19 21:59;
     column `exclusive_group`), and its superseded predecessor gene_set_comparison_s0_archive.tsv
     (21:57 the same day). Groups are assigned by H17_motif_division_of_labor.py::assign_exclusive_group
     -> the three sets are MUTUALLY EXCLUSIVE by construction ('GCCGGC-only', 'AAGCCCG-only',
     'Dual-targeted'). 62_GO_KEGG_enrichment/scripts/GO_KEGG_enrichment.py builds its query sets by
     substring match on `exclusive_group` ('GCCGGC' -> GCCGGC-only, 'AAGCCCG' -> AAGCCCG-only,
     == 'Dual-targeted'), i.e. the enrichment was run on the exclusive sets.
"""
import pandas as pd
from pathlib import Path
B = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis"); H = Path(__file__).resolve().parent
reg = pd.read_csv(B/"52_shielded_exposed_boundary/tables/SuppTable_regulatory_gene_classification_n1051.tsv", sep="\t")
rows = [dict(universe="all TSS-assignable regulatory genes", n=len(reg), **reg["class"].value_counts().to_dict())]
d = reg.dropna(subset=["LFC_T2vsT1"]); rows.append(dict(universe="with defined LFC_T2vsT1", n=len(d), **d["class"].value_counts().to_dict()))
d3 = reg.dropna(subset=["LFC_T3vsT1"]); rows.append(dict(universe="with defined LFC_T3vsT1", n=len(d3), **d3["class"].value_counts().to_dict()))
split = pd.DataFrame(rows); split.to_csv(H/"tables/B10_exposed_shielded_split.tsv", sep="\t", index=False); print(split.to_string(index=False))
out = []
for tag, f in (("current", "gene_set_comparison.tsv"), ("s0_archive (superseded)", "gene_set_comparison_s0_archive.tsv")):
    g = pd.read_csv(B/"40_motif_division_of_labor/tables"/f, sep="\t")
    G = set(g.gene_id[g.exclusive_group == "GCCGGC-only"]); A = set(g.gene_id[g.exclusive_group == "AAGCCCG-only"]); D = set(g.gene_id[g.exclusive_group == "Dual-targeted"])
    out.append(dict(table=tag, n_genes_in_table=len(g), GCCGGC_only=len(G), AAGCCCG_only=len(A), Dual_targeted=len(D), Background=int((g.exclusive_group == "Background").sum()),
                    GCCGGC_only_and_Dual=len(G & D), AAGCCCG_only_and_Dual=len(A & D), GCCGGC_only_and_AAGCCCG_only=len(G & A),
                    GCCGGC_proximal_inclusive=len(G | D), AAGCCCG_proximal_inclusive=len(A | D),
                    Dual_subset_of_GCCGGC_inclusive=D <= (G | D), Dual_subset_of_AAGCCCG_inclusive=D <= (A | D)))
nest = pd.DataFrame(out); nest.to_csv(H/"tables/B10_geneset_nesting.tsv", sep="\t", index=False); print(nest.T.to_string())
