#!/usr/bin/env python3
"""Recompute tables/expression_quintile.tsv with the CANONICAL Exposed set (n = 62).

2026-09-13 (review finding FIG-06 / S22 panel c). The stored table dated
2026-05-23 was produced by H29_shielded_exposed_boundary.py §10 while the
unified table's `is_exposed` still held the expression-selected n = 57 set
(now kept as `is_exposed_STALE_n57`); its n_exposed column summed to 57.

Method is H29 §10 verbatim: baseMean quintiles by pd.qcut(q=5) over the
regulatory genes with a baseMean, Exposed fraction per quintile, and H29's
Jonckheere–Terpstra Z approximation (no tie correction). Exposed membership is
taken from SuppTable1_Exposed62_identity.tsv and cross-checked against the
unified table's `is_exposed` column (must be identical).
"""
import numpy as np, pandas as pd
from itertools import combinations
from scipy import stats
from pathlib import Path

A = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/52_shielded_exposed_boundary")
TBL = A / "tables"

reg = pd.read_csv(TBL / "all_genes_features_unified_n57.tsv", sep="\t")
can = set(pd.read_csv(TBL / "SuppTable1_Exposed62_identity.tsv", sep="\t")["locus_tag"])
reg["is_exposed_canonical"] = reg["locus_tag"].isin(can).astype(int)
assert reg["is_exposed_canonical"].sum() == 62, reg["is_exposed_canonical"].sum()
assert (reg["is_exposed_canonical"] == reg["is_exposed"].astype(int)).all(), \
    "unified-table is_exposed differs from the canonical 62 — stop and report"

df = reg.dropna(subset=["baseMean"]).copy()
df["baseMean_quintile"] = pd.qcut(df["baseMean"], q=5, labels=False, duplicates="drop") + 1

rows = []
for q in sorted(df["baseMean_quintile"].unique()):
    sub = df[df["baseMean_quintile"] == q]
    rows.append({
        "quintile": int(q), "n_genes": len(sub),
        "n_exposed": int(sub["is_exposed_canonical"].sum()),
        "frac_exposed": round(sub["is_exposed_canonical"].mean(), 4),
        "mean_baseMean": round(sub["baseMean"].mean(), 1),
        "median_baseMean": round(sub["baseMean"].median(), 1),
        "mean_nearest_methyl_distance": round(sub["nearest_methyl_distance"].mean(), 1),
        "median_nearest_methyl_distance": round(sub["nearest_methyl_distance"].median(), 1),
        "mean_n_methyl_sites_2kb": round(sub["n_methyl_sites_2kb"].mean(), 2),
    })
quint = pd.DataFrame(rows)


def jonckheere_terpstra(groups, values):
    """H29 §10 implementation (Z approximation, no tie correction)."""
    ug = sorted(set(groups)); jt = 0.0
    for i, j in combinations(range(len(ug)), 2):
        vi = values[groups == ug[i]]; vj = values[groups == ug[j]]
        gt = (vj[None, :] > vi[:, None]).sum(); eq = (vj[None, :] == vi[:, None]).sum()
        jt += gt + 0.5 * eq
    n = len(values); sizes = [int(np.sum(groups == g)) for g in ug]
    E = (n ** 2 - sum(s ** 2 for s in sizes)) / 4
    var = (n ** 2 * (2 * n + 3) - sum(s ** 2 * (2 * s + 3) for s in sizes)) / 72
    z = (jt - E) / np.sqrt(var)
    return jt, z, 2 * (1 - stats.norm.cdf(abs(z)))


jt, z, p = jonckheere_terpstra(df["baseMean_quintile"].values, df["is_exposed_canonical"].values)

old = TBL / "expression_quintile.tsv"
if old.exists():
    prev = pd.read_csv(old, sep="\t")
    arch = TBL / "expression_quintile_ARCHIVED_260913_n57.tsv"
    if not arch.exists():
        prev.to_csv(arch, sep="\t", index=False)
    print("previous table n_exposed:", prev["n_exposed"].tolist(), "sum", int(prev["n_exposed"].sum()))
quint.to_csv(old, sep="\t", index=False)
pd.DataFrame([{"test": "Jonckheere-Terpstra (H29 Z approx)", "outcome": "is_exposed_canonical(62)",
               "groups": "baseMean quintile", "n": len(df), "JT": jt, "z": z, "p": p}]).to_csv(
    TBL / "expression_quintile_JT_n62.tsv", sep="\t", index=False)
print(quint[["quintile", "n_genes", "n_exposed", "frac_exposed"]].to_string(index=False))
print(f"n genes with baseMean = {len(df)}; Exposed with baseMean = {int(df['is_exposed_canonical'].sum())} of 62")
print(f"JT stat={jt:.1f} z={z:.3f} p={p:.3f}")
