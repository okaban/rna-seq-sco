#!/usr/bin/env python3
"""
P1-3  Spatial-autocorrelation-aware re-test of the genome-wide GCCGGC
      dose-response (methylation density vs expression change).

Manuscript reports a Kruskal-Wallis association between GCCGGC site count
(dose group) and T2-vs-T1 log2FC with an extremely small naive p
(~4.1e-10; current pipeline 1.76e-12). Because both methylation density and
expression change are spatially autocorrelated along the chromosome, the naive
p is anti-conservative. We recompute significance with a circular-shift
(block) permutation that PRESERVES the spatial autocorrelation of LFC while
breaking its association with dose group.

Also reports the Simpson's-paradox structure: pooled vs core-only vs arm-only.

Inputs (existing canonical tables):
  41_GCCGGC_dose_response/tables/GCCGGC_gene_site_counts.tsv  (gene_start, region, dose_group)
  42_GCCGGC_temporal_derepression/tables/gene_methylation_transitions.tsv (LFC_T2vsT1)
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

np.random.seed(42)
BASE = os.path.expanduser("~/bioinfo/rna-seq/11_epigenome_integration/analysis")
DOSE = os.path.join(BASE, "41_GCCGGC_dose_response/tables/GCCGGC_gene_site_counts.tsv")
TRANS = os.path.join(BASE, "42_GCCGGC_temporal_derepression/tables/gene_methylation_transitions.tsv")
OUT = os.path.join(BASE, "78_reviewer_robustness/tables")

dose = pd.read_csv(DOSE, sep="\t")
trans = pd.read_csv(TRANS, sep="\t")[["gene_id", "LFC_T2vsT1", "LFC_T3vsT1"]]
df = dose.merge(trans, left_on="locus_tag", right_on="gene_id", how="inner")
df = df.dropna(subset=["LFC_T2vsT1", "gene_start", "dose_group", "region"]).copy()
df = df.sort_values("gene_start").reset_index(drop=True)
print(f"Merged genes with dose + LFC + position: {len(df)}")
print("dose_group counts:\n", df["dose_group"].value_counts().sort_index())
print("region counts:\n", df["region"].value_counts())


def kw_by_dose(sub):
    groups = [g["LFC_T2vsT1"].values for _, g in sub.groupby("dose_group")]
    groups = [g for g in groups if len(g) > 0]
    H, p = stats.kruskal(*groups)
    return H, p


def spearman_dose(sub):
    # dose_group as ordinal (0,1,2,3,4+ -> 0..4)
    order = {"0": 0, "1": 1, "2": 2, "3": 3, "4+": 4}
    x = sub["dose_group"].astype(str).map(order)
    return stats.spearmanr(x, sub["LFC_T2vsT1"])


# ---- naive KW (pooled / core / arm) ---------------------------------------
print("\n" + "=" * 66)
print("Simpson's-paradox structure (naive KW across dose groups)")
print("=" * 66)
H_all, p_all = kw_by_dose(df)
rho_all, prho_all = spearman_dose(df)
print(f"  POOLED : KW H={H_all:.2f}  p={p_all:.2e}   Spearman rho={rho_all:+.3f} p={prho_all:.2e}  (n={len(df)})")
for reg in ["core", "arm"]:
    sub = df[df["region"] == reg]
    if sub["dose_group"].nunique() > 1:
        H, p = kw_by_dose(sub)
        rho, prho = spearman_dose(sub)
        print(f"  {reg:6}: KW H={H:.2f}  p={p:.2e}   Spearman rho={rho:+.3f} p={prho:.2e}  (n={len(sub)})")


# ---- circular-shift permutation (preserves LFC autocorrelation) -----------
print("\n" + "=" * 66)
print("Circular-shift (block) permutation null  [preserves spatial autocorr]")
print("=" * 66)
N = 10000
lfc = df["LFC_T2vsT1"].values            # ordered by gene_start
dose_codes = df["dose_group"].astype(str).map(
    {"0": 0, "1": 1, "2": 2, "3": 3, "4+": 4}).values
obs_rho, _ = stats.spearmanr(dose_codes, lfc)
# observed KW
H_obs, _ = kw_by_dose(df)

rng = np.default_rng if hasattr(np, "default_rng") else None
rng = np.random.default_rng(42)
null_rho = np.empty(N)
null_H = np.empty(N)
n = len(lfc)
for i in range(N):
    shift = rng.integers(1, n)
    lfc_shift = np.roll(lfc, shift)       # circular shift keeps autocorrelation
    null_rho[i], _ = stats.spearmanr(dose_codes, lfc_shift)
    # KW with shifted LFC grouped by fixed dose
    groups = [lfc_shift[dose_codes == k] for k in np.unique(dose_codes)]
    null_H[i], _ = stats.kruskal(*[g for g in groups if len(g) > 0])

p_perm_rho = (np.sum(np.abs(null_rho) >= abs(obs_rho)) + 1) / (N + 1)
p_perm_H = (np.sum(null_H >= H_obs) + 1) / (N + 1)
print(f"  Observed Spearman rho (dose vs LFC) = {obs_rho:+.3f}")
print(f"  Block-permutation p (rho)           = {p_perm_rho:.4f}   "
      f"[null |rho| mean={np.mean(np.abs(null_rho)):.3f}, 95th={np.percentile(np.abs(null_rho),95):.3f}]")
print(f"  Observed KW H = {H_obs:.2f}")
print(f"  Block-permutation p (KW H)          = {p_perm_H:.4f}   "
      f"[null H 95th={np.percentile(null_H,95):.2f}]")

pd.DataFrame([dict(
    n=len(df), kw_H_pooled=H_all, kw_p_pooled_naive=p_all,
    spearman_rho_pooled=rho_all, spearman_p_pooled_naive=prho_all,
    block_perm_p_rho=p_perm_rho, block_perm_p_KW=p_perm_H,
    null_rho_mean_abs=float(np.mean(np.abs(null_rho))),
    null_rho_95=float(np.percentile(np.abs(null_rho), 95)),
)]).to_csv(os.path.join(OUT, "P1_3_spatial_autocorr.tsv"), sep="\t", index=False)

# region-stratified table
rows = []
for reg in ["pooled", "core", "arm"]:
    sub = df if reg == "pooled" else df[df["region"] == reg]
    if sub["dose_group"].nunique() > 1:
        H, p = kw_by_dose(sub); rho, prho = spearman_dose(sub)
        rows.append(dict(stratum=reg, n=len(sub), kw_H=H, kw_p_naive=p,
                         spearman_rho=rho, spearman_p_naive=prho))
pd.DataFrame(rows).to_csv(os.path.join(OUT, "P1_3_simpson_strata.tsv"),
                          sep="\t", index=False)
print("\nDONE -> 78_reviewer_robustness/tables/P1_3_*.tsv")
