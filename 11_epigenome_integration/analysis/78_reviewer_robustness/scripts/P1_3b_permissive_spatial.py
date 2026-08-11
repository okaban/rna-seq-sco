#!/usr/bin/env python3
"""
P1-3b  Does the central 'weak modulatory bias' permissive correlation survive
       spatial-autocorrelation correction?

The manuscript's load-bearing permissive number is a geography-controlled
partial correlation between promoter GCCGGC methylation and T1->T2 expression
change in regulatory genes (partial Spearman r = -0.09, p = 0.004). P1-3 showed
the genome-wide dose-response association vanishes under a circular-shift
(spatial-autocorrelation-preserving) permutation. Here we apply the same test
to the regulatory-gene permissive correlation, using nearest-GCCGGC distance at
T1 as the (continuous) methylation proxy available in the canonical table.

If it does NOT survive -> the honest statement strengthens from 'weak bias
r=-0.09' to 'no bias detectable after autocorrelation correction' (= fully
permissive). This touches a locked number, so it is reported for decision only.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

np.random.seed(42)
BASE = os.path.expanduser("~/bioinfo/rna-seq/11_epigenome_integration/analysis")
TBL = os.path.join(BASE, "52_shielded_exposed_boundary/tables/"
                   "SuppTable_regulatory_gene_classification_n1051.tsv")
OUT = os.path.join(BASE, "78_reviewer_robustness/tables")

df = pd.read_csv(TBL, sep="\t").dropna(subset=["LFC_T2vsT1", "nearest_GCCGGC_T1", "tss"]).copy()
df = df.sort_values("tss").reset_index(drop=True)
df["logdist"] = np.log10(df["nearest_GCCGGC_T1"] + 1)
df["arm"] = (df["region_label"] == "arm").astype(int)
n = len(df)
print(f"Regulatory genes ordered by TSS: {n}")


def partial_spearman_resid(x, y, z):
    rx = pd.Series(x).rank().values
    ry = pd.Series(y).rank().values
    rz = pd.Series(z).rank().values
    ex = rx - np.polyval(np.polyfit(rz, rx, 1), rz)
    ey = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
    r, p = stats.pearsonr(ex, ey)
    return r, p


x = df["logdist"].values
y = df["LFC_T2vsT1"].values
z = df["arm"].values

# observed
rho_raw, p_raw = stats.spearmanr(x, y)
rho_par, p_par = partial_spearman_resid(x, y, z)
print(f"\nObserved (continuous distance proxy):")
print(f"  raw     Spearman(logdist, signed LFC)        = {rho_raw:+.3f}  naive p={p_raw:.2e}")
print(f"  partial Spearman(logdist, signed LFC|region) = {rho_par:+.3f}  naive p={p_par:.2e}")

# circular-shift permutation: shift LFC along TSS order, keep (dist, region) fixed
N = 10000
rng = np.random.default_rng(42)
null_raw = np.empty(N)
null_par = np.empty(N)
for i in range(N):
    s = rng.integers(1, n)
    ys = np.roll(y, s)
    null_raw[i], _ = stats.spearmanr(x, ys)
    null_par[i], _ = partial_spearman_resid(x, ys, z)

p_perm_raw = (np.sum(np.abs(null_raw) >= abs(rho_raw)) + 1) / (N + 1)
p_perm_par = (np.sum(np.abs(null_par) >= abs(rho_par)) + 1) / (N + 1)
print(f"\nCircular-shift (spatial-autocorrelation-preserving) permutation:")
print(f"  raw     block-perm p = {p_perm_raw:.4f}  "
      f"[null |rho| 95th = {np.percentile(np.abs(null_raw),95):.3f}]")
print(f"  partial block-perm p = {p_perm_par:.4f}  "
      f"[null |rho| 95th = {np.percentile(np.abs(null_par),95):.3f}]")
verdict = ("SURVIVES (weak bias is real)" if p_perm_par < 0.05
           else "DOES NOT SURVIVE (bias is an autocorrelation artefact -> fully permissive)")
print(f"\n  >>> Permissive weak-bias correlation: {verdict}")

pd.DataFrame([dict(
    n=n, rho_raw=rho_raw, p_raw_naive=p_raw, rho_partial=rho_par, p_partial_naive=p_par,
    block_perm_p_raw=p_perm_raw, block_perm_p_partial=p_perm_par,
    null_partial_95=float(np.percentile(np.abs(null_par), 95)),
    verdict=verdict)]).to_csv(
    os.path.join(OUT, "P1_3b_permissive_spatial.tsv"), sep="\t", index=False)
print("\nDONE -> P1_3b_permissive_spatial.tsv")
