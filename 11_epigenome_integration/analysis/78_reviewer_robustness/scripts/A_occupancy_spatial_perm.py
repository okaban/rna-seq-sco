"""A: spatial-autocorrelation permutation on the CANONICAL occupancy-based
permissive bias (the central r = -0.09, p = 0.004).

P1-3b tested the weak bias with a distance proxy. Here we reproduce the exact
canonical statistic from Exposed_dynamicA_define_and_verify.py — region-
controlled partial correlation of T1 promoter GCCGGC m4C occupancy (summed
modification frequency within +-2 kb of TSS) versus T2-vs-T1 log2FC — and
attach a circular-shift (spatial-autocorrelation-preserving) permutation p
directly to it.
"""
import numpy as np, pandas as pd
from scipy.stats import pearsonr
from pathlib import Path

np.random.seed(42)
B = Path('/Users/okaban/bioinfo/rna-seq')
reg = pd.read_csv(B/'11_epigenome_integration/analysis/52_shielded_exposed_boundary/tables/all_genes_features_unified_n57.tsv', sep='\t')
g = pd.read_csv(B/'11_epigenome_integration/analysis/37_defense_island_GCCGGC/tables/GCCGGC_sites_by_timepoint.tsv', sep='\t')
OUT = B/'11_epigenome_integration/analysis/78_reviewer_robustness/tables'

r = reg.dropna(subset=['tss']).copy(); r['tss'] = r['tss'].astype(int)
gT1 = g[g.timepoint == 'T1']
P = np.sort(gT1.position.values)
Fq = gT1.sort_values('position').frequency.values

def occ(tss, w=2000):
    m = (P >= tss - w) & (P <= tss + w)
    return Fq[m].sum() if m.any() else 0.0

r['occ2k'] = r['tss'].apply(occ)
d = r[['occ2k', 'LFC_T2vsT1', 'region', 'tss']].dropna().sort_values('tss').reset_index(drop=True)
n = len(d)

def partial_r(occ, lfc, region):
    rocc = pd.Series(occ).rank().values
    rlfc = pd.Series(lfc).rank().values
    rreg = pd.Series(region).rank().values
    rx = rocc - np.polyval(np.polyfit(rreg, rocc, 1), rreg)
    ry = rlfc - np.polyval(np.polyfit(rreg, rlfc, 1), rreg)
    return pearsonr(rx, ry)

occ_v = d.occ2k.values; lfc_v = d.LFC_T2vsT1.values; reg_v = d.region.values
r_obs, p_obs = partial_r(occ_v, lfc_v, reg_v)
print(f"n={n}")
print(f"CANONICAL occupancy x LFC_T2, region-controlled partial r = {r_obs:+.3f}  naive p = {p_obs:.4f}")

# circular-shift permutation: shift LFC along TSS order, keep (occ, region) fixed
N = 10000
rng = np.random.default_rng(42)
null = np.empty(N)
for i in range(N):
    s = rng.integers(1, n)
    null[i], _ = partial_r(occ_v, np.roll(lfc_v, s), reg_v)
p_perm = (np.sum(np.abs(null) >= abs(r_obs)) + 1) / (N + 1)
print(f"circular-shift block-permutation p = {p_perm:.4f}  "
      f"[null |r| mean={np.mean(np.abs(null)):.3f}, 95th={np.percentile(np.abs(null),95):.3f}]")
verdict = "SURVIVES" if p_perm < 0.05 else "does NOT survive"
print(f">>> canonical occupancy-based weak bias {verdict} spatial-autocorrelation correction")

pd.DataFrame([dict(n=n, partial_r=r_obs, naive_p=p_obs, block_perm_p=p_perm,
                   null_abs_mean=float(np.mean(np.abs(null))),
                   null_abs_95=float(np.percentile(np.abs(null), 95)),
                   verdict=verdict)]).to_csv(OUT/'A_occupancy_spatial_perm.tsv', sep='\t', index=False)
print("saved -> A_occupancy_spatial_perm.tsv")
