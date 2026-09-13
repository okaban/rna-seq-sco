#!/usr/bin/env python3
"""80 — Partial Spearman correlation of promoter GCCGGC 4mC occupancy (T1) vs
LFC_T2vsT1 in regulatory genes, under three position covariates (ledger B04).

Canon (EN L85): partial r = -0.090, 95% CI [-0.15,-0.03], p = 0.0041, n = 1,019.
Manuscript also quotes (no producing script existed): continuous oriC-distance
covariate r = -0.09, p = 0.0028; both jointly r = -0.10, p = 0.0018.

Definitions (identical to 52_shielded_exposed_boundary/scripts/
Exposed_dynamicA_define_and_verify.py and 78_reviewer_robustness/scripts/
A_occupancy_spatial_perm.py, which reproduce the canon):
  x  = occ2k  : sum of T1 GCCGGC 4mC site frequencies (%) within TSS +/- 2,000 bp
                (sites: 37_defense_island_GCCGGC/tables/GCCGGC_sites_by_timepoint.tsv, T1 rows)
  y  = LFC_T2vsT1 (DESeq2)
  z1 = core/arm binary (region_label; core = 1.5-7.17 Mb)
  z2 = |TSS - oriC|, oriC = 4,271,763 bp (dnaA SCO3879 midpoint; same constant as
       83_oriC_confound_FIRE/scripts/oric_confound.py)
  z3 = (z1, z2) jointly
Partial Spearman = Pearson r between the residuals of rank(x) and rank(y) after
OLS on rank(z) (ranks of each covariate; multi-covariate OLS for z3). p = naive
Pearson p on the residuals (as in the canon script). 95% CI = percentile
bootstrap over genes, 1,000 resamples, numpy default_rng(seed=1).

Gene universe: 52_shielded_exposed_boundary/tables/
SuppTable_regulatory_gene_classification_n1051.tsv (canonical identity table,
1,051 genes; 1,019 with LFC_T2vsT1). all_genes_features_unified_n57.tsv is NOT
used (it carries is_exposed_STALE_n57).
"""
import numpy as np, pandas as pd
from pathlib import Path
from scipy.stats import pearsonr, spearmanr
B = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis")
H = Path(__file__).resolve().parent
ORIC = 4_271_763.0
reg = pd.read_csv(B/"52_shielded_exposed_boundary/tables/SuppTable_regulatory_gene_classification_n1051.tsv", sep="\t")
g = pd.read_csv(B/"37_defense_island_GCCGGC/tables/GCCGGC_sites_by_timepoint.tsv", sep="\t")
gT1 = g[g.timepoint == "T1"].sort_values("position")
P, Fq = gT1.position.values, gT1.frequency.values
def occ(tss, w=2000):
    m = (P >= tss - w) & (P <= tss + w); return Fq[m].sum() if m.any() else 0.0
d = reg.dropna(subset=["tss", "LFC_T2vsT1"]).copy()
d["tss"] = d.tss.astype(int); d["occ2k"] = d.tss.apply(occ)
d["arm"] = (d.region_label == "arm").astype(int); d["dist_oriC"] = (d.tss - ORIC).abs()
n = len(d); assert n == 1019, n
x, y = d.occ2k.values, d.LFC_T2vsT1.values
Z = {"core_arm_binary": d[["arm"]].values.astype(float),
     "dist_oriC_continuous": d[["dist_oriC"]].values,
     "both_jointly": d[["arm", "dist_oriC"]].values.astype(float)}
def rank(v): return pd.Series(v).rank().values
def partial_spearman(x, y, Zm):
    rx, ry = rank(x), rank(y)
    A = np.column_stack([np.ones(len(x))] + [rank(Zm[:, j]) for j in range(Zm.shape[1])])
    ex = rx - A @ np.linalg.lstsq(A, rx, rcond=None)[0]
    ey = ry - A @ np.linalg.lstsq(A, ry, rcond=None)[0]
    return pearsonr(ex, ey)
rows = []
rho0, p0 = spearmanr(x, y)
rows.append(dict(covariate="none (raw Spearman)", n=n, r=rho0, p=p0, ci_lower=np.nan, ci_upper=np.nan))
rng = np.random.default_rng(1); NB = 1000
for name, Zm in Z.items():
    r, p = partial_spearman(x, y, Zm)
    boots = np.empty(NB)
    for b in range(NB):
        idx = rng.integers(0, n, n); boots[b] = partial_spearman(x[idx], y[idx], Zm[idx])[0]
    lo, hi = np.percentile(boots, [2.5, 97.5])
    rows.append(dict(covariate=name, n=n, r=r, p=p, ci_lower=lo, ci_upper=hi))
out = pd.DataFrame(rows); out.to_csv(H/"tables/partial_corr_covariates.tsv", sep="\t", index=False, float_format="%.6g")
print(out.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
# Exposed/Shielded split among the n=1,019 (B10)
print("class split among n with LFC:", d["class"].value_counts().to_dict())
print("core/arm among n:", d.region_label.value_counts().to_dict(), "| Spearman(dist_oriC, arm) =", round(spearmanr(d.dist_oriC, d.arm)[0], 3))

# ---- sensitivity: alternative oriC-distance transforms (which definition did the manuscript use?) ----
# The quoted oriC values (r = -0.09, p = 0.0028; joint r = -0.10, p = 0.0018) are not reproduced by the
# rank-based |TSS - oriC| covariate above. Alternatives are tabulated so the editor can see what comes closest.
def ps_lin(x, y, Zm):   # covariates entered untransformed (not ranked)
    rx, ry = rank(x), rank(y); A = np.column_stack([np.ones(len(x)), Zm])
    ex = rx - A @ np.linalg.lstsq(A, rx, rcond=None)[0]; ey = ry - A @ np.linalg.lstsq(A, ry, rcond=None)[0]
    return pearsonr(ex, ey)
tss = d.tss.values.astype(float); arm = d.arm.values.astype(float); dist = d.dist_oriC.values
sens = {"|d| ranked": (partial_spearman, dist[:, None]), "|d| linear": (ps_lin, dist[:, None]),
        "log10(|d|+1) linear": (ps_lin, np.log10(dist + 1)[:, None]), "log10(|d|+1) ranked (= |d| ranked)": (partial_spearman, np.log10(dist + 1)[:, None]),
        "signed (TSS-oriC) linear": (ps_lin, (tss - ORIC)[:, None]),
        "|d| linear + arm": (ps_lin, np.column_stack([dist, arm])), "log10(|d|+1) linear + arm": (ps_lin, np.column_stack([np.log10(dist + 1), arm])),
        "signed linear + arm": (ps_lin, np.column_stack([tss - ORIC, arm]))}
srows = [dict(definition=k, r=f(x, y, Zm)[0], p=f(x, y, Zm)[1]) for k, (f, Zm) in sens.items()]
sdf = pd.DataFrame(srows); sdf.to_csv(H/"tables/partial_corr_oriC_definition_sensitivity.tsv", sep="\t", index=False, float_format="%.6g")
print(sdf.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
