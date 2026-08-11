#!/usr/bin/env python3
"""
P1-4 (TOST equivalence: permissive null -> positive equivalence claim)
P1-5 (continuous treatment + TSS-threshold sensitivity of Exposed/Shielded)

Input: non-circular regulatory-gene classification table (62 Exposed / 989 Shielded).
All from the methylation map alone (nearest GCCGGC distance at T1), with DESeq2 LFCs.

No new data; re-analysis of existing canonical table.
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

df = pd.read_csv(TBL, sep="\t")
df = df.dropna(subset=["LFC_T2vsT1"])
exp = df[df["class"] == "Exposed"]
shi = df[df["class"] == "Shielded"]
print(f"Loaded {len(df)} regulators with LFC_T2vsT1: "
      f"{len(exp)} Exposed / {len(shi)} Shielded")


# ---------------------------------------------------------------------------
def cliffs_delta(a, b):
    """Cliff's delta = P(a>b) - P(a<b), via rank method (exact, vectorised)."""
    a = np.asarray(a); b = np.asarray(b)
    n_gt = sum((a[:, None] > b[None, :]).sum(axis=1))
    n_lt = sum((a[:, None] < b[None, :]).sum(axis=1))
    return (n_gt - n_lt) / (len(a) * len(b))


def boot_cliffs_ci(a, b, n=10000, seed=42):
    rng = np.random.default_rng(seed)
    a = np.asarray(a); b = np.asarray(b)
    deltas = np.empty(n)
    for i in range(n):
        aa = rng.choice(a, size=len(a), replace=True)
        bb = rng.choice(b, size=len(b), replace=True)
        deltas[i] = cliffs_delta(aa, bb)
    return np.percentile(deltas, [2.5, 97.5])


def tost_welch(a, b, low, high):
    """Two one-sided Welch t-tests for equivalence of means within (low, high)
    of the (a-b) difference. Returns (p_tost, mean_diff, ci90)."""
    a = np.asarray(a); b = np.asarray(b)
    ma, mb = a.mean(), b.mean()
    diff = ma - mb
    se = np.sqrt(a.var(ddof=1)/len(a) + b.var(ddof=1)/len(b))
    # Welch df
    df_w = (se**2)**2 / ((a.var(ddof=1)/len(a))**2/(len(a)-1)
                         + (b.var(ddof=1)/len(b))**2/(len(b)-1))
    t_low = (diff - low) / se        # H0: diff <= low
    t_high = (diff - high) / se      # H0: diff >= high
    p_low = stats.t.sf(t_low, df_w)          # one-sided upper
    p_high = stats.t.cdf(t_high, df_w)       # one-sided lower
    p_tost = max(p_low, p_high)
    tcrit = stats.t.ppf(0.95, df_w)          # 90% CI
    ci90 = (diff - tcrit*se, diff + tcrit*se)
    return p_tost, diff, ci90


# ===========================================================================
# P1-4  EQUIVALENCE (permissive)
# ===========================================================================
print("\n" + "="*70)
print("P1-4  EQUIVALENCE TESTING (permissive: Exposed ~ Shielded)")
print("="*70)

rows = []
for lab, col in [("|LFC| T2vsT1", "LFC_T2vsT1"), ("|LFC| T3vsT1", "LFC_T3vsT1")]:
    e = exp[col].dropna().abs().values
    s = shi[col].dropna().abs().values
    U, p_mwu = stats.mannwhitneyu(e, s, alternative="two-sided")
    d = cliffs_delta(e, s)
    ci = boot_cliffs_ci(e, s)
    # parametric TOST on |LFC| means, margin +/-0.5 log2 units
    margin = 0.5
    p_tost, diff, ci90 = tost_welch(e, s, -margin, margin)
    delta_equiv = (ci[0] > -0.147) and (ci[1] < 0.147)
    print(f"\n[{lab}]  Exposed n={len(e)}  Shielded n={len(s)}")
    print(f"  median Exposed={np.median(e):.3f}  Shielded={np.median(s):.3f}")
    print(f"  MWU p={p_mwu:.3f}")
    print(f"  Cliff delta={d:+.3f}  boot95%CI=[{ci[0]:+.3f},{ci[1]:+.3f}]  "
          f"delta-equivalent(|d|<0.147)? {delta_equiv}")
    print(f"  TOST mean-diff={diff:+.3f}  90%CI=[{ci90[0]:+.3f},{ci90[1]:+.3f}]  "
          f"margin=+/-{margin}  p_TOST={p_tost:.4f}  "
          f"equivalent? {p_tost < 0.05}")
    rows.append(dict(contrast=lab, n_exp=len(e), n_shi=len(s),
                     med_exp=np.median(e), med_shi=np.median(s), mwu_p=p_mwu,
                     cliff_delta=d, cliff_lo=ci[0], cliff_hi=ci[1],
                     delta_equiv=delta_equiv, tost_diff=diff,
                     tost_ci90_lo=ci90[0], tost_ci90_hi=ci90[1],
                     tost_margin=margin, tost_p=p_tost, tost_equiv=p_tost < 0.05))

# DEG proportion equivalence (|LFC|>=1)
print("\n[DEG proportion |LFC_T2vsT1|>=1]")
e_deg = (exp["LFC_T2vsT1"].abs() >= 1)
s_deg = (shi["LFC_T2vsT1"].abs() >= 1)
pe, ps = e_deg.mean(), s_deg.mean()
# Newcombe-style 2-prop test + TOST on proportion diff, margin +/-0.10
n1, n2 = len(e_deg), len(s_deg)
x1, x2 = e_deg.sum(), s_deg.sum()
pdiff = pe - ps
se_p = np.sqrt(pe*(1-pe)/n1 + ps*(1-ps)/n2)
z_low = (pdiff - (-0.10))/se_p
z_high = (pdiff - 0.10)/se_p
p_tost_prop = max(stats.norm.sf(z_low), stats.norm.cdf(z_high))
ci90_p = (pdiff - 1.645*se_p, pdiff + 1.645*se_p)
print(f"  Exposed {x1}/{n1}={pe:.3f}  Shielded {x2}/{n2}={ps:.3f}  diff={pdiff:+.3f}")
print(f"  TOST 90%CI=[{ci90_p[0]:+.3f},{ci90_p[1]:+.3f}] margin=+/-0.10 "
      f"p_TOST={p_tost_prop:.4f} equivalent? {p_tost_prop < 0.05}")
rows.append(dict(contrast="DEG-prop T2vsT1", n_exp=n1, n_shi=n2,
                 med_exp=pe, med_shi=ps, mwu_p=np.nan,
                 cliff_delta=pdiff, cliff_lo=ci90_p[0], cliff_hi=ci90_p[1],
                 delta_equiv=np.nan, tost_diff=pdiff, tost_ci90_lo=ci90_p[0],
                 tost_ci90_hi=ci90_p[1], tost_margin=0.10, tost_p=p_tost_prop,
                 tost_equiv=p_tost_prop < 0.05))

pd.DataFrame(rows).to_csv(os.path.join(OUT, "P1_4_equivalence.tsv"),
                          sep="\t", index=False)

# ===========================================================================
# P1-5  CONTINUOUS treatment + TSS-threshold sensitivity
# ===========================================================================
print("\n" + "="*70)
print("P1-5  CONTINUOUS distance vs expression  (does binary lose signal?)")
print("="*70)

d2 = df.dropna(subset=["nearest_GCCGGC_T1", "LFC_T2vsT1"]).copy()
d2["logdist"] = np.log10(d2["nearest_GCCGGC_T1"] + 1)
d2["absLFC"] = d2["LFC_T2vsT1"].abs()
d2["arm"] = (d2["region_label"] == "arm").astype(int)

# raw Spearman (continuous distance vs signed and |LFC|)
for yl, yc in [("signed LFC_T2vsT1", "LFC_T2vsT1"), ("|LFC_T2vsT1|", "absLFC")]:
    rho, p = stats.spearmanr(d2["logdist"], d2[yc])
    print(f"  Spearman(logdist, {yl}) = {rho:+.3f}  p={p:.3f}  (n={len(d2)})")

# partial Spearman controlling core/arm (residual-rank method)
def partial_spearman(x, y, z):
    rx = pd.Series(x).rank(); ry = pd.Series(y).rank(); rz = pd.Series(z).rank()
    ex = rx - np.polyval(np.polyfit(rz, rx, 1), rz)
    ey = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
    r, p = stats.pearsonr(ex, ey)
    return r, p

rp, pp = partial_spearman(d2["logdist"].values, d2["LFC_T2vsT1"].values, d2["arm"].values)
print(f"  partial Spearman(logdist, signed LFC | core/arm) = {rp:+.3f}  p={pp:.3f}")
rp2, pp2 = partial_spearman(d2["logdist"].values, d2["absLFC"].values, d2["arm"].values)
print(f"  partial Spearman(logdist, |LFC|     | core/arm) = {rp2:+.3f}  p={pp2:.3f}")

# modality of the distance distribution (Hartigan dip)
try:
    import diptest
    dipstat, dip_p = diptest.diptest(d2["logdist"].values)
    print(f"  Hartigan dip (log10 nearest-distance): D={dipstat:.4f}  p={dip_p:.3f}")
except Exception as ex:
    dipstat, dip_p = np.nan, np.nan
    print(f"  diptest unavailable: {ex}")

# threshold sensitivity: vary Exposed window, re-test permissive null
print("\n  Threshold sensitivity (Exposed = nearest_GCCGGC_T1 <= W):")
sens = []
for W in range(100, 501, 50):
    is_exp = d2["nearest_GCCGGC_T1"] <= W
    e = d2.loc[is_exp, "absLFC"].values
    s = d2.loc[~is_exp, "absLFC"].values
    if len(e) < 5:
        continue
    U, p = stats.mannwhitneyu(e, s, alternative="two-sided")
    dd = cliffs_delta(e, s)
    print(f"    W={W:>3}bp  n_Exposed={len(e):>3}  MWU p={p:.3f}  Cliff d={dd:+.3f}")
    sens.append(dict(window_bp=W, n_exposed=int(is_exp.sum()),
                     mwu_p=p, cliff_delta=dd))
pd.DataFrame(sens).to_csv(os.path.join(OUT, "P1_5_threshold_sensitivity.tsv"),
                          sep="\t", index=False)

# save continuous summary
pd.DataFrame([dict(
    n=len(d2),
    spearman_logdist_signedLFC=stats.spearmanr(d2["logdist"], d2["LFC_T2vsT1"])[0],
    spearman_logdist_signedLFC_p=stats.spearmanr(d2["logdist"], d2["LFC_T2vsT1"])[1],
    partial_logdist_signedLFC=rp, partial_logdist_signedLFC_p=pp,
    partial_logdist_absLFC=rp2, partial_logdist_absLFC_p=pp2,
    dip_D=dipstat, dip_p=dip_p)]).to_csv(
    os.path.join(OUT, "P1_5_continuous_summary.tsv"), sep="\t", index=False)

print("\nDONE. Tables -> 78_reviewer_robustness/tables/")
