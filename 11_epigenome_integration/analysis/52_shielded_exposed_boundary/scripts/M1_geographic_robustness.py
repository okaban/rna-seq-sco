#!/usr/bin/env python3
"""
M1 (peer-review pre-emption): Is the Shielded/Exposed partition (nearest methylation
distance -> is_exposed) explained by core/arm geographic confounding?

- Recompute region cleanly for ALL genes from gene midpoint (fixes mixed 0/1 + arm/core coding).
- Overall AUC, per-stratum AUC, Mann-Whitney within strata.
- Logistic regression: is_exposed ~ log10(distance) + region  -> distance effect adjusted for geography.
Outputs: tables/M1_geographic_robustness.tsv, figures/M1_stratified_ROC.png/.pdf
"""
import csv, numpy as np
from pathlib import Path
from scipy.stats import rankdata, mannwhitneyu, fisher_exact
import statsmodels.api as sm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve

ARM_LEFT_MAX = 1_500_000
ARM_RIGHT_MIN = 7_167_508
HERE = Path(__file__).resolve().parent.parent
TBL = HERE / "tables"; FIG = HERE / "figures"
src = TBL / "all_genes_features_unified_n57.tsv"

rows = list(csv.DictReader(open(src), delimiter="\t"))
def fnum(r, k):
    try: return float(r[k])
    except: return np.nan
y   = np.array([int(float(r["is_exposed"])) for r in rows])
mid = np.array([(fnum(r,"start")+fnum(r,"end"))/2 for r in rows])
dist= np.array([fnum(r,"nearest_methyl_distance") for r in rows])
# clean region: arm=1 if midpoint in chromosome arms, else core=0
region = np.where((mid <= ARM_LEFT_MAX) | (mid >= ARM_RIGHT_MIN), "arm", "core")

def auc(yy, xx):
    m = ~np.isnan(xx); yy, xx = yy[m], xx[m]
    n1=(yy==1).sum(); n0=(yy==0).sum()
    if n1==0 or n0==0: return np.nan, int(n1), int(n0)
    r = rankdata(xx); U = r[yy==1].sum()-n1*(n1+1)/2
    return U/(n1*n0), int(n1), int(n0)

out = []
# overall (smaller distance -> exposed => rank on -dist)
a, n1, n0 = auc(y, -dist)
out.append(("overall", a, n1, n0, np.nan))
for strat in ("core", "arm"):
    mk = region == strat
    a_s, e_s, s_s = auc(y[mk], -dist[mk])
    md = ~np.isnan(dist) & mk
    try:
        U, p = mannwhitneyu(dist[md & (y==1)], dist[md & (y==0)], alternative="less")
    except ValueError:
        p = np.nan
    out.append((strat, a_s, e_s, s_s, p))

# exposed enrichment by region
tab = [[int(((region=="arm")&(y==1)).sum()), int(((region=="arm")&(y==0)).sum())],
       [int(((region=="core")&(y==1)).sum()), int(((region=="core")&(y==0)).sum())]]
OR, p_fisher = fisher_exact(tab)

# logistic: is_exposed ~ log10(dist+1) + region(arm=1)
m = ~np.isnan(dist)
X = np.column_stack([np.log10(dist[m]+1), (region[m]=="arm").astype(float)])
X = sm.add_constant(X)
res = sm.Logit(y[m], X).fit(disp=0)
coef = res.params; pvals = res.pvalues

with open(TBL/"M1_geographic_robustness.tsv","w") as f:
    f.write("stratum\tAUC_distance\tn_exposed\tn_shielded\tMWU_p\n")
    for s,a_,e,sh,p in out:
        f.write(f"{s}\t{a_:.4f}\t{e}\t{sh}\t{'' if np.isnan(p) else f'{p:.2e}'}\n")
    f.write("\n# Exposed enrichment in arm vs core\n")
    f.write(f"# contingency [arm:exp,shi / core:exp,shi] = {tab}\n")
    f.write(f"# Fisher OR(arm vs core) = {OR:.3f}, p = {p_fisher:.3f}\n")
    f.write("\n# Logistic: is_exposed ~ log10(dist+1) + region(arm=1)\n")
    f.write(f"# beta_log10dist = {coef[1]:.3f} (p={pvals[1]:.2e})  <- distance effect ADJUSTED for geography\n")
    f.write(f"# beta_region    = {coef[2]:.3f} (p={pvals[2]:.2e})\n")

# figure: ROC per stratum
fig, ax = plt.subplots(1,2, figsize=(9,4))
for strat,c in (("core","#3C3489"),("arm","#D85A30")):
    mk = (region==strat) & (~np.isnan(dist))
    fpr,tpr,_ = roc_curve(y[mk], -dist[mk])
    a_s,_,_ = auc(y[mk], -dist[mk])
    ax[0].plot(fpr,tpr,color=c,lw=2,label=f"{strat} (AUC={a_s:.3f})")
ax[0].plot([0,1],[0,1],"--",color="#bbb")
ax[0].set_xlabel("FPR"); ax[0].set_ylabel("TPR"); ax[0].set_title("Partition holds within strata")
ax[0].legend(loc="lower right", fontsize=9)
for strat,c in (("core","#3C3489"),("arm","#D85A30")):
    mk=(region==strat)&(~np.isnan(dist))
    ax[1].hist(np.clip(dist[mk&(y==1)],0,3000),bins=20,alpha=.6,color=c,density=True,label=f"{strat} exposed")
ax[1].set_xlabel("nearest methylation distance (bp)"); ax[1].set_ylabel("density")
ax[1].set_title("Exposed = methylation-proximal in both strata"); ax[1].legend(fontsize=8)
plt.tight_layout()
fig.savefig(FIG/"M1_stratified_ROC.png", dpi=150)
fig.savefig(FIG/"M1_stratified_ROC.pdf")
print("OVERALL/stratified AUC and adjusted logistic written.")
for s,a_,e,sh,p in out: print(f"  {s}: AUC={a_:.4f} exp={e} shi={sh}")
print(f"  Fisher arm-vs-core OR={OR:.3f} p={p_fisher:.3f}")
print(f"  Logistic beta_log10dist={coef[1]:.3f} p={pvals[1]:.2e}; beta_region={coef[2]:.3f} p={pvals[2]:.2e}")
