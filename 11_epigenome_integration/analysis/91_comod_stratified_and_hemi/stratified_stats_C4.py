#!/usr/bin/env python3
"""REFA-02(b,c): per-read AAGCCCG 4mC(C4)/6mA(A0/A1) coupling, T1, from tables/perread_pairs_T1_C4.tsv.gz.
(i) Mantel-Haenszel OR stratified by motif instance (statsmodels StratifiedTable, RGB CI);
(ii) per-replicate pooled ORs; (iii) fraction of instances with individual OR>1 among instances
with >=5 in every 2x2 margin; (c) within-read cross-instance null: 4mC at instance i vs 6mA at a
different instance j on the same read (all ordered pairs i!=j among instances covered by the read).
"""
import pandas as pd, numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.contingency_tables import StratifiedTable, Table2x2
d = pd.read_csv("tables/perread_pairs_T1_C4.tsv.gz", sep="\t")
d["inst"] = d.strand + d.inst_start.astype(str)
print("pairs", len(d), "instances", d.inst.nunique(), "reads", d.read.nunique(), d.rep.value_counts().to_dict())
def t2x2(a, c):
    a = np.asarray(a, bool); c = np.asarray(c, bool)
    return np.array([[np.sum(a & c), np.sum(a & ~c)], [np.sum(~a & c), np.sum(~a & ~c)]])
rows = []; inst_rows = []
for t in (0.5, 0.75, 0.9):
    a = d.pA >= t; c = d.pC >= t
    T = t2x2(a, c); OR, p = fisher_exact(T); ci = Table2x2(T).oddsratio_confint()
    rows.append(dict(threshold=t, stratum="pooled (79_ reproduction)", n_both=T[0,0], n_A_only=T[0,1], n_C_only=T[1,0], n_neither=T[1,1], OR=OR, CI_lo=ci[0], CI_hi=ci[1], p=p))
    # per replicate
    for rep, g in d.groupby("rep"):
        T = t2x2(g.pA >= t, g.pC >= t); OR, p = fisher_exact(T); ci = Table2x2(T).oddsratio_confint()
        rows.append(dict(threshold=t, stratum=f"replicate {rep}", n_both=T[0,0], n_A_only=T[0,1], n_C_only=T[1,0], n_neither=T[1,1], OR=OR, CI_lo=ci[0], CI_hi=ci[1], p=p))
    # MH by instance
    tabs = []; n_or1 = n_elig = 0; n_inf = 0
    for inst, g in d.groupby("inst"):
        T = t2x2(g.pA >= t, g.pC >= t)
        if T.sum() < 2: continue
        tabs.append(T)
        r1, r2 = T[0].sum(), T[1].sum(); c1, c2 = T[:,0].sum(), T[:,1].sum()
        if min(r1, r2, c1, c2) >= 5:
            n_elig += 1; o = (T[0,0]*T[1,1]) / (T[0,1]*T[1,0]) if T[0,1]*T[1,0] > 0 else np.inf
            n_or1 += o > 1
            inst_rows.append(dict(threshold=t, inst=inst, n=T.sum(), n_both=T[0,0], n_A_only=T[0,1], n_C_only=T[1,0], n_neither=T[1,1], OR=o))
    st = StratifiedTable([x for x in tabs]); mh = st.oddsratio_pooled; mhci = st.oddsratio_pooled_confint()
    tst = st.test_null_odds(); hom = st.test_equal_odds()
    rows.append(dict(threshold=t, stratum=f"Mantel-Haenszel by instance (k={len(tabs)} strata)", n_both=np.nan, n_A_only=np.nan, n_C_only=np.nan, n_neither=np.nan, OR=mh, CI_lo=mhci[0], CI_hi=mhci[1], p=tst.pvalue,
                     note=f"Breslow-Day homogeneity p={hom.pvalue:.3g}; instances with all margins>=5: {n_elig}; of which OR>1: {n_or1} ({n_or1/n_elig:.1%})"))
    # informative strata: instances with at least one A+ and one C+ pair
    # (c) cross-instance null within reads covering >=2 instances
    multi = d.groupby("read").filter(lambda g: len(g) >= 2)
    xs = []
    for read, g in multi.groupby("read"):
        A = (g.pA.values >= t); C = (g.pC.values >= t); n = len(g)
        for i in range(n):
            for j in range(n):
                if i != j: xs.append((C[i], A[j]))
    xs = np.array(xs); T = t2x2(xs[:,1], xs[:,0]); OR, p = fisher_exact(T); ci = Table2x2(T).oddsratio_confint()
    rows.append(dict(threshold=t, stratum=f"cross-instance same-read null (reads>=2 inst: {multi.read.nunique()}; ordered pairs)", n_both=T[0,0], n_A_only=T[0,1], n_C_only=T[1,0], n_neither=T[1,1], OR=OR, CI_lo=ci[0], CI_hi=ci[1], p=p))
    # same-instance OR restricted to the same multi-instance reads, for like-for-like comparison
    T = t2x2(multi.pA >= t, multi.pC >= t); OR, p = fisher_exact(T); ci = Table2x2(T).oddsratio_confint()
    rows.append(dict(threshold=t, stratum="same-instance OR restricted to reads>=2 inst", n_both=T[0,0], n_A_only=T[0,1], n_C_only=T[1,0], n_neither=T[1,1], OR=OR, CI_lo=ci[0], CI_hi=ci[1], p=p))
    # MH stratified by read (within-read) is the cleanest read-quality control: not identifiable for single-instance reads; skip.
R = pd.DataFrame(rows); R.to_csv("tables/REFA02b_perread_OR_stratified_T1_C4.tsv", sep="\t", index=False)
pd.DataFrame(inst_rows).to_csv("tables/REFA02b_per_instance_OR_margins_ge5.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); print(R[["threshold","stratum","n_both","n_A_only","n_C_only","n_neither","OR","CI_lo","CI_hi","p"]].round(3).to_string()); print(R.dropna(subset=["note"]).note.to_string())
