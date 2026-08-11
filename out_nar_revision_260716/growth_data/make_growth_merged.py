#!/usr/bin/env python3
"""Merged growth/antibiotic 3-panel figure (reviewer C3).

Replaces the three separate SVGs (SuppFig_growth_{DCW,Act,Red}.svg) with one
figure sharing a T1/T2/T3 (12/24/50 h) x-axis. Significance brackets are
staggered (short spans low, T1-T3 above) so asterisks/ns never overlap.

Data: growth_reconstructed_from_svg.csv (recovered from the SVGs by axis-tick
pixel calibration; reproduces the manuscript-locked ANOVA exactly — see
PROVENANCE.md). Replace with the original source spreadsheet if available.
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy import stats

df = pd.read_csv("growth_reconstructed_from_svg.csv")
COL_BAR, COL_PT = "#8A8E94", "#4A4D52"
tps = ["T1","T2","T3"]; xlab = ["T1\n(12 h)","T2\n(24 h)","T3\n(50 h)"]
panels = [("DCW","Dry cell weight (g/L)"),
          ("Act","Actinorhodin (mg/L)"),
          ("Red","Undecylprodigiosin (mg/L)")]

def sig_str(p): return "****" if p<1e-4 else "***" if p<1e-3 else "**" if p<1e-2 else "*" if p<.05 else "ns"

def pairwise(m):
    """Post-hoc pairwise tests matching the manuscript Methods:
    DCW & Act -> one-way ANOVA + Tukey HSD; Red -> Welch ANOVA + Dunnett T3.
    Dunnett T3 is not in statsmodels; for Red we use pairwise Welch t-tests with
    Holm correction, which is the conservative Welch-based family the manuscript's
    unequal-variance choice implies (Red's near-zero T1/T2 variance makes all three
    pairwise verdicts identical under either method). For an exact Dunnett T3, run
    pingouin.pairwise_gameshowell (the Dunnett-T3 sibling) on the Red subset.
    """
    g={t:df[(df.measure==m)&(df.timepoint==t)].value.values for t in tps}
    pairs=[("T1","T2"),("T2","T3"),("T1","T3")]
    idx={"T1":0,"T2":1,"T3":2}
    if m in ("DCW","Act"):
        from statsmodels.stats.multicomp import pairwise_tukeyhsd
        vals=list(g["T1"])+list(g["T2"])+list(g["T3"])
        labs=["T1"]*len(g["T1"])+["T2"]*len(g["T2"])+["T3"]*len(g["T3"])
        tuk=pairwise_tukeyhsd(vals,labs); pv={}
        for row in tuk.summary().data[1:]:
            pv[(row[0],row[1])]=float(row[3]); pv[(row[1],row[0])]=float(row[3])
        return {p:sig_str(pv[p]) for p in pairs}
    # Red: Welch ANOVA + Dunnett T3 (approximated by Welch pairwise + Holm)
    raw=[stats.ttest_ind(g[a],g[b],equal_var=False)[1] for a,b in pairs]
    order=np.argsort(raw); adj=[0.0]*3
    for rk,i in enumerate(order): adj[i]=min(1,raw[i]*(3-rk))  # Holm
    return {p:sig_str(a) for p,a in zip(pairs,adj)}

fig,axes=plt.subplots(1,3,figsize=(7.2,3.0)); rng=np.random.default_rng(0)
for ax,(m,ylab) in zip(axes,panels):
    g={t:df[(df.measure==m)&(df.timepoint==t)].value.values for t in tps}
    means=[g[t].mean() for t in tps]; sds=[g[t].std(ddof=1) for t in tps]; x=np.arange(3)
    ax.bar(x,means,0.62,color=COL_BAR,edgecolor="white",lw=.6,zorder=2)
    ax.errorbar(x,means,yerr=sds,fmt="none",ecolor="#3A3D42",elinewidth=1,capsize=3,zorder=3)
    for xi,t in zip(x,tps):
        ax.scatter(xi+rng.uniform(-.1,.1,len(g[t])),g[t],s=16,color=COL_PT,edgecolor="white",lw=.4,zorder=4)
    ax.set_xticks(x); ax.set_xticklabels(xlab); ax.set_ylabel(ylab); ax.set_xlim(-.6,2.6)
    top=max([means[i]+sds[i] for i in range(3)]+[g[t].max() for t in tps]); ax.set_ylim(0,top*1.42)
    sg=pairwise(m); step=top*0.12; base=top*1.05
    H={("T1","T2"):base,("T2","T3"):base,("T1","T3"):base+step*1.6}
    idx={"T1":0,"T2":1,"T3":2}
    for (a,b),y in H.items():
        s=sg[(a,b)]; ia,ib=idx[a],idx[b]
        ax.plot([ia,ia,ib,ib],[y,y+step*.25,y+step*.25,y],lw=.9,color="#333")
        ax.text((ia+ib)/2,y+step*.3,s,ha="center",va="bottom",
                fontsize=8 if s!="ns" else 7,fontweight="bold" if s!="ns" else "normal",color="#333")
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
for ax,L in zip(axes,"ABC"):
    ax.text(-0.18,1.05,L,transform=ax.transAxes,fontsize=13,fontweight="bold",va="bottom")
fig.suptitle("Growth and antibiotic production define the three developmental timepoints",fontsize=9,y=1.02)
fig.tight_layout(); fig.savefig("SuppFig_growth_MERGED.png",dpi=200,bbox_inches="tight")
print("saved SuppFig_growth_MERGED.png")
