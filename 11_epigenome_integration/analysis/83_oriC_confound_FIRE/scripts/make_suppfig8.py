"""Revised Supplementary Figure 8 — methylome-FIRE is compartment-scale only:
(a) genome-wide meth vs FIRE coloured core/arm; (b) within-core correlation null
across 5-100 kb; (c) developmental co-change Δmeth vs ΔFIRE (G2, rho=0.43).
Canonical muted palette. 2026-06-30."""
import csv, numpy as np
from collections import defaultdict
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
np.random.seed(42)
C_CORE="#C26B6B"; C_ARM="#BBBBBB"; C_ACC="#4477AA"
B="/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/76_FIRE_methylation_crossref/tables/fire_methyl_bins.tsv"
bn=list(csv.DictReader(open(B),delimiter="\t"))
def col(k): return np.array([float(r[k]) for r in bn])
core=np.array([str(r["core"]).lower() in ("true","1","yes") for r in bn])
fireM=col("FIRE_M"); methT1=col("m_n_T1"); cen=col("center")
dmeth=col("m_n_T2")-methT1; dfire=col("FIRE_L")-fireM

fig,ax=plt.subplots(1,3,figsize=(13,4.0))
plt.rcParams.update({"axes.spines.top":False,"axes.spines.right":False})
# (a) genome-wide scatter coloured core/arm
ax[0].scatter(methT1[~core],fireM[~core],s=8,c=C_ARM,alpha=0.5,edgecolors="none",label="arm")
ax[0].scatter(methT1[core],fireM[core],s=8,c=C_CORE,alpha=0.5,edgecolors="none",label="core")
rho=stats.spearmanr(methT1,fireM).correlation
ax[0].set_xlabel("T1 GCCGGC m4C sites / 5 kb"); ax[0].set_ylabel("FIRE (M-phase)")
ax[0].set_title(f"(a) Genome-wide co-localisation\nSpearman ρ = {rho:.2f} (= core/arm contrast)",fontsize=10)
ax[0].legend(frameon=False,fontsize=8,markerscale=2)
# (b) within-core correlation by resolution (null)
def agg_core_rho(res):
    g=defaultdict(lambda:[0.0,[],0])
    for c,f,m,cr in zip(cen,fireM,methT1,core):
        k=int(c//res); g[k][0]+=m; g[k][1].append(f); g[k][2]+=cr
    ms=[];fs=[]
    for k,(msum,fl,crsum) in g.items():
        if crsum>len(fl)/2: ms.append(msum); fs.append(np.mean(fl))
    if len(ms)<8: return np.nan,np.nan,len(ms)
    r=stats.spearmanr(ms,fs); return r.correlation,r.pvalue,len(ms)
res_kb=[5,20,50,100]; rhos=[];ps=[]
for r in res_kb:
    rr,pp,nn=agg_core_rho(r*1000); rhos.append(rr); ps.append(pp)
bars=ax[1].bar(range(len(res_kb)),rhos,color=C_ACC,width=0.6)
for i,(rr,pp) in enumerate(zip(rhos,ps)):
    ax[1].text(i,rr+0.01,f"ρ={rr:.2f}\n(ns)" if pp>=0.05 else f"ρ={rr:.2f}\np={pp:.0e}",ha="center",va="bottom",fontsize=7.5)
ax[1].axhline(0,color="k",lw=0.6)
ax[1].set_xticks(range(len(res_kb))); ax[1].set_xticklabels([f"{r} kb" for r in res_kb])
ax[1].set_ylim(-0.05,0.45); ax[1].set_xlabel("aggregation window")
ax[1].set_ylabel("within-core Spearman ρ (meth vs FIRE)")
ax[1].set_title("(b) No within-compartment resolution\n(all p > 0.1, 5–100 kb)",fontsize=10)
# (c) dynamic co-change Δmeth vs ΔFIRE
ax[2].scatter(dmeth,dfire,s=8,c=C_CORE,alpha=0.4,edgecolors="none")
rd=stats.spearmanr(dmeth,dfire)
ax[2].axhline(0,color="k",lw=0.5); ax[2].axvline(0,color="k",lw=0.5)
ax[2].set_xlabel("Δ methylation (T1→T2), sites/5 kb"); ax[2].set_ylabel("Δ FIRE (M→L)")
ax[2].set_title(f"(c) Developmental co-change\nSpearman ρ = {rd.correlation:.2f} (p = {rd.pvalue:.0e})",fontsize=10)
fig.suptitle("Supplementary Figure 8 | Methylome–FIRE relationship is compartment-scale and developmentally co-dynamic",
             y=1.03,fontsize=11,fontweight="bold")
fig.tight_layout()
O="/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/83_oriC_confound_FIRE/figures"
import os; os.makedirs(O,exist_ok=True)
fig.savefig(O+"/SuppFig8_compartment_scale.png",dpi=150,bbox_inches="tight")
fig.savefig(O+"/SuppFig8_compartment_scale.pdf",bbox_inches="tight"); plt.close(fig)
print("SuppFig8 done. within-core rho by res:",list(zip(res_kb,[round(r,3) for r in rhos])),"| genome rho",round(rho,3),"| dynamic rho",round(rd.correlation,3))
