"""C18: is the methylome-FIRE / Exposed-FIRE signal just oriC proximity?
oriC ~ dnaA (SCO3879) midpoint = 4,271,763 bp on NC_003888.3. Seed 42. 2026-06-29."""
import csv, numpy as np
from scipy import stats
np.random.seed(42)
ORIC = 4271763.0
B="/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/76_FIRE_methylation_crossref/tables/"
OUT="/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/83_oriC_confound_FIRE/tables/oric_confound_results.tsv"

def sp(x,y): return stats.spearmanr(x,y)
def partial_spearman(x,y,z):
    rxy=stats.spearmanr(x,y).correlation; rxz=stats.spearmanr(x,z).correlation; ryz=stats.spearmanr(y,z).correlation
    return (rxy-rxz*ryz)/np.sqrt((1-rxz**2)*(1-ryz**2))

# windows
bn=list(csv.DictReader(open(B+"fire_methyl_bins.tsv"),delimiter="\t"))
fire=np.array([float(r["FIRE_M"]) for r in bn]); meth=np.array([float(r["m_n_T1"]) for r in bn])
dor=np.abs(np.array([float(r["center"]) for r in bn])-ORIC)
out=[]
r0=sp(meth,fire); out.append(("windows: Spearman(meth_T1, FIRE) raw", f"rho={r0.correlation:.3f}", f"p={r0.pvalue:.2e}", f"n={len(fire)}"))
pp=partial_spearman(meth,fire,dor); out.append(("windows: partial Spearman(meth,FIRE | dist_oriC)", f"rho={pp:.3f}", "", ""))
rf=sp(fire,dor); rm=sp(meth,dor)
out.append(("windows: FIRE vs dist_oriC (gradient)", f"rho={rf.correlation:.3f}", f"p={rf.pvalue:.2e}", ""))
out.append(("windows: meth_T1 vs dist_oriC (gradient)", f"rho={rm.correlation:.3f}", f"p={rm.pvalue:.2e}", ""))

# genes: Exposed vs Shielded FIRE, overall + stratified by oriC-distance quartile
gn=list(csv.DictReader(open(B+"regulatory_FIRE_at_TSS.tsv"),delimiter="\t"))
gn=[r for r in gn if r["FIRE_tss"] not in ("","NA")]
for r in gn: r["_d"]=abs(float(r["tss"])-ORIC); r["_f"]=float(r["FIRE_tss"]); r["_e"]=(str(r["Exposed"]).lower() in ("1","true","yes","exposed"))
E=[r["_f"] for r in gn if r["_e"]]; S=[r["_f"] for r in gn if not r["_e"]]
U,p=stats.mannwhitneyu(E,S,alternative="two-sided")
out.append(("genes: Exposed vs Shielded FIRE_tss (overall)", f"medE={np.median(E):.2f} medS={np.median(S):.2f}", f"MWU p={p:.2e}", f"nE={len(E)} nS={len(S)}"))
ds=np.array([r["_d"] for r in gn]); q=np.quantile(ds,[.25,.5,.75])
for i,(lo,hi) in enumerate([(0,q[0]),(q[0],q[1]),(q[1],q[2]),(q[2],1e12)]):
    sub=[r for r in gn if lo<=r["_d"]<hi]
    e=[r["_f"] for r in sub if r["_e"]]; s=[r["_f"] for r in sub if not r["_e"]]
    if len(e)>=3 and len(s)>=3:
        u,pp2=stats.mannwhitneyu(e,s,alternative="two-sided")
        out.append((f"genes: Exposed>Shielded FIRE within oriC-dist Q{i+1}", f"medE={np.median(e):.2f} medS={np.median(s):.2f}", f"MWU p={pp2:.2e}", f"nE={len(e)} nS={len(s)}"))
    else:
        out.append((f"genes: oriC-dist Q{i+1}", f"nE={len(e)} nS={len(s)} (too few)", "", ""))

with open(OUT,"w") as f:
    f.write("analysis\teffect\tp\tn\n")
    for a,b,c,d in out: f.write(f"{a}\t{b}\t{c}\t{d}\n")
print("oriC =",ORIC,"(dnaA SCO3879 midpoint)\n")
for a,b,c,d in out: print(f"{a:52s} {b:28s} {c:16s} {d}")
