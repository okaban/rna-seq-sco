"""C18 (corrected): test methylome-FIRE WITHIN the active compartment, not by
over-controlling genome-wide oriC-distance (which == the compartment axis Deng
identifies as functional). oriC=4,271,763 (dnaA). Seed 42. 2026-06-29."""
import csv, numpy as np
from scipy import stats
np.random.seed(42)
ORIC=4271763.0
B="/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/76_FIRE_methylation_crossref/tables/"
bn=list(csv.DictReader(open(B+"fire_methyl_bins.tsv"),delimiter="\t"))
def col(r,k): 
    try:return float(r[k])
    except:return None
fire=np.array([col(r,"FIRE_M") for r in bn]); meth=np.array([col(r,"m_n_T1") for r in bn])
core=np.array([str(r["core"]).lower() in ("true","1","yes") for r in bn])
cen=np.array([col(r,"center") for r in bn]); dor=np.abs(cen-ORIC)
out=[]
def sr(x,y,m=None):
    if m is not None: x,y=x[m],y[m]
    r=stats.spearmanr(x,y); return r.correlation,r.pvalue,len(x)
# genome-wide (for reference)
r,p,n=sr(meth,fire); out.append(("genome-wide meth–FIRE",f"rho={r:.3f}",f"p={p:.1e}",f"n={n}"))
# WITHIN CORE (active compartment) — the appropriate test
r,p,n=sr(meth,fire,core); out.append(("WITHIN-core meth–FIRE",f"rho={r:.3f}",f"p={p:.1e}",f"n={n}"))
# WITHIN ARM (control: do they also co-vary off the central compartment?)
r,p,n=sr(meth,fire,~core); out.append(("within-arm meth–FIRE",f"rho={r:.3f}",f"p={p:.1e}",f"n={n}"))
# within-core, is FIRE still just oriC distance? (if not deterministic within core, the link is real)
r,p,n=sr(fire[core],dor[core]); out.append(("within-core FIRE vs oriC-dist",f"rho={r:.3f}",f"p={p:.1e}",f"n={n}"))
# within-core partial meth-FIRE controlling oriC-dist (now legitimate: within compartment)
def psp(x,y,z):
    rxy=stats.spearmanr(x,y).correlation;rxz=stats.spearmanr(x,z).correlation;ryz=stats.spearmanr(y,z).correlation
    return (rxy-rxz*ryz)/np.sqrt((1-rxz**2)*(1-ryz**2))
pp=psp(meth[core],fire[core],dor[core]); out.append(("within-core partial meth–FIRE | oriC-dist",f"rho={pp:.3f}","",""))
print("oriC=4,271,763 (dnaA SCO3879)\n")
for a,b,c,d in out: print(f"{a:42s} {b:14s} {c:12s} {d}")
open("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/83_oriC_confound_FIRE/tables/within_compartment.tsv","w").write(
  "analysis\teffect\tp\tn\n"+"\n".join(f"{a}\t{b}\t{c}\t{d}" for a,b,c,d in out)+"\n")
