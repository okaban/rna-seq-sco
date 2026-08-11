"""C18 final: does methylome predict FIRE WITHIN the active core at the CLAIMED
20-100 kb resolution? Aggregate 5-kb windows; Spearman + top-quintile concordance
(precision/recall) within core, per resolution. oriC=4,271,763. Seed 42. 2026-06-29."""
import csv, numpy as np
from collections import defaultdict
from scipy import stats
np.random.seed(42)
B="/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/76_FIRE_methylation_crossref/tables/fire_methyl_bins.tsv"
bn=list(csv.DictReader(open(B),delimiter="\t"))
def f(r,k):
    try:return float(r[k])
    except:return None
rows=[(f(r,"center"),f(r,"FIRE_M"),f(r,"m_n_T1"),str(r["core"]).lower() in ("true","1","yes")) for r in bn]
rows=[r for r in rows if None not in r[:3]]

def agg(res):
    g=defaultdict(lambda:[0.0,[],[]])  # key -> [meth_sum, fire_list, core_list]
    for c,fire,meth,core in rows:
        k=int(c//res); g[k][0]+=meth; g[k][1].append(fire); g[k][2].append(core)
    out=[]
    for k,(ms,fl,cl) in g.items():
        out.append((ms, np.mean(fl), sum(cl)>len(cl)/2))  # meth_sum, mean FIRE, majority-core
    return out

def concordance(meth, fire):
    # precision/recall of top-quintile methylation predicting top-quintile FIRE
    mt=np.quantile(meth,0.8); ft=np.quantile(fire,0.8)
    pm=meth>=mt; pf=fire>=ft
    tp=np.sum(pm&pf)
    prec=tp/np.sum(pm) if np.sum(pm) else float('nan')
    rec=tp/np.sum(pf) if np.sum(pf) else float('nan')
    return prec,rec

print("oriC=4,271,763 (dnaA). Within-CORE methylome–FIRE by resolution:\n")
print(f"{'res':>7} {'scope':12} {'n':>5} {'Spearman rho':>13} {'p':>10} {'prec':>6} {'rec':>6}")
for res in [5000,20000,50000,100000]:
    a=agg(res)
    for scope,sub in [("within-core",[x for x in a if x[2]]),("genome-wide",a)]:
        meth=np.array([x[0] for x in sub]); fire=np.array([x[1] for x in sub])
        if len(meth)<8: 
            print(f"{res//1000:>5}kb {scope:12} {len(meth):>5}  (too few)"); continue
        rr=stats.spearmanr(meth,fire); pr,rc=concordance(meth,fire)
        print(f"{res//1000:>5}kb {scope:12} {len(meth):>5} {rr.correlation:>13.3f} {rr.pvalue:>10.1e} {pr:>6.2f} {rc:>6.2f}")
