"""G2: does methylation CHANGE (T1->T2) co-vary with the 3D refolding
(compartment PC and FIRE change M->L)? Tests the central 'mirrors refolding'
claim dynamically, not by static overlap. Seed 42. 2026-06-29."""
import csv, numpy as np
from scipy import stats
np.random.seed(42)
bn=list(csv.DictReader(open("76_FIRE_methylation_crossref/tables/fire_methyl_bins.tsv"),delimiter="\t"))
def col(k): return np.array([float(r[k]) for r in bn])
core=np.array([str(r["core"]).lower() in ("true","1","yes") for r in bn])
dmeth=col("m_n_T2")-col("m_n_T1")
dfire=col("FIRE_L")-col("FIRE_M")
dpc=col("PC_L")-col("PC_M")
def rep(lbl,x,y,m=None):
    if m is not None: x,y=x[m],y[m]
    r=stats.spearmanr(x,y); print(f"{lbl:46s} rho={r.correlation:+.3f}  p={r.pvalue:.2e}  n={len(x)}")
print("G2 dynamic co-change (per 5-kb window):\n")
rep("Δmeth(T1→T2) vs ΔFIRE(M→L) [genome]", dmeth, dfire)
rep("Δmeth(T1→T2) vs ΔPC/compartment(M→L) [genome]", dmeth, dpc)
rep("Δmeth vs ΔFIRE [within core]", dmeth, dfire, core)
rep("Δmeth vs ΔFIRE [within arm]", dmeth, dfire, ~core)
# directional check: do windows that LOSE methylation also lose FIRE/compartment?
lost=dmeth<0
print(f"\nwindows losing meth (n={lost.sum()}): median ΔFIRE={np.median(dfire[lost]):+.3f}; "
      f"gaining meth (n={(~lost & (dmeth>0)).sum()}): median ΔFIRE={np.median(dfire[dmeth>0]):+.3f}")
open("85_deng_integration_gap_analysis/tables/G2_dynamic.tsv","w").write(
 f"test\trho\tp\n"
 f"dmeth_vs_dFIRE_genome\t{stats.spearmanr(dmeth,dfire).correlation:.3f}\t{stats.spearmanr(dmeth,dfire).pvalue:.2e}\n"
 f"dmeth_vs_dPC_genome\t{stats.spearmanr(dmeth,dpc).correlation:.3f}\t{stats.spearmanr(dmeth,dpc).pvalue:.2e}\n")
