"""G1: do methylation / Exposed promoters align with Deng CID boundaries?
Boundaries = midpoints of the 81 CID-boundary-located genes (sd02, M-phase).
Seed 42. 2026-06-29."""
import pandas as pd, numpy as np, csv
from scipy import stats
np.random.seed(42)
D="76_FIRE_methylation_crossref/"
# CID boundaries
b=pd.ExcelFile(D+"data/pnas.2222045120.sd02.xlsx").parse("M-phase",header=1)
b=b.dropna(subset=["Start(bp)","End(bp)"])
bnd=((b["Start(bp)"]+b["End(bp)"])/2).values
print(f"CID boundaries (M-phase): n={len(bnd)}")
def dmin(x): return np.min(np.abs(bnd-x))
# windows: methylation density vs distance to nearest boundary
bn=list(csv.DictReader(open(D+"tables/fire_methyl_bins.tsv"),delimiter="\t"))
cen=np.array([float(r["center"]) for r in bn]); meth=np.array([float(r["m_n_T1"]) for r in bn])
dwin=np.array([dmin(c) for c in cen])
# is methylation higher in windows NEAR a boundary (<=10kb) vs far?
near=dwin<=10000
U,p=stats.mannwhitneyu(meth[near],meth[~near],alternative="two-sided")
print(f"\n[windows] meth density near-boundary(<=10kb, n={near.sum()}) vs far(n={(~near).sum()}): "
      f"med {np.median(meth[near]):.1f} vs {np.median(meth[~near]):.1f}, MWU p={p:.2e}")
rho=stats.spearmanr(meth,dwin); print(f"[windows] Spearman(meth, dist_to_boundary) rho={rho.correlation:.3f} p={rho.pvalue:.2e}")
# genes: Exposed vs Shielded distance to nearest boundary
g=pd.read_csv("52_shielded_exposed_boundary/tables/SuppTable_regulatory_gene_classification_n1051.tsv",sep="\t")
g=g.dropna(subset=["tss"])
g["d_bnd"]=g["tss"].apply(dmin)
E=g[g["class"]=="Exposed"]["d_bnd"].values; S=g[g["class"]=="Shielded"]["d_bnd"].values
U,p=stats.mannwhitneyu(E,S,alternative="two-sided")
print(f"\n[genes] dist-to-CID-boundary Exposed(n={len(E)}) vs Shielded(n={len(S)}): "
      f"med {np.median(E):.0f} vs {np.median(S):.0f} bp, MWU p={p:.2e}")
# fraction within 20kb of a boundary
fE=np.mean(E<=20000); fS=np.mean(S<=20000)
from scipy.stats import fisher_exact
orr,pf=fisher_exact([[sum(E<=20000),sum(E>20000)],[sum(S<=20000),sum(S>20000)]])
print(f"[genes] frac within 20kb of boundary: Exposed {fE:.2f} vs Shielded {fS:.2f}, OR={orr:.2f} p={pf:.2e}")
open("85_deng_integration_gap_analysis/tables/G1_cid_boundary.tsv","w").write(
 f"metric\tvalue\nn_boundaries\t{len(bnd)}\nwin_meth_near_vs_far_p\t{p:.2e}\n"
 f"gene_dist_Exposed_med\t{np.median(E):.0f}\ngene_dist_Shielded_med\t{np.median(S):.0f}\n")
