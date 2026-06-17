"""Regenerate main Figures 4 and 5 under the LOCKED reframe (2026-06-15).

F4 = two methylation systems mark distinct regulator classes (vegetative):
     GCCGGC-Exposed (62) vs AAGCCCG-promoter (22) by position, TF family, demethylation.
F5 = the 62 Exposed promoters are synchronously demethylated at T2 and bias
     expression only weakly (geography-controlled r=-0.09); bldD exemplar.

Non-circular definition: Exposed = regulatory gene with TSS<=293bp GCCGGC m4C at T1.
region: 0=core, 1=arm.
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact, mannwhitneyu, pearsonr
from pathlib import Path

B = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis")
A = B/"52_shielded_exposed_boundary"
FIG = A/"figures"; FIG.mkdir(exist_ok=True)
W = 293
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})

reg = pd.read_csv(A/"tables/all_genes_features_unified_n57.tsv", sep="\t").dropna(subset=['tss']).copy()
reg['tss'] = reg['tss'].astype(int)
g = pd.read_csv(B/"37_defense_island_GCCGGC/tables/GCCGGC_sites_by_timepoint.tsv", sep="\t")
aag = pd.read_csv(B/"36_AAGCCCG_distribution/tables/AAGCCCG_site_gene_mapping.tsv", sep="\t")

def nd(tss, pos):
    if len(pos) == 0: return np.nan
    i = np.clip(np.searchsorted(pos, tss), 1, len(pos)-1); return min(abs(tss-pos[i-1]), abs(tss-pos[i]))
for tp in ['T1','T2','T3']:
    pos = np.sort(g[g.timepoint==tp].position.values); reg['d'+tp] = reg['tss'].apply(lambda t: nd(t, pos))
reg['Exposed'] = reg['dT1'] <= W
aagprom = set(aag[(aag.timepoint=='T1') & (aag['distance']<=W)]['locus_tag'])
reg['AAG'] = reg['locus_tag'].isin(aagprom)
E = reg[reg.Exposed]; AAGp = reg[reg.AAG]
print(f"GCCGGC-Exposed={len(E)}  AAGCCCG-promoter={len(AAGp)}  overlap={(reg.Exposed&reg.AAG).sum()}")
print(f"Exposed core frac={(E.region==0).mean():.2f}  AAGp core frac={(AAGp.region==0).mean():.2f}")

# ============ FIGURE 4 — two-system marking ============
fig, ax = plt.subplots(1, 3, figsize=(12, 3.8))
# (a) core-enrichment is robust across boundary definitions (not an arbitrary cut)
#     three independent central-region definitions, applied to the same gene sets:
#       - Compartment A (Deng 2023 Hi-C):           2.30-6.20 Mb
#       - chromosomal core (conserved-gene span):   1.50-7.17 Mb
#       - geometry-only central half of replicon:   25-75% of 8.67 Mb (2.17-6.50 Mb)
GENOME = 8_667_507
CUTS = [("Compartment A\n(2.3-6.2 Mb)", 2_300_000, 6_200_000),
        ("Chromosomal core\n(1.5-7.17 Mb)", 1_500_000, 7_170_000),
        ("Central half\n(geometry only)", 0.25*GENOME, 0.75*GENOME)]
def core_frac(df, lo, hi):
    p = df['tss'].astype(float)
    return ((p >= lo) & (p <= hi)).mean() * 100
CLASSES = [("GCCGGC-Exposed", E, "#C26B6B"),
           ("AAGCCCG-promoter", AAGp, "#4477AA"),
           ("all regulators", reg, "#BBBBBB")]
xc = np.arange(len(CUTS)); bw = 0.26
for j, (cname, df, col) in enumerate(CLASSES):
    vals = [core_frac(df, lo, hi) for _, lo, hi in CUTS]
    ax[0].bar(xc + (j-1)*bw, vals, bw, color=col, label=f"{cname} (n={len(df)})")
    for k, v in enumerate(vals):
        ax[0].text(xc[k] + (j-1)*bw, v+1.5, f"{v:.0f}", ha="center", fontsize=7)
ax[0].set_xticks(xc); ax[0].set_xticklabels([c[0] for c in CUTS], fontsize=8)
ax[0].set_ylabel("% of class in central region"); ax[0].set_ylim(0, 100)
ax[0].set_title("(a) Core-enrichment is robust to the boundary cut")
ax[0].legend(frameon=False, fontsize=7, loc="upper right")
# (b) TF family enrichment in Exposed vs rest (Fisher OR)
fams = ['MerR','LysR','LacI','TetR','Sigma factor','Sensor kinase']
ors, ps = [], []
for fam in fams:
    a=((reg.Exposed)&(reg.tf_family==fam)).sum(); c=((~reg.Exposed)&(reg.tf_family==fam)).sum()
    od,p=fisher_exact([[a,len(E)-a],[c,len(reg)-len(E)-c]]); ors.append(od); ps.append(p)
cols=["#009E73" if p<0.05 else "#999999" for p in ps]
ax[1].barh(range(len(fams)), ors, color=cols)
ax[1].axvline(1, color="k", ls="--", lw=0.8)
ax[1].set_yticks(range(len(fams))); ax[1].set_yticklabels(fams)
for i,(o,p) in enumerate(zip(ors,ps)):
    ax[1].text(o+0.1, i, f"OR={o:.1f}{'*' if p<0.05 else ''}", va="center", fontsize=8)
ax[1].set_xlabel("odds ratio (Exposed vs other regulators)")
ax[1].set_title("(b) TF family enrichment (GCCGGC-Exposed)")
# (c) demethylation trajectory: fraction with promoter mark <=293bp at T1/T2/T3
def near_frac(df):
    return [ (df['dT1']<=W).mean(), (df['dT2']<=W).mean(), (df['dT3']<=W).mean() ]
ax[2].plot([1,2,3], np.array(near_frac(E))*100, '-o', color="#D55E00", label=f"GCCGGC-Exposed (n={len(E)})")
ax[2].plot([1,2,3], np.array(near_frac(AAGp))*100, '-s', color="#CC79A7", label=f"AAGCCCG-prom (n={len(AAGp)})")
ax[2].set_xticks([1,2,3]); ax[2].set_xticklabels(['T1\n(12h)','T2\n(24h)','T3\n(50h)'])
ax[2].set_ylabel("% with promoter mark (≤293 bp)"); ax[2].set_ylim(-3,103)
ax[2].set_title("(c) Synchronous demethylation"); ax[2].legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(FIG/"Figure4_two_system_marking.pdf"); fig.savefig(FIG/"Figure4_two_system_marking.png", dpi=150); plt.close(fig)

# ============ FIGURE 5 — synchronized demethylation + weak bias ============
fig, ax = plt.subplots(1, 3, figsize=(12, 3.8))
# (a) count of Exposed promoters still methylated at each timepoint
cnt = [ (E['dT1']<=W).sum(), (E['dT2']<=W).sum(), (E['dT3']<=W).sum() ]
ax[0].bar(['T1\n(12h)','T2\n(24h)','T3\n(50h)'], cnt, color=["#D55E00","#BBBBBB","#BBBBBB"])
for i,c in enumerate(cnt): ax[0].text(i, c+0.8, str(int(c)), ha="center", fontsize=10)
ax[0].set_ylabel("Exposed promoters methylated (≤293 bp)")
ax[0].set_title(f"(a) Synchronous erasure ({len(E)}→0 at T2)")
# (b) weak bias: promoter GCCGGC occupancy (+-2kb, T1) vs LFC_T2, region-controlled
gT1 = g[g.timepoint=='T1']; P=np.sort(gT1.position.values); Fq=gT1.sort_values('position').frequency.values
def occ(tss,w):
    m=(P>=tss-w)&(P<=tss+w); return Fq[m].sum() if m.any() else 0.0
reg['occ2k']=reg['tss'].apply(lambda t:occ(t,2000))
d=reg[['occ2k','LFC_T2vsT1','region']].dropna()
rx=d.occ2k.rank()-np.polyval(np.polyfit(d.region.rank(),d.occ2k.rank(),1),d.region.rank())
ry=d.LFC_T2vsT1.rank()-np.polyval(np.polyfit(d.region.rank(),d.LFC_T2vsT1.rank(),1),d.region.rank())
rr,pp=pearsonr(rx,ry)
ax[1].scatter(d.occ2k, d.LFC_T2vsT1, s=8, alpha=0.4, c="#555555", edgecolors="none")
ax[1].set_xlabel("promoter GCCGGC occupancy (±2 kb, T1)"); ax[1].set_ylabel("log2 FC (T2 vs T1)")
ax[1].set_title(f"(b) Weak modulatory bias\nregion-controlled r={rr:.2f} (p={pp:.3f})")
ax[1].axhline(0, color="k", lw=0.5)
# (c) bidirectional exemplars: same promoter mark, opposite outcomes (permissive)
def lab(r):
    n = str(r.get('gene_name','')).strip()
    return n if n and n.lower()!='nan' else (str(r.get('old_locus_tag','')).strip() or r['locus_tag'])
ex = E.dropna(subset=['LFC_T2vsT1']).copy(); ex['name'] = ex.apply(lab, axis=1)
up = ex.sort_values('LFC_T2vsT1', ascending=False).head(4)
dn = ex.sort_values('LFC_T2vsT1').head(4)
sel = pd.concat([up, dn]).drop_duplicates('locus_tag').sort_values('LFC_T2vsT1')
colors = ["#D55E00" if v < 0 else "#009E73" for v in sel['LFC_T2vsT1']]  # Okabe-Ito: vermillion=repressed, bluish-green=induced
ax[2].barh(range(len(sel)), sel['LFC_T2vsT1'].values, color=colors)
ax[2].set_yticks(range(len(sel))); ax[2].set_yticklabels(sel['name'], fontsize=8)
ax[2].axvline(0, color="k", lw=0.6)
ax[2].set_xlabel("log2 FC (T2 vs T1)")
ax[2].set_title("(c) Same mark, opposite outcomes\n(Exposed regulators; permissive)")
fig.suptitle("Methylation marks the Exposed regulators but does NOT direct their expression (permissive)", y=1.04, fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(FIG/"Figure5_synchronized_demethylation.pdf", bbox_inches="tight"); fig.savefig(FIG/"Figure5_synchronized_demethylation.png", dpi=150); plt.close(fig)

print("F4/F5 saved. weak-bias region-controlled r=%.3f p=%.4f" % (rr,pp))
print("MerR/LysR/LacI ORs:", {f:round(o,1) for f,o in zip(fams,ors)})
print("bldD LFC_T2:", reg.loc[reg.locus_tag=='SC_RS09420','LFC_T2vsT1'].values)
