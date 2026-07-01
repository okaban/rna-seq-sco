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

# ---- canonical muted palette (matches Fig7/Fig9 + figure-revision-spec_260623) ----
C_GCC = "#C26B6B"   # GCCGGC m4C system (muted red)
C_AAG = "#9970AB"   # AAGCCCG m4C/6mA system (muted purple)
C_ALL = "#BBBBBB"   # all regulators / background (grey)
CORE_LO, CORE_HI = 1_500_000, 7_170_000   # chromosomal core (single canonical definition)

# ============ FIGURE 4 — two-system marking ============
fig, ax = plt.subplots(1, 3, figsize=(12, 4.0))
# (a) the two systems mark POSITIONALLY DISTINCT regulator sets (core localisation + 2-gene overlap)
def core_frac(df, lo=CORE_LO, hi=CORE_HI):
    p = df['tss'].astype(float)
    return ((p >= lo) & (p <= hi)).mean() * 100
labels = [f"GCCGGC-\nExposed\n(n={len(E)})", f"AAGCCCG-\npromoter\n(n={len(AAGp)})", f"all\nregulators\n(n={len(reg)})"]
vals = [core_frac(E), core_frac(AAGp), core_frac(reg)]
cols = [C_GCC, C_AAG, C_ALL]
ax[0].bar(range(3), vals, color=cols, width=0.62)
for k, v in enumerate(vals):
    ax[0].text(k, v+1.5, f"{v:.0f}%", ha="center", fontsize=9)
ax[0].set_xticks(range(3)); ax[0].set_xticklabels(labels, fontsize=8)
ax[0].set_ylabel("% of class in chromosomal core"); ax[0].set_ylim(0, 100)
ax[0].set_title("(a) Distinct chromosomal positioning")
ax[0].annotate(f"overlap = {(reg.Exposed&reg.AAG).sum()}/{len(E)} genes\n(largely distinct sets)",
               xy=(0.5, 0.93), xycoords="axes fraction", ha="center", va="top",
               fontsize=8, style="italic", color="#333333")
# (b) TF family enrichment in Exposed vs rest (Fisher OR) — significant = GCCGGC colour
fams = ['MerR','LysR','LacI','TetR','Sigma factor','Sensor kinase']
ors, ps = [], []
for fam in fams:
    a=((reg.Exposed)&(reg.tf_family==fam)).sum(); c=((~reg.Exposed)&(reg.tf_family==fam)).sum()
    od,p=fisher_exact([[a,len(E)-a],[c,len(reg)-len(E)-c]]); ors.append(od); ps.append(p)
cols=[C_GCC if p<0.05 else C_ALL for p in ps]
ax[1].barh(range(len(fams)), ors, color=cols)
ax[1].axvline(1, color="k", ls="--", lw=0.8)
ax[1].set_yticks(range(len(fams))); ax[1].set_yticklabels(fams)
for i,(o,p) in enumerate(zip(ors,ps)):
    ax[1].text(o+0.1, i, f"OR={o:.1f}{'*' if p<0.05 else ''}", va="center", fontsize=8)
ax[1].set_xlim(0, max(ors)*1.25)
ax[1].set_xlabel("odds ratio (Exposed vs other regulators)")
ax[1].set_title("(b) GCCGGC-Exposed: family enrichment")
# (c) demethylation trajectory: fraction with promoter mark <=293bp at T1/T2/T3
def near_frac(df):
    return [ (df['dT1']<=W).mean(), (df['dT2']<=W).mean(), (df['dT3']<=W).mean() ]
ax[2].plot([1,2,3], np.array(near_frac(E))*100, '-o', color=C_GCC, label=f"GCCGGC-Exposed (n={len(E)})")
ax[2].plot([1,2,3], np.array(near_frac(AAGp))*100, '-s', color=C_AAG, label=f"AAGCCCG-prom (n={len(AAGp)})")
ax[2].set_xticks([1,2,3]); ax[2].set_xticklabels(['T1\n(12h)','T2\n(24h)','T3\n(50h)'])
ax[2].set_ylabel("% with promoter mark (≤293 bp)"); ax[2].set_ylim(-3,103)
ax[2].set_title("(c) Synchronous demethylation at T2"); ax[2].legend(frameon=False, fontsize=8)
fig.suptitle("Two methylation systems mark largely distinct regulator classes (vegetative growth)",
             y=1.02, fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(FIG/"Figure4_two_system_marking.pdf", bbox_inches="tight"); fig.savefig(FIG/"Figure4_two_system_marking.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ============ FIGURE 5 — synchronized demethylation + weak bias ============
fig, ax = plt.subplots(1, 3, figsize=(12, 3.8))
# (a) count of Exposed promoters still methylated at each timepoint
cnt = [ (E['dT1']<=W).sum(), (E['dT2']<=W).sum(), (E['dT3']<=W).sum() ]
ax[0].bar(['T1\n(12h)','T2\n(24h)','T3\n(50h)'], cnt, color=["#C26B6B","#BBBBBB","#BBBBBB"])
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
