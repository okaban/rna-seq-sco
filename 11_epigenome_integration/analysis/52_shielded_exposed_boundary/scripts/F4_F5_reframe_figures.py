"""Regenerate main Figures 4 and 5 under the LOCKED reframe (2026-06-15).

F4 = two methylation systems mark distinct regulator classes (vegetative):
     GCCGGC-Exposed (62) vs AAGCCCG-promoter (22) by position, TF family, demethylation.
F5 = the 62 Exposed promoters are synchronously demethylated at T2 and bias
     expression only weakly (geography-controlled r=-0.09); bldD exemplar.

Non-circular definition: Exposed = regulatory gene with TSS<=293bp GCCGGC 4mC at T1.
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
# NAR full-width figures: 174 mm = 6.85 in. Fonts sized for the compressed panels.
plt.rcParams.update({"font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7.5,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6.5,
                     "axes.spines.top": False, "axes.spines.right": False})
FULL_W = 6.85   # 174 mm
NAR_FULL_PX = 2055   # 174 mm at 300 dpi

def _panel_letter(ax, lab, fontsize=13):
    """Uniform panel letter: anchored in the figure's outer margin at the panel's
    top-left, just left of the widest left-side text (y-label/ticks), so every
    panel gets an identically-sized, identically-placed bold lowercase letter
    regardless of y-label width. Mirrors 00_shared_utils.add_panel_label."""
    fig = ax.figure
    try:
        fig.canvas.draw()
        rend = fig.canvas.get_renderer()
        tb = ax.get_tightbbox(rend).transformed(fig.transFigure.inverted())
        xf = max(0.002, tb.x0 - 0.006)
        # Sit just above the axes top, but stay clear of the suptitle band (these
        # compact 3-panel figures reserve the top 10% for the suptitle).
        # 2026-09-13: anchor ABOVE the panel's tight bbox (which includes its
        # left-aligned title) so the letter never prints over the title text.
        # Common baseline: the tallest panel's tight-bbox top, so a/b/c align.
        top = max(a_.get_tightbbox(rend).transformed(fig.transFigure.inverted()).y1
                  for a_ in fig.axes)
        yf = min(top + 0.004, 0.975)
        fig.text(xf, yf, lab.lower(), fontsize=fontsize, fontweight="bold",
                 va="bottom", ha="left")
    except Exception:
        ax.annotate(lab.lower(), xy=(0, 1), xycoords="axes fraction",
                    xytext=(-26, 12), textcoords="offset points",
                    fontsize=fontsize, fontweight="bold", va="bottom",
                    ha="left", annotation_clip=False)


def _clamp_png_width(png_path, target_px=NAR_FULL_PX):
    """bbox_inches='tight' re-expands past the declared figsize; downscale the
    saved PNG to the NAR full-width pixel ceiling (aspect preserved) if it came
    out wider. Never upscales, so no content is ever clipped."""
    from PIL import Image
    im = Image.open(png_path)
    if im.size[0] > target_px:
        h = round(im.size[1] * target_px / im.size[0])
        im.resize((target_px, h), Image.LANCZOS).save(png_path, dpi=(300, 300))

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

# ---- unified palette "Rich & Calm" (2026-07, matches all main + supp figures) ----
C_GCC = "#A64B44"   # GCCGGC 4mC system → unified rich muted red (4mC)
C_AAG = "#6E5495"   # AAGCCCG dual 4mC/6mA system → unified muted purple (dual)
C_ALL = "#9AA7B0"   # all regulators / background → unified neutral grey
C_BAR = "#293039"   # single-category / count bars → unified grayscale dark (SuppFig15)
C_INDUCED = "#3E7256"   # induced (LFC>0) → unified deep muted green
C_REPRESSED = "#A64B44" # repressed (LFC<0) → unified rich muted red
CORE_LO, CORE_HI = 1_500_000, 7_170_000   # chromosomal core (single canonical definition)

# ============ FIGURE 4 — two-system marking ============
fig, ax = plt.subplots(1, 3, figsize=(FULL_W, 2.7))
# (a) the two systems mark POSITIONALLY DISTINCT regulator sets (core localisation + 2-gene overlap)
def core_frac(df, lo=CORE_LO, hi=CORE_HI):
    p = df['tss'].astype(float)
    return ((p >= lo) & (p <= hi)).mean() * 100
labels = [f"GCCGGC-\nExposed\n(n={len(E)})", f"AAGCCCG-\nprom.\n(n={len(AAGp)})", f"all\nreg.\n(n={len(reg)})"]
vals = [core_frac(E), core_frac(AAGp), core_frac(reg)]
cols = [C_GCC, C_AAG, C_ALL]
ax[0].bar(range(3), vals, color=cols, width=0.6)
for k, v in enumerate(vals):
    ax[0].text(k, v+2, f"{v:.0f}%", ha="center", fontsize=7, fontweight="bold")
ax[0].set_xticks(range(3)); ax[0].set_xticklabels(labels, fontsize=6)
ax[0].set_ylabel("% of class in chromosomal core"); ax[0].set_ylim(0, 108)
ax[0].set_title("Distinct chromosomal positioning", fontsize=8, pad=6, loc="left")
# overlap note in the free upper-right corner, clear of the bars/percent labels
ax[0].annotate(f"overlap = {(reg.Exposed&reg.AAG).sum()}/{len(E)}\n(largely distinct)",
               xy=(0.97, 0.97), xycoords="axes fraction", ha="right", va="top",
               fontsize=6, style="italic", color="#333333")
# (b) TF family enrichment in Exposed vs rest (Fisher OR) — significant = GCCGGC colour
fams = ['MerR','LysR','LacI','TetR','Sigma factor','Sensor kinase']
ors, ps = [], []
for fam in fams:
    a=((reg.Exposed)&(reg.tf_family==fam)).sum(); c=((~reg.Exposed)&(reg.tf_family==fam)).sum()
    od,p=fisher_exact([[a,len(E)-a],[c,len(reg)-len(E)-c]]); ors.append(od); ps.append(p)
# 2026-09-13 (FIG-07): the legend/body report Benjamini-Hochberg-corrected p across
# the six families (MerR/LysR significant; LacI p_BH = 0.10 not). The asterisk and
# the bar colour now follow the SAME BH rule instead of the raw Fisher p.
from statsmodels.stats.multitest import multipletests
p_bh = multipletests(ps, method="fdr_bh")[1]
print("TF-family Fisher:", {f:(round(o,2), round(p,4), round(q,4)) for f,o,p,q in zip(fams,ors,ps,p_bh)})
cols=[C_GCC if q<0.05 else C_ALL for q in p_bh]
ax[1].barh(range(len(fams)), ors, color=cols)
ax[1].axvline(1, color="k", ls="--", lw=0.8)
ax[1].set_yticks(range(len(fams))); ax[1].set_yticklabels(fams)
for i,(o,q) in enumerate(zip(ors,p_bh)):
    ax[1].text(o+max(ors)*0.03, i, f"OR={o:.1f}{'*' if q<0.05 else ''}",
               va="center", fontsize=6, color="#333")
ax[1].text(0.98, 0.98, "* $p_{BH}$ < 0.05", transform=ax[1].transAxes,
           ha="right", va="top", fontsize=5.5, color="#333")
# extra right headroom so the 'OR=4.7*' label on the longest bar stays inside the axis
ax[1].set_xlim(0, max(ors)*1.45)
ax[1].set_xlabel("odds ratio (Exposed vs other regulators)")
ax[1].set_title("Exposed: TF-family enrichment", fontsize=8, pad=6, loc="left")
# (c) demethylation trajectory: fraction with promoter mark <=293bp at T1/T2/T3
def near_frac(df):
    return [ (df['dT1']<=W).mean(), (df['dT2']<=W).mean(), (df['dT3']<=W).mean() ]
ax[2].plot([1,2,3], np.array(near_frac(E))*100, '-o', color=C_GCC, label=f"GCCGGC-Exposed (n={len(E)})")
ax[2].plot([1,2,3], np.array(near_frac(AAGp))*100, '-s', color=C_AAG, label=f"AAGCCCG-prom (n={len(AAGp)})")
ax[2].set_xticks([1,2,3]); ax[2].set_xticklabels(['T1\n(12 h)','T2\n(24 h)','T3\n(50 h)'])
ax[2].set_ylabel("% with promoter mark (≤293 bp)"); ax[2].set_ylim(-3,103)
ax[2].set_title("Synchronous demethylation at T2", fontsize=8, pad=6, loc="left")
# both series start high at T1 (left) and collapse by T2; the upper-right is the
# only clear space, so anchor the legend there away from the T1 markers.
ax[2].legend(frameon=False, fontsize=6, loc="upper right", bbox_to_anchor=(1.0, 1.0))
for _i, _lab in enumerate(["a","b","c"]):
    _panel_letter(ax[_i], _lab)
# 2026-09-13 (FIG-07): banner suptitle removed — it asserted the conclusion in the
# figure and collided with the panel letters; the title now lives in the legend.
fig.tight_layout(rect=[0, 0, 1, 0.93], w_pad=2.6)
fig.savefig(FIG/"Figure4_two_system_marking.pdf", bbox_inches="tight"); fig.savefig(FIG/"Figure4_two_system_marking.png", dpi=300, bbox_inches="tight"); plt.close(fig)
_clamp_png_width(FIG/"Figure4_two_system_marking.png")

# ============ FIGURE 5 — synchronized demethylation + weak bias ============
fig, ax = plt.subplots(1, 3, figsize=(FULL_W, 2.7))
# (a) count of Exposed promoters still methylated at each timepoint
cnt = [ (E['dT1']<=W).sum(), (E['dT2']<=W).sum(), (E['dT3']<=W).sum() ]
ax[0].bar(['T1\n(12 h)','T2\n(24 h)','T3\n(50 h)'], cnt, color=[C_BAR, C_ALL, C_ALL])
for i,c in enumerate(cnt): ax[0].text(i, c+0.8, str(int(c)), ha="center", fontsize=10)
ax[0].set_ylabel("Exposed promoters\nmethylated (≤293 bp)")
ax[0].set_title(f"Synchronous erasure ({len(E)}→0 at T2)", fontsize=8, pad=6, loc="left")
# (b) weak bias: promoter GCCGGC occupancy (+-2kb, T1) vs LFC_T2, region-controlled
gT1 = g[g.timepoint=='T1']; P=np.sort(gT1.position.values); Fq=gT1.sort_values('position').frequency.values
def occ(tss,w):
    m=(P>=tss-w)&(P<=tss+w); return Fq[m].sum() if m.any() else 0.0
reg['occ2k']=reg['tss'].apply(lambda t:occ(t,2000))
d=reg[['occ2k','LFC_T2vsT1','region']].dropna()
rx=d.occ2k.rank()-np.polyval(np.polyfit(d.region.rank(),d.occ2k.rank(),1),d.region.rank())
ry=d.LFC_T2vsT1.rank()-np.polyval(np.polyfit(d.region.rank(),d.LFC_T2vsT1.rank(),1),d.region.rank())
rr,pp=pearsonr(rx,ry)
ax[1].scatter(d.occ2k, d.LFC_T2vsT1, s=8, alpha=0.4, c="#6B7783", edgecolors="none")
ax[1].set_xlabel("promoter GCCGGC occupancy (±2 kb, T1)"); ax[1].set_ylabel("log$_2$ FC (T2 vs T1)")
ax[1].set_title(f"Weak modulatory bias\nregion-controlled r = {rr:.2f} (p = {pp:.3f})", fontsize=8, pad=6, loc="left")
ax[1].axhline(0, color="k", lw=0.5)
# (c) bidirectional exemplars: same promoter mark, opposite outcomes (permissive)
def lab(r):
    n = str(r.get('gene_name','')).strip()
    return n if n and n.lower()!='nan' else (str(r.get('old_locus_tag','')).strip() or r['locus_tag'])
ex = E.dropna(subset=['LFC_T2vsT1']).copy(); ex['name'] = ex.apply(lab, axis=1)
up = ex.sort_values('LFC_T2vsT1', ascending=False).head(4)
dn = ex.sort_values('LFC_T2vsT1').head(4)
sel = pd.concat([up, dn]).drop_duplicates('locus_tag').sort_values('LFC_T2vsT1')
colors = [C_REPRESSED if v < 0 else C_INDUCED for v in sel['LFC_T2vsT1']]  # unified: red=repressed, green=induced
ax[2].barh(range(len(sel)), sel['LFC_T2vsT1'].values, color=colors)
ax[2].set_yticks(range(len(sel))); ax[2].set_yticklabels(sel['name'], fontsize=8)
ax[2].axvline(0, color="k", lw=0.6)
ax[2].set_xlabel("log$_2$ FC (T2 vs T1)")
ax[2].set_title("Same mark, opposite outcomes\n(Exposed regulators; permissive)", fontsize=8, pad=6, loc="left")
for _i, _lab in enumerate(["a","b","c"]):
    _panel_letter(ax[_i], _lab)
# 2026-09-13 (FIG-03): banner suptitle removed — interpretive in-figure title
# ("...does NOT direct...") is not NAR style and overlapped the panel letters.
fig.tight_layout(rect=[0, 0, 1, 0.93], w_pad=2.6)
fig.savefig(FIG/"Figure5_synchronized_demethylation.pdf", bbox_inches="tight"); fig.savefig(FIG/"Figure5_synchronized_demethylation.png", dpi=300, bbox_inches="tight"); plt.close(fig)
_clamp_png_width(FIG/"Figure5_synchronized_demethylation.png")

# sync the two PNGs into the manuscript image slots (Figure4.png / Figure5.png)
import shutil, os
_slotdir = Path.home()/"obsidian"/"Research"/"rna-seq"/"Writing"/"fig_images"
# 2026-09-13: sync is opt-in (SYNC_FIG_SLOTS=1). Manuscript slots are now
# Figure5.png (main Fig 4) and SuppFigure20.png (S20); back up before overwriting.
if os.environ.get("SYNC_FIG_SLOTS") == "1" and _slotdir.is_dir():
    shutil.copyfile(FIG/"Figure4_two_system_marking.png", _slotdir/"Figure4.png")
    shutil.copyfile(FIG/"Figure5_synchronized_demethylation.png", _slotdir/"Figure5.png")
    print(f"  Synced → {_slotdir}/Figure4.png, Figure5.png")

print("F4/F5 saved. weak-bias region-controlled r=%.3f p=%.4f" % (rr,pp))
print("MerR/LysR/LacI ORs:", {f:round(o,1) for f,o in zip(fams,ors)})
print("bldD LFC_T2:", reg.loc[reg.locus_tag=='SC_RS09420','LFC_T2vsT1'].values)
