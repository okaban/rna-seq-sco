"""Genus-wide conservation / selection figure (candidate MAIN figure).
O/E of the identified motifs across 833 Streptomyces genomes:
AAGCCCG is strongly avoided (purifying selection genus-wide); GCCGGC neutral.
"""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
B = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis")
OUT = B/"77_reviewer_figures/figures"
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})

d = pd.read_csv(B/"21_genuswide_motif_conservation/motif_site_density_genuswide.csv")
# motif label, column, our-motif?
motifs = [("AAGCCCG\n(this study)","AAGCCCG_oe",True),
          ("GCCGGC\n(this study)","GGCCGG_oe",True),
          ("GATC\n(Dam-type)","GATC_oe",False),
          ("CCGG","CCGG_oe",False),
          ("CGACNNNCTCC","CGACNNNCTCC_oe",False)]
fig, ax = plt.subplots(figsize=(8, 5.0))
data = [d[c].dropna().values for _,c,_ in motifs]
parts = ax.violinplot(data, showmedians=True, widths=0.85)
for i,(lab,c,ours) in enumerate(motifs):
    parts['bodies'][i].set_facecolor("#C26B6B" if ours else "#AAAAAA")
    parts['bodies'][i].set_alpha(0.75)
ax.axhline(1.0, color="k", ls="--", lw=1, label="O/E = 1 (neutral expectation)")
ax.set_xticks(range(1,len(motifs)+1)); ax.set_xticklabels([m[0] for m in motifs], fontsize=9)
ax.set_ylabel("Observed/Expected motif frequency\n(per genome, GC-corrected)")
# Y-limit covers the data (incl. the tall CGACNNNCTCC violin) with headroom; the
# title sits above the axes via pad so it never overlaps a violin.
YTOP = 4.2
ax.set_ylim(0, YTOP)
ax.set_title(r"Genus-wide selection on methylation motifs across 833 $\it{Streptomyces}$ genomes",
             pad=12)
# annotate medians just above each violin's bulk (clamped inside the axes so the
# label stays attached to its violin and never floats out at the top)
for i,(lab,c,ours) in enumerate(motifs):
    v=d[c].dropna(); med=v.median()
    ylab = min(np.percentile(v, 95) + 0.12, YTOP - 0.18)
    ax.text(i+1, ylab, f"med {med:.2f}", ha="center", fontsize=8,
            color="#A0484A" if ours else "#555", fontweight="bold" if ours else "normal")
# callout placed in empty space above the low GCCGGC/GATC violins (no leader over data)
ax.text(3.45, 4.08, "AAGCCCG strongly avoided\n(median O/E 0.70; 73% of\ngenomes < 0.75) = long-term\npurifying selection",
        ha="center", va="top", fontsize=8, color="#A0484A",
        bbox=dict(boxstyle="round,pad=0.3", fc="#f7eeee", ec="#C26B6B", lw=0.8))
ax.legend(frameon=False, fontsize=8, loc="upper left")
ax.text(0.5,-0.30,"Claim: AAGCCCG is under genus-wide purifying selection (avoided in 73% of genomes) — a universal, not strain-specific, feature",transform=ax.transAxes,ha="center",fontsize=8,style="italic",color="#444")
ax.text(0.012, 0.80, "red = motifs identified in this study", transform=ax.transAxes,
        fontsize=7.5, color="#A0484A", style="italic")
fig.tight_layout()
fig.savefig(OUT/"Figure_genuswide_conservation.pdf"); fig.savefig(OUT/"Figure_genuswide_conservation.png", dpi=150)
# sync into the Obsidian manuscript slot (Figure 9 previously had no image slot)
import shutil
slot = Path.home()/"obsidian"/"Research"/"rna-seq"/"Writing"/"fig_images"/"Figure9.png"
if slot.parent.is_dir():
    shutil.copyfile(OUT/"Figure_genuswide_conservation.png", slot); print(f"  Synced → {slot}")
plt.close(fig)
print("saved conservation figure")
for lab,c,_ in motifs:
    v=d[c].dropna(); print(f"  {c}: median={v.median():.3f}, %<0.75={100*(v<0.75).mean():.0f}%")
