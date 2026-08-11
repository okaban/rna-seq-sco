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
plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
                     "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7,
                     "axes.spines.top": False, "axes.spines.right": False})

d = pd.read_csv(B/"21_genuswide_motif_conservation/motif_site_density_genuswide.csv")
# motif label, column, our-motif?
motifs = [("AAGCCCG\n(this study)","AAGCCCG_oe",True),
          ("GCCGGC\n(this study)","GGCCGG_oe",True),
          ("GATC\n(Dam-type)","GATC_oe",False),
          ("CCGG","CCGG_oe",False),
          ("CGACNNNCTCC","CGACNNNCTCC_oe",False)]
fig, ax = plt.subplots(figsize=(6.85, 4.4))   # NAR full width (174 mm)
data = [d[c].dropna().values for _,c,_ in motifs]
parts = ax.violinplot(data, showmedians=True, widths=0.85)
for i,(lab,c,ours) in enumerate(motifs):
    parts['bodies'][i].set_facecolor("#A64B44" if ours else "#9AA7B0")
    parts['bodies'][i].set_alpha(0.75)
ax.axhline(1.0, color="k", ls="--", lw=1, label="O/E = 1 (neutral expectation)")
ax.set_xticks(range(1,len(motifs)+1)); ax.set_xticklabels([m[0] for m in motifs], fontsize=7.5)
ax.set_xlim(0.4, len(motifs)+0.9)   # right headroom for the 'med' labels beside each violin
ax.set_ylabel("Observed/Expected motif frequency\n(per genome, GC-corrected)")
# Y-limit covers the data (incl. the tall CGACNNNCTCC violin) with headroom; the
# title sits above the axes via pad so it never overlaps a violin.
YTOP = 4.2
ax.set_ylim(0, YTOP)
ax.set_title(r"Genus-wide selection on methylation motifs across 833 $\it{Streptomyces}$ genomes",
             pad=12)
# annotate medians to the RIGHT of each violin, at the median height, so the label
# never crosses the central median line drawn by showmedians=True
for i,(lab,c,ours) in enumerate(motifs):
    v=d[c].dropna(); med=v.median()
    # nudge the label off the O/E=1 dashed reference line when the median sits on it
    y_lab = med + 0.14 if abs(med - 1.0) < 0.08 else med
    ax.text(i+1+0.46, y_lab, f"med {med:.2f}", ha="left", va="center", fontsize=7,
            color="#A64B44" if ours else "#555", fontweight="bold" if ours else "normal")
# Legend: reference line + a colour key for the two motifs from this study.
# Interpretive statements (purifying-selection claim, %-avoided) live in the
# figure legend text, not inside the plot — journal convention.
from matplotlib.patches import Patch
handles = [plt.Line2D([0], [0], color="k", ls="--", lw=1, label="O/E = 1 (neutral expectation)"),
           Patch(facecolor="#A64B44", alpha=0.75, label="motifs identified in this study")]
ax.legend(handles=handles, frameon=False, fontsize=7, loc="upper left")
fig.tight_layout()
fig.savefig(OUT/"Figure_genuswide_conservation.pdf"); fig.savefig(OUT/"Figure_genuswide_conservation.png", dpi=300)
# sync into the Obsidian manuscript slot (Figure 9 previously had no image slot)
import shutil
slot = Path.home()/"obsidian"/"Research"/"rna-seq"/"Writing"/"fig_images"/"Figure9.png"
if slot.parent.is_dir():
    shutil.copyfile(OUT/"Figure_genuswide_conservation.png", slot); print(f"  Synced → {slot}")
plt.close(fig)
print("saved conservation figure")
for lab,c,_ in motifs:
    v=d[c].dropna(); print(f"  {c}: median={v.median():.3f}, %<0.75={100*(v<0.75).mean():.0f}%")
