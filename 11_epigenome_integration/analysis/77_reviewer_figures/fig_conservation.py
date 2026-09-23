"""Genus-wide conservation / selection figure (candidate MAIN figure).
O/E of the identified motifs across 833 Streptomyces genomes:
AAGCCCG is strongly avoided (purifying selection genus-wide); GCCGGC neutral.
"""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
import importlib.util, sys
B = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis")
OUT = B/"77_reviewer_figures/figures"

# 2026-09-22: go through the approved style module instead of local rcParams and
# hard-coded hexes, so this figure obeys the same palette, font ladder and NAR
# width clamp as the rest of the deck (Writing/check_figures.py enforces it).
_sp = importlib.util.spec_from_file_location(
    "shared_utils", B.parent.parent/"15_paper_figures/scripts/00_shared_utils.py")
su = importlib.util.module_from_spec(_sp); sys.modules["shared_utils"] = su
_sp.loader.exec_module(su)
su.apply_unified_style()

d = pd.read_csv(B/"21_genuswide_motif_conservation/motif_site_density_genuswide.csv")
# motif label, column, our-motif?
# Colour encodes WHICH MOTIF, the same meaning it carries in every other figure
# (author decision 2026-09-22), not "ours versus reference" — the x-axis labels
# already say that, so colour would have been a redundant second encoding.
motifs = [("AAGCCCG\n(this study)","AAGCCCG_oe",su.COL_BOTH),
          ("GCCGGC\n(this study)","GGCCGG_oe",su.COL_4mC),
          ("GATC\n(Dam-type)","GATC_oe",su.COL_GRAY),
          ("CCGG","CCGG_oe",su.COL_GRAY),
          ("CGACNNNCTCC","CGACNNNCTCC_oe",su.COL_GRAY)]
fig, ax = plt.subplots(figsize=(su.mm_to_inch(174), su.mm_to_inch(72)))
data = [d[c].dropna().values for _,c,_ in motifs]
parts = ax.violinplot(data, showmedians=True, widths=0.85)
for i,(lab,c,col) in enumerate(motifs):
    parts['bodies'][i].set_facecolor(col)
    parts['bodies'][i].set_alpha(0.8)
    parts['bodies'][i].set_edgecolor(su.COL_DARK)
    parts['bodies'][i].set_linewidth(0.6)
# violinplot leaves cbars/cmins/cmaxes/cmedians at the default cycle blue, which
# denotes the 6mA mark everywhere else in the deck
for k in ("cbars", "cmins", "cmaxes", "cmedians"):
    if k in parts:
        parts[k].set_color(su.COL_DARK)
        parts[k].set_linewidth(0.9 if k == "cmedians" else 0.6)
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
for i,(lab,c,col) in enumerate(motifs):
    v=d[c].dropna(); med=v.median()
    # nudge the label off the O/E=1 dashed reference line when the median sits on it
    y_lab = med + 0.14 if abs(med - 1.0) < 0.08 else med
    ax.text(i+1+0.46, y_lab, f"med {med:.2f}", ha="left", va="center",
            fontsize=su.FONT["annot"], color=su.COL_DARK)
# Legend: reference line + a colour key for the two motifs from this study.
# Interpretive statements (purifying-selection claim, %-avoided) live in the
# figure legend text, not inside the plot — journal convention.
from matplotlib.patches import Patch
handles = [plt.Line2D([0], [0], color=su.COL_DARK, ls="--", lw=1,
                      label="O/E = 1 (neutral expectation)"),
           Patch(facecolor=su.COL_BOTH, alpha=0.8, label="AAGCCCG"),
           Patch(facecolor=su.COL_4mC, alpha=0.8, label="GCCGGC"),
           Patch(facecolor=su.COL_GRAY, alpha=0.8, label="reference motifs")]
ax.legend(handles=handles, frameon=False, fontsize=su.FONT["legend"],
          loc="upper left", ncol=2, columnspacing=1.2)
# explicit margins: tight_layout is incompatible with the fixed-width save,
# and the default box leaves no room for the two-line x tick labels
fig.subplots_adjust(left=0.085, right=0.985, top=0.88, bottom=0.20)
su.assert_no_text_collisions(fig, "Figure9_conservation")
su.save_figure(fig, OUT/"Figure_genuswide_conservation",
               formats=("pdf", "svg", "png"), width_class="full")
# sync into the Obsidian manuscript slot (Figure 9 previously had no image slot)
import shutil
slot = Path.home()/"obsidian"/"Research"/"rna-seq"/"Writing"/"fig_images"/"Figure9.png"
if slot.parent.is_dir():
    shutil.copyfile(OUT/"Figure_genuswide_conservation.png", slot); print(f"  Synced → {slot}")
plt.close(fig)
print("saved conservation figure")
for lab,c,_ in motifs:
    v=d[c].dropna(); print(f"  {c}: median={v.median():.3f}, %<0.75={100*(v<0.75).mean():.0f}%")
