#!/usr/bin/env python3
"""
Figure 1 panel: GCCGGC methylation site distribution across genome regions
(left arm / core / right arm) at T1, T2, T3.

Genome boundaries (S. coelicolor M145, NC_003888.3, ~8.7 Mb):
  Left arm : 0 – 1.5 Mb
  Core     : 1.5 – 6.5 Mb
  Right arm: 6.5 – 8.7 Mb

Input : 11_epigenome_integration/analysis/37_defense_island_GCCGGC/
        tables/GCCGGC_sites_by_timepoint.tsv
Output: 15_paper_figures/figures/main/Figure1_GCCGGC_arm_distribution.png
        Writing/fig_images/Figure1_arm_panel.png
"""

import os
import shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

# ── Paths ────────────────────────────────────────────────────────────────────
BASE = Path(__file__).resolve().parents[2]          # rna-seq/
DATA = (BASE / "11_epigenome_integration" / "analysis"
        / "37_defense_island_GCCGGC" / "tables"
        / "GCCGGC_sites_by_timepoint.tsv")
OUT_MAIN = BASE / "15_paper_figures" / "figures" / "main" / "Figure1_GCCGGC_arm_distribution.png"
OUT_WRITE = BASE / "Writing" / "fig_images" / "Figure1_arm_panel.png"

# ── Genome boundaries (bp) ────────────────────────────────────────────────────
LEFT_ARM_END  = 1_500_000
RIGHT_ARM_START = 6_500_000

# ── Colours ──────────────────────────────────────────────────────────────────
COL_LEFT  = "#FF5722"   # left arm  – orange-red
COL_CORE  = "#2196F3"   # core      – blue
COL_RIGHT = "#4CAF50"   # right arm – green

# ── Load & classify ──────────────────────────────────────────────────────────
df = pd.read_csv(DATA, sep="\t")

def classify(pos):
    if pos < LEFT_ARM_END:
        return "left_arm"
    elif pos <= RIGHT_ARM_START:
        return "core"
    else:
        return "right_arm"

df["arm_region"] = df["position"].apply(classify)

# ── Count sites per timepoint × region ───────────────────────────────────────
timepoints = ["T1", "T2", "T3"]
regions    = ["left_arm", "core", "right_arm"]

counts = (df.groupby(["timepoint", "arm_region"])
            .size()
            .unstack(fill_value=0)
            .reindex(index=timepoints, columns=regions, fill_value=0))

totals = counts.sum(axis=1)
fracs  = counts.div(totals, axis=0) * 100   # percent

# ── Figure layout ─────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(9, 4.2),
                          gridspec_kw={"wspace": 0.38})

region_colors  = [COL_LEFT, COL_CORE, COL_RIGHT]
region_labels  = ["Left arm\n(<1.5 Mb)", "Core\n(1.5–6.5 Mb)", "Right arm\n(>6.5 Mb)"]
x = np.arange(len(timepoints))
bar_w = 0.55

for ax, data, ylabel, title_tag in [
    (axes[0], counts, "Number of methylated sites", "(i)"),
    (axes[1], fracs,  "Fraction of methylated sites (%)", "(ii)"),
]:
    bottom = np.zeros(len(timepoints))
    for col, color in zip(regions, region_colors):
        vals = data[col].values.astype(float)
        ax.bar(x, vals, bar_w, bottom=bottom, color=color,
               edgecolor="white", linewidth=0.6, zorder=3)
        # Annotate bars for counts panel (skip very small slices)
        if data is counts:
            for xi, (v, b) in enumerate(zip(vals, bottom)):
                if v >= 15:
                    ax.text(xi, b + v / 2, f"{int(v)}",
                            ha="center", va="center",
                            fontsize=7.5, color="white", fontweight="bold")
        bottom = bottom + vals

    ax.set_xticks(x)
    ax.set_xticklabels(timepoints, fontsize=11)
    ax.set_xlabel("Timepoint", fontsize=11, labelpad=6)
    ax.set_ylabel(ylabel, fontsize=10, labelpad=6)
    ax.tick_params(axis="y", labelsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.yaxis.grid(True, linestyle="--", linewidth=0.5, alpha=0.6, zorder=0)
    ax.set_axisbelow(True)

    # Total-count annotation above each bar
    if data is counts:
        for xi, tot in enumerate(totals):
            ax.text(xi, tot + totals.max() * 0.02, f"n={int(tot)}",
                    ha="center", va="bottom", fontsize=8, color="#333333")
    else:
        ax.set_ylim(0, 115)
        ax.yaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:.0f}%"))

# ── Legend ────────────────────────────────────────────────────────────────────
patches = [mpatches.Patch(color=c, label=l)
           for c, l in zip(region_colors, region_labels)]
axes[1].legend(handles=patches, loc="upper right", fontsize=8.5,
               frameon=True, framealpha=0.85, edgecolor="#cccccc",
               title="Genomic region", title_fontsize=8.5)

# ── Panel label ───────────────────────────────────────────────────────────────
fig.text(0.01, 0.97, "A", fontsize=16, fontweight="bold",
         va="top", ha="left")

fig.suptitle("GCCGGC methylation site distribution across genome regions",
             fontsize=11, y=1.01, ha="center", color="#222222")

# ── Save ──────────────────────────────────────────────────────────────────────
for out_path in [OUT_MAIN, OUT_WRITE]:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight",
                facecolor="white", transparent=False)
    print(f"Saved: {out_path}")

plt.close(fig)
print("Done.")
