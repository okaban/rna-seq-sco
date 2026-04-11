#!/usr/bin/env python3
"""
Visualize temporal dynamics of methylation motif sites across three timepoints
in Streptomyces coelicolor M145.

Panels:
  A) 6mA stacked bar chart (site counts by motif, T1/T2/T3)
  B) 4mC stacked bar chart
  C) Fold-change line plot for major motif-site combinations
  D) Heatmap of motif prevalence (%) across timepoints

Input:  temporal_motif_enrichment.csv
Output: motif_temporal_dynamics.{pdf,svg,png}
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.patches import FancyBboxPatch
import matplotlib.colors as mcolors
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration")
INPUT_CSV = BASE / "analysis/18_tss_analyses/temporal_motif_enrichment.csv"
OUT_DIR = BASE / "analysis/02_publication_figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 10,
    "figure.dpi": 300,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "pdf.fonttype": 42,  # TrueType for publication
    "ps.fonttype": 42,
})

# Motif colour palette
MOTIF_COLORS = {
    "CCGG":    "#E53935",
    "AAGCCCG": "#1E88E5",
    "GATC":    "#43A047",
    "GCGC":    "#FB8C00",
    "TCGA":    "#8E24AA",
    "Other":   "#B0BEC5",
}

TIMEPOINTS = ["T1", "T2", "T3"]
# Desired stacking order (bottom to top) for 6mA and 4mC
# Note: motifs overlap (a single site can match multiple motifs), so no "Other"
MOTIF_ORDER_6mA = ["AAGCCCG", "GATC", "GCGC", "CCGG", "TCGA"]
MOTIF_ORDER_4mC = ["CCGG", "AAGCCCG", "GCGC", "GATC", "TCGA"]

# ---------------------------------------------------------------------------
# Load and compute
# ---------------------------------------------------------------------------
df = pd.read_csv(INPUT_CSV)

# Compute absolute site counts per motif per timepoint
for tp in TIMEPOINTS:
    df[f"{tp}_sites"] = (df[f"{tp}_pct"] / 100.0) * df[f"{tp}_n"]

# Build a tidy table of site counts
records = []
for _, row in df.iterrows():
    for tp in TIMEPOINTS:
        records.append({
            "mod_type": row["mod_type"],
            "motif": row["motif"],
            "timepoint": tp,
            "pct": row[f"{tp}_pct"],
            "total_n": row[f"{tp}_n"],
            "site_count": row[f"{tp}_sites"],
        })
tidy = pd.DataFrame(records)

# Note: motifs overlap (a single site can match multiple motifs).
# For 6mA, pct sums to ~104%; for 4mC, pct sums to ~141%.
# Stacked bars therefore represent motif-hit counts, not exclusive partitions.
tidy_full = tidy.copy()

# ---------------------------------------------------------------------------
# Summary statistics
# ---------------------------------------------------------------------------
print("=" * 78)
print("TEMPORAL MOTIF DYNAMICS  --  Streptomyces coelicolor M145")
print("=" * 78)

for mod in ["6mA", "4mC"]:
    print(f"\n{'─' * 40}")
    print(f"  {mod}  --  Absolute site counts")
    print(f"{'─' * 40}")
    sub = tidy_full[tidy_full["mod_type"] == mod].copy()
    pivot = sub.pivot_table(index="motif", columns="timepoint",
                            values="site_count", aggfunc="first")
    pivot = pivot[TIMEPOINTS]
    # reorder rows
    order = MOTIF_ORDER_6mA if mod == "6mA" else MOTIF_ORDER_4mC
    pivot = pivot.reindex(order)
    print(pivot.round(1).to_string())
    for tp in TIMEPOINTS:
        print(f"  Total {tp}: {sub[sub['timepoint']==tp]['total_n'].iloc[0]:.0f}")

print(f"\n{'─' * 40}")
print("  Fold changes from T1 (major motifs)")
print(f"{'─' * 40}")

major_combos = [
    ("6mA", "AAGCCCG"),
    ("6mA", "GATC"),
    ("4mC", "CCGG"),
    ("4mC", "AAGCCCG"),
]

fc_data = {}
for mod, motif in major_combos:
    sub = tidy[(tidy["mod_type"] == mod) & (tidy["motif"] == motif)]
    t1_val = sub.loc[sub["timepoint"] == "T1", "site_count"].values[0]
    vals = {}
    for tp in TIMEPOINTS:
        v = sub.loc[sub["timepoint"] == tp, "site_count"].values[0]
        vals[tp] = v / t1_val if t1_val > 0 else np.nan
    fc_data[(mod, motif)] = vals
    print(f"  {mod}-{motif:>8s}:  "
          f"T1={vals['T1']:.2f}  T2={vals['T2']:.2f}  T3={vals['T3']:.2f}")

print(f"\n{'─' * 40}")
print("  Key observations")
print(f"{'─' * 40}")

# 4mC total sites drop at T3
t1_4mc = df.loc[df["mod_type"] == "4mC", "T1_n"].iloc[0]
t3_4mc = df.loc[df["mod_type"] == "4mC", "T3_n"].iloc[0]
print(f"  * 4mC total sites drop from {t1_4mc:.0f} (T1) to {t3_4mc:.0f} (T3) "
      f"({(t3_4mc/t1_4mc - 1)*100:.1f}%)")

# 4mC-CCGG absolute count
ccgg_t1 = tidy.loc[(tidy["mod_type"]=="4mC") & (tidy["motif"]=="CCGG") & (tidy["timepoint"]=="T1"), "site_count"].values[0]
ccgg_t3 = tidy.loc[(tidy["mod_type"]=="4mC") & (tidy["motif"]=="CCGG") & (tidy["timepoint"]=="T3"), "site_count"].values[0]
print(f"  * 4mC-CCGG sites: {ccgg_t1:.0f} (T1) -> {ccgg_t3:.0f} (T3), "
      f"fold change = {ccgg_t3/ccgg_t1:.2f}")

# 6mA total sites increase
t1_6ma = df.loc[df["mod_type"] == "6mA", "T1_n"].iloc[0]
t3_6ma = df.loc[df["mod_type"] == "6mA", "T3_n"].iloc[0]
print(f"  * 6mA total sites increase from {t1_6ma:.0f} (T1) to {t3_6ma:.0f} (T3) "
      f"(+{(t3_6ma/t1_6ma - 1)*100:.1f}%)")

# 6mA-AAGCCCG
aag_t1 = tidy.loc[(tidy["mod_type"]=="6mA") & (tidy["motif"]=="AAGCCCG") & (tidy["timepoint"]=="T1"), "site_count"].values[0]
aag_t3 = tidy.loc[(tidy["mod_type"]=="6mA") & (tidy["motif"]=="AAGCCCG") & (tidy["timepoint"]=="T3"), "site_count"].values[0]
print(f"  * 6mA-AAGCCCG sites: {aag_t1:.0f} (T1) -> {aag_t3:.0f} (T3), "
      f"fold change = {aag_t3/aag_t1:.2f}")

print()

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(12, 9))
fig.subplots_adjust(hspace=0.38, wspace=0.30, top=0.93, bottom=0.07,
                    left=0.08, right=0.96)

# ── Helper: stacked bar chart ──
def stacked_bar(ax, mod_type, motif_order, panel_label):
    sub = tidy_full[tidy_full["mod_type"] == mod_type].copy()
    x_pos = np.arange(len(TIMEPOINTS))
    bar_width = 0.55
    bottoms = np.zeros(len(TIMEPOINTS))

    for motif in motif_order:
        vals = []
        for tp in TIMEPOINTS:
            row = sub[(sub["motif"] == motif) & (sub["timepoint"] == tp)]
            vals.append(row["site_count"].values[0] if len(row) else 0)
        vals = np.array(vals)
        ax.bar(x_pos, vals, bar_width, bottom=bottoms,
               color=MOTIF_COLORS[motif], edgecolor="white", linewidth=0.5,
               label=motif)
        bottoms += vals

    # Total count labels above bars
    for i, tp in enumerate(TIMEPOINTS):
        total = sub.loc[sub["timepoint"] == tp, "total_n"].iloc[0]
        ax.text(x_pos[i], bottoms[i] + 30, f"n={int(total):,}",
                ha="center", va="bottom", fontsize=8.5, fontweight="bold")

    ax.set_xticks(x_pos)
    ax.set_xticklabels(TIMEPOINTS, fontsize=10)
    ax.set_ylabel("Motif-containing sites", fontsize=10)
    ax.set_title(f"{panel_label}) {mod_type} methylation sites by motif",
                 fontsize=11, fontweight="bold", loc="left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(fontsize=8, frameon=False, ncol=2, loc="upper right",
              bbox_to_anchor=(1.0, 1.0))
    # nice y-limit
    ax.set_ylim(0, bottoms.max() * 1.18)
    # Footnote about overlap
    ax.annotate("*motifs may overlap", xy=(0.98, 0.02),
                xycoords="axes fraction", fontsize=6.5, color="grey",
                ha="right", va="bottom")


# Panel A: 6mA
stacked_bar(axes[0, 0], "6mA", MOTIF_ORDER_6mA, "A")

# Panel B: 4mC
stacked_bar(axes[0, 1], "4mC", MOTIF_ORDER_4mC, "B")

# ── Panel C: Fold-change line plot ──
ax_c = axes[1, 0]
line_combos = [
    ("6mA", "AAGCCCG", "#1E88E5", "o", "-"),
    ("6mA", "GATC",    "#43A047", "s", "-"),
    ("4mC", "CCGG",    "#E53935", "D", "--"),
    ("4mC", "AAGCCCG", "#1E88E5", "^", "--"),
]

x_idx = np.arange(len(TIMEPOINTS))
for mod, motif, color, marker, ls in line_combos:
    sub = tidy[(tidy["mod_type"] == mod) & (tidy["motif"] == motif)]
    t1_val = sub.loc[sub["timepoint"] == "T1", "site_count"].values[0]
    fc_vals = []
    for tp in TIMEPOINTS:
        v = sub.loc[sub["timepoint"] == tp, "site_count"].values[0]
        fc_vals.append(v / t1_val)
    ax_c.plot(x_idx, fc_vals, color=color, marker=marker, markersize=7,
              linewidth=2, linestyle=ls, label=f"{mod}-{motif}", zorder=3)

ax_c.axhline(1.0, color="grey", linewidth=0.8, linestyle=":", zorder=1)
ax_c.set_xticks(x_idx)
ax_c.set_xticklabels(TIMEPOINTS, fontsize=10)
ax_c.set_ylabel("Fold change (relative to T1)", fontsize=10)
ax_c.set_title("C) Temporal fold change of major motif sites",
                fontsize=11, fontweight="bold", loc="left")
ax_c.spines["top"].set_visible(False)
ax_c.spines["right"].set_visible(False)
ax_c.legend(fontsize=8.5, frameon=False, loc="best")
ax_c.set_ylim(0.4, 1.65)
ax_c.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.1f"))

# ── Panel D: Heatmap ──
ax_d = axes[1, 1]

heatmap_combos = [
    ("CCGG",    "4mC"),
    ("AAGCCCG", "6mA"),
    ("AAGCCCG", "4mC"),
    ("GATC",    "6mA"),
    ("GCGC",    "6mA"),
    ("GCGC",    "4mC"),
]

hm_matrix = []
hm_labels = []
for motif, mod in heatmap_combos:
    row_vals = []
    for tp in TIMEPOINTS:
        sub = tidy[(tidy["mod_type"] == mod) & (tidy["motif"] == motif)
                   & (tidy["timepoint"] == tp)]
        row_vals.append(sub["pct"].values[0] if len(sub) else 0)
    hm_matrix.append(row_vals)
    hm_labels.append(f"{motif} ({mod})")

hm_array = np.array(hm_matrix)

# Custom diverging colormap: light -> medium -> dark blue
cmap = plt.cm.YlOrRd

im = ax_d.imshow(hm_array, aspect="auto", cmap=cmap, vmin=0,
                 vmax=hm_array.max() * 1.05)

# Annotate cells
for i in range(hm_array.shape[0]):
    for j in range(hm_array.shape[1]):
        val = hm_array[i, j]
        # Choose text colour for readability
        text_color = "white" if val > 45 else "black"
        ax_d.text(j, i, f"{val:.1f}%", ha="center", va="center",
                  fontsize=9, fontweight="bold", color=text_color)

ax_d.set_xticks(np.arange(len(TIMEPOINTS)))
ax_d.set_xticklabels(TIMEPOINTS, fontsize=10)
ax_d.set_yticks(np.arange(len(hm_labels)))
ax_d.set_yticklabels(hm_labels, fontsize=9)
ax_d.set_title("D) Motif prevalence (%) across timepoints",
                fontsize=11, fontweight="bold", loc="left")

# Colorbar
cbar = fig.colorbar(im, ax=ax_d, shrink=0.75, pad=0.04)
cbar.set_label("Prevalence (%)", fontsize=9)
cbar.ax.tick_params(labelsize=8)

# ── Save ──
for ext in ["pdf", "svg", "png"]:
    out_path = OUT_DIR / f"motif_temporal_dynamics.{ext}"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"Saved: {out_path}")

plt.close(fig)
print("\nDone.")
