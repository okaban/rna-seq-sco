#!/usr/bin/env python3
"""Main Figure v2: T1 4mC density at exposed-TF promoters predicts T1->T2 LFC.

Differences vs. v1:
  * Primary inferential test changed to Mann-Whitney U on LFC between
    count_T1=0 (gray) vs count_T1>=1 (blue) groups, since the underlying
    count distribution is effectively binary (44 zero vs 13 nonzero).
  * Spearman correlation kept as a secondary statistic.
  * Category A genes are highlighted with a dark-orange edge ring (no fill change).

Outputs:
  analysis/figures/T1_4mC_vs_T1T2_LFC_main_v2.pdf
  analysis/figures/T1_4mC_vs_T1T2_LFC_main_v2.png
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy import stats

BASE = Path("/Users/okaban/bioinfo/rna-seq")
ANA = BASE / "11_epigenome_integration/analysis"

EXPOSED_TF = ANA / "58_AAGCCCG_exposed_TF_causal/tables/exposed_TF_methylation_status.tsv"
HC_SITES = ANA / "01_integration/high_confidence_sites_weighted.csv"
CATA_PATH = ANA / "57_temporal_dynamics_exposed_TF/data/categoryA_genes.tsv"
FIG_DIR = ANA / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

WINDOW = 293
LABEL_GENES = {"ramR", "SCO4122", "SCO7314", "tcrA"}
COLOR_ZERO = "#9aa3ad"   # gray
COLOR_NONZERO = "#1f6fb4"  # blue
COLOR_HIGHLIGHT = "darkorange"

# ------------------------------------------------------------------ #
# Load & build per-TF table
# ------------------------------------------------------------------ #
tf = pd.read_csv(EXPOSED_TF, sep="\t")
sites = pd.read_csv(HC_SITES).rename(columns={"weighted_mod_freq": "freq"})
catA = pd.read_csv(CATA_PATH, sep="\t")
catA_set = set(catA["locus_tag"])

sites_4mC_T1 = sites[(sites["mod_type"] == "4mC") & (sites["timepoint"] == "T1")]
pos = sites_4mC_T1["position"].to_numpy()

records = []
for _, r in tf.iterrows():
    tss = int(r["tss"])
    cnt = int(np.sum(np.abs(pos - tss) <= WINDOW))
    records.append({"locus_tag": r["locus_tag"], "count_T1": cnt})
cnts = pd.DataFrame(records)

df = tf.merge(cnts, on="locus_tag", how="left")
df = df.dropna(subset=["LFC_T2vsT1", "count_T1"]).reset_index(drop=True)
df["density_per_100bp"] = df["count_T1"] / (2 * WINDOW) * 100

n = len(df)

# ------------------------------------------------------------------ #
# Statistics
# ------------------------------------------------------------------ #
rho, p_rho = stats.spearmanr(df["count_T1"], df["LFC_T2vsT1"])

zero_mask = df["count_T1"] == 0
nonzero_mask = df["count_T1"] >= 1
n_zero = int(zero_mask.sum())
n_nonzero = int(nonzero_mask.sum())

lfc_zero = df.loc[zero_mask, "LFC_T2vsT1"].to_numpy()
lfc_nonzero = df.loc[nonzero_mask, "LFC_T2vsT1"].to_numpy()

U, p_mw = stats.mannwhitneyu(lfc_zero, lfc_nonzero, alternative="two-sided")
# rank-biserial effect size: r = 1 - 2U / (n1 * n2)
rb_r = 1.0 - 2.0 * U / (n_zero * n_nonzero)

print(f"[stat] n={n}  count_T1==0: {n_zero},  count_T1>=1: {n_nonzero}")
print(f"[stat] Spearman rho={rho:+.4f}  p={p_rho:.4g}")
print(f"[stat] Mann-Whitney U={U:.0f}  p={p_mw:.4g}  rank-biserial r={rb_r:+.3f}")
print(f"[stat]   median LFC (count=0):  {np.median(lfc_zero):+.3f}")
print(f"[stat]   median LFC (count>=1): {np.median(lfc_nonzero):+.3f}")

# Linear regression (kept for the visual fit & 95% band)
x = df["density_per_100bp"].to_numpy()
y = df["LFC_T2vsT1"].to_numpy()
slope, intercept, r_pearson, p_pearson, se_slope = stats.linregress(x, y)
x_range = float(np.ptp(x))
xx = np.linspace(x.min() - 0.05 * x_range, x.max() + 0.08 * x_range, 200)
yy = slope * xx + intercept

# Horizontal jitter for plotting only (counts are 0/1/2, points stack)
rng = np.random.default_rng(seed=42)
jitter_amp = 0.012  # density units (per 100 bp)
df["density_plot"] = df["density_per_100bp"] + rng.uniform(
    -jitter_amp, jitter_amp, size=len(df))

x_mean = x.mean()
ssx = np.sum((x - x_mean) ** 2)
resid = y - (slope * x + intercept)
df_resid = n - 2
sigma2 = np.sum(resid ** 2) / df_resid
se_mean = np.sqrt(sigma2 * (1.0 / n + (xx - x_mean) ** 2 / ssx))
t_crit = stats.t.ppf(0.975, df_resid)
ci_lo = yy - t_crit * se_mean
ci_hi = yy + t_crit * se_mean

# ------------------------------------------------------------------ #
# Plot
# ------------------------------------------------------------------ #
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 10,
    "axes.linewidth": 0.8,
    "axes.labelsize": 11,
    "axes.titlesize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 8.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

fig, ax = plt.subplots(figsize=(4.8, 4.4))

# Regression band + line
ax.fill_between(xx, ci_lo, ci_hi, color="#444444", alpha=0.10, linewidth=0,
                zorder=1, label="95% CI")
ax.plot(xx, yy, color="#444444", lw=1.2, zorder=2, label="linear fit")

# --------------------------- scatter ----------------------------- #
is_A = df["locus_tag"].isin(catA_set)


def _scatter(group_df, *, facecolor, label):
    """Plot a group with Category A genes ringed in dark-orange."""
    if len(group_df) == 0:
        return
    nonA = group_df.loc[~group_df["locus_tag"].isin(catA_set)]
    A = group_df.loc[group_df["locus_tag"].isin(catA_set)]

    if len(nonA) > 0:
        ax.scatter(nonA["density_plot"], nonA["LFC_T2vsT1"],
                   s=30, facecolor=facecolor, edgecolor="black",
                   linewidth=0.4, alpha=0.9, zorder=3, label=label)
    if len(A) > 0:
        ax.scatter(A["density_plot"], A["LFC_T2vsT1"],
                   s=46, facecolor=facecolor, edgecolor=COLOR_HIGHLIGHT,
                   linewidth=1.5, alpha=0.95, zorder=4)


_scatter(df.loc[zero_mask], facecolor=COLOR_ZERO,
         label=f"count$_{{T1}}$ = 0  (n={n_zero})")
_scatter(df.loc[nonzero_mask], facecolor=COLOR_NONZERO,
         label=f"count$_{{T1}}$ ≥ 1  (n={n_nonzero})")

# Reference axes
ax.axhline(0, color="black", lw=0.5, ls=":", zorder=0)


# Labels for the four highlighted genes
def resolve_name(row) -> str:
    name = str(row.get("gene_name", "")).strip()
    if name and name.lower() != "nan":
        return name
    old = str(row.get("old_locus_tag", "")).strip()
    if old and old.lower() != "nan":
        return old
    return str(row["locus_tag"])


for _, row in df.iterrows():
    label = resolve_name(row)
    if label not in LABEL_GENES:
        continue
    xv = row["density_plot"]
    yv = row["LFC_T2vsT1"]
    offsets = {
        "ramR":    (22, 0),
        "SCO7314": (22, 0),
        "tcrA":    (22, 0),
        "SCO4122": (22, 0),
    }
    dx, dy = offsets.get(label, (16, 0))
    ax.annotate(
        label,
        xy=(xv, yv),
        xytext=(dx, dy),
        textcoords="offset points",
        fontsize=9,
        fontstyle="italic",
        ha="left",
        va="center",
        arrowprops=dict(arrowstyle="-", lw=0.5, color="black",
                        shrinkA=0, shrinkB=2),
    )

# Stats annotation (top-right)
stats_str = (
    f"Spearman ρ = {rho:+.2f}  (p = {p_rho:.4f})\n"
    f"Mann–Whitney U  p = {p_mw:.4f}\n"
    f"  (count$_{{T1}}$=0 vs ≥1; r = {rb_r:+.2f})\n"
    f"n = {n}"
)
ax.text(
    0.97, 0.97, stats_str,
    transform=ax.transAxes,
    ha="right", va="top",
    fontsize=9,
    bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
              edgecolor="#888888", linewidth=0.6, alpha=0.95),
)

ax.set_xlabel("T1 all-4mC site density (sites per 100 bp, TSS ±293 bp)")
ax.set_ylabel("log$_2$ fold-change, T1 → T2")
ax.set_xlim(x.min() - 0.06 * x_range, x.max() + 0.06 * x_range)

# Legend (manual handles to keep the order tidy and include Cat-A explanation)
n_A = int(is_A.sum())
handles = [
    Line2D([0], [0], color="#444444", lw=1.2, label="linear fit"),
    plt.Rectangle((0, 0), 1, 1, fc="#444444", alpha=0.10, ec="none", label="95% CI"),
    Line2D([0], [0], marker="o", color="none",
           markerfacecolor=COLOR_ZERO, markeredgecolor="black", markersize=6,
           markeredgewidth=0.4, label=f"count$_{{T1}}$ = 0  (n={n_zero})"),
    Line2D([0], [0], marker="o", color="none",
           markerfacecolor=COLOR_NONZERO, markeredgecolor="black", markersize=6,
           markeredgewidth=0.4, label=f"count$_{{T1}}$ ≥ 1  (n={n_nonzero})"),
    Line2D([0], [0], marker="o", color="none",
           markerfacecolor="white", markeredgecolor=COLOR_HIGHLIGHT, markersize=7,
           markeredgewidth=1.5, label=f"Category A  (n={n_A})"),
]
ax.legend(handles=handles, loc="lower right", frameon=False,
          handlelength=1.6, handletextpad=0.6, borderpad=0.2)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.tight_layout()

# ------------------------------------------------------------------ #
# Save
# ------------------------------------------------------------------ #
PDF = FIG_DIR / "T1_4mC_vs_T1T2_LFC_main_v2.pdf"
PNG = FIG_DIR / "T1_4mC_vs_T1T2_LFC_main_v2.png"
fig.savefig(PDF)
fig.savefig(PNG, dpi=600)
plt.close(fig)

print(f"[write] {PDF}")
print(f"[write] {PNG}  (600 dpi)")
