#!/usr/bin/env python3
"""Main Figure: T1 4mC site density at exposed-TF promoters predicts T1->T2 repression.

Recomputes the headline statistic ρ=-0.40, p=0.0019 (n=57 exposed TFs;
all_4mC sites within TSS ±293 bp; count_T1 vs LFC_T2vsT1) and renders a
publication-quality scatter with regression line, 95 % CI band, and Category A
gene highlights.

Outputs:
  analysis/figures/T1_4mC_vs_T1T2_LFC_main.pdf
  analysis/figures/T1_4mC_vs_T1T2_LFC_main.png
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

# ------------------------------------------------------------------ #
# Load
# ------------------------------------------------------------------ #
tf = pd.read_csv(EXPOSED_TF, sep="\t")
sites = pd.read_csv(HC_SITES).rename(columns={"weighted_mod_freq": "freq"})
catA = pd.read_csv(CATA_PATH, sep="\t")
catA_set = set(catA["locus_tag"])

# Per-TF count of all 4mC sites at T1 within TSS ±293 bp
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

# Density (sites per 100 bp): just a linear rescale; Spearman ρ unchanged
df["density_per_100bp"] = df["count_T1"] / (2 * WINDOW) * 100

n = len(df)
rho, p = stats.spearmanr(df["count_T1"], df["LFC_T2vsT1"])
print(f"[stat] n={n}  Spearman rho={rho:.4f}  p={p:.4g}")

# Linear regression on (density, LFC) for the visual fit + 95 % CI
x = df["density_per_100bp"].to_numpy()
y = df["LFC_T2vsT1"].to_numpy()

slope, intercept, r_pearson, p_pearson, se_slope = stats.linregress(x, y)
x_range = float(np.ptp(x))
xx = np.linspace(x.min() - 0.05 * x_range, x.max() + 0.08 * x_range, 200)
yy = slope * xx + intercept

# Horizontal jitter for plotting only (counts are 0/1/2, so points stack)
rng = np.random.default_rng(seed=42)
jitter_amp = 0.012  # in density units (per 100 bp)
df["density_plot"] = df["density_per_100bp"] + rng.uniform(
    -jitter_amp, jitter_amp, size=len(df))

# 95 % CI band for the mean response (E[y|x])
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
    "legend.fontsize": 9,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

fig, ax = plt.subplots(figsize=(4.6, 4.2))

# Regression band + line
ax.fill_between(xx, ci_lo, ci_hi, color="#444444", alpha=0.12, linewidth=0,
                zorder=1, label="95% CI")
ax.plot(xx, yy, color="#444444", lw=1.2, zorder=2,
        label=f"linear fit (β={slope:+.2f})")

# Scatter: non-Category A
is_A = df["locus_tag"].isin(catA_set)
non_A = df.loc[~is_A]
A = df.loc[is_A]

ax.scatter(non_A["density_plot"], non_A["LFC_T2vsT1"],
           s=28, facecolor="#9aa3ad", edgecolor="black", linewidth=0.4,
           alpha=0.85, zorder=3, label=f"other exposed TF (n={len(non_A)})")

ax.scatter(A["density_plot"], A["LFC_T2vsT1"],
           s=42, facecolor="#e15a1f", edgecolor="black", linewidth=0.5,
           alpha=0.95, zorder=4, label=f"Category A (n={len(A)})")

# Reference axes
ax.axhline(0, color="black", lw=0.5, ls=":", zorder=0)

# Gene labels: name resolution
def resolve_name(row) -> str:
    name = str(row.get("gene_name", "")).strip()
    if name and name.lower() != "nan":
        return name
    old = str(row.get("old_locus_tag", "")).strip()
    if old and old.lower() != "nan":
        return old
    return str(row["locus_tag"])

# Label requested genes (and only those — keep figure clean)
for _, row in df.iterrows():
    label = resolve_name(row)
    if label not in LABEL_GENES:
        continue
    xv = row["density_plot"]
    yv = row["LFC_T2vsT1"]
    # offset directions chosen per-gene to avoid overlap; all four sit at count=0
    offsets = {
        "ramR":    ( 22,   0),   # LFC=+8.85, far above
        "SCO7314": ( 22,   0),   # LFC=+1.91
        "tcrA":    ( 22,   0),   # LFC=+1.19
        "SCO4122": ( 22,   0),   # LFC=-3.86, far below
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

# Stats annotation (top-right corner)
stats_str = (
    f"Spearman ρ = {rho:+.2f}\n"
    f"p = {p:.4f}\n"
    f"n = {n}"
)
ax.text(
    0.97, 0.97, stats_str,
    transform=ax.transAxes,
    ha="right", va="top",
    fontsize=10,
    bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
              edgecolor="#888888", linewidth=0.6, alpha=0.95),
)

# Axes labels & title
ax.set_xlabel("T1 4mC site density (sites per 100 bp, TSS ±293 bp)")
ax.set_ylabel("log$_2$ fold-change, T1 → T2")

# Tight x-limits with a small margin
ax.set_xlim(x.min() - 0.06 * x_range, x.max() + 0.06 * x_range)

# Build a compact legend (drop autogenerated "linear fit (β=...)" handle for clarity)
handles = [
    Line2D([0], [0], color="#444444", lw=1.2, label="linear fit"),
    plt.Rectangle((0, 0), 1, 1, fc="#444444", alpha=0.12, ec="none", label="95% CI"),
    Line2D([0], [0], marker="o", color="none",
           markerfacecolor="#9aa3ad", markeredgecolor="black", markersize=6,
           markeredgewidth=0.4, label=f"other exposed TF (n={len(non_A)})"),
    Line2D([0], [0], marker="o", color="none",
           markerfacecolor="#e15a1f", markeredgecolor="black", markersize=7,
           markeredgewidth=0.5, label=f"Category A (n={len(A)})"),
]
ax.legend(handles=handles, loc="lower right", frameon=False,
          handlelength=1.6, handletextpad=0.6, borderpad=0.2)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.tight_layout()

# ------------------------------------------------------------------ #
# Save
# ------------------------------------------------------------------ #
PDF = FIG_DIR / "T1_4mC_vs_T1T2_LFC_main.pdf"
PNG = FIG_DIR / "T1_4mC_vs_T1T2_LFC_main.png"
fig.savefig(PDF)
fig.savefig(PNG, dpi=600)
plt.close(fig)

print(f"[write] {PDF}")
print(f"[write] {PNG}  (600 dpi)")

# Persist underlying numbers for record
TSV = ANA / "57_temporal_dynamics_exposed_TF/results/main_fig_T1_4mC_vs_T1T2_LFC_data.tsv"
out_cols = [
    "locus_tag", "gene_name", "old_locus_tag", "tf_family", "region",
    "tss", "count_T1", "density_per_100bp", "LFC_T2vsT1", "padj_T2",
]
out_cols = [c for c in out_cols if c in df.columns]
df["is_categoryA"] = df["locus_tag"].isin(catA_set)
df[out_cols + ["is_categoryA"]].sort_values("density_per_100bp", ascending=False)\
    .to_csv(TSV, sep="\t", index=False)
print(f"[write] {TSV}")
