#!/usr/bin/env python3
"""Split the exposed-TF T1 4mC density × T1→T2 LFC correlation by modification type.

Splits the n=57 coordinated exposed TFs into:
  * 4mC-selected: 4mC counts vary across T1/T2/T3 (drives methyl_change)
  * 6mA-selected: 6mA counts vary across T1/T2/T3
  * Both:         both modification counts vary

For each subgroup, computes Spearman rho between:
  * T1 4mC density (TSS ± 293 bp, sites/100 bp) vs T1→T2 LFC
  * T1 6mA density (TSS ± 293 bp, sites/100 bp) vs T1→T2 LFC

Outputs a 2x2 figure (rows: 4mC-selected vs 6mA-selected;
columns: x = 4mC density vs 6mA density), plus a TSV summary.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

BASE = Path("/Users/okaban/bioinfo/rna-seq")
ANA = BASE / "11_epigenome_integration/analysis"

COORDINATED = ANA / "29_genomewide_TF_screen/tables/coordinated_regulatory_genes.tsv"
EXPOSED_TF = ANA / "58_AAGCCCG_exposed_TF_causal/tables/exposed_TF_methylation_status.tsv"
HC_SITES = ANA / "01_integration/high_confidence_sites_weighted.csv"
FIG_DIR = ANA / "figures"
RESULTS_DIR = ANA / "57_temporal_dynamics_exposed_TF/results"
OBS_FIG = Path("/Users/okaban/obsidian/Research/rna-seq/Writing/Fig_exposed_4mC_6mA_split.png")

WINDOW = 293  # bp

FIG_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def classify_modification(row: pd.Series) -> str:
    """Classify whether a coordinated gene is driven by 4mC, 6mA, or both."""
    counts_4mc = np.array([row["4mC_T1_count"], row["4mC_T2_count"], row["4mC_T3_count"]], dtype=float)
    counts_6ma = np.array([row["6mA_T1_count"], row["6mA_T2_count"], row["6mA_T3_count"]], dtype=float)
    counts_4mc = np.nan_to_num(counts_4mc)
    counts_6ma = np.nan_to_num(counts_6ma)
    var_4mc = counts_4mc.max() != counts_4mc.min()
    var_6ma = counts_6ma.max() != counts_6ma.min()
    if var_4mc and var_6ma:
        return "both"
    if var_4mc:
        return "4mC"
    if var_6ma:
        return "6mA"
    return "neither"


# ------------------------------------------------------------------ #
# 1. Load coordinated TF table and classify by modification type
# ------------------------------------------------------------------ #
coord = pd.read_csv(COORDINATED, sep="\t")
coord["mod_class"] = coord.apply(classify_modification, axis=1)
print("[classification] mod_class counts among coordinated TFs:")
print(coord["mod_class"].value_counts())

# ------------------------------------------------------------------ #
# 2. Load exposed-TF status (TSS, LFC, region) and merge classification
# ------------------------------------------------------------------ #
tf = pd.read_csv(EXPOSED_TF, sep="\t")
df = tf.merge(coord[["locus_tag", "mod_class"]], on="locus_tag", how="inner")
df = df.dropna(subset=["LFC_T2vsT1"]).reset_index(drop=True)
print(f"\n[merge] coordinated × exposed_TF_methylation_status n={len(df)}")

# ------------------------------------------------------------------ #
# 3. Compute T1 4mC and T1 6mA densities at TSS±WINDOW
# ------------------------------------------------------------------ #
sites = pd.read_csv(HC_SITES).rename(columns={"weighted_mod_freq": "freq"})
pos_4mc_t1 = sites.loc[(sites["mod_type"] == "4mC") & (sites["timepoint"] == "T1"), "position"].to_numpy()
pos_6ma_t1 = sites.loc[(sites["mod_type"] == "6mA") & (sites["timepoint"] == "T1"), "position"].to_numpy()


def count_in_window(tss: int, positions: np.ndarray) -> int:
    return int(np.sum(np.abs(positions - tss) <= WINDOW))


df["count_4mC_T1"] = df["tss"].astype(int).apply(lambda t: count_in_window(t, pos_4mc_t1))
df["count_6mA_T1"] = df["tss"].astype(int).apply(lambda t: count_in_window(t, pos_6ma_t1))
df["density_4mC_T1"] = df["count_4mC_T1"] / (2 * WINDOW) * 100
df["density_6mA_T1"] = df["count_6mA_T1"] / (2 * WINDOW) * 100

# ------------------------------------------------------------------ #
# 4. Spearman rho per group
# ------------------------------------------------------------------ #
groups = {
    "4mC-selected": df[df["mod_class"] == "4mC"],
    "6mA-selected": df[df["mod_class"] == "6mA"],
    "Both": df[df["mod_class"] == "both"],
    "All exposed": df,
}

rows = []
for gname, g in groups.items():
    n = len(g)
    for x_label, x_col in (("T1_4mC_density", "density_4mC_T1"), ("T1_6mA_density", "density_6mA_T1")):
        if n >= 3 and g[x_col].nunique() >= 2:
            rho, p = stats.spearmanr(g[x_col], g["LFC_T2vsT1"])
        else:
            rho, p = np.nan, np.nan
        rows.append(
            {
                "group": gname,
                "n": n,
                "x_var": x_label,
                "y_var": "LFC_T2vsT1",
                "spearman_rho": rho,
                "p_value": p,
            }
        )
res = pd.DataFrame(rows)
print("\n[stats] split correlations:")
print(res.to_string(index=False))

out_tsv = RESULTS_DIR / "exposed_4mC_6mA_split_correlation.tsv"
res.to_csv(out_tsv, sep="\t", index=False)
print(f"\n[write] {out_tsv}")

# ------------------------------------------------------------------ #
# 5. Per-gene table for reproducibility
# ------------------------------------------------------------------ #
per_gene_cols = [
    "locus_tag",
    "gene_name",
    "old_locus_tag",
    "tf_family",
    "region",
    "tss",
    "mod_class",
    "count_4mC_T1",
    "density_4mC_T1",
    "count_6mA_T1",
    "density_6mA_T1",
    "LFC_T2vsT1",
    "padj_T2",
]
df[per_gene_cols].to_csv(RESULTS_DIR / "exposed_4mC_6mA_split_per_gene.tsv", sep="\t", index=False)
print(f"[write] {RESULTS_DIR / 'exposed_4mC_6mA_split_per_gene.tsv'}")

# ------------------------------------------------------------------ #
# 6. 2x2 panel figure (rows: 4mC-selected / 6mA-selected;
#                     cols: x=4mC density / 6mA density)
# ------------------------------------------------------------------ #
COLOR_4MC = "#1f6fb4"
COLOR_6MA = "#c0504d"
JITTER = 0.08

panels = [
    ("4mC-selected", "density_4mC_T1", "T1 4mC density (sites / 100 bp, TSS ±293 bp)", COLOR_4MC),
    ("4mC-selected", "density_6mA_T1", "T1 6mA density (sites / 100 bp, TSS ±293 bp)", COLOR_6MA),
    ("6mA-selected", "density_4mC_T1", "T1 4mC density (sites / 100 bp, TSS ±293 bp)", COLOR_4MC),
    ("6mA-selected", "density_6mA_T1", "T1 6mA density (sites / 100 bp, TSS ±293 bp)", COLOR_6MA),
]

fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.6), sharey=True)
rng = np.random.default_rng(0)

for ax, (group_label, x_col, x_title, color) in zip(axes.flatten(), panels):
    g = groups[group_label]
    n = len(g)
    if n == 0:
        ax.set_visible(False)
        continue
    x = g[x_col].to_numpy()
    y = g["LFC_T2vsT1"].to_numpy()
    x_jit = x + rng.uniform(-JITTER, JITTER, size=len(x))
    ax.scatter(x_jit, y, s=42, c=color, alpha=0.85, edgecolor="black", linewidth=0.4)
    ax.axhline(0, color="grey", linewidth=0.6, linestyle="--", alpha=0.7)

    if n >= 3 and g[x_col].nunique() >= 2:
        rho, p = stats.spearmanr(x, y)
        rho_str = f"Spearman ρ = {rho:.3f}\np = {p:.3g}\nn = {n}"
    else:
        rho_str = f"n = {n}\n(insufficient variance)"

    ax.text(
        0.04,
        0.04,
        rho_str,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="grey", alpha=0.85),
    )
    ax.set_title(f"{group_label}  (n={n})", fontsize=11)
    ax.set_xlabel(x_title, fontsize=9)
    if ax in axes[:, 0]:
        ax.set_ylabel("log2FC (T2 vs T1)", fontsize=10)
    ax.tick_params(labelsize=9)
    ax.grid(alpha=0.25, linewidth=0.4)

fig.suptitle(
    "Exposed TF: T1 methylation density vs T1→T2 LFC, split by modification type",
    fontsize=12,
)
fig.tight_layout(rect=(0, 0, 1, 0.96))

pdf_out = FIG_DIR / "Fig_exposed_4mC_vs_6mA_split_correlation.pdf"
png_out = FIG_DIR / "Fig_exposed_4mC_vs_6mA_split_correlation.png"
fig.savefig(pdf_out)
fig.savefig(png_out, dpi=300)
plt.close(fig)
print(f"[fig] {pdf_out}")
print(f"[fig] {png_out}")

OBS_FIG.parent.mkdir(parents=True, exist_ok=True)
import shutil

shutil.copy(png_out, OBS_FIG)
print(f"[copy] {OBS_FIG}")
