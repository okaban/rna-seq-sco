#!/usr/bin/env python3
"""Shielded TF control: same T1-4mC density vs T1->T2 LFC test.

Hypothesis: the negative correlation between absolute T1 4mC density
(TSS ±293 bp) and T1->T2 log2FC is exposed-TF-specific. If so, the same
test applied to *shielded* TF promoters (all regulatory genes - exposed)
should be null.

Outputs:
  analysis/figures/T1_4mC_vs_T1T2_LFC_exposed_vs_shielded.pdf
  analysis/figures/T1_4mC_vs_T1T2_LFC_exposed_vs_shielded.png
  analysis/57_temporal_dynamics_exposed_TF/results/exposed_vs_shielded_stats.tsv
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

EXPOSED_TF_FILE = ANA / "58_AAGCCCG_exposed_TF_causal/tables/exposed_TF_methylation_status.tsv"
ALL_REG_FILE = ANA / "29_genomewide_TF_screen/tables/all_regulatory_genes.tsv"
HC_SITES = ANA / "01_integration/high_confidence_sites_weighted.csv"
FIG_DIR = ANA / "figures"
RES_DIR = ANA / "57_temporal_dynamics_exposed_TF/results"

WINDOW = 293
COLOR_ZERO = "#9aa3ad"     # gray
COLOR_NONZERO = "#1f6fb4"  # blue


# ------------------------------------------------------------------ #
# Load data
# ------------------------------------------------------------------ #
exposed = pd.read_csv(EXPOSED_TF_FILE, sep="\t")
all_reg = pd.read_csv(ALL_REG_FILE, sep="\t")
sites = pd.read_csv(HC_SITES)

exposed_loci = set(exposed["locus_tag"])

# Compute TSS for shielded TFs
all_reg["tss"] = all_reg.apply(
    lambda r: r["start"] if r["strand"] == "+" else r["end"], axis=1
)

shielded = all_reg[~all_reg["locus_tag"].isin(exposed_loci)].copy()
shielded = shielded.rename(columns={"log2FC_T2vsT1_deseq": "LFC_T2vsT1",
                                    "padj_T2vsT1_deseq": "padj_T2"})
shielded = shielded.dropna(subset=["LFC_T2vsT1", "tss"]).reset_index(drop=True)

# T1 4mC positions
sites_4mC_T1 = sites[(sites["mod_type"] == "4mC") & (sites["timepoint"] == "T1")]
pos = sites_4mC_T1["position"].to_numpy()


def add_count(df, tss_col="tss"):
    counts = []
    for tss in df[tss_col].astype(int):
        counts.append(int(np.sum(np.abs(pos - tss) <= WINDOW)))
    df = df.copy()
    df["count_T1"] = counts
    df["density_per_100bp"] = df["count_T1"] / (2 * WINDOW) * 100
    return df


# Exposed group (use the same construction as v2 main figure)
exposed_full = exposed.dropna(subset=["LFC_T2vsT1"]).copy()
exposed_full = add_count(exposed_full)
shielded = add_count(shielded)


def run_stats(df, label):
    n = len(df)
    rho, p_rho = stats.spearmanr(df["count_T1"], df["LFC_T2vsT1"])
    zero = df["count_T1"] == 0
    nonz = df["count_T1"] >= 1
    n_zero = int(zero.sum())
    n_nonz = int(nonz.sum())
    if n_zero > 0 and n_nonz > 0:
        U, p_mw = stats.mannwhitneyu(
            df.loc[zero, "LFC_T2vsT1"], df.loc[nonz, "LFC_T2vsT1"],
            alternative="two-sided")
        rb_r = 1.0 - 2.0 * U / (n_zero * n_nonz)
        med_zero = float(df.loc[zero, "LFC_T2vsT1"].median())
        med_nonz = float(df.loc[nonz, "LFC_T2vsT1"].median())
    else:
        U = p_mw = rb_r = med_zero = med_nonz = np.nan

    print(f"[stat] {label}: n={n}  zero={n_zero}  nonzero={n_nonz}")
    print(f"        Spearman ρ={rho:+.4f}  p={p_rho:.4g}")
    print(f"        Mann-Whitney U={U}  p={p_mw}  rb-r={rb_r}")
    print(f"        median LFC: zero={med_zero}, nonzero={med_nonz}")
    return dict(group=label, n=n, n_zero=n_zero, n_nonzero=n_nonz,
                spearman_rho=rho, spearman_p=p_rho,
                mw_U=U, mw_p=p_mw, rank_biserial_r=rb_r,
                median_LFC_zero=med_zero, median_LFC_nonzero=med_nonz)


print("=" * 70)
print("Statistics")
print("=" * 70)
res_exposed = run_stats(exposed_full, "exposed")
res_shielded = run_stats(shielded, "shielded")

stats_df = pd.DataFrame([res_exposed, res_shielded])
RES_DIR.mkdir(parents=True, exist_ok=True)
stats_path = RES_DIR / "exposed_vs_shielded_stats.tsv"
stats_df.to_csv(stats_path, sep="\t", index=False)
print(f"[write] {stats_path}")


# ------------------------------------------------------------------ #
# Plot: 2 panels (exposed | shielded)
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

fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.4), sharey=False)

rng = np.random.default_rng(seed=42)


def panel(ax, df, title, jitter_amp_local):
    n = len(df)
    df = df.copy()
    df["density_plot"] = df["density_per_100bp"] + rng.uniform(
        -jitter_amp_local, jitter_amp_local, size=n)

    rho, p_rho = stats.spearmanr(df["count_T1"], df["LFC_T2vsT1"])

    zero = df["count_T1"] == 0
    nonz = df["count_T1"] >= 1
    n_zero = int(zero.sum())
    n_nonz = int(nonz.sum())

    if n_zero > 0 and n_nonz > 0:
        U, p_mw = stats.mannwhitneyu(
            df.loc[zero, "LFC_T2vsT1"], df.loc[nonz, "LFC_T2vsT1"],
            alternative="two-sided")
    else:
        p_mw = np.nan

    # Linear regression for visual fit
    x = df["density_per_100bp"].to_numpy()
    y = df["LFC_T2vsT1"].to_numpy()
    if np.ptp(x) > 0:
        slope, intercept, *_ = stats.linregress(x, y)
        x_range = float(np.ptp(x))
        xx = np.linspace(x.min() - 0.05 * x_range, x.max() + 0.08 * x_range, 200)
        yy = slope * xx + intercept
        # 95 % CI band
        x_mean = x.mean()
        ssx = np.sum((x - x_mean) ** 2)
        resid = y - (slope * x + intercept)
        df_resid = max(n - 2, 1)
        sigma2 = np.sum(resid ** 2) / df_resid
        if ssx > 0:
            se_mean = np.sqrt(sigma2 * (1.0 / n + (xx - x_mean) ** 2 / ssx))
            t_crit = stats.t.ppf(0.975, df_resid)
            ci_lo = yy - t_crit * se_mean
            ci_hi = yy + t_crit * se_mean
            ax.fill_between(xx, ci_lo, ci_hi, color="#444444",
                            alpha=0.10, linewidth=0, zorder=1)
        ax.plot(xx, yy, color="#444444", lw=1.2, zorder=2)

    ax.scatter(df.loc[zero, "density_plot"], df.loc[zero, "LFC_T2vsT1"],
               s=22, facecolor=COLOR_ZERO, edgecolor="black", linewidth=0.3,
               alpha=0.75, zorder=3,
               label=f"count$_{{T1}}$ = 0  (n={n_zero})")
    ax.scatter(df.loc[nonz, "density_plot"], df.loc[nonz, "LFC_T2vsT1"],
               s=26, facecolor=COLOR_NONZERO, edgecolor="black", linewidth=0.3,
               alpha=0.85, zorder=4,
               label=f"count$_{{T1}}$ ≥ 1  (n={n_nonz})")

    ax.axhline(0, color="black", lw=0.5, ls=":", zorder=0)
    ax.set_xlabel("T1 all-4mC site density (per 100 bp, TSS ±293 bp)")
    ax.set_ylabel("log$_2$ fold-change, T1 → T2")
    ax.set_title(title)

    stats_str = (
        f"Spearman ρ = {rho:+.2f}\n"
        f"  p$_{{ρ}}$ = {p_rho:.4f}\n"
        f"Mann–Whitney p = {p_mw:.4f}\n"
        f"n = {n}"
    )
    ax.text(
        0.97, 0.97, stats_str,
        transform=ax.transAxes,
        ha="right", va="top",
        fontsize=8.5,
        bbox=dict(boxstyle="round,pad=0.32", facecolor="white",
                  edgecolor="#888888", linewidth=0.6, alpha=0.95),
    )
    ax.legend(loc="lower right", frameon=False, fontsize=8,
              handlelength=1.2, handletextpad=0.5, borderpad=0.2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


panel(axes[0], exposed_full,
      f"(a) Exposed TFs",
      jitter_amp_local=0.012)
panel(axes[1], shielded,
      f"(b) Shielded TFs (control)",
      jitter_amp_local=0.012)

fig.tight_layout()

PDF = FIG_DIR / "T1_4mC_vs_T1T2_LFC_exposed_vs_shielded.pdf"
PNG = FIG_DIR / "T1_4mC_vs_T1T2_LFC_exposed_vs_shielded.png"
fig.savefig(PDF)
fig.savefig(PNG, dpi=600)
plt.close(fig)
print(f"[write] {PDF}")
print(f"[write] {PNG}  (600 dpi)")
