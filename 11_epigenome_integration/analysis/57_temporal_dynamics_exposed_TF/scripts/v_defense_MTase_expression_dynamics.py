#!/usr/bin/env python3
"""
V-Defense MTase expression dynamics + global GCCGGC-m4C density correlation.

Two analyses:
  1. Per-gene T1/T2/T3 expression trajectories for the 12 V-Defense 4mC-hypo genes
     (10/12 are SAM-dependent methyltransferases).
  2. Mean MTase expression vs whole-genome GCCGGC-m4C density across T1/T2/T3.

Inputs:
  - v_defense_hypo_genes.tsv (12 genes)
  - normalized_counts_M145.tsv (DESeq2 size-factor normalized counts, 9 samples)
  - 4mC_final_census.csv (per-site SMRT 4mC calls per timepoint; freq in percent)
  - NC_003888.3.fna (genome reference, used to count total GCCGGC palindromes)

Outputs:
  - analysis/figures/v_defense_MTase_expression_dynamics.{pdf,png}
  - analysis/figures/v_defense_MTase_vs_global_4mC.{pdf,png}
  - analysis/57_*/results/v_defense_MTase_expression_long.tsv
  - analysis/57_*/results/v_defense_global_4mC_density.tsv

Author: Claude Scholar
Date: 2026-05-04
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
BASE = Path("/Users/okaban/bioinfo/rna-seq")
ANALYSIS_DIR = BASE / "11_epigenome_integration/analysis/57_temporal_dynamics_exposed_TF"
GENES_TSV = ANALYSIS_DIR / "results/v_defense_hypo_genes.tsv"
COUNTS_TSV = BASE / "04_deseq2/results/results/normalized_counts_M145.tsv"
CENSUS_CSV = BASE / "11_epigenome_integration/analysis/23_expanded_motif_search/4mC_final_census.csv"
GENOME_FNA = BASE / "11_epigenome_integration/data/NC_003888.3.fna"

FIG_DIR = BASE / "11_epigenome_integration/analysis/figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR = ANALYSIS_DIR / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# DESeq2 sample naming: M145_<timepoint>_<replicate>
TIMEPOINTS = ["T1", "T2", "T3"]
SAMPLE_GROUPS = {
    "T1": ["M145_1_1", "M145_1_2", "M145_1_3"],
    "T2": ["M145_2_1", "M145_2_3", "M145_2_4"],
    "T3": ["M145_3_2", "M145_3_3", "M145_3_4"],
}
TIMEPOINT_COLORS = {"T1": "#4477AA", "T2": "#EE6677", "T3": "#228833"}

# Methylation density: fraction of detected GCCGGC-4mC sites (frequency >=50%) / total
# strand-resolved GCCGGC palindromic positions in the genome.
FREQ_THRESHOLD = 50.0  # percent

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 9,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.1,
})


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def load_genome_total_gccggc(fna: Path) -> int:
    """Return number of strand-resolved GCCGGC sites in the genome.

    GCCGGC is palindromic, so each motif occurrence contributes one
    methylatable C on the forward strand and one on the reverse strand.
    Returns 2 * (number of motif occurrences in the assembled sequence).
    """
    seq = "".join(line.strip() for line in fna.read_text().splitlines() if not line.startswith(">"))
    seq = seq.upper()
    n_fwd = sum(1 for _ in re.finditer("GCCGGC", seq))
    return 2 * n_fwd


def load_expression_long(genes_tsv: Path, counts_tsv: Path) -> pd.DataFrame:
    genes = pd.read_csv(genes_tsv, sep="\t")
    counts = pd.read_csv(counts_tsv, sep="\t").set_index("gene_id")
    keep = counts.index.intersection(genes["gene_id"])
    missing = sorted(set(genes["gene_id"]) - set(keep))
    if missing:
        logger.warning("Missing genes in counts: %s", missing)
    sub = counts.loc[keep]
    rows = []
    for gene_id in genes["gene_id"]:
        if gene_id not in sub.index:
            continue
        meta = genes[genes["gene_id"] == gene_id].iloc[0]
        for tp, samples in SAMPLE_GROUPS.items():
            for sample in samples:
                rows.append({
                    "gene_id": gene_id,
                    "old_locus_tag": meta["old_locus_tag"],
                    "product": meta["product"],
                    "timepoint": tp,
                    "sample": sample,
                    "norm_count": float(sub.at[gene_id, sample]),
                })
    return pd.DataFrame(rows)


def per_timepoint_mean_log2(long_df: pd.DataFrame) -> pd.DataFrame:
    """Return per-(gene, timepoint) mean log2(norm_count + 1)."""
    df = long_df.copy()
    df["log2_count"] = np.log2(df["norm_count"] + 1.0)
    agg = (
        df.groupby(["gene_id", "old_locus_tag", "product", "timepoint"], dropna=False)
          .agg(mean_log2=("log2_count", "mean"),
               sd_log2=("log2_count", "std"),
               n=("log2_count", "size"))
          .reset_index()
    )
    return agg


def compute_global_4mC_density(census: Path, genome_fna: Path) -> pd.DataFrame:
    df = pd.read_csv(census)
    df = df[(df["mod_type"] == "4mC") & (df["final_motif"].str.contains("GCCGGC", na=False))].copy()
    total_sites = load_genome_total_gccggc(genome_fna)
    rows = []
    for tp in TIMEPOINTS:
        sub = df[df["timepoint"] == tp]
        n_called = len(sub)
        n_ge_thr = (sub["frequency"] >= FREQ_THRESHOLD).sum()
        rows.append({
            "timepoint": tp,
            "n_sites_called": n_called,
            "n_sites_ge50pct": int(n_ge_thr),
            "total_GCCGGC_strand_sites": total_sites,
            "density_ge50pct": n_ge_thr / total_sites,
            "mean_freq": sub["frequency"].mean() if n_called else float("nan"),
        })
    return pd.DataFrame(rows)


def compute_mean_mtase_expression(long_df: pd.DataFrame) -> pd.DataFrame:
    """Mean across the 12 genes of per-sample log2(norm_count+1), summarised per timepoint."""
    df = long_df.copy()
    df["log2_count"] = np.log2(df["norm_count"] + 1.0)
    per_sample_mean = df.groupby(["timepoint", "sample"])["log2_count"].mean().reset_index()
    summary = (
        per_sample_mean
        .groupby("timepoint")
        .agg(mean_log2_mtase=("log2_count", "mean"),
             sd_log2_mtase=("log2_count", "std"),
             n_replicates=("log2_count", "size"))
        .reset_index()
    )
    return per_sample_mean, summary


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_expression_dynamics(per_gene: pd.DataFrame, out_stem: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    cmap = plt.colormaps.get_cmap("tab20")
    pivot = per_gene.pivot_table(index="gene_id", columns="timepoint", values="mean_log2")[TIMEPOINTS]

    label_lookup = (
        per_gene.drop_duplicates("gene_id")
                .set_index("gene_id")[["old_locus_tag", "product"]]
    )

    n = len(pivot)
    for i, gid in enumerate(pivot.index):
        y = pivot.loc[gid].values
        sco = label_lookup.at[gid, "old_locus_tag"]
        prod_full = str(label_lookup.at[gid, "product"])
        prod_short = prod_full[:40] + ("..." if len(prod_full) > 40 else "")
        label = f"{gid} ({sco}) — {prod_short}" if isinstance(sco, str) and sco else f"{gid} — {prod_short}"
        ax.plot([0, 1, 2], y, marker="o", lw=1.4, ms=4.5, color=cmap(i / max(1, n - 1)), label=label)

    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(TIMEPOINTS)
    ax.set_xlabel("Developmental timepoint", fontsize=10)
    ax.set_ylabel("log2(normalized count + 1)", fontsize=10)
    ax.set_title(
        "V-Defense 4mC-hypomethylated genes (n=12) — expression trajectory\n"
        "10/12 are SAM-dependent methyltransferases",
        fontsize=10,
    )
    ax.grid(True, alpha=0.25)

    # Count downregulated T1->T2 (mean-of-means)
    t1 = pivot["T1"].values
    t2 = pivot["T2"].values
    n_down = int(np.sum(t2 < t1))
    n_up = int(np.sum(t2 > t1))
    ax.text(0.02, 0.98,
            f"T1->T2 down: {n_down}/{len(pivot)}    T1->T2 up: {n_up}/{len(pivot)}",
            transform=ax.transAxes, va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.6", alpha=0.85),
            fontsize=8)

    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5),
              fontsize=7, frameon=False, title="gene (SCO) — product", title_fontsize=7)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(f"{out_stem}.{ext}")
    plt.close(fig)
    logger.info("Saved %s.{pdf,png}", out_stem)


def plot_mtase_vs_global_4mC(mtase_summary: pd.DataFrame,
                             mtase_per_sample: pd.DataFrame,
                             m4c_density: pd.DataFrame,
                             out_stem: Path) -> tuple[float, float, float, float]:
    fig, ax_left = plt.subplots(figsize=(6.0, 4.2))
    ax_right = ax_left.twinx()

    x = np.arange(len(TIMEPOINTS))
    mtase_summary = mtase_summary.set_index("timepoint").loc[TIMEPOINTS].reset_index()
    m4c_density = m4c_density.set_index("timepoint").loc[TIMEPOINTS].reset_index()

    color_left = "#1f77b4"
    color_right = "#d62728"

    # left: mean MTase expression with replicate dots and SD error bars
    rep = mtase_per_sample.set_index("sample")
    for i, tp in enumerate(TIMEPOINTS):
        vals = rep.loc[rep["timepoint"] == tp, "log2_count"].values
        ax_left.scatter([i] * len(vals), vals, color=color_left, alpha=0.5, s=22, zorder=3)
    ax_left.errorbar(
        x, mtase_summary["mean_log2_mtase"],
        yerr=mtase_summary["sd_log2_mtase"],
        marker="o", ms=7, lw=1.8, color=color_left, capsize=3, zorder=4,
        label="Mean log2 expression of 12 V-Defense MTases (per-sample mean ± SD across n=3 reps)",
    )

    # right: global GCCGGC-m4C density as percentage
    density_pct = m4c_density["density_ge50pct"] * 100.0
    ax_right.plot(
        x, density_pct,
        marker="s", ms=7, lw=1.8, color=color_right,
        label=f"GCCGGC-4mC density (% of {m4c_density['total_GCCGGC_strand_sites'].iloc[0]} strand-sites with freq>=50%)",
    )
    for i, v in enumerate(density_pct):
        ax_right.annotate(f"{v:.2f}%", (x[i], v),
                          textcoords="offset points", xytext=(6, 4),
                          fontsize=8, color=color_right)

    ax_left.set_xticks(x)
    ax_left.set_xticklabels(TIMEPOINTS)
    ax_left.set_xlabel("Developmental timepoint", fontsize=10)
    ax_left.set_ylabel("Mean log2(norm count + 1) of 12 V-Defense MTases",
                      color=color_left, fontsize=10)
    ax_right.set_ylabel("Whole-genome GCCGGC-m4C density (%)", color=color_right, fontsize=10)
    ax_left.tick_params(axis="y", colors=color_left)
    ax_right.tick_params(axis="y", colors=color_right)
    for spine, col in (("left", color_left), ("right", color_right)):
        ax_right.spines[spine].set_color(col)
    ax_left.set_title(
        "V-Defense MTase expression vs whole-genome GCCGGC-m4C density",
        fontsize=10,
    )
    ax_left.grid(True, alpha=0.25)

    # correlations (n=3, descriptive only)
    expr = mtase_summary["mean_log2_mtase"].values
    dens = m4c_density["density_ge50pct"].values
    pearson_r, pearson_p = stats.pearsonr(expr, dens)
    spearman_r, spearman_p = stats.spearmanr(expr, dens)
    ax_left.text(
        0.02, 0.02,
        f"Pearson r = {pearson_r:+.3f} (p={pearson_p:.3f}, n=3)\n"
        f"Spearman rho = {spearman_r:+.3f} (p={spearman_p:.3f}, n=3)\n"
        "[n=3, descriptive only]",
        transform=ax_left.transAxes, va="bottom", ha="left", fontsize=8,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.6", alpha=0.85),
    )

    # combined legend
    h1, l1 = ax_left.get_legend_handles_labels()
    h2, l2 = ax_right.get_legend_handles_labels()
    ax_left.legend(h1 + h2, l1 + l2, loc="upper center",
                   bbox_to_anchor=(0.5, -0.16), ncol=1, frameon=False, fontsize=7.5)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(f"{out_stem}.{ext}")
    plt.close(fig)
    logger.info("Saved %s.{pdf,png}", out_stem)
    return pearson_r, pearson_p, spearman_r, spearman_p


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    logger.info("Loading expression for V-Defense 4mC-hypo genes (n=12)")
    long_df = load_expression_long(GENES_TSV, COUNTS_TSV)
    long_path = RESULTS_DIR / "v_defense_MTase_expression_long.tsv"
    long_df.to_csv(long_path, sep="\t", index=False)
    logger.info("Wrote %s (%d rows, %d unique genes)",
                long_path, len(long_df), long_df["gene_id"].nunique())

    per_gene = per_timepoint_mean_log2(long_df)
    per_gene_path = RESULTS_DIR / "v_defense_MTase_expression_per_gene_timepoint.tsv"
    per_gene.to_csv(per_gene_path, sep="\t", index=False)

    plot_expression_dynamics(per_gene, FIG_DIR / "v_defense_MTase_expression_dynamics")

    logger.info("Computing whole-genome GCCGGC-4mC density per timepoint")
    density = compute_global_4mC_density(CENSUS_CSV, GENOME_FNA)
    density_path = RESULTS_DIR / "v_defense_global_4mC_density.tsv"
    density.to_csv(density_path, sep="\t", index=False)
    logger.info("\n%s", density.to_string(index=False))

    per_sample, summary = compute_mean_mtase_expression(long_df)
    summary_path = RESULTS_DIR / "v_defense_MTase_mean_expression_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    logger.info("\n%s", summary.to_string(index=False))

    pearson_r, pearson_p, spearman_r, spearman_p = plot_mtase_vs_global_4mC(
        summary, per_sample, density, FIG_DIR / "v_defense_MTase_vs_global_4mC",
    )
    correlation_path = RESULTS_DIR / "v_defense_MTase_vs_global_4mC_correlation.tsv"
    pd.DataFrame([{
        "metric": "mean_log2_mtase_expression vs density_ge50pct",
        "n": 3,
        "pearson_r": pearson_r,
        "pearson_p": pearson_p,
        "spearman_r": spearman_r,
        "spearman_p": spearman_p,
        "note": "n=3 timepoints; descriptive only",
    }]).to_csv(correlation_path, sep="\t", index=False)

    # Quick T1->T2 down summary
    pivot = per_gene.pivot_table(index="gene_id", columns="timepoint", values="mean_log2")[TIMEPOINTS]
    n_down_t1t2 = int((pivot["T2"] < pivot["T1"]).sum())
    n_up_t1t2 = int((pivot["T2"] > pivot["T1"]).sum())
    n_down_t1t3 = int((pivot["T3"] < pivot["T1"]).sum())
    logger.info(
        "T1->T2 down: %d/%d  | T1->T2 up: %d/%d  | T1->T3 down: %d/%d",
        n_down_t1t2, len(pivot), n_up_t1t2, len(pivot), n_down_t1t3, len(pivot),
    )


if __name__ == "__main__":
    main()
