#!/usr/bin/env python3
"""
Motif decomposition + window sensitivity for the
Exposed-TF "TSS-proximal methylation change vs T1->T2 LFC" scatter.

Tasks
-----
1. Decompose by motif:
   - GCCGGC (4mC, R-M defense)
   - AAGCCCG (6mA, candidate regulatory)
   Per-motif Spearman on Exposed (n=57) and Shielded (n=355).

2. Window sensitivity:
   - Re-run combined (both motifs) at TSS +/- {200, 293, 500, 1000} bp.
   - Print a comparison table (rho_exposed, p, rho_shielded, p) for each window.

Outputs
-------
- exposed_TF_motif_decomposition.{pdf,svg,png}: 2-panel scatter
- exposed_TF_window_sensitivity.tsv: comparison table
- exposed_TF_window_sensitivity.{pdf,svg,png}: rho vs window plot
"""

from pathlib import Path
import logging

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

BASE = Path("/Users/okaban/bioinfo/rna-seq")
EPI = BASE / "11_epigenome_integration/analysis"

EXPOSED_FILE = EPI / "51_exposed_regulators_characteristics/tables/exposed_regulators_full_table.tsv"
ALL_FEATURES = EPI / "52_shielded_exposed_boundary/tables/all_genes_features.tsv"
METHYL_FILE = EPI / "01_integration/high_confidence_sites_weighted.csv"
DESEQ_T2 = BASE / "04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_2_vs_1.tsv"
AAGCCCG_FILE = EPI / "36_AAGCCCG_distribution/tables/AAGCCCG_site_gene_mapping.tsv"
GCCGGC_FILE = EPI / "37_defense_island_GCCGGC/tables/GCCGGC_sites_by_timepoint.tsv"

OUT_DIR_LOCAL = EPI / "57_temporal_dynamics_exposed_TF/figures"
OUT_DIR_TOP = EPI / "figures"
TABLE_DIR = EPI / "57_temporal_dynamics_exposed_TF/tables"


def load_universe() -> pd.DataFrame:
    """Build TF universe (57 Exposed + 355 Shielded) with TSS."""
    exposed_df = pd.read_csv(EXPOSED_FILE, sep="\t")
    all_features = pd.read_csv(ALL_FEATURES, sep="\t")
    exposed_set = set(exposed_df["locus_tag"].dropna())
    h29_tss_by_locus = dict(zip(all_features["locus_tag"], all_features["tss"]))

    rows = []
    for _, r in exposed_df.iterrows():
        lt = r["locus_tag"]
        if lt in h29_tss_by_locus and not pd.isna(h29_tss_by_locus[lt]):
            tss = int(h29_tss_by_locus[lt])
            src = "experimental"
        else:
            tss = int(r["start"]) if r["strand"] == "+" else int(r["end"])
            src = "gene_start"
        rows.append({"locus_tag": lt, "tss": tss, "tss_source": src, "is_exposed": 1})
    for _, r in all_features.iterrows():
        if r["locus_tag"] in exposed_set or pd.isna(r["tss"]):
            continue
        rows.append({"locus_tag": r["locus_tag"], "tss": int(r["tss"]),
                     "tss_source": "experimental", "is_exposed": 0})
    return pd.DataFrame(rows)


def load_methyl_with_motif() -> pd.DataFrame:
    """Load high-confidence sites + tag motif annotation."""
    hc = pd.read_csv(METHYL_FILE)
    hc = hc.loc[hc["timepoint"].isin(["T1", "T2"])].copy()

    aagc = pd.read_csv(AAGCCCG_FILE, sep="\t")
    gccg = pd.read_csv(GCCGGC_FILE, sep="\t")
    aagcccg_pos = set(aagc["position"].astype(int))
    gccggc_pos = set(gccg["position"].astype(int))

    def motif_label(row):
        p = int(row["position"])
        m = row["mod_type"]
        if m == "6mA" and p in aagcccg_pos:
            return "AAGCCCG"
        if m == "4mC" and p in gccggc_pos:
            return "GCCGGC"
        return "other"

    hc["motif"] = hc.apply(motif_label, axis=1)
    logger.info(f"  Methyl-site motif tagging: {hc['motif'].value_counts().to_dict()}")
    return hc


def methyl_score(methyl_df: pd.DataFrame, tss: int, window: int) -> tuple[float, float]:
    """Sum of weighted_mod_freq (fraction) for sites in TSS +/- window, T1 and T2."""
    sub = methyl_df.loc[
        (methyl_df["position"] >= tss - window) & (methyl_df["position"] <= tss + window)
    ]
    t1 = sub.loc[sub["timepoint"] == "T1", "weighted_mod_freq"].sum() / 100.0
    t2 = sub.loc[sub["timepoint"] == "T2", "weighted_mod_freq"].sum() / 100.0
    return t1, t2


def build_scores(universe: pd.DataFrame, methyl: pd.DataFrame, window: int) -> pd.DataFrame:
    """Per-gene methylation change at given window; returns df with methyl_change."""
    rows = []
    for _, r in universe.iterrows():
        t1, t2 = methyl_score(methyl, int(r["tss"]), window)
        rows.append({**r.to_dict(), "methyl_change": t2 - t1, "T1": t1, "T2": t2})
    return pd.DataFrame(rows)


def add_lfc(df: pd.DataFrame, deseq: pd.DataFrame) -> pd.DataFrame:
    return df.merge(
        deseq.rename(columns={"gene_id": "locus_tag", "log2FoldChange": "LFC_T2vsT1"})[
            ["locus_tag", "LFC_T2vsT1", "padj"]
        ],
        on="locus_tag",
        how="left",
    ).dropna(subset=["LFC_T2vsT1"])


def spearman(df: pd.DataFrame) -> tuple[float, float, int]:
    if len(df) < 5:
        return np.nan, np.nan, len(df)
    rho, p = stats.spearmanr(df["methyl_change"], df["LFC_T2vsT1"])
    return rho, p, len(df)


# ============================================================================
# Task 1: Motif decomposition
# ============================================================================
def task1_motif_decomposition(universe: pd.DataFrame, methyl_full: pd.DataFrame,
                               deseq: pd.DataFrame, window: int = 500) -> dict:
    logger.info("\n=== Task 1: Motif decomposition (window = +/- %d bp) ===" % window)

    fig, axes = plt.subplots(1, 2, figsize=(13, 6), sharey=True)

    summary = {}
    for ax, motif, color in zip(axes, ["GCCGGC", "AAGCCCG"], ["#1F77B4", "#2CA02C"]):
        m = methyl_full[methyl_full["motif"] == motif]
        df = build_scores(universe, m, window)
        df = add_lfc(df, deseq)
        exp = df[df.is_exposed == 1]
        shld = df[df.is_exposed == 0]

        rho_e, p_e, n_e = spearman(exp)
        rho_s, p_s, n_s = spearman(shld)
        logger.info(f"  [{motif}] Exposed:  rho={rho_e:.4f}, p={p_e:.4g}, n={n_e}")
        logger.info(f"  [{motif}] Shielded: rho={rho_s:.4f}, p={p_s:.4g}, n={n_s}")

        # Quadrant counts (Exposed only)
        eps = 0.005
        q_gain_up = ((exp["methyl_change"] > eps) & (exp["LFC_T2vsT1"] > 0)).sum()
        q_loss_up = ((exp["methyl_change"] < -eps) & (exp["LFC_T2vsT1"] > 0)).sum()
        q_loss_down = ((exp["methyl_change"] < -eps) & (exp["LFC_T2vsT1"] < 0)).sum()
        q_gain_down = ((exp["methyl_change"] > eps) & (exp["LFC_T2vsT1"] < 0)).sum()
        n_zero = (exp["methyl_change"].abs() <= eps).sum()

        summary[motif] = {
            "rho_exposed": rho_e, "p_exposed": p_e, "n_exposed": n_e,
            "rho_shielded": rho_s, "p_shielded": p_s, "n_shielded": n_s,
            "exp_gain_up": int(q_gain_up), "exp_loss_up": int(q_loss_up),
            "exp_loss_down": int(q_loss_down), "exp_gain_down": int(q_gain_down),
            "exp_zero_change": int(n_zero),
        }

        ax.scatter(shld["methyl_change"], shld["LFC_T2vsT1"], c="lightgray", s=12,
                   alpha=0.55, edgecolors="none", label=f"Shielded (n={len(shld)})", zorder=1)
        ax.scatter(exp["methyl_change"], exp["LFC_T2vsT1"], c=color, s=40,
                   alpha=0.85, edgecolors="black", linewidths=0.4,
                   label=f"Exposed (n={len(exp)})", zorder=3)
        ax.axhline(0, color="gray", lw=0.6)
        ax.axvline(0, color="gray", lw=0.6)

        if len(exp) >= 5 and exp["methyl_change"].std() > 0:
            slope, intercept = np.polyfit(exp["methyl_change"], exp["LFC_T2vsT1"], 1)
            xs = np.linspace(exp["methyl_change"].min(), exp["methyl_change"].max(), 100)
            ax.plot(xs, slope * xs + intercept, color=color, lw=1.4, ls="--", alpha=0.8)

        ax.set_xlabel(f"{motif} methylation change\n(T2 - T1, TSS +/- {window} bp)")
        if ax is axes[0]:
            ax.set_ylabel("Expression log$_2$FC ($T_2$ vs $T_1$)")

        # Subtitle with stats
        ttl = (f"{motif}\n"
               f"Exposed $\\rho$={rho_e:.3f}, p={p_e:.3g}, n={n_e}"
               f"  |  Shielded $\\rho$={rho_s:.3f}, p={p_s:.3g}, n={n_s}\n"
               f"Exposed quadrants  gain+up={q_gain_up}, loss+up={q_loss_up},"
               f" loss+down={q_loss_down}, gain+down={q_gain_down},"
               f" |$\\Delta$|<={eps}: {n_zero}")
        ax.set_title(ttl, fontsize=9, loc="left")
        ax.legend(loc="lower right", fontsize=8)
        ax.grid(True, linestyle=":", alpha=0.4)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.suptitle("Motif decomposition: TSS-proximal methylation change vs T1->T2 LFC",
                 fontsize=11, y=1.0)
    fig.tight_layout()

    stem = "exposed_TF_motif_decomposition"
    for d in (OUT_DIR_LOCAL, OUT_DIR_TOP):
        for ext in ("pdf", "svg", "png"):
            p = d / f"{stem}.{ext}"
            fig.savefig(p, dpi=200 if ext == "png" else None, bbox_inches="tight")
            logger.info(f"  Saved: {p}")
    plt.close(fig)

    return summary


# ============================================================================
# Task 2: Window sensitivity
# ============================================================================
def task2_window_sensitivity(universe: pd.DataFrame, methyl_full: pd.DataFrame,
                              deseq: pd.DataFrame, windows: list[int]) -> pd.DataFrame:
    logger.info("\n=== Task 2: Window sensitivity ===")

    rows = []
    methyl_combined = methyl_full[methyl_full["motif"].isin(["GCCGGC", "AAGCCCG"])]
    methyl_all = methyl_full  # includes "other" sites too
    for window in windows:
        for label, m in (("GCCGGC+AAGCCCG", methyl_combined), ("All_4mC+6mA", methyl_all)):
            df = build_scores(universe, m, window)
            df = add_lfc(df, deseq)
            exp = df[df.is_exposed == 1]
            shld = df[df.is_exposed == 0]
            rho_e, p_e, n_e = spearman(exp)
            rho_s, p_s, n_s = spearman(shld)
            rows.append({
                "window_bp": window, "site_set": label,
                "rho_exposed": rho_e, "p_exposed": p_e, "n_exposed": n_e,
                "rho_shielded": rho_s, "p_shielded": p_s, "n_shielded": n_s,
            })
            logger.info(f"  +/-{window}bp [{label}] Exposed rho={rho_e:.4f} (p={p_e:.4g}) | "
                        f"Shielded rho={rho_s:.4f} (p={p_s:.4g})")

    out = pd.DataFrame(rows)
    out.to_csv(TABLE_DIR / "exposed_TF_window_sensitivity.tsv", sep="\t", index=False)
    logger.info(f"  Wrote table: {TABLE_DIR / 'exposed_TF_window_sensitivity.tsv'}")

    # Plot rho vs window
    fig, ax = plt.subplots(figsize=(8, 5))
    for label, marker in (("GCCGGC+AAGCCCG", "o"), ("All_4mC+6mA", "s")):
        sub = out[out.site_set == label]
        ax.errorbar(sub.window_bp, sub.rho_exposed, fmt=f"{marker}-", color="#C0392B",
                    label=f"Exposed [{label}]", lw=1.6, ms=8)
        ax.errorbar(sub.window_bp, sub.rho_shielded, fmt=f"{marker}--", color="gray",
                    label=f"Shielded [{label}]", lw=1.2, ms=6)

    ax.axhline(0, color="black", lw=0.5)
    ax.set_xlabel("TSS window (+/- bp)")
    ax.set_ylabel("Spearman $\\rho$ (methyl change vs LFC T2-T1)")
    ax.set_title("Window sensitivity of methyl-LFC correlation", fontsize=11)
    ax.set_xticks(windows)
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, linestyle=":", alpha=0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()

    stem = "exposed_TF_window_sensitivity"
    for d in (OUT_DIR_LOCAL, OUT_DIR_TOP):
        for ext in ("pdf", "svg", "png"):
            p = d / f"{stem}.{ext}"
            fig.savefig(p, dpi=200 if ext == "png" else None, bbox_inches="tight")
            logger.info(f"  Saved: {p}")
    plt.close(fig)

    return out


def main() -> None:
    OUT_DIR_LOCAL.mkdir(parents=True, exist_ok=True)
    OUT_DIR_TOP.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    logger.info("Loading universe and methylation ...")
    universe = load_universe()
    logger.info(f"  TF universe: n={len(universe)} (Exposed={universe.is_exposed.sum()})")

    methyl = load_methyl_with_motif()
    deseq = pd.read_csv(DESEQ_T2, sep="\t")

    # Task 1: motif decomposition (default window 500)
    s1 = task1_motif_decomposition(universe, methyl, deseq, window=500)
    logger.info(f"\nTask 1 summary: {s1}")

    # Task 2: window sensitivity
    s2 = task2_window_sensitivity(universe, methyl, deseq, windows=[200, 293, 500, 1000])
    logger.info("\nTask 2 summary table:")
    logger.info(s2.to_string(index=False))


if __name__ == "__main__":
    main()
