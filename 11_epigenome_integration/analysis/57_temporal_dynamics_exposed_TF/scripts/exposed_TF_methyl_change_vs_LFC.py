#!/usr/bin/env python3
"""
Exposed TF: TSS-proximal methylation change vs T1->T2 log2FC scatter.

Question: Is protection-zone collapse (loss of TSS-proximal methylation) at T2
quantitatively associated with transcriptional derepression of the 57 Exposed
TFs?

Design:
  - X = TSS-proximal methylation score (T2) - score (T1)
        score = sum of weighted_mod_freq (units of fraction, 0-1) for all
        4mC + 6mA high-confidence sites within TSS +/- 500 bp
  - Y = DESeq2 log2FC, T2 vs T1
  - Foreground: 57 Exposed TFs (red)
  - Background: remaining TFs (Shielded, gray)
  - Spearman rho + p reported for Exposed and Shielded separately

Outputs (dir 57 figures + top-level figures):
  - exposed_TF_methyl_change_vs_LFC.pdf / .svg / .png
  - exposed_TF_methyl_change_vs_LFC_table.tsv
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

OUT_DIR_LOCAL = EPI / "57_temporal_dynamics_exposed_TF/figures"
OUT_DIR_TOP = EPI / "figures"
OUT_TABLE = EPI / "57_temporal_dynamics_exposed_TF/tables/exposed_TF_methyl_change_vs_LFC.tsv"

WINDOW = 500  # +/- bp around TSS


def methyl_score(methyl_df: pd.DataFrame, tss: int, window: int = WINDOW) -> tuple[float, float, int, int]:
    """Sum of weighted_mod_freq (as fraction) within TSS +/- window for T1, T2.

    Returns (score_T1, score_T2, n_sites_T1, n_sites_T2).
    """
    in_window = (methyl_df["position"] >= tss - window) & (methyl_df["position"] <= tss + window)
    sub = methyl_df.loc[in_window]
    t1 = sub.loc[sub["timepoint"] == "T1", "weighted_mod_freq"]
    t2 = sub.loc[sub["timepoint"] == "T2", "weighted_mod_freq"]
    return t1.sum() / 100.0, t2.sum() / 100.0, len(t1), len(t2)


def _tss_from_start_end(start: int, end: int, strand: str) -> int:
    """Gene-start-based TSS proxy: start (+ strand) or end (- strand)."""
    return int(start) if strand == "+" else int(end)


def main() -> None:
    logger.info("Loading inputs ...")
    exposed_df = pd.read_csv(EXPOSED_FILE, sep="\t")
    all_features = pd.read_csv(ALL_FEATURES, sep="\t")
    methyl = pd.read_csv(METHYL_FILE)
    deseq_t2 = pd.read_csv(DESEQ_T2, sep="\t")

    exposed_set = set(exposed_df["locus_tag"].dropna().tolist())
    logger.info(f"  Exposed TF set (H28): n={len(exposed_set)}")
    logger.info(f"  All TFs in H29 features table: n={len(all_features)}")
    logger.info(f"  Methylation sites (rows): n={len(methyl)}")

    # Restrict methyl to T1/T2
    methyl = methyl.loc[methyl["timepoint"].isin(["T1", "T2"])].copy()

    # --- Build TF universe with TSS ---
    # Foreground: 57 Exposed from H28. Use experimental TSS if present in H29,
    # else fall back to gene-start (strand-aware) as TSS proxy.
    h29_tss_by_locus = dict(zip(all_features["locus_tag"], all_features["tss"]))

    universe_rows = []
    # Add 57 exposed TFs first
    for _, r in exposed_df.iterrows():
        lt = r["locus_tag"]
        if lt in h29_tss_by_locus and not pd.isna(h29_tss_by_locus[lt]):
            tss = int(h29_tss_by_locus[lt])
            tss_source = "experimental"
        else:
            tss = _tss_from_start_end(r["start"], r["end"], r["strand"])
            tss_source = "gene_start"
        universe_rows.append(
            {
                "locus_tag": lt,
                "gene_name": r.get("gene_name", ""),
                "tf_family": r.get("tf_family", ""),
                "tss": tss,
                "tss_source": tss_source,
                "is_exposed": 1,
            }
        )

    # Add shielded TFs from H29 (those NOT in exposed set)
    for _, r in all_features.iterrows():
        lt = r["locus_tag"]
        if lt in exposed_set:
            continue  # already added as exposed
        if pd.isna(r["tss"]):
            continue
        universe_rows.append(
            {
                "locus_tag": lt,
                "gene_name": r.get("gene_name", ""),
                "tf_family": r.get("tf_family", ""),
                "tss": int(r["tss"]),
                "tss_source": "experimental",
                "is_exposed": 0,
            }
        )

    universe = pd.DataFrame(universe_rows)
    logger.info(
        f"  TF universe: n={len(universe)} "
        f"(Exposed={universe.is_exposed.sum()}, Shielded={(universe.is_exposed == 0).sum()})"
    )
    logger.info(
        f"  Exposed TSS sources: {universe.loc[universe.is_exposed == 1, 'tss_source'].value_counts().to_dict()}"
    )

    # --- Methylation score per gene ---
    rows = []
    for _, r in universe.iterrows():
        t1_score, t2_score, n_t1, n_t2 = methyl_score(methyl, int(r["tss"]))
        rows.append(
            {
                "locus_tag": r["locus_tag"],
                "gene_name": r["gene_name"],
                "tf_family": r["tf_family"],
                "tss": r["tss"],
                "tss_source": r["tss_source"],
                "is_exposed": r["is_exposed"],
                "methyl_score_T1": t1_score,
                "methyl_score_T2": t2_score,
                "methyl_change": t2_score - t1_score,
                "n_sites_T1": n_t1,
                "n_sites_T2": n_t2,
            }
        )
    df = pd.DataFrame(rows)
    logger.info(f"  Built per-TF methylation scores: n={len(df)}")

    # Merge canonical DESeq2 LFC
    deseq_t2 = deseq_t2.rename(columns={"gene_id": "locus_tag", "log2FoldChange": "LFC_T2vsT1"})[
        ["locus_tag", "LFC_T2vsT1", "padj"]
    ]
    df = df.merge(deseq_t2, on="locus_tag", how="left")

    # Drop rows missing LFC
    df = df.dropna(subset=["LFC_T2vsT1"]).reset_index(drop=True)
    logger.info(f"  Rows with LFC available: n={len(df)}")

    # Save full table
    OUT_TABLE.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_TABLE, sep="\t", index=False)
    logger.info(f"  Wrote table: {OUT_TABLE}")

    # Split groups
    exp = df[df["is_exposed"] == 1].copy()
    shld = df[df["is_exposed"] == 0].copy()
    logger.info(f"  Exposed n={len(exp)}, Shielded n={len(shld)}")

    # Spearman rho
    rho_e, p_e = stats.spearmanr(exp["methyl_change"], exp["LFC_T2vsT1"])
    rho_s, p_s = stats.spearmanr(shld["methyl_change"], shld["LFC_T2vsT1"])
    logger.info(f"  Spearman Exposed:  rho = {rho_e:.4f}, p = {p_e:.4g}, n = {len(exp)}")
    logger.info(f"  Spearman Shielded: rho = {rho_s:.4f}, p = {p_s:.4g}, n = {len(shld)}")

    # Among exposed-with-methyl-change != 0 (i.e., genes where TSS-proximal methylation actually changed)
    exp_nz = exp[exp["methyl_change"].abs() > 1e-9]
    if len(exp_nz) >= 5:
        rho_e_nz, p_e_nz = stats.spearmanr(exp_nz["methyl_change"], exp_nz["LFC_T2vsT1"])
        logger.info(f"  Spearman Exposed (|delta|>0): rho = {rho_e_nz:.4f}, p = {p_e_nz:.4g}, n = {len(exp_nz)}")
    else:
        rho_e_nz, p_e_nz = np.nan, np.nan

    # ---- Plot ----
    fig, ax = plt.subplots(figsize=(7.5, 6.5))

    ax.scatter(
        shld["methyl_change"],
        shld["LFC_T2vsT1"],
        c="lightgray",
        s=14,
        alpha=0.55,
        edgecolors="none",
        label=f"Shielded TFs (n={len(shld)})",
        zorder=1,
    )
    ax.scatter(
        exp["methyl_change"],
        exp["LFC_T2vsT1"],
        c="#C0392B",
        s=46,
        alpha=0.85,
        edgecolors="black",
        linewidths=0.5,
        label=f"Exposed TFs (n={len(exp)})",
        zorder=3,
    )

    ax.axhline(0, color="gray", lw=0.6, zorder=0)
    ax.axvline(0, color="gray", lw=0.6, zorder=0)

    # Regression line for exposed (least-squares on ranks via OLS on values for visualization)
    if len(exp) >= 5:
        x = exp["methyl_change"].values
        y = exp["LFC_T2vsT1"].values
        m, b = np.polyfit(x, y, 1)
        xs = np.linspace(x.min(), x.max(), 100)
        ax.plot(xs, m * xs + b, color="#C0392B", lw=1.6, ls="--", alpha=0.8, zorder=2,
                label=f"Exposed OLS (slope={m:.2f})")

    # Quadrant counts (Exposed only, |dx|, |dy| > small epsilon for cleaner counts)
    eps = 0.01
    q1 = ((exp["methyl_change"] > eps) & (exp["LFC_T2vsT1"] > 0)).sum()  # gain x up
    q2 = ((exp["methyl_change"] < -eps) & (exp["LFC_T2vsT1"] > 0)).sum()  # loss x up
    q3 = ((exp["methyl_change"] < -eps) & (exp["LFC_T2vsT1"] < 0)).sum()  # loss x down
    q4 = ((exp["methyl_change"] > eps) & (exp["LFC_T2vsT1"] < 0)).sum()  # gain x down

    # Place quadrant labels in axes-fraction coordinates
    ax.text(0.97, 0.97, f"Methyl gain + LFC up\nn = {q1}", transform=ax.transAxes,
            ha="right", va="top", fontsize=9, color="#7B241C",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#C0392B", lw=0.6, alpha=0.85))
    ax.text(0.03, 0.97, f"Methyl loss + LFC up\nn = {q2}", transform=ax.transAxes,
            ha="left", va="top", fontsize=9, color="#7B241C",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#C0392B", lw=0.6, alpha=0.85))
    ax.text(0.03, 0.03, f"Methyl loss + LFC down\nn = {q3}", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=9, color="#7B241C",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#C0392B", lw=0.6, alpha=0.85))
    ax.text(0.97, 0.20, f"Methyl gain + LFC down\nn = {q4}", transform=ax.transAxes,
            ha="right", va="top", fontsize=9, color="#7B241C",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#C0392B", lw=0.6, alpha=0.85))

    ax.set_xlabel("TSS-proximal methylation change\n(score $T_2 - T_1$, sites within TSS $\\pm$500 bp)")
    ax.set_ylabel("Expression log$_2$FC ($T_2$ vs $T_1$, DESeq2)")

    title_lines = [
        f"Exposed TFs: Spearman $\\rho$ = {rho_e:.3f}, p = {p_e:.3g}, n = {len(exp)}",
        f"Shielded TFs: Spearman $\\rho$ = {rho_s:.3f}, p = {p_s:.3g}, n = {len(shld)}",
    ]
    ax.set_title("\n".join(title_lines), fontsize=10, loc="left")

    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 0.0), frameon=True,
              fontsize=8.5, framealpha=0.95)

    ax.grid(True, linestyle=":", alpha=0.4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()

    # Save
    OUT_DIR_LOCAL.mkdir(parents=True, exist_ok=True)
    OUT_DIR_TOP.mkdir(parents=True, exist_ok=True)
    stem = "exposed_TF_methyl_change_vs_LFC"
    for d in (OUT_DIR_LOCAL, OUT_DIR_TOP):
        for ext in ("pdf", "svg", "png"):
            p = d / f"{stem}.{ext}"
            fig.savefig(p, dpi=200 if ext == "png" else None, bbox_inches="tight")
            logger.info(f"  Saved: {p}")

    plt.close(fig)
    logger.info("Done.")


if __name__ == "__main__":
    main()
