#!/usr/bin/env python3
"""
Exposed TF: TSS-proximal methylation change vs T1->T2 log2FC, stratified by motif.

Splits the 57 Exposed TFs into three motif strata based on TSS +/- 500 bp
high-confidence methylation sites:
  - GCCGGC      : 4mC sites whose surrounding sequence contains GCCGGC
  - AAGCCCG     : sites (4mC or 6mA) whose surrounding sequence contains AAGCCCG
  - unassigned_6mA : 6mA sites not in AAGCCCG context

Per-gene methylation change for each stratum:
  delta = mean(weighted_mod_freq at T2) - mean(weighted_mod_freq at T1)
  (computed over sites in that motif within TSS +/- 500 bp; NaN if no sites)

Outputs:
  - analysis/figures/exposed_TF_methyl_motif_stratified.{pdf,png}
  - analysis/57_temporal_dynamics_exposed_TF/figures/exposed_TF_methyl_motif_stratified.{pdf,png}
  - analysis/57_temporal_dynamics_exposed_TF/results/motif_stratified_stats.txt
  - analysis/57_temporal_dynamics_exposed_TF/tables/exposed_TF_methyl_motif_stratified.tsv
"""

from __future__ import annotations

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
SEQ_FILE = EPI / "07_motif_analysis/methylation_site_sequences.csv"
DESEQ_T2 = BASE / "04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_2_vs_1.tsv"

OUT_DIR_LOCAL = EPI / "57_temporal_dynamics_exposed_TF/figures"
OUT_DIR_TOP = EPI / "figures"
OUT_TABLE = EPI / "57_temporal_dynamics_exposed_TF/tables/exposed_TF_methyl_motif_stratified.tsv"
OUT_STATS = EPI / "57_temporal_dynamics_exposed_TF/results/motif_stratified_stats.txt"

WINDOW = 500

MOTIF_GCCGGC = "GCCGGC"
MOTIF_AAGCCCG = "AAGCCCG"
MOTIF_UNASSIGNED = "unassigned_6mA"

PANEL_ORDER = [MOTIF_GCCGGC, MOTIF_AAGCCCG, MOTIF_UNASSIGNED]
PANEL_TITLES = {
    MOTIF_GCCGGC: "GCCGGC (4mC)",
    MOTIF_AAGCCCG: "AAGCCCG (4mC + 6mA)",
    MOTIF_UNASSIGNED: "Unassigned 6mA",
}
PANEL_COLORS = {
    MOTIF_GCCGGC: "#1F77B4",
    MOTIF_AAGCCCG: "#C0392B",
    MOTIF_UNASSIGNED: "#2CA02C",
}


def classify_motif(row: pd.Series) -> str:
    """Assign each site to one of three motif strata or 'other'."""
    seq = row.get("sequence", "")
    if not isinstance(seq, str) or len(seq) == 0:
        return "other"
    has_aagcccg = "AAGCCCG" in seq
    has_gccggc = "GCCGGC" in seq
    if has_aagcccg:
        return MOTIF_AAGCCCG
    if row["mod_type"] == "4mC" and has_gccggc:
        return MOTIF_GCCGGC
    if row["mod_type"] == "6mA":
        return MOTIF_UNASSIGNED
    return "other"


def annotate_motif(methyl: pd.DataFrame, seq: pd.DataFrame) -> pd.DataFrame:
    """Join motif label onto methylation sites using site-level sequence file."""
    seq_unique = (
        seq.drop_duplicates(subset=["chrom", "position", "strand", "mod_type"])[
            ["chrom", "position", "strand", "mod_type", "sequence"]
        ]
        .copy()
    )
    seq_unique["motif"] = seq_unique.apply(classify_motif, axis=1)
    merged = methyl.merge(
        seq_unique[["chrom", "position", "strand", "mod_type", "motif"]],
        on=["chrom", "position", "strand", "mod_type"],
        how="left",
    )
    return merged


def per_gene_mean_change(
    methyl: pd.DataFrame, tss: int, motif: str, window: int = WINDOW
) -> tuple[float, int, int]:
    """Mean-based methylation change for a single gene within a motif stratum.

    Returns (delta, n_sites_T1, n_sites_T2). delta is np.nan if no sites in
    either timepoint for that motif within the TSS window.
    """
    sub = methyl.loc[
        (methyl["motif"] == motif)
        & (methyl["position"] >= tss - window)
        & (methyl["position"] <= tss + window)
    ]
    t1 = sub.loc[sub["timepoint"] == "T1", "weighted_mod_freq"]
    t2 = sub.loc[sub["timepoint"] == "T2", "weighted_mod_freq"]
    if len(t1) == 0 and len(t2) == 0:
        return np.nan, 0, 0
    mean_t1 = t1.mean() if len(t1) > 0 else 0.0
    mean_t2 = t2.mean() if len(t2) > 0 else 0.0
    return float(mean_t2 - mean_t1), int(len(t1)), int(len(t2))


def _tss_from_start_end(start: int, end: int, strand: str) -> int:
    return int(start) if strand == "+" else int(end)


def main() -> None:
    logger.info("Loading inputs ...")
    exposed_df = pd.read_csv(EXPOSED_FILE, sep="\t")
    all_features = pd.read_csv(ALL_FEATURES, sep="\t")
    methyl = pd.read_csv(METHYL_FILE)
    seq = pd.read_csv(SEQ_FILE)
    deseq = pd.read_csv(DESEQ_T2, sep="\t")

    logger.info(f"  Exposed TFs: n={len(exposed_df)}")
    logger.info(f"  Methylation sites (rows): n={len(methyl)}")
    logger.info(f"  Sequence-annotated sites (rows): n={len(seq)}")

    methyl = methyl.loc[methyl["timepoint"].isin(["T1", "T2"])].copy()
    methyl = annotate_motif(methyl, seq)

    motif_counts = methyl["motif"].value_counts(dropna=False).to_dict()
    logger.info(f"  Methylation site motif distribution (T1+T2 rows): {motif_counts}")

    # Build TSS table for the 57 exposed TFs (prefer experimental TSS)
    h29_tss = dict(zip(all_features["locus_tag"], all_features["tss"]))
    tss_rows = []
    for _, r in exposed_df.iterrows():
        lt = r["locus_tag"]
        if lt in h29_tss and not pd.isna(h29_tss[lt]):
            tss_val = int(h29_tss[lt])
            src = "experimental"
        else:
            tss_val = _tss_from_start_end(r["start"], r["end"], r["strand"])
            src = "gene_start"
        tss_rows.append(
            {
                "locus_tag": lt,
                "gene_name": r.get("gene_name", ""),
                "tf_family": r.get("tf_family", ""),
                "tss": tss_val,
                "tss_source": src,
            }
        )
    universe = pd.DataFrame(tss_rows)
    logger.info(f"  Built TSS table for {len(universe)} exposed TFs")

    # Per-gene per-motif methylation change
    rows = []
    for _, g in universe.iterrows():
        rec = {
            "locus_tag": g["locus_tag"],
            "gene_name": g["gene_name"],
            "tf_family": g["tf_family"],
            "tss": int(g["tss"]),
            "tss_source": g["tss_source"],
        }
        for motif in PANEL_ORDER:
            d, n1, n2 = per_gene_mean_change(methyl, int(g["tss"]), motif)
            rec[f"{motif}_delta"] = d
            rec[f"{motif}_n_T1"] = n1
            rec[f"{motif}_n_T2"] = n2
        rows.append(rec)
    df = pd.DataFrame(rows)

    # Merge T2 vs T1 LFC
    deseq = deseq.rename(
        columns={"gene_id": "locus_tag", "log2FoldChange": "LFC_T2vsT1"}
    )[["locus_tag", "LFC_T2vsT1", "padj"]]
    df = df.merge(deseq, on="locus_tag", how="left")
    df_with_lfc = df.dropna(subset=["LFC_T2vsT1"]).reset_index(drop=True)
    logger.info(
        f"  Exposed TFs with DESeq2 LFC available: n={len(df_with_lfc)} / {len(df)}"
    )

    OUT_TABLE.parent.mkdir(parents=True, exist_ok=True)
    df_with_lfc.to_csv(OUT_TABLE, sep="\t", index=False)
    logger.info(f"  Wrote table: {OUT_TABLE}")

    # ---- Stats per motif ----
    stat_lines: list[str] = []
    stat_lines.append("Exposed TF (n=57) motif-stratified methylation change vs LFC (T2 vs T1)")
    stat_lines.append("=" * 78)
    stat_lines.append(
        f"Window: TSS +/- {WINDOW} bp | delta = mean(weighted_mod_freq, T2) - mean(T1)"
    )
    stat_lines.append("Genes lacking sites for a given motif are excluded from that panel.")
    stat_lines.append("")

    panel_data: dict[str, pd.DataFrame] = {}
    for motif in PANEL_ORDER:
        col = f"{motif}_delta"
        sub = df_with_lfc.dropna(subset=[col]).copy()
        panel_data[motif] = sub
        if len(sub) >= 3:
            rho, pval = stats.spearmanr(sub[col], sub["LFC_T2vsT1"])
        else:
            rho, pval = np.nan, np.nan
        stat_lines.append(f"[{PANEL_TITLES[motif]}]")
        stat_lines.append(f"  n_genes_with_sites = {len(sub)}")
        stat_lines.append(f"  Spearman rho       = {rho:.4f}" if not np.isnan(rho) else "  Spearman rho       = NA")
        stat_lines.append(f"  p-value            = {pval:.4g}" if not np.isnan(pval) else "  p-value            = NA")
        if len(sub) > 0:
            stat_lines.append(
                f"  delta range        = [{sub[col].min():.3f}, {sub[col].max():.3f}] (units: %% mod freq)"
            )
            stat_lines.append(
                f"  LFC range          = [{sub['LFC_T2vsT1'].min():.3f}, {sub['LFC_T2vsT1'].max():.3f}]"
            )
        stat_lines.append("")

        logger.info(
            f"  {motif}: n={len(sub)}, rho={rho:.4f}, p={pval:.4g}"
            if not np.isnan(rho)
            else f"  {motif}: n={len(sub)} (insufficient for correlation)"
        )

    OUT_STATS.parent.mkdir(parents=True, exist_ok=True)
    OUT_STATS.write_text("\n".join(stat_lines))
    logger.info(f"  Wrote stats: {OUT_STATS}")

    # ---- Plot: 3 panels ----
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    for ax, motif in zip(axes, PANEL_ORDER):
        col = f"{motif}_delta"
        sub = panel_data[motif]
        color = PANEL_COLORS[motif]

        if len(sub) >= 3:
            rho, pval = stats.spearmanr(sub[col], sub["LFC_T2vsT1"])
        else:
            rho, pval = np.nan, np.nan

        ax.scatter(
            sub[col],
            sub["LFC_T2vsT1"],
            c=color,
            s=46,
            alpha=0.85,
            edgecolors="black",
            linewidths=0.4,
            zorder=3,
        )

        ax.axhline(0, color="gray", lw=0.8, ls=":", zorder=1)
        ax.axvline(0, color="gray", lw=0.8, ls=":", zorder=1)

        # Stat annotation
        if not np.isnan(rho):
            stat_text = f"Spearman $\\rho$ = {rho:.3f}\np = {pval:.3g}\nn = {len(sub)}"
        else:
            stat_text = f"n = {len(sub)} (too few)"
        ax.text(
            0.03,
            0.97,
            stat_text,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=10,
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=color, lw=0.8, alpha=0.9),
        )

        ax.set_title(PANEL_TITLES[motif], fontsize=12, color=color, fontweight="bold")
        ax.set_xlabel(
            r"Methylation change $\Delta$ (% mod freq, $T_2 - T_1$)" + "\nmean over motif sites in TSS $\\pm$500 bp",
            fontsize=10,
        )
        ax.set_ylabel("Expression log$_2$FC ($T_2$ vs $T_1$, DESeq2)", fontsize=10)
        ax.grid(True, linestyle=":", alpha=0.35)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.suptitle(
        "Exposed TFs (n=57): TSS-proximal methylation change vs expression change, stratified by motif",
        fontsize=13,
        y=1.02,
    )
    fig.tight_layout()

    # Save
    OUT_DIR_LOCAL.mkdir(parents=True, exist_ok=True)
    OUT_DIR_TOP.mkdir(parents=True, exist_ok=True)
    stem = "exposed_TF_methyl_motif_stratified"
    for d in (OUT_DIR_TOP, OUT_DIR_LOCAL):
        for ext in ("pdf", "png"):
            p = d / f"{stem}.{ext}"
            fig.savefig(p, dpi=200 if ext == "png" else None, bbox_inches="tight")
            logger.info(f"  Saved: {p}")

    plt.close(fig)
    logger.info("Done.")


if __name__ == "__main__":
    main()
