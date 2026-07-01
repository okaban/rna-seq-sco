#!/usr/bin/env python3
"""Standalone heatmap of the 20 Category A exposed-TF cohort.

Rows  : Cat A genes, sorted by T1->T2 LFC (descending; up genes on top,
        down genes on the bottom).
Cols  : T1 / T2 / T3.
Values: row-wise z-score of mean log2(normalized count + 1).
Cmap  : RdBu_r (blue = low, white = mid, red = high).

A thin LFC_T2vsT1 colour-bar is drawn on the right of the main panel as a
visual annotation; this is optional and can be disabled with DRAW_LFC_BAR.

Outputs (300+ dpi):
  analysis/figures/FigCatA_heatmap_standalone.{pdf,png}
  Mirrored into /Users/okaban/obsidian/Research/rna-seq/Writing/

Window for site density / sorting matches the v2 main scatter (TSS +/- 293 bp),
though the heatmap itself does not consume the methylation signal.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
BASE = Path("/Users/okaban/bioinfo/rna-seq")
ANA = BASE / "11_epigenome_integration/analysis"

EXPOSED_TF = ANA / "58_AAGCCCG_exposed_TF_causal/tables/exposed_TF_methylation_status.tsv"
CATA_PATH = ANA / "57_temporal_dynamics_exposed_TF/data/categoryA_genes.tsv"
NOTABLE_PATH = ANA / "57_temporal_dynamics_exposed_TF/tables/exposed_TF_notable_list.tsv"
COUNTS_PATH = BASE / "04_deseq2/results/results/normalized_counts_M145.tsv"

FIG_DIR = ANA / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

OBSIDIAN_DIR = Path("/Users/okaban/obsidian/Research/rna-seq/Writing")

# Promoter window referenced in the figure caption (same as scatter v2).
WINDOW_BP = 293

# Toggle the right-hand LFC sidebar.
DRAW_LFC_BAR = False

# Sample -> timepoint assignment (3 reps per timepoint).
SAMPLE_TP = {
    "M145_1_1": "T1", "M145_1_2": "T1", "M145_1_3": "T1",
    "M145_2_1": "T2", "M145_2_3": "T2", "M145_2_4": "T2",
    "M145_3_2": "T3", "M145_3_3": "T3", "M145_3_4": "T3",
}

# --------------------------------------------------------------------------- #
# Style
# --------------------------------------------------------------------------- #
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

EXPR_CMAP = "RdBu_r"
LFC_CMAP = "PuOr_r"


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def resolve_label(row: pd.Series) -> str:
    name = str(row.get("gene_name", "")).strip()
    if name and name.lower() != "nan":
        return name
    old = str(row.get("old_locus_tag", "")).strip()
    if old and old.lower() != "nan":
        return old
    return str(row["locus_tag"])


def build_catA_expression(catA_set: set[str]) -> pd.DataFrame:
    """Long-form per-gene per-rep log2(normalized count + 1)."""
    counts = pd.read_csv(COUNTS_PATH, sep="\t")
    keep = counts[counts["gene_id"].isin(catA_set)].copy()
    rep_cols = [c for c in keep.columns if c.startswith("M145_")]
    log2 = np.log2(keep[rep_cols].astype(float) + 1.0)
    log2.insert(0, "locus_tag", keep["gene_id"].values)

    rows = []
    for col in rep_cols:
        tp = SAMPLE_TP[col]
        for _, row in log2.iterrows():
            rows.append({
                "locus_tag": row["locus_tag"],
                "timepoint": tp,
                "log2_norm": float(row[col]),
            })
    return pd.DataFrame(rows)


def build_zscore_matrix(expr_long: pd.DataFrame) -> pd.DataFrame:
    pivot = (
        expr_long.groupby(["locus_tag", "timepoint"])["log2_norm"]
        .mean()
        .unstack("timepoint")[["T1", "T2", "T3"]]
    )
    z = pivot.sub(pivot.mean(axis=1), axis=0).div(
        pivot.std(axis=1).replace(0, np.nan), axis=0
    ).fillna(0.0)
    return z


def order_rows_by_lfc(z: pd.DataFrame, lfc_map: dict[str, float]) -> list[str]:
    items = [(lt, lfc_map.get(lt, np.nan)) for lt in z.index]
    items.sort(key=lambda kv: (np.nan if pd.isna(kv[1]) else -kv[1], kv[0]))
    return [lt for lt, _ in items]


# --------------------------------------------------------------------------- #
# Plot
# --------------------------------------------------------------------------- #
def draw(ax_heat: plt.Axes,
         ax_cbar: plt.Axes,
         z_ord: pd.DataFrame,
         row_labels: list[str],
         lfc_values: np.ndarray,
         ax_lfc: plt.Axes | None,
         ax_lfc_cbar: plt.Axes | None) -> None:
    vmax = float(np.nanmax(np.abs(z_ord.values)))
    vmax = max(vmax, 1.5)
    im = ax_heat.imshow(
        z_ord.values, aspect="auto", cmap=EXPR_CMAP,
        vmin=-vmax, vmax=vmax, interpolation="nearest",
    )
    ax_heat.set_xticks([0, 1, 2])
    ax_heat.set_xticklabels(["T1", "T2", "T3"])
    ax_heat.set_yticks(np.arange(len(row_labels)))
    ax_heat.set_yticklabels(row_labels, fontsize=8.5, fontstyle="italic")
    ax_heat.tick_params(axis="x", length=0, pad=2)
    ax_heat.tick_params(axis="y", length=0, pad=2)
    for spine in ax_heat.spines.values():
        spine.set_visible(False)
    ax_heat.set_xlabel("Timepoint")
    ax_heat.set_title(
        f"Cat A expression dynamics (n={len(row_labels)}, z-score per gene)",
        fontsize=10,
    )

    cb = plt.colorbar(im, cax=ax_cbar, orientation="vertical")
    cb.set_label("z-score", fontsize=9)
    cb.ax.tick_params(labelsize=8, length=2)

    if ax_lfc is not None and ax_lfc_cbar is not None:
        lfc_col = lfc_values.reshape(-1, 1)
        finite = lfc_col[np.isfinite(lfc_col)]
        if finite.size:
            lfc_lim = float(np.nanmax(np.abs(finite)))
        else:
            lfc_lim = 1.0
        lfc_lim = max(lfc_lim, 1.0)
        im2 = ax_lfc.imshow(
            lfc_col, aspect="auto", cmap=LFC_CMAP,
            vmin=-lfc_lim, vmax=lfc_lim, interpolation="nearest",
        )
        ax_lfc.set_xticks([0])
        ax_lfc.set_xticklabels(["LFC\nT1→T2"], fontsize=8)
        ax_lfc.set_yticks([])
        ax_lfc.tick_params(axis="x", length=0, pad=2)
        for spine in ax_lfc.spines.values():
            spine.set_visible(False)

        cb2 = plt.colorbar(im2, cax=ax_lfc_cbar, orientation="vertical")
        cb2.set_label("LFC", fontsize=9)
        cb2.ax.tick_params(labelsize=8, length=2)


# --------------------------------------------------------------------------- #
# Save / mirror
# --------------------------------------------------------------------------- #
def save_pair(fig: plt.Figure, stem: str) -> tuple[Path, Path]:
    pdf = FIG_DIR / f"{stem}.pdf"
    png = FIG_DIR / f"{stem}.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {pdf}")
    print(f"  wrote {png}")
    return pdf, png


def mirror_to_obsidian(paths: list[Path]) -> None:
    if not OBSIDIAN_DIR.exists():
        print(f"  [skip] obsidian dir not present: {OBSIDIAN_DIR}")
        return
    for src in paths:
        dst = OBSIDIAN_DIR / src.name
        dst.write_bytes(src.read_bytes())
        print(f"  copied to {dst}")


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main() -> None:
    print("[load] Cat A gene set + LFC")
    catA = pd.read_csv(CATA_PATH, sep="\t")
    catA_set = set(catA["locus_tag"])
    print(f"  Cat A genes: {len(catA_set)}")

    tf = pd.read_csv(EXPOSED_TF, sep="\t")
    tf_catA = tf[tf["locus_tag"].isin(catA_set)][["locus_tag", "LFC_T2vsT1"]]
    lfc_map = dict(zip(tf_catA["locus_tag"], tf_catA["LFC_T2vsT1"]))

    print("[load] Cat A expression long-form")
    expr_long = build_catA_expression(catA_set)
    print(f"  expr rows: {len(expr_long)}")

    print("[build] z-score matrix and order rows by T1->T2 LFC (desc)")
    z = build_zscore_matrix(expr_long)
    order = order_rows_by_lfc(z, lfc_map)
    z_ord = z.loc[order]

    notable = pd.read_csv(NOTABLE_PATH, sep="\t")
    catA_meta = notable[notable["category"].fillna("").str.contains("A")].copy()
    label_lookup = catA_meta.set_index("locus_tag")
    # Locus tags with no legacy SCO equivalent are marked with '*'.
    NO_SCO_LOCUS = {"SC_RS36220"}
    row_labels = []
    for lt in z_ord.index:
        if lt in label_lookup.index:
            meta_row = label_lookup.loc[lt].copy()
            meta_row["locus_tag"] = lt
            label = resolve_label(meta_row)
        else:
            label = lt
        if lt in NO_SCO_LOCUS:
            label = f"{label}*"
        row_labels.append(label)

    lfc_values = np.array([lfc_map.get(lt, np.nan) for lt in z_ord.index],
                          dtype=float)
    n_up = int(np.sum(lfc_values > 0))
    n_dn = int(np.sum(lfc_values < 0))
    print(f"  ordered: up={n_up}  down={n_dn}  total={len(row_labels)}")

    print("[plot] standalone Cat A heatmap")
    if DRAW_LFC_BAR:
        fig = plt.figure(figsize=(4.6, 6.4))
        gs = GridSpec(
            nrows=1, ncols=4,
            width_ratios=[1.0, 0.10, 0.06, 0.06],
            wspace=0.18,
            figure=fig,
        )
        ax_heat = fig.add_subplot(gs[0, 0])
        ax_lfc = fig.add_subplot(gs[0, 1])
        ax_cbar = fig.add_subplot(gs[0, 2])
        ax_lfc_cbar = fig.add_subplot(gs[0, 3])
    else:
        fig = plt.figure(figsize=(4.0, 6.4))
        gs = GridSpec(
            nrows=1, ncols=2,
            width_ratios=[1.0, 0.06],
            wspace=0.10,
            figure=fig,
        )
        ax_heat = fig.add_subplot(gs[0, 0])
        ax_cbar = fig.add_subplot(gs[0, 1])
        ax_lfc = None
        ax_lfc_cbar = None

    draw(ax_heat, ax_cbar, z_ord, row_labels, lfc_values, ax_lfc, ax_lfc_cbar)

    pdf, png = save_pair(fig, "FigCatA_heatmap_standalone")
    mirror_to_obsidian([pdf, png])

    print("[done] standalone Cat A heatmap written under:", FIG_DIR)
    print(f"  promoter window referenced in caption: TSS +/- {WINDOW_BP} bp")


if __name__ == "__main__":
    main()
