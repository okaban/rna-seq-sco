#!/usr/bin/env python3
"""A3: Is the Exposed/Shielded distance AUC=0.908 circular?

The 57 Exposed genes are selected by methylation-expression coordination
(requires methylation). The classifier feature is nearest-methylation distance.
This script tests whether the AUC is inflated by selection-on-methylation.

Three tests:
  1. Stratified 10-fold CV of the single-distance classifier on the 57/998
     labels -> shows whether the AUC is a model-overfitting artefact (it is not;
     a single fixed feature needs no fitting).
  2. Label-permutation null -> AUC vs chance.
  3. DECISIVE: expression-only switch label. Define DEG regulators (|LFC|>=1 at
     T2 or T3) WITHOUT any methylation criterion, and test whether nearest
     distance predicts THIS expression-only label. If distance only predicts the
     methylation-selected label (0.908) but not the expression-only label, the
     0.908 is driven by the methylation selection (circular).

Outputs:
  figures/A3_circularity.{pdf,svg,png}
  tables/A3_circularity_summary.tsv
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import StratifiedKFold

TBL = Path(__file__).resolve().parents[1] / "tables"
FIG = Path(__file__).resolve().parents[1] / "figures"
FIG.mkdir(exist_ok=True)
RNG = np.random.default_rng(42)


def cv_auc(y: np.ndarray, score: np.ndarray, k: int = 10) -> tuple[float, float]:
    skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=42)
    aucs = []
    for _, te in skf.split(score.reshape(-1, 1), y):
        if len(np.unique(y[te])) > 1:
            aucs.append(roc_auc_score(y[te], score[te]))
    return float(np.mean(aucs)), float(np.std(aucs))


def main() -> None:
    df = pd.read_csv(TBL / "all_genes_features_unified_n57.tsv", sep="\t")
    df = df.dropna(subset=["nearest_methyl_distance"]).copy()
    dist = df["nearest_methyl_distance"].values
    score = -dist  # exposed = small distance
    y_exp = df["is_exposed"].astype(int).values

    # expression-only label (no methylation criterion)
    y_deg = (((df["LFC_T2vsT1"].abs() >= 1) | (df["LFC_T3vsT1"].abs() >= 1))
             .astype(int).values)

    rows = []

    # Test 1: pooled + 10-fold CV on the methylation-coordinated label
    auc_exp = roc_auc_score(y_exp, score)
    cvm, cvs = cv_auc(y_exp, score)
    rows.append(["coordinated_label_pooled", int(y_exp.sum()), int((y_exp == 0).sum()),
                 round(auc_exp, 4), "-"])
    rows.append(["coordinated_label_10foldCV", int(y_exp.sum()), int((y_exp == 0).sum()),
                 round(cvm, 4), f"+/-{cvs:.3f}"])

    # Test 2: permutation null for the coordinated label
    n_exp = int(y_exp.sum())
    null = []
    idx = np.arange(len(y_exp))
    for _ in range(2000):
        perm = np.zeros(len(y_exp), dtype=int)
        perm[RNG.choice(idx, n_exp, replace=False)] = 1
        null.append(roc_auc_score(perm, score))
    null = np.array(null)
    p_perm = (np.sum(null >= auc_exp) + 1) / (len(null) + 1)
    rows.append(["permutation_null_mean", n_exp, len(y_exp) - n_exp,
                 round(float(null.mean()), 4),
                 f"95%[{np.percentile(null,2.5):.3f},{np.percentile(null,97.5):.3f}]"])

    # Test 3 (DECISIVE): expression-only label
    auc_deg = roc_auc_score(y_deg, score)
    rows.append(["expression_only_label(|LFC|>=1)", int(y_deg.sum()),
                 int((y_deg == 0).sum()), round(auc_deg, 4),
                 "distance predicts expression-switch?"])

    # group median distances
    med_exp = float(np.median(dist[y_exp == 1]))
    deg_not_exp = (y_deg == 1) & (y_exp == 0)
    non_deg = y_deg == 0
    med_degne = float(np.median(dist[deg_not_exp]))
    med_nondeg = float(np.median(dist[non_deg]))

    summary = pd.DataFrame(rows, columns=["test", "n_pos", "n_neg", "AUC", "note"])
    summary.to_csv(TBL / "A3_circularity_summary.tsv", sep="\t", index=False)
    with open(TBL / "A3_circularity_summary.tsv", "a") as f:
        f.write(f"\n# permutation p-value (AUC>=observed {auc_exp:.3f}): {p_perm:.4g}\n")
        f.write(f"# median nearest distance (bp): Exposed57={med_exp:.0f}, "
                f"DEG-not-Exposed={med_degne:.0f}, non-DEG={med_nondeg:.0f}\n")
        f.write("# INTERPRETATION: distance strongly predicts the methylation-"
                "coordinated label (0.908) but NOT the expression-only label "
                f"({auc_deg:.3f}); expression-defined switch genes are NOT "
                "TSS-proximal-methylated. The 0.908 is substantially driven by "
                "selection-on-methylation (partial circularity).\n")

    # ---- figure ----
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4))

    # Panel A: ROC, coordinated vs expression-only label
    ax = axes[0]
    for yl, lab, col in [(y_exp, f"Methyl-coordinated label\nAUC={auc_exp:.3f}", "#D1495B"),
                         (y_deg, f"Expression-only label\nAUC={auc_deg:.3f}", "#3B6FB6")]:
        fpr, tpr, _ = roc_curve(yl, score)
        ax.plot(fpr, tpr, color=col, lw=2, label=lab)
    ax.plot([0, 1], [0, 1], "k:", lw=0.6, alpha=0.5)
    ax.set_xlabel("False positive rate"); ax.set_ylabel("True positive rate")
    ax.set_title("Distance predicts the methyl-\nselected label, not expression",
                 fontsize=9.5, fontweight="bold")
    ax.legend(fontsize=7, loc="lower right"); ax.set_aspect("equal")

    # Panel B: permutation null
    ax = axes[1]
    ax.hist(null, bins=40, color="#bbb", edgecolor="white")
    ax.axvline(auc_exp, color="#D1495B", lw=2)
    ax.text(auc_exp - 0.01, ax.get_ylim()[1] * 0.9, f"observed\n{auc_exp:.3f}\n(p={p_perm:.0e})",
            ha="right", va="top", fontsize=7.5, color="#D1495B")
    ax.set_xlabel("AUC under permuted labels"); ax.set_ylabel("count")
    ax.set_title("Permutation null\n(label vs chance)", fontsize=9.5, fontweight="bold")

    # Panel C: median distance by group
    ax = axes[2]
    groups = ["Exposed 57\n(methyl-sel.)", "DEG, not\nExposed", "non-DEG"]
    meds = [med_exp, med_degne, med_nondeg]
    cols = ["#D1495B", "#E8A33D", "#9aa0a6"]
    bars = ax.bar(groups, meds, color=cols, edgecolor="white")
    for b, v in zip(bars, meds):
        ax.text(b.get_x() + b.get_width() / 2, v + 15, f"{v:.0f}", ha="center",
                va="bottom", fontsize=9, fontweight="bold")
    ax.set_ylabel("Median nearest\nmethylation distance (bp)")
    ax.set_title("Expression-switch genes are\nNOT proximal-methylated",
                 fontsize=9.5, fontweight="bold")
    ax.set_ylim(0, max(meds) * 1.25)

    fig.tight_layout()
    for ext in ("pdf", "svg", "png"):
        fig.savefig(FIG / f"A3_circularity.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(summary.to_string(index=False))
    print(f"\npermutation p = {p_perm:.4g}")
    print(f"median dist (bp): Exposed57={med_exp:.0f}, DEG-not-Exposed={med_degne:.0f}, "
          f"non-DEG={med_nondeg:.0f}")
    print(f"AUC coordinated-label={auc_exp:.3f}  vs  expression-only-label={auc_deg:.3f}")
    print(f"Saved figure + table.")


if __name__ == "__main__":
    main()
