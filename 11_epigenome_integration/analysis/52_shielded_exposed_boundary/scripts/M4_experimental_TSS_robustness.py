#!/usr/bin/env python3
"""M4 (reviewer self-review): Shielded/Exposed robustness under experimental TSS.

Question: Does the 293 bp distance classifier survive when the TSS proxy is
swapped from genome-annotation gene starts (predicted TSS, n=57/998) to
Jeong et al. (2016) dRNA-seq experimental TSSs (n=11/353)?

Two outputs:
  - figures/M4_experimental_TSS_robustness.{pdf,svg,png}
      Panel A: overlaid ROC (predicted vs experimental TSS)
      Panel B: canonical-57 member retention among testable genes
  - tables/M4_experimental_TSS_robustness.tsv  (exact numbers)

Reads:
  tables/all_genes_features_unified_n57.tsv  (predicted TSS, canonical 57/998)
  tables/all_genes_features.tsv              (experimental TSS, Jeong2016)
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score

TBL = Path(__file__).resolve().parents[1] / "tables"
FIG = Path(__file__).resolve().parents[1] / "figures"
FIG.mkdir(exist_ok=True)

RNG = np.random.default_rng(42)


def boot_ci(y: np.ndarray, score: np.ndarray, n: int = 2000) -> tuple[float, float]:
    """Bootstrap 95% CI of ROC AUC."""
    idx = np.arange(len(y))
    vals = []
    for _ in range(n):
        b = RNG.choice(idx, len(idx), replace=True)
        if len(np.unique(y[b])) > 1:
            vals.append(roc_auc_score(y[b], score[b]))
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def load(name: str) -> pd.DataFrame:
    df = pd.read_csv(TBL / name, sep="\t")
    return df.dropna(subset=["nearest_methyl_distance"]).copy()


def main() -> None:
    pred = load("all_genes_features_unified_n57.tsv")
    exp = load("all_genes_features.tsv")

    sets = {}
    for label, df in [("Predicted TSS", pred), ("Experimental TSS (Jeong2016)", exp)]:
        y = df["is_exposed"].astype(int).values
        score = -df["nearest_methyl_distance"].values
        auc = roc_auc_score(y, score)
        lo, hi = boot_ci(y, score)
        fpr, tpr, _ = roc_curve(y, score)
        sets[label] = dict(y=y, score=score, auc=auc, lo=lo, hi=hi,
                           fpr=fpr, tpr=tpr,
                           n_exp=int(y.sum()), n_shi=int((y == 0).sum()))

    # Member retention
    canon57 = set(pred.loc[pred["is_exposed"] == 1, "locus_tag"])
    exp_cov = set(exp["locus_tag"])
    exp_pos = set(exp.loc[exp["is_exposed"] == 1, "locus_tag"])
    covered = canon57 & exp_cov
    retained = covered & exp_pos
    n_cov, n_ret = len(covered), len(retained)

    # ---- table ----
    rows = [
        ["predicted_TSS", sets["Predicted TSS"]["n_exp"], sets["Predicted TSS"]["n_shi"],
         f"{sets['Predicted TSS']['auc']:.4f}",
         f"[{sets['Predicted TSS']['lo']:.3f}, {sets['Predicted TSS']['hi']:.3f}]"],
        ["experimental_TSS", sets["Experimental TSS (Jeong2016)"]["n_exp"],
         sets["Experimental TSS (Jeong2016)"]["n_shi"],
         f"{sets['Experimental TSS (Jeong2016)']['auc']:.4f}",
         f"[{sets['Experimental TSS (Jeong2016)']['lo']:.3f}, "
         f"{sets['Experimental TSS (Jeong2016)']['hi']:.3f}]"],
    ]
    tbl = pd.DataFrame(rows, columns=["tss_set", "n_exposed", "n_shielded", "AUC", "AUC_95CI"])
    tbl.to_csv(TBL / "M4_experimental_TSS_robustness.tsv", sep="\t", index=False)
    with open(TBL / "M4_experimental_TSS_robustness.tsv", "a") as f:
        f.write(f"\n# canonical Exposed (predicted TSS): {len(canon57)}\n")
        f.write(f"# canonical Exposed covered by experimental-TSS set ({len(exp_cov)} genes): {n_cov}\n")
        f.write(f"# of covered, retained Exposed under experimental TSS: {n_ret}/{n_cov}\n")
        f.write(f"# experimental-TSS-only Exposed (not in canonical 57): "
                f"{sorted(exp_pos - canon57)}\n")

    # ---- figure ----
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.3))

    ax = axes[0]
    colors = {"Predicted TSS": "#3B6FB6",
              "Experimental TSS (Jeong2016)": "#D1495B"}
    for label, d in sets.items():
        ax.plot(d["fpr"], d["tpr"], color=colors[label], lw=2.0,
                label=f"{label}\nAUC={d['auc']:.3f} [{d['lo']:.2f}, {d['hi']:.2f}]\n"
                      f"(n={d['n_exp']}/{d['n_shi']})")
    ax.plot([0, 1], [0, 1], "k:", lw=0.6, alpha=0.5)
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title("Classifier is TSS-definition robust", fontsize=10, fontweight="bold")
    ax.legend(fontsize=6.5, loc="lower right", frameon=True, edgecolor="#ccc")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_aspect("equal")

    ax = axes[1]
    bars = ax.bar([0, 1], [n_cov, n_ret], color=["#9CC3E6", "#D1495B"],
                  width=0.6, edgecolor="white")
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"Canonical Exposed\ntestable\n(n={n_cov})",
                        f"Retained Exposed\nunder exp. TSS\n(n={n_ret})"], fontsize=8)
    for b, v in zip(bars, [n_cov, n_ret]):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.15, str(v),
                ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.set_ylabel("Genes")
    ax.set_ylim(0, max(n_cov, n_ret) * 1.35 + 1)
    ax.set_title(f"{n_ret}/{n_cov} canonical members retained",
                 fontsize=10, fontweight="bold")
    ax.text(0.5, 0.92,
            f"0 reclassified to Shielded\n({len(canon57) - n_cov} of 57 untestable:\nno Jeong2016 dRNA-seq TSS)",
            transform=ax.transAxes, ha="center", va="top", fontsize=7,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#FFF8E1", edgecolor="#E0C200"))

    fig.tight_layout()
    for ext in ("pdf", "svg", "png"):
        fig.savefig(FIG / f"M4_experimental_TSS_robustness.{ext}",
                    dpi=200, bbox_inches="tight")
    plt.close(fig)

    print("=== M4 experimental-TSS robustness ===")
    print(tbl.to_string(index=False))
    print(f"\ncanonical Exposed testable under experimental TSS: {n_cov}/{len(canon57)}")
    print(f"retained Exposed: {n_ret}/{n_cov}  (reclassified to Shielded: {n_cov - n_ret})")
    print(f"experimental-only Exposed: {sorted(exp_pos - canon57)}")
    print(f"\nSaved: {FIG / 'M4_experimental_TSS_robustness.pdf'}")
    print(f"Saved: {TBL / 'M4_experimental_TSS_robustness.tsv'}")


if __name__ == "__main__":
    main()
