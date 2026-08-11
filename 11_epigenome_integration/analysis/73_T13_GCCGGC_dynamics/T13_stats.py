#!/usr/bin/env python3
"""
T13 GCCGGC m4C temporal dynamics: formal statistical tests across T1/T2/T3.

Input
-----
GCCGGC_sites_by_timepoint.tsv (from analysis/37_defense_island_GCCGGC):
    columns include: chrom, position, strand, timepoint (T1|T2|T3),
                     frequency (per-site methylation fraction, 0-100),
                     final_motif, region, dist_to_DI, dist_to_RS19770.
    Sites are pre-filtered by the upstream modkit pileup pipeline
    (mod_code 'm', motif GCCGGC ±3bp context).

Tests
-----
1. Kruskal-Wallis H across all three timepoints (T1, T2, T3).
2. Pairwise Mann-Whitney U (two-sided) for T1-T2, T1-T3, T2-T3
   with rank-biserial correlation r as effect size.
3. Bonferroni correction across the three pairwise comparisons.

Outputs
-------
T13_results.tsv : tidy table of all test statistics and adjusted p-values.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

BASE = Path("/Users/okaban/bioinfo/rna-seq")
SITES = (
    BASE
    / "11_epigenome_integration/analysis/37_defense_island_GCCGGC/tables/"
    "GCCGGC_sites_by_timepoint.tsv"
)
OUT_DIR = BASE / "11_epigenome_integration/analysis/73_T13_GCCGGC_dynamics"
OUT_TSV = OUT_DIR / "T13_results.tsv"


def rank_biserial(u_stat: float, n1: int, n2: int) -> float:
    """Rank-biserial correlation derived from Mann-Whitney U.

    r = 1 - 2U / (n1 * n2). Range [-1, 1]; positive r means group1 > group2.
    """
    return 1.0 - (2.0 * u_stat) / (n1 * n2)


def main() -> None:
    df = pd.read_csv(SITES, sep="\t")
    needed = {"timepoint", "frequency"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"Input TSV missing columns: {missing}")

    df = df.dropna(subset=["frequency"]).copy()
    groups = {tp: df.loc[df["timepoint"] == tp, "frequency"].to_numpy() for tp in ("T1", "T2", "T3")}
    n = {tp: int(len(v)) for tp, v in groups.items()}
    medians = {tp: float(np.median(v)) if len(v) else float("nan") for tp, v in groups.items()}
    means = {tp: float(np.mean(v)) if len(v) else float("nan") for tp, v in groups.items()}

    print("=== Per-timepoint sample summary ===")
    for tp in ("T1", "T2", "T3"):
        v = groups[tp]
        print(
            f"  {tp}: n={n[tp]}, median={medians[tp]:.2f}, mean={means[tp]:.2f}, "
            f"min={float(np.min(v)):.2f}, max={float(np.max(v)):.2f}"
        )

    # --- Kruskal-Wallis ----------------------------------------------------
    h_stat, kw_p = stats.kruskal(groups["T1"], groups["T2"], groups["T3"])
    print("\n=== Kruskal-Wallis (T1 vs T2 vs T3) ===")
    print(f"  H = {h_stat:.4f}, p = {kw_p:.3e}, df = 2")

    # --- Pairwise Mann-Whitney U ------------------------------------------
    pairs = [("T1", "T2"), ("T1", "T3"), ("T2", "T3")]
    pairwise = []
    for a, b in pairs:
        x, y = groups[a], groups[b]
        u, p = stats.mannwhitneyu(x, y, alternative="two-sided")
        r = rank_biserial(u, len(x), len(y))
        pairwise.append({"a": a, "b": b, "U": float(u), "raw_p": float(p), "rb_r": float(r)})

    raw_ps = np.array([row["raw_p"] for row in pairwise])
    bonf_ps = np.minimum(raw_ps * len(pairs), 1.0)
    for row, padj in zip(pairwise, bonf_ps):
        row["bonferroni_p"] = float(padj)

    print("\n=== Pairwise Mann-Whitney U (two-sided, Bonferroni x3) ===")
    for row in pairwise:
        print(
            f"  {row['a']} vs {row['b']}: "
            f"U={row['U']:.1f}, raw_p={row['raw_p']:.3e}, "
            f"adj_p={row['bonferroni_p']:.3e}, rank-biserial r={row['rb_r']:+.3f}"
        )

    # --- Tidy results table ------------------------------------------------
    rows = [
        {
            "test": "Kruskal-Wallis",
            "comparison": "T1 vs T2 vs T3",
            "n_T1": n["T1"],
            "n_T2": n["T2"],
            "n_T3": n["T3"],
            "median_T1": medians["T1"],
            "median_T2": medians["T2"],
            "median_T3": medians["T3"],
            "statistic": h_stat,
            "raw_p": kw_p,
            "bonferroni_p": kw_p,  # KW is the single omnibus test, not corrected
            "effect_size": np.nan,
            "effect_size_label": "",
        }
    ]
    for row in pairwise:
        a, b = row["a"], row["b"]
        rows.append(
            {
                "test": "Mann-Whitney U (two-sided)",
                "comparison": f"{a} vs {b}",
                "n_T1": n["T1"] if a == "T1" or b == "T1" else np.nan,
                "n_T2": n["T2"] if a == "T2" or b == "T2" else np.nan,
                "n_T3": n["T3"] if a == "T3" or b == "T3" else np.nan,
                "median_T1": medians["T1"] if a == "T1" or b == "T1" else np.nan,
                "median_T2": medians["T2"] if a == "T2" or b == "T2" else np.nan,
                "median_T3": medians["T3"] if a == "T3" or b == "T3" else np.nan,
                "statistic": row["U"],
                "raw_p": row["raw_p"],
                "bonferroni_p": row["bonferroni_p"],
                "effect_size": row["rb_r"],
                "effect_size_label": "rank-biserial r",
            }
        )

    out = pd.DataFrame(rows)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_TSV, sep="\t", index=False, float_format="%.6g")
    print(f"\nWrote {OUT_TSV}")


if __name__ == "__main__":
    main()
