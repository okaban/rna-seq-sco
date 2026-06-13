#!/usr/bin/env python3
"""Effect sizes for Exposed vs Shielded nearest-methylation distance.

Cliff's delta (non-parametric) and Cohen's d on log10(distance+1), with
non-parametric bootstrap 95% CIs. Canonical 57/998 set.
"""

from pathlib import Path

import numpy as np
import pandas as pd

TBL = Path(__file__).resolve().parents[1] / "tables"
RNG = np.random.default_rng(42)


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """Delta = P(a>b) - P(a<b). a=Exposed (smaller distance) -> negative."""
    # memory-safe pairwise via sorting
    n_a, n_b = len(a), len(b)
    more = sum((a[:, None] > b[None, :]).sum() for _ in [0])
    less = (a[:, None] < b[None, :]).sum()
    return float((more - less) / (n_a * n_b))


def cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    na, nb = len(a), len(b)
    sp = np.sqrt(((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2))
    return float((a.mean() - b.mean()) / sp)


def boot_ci(a: np.ndarray, b: np.ndarray, fn, n: int = 10000) -> tuple[float, float]:
    vals = []
    for _ in range(n):
        sa = a[RNG.integers(0, len(a), len(a))]
        sb = b[RNG.integers(0, len(b), len(b))]
        vals.append(fn(sa, sb))
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def main() -> None:
    df = pd.read_csv(TBL / "all_genes_features_unified_n57.tsv", sep="\t")
    df = df.dropna(subset=["nearest_methyl_distance"])
    exp = df.loc[df.is_exposed == 1, "nearest_methyl_distance"].values.astype(float)
    shi = df.loc[df.is_exposed == 0, "nearest_methyl_distance"].values.astype(float)
    le = np.log10(exp + 1)
    ls = np.log10(shi + 1)

    delta = cliffs_delta(exp, shi)
    dlo, dhi = boot_ci(exp, shi, cliffs_delta, n=2000)
    d = cohens_d(le, ls)
    clo, chi = boot_ci(le, ls, cohens_d, n=10000)

    print(f"n_Exposed={len(exp)}  n_Shielded={len(shi)}")
    print(f"Cliff's delta = {delta:.3f}  95% CI [{dlo:.3f}, {dhi:.3f}]")
    print(f"Cohen's d (log10 dist) = {d:.3f}  95% CI [{clo:.3f}, {chi:.3f}]")


if __name__ == "__main__":
    main()
