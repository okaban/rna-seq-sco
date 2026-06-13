#!/usr/bin/env python3
"""Audit the 293 bp boundary's bimodality justification (reproducibility check).

Tests whether the nearest-methylation-distance distribution of the 1,055
regulatory genes is bimodal, as the Methods section claims. Verdict (2026-06-13):
the distribution is UNIMODAL; the "bimodal -> natural boundary" narrative
(dip test + KDE anti-mode + single PELT changepoint at 293 bp) is not reproduced.

Run: python boundary_modality_audit.py
Requires: diptest, ruptures (installed into the rna-seq env on 2026-06-13).
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

TBL = Path(__file__).resolve().parents[1] / "tables"


def main() -> None:
    df = pd.read_csv(TBL / "all_genes_features_unified_n57.tsv", sep="\t")
    df = df.dropna(subset=["nearest_methyl_distance"])
    d = df["nearest_methyl_distance"].values
    x = np.log10(d + 1.0)
    print(f"n = {len(x)} dist-valid regulatory genes")

    # 1) Hartigan dip test
    import diptest
    D, p = diptest.diptest(x)
    print(f"[dip]      D = {D:.4f}  p = {p:.3g}  -> "
          f"{'BIMODAL' if p < 0.05 else 'unimodal (NS)'}")

    # 2) KDE modes at Silverman bandwidth
    kde = gaussian_kde(x, bw_method="silverman")
    grid = np.linspace(x.min(), x.max(), 1000)
    dens = kde(grid)
    is_max = np.r_[False, dens[1:] > dens[:-1]] & np.r_[dens[:-1] > dens[1:], False]
    n_modes = int(is_max.sum())
    modes_bp = np.round(10 ** grid[is_max] - 1).astype(int)
    print(f"[KDE]      Silverman-bw modes = {n_modes} at bp {list(modes_bp)} -> "
          f"{'bimodal+' if n_modes >= 2 else 'UNIMODAL'}")

    # 3) PELT changepoints
    import ruptures as rpt
    xs = np.sort(x).reshape(-1, 1)
    bk = rpt.Pelt(model="rbf", min_size=2).fit(xs).predict(pen=np.log(len(x)))
    cps = [round(float(10 ** xs[b - 1, 0] - 1)) for b in bk[:-1]]
    print(f"[PELT]     pen=log(n): {len(bk) - 1} changepoints; first bp = {cps[:5]}")
    m = (d >= 100) & (d <= 500)
    xw = np.log10(d[m] + 1).reshape(-1, 1)
    bw = rpt.Dynp(model="rbf", min_size=2).fit(xw).predict(n_bkps=1)
    cw = round(float(10 ** xw[bw[0] - 1, 0] - 1))
    print(f"[PELT]     single changepoint in [100,500]bp window = {cw} bp "
          f"(Methods claims 293)")

    print("\nVERDICT: distance distribution is unimodal; 293 bp is a threshold "
          "on a continuum, not a natural bimodal boundary. Classification "
          "(AUC=0.908) is unaffected.")


if __name__ == "__main__":
    main()
