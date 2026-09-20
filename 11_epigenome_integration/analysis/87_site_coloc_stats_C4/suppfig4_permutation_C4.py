#!/usr/bin/env python3
"""Supplementary Figure 4 — permutation null for AAGCCCG dual modification, corrected offsets.

Replaces 68_or_permutation/A3_OR_CI_permutation.py, which plotted log10(OR) with an
observed OR of 138,440 computed through the 1-bp-misaligned sequence window.

With the corrected C4 census no co-modified pair occurs outside the motif (c = 0), so
the observed odds ratio is saturated and a log10(OR) axis has no finite observed value
to mark. The informative quantity is the one the permutation actually resamples: the
number of in-motif co-modified candidate pairs. Null = hypergeometric label permutation
over the 1,696,280 candidate pairs (same design as 87_/68_), 10,000 draws, seed 42.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from pathlib import Path
import sys; sys.path.insert(0, str(Path.home()/".claude-science"))
H = Path(__file__).resolve().parent
t = pd.read_csv(H/"tables/site_coloc_2x2_C4.tsv", sep="\t").set_index("scope")
FIG = "/Users/okaban/obsidian/Research/rna-seq/Writing/fig_images/SuppFigure4.png"

rng = np.random.default_rng(42); N = 10_000
panels = []
for scope, label in [("pooled", "All timepoints pooled"), ("T1", "T1 only")]:
    r = t.loc[scope]
    a, b, c, d = int(r.in_comod), int(r.in_notcomod), int(r.out_comod), int(r.out_notcomod)
    n_tot, n_in, n_comod = a+b+c+d, a+b, a+c
    null = rng.hypergeometric(n_comod, n_tot-n_comod, n_in, N)
    panels.append((label, null, a, n_comod*n_in/n_tot))

fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)
for ax, (label, null, obs, exp) in zip(axes, panels):
    bins = np.arange(-0.5, max(null.max(), 6)+1.5)
    ax.hist(null, bins=bins, color="#6B7783", edgecolor="white", linewidth=0.4)
    ax.axvline(exp, color="#3C3C3C", linestyle=":", linewidth=1.2)
    ax.annotate(f"observed\n{obs:,}", xy=(0.97, 0.62), xycoords="axes fraction",
                ha="right", va="center", color="#A64B44", fontsize=7)
    ax.annotate("", xy=(0.995, 0.40), xytext=(0.80, 0.40), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color="#A64B44", linewidth=1.6))
    ax.set_title(f"{label} (expected {exp:.2f}, max of null {null.max()})", fontsize=8)
    ax.set_xlabel("Co-modified pairs inside AAGCCCG\n(permutation null)", fontsize=8)
    ax.tick_params(labelsize=7, direction="in", which="both")
    ax.margins(0.04)
axes[0].set_ylabel("Permutation replicates", fontsize=8)
fig.suptitle("No permutation replicate reaches 2 % of the observed count (null maxima 5 of 406; 4 of 214)", fontsize=9, y=1.02)
fig.tight_layout()
fig.savefig(FIG, dpi=300, bbox_inches="tight")
r = fig.canvas.get_renderer()
texts = [(x, x.get_window_extent(r)) for x in fig.findobj(mpl.text.Text) if x.get_text().strip() and x.get_visible()]
print("overlaps:", [(a.get_text()[:18], b.get_text()[:18]) for i,(a,ba) in enumerate(texts) for b,bb in texts[i+1:] if ba.overlaps(bb)])
print("saved", FIG, [(l, int(o), round(e,2), int(n.max())) for l,n,o,e in panels])
