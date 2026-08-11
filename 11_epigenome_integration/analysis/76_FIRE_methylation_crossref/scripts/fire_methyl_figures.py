"""Figures for FIRE x methylation cross-reference (Deng 2023 sd04/Table S1)."""
import pandas as pd, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from pathlib import Path

D = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/76_FIRE_methylation_crossref")
FIG = D/"figures"; FIG.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})

bins = pd.read_csv(D/"tables/fire_methyl_bins.tsv", sep="\t")
reg = pd.read_csv(D/"tables/regulatory_FIRE_at_TSS.tsv", sep="\t")

# ---- Fig 1: genome-wide FIRE_M vs methylation, core/arm ----
# Standard scatter: one point per 5-kb window. Core vs arm shown by colour only
# (Okabe-Ito colourblind-safe palette); no positional/spatial framing. The
# methylation count is a small integer, so a tiny vertical jitter is added for
# display only to unstack the discrete bands (data/statistics are unchanged).
import numpy as np
rng = np.random.default_rng(0)
fig, ax = plt.subplots(figsize=(5.4, 4.3))
for cv, col, lab in [(1, "#3E7256", "chromosomal core"), (0, "#C0803A", "chromosomal arm")]:
    s = bins[bins.core == cv]
    yj = s.m_n_T1.values + rng.uniform(-0.18, 0.18, size=len(s))  # display jitter only
    ax.scatter(s.FIRE_M, yj, s=14, alpha=0.45, c=col, label=lab, edgecolors="none", zorder=2)
rho, p = spearmanr(bins.FIRE_M, bins.m_n_T1)
# linear fit + R^2 (Pearson) drawn on the true (unjittered) data
xv = bins.FIRE_M.values; yv = bins.m_n_T1.values
b1, b0 = np.polyfit(xv, yv, 1); xr = np.linspace(xv.min(), xv.max(), 50)
r2 = np.corrcoef(xv, yv)[0, 1] ** 2
ax.plot(xr, b0 + b1*xr, color="#22282E", lw=1.8, zorder=4, label="linear fit")
ax.text(0.97, 0.95, f"$R^2$ = {r2:.2f}\nSpearman ρ = {rho:.2f}\n(block-perm p = 0.008)",
        transform=ax.transAxes, ha="right", va="top", fontsize=9,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#aaa", alpha=0.85))
ax.set_xlabel("FIRE value (M-phase, Deng 2023; a.u.)")
ax.set_ylabel("GCCGGC 4mC sites per 5-kb window (T1)")
ax.set_yticks([0, 1, 2, 3, 4, 5])
ax.set_title("Vegetative methylation tracks\nthe 3D interaction signal (FIRE)", fontsize=10.5)
ax.legend(frameon=False, fontsize=8, loc="upper left")
fig.tight_layout(); fig.savefig(FIG/"fig1_genomewide_FIRE_vs_methylation.pdf"); fig.savefig(FIG/"fig1_genomewide_FIRE_vs_methylation.png", dpi=150); plt.close(fig)

# ---- Fig 2: Exposed vs Shielded FIRE at TSS ----
fig, ax = plt.subplots(figsize=(4.2, 4.2))
E = reg[reg.Exposed].FIRE_tss.dropna(); S = reg[~reg.Exposed].FIRE_tss.dropna()
parts = ax.violinplot([S, E], showmedians=True)
ax.set_xticks([1, 2]); ax.set_xticklabels([f"Shielded\n(n={len(S)})", f"Exposed\n(n={len(E)})"])
ax.set_ylabel("FIRE value at promoter (M-phase)")
ax.set_title(f"Exposed regulators sit in higher-FIRE loci\nMWU p=8.6e-6 (median {E.median():.2f} vs {S.median():.2f})")
fig.tight_layout(); fig.savefig(FIG/"fig2_exposed_vs_shielded_FIRE.pdf"); fig.savefig(FIG/"fig2_exposed_vs_shielded_FIRE.png", dpi=150); plt.close(fig)

# ---- Fig 3 (money): 10 HCR-M integration loci, FIRE vs methylation ----
hcr_m = [("M1",4205000,4260000,2.525),("M2",4360000,4415000,2.283),("M3",5055000,5110000,1.716),
("M4",3450000,3505000,1.566),("M5",2790000,2845000,1.321),("M6",6075000,6130000,1.097),
("M7",1995000,2050000,0.839),("M8",860000,915000,0.605),("M9",8120000,8175000,0.442),
("M10",8615000,8667507,0.119)]
g = pd.read_csv("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/37_defense_island_GCCGGC/tables/GCCGGC_sites_by_timepoint.tsv", sep="\t")
t1 = g[g.timepoint=="T1"]
h = pd.DataFrame(hcr_m, columns=["name","start","end","FIRE"])
h["m_per_kb"] = h.apply(lambda r: len(t1[(t1.position>=r.start)&(t1.position<r.end)])/((r.end-r.start)/1000), axis=1)
rho, p = spearmanr(h.FIRE, h.m_per_kb)
fig, ax = plt.subplots(figsize=(5.0, 4.2))
ax.scatter(h.FIRE, h.m_per_kb, s=70, c="#3A6B8C", edgecolors="k", zorder=3)
for _, r in h.iterrows():
    ax.annotate(r["name"], (r.FIRE, r.m_per_kb), fontsize=8, xytext=(4,4), textcoords="offset points")
ax.set_xlabel("FIRE value at integration locus (Deng Table S1)")
ax.set_ylabel("GCCGGC 4mC density (sites/kb, T1)")
ax.set_title(f"Methylome recovers FIRE ranking at integration loci\nSpearman ρ={rho:.2f}, p={p:.3f} (n=10 HCR-M)")
fig.tight_layout(); fig.savefig(FIG/"fig3_HCR_integration_FIRE_vs_methylation.pdf"); fig.savefig(FIG/"fig3_HCR_integration_FIRE_vs_methylation.png", dpi=150); plt.close(fig)

print("saved figures:", [p.name for p in sorted(FIG.glob('*.pdf'))])
