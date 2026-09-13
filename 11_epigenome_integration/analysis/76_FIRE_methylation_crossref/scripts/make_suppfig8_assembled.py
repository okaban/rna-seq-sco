#!/usr/bin/env python3
"""Supplementary Figure 8 (manuscript slot fig_images/SuppFigure8.png) — assembled
2026-09-13 (review finding FIG-04) to match the legend panel-for-panel:

  (a) genome-wide, per 5-kb window (n = 1,724): T1 GCCGGC 4mC site count vs
      M-phase FIRE, coloured by Deng compartment-A core (2.3–6.2 Mb) vs arm;
  (b) FIRE at the TSS of Exposed (n = 62) vs Shielded (n = 985 of 989) regulators;
  (c) developmental co-change: Δ methylation (T1→T2) vs Δ FIRE (M→L) per window.

The previously embedded PNG carried the withdrawn integration-locus ranking
panel (former S12) in slot (c); that panel is dropped here.

Inputs (unchanged): tables/fire_methyl_bins.tsv, tables/regulatory_FIRE_at_TSS.tsv
(the `Exposed` column there equals the canonical 62-gene identity table —
verified 2026-09-13, symmetric difference 0). Only Spearman ρ / n / medians are
printed on the panels; the spatial-block permutation p values quoted in the
legend are not recomputed here.
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, mannwhitneyu
from pathlib import Path

D = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/76_FIRE_methylation_crossref")
FIG = D / "figures"; FIG.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7.5,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6.5,
                     "axes.spines.top": False, "axes.spines.right": False})
C_CORE, C_ARM = "#3E7256", "#C0803A"          # Okabe-Ito-like, matches fig1
C_SHIELDED, C_EXPOSED = "#9AA7B0", "#A64B44"  # unified palette (Shielded grey / 4mC red)
FULL_W = 6.85  # NAR full width, 174 mm

bins = pd.read_csv(D / "tables/fire_methyl_bins.tsv", sep="\t").dropna(subset=["FIRE_M"])
reg = pd.read_csv(D / "tables/regulatory_FIRE_at_TSS.tsv", sep="\t")

fig, ax = plt.subplots(1, 3, figsize=(FULL_W, 2.6))
rng = np.random.default_rng(0)

# ---- (a) genome-wide scatter, core/arm ----
for cv, col, lab in [(0, C_ARM, "arm"), (1, C_CORE, "core (2.3–6.2 Mb)")]:
    s = bins[bins.core == cv]
    yj = s.m_n_T1.values + rng.uniform(-0.18, 0.18, size=len(s))  # display jitter only
    ax[0].scatter(s.FIRE_M, yj, s=6, alpha=0.45, c=col, label=lab, edgecolors="none")
rho_a, p_a = spearmanr(bins.FIRE_M, bins.m_n_T1)
ax[0].set_xlabel("FIRE (M-phase, Deng 2023; a.u.)")
ax[0].set_ylabel("GCCGGC 4mC sites per 5-kb window (T1)")
ax[0].set_yticks([0, 1, 2, 3, 4, 5])
ax[0].set_title(f"Genome-wide, per 5-kb window\nSpearman ρ = {rho_a:.2f}, n = {len(bins):,}", loc="left")
ax[0].set_ylim(-0.45, 6.6)  # headroom so the legend sits above the y=5 points
ax[0].legend(frameon=False, loc="upper right", markerscale=2, ncol=2, columnspacing=0.8)

# ---- (b) FIRE at TSS: Shielded vs Exposed ----
d = reg.dropna(subset=["FIRE_tss"])
S = d.loc[~d.Exposed, "FIRE_tss"].values; E = d.loc[d.Exposed, "FIRE_tss"].values
n_sh_total = int((~reg.Exposed).sum()); n_ex_total = int(reg.Exposed.sum())
parts = ax[1].violinplot([S, E], positions=[0, 1], showextrema=False, widths=0.8)
for body, col in zip(parts["bodies"], [C_SHIELDED, C_EXPOSED]):
    body.set_facecolor(col); body.set_edgecolor(col); body.set_alpha(0.6)
for x, v in [(0, S), (1, E)]:
    ax[1].hlines(np.median(v), x - 0.22, x + 0.22, color="k", lw=1.4, zorder=5)
    ax[1].text(x, np.median(v) + 0.06, f"median {np.median(v):.2f}", ha="center",
               va="bottom", fontsize=6.5)
U, p_b = mannwhitneyu(E, S, alternative="two-sided")
ax[1].set_xticks([0, 1])
ax[1].set_xticklabels([f"Shielded\n(n = {len(S)} of {n_sh_total})", f"Exposed\n(n = {len(E)} of {n_ex_total})"])
ax[1].set_ylabel("FIRE at promoter TSS (M-phase)")
ax[1].set_title(f"Promoter FIRE by class\nMann–Whitney p = {p_b:.1e}", loc="left")

# ---- (c) developmental co-change ----
dm = bins.m_n_T2 - bins.m_n_T1; dfire = bins.FIRE_L - bins.FIRE_M
rho_c, p_c = spearmanr(dm, dfire)
ax[2].scatter(dm + rng.uniform(-0.15, 0.15, size=len(dm)), dfire, s=6, alpha=0.4,
              c="#6B7783", edgecolors="none")
ax[2].axhline(0, color="k", lw=0.5); ax[2].axvline(0, color="k", lw=0.5)
ax[2].set_xlabel("Δ GCCGGC 4mC sites / 5 kb (T1→T2)")
ax[2].set_ylabel("Δ FIRE (M→L phase)")
ax[2].set_title(f"Developmental co-change\nSpearman ρ = {rho_c:.2f}, n = {len(bins):,}", loc="left")

fig.tight_layout(w_pad=2.2, rect=[0, 0, 0.99, 0.94])
# panel letters on one baseline above the tallest tight bbox
fig.canvas.draw(); rend = fig.canvas.get_renderer()
tbs = [a.get_tightbbox(rend).transformed(fig.transFigure.inverted()) for a in ax]
top = max(tb.y1 for tb in tbs)
for a, tb, lab in zip(ax, tbs, "abc"):
    fig.text(max(0.002, tb.x0 - 0.006), min(top + 0.004, 0.975), lab,
             fontsize=13, fontweight="bold", va="bottom", ha="left")

out = FIG / "SuppFig8_compartment_scale_assembled"
fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")
# clamp to NAR full-width pixel ceiling (bbox='tight' can re-expand)
from PIL import Image
im = Image.open(out.with_suffix(".png"))
if im.size[0] > 2055:
    im.resize((2055, round(im.size[1] * 2055 / im.size[0])), Image.LANCZOS).save(out.with_suffix(".png"), dpi=(300, 300))
print(f"(a) rho={rho_a:.4f} p={p_a:.2e} n={len(bins)} | (b) medS={np.median(S):.4f} medE={np.median(E):.4f} "
      f"nS={len(S)}/{n_sh_total} nE={len(E)}/{n_ex_total} MWU p={p_b:.2e} | (c) rho={rho_c:.4f} p={p_c:.2e} n={len(bins)}")
print("[saved]", out.with_suffix(".png"))
