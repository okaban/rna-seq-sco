#!/usr/bin/env python3
"""Real figures for the P1 reviewer-robustness re-analysis bundle.
Reads the TSVs written by the P1_* scripts; muted Okabe-Ito palette."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.expanduser("~/bioinfo/rna-seq/11_epigenome_integration/analysis/78_reviewer_robustness")
T = os.path.join(BASE, "tables"); F = os.path.join(BASE, "figures")
OK = dict(blue="#0072B2", vermillion="#D55E00", green="#009E73",
          grey="#999999", purple="#CC79A7", orange="#E69F00")
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})


# ---- Fig A: cross-talk specificity + probability QC ------------------------
# 2026-09-20: read the C4 tables (4mC is at C4; the C3/C5 QC panel measured
# positions that carry ~0 modification). The near/far test itself is unchanged
# because it scans offsets around each 6mA call, not fixed motif positions.
rates = pd.read_csv(os.path.join(T, "P1_1_crosstalk_rates_C4.tsv"), sep="\t")
probs = pd.read_csv(os.path.join(T, "P1_1_prob_distributions_C4.tsv"), sep="\t")
fig, ax = plt.subplots(1, 3, figsize=(10, 3.2))
aag = rates[rates.group == "AAGCCCG_6mA"].iloc[0]
non = rates[rates.group == "nonAAGCCCG_6mA_generic"].iloc[0]
# panel a: near vs far co-call rate
x = np.arange(2); w = 0.35
ax[0].bar(x - w/2, [aag.near_rate, non.near_rate], w, label="near (±3,4 bp)", color=OK["vermillion"])
ax[0].bar(x + w/2, [aag.far_rate, non.far_rate], w, label="far (20–50 bp)", color=OK["grey"])
ax[0].set_yscale("log"); ax[0].set_xticks(x)
ax[0].set_xticklabels(["AAGCCCG\n6mA", "non-AAGCCCG\n6mA (generic)"])
ax[0].set_ylabel("4mC co-call rate / offset (log)")
ax[0].set_title("4mC near a 6mA call", loc="center", fontsize=9)
ax[0].legend(fontsize=7, frameon=False, loc="upper right")
# fold annotations sit just above each bar-pair. Headroom is opened so the
# title band above stays clear (R2-11: title/letter must not touch the label).
ax[0].set_ylim(top=ax[0].get_ylim()[1]*3)
ax[0].text(0 - w/2, aag.near_rate*1.35, f"{aag.near_over_far:.0f}×", ha="center", fontsize=8)
ax[0].text(1 - w/2, non.near_rate*1.35, f"{non.near_over_far:.2f}×", ha="center", fontsize=8)
# panel b: motif specificity
spec = aag.near_rate / non.near_rate
ax[1].bar([0, 1], [non.near_rate, aag.near_rate], color=[OK["grey"], OK["vermillion"]])
ax[1].set_yscale("log"); ax[1].set_xticks([0, 1])
ax[1].set_xticklabels(["generic\n6mA", "AAGCCCG\n6mA"])
ax[1].set_ylabel("near 4mC co-call rate (log)")
ax[1].set_title(f"motif specificity = {spec:.0f}×", loc="center", fontsize=9)
# panel c: probability QC
ax[2].plot(probs.prob_bin_center, probs.AAGCCCG_A_6mA/probs.AAGCCCG_A_6mA.sum(),
           "-o", ms=3, color=OK["blue"], label="A₀/A₁ (6mA)")
ax[2].plot(probs.prob_bin_center, probs.AAGCCCG_C_4mC/probs.AAGCCCG_C_4mC.sum(),
           "-s", ms=3, color=OK["vermillion"], label="C₄ (4mC)")
ax[2].set_xlabel("per-read modification probability")
ax[2].set_ylabel("fraction of calls")
ax[2].set_title("call confidence at AAGCCCG", loc="center", fontsize=9)
ax[2].legend(fontsize=7, frameon=False)
# uniform panel letters in the outer margin (R2-11: letter must not collide
# with the in-plot fold annotation). qc helper anchors a/b/c identically.
try:
    for _ax, _L in zip(ax, "abc"):
        qc_add_panel_label(_ax, _L)
except NameError:
    for _ax, _L in zip(ax, "abc"):
        _ax.text(-0.12, 1.06, _L, transform=_ax.transAxes,
                 fontsize=13, fontweight="bold", va="bottom", ha="right")
fig.tight_layout(); fig.savefig(os.path.join(F, "P1_1_crosstalk.png"), dpi=200); plt.close(fig)

# ---- Fig B: 5mC vs 4mC (P1-2) ---------------------------------------------
fig, ax = plt.subplots(figsize=(3.4, 3.2))
ax.bar([0, 1], [0.0014, 0.0226], color=[OK["grey"], OK["vermillion"]])
ax.axhline(0.01, ls="--", color="k", lw=0.8)
ax.text(1.4, 0.0105, "0.01% claim", fontsize=7, va="bottom", ha="right")
ax.set_xticks([0, 1]); ax.set_xticklabels(["5mC (m)", "4mC (21839)"])
ax.set_ylabel("genome-wide modified fraction (%)\ncytosine calls, cov≥10")
ax.set_title("Cytosine modification is m4C, not 5mC", loc="left", fontweight="bold", fontsize=9)
for i, v in enumerate([0.0014, 0.0226]):
    ax.text(i, v + 0.0008, f"{v:.4f}%", ha="center", fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(F, "P1_2_5mC_vs_4mC.png"), dpi=200); plt.close(fig)

# ---- Fig C: spatial autocorrelation (P1-3) --------------------------------
# rebuild a quick null for the figure (cheap)
from scipy import stats
A = os.path.expanduser("~/bioinfo/rna-seq/11_epigenome_integration/analysis")
dose = pd.read_csv(os.path.join(A, "41_GCCGGC_dose_response/tables/GCCGGC_gene_site_counts.tsv"), sep="\t")
trans = pd.read_csv(os.path.join(A, "42_GCCGGC_temporal_derepression/tables/gene_methylation_transitions.tsv"), sep="\t")[["gene_id", "LFC_T2vsT1"]]
d = dose.merge(trans, left_on="locus_tag", right_on="gene_id").dropna(subset=["LFC_T2vsT1", "gene_start", "dose_group"]).sort_values("gene_start").reset_index(drop=True)
code = d.dose_group.astype(str).map({"0":0,"1":1,"2":2,"3":3,"4+":4}).values
lfc = d.LFC_T2vsT1.values
obs, _ = stats.spearmanr(code, lfc)
rng = np.random.default_rng(42)
null = np.array([stats.spearmanr(code, np.roll(lfc, rng.integers(1, len(lfc))))[0] for _ in range(3000)])
fig, ax = plt.subplots(1, 2, figsize=(7.6, 3.2))
ax[0].hist(null, bins=40, color=OK["grey"], alpha=0.8)
ax[0].axvline(obs, color=OK["vermillion"], lw=2, label=f"observed ρ = {obs:+.3f}\n(block-perm p = 0.22)")
ax[0].axvline(-np.percentile(np.abs(null),95), color="k", ls=":", lw=0.8)
ax[0].axvline(np.percentile(np.abs(null),95), color="k", ls=":", lw=0.8, label="null 95%")
ax[0].set_xlabel("Spearman ρ (dose vs LFC)"); ax[0].set_ylabel("permutations")
# 2026-09-13 (FIG-05): neutral titles, no editorial labels ("→ artefact",
# "SURVIVES"); panel letters a-d run across the two stacked rows (S18).
ax[0].set_title("a  Genome-wide dose vs log$_2$FC (n = 7,496)\n    spatial-block permutation null",
                loc="left", fontweight="bold", fontsize=9)
ax[0].set_ylim(top=ax[0].get_ylim()[1] * 1.28)  # headroom for the legend
ax[0].legend(fontsize=7, frameon=False, loc="upper center")
# panel b: regulatory survives
reg = pd.read_csv(os.path.join(T, "P1_3b_permissive_spatial.tsv"), sep="\t").iloc[0]
ax[1].bar([0, 1], [reg.p_partial_naive, reg.block_perm_p_partial], color=[OK["grey"], OK["green"]])
ax[1].axhline(0.05, ls="--", color="k", lw=0.8); ax[1].set_yscale("log")
ax[1].set_xticks([0, 1]); ax[1].set_xticklabels(["naive p", "block-perm p"])
ax[1].set_ylabel("p (partial Spearman | region)")
ax[1].set_title(f"b  Regulatory genes (n = {int(reg.n):,})\n    partial Spearman ρ = {reg.rho_partial:+.2f}\n    (TSS–GCCGGC distance vs log$_2$FC | region)",
                loc="left", fontweight="bold", fontsize=8)
ax[1].text(0.5, 0.95, f"naive p = {reg.p_partial_naive:.1e}\nblock-perm p = {reg.block_perm_p_partial:.4f}",
           transform=ax[1].transAxes, fontsize=7.5, ha="center", va="top")
fig.tight_layout(); fig.savefig(os.path.join(F, "P1_3_spatial_autocorr.png"), dpi=200); plt.close(fig)

# ---- Fig D: equivalence + threshold sensitivity (P1-4/5) ------------------
eq = pd.read_csv(os.path.join(T, "P1_4_equivalence.tsv"), sep="\t")
sens = pd.read_csv(os.path.join(T, "P1_5_threshold_sensitivity.tsv"), sep="\t")
fig, ax = plt.subplots(1, 2, figsize=(7.6, 3.2))
row = eq[eq.contrast == "|LFC| T2vsT1"].iloc[0]
ax[0].axvspan(-0.147, 0.147, color=OK["green"], alpha=0.15, label="negligible band")
ax[0].errorbar([row.cliff_delta], [0], xerr=[[row.cliff_delta-row.cliff_lo],[row.cliff_hi-row.cliff_delta]],
               fmt="o", color=OK["vermillion"], capsize=4)
ax[0].axvline(0, color="k", lw=0.6)
ax[0].set_yticks([]); ax[0].set_xlim(-0.5, 0.5)
ax[0].set_xlabel("Cliff's δ (Exposed − Shielded |LFC|)")
ax[0].set_title(f"c  Permissive equivalence (TOST)\n    Cliff's δ = {row.cliff_delta:+.2f}, 95% CI [{row.cliff_lo:+.2f}, {row.cliff_hi:+.2f}]",
                loc="left", fontweight="bold", fontsize=9)
ax[0].legend(fontsize=7, frameon=False, loc="upper right")
ax[1].plot(sens.window_bp, sens.cliff_delta, "-o", color=OK["blue"], ms=4)
ax[1].axhspan(-0.147, 0.147, color=OK["green"], alpha=0.15)
ax[1].axvline(293, ls="--", color=OK["vermillion"], lw=1, label="293 bp")
ax[1].set_xlabel("Exposed window (bp)"); ax[1].set_ylabel("Cliff's δ (|LFC|)")
ax[1].set_ylim(-0.4, 0.4)
ax[1].set_title("d  Cliff's δ across Exposed-window\n    definitions (100–500 bp)", loc="left", fontweight="bold", fontsize=9)
ax[1].legend(fontsize=7, frameon=False)
fig.tight_layout(); fig.savefig(os.path.join(F, "P1_4_5_equivalence.png"), dpi=200); plt.close(fig)

# ---- S18 slot: rows C (a,b) and D (c,d) stacked into one PNG -------------
from PIL import Image
_top = Image.open(os.path.join(F, "P1_3_spatial_autocorr.png"))
_bot = Image.open(os.path.join(F, "P1_4_5_equivalence.png"))
_W = max(_top.width, _bot.width); _gap = 40
_st = Image.new("RGB", (_W, _top.height + _bot.height + _gap), "white")
_st.paste(_top, (0, 0)); _st.paste(_bot, (0, _top.height + _gap))
_st.save(os.path.join(F, "P1_3_4_5_stacked_S18.png"), dpi=(200, 200))

print("Figures ->", F)
for fn in ["P1_1_crosstalk.png","P1_2_5mC_vs_4mC.png","P1_3_spatial_autocorr.png","P1_4_5_equivalence.png"]:
    p = os.path.join(F, fn); print("  ", fn, "OK" if os.path.exists(p) else "MISSING")
