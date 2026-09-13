"""C5/C6/C7 summary figure.

Panel A  Genome track: FIRE (5-kb) and T1 GCCGGC methylation density along the
         chromosome, hotspot windows marked, BGCs + Compartment A annotated.
Panel B  FIRE vs methylation scatter (all bins), hotspot quadrant highlighted.
Panel C  C7 resolution sweep: candidate-window count and FIRE concordance
         (precision / Jaccard) vs window size.

Reads tables produced by C5C6C7_fire_hotspots.py. Output PNG only.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from pathlib import Path

BASE = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis")
OUT = BASE / "78_reviewer_robustness"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)
BIN = 5000
COMP_A_LO, COMP_A_HI = 2_300_000, 6_200_000
BGCS = [("cda", 3_519_449, 3_602_320), ("act", 5_513_809, 5_535_091),
        ("red", 6_432_812, 6_464_206), ("cpk", 6_900_898, 6_948_414)]

# palette (manuscript-consistent muted tones)
C_FIRE = "#2c6e9c"
C_METH = "#c0504d"
C_HOT = "#e8a33d"
C_COMPA = "#dfe7ef"
C_BGC = "#6a6a6a"

bins = pd.read_csv(BASE / "76_FIRE_methylation_crossref/tables/fire_methyl_bins.tsv",
                   sep="\t")
bins = bins[bins["FIRE_M"].notna()].copy()
bins["center"] = bins["bin"] * BIN + BIN / 2

# recompute hotspot flags (same thresholds as analysis)
fire_thr = bins["FIRE_M"].quantile(0.80)
methyl_thr = bins.loc[bins["m_n_T1"] > 0, "m_n_T1"].quantile(0.80)
bins["hot"] = ((bins["FIRE_M"] >= fire_thr) & (bins["m_n_T1"] >= methyl_thr)).astype(int)

loci = pd.read_csv(OUT / "tables/C5_hotspot_loci_merged.tsv", sep="\t")
sweep = pd.read_csv(OUT / "tables/C7_resolution_sweep.tsv", sep="\t")

fig = plt.figure(figsize=(13, 9))
gs = fig.add_gridspec(3, 2, height_ratios=[1.15, 1.0, 1.15],
                      hspace=0.55, wspace=0.32)

# ---------------- Panel A: genome track ----------------
axA = fig.add_subplot(gs[0, :])
mb = bins["center"] / 1e6
# Compartment A shading
axA.axvspan(COMP_A_LO / 1e6, COMP_A_HI / 1e6, color=C_COMPA, zorder=0,
            label="Compartment A")
# FIRE line
axA.plot(mb, bins["FIRE_M"], color=C_FIRE, lw=0.7, alpha=0.85, label="FIRE (5 kb)")
axA.axhline(fire_thr, color=C_FIRE, ls=":", lw=0.8, alpha=0.6)
axA.set_ylabel("FIRE score", color=C_FIRE)
axA.tick_params(axis="y", labelcolor=C_FIRE)
axA.set_xlim(0, 8.7)
axA.set_ylim(bottom=min(0, bins["FIRE_M"].min()))
# methylation on twin axis
axM = axA.twinx()
axM.bar(mb, bins["m_n_T1"], width=BIN / 1e6, color=C_METH, alpha=0.45,
        label="T1 GCCGGC sites / 5 kb")
axM.set_ylabel("GCCGGC sites / 5 kb", color=C_METH)
axM.tick_params(axis="y", labelcolor=C_METH)
axM.set_ylim(0, bins["m_n_T1"].max() * 1.05)
# hotspot ticks
hb = bins[bins["hot"] == 1]
axA.scatter(hb["center"] / 1e6, [axA.get_ylim()[1] * 0.97] * len(hb),
            marker="v", s=18, color=C_HOT, zorder=5,
            label="hotspot bin (n=%d)" % len(hb))
# BGCs
for name, s, e in BGCS:
    axA.axvspan(s / 1e6, e / 1e6, color=C_BGC, alpha=0.85, lw=0)
    axA.text((s + e) / 2 / 1e6, axA.get_ylim()[1] * 0.60, name, ha="center",
             va="top", fontsize=7, color="#222222", fontweight="bold",
             bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
axA.set_xlabel("Chromosome position (Mb)")
axA.set_title("a  FIRE and T1 GCCGGC methylation along the chromosome; "
              "co-high hotspots concentrate in the central core (Compartment A; IQR 3.9–4.8 Mb)",
              fontsize=10, loc="left")
h1, l1 = axA.get_legend_handles_labels()
h2, l2 = axM.get_legend_handles_labels()
# 2026-09-13 (FIG-10): legend moved out of the data area (it covered the FIRE
# trace at 6-8 Mb and the hotspot triangles) to a single row under the x-axis.
axA.legend(h1 + h2, l1 + l2, fontsize=7, loc="upper center", ncol=4,
           bbox_to_anchor=(0.5, -0.22), frameon=False)

# ---------------- Panel B: FIRE vs methylation scatter ----------------
axB = fig.add_subplot(gs[1, 0])
bg = bins[bins["hot"] == 0]
axB.scatter(bg["m_n_T1"] + np.random.uniform(-0.15, 0.15, len(bg)),
            bg["FIRE_M"], s=6, color="#b9c2cb", alpha=0.4, label="other bins")
axB.scatter(hb["m_n_T1"] + np.random.uniform(-0.15, 0.15, len(hb)),
            hb["FIRE_M"], s=14, color=C_HOT, edgecolor="#9a6a14", lw=0.3,
            label="hotspot")
axB.axhline(fire_thr, color=C_FIRE, ls=":", lw=1)
axB.axvline(methyl_thr - 0.5, color=C_METH, ls=":", lw=1)
axB.set_xlabel("T1 GCCGGC sites / 5-kb bin")
axB.set_ylabel("FIRE score")
axB.set_title("b  Co-high quadrant (top-20% × top-20%)", fontsize=10, loc="left")
axB.legend(fontsize=7, loc="upper right")
axB.set_xlim(-0.5, bins["m_n_T1"].max() + 0.5)

# ---------------- Panel C: resolution sweep ----------------
axC = fig.add_subplot(gs[1, 1])
sw = sweep.dropna(subset=["jaccard"]).copy()
x = np.arange(len(sw))
axC.plot(x, sw["precision_methyl->FIRE"], "-o", color=C_FIRE, ms=5,
         label="precision (top-methyl→FIRE)")
axC.plot(x, sw["recall_FIRE_captured"], "-s", color=C_METH, ms=5,
         label="recall (FIRE captured)")
axC.plot(x, sw["jaccard"], "-^", color=C_HOT, ms=5, label="Jaccard")
axC.set_xticks(x)
axC.set_xticklabels(["%dkb" % (w // 1000) for w in sw["window_bp"]])
axC.set_xlabel("Window size")
axC.set_ylabel("Concordance")
axC.set_ylim(0, 0.65)
axC.set_title("c  Methylation↔FIRE concordance vs resolution", fontsize=10,
              loc="left")
axC.legend(fontsize=7, loc="upper left")
axC.grid(alpha=0.25)

# ---------------- Panel D: candidate window counts ----------------
axD = fig.add_subplot(gs[2, 0])
allsw = sweep.copy()
xx = np.arange(len(allsw))
axD.bar(xx - 0.2, allsw["n_methyl_windows"], width=0.4, color="#9bb7cd",
        label="methylated windows")
axD.bar(xx + 0.2, allsw["n_top_methyl_windows"], width=0.4, color=C_METH,
        label="top-methyl candidate windows")
axD.set_xticks(xx)
axD.set_xticklabels(["%dkb" % (w // 1000) for w in allsw["window_bp"]])
axD.set_xlabel("Window size")
axD.set_ylabel("Number of windows")
axD.set_title("d  Candidate-window count by resolution", fontsize=10, loc="left")
axD.set_ylim(0, allsw["n_methyl_windows"].max() * 1.30)  # headroom for count labels
axD.legend(fontsize=7, loc="upper right")
for xi, n in zip(xx + 0.2, allsw["n_top_methyl_windows"]):
    axD.text(xi, n + 8, str(int(n)), ha="center", fontsize=6.5, color=C_METH)

# ---------------- Panel E: top loci table-as-bars ----------------
axE = fig.add_subplot(gs[2, 1])
top = loci.head(10).iloc[::-1]
ypos = np.arange(len(top))
axE.barh(ypos, top["sum_methyl_T1"], color=C_HOT, edgecolor="#9a6a14", lw=0.4)
labels = ["%.2f–%.2f Mb (%dkb, %s+%dkb)" %
          (r.start_bp / 1e6, r.end_bp / 1e6, r.width_kb, r.nearest_BGC,
           r.dist_to_BGC_bp // 1000) for r in top.itertuples()]
axE.set_yticks(ypos)
axE.set_yticklabels(labels, fontsize=5.8)
axE.set_xlabel("Σ T1 GCCGGC sites in locus")
axE.set_title("e  Top-10 discrete hotspot loci", fontsize=10, loc="left")
for yi, (m, f) in enumerate(zip(top["sum_methyl_T1"], top["max_FIRE"])):
    axE.text(m + 0.1, yi, "FIRE %.1f" % f, va="center", fontsize=6)

# 2026-09-13 (FIG-10): suptitle removed — it carried internal reviewer-comment
# codes "(C5/C6/C7)"; the figure title belongs in the manuscript legend.
out = FIG / "C5C6C7_hotspots_summary.png"
fig.savefig(out, dpi=200, bbox_inches="tight")
print("[saved]", out)
