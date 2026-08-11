#!/usr/bin/env python3
"""Supplementary figure: m4C calling specificity (GAP-3) and GCCGGC strand
structure / hemimethylation (GAP-5).  Numbers verified from
G3G5_background_and_symmetry.py on T1 replicate 1-3 (reproduced in 1-2)."""
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C_4MC = "#D55E00"     # Okabe-Ito vermillion (4mC)
C_GREY = "#999999"

# --- verified GAP-3 numbers (T1 rep 1-3) ---
BG_N = 6_140_891
GC_N, GC_MET = 31_679, 1_544          # GCCGGC modified C; methylated >=50%
# non-motif background positive rate (%) at freq thresholds 50/25/10
bg_rate = {50: 0.0, 25: 100*2/BG_N, 10: 100*102/BG_N}
gc_rate50 = 100*GC_MET/GC_N           # 4.87%

# --- GAP-5 active pairs ---
PAIRS = "tables/G3G5_T1rep13_active_pairs.tsv"
xp, yp = [], []
with open(PAIRS) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        xp.append(float(row["plus_frac"])); yp.append(float(row["minus_frac"]))
act = [(a, b) for a, b in zip(xp, yp) if max(a, b) >= 50]
hemi = sum(1 for a, b in act if (a >= 50) ^ (b >= 50))

fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.6))

# Panel A: specificity at canonical threshold
labels = ["non-motif C\n(n=6.14M)", "GCCGGC\nmodified C\n(n=31,679)"]
vals = [bg_rate[50], gc_rate50]
bars = axA.bar(labels, vals, color=[C_GREY, C_4MC], width=0.55)
axA.set_ylabel("cytosines called 4mC-positive (freq ≥ 50%), %")
axA.set_title("4mC calling specificity (empirical FDR proxy)", loc="center",
              fontsize=9)
axA.set_ylim(0, gc_rate50 * 1.25)
axA.bar_label(bars, labels=["0 / 6.14M\n(0.0000%)", f"{gc_rate50:.2f}%"],
              padding=3, fontsize=9)
# (relaxed-background counts moved to the figure caption per reviewer R2-12)
axA.spines[["top", "right"]].set_visible(False)

# Panel B: strand symmetry / hemimethylation
axB.scatter(xp, yp, s=6, alpha=0.25, color=C_4MC, edgecolors="none")
axB.axhline(50, ls=":", c=C_GREY, lw=0.8); axB.axvline(50, ls=":", c=C_GREY, lw=0.8)
axB.set_xlabel("+ strand 4mC frequency (%)")
axB.set_ylabel("− strand 4mC frequency (%)")
axB.set_title("GCCGGC palindrome strand structure", loc="center", fontsize=9)
axB.set_xlim(-3, 103); axB.set_ylim(-3, 103)
# minimal quadrant orientation labels only (statistics → caption, R2-12):
# points near an axis = one strand methylated (hemimethylated); points in the
# upper-right = both strands (fully methylated).
axB.text(52, 8, "hemimethylated", fontsize=8, color="#222")
axB.text(62, 62, "full", fontsize=8, color=C_GREY, ha="center")
axB.spines[["top", "right"]].set_visible(False)

# panel letters via qc helper (uniform, outer margin)
try:
    qc_add_panel_label(axA, "a"); qc_add_panel_label(axB, "b")
except NameError:
    axA.text(-0.14, 1.04, "a", transform=axA.transAxes, fontsize=13,
             fontweight="bold", va="bottom", ha="right")
    axB.text(-0.14, 1.04, "b", transform=axB.transAxes, fontsize=13,
             fontweight="bold", va="bottom", ha="right")
fig.tight_layout()
out = "figures/G3G5_specificity_symmetry.png"
fig.savefig(out, dpi=200)
print("[saved]", out)
