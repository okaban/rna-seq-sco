"""Reviewer-requested figures (2026-06-16).
Q10: modification-position schematic (which base of each motif is modified).
Q14: Fig2 redistribution — linear chromosome landscape + quantitative core/arm,
     for BOTH GCCGGC and AAGCCCG (GCCGGC strongly dynamic, AAGCCCG as contrast).
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from pathlib import Path

B = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis")
OUT = B/"77_reviewer_figures/figures"; OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
GLEN = 8_667_507; BIN = 50_000
# Match Figure 1a exactly: core 1.5–7.17 Mb, arms outer 1.5 Mb, oriC at dnaA.
ARM_LEFT = 1_500_000; ARM_RIGHT = 7_167_507; ORIC = 4_270_777
COL_CORE = "#88CCEE"; COL_ARM = "#DDCC77"  # muted Tol, matches Fig1
CORE_LO, CORE_HI = ARM_LEFT, ARM_RIGHT
C4 = "#C26B6B"; C6 = "#4477AA"  # muted Tol (4mC rose / 6mA blue), matches Fig1

# ---------- data ----------
# 2026-09-21 (BLOCKER-0): both former sources (37_ GCCGGC table; 07_ master site
# table) are position-DEDUPLICATED across timepoints (first-appearance), which gave
# the retracted 83/18/38 % and 82/38/61 % core series. Read the canonical
# per-timepoint file through 90_/canonical_sites.py: GCCGGC 4mC 1,289/1,595/1,073
# (83.0/66.0/68.7 % core); AAGCCCG with 4mC(C4) and 6mA(A0/A1) pooled
# 1,116/1,302/1,017 (82.2/70.4/73.0 % core; motif called from the reference, not from
# the 1-bp-misaligned `sequence` window).
import importlib.util as _ilu
_spec = _ilu.spec_from_file_location("canonical_sites", B/"90_per_timepoint_census_audit/canonical_sites.py")
_cs = _ilu.module_from_spec(_spec); _spec.loader.exec_module(_cs)
g = _cs.gccggc_by_timepoint()
a = _cs.aagcccg_pooled_by_timepoint()[['position','timepoint']].copy()

# ============================================================
# FIGURE A — modification position within each motif (Q10)
#   left: schematic of the modified base; right: quantitative per-position
#         mean modification frequency (empirical evidence for which base).
# ============================================================
_COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}

def draw_motif(ax, seq, mods, title, bot_mods=None):
    """Double-stranded schematic.

    seq      : top strand 5'->3' (list of bases).
    mods     : dict idx -> (label, color) for TOP-strand modifications.
    bot_mods : optional dict idx -> (label, color) for BOTTOM-strand
               modifications, keyed by the SAME top-strand index (the box sits
               under that column). Use this to show the unmodified / hemi
               partner on the complementary strand.
    """
    bot_mods = bot_mods or {}
    n = len(seq)
    comp = [_COMP.get(b, b) for b in seq]  # complementary base under each column
    TOP_Y, BOT_Y, BH = 1.15, 0.0, 0.9       # row y-origins and box height

    def _box(i, y, base, mod):
        modded = mod is not None
        fc = mod[1] if modded else "#ecf0f1"
        ax.add_patch(FancyBboxPatch((i, y), 0.9, BH, boxstyle="round,pad=0.02",
                     fc=fc, ec="#2c3e50", lw=1.2))
        ax.text(i + 0.45, y + BH / 2, base, ha="center", va="center",
                fontsize=14, fontweight="bold",
                color="white" if modded else "#2c3e50")

    # top strand (5'->3') with arrows pointing DOWN onto the modified base
    for i, base in enumerate(seq):
        _box(i, TOP_Y, base, mods.get(i))
        if i in mods:
            ax.annotate(mods[i][0], (i + 0.45, TOP_Y + BH), ha="center", va="bottom",
                        fontsize=9, fontweight="bold", color=mods[i][1],
                        xytext=(i + 0.45, TOP_Y + BH + 0.75), textcoords="data",
                        arrowprops=dict(arrowstyle="->", color=mods[i][1], lw=1.4))
    # bottom strand (3'->5') with arrows pointing UP onto the modified base
    for i, base in enumerate(comp):
        _box(i, BOT_Y, base, bot_mods.get(i))
        if i in bot_mods:
            ax.annotate(bot_mods[i][0], (i + 0.45, BOT_Y), ha="center", va="top",
                        fontsize=9, fontweight="bold", color=bot_mods[i][1],
                        xytext=(i + 0.45, BOT_Y - 0.75), textcoords="data",
                        arrowprops=dict(arrowstyle="->", color=bot_mods[i][1], lw=1.4))
    # base-pair rungs between the two strands
    for i in range(n):
        ax.plot([i + 0.45, i + 0.45], [BOT_Y + BH, TOP_Y], color="#bbb", lw=0.8, zorder=0)
    # strand-end labels
    ax.text(-0.15, TOP_Y + BH / 2, "5′", ha="right", va="center", fontsize=10)
    ax.text(n + 0.0, TOP_Y + BH / 2, "3′", ha="left", va="center", fontsize=10)
    ax.text(-0.15, BOT_Y + BH / 2, "3′", ha="right", va="center", fontsize=10)
    ax.text(n + 0.0, BOT_Y + BH / 2, "5′", ha="left", va="center", fontsize=10)
    ax.set_xlim(-0.8, n + 0.7); ax.set_ylim(-1.05, TOP_Y + BH + 1.15)
    ax.axis("off"); ax.set_title(title, fontsize=11, fontweight="bold")

def quant_positions(ax, seq, freqs, title, ylab="mean modification\nfrequency (%)"):
    """freqs: dict idx->(value,color,label). Bar per motif position."""
    n=len(seq)
    for i in range(n):
        if i in freqs:
            v,col,_=freqs[i]; ax.bar(i, v, color=col, width=0.8, zorder=3)
            ax.text(i, v+1.5, f"{v:.0f}%", ha="center", fontsize=8, color=col, fontweight="bold")
        else:
            ax.bar(i, 0, width=0.8)
    ax.set_xticks(range(n)); ax.set_xticklabels(list(seq), fontsize=11, fontweight="bold")
    ax.set_ylim(0, 118); ax.set_ylabel(ylab, fontsize=8.5)
    ax.set_title(title, fontsize=9.5)

import matplotlib.patches as mpatches
fig, axes = plt.subplots(2, 3, figsize=(15.5, 7.2), gridspec_kw={"width_ratios":[1,1.15,0.85]})
# GCCGGC (palindrome): 4mC at internal C (C2) on the TOP strand only.
# Hemi-methylated: the palindromic partner C on the BOTTOM strand (under
# column 3, where the complementary base is C) is ~0% methylated.
GREY0 = "#9aa0a6"  # marks the unmodified hemi partner
draw_motif(axes[0,0], list("GCCGGC"), {2:("4mC", C4)},
           "GCCGGC  (palindromic; hemi-methylated 4mC)",
           bot_mods={3:("", GREY0)})
quant_positions(axes[0,1], "GCCGGC", {2:(82,C4,"4mC")},
                "Per-position modification frequency — GCCGGC")
# (per-position occupancy detail moved to the figure caption)
# AAGCCCG: 6mA at A0/A1 (~23%), 4mC at C3/C5 (~35%); dominant co-modified
# pair A1<->C5. Both modification types sit on the SAME (top) strand; the
# bottom strand carries no modification.
# 2026-09-20: the 4mC sits at C4, not C3/C5. The earlier offsets came from the
# 1-bp-misaligned `sequence` window of methylation_site_sequences.csv; the T1
# pileups put 4mC at C4 (mean 74.8% at canonical sites, 65.0% of instances
# called) and ~0% at C3/C5 (0.15% / 0.05%). 6mA: A0 62.1%, A1 59.6% at canonical
# sites. Values from 89_supp_deliverables/per_position_motif_frequency_T1.tsv.
draw_motif(axes[1,0], list("AAGCCCG"), {0:("6mA", C6), 1:("6mA", C6), 4:("4mC", C4)},
           "AAGCCCG  (same-strand dual modification)")
quant_positions(axes[1,1], "AAGCCCG", {0:(62,C6,"6mA"),1:(60,C6,"6mA"),4:(75,C4,"4mC")},
                "Per-position modification frequency — AAGCCCG")
# (per-position occupancy + dominant-pair detail moved to the figure caption)
# (third column, top) GCCGGC summary. Reviewer C7: the previous version placed
# a red-bordered text box floating in an otherwise-empty axis, which reads as a
# detached, non-standard figure element. Replace it with a titled, borderless
# summary panel that mirrors the AAGCCCG summary panel directly below it, so the
# two right-column cells form a coherent "modification-type summary" column.
# GCCGGC breakdown bar — mirrors the AAGCCCG breakdown below so both right-column
# cells are parallel DATA plots (not a text panel). GCCGGC carries 4mC essentially
# only: T1 4mC = 1289 sites vs 6mA = 10 (~0).
axa=axes[0,2]
axa.bar([0,1],[1289,10],color=[C4,C6],width=0.6,zorder=3)
axa.text(0,1289+15,"1289",ha="center",fontsize=9,color=C4,fontweight="bold",zorder=4)
axa.text(1,10+15,"10",ha="center",fontsize=9,color=C6,fontweight="bold",zorder=4)
axa.set_xticks([0,1]);axa.set_xticklabels(["4mC","6mA"],fontsize=9)
axa.set_ylabel("GCCGGC sites (T1)",fontsize=8.5)
axa.set_title("GCCGGC",fontsize=9)
axa.set_ylim(0,1450)
for sp in ("top","right"): axa.spines[sp].set_visible(False)
# AAGCCCG breakdown at the corrected offsets (T1 canonical sites, cov>=10, >=50%):
# 4mC at C4 = 867 sites; 6mA at A0/A1 = 410 sites (312 + 98).
axb=axes[1,2]
axb.bar([0,1],[867,410],color=[C4,C6],width=0.6,zorder=3)
axb.text(0,867+12,"867",ha="center",fontsize=9,color=C4,fontweight="bold",zorder=4)
axb.text(1,410+12,"410",ha="center",fontsize=9,color=C6,fontweight="bold",zorder=4)
axb.set_xticks([0,1]);axb.set_xticklabels(["4mC","6mA"],fontsize=9)
axb.set_ylabel("AAGCCCG sites (T1)",fontsize=8.5);axb.set_title("AAGCCCG",fontsize=9)
axb.set_ylim(0, 980)
for sp in ("top","right"): axb.spines[sp].set_visible(False)
# (per-molecule dual co-occurrence 31.6%→7.2%→4.3% across T1/T2/T3 moved to the caption)
fig.legend(handles=[mpatches.Patch(color=C4,label="4mC"),mpatches.Patch(color=C6,label="6mA")],
           loc="lower center", frameon=False, fontsize=9, ncol=2, bbox_to_anchor=(0.5,-0.02))
fig.suptitle("Position of the methyl modification within each recognition motif", y=1.0, fontsize=12)
fig.tight_layout(rect=[0,0.03,1,1])
fig.savefig(OUT/"Figure_modification_position.pdf", bbox_inches="tight"); fig.savefig(OUT/"Figure_modification_position.png", dpi=150, bbox_inches="tight"); plt.close(fig)
print("[saved] Figure_modification_position (schematic + quantitative)")

# ============================================================
# FIGURE B — geographic distribution: landscape + quantitative, both motifs (Q14)
# ============================================================
def landscape(ax, df, tps, tp_lbls, color, title):
    """Fig1a-identical design: arm/core shading, oriC marker, region labels,
       stacked per-timepoint density tracks (filled bars), labelled axes."""
    edges = np.arange(0, GLEN+BIN, BIN); centers = (edges[:-1]+edges[1:])/2/1e6
    hmax = max((np.histogram(df[df.timepoint==tp]["position"].values, bins=edges)[0].max() or 1) for tp in tps)
    band_h = 0.80; ntp = len(tps); top = ntp
    # genome-region shading spanning all tracks (Fig1a colours)
    ax.add_patch(plt.Rectangle((0, 0), ARM_LEFT/1e6, top, fc=COL_ARM, alpha=0.18, ec='none', zorder=0))
    ax.add_patch(plt.Rectangle((ARM_LEFT/1e6, 0), (ARM_RIGHT-ARM_LEFT)/1e6, top, fc=COL_CORE, alpha=0.10, ec='none', zorder=0))
    ax.add_patch(plt.Rectangle((ARM_RIGHT/1e6, 0), (GLEN-ARM_RIGHT)/1e6, top, fc=COL_ARM, alpha=0.18, ec='none', zorder=0))
    ax.axvline(ORIC/1e6, color="#2E7D32", lw=1.2, ls='-', zorder=4)
    ax.text(ORIC/1e6, top+0.02, "oriC", ha="center", va="bottom", fontsize=7.5, color="#2E7D32", fontweight="bold")
    for bnd in (ARM_LEFT, ARM_RIGHT):
        ax.axvline(bnd/1e6, color="gray", ls=":", lw=0.6, zorder=1)
    for k, (tp, lab) in enumerate(zip(tps, tp_lbls)):
        base = ntp-1-k
        pos = df[df.timepoint==tp]["position"].values
        h,_ = np.histogram(pos, bins=edges)
        ax.bar(centers, h/hmax*band_h, width=BIN/1e6, bottom=base, color=color, alpha=0.95, linewidth=0, zorder=2)
        ax.axhline(base, color="#aaa", lw=0.5, zorder=1)
        # Track label INSIDE the panel (top-left of each track) so it never
        # collides with the y-axis title on the outer margin (figure-legibility-qc §3).
        ax.text(0.15, base+band_h-0.03, f"{lab}  (n={len(pos)})", ha="left", va="top",
                fontsize=7.2, color="#222", zorder=6)
    # region labels under the axis
    ax.text(ARM_LEFT/2/1e6, -0.34, "Left arm", ha="center", fontsize=7, color=COL_ARM, style="italic")
    ax.text((ARM_LEFT+ARM_RIGHT)/2/1e6, -0.34, "Core", ha="center", fontsize=7, color=COL_CORE, style="italic")
    ax.text((ARM_RIGHT+GLEN)/2/1e6, -0.34, "Right arm", ha="center", fontsize=7, color=COL_ARM, style="italic")
    ax.set_xlim(0, GLEN/1e6); ax.set_ylim(-0.5, top+0.35)
    ax.set_yticks([]); ax.spines['left'].set_visible(False)
    ax.set_xlabel("Chromosome position (Mb)", fontsize=9)
    ax.set_ylabel("HC methylation sites\n(per 50 kb, by timepoint)", fontsize=8)
    ax.set_title(title, fontsize=10)

def quant_core_arm(ax, df, tps, color, title):
    """Core fraction computed geographically with Fig1a's core (1.5–7.17 Mb)."""
    fr = []
    for tp in tps:
        p = df[df.timepoint==tp]["position"].values
        core = ((p>=ARM_LEFT)&(p<=ARM_RIGHT)).mean()*100 if len(p) else np.nan
        fr.append(core)
    ax.bar(range(len(tps)), fr, color=color, alpha=0.9)
    ax.axhline(((ARM_RIGHT-ARM_LEFT)/GLEN)*100, color="k", ls="--", lw=0.8, label="core size expectation (65%)")
    for i,v in enumerate(fr):
        if not np.isnan(v): ax.text(i, v+1, f"{v:.0f}%", ha="center", fontsize=9)
    ax.set_xticks(range(len(tps))); ax.set_xticklabels(tps)
    ax.set_ylabel("% of sites in core (1.5–7.17 Mb)"); ax.set_ylim(0,100)
    ax.set_title(title, fontsize=10); ax.legend(frameon=False, fontsize=7.5)

fig, ax = plt.subplots(2, 2, figsize=(13, 7), gridspec_kw={"width_ratios":[1.9,1]})
# GCCGGC (T1/T2/T3)
landscape(ax[0,0], g, ["T1","T2","T3"], ["T1 (12 h)","T2 (24 h)","T3 (50 h)"], C4,
          "(a) GCCGGC 4mC — linear chromosome landscape")
quant_core_arm(ax[0,1], g, ["T1","T2","T3"], C4,
          "(b) GCCGGC 4mC — core fraction")
# AAGCCCG (T1/T2/T3 from master site table)
landscape(ax[1,0], a, ["T1","T2","T3"], ["T1 (12 h)","T2 (24 h)","T3 (50 h)"], C6,
          "(c) AAGCCCG (4mC + 6mA pooled) — linear chromosome landscape")
quant_core_arm(ax[1,1], a, ["T1","T2","T3"], C6,
          "(d) AAGCCCG (4mC + 6mA pooled) — core fraction")
fig.suptitle("Geographic distribution of GCCGGC and AAGCCCG methylation across development (T1–T3)", y=1.0, fontsize=12)
fig.tight_layout(); fig.subplots_adjust(left=0.10, wspace=0.28)
fig.savefig(OUT/"Figure2_redistribution_4panel.pdf", bbox_inches="tight"); fig.savefig(OUT/"Figure2_redistribution_4panel.png", dpi=150); plt.close(fig)
print("[saved] Figure2_redistribution_4panel")
# report quant numbers
def geocore(p): return ((p>=ARM_LEFT)&(p<=ARM_RIGHT)).mean()*100 if len(p) else float('nan')
for name,df,tps in [("GCCGGC",g,["T1","T2","T3"]),("AAGCCCG",a,["T1","T2","T3"])]:
    print(name, {tp: f"{geocore(df[df.timepoint==tp]['position'].values):.0f}% core (n={len(df[df.timepoint==tp])})" for tp in tps})
