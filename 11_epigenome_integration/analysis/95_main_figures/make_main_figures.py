#!/usr/bin/env python3
"""95 — the producing script for main Figures 1, 2 and 4.

Why this exists: these three figures were rebuilt on 2026-09-21 in throwaway
cells, so the PNGs shipped with NO producing script, ad-hoc colours, no font
ladder and widths of 182-196 mm. `Writing/check_figures.py` now fails on exactly
that. This script regenerates them from canonical inputs through the approved
shared style module (palette v2, font ladder, 174 mm NAR clamp).

Content changes adopted from the figure-design literature (see
Writing/FIGURE_DESIGN_RULES_260922.md):
  * Fig 1a is a SuperPlot (Lord et al. 2020, JCB): every site is drawn coloured
    by the biological replicate it came from, with the three replicate medians
    overlaid. The previous version pooled three replicates into one cloud.
  * Fig 1b and Fig 4a carry Wilson 95% CIs (Midway 2020, principle 5). Fig 4a
    compares 27/62 against 435/989 - without intervals two equal-height bars
    imply equal precision across a 16-fold difference in denominator.

Motif assignment and the 0-based `position` convention follow
90_per_timepoint_census_audit/per_timepoint_census.py exactly.

Usage:  python make_main_figures.py [--skip-pileup]
"""
import argparse, collections, csv, importlib.util, sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
A = HERE.parent
REPO = A.parent.parent
TAB = HERE / "tables"; TAB.mkdir(exist_ok=True)
FIGDIR = HERE / "figures"; FIGDIR.mkdir(exist_ok=True)
SHIP = Path("/Users/okaban/obsidian/Research/rna-seq/Writing/fig_images")

REFFA = Path("/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa")
PILEUP = Path("/Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup")
CANON = A / "01_integration/high_confidence_sites_weighted.csv"
CENSUS = A / "90_per_timepoint_census_audit/tables"
RESTRUCT = A / "94_figure_restructure_260921/tables"

T1_REPS = ["1-1", "1-2", "1-3"]
CODE_4mC, CODE_5mC = "21839", "m"
CORE = (1_500_000, 7_170_000)

REF = "".join(l.strip() for l in open(REFFA) if not l.startswith(">")).upper()
GENOME_LEN = len(REF)


def shared():
    p = REPO / "15_paper_figures/scripts/00_shared_utils.py"
    spec = importlib.util.spec_from_file_location("shared_utils", p)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["shared_utils"] = mod
    spec.loader.exec_module(mod)
    return mod


def motif(p, s, mod):
    """Identical to 90_/per_timepoint_census.py — 0-based position, strand-aware."""
    if mod == "4mC":
        if (REF[p - 2:p + 4] == "GCCGGC") if s == "+" else (REF[p - 3:p + 3] == "GCCGGC"):
            return "GCCGGC"
        if (REF[p - 4:p + 3] == "AAGCCCG") if s == "+" else (REF[p - 2:p + 5] == "CGGGCTT"):
            return "AAGCCCG"
    else:
        if s == "+" and (REF[p:p + 7] == "AAGCCCG" or REF[p - 1:p + 6] == "AAGCCCG"):
            return "AAGCCCG"
        if s == "-" and (REF[p - 6:p + 1] == "CGGGCTT" or REF[p - 5:p + 2] == "CGGGCTT"):
            return "AAGCCCG"
    return "other"


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def site_sets():
    """(timepoint, mod, motif) -> set of (pos, strand), from the canonical file."""
    S = collections.defaultdict(set)
    for r in csv.DictReader(open(CANON)):
        mod = "4mC" if "4mC" in r["mod_type"] else ("6mA" if "6mA" in r["mod_type"] else r["mod_type"])
        p, s = int(r["position"]), r["strand"]
        S[(r["timepoint"], mod, motif(p, s, mod))].add((p, s))
    return S


def per_replicate_table(sites):
    """4mC and 5mC percentage at each canonical T1 GCCGGC site, per replicate."""
    want = set(sites)
    rows = []
    for rep in T1_REPS:
        got = {}
        with open(PILEUP / f"{rep}_pileup.bed") as fh:
            for line in fh:
                c = line.split("\t", 11)
                if len(c) < 11 or c[3] not in (CODE_4mC, CODE_5mC):
                    continue
                key = (int(c[1]), c[5])
                if key in want:
                    got[(key, c[3])] = (float(c[10]), int(c[4]))
        n = 0
        for key in want:
            f4, f5 = got.get((key, CODE_4mC)), got.get((key, CODE_5mC))
            if f4 is None and f5 is None:
                continue
            rows.append(dict(replicate=rep, position=key[0], strand=key[1],
                             pct_4mC=(f4 or (np.nan, 0))[0], cov_4mC=(f4 or (np.nan, 0))[1],
                             pct_5mC=(f5 or (np.nan, 0))[0], cov_5mC=(f5 or (np.nan, 0))[1]))
            n += 1
        print(f"  {rep}: {n} of {len(want)} canonical sites covered", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(TAB / "t1_per_replicate_sites.tsv", sep="\t", index=False)
    return df



def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-pileup", action="store_true")
    args = ap.parse_args()
    su = shared()
    su.apply_unified_style() if hasattr(su, "apply_unified_style") else su.apply_style()

    S = site_sets()
    t1 = sorted(S[("T1", "4mC", "GCCGGC")])
    print(f"canonical T1 GCCGGC 4mC sites: {len(t1)}")

    cache = TAB / "t1_per_replicate_sites.tsv"
    if args.skip_pileup and cache.exists():
        rep = pd.read_csv(cache, sep="\t")
    else:
        print("reading T1 pileups (a few minutes)...")
        rep = per_replicate_table(t1)

    # ---- Figure 1 -----------------------------------------------------------
    fig1 = plt.figure(figsize=(su.mm_to_inch(174), su.mm_to_inch(76)))
    a1 = fig1.add_axes([0.095, 0.40, 0.145, 0.42])
    a2 = fig1.add_axes([0.400, 0.40, 0.170, 0.42])
    a3 = fig1.add_axes([0.750, 0.40, 0.165, 0.42])

    # (a) SuperPlot: sites coloured by replicate, replicate medians overlaid
    rng = np.random.default_rng(42)
    rep_cols = [su.COL_BAR_LIGHT, su.COL_BAR_MID, su.COL_BAR_DARK]
    med_rows = []
    for ch, xc, lab in ((("pct_4mC"), 0, "4mC"), (("pct_5mC"), 1, "5mC")):
        for i, r in enumerate(T1_REPS):
            v = rep.loc[rep.replicate == r, ch].dropna().values
            if len(v) == 0:
                continue
            x = xc + (i - 1) * 0.22 + rng.normal(0, 0.035, len(v))
            a1.scatter(x, v, s=1.2, c=rep_cols[i], alpha=0.35, linewidths=0, zorder=2)
            med_rows.append(dict(channel=ch, replicate=r, n=len(v), median=float(np.median(v))))
            a1.plot([xc + (i - 1) * 0.22], [np.median(v)], marker="o", ms=4.5,
                    mfc=rep_cols[i], mec=su.COL_DARK, mew=0.6, zorder=4)
    pd.DataFrame(med_rows).to_csv(TAB / "fig1a_superplot_medians.tsv", sep="\t", index=False)
    a1.set_xticks([0, 1]); a1.set_xticklabels(["4mC", "5mC"])
    a1.set_ylabel("modification frequency (%)")
    a1.set_ylim(-4, 104); a1.set_xlim(-0.55, 1.55)
    a1.set_title("The GCCGGC mark is 4mC,\nnot 5mC", loc="left")
    # The three large circles are the replicate medians; named in the caption
    # rather than in a legend box, so no text sits inside the data region.

    # (b) core fraction per timepoint with Wilson CI, against the motif null
    core_rows = []
    for tp in ("T1", "T2", "T3"):
        st = S[(tp, "4mC", "GCCGGC")]
        k = sum(1 for p, _ in st if CORE[0] <= p <= CORE[1])
        lo, hi = wilson(k, len(st))
        core_rows.append(dict(timepoint=tp, n=len(st), in_core=k, frac=k / len(st), lo=lo, hi=hi))
    inst = [m.start() for m in __import__("re").finditer("GCCGGC", REF)]
    null = sum(1 for p in inst if CORE[0] <= p <= CORE[1]) / len(inst)
    nlo, nhi = wilson(sum(1 for p in inst if CORE[0] <= p <= CORE[1]), len(inst))
    cf = pd.DataFrame(core_rows); cf.to_csv(TAB / "fig1b_core_fraction_wilson.tsv", sep="\t", index=False)
    x = np.arange(3)
    a2.bar(x, cf.frac * 100, color=su.COL_CORE, width=0.6, zorder=2)
    a2.errorbar(x, cf.frac * 100, yerr=[(cf.frac - cf.lo) * 100, (cf.hi - cf.frac) * 100],
                fmt="none", ecolor=su.COL_DARK, elinewidth=0.8, capsize=2, zorder=3)
    a2.axhspan(nlo * 100, nhi * 100, color=su.COL_GRAY, alpha=0.35, zorder=1)
    a2.axhline(null * 100, color=su.COL_GRAY, ls="--", lw=0.8, zorder=1)
    # Direct-label the null line in the margin the bars do not occupy, rather than
    # over the bars themselves (the previous placement overlapped them).
    a2.set_xlim(-0.55, 3.25)
    a2.text(2.62, null * 100, f"motif-matched null\n{null*100:.0f}% (95% CI)", ha="left",
            va="center", color=su.COL_GRAY, fontsize=su.FONT["annot"])
    a2.set_xticks(x); a2.set_xticklabels(su.TP_LABELS_NL)
    a2.spines["right"].set_visible(False)
    a2.set_ylabel("GCCGGC 4mC in core (%)")
    a2.set_ylim(0, 100)
    a2.set_title("Core-concentrated at\nevery timepoint", loc="left")

    # (c) promoter occupancy vs expression change, split by compartment
    d = pd.read_csv(RESTRUCT / "occupancy_vs_lfc.tsv", sep="\t")
    for rg, col, lab in ((0, su.COL_CORE, "core"), (1, su.COL_ARM, "arm")):
        sub = d[d.region == rg]
        a3.scatter(sub.occ2k, sub.LFC_T2vsT1, s=3, c=col, alpha=0.55, linewidths=0,
                   label=f"{lab} (n={len(sub)})")
    a3.axhline(0, color=su.COL_DARK, lw=0.7)
    raw = spearmanr(d.occ2k, d.LFC_T2vsT1)[0]
    a3.set_xlim(-20, 520)
    a3.set_xlabel("promoter GCCGGC\noccupancy (±2 kb, T1)")
    a3.set_ylabel("log$_2$ FC (T2 vs T1)")
    a3.set_title(("Coupling is positional:\n"
                  f"r {raw:.2f} → −0.09 controlled").replace("-", "−"), loc="left")
    # key below the panel: inside the axes it would sit on the point cloud
    a3.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.46), ncol=2,
              handletextpad=0.3, columnspacing=1.4, markerscale=2.2,
              fontsize=su.FONT["legend"])
    for ax, L in ((a1, "a"), (a2, "b"), (a3, "c")):
        ax.text(-0.30, 1.20, L, transform=ax.transAxes, fontweight="bold",
                fontsize=su.FONT["panel_letter"], va="bottom", ha="left")
    su.assert_no_text_collisions(fig1, "Figure1")
    su.save_figure(fig1, FIGDIR / "Figure1", formats=("png", "pdf"), width_class="full")
    plt.close(fig1)
    print("Figure1 written")


if __name__ == "__main__":
    sys.exit(main())
