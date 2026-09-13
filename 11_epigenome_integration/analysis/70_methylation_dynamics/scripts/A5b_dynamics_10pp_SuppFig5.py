#!/usr/bin/env python3
"""A5b — per-site methylation dynamics with a 10 PERCENTAGE-POINT 'dynamic' threshold
and regeneration of Supplementary Figure 5 (2026-09-13, ledger item B03).

Background: A5_methylation_dynamics.py stores weighted_mod_freq in PERCENT and used
MAIN_THRESHOLD = 0.3, i.e. a range of 0.3 percentage points — which classifies
853/855 GCCGGC and 147/147 AAGCCCG sites as 'dynamic' (manuscript EN L49 '> 0.3').
Here the same wide tables are rebuilt from the same input and re-classified with
delta = max(T1,T2,T3) - min(T1,T2,T3) >= 10 percentage points (also > 10, and 5/20
pp for sensitivity). Non-monotonic share ('Valley (T2 min)' + 'Peak (T2 max)') is
computed among dynamic sites with A5's classify_pattern.

Figure: same three-block layout as the embedded fig_images/SuppFigure5.png
(violin + delta histogram; pattern classification; chromosomal distribution),
axis/legend text states 'percentage points'.

Outputs (tables/): A5b_dynamics_summary_pp.tsv, A5b_pattern_statistics_10pp.tsv,
A5b_GCCGGC_site_dynamics_10pp.tsv, A5b_AAGCCCG_site_dynamics_10pp.tsv
Figure (figures/): SuppFigure5_10pp.png  (copied to Writing/fig_images/SuppFigure5.png
by the caller after backing up the previous PNG).
"""
import importlib.util, sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from PIL import Image

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("A5", HERE / "A5_methylation_dynamics.py")
A5 = importlib.util.module_from_spec(spec); spec.loader.exec_module(A5)
OUT_T = A5.ANALYSIS_DIR / "tables"; OUT_F = A5.ANALYSIS_DIR / "figures"
THRESHOLDS_PP = [5, 10, 20]; MAIN_PP = 10
GENOME_LEN = 8_667_507

# ---- data (identical pipeline to A5.main steps 1-4) ----
genome = A5.load_genome(A5.GENOME_PATH)
df = pd.read_csv(A5.DATA_PATH)
uniq = df[["chrom","position","strand","mod_type"]].drop_duplicates().copy()
uniq["motif"] = A5.assign_motifs_fast(uniq, genome)
mm = uniq.set_index(["chrom","position","strand","mod_type"])["motif"]
df["motif"] = df.set_index(["chrom","position","strand","mod_type"]).index.map(mm).values
g = A5.build_wide_table(df, "4mC"); g = g[g.motif=="GCCGGC"].copy()
a = A5.build_wide_table(df, "6mA"); a = a[a.motif=="AAGCCCG"].copy()
assert len(g) == 855 and len(a) == 147, (len(g), len(a))
for d in (g, a):
    d["delta"] = d[["T1_freq","T2_freq","T3_freq"]].max(axis=1) - d[["T1_freq","T2_freq","T3_freq"]].min(axis=1)
    assert d[["T1_freq","T2_freq","T3_freq"]].max().max() > 1.5, "frequencies are not in percent"

# ---- sensitivity summary (>= and >) ----
rows = []
for thr in THRESHOLDS_PP:
    for name, d in (("GCCGGC", g), ("AAGCCCG", a)):
        rows.append(dict(threshold_pp=thr, motif=name, n_total=len(d),
                         n_dynamic_ge=int((d.delta >= thr).sum()), n_dynamic_gt=int((d.delta > thr).sum()),
                         pct_dynamic_ge=round((d.delta >= thr).mean()*100, 1)))
summ = pd.DataFrame(rows); summ.to_csv(OUT_T / "A5b_dynamics_summary_pp.tsv", sep="\t", index=False)

# ---- main classification at 10 pp (>=) ----
pat_rows = []
for name, d in (("GCCGGC", g), ("AAGCCCG", a)):
    d["category"] = np.where(d.delta >= MAIN_PP, "dynamic", "stable")
    dyn = d.category == "dynamic"
    d.loc[dyn, "pattern"] = d[dyn].apply(A5.classify_pattern, axis=1); d.loc[~dyn, "pattern"] = "stable"
    d["region"] = d.position.apply(A5.assign_region)
    d.sort_values("delta", ascending=False)[["chrom","position","strand","T1_freq","T2_freq","T3_freq","delta","category","pattern","region"]]\
        .to_csv(OUT_T / f"A5b_{name}_site_dynamics_10pp.tsv", sep="\t", index=False)
    vc = d.loc[dyn, "pattern"].value_counts(); n_dyn = int(dyn.sum())
    for p, c in vc.items():
        pat_rows.append(dict(motif=name, pattern=p, count=int(c), pct_of_dynamic=round(c/n_dyn*100, 1)))
    nm = int(vc.get("Valley (T2 min)", 0) + vc.get("Peak (T2 max)", 0))
    pat_rows.append(dict(motif=name, pattern="NON-MONOTONIC (Valley+Peak)", count=nm, pct_of_dynamic=round(nm/n_dyn*100, 1)))
    pat_rows.append(dict(motif=name, pattern="DYNAMIC TOTAL", count=n_dyn, pct_of_dynamic=100.0))
pats = pd.DataFrame(pat_rows); pats.to_csv(OUT_T / "A5b_pattern_statistics_10pp.tsv", sep="\t", index=False)
print(summ.to_string(index=False)); print(pats.to_string(index=False))

# ---- figure blocks (same layout as the embedded SuppFigure5.png) ----
C = A5.COLORS; PC = A5.PATTERN_COLORS
tp_cols, tp_lab = ["T1_freq","T2_freq","T3_freq"], ["T1","T2","T3"]
def style(ax): ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

def block1():
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("Methylation frequency dynamics: per-site trajectories", fontsize=13, fontweight="bold")
    for ax, d, lab in ((axes[0,0], g, "GCCGGC 4mC\n(T1→T2→T3 methylation frequency)"), (axes[0,1], a, "AAGCCCG 6mA\n(T1→T2→T3 methylation frequency)")):
        vp = ax.violinplot([d[c].values for c in tp_cols], positions=[1,2,3], showmedians=True, showextrema=True, widths=0.6)
        for body, tp in zip(vp["bodies"], tp_lab): body.set_facecolor(C[tp]); body.set_alpha(0.7)
        vp["cmedians"].set_color("black"); vp["cmedians"].set_linewidth(2)
        ax.set_xticks([1,2,3]); ax.set_xticklabels(tp_lab); ax.set_ylabel("Methylation frequency (%)", fontsize=10)
        ax.set_title(lab, fontsize=10); ax.set_ylim(0, 105)
        ax.text(0.98, 0.98, f"n={len(d)} sites\n(all 3 TPs present)", transform=ax.transAxes, ha="right", va="top", fontsize=8,
                bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5)); style(ax)
    for ax, d, lab in ((axes[1,0], g, "GCCGGC 4mC — range distribution"), (axes[1,1], a, "AAGCCCG 6mA — range distribution")):
        ax.hist(d.delta.values, bins=40, color="#6B7783", alpha=0.8, edgecolor="white")
        for thr, col, ls in ((5, "#C0803A", "--"), (10, "#A64B44", "-"), (20, "#6E5495", ":")):
            n_ab = int((d.delta >= thr).sum())
            ax.axvline(thr, color=col, linewidth=1.5, linestyle=ls, label=f"≥ {thr} pp: {n_ab} ({n_ab/len(d)*100:.0f}%)")
        ax.set_xlabel("Range of methylation frequency across T1–T3, max − min (percentage points)", fontsize=10)
        ax.set_ylabel("Number of sites", fontsize=10); ax.set_title(lab, fontsize=10)
        ax.legend(fontsize=8, title="Dynamic threshold (percentage points)"); style(ax)
    fig.tight_layout(); out = OUT_F / "A5b_dynamics_violin_pp.png"; fig.savefig(out, dpi=300, bbox_inches="tight"); plt.close(fig); return out

def block2():
    gd, ad = g[g.category=="dynamic"], a[a.category=="dynamic"]
    fig = plt.figure(figsize=(14, 8)); gs = GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.35)
    order = ["Monotone increase","Monotone decrease","Peak (T2 max)","Valley (T2 min)","Other"]
    for ax, d, lab in ((fig.add_subplot(gs[0,0]), gd, f"GCCGGC (4mC), dynamic n={len(gd)}"), (fig.add_subplot(gs[1,0]), ad, f"AAGCCCG (6mA), dynamic n={len(ad)}")):
        cnt = d.pattern.value_counts(); y = [int(cnt.get(p, 0)) for p in order]
        bars = ax.barh(order, y, color=[PC[p] for p in order], edgecolor="white")
        for b, v in zip(bars, y): ax.text(b.get_width()+0.3, b.get_y()+b.get_height()/2, f"{v} ({v/max(sum(y),1)*100:.0f}%)", va="center", fontsize=8)
        ax.set_title(lab, fontsize=10); ax.set_xlabel("Count", fontsize=9); style(ax)
    for ax, pat in zip((fig.add_subplot(gs[0,1]), fig.add_subplot(gs[0,2]), fig.add_subplot(gs[1,1]), fig.add_subplot(gs[1,2])), order[:4]):
        plotted = False
        for d, lab, ls in ((gd, "GCCGGC", "-"), (ad, "AAGCCCG", "--")):
            sub = d[d.pattern == pat]
            if len(sub):
                rep = sub.nlargest(1, "delta").iloc[0]
                ax.plot(tp_lab, [rep[c] for c in tp_cols], marker="o", linestyle=ls, color=PC[pat], label=lab, linewidth=1.8); plotted = True
        ax.set_title(pat, fontsize=9, color=PC[pat], fontweight="bold"); ax.set_ylabel("Methylation frequency (%)", fontsize=8); ax.set_ylim(0, 105)
        ax.legend(fontsize=7); style(ax)
        if not plotted: ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center", va="center", fontsize=10, color="grey")
    fig.suptitle(f"Dynamic site pattern classification (dynamic = range ≥ {MAIN_PP} percentage points)", fontsize=12, fontweight="bold")
    out = OUT_F / "A5b_pattern_classification_10pp.png"; fig.savefig(out, dpi=300, bbox_inches="tight"); plt.close(fig); return out

def block3():
    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
    fig.suptitle(f"Chromosomal distribution of methylation sites\n(dynamic = range ≥ {MAIN_PP} percentage points, red; stable, grey)", fontsize=12, fontweight="bold")
    bins = np.arange(0, GENOME_LEN + 100_000, 100_000)
    for ax, d, lab in ((axes[0], g, "GCCGGC (4mC)"), (axes[1], a, "AAGCCCG (6mA)")):
        for cat, col, al in (("stable", C["stable"], 0.7), ("dynamic", C["dynamic"], 0.9)):
            sub = d[d.category == cat]; ax.hist(sub.position, bins=bins, color=col, alpha=al, label=f"{cat} (n={len(sub)})")
        for b in (A5.ARM_LEFT_END, A5.ARM_RIGHT_START): ax.axvline(b, color="navy", linewidth=1.2, linestyle="--", alpha=0.7)
        yt = ax.get_ylim()[1]*0.9
        for x, t in ((A5.ARM_LEFT_END/2, "arm"), ((A5.ARM_LEFT_END+A5.ARM_RIGHT_START)/2, "core"), ((A5.ARM_RIGHT_START+GENOME_LEN)/2, "arm")):
            ax.text(x, yt, t, ha="center", fontsize=9, color="navy")
        ax.set_ylabel("Site count per 100 kb", fontsize=9); ax.set_title(lab, fontsize=10); ax.legend(fontsize=8); style(ax)
    axes[1].set_xlabel("Chromosomal position (Mb)", fontsize=10)
    axes[1].xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, _: f"{x/1e6:.1f}"))
    fig.tight_layout(); out = OUT_F / "A5b_chromosomal_distribution_10pp.png"; fig.savefig(out, dpi=300, bbox_inches="tight"); plt.close(fig); return out

parts = [Image.open(p).convert("RGB") for p in (block1(), block2(), block3())]
W = max(im.width for im in parts); H = sum(int(im.height * W / im.width) for im in parts)
canvas = Image.new("RGB", (W, H), "white"); y = 0
for im in parts:
    im2 = im.resize((W, int(im.height * W / im.width)), Image.LANCZOS); canvas.paste(im2, (0, y)); y += im2.height
final = OUT_F / "SuppFigure5_10pp.png"; canvas.save(final, dpi=(300, 300)); print("FIGURE", final, canvas.size)
