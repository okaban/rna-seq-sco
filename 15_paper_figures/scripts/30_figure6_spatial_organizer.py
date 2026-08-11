#!/usr/bin/env python3
"""
Figure 6 (real-data redesign, 2026-07-06): Methylation as a permissive spatial
organiser of the regulatory genome.

Replaces the earlier random-lollipop schematic (site positions were drawn from
RNG.normal and carried no quantitative content). This version renders the
paper's capstone thesis directly from data:

  Panel A  GCCGGC 4mC sites reverse their core/arm partition across development:
           T1  1,289 sites, 83% core   (vegetative growth, 12 h)
           T2    407 sites, 18% core / 82% arm   (developmental switch, 24 h)
           T3     21 sites, 38% core
           (37_defense_island_GCCGGC/tables/GCCGGC_sites_by_timepoint.tsv;
            core = 1.5-7.17 Mb, the manuscript definition)

  Panel B  The redistribution is PERMISSIVE: promoter GCCGGC occupancy (+-2 kb,
           T1) is only weakly, region-controlled, associated with the T1->T2
           expression change (region-controlled partial rank r = -0.09,
           p = 0.004, n = 1019;
           52_shielded_exposed_boundary/tables/all_genes_features_unified_n57.tsv).

A compact core->arm switch schematic sits above Panel A, and the
"permissive, not instructive" statement is retained as a bottom banner.

Outputs:
  15_paper_figures/figures/main/Figure6_spatial_organizer.{pdf,svg,png}
  Writing/fig_images/Figure6.png   (manuscript image slot; no renumbering)
"""
from pathlib import Path
import importlib
import sys

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch

sys.path.insert(0, str(Path(__file__).parent))
_utils = importlib.import_module('00_shared_utils')
for _attr in dir(_utils):
    if not _attr.startswith('_'):
        globals()[_attr] = getattr(_utils, _attr)

# ── Data locations ────────────────────────────────────────────────────────────
BASE_ANALYSIS = Path.home() / 'bioinfo' / 'rna-seq' / '11_epigenome_integration' / 'analysis'
SITES_TSV = BASE_ANALYSIS / '37_defense_island_GCCGGC' / 'tables' / 'GCCGGC_sites_by_timepoint.tsv'
REG_TSV = BASE_ANALYSIS / '52_shielded_exposed_boundary' / 'tables' / 'all_genes_features_unified_n57.tsv'

GENOME = 8_667_507                 # NC_003888.3 length (bp)
CORE_LO, CORE_HI = 1.5e6, 7.17e6   # manuscript core definition (Methods)

# ── Colours (threaded from shared palette) ────────────────────────────────────
COL_CORE = '#3E7256'      # core fraction — unified deep muted green
COL_ARM  = '#C0803A'      # arm fraction — unified calm amber


def _load_panelA():
    """Core/arm site counts per timepoint from real GCCGGC calls."""
    g = pd.read_csv(SITES_TSV, sep='\t')
    g['core'] = (g['position'] >= CORE_LO) & (g['position'] <= CORE_HI)
    rows = {}
    for tp in ('T1', 'T2', 'T3'):
        sub = g[g['timepoint'] == tp]
        n = len(sub); nc = int(sub['core'].sum())
        rows[tp] = dict(n=n, core=nc, arm=n - nc, pct_core=100 * nc / n)
    return g, rows


def _load_panelB(g):
    """Promoter GCCGGC occupancy (+-2 kb, T1) vs LFC_T2vsT1, region-controlled."""
    reg = pd.read_csv(REG_TSV, sep='\t').dropna(subset=['tss']).copy()
    gT1 = g[g['timepoint'] == 'T1']
    P = np.sort(gT1['position'].values)
    Fq = gT1.sort_values('position')['frequency'].values

    def occ(tss, w=2000):
        m = (P >= tss - w) & (P <= tss + w)
        return Fq[m].sum() if m.any() else 0.0

    reg['occ2k'] = reg['tss'].apply(occ)
    d = reg[['occ2k', 'LFC_T2vsT1', 'region']].dropna()
    rx = d['occ2k'].rank() - np.polyval(np.polyfit(d['region'].rank(), d['occ2k'].rank(), 1), d['region'].rank())
    ry = d['LFC_T2vsT1'].rank() - np.polyval(np.polyfit(d['region'].rank(), d['LFC_T2vsT1'].rank(), 1), d['region'].rank())
    rr, pp = pearsonr(rx, ry)
    return d, rr, pp


# ── Panel A: core/arm reversal ────────────────────────────────────────────────
def draw_panelA(ax, rows):
    tps = ['T1', 'T2', 'T3']
    xlabels = ['T1\n(12 h)', 'T2\n(24 h)', 'T3\n(50 h)']
    core_pct = [rows[t]['pct_core'] for t in tps]
    arm_pct = [100 - c for c in core_pct]
    x = np.arange(3)
    ax.bar(x, core_pct, width=0.62, color=COL_CORE, label='Core (1.5–7.17 Mb)')
    ax.bar(x, arm_pct, width=0.62, bottom=core_pct, color=COL_ARM, label='Arms')
    # site-count annotation above each bar
    for i, t in enumerate(tps):
        ax.text(i, 103, f"n = {rows[t]['n']:,}", ha='center', va='bottom',
                fontsize=6, color=COL_DARK)
    # core-% value on the core segment (headline number)
    for i, t in enumerate(tps):
        cp = rows[t]['pct_core']
        ax.text(i, cp / 2, f"{cp:.0f}%", ha='center', va='center',
                fontsize=7, color='white', fontweight='bold')
    ax.set_xticks(x); ax.set_xticklabels(xlabels, fontsize=7)
    ax.set_ylim(0, 112)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylabel('GCCGGC 4mC sites (%)', fontsize=8)
    ax.set_title('Core→arm redistribution of GCCGGC 4mC', fontsize=8, pad=14)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(frameon=False, fontsize=6, loc='lower center',
              bbox_to_anchor=(0.5, -0.30), ncol=2, handlelength=1.1,
              columnspacing=1.2, borderpad=0.2)


# ── Panel B: permissive scatter ───────────────────────────────────────────────
def draw_panelB(ax, d, rr, pp):
    ax.scatter(d['occ2k'], d['LFC_T2vsT1'], s=7, alpha=0.35,
               c=COL_DARK, edgecolors='none', zorder=2)
    ax.axhline(0, color='k', lw=0.5, zorder=1)
    ax.set_xlabel('Promoter GCCGGC occupancy (±2 kb, T1)', fontsize=8)
    ax.set_ylabel('log$_2$ FC (T2 vs T1)', fontsize=8)
    ax.set_title('Methylation biases expression only weakly', fontsize=8, pad=14)
    ax.spines[['top', 'right']].set_visible(False)
    ax.margins(x=0.04)
    # stat annotation in clear upper-right whitespace, boxed, no data overlap
    ax.text(0.97, 0.96,
            f"region-controlled\npartial rank $r$ = {rr:.2f}\n$p$ = {pp:.3f}  ($n$ = {len(d):,})",
            transform=ax.transAxes, ha='right', va='top', fontsize=6.2,
            color=COL_DARK,
            bbox=dict(boxstyle='round,pad=0.35', facecolor='white',
                      edgecolor=COL_GRAY, linewidth=0.7))


# ── Top schematic strip: core → arm switch ────────────────────────────────────
def draw_schematic(ax):
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis('off')

    def chromosome(cx, concentrate_core, label, sub):
        w, h, y = 3.6, 0.55, 5.4
        x0 = cx - w / 2
        # backbone
        ax.add_patch(Rectangle((x0, y), w, h, facecolor='#ECEFF1',
                               edgecolor=COL_DARK, linewidth=0.9, zorder=2))
        # core shade (central ~65% of the bar)
        core_frac = (CORE_HI - CORE_LO) / GENOME
        cw = w * core_frac
        ax.add_patch(Rectangle((cx - cw / 2, y), cw, h, facecolor='#CFE3D6',
                               edgecolor='none', zorder=1.5))
        # methyl ticks: core-concentrated (T1) or arm-concentrated (T2)
        rng = np.random.default_rng(3)
        if concentrate_core:
            xs = np.clip(rng.normal(cx, cw / 3.2, 22), x0 + 0.05, x0 + w - 0.05)
        else:
            left = rng.normal(x0 + (w - cw) / 4, 0.28, 9)
            right = rng.normal(x0 + w - (w - cw) / 4, 0.28, 9)
            xs = np.clip(np.concatenate([left, right]), x0 + 0.05, x0 + w - 0.05)
        for xt in xs:
            ax.vlines(xt, y + h, y + h + 0.55, color=COL_4mC, lw=0.6, alpha=0.85, zorder=3)
        ax.text(cx, y - 0.35, label, ha='center', va='top', fontsize=6.6,
                fontweight='bold', color=COL_DARK)
        ax.text(cx, y - 1.15, sub, ha='center', va='top', fontsize=5.6, color=COL_DARK)

    chromosome(2.4, True, 'T1 — core-methylated', 'vegetative growth (12 h)')
    chromosome(7.6, False, 'T2 — arm-redistributed', 'developmental switch (24 h)')
    # arrow between the two states
    arr = FancyArrowPatch((4.35, 5.68), (5.65, 5.68), arrowstyle='-|>',
                          mutation_scale=13, lw=1.8, color=COL_DARK, zorder=5)
    ax.add_patch(arr)
    ax.text(5.0, 6.55, 'developmental\nswitch', ha='center', va='bottom',
            fontsize=5.8, fontweight='bold', color=COL_DARK)
    # mark legend, top-left whitespace
    ax.text(0.1, 9.4, '│ GCCGGC 4mC site', fontsize=5.8, color=COL_4mC,
            ha='left', va='top')


def main():
    apply_style()
    print('=== Figure 6: permissive spatial organiser (real-data redesign) ===')

    g, rowsA = _load_panelA()
    d, rr, pp = _load_panelB(g)
    print(f'  Panel A: ' + ', '.join(f"{t} n={rowsA[t]['n']} core={rowsA[t]['pct_core']:.0f}%" for t in ('T1','T2','T3')))
    print(f'  Panel B: region-controlled r={rr:.3f} p={pp:.4f} n={len(d)}')

    fig = plt.figure(figsize=(mm_to_inch(174), mm_to_inch(100)))

    # layout: top schematic strip, two data panels, bottom banner
    ax_s = fig.add_axes([0.02, 0.74, 0.96, 0.24]); 
    axA = fig.add_axes([0.085, 0.155, 0.37, 0.50])
    axB = fig.add_axes([0.60, 0.155, 0.37, 0.50])

    draw_schematic(ax_s)
    draw_panelA(axA, rowsA)
    draw_panelB(axB, d, rr, pp)

    # panel letters — uniform offset-points placement (shared helper)
    for ax, L in ((axA, 'a'), (axB, 'b')):
        _utils.add_panel_label(ax, L)

    # The former yellow "Permissive, not instructive" interpretation banner was
    # removed (figure-legibility-qc §5: in-plot interpretation prose belongs in the
    # caption). Its content is fully covered by the Figure 6 legend text; the
    # load-bearing r = -0.09 statistic remains boxed in panel b.

    out = FIG_DIR / 'Figure6_spatial_organizer'
    save_figure(fig, out, formats=('pdf', 'svg', 'png'))

    # write the manuscript image slot (keep Figure6.png; no renumbering)
    import shutil
    png = out.with_suffix('.png')
    obsidian_slot = Path.home() / 'obsidian' / 'Research' / 'rna-seq' / 'Writing' / 'fig_images' / 'Figure6.png'
    targets = [obsidian_slot] if obsidian_slot.parent.is_dir() else []
    targets.append(BASE / 'Writing' / 'fig_images' / 'Figure6.png')
    for wr in targets:
        wr.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(png, wr)
        print(f'  Copied → {wr}')
    print('=== Done ===')


if __name__ == '__main__':
    main()
