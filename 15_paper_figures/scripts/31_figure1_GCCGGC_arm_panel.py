#!/usr/bin/env python3
"""
Figure 1 supplementary panel: GCCGGC core->arm geographic redistribution
across T1/T2/T3, contrasted with AAGCCCG (stable distribution).

Two sub-panels share the y-axis (sites per Mb):
  Left:  GCCGGC (4mC) — dramatic core-to-arm shift (T1 core dominant ->
                       T2 arm dominant -> T3 collapsed)
  Right: AAGCCCG (6mA) — flat core/arm ratio across all timepoints
                        (negative control / contrast)

Outputs:
  15_paper_figures/figures/main/Figure1_panel_GCCGGC_distribution.{pdf,svg,png}
  Writing/fig_images/Figure1_GCCGGC_arm_panel.png

Data sources:
  GCCGGC sites: 11_epigenome_integration/analysis/37_defense_island_GCCGGC/
                tables/GCCGGC_sites_by_timepoint.tsv  (n=1718, region pre-assigned)
  AAGCCCG sites: 11_epigenome_integration/analysis/23_expanded_motif_search/
                 6mA_final_census.csv  (filtered final_motif == 'AAGCCCG')

Core/arm definition (consistent with 01_figure1_landscape.py and 10_new_figure1_overview.py):
  Genome:    8,667,507 bp
  Left arm:  0 – 1,500,000 bp           (1.5 Mb)
  Core:      1,500,000 – 7,167,507 bp   (5.667 Mb)
  Right arm: 7,167,507 – 8,667,507 bp   (1.5 Mb)
  Total arm length: 3.0 Mb; total core length: 5.667 Mb
"""

import sys
import importlib
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).parent))
_utils = importlib.import_module('00_shared_utils')
for _attr in dir(_utils):
    if not _attr.startswith('_'):
        globals()[_attr] = getattr(_utils, _attr)

# ── Genome landmarks (must match Figure 1 main panel) ──────────────────────
GENOME_LEN = 8_667_507
ARM_LEFT = 1_500_000
ARM_RIGHT = 7_167_507
ARM_LEN_MB = (ARM_LEFT + (GENOME_LEN - ARM_RIGHT)) / 1e6        # 3.0 Mb
CORE_LEN_MB = (ARM_RIGHT - ARM_LEFT) / 1e6                       # 5.667 Mb

COL_CORE = '#4169E1'   # Royal blue
COL_ARM = '#DAA520'    # Goldenrod
TIMEPOINTS = ['T1', 'T2', 'T3']

# Output paths
WRITING_FIG_DIR = BASE / 'Writing' / 'fig_images'


def load_gccggc_region_counts():
    """Per-timepoint, per-region GCCGGC site counts."""
    path = (EPIGENOME / '37_defense_island_GCCGGC' / 'tables' /
            'GCCGGC_sites_by_timepoint.tsv')
    df = pd.read_csv(path, sep='\t')
    counts = (df.groupby(['timepoint', 'region']).size()
              .unstack(fill_value=0)
              .reindex(index=TIMEPOINTS, columns=['core', 'arm'], fill_value=0))
    print(f'  GCCGGC counts (core / arm):')
    for tp in TIMEPOINTS:
        c, a = counts.loc[tp, 'core'], counts.loc[tp, 'arm']
        tot = c + a
        pct = 100 * c / tot if tot else 0
        print(f'    {tp}: core={c:>5}  arm={a:>5}  total={tot:>5}  '
              f'(core {pct:.1f}%)')
    return counts


def load_aagcccg_region_counts():
    """Per-timepoint, per-region AAGCCCG site counts (computed from census)."""
    path = (EPIGENOME / '23_expanded_motif_search' / '6mA_final_census.csv')
    df = pd.read_csv(path)
    df = df[df['final_motif'] == 'AAGCCCG'].copy()
    df['region'] = np.where(
        (df['position'] >= ARM_LEFT) & (df['position'] <= ARM_RIGHT),
        'core', 'arm')
    counts = (df.groupby(['timepoint', 'region']).size()
              .unstack(fill_value=0)
              .reindex(index=TIMEPOINTS, columns=['core', 'arm'], fill_value=0))
    print(f'  AAGCCCG counts (core / arm):')
    for tp in TIMEPOINTS:
        c, a = counts.loc[tp, 'core'], counts.loc[tp, 'arm']
        tot = c + a
        pct = 100 * c / tot if tot else 0
        print(f'    {tp}: core={c:>5}  arm={a:>5}  total={tot:>5}  '
              f'(core {pct:.1f}%)')
    return counts


def counts_to_density(counts):
    """Convert {tp×region: n_sites} counts to sites per Mb of region length."""
    density = counts.copy().astype(float)
    density['core'] = counts['core'] / CORE_LEN_MB
    density['arm'] = counts['arm'] / ARM_LEN_MB
    return density


def draw_panel(ax, density, counts, title, motif, mod_label):
    """Grouped bar plot: T1/T2/T3 on x-axis; core vs arm bars per timepoint."""
    x = np.arange(len(TIMEPOINTS))
    width = 0.36

    bars_core = ax.bar(x - width / 2, density['core'].values, width,
                       color=COL_CORE, alpha=0.85, edgecolor='white',
                       linewidth=0.6, label=f'Core ({CORE_LEN_MB:.2f} Mb)')
    bars_arm = ax.bar(x + width / 2, density['arm'].values, width,
                      color=COL_ARM, alpha=0.85, edgecolor='white',
                      linewidth=0.6, label=f'Arm ({ARM_LEN_MB:.1f} Mb total)')

    # Per-bar count annotation (sites/Mb on top, raw n in parentheses)
    y_max = density.values.max()
    for region, bars in [('core', bars_core), ('arm', bars_arm)]:
        for i, bar in enumerate(bars):
            h = bar.get_height()
            n = counts.loc[TIMEPOINTS[i], region]
            if h > 0:
                ax.text(bar.get_x() + bar.get_width() / 2,
                        h + y_max * 0.015,
                        f'{h:.1f}\n(n={int(n)})',
                        ha='center', va='bottom', fontsize=6.5,
                        color='#333', linespacing=1.0)

    # Per-timepoint % core annotation (placed above the taller bar of the pair)
    for i, tp in enumerate(TIMEPOINTS):
        tot = counts.loc[tp].sum()
        pct_core = 100 * counts.loc[tp, 'core'] / tot if tot else 0
        pair_top = max(density.loc[tp, 'core'], density.loc[tp, 'arm'])
        ax.text(i, pair_top + y_max * 0.18,
                f'{pct_core:.0f}% core',
                ha='center', va='bottom', fontsize=7.5,
                color=COL_CORE if pct_core >= 50 else COL_ARM,
                fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                          edgecolor=COL_CORE if pct_core >= 50 else COL_ARM,
                          alpha=0.9, linewidth=0.6))

    ax.set_xticks(x)
    ax.set_xticklabels(TP_LABELS_NL, fontsize=8)
    ax.set_ylabel(f'{motif} {mod_label} sites per Mb', fontsize=9)
    ax.set_title(title, fontsize=10, fontweight='bold')
    ax.set_ylim(0, density.values.max() * 1.45)
    ax.grid(axis='y', alpha=0.3, lw=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def main():
    apply_style()
    print('=== Figure 1 panel: GCCGGC core->arm distribution ===')
    print()
    print(f'  Core length: {CORE_LEN_MB:.3f} Mb')
    print(f'  Arm length:  {ARM_LEN_MB:.3f} Mb (left + right)')
    print()

    print('Loading GCCGGC sites...')
    g_counts = load_gccggc_region_counts()
    g_density = counts_to_density(g_counts)
    print()
    print('Loading AAGCCCG sites...')
    a_counts = load_aagcccg_region_counts()
    a_density = counts_to_density(a_counts)
    print()

    # Two sub-panels (independent y-scales — densities differ by ~5x)
    fig, axes = plt.subplots(1, 2, figsize=(mm_to_inch(180), mm_to_inch(80)),
                             gridspec_kw={'wspace': 0.35})

    print('Drawing GCCGGC sub-panel...')
    draw_panel(axes[0], g_density, g_counts,
               'GCCGGC: core → arm redistribution',
               'GCCGGC', '(4mC)')
    axes[0].legend(fontsize=7, loc='upper right', frameon=False)

    # Highlight the T1→T2 inversion with a horizontal bracket above bars
    g_max = g_density.values.max()
    bracket_y = g_max * 1.32
    axes[0].annotate(
        '', xy=(0, bracket_y), xytext=(1, bracket_y),
        arrowprops=dict(arrowstyle='->', color='#B71C1C', lw=1.2))
    axes[0].text(0.5, bracket_y + g_max * 0.02,
                 'T1 → T2 inversion',
                 ha='center', va='bottom', fontsize=7.5,
                 color='#B71C1C', fontweight='bold')

    print('Drawing AAGCCCG sub-panel...')
    draw_panel(axes[1], a_density, a_counts,
               'AAGCCCG: stable core/arm ratio',
               'AAGCCCG', '(6mA)')
    axes[1].legend(fontsize=7, loc='upper right', frameon=False)

    # Annotate AAGCCCG as negative control
    a_max = a_density.values.max()
    axes[1].text(0.98, 0.50,
                 'No geographic shift\n(core ≈ 60–70% throughout)',
                 transform=axes[1].transAxes, fontsize=7,
                 color='#1565C0', ha='right', va='center',
                 fontstyle='italic',
                 bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                           edgecolor='#1565C0', alpha=0.85, linewidth=0.6))

    # Panel labels
    add_panel_label(axes[0], 'a', x=-0.18, y=1.10)
    add_panel_label(axes[1], 'b', x=-0.18, y=1.10)

    # Suptitle
    fig.suptitle(
        '4mC GCCGGC undergoes T1→T2 core-to-arm redistribution; '
        '6mA AAGCCCG distribution is timepoint-invariant',
        fontsize=10.5, y=1.02)

    out_main = FIG_DIR / 'Figure1_panel_GCCGGC_distribution'
    save_figure(fig, out_main, formats=('pdf', 'svg', 'png'))

    # Mirror PNG into Writing/fig_images for the manuscript
    png_src = out_main.with_suffix('.png')
    if png_src.exists():
        WRITING_FIG_DIR.mkdir(parents=True, exist_ok=True)
        dest = WRITING_FIG_DIR / 'Figure1_GCCGGC_arm_panel.png'
        dest.write_bytes(png_src.read_bytes())
        print(f'  Mirrored: {dest}')

    print()
    print('=== Done ===')


if __name__ == '__main__':
    main()
