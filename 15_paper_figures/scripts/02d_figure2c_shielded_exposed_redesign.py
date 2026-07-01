#!/usr/bin/env python3
"""
Figure 2c redesign: Shielded/Exposed classification — threshold justification.

Single-panel layout:
  (A) Distance histogram — nearest GCCGGC site from TSS, log-scale x-axis
      • All shielded regulatory genes (gray fill)
      • Exposed TF positions (rug plot, purple)
      • 293 bp operating point line (red dashed)

NOTE (2026-07-01): the ROC/AUC panel was REMOVED. Per the locked paper
direction, the AUC of a distance-based classifier is circular (the 293 bp
threshold is itself defined from the distance distribution), so AUC is not
reported. The distance histogram alone justifies the 293 bp operating point
without the circular metric.

Data: all_genes_features_unified_n57.tsv  (is_exposed rebuilt to canonical 62)
  → 62 Exposed TFs + 989 Shielded regulatory genes  (n_total = 1,051)

Output:
  15_paper_figures/figures/main/Figure2c_shielded_exposed.png  (300 dpi)
  Writing/fig_images/Figure2c_shielded_exposed.png
"""

import importlib
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# ── Load shared utilities ───────────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).parent))
_utils = importlib.import_module('00_shared_utils')
for _attr in dir(_utils):
    if not _attr.startswith('_'):
        globals()[_attr] = getattr(_utils, _attr)

# ── Constants ───────────────────────────────────────────────────────────────
THRESHOLD_BP = 293
COL_SHIELDED_HIST = '#B0BEC5'   # gray (shielded)
COL_EXPOSED_RUG  = '#7E57C2'    # purple (exposed)
COL_THRESHOLD    = '#E53935'    # red (293 bp line / operating point)
COL_DIAGONAL     = '#90A4AE'    # gray (random classifier)


# ── Data loading ─────────────────────────────────────────────────────────────

def load_data():
    # Primary path (Mac native); fall back to bash sandbox mount
    path = (EPIGENOME / '52_shielded_exposed_boundary' /
            'tables' / 'all_genes_features_unified_n57.tsv')
    if not path.exists():
        path = Path('/sessions/clever-sweet-darwin/mnt/rna-seq/11_epigenome_integration'
                    '/analysis/52_shielded_exposed_boundary/tables'
                    '/all_genes_features_unified_n57.tsv')
    df = pd.read_csv(path, sep='\t')
    df = df.dropna(subset=['nearest_methyl_distance']).copy()
    n_exp = int(df['is_exposed'].sum())
    n_shi = int((df['is_exposed'] == 0).sum())
    print(f'  Loaded: {len(df)} genes  ({n_exp} Exposed / {n_shi} Shielded)')
    return df


# ── Panel A: Distance histogram ───────────────────────────────────────────────

def panel_a(ax, df):
    shielded_d = df.loc[df['is_exposed'] == 0, 'nearest_methyl_distance'].values
    exposed_d  = df.loc[df['is_exposed'] == 1, 'nearest_methyl_distance'].values

    # Log-spaced bins from 1 bp to 7000 bp
    bins = np.logspace(np.log10(1), np.log10(7000), 45)

    # Shielded histogram (gray, normalized to density)
    ax.hist(shielded_d, bins=bins, color=COL_SHIELDED_HIST, alpha=0.85,
            edgecolor='white', linewidth=0.3,
            label=f'Shielded (n={len(shielded_d)})',
            density=False, zorder=2)

    # Exposed histogram stacked on top (purple, small)
    ax.hist(exposed_d, bins=bins, color=COL_EXPOSED_RUG, alpha=0.85,
            edgecolor='white', linewidth=0.3,
            label=f'Exposed TF (n={len(exposed_d)})',
            density=False, zorder=3)

    # 293 bp threshold vertical line
    ax.axvline(THRESHOLD_BP, color=COL_THRESHOLD, linewidth=1.5,
               linestyle='--', zorder=5, label=f'{THRESHOLD_BP} bp threshold')

    # Annotation for threshold
    y_top = ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 300
    ax.text(THRESHOLD_BP * 1.07, ax.get_ylim()[1] * 0.92 if ax.get_ylim()[1] > 0 else 280,
            f'{THRESHOLD_BP} bp',
            color=COL_THRESHOLD, fontsize=8, va='top', fontweight='bold')

    # Rug plot for exposed TFs at bottom
    rug_y = -ax.get_ylim()[1] * 0.04 if ax.get_ylim()[1] > 0 else -8
    for d in exposed_d:
        ax.plot([d, d], [0, rug_y], color=COL_EXPOSED_RUG,
                alpha=0.7, linewidth=0.8, zorder=4, clip_on=False)

    ax.set_xscale('log')
    ax.set_xlim(3, 8000)

    # x-axis ticks: 10, 100, 1000 bp
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(
        lambda x, _: f'{int(x):,}'))
    ax.xaxis.set_major_locator(mticker.LogLocator(base=10, numticks=6))
    ax.xaxis.set_minor_locator(mticker.NullLocator())

    ax.set_xlabel('Nearest GCCGGC site to TSS (bp)', fontsize=9)
    ax.set_ylabel('Number of genes', fontsize=9)
    ax.set_title('Distance distribution', fontsize=9, fontweight='bold')

    leg = ax.legend(fontsize=7.5, loc='upper right', frameon=True,
                    fancybox=False, edgecolor='#ccc',
                    handlelength=1.2, handletextpad=0.5,
                    borderpad=0.5, labelspacing=0.3)

    # Spine cleanup
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def _fix_panel_a_rug(ax, df):
    """Re-draw rug after final y-limits are set."""
    exposed_d = df.loc[df['is_exposed'] == 1, 'nearest_methyl_distance'].values
    ymin, ymax = ax.get_ylim()
    rug_h = ymax * 0.03
    for d in exposed_d:
        ax.plot([d, d], [-rug_h, 0], color=COL_EXPOSED_RUG,
                alpha=0.7, linewidth=0.8, zorder=4, clip_on=False)
    # Re-apply threshold label at final height
    for t in list(ax.texts):
        t.remove()
    ax.text(THRESHOLD_BP * 1.07, ymax * 0.92,
            f'{THRESHOLD_BP} bp',
            color=COL_THRESHOLD, fontsize=8, va='top', fontweight='bold')


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    apply_style()
    print('=== Figure 2c redesign: Shielded/Exposed threshold justification ===')
    print()

    print('Loading data...')
    df = load_data()
    print()

    # Figure: single panel (ROC/AUC panel removed — AUC is circular, not reported),
    # 70 × 65 mm
    fig, ax_a = plt.subplots(
        1, 1,
        figsize=(mm_to_inch(70), mm_to_inch(65)),
    )
    fig.subplots_adjust(left=0.18, right=0.95, top=0.88, bottom=0.20)

    print('Drawing distance histogram...')
    panel_a(ax_a, df)
    # Fix rug after first draw sets y-limits
    _fix_panel_a_rug(ax_a, df)

    # ── Save ──────────────────────────────────────────────────────────────────
    out_dir = FIG_DIR
    if not out_dir.exists():
        # bash sandbox fallback
        out_dir = Path('/sessions/clever-sweet-darwin/mnt/rna-seq/15_paper_figures/figures/main')
    out_stem = out_dir / 'Figure2c_shielded_exposed'
    out_stem.parent.mkdir(parents=True, exist_ok=True)

    # PNG (primary deliverable — 300 dpi)
    png_path = out_stem.with_suffix('.png')
    fig.savefig(png_path, format='png', dpi=300, bbox_inches='tight')
    print(f'  Saved PNG: {png_path}')

    # PDF + SVG
    for fmt in ('pdf', 'svg'):
        p = out_stem.with_suffix(f'.{fmt}')
        fig.savefig(p, format=fmt, dpi=300, bbox_inches='tight')
        print(f'  Saved {fmt.upper()}: {p}')

    plt.close(fig)

    # Copy PNG to Writing/fig_images/
    writing_dir = BASE / 'Writing' / 'fig_images'
    if not (BASE / 'Writing').exists():
        writing_dir = Path('/sessions/clever-sweet-darwin/mnt/rna-seq/Writing/fig_images')
    writing_dir.mkdir(parents=True, exist_ok=True)
    dst = writing_dir / 'Figure2c_shielded_exposed.png'
    shutil.copy2(png_path, dst)
    print(f'  Copied to: {dst}')

    print()
    print('=== Done ===')


if __name__ == '__main__':
    main()
