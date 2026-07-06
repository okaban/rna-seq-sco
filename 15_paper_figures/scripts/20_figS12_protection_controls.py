#!/usr/bin/env python3
"""
FigS12: Protection Zone Controls — Regulatory vs Non-regulatory Specificity

Panel A: TSS methylation gradient — regulatory vs non-regulatory (all methylation types, T1)
Panel B: Protection zone by methylation type (width / minimum depletion ratio)
Panel C: TF subcategory methylation profiles near TSS
Panel D: T1 vs T2 temporal comparison — protection zone stability
"""

import sys, importlib
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
_utils = importlib.import_module('00_shared_utils')
for _attr in dir(_utils):
    if not _attr.startswith('_'):
        globals()[_attr] = getattr(_utils, _attr)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches

GRAD_DIR = EPIGENOME / '48_TSS_methylation_gradient' / 'tables'


def _style(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


# ── Panel A: Regulatory vs non-regulatory spatial profiles ───────────────────

def panel_a(ax, df_spatial):
    """Line plot — regulatory vs non-regulatory TSS methylation gradient."""
    site_types = [
        ('All_methylation', 'All methylation', '#333333', '-'),
        ('All_4mC',         'All 4mC',         COL_4mC,  '-'),
        ('All_6mA',         'All 6mA',         COL_6mA,  '-'),
    ]

    for site_type, label, color, ls in site_types:
        for cat, lw, alpha, suffix in [
            ('regulatory',     1.8, 0.9, ' (regulatory)'),
            ('non_regulatory', 1.0, 0.5, ' (non-reg)'),
        ]:
            sub = df_spatial[
                (df_spatial['site_type'] == site_type) &
                (df_spatial['timepoint'] == 'T1') &
                (df_spatial['gene_category'] == cat)
            ].sort_values('bin_center_bp')
            if len(sub) == 0:
                continue
            x = sub['bin_center_bp'].values / 1000
            y = sub['density_per_kb_per_gene'].values
            ax.plot(x, y, color=color, linewidth=lw, alpha=alpha,
                    linestyle=ls if cat == 'regulatory' else '--',
                    label=label + suffix if cat == 'regulatory' else None)

    # Shade protection zone region
    ax.axvspan(-0.3, 0.5, alpha=0.07, color='#1565C0', zorder=0)
    ax.axvline(0, color='gray', linewidth=0.8, linestyle=':')
    ax.set_xlabel('Distance from TSS (kb)', fontsize=8)
    ax.set_ylabel('Methylation density\n(sites / kb / gene)', fontsize=8)
    ax.set_title('TSS methylation gradient:\nregulatory vs non-regulatory', fontsize=9,
                 fontweight='bold', linespacing=1.3)
    ax.legend(fontsize=6.5, loc='upper right', frameon=True, fancybox=False,
              edgecolor='#ccc', ncol=1)
    _style(ax)

    # Annotation: protection zone
    ax.text(-0.3, ax.get_ylim()[1] * 0.95, 'Protection\nzone', fontsize=6.5,
            color='#1565C0', style='italic', va='top')


# ── Panel B: Protection zone metrics ─────────────────────────────────────────

def panel_b(ax_w, ax_r, df_metrics):
    """Bar charts: protection zone width and minimum depletion ratio by site type."""
    # Exclude AAGCCCG_6mA (no zone)
    df = df_metrics[df_metrics['site_type'] != 'AAGCCCG_6mA'].copy()
    labels_map = {
        'All_methylation': 'All\nmethyl',
        'All_4mC':         'All\n4mC',
        'All_6mA':         'All\n6mA',
        'GCCGGC_4mC':      'GCCGGC\n4mC',
    }
    color_map = {
        'All_methylation': '#333333',
        'All_4mC':         COL_4mC,
        'All_6mA':         COL_6mA,
        'GCCGGC_4mC':      COL_4mC,
    }
    labels = [labels_map.get(s, s) for s in df['site_type']]
    colors = [color_map.get(s, '#90A4AE') for s in df['site_type']]
    x = np.arange(len(df))

    # Width
    ax_w.bar(x, df['protection_zone_width_bp'].values, color=colors, alpha=0.8,
             edgecolor='white', linewidth=0.5, width=0.6)
    for xi, v in zip(x, df['protection_zone_width_bp'].values):
        ax_w.text(xi, v + 10, f'{int(v)}', ha='center', va='bottom', fontsize=7)
    ax_w.set_xticks(x)
    ax_w.set_xticklabels(labels, fontsize=7)
    ax_w.set_ylabel('Protection zone width (bp)', fontsize=8)
    ax_w.set_title('Zone width', fontsize=8, fontweight='bold')
    _style(ax_w)

    # Min ratio
    ax_r.bar(x, df['min_ratio_reg_vs_nonreg'].values, color=colors, alpha=0.8,
             edgecolor='white', linewidth=0.5, width=0.6)
    ax_r.axhline(1.0, color='gray', linewidth=0.8, linestyle='--')
    for xi, v in zip(x, df['min_ratio_reg_vs_nonreg'].values):
        ax_r.text(xi, v + 0.01, f'{v:.2f}', ha='center', va='bottom', fontsize=7)
    ax_r.set_xticks(x)
    ax_r.set_xticklabels(labels, fontsize=7)
    ax_r.set_ylabel('Min ratio (reg / non-reg)', fontsize=8)
    ax_r.set_title('Depletion magnitude', fontsize=8, fontweight='bold')
    _style(ax_r)


# ── Panel C: TF subcategory profiles ─────────────────────────────────────────

def panel_c(ax, df_subcat):
    """TSS methylation profiles by TF subcategory."""
    subtype_config = [
        ('sigma_factor',    'σ factor (n=36)',          '#E53935', '-',  1.6),
        ('TCS_response_reg','TCS response reg (n=29)',  '#1E88E5', '-',  1.6),
        ('sensor_kinase',   'Sensor kinase (n=16)',     '#43A047', '-',  1.6),
        ('other_TF',        'Other TF (n=365)',         '#9E9E9E', '--', 1.0),
        ('non_regulatory',  'Non-regulatory (n=2196)',  '#CFD8DC', ':',  0.8),
    ]

    for subtype, label, color, ls, lw in subtype_config:
        sub = df_subcat[df_subcat['reg_subtype'] == subtype].sort_values('bin_center_bp')
        if len(sub) == 0:
            continue
        x = sub['bin_center_bp'].values / 1000
        y = sub['density_per_kb_per_gene'].values
        ci_lo = sub.get('ci_low', pd.Series([np.nan] * len(sub))).values
        ci_hi = sub.get('ci_high', pd.Series([np.nan] * len(sub))).values

        ax.plot(x, y, color=color, linewidth=lw, linestyle=ls, label=label)
        if not np.all(np.isnan(ci_lo)):
            ax.fill_between(x, ci_lo, ci_hi, color=color, alpha=0.10)

    ax.axvline(0, color='gray', linewidth=0.8, linestyle=':')
    ax.axvspan(-0.3, 0.5, alpha=0.05, color='#1565C0', zorder=0)
    ax.set_xlabel('Distance from TSS (kb)', fontsize=8)
    ax.set_ylabel('Methylation density\n(sites / kb / gene)', fontsize=8)
    ax.set_title('Protection zone by TF subcategory\n(All methylation, T1)', fontsize=9,
                 fontweight='bold', linespacing=1.3)
    ax.legend(fontsize=6.5, loc='upper right', frameon=True, fancybox=False,
              edgecolor='#ccc')
    _style(ax)


# ── Panel D: T1 vs T2 stability ──────────────────────────────────────────────

def panel_d(ax, df_spatial):
    """Compare regulatory TSS profiles T1 vs T2 — test if protection is stable."""
    for tp, color, label in [('T1', COL_EXPOSED, 'T1 (growth)'),
                               ('T2', '#FB8C00',  'T2 (transition)')]:
        sub = df_spatial[
            (df_spatial['site_type'] == 'All_methylation') &
            (df_spatial['timepoint'] == tp) &
            (df_spatial['gene_category'] == 'regulatory')
        ].sort_values('bin_center_bp')
        if len(sub) == 0:
            continue
        x = sub['bin_center_bp'].values / 1000
        y = sub['density_per_kb_per_gene'].values
        ax.plot(x, y, color=color, linewidth=1.8, label=label)

    ax.axvline(0, color='gray', linewidth=0.8, linestyle=':')
    ax.axvspan(-0.3, 0.5, alpha=0.07, color='#1565C0', zorder=0)
    ax.set_xlabel('Distance from TSS (kb)', fontsize=8)
    ax.set_ylabel('Methylation density\n(sites / kb / gene)', fontsize=8)
    ax.set_title('Protection zone temporal stability\n(regulatory genes)', fontsize=9,
                 fontweight='bold', linespacing=1.3)
    ax.legend(fontsize=7, frameon=False)
    _style(ax)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    apply_style()
    print('=== FigS12: Protection Zone Controls ===')

    df_spatial = pd.read_csv(GRAD_DIR / 'spatial_profile_data.tsv', sep='\t')
    df_metrics = pd.read_csv(GRAD_DIR / 'protection_zone_metrics.tsv', sep='\t')
    df_subcat = pd.read_csv(GRAD_DIR / 'subcategory_profiles.tsv', sep='\t')
    print(f'Spatial: {len(df_spatial)} rows')
    print(f'Metrics: {len(df_metrics)} site types')
    print(f'Subcategories: {df_subcat["reg_subtype"].unique().tolist()}')

    # ── Figure layout ─────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(mm_to_inch(174), mm_to_inch(190)))
    gs = gridspec.GridSpec(2, 4, figure=fig,
                           hspace=0.50, wspace=0.50,
                           left=0.09, right=0.97, top=0.96, bottom=0.07)

    ax_a = fig.add_subplot(gs[0, 0:2])
    ax_bw = fig.add_subplot(gs[0, 2])
    ax_br = fig.add_subplot(gs[0, 3])
    ax_c = fig.add_subplot(gs[1, 0:2])
    ax_d = fig.add_subplot(gs[1, 2:4])

    add_panel_label(ax_a, 'a', x=-0.10, y=1.08)
    add_panel_label(ax_bw, 'b', x=-0.18, y=1.08)
    add_panel_label(ax_c, 'c', x=-0.10, y=1.08)
    add_panel_label(ax_d, 'd', x=-0.10, y=1.08)

    print('Drawing Panel A: Regulatory vs non-regulatory profiles...')
    panel_a(ax_a, df_spatial)

    print('Drawing Panel B: Protection zone metrics...')
    panel_b(ax_bw, ax_br, df_metrics)

    print('Drawing Panel C: TF subcategory profiles...')
    panel_c(ax_c, df_subcat)

    print('Drawing Panel D: T1 vs T2 stability...')
    panel_d(ax_d, df_spatial)

    # ── Save ─────────────────────────────────────────────────────────────────
    out_path = FIG_SUP_DIR / 'FigS12_protection_controls'
    save_figure(fig, out_path)
    print(f'\nSaved: {out_path}.pdf / .svg')
    print('=== Done ===')


if __name__ == '__main__':
    main()
