#!/usr/bin/env python3
"""
FigS10: Per-motif Temporal Dynamics × Regional Stratification

Panel A: GCCGGC — T1→T2 transition groups × expression change
Panel B: GCCGGC — geographic stratification (core vs arm)
Panel C: AAGCCCG — T1→T2 transition groups × expression change
Panel D: AAGCCCG — geographic stratification (core vs arm)
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
from scipy import stats


def _style(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def _p_label(p):
    if p < 0.001:
        return '***'
    if p < 0.01:
        return '**'
    if p < 0.05:
        return '*'
    return 'n.s.'


# ── Panel A/C: Transition group × LFC violin ─────────────────────────────────

def draw_transition_violins(ax, df_trans, motif, transition_col, lfc_col, padj_col,
                            title, group_order, group_colors):
    """Violin plot of LFC by methylation transition group, overlaid with jitter."""
    df = df_trans[df_trans[transition_col].notna()].copy()
    df = df[df[transition_col].isin(group_order)]

    positions = np.arange(len(group_order))
    parts = ax.violinplot(
        [df[df[transition_col] == g][lfc_col].dropna().values for g in group_order],
        positions=positions,
        widths=0.7,
        showmedians=True,
        showextrema=False,
    )
    for i, (pc, g) in enumerate(zip(parts['bodies'], group_order)):
        pc.set_facecolor(group_colors[g])
        pc.set_alpha(0.6)
        pc.set_edgecolor('none')
    parts['cmedians'].set_color('#333333')
    parts['cmedians'].set_linewidth(1.5)

    # DEG fraction annotation
    for i, g in enumerate(group_order):
        sub = df[df[transition_col] == g]
        n_total = len(sub)
        n_sig = (sub[padj_col] < 0.05).sum() if padj_col in sub.columns else 0
        median = sub[lfc_col].median()
        frac = n_sig / n_total * 100 if n_total > 0 else 0
        ax.text(i, ax.get_ylim()[0] if ax.get_ylim()[0] > -5 else -4.5,
                f'n={n_total}', ha='center', va='bottom', fontsize=6, color='#555')

    ax.axhline(0, color='gray', linewidth=0.8, linestyle='--', alpha=0.5)
    ax.set_xticks(positions)
    ax.set_xticklabels([f'{g}\n({len(df[df[transition_col]==g]):,})'
                        for g in group_order], fontsize=7)
    ax.set_ylabel('log₂FC (vs T1)', fontsize=8)
    ax.set_title(title, fontsize=9, fontweight='bold')
    _style(ax)


# ── Panel B/D: Geographic stratified forest plot ──────────────────────────────

def draw_geo_forest(ax, df_geo, motif, color):
    """Forest plot of rank-biserial correlation by region × comparison."""
    df = df_geo.copy()

    # Create compound label
    df['label'] = df['region'].str.capitalize() + ': ' + df['comparison']

    # Sort: core first
    df['sort_key'] = df['region'].map({'core': 0, 'arm': 1, 'unknown': 2}).fillna(3)
    df = df.sort_values(['sort_key', 'comparison']).reset_index(drop=True)

    y = np.arange(len(df))
    for i, row in df.iterrows():
        col = color if row['p_value'] < 0.05 else '#BDBDBD'
        ax.scatter(row['rank_biserial'], i, color=col, s=40, zorder=3)

        # significance text
        lbl = _p_label(row['p_value'])
        ax.text(row['rank_biserial'] + 0.02, i, lbl, va='center', fontsize=7,
                color=col)

    ax.axvline(0, color='gray', linewidth=0.8, linestyle=':')
    ax.set_yticks(y)
    ax.set_yticklabels(df['label'], fontsize=7)
    ax.set_xlabel('Rank-biserial correlation', fontsize=8)
    ax.set_xlim(-0.35, 0.45)
    ax.set_title(f'{motif}: geographic stratification', fontsize=9, fontweight='bold')
    _style(ax)


# ── Panel for transition group DEG fraction ───────────────────────────────────

def draw_deg_bars(ax, df_expr, motif, color, title):
    """Bar chart of DEG fraction by transition group."""
    # df_expr has: group, n_genes, median_LFC, mean_LFC, n_DEG, frac_DEG, n_up, frac_up, n_down
    groups = df_expr['group'].values
    frac_up = df_expr['frac_up'].values * 100
    frac_down = df_expr['frac_down'].values * 100
    frac_ns = (1 - df_expr['frac_DEG'].values) * 100

    x = np.arange(len(groups))
    w = 0.6
    ax.bar(x, frac_up, width=w, color='#E53935', alpha=0.8, label='Up-regulated')
    ax.bar(x, frac_down, width=w, bottom=frac_up, color='#1E88E5', alpha=0.8,
           label='Down-regulated')
    ax.bar(x, frac_ns, width=w, bottom=frac_up + frac_down, color='#E0E0E0',
           alpha=0.8, label='Non-DEG')

    for i, (g, row) in enumerate(df_expr.iterrows()):
        ax.text(x[i], 102, f"n={df_expr.iloc[i]['n_genes']:,}", ha='center',
                fontsize=6, color='#444')

    ax.set_xticks(x)
    ax.set_xticklabels(groups, fontsize=7)
    ax.set_ylabel('% genes', fontsize=8)
    ax.set_ylim(0, 115)
    ax.set_title(title, fontsize=9, fontweight='bold')
    ax.legend(fontsize=6, loc='upper right', frameon=False)
    _style(ax)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    apply_style()
    print('=== FigS10: Per-motif Temporal × Regional ===')

    # Load data
    print('Loading GCCGGC data...')
    gccggc_trans, gccggc_geo = load_gene_methylation_transitions('GCCGGC')
    gccggc_expr = pd.read_csv(
        EPIGENOME / '42_GCCGGC_temporal_derepression' / 'tables' /
        'transition_group_expression.tsv', sep='\t')

    print('Loading AAGCCCG data...')
    aagcccg_trans, aagcccg_geo = load_gene_methylation_transitions('AAGCCCG')
    aagcccg_expr = pd.read_csv(
        EPIGENOME / '44_AAGCCCG_temporal_causality' / 'tables' /
        'transition_group_expression.tsv', sep='\t')

    # ── Figure layout ─────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(mm_to_inch(180), mm_to_inch(200)))
    outer = gridspec.GridSpec(2, 1, figure=fig, hspace=0.55,
                              left=0.09, right=0.97, top=0.96, bottom=0.06)

    # Row 1: GCCGGC
    gs_top = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[0],
                                              wspace=0.45, width_ratios=[1.4, 1.4, 1])
    ax_gv = fig.add_subplot(gs_top[0])
    ax_gd = fig.add_subplot(gs_top[1])
    ax_gg = fig.add_subplot(gs_top[2])

    # Row 2: AAGCCCG
    gs_bot = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[1],
                                              wspace=0.45, width_ratios=[1.4, 1.4, 1])
    ax_av = fig.add_subplot(gs_bot[0])
    ax_ad = fig.add_subplot(gs_bot[1])
    ax_ag = fig.add_subplot(gs_bot[2])

    # ── GCCGGC panels ────────────────────────────────────────────────────────
    print('Drawing GCCGGC panels...')
    add_panel_label(ax_gv, 'a', x=-0.12, y=1.08)

    # Group colors: Lost→blue (methylation lost=exposed), Gained→orange, Both→purple, Never→gray
    grp_colors_gcc = {
        'Lost': '#1E88E5', 'Gained': '#E53935',
        'Both': '#8E24AA', 'Never': '#9E9E9E',
    }
    grp_order_gcc = ['Lost', 'Gained', 'Both', 'Never']

    # Violin: T1→T2 LFC
    draw_transition_violins(
        ax_gv, gccggc_trans, 'GCCGGC', 'transition_T1T2', 'LFC_T2vsT1', 'padj_T2vsT1',
        'GCCGGC: T1→T2 methylation × expression',
        grp_order_gcc, grp_colors_gcc,
    )

    # DEG fraction bars: T2vsT1
    gccggc_expr_sorted = gccggc_expr.set_index('group').loc[
        [g for g in grp_order_gcc if g in gccggc_expr['group'].values]
    ].reset_index()
    draw_deg_bars(ax_gd, gccggc_expr_sorted, 'GCCGGC', COL_4mC,
                  'GCCGGC: DEG composition by group')

    # Geographic forest plot
    draw_geo_forest(ax_gg, gccggc_geo, 'GCCGGC', COL_4mC)

    # ── AAGCCCG panels ───────────────────────────────────────────────────────
    print('Drawing AAGCCCG panels...')
    add_panel_label(ax_av, 'b', x=-0.12, y=1.08)

    grp_colors_6ma = {
        'Lost': '#1E88E5', 'Gained': '#E53935',
        'Both': '#8E24AA', 'Never': '#9E9E9E',
    }
    grp_order_6ma = ['Lost', 'Gained', 'Both', 'Never']

    draw_transition_violins(
        ax_av, aagcccg_trans, 'AAGCCCG', 'transition', 'LFC_T2vsT1', 'padj_T2vsT1',
        'AAGCCCG: T1→T2 methylation × expression',
        grp_order_6ma, grp_colors_6ma,
    )

    aagcccg_expr_sorted = aagcccg_expr.set_index('group').loc[
        [g for g in grp_order_6ma if g in aagcccg_expr['group'].values]
    ].reset_index()
    draw_deg_bars(ax_ad, aagcccg_expr_sorted, 'AAGCCCG', COL_6mA,
                  'AAGCCCG: DEG composition by group')

    draw_geo_forest(ax_ag, aagcccg_geo, 'AAGCCCG', COL_6mA)

    # ── Shared legend for significance ───────────────────────────────────────
    sig_legend = [
        mpatches.Patch(color=COL_4mC, label='GCCGGC (4mC)'),
        mpatches.Patch(color=COL_6mA, label='AAGCCCG (6mA)'),
        mpatches.Patch(color='#BDBDBD', label='p ≥ 0.05 (n.s.)'),
    ]
    fig.legend(handles=sig_legend, loc='lower center', ncol=3, fontsize=7,
               frameon=False, bbox_to_anchor=(0.5, 0.0))

    # ── Save ─────────────────────────────────────────────────────────────────
    out_path = FIG_SUP_DIR / 'FigS10_per_motif_temporal_regional'
    save_figure(fig, out_path)
    print(f'\nSaved: {out_path}.pdf / .svg')
    print('=== Done ===')


if __name__ == '__main__':
    main()
