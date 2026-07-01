#!/usr/bin/env python3
"""
Figure 5 (exposed TFs, n=57 update) — standalone runner.
Panels: (a) Co-expression heatmap 57×57
        (b) Temporal trajectories (bloc)
        (c) Early vs Late scatter
        (d) TF family by functional category
        (e) TCS pair asymmetry
Output: figures/main/Figure5_exposed_TFs.{pdf,svg,png}
"""
import sys
import warnings
warnings.filterwarnings('ignore')

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from scipy import stats

# ── Paths ────────────────────────────────────────────────────────────────────
BASE     = Path(__file__).resolve().parents[2]          # …/rna-seq
EPIGENOME = BASE / '11_epigenome_integration' / 'analysis'
FIG_DIR  = BASE / '15_paper_figures' / 'figures' / 'main'

# ── Style constants ───────────────────────────────────────────────────────────
COL_ACTIVATION = '#43A047'
COL_REPRESSION = '#FB8C00'
COL_EXPOSED    = '#7E57C2'
COL_SHIELDED   = '#B0BEC5'
COL_GRAY       = '#B0BEC5'
COL_DARK       = '#37474F'

STYLE = {
    'font.family':      'DejaVu Sans',
    'font.size':        10,
    'axes.titlesize':   11,
    'axes.labelsize':   10,
    'xtick.labelsize':  9,
    'ytick.labelsize':  9,
    'legend.fontsize':  8,
    'figure.dpi':       300,
    'savefig.dpi':      300,
    'axes.linewidth':   1.0,
    'axes.spines.top':  False,
    'axes.spines.right': False,
    'pdf.fonttype':     42,
    'ps.fonttype':      42,
    'svg.fonttype':     'none',
}
plt.rcParams.update(STYLE)


def mm(x):
    return x / 25.4


def add_label(ax, label, x=-0.12, y=1.08):
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=13, fontweight='bold', va='top', ha='left')


# ── Data loaders ──────────────────────────────────────────────────────────────

def load_coexpr():
    p = EPIGENOME / '55_exposed_regulatory_module' / 'tables' / 'coexpression_matrix.tsv'
    df = pd.read_csv(p, sep='\t', index_col=0)
    print(f'  Co-expression matrix: {df.shape}')
    return df


def load_temporal():
    p = EPIGENOME / '57_temporal_dynamics_exposed_TF' / 'tables' / 'temporal_classification.tsv'
    df = pd.read_csv(p, sep='\t')
    print(f'  Temporal classification: {len(df)} TFs  '
          f'(act={sum(df.bloc=="activation")}, '
          f'rep={sum(df.bloc=="repression")}, '
          f'una={sum(df.bloc=="unassigned")})')
    return df


def load_bloc_comp():
    p = EPIGENOME / '57_temporal_dynamics_exposed_TF' / 'tables' / 'bloc_comparison.tsv'
    return pd.read_csv(p, sep='\t')


def load_family():
    p = EPIGENOME / '58_exposed_TF_functional_prediction' / 'tables' / 'bloc_family_enrichment.tsv'
    df = pd.read_csv(p, sep='\t').dropna(subset=['tf_family'])
    print(f'  Family enrichment: {len(df)} families')
    return df


def load_tcs():
    p = EPIGENOME / '55_exposed_regulatory_module' / 'tables' / 'TCS_pairs_analysis.tsv'
    df = pd.read_csv(p, sep='\t')
    print(f'  TCS pairs: {len(df)}')
    return df


# ── Panels ────────────────────────────────────────────────────────────────────

def panel_a(ax, df_coexpr, df_temporal):
    """Co-expression heatmap, 57×57, ordered by direction then module."""
    # Use 'direction' (up/down) if repression bloc is absent
    has_repression = (df_temporal['bloc'] == 'repression').any()
    if has_repression:
        bloc_order = {'activation': 0, 'repression': 1, 'unassigned': 2}
        sort_key   = 'bloc'
    else:
        # Fall back to direction: down-regulated TFs map to "repression-like" group
        df_temporal = df_temporal.copy()
        df_temporal['display_bloc'] = df_temporal.apply(
            lambda r: 'repression' if r['direction'] == 'down' else r['bloc'], axis=1)
        bloc_order = {'activation': 0, 'repression': 1, 'unassigned': 2}
        sort_key   = 'display_bloc'

    gene_order = (df_temporal
                  .assign(bloc_rank=df_temporal[sort_key].map(bloc_order))
                  .sort_values(['bloc_rank', 'module', 'locus_tag'])
                  ['locus_tag'].tolist())
    gene_order = [g for g in gene_order if g in df_coexpr.columns]

    mat = df_coexpr.loc[gene_order, gene_order].values

    n_act = int((df_temporal.set_index('locus_tag')
                 .loc[[g for g in gene_order if g in df_temporal['locus_tag'].values]]
                 [sort_key] == 'activation').sum())
    n_rep = int((df_temporal.set_index('locus_tag')
                 .loc[[g for g in gene_order if g in df_temporal['locus_tag'].values]]
                 [sort_key] == 'repression').sum())

    im = ax.imshow(mat, cmap='RdBu_r', vmin=-1, vmax=1, aspect='equal')

    if n_act > 0:
        ax.axhline(n_act - 0.5, color='black', linewidth=1.0)
        ax.axvline(n_act - 0.5, color='black', linewidth=1.0)
        ax.text(-3, n_act / 2, 'Up', ha='right', va='center', fontsize=7,
                color=COL_ACTIVATION, fontweight='bold', rotation=90)
    if n_rep > 0:
        ax.text(-3, n_act + n_rep / 2, 'Down', ha='right', va='center', fontsize=7,
                color=COL_REPRESSION, fontweight='bold', rotation=90)
        if n_act + n_rep < len(gene_order):
            ax.axhline(n_act + n_rep - 0.5, color='black', linewidth=0.5, linestyle='--')
            ax.axvline(n_act + n_rep - 0.5, color='black', linewidth=0.5, linestyle='--')

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('Co-expression matrix\n(57 exposed TFs)', fontsize=9, fontweight='bold')

    cbar = plt.colorbar(im, ax=ax, shrink=0.7, pad=0.03, aspect=15)
    cbar.set_label('Spearman ρ', fontsize=7)
    cbar.ax.tick_params(labelsize=6)


def panel_b(ax, df_temporal):
    """Temporal trajectories grouped by direction (up/down)."""
    z_cols = ['T1_z', 'T2_z', 'T3_z']
    x = [1, 2, 3]

    for direction, color, marker in [
        ('up',   COL_ACTIVATION, 'o'),
        ('down', COL_REPRESSION, 's'),
    ]:
        sub = df_temporal[df_temporal['direction'] == direction]
        if len(sub) == 0:
            continue
        for _, row in sub.iterrows():
            ax.plot(x, [row[c] for c in z_cols],
                    color=color, alpha=0.08, linewidth=0.5)
        means = [sub[c].mean() for c in z_cols]
        sems  = [sub[c].std() / np.sqrt(len(sub)) for c in z_cols]
        lbl   = f'Up (n={len(sub)})' if direction == 'up' else f'Down (n={len(sub)})'
        ax.plot(x, means, color=color, linewidth=2.5, marker=marker,
                markersize=6, label=lbl, zorder=5)
        ax.fill_between(x,
                        [m - s for m, s in zip(means, sems)],
                        [m + s for m, s in zip(means, sems)],
                        color=color, alpha=0.2)

    ax.axhline(0, color='gray', linewidth=0.5, linestyle=':')
    ax.set_xticks(x)
    ax.set_xticklabels(['T1', 'T2', 'T3'])
    ax.set_xlabel('Timepoint', fontsize=9)
    ax.set_ylabel('Expression (z-score)', fontsize=9)
    ax.set_title('Temporal trajectories\n(antagonistic groups)', fontsize=9, fontweight='bold')
    ax.legend(fontsize=7, loc='upper left', frameon=False)


def panel_c(ax, df_temporal):
    """Early vs Late scatter: |LFC_T2vsT1| vs |LFC_T3vsT2|."""
    df = df_temporal.copy()
    df['abs_t1t2'] = df['LFC_T2vsT1'].abs()
    df['abs_t2t3'] = df['LFC_T3vsT2'].abs()

    colors = {'up': COL_ACTIVATION, 'down': COL_REPRESSION}
    markers = {'up': 'o', 'down': 's'}

    for direction in ['up', 'down']:
        sub = df[df['direction'] == direction]
        if len(sub) == 0:
            continue
        label = f'Up (n={len(sub)})' if direction == 'up' else f'Down (n={len(sub)})'
        ax.scatter(sub['abs_t1t2'], sub['abs_t2t3'],
                   color=colors[direction], marker=markers[direction],
                   s=28, alpha=0.7, edgecolors='white', linewidths=0.4,
                   label=label, zorder=4)

    x_max = df['abs_t1t2'].max() * 1.1 if len(df) else 5
    ax.plot([0, x_max], [0, x_max * (2/3)], color='gray', linewidth=0.9,
            linestyle='--', zorder=2, label='Phase ratio = 0.6')
    ax.text(0.03, 0.70, 'Late-\ndominant', transform=ax.transAxes,
            fontsize=6, color='gray', va='center', style='italic')
    ax.text(0.75, 0.12, 'Early-\ndominant', transform=ax.transAxes,
            fontsize=6, color='gray', va='bottom', style='italic')
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel('T1→T2 |log2FC|', fontsize=8)
    ax.set_ylabel('T2→T3 |log2FC|', fontsize=8)
    n_early = (df['phase_ratio'] > 0.6).sum() if 'phase_ratio' in df.columns else 0
    ax.set_title(f'Response timing\n({n_early}/{len(df)} early-dominant)',
                 fontsize=9, fontweight='bold')
    ax.legend(fontsize=6, loc='upper right', frameon=False)


def panel_d(ax, df_family):
    """TF family by functional category."""
    FUNC_COLORS = {
        'Signal transduction': '#5C6BC0',
        'Stress/development':  '#26A69A',
        'Defense/resistance':  '#EF5350',
        'Metabolic':           '#8D6E63',
        'General':             '#78909C',
    }
    FUNC_ORDER = ['Signal transduction', 'Stress/development',
                  'Defense/resistance', 'Metabolic', 'General']

    df = df_family.copy()
    df['total'] = df['activation_count'] + df['repression_count']
    df = df[df['total'] >= 1].copy()
    df['func_rank'] = df['functional_category'].map(
        {c: i for i, c in enumerate(FUNC_ORDER)}).fillna(99)
    df = df.sort_values(['func_rank', 'total'], ascending=[True, True]).reset_index(drop=True)

    gap_offset = 0
    y_positions = []
    prev_cat = None
    for i, row in df.iterrows():
        fcat = row.get('functional_category', 'General')
        if fcat != prev_cat and prev_cat is not None:
            gap_offset += 0.7
        y_positions.append(i + gap_offset)
        prev_cat = fcat

    y_arr = np.array(y_positions)
    bar_h = 0.32

    for i, row in df.iterrows():
        fcat = row.get('functional_category', 'General')
        fc = FUNC_COLORS.get(fcat, '#78909C')
        y = y_arr[i]
        ax.barh(y - bar_h/2, row['activation_count'], bar_h, color=fc, alpha=0.85)
        ax.barh(y + bar_h/2, row['repression_count'], bar_h, color=fc, alpha=0.35,
                hatch='////', edgecolor=fc)

    ax.set_yticks(y_arr)
    ax.set_yticklabels(df['tf_family'], fontsize=7)
    ax.set_xlabel('Number of TFs', fontsize=9)
    ax.set_title('TF family by functional category', fontsize=9, fontweight='bold')
    legend_els = [
        Patch(fc='gray', alpha=0.85, label='Up-regulated'),
        Patch(fc='gray', alpha=0.35, hatch='////', ec='gray', label='Down-regulated'),
    ]
    ax.legend(handles=legend_els, fontsize=7, loc='lower right', frameon=False)


def panel_e(ax, df_tcs):
    """TCS pair asymmetry — 7 identified pairs."""
    n_pairs = len(df_tcs)

    for i, (_, row) in enumerate(df_tcs.iterrows()):
        y = n_pairs - 1 - i
        sk_col = COL_EXPOSED if row['sk_exposed'] else COL_SHIELDED
        rr_col = COL_EXPOSED if row['rr_exposed'] else COL_SHIELDED

        ax.scatter(0.30, y, color=sk_col, s=80, zorder=5,
                   edgecolors='black', linewidths=0.5, marker='D')
        ax.text(0.05, y, str(row.get('sk_old_locus', '')),
                ha='left', va='center', fontsize=6)
        ax.scatter(0.70, y, color=rr_col, s=80, zorder=5,
                   edgecolors='black', linewidths=0.5, marker='o')
        ax.text(0.95, y, str(row.get('rr_old_locus', '')),
                ha='right', va='center', fontsize=6)
        ax.plot([0.34, 0.66], [y, y], color='#999', linewidth=0.8)

    ax.text(0.30, n_pairs - 0.2, 'SK', ha='center', va='bottom',
            fontsize=8, fontweight='bold')
    ax.text(0.70, n_pairs - 0.2, 'RR', ha='center', va='bottom',
            fontsize=8, fontweight='bold')

    legend_els = [
        Line2D([0], [0], marker='D', color='w', markerfacecolor=COL_EXPOSED,
               markersize=8, label='Exposed', markeredgecolor='black', markeredgewidth=0.5),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=COL_SHIELDED,
               markersize=8, label='Shielded', markeredgecolor='black', markeredgewidth=0.5),
    ]
    ax.legend(handles=legend_els, fontsize=7, loc='lower center',
              frameon=False, ncol=2)
    ax.text(0.50, n_pairs + 0.15,
            '7/7 identified pairs: exactly one partner exposed',
            ha='center', va='bottom', fontsize=6, fontstyle='italic', color=COL_DARK)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-1.2, n_pairs + 1.2)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_title('TCS pair asymmetry', fontsize=9, fontweight='bold')


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print('=== Figure 5: 57 Exposed TFs (n=57 update) ===')
    print()
    print('Loading data...')
    df_coexpr   = load_coexpr()
    df_temporal = load_temporal()
    df_bloc     = load_bloc_comp()
    df_family   = load_family()
    df_tcs      = load_tcs()
    print()

    fig = plt.figure(figsize=(mm(180), mm(210)))
    gs = fig.add_gridspec(2, 12, hspace=0.50, wspace=1.2,
                          left=0.08, right=0.97, top=0.95, bottom=0.05,
                          height_ratios=[1, 1])
    ax_a = fig.add_subplot(gs[0, 0:5])
    ax_b = fig.add_subplot(gs[0, 5:9])
    ax_c = fig.add_subplot(gs[0, 9:12])
    ax_d = fig.add_subplot(gs[1, 0:6])
    ax_e = fig.add_subplot(gs[1, 6:12])

    print('Panel A: Co-expression heatmap...')
    panel_a(ax_a, df_coexpr, df_temporal)
    add_label(ax_a, 'a', x=-0.12, y=1.08)

    print('Panel B: Temporal trajectories...')
    panel_b(ax_b, df_temporal)
    add_label(ax_b, 'b', x=-0.18, y=1.08)

    print('Panel C: Early vs Late scatter...')
    panel_c(ax_c, df_temporal)
    add_label(ax_c, 'c', x=-0.20, y=1.08)

    print('Panel D: TF family...')
    panel_d(ax_d, df_family)
    add_label(ax_d, 'd', x=-0.10, y=1.08)

    print('Panel E: TCS pairs...')
    panel_e(ax_e, df_tcs)
    add_label(ax_e, 'e', x=-0.06, y=1.08)

    out = FIG_DIR / 'Figure5_exposed_TFs'
    for fmt in ('pdf', 'svg', 'png'):
        fpath = out.with_suffix(f'.{fmt}')
        fig.savefig(fpath, format=fmt, dpi=300, bbox_inches='tight')
        print(f'  Saved: {fpath}')
    plt.close(fig)
    print()
    print('=== Done ===')


if __name__ == '__main__':
    main()
