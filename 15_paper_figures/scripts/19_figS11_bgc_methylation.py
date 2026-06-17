#!/usr/bin/env python3
"""
FigS11: BGC × Methylation

Panel A: BGC gene GCCGGC/AAGCCCG enrichment — fold enrichment vs non-BGC
         (all genome and core-only subsets)
Panel B: Cluster-specific GCCGGC methylation density
Panel C: DNA sequence content vs actual methylation density in BGC vs non-BGC
Panel D: GCCGGC fold enrichment per cluster vs core baseline
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

BGC_DIR = EPIGENOME / '47_BGC_methylation_geographic_test' / 'tables'

CLUSTER_COLORS = {  # muted (Okabe-Ito/Tol) for palette consistency with main figures
    'act': '#C26B6B',  # muted rose — actinorhodin
    'cda': '#4477AA',  # muted blue — calcium-dependent antibiotic
    'cpk': '#009E73',  # muted green — coelimycin
    'red': '#E69F00',  # muted orange — undecylprodigiosin
}

def _style(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def _p_label(p):
    if p < 0.001: return '***'
    if p < 0.01:  return '**'
    if p < 0.05:  return '*'
    return 'n.s.'


# ── Panel A: Enrichment fold change ──────────────────────────────────────────

def panel_a(ax, df_enrich):
    """Grouped bar chart: fold enrichment in BGC vs non-BGC, by motif and subset."""
    # Focus on GCCGGC and AAGCCCG, all_genes subset only
    motifs = ['GCCGGC', 'AAGCCCG', 'All4mC']
    labels = ['GCCGGC\n(4mC)', 'AAGCCCG\n(6mA)', 'All 4mC']
    colors = [COL_4mC, COL_6mA, '#90A4AE']

    all_sub = df_enrich[df_enrich['subset'] == 'all_genes']
    core_sub = df_enrich[df_enrich['subset'] == 'core_only']

    x = np.arange(len(motifs))
    w = 0.35

    for i, (motif, label, color) in enumerate(zip(motifs, labels, colors)):
        ar = all_sub[all_sub['motif'] == motif].iloc[0]
        cr = core_sub[core_sub['motif'] == motif].iloc[0]

        # All-genome bar
        ax.bar(x[i] - w/2, ar['fold_enrichment'], width=w, color=color,
               alpha=0.85, edgecolor='white', linewidth=0.5, label='All genome' if i == 0 else '')
        # Core bar
        ax.bar(x[i] + w/2, cr['fold_enrichment'], width=w, color=color,
               alpha=0.45, edgecolor='white', linewidth=0.5,
               hatch='///', label='Core only' if i == 0 else '')

        # Significance annotation
        for j, (row, offset) in enumerate([(ar, -w/2), (cr, +w/2)]):
            lbl = _p_label(row['pvalue_bonferroni'])
            ypos = max(row['fold_enrichment'], 1.0) + 0.05
            ax.text(x[i] + offset, ypos, lbl, ha='center', va='bottom',
                    fontsize=7, fontweight='bold',
                    color='#C62828' if lbl != 'n.s.' else '#9E9E9E')

    ax.axhline(1.0, color='gray', linewidth=0.8, linestyle='--')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel('Fold enrichment (BGC vs non-BGC)', fontsize=8)
    ax.set_title('Methylation enrichment in BGC genes', fontsize=9, fontweight='bold')
    ax.legend(fontsize=7, loc='upper right', frameon=False)
    _style(ax)


# ── Panel B: Per-cluster density ──────────────────────────────────────────────

def panel_b(ax, df_cluster):
    """Bar chart of per-cluster GCCGGC methylation density."""
    clusters = df_cluster['cluster'].values
    density = df_cluster['GCCGGC_density_per_kb'].values
    colors = [CLUSTER_COLORS.get(c, '#BDBDBD') for c in clusters]

    x = np.arange(len(clusters))
    ax.bar(x, density, color=colors, alpha=0.85, edgecolor='white', linewidth=0.5,
           width=0.6)

    # Core genome baseline
    core_baseline = 0.2333  # mean GCCGGC density in core non-BGC genes (from DNA_vs_methylation)
    ax.axhline(core_baseline, color='#333', linewidth=1.0, linestyle='--',
               label=f'Core genome baseline ({core_baseline:.3f})')

    # Value labels
    for xi, d in zip(x, density):
        ax.text(xi, d + 0.005, f'{d:.3f}', ha='center', va='bottom', fontsize=7)

    ax.set_xticks(x)
    ax.set_xticklabels([c.upper() for c in clusters], fontsize=8)
    ax.set_ylabel('GCCGGC density (sites / kb)', fontsize=8)
    ax.set_title('Per-cluster GCCGGC methylation density', fontsize=9, fontweight='bold')
    ax.legend(fontsize=7, loc='upper right', frameon=False)
    _style(ax)


# ── Panel C: DNA sequence vs methylation ──────────────────────────────────────

def panel_c(ax, df_dna):
    """Compare DNA GCCGGC content vs actual methylation density, BGC vs non-BGC."""
    # Rows of interest: DNA_TGGCCGGC_extended, Methylation_GCCGGC_density
    dna_row = df_dna[df_dna['comparison'] == 'DNA_TGGCCGGC_body'].iloc[0]
    meth_row = df_dna[df_dna['comparison'] == 'Methylation_GCCGGC_density'].iloc[0]

    categories = ['DNA\n(TGGCCGGC seq)', 'Methylation\n(GCCGGC 4mC)']
    bgc_vals = [dna_row['BGC_core_mean'], meth_row['BGC_core_mean']]
    nonbgc_vals = [dna_row['nonBGC_core_mean'], meth_row['nonBGC_core_mean']]
    pvals = [dna_row['MWU_pvalue'], meth_row['MWU_pvalue']]
    fold = [dna_row['fold'], meth_row['fold']]

    x = np.arange(len(categories))
    w = 0.35
    ax.bar(x - w/2, bgc_vals, width=w, color=COL_4mC, alpha=0.8,
           edgecolor='white', linewidth=0.5, label='BGC genes')
    ax.bar(x + w/2, nonbgc_vals, width=w, color='#90A4AE', alpha=0.8,
           edgecolor='white', linewidth=0.5, label='Non-BGC genes')

    for i, (p, f) in enumerate(zip(pvals, fold)):
        ymax = max(bgc_vals[i], nonbgc_vals[i])
        ax.text(x[i], ymax + 0.01, _p_label(p),
                ha='center', va='bottom', fontsize=8, fontweight='bold',
                color='#C62828' if _p_label(p) != 'n.s.' else '#9E9E9E')
        ax.text(x[i], ymax + 0.03, f'×{f:.2f}',
                ha='center', va='bottom', fontsize=6.5, color='#444')

    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=8)
    ax.set_ylabel('Mean density (proportion)', fontsize=8)
    ax.set_title('DNA sequence vs methylation in BGC\n(core genes only)', fontsize=9,
                 fontweight='bold', linespacing=1.3)
    ax.legend(fontsize=7, frameon=False)
    _style(ax)


# ── Panel D: Per-cluster fold vs core baseline ────────────────────────────────

def panel_d(ax, df_cluster):
    """Show fold enrichment per cluster vs core genome baseline."""
    clusters = df_cluster['cluster'].values
    fold_gcc = df_cluster['GCCGGC_fold_vs_core'].values
    fold_6ma = df_cluster['AAGCCCG_fold_vs_core'].values

    x = np.arange(len(clusters))
    w = 0.35

    ax.bar(x - w/2, fold_gcc, width=w, color=COL_4mC, alpha=0.8,
           edgecolor='white', linewidth=0.5, label='GCCGGC (4mC)')
    ax.bar(x + w/2, fold_6ma, width=w, color=COL_6mA, alpha=0.8,
           edgecolor='white', linewidth=0.5, label='AAGCCCG (6mA)')

    ax.axhline(1.0, color='gray', linewidth=0.8, linestyle='--',
               label='Core genome baseline')

    ax.set_xticks(x)
    ax.set_xticklabels([c.upper() for c in clusters], fontsize=8)
    ax.set_ylabel('Fold enrichment vs core baseline', fontsize=8)
    ax.set_title('BGC-specific enrichment per cluster', fontsize=9, fontweight='bold')
    ax.legend(fontsize=7, frameon=False)
    _style(ax)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    apply_style()
    print('=== FigS11: BGC × Methylation ===')

    df_enrich = pd.read_csv(BGC_DIR / 'enrichment_tests.tsv', sep='\t')
    df_cluster = pd.read_csv(BGC_DIR / 'cluster_specific_methylation.tsv', sep='\t')
    df_dna = pd.read_csv(BGC_DIR / 'DNA_vs_methylation_BGC.tsv', sep='\t')
    print(f'Enrichment tests: {len(df_enrich)} rows')
    print(f'Clusters: {len(df_cluster)}')

    # ── Figure layout ─────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(mm_to_inch(180), mm_to_inch(160)))
    gs = gridspec.GridSpec(2, 2, figure=fig,
                           hspace=0.50, wspace=0.38,
                           left=0.09, right=0.97, top=0.95, bottom=0.08)

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    add_panel_label(ax_a, 'a', x=-0.14, y=1.10)
    add_panel_label(ax_b, 'b', x=-0.14, y=1.10)
    add_panel_label(ax_c, 'c', x=-0.14, y=1.10)
    add_panel_label(ax_d, 'd', x=-0.14, y=1.10)

    print('Drawing Panel A: Enrichment fold change...')
    panel_a(ax_a, df_enrich)

    print('Drawing Panel B: Per-cluster density...')
    panel_b(ax_b, df_cluster)

    print('Drawing Panel C: DNA vs methylation...')
    panel_c(ax_c, df_dna)

    print('Drawing Panel D: Per-cluster fold enrichment...')
    panel_d(ax_d, df_cluster)

    # ── Save ─────────────────────────────────────────────────────────────────
    out_path = FIG_SUP_DIR / 'FigS11_bgc_methylation'
    save_figure(fig, out_path)
    print(f'\nSaved: {out_path}.pdf / .svg')
    print('=== Done ===')


if __name__ == '__main__':
    main()
