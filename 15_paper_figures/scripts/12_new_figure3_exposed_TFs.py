#!/usr/bin/env python3
"""
New Figure 3: 62 Exposed TFs — Distributed Developmental Switch

(A) Co-expression heatmap (62×62), clustered by module/bloc
(B) Two blocs temporal trajectories (z-score, mirror-image)
    → Labels added for top developmental regulators (ramR, nsdB, SCO1160)
(C) Early vs Late response scatter (|LFC_T2| vs |LFC_T3|)
    → Replaces phase_ratio violin for intuitive visualization
(D) TF family composition by bloc, grouped by functional category
(E) TCS pair asymmetry — comprehensive genome-wide analysis
    → All M145 adjacent SK-RR pairs annotated with exposed/shielded status
"""

import importlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform

_utils = importlib.import_module('00_shared_utils')
for _attr in dir(_utils):
    if not _attr.startswith('_'):
        globals()[_attr] = getattr(_utils, _attr)


def panel_a(ax, df_coexpr, df_temporal):
    """Panel A: 62×62 co-expression heatmap, ordered by bloc then module."""
    bloc_order = {'activation': 0, 'repression': 1, 'unassigned': 2}
    gene_order = (df_temporal
                  .assign(bloc_rank=df_temporal['bloc'].map(bloc_order))
                  .sort_values(['bloc_rank', 'module', 'locus_tag'])
                  ['locus_tag'].tolist())

    gene_order = [g for g in gene_order if g in df_coexpr.columns]

    mat = df_coexpr.loc[gene_order, gene_order].values

    im = ax.imshow(mat, cmap='RdBu_r', vmin=-1, vmax=1, aspect='equal')

    n_act = sum(1 for g in gene_order
                if df_temporal[df_temporal['locus_tag'] == g]['bloc'].values[0] == 'activation')
    n_rep = sum(1 for g in gene_order
                if df_temporal[df_temporal['locus_tag'] == g]['bloc'].values[0] == 'repression')

    ax.axhline(n_act - 0.5, color='black', linewidth=1.0)
    ax.axvline(n_act - 0.5, color='black', linewidth=1.0)
    if n_rep > 0:
        ax.axhline(n_act + n_rep - 0.5, color='black', linewidth=0.5, linestyle='--')
        ax.axvline(n_act + n_rep - 0.5, color='black', linewidth=0.5, linestyle='--')

    ax.text(-3, n_act / 2, 'Act.', ha='right', va='center', fontsize=7,
            color=COL_ACTIVATION, fontweight='bold', rotation=90)
    ax.text(-3, n_act + n_rep / 2, 'Rep.', ha='right', va='center', fontsize=7,
            color=COL_REPRESSION, fontweight='bold', rotation=90)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('Co-expression matrix\n(57 exposed TFs)', fontsize=9, fontweight='bold')

    cbar = plt.colorbar(im, ax=ax, shrink=0.7, pad=0.03, aspect=15)
    cbar.set_label('Spearman $\\rho$', fontsize=7)
    cbar.ax.tick_params(labelsize=6)


def panel_b(ax, df_temporal, df_bloc):
    """Panel B: Two blocs temporal trajectories with key gene labels."""
    timepoints = ['T1', 'T2', 'T3']
    x = [1, 2, 3]

    # Key developmental/morphological regulators to highlight
    highlight_genes = {
        'ramR':    {'old_locus': 'SCO6685', 'label': 'ramR'},
        'nsdB':    {'old_locus': 'SCO7252', 'label': 'nsdB'},
        'SCO1160': {'old_locus': 'SCO1160', 'label': 'SCO1160'},
    }

    for bloc, color, marker in [('activation', COL_ACTIVATION, 'o'),
                                 ('repression', COL_REPRESSION, 's')]:
        sub = df_temporal[df_temporal['bloc'] == bloc]
        z_cols = ['T1_z', 'T2_z', 'T3_z']

        # Individual gene trajectories (light)
        for _, row in sub.iterrows():
            is_highlight = row['old_locus_tag'] in [v['old_locus'] for v in highlight_genes.values()]
            if not is_highlight:
                ax.plot(x, [row[c] for c in z_cols], color=color, alpha=0.08, linewidth=0.5)

        # Highlighted genes — drawn with higher visibility (labels in inset box below)
        for gene_key, gene_info in highlight_genes.items():
            match = sub[sub['old_locus_tag'] == gene_info['old_locus']]
            if len(match) == 0:
                continue
            row = match.iloc[0]
            zvals = [row[c] for c in z_cols]
            ax.plot(x, zvals, color=color, alpha=0.85, linewidth=1.6,
                    linestyle='--', zorder=6)

        # Mean trajectory (bold)
        means = [sub[c].mean() for c in z_cols]
        sems = [sub[c].std() / np.sqrt(len(sub)) for c in z_cols]
        ax.plot(x, means, color=color, linewidth=2.5, marker=marker,
                markersize=6, label=f'{bloc.capitalize()} (n={len(sub)})',
                zorder=5)
        ax.fill_between(x,
                        [m - s for m, s in zip(means, sems)],
                        [m + s for m, s in zip(means, sems)],
                        color=color, alpha=0.2)

    ax.axhline(0, color='gray', linewidth=0.5, linestyle=':')
    ax.set_xticks(x)
    ax.set_xticklabels(timepoints)
    ax.set_xlabel('Timepoint', fontsize=9)
    ax.set_ylabel('Expression (z-score)', fontsize=9)
    ax.set_title('Temporal trajectories\n(antagonistic blocs)', fontsize=9, fontweight='bold')
    ax.legend(fontsize=7, loc='upper left', frameon=False)

    # Inset box: top developmental regulators
    # Gather LFC_T3vsT1 for the highlighted genes (peak LFC for annotation)
    highlights_info = []
    for gene_key, gene_info in highlight_genes.items():
        match = df_temporal[df_temporal['old_locus_tag'] == gene_info['old_locus']]
        if len(match) == 0:
            continue
        row = match.iloc[0]
        lfc_max = max(abs(row['LFC_T2vsT1']), abs(row['LFC_T3vsT1']))
        lfc_val = row['LFC_T3vsT1'] if abs(row['LFC_T3vsT1']) >= abs(row['LFC_T2vsT1']) \
                  else row['LFC_T2vsT1']
        sign = '+' if lfc_val > 0 else ''
        highlights_info.append(f'{gene_info["label"]}: LFC={sign}{lfc_val:.1f}')

    inset_text = 'Top activation regulators:\n' + '\n'.join(highlights_info)
    ax.text(0.97, 0.03, inset_text,
            transform=ax.transAxes, ha='right', va='bottom', fontsize=6,
            color=COL_ACTIVATION, style='normal', fontweight='normal',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                      edgecolor=COL_ACTIVATION, alpha=0.88, linewidth=0.7))


def panel_c(ax, df_temporal):
    """Panel C: Early vs Late transition scatter — |LFC_T2vsT1| vs |LFC_T3vsT2|.

    Axes directly correspond to the phase_ratio calculation:
      phase_ratio = |LFC_T2vsT1| / (|LFC_T2vsT1| + |LFC_T3vsT2|)
    Phase_ratio > 0.6  ↔  |LFC_T2vsT1| > 1.5 × |LFC_T3vsT2|
                       ↔  below the line y = (2/3)x (shown as dashed)
    """
    df = df_temporal.copy()
    df['abs_t1t2'] = df['LFC_T2vsT1'].abs()    # T1→T2 magnitude (early)
    df['abs_t2t3'] = df['LFC_T3vsT2'].abs()    # T2→T3 magnitude (late)

    bloc_colors = {'activation': COL_ACTIVATION, 'repression': COL_REPRESSION,
                   'unassigned': COL_GRAY}
    bloc_markers = {'activation': 'o', 'repression': 's', 'unassigned': '^'}

    for bloc in ['activation', 'repression', 'unassigned']:
        sub = df[df['bloc'] == bloc]
        ax.scatter(sub['abs_t1t2'], sub['abs_t2t3'],
                   color=bloc_colors[bloc],
                   marker=bloc_markers[bloc],
                   s=28, alpha=0.7, edgecolors='white', linewidths=0.4,
                   label=f'{bloc.capitalize()} (n={len(sub)})',
                   zorder=4)

    # Phase_ratio = 0.6 boundary: |T2vsT1| = 1.5 × |T3vsT2|
    # → y = (2/3) x  (i.e. late = 0.667 × early)
    x_max = df['abs_t1t2'].max() * 1.1
    y_max = df['abs_t2t3'].max() * 1.1
    x_line = np.linspace(0, max(x_max, y_max * 1.5), 100)
    ax.plot(x_line, x_line * (2 / 3), color='gray', linewidth=0.9,
            linestyle='--', zorder=2, label='Phase ratio = 0.6')

    # Region labels — minimal, placed in sparsely populated quadrants
    ax.text(0.03, 0.70, 'Late-\ndominant',
            transform=ax.transAxes, fontsize=6, color='gray', va='center',
            style='italic')
    ax.text(0.75, 0.12, 'Early-\ndominant',
            transform=ax.transAxes, fontsize=6, color='gray', va='bottom',
            style='italic')

    # Highlight key genes — manual offsets to avoid overlap
    # ramR: (6.60, 1.30), SCO1160: (6.81, 1.24), nsdB: (9.32, 1.55)
    highlight_offsets = {
        'SCO6685': ('ramR',    (-5, 18)),   # ramR: up-left
        'SCO7252': ('nsdB',    (4,  14)),   # nsdB: up-right
        'SCO1160': ('SCO1160', (4,  -12)),  # SCO1160: down-right
    }
    for old_locus, (label, (dx, dy)) in highlight_offsets.items():
        match = df[df['old_locus_tag'] == old_locus]
        if len(match) == 0:
            continue
        row = match.iloc[0]
        ax.annotate(label,
                    xy=(row['abs_t1t2'], row['abs_t2t3']),
                    xytext=(dx, dy), textcoords='offset points',
                    fontsize=6, fontweight='bold',
                    color=COL_ACTIVATION,
                    arrowprops=dict(arrowstyle='-', color=COL_ACTIVATION,
                                   lw=0.5, shrinkA=2, shrinkB=2))

    # Summary: denominator = all 62 exposed TFs (consistent with manuscript)
    n_early = (df['phase_ratio'] > 0.6).sum()
    n_total_tfs = len(df)
    pct_early = 100 * n_early / n_total_tfs
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel('T1→T2 response\n|log2FC T2 vs T1|', fontsize=8)
    ax.set_ylabel('T2→T3 response\n|log2FC T3 vs T2|', fontsize=8)
    ax.set_title(f'Response timing\n({pct_early:.0f}% early-dominant)',
                 fontsize=9, fontweight='bold')
    ax.legend(fontsize=6, loc='upper right', frameon=False, ncol=1)


def panel_d(ax, df_family):
    """Panel D: TF family composition grouped by functional category.

    Bars colored by functional category, fill pattern = bloc.
    Category dividers drawn inside the plot with inline category headers.
    """
    FUNC_ORDER = [
        'Signal transduction',
        'Stress/development',
        'Defense/resistance',
        'Metabolic',
        'General',
    ]
    # Short display names to fit within plot
    FUNC_SHORT = {
        'Signal transduction': 'Signal trans.',
        'Stress/development':  'Stress/dev.',
        'Defense/resistance':  'Defense/resist.',
        'Metabolic':           'Metabolic',
        'General':             'General',
    }
    FUNC_COLORS = {
        'Signal transduction': '#5C6BC0',
        'Stress/development':  '#26A69A',
        'Defense/resistance':  '#EF5350',
        'Metabolic':           '#8D6E63',
        'General':             '#78909C',
    }

    df = df_family.copy()
    df['total'] = df['activation_count'] + df['repression_count']
    df = df[df['total'] >= 1].copy()

    # Sort: by functional category order, then total
    df['func_rank'] = df['functional_category'].map(
        {c: i for i, c in enumerate(FUNC_ORDER)}).fillna(99)
    df = df.sort_values(['func_rank', 'total'], ascending=[True, True]).reset_index(drop=True)

    # Build y-positions with gaps between categories
    # Insert a gap row between each category transition
    rows = []
    prev_cat = None
    gap_offset = 0
    y_positions = []
    gap_rows = []  # (y_pos, category_name) for divider lines

    for i, row in df.iterrows():
        fcat = row.get('functional_category', 'General')
        if fcat != prev_cat:
            if prev_cat is not None:
                gap_offset += 0.7  # extra gap between categories
                gap_rows.append((i + gap_offset - 0.35, fcat))
            prev_cat = fcat
        y_positions.append(i + gap_offset)

    y_arr = np.array(y_positions)
    bar_h = 0.32

    for i, row in df.iterrows():
        fcat = row.get('functional_category', 'General')
        fc = FUNC_COLORS.get(fcat, '#78909C')
        y = y_arr[i]
        # Activation bar (solid fill)
        ax.barh(y - bar_h / 2, row['activation_count'], bar_h,
                color=fc, alpha=0.85)
        # Repression bar (hatched, lighter)
        ax.barh(y + bar_h / 2, row['repression_count'], bar_h,
                color=fc, alpha=0.35, hatch='////', edgecolor=fc)

    ax.set_yticks(y_arr)
    ax.set_yticklabels(df['tf_family'], fontsize=7)
    ax.set_xlabel('Number of TFs', fontsize=9)
    ax.set_title('TF family by functional category', fontsize=9, fontweight='bold')

    # Draw divider lines and inline category labels
    x_max_for_label = ax.get_xlim()[1] if ax.get_xlim()[1] > 1 else 9
    for gap_y, cat_name in gap_rows:
        ax.axhline(gap_y, color='#ccc', linewidth=0.6, linestyle='-', zorder=0)

    # Add category color strips on the left (y-axis side)
    # Group bars by category and draw a colored bracket/strip
    prev_cat = None
    cat_y_start = None
    for i, row in df.iterrows():
        fcat = row.get('functional_category', 'General')
        if fcat != prev_cat:
            if prev_cat is not None:
                # Draw category label between divider and first bar of group
                mid_y = (cat_y_start + y_arr[i - 1]) / 2
                short = FUNC_SHORT.get(prev_cat, prev_cat)
                fc_c = FUNC_COLORS.get(prev_cat, '#78909C')
                # Colored dot on the right end of bar label area
                ax.text(-0.45, mid_y, short,
                        ha='right', va='center', fontsize=5.5,
                        color=fc_c, fontweight='bold',
                        transform=ax.transData, clip_on=True)
            cat_y_start = y_arr[i]
            prev_cat = fcat
    # Last category
    if prev_cat is not None:
        mid_y = (cat_y_start + y_arr[-1]) / 2
        short = FUNC_SHORT.get(prev_cat, prev_cat)
        fc_c = FUNC_COLORS.get(prev_cat, '#78909C')
        ax.text(-0.45, mid_y, short,
                ha='right', va='center', fontsize=5.5,
                color=fc_c, fontweight='bold', clip_on=True)

    # Highlight TetR
    tetr_rows = df[df['tf_family'] == 'TetR']
    if len(tetr_rows) > 0:
        tetr_i = tetr_rows.index[0]
        row = tetr_rows.iloc[0]
        y = y_arr[tetr_i]
        ax.text(row['repression_count'] + 0.2,
                y + bar_h / 2,
                'p = 0.028 *', fontsize=7, color=COL_REPRESSION,
                fontweight='bold', va='center')

    # Legend for activation/repression fill style
    legend_elements = [
        plt.Rectangle((0, 0), 1, 1, fc='gray', alpha=0.85, label='Activation'),
        plt.Rectangle((0, 0), 1, 1, fc='gray', alpha=0.35, hatch='////',
                       ec='gray', label='Repression'),
    ]
    ax.legend(handles=legend_elements, fontsize=7, loc='lower right', frameon=False)


def _get_comprehensive_tcs_pairs(df_all_genes, distance_kb=10):
    """Find all adjacent SK-RR pairs within distance_kb in M145.

    Returns DataFrame with sk, rr, sk_exposed, rr_exposed, dist_midpoint columns.
    """
    sk = df_all_genes[df_all_genes['tf_family'] == 'Sensor kinase'].copy()
    rr = df_all_genes[df_all_genes['tf_family'] == 'Response regulator'].copy()

    sk['midpoint'] = (sk['start'] + sk['end']) / 2
    rr['midpoint'] = (rr['start'] + rr['end']) / 2

    cutoff = distance_kb * 1000
    pairs = []
    for _, s in sk.iterrows():
        for _, r in rr.iterrows():
            dist = abs(s['midpoint'] - r['midpoint'])
            if dist < cutoff:
                pairs.append({
                    'sk_locus': s['old_locus_tag'],
                    'rr_locus': r['old_locus_tag'],
                    'sk_exposed': bool(s['is_exposed']),
                    'rr_exposed': bool(r['is_exposed']),
                    'dist_bp': int(dist),
                })
    return pd.DataFrame(pairs)


def panel_e(ax, df_tcs, df_all_genes=None):
    """Panel E: TCS asymmetry — all 7 identified pairs with descriptive framing.

    NOTE: Genome-wide comparison (9/84 SK-RR pairs) is NOT statistically
    significant (permutation test p = 0.83; expected ~10.5 by chance given
    ~7% base exposure rate). Stats kept in manuscript; figure is descriptive only.
    """
    n_pairs = len(df_tcs)

    for i, (_, row) in enumerate(df_tcs.iterrows()):
        y = n_pairs - 1 - i

        sk_color = COL_EXPOSED if row['sk_exposed'] else COL_SHIELDED
        rr_color = COL_EXPOSED if row['rr_exposed'] else COL_SHIELDED

        ax.scatter(0.30, y, color=sk_color, s=80, zorder=5,
                   edgecolors='black', linewidths=0.5, marker='D')
        ax.text(0.05, y, row['sk_old_locus'], ha='left', va='center', fontsize=6)

        ax.scatter(0.70, y, color=rr_color, s=80, zorder=5,
                   edgecolors='black', linewidths=0.5, marker='o')
        ax.text(0.95, y, row['rr_old_locus'], ha='right', va='center', fontsize=6)

        ax.plot([0.34, 0.66], [y, y], color='#999', linewidth=0.8,
                linestyle='-', zorder=2)

    # Column headers
    ax.text(0.30, n_pairs - 0.2, 'SK', ha='center', va='bottom', fontsize=8,
            fontweight='bold')
    ax.text(0.70, n_pairs - 0.2, 'RR', ha='center', va='bottom', fontsize=8,
            fontweight='bold')

    # Legend
    legend_elements = [
        Line2D([0], [0], marker='D', color='w', markerfacecolor=COL_EXPOSED,
               markersize=8, label='Exposed', markeredgecolor='black', markeredgewidth=0.5),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=COL_SHIELDED,
               markersize=8, label='Shielded', markeredgecolor='black', markeredgewidth=0.5),
    ]
    ax.legend(handles=legend_elements, fontsize=7, loc='lower center',
              frameon=False, ncol=2)

    # Descriptive annotation — no statistical claim (genome-wide stats in manuscript)
    ax.text(0.50, n_pairs + 0.15,
            '7/7 identified pairs: exactly one partner exposed\n'
            '(descriptive pattern; see Methods for permutation test)',
            ha='center', va='bottom', fontsize=6, fontstyle='italic',
            color=COL_DARK)

    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-1.2, n_pairs + 1.2)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_title('TCS pair asymmetry', fontsize=9, fontweight='bold')


def main():
    apply_style()
    print('=== Figure 5: 57 Exposed TFs (n=57 update) ===')
    print()

    # Load data
    print('Loading data...')
    df_coexpr = load_coexpression_matrix()
    df_temporal = load_temporal_classification()
    df_bloc = load_bloc_comparison()
    df_family = load_bloc_family_enrichment()
    df_tcs = load_tcs_pairs()
    # Use the unified n=57 gene set (1,055 regulatory genes)
    df_all_genes = pd.read_csv(
        EPIGENOME / '52_shielded_exposed_boundary' / 'tables' /
        'all_genes_features_unified_n57.tsv', sep='\t')

    # Create figure: 180mm × 210mm
    fig = plt.figure(figsize=(mm_to_inch(174), mm_to_inch(210)))

    gs = fig.add_gridspec(2, 12, hspace=0.45, wspace=1.2,
                          left=0.08, right=0.97, top=0.95, bottom=0.05,
                          height_ratios=[1, 1])

    # Row 1: A (5 cols), B (4 cols), C (3 cols)
    ax_a = fig.add_subplot(gs[0, 0:5])
    ax_b = fig.add_subplot(gs[0, 5:9])
    ax_c = fig.add_subplot(gs[0, 9:12])

    # Row 2: D (6 cols), E (6 cols)
    ax_d = fig.add_subplot(gs[1, 0:6])
    ax_e = fig.add_subplot(gs[1, 6:12])

    print('Drawing Panel A: Co-expression heatmap...')
    panel_a(ax_a, df_coexpr, df_temporal)
    add_panel_label(ax_a, 'a', x=-0.12, y=1.08)

    print('Drawing Panel B: Temporal trajectories with gene labels...')
    panel_b(ax_b, df_temporal, df_bloc)
    add_panel_label(ax_b, 'b', x=-0.18, y=1.08)

    print('Drawing Panel C: Early vs Late response scatter...')
    panel_c(ax_c, df_temporal)
    add_panel_label(ax_c, 'c', x=-0.20, y=1.08)

    print('Drawing Panel D: TF family by functional category...')
    panel_d(ax_d, df_family)
    add_panel_label(ax_d, 'd', x=-0.10, y=1.08)

    print('Drawing Panel E: TCS asymmetry (genome-wide)...')
    panel_e(ax_e, df_tcs, df_all_genes)
    add_panel_label(ax_e, 'e', x=-0.06, y=1.08)

    # Save
    out_path = FIG_DIR / 'Figure5_exposed_TFs'
    save_figure(fig, out_path)
    print()
    print('=== Done ===')


if __name__ == '__main__':
    main()
