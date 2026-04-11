#!/usr/bin/env python3
"""
New Figure 3: Protection Zone Characterization (Layer 2 quantification)

(A) TSS methylation gradient heatmap (site_type × distance bins)
(B) 293bp boundary ROC curve (nearest_methyl_distance AUC=0.917 vs baseMean AUC=0.547)
(C) Sequence vs Protein occupancy decomposition (CV AUC comparison)
(D) Individual TF BS is NOT protective (TFBS-centered methylation profile)
(E) Expression-independence (exposed fraction by expression quintile)
"""

import importlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from sklearn.metrics import roc_curve, auc

_utils = importlib.import_module('00_shared_utils')
for _attr in dir(_utils):
    if not _attr.startswith('_'):
        globals()[_attr] = getattr(_utils, _attr)


def panel_a(ax, df_spatial):
    """Panel A: TSS methylation gradient heatmap — site type × distance."""
    # Show regulatory genes, T1, different methylation types
    categories = [
        ('GCCGGC_4mC', 'GCCGGC (4mC)'),
        ('All_4mC', 'All 4mC'),
        ('AAGCCCG_6mA', 'AAGCCCG (6mA)'),
        ('All_6mA', 'All 6mA'),
        ('All_methylation', 'All methylation'),
    ]

    matrix_rows = []
    row_labels = []
    for site_type, label in categories:
        sub = df_spatial[(df_spatial['site_type'] == site_type) &
                         (df_spatial['timepoint'] == 'T1') &
                         (df_spatial['gene_category'] == 'regulatory')]
        sub = sub.sort_values('bin_center_bp')
        if len(sub) == 0:
            continue
        matrix_rows.append(sub['density_per_kb_per_gene'].values)
        row_labels.append(label)

    if not matrix_rows:
        ax.text(0.5, 0.5, 'No data', transform=ax.transAxes, ha='center')
        return

    matrix = np.array(matrix_rows)
    bin_centers = df_spatial[(df_spatial['site_type'] == 'All_methylation') &
                             (df_spatial['timepoint'] == 'T1') &
                             (df_spatial['gene_category'] == 'regulatory')].sort_values('bin_center_bp')['bin_center_bp'].values

    # Normalize each row by its flanking mean (±3-5 kb) to show relative depletion
    # Apply 5-bin running mean to smooth out noise (especially AAGCCCG with sparse sites)
    norm_matrix = np.zeros_like(matrix)
    for i in range(len(matrix)):
        smoothed = np.convolve(matrix[i], np.ones(5) / 5, mode='same')
        flank_mask = (np.abs(bin_centers) >= 3000) & (np.abs(bin_centers) <= 5000)
        flank_mean = smoothed[flank_mask].mean() if flank_mask.sum() > 0 else smoothed.mean()
        norm_matrix[i] = smoothed / flank_mean if flank_mean > 0 else 1.0

    # Clip extreme values to prevent saturated patches
    norm_matrix = np.clip(norm_matrix, 0.5, 1.5)

    im = ax.imshow(norm_matrix, aspect='auto', cmap='RdBu_r',
                   vmin=0.6, vmax=1.4,
                   extent=[bin_centers[0] / 1000, bin_centers[-1] / 1000, len(row_labels) - 0.5, -0.5])

    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=7)
    ax.set_xlabel('Distance from TSS (kb)', fontsize=9)
    ax.set_title('Methylation at TSS\n(regulatory, T1)', fontsize=9, fontweight='bold')

    # Protection zone boundary lines
    ax.axvline(-1.3, color='white', linewidth=0.8, linestyle='--')
    ax.axvline(0.7, color='white', linewidth=0.8, linestyle='--')
    ax.axvline(0, color='white', linewidth=0.5, linestyle=':')

    # Colorbar — on the right, small and tight
    cbar = plt.colorbar(im, ax=ax, shrink=0.7, pad=0.03, aspect=15)
    cbar.set_label('Fold vs\nflanking', fontsize=6, labelpad=2)
    cbar.ax.tick_params(labelsize=6)


def panel_b(ax, df_genes):
    """Panel B: ROC curves — nearest_methyl_distance (AUC=0.917) vs baseMean (AUC=0.547)."""
    y_true = df_genes['is_exposed'].values
    valid = ~np.isnan(df_genes['nearest_methyl_distance'].values)

    # nearest_methyl_distance — lower = exposed
    dist = df_genes.loc[valid, 'nearest_methyl_distance'].values
    y_val = y_true[valid]
    fpr_d, tpr_d, _ = roc_curve(y_val, -dist)  # negate for lower=positive
    auc_d = auc(fpr_d, tpr_d)

    # baseMean — lower = exposed
    bm = df_genes.loc[valid, 'baseMean'].values
    fpr_b, tpr_b, _ = roc_curve(y_val, -bm)
    auc_b = auc(fpr_b, tpr_b)

    ax.plot(fpr_d, tpr_d, color=COL_EXPOSED, linewidth=2,
            label=f'Nearest distance\nAUC = {auc_d:.3f}')
    ax.plot(fpr_b, tpr_b, color=COL_GRAY, linewidth=1.5, linestyle='--',
            label=f'Expression level\nAUC = {auc_b:.3f}')
    ax.plot([0, 1], [0, 1], 'k:', linewidth=0.5, alpha=0.5)

    # 293bp threshold marker
    # Find the point on the ROC closest to threshold=293
    thresh_idx = np.argmin(np.abs(-dist[:, np.newaxis] - np.array([[-293]])), axis=0)

    ax.set_xlabel('False positive rate', fontsize=9)
    ax.set_ylabel('True positive rate', fontsize=9)
    ax.set_title('Exposed TF classification', fontsize=9, fontweight='bold')
    ax.legend(fontsize=7, loc='lower right', frameon=True, fancybox=False,
              edgecolor='#ccc')
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_aspect('equal')


def panel_c(ax, df_roc_boundary, df_roc_sequence):
    """Panel C: Sequence vs Protein occupancy — AUC comparison bar chart."""
    # Key features to compare
    features = [
        ('Distance\n(observed)', df_roc_boundary[df_roc_boundary['feature'] == 'nearest_methyl_distance'].iloc[0]),
        ('Sequence\n(5-fold CV)', df_roc_sequence[df_roc_sequence['feature'] == 'combined_LR_sequence_CV'].iloc[0]),
        ('GC%\n(±300 bp)', df_roc_sequence[df_roc_sequence['feature'] == 'GC_300bp'].iloc[0]),
        ('Expression\nlevel', df_roc_boundary[df_roc_boundary['feature'] == 'baseMean'].iloc[0]),
    ]

    labels = [f[0] for f in features]
    aucs = [f[1]['AUC'] for f in features]
    ci_lo = [f[1].get('AUC_CI_lo', f[1].get('AUC_95CI_lower', np.nan)) for f in features]
    ci_hi = [f[1].get('AUC_CI_hi', f[1].get('AUC_95CI_upper', np.nan)) for f in features]
    errors = [[max(0, a - lo) if not np.isnan(lo) else 0 for a, lo in zip(aucs, ci_lo)],
              [max(0, hi - a) if not np.isnan(hi) else 0 for a, hi in zip(aucs, ci_hi)]]

    colors = [COL_EXPOSED, COL_SIGNAL, COL_GRAY, COL_GRAY]
    bars = ax.bar(range(len(labels)), aucs, color=colors, alpha=0.8,
                  edgecolor='white', linewidth=0.5, width=0.6)
    ax.errorbar(range(len(labels)), aucs, yerr=errors,
                fmt='none', ecolor='#333', capsize=3, linewidth=1)

    # Value labels
    for i, (v, bar) in enumerate(zip(aucs, bars)):
        ax.text(i, v + errors[1][i] + 0.02, f'{v:.3f}',
                ha='center', va='bottom', fontsize=7, fontweight='bold')

    ax.axhline(0.5, color='gray', linewidth=0.5, linestyle=':')
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel('AUC', fontsize=9)
    ax.set_ylim(0.35, 1.05)
    ax.set_title('Feature discriminative power', fontsize=9, fontweight='bold')

    # Annotation: sequence ~33%, occupancy ~67%
    ax.annotate('~33% sequence\n~67% protein occupancy',
                xy=(0.5, 0.03), xycoords='axes fraction',
                fontsize=7, ha='center', va='bottom',
                fontstyle='italic', color=COL_DARK)


def panel_d(ax, df_tfbs):
    """Panel D: Individual TF BS is NOT protective — TFBS-centered profile."""
    # Plot observed vs random for different methylation types
    types_to_plot = [
        ('GCCGGC_T1', 'GCCGGC', COL_4mC),
        ('All_T1', 'All methylation', COL_DARK),
    ]

    for mtype, label, color in types_to_plot:
        sub = df_tfbs[df_tfbs['methylation_type'] == mtype].sort_values('bin_center_bp')
        if len(sub) == 0:
            continue
        x = sub['bin_center_bp'].values
        obs_ratio = sub['obs_over_random'].values

        ax.plot(x, obs_ratio, color=color, linewidth=1.5, label=label)

    # Reference line at 1.0 (no enrichment/depletion)
    ax.axhline(1.0, color='gray', linewidth=0.8, linestyle='--')

    # GCCGGC enrichment annotation
    ax.text(0.97, 0.95, 'No depletion at TFBS\n(GCCGGC fold = 1.16)',
            transform=ax.transAxes, ha='right', va='top', fontsize=7,
            color=COL_4mC, fontstyle='italic',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

    ax.set_xlabel('Distance from TFBS center (bp)', fontsize=9)
    ax.set_ylabel('Observed / Random', fontsize=9)
    ax.set_title('TFBS methylation profile', fontsize=9, fontweight='bold')
    ax.legend(fontsize=7, loc='lower left', frameon=False)


def panel_e(ax, df_quintile):
    """Panel E: Expression-independence — exposed fraction by expression quintile."""
    q = df_quintile['quintile'].values
    frac = df_quintile['frac_exposed'].values * 100  # Convert to percentage

    ax.bar(q, frac, color=COL_EXPOSED, alpha=0.7, edgecolor='white',
           linewidth=0.5, width=0.6)

    # Value labels
    for i, (qi, fi) in enumerate(zip(q, frac)):
        ax.text(qi, fi + 0.3, f'{fi:.1f}%', ha='center', va='bottom', fontsize=7)

    # Trend line
    slope, intercept, r, p, se = __import__('scipy').stats.linregress(q, frac)
    x_fit = np.array([0.5, 5.5])
    ax.plot(x_fit, slope * x_fit + intercept, 'k--', linewidth=0.8, alpha=0.5)

    ax.text(0.97, 0.95, f'baseMean AUC = 0.547\nJT p = 0.730',
            transform=ax.transAxes, ha='right', va='top', fontsize=7,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='#ccc', alpha=0.9))

    ax.set_xlabel('Expression quintile\n(1 = lowest, 5 = highest)', fontsize=9)
    ax.set_ylabel('% Exposed', fontsize=9)
    ax.set_ylim(0, max(frac) * 1.4)
    ax.set_xticks(q)
    ax.set_title('Expression-independence', fontsize=9, fontweight='bold')


def main():
    apply_style()
    print('=== New Figure 3: Protection Zone Characterization ===')
    print()

    # Load data
    print('Loading data...')
    df_spatial = load_spatial_profile()
    df_genes = load_all_genes_features()
    df_roc_boundary = load_roc_analysis('boundary')
    df_roc_sequence = load_roc_analysis('sequence')
    df_tfbs = load_tfbs_spatial_profile()
    df_quintile = load_expression_quintile()

    # Create figure: 180mm × 180mm, 2 rows (3 + 2)
    fig = plt.figure(figsize=(mm_to_inch(180), mm_to_inch(185)))

    gs = fig.add_gridspec(2, 12, hspace=0.55, wspace=1.2,
                          left=0.07, right=0.97, top=0.93, bottom=0.07,
                          height_ratios=[1, 1])

    # Row 1: A (5 cols with colorbar space), B (3 cols), C (4 cols)
    ax_a = fig.add_subplot(gs[0, 0:4])
    ax_b = fig.add_subplot(gs[0, 5:8])
    ax_c = fig.add_subplot(gs[0, 9:12])

    # Row 2: D (6 cols), E (6 cols)
    ax_d = fig.add_subplot(gs[1, 0:6])
    ax_e = fig.add_subplot(gs[1, 6:12])

    print('Drawing Panel A: TSS methylation gradient heatmap...')
    panel_a(ax_a, df_spatial)
    add_panel_label(ax_a, 'a', x=-0.20, y=1.12)

    print('Drawing Panel B: ROC curves...')
    panel_b(ax_b, df_genes)
    add_panel_label(ax_b, 'b', x=-0.20, y=1.12)

    print('Drawing Panel C: Sequence vs occupancy...')
    panel_c(ax_c, df_roc_boundary, df_roc_sequence)
    add_panel_label(ax_c, 'c', x=-0.20, y=1.12)

    print('Drawing Panel D: TFBS not protective...')
    panel_d(ax_d, df_tfbs)
    add_panel_label(ax_d, 'd', x=-0.12, y=1.12)

    print('Drawing Panel E: Expression-independence...')
    panel_e(ax_e, df_quintile)
    add_panel_label(ax_e, 'e', x=-0.12, y=1.12)

    # Save
    out_path = FIG_DIR / 'new_Figure3_protection_zone'
    save_figure(fig, out_path)
    print()
    print('=== Done ===')


if __name__ == '__main__':
    main()
