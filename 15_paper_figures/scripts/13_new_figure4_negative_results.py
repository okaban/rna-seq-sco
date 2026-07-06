#!/usr/bin/env python3
"""
New Figure 4: Negative Results Panel (Simpson's Paradox + Null Results)

(A) GCCGGC de-repression — unstratified vs core-only (Simpson's Paradox)
(B) AAGCCCG clean null (r=0.003, p=0.91)
(C) Sequence ≠ Methylation enrichment at exposed TF TSS
(D) No neighborhood effect (permutation p=0.857, no distance decay rho=+0.026)
"""

import importlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D
from scipy import stats

_utils = importlib.import_module('00_shared_utils')
for _attr in dir(_utils):
    if not _attr.startswith('_'):
        globals()[_attr] = getattr(_utils, _attr)


def panel_a(ax, df_trans, df_geo):
    """Panel A: Simpson's Paradox — GCCGGC unstratified vs core-only."""
    # Left side: unstratified (all regions)
    lost_all = df_trans[df_trans['transition_T1T2'] == 'Lost']['LFC_T2vsT1'].dropna()
    never_all = df_trans[df_trans['transition_T1T2'] == 'Never']['LFC_T2vsT1'].dropna()

    # Right side: core-only
    lost_core = df_trans[(df_trans['transition_T1T2'] == 'Lost') &
                         (df_trans['region'] == 'core')]['LFC_T2vsT1'].dropna()
    never_core = df_trans[(df_trans['transition_T1T2'] == 'Never') &
                          (df_trans['region'] == 'core')]['LFC_T2vsT1'].dropna()

    # Two side-by-side violin plots
    positions = [1, 2, 4, 5]
    data = [lost_all.values, never_all.values, lost_core.values, never_core.values]
    colors = [COL_ARTIFACT, COL_GRAY, COL_SIGNAL, COL_GRAY]
    labels_inner = ['Lost', 'Never', 'Lost', 'Never']

    vp = ax.violinplot(data, positions=positions, showextrema=False,
                       showmedians=False, widths=0.7)
    for i, body in enumerate(vp['bodies']):
        body.set_facecolor(colors[i])
        body.set_alpha(0.4)

    bp = ax.boxplot(data, positions=positions, widths=0.2, patch_artist=True,
                    showfliers=False,
                    medianprops=dict(color='white', linewidth=1.5),
                    whiskerprops=dict(linewidth=0.8),
                    capprops=dict(linewidth=0.8))
    for i, patch in enumerate(bp['boxes']):
        patch.set_facecolor(colors[i])
        patch.set_alpha(0.8)

    # Significance brackets
    # Unstratified: p = 8.3e-8
    y_max = 6
    ax.plot([1, 1, 2, 2], [y_max, y_max + 0.3, y_max + 0.3, y_max], 'k-', linewidth=0.8)
    ax.text(1.5, y_max + 0.4, 'p = 8.3 $\\times$ 10$^{-8}$', ha='center', fontsize=6.5,
            color=COL_ARTIFACT, fontweight='bold')

    # Core-only: p = 0.87
    ax.plot([4, 4, 5, 5], [y_max, y_max + 0.3, y_max + 0.3, y_max], 'k-', linewidth=0.8)
    ax.text(4.5, y_max + 0.4, 'p = 0.87', ha='center', fontsize=6.5,
            color=COL_SIGNAL, fontweight='bold')

    ax.set_xticks(positions)
    ax.set_xticklabels(labels_inner, fontsize=7)
    ax.set_ylabel('$\\log_2$FC (T2 vs T1)', fontsize=9)
    ax.set_ylim(-8, 9)

    # Group labels
    ax.text(1.5, -7.5, 'Unstratified\n(all regions)', ha='center', fontsize=7,
            fontweight='bold', color=COL_ARTIFACT)
    ax.text(4.5, -7.5, 'Core only', ha='center', fontsize=7,
            fontweight='bold', color=COL_SIGNAL)

    ax.set_title("Simpson's Paradox\n(GCCGGC de-repression)", fontsize=9, fontweight='bold')


def panel_b(ax, df_trans_aag):
    """Panel B: AAGCCCG clean null — no de-repression effect."""
    lost = df_trans_aag[df_trans_aag['transition'] == 'Lost']['LFC_T2vsT1'].dropna()
    never = df_trans_aag[df_trans_aag['transition'] == 'Never']['LFC_T2vsT1'].dropna()

    positions = [1, 2]
    data = [lost.values, never.values]
    colors = [COL_SIGNAL, COL_GRAY]

    vp = ax.violinplot(data, positions=positions, showextrema=False,
                       showmedians=False, widths=0.6)
    for i, body in enumerate(vp['bodies']):
        body.set_facecolor(colors[i])
        body.set_alpha(0.4)

    bp = ax.boxplot(data, positions=positions, widths=0.15, patch_artist=True,
                    showfliers=False,
                    medianprops=dict(color='white', linewidth=1.5),
                    whiskerprops=dict(linewidth=0.8),
                    capprops=dict(linewidth=0.8))
    for i, patch in enumerate(bp['boxes']):
        patch.set_facecolor(colors[i])
        patch.set_alpha(0.8)

    # Test result
    U, p = stats.mannwhitneyu(lost, never, alternative='two-sided')
    ax.text(0.5, 0.97, f'MW p = {p:.2f}\nr = 0.003',
            transform=ax.transAxes, ha='center', va='top', fontsize=7,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='#ccc', alpha=0.9))

    ax.set_xticks(positions)
    ax.set_xticklabels([f'Lost\n(n={len(lost)})', f'Never\n(n={len(never)})'], fontsize=7)
    ax.set_ylabel('$\\log_2$FC (T2 vs T1)', fontsize=9)
    ax.set_title('AAGCCCG null result\n(no geographic confound)', fontsize=9, fontweight='bold')


def panel_c(ax, df_seq_tests, df_methyl_status):
    """Panel C: Sequence ≠ Methylation — AAGCCCG sequence enrichment vs methylation presence."""
    # AAGCCCG sequence is 5.0x enriched at exposed TF TSS (from sequence analysis)
    aag_row = df_seq_tests[(df_seq_tests['motif'] == 'AAGCCCG') &
                            (df_seq_tests['window'] == '±300bp')]

    # How many exposed TFs actually have methylated AAGCCCG?
    n_total = len(df_methyl_status)
    n_has_sequence = (df_methyl_status['n_AAGCCCG_T1'] > 0).sum() + \
                     (df_methyl_status['n_AAGCCCG_T2'] > 0).sum()
    # Use AAGCCCG_status column
    n_lost = (df_methyl_status['AAGCCCG_status'] == 'Lost').sum()
    n_never = (df_methyl_status['AAGCCCG_status'] == 'Never').sum()
    n_has_any_aagcccg = n_total - n_never  # Those with any AAGCCCG methylation

    pct_seq = 25.8   # from statistical_tests (exposed_has_motif_pct for AAGCCCG ±300bp)
    pct_shielded_seq = 6.3  # shielded_has_motif_pct
    pct_methyl = 100 * n_has_any_aagcccg / n_total

    # Bar chart
    x = [0, 1, 2]
    values = [pct_seq, pct_shielded_seq, pct_methyl]
    colors = [COL_EXPOSED, COL_SHIELDED, COL_4mC]
    labels = ['Exposed\nsequence', 'Shielded\nsequence', 'Exposed\nmethylated']

    bars = ax.bar(x, values, color=colors, alpha=0.8, width=0.6,
                  edgecolor='white', linewidth=0.5)

    for i, (v, bar) in enumerate(zip(values, bars)):
        ax.text(i, v + 1, f'{v:.1f}%', ha='center', va='bottom', fontsize=7,
                fontweight='bold')

    # Fold enrichment annotation
    if pct_shielded_seq > 0:
        fold = pct_seq / pct_shielded_seq
        ax.annotate(f'{fold:.1f}x enriched\n(sequence)', xy=(0, pct_seq),
                    xytext=(0.8, pct_seq + 8), fontsize=7,
                    arrowprops=dict(arrowstyle='->', color=COL_DARK, lw=0.8),
                    ha='center', color=COL_DARK)

    ax.text(2, pct_methyl + 5, f'Only {pct_methyl:.1f}%\nmethylated',
            ha='center', fontsize=7, color=COL_4mC, fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel('% with AAGCCCG (±300 bp)', fontsize=9)
    ax.set_ylim(0, 45)
    ax.set_title('Sequence $\\neq$ Methylation', fontsize=9, fontweight='bold')


def panel_d(ax, df_perm, df_stats):
    """Panel D: No neighborhood effect — permutation and distance decay."""
    # Permutation test
    perm_p = df_perm[df_perm['metric'] == 'mean_absLFC_T3']['empirical_p'].values[0]
    perm_obs = df_perm[df_perm['metric'] == 'mean_absLFC_T3']['observed'].values[0]
    perm_mean = df_perm[df_perm['metric'] == 'mean_absLFC_T3']['perm_mean'].values[0]
    perm_sd = df_perm[df_perm['metric'] == 'mean_absLFC_T3']['perm_sd'].values[0]

    # Distance decay rho
    rho_row = df_stats[df_stats['metric'] == 'Spearman_rho_distance_absLFC']
    rho = rho_row['effect_size_r'].values[0] if len(rho_row) > 0 else 0.026
    rho_p = rho_row['p_value'].values[0] if len(rho_row) > 0 else 0.197

    # Draw permutation distribution
    perm_x = np.linspace(perm_mean - 4 * perm_sd, perm_mean + 4 * perm_sd, 200)
    perm_y = stats.norm.pdf(perm_x, perm_mean, perm_sd)
    ax.fill_between(perm_x, perm_y, color=COL_GRAY, alpha=0.3)
    ax.plot(perm_x, perm_y, color=COL_GRAY, linewidth=1.5, label='Permutation\ndistribution')

    # Observed value
    ax.axvline(perm_obs, color=COL_EXPOSED, linewidth=2, linestyle='-',
               label=f'Observed\n(p = {perm_p:.3f})')
    ax.axvline(perm_mean, color='gray', linewidth=0.8, linestyle=':')

    ax.set_xlabel('Mean |$\\log_2$FC| in ±20 kb neighborhood', fontsize=9)
    ax.set_ylabel('Density', fontsize=9)
    ax.set_title('No neighborhood effect', fontsize=9, fontweight='bold')
    ax.legend(fontsize=7, loc='upper right', frameon=False)

    # Distance decay annotation
    ax.text(0.03, 0.95, f'Distance decay:\n$\\rho$ = +{rho:.3f}, p = {rho_p:.2f}',
            transform=ax.transAxes, ha='left', va='top', fontsize=7,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='#ccc', alpha=0.9))


def main():
    apply_style()
    print('=== Fig. S9: Negative Results (Supplementary) ===')
    print()

    # Load data
    print('Loading data...')
    df_trans_gcc, df_geo_gcc = load_gene_methylation_transitions('GCCGGC')
    df_trans_aag, df_geo_aag = load_gene_methylation_transitions('AAGCCCG')
    df_perm, df_stats = load_neighborhood_results()
    df_methyl_status = load_exposed_methylation_status()

    # Load sequence test results
    path_seq_tests = EPIGENOME / '53_TSS_sequence_determinants' / 'tables' / 'statistical_tests.tsv'
    df_seq_tests = pd.read_csv(path_seq_tests, sep='\t')

    # Create figure: 180mm × 160mm
    fig = plt.figure(figsize=(mm_to_inch(174), mm_to_inch(170)))

    gs = fig.add_gridspec(2, 2, hspace=0.50, wspace=0.40,
                          left=0.10, right=0.95, top=0.93, bottom=0.07)

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    print('Drawing Panel A: Simpson\'s Paradox...')
    panel_a(ax_a, df_trans_gcc, df_geo_gcc)
    add_panel_label(ax_a, 'a', x=-0.15, y=1.10)

    print('Drawing Panel B: AAGCCCG null...')
    panel_b(ax_b, df_trans_aag)
    add_panel_label(ax_b, 'b', x=-0.15, y=1.10)

    print('Drawing Panel C: Sequence ≠ Methylation...')
    panel_c(ax_c, df_seq_tests, df_methyl_status)
    add_panel_label(ax_c, 'c', x=-0.15, y=1.10)

    print('Drawing Panel D: No neighborhood effect...')
    panel_d(ax_d, df_perm, df_stats)
    add_panel_label(ax_d, 'd', x=-0.15, y=1.10)

    # Save to supplementary figures
    out_path = FIG_SUP_DIR / 'FigS9_negative_results'
    save_figure(fig, out_path)
    print()
    print('=== Done ===')


if __name__ == '__main__':
    main()
