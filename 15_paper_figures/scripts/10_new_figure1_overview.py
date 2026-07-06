#!/usr/bin/env python3
"""
New Figure 2: Gatekeeper Model Overview
4 panels showing the 4-layer architecture.

(A) Chromosome ideogram: GCCGGC vs AAGCCCG geographic redistribution (T1→T2→T3)
(B) Methylation density around TSS (visual trend; regulatory vs non-regulatory n.s.)
(C) Shielded/Exposed distance distribution (293bp ROC threshold, AUC=0.917; n=1055 full TF dataset)
(D) 15 exposed TF expression heatmap (activation vs repression blocs)
"""

import sys
import importlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
from scipy import stats

_utils = importlib.import_module('00_shared_utils')
for _attr in dir(_utils):
    if not _attr.startswith('_'):
        globals()[_attr] = getattr(_utils, _attr)

# ── Constants ──────────────────────────────────────────────────────────────
GENOME_LEN = 8_667_507
ARM_LEFT = 1_500_000
ARM_RIGHT = 7_167_507

COL_CORE = '#4169E1'   # Royal blue
COL_ARM = '#DAA520'    # Goldenrod


def _load_aagcccg_sites():
    """Load AAGCCCG 6mA sites per timepoint from full HC data."""
    df_hc = load_methylation_hc_all()
    df_6ma = df_hc[df_hc['mod_type'] == '6mA']

    # Get AAGCCCG positions from census
    df_census = pd.read_csv(
        EPIGENOME / '23_expanded_motif_search' / '6mA_final_census.csv')
    aagcccg_pos = set(
        df_census[df_census['final_motif'] == 'AAGCCCG']['position'])

    df_aag = df_6ma[df_6ma['position'].isin(aagcccg_pos)].copy()
    print(f'  AAGCCCG sites: T1={len(df_aag[df_aag.timepoint=="T1"])}, '
          f'T2={len(df_aag[df_aag.timepoint=="T2"])}, '
          f'T3={len(df_aag[df_aag.timepoint=="T3"])}')
    return df_aag


def _draw_chromosome_row(ax, y, tp_data, chrom_h, color, max_dens,
                         bar_h, bin_size_bp, bins_bp, bin_centers_mb,
                         bin_size_mb, above=True):
    """Draw density bars for a single motif on one timepoint row."""
    if len(tp_data) == 0:
        return
    h, _ = np.histogram(tp_data['position'].values, bins=bins_bp)
    d = h / (bin_size_bp / 1000)
    heights = d / max_dens * bar_h if max_dens > 0 else d * 0

    if above:
        ax.bar(bin_centers_mb, heights, width=bin_size_mb,
               bottom=y + chrom_h / 2,
               color=color, alpha=0.7, edgecolor='none')
    else:
        ax.bar(bin_centers_mb, -heights, width=bin_size_mb,
               bottom=y - chrom_h / 2,
               color=color, alpha=0.7, edgecolor='none')


def panel_a(ax, df_gccggc, df_aagcccg):
    """Panel A: GCCGGC (4mC) vs AAGCCCG (6mA) geographic redistribution."""
    timepoints = ['T1', 'T2', 'T3']
    tp_labels = TP_LABELS
    y_positions = [2.5, 1.5, 0.5]
    chrom_h = 0.18
    bar_h = 0.40
    genome_mb = GENOME_LEN / 1e6

    bin_size_bp = 100_000
    bin_size_mb = bin_size_bp / 1e6
    bins_bp = np.arange(0, GENOME_LEN + bin_size_bp, bin_size_bp)
    bin_centers_mb = (bins_bp[:-1] + bins_bp[1:]) / 2 / 1e6

    # Global max density for consistent scaling
    max_dens = 0
    for tp in timepoints:
        for df_mod in [df_gccggc, df_aagcccg]:
            sub = df_mod[df_mod['timepoint'] == tp]
            if len(sub) > 0:
                h, _ = np.histogram(sub['position'].values, bins=bins_bp)
                d = h / (bin_size_bp / 1000)
                max_dens = max(max_dens, d.max())

    for tp, tp_label, y in zip(timepoints, tp_labels, y_positions):
        tp_gccggc = df_gccggc[df_gccggc['timepoint'] == tp]
        tp_aagcccg = df_aagcccg[df_aagcccg['timepoint'] == tp]
        n_gcc = len(tp_gccggc)
        n_aag = len(tp_aagcccg)

        n_core_gcc = ((tp_gccggc['position'] >= ARM_LEFT) &
                      (tp_gccggc['position'] <= ARM_RIGHT)).sum()
        pct_core_gcc = 100 * n_core_gcc / n_gcc if n_gcc > 0 else 0

        n_core_aag = ((tp_aagcccg['position'] >= ARM_LEFT) &
                      (tp_aagcccg['position'] <= ARM_RIGHT)).sum()
        pct_core_aag = 100 * n_core_aag / n_aag if n_aag > 0 else 0

        # Chromosome
        ax.add_patch(Rectangle((0, y - chrom_h / 2), ARM_LEFT / 1e6, chrom_h,
                               facecolor=COL_ARM, alpha=0.25, edgecolor='none'))
        ax.add_patch(Rectangle((ARM_LEFT / 1e6, y - chrom_h / 2),
                               (ARM_RIGHT - ARM_LEFT) / 1e6, chrom_h,
                               facecolor=COL_CORE, alpha=0.12, edgecolor='none'))
        ax.add_patch(Rectangle((ARM_RIGHT / 1e6, y - chrom_h / 2),
                               (GENOME_LEN - ARM_RIGHT) / 1e6, chrom_h,
                               facecolor=COL_ARM, alpha=0.25, edgecolor='none'))
        # Outline (linear ends)
        ax.plot([0, genome_mb], [y - chrom_h / 2] * 2, color='#37474F', lw=0.8)
        ax.plot([0, genome_mb], [y + chrom_h / 2] * 2, color='#37474F', lw=0.8)
        ax.plot([0, 0], [y - chrom_h / 2, y + chrom_h / 2], color='#37474F', lw=1.0)
        ax.plot([genome_mb, genome_mb], [y - chrom_h / 2, y + chrom_h / 2],
                color='#37474F', lw=1.0)

        # GCCGGC density (above, red)
        _draw_chromosome_row(ax, y, tp_gccggc, chrom_h, COL_4mC,
                             max_dens, bar_h, bin_size_bp, bins_bp,
                             bin_centers_mb, bin_size_mb, above=True)
        # AAGCCCG density (below, blue)
        _draw_chromosome_row(ax, y, tp_aagcccg, chrom_h, COL_6mA,
                             max_dens, bar_h, bin_size_bp, bins_bp,
                             bin_centers_mb, bin_size_mb, above=False)

        # Labels
        ax.text(-0.3, y + 0.15, tp_label, ha='right', va='center',
                fontsize=8, fontweight='bold')
        # GCCGGC stats (above)
        gcc_pct_label = (f'{pct_core_gcc:.0f}% core'
                         if pct_core_gcc > 50
                         else f'{100 - pct_core_gcc:.0f}% arm')
        ax.text(genome_mb + 0.2, y + 0.15,
                f'GCCGGC: {n_gcc:,} ({gcc_pct_label})',
                ha='left', va='center', fontsize=6.5, color=COL_4mC)
        # AAGCCCG stats (below)
        ax.text(genome_mb + 0.2, y - 0.15,
                f'AAGCCCG: {n_aag:,} ({pct_core_aag:.0f}% core)',
                ha='left', va='center', fontsize=6.5, color=COL_6mA)

    # Boundaries
    for bnd in [ARM_LEFT, ARM_RIGHT]:
        ax.axvline(bnd / 1e6, color='gray', linestyle=':', linewidth=0.5,
                   ymin=0.02, ymax=0.98)

    # Region labels
    ax.text(ARM_LEFT / 2 / 1e6, -0.15, 'Left arm', ha='center', va='center',
            fontsize=7, color=COL_ARM, fontstyle='italic')
    ax.text((ARM_LEFT + ARM_RIGHT) / 2 / 1e6, -0.15, 'Core',
            ha='center', va='center', fontsize=7, color=COL_CORE,
            fontstyle='italic')
    ax.text((ARM_RIGHT + GENOME_LEN) / 2 / 1e6, -0.15, 'Right arm',
            ha='center', va='center', fontsize=7, color=COL_ARM,
            fontstyle='italic')

    ax.set_xlim(-0.3, genome_mb + 4.5)
    ax.set_ylim(-0.4, 3.5)
    ax.set_xlabel('Chromosome position (Mb)', fontsize=9)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)

    # Legend
    legend_elements = [
        Line2D([0], [0], color=COL_4mC, lw=6, alpha=0.7,
               label='GCCGGC (4mC)'),
        Line2D([0], [0], color=COL_6mA, lw=6, alpha=0.7,
               label='AAGCCCG (6mA)'),
    ]
    ax.legend(handles=legend_elements, fontsize=7, loc='upper right',
              frameon=True, framealpha=0.9, edgecolor='#CCC', ncol=2)

    ax.set_title('Motif-specific geographic redistribution',
                 fontsize=10, fontweight='bold')


def panel_b(ax, df_spatial):
    """Panel B: TSS protection zone — methylation density ±5kb of TSS."""
    mask = ((df_spatial['site_type'] == 'All_methylation') &
            (df_spatial['timepoint'] == 'T1'))

    for cat, color, label, ls in [
        ('regulatory', COL_EXPOSED, 'Regulatory (n = 450)', '-'),
        ('all', COL_GRAY, 'All genes with exp. TSS (n = 2,646)', '--'),
    ]:
        sub = df_spatial[mask & (df_spatial['gene_category'] == cat)].sort_values('bin_center_bp')
        if len(sub) == 0:
            continue
        x = sub['bin_center_bp'].values
        y = sub['density_per_kb_per_gene'].values
        ax.plot(x / 1000, y, color=color, linewidth=1.5, label=label, linestyle=ls)
        if 'ci_low' in sub.columns:
            ax.fill_between(x / 1000, sub['ci_low'].values, sub['ci_high'].values,
                            color=color, alpha=0.15)

    # Vertical TSS reference line
    ax.axvline(0, color='#AAAAAA', linewidth=0.8, linestyle=':', zorder=1)

    # Note: regulatory vs non-regulatory difference not statistically significant
    # by permutation test (p_adj > 0.05 at all bins near TSS). Shown as visual trend only.
    ax.text(0.98, 0.05,
            'Regulatory vs. non-regulatory\ndifference n.s. (permutation test)',
            transform=ax.transAxes, fontsize=6, ha='right', va='bottom',
            color='#888888', fontstyle='italic')

    ax.set_xlabel('Distance from TSS (kb)', fontsize=9)
    ax.set_ylabel('Methylation density\n(sites/kb/gene)', fontsize=9)
    ax.set_title('Methylation density around TSS\n(experimental TSSs, Jeong et al. 2016)',
                 fontsize=9, fontweight='bold')
    ax.legend(fontsize=7, loc='upper left', frameon=False)
    ax.set_xlim(-5, 5)


def panel_c(ax, df_genes):
    """Panel C: Shielded vs Exposed distance distribution — violin + box (full TF dataset, 293bp ROC threshold, AUC=0.917)."""
    exposed = df_genes[df_genes['is_exposed'] == 1]['nearest_methyl_distance'].values
    shielded = df_genes[df_genes['is_exposed'] == 0]['nearest_methyl_distance'].values

    log_exp = np.log10(exposed + 1)
    log_shi = np.log10(shielded + 1)

    positions = [1, 2]
    vp = ax.violinplot([log_exp, log_shi], positions=positions,
                       showextrema=False, showmedians=False, widths=0.7)
    for i, body in enumerate(vp['bodies']):
        body.set_facecolor([COL_EXPOSED, COL_SHIELDED][i])
        body.set_alpha(0.4)
        body.set_edgecolor([COL_EXPOSED, COL_SHIELDED][i])

    bp = ax.boxplot([log_exp, log_shi], positions=positions,
                    widths=0.15, patch_artist=True, showfliers=False,
                    medianprops=dict(color='white', linewidth=1.5),
                    whiskerprops=dict(linewidth=0.8),
                    capprops=dict(linewidth=0.8))
    for i, patch in enumerate(bp['boxes']):
        patch.set_facecolor([COL_EXPOSED, COL_SHIELDED][i])
        patch.set_alpha(0.8)

    # 293bp boundary (Youden threshold from full TF dataset ROC analysis, AUC=0.917)
    log_boundary = np.log10(293)
    ax.axhline(log_boundary, color=COL_4mC, linewidth=1.5, linestyle='--', zorder=3)
    ax.text(2.65, log_boundary + 0.08, '293 bp\n(ROC optimal)',
            fontsize=6.5, color=COL_4mC, fontweight='bold', va='bottom')

    # Statistics
    U, p = stats.mannwhitneyu(exposed, shielded, alternative='less')
    ax.text(0.50, 0.97,
            f'Mann-Whitney p = {p:.1e}\nROC AUC = 0.917\nYouden\'s J threshold',
            transform=ax.transAxes, ha='center', va='top', fontsize=6.5,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='#999', alpha=0.9))

    ax.set_xticks(positions)
    ax.set_xticklabels([f'Exposed\n(n = {len(exposed)})',
                        f'Shielded\n(n = {len(shielded)})'], fontsize=8)

    yticks_bp = [1, 10, 100, 293, 1000, 10000, 100000]
    yticks_log = [np.log10(v) for v in yticks_bp]
    ax.set_yticks(yticks_log)
    ax.set_yticklabels([str(v) if v != 293 else '' for v in yticks_bp])
    ax.set_ylabel('Nearest methylation\ndistance (bp)', fontsize=9)
    ax.set_title('Shielded / Exposed dichotomy', fontsize=10, fontweight='bold')


def panel_d(ax, df_genes):
    """Panel D: Shielded vs Exposed expression variability (signed LFC violin)."""
    from scipy.stats import levene

    df = df_genes.dropna(subset=['LFC_T2vsT1', 'LFC_T3vsT1']).copy()
    exposed = df[df['is_exposed'] == 1]
    shielded = df[df['is_exposed'] == 0]

    # Colors
    col_shielded = '#999999'
    col_exposed = COL_EXPOSED  # purple #7E57C2 — unified with shared_utils

    # Positions: T2vsT1 at x=0,1; T3vsT1 at x=3,4
    comparisons = [
        (0, 1, 'LFC_T2vsT1', 'T2 vs T1\n(24 h vs 12 h)'),
        (3, 4, 'LFC_T3vsT1', 'T3 vs T1\n(50 h vs 12 h)'),
    ]

    for pos_s, pos_e, col_name, label in comparisons:
        data_s = shielded[col_name].values
        data_e = exposed[col_name].values

        # Violin - shielded
        vp_s = ax.violinplot([data_s], positions=[pos_s], widths=0.8,
                             showmedians=False, showextrema=False)
        for body in vp_s['bodies']:
            body.set_facecolor(col_shielded)
            body.set_alpha(0.5)

        # Violin - exposed
        vp_e = ax.violinplot([data_e], positions=[pos_e], widths=0.8,
                             showmedians=False, showextrema=False)
        for body in vp_e['bodies']:
            body.set_facecolor(col_exposed)
            body.set_alpha(0.5)

        # Box plot overlay for both groups
        for pos, data, col in [(pos_s, data_s, col_shielded),
                                (pos_e, data_e, col_exposed)]:
            q1, med, q3 = np.percentile(data, [25, 50, 75])
            iqr = q3 - q1
            whisker_lo = max(data[data >= q1 - 1.5 * iqr].min(), data.min())
            whisker_hi = min(data[data <= q3 + 1.5 * iqr].max(), data.max())
            # Whiskers
            ax.vlines(pos, whisker_lo, whisker_hi, color='black',
                      linewidth=0.8, zorder=5)
            # IQR box
            box_w = 0.25
            box = plt.Rectangle((pos - box_w / 2, q1), box_w, q3 - q1,
                                facecolor='white', edgecolor='black',
                                linewidth=0.8, zorder=6)
            ax.add_patch(box)
            # Median line
            ax.hlines(med, pos - box_w / 2, pos + box_w / 2,
                      color='black', linewidth=1.5, zorder=7)

        # Strip plot for both groups (unified design)
        rng = np.random.default_rng(42)
        # Shielded: subsample for readability (n=955 → 200 random points)
        n_show = min(200, len(data_s))
        idx_sub = rng.choice(len(data_s), n_show, replace=False)
        jitter_s = rng.uniform(-0.15, 0.15, n_show)
        ax.scatter(pos_s + jitter_s, data_s[idx_sub], s=4, alpha=0.25,
                   color=col_shielded, edgecolors='none', zorder=4)
        # Exposed: all points
        jitter_e = rng.uniform(-0.15, 0.15, len(data_e))
        ax.scatter(pos_e + jitter_e, data_e, s=6, alpha=0.5,
                   color=col_exposed, edgecolors='none', zorder=4)

        # Resampling Levene's test (size-matched permutation)
        n_exposed = len(data_e)
        n_perm = 1000
        rng_perm = np.random.default_rng(123)
        perm_stats = np.empty(n_perm)
        for i in range(n_perm):
            sub_s = rng_perm.choice(data_s, n_exposed, replace=False)
            perm_stats[i], _ = levene(data_e, sub_s)
        # Empirical p-value: fraction of resampled stats >= observed
        obs_stat, _ = levene(data_e, data_s)
        p_resamp = (np.sum(perm_stats >= obs_stat) + 1) / (n_perm + 1)

        # Variance ratio (effect size)
        var_ratio = np.var(data_e) / np.var(data_s)

        stars = ('***' if p_resamp < 0.001 else '**' if p_resamp < 0.01
                 else '*' if p_resamp < 0.05 else 'n.s.')
        y_absmax = max(np.abs(data_s).max(), np.abs(data_e).max())
        bracket_y = y_absmax * 1.10
        ax.plot([pos_s, pos_s, pos_e, pos_e],
                [bracket_y, bracket_y + 0.3, bracket_y + 0.3, bracket_y],
                color='black', linewidth=0.8)
        ax.text((pos_s + pos_e) / 2, bracket_y + 0.4,
                f'{stars}\nVR = {var_ratio:.1f}',
                ha='center', va='bottom', fontsize=7, fontweight='bold',
                linespacing=1.2)
        print(f'  {col_name}: Levene resampling p={p_resamp:.4f}, '
              f'VR={var_ratio:.2f} (Var_exposed/Var_shielded)')

    # Reference line at y=0
    ax.axhline(0, color='black', linewidth=0.5, linestyle='--', alpha=0.5)

    # X-axis labels
    ax.set_xticks([0.5, 3.5])
    ax.set_xticklabels([c[3] for c in comparisons], fontsize=8)
    ax.set_ylabel(r'log$_2$ fold change', fontsize=9)
    ax.set_xlim(-0.8, 5.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Legend (bottom-left to avoid bracket overlap)
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=col_shielded, alpha=0.5,
              label=f'Shielded (n={len(shielded)})'),
        Patch(facecolor=col_exposed, alpha=0.5,
              label=f'Exposed (n={len(exposed)})'),
    ]
    ax.legend(handles=legend_elements, loc='lower left', fontsize=7,
              framealpha=0.8)

    ax.set_title('Expression variability: Shielded vs Exposed TFs',
                 fontsize=10, fontweight='bold')


def main():
    apply_style()
    print('=== New Figure 2: Gatekeeper Model Overview ===')
    print()

    # Load data
    print('Loading data...')
    _, df_gccggc = load_geographic_redistribution()
    df_aagcccg = _load_aagcccg_sites()
    df_spatial = load_spatial_profile()
    # Load full TF dataset (1,055 genes, 57 Exposed, 998 Shielded; 293bp ROC threshold)
    _full_path = (BASE /
                  '11_epigenome_integration/analysis/52_shielded_exposed_boundary/tables' /
                  'all_genes_features_unified_n57.tsv')
    df_genes = pd.read_csv(_full_path, sep='\t')
    print(f'  Full TF dataset: {len(df_genes)} genes ({df_genes["is_exposed"].sum()} exposed)')


    # Layout: 2×2 with panel A full-width top
    fig = plt.figure(figsize=(mm_to_inch(174), mm_to_inch(240)))

    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 1.2],
                          hspace=0.45, wspace=0.40,
                          left=0.10, right=0.92, top=0.96, bottom=0.04)

    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])
    ax_d = fig.add_subplot(gs[2, :])

    print('Drawing Panel A: GCCGGC vs AAGCCCG redistribution...')
    panel_a(ax_a, df_gccggc, df_aagcccg)
    add_panel_label(ax_a, 'a', x=-0.05, y=1.12)

    print('Drawing Panel B: TSS protection zone...')
    panel_b(ax_b, df_spatial)
    add_panel_label(ax_b, 'b', x=-0.15, y=1.12)

    print('Drawing Panel C: Shielded/Exposed distribution...')
    panel_c(ax_c, df_genes)
    add_panel_label(ax_c, 'c', x=-0.15, y=1.12)

    print('Drawing Panel D: Shielded vs Exposed variability...')
    panel_d(ax_d, df_genes)
    add_panel_label(ax_d, 'd', x=-0.05, y=1.08)

    # Save
    out_path = FIG_DIR / 'Figure1_overview'
    save_figure(fig, out_path)
    print()
    print('=== Done ===')


if __name__ == '__main__':
    main()
