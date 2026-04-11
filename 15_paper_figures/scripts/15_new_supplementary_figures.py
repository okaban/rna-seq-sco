#!/usr/bin/env python3
"""
Generate new Supplementary Figures for the Gatekeeper Model paper.

S1: Reuse existing old figure (02+03 reclassification)
S2: Simpson's paradox — 6-panel detailed
S3: Regulatory avoidance forest plot (CMH)
S4: Sequence-level motif depletion — reuse existing analysis
S5: 62 exposed TF annotation heatmap
S6: TCS pair detailed — reuse existing
S7: Conservation metrics (H36)
S8: Protection zone characterization (moved from main Fig.3)

This script generates S2, S3, S4, S5, S6, S7, and S8.
S1 is a copy of the old main Figure 2 (done via file rename).
"""

import sys, importlib
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent))
_utils = importlib.import_module('00_shared_utils')
for _attr in dir(_utils):
    if not _attr.startswith('_'):
        globals()[_attr] = getattr(_utils, _attr)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats

FIG_SUP_DIR.mkdir(parents=True, exist_ok=True)
apply_style()

def _style(ax):
    """Remove top/right spines."""
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


# ==============================================================
# S2: Simpson's Paradox — 6-panel detailed
# ==============================================================
def make_fig_s2():
    print("Generating Fig S2: Simpson's Paradox...")

    fig = plt.figure(figsize=(mm_to_inch(180), mm_to_inch(220)))
    gs = gridspec.GridSpec(3, 2, hspace=0.45, wspace=0.35,
                           left=0.10, right=0.95, top=0.95, bottom=0.06)

    # Load data (returns tuple: transitions_df, geographic_df)
    gccggc, geo_gccggc = load_gene_methylation_transitions('GCCGGC')
    aagcccg, geo_aagcccg = load_gene_methylation_transitions('AAGCCCG')

    # --- Helper: violin plot for Lost vs Never ---
    def plot_lost_vs_never(ax, df, title, region_filter=None):
        d = df.copy()
        if region_filter and 'region' in d.columns:
            d = d[d['region'] == region_filter]

        if 'transition_T1T2' not in d.columns:
            # Try alternate column names
            for col in d.columns:
                if 'transition' in col.lower() or 'status' in col.lower():
                    d['transition_T1T2'] = d[col]
                    break

        lfc_col = None
        for c in ['LFC_T2vsT1', 'LFC', 'log2FC_T2']:
            if c in d.columns:
                lfc_col = c
                break
        if lfc_col is None:
            ax.text(0.5, 0.5, 'LFC column not found', transform=ax.transAxes, ha='center')
            ax.set_title(title, fontsize=9, fontweight='bold')
            return

        lost = d.loc[d['transition_T1T2'].str.contains('lost|Lost', na=False), lfc_col].dropna()
        never = d.loc[d['transition_T1T2'].str.contains('never|Never', na=False), lfc_col].dropna()

        if len(lost) < 5 or len(never) < 5:
            ax.text(0.5, 0.5, f'n_lost={len(lost)}, n_never={len(never)}\n(insufficient)',
                    transform=ax.transAxes, ha='center', fontsize=8)
            ax.set_title(title, fontsize=9, fontweight='bold')
            return

        parts = ax.violinplot([lost.values, never.values], positions=[1, 2], showmedians=True)
        for i, pc in enumerate(parts['bodies']):
            pc.set_facecolor(COL_ARTIFACT if i == 0 else COL_SIGNAL)
            pc.set_alpha(0.6)
        for key in ['cmins', 'cmaxes', 'cbars', 'cmedians']:
            if key in parts:
                parts[key].set_color(COL_DARK)

        stat, pval = stats.mannwhitneyu(lost, never, alternative='two-sided')
        ax.set_xticks([1, 2])
        ax.set_xticklabels([f'Lost\n(n={len(lost)})', f'Never\n(n={len(never)})'], fontsize=7)
        ax.set_ylabel('LFC (T2 vs T1)', fontsize=8)
        ax.set_title(title, fontsize=9, fontweight='bold')

        p_str = f'p = {pval:.2e}' if pval < 0.001 else f'p = {pval:.3f}'
        ax.text(0.5, 0.97, p_str, transform=ax.transAxes, ha='center', va='top', fontsize=8,
                color='red' if pval < 0.05 else COL_DARK)
        _style(ax)

    # Panel A: GCCGGC unstratified
    ax_a = fig.add_subplot(gs[0, 0])
    plot_lost_vs_never(ax_a, gccggc, 'A  GCCGGC — Unstratified')

    # Panel B: GCCGGC core-only
    ax_b = fig.add_subplot(gs[0, 1])
    plot_lost_vs_never(ax_b, gccggc, 'B  GCCGGC — Core only', region_filter='core')

    # Panel C: GCCGGC arm-only
    ax_c = fig.add_subplot(gs[1, 0])
    plot_lost_vs_never(ax_c, gccggc, 'C  GCCGGC — Arm only', region_filter='arm')

    # Panel D: AAGCCCG unstratified
    ax_d = fig.add_subplot(gs[1, 1])
    plot_lost_vs_never(ax_d, aagcccg, 'D  AAGCCCG — Unstratified')

    # Panel E: AAGCCCG core-only
    ax_e = fig.add_subplot(gs[2, 0])
    plot_lost_vs_never(ax_e, aagcccg, 'E  AAGCCCG — Core only', region_filter='core')

    # Panel F: AAGCCCG arm-only
    ax_f = fig.add_subplot(gs[2, 1])
    plot_lost_vs_never(ax_f, aagcccg, 'F  AAGCCCG — Arm only', region_filter='arm')

    save_figure(fig, FIG_SUP_DIR / 'FigS2_simpsons_paradox_detail')
    print("  Fig S2 saved.")


# ==============================================================
# S3: Regulatory avoidance forest plot (CMH)
# ==============================================================
def make_fig_s3():
    print("Generating Fig S3: Regulatory avoidance forest plot...")

    cmh = pd.read_csv(
        EPIGENOME / '43_regulatory_avoidance_geographic_test' / 'tables' / 'CMH_test_results.tsv',
        sep='\t')

    fig, ax = plt.subplots(figsize=(mm_to_inch(180), mm_to_inch(120)))

    # Filter to Regulatory/TF rows
    reg = cmh[cmh['category'].str.contains('Regulatory', na=False)].copy()
    if len(reg) == 0:
        reg = cmh.head(8)

    reg = reg.sort_values('CMH_OR').reset_index(drop=True)
    y_pos = np.arange(len(reg))

    # Plot forest
    ax.errorbar(reg['CMH_OR'], y_pos,
                xerr=[reg['CMH_OR'] - reg['CMH_OR_CI_low'],
                      reg['CMH_OR_CI_high'] - reg['CMH_OR']],
                fmt='o', color=COL_EXPOSED, markersize=6, capsize=3, linewidth=1.5)

    ax.axvline(1.0, color=COL_DARK, linestyle='--', linewidth=0.8, alpha=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(reg['motif'], fontsize=8)
    ax.set_xlabel('CMH-adjusted Odds Ratio', fontsize=9)
    ax.set_title('Regulatory gene methylation avoidance\n(CMH-adjusted for core/arm)', fontsize=10, fontweight='bold')

    # Add p-values
    for i, row in reg.iterrows():
        idx = reg.index.get_loc(i)
        p = row['CMH_p']
        p_str = f'p = {p:.1e}' if p < 0.001 else f'p = {p:.3f}'
        ax.text(max(reg['CMH_OR_CI_high']) * 1.05, idx, p_str, fontsize=7, va='center')

    ax.set_xlim(0, max(reg['CMH_OR_CI_high']) * 1.3)
    _style(ax)
    fig.tight_layout()
    save_figure(fig, FIG_SUP_DIR / 'FigS3_regulatory_avoidance_CMH')
    print("  Fig S3 saved.")


# ==============================================================
# S4: Sequence-level motif depletion (gene body + promoter)
# ==============================================================
def make_fig_s4():
    print("Generating Fig S4: Sequence-level motif depletion...")

    dep = pd.read_csv(
        EPIGENOME / '45_sequence_level_motif_depletion' / 'tables' / 'depletion_statistics.tsv',
        sep='\t')

    fig, axes = plt.subplots(1, 2, figsize=(mm_to_inch(180), mm_to_inch(100)),
                              gridspec_kw={'wspace': 0.35})

    for idx, (motif, color) in enumerate([('TGGCCGGC', COL_4mC), ('AAGCCCG', COL_6mA)]):
        ax = axes[idx]
        sub = dep[(dep['motif'] == motif) & (dep['region'] == 'combined')].copy()
        if len(sub) == 0:
            ax.text(0.5, 0.5, 'No data', transform=ax.transAxes, ha='center')
            continue

        zones = ['promoter', 'gene_body', 'extended_2kb']
        zone_labels = ['Promoter\n(±300 bp)', 'Gene body', 'Extended\n(±2 kb)']
        folds = []
        pvals = []
        for z in zones:
            row = sub[sub['zone'] == z]
            if len(row) > 0:
                folds.append(row.iloc[0]['mean_density_fold'])
                pvals.append(row.iloc[0]['wilcoxon_p'])
            else:
                folds.append(1.0)
                pvals.append(1.0)

        x = np.arange(len(zones))
        bars = ax.bar(x, folds, color=color, alpha=0.7, width=0.5, edgecolor='white')

        # Reference line at 1.0
        ax.axhline(1.0, color='gray', linewidth=0.8, linestyle='--')

        # P-value annotations
        for i, (f, p) in enumerate(zip(folds, pvals)):
            p_str = f'p = {p:.1e}' if p < 0.001 else f'p = {p:.3f}'
            sig_color = 'red' if p < 0.05 else COL_DARK
            ax.text(i, f - 0.03, p_str, ha='center', va='top', fontsize=6, color=sig_color)

        ax.set_xticks(x)
        ax.set_xticklabels(zone_labels, fontsize=7)
        ax.set_ylabel('Fold (regulatory / non-regulatory)', fontsize=8)
        ax.set_title(f'{chr(65 + idx)}  {motif} motif density', fontsize=9, fontweight='bold')
        ax.set_ylim(0.5, 1.2)
        _style(ax)

    fig.tight_layout()
    save_figure(fig, FIG_SUP_DIR / 'FigS4_sequence_motif_depletion')
    print("  Fig S4 saved.")


# ==============================================================
# S5: 62 exposed TF annotation summary
# ==============================================================
def make_fig_s5():
    print("Generating Fig S5: Exposed TF annotation...")

    tc = load_temporal_classification()
    exp = tc.copy()

    fig, axes = plt.subplots(1, 3, figsize=(mm_to_inch(180), mm_to_inch(200)),
                              gridspec_kw={'width_ratios': [1.5, 1, 1], 'wspace': 0.3})

    # Panel A: Waterfall plot sorted by LFC
    ax = axes[0]
    exp_sorted = exp.sort_values('LFC_T3vsT1')
    colors = [COL_ACTIVATION if b == 'activation' else COL_REPRESSION
              for b in exp_sorted['bloc']]
    y_pos = np.arange(len(exp_sorted))
    ax.barh(y_pos, exp_sorted['LFC_T3vsT1'], color=colors, height=0.8, edgecolor='none')
    ax.axvline(0, color=COL_DARK, linewidth=0.5)

    # Labels for top/bottom genes
    labels = []
    for _, row in exp_sorted.iterrows():
        name = row.get('gene_name', '')
        if pd.isna(name) or name == '':
            name = row.get('old_locus_tag', row.get('locus_tag', ''))
        labels.append(str(name))
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=4.5)
    ax.set_xlabel('LFC (T3 vs T1)', fontsize=8)
    ax.set_title('A  Expression change', fontsize=9, fontweight='bold')
    _style(ax)

    # Panel B: TF family bar chart
    ax = axes[1]
    fam_counts = exp['tf_family'].value_counts().head(10)
    y = np.arange(len(fam_counts))
    ax.barh(y, fam_counts.values, color=COL_EXPOSED, height=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(fam_counts.index, fontsize=7)
    ax.set_xlabel('Count', fontsize=8)
    ax.set_title('B  TF family', fontsize=9, fontweight='bold')
    _style(ax)

    # Panel C: Phase ratio distribution
    ax = axes[2]
    pr = exp['phase_ratio'].dropna()
    act_pr = exp.loc[exp['bloc'] == 'activation', 'phase_ratio'].dropna()
    rep_pr = exp.loc[exp['bloc'] == 'repression', 'phase_ratio'].dropna()
    bins = np.linspace(-0.5, 2.0, 20)
    ax.hist(act_pr, bins=bins, alpha=0.6, color=COL_ACTIVATION, label='Activation')
    ax.hist(rep_pr, bins=bins, alpha=0.6, color=COL_REPRESSION, label='Repression')
    ax.axvline(0.6, color=COL_DARK, linestyle='--', linewidth=0.8)
    ax.set_xlabel('Phase ratio', fontsize=8)
    ax.set_ylabel('Count', fontsize=8)
    ax.set_title('C  Timing', fontsize=9, fontweight='bold')
    ax.legend(fontsize=6, frameon=False)
    _style(ax)

    fig.tight_layout()
    save_figure(fig, FIG_SUP_DIR / 'FigS5_exposed_TF_annotation')
    print("  Fig S5 saved.")


# ==============================================================
# S6: TCS pair detailed analysis
# ==============================================================
def make_fig_s6():
    print("Generating Fig S6: TCS pair detailed analysis...")

    tcs = load_tcs_pairs()
    tc = load_temporal_classification()

    fig, axes = plt.subplots(1, 2, figsize=(mm_to_inch(180), mm_to_inch(120)),
                              gridspec_kw={'wspace': 0.40, 'width_ratios': [1.2, 1]})

    # Panel A: TCS pair diagram — exposed/shielded status with labels
    ax = axes[0]
    n_pairs = len(tcs)
    for i, (_, row) in enumerate(tcs.iterrows()):
        y = n_pairs - 1 - i

        # Sensor kinase
        sk_color = COL_EXPOSED if row['sk_exposed'] else COL_SHIELDED
        ax.scatter(0.3, y, color=sk_color, s=100, zorder=5,
                   edgecolors='black', linewidths=0.5, marker='D')
        ax.text(0.02, y, str(row['sk_old_locus']), ha='right', va='center', fontsize=7)

        # Response regulator
        rr_color = COL_EXPOSED if row['rr_exposed'] else COL_SHIELDED
        ax.scatter(0.7, y, color=rr_color, s=100, zorder=5,
                   edgecolors='black', linewidths=0.5, marker='o')
        ax.text(0.98, y, str(row['rr_old_locus']), ha='left', va='center', fontsize=7)

        # Connecting line
        ax.plot([0.35, 0.65], [y, y], color='#999', linewidth=0.8, zorder=2)

        # Distance annotation
        dist = row['distance_bp']
        ax.text(0.50, y + 0.25, f'{dist} bp', ha='center', va='bottom', fontsize=5.5,
                color='gray')

    ax.text(0.3, n_pairs - 0.3, 'SK', ha='center', va='bottom', fontsize=9, fontweight='bold')
    ax.text(0.7, n_pairs - 0.3, 'RR', ha='center', va='bottom', fontsize=9, fontweight='bold')

    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='D', color='w', markerfacecolor=COL_EXPOSED,
               markersize=8, label='Exposed', markeredgecolor='black', markeredgewidth=0.5),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=COL_SHIELDED,
               markersize=8, label='Shielded', markeredgecolor='black', markeredgewidth=0.5),
    ]
    ax.legend(handles=legend_elements, fontsize=7, loc='lower center', frameon=False, ncol=2)
    ax.set_xlim(-0.3, 1.3)
    ax.set_ylim(-1.2, n_pairs + 0.5)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title('A  7/7 TCS pairs:\none partner exposed', fontsize=9, fontweight='bold')

    # Panel B: Expression correlation of TCS pairs
    ax = axes[1]
    rhos = tcs['expression_rho'].values
    pvals = tcs['expression_pval'].values
    pair_labels = [f"{row['sk_old_locus']}\n{row['rr_old_locus']}" for _, row in tcs.iterrows()]

    y_pos = np.arange(n_pairs)
    colors_bar = ['#43A047' if p < 0.05 else '#B0BEC5' for p in pvals]
    ax.barh(y_pos, rhos, color=colors_bar, height=0.6, edgecolor='white')
    ax.axvline(0, color='gray', linewidth=0.5, linestyle=':')

    ax.set_yticks(y_pos)
    ax.set_yticklabels(pair_labels, fontsize=6)
    ax.set_xlabel('Spearman $\\rho$', fontsize=8)
    ax.set_title('B  SK-RR expression\ncorrelation', fontsize=9, fontweight='bold')
    ax.set_xlim(-1, 1)
    _style(ax)

    fig.tight_layout()
    save_figure(fig, FIG_SUP_DIR / 'FigS6_TCS_pair_analysis')
    print("  Fig S6 saved.")


# ==============================================================
# S7: Conservation metrics (H36)
# ==============================================================
def make_fig_s7():
    print("Generating Fig S7: Conservation metrics...")

    features = load_all_genes_features()
    exp = features[features['is_exposed'] == 1]
    shi = features[features['is_exposed'] == 0]

    fig, axes = plt.subplots(2, 2, figsize=(mm_to_inch(180), mm_to_inch(150)),
                              gridspec_kw={'hspace': 0.45, 'wspace': 0.35})

    # Conservation proxy data from H36 — use available columns
    metrics = [
        ('baseMean', 'Mean expression (baseMean)', 'log'),
        ('nearest_methyl_distance', 'Nearest methylation (bp)', 'log'),
        ('LFC_T3vsT1', 'LFC (T3 vs T1)', 'linear'),
        ('n_FIMO_hits', 'FIMO hits (known TF targets)', 'linear'),
    ]

    for idx, (col, label, scale) in enumerate(metrics):
        ax = axes[idx // 2, idx % 2]
        if col not in features.columns:
            ax.text(0.5, 0.5, f'{col} not available', transform=ax.transAxes, ha='center')
            continue

        e_vals = exp[col].dropna()
        s_vals = shi[col].dropna()

        if scale == 'log':
            e_vals = e_vals[e_vals > 0]
            s_vals = s_vals[s_vals > 0]
            e_vals = np.log10(e_vals)
            s_vals = np.log10(s_vals)
            label = f'log10({label})'

        parts = ax.violinplot([s_vals.values, e_vals.values], positions=[1, 2], showmedians=True)
        for i, pc in enumerate(parts['bodies']):
            pc.set_facecolor(COL_SHIELDED if i == 0 else COL_EXPOSED)
            pc.set_alpha(0.6)
        for key in ['cmins', 'cmaxes', 'cbars', 'cmedians']:
            if key in parts:
                parts[key].set_color(COL_DARK)

        stat, pval = stats.mannwhitneyu(e_vals, s_vals, alternative='two-sided')
        p_str = f'p = {pval:.2e}' if pval < 0.001 else f'p = {pval:.3f}'

        ax.set_xticks([1, 2])
        ax.set_xticklabels(['Shielded\n(n=993)', 'Exposed\n(n=62)'], fontsize=7)
        ax.set_ylabel(label, fontsize=8)
        panel_label = chr(65 + idx)
        ax.set_title(f'{panel_label}  {col}', fontsize=9, fontweight='bold')
        ax.text(0.5, 0.97, p_str, transform=ax.transAxes, ha='center', va='top', fontsize=8,
                color='red' if pval < 0.05 else COL_DARK)
        _style(ax)

    fig.tight_layout()
    save_figure(fig, FIG_SUP_DIR / 'FigS7_conservation_metrics')
    print("  Fig S7 saved.")


# ==============================================================
# S8: Protection Zone Characterization (moved from main Fig.3)
# ==============================================================
def make_fig_s8():
    print("Generating Fig S8: Protection zone characterization...")
    from sklearn.metrics import roc_curve, auc as sk_auc

    df_spatial = load_spatial_profile()
    df_genes = load_all_genes_features()
    df_roc_boundary = load_roc_analysis('boundary')
    df_roc_sequence = load_roc_analysis('sequence')
    df_tfbs = load_tfbs_spatial_profile()
    df_quintile = load_expression_quintile()

    fig = plt.figure(figsize=(mm_to_inch(180), mm_to_inch(185)))
    gs = gridspec.GridSpec(2, 12, hspace=0.55, wspace=1.2,
                           left=0.07, right=0.97, top=0.93, bottom=0.07,
                           height_ratios=[1, 1])

    ax_a = fig.add_subplot(gs[0, 0:4])
    ax_b = fig.add_subplot(gs[0, 5:8])
    ax_c = fig.add_subplot(gs[0, 9:12])
    ax_d = fig.add_subplot(gs[1, 0:6])
    ax_e = fig.add_subplot(gs[1, 6:12])

    # --- Panel A: TSS methylation gradient heatmap (all 3 timepoints) ---
    categories = [
        ('GCCGGC_4mC', 'GCCGGC (4mC)'),
        ('All_4mC', 'All 4mC'),
        ('AAGCCCG_6mA', 'AAGCCCG (6mA)'),
        ('All_6mA', 'All 6mA'),
        ('All_methylation', 'All methylation'),
    ]
    timepoints = ['T1', 'T2', 'T3']
    matrix_rows = []
    row_labels = []
    for tp in timepoints:
        for site_type, label in categories:
            sub = df_spatial[(df_spatial['site_type'] == site_type) &
                             (df_spatial['timepoint'] == tp) &
                             (df_spatial['gene_category'] == 'regulatory')]
            sub = sub.sort_values('bin_center_bp')
            if len(sub) == 0:
                continue
            matrix_rows.append(sub['density_per_kb_per_gene'].values)
            row_labels.append(f'{tp} {label}')

    if matrix_rows:
        matrix = np.array(matrix_rows)
        bin_centers = df_spatial[
            (df_spatial['site_type'] == 'All_methylation') &
            (df_spatial['timepoint'] == 'T1') &
            (df_spatial['gene_category'] == 'regulatory')
        ].sort_values('bin_center_bp')['bin_center_bp'].values

        norm_matrix = np.zeros_like(matrix)
        for i in range(len(matrix)):
            smoothed = np.convolve(matrix[i], np.ones(5) / 5, mode='same')
            flank_mask = (np.abs(bin_centers) >= 3000) & (np.abs(bin_centers) <= 5000)
            flank_mean = smoothed[flank_mask].mean() if flank_mask.sum() > 0 else smoothed.mean()
            norm_matrix[i] = smoothed / flank_mean if flank_mean > 0 else 1.0
        norm_matrix = np.clip(norm_matrix, 0.5, 1.5)

        im = ax_a.imshow(norm_matrix, aspect='auto', cmap='RdBu_r',
                         vmin=0.6, vmax=1.4,
                         extent=[bin_centers[0]/1000, bin_centers[-1]/1000,
                                 len(row_labels)-0.5, -0.5])
        ax_a.set_yticks(range(len(row_labels)))
        ax_a.set_yticklabels(row_labels, fontsize=5)
        ax_a.set_xlabel('Distance from TSS (kb)', fontsize=8)

        # T1/T2/T3 boundary lines
        n_cat = len(categories)
        for boundary in range(1, len(timepoints)):
            ax_a.axhline(boundary * n_cat - 0.5, color='white', linewidth=1.5)

        ax_a.axvline(-1.3, color='white', linewidth=0.8, linestyle='--')
        ax_a.axvline(0.7, color='white', linewidth=0.8, linestyle='--')
        ax_a.axvline(0, color='white', linewidth=0.5, linestyle=':')

        cbar = plt.colorbar(im, ax=ax_a, shrink=0.6, pad=0.03, aspect=15)
        cbar.set_label('Fold vs flanking', fontsize=6, labelpad=2)
        cbar.ax.tick_params(labelsize=5)
    ax_a.set_title('A  Methylation at TSS\n(regulatory, all timepoints)',
                    fontsize=8, fontweight='bold')

    # --- Panel B: ROC curves ---
    y_true = df_genes['is_exposed'].values
    valid = ~np.isnan(df_genes['nearest_methyl_distance'].values)
    dist = df_genes.loc[valid, 'nearest_methyl_distance'].values
    y_val = y_true[valid]
    fpr_d, tpr_d, _ = roc_curve(y_val, -dist)
    auc_d = sk_auc(fpr_d, tpr_d)
    bm = df_genes.loc[valid, 'baseMean'].values
    fpr_b, tpr_b, _ = roc_curve(y_val, -bm)
    auc_b = sk_auc(fpr_b, tpr_b)

    ax_b.plot(fpr_d, tpr_d, color=COL_EXPOSED, linewidth=2,
              label=f'Nearest dist.\nAUC = {auc_d:.3f}')
    ax_b.plot(fpr_b, tpr_b, color=COL_GRAY, linewidth=1.5, linestyle='--',
              label=f'Expression\nAUC = {auc_b:.3f}')
    ax_b.plot([0, 1], [0, 1], 'k:', linewidth=0.5, alpha=0.5)
    ax_b.set_xlabel('FPR', fontsize=8)
    ax_b.set_ylabel('TPR', fontsize=8)
    ax_b.set_title('B  Exposed classification', fontsize=8, fontweight='bold')
    ax_b.legend(fontsize=6, loc='lower right', frameon=True, edgecolor='#ccc')
    ax_b.set_aspect('equal')
    _style(ax_b)

    # --- Panel C: Feature discriminative power ---
    features_list = [
        ('Distance\n(observed)', df_roc_boundary[df_roc_boundary['feature'] == 'nearest_methyl_distance'].iloc[0]),
        ('Sequence\n(5-fold CV)', df_roc_sequence[df_roc_sequence['feature'] == 'combined_LR_sequence_CV'].iloc[0]),
        ('GC%\n(±300 bp)', df_roc_sequence[df_roc_sequence['feature'] == 'GC_300bp'].iloc[0]),
        ('Expression', df_roc_boundary[df_roc_boundary['feature'] == 'baseMean'].iloc[0]),
    ]
    labels_c = [f[0] for f in features_list]
    aucs_c = [f[1]['AUC'] for f in features_list]
    ci_lo = [f[1].get('AUC_CI_lo', f[1].get('AUC_95CI_lower', 0)) for f in features_list]
    ci_hi = [f[1].get('AUC_CI_hi', f[1].get('AUC_95CI_upper', 0)) for f in features_list]
    errors_c = [[a - lo for a, lo in zip(aucs_c, ci_lo)],
                [hi - a for a, hi in zip(aucs_c, ci_hi)]]
    colors_c = [COL_EXPOSED, COL_SIGNAL, COL_GRAY, COL_GRAY]

    ax_c.bar(range(len(labels_c)), aucs_c, color=colors_c, alpha=0.8, width=0.6)
    ax_c.errorbar(range(len(labels_c)), aucs_c, yerr=errors_c,
                  fmt='none', ecolor='#333', capsize=3, linewidth=1)
    for i, v in enumerate(aucs_c):
        ax_c.text(i, v + errors_c[1][i] + 0.02, f'{v:.3f}',
                  ha='center', va='bottom', fontsize=6, fontweight='bold')
    ax_c.axhline(0.5, color='gray', linewidth=0.5, linestyle=':')
    ax_c.set_xticks(range(len(labels_c)))
    ax_c.set_xticklabels(labels_c, fontsize=6)
    ax_c.set_ylabel('AUC', fontsize=8)
    ax_c.set_ylim(0.35, 1.05)
    ax_c.set_title('C  Feature comparison', fontsize=8, fontweight='bold')
    _style(ax_c)

    # --- Panel D: TFBS not protective ---
    for mtype, label, color in [('GCCGGC_T1', 'GCCGGC', COL_4mC),
                                 ('All_T1', 'All methylation', COL_DARK)]:
        sub = df_tfbs[df_tfbs['methylation_type'] == mtype].sort_values('bin_center_bp')
        if len(sub) == 0:
            continue
        ax_d.plot(sub['bin_center_bp'], sub['obs_over_random'],
                  color=color, linewidth=1.5, label=label)
    ax_d.axhline(1.0, color='gray', linewidth=0.8, linestyle='--')
    ax_d.set_xlabel('Distance from TFBS center (bp)', fontsize=8)
    ax_d.set_ylabel('Observed / Random', fontsize=8)
    ax_d.set_title('D  TFBS methylation profile', fontsize=8, fontweight='bold')
    ax_d.legend(fontsize=6, loc='lower left', frameon=False)
    _style(ax_d)

    # --- Panel E: Expression-independence ---
    q = df_quintile['quintile'].values
    frac = df_quintile['frac_exposed'].values * 100
    ax_e.bar(q, frac, color=COL_EXPOSED, alpha=0.7, width=0.6)
    for qi, fi in zip(q, frac):
        ax_e.text(qi, fi + 0.3, f'{fi:.1f}%', ha='center', va='bottom', fontsize=6)
    slope, intercept, r, p, se = stats.linregress(q, frac)
    x_fit = np.array([0.5, 5.5])
    ax_e.plot(x_fit, slope * x_fit + intercept, 'k--', linewidth=0.8, alpha=0.5)
    ax_e.text(0.97, 0.95, f'baseMean AUC = 0.547\nJT p = 0.730',
              transform=ax_e.transAxes, ha='right', va='top', fontsize=6,
              bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                        edgecolor='#ccc', alpha=0.9))
    ax_e.set_xlabel('Expression quintile\n(1 = lowest, 5 = highest)', fontsize=8)
    ax_e.set_ylabel('% Exposed', fontsize=8)
    ax_e.set_ylim(0, max(frac) * 1.4)
    ax_e.set_xticks(q)
    ax_e.set_title('E  Expression-independence', fontsize=8, fontweight='bold')
    _style(ax_e)

    save_figure(fig, FIG_SUP_DIR / 'FigS8_protection_zone')
    print("  Fig S8 saved.")


# ==============================================================
# Main
# ==============================================================
if __name__ == '__main__':
    make_fig_s2()
    make_fig_s3()
    make_fig_s4()
    make_fig_s5()
    make_fig_s6()
    make_fig_s7()
    make_fig_s8()
    print(f"\nAll supplementary figures saved to {FIG_SUP_DIR}")
