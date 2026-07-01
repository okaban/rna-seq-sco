#!/usr/bin/env python3
"""
New Figure 3: Promoter protection zone architecture and Shielded/Exposed classification.

[2026-06-17 non-circular reframe] The former ROC/feature-AUC panels were REMOVED:
under the non-circular definition the Exposed label is itself a function of
promoter-proximal distance, so a distance classifier is tautological and no AUC is
reported. Three panels remain, matching the manuscript legend:

(a) TSS methylation gradient heatmap (site_type × distance bins)
(b) TFBS-centred methylation profile (individual TFBS is not protective)
(c) Expression-independence (Exposed fraction by expression quintile; JT trend n.s.)
"""

import importlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D

_utils = importlib.import_module('00_shared_utils')
for _attr in dir(_utils):
    if not _attr.startswith('_'):
        globals()[_attr] = getattr(_utils, _attr)


def panel_a(ax, df_spatial):
    """Panel A: TSS methylation gradient heatmap — site type × distance."""
    # Show regulatory genes, T1, different methylation types
    # NOTE: AAGCCCG_6mA is omitted from this per-gene metagene. With only ~260
    # AAGCCCG 6mA sites across all regulatory genes, the per-gene density profile
    # is dominated by the minority of promoter-marked regulators and shows a
    # spurious near-TSS 'enrichment' that contradicts the statistically-correct
    # depletion (CMH-adjusted OR = 0.409; reported in the text). 'All 6mA'
    # (n=1,934) is well-sampled and correctly shows promoter depletion.
    categories = [
        ('GCCGGC_4mC', 'GCCGGC (4mC)'),
        ('All_4mC', 'All 4mC'),
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


def panel_b(ax, df_tfbs):
    """Panel B: Individual TF BS is NOT protective — TFBS-centered profile."""
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


def panel_c(ax, df_quintile):
    """Panel C: Expression-independence — Exposed fraction by expression quintile."""
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

    # Trend-test only — no classifier AUC (tautological under the non-circular definition)
    ax.text(0.97, 0.95, 'Jonckheere–Terpstra\np = 0.730 (n.s.)',
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

    # Load data (no ROC/classifier data — AUC removed under the non-circular reframe)
    print('Loading data...')
    df_spatial = load_spatial_profile()
    df_tfbs = load_tfbs_spatial_profile()
    df_quintile = load_expression_quintile()

    # Create figure: 180mm × 150mm. Row 1 = (a) heatmap (full width);
    # Row 2 = (b) TFBS profile + (c) expression-independence.
    fig = plt.figure(figsize=(mm_to_inch(180), mm_to_inch(150)))

    gs = fig.add_gridspec(2, 12, hspace=0.55, wspace=1.4,
                          left=0.08, right=0.95, top=0.92, bottom=0.10,
                          height_ratios=[1.05, 1])

    ax_a = fig.add_subplot(gs[0, 1:11])     # (a) heatmap, centred full width
    ax_b = fig.add_subplot(gs[1, 0:6])      # (b) TFBS metagene
    ax_c = fig.add_subplot(gs[1, 6:12])     # (c) expression-independence

    print('Drawing Panel A: TSS methylation gradient heatmap...')
    panel_a(ax_a, df_spatial)
    add_panel_label(ax_a, 'a', x=-0.10, y=1.12)

    print('Drawing Panel B: TFBS not protective...')
    panel_b(ax_b, df_tfbs)
    add_panel_label(ax_b, 'b', x=-0.12, y=1.12)

    print('Drawing Panel C: Expression-independence...')
    panel_c(ax_c, df_quintile)
    add_panel_label(ax_c, 'c', x=-0.12, y=1.12)

    # Save (and sync into the Obsidian manuscript slot Figure3.png)
    out_path = FIG_DIR / 'Figure3_protection_zone'
    save_figure(fig, out_path, formats=('pdf', 'svg', 'png'))
    import shutil
    slot = Path.home() / 'obsidian' / 'Research' / 'rna-seq' / 'Writing' / 'fig_images' / 'Figure3.png'
    if slot.parent.is_dir():
        shutil.copyfile(out_path.with_suffix('.png'), slot)
        print(f'  Synced → {slot}')
    print()
    print('=== Done ===')


if __name__ == '__main__':
    main()
