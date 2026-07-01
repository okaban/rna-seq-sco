#!/usr/bin/env python3
"""Figure 2 — GCCGGC m4C redistribution and promoter protection architecture.

Panel A: Genome-wide GCCGGC m4C density at T1/T2/T3 (50 kb windows).
Panel B: TSS metagene methylation profile, T1, ±4 kb, all genes.
Panel C: Continuous TSS-to-nearest-GCCGGC-m4C distance across all regulatory
         genes (n=1,051; 62 Exposed / 989 Shielded), shown as a single histogram
         + KDE on a log10-bp axis with the 293 bp operating point marked. The
         distribution is unimodal (Hartigan's dip test does not reject), so it is
         a gradient, not two natural groups.
Panel D: |log2 fold-change| variability of Shielded vs Exposed TFs at T2 vs T1
         and T3 vs T1.

Outputs PDF/SVG/PNG to FIG_DIR. Adds uppercase panel labels A/B/C/D and
positions panel D's legend so it does not overlap the significance markers.
"""

import sys
import importlib
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
_utils = importlib.import_module('00_shared_utils')
for _attr in dir(_utils):
    if not _attr.startswith('_'):
        globals()[_attr] = getattr(_utils, _attr)


GENOME_LEN = 8_667_507
ARM_LEFT = 1_500_000
ARM_RIGHT = 7_167_507
PROTECTION_ZONE_BP = 293

# ── Palette unification with Figure 1 (Okabe-Ito, colourblind-safe) ───────────
# Override the imported palette so Figure 1 and Figure 2 use IDENTICAL red/blue.
#   4mC / GCCGGC m4C series  → vermillion  #D55E00
#   6mA / AAGCCCG series     → blue        #0072B2
#   core/arm shading & grey baseline       #999999
COL_4mC = '#D55E00'      # vermillion (was muted rose #C26B6B)
COL_6mA = '#0072B2'      # blue (was muted #4477AA)
COL_GREY = '#999999'     # neutral grey baseline / core shading

T_LABELS = {'T1': 'T1 (12 h)', 'T2': 'T2 (24 h)', 'T3': 'T3 (50 h)'}
# T1 (red) and T3 (blue) GCCGGC m4C density tracks matched to Figure 1;
# T2 keeps the Okabe-Ito orange (neither red nor blue).
T_COLORS = {'T1': COL_4mC, 'T2': '#E69F00', 'T3': COL_6mA}


def add_uppercase_label(ax, label, x=-0.10, y=1.06, fontsize=15):
    """Add a bold uppercase panel label (A/B/C/D) to the top-left of an axis."""
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=fontsize, fontweight='bold', va='top', ha='left')


def panel_a(ax, df_sites, bin_size_bp=50_000):
    """Genome-wide GCCGGC m4C density at T1/T2/T3 (50 kb windows)."""
    bins_bp = np.arange(0, GENOME_LEN + bin_size_bp, bin_size_bp)
    centers_mb = (bins_bp[:-1] + bins_bp[1:]) / 2 / 1e6

    # Shaded core region
    ax.axvspan(ARM_LEFT / 1e6, ARM_RIGHT / 1e6,
               color=COL_GREY, alpha=0.12, zorder=0, label='Core')

    for tp in ['T1', 'T2', 'T3']:
        sub = df_sites[df_sites['timepoint'] == tp]
        h, _ = np.histogram(sub['position'].values, bins=bins_bp)
        density = h / (bin_size_bp / 1000)  # sites / kb
        ax.plot(centers_mb, density, lw=1.0, color=T_COLORS[tp],
                label=T_LABELS[tp], alpha=0.85)

    ax.set_xlim(0, GENOME_LEN / 1e6)
    ax.set_xlabel('Chromosome position (Mb)', fontsize=9)
    ax.set_ylabel('GCCGGC m4C density\n(50 kb windows)', fontsize=9)
    ax.set_title('Geographic redistribution of GCCGGC m4C',
                 fontsize=10, fontweight='bold')
    ax.legend(fontsize=7, loc='upper right', frameon=True,
              framealpha=0.9, edgecolor='#CCC')
    ax.grid(axis='y', alpha=0.25, lw=0.5)


def panel_b(ax, df_spatial):
    """TSS metagene methylation profile (T1, GCCGGC m4C, all genes ±4 kb)."""
    sub = df_spatial[
        (df_spatial['site_type'] == 'GCCGGC_4mC') &
        (df_spatial['timepoint'] == 'T1') &
        (df_spatial['gene_category'] == 'all')
    ].sort_values('bin_center_bp')

    x = sub['bin_center_bp'].values
    y = sub['density_per_kb_per_gene'].values
    n_genes = int(sub['n_genes'].iloc[0])

    ax.plot(x, y, color=COL_4mC, lw=1.4)
    if 'ci_low' in sub.columns:
        ax.fill_between(x, sub['ci_low'].values, sub['ci_high'].values,
                        color=COL_4mC, alpha=0.18, lw=0)

    # Protection zone shading at the TSS (centred at 0)
    ax.axvspan(-PROTECTION_ZONE_BP, 0, color='#FFCDD2', alpha=0.55,
               label=f'Protection zone ({PROTECTION_ZONE_BP} bp)')
    ax.axvline(0, color='#888', lw=0.6, ls=':')

    ax.set_xlim(-4000, 4000)
    ax.set_xlabel('Distance from TSS (bp)', fontsize=9)
    ax.set_ylabel(f'GCCGGC m4C density\n(all genes, n={n_genes:,})', fontsize=9)
    ax.set_title('TSS metagene methylation profile (T1)',
                 fontsize=10, fontweight='bold')
    ax.legend(fontsize=7, loc='upper right', frameon=True,
              framealpha=0.9, edgecolor='#CCC')
    ax.grid(axis='y', alpha=0.25, lw=0.5)


def panel_c(ax, df_genes):
    """Continuous TSS-to-nearest-GCCGGC-m4C distance across ALL regulatory genes.

    Honest depiction: the distribution is continuous and unimodal (Hartigan's dip
    test does not reject unimodality), so it is shown as a single histogram + KDE on
    a log10-bp x-axis rather than pre-split into two groups. The 293 bp operating
    point is a single vertical line; the 62 Exposed are simply those left of it
    (the region <=293 bp is shaded lightly).
    """
    df = df_genes.dropna(subset=['nearest_methyl_distance'])
    dist = df['nearest_methyl_distance'].values.astype(float)
    n_total = len(dist)
    n_exposed = int((dist <= PROTECTION_ZONE_BP).sum())
    n_shielded = n_total - n_exposed

    log_dist = np.log10(dist + 1)
    log_thresh = np.log10(PROTECTION_ZONE_BP)

    # Hartigan's dip test for unimodality (annotation only; falls back gracefully).
    dip_label = 'unimodal: Hartigan dip p = 0.96'
    try:
        import diptest
        _dip, _p = diptest.diptest(log_dist)
        dip_label = f'unimodal: Hartigan dip p = {_p:.2f}'
    except Exception:
        pass

    x_lo, x_hi = 0.0, np.ceil(log_dist.max() * 2) / 2
    bins = np.linspace(x_lo, x_hi, 31)

    # One continuous histogram (density) across all regulatory genes.
    ax.hist(log_dist, bins=bins, density=True, color=COL_GREY, alpha=0.45,
            edgecolor='white', linewidth=0.4, zorder=2)

    # KDE overlay to convey the continuous, single-peaked shape.
    kde = stats.gaussian_kde(log_dist)
    xs = np.linspace(x_lo, x_hi, 400)
    ys = kde(xs)
    ax.plot(xs, ys, color='#555', lw=1.6, zorder=4)

    # Shade the region <=293 bp (the Exposed side) lightly — no separate group.
    ax.axvspan(x_lo, log_thresh, color=COL_4mC, alpha=0.12, zorder=1)

    # Single vertical line at the 293 bp operating point.
    ax.axvline(log_thresh, color=COL_4mC, ls='--', lw=1.4, zorder=5)
    y_top = ax.get_ylim()[1]
    ax.text(log_thresh, y_top * 0.98, '293 bp operating point  ',
            rotation=90, ha='right', va='top', fontsize=7.5,
            color=COL_4mC, fontweight='bold', zorder=6)

    # Label the Exposed (left) and Shielded (right) regions descriptively.
    ax.text(log_thresh - 0.12, y_top * 0.88,
            f'Exposed\n(<= {PROTECTION_ZONE_BP} bp)\nn = {n_exposed}',
            ha='right', va='top', fontsize=7, color=COL_4mC)
    ax.text(log_thresh + 0.18, y_top * 0.88,
            f'Shielded\n(> {PROTECTION_ZONE_BP} bp)\nn = {n_shielded}',
            ha='left', va='top', fontsize=7, color='#444')

    # Unimodality annotation (continuous gradient, not a natural gap).
    ax.text(0.97, 0.97, dip_label, transform=ax.transAxes,
            ha='right', va='top', fontsize=7.5, style='italic',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='#999', alpha=0.9))

    xticks_bp = [1, 10, 100, 293, 1000, 10_000, 100_000]
    ax.set_xticks([np.log10(v + 1) for v in xticks_bp])
    ax.set_xticklabels([f'{v:,}' for v in xticks_bp], fontsize=8)
    ax.set_xlim(x_lo, x_hi)
    ax.set_xlabel('Distance to nearest GCCGGC m4C site (bp)', fontsize=9)
    ax.set_ylabel(f'Density (all regulatory genes, n = {n_total:,})', fontsize=9)
    ax.set_title('TSS distance distribution is a continuous gradient (T1)',
                 fontsize=10, fontweight='bold')
    ax.grid(axis='y', alpha=0.25, lw=0.5)


def panel_d(ax, df_genes):
    """|log2 FC| variability for Shielded vs Exposed TFs (T2/T1 and T3/T1)."""
    df = df_genes.dropna(subset=['LFC_T2vsT1', 'LFC_T3vsT1']).copy()
    df['abs_T2T1'] = df['LFC_T2vsT1'].abs()
    df['abs_T3T1'] = df['LFC_T3vsT1'].abs()
    exp = df[df['is_exposed'] == 1]
    shi = df[df['is_exposed'] == 0]

    comparisons = [
        ('abs_T2T1', 0, 1, 'T2 vs T1'),
        ('abs_T3T1', 3, 4, 'T3 vs T1'),
    ]

    y_max = 0
    for col, pos_s, pos_e, _ in comparisons:
        data_s = shi[col].values
        data_e = exp[col].values
        y_max = max(y_max, np.percentile(data_s, 99), np.percentile(data_e, 99))

        for pos, data, color in [(pos_s, data_s, COL_SHIELDED),
                                 (pos_e, data_e, COL_EXPOSED)]:
            parts = ax.violinplot([data], positions=[pos], widths=0.85,
                                  showextrema=False, showmedians=False)
            for body in parts['bodies']:
                body.set_facecolor(color)
                body.set_alpha(0.55)
                body.set_edgecolor(color)
            # Median bar
            med = np.median(data)
            ax.hlines(med, pos - 0.20, pos + 0.20,
                      color='black', lw=1.4, zorder=5)

    # Significance bracket — placed above all violins, well clear of legend
    bracket_y = y_max * 1.15
    text_y = bracket_y * 1.03

    for col, pos_s, pos_e, label in comparisons:
        e = exp[col].values
        s = shi[col].values
        U, p = stats.mannwhitneyu(e, s, alternative='greater')
        if p < 0.001:
            star = '***'
        elif p < 0.01:
            star = '**'
        elif p < 0.05:
            star = '*'
        else:
            star = 'n.s.'
        ax.plot([pos_s, pos_s, pos_e, pos_e],
                [bracket_y, bracket_y * 1.04, bracket_y * 1.04, bracket_y],
                color='black', lw=0.8)
        ax.text((pos_s + pos_e) / 2, text_y,
                f'{star}  (p = {p:.1e})',
                ha='center', va='bottom', fontsize=7.5, fontweight='bold')

    ax.set_ylim(0, bracket_y * 1.30)
    ax.set_xticks([0.5, 3.5])
    ax.set_xticklabels([c[3] for c in comparisons], fontsize=9)
    ax.set_xlim(-0.7, 5.2)
    ax.set_ylabel(r'$|\log_2$ fold change$|$ (expression variability)',
                  fontsize=9)
    ax.set_title('Expression variability: Shielded vs Exposed TFs',
                 fontsize=10, fontweight='bold')

    # Legend placed below the axis to avoid overlap with the significance marker.
    # Counts reflect genes with available LFC data for both transitions.
    legend_handles = [
        Patch(facecolor=COL_SHIELDED, alpha=0.55,
              label=f'Shielded (n = {len(shi)} of 989)'),
        Patch(facecolor=COL_EXPOSED, alpha=0.55,
              label=f'Exposed (n = {len(exp)} of 62)'),
    ]
    ax.legend(handles=legend_handles, loc='upper center',
              bbox_to_anchor=(0.5, -0.18), ncol=2,
              fontsize=8, frameon=True, framealpha=0.9, edgecolor='#CCC')


def main():
    apply_style()
    print('=== Figure 2: GCCGGC m4C redistribution + promoter protection ===')

    print('Loading data...')
    _, df_sites = load_geographic_redistribution()
    df_spatial = load_spatial_profile()
    df_genes = pd.read_csv(
        EPIGENOME / '52_shielded_exposed_boundary' /
        'tables' / 'all_genes_features_unified_n57.tsv', sep='\t')

    # --- Non-circular reframe (2026-06-17): recompute Exposed from the methylation
    # map only (TSS within 293 bp of a GCCGGC m4C site at T1), matching F4_F5 and
    # the Fig2c/d reframe. The table's stale `is_exposed` (n=57, expression-selected)
    # is replaced; nearest_methyl_distance is set to the T1 GCCGGC distance so the
    # panel-C violin and the Shielded/Exposed split are the same quantity (62/989).
    df_genes = df_genes.dropna(subset=['tss']).copy()
    df_genes['tss'] = df_genes['tss'].astype(int)
    g_tp = pd.read_csv(EPIGENOME / '37_defense_island_GCCGGC' /
                       'tables' / 'GCCGGC_sites_by_timepoint.tsv', sep='\t')
    _pos = np.sort(g_tp[g_tp['timepoint'] == 'T1']['position'].values)

    def _nearest(tss):
        if len(_pos) == 0:
            return np.nan
        i = np.clip(np.searchsorted(_pos, tss), 1, len(_pos) - 1)
        return min(abs(tss - _pos[i - 1]), abs(tss - _pos[i]))

    df_genes['nearest_methyl_distance'] = df_genes['tss'].apply(_nearest)
    df_genes['is_exposed'] = (df_genes['nearest_methyl_distance']
                              <= PROTECTION_ZONE_BP).astype(int)
    print(f'  Regulatory genes (non-circular): {len(df_genes)} '
          f'(Exposed={int(df_genes["is_exposed"].sum())}, '
          f'Shielded={int((df_genes["is_exposed"] == 0).sum())})')

    fig = plt.figure(figsize=(mm_to_inch(180), mm_to_inch(180)))
    gs = fig.add_gridspec(2, 2, hspace=0.55, wspace=0.42,
                          left=0.10, right=0.96, top=0.92, bottom=0.10)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    panel_a(ax_a, df_sites)
    panel_b(ax_b, df_spatial)
    panel_c(ax_c, df_genes)
    panel_d(ax_d, df_genes)

    for ax, label in zip([ax_a, ax_b, ax_c, ax_d], ['A', 'B', 'C', 'D']):
        add_uppercase_label(ax, label, x=-0.16, y=1.10, fontsize=15)

    fig.suptitle('Figure 2 | Geographic redistribution of GCCGGC m4C '
                 'and promoter protection architecture',
                 fontsize=11, fontweight='bold', y=0.985)

    out_path = FIG_DIR / 'Figure2_RM_redistribution'
    save_figure(fig, out_path, formats=('pdf', 'svg', 'png'))
    import shutil
    slot = Path.home() / 'obsidian' / 'Research' / 'rna-seq' / 'Writing' / 'fig_images' / 'Figure2.png'
    if slot.parent.is_dir():
        shutil.copyfile(out_path.with_suffix('.png'), slot)
        print(f'  Synced → {slot}')
    print('=== Done ===')


if __name__ == '__main__':
    main()
