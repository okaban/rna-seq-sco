#!/usr/bin/env python3
"""
Figure 5: Temporal Dynamics of DNA Methylation
  Panel A: HC site counts by timepoint and modification type
  Panel B: Motif-specific temporal loss patterns (heatmap)
  Panel C: Genome browser tracks for representative loci
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import importlib
_utils = importlib.import_module('00_shared_utils')
for _attr in dir(_utils):
    if not _attr.startswith('_'):
        globals()[_attr] = getattr(_utils, _attr)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
from Bio import SeqIO


def panel_a_temporal_counts(ax):
    """HC methylation site counts by timepoint and modification type."""
    df_4mc, df_6ma = load_methylation_census()

    counts_4mc = df_4mc.groupby('timepoint').size().reindex(['T1', 'T2', 'T3'], fill_value=0)
    counts_6ma = df_6ma.groupby('timepoint').size().reindex(['T1', 'T2', 'T3'], fill_value=0)

    x = np.arange(3)
    width = 0.35

    bars1 = ax.bar(x - width/2, counts_4mc.values, width, color=COL_4mC,
                   edgecolor='white', linewidth=0.5, label='4mC', alpha=0.9,
                   zorder=3)
    bars2 = ax.bar(x + width/2, counts_6ma.values, width, color=COL_6mA,
                   edgecolor='white', linewidth=0.5, label='6mA', alpha=0.9,
                   zorder=3)

    # Value labels
    for bars in [bars1, bars2]:
        for bar in bars:
            h = bar.get_height()
            if h > 0:
                ax.text(bar.get_x() + bar.get_width()/2, h + 30,
                        f'{int(h):,}', ha='center', va='bottom', fontsize=7.5,
                        fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels(['T1 (24 h)', 'T2 (48 h)', 'T3 (72 h)'])
    ax.set_ylabel('HC methylation sites')
    ax.set_title('Temporal distribution', fontsize=11, fontweight='bold')
    ax.legend(fontsize=8, frameon=False, loc='upper right')
    ax.set_ylim(0, max(counts_4mc.max(), counts_6ma.max()) * 1.15)
    ax.grid(axis='y', alpha=0.3, lw=0.5)

    # Percentage loss annotation
    loss_4mc = (1 - counts_4mc['T3'] / counts_4mc['T1']) * 100
    loss_6ma = (1 - counts_6ma['T3'] / counts_6ma['T1']) * 100
    ax.text(0.02, 0.95,
            f'T1\u2192T3 loss:\n4mC: {loss_4mc:.0f}%\n6mA: {loss_6ma:.0f}%',
            transform=ax.transAxes, fontsize=7.5, va='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#FFF3E0',
                      edgecolor='#FB8C00', alpha=0.9, linewidth=0.8))


def panel_b_motif_temporal(ax):
    """Motif-specific temporal patterns as heatmap."""
    df_4mc, df_6ma = load_methylation_census()
    df_all = pd.concat([df_4mc, df_6ma])

    # Group CCGG family
    ccgg_fam = ['TGGCCGGC', 'GGCCGG', 'CCGG (other)', 'CCGG']
    df_all['motif_group'] = df_all['final_motif'].apply(
        lambda x: 'CCGG family' if x in ccgg_fam else
                  ('Unassigned' if x == 'unassigned' else x))

    # Keep only main motifs
    main_motifs = ['CCGG family', 'AAGCCCG', 'GATC', 'CCGKCA', 'Unassigned']
    df_main = df_all[df_all['motif_group'].isin(main_motifs)].copy()

    # Counts per motif × timepoint
    ct = df_main.groupby(['motif_group', 'timepoint']).size().unstack(fill_value=0)
    ct = ct.reindex(columns=['T1', 'T2', 'T3'], fill_value=0)
    ct = ct.reindex(main_motifs)

    # Normalize: percentage of T1 count (T1 = 100%)
    ct_pct = ct.div(ct['T1'], axis=0) * 100
    ct_pct = ct_pct.fillna(0)

    # Custom colormap: red (low) -> white (mid) -> blue (high)
    cmap = LinearSegmentedColormap.from_list('temporal',
        ['#FFCDD2', '#FFFFFF', '#BBDEFB'], N=256)
    # Actually better: white to dark, where darker = more sites retained
    cmap = LinearSegmentedColormap.from_list('retention',
        ['#FBE9E7', '#FF8A65', '#BF360C'], N=256)

    im = ax.imshow(ct_pct.values, cmap=cmap, aspect='auto', vmin=0, vmax=100)

    # Text annotations: show both count and percentage
    for i in range(len(main_motifs)):
        for j, tp in enumerate(['T1', 'T2', 'T3']):
            count = ct.iloc[i, j]
            pct = ct_pct.iloc[i, j]
            color = 'white' if pct > 60 else 'black'
            ax.text(j, i, f'{count:,}\n({pct:.0f}%)',
                    ha='center', va='center', fontsize=7, color=color,
                    fontweight='bold' if j == 0 else 'normal')

    ax.set_xticks(range(3))
    ax.set_xticklabels(['T1 (24 h)', 'T2 (48 h)', 'T3 (72 h)'])
    ax.set_yticks(range(len(main_motifs)))
    ax.set_yticklabels(main_motifs, fontsize=9)
    ax.set_title('Motif-specific temporal patterns', fontsize=11,
                 fontweight='bold')

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label('% of T1 sites', fontsize=8)
    cbar.ax.tick_params(labelsize=7)


def _draw_gene_arrow(ax, start, end, y, h, strand, fc, ec, lw):
    """Draw a gene as an arrow-shaped polygon."""
    arrow_len = min(300, (end - start) * 0.15)  # arrow tip length
    if strand == 1:  # forward
        xs = [start, end - arrow_len, end, end - arrow_len, start]
        ys = [y - h/2, y - h/2, y, y + h/2, y + h/2]
    else:  # reverse
        xs = [start + arrow_len, end, end, start + arrow_len, start]
        ys = [y - h/2, y - h/2, y + h/2, y + h/2, y]
    ax.fill(xs, ys, facecolor=fc, edgecolor=ec, linewidth=lw, zorder=2)


def panel_c_browser_tracks(ax):
    """Genome browser-style tracks for SC_RS17645 region."""
    gbk = load_reference_gbk()
    df_4mc, df_6ma = load_methylation_census()

    region_start = 3396000
    region_end = 3409000
    region_span = region_end - region_start

    ax.set_xlim(region_start, region_end)
    ax.set_ylim(-0.5, 9.5)
    ax.set_title('SC_RS17645 region (3.396\u20133.409 Mb)', fontsize=10,
                 fontweight='bold')

    # ── Track 1: Gene annotations (y ~ 7-8.5) ──
    gene_y = 7.8
    gene_h = 0.8
    ax.axhline(gene_y, color='#EEEEEE', lw=0.3, zorder=0)

    genes_in_region = []
    for feat in gbk.features:
        if feat.type != 'CDS':
            continue
        s = int(feat.location.start)
        e = int(feat.location.end)
        if e < region_start or s > region_end:
            continue
        strand = feat.location.strand
        name = ''
        for key in ['locus_tag', 'gene']:
            if key in feat.qualifiers:
                name = feat.qualifiers[key][0]
                break
        genes_in_region.append((s, e, strand, name))

    # Draw genes
    for s, e, strand, name in genes_in_region:
        s_clip = max(s, region_start)
        e_clip = min(e, region_end)

        is_target = 'SC_RS17645' in name
        if is_target:
            fc, ec, lw = '#EDE7F6', COL_BOTH, 1.8
        elif strand == 1:
            fc, ec, lw = '#E3F2FD', '#64B5F6', 0.8
        else:
            fc, ec, lw = '#FFF3E0', '#FFB74D', 0.8

        _draw_gene_arrow(ax, s_clip, e_clip, gene_y, gene_h, strand,
                         fc, ec, lw)

    # Label only key genes (every other + SC_RS17645 always)
    labeled_x = []
    for s, e, strand, name in genes_in_region:
        s_clip = max(s, region_start)
        e_clip = min(e, region_end)
        mid = (s_clip + e_clip) / 2
        is_target = 'SC_RS17645' in name

        # Skip if too close to previous label
        too_close = any(abs(mid - lx) < region_span * 0.06 for lx in labeled_x)
        if too_close and not is_target:
            continue

        fontsize = 7 if is_target else 5.5
        fw = 'bold' if is_target else 'normal'
        color = COL_BOTH if is_target else '#616161'
        # Shorten label
        short = name.replace('SC_RS', '')
        if is_target:
            short = 'SC_RS17645'

        ax.text(mid, gene_y + gene_h/2 + 0.2, short,
                ha='center', va='bottom', fontsize=fontsize,
                fontweight=fw, color=color, fontstyle='italic')
        labeled_x.append(mid)

    # Highlight SC_RS17645 with background shading
    for s, e, strand, name in genes_in_region:
        if 'SC_RS17645' in name:
            ax.axvspan(s, e, alpha=0.08, color=COL_BOTH, zorder=0)

    # ── Track 2: 4mC lollipop plot (y ~ 4-6) ──
    base_4mc = 5.2
    ax.axhline(base_4mc, color='#FFCDD2', lw=0.5, zorder=0)
    ax.text(region_start + 100, base_4mc + 1.0, '4mC', fontsize=8,
            fontweight='bold', color=COL_4mC, va='bottom')

    meth_4mc = df_4mc[(df_4mc['position'] >= region_start) &
                       (df_4mc['position'] <= region_end)]
    tp_offsets = {'T1': 0, 'T2': 0.33, 'T3': 0.66}
    tp_alphas = {'T1': 1.0, 'T2': 0.6, 'T3': 0.35}
    for _, row in meth_4mc.iterrows():
        tp = row['timepoint']
        freq = row['frequency'] / 100
        x_off = tp_offsets[tp] * 120
        stem_h = freq * 1.2
        x = row['position'] + x_off
        # Stem
        ax.plot([x, x], [base_4mc, base_4mc + stem_h],
                color=COL_4mC, alpha=tp_alphas[tp], lw=1.0, zorder=3)
        # Head
        ax.plot(x, base_4mc + stem_h, 'o', color=COL_4mC,
                alpha=tp_alphas[tp], markersize=3, zorder=4)

    # ── Track 3: 6mA lollipop plot (y ~ 1-3) ──
    base_6ma = 2.0
    ax.axhline(base_6ma, color='#BBDEFB', lw=0.5, zorder=0)
    ax.text(region_start + 100, base_6ma + 1.0, '6mA', fontsize=8,
            fontweight='bold', color=COL_6mA, va='bottom')

    meth_6ma = df_6ma[(df_6ma['position'] >= region_start) &
                       (df_6ma['position'] <= region_end)]
    for _, row in meth_6ma.iterrows():
        tp = row['timepoint']
        freq = row['frequency'] / 100
        x_off = tp_offsets[tp] * 120
        stem_h = freq * 1.2
        x = row['position'] + x_off
        ax.plot([x, x], [base_6ma, base_6ma + stem_h],
                color=COL_6mA, alpha=tp_alphas[tp], lw=1.0, zorder=3)
        ax.plot(x, base_6ma + stem_h, 'o', color=COL_6mA,
                alpha=tp_alphas[tp], markersize=3, zorder=4)

    # ── Legend for timepoints ──
    for tp, alpha in [('T1', 1.0), ('T2', 0.6), ('T3', 0.35)]:
        ax.plot([], [], 'o-', color='#616161', alpha=alpha, markersize=4,
                lw=1.5, label=tp)
    ax.legend(fontsize=7, frameon=True, loc='lower right', ncol=3,
              fancybox=True, framealpha=0.8, edgecolor='#E0E0E0')

    # ── X axis ──
    ticks = np.arange(3396000, 3410000, 2000)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{t/1e6:.3f}' for t in ticks], fontsize=7)
    ax.set_xlabel('Genome position (Mb)')

    # Clean y axis
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)


def main():
    apply_style()
    print('=== Figure 5: Temporal Dynamics ===')

    fig = plt.figure(figsize=(mm_to_inch(174), mm_to_inch(180)))

    # Layout: top row 2 panels, bottom row full-width browser track
    gs = GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.45,
                  height_ratios=[1, 1])

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])  # full width

    print('Panel A: Temporal site counts...')
    panel_a_temporal_counts(ax_a)
    add_panel_label(ax_a, 'a')

    print('Panel B: Motif-specific temporal patterns...')
    panel_b_motif_temporal(ax_b)
    add_panel_label(ax_b, 'b')

    print('Panel C: Genome browser tracks...')
    panel_c_browser_tracks(ax_c)
    add_panel_label(ax_c, 'c', x=-0.05, y=1.05)

    # Save
    out_path = FIG_DIR / 'Figure5_temporal_dynamics'
    save_figure(fig, out_path)

    print('\n=== Figure 5 complete ===')


if __name__ == '__main__':
    main()
