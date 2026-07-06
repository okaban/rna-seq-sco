#!/usr/bin/env python3
"""
FigS1 (revised): Motif Discovery Landscape — All Candidates and Attribution

Panel A: 4mC MEME results — MEME-1 logo (reference to main Fig 1d) + all-motif table
Panel B: 6mA MEME results — MEME-1 logo + MEME-2 logo (secondary candidate) + all-motif table
Panel C: Final motif attribution proportions (4mC and 6mA donut charts)

Replaces the previous FigS1 (4mC reclassification) with a motif-focused overview
that transparently shows all discovered candidate motifs and their quantitative
contributions to the methylome.
"""

import sys
import importlib
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
import logomaker

# ── Paths ─────────────────────────────────────────────────────────────────────
MEME_4mC = (EPIGENOME / 'archive' / 'v1_weighted_minreps2' /
            '07_motif_analysis' / 'meme_4mC' / 'meme.txt')
MEME_6mA = (EPIGENOME / 'archive' / 'v1_weighted_minreps2' /
            '07_motif_analysis' / 'meme_6mA' / 'meme.txt')

# Nucleotide colors matching main figure style
NUC_COLORS = {'A': '#2ecc71', 'C': '#3498db', 'G': '#e67e22', 'T': '#e74c3c'}


# ── Helper functions ───────────────────────────────────────────────────────────

def sig_label(evalue_str):
    """Convert E-value string to significance symbol."""
    try:
        e = float(evalue_str)
        if e < 1e-100:
            return '***'
        if e < 1e-10:
            return '**'
        if e < 0.05:
            return '*'
        return 'n.s.'
    except Exception:
        return '-'


def draw_logo(ax, pwm, title='', subtitle=''):
    """Draw a sequence logo from a probability matrix."""
    info = logomaker.transform_matrix(
        pwm, from_type='probability', to_type='information')
    logomaker.Logo(info, ax=ax, color_scheme=NUC_COLORS, show_spines=False)
    ax.set_xticks([])
    ax.set_yticks([0, 1, 2])
    ax.set_yticklabels(['0', '1', '2'], fontsize=6)
    ax.set_ylabel('bits', fontsize=6)
    ax.spines['left'].set_visible(True)
    ax.spines['bottom'].set_visible(False)
    if title:
        ax.set_title(title, fontsize=8, fontweight='bold', linespacing=1.3)
    if subtitle:
        ax.text(0.5, -0.28, subtitle, transform=ax.transAxes, fontsize=7,
                ha='center', va='top', color='#555555', style='italic')


def draw_meme_table(ax, motifs, highlight_ranks=None):
    """Draw a styled summary table of all MEME motifs.

    Rows corresponding to highlight_ranks are shaded blue.
    """
    highlight_ranks = set(highlight_ranks or [])
    ax.axis('off')

    rows = []
    for rank_key in sorted(motifs.keys(), key=lambda x: int(x.split('-')[1])):
        m = motifs[rank_key]
        rank_num = int(rank_key.split('-')[1])
        sig = sig_label(m['evalue'])
        rows.append((rank_num, rank_key, m['name'], f"{m['nsites']:,}",
                     m['evalue'], sig))

    col_labels = ['Rank', 'Consensus', 'Sites', 'E-value', 'Sig.']
    col_x = [0.00, 0.14, 0.56, 0.72, 0.90]
    col_w = [0.14, 0.42, 0.16, 0.18, 0.10]
    row_h = 0.85
    n_rows = len(rows)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, n_rows + 1.1)

    # Header row
    for label, cx, cw in zip(col_labels, col_x, col_w):
        ax.add_patch(plt.Rectangle((cx, n_rows + 0.1), cw, row_h,
                                   facecolor=COL_DARK, edgecolor='none'))
        ax.text(cx + cw / 2, n_rows + 0.1 + row_h / 2, label,
                ha='center', va='center', fontsize=7,
                fontweight='bold', color='white')

    # Data rows (bottom-to-top indexing so rank 1 is at top)
    for i, (rank_num, rank_key, consensus, nsites, evalue, sig) in enumerate(rows):
        y = n_rows - i - 1 + 0.1
        highlight = rank_num in highlight_ranks
        bg = '#E3F2FD' if highlight else ('#F5F5F5' if i % 2 == 0 else 'white')

        vals = [rank_key, consensus, nsites, evalue, sig]
        for j, (val, cx, cw) in enumerate(zip(vals, col_x, col_w)):
            ax.add_patch(plt.Rectangle((cx, y), cw, row_h,
                                       facecolor=bg, edgecolor='#DDDDDD',
                                       linewidth=0.4))
            if j == 1:  # consensus — monospace, left-aligned
                ax.text(cx + 0.02, y + row_h / 2, val,
                        ha='left', va='center', fontsize=6.5,
                        fontfamily='monospace',
                        fontweight='bold' if highlight else 'normal',
                        color=COL_DARK)
            elif j == 4:  # significance — coloured
                sig_color = {
                    '***': '#C62828', '**': '#E64A19',
                    '*': '#F57C00', 'n.s.': '#9E9E9E'
                }.get(val, COL_DARK)
                ax.text(cx + cw / 2, y + row_h / 2, val,
                        ha='center', va='center', fontsize=7,
                        color=sig_color, fontweight='bold')
            else:
                ax.text(cx + cw / 2, y + row_h / 2, val,
                        ha='center', va='center', fontsize=6.5,
                        color=COL_DARK,
                        fontweight='bold' if (highlight and j == 0) else 'normal')

    # Legend note for highlight colour
    ax.text(0.0, -0.08, '*** E < 1e-100   ** E < 1e-10   * E < 0.05   n.s. not significant',
            transform=ax.transAxes, fontsize=6, color='#666666', va='top')
    if highlight_ranks:
        ax.add_patch(mpatches.FancyBboxPatch(
            (0.0, -0.22), 0.22, 0.10, transform=ax.transAxes,
            boxstyle='round,pad=0.01', facecolor='#E3F2FD', edgecolor='#90CAF9'))
        ax.text(0.23, -0.17, '= shown as logo', transform=ax.transAxes,
                fontsize=6, color='#1565C0', va='center')


def draw_attribution_donut(ax, census_df, motif_groups, color_map, title):
    """Draw a donut chart of final motif attribution proportions.

    motif_groups: OrderedDict {display_label: [list of final_motif values]}
    """
    df_uniq = census_df.drop_duplicates(subset=['chrom', 'position', 'strand'])
    total = len(df_uniq)

    counts = {}
    for label, motif_vals in motif_groups.items():
        n = df_uniq['final_motif'].isin(motif_vals).sum()
        counts[label] = n

    labels = list(counts.keys())
    sizes = [counts[k] for k in labels]
    colors = [color_map.get(k, '#BDBDBD') for k in labels]

    wedges, _, autotexts = ax.pie(
        sizes,
        labels=None,
        colors=colors,
        startangle=90,
        counterclock=False,
        autopct=lambda p: f'{p:.1f}%' if p > 3 else '',
        pctdistance=0.70,
        wedgeprops={'edgecolor': 'white', 'linewidth': 1.2, 'width': 0.58},
    )
    for at in autotexts:
        at.set_fontsize(6.5)
        at.set_color('white')
        at.set_fontweight('bold')

    # Legend with n and %
    legend_labels = [
        f'{k}\n{counts[k]:,} sites ({counts[k] / total * 100:.1f}%)'
        for k in labels
    ]
    ax.legend(wedges, legend_labels,
              loc='center left', bbox_to_anchor=(0.92, 0.5),
              fontsize=7, frameon=False,
              handlelength=1.0, handletextpad=0.5, labelspacing=1.1)
    ax.set_title(title, fontsize=9, fontweight='bold', pad=6)


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    apply_style()
    print('=== FigS1 (revised): Motif Discovery Landscape ===')

    # Load MEME motifs
    print('Loading MEME motifs...')
    motifs_4mc = parse_meme_pwm(MEME_4mC)
    motifs_6ma = parse_meme_pwm(MEME_6mA)

    # Load census
    print('Loading methylation census...')
    df_4mc, df_6ma = load_methylation_census()

    # ── Figure layout ────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(mm_to_inch(174), mm_to_inch(235)))

    outer = gridspec.GridSpec(
        3, 1, figure=fig,
        hspace=0.55,
        left=0.07, right=0.97, top=0.97, bottom=0.04,
        height_ratios=[1.0, 1.2, 1.1],
    )

    # Row 1 — 4mC: logo + table
    gs_a = gridspec.GridSpecFromSubplotSpec(
        1, 2, subplot_spec=outer[0], wspace=0.40, width_ratios=[1, 1.9])
    ax_logo_4mc = fig.add_subplot(gs_a[0])
    ax_tbl_4mc = fig.add_subplot(gs_a[1])

    # Row 2 — 6mA: logo1 + logo2 + table
    gs_b = gridspec.GridSpecFromSubplotSpec(
        1, 3, subplot_spec=outer[1], wspace=0.40, width_ratios=[1, 1, 1.6])
    ax_logo_6ma1 = fig.add_subplot(gs_b[0])
    ax_logo_6ma2 = fig.add_subplot(gs_b[1])
    ax_tbl_6ma = fig.add_subplot(gs_b[2])

    # Row 3 — Attribution donuts
    gs_c = gridspec.GridSpecFromSubplotSpec(
        1, 2, subplot_spec=outer[2], wspace=0.70)
    ax_pie4 = fig.add_subplot(gs_c[0])
    ax_pie6 = fig.add_subplot(gs_c[1])

    # ── Panel A: 4mC ─────────────────────────────────────────────────────────
    print('Panel A: 4mC MEME motifs...')
    add_panel_label(ax_logo_4mc, 'a', x=-0.10, y=1.15)

    m1_4mc = motifs_4mc['MEME-1']
    draw_logo(
        ax_logo_4mc, m1_4mc['pwm'],
        title=f"4mC MEME-1  (= Main Fig. 1d)\n{m1_4mc['name']}",
        subtitle=f"n = {m1_4mc['nsites']:,} sites   E = {m1_4mc['evalue']}",
    )

    draw_meme_table(ax_tbl_4mc, motifs_4mc, highlight_ranks=[1])
    ax_tbl_4mc.set_title('All 4mC MEME motifs (nmotifs = 5)', fontsize=8,
                          fontweight='bold', pad=4)

    # ── Panel B: 6mA ─────────────────────────────────────────────────────────
    print('Panel B: 6mA MEME motifs...')
    add_panel_label(ax_logo_6ma1, 'b', x=-0.10, y=1.15)

    m1_6ma = motifs_6ma['MEME-1']
    draw_logo(
        ax_logo_6ma1, m1_6ma['pwm'],
        title=f"6mA MEME-1  (= Main Fig. 1d)\n{m1_6ma['name']}",
        subtitle=f"n = {m1_6ma['nsites']:,} sites   E = {m1_6ma['evalue']}",
    )

    m2_6ma = motifs_6ma['MEME-2']
    draw_logo(
        ax_logo_6ma2, m2_6ma['pwm'],
        title=f"6mA MEME-2  (secondary)\n{m2_6ma['name']}",
        subtitle=f"n = {m2_6ma['nsites']:,} sites   E = {m2_6ma['evalue']}",
    )

    draw_meme_table(ax_tbl_6ma, motifs_6ma, highlight_ranks=[1, 2])
    ax_tbl_6ma.set_title('All 6mA MEME motifs (nmotifs = 5)', fontsize=8,
                          fontweight='bold', pad=4)

    # ── Panel C: Attribution donuts ───────────────────────────────────────────
    # Unified presentation: AAGCCCG receives COL_BOTH (purple) in BOTH charts
    # to highlight dual-modification across modification types.
    print('Panel C: Motif attribution...')
    add_panel_label(ax_pie4, 'c', x=-0.08, y=1.15)

    # 4mC: TGGCCGGC / GGCCGG / CCGG are all sub-contexts of the GCCGGC palindrome;
    # AAGCCCG is the dual-modification motif shared with 6mA.
    motif_groups_4mc = {
        'GCCGGC': ['TGGCCGGC', 'GGCCGG', 'CCGG (other)'],
        'AAGCCCG': ['AAGCCCG'],
    }
    color_4mc = {
        'GCCGGC': COL_4mC,
        'AAGCCCG': COL_BOTH,
    }
    draw_attribution_donut(ax_pie4, df_4mc, motif_groups_4mc, color_4mc,
                           '4mC motif attribution')
    # Annotate 100% assignment
    ax_pie4.text(0.5, -0.06, '100% of sites assigned',
                 transform=ax_pie4.transAxes, fontsize=7,
                 ha='center', color='#388E3C', style='italic')

    # 6mA: AAGCCCG context uses same COL_BOTH (purple) to link dual-modification theme.
    # CCGKCA-related: MEME-2/STREME secondary cluster.
    # "Other assigned" (GAACCGG/CGGCAACC/CTGCTCGCCG; 3.8%) merged into Unassigned
    # because the distinction between minor STREME candidates and truly unassigned
    # is not meaningful at this resolution.
    motif_groups_6ma = {
        'AAGCCCG': ['AAGCCCG'],
        'CCGKCA': ['CCGKCA', 'GCCG', 'CCGC', 'CCGG', 'CCSGG'],
        'GATC': ['GATC'],
        'Unassigned': ['GAACCGG', 'CGGCAACC', 'CTGCTCGCCG', 'unassigned'],
    }
    color_6ma = {
        'AAGCCCG': COL_BOTH,
        'CCGKCA': '#42A5F5',
        'GATC': '#0D47A1',
        'Unassigned': '#CFD8DC',
    }
    draw_attribution_donut(ax_pie6, df_6ma, motif_groups_6ma, color_6ma,
                           '6mA motif attribution')

    # ── Save ─────────────────────────────────────────────────────────────────
    out_path = FIG_SUP_DIR / 'FigS1_motif_landscape'
    save_figure(fig, out_path)
    print(f'\nSaved: {out_path}.pdf / .svg')
    print('=== Done ===')


if __name__ == '__main__':
    main()
