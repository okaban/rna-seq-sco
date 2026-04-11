#!/usr/bin/env python3
"""
Figure 3: AAGCCCG Novel Motif
  Panel A: REBASE conservation (bar chart)
  Panel B: O/E ratio comparison (bar chart)
  Panel C: SC_RS17645 domain architecture + logic
  Panel D: SC_RS17645 expression vs AAGCCCG methylation timeline
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


def panel_a_rebase(ax):
    """REBASE conservation across Streptomyces species."""
    df = load_rebase_conservation()

    # Filter to our motifs only (from this study)
    our_motifs = df[df['source'] == 'This study'].copy()
    our_motifs = our_motifs.sort_values('rate_pct', ascending=True)

    motifs = our_motifs['motif'].values
    rates = our_motifs['rate_pct'].values
    ci_lo = our_motifs['ci_low'].values
    ci_hi = our_motifs['ci_high'].values
    n_matches = our_motifs['n_match'].values.astype(int)

    y = np.arange(len(motifs))
    colors = []
    for m in motifs:
        if m == 'AAGCCCG':
            colors.append(COL_BOTH)
        elif m == 'CCGG':
            colors.append(COL_4mC)
        else:
            colors.append(COL_6mA)

    xerr_lo = rates - ci_lo
    xerr_hi = ci_hi - rates
    xerr = np.array([xerr_lo, xerr_hi])

    ax.barh(y, rates, height=0.55, color=colors, edgecolor='white',
            linewidth=0.5, alpha=0.9, zorder=3)
    ax.errorbar(rates, y, xerr=xerr, fmt='none', ecolor='black',
                capsize=3, capthick=0.8, elinewidth=0.8, zorder=4)

    # Add n/82 labels to the right of error bars
    for i in range(len(motifs)):
        label_x = max(rates[i] + xerr_hi[i], rates[i]) + 1.5
        if motifs[i] == 'AAGCCCG':
            ax.text(label_x, y[i], f'{n_matches[i]}/82 (novel)',
                    fontsize=8, fontweight='bold', color=COL_BOTH, va='center')
        else:
            ax.text(label_x, y[i], f'{n_matches[i]}/82',
                    fontsize=8, va='center', color='#424242')

    ax.set_yticks(y)
    ax.set_yticklabels(motifs, fontfamily='monospace', fontweight='bold',
                       fontsize=10)
    ax.set_xlabel('Conservation in Streptomyces\nREBASE entries (%)')
    ax.set_xlim(-0.5, 42)
    ax.set_title('REBASE conservation', fontsize=11, fontweight='bold')
    ax.grid(axis='x', alpha=0.3, lw=0.5)


def panel_b_oe_ratio(ax):
    """O/E ratio comparison for M145 motifs."""
    df = load_motif_summary()

    oe_data = []
    for _, row in df.iterrows():
        motif = row['motif']
        oe_str = str(row['m145_oe'])
        try:
            oe_val = float(oe_str.split('(')[0].strip().split()[0])
        except (ValueError, IndexError):
            oe_val = float(oe_str)
        if '/' in motif:
            motif = 'CCGG family'
        oe_data.append({'motif': motif, 'oe': oe_val})

    df_oe = pd.DataFrame(oe_data)
    df_oe = df_oe.sort_values('oe', ascending=True)

    y = np.arange(len(df_oe))
    colors = []
    for m in df_oe['motif']:
        if m == 'AAGCCCG':
            colors.append(COL_BOTH)
        elif m == 'CCGG family':
            colors.append(COL_4mC)
        else:
            colors.append(COL_6mA)

    ax.barh(y, df_oe['oe'].values, height=0.55, color=colors,
            edgecolor='white', linewidth=0.5, alpha=0.9, zorder=3)

    # Reference line at O/E = 1.0
    ax.axvline(1.0, color='black', linestyle='--', linewidth=0.8,
               alpha=0.5, zorder=2)

    ax.set_yticks(y)
    ax.set_yticklabels(df_oe['motif'].values, fontfamily='monospace',
                       fontsize=8.5)
    ax.set_xlabel('Observed / Expected ratio')
    ax.set_title('Motif O/E in M145 genome', fontsize=11, fontweight='bold')
    ax.grid(axis='x', alpha=0.3, lw=0.5)

    # Annotate AAGCCCG
    aag_idx = list(df_oe['motif'].values).index('AAGCCCG')
    ax.text(df_oe.iloc[aag_idx]['oe'] + 0.05, y[aag_idx],
            'avoidance', fontsize=8, fontweight='bold',
            color=COL_BOTH, va='center')

    # Annotate the expected=1 line at the bottom
    ax.text(1.02, -0.7, 'O/E = 1', fontsize=7, color='#757575',
            fontstyle='italic')


def panel_c_domain(ax):
    """SC_RS17645 domain architecture and candidate logic."""
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis('off')

    # Title
    ax.text(50, 97, 'SC_RS17645 (SCO3104)', ha='center', fontsize=10,
            fontweight='bold')
    ax.text(50, 92, 'Candidate N-6 adenine methyltransferase', ha='center',
            fontsize=9, fontstyle='italic', color='#616161')

    # Logic flow — vertical chain of reasoning
    steps = [
        ('AAGCCCG = N6-methyladenine motif', COL_6mA),
        ('Requires N-6 adenine MTase', '#424242'),
        ('22 MTase genes in M145 genome', '#424242'),
        ('Only 1 annotated as "N-6 DNA methylase"', COL_BOTH),
    ]
    box_w = 72
    box_h = 7
    x_center = 42
    y_start = 84

    for i, (text, color) in enumerate(steps):
        y = y_start - i * 12
        # Box
        rect = mpatches.FancyBboxPatch(
            (x_center - box_w/2, y - box_h/2), box_w, box_h,
            boxstyle='round,pad=1.5', facecolor='white',
            edgecolor=color, linewidth=1.2)
        ax.add_patch(rect)
        ax.text(x_center, y, text, ha='center', va='center',
                fontsize=7.5, color=color)
        # Arrow to next
        if i < len(steps) - 1:
            ax.annotate('', xy=(x_center, y - box_h/2 - 1),
                        xytext=(x_center, y - box_h/2 - 4),
                        arrowprops=dict(arrowstyle='<-', lw=1.0,
                                        color='#BDBDBD'))

    # Result box
    result_y = y_start - len(steps) * 12 + 2
    rect = mpatches.FancyBboxPatch(
        (x_center - box_w/2, result_y - box_h/2), box_w, box_h,
        boxstyle='round,pad=1.5', facecolor='#EDE7F6',
        edgecolor=COL_BOTH, linewidth=1.5)
    ax.add_patch(rect)
    ax.text(x_center, result_y,
            'SC_RS17645: HsdM (Type I R-M) with TRD domain',
            ha='center', va='center', fontsize=7.5, fontweight='bold',
            color=COL_BOTH)

    # Domain architecture
    dom_y = 18
    protein_len = 679
    x_left = 8
    x_right = 92
    bar_w = x_right - x_left
    scale = bar_w / protein_len

    # Full protein
    rect = mpatches.FancyBboxPatch(
        (x_left, dom_y - 3), bar_w, 6,
        boxstyle='round,pad=0.5', facecolor='#ECEFF1',
        edgecolor='#90A4AE', linewidth=0.8)
    ax.add_patch(rect)

    # Domains
    domains = [
        (1, 484, 'IPR052916', '#C5CAE9', '#3F51B5', 5.5),
        (167, 393, 'N6_Mtase (PF02384)', '#D1C4E9', '#7E57C2', 7),
        (539, 669, 'TRD', '#BBDEFB', '#1565C0', 5.5),
    ]
    for start, end, name, fc, ec, h in domains:
        x = x_left + start * scale
        w = (end - start) * scale
        rect = mpatches.FancyBboxPatch(
            (x, dom_y - h/2), w, h,
            boxstyle='round,pad=0.3', facecolor=fc,
            edgecolor=ec, linewidth=1.5, alpha=0.95)
        ax.add_patch(rect)
        ax.text(x + w/2, dom_y, name, ha='center', va='center',
                fontsize=6.5, fontweight='bold', color=ec)

    # Position markers
    markers = [(1, 'left'), (167, 'center'), (393, 'center'),
               (484, 'center'), (539, 'center'), (669, 'center'),
               (679, 'right')]
    for pos, ha in markers:
        x = x_left + pos * scale
        ax.plot([x, x], [dom_y - 4.5, dom_y - 6], color='#9E9E9E',
                lw=0.5)
        ax.text(x, dom_y - 7.5, str(pos), ha=ha, va='top', fontsize=5.5,
                color='#757575')

    ax.text(50, dom_y - 12, f'{protein_len} aa  |  REBASE: M.ScoA3ORF3104P',
            ha='center', fontsize=7.5, color='#757575', fontstyle='italic')


def panel_d_expression_timeline(ax):
    """SC_RS17645 expression vs AAGCCCG methylation over time."""
    df_mtase = load_mtase_expression()
    sc17645 = df_mtase[df_mtase['locus_tag'] == 'SC_RS17645'].iloc[0]

    expr_vals = [0, sc17645['log2FC_T2vsT1'], sc17645['log2FC_T3vsT1']]

    # Load AAGCCCG methylation
    df_4mc, df_6ma = load_methylation_census()
    df_all = pd.concat([df_4mc, df_6ma], ignore_index=True)
    aag = df_all[df_all['final_motif'] == 'AAGCCCG'].copy()
    aag_freq = aag.groupby('timepoint')['frequency'].mean().reindex(['T1', 'T2', 'T3'])

    x = np.arange(3)

    # Primary axis: expression (green)
    color_expr = COL_GREEN
    ln1 = ax.plot(x, expr_vals, 'o-', color=color_expr, linewidth=2.5,
                  markersize=9, markerfacecolor=color_expr,
                  markeredgecolor='white', markeredgewidth=1.5,
                  zorder=4, label='SC_RS17645 expression')
    ax.axhline(0, color='gray', linestyle=':', linewidth=0.5, alpha=0.5)
    ax.set_ylabel('SC_RS17645\nlog$_2$FC vs T1', color=color_expr, fontsize=9)
    ax.tick_params(axis='y', labelcolor=color_expr)

    # Secondary axis: methylation (purple)
    ax2 = ax.twinx()
    color_meth = COL_BOTH
    ln2 = ax2.plot(x, aag_freq.values, 's--', color=color_meth, linewidth=2.0,
                   markersize=9, markerfacecolor=color_meth,
                   markeredgecolor='white', markeredgewidth=1.5,
                   zorder=3, label='AAGCCCG methylation')
    ax2.set_ylabel('AAGCCCG mean\nfrequency (%)', color=color_meth, fontsize=9)
    ax2.tick_params(axis='y', labelcolor=color_meth)

    ax.set_xticks(x)
    ax.set_xticklabels(['T1 (24 h)', 'T2 (48 h)', 'T3 (72 h)'])
    ax.set_title('MTase expression vs\nAAGCCCG methylation', fontsize=11,
                 fontweight='bold')

    # Combined legend
    lns = ln1 + ln2
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, fontsize=7.5, frameon=False, loc='lower left',
              bbox_to_anchor=(0.0, 0.0))

    # Key stat annotation
    padj = sc17645['padj_T2vsT1']
    ax.text(0.97, 0.95,
            f'T2 vs T1:\nlog$_2$FC = {expr_vals[1]:.2f}\np$_{{adj}}$ = {padj:.1e}',
            transform=ax.transAxes, fontsize=7, ha='right', va='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#E8F5E9',
                      edgecolor=color_expr, alpha=0.9, linewidth=0.8))


def main():
    apply_style()
    print('=== Figure 3: AAGCCCG Novel Motif ===')

    fig = plt.figure(figsize=(mm_to_inch(180), mm_to_inch(180)))

    # Layout: 2×2 grid, bottom row taller for domain panel
    gs = GridSpec(2, 2, figure=fig, hspace=0.5, wspace=0.5,
                  height_ratios=[0.8, 1.2])

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    print('Panel A: REBASE conservation...')
    panel_a_rebase(ax_a)
    add_panel_label(ax_a, 'a')

    print('Panel B: O/E ratio...')
    panel_b_oe_ratio(ax_b)
    add_panel_label(ax_b, 'b')

    print('Panel C: SC_RS17645 domain architecture...')
    panel_c_domain(ax_c)
    add_panel_label(ax_c, 'c', x=-0.05, y=1.02)

    print('Panel D: Expression vs methylation timeline...')
    panel_d_expression_timeline(ax_d)
    add_panel_label(ax_d, 'd')

    # Save
    out_path = FIG_DIR / 'Figure3_AAGCCCG_novel_motif'
    save_figure(fig, out_path)

    print('\n=== Figure 3 complete ===')


if __name__ == '__main__':
    main()
