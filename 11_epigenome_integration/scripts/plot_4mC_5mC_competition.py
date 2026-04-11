#!/usr/bin/env python3
"""
4mC/5mC Competition Hypothesis at CCGG Sites
Publication-quality figure showing evidence for dual cytosine modification.

Streptomyces coelicolor A3(2) M145 Project
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
from pathlib import Path

# === Paths ===
BASE = Path('/Users/okaban/bioinfo/rna-seq/11_epigenome_integration')
OUTPUT_DIR = BASE / 'analysis' / '22_4mC_5mC_competition'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# === Publication style ===
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 8,
    'axes.titlesize': 10,
    'axes.labelsize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 0.6,
})

# ====================================================================
# Panel A: CCGG 4mC dynamics + Dcm-like MTase expression
# ====================================================================
def draw_4mc_dynamics(ax):
    """4mC-CCGG site counts + Dcm-like MTase expression overlay."""
    timepoints = ['T1', 'T2', 'T3']
    x = np.arange(len(timepoints))

    # 4mC-CCGG site counts (from analysis)
    mc4_counts = [1576, 1961, 866]

    # Dcm-like MTase expression (log2FC relative to T1, 0 = T1 baseline)
    # SC_RS19770: T2vsT1=-0.94(ns), T3vsT1=+2.32
    # SC_RS36410: T2vsT1=+1.06(ns), T3vsT1=+5.34
    # Use mean of the two
    dcm_log2fc = [0, 0.06, 3.83]  # mean of (-0.94,1.06)=0.06, mean of (2.32,5.34)=3.83

    # Left axis: 4mC counts
    color1 = '#E53935'
    bars = ax.bar(x, mc4_counts, width=0.45, color=color1, alpha=0.7, edgecolor='black', lw=0.6)
    ax.set_ylabel('4mC-CCGG sites (count)', color=color1, fontsize=9)
    ax.tick_params(axis='y', labelcolor=color1)
    ax.set_ylim(0, 2500)

    # Annotate bars
    for i, (bar, val) in enumerate(zip(bars, mc4_counts)):
        ax.text(i, val + 50, str(val), ha='center', va='bottom', fontsize=8,
                fontweight='bold', color=color1)

    # Collapse annotation
    ax.annotate('56% collapse', xy=(2, 866), xytext=(2.3, 1500),
                fontsize=8, color=color1, fontweight='bold',
                arrowprops=dict(arrowstyle='->', color=color1, lw=1.0))

    # Right axis: Dcm-like expression
    ax2 = ax.twinx()
    color2 = '#1565C0'
    ax2.plot(x, dcm_log2fc, 'D-', color=color2, lw=2, markersize=8, zorder=5)
    ax2.set_ylabel('Dcm-like MTase\nlog$_2$FC (mean of 2)', color=color2, fontsize=8)
    ax2.tick_params(axis='y', labelcolor=color2)
    ax2.set_ylim(-1, 6)

    # Annotate Dcm expression
    for i, val in enumerate(dcm_log2fc):
        if i > 0:
            label = f'+{val:.1f}' if val > 0 else f'{val:.1f}'
            ax2.text(i + 0.12, val + 0.3, label, fontsize=7, color=color2, fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels(timepoints, fontsize=10, fontweight='bold')
    ax.set_title('A. 4mC-CCGG collapse despite Dcm-like MTase upregulation',
                 fontsize=9.5, fontweight='bold', loc='left')
    ax.spines['top'].set_visible(False)

    # Paradox box
    ax.text(0.02, 0.98, 'PARADOX: If Dcm-like produces 4mC,\n'
            'why does 4mC collapse when Dcm-like\n'
            'is upregulated 12-fold at T3?',
            transform=ax.transAxes, fontsize=6.5, va='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#FFEBEE', edgecolor='#E53935', lw=0.6))


# ====================================================================
# Panel B: GGCCGG enrichment in 4mC sites
# ====================================================================
def draw_ggccgg_enrichment(ax):
    """Show GGCCGG enrichment: genome vs methylated sites."""
    categories = ['Genome\n(all CCGGs)', '4mC methylated\nCCGGs (top 20)']
    ggccgg_pct = [24.4, 100.0]
    bare_ccgg_pct = [75.6, 0.0]

    x = np.arange(len(categories))
    width = 0.5

    bars1 = ax.bar(x, ggccgg_pct, width, label='Within GGCCGG', color='#E53935',
                   edgecolor='black', lw=0.6)
    bars2 = ax.bar(x, bare_ccgg_pct, width, bottom=ggccgg_pct, label='Bare CCGG only',
                   color='#BBDEFB', edgecolor='black', lw=0.6)

    # Annotations
    ax.text(0, 12, '24.4%\nGGCCGG', ha='center', fontsize=8, fontweight='bold', color='white')
    ax.text(0, 62, '75.6%\nbare CCGG', ha='center', fontsize=8, fontweight='bold', color='#1565C0')
    ax.text(1, 50, '100%\nGGCCGG', ha='center', fontsize=11, fontweight='bold', color='white')

    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=8)
    ax.set_ylabel('Fraction (%)', fontsize=9)
    ax.set_ylim(0, 115)
    ax.set_title('B. 4mC is enriched at GGCCGG = Pisciotta m5C motif',
                 fontsize=9.5, fontweight='bold', loc='left')
    ax.legend(fontsize=7, loc='upper right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Statistical note
    ax.text(0.5, -0.18, 'Pisciotta 2023 BS-seq: 5mC at GG$\\bf{C^{m5}}$CGG\n'
            'This study Nanopore: 4mC at GG$\\bf{C^{m4}}$CGG',
            transform=ax.transAxes, fontsize=7, ha='center', va='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#E3F2FD', edgecolor='#1565C0', lw=0.5))


# ====================================================================
# Panel C: Dual modification model
# ====================================================================
def draw_model(ax):
    """Draw the m4C/m5C dual modification model."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('C. Proposed dual modification model at CCGG',
                 fontsize=9.5, fontweight='bold', loc='left')

    # DNA sequence
    dna_y = 7.5
    ax.text(5, dna_y + 0.7, 'Same cytosine in CCGG:', ha='center', fontsize=8, fontweight='bold')
    ax.text(5, dna_y, '5\'-...GG-C-CGG-...-3\'', ha='center', fontsize=10,
            fontfamily='monospace', fontweight='bold')
    ax.annotate('', xy=(5.15, dna_y - 0.3), xytext=(5.15, dna_y - 1.0),
                arrowprops=dict(arrowstyle='->', color='black', lw=1.5))

    # T1-T2 box (left)
    t1_x, t1_y = 2.5, 4.5
    rect1 = FancyBboxPatch((t1_x - 1.8, t1_y - 1.2), 3.6, 2.4,
                           boxstyle="round,pad=0.2", fc='#FFCDD2', ec='#E53935', lw=1.0)
    ax.add_patch(rect1)
    ax.text(t1_x, t1_y + 0.7, 'T1-T2: m4C dominant', ha='center', fontsize=8,
            fontweight='bold', color='#C62828')
    ax.text(t1_x, t1_y + 0.1, 'Nanopore detection', ha='center', fontsize=7, color='#C62828')
    ax.text(t1_x, t1_y - 0.5, '~2,000 sites', ha='center', fontsize=7, color='#C62828')
    ax.text(t1_x, t1_y - 0.9, 'Unknown MTase', ha='center', fontsize=6.5, fontstyle='italic',
            color='#C62828')

    # T3 box (right)
    t3_x, t3_y = 7.5, 4.5
    rect3 = FancyBboxPatch((t3_x - 1.8, t3_y - 1.2), 3.6, 2.4,
                           boxstyle="round,pad=0.2", fc='#BBDEFB', ec='#1565C0', lw=1.0)
    ax.add_patch(rect3)
    ax.text(t3_x, t3_y + 0.7, 'T3: m5C dominant', ha='center', fontsize=8,
            fontweight='bold', color='#0D47A1')
    ax.text(t3_x, t3_y + 0.1, 'BS-seq detection', ha='center', fontsize=7, color='#0D47A1')
    ax.text(t3_x, t3_y - 0.5, '~866 4mC sites remain', ha='center', fontsize=7, color='#0D47A1')
    ax.text(t3_x, t3_y - 0.9, 'Dcm-like MTases\n(12x upregulated)', ha='center',
            fontsize=6.5, fontstyle='italic', color='#0D47A1', linespacing=1.0)

    # Arrow between boxes
    ax.annotate('', xy=(t3_x - 2.0, t1_y), xytext=(t1_x + 2.0, t1_y),
                arrowprops=dict(arrowstyle='->', color='black', lw=2.0,
                                connectionstyle='arc3,rad=0'))
    ax.text(5, t1_y + 0.4, 'Phase switch', ha='center', fontsize=7, fontweight='bold')

    # Evidence box at bottom
    evidence_y = 1.5
    rect_ev = FancyBboxPatch((0.5, evidence_y - 1.0), 9.0, 2.0,
                             boxstyle="round,pad=0.2", fc='#F5F5F5', ec='#666666', lw=0.6)
    ax.add_patch(rect_ev)
    ax.text(5, evidence_y + 0.6, 'Key evidence for competition:', ha='center', fontsize=7.5,
            fontweight='bold')
    ax.text(5, evidence_y + 0.0, '1. Same target: 4mC and 5mC both at inner C of GGCCGG',
            ha='center', fontsize=6.5)
    ax.text(5, evidence_y - 0.4, '2. Temporal anticorrelation: Dcm-like up + 4mC down at T3',
            ha='center', fontsize=6.5)
    ax.text(5, evidence_y - 0.8, '3. BS-seq blindness: Bisulfite cannot distinguish unmethylated C from m4C',
            ha='center', fontsize=6.5)


# ====================================================================
# Panel D: Evidence summary table
# ====================================================================
def draw_evidence_table(ax):
    """Draw evidence summary as a formatted table."""
    ax.axis('off')
    ax.set_title('D. Lines of evidence for m4C/m5C competition',
                 fontsize=9.5, fontweight='bold', loc='left')

    table_data = [
        ['#', 'Evidence', 'Source', 'Strength'],
        ['1', '4mC at CCGG (Nanopore)', 'This study', 'Direct'],
        ['2', '5mC at GGCCGG (BS-seq)', 'Pisciotta 2023', 'Direct'],
        ['3', '4mC consensus = GGCCGG\n(100% of top 20 motifs)', 'This study', 'Strong'],
        ['4', 'Dcm-like MTase +12x at T3\nyet 4mC collapses 56%', 'This study', 'Strong'],
        ['5', 'BS-seq cannot detect m4C\n(treated as unmethylated)', 'Methodology', 'Logical'],
        ['6', 'Nanopore has low m5C\nsensitivity (<0.01%)', 'Methodology', 'Logical'],
    ]

    # Draw table
    cell_h = 0.14
    cell_w = [0.04, 0.40, 0.28, 0.12]
    y_start = 0.88

    colors_row = ['#E8EAF6', '#FFFFFF', '#F5F5F5', '#FFFFFF', '#F5F5F5', '#FFFFFF', '#F5F5F5']

    for i, row in enumerate(table_data):
        y = y_start - i * cell_h
        x = 0.02
        for j, (cell, cw) in enumerate(zip(row, cell_w)):
            weight = 'bold' if i == 0 else 'normal'
            fontsize = 6.5 if i > 0 else 7
            color = '#333333' if i == 0 else 'black'
            bgcolor = '#C5CAE9' if i == 0 else colors_row[i]

            rect = FancyBboxPatch((x, y - cell_h * 0.4), cw - 0.01, cell_h * 0.9,
                                  boxstyle="square,pad=0", fc=bgcolor, ec='#BDBDBD', lw=0.3,
                                  transform=ax.transAxes)
            ax.add_patch(rect)
            ax.text(x + cw/2, y + cell_h * 0.05, cell, transform=ax.transAxes,
                    ha='center', va='center', fontsize=fontsize, fontweight=weight,
                    color=color, linespacing=0.9)
            x += cw


# ====================================================================
# Compose Figure
# ====================================================================
def main():
    fig = plt.figure(figsize=(180/25.4, 240/25.4))

    gs = fig.add_gridspec(2, 2, hspace=0.40, wspace=0.35,
                          height_ratios=[1.0, 1.2])

    ax_a = fig.add_subplot(gs[0, 0])    # 4mC dynamics
    ax_b = fig.add_subplot(gs[0, 1])    # GGCCGG enrichment
    ax_c = fig.add_subplot(gs[1, 0])    # Model
    ax_d = fig.add_subplot(gs[1, 1])    # Evidence table

    draw_4mc_dynamics(ax_a)
    draw_ggccgg_enrichment(ax_b)
    draw_model(ax_c)
    draw_evidence_table(ax_d)

    # Save
    out_base = OUTPUT_DIR / 'fig_4mC_5mC_competition'
    fig.savefig(f'{out_base}.pdf')
    fig.savefig(f'{out_base}.svg')
    fig.savefig(f'{out_base}.png', dpi=300)
    plt.close()
    print(f'Saved: {out_base}.pdf/.svg/.png')


if __name__ == '__main__':
    main()
