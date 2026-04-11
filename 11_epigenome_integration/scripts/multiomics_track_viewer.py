#!/usr/bin/env python3
"""
Multi-omics Genomic Track Viewer
Publication-quality genome browser style visualizations:
  A: 6mA methylation frequency by timepoint
  B: 4mC methylation frequency by timepoint
  C: Gene expression changes (log2FC)
  D: Gene annotation with strand arrows

Streptomyces coelicolor A3(2) M145 Project
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from matplotlib.patches import FancyArrowPatch
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# === Paths ===
BASE = Path('/Users/okaban/bioinfo/rna-seq/11_epigenome_integration')
INTEGRATION_DIR = BASE / 'analysis' / '01_integration'
OUTPUT_DIR = BASE / 'analysis' / '04_multiomics_tracks'
GFF_PATH = Path('/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/'
                'M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/'
                'GCF_000203835.1/genomic.gff')
METHYL_PATH = INTEGRATION_DIR / 'high_confidence_sites_weighted.csv'
EXPRESSION_PATH = INTEGRATION_DIR / 'integrated_methyl_expression_weighted.csv'

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# === Publication style ===
# Colorblind-friendly palette (consistent with other project figures)
COL_T1 = '#4DAF4A'   # Green
COL_T2 = '#E66100'   # Orange
COL_T3 = '#5D3A9B'   # Purple
COL_T3T2 = '#808080'  # Gray
COL_TARGET = '#D32F2F'  # Red for target gene
COL_OTHER = '#607D8B'   # Blue-gray for other genes
COL_PROMOTER = '#FFF9C4'  # Light yellow for promoter

plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 10,
    'axes.titlesize': 11,
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'axes.spines.top': False,
    'axes.spines.right': False,
})

# Figure size: 180mm width (journal standard), ~200mm height
FIG_W_MM = 180
FIG_H_MM = 195
FIG_W = FIG_W_MM / 25.4
FIG_H = FIG_H_MM / 25.4

# === Target genes ===
TARGETS = {
    'RamR': {
        'new_locus': 'SC_RS35610',
        'old_locus': 'SCO6685',
        'start': 7426396,
        'end': 7427004,
        'strand': '-',
        'function': 'Response regulator (SapB activator)',
        'region_start': 7425000,
        'region_end': 7428500,
    },
    'NsdB': {
        'new_locus': 'SC_RS38475',
        'old_locus': 'SCO7252',
        'start': 8062135,
        'end': 8063643,
        'strand': '+',
        'function': 'DNA-binding protein (master regulator)',
        'region_start': 8060500,
        'region_end': 8065000,
    },
    'Red_SCO5897': {
        'new_locus': 'SC_RS31730',
        'old_locus': 'SCO5897',
        'start': 6462296,
        'end': 6463483,
        'strand': '+',
        'function': 'Red cluster oxygenase',
        'region_start': 6460500,
        'region_end': 6465000,
    },
    'Act_SCO5079': {
        'new_locus': 'SC_RS27555',
        'old_locus': 'SCO5079',
        'start': 5520857,
        'end': 5521741,
        'strand': '+',
        'function': 'Act cluster NmrA protein',
        'region_start': 5519000,
        'region_end': 5523500,
    },
    'Cpk_SCO6284': {
        'new_locus': 'SC_RS33670',
        'old_locus': 'SCO6284',
        'start': 6943673,
        'end': 6945265,
        'strand': '+',
        'function': 'Cpk cluster carboxylase',
        'region_start': 6942000,
        'region_end': 6947000,
    },
}


def load_methylation_data():
    return pd.read_csv(METHYL_PATH)


def load_expression_data():
    return pd.read_csv(EXPRESSION_PATH)


def parse_gff_region(gff_path, chrom, start, end):
    genes = []
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            if parts[0] != chrom or parts[2] != 'gene':
                continue

            gene_start = int(parts[3])
            gene_end = int(parts[4])
            if gene_end < start or gene_start > end:
                continue

            strand = parts[6]
            attrs = {}
            for attr in parts[8].split(';'):
                if '=' in attr:
                    key, val = attr.split('=', 1)
                    attrs[key] = val

            genes.append({
                'start': gene_start,
                'end': gene_end,
                'strand': strand,
                'locus_tag': attrs.get('locus_tag', ''),
                'old_locus_tag': attrs.get('old_locus_tag', ''),
                'name': attrs.get('Name', ''),
                'gene': attrs.get('gene', ''),
            })
    return genes


def get_methylation_in_region(methyl_df, start, end):
    return methyl_df[
        (methyl_df['position'] >= start) & (methyl_df['position'] <= end)
    ].copy()


def format_kb_position(x, pos):
    """Format genomic position as X,XXX kb."""
    return f'{x / 1000:,.1f}'


def assign_gene_levels(genes, region_start, region_end):
    """Assign y-levels to genes to avoid overlapping arrows.

    Uses a greedy interval scheduling approach: for each gene (sorted by start),
    place it on the lowest level where it doesn't overlap with any existing gene
    (with a small gap buffer).
    """
    gap = (region_end - region_start) * 0.02  # 2% buffer
    sorted_genes = sorted(genes, key=lambda g: g['start'])
    levels = []  # list of lists, each sub-list tracks occupied intervals

    for gene in sorted_genes:
        placed = False
        for level_idx, occupied in enumerate(levels):
            # Check if gene fits on this level
            fits = True
            for occ_start, occ_end in occupied:
                if gene['start'] - gap < occ_end and gene['end'] + gap > occ_start:
                    fits = False
                    break
            if fits:
                occupied.append((gene['start'], gene['end']))
                gene['level'] = level_idx
                placed = True
                break
        if not placed:
            levels.append([(gene['start'], gene['end'])])
            gene['level'] = len(levels) - 1

    return genes, len(levels)


def draw_gene_arrow(ax, gene, y_center, height, is_target, region_span):
    """Draw a clean gene arrow with proper proportions."""
    gene_start = gene['start']
    gene_end = gene['end']
    gene_len = gene_end - gene_start
    head_len = min(region_span * 0.015, gene_len * 0.2)

    fc = COL_TARGET if is_target else COL_OTHER
    alpha = 1.0 if is_target else 0.7

    if gene['strand'] == '+':
        body_end = gene_end - head_len
        # Body rectangle
        body = mpatches.FancyBboxPatch(
            (gene_start, y_center - height / 2), body_end - gene_start, height,
            boxstyle='round,pad=0', fc=fc, ec='none', alpha=alpha, zorder=3,
        )
        ax.add_patch(body)
        # Arrowhead
        arrow_xs = [body_end, gene_end, body_end]
        arrow_ys = [y_center - height * 0.7, y_center, y_center + height * 0.7]
        ax.fill(arrow_xs, arrow_ys, fc=fc, ec='none', alpha=alpha, zorder=3)
    else:
        body_start = gene_start + head_len
        body = mpatches.FancyBboxPatch(
            (body_start, y_center - height / 2), gene_end - body_start, height,
            boxstyle='round,pad=0', fc=fc, ec='none', alpha=alpha, zorder=3,
        )
        ax.add_patch(body)
        arrow_xs = [body_start, gene_start, body_start]
        arrow_ys = [y_center - height * 0.7, y_center, y_center + height * 0.7]
        ax.fill(arrow_xs, arrow_ys, fc=fc, ec='none', alpha=alpha, zorder=3)


def create_multiomics_track(target_name, target_info, methyl_df, expr_df):
    """Create publication-quality multi-omics track visualization."""

    region_start = target_info['region_start']
    region_end = target_info['region_end']
    region_span = region_end - region_start
    gene_start = target_info['start']
    gene_end = target_info['end']

    # Load data
    methyl_region = get_methylation_in_region(methyl_df, region_start, region_end)
    expr_data = expr_df[expr_df['gene_id'] == target_info['new_locus']]
    genes = parse_gff_region(GFF_PATH, 'NC_003888.3', region_start, region_end)

    # Assign gene levels to avoid overlap
    genes, n_levels = assign_gene_levels(genes, region_start, region_end)

    # Adjust gene panel height ratio based on number of levels
    gene_ratio = max(0.5, 0.3 * n_levels)

    fig, axes = plt.subplots(
        4, 1, figsize=(FIG_W, FIG_H),
        height_ratios=[1, 1, 1, gene_ratio],
        gridspec_kw={'hspace': 0.35},
    )

    # === Promoter highlight (all tracks except gene) ===
    if target_info['strand'] == '+':
        prom_start = gene_start - 300
        prom_end = gene_start + 50
    else:
        prom_start = gene_end - 50
        prom_end = gene_end + 300

    for ax in axes[:3]:
        ax.axvspan(prom_start, prom_end, alpha=0.25, color=COL_PROMOTER,
                   zorder=0, ec='none')

    # === Panel labels ===
    panel_labels = ['A', 'B', 'C', 'D']
    for i, (ax, label) in enumerate(zip(axes, panel_labels)):
        ax.text(-0.08, 1.05, label, transform=ax.transAxes,
                fontsize=13, fontweight='bold', va='top', ha='left')

    # Timepoint plotting config
    tp_config = [
        ('T1', COL_T1, 'o'),
        ('T2', COL_T2, 's'),
        ('T3', COL_T3, 'D'),
    ]

    # === Track A: 6mA methylation ===
    ax1 = axes[0]
    for tp, color, marker in tp_config:
        tp_data = methyl_region[
            (methyl_region['timepoint'] == tp) &
            (methyl_region['mod_type'] == '6mA')
        ]
        if len(tp_data) > 0:
            ax1.scatter(tp_data['position'], tp_data['weighted_mod_freq'],
                        c=color, alpha=0.8, s=45, label=tp, marker=marker,
                        edgecolors='white', linewidth=0.4, zorder=5)

    ax1.axhline(50, color='#BDBDBD', linestyle=':', linewidth=0.7, zorder=1)
    ax1.axvline(gene_start, color='#9E9E9E', linestyle='--', linewidth=0.7, alpha=0.6)
    ax1.axvline(gene_end, color='#9E9E9E', linestyle='--', linewidth=0.7, alpha=0.6)
    ax1.set_ylabel('6mA frequency (%)')
    ax1.set_xlim(region_start, region_end)
    ax1.set_ylim(0, 105)
    ax1.set_title(f'6mA methylation', fontsize=11, fontweight='bold', loc='left')
    ax1.tick_params(labelbottom=False)
    ax1.legend(loc='upper right', frameon=True, fancybox=False,
               edgecolor='#BDBDBD', framealpha=0.9)

    # === Track B: 4mC methylation ===
    ax2 = axes[1]
    for tp, color, marker in tp_config:
        tp_data = methyl_region[
            (methyl_region['timepoint'] == tp) &
            (methyl_region['mod_type'] == '4mC')
        ]
        if len(tp_data) > 0:
            ax2.scatter(tp_data['position'], tp_data['weighted_mod_freq'],
                        c=color, alpha=0.8, s=45, label=tp, marker=marker,
                        edgecolors='white', linewidth=0.4, zorder=5)

    ax2.axhline(50, color='#BDBDBD', linestyle=':', linewidth=0.7, zorder=1)
    ax2.axvline(gene_start, color='#9E9E9E', linestyle='--', linewidth=0.7, alpha=0.6)
    ax2.axvline(gene_end, color='#9E9E9E', linestyle='--', linewidth=0.7, alpha=0.6)
    ax2.set_ylabel('4mC frequency (%)')
    ax2.set_xlim(region_start, region_end)
    ax2.set_ylim(0, 105)
    ax2.set_title('4mC methylation', fontsize=11, fontweight='bold', loc='left')
    ax2.tick_params(labelbottom=False)
    ax2.legend(loc='upper right', frameon=True, fancybox=False,
               edgecolor='#BDBDBD', framealpha=0.9)

    # === Track C: Expression bar chart ===
    ax3 = axes[2]
    if len(expr_data) > 0:
        row = expr_data.iloc[0]
        log2fc_t2t1 = row.get('log2FC_T2_vs_T1', 0)
        log2fc_t3t1 = row.get('log2FC_T3_vs_T1', 0)
        log2fc_t3t2 = row.get('log2FC_T3_vs_T2', 0)

        if pd.notna(log2fc_t2t1):
            bar_width = region_span / 14
            gene_center = (gene_start + gene_end) / 2

            bar_positions = [
                gene_center - bar_width * 1.1,
                gene_center,
                gene_center + bar_width * 1.1,
            ]
            bar_values = [
                log2fc_t2t1 if pd.notna(log2fc_t2t1) else 0,
                log2fc_t3t1 if pd.notna(log2fc_t3t1) else 0,
                log2fc_t3t2 if pd.notna(log2fc_t3t2) else 0,
            ]
            bar_colors = [COL_T2, COL_T3, COL_T3T2]

            bars = ax3.bar(bar_positions, bar_values, width=bar_width * 0.85,
                           color=bar_colors, alpha=0.85, edgecolor='white',
                           linewidth=0.5, zorder=5)

            for bar, val in zip(bars, bar_values):
                if val != 0:
                    offset = 0.25 if val >= 0 else -0.25
                    va = 'bottom' if val >= 0 else 'top'
                    ax3.text(bar.get_x() + bar.get_width() / 2,
                             bar.get_height() + offset,
                             f'{val:.1f}', ha='center', va=va,
                             fontsize=9, fontweight='bold', zorder=6)

    ax3.axhline(0, color='black', linewidth=0.6, zorder=2)
    ax3.axvline(gene_start, color='#9E9E9E', linestyle='--', linewidth=0.7, alpha=0.6)
    ax3.axvline(gene_end, color='#9E9E9E', linestyle='--', linewidth=0.7, alpha=0.6)
    ax3.set_ylabel('log$_2$ fold change')
    ax3.set_xlim(region_start, region_end)
    ax3.set_title('Gene expression', fontsize=11, fontweight='bold', loc='left')
    ax3.tick_params(labelbottom=False)

    legend_elements = [
        mpatches.Patch(facecolor=COL_T2, label='T2 vs T1'),
        mpatches.Patch(facecolor=COL_T3, label='T3 vs T1'),
        mpatches.Patch(facecolor=COL_T3T2, label='T3 vs T2'),
    ]
    ax3.legend(handles=legend_elements, loc='upper right', frameon=True,
               fancybox=False, edgecolor='#BDBDBD', framealpha=0.9)

    # === Track D: Gene annotation ===
    ax4 = axes[3]
    ax4.set_xlim(region_start, region_end)

    # Calculate y-limits based on gene levels
    level_spacing = 1.0
    arrow_h = 0.25
    y_min = -0.5
    y_max = (n_levels - 1) * level_spacing + 0.8

    ax4.set_ylim(y_min, y_max)

    for gene in genes:
        is_target = gene['locus_tag'] == target_info['new_locus']
        y_center = gene['level'] * level_spacing

        draw_gene_arrow(ax4, gene, y_center, arrow_h, is_target, region_span)

        # Gene label
        gene_center = (gene['start'] + gene['end']) / 2
        if gene['gene']:
            label = gene['gene']
            if is_target:
                label = f"{gene['gene']} ({gene['old_locus_tag']})"
        else:
            label = gene['old_locus_tag'] or gene['locus_tag']

        label_y = y_center + arrow_h * 0.7 + 0.15
        ax4.text(gene_center, label_y, label, ha='center', va='bottom',
                 fontsize=9, fontstyle='italic' if gene['gene'] else 'normal',
                 fontweight='bold' if is_target else 'normal',
                 color=COL_TARGET if is_target else '#37474F',
                 zorder=10)

    # Promoter region indicator in gene track
    ax4.axvspan(prom_start, prom_end, alpha=0.25, color=COL_PROMOTER,
                zorder=0, ec='none')

    ax4.set_ylabel('Genes', fontsize=10)
    ax4.set_yticks([])
    ax4.spines['left'].set_visible(False)

    # X-axis: position in kb (clean, no scientific notation)
    ax4.xaxis.set_major_formatter(ticker.FuncFormatter(format_kb_position))
    ax4.set_xlabel('Genomic position (kb)', fontsize=10)

    # Figure title
    fig.suptitle(
        f'{target_name} ({target_info["old_locus"]})',
        fontsize=13, fontweight='bold', y=0.98,
    )

    plt.tight_layout(rect=[0, 0, 1, 0.96])

    # Save
    out_base = OUTPUT_DIR / f'multiomics_track_{target_name}'
    fig.savefig(f'{out_base}.pdf')
    fig.savefig(f'{out_base}.svg')
    fig.savefig(f'{out_base}.png', dpi=300)
    plt.close()
    print(f'  Saved: {out_base}.pdf/.svg/.png')
    return str(out_base)


def main():
    print('=' * 60)
    print('Multi-omics Genomic Track Viewer (Publication Quality)')
    print('=' * 60)

    methyl_df = load_methylation_data()
    expr_df = load_expression_data()
    print(f'Loaded {len(methyl_df)} methylation sites')
    print(f'Loaded {len(expr_df)} genes with expression data')

    for target_name, target_info in TARGETS.items():
        print(f'\nProcessing {target_name} ({target_info["old_locus"]})...')
        create_multiomics_track(target_name, target_info, methyl_df, expr_df)

    print('\n' + '=' * 60)
    print(f'Output: {OUTPUT_DIR}')
    print('=' * 60)


if __name__ == '__main__':
    main()
