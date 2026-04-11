#!/usr/bin/env python3
"""
Red BGC Gene Expression Heatmap
Publication-quality heatmap showing log2FC expression changes
for all 22 genes in the undecylprodigiosin (Red) biosynthetic gene cluster.

Streptomyces coelicolor A3(2) M145 Project
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# === Paths ===
BASE = Path('/Users/okaban/bioinfo/rna-seq/11_epigenome_integration')
EXPR_PATH = BASE / 'analysis' / '01_integration' / 'integrated_methyl_expression_weighted.csv'
OUTPUT_DIR = BASE / 'analysis' / '02_publication_figures'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# === Publication style ===
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 9,
    'axes.titlesize': 11,
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 0.7,
})

# === Red BGC gene definitions (SCO5877-SCO5898, genomic order) ===
RED_GENES = [
    ('SC_RS31630', 'SCO5877', 'redD',  'Pathway-specific activator (SARP)', '+'),
    ('SC_RS31635', 'SCO5878', '',       'Proline iminopeptidase', '+'),
    ('SC_RS31640', 'SCO5879', 'redW',   'Oxidoreductase', '+'),
    ('SC_RS31645', 'SCO5880', '',       'Acyl carrier protein', '+'),
    ('SC_RS31650', 'SCO5881', 'redZ',   'Response regulator', '-'),
    ('SC_RS31655', 'SCO5882', '',       'Thioesterase', '-'),
    ('SC_RS31660', 'SCO5883', '',       'Oxidoreductase', '-'),
    ('SC_RS31665', 'SCO5884', '',       'Condensation domain', '-'),
    ('SC_RS31670', 'SCO5885', '',       'Peptidase', '-'),
    ('SC_RS31675', 'SCO5886', '',       'Phosphoenolpyruvate mutase', '-'),
    ('SC_RS31680', 'SCO5887', '',       'Hypothetical protein', '-'),
    ('SC_RS31685', 'SCO5888', '',       'Methyltransferase', '+'),
    ('SC_RS31690', 'SCO5889', '',       'Oxidoreductase', '+'),
    ('SC_RS31695', 'SCO5890', '',       'KS-AT didomain', '+'),
    ('SC_RS31700', 'SCO5891', 'redM',   'Methyltransferase', '+'),
    ('SC_RS31705', 'SCO5892', '',       'PKS (mega-enzyme)', '+'),
    ('SC_RS31710', 'SCO5893', '',       'Condensation', '+'),
    ('SC_RS31715', 'SCO5894', '',       'Acyltransferase', '+'),
    ('SC_RS31720', 'SCO5895', '',       'Peptidyl carrier protein', '+'),
    ('SC_RS31725', 'SCO5896', '',       'Thioester reductase', '+'),
    ('SC_RS31730', 'SCO5897', '',       'Aromatic ring hydroxylase', '+'),
    ('SC_RS31735', 'SCO5898', '',       'Hypothetical protein', '+'),
]


def main():
    df = pd.read_csv(EXPR_PATH)
    gene_ids = [g[0] for g in RED_GENES]
    red_df = df[df['gene_id'].isin(gene_ids)].set_index('gene_id')

    comparisons = [
        ('log2FC_T2_vs_T1', 'padj_T2_vs_T1', 'T2 vs T1'),
        ('log2FC_T3_vs_T1', 'padj_T3_vs_T1', 'T3 vs T1'),
        ('log2FC_T3_vs_T2', 'padj_T3_vs_T2', 'T3 vs T2'),
    ]

    n_genes = len(RED_GENES)
    fc_matrix = np.full((n_genes, 3), np.nan)
    sig_matrix = np.full((n_genes, 3), False)

    gene_labels = []
    for i, (gid, sco, name, func, strand) in enumerate(RED_GENES):
        if name:
            gene_labels.append(f'{name} ({sco})')
        else:
            gene_labels.append(sco)

        if gid in red_df.index:
            row = red_df.loc[gid]
            for j, (fc_col, padj_col, _) in enumerate(comparisons):
                fc_matrix[i, j] = row.get(fc_col, np.nan)
                padj = row.get(padj_col, 1.0)
                sig_matrix[i, j] = pd.notna(padj) and padj < 0.05

    # === Create figure (single heatmap, no methylation sidebar) ===
    fig_w = 100 / 25.4  # 100mm
    fig_h = 160 / 25.4  # 160mm

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # Diverging colormap centered at 0
    vmax = np.nanmax(np.abs(fc_matrix))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    im = ax.imshow(
        fc_matrix, aspect='auto', cmap='RdBu_r', norm=norm,
        interpolation='nearest',
    )

    # Annotate cells
    for i in range(n_genes):
        for j in range(3):
            val = fc_matrix[i, j]
            if np.isnan(val):
                continue
            sig = '**' if sig_matrix[i, j] and abs(val) > 1 else ('*' if sig_matrix[i, j] else '')
            text_color = 'white' if abs(val) > vmax * 0.55 else 'black'
            ax.text(j, i, f'{val:.1f}{sig}', ha='center', va='center',
                    fontsize=7.5, color=text_color, fontweight='bold')

    ax.set_xticks(range(3))
    ax.set_xticklabels([c[2] for c in comparisons], fontsize=9, fontweight='bold')
    ax.set_yticks(range(n_genes))
    ax.set_yticklabels(gene_labels, fontsize=8)

    # Highlight regulatory genes in red
    for i, (_, _, name, _, _) in enumerate(RED_GENES):
        if name in ('redD', 'redZ'):
            ax.get_yticklabels()[i].set_fontweight('bold')
            ax.get_yticklabels()[i].set_color('#D32F2F')

    ax.set_title('Red BGC gene expression (log$_2$FC)',
                 fontsize=11, fontweight='bold', pad=10)
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)

    # Colorbar
    cbar = fig.colorbar(im, ax=ax, shrink=0.35, pad=0.03, aspect=12)
    cbar.ax.set_ylabel('log$_2$FC', fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    plt.tight_layout()

    # Save
    out_base = OUTPUT_DIR / 'red_bgc_expression_heatmap'
    fig.savefig(f'{out_base}.pdf')
    fig.savefig(f'{out_base}.svg')
    fig.savefig(f'{out_base}.png', dpi=300)
    plt.close()
    print(f'Saved: {out_base}.pdf/.svg/.png')

    # Summary
    print(f'\nRed BGC expression summary:')
    print(f'  Genes: {n_genes}')
    up_t2 = np.sum((fc_matrix[:, 0] > 1) & sig_matrix[:, 0])
    up_t3 = np.sum((fc_matrix[:, 1] > 1) & sig_matrix[:, 1])
    print(f'  Significantly upregulated (|LFC|>1, padj<0.05):')
    print(f'    T2 vs T1: {up_t2}/{n_genes}')
    print(f'    T3 vs T1: {up_t3}/{n_genes}')
    mean_t2 = np.nanmean(fc_matrix[:, 0])
    mean_t3 = np.nanmean(fc_matrix[:, 1])
    print(f'  Mean LFC: T2/T1={mean_t2:.2f}, T3/T1={mean_t3:.2f}')
    print(f'  Note: Red BGC genes have 0 promoter methylation sites')
    print(f'        (except SCO5897 with 1 6mA + 1 4mC)')
    print(f'  Epigenetic control is indirect via redZ promoter AAGCCCG methylation')


if __name__ == '__main__':
    main()
