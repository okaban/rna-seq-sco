#!/usr/bin/env python3
"""
Act BGC Gene Expression Heatmap
Publication-quality heatmap showing log2FC expression changes
for all 22 genes in the actinorhodin (Act) biosynthetic gene cluster.

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

# === Act BGC gene definitions (SCO5071-SCO5092, genomic order) ===
ACT_GENES = [
    ('SC_RS27515', 'SCO5071', 'actII-ORF1', 'KR (ketoreductase)', '+'),
    ('SC_RS27520', 'SCO5072', 'actII-ORF2', 'ActII-2 (export)', '+'),
    ('SC_RS27525', 'SCO5073', 'actII-ORF3', 'ActII-3 (export)', '+'),
    ('SC_RS27530', 'SCO5074', 'actII-ORF4', 'ActII-ORF4 (SARP)', '+'),
    ('SC_RS27535', 'SCO5075', 'actIII', 'KR (ketoreductase)', '+'),
    ('SC_RS27540', 'SCO5076', 'actI-ORF1', 'KSa (minimal PKS)', '+'),
    ('SC_RS27545', 'SCO5077', 'actI-ORF2', 'CLF (chain length factor)', '+'),
    ('SC_RS27550', 'SCO5078', 'actI-ORF3', 'ACP (acyl carrier protein)', '+'),
    ('SC_RS27555', 'SCO5079', 'actVI-ORF1', 'Cyclase/oxygenase', '+'),
    ('SC_RS27560', 'SCO5080', 'actVI-ORF2', 'Oxygenase', '+'),
    ('SC_RS27565', 'SCO5081', 'actVI-ORFA', 'Cyclase', '+'),
    ('SC_RS27570', 'SCO5082', 'actVA-ORF1', 'Hydroxylase', '+'),
    ('SC_RS27575', 'SCO5083', 'actVA-ORF2', 'Hypothetical', '+'),
    ('SC_RS27580', 'SCO5084', 'actVA-ORF3', 'Hydroxylase', '+'),
    ('SC_RS27585', 'SCO5085', 'actVA-ORF4', 'Hypothetical', '+'),
    ('SC_RS27590', 'SCO5086', 'actVA-ORF5', 'Oxygenase', '+'),
    ('SC_RS27595', 'SCO5087', 'actVA-ORF6', 'Dehydratase', '+'),
    ('SC_RS27600', 'SCO5088', 'actVB', 'Oxidoreductase', '+'),
    ('SC_RS27605', 'SCO5089', '', 'Hypothetical', '+'),
    ('SC_RS27610', 'SCO5090', 'actVI-ORF3', 'Oxidoreductase', '+'),
    ('SC_RS27615', 'SCO5091', 'actVI-ORF4', 'Dehydratase', '+'),
    ('SC_RS27620', 'SCO5092', 'actVII', 'ActVII (export/resistance)', '+'),
]


def main():
    df = pd.read_csv(EXPR_PATH)
    gene_ids = [g[0] for g in ACT_GENES]
    act_df = df[df['gene_id'].isin(gene_ids)].set_index('gene_id')

    comparisons = [
        ('log2FC_T2_vs_T1', 'padj_T2_vs_T1', 'T2 vs T1'),
        ('log2FC_T3_vs_T1', 'padj_T3_vs_T1', 'T3 vs T1'),
        ('log2FC_T3_vs_T2', 'padj_T3_vs_T2', 'T3 vs T2'),
    ]

    n_genes = len(ACT_GENES)
    fc_matrix = np.full((n_genes, 3), np.nan)
    sig_matrix = np.full((n_genes, 3), False)

    gene_labels = []
    for i, (gid, sco, name, func, strand) in enumerate(ACT_GENES):
        if name:
            gene_labels.append(f'{name} ({sco})')
        else:
            gene_labels.append(sco)

        if gid in act_df.index:
            row = act_df.loc[gid]
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

    # Highlight regulatory gene in red
    for i, (_, _, name, _, _) in enumerate(ACT_GENES):
        if name == 'actII-ORF4':
            ax.get_yticklabels()[i].set_fontweight('bold')
            ax.get_yticklabels()[i].set_color('#D32F2F')

    ax.set_title('Act BGC gene expression (log$_2$FC)',
                 fontsize=11, fontweight='bold', pad=10)
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)

    # Colorbar
    cbar = fig.colorbar(im, ax=ax, shrink=0.35, pad=0.03, aspect=12)
    cbar.ax.set_ylabel('log$_2$FC', fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    plt.tight_layout()

    # Save
    out_base = OUTPUT_DIR / 'act_bgc_expression_heatmap'
    fig.savefig(f'{out_base}.pdf')
    fig.savefig(f'{out_base}.svg')
    fig.savefig(f'{out_base}.png', dpi=300)
    plt.close()
    print(f'Saved: {out_base}.pdf/.svg/.png')

    # Summary
    print(f'\nAct BGC expression summary:')
    print(f'  Genes: {n_genes}')
    up_t2 = np.sum((fc_matrix[:, 0] > 1) & sig_matrix[:, 0])
    up_t3 = np.sum((fc_matrix[:, 1] > 1) & sig_matrix[:, 1])
    print(f'  Significantly upregulated (|LFC|>1, padj<0.05):')
    print(f'    T2 vs T1: {up_t2}/{n_genes}')
    print(f'    T3 vs T1: {up_t3}/{n_genes}')
    mean_t2 = np.nanmean(fc_matrix[:, 0])
    mean_t3 = np.nanmean(fc_matrix[:, 1])
    print(f'  Mean LFC: T2/T1={mean_t2:.2f}, T3/T1={mean_t3:.2f}')


if __name__ == '__main__':
    main()
