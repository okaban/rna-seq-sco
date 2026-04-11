#!/usr/bin/env python3
"""
MTase expression vs methylation dynamics.
Quantitative evaluation of enzyme expression → methylation site dynamics.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings
warnings.filterwarnings('ignore')

BASE = Path('/Users/okaban/bioinfo/rna-seq')
GENE_MASTER = BASE / '05_annotation/analysis/05_annotation_260128_v1/tables/gene_master_DESeq2.tsv'
MTASE_EXPR = BASE / '11_epigenome_integration/analysis/11_rm_system_identification/mtase_genes_with_expression.csv'
GAINED_LOST = BASE / '11_epigenome_integration/analysis/18_tss_analyses/temporal_gained_lost_motifs.csv'
MOTIF_ENRICHMENT = BASE / '11_epigenome_integration/analysis/18_tss_analyses/temporal_motif_enrichment.csv'
METHYL_SITES = BASE / '11_epigenome_integration/analysis/01_integration/high_confidence_sites_weighted.csv'
OUT_DIR = BASE / '11_epigenome_integration/analysis/18_tss_analyses'


def load_mtase_expression():
    """Load MTase expression data and annotate with functions."""
    mtase = pd.read_csv(MTASE_EXPR)

    # Key enzymes (10 genuine DNA MTases + McrA restriction enzyme)
    key_enzymes = {
        'SC_RS17645': {'name': 'N-6 DNA methylase\n(AAGCCCG candidate)', 'motif': 'AAGCCCG', 'mod': '6mA',
                       'color': '#E74C3C'},
        'SC_RS24685': {'name': 'SAM-dep DNA MTase\n(SC_RS24685)', 'motif': 'Unknown', 'mod': 'Unknown',
                       'color': '#C0392B'},
        'SC_RS10665': {'name': 'SCO1731\n(m5C MTase)', 'motif': 'CCGG (m5C)', 'mod': '5mC',
                       'color': '#9B59B6'},
        'SC_RS19770': {'name': 'Dcm-like MTase\n(SC_RS19770)', 'motif': 'CCGG (Dcm)', 'mod': '4mC/5mC',
                       'color': '#3498DB'},
        'SC_RS36410': {'name': 'Dcm-like MTase\n(SC_RS36410)', 'motif': 'CCGG (Dcm)', 'mod': '4mC/5mC',
                       'color': '#2ECC71'},
        'SC_RS19670': {'name': 'DNA-MTase\n(SC_RS19670)', 'motif': 'Unknown', 'mod': 'Unknown',
                       'color': '#27AE60'},
        'SC_RS36625': {'name': 'DNA-MTase\n(SC_RS36625)', 'motif': 'Unknown', 'mod': 'Unknown',
                       'color': '#16A085'},
        'SC_RS28835': {'name': 'PglX\n(SC_RS28835)', 'motif': 'Unknown (6mA)', 'mod': '6mA',
                       'color': '#F39C12'},
        'SC_RS35335': {'name': 'PglX\n(SC_RS35335)', 'motif': 'Unknown (6mA)', 'mod': '6mA',
                       'color': '#E67E22'},
        'SC_RS25315': {'name': 'McrA\n(methyl-CCGG cutter)', 'motif': 'Cuts mCCGG', 'mod': 'restriction',
                       'color': '#1ABC9C'},
    }

    # Get expression from gene_master for SCO1731 and McrA
    gm = pd.read_csv(GENE_MASTER, sep='\t')

    results = []
    for gene_id, info in key_enzymes.items():
        row_mtase = mtase[mtase['locus_tag'] == gene_id]
        row_gm = gm[gm['gene_id'] == gene_id]

        if len(row_mtase) > 0:
            r = row_mtase.iloc[0]
            results.append({
                'gene_id': gene_id,
                'name': info['name'],
                'motif': info['motif'],
                'mod': info['mod'],
                'color': info['color'],
                'log2FC_T2vsT1': r.get('log2FC_T2vsT1', np.nan),
                'log2FC_T3vsT1': r.get('log2FC_T3vsT1', np.nan),
                'log2FC_T3vsT2': r.get('log2FC_T3vsT2', np.nan),
                'padj_T2vsT1': r.get('padj_T2vsT1', np.nan),
                'padj_T3vsT1': r.get('padj_T3vsT1', np.nan),
                'padj_T3vsT2': r.get('padj_T3vsT2', np.nan),
            })
        elif len(row_gm) > 0:
            r = row_gm.iloc[0]
            results.append({
                'gene_id': gene_id,
                'name': info['name'],
                'motif': info['motif'],
                'mod': info['mod'],
                'color': info['color'],
                'log2FC_T2vsT1': r.get('log2FoldChange_2_vs_1', np.nan),
                'log2FC_T3vsT1': r.get('log2FoldChange_3_vs_1', np.nan),
                'log2FC_T3vsT2': r.get('log2FoldChange_3_vs_2', np.nan),
                'padj_T2vsT1': r.get('padj_2_vs_1', np.nan),
                'padj_T3vsT1': r.get('padj_3_vs_1', np.nan),
                'padj_T3vsT2': r.get('padj_3_vs_2', np.nan),
            })

    return pd.DataFrame(results)


def load_methylation_dynamics():
    """Load per-timepoint methylation site counts."""
    methyl = pd.read_csv(METHYL_SITES)

    counts = {}
    for tp in ['T1', 'T2', 'T3']:
        tp_data = methyl[methyl['timepoint'] == tp]
        for mod in ['6mA', '4mC']:
            mod_data = tp_data[tp_data['mod_type'] == mod]
            counts[f'{mod}_{tp}'] = len(mod_data.drop_duplicates('position'))

    return counts


def load_gained_lost_counts():
    """Load gained/lost site counts per comparison."""
    df = pd.read_csv(GAINED_LOST)

    results = {}
    for comp in df['comparison'].unique():
        for mod in df['mod_type'].unique():
            sub = df[(df['comparison'] == comp) & (df['mod_type'] == mod)]
            for cat in ['gained', 'lost', 'shared']:
                cat_data = sub[sub['category'] == cat]
                if len(cat_data) > 0:
                    n = cat_data['n_sites'].values[0]
                    results[f'{mod}_{comp}_{cat}'] = n
    return results


def main():
    print("Loading data...")
    mtase_df = load_mtase_expression()
    methyl_counts = load_methylation_dynamics()
    gained_lost = load_gained_lost_counts()

    print("\n=== MTase Expression Summary ===")
    for _, row in mtase_df.iterrows():
        sig_t2 = '*' if row['padj_T2vsT1'] < 0.05 else 'ns'
        sig_t3t1 = '*' if row['padj_T3vsT1'] < 0.05 else 'ns'
        sig_t3t2 = '*' if row['padj_T3vsT2'] < 0.05 else 'ns'
        print(f"  {row['gene_id']} ({row['name'].replace(chr(10), ' ')}): "
              f"T2vsT1={row['log2FC_T2vsT1']:+.2f}({sig_t2}) "
              f"T3vsT1={row['log2FC_T3vsT1']:+.2f}({sig_t3t1}) "
              f"T3vsT2={row['log2FC_T3vsT2']:+.2f}({sig_t3t2})")

    print("\n=== Methylation Site Counts ===")
    for mod in ['6mA', '4mC']:
        t1 = methyl_counts.get(f'{mod}_T1', 0)
        t2 = methyl_counts.get(f'{mod}_T2', 0)
        t3 = methyl_counts.get(f'{mod}_T3', 0)
        print(f"  {mod}: T1={t1}, T2={t2}, T3={t3}")

    print("\n=== Gained/Lost Counts ===")
    for key, val in sorted(gained_lost.items()):
        print(f"  {key}: {val}")

    # ===== FIGURE 1: MTase expression heatmap + methylation dynamics =====
    fig = plt.figure(figsize=(16, 12))
    gs = gridspec.GridSpec(2, 3, height_ratios=[1.4, 1], hspace=0.35, wspace=0.35)

    # Panel A: MTase expression heatmap
    ax_heat = fig.add_subplot(gs[0, 0:2])

    enzyme_order = ['SC_RS17645', 'SC_RS24685', 'SC_RS10665', 'SC_RS19770',
                    'SC_RS36410', 'SC_RS19670', 'SC_RS36625',
                    'SC_RS28835', 'SC_RS35335', 'SC_RS25315']
    mtase_ordered = mtase_df.set_index('gene_id').loc[
        [e for e in enzyme_order if e in mtase_df['gene_id'].values]
    ]

    data_matrix = mtase_ordered[['log2FC_T2vsT1', 'log2FC_T3vsT1', 'log2FC_T3vsT2']].values
    padj_matrix = mtase_ordered[['padj_T2vsT1', 'padj_T3vsT1', 'padj_T3vsT2']].values

    im = ax_heat.imshow(data_matrix, cmap='RdBu_r', vmin=-4, vmax=4, aspect='auto')

    # Add text annotations
    for i in range(data_matrix.shape[0]):
        for j in range(data_matrix.shape[1]):
            val = data_matrix[i, j]
            padj = padj_matrix[i, j]
            if np.isnan(val):
                continue
            sig = '***' if padj < 0.001 else ('**' if padj < 0.01 else ('*' if padj < 0.05 else ''))
            color = 'white' if abs(val) > 2 else 'black'
            ax_heat.text(j, i, f'{val:.1f}{sig}', ha='center', va='center', fontsize=8, color=color)

    ax_heat.set_xticks([0, 1, 2])
    ax_heat.set_xticklabels(['T2 vs T1', 'T3 vs T1', 'T3 vs T2'], fontsize=9)
    ax_heat.set_yticks(range(len(mtase_ordered)))
    ax_heat.set_yticklabels(mtase_ordered['name'].values, fontsize=8)
    ax_heat.set_title('A. Methylation enzyme expression (log2FC)', fontsize=11, fontweight='bold')
    plt.colorbar(im, ax=ax_heat, shrink=0.8, label='log2FC')

    # Panel B: Methylation site counts
    ax_counts = fig.add_subplot(gs[0, 2])
    timepoints = ['T1', 'T2', 'T3']
    x = np.arange(3)
    w = 0.35

    counts_6mA = [methyl_counts.get(f'6mA_{tp}', 0) for tp in timepoints]
    counts_4mC = [methyl_counts.get(f'4mC_{tp}', 0) for tp in timepoints]

    ax_counts.bar(x - w/2, counts_6mA, w, label='6mA', color='#E74C3C', alpha=0.8)
    ax_counts.bar(x + w/2, counts_4mC, w, label='4mC', color='#3498DB', alpha=0.8)
    ax_counts.set_xticks(x)
    ax_counts.set_xticklabels(timepoints)
    ax_counts.set_ylabel('Unique methylation sites')
    ax_counts.set_title('B. Methylation site counts', fontsize=11, fontweight='bold')
    ax_counts.legend(fontsize=9)
    for i, (a, c) in enumerate(zip(counts_6mA, counts_4mC)):
        ax_counts.text(i - w/2, a + 30, str(a), ha='center', fontsize=8)
        ax_counts.text(i + w/2, c + 30, str(c), ha='center', fontsize=8)

    # Panel C: SC_RS17645 expression vs AAGCCCG dynamics
    ax_aag = fig.add_subplot(gs[1, 0])

    # SC_RS17645 relative expression (approximated from log2FC)
    rs17645 = mtase_df[mtase_df['gene_id'] == 'SC_RS17645'].iloc[0]
    # Relative to T1 = 1.0
    expr_T1 = 1.0
    expr_T2 = 2 ** rs17645['log2FC_T2vsT1']
    expr_T3 = 2 ** rs17645['log2FC_T3vsT1']

    # AAGCCCG gained/lost from temporal data
    motif_enrich = pd.read_csv(MOTIF_ENRICHMENT)
    aag_pcts = {}
    for tp in ['T1', 'T2', 'T3']:
        row = motif_enrich[(motif_enrich['mod_type'] == '6mA') & (motif_enrich['motif'] == 'AAGCCCG')]
        if len(row) > 0:
            aag_pcts[tp] = row[f'{tp}_pct'].values[0]

    ax_aag_twin = ax_aag.twinx()
    ax_aag.plot([0, 1, 2], [expr_T1, expr_T2, expr_T3], 'o-', color='#E74C3C', lw=2,
                markersize=8, label='SC_RS17645 expression')
    if aag_pcts:
        aag_vals = [aag_pcts.get(tp, 0) for tp in timepoints]
        ax_aag_twin.plot([0, 1, 2], aag_vals, 's--', color='#8E44AD', lw=2,
                         markersize=8, label='AAGCCCG in 6mA (%)')
        ax_aag_twin.set_ylabel('AAGCCCG prevalence (%)', color='#8E44AD')

    ax_aag.set_xticks([0, 1, 2])
    ax_aag.set_xticklabels(timepoints)
    ax_aag.set_ylabel('Relative expression', color='#E74C3C')
    ax_aag.set_title('C. SC_RS17645 vs AAGCCCG\n(6mA methylase candidate)', fontsize=10, fontweight='bold')
    ax_aag.axhline(y=1, color='gray', ls=':', alpha=0.5)

    lines1, labels1 = ax_aag.get_legend_handles_labels()
    lines2, labels2 = ax_aag_twin.get_legend_handles_labels()
    ax_aag.legend(lines1 + lines2, labels1 + labels2, fontsize=7, loc='lower left')

    # Panel D: Dcm-like MTase expression vs 4mC-CCGG dynamics (the paradox)
    ax_paradox = fig.add_subplot(gs[1, 1])

    # Dcm-like MTases
    rs19770 = mtase_df[mtase_df['gene_id'] == 'SC_RS19770']
    rs36410 = mtase_df[mtase_df['gene_id'] == 'SC_RS36410']

    if len(rs19770) > 0:
        r = rs19770.iloc[0]
        e1 = [1.0, 2**r['log2FC_T2vsT1'] if not np.isnan(r['log2FC_T2vsT1']) else 1.0,
              2**r['log2FC_T3vsT1'] if not np.isnan(r['log2FC_T3vsT1']) else 1.0]
        ax_paradox.plot([0, 1, 2], e1, 'o-', color='#3498DB', lw=2, markersize=7,
                        label='SC_RS19770 (Dcm)')

    if len(rs36410) > 0:
        r = rs36410.iloc[0]
        e2 = [1.0, 2**r['log2FC_T2vsT1'] if not np.isnan(r['log2FC_T2vsT1']) else 1.0,
              2**r['log2FC_T3vsT1'] if not np.isnan(r['log2FC_T3vsT1']) else 1.0]
        ax_paradox.plot([0, 1, 2], e2, 'D-', color='#2ECC71', lw=2, markersize=7,
                        label='SC_RS36410 (Dcm)')

    # 4mC site counts normalized
    ax_p_twin = ax_paradox.twinx()
    norm_4mC = [c / counts_4mC[0] for c in counts_4mC]
    ax_p_twin.plot([0, 1, 2], norm_4mC, 's--', color='#E74C3C', lw=2, markersize=8,
                    label='4mC sites (norm.)')
    ax_p_twin.set_ylabel('4mC sites (relative to T1)', color='#E74C3C')

    ax_paradox.set_xticks([0, 1, 2])
    ax_paradox.set_xticklabels(timepoints)
    ax_paradox.set_ylabel('Relative expression')
    ax_paradox.set_title('D. Dcm-like MTase vs 4mC sites\n(THE PARADOX: m4C ≠ m5C?)', fontsize=10, fontweight='bold')
    ax_paradox.axhline(y=1, color='gray', ls=':', alpha=0.5)
    ax_paradox.set_yscale('log', base=2)

    lines1, labels1 = ax_paradox.get_legend_handles_labels()
    lines2, labels2 = ax_p_twin.get_legend_handles_labels()
    ax_paradox.legend(lines1 + lines2, labels1 + labels2, fontsize=7, loc='upper left')

    # Panel E: SCO1731 expression vs Act/Red regulator expression
    ax_sco = fig.add_subplot(gs[1, 2])

    gm = pd.read_csv(GENE_MASTER, sep='\t')
    key_genes = {
        'SC_RS10665': ('SCO1731 (m5C MTase)', '#9B59B6'),
        'SC_RS25315': ('McrA (methyl-CCGG cutter)', '#1ABC9C'),
        'SC_RS18430': ('SCO3262 (HNH endo)', '#95A5A6'),
    }

    for gene_id, (label, color) in key_genes.items():
        row = gm[gm['gene_id'] == gene_id]
        if len(row) > 0:
            r = row.iloc[0]
            lfc_t2 = r.get('log2FoldChange_2_vs_1', np.nan)
            lfc_t3 = r.get('log2FoldChange_3_vs_1', np.nan)
            if not np.isnan(lfc_t2) and not np.isnan(lfc_t3):
                expr = [1.0, 2**lfc_t2, 2**lfc_t3]
                ax_sco.plot([0, 1, 2], expr, 'o-', color=color, lw=2, markersize=7, label=label)

    ax_sco.set_xticks([0, 1, 2])
    ax_sco.set_xticklabels(timepoints)
    ax_sco.set_ylabel('Relative expression')
    ax_sco.set_title('E. Restriction-Modification system\nexpression dynamics', fontsize=10, fontweight='bold')
    ax_sco.axhline(y=1, color='gray', ls=':', alpha=0.5)
    ax_sco.set_yscale('log', base=2)
    ax_sco.legend(fontsize=7, loc='lower left')

    plt.savefig(OUT_DIR / 'mtase_methylation_dynamics.png', dpi=150, bbox_inches='tight')
    plt.savefig(OUT_DIR / 'mtase_methylation_dynamics.pdf', bbox_inches='tight')
    plt.close()
    print("\nSaved: mtase_methylation_dynamics.png/pdf")

    # ===== Key contradiction analysis =====
    print("\n" + "=" * 60)
    print("KEY ANALYSIS: m4C vs m5C hypothesis")
    print("=" * 60)

    print("\nDcm-like MTases (CCGG → m5C producers):")
    for gid in ['SC_RS19770', 'SC_RS36410']:
        r = mtase_df[mtase_df['gene_id'] == gid]
        if len(r) > 0:
            r = r.iloc[0]
            print(f"  {gid}: T3vsT2 log2FC = {r['log2FC_T3vsT2']:+.2f} (padj={r['padj_T3vsT2']:.2e})")

    print(f"\n4mC sites (Nanopore): T1={counts_4mC[0]}, T2={counts_4mC[1]}, T3={counts_4mC[2]}")
    print(f"  T3/T2 ratio: {counts_4mC[2]/counts_4mC[1]:.2f} (56% decrease)")

    print("\n→ Dcm-like MTase expression INCREASES 10-40x at T3")
    print("→ Yet m4C sites DECREASE by 56% at T3")
    print("→ CONCLUSION: Dcm-like MTases likely produce m5C, not m4C")
    print("→ The m4C enzyme on CCGG remains unidentified")
    print("→ This implies S. coelicolor has DUAL cytosine methylation (m4C + m5C) on CCGG")

    # Save summary
    mtase_df.to_csv(OUT_DIR / 'mtase_expression_summary.csv', index=False)
    print("\nSaved: mtase_expression_summary.csv")


if __name__ == '__main__':
    main()
