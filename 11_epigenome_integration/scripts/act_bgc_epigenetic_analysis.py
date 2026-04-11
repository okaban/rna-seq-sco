#!/usr/bin/env python3
"""
Act BGC Detailed Epigenetic Analysis
Comprehensive analysis of actinorhodin BGC methylation and expression
Parallel to Red BGC analysis
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

BASE_DIR = Path('/Users/okaban/bioinfo/rna-seq/11_epigenome_integration')
OUTPUT_DIR = BASE_DIR / 'analysis/15_act_bgc_epigenetic'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Act BGC coordinates (from antiSMASH/literature)
ACT_BGC = {
    'name': 'Actinorhodin',
    'start': 5532640,
    'end': 5558199,
    'size_kb': 25.6,
    'genes': {
        'SC_RS27500': {'name': 'actVA-ORF5', 'function': 'Hydroxylase'},
        'SC_RS27505': {'name': 'actVA-ORF4', 'function': 'Monooxygenase'},
        'SC_RS27510': {'name': 'actVA-ORF3', 'function': 'Dehydratase'},
        'SC_RS27515': {'name': 'actVA-ORF2', 'function': 'Unknown'},
        'SC_RS27520': {'name': 'actVA-ORF1', 'function': 'Unknown'},
        'SC_RS27525': {'name': 'actVI-ORF4', 'function': 'Cyclase'},
        'SC_RS27530': {'name': 'actVI-ORF3', 'function': 'Cyclase'},
        'SC_RS27535': {'name': 'actVI-ORF2', 'function': 'Cyclase'},
        'SC_RS27540': {'name': 'actVI-ORF1', 'function': 'Cyclase'},
        'SC_RS27545': {'name': 'actVII', 'function': 'Aromatase'},
        'SC_RS27550': {'name': 'actIV', 'function': 'Cyclase'},
        'SC_RS27555': {'name': 'actI-ORF3', 'function': 'KS'},
        'SC_RS27560': {'name': 'actI-ORF2', 'function': 'CLF'},
        'SC_RS27565': {'name': 'actI-ORF1', 'function': 'ACP'},
        'SC_RS27570': {'name': 'actII-ORF4', 'function': 'SARP regulator'},
        'SC_RS27575': {'name': 'actII-ORF3', 'function': 'Dehydrase'},
        'SC_RS27580': {'name': 'actII-ORF2', 'function': 'ActR (repressor)'},
        'SC_RS27585': {'name': 'actII-ORF1', 'function': 'Transporter'},
        'SC_RS27590': {'name': 'actIII', 'function': 'Ketoreductase'},
        'SC_RS27595': {'name': 'actVB', 'function': 'Hydroxylase'},
    }
}

# Red BGC for comparison
RED_BGC = {
    'name': 'Undecylprodigiosin',
    'start': 5476648,
    'end': 5511923,
    'size_kb': 35.3,
    'csr': 'redD',
    'regulator2': 'redZ',
}

def load_data():
    """Load all required data"""
    # Methylation-expression integration
    methyl_expr = pd.read_csv(BASE_DIR / 'analysis/01_integration/integrated_methyl_expression_weighted.csv')

    # DESeq2 results
    deseq_dir = Path('/Users/okaban/bioinfo/rna-seq/04_deseq2/analysis/04_deseq2_260128_v1/results')
    t2vst1 = pd.read_csv(deseq_dir / 'DESeq2_M145_2_vs_1.tsv', sep='\t')
    t3vst1 = pd.read_csv(deseq_dir / 'DESeq2_M145_3_vs_1.tsv', sep='\t')

    # TF methylation
    tf_methyl = pd.read_csv(BASE_DIR / 'analysis/12_grn_tf_methylation/tf_promoter_methylation/tf_methylation_summary.csv')

    # Coordinated genes
    coord_genes = pd.read_csv(BASE_DIR / 'analysis/01_integration/T2vsT1_coordinated_genes.csv')

    return methyl_expr, t2vst1, t3vst1, tf_methyl, coord_genes

def analyze_act_bgc_genes(methyl_expr, t2vst1, t3vst1):
    """Detailed analysis of Act BGC genes"""
    results = []

    for locus, info in ACT_BGC['genes'].items():
        # Get methylation data
        methyl = methyl_expr[methyl_expr['gene_id'] == locus]

        # Get expression data
        expr_t2 = t2vst1[t2vst1['gene_id'] == locus]
        expr_t3 = t3vst1[t3vst1['gene_id'] == locus]

        result = {
            'locus_tag': locus,
            'gene_name': info['name'],
            'function': info['function'],
        }

        # Methylation
        if len(methyl) > 0:
            m = methyl.iloc[0]
            result['6mA_T1'] = m['6mA_T1_count']
            result['6mA_T2'] = m['6mA_T2_count']
            result['4mC_T1'] = m['4mC_T1_count']
            result['4mC_T2'] = m['4mC_T2_count']
            result['total_T1'] = m['6mA_T1_count'] + m['4mC_T1_count']
            result['total_T2'] = m['6mA_T2_count'] + m['4mC_T2_count']
        else:
            result['6mA_T1'] = result['6mA_T2'] = 0
            result['4mC_T1'] = result['4mC_T2'] = 0
            result['total_T1'] = result['total_T2'] = 0

        # Expression T2vsT1
        if len(expr_t2) > 0:
            e = expr_t2.iloc[0]
            result['log2FC_T2vsT1'] = e['log2FoldChange']
            result['padj_T2vsT1'] = e['padj']
        else:
            result['log2FC_T2vsT1'] = np.nan
            result['padj_T2vsT1'] = np.nan

        # Expression T3vsT1
        if len(expr_t3) > 0:
            e = expr_t3.iloc[0]
            result['log2FC_T3vsT1'] = e['log2FoldChange']
            result['padj_T3vsT1'] = e['padj']
        else:
            result['log2FC_T3vsT1'] = np.nan
            result['padj_T3vsT1'] = np.nan

        # Classify methylation change
        if result['total_T2'] > result['total_T1']:
            result['methyl_change'] = 'Gained'
        elif result['total_T2'] < result['total_T1']:
            result['methyl_change'] = 'Lost'
        elif result['total_T2'] == 0:
            result['methyl_change'] = 'None'
        else:
            result['methyl_change'] = 'Stable'

        results.append(result)

    return pd.DataFrame(results)

def analyze_actii_orf4(tf_methyl, t2vst1):
    """Detailed analysis of actII-ORF4 (SARP regulator)"""
    # Get TF data
    actii = tf_methyl[tf_methyl['name'] == 'actII-ORF4']

    # Get expression
    expr = t2vst1[t2vst1['gene_id'] == 'SC_RS27570']

    analysis = {
        'locus_tag': 'SC_RS27570',
        'gene_name': 'actII-ORF4',
        'function': 'SARP family transcriptional activator',
        'target': 'Act biosynthetic genes',
    }

    if len(actii) > 0:
        a = actii.iloc[0]
        analysis['methyl_T1'] = a['total_T1_sites']
        analysis['methyl_T2'] = a['total_T2_sites']
        analysis['methyl_status'] = a['methyl_status_T2vsT1']
        analysis['corr_act'] = a['corr_act']
        analysis['corr_red'] = a['corr_red']

    if len(expr) > 0:
        e = expr.iloc[0]
        analysis['log2FC'] = e['log2FoldChange']
        analysis['padj'] = e['padj']

    return analysis

def compare_act_vs_red(act_results, methyl_expr, t2vst1, coord_genes):
    """Compare Act and Red BGC epigenetic control"""
    # Red BGC genes
    red_genes = [f'SC_RS272{i:02d}' for i in range(10, 65, 5)]
    red_results = []

    for locus in red_genes:
        methyl = methyl_expr[methyl_expr['gene_id'] == locus]
        expr = t2vst1[t2vst1['gene_id'] == locus]

        result = {'locus_tag': locus, 'bgc': 'Red'}

        if len(methyl) > 0:
            m = methyl.iloc[0]
            result['total_T1'] = m['6mA_T1_count'] + m['4mC_T1_count']
            result['total_T2'] = m['6mA_T2_count'] + m['4mC_T2_count']
        else:
            result['total_T1'] = result['total_T2'] = 0

        if len(expr) > 0:
            result['log2FC'] = expr.iloc[0]['log2FoldChange']
        else:
            result['log2FC'] = np.nan

        red_results.append(result)

    red_df = pd.DataFrame(red_results)

    # Add BGC label to Act
    act_df = act_results[['locus_tag', 'total_T1', 'total_T2', 'log2FC_T2vsT1']].copy()
    act_df.columns = ['locus_tag', 'total_T1', 'total_T2', 'log2FC']
    act_df['bgc'] = 'Act'

    # Combine
    comparison = pd.concat([act_df, red_df], ignore_index=True)

    # Statistics
    act_methylated = (act_df['total_T2'] > 0).sum()
    red_methylated = (red_df['total_T2'] > 0).sum()

    act_mean_fc = act_df['log2FC'].mean()
    red_mean_fc = red_df['log2FC'].mean()

    # Coordinated genes in each BGC
    coord_locus = set(coord_genes['gene_id'].values)
    act_coord = act_df['locus_tag'].isin(coord_locus).sum()
    red_coord = red_df['locus_tag'].isin(coord_locus).sum()

    return {
        'act_genes': len(act_df),
        'red_genes': len(red_df),
        'act_methylated': act_methylated,
        'red_methylated': red_methylated,
        'act_mean_fc': act_mean_fc,
        'red_mean_fc': red_mean_fc,
        'act_coordinated': act_coord,
        'red_coordinated': red_coord,
        'comparison_df': comparison,
    }

def create_visualizations(act_results, actii_analysis, comparison):
    """Create comprehensive visualizations"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # 1. Act BGC gene expression
    ax1 = axes[0, 0]
    act_sorted = act_results.sort_values('log2FC_T2vsT1')
    colors = ['#59a14f' if x > 0 else '#e15759' for x in act_sorted['log2FC_T2vsT1'].fillna(0)]
    ax1.barh(act_sorted['gene_name'], act_sorted['log2FC_T2vsT1'].fillna(0), color=colors)
    ax1.axvline(x=0, color='black', linewidth=0.5)
    ax1.set_xlabel('log2 Fold Change (T2 vs T1)')
    ax1.set_title('Act BGC Gene Expression')

    # Highlight actII-ORF4
    actii_idx = act_sorted[act_sorted['gene_name'] == 'actII-ORF4'].index
    if len(actii_idx) > 0:
        ax1.get_children()[list(act_sorted['gene_name']).index('actII-ORF4')].set_edgecolor('blue')
        ax1.get_children()[list(act_sorted['gene_name']).index('actII-ORF4')].set_linewidth(2)

    # 2. Methylation status
    ax2 = axes[0, 1]
    methyl_counts = act_results['methyl_change'].value_counts()
    colors = {'None': '#e0e0e0', 'Gained': '#59a14f', 'Lost': '#e15759', 'Stable': '#76b7b2'}
    ax2.pie(methyl_counts.values, labels=methyl_counts.index, autopct='%1.1f%%',
            colors=[colors.get(x, '#999999') for x in methyl_counts.index])
    ax2.set_title('Act BGC Methylation Status (T2 vs T1)')

    # 3. Act vs Red comparison
    ax3 = axes[1, 0]
    comp_df = comparison['comparison_df']

    act_fc = comp_df[comp_df['bgc'] == 'Act']['log2FC'].dropna()
    red_fc = comp_df[comp_df['bgc'] == 'Red']['log2FC'].dropna()

    bp = ax3.boxplot([act_fc, red_fc], labels=['Act BGC', 'Red BGC'], patch_artist=True)
    bp['boxes'][0].set_facecolor('#4e79a7')
    bp['boxes'][1].set_facecolor('#f28e2b')
    ax3.axhline(y=0, color='gray', linestyle='--')
    ax3.set_ylabel('log2 Fold Change (T2 vs T1)')
    ax3.set_title('Act vs Red BGC Expression Comparison')

    # Add stats
    if len(act_fc) > 0 and len(red_fc) > 0:
        stat, pval = stats.mannwhitneyu(act_fc, red_fc, alternative='two-sided')
        ax3.annotate(f'Mann-Whitney p = {pval:.3f}', (0.5, 0.95), xycoords='axes fraction', ha='center')

    # 4. Summary table
    ax4 = axes[1, 1]
    ax4.axis('off')

    summary_data = [
        ['Metric', 'Act BGC', 'Red BGC'],
        ['Total genes', str(comparison['act_genes']), str(comparison['red_genes'])],
        ['Methylated (T2)', str(comparison['act_methylated']), str(comparison['red_methylated'])],
        ['Coordinated genes', str(comparison['act_coordinated']), str(comparison['red_coordinated'])],
        ['Mean log2FC', f"{comparison['act_mean_fc']:.2f}", f"{comparison['red_mean_fc']:.2f}"],
        ['SARP methylation', actii_analysis['methyl_status'], 'redZ: Lost (coordinated)'],
    ]

    table = ax4.table(cellText=summary_data, loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.8)

    # Style header
    for i in range(3):
        table[(0, i)].set_facecolor('#4472C4')
        table[(0, i)].set_text_props(color='white', fontweight='bold')

    ax4.set_title('Act vs Red BGC Summary')

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'act_bgc_epigenetic_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: act_bgc_epigenetic_analysis.png")

def main():
    print("="*60)
    print("ACT BGC DETAILED EPIGENETIC ANALYSIS")
    print("="*60)

    # Load data
    print("\n1. Loading data...")
    methyl_expr, t2vst1, t3vst1, tf_methyl, coord_genes = load_data()

    # Analyze Act BGC genes
    print("\n2. Analyzing Act BGC genes...")
    act_results = analyze_act_bgc_genes(methyl_expr, t2vst1, t3vst1)
    print(f"   Found {len(act_results)} Act BGC genes")
    act_results.to_csv(OUTPUT_DIR / 'act_bgc_gene_analysis.csv', index=False)

    # Analyze actII-ORF4
    print("\n3. Analyzing actII-ORF4 (SARP regulator)...")
    actii_analysis = analyze_actii_orf4(tf_methyl, t2vst1)
    print(f"   Methylation: T1={actii_analysis['methyl_T1']:.0f}, T2={actii_analysis['methyl_T2']:.0f}")
    print(f"   Status: {actii_analysis['methyl_status']}")
    print(f"   Expression: log2FC = {actii_analysis['log2FC']:.2f}")
    print(f"   Correlation with Act: r = {actii_analysis['corr_act']:.3f}")

    # Compare Act vs Red
    print("\n4. Comparing Act vs Red BGC...")
    comparison = compare_act_vs_red(act_results, methyl_expr, t2vst1, coord_genes)
    print(f"   Act: {comparison['act_methylated']}/{comparison['act_genes']} methylated, mean FC = {comparison['act_mean_fc']:.2f}")
    print(f"   Red: {comparison['red_methylated']}/{comparison['red_genes']} methylated, mean FC = {comparison['red_mean_fc']:.2f}")

    # Save comparison
    comparison['comparison_df'].to_csv(OUTPUT_DIR / 'act_vs_red_comparison.csv', index=False)

    # Create visualizations
    print("\n5. Creating visualizations...")
    create_visualizations(act_results, actii_analysis, comparison)

    # Summary
    print("\n" + "="*60)
    print("KEY FINDINGS")
    print("="*60)
    print(f"\n1. actII-ORF4 (SARP regulator):")
    print(f"   - Methylation: {actii_analysis['methyl_status']}")
    print(f"   - Expression: log2FC = {actii_analysis['log2FC']:.2f}")
    print(f"   - Strong correlation with Act BGC (r = {actii_analysis['corr_act']:.3f})")

    print(f"\n2. Act vs Red epigenetic control:")
    print(f"   - Act: {comparison['act_methylated']} genes with methylation, {comparison['act_coordinated']} coordinated")
    print(f"   - Red: {comparison['red_methylated']} genes with methylation, {comparison['red_coordinated']} coordinated")
    print(f"   - Key difference: Red has redZ (coordinated methylation), Act has no equivalent")

    print(f"\n3. Act BGC expression:")
    print(f"   - Mean log2FC = {comparison['act_mean_fc']:.2f} (T2 vs T1)")
    print(f"   - Mixed pattern (both up and down regulated genes)")

    print(f"\nOutput directory: {OUTPUT_DIR}")

if __name__ == '__main__':
    main()
