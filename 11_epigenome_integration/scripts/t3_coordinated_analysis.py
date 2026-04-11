#!/usr/bin/env python3
"""
T3 Coordinated Gene Analysis
Generate T3vsT1 and T3vsT2 coordinated gene lists
Compare with T2vsT1 patterns
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
OUTPUT_DIR = BASE_DIR / 'analysis/16_t3_coordinated'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def load_data():
    """Load methylation and expression data"""
    # Integrated data (weighted method for T3)
    methyl_expr = pd.read_csv(BASE_DIR / 'analysis/01_integration/integrated_methyl_expression_weighted.csv')

    # DESeq2 results
    deseq_dir = Path('/Users/okaban/bioinfo/rna-seq/04_deseq2/analysis/04_deseq2_260128_v1/results')
    t2vst1 = pd.read_csv(deseq_dir / 'DESeq2_M145_2_vs_1.tsv', sep='\t')
    t3vst1 = pd.read_csv(deseq_dir / 'DESeq2_M145_3_vs_1.tsv', sep='\t')
    t3vst2 = pd.read_csv(deseq_dir / 'DESeq2_M145_3_vs_2.tsv', sep='\t')

    # T2vsT1 coordinated genes (for comparison)
    t2_coord = pd.read_csv(BASE_DIR / 'analysis/01_integration/T2vsT1_coordinated_genes.csv')

    return methyl_expr, t2vst1, t3vst1, t3vst2, t2_coord

def identify_coordinated_genes(methyl_expr, deseq_results, comparison_name, methyl_col_t1, methyl_col_t2):
    """Identify genes with coordinated methylation-expression changes"""
    # Merge methylation and expression
    merged = methyl_expr.merge(
        deseq_results[['gene_id', 'log2FoldChange', 'padj']],
        on='gene_id',
        how='inner'
    )

    # Calculate methylation change
    merged['6mA_change'] = merged[f'6mA_{methyl_col_t2}_mean_freq'] - merged[f'6mA_{methyl_col_t1}_mean_freq']
    merged['4mC_change'] = merged[f'4mC_{methyl_col_t2}_mean_freq'] - merged[f'4mC_{methyl_col_t1}_mean_freq']

    # Define coordinated criteria
    # 1. Significant expression change (|log2FC| > 1, padj < 0.05)
    # 2. Methylation change in same direction

    coordinated = []

    for _, row in merged.iterrows():
        if pd.isna(row['log2FoldChange']) or pd.isna(row['padj']):
            continue

        expr_sig = abs(row['log2FoldChange']) > 1 and row['padj'] < 0.05
        if not expr_sig:
            continue

        expr_up = row['log2FoldChange'] > 1
        expr_down = row['log2FoldChange'] < -1

        # Check 6mA
        if row['6mA_change'] > 10:  # Gained
            if expr_up:
                coordinated.append({
                    'gene_id': row['gene_id'],
                    'mod_type': '6mA',
                    'methyl_change': row['6mA_change'],
                    'methyl_category': 'Gained',
                    'log2FC': row['log2FoldChange'],
                    'padj': row['padj'],
                    'coordination': 'Gained_Up',
                    'correlation_type': 'Positive',
                })
        elif row['6mA_change'] < -10:  # Lost
            if expr_down:
                coordinated.append({
                    'gene_id': row['gene_id'],
                    'mod_type': '6mA',
                    'methyl_change': row['6mA_change'],
                    'methyl_category': 'Lost',
                    'log2FC': row['log2FoldChange'],
                    'padj': row['padj'],
                    'coordination': 'Lost_Down',
                    'correlation_type': 'Positive',
                })
            elif expr_up:
                coordinated.append({
                    'gene_id': row['gene_id'],
                    'mod_type': '6mA',
                    'methyl_change': row['6mA_change'],
                    'methyl_category': 'Lost',
                    'log2FC': row['log2FoldChange'],
                    'padj': row['padj'],
                    'coordination': 'Lost_Up',
                    'correlation_type': 'Negative',
                })

        # Check 4mC
        if row['4mC_change'] > 10:  # Gained
            if expr_up:
                coordinated.append({
                    'gene_id': row['gene_id'],
                    'mod_type': '4mC',
                    'methyl_change': row['4mC_change'],
                    'methyl_category': 'Gained',
                    'log2FC': row['log2FoldChange'],
                    'padj': row['padj'],
                    'coordination': 'Gained_Up',
                    'correlation_type': 'Positive',
                })
        elif row['4mC_change'] < -10:  # Lost
            if expr_down:
                coordinated.append({
                    'gene_id': row['gene_id'],
                    'mod_type': '4mC',
                    'methyl_change': row['4mC_change'],
                    'methyl_category': 'Lost',
                    'log2FC': row['log2FoldChange'],
                    'padj': row['padj'],
                    'coordination': 'Lost_Down',
                    'correlation_type': 'Positive',
                })

    return pd.DataFrame(coordinated)

def compare_timepoints(t2_coord, t3vst1_coord, t3vst2_coord):
    """Compare coordinated genes across timepoints"""
    t2_genes = set(t2_coord['gene_id'].unique())
    t3vst1_genes = set(t3vst1_coord['gene_id'].unique()) if len(t3vst1_coord) > 0 else set()
    t3vst2_genes = set(t3vst2_coord['gene_id'].unique()) if len(t3vst2_coord) > 0 else set()

    comparison = {
        't2vst1_total': len(t2_genes),
        't3vst1_total': len(t3vst1_genes),
        't3vst2_total': len(t3vst2_genes),
        't2_and_t3vst1': len(t2_genes & t3vst1_genes),
        't2_only': len(t2_genes - t3vst1_genes),
        't3vst1_only': len(t3vst1_genes - t2_genes),
        'all_three': len(t2_genes & t3vst1_genes & t3vst2_genes),
    }

    return comparison

def create_visualizations(t2_coord, t3vst1_coord, t3vst2_coord, comparison):
    """Create comparison visualizations"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # 1. Coordinated gene counts
    ax1 = axes[0, 0]
    timepoints = ['T2 vs T1', 'T3 vs T1', 'T3 vs T2']
    counts = [comparison['t2vst1_total'], comparison['t3vst1_total'], comparison['t3vst2_total']]
    colors = ['#4e79a7', '#f28e2b', '#59a14f']
    ax1.bar(timepoints, counts, color=colors)
    ax1.set_ylabel('Number of coordinated genes')
    ax1.set_title('Coordinated Genes by Comparison')
    for i, v in enumerate(counts):
        ax1.text(i, v + 5, str(v), ha='center')

    # 2. Correlation type distribution
    ax2 = axes[0, 1]

    data_for_plot = []
    for name, df in [('T2vsT1', t2_coord), ('T3vsT1', t3vst1_coord), ('T3vsT2', t3vst2_coord)]:
        if len(df) > 0:
            # Handle different column names
            if 'correlation_type' in df.columns:
                counts = df['correlation_type'].value_counts()
            elif 'coordination' in df.columns:
                # Derive correlation type from coordination
                df_copy = df.copy()
                df_copy['correlation_type'] = df_copy['coordination'].apply(
                    lambda x: 'Positive' if x in ['Gained_Up', 'Lost_Down', 'Increased_Up', 'Decreased_Down']
                    else 'Negative'
                )
                counts = df_copy['correlation_type'].value_counts()
            else:
                continue
            for corr_type in ['Positive', 'Negative']:
                data_for_plot.append({
                    'Comparison': name,
                    'Type': corr_type,
                    'Count': counts.get(corr_type, 0)
                })

    if data_for_plot:
        plot_df = pd.DataFrame(data_for_plot)
        plot_pivot = plot_df.pivot(index='Comparison', columns='Type', values='Count').fillna(0)
        plot_pivot.plot(kind='bar', ax=ax2, color=['#59a14f', '#e15759'])
        ax2.set_ylabel('Number of genes')
        ax2.set_title('Correlation Type Distribution')
        ax2.legend(title='Correlation')
        ax2.set_xticklabels(ax2.get_xticklabels(), rotation=0)

    # 3. Venn-like overlap
    ax3 = axes[1, 0]

    # Simple bar for overlap
    overlap_data = {
        'T2 only': comparison['t2_only'],
        'Shared T2-T3': comparison['t2_and_t3vst1'],
        'T3vsT1 only': comparison['t3vst1_only'],
    }
    ax3.bar(overlap_data.keys(), overlap_data.values(), color=['#4e79a7', '#9467bd', '#f28e2b'])
    ax3.set_ylabel('Number of genes')
    ax3.set_title('Overlap: T2vsT1 and T3vsT1 Coordinated Genes')

    # 4. Methylation type breakdown
    ax4 = axes[1, 1]

    methyl_data = []
    for name, df in [('T2vsT1', t2_coord), ('T3vsT1', t3vst1_coord), ('T3vsT2', t3vst2_coord)]:
        if len(df) > 0 and 'mod_type' in df.columns:
            counts = df['mod_type'].value_counts()
            for mod in ['4mC', '6mA']:
                methyl_data.append({
                    'Comparison': name,
                    'Type': mod,
                    'Count': counts.get(mod, 0)
                })

    if methyl_data:
        methyl_df = pd.DataFrame(methyl_data)
        methyl_pivot = methyl_df.pivot(index='Comparison', columns='Type', values='Count').fillna(0)
        methyl_pivot.plot(kind='bar', ax=ax4, color=['#f28e2b', '#4e79a7'])
        ax4.set_ylabel('Number of genes')
        ax4.set_title('Methylation Type Distribution')
        ax4.legend(title='Modification')
        ax4.set_xticklabels(ax4.get_xticklabels(), rotation=0)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 't3_coordinated_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: t3_coordinated_comparison.png")

def main():
    print("="*60)
    print("T3 COORDINATED GENE ANALYSIS")
    print("="*60)

    # Load data
    print("\n1. Loading data...")
    methyl_expr, t2vst1, t3vst1, t3vst2, t2_coord = load_data()
    print(f"   Loaded {len(methyl_expr)} genes with methylation")
    print(f"   T2vsT1 coordinated: {len(t2_coord)} genes")

    # Identify T3vsT1 coordinated genes
    print("\n2. Identifying T3vsT1 coordinated genes...")
    t3vst1_coord = identify_coordinated_genes(methyl_expr, t3vst1, 'T3vsT1', 'T1', 'T3')
    print(f"   Found {len(t3vst1_coord)} T3vsT1 coordinated genes")
    t3vst1_coord.to_csv(OUTPUT_DIR / 'T3vsT1_coordinated_genes.csv', index=False)

    # Identify T3vsT2 coordinated genes
    print("\n3. Identifying T3vsT2 coordinated genes...")
    t3vst2_coord = identify_coordinated_genes(methyl_expr, t3vst2, 'T3vsT2', 'T2', 'T3')
    print(f"   Found {len(t3vst2_coord)} T3vsT2 coordinated genes")
    t3vst2_coord.to_csv(OUTPUT_DIR / 'T3vsT2_coordinated_genes.csv', index=False)

    # Compare timepoints
    print("\n4. Comparing timepoints...")
    comparison = compare_timepoints(t2_coord, t3vst1_coord, t3vst2_coord)

    print(f"   T2vsT1: {comparison['t2vst1_total']} genes")
    print(f"   T3vsT1: {comparison['t3vst1_total']} genes")
    print(f"   T3vsT2: {comparison['t3vst2_total']} genes")
    print(f"   Shared T2-T3vsT1: {comparison['t2_and_t3vst1']} genes")

    # Save comparison
    pd.DataFrame([comparison]).to_csv(OUTPUT_DIR / 'timepoint_comparison.csv', index=False)

    # Create visualizations
    print("\n5. Creating visualizations...")
    create_visualizations(t2_coord, t3vst1_coord, t3vst2_coord, comparison)

    # Analyze patterns
    print("\n6. Analyzing patterns...")

    if len(t3vst1_coord) > 0:
        # Correlation type
        pos_t3 = (t3vst1_coord['correlation_type'] == 'Positive').sum()
        neg_t3 = (t3vst1_coord['correlation_type'] == 'Negative').sum()
        print(f"   T3vsT1 correlation: {pos_t3} positive, {neg_t3} negative")

        # Methylation type
        mC4_t3 = (t3vst1_coord['mod_type'] == '4mC').sum()
        mA6_t3 = (t3vst1_coord['mod_type'] == '6mA').sum()
        print(f"   T3vsT1 methylation: {mC4_t3} 4mC, {mA6_t3} 6mA")

    # Summary
    print("\n" + "="*60)
    print("KEY FINDINGS")
    print("="*60)
    print(f"\n1. T3vsT1 coordinated genes: {comparison['t3vst1_total']}")
    print(f"   (vs T2vsT1: {comparison['t2vst1_total']})")

    print(f"\n2. Overlap with T2vsT1: {comparison['t2_and_t3vst1']} genes")
    print(f"   - Persistent coordination across timepoints")

    print(f"\n3. T3-specific: {comparison['t3vst1_only']} genes")
    print(f"   - New methylation-expression coordination in late phase")

    print(f"\n4. T3vsT2 (transition): {comparison['t3vst2_total']} genes")
    print(f"   - Changes between T2 and T3")

    print(f"\nOutput directory: {OUTPUT_DIR}")

if __name__ == '__main__':
    main()
