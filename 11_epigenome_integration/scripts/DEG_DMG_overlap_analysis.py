#!/usr/bin/env python3
"""
DEGs vs DMGs Overlap Analysis
- Fisher's exact test for DEG-DMG overlap
- 2x2 contingency table
- Venn diagram visualization

Streptomyces coelicolor A3(2) M145 Project
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import os

# Paths
DESEQ2_DIR = "/Users/okaban/bioinfo/rna-seq/04_deseq2/analysis/04_deseq2_260128_v1/results"
INTEGRATION_DIR = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis"
OUTPUT_DIR = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/overlap_analysis"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Thresholds
DEG_PADJ_THRESHOLD = 0.05
DEG_LOG2FC_THRESHOLD = 1.0
DMG_CHANGE_THRESHOLD = 10.0  # 10% methylation change

def load_deseq2_results():
    """Load DESeq2 results for all comparisons"""
    comparisons = {
        'T2vsT1': 'DESeq2_M145_2_vs_1.tsv',
        'T3vsT1': 'DESeq2_M145_3_vs_1.tsv',
        'T3vsT2': 'DESeq2_M145_3_vs_2.tsv'
    }

    results = {}
    for comp, filename in comparisons.items():
        filepath = os.path.join(DESEQ2_DIR, filename)
        df = pd.read_csv(filepath, sep='\t')
        results[comp] = df

    return results

def load_methylation_data():
    """Load integrated methylation-expression data"""
    filepath = os.path.join(INTEGRATION_DIR, "integrated_methyl_expression_weighted.csv")
    return pd.read_csv(filepath)

def identify_DEGs(deseq2_df, padj_thresh=DEG_PADJ_THRESHOLD, log2fc_thresh=DEG_LOG2FC_THRESHOLD):
    """Identify differentially expressed genes"""
    degs = deseq2_df[
        (deseq2_df['padj'] < padj_thresh) &
        (deseq2_df['log2FoldChange'].abs() > log2fc_thresh)
    ]['gene_id'].tolist()

    degs_up = deseq2_df[
        (deseq2_df['padj'] < padj_thresh) &
        (deseq2_df['log2FoldChange'] > log2fc_thresh)
    ]['gene_id'].tolist()

    degs_down = deseq2_df[
        (deseq2_df['padj'] < padj_thresh) &
        (deseq2_df['log2FoldChange'] < -log2fc_thresh)
    ]['gene_id'].tolist()

    return set(degs), set(degs_up), set(degs_down)

def identify_DMGs(methyl_df, comparison, change_thresh=DMG_CHANGE_THRESHOLD):
    """
    Identify differentially methylated genes (DMGs)
    DMG: gene with methylation site count change or frequency change > threshold
    """
    col_map = {
        'T2vsT1': ('6mA_change_T2_vs_T1', '4mC_change_T2_vs_T1'),
        'T3vsT1': ('6mA_change_T3_vs_T1', '4mC_change_T3_vs_T1'),
        'T3vsT2': ('6mA_change_T3_vs_T2', '4mC_change_T3_vs_T2')
    }

    col_6mA, col_4mC = col_map[comparison]

    # DMG: absolute methylation change > threshold for either 6mA or 4mC
    dmgs = methyl_df[
        (methyl_df[col_6mA].abs() > change_thresh) |
        (methyl_df[col_4mC].abs() > change_thresh)
    ]['gene_id'].tolist()

    # Hypermethylated: methylation increased
    dmgs_hyper = methyl_df[
        (methyl_df[col_6mA] > change_thresh) |
        (methyl_df[col_4mC] > change_thresh)
    ]['gene_id'].tolist()

    # Hypomethylated: methylation decreased
    dmgs_hypo = methyl_df[
        (methyl_df[col_6mA] < -change_thresh) |
        (methyl_df[col_4mC] < -change_thresh)
    ]['gene_id'].tolist()

    return set(dmgs), set(dmgs_hyper), set(dmgs_hypo)

def create_contingency_table(all_genes, degs, dmgs):
    """Create 2x2 contingency table"""
    deg_and_dmg = len(degs & dmgs)
    deg_not_dmg = len(degs - dmgs)
    dmg_not_deg = len(dmgs - degs)
    neither = len(all_genes - degs - dmgs)

    table = np.array([
        [deg_and_dmg, deg_not_dmg],
        [dmg_not_deg, neither]
    ])

    return table

def fisher_exact_test(table):
    """Perform Fisher's exact test"""
    odds_ratio, p_value = stats.fisher_exact(table)
    return odds_ratio, p_value

def create_venn_diagram(degs, dmgs, comparison, output_path, title_suffix=""):
    """Create Venn diagram for DEGs vs DMGs overlap"""
    fig, ax = plt.subplots(figsize=(8, 6))

    v = venn2([degs, dmgs], set_labels=('DEGs', 'DMGs'), ax=ax)

    # Customize colors
    if v.get_patch_by_id('10'):
        v.get_patch_by_id('10').set_color('#FF6B6B')
        v.get_patch_by_id('10').set_alpha(0.7)
    if v.get_patch_by_id('01'):
        v.get_patch_by_id('01').set_color('#4ECDC4')
        v.get_patch_by_id('01').set_alpha(0.7)
    if v.get_patch_by_id('11'):
        v.get_patch_by_id('11').set_color('#95E1D3')
        v.get_patch_by_id('11').set_alpha(0.7)

    plt.title(f'DEGs vs DMGs Overlap ({comparison}){title_suffix}\n'
              f'DEGs: |log2FC|>1, padj<0.05 | DMGs: |Δmethyl|>10%',
              fontsize=12)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()

def main():
    print("=" * 60)
    print("DEGs vs DMGs Overlap Analysis")
    print("=" * 60)

    # Load data
    print("\nLoading data...")
    deseq2_results = load_deseq2_results()
    methyl_df = load_methylation_data()

    # Get all genes (background)
    all_genes = set(methyl_df['gene_id'].tolist())
    print(f"Total genes in analysis: {len(all_genes)}")

    # Results storage
    results = []

    # Analyze each comparison
    for comparison in ['T2vsT1', 'T3vsT1', 'T3vsT2']:
        print(f"\n{'='*60}")
        print(f"Analyzing {comparison}")
        print("=" * 60)

        # Get DEGs
        degs, degs_up, degs_down = identify_DEGs(deseq2_results[comparison])
        print(f"DEGs: {len(degs)} (Up: {len(degs_up)}, Down: {len(degs_down)})")

        # Get DMGs
        dmgs, dmgs_hyper, dmgs_hypo = identify_DMGs(methyl_df, comparison)
        print(f"DMGs: {len(dmgs)} (Hyper: {len(dmgs_hyper)}, Hypo: {len(dmgs_hypo)})")

        # Overall overlap
        overlap = degs & dmgs
        print(f"Overlap (DEG ∩ DMG): {len(overlap)}")

        # Create contingency table
        table = create_contingency_table(all_genes, degs, dmgs)
        print(f"\nContingency Table:")
        print(f"                    DMG      Not-DMG")
        print(f"  DEG              {table[0,0]:5d}      {table[0,1]:5d}")
        print(f"  Not-DEG          {table[1,0]:5d}      {table[1,1]:5d}")

        # Fisher's exact test
        odds_ratio, p_value = fisher_exact_test(table)
        print(f"\nFisher's Exact Test:")
        print(f"  Odds Ratio: {odds_ratio:.4f}")
        print(f"  P-value: {p_value:.4e}")

        # Interpretation
        if p_value < 0.05:
            if odds_ratio > 1:
                interpretation = "Significant positive association (DEGs enriched in DMGs)"
            else:
                interpretation = "Significant negative association (DEGs depleted in DMGs)"
        else:
            interpretation = "No significant association"
        print(f"  Interpretation: {interpretation}")

        # Store results
        results.append({
            'Comparison': comparison,
            'Total_genes': len(all_genes),
            'DEGs': len(degs),
            'DEGs_up': len(degs_up),
            'DEGs_down': len(degs_down),
            'DMGs': len(dmgs),
            'DMGs_hyper': len(dmgs_hyper),
            'DMGs_hypo': len(dmgs_hypo),
            'Overlap_DEG_DMG': len(overlap),
            'DEG_only': len(degs - dmgs),
            'DMG_only': len(dmgs - degs),
            'Neither': len(all_genes - degs - dmgs),
            'Odds_ratio': odds_ratio,
            'P_value': p_value,
            'Significant': p_value < 0.05,
            'Interpretation': interpretation
        })

        # Create Venn diagram
        venn_path = os.path.join(OUTPUT_DIR, f"venn_DEG_DMG_{comparison}.png")
        create_venn_diagram(degs, dmgs, comparison, venn_path)
        print(f"\nVenn diagram saved: {venn_path}")

        # Detailed subgroup analysis
        print(f"\n--- Subgroup Analysis ---")

        # Up-regulated DEGs vs Hypermethylated DMGs (positive correlation)
        up_hyper = degs_up & dmgs_hyper
        print(f"DEG↑ ∩ DMG↑ (hyper): {len(up_hyper)}")

        # Down-regulated DEGs vs Hypomethylated DMGs (positive correlation)
        down_hypo = degs_down & dmgs_hypo
        print(f"DEG↓ ∩ DMG↓ (hypo): {len(down_hypo)}")

        # Up-regulated DEGs vs Hypomethylated DMGs (negative correlation)
        up_hypo = degs_up & dmgs_hypo
        print(f"DEG↑ ∩ DMG↓ (hypo): {len(up_hypo)}")

        # Down-regulated DEGs vs Hypermethylated DMGs (negative correlation)
        down_hyper = degs_down & dmgs_hyper
        print(f"DEG↓ ∩ DMG↑ (hyper): {len(down_hyper)}")

        # Positive vs Negative correlation
        positive_corr = len(up_hyper) + len(down_hypo)
        negative_corr = len(up_hypo) + len(down_hyper)
        print(f"\nPositive correlation (methyl↔expr same direction): {positive_corr}")
        print(f"Negative correlation (methyl↔expr opposite direction): {negative_corr}")

    # Save results table
    results_df = pd.DataFrame(results)
    results_path = os.path.join(OUTPUT_DIR, "DEG_DMG_overlap_statistics.csv")
    results_df.to_csv(results_path, index=False)
    print(f"\n\nResults saved: {results_path}")

    # Create summary figure with all three comparisons
    create_summary_figure(results, OUTPUT_DIR)

    # Create contingency table visualization
    create_contingency_table_figure(results, OUTPUT_DIR)

    print("\n" + "=" * 60)
    print("Analysis Complete")
    print("=" * 60)

def create_summary_figure(results, output_dir):
    """Create summary bar chart of overlap statistics"""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    comparisons = [r['Comparison'] for r in results]

    # Left panel: DEGs and DMGs counts
    ax1 = axes[0]
    x = np.arange(len(comparisons))
    width = 0.35

    degs = [r['DEGs'] for r in results]
    dmgs = [r['DMGs'] for r in results]
    overlaps = [r['Overlap_DEG_DMG'] for r in results]

    bars1 = ax1.bar(x - width/2, degs, width, label='DEGs', color='#FF6B6B', alpha=0.8)
    bars2 = ax1.bar(x + width/2, dmgs, width, label='DMGs', color='#4ECDC4', alpha=0.8)

    ax1.set_xlabel('Comparison')
    ax1.set_ylabel('Number of Genes')
    ax1.set_title('DEGs and DMGs by Comparison')
    ax1.set_xticks(x)
    ax1.set_xticklabels(comparisons)
    ax1.legend()

    # Add overlap numbers as text
    for i, (d, m, o) in enumerate(zip(degs, dmgs, overlaps)):
        ax1.annotate(f'Overlap: {o}', xy=(i, max(d, m) + 50), ha='center', fontsize=10)

    # Right panel: Odds ratios with significance
    ax2 = axes[1]
    odds_ratios = [r['Odds_ratio'] for r in results]
    p_values = [r['P_value'] for r in results]

    colors = ['#27AE60' if p < 0.05 else '#E74C3C' for p in p_values]
    bars = ax2.bar(comparisons, odds_ratios, color=colors, alpha=0.8)

    ax2.axhline(y=1, color='black', linestyle='--', linewidth=1, label='No association (OR=1)')
    ax2.set_xlabel('Comparison')
    ax2.set_ylabel('Odds Ratio')
    ax2.set_title('Fisher\'s Exact Test: DEG-DMG Association')
    ax2.legend()

    # Add p-values as text
    for i, (bar, p) in enumerate(zip(bars, p_values)):
        sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
        ax2.annotate(f'p={p:.2e}\n{sig}',
                    xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                    ha='center', va='bottom', fontsize=9)

    plt.tight_layout()

    output_path = os.path.join(output_dir, "DEG_DMG_summary_statistics.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()
    print(f"Summary figure saved: {output_path}")

def create_contingency_table_figure(results, output_dir):
    """Create contingency table visualization"""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    for i, r in enumerate(results):
        ax = axes[i]

        table = np.array([
            [r['Overlap_DEG_DMG'], r['DEG_only']],
            [r['DMG_only'], r['Neither']]
        ])

        im = ax.imshow(table, cmap='Blues', aspect='auto')

        # Add text annotations
        for row in range(2):
            for col in range(2):
                ax.text(col, row, f'{table[row, col]}',
                       ha='center', va='center', fontsize=14, fontweight='bold',
                       color='white' if table[row, col] > table.max()/2 else 'black')

        ax.set_xticks([0, 1])
        ax.set_yticks([0, 1])
        ax.set_xticklabels(['DMG', 'Not-DMG'])
        ax.set_yticklabels(['DEG', 'Not-DEG'])
        ax.set_title(f"{r['Comparison']}\nOR={r['Odds_ratio']:.2f}, p={r['P_value']:.2e}")

        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    plt.suptitle('2×2 Contingency Tables: DEG vs DMG Overlap', fontsize=14, y=1.02)
    plt.tight_layout()

    output_path = os.path.join(output_dir, "contingency_tables.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()
    print(f"Contingency tables figure saved: {output_path}")

if __name__ == "__main__":
    main()
