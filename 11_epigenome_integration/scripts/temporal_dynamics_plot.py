#!/usr/bin/env python3
"""
Temporal Dynamics Plot
Creates line charts showing methylation frequency and expression changes across timepoints (T1→T2→T3)
Reference: Zhang et al., 2023 (Salmonella) - Fig showing ahpCF promoter methylation and expression

Streptomyces coelicolor A3(2) M145 Project
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os

# Configuration
INTEGRATION_DIR = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis"
OUTPUT_DIR = os.path.join(INTEGRATION_DIR, "temporal_dynamics")
METHYL_PATH = os.path.join(INTEGRATION_DIR, "integrated_methyl_expression_weighted.csv")
DESEQ2_DIR = "/Users/okaban/bioinfo/rna-seq/04_deseq2/analysis/04_deseq2_260128_v1/results"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Key genes for temporal analysis
KEY_GENES = {
    'RamR': {
        'gene_id': 'SC_RS35610',
        'old_locus': 'SCO6685',
        'function': 'SapB pathway activator',
        'pattern': 'Demethylation → Activation'
    },
    'NsdB': {
        'gene_id': 'SC_RS38475',
        'old_locus': 'SCO7252',
        'function': 'Master developmental regulator',
        'pattern': 'Methylation gain → Strong activation'
    },
    'Red_oxygenase': {
        'gene_id': 'SC_RS31730',
        'old_locus': 'SCO5897',
        'function': 'Red cluster biosynthesis',
        'pattern': '4mC gain → Upregulation'
    },
    'Act_NmrA': {
        'gene_id': 'SC_RS27555',
        'old_locus': 'SCO5079',
        'function': 'Act cluster protein',
        'pattern': '4mC loss → Upregulation'
    },
    'Cpk_carboxylase': {
        'gene_id': 'SC_RS33670',
        'old_locus': 'SCO6284',
        'function': 'Cpk cluster enzyme',
        'pattern': '6mA gain → Strong upregulation'
    },
    'CDA_trpC': {
        'gene_id': 'SC_RS18170',
        'old_locus': 'SCO3211',
        'function': 'CDA cluster TrpC',
        'pattern': 'Stable methylation, strong upregulation'
    }
}

def load_data():
    """Load methylation and expression data"""
    # Load integrated data
    expr_methyl = pd.read_csv(METHYL_PATH)

    # Load normalized counts for expression levels
    counts_path = os.path.join(DESEQ2_DIR, "normalized_counts_M145.tsv")
    counts = pd.read_csv(counts_path, sep='\t')

    return expr_methyl, counts

def calculate_mean_expression(counts_df, gene_id):
    """Calculate mean normalized expression for each timepoint"""
    gene_row = counts_df[counts_df['gene_id'] == gene_id]
    if len(gene_row) == 0:
        return {'T1': np.nan, 'T2': np.nan, 'T3': np.nan}

    gene_row = gene_row.iloc[0]

    # Calculate means by timepoint
    t1_cols = [c for c in counts_df.columns if c.startswith('M145_1_')]
    t2_cols = [c for c in counts_df.columns if c.startswith('M145_2_')]
    t3_cols = [c for c in counts_df.columns if c.startswith('M145_3_')]

    return {
        'T1': gene_row[t1_cols].mean() if t1_cols else np.nan,
        'T2': gene_row[t2_cols].mean() if t2_cols else np.nan,
        'T3': gene_row[t3_cols].mean() if t3_cols else np.nan
    }

def create_temporal_plot_single(gene_name, gene_info, expr_methyl, counts_df, output_dir):
    """Create temporal dynamics plot for a single gene"""

    gene_id = gene_info['gene_id']
    gene_data = expr_methyl[expr_methyl['gene_id'] == gene_id]

    if len(gene_data) == 0:
        print(f"    No data for {gene_name}")
        return None

    row = gene_data.iloc[0]

    # Get methylation data
    timepoints = ['T1', 'T2', 'T3']
    x_pos = [1, 2, 3]

    # 6mA methylation frequency (mean across sites in promoter)
    methyl_6mA = [
        row['6mA_T1_mean_freq'],
        row['6mA_T2_mean_freq'],
        row['6mA_T3_mean_freq']
    ]

    # 4mC methylation frequency
    methyl_4mC = [
        row['4mC_T1_mean_freq'],
        row['4mC_T2_mean_freq'],
        row['4mC_T3_mean_freq']
    ]

    # Get expression levels
    expr_levels = calculate_mean_expression(counts_df, gene_id)
    expr_values = [expr_levels['T1'], expr_levels['T2'], expr_levels['T3']]

    # Create figure
    fig = plt.figure(figsize=(10, 8))
    gs = GridSpec(2, 1, height_ratios=[1, 1], hspace=0.3)

    colors = {
        '6mA': '#E74C3C',
        '4mC': '#3498DB',
        'expr': '#27AE60'
    }

    # Top panel: Methylation frequency
    ax1 = fig.add_subplot(gs[0])

    # Plot 6mA if present
    if any(v > 0 for v in methyl_6mA):
        ax1.plot(x_pos, methyl_6mA, 'o-', color=colors['6mA'], linewidth=2.5,
                markersize=10, label='6mA', markerfacecolor='white', markeredgewidth=2)
        for x, y in zip(x_pos, methyl_6mA):
            if y > 0:
                ax1.annotate(f'{y:.1f}%', (x, y), textcoords="offset points",
                            xytext=(0, 10), ha='center', fontsize=9, color=colors['6mA'])

    # Plot 4mC if present
    if any(v > 0 for v in methyl_4mC):
        ax1.plot(x_pos, methyl_4mC, 's-', color=colors['4mC'], linewidth=2.5,
                markersize=10, label='4mC', markerfacecolor='white', markeredgewidth=2)
        for x, y in zip(x_pos, methyl_4mC):
            if y > 0:
                ax1.annotate(f'{y:.1f}%', (x, y), textcoords="offset points",
                            xytext=(0, -15), ha='center', fontsize=9, color=colors['4mC'])

    ax1.set_xlim(0.5, 3.5)
    ax1.set_ylim(0, 100)
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(timepoints, fontsize=12)
    ax1.set_ylabel('Methylation Frequency (%)', fontsize=12)
    ax1.set_title(f'Promoter Methylation Dynamics', fontsize=12, fontweight='bold')
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.3)

    # Add shading for methylation change phases
    if any(v > 0 for v in methyl_6mA + methyl_4mC):
        ax1.axvspan(1, 2, alpha=0.1, color='orange', label='T1→T2')
        ax1.axvspan(2, 3, alpha=0.1, color='purple', label='T2→T3')

    # Bottom panel: Expression level
    ax2 = fig.add_subplot(gs[1])

    # Log transform expression for better visualization
    expr_log = [np.log2(v + 1) if not np.isnan(v) else 0 for v in expr_values]

    ax2.plot(x_pos, expr_log, 'D-', color=colors['expr'], linewidth=2.5,
            markersize=10, label='Expression', markerfacecolor='white', markeredgewidth=2)

    for x, y, raw in zip(x_pos, expr_log, expr_values):
        if not np.isnan(raw):
            ax2.annotate(f'{raw:.0f}', (x, y), textcoords="offset points",
                        xytext=(0, 10), ha='center', fontsize=9, color=colors['expr'])

    ax2.set_xlim(0.5, 3.5)
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(timepoints, fontsize=12)
    ax2.set_ylabel('log2(Normalized Counts + 1)', fontsize=12)
    ax2.set_xlabel('Timepoint', fontsize=12)
    ax2.set_title('Gene Expression Dynamics', fontsize=12, fontweight='bold')
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3)

    # Add fold change annotations
    if not np.isnan(expr_values[0]) and not np.isnan(expr_values[1]) and expr_values[0] > 0:
        fc_t2t1 = expr_values[1] / expr_values[0]
        ax2.annotate(f'{fc_t2t1:.1f}×', xy=(1.5, (expr_log[0] + expr_log[1])/2),
                    fontsize=10, ha='center', color='#E67E22', fontweight='bold')

    if not np.isnan(expr_values[1]) and not np.isnan(expr_values[2]) and expr_values[1] > 0:
        fc_t3t2 = expr_values[2] / expr_values[1]
        ax2.annotate(f'{fc_t3t2:.1f}×', xy=(2.5, (expr_log[1] + expr_log[2])/2),
                    fontsize=10, ha='center', color='#8E44AD', fontweight='bold')

    # Main title
    fig.suptitle(f'{gene_name} ({gene_info["old_locus"]})\n{gene_info["function"]}\n'
                 f'Pattern: {gene_info["pattern"]}',
                 fontsize=14, fontweight='bold', y=0.98)

    plt.tight_layout()
    plt.subplots_adjust(top=0.85)

    # Save
    output_path = os.path.join(output_dir, f'temporal_dynamics_{gene_name}.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()

    return output_path

def create_combined_temporal_plot(expr_methyl, counts_df, output_dir):
    """Create combined temporal dynamics plot for all key genes"""

    fig, axes = plt.subplots(len(KEY_GENES), 2, figsize=(14, 4*len(KEY_GENES)))

    timepoints = ['T1', 'T2', 'T3']
    x_pos = [1, 2, 3]

    colors = {
        '6mA': '#E74C3C',
        '4mC': '#3498DB',
        'expr': '#27AE60'
    }

    for idx, (gene_name, gene_info) in enumerate(KEY_GENES.items()):
        gene_id = gene_info['gene_id']
        gene_data = expr_methyl[expr_methyl['gene_id'] == gene_id]

        if len(gene_data) == 0:
            continue

        row = gene_data.iloc[0]

        # Get methylation data
        methyl_6mA = [row['6mA_T1_mean_freq'], row['6mA_T2_mean_freq'], row['6mA_T3_mean_freq']]
        methyl_4mC = [row['4mC_T1_mean_freq'], row['4mC_T2_mean_freq'], row['4mC_T3_mean_freq']]

        # Get expression levels
        expr_levels = calculate_mean_expression(counts_df, gene_id)
        expr_values = [expr_levels['T1'], expr_levels['T2'], expr_levels['T3']]
        expr_log = [np.log2(v + 1) if not np.isnan(v) else 0 for v in expr_values]

        # Left panel: Methylation
        ax_left = axes[idx, 0]

        if any(v > 0 for v in methyl_6mA):
            ax_left.plot(x_pos, methyl_6mA, 'o-', color=colors['6mA'], linewidth=2,
                        markersize=8, label='6mA')

        if any(v > 0 for v in methyl_4mC):
            ax_left.plot(x_pos, methyl_4mC, 's-', color=colors['4mC'], linewidth=2,
                        markersize=8, label='4mC')

        ax_left.set_xlim(0.5, 3.5)
        ax_left.set_ylim(0, 100)
        ax_left.set_xticks(x_pos)
        ax_left.set_xticklabels(timepoints)
        ax_left.set_ylabel('Methyl Freq (%)')
        ax_left.set_title(f'{gene_name} ({gene_info["old_locus"]})\n{gene_info["pattern"]}',
                         fontsize=10, fontweight='bold')
        ax_left.legend(loc='upper right', fontsize=8)
        ax_left.grid(True, alpha=0.3)

        # Right panel: Expression
        ax_right = axes[idx, 1]

        ax_right.plot(x_pos, expr_log, 'D-', color=colors['expr'], linewidth=2,
                     markersize=8, label='Expression')

        # Add actual values
        for x, y, raw in zip(x_pos, expr_log, expr_values):
            if not np.isnan(raw):
                ax_right.annotate(f'{raw:.0f}', (x, y), textcoords="offset points",
                                 xytext=(0, 8), ha='center', fontsize=8)

        ax_right.set_xlim(0.5, 3.5)
        ax_right.set_xticks(x_pos)
        ax_right.set_xticklabels(timepoints)
        ax_right.set_ylabel('log2(Counts+1)')
        ax_right.set_title(f'Expression: {gene_info["function"]}', fontsize=10)
        ax_right.grid(True, alpha=0.3)

    # Add main title
    fig.suptitle('Temporal Dynamics of Methylation and Expression\n'
                 'Key Epigenetically Regulated Genes in S. coelicolor M145',
                 fontsize=14, fontweight='bold', y=1.01)

    plt.tight_layout()

    output_path = os.path.join(output_dir, 'temporal_dynamics_all_genes.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()

    return output_path

def create_correlation_trajectory_plot(expr_methyl, counts_df, output_dir):
    """Create a plot showing methylation vs expression trajectories"""

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    colors_tp = {'T1→T2': '#F39C12', 'T2→T3': '#9B59B6'}

    for idx, (gene_name, gene_info) in enumerate(list(KEY_GENES.items())[:6]):
        ax = axes[idx]
        gene_id = gene_info['gene_id']
        gene_data = expr_methyl[expr_methyl['gene_id'] == gene_id]

        if len(gene_data) == 0:
            continue

        row = gene_data.iloc[0]

        # Get methylation changes
        methyl_change_t2t1 = row['6mA_change_T2_vs_T1'] if row['6mA_T1_mean_freq'] > 0 or row['6mA_T2_mean_freq'] > 0 else row['4mC_change_T2_vs_T1']
        methyl_change_t3t2 = row['6mA_change_T3_vs_T2'] if row['6mA_T2_mean_freq'] > 0 or row['6mA_T3_mean_freq'] > 0 else row['4mC_change_T3_vs_T2']

        # Get expression changes
        log2fc_t2t1 = row['log2FC_T2_vs_T1'] if pd.notna(row['log2FC_T2_vs_T1']) else 0
        log2fc_t3t2 = row['log2FC_T3_vs_T2'] if pd.notna(row['log2FC_T3_vs_T2']) else 0

        # Plot trajectory
        ax.annotate('', xy=(methyl_change_t2t1, log2fc_t2t1), xytext=(0, 0),
                   arrowprops=dict(arrowstyle='->', color=colors_tp['T1→T2'], lw=2))
        ax.annotate('', xy=(methyl_change_t2t1 + methyl_change_t3t2, log2fc_t2t1 + log2fc_t3t2),
                   xytext=(methyl_change_t2t1, log2fc_t2t1),
                   arrowprops=dict(arrowstyle='->', color=colors_tp['T2→T3'], lw=2))

        # Mark points
        ax.scatter([0], [0], c='green', s=100, zorder=5, label='T1')
        ax.scatter([methyl_change_t2t1], [log2fc_t2t1], c='orange', s=100, zorder=5, label='T2')
        ax.scatter([methyl_change_t2t1 + methyl_change_t3t2], [log2fc_t2t1 + log2fc_t3t2],
                  c='purple', s=100, zorder=5, label='T3')

        ax.axhline(0, color='gray', linestyle='--', alpha=0.5)
        ax.axvline(0, color='gray', linestyle='--', alpha=0.5)

        ax.set_xlabel('ΔMethylation (%)')
        ax.set_ylabel('log2FC Expression')
        ax.set_title(f'{gene_name}\n({gene_info["old_locus"]})', fontweight='bold')
        ax.legend(loc='best', fontsize=8)
        ax.grid(True, alpha=0.3)

    fig.suptitle('Methylation-Expression Trajectories Across Timepoints\n'
                 'T1 (green) → T2 (orange) → T3 (purple)',
                 fontsize=14, fontweight='bold', y=1.02)

    plt.tight_layout()

    output_path = os.path.join(output_dir, 'methylation_expression_trajectories.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()

    return output_path

def main():
    print("=" * 60)
    print("Temporal Dynamics Plot Generator")
    print("=" * 60)

    # Load data
    print("\nLoading data...")
    expr_methyl, counts_df = load_data()

    print(f"Loaded expression-methylation data: {len(expr_methyl)} genes")
    print(f"Loaded normalized counts: {len(counts_df)} genes")

    # Generate individual plots
    print("\nGenerating individual temporal dynamics plots...")
    for gene_name, gene_info in KEY_GENES.items():
        print(f"  Processing {gene_name}...")
        output_path = create_temporal_plot_single(gene_name, gene_info, expr_methyl, counts_df, OUTPUT_DIR)
        if output_path:
            print(f"    Saved: {output_path}")

    # Generate combined plot
    print("\nGenerating combined temporal dynamics plot...")
    output_path = create_combined_temporal_plot(expr_methyl, counts_df, OUTPUT_DIR)
    print(f"  Saved: {output_path}")

    # Generate trajectory plot
    print("\nGenerating methylation-expression trajectory plot...")
    output_path = create_correlation_trajectory_plot(expr_methyl, counts_df, OUTPUT_DIR)
    print(f"  Saved: {output_path}")

    print("\n" + "=" * 60)
    print("Temporal Dynamics Generation Complete")
    print(f"Output directory: {OUTPUT_DIR}")
    print("=" * 60)

if __name__ == "__main__":
    main()
