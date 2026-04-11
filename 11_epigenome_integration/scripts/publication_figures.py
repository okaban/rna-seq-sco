#!/usr/bin/env python3
"""
Publication-quality figures for Methylation-Expression Correlation Analysis
Streptomyces coelicolor A3(2) M145: T2 vs T1 comparison
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from scipy import stats
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set publication-quality defaults
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 1.0,
    'xtick.major.width': 1.0,
    'ytick.major.width': 1.0,
})

# Paths
DATA_PATH = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis")
OUTPUT_PATH = DATA_PATH / "figures"
OUTPUT_PATH.mkdir(exist_ok=True)

# Color palette (colorblind-friendly)
COLORS = {
    '6mA': '#E64B35',  # Red
    '4mC': '#4DBBD5',  # Blue
    'up': '#00A087',   # Green
    'down': '#3C5488', # Dark blue
    'ns': '#B8B8B8',   # Gray
    'positive': '#E64B35',
    'negative': '#4DBBD5',
}

def load_data():
    """Load analysis results."""
    coord_df = pd.read_csv(DATA_PATH / 'T2vsT1_coordinated_genes.csv')
    pos_df = pd.read_csv(DATA_PATH / 'T2vsT1_positive_correlation.csv')
    neg_df = pd.read_csv(DATA_PATH / 'T2vsT1_negative_correlation.csv')
    return coord_df, pos_df, neg_df

def figure1_scatter_correlation(coord_df):
    """
    Figure 1: Scatter plot of methylation change vs expression change
    Main correlation figure for the paper
    """
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

    for idx, mod_type in enumerate(['6mA', '4mC']):
        ax = axes[idx]
        df = coord_df[coord_df['mod_type'] == mod_type].copy()

        # Calculate correlation
        valid = df.dropna(subset=['methyl_change', 'log2FC'])
        r, p = stats.spearmanr(valid['methyl_change'], valid['log2FC'])

        # Color by expression direction
        colors = []
        for _, row in df.iterrows():
            if row['expr_category'] == 'Up':
                colors.append(COLORS['up'])
            elif row['expr_category'] == 'Down':
                colors.append(COLORS['down'])
            else:
                colors.append(COLORS['ns'])

        # Scatter plot
        scatter = ax.scatter(df['methyl_change'], df['log2FC'],
                           c=colors, alpha=0.7, s=40, edgecolors='white', linewidth=0.5)

        # Add regression line
        z = np.polyfit(valid['methyl_change'], valid['log2FC'], 1)
        p_line = np.poly1d(z)
        x_line = np.linspace(df['methyl_change'].min(), df['methyl_change'].max(), 100)
        ax.plot(x_line, p_line(x_line), 'k--', linewidth=1.5, alpha=0.7)

        # Reference lines
        ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.5, alpha=0.5)
        ax.axvline(x=0, color='gray', linestyle='-', linewidth=0.5, alpha=0.5)

        # Labels
        ax.set_xlabel('Δ Methylation frequency (%)', fontweight='bold')
        ax.set_ylabel('log₂ Fold Change (T2/T1)', fontweight='bold')
        ax.set_title(f'{mod_type} Promoter Methylation', fontweight='bold', fontsize=12)

        # Statistics annotation
        stats_text = f'ρ = {r:.3f}\np = {p:.2e}\nn = {len(valid)}'
        ax.text(0.05, 0.95, stats_text, transform=ax.transAxes,
               fontsize=9, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))

        # Set axis limits with padding
        xlim = max(abs(df['methyl_change'].min()), abs(df['methyl_change'].max())) * 1.1
        ax.set_xlim(-xlim, xlim)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # Legend
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS['up'],
               markersize=8, label='Upregulated'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS['down'],
               markersize=8, label='Downregulated'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS['ns'],
               markersize=8, label='Mild change'),
    ]
    fig.legend(handles=legend_elements, loc='upper center', ncol=3,
              bbox_to_anchor=(0.5, 0.02), frameon=False)

    plt.tight_layout(rect=[0, 0.05, 1, 1])
    fig.savefig(OUTPUT_PATH / 'Fig1_methylation_expression_correlation.pdf')
    fig.savefig(OUTPUT_PATH / 'Fig1_methylation_expression_correlation.png', dpi=300)
    plt.close()
    print("  Saved: Fig1_methylation_expression_correlation")

def figure2_coordination_patterns(coord_df):
    """
    Figure 2: Bar chart of coordination patterns
    """
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    for idx, mod_type in enumerate(['6mA', '4mC']):
        ax = axes[idx]
        df = coord_df[coord_df['mod_type'] == mod_type]

        # Define categories
        categories = {
            'Gained + Up': len(df[(df['methyl_category'] == 'Gained') & (df['expr_category'] == 'Up')]),
            'Gained + Down': len(df[(df['methyl_category'] == 'Gained') & (df['expr_category'] == 'Down')]),
            'Lost + Up': len(df[(df['methyl_category'] == 'Lost') & (df['expr_category'] == 'Up')]),
            'Lost + Down': len(df[(df['methyl_category'] == 'Lost') & (df['expr_category'] == 'Down')]),
            'Increased + Up': len(df[(df['methyl_category'] == 'Increased') & (df['expr_category'] == 'Up')]),
            'Increased + Down': len(df[(df['methyl_category'] == 'Increased') & (df['expr_category'] == 'Down')]),
            'Decreased + Up': len(df[(df['methyl_category'] == 'Decreased') & (df['expr_category'] == 'Up')]),
            'Decreased + Down': len(df[(df['methyl_category'] == 'Decreased') & (df['expr_category'] == 'Down')]),
        }

        # Filter non-zero
        categories = {k: v for k, v in categories.items() if v > 0}

        # Colors: positive correlation = red tones, negative = blue tones
        bar_colors = []
        for cat in categories.keys():
            if ('Gained' in cat and 'Up' in cat) or ('Lost' in cat and 'Down' in cat) or \
               ('Increased' in cat and 'Up' in cat) or ('Decreased' in cat and 'Down' in cat):
                bar_colors.append(COLORS['positive'])
            else:
                bar_colors.append(COLORS['negative'])

        bars = ax.barh(list(categories.keys()), list(categories.values()),
                      color=bar_colors, edgecolor='white', linewidth=0.5)

        # Add value labels
        for bar, val in zip(bars, categories.values()):
            ax.text(val + 0.5, bar.get_y() + bar.get_height()/2,
                   str(val), va='center', fontsize=9)

        ax.set_xlabel('Number of genes', fontweight='bold')
        ax.set_title(f'{mod_type}', fontweight='bold', fontsize=12)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor=COLORS['positive'], label='Positive correlation'),
        mpatches.Patch(facecolor=COLORS['negative'], label='Negative correlation'),
    ]
    fig.legend(handles=legend_elements, loc='upper center', ncol=2,
              bbox_to_anchor=(0.5, 0.02), frameon=False)

    plt.tight_layout(rect=[0, 0.05, 1, 1])
    fig.savefig(OUTPUT_PATH / 'Fig2_coordination_patterns.pdf')
    fig.savefig(OUTPUT_PATH / 'Fig2_coordination_patterns.png', dpi=300)
    plt.close()
    print("  Saved: Fig2_coordination_patterns")

def figure3_volcano_style(coord_df):
    """
    Figure 3: Volcano-style plot showing methylation gain/loss with expression
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    # Combine both modification types
    df = coord_df.copy()

    # Calculate -log10(padj)
    df['neg_log_padj'] = -np.log10(df['padj'].clip(lower=1e-300))

    # Color scheme
    colors = []
    for _, row in df.iterrows():
        if row['methyl_category'] in ['Gained', 'Increased']:
            if row['log2FC'] > 0:
                colors.append('#E64B35')  # Red: methyl up, expr up
            else:
                colors.append('#F39B7F')  # Light red: methyl up, expr down
        elif row['methyl_category'] in ['Lost', 'Decreased']:
            if row['log2FC'] < 0:
                colors.append('#4DBBD5')  # Blue: methyl down, expr down
            else:
                colors.append('#91D1C2')  # Light blue: methyl down, expr up
        else:
            colors.append('#B8B8B8')  # Gray: stable methylation

    # Size by methylation change magnitude
    sizes = np.abs(df['methyl_change']) * 0.8 + 20

    scatter = ax.scatter(df['log2FC'], df['neg_log_padj'],
                        c=colors, s=sizes, alpha=0.7,
                        edgecolors='white', linewidth=0.3)

    # Reference lines
    ax.axhline(y=-np.log10(0.05), color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    ax.axvline(x=1, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    ax.axvline(x=-1, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    ax.axvline(x=0, color='gray', linestyle='-', linewidth=0.5, alpha=0.5)

    # Labels
    ax.set_xlabel('log₂ Fold Change (T2/T1)', fontweight='bold', fontsize=11)
    ax.set_ylabel('-log₁₀(adjusted p-value)', fontweight='bold', fontsize=11)
    ax.set_title('Methylation-Associated Differential Expression\n(T2 vs T1)',
                fontweight='bold', fontsize=12)

    # Annotate significance threshold
    ax.text(ax.get_xlim()[1] * 0.95, -np.log10(0.05) + 2,
           'padj = 0.05', fontsize=8, ha='right', color='gray')

    # Custom legend
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#E64B35',
               markersize=10, label='Methyl ↑, Expr ↑'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#F39B7F',
               markersize=10, label='Methyl ↑, Expr ↓'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#91D1C2',
               markersize=10, label='Methyl ↓, Expr ↑'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#4DBBD5',
               markersize=10, label='Methyl ↓, Expr ↓'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#B8B8B8',
               markersize=10, label='Methyl stable'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', frameon=True,
             fancybox=True, framealpha=0.9)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    fig.savefig(OUTPUT_PATH / 'Fig3_volcano_methylation.pdf')
    fig.savefig(OUTPUT_PATH / 'Fig3_volcano_methylation.png', dpi=300)
    plt.close()
    print("  Saved: Fig3_volcano_methylation")

def figure4_summary_statistics(coord_df, pos_df, neg_df):
    """
    Figure 4: Summary statistics panel
    """
    fig = plt.figure(figsize=(10, 8))

    # Panel A: Pie chart of correlation directions
    ax1 = fig.add_subplot(2, 2, 1)

    positive = len(pos_df)
    negative = len(neg_df)
    other = len(coord_df) - positive - negative

    sizes = [positive, negative, other]
    labels = [f'Positive\n(n={positive})', f'Negative\n(n={negative})', f'Other\n(n={other})']
    colors_pie = [COLORS['positive'], COLORS['negative'], COLORS['ns']]
    explode = (0.02, 0.02, 0)

    ax1.pie(sizes, labels=labels, colors=colors_pie, explode=explode,
           autopct='%1.1f%%', startangle=90, pctdistance=0.6,
           textprops={'fontsize': 9})
    ax1.set_title('A. Correlation Direction', fontweight='bold', fontsize=11)

    # Panel B: Modification type breakdown
    ax2 = fig.add_subplot(2, 2, 2)

    mod_data = {
        '6mA': {
            'Positive': len(pos_df[pos_df['mod_type'] == '6mA']),
            'Negative': len(neg_df[neg_df['mod_type'] == '6mA']),
        },
        '4mC': {
            'Positive': len(pos_df[pos_df['mod_type'] == '4mC']),
            'Negative': len(neg_df[neg_df['mod_type'] == '4mC']),
        }
    }

    x = np.arange(2)
    width = 0.35

    bars1 = ax2.bar(x - width/2, [mod_data['6mA']['Positive'], mod_data['4mC']['Positive']],
                   width, label='Positive', color=COLORS['positive'])
    bars2 = ax2.bar(x + width/2, [mod_data['6mA']['Negative'], mod_data['4mC']['Negative']],
                   width, label='Negative', color=COLORS['negative'])

    ax2.set_ylabel('Number of genes', fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(['6mA', '4mC'])
    ax2.legend(frameon=False)
    ax2.set_title('B. By Modification Type', fontweight='bold', fontsize=11)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # Add value labels
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax2.annotate(f'{int(height)}',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3), textcoords="offset points",
                        ha='center', va='bottom', fontsize=9)

    # Panel C: Expression magnitude distribution
    ax3 = fig.add_subplot(2, 2, 3)

    pos_fc = pos_df['log2FC'].abs()
    neg_fc = neg_df['log2FC'].abs()

    bp = ax3.boxplot([pos_fc, neg_fc], labels=['Positive\ncorrelation', 'Negative\ncorrelation'],
                    patch_artist=True, widths=0.6)

    bp['boxes'][0].set_facecolor(COLORS['positive'])
    bp['boxes'][1].set_facecolor(COLORS['negative'])
    for box in bp['boxes']:
        box.set_alpha(0.7)

    # Statistical test
    stat, pval = stats.mannwhitneyu(pos_fc, neg_fc, alternative='two-sided')
    ax3.text(0.5, 0.95, f'Mann-Whitney U\np = {pval:.3f}', transform=ax3.transAxes,
            fontsize=9, ha='center', va='top')

    ax3.set_ylabel('|log₂ Fold Change|', fontweight='bold')
    ax3.set_title('C. Expression Magnitude', fontweight='bold', fontsize=11)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)

    # Panel D: Methylation change magnitude
    ax4 = fig.add_subplot(2, 2, 4)

    pos_methyl = pos_df['methyl_change'].abs()
    neg_methyl = neg_df['methyl_change'].abs()

    bp2 = ax4.boxplot([pos_methyl, neg_methyl], labels=['Positive\ncorrelation', 'Negative\ncorrelation'],
                     patch_artist=True, widths=0.6)

    bp2['boxes'][0].set_facecolor(COLORS['positive'])
    bp2['boxes'][1].set_facecolor(COLORS['negative'])
    for box in bp2['boxes']:
        box.set_alpha(0.7)

    stat2, pval2 = stats.mannwhitneyu(pos_methyl, neg_methyl, alternative='two-sided')
    ax4.text(0.5, 0.95, f'Mann-Whitney U\np = {pval2:.3f}', transform=ax4.transAxes,
            fontsize=9, ha='center', va='top')

    ax4.set_ylabel('|Δ Methylation frequency| (%)', fontweight='bold')
    ax4.set_title('D. Methylation Magnitude', fontweight='bold', fontsize=11)
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)

    plt.tight_layout()
    fig.savefig(OUTPUT_PATH / 'Fig4_summary_statistics.pdf')
    fig.savefig(OUTPUT_PATH / 'Fig4_summary_statistics.png', dpi=300)
    plt.close()
    print("  Saved: Fig4_summary_statistics")

def figure5_top_genes_heatmap(pos_df, neg_df):
    """
    Figure 5: Heatmap of top coordinated genes
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 8))

    for idx, (df, title, cmap) in enumerate([
        (pos_df.head(20), 'Positive Correlation (Top 20)', 'Reds'),
        (neg_df.head(20), 'Negative Correlation (Top 20)', 'Blues')
    ]):
        ax = axes[idx]

        if len(df) == 0:
            continue

        # Prepare data
        df = df.copy()
        df['label'] = df.apply(
            lambda x: f"{x['gene_id']} ({x['old_locus_tag']})" if pd.notna(x['old_locus_tag']) and x['old_locus_tag'] else x['gene_id'],
            axis=1
        )

        # Create matrix
        data = df[['methyl_change', 'log2FC']].values
        labels = df['label'].values

        # Normalize for visualization
        norm_data = np.zeros_like(data)
        norm_data[:, 0] = data[:, 0] / 100  # Methylation to 0-1 scale
        norm_data[:, 1] = data[:, 1] / data[:, 1].max()  # Expression to 0-1 scale

        im = ax.imshow(np.abs(norm_data), cmap=cmap, aspect='auto', vmin=0, vmax=1)

        # Labels
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['Δ Methylation', 'log₂FC'], fontsize=10)
        ax.set_title(title, fontweight='bold', fontsize=11)

        # Add text annotations
        for i in range(len(df)):
            ax.text(0, i, f'{df.iloc[i]["methyl_change"]:+.0f}%',
                   ha='center', va='center', fontsize=7, color='white' if abs(norm_data[i, 0]) > 0.5 else 'black')
            ax.text(1, i, f'{df.iloc[i]["log2FC"]:+.1f}',
                   ha='center', va='center', fontsize=7, color='white' if abs(norm_data[i, 1]) > 0.5 else 'black')

    plt.tight_layout()
    fig.savefig(OUTPUT_PATH / 'Fig5_top_genes_heatmap.pdf')
    fig.savefig(OUTPUT_PATH / 'Fig5_top_genes_heatmap.png', dpi=300)
    plt.close()
    print("  Saved: Fig5_top_genes_heatmap")

def main():
    print("=" * 70)
    print("Generating Publication-Quality Figures")
    print("=" * 70)

    # Load data
    print("\nLoading data...")
    coord_df, pos_df, neg_df = load_data()
    print(f"  Total coordinated genes: {len(coord_df)}")
    print(f"  Positive correlation: {len(pos_df)}")
    print(f"  Negative correlation: {len(neg_df)}")

    # Generate figures
    print("\nGenerating figures...")
    figure1_scatter_correlation(coord_df)
    figure2_coordination_patterns(coord_df)
    figure3_volcano_style(coord_df)
    figure4_summary_statistics(coord_df, pos_df, neg_df)
    figure5_top_genes_heatmap(pos_df, neg_df)

    print(f"\nAll figures saved to: {OUTPUT_PATH}")
    print("\nFigure descriptions:")
    print("  Fig1: Main scatter plot - methylation vs expression correlation")
    print("  Fig2: Coordination patterns breakdown by category")
    print("  Fig3: Volcano-style plot with methylation direction")
    print("  Fig4: Summary statistics panel (4 subpanels)")
    print("  Fig5: Top genes heatmap")

if __name__ == "__main__":
    main()
