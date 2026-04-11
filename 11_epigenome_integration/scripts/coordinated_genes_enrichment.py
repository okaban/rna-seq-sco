#!/usr/bin/env python3
"""
Coordinated Genes Functional Enrichment Analysis
- GO/KEGG enrichment for 558 coordinated methylation-expression genes
- Compare positive vs negative correlation genes
- Compare Gained vs Lost methylation patterns

Streptomyces coelicolor A3(2) M145 Project
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import os

# Configuration
INTEGRATION_DIR = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis"
ANNOTATION_DIR = "/Users/okaban/bioinfo/rna-seq/05_annotation/analysis"
SUPP_DIR = "/Users/okaban/bioinfo/rna-seq/12_supplementary_figures/analysis/12_supplementary_260202_v1/tables"
OUTPUT_DIR = os.path.join(INTEGRATION_DIR, "coordinated_enrichment")
COORDINATED_PATH = os.path.join(INTEGRATION_DIR, "T2vsT1_coordinated_genes.csv")

os.makedirs(OUTPUT_DIR, exist_ok=True)

def load_data():
    """Load coordinated genes and annotation data"""
    # Load coordinated genes
    coord_df = pd.read_csv(COORDINATED_PATH)
    print(f"Loaded {len(coord_df)} coordinated gene-modification pairs")

    # Load COG classification if available
    cog_path = os.path.join(SUPP_DIR, "gene_COG_classification.tsv")
    if os.path.exists(cog_path):
        cog_df = pd.read_csv(cog_path, sep='\t')
        print(f"Loaded COG classification for {len(cog_df)} genes")
    else:
        cog_df = None

    return coord_df, cog_df

def classify_coordination_patterns(coord_df):
    """Classify genes by their coordination pattern"""
    # Positive correlation: methyl and expr change in same direction
    # Negative correlation: methyl and expr change in opposite direction

    positive_patterns = ['Gained_Up', 'Lost_Down', 'Increased_Up', 'Decreased_Down']
    negative_patterns = ['Gained_Down', 'Lost_Up', 'Increased_Down', 'Decreased_Up']

    coord_df['correlation_type'] = 'Other'
    coord_df.loc[coord_df['coordination'].isin(positive_patterns), 'correlation_type'] = 'Positive'
    coord_df.loc[coord_df['coordination'].isin(negative_patterns), 'correlation_type'] = 'Negative'

    # Methylation direction
    coord_df['methyl_direction'] = 'Stable'
    coord_df.loc[coord_df['coordination'].str.startswith('Gained'), 'methyl_direction'] = 'Gained'
    coord_df.loc[coord_df['coordination'].str.startswith('Lost'), 'methyl_direction'] = 'Lost'
    coord_df.loc[coord_df['coordination'].str.startswith('Increased'), 'methyl_direction'] = 'Increased'
    coord_df.loc[coord_df['coordination'].str.startswith('Decreased'), 'methyl_direction'] = 'Decreased'

    return coord_df

def cog_enrichment_analysis(coord_df, cog_df, group_col, group_name, output_dir):
    """Perform COG category enrichment analysis"""
    if cog_df is None:
        print("COG data not available, skipping COG enrichment")
        return None

    # Merge with COG data
    merged = coord_df.merge(cog_df, left_on='gene_id', right_on='gene_id', how='left')

    # Get COG categories for the group
    group_genes = merged[merged[group_col] == group_name]
    all_genes = merged

    results = []

    # COG category descriptions
    cog_descriptions = {
        'J': 'Translation, ribosomal structure',
        'A': 'RNA processing and modification',
        'K': 'Transcription',
        'L': 'Replication, recombination and repair',
        'B': 'Chromatin structure and dynamics',
        'D': 'Cell cycle control, cell division',
        'Y': 'Nuclear structure',
        'V': 'Defense mechanisms',
        'T': 'Signal transduction mechanisms',
        'M': 'Cell wall/membrane/envelope biogenesis',
        'N': 'Cell motility',
        'Z': 'Cytoskeleton',
        'W': 'Extracellular structures',
        'U': 'Intracellular trafficking and secretion',
        'O': 'PTM, protein turnover, chaperones',
        'C': 'Energy production and conversion',
        'G': 'Carbohydrate transport and metabolism',
        'E': 'Amino acid transport and metabolism',
        'F': 'Nucleotide transport and metabolism',
        'H': 'Coenzyme transport and metabolism',
        'I': 'Lipid transport and metabolism',
        'P': 'Inorganic ion transport and metabolism',
        'Q': 'Secondary metabolites biosynthesis',
        'R': 'General function prediction only',
        'S': 'Function unknown',
        '-': 'Not assigned'
    }

    # Count COG categories
    if 'COG_category' not in merged.columns:
        print("COG_category column not found")
        return None

    for cog_cat in cog_descriptions.keys():
        group_count = len(group_genes[group_genes['COG_category'].str.contains(cog_cat, na=False)])
        total_count = len(all_genes[all_genes['COG_category'].str.contains(cog_cat, na=False)])

        group_size = len(group_genes)
        total_size = len(all_genes)

        # Fisher's exact test
        # [[group_with_cog, group_without_cog], [other_with_cog, other_without_cog]]
        other_count = total_count - group_count
        table = [
            [group_count, group_size - group_count],
            [other_count, (total_size - group_size) - other_count]
        ]

        if min(table[0]) >= 0 and min(table[1]) >= 0:
            odds_ratio, p_value = stats.fisher_exact(table)
        else:
            odds_ratio, p_value = 1.0, 1.0

        group_pct = group_count / group_size * 100 if group_size > 0 else 0
        total_pct = total_count / total_size * 100 if total_size > 0 else 0

        results.append({
            'COG_category': cog_cat,
            'Description': cog_descriptions.get(cog_cat, 'Unknown'),
            'Group_count': group_count,
            'Group_total': group_size,
            'Group_pct': group_pct,
            'Background_count': total_count,
            'Background_total': total_size,
            'Background_pct': total_pct,
            'Fold_enrichment': group_pct / total_pct if total_pct > 0 else 0,
            'Odds_ratio': odds_ratio,
            'P_value': p_value
        })

    results_df = pd.DataFrame(results)
    results_df['P_adjusted'] = results_df['P_value'] * len(results_df)  # Bonferroni
    results_df['P_adjusted'] = results_df['P_adjusted'].clip(upper=1.0)

    # Sort by p-value
    results_df = results_df.sort_values('P_value')

    return results_df

def plot_cog_comparison(positive_cog, negative_cog, output_path):
    """Plot COG category comparison between positive and negative correlation genes"""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(16, 8))

    # Filter to categories with at least some genes
    pos_df = positive_cog[positive_cog['Group_count'] > 0].head(15)
    neg_df = negative_cog[negative_cog['Group_count'] > 0].head(15)

    # Positive correlation
    ax1 = axes[0]
    y_pos = range(len(pos_df))
    colors1 = ['#27AE60' if p < 0.05 else '#BDC3C7' for p in pos_df['P_value']]
    ax1.barh(y_pos, pos_df['Fold_enrichment'], color=colors1, alpha=0.8)
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels([f"{row['COG_category']}: {row['Description'][:30]}" for _, row in pos_df.iterrows()])
    ax1.axvline(1, color='red', linestyle='--', linewidth=1)
    ax1.set_xlabel('Fold Enrichment')
    ax1.set_title(f'Positive Correlation Genes (n={pos_df.iloc[0]["Group_total"]})\nMethyl↔Expr same direction')
    ax1.invert_yaxis()

    # Negative correlation
    ax2 = axes[1]
    y_pos = range(len(neg_df))
    colors2 = ['#E74C3C' if p < 0.05 else '#BDC3C7' for p in neg_df['P_value']]
    ax2.barh(y_pos, neg_df['Fold_enrichment'], color=colors2, alpha=0.8)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels([f"{row['COG_category']}: {row['Description'][:30]}" for _, row in neg_df.iterrows()])
    ax2.axvline(1, color='red', linestyle='--', linewidth=1)
    ax2.set_xlabel('Fold Enrichment')
    ax2.set_title(f'Negative Correlation Genes (n={neg_df.iloc[0]["Group_total"]})\nMethyl↔Expr opposite direction')
    ax2.invert_yaxis()

    plt.suptitle('COG Category Enrichment: Positive vs Negative Correlation Genes',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()

    print(f"Saved COG comparison plot: {output_path}")

def analyze_product_keywords(coord_df, group_col, output_dir):
    """Analyze gene product keywords"""
    from collections import Counter

    results = {}

    for group_name in coord_df[group_col].unique():
        group_df = coord_df[coord_df[group_col] == group_name]

        # Extract keywords from product descriptions
        keywords = []
        for product in group_df['product'].dropna():
            # Split and clean keywords
            words = product.lower().replace(',', ' ').replace('/', ' ').split()
            # Filter meaningful words
            meaningful = [w for w in words if len(w) > 3 and w not in
                         ['protein', 'family', 'domain', 'containing', 'like', 'type', 'subunit']]
            keywords.extend(meaningful)

        keyword_counts = Counter(keywords)
        results[group_name] = keyword_counts.most_common(20)

    return results

def create_coordination_summary_plot(coord_df, output_path):
    """Create summary plot of coordination patterns"""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # 1. Coordination pattern distribution
    ax1 = axes[0, 0]
    pattern_counts = coord_df['coordination'].value_counts()
    colors = plt.cm.Set3(np.linspace(0, 1, len(pattern_counts)))
    bars = ax1.bar(range(len(pattern_counts)), pattern_counts.values, color=colors)
    ax1.set_xticks(range(len(pattern_counts)))
    ax1.set_xticklabels(pattern_counts.index, rotation=45, ha='right')
    ax1.set_ylabel('Number of genes')
    ax1.set_title('Coordination Pattern Distribution')

    # Add value labels
    for bar, val in zip(bars, pattern_counts.values):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                str(val), ha='center', fontsize=9)

    # 2. Positive vs Negative correlation
    ax2 = axes[0, 1]
    corr_counts = coord_df['correlation_type'].value_counts()
    colors_pie = ['#27AE60', '#E74C3C', '#BDC3C7']
    wedges, texts, autotexts = ax2.pie(corr_counts.values, labels=corr_counts.index,
                                        autopct='%1.1f%%', colors=colors_pie[:len(corr_counts)])
    ax2.set_title('Correlation Type Distribution')

    # 3. Methylation direction
    ax3 = axes[1, 0]
    methyl_counts = coord_df['methyl_direction'].value_counts()
    colors_methyl = ['#3498DB', '#E74C3C', '#F39C12', '#9B59B6', '#95A5A6']
    ax3.bar(range(len(methyl_counts)), methyl_counts.values, color=colors_methyl[:len(methyl_counts)])
    ax3.set_xticks(range(len(methyl_counts)))
    ax3.set_xticklabels(methyl_counts.index)
    ax3.set_ylabel('Number of genes')
    ax3.set_title('Methylation Direction Distribution')

    # 4. Modification type
    ax4 = axes[1, 1]
    mod_counts = coord_df['mod_type'].value_counts()
    colors_mod = ['#E74C3C', '#3498DB']
    ax4.bar(mod_counts.index, mod_counts.values, color=colors_mod[:len(mod_counts)])
    ax4.set_ylabel('Number of genes')
    ax4.set_title('Modification Type Distribution')

    plt.suptitle('Coordinated Methylation-Expression Genes Summary\n'
                 f'Total: {len(coord_df)} gene-modification pairs',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()

    print(f"Saved coordination summary: {output_path}")

def create_enrichment_barplot(results_df, title, output_path, top_n=15):
    """Create enrichment bar plot"""
    import matplotlib.pyplot as plt

    # Filter significant and sort
    df = results_df.head(top_n)

    fig, ax = plt.subplots(figsize=(12, 8))

    y_pos = range(len(df))
    colors = ['#27AE60' if p < 0.05 else '#E74C3C' if p < 0.1 else '#BDC3C7'
              for p in df['P_value']]

    bars = ax.barh(y_pos, df['Fold_enrichment'], color=colors, alpha=0.8)

    ax.set_yticks(y_pos)
    ax.set_yticklabels([f"{row['COG_category']}: {row['Description']}" for _, row in df.iterrows()])
    ax.axvline(1, color='black', linestyle='--', linewidth=1, label='No enrichment')
    ax.set_xlabel('Fold Enrichment')
    ax.set_title(title)
    ax.invert_yaxis()

    # Add p-value annotations
    for i, (bar, pval) in enumerate(zip(bars, df['P_value'])):
        sig = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else ''
        ax.text(bar.get_width() + 0.05, bar.get_y() + bar.get_height()/2,
               f'{sig} p={pval:.3f}', va='center', fontsize=8)

    # Legend
    legend_elements = [
        plt.Rectangle((0,0), 1, 1, fc='#27AE60', label='p < 0.05'),
        plt.Rectangle((0,0), 1, 1, fc='#E74C3C', label='0.05 ≤ p < 0.1'),
        plt.Rectangle((0,0), 1, 1, fc='#BDC3C7', label='p ≥ 0.1')
    ]
    ax.legend(handles=legend_elements, loc='lower right')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()

    print(f"Saved enrichment plot: {output_path}")

def main():
    print("=" * 60)
    print("Coordinated Genes Functional Enrichment Analysis")
    print("=" * 60)

    # Load data
    coord_df, cog_df = load_data()

    # Classify coordination patterns
    coord_df = classify_coordination_patterns(coord_df)

    # Print summary
    print("\n--- Coordination Pattern Summary ---")
    print(f"Total gene-modification pairs: {len(coord_df)}")
    print(f"\nCorrelation types:")
    print(coord_df['correlation_type'].value_counts().to_string())
    print(f"\nMethylation directions:")
    print(coord_df['methyl_direction'].value_counts().to_string())

    # Save classified data
    coord_df.to_csv(os.path.join(OUTPUT_DIR, 'coordinated_genes_classified.csv'), index=False)

    # Create summary plot
    create_coordination_summary_plot(coord_df,
                                     os.path.join(OUTPUT_DIR, 'coordination_summary.png'))

    # COG enrichment analysis
    if cog_df is not None:
        print("\n--- COG Enrichment Analysis ---")

        # All coordinated genes vs background
        print("\nAnalyzing all coordinated genes...")
        all_cog = cog_enrichment_analysis(coord_df, cog_df, 'correlation_type', 'Positive', OUTPUT_DIR)
        # Actually analyze all by setting a dummy column
        coord_df['all'] = 'All'
        all_cog = cog_enrichment_analysis(coord_df, cog_df, 'all', 'All', OUTPUT_DIR)
        if all_cog is not None:
            all_cog.to_csv(os.path.join(OUTPUT_DIR, 'COG_enrichment_all_coordinated.csv'), index=False)
            create_enrichment_barplot(all_cog, 'COG Enrichment: All Coordinated Genes (n=558)',
                                     os.path.join(OUTPUT_DIR, 'COG_enrichment_all.png'))

        # Positive correlation genes
        print("\nAnalyzing positive correlation genes...")
        positive_cog = cog_enrichment_analysis(coord_df, cog_df, 'correlation_type', 'Positive', OUTPUT_DIR)
        if positive_cog is not None:
            positive_cog.to_csv(os.path.join(OUTPUT_DIR, 'COG_enrichment_positive.csv'), index=False)

        # Negative correlation genes
        print("\nAnalyzing negative correlation genes...")
        negative_cog = cog_enrichment_analysis(coord_df, cog_df, 'correlation_type', 'Negative', OUTPUT_DIR)
        if negative_cog is not None:
            negative_cog.to_csv(os.path.join(OUTPUT_DIR, 'COG_enrichment_negative.csv'), index=False)

        # Plot comparison
        if positive_cog is not None and negative_cog is not None:
            plot_cog_comparison(positive_cog, negative_cog,
                               os.path.join(OUTPUT_DIR, 'COG_comparison_pos_vs_neg.png'))

        # Gained vs Lost
        print("\nAnalyzing Gained vs Lost methylation...")
        gained_cog = cog_enrichment_analysis(coord_df, cog_df, 'methyl_direction', 'Gained', OUTPUT_DIR)
        lost_cog = cog_enrichment_analysis(coord_df, cog_df, 'methyl_direction', 'Lost', OUTPUT_DIR)

        if gained_cog is not None:
            gained_cog.to_csv(os.path.join(OUTPUT_DIR, 'COG_enrichment_gained.csv'), index=False)
        if lost_cog is not None:
            lost_cog.to_csv(os.path.join(OUTPUT_DIR, 'COG_enrichment_lost.csv'), index=False)

    # Analyze product keywords
    print("\n--- Product Keyword Analysis ---")
    keywords = analyze_product_keywords(coord_df, 'correlation_type', OUTPUT_DIR)

    with open(os.path.join(OUTPUT_DIR, 'product_keywords.txt'), 'w') as f:
        for group, kw_list in keywords.items():
            f.write(f"\n=== {group} correlation genes ===\n")
            for kw, count in kw_list:
                f.write(f"  {kw}: {count}\n")

    print("Top keywords by correlation type:")
    for group, kw_list in keywords.items():
        print(f"\n{group}: {[kw for kw, _ in kw_list[:5]]}")

    # Create summary report
    create_summary_report(coord_df, OUTPUT_DIR)

    print("\n" + "=" * 60)
    print("Enrichment Analysis Complete")
    print(f"Output directory: {OUTPUT_DIR}")
    print("=" * 60)

def create_summary_report(coord_df, output_dir):
    """Create summary report"""
    report_path = os.path.join(output_dir, 'ENRICHMENT_ANALYSIS_REPORT.md')

    with open(report_path, 'w') as f:
        f.write("# Coordinated Genes Functional Enrichment Report\n\n")
        f.write(f"**Date:** 2026-02-02\n")
        f.write(f"**Project:** *Streptomyces coelicolor* A3(2) M145\n\n")

        f.write("## Summary\n\n")
        f.write(f"Total coordinated gene-modification pairs: {len(coord_df)}\n\n")

        f.write("### Correlation Types\n\n")
        f.write("| Type | Count | Percentage |\n")
        f.write("|------|-------|------------|\n")
        for corr_type, count in coord_df['correlation_type'].value_counts().items():
            pct = count / len(coord_df) * 100
            f.write(f"| {corr_type} | {count} | {pct:.1f}% |\n")

        f.write("\n### Methylation Directions\n\n")
        f.write("| Direction | Count | Percentage |\n")
        f.write("|-----------|-------|------------|\n")
        for direction, count in coord_df['methyl_direction'].value_counts().items():
            pct = count / len(coord_df) * 100
            f.write(f"| {direction} | {count} | {pct:.1f}% |\n")

        f.write("\n## Output Files\n\n")
        f.write("| File | Description |\n")
        f.write("|------|-------------|\n")
        f.write("| `coordinated_genes_classified.csv` | Classified coordinated genes |\n")
        f.write("| `coordination_summary.png` | Summary visualization |\n")
        f.write("| `COG_enrichment_*.csv` | COG enrichment results |\n")
        f.write("| `COG_comparison_pos_vs_neg.png` | Positive vs negative comparison |\n")

    print(f"Saved report: {report_path}")

if __name__ == "__main__":
    main()
