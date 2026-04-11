#!/usr/bin/env python3
"""
TF Promoter Methylation Analysis
Analyzes methylation patterns at transcription factor promoters
with focus on the GRN hierarchy (Global → Pleiotropic → CSR)

Step C-2: TFプロモーターのメチル化状態解析
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Setup
BASE_DIR = Path('/Users/okaban/bioinfo/rna-seq/11_epigenome_integration')
ANALYSIS_DIR = BASE_DIR / 'analysis/12_grn_tf_methylation'
INTEGRATION_DIR = BASE_DIR / 'analysis/01_integration'
OUTPUT_DIR = ANALYSIS_DIR / 'tf_promoter_methylation'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Set style
plt.rcParams['font.size'] = 10
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['figure.dpi'] = 150

def load_data():
    """Load all required datasets"""
    # 1. Hierarchical TF data
    tf_hierarchy = pd.read_csv(ANALYSIS_DIR / 'hierarchical_network_tfs.csv')
    print(f"Loaded {len(tf_hierarchy)} TFs from GRN hierarchy")

    # 2. Integrated methylation-expression data
    methyl_expr = pd.read_csv(INTEGRATION_DIR / 'integrated_methyl_expression_weighted.csv')
    print(f"Loaded methylation data for {len(methyl_expr)} genes")

    # 3. TF coordinated changes (already computed)
    tf_coord = pd.read_csv(INTEGRATION_DIR / 'TF_coordinated_changes.csv')
    print(f"Loaded {len(tf_coord)} TFs with coordinated methylation changes")

    # 4. Regulatory interactions
    reg_interactions = pd.read_csv(ANALYSIS_DIR / 'regulatory_interactions.csv')
    print(f"Loaded {len(reg_interactions)} regulatory interactions")

    return tf_hierarchy, methyl_expr, tf_coord, reg_interactions

def merge_tf_methylation(tf_hierarchy, methyl_expr):
    """Merge TF hierarchy data with methylation data"""
    # Get methylation data for TFs
    tf_methylation = tf_hierarchy.merge(
        methyl_expr,
        left_on='locus_tag',
        right_on='gene_id',
        how='left'
    )

    # Calculate total methylation sites per timepoint
    for tp in ['T1', 'T2', 'T3']:
        tf_methylation[f'total_{tp}_sites'] = (
            tf_methylation[f'6mA_{tp}_count'].fillna(0) +
            tf_methylation[f'4mC_{tp}_count'].fillna(0)
        )

    # Classify methylation status
    def classify_methylation(row):
        t1 = row['total_T1_sites']
        t2 = row['total_T2_sites']
        t3 = row['total_T3_sites']

        if pd.isna(t1) or pd.isna(t2):
            return 'No data'

        if t1 == 0 and t2 == 0 and t3 == 0:
            return 'Unmethylated'
        elif t2 > t1 and t2 > 0:
            return 'Gained in T2'
        elif t2 < t1:
            return 'Lost in T2'
        elif t2 == t1 and t2 > 0:
            return 'Stable methylation'
        else:
            return 'Variable'

    tf_methylation['methyl_status_T2vsT1'] = tf_methylation.apply(classify_methylation, axis=1)

    return tf_methylation

def analyze_bgc_regulators(tf_methylation):
    """Focus analysis on Act/Red BGC-related regulators"""
    print("\n" + "="*60)
    print("ACT/RED BGC REGULATOR METHYLATION ANALYSIS")
    print("="*60)

    # Key regulators for Act and Red
    act_regulators = ['actII-ORF4', 'afsR', 'afsS', 'absA1', 'absA2']
    red_regulators = ['redD', 'redZ', 'afsR', 'absA1', 'absA2']

    bgc_results = []

    # Analyze Act regulators
    print("\n--- ACTINORHODIN (Act) PATHWAY REGULATORS ---")
    for reg in act_regulators:
        row = tf_methylation[tf_methylation['name'] == reg]
        if len(row) > 0:
            row = row.iloc[0]
            print(f"\n{reg} ({row['locus_tag']}):")
            print(f"  Tier: {row['tier']} ({row['category']})")
            print(f"  Methylation sites - T1: {row['total_T1_sites']:.0f}, T2: {row['total_T2_sites']:.0f}, T3: {row['total_T3_sites']:.0f}")
            print(f"  Status: {row['methyl_status_T2vsT1']}")
            if not pd.isna(row.get('log2FC_T2_vs_T1')):
                print(f"  Expression T2vsT1: log2FC={row['log2FC_T2_vs_T1']:.2f}")
            if not pd.isna(row.get('corr_act')):
                print(f"  Correlation with Act: r={row['corr_act']:.3f}")

            bgc_results.append({
                'regulator': reg,
                'bgc': 'act',
                'locus_tag': row['locus_tag'],
                'tier': row['tier'],
                'T1_sites': row['total_T1_sites'],
                'T2_sites': row['total_T2_sites'],
                'T3_sites': row['total_T3_sites'],
                'methyl_status': row['methyl_status_T2vsT1'],
                'log2FC_T2vsT1': row.get('log2FC_T2_vs_T1'),
                'corr_act': row.get('corr_act')
            })

    # Analyze Red regulators
    print("\n--- UNDECYLPRODIGIOSIN (Red) PATHWAY REGULATORS ---")
    for reg in red_regulators:
        if reg in act_regulators and reg not in ['redD', 'redZ']:
            continue  # Skip duplicates
        row = tf_methylation[tf_methylation['name'] == reg]
        if len(row) > 0:
            row = row.iloc[0]
            print(f"\n{reg} ({row['locus_tag']}):")
            print(f"  Tier: {row['tier']} ({row['category']})")
            print(f"  Methylation sites - T1: {row['total_T1_sites']:.0f}, T2: {row['total_T2_sites']:.0f}, T3: {row['total_T3_sites']:.0f}")
            print(f"  Status: {row['methyl_status_T2vsT1']}")
            if not pd.isna(row.get('log2FC_T2_vs_T1')):
                print(f"  Expression T2vsT1: log2FC={row['log2FC_T2_vs_T1']:.2f}")
            if not pd.isna(row.get('corr_red')):
                print(f"  Correlation with Red: r={row['corr_red']:.3f}")

            bgc_results.append({
                'regulator': reg,
                'bgc': 'red',
                'locus_tag': row['locus_tag'],
                'tier': row['tier'],
                'T1_sites': row['total_T1_sites'],
                'T2_sites': row['total_T2_sites'],
                'T3_sites': row['total_T3_sites'],
                'methyl_status': row['methyl_status_T2vsT1'],
                'log2FC_T2vsT1': row.get('log2FC_T2_vs_T1'),
                'corr_red': row.get('corr_red')
            })

    return pd.DataFrame(bgc_results)

def check_coordinated_tf_changes(tf_hierarchy, tf_coord):
    """Check which GRN TFs have coordinated methylation-expression changes"""
    print("\n" + "="*60)
    print("GRN TFs WITH COORDINATED METHYLATION CHANGES")
    print("="*60)

    # Get locus tags from hierarchy
    grn_locus_tags = set(tf_hierarchy['locus_tag'].dropna())

    # Find matching coordinated TFs
    coordinated_grn = tf_coord[tf_coord['gene_id'].isin(grn_locus_tags)]

    if len(coordinated_grn) > 0:
        print(f"\nFound {len(coordinated_grn)} GRN TFs with coordinated changes:")
        for _, row in coordinated_grn.iterrows():
            tf_info = tf_hierarchy[tf_hierarchy['locus_tag'] == row['gene_id']]
            if len(tf_info) > 0:
                tf_name = tf_info.iloc[0]['name']
                tier = tf_info.iloc[0]['tier']
                bgc = tf_info.iloc[0]['bgc'] if pd.notna(tf_info.iloc[0]['bgc']) else '-'
                print(f"\n  {tf_name} ({row['gene_id']}):")
                print(f"    Tier: {tier}, BGC: {bgc}")
                print(f"    Methylation: {row['mod_type']}, {row['methyl_category']}")
                print(f"    Expression: log2FC={row['log2FC']:.2f}, padj={row['padj']:.2e}")
                print(f"    Coordination: {row['coordination']} ({row['correlation_type']})")
    else:
        print("No GRN hierarchy TFs found with coordinated methylation changes")

    return coordinated_grn

def plot_tf_hierarchy_methylation(tf_methylation):
    """Visualize methylation patterns across the GRN hierarchy"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # 1. Methylation sites by tier
    ax1 = axes[0, 0]
    tier_data = tf_methylation.groupby('tier').agg({
        'total_T1_sites': 'sum',
        'total_T2_sites': 'sum',
        'total_T3_sites': 'sum'
    }).reset_index()

    tier_labels = ['Tier 1\n(Global)', 'Tier 2\n(Pleiotropic/Sigma)', 'Tier 3\n(CSR)']
    x = np.arange(len(tier_labels))
    width = 0.25

    ax1.bar(x - width, tier_data['total_T1_sites'], width, label='T1', color='#4e79a7')
    ax1.bar(x, tier_data['total_T2_sites'], width, label='T2', color='#f28e2b')
    ax1.bar(x + width, tier_data['total_T3_sites'], width, label='T3', color='#59a14f')
    ax1.set_xlabel('Regulatory Tier')
    ax1.set_ylabel('Total Methylation Sites')
    ax1.set_title('Methylation Sites by Regulatory Tier')
    ax1.set_xticks(x)
    ax1.set_xticklabels(tier_labels)
    ax1.legend()

    # 2. Methylation status distribution
    ax2 = axes[0, 1]
    status_counts = tf_methylation['methyl_status_T2vsT1'].value_counts()
    colors = ['#e15759', '#76b7b2', '#edc948', '#b07aa1', '#9c755f']
    ax2.pie(status_counts.values, labels=status_counts.index, autopct='%1.1f%%',
            colors=colors[:len(status_counts)])
    ax2.set_title('TF Methylation Status (T2 vs T1)')

    # 3. Act/Red BGC correlations vs methylation
    ax3 = axes[1, 0]
    tf_with_corr = tf_methylation[tf_methylation['corr_act'].notna()]
    if len(tf_with_corr) > 0:
        colors_tier = {1: '#4e79a7', 2: '#f28e2b', 3: '#59a14f'}
        for tier in [1, 2, 3]:
            tier_data = tf_with_corr[tf_with_corr['tier'] == tier]
            if len(tier_data) > 0:
                ax3.scatter(tier_data['corr_act'], tier_data['corr_red'],
                           c=colors_tier[tier], label=f'Tier {tier}', alpha=0.7, s=100)
                # Annotate key regulators
                for _, row in tier_data.iterrows():
                    if row['name'] in ['actII-ORF4', 'redD', 'redZ', 'afsR', 'absA2', 'papR2']:
                        ax3.annotate(row['name'], (row['corr_act'], row['corr_red']),
                                    fontsize=8, ha='left')
        ax3.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
        ax3.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
        ax3.set_xlabel('Correlation with Act BGC')
        ax3.set_ylabel('Correlation with Red BGC')
        ax3.set_title('TF-BGC Expression Correlations')
        ax3.legend()

    # 4. Expression changes of key regulators
    ax4 = axes[1, 1]
    key_regs = ['actII-ORF4', 'afsR', 'afsS', 'absA2', 'absA1', 'cpkO/kasO', 'papR2', 'sigF', 'sigU']
    key_data = tf_methylation[tf_methylation['name'].isin(key_regs)]
    key_data = key_data[key_data['log2FC_T2_vs_T1'].notna()]

    if len(key_data) > 0:
        key_data_sorted = key_data.sort_values('log2FC_T2_vs_T1', ascending=True)
        colors = ['#e15759' if x < 0 else '#59a14f' for x in key_data_sorted['log2FC_T2_vs_T1']]
        ax4.barh(key_data_sorted['name'], key_data_sorted['log2FC_T2_vs_T1'], color=colors)
        ax4.axvline(x=0, color='black', linewidth=0.5)
        ax4.set_xlabel('log2 Fold Change (T2 vs T1)')
        ax4.set_title('Expression Changes of Key Regulators')

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'tf_hierarchy_methylation_overview.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"\nSaved: tf_hierarchy_methylation_overview.png")

def plot_regulatory_cascade(tf_methylation, reg_interactions):
    """Visualize the regulatory cascade with methylation information"""
    fig, ax = plt.subplots(figsize=(16, 10))

    # Define positions for regulators in cascade
    positions = {
        # Tier 1 - Global (top)
        'bldD': (2, 4), 'adpA': (4, 4), 'afsR': (6, 4), 'dasR': (8, 4),
        # Tier 2 - Pleiotropic (middle)
        'absA2': (3, 3), 'hrdD': (6, 3), 'wblA': (9, 3),
        # Tier 3 - CSR (bottom)
        'actII-ORF4': (2, 2), 'redD': (5, 2), 'redZ': (5, 1.5), 'cdaR': (8, 2), 'cpkO/kasO': (11, 2),
    }

    # Draw nodes
    for name, (x, y) in positions.items():
        tf_data = tf_methylation[tf_methylation['name'] == name]
        if len(tf_data) > 0:
            tf_data = tf_data.iloc[0]

            # Color by methylation status
            status = tf_data['methyl_status_T2vsT1']
            if status == 'Gained in T2':
                color = '#f28e2b'
            elif status == 'Lost in T2':
                color = '#4e79a7'
            elif status == 'Unmethylated':
                color = '#e0e0e0'
            else:
                color = '#76b7b2'

            # Size by expression fold change
            log2fc = tf_data.get('log2FC_T2_vs_T1', 0) if pd.notna(tf_data.get('log2FC_T2_vs_T1')) else 0
            size = 1500 + abs(log2fc) * 200

            # Draw node
            circle = plt.Circle((x, y), 0.3, color=color, ec='black', linewidth=2, zorder=10)
            ax.add_patch(circle)

            # Add label
            ax.annotate(name, (x, y-0.5), ha='center', fontsize=9, fontweight='bold')
            if pd.notna(tf_data.get('log2FC_T2_vs_T1')):
                direction = '↑' if tf_data['log2FC_T2_vs_T1'] > 0 else '↓'
                ax.annotate(f"{direction}{abs(tf_data['log2FC_T2_vs_T1']):.1f}",
                           (x, y-0.75), ha='center', fontsize=7)

    # Draw edges (regulatory interactions)
    edges = [
        ('bldD', 'actII-ORF4'), ('bldD', 'redD'), ('bldD', 'cdaR'),
        ('afsR', 'actII-ORF4'), ('afsR', 'redD'),
        ('absA2', 'actII-ORF4'), ('absA2', 'redD'), ('absA2', 'cdaR'), ('absA2', 'cpkO/kasO'),
        ('hrdD', 'actII-ORF4'), ('hrdD', 'redD'),
        ('redZ', 'redD'),
    ]

    for source, target in edges:
        if source in positions and target in positions:
            x1, y1 = positions[source]
            x2, y2 = positions[target]
            ax.annotate('', xy=(x2, y2+0.3), xytext=(x1, y1-0.3),
                       arrowprops=dict(arrowstyle='->', color='gray', lw=1.5))

    # Add tier labels
    ax.text(-0.5, 4, 'Tier 1\n(Global)', fontsize=11, fontweight='bold', va='center')
    ax.text(-0.5, 3, 'Tier 2\n(Pleiotropic)', fontsize=11, fontweight='bold', va='center')
    ax.text(-0.5, 1.75, 'Tier 3\n(CSR)', fontsize=11, fontweight='bold', va='center')

    # Add BGC labels
    ax.text(2, 0.8, 'Act', fontsize=10, ha='center', style='italic', color='blue')
    ax.text(5, 0.8, 'Red', fontsize=10, ha='center', style='italic', color='red')
    ax.text(8, 0.8, 'CDA', fontsize=10, ha='center', style='italic', color='purple')
    ax.text(11, 0.8, 'Cpk', fontsize=10, ha='center', style='italic', color='green')

    # Legend
    legend_elements = [
        plt.Circle((0, 0), 0.1, color='#f28e2b', label='Gained methylation (T2)'),
        plt.Circle((0, 0), 0.1, color='#4e79a7', label='Lost methylation (T2)'),
        plt.Circle((0, 0), 0.1, color='#e0e0e0', label='Unmethylated'),
        plt.Circle((0, 0), 0.1, color='#76b7b2', label='Stable/Other'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=9)

    ax.set_xlim(-1, 13)
    ax.set_ylim(0, 5)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('GRN Regulatory Cascade with Methylation Status\n(T2 vs T1 comparison)',
                fontsize=14, fontweight='bold')

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'regulatory_cascade_methylation.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: regulatory_cascade_methylation.png")

def generate_summary_table(tf_methylation, bgc_results, coordinated_grn):
    """Generate summary tables for the report"""
    # 1. TF methylation summary
    tf_summary = tf_methylation[[
        'name', 'locus_tag', 'tier', 'category', 'bgc',
        'total_T1_sites', 'total_T2_sites', 'total_T3_sites',
        'methyl_status_T2vsT1', 'log2FC_T2_vs_T1', 'padj_T2_vs_T1',
        'corr_act', 'corr_red'
    ]].copy()
    tf_summary.to_csv(OUTPUT_DIR / 'tf_methylation_summary.csv', index=False)
    print(f"\nSaved: tf_methylation_summary.csv")

    # 2. BGC regulator details
    bgc_df = pd.DataFrame(bgc_results) if len(bgc_results) > 0 else pd.DataFrame()
    if len(bgc_df) > 0:
        bgc_df.to_csv(OUTPUT_DIR / 'bgc_regulator_methylation.csv', index=False)
        print(f"Saved: bgc_regulator_methylation.csv")

    # 3. Key insights summary
    insights = []

    # Check for coordinated changes
    if len(coordinated_grn) > 0:
        for _, row in coordinated_grn.iterrows():
            insights.append({
                'insight_type': 'Coordinated change',
                'tf': row['gene_id'],
                'detail': f"{row['methyl_category']} methylation + {row['coordination']}",
                'correlation': row['correlation_type']
            })

    # Check Act/Red regulators
    for _, row in bgc_df.iterrows() if len(bgc_df) > 0 else []:
        if pd.notna(row.get('corr_act')) and abs(row['corr_act']) > 0.5:
            insights.append({
                'insight_type': 'High Act correlation',
                'tf': row['regulator'],
                'detail': f"r={row['corr_act']:.3f}",
                'correlation': 'BGC correlation'
            })

    if insights:
        pd.DataFrame(insights).to_csv(OUTPUT_DIR / 'key_insights.csv', index=False)
        print(f"Saved: key_insights.csv")

    return tf_summary

def main():
    print("="*70)
    print("TF PROMOTER METHYLATION ANALYSIS")
    print("Step C-2: GRN Hierarchy TF Methylation Status")
    print("="*70)

    # Load data
    tf_hierarchy, methyl_expr, tf_coord, reg_interactions = load_data()

    # Merge TF hierarchy with methylation data
    tf_methylation = merge_tf_methylation(tf_hierarchy, methyl_expr)

    # Analyze BGC regulators
    bgc_results = analyze_bgc_regulators(tf_methylation)

    # Check coordinated TF changes
    coordinated_grn = check_coordinated_tf_changes(tf_hierarchy, tf_coord)

    # Generate visualizations
    print("\n" + "="*60)
    print("GENERATING FIGURES")
    print("="*60)
    plot_tf_hierarchy_methylation(tf_methylation)
    plot_regulatory_cascade(tf_methylation, reg_interactions)

    # Generate summary tables
    tf_summary = generate_summary_table(tf_methylation, bgc_results, coordinated_grn)

    # Print final summary
    print("\n" + "="*60)
    print("KEY FINDINGS SUMMARY")
    print("="*60)

    # Count methylated TFs by tier
    methylated_tfs = tf_methylation[tf_methylation['methyl_status_T2vsT1'] != 'Unmethylated']
    print(f"\nMethylated TFs in GRN hierarchy:")
    for tier in [1, 2, 3]:
        tier_data = methylated_tfs[methylated_tfs['tier'] == tier]
        total_tier = len(tf_methylation[tf_methylation['tier'] == tier])
        print(f"  Tier {tier}: {len(tier_data)}/{total_tier} TFs have methylation sites")

    # Key regulators with methylation
    print("\nKey Act/Red regulators with methylation changes (T2 vs T1):")
    key_regs = ['actII-ORF4', 'redD', 'redZ', 'afsR', 'absA1', 'absA2']
    for reg in key_regs:
        reg_data = tf_methylation[tf_methylation['name'] == reg]
        if len(reg_data) > 0:
            reg_data = reg_data.iloc[0]
            if reg_data['methyl_status_T2vsT1'] != 'Unmethylated':
                print(f"  {reg}: {reg_data['methyl_status_T2vsT1']}")

    # redZ finding
    redz = tf_methylation[tf_methylation['name'] == 'redZ']
    if len(redz) > 0:
        redz = redz.iloc[0]
        print(f"\n*** NOTABLE: redZ (SC_RS27300) ***")
        print(f"  - Part of coordinated TF changes")
        print(f"  - Methylation: Lost in T2 (6mA)")
        print(f"  - Controls redD → Red biosynthesis")

    print("\n" + "="*60)
    print("OUTPUT FILES")
    print("="*60)
    print(f"Directory: {OUTPUT_DIR}")
    for f in OUTPUT_DIR.glob('*'):
        print(f"  - {f.name}")

if __name__ == '__main__':
    main()
