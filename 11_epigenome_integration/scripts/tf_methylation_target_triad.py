#!/usr/bin/env python3
"""
TF-Methylation-Target Triad Analysis
Step C-3: Connects TFs, their methylation status, and downstream targets

Key question: Does TF methylation affect target gene expression?
Focus: Act and Red BGC pathways
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import seaborn as sns
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Setup
BASE_DIR = Path('/Users/okaban/bioinfo/rna-seq/11_epigenome_integration')
ANALYSIS_DIR = BASE_DIR / 'analysis/12_grn_tf_methylation'
INTEGRATION_DIR = BASE_DIR / 'analysis/01_integration'
OUTPUT_DIR = ANALYSIS_DIR / 'triad_analysis'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# BGC gene coordinates (from previous analysis)
BGC_GENES = {
    'act': {
        'range': (5532640, 5558199),
        'key_genes': ['actI-ORF1', 'actI-ORF2', 'actI-ORF3', 'actII-ORF4', 'actVI-ORF1', 'actVII'],
        'csr': 'actII-ORF4',
        'locus_prefix': 'SC_RS2755', # approximate
    },
    'red': {
        'range': (5476648, 5511923),
        'key_genes': ['redD', 'redZ', 'redL', 'redK', 'redM', 'redN', 'redO', 'redP', 'redQ', 'redR'],
        'csr': 'redD',
        'locus_prefix': 'SC_RS2720',
    },
    'cda': {
        'range': (3495761, 3534082),
        'key_genes': ['cdaR', 'cdaPS1', 'cdaPS2', 'cdaPS3'],
        'csr': 'cdaR',
        'locus_prefix': 'SC_RS1783',
    },
    'cpk': {
        'range': (6410614, 6462648),
        'key_genes': ['cpkO', 'cpkA', 'cpkB', 'cpkC'],
        'csr': 'cpkO/kasO',
        'locus_prefix': 'SC_RS3155',
    }
}

# Known regulatory relationships from literature
REGULATORY_NETWORK = {
    # TF -> list of known targets
    'redZ': ['redD'],  # redZ directly activates redD
    'redD': ['redL', 'redK', 'redM', 'redN', 'redO', 'redP', 'redQ'],  # redD activates red biosynthetic genes
    'actII-ORF4': ['actI', 'actVI', 'actVII', 'actVA', 'actVB'],  # SARP activates act genes
    'absA2': ['actII-ORF4', 'redD', 'cdaR', 'cpkO'],  # negative regulator of multiple CSRs
    'afsR': ['actII-ORF4', 'redD', 'afsS'],  # global regulator
    'bldD': ['actII-ORF4', 'redD', 'cdaR'],  # developmental control of BGCs
}

def load_data():
    """Load all required datasets"""
    # TF methylation summary
    tf_methyl = pd.read_csv(ANALYSIS_DIR / 'tf_promoter_methylation/tf_methylation_summary.csv')
    print(f"Loaded {len(tf_methyl)} TF methylation records")

    # Full expression data
    deseq_dir = Path('/Users/okaban/bioinfo/rna-seq/04_deseq2/analysis/04_deseq2_260128_v1/results')
    t2vst1 = pd.read_csv(deseq_dir / 'DESeq2_M145_2_vs_1.tsv', sep='\t')
    t3vst1 = pd.read_csv(deseq_dir / 'DESeq2_M145_3_vs_1.tsv', sep='\t')
    print(f"Loaded DESeq2 results: T2vsT1 ({len(t2vst1)}), T3vsT1 ({len(t3vst1)})")

    # Methylation data
    methyl_expr = pd.read_csv(INTEGRATION_DIR / 'integrated_methyl_expression_weighted.csv')
    print(f"Loaded methylation data for {len(methyl_expr)} genes")

    # Regulatory interactions
    reg_interactions = pd.read_csv(ANALYSIS_DIR / 'regulatory_interactions.csv')
    print(f"Loaded {len(reg_interactions)} regulatory interactions")

    return tf_methyl, t2vst1, t3vst1, methyl_expr, reg_interactions

def identify_bgc_genes(t2vst1):
    """Identify genes in each BGC"""
    bgc_gene_lists = {}

    for bgc_name, bgc_info in BGC_GENES.items():
        # Find genes by locus prefix (approximate method)
        prefix = bgc_info['locus_prefix']
        bgc_genes = t2vst1[t2vst1['gene_id'].str.startswith(prefix, na=False)]

        # Also filter by genomic coordinates if available
        start, end = bgc_info['range']

        bgc_gene_lists[bgc_name] = {
            'genes': bgc_genes['gene_id'].tolist(),
            'n_genes': len(bgc_genes),
            'mean_log2fc_t2': bgc_genes['log2FoldChange'].mean() if len(bgc_genes) > 0 else np.nan,
            'deg_up': (bgc_genes['log2FoldChange'] > 1).sum() if len(bgc_genes) > 0 else 0,
            'deg_down': (bgc_genes['log2FoldChange'] < -1).sum() if len(bgc_genes) > 0 else 0,
        }

    return bgc_gene_lists

def build_triad_table(tf_methyl, t2vst1, methyl_expr):
    """Build the TF-Methylation-Target triad table"""
    triads = []

    for tf_name, targets in REGULATORY_NETWORK.items():
        # Get TF methylation info
        tf_info = tf_methyl[tf_methyl['name'] == tf_name]
        if len(tf_info) == 0:
            continue
        tf_info = tf_info.iloc[0]

        # Get TF expression
        tf_expr = t2vst1[t2vst1['gene_id'] == tf_info['locus_tag']]
        tf_log2fc = tf_expr['log2FoldChange'].iloc[0] if len(tf_expr) > 0 else np.nan
        tf_padj = tf_expr['padj'].iloc[0] if len(tf_expr) > 0 else np.nan

        # Get TF methylation sites
        tf_methyl_t1 = tf_info['total_T1_sites']
        tf_methyl_t2 = tf_info['total_T2_sites']
        tf_methyl_status = tf_info['methyl_status_T2vsT1']

        for target in targets:
            # Try to find target in expression data
            # First try exact locus_tag match
            target_locus = None
            target_expr = None

            # Check if target is a locus tag directly
            if target.startswith('SC_RS'):
                target_locus = target
            else:
                # Try to find by gene name (approximate)
                for bgc_name, bgc_info in BGC_GENES.items():
                    if target in bgc_info['key_genes'] or target == bgc_info['csr']:
                        # Find in TF list first
                        tf_match = tf_methyl[tf_methyl['name'] == target]
                        if len(tf_match) > 0:
                            target_locus = tf_match.iloc[0]['locus_tag']
                            break

            if target_locus:
                target_expr = t2vst1[t2vst1['gene_id'] == target_locus]

            target_log2fc = target_expr['log2FoldChange'].iloc[0] if target_expr is not None and len(target_expr) > 0 else np.nan
            target_padj = target_expr['padj'].iloc[0] if target_expr is not None and len(target_expr) > 0 else np.nan

            # Get target methylation
            target_methyl = methyl_expr[methyl_expr['gene_id'] == target_locus] if target_locus else None
            target_methyl_t1 = target_methyl['6mA_T1_count'].iloc[0] + target_methyl['4mC_T1_count'].iloc[0] if target_methyl is not None and len(target_methyl) > 0 else 0
            target_methyl_t2 = target_methyl['6mA_T2_count'].iloc[0] + target_methyl['4mC_T2_count'].iloc[0] if target_methyl is not None and len(target_methyl) > 0 else 0

            triads.append({
                'tf_name': tf_name,
                'tf_locus': tf_info['locus_tag'],
                'tf_tier': tf_info['tier'],
                'tf_methyl_T1': tf_methyl_t1,
                'tf_methyl_T2': tf_methyl_t2,
                'tf_methyl_status': tf_methyl_status,
                'tf_log2fc': tf_log2fc,
                'tf_padj': tf_padj,
                'target_name': target,
                'target_locus': target_locus,
                'target_methyl_T1': target_methyl_t1,
                'target_methyl_T2': target_methyl_t2,
                'target_log2fc': target_log2fc,
                'target_padj': target_padj,
                'regulation_type': 'activation' if tf_name not in ['absA2', 'nsdA', 'nsdB', 'wblA'] else 'repression'
            })

    return pd.DataFrame(triads)

def analyze_redz_cascade(tf_methyl, t2vst1, methyl_expr):
    """Detailed analysis of the redZ → redD → Red BGC cascade"""
    print("\n" + "="*60)
    print("redZ → redD → Red BGC CASCADE ANALYSIS")
    print("="*60)

    cascade_results = {}

    # 1. redZ (methylated TF)
    redz_info = tf_methyl[tf_methyl['name'] == 'redZ'].iloc[0]
    redz_expr = t2vst1[t2vst1['gene_id'] == 'SC_RS27300']

    cascade_results['redZ'] = {
        'locus_tag': 'SC_RS27300',
        'methyl_T1': redz_info['total_T1_sites'],
        'methyl_T2': redz_info['total_T2_sites'],
        'methyl_status': redz_info['methyl_status_T2vsT1'],
        'log2fc': redz_expr['log2FoldChange'].iloc[0] if len(redz_expr) > 0 else np.nan,
        'padj': redz_expr['padj'].iloc[0] if len(redz_expr) > 0 else np.nan,
    }

    print(f"\n1. redZ (Response Regulator):")
    print(f"   Locus: SC_RS27300")
    print(f"   Methylation: T1={cascade_results['redZ']['methyl_T1']:.0f}, T2={cascade_results['redZ']['methyl_T2']:.0f} sites")
    print(f"   Status: {cascade_results['redZ']['methyl_status']}")
    print(f"   Expression: log2FC={cascade_results['redZ']['log2fc']:.2f}, padj={cascade_results['redZ']['padj']:.2e}")
    print(f"   Interpretation: 6mA methylation LOST → expression DOWN")

    # 2. redD (direct target of redZ)
    redd_info = tf_methyl[tf_methyl['name'] == 'redD'].iloc[0]
    redd_expr = t2vst1[t2vst1['gene_id'] == 'SC_RS27225']

    cascade_results['redD'] = {
        'locus_tag': 'SC_RS27225',
        'methyl_T1': redd_info['total_T1_sites'],
        'methyl_T2': redd_info['total_T2_sites'],
        'methyl_status': redd_info['methyl_status_T2vsT1'],
        'log2fc': redd_expr['log2FoldChange'].iloc[0] if len(redd_expr) > 0 else np.nan,
        'padj': redd_expr['padj'].iloc[0] if len(redd_expr) > 0 else np.nan,
    }

    print(f"\n2. redD (SARP - Red Activator):")
    print(f"   Locus: SC_RS27225")
    print(f"   Methylation: T1={cascade_results['redD']['methyl_T1']:.0f}, T2={cascade_results['redD']['methyl_T2']:.0f} sites")
    print(f"   Status: {cascade_results['redD']['methyl_status']}")
    print(f"   Expression: log2FC={cascade_results['redD']['log2fc']:.2f}, padj={cascade_results['redD']['padj']:.2e}")
    print(f"   Regulated by: redZ (positive), AbsA2 (negative)")

    # 3. Red BGC genes
    print(f"\n3. Red BGC Biosynthetic Genes:")
    red_genes = ['SC_RS27210', 'SC_RS27215', 'SC_RS27220', 'SC_RS27230', 'SC_RS27235',
                 'SC_RS27240', 'SC_RS27245', 'SC_RS27250', 'SC_RS27255', 'SC_RS27260']

    cascade_results['red_bgc'] = {'genes': [], 'mean_log2fc': 0, 'n_up': 0, 'n_down': 0}

    for gene in red_genes:
        expr = t2vst1[t2vst1['gene_id'] == gene]
        if len(expr) > 0:
            log2fc = expr['log2FoldChange'].iloc[0]
            cascade_results['red_bgc']['genes'].append({
                'locus': gene,
                'log2fc': log2fc,
                'padj': expr['padj'].iloc[0]
            })
            if log2fc > 1:
                cascade_results['red_bgc']['n_up'] += 1
            elif log2fc < -1:
                cascade_results['red_bgc']['n_down'] += 1

    if cascade_results['red_bgc']['genes']:
        mean_fc = np.mean([g['log2fc'] for g in cascade_results['red_bgc']['genes']])
        cascade_results['red_bgc']['mean_log2fc'] = mean_fc
        print(f"   Found {len(cascade_results['red_bgc']['genes'])} genes in region")
        print(f"   Mean log2FC: {mean_fc:.2f}")
        print(f"   DEGs: {cascade_results['red_bgc']['n_up']} up, {cascade_results['red_bgc']['n_down']} down")

    return cascade_results

def plot_triad_network(triads, cascade_results):
    """Visualize the TF-Methylation-Target network"""
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))

    # Left panel: Triad network
    ax1 = axes[0]

    # Position nodes
    positions = {
        # TFs with methylation
        'redZ': (0, 2),
        'afsR': (2, 2),
        'absA2': (4, 2),
        # CSR targets
        'redD': (0, 1),
        'actII-ORF4': (2, 1),
        # BGC endpoints
        'Red BGC': (0, 0),
        'Act BGC': (2, 0),
    }

    # Draw nodes
    node_info = {
        'redZ': {'color': '#e15759', 'methyl': 'Lost 6mA'},
        'afsR': {'color': '#76b7b2', 'methyl': 'Stable'},
        'absA2': {'color': '#e0e0e0', 'methyl': 'None'},
        'redD': {'color': '#e0e0e0', 'methyl': 'None'},
        'actII-ORF4': {'color': '#e0e0e0', 'methyl': 'None'},
        'Red BGC': {'color': '#f28e2b', 'methyl': 'Pathway'},
        'Act BGC': {'color': '#4e79a7', 'methyl': 'Pathway'},
    }

    for node, (x, y) in positions.items():
        info = node_info[node]
        circle = plt.Circle((x, y), 0.3, color=info['color'], ec='black', linewidth=2)
        ax1.add_patch(circle)
        ax1.annotate(node, (x, y+0.4), ha='center', fontsize=9, fontweight='bold')
        ax1.annotate(info['methyl'], (x, y-0.45), ha='center', fontsize=7, style='italic')

    # Draw edges
    edges = [
        ('redZ', 'redD', '+', 'green'),
        ('afsR', 'redD', '+', 'green'),
        ('afsR', 'actII-ORF4', '+', 'green'),
        ('absA2', 'redD', '-', 'red'),
        ('absA2', 'actII-ORF4', '-', 'red'),
        ('redD', 'Red BGC', '+', 'green'),
        ('actII-ORF4', 'Act BGC', '+', 'green'),
    ]

    for src, tgt, sign, color in edges:
        x1, y1 = positions[src]
        x2, y2 = positions[tgt]
        ax1.annotate('', xy=(x2, y2+0.3), xytext=(x1, y1-0.3),
                    arrowprops=dict(arrowstyle='->', color=color, lw=2))
        # Add +/- label
        mid_x, mid_y = (x1+x2)/2, (y1+y2)/2
        ax1.annotate(sign, (mid_x+0.1, mid_y), fontsize=10, color=color, fontweight='bold')

    ax1.set_xlim(-1, 5)
    ax1.set_ylim(-0.8, 3)
    ax1.set_aspect('equal')
    ax1.axis('off')
    ax1.set_title('TF-Methylation-Target Regulatory Network', fontsize=12, fontweight='bold')

    # Legend
    legend_elements = [
        mpatches.Patch(color='#e15759', label='Methylation Lost'),
        mpatches.Patch(color='#76b7b2', label='Stable Methylation'),
        mpatches.Patch(color='#e0e0e0', label='No Methylation'),
        Line2D([0], [0], color='green', lw=2, label='Activation (+)'),
        Line2D([0], [0], color='red', lw=2, label='Repression (-)'),
    ]
    ax1.legend(handles=legend_elements, loc='lower right', fontsize=8)

    # Right panel: redZ cascade detail
    ax2 = axes[1]

    # Draw cascade flowchart
    cascade_positions = {
        'SC_RS17645\n(N-6 MTase)': (1, 4),
        '6mA\nmethylation': (1, 3),
        'redZ\npromoter': (1, 2),
        'redZ\nexpression': (1, 1),
        'redD\nexpression': (1, 0),
    }

    # Expression changes
    expr_colors = {
        'SC_RS17645\n(N-6 MTase)': '#e15759',  # Down
        '6mA\nmethylation': '#e15759',  # Lost
        'redZ\npromoter': '#4e79a7',  # Demethylated
        'redZ\nexpression': '#e15759',  # Down
        'redD\nexpression': '#59a14f',  # Up
    }

    for node, (x, y) in cascade_positions.items():
        rect = plt.Rectangle((x-0.4, y-0.25), 0.8, 0.5,
                             color=expr_colors[node], ec='black', linewidth=2)
        ax2.add_patch(rect)
        ax2.annotate(node, (x, y), ha='center', va='center', fontsize=9)

    # Draw arrows between cascade steps
    for i, (node1, node2) in enumerate(zip(list(cascade_positions.keys())[:-1],
                                           list(cascade_positions.keys())[1:])):
        x1, y1 = cascade_positions[node1]
        x2, y2 = cascade_positions[node2]
        ax2.annotate('', xy=(x2, y2+0.25), xytext=(x1, y1-0.25),
                    arrowprops=dict(arrowstyle='->', color='black', lw=2))

    # Add annotations
    ax2.annotate('log2FC = -2.19\n(T2 vs T1)', (2, 4), fontsize=8)
    ax2.annotate('Lost\n(T2 vs T1)', (2, 3), fontsize=8)
    ax2.annotate('6mA site\ndemethylated', (2, 2), fontsize=8)
    ax2.annotate('log2FC = -2.25\n(coordinated)', (2, 1), fontsize=8)
    ax2.annotate('log2FC = +4.77\n(paradox?)', (2, 0), fontsize=8)

    ax2.set_xlim(-0.5, 3)
    ax2.set_ylim(-1, 5)
    ax2.axis('off')
    ax2.set_title('Proposed Epigenetic Cascade:\nMTase → Methylation → TF → Target',
                 fontsize=12, fontweight='bold')

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'triad_network_visualization.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"\nSaved: triad_network_visualization.png")

def plot_bgc_expression_summary(t2vst1):
    """Summarize expression changes in all BGCs"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    bgc_colors = {'act': '#4e79a7', 'red': '#f28e2b', 'cda': '#59a14f', 'cpk': '#e15759'}

    for idx, (bgc_name, bgc_info) in enumerate(BGC_GENES.items()):
        ax = axes[idx // 2, idx % 2]

        # Get genes in this BGC region
        prefix = bgc_info['locus_prefix']
        bgc_genes = t2vst1[t2vst1['gene_id'].str.startswith(prefix, na=False)].copy()

        if len(bgc_genes) > 0:
            bgc_genes = bgc_genes.sort_values('log2FoldChange')

            colors = ['#59a14f' if x > 0 else '#e15759' for x in bgc_genes['log2FoldChange']]
            ax.barh(range(len(bgc_genes)), bgc_genes['log2FoldChange'], color=colors)
            ax.axvline(x=0, color='black', linewidth=0.5)
            ax.set_yticks(range(len(bgc_genes)))
            ax.set_yticklabels(bgc_genes['gene_id'].str.replace('SC_RS', ''), fontsize=7)
            ax.set_xlabel('log2 Fold Change (T2 vs T1)')
            ax.set_title(f'{bgc_name.upper()} BGC Expression')

            # Add summary stats
            mean_fc = bgc_genes['log2FoldChange'].mean()
            n_up = (bgc_genes['log2FoldChange'] > 1).sum()
            n_down = (bgc_genes['log2FoldChange'] < -1).sum()
            ax.text(0.95, 0.95, f'Mean: {mean_fc:.2f}\nUp: {n_up}, Down: {n_down}',
                   transform=ax.transAxes, ha='right', va='top', fontsize=9,
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'bgc_expression_summary.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: bgc_expression_summary.png")

def generate_report_data(triads, cascade_results, tf_methyl):
    """Generate data for the final report"""
    # Key insights
    insights = []

    # Insight 1: redZ coordinated methylation-expression
    insights.append({
        'number': 1,
        'insight': 'redZはメチル化-発現協調変動を示す唯一のGRN TF',
        'evidence': '6mA Lost + expression DOWN (log2FC=-2.25)',
        'figure': 'triad_network_visualization.png',
        'significance': 'Red BGCへのエピジェネティック制御の直接的証拠'
    })

    # Insight 2: MTase expression correlates with methylation
    insights.append({
        'number': 2,
        'insight': 'SC_RS17645 (N-6 MTase) 発現低下がredZメチル化消失と一致',
        'evidence': 'MTase log2FC=-2.19 (T2vsT1)',
        'figure': 'regulatory_cascade_methylation.png',
        'significance': 'AAGCCCGモチーフのメチル化機構解明'
    })

    # Insight 3: redD paradox
    insights.append({
        'number': 3,
        'insight': 'redDはredZ発現低下にもかかわらず高発現',
        'evidence': 'redZ: -2.25, redD: +4.77 (log2FC)',
        'figure': 'bgc_expression_summary.png',
        'significance': 'AbsA2抑制解除の影響が支配的'
    })

    # Save insights
    pd.DataFrame(insights).to_csv(OUTPUT_DIR / 'key_insights_triad.csv', index=False)

    # Save triads table
    triads.to_csv(OUTPUT_DIR / 'tf_target_triads.csv', index=False)

    print(f"\nSaved: key_insights_triad.csv")
    print(f"Saved: tf_target_triads.csv")

    return insights

def main():
    print("="*70)
    print("TF-METHYLATION-TARGET TRIAD ANALYSIS")
    print("Step C-3: Connecting TFs, Methylation, and Downstream Targets")
    print("="*70)

    # Load data
    tf_methyl, t2vst1, t3vst1, methyl_expr, reg_interactions = load_data()

    # Build triad table
    print("\n" + "="*60)
    print("BUILDING TRIAD RELATIONSHIPS")
    print("="*60)
    triads = build_triad_table(tf_methyl, t2vst1, methyl_expr)
    print(f"\nBuilt {len(triads)} TF-target triads")

    # Analyze redZ cascade in detail
    cascade_results = analyze_redz_cascade(tf_methyl, t2vst1, methyl_expr)

    # Generate visualizations
    print("\n" + "="*60)
    print("GENERATING VISUALIZATIONS")
    print("="*60)
    plot_triad_network(triads, cascade_results)
    plot_bgc_expression_summary(t2vst1)

    # Generate report data
    insights = generate_report_data(triads, cascade_results, tf_methyl)

    # Print summary
    print("\n" + "="*60)
    print("KEY FINDINGS SUMMARY")
    print("="*60)

    print("\n1. EPIGENETIC CONTROL OF RED PATHWAY:")
    print("   redZ (SC_RS27300) shows coordinated methylation-expression:")
    print("   - 6mA methylation LOST at promoter (T2 vs T1)")
    print("   - Expression DOWN (log2FC = -2.25)")
    print("   - This is the ONLY GRN TF with such coordination")

    print("\n2. PROPOSED MECHANISM:")
    print("   SC_RS17645 (N-6 MTase) expression DOWN (T2)")
    print("   → AAGCCCG methylation reduced")
    print("   → redZ promoter demethylated")
    print("   → redZ expression altered")
    print("   → Downstream effects on Red BGC")

    print("\n3. PARADOX - redD activation despite redZ down:")
    print("   - AbsA2: log2FC = +6.29 (massive increase)")
    print("   - AbsA2 is a NEGATIVE regulator of redD")
    print("   - Yet redD is UP (log2FC = +4.77)")
    print("   - Suggests other activators or AbsA2 saturation")

    print("\n" + "="*60)
    print("OUTPUT FILES")
    print("="*60)
    print(f"Directory: {OUTPUT_DIR}")
    for f in OUTPUT_DIR.glob('*'):
        print(f"  - {f.name}")

if __name__ == '__main__':
    main()
