#!/usr/bin/env python3
"""
Epigenetic Regulatory Network Visualization
Step C-4: Comprehensive network visualization integrating:
- GRN hierarchy (Global → Pleiotropic → CSR)
- Methylation status (4mC, 6mA)
- Expression changes
- BGC targets
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Setup
BASE_DIR = Path('/Users/okaban/bioinfo/rna-seq/11_epigenome_integration')
ANALYSIS_DIR = BASE_DIR / 'analysis/12_grn_tf_methylation'
OUTPUT_DIR = ANALYSIS_DIR / 'network_visualization'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def load_all_data():
    """Load all required data for visualization"""
    # TF methylation summary
    tf_methyl = pd.read_csv(ANALYSIS_DIR / 'tf_promoter_methylation/tf_methylation_summary.csv')

    # Regulatory interactions
    reg_interactions = pd.read_csv(ANALYSIS_DIR / 'regulatory_interactions.csv')

    # Triad data
    triads = pd.read_csv(ANALYSIS_DIR / 'triad_analysis/tf_target_triads.csv')

    # MTase data
    mtase = pd.read_csv(BASE_DIR / 'analysis/11_rm_system_identification/mtase_genes_with_expression.csv')

    return tf_methyl, reg_interactions, triads, mtase

def create_comprehensive_network(tf_methyl, reg_interactions):
    """Create comprehensive multi-layer network visualization"""
    fig = plt.figure(figsize=(20, 16))

    # Main network panel
    ax_main = fig.add_axes([0.05, 0.25, 0.65, 0.70])

    # Define node positions (hierarchical layout)
    positions = {
        # Layer 1: MTase (top) - Epigenetic machinery
        'SC_RS17645': (10, 8.5),  # N-6 DNA methylase

        # Layer 2: Global regulators (Tier 1)
        'bldA': (2, 7), 'bldD': (5, 7), 'adpA': (8, 7), 'afsR': (11, 7),
        'dasR': (14, 7), 'afsS': (17, 7),

        # Layer 3: Pleiotropic regulators (Tier 2)
        'absA1': (4, 5.5), 'absA2': (7, 5.5), 'afsQ1': (10, 5.5), 'wblA': (13, 5.5),
        'hrdD': (16, 5.5),

        # Layer 4: CSR (Tier 3)
        'actII-ORF4': (3, 4), 'redZ': (7, 4), 'redD': (10, 4), 'cdaR': (13, 4),
        'cpkO/kasO': (16, 4), 'papR2': (19, 4),

        # Layer 5: BGC endpoints (bottom)
        'Act': (3, 2), 'Red': (8.5, 2), 'CDA': (13, 2), 'Cpk': (16, 2),
    }

    # Draw tier separators
    tier_labels = [
        (8.2, 'MTase'),
        (7, 'Tier 1: Global'),
        (5.5, 'Tier 2: Pleiotropic/Sigma'),
        (4, 'Tier 3: CSR'),
        (2, 'BGC Products'),
    ]
    for y, label in tier_labels:
        ax_main.axhline(y=y-0.6, color='lightgray', linestyle='--', alpha=0.5, zorder=0)
        ax_main.text(0.3, y, label, fontsize=10, fontweight='bold', va='center', color='gray')

    # Draw nodes with methylation status
    def get_node_style(name, tf_methyl):
        """Determine node color based on methylation and expression"""
        if name in ['Act', 'Red', 'CDA', 'Cpk']:
            return {'color': '#f0f0f0', 'size': 0.8, 'shape': 'rect'}

        if name == 'SC_RS17645':
            return {'color': '#9467bd', 'size': 0.5, 'shape': 'circle'}  # MTase purple

        tf_data = tf_methyl[tf_methyl['name'] == name]
        if len(tf_data) == 0:
            return {'color': '#e0e0e0', 'size': 0.4, 'shape': 'circle'}

        tf_data = tf_data.iloc[0]
        methyl_status = tf_data['methyl_status_T2vsT1']
        log2fc = tf_data.get('log2FC_T2_vs_T1', 0)
        log2fc = 0 if pd.isna(log2fc) else log2fc

        # Color by methylation status
        if methyl_status == 'Lost in T2':
            color = '#e15759'  # Red for lost
        elif methyl_status == 'Gained in T2':
            color = '#59a14f'  # Green for gained
        elif methyl_status == 'Stable methylation':
            color = '#76b7b2'  # Teal for stable
        else:
            color = '#e0e0e0'  # Gray for unmethylated

        # Size by expression change
        size = 0.35 + min(abs(log2fc), 5) * 0.08

        return {'color': color, 'size': size, 'shape': 'circle'}

    # Draw all nodes
    for name, (x, y) in positions.items():
        style = get_node_style(name, tf_methyl)

        if style['shape'] == 'circle':
            circle = plt.Circle((x, y), style['size'], color=style['color'],
                               ec='black', linewidth=1.5, zorder=10)
            ax_main.add_patch(circle)
        else:
            rect = mpatches.FancyBboxPatch((x-0.8, y-0.3), 1.6, 0.6,
                                           boxstyle="round,pad=0.05",
                                           facecolor=style['color'], ec='black',
                                           linewidth=2, zorder=10)
            ax_main.add_patch(rect)

        # Add label
        if name in ['Act', 'Red', 'CDA', 'Cpk']:
            ax_main.text(x, y, name, ha='center', va='center', fontsize=11,
                        fontweight='bold', style='italic')
        elif name == 'SC_RS17645':
            ax_main.text(x, y+0.7, 'N-6 MTase', ha='center', fontsize=8, fontweight='bold')
            ax_main.text(x, y-0.7, '(SC_RS17645)', ha='center', fontsize=7)
        else:
            ax_main.text(x, y-0.55, name, ha='center', fontsize=8, fontweight='bold')

            # Add expression direction arrow for key regulators
            tf_data = tf_methyl[tf_methyl['name'] == name]
            if len(tf_data) > 0:
                log2fc = tf_data.iloc[0].get('log2FC_T2_vs_T1')
                if pd.notna(log2fc) and abs(log2fc) > 1:
                    direction = '↑' if log2fc > 0 else '↓'
                    color = '#59a14f' if log2fc > 0 else '#e15759'
                    ax_main.text(x+0.4, y+0.15, direction, fontsize=12,
                                color=color, fontweight='bold')

    # Draw regulatory edges
    edges = [
        # MTase → methylation targets
        ('SC_RS17645', 'redZ', '#9467bd', '--', 'Methylates'),

        # Global → Pleiotropic
        ('afsR', 'afsS', 'green', '-', ''),
        ('bldD', 'absA2', 'red', '-', ''),

        # Global → CSR
        ('bldD', 'actII-ORF4', 'gray', '-', ''),
        ('bldD', 'redD', 'gray', '-', ''),
        ('afsR', 'actII-ORF4', 'green', '-', ''),

        # Pleiotropic → CSR
        ('absA2', 'actII-ORF4', 'red', '-', ''),
        ('absA2', 'redD', 'red', '-', ''),
        ('absA2', 'cdaR', 'red', '-', ''),
        ('absA2', 'cpkO/kasO', 'red', '-', ''),
        ('hrdD', 'actII-ORF4', 'green', '-', ''),
        ('hrdD', 'redD', 'green', '-', ''),

        # CSR → CSR (cascade)
        ('redZ', 'redD', 'green', '-', ''),

        # CSR → BGC
        ('actII-ORF4', 'Act', 'green', '-', ''),
        ('redD', 'Red', 'green', '-', ''),
        ('cdaR', 'CDA', 'green', '-', ''),
        ('cpkO/kasO', 'Cpk', 'green', '-', ''),
    ]

    for src, tgt, color, style, label in edges:
        if src in positions and tgt in positions:
            x1, y1 = positions[src]
            x2, y2 = positions[tgt]

            # Calculate offset for source and target
            dx = x2 - x1
            dy = y2 - y1
            dist = np.sqrt(dx**2 + dy**2)
            if dist > 0:
                # Offset start and end points
                src_offset = 0.4 if src != 'SC_RS17645' else 0.5
                tgt_offset = 0.4 if tgt not in ['Act', 'Red', 'CDA', 'Cpk'] else 0.3

                x1_adj = x1 + (dx/dist) * src_offset
                y1_adj = y1 + (dy/dist) * src_offset
                x2_adj = x2 - (dx/dist) * tgt_offset
                y2_adj = y2 - (dy/dist) * tgt_offset

                linestyle = '--' if style == '--' else '-'
                arrow = mpatches.FancyArrowPatch(
                    (x1_adj, y1_adj), (x2_adj, y2_adj),
                    connectionstyle="arc3,rad=0.1",
                    arrowstyle='-|>',
                    mutation_scale=12,
                    color=color,
                    linewidth=1.5,
                    linestyle=linestyle,
                    alpha=0.7,
                    zorder=5
                )
                ax_main.add_patch(arrow)

                # Add label for special edges
                if label:
                    mid_x = (x1 + x2) / 2
                    mid_y = (y1 + y2) / 2
                    ax_main.text(mid_x+0.3, mid_y, label, fontsize=7, style='italic')

    # Set axis limits and properties
    ax_main.set_xlim(0, 21)
    ax_main.set_ylim(1, 9.5)
    ax_main.set_aspect('equal')
    ax_main.axis('off')
    ax_main.set_title('Epigenetic Regulatory Network: MTase → TF → BGC',
                     fontsize=14, fontweight='bold', y=0.98)

    # Add legend panel
    ax_legend = fig.add_axes([0.72, 0.55, 0.26, 0.40])
    ax_legend.axis('off')
    ax_legend.set_title('Legend', fontsize=12, fontweight='bold', loc='left')

    legend_items = [
        ('Methylation Status (T2 vs T1)', None, None),
        ('', plt.Circle((0, 0), 0.1, color='#e15759'), 'Lost'),
        ('', plt.Circle((0, 0), 0.1, color='#59a14f'), 'Gained'),
        ('', plt.Circle((0, 0), 0.1, color='#76b7b2'), 'Stable'),
        ('', plt.Circle((0, 0), 0.1, color='#e0e0e0'), 'Unmethylated'),
        ('', plt.Circle((0, 0), 0.1, color='#9467bd'), 'MTase enzyme'),
        ('', None, ''),
        ('Regulatory Interaction', None, None),
        ('', Line2D([0], [0], color='green', lw=2), 'Activation'),
        ('', Line2D([0], [0], color='red', lw=2), 'Repression'),
        ('', Line2D([0], [0], color='gray', lw=2), 'Regulation'),
        ('', Line2D([0], [0], color='#9467bd', lw=2, linestyle='--'), 'Methylation'),
    ]

    y_pos = 0.95
    for label, handle, text in legend_items:
        if handle is None and text is None:
            ax_legend.text(0.05, y_pos, label, fontsize=10, fontweight='bold', transform=ax_legend.transAxes)
            y_pos -= 0.08
        elif handle is None:
            y_pos -= 0.05
        else:
            ax_legend.text(0.15, y_pos, text, fontsize=9, transform=ax_legend.transAxes, va='center')
            y_pos -= 0.07

    # Add key insight panel
    ax_insight = fig.add_axes([0.72, 0.25, 0.26, 0.28])
    ax_insight.axis('off')
    ax_insight.set_title('Key Insights', fontsize=12, fontweight='bold', loc='left')

    insights_text = """
1. redZ is the ONLY GRN TF with
   coordinated methylation-expression
   (6mA lost → expression DOWN)

2. SC_RS17645 (N-6 MTase)
   expression DOWN in T2
   → Reduced AAGCCCG methylation

3. PARADOX: redD UP despite redZ DOWN
   → AbsA2 repression relief or
   other activators compensate
"""
    ax_insight.text(0.05, 0.9, insights_text, fontsize=9, transform=ax_insight.transAxes,
                   va='top', family='monospace')

    # Add summary statistics panel
    ax_stats = fig.add_axes([0.05, 0.05, 0.90, 0.18])
    ax_stats.axis('off')

    # Create summary table
    stats_data = [
        ['Category', 'Total', 'Methylated', 'Key Findings'],
        ['Tier 1 (Global)', '14', '3', 'bldN lost, afsR/afsS stable'],
        ['Tier 2 (Pleiotropic)', '17', '0', 'absA2 unmethylated but UP (6.3x)'],
        ['Tier 3 (CSR)', '6', '1', 'redZ coordinated (6mA lost + DOWN)'],
        ['BGC Expression', '-', '-', 'Act: mixed, Red: 9/10 UP, CDA: mixed, Cpk: mixed'],
    ]

    table = ax_stats.table(cellText=stats_data, loc='center', cellLoc='left',
                          colWidths=[0.15, 0.08, 0.10, 0.50])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.2, 1.5)

    # Style header row
    for i in range(4):
        table[(0, i)].set_facecolor('#4472C4')
        table[(0, i)].set_text_props(color='white', fontweight='bold')

    plt.savefig(OUTPUT_DIR / 'comprehensive_epigenetic_network.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: comprehensive_epigenetic_network.png")

def create_cascade_diagram():
    """Create detailed cascade diagram for the paper"""
    fig, ax = plt.subplots(figsize=(14, 10))

    # Define cascade steps
    steps = [
        {'y': 9, 'label': 'Environmental Signal\n(T2 time point)', 'color': '#e0e0e0'},
        {'y': 7.5, 'label': 'SC_RS17645 (N-6 MTase)\nExpression DOWN (log2FC = -2.19)', 'color': '#9467bd'},
        {'y': 6, 'label': 'AAGCCCG (6mA) Methylation\nReduced', 'color': '#f28e2b'},
        {'y': 4.5, 'label': 'redZ Promoter\n6mA Demethylated', 'color': '#4e79a7'},
        {'y': 3, 'label': 'redZ Expression\nDOWN (log2FC = -2.25)', 'color': '#e15759'},
        {'y': 1.5, 'label': 'redD + Red BGC\n(Complex regulation)', 'color': '#59a14f'},
    ]

    # Draw boxes and arrows
    box_width = 5
    box_height = 1

    for i, step in enumerate(steps):
        # Draw box
        rect = mpatches.FancyBboxPatch((5-box_width/2, step['y']-box_height/2),
                                       box_width, box_height,
                                       boxstyle="round,pad=0.05",
                                       facecolor=step['color'], ec='black',
                                       linewidth=2, alpha=0.8)
        ax.add_patch(rect)

        # Add label
        ax.text(5, step['y'], step['label'], ha='center', va='center',
               fontsize=11, fontweight='bold')

        # Draw arrow to next step
        if i < len(steps) - 1:
            next_y = steps[i+1]['y']
            ax.annotate('', xy=(5, next_y + box_height/2 + 0.1),
                       xytext=(5, step['y'] - box_height/2 - 0.1),
                       arrowprops=dict(arrowstyle='-|>', color='black', lw=2))

    # Add side annotations
    annotations = [
        (10, 7.5, 'REBASE: AAGCCCG is novel\n(not registered)', '#9467bd'),
        (10, 6, 'MEME analysis:\n75.6% 4mC at CCGG\nAAGCCCG for 6mA', '#f28e2b'),
        (10, 4.5, 'Coordinated change:\nMethylation + Expression', '#4e79a7'),
        (10, 3, 'Only GRN TF with\nthis coordination', '#e15759'),
        (10, 1.5, 'Paradox:\nredD UP (+4.77)\ndespite redZ DOWN', '#59a14f'),
    ]

    for x, y, text, color in annotations:
        ax.text(x, y, text, fontsize=9, va='center',
               bbox=dict(boxstyle='round', facecolor='white', edgecolor=color, alpha=0.8))
        ax.annotate('', xy=(7.7, y), xytext=(x-0.2, y),
                   arrowprops=dict(arrowstyle='-', color=color, lw=1, linestyle='--'))

    # Add title
    ax.set_title('Proposed Epigenetic Cascade for Red BGC Regulation\n'
                'MTase → Methylation → TF → Target',
                fontsize=14, fontweight='bold')

    ax.set_xlim(0, 14)
    ax.set_ylim(0, 10)
    ax.axis('off')

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'epigenetic_cascade_diagram.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: epigenetic_cascade_diagram.png")

def create_paper_figure_summary():
    """Create a summary figure suitable for publication"""
    fig = plt.figure(figsize=(18, 14))

    # Panel A: Network overview
    ax_a = fig.add_subplot(2, 2, 1)
    ax_a.set_title('A. Regulatory Network Overview', fontsize=12, fontweight='bold', loc='left')

    # Simple network diagram
    positions = {
        'MTase': (0.5, 0.9),
        'Global TF': (0.3, 0.6),
        'CSR': (0.7, 0.6),
        'redZ': (0.5, 0.4),
        'BGC': (0.5, 0.15),
    }

    colors = {
        'MTase': '#9467bd',
        'Global TF': '#e0e0e0',
        'CSR': '#e0e0e0',
        'redZ': '#e15759',
        'BGC': '#f28e2b',
    }

    for name, (x, y) in positions.items():
        circle = plt.Circle((x, y), 0.08, color=colors[name], ec='black', lw=2)
        ax_a.add_patch(circle)
        ax_a.text(x, y-0.15, name, ha='center', fontsize=9, fontweight='bold')

    # Arrows
    arrows = [('MTase', 'redZ'), ('Global TF', 'CSR'), ('CSR', 'BGC'), ('redZ', 'BGC')]
    for src, tgt in arrows:
        x1, y1 = positions[src]
        x2, y2 = positions[tgt]
        ax_a.annotate('', xy=(x2, y2+0.08), xytext=(x1, y1-0.08),
                     arrowprops=dict(arrowstyle='->', color='gray', lw=1.5))

    ax_a.set_xlim(0, 1)
    ax_a.set_ylim(0, 1)
    ax_a.axis('off')

    # Panel B: Methylation heatmap
    ax_b = fig.add_subplot(2, 2, 2)
    ax_b.set_title('B. TF Methylation Status', fontsize=12, fontweight='bold', loc='left')

    tf_methyl = pd.read_csv(ANALYSIS_DIR / 'tf_promoter_methylation/tf_methylation_summary.csv')

    # Prepare data for heatmap
    key_tfs = ['bldD', 'bldN', 'afsR', 'afsS', 'absA2', 'redZ', 'redD', 'actII-ORF4']
    heatmap_data = []
    for tf in key_tfs:
        tf_data = tf_methyl[tf_methyl['name'] == tf]
        if len(tf_data) > 0:
            row = tf_data.iloc[0]
            heatmap_data.append([
                row['total_T1_sites'],
                row['total_T2_sites'],
                row['total_T3_sites']
            ])
        else:
            heatmap_data.append([0, 0, 0])

    heatmap_df = pd.DataFrame(heatmap_data, index=key_tfs, columns=['T1', 'T2', 'T3'])
    im = ax_b.imshow(heatmap_df.values, cmap='Blues', aspect='auto')
    ax_b.set_xticks(range(3))
    ax_b.set_xticklabels(['T1', 'T2', 'T3'])
    ax_b.set_yticks(range(len(key_tfs)))
    ax_b.set_yticklabels(key_tfs)
    ax_b.set_xlabel('Time point')
    ax_b.set_ylabel('Transcription Factor')
    plt.colorbar(im, ax=ax_b, label='Methylation sites')

    # Panel C: Expression changes
    ax_c = fig.add_subplot(2, 2, 3)
    ax_c.set_title('C. Key Regulator Expression (T2 vs T1)', fontsize=12, fontweight='bold', loc='left')

    expr_data = []
    for tf in key_tfs:
        tf_data = tf_methyl[tf_methyl['name'] == tf]
        if len(tf_data) > 0 and pd.notna(tf_data.iloc[0]['log2FC_T2_vs_T1']):
            expr_data.append({'TF': tf, 'log2FC': tf_data.iloc[0]['log2FC_T2_vs_T1']})

    expr_df = pd.DataFrame(expr_data).sort_values('log2FC')
    colors = ['#59a14f' if x > 0 else '#e15759' for x in expr_df['log2FC']]
    ax_c.barh(expr_df['TF'], expr_df['log2FC'], color=colors)
    ax_c.axvline(x=0, color='black', lw=0.5)
    ax_c.set_xlabel('log2 Fold Change')

    # Panel D: BGC correlation
    ax_d = fig.add_subplot(2, 2, 4)
    ax_d.set_title('D. TF-BGC Expression Correlation', fontsize=12, fontweight='bold', loc='left')

    corr_data = tf_methyl[['name', 'corr_act', 'corr_red']].dropna(subset=['corr_act'])
    ax_d.scatter(corr_data['corr_act'], corr_data['corr_red'], s=80, alpha=0.7)
    for _, row in corr_data.iterrows():
        if row['name'] in ['actII-ORF4', 'afsR', 'absA2', 'sigU']:
            ax_d.annotate(row['name'], (row['corr_act'], row['corr_red']), fontsize=8)

    ax_d.axhline(y=0, color='gray', ls='--', alpha=0.5)
    ax_d.axvline(x=0, color='gray', ls='--', alpha=0.5)
    ax_d.set_xlabel('Correlation with Act BGC')
    ax_d.set_ylabel('Correlation with Red BGC')
    ax_d.set_xlim(-1, 1)
    ax_d.set_ylim(-1, 1)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'paper_figure_summary.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: paper_figure_summary.png")

def main():
    print("="*70)
    print("EPIGENETIC REGULATORY NETWORK VISUALIZATION")
    print("Step C-4: Comprehensive Network Figures")
    print("="*70)

    # Load data
    tf_methyl, reg_interactions, triads, mtase = load_all_data()
    print(f"Loaded {len(tf_methyl)} TFs, {len(reg_interactions)} interactions")

    # Create visualizations
    print("\nGenerating figures...")
    create_comprehensive_network(tf_methyl, reg_interactions)
    create_cascade_diagram()
    create_paper_figure_summary()

    print("\n" + "="*60)
    print("OUTPUT FILES")
    print("="*60)
    print(f"Directory: {OUTPUT_DIR}")
    for f in sorted(OUTPUT_DIR.glob('*.png')):
        print(f"  - {f.name}")

    print("\n" + "="*60)
    print("FIGURE DESCRIPTIONS FOR REPORT")
    print("="*60)

    descriptions = """
1. comprehensive_epigenetic_network.png
   - Full regulatory network from MTase to BGC
   - Shows methylation status by color (lost=red, gained=green, stable=teal)
   - Includes expression direction (↑/↓) for key TFs
   - Summary statistics table at bottom

2. epigenetic_cascade_diagram.png
   - Linear cascade: MTase → Methylation → TF → Target
   - Side annotations with supporting evidence
   - Highlights the redZ → redD → Red pathway

3. paper_figure_summary.png
   - Four-panel summary for publication
   - A: Network overview, B: Methylation heatmap
   - C: Expression changes, D: TF-BGC correlations
"""
    print(descriptions)

if __name__ == '__main__':
    main()
