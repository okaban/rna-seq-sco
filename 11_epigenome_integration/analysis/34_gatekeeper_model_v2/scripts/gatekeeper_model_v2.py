#!/usr/bin/env python3
"""
H11: Gatekeeper Model v2 - Quantitative Epigenome-Transcriptome Integration
S. coelicolor M145 RNA-seq + Epigenome Integration Study

Synthesizes results from 9 hypotheses (H1-H9) into a 3-layer quantitative model.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, ArrowStyle
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd
import os

# Paths
BASE = '/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/34_gatekeeper_model_v2'
FIG_DIR = os.path.join(BASE, 'figures')
TBL_DIR = os.path.join(BASE, 'tables')

os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TBL_DIR, exist_ok=True)

# ============================================================
# Color palette
# ============================================================
COL_L1 = '#2E86AB'   # Layer 1: Landscape Remodeling (blue)
COL_L2 = '#A23B72'   # Layer 2: Protection/Depletion (magenta)
COL_L3 = '#F18F01'   # Layer 3: Signal Gating (orange)
COL_EXCLUDED = '#E63946'  # Excluded pathway (red)
COL_BG = '#F8F9FA'
COL_DARK = '#2D3436'
COL_SUPPORT = '#27AE60'
COL_REJECT = '#E74C3C'
COL_PARTIAL = '#F39C12'

# ============================================================
# FIGURE 1: Gatekeeper Model v2 Conceptual Diagram
# ============================================================
def create_model_diagram():
    fig, ax = plt.subplots(figsize=(14, 10))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 10)
    ax.set_aspect('equal')
    ax.axis('off')
    fig.patch.set_facecolor('white')

    # Title
    ax.text(7, 9.6, 'Gatekeeper Model v2: Epigenome-Transcriptome Integration',
            ha='center', va='center', fontsize=16, fontweight='bold', color=COL_DARK,
            fontstyle='italic')
    ax.text(7, 9.2, 'S. coelicolor A3(2) M145 | 3 Timepoints | 9 Hypotheses Integrated',
            ha='center', va='center', fontsize=10, color='#636E72')

    # --- Layer 1: Landscape Remodeling ---
    l1_box = FancyBboxPatch((0.5, 6.8), 13, 2.1,
                             boxstyle="round,pad=0.1",
                             facecolor=COL_L1, alpha=0.12, edgecolor=COL_L1, linewidth=2.5)
    ax.add_patch(l1_box)
    ax.text(0.8, 8.65, 'LAYER 1', fontsize=9, fontweight='bold', color=COL_L1, va='center')
    ax.text(2.4, 8.65, 'LANDSCAPE REMODELING', fontsize=12, fontweight='bold',
            color=COL_L1, va='center')

    # Layer 1 content boxes
    items_l1 = [
        ('MTase expression\nchange', 'SC_RS17645 \u2193\n\u2192 260 sites lost', '7.0'),
        ('CCGG/GCCGGC\nturnover', 'Jaccard = 0.000\nComplete positional shift', '10.0'),
        ('4mC-expression\ncorrelation', '\u03c1 = 0.139\nperm p = 0.002', '4.0'),
    ]
    for i, (title, detail, x_pos) in enumerate(items_l1):
        x = float(x_pos)
        box = FancyBboxPatch((x-1.3, 7.0), 2.6, 1.45,
                              boxstyle="round,pad=0.08",
                              facecolor='white', edgecolor=COL_L1, linewidth=1.5, alpha=0.95)
        ax.add_patch(box)
        ax.text(x, 8.15, title, ha='center', va='center', fontsize=8.5, fontweight='bold',
                color=COL_L1)
        ax.text(x, 7.45, detail, ha='center', va='center', fontsize=7.5, color=COL_DARK)

    # Arrows L1 -> L2
    for x_pos in [4.0, 7.0, 10.0]:
        ax.annotate('', xy=(x_pos, 6.65), xytext=(x_pos, 6.85),
                     arrowprops=dict(arrowstyle='->', color=COL_DARK, lw=1.5))

    # --- Layer 2: Protection / Depletion ---
    l2_box = FancyBboxPatch((0.5, 4.3), 13, 2.1,
                             boxstyle="round,pad=0.1",
                             facecolor=COL_L2, alpha=0.12, edgecolor=COL_L2, linewidth=2.5)
    ax.add_patch(l2_box)
    ax.text(0.8, 6.15, 'LAYER 2', fontsize=9, fontweight='bold', color=COL_L2, va='center')
    ax.text(2.4, 6.15, 'PROTECTION / DEPLETION', fontsize=12, fontweight='bold',
            color=COL_L2, va='center')

    items_l2 = [
        ('TF binding sites\ndepleted', 'fold = 0.66\np = 6.8\u00d710\u207b\u2075', '3.0'),
        ('\u03c3 factor -10 box\ndepleted', 'p = 1.73\u00d710\u207b\xb9\xb2', '6.0'),
        ('Literature 37 TFs\nunmethylated', '92% methyl-free\nSARP 100%', '9.0'),
        ('DMGs: no\nfunctional enrich.', 'Position-dependent\nFunction-independent', '12.0'),
    ]
    for i, (title, detail, x_pos) in enumerate(items_l2):
        x = float(x_pos)
        box = FancyBboxPatch((x-1.3, 4.5), 2.6, 1.45,
                              boxstyle="round,pad=0.08",
                              facecolor='white', edgecolor=COL_L2, linewidth=1.5, alpha=0.95)
        ax.add_patch(box)
        ax.text(x, 5.65, title, ha='center', va='center', fontsize=8.5, fontweight='bold',
                color=COL_L2)
        ax.text(x, 4.95, detail, ha='center', va='center', fontsize=7.5, color=COL_DARK)

    # Arrows L2 -> L3
    for x_pos in [4.5, 7.5, 10.5]:
        ax.annotate('', xy=(x_pos, 4.15), xytext=(x_pos, 4.35),
                     arrowprops=dict(arrowstyle='->', color=COL_DARK, lw=1.5))

    # --- Layer 3: Signal Gating (Revised) ---
    l3_box = FancyBboxPatch((0.5, 1.6), 13, 2.3,
                             boxstyle="round,pad=0.1",
                             facecolor=COL_L3, alpha=0.12, edgecolor=COL_L3, linewidth=2.5)
    ax.add_patch(l3_box)
    ax.text(0.8, 3.65, 'LAYER 3', fontsize=9, fontweight='bold', color=COL_L3, va='center')
    ax.text(2.4, 3.65, 'SIGNAL GATING (Revised)', fontsize=12, fontweight='bold',
            color=COL_L3, va='center')

    items_l3 = [
        ('57 non-literature\nregulators', 'All outside 37 list\n15.6% of 1,055 methylated', '2.5'),
        ('TCS asymmetric\nmethylation', '7 pairs identified\nSK methylated / RR free', '5.5'),
        ('T3 arm enrichment\n(gain + up)', 'OR = 8.05, p = 0.001\n87% at arms', '8.5'),
        ('COG T enriched\n(signal transd.)', 'OR = 2.47, p = 0.045\nCOG Q = 0', '11.5'),
    ]
    for i, (title, detail, x_pos) in enumerate(items_l3):
        x = float(x_pos)
        box = FancyBboxPatch((x-1.3, 1.8), 2.6, 1.55,
                              boxstyle="round,pad=0.08",
                              facecolor='white', edgecolor=COL_L3, linewidth=1.5, alpha=0.95)
        ax.add_patch(box)
        ax.text(x, 3.05, title, ha='center', va='center', fontsize=8.5, fontweight='bold',
                color=COL_L3)
        ax.text(x, 2.3, detail, ha='center', va='center', fontsize=7.5, color=COL_DARK)

    # --- Excluded Pathway (bottom bar) ---
    excl_box = FancyBboxPatch((0.5, 0.2), 13, 1.1,
                               boxstyle="round,pad=0.1",
                               facecolor=COL_EXCLUDED, alpha=0.08, edgecolor=COL_EXCLUDED,
                               linewidth=2.0, linestyle='--')
    ax.add_patch(excl_box)
    ax.text(7, 1.0, 'EXCLUDED: Methylation \u2192 TF cascade \u2192 BGC activation',
            ha='center', va='center', fontsize=11, fontweight='bold', color=COL_EXCLUDED)
    ax.text(7, 0.55, 'H4: 0/37 literature TFs show methylation-expression coordination  |  '
            'SARP activators (redD, actII-ORF4, cdaR, cpkO) ALL unmethylated',
            ha='center', va='center', fontsize=7.5, color=COL_EXCLUDED, alpha=0.85)

    # Big X over excluded
    ax.plot([0.8, 2.0], [0.35, 1.15], color=COL_EXCLUDED, lw=3, alpha=0.4)
    ax.plot([0.8, 2.0], [1.15, 0.35], color=COL_EXCLUDED, lw=3, alpha=0.4)

    # Flow arrows on the side
    ax.annotate('', xy=(0.2, 4.4), xytext=(0.2, 8.8),
                 arrowprops=dict(arrowstyle='->', color=COL_DARK, lw=2.5,
                                  connectionstyle='arc3,rad=0'))
    ax.text(0.15, 6.6, 'Flow', rotation=90, ha='center', va='center',
            fontsize=9, fontweight='bold', color=COL_DARK)

    plt.tight_layout()
    for fmt in ['pdf', 'svg']:
        fig.savefig(os.path.join(FIG_DIR, f'gatekeeper_model_v2.{fmt}'),
                    dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print('  [OK] Figure 1: Gatekeeper Model v2 diagram')


# ============================================================
# FIGURE 2: Evidence Integration Matrix (H1-H9 x Layers)
# ============================================================
def create_evidence_matrix():
    hypotheses = ['H1', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9']
    h_labels = [
        'H1: 4mC/6mA distribution\n& expression',
        'H3: redZ paradox\nresolution',
        'H4: TF methylation-\nexpression cascade',
        'H5: MTase expression\n\u2192 site dynamics',
        'H6: Genome-wide\nregulator screen',
        'H7: CCGG MTase\nparadox',
        'H8: Coordinated reg.\ncharacterization',
        'H9: CCGG 5mC\nmisclassification',
    ]
    layers = ['Layer 1:\nLandscape\nRemodeling',
              'Layer 2:\nProtection/\nDepletion',
              'Layer 3:\nSignal\nGating',
              'Excluded:\nTF Cascade']

    verdicts = ['Partial', 'Rejected\n(premise)', 'Rejected', 'Partial',
                'Supported', 'Rejected\n(passive)', 'Partial', 'Rejected']
    verdict_colors = [COL_PARTIAL, COL_SUPPORT, COL_REJECT, COL_PARTIAL,
                      COL_SUPPORT, COL_REJECT, COL_PARTIAL, COL_REJECT]

    # Matrix: rows=hypotheses, cols=layers
    # Values: 0=no link, 1=supports, 2=strong support, -1=rejects
    matrix = np.array([
        [2,  0,  0,  0],   # H1: strong support L1
        [0,  0,  0, -1],   # H3: rejects excluded (premise error)
        [0,  1,  0, -1],   # H4: supports L2, rejects excluded
        [2,  0,  0,  0],   # H5: strong support L1
        [0,  0,  2,  0],   # H6: strong support L3
        [2,  0,  1,  0],   # H7: strong support L1, supports L3
        [0,  0,  2,  0],   # H8: strong support L3
        [1,  0,  0,  0],   # H9: supports L1 (4mC identity confirmed)
    ], dtype=float)

    fig, (ax_main, ax_verdict) = plt.subplots(1, 2, figsize=(12, 6),
                                               gridspec_kw={'width_ratios': [3, 1]})
    fig.patch.set_facecolor('white')

    # Color map
    cmap_colors = ['#E74C3C', '#FFFFFF', '#F0E68C', '#27AE60', '#1B7A3D']
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list('evidence',
           [(0, '#E74C3C'), (0.25, '#F5B7B1'), (0.5, '#F8F9FA'),
            (0.75, '#ABEBC6'), (1.0, '#27AE60')])

    im = ax_main.imshow(matrix, cmap=cmap, aspect='auto', vmin=-1.5, vmax=2.5,
                         interpolation='nearest')

    ax_main.set_xticks(range(4))
    ax_main.set_xticklabels(layers, fontsize=9, fontweight='bold')
    ax_main.set_yticks(range(8))
    ax_main.set_yticklabels(h_labels, fontsize=8.5)

    # Annotate cells
    labels_map = {-1: 'Rejects', 0: '--', 1: 'Supports', 2: 'Strong\nSupport'}
    for i in range(8):
        for j in range(4):
            val = matrix[i, j]
            txt = labels_map.get(val, '')
            color = 'white' if abs(val) >= 1.5 else COL_DARK
            ax_main.text(j, i, txt, ha='center', va='center', fontsize=7.5,
                         fontweight='bold' if val != 0 else 'normal', color=color)

    # Grid
    for i in range(9):
        ax_main.axhline(i - 0.5, color='white', lw=2)
    for j in range(5):
        ax_main.axvline(j - 0.5, color='white', lw=2)

    ax_main.set_title('Evidence Integration Matrix\n(Hypotheses H1-H9 \u00d7 Model Layers)',
                       fontsize=13, fontweight='bold', pad=15)

    # Verdict column
    ax_verdict.set_xlim(0, 2)
    ax_verdict.set_ylim(-0.5, 7.5)
    ax_verdict.invert_yaxis()
    ax_verdict.axis('off')
    ax_verdict.set_title('Verdict', fontsize=12, fontweight='bold', pad=15)

    for i, (v, c) in enumerate(zip(verdicts, verdict_colors)):
        box = FancyBboxPatch((0.1, i - 0.35), 1.8, 0.7,
                              boxstyle="round,pad=0.05",
                              facecolor=c, alpha=0.2, edgecolor=c, linewidth=1.5)
        ax_verdict.add_patch(box)
        ax_verdict.text(1.0, i, v, ha='center', va='center', fontsize=8.5,
                         fontweight='bold', color=c)

    plt.tight_layout()
    for fmt in ['pdf', 'svg']:
        fig.savefig(os.path.join(FIG_DIR, f'evidence_integration_matrix.{fmt}'),
                    dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print('  [OK] Figure 2: Evidence Integration Matrix')


# ============================================================
# FIGURE 3: Layer-Specific Quantitative Summaries
# ============================================================
def create_layer_summaries():
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.patch.set_facecolor('white')
    fig.suptitle('Gatekeeper Model v2: Quantitative Parameters by Layer',
                 fontsize=14, fontweight='bold', y=0.98)

    # --- Panel A: Layer 1 - Key statistics ---
    ax = axes[0, 0]
    params = ['4mC-expr\ncorrelation\n(\u03c1)', 'Permutation\np-value',
              'CCGG site\noverlap\n(Jaccard)', 'AAGCCCG\nsites lost\n(n)']
    values = [0.139, 0.002, 0.000, 260]
    colors = [COL_L1] * 4
    # Normalize for display
    display_vals = [0.139, 0.002, 0.000, 260]

    ax.barh(range(4), [0.139, -np.log10(0.002), 0.001, 260/260], color=COL_L1, alpha=0.7,
            edgecolor=COL_L1, linewidth=1.5, height=0.6)
    for i, (p, v) in enumerate(zip(params, values)):
        ax.text(0.02, i, f'{p}: {v}', va='center', fontsize=9, fontweight='bold', color='white',
                path_effects=[pe.withStroke(linewidth=2, foreground=COL_L1)])
    ax.set_xlim(-0.1, 1.2)
    ax.set_yticks([])
    ax.set_title('Layer 1: Landscape Remodeling', fontsize=11, fontweight='bold', color=COL_L1)
    ax.set_xlabel('Normalized Scale', fontsize=9)

    # Redraw as a proper table-like display
    ax.clear()
    ax.axis('off')
    ax.set_title('Layer 1: Landscape Remodeling', fontsize=12, fontweight='bold',
                  color=COL_L1, pad=10)
    table_data_l1 = [
        ['4mC promoter-expression correlation', '\u03c1 = 0.139, p = 8.26\u00d710\u207b\u2074'],
        ['Permutation test (4mC T2vsT1)', 'r = 0.125, p = 0.002'],
        ['CCGG 4mC positional overlap', 'Jaccard = 0.000'],
        ['AAGCCCG site loss (SC_RS17645\u2193)', '260 sites, 0% retention'],
        ['T1 CCGG core genome fraction', '77%'],
        ['T2 CCGG chromosomal arms', '89%'],
        ['6mA T2vsT1 correlation', 'Lost after T3 resequencing'],
    ]
    table = ax.table(cellText=table_data_l1,
                      colLabels=['Parameter', 'Value'],
                      loc='center', cellLoc='left',
                      colColours=[COL_L1 + '33', COL_L1 + '33'])
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1.0, 1.4)
    for key, cell in table.get_celld().items():
        if key[0] == 0:
            cell.set_text_props(fontweight='bold', color=COL_L1)
        cell.set_edgecolor('#BDC3C7')

    # --- Panel B: Layer 2 - Protection ---
    ax = axes[0, 1]
    ax.axis('off')
    ax.set_title('Layer 2: Protection / Depletion', fontsize=12, fontweight='bold',
                  color=COL_L2, pad=10)
    table_data_l2 = [
        ['TF binding site methylation', 'fold = 0.66, p = 6.8\u00d710\u207b\u2075'],
        ['\u03c3-10 box methylation depletion', 'p = 1.73\u00d710\u207b\xb9\xb2'],
        ['Literature 37 TFs methylated', '3/37 (8%), all constitutive'],
        ['SARP activators methylated', '0/4 (0%)'],
        ['SARP activator LFC range', '+2.0 to +7.2 (unmethylated)'],
        ['DMG functional enrichment', 'NONE (GO/KEGG/COG)'],
        ['TF BS methylation: Lost T2T3', '75%'],
    ]
    table = ax.table(cellText=table_data_l2,
                      colLabels=['Parameter', 'Value'],
                      loc='center', cellLoc='left',
                      colColours=[COL_L2 + '33', COL_L2 + '33'])
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1.0, 1.4)
    for key, cell in table.get_celld().items():
        if key[0] == 0:
            cell.set_text_props(fontweight='bold', color=COL_L2)
        cell.set_edgecolor('#BDC3C7')

    # --- Panel C: Layer 3 - Signal Gating ---
    ax = axes[1, 0]
    ax.axis('off')
    ax.set_title('Layer 3: Signal Gating (Revised)', fontsize=12, fontweight='bold',
                  color=COL_L3, pad=10)
    table_data_l3 = [
        ['Total regulatory genes screened', '1,055 (25 families)'],
        ['Methylated regulatory genes', '165/1,055 (15.6%)'],
        ['Coordinated methyl+expr genes', '57 (all outside literature 37)'],
        ['T3 discordant_gain_up at arms', '87% (OR=8.05, p=0.001)'],
        ['COG T (signal transduction)', 'OR=2.47, p=0.045'],
        ['COG Q (secondary metabolism)', '0 genes'],
        ['TCS asymmetric pairs', '7 identified'],
        ['Top candidate: SC_RS10435', 'Chaplin-adjacent TetR'],
        ['Top candidate: SC_RS35525', 'PPTase-adjacent SK'],
    ]
    table = ax.table(cellText=table_data_l3,
                      colLabels=['Parameter', 'Value'],
                      loc='center', cellLoc='left',
                      colColours=[COL_L3 + '33', COL_L3 + '33'])
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1.0, 1.3)
    for key, cell in table.get_celld().items():
        if key[0] == 0:
            cell.set_text_props(fontweight='bold', color=COL_L3)
        cell.set_edgecolor('#BDC3C7')

    # --- Panel D: Family methylation rates ---
    ax = axes[1, 1]
    # Read family data
    fam_df = pd.read_csv('/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/'
                          '29_genomewide_TF_screen/tables/regulatory_methylation_summary.tsv',
                          sep='\t')
    fam_df = fam_df.sort_values('methylation_rate', ascending=True)

    colors_bar = []
    for _, row in fam_df.iterrows():
        if row['tf_family'] == 'SARP':
            colors_bar.append(COL_EXCLUDED)
        elif row['methylation_rate'] >= 20:
            colors_bar.append(COL_L3)
        elif row['methylation_rate'] > 0:
            colors_bar.append(COL_L1)
        else:
            colors_bar.append('#BDC3C7')

    ax.barh(range(len(fam_df)), fam_df['methylation_rate'], color=colors_bar, alpha=0.8,
            edgecolor='white', linewidth=0.5, height=0.7)
    ax.set_yticks(range(len(fam_df)))
    ax.set_yticklabels([f"{row['tf_family']} (n={row['n_genes']})"
                         for _, row in fam_df.iterrows()], fontsize=7.5)
    ax.set_xlabel('Methylation Rate (%)', fontsize=10)
    ax.set_title('Regulatory Family Methylation Rates', fontsize=12, fontweight='bold',
                  color=COL_DARK, pad=10)
    ax.axvline(x=0, color='black', linewidth=0.5)

    # Highlight SARP
    sarp_idx = list(fam_df['tf_family']).index('SARP')
    ax.get_yticklabels()[sarp_idx].set_color(COL_EXCLUDED)
    ax.get_yticklabels()[sarp_idx].set_fontweight('bold')

    # Add MerR label
    merr_idx = list(fam_df['tf_family']).index('MerR')
    ax.get_yticklabels()[merr_idx].set_fontweight('bold')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    for fmt in ['pdf', 'svg']:
        fig.savefig(os.path.join(FIG_DIR, f'layer_quantitative_summaries.{fmt}'),
                    dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print('  [OK] Figure 3: Layer Quantitative Summaries')


# ============================================================
# FIGURE 4: Excluded Pathway Diagram
# ============================================================
def create_excluded_pathway():
    fig, ax = plt.subplots(figsize=(14, 6))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 6)
    ax.set_aspect('equal')
    ax.axis('off')
    fig.patch.set_facecolor('white')

    ax.text(7, 5.7, 'Excluded Pathway: Methylation-Dependent TF Cascade',
            ha='center', va='center', fontsize=14, fontweight='bold', color=COL_EXCLUDED)

    # --- Top row: REJECTED pathway ---
    ax.text(0.5, 4.8, 'REJECTED (H4)', fontsize=11, fontweight='bold', color=COL_EXCLUDED)

    boxes_top = [
        (1.5, 3.8, 'DNA\nMethylation', COL_EXCLUDED),
        (4.5, 3.8, 'Literature 37\nTFs Activation', COL_EXCLUDED),
        (7.5, 3.8, 'SARP\nActivators', COL_EXCLUDED),
        (10.5, 3.8, 'BGC\nExpression', COL_EXCLUDED),
    ]
    for x, y, label, col in boxes_top:
        box = FancyBboxPatch((x-1.0, y), 2.0, 1.0,
                              boxstyle="round,pad=0.08",
                              facecolor=col, alpha=0.1, edgecolor=col,
                              linewidth=2, linestyle='--')
        ax.add_patch(box)
        ax.text(x, y+0.5, label, ha='center', va='center', fontsize=9,
                fontweight='bold', color=col, alpha=0.6)

    # Crossed arrows
    for x1, x2 in [(2.5, 3.5), (5.5, 6.5), (8.5, 9.5)]:
        ax.annotate('', xy=(x2, 4.3), xytext=(x1, 4.3),
                     arrowprops=dict(arrowstyle='->', color=COL_EXCLUDED, lw=2,
                                      linestyle='--', alpha=0.4))
        # X marks
        mid = (x1 + x2) / 2
        ax.plot([mid-0.15, mid+0.15], [4.1, 4.5], color=COL_EXCLUDED, lw=3, alpha=0.7)
        ax.plot([mid-0.15, mid+0.15], [4.5, 4.1], color=COL_EXCLUDED, lw=3, alpha=0.7)

    # --- Bottom row: ACTUAL pathway ---
    ax.text(0.5, 2.6, 'ACTUAL (H6, H8)', fontsize=11, fontweight='bold', color=COL_SUPPORT)

    boxes_bot = [
        (1.5, 1.4, 'Landscape\nRemodeling', COL_L1),
        (4.5, 1.4, 'Protected\nZones (TF BS)', COL_L2),
        (7.5, 1.4, '57 Non-Lit.\nRegulators', COL_L3),
        (10.5, 1.4, 'Signal\nTransduction', COL_L3),
    ]
    for x, y, label, col in boxes_bot:
        box = FancyBboxPatch((x-1.0, y), 2.0, 1.0,
                              boxstyle="round,pad=0.08",
                              facecolor=col, alpha=0.2, edgecolor=col, linewidth=2)
        ax.add_patch(box)
        ax.text(x, y+0.5, label, ha='center', va='center', fontsize=9,
                fontweight='bold', color=col)

    for x1, x2 in [(2.5, 3.5), (5.5, 6.5), (8.5, 9.5)]:
        ax.annotate('', xy=(x2, 1.9), xytext=(x1, 1.9),
                     arrowprops=dict(arrowstyle='->', color=COL_SUPPORT, lw=2.5))

    # Evidence annotations
    evidence_items = [
        (3.0, 3.2, '3/37 methylated\n(8%, all constitutive)', COL_EXCLUDED),
        (6.0, 3.2, 'redD +5.79, actII-ORF4\ncpkO, cdaR: 0 methylation', COL_EXCLUDED),
        (9.0, 3.2, 'SARP 0% methylation\nrate (n=7)', COL_EXCLUDED),

        (3.0, 0.9, 'TF BS fold=0.66\n\u03c3-10 depleted', COL_L2),
        (6.0, 0.9, '57/57 outside\nliterature 37', COL_L3),
        (9.0, 0.9, 'COG T: OR=2.47\nTCS asymmetric', COL_L3),
    ]
    for x, y, txt, col in evidence_items:
        ax.text(x, y, txt, ha='center', va='center', fontsize=7, color=col,
                fontstyle='italic')

    # Output box
    out_box = FancyBboxPatch((12.0, 1.4), 1.5, 1.0,
                              boxstyle="round,pad=0.08",
                              facecolor=COL_SUPPORT, alpha=0.15, edgecolor=COL_SUPPORT, linewidth=2)
    ax.add_patch(out_box)
    ax.text(12.75, 1.9, 'Secondary\nMetabolism\n(Indirect)', ha='center', va='center',
            fontsize=8.5, fontweight='bold', color=COL_SUPPORT)
    ax.annotate('', xy=(12.0, 1.9), xytext=(11.5, 1.9),
                 arrowprops=dict(arrowstyle='->', color=COL_SUPPORT, lw=2.5))

    plt.tight_layout()
    for fmt in ['pdf', 'svg']:
        fig.savefig(os.path.join(FIG_DIR, f'excluded_pathway_diagram.{fmt}'),
                    dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print('  [OK] Figure 4: Excluded Pathway Diagram')


# ============================================================
# TABLE 1: H1-H9 Evidence Integration
# ============================================================
def create_evidence_table():
    data = [
        {
            'Hypothesis': 'H1',
            'Description': '4mC and 6mA have different genomic distributions and expression correlations',
            'Verdict': 'Partial support',
            'Key_Statistic': '4mC/6mA distribution p=1.56e-14; 4mC promoter rho=0.139 p=8.26e-4; 6mA: no correlation',
            'Layer_Supported': 'Layer 1',
            'Layer_Impact': 'Strong support - 4mC is the dominant modification linked to expression ("quality not quantity")',
            'Evidence_Strength': 'Strong',
        },
        {
            'Hypothesis': 'H3',
            'Description': 'redZ paradox (apparent downregulation despite Red pathway activation)',
            'Verdict': 'Premise rejected',
            'Key_Statistic': 'True redZ (SC_RS31650) LFC=+0.91; redD LFC=+5.79; locus_tag error 36/37',
            'Layer_Supported': 'Excluded pathway',
            'Layer_Impact': 'Literature locus_tag error discovered; paradox does not exist; redD is main driver',
            'Evidence_Strength': 'Definitive',
        },
        {
            'Hypothesis': 'H4',
            'Description': '37 literature TFs have methylation-expression coordination',
            'Verdict': 'Rejected',
            'Key_Statistic': '3/37 (8%) methylated, all constitutive; SARP 0/4; 0 coordinated',
            'Layer_Supported': 'Layer 2 + Excluded pathway',
            'Layer_Impact': 'TF network is methylation-INDEPENDENT; Excluded pathway definitively rejected',
            'Evidence_Strength': 'Definitive',
        },
        {
            'Hypothesis': 'H5',
            'Description': 'MTase expression changes drive methylation site dynamics',
            'Verdict': 'Partial support',
            'Key_Statistic': 'SC_RS17645 down -> 260 AAGCCCG sites lost; gene-level |rho|<0.16 for all motifs',
            'Layer_Supported': 'Layer 1',
            'Layer_Impact': 'Strong support - methylation operates at global/wholesale level not gene-specific cis',
            'Evidence_Strength': 'Strong (site level) / Weak (gene level)',
        },
        {
            'Hypothesis': 'H6',
            'Description': 'Genome-wide regulatory gene methylation screen (1055 genes)',
            'Verdict': 'Supported',
            'Key_Statistic': '165/1055 (15.6%) methylated; 57 coordinated; ALL outside literature 37; MerR 28.6%, SARP 0%',
            'Layer_Supported': 'Layer 3',
            'Layer_Impact': 'Core evidence for revised Signal Gating layer; non-literature regulators are the targets',
            'Evidence_Strength': 'Strong',
        },
        {
            'Hypothesis': 'H7',
            'Description': 'CCGG 4mC passive dilution vs active remodeling',
            'Verdict': 'Passive dilution rejected',
            'Key_Statistic': 'Jaccard=0.000; T1 core 77% -> T2 arms 89%; defense island co-induction T3',
            'Layer_Supported': 'Layer 1 + Layer 3',
            'Layer_Impact': 'Active geographic remodeling; CCGG system has defense function',
            'Evidence_Strength': 'Definitive',
        },
        {
            'Hypothesis': 'H8',
            'Description': '57 coordinated regulators: geographic and functional characterization',
            'Verdict': 'Partial support',
            'Key_Statistic': 'Overall arm p=0.17; T3 gain_up arms 87% OR=8.05 p=0.001; COG T OR=2.47 p=0.045',
            'Layer_Supported': 'Layer 3',
            'Layer_Impact': 'Asymmetric Signal Gating model; TCS pairs show one-sided methylation',
            'Evidence_Strength': 'Moderate-Strong',
        },
        {
            'Hypothesis': 'H9',
            'Description': 'CCGG "4mC" is 5mC misclassification?',
            'Verdict': 'Rejected',
            'Key_Statistic': '5mC signal=0.01%; 4mC freq=82.7%; 0% CCWGG; 91.4% GCCGGC palindrome',
            'Layer_Supported': 'Layer 1',
            'Layer_Impact': 'Confirms 4mC identity; GCCGGC is true recognition motif for novel N4-C MTase',
            'Evidence_Strength': 'Definitive',
        },
    ]

    df = pd.DataFrame(data)
    df.to_csv(os.path.join(TBL_DIR, 'H1_H9_evidence_integration.tsv'), sep='\t', index=False)
    print('  [OK] Table 1: H1-H9 Evidence Integration')
    return df


# ============================================================
# TABLE 2: Layer Parameters
# ============================================================
def create_layer_params_table():
    data = [
        # Layer 1
        {'Layer': 'Layer 1: Landscape Remodeling', 'Parameter': '4mC promoter-expression correlation',
         'Value': 'rho=0.139', 'P_value': '8.26e-4', 'Source_Hypothesis': 'H1',
         'Interpretation': '4mC is the dominant modification linked to expression'},
        {'Layer': 'Layer 1: Landscape Remodeling', 'Parameter': '4mC T2vsT1 robust correlation',
         'Value': 'r=0.125', 'P_value': '0.002 (perm)', 'Source_Hypothesis': 'Pre-loop',
         'Interpretation': 'Confirmed by permutation test'},
        {'Layer': 'Layer 1: Landscape Remodeling', 'Parameter': 'CCGG site positional overlap',
         'Value': 'Jaccard=0.000', 'P_value': 'N/A', 'Source_Hypothesis': 'H7',
         'Interpretation': 'Complete positional turnover between timepoints'},
        {'Layer': 'Layer 1: Landscape Remodeling', 'Parameter': 'AAGCCCG site loss',
         'Value': '260 sites', 'P_value': 'N/A', 'Source_Hypothesis': 'H5',
         'Interpretation': 'MTase downregulation leads to wholesale site loss'},
        {'Layer': 'Layer 1: Landscape Remodeling', 'Parameter': '6mA T2vsT1 correlation',
         'Value': 'Lost', 'P_value': 'NS after T3 reseq', 'Source_Hypothesis': 'Pre-loop',
         'Interpretation': '6mA correlation was artifactual (coverage bias)'},
        {'Layer': 'Layer 1: Landscape Remodeling', 'Parameter': 'T1 CCGG core genome fraction',
         'Value': '77%', 'P_value': 'N/A', 'Source_Hypothesis': 'H7',
         'Interpretation': 'CCGG sites concentrated in core at T1'},
        {'Layer': 'Layer 1: Landscape Remodeling', 'Parameter': 'T2 CCGG arms fraction',
         'Value': '89%', 'P_value': 'N/A', 'Source_Hypothesis': 'H7',
         'Interpretation': 'Geographic shift to chromosomal arms at T2'},
        {'Layer': 'Layer 1: Landscape Remodeling', 'Parameter': '5mC misclassification excluded',
         'Value': '5mC=0.01%', 'P_value': 'N/A', 'Source_Hypothesis': 'H9',
         'Interpretation': 'CCGG 4mC is genuine N4-cytosine modification'},

        # Layer 2
        {'Layer': 'Layer 2: Protection/Depletion', 'Parameter': 'TF binding site methylation fold',
         'Value': '0.66', 'P_value': '6.8e-5', 'Source_Hypothesis': 'Pre-loop',
         'Interpretation': 'Methylation depleted at TF binding sites'},
        {'Layer': 'Layer 2: Protection/Depletion', 'Parameter': 'Sigma -10 box depletion',
         'Value': 'Depleted', 'P_value': '1.73e-12', 'Source_Hypothesis': 'Pre-loop',
         'Interpretation': 'Methylation actively excluded from promoter elements'},
        {'Layer': 'Layer 2: Protection/Depletion', 'Parameter': 'Literature TFs methylation rate',
         'Value': '8% (3/37)', 'P_value': 'N/A', 'Source_Hypothesis': 'H4',
         'Interpretation': 'Key TFs are in the protected zone'},
        {'Layer': 'Layer 2: Protection/Depletion', 'Parameter': 'SARP activator methylation',
         'Value': '0% (0/4)', 'P_value': 'N/A', 'Source_Hypothesis': 'H4',
         'Interpretation': 'Pathway-specific activators completely unmethylated'},
        {'Layer': 'Layer 2: Protection/Depletion', 'Parameter': 'DMG functional enrichment',
         'Value': 'NONE', 'P_value': 'All NS', 'Source_Hypothesis': 'Pre-loop',
         'Interpretation': 'Methylation is position-dependent, function-independent'},
        {'Layer': 'Layer 2: Protection/Depletion', 'Parameter': 'TF BS methylation lost T2T3',
         'Value': '75%', 'P_value': 'N/A', 'Source_Hypothesis': 'Pre-loop',
         'Interpretation': 'Methylation at regulatory regions is transient'},

        # Layer 3
        {'Layer': 'Layer 3: Signal Gating', 'Parameter': 'Total regulatory genes screened',
         'Value': '1,055', 'P_value': 'N/A', 'Source_Hypothesis': 'H6',
         'Interpretation': '25 families from GFF annotation'},
        {'Layer': 'Layer 3: Signal Gating', 'Parameter': 'Methylated regulatory genes',
         'Value': '165 (15.6%)', 'P_value': 'N/A', 'Source_Hypothesis': 'H6',
         'Interpretation': '2.3x higher than literature 37 rate (6.7%)'},
        {'Layer': 'Layer 3: Signal Gating', 'Parameter': 'Coordinated regulators',
         'Value': '57', 'P_value': 'N/A', 'Source_Hypothesis': 'H6',
         'Interpretation': 'All outside literature 37 list'},
        {'Layer': 'Layer 3: Signal Gating', 'Parameter': 'T3 discordant_gain_up arm enrichment',
         'Value': '87% at arms', 'P_value': '0.001 (OR=8.05)', 'Source_Hypothesis': 'H8',
         'Interpretation': 'Geographic specificity at late timepoint'},
        {'Layer': 'Layer 3: Signal Gating', 'Parameter': 'COG T enrichment',
         'Value': 'OR=2.47', 'P_value': '0.045', 'Source_Hypothesis': 'H8',
         'Interpretation': 'Signal transduction enriched'},
        {'Layer': 'Layer 3: Signal Gating', 'Parameter': 'TCS asymmetric pairs',
         'Value': '7', 'P_value': 'N/A', 'Source_Hypothesis': 'H8',
         'Interpretation': 'Sensor kinase methylated, response regulator free'},
        {'Layer': 'Layer 3: Signal Gating', 'Parameter': 'Highest methylation family',
         'Value': 'MerR 28.6%', 'P_value': 'N/A', 'Source_Hypothesis': 'H6',
         'Interpretation': 'Metal-sensing regulators most targeted'},
        {'Layer': 'Layer 3: Signal Gating', 'Parameter': 'SARP methylation rate',
         'Value': '0%', 'P_value': 'N/A', 'Source_Hypothesis': 'H6',
         'Interpretation': 'Pathway-specific activators completely excluded'},
    ]

    df = pd.DataFrame(data)
    df.to_csv(os.path.join(TBL_DIR, 'layer_parameters.tsv'), sep='\t', index=False)
    print('  [OK] Table 2: Layer Parameters')
    return df


# ============================================================
# TABLE 3: Testable Predictions
# ============================================================
def create_predictions_table():
    data = [
        {
            'Prediction_ID': 'P1',
            'Layer': 'Layer 1',
            'Prediction': 'SC_RS17645 knockout will eliminate all AAGCCCG methylation sites genome-wide',
            'Test_Method': 'Gene knockout + PacBio/Nanopore sequencing',
            'Expected_Outcome': 'Complete loss of AAGCCCG methylation; no compensatory MTase activity',
            'Priority': 'High',
            'Derived_From': 'H5',
        },
        {
            'Prediction_ID': 'P2',
            'Layer': 'Layer 1',
            'Prediction': 'GCCGGC-recognizing MTase exists but is not SC_RS36410 alone; a separate unidentified enzyme creates T1 core CCGG sites',
            'Test_Method': 'SC_RS36410 knockout + methylome profiling; search for GCCGGC MTase in un-annotated ORFs',
            'Expected_Outcome': 'SC_RS36410 KO eliminates T2/T3 arm CCGG sites but T1-type core sites persist or are lost by different MTase',
            'Priority': 'High',
            'Derived_From': 'H7',
        },
        {
            'Prediction_ID': 'P3',
            'Layer': 'Layer 2',
            'Prediction': 'Artificial methylation of SARP activator promoters (redD, actII-ORF4) will NOT reduce BGC expression',
            'Test_Method': 'Engineered MTase targeting SARP promoters + RT-qPCR',
            'Expected_Outcome': 'No significant change in BGC expression (methylation is inert at these loci)',
            'Priority': 'Medium',
            'Derived_From': 'H4',
        },
        {
            'Prediction_ID': 'P4',
            'Layer': 'Layer 3',
            'Prediction': 'SC_RS10435 (chaplin-adjacent TetR) methylation loss causes derepression of chaplin/rodlin morphogenesis genes',
            'Test_Method': 'SC_RS10435 KO + RNA-seq during aerial mycelium development',
            'Expected_Outcome': 'Chaplin genes upregulated; premature aerial mycelium formation',
            'Priority': 'High',
            'Derived_From': 'H6, H8',
        },
        {
            'Prediction_ID': 'P5',
            'Layer': 'Layer 3',
            'Prediction': 'SC_RS35525 (PPTase-adjacent SK) demethylation activates its kinase activity, phosphorylating its cognate RR',
            'Test_Method': 'In vitro phosphorylation assay +/- methylated vs unmethylated SC_RS35525 promoter',
            'Expected_Outcome': 'Higher SK expression leads to RR phosphorylation and downstream activation',
            'Priority': 'High',
            'Derived_From': 'H6, H8',
        },
        {
            'Prediction_ID': 'P6',
            'Layer': 'Layer 3',
            'Prediction': 'TCS asymmetric methylation pattern will be conserved across Streptomyces species with similar R-M systems',
            'Test_Method': 'Comparative methylome analysis of close Streptomyces relatives',
            'Expected_Outcome': 'SK preferentially methylated vs RR in species with same MTase families',
            'Priority': 'Medium',
            'Derived_From': 'H8',
        },
        {
            'Prediction_ID': 'P7',
            'Layer': 'Layer 1',
            'Prediction': 'T1->T2 CCGG geographic shift (core->arms) is reproducible across biological replicates',
            'Test_Method': 'Repeat PacBio sequencing of independent M145 cultures at equivalent growth stages',
            'Expected_Outcome': 'Consistent core-to-arm shift pattern; stochastic at individual site level but deterministic at regional level',
            'Priority': 'High',
            'Derived_From': 'H7',
        },
        {
            'Prediction_ID': 'P8',
            'Layer': 'Layer 2',
            'Prediction': 'TF binding site methylation depletion is active (not passive): sites near but outside BS are methylated at normal rates',
            'Test_Method': 'High-resolution methylome analysis at known BS positions +/- flanking regions',
            'Expected_Outcome': 'Sharp boundary of methylation depletion at BS edges; flanking regions show normal methylation',
            'Priority': 'Medium',
            'Derived_From': 'H4, Pre-loop',
        },
        {
            'Prediction_ID': 'P9',
            'Layer': 'Layer 3',
            'Prediction': 'MerR family regulators (28.6% methylated) modulate metal-responsive gene expression via methylation gating',
            'Test_Method': 'Metal stress experiments + methylome profiling',
            'Expected_Outcome': 'MerR methylation changes in response to metal stress; correlated with metal-responsive gene expression',
            'Priority': 'Medium',
            'Derived_From': 'H6',
        },
        {
            'Prediction_ID': 'P10',
            'Layer': 'All layers',
            'Prediction': 'The 3-layer model applies to other Streptomyces species with complex R-M systems',
            'Test_Method': 'Multi-species epigenome-transcriptome integration study',
            'Expected_Outcome': 'Landscape remodeling + protection zones + signal gating observed in >= 2 additional species',
            'Priority': 'Low',
            'Derived_From': 'All hypotheses',
        },
    ]

    df = pd.DataFrame(data)
    df.to_csv(os.path.join(TBL_DIR, 'testable_predictions.tsv'), sep='\t', index=False)
    print('  [OK] Table 3: Testable Predictions')
    return df


# ============================================================
# FIGURE 5: Comprehensive 4-panel Figure (publication quality)
# ============================================================
def create_comprehensive_figure():
    fig = plt.figure(figsize=(16, 18))
    fig.patch.set_facecolor('white')

    # Use gridspec for fine control
    from matplotlib.gridspec import GridSpec
    gs = GridSpec(4, 2, figure=fig, hspace=0.35, wspace=0.3,
                  height_ratios=[1.2, 1, 1, 0.8])

    # --- Panel A: Model Overview (top, spans both columns) ---
    ax_a = fig.add_subplot(gs[0, :])
    ax_a.set_xlim(0, 16)
    ax_a.set_ylim(0, 4)
    ax_a.axis('off')
    ax_a.text(0.0, 3.8, 'A', fontsize=18, fontweight='bold', transform=ax_a.transAxes,
              va='top')

    # Three layer boxes
    layer_info = [
        (2.5, COL_L1, 'Layer 1:\nLandscape\nRemodeling',
         'Global site turnover\nJaccard=0.000\n4mC \u03c1=0.139'),
        (7.0, COL_L2, 'Layer 2:\nProtection /\nDepletion',
         'TF BS fold=0.66\n\u03c3-10 p=1.73e-12\n37 TFs: 92% free'),
        (11.5, COL_L3, 'Layer 3:\nSignal\nGating',
         '57 regulators\nTCS asymmetric\nCOG T OR=2.47'),
    ]

    for x, col, title, details in layer_info:
        box = FancyBboxPatch((x-1.8, 0.5), 3.6, 3.0,
                              boxstyle="round,pad=0.15",
                              facecolor=col, alpha=0.15, edgecolor=col, linewidth=2.5)
        ax_a.add_patch(box)
        ax_a.text(x, 3.0, title, ha='center', va='center', fontsize=12,
                  fontweight='bold', color=col)
        ax_a.text(x, 1.5, details, ha='center', va='center', fontsize=9, color=COL_DARK)

    # Arrows between layers
    ax_a.annotate('', xy=(5.2, 2.0), xytext=(4.3, 2.0),
                   arrowprops=dict(arrowstyle='->', color=COL_DARK, lw=2.5))
    ax_a.annotate('', xy=(9.7, 2.0), xytext=(8.8, 2.0),
                   arrowprops=dict(arrowstyle='->', color=COL_DARK, lw=2.5))

    # Excluded pathway below
    excl_box = FancyBboxPatch((5.2, -0.1), 5.6, 0.5,
                               boxstyle="round,pad=0.05",
                               facecolor=COL_EXCLUDED, alpha=0.1, edgecolor=COL_EXCLUDED,
                               linewidth=1.5, linestyle='--')
    ax_a.add_patch(excl_box)
    ax_a.text(8.0, 0.15, 'EXCLUDED: Methylation \u2192 TF cascade (H4 rejected)',
              ha='center', va='center', fontsize=9, color=COL_EXCLUDED, fontweight='bold')
    ax_a.plot([5.5, 5.9], [-0.05, 0.35], color=COL_EXCLUDED, lw=2, alpha=0.5)
    ax_a.plot([5.5, 5.9], [0.35, -0.05], color=COL_EXCLUDED, lw=2, alpha=0.5)

    # --- Panel B: Hypothesis verdict pie ---
    ax_b = fig.add_subplot(gs[1, 0])
    ax_b.text(-0.1, 1.05, 'B', fontsize=18, fontweight='bold', transform=ax_b.transAxes,
              va='top')

    verdict_counts = {'Supported': 1, 'Partial Support': 3, 'Rejected\n(hypothesis)': 2,
                      'Rejected\n(premise)': 2}
    colors_pie = [COL_SUPPORT, COL_PARTIAL, COL_REJECT, '#C0392B']
    wedges, texts, autotexts = ax_b.pie(
        verdict_counts.values(), labels=verdict_counts.keys(),
        colors=colors_pie, autopct='%1.0f%%', startangle=90,
        textprops={'fontsize': 9}, pctdistance=0.75)
    for at in autotexts:
        at.set_fontweight('bold')
        at.set_fontsize(10)
    ax_b.set_title('Hypothesis Verdicts (n=8)', fontsize=12, fontweight='bold')

    # --- Panel C: Family methylation rates bar chart ---
    ax_c = fig.add_subplot(gs[1, 1])
    ax_c.text(-0.1, 1.05, 'C', fontsize=18, fontweight='bold', transform=ax_c.transAxes,
              va='top')

    fam_df = pd.read_csv('/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/'
                          '29_genomewide_TF_screen/tables/regulatory_methylation_summary.tsv',
                          sep='\t')
    # Top 10 by n_genes
    fam_top = fam_df.nlargest(10, 'n_genes').sort_values('methylation_rate', ascending=True)

    bars_colors = []
    for _, row in fam_top.iterrows():
        if row['tf_family'] == 'SARP':
            bars_colors.append(COL_EXCLUDED)
        elif row['methylation_rate'] >= 20:
            bars_colors.append(COL_L3)
        else:
            bars_colors.append(COL_L1)

    ax_c.barh(range(len(fam_top)), fam_top['methylation_rate'],
              color=bars_colors, alpha=0.8, edgecolor='white', height=0.7)
    ax_c.set_yticks(range(len(fam_top)))
    ax_c.set_yticklabels([f"{r['tf_family']} (n={r['n_genes']})"
                           for _, r in fam_top.iterrows()], fontsize=9)
    ax_c.set_xlabel('Methylation Rate (%)', fontsize=10)
    ax_c.set_title('Regulatory Family Methylation Rates\n(Top 10 by gene count)',
                    fontsize=12, fontweight='bold')
    ax_c.spines['top'].set_visible(False)
    ax_c.spines['right'].set_visible(False)

    # --- Panel D: Coordination types bar ---
    ax_d = fig.add_subplot(gs[2, 0])
    ax_d.text(-0.1, 1.05, 'D', fontsize=18, fontweight='bold', transform=ax_d.transAxes,
              va='top')

    # Count coordination types from 57 genes
    coord_df = pd.read_csv('/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/'
                            '29_genomewide_TF_screen/tables/coordinated_regulatory_genes.tsv',
                            sep='\t')

    # T2 and T3 coordination types
    coord_types_t2 = coord_df['coordination_T2'].value_counts()
    coord_types_t3 = coord_df['coordination_T3'].value_counts()

    all_types = ['concordant_derepression', 'concordant_repression',
                 'discordant_gain_up', 'discordant_loss_down',
                 'methyl_change_no_expr_change', 'ambiguous']
    type_labels = ['Concordant\nDerepression', 'Concordant\nRepression',
                   'Discordant\nGain+Up', 'Discordant\nLoss+Down',
                   'Methyl only\n(no expr change)', 'Ambiguous']

    t2_vals = [coord_types_t2.get(t, 0) for t in all_types]
    t3_vals = [coord_types_t3.get(t, 0) for t in all_types]

    x = np.arange(len(all_types))
    width = 0.35
    ax_d.bar(x - width/2, t2_vals, width, label='T2 vs T1', color=COL_L1, alpha=0.8)
    ax_d.bar(x + width/2, t3_vals, width, label='T3 vs T1', color=COL_L3, alpha=0.8)
    ax_d.set_xticks(x)
    ax_d.set_xticklabels(type_labels, fontsize=7.5, rotation=0, ha='center')
    ax_d.set_ylabel('Number of Genes', fontsize=10)
    ax_d.set_title('Coordination Types in 57 Regulators', fontsize=12, fontweight='bold')
    ax_d.legend(fontsize=9)
    ax_d.spines['top'].set_visible(False)
    ax_d.spines['right'].set_visible(False)

    # --- Panel E: Geographic distribution ---
    ax_e = fig.add_subplot(gs[2, 1])
    ax_e.text(-0.1, 1.05, 'E', fontsize=18, fontweight='bold', transform=ax_e.transAxes,
              va='top')

    groups = ['Genome-wide', 'All regulators\n(1,055)', 'Coordinated 57']
    arm_pct = [44.7, 45.0, 51.6]
    core_pct = [55.3, 55.0, 48.4]

    x = np.arange(len(groups))
    ax_e.bar(x, arm_pct, 0.5, label='Chromosomal Arms', color=COL_L3, alpha=0.8)
    ax_e.bar(x, core_pct, 0.5, bottom=arm_pct, label='Core', color=COL_L1, alpha=0.8)
    ax_e.set_xticks(x)
    ax_e.set_xticklabels(groups, fontsize=9)
    ax_e.set_ylabel('Percentage (%)', fontsize=10)
    ax_e.set_title('Geographic Distribution\n(Core vs Arms)', fontsize=12, fontweight='bold')
    ax_e.legend(fontsize=9, loc='upper right')
    ax_e.set_ylim(0, 105)
    ax_e.axhline(y=50, color='gray', linestyle='--', alpha=0.5, lw=0.8)
    ax_e.spines['top'].set_visible(False)
    ax_e.spines['right'].set_visible(False)

    # Annotation
    ax_e.text(2, 25, 'p=0.17 (NS)\nOverall', ha='center', fontsize=9, fontweight='bold',
              color='white')
    ax_e.text(2, 75, 'T3 gain+up:\n87% arms\np=0.001', ha='center', fontsize=8,
              fontweight='bold', color='white', fontstyle='italic')

    # --- Panel F: Key numbers summary ---
    ax_f = fig.add_subplot(gs[3, :])
    ax_f.axis('off')
    ax_f.text(0.0, 1.05, 'F', fontsize=18, fontweight='bold', transform=ax_f.transAxes,
              va='top')

    summary_text = (
        'Model Summary  |  '
        'Layer 1: 4mC dominates (not 6mA), complete positional turnover (Jaccard=0.000), wholesale MTase-driven dynamics  |  '
        'Layer 2: TF BS depleted (fold=0.66), \u03c3-10 depleted (p=1.73e-12), 92% literature TFs methylation-free  |  '
        'Layer 3: 57 non-literature regulators, TCS asymmetric methylation, COG T enriched (OR=2.47), arms-enriched at T3  |  '
        'EXCLUDED: Direct methylation-TF cascade-BGC pathway (0/37 coordination, SARP 0% methylated)'
    )
    ax_f.text(0.5, 0.5, summary_text, ha='center', va='center', fontsize=9,
              color=COL_DARK, wrap=True,
              bbox=dict(boxstyle='round,pad=0.5', facecolor='#F8F9FA', edgecolor='#BDC3C7'))

    plt.savefig(os.path.join(FIG_DIR, 'comprehensive_model_figure.pdf'),
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(os.path.join(FIG_DIR, 'comprehensive_model_figure.svg'),
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print('  [OK] Figure 5: Comprehensive Model Figure')


# ============================================================
# Main execution
# ============================================================
if __name__ == '__main__':
    print('='*60)
    print('H11: Gatekeeper Model v2 - Figure & Table Generation')
    print('='*60)
    print()

    print('Generating figures...')
    create_model_diagram()
    create_evidence_matrix()
    create_layer_summaries()
    create_excluded_pathway()
    create_comprehensive_figure()

    print()
    print('Generating tables...')
    ev_df = create_evidence_table()
    lp_df = create_layer_params_table()
    pred_df = create_predictions_table()

    print()
    print('='*60)
    print('All outputs generated successfully.')
    print(f'Figures: {FIG_DIR}/')
    print(f'Tables:  {TBL_DIR}/')
    print('='*60)
