#!/usr/bin/env python3
"""
Methylation Motif Position Schematic
Creates schematic diagrams showing:
- Gene structure with TSS, promoter elements (-10/-35)
- Methylation site positions relative to TSS
- Methylation frequency at each site

Reference: Schmidt et al., 2024; Huang et al., 2025

Streptomyces coelicolor A3(2) M145 Project
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle, FancyArrowPatch, FancyBboxPatch
from matplotlib.collections import PatchCollection
import os

# Configuration
INTEGRATION_DIR = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis"
OUTPUT_DIR = os.path.join(INTEGRATION_DIR, "motif_schematics")
METHYL_PATH = os.path.join(INTEGRATION_DIR, "high_confidence_sites_weighted.csv")
GFF_PATH = "/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/genomic.gff"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Key genes with known methylation patterns
KEY_GENES = {
    'RamR': {
        'gene_id': 'SC_RS35610',
        'old_locus': 'SCO6685',
        'start': 7426396,
        'end': 7427004,
        'strand': '-',
        'function': 'Response regulator (SapB activator)',
        'pattern': '6mA loss → Activation',
        'mod_type': '6mA'
    },
    'NsdB': {
        'gene_id': 'SC_RS38475',
        'old_locus': 'SCO7252',
        'start': 8062135,
        'end': 8063643,
        'strand': '+',
        'function': 'DNA-binding protein (master regulator)',
        'pattern': '4mC gain → Strong activation',
        'mod_type': '4mC'
    },
    'Red_SCO5897': {
        'gene_id': 'SC_RS31730',
        'old_locus': 'SCO5897',
        'start': 6462296,
        'end': 6463483,
        'strand': '+',
        'function': 'Red cluster oxygenase',
        'pattern': '4mC gain → Upregulation',
        'mod_type': '4mC'
    },
    'Act_SCO5079': {
        'gene_id': 'SC_RS27555',
        'old_locus': 'SCO5079',
        'start': 5520857,
        'end': 5521741,
        'strand': '+',
        'function': 'Act cluster NmrA protein',
        'pattern': '4mC loss → Upregulation',
        'mod_type': '4mC'
    },
    'Cpk_SCO6284': {
        'gene_id': 'SC_RS33670',
        'old_locus': 'SCO6284',
        'start': 6943673,
        'end': 6945265,
        'strand': '+',
        'function': 'Cpk cluster carboxylase',
        'pattern': '6mA gain → Strong upregulation',
        'mod_type': '6mA'
    }
}

def load_methylation_data():
    """Load high-confidence methylation sites"""
    return pd.read_csv(METHYL_PATH)

def get_promoter_methylation_sites(methyl_df, gene_info):
    """Get methylation sites in promoter region (-300 to +50 from TSS)"""

    if gene_info['strand'] == '+':
        tss = gene_info['start']
        promoter_start = tss - 300
        promoter_end = tss + 50
    else:
        tss = gene_info['end']
        promoter_start = tss - 50
        promoter_end = tss + 300

    sites = methyl_df[
        (methyl_df['position'] >= min(promoter_start, promoter_end)) &
        (methyl_df['position'] <= max(promoter_start, promoter_end))
    ].copy()

    # Calculate position relative to TSS
    if gene_info['strand'] == '+':
        sites['rel_position'] = sites['position'] - tss
    else:
        sites['rel_position'] = tss - sites['position']

    return sites, tss

def create_gene_schematic(gene_name, gene_info, methyl_df):
    """Create schematic diagram for a single gene"""

    sites, tss = get_promoter_methylation_sites(methyl_df, gene_info)

    # Filter by modification type of interest
    sites_filtered = sites[sites['mod_type'] == gene_info['mod_type']]

    fig, ax = plt.subplots(figsize=(14, 8))

    # Y-axis positions
    y_gene = 0.3
    y_promoter = 0.5
    y_methyl = 0.7
    y_freq = 0.85

    # X-axis range: -350 to +100 from TSS
    x_min, x_max = -350, 150

    # Colors
    colors = {
        '6mA': '#E74C3C',
        '4mC': '#3498DB',
        'gene': '#2C3E50',
        'promoter': '#F39C12',
        'tss': '#27AE60'
    }

    # Draw scale bar
    ax.plot([x_min, x_max], [0.1, 0.1], 'k-', linewidth=1)
    for x in range(-300, 151, 50):
        ax.plot([x, x], [0.08, 0.12], 'k-', linewidth=1)
        ax.text(x, 0.03, str(x), ha='center', fontsize=9)

    ax.text(0, 0.0, 'Position relative to TSS (bp)', ha='center', fontsize=10)

    # Draw promoter elements (-35 and -10 boxes)
    # Typical bacterial promoter elements
    box_35 = Rectangle((-38, y_promoter - 0.03), 6, 0.06, facecolor='#9B59B6',
                       edgecolor='black', linewidth=1, alpha=0.8)
    box_10 = Rectangle((-13, y_promoter - 0.03), 6, 0.06, facecolor='#9B59B6',
                       edgecolor='black', linewidth=1, alpha=0.8)
    ax.add_patch(box_35)
    ax.add_patch(box_10)
    ax.text(-35, y_promoter + 0.06, '-35', ha='center', fontsize=9, fontweight='bold')
    ax.text(-10, y_promoter + 0.06, '-10', ha='center', fontsize=9, fontweight='bold')

    # Draw TSS arrow
    ax.annotate('', xy=(5, y_promoter), xytext=(0, y_promoter),
                arrowprops=dict(arrowstyle='->', color=colors['tss'], lw=2))
    ax.plot([0, 0], [y_promoter - 0.08, y_promoter + 0.08], color=colors['tss'],
            linewidth=2, linestyle='-')
    ax.text(0, y_promoter + 0.10, 'TSS', ha='center', fontsize=10, fontweight='bold',
            color=colors['tss'])

    # Draw gene body
    gene_arrow = FancyArrowPatch((10, y_gene), (100, y_gene),
                                  arrowstyle='->', mutation_scale=20,
                                  color=colors['gene'], linewidth=3)
    ax.add_patch(gene_arrow)
    ax.text(55, y_gene + 0.06, f'{gene_info["old_locus"]}', ha='center',
            fontsize=11, fontweight='bold', color=colors['gene'])

    # Draw promoter region box
    promoter_box = Rectangle((-300, y_promoter - 0.05), 350, 0.10,
                              facecolor=colors['promoter'], alpha=0.2,
                              edgecolor=colors['promoter'], linewidth=1, linestyle='--')
    ax.add_patch(promoter_box)
    ax.text(-150, y_promoter - 0.08, 'Promoter region (-300 to +50)',
            ha='center', fontsize=9, style='italic', color=colors['promoter'])

    # Draw methylation sites
    mod_color = colors[gene_info['mod_type']]

    # Get unique positions and their methylation frequencies by timepoint
    unique_positions = sites_filtered['rel_position'].unique()

    for pos in unique_positions:
        pos_data = sites_filtered[sites_filtered['rel_position'] == pos]

        # Draw methylation site marker
        ax.plot([pos, pos], [y_methyl - 0.02, y_methyl + 0.02],
                color=mod_color, linewidth=3)
        ax.scatter([pos], [y_methyl], c=mod_color, s=100, zorder=5, marker='o')

        # Add frequency values by timepoint
        for tp_idx, tp in enumerate(['T1', 'T2', 'T3']):
            tp_data = pos_data[pos_data['timepoint'] == tp]
            if len(tp_data) > 0:
                freq = tp_data['weighted_mod_freq'].values[0]
                y_offset = y_freq + tp_idx * 0.04
                color_tp = ['#2ECC71', '#F39C12', '#9B59B6'][tp_idx]
                ax.text(pos, y_offset, f'{tp}: {freq:.0f}%',
                       ha='center', fontsize=8, color=color_tp)

    # Draw label for methylation type
    ax.text(x_min + 20, y_methyl, f'{gene_info["mod_type"]} sites',
            ha='left', fontsize=10, fontweight='bold', color=mod_color)

    # Add legend
    legend_elements = [
        mpatches.Patch(facecolor='#9B59B6', label='-35/-10 boxes'),
        plt.Line2D([0], [0], color=colors['tss'], linewidth=2, label='TSS'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=mod_color,
                   markersize=10, label=f'{gene_info["mod_type"]} methylation'),
        mpatches.Patch(facecolor=colors['promoter'], alpha=0.3, label='Promoter region')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=9)

    # Set axis properties
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(-0.05, 1.0)
    ax.set_aspect('auto')
    ax.axis('off')

    # Title
    ax.set_title(f'{gene_name} ({gene_info["old_locus"]}) - Promoter Methylation Schematic\n'
                 f'{gene_info["function"]}\n'
                 f'Pattern: {gene_info["pattern"]}',
                 fontsize=12, fontweight='bold', pad=20)

    plt.tight_layout()

    # Save
    output_path = os.path.join(OUTPUT_DIR, f'motif_schematic_{gene_name}.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight', facecolor='white')
    plt.close()

    return output_path, len(unique_positions)

def create_combined_schematic(key_genes, methyl_df):
    """Create combined schematic showing all key genes"""

    n_genes = len(key_genes)
    fig, axes = plt.subplots(n_genes, 1, figsize=(14, 4*n_genes))

    if n_genes == 1:
        axes = [axes]

    colors = {
        '6mA': '#E74C3C',
        '4mC': '#3498DB',
        'gene': '#2C3E50',
        'promoter': '#F39C12',
        'tss': '#27AE60'
    }

    x_min, x_max = -350, 150

    for idx, (gene_name, gene_info) in enumerate(key_genes.items()):
        ax = axes[idx]

        sites, tss = get_promoter_methylation_sites(methyl_df, gene_info)
        sites_filtered = sites[sites['mod_type'] == gene_info['mod_type']]

        y_base = 0.2
        y_promoter = 0.4
        y_methyl = 0.6

        # Scale bar
        ax.plot([x_min, x_max], [0.05, 0.05], 'k-', linewidth=0.5)
        for x in range(-300, 151, 100):
            ax.plot([x, x], [0.03, 0.07], 'k-', linewidth=0.5)
            ax.text(x, 0.0, str(x), ha='center', fontsize=7)

        # Promoter boxes
        box_35 = Rectangle((-38, y_promoter - 0.05), 6, 0.1, facecolor='#9B59B6',
                           edgecolor='black', linewidth=0.5, alpha=0.7)
        box_10 = Rectangle((-13, y_promoter - 0.05), 6, 0.1, facecolor='#9B59B6',
                           edgecolor='black', linewidth=0.5, alpha=0.7)
        ax.add_patch(box_35)
        ax.add_patch(box_10)

        # TSS
        ax.plot([0, 0], [y_promoter - 0.1, y_promoter + 0.1], color=colors['tss'],
                linewidth=2)
        ax.annotate('', xy=(10, y_promoter), xytext=(0, y_promoter),
                    arrowprops=dict(arrowstyle='->', color=colors['tss'], lw=1.5))

        # Gene body
        gene_arrow = FancyArrowPatch((15, y_base), (100, y_base),
                                      arrowstyle='->', mutation_scale=15,
                                      color=colors['gene'], linewidth=2)
        ax.add_patch(gene_arrow)

        # Methylation sites
        mod_color = colors[gene_info['mod_type']]
        unique_positions = sites_filtered['rel_position'].unique()

        for pos in unique_positions:
            pos_data = sites_filtered[sites_filtered['rel_position'] == pos]

            ax.plot([pos, pos], [y_methyl - 0.05, y_methyl + 0.05],
                    color=mod_color, linewidth=2)
            ax.scatter([pos], [y_methyl], c=mod_color, s=60, zorder=5, marker='o')

            # Show T1 and T2 frequencies
            t1_data = pos_data[pos_data['timepoint'] == 'T1']
            t2_data = pos_data[pos_data['timepoint'] == 'T2']

            freq_t1 = t1_data['weighted_mod_freq'].values[0] if len(t1_data) > 0 else 0
            freq_t2 = t2_data['weighted_mod_freq'].values[0] if len(t2_data) > 0 else 0

            if freq_t1 > 0 or freq_t2 > 0:
                ax.text(pos, y_methyl + 0.12, f'T1:{freq_t1:.0f}%\nT2:{freq_t2:.0f}%',
                       ha='center', fontsize=7, linespacing=0.8)

        # Labels
        ax.text(x_min + 10, y_methyl, gene_info['mod_type'], ha='left', fontsize=9,
               fontweight='bold', color=mod_color)
        ax.text(120, y_base, gene_info['old_locus'], ha='left', fontsize=10,
               fontweight='bold', color=colors['gene'])

        ax.set_xlim(x_min, x_max)
        ax.set_ylim(-0.1, 0.85)
        ax.axis('off')

        ax.set_title(f'{gene_name}: {gene_info["pattern"]}', fontsize=10,
                    fontweight='bold', loc='left', pad=5)

    fig.suptitle('Methylation Motif Positions in Promoter Regions\n'
                 'Key Epigenetically Regulated Genes - S. coelicolor M145',
                 fontsize=13, fontweight='bold', y=1.01)

    plt.tight_layout()

    output_path = os.path.join(OUTPUT_DIR, 'motif_schematic_all_genes.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight', facecolor='white')
    plt.close()

    return output_path

def create_promoter_architecture_summary(key_genes, methyl_df):
    """Create summary figure showing promoter architecture with methylation"""

    fig, ax = plt.subplots(figsize=(16, 10))

    colors = {
        '6mA': '#E74C3C',
        '4mC': '#3498DB',
        'T1': '#2ECC71',
        'T2': '#F39C12',
        'T3': '#9B59B6'
    }

    y_spacing = 1.5
    x_min, x_max = -350, 200

    for idx, (gene_name, gene_info) in enumerate(key_genes.items()):
        y_base = (len(key_genes) - idx - 1) * y_spacing

        sites, tss = get_promoter_methylation_sites(methyl_df, gene_info)

        # Draw DNA backbone
        ax.plot([x_min, x_max], [y_base, y_base], 'k-', linewidth=2, alpha=0.5)

        # Draw promoter region highlight
        ax.axvspan(-300, 50, ymin=(y_base - 0.3) / (len(key_genes) * y_spacing),
                   ymax=(y_base + 0.3) / (len(key_genes) * y_spacing),
                   alpha=0.1, color='yellow')

        # Draw -35 and -10 boxes
        ax.add_patch(Rectangle((-38, y_base - 0.15), 6, 0.3, facecolor='#9B59B6',
                               edgecolor='black', linewidth=1))
        ax.add_patch(Rectangle((-13, y_base - 0.15), 6, 0.3, facecolor='#9B59B6',
                               edgecolor='black', linewidth=1))

        # Draw TSS
        ax.annotate('', xy=(15, y_base), xytext=(0, y_base),
                    arrowprops=dict(arrowstyle='->', color='green', lw=2))
        ax.plot([0, 0], [y_base - 0.25, y_base + 0.25], 'g-', linewidth=2)

        # Draw gene arrow
        ax.annotate('', xy=(180, y_base), xytext=(20, y_base),
                    arrowprops=dict(arrowstyle='->', color='#2C3E50', lw=3))

        # Draw methylation sites for both 6mA and 4mC
        for mod_type in ['6mA', '4mC']:
            mod_sites = sites[sites['mod_type'] == mod_type]
            unique_positions = mod_sites['rel_position'].unique()

            y_offset = 0.4 if mod_type == '6mA' else -0.4

            for pos in unique_positions:
                pos_data = mod_sites[mod_sites['rel_position'] == pos]

                # Get T1 and T2 frequencies
                t1_freq = pos_data[pos_data['timepoint'] == 'T1']['weighted_mod_freq'].values
                t2_freq = pos_data[pos_data['timepoint'] == 'T2']['weighted_mod_freq'].values

                t1_freq = t1_freq[0] if len(t1_freq) > 0 else 0
                t2_freq = t2_freq[0] if len(t2_freq) > 0 else 0

                # Draw site marker
                marker = 'o' if mod_type == '6mA' else 's'
                ax.scatter([pos], [y_base + y_offset], c=colors[mod_type], s=80,
                          marker=marker, zorder=5, edgecolors='black', linewidths=0.5)

                # Draw connection to backbone
                ax.plot([pos, pos], [y_base, y_base + y_offset], '--',
                       color=colors[mod_type], alpha=0.5, linewidth=1)

                # Add frequency change annotation
                if t1_freq > 0 or t2_freq > 0:
                    change = t2_freq - t1_freq
                    change_str = f'+{change:.0f}' if change > 0 else f'{change:.0f}'
                    change_color = colors['T2'] if change > 0 else colors['T3']
                    ax.text(pos, y_base + y_offset + 0.15 * (1 if mod_type == '6mA' else -1),
                           f'{change_str}%', ha='center', fontsize=7, color=change_color)

        # Gene name and pattern
        ax.text(x_min + 5, y_base + 0.6, f'{gene_name}', ha='left', fontsize=11,
               fontweight='bold')
        ax.text(x_min + 5, y_base - 0.6, f'{gene_info["old_locus"]}: {gene_info["pattern"]}',
               ha='left', fontsize=9, style='italic')

    # Add legend
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=colors['6mA'],
                   markersize=10, label='6mA site'),
        plt.Line2D([0], [0], marker='s', color='w', markerfacecolor=colors['4mC'],
                   markersize=10, label='4mC site'),
        mpatches.Patch(facecolor='#9B59B6', label='-35/-10 boxes'),
        plt.Line2D([0], [0], color='green', linewidth=2, label='TSS'),
        mpatches.Patch(facecolor='yellow', alpha=0.3, label='Promoter region')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)

    # Scale bar
    ax.text(0, -1.0, '0', ha='center', fontsize=9)
    ax.text(-300, -1.0, '-300', ha='center', fontsize=9)
    ax.text(100, -1.0, '+100', ha='center', fontsize=9)
    ax.text(-100, -1.3, 'Position relative to TSS (bp)', ha='center', fontsize=10)

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(-1.5, len(key_genes) * y_spacing)
    ax.axis('off')

    ax.set_title('Promoter Architecture and Methylation Site Positions\n'
                 'Key Epigenetically Regulated Genes in S. coelicolor M145',
                 fontsize=14, fontweight='bold', pad=20)

    plt.tight_layout()

    output_path = os.path.join(OUTPUT_DIR, 'promoter_architecture_summary.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight', facecolor='white')
    plt.close()

    return output_path

def main():
    print("=" * 60)
    print("Methylation Motif Position Schematic Generator")
    print("=" * 60)

    # Load data
    print("\nLoading methylation data...")
    methyl_df = load_methylation_data()
    print(f"Loaded {len(methyl_df)} methylation sites")

    # Generate individual schematics
    print("\nGenerating individual gene schematics...")
    for gene_name, gene_info in KEY_GENES.items():
        print(f"  Processing {gene_name}...")
        output_path, n_sites = create_gene_schematic(gene_name, gene_info, methyl_df)
        print(f"    Found {n_sites} {gene_info['mod_type']} sites in promoter")
        print(f"    Saved: {output_path}")

    # Generate combined schematic
    print("\nGenerating combined schematic...")
    output_path = create_combined_schematic(KEY_GENES, methyl_df)
    print(f"  Saved: {output_path}")

    # Generate promoter architecture summary
    print("\nGenerating promoter architecture summary...")
    output_path = create_promoter_architecture_summary(KEY_GENES, methyl_df)
    print(f"  Saved: {output_path}")

    print("\n" + "=" * 60)
    print("Schematic Generation Complete")
    print(f"Output directory: {OUTPUT_DIR}")
    print("=" * 60)

if __name__ == "__main__":
    main()
