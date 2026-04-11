#!/usr/bin/env python3
"""
AAGCCCG Methylation Cascade: Promoter schematics + MTase dynamics + Genome-wide context
Publication-quality figure for peer-reviewed journal

Panels:
  A-1: afsS promoter architecture (4mC site on AAGCCCG)
  A-2: redZ promoter architecture (6mA site on AAGCCCG)
  B:   SC_RS17645 MTase expression vs methylation loss
  C:   Genome-wide AAGCCCG methylation summary
  D:   BGC regulator AAGCCCG status matrix
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# === Paths ===
BASE = Path('/Users/okaban/bioinfo/rna-seq/11_epigenome_integration')
OUTDIR = BASE / 'analysis/19_bgc_regulator_overview'
OUTDIR.mkdir(parents=True, exist_ok=True)

# === Publication style ===
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 9,
    'axes.titlesize': 11,
    'axes.labelsize': 10,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 7.5,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'xtick.major.size': 3,
    'ytick.major.size': 3,
})

# === Color palette (colorblind-friendly) ===
C_TSS = '#2166AC'
C_GENE = '#2E7D32'
C_10BOX = '#E65100'
C_35BOX = '#1565C0'
C_MOTIF_BG = '#FFF3E0'
C_MOTIF_BORDER = '#E65100'
C_MTASE = '#6A1B9A'
C_4MC = '#C62828'
C_6MA = '#1565C0'
C_T1 = '#C62828'
C_T2 = '#E65100'
C_T3 = '#2E7D32'
C_LOST_BG = '#FFEBEE'


def draw_promoter_schematic(ax, gene_name, gene_product, methyl_type,
                            sigma10_seq, sigma10_dist, aagcccg_tss_dist,
                            methyl_freqs, lfc_values, padj_values):
    """Draw publication-quality promoter schematic."""

    ax.set_xlim(-300, 160)
    ax.set_ylim(-1.8, 3.8)
    ax.set_aspect('auto')

    # --- DNA backbone ---
    ax.plot([-280, 140], [0, 0], color='#333333', linewidth=2.5, zorder=2,
            solid_capstyle='round')

    # --- TSS bent arrow ---
    ax.plot([0, 0], [0, 0.7], color=C_TSS, linewidth=1.8, zorder=3)
    ax.annotate('', xy=(18, 0.7), xytext=(0, 0.7),
                arrowprops=dict(arrowstyle='-|>', color=C_TSS, lw=1.8,
                                mutation_scale=12))
    ax.text(0, 0.95, 'TSS', ha='center', va='bottom', fontsize=9,
            fontweight='bold', color=C_TSS)

    # --- -10 box ---
    pos_10 = -abs(sigma10_dist)
    bw = 28
    rect_10 = FancyBboxPatch((pos_10 - bw/2, -0.28), bw, 0.56,
                              boxstyle='round,pad=2', facecolor='#FFE0B2',
                              edgecolor=C_10BOX, linewidth=1.2, zorder=3)
    ax.add_patch(rect_10)
    ax.text(pos_10, 0, f'−10\n{sigma10_seq}', ha='center', va='center',
            fontsize=7, fontweight='bold', color=C_10BOX, linespacing=1.1)

    # --- -35 box ---
    pos_35 = pos_10 - 17
    rect_35 = FancyBboxPatch((pos_35 - bw/2, -0.28), bw, 0.56,
                              boxstyle='round,pad=2', facecolor='#E3F2FD',
                              edgecolor=C_35BOX, linewidth=1.2, zorder=3)
    ax.add_patch(rect_35)
    ax.text(pos_35, 0, '−35', ha='center', va='center',
            fontsize=8, fontweight='bold', color=C_35BOX)

    # --- AAGCCCG highlight region ---
    mx = -abs(aagcccg_tss_dist)
    rect_motif = FancyBboxPatch((mx - 22, -0.65), 44, 1.3,
                                 boxstyle='round,pad=2',
                                 facecolor=C_MOTIF_BG, alpha=0.7,
                                 edgecolor=C_MOTIF_BORDER, linewidth=1.5,
                                 linestyle='--', zorder=1)
    ax.add_patch(rect_motif)

    # --- Methylation lollipops ---
    colors_tp = [C_T1, C_T2, C_T3]
    labels_tp = ['T1', 'T2', 'T3']
    x_offsets = [-7, 0, 7]

    for freq, col, lab, xo in zip(methyl_freqs, colors_tp, labels_tp, x_offsets):
        x = mx + xo
        if freq > 0:
            ax.plot([x, x], [0.4, 2.0], color=col, linewidth=1.3, zorder=4)
            ax.scatter(x, 2.0, s=90, c=col, edgecolors='#333333',
                       linewidth=0.8, zorder=5)
            ax.text(x, 2.25, f'{freq:.0f}%', ha='center', va='bottom',
                    fontsize=7.5, color=col, fontweight='bold')
        else:
            ax.plot([x, x], [0.4, 2.0], color=col, linewidth=0.8,
                    linestyle=':', zorder=4, alpha=0.5)
            ax.scatter(x, 2.0, s=90, c='white', edgecolors=col,
                       linewidth=1.0, zorder=5)
            ax.text(x, 2.25, '0%', ha='center', va='bottom',
                    fontsize=7.5, color=col, alpha=0.7)
        ax.text(x, 0.25, lab, ha='center', va='bottom', fontsize=6.5,
                color=col, fontweight='bold')

    # --- Motif label ---
    motif_color = C_4MC if methyl_type == '4mC' else C_6MA
    ax.text(mx, -1.1, f'AAGCCCG ({methyl_type})',
            ha='center', va='top', fontsize=8.5, fontweight='bold',
            color=motif_color,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#FFF9C4',
                      edgecolor=motif_color, alpha=0.9, linewidth=1.0))

    # --- Distance annotation ---
    ax.annotate('', xy=(0, -0.5), xytext=(mx, -0.5),
                arrowprops=dict(arrowstyle='<->', color='#666666', lw=0.8))
    ax.text(mx / 2, -0.7, f'{abs(aagcccg_tss_dist)} bp',
            ha='center', va='top', fontsize=7.5, color='#666666')

    # --- Gene body arrow ---
    ax.annotate('', xy=(140, 0), xytext=(50, 0),
                arrowprops=dict(arrowstyle='-|>', color=C_GENE,
                                lw=6, mutation_scale=12))
    ax.text(95, 0.45, gene_name, ha='center', va='bottom', fontsize=11,
            fontweight='bold', fontstyle='italic', color=C_GENE)

    # --- Expression info box ---
    comparisons = ['T2 vs T1', 'T3 vs T1', 'T3 vs T2']
    lines = []
    for comp, lfc, padj in zip(comparisons, lfc_values, padj_values):
        sig = '***' if padj < 0.001 else ('**' if padj < 0.01 else
              ('*' if padj < 0.05 else 'n.s.'))
        lines.append(f'{comp}: {lfc:+.2f} {sig}')
    expr_text = 'log₂FC\n' + '\n'.join(lines)

    ax.text(155, 3.2, expr_text, ha='right', va='top', fontsize=7,
            fontfamily='monospace', linespacing=1.3,
            bbox=dict(boxstyle='round,pad=0.4', facecolor='#F5F5F5',
                      edgecolor='#BDBDBD', alpha=0.95, linewidth=0.8))

    # --- Product annotation ---
    ax.text(-278, 3.2, f'{gene_name} ({gene_product})',
            ha='left', va='top', fontsize=8.5, color='#555555',
            fontstyle='italic')

    # --- Axis cleanup ---
    ax.set_yticks([])
    for spine in ['top', 'right', 'left']:
        ax.spines[spine].set_visible(False)
    ax.set_xlabel('Distance from TSS (bp)', fontsize=9)
    ax.tick_params(axis='x', labelsize=8)


def plot_mtase_dynamics(ax):
    """Panel B: SC_RS17645 expression vs methylation loss."""

    tp_labels = ['T1', 'T2', 'T3']
    tp_x = [0, 1, 2]

    # SC_RS17645 normalized counts
    mtase_counts = {
        'T1': [239.69, 233.97, 261.65],
        'T2': [56.42, 46.54, 52.25],
        'T3': [182.20, 93.93, 159.28],
    }
    mtase_mean = [np.mean(v) for v in mtase_counts.values()]
    mtase_std = [np.std(v) for v in mtase_counts.values()]

    # Methylation frequency
    afsS_methyl = [83.19, 0, 0]
    redZ_methyl = [58.98, 0, 0]

    # Left axis: MTase expression
    ax.errorbar(tp_x, mtase_mean, yerr=mtase_std, color=C_MTASE,
                linewidth=2.0, marker='s', markersize=7, capsize=4,
                capthick=1.2, zorder=5, label='SC_RS17645 (counts)')
    ax.fill_between(tp_x, [m - s for m, s in zip(mtase_mean, mtase_std)],
                    [m + s for m, s in zip(mtase_mean, mtase_std)],
                    color=C_MTASE, alpha=0.08)
    ax.set_ylabel('SC_RS17645\nnormalized counts', color=C_MTASE,
                  fontsize=9, fontweight='bold')
    ax.tick_params(axis='y', labelcolor=C_MTASE, labelsize=8)
    ax.set_ylim(0, 320)

    # LFC annotations
    ax.text(0.5, 75, 'LFC = −2.19***', fontsize=7, color=C_MTASE,
            ha='center', fontstyle='italic')
    ax.text(1.5, 155, 'LFC = +1.45***', fontsize=7, color=C_MTASE,
            ha='center', fontstyle='italic')

    # Right axis: methylation
    ax2 = ax.twinx()
    ax2.plot(tp_x, afsS_methyl, color=C_4MC, linewidth=1.8, marker='o',
             markersize=6, linestyle='--', label='afsS 4mC (%)', zorder=4)
    ax2.plot(tp_x, redZ_methyl, color=C_6MA, linewidth=1.8, marker='^',
             markersize=6, linestyle='--', label='redZ 6mA (%)', zorder=4)

    # "Lost" markers
    for tp_i in [1, 2]:
        ax2.scatter(tp_i, 0, marker='x', s=60, c=C_4MC, linewidth=1.5, zorder=6)
        ax2.scatter(tp_i, 0, marker='x', s=60, c=C_6MA, linewidth=1.5, zorder=6)

    ax2.set_ylabel('Methylation (%)', fontsize=9, fontweight='bold', color='#333333')
    ax2.set_ylim(-5, 105)
    ax2.tick_params(axis='y', labelsize=8)

    ax.set_xticks(tp_x)
    ax.set_xticklabels(tp_labels, fontsize=9, fontweight='bold')
    ax.set_xlim(-0.3, 2.3)

    # Legend
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    leg = ax.legend(lines1 + lines2, labels1 + labels2, loc='upper right',
                    fontsize=7, frameon=True, framealpha=0.95,
                    edgecolor='#CCCCCC', fancybox=False)
    leg.get_frame().set_linewidth(0.6)

    # Cascade annotation
    ax.text(0.7, 30, 'MTase ↓ → methylation lost',
            fontsize=7.5, ha='center', fontweight='bold', color=C_4MC,
            bbox=dict(boxstyle='round,pad=0.3', facecolor=C_LOST_BG,
                      edgecolor=C_4MC, alpha=0.9, linewidth=0.8))

    for spine in ['top']:
        ax.spines[spine].set_visible(False)
        ax2.spines[spine].set_visible(False)


def plot_aagcccg_genome_wide(ax_left, ax_right):
    """Panels C and D: Genome-wide AAGCCCG context."""

    # --- Panel C: Methylation summary ---
    methylated = 669
    unmethylated = 180
    total = methylated + unmethylated

    cats = ['Meth.', 'Unmeth.']
    vals = [methylated, unmethylated]
    colors = [C_MOTIF_BORDER, '#BDBDBD']

    bars = ax_left.bar(cats, vals, color=colors, edgecolor='#333333',
                       linewidth=0.8, width=0.55)

    for bar, val in zip(bars, vals):
        pct = val / total * 100
        ax_left.text(bar.get_x() + bar.get_width() / 2,
                     bar.get_height() + 12,
                     f'{val}\n({pct:.1f}%)',
                     ha='center', fontsize=8, fontweight='bold',
                     linespacing=1.2)

    ax_left.set_ylabel('AAGCCCG sites\n(TSS ± 500 bp)', fontsize=9)
    ax_left.set_ylim(0, 820)
    ax_left.tick_params(axis='x', labelsize=7.5, rotation=0)
    for spine in ['top', 'right']:
        ax_left.spines[spine].set_visible(False)

    # Key genes annotation inside bar
    ax_left.text(0, methylated * 0.40,
                 'incl. afsS (4mC)\n    redZ (6mA)',
                 ha='center', va='center', fontsize=6, color='white',
                 fontweight='bold')

    # --- Panel D: BGC regulator matrix ---
    regulators = {
        'afsS':       {'motif': True,  'meth': True,  'lost': True},
        'redZ':       {'motif': True,  'meth': True,  'lost': True},
        'afsR':       {'motif': False, 'meth': True,  'lost': False},
        'bldN':       {'motif': False, 'meth': True,  'lost': True},
        'actII-ORF4': {'motif': False, 'meth': False, 'lost': False},
        'redD':       {'motif': False, 'meth': False, 'lost': False},
        'cdaR':       {'motif': False, 'meth': False, 'lost': False},
        'cpkO':       {'motif': False, 'meth': False, 'lost': False},
        'bldA':       {'motif': False, 'meth': False, 'lost': False},
        'bldD':       {'motif': False, 'meth': False, 'lost': False},
        'adpA':       {'motif': False, 'meth': False, 'lost': False},
        'papR2':      {'motif': False, 'meth': False, 'lost': False},
    }

    names = list(regulators.keys())
    cols_label = ['AAGCCCG\nprom.', 'Methyl.', 'Lost']
    data = np.array([[int(v['motif']), int(v['meth']), int(v['lost'])]
                     for v in regulators.values()])

    col_colors = ['#FF9800', C_MOTIF_BORDER, C_4MC]

    for j in range(3):
        for i in range(len(names)):
            if data[i, j] == 1:
                ax_right.scatter(j, i, s=120, c=col_colors[j],
                                 edgecolors='#333333', linewidth=0.6,
                                 marker='o', zorder=5)
                ax_right.text(j, i, '✓', ha='center', va='center',
                              fontsize=8, color='white', fontweight='bold',
                              zorder=6)
            else:
                ax_right.scatter(j, i, s=60, c='#F5F5F5', edgecolors='#CCCCCC',
                                 linewidth=0.5, marker='o', zorder=3)

    ax_right.set_xticks(range(3))
    ax_right.set_xticklabels(cols_label, fontsize=6.5, fontweight='bold',
                              linespacing=1.0)
    ax_right.set_yticks(range(len(names)))
    ax_right.set_yticklabels(names, fontsize=8, fontstyle='italic')
    ax_right.invert_yaxis()
    ax_right.set_xlim(-0.5, 2.5)
    for spine in ['top', 'right', 'left', 'bottom']:
        ax_right.spines[spine].set_visible(False)
    ax_right.tick_params(axis='both', length=0)

    # Highlight rows with all three positive
    for i, name in enumerate(names):
        if regulators[name]['lost'] and regulators[name]['motif']:
            ax_right.axhspan(i - 0.4, i + 0.4, color='#FFF9C4', alpha=0.5,
                             zorder=0)


# ============================================================
# MAIN FIGURE
# ============================================================

fig = plt.figure(figsize=(200 / 25.4, 240 / 25.4))  # ~200mm × 240mm

gs = gridspec.GridSpec(4, 2,
                       height_ratios=[1.0, 1.0, 0.08, 1.4],
                       width_ratios=[1.2, 1],
                       hspace=0.50, wspace=0.40,
                       left=0.10, right=0.95, top=0.94, bottom=0.04)

# --- Panel A-1: afsS ---
ax_a1 = fig.add_subplot(gs[0, :])
draw_promoter_schematic(
    ax_a1,
    gene_name='afsS', gene_product='σ-like factor',
    methyl_type='4mC',
    sigma10_seq='TAGACT', sigma10_dist=12,
    aagcccg_tss_dist=172,
    methyl_freqs=[83.19, 0, 0],
    lfc_values=[-0.03, -1.39, -1.36],
    padj_values=[0.91, 3.4e-12, 2.3e-11],
)
ax_a1.text(-0.02, 1.12, 'A', transform=ax_a1.transAxes,
           fontsize=14, fontweight='bold', va='top')
ax_a1.set_title('afsS promoter — AAGCCCG 4mC lost at T2',
                fontsize=10, fontweight='bold', pad=8, loc='left',
                x=0.04)

# --- Panel A-2: redZ ---
ax_a2 = fig.add_subplot(gs[1, :])
draw_promoter_schematic(
    ax_a2,
    gene_name='redZ', gene_product='response regulator (Red CSR)',
    methyl_type='6mA',
    sigma10_seq='TAACGT', sigma10_dist=12,
    aagcccg_tss_dist=108,
    methyl_freqs=[58.98, 0, 0],
    lfc_values=[-2.25, -1.10, 1.11],
    padj_values=[1.3e-17, 1.7e-6, 3.5e-5],
)
ax_a2.text(-0.02, 1.12, 'B', transform=ax_a2.transAxes,
           fontsize=14, fontweight='bold', va='top')
ax_a2.set_title('redZ promoter — AAGCCCG 6mA lost at T2',
                fontsize=10, fontweight='bold', pad=8, loc='left',
                x=0.04)

# --- Panel B: MTase dynamics ---
ax_b = fig.add_subplot(gs[3, 0])
plot_mtase_dynamics(ax_b)
ax_b.text(-0.12, 1.08, 'C', transform=ax_b.transAxes,
          fontsize=14, fontweight='bold', va='top')
ax_b.set_title('SC_RS17645 expression vs\ntarget methylation',
               fontsize=9, fontweight='bold', pad=10)

# --- Panel C+D: Genome-wide ---
gs_cd = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[3, 1],
                                         width_ratios=[1, 1.2], wspace=0.55)
ax_c = fig.add_subplot(gs_cd[0])
ax_d = fig.add_subplot(gs_cd[1])
plot_aagcccg_genome_wide(ax_c, ax_d)

ax_c.text(-0.15, 1.12, 'D', transform=ax_c.transAxes,
          fontsize=14, fontweight='bold', va='top')
ax_c.set_title('AAGCCCG sites\n(genome-wide)', fontsize=8.5,
               fontweight='bold', pad=12)

ax_d.text(-0.10, 1.12, 'E', transform=ax_d.transAxes,
          fontsize=14, fontweight='bold', va='top')
ax_d.set_title('BGC regulator\nAAGCCCG status', fontsize=8.5,
               fontweight='bold', pad=12)

# --- Save ---
out = OUTDIR / 'aagcccg_cascade_analysis'
fig.savefig(f'{out}.pdf')
fig.savefig(f'{out}.svg')
fig.savefig(f'{out}.png', dpi=300)
plt.close()
print(f'Saved: {out}.pdf / .svg / .png')


# ============================================================
# Summary statistics (unchanged)
# ============================================================

aagcccg_df = pd.read_csv(
    BASE / 'analysis/14_aagcccg_promoter_analysis/aagcccg_promoter_counts.csv'
)
aagcccg_with = aagcccg_df[aagcccg_df['has_motif'] == True]

spatial_df = pd.read_csv(
    BASE / 'analysis/18_tss_analyses/motif_AAGCCCG_spatial_detail.csv'
)

total_sites = len(spatial_df)
methylated_sites = spatial_df['is_methylated'].sum()
unmethylated_sites = total_sites - methylated_sites
proximal = spatial_df[spatial_df['rel_pos'].abs() <= 200]
proximal_methylated = proximal['is_methylated'].sum()
proximal_total = len(proximal)

summary_stats = {
    'Metric': [
        'Total promoters analyzed',
        'Promoters with AAGCCCG motif',
        'AAGCCCG motif prevalence',
        'Total AAGCCCG sites (TSS ±500bp)',
        'Methylated AAGCCCG sites',
        'Unmethylated AAGCCCG sites',
        'Methylation rate',
        'TSS-proximal AAGCCCG sites (±200bp)',
        'TSS-proximal methylated',
        'TSS-proximal methylation rate',
        'Coordinated (methylation + expression)',
        'SC_RS17645 (MTase) LFC T2vT1',
        'SC_RS17645 (MTase) LFC T3vT1',
        'afsS 4mC site TSS distance',
        'redZ 6mA site TSS distance',
    ],
    'Value': [
        f'{len(aagcccg_df):,}',
        f'{len(aagcccg_with):,}',
        f'{len(aagcccg_with)/len(aagcccg_df)*100:.1f}%',
        f'{total_sites}',
        f'{methylated_sites}',
        f'{unmethylated_sites}',
        f'{methylated_sites/total_sites*100:.1f}%',
        f'{proximal_total}',
        f'{proximal_methylated}',
        f'{proximal_methylated/proximal_total*100:.1f}%' if proximal_total > 0 else 'N/A',
        f"{len(aagcccg_with[aagcccg_with['is_coordinated']==True])}",
        '-2.19 (padj = 6.5e-16)',
        '-0.71 (padj = 3.3e-3)',
        '172 bp upstream',
        '108 bp upstream',
    ],
}
summary_table = pd.DataFrame(summary_stats)
summary_table.to_csv(OUTDIR / 'aagcccg_cascade_summary.tsv', sep='\t', index=False)
print(f'Saved summary to {OUTDIR}/aagcccg_cascade_summary.tsv')

print('\n=== KEY FINDING ===')
print(f'AAGCCCG motif prevalence: {len(aagcccg_with)}/{len(aagcccg_df)} '
      f'promoters ({len(aagcccg_with)/len(aagcccg_df)*100:.1f}%)')
print(f'AAGCCCG methylation rate: {methylated_sites}/{total_sites} '
      f'sites ({methylated_sites/total_sites*100:.1f}%)')
coord = len(aagcccg_with[aagcccg_with['is_coordinated'] == True])
print(f'Coordinated genes: {coord}/{len(aagcccg_with)} '
      f'({coord/len(aagcccg_with)*100:.1f}%)')
