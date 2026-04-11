#!/usr/bin/env python3
"""
SC_RS17645 (SCO3104) Sequence & Structure Homology Analysis Figure
Publication-quality composite figure showing domain architecture,
BLAST conservation, structural homologs, and expression pattern.

Streptomyces coelicolor A3(2) M145 Project
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np
from pathlib import Path

# === Paths ===
BASE = Path('/Users/okaban/bioinfo/rna-seq/11_epigenome_integration')
OUTPUT_DIR = BASE / 'analysis' / '13_sc_rs17645_analysis'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# === Publication style ===
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 8,
    'axes.titlesize': 10,
    'axes.labelsize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 0.6,
})

# ====================================================================
# Panel A: Domain Architecture
# ====================================================================
def draw_domain_architecture(ax):
    """Draw protein domain architecture with annotations."""
    protein_len = 679
    bar_y = 0.5
    bar_h = 0.22

    # Full protein bar (gray background)
    ax.add_patch(FancyBboxPatch((0, bar_y - bar_h/2), protein_len, bar_h,
                                boxstyle="round,pad=2", fc='#E0E0E0', ec='black', lw=0.8))

    # Domain definitions: (start, end, color, label, label_pos)
    domains = [
        (1, 166, '#90CAF9', 'N-terminal\nregion', 'above'),
        (167, 393, '#C62828', 'N6_Mtase\n(PF02384)', 'inside'),
        (394, 538, '#FFE082', 'Linker', 'above'),
        (539, 669, '#0D47A1', 'TRD\n(specificity)', 'inside'),
    ]

    for start, end, color, label, pos in domains:
        width = end - start
        rect = FancyBboxPatch((start, bar_y - bar_h/2), width, bar_h,
                              boxstyle="round,pad=1", fc=color, ec='black', lw=0.6,
                              alpha=0.90)
        ax.add_patch(rect)

        mid_x = start + width / 2
        if pos == 'inside':
            ax.text(mid_x, bar_y, label, ha='center', va='center',
                    fontsize=7.5, fontweight='bold', color='white', linespacing=1.1)
        else:
            ax.text(mid_x, bar_y + bar_h/2 + 0.08, label, ha='center', va='bottom',
                    fontsize=7, fontweight='bold', color='#212121', linespacing=1.1)

    # Superfamily annotations (brackets below) — with white background for readability
    bracket_y = bar_y - bar_h/2 - 0.06
    # Type I RE MTase Subunit (IPR052916): 1-484
    ax.annotate('', xy=(1, bracket_y), xytext=(484, bracket_y),
                arrowprops=dict(arrowstyle='|-|', color='#7B1FA2', lw=1.2))
    ax.text(242, bracket_y - 0.05, 'Type I RE MTase Subunit (IPR052916)',
            ha='center', va='top', fontsize=7, fontweight='bold', color='#7B1FA2',
            fontstyle='italic',
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=1))

    # SAM-dependent MTase SF (IPR029063): 166-471
    sf_y = bracket_y - 0.15
    ax.annotate('', xy=(166, sf_y), xytext=(471, sf_y),
                arrowprops=dict(arrowstyle='|-|', color='#E65100', lw=1.0))
    ax.text(318, sf_y - 0.05, 'SAM-dep. MTase SF (IPR029063)',
            ha='center', va='top', fontsize=6.5, fontweight='bold', color='#E65100',
            fontstyle='italic',
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=1))

    # Motif I annotation (SAM binding) — place above the Linker, clear of domain bar
    ax.annotate('Motif I (FXGXG)\nSAM binding', xy=(394, bar_y + bar_h/2),
                xytext=(520, bar_y + 0.35),
                fontsize=6, ha='center', va='bottom', color='#C62828', fontweight='bold',
                arrowprops=dict(arrowstyle='->', color='#C62828', lw=0.8))

    # Scale bar
    for pos in [0, 100, 200, 300, 400, 500, 600, 679]:
        ax.plot([pos, pos], [bar_y - bar_h/2 - 0.28, bar_y - bar_h/2 - 0.30], 'k-', lw=0.5)
        if pos % 200 == 0 or pos == 679:
            ax.text(pos, bar_y - bar_h/2 - 0.33, str(pos), ha='center', va='top', fontsize=6)

    ax.set_xlim(-20, 720)
    ax.set_ylim(-0.15, 1.05)
    ax.set_xlabel('Amino acid position', fontsize=8)
    ax.set_title('A. Domain architecture of SC_RS17645 (SCO3104, 679 aa)',
                 fontsize=10, fontweight='bold', loc='left')
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

# ====================================================================
# Panel B: BLAST Top Hits (Sequence Conservation)
# ====================================================================
def draw_blast_results(ax):
    """Draw BLAST top hits as horizontal bar chart."""
    # Representative species from BLAST (best hit per species)
    species = [
        ('S. coelicolor\n(self)', 100.0, '#E53935'),
        ('S. lividans TK24', 99.1, '#EF5350'),
        ('S. anthocyanicus', 99.7, '#EF5350'),
        ('S. violaceoruber', 99.6, '#EF5350'),
        ('S. sp. RK76', 99.7, '#EF5350'),
        ('S. rubrogriseus', 98.5, '#FF7043'),
        ('S. sp. NPDC127039', 97.8, '#FF7043'),
        ('S. sp. NPDC059906', 97.6, '#FF7043'),
        ('S. sp. NPDC001599', 96.8, '#FFA726'),
        ('S. coelicoflavus', 93.4, '#FFA726'),
        ('S. tendae', 93.4, '#FFA726'),
        ('S. sp. ARC32', 92.6, '#FFB74D'),
    ]

    names = [s[0] for s in species]
    identities = [s[1] for s in species]
    colors = [s[2] for s in species]

    y_pos = np.arange(len(names))
    bars = ax.barh(y_pos, identities, height=0.65, color=colors, edgecolor='black', lw=0.4)

    # Value labels
    for i, (bar, val) in enumerate(zip(bars, identities)):
        ax.text(val - 0.5, i, f'{val:.1f}%', ha='right', va='center',
                fontsize=7, fontweight='bold', color='white')

    ax.set_yticks(y_pos)
    ax.set_yticklabels(names, fontsize=7)
    ax.set_xlabel('Sequence identity (%)', fontsize=8)
    ax.set_xlim(88, 102)
    ax.invert_yaxis()
    ax.set_title('B. NCBI BLAST top hits (nr, E = 0.0)',
                 fontsize=10, fontweight='bold', loc='left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Add annotation — positioned in lower-right to avoid overlap with bars
    ax.text(0.98, 0.02, 'All 50 hits:\nStreptomyces\n(>92% identity)',
            transform=ax.transAxes, fontsize=7, va='bottom', ha='right',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#FFF3E0', edgecolor='#F57C00', lw=0.6))

# ====================================================================
# Panel C: Foldseek Structural Homologs (PDB)
# ====================================================================
def draw_structural_homologs(ax):
    """Draw Foldseek PDB structural homolog comparison."""
    # Top PDB structural homologs from Foldseek
    pdb_hits = [
        ('7VS4_A', 'PacII M1M2S\n(DNA-m6A-SAH)', 24.1, 4.2e-22, '103-482', 'Type I R-M'),
        ('7VRU_A', 'PacII M1M2S\n(DNA-SAH)', 21.2, 1.4e-21, '67-482', 'Type I R-M'),
        ('7EEW_A', 'V. vulnificus\nMTase-Ocr-SAH', 16.1, 1.7e-21, '67-676', 'Type I R-M'),
        ('3KHK_B', 'M. mazei\nMM_0429', 20.6, 3.1e-17, '111-453', 'Type I R-M'),
        ('3LKD_A', 'S. thermophilus\nHsdM', 19.0, 9.6e-18, '48-453', 'Type I R-M'),
        ('2OKC_A', 'B. thetaiotaomicron\nStySJI M', 20.3, 2.9e-15, '67-447', 'Type I R-M'),
        ('7BTO_A', 'EcoR124I\nTranslocation', 18.0, 5.0e-14, '69-484', 'Type I R-M'),
        ('2ADM_B', 'M.TaqI\n(adenine MTase)', 15.9, 3.0e-11, '187-630', 'Type II'),
    ]

    # Create table-like visualization
    names = [h[0] for h in pdb_hits]
    descs = [h[1] for h in pdb_hits]
    seqids = [h[2] for h in pdb_hits]
    evalues = [h[3] for h in pdb_hits]
    ranges = [h[4] for h in pdb_hits]
    types = [h[5] for h in pdb_hits]

    y_pos = np.arange(len(names))
    colors = ['#1565C0' if t == 'Type I R-M' else '#F57C00' for t in types]

    bars = ax.barh(y_pos, seqids, height=0.65, color=colors, edgecolor='black', lw=0.4, alpha=0.8)

    for i, (bar, val, ev, rng) in enumerate(zip(bars, seqids, evalues, ranges)):
        ax.text(val + 0.3, i, f'{val:.1f}%  E={ev:.0e}  [{rng}]',
                ha='left', va='center', fontsize=6.5)

    ax.set_yticks(y_pos)
    labels = [f'{n}\n{d}' for n, d in zip(names, descs)]
    ax.set_yticklabels(labels, fontsize=6.5, linespacing=0.9)
    ax.set_xlabel('Sequence identity (%)', fontsize=8)
    ax.set_xlim(0, 38)
    ax.invert_yaxis()
    ax.set_title('C. Foldseek structural homologs (PDB, prob = 1.0)',
                 fontsize=10, fontweight='bold', loc='left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Legend
    type1_patch = mpatches.Patch(color='#1565C0', label='Type I R-M', alpha=0.8)
    type2_patch = mpatches.Patch(color='#F57C00', label='Type II', alpha=0.8)
    ax.legend(handles=[type1_patch, type2_patch], fontsize=7, loc='lower right',
              framealpha=0.9, edgecolor='gray')

# ====================================================================
# Panel D: Expression Pattern
# ====================================================================
def draw_expression(ax):
    """Draw expression pattern across timepoints with methylation correlation."""
    comparisons = ['T2 vs T1', 'T3 vs T1', 'T3 vs T2']
    log2fc = [-2.19, -0.71, 1.45]
    padj = [6.48e-16, 3.33e-3, 1.48e-7]
    sig = ['***', '**', '***']

    colors_bar = ['#1565C0', '#42A5F5', '#E53935']

    bars = ax.bar(comparisons, log2fc, width=0.55, color=colors_bar, edgecolor='black', lw=0.6)

    # Significance annotations
    for i, (bar, fc, p, s) in enumerate(zip(bars, log2fc, padj, sig)):
        y = fc + (0.15 if fc > 0 else -0.25)
        ax.text(i, y, f'{fc:.2f}\n({s})', ha='center', va='bottom' if fc > 0 else 'top',
                fontsize=8, fontweight='bold')

    ax.axhline(0, color='black', lw=0.5)
    ax.axhline(1, color='gray', lw=0.5, ls='--', alpha=0.5)
    ax.axhline(-1, color='gray', lw=0.5, ls='--', alpha=0.5)
    ax.set_ylabel('log$_2$FC', fontsize=9)
    ax.set_ylim(-3.0, 2.5)
    ax.set_title('D. SC_RS17645 expression changes',
                 fontsize=10, fontweight='bold', loc='left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Interpretation annotation — positioned in upper-right to avoid overlap with bars
    ax.text(0.98, 0.98, 'T2: MTase expression $\\downarrow$ $\\rightarrow$ AAGCCCG methylation $\\downarrow$\n'
            'T3: Partial recovery (log$_2$FC = +1.45 vs T2)',
            transform=ax.transAxes, fontsize=6.5, va='top', ha='right',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#E3F2FD', edgecolor='#1565C0', lw=0.5))


# ====================================================================
# Compose Figure
# ====================================================================
def main():
    fig = plt.figure(figsize=(180/25.4, 240/25.4))

    # Create grid: Panel A (top full width), B+C (middle row), D (bottom)
    gs = fig.add_gridspec(3, 2, hspace=0.45, wspace=0.40,
                          height_ratios=[1.0, 1.6, 0.9])

    ax_a = fig.add_subplot(gs[0, :])    # Domain architecture, full width
    ax_b = fig.add_subplot(gs[1, 0])    # BLAST hits
    ax_c = fig.add_subplot(gs[1, 1])    # Structural homologs
    ax_d = fig.add_subplot(gs[2, 0])    # Expression pattern

    draw_domain_architecture(ax_a)
    draw_blast_results(ax_b)
    draw_structural_homologs(ax_c)
    draw_expression(ax_d)

    # Add summary text box in bottom right
    ax_summary = fig.add_subplot(gs[2, 1])
    ax_summary.axis('off')

    summary_text = (
        "Summary: SC_RS17645 / SCO3104\n"
        "━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
        "UniProt: Q9F2P6 (679 aa, 72.2 kDa)\n"
        "REBASE: M.ScoA3ORF3104P\n"
        "AlphaFold: AF-Q9F2P6-F1 (pLDDT=87.1)\n"
        "\n"
        "Classification:\n"
        "  Type I R-M system M subunit\n"
        "  N-6 adenine methyltransferase\n"
        "\n"
        "Predicted motif: AAGCCCG (6mA)\n"
        "BLAST: Streptomyces-specific\n"
        "  (top 50 hits >92% identity)\n"
        "Closest PDB: 7VS4 (PacII, 24.1%)"
    )
    ax_summary.text(0.05, 0.95, summary_text, transform=ax_summary.transAxes,
                    fontsize=7.5, va='top', ha='left', family='monospace',
                    bbox=dict(boxstyle='round,pad=0.5', facecolor='#F5F5F5',
                              edgecolor='#666666', lw=0.6))

    # Save
    out_base = OUTPUT_DIR / 'fig_sc_rs17645_homology_composite'
    fig.savefig(f'{out_base}.pdf')
    fig.savefig(f'{out_base}.svg')
    fig.savefig(f'{out_base}.png', dpi=300)
    plt.close()
    print(f'Saved: {out_base}.pdf/.svg/.png')

    # Also save individual BLAST and Foldseek tables as CSV
    import csv

    # BLAST table
    blast_csv = OUTPUT_DIR / 'blast_top_hits.csv'
    with open(blast_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['rank', 'accession', 'organism', 'description', 'pct_identity', 'qcov_pct', 'evalue', 'bit_score'])
        blast_data = [
            (1, 'WP_011028768', 'Streptomyces (MULTISPECIES)', 'N-6 DNA methylase', 100.0, 100, 0.0, 1302),
            (2, 'WP_093456308', 'Streptomyces sp. 2114.2', 'N-6 DNA methylase', 99.9, 100, 0.0, 1302),
            (3, 'WP_359571097', 'Streptomyces anthocyanicus', 'N-6 DNA methylase', 99.7, 100, 0.0, 1300),
            (4, 'WP_210984784', 'Streptomyces sp. RK76', 'N-6 DNA methylase', 99.7, 100, 0.0, 1298),
            (5, 'WP_189284310', 'Streptomyces violaceoruber group', 'N-6 DNA methylase', 99.7, 100, 0.0, 1298),
            (6, 'WP_332056517', 'Streptomyces violaceoruber', 'N-6 DNA methylase', 99.4, 100, 0.0, 1295),
            (7, 'EFD68789', 'Streptomyces lividans TK24', 'Type II R-M system DNA adenine-specific methylase', 99.1, 100, 0.0, 1290),
            (8, 'WP_382820331', 'Streptomyces rubrogriseus', 'N-6 DNA methylase', 98.5, 100, 0.0, 1283),
            (9, 'WP_042799561', 'Streptomyces lividans', 'N-6 DNA methylase', 98.1, 100, 0.0, 1272),
            (10, 'EOY48171', 'Streptomyces lividans 1326', 'Type I R-M system, DNA-MTase subunit M', 97.9, 99, 0.0, 1254),
            (11, 'WP_386904923', 'Streptomyces sp. NPDC127039', 'N-6 DNA methylase', 97.8, 100, 0.0, 1273),
            (12, 'WP_385359592', 'Streptomyces sp. NPDC059906', 'N-6 DNA methylase', 97.6, 100, 0.0, 1272),
            (13, 'WP_356859319', 'Streptomyces sp. NPDC006617', 'N-6 DNA methylase', 97.2, 100, 0.0, 1267),
            (14, 'WP_398328706', 'Streptomyces sp. NPDC018026', 'N-6 DNA methylase', 96.8, 100, 0.0, 1265),
            (15, 'WP_388523046', 'Streptomyces sp. NPDC001599', 'N-6 DNA methylase', 96.8, 100, 0.0, 1261),
            (16, 'WP_381565166', 'Streptomyces coelicoflavus', 'N-6 DNA methylase', 93.4, 100, 0.0, 1198),
            (17, 'WP_385935194', 'Streptomyces tendae', 'N-6 DNA methylase', 93.4, 100, 0.0, 1193),
            (18, 'XKK58901', 'Streptomyces sp. ARC32', 'N-6 DNA methylase', 92.6, 100, 0.0, 1193),
        ]
        for row in blast_data:
            w.writerow(row)
    print(f'Saved: {blast_csv}')

    # Foldseek table
    foldseek_csv = OUTPUT_DIR / 'foldseek_structural_homologs.csv'
    with open(foldseek_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['rank', 'pdb_id', 'chain', 'description', 'organism', 'rm_type',
                     'seq_identity_pct', 'evalue', 'prob', 'query_range', 'target_len'])
        foldseek_data = [
            (1, '7VS4', 'A', 'PacII M1M2S-DNA(m6A)-SAH complex', 'Providencia alcalifaciens', 'Type I',
             24.1, 4.2e-22, 1.0, '103-482', 495),
            (2, '7VRU', 'A', 'PacII M1M2S-DNA-SAH complex', 'Providencia alcalifaciens', 'Type I',
             21.2, 1.4e-21, 1.0, '67-482', 497),
            (3, '7EEW', 'A', 'Intact MTase with Ocr and SAH', 'Vibrio vulnificus', 'Type I',
             16.1, 1.7e-21, 1.0, '67-676', 610),
            (4, '3KHK', 'B', 'Type I R-M system M subunit (MM_0429)', 'Methanosarchina mazei', 'Type I',
             20.6, 3.1e-17, 1.0, '111-453', 487),
            (5, '3LKD', 'A', 'Type I R-M system MTase subunit', 'Streptococcus thermophilus', 'Type I',
             19.0, 9.6e-18, 1.0, '48-453', 468),
            (6, '2OKC', 'A', 'Type I RE StySJI M protein', 'Bacteroides thetaiotaomicron', 'Type I',
             20.3, 2.9e-15, 1.0, '67-447', 424),
            (7, '7BTO', 'A', 'EcoR124I in Translocation State', 'Escherichia coli', 'Type I',
             18.0, 5.0e-14, 1.0, '69-484', 487),
            (8, '5YBB', 'B', 'Type I R-M complex assembly', 'Escherichia coli', 'Type I',
             21.1, 4.3e-14, 1.0, '82-482', 477),
            (9, '3S1S', 'A', 'Type IIG RE BpuSI', 'Bacillus pumilus', 'Type IIG',
             11.3, 1.4e-14, 1.0, '51-659', 791),
            (10, '2ADM', 'B', 'Adenine-N6-DNA-MTase M.TaqI', 'Thermus aquaticus', 'Type II',
             15.9, 3.0e-11, 1.0, '187-630', 496),
        ]
        for row in foldseek_data:
            w.writerow(row)
    print(f'Saved: {foldseek_csv}')


if __name__ == '__main__':
    main()
