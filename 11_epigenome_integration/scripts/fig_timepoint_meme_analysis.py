#!/usr/bin/env python3
"""
Timepoint-specific Motif Analysis and Sequence Logo Generation
for Streptomyces coelicolor M145 Methylation Data

This script:
1. Extracts sequences around methylation sites for each timepoint
2. Constructs Position Weight Matrices (PWMs) from known motif patterns
3. Generates publication-quality sequence logos
4. Analyzes temporal dynamics of motif prevalence

Output:
- Figure A: Sequence logos for top 3 motifs per modification type per timepoint
- Figure B: Temporal dynamics of motif prevalence

Author: Generated with Claude Code
Date: 2026-02-06
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

# Try to import logomaker, install if not available
try:
    import logomaker
    LOGOMAKER_AVAILABLE = True
except ImportError:
    LOGOMAKER_AVAILABLE = False
    print("Warning: logomaker not installed. Using custom logo implementation.")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration")
METHYL_PATH = BASE / "analysis/01_integration/high_confidence_sites_weighted.csv"
GENOME_PATH = Path("/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/GCF_000203835.1_ASM20383v1_genomic.fna")
OUT_DIR = BASE / "analysis/02_publication_figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Style - Publication Quality
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 10,
    "figure.dpi": 300,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.dpi": 300,
})

# Timepoints
TIMEPOINTS = ["T1", "T2", "T3"]
MOD_TYPES = ["4mC", "6mA"]

# Known motifs to analyze
KNOWN_MOTIFS = {
    "4mC": ["CCGG", "AAGCCCG", "TGGCCGGC"],  # Top motifs for 4mC (TGGCCGGC is extended CCGG)
    "6mA": ["CCGG", "AAGCCCG", "GCGC"],  # Top 3 for 6mA
}

# Motif colors
MOTIF_COLORS = {
    "CCGG":     "#E53935",
    "AAGCCCG":  "#1E88E5",
    "GCGC":     "#FB8C00",
    "GATC":     "#43A047",
    "TCGA":     "#8E24AA",
    "TGGCCGGC": "#7B1FA2",  # Purple for extended CCGG
}

# Timepoint colors
TP_COLORS = {
    "T1": "#4CAF50",  # Green
    "T2": "#2196F3",  # Blue
    "T3": "#FF5722",  # Deep Orange
}

# DNA base colors for sequence logos
BASE_COLORS = {
    'A': '#109648',  # Green
    'C': '#255C99',  # Blue
    'G': '#F7B32B',  # Yellow/Gold
    'T': '#D62839',  # Red
}

FLANK_SIZE = 10  # Extract +-10bp around methylation site


# ---------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------
def load_genome(genome_path):
    """Load genome sequence from FASTA"""
    print("Loading genome sequence...")
    sequences = {}
    current_id = None
    current_seq = []

    with open(genome_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
        if current_id:
            sequences[current_id] = ''.join(current_seq)

    print(f"  Loaded {len(sequences)} chromosome(s)")
    for chrom, seq in sequences.items():
        print(f"    {chrom}: {len(seq):,} bp")

    return sequences


def reverse_complement(seq):
    """Get reverse complement of DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


def extract_sequences_by_timepoint(methyl_df, genome, flank_size=FLANK_SIZE):
    """Extract flanking sequences for each timepoint and mod_type"""
    results = {}

    for tp in TIMEPOINTS:
        for mod in MOD_TYPES:
            key = (tp, mod)
            subset = methyl_df[(methyl_df['timepoint'] == tp) &
                               (methyl_df['mod_type'] == mod)]

            sequences = []
            for _, row in subset.iterrows():
                chrom = row['chrom']
                pos = int(row['position'])
                strand = row['strand']

                if chrom not in genome:
                    continue

                start = pos - flank_size - 1  # 0-based
                end = pos + flank_size

                if start < 0 or end > len(genome[chrom]):
                    continue

                seq = genome[chrom][start:end]

                if strand == '-':
                    seq = reverse_complement(seq)

                sequences.append(seq)

            results[key] = sequences
            print(f"  {tp} {mod}: {len(sequences)} sequences extracted")

    return results


def find_motif_containing_sequences(sequences, motif, center_window=5):
    """Find sequences containing the motif near center"""
    matching_seqs = []
    center = FLANK_SIZE

    for seq in sequences:
        # Check if motif is near center
        for offset in range(-center_window, center_window + 1):
            start = center + offset
            end = start + len(motif)
            if 0 <= start and end <= len(seq):
                if seq[start:end] == motif or reverse_complement(seq[start:end]) == motif:
                    # Align sequence to motif center
                    aligned_start = start
                    aligned_end = end
                    # Extract fixed window around motif
                    motif_center = start + len(motif) // 2
                    extract_start = max(0, motif_center - len(motif) - 2)
                    extract_end = min(len(seq), motif_center + len(motif) + 3)
                    matching_seqs.append(seq[extract_start:extract_end])
                    break

    return matching_seqs


def compute_pwm(sequences):
    """Compute Position Weight Matrix from sequences"""
    if not sequences:
        return None

    # Ensure all sequences have same length
    seq_len = min(len(s) for s in sequences)
    sequences = [s[:seq_len] for s in sequences]

    # Count bases at each position
    counts = {base: [0] * seq_len for base in 'ACGT'}

    for seq in sequences:
        for pos, base in enumerate(seq):
            if base in counts:
                counts[base][pos] += 1

    # Convert to DataFrame
    pwm = pd.DataFrame(counts)

    # Normalize to frequencies
    pwm = pwm.div(pwm.sum(axis=1), axis=0)

    # Add small pseudocount to avoid log(0)
    pwm = pwm + 0.0001
    pwm = pwm.div(pwm.sum(axis=1), axis=0)

    return pwm


def compute_information_content(pwm):
    """Compute information content matrix for sequence logo"""
    if pwm is None:
        return None

    # Background frequency (assume uniform)
    bg = 0.25

    # Compute information content
    ic = pwm.copy()
    for col in ic.columns:
        ic[col] = pwm[col] * np.log2(pwm[col] / bg)

    # Total IC at each position
    total_ic = ic.sum(axis=1)

    # Scale heights by IC
    heights = pwm.copy()
    for col in heights.columns:
        heights[col] = pwm[col] * total_ic

    return heights


def draw_logo_custom(ax, pwm, title="", fontsize=10):
    """Draw sequence logo without logomaker"""
    if pwm is None:
        ax.text(0.5, 0.5, "No data", ha='center', va='center', fontsize=fontsize)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 2)
        return

    heights = compute_information_content(pwm)

    positions = range(len(heights))

    for pos in positions:
        y_bottom = 0
        # Sort bases by height at this position
        pos_data = [(base, heights.loc[pos, base]) for base in 'ACGT']
        pos_data.sort(key=lambda x: x[1])

        for base, height in pos_data:
            if height > 0.01:
                ax.text(pos + 0.5, y_bottom + height/2, base,
                       fontsize=int(fontsize * min(1, height * 1.5)),
                       fontweight='bold',
                       ha='center', va='center',
                       color=BASE_COLORS[base])
                y_bottom += height

    ax.set_xlim(-0.5, len(heights) - 0.5)
    ax.set_ylim(0, 2)
    ax.set_xticks(range(len(heights)))
    ax.set_xticklabels(range(1, len(heights) + 1), fontsize=8)
    ax.set_ylabel("bits", fontsize=9)
    ax.set_title(title, fontsize=fontsize, fontweight='bold', loc='left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def draw_logo(ax, pwm, title="", n_sites=0, e_value=None):
    """Draw sequence logo using logomaker or custom implementation"""
    if pwm is None or len(pwm) == 0:
        ax.text(0.5, 0.5, "Insufficient data", ha='center', va='center',
                fontsize=10, color='gray')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 2)
        ax.axis('off')
        return

    heights = compute_information_content(pwm)

    if LOGOMAKER_AVAILABLE:
        try:
            logomaker.Logo(heights, ax=ax, color_scheme={
                'A': BASE_COLORS['A'],
                'C': BASE_COLORS['C'],
                'G': BASE_COLORS['G'],
                'T': BASE_COLORS['T'],
            })
        except Exception:
            draw_logo_custom(ax, pwm, title)
    else:
        draw_logo_custom(ax, pwm, title)

    ax.set_ylim(0, 2)
    ax.set_ylabel("bits", fontsize=8)
    ax.set_xlabel("")

    # Title with statistics
    title_text = title
    if n_sites > 0:
        title_text += f"\n(n={n_sites:,})"
    ax.set_title(title_text, fontsize=9, fontweight='bold', loc='left', pad=2)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def analyze_motif_prevalence(methyl_df, genome, flank_size=FLANK_SIZE):
    """Analyze motif prevalence across timepoints"""
    results = []

    for tp in TIMEPOINTS:
        for mod in MOD_TYPES:
            subset = methyl_df[(methyl_df['timepoint'] == tp) &
                               (methyl_df['mod_type'] == mod)]

            # Extract sequences
            sequences = []
            for _, row in subset.iterrows():
                chrom = row['chrom']
                pos = int(row['position'])
                strand = row['strand']

                if chrom not in genome:
                    continue

                start = pos - flank_size - 1
                end = pos + flank_size

                if start < 0 or end > len(genome[chrom]):
                    continue

                seq = genome[chrom][start:end]
                if strand == '-':
                    seq = reverse_complement(seq)
                sequences.append(seq)

            total = len(sequences)

            for motif in ["CCGG", "AAGCCCG", "GCGC", "GATC", "TCGA", "TGGCCGGC"]:
                matching = find_motif_containing_sequences(sequences, motif)
                count = len(matching)
                pct = (count / total * 100) if total > 0 else 0

                results.append({
                    'timepoint': tp,
                    'mod_type': mod,
                    'motif': motif,
                    'count': count,
                    'total': total,
                    'percentage': pct
                })

    return pd.DataFrame(results)


def create_figure_a(seq_data, genome, output_dir):
    """
    Create Figure A: Sequence Logo Comparison
    Layout: 3 rows (T1, T2, T3) x 6 columns (4mC x 3 motifs, 6mA x 3 motifs)
    Total 18 logos for publication
    """
    print("\nCreating Figure A: Sequence Logo Comparison (18 logos)...")

    # Create figure with 3 rows x 6 columns
    fig, axes = plt.subplots(3, 6, figsize=(18, 10))
    fig.subplots_adjust(hspace=0.45, wspace=0.3, top=0.88, bottom=0.06,
                        left=0.06, right=0.98)

    # Super-titles for modification types
    fig.text(0.27, 0.93, "4mC Methylation", ha='center', fontsize=14, fontweight='bold')
    fig.text(0.73, 0.93, "6mA Methylation", ha='center', fontsize=14, fontweight='bold')

    # Add separation line between 4mC and 6mA
    line = plt.Line2D([0.50, 0.50], [0.05, 0.90], transform=fig.transFigure,
                      color='gray', linewidth=1, linestyle='--', alpha=0.5)
    fig.add_artist(line)

    # Column labels for motifs
    motif_labels = ["CCGG", "AAGCCCG", "GCGC", "CCGG", "AAGCCCG", "GCGC"]
    for col, label in enumerate(motif_labels):
        x_pos = 0.09 + col * 0.155
        fig.text(x_pos, 0.90, label, ha='center', fontsize=10, fontweight='bold',
                color=MOTIF_COLORS[label])

    for row_idx, tp in enumerate(TIMEPOINTS):
        # Row label (timepoint)
        fig.text(0.02, 0.78 - row_idx * 0.28, tp, fontsize=14, fontweight='bold',
                va='center', ha='center')

        for col_idx in range(6):
            ax = axes[row_idx, col_idx]

            # Determine mod_type and motif
            if col_idx < 3:
                mod = "4mC"
                motif = KNOWN_MOTIFS["4mC"][col_idx]
            else:
                mod = "6mA"
                motif = KNOWN_MOTIFS["6mA"][col_idx - 3]

            # Get sequences
            key = (tp, mod)
            sequences = seq_data.get(key, [])

            if not sequences:
                ax.text(0.5, 0.5, "No data", ha='center', va='center',
                       fontsize=10, color='gray')
                ax.axis('off')
                continue

            # Find matching sequences
            matching_seqs = find_motif_containing_sequences(sequences, motif)

            if matching_seqs and len(matching_seqs) >= 5:
                pwm = compute_pwm(matching_seqs)
                n_sites = len(matching_seqs)

                # Draw logo
                if pwm is not None and len(pwm) > 0:
                    heights = compute_information_content(pwm)
                    if LOGOMAKER_AVAILABLE:
                        try:
                            logomaker.Logo(heights, ax=ax, color_scheme={
                                'A': BASE_COLORS['A'],
                                'C': BASE_COLORS['C'],
                                'G': BASE_COLORS['G'],
                                'T': BASE_COLORS['T'],
                            })
                        except Exception:
                            draw_logo_custom(ax, pwm)

                    ax.set_ylim(0, 2.2)
                    ax.set_ylabel("bits", fontsize=7)
                    ax.set_xlabel("")
                    ax.set_title(f"n={n_sites:,}", fontsize=8, loc='right', pad=2)
                    ax.spines['top'].set_visible(False)
                    ax.spines['right'].set_visible(False)
                    ax.tick_params(axis='x', labelsize=6)
                    ax.tick_params(axis='y', labelsize=6)
                else:
                    ax.text(0.5, 0.5, f"<5 sites", ha='center', va='center',
                           fontsize=9, color='gray')
                    ax.axis('off')
            else:
                ax.text(0.5, 0.5, f"<5 sites", ha='center', va='center',
                       fontsize=9, color='gray')
                ax.axis('off')

    # Save
    for ext in ["pdf", "svg", "png"]:
        out_path = output_dir / f"timepoint_meme_logos.{ext}"
        fig.savefig(out_path, dpi=300, bbox_inches='tight')
        print(f"  Saved: {out_path}")

    plt.close(fig)


def create_figure_b(prevalence_df, output_dir):
    """
    Create Figure B: Temporal Dynamics of Motif Prevalence
    Two panels: 4mC and 6mA, showing motif site counts across timepoints
    """
    print("\nCreating Figure B: Temporal Dynamics...")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.subplots_adjust(wspace=0.3, top=0.88, bottom=0.12, left=0.08, right=0.95)

    motif_order = ["CCGG", "AAGCCCG", "GCGC", "GATC", "TCGA"]

    for ax_idx, mod in enumerate(MOD_TYPES):
        ax = axes[ax_idx]
        mod_data = prevalence_df[prevalence_df['mod_type'] == mod]

        x_pos = np.arange(len(TIMEPOINTS))
        width = 0.15

        for motif_idx, motif in enumerate(motif_order):
            motif_data = mod_data[mod_data['motif'] == motif]
            counts = [motif_data[motif_data['timepoint'] == tp]['count'].values[0]
                     if len(motif_data[motif_data['timepoint'] == tp]) > 0 else 0
                     for tp in TIMEPOINTS]

            offset = (motif_idx - 2) * width
            bars = ax.bar(x_pos + offset, counts, width,
                         label=motif, color=MOTIF_COLORS[motif],
                         edgecolor='white', linewidth=0.5)

        # Add total site count annotations
        for i, tp in enumerate(TIMEPOINTS):
            tp_data = mod_data[mod_data['timepoint'] == tp]
            if len(tp_data) > 0:
                total = tp_data.iloc[0]['total']
                ax.text(i, ax.get_ylim()[1] * 0.95, f"n={total:,}",
                       ha='center', fontsize=9, fontweight='bold')

        ax.set_xticks(x_pos)
        ax.set_xticklabels(TIMEPOINTS, fontsize=11)
        ax.set_ylabel("Motif-containing sites", fontsize=11)
        ax.set_title(f"{mod} Methylation", fontsize=12, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.legend(fontsize=9, frameon=False, loc='upper right')

    fig.suptitle("Temporal Dynamics of Methylation Motifs",
                 fontsize=14, fontweight='bold', y=0.98)

    # Save
    for ext in ["pdf", "svg", "png"]:
        out_path = output_dir / f"timepoint_motif_dynamics.{ext}"
        fig.savefig(out_path, dpi=300, bbox_inches='tight')
        print(f"  Saved: {out_path}")

    plt.close(fig)


def create_figure_c(prevalence_df, output_dir):
    """
    Create Figure C: Fold Change of Major Motifs Relative to T1
    """
    print("\nCreating Figure C: Fold Change Analysis...")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.subplots_adjust(wspace=0.3, top=0.88, bottom=0.12, left=0.08, right=0.95)

    major_motifs = ["CCGG", "AAGCCCG", "GCGC"]
    markers = ['o', 's', 'D']

    for ax_idx, mod in enumerate(MOD_TYPES):
        ax = axes[ax_idx]
        mod_data = prevalence_df[prevalence_df['mod_type'] == mod]

        x_pos = np.arange(len(TIMEPOINTS))

        for motif_idx, motif in enumerate(major_motifs):
            motif_data = mod_data[mod_data['motif'] == motif]

            # Get T1 baseline
            t1_count = motif_data[motif_data['timepoint'] == 'T1']['count'].values
            if len(t1_count) == 0 or t1_count[0] == 0:
                continue
            t1_count = t1_count[0]

            fc_values = []
            for tp in TIMEPOINTS:
                tp_count = motif_data[motif_data['timepoint'] == tp]['count'].values
                if len(tp_count) > 0:
                    fc_values.append(tp_count[0] / t1_count)
                else:
                    fc_values.append(1.0)

            ax.plot(x_pos, fc_values, marker=markers[motif_idx],
                   markersize=8, linewidth=2,
                   color=MOTIF_COLORS[motif], label=motif)

        ax.axhline(1.0, color='gray', linestyle='--', linewidth=1, alpha=0.7)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(TIMEPOINTS, fontsize=11)
        ax.set_ylabel("Fold change (relative to T1)", fontsize=11)
        ax.set_title(f"{mod} Methylation", fontsize=12, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.legend(fontsize=9, frameon=False, loc='best')
        ax.set_ylim(0.4, 1.8)

    fig.suptitle("Temporal Fold Change of Major Methylation Motifs",
                 fontsize=14, fontweight='bold', y=0.98)

    # Save
    for ext in ["pdf", "svg", "png"]:
        out_path = output_dir / f"timepoint_motif_foldchange.{ext}"
        fig.savefig(out_path, dpi=300, bbox_inches='tight')
        print(f"  Saved: {out_path}")

    plt.close(fig)


def create_combined_figure(seq_data, prevalence_df, output_dir):
    """
    Create combined publication figure with all 18 logos and dynamics
    """
    print("\nCreating Combined Publication Figure...")

    # Create large figure for comprehensive view
    fig = plt.figure(figsize=(18, 18))

    # Define grid: top 3 rows for logos, bottom 2 rows for dynamics
    # Increased vertical spacing and margins
    gs = fig.add_gridspec(5, 6, hspace=0.55, wspace=0.35,
                          top=0.82, bottom=0.06, left=0.07, right=0.98)

    # ---- Panel A: Sequence Logos (top, 3x6 grid) ----
    fig.text(0.02, 0.95, "A", fontsize=18, fontweight='bold')
    fig.text(0.06, 0.95, "Sequence Logos by Timepoint and Motif", fontsize=13, fontweight='bold')

    # Modification type labels
    fig.text(0.28, 0.91, "4mC Methylation", ha='center', fontsize=12, fontweight='bold',
            color='#424242')
    fig.text(0.74, 0.91, "6mA Methylation", ha='center', fontsize=12, fontweight='bold',
            color='#424242')

    # Add vertical separator between 4mC and 6mA
    line = plt.Line2D([0.51, 0.51], [0.36, 0.90], transform=fig.transFigure,
                      color='gray', linewidth=1, linestyle='--', alpha=0.5)
    fig.add_artist(line)

    # Motif labels above each column (4mC has TGGCCGGC instead of GCGC)
    motif_labels_4mc = KNOWN_MOTIFS["4mC"]
    motif_labels_6ma = KNOWN_MOTIFS["6mA"]
    motif_labels = motif_labels_4mc + motif_labels_6ma
    for col, label in enumerate(motif_labels):
        x_pos = 0.11 + col * 0.152
        # Shorten TGGCCGGC display label to fit
        display_label = label if len(label) <= 7 else label[:4] + "..."
        fig.text(x_pos, 0.88, display_label, ha='center', fontsize=9, fontweight='bold',
                color=MOTIF_COLORS[label])

    # Calculate row positions for timepoint labels
    # Grid spans from bottom=0.06 to top=0.82, with 5 rows
    # Logo rows (0, 1, 2) occupy top 3/5 of grid
    grid_height = 0.82 - 0.06  # 0.76
    row_height = grid_height / 5  # 0.152 per row

    for row_idx, tp in enumerate(TIMEPOINTS):
        # Row label aligned with actual row center
        row_center = 0.82 - (row_idx + 0.5) * row_height
        fig.text(0.025, row_center, tp, fontsize=12, fontweight='bold',
                va='center', ha='center', color=TP_COLORS[tp])

        for col_idx in range(6):
            ax = fig.add_subplot(gs[row_idx, col_idx])

            if col_idx < 3:
                mod = "4mC"
                motif = KNOWN_MOTIFS["4mC"][col_idx]
            else:
                mod = "6mA"
                motif = KNOWN_MOTIFS["6mA"][col_idx - 3]

            key = (tp, mod)
            sequences = seq_data.get(key, [])

            if not sequences:
                ax.text(0.5, 0.5, "No data", ha='center', va='center',
                       fontsize=9, color='gray')
                ax.axis('off')
                continue

            matching_seqs = find_motif_containing_sequences(sequences, motif)

            if matching_seqs and len(matching_seqs) >= 5:
                pwm = compute_pwm(matching_seqs)
                n_sites = len(matching_seqs)

                if pwm is not None and len(pwm) > 0:
                    heights = compute_information_content(pwm)
                    if LOGOMAKER_AVAILABLE:
                        try:
                            logomaker.Logo(heights, ax=ax, color_scheme={
                                'A': BASE_COLORS['A'],
                                'C': BASE_COLORS['C'],
                                'G': BASE_COLORS['G'],
                                'T': BASE_COLORS['T'],
                            })
                        except Exception:
                            draw_logo_custom(ax, pwm)

                    ax.set_ylim(0, 2.2)
                    ax.set_ylabel("bits", fontsize=6)
                    ax.set_title(f"n={n_sites:,}", fontsize=7, loc='right', pad=1)
                    ax.spines['top'].set_visible(False)
                    ax.spines['right'].set_visible(False)
                    ax.tick_params(axis='both', labelsize=5)
            else:
                ax.text(0.5, 0.5, f"<5", ha='center', va='center',
                       fontsize=8, color='gray')
                ax.axis('off')

    # ---- Panel B: Temporal Dynamics (bottom left, spanning 2 columns each) ----
    # Panel B label positioned above the bar charts
    panel_b_top = 0.06 + 2 * row_height + 0.02  # Top of rows 3-4
    fig.text(0.02, panel_b_top, "B", fontsize=18, fontweight='bold')
    fig.text(0.06, panel_b_top, "Motif Site Counts", fontsize=13, fontweight='bold')

    for ax_idx, mod in enumerate(MOD_TYPES):
        ax = fig.add_subplot(gs[3:5, ax_idx * 2:(ax_idx + 1) * 2])
        mod_data = prevalence_df[prevalence_df['mod_type'] == mod]

        # Use modification-specific motif list
        motif_order = KNOWN_MOTIFS[mod]
        x_pos = np.arange(len(TIMEPOINTS))
        width = 0.25

        for motif_idx, motif in enumerate(motif_order):
            motif_data = mod_data[mod_data['motif'] == motif]
            counts = [motif_data[motif_data['timepoint'] == tp]['count'].values[0]
                     if len(motif_data[motif_data['timepoint'] == tp]) > 0 else 0
                     for tp in TIMEPOINTS]

            offset = (motif_idx - 1) * width
            # Shorten label for display
            display_label = motif if len(motif) <= 7 else motif[:6] + ".."
            ax.bar(x_pos + offset, counts, width,
                  label=display_label, color=MOTIF_COLORS[motif],
                  edgecolor='white', linewidth=0.5)

        # Total site annotations at bottom of bars
        for i, tp in enumerate(TIMEPOINTS):
            tp_data = mod_data[mod_data['timepoint'] == tp]
            if len(tp_data) > 0:
                total = tp_data.iloc[0]['total']
                ax.text(i, -0.05, f"n={total:,}",
                       ha='center', va='top', fontsize=8, fontweight='bold',
                       transform=ax.get_xaxis_transform())

        ax.set_xticks(x_pos)
        ax.set_xticklabels(TIMEPOINTS, fontsize=10)
        ax.set_ylabel("Sites", fontsize=10)
        ax.set_title(f"{mod}", fontsize=11, fontweight='bold', pad=8)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.legend(fontsize=8, frameon=False, loc='upper right')

    # ---- Panel C: Fold Change (bottom right) ----
    fig.text(0.68, panel_b_top, "C", fontsize=18, fontweight='bold')
    fig.text(0.72, panel_b_top, "Fold Change (vs T1)", fontsize=13, fontweight='bold')

    ax_fc = fig.add_subplot(gs[3:5, 4:6])

    # Distinct color palettes for 4mC (warm/red tones) and 6mA (cool/blue tones)
    mod_colors = {
        '4mC': {
            'CCGG': '#D32F2F',      # Dark red
            'AAGCCCG': '#F57C00',   # Orange
            'TGGCCGGC': '#C2185B',  # Pink
        },
        '6mA': {
            'CCGG': '#1976D2',      # Blue
            'AAGCCCG': '#00796B',   # Teal
            'GCGC': '#512DA8',      # Purple
        }
    }

    # Different markers: 4mC = filled, 6mA = hollow
    markers_4mc = {'CCGG': 'o', 'AAGCCCG': 's', 'TGGCCGGC': '^'}
    markers_6ma = {'CCGG': 'o', 'AAGCCCG': 's', 'GCGC': 'D'}

    for mod in MOD_TYPES:
        mod_data = prevalence_df[prevalence_df['mod_type'] == mod]

        # Use modification-specific motif list
        for motif in KNOWN_MOTIFS[mod]:
            motif_data = mod_data[mod_data['motif'] == motif]

            t1_count = motif_data[motif_data['timepoint'] == 'T1']['count'].values
            if len(t1_count) == 0 or t1_count[0] == 0:
                continue
            t1_count = t1_count[0]

            fc_values = []
            for tp in TIMEPOINTS:
                tp_count = motif_data[motif_data['timepoint'] == tp]['count'].values
                fc_values.append(tp_count[0] / t1_count if len(tp_count) > 0 else 1.0)

            # Get color and marker based on modification type
            color = mod_colors[mod].get(motif, '#666666')
            if mod == '4mC':
                marker = markers_4mc.get(motif, 'o')
                linestyle = '-'
                linewidth = 2.5
                markersize = 9
                markerfacecolor = color
                markeredgecolor = 'white'
                markeredgewidth = 1
            else:  # 6mA
                marker = markers_6ma.get(motif, 'o')
                linestyle = '--'
                linewidth = 2
                markersize = 8
                markerfacecolor = 'white'
                markeredgecolor = color
                markeredgewidth = 2

            # Shorten label for display
            display_motif = motif if len(motif) <= 7 else motif[:6] + ".."
            ax_fc.plot(range(len(TIMEPOINTS)), fc_values,
                      marker=marker, markersize=markersize, linewidth=linewidth,
                      linestyle=linestyle, color=color,
                      markerfacecolor=markerfacecolor,
                      markeredgecolor=markeredgecolor,
                      markeredgewidth=markeredgewidth,
                      label=f"{mod}-{display_motif}")

    ax_fc.axhline(1.0, color='gray', linestyle=':', linewidth=1, alpha=0.7)
    ax_fc.set_xticks(range(len(TIMEPOINTS)))
    ax_fc.set_xticklabels(TIMEPOINTS, fontsize=10)
    ax_fc.set_ylabel("Fold change", fontsize=10)
    ax_fc.set_title("", pad=8)  # Empty title for consistent padding
    ax_fc.spines['top'].set_visible(False)
    ax_fc.spines['right'].set_visible(False)

    # Two-column legend with modification type grouping
    ax_fc.legend(fontsize=7, frameon=True, fancybox=True, framealpha=0.9,
                ncol=2, loc='upper left', bbox_to_anchor=(0.0, 1.0),
                columnspacing=1.0, handlelength=2.5)
    ax_fc.set_ylim(0.35, 1.65)

    # Add shaded regions to emphasize 4mC decrease vs 6mA increase
    ax_fc.axhspan(0.35, 1.0, alpha=0.05, color='red', zorder=0)
    ax_fc.axhspan(1.0, 1.65, alpha=0.05, color='blue', zorder=0)

    # Add annotations for key observations (adjusted positions)
    ax_fc.annotate("4mC decreases\nat T3", xy=(2, 0.55), xytext=(2.3, 0.70),
                  fontsize=8, ha='left', color='#B71C1C', fontweight='bold',
                  arrowprops=dict(arrowstyle='->', color='#B71C1C', lw=1))
    ax_fc.annotate("6mA stable/\nincreases", xy=(2, 1.20), xytext=(2.3, 1.40),
                  fontsize=8, ha='left', color='#0D47A1', fontweight='bold',
                  arrowprops=dict(arrowstyle='->', color='#0D47A1', lw=1))

    # Save
    for ext in ["pdf", "svg", "png"]:
        out_path = output_dir / f"timepoint_meme_combined.{ext}"
        fig.savefig(out_path, dpi=300, bbox_inches='tight')
        print(f"  Saved: {out_path}")

    plt.close(fig)


def print_summary_statistics(prevalence_df):
    """Print summary statistics of motif analysis"""
    print("\n" + "=" * 70)
    print("TIMEPOINT-SPECIFIC MOTIF ANALYSIS SUMMARY")
    print("=" * 70)

    for mod in MOD_TYPES:
        print(f"\n{'─' * 40}")
        print(f"  {mod} Methylation")
        print(f"{'─' * 40}")

        mod_data = prevalence_df[prevalence_df['mod_type'] == mod]

        # Pivot table
        pivot = mod_data.pivot_table(index='motif', columns='timepoint',
                                     values='count', aggfunc='first')
        pivot = pivot[TIMEPOINTS]
        print("\nSite counts by motif and timepoint:")
        print(pivot.to_string())

        # Total sites
        print("\nTotal sites:")
        for tp in TIMEPOINTS:
            tp_data = mod_data[mod_data['timepoint'] == tp]
            if len(tp_data) > 0:
                print(f"  {tp}: {tp_data.iloc[0]['total']:,}")

        # Fold changes
        print("\nFold change (T3/T1):")
        for motif in ["CCGG", "AAGCCCG", "GCGC"]:
            motif_data = mod_data[mod_data['motif'] == motif]
            t1 = motif_data[motif_data['timepoint'] == 'T1']['count'].values
            t3 = motif_data[motif_data['timepoint'] == 'T3']['count'].values
            if len(t1) > 0 and len(t3) > 0 and t1[0] > 0:
                fc = t3[0] / t1[0]
                print(f"  {motif}: {fc:.2f}")


def main():
    """Main function"""
    print("=" * 70)
    print("Timepoint-Specific MEME Motif Analysis")
    print("Streptomyces coelicolor M145")
    print("=" * 70)

    # Load genome
    genome = load_genome(GENOME_PATH)

    # Load methylation data
    print("\nLoading methylation sites...")
    methyl_df = pd.read_csv(METHYL_PATH)
    print(f"  Total sites: {len(methyl_df):,}")

    # Summary by timepoint and mod_type
    print("\nSites by timepoint and modification:")
    summary = methyl_df.groupby(['timepoint', 'mod_type']).size().unstack(fill_value=0)
    print(summary)

    # Extract sequences by timepoint
    print("\nExtracting flanking sequences...")
    seq_data = extract_sequences_by_timepoint(methyl_df, genome)

    # Analyze motif prevalence
    print("\nAnalyzing motif prevalence...")
    prevalence_df = analyze_motif_prevalence(methyl_df, genome)

    # Save prevalence data
    prevalence_path = OUT_DIR / "timepoint_motif_prevalence.csv"
    prevalence_df.to_csv(prevalence_path, index=False)
    print(f"  Saved: {prevalence_path}")

    # Print summary statistics
    print_summary_statistics(prevalence_df)

    # Create figures
    create_figure_a(seq_data, genome, OUT_DIR)
    create_figure_b(prevalence_df, OUT_DIR)
    create_figure_c(prevalence_df, OUT_DIR)
    create_combined_figure(seq_data, prevalence_df, OUT_DIR)

    print("\n" + "=" * 70)
    print("Analysis Complete")
    print(f"Output directory: {OUT_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()
