#!/usr/bin/env python3
"""
Methylation Motif – TSS – σ Factor Binding Site Spatial Relationship Analysis

Maps AAGCCCG (6mA) and CCGG (4mC) motif occurrences relative to TSS,
and quantifies their spatial relationship with -35/-10 box elements.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import re
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow
import matplotlib.ticker as ticker
import warnings
warnings.filterwarnings('ignore')

# ── Paths ──
TSS_TABLE = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/18_tss_analyses/comprehensive_tss_table.csv")
SIGMA_MOTIFS = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/18_tss_analyses/D1_sigma_motifs.csv")
METHYL_SITES = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/01_integration/high_confidence_sites_weighted.csv")
GENOME_FASTA = "/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/GCF_000203835.1_ASM20383v1_genomic.fna"
OUTPUT_DIR = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/18_tss_analyses")

plt.rcParams.update({
    'font.size': 10, 'axes.labelsize': 12, 'axes.titlesize': 13,
    'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
    'font.family': 'sans-serif',
})


def load_genome(fasta_path):
    """Load genome sequence."""
    seq_parts = []
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith('>'):
                seq_parts.append(line.strip())
    return ''.join(seq_parts).upper()


def revcomp(seq):
    comp = str.maketrans('ACGT', 'TGCA')
    return seq.translate(comp)[::-1]


def find_motif_in_promoters(genome_seq, genes_df, motif_pattern, motif_name, window_up=500, window_down=200):
    """Find all occurrences of a motif in promoter regions relative to TSS."""
    records = []
    motif_re = re.compile(motif_pattern)
    motif_rc = revcomp(motif_pattern) if not any(c in motif_pattern for c in '[]().+*?') else None

    for _, gene in genes_df.iterrows():
        tss = gene['tss']
        strand = gene['strand']
        gid = gene['gene_id']

        # Extract region around TSS
        if strand == '+':
            region_start = max(0, tss - window_up - 1)  # 0-indexed
            region_end = min(len(genome_seq), tss + window_down)
            region_seq = genome_seq[region_start:region_end]

            # Search forward strand
            for m in re.finditer(motif_pattern, region_seq):
                abs_pos = region_start + m.start()
                rel_pos = abs_pos - (tss - 1)  # TSS-relative (1-based to 0-based adj)
                records.append({
                    'gene_id': gid, 'motif': motif_name,
                    'abs_position': abs_pos + 1,  # 1-based
                    'rel_pos': rel_pos, 'match_strand': '+',
                    'gene_strand': strand,
                })
            # Search reverse complement
            if motif_rc and motif_rc != motif_pattern:
                for m in re.finditer(motif_rc, region_seq):
                    abs_pos = region_start + m.start()
                    rel_pos = abs_pos - (tss - 1)
                    records.append({
                        'gene_id': gid, 'motif': motif_name,
                        'abs_position': abs_pos + 1,
                        'rel_pos': rel_pos, 'match_strand': '-',
                        'gene_strand': strand,
                    })
        else:
            region_start = max(0, tss - window_down - 1)
            region_end = min(len(genome_seq), tss + window_up)
            region_seq = genome_seq[region_start:region_end]

            # For minus strand genes, we reverse the relative position
            for m in re.finditer(motif_pattern, region_seq):
                abs_pos = region_start + m.start()
                rel_pos = (tss - 1) - abs_pos
                records.append({
                    'gene_id': gid, 'motif': motif_name,
                    'abs_position': abs_pos + 1,
                    'rel_pos': rel_pos, 'match_strand': '+',
                    'gene_strand': strand,
                })
            if motif_rc and motif_rc != motif_pattern:
                for m in re.finditer(motif_rc, region_seq):
                    abs_pos = region_start + m.start()
                    rel_pos = (tss - 1) - abs_pos
                    records.append({
                        'gene_id': gid, 'motif': motif_name,
                        'abs_position': abs_pos + 1,
                        'rel_pos': rel_pos, 'match_strand': '-',
                        'gene_strand': strand,
                    })

    return pd.DataFrame(records)


def find_nearest_sigma(motif_row, sigma_df):
    """Find nearest -35 and -10 box for a given motif occurrence."""
    gid = motif_row['gene_id']
    rel = motif_row['rel_pos']
    gene_sigma = sigma_df[sigma_df['gene_id'] == gid]

    result = {}
    for box_type in ['-35_box', '-10_box']:
        sub = gene_sigma[gene_sigma['motif_type'] == box_type]
        if len(sub) > 0:
            distances = sub['rel_pos'].values - rel
            closest_idx = np.argmin(np.abs(distances))
            result[f'nearest_{box_type}_dist'] = distances[closest_idx]
            result[f'nearest_{box_type}_pos'] = sub['rel_pos'].values[closest_idx]
        else:
            result[f'nearest_{box_type}_dist'] = np.nan
            result[f'nearest_{box_type}_pos'] = np.nan
    return result


def main():
    print("=" * 60)
    print("Motif-TSS-σ Factor Spatial Relationship Analysis")
    print("=" * 60)

    # Load data
    genes_df = pd.read_csv(TSS_TABLE)
    sigma_df = pd.read_csv(SIGMA_MOTIFS)
    methyl_df = pd.read_csv(METHYL_SITES)
    genome_seq = load_genome(GENOME_FASTA)
    print(f"  Genes: {len(genes_df)}, σ motifs: {len(sigma_df)}")

    # ── 1. Find AAGCCCG and CCGG in promoter regions ──
    print("\n[1] Scanning motifs in promoter regions (-500 to +200 bp from TSS)...")

    aagcccg_df = find_motif_in_promoters(genome_seq, genes_df, 'AAGCCCG', 'AAGCCCG')
    ccgg_df = find_motif_in_promoters(genome_seq, genes_df, 'CCGG', 'CCGG')

    print(f"  AAGCCCG occurrences: {len(aagcccg_df)} (in {aagcccg_df['gene_id'].nunique()} genes)")
    print(f"  CCGG occurrences:    {len(ccgg_df)} (in {ccgg_df['gene_id'].nunique()} genes)")

    all_motifs = pd.concat([aagcccg_df, ccgg_df], ignore_index=True)

    # ── 2. Compute distance to nearest σ motifs ──
    print("\n[2] Computing distances to nearest -35/-10 boxes...")

    sigma_records = []
    for _, row in all_motifs.iterrows():
        dist_info = find_nearest_sigma(row, sigma_df)
        sigma_records.append(dist_info)

    sigma_dist_df = pd.DataFrame(sigma_records)
    all_motifs = pd.concat([all_motifs.reset_index(drop=True), sigma_dist_df], axis=1)
    all_motifs.to_csv(OUTPUT_DIR / 'motif_tss_sigma_spatial.csv', index=False)

    # ── 3. Summary statistics ──
    print("\n[3] Summary statistics:")

    for motif_name in ['AAGCCCG', 'CCGG']:
        sub = all_motifs[all_motifs['motif'] == motif_name]
        print(f"\n  === {motif_name} ===")
        print(f"  Total occurrences: {len(sub)}")
        print(f"  TSS-relative position: median={sub['rel_pos'].median():.0f}, mean={sub['rel_pos'].mean():.1f}")

        for box in ['-35_box', '-10_box']:
            col = f'nearest_{box}_dist'
            valid = sub[col].dropna()
            if len(valid) > 0:
                print(f"  Distance to nearest {box}: median={valid.median():.0f}, mean={valid.mean():.1f} bp")
                print(f"    Within ±5 bp:  {(valid.abs() <= 5).sum()} ({(valid.abs() <= 5).mean()*100:.1f}%)")
                print(f"    Within ±10 bp: {(valid.abs() <= 10).sum()} ({(valid.abs() <= 10).mean()*100:.1f}%)")
                print(f"    Overlapping (within motif length): {(valid.abs() <= len(motif_name)).sum()}")

    # ── 4. Methylated vs unmethylated motif positions ──
    print("\n[4] Comparing methylated vs unmethylated motif positions...")

    # Map actual methylation sites to motif occurrences
    methyl_positions = set(methyl_df['position'].values)

    for motif_name, motif_len in [('AAGCCCG', 7), ('CCGG', 4)]:
        sub = all_motifs[all_motifs['motif'] == motif_name].copy()
        # Check if any methylation site falls within the motif
        methylated = []
        for _, row in sub.iterrows():
            abs_pos = row['abs_position']
            is_meth = any(p in methyl_positions for p in range(abs_pos, abs_pos + motif_len))
            methylated.append(is_meth)
        sub['is_methylated'] = methylated

        n_meth = sum(methylated)
        n_unmeth = len(methylated) - n_meth
        print(f"\n  {motif_name}: methylated={n_meth}, unmethylated={n_unmeth}")

        if n_meth > 10 and n_unmeth > 10:
            meth_pos = sub[sub['is_methylated']]['rel_pos']
            unmeth_pos = sub[~sub['is_methylated']]['rel_pos']
            stat, p = stats.mannwhitneyu(meth_pos, unmeth_pos, alternative='two-sided')
            print(f"    Methylated median pos: {meth_pos.median():.0f} bp from TSS")
            print(f"    Unmethylated median pos: {unmeth_pos.median():.0f} bp from TSS")
            print(f"    Mann-Whitney U p-value: {p:.2e}")

        # Store for plotting
        sub.to_csv(OUTPUT_DIR / f'motif_{motif_name}_spatial_detail.csv', index=False)

    # ══════════════════════════════════════════════
    # FIGURES
    # ══════════════════════════════════════════════

    # ── Figure 1: Motif TSS-relative position distribution ──
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)

    for ax, (motif_name, color) in zip(axes, [('AAGCCCG', '#d6604d'), ('CCGG', '#2166ac')]):
        sub = all_motifs[all_motifs['motif'] == motif_name]
        ax.hist(sub['rel_pos'], bins=np.arange(-500, 201, 10), color=color,
                edgecolor='none', alpha=0.7, label=f'{motif_name} (n={len(sub)})')

        # Mark σ regions
        ax.axvline(0, color='green', ls='-', lw=1.5, label='TSS')
        ax.axvspan(-45, -25, alpha=0.15, color='purple', label='-35 box region')
        ax.axvspan(-15, -5, alpha=0.15, color='orange', label='-10 box region')
        ax.set_ylabel('Count')
        ax.set_title(f'{motif_name} motif distribution relative to TSS')
        ax.legend(fontsize=8, loc='upper left')

    axes[1].set_xlabel('Distance from TSS (bp)')
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / 'motif_tss_position_distribution.png')
    fig.savefig(OUTPUT_DIR / 'motif_tss_position_distribution.pdf')
    plt.close(fig)

    # ── Figure 2: Distance to nearest -35/-10 box ──
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))

    for col_idx, (motif_name, color) in enumerate([('AAGCCCG', '#d6604d'), ('CCGG', '#2166ac')]):
        sub = all_motifs[all_motifs['motif'] == motif_name]
        for row_idx, box_type in enumerate(['-35_box', '-10_box']):
            ax = axes[row_idx, col_idx]
            dist_col = f'nearest_{box_type}_dist'
            valid = sub[dist_col].dropna()

            if len(valid) > 0:
                ax.hist(valid.clip(-100, 100), bins=np.arange(-100, 101, 5),
                        color=color, edgecolor='none', alpha=0.7)
                ax.axvline(0, color='red', ls='--', lw=1, label='Co-localized')
                ax.axvspan(-len(motif_name), len(motif_name), alpha=0.1, color='red', label=f'±{len(motif_name)} bp (overlap)')

                median_d = valid.median()
                ax.axvline(median_d, color='orange', ls='--', lw=1, label=f'Median={median_d:.0f} bp')

            ax.set_title(f'{motif_name} → nearest {box_type}')
            ax.set_ylabel('Count')
            ax.legend(fontsize=7)

        axes[1, col_idx].set_xlabel('Distance to nearest σ70 box (bp)')

    fig.suptitle('Methylation motif distance to σ70 binding elements', fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / 'motif_sigma_distance_histogram.png')
    fig.savefig(OUTPUT_DIR / 'motif_sigma_distance_histogram.pdf')
    plt.close(fig)

    # ── Figure 3: Comprehensive spatial map (schematic) ──
    fig, ax = plt.subplots(figsize=(14, 6))

    # Draw promoter architecture
    y_base = 3.0
    ax.add_patch(Rectangle((-500, y_base - 0.15), 700, 0.3, fc='#f0f0f0', ec='grey', lw=0.5))
    ax.annotate('', xy=(200, y_base), xytext=(-500, y_base),
                arrowprops=dict(arrowstyle='->', color='grey', lw=1))

    # -35 and -10 boxes
    ax.add_patch(Rectangle((-42, y_base - 0.3), 15, 0.6, fc='purple', alpha=0.3, ec='purple', lw=1))
    ax.text(-35, y_base + 0.45, '-35 box', ha='center', fontsize=9, color='purple', fontweight='bold')
    ax.add_patch(Rectangle((-15, y_base - 0.3), 10, 0.6, fc='orange', alpha=0.3, ec='orange', lw=1))
    ax.text(-10, y_base + 0.45, '-10 box', ha='center', fontsize=9, color='orange', fontweight='bold')
    ax.axvline(0, color='green', ls='-', lw=2)
    ax.text(2, y_base + 0.45, 'TSS', ha='left', fontsize=9, color='green', fontweight='bold')

    # Plot AAGCCCG density as kernel density
    aag_sub = all_motifs[all_motifs['motif'] == 'AAGCCCG']
    ccgg_sub = all_motifs[all_motifs['motif'] == 'CCGG']

    from scipy.stats import gaussian_kde

    x_range = np.linspace(-500, 200, 700)

    if len(aag_sub) > 10:
        kde_aag = gaussian_kde(aag_sub['rel_pos'].values, bw_method=0.05)
        density_aag = kde_aag(x_range)
        density_aag = density_aag / density_aag.max() * 1.2  # normalize
        ax.fill_between(x_range, y_base - 1.5, y_base - 1.5 + density_aag,
                        alpha=0.5, color='#d6604d', label='AAGCCCG density')
        ax.plot(x_range, y_base - 1.5 + density_aag, color='#d6604d', lw=1.5)

    if len(ccgg_sub) > 10:
        kde_ccgg = gaussian_kde(ccgg_sub['rel_pos'].values, bw_method=0.05)
        density_ccgg = kde_ccgg(x_range)
        density_ccgg = density_ccgg / density_ccgg.max() * 1.2
        ax.fill_between(x_range, y_base - 3.0, y_base - 3.0 + density_ccgg,
                        alpha=0.5, color='#2166ac', label='CCGG density')
        ax.plot(x_range, y_base - 3.0 + density_ccgg, color='#2166ac', lw=1.5)

    # Add labels
    ax.text(-490, y_base - 1.0, 'AAGCCCG\n(6mA target)', fontsize=10, color='#d6604d', fontweight='bold')
    ax.text(-490, y_base - 2.5, 'CCGG\n(4mC target)', fontsize=10, color='#2166ac', fontweight='bold')

    # Extend σ region lines
    for box_pos, box_color in [(-35, 'purple'), (-10, 'orange')]:
        ax.axvline(box_pos, color=box_color, ls=':', lw=0.8, alpha=0.5, ymin=0, ymax=1)

    ax.set_xlim(-510, 210)
    ax.set_ylim(-1, 4.5)
    ax.set_xlabel('Distance from TSS (bp)', fontsize=12)
    ax.set_yticks([])
    ax.set_title('Methylation motif spatial distribution relative to promoter architecture', fontsize=13)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / 'motif_promoter_architecture_map.png')
    fig.savefig(OUTPUT_DIR / 'motif_promoter_architecture_map.pdf')
    plt.close(fig)

    # ── Figure 4: Methylated vs unmethylated motif TSS position ──
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for ax, motif_name in zip(axes, ['AAGCCCG', 'CCGG']):
        detail_df = pd.read_csv(OUTPUT_DIR / f'motif_{motif_name}_spatial_detail.csv')
        if 'is_methylated' not in detail_df.columns:
            continue

        meth = detail_df[detail_df['is_methylated'] == True]
        unmeth = detail_df[detail_df['is_methylated'] == False]

        bins = np.arange(-500, 201, 15)
        ax.hist(unmeth['rel_pos'], bins=bins, alpha=0.5, color='grey', label=f'Unmethylated (n={len(unmeth)})', density=True)
        if len(meth) > 0:
            ax.hist(meth['rel_pos'], bins=bins, alpha=0.7, color='red', label=f'Methylated (n={len(meth)})', density=True)

        ax.axvline(0, color='green', ls='-', lw=1.5)
        ax.axvspan(-45, -25, alpha=0.1, color='purple')
        ax.axvspan(-15, -5, alpha=0.1, color='orange')
        ax.set_xlabel('Distance from TSS (bp)')
        ax.set_ylabel('Density')
        ax.set_title(f'{motif_name}: methylated vs unmethylated')
        ax.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / 'motif_methylated_vs_unmethylated_position.png')
    fig.savefig(OUTPUT_DIR / 'motif_methylated_vs_unmethylated_position.pdf')
    plt.close(fig)

    # ── Figure 5: Zoomed view of motif-sigma overlap zone ──
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)

    for ax, (motif_name, color) in zip(axes, [('AAGCCCG', '#d6604d'), ('CCGG', '#2166ac')]):
        sub = all_motifs[all_motifs['motif'] == motif_name]
        zoom = sub[(sub['rel_pos'] >= -80) & (sub['rel_pos'] <= 30)]

        ax.hist(zoom['rel_pos'], bins=np.arange(-80, 31, 2), color=color,
                edgecolor='white', alpha=0.8, label=f'{motif_name} (n={len(zoom)})')
        ax.axvline(0, color='green', ls='-', lw=2, label='TSS')
        ax.axvspan(-45, -25, alpha=0.2, color='purple', label='-35 box')
        ax.axvspan(-15, -5, alpha=0.2, color='orange', label='-10 box')
        ax.set_ylabel('Count')
        ax.legend(fontsize=8)

    axes[1].set_xlabel('Distance from TSS (bp)')
    fig.suptitle('Motif distribution in core promoter region (zoomed)', fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / 'motif_core_promoter_zoom.png')
    fig.savefig(OUTPUT_DIR / 'motif_core_promoter_zoom.pdf')
    plt.close(fig)

    print("\n  All figures saved to", OUTPUT_DIR)
    print("=" * 60)
    print("Analysis complete.")


if __name__ == '__main__':
    main()
