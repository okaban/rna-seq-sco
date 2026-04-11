#!/usr/bin/env python3
"""
06: REBASE & RefSeq Conservation Figure for Novel Motif Candidates
===================================================================
Compute O/E ratios for novel candidate motifs (CCGKCA, GAACCGG, CGGCAACC)
across 833 RefSeq Streptomyces genomes, and generate a comprehensive figure
combining REBASE conservation rates with genuswide O/E distributions.

Panel A: REBASE conservation rate (82 species) for all motifs
Panel B: 833-genome O/E distribution (boxplot) with M145 highlighted
"""

import re
import gzip
from pathlib import Path
from multiprocessing import Pool, cpu_count

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings('ignore')

# ── Paths ──
BASE = Path('/Users/okaban/bioinfo/rna-seq/11_epigenome_integration')
GENOME_DIR = BASE / 'data' / 'ncbi_genomes_all_species'
OUTPUT_DIR = BASE / 'analysis' / '23_expanded_motif_search'
FIG_DIR = OUTPUT_DIR / 'figures'
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Existing genuswide data (7 original motifs)
GENUSWIDE_CSV = BASE / 'analysis' / '21_genuswide_motif_conservation' / 'motif_site_density_genuswide.csv'
REBASE_CONSERVATION = BASE / 'analysis' / '21_genuswide_motif_conservation' / 'rebase_motif_conservation_matrix.csv'
M145_OE = OUTPUT_DIR / 'm145_novel_motif_oe.csv'

# ── Style ──
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 9,
    'axes.titlesize': 11,
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 8,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 0.7,
    'axes.spines.top': False,
    'axes.spines.right': False,
})

COL_NOVEL = '#D32F2F'       # Red – novel (this study)
COL_ESTABLISHED = '#1565C0' # Blue – established/REBASE-known
COL_M145 = '#D32F2F'        # M145 marker color

# ── IUPAC & motif utilities ──
IUPAC = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'},
    'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'G', 'C'},
    'W': {'A', 'T'}, 'K': {'G', 'T'}, 'M': {'A', 'C'},
    'B': {'C', 'G', 'T'}, 'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'},
    'N': {'A', 'C', 'G', 'T'},
}

COMPLEMENT = {
    'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
    'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
    'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
    'D': 'H', 'H': 'D', 'N': 'N',
}


def reverse_complement(seq):
    return ''.join(COMPLEMENT.get(b, 'N') for b in reversed(seq.upper()))


def iupac_to_regex(motif):
    base_map = {
        'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
        'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
        'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
        'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]',
    }
    return ''.join(base_map.get(b, b) for b in motif.upper())


def expected_count_gc(motif_seq, genome_size, gc_pct):
    gc = gc_pct / 100.0
    at = 1.0 - gc
    base_prob = {"A": at / 2, "T": at / 2, "C": gc / 2, "G": gc / 2}
    iupac_bases = {
        "A": ["A"], "C": ["C"], "G": ["G"], "T": ["T"],
        "R": ["A", "G"], "Y": ["C", "T"], "S": ["C", "G"], "W": ["A", "T"],
        "K": ["G", "T"], "M": ["A", "C"], "B": ["C", "G", "T"],
        "D": ["A", "G", "T"], "H": ["A", "C", "T"], "V": ["A", "C", "G"],
        "N": ["A", "C", "G", "T"],
    }
    prob = 1.0
    for ch in motif_seq.upper():
        bases = iupac_bases.get(ch, [ch])
        pos_prob = sum(base_prob[b] for b in bases)
        prob *= pos_prob

    rc = reverse_complement(motif_seq)
    if rc == motif_seq:
        return prob * genome_size
    else:
        return prob * genome_size * 2


def load_genome(fasta_path):
    parts = []
    opener = gzip.open if str(fasta_path).endswith('.gz') else open
    with opener(fasta_path, 'rt') as f:
        for line in f:
            if not line.startswith('>'):
                parts.append(line.strip().upper())
    return ''.join(parts)


def find_genome_fastas():
    """Find all genome FASTA files."""
    extracted = GENOME_DIR / 'extracted' / 'ncbi_dataset' / 'data'
    fastas = {}
    if extracted.exists():
        for d in extracted.iterdir():
            if d.is_dir() and d.name.startswith('GCF_'):
                for f in d.glob('*.fna*'):
                    fastas[d.name] = f
                    break
    return fastas


# ── Novel motifs to compute ──
NOVEL_MOTIFS = {
    'CCGKCA': iupac_to_regex('CCGKCA'),
    'GAACCGG': 'GAACCGG',
    'CGGCAACC': 'CGGCAACC',
}


def process_genome(args):
    """Count novel motifs in one genome."""
    accession, fasta_path = args
    try:
        genome = load_genome(fasta_path)
        genome_len = len(genome)
        gc_count = genome.count('G') + genome.count('C')
        gc_pct = gc_count / genome_len * 100 if genome_len > 0 else 0

        result = {
            'accession': accession,
            'size_mb': round(genome_len / 1e6, 2),
            'gc_pct': round(gc_pct, 1),
        }

        for motif_name, regex in NOVEL_MOTIFS.items():
            compiled = re.compile(regex)
            fwd_count = len(compiled.findall(genome))

            # Reverse complement count
            rc_motif = reverse_complement(motif_name)
            is_palindrome = (rc_motif == motif_name)

            if is_palindrome:
                total = fwd_count
            else:
                rc_regex = iupac_to_regex(rc_motif)
                rc_compiled = re.compile(rc_regex)
                rev_count = len(rc_compiled.findall(genome))
                total = fwd_count + rev_count

            density = total / (genome_len / 1e6) if genome_len > 0 else 0
            expected = expected_count_gc(motif_name, genome_len, gc_pct)
            oe = total / expected if expected > 0 else float('inf')

            result[f'{motif_name}_count'] = total
            result[f'{motif_name}_density'] = round(density, 1)
            result[f'{motif_name}_expected'] = round(expected, 1)
            result[f'{motif_name}_oe'] = round(oe, 3)

        return result
    except Exception as e:
        print(f"  Error processing {accession}: {e}")
        return None


def main():
    # ============================================================
    # Main analysis
    # ============================================================
    print("=" * 70)
    print("REBASE & RefSeq CONSERVATION FIGURE")
    print("=" * 70)

    # 1. Load existing data
    print("\n1. Loading existing genuswide data...")
    genuswide = pd.read_csv(GENUSWIDE_CSV)
    rebase_cons = pd.read_csv(REBASE_CONSERVATION)
    m145_oe = pd.read_csv(M145_OE)

    print(f"   Existing genuswide: {len(genuswide)} species, {len(genuswide.columns)} columns")
    print(f"   REBASE conservation: {len(rebase_cons)} motifs")

    # 2. Compute novel motif O/E across 833 genomes
    novel_csv = OUTPUT_DIR / 'novel_motif_genuswide_oe.csv'

    print("\n2. Computing novel motif O/E across 833 genomes...")
    fastas = find_genome_fastas()
    print(f"   Found {len(fastas)} genome FASTA files")

    n_workers = min(cpu_count(), 8)
    print(f"   Using {n_workers} workers for parallel processing...")

    tasks = list(fastas.items())
    with Pool(n_workers) as pool:
        results = pool.map(process_genome, tasks)

    novel_df = pd.DataFrame([r for r in results if r is not None])
    novel_df.to_csv(novel_csv, index=False)
    print(f"   Processed {len(novel_df)} genomes → {novel_csv}")

    # Merge with existing genuswide data for species names
    species_map = dict(zip(genuswide['accession'], genuswide['species']))
    novel_df['species'] = novel_df['accession'].map(species_map)

    # 3. Compute summary statistics
    print("\n3. Novel motif O/E summary across genus:")
    print(f"{'Motif':<12} {'Median O/E':>10} {'IQR':>14} {'Min':>8} {'Max':>8} {'M145 O/E':>10}")
    print("-" * 70)

    m145_oe_dict = dict(zip(m145_oe['motif'], m145_oe['O_E_ratio']))

    for motif in NOVEL_MOTIFS:
        col = f'{motif}_oe'
        vals = novel_df[col].replace([np.inf, -np.inf], np.nan).dropna()
        q25, q50, q75 = vals.quantile([0.25, 0.5, 0.75])
        m145_val = m145_oe_dict.get(motif, 'N/A')
        print(f"{motif:<12} {q50:>10.3f} [{q25:.3f}–{q75:.3f}] {vals.min():>8.3f} {vals.max():>8.3f} {m145_val:>10}")


    # ============================================================
    # 4. Generate figure
    # ============================================================
    print("\n4. Generating figure...")

    # Define motifs for figure (ordered: established first, then novel)
    all_motifs_ordered = [
        # Established
        ('CCGG', '4mC', 'established'),
        ('GATC', '6mA', 'established'),
        ('AAGCCCG', '4mC+6mA', 'novel'),
        ('CCGKCA', '6mA', 'novel'),
        ('GAACCGG', '6mA', 'novel'),
        ('CGGCAACC', '6mA', 'novel'),
    ]

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={'height_ratios': [1, 1.5]})

    # ── Panel A: REBASE conservation rate ──
    ax_a = axes[0]

    # Existing REBASE data (from analysis 21)
    rebase_rates = {}
    rebase_ci = {}
    for _, row in rebase_cons.iterrows():
        motif = row['motif']
        rebase_rates[motif] = row['rate_pct']
        rebase_ci[motif] = (row['ci_low'], row['ci_high'])

    # Novel motifs: all 0% in REBASE
    for m, _, status in all_motifs_ordered:
        if m not in rebase_rates:
            rebase_rates[m] = 0.0
            rebase_ci[m] = (0.0, 0.0)

    labels_a = [f"{m}\n({mod})" for m, mod, _ in all_motifs_ordered]
    rates = [rebase_rates.get(m, 0) for m, _, _ in all_motifs_ordered]
    colors_a = [COL_ESTABLISHED if s == 'established' else COL_NOVEL for _, _, s in all_motifs_ordered]
    yerr_low = [rebase_rates.get(m, 0) - rebase_ci.get(m, (0, 0))[0] for m, _, _ in all_motifs_ordered]
    yerr_high = [rebase_ci.get(m, (0, 0))[1] - rebase_rates.get(m, 0) for m, _, _ in all_motifs_ordered]

    bars = ax_a.bar(range(len(labels_a)), rates, color=colors_a, edgecolor='black', linewidth=0.5,
                    yerr=[yerr_low, yerr_high], capsize=4, error_kw={'linewidth': 0.8})

    # Annotate
    for i, (m, _, s) in enumerate(all_motifs_ordered):
        rate = rates[i]
        if rate == 0:
            ax_a.text(i, 1.5, "0%\n(NOT IN\nREBASE)", ha='center', va='bottom',
                      fontsize=7, color=COL_NOVEL, fontweight='bold')
        else:
            n_match = int(round(rate * 82 / 100))
            ax_a.text(i, rate + 3, f"{rate:.1f}%\n({n_match}/82)",
                      ha='center', va='bottom', fontsize=7)

    ax_a.set_xticks(range(len(labels_a)))
    ax_a.set_xticklabels(labels_a, fontsize=8)
    ax_a.set_ylabel('REBASE conservation\nrate (%, 82 species)', fontsize=10)
    ax_a.set_title('A. REBASE R-M System Conservation in Streptomyces', fontsize=11, fontweight='bold', loc='left')
    ax_a.set_ylim(0, 35)

    # Legend
    legend_a = [
        mpatches.Patch(facecolor=COL_ESTABLISHED, edgecolor='black', linewidth=0.5, label='REBASE-known'),
        mpatches.Patch(facecolor=COL_NOVEL, edgecolor='black', linewidth=0.5, label='REBASE-novel (this study)'),
    ]
    ax_a.legend(handles=legend_a, loc='upper right', fontsize=8)

    # ── Panel B: O/E ratio distributions (833 genomes) ──
    ax_b = axes[1]

    # Collect O/E data for each motif
    oe_data = []
    oe_labels = []
    oe_colors = []

    for m, mod, status in all_motifs_ordered:
        col = f'{m}_oe'
        if col in genuswide.columns:
            vals = genuswide[col].replace([np.inf, -np.inf], np.nan).dropna().values
        elif col in novel_df.columns:
            vals = novel_df[col].replace([np.inf, -np.inf], np.nan).dropna().values
        else:
            vals = np.array([])
        oe_data.append(vals)
        oe_labels.append(f"{m}\n({mod})")
        oe_colors.append(COL_ESTABLISHED if status == 'established' else COL_NOVEL)

    # Create boxplots
    positions = range(len(oe_data))
    bps = ax_b.boxplot(oe_data, positions=positions, widths=0.6,
                        patch_artist=True, showfliers=False,
                        medianprops={'color': 'black', 'linewidth': 1.2},
                        whiskerprops={'linewidth': 0.7},
                        capprops={'linewidth': 0.7},
                        boxprops={'linewidth': 0.7})

    for patch, color in zip(bps['boxes'], oe_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.4)

    # M145 values (red star)
    for i, (m, _, _) in enumerate(all_motifs_ordered):
        m145_val = m145_oe_dict.get(m, None)
        if m145_val is not None:
            ax_b.plot(i, m145_val, '*', color=COL_M145, markersize=12,
                      markeredgecolor='black', markeredgewidth=0.5, zorder=5)

    # Reference lines
    ax_b.axhline(y=1.0, color='black', linestyle='--', linewidth=0.8, alpha=0.7)
    ax_b.axhline(y=0.75, color='gray', linestyle=':', linewidth=0.7, alpha=0.5)
    ax_b.axhline(y=1.5, color='gray', linestyle=':', linewidth=0.7, alpha=0.5)

    # Annotations for reference lines
    ax_b.text(len(oe_data) - 0.5, 1.02, 'O/E = 1.0', fontsize=7, ha='right', color='gray')
    ax_b.text(len(oe_data) - 0.5, 0.77, 'Avoidance', fontsize=7, ha='right', color='gray')
    ax_b.text(len(oe_data) - 0.5, 1.52, 'Maintenance', fontsize=7, ha='right', color='gray')

    ax_b.set_xticks(positions)
    ax_b.set_xticklabels(oe_labels, fontsize=8)
    ax_b.set_ylabel('Observed / Expected ratio\n(GC-content model)', fontsize=10)
    ax_b.set_title('B. O/E Ratio Distribution Across 833 Streptomyces Genomes', fontsize=11, fontweight='bold', loc='left')

    # Dynamic y-limit based on data
    all_vals = np.concatenate([d for d in oe_data if len(d) > 0])
    y_max = min(np.percentile(all_vals, 99), 5.0)
    ax_b.set_ylim(0, y_max)

    # Legend
    legend_b = [
        mpatches.Patch(facecolor=COL_ESTABLISHED, edgecolor='black', linewidth=0.5, alpha=0.4,
                       label='REBASE-known (833 genomes)'),
        mpatches.Patch(facecolor=COL_NOVEL, edgecolor='black', linewidth=0.5, alpha=0.4,
                       label='REBASE-novel (833 genomes)'),
        plt.Line2D([0], [0], marker='*', color='w', markerfacecolor=COL_M145,
                   markersize=10, markeredgecolor='black', label='M145'),
    ]
    ax_b.legend(handles=legend_b, loc='upper right', fontsize=8)

    plt.tight_layout()
    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(FIG_DIR / f'rebase_refseq_conservation.{ext}', dpi=300, bbox_inches='tight')
    print(f"   Saved: figures/rebase_refseq_conservation.{{pdf,svg,png}}")
    plt.close()

    # ============================================================
    # 5. Summary table
    # ============================================================
    print("\n" + "=" * 70)
    print("5. COMPREHENSIVE CONSERVATION SUMMARY")
    print("=" * 70)

    print(f"\n{'Motif':<12} {'Mod':<8} {'REBASE':>10} {'Genus O/E':>10} {'M145 O/E':>10} {'Novel?':>7}")
    print("-" * 65)

    for m, mod, status in all_motifs_ordered:
        rate = rebase_rates.get(m, 0)
        rebase_str = f"{rate:.1f}%" if rate > 0 else "0% (NOVEL)"

        col = f'{m}_oe'
        if col in genuswide.columns:
            genus_vals = genuswide[col].replace([np.inf, -np.inf], np.nan).dropna()
        elif col in novel_df.columns:
            genus_vals = novel_df[col].replace([np.inf, -np.inf], np.nan).dropna()
        else:
            genus_vals = pd.Series([])

        genus_med = f"{genus_vals.median():.3f}" if len(genus_vals) > 0 else "N/A"
        m145_val = m145_oe_dict.get(m, None)
        m145_str = f"{m145_val:.3f}" if m145_val else "N/A"
        novel_str = "YES" if status == 'novel' else "No"

        print(f"{m:<12} {mod:<8} {rebase_str:>10} {genus_med:>10} {m145_str:>10} {novel_str:>7}")

    print(f"\nAll results saved to: {OUTPUT_DIR}/")
    print("Done.")


if __name__ == '__main__':
    main()
