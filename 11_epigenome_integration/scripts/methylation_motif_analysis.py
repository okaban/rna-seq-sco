#!/usr/bin/env python3
"""
Methylation Motif Analysis
- Extract sequences around methylation sites
- Prepare FASTA for MEME motif discovery
- Analyze motif composition

Streptomyces coelicolor A3(2) M145 Project
"""

import pandas as pd
import numpy as np
from collections import Counter
import os
import subprocess

# Configuration
INTEGRATION_DIR = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis"
OUTPUT_DIR = os.path.join(INTEGRATION_DIR, "motif_analysis")
GENOME_PATH = "/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/GCF_000203835.1_ASM20383v1_genomic.fna"
METHYL_PATH = os.path.join(INTEGRATION_DIR, "high_confidence_sites_weighted.csv")

os.makedirs(OUTPUT_DIR, exist_ok=True)

FLANK_SIZE = 15  # Extract ±15bp around methylation site

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

    print(f"Loaded {len(sequences)} chromosome(s)")
    for chrom, seq in sequences.items():
        print(f"  {chrom}: {len(seq):,} bp")

    return sequences

def reverse_complement(seq):
    """Get reverse complement of DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))

def extract_flanking_sequences(methyl_df, genome, flank_size=FLANK_SIZE):
    """Extract sequences around methylation sites"""
    sequences = []

    for _, row in methyl_df.iterrows():
        chrom = row['chrom']
        pos = int(row['position'])
        strand = row['strand']
        mod_type = row['mod_type']
        timepoint = row['timepoint']
        freq = row['weighted_mod_freq']

        if chrom not in genome:
            continue

        # Extract flanking sequence
        start = pos - flank_size - 1  # 0-based
        end = pos + flank_size

        if start < 0 or end > len(genome[chrom]):
            continue

        seq = genome[chrom][start:end]

        # Reverse complement if on minus strand
        if strand == '-':
            seq = reverse_complement(seq)

        sequences.append({
            'chrom': chrom,
            'position': pos,
            'strand': strand,
            'mod_type': mod_type,
            'timepoint': timepoint,
            'frequency': freq,
            'sequence': seq,
            'center_base': seq[flank_size]  # The methylated base
        })

    return pd.DataFrame(sequences)

def write_fasta(sequences_df, output_path, mod_type=None):
    """Write sequences to FASTA file"""
    if mod_type:
        df = sequences_df[sequences_df['mod_type'] == mod_type]
    else:
        df = sequences_df

    with open(output_path, 'w') as f:
        for idx, row in df.iterrows():
            header = f">{row['chrom']}_{row['position']}_{row['strand']}_{row['mod_type']}_{row['timepoint']}"
            f.write(f"{header}\n{row['sequence']}\n")

    print(f"Wrote {len(df)} sequences to {output_path}")
    return len(df)

def analyze_base_composition(sequences_df, mod_type):
    """Analyze base composition around methylation sites"""
    df = sequences_df[sequences_df['mod_type'] == mod_type]

    if len(df) == 0:
        return None

    # Initialize position-wise counts
    seq_len = len(df.iloc[0]['sequence'])
    position_counts = {pos: Counter() for pos in range(seq_len)}

    for _, row in df.iterrows():
        seq = row['sequence']
        for pos, base in enumerate(seq):
            position_counts[pos][base] += 1

    # Convert to DataFrame
    positions = list(range(-FLANK_SIZE, FLANK_SIZE + 1))
    composition = pd.DataFrame(index=positions, columns=['A', 'T', 'G', 'C', 'N'])

    for i, pos in enumerate(positions):
        total = sum(position_counts[i].values())
        for base in ['A', 'T', 'G', 'C', 'N']:
            composition.loc[pos, base] = position_counts[i][base] / total * 100 if total > 0 else 0

    return composition

def find_consensus_motif(sequences_df, mod_type, window=5):
    """Find consensus motif around methylation site"""
    df = sequences_df[sequences_df['mod_type'] == mod_type]

    if len(df) == 0:
        return None

    # Get sequences centered on methylation site
    center = FLANK_SIZE
    start = center - window
    end = center + window + 1

    motifs = [row['sequence'][start:end] for _, row in df.iterrows()]

    # Count k-mers
    kmer_counts = Counter(motifs)

    return kmer_counts.most_common(20)

def plot_base_composition(composition_df, mod_type, output_path):
    """Plot base composition around methylation site"""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(14, 6))

    positions = composition_df.index

    colors = {'A': '#2ECC71', 'T': '#E74C3C', 'G': '#F39C12', 'C': '#3498DB'}

    for base in ['A', 'T', 'G', 'C']:
        ax.plot(positions, composition_df[base], label=base, color=colors[base], linewidth=2)

    ax.axvline(0, color='black', linestyle='--', linewidth=1, label='Methylation site')
    ax.axhline(25, color='gray', linestyle=':', alpha=0.5)  # Expected random frequency

    ax.set_xlabel('Position relative to methylation site (bp)', fontsize=12)
    ax.set_ylabel('Base frequency (%)', fontsize=12)
    ax.set_title(f'{mod_type} Methylation Site Flanking Sequence Composition\n'
                 f'(n = {len(composition_df)} sites, ±{FLANK_SIZE}bp)', fontsize=14)
    ax.legend(loc='upper right')
    ax.set_xlim(-FLANK_SIZE, FLANK_SIZE)
    ax.set_ylim(0, 60)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()

    print(f"Saved composition plot: {output_path}")

def create_motif_logo_data(sequences_df, mod_type, output_path):
    """Create data for sequence logo visualization"""
    df = sequences_df[sequences_df['mod_type'] == mod_type]

    # Write centered sequences (±5bp) for logo
    center = FLANK_SIZE
    with open(output_path, 'w') as f:
        for idx, row in df.iterrows():
            seq = row['sequence'][center-5:center+6]  # 11bp centered
            f.write(f">{row['chrom']}_{row['position']}\n{seq}\n")

    print(f"Wrote logo sequences: {output_path}")

def run_meme_analysis(fasta_path, output_dir, mod_type):
    """Run MEME for de novo motif discovery (if installed)"""
    meme_output = os.path.join(output_dir, f"meme_{mod_type}")

    # Check if MEME is available
    try:
        result = subprocess.run(['which', 'meme'], capture_output=True, text=True)
        if result.returncode != 0:
            print("MEME not found in PATH. Skipping MEME analysis.")
            print("To install: conda install -c bioconda meme")
            return None

        # Run MEME
        cmd = [
            'meme', fasta_path,
            '-dna',
            '-oc', meme_output,
            '-mod', 'zoops',
            '-nmotifs', '5',
            '-minw', '4',
            '-maxw', '12',
            '-revcomp'
        ]

        print(f"Running MEME for {mod_type}...")
        subprocess.run(cmd, capture_output=True, text=True)
        print(f"MEME output: {meme_output}")
        return meme_output

    except Exception as e:
        print(f"MEME analysis failed: {e}")
        return None

def analyze_known_motifs(sequences_df, mod_type):
    """Check for known bacterial methylation motifs"""
    # Common bacterial methylation motifs
    known_motifs = {
        '6mA': {
            'GATC': 'Dam methylase (E. coli)',
            'GAATTC': 'EcoRI-like',
            'GANTC': 'HinfI-like',
            'CTAG': 'Type II R-M',
            'CAGCTG': 'PvuII-like',
        },
        '4mC': {
            'CCGG': 'MspI-like',
            'GCGC': 'HhaI-like',
            'CCWGG': 'EcoRII-like (W=A/T)',
            'GATC': 'Dcm-like',
            'CCSGG': 'Type II R-M',
        }
    }

    df = sequences_df[sequences_df['mod_type'] == mod_type]
    results = []

    motifs_to_check = known_motifs.get(mod_type, {})

    for motif, description in motifs_to_check.items():
        # Check centered region
        center = FLANK_SIZE
        count = 0

        for _, row in df.iterrows():
            seq = row['sequence']
            # Check if motif is present near center
            for offset in range(-len(motif)//2, len(motif)//2 + 1):
                start = center + offset
                end = start + len(motif)
                if 0 <= start and end <= len(seq):
                    if seq[start:end] == motif or \
                       reverse_complement(seq[start:end]) == motif:
                        count += 1
                        break

        pct = count / len(df) * 100 if len(df) > 0 else 0
        results.append({
            'motif': motif,
            'description': description,
            'count': count,
            'total': len(df),
            'percentage': pct
        })

    return pd.DataFrame(results)

def main():
    print("=" * 60)
    print("Methylation Motif Analysis")
    print("=" * 60)

    # Load genome
    genome = load_genome(GENOME_PATH)

    # Load methylation sites
    print("\nLoading methylation sites...")
    methyl_df = pd.read_csv(METHYL_PATH)
    print(f"Loaded {len(methyl_df)} methylation sites")

    # Get unique sites (collapse across timepoints)
    unique_sites = methyl_df.drop_duplicates(subset=['chrom', 'position', 'strand', 'mod_type'])
    print(f"Unique methylation sites: {len(unique_sites)}")

    # Extract flanking sequences
    print(f"\nExtracting ±{FLANK_SIZE}bp flanking sequences...")
    sequences_df = extract_flanking_sequences(unique_sites, genome)
    print(f"Extracted {len(sequences_df)} sequences")

    # Save sequences
    sequences_df.to_csv(os.path.join(OUTPUT_DIR, 'methylation_site_sequences.csv'), index=False)

    # Analyze each modification type
    for mod_type in ['6mA', '4mC']:
        print(f"\n{'='*60}")
        print(f"Analyzing {mod_type} sites")
        print("=" * 60)

        mod_df = sequences_df[sequences_df['mod_type'] == mod_type]
        print(f"Total {mod_type} sites: {len(mod_df)}")

        if len(mod_df) == 0:
            continue

        # Write FASTA
        fasta_path = os.path.join(OUTPUT_DIR, f'{mod_type}_sites.fasta')
        write_fasta(sequences_df, fasta_path, mod_type)

        # Write logo sequences
        logo_path = os.path.join(OUTPUT_DIR, f'{mod_type}_logo_sequences.fasta')
        create_motif_logo_data(sequences_df, logo_path, mod_type)

        # Analyze base composition
        print("\nAnalyzing base composition...")
        composition = analyze_base_composition(sequences_df, mod_type)
        if composition is not None:
            composition.to_csv(os.path.join(OUTPUT_DIR, f'{mod_type}_base_composition.csv'))
            plot_base_composition(composition, mod_type,
                                 os.path.join(OUTPUT_DIR, f'{mod_type}_base_composition.png'))

        # Find consensus motifs
        print("\nFinding consensus motifs (±5bp)...")
        consensus = find_consensus_motif(sequences_df, mod_type, window=5)
        if consensus:
            print(f"Top 10 most common {mod_type} motifs:")
            for motif, count in consensus[:10]:
                pct = count / len(mod_df) * 100
                print(f"  {motif}: {count} ({pct:.1f}%)")

            # Save to file
            with open(os.path.join(OUTPUT_DIR, f'{mod_type}_consensus_motifs.txt'), 'w') as f:
                f.write(f"# {mod_type} consensus motifs (±5bp around methylation site)\n")
                f.write(f"# Total sites: {len(mod_df)}\n\n")
                f.write("Motif\tCount\tPercentage\n")
                for motif, count in consensus:
                    pct = count / len(mod_df) * 100
                    f.write(f"{motif}\t{count}\t{pct:.2f}%\n")

        # Check known motifs
        print("\nChecking for known R-M system motifs...")
        known_results = analyze_known_motifs(sequences_df, mod_type)
        known_results.to_csv(os.path.join(OUTPUT_DIR, f'{mod_type}_known_motifs.csv'), index=False)
        print(known_results.to_string(index=False))

        # Try MEME analysis
        meme_result = run_meme_analysis(fasta_path, OUTPUT_DIR, mod_type)

    # Create summary report
    create_summary_report(sequences_df, OUTPUT_DIR)

    print("\n" + "=" * 60)
    print("Motif Analysis Complete")
    print(f"Output directory: {OUTPUT_DIR}")
    print("=" * 60)

def create_summary_report(sequences_df, output_dir):
    """Create summary report of motif analysis"""
    report_path = os.path.join(output_dir, 'MOTIF_ANALYSIS_REPORT.md')

    with open(report_path, 'w') as f:
        f.write("# Methylation Motif Analysis Report\n\n")
        f.write(f"**Date:** 2026-02-02\n")
        f.write(f"**Project:** *Streptomyces coelicolor* A3(2) M145\n\n")

        f.write("## Summary\n\n")
        f.write(f"| Modification | Total Sites | Unique Sequences |\n")
        f.write(f"|--------------|-------------|------------------|\n")

        for mod_type in ['6mA', '4mC']:
            mod_df = sequences_df[sequences_df['mod_type'] == mod_type]
            f.write(f"| {mod_type} | {len(mod_df)} | {len(mod_df.drop_duplicates('sequence'))} |\n")

        f.write("\n## Output Files\n\n")
        f.write("| File | Description |\n")
        f.write("|------|-------------|\n")
        f.write("| `*_sites.fasta` | Flanking sequences for MEME |\n")
        f.write("| `*_base_composition.csv` | Position-wise base frequencies |\n")
        f.write("| `*_base_composition.png` | Base composition plot |\n")
        f.write("| `*_consensus_motifs.txt` | Top k-mer motifs |\n")
        f.write("| `*_known_motifs.csv` | Known R-M motif matches |\n")

    print(f"\nSaved report: {report_path}")

if __name__ == "__main__":
    main()
