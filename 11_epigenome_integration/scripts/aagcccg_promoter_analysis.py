#!/usr/bin/env python3
"""
AAGCCCG Promoter Distribution Analysis
Quantitative analysis of AAGCCCG motif in promoter regions
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
from Bio import SeqIO
import re
import warnings
warnings.filterwarnings('ignore')

BASE_DIR = Path('/Users/okaban/bioinfo/rna-seq/11_epigenome_integration')
OUTPUT_DIR = BASE_DIR / 'analysis/14_aagcccg_promoter_analysis'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Motifs
AAGCCCG_MOTIF = 'AAGCCCG'  # 6mA motif
CCGG_MOTIF = 'CCGG'  # 4mC motif

def load_data():
    """Load methylation and expression data"""
    # Coordinated genes
    coord_genes = pd.read_csv(BASE_DIR / 'analysis/01_integration/T2vsT1_coordinated_genes.csv')
    print(f"Loaded {len(coord_genes)} T2vsT1 coordinated genes")

    # Full integrated data
    methyl_expr = pd.read_csv(BASE_DIR / 'analysis/01_integration/integrated_methyl_expression_weighted.csv')
    print(f"Loaded {len(methyl_expr)} genes with methylation data")

    # DESeq2 results
    deseq = pd.read_csv(Path('/Users/okaban/bioinfo/rna-seq/04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_2_vs_1.tsv'), sep='\t')
    print(f"Loaded {len(deseq)} genes with expression data")

    return coord_genes, methyl_expr, deseq

def load_genome_and_gff():
    """Load genome sequence and gene annotations"""
    genome_path = Path('/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/GCF_000203835.1_ASM20383v1_genomic.fna')
    gff_path = Path('/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/genomic.gff')

    # Load genome
    genome_seq = None
    for record in SeqIO.parse(genome_path, 'fasta'):
        genome_seq = str(record.seq)
        break
    print(f"Loaded genome: {len(genome_seq)} bp")

    # Load gene coordinates
    genes = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'gene':
                continue

            attrs = dict(item.split('=') for item in parts[8].split(';') if '=' in item)
            genes.append({
                'locus_tag': attrs.get('locus_tag', ''),
                'start': int(parts[3]),
                'end': int(parts[4]),
                'strand': parts[6],
            })

    genes_df = pd.DataFrame(genes)
    print(f"Loaded {len(genes_df)} genes from GFF")

    return genome_seq, genes_df

def extract_promoter_sequences(genome_seq, genes_df, upstream=300, downstream=50):
    """Extract promoter sequences for all genes"""
    promoters = []

    for _, gene in genes_df.iterrows():
        if gene['strand'] == '+':
            prom_start = max(0, gene['start'] - upstream - 1)
            prom_end = gene['start'] + downstream - 1
        else:
            prom_start = gene['end'] - downstream
            prom_end = min(len(genome_seq), gene['end'] + upstream)

        prom_seq = genome_seq[prom_start:prom_end]

        if gene['strand'] == '-':
            # Reverse complement
            complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
            prom_seq = ''.join(complement.get(b, 'N') for b in reversed(prom_seq))

        promoters.append({
            'locus_tag': gene['locus_tag'],
            'strand': gene['strand'],
            'prom_start': prom_start,
            'prom_end': prom_end,
            'sequence': prom_seq,
        })

    return pd.DataFrame(promoters)

def count_motifs_in_promoters(promoters_df, motif):
    """Count motif occurrences in promoter sequences"""
    results = []

    for _, prom in promoters_df.iterrows():
        seq = prom['sequence'].upper()

        # Find all occurrences
        positions = []
        for match in re.finditer(motif, seq):
            pos = match.start() - 300  # Relative to TSS (upstream = negative)
            positions.append(pos)

        # Also check reverse complement
        rev_motif = motif[::-1].translate(str.maketrans('ATGC', 'TACG'))
        for match in re.finditer(rev_motif, seq):
            pos = match.start() - 300
            positions.append(pos)

        results.append({
            'locus_tag': prom['locus_tag'],
            'motif_count': len(positions),
            'has_motif': len(positions) > 0,
            'positions': positions,
        })

    return pd.DataFrame(results)

def analyze_enrichment(motif_counts, coord_genes, all_genes):
    """Analyze motif enrichment in coordinated vs non-coordinated genes"""
    # Merge with coordination status
    coord_locus = set(coord_genes['gene_id'].values)
    motif_counts['is_coordinated'] = motif_counts['locus_tag'].isin(coord_locus)

    # Calculate enrichment
    coord_with_motif = motif_counts[motif_counts['is_coordinated'] & motif_counts['has_motif']]
    coord_without_motif = motif_counts[motif_counts['is_coordinated'] & ~motif_counts['has_motif']]
    noncoord_with_motif = motif_counts[~motif_counts['is_coordinated'] & motif_counts['has_motif']]
    noncoord_without_motif = motif_counts[~motif_counts['is_coordinated'] & ~motif_counts['has_motif']]

    # Contingency table
    contingency = [
        [len(coord_with_motif), len(coord_without_motif)],
        [len(noncoord_with_motif), len(noncoord_without_motif)]
    ]

    # Fisher's exact test
    odds_ratio, pvalue = stats.fisher_exact(contingency)

    return {
        'coord_with': len(coord_with_motif),
        'coord_without': len(coord_without_motif),
        'noncoord_with': len(noncoord_with_motif),
        'noncoord_without': len(noncoord_without_motif),
        'coord_freq': len(coord_with_motif) / (len(coord_with_motif) + len(coord_without_motif)) if (len(coord_with_motif) + len(coord_without_motif)) > 0 else 0,
        'noncoord_freq': len(noncoord_with_motif) / (len(noncoord_with_motif) + len(noncoord_without_motif)) if (len(noncoord_with_motif) + len(noncoord_without_motif)) > 0 else 0,
        'odds_ratio': odds_ratio,
        'pvalue': pvalue,
    }

def analyze_position_distribution(motif_counts):
    """Analyze position distribution of motifs relative to TSS"""
    all_positions = []

    for _, row in motif_counts.iterrows():
        for pos in row['positions']:
            all_positions.append(pos)

    positions_df = pd.DataFrame({'position': all_positions})

    # Binning
    bins = [-300, -200, -100, -50, -35, -10, 0, 50]
    labels = ['-300 to -200', '-200 to -100', '-100 to -50', '-50 to -35', '-35 to -10', '-10 to TSS', 'TSS to +50']
    positions_df['region'] = pd.cut(positions_df['position'], bins=bins, labels=labels)

    return positions_df

def create_visualizations(aagcccg_counts, ccgg_counts, aagcccg_positions, ccgg_positions, aagcccg_enrich, ccgg_enrich):
    """Create visualization figures"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # 1. Motif frequency comparison
    ax1 = axes[0, 0]
    motifs = ['AAGCCCG', 'CCGG']
    coord_freq = [aagcccg_enrich['coord_freq']*100, ccgg_enrich['coord_freq']*100]
    noncoord_freq = [aagcccg_enrich['noncoord_freq']*100, ccgg_enrich['noncoord_freq']*100]

    x = np.arange(len(motifs))
    width = 0.35
    ax1.bar(x - width/2, coord_freq, width, label='Coordinated genes', color='#e15759')
    ax1.bar(x + width/2, noncoord_freq, width, label='Non-coordinated genes', color='#4e79a7')
    ax1.set_ylabel('% promoters with motif')
    ax1.set_title('Motif Frequency in Promoters')
    ax1.set_xticks(x)
    ax1.set_xticklabels(motifs)
    ax1.legend()

    # Add p-values
    for i, (p_aag, p_ccgg) in enumerate([(aagcccg_enrich['pvalue'], ccgg_enrich['pvalue'])]):
        pass  # Will annotate below
    ax1.annotate(f"p={aagcccg_enrich['pvalue']:.3f}", (0, max(coord_freq[0], noncoord_freq[0])+1), ha='center')
    ax1.annotate(f"p={ccgg_enrich['pvalue']:.3f}", (1, max(coord_freq[1], noncoord_freq[1])+1), ha='center')

    # 2. Position distribution - AAGCCCG
    ax2 = axes[0, 1]
    if len(aagcccg_positions) > 0:
        aagcccg_positions['position'].hist(bins=30, ax=ax2, color='#f28e2b', edgecolor='black', alpha=0.7)
        ax2.axvline(x=-35, color='red', linestyle='--', label='-35 box')
        ax2.axvline(x=-10, color='blue', linestyle='--', label='-10 box')
        ax2.axvline(x=0, color='green', linestyle='-', label='TSS')
    ax2.set_xlabel('Position relative to TSS (bp)')
    ax2.set_ylabel('Count')
    ax2.set_title('AAGCCCG Position Distribution')
    ax2.legend()

    # 3. Position distribution - CCGG
    ax3 = axes[1, 0]
    if len(ccgg_positions) > 0:
        ccgg_positions['position'].hist(bins=30, ax=ax3, color='#4e79a7', edgecolor='black', alpha=0.7)
        ax3.axvline(x=-35, color='red', linestyle='--', label='-35 box')
        ax3.axvline(x=-10, color='blue', linestyle='--', label='-10 box')
        ax3.axvline(x=0, color='green', linestyle='-', label='TSS')
    ax3.set_xlabel('Position relative to TSS (bp)')
    ax3.set_ylabel('Count')
    ax3.set_title('CCGG Position Distribution')
    ax3.legend()

    # 4. Enrichment summary
    ax4 = axes[1, 1]
    enrichment_data = pd.DataFrame({
        'Motif': ['AAGCCCG', 'CCGG'],
        'Odds Ratio': [aagcccg_enrich['odds_ratio'], ccgg_enrich['odds_ratio']],
        'p-value': [aagcccg_enrich['pvalue'], ccgg_enrich['pvalue']],
    })

    colors = ['#f28e2b' if or_val > 1 else '#4e79a7' for or_val in enrichment_data['Odds Ratio']]
    ax4.barh(enrichment_data['Motif'], enrichment_data['Odds Ratio'], color=colors)
    ax4.axvline(x=1, color='black', linestyle='--')
    ax4.set_xlabel('Odds Ratio (coordinated vs non-coordinated)')
    ax4.set_title('Motif Enrichment in Coordinated Genes')

    for i, (or_val, pval) in enumerate(zip(enrichment_data['Odds Ratio'], enrichment_data['p-value'])):
        ax4.annotate(f'OR={or_val:.2f}, p={pval:.3f}', (or_val + 0.05, i), va='center')

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'aagcccg_promoter_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: aagcccg_promoter_analysis.png")

def main():
    print("="*60)
    print("AAGCCCG PROMOTER DISTRIBUTION ANALYSIS")
    print("="*60)

    # Load data
    print("\n1. Loading data...")
    coord_genes, methyl_expr, deseq = load_data()
    genome_seq, genes_df = load_genome_and_gff()

    # Extract promoter sequences
    print("\n2. Extracting promoter sequences...")
    promoters = extract_promoter_sequences(genome_seq, genes_df)
    print(f"   Extracted {len(promoters)} promoter sequences")

    # Count AAGCCCG motifs
    print("\n3. Counting AAGCCCG motifs...")
    aagcccg_counts = count_motifs_in_promoters(promoters, AAGCCCG_MOTIF)
    n_with_aagcccg = aagcccg_counts['has_motif'].sum()
    print(f"   {n_with_aagcccg} promoters ({n_with_aagcccg/len(aagcccg_counts)*100:.1f}%) contain AAGCCCG")

    # Count CCGG motifs
    print("\n4. Counting CCGG motifs...")
    ccgg_counts = count_motifs_in_promoters(promoters, CCGG_MOTIF)
    n_with_ccgg = ccgg_counts['has_motif'].sum()
    print(f"   {n_with_ccgg} promoters ({n_with_ccgg/len(ccgg_counts)*100:.1f}%) contain CCGG")

    # Enrichment analysis
    print("\n5. Analyzing enrichment in coordinated genes...")
    aagcccg_enrich = analyze_enrichment(aagcccg_counts, coord_genes, genes_df)
    ccgg_enrich = analyze_enrichment(ccgg_counts, coord_genes, genes_df)

    print(f"\n   AAGCCCG enrichment:")
    print(f"     Coordinated: {aagcccg_enrich['coord_freq']*100:.1f}%")
    print(f"     Non-coordinated: {aagcccg_enrich['noncoord_freq']*100:.1f}%")
    print(f"     Odds ratio: {aagcccg_enrich['odds_ratio']:.2f}")
    print(f"     p-value: {aagcccg_enrich['pvalue']:.4f}")

    print(f"\n   CCGG enrichment:")
    print(f"     Coordinated: {ccgg_enrich['coord_freq']*100:.1f}%")
    print(f"     Non-coordinated: {ccgg_enrich['noncoord_freq']*100:.1f}%")
    print(f"     Odds ratio: {ccgg_enrich['odds_ratio']:.2f}")
    print(f"     p-value: {ccgg_enrich['pvalue']:.4f}")

    # Position distribution
    print("\n6. Analyzing position distributions...")
    aagcccg_positions = analyze_position_distribution(aagcccg_counts)
    ccgg_positions = analyze_position_distribution(ccgg_counts)

    # Save results
    print("\n7. Saving results...")
    aagcccg_counts.to_csv(OUTPUT_DIR / 'aagcccg_promoter_counts.csv', index=False)
    ccgg_counts.to_csv(OUTPUT_DIR / 'ccgg_promoter_counts.csv', index=False)

    enrichment_summary = pd.DataFrame([
        {'motif': 'AAGCCCG', **aagcccg_enrich},
        {'motif': 'CCGG', **ccgg_enrich},
    ])
    enrichment_summary.to_csv(OUTPUT_DIR / 'motif_enrichment_summary.csv', index=False)

    # Create visualizations
    print("\n8. Creating visualizations...")
    create_visualizations(aagcccg_counts, ccgg_counts, aagcccg_positions, ccgg_positions, aagcccg_enrich, ccgg_enrich)

    # Summary
    print("\n" + "="*60)
    print("KEY FINDINGS")
    print("="*60)
    print(f"- AAGCCCG found in {aagcccg_enrich['coord_freq']*100:.1f}% of coordinated gene promoters")
    print(f"- CCGG found in {ccgg_enrich['coord_freq']*100:.1f}% of coordinated gene promoters")
    print(f"- AAGCCCG enrichment OR={aagcccg_enrich['odds_ratio']:.2f} (p={aagcccg_enrich['pvalue']:.4f})")
    print(f"- CCGG enrichment OR={ccgg_enrich['odds_ratio']:.2f} (p={ccgg_enrich['pvalue']:.4f})")

    print(f"\nOutput directory: {OUTPUT_DIR}")

if __name__ == '__main__':
    main()
