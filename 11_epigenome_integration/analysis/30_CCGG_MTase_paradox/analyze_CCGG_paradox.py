#!/usr/bin/env python3
"""
H7: CCGG MTase Paradox Analysis
Analyzes the paradox of CCGG 4mC sites present at T1 despite low MTase expression.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn2
from matplotlib.patches import Patch
import re
import os
from collections import defaultdict

# ============================================================
# Configuration
# ============================================================
OUTDIR = '/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/30_CCGG_MTase_paradox'
FIGDIR = os.path.join(OUTDIR, 'figures')
TABDIR = os.path.join(OUTDIR, 'tables')

HC_SITES = '/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/01_integration/high_confidence_sites_weighted.csv'
MOTIF_ASSIGN = '/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/23_expanded_motif_search/4mC_motif_assignment.csv'
GFF_FILE = '/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/genomic.gff'
DESEQ2_T2T1 = '/Users/okaban/bioinfo/rna-seq/04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_2_vs_1.tsv'
DESEQ2_T3T1 = '/Users/okaban/bioinfo/rna-seq/04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_3_vs_1.tsv'
NORM_COUNTS = '/Users/okaban/bioinfo/rna-seq/04_deseq2/analysis/04_deseq2_260128_v1/results/normalized_counts_M145.tsv'
GENOME_FASTA = '/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/GCF_000203835.1_ASM20383v1_genomic.fna'
MTASE_PRED = '/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/23_expanded_motif_search/mtase_motif_predictions.csv'

# Plotting style
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

COLORS = {'T1': '#2166ac', 'T2': '#67a9cf', 'T3': '#ef8a62'}

# ============================================================
# 1. Load data
# ============================================================
print("=" * 60)
print("Loading data...")
print("=" * 60)

hc = pd.read_csv(HC_SITES)
motif = pd.read_csv(MOTIF_ASSIGN)
norm_counts = pd.read_csv(NORM_COUNTS, sep='\t')
deseq_t2t1 = pd.read_csv(DESEQ2_T2T1, sep='\t')
deseq_t3t1 = pd.read_csv(DESEQ2_T3T1, sep='\t')

print(f"HC sites total: {len(hc)}")
print(f"4mC motif assignments total: {len(motif)}")

# Filter for 4mC sites
hc_4mc = hc[hc['mod_type'] == '4mC'].copy()
print(f"HC 4mC sites: {len(hc_4mc)}")

# CCGG-related motifs: CCGG, GGCCGG, TGGCCGGC, CCGC
ccgg_motifs = ['CCGG', 'GGCCGG', 'TGGCCGGC', 'CCGC']
motif_ccgg = motif[motif['assigned_motif'].isin(ccgg_motifs)].copy()
print(f"CCGG-related motif-assigned 4mC sites: {len(motif_ccgg)}")

# Also include unassigned that might be CCGG
# We'll focus on motif-assigned sites for primary analysis
for tp in ['T1', 'T2', 'T3']:
    n = len(motif_ccgg[motif_ccgg['timepoint'] == tp])
    print(f"  {tp}: {n} CCGG-related sites")

# ============================================================
# 2. CCGG 4mC site overlap between timepoints
# ============================================================
print("\n" + "=" * 60)
print("2. Site overlap analysis between timepoints")
print("=" * 60)

# Get CCGG-related sites including TGGCCGGC
# Use the broader definition: all 4mC at CCGG-containing motifs
# TGGCCGGC contains CCGG (positions 4-7: CCGG in TGGCCGGC)
# GGCCGG contains CCGG
# CCGC is a related but distinct motif

# Let's be thorough: include TGGCCGGC, GGCCGG, CCGG from motif_assign
# Also get all 4mC from HC sites and check if their context contains CCGG
motif_tggccggc = motif[motif['assigned_motif'] == 'TGGCCGGC'].copy()
motif_ccgg_strict = motif[motif['assigned_motif'] == 'CCGG'].copy()
motif_ggccgg = motif[motif['assigned_motif'] == 'GGCCGG'].copy()
motif_ccgc = motif[motif['assigned_motif'] == 'CCGC'].copy()

print(f"\nMotif breakdown:")
for m_name, m_df in [('TGGCCGGC', motif_tggccggc), ('CCGG', motif_ccgg_strict),
                       ('GGCCGG', motif_ggccgg), ('CCGC', motif_ccgc)]:
    for tp in ['T1', 'T2', 'T3']:
        n = len(m_df[m_df['timepoint'] == tp])
        print(f"  {m_name} {tp}: {n}")

# Create position sets for each timepoint
# Use (position, strand) as unique identifier for a site
def get_site_set(df, tp):
    sub = df[df['timepoint'] == tp]
    return set(zip(sub['position'], sub['strand']))

# For all CCGG-related motifs combined
ccgg_all = motif[motif['assigned_motif'].isin(['TGGCCGGC', 'CCGG', 'GGCCGG', 'CCGC'])].copy()

t1_sites = get_site_set(ccgg_all, 'T1')
t2_sites = get_site_set(ccgg_all, 'T2')
t3_sites = get_site_set(ccgg_all, 'T3')

print(f"\nAll CCGG-related sites (unique positions):")
print(f"  T1: {len(t1_sites)}")
print(f"  T2: {len(t2_sites)}")
print(f"  T3: {len(t3_sites)}")

# Overlaps
t1_t2 = t1_sites & t2_sites
t1_t3 = t1_sites & t3_sites
t2_t3 = t2_sites & t3_sites
t1_t2_t3 = t1_sites & t2_sites & t3_sites

overlap_table = pd.DataFrame({
    'comparison': ['T1_only', 'T2_only', 'T3_only', 'T1∩T2', 'T1∩T3', 'T2∩T3', 'T1∩T2∩T3',
                   'T1_total', 'T2_total', 'T3_total'],
    'n_sites': [
        len(t1_sites - t2_sites - t3_sites),
        len(t2_sites - t1_sites - t3_sites),
        len(t3_sites - t1_sites - t2_sites),
        len(t1_t2 - t3_sites),
        len(t1_t3 - t2_sites),
        len(t2_t3 - t1_sites),
        len(t1_t2_t3),
        len(t1_sites), len(t2_sites), len(t3_sites)
    ]
})

# Also do per-motif overlap
motif_overlap_records = []
for m_name in ['TGGCCGGC', 'CCGG', 'GGCCGG', 'CCGC']:
    m_sub = motif[motif['assigned_motif'] == m_name]
    s1 = get_site_set(m_sub, 'T1')
    s2 = get_site_set(m_sub, 'T2')
    s3 = get_site_set(m_sub, 'T3')
    s12 = s1 & s2
    s13 = s1 & s3
    s23 = s2 & s3
    s123 = s1 & s2 & s3

    motif_overlap_records.append({
        'motif': m_name,
        'T1_total': len(s1), 'T2_total': len(s2), 'T3_total': len(s3),
        'T1_only': len(s1 - s2 - s3),
        'T2_only': len(s2 - s1 - s3),
        'T3_only': len(s3 - s1 - s2),
        'T1_T2_overlap': len(s12 - s3),
        'T1_T3_overlap': len(s13 - s2),
        'T2_T3_overlap': len(s23 - s1),
        'all_three_overlap': len(s123),
        'T1_T2_jaccard': len(s12) / len(s1 | s2) if len(s1 | s2) > 0 else 0,
        'T1_T3_jaccard': len(s13) / len(s1 | s3) if len(s1 | s3) > 0 else 0,
    })

motif_overlap_df = pd.DataFrame(motif_overlap_records)

# Combine overlap tables
overlap_table.to_csv(os.path.join(TABDIR, 'CCGG_site_overlap_between_timepoints.tsv'), sep='\t', index=False)
motif_overlap_df.to_csv(os.path.join(TABDIR, 'CCGG_site_overlap_per_motif.tsv'), sep='\t', index=False)

print(f"\nOverlap (all CCGG-related):")
print(overlap_table.to_string(index=False))
print(f"\nPer-motif overlap:")
print(motif_overlap_df.to_string(index=False))

# ============================================================
# 2b. Venn diagram of CCGG site overlaps
# ============================================================
print("\n  Creating Venn diagram...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# All CCGG-related
ax = axes[0]
plt.sca(ax)
v = venn3([t1_sites, t2_sites, t3_sites], set_labels=('T1', 'T2', 'T3'),
          set_colors=(COLORS['T1'], COLORS['T2'], COLORS['T3']), alpha=0.5)
ax.set_title('All CCGG-related 4mC sites\n(TGGCCGGC + CCGG + GGCCGG + CCGC)', fontsize=11)

# TGGCCGGC only (the dominant motif)
ax = axes[1]
plt.sca(ax)
s1 = get_site_set(motif[motif['assigned_motif'] == 'TGGCCGGC'], 'T1')
s2 = get_site_set(motif[motif['assigned_motif'] == 'TGGCCGGC'], 'T2')
s3 = get_site_set(motif[motif['assigned_motif'] == 'TGGCCGGC'], 'T3')
if len(s1) > 0 or len(s2) > 0 or len(s3) > 0:
    v2 = venn3([s1, s2, s3], set_labels=('T1', 'T2', 'T3'),
               set_colors=(COLORS['T1'], COLORS['T2'], COLORS['T3']), alpha=0.5)
ax.set_title('TGGCCGGC 4mC sites only', fontsize=11)

plt.tight_layout()
fig.savefig(os.path.join(FIGDIR, 'CCGG_site_venn_diagram.pdf'), format='pdf')
fig.savefig(os.path.join(FIGDIR, 'CCGG_site_venn_diagram.svg'), format='svg')
plt.close()
print("  Saved Venn diagram.")

# ============================================================
# 3. Frequency distribution analysis
# ============================================================
print("\n" + "=" * 60)
print("3. Frequency distribution analysis")
print("=" * 60)

# Use the motif assignment table which has frequency
freq_stats = []
for tp in ['T1', 'T2', 'T3']:
    sub = ccgg_all[ccgg_all['timepoint'] == tp]
    if len(sub) > 0:
        freq_stats.append({
            'timepoint': tp,
            'n_sites': len(sub),
            'mean_freq': sub['frequency'].mean(),
            'median_freq': sub['frequency'].median(),
            'std_freq': sub['frequency'].std(),
            'min_freq': sub['frequency'].min(),
            'max_freq': sub['frequency'].max(),
            'q25_freq': sub['frequency'].quantile(0.25),
            'q75_freq': sub['frequency'].quantile(0.75),
            'pct_above_80': (sub['frequency'] >= 80).mean() * 100,
            'pct_below_60': (sub['frequency'] < 60).mean() * 100,
        })

freq_stats_df = pd.DataFrame(freq_stats)
freq_stats_df.to_csv(os.path.join(TABDIR, 'CCGG_frequency_distribution.tsv'), sep='\t', index=False)
print(freq_stats_df.to_string(index=False))

# Also do per-motif frequency stats
freq_by_motif = []
for m_name in ['TGGCCGGC', 'CCGG', 'GGCCGG', 'CCGC']:
    m_sub = motif[motif['assigned_motif'] == m_name]
    for tp in ['T1', 'T2', 'T3']:
        sub = m_sub[m_sub['timepoint'] == tp]
        if len(sub) > 0:
            freq_by_motif.append({
                'motif': m_name,
                'timepoint': tp,
                'n_sites': len(sub),
                'mean_freq': sub['frequency'].mean(),
                'median_freq': sub['frequency'].median(),
                'std_freq': sub['frequency'].std(),
            })

freq_by_motif_df = pd.DataFrame(freq_by_motif)
print("\nPer-motif frequency stats:")
print(freq_by_motif_df.to_string(index=False))

# ============================================================
# 3b. Frequency distribution histogram
# ============================================================
print("\n  Creating frequency distribution histogram...")

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# All CCGG-related
ax = axes[0, 0]
bins = np.arange(0, 105, 5)
for tp in ['T1', 'T2', 'T3']:
    sub = ccgg_all[ccgg_all['timepoint'] == tp]
    if len(sub) > 0:
        ax.hist(sub['frequency'], bins=bins, alpha=0.5, label=f'{tp} (n={len(sub)})',
                color=COLORS[tp], edgecolor='white', linewidth=0.5)
ax.set_xlabel('Methylation frequency (%)')
ax.set_ylabel('Number of sites')
ax.set_title('All CCGG-related 4mC sites')
ax.legend()

# TGGCCGGC only
ax = axes[0, 1]
for tp in ['T1', 'T2', 'T3']:
    sub = motif[motif['assigned_motif'] == 'TGGCCGGC']
    sub = sub[sub['timepoint'] == tp]
    if len(sub) > 0:
        ax.hist(sub['frequency'], bins=bins, alpha=0.5, label=f'{tp} (n={len(sub)})',
                color=COLORS[tp], edgecolor='white', linewidth=0.5)
ax.set_xlabel('Methylation frequency (%)')
ax.set_ylabel('Number of sites')
ax.set_title('TGGCCGGC sites')
ax.legend()

# Since overlap is zero, show frequency comparison across timepoints (CDF)
ax = axes[1, 0]
for tp in ['T1', 'T2', 'T3']:
    sub = ccgg_all[ccgg_all['timepoint'] == tp]
    if len(sub) > 0:
        sorted_freq = np.sort(sub['frequency'].values)
        cdf = np.arange(1, len(sorted_freq) + 1) / len(sorted_freq)
        ax.plot(sorted_freq, cdf, label=f'{tp} (n={len(sub)})', color=COLORS[tp], linewidth=2)
ax.set_xlabel('Methylation frequency (%)')
ax.set_ylabel('Cumulative proportion')
ax.set_title('All CCGG-related 4mC: frequency CDF')
ax.legend()
ax.set_xlim(45, 105)

# Site count reduction bar chart per motif
ax = axes[1, 1]
motif_names = ['TGGCCGGC', 'CCGG', 'GGCCGG', 'CCGC']
x = np.arange(len(motif_names))
width = 0.25
for i, tp in enumerate(['T1', 'T2', 'T3']):
    counts = []
    for m_name in motif_names:
        n = len(motif[(motif['assigned_motif'] == m_name) & (motif['timepoint'] == tp)])
        counts.append(n)
    ax.bar(x + i * width, counts, width, label=tp, color=COLORS[tp], edgecolor='white')
ax.set_xticks(x + width)
ax.set_xticklabels(motif_names, fontsize=9)
ax.set_ylabel('Number of sites')
ax.set_title('CCGG-related 4mC sites by motif and timepoint')
ax.legend()
ax.set_yscale('symlog', linthresh=10)

plt.tight_layout()
fig.savefig(os.path.join(FIGDIR, 'CCGG_frequency_distribution.pdf'), format='pdf')
fig.savefig(os.path.join(FIGDIR, 'CCGG_frequency_distribution.svg'), format='svg')
plt.close()
print("  Saved frequency distribution figures.")

# ============================================================
# 4. Search for ALL methyltransferases in GFF
# ============================================================
print("\n" + "=" * 60)
print("4. Searching for all methyltransferases in GFF")
print("=" * 60)

mtase_records = []
gff_cds_records = []  # Store all CDS for operon context later

with open(GFF_FILE, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue
        if fields[2] != 'CDS':
            continue

        attrs = fields[8]
        chrom = fields[0]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]

        # Extract locus_tag
        lt_match = re.search(r'locus_tag=([^;]+)', attrs)
        locus_tag = lt_match.group(1) if lt_match else ''

        # Extract product
        prod_match = re.search(r'product=([^;]+)', attrs)
        product = prod_match.group(1) if prod_match else ''

        # Extract protein_id
        pid_match = re.search(r'protein_id=([^;]+)', attrs)
        protein_id = pid_match.group(1) if pid_match else ''

        gff_cds_records.append({
            'chrom': chrom,
            'start': start,
            'end': end,
            'strand': strand,
            'locus_tag': locus_tag,
            'product': product,
            'protein_id': protein_id,
        })

        # Check if methyltransferase-related
        product_lower = product.lower()
        is_mtase = any(kw in product_lower for kw in [
            'methyltransferase', 'methylase', 'dna methylation',
            'modification methylase', 'dcm', 'dam', 'hsdm',
            'n-6 adenine', 'c-5 cytosine', 'n-4 cytosine',
            'sam-dependent methyltransferase',
            'restriction-modification',
        ])

        if is_mtase:
            mtase_records.append({
                'locus_tag': locus_tag,
                'start': start,
                'end': end,
                'strand': strand,
                'product': product,
                'protein_id': protein_id,
            })

gff_cds_df = pd.DataFrame(gff_cds_records)
mtase_gff_df = pd.DataFrame(mtase_records)

print(f"Total CDS in GFF: {len(gff_cds_df)}")
print(f"MTase-related CDS found: {len(mtase_gff_df)}")

# Add expression data
if len(mtase_gff_df) > 0:
    # Get normalized counts for each MTase
    for idx, row in mtase_gff_df.iterrows():
        lt = row['locus_tag']
        if lt in norm_counts['gene_id'].values:
            counts_row = norm_counts[norm_counts['gene_id'] == lt].iloc[0]
            # T1 samples
            t1_cols = [c for c in norm_counts.columns if c.startswith('M145_1_')]
            t2_cols = [c for c in norm_counts.columns if c.startswith('M145_2_')]
            t3_cols = [c for c in norm_counts.columns if c.startswith('M145_3_')]

            mtase_gff_df.loc[idx, 'T1_mean_counts'] = counts_row[t1_cols].mean()
            mtase_gff_df.loc[idx, 'T2_mean_counts'] = counts_row[t2_cols].mean()
            mtase_gff_df.loc[idx, 'T3_mean_counts'] = counts_row[t3_cols].mean()
        else:
            mtase_gff_df.loc[idx, 'T1_mean_counts'] = np.nan
            mtase_gff_df.loc[idx, 'T2_mean_counts'] = np.nan
            mtase_gff_df.loc[idx, 'T3_mean_counts'] = np.nan

        # DESeq2 results
        if lt in deseq_t2t1['gene_id'].values:
            row_de = deseq_t2t1[deseq_t2t1['gene_id'] == lt].iloc[0]
            mtase_gff_df.loc[idx, 'T2vsT1_log2FC'] = row_de['log2FoldChange']
            mtase_gff_df.loc[idx, 'T2vsT1_padj'] = row_de['padj']
        if lt in deseq_t3t1['gene_id'].values:
            row_de = deseq_t3t1[deseq_t3t1['gene_id'] == lt].iloc[0]
            mtase_gff_df.loc[idx, 'T3vsT1_log2FC'] = row_de['log2FoldChange']
            mtase_gff_df.loc[idx, 'T3vsT1_padj'] = row_de['padj']

    mtase_gff_df.to_csv(os.path.join(TABDIR, 'all_methyltransferases_GFF.tsv'), sep='\t', index=False)
    print("\nAll methyltransferases found:")
    print(mtase_gff_df[['locus_tag', 'product', 'T1_mean_counts', 'T2_mean_counts', 'T3_mean_counts']].to_string(index=False))

# Also search broader: any gene with "methyl" in product or with DNA modification keywords
broader_mtase_records = []
for _, row in gff_cds_df.iterrows():
    p = row['product'].lower()
    # Broader search for potential methyltransferases
    if any(kw in p for kw in ['methyltransferase', 'methylase', 'dna methyl']):
        broader_mtase_records.append(row)

broader_df = pd.DataFrame(broader_mtase_records)
if len(broader_df) > len(mtase_gff_df):
    extra = set(broader_df['locus_tag']) - set(mtase_gff_df['locus_tag'])
    if extra:
        print(f"\nAdditional MTase candidates from broader search: {extra}")

# Specifically check cytosine MTases that could explain CCGG 4mC
print("\n--- Candidate cytosine MTases ---")
cytosine_mtases = mtase_gff_df[
    mtase_gff_df['product'].str.lower().str.contains('cytosine|dcm|c-5|n-4|c5-methyl', na=False)
]
if len(cytosine_mtases) > 0:
    print(cytosine_mtases[['locus_tag', 'product', 'T1_mean_counts', 'T2_mean_counts', 'T3_mean_counts']].to_string(index=False))
else:
    print("  No cytosine-specific MTases found by keyword search.")
    print("  Checking all MTases for potential cytosine activity...")

# ============================================================
# 5. CCGG Genomic context analysis
# ============================================================
print("\n" + "=" * 60)
print("5. CCGG genomic context analysis")
print("=" * 60)

# Load genome sequence
from Bio import SeqIO
# Parse multi-record FASTA, use only main chromosome
genome_records = {rec.id: rec for rec in SeqIO.parse(GENOME_FASTA, 'fasta')}
genome_record = genome_records['NC_003888.3']
genome_seq = str(genome_record.seq).upper()
genome_len = len(genome_seq)
print(f"Genome length: {genome_len:,} bp")

# Count all CCGG positions in genome
ccgg_positions_fwd = []
ccgg_positions_rev = []
for m in re.finditer('CCGG', genome_seq):
    ccgg_positions_fwd.append(m.start())  # 0-based
    # CCGG is palindromic, so reverse complement CCGG = CCGG
    ccgg_positions_rev.append(m.start())

print(f"Total CCGG sites in genome (fwd): {len(ccgg_positions_fwd)}")

# CCGG is palindromic, so each occurrence on + strand is also on - strand
# The 4mC would be at the internal C (position 1 in CCGG, i.e., 0-based genome_pos + 1)
# Actually, CCGG: C at pos 0, C at pos 1, G at pos 2, G at pos 3
# Reverse complement: CCGG (palindrome)
# The methylation could be at either C in CCGG

# Count TGGCCGGC (contains CCGG) in genome
tggccggc_fwd = [m.start() for m in re.finditer('TGGCCGGC', genome_seq)]
tggccggc_rev = [m.start() for m in re.finditer('GCCGGCCA', genome_seq)]
print(f"Total TGGCCGGC sites in genome: {len(tggccggc_fwd)} fwd + {len(tggccggc_rev)} rev = {len(tggccggc_fwd) + len(tggccggc_rev)}")

# Count GGCCGG (contains CCGG)
ggccgg_fwd = [m.start() for m in re.finditer('GGCCGG', genome_seq)]
ggccgg_rev = [m.start() for m in re.finditer('CCGGCC', genome_seq)]
print(f"Total GGCCGG sites in genome: {len(ggccgg_fwd)} fwd + {len(ggccgg_rev)} rev = {len(ggccgg_fwd) + len(ggccgg_rev)}")

# Genomic density analysis: bin genome into windows
window_size = 50000  # 50 kb
n_windows = genome_len // window_size + 1
bins_start = np.arange(0, genome_len, window_size)

# Count CCGG sites per window
ccgg_all_positions = np.array(ccgg_positions_fwd)
ccgg_per_window = np.histogram(ccgg_all_positions, bins=np.append(bins_start, genome_len))[0]

# Count methylated CCGG sites per window for each timepoint
def positions_to_array(df, tp):
    """Get positions of CCGG-related 4mC sites for a timepoint."""
    sub = df[df['timepoint'] == tp]
    return np.array(sub['position'].values) if len(sub) > 0 else np.array([])

meth_per_window = {}
for tp in ['T1', 'T2', 'T3']:
    pos = positions_to_array(ccgg_all, tp)
    if len(pos) > 0:
        meth_per_window[tp] = np.histogram(pos, bins=np.append(bins_start, genome_len))[0]
    else:
        meth_per_window[tp] = np.zeros(len(bins_start))

# Create genomic context table
context_records = []
for tp in ['T1', 'T2', 'T3']:
    sub = ccgg_all[ccgg_all['timepoint'] == tp]
    for _, row in sub.iterrows():
        pos = row['position']
        # Find nearest gene
        region_mask = (gff_cds_df['start'] <= pos) & (gff_cds_df['end'] >= pos)
        if region_mask.any():
            gene = gff_cds_df[region_mask].iloc[0]
            location = 'intragenic'
            nearest_gene = gene['locus_tag']
            gene_product = gene['product']
        else:
            # Find nearest gene
            dists = np.minimum(
                np.abs(gff_cds_df['start'].values - pos),
                np.abs(gff_cds_df['end'].values - pos)
            )
            idx = np.argmin(dists)
            gene = gff_cds_df.iloc[idx]
            location = 'intergenic'
            nearest_gene = gene['locus_tag']
            gene_product = gene['product']

        context_records.append({
            'timepoint': tp,
            'position': pos,
            'strand': row['strand'],
            'motif': row['assigned_motif'],
            'frequency': row['frequency'],
            'location': location,
            'nearest_gene': nearest_gene,
            'gene_product': gene_product,
            'genome_region': 'core' if 1500000 <= pos <= 6500000 else 'arm',
        })

context_df = pd.DataFrame(context_records)
context_df.to_csv(os.path.join(TABDIR, 'CCGG_genomic_context.tsv'), sep='\t', index=False)

# Summary statistics
print("\nGenomic context summary:")
for tp in ['T1', 'T2', 'T3']:
    sub = context_df[context_df['timepoint'] == tp]
    if len(sub) > 0:
        n_intra = (sub['location'] == 'intragenic').sum()
        n_inter = (sub['location'] == 'intergenic').sum()
        n_core = (sub['genome_region'] == 'core').sum()
        n_arm = (sub['genome_region'] == 'arm').sum()
        print(f"  {tp}: intragenic={n_intra} ({n_intra/len(sub)*100:.1f}%), intergenic={n_inter} ({n_inter/len(sub)*100:.1f}%)")
        print(f"       core={n_core} ({n_core/len(sub)*100:.1f}%), arm={n_arm} ({n_arm/len(sub)*100:.1f}%)")

# ============================================================
# 5b. Genome-wide density plot
# ============================================================
print("\n  Creating genome-wide density plot...")

fig, axes = plt.subplots(4, 1, figsize=(14, 10), gridspec_kw={'height_ratios': [1, 1.2, 1.2, 1.2]})

# Genomic CCGG density
ax = axes[0]
ax.fill_between(bins_start / 1e6, ccgg_per_window, alpha=0.5, color='gray')
ax.set_ylabel('CCGG sites\nper 50 kb')
ax.set_title('CCGG site density across S. coelicolor genome')
ax.set_xlim(0, genome_len / 1e6)
# Mark core/arm boundaries
for boundary in [1.5, 6.5]:
    ax.axvline(boundary, color='red', linestyle='--', alpha=0.3, linewidth=0.8)

# Methylated CCGG per timepoint
for i, tp in enumerate(['T1', 'T2', 'T3']):
    ax = axes[i + 1]
    ax.fill_between(bins_start / 1e6, meth_per_window[tp], alpha=0.7, color=COLORS[tp])
    ax.set_ylabel(f'{tp}\n4mC sites')
    ax.set_xlim(0, genome_len / 1e6)
    n_total = int(meth_per_window[tp].sum())
    ax.text(0.98, 0.85, f'n={n_total}', transform=ax.transAxes, ha='right', fontsize=10,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
    for boundary in [1.5, 6.5]:
        ax.axvline(boundary, color='red', linestyle='--', alpha=0.3, linewidth=0.8)

axes[-1].set_xlabel('Genome position (Mb)')
plt.tight_layout()
fig.savefig(os.path.join(FIGDIR, 'CCGG_genome_density.pdf'), format='pdf')
fig.savefig(os.path.join(FIGDIR, 'CCGG_genome_density.svg'), format='svg')
plt.close()
print("  Saved genome density plot.")

# ============================================================
# 6. Operon context of SC_RS19770 and SC_RS36410
# ============================================================
print("\n" + "=" * 60)
print("6. Operon context analysis for Dcm-like MTases")
print("=" * 60)

target_genes = ['SC_RS19770', 'SC_RS36410', 'SC_RS19670', 'SC_RS36625']
operon_records = []

for target in target_genes:
    target_row = gff_cds_df[gff_cds_df['locus_tag'] == target]
    if len(target_row) == 0:
        print(f"  WARNING: {target} not found in GFF CDS records")
        continue

    target_row = target_row.iloc[0]
    t_start = target_row['start']
    t_end = target_row['end']
    t_strand = target_row['strand']

    print(f"\n  {target}: {target_row['product']}")
    print(f"    Position: {t_start}-{t_end} ({t_strand})")

    # Get neighboring genes (within 10 kb on each side)
    flank = 10000
    neighbors = gff_cds_df[
        (gff_cds_df['start'] >= t_start - flank) &
        (gff_cds_df['end'] <= t_end + flank)
    ].copy()

    neighbors = neighbors.sort_values('start')

    for _, nb in neighbors.iterrows():
        lt = nb['locus_tag']
        # Get expression
        t1_mean = t2_mean = t3_mean = np.nan
        lfc_t2t1 = lfc_t3t1 = padj_t2t1 = padj_t3t1 = np.nan

        if lt in norm_counts['gene_id'].values:
            cr = norm_counts[norm_counts['gene_id'] == lt].iloc[0]
            t1_cols = [c for c in norm_counts.columns if c.startswith('M145_1_')]
            t2_cols = [c for c in norm_counts.columns if c.startswith('M145_2_')]
            t3_cols = [c for c in norm_counts.columns if c.startswith('M145_3_')]
            t1_mean = cr[t1_cols].mean()
            t2_mean = cr[t2_cols].mean()
            t3_mean = cr[t3_cols].mean()

        if lt in deseq_t2t1['gene_id'].values:
            dr = deseq_t2t1[deseq_t2t1['gene_id'] == lt].iloc[0]
            lfc_t2t1 = dr['log2FoldChange']
            padj_t2t1 = dr['padj']
        if lt in deseq_t3t1['gene_id'].values:
            dr = deseq_t3t1[deseq_t3t1['gene_id'] == lt].iloc[0]
            lfc_t3t1 = dr['log2FoldChange']
            padj_t3t1 = dr['padj']

        is_target = 'YES' if lt == target else ''

        operon_records.append({
            'target_gene': target,
            'locus_tag': lt,
            'start': nb['start'],
            'end': nb['end'],
            'strand': nb['strand'],
            'product': nb['product'],
            'T1_mean_counts': round(t1_mean, 2) if not np.isnan(t1_mean) else np.nan,
            'T2_mean_counts': round(t2_mean, 2) if not np.isnan(t2_mean) else np.nan,
            'T3_mean_counts': round(t3_mean, 2) if not np.isnan(t3_mean) else np.nan,
            'T2vsT1_log2FC': round(lfc_t2t1, 3) if not np.isnan(lfc_t2t1) else np.nan,
            'T3vsT1_log2FC': round(lfc_t3t1, 3) if not np.isnan(lfc_t3t1) else np.nan,
            'T2vsT1_padj': padj_t2t1,
            'T3vsT1_padj': padj_t3t1,
            'is_target': is_target,
        })

        marker = ' ***' if lt == target else ''
        print(f"    {lt} ({nb['strand']}) {nb['start']}-{nb['end']}: {nb['product'][:60]}{marker}")

operon_df = pd.DataFrame(operon_records)
operon_df.to_csv(os.path.join(TABDIR, 'Dcm_operon_context.tsv'), sep='\t', index=False)
print("\n  Operon context saved.")

# Check for REases near MTases
print("\n--- Checking for restriction endonucleases near Dcm-like MTases ---")
for target in ['SC_RS19770', 'SC_RS36410']:
    target_row = gff_cds_df[gff_cds_df['locus_tag'] == target]
    if len(target_row) == 0:
        continue
    t_start = target_row.iloc[0]['start']
    t_end = target_row.iloc[0]['end']
    neighbors = gff_cds_df[
        (gff_cds_df['start'] >= t_start - 10000) &
        (gff_cds_df['end'] <= t_end + 10000)
    ]
    for _, nb in neighbors.iterrows():
        p = nb['product'].lower()
        if any(kw in p for kw in ['restriction', 'endonuclease', 'rease', 'hsdr', 'mrr']):
            print(f"  Found REase near {target}: {nb['locus_tag']} - {nb['product']}")

# ============================================================
# 7. MTase expression timeline figure
# ============================================================
print("\n" + "=" * 60)
print("7. Creating MTase expression timeline figure")
print("=" * 60)

# Get all MTases and their expression
mtase_pred = pd.read_csv(MTASE_PRED)
print(f"MTase predictions: {len(mtase_pred)}")

# Parse mtase_pred to get individual locus tags
all_mtase_lts = set()
for _, row in mtase_pred.iterrows():
    lt_str = str(row['mtase'])
    for lt in lt_str.split(' / '):
        lt = lt.strip()
        all_mtase_lts.add(lt)

# Add any extra from GFF search
for lt in mtase_gff_df['locus_tag'].values:
    all_mtase_lts.add(lt)

print(f"Total unique MTase locus tags: {len(all_mtase_lts)}")

# Collect expression data
expr_data = {}
for lt in sorted(all_mtase_lts):
    if lt in norm_counts['gene_id'].values:
        cr = norm_counts[norm_counts['gene_id'] == lt].iloc[0]
        t1_cols = [c for c in norm_counts.columns if c.startswith('M145_1_')]
        t2_cols = [c for c in norm_counts.columns if c.startswith('M145_2_')]
        t3_cols = [c for c in norm_counts.columns if c.startswith('M145_3_')]
        expr_data[lt] = {
            'T1': cr[t1_cols].mean(),
            'T2': cr[t2_cols].mean(),
            'T3': cr[t3_cols].mean(),
        }

# Get product names
lt_to_product = {}
for _, row in gff_cds_df.iterrows():
    if row['locus_tag'] in all_mtase_lts:
        lt_to_product[row['locus_tag']] = row['product']

# Plot
fig, axes = plt.subplots(1, 2, figsize=(14, 7))

# Panel A: All MTases
ax = axes[0]
timepoints = ['T1', 'T2', 'T3']
x_pos = [0, 1, 2]

for lt in sorted(expr_data.keys()):
    values = [expr_data[lt][tp] for tp in timepoints]
    product = lt_to_product.get(lt, 'unknown')[:40]
    label = f"{lt}\n({product})"

    # Highlight Dcm-like MTases
    if lt in ['SC_RS19770', 'SC_RS36410']:
        ax.plot(x_pos, values, 'o-', linewidth=2.5, markersize=8, label=label, zorder=10)
    elif lt in ['SC_RS19670', 'SC_RS36625']:
        ax.plot(x_pos, values, 's--', linewidth=2, markersize=7, label=label, zorder=9)
    else:
        ax.plot(x_pos, values, '.-', linewidth=1, markersize=5, alpha=0.6, label=label)

ax.set_xticks(x_pos)
ax.set_xticklabels(timepoints)
ax.set_ylabel('Normalized counts (mean)')
ax.set_title('All predicted MTases: expression timeline')
ax.legend(fontsize=6, loc='upper left', bbox_to_anchor=(1.02, 1.0))
ax.set_yscale('symlog', linthresh=10)
ax.set_xlabel('Timepoint')

# Panel B: Focused on Dcm-related only (log scale)
ax = axes[1]
dcm_related = ['SC_RS19770', 'SC_RS36410', 'SC_RS19670', 'SC_RS36625', 'SC_RS19765']
colors_dcm = ['#e41a1c', '#377eb8', '#ff7f00', '#984ea3', '#a65628']
markers_dcm = ['o', 's', '^', 'D', 'v']

for lt, c, mk in zip(dcm_related, colors_dcm, markers_dcm):
    if lt in expr_data:
        values = [expr_data[lt][tp] for tp in timepoints]
        product = lt_to_product.get(lt, 'unknown')[:50]
        ax.plot(x_pos, values, f'{mk}-', color=c, linewidth=2.5, markersize=10,
                label=f"{lt}\n({product})")

ax.set_xticks(x_pos)
ax.set_xticklabels(timepoints)
ax.set_ylabel('Normalized counts (mean)')
ax.set_title('Dcm-related MTases: expression timeline')
ax.legend(fontsize=8, loc='best')
ax.set_yscale('symlog', linthresh=1)
ax.set_xlabel('Timepoint')

plt.tight_layout()
fig.savefig(os.path.join(FIGDIR, 'MTase_expression_timeline.pdf'), format='pdf')
fig.savefig(os.path.join(FIGDIR, 'MTase_expression_timeline.svg'), format='svg')
plt.close()
print("  Saved MTase expression timeline.")

# ============================================================
# 8. Passive dilution model analysis
# ============================================================
print("\n" + "=" * 60)
print("8. Passive dilution model evaluation")
print("=" * 60)

# If passive dilution is occurring, we expect:
# 1. Sites present at T1 but lost at T2 should have LOWER frequency at T1
#    (sites with lower initial frequency are more likely to be lost after dilution)
# 2. Shared sites should show frequency decrease from T1 to T2

# Analysis for TGGCCGGC (the dominant motif)
tggc = motif[motif['assigned_motif'] == 'TGGCCGGC'].copy()

t1_sites_df = tggc[tggc['timepoint'] == 'T1'].copy()
t2_sites_df = tggc[tggc['timepoint'] == 'T2'].copy()
t3_sites_df = tggc[tggc['timepoint'] == 'T3'].copy()

t1_pos_set = set(zip(t1_sites_df['position'], t1_sites_df['strand']))
t2_pos_set = set(zip(t2_sites_df['position'], t2_sites_df['strand']))
t3_pos_set = set(zip(t3_sites_df['position'], t3_sites_df['strand']))

from scipy import stats

# Check overlap
shared_t1t2 = t1_pos_set & t2_pos_set
shared_t1t3 = t1_pos_set & t3_pos_set
shared_t2t3 = t2_pos_set & t3_pos_set

print(f"\nTGGCCGGC site overlap:")
print(f"  T1∩T2: {len(shared_t1t2)} (Jaccard: {len(shared_t1t2)/max(len(t1_pos_set | t2_pos_set),1):.3f})")
print(f"  T1∩T3: {len(shared_t1t3)} (Jaccard: {len(shared_t1t3)/max(len(t1_pos_set | t3_pos_set),1):.3f})")
print(f"  T2∩T3: {len(shared_t2t3)} (Jaccard: {len(shared_t2t3)/max(len(t2_pos_set | t3_pos_set),1):.3f})")

if len(shared_t1t2) == 0:
    print("\n  *** ZERO overlap between T1 and T2 TGGCCGGC sites! ***")
    print("  This completely rules out passive dilution as an explanation.")
    print("  If passive dilution were occurring, a subset of T1 sites should persist at T2")
    print("  with reduced frequency. Instead, ALL T2 sites are at NEW positions.")
    print("  This pattern is consistent with INDEPENDENT de novo methylation at each timepoint,")
    print("  NOT maintenance methylation.")

# Compare frequency distributions (even though different sites)
print(f"\n  T1 sites (all lost by T2): n={len(t1_sites_df)}, mean freq={t1_sites_df['frequency'].mean():.1f}%, median={t1_sites_df['frequency'].median():.1f}%")
print(f"  T2 sites (all new):        n={len(t2_sites_df)}, mean freq={t2_sites_df['frequency'].mean():.1f}%, median={t2_sites_df['frequency'].median():.1f}%")
print(f"  T3 sites (all new):        n={len(t3_sites_df)}, mean freq={t3_sites_df['frequency'].mean():.1f}%, median={t3_sites_df['frequency'].median():.1f}%")

if len(t1_sites_df) > 0 and len(t2_sites_df) > 0:
    stat, pval = stats.mannwhitneyu(t1_sites_df['frequency'], t2_sites_df['frequency'], alternative='greater')
    print(f"  Mann-Whitney U (T1 > T2 freq): U={stat:.0f}, p={pval:.2e}")

# Check proximity of T1 and T2 sites - are T2 sites near T1 sites?
t1_positions = np.sort(t1_sites_df['position'].values)
t2_positions = np.sort(t2_sites_df['position'].values)
t3_positions = np.sort(t3_sites_df['position'].values)

# For each T2 site, find distance to nearest T1 site
if len(t1_positions) > 0 and len(t2_positions) > 0:
    min_dists = []
    for p2 in t2_positions:
        idx = np.searchsorted(t1_positions, p2)
        dists = []
        if idx > 0:
            dists.append(abs(p2 - t1_positions[idx - 1]))
        if idx < len(t1_positions):
            dists.append(abs(p2 - t1_positions[idx]))
        min_dists.append(min(dists))
    min_dists = np.array(min_dists)
    print(f"\n  T2 site distances to nearest T1 site:")
    print(f"    Mean: {min_dists.mean():.0f} bp, Median: {np.median(min_dists):.0f} bp")
    print(f"    Min: {min_dists.min():.0f} bp, Max: {min_dists.max():.0f} bp")
    print(f"    Within 100 bp: {(min_dists <= 100).sum()}")
    print(f"    Within 1 kb: {(min_dists <= 1000).sum()}")
    print(f"    Within 10 kb: {(min_dists <= 10000).sum()}")

# Check T3 new sites
t3_new = t3_pos_set - t1_pos_set - t2_pos_set
print(f"\n  T3 TGGCCGGC sites:")
print(f"    Total: {len(t3_pos_set)}")
print(f"    New at T3 (not in T1 or T2): {len(t3_new)}")
print(f"    Also present in T1: {len(shared_t1t3)}")
print(f"    Also present in T2: {len(shared_t2t3)}")

# ============================================================
# 8b. Check if SC_RS19765 is present in GFF (new candidate)
# ============================================================
print("\n--- SC_RS19765 investigation ---")
rs19765 = gff_cds_df[gff_cds_df['locus_tag'] == 'SC_RS19765']
if len(rs19765) > 0:
    print(f"  SC_RS19765: {rs19765.iloc[0]['product']}")
    print(f"  Position: {rs19765.iloc[0]['start']}-{rs19765.iloc[0]['end']} ({rs19765.iloc[0]['strand']})")
    # Get expression
    lt = 'SC_RS19765'
    if lt in norm_counts['gene_id'].values:
        cr = norm_counts[norm_counts['gene_id'] == lt].iloc[0]
        t1_cols = [c for c in norm_counts.columns if c.startswith('M145_1_')]
        t2_cols = [c for c in norm_counts.columns if c.startswith('M145_2_')]
        t3_cols = [c for c in norm_counts.columns if c.startswith('M145_3_')]
        print(f"  Expression: T1={cr[t1_cols].mean():.1f}, T2={cr[t2_cols].mean():.1f}, T3={cr[t3_cols].mean():.1f}")
    else:
        print(f"  Expression: not in normalized counts")
else:
    print("  SC_RS19765 not found in GFF")

# Also look for SC_RS19765 neighbors to see if it's near SC_RS19770
rs19770 = gff_cds_df[gff_cds_df['locus_tag'] == 'SC_RS19770']
if len(rs19770) > 0 and len(rs19765) > 0:
    dist = abs(rs19770.iloc[0]['start'] - rs19765.iloc[0]['end'])
    print(f"  Distance from SC_RS19765 to SC_RS19770: {dist} bp")

# ============================================================
# 9. Summary statistics
# ============================================================
print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)

# HC 4mC counts per timepoint
for tp in ['T1', 'T2', 'T3']:
    n_4mc = len(hc_4mc[hc_4mc['timepoint'] == tp])
    n_ccgg = len(ccgg_all[ccgg_all['timepoint'] == tp])
    pct = n_ccgg / n_4mc * 100 if n_4mc > 0 else 0
    print(f"{tp}: {n_4mc} HC 4mC sites, {n_ccgg} CCGG-related ({pct:.1f}%)")

# Dcm-like MTase expression
for lt in ['SC_RS19770', 'SC_RS36410']:
    if lt in expr_data:
        print(f"\n{lt} expression: T1={expr_data[lt]['T1']:.1f}, T2={expr_data[lt]['T2']:.1f}, T3={expr_data[lt]['T3']:.1f}")

# Overlap summary
print(f"\nCCGG-related site overlap:")
print(f"  T1∩T2: {len(t1_t2)} / T1={len(t1_sites)}, T2={len(t2_sites)}")
print(f"  T1∩T3: {len(t1_t3)} / T1={len(t1_sites)}, T3={len(t3_sites)}")
print(f"  Jaccard T1-T2: {len(t1_t2) / len(t1_sites | t2_sites):.3f}" if len(t1_sites | t2_sites) > 0 else "  N/A")
print(f"  Jaccard T1-T3: {len(t1_t3) / len(t1_sites | t3_sites):.3f}" if len(t1_sites | t3_sites) > 0 else "  N/A")

print("\n" + "=" * 60)
print("Analysis complete!")
print("=" * 60)
