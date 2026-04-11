#!/usr/bin/env python3
"""
H9: Systematic verification of whether CCGG "4mC" calls are actually 5mC
misclassified by Nanopore basecalling.

Hypothesis: The CCGG "4mC" detections are actually 5mC produced by Dcm-like
enzymes, misclassified by the Nanopore basecaller.

Evidence lines:
(a) Modification at C2 position (internal cytosine of CCGG = Dcm target)
(b) Systematic differences in detection scores between CCGG-context 4mC
    and other true 4mC sites
(c) CCWGG context enrichment (Dcm signature)
(d) Lower replicate consistency for CCGG 4mC vs other 4mC
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import stats
from collections import defaultdict, Counter
from pathlib import Path
import re
import os
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Configuration
# ============================================================
OUTDIR = '/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/32_CCGG_5mC_misclassification'
FIGDIR = os.path.join(OUTDIR, 'figures')
TABDIR = os.path.join(OUTDIR, 'tables')

HC_SITES = '/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/01_integration/high_confidence_sites_weighted.csv'
MOTIF_ASSIGN = '/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/23_expanded_motif_search/4mC_motif_assignment.csv'
GENOME_FASTA = '/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/GCF_000203835.1_ASM20383v1_genomic.fna'
PILEUP_DIR = '/Users/okaban/bioinfo/methyl/260102_M145/analysis/cursor_results/20260108/pileup'

# Sample mapping
SAMPLE_TIMEPOINT = {
    '1-1': 'T1', '1-2': 'T1', '1-3': 'T1',
    '2-1': 'T2', '2-3': 'T2', '2-4': 'T2',
    '3-2': 'T3', '3-3': 'T3', '3-4': 'T3'
}

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

COLORS = {
    'T1': '#2166ac', 'T2': '#67a9cf', 'T3': '#ef8a62',
    'CCGG_4mC': '#e41a1c', 'non_CCGG_4mC': '#377eb8',
    '6mA': '#4daf4a', 'AAGCCCG_4mC': '#984ea3',
    'C1': '#ff7f00', 'C2': '#984ea3',
}

# ============================================================
# 1. Load data
# ============================================================
print("=" * 70)
print("H9: CCGG 5mC Misclassification Analysis")
print("=" * 70)

print("\n[1] Loading data...")

hc = pd.read_csv(HC_SITES)
motif = pd.read_csv(MOTIF_ASSIGN)

print(f"  HC sites total: {len(hc)}")
print(f"  4mC motif assignments: {len(motif)}")

# Load genome
from Bio import SeqIO
genome_records = {rec.id: rec for rec in SeqIO.parse(GENOME_FASTA, 'fasta')}
genome_record = genome_records['NC_003888.3']
genome_seq = str(genome_record.seq).upper()
genome_len = len(genome_seq)
print(f"  Genome length: {genome_len:,} bp")

# Separate modification types
hc_4mc = hc[hc['mod_type'] == '4mC'].copy()
hc_6ma = hc[hc['mod_type'] == '6mA'].copy()
hc_5mc = hc[hc['mod_type'] == '5mC'].copy() if '5mC' in hc['mod_type'].values else pd.DataFrame()

print(f"  HC 4mC: {len(hc_4mc)}, HC 6mA: {len(hc_6ma)}, HC 5mC: {len(hc_5mc)}")

# CCGG-related motifs
ccgg_motifs = ['CCGG', 'GGCCGG', 'TGGCCGGC', 'CCGC']
motif_ccgg = motif[motif['assigned_motif'].isin(ccgg_motifs)].copy()
motif_non_ccgg = motif[~motif['assigned_motif'].isin(ccgg_motifs + ['unassigned'])].copy()
motif_unassigned = motif[motif['assigned_motif'] == 'unassigned'].copy()

print(f"  CCGG-related 4mC motif sites: {len(motif_ccgg)}")
print(f"  Non-CCGG 4mC sites (assigned): {len(motif_non_ccgg)}")
print(f"  Unassigned 4mC sites: {len(motif_unassigned)}")

# ============================================================
# 2. STEP 1 & 2: C position within CCGG analysis
# ============================================================
print("\n" + "=" * 70)
print("[2] C position within CCGG motif analysis")
print("=" * 70)

def find_ccgg_position(position, strand, genome_seq):
    """
    For a given methylation site position, determine which C in CCGG it is.

    CCGG on + strand: C1=pos0, C2=pos1, G=pos2, G=pos3
    CCGG on - strand (complement CCGG): positions are reversed

    Returns: 'C1', 'C2', 'not_in_CCGG', or 'ambiguous'
    Also returns the CCGG start position.
    """
    # The position in the HC sites is 0-based genome coordinate
    # Check if this position falls within a CCGG on either strand

    results = []

    # Check + strand CCGG occurrences near this position
    for offset in range(-3, 1):  # CCGG could start at pos-3 to pos
        start = position + offset
        if 0 <= start <= len(genome_seq) - 4:
            context = genome_seq[start:start+4]
            if context == 'CCGG':
                rel_pos = position - start
                if strand == '+':
                    if rel_pos == 0:
                        results.append(('C1', start, '+'))
                    elif rel_pos == 1:
                        results.append(('C2', start, '+'))
                elif strand == '-':
                    # On - strand, the CCGG complement is read as CCGG
                    # Position 3 = C1 of reverse complement
                    # Position 2 = C2 of reverse complement
                    if rel_pos == 3:
                        results.append(('C1', start, '-'))
                    elif rel_pos == 2:
                        results.append(('C2', start, '-'))

    # Also check reverse complement CCGG on - strand
    # Reverse complement of CCGG is CCGG (palindrome)
    # So same search applies

    if len(results) == 1:
        return results[0][0], results[0][1], results[0][2]
    elif len(results) > 1:
        return 'ambiguous', results[0][1], results[0][2]
    else:
        return 'not_in_CCGG', -1, ''


# Analyze all CCGG-related 4mC sites
ccgg_position_records = []
for _, row in motif_ccgg.iterrows():
    pos = row['position']
    strand = row['strand']
    seq_context = row['sequence']

    c_pos, ccgg_start, ccgg_strand = find_ccgg_position(pos, strand, genome_seq)

    # Also extract the broader context (5bp on each side of CCGG)
    if ccgg_start >= 0 and ccgg_start + 4 <= len(genome_seq):
        broad_start = max(0, ccgg_start - 5)
        broad_end = min(len(genome_seq), ccgg_start + 9)
        broad_context = genome_seq[broad_start:broad_end]
    else:
        broad_context = ''

    ccgg_position_records.append({
        'position': pos,
        'strand': strand,
        'timepoint': row['timepoint'],
        'motif': row['assigned_motif'],
        'frequency': row['frequency'],
        'C_position_in_CCGG': c_pos,
        'ccgg_start': ccgg_start,
        'broad_context': broad_context,
        'sequence': seq_context,
    })

ccgg_pos_df = pd.DataFrame(ccgg_position_records)

# Summary
print("\nC position within CCGG distribution:")
pos_counts = ccgg_pos_df['C_position_in_CCGG'].value_counts()
for p, n in pos_counts.items():
    print(f"  {p}: {n} ({n/len(ccgg_pos_df)*100:.1f}%)")

# By timepoint
print("\nC position by timepoint:")
for tp in ['T1', 'T2', 'T3']:
    sub = ccgg_pos_df[ccgg_pos_df['timepoint'] == tp]
    if len(sub) > 0:
        pc = sub['C_position_in_CCGG'].value_counts()
        print(f"  {tp} (n={len(sub)}):")
        for p in ['C1', 'C2', 'ambiguous', 'not_in_CCGG']:
            n = pc.get(p, 0)
            print(f"    {p}: {n} ({n/len(sub)*100:.1f}%)")

# By motif
print("\nC position by assigned motif:")
for m in ['TGGCCGGC', 'CCGG', 'GGCCGG', 'CCGC']:
    sub = ccgg_pos_df[ccgg_pos_df['motif'] == m]
    if len(sub) > 0:
        pc = sub['C_position_in_CCGG'].value_counts()
        print(f"  {m} (n={len(sub)}):")
        for p in ['C1', 'C2', 'ambiguous', 'not_in_CCGG']:
            n = pc.get(p, 0)
            print(f"    {p}: {n} ({n/len(sub)*100:.1f}%)")

# Save position table
ccgg_pos_df.to_csv(os.path.join(TABDIR, 'CCGG_C_position_analysis.tsv'), sep='\t', index=False)

# ============================================================
# 3. STEP 3: Compare detection scores from raw pileup data
# ============================================================
print("\n" + "=" * 70)
print("[3] Detection score comparison from raw pileup data")
print("=" * 70)

def parse_pileup_for_sites(filepath, target_sites):
    """
    Parse a pileup file and extract data for specific target sites.
    target_sites: set of (chrom, pos, strand, mod_type) tuples
    Returns dict: site_key -> {coverage, mod_freq, n_modified, n_canonical, ...}
    """
    results = {}
    mod_type_map = {'a': '6mA', 'm': '5mC', '21839': '4mC'}

    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 17:
                continue

            chrom = parts[0]
            pos = int(parts[1])  # 0-based
            mod_code = parts[3]
            strand = parts[5]
            coverage = int(parts[9])
            mod_freq = float(parts[10])

            if mod_code not in mod_type_map:
                continue

            mod_type = mod_type_map[mod_code]
            key = (chrom, pos, strand, mod_type)

            if key in target_sites:
                # Parse detailed counts
                n_modified = int(parts[11])
                n_canonical = int(parts[12])
                n_other_mod = int(parts[13])
                n_delete = int(parts[14])
                n_fail = int(parts[15])
                n_diff = int(parts[16])
                n_nocall = int(parts[17]) if len(parts) > 17 else 0

                results[key] = {
                    'coverage': coverage,
                    'mod_freq': mod_freq,
                    'n_modified': n_modified,
                    'n_canonical': n_canonical,
                    'n_other_mod': n_other_mod,
                    'n_delete': n_delete,
                    'n_fail': n_fail,
                    'n_diff': n_diff,
                    'n_nocall': n_nocall,
                }

    return results


# Build target site sets
# Group A: CCGG-context 4mC
ccgg_4mc_sites = set()
for _, row in motif_ccgg.iterrows():
    ccgg_4mc_sites.add(('NC_003888.3', row['position'], row['strand'], '4mC'))

# Group B: Non-CCGG 4mC (true 4mC controls - CCGCGG motif)
non_ccgg_4mc_sites = set()
for _, row in motif_non_ccgg.iterrows():
    non_ccgg_4mc_sites.add(('NC_003888.3', row['position'], row['strand'], '4mC'))

# Group C: 6mA sites
sixma_sites = set()
for _, row in hc_6ma.iterrows():
    sixma_sites.add(('NC_003888.3', row['position'], row['strand'], '6mA'))

# Also get 5mC calls at the SAME positions as CCGG 4mC
# This tells us if modkit also detected 5mC at these positions
ccgg_5mc_check = set()
for _, row in motif_ccgg.iterrows():
    ccgg_5mc_check.add(('NC_003888.3', row['position'], row['strand'], '5mC'))

all_target_sites = ccgg_4mc_sites | non_ccgg_4mc_sites | sixma_sites | ccgg_5mc_check

print(f"  Target sites: CCGG 4mC={len(ccgg_4mc_sites)}, non-CCGG 4mC={len(non_ccgg_4mc_sites)}, 6mA={len(sixma_sites)}")
print(f"  Also checking 5mC at CCGG positions: {len(ccgg_5mc_check)}")

# Parse all pileup files
pileup_files = sorted(Path(PILEUP_DIR).glob('*_pileup.bed'))
# Exclude _v1 files
pileup_files = [f for f in pileup_files if '_v1' not in f.name]

print(f"\n  Parsing {len(pileup_files)} pileup files...")

sample_data = {}
for pf in pileup_files:
    sample_id = pf.name.replace('_pileup.bed', '')
    if sample_id not in SAMPLE_TIMEPOINT:
        continue

    tp = SAMPLE_TIMEPOINT[sample_id]
    print(f"    {sample_id} ({tp})...")
    data = parse_pileup_for_sites(str(pf), all_target_sites)
    sample_data[sample_id] = {'timepoint': tp, 'data': data}
    print(f"      Found {len(data)} target sites")

# Aggregate per-replicate data for score comparison
print("\n  Aggregating per-replicate scores...")

score_records = []
for sample_id, sdata in sample_data.items():
    tp = sdata['timepoint']
    for site_key, site_data in sdata['data'].items():
        chrom, pos, strand, mod_type = site_key

        # Determine group
        if site_key in ccgg_4mc_sites:
            group = 'CCGG_4mC'
        elif site_key in non_ccgg_4mc_sites:
            group = 'non_CCGG_4mC'
        elif mod_type == '6mA' and site_key in sixma_sites:
            group = '6mA'
        elif site_key in ccgg_5mc_check:
            group = 'CCGG_pos_5mC'
        else:
            continue

        score_records.append({
            'sample': sample_id,
            'timepoint': tp,
            'chrom': chrom,
            'position': pos,
            'strand': strand,
            'mod_type': mod_type,
            'group': group,
            'coverage': site_data['coverage'],
            'mod_freq': site_data['mod_freq'],
            'n_modified': site_data['n_modified'],
            'n_canonical': site_data['n_canonical'],
            'n_other_mod': site_data['n_other_mod'],
            'n_delete': site_data['n_delete'],
            'n_fail': site_data['n_fail'],
            'n_diff': site_data['n_diff'],
        })

score_df = pd.DataFrame(score_records)
print(f"\n  Total replicate-level records: {len(score_df)}")
print(f"  Group breakdown:")
for g, n in score_df['group'].value_counts().items():
    print(f"    {g}: {n}")

# Save raw score data
score_df.to_csv(os.path.join(TABDIR, 'replicate_level_scores.tsv'), sep='\t', index=False)

# ============================================================
# 3b. Statistical comparison of detection frequencies
# ============================================================
print("\n  Statistical comparison of detection frequencies...")

# Compare per-replicate modification frequencies
# For fair comparison, use sites that are detected (mod_freq > 0)
# and with adequate coverage (>= 5)
min_cov = 5

groups_to_compare = ['CCGG_4mC', 'non_CCGG_4mC', '6mA']
for tp in ['T1', 'T2', 'T3']:
    tp_data = score_df[(score_df['timepoint'] == tp) & (score_df['coverage'] >= min_cov)]
    print(f"\n  --- {tp} ---")
    for g in groups_to_compare:
        g_data = tp_data[tp_data['group'] == g]
        if len(g_data) > 0:
            print(f"  {g}: n={len(g_data)}, mean_freq={g_data['mod_freq'].mean():.1f}%, "
                  f"median={g_data['mod_freq'].median():.1f}%, "
                  f"mean_cov={g_data['coverage'].mean():.1f}")

# Pairwise comparisons
comparison_stats = []

for tp in ['T1', 'T2', 'T3']:
    tp_data = score_df[(score_df['timepoint'] == tp) & (score_df['coverage'] >= min_cov)]

    ccgg_freq = tp_data[tp_data['group'] == 'CCGG_4mC']['mod_freq'].values
    non_ccgg_freq = tp_data[tp_data['group'] == 'non_CCGG_4mC']['mod_freq'].values
    sixma_freq = tp_data[tp_data['group'] == '6mA']['mod_freq'].values

    # CCGG 4mC vs non-CCGG 4mC
    if len(ccgg_freq) > 5 and len(non_ccgg_freq) > 5:
        u_stat, u_pval = stats.mannwhitneyu(ccgg_freq, non_ccgg_freq, alternative='two-sided')
        ks_stat, ks_pval = stats.ks_2samp(ccgg_freq, non_ccgg_freq)
        comparison_stats.append({
            'timepoint': tp,
            'comparison': 'CCGG_4mC vs non_CCGG_4mC',
            'n1': len(ccgg_freq), 'n2': len(non_ccgg_freq),
            'mean1': ccgg_freq.mean(), 'mean2': non_ccgg_freq.mean(),
            'median1': np.median(ccgg_freq), 'median2': np.median(non_ccgg_freq),
            'MannWhitney_U': u_stat, 'MannWhitney_p': u_pval,
            'KS_stat': ks_stat, 'KS_p': ks_pval,
        })

    # CCGG 4mC vs 6mA
    if len(ccgg_freq) > 5 and len(sixma_freq) > 5:
        u_stat, u_pval = stats.mannwhitneyu(ccgg_freq, sixma_freq, alternative='two-sided')
        ks_stat, ks_pval = stats.ks_2samp(ccgg_freq, sixma_freq)
        comparison_stats.append({
            'timepoint': tp,
            'comparison': 'CCGG_4mC vs 6mA',
            'n1': len(ccgg_freq), 'n2': len(sixma_freq),
            'mean1': ccgg_freq.mean(), 'mean2': sixma_freq.mean(),
            'median1': np.median(ccgg_freq), 'median2': np.median(sixma_freq),
            'MannWhitney_U': u_stat, 'MannWhitney_p': u_pval,
            'KS_stat': ks_stat, 'KS_p': ks_pval,
        })

comparison_stats_df = pd.DataFrame(comparison_stats)
comparison_stats_df.to_csv(os.path.join(TABDIR, 'score_comparison_statistics.tsv'), sep='\t', index=False)

print("\n  Score comparison statistics:")
print(comparison_stats_df.to_string(index=False))

# ============================================================
# 3c. Check 5mC signal at CCGG 4mC positions
# ============================================================
print("\n\n  Checking 5mC signal at CCGG 4mC positions...")

ccgg_5mc_data = score_df[score_df['group'] == 'CCGG_pos_5mC'].copy()
ccgg_4mc_data = score_df[score_df['group'] == 'CCGG_4mC'].copy()

print(f"  CCGG positions with 5mC signal: {len(ccgg_5mc_data)} records")
print(f"  CCGG positions with 4mC signal: {len(ccgg_4mc_data)} records")

if len(ccgg_5mc_data) > 0:
    # Sites where both 5mC and 4mC are detected
    ccgg_5mc_positions = set(zip(ccgg_5mc_data['position'], ccgg_5mc_data['strand']))
    ccgg_4mc_positions = set(zip(ccgg_4mc_data['position'], ccgg_4mc_data['strand']))
    both = ccgg_5mc_positions & ccgg_4mc_positions

    print(f"  Positions with both 5mC and 4mC signal: {len(both)}")

    # Compare frequencies at dual-signal sites
    if len(both) > 0:
        dual_records = []
        for pos, strand in both:
            fmc_data = ccgg_4mc_data[(ccgg_4mc_data['position'] == pos) & (ccgg_4mc_data['strand'] == strand)]
            fivemc_data = ccgg_5mc_data[(ccgg_5mc_data['position'] == pos) & (ccgg_5mc_data['strand'] == strand)]

            for _, r4 in fmc_data.iterrows():
                # Find matching 5mC record from same sample
                match_5 = fivemc_data[fivemc_data['sample'] == r4['sample']]
                if len(match_5) > 0:
                    dual_records.append({
                        'position': pos,
                        'strand': strand,
                        'sample': r4['sample'],
                        'timepoint': r4['timepoint'],
                        '4mC_freq': r4['mod_freq'],
                        '5mC_freq': match_5.iloc[0]['mod_freq'],
                        '4mC_coverage': r4['coverage'],
                        '5mC_coverage': match_5.iloc[0]['coverage'],
                    })

        dual_df = pd.DataFrame(dual_records)
        if len(dual_df) > 0:
            print(f"\n  Dual-signal site comparison ({len(dual_df)} records):")
            print(f"    Mean 4mC freq: {dual_df['4mC_freq'].mean():.1f}%")
            print(f"    Mean 5mC freq: {dual_df['5mC_freq'].mean():.1f}%")
            dual_df.to_csv(os.path.join(TABDIR, 'dual_4mC_5mC_signal.tsv'), sep='\t', index=False)

# ============================================================
# 4. STEP 4: AAGCCCG dual-modification comparison
# ============================================================
print("\n" + "=" * 70)
print("[4] AAGCCCG dual-modification comparison")
print("=" * 70)

# AAGCCCG has both 6mA (at position A1) and 4mC (at the C within the GG-CC-GG overlap)
# Find AAGCCCG 4mC sites from motif assignment
# AAGCCCG contains the CCGG subsequence? Let's check: AAGCCCG -> no direct CCGG
# But within AAGCCCG context, the 4mC is at the internal C of a CC motif
# The GGCCGG motif overlaps with AAGCCCG

# From the expanded motif search, AAGCCCG was identified as a dual-modification motif
# Let's check which 4mC sites are in AAGCCCG context

aagcccg_4mc = []
for _, row in motif.iterrows():
    seq = row['sequence']
    if 'AAGCCCG' in seq or 'CGGGCTT' in seq:  # Forward or reverse complement
        aagcccg_4mc.append(row)

aagcccg_4mc_df = pd.DataFrame(aagcccg_4mc)
print(f"  4mC sites in AAGCCCG context: {len(aagcccg_4mc_df)}")

if len(aagcccg_4mc_df) > 0:
    print(f"  Motif assignments within AAGCCCG context:")
    print(aagcccg_4mc_df['assigned_motif'].value_counts().to_string())

    # Compare frequency distributions
    aagcccg_freq = aagcccg_4mc_df['frequency'].values
    standalone_ccgg_freq = motif_ccgg[~motif_ccgg.index.isin(aagcccg_4mc_df.index)]['frequency'].values

    print(f"\n  AAGCCCG-context 4mC: n={len(aagcccg_freq)}, mean={aagcccg_freq.mean():.1f}%, median={np.median(aagcccg_freq):.1f}%")
    print(f"  Standalone CCGG 4mC: n={len(standalone_ccgg_freq)}, mean={standalone_ccgg_freq.mean():.1f}%, median={np.median(standalone_ccgg_freq):.1f}%")

    if len(aagcccg_freq) > 5 and len(standalone_ccgg_freq) > 5:
        u, p = stats.mannwhitneyu(aagcccg_freq, standalone_ccgg_freq, alternative='two-sided')
        print(f"  Mann-Whitney U: U={u:.0f}, p={p:.2e}")

# ============================================================
# 5. STEP 5: CCWGG context analysis (Dcm signature)
# ============================================================
print("\n" + "=" * 70)
print("[5] CCWGG context analysis (Dcm specificity)")
print("=" * 70)

# E. coli Dcm methylates the internal C of CC(A/T)GG = CCWGG
# Check what proportion of CCGG 4mC sites are in CCWGG context

ccwgg_records = []
for _, row in motif_ccgg.iterrows():
    pos = row['position']
    strand = row['strand']

    # Get 5bp context centered on the nearest CCGG
    # Find nearest CCGG
    for offset in range(-4, 2):
        start = pos + offset
        if 0 <= start <= len(genome_seq) - 4:
            if genome_seq[start:start+4] == 'CCGG':
                # Check if this is part of CCWGG (where W = A or T)
                # CCWGG: CC at positions 0-1, W at position 2, GG at 3-4
                # But our match is CCGG (4bp). Check if there's a W before the GG
                # Wait - CCWGG is 5bp: C-C-W-G-G
                # CCGG is 4bp: C-C-G-G
                # CCWGG != CCGG; CCWGG has an A/T between CC and GG
                # So CCGG sites are NOT in CCWGG context by definition
                #
                # Actually, let me re-read: CCWGG = CC(A/T)GG
                # CCGG = C-C-G-G (no W in between)
                # These are DIFFERENT motifs!
                #
                # Dcm targets CCWGG (5bp), not CCGG (4bp)
                #
                # But the question asks: of the sites called as "CCGG" context 4mC,
                # how many are actually in a broader CCWGG context?
                # CCGG is a subset if we consider extended context
                # Let me check the broader context for each site

                # Get extended context: 2bp before CCGG start, CCGG, 2bp after
                ext_start = max(0, start - 1)
                ext_end = min(len(genome_seq), start + 5)
                ext_context = genome_seq[ext_start:ext_end]

                # Check for CCWGG pattern
                # Pattern: CC[AT]GG
                is_ccwgg = False
                is_ccagg = False
                is_cctgg = False
                is_ccggg = False
                is_cccgg = False

                # The CCGG is at positions start:start+4
                # Check 5bp starting at start: CC_GG where _ is the 3rd base
                if start + 4 < len(genome_seq):
                    five_mer = genome_seq[start:start+5]
                    # Wait, CCGG is only 4 bases
                    # I need to check if the CCGG is part of a 5bp context like CC_GG
                    pass

                # Actually, CCGG (4bp) and CCWGG (5bp) are completely different motifs.
                # CCWGG = CC-A/T-GG  (5bp)
                # CCGG = CC-GG (4bp)
                # They share CC...GG but CCWGG has an extra base in between

                # So the real question is: are the Dcm targets in S. coelicolor
                # actually CCGG or CCWGG?
                # Let's just check the extended contexts

                # Get the extended context around the methylated C
                ext5 = genome_seq[max(0,pos-5):min(len(genome_seq),pos+6)]

                ccwgg_records.append({
                    'position': pos,
                    'strand': strand,
                    'timepoint': row['timepoint'],
                    'motif': row['assigned_motif'],
                    'frequency': row['frequency'],
                    'ccgg_at': start,
                    'extended_context': ext5,
                })
                break

ccwgg_df = pd.DataFrame(ccwgg_records)

# Now check for CCWGG pattern in the genome near CCGG sites
# Actually, CCGG and CCWGG are distinct. Let me instead check:
# 1. How many CCWGG sites exist in the genome
# 2. How many CCGG sites exist
# 3. Do any 4mC sites fall at CCWGG rather than CCGG

ccwgg_fwd = [m.start() for m in re.finditer('CC[AT]GG', genome_seq)]
ccgg_fwd = [m.start() for m in re.finditer('CCGG', genome_seq)]

print(f"  Genome-wide motif counts:")
print(f"    CCGG sites: {len(ccgg_fwd)}")
print(f"    CCWGG (CCAGG+CCTGG) sites: {len(ccwgg_fwd)}")
print(f"    CCAGG sites: {len([m.start() for m in re.finditer('CCAGG', genome_seq)])}")
print(f"    CCTGG sites: {len([m.start() for m in re.finditer('CCTGG', genome_seq)])}")

# Check if any of the "CCGG" 4mC sites are actually at CCWGG positions
# This would require the methylated C to be near a CCWGG
ccwgg_positions_fwd = set(ccwgg_fwd)
ccwgg_positions_all = set()
for p in ccwgg_fwd:
    ccwgg_positions_all.add(p)       # C1 of CCWGG
    ccwgg_positions_all.add(p + 1)   # C2 of CCWGG

# Check overlap
n_at_ccwgg = 0
n_at_ccgg_only = 0
n_neither = 0

for _, row in motif_ccgg.iterrows():
    pos = row['position']
    # Check if pos is within a CCWGG context
    in_ccwgg = False
    for offset in range(-4, 2):
        s = pos + offset
        if s in ccwgg_positions_fwd:
            # pos is within a CCWGG site
            in_ccwgg = True
            break

    if in_ccwgg:
        n_at_ccwgg += 1
    else:
        # Check if at CCGG
        in_ccgg = False
        for offset in range(-3, 1):
            s = pos + offset
            if 0 <= s <= len(genome_seq) - 4:
                if genome_seq[s:s+4] == 'CCGG':
                    in_ccgg = True
                    break
        if in_ccgg:
            n_at_ccgg_only += 1
        else:
            n_neither += 1

print(f"\n  CCGG 4mC sites context breakdown:")
print(f"    At CCWGG (Dcm target): {n_at_ccwgg} ({n_at_ccwgg/len(motif_ccgg)*100:.1f}%)")
print(f"    At CCGG only (not CCWGG): {n_at_ccgg_only} ({n_at_ccgg_only/len(motif_ccgg)*100:.1f}%)")
print(f"    Neither: {n_neither} ({n_neither/len(motif_ccgg)*100:.1f}%)")

# Deeper analysis: for each CCGG 4mC site, what is the 5bp context?
# Look at base at position +2 relative to CCGG start (the base between CC and GG)
# Wait - CCGG doesn't have a base between CC and GG. It's CC immediately followed by GG.
# The question should be: what flanking bases surround CCGG?
# Specifically, the base BEFORE CCGG (position -1) and after CCGG (position +4)

flanking_analysis = []
for _, row in ccwgg_df.iterrows():
    ccgg_at = row['ccgg_at']
    if ccgg_at > 0 and ccgg_at + 4 < len(genome_seq):
        base_before = genome_seq[ccgg_at - 1]
        base_after = genome_seq[ccgg_at + 4]
        five_bp_context = genome_seq[ccgg_at-1:ccgg_at+5]  # [-1]CCGG[+4]
        six_bp_context = genome_seq[max(0,ccgg_at-1):min(len(genome_seq),ccgg_at+5)]

        flanking_analysis.append({
            'position': row['position'],
            'strand': row['strand'],
            'timepoint': row['timepoint'],
            'motif': row['motif'],
            'frequency': row['frequency'],
            'base_before': base_before,
            'base_after': base_after,
            'six_bp_context': six_bp_context,
        })

flanking_df = pd.DataFrame(flanking_analysis)

if len(flanking_df) > 0:
    print(f"\n  Flanking base analysis around CCGG:")
    print(f"  Base before CCGG:")
    for b, n in flanking_df['base_before'].value_counts().items():
        print(f"    {b}: {n} ({n/len(flanking_df)*100:.1f}%)")

    print(f"  Base after CCGG:")
    for b, n in flanking_df['base_after'].value_counts().items():
        print(f"    {b}: {n} ({n/len(flanking_df)*100:.1f}%)")

    # 6-bp context counts
    print(f"\n  Top 6-bp contexts (nCCGGn):")
    ctx_counts = flanking_df['six_bp_context'].value_counts()
    for ctx, n in ctx_counts.head(20).items():
        print(f"    {ctx}: {n} ({n/len(flanking_df)*100:.1f}%)")

    flanking_df.to_csv(os.path.join(TABDIR, 'CCGG_flanking_context.tsv'), sep='\t', index=False)

# ============================================================
# 6. STEP 6: Replicate consistency comparison
# ============================================================
print("\n" + "=" * 70)
print("[6] Replicate consistency comparison")
print("=" * 70)

# The high-confidence sites have n_reps column
# Compare n_reps distribution for CCGG 4mC vs other 4mC vs 6mA

# First, merge HC sites with motif assignments to get CCGG annotation
hc_4mc_annotated = hc_4mc.copy()

# Create a lookup from motif assignments
motif_lookup = {}
for _, row in motif.iterrows():
    key = (row['position'], row['strand'], row['timepoint'])
    motif_lookup[key] = row['assigned_motif']

hc_4mc_annotated['motif_group'] = 'other'
for idx, row in hc_4mc_annotated.iterrows():
    key = (row['position'], row['strand'], row['timepoint'])
    if key in motif_lookup:
        assigned = motif_lookup[key]
        if assigned in ccgg_motifs:
            hc_4mc_annotated.loc[idx, 'motif_group'] = 'CCGG_related'
        elif assigned == 'unassigned':
            hc_4mc_annotated.loc[idx, 'motif_group'] = 'unassigned'
        else:
            hc_4mc_annotated.loc[idx, 'motif_group'] = 'non_CCGG'

# n_reps distribution
print("\nReplicate consistency (n_reps) by group:")
for tp in ['T1', 'T2', 'T3']:
    print(f"\n  --- {tp} ---")
    for group in ['CCGG_related', 'non_CCGG', 'unassigned']:
        sub = hc_4mc_annotated[(hc_4mc_annotated['timepoint'] == tp) &
                                (hc_4mc_annotated['motif_group'] == group)]
        if len(sub) > 0:
            rep_counts = sub['n_reps'].value_counts().sort_index()
            mean_reps = sub['n_reps'].mean()
            pct_3reps = (sub['n_reps'] == 3).mean() * 100
            print(f"  4mC {group}: n={len(sub)}, mean_reps={mean_reps:.2f}, 3-rep%={pct_3reps:.1f}%")
            for r, n in rep_counts.items():
                print(f"    {r} reps: {n} ({n/len(sub)*100:.1f}%)")

    # 6mA for comparison
    sub_6ma = hc_6ma[hc_6ma['timepoint'] == tp]
    if len(sub_6ma) > 0:
        mean_reps = sub_6ma['n_reps'].mean()
        pct_3reps = (sub_6ma['n_reps'] == 3).mean() * 100
        print(f"  6mA: n={len(sub_6ma)}, mean_reps={mean_reps:.2f}, 3-rep%={pct_3reps:.1f}%")

# Statistical test: compare n_reps distribution
rep_comparison_stats = []
for tp in ['T1', 'T2', 'T3']:
    ccgg_reps = hc_4mc_annotated[(hc_4mc_annotated['timepoint'] == tp) &
                                  (hc_4mc_annotated['motif_group'] == 'CCGG_related')]['n_reps'].values
    non_ccgg_reps = hc_4mc_annotated[(hc_4mc_annotated['timepoint'] == tp) &
                                      (hc_4mc_annotated['motif_group'] == 'non_CCGG')]['n_reps'].values
    sixma_reps = hc_6ma[hc_6ma['timepoint'] == tp]['n_reps'].values

    if len(ccgg_reps) > 5 and len(non_ccgg_reps) > 5:
        u, p = stats.mannwhitneyu(ccgg_reps, non_ccgg_reps, alternative='two-sided')
        rep_comparison_stats.append({
            'timepoint': tp,
            'comparison': 'CCGG_4mC vs non_CCGG_4mC',
            'n1': len(ccgg_reps), 'n2': len(non_ccgg_reps),
            'mean1': ccgg_reps.mean(), 'mean2': non_ccgg_reps.mean(),
            'MannWhitney_U': u, 'MannWhitney_p': p,
        })

    if len(ccgg_reps) > 5 and len(sixma_reps) > 5:
        u, p = stats.mannwhitneyu(ccgg_reps, sixma_reps, alternative='two-sided')
        rep_comparison_stats.append({
            'timepoint': tp,
            'comparison': 'CCGG_4mC vs 6mA',
            'n1': len(ccgg_reps), 'n2': len(sixma_reps),
            'mean1': ccgg_reps.mean(), 'mean2': sixma_reps.mean(),
            'MannWhitney_U': u, 'MannWhitney_p': p,
        })

rep_stats_df = pd.DataFrame(rep_comparison_stats)
rep_stats_df.to_csv(os.path.join(TABDIR, 'replicate_consistency_statistics.tsv'), sep='\t', index=False)
print("\nReplicate consistency statistics:")
print(rep_stats_df.to_string(index=False))

# ============================================================
# 7. STEP 3 continued: weighted frequency comparison (HC level)
# ============================================================
print("\n" + "=" * 70)
print("[7] Weighted frequency comparison (high-confidence level)")
print("=" * 70)

# Compare weighted_mod_freq distributions
freq_comparison_stats = []
for tp in ['T1', 'T2', 'T3']:
    ccgg_freq = hc_4mc_annotated[(hc_4mc_annotated['timepoint'] == tp) &
                                  (hc_4mc_annotated['motif_group'] == 'CCGG_related')]['weighted_mod_freq'].values
    non_ccgg_freq = hc_4mc_annotated[(hc_4mc_annotated['timepoint'] == tp) &
                                      (hc_4mc_annotated['motif_group'] == 'non_CCGG')]['weighted_mod_freq'].values
    sixma_freq = hc_6ma[hc_6ma['timepoint'] == tp]['weighted_mod_freq'].values

    print(f"\n  --- {tp} ---")
    for name, vals in [('CCGG_4mC', ccgg_freq), ('non_CCGG_4mC', non_ccgg_freq), ('6mA', sixma_freq)]:
        if len(vals) > 0:
            print(f"  {name}: n={len(vals)}, mean={vals.mean():.1f}%, median={np.median(vals):.1f}%, "
                  f"std={vals.std():.1f}%, IQR=[{np.percentile(vals,25):.1f}, {np.percentile(vals,75):.1f}]")

    if len(ccgg_freq) > 5 and len(non_ccgg_freq) > 5:
        u, p = stats.mannwhitneyu(ccgg_freq, non_ccgg_freq, alternative='two-sided')
        ks, ksp = stats.ks_2samp(ccgg_freq, non_ccgg_freq)
        freq_comparison_stats.append({
            'timepoint': tp,
            'comparison': 'CCGG vs non_CCGG (weighted_freq)',
            'n1': len(ccgg_freq), 'n2': len(non_ccgg_freq),
            'mean1': ccgg_freq.mean(), 'mean2': non_ccgg_freq.mean(),
            'MannWhitney_U': u, 'MannWhitney_p': p,
            'KS_stat': ks, 'KS_p': ksp,
        })

    if len(ccgg_freq) > 5 and len(sixma_freq) > 5:
        u, p = stats.mannwhitneyu(ccgg_freq, sixma_freq, alternative='two-sided')
        ks, ksp = stats.ks_2samp(ccgg_freq, sixma_freq)
        freq_comparison_stats.append({
            'timepoint': tp,
            'comparison': 'CCGG_4mC vs 6mA (weighted_freq)',
            'n1': len(ccgg_freq), 'n2': len(sixma_freq),
            'mean1': ccgg_freq.mean(), 'mean2': sixma_freq.mean(),
            'MannWhitney_U': u, 'MannWhitney_p': p,
            'KS_stat': ks, 'KS_p': ksp,
        })

freq_stats_df = pd.DataFrame(freq_comparison_stats)
freq_stats_df.to_csv(os.path.join(TABDIR, 'weighted_freq_comparison_statistics.tsv'), sep='\t', index=False)
print("\nWeighted frequency comparison statistics:")
print(freq_stats_df.to_string(index=False))

# Also compare coverage
print("\n  Coverage comparison:")
for tp in ['T1', 'T2', 'T3']:
    ccgg_cov = hc_4mc_annotated[(hc_4mc_annotated['timepoint'] == tp) &
                                 (hc_4mc_annotated['motif_group'] == 'CCGG_related')]['total_coverage'].values
    non_ccgg_cov = hc_4mc_annotated[(hc_4mc_annotated['timepoint'] == tp) &
                                     (hc_4mc_annotated['motif_group'] == 'non_CCGG')]['total_coverage'].values
    sixma_cov = hc_6ma[hc_6ma['timepoint'] == tp]['total_coverage'].values

    print(f"  {tp}: CCGG mean_cov={ccgg_cov.mean():.0f} (n={len(ccgg_cov)}), "
          f"non-CCGG mean_cov={non_ccgg_cov.mean():.0f} (n={len(non_ccgg_cov)}), "
          f"6mA mean_cov={sixma_cov.mean():.0f} (n={len(sixma_cov)})")

# ============================================================
# 8. FIGURES
# ============================================================
print("\n" + "=" * 70)
print("[8] Generating figures")
print("=" * 70)

# ----- Figure 1: C position within CCGG -----
print("  Figure 1: C position within CCGG...")

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Panel A: Overall C position distribution
ax = axes[0]
pos_data = ccgg_pos_df['C_position_in_CCGG'].value_counts()
positions = ['C1', 'C2', 'ambiguous', 'not_in_CCGG']
counts = [pos_data.get(p, 0) for p in positions]
colors_pos = [COLORS['C1'], COLORS['C2'], '#999999', '#cccccc']
bars = ax.bar(range(len(positions)), counts, color=colors_pos, edgecolor='black', linewidth=0.5)
ax.set_xticks(range(len(positions)))
ax.set_xticklabels(positions, rotation=45, ha='right')
ax.set_ylabel('Number of sites')
ax.set_title('A. C position in CCGG\n(all CCGG-related 4mC sites)')
for bar, count in zip(bars, counts):
    if count > 0:
        ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 5,
                f'{count}\n({count/sum(counts)*100:.1f}%)', ha='center', va='bottom', fontsize=8)

# Panel B: C position by timepoint
ax = axes[1]
width = 0.25
x = np.arange(len(positions))
for i, tp in enumerate(['T1', 'T2', 'T3']):
    sub = ccgg_pos_df[ccgg_pos_df['timepoint'] == tp]
    tp_counts = [len(sub[sub['C_position_in_CCGG'] == p]) for p in positions]
    ax.bar(x + i*width, tp_counts, width, label=f'{tp} (n={len(sub)})',
           color=COLORS[tp], edgecolor='black', linewidth=0.5, alpha=0.8)
ax.set_xticks(x + width)
ax.set_xticklabels(positions, rotation=45, ha='right')
ax.set_ylabel('Number of sites')
ax.set_title('B. C position by timepoint')
ax.legend()

# Panel C: C position by motif
ax = axes[2]
motif_names = ['TGGCCGGC', 'CCGG', 'GGCCGG', 'CCGC']
width = 0.2
x = np.arange(len(positions))
motif_colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3']
for i, (m_name, mc) in enumerate(zip(motif_names, motif_colors)):
    sub = ccgg_pos_df[ccgg_pos_df['motif'] == m_name]
    m_counts = [len(sub[sub['C_position_in_CCGG'] == p]) for p in positions]
    ax.bar(x + i*width, m_counts, width, label=f'{m_name} (n={len(sub)})',
           color=mc, edgecolor='black', linewidth=0.5, alpha=0.8)
ax.set_xticks(x + 1.5*width)
ax.set_xticklabels(positions, rotation=45, ha='right')
ax.set_ylabel('Number of sites')
ax.set_title('C. C position by motif')
ax.legend(fontsize=7)

plt.tight_layout()
fig.savefig(os.path.join(FIGDIR, 'fig1_C_position_in_CCGG.pdf'), format='pdf')
fig.savefig(os.path.join(FIGDIR, 'fig1_C_position_in_CCGG.svg'), format='svg')
plt.close()
print("    Saved fig1_C_position_in_CCGG")

# ----- Figure 2: Score/frequency comparison (violin + box) -----
print("  Figure 2: Score/frequency comparison...")

fig, axes = plt.subplots(2, 3, figsize=(16, 10))

for col, tp in enumerate(['T1', 'T2', 'T3']):
    # Top row: Weighted modification frequency from HC sites
    ax = axes[0, col]
    plot_data = []
    plot_labels = []

    ccgg_f = hc_4mc_annotated[(hc_4mc_annotated['timepoint'] == tp) &
                               (hc_4mc_annotated['motif_group'] == 'CCGG_related')]['weighted_mod_freq'].values
    non_ccgg_f = hc_4mc_annotated[(hc_4mc_annotated['timepoint'] == tp) &
                                   (hc_4mc_annotated['motif_group'] == 'non_CCGG')]['weighted_mod_freq'].values
    sixma_f = hc_6ma[hc_6ma['timepoint'] == tp]['weighted_mod_freq'].values

    if len(ccgg_f) > 0:
        plot_data.append(ccgg_f)
        plot_labels.append(f'CCGG 4mC\n(n={len(ccgg_f)})')
    if len(non_ccgg_f) > 0:
        plot_data.append(non_ccgg_f)
        plot_labels.append(f'non-CCGG 4mC\n(n={len(non_ccgg_f)})')
    if len(sixma_f) > 0:
        plot_data.append(sixma_f)
        plot_labels.append(f'6mA\n(n={len(sixma_f)})')

    if len(plot_data) > 0:
        parts = ax.violinplot(plot_data, showmeans=True, showmedians=True)
        for i, pc in enumerate(parts['bodies']):
            color = [COLORS['CCGG_4mC'], COLORS['non_CCGG_4mC'], COLORS['6mA']][i] if i < 3 else 'gray'
            pc.set_facecolor(color)
            pc.set_alpha(0.6)
        ax.set_xticks(range(1, len(plot_labels) + 1))
        ax.set_xticklabels(plot_labels, fontsize=8)
    ax.set_ylabel('Weighted mod frequency (%)')
    ax.set_title(f'{tp} - HC weighted frequency')
    ax.set_ylim(45, 105)

    # Bottom row: Per-replicate frequencies from pileup
    ax = axes[1, col]
    tp_scores = score_df[(score_df['timepoint'] == tp) & (score_df['coverage'] >= min_cov)]

    plot_data2 = []
    plot_labels2 = []
    for g, color in [('CCGG_4mC', COLORS['CCGG_4mC']),
                      ('non_CCGG_4mC', COLORS['non_CCGG_4mC']),
                      ('6mA', COLORS['6mA'])]:
        vals = tp_scores[tp_scores['group'] == g]['mod_freq'].values
        if len(vals) > 0:
            plot_data2.append(vals)
            plot_labels2.append(f'{g}\n(n={len(vals)})')

    if len(plot_data2) > 0:
        parts = ax.violinplot(plot_data2, showmeans=True, showmedians=True)
        for i, pc in enumerate(parts['bodies']):
            color = [COLORS['CCGG_4mC'], COLORS['non_CCGG_4mC'], COLORS['6mA']][i] if i < 3 else 'gray'
            pc.set_facecolor(color)
            pc.set_alpha(0.6)
        ax.set_xticks(range(1, len(plot_labels2) + 1))
        ax.set_xticklabels(plot_labels2, fontsize=8)
    ax.set_ylabel('Per-replicate mod frequency (%)')
    ax.set_title(f'{tp} - Per-replicate frequency')

plt.suptitle('Detection Score Comparison: CCGG 4mC vs non-CCGG 4mC vs 6mA', fontsize=14, y=1.02)
plt.tight_layout()
fig.savefig(os.path.join(FIGDIR, 'fig2_score_frequency_comparison.pdf'), format='pdf')
fig.savefig(os.path.join(FIGDIR, 'fig2_score_frequency_comparison.svg'), format='svg')
plt.close()
print("    Saved fig2_score_frequency_comparison")

# ----- Figure 3: CCWGG context analysis -----
print("  Figure 3: Flanking context analysis...")

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Panel A: Base before CCGG
ax = axes[0]
if len(flanking_df) > 0:
    bases = ['A', 'C', 'G', 'T']
    base_colors = {'A': '#e41a1c', 'C': '#377eb8', 'G': '#4daf4a', 'T': '#984ea3'}

    # Expected from genome composition (random)
    # Count all CCGG flanking bases in genome
    genome_flanking_before = Counter()
    genome_flanking_after = Counter()
    for m in re.finditer('CCGG', genome_seq):
        s = m.start()
        if s > 0:
            genome_flanking_before[genome_seq[s-1]] += 1
        if s + 4 < len(genome_seq):
            genome_flanking_after[genome_seq[s+4]] += 1

    total_genome_before = sum(genome_flanking_before.values())
    total_genome_after = sum(genome_flanking_after.values())

    # Observed in methylated sites
    obs_before = flanking_df['base_before'].value_counts()
    total_obs = len(flanking_df)

    x = np.arange(len(bases))
    width = 0.35

    obs_vals = [obs_before.get(b, 0)/total_obs*100 for b in bases]
    exp_vals = [genome_flanking_before.get(b, 0)/total_genome_before*100 for b in bases]

    ax.bar(x - width/2, obs_vals, width, label='Methylated CCGG',
           color=[base_colors[b] for b in bases], edgecolor='black', linewidth=0.5, alpha=0.8)
    ax.bar(x + width/2, exp_vals, width, label='All CCGG (genome)',
           color=[base_colors[b] for b in bases], edgecolor='black', linewidth=0.5, alpha=0.3)

    ax.set_xticks(x)
    ax.set_xticklabels([f'{b}CCGG' for b in bases])
    ax.set_ylabel('Frequency (%)')
    ax.set_title('A. Base before CCGG')
    ax.legend()

# Panel B: Base after CCGG
ax = axes[1]
if len(flanking_df) > 0:
    obs_after = flanking_df['base_after'].value_counts()

    obs_vals = [obs_after.get(b, 0)/total_obs*100 for b in bases]
    exp_vals = [genome_flanking_after.get(b, 0)/total_genome_after*100 for b in bases]

    ax.bar(x - width/2, obs_vals, width, label='Methylated CCGG',
           color=[base_colors[b] for b in bases], edgecolor='black', linewidth=0.5, alpha=0.8)
    ax.bar(x + width/2, exp_vals, width, label='All CCGG (genome)',
           color=[base_colors[b] for b in bases], edgecolor='black', linewidth=0.5, alpha=0.3)

    ax.set_xticks(x)
    ax.set_xticklabels([f'CCGG{b}' for b in bases])
    ax.set_ylabel('Frequency (%)')
    ax.set_title('B. Base after CCGG')
    ax.legend()

# Panel C: Top 6-bp contexts
ax = axes[2]
if len(flanking_df) > 0:
    ctx_counts = flanking_df['six_bp_context'].value_counts().head(15)
    y_pos = range(len(ctx_counts))
    bars = ax.barh(y_pos, ctx_counts.values, color='steelblue', edgecolor='black', linewidth=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(ctx_counts.index, fontfamily='monospace', fontsize=8)
    ax.set_xlabel('Number of sites')
    ax.set_title('C. Top 6-bp contexts (nCCGGn)')
    ax.invert_yaxis()

    # Annotate with percentages
    for bar, count in zip(bars, ctx_counts.values):
        ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2,
                f'{count/total_obs*100:.1f}%', va='center', fontsize=7)

plt.tight_layout()
fig.savefig(os.path.join(FIGDIR, 'fig3_CCGG_flanking_context.pdf'), format='pdf')
fig.savefig(os.path.join(FIGDIR, 'fig3_CCGG_flanking_context.svg'), format='svg')
plt.close()
print("    Saved fig3_CCGG_flanking_context")

# ----- Figure 4: Replicate consistency comparison -----
print("  Figure 4: Replicate consistency comparison...")

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

for col, tp in enumerate(['T1', 'T2', 'T3']):
    ax = axes[col]

    groups = ['CCGG_related', 'non_CCGG', 'unassigned']
    group_labels = ['CCGG 4mC', 'non-CCGG 4mC', 'Unassigned 4mC']
    group_colors = [COLORS['CCGG_4mC'], COLORS['non_CCGG_4mC'], '#999999']

    x = np.arange(3)  # n_reps = 1, 2, 3 but HC sites have min 2
    width = 0.2

    for i, (g, gl, gc) in enumerate(zip(groups, group_labels, group_colors)):
        sub = hc_4mc_annotated[(hc_4mc_annotated['timepoint'] == tp) &
                                (hc_4mc_annotated['motif_group'] == g)]
        if len(sub) > 0:
            rep_vals = sub['n_reps'].value_counts().sort_index()
            vals = [rep_vals.get(r, 0)/len(sub)*100 for r in [1, 2, 3]]
            ax.bar(x + i*width, vals, width, label=f'{gl} (n={len(sub)})',
                   color=gc, edgecolor='black', linewidth=0.5, alpha=0.8)

    # Add 6mA
    sub_6ma = hc_6ma[hc_6ma['timepoint'] == tp]
    if len(sub_6ma) > 0:
        rep_vals = sub_6ma['n_reps'].value_counts().sort_index()
        vals = [rep_vals.get(r, 0)/len(sub_6ma)*100 for r in [1, 2, 3]]
        ax.bar(x + len(groups)*width, vals, width, label=f'6mA (n={len(sub_6ma)})',
               color=COLORS['6mA'], edgecolor='black', linewidth=0.5, alpha=0.8)

    ax.set_xticks(x + 1.5*width)
    ax.set_xticklabels(['1 rep', '2 reps', '3 reps'])
    ax.set_ylabel('Percentage of sites (%)')
    ax.set_title(f'{tp}')
    ax.legend(fontsize=7)

plt.suptitle('Replicate Consistency: CCGG 4mC vs non-CCGG 4mC vs 6mA', fontsize=14, y=1.02)
plt.tight_layout()
fig.savefig(os.path.join(FIGDIR, 'fig4_replicate_consistency.pdf'), format='pdf')
fig.savefig(os.path.join(FIGDIR, 'fig4_replicate_consistency.svg'), format='svg')
plt.close()
print("    Saved fig4_replicate_consistency")

# ----- Figure 5: Comprehensive summary figure -----
print("  Figure 5: Comprehensive summary...")

fig = plt.figure(figsize=(16, 12))
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.35)

# Panel A: C position pie chart
ax = fig.add_subplot(gs[0, 0])
pos_data_summary = ccgg_pos_df['C_position_in_CCGG'].value_counts()
labels = []
sizes = []
colors_pie = []
for p in ['C2', 'C1', 'ambiguous', 'not_in_CCGG']:
    n = pos_data_summary.get(p, 0)
    if n > 0:
        labels.append(f'{p}\n({n}, {n/len(ccgg_pos_df)*100:.1f}%)')
        sizes.append(n)
        colors_pie.append({'C1': COLORS['C1'], 'C2': COLORS['C2'],
                          'ambiguous': '#999999', 'not_in_CCGG': '#cccccc'}[p])

wedges, texts = ax.pie(sizes, labels=labels, colors=colors_pie, startangle=90,
                        textprops={'fontsize': 8})
ax.set_title('A. Methylated C position\nwithin CCGG motif', fontsize=11)

# Panel B: Frequency comparison boxplot (T1 only, largest sample)
ax = fig.add_subplot(gs[0, 1])
t1_ccgg = hc_4mc_annotated[(hc_4mc_annotated['timepoint'] == 'T1') &
                             (hc_4mc_annotated['motif_group'] == 'CCGG_related')]['weighted_mod_freq'].values
t1_non = hc_4mc_annotated[(hc_4mc_annotated['timepoint'] == 'T1') &
                           (hc_4mc_annotated['motif_group'] == 'non_CCGG')]['weighted_mod_freq'].values
t1_6ma = hc_6ma[hc_6ma['timepoint'] == 'T1']['weighted_mod_freq'].values

box_data = [t1_ccgg, t1_non, t1_6ma]
box_labels = [f'CCGG 4mC\n(n={len(t1_ccgg)})', f'non-CCGG 4mC\n(n={len(t1_non)})', f'6mA\n(n={len(t1_6ma)})']
box_colors = [COLORS['CCGG_4mC'], COLORS['non_CCGG_4mC'], COLORS['6mA']]

bp = ax.boxplot(box_data, labels=box_labels, patch_artist=True, widths=0.6)
for patch, color in zip(bp['boxes'], box_colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.6)
ax.set_ylabel('Weighted mod frequency (%)')
ax.set_title('B. Modification frequency\n(T1, high-confidence)', fontsize=11)

# Panel C: n_reps comparison (T1)
ax = fig.add_subplot(gs[0, 2])
for g, gl, gc in [('CCGG_related', 'CCGG 4mC', COLORS['CCGG_4mC']),
                    ('non_CCGG', 'non-CCGG 4mC', COLORS['non_CCGG_4mC'])]:
    sub = hc_4mc_annotated[(hc_4mc_annotated['timepoint'] == 'T1') &
                            (hc_4mc_annotated['motif_group'] == g)]
    if len(sub) > 0:
        pct_3rep = (sub['n_reps'] == 3).mean() * 100
        pct_2rep = (sub['n_reps'] == 2).mean() * 100
        ax.bar([gl], [pct_3rep], color=gc, alpha=0.8, label='3 reps', edgecolor='black', linewidth=0.5)
        ax.bar([gl], [pct_2rep], bottom=[pct_3rep], color=gc, alpha=0.4, label='2 reps', edgecolor='black', linewidth=0.5)

sub_6ma = hc_6ma[hc_6ma['timepoint'] == 'T1']
if len(sub_6ma) > 0:
    pct_3rep = (sub_6ma['n_reps'] == 3).mean() * 100
    pct_2rep = (sub_6ma['n_reps'] == 2).mean() * 100
    ax.bar(['6mA'], [pct_3rep], color=COLORS['6mA'], alpha=0.8, edgecolor='black', linewidth=0.5)
    ax.bar(['6mA'], [pct_2rep], bottom=[pct_3rep], color=COLORS['6mA'], alpha=0.4, edgecolor='black', linewidth=0.5)

ax.set_ylabel('Percentage of sites (%)')
ax.set_title('C. Replicate consistency\n(T1)', fontsize=11)
ax.set_ylim(0, 105)

# Panel D: 6-bp context analysis
ax = fig.add_subplot(gs[1, 0])
if len(flanking_df) > 0:
    ctx_counts = flanking_df['six_bp_context'].value_counts().head(10)
    y_pos = range(len(ctx_counts))
    ax.barh(y_pos, ctx_counts.values, color='steelblue', edgecolor='black', linewidth=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(ctx_counts.index, fontfamily='monospace', fontsize=8)
    ax.set_xlabel('Number of sites')
    ax.set_title('D. Top 6-bp contexts', fontsize=11)
    ax.invert_yaxis()

# Panel E: CDF comparison
ax = fig.add_subplot(gs[1, 1])
for name, vals, color in [('CCGG 4mC', t1_ccgg, COLORS['CCGG_4mC']),
                           ('non-CCGG 4mC', t1_non, COLORS['non_CCGG_4mC']),
                           ('6mA', t1_6ma, COLORS['6mA'])]:
    if len(vals) > 0:
        sorted_v = np.sort(vals)
        cdf = np.arange(1, len(sorted_v) + 1) / len(sorted_v)
        ax.plot(sorted_v, cdf, label=f'{name} (n={len(vals)})', color=color, linewidth=2)

ax.set_xlabel('Weighted mod frequency (%)')
ax.set_ylabel('Cumulative proportion')
ax.set_title('E. Frequency CDF\n(T1)', fontsize=11)
ax.legend(fontsize=8)
ax.set_xlim(45, 105)

# Panel F: Evidence summary table
ax = fig.add_subplot(gs[1, 2])
ax.axis('off')

# Create a summary text
c2_pct = pos_data_summary.get('C2', 0) / len(ccgg_pos_df) * 100 if len(ccgg_pos_df) > 0 else 0
c1_pct = pos_data_summary.get('C1', 0) / len(ccgg_pos_df) * 100 if len(ccgg_pos_df) > 0 else 0

evidence_text = "H9 Evidence Summary\n"
evidence_text += "=" * 35 + "\n\n"
evidence_text += f"C2 position (Dcm target): {c2_pct:.1f}%\n"
evidence_text += f"C1 position: {c1_pct:.1f}%\n\n"

if len(t1_ccgg) > 0 and len(t1_non) > 0:
    evidence_text += f"T1 CCGG freq: {t1_ccgg.mean():.1f}% (mean)\n"
    evidence_text += f"T1 non-CCGG freq: {t1_non.mean():.1f}% (mean)\n\n"

ccgg_3rep = hc_4mc_annotated[(hc_4mc_annotated['timepoint'] == 'T1') &
                              (hc_4mc_annotated['motif_group'] == 'CCGG_related')]
if len(ccgg_3rep) > 0:
    evidence_text += f"CCGG 3-rep: {(ccgg_3rep['n_reps']==3).mean()*100:.1f}%\n"
non_ccgg_3rep = hc_4mc_annotated[(hc_4mc_annotated['timepoint'] == 'T1') &
                                  (hc_4mc_annotated['motif_group'] == 'non_CCGG')]
if len(non_ccgg_3rep) > 0:
    evidence_text += f"non-CCGG 3-rep: {(non_ccgg_3rep['n_reps']==3).mean()*100:.1f}%\n"

ax.text(0.05, 0.95, evidence_text, transform=ax.transAxes, fontsize=9,
        verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='lightyellow', alpha=0.8))
ax.set_title('F. Evidence summary', fontsize=11)

plt.suptitle('H9: Is CCGG "4mC" actually 5mC misclassified by Nanopore?',
             fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()
fig.savefig(os.path.join(FIGDIR, 'fig5_comprehensive_summary.pdf'), format='pdf')
fig.savefig(os.path.join(FIGDIR, 'fig5_comprehensive_summary.svg'), format='svg')
plt.close()
print("    Saved fig5_comprehensive_summary")

# ============================================================
# 9. Extended analysis: TGGCCGGC motif deep dive
# ============================================================
print("\n" + "=" * 70)
print("[9] TGGCCGGC motif deep dive")
print("=" * 70)

# TGGCCGGC is the dominant CCGG-containing motif (n=1717)
# Within TGGCCGGC, where exactly is the methylated C?
# TGGCCGGC: T(0) G(1) G(2) C(3) C(4) G(5) G(6) C(7)
# Contains CCGG at positions 3-6
# The C at position 4 = C2 of CCGG
# The C at position 3 = C1 of CCGG
# Also C at position 7 is not part of CCGG

tggccggc_sites = motif[motif['assigned_motif'] == 'TGGCCGGC'].copy()
print(f"  TGGCCGGC sites: {len(tggccggc_sites)}")

# For each TGGCCGGC site, the center_base should tell us which base is modified
print(f"  Center base distribution:")
print(tggccggc_sites['center_base'].value_counts().to_string())

# Analyze the sequence context more carefully
# The 'sequence' column has 31bp with the methylated C at center (position 15)
# We need to find where the TGGCCGGC motif falls relative to the methylated position

def find_motif_in_sequence(seq, motif_str, center_pos=15):
    """
    Find motif_str in sequence and return the position of the methylated base
    relative to the motif.
    center_pos: position of the methylated base in the sequence (0-based)
    """
    # Search for motif in forward direction
    for m in re.finditer(motif_str, seq):
        rel_pos = center_pos - m.start()
        return rel_pos, m.start(), '+'

    # Search for reverse complement
    rc_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    rc_motif = ''.join(rc_map.get(b, b) for b in reversed(motif_str))
    for m in re.finditer(rc_motif, seq):
        rel_pos = center_pos - m.start()
        # In RC, position is reversed
        motif_len = len(motif_str)
        rc_rel_pos = motif_len - 1 - rel_pos
        return rc_rel_pos, m.start(), '-'

    return None, None, None

tggccggc_pos_in_motif = []
for _, row in tggccggc_sites.iterrows():
    seq = row['sequence']
    rel_pos, motif_start, orientation = find_motif_in_sequence(seq, 'TGGCCGGC')

    if rel_pos is not None:
        # Map to TGGCCGGC positions
        # TGGCCGGC: T=0, G=1, G=2, C=3, C=4, G=5, G=6, C=7
        # The CCGG within is at positions 3,4,5,6
        # C1 of CCGG = position 3
        # C2 of CCGG = position 4

        tggccggc_pos_in_motif.append({
            'position': row['position'],
            'strand': row['strand'],
            'timepoint': row['timepoint'],
            'frequency': row['frequency'],
            'rel_pos_in_motif': rel_pos,
            'motif_base': 'TGGCCGGC'[rel_pos] if 0 <= rel_pos < 8 else '?',
            'in_CCGG_C1': rel_pos == 3,
            'in_CCGG_C2': rel_pos == 4,
            'in_CCGG_other_C': rel_pos == 7,
            'orientation': orientation,
        })
    else:
        tggccggc_pos_in_motif.append({
            'position': row['position'],
            'strand': row['strand'],
            'timepoint': row['timepoint'],
            'frequency': row['frequency'],
            'rel_pos_in_motif': -1,
            'motif_base': '?',
            'in_CCGG_C1': False,
            'in_CCGG_C2': False,
            'in_CCGG_other_C': False,
            'orientation': '?',
        })

tggccggc_pos_df = pd.DataFrame(tggccggc_pos_in_motif)

print(f"\nTGGCCGGC methylated position distribution:")
for pos_val in sorted(tggccggc_pos_df['rel_pos_in_motif'].unique()):
    n = len(tggccggc_pos_df[tggccggc_pos_df['rel_pos_in_motif'] == pos_val])
    base = 'TGGCCGGC'[pos_val] if 0 <= pos_val < 8 else '?'
    ccgg_label = ''
    if pos_val == 3:
        ccgg_label = ' (C1 of CCGG)'
    elif pos_val == 4:
        ccgg_label = ' (C2 of CCGG = Dcm target)'
    elif pos_val == 7:
        ccgg_label = ' (terminal C)'
    print(f"  Position {pos_val} ({base}){ccgg_label}: {n} ({n/len(tggccggc_pos_df)*100:.1f}%)")

print(f"\n  Summary:")
n_c1 = tggccggc_pos_df['in_CCGG_C1'].sum()
n_c2 = tggccggc_pos_df['in_CCGG_C2'].sum()
n_c7 = tggccggc_pos_df['in_CCGG_other_C'].sum()
print(f"    CCGG C1 (pos 3): {n_c1} ({n_c1/len(tggccggc_pos_df)*100:.1f}%)")
print(f"    CCGG C2/Dcm (pos 4): {n_c2} ({n_c2/len(tggccggc_pos_df)*100:.1f}%)")
print(f"    Terminal C (pos 7): {n_c7} ({n_c7/len(tggccggc_pos_df)*100:.1f}%)")

tggccggc_pos_df.to_csv(os.path.join(TABDIR, 'TGGCCGGC_position_analysis.tsv'), sep='\t', index=False)

# ============================================================
# 10. Comprehensive summary statistics
# ============================================================
print("\n" + "=" * 70)
print("[10] COMPREHENSIVE SUMMARY")
print("=" * 70)

# Compile all evidence
summary_records = []

# Evidence 1: C position
c2_total = len(ccgg_pos_df[ccgg_pos_df['C_position_in_CCGG'] == 'C2'])
c1_total = len(ccgg_pos_df[ccgg_pos_df['C_position_in_CCGG'] == 'C1'])
c2_pct = c2_total / len(ccgg_pos_df) * 100 if len(ccgg_pos_df) > 0 else 0
c1_pct = c1_total / len(ccgg_pos_df) * 100 if len(ccgg_pos_df) > 0 else 0

summary_records.append({
    'evidence': 'C position in CCGG',
    'observation': f'C2={c2_total} ({c2_pct:.1f}%), C1={c1_total} ({c1_pct:.1f}%)',
    'supports_H9': 'YES' if c2_pct > c1_pct else 'PARTIAL' if c2_pct > 0 else 'NO',
    'note': 'C2 is the Dcm target position',
})

# Evidence 2: Frequency comparison
if len(freq_stats_df) > 0:
    t1_row = freq_stats_df[(freq_stats_df['timepoint'] == 'T1') &
                            (freq_stats_df['comparison'].str.contains('non_CCGG'))]
    if len(t1_row) > 0:
        r = t1_row.iloc[0]
        lower = r['mean1'] < r['mean2']
        summary_records.append({
            'evidence': 'Mod frequency (CCGG vs non-CCGG 4mC)',
            'observation': f'CCGG={r["mean1"]:.1f}% vs non-CCGG={r["mean2"]:.1f}%, p={r["MannWhitney_p"]:.2e}',
            'supports_H9': 'YES' if lower and r['MannWhitney_p'] < 0.05 else 'PARTIAL' if lower else 'NO',
            'note': 'Lower frequency suggests weaker/misclassified signal',
        })

# Evidence 3: Replicate consistency
if len(rep_stats_df) > 0:
    t1_rep = rep_stats_df[(rep_stats_df['timepoint'] == 'T1') &
                           (rep_stats_df['comparison'].str.contains('non_CCGG'))]
    if len(t1_rep) > 0:
        r = t1_rep.iloc[0]
        lower = r['mean1'] < r['mean2']
        summary_records.append({
            'evidence': 'Replicate consistency (n_reps)',
            'observation': f'CCGG={r["mean1"]:.2f} vs non-CCGG={r["mean2"]:.2f} mean reps, p={r["MannWhitney_p"]:.2e}',
            'supports_H9': 'YES' if lower and r['MannWhitney_p'] < 0.05 else 'PARTIAL' if lower else 'NO',
            'note': 'Lower consistency suggests stochastic/artifactual detection',
        })

# Evidence 4: Temporal instability (zero overlap from H7)
summary_records.append({
    'evidence': 'Temporal site overlap',
    'observation': 'Jaccard=0.000 between T1 and T2 (from H7)',
    'supports_H9': 'YES',
    'note': 'Complete lack of positional overlap is inconsistent with enzymatic methylation',
})

# Evidence 5: No adequate MTase expression
summary_records.append({
    'evidence': 'MTase expression',
    'observation': 'No identified cytosine MTase with adequate T1 expression (from H7)',
    'supports_H9': 'YES',
    'note': 'If sites are artifacts, no MTase needed',
})

# TGGCCGGC specific evidence
if len(tggccggc_pos_df) > 0:
    n_c2_tggc = tggccggc_pos_df['in_CCGG_C2'].sum()
    pct_c2_tggc = n_c2_tggc / len(tggccggc_pos_df) * 100
    summary_records.append({
        'evidence': 'TGGCCGGC C position',
        'observation': f'C2 of CCGG (Dcm target): {n_c2_tggc} ({pct_c2_tggc:.1f}%)',
        'supports_H9': 'YES' if pct_c2_tggc > 30 else 'PARTIAL' if pct_c2_tggc > 10 else 'NO',
        'note': 'High C2 fraction supports 5mC at Dcm target misread as 4mC',
    })

summary_df = pd.DataFrame(summary_records)
summary_df.to_csv(os.path.join(TABDIR, 'H9_evidence_summary.tsv'), sep='\t', index=False)

print("\nH9 Evidence Summary:")
print(summary_df.to_string(index=False))

# Count support
n_yes = (summary_df['supports_H9'] == 'YES').sum()
n_partial = (summary_df['supports_H9'] == 'PARTIAL').sum()
n_no = (summary_df['supports_H9'] == 'NO').sum()
total = len(summary_df)

print(f"\n  YES: {n_yes}/{total}, PARTIAL: {n_partial}/{total}, NO: {n_no}/{total}")

verdict = ""
if n_yes >= total * 0.5:
    verdict = "SUPPORTED - Majority of evidence supports CCGG 4mC as misclassified 5mC"
elif n_yes + n_partial >= total * 0.5:
    verdict = "PARTIALLY SUPPORTED - Mixed evidence, some consistent with misclassification"
else:
    verdict = "NOT SUPPORTED - Evidence does not support systematic 5mC misclassification"

print(f"\n  VERDICT: {verdict}")

# ============================================================
# DONE
# ============================================================
print("\n" + "=" * 70)
print("Analysis complete!")
print(f"  Figures: {FIGDIR}/")
print(f"  Tables: {TABDIR}/")
print("=" * 70)
