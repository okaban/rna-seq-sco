#!/usr/bin/env python3
"""
H30: TSS-Proximal DNA Sequence Features as Determinants of Protection Zone

Tests whether exposed regulators (no protection zone) have higher density of
R-M recognition motifs in the TSS-proximal region compared to shielded regulators.

Background:
- H27: 57 "exposed" regulators have no protection zone (methylation at median 114bp from TSS)
         998 "shielded" regulators have a 1,200bp protection zone
- H29: Expression level does NOT predict this classification (baseMean AUC=0.547)
- H22: R-M motifs depleted in regulatory gene bodies (fold=0.74-0.76)

Hypothesis: Exposed regulators have higher density of R-M recognition motifs
in the TSS +/-300bp region, and this is a structural determinant of protection.
"""

import re
import numpy as np
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from scipy import stats
from collections import defaultdict, Counter
from itertools import product
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc, roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Configuration
# ============================================================
BASE_DIR = Path('/Users/okaban/bioinfo/rna-seq')
ANALYSIS_DIR = BASE_DIR / '11_epigenome_integration/analysis/53_TSS_sequence_determinants'
FIGURES_DIR = ANALYSIS_DIR / 'figures'
TABLES_DIR = ANALYSIS_DIR / 'tables'

GENOME_FASTA = Path('/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/GCF_000203835.1_ASM20383v1_genomic.fna')
GENOME_FASTA_ALT = Path('/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa')

GENE_ANNOTATION = BASE_DIR / '05_annotation/analysis/05_annotation_260128_v1/tables/gene_annotation_basic.tsv'
ALL_REGULATORY = BASE_DIR / '11_epigenome_integration/analysis/29_genomewide_TF_screen/tables/all_regulatory_genes.tsv'
COORDINATED_REGULATORY = BASE_DIR / '11_epigenome_integration/analysis/29_genomewide_TF_screen/tables/coordinated_regulatory_genes.tsv'
H29_FEATURES = BASE_DIR / '11_epigenome_integration/analysis/52_shielded_exposed_boundary/tables/all_genes_features.tsv'
H27_METRICS = BASE_DIR / '11_epigenome_integration/analysis/50_coordinated_regulators_protection/tables/gene_level_metrics.tsv'
H22_MOTIF_COUNTS = BASE_DIR / '11_epigenome_integration/analysis/45_sequence_level_motif_depletion/tables/gene_motif_counts.tsv'

# Genomic parameters
CHROMOSOME_LENGTH = 8_667_507
ARM_BOUNDARY = 1_500_000

# Window sizes to analyze
WINDOWS = [300, 500, 1000]

# R-M motifs to search
MOTIFS = {
    'TGGCCGGC': 'GCCGGCCA',   # 8-mer with RC
    'AAGCCCG': 'CGGGCTT',      # 7-mer with RC
    'GCCGGC': 'GCCGGC',        # Core 6-mer (palindrome)
    'CCGG': 'CCGG',             # Common methylation target (palindrome)
}

print("=" * 70)
print("H30: TSS-Proximal DNA Sequence Features as Determinants of Protection Zone")
print("=" * 70)

# ============================================================
# 1. Load genome sequence
# ============================================================
print("\n[1] Loading genome sequence...")
genome_path = GENOME_FASTA if GENOME_FASTA.exists() else GENOME_FASTA_ALT
genome = {}
for record in SeqIO.parse(genome_path, 'fasta'):
    genome[record.id] = str(record.seq).upper()
    print(f"  Loaded: {record.id} ({len(record.seq):,} bp)")

# Get the main chromosome
chrom_id = [k for k in genome.keys() if 'NC_003888' in k][0]
chrom_seq = genome[chrom_id]
print(f"  Chromosome: {chrom_id}, length={len(chrom_seq):,} bp")

# ============================================================
# 2. Load gene data and classify exposed/shielded
# ============================================================
print("\n[2] Loading regulatory gene data...")

# H29 feature matrix has exposed/shielded classification
h29 = pd.read_csv(H29_FEATURES, sep='\t')
print(f"  H29 feature matrix: {len(h29)} genes")
print(f"  Exposed: {h29['is_exposed'].sum()}, Shielded: {(~h29['is_exposed'].astype(bool)).sum()}")

# All regulatory genes
all_reg = pd.read_csv(ALL_REGULATORY, sep='\t')
print(f"  All regulatory genes: {len(all_reg)}")

# Coordinated regulatory genes (exposed)
coord_reg = pd.read_csv(COORDINATED_REGULATORY, sep='\t')
print(f"  Coordinated (exposed): {len(coord_reg)}")

# H27 metrics
h27 = pd.read_csv(H27_METRICS, sep='\t')
print(f"  H27 gene-level metrics: {len(h27)}")

# Use H29 feature matrix as the primary source
# is_exposed: 1=exposed, 0=shielded
df = h29.copy()
df['is_exposed'] = df['is_exposed'].astype(int)
df['class'] = df['is_exposed'].map({1: 'Exposed', 0: 'Shielded'})

n_exposed = (df['is_exposed'] == 1).sum()
n_shielded = (df['is_exposed'] == 0).sum()
print(f"\n  Classification: {n_exposed} exposed, {n_shielded} shielded")

# ============================================================
# 3. Extract TSS-proximal sequences
# ============================================================
print("\n[3] Extracting TSS-proximal sequences...")

def reverse_complement(seq):
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq))

def extract_tss_region(tss, strand, chrom_seq, window):
    """Extract sequence around TSS. For + strand genes, TSS=start; for - strand, TSS=end."""
    chrom_len = len(chrom_seq)
    left = max(0, tss - window)
    right = min(chrom_len, tss + window)
    seq = chrom_seq[left:right]
    if strand == '-':
        seq = reverse_complement(seq)
    return seq

# Extract sequences for each window size
for w in WINDOWS:
    col = f'seq_{w}bp'
    seqs = []
    for _, row in df.iterrows():
        tss = int(row['tss'])
        strand = row['strand']
        seq = extract_tss_region(tss, strand, chrom_seq, w)
        seqs.append(seq)
    df[col] = seqs
    print(f"  Extracted {w}bp window: median length = {np.median([len(s) for s in seqs]):.0f}")

# ============================================================
# 4. Count R-M recognition motifs in TSS region
# ============================================================
print("\n[4] Counting R-M recognition motifs...")

def count_motif(seq, motif, rc_motif):
    """Count occurrences of motif and its reverse complement in sequence."""
    count = 0
    for m in [motif, rc_motif]:
        start = 0
        while True:
            pos = seq.find(m, start)
            if pos == -1:
                break
            count += 1
            start = pos + 1
    # If motif is a palindrome, we counted it once already per occurrence
    if motif == rc_motif:
        # Actually for palindromes, find(motif) already catches both strands
        pass
    return count

for motif_name, rc in MOTIFS.items():
    for w in WINDOWS:
        col_seq = f'seq_{w}bp'
        col_count = f'{motif_name}_count_{w}bp'
        col_density = f'{motif_name}_density_{w}bp'
        counts = []
        for seq in df[col_seq]:
            c = count_motif(seq, motif_name, rc)
            counts.append(c)
        df[col_count] = counts
        df[col_density] = [c / (len(s) / 1000) if len(s) > 0 else 0
                            for c, s in zip(counts, df[col_seq])]

    # Summary
    exp_mean = df.loc[df['is_exposed'] == 1, f'{motif_name}_count_300bp'].mean()
    shi_mean = df.loc[df['is_exposed'] == 0, f'{motif_name}_count_300bp'].mean()
    print(f"  {motif_name} ±300bp: exposed mean={exp_mean:.3f}, shielded mean={shi_mean:.3f}")

# ============================================================
# 5. Statistical comparisons: exposed vs shielded
# ============================================================
print("\n[5] Statistical comparisons (exposed vs shielded)...")

stat_results = []

for motif_name in MOTIFS.keys():
    for w in WINDOWS:
        col = f'{motif_name}_count_{w}bp'

        exp_vals = df.loc[df['is_exposed'] == 1, col].values
        shi_vals = df.loc[df['is_exposed'] == 0, col].values

        # Wilcoxon rank-sum test
        if np.std(exp_vals) > 0 or np.std(shi_vals) > 0:
            u_stat, p_wilcox = stats.mannwhitneyu(exp_vals, shi_vals, alternative='two-sided')
            # Rank-biserial correlation
            n1, n2 = len(exp_vals), len(shi_vals)
            r_rb = 1 - (2 * u_stat) / (n1 * n2)
        else:
            p_wilcox = 1.0
            r_rb = 0.0

        # Fisher exact: has_motif vs no_motif
        exp_has = np.sum(exp_vals > 0)
        exp_no = np.sum(exp_vals == 0)
        shi_has = np.sum(shi_vals > 0)
        shi_no = np.sum(shi_vals == 0)

        if (exp_has + shi_has) > 0:
            table = [[exp_has, exp_no], [shi_has, shi_no]]
            or_fisher, p_fisher = stats.fisher_exact(table)
        else:
            or_fisher, p_fisher = np.nan, 1.0

        stat_results.append({
            'motif': motif_name,
            'window': f'±{w}bp',
            'exposed_mean': np.mean(exp_vals),
            'exposed_median': np.median(exp_vals),
            'shielded_mean': np.mean(shi_vals),
            'shielded_median': np.median(shi_vals),
            'fold_change': np.mean(exp_vals) / np.mean(shi_vals) if np.mean(shi_vals) > 0 else np.inf,
            'wilcox_U': u_stat if np.std(exp_vals) > 0 or np.std(shi_vals) > 0 else np.nan,
            'wilcox_p': p_wilcox,
            'rank_biserial_r': r_rb,
            'fisher_OR': or_fisher,
            'fisher_p': p_fisher,
            'exposed_has_motif_pct': exp_has / len(exp_vals) * 100 if len(exp_vals) > 0 else 0,
            'shielded_has_motif_pct': shi_has / len(shi_vals) * 100 if len(shi_vals) > 0 else 0,
        })

        if w == 300:
            sig = "***" if p_wilcox < 0.001 else "**" if p_wilcox < 0.01 else "*" if p_wilcox < 0.05 else "NS"
            print(f"  {motif_name} ±{w}bp: exp={np.mean(exp_vals):.3f} vs shi={np.mean(shi_vals):.3f}, "
                  f"fold={np.mean(exp_vals)/np.mean(shi_vals):.2f}, p={p_wilcox:.3e} {sig}")

stat_df = pd.DataFrame(stat_results)
stat_df.to_csv(TABLES_DIR / 'statistical_tests.tsv', sep='\t', index=False)
print(f"\n  Saved: {TABLES_DIR / 'statistical_tests.tsv'}")

# ============================================================
# 6. GC content and dinucleotide composition
# ============================================================
print("\n[6] GC content and dinucleotide composition...")

def calc_gc(seq):
    if len(seq) == 0:
        return 0
    gc = seq.count('G') + seq.count('C')
    return gc / len(seq) * 100

def calc_cpg_oe(seq):
    """CpG observed/expected ratio."""
    if len(seq) < 2:
        return 0
    cpg = seq.count('CG')
    c = seq.count('C')
    g = seq.count('G')
    expected = (c * g) / len(seq) if len(seq) > 0 else 0
    return cpg / expected if expected > 0 else 0

def calc_dinucleotide_freq(seq):
    """Calculate all 16 dinucleotide frequencies."""
    if len(seq) < 2:
        return {}
    total = len(seq) - 1
    freqs = {}
    for d1, d2 in product('ACGT', repeat=2):
        di = d1 + d2
        count = sum(1 for i in range(total) if seq[i:i+2] == di)
        freqs[di] = count / total
    return freqs

# GC content for each window
for w in WINDOWS:
    col_seq = f'seq_{w}bp'
    df[f'GC_{w}bp'] = df[col_seq].apply(calc_gc)
    df[f'CpG_OE_{w}bp'] = df[col_seq].apply(calc_cpg_oe)

# Compare GC content
for w in WINDOWS:
    exp_gc = df.loc[df['is_exposed'] == 1, f'GC_{w}bp'].values
    shi_gc = df.loc[df['is_exposed'] == 0, f'GC_{w}bp'].values
    u, p = stats.mannwhitneyu(exp_gc, shi_gc, alternative='two-sided')
    n1, n2 = len(exp_gc), len(shi_gc)
    r_rb = 1 - (2 * u) / (n1 * n2)

    stat_results.append({
        'motif': 'GC_content',
        'window': f'±{w}bp',
        'exposed_mean': np.mean(exp_gc),
        'exposed_median': np.median(exp_gc),
        'shielded_mean': np.mean(shi_gc),
        'shielded_median': np.median(shi_gc),
        'fold_change': np.mean(exp_gc) / np.mean(shi_gc) if np.mean(shi_gc) > 0 else np.inf,
        'wilcox_U': u,
        'wilcox_p': p,
        'rank_biserial_r': r_rb,
        'fisher_OR': np.nan,
        'fisher_p': np.nan,
        'exposed_has_motif_pct': np.nan,
        'shielded_has_motif_pct': np.nan,
    })

    sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "NS"
    print(f"  GC ±{w}bp: exp={np.mean(exp_gc):.2f}% vs shi={np.mean(shi_gc):.2f}%, p={p:.3e} {sig}")

# CpG O/E
for w in WINDOWS:
    exp_cpg = df.loc[df['is_exposed'] == 1, f'CpG_OE_{w}bp'].values
    shi_cpg = df.loc[df['is_exposed'] == 0, f'CpG_OE_{w}bp'].values
    u, p = stats.mannwhitneyu(exp_cpg, shi_cpg, alternative='two-sided')
    n1, n2 = len(exp_cpg), len(shi_cpg)
    r_rb = 1 - (2 * u) / (n1 * n2)

    stat_results.append({
        'motif': 'CpG_OE',
        'window': f'±{w}bp',
        'exposed_mean': np.mean(exp_cpg),
        'exposed_median': np.median(exp_cpg),
        'shielded_mean': np.mean(shi_cpg),
        'shielded_median': np.median(shi_cpg),
        'fold_change': np.mean(exp_cpg) / np.mean(shi_cpg) if np.mean(shi_cpg) > 0 else np.inf,
        'wilcox_U': u,
        'wilcox_p': p,
        'rank_biserial_r': r_rb,
        'fisher_OR': np.nan,
        'fisher_p': np.nan,
        'exposed_has_motif_pct': np.nan,
        'shielded_has_motif_pct': np.nan,
    })

    sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "NS"
    print(f"  CpG O/E ±{w}bp: exp={np.mean(exp_cpg):.4f} vs shi={np.mean(shi_cpg):.4f}, p={p:.3e} {sig}")

# Dinucleotide frequencies for ±300bp
print("\n  Calculating dinucleotide frequencies (±300bp)...")
dinuc_data = []
for idx, row in df.iterrows():
    seq = row['seq_300bp']
    freqs = calc_dinucleotide_freq(seq)
    freqs['locus_tag'] = row['locus_tag']
    freqs['is_exposed'] = row['is_exposed']
    freqs['class'] = row['class']
    dinuc_data.append(freqs)

dinuc_df = pd.DataFrame(dinuc_data)
dinuc_df.to_csv(TABLES_DIR / 'dinucleotide_frequencies.tsv', sep='\t', index=False)

# Compare each dinucleotide
dinuc_results = []
dinucleotides = [d1+d2 for d1, d2 in product('ACGT', repeat=2)]
for di in dinucleotides:
    exp_vals = dinuc_df.loc[dinuc_df['is_exposed'] == 1, di].values
    shi_vals = dinuc_df.loc[dinuc_df['is_exposed'] == 0, di].values
    u, p = stats.mannwhitneyu(exp_vals, shi_vals, alternative='two-sided')
    n1, n2 = len(exp_vals), len(shi_vals)
    r_rb = 1 - (2 * u) / (n1 * n2)
    dinuc_results.append({
        'dinucleotide': di,
        'exposed_mean': np.mean(exp_vals),
        'shielded_mean': np.mean(shi_vals),
        'diff': np.mean(exp_vals) - np.mean(shi_vals),
        'wilcox_p': p,
        'rank_biserial_r': r_rb,
    })

dinuc_stat = pd.DataFrame(dinuc_results).sort_values('wilcox_p')
print(f"  Top 5 differential dinucleotides (±300bp):")
for _, row in dinuc_stat.head(5).iterrows():
    sig = "***" if row['wilcox_p'] < 0.001 else "**" if row['wilcox_p'] < 0.01 else "*" if row['wilcox_p'] < 0.05 else "NS"
    print(f"    {row['dinucleotide']}: exp={row['exposed_mean']:.5f} vs shi={row['shielded_mean']:.5f}, "
          f"diff={row['diff']:.5f}, p={row['wilcox_p']:.3e} {sig}")

# ============================================================
# 7. Palindrome density
# ============================================================
print("\n[7] Palindrome density analysis...")

def count_palindromes(seq, min_len=6, max_len=12):
    """Count palindromic sequences (potential R-M targets)."""
    count = 0
    palindromes = []
    for length in range(min_len, max_len + 1, 2):  # Only even lengths are palindromic
        for i in range(len(seq) - length + 1):
            subseq = seq[i:i+length]
            if 'N' in subseq:
                continue
            rc = reverse_complement(subseq)
            if subseq == rc:
                count += 1
                palindromes.append(subseq)
    return count, palindromes

for w in WINDOWS:
    col_seq = f'seq_{w}bp'
    counts = []
    for seq in df[col_seq]:
        c, _ = count_palindromes(seq)
        counts.append(c)
    df[f'palindrome_count_{w}bp'] = counts
    df[f'palindrome_density_{w}bp'] = [c / (len(s) / 1000) if len(s) > 0 else 0
                                         for c, s in zip(counts, df[col_seq])]

# Compare palindrome counts
for w in WINDOWS:
    col = f'palindrome_count_{w}bp'
    exp_vals = df.loc[df['is_exposed'] == 1, col].values
    shi_vals = df.loc[df['is_exposed'] == 0, col].values
    u, p = stats.mannwhitneyu(exp_vals, shi_vals, alternative='two-sided')
    n1, n2 = len(exp_vals), len(shi_vals)
    r_rb = 1 - (2 * u) / (n1 * n2)

    stat_results.append({
        'motif': 'palindrome_6-12bp',
        'window': f'±{w}bp',
        'exposed_mean': np.mean(exp_vals),
        'exposed_median': np.median(exp_vals),
        'shielded_mean': np.mean(shi_vals),
        'shielded_median': np.median(shi_vals),
        'fold_change': np.mean(exp_vals) / np.mean(shi_vals) if np.mean(shi_vals) > 0 else np.inf,
        'wilcox_U': u,
        'wilcox_p': p,
        'rank_biserial_r': r_rb,
        'fisher_OR': np.nan,
        'fisher_p': np.nan,
        'exposed_has_motif_pct': np.nan,
        'shielded_has_motif_pct': np.nan,
    })

    sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "NS"
    print(f"  Palindrome ±{w}bp: exp={np.mean(exp_vals):.2f} vs shi={np.mean(shi_vals):.2f}, p={p:.3e} {sig}")

# Update statistical tests
stat_df = pd.DataFrame(stat_results)
stat_df.to_csv(TABLES_DIR / 'statistical_tests.tsv', sep='\t', index=False)

# ============================================================
# 8. Local motif density gradient (sliding window)
# ============================================================
print("\n[8] Local motif density gradient...")

def count_motif_in_region(chrom_seq, center, start_offset, end_offset, motif, rc_motif):
    """Count motif occurrences in a region relative to center."""
    left = max(0, center + start_offset)
    right = min(len(chrom_seq), center + end_offset)
    region = chrom_seq[left:right]
    count = 0
    for m in [motif, rc_motif]:
        start = 0
        while True:
            pos = region.find(m, start)
            if pos == -1:
                break
            count += 1
            start = pos + 1
    if motif == rc_motif:
        pass  # palindrome already counted correctly
    return count

# Sliding window: 100bp windows from -1000 to +1000
window_size = 100
step = 50
positions = list(range(-1000, 1001 - window_size, step))
gradient_motif = 'GCCGGC'  # Core 6-mer for gradient analysis
gradient_rc = 'GCCGGC'     # Self-complementary palindrome

# Also do AAGCCCG gradient
gradient_motifs = {
    'GCCGGC': 'GCCGGC',
    'AAGCCCG': 'CGGGCTT',
}

gradient_profiles = {}
for gm_name, gm_rc in gradient_motifs.items():
    print(f"\n  Computing gradient for {gm_name}...")

    exposed_profiles = []
    shielded_profiles = []

    for _, row in df.iterrows():
        tss = int(row['tss'])
        strand = row['strand']
        profile = []

        for pos in positions:
            if strand == '+':
                center = tss
                start_off = pos
                end_off = pos + window_size
            else:
                center = tss
                start_off = -(pos + window_size)
                end_off = -pos

            count = count_motif_in_region(chrom_seq, center, start_off, end_off, gm_name, gm_rc)
            profile.append(count)

        if row['is_exposed'] == 1:
            exposed_profiles.append(profile)
        else:
            shielded_profiles.append(profile)

    exposed_mean = np.mean(exposed_profiles, axis=0)
    shielded_mean = np.mean(shielded_profiles, axis=0)
    exposed_se = np.std(exposed_profiles, axis=0) / np.sqrt(len(exposed_profiles))
    shielded_se = np.std(shielded_profiles, axis=0) / np.sqrt(len(shielded_profiles))

    gradient_profiles[gm_name] = {
        'positions': np.array(positions) + window_size / 2,  # Center of window
        'exposed_mean': exposed_mean,
        'shielded_mean': shielded_mean,
        'exposed_se': exposed_se,
        'shielded_se': shielded_se,
    }

    # Find peak difference position
    diff = exposed_mean - shielded_mean
    peak_idx = np.argmax(np.abs(diff))
    peak_pos = positions[peak_idx] + window_size / 2
    print(f"    Peak difference at position {peak_pos:.0f}: exp={exposed_mean[peak_idx]:.4f}, shi={shielded_mean[peak_idx]:.4f}")

# ============================================================
# 9. ROC analysis for sequence features
# ============================================================
print("\n[9] ROC analysis for sequence features...")

y = df['is_exposed'].values

feature_cols = {
    'GCCGGC_count_300bp': 'GCCGGC count ±300bp',
    'AAGCCCG_count_300bp': 'AAGCCCG count ±300bp',
    'TGGCCGGC_count_300bp': 'TGGCCGGC count ±300bp',
    'CCGG_count_300bp': 'CCGG count ±300bp',
    'GC_300bp': 'GC% ±300bp',
    'CpG_OE_300bp': 'CpG O/E ±300bp',
    'palindrome_count_300bp': 'Palindrome count ±300bp',
}

# Add H29 features
if 'nearest_methyl_distance' in df.columns:
    feature_cols['nearest_methyl_distance'] = 'Nearest methyl distance (H29)'
if 'n_methyl_sites_2kb' in df.columns:
    feature_cols['n_methyl_sites_2kb'] = 'Methyl sites 2kb (H29)'

roc_results = []
roc_curves = {}

for col, label in feature_cols.items():
    if col not in df.columns:
        continue
    x = df[col].values
    valid = ~np.isnan(x)
    if valid.sum() < 10:
        continue

    x_valid = x[valid]
    y_valid = y[valid]

    # For distance, lower = more likely exposed, so negate for AUC
    if 'distance' in col:
        x_valid_roc = -x_valid
    elif 'n_methyl' in col:
        x_valid_roc = x_valid  # More sites = more likely exposed
    else:
        x_valid_roc = x_valid

    try:
        auc_val = roc_auc_score(y_valid, x_valid_roc)
        fpr, tpr, thresholds = roc_curve(y_valid, x_valid_roc)
        roc_curves[label] = (fpr, tpr, auc_val)

        # Bootstrap 95% CI
        n_boot = 1000
        boot_aucs = []
        rng = np.random.RandomState(42)
        for _ in range(n_boot):
            idx = rng.choice(len(y_valid), len(y_valid), replace=True)
            if len(np.unique(y_valid[idx])) < 2:
                continue
            try:
                ba = roc_auc_score(y_valid[idx], x_valid_roc[idx])
                boot_aucs.append(ba)
            except:
                pass
        ci_lo = np.percentile(boot_aucs, 2.5) if boot_aucs else np.nan
        ci_hi = np.percentile(boot_aucs, 97.5) if boot_aucs else np.nan

        roc_results.append({
            'feature': col,
            'label': label,
            'AUC': auc_val,
            'AUC_95CI_lower': ci_lo,
            'AUC_95CI_upper': ci_hi,
            'n_samples': len(y_valid),
            'n_exposed': sum(y_valid),
            'n_shielded': sum(1 - y_valid),
        })

        sig_str = "***" if auc_val > 0.7 else "**" if auc_val > 0.6 else "*" if auc_val > 0.55 else "~random"
        print(f"  {label}: AUC={auc_val:.3f} ({ci_lo:.3f}-{ci_hi:.3f}) {sig_str}")
    except Exception as e:
        print(f"  {label}: ERROR - {e}")

roc_df = pd.DataFrame(roc_results).sort_values('AUC', ascending=False)
roc_df.to_csv(TABLES_DIR / 'ROC_analysis.tsv', sep='\t', index=False)

# ============================================================
# 10. Combined logistic regression model
# ============================================================
print("\n[10] Combined logistic regression model...")

# Sequence-only features
seq_features = ['GCCGGC_count_300bp', 'AAGCCCG_count_300bp', 'TGGCCGGC_count_300bp',
                'CCGG_count_300bp', 'GC_300bp', 'palindrome_count_300bp']

# Check which features exist
available_features = [f for f in seq_features if f in df.columns]

X_seq = df[available_features].values
y = df['is_exposed'].values

# Remove NaN rows
valid = ~np.any(np.isnan(X_seq), axis=1)
X_seq_valid = X_seq[valid]
y_valid = y[valid]

# Standardize
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_seq_valid)

# Logistic regression
lr = LogisticRegression(max_iter=1000, random_state=42)
lr.fit(X_scaled, y_valid)

# Predictions
y_pred_proba = lr.predict_proba(X_scaled)[:, 1]
combined_auc = roc_auc_score(y_valid, y_pred_proba)
fpr_comb, tpr_comb, _ = roc_curve(y_valid, y_pred_proba)
roc_curves['Combined LR (sequence)'] = (fpr_comb, tpr_comb, combined_auc)

# Cross-validated AUC
cv_aucs = cross_val_score(lr, X_scaled, y_valid, cv=5, scoring='roc_auc')

print(f"  Combined LR AUC (train): {combined_auc:.3f}")
print(f"  Cross-validated AUC: {np.mean(cv_aucs):.3f} ± {np.std(cv_aucs):.3f}")
print(f"  Feature coefficients:")
for feat, coef in zip(available_features, lr.coef_[0]):
    print(f"    {feat}: {coef:.4f}")

# Bootstrap CI for combined model
n_boot = 1000
boot_aucs_comb = []
rng = np.random.RandomState(42)
for _ in range(n_boot):
    idx = rng.choice(len(y_valid), len(y_valid), replace=True)
    if len(np.unique(y_valid[idx])) < 2:
        continue
    try:
        lr_boot = LogisticRegression(max_iter=1000, random_state=42)
        scaler_boot = StandardScaler()
        X_boot = scaler_boot.fit_transform(X_seq_valid[idx])
        lr_boot.fit(X_boot, y_valid[idx])
        pred_boot = lr_boot.predict_proba(X_boot)[:, 1]
        ba = roc_auc_score(y_valid[idx], pred_boot)
        boot_aucs_comb.append(ba)
    except:
        pass
ci_lo_comb = np.percentile(boot_aucs_comb, 2.5)
ci_hi_comb = np.percentile(boot_aucs_comb, 97.5)
print(f"  Bootstrap 95% CI: {ci_lo_comb:.3f}-{ci_hi_comb:.3f}")

roc_results.append({
    'feature': 'combined_LR_sequence',
    'label': 'Combined LR (sequence features)',
    'AUC': combined_auc,
    'AUC_95CI_lower': ci_lo_comb,
    'AUC_95CI_upper': ci_hi_comb,
    'n_samples': len(y_valid),
    'n_exposed': sum(y_valid),
    'n_shielded': int(sum(1 - y_valid)),
})

roc_results.append({
    'feature': 'combined_LR_sequence_CV',
    'label': 'Combined LR (5-fold CV)',
    'AUC': np.mean(cv_aucs),
    'AUC_95CI_lower': np.mean(cv_aucs) - 1.96 * np.std(cv_aucs),
    'AUC_95CI_upper': np.mean(cv_aucs) + 1.96 * np.std(cv_aucs),
    'n_samples': len(y_valid),
    'n_exposed': sum(y_valid),
    'n_shielded': int(sum(1 - y_valid)),
})

# Update ROC table
roc_df = pd.DataFrame(roc_results).sort_values('AUC', ascending=False)
roc_df.to_csv(TABLES_DIR / 'ROC_analysis.tsv', sep='\t', index=False)

# ============================================================
# 11. Sequence composition at TSS ±50bp (exposed vs shielded)
# ============================================================
print("\n[11] Sequence composition at TSS ±50bp...")

# Extract ±50bp sequences
tss_50bp_seqs = {}
for cls in ['Exposed', 'Shielded']:
    seqs = []
    for _, row in df[df['class'] == cls].iterrows():
        tss = int(row['tss'])
        strand = row['strand']
        seq = extract_tss_region(tss, strand, chrom_seq, 50)
        if len(seq) == 100:  # Only use full-length sequences
            seqs.append(seq)
    tss_50bp_seqs[cls] = seqs
    print(f"  {cls}: {len(seqs)} full-length ±50bp sequences")

# Position-specific nucleotide frequencies
for cls in ['Exposed', 'Shielded']:
    seqs = tss_50bp_seqs[cls]
    n = len(seqs)
    if n == 0:
        continue
    freq_matrix = np.zeros((100, 4))  # 100 positions, 4 nucleotides
    nuc_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for seq in seqs:
        for i, base in enumerate(seq):
            if base in nuc_idx:
                freq_matrix[i, nuc_idx[base]] += 1
    freq_matrix /= n

    gc_profile = freq_matrix[:, 1] + freq_matrix[:, 2]  # C + G

    if cls == 'Exposed':
        exposed_gc_profile = gc_profile
        exposed_freq = freq_matrix
    else:
        shielded_gc_profile = gc_profile
        shielded_freq = freq_matrix

# ============================================================
# 12. Save per-gene feature table
# ============================================================
print("\n[12] Saving per-gene feature table...")

# Select columns for output
output_cols = ['locus_tag', 'gene_name', 'old_locus_tag', 'product', 'tf_family',
               'start', 'end', 'strand', 'tss', 'region', 'is_exposed', 'class']

# Add motif counts
for motif_name in MOTIFS.keys():
    for w in WINDOWS:
        output_cols.append(f'{motif_name}_count_{w}bp')
        output_cols.append(f'{motif_name}_density_{w}bp')

# Add GC and palindrome
for w in WINDOWS:
    output_cols.extend([f'GC_{w}bp', f'CpG_OE_{w}bp', f'palindrome_count_{w}bp', f'palindrome_density_{w}bp'])

# Add H29 features if present
for col in ['baseMean', 'nearest_methyl_distance', 'n_methyl_sites_2kb', 'gene_length']:
    if col in df.columns:
        output_cols.append(col)

# Filter to existing columns
output_cols = [c for c in output_cols if c in df.columns]
output_df = df[output_cols].copy()
output_df.to_csv(TABLES_DIR / 'TSS_sequence_features.tsv', sep='\t', index=False)
print(f"  Saved: {TABLES_DIR / 'TSS_sequence_features.tsv'} ({len(output_df)} genes, {len(output_cols)} columns)")

# ============================================================
# 13. FIGURES
# ============================================================
print("\n[13] Generating figures...")

# Color scheme
COLOR_EXPOSED = '#E74C3C'
COLOR_SHIELDED = '#3498DB'
COLOR_BG = '#95A5A6'

# --- Figure 1: Motif density comparison ---
print("  Figure 1: Motif density comparison...")
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

for idx, motif_name in enumerate(MOTIFS.keys()):
    ax = axes[idx // 2, idx % 2]
    col = f'{motif_name}_count_300bp'

    exp_vals = df.loc[df['is_exposed'] == 1, col].values
    shi_vals = df.loc[df['is_exposed'] == 0, col].values

    # Boxplot comparison
    bp = ax.boxplot([exp_vals, shi_vals],
                     labels=['Exposed\n(n={})'.format(n_exposed), 'Shielded\n(n={})'.format(n_shielded)],
                     patch_artist=True,
                     widths=0.5)
    bp['boxes'][0].set_facecolor(COLOR_EXPOSED)
    bp['boxes'][0].set_alpha(0.5)
    bp['boxes'][1].set_facecolor(COLOR_SHIELDED)
    bp['boxes'][1].set_alpha(0.5)

    # Overlay points
    jitter_exp = np.random.RandomState(42).normal(0, 0.05, len(exp_vals))
    jitter_shi = np.random.RandomState(42).normal(0, 0.05, len(shi_vals))
    ax.scatter(1 + jitter_exp, exp_vals, c=COLOR_EXPOSED, alpha=0.3, s=10, zorder=5)
    ax.scatter(2 + jitter_shi, shi_vals, c=COLOR_SHIELDED, alpha=0.1, s=5, zorder=5)

    # p-value annotation
    u, p = stats.mannwhitneyu(exp_vals, shi_vals, alternative='two-sided')
    fold = np.mean(exp_vals) / np.mean(shi_vals) if np.mean(shi_vals) > 0 else float('inf')
    ax.set_title(f'{motif_name} (±300bp)\nfold={fold:.2f}, p={p:.2e}', fontsize=11)
    ax.set_ylabel('Motif count')

fig.suptitle('H30: R-M Motif Counts in TSS ±300bp Region\nExposed vs Shielded Regulators',
             fontsize=13, fontweight='bold')
plt.tight_layout()
for fmt in ['pdf', 'svg']:
    fig.savefig(FIGURES_DIR / f'motif_density_comparison.{fmt}', dpi=300, bbox_inches='tight')
plt.close()

# --- Figure 2: Motif gradient profile ---
print("  Figure 2: Motif gradient profile...")
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

for idx, (gm_name, prof) in enumerate(gradient_profiles.items()):
    ax = axes[idx]
    pos = prof['positions']

    ax.plot(pos, prof['exposed_mean'], color=COLOR_EXPOSED, linewidth=2, label=f'Exposed (n={n_exposed})')
    ax.fill_between(pos,
                     prof['exposed_mean'] - 1.96 * prof['exposed_se'],
                     prof['exposed_mean'] + 1.96 * prof['exposed_se'],
                     color=COLOR_EXPOSED, alpha=0.2)

    ax.plot(pos, prof['shielded_mean'], color=COLOR_SHIELDED, linewidth=2, label=f'Shielded (n={n_shielded})')
    ax.fill_between(pos,
                     prof['shielded_mean'] - 1.96 * prof['shielded_se'],
                     prof['shielded_mean'] + 1.96 * prof['shielded_se'],
                     color=COLOR_SHIELDED, alpha=0.2)

    ax.axvline(0, color='black', linestyle='--', alpha=0.5, label='TSS')
    ax.set_xlabel('Distance from TSS (bp)')
    ax.set_ylabel(f'Mean {gm_name} count per 100bp window')
    ax.set_title(f'{gm_name} density gradient around TSS')
    ax.legend(fontsize=9)

fig.suptitle('H30: Motif Density Gradient Around TSS\n(Sliding 100bp windows, step=50bp)',
             fontsize=13, fontweight='bold')
plt.tight_layout()
for fmt in ['pdf', 'svg']:
    fig.savefig(FIGURES_DIR / f'motif_gradient_profile.{fmt}', dpi=300, bbox_inches='tight')
plt.close()

# --- Figure 3: GC content comparison ---
print("  Figure 3: GC content comparison...")
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

for idx, w in enumerate(WINDOWS):
    ax = axes[idx]
    col = f'GC_{w}bp'

    exp_vals = df.loc[df['is_exposed'] == 1, col].values
    shi_vals = df.loc[df['is_exposed'] == 0, col].values

    ax.hist(shi_vals, bins=30, alpha=0.5, color=COLOR_SHIELDED, label=f'Shielded (n={n_shielded})', density=True)
    ax.hist(exp_vals, bins=20, alpha=0.5, color=COLOR_EXPOSED, label=f'Exposed (n={n_exposed})', density=True)

    ax.axvline(np.median(exp_vals), color=COLOR_EXPOSED, linestyle='--', linewidth=2)
    ax.axvline(np.median(shi_vals), color=COLOR_SHIELDED, linestyle='--', linewidth=2)

    u, p = stats.mannwhitneyu(exp_vals, shi_vals, alternative='two-sided')
    ax.set_title(f'±{w}bp window\nexp median={np.median(exp_vals):.1f}%, shi median={np.median(shi_vals):.1f}%\np={p:.2e}')
    ax.set_xlabel('GC content (%)')
    ax.set_ylabel('Density')
    ax.legend(fontsize=9)

fig.suptitle('H30: GC Content in TSS-Proximal Regions\nExposed vs Shielded Regulators',
             fontsize=13, fontweight='bold')
plt.tight_layout()
for fmt in ['pdf', 'svg']:
    fig.savefig(FIGURES_DIR / f'GC_content_comparison.{fmt}', dpi=300, bbox_inches='tight')
plt.close()

# --- Figure 4: ROC curves ---
print("  Figure 4: ROC curves...")
fig, ax = plt.subplots(1, 1, figsize=(8, 8))

# Sort by AUC for legend ordering
sorted_rocs = sorted(roc_curves.items(), key=lambda x: x[1][2], reverse=True)

colors = plt.cm.tab10(np.linspace(0, 1, len(sorted_rocs)))
for i, (label, (fpr, tpr, auc_val)) in enumerate(sorted_rocs):
    ls = '-' if auc_val > 0.6 else '--'
    lw = 2.5 if 'distance' in label.lower() or 'Combined' in label else 1.5
    ax.plot(fpr, tpr, color=colors[i], linestyle=ls, linewidth=lw,
            label=f'{label} (AUC={auc_val:.3f})')

ax.plot([0, 1], [0, 1], 'k--', alpha=0.5, label='Random (AUC=0.500)')
ax.set_xlabel('False Positive Rate', fontsize=12)
ax.set_ylabel('True Positive Rate', fontsize=12)
ax.set_title('H30: ROC Analysis - Sequence Features\nPredicting Exposed vs Shielded Regulators',
             fontsize=13, fontweight='bold')
ax.legend(loc='lower right', fontsize=8)
ax.set_xlim([0, 1])
ax.set_ylim([0, 1])

plt.tight_layout()
for fmt in ['pdf', 'svg']:
    fig.savefig(FIGURES_DIR / f'ROC_sequence_features.{fmt}', dpi=300, bbox_inches='tight')
plt.close()

# --- Figure 5: Comprehensive summary ---
print("  Figure 5: Comprehensive summary...")
fig = plt.figure(figsize=(18, 14))
gs = gridspec.GridSpec(3, 3, hspace=0.4, wspace=0.35)

# Panel A: Motif counts (GCCGGC ±300bp)
ax1 = fig.add_subplot(gs[0, 0])
col = 'GCCGGC_count_300bp'
exp_vals = df.loc[df['is_exposed'] == 1, col].values
shi_vals = df.loc[df['is_exposed'] == 0, col].values
bp = ax1.boxplot([exp_vals, shi_vals], labels=['Exposed', 'Shielded'], patch_artist=True, widths=0.5)
bp['boxes'][0].set_facecolor(COLOR_EXPOSED); bp['boxes'][0].set_alpha(0.5)
bp['boxes'][1].set_facecolor(COLOR_SHIELDED); bp['boxes'][1].set_alpha(0.5)
u, p = stats.mannwhitneyu(exp_vals, shi_vals, alternative='two-sided')
ax1.set_title(f'A) GCCGGC ±300bp\np={p:.2e}', fontsize=10)
ax1.set_ylabel('Motif count')

# Panel B: GC content
ax2 = fig.add_subplot(gs[0, 1])
col = 'GC_300bp'
exp_vals = df.loc[df['is_exposed'] == 1, col].values
shi_vals = df.loc[df['is_exposed'] == 0, col].values
ax2.hist(shi_vals, bins=30, alpha=0.5, color=COLOR_SHIELDED, label='Shielded', density=True)
ax2.hist(exp_vals, bins=20, alpha=0.5, color=COLOR_EXPOSED, label='Exposed', density=True)
u, p = stats.mannwhitneyu(exp_vals, shi_vals, alternative='two-sided')
ax2.set_title(f'B) GC% ±300bp\np={p:.2e}', fontsize=10)
ax2.set_xlabel('GC%')
ax2.legend(fontsize=8)

# Panel C: Palindrome count
ax3 = fig.add_subplot(gs[0, 2])
col = 'palindrome_count_300bp'
exp_vals = df.loc[df['is_exposed'] == 1, col].values
shi_vals = df.loc[df['is_exposed'] == 0, col].values
bp = ax3.boxplot([exp_vals, shi_vals], labels=['Exposed', 'Shielded'], patch_artist=True, widths=0.5)
bp['boxes'][0].set_facecolor(COLOR_EXPOSED); bp['boxes'][0].set_alpha(0.5)
bp['boxes'][1].set_facecolor(COLOR_SHIELDED); bp['boxes'][1].set_alpha(0.5)
u, p = stats.mannwhitneyu(exp_vals, shi_vals, alternative='two-sided')
ax3.set_title(f'C) Palindromes ±300bp\np={p:.2e}', fontsize=10)
ax3.set_ylabel('Palindrome count (>=6bp)')

# Panel D: GCCGGC gradient
ax4 = fig.add_subplot(gs[1, 0:2])
prof = gradient_profiles['GCCGGC']
pos = prof['positions']
ax4.plot(pos, prof['exposed_mean'], color=COLOR_EXPOSED, linewidth=2, label=f'Exposed (n={n_exposed})')
ax4.fill_between(pos, prof['exposed_mean'] - 1.96*prof['exposed_se'],
                 prof['exposed_mean'] + 1.96*prof['exposed_se'], color=COLOR_EXPOSED, alpha=0.2)
ax4.plot(pos, prof['shielded_mean'], color=COLOR_SHIELDED, linewidth=2, label=f'Shielded (n={n_shielded})')
ax4.fill_between(pos, prof['shielded_mean'] - 1.96*prof['shielded_se'],
                 prof['shielded_mean'] + 1.96*prof['shielded_se'], color=COLOR_SHIELDED, alpha=0.2)
ax4.axvline(0, color='black', linestyle='--', alpha=0.5)
ax4.set_xlabel('Distance from TSS (bp)')
ax4.set_ylabel('Mean GCCGGC count / 100bp')
ax4.set_title('D) GCCGGC motif density gradient around TSS', fontsize=10)
ax4.legend(fontsize=9)

# Panel E: ROC curves (top 5)
ax5 = fig.add_subplot(gs[1, 2])
top_rocs = sorted_rocs[:6]  # Top 6 features
colors_roc = plt.cm.tab10(np.linspace(0, 1, len(top_rocs)))
for i, (label, (fpr, tpr, auc_val)) in enumerate(top_rocs):
    short_label = label.split('(')[0].strip()[:30]
    ax5.plot(fpr, tpr, color=colors_roc[i], linewidth=1.5,
             label=f'{short_label} ({auc_val:.3f})')
ax5.plot([0, 1], [0, 1], 'k--', alpha=0.5)
ax5.set_xlabel('FPR')
ax5.set_ylabel('TPR')
ax5.set_title('E) ROC: sequence features', fontsize=10)
ax5.legend(fontsize=7, loc='lower right')

# Panel F: GC profile at TSS ±50bp
ax6 = fig.add_subplot(gs[2, 0:2])
positions_50 = np.arange(-50, 50)
if 'exposed_gc_profile' in dir() and 'shielded_gc_profile' in dir():
    ax6.plot(positions_50, exposed_gc_profile * 100, color=COLOR_EXPOSED, linewidth=1.5,
             label=f'Exposed (n={len(tss_50bp_seqs["Exposed"])})')
    ax6.plot(positions_50, shielded_gc_profile * 100, color=COLOR_SHIELDED, linewidth=1.5,
             label=f'Shielded (n={len(tss_50bp_seqs["Shielded"])})')
    ax6.axvline(0, color='black', linestyle='--', alpha=0.5, label='TSS')
    ax6.set_xlabel('Position relative to TSS (bp)')
    ax6.set_ylabel('GC frequency (%)')
    ax6.set_title('F) Position-specific GC content at TSS ±50bp', fontsize=10)
    ax6.legend(fontsize=9)

# Panel G: Dinucleotide differences (top differential)
ax7 = fig.add_subplot(gs[2, 2])
top_dinuc = dinuc_stat.head(8)
colors_dinuc = ['#E74C3C' if d > 0 else '#3498DB' for d in top_dinuc['diff']]
ax7.barh(range(len(top_dinuc)), top_dinuc['diff'].values, color=colors_dinuc, alpha=0.7)
ax7.set_yticks(range(len(top_dinuc)))
ax7.set_yticklabels(top_dinuc['dinucleotide'].values)
for i, (_, row) in enumerate(top_dinuc.iterrows()):
    sig = "*" if row['wilcox_p'] < 0.05 else ""
    ax7.text(row['diff'] + 0.0001 * np.sign(row['diff']), i,
             f"p={row['wilcox_p']:.1e}{sig}", fontsize=7, va='center')
ax7.set_xlabel('Freq difference (Exposed - Shielded)')
ax7.set_title('G) Top dinucleotide differences ±300bp', fontsize=10)
ax7.axvline(0, color='black', linewidth=0.5)

fig.suptitle('H30: TSS-Proximal DNA Sequence Features as Determinants of Protection Zone',
             fontsize=14, fontweight='bold', y=0.98)
plt.tight_layout(rect=[0, 0, 1, 0.96])
for fmt in ['pdf', 'svg']:
    fig.savefig(FIGURES_DIR / f'H30_comprehensive_summary.{fmt}', dpi=300, bbox_inches='tight')
plt.close()

# ============================================================
# 14. Final summary
# ============================================================
print("\n" + "=" * 70)
print("H30 ANALYSIS COMPLETE - SUMMARY")
print("=" * 70)

# Best sequence feature
best_seq = roc_df[~roc_df['feature'].str.contains('combined|CV|nearest|n_methyl')].iloc[0] if len(roc_df) > 0 else None
best_overall = roc_df.iloc[0] if len(roc_df) > 0 else None

print(f"\nKey Results:")
print(f"  Genes analyzed: {len(df)} ({n_exposed} exposed, {n_shielded} shielded)")

# Motif comparison summary
for motif_name in MOTIFS.keys():
    col = f'{motif_name}_count_300bp'
    exp_mean = df.loc[df['is_exposed'] == 1, col].mean()
    shi_mean = df.loc[df['is_exposed'] == 0, col].mean()
    fold = exp_mean / shi_mean if shi_mean > 0 else float('inf')
    u, p = stats.mannwhitneyu(
        df.loc[df['is_exposed'] == 1, col].values,
        df.loc[df['is_exposed'] == 0, col].values,
        alternative='two-sided'
    )
    sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "NS"
    print(f"  {motif_name} ±300bp: fold={fold:.2f}, p={p:.2e} {sig}")

# GC content
exp_gc = df.loc[df['is_exposed'] == 1, 'GC_300bp'].mean()
shi_gc = df.loc[df['is_exposed'] == 0, 'GC_300bp'].mean()
u, p_gc = stats.mannwhitneyu(
    df.loc[df['is_exposed'] == 1, 'GC_300bp'].values,
    df.loc[df['is_exposed'] == 0, 'GC_300bp'].values,
    alternative='two-sided'
)
print(f"  GC% ±300bp: exp={exp_gc:.1f}%, shi={shi_gc:.1f}%, p={p_gc:.2e}")

# Palindromes
exp_pal = df.loc[df['is_exposed'] == 1, 'palindrome_count_300bp'].mean()
shi_pal = df.loc[df['is_exposed'] == 0, 'palindrome_count_300bp'].mean()
u, p_pal = stats.mannwhitneyu(
    df.loc[df['is_exposed'] == 1, 'palindrome_count_300bp'].values,
    df.loc[df['is_exposed'] == 0, 'palindrome_count_300bp'].values,
    alternative='two-sided'
)
print(f"  Palindromes ±300bp: exp={exp_pal:.1f}, shi={shi_pal:.1f}, p={p_pal:.2e}")

print(f"\nROC Analysis (top features):")
for _, row in roc_df.head(8).iterrows():
    print(f"  {row['label']}: AUC={row['AUC']:.3f} ({row['AUC_95CI_lower']:.3f}-{row['AUC_95CI_upper']:.3f})")

print(f"\nCombined LR model:")
print(f"  Training AUC: {combined_auc:.3f} ({ci_lo_comb:.3f}-{ci_hi_comb:.3f})")
print(f"  CV AUC: {np.mean(cv_aucs):.3f} +/- {np.std(cv_aucs):.3f}")

# H29 comparison
h29_ref_auc = 0.917
print(f"\nComparison to H29 nearest_methyl_distance AUC={h29_ref_auc:.3f}:")
if best_seq is not None:
    print(f"  Best sequence feature: {best_seq['label']}, AUC={best_seq['AUC']:.3f}")
    gap = h29_ref_auc - best_seq['AUC']
    print(f"  AUC gap to H29: {gap:.3f}")
    print(f"  Sequence features explain {best_seq['AUC']/h29_ref_auc*100:.1f}% of H29 discrimination")

print(f"\nConclusion:")
if best_seq is not None and best_seq['AUC'] > 0.7:
    print(f"  SUPPORTED: Sequence features substantially predict exposed/shielded classification")
elif best_seq is not None and best_seq['AUC'] > 0.6:
    print(f"  PARTIAL: Sequence features contribute modestly to exposed/shielded classification")
else:
    print(f"  REJECTED: TSS-proximal sequence features do NOT substantially predict protection zone")
    print(f"  Protection zone is determined by factors beyond DNA sequence composition")

print(f"\nOutput files:")
print(f"  {TABLES_DIR / 'TSS_sequence_features.tsv'}")
print(f"  {TABLES_DIR / 'statistical_tests.tsv'}")
print(f"  {TABLES_DIR / 'ROC_analysis.tsv'}")
print(f"  {TABLES_DIR / 'dinucleotide_frequencies.tsv'}")
for f in FIGURES_DIR.glob('*.pdf'):
    print(f"  {f}")

print("\nDone.")
