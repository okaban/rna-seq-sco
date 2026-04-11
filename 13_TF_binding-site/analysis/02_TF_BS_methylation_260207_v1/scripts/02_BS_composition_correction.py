#!/usr/bin/env python3
"""
TF Binding Site methylation depletion — base composition correction

The original analysis (01 script, fold=0.66, p=6.8e-5 with 782 Tier 1 BS)
assumed uniform methylation distribution across the genome.

This script tests whether the depletion survives correction for:
  M0: Uniform (original, for reference)
  M1: GC-content adjusted (base-specific methylation rates)
  M2: Motif-count adjusted (CCGG, AAGCCCG, GATC occurrence rates in BS)
  M3: Permutation with random + GC-matched random regions

NOTE: Coordinate resolution from 01 script is replicated here
  - FIMO: absolute genomic coordinates
  - RegPrecise: relative positions -> absolute using target gene TSS
  - ZorroAranda Curated Strong: target gene promoter region (-300 to +50)
"""

import pandas as pd
import numpy as np
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# === Paths ===
PROJECT = Path("/Users/okaban/bioinfo/rna-seq")
BS_DIR = PROJECT / "13_TF_binding-site/analysis/01_master_TF_list_260206_v1"
BS_FILE = BS_DIR / "master_TF_binding_sites_M145_ARCHIVED_260224.tsv"
EPIGENOME = PROJECT / "11_epigenome_integration"
METHYL_4mC = EPIGENOME / "analysis/23_expanded_motif_search/4mC_final_census.csv"
METHYL_6mA = EPIGENOME / "analysis/23_expanded_motif_search/6mA_final_census.csv"
REF_FASTA = Path("/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa")
OUT_DIR = PROJECT / "13_TF_binding-site/analysis/02_TF_BS_methylation_260207_v1"

GENOME_SIZE = 8_667_507
TIER1_SOURCES = {'FIMO_curated_motif', 'RegPrecise', 'ZorroAranda2022_Curated_Strong'}
CONTIG_MAP = {'chromosome': 'NC_003888.3'}
N_PERMUTATIONS = 10000
np.random.seed(42)

# === Load reference genome ===
print("Loading reference genome...")
ref_record = next(SeqIO.parse(REF_FASTA, "fasta"))
genome_seq = str(ref_record.seq).upper()
print(f"  Genome: {ref_record.id}, {len(genome_seq):,} bp")

# === Load and resolve BS coordinates (same logic as 01 script) ===
print("\nLoading and resolving binding site coordinates...")
bs_raw = pd.read_csv(BS_FILE, sep='\t')

# Gene info for coordinate conversion
gene_info = pd.read_csv(BS_DIR / 'intermediate' / 'M145_gene_basic_info.tsv', sep='\t')
gene_coord = {}
sco_to_geneid = {}
for _, g in gene_info.iterrows():
    gid = g.get('gene_id', '')
    sco = g.get('old_locus_tag', g.get('SCO_ID', ''))
    start = g.get('start', np.nan)
    end = g.get('end', np.nan)
    strand = g.get('strand', '+')
    if pd.notna(start) and pd.notna(end):
        gene_coord[gid] = (int(start), int(end), strand)
        if sco:
            sco_to_geneid[sco] = gid
            gene_coord[sco] = (int(start), int(end), strand)

# TSS table
tss_df = pd.read_csv(EPIGENOME / 'analysis/18_tss_analyses/comprehensive_tss_table.csv')
tss_map = dict(zip(tss_df['gene_id'], tss_df['tss']))
for _, t in tss_df.iterrows():
    olt = t.get('old_locus_tag', '')
    if pd.notna(olt) and olt:
        tss_map[olt] = t['tss']

# Resolve absolute coordinates
resolved_bs = []
for _, row in bs_raw.iterrows():
    src = row['BS_source']
    if src not in TIER1_SOURCES:
        continue

    rec = {
        'TF_name': row['TF_name'],
        'BS_source': src,
        'abs_start': np.nan,
        'abs_end': np.nan,
    }

    tg_sco = row.get('TG_SCO_ID', '')
    tg_gid = row.get('TG_gene_id', '')
    bs_seq = row.get('BS_sequence', '')

    if src == 'FIMO_curated_motif' and pd.notna(row['BS_start']):
        rec['abs_start'] = int(row['BS_start'])
        rec['abs_end'] = int(row['BS_end'])

    elif src == 'RegPrecise' and pd.notna(row['BS_start']):
        rel_pos = int(row['BS_start'])
        tg_key = tg_gid if tg_gid else tg_sco
        tss_val = tss_map.get(tg_key) or tss_map.get(tg_sco)
        if tss_val and pd.notna(tss_val):
            tss_val = int(tss_val)
            coord = gene_coord.get(tg_key) or gene_coord.get(tg_sco)
            if coord:
                strand = coord[2]
                seq_len = len(bs_seq) if pd.notna(bs_seq) and bs_seq else 18
                if strand == '+':
                    abs_s = tss_val + rel_pos
                    abs_e = abs_s + seq_len - 1
                else:
                    abs_e = tss_val - rel_pos
                    abs_s = abs_e - seq_len + 1
                rec['abs_start'] = max(1, abs_s)
                rec['abs_end'] = abs_e

    elif src == 'ZorroAranda2022_Curated_Strong':
        tg_key = tg_gid if tg_gid else tg_sco
        tss_val = tss_map.get(tg_key) or tss_map.get(tg_sco)
        coord = gene_coord.get(tg_key) or gene_coord.get(tg_sco)
        if tss_val and pd.notna(tss_val) and coord:
            tss_val = int(tss_val)
            strand = coord[2]
            if strand == '+':
                rec['abs_start'] = max(1, tss_val - 300)
                rec['abs_end'] = tss_val + 50
            else:
                rec['abs_start'] = max(1, tss_val - 50)
                rec['abs_end'] = tss_val + 300

    resolved_bs.append(rec)

bs_df = pd.DataFrame(resolved_bs)
bs_df = bs_df.dropna(subset=['abs_start', 'abs_end'])
bs_df['abs_start'] = bs_df['abs_start'].astype(int)
bs_df['abs_end'] = bs_df['abs_end'].astype(int)
bs_df = bs_df[(bs_df['abs_start'] >= 0) & (bs_df['abs_end'] <= GENOME_SIZE) & (bs_df['abs_start'] < bs_df['abs_end'])]
print(f"  Resolved Tier 1 BS: {len(bs_df)} sites")
for src, cnt in bs_df['BS_source'].value_counts().items():
    print(f"    {src}: {cnt}")

# === Load methylation data ===
print("\nLoading methylation sites...")
df_4mC = pd.read_csv(METHYL_4mC)
df_6mA = pd.read_csv(METHYL_6mA)
pos_4mC = df_4mC.drop_duplicates(subset=['chrom', 'position'])
pos_6mA = df_6mA.drop_duplicates(subset=['chrom', 'position'])
pos_4mC_set = set(pos_4mC[pos_4mC['chrom'] == 'NC_003888.3']['position'].values)
pos_6mA_set = set(pos_6mA[pos_6mA['chrom'] == 'NC_003888.3']['position'].values)
pos_all_set = pos_4mC_set | pos_6mA_set
print(f"  4mC: {len(pos_4mC_set)}, 6mA: {len(pos_6mA_set)}, all: {len(pos_all_set)}")

# === Helper functions ===
def count_methyl_fast(regions, sorted_arr):
    """Fast counting using sorted array and searchsorted."""
    count = 0
    for s, e in regions:
        left = np.searchsorted(sorted_arr, s, side='left')
        right = np.searchsorted(sorted_arr, e, side='right')
        count += (right - left)
    return count

def extract_sequences(regions, genome):
    """Extract sequences for regions."""
    seqs = []
    for s, e in regions:
        s = max(0, s)
        e = min(len(genome), e)
        seqs.append(genome[s:e])
    return seqs

def gc_content(seq):
    if len(seq) == 0:
        return 0
    return (seq.count('G') + seq.count('C')) / len(seq)

def count_motif(seq, motif):
    """Count occurrences of motif in sequence (both strands). Uses str.count()."""
    motif_rc = str(Seq(motif).reverse_complement())
    if motif == motif_rc:
        return seq.count(motif)
    return seq.count(motif) + seq.count(motif_rc)

# === Prepare BS regions ===
bs_regions = list(zip(bs_df['abs_start'].values, bs_df['abs_end'].values))
bs_total_width = sum(e - s for s, e in bs_regions)
print(f"\nBS total width: {bs_total_width:,} bp ({len(bs_regions)} regions)")

# Sort methylation positions for fast lookup
arr_4mC = np.sort(list(pos_4mC_set))
arr_6mA = np.sort(list(pos_6mA_set))
arr_all = np.sort(list(pos_all_set))

# Observed counts
obs_4mC = count_methyl_fast(bs_regions, arr_4mC)
obs_6mA = count_methyl_fast(bs_regions, arr_6mA)
obs_all = count_methyl_fast(bs_regions, arr_all)

print(f"\nObserved methylation in Tier 1 BS regions:")
print(f"  4mC: {obs_4mC}")
print(f"  6mA: {obs_6mA}")
print(f"  Any: {obs_all}")

# === Method 0: Original (uniform expectation) ===
print("\n" + "=" * 70)
print("METHOD 0: Original uniform expectation (for reference)")
print("=" * 70)

exp_4mC_uniform = len(pos_4mC_set) * bs_total_width / GENOME_SIZE
exp_6mA_uniform = len(pos_6mA_set) * bs_total_width / GENOME_SIZE
exp_all_uniform = len(pos_all_set) * bs_total_width / GENOME_SIZE

for label, obs, exp in [("4mC", obs_4mC, exp_4mC_uniform),
                         ("6mA", obs_6mA, exp_6mA_uniform),
                         ("Any", obs_all, exp_all_uniform)]:
    fold = obs / exp if exp > 0 else float('inf')
    p = stats.poisson.cdf(obs, exp)
    print(f"  {label}: obs={obs}, exp={exp:.1f}, fold={fold:.3f}, p(depletion)={p:.2e}")

# === Method 1: GC-content adjusted expectation ===
print("\n" + "=" * 70)
print("METHOD 1: GC-content adjusted expectation")
print("=" * 70)

genome_gc = gc_content(genome_seq)
print(f"  Genome GC: {genome_gc:.4f}")

bs_seqs = extract_sequences(bs_regions, genome_seq)
bs_gc_values = [gc_content(s) for s in bs_seqs]
bs_gc_mean = np.mean(bs_gc_values)
print(f"  BS mean GC: {bs_gc_mean:.4f}")
print(f"  BS GC range: [{min(bs_gc_values):.3f}, {max(bs_gc_values):.3f}]")

# Base-specific methylation rates
genome_C_count = genome_seq.count('C') + genome_seq.count('G')
genome_A_count = genome_seq.count('A') + genome_seq.count('T')
rate_4mC_per_C = len(pos_4mC_set) / genome_C_count
rate_6mA_per_A = len(pos_6mA_set) / genome_A_count
print(f"  Genome C+G: {genome_C_count:,}, A+T: {genome_A_count:,}")
print(f"  4mC rate per C/G: {rate_4mC_per_C:.6f}")
print(f"  6mA rate per A/T: {rate_6mA_per_A:.6f}")

# GC-adjusted expectation per BS
exp_4mC_gc = 0
exp_6mA_gc = 0
for seq in bs_seqs:
    n_cg = sum(1 for c in seq if c in 'CG')
    n_at = sum(1 for c in seq if c in 'AT')
    exp_4mC_gc += n_cg * rate_4mC_per_C
    exp_6mA_gc += n_at * rate_6mA_per_A
exp_all_gc = exp_4mC_gc + exp_6mA_gc

for label, obs, exp in [("4mC", obs_4mC, exp_4mC_gc),
                         ("6mA", obs_6mA, exp_6mA_gc),
                         ("Any", obs_all, exp_all_gc)]:
    fold = obs / exp if exp > 0 else float('inf')
    p = stats.poisson.cdf(obs, exp)
    print(f"  {label}: obs={obs}, exp={exp:.1f}, fold={fold:.3f}, p(depletion)={p:.2e}")

# === Method 2: Motif-count adjusted expectation ===
print("\n" + "=" * 70)
print("METHOD 2: Motif-count adjusted expectation")
print("=" * 70)

MOTIFS = ['CCGG', 'AAGCCCG', 'GATC']
genome_motif_counts = {}
for motif in MOTIFS:
    genome_motif_counts[motif] = count_motif(genome_seq, motif)
    print(f"  Genome {motif}: {genome_motif_counts[motif]:,}")

# Methylation rate per motif from census data
print("  Loading motif assignments from census...")
motif_4mC_counts = {}
if 'final_motif' in df_4mC.columns:
    for motif in MOTIFS:
        n = len(df_4mC[df_4mC['final_motif'].str.contains(motif, na=False)].drop_duplicates(subset=['chrom', 'position']))
        motif_4mC_counts[motif] = n
        print(f"    4mC at {motif}: {n} sites")

motif_6mA_counts = {}
if 'final_motif' in df_6mA.columns:
    for motif in MOTIFS:
        n = len(df_6mA[df_6mA['final_motif'].str.contains(motif, na=False)].drop_duplicates(subset=['chrom', 'position']))
        motif_6mA_counts[motif] = n
        print(f"    6mA at {motif}: {n} sites")

rates = {}
for motif in MOTIFS:
    gcount = genome_motif_counts[motif]
    if gcount > 0:
        r4 = motif_4mC_counts.get(motif, 0) / gcount
        r6 = motif_6mA_counts.get(motif, 0) / gcount
        rates[motif] = (r4, r6)
        print(f"  Rate {motif}: 4mC={r4:.4f}, 6mA={r6:.4f} per occurrence")

# Count motifs in BS regions
bs_motif_counts = {m: 0 for m in MOTIFS}
for seq in bs_seqs:
    for motif in MOTIFS:
        bs_motif_counts[motif] += count_motif(seq, motif)

for motif in MOTIFS:
    print(f"  BS {motif}: {bs_motif_counts[motif]} occurrences")

# Expected methylation from motif counts
exp_4mC_motif = sum(bs_motif_counts[m] * rates[m][0] for m in MOTIFS if m in rates)
exp_6mA_motif = sum(bs_motif_counts[m] * rates[m][1] for m in MOTIFS if m in rates)
# Add residual (non-motif methylation) proportional to width
residual_4mC = len(pos_4mC_set) - sum(motif_4mC_counts.get(m, 0) for m in MOTIFS)
residual_6mA = len(pos_6mA_set) - sum(motif_6mA_counts.get(m, 0) for m in MOTIFS)
exp_4mC_motif += residual_4mC * bs_total_width / GENOME_SIZE
exp_6mA_motif += residual_6mA * bs_total_width / GENOME_SIZE
exp_all_motif = exp_4mC_motif + exp_6mA_motif

for label, obs, exp in [("4mC", obs_4mC, exp_4mC_motif),
                         ("6mA", obs_6mA, exp_6mA_motif),
                         ("Any", obs_all, exp_all_motif)]:
    fold = obs / exp if exp > 0 else float('inf')
    p = stats.poisson.cdf(obs, exp)
    print(f"  {label}: obs={obs}, exp={exp:.1f}, fold={fold:.3f}, p(depletion)={p:.2e}")

# === Method 3: Permutation with random regions ===
print("\n" + "=" * 70)
print(f"METHOD 3: Permutation test ({N_PERMUTATIONS:,} iterations)")
print("=" * 70)

bs_lengths = [e - s for s, e in bs_regions]

def generate_random_regions(lengths, genome_size, margin=200):
    """Generate random genomic regions with same length distribution."""
    regions = []
    for l in lengths:
        start = np.random.randint(margin, genome_size - l - margin)
        regions.append((start, start + l))
    return regions

print("  Running random permutations...")
perm_counts_4mC = []
perm_counts_6mA = []
perm_counts_all = []

for i in range(N_PERMUTATIONS):
    rand_regions = generate_random_regions(bs_lengths, GENOME_SIZE)
    c4 = count_methyl_fast(rand_regions, arr_4mC)
    c6 = count_methyl_fast(rand_regions, arr_6mA)
    ca = count_methyl_fast(rand_regions, arr_all)
    perm_counts_4mC.append(c4)
    perm_counts_6mA.append(c6)
    perm_counts_all.append(ca)
    if (i + 1) % 2000 == 0:
        print(f"    {i + 1}/{N_PERMUTATIONS}")

perm_counts_4mC = np.array(perm_counts_4mC)
perm_counts_6mA = np.array(perm_counts_6mA)
perm_counts_all = np.array(perm_counts_all)

for label, obs, perm in [("4mC", obs_4mC, perm_counts_4mC),
                          ("6mA", obs_6mA, perm_counts_6mA),
                          ("Any", obs_all, perm_counts_all)]:
    p_depletion = np.mean(perm <= obs)
    fold = obs / np.mean(perm) if np.mean(perm) > 0 else float('inf')
    print(f"  {label}: obs={obs}, perm_mean={np.mean(perm):.1f} ± {np.std(perm):.1f}, "
          f"fold={fold:.3f}, p_depletion={p_depletion:.4f}")

# === GC-matched permutation ===
print("\n  Running GC-MATCHED permutations...")
bs_gc_target = bs_gc_mean
gc_tolerance = 0.02  # ±2%

perm_gc_counts_all = []
perm_gc_counts_4mC = []
perm_gc_counts_6mA = []
attempts = 0
matched = 0

while matched < N_PERMUTATIONS and attempts < N_PERMUTATIONS * 100:
    rand_regions = generate_random_regions(bs_lengths, GENOME_SIZE)
    rand_seqs = extract_sequences(rand_regions, genome_seq)
    rand_gc = np.mean([gc_content(s) for s in rand_seqs])
    attempts += 1

    if abs(rand_gc - bs_gc_target) <= gc_tolerance:
        c4 = count_methyl_fast(rand_regions, arr_4mC)
        c6 = count_methyl_fast(rand_regions, arr_6mA)
        ca = count_methyl_fast(rand_regions, arr_all)
        perm_gc_counts_4mC.append(c4)
        perm_gc_counts_6mA.append(c6)
        perm_gc_counts_all.append(ca)
        matched += 1
        if matched % 2000 == 0:
            print(f"    {matched}/{N_PERMUTATIONS} (attempts: {attempts})")

if matched < N_PERMUTATIONS:
    print(f"  WARNING: Only {matched}/{N_PERMUTATIONS} GC-matched permutations achieved ({attempts} attempts)")

perm_gc_counts_all = np.array(perm_gc_counts_all)
perm_gc_counts_4mC = np.array(perm_gc_counts_4mC)
perm_gc_counts_6mA = np.array(perm_gc_counts_6mA)

for label, obs, perm in [("4mC", obs_4mC, perm_gc_counts_4mC),
                          ("6mA", obs_6mA, perm_gc_counts_6mA),
                          ("Any", obs_all, perm_gc_counts_all)]:
    p_gc = np.mean(perm <= obs) if len(perm) > 0 else float('nan')
    fold_gc = obs / np.mean(perm) if len(perm) > 0 and np.mean(perm) > 0 else float('nan')
    print(f"  GC-matched {label}: obs={obs}, perm_mean={np.mean(perm):.1f} ± {np.std(perm):.1f}, "
          f"fold={fold_gc:.3f}, p_depletion={p_gc:.4f}")

# === Summary table ===
print("\n" + "=" * 70)
print("SUMMARY: All methods comparison")
print("=" * 70)

summary_data = {
    'Method': [
        'M0: Uniform (original)',
        'M1: GC-adjusted',
        'M2: Motif-adjusted',
        'M3a: Permutation (random)',
        'M3b: Permutation (GC-matched)',
    ],
    'Obs_4mC': [obs_4mC] * 5,
    'Exp_4mC': [
        exp_4mC_uniform, exp_4mC_gc, exp_4mC_motif,
        np.mean(perm_counts_4mC), np.mean(perm_gc_counts_4mC),
    ],
    'Fold_4mC': [
        obs_4mC / exp_4mC_uniform if exp_4mC_uniform > 0 else np.nan,
        obs_4mC / exp_4mC_gc if exp_4mC_gc > 0 else np.nan,
        obs_4mC / exp_4mC_motif if exp_4mC_motif > 0 else np.nan,
        obs_4mC / np.mean(perm_counts_4mC) if np.mean(perm_counts_4mC) > 0 else np.nan,
        obs_4mC / np.mean(perm_gc_counts_4mC) if len(perm_gc_counts_4mC) > 0 and np.mean(perm_gc_counts_4mC) > 0 else np.nan,
    ],
    'Obs_6mA': [obs_6mA] * 5,
    'Exp_6mA': [
        exp_6mA_uniform, exp_6mA_gc, exp_6mA_motif,
        np.mean(perm_counts_6mA), np.mean(perm_gc_counts_6mA),
    ],
    'Fold_6mA': [
        obs_6mA / exp_6mA_uniform if exp_6mA_uniform > 0 else np.nan,
        obs_6mA / exp_6mA_gc if exp_6mA_gc > 0 else np.nan,
        obs_6mA / exp_6mA_motif if exp_6mA_motif > 0 else np.nan,
        obs_6mA / np.mean(perm_counts_6mA) if np.mean(perm_counts_6mA) > 0 else np.nan,
        obs_6mA / np.mean(perm_gc_counts_6mA) if len(perm_gc_counts_6mA) > 0 and np.mean(perm_gc_counts_6mA) > 0 else np.nan,
    ],
    'Obs_any': [obs_all] * 5,
    'Exp_any': [
        exp_all_uniform, exp_all_gc, exp_all_motif,
        np.mean(perm_counts_all), np.mean(perm_gc_counts_all),
    ],
    'Fold_any': [
        obs_all / exp_all_uniform if exp_all_uniform > 0 else np.nan,
        obs_all / exp_all_gc if exp_all_gc > 0 else np.nan,
        obs_all / exp_all_motif if exp_all_motif > 0 else np.nan,
        obs_all / np.mean(perm_counts_all) if np.mean(perm_counts_all) > 0 else np.nan,
        obs_all / np.mean(perm_gc_counts_all) if len(perm_gc_counts_all) > 0 and np.mean(perm_gc_counts_all) > 0 else np.nan,
    ],
    'P_depletion_any': [
        stats.poisson.cdf(obs_all, exp_all_uniform),
        stats.poisson.cdf(obs_all, exp_all_gc),
        stats.poisson.cdf(obs_all, exp_all_motif),
        np.mean(perm_counts_all <= obs_all),
        np.mean(perm_gc_counts_all <= obs_all) if len(perm_gc_counts_all) > 0 else np.nan,
    ]
}

df_summary = pd.DataFrame(summary_data)
for col in df_summary.columns:
    if 'Fold' in col or 'Exp' in col:
        df_summary[col] = df_summary[col].round(3)
print(df_summary.to_string(index=False))

df_summary.to_csv(OUT_DIR / "tables" / "BS_methylation_composition_correction.csv", index=False)

# === Figure: Permutation distribution ===
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Panel A: Random permutation
ax = axes[0]
ax.hist(perm_counts_all, bins=50, color='gray', alpha=0.7, edgecolor='black', linewidth=0.5,
        label=f'Random permutation\n(n={N_PERMUTATIONS:,})')
ax.axvline(obs_all, color='red', linewidth=2, linestyle='--',
           label=f'Observed = {obs_all}')
ax.axvline(np.mean(perm_counts_all), color='blue', linewidth=1.5, linestyle=':',
           label=f'Perm mean = {np.mean(perm_counts_all):.0f}')
ax.set_xlabel('Methylation sites in regions')
ax.set_ylabel('Count')
fold_r = obs_all / np.mean(perm_counts_all) if np.mean(perm_counts_all) > 0 else 0
p_r = np.mean(perm_counts_all <= obs_all)
ax.set_title(f'A. Random permutation\nfold={fold_r:.3f}, p={p_r:.4f}')
ax.legend(fontsize=9)

# Panel B: GC-matched permutation
ax = axes[1]
if len(perm_gc_counts_all) > 0:
    ax.hist(perm_gc_counts_all, bins=50, color='#2ca02c', alpha=0.7, edgecolor='black', linewidth=0.5,
            label=f'GC-matched permutation\n(n={len(perm_gc_counts_all):,}, GC={bs_gc_target:.3f}±{gc_tolerance})')
    ax.axvline(obs_all, color='red', linewidth=2, linestyle='--',
               label=f'Observed = {obs_all}')
    ax.axvline(np.mean(perm_gc_counts_all), color='blue', linewidth=1.5, linestyle=':',
               label=f'Perm mean = {np.mean(perm_gc_counts_all):.0f}')
    fold_g = obs_all / np.mean(perm_gc_counts_all) if np.mean(perm_gc_counts_all) > 0 else 0
    p_g = np.mean(perm_gc_counts_all <= obs_all)
    ax.set_title(f'B. GC-matched permutation\nfold={fold_g:.3f}, p={p_g:.4f}')
else:
    ax.text(0.5, 0.5, 'Insufficient GC-matched samples', ha='center', va='center',
            transform=ax.transAxes, fontsize=12, color='gray')
    ax.set_title('B. GC-matched permutation')
ax.set_xlabel('Methylation sites in regions')
ax.set_ylabel('Count')
ax.legend(fontsize=9)

plt.suptitle(f'TF Binding Site methylation depletion — composition correction\n'
             f'Tier 1 BS (n={len(bs_df)}), {bs_total_width:,} bp total, direct overlap',
             fontsize=12, fontweight='bold')
plt.tight_layout()

for ext in ['pdf', 'png', 'svg']:
    fig.savefig(OUT_DIR / "figures" / f"BS_composition_correction_permutation.{ext}",
                dpi=200, bbox_inches='tight')
plt.close()

print(f"\nFigure saved to {OUT_DIR / 'figures'}/")
print(f"Table saved to {OUT_DIR / 'tables'}/")
print("\n=== Composition correction analysis complete ===")
