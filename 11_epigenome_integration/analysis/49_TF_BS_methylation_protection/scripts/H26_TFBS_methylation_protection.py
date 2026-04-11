#!/usr/bin/env python3
"""
H26: TF Binding Site-Level Methylation Protection
==================================================
Tests whether predicted TF binding sites show methylation depletion,
and whether the effect is stronger at TF BS near regulatory gene promoters.

Builds on H25 (TSS protection zone) and H22 (protein occupancy model).
"""

import pandas as pd
import numpy as np
import re
import os
import sys
from collections import defaultdict
from scipy import stats
from scipy.stats import fisher_exact, binomtest, mannwhitneyu
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FuncFormatter
import warnings
warnings.filterwarnings('ignore')

# ─── Paths ───────────────────────────────────────────────────────────────────
BASE = "/Users/okaban/bioinfo/rna-seq"
OUT = f"{BASE}/11_epigenome_integration/analysis/49_TF_BS_methylation_protection"
FIG = f"{OUT}/figures"
TAB = f"{OUT}/tables"

FIMO_GFF = f"{BASE}/13_TF_binding-site/analysis/01_master_TF_list_260206_v1/intermediate/fimo_results/fimo.gff"
GCCGGC_FILE = f"{BASE}/11_epigenome_integration/analysis/37_defense_island_GCCGGC/tables/GCCGGC_sites_by_timepoint.tsv"
ALL_METH_FILE = f"{BASE}/11_epigenome_integration/analysis/01_integration/high_confidence_sites_weighted.csv"
GENE_ANN_FILE = f"{BASE}/05_annotation/analysis/05_annotation_260128_v1/tables/gene_annotation_basic.tsv"
REG_GENES_FILE = f"{BASE}/11_epigenome_integration/analysis/29_genomewide_TF_screen/tables/all_regulatory_genes.tsv"
AAGCCCG_MAP_FILE = f"{BASE}/11_epigenome_integration/analysis/36_AAGCCCG_distribution/tables/AAGCCCG_site_gene_mapping.tsv"

# Chromosome parameters
CHROM_LEN = 8_667_507
ARM_LEFT_MAX = 1_500_000
ARM_RIGHT_MIN = 7_167_508
CHROM_ID = "NC_003888.3"

np.random.seed(42)

# ─── 1. Parse FIMO GFF ──────────────────────────────────────────────────────
print("=" * 70)
print("STEP 1: Parsing FIMO GFF")
print("=" * 70)

fimo_records = []
with open(FIMO_GFF) as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9:
            continue
        seqid, source, ftype, start, end, score, strand, phase, attrs = parts
        if seqid != CHROM_ID:
            continue

        # Parse attributes
        attr_dict = {}
        for kv in attrs.split(';'):
            if '=' in kv:
                k, v = kv.split('=', 1)
                attr_dict[k] = v

        # Extract TF name (before _NC_003888.3)
        name_raw = attr_dict.get('Name', '')
        tf_name = re.sub(r'_NC_003888\.3[+-]?$', '', name_raw)

        pval = float(attr_dict.get('pvalue', 1))
        qval = float(attr_dict.get('qvalue', 1))

        fimo_records.append({
            'chrom': seqid,
            'start': int(start),
            'end': int(end),
            'strand': strand,
            'tf_name': tf_name,
            'score': float(score),
            'pvalue': pval,
            'qvalue': qval,
            'sequence': attr_dict.get('sequence', '')
        })

fimo_df = pd.DataFrame(fimo_records)
print(f"Total FIMO hits on {CHROM_ID}: {len(fimo_df)}")

# Calculate BS center and width
fimo_df['center'] = ((fimo_df['start'] + fimo_df['end']) / 2).astype(int)
fimo_df['width'] = fimo_df['end'] - fimo_df['start'] + 1

# Region classification
fimo_df['region'] = fimo_df['center'].apply(
    lambda x: 'arm' if x <= ARM_LEFT_MAX or x >= ARM_RIGHT_MIN else 'core')

print(f"\nTF family distribution:")
tf_counts = fimo_df['tf_name'].value_counts()
for tf, cnt in tf_counts.items():
    print(f"  {tf}: {cnt}")

print(f"\nMedian BS width: {fimo_df['width'].median():.0f} bp")
print(f"Mean BS width: {fimo_df['width'].mean():.1f} bp")
print(f"Core: {(fimo_df['region']=='core').sum()}, Arm: {(fimo_df['region']=='arm').sum()}")

# p-value distribution
print(f"\np-value distribution:")
for thresh in [1e-4, 1e-3, 0.01, 0.1]:
    n = (fimo_df['pvalue'] < thresh).sum()
    print(f"  p < {thresh}: {n} ({100*n/len(fimo_df):.1f}%)")

# Filter for high-confidence (use all since these are already FIMO-filtered)
# Check if p-value < 1e-4 gives enough sites
n_hc = (fimo_df['pvalue'] < 1e-4).sum()
print(f"\nHigh-confidence (p < 1e-4): {n_hc}")
# Use all sites as primary analysis, high-confidence as sensitivity check
fimo_all = fimo_df.copy()

# ─── 2. Load Methylation Sites ──────────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 2: Loading methylation sites")
print("=" * 70)

# GCCGGC sites
gccggc_df = pd.read_csv(GCCGGC_FILE, sep='\t')
gccggc_df = gccggc_df[gccggc_df['chrom'] == CHROM_ID].copy()
gccggc_t1 = gccggc_df[gccggc_df['timepoint'] == 'T1']['position'].unique()
gccggc_t2 = gccggc_df[gccggc_df['timepoint'] == 'T2']['position'].unique()
print(f"GCCGGC T1: {len(gccggc_t1)} sites")
print(f"GCCGGC T2: {len(gccggc_t2)} sites")

# All methylation data
all_meth = pd.read_csv(ALL_METH_FILE)
all_meth = all_meth[all_meth['chrom'] == CHROM_ID].copy()

# 4mC T1
all_4mc_t1 = all_meth[(all_meth['mod_type'] == '4mC') & (all_meth['timepoint'] == 'T1')]['position'].unique()
print(f"All 4mC T1: {len(all_4mc_t1)} sites")

# 6mA T1
all_6ma_t1 = all_meth[(all_meth['mod_type'] == '6mA') & (all_meth['timepoint'] == 'T1')]['position'].unique()
print(f"All 6mA T1: {len(all_6ma_t1)} sites")

# AAGCCCG sites from gene mapping file (get unique positions per timepoint)
aagcccg_map = pd.read_csv(AAGCCCG_MAP_FILE, sep='\t')
aagcccg_t1_positions = aagcccg_map[aagcccg_map['timepoint'] == 'T1']['position'].unique()
aagcccg_t2_positions = aagcccg_map[aagcccg_map['timepoint'] == 'T2']['position'].unique()
print(f"AAGCCCG T1: {len(aagcccg_t1_positions)} sites")
print(f"AAGCCCG T2: {len(aagcccg_t2_positions)} sites")

# All T1 methylation
all_t1 = all_meth[all_meth['timepoint'] == 'T1']['position'].unique()
print(f"All T1 methylation: {len(all_t1)} sites")

# Create sorted arrays for efficient searching
gccggc_t1_sorted = np.sort(gccggc_t1)
all_4mc_t1_sorted = np.sort(all_4mc_t1)
all_6ma_t1_sorted = np.sort(all_6ma_t1)
aagcccg_t1_sorted = np.sort(aagcccg_t1_positions)
all_t1_sorted = np.sort(all_t1)

# Methylation datasets for analysis
meth_datasets = {
    'GCCGGC_T1': gccggc_t1_sorted,
    'AAGCCCG_T1': aagcccg_t1_sorted,
    'All_4mC_T1': all_4mc_t1_sorted,
    'All_6mA_T1': all_6ma_t1_sorted,
    'All_T1': all_t1_sorted,
}

# ─── 3. Load Gene Annotations & Regulatory Genes ────────────────────────────
print("\n" + "=" * 70)
print("STEP 3: Loading gene annotations")
print("=" * 70)

genes = pd.read_csv(GENE_ANN_FILE, sep='\t')
genes = genes[genes['contig'] == CHROM_ID].copy()
print(f"Total genes on {CHROM_ID}: {len(genes)}")

reg_genes = pd.read_csv(REG_GENES_FILE, sep='\t')
reg_locus_tags = set(reg_genes['locus_tag'].unique())
print(f"Regulatory genes: {len(reg_locus_tags)}")

# Classify genes
genes['is_regulatory'] = genes['gene_id'].isin(reg_locus_tags)
genes['tss'] = genes.apply(lambda r: r['start'] if r['strand'] == '+' else r['end'], axis=1)

print(f"Regulatory genes in annotation: {genes['is_regulatory'].sum()}")
print(f"Non-regulatory genes: {(~genes['is_regulatory']).sum()}")


# ─── 4. Classify TF BS by genomic context ───────────────────────────────────
print("\n" + "=" * 70)
print("STEP 4: Classifying TF BS by genomic context")
print("=" * 70)

# Build arrays for efficient vectorized classification
gene_starts = genes['start'].values
gene_ends = genes['end'].values
gene_strands = genes['strand'].values
gene_ids = genes['gene_id'].values
gene_is_reg = genes['is_regulatory'].values

# Compute TSS and promoter regions
gene_tss = np.where(gene_strands == '+', gene_starts, gene_ends)
prom_starts = np.where(gene_strands == '+', gene_tss - 500, gene_tss - 100)
prom_ends = np.where(gene_strands == '+', gene_tss + 100, gene_tss + 500)

# Sort genes by start for binary search
gene_sort_idx = np.argsort(gene_starts)

# Classify each TFBS using vectorized operations
tfbs_centers = fimo_all['center'].values
n_tfbs = len(tfbs_centers)

contexts = np.full(n_tfbs, 'intergenic', dtype=object)
associated_genes = np.full(n_tfbs, None, dtype=object)
near_reg_promoter = np.zeros(n_tfbs, dtype=bool)
near_nonreg_promoter = np.zeros(n_tfbs, dtype=bool)

print(f"Classifying {n_tfbs} TF BS against {len(gene_ids)} genes...")

# Vectorized: broadcast TFBS centers against promoter regions
# Process in chunks to avoid memory explosion
CHUNK = 5000
for chunk_start in range(0, n_tfbs, CHUNK):
    chunk_end = min(chunk_start + CHUNK, n_tfbs)
    centers_chunk = tfbs_centers[chunk_start:chunk_end]

    # Check promoter overlap (N_tfbs_chunk x N_genes)
    in_prom = ((centers_chunk[:, None] >= prom_starts[None, :]) &
               (centers_chunk[:, None] <= prom_ends[None, :]))

    for i_local in range(len(centers_chunk)):
        i_global = chunk_start + i_local
        prom_hits = np.where(in_prom[i_local])[0]
        if len(prom_hits) > 0:
            gi = prom_hits[0]
            contexts[i_global] = 'promoter'
            associated_genes[i_global] = gene_ids[gi]
            if gene_is_reg[gi]:
                near_reg_promoter[i_global] = True
            else:
                near_nonreg_promoter[i_global] = True
            continue

        # Check gene body
        c = centers_chunk[i_local]
        body_mask = (gene_starts <= c) & (gene_ends >= c)
        if body_mask.any():
            gi = np.where(body_mask)[0][0]
            contexts[i_global] = 'intragenic'
            associated_genes[i_global] = gene_ids[gi]

    if (chunk_start + CHUNK) % 10000 == 0:
        print(f"  Processed {chunk_end}/{n_tfbs} TF BS...")

fimo_all['context'] = contexts
fimo_all['associated_gene'] = associated_genes
fimo_all['in_reg_promoter'] = near_reg_promoter
fimo_all['in_nonreg_promoter'] = near_nonreg_promoter

print("\nGenomic context distribution:")
ctx_counts = fimo_all['context'].value_counts()
for ctx, cnt in ctx_counts.items():
    print(f"  {ctx}: {cnt} ({100*cnt/len(fimo_all):.1f}%)")

print(f"\nPromoter TF BS near regulatory genes: {sum(near_reg_promoter)}")
print(f"Promoter TF BS near non-regulatory genes: {sum(near_nonreg_promoter)}")


# ─── 5. Overlap Analysis ────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 5: Overlap analysis (TF BS vs methylation)")
print("=" * 70)

def count_overlaps(tfbs_df, meth_positions, extension=0):
    """Count TF BS that overlap with at least one methylation site."""
    n_overlap = 0
    for _, row in tfbs_df.iterrows():
        bs_start = row['start'] - extension
        bs_end = row['end'] + extension
        # Binary search for methylation sites in range
        left = np.searchsorted(meth_positions, bs_start, side='left')
        right = np.searchsorted(meth_positions, bs_end, side='right')
        if right > left:
            n_overlap += 1
    return n_overlap

# Calculate overlap for each methylation dataset
overlap_results = []
total_bs = len(fimo_all)
avg_width = fimo_all['width'].mean()

print(f"\nTotal TF BS: {total_bs}")
print(f"Average BS width: {avg_width:.1f} bp")

for ext in [0, 10, 25, 50]:
    print(f"\n--- Extension: ±{ext} bp ---")
    effective_width = avg_width + 2 * ext

    for meth_name, meth_pos in meth_datasets.items():
        n_meth = len(meth_pos)
        n_overlap = count_overlaps(fimo_all, meth_pos, extension=ext)
        obs_frac = n_overlap / total_bs

        # Expected by chance:
        # P(at least 1 meth site in window) = 1 - (1 - effective_width/CHROM_LEN)^n_meth
        # For small probabilities, ~= n_meth * effective_width / CHROM_LEN
        p_single = effective_width / CHROM_LEN
        exp_frac = 1 - (1 - p_single) ** n_meth
        exp_overlap = exp_frac * total_bs

        fold = obs_frac / exp_frac if exp_frac > 0 else np.nan

        # Binomial test
        pval = binomtest(n_overlap, total_bs, exp_frac, alternative='two-sided').pvalue

        print(f"  {meth_name}: {n_overlap}/{total_bs} overlap ({100*obs_frac:.2f}%), "
              f"expected {exp_overlap:.1f} ({100*exp_frac:.2f}%), "
              f"fold={fold:.3f}, p={pval:.2e}")

        overlap_results.append({
            'methylation_type': meth_name,
            'extension_bp': ext,
            'n_TFBS_total': total_bs,
            'n_TFBS_overlap': n_overlap,
            'obs_fraction': obs_frac,
            'expected_overlap': exp_overlap,
            'expected_fraction': exp_frac,
            'fold_enrichment': fold,
            'binomial_pvalue': pval,
            'n_meth_sites': n_meth,
            'effective_width': effective_width
        })

overlap_df = pd.DataFrame(overlap_results)

# ─── 6. Spatial Profile Around TF BS Center ─────────────────────────────────
print("\n" + "=" * 70)
print("STEP 6: Spatial methylation density profile around TF BS center")
print("=" * 70)

WINDOW = 2000
BIN_SIZE = 50
bins = np.arange(-WINDOW, WINDOW + BIN_SIZE, BIN_SIZE)
bin_centers = (bins[:-1] + bins[1:]) / 2

def compute_density_profile(centers, meth_positions, window=WINDOW, bin_size=BIN_SIZE):
    """Compute methylation density profile around a set of centers."""
    bins_local = np.arange(-window, window + bin_size, bin_size)
    counts = np.zeros(len(bins_local) - 1)
    n_valid = 0

    for c in centers:
        # Find methylation sites within window
        left = np.searchsorted(meth_positions, c - window, side='left')
        right = np.searchsorted(meth_positions, c + window, side='right')

        if c - window >= 1 and c + window <= CHROM_LEN:
            n_valid += 1
            if right > left:
                dists = meth_positions[left:right] - c
                hist, _ = np.histogram(dists, bins=bins_local)
                counts += hist

    # Density: sites per kb per center
    density = counts / (n_valid * bin_size / 1000) if n_valid > 0 else counts
    return density, n_valid

# Observed profiles for key methylation types
profile_results = {}
key_meth = ['GCCGGC_T1', 'All_4mC_T1', 'All_6mA_T1', 'All_T1']

tfbs_centers = fimo_all['center'].values

for mname in key_meth:
    mpos = meth_datasets[mname]
    density, n_valid = compute_density_profile(tfbs_centers, mpos)
    profile_results[mname] = {'density': density, 'n_valid': n_valid}
    print(f"{mname}: computed profile over {n_valid} valid TF BS")

# Random control: generate random positions and compute profiles
print("\nComputing random control profiles (1000 permutations)...")
N_PERM = 1000
random_profiles = {mname: [] for mname in key_meth}

for i in range(N_PERM):
    if (i + 1) % 100 == 0:
        print(f"  Permutation {i+1}/{N_PERM}")
    rand_centers = np.random.randint(WINDOW + 1, CHROM_LEN - WINDOW, size=len(tfbs_centers))

    for mname in key_meth:
        mpos = meth_datasets[mname]
        density, _ = compute_density_profile(rand_centers, mpos)
        random_profiles[mname].append(density)

# Calculate CI
for mname in key_meth:
    rand_arr = np.array(random_profiles[mname])
    profile_results[mname]['rand_mean'] = rand_arr.mean(axis=0)
    profile_results[mname]['rand_ci_lo'] = np.percentile(rand_arr, 2.5, axis=0)
    profile_results[mname]['rand_ci_hi'] = np.percentile(rand_arr, 97.5, axis=0)

# Save spatial profile data
profile_data_rows = []
for mname in key_meth:
    pr = profile_results[mname]
    for i, bc in enumerate(bin_centers):
        profile_data_rows.append({
            'methylation_type': mname,
            'bin_center_bp': int(bc),
            'observed_density': pr['density'][i],
            'random_mean': pr['rand_mean'][i],
            'random_ci_lo': pr['rand_ci_lo'][i],
            'random_ci_hi': pr['rand_ci_hi'][i],
            'obs_over_random': pr['density'][i] / pr['rand_mean'][i] if pr['rand_mean'][i] > 0 else np.nan,
            'n_valid_centers': pr['n_valid']
        })

profile_data_df = pd.DataFrame(profile_data_rows)
profile_data_df.to_csv(f"{TAB}/spatial_profile_data.tsv", sep='\t', index=False)
print("Saved spatial_profile_data.tsv")


# ─── 7. Regulatory vs Non-regulatory Promoter TF BS ─────────────────────────
print("\n" + "=" * 70)
print("STEP 7: Regulatory vs non-regulatory promoter TF BS")
print("=" * 70)

reg_promoter_bs = fimo_all[fimo_all['in_reg_promoter']].copy()
nonreg_promoter_bs = fimo_all[fimo_all['in_nonreg_promoter']].copy()
intergenic_bs = fimo_all[fimo_all['context'] == 'intergenic'].copy()
intragenic_bs = fimo_all[fimo_all['context'] == 'intragenic'].copy()

print(f"Regulatory promoter TF BS: {len(reg_promoter_bs)}")
print(f"Non-regulatory promoter TF BS: {len(nonreg_promoter_bs)}")
print(f"Intergenic TF BS: {len(intergenic_bs)}")
print(f"Intragenic TF BS: {len(intragenic_bs)}")

# Compute profiles for regulatory vs non-regulatory promoter TF BS
reg_nonreg_profiles = {}
for label, subset in [('reg_promoter', reg_promoter_bs),
                       ('nonreg_promoter', nonreg_promoter_bs),
                       ('intergenic', intergenic_bs),
                       ('intragenic', intragenic_bs),
                       ('all_TFBS', fimo_all)]:
    if len(subset) < 10:
        print(f"  {label}: too few sites ({len(subset)}), skipping profile")
        continue
    centers = subset['center'].values
    for mname in ['All_4mC_T1', 'All_T1']:
        mpos = meth_datasets[mname]
        density, n_valid = compute_density_profile(centers, mpos)
        key = f"{label}_{mname}"
        reg_nonreg_profiles[key] = {'density': density, 'n_valid': n_valid, 'n_bs': len(subset)}
        print(f"  {key}: {n_valid} valid centers")

# Overlap rates by context
context_results = []
for ctx_label, ctx_df in [('promoter_reg', reg_promoter_bs),
                           ('promoter_nonreg', nonreg_promoter_bs),
                           ('intergenic', intergenic_bs),
                           ('intragenic', intragenic_bs),
                           ('all', fimo_all)]:
    if len(ctx_df) == 0:
        continue
    for mname, mpos in meth_datasets.items():
        n_ov = count_overlaps(ctx_df, mpos, extension=10)
        n_total = len(ctx_df)
        frac = n_ov / n_total if n_total > 0 else 0

        # Expected
        eff_w = ctx_df['width'].mean() + 20  # ±10 extension
        p_s = eff_w / CHROM_LEN
        exp_f = 1 - (1 - p_s) ** len(mpos)
        fold = frac / exp_f if exp_f > 0 else np.nan

        pval = binomtest(n_ov, n_total, exp_f, alternative='two-sided').pvalue if n_total > 0 else 1

        context_results.append({
            'context': ctx_label,
            'methylation_type': mname,
            'n_TFBS': n_total,
            'n_overlap': n_ov,
            'overlap_fraction': frac,
            'expected_fraction': exp_f,
            'fold': fold,
            'pvalue': pval
        })

context_df = pd.DataFrame(context_results)
context_df.to_csv(f"{TAB}/context_analysis.tsv", sep='\t', index=False)
print("\nSaved context_analysis.tsv")

# Print key results
print("\n--- Overlap rates by context (±10 bp, All_4mC_T1) ---")
sub = context_df[context_df['methylation_type'] == 'All_4mC_T1']
for _, row in sub.iterrows():
    print(f"  {row['context']}: {row['n_overlap']}/{row['n_TFBS']} "
          f"({100*row['overlap_fraction']:.2f}%), fold={row['fold']:.3f}, p={row['pvalue']:.2e}")


# ─── 8. TF Family-Specific Analysis ─────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 8: TF family-specific methylation depletion")
print("=" * 70)

tf_family_results = []
for tf_name in fimo_all['tf_name'].unique():
    tf_subset = fimo_all[fimo_all['tf_name'] == tf_name]
    n_tf = len(tf_subset)

    for mname in ['All_4mC_T1', 'All_T1']:
        mpos = meth_datasets[mname]
        n_ov = count_overlaps(tf_subset, mpos, extension=10)
        frac = n_ov / n_tf

        eff_w = tf_subset['width'].mean() + 20
        p_s = eff_w / CHROM_LEN
        exp_f = 1 - (1 - p_s) ** len(mpos)
        fold = frac / exp_f if exp_f > 0 else np.nan

        pval = binomtest(n_ov, n_tf, exp_f, alternative='two-sided').pvalue

        tf_family_results.append({
            'tf_name': tf_name,
            'methylation_type': mname,
            'n_TFBS': n_tf,
            'n_overlap': n_ov,
            'overlap_fraction': frac,
            'expected_fraction': exp_f,
            'fold': fold,
            'pvalue': pval
        })

tf_family_df = pd.DataFrame(tf_family_results)

print("\nTF family overlap (All_4mC_T1, ±10 bp):")
sub = tf_family_df[tf_family_df['methylation_type'] == 'All_4mC_T1'].sort_values('fold')
for _, row in sub.iterrows():
    sig = '*' if row['pvalue'] < 0.05 else ''
    print(f"  {row['tf_name']:12s}: {row['n_overlap']:4d}/{row['n_TFBS']:5d} "
          f"({100*row['overlap_fraction']:5.2f}%), fold={row['fold']:.3f}, p={row['pvalue']:.2e} {sig}")


# ─── 9. Core/Arm Stratification ─────────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 9: Core vs Arm stratification")
print("=" * 70)

region_results = []
for region in ['core', 'arm']:
    region_bs = fimo_all[fimo_all['region'] == region]
    n_bs = len(region_bs)

    for mname in ['GCCGGC_T1', 'All_4mC_T1', 'All_T1']:
        mpos = meth_datasets[mname]
        # Filter methylation to same region
        if region == 'core':
            region_meth = mpos[(mpos > ARM_LEFT_MAX) & (mpos < ARM_RIGHT_MIN)]
            region_len = ARM_RIGHT_MIN - ARM_LEFT_MAX
        else:
            region_meth = mpos[(mpos <= ARM_LEFT_MAX) | (mpos >= ARM_RIGHT_MIN)]
            region_len = ARM_LEFT_MAX + (CHROM_LEN - ARM_RIGHT_MIN)

        n_ov = count_overlaps(region_bs, region_meth, extension=10)
        frac = n_ov / n_bs if n_bs > 0 else 0

        eff_w = region_bs['width'].mean() + 20
        p_s = eff_w / region_len
        exp_f = 1 - (1 - p_s) ** len(region_meth)
        fold = frac / exp_f if exp_f > 0 else np.nan

        pval = binomtest(n_ov, n_bs, exp_f, alternative='two-sided').pvalue if n_bs > 0 else 1.0

        region_results.append({
            'region': region,
            'methylation_type': mname,
            'n_TFBS': n_bs,
            'n_overlap': n_ov,
            'overlap_fraction': frac,
            'expected_fraction': exp_f,
            'fold': fold,
            'pvalue': pval,
            'n_region_meth': len(region_meth)
        })

region_df = pd.DataFrame(region_results)

print("\nCore/Arm stratified overlap:")
for _, row in region_df.iterrows():
    sig = '*' if row['pvalue'] < 0.05 else ''
    print(f"  {row['region']:4s} {row['methylation_type']:12s}: "
          f"{row['n_overlap']}/{row['n_TFBS']} ({100*row['overlap_fraction']:.2f}%), "
          f"fold={row['fold']:.3f}, p={row['pvalue']:.2e} {sig}")


# ─── 10. Per-TFBS Overlap Table ─────────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 10: Generating per-TFBS overlap table")
print("=" * 70)

# For each TFBS, record overlap with each methylation type
per_bs_records = []
for idx, row in fimo_all.iterrows():
    rec = {
        'tfbs_id': idx,
        'tf_name': row['tf_name'],
        'start': row['start'],
        'end': row['end'],
        'center': row['center'],
        'width': row['width'],
        'strand': row['strand'],
        'pvalue': row['pvalue'],
        'context': row['context'],
        'in_reg_promoter': row['in_reg_promoter'],
        'in_nonreg_promoter': row['in_nonreg_promoter'],
        'region': row['region'],
    }

    # Check overlap with ±10 bp extension
    for mname, mpos in meth_datasets.items():
        left = np.searchsorted(mpos, row['start'] - 10, side='left')
        right = np.searchsorted(mpos, row['end'] + 10, side='right')
        n_sites = right - left
        rec[f'{mname}_overlap'] = n_sites > 0
        rec[f'{mname}_count'] = n_sites

        # Nearest distance
        if len(mpos) > 0:
            nearest_idx = np.searchsorted(mpos, row['center'])
            candidates = []
            if nearest_idx > 0:
                candidates.append(abs(mpos[nearest_idx - 1] - row['center']))
            if nearest_idx < len(mpos):
                candidates.append(abs(mpos[nearest_idx] - row['center']))
            rec[f'{mname}_nearest_dist'] = min(candidates) if candidates else np.nan
        else:
            rec[f'{mname}_nearest_dist'] = np.nan

    per_bs_records.append(rec)

per_bs_df = pd.DataFrame(per_bs_records)
per_bs_df.to_csv(f"{TAB}/TFBS_methylation_overlap.tsv", sep='\t', index=False)
print(f"Saved TFBS_methylation_overlap.tsv ({len(per_bs_df)} rows)")


# ─── 11. Statistical Tests Summary ──────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 11: Compiling statistical tests")
print("=" * 70)

stat_tests = []

# Test 1: Overall depletion at TFBS (binomial tests from step 5)
for _, row in overlap_df[overlap_df['extension_bp'] == 10].iterrows():
    stat_tests.append({
        'test_id': f"T1_overlap_{row['methylation_type']}",
        'description': f"TFBS overlap with {row['methylation_type']} (±10bp)",
        'test_type': 'binomial',
        'observed': row['n_TFBS_overlap'],
        'expected': row['expected_overlap'],
        'n': row['n_TFBS_total'],
        'fold': row['fold_enrichment'],
        'pvalue': row['binomial_pvalue'],
        'significant': row['binomial_pvalue'] < 0.05,
        'direction': 'depleted' if row['fold_enrichment'] < 1 else 'enriched'
    })

# Test 2: Spatial depletion (central bin vs flanks)
for mname in key_meth:
    pr = profile_results[mname]
    obs = pr['density']
    rand_mean = pr['rand_mean']

    # Central region (-200 to +200)
    central_mask = (bin_centers >= -200) & (bin_centers <= 200)
    flank_mask = ((bin_centers >= -2000) & (bin_centers <= -1000)) | ((bin_centers >= 1000) & (bin_centers <= 2000))

    if central_mask.sum() > 0 and flank_mask.sum() > 0:
        central_obs = obs[central_mask].mean()
        flank_obs = obs[flank_mask].mean()
        central_rand = rand_mean[central_mask].mean()

        # Ratio of central to flank
        ratio = central_obs / flank_obs if flank_obs > 0 else np.nan
        depletion_pct = (1 - central_obs / central_rand) * 100 if central_rand > 0 else 0

        stat_tests.append({
            'test_id': f"T2_spatial_{mname}",
            'description': f"Central (-200 to +200) vs flank density for {mname}",
            'test_type': 'ratio',
            'observed': central_obs,
            'expected': central_rand,
            'n': pr['n_valid'],
            'fold': central_obs / central_rand if central_rand > 0 else np.nan,
            'pvalue': np.nan,  # Computed from permutations below
            'significant': np.nan,
            'direction': 'depleted' if central_obs < central_rand else 'enriched'
        })

# Test 3: Permutation p-value for central depletion
for mname in key_meth:
    pr = profile_results[mname]
    obs = pr['density']
    central_mask = (bin_centers >= -200) & (bin_centers <= 200)

    obs_central_mean = obs[central_mask].mean()
    rand_central_means = [rp[central_mask].mean() for rp in random_profiles[mname]]
    perm_p = np.mean([r <= obs_central_mean for r in rand_central_means])

    stat_tests.append({
        'test_id': f"T3_perm_{mname}",
        'description': f"Permutation test: central depletion for {mname}",
        'test_type': 'permutation',
        'observed': obs_central_mean,
        'expected': np.mean(rand_central_means),
        'n': N_PERM,
        'fold': obs_central_mean / np.mean(rand_central_means) if np.mean(rand_central_means) > 0 else np.nan,
        'pvalue': perm_p,
        'significant': perm_p < 0.05,
        'direction': 'depleted' if obs_central_mean < np.mean(rand_central_means) else 'enriched'
    })

# Test 4: Regulatory vs non-regulatory promoter comparison
for mname in ['All_4mC_T1', 'All_T1']:
    reg_key = f"reg_promoter_{mname}"
    nonreg_key = f"nonreg_promoter_{mname}"

    if reg_key in reg_nonreg_profiles and nonreg_key in reg_nonreg_profiles:
        reg_d = reg_nonreg_profiles[reg_key]['density']
        nonreg_d = reg_nonreg_profiles[nonreg_key]['density']

        central_mask = (bin_centers >= -200) & (bin_centers <= 200)
        reg_central = reg_d[central_mask].mean()
        nonreg_central = nonreg_d[central_mask].mean()

        stat_tests.append({
            'test_id': f"T4_reg_vs_nonreg_{mname}",
            'description': f"Regulatory vs non-reg promoter TFBS central density ({mname})",
            'test_type': 'comparison',
            'observed': reg_central,
            'expected': nonreg_central,
            'n': reg_nonreg_profiles[reg_key]['n_bs'],
            'fold': reg_central / nonreg_central if nonreg_central > 0 else np.nan,
            'pvalue': np.nan,
            'significant': np.nan,
            'direction': 'reg < nonreg' if reg_central < nonreg_central else 'reg > nonreg'
        })

# Test 5: Context comparison (Fisher exact: promoter vs intragenic overlap)
for mname in ['All_4mC_T1']:
    prom_sub = context_df[(context_df['context'].isin(['promoter_reg', 'promoter_nonreg'])) &
                           (context_df['methylation_type'] == mname)]
    intra_sub = context_df[(context_df['context'] == 'intragenic') &
                            (context_df['methylation_type'] == mname)]

    if len(prom_sub) > 0 and len(intra_sub) > 0:
        prom_ov = prom_sub['n_overlap'].sum()
        prom_total = prom_sub['n_TFBS'].sum()
        intra_ov = intra_sub['n_overlap'].sum()
        intra_total = intra_sub['n_TFBS'].sum()

        table = [[prom_ov, prom_total - prom_ov],
                 [intra_ov, intra_total - intra_ov]]

        or_val, fp = fisher_exact(table, alternative='two-sided')

        stat_tests.append({
            'test_id': f"T5_promoter_vs_intragenic_{mname}",
            'description': f"Promoter vs intragenic TFBS overlap ({mname})",
            'test_type': 'fisher_exact',
            'observed': prom_ov / prom_total if prom_total > 0 else 0,
            'expected': intra_ov / intra_total if intra_total > 0 else 0,
            'n': prom_total + intra_total,
            'fold': or_val,
            'pvalue': fp,
            'significant': fp < 0.05,
            'direction': 'promoter < intragenic' if or_val < 1 else 'promoter > intragenic'
        })

stat_df = pd.DataFrame(stat_tests)
stat_df.to_csv(f"{TAB}/statistical_tests.tsv", sep='\t', index=False)
print(f"Saved statistical_tests.tsv ({len(stat_df)} tests)")

print("\nKey statistical results:")
for _, row in stat_df.iterrows():
    sig = '***' if row['pvalue'] < 0.001 else '**' if row['pvalue'] < 0.01 else '*' if row['pvalue'] < 0.05 else 'NS'
    print(f"  {row['test_id']}: fold={row['fold']:.3f}, p={row['pvalue']:.2e}, "
          f"{row['direction']} [{sig}]")


# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATIONS
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("VISUALIZATIONS")
print("=" * 70)

plt.rcParams.update({
    'font.size': 10,
    'axes.titlesize': 11,
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 8,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

# ─── Panel A: TF BS-centered methylation density profile ────────────────────
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

for ax_idx, mname in enumerate(key_meth):
    ax = axes.flat[ax_idx]
    pr = profile_results[mname]

    ax.fill_between(bin_centers, pr['rand_ci_lo'], pr['rand_ci_hi'],
                     alpha=0.3, color='gray', label='95% CI (random)')
    ax.plot(bin_centers, pr['rand_mean'], '--', color='gray', lw=1, label='Random mean')
    ax.plot(bin_centers, pr['density'], '-', color='#e74c3c', lw=1.5, label='Observed')

    ax.axvline(0, color='black', ls=':', lw=0.5)
    ax.set_xlabel('Distance from TF BS center (bp)')
    ax.set_ylabel('Methylation density\n(sites/kb/TFBS)')
    ax.set_title(mname.replace('_', ' '))
    ax.legend(loc='upper right', framealpha=0.8)

    # Add depletion annotation
    central_mask = (bin_centers >= -200) & (bin_centers <= 200)
    obs_central = pr['density'][central_mask].mean()
    rand_central = pr['rand_mean'][central_mask].mean()
    if rand_central > 0:
        depl = (1 - obs_central / rand_central) * 100
        ax.text(0.02, 0.95, f'Central depletion: {depl:.1f}%',
                transform=ax.transAxes, fontsize=8, va='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

fig.suptitle('H26: Methylation Density Profile Around TF Binding Sites', fontsize=13, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.96])
for ext in ['pdf', 'svg']:
    fig.savefig(f"{FIG}/TFBS_methylation_profile.{ext}")
plt.close()
print("Saved TFBS_methylation_profile.pdf/svg")


# ─── Panel B: Regulatory vs Non-regulatory Promoter TF BS ───────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for ax_idx, mname in enumerate(['All_4mC_T1', 'All_T1']):
    ax = axes[ax_idx]

    # Plot all TFBS profile as reference
    all_key = f"all_TFBS_{mname}"
    if all_key in reg_nonreg_profiles:
        ax.plot(bin_centers, reg_nonreg_profiles[all_key]['density'], '-',
                color='gray', lw=1, alpha=0.5, label=f'All TFBS (n={reg_nonreg_profiles[all_key]["n_bs"]})')

    reg_key = f"reg_promoter_{mname}"
    nonreg_key = f"nonreg_promoter_{mname}"
    inter_key = f"intergenic_{mname}"
    intra_key = f"intragenic_{mname}"

    colors = {'reg_promoter': '#e74c3c', 'nonreg_promoter': '#3498db',
              'intergenic': '#2ecc71', 'intragenic': '#f39c12'}
    labels_map = {'reg_promoter': 'Reg. promoter', 'nonreg_promoter': 'Non-reg. promoter',
                  'intergenic': 'Intergenic', 'intragenic': 'Intragenic'}

    for ctx_key, ctx_label in [('reg_promoter', reg_key), ('nonreg_promoter', nonreg_key),
                                ('intergenic', inter_key), ('intragenic', intra_key)]:
        if ctx_label in reg_nonreg_profiles:
            pr = reg_nonreg_profiles[ctx_label]
            ax.plot(bin_centers, pr['density'], '-', color=colors[ctx_key], lw=1.5,
                    label=f"{labels_map[ctx_key]} (n={pr['n_bs']})")

    ax.axvline(0, color='black', ls=':', lw=0.5)
    ax.set_xlabel('Distance from TF BS center (bp)')
    ax.set_ylabel('Methylation density (sites/kb/TFBS)')
    ax.set_title(mname.replace('_', ' '))
    ax.legend(loc='upper right', framealpha=0.8, fontsize=7)

fig.suptitle('H26: Regulatory vs Non-regulatory Promoter TF BS Methylation', fontsize=13, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.95])
for ext in ['pdf', 'svg']:
    fig.savefig(f"{FIG}/regulatory_vs_nonreg_TFBS.{ext}")
plt.close()
print("Saved regulatory_vs_nonreg_TFBS.pdf/svg")


# ─── Panel C: Overlap Rates by Genomic Context ──────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for ax_idx, mname in enumerate(['All_4mC_T1', 'All_T1']):
    ax = axes[ax_idx]
    sub = context_df[context_df['methylation_type'] == mname].copy()

    contexts_order = ['promoter_reg', 'promoter_nonreg', 'intergenic', 'intragenic', 'all']
    ctx_labels = ['Reg.\npromoter', 'Non-reg.\npromoter', 'Intergenic', 'Intragenic', 'All\nTFBS']

    folds = []
    pvals = []
    ns = []
    for ctx in contexts_order:
        row = sub[sub['context'] == ctx]
        if len(row) > 0:
            folds.append(row['fold'].values[0])
            pvals.append(row['pvalue'].values[0])
            ns.append(row['n_TFBS'].values[0])
        else:
            folds.append(np.nan)
            pvals.append(1)
            ns.append(0)

    colors_bar = ['#e74c3c', '#3498db', '#2ecc71', '#f39c12', '#95a5a6']
    bars = ax.bar(range(len(contexts_order)), folds, color=colors_bar, alpha=0.8, edgecolor='black', lw=0.5)

    ax.axhline(1.0, color='black', ls='--', lw=1, label='Expected (random)')

    # Add significance stars
    for i, (f, p, n) in enumerate(zip(folds, pvals, ns)):
        if not np.isnan(f):
            sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'NS'
            ax.text(i, f + 0.02, sig, ha='center', va='bottom', fontsize=8, fontweight='bold')
            ax.text(i, f - 0.05, f'n={n}', ha='center', va='top', fontsize=7, color='white')

    ax.set_xticks(range(len(contexts_order)))
    ax.set_xticklabels(ctx_labels, fontsize=8)
    ax.set_ylabel('Fold enrichment\n(observed/expected)')
    ax.set_title(mname.replace('_', ' '))
    ax.set_ylim(0, max([f for f in folds if not np.isnan(f)]) * 1.3 if any(not np.isnan(f) for f in folds) else 2)

fig.suptitle('H26: Methylation Overlap at TF BS by Genomic Context', fontsize=13, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.95])
for ext in ['pdf', 'svg']:
    fig.savefig(f"{FIG}/overlap_by_context.{ext}")
plt.close()
print("Saved overlap_by_context.pdf/svg")


# ─── Panel D: Comprehensive Summary ─────────────────────────────────────────
fig = plt.figure(figsize=(16, 14))
gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.45, wspace=0.35)

# D1: Spatial profile (All_4mC_T1) — main result
ax1 = fig.add_subplot(gs[0, 0:2])
pr = profile_results['All_4mC_T1']
ax1.fill_between(bin_centers, pr['rand_ci_lo'], pr['rand_ci_hi'],
                  alpha=0.3, color='gray', label='95% CI (random)')
ax1.plot(bin_centers, pr['rand_mean'], '--', color='gray', lw=1, label='Random mean')
ax1.plot(bin_centers, pr['density'], '-', color='#e74c3c', lw=2, label='Observed')
ax1.axvline(0, color='black', ls=':', lw=0.5)
ax1.set_xlabel('Distance from TF BS center (bp)')
ax1.set_ylabel('4mC density (sites/kb/TFBS)')
ax1.set_title('A. Spatial Methylation Profile Around TF Binding Sites (All 4mC T1)', fontweight='bold')
ax1.legend(loc='upper right', framealpha=0.8, fontsize=7)
central_mask = (bin_centers >= -200) & (bin_centers <= 200)
obs_c = pr['density'][central_mask].mean()
rand_c = pr['rand_mean'][central_mask].mean()
if rand_c > 0:
    depl = (1 - obs_c / rand_c) * 100
    ax1.text(0.02, 0.95, f'Central depletion: {depl:.1f}%',
             transform=ax1.transAxes, fontsize=9, va='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

# D2: TF family comparison
ax2 = fig.add_subplot(gs[0, 2])
tf_sub = tf_family_df[tf_family_df['methylation_type'] == 'All_4mC_T1'].sort_values('fold')
colors_tf = ['#e74c3c' if p < 0.05 else '#bdc3c7' for p in tf_sub['pvalue']]
ax2.barh(range(len(tf_sub)), tf_sub['fold'], color=colors_tf, edgecolor='black', lw=0.5)
ax2.set_yticks(range(len(tf_sub)))
ax2.set_yticklabels(tf_sub['tf_name'], fontsize=7)
ax2.axvline(1.0, color='black', ls='--', lw=1)
ax2.set_xlabel('Fold enrichment')
ax2.set_title('B. TF Family-Specific\nMethylation at BS', fontweight='bold')

# D3: Regulatory vs non-regulatory profile
ax3 = fig.add_subplot(gs[1, 0:2])
mname = 'All_4mC_T1'
for ctx_key, color, label in [('reg_promoter', '#e74c3c', 'Regulatory promoter'),
                               ('nonreg_promoter', '#3498db', 'Non-reg. promoter'),
                               ('intergenic', '#2ecc71', 'Intergenic'),
                               ('intragenic', '#f39c12', 'Intragenic')]:
    k = f"{ctx_key}_{mname}"
    if k in reg_nonreg_profiles:
        pr = reg_nonreg_profiles[k]
        ax3.plot(bin_centers, pr['density'], '-', color=color, lw=1.5,
                 label=f"{label} (n={pr['n_bs']})")

ax3.axvline(0, color='black', ls=':', lw=0.5)
ax3.set_xlabel('Distance from TF BS center (bp)')
ax3.set_ylabel('4mC density (sites/kb/TFBS)')
ax3.set_title('C. Context-Specific Methylation Profiles (All 4mC T1)', fontweight='bold')
ax3.legend(loc='upper right', framealpha=0.8, fontsize=7)

# D4: Core/Arm comparison
ax4 = fig.add_subplot(gs[1, 2])
reg_sub = region_df[region_df['methylation_type'] == 'All_4mC_T1']
x_pos = [0, 1]
for i, (_, row) in enumerate(reg_sub.iterrows()):
    color = '#2980b9' if row['region'] == 'core' else '#e67e22'
    ax4.bar(i, row['fold'], color=color, edgecolor='black', lw=0.5, alpha=0.8)
    sig = '***' if row['pvalue'] < 0.001 else '**' if row['pvalue'] < 0.01 else '*' if row['pvalue'] < 0.05 else 'NS'
    ax4.text(i, row['fold'] + 0.02, sig, ha='center', va='bottom', fontsize=10, fontweight='bold')
ax4.set_xticks(x_pos)
ax4.set_xticklabels(['Core', 'Arm'])
ax4.axhline(1.0, color='black', ls='--', lw=1)
ax4.set_ylabel('Fold enrichment')
ax4.set_title('D. Core vs Arm\n(All 4mC T1)', fontweight='bold')

# D5: Overlap by context bar chart
ax5 = fig.add_subplot(gs[2, 0])
sub = context_df[context_df['methylation_type'] == 'All_4mC_T1']
contexts_order = ['promoter_reg', 'promoter_nonreg', 'intergenic', 'intragenic']
ctx_labels_short = ['Reg.\nprom.', 'Non-reg.\nprom.', 'Inter-\ngenic', 'Intra-\ngenic']
colors_ctx = ['#e74c3c', '#3498db', '#2ecc71', '#f39c12']
for i, ctx in enumerate(contexts_order):
    row = sub[sub['context'] == ctx]
    if len(row) > 0:
        ax5.bar(i, row['fold'].values[0], color=colors_ctx[i], edgecolor='black', lw=0.5, alpha=0.8)
        p = row['pvalue'].values[0]
        sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'NS'
        ax5.text(i, row['fold'].values[0] + 0.02, sig, ha='center', fontsize=8, fontweight='bold')
ax5.set_xticks(range(len(contexts_order)))
ax5.set_xticklabels(ctx_labels_short, fontsize=7)
ax5.axhline(1.0, color='black', ls='--', lw=1)
ax5.set_ylabel('Fold enrichment')
ax5.set_title('E. Overlap by Context\n(All 4mC T1)', fontweight='bold')

# D6: H25 vs H26 comparison summary (text panel)
ax6 = fig.add_subplot(gs[2, 1:3])
ax6.axis('off')

# Compile summary text
summary_lines = [
    "F. Comparison with H25 TSS-Level Protection",
    "",
    "H25 (TSS-level):                    H26 (TF BS-level):",
    f"  Protection zone: 2,200 bp           TF BS width: {fimo_all['width'].mean():.0f} bp (median)",
    f"  Max depletion: 17.3%                Central depletion:",
]

for mname in ['All_4mC_T1', 'All_T1']:
    pr = profile_results[mname]
    central_mask_inner = (bin_centers >= -200) & (bin_centers <= 200)
    obs_central = pr['density'][central_mask_inner].mean()
    rand_central = pr['rand_mean'][central_mask_inner].mean()
    depl = (1 - obs_central / rand_central) * 100 if rand_central > 0 else 0
    summary_lines.append(f"                                        {mname}: {depl:.1f}%")

# Add overlap summary
ov_all_4mc = overlap_df[(overlap_df['methylation_type'] == 'All_4mC_T1') & (overlap_df['extension_bp'] == 10)]
if len(ov_all_4mc) > 0:
    r = ov_all_4mc.iloc[0]
    summary_lines.append(f"")
    summary_lines.append(f"  Deepest at +300 bp                  Overall overlap fold: {r['fold_enrichment']:.3f}")
    summary_lines.append(f"  Reg. genes: 4.6x stronger           p = {r['binomial_pvalue']:.2e}")

summary_lines.append("")
summary_lines.append("Interpretation:")
summary_lines.append("  H25 protection zone and H26 TF BS-level protection provide")
summary_lines.append("  complementary evidence for the protein occupancy model.")

text = '\n'.join(summary_lines)
ax6.text(0.02, 0.95, text, transform=ax6.transAxes, fontsize=8,
         fontfamily='monospace', va='top',
         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

fig.suptitle('H26: TF Binding Site-Level Methylation Protection — Comprehensive Summary',
             fontsize=14, fontweight='bold')
for ext in ['pdf', 'svg']:
    fig.savefig(f"{FIG}/H26_comprehensive_summary.{ext}")
plt.close()
print("Saved H26_comprehensive_summary.pdf/svg")


# ═══════════════════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("FINAL SUMMARY")
print("=" * 70)

print(f"\nTotal TF BS analyzed: {len(fimo_all)} (from 16 TF families)")
print(f"Median BS width: {fimo_all['width'].median():.0f} bp")

print("\n--- Overlap with methylation (±10 bp extension) ---")
for _, row in overlap_df[overlap_df['extension_bp'] == 10].iterrows():
    sig = '***' if row['binomial_pvalue'] < 0.001 else '**' if row['binomial_pvalue'] < 0.01 else '*' if row['binomial_pvalue'] < 0.05 else 'NS'
    print(f"  {row['methylation_type']:15s}: fold={row['fold_enrichment']:.3f}, p={row['binomial_pvalue']:.2e} [{sig}]")

print("\n--- Spatial central depletion (±200 bp) ---")
for mname in key_meth:
    pr = profile_results[mname]
    central_mask = (bin_centers >= -200) & (bin_centers <= 200)
    obs_c = pr['density'][central_mask].mean()
    rand_c = pr['rand_mean'][central_mask].mean()
    depl = (1 - obs_c / rand_c) * 100 if rand_c > 0 else 0

    rand_central_means = [rp[central_mask].mean() for rp in random_profiles[mname]]
    perm_p = np.mean([r <= obs_c for r in rand_central_means])
    print(f"  {mname:15s}: {depl:.1f}% depletion, permutation p={perm_p:.3f}")

print("\n--- Context comparison (All_4mC_T1) ---")
sub = context_df[(context_df['methylation_type'] == 'All_4mC_T1')]
for _, row in sub.iterrows():
    sig = '*' if row['pvalue'] < 0.05 else 'NS'
    print(f"  {row['context']:20s}: fold={row['fold']:.3f}, p={row['pvalue']:.2e} [{sig}]")

print("\n--- Core/Arm (All_4mC_T1) ---")
for _, row in region_df[region_df['methylation_type'] == 'All_4mC_T1'].iterrows():
    sig = '*' if row['pvalue'] < 0.05 else 'NS'
    print(f"  {row['region']:6s}: fold={row['fold']:.3f}, p={row['pvalue']:.2e} [{sig}]")

print("\nDone. All outputs saved to:")
print(f"  Figures: {FIG}/")
print(f"  Tables:  {TAB}/")
