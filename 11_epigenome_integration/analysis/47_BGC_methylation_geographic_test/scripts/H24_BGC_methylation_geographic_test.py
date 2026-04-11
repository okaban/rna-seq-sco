#!/usr/bin/env python3
"""
H24: BGC Methylation Enrichment -- Geographic Confound Test
============================================================
H23 found BGC genes are enriched near GCCGGC 4mC methylation sites
(fold=1.66, p=4.2e-07). However, BGC genes are 100% in the core genome
and GCCGGC T1 sites are 83% core. This co-location could create a
geographic confound (Simpson's paradox, cf. H19).

Key question: Does the enrichment survive core-only stratification?
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import fisher_exact, mannwhitneyu
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.collections import PatchCollection
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Configuration
# ============================================================
BASE = '/Users/okaban/bioinfo/rna-seq'
ANALYSIS_DIR = f'{BASE}/11_epigenome_integration/analysis/47_BGC_methylation_geographic_test'
FIG_DIR = f'{ANALYSIS_DIR}/figures'
TBL_DIR = f'{ANALYSIS_DIR}/tables'

CHROM_SIZE = 8_667_507
ARM_LEFT_END = 1_500_000
ARM_RIGHT_START = 7_167_508
PROXIMITY_KB = 2  # 2 kb proximity window

# ============================================================
# 1. Load all data
# ============================================================
print("=" * 70)
print("H24: BGC Methylation Enrichment -- Geographic Confound Test")
print("=" * 70)

# Gene annotations (basic)
gene_annot = pd.read_csv(
    f'{BASE}/05_annotation/analysis/05_annotation_260128_v1/tables/gene_annotation_basic.tsv',
    sep='\t'
)
print(f"Gene annotations: {len(gene_annot)} genes")

# BGC master table
bgc_master = pd.read_csv(
    f'{BASE}/05_annotation/analysis/05_annotation_260128_v1/tables/gene_master_with_BGC.tsv',
    sep='\t'
)
bgc_info = bgc_master[bgc_master['bgc_name'].notna()][['gene_id', 'bgc_name', 'bgc_role']].copy()
bgc_gene_ids = set(bgc_info['gene_id'].values)
print(f"BGC genes: {len(bgc_gene_ids)}")

# GCCGGC sites (T1 timepoint)
gccggc_all = pd.read_csv(
    f'{BASE}/11_epigenome_integration/analysis/37_defense_island_GCCGGC/tables/GCCGGC_sites_by_timepoint.tsv',
    sep='\t'
)
gccggc_t1 = gccggc_all[gccggc_all['timepoint'] == 'T1'].copy()
gccggc_t1_positions = np.sort(gccggc_t1['position'].unique())
print(f"GCCGGC 4mC T1 sites: {len(gccggc_t1_positions)}")

# AAGCCCG sites (T1)
aagcccg_map = pd.read_csv(
    f'{BASE}/11_epigenome_integration/analysis/36_AAGCCCG_distribution/tables/AAGCCCG_site_gene_mapping.tsv',
    sep='\t'
)
aagcccg_t1_positions = np.sort(aagcccg_map[aagcccg_map['timepoint'] == 'T1']['position'].unique())
print(f"AAGCCCG 6mA T1 sites: {len(aagcccg_t1_positions)}")

# All 4mC sites (T1)
all_4mc = pd.read_csv(
    f'{BASE}/11_epigenome_integration/analysis/23_expanded_motif_search/4mC_final_census.csv'
)
all_4mc_t1 = all_4mc[all_4mc['timepoint'] == 'T1']
all_4mc_t1_positions = np.sort(all_4mc_t1['position'].unique())
print(f"All 4mC T1 sites: {len(all_4mc_t1_positions)}")

# Motif counts from H22
motif_counts = pd.read_csv(
    f'{BASE}/11_epigenome_integration/analysis/45_sequence_level_motif_depletion/tables/gene_motif_counts.tsv',
    sep='\t'
)
print(f"Motif count data: {len(motif_counts)} genes")

# Regulatory genes
reg_genes = pd.read_csv(
    f'{BASE}/11_epigenome_integration/analysis/29_genomewide_TF_screen/tables/all_regulatory_genes.tsv',
    sep='\t'
)
reg_locus_tags = set(reg_genes['locus_tag'].values)
print(f"Regulatory genes: {len(reg_locus_tags)}")

# ============================================================
# 2. Annotate genes with BGC status and region
# ============================================================
print("\n" + "=" * 70)
print("Step 1: BGC Gene Geography")
print("=" * 70)

gene_annot['midpoint'] = (gene_annot['start'] + gene_annot['end']) / 2
gene_annot['region'] = gene_annot['midpoint'].apply(
    lambda x: 'arm' if x <= ARM_LEFT_END or x >= ARM_RIGHT_START else 'core'
)
gene_annot['is_BGC'] = gene_annot['gene_id'].isin(bgc_gene_ids)

# Merge BGC names
gene_annot = gene_annot.merge(bgc_info, on='gene_id', how='left')

# BGC geography
bgc_df = gene_annot[gene_annot['is_BGC']].copy()
non_bgc_df = gene_annot[~gene_annot['is_BGC']].copy()

print(f"\nBGC gene geography:")
bgc_region_counts = bgc_df['region'].value_counts()
for reg, cnt in bgc_region_counts.items():
    print(f"  {reg}: {cnt} ({cnt/len(bgc_df)*100:.1f}%)")

print(f"\nNon-BGC gene geography:")
non_bgc_region_counts = non_bgc_df['region'].value_counts()
for reg, cnt in non_bgc_region_counts.items():
    print(f"  {reg}: {cnt} ({cnt/len(non_bgc_df)*100:.1f}%)")

# Save BGC gene geography table
bgc_geo = gene_annot[gene_annot['is_BGC']][['gene_id', 'start', 'end', 'strand',
                                              'product', 'bgc_name', 'bgc_role',
                                              'region', 'midpoint']].copy()
bgc_geo = bgc_geo.sort_values('start')
bgc_geo.to_csv(f'{TBL_DIR}/BGC_gene_geography.tsv', sep='\t', index=False)
print(f"\nSaved: BGC_gene_geography.tsv ({len(bgc_geo)} genes)")

# Per-cluster summary
print(f"\nPer-cluster geography:")
for cluster in sorted(bgc_df['bgc_name'].unique()):
    cdf = bgc_df[bgc_df['bgc_name'] == cluster]
    core_n = (cdf['region'] == 'core').sum()
    arm_n = (cdf['region'] == 'arm').sum()
    print(f"  {cluster}: {len(cdf)} genes, core={core_n}, arm={arm_n}, "
          f"range={cdf['start'].min():,}-{cdf['end'].max():,}")

# ============================================================
# 3. Proximity calculation
# ============================================================
print("\n" + "=" * 70)
print("Step 2: Methylation Proximity Analysis")
print("=" * 70)


def calculate_proximity(gene_df, site_positions, window_bp=2000):
    """For each gene, check if any methylation site is within window_bp of gene boundaries."""
    sites = np.array(sorted(site_positions))
    proximal = np.zeros(len(gene_df), dtype=bool)
    for i, (_, row) in enumerate(gene_df.iterrows()):
        g_start = row['start']
        g_end = row['end']
        left = np.searchsorted(sites, g_start - window_bp)
        right = np.searchsorted(sites, g_end + window_bp, side='right')
        if left < right:
            proximal[i] = True
    return proximal


def count_proximal_sites(gene_df, site_positions, window_bp=2000):
    """Count number of methylation sites within window_bp of each gene."""
    sites = np.array(sorted(site_positions))
    counts = np.zeros(len(gene_df), dtype=int)
    for i, (_, row) in enumerate(gene_df.iterrows()):
        g_start = row['start']
        g_end = row['end']
        left = np.searchsorted(sites, g_start - window_bp)
        right = np.searchsorted(sites, g_end + window_bp, side='right')
        counts[i] = right - left
    return counts


def calculate_density_per_kb(gene_df, site_positions, window_bp=2000):
    """Calculate methylation site density (sites per kb) in each gene's extended region."""
    sites = np.array(sorted(site_positions))
    densities = np.zeros(len(gene_df), dtype=float)
    for i, (_, row) in enumerate(gene_df.iterrows()):
        g_start = row['start']
        g_end = row['end']
        region_start = g_start - window_bp
        region_end = g_end + window_bp
        region_len_kb = (region_end - region_start) / 1000.0
        left = np.searchsorted(sites, region_start)
        right = np.searchsorted(sites, region_end, side='right')
        n_sites = right - left
        densities[i] = n_sites / region_len_kb if region_len_kb > 0 else 0
    return densities


window = PROXIMITY_KB * 1000

# Calculate proximity for all genes
gene_annot['prox_GCCGGC'] = calculate_proximity(gene_annot, gccggc_t1_positions, window)
gene_annot['prox_AAGCCCG'] = calculate_proximity(gene_annot, aagcccg_t1_positions, window)
gene_annot['prox_All4mC'] = calculate_proximity(gene_annot, all_4mc_t1_positions, window)

# Calculate site counts and density
gene_annot['count_GCCGGC'] = count_proximal_sites(gene_annot, gccggc_t1_positions, window)
gene_annot['count_AAGCCCG'] = count_proximal_sites(gene_annot, aagcccg_t1_positions, window)
gene_annot['count_All4mC'] = count_proximal_sites(gene_annot, all_4mc_t1_positions, window)
gene_annot['density_GCCGGC'] = calculate_density_per_kb(gene_annot, gccggc_t1_positions, window)
gene_annot['density_AAGCCCG'] = calculate_density_per_kb(gene_annot, aagcccg_t1_positions, window)
gene_annot['density_All4mC'] = calculate_density_per_kb(gene_annot, all_4mc_t1_positions, window)

for methyl_type in ['GCCGGC', 'AAGCCCG', 'All4mC']:
    col = f'prox_{methyl_type}'
    n_prox = gene_annot[col].sum()
    print(f"  {methyl_type}: {n_prox}/{len(gene_annot)} genes proximal "
          f"({n_prox/len(gene_annot)*100:.1f}%)")

# ============================================================
# 4. Replicate H23's enrichment (ALL genes)
# ============================================================
print("\n" + "=" * 70)
print("Step 3: Replicate H23 Enrichment (All Genes)")
print("=" * 70)

results = []


def fisher_enrichment(group_proximal, group_total, bg_proximal, bg_total, label, motif, subset):
    """Perform Fisher exact test for enrichment."""
    group_not_prox = group_total - group_proximal
    bg_not_prox = bg_total - bg_proximal
    table = [[group_proximal, group_not_prox],
             [bg_proximal, bg_not_prox]]
    odds_ratio, pval = fisher_exact(table, alternative='two-sided')
    # Calculate fold enrichment
    if bg_total > 0 and bg_proximal > 0:
        group_rate = group_proximal / group_total
        bg_rate = bg_proximal / bg_total
        fold = group_rate / bg_rate if bg_rate > 0 else np.inf
    else:
        fold = np.nan

    result = {
        'label': label,
        'motif': motif,
        'subset': subset,
        'group_proximal': group_proximal,
        'group_total': group_total,
        'group_rate': group_proximal / group_total if group_total > 0 else 0,
        'bg_proximal': bg_proximal,
        'bg_total': bg_total,
        'bg_rate': bg_proximal / bg_total if bg_total > 0 else 0,
        'fold_enrichment': fold,
        'odds_ratio': odds_ratio,
        'pvalue': pval,
    }
    return result


# ALL genes: BGC vs non-BGC
for motif in ['GCCGGC', 'AAGCCCG', 'All4mC']:
    col = f'prox_{motif}'
    bgc_prox = gene_annot.loc[gene_annot['is_BGC'], col].sum()
    bgc_tot = gene_annot['is_BGC'].sum()
    non_bgc_prox = gene_annot.loc[~gene_annot['is_BGC'], col].sum()
    non_bgc_tot = (~gene_annot['is_BGC']).sum()

    res = fisher_enrichment(bgc_prox, bgc_tot, non_bgc_prox, non_bgc_tot,
                            'BGC_vs_nonBGC', motif, 'all_genes')
    results.append(res)
    print(f"\n  {motif} - All genes:")
    print(f"    BGC: {bgc_prox}/{bgc_tot} ({res['group_rate']*100:.1f}%)")
    print(f"    non-BGC: {non_bgc_prox}/{non_bgc_tot} ({res['bg_rate']*100:.1f}%)")
    print(f"    Fold={res['fold_enrichment']:.3f}, OR={res['odds_ratio']:.3f}, p={res['pvalue']:.2e}")

# ============================================================
# 5. Core-only geographic stratification (THE KEY TEST)
# ============================================================
print("\n" + "=" * 70)
print("Step 4: Core-Only Stratification (Key Test)")
print("=" * 70)

core_genes = gene_annot[gene_annot['region'] == 'core'].copy()
print(f"Core genes: {len(core_genes)}")
print(f"  BGC core: {core_genes['is_BGC'].sum()}")
print(f"  non-BGC core: {(~core_genes['is_BGC']).sum()}")

for motif in ['GCCGGC', 'AAGCCCG', 'All4mC']:
    col = f'prox_{motif}'
    bgc_prox = core_genes.loc[core_genes['is_BGC'], col].sum()
    bgc_tot = core_genes['is_BGC'].sum()
    non_bgc_prox = core_genes.loc[~core_genes['is_BGC'], col].sum()
    non_bgc_tot = (~core_genes['is_BGC']).sum()

    res = fisher_enrichment(bgc_prox, bgc_tot, non_bgc_prox, non_bgc_tot,
                            'BGC_vs_nonBGC', motif, 'core_only')
    results.append(res)
    print(f"\n  {motif} - Core only:")
    print(f"    BGC core: {bgc_prox}/{bgc_tot} ({res['group_rate']*100:.1f}%)")
    print(f"    non-BGC core: {non_bgc_prox}/{non_bgc_tot} ({res['bg_rate']*100:.1f}%)")
    print(f"    Fold={res['fold_enrichment']:.3f}, OR={res['odds_ratio']:.3f}, p={res['pvalue']:.2e}")

# Also test arm genes (for completeness)
arm_genes = gene_annot[gene_annot['region'] == 'arm'].copy()
print(f"\nArm genes: {len(arm_genes)}")
print(f"  BGC arm: {arm_genes['is_BGC'].sum()}")
if arm_genes['is_BGC'].sum() > 0:
    for motif in ['GCCGGC', 'AAGCCCG', 'All4mC']:
        col = f'prox_{motif}'
        bgc_prox = arm_genes.loc[arm_genes['is_BGC'], col].sum()
        bgc_tot = arm_genes['is_BGC'].sum()
        non_bgc_prox = arm_genes.loc[~arm_genes['is_BGC'], col].sum()
        non_bgc_tot = (~arm_genes['is_BGC']).sum()
        res = fisher_enrichment(bgc_prox, bgc_tot, non_bgc_prox, non_bgc_tot,
                                'BGC_vs_nonBGC', motif, 'arm_only')
        results.append(res)
else:
    print("  (No BGC genes in arms -- test skipped)")

# Bonferroni correction
n_tests = len(results)
for r in results:
    r['pvalue_bonferroni'] = min(r['pvalue'] * n_tests, 1.0)
    r['significant_bonf'] = r['pvalue_bonferroni'] < 0.05

# Summary comparison
print("\n" + "-" * 50)
print("COMPARISON: All genes vs Core-only")
print("-" * 50)
for motif in ['GCCGGC', 'AAGCCCG', 'All4mC']:
    all_res = [r for r in results if r['motif'] == motif and r['subset'] == 'all_genes'][0]
    core_res = [r for r in results if r['motif'] == motif and r['subset'] == 'core_only'][0]
    print(f"\n  {motif}:")
    print(f"    All genes:  fold={all_res['fold_enrichment']:.3f}, "
          f"p={all_res['pvalue']:.2e}, p_bonf={all_res['pvalue_bonferroni']:.2e}")
    print(f"    Core only:  fold={core_res['fold_enrichment']:.3f}, "
          f"p={core_res['pvalue']:.2e}, p_bonf={core_res['pvalue_bonferroni']:.2e}")
    fold_change = core_res['fold_enrichment'] / all_res['fold_enrichment'] if all_res['fold_enrichment'] > 0 else np.nan
    print(f"    Fold change ratio (core/all): {fold_change:.3f}")

# Save enrichment results
results_df = pd.DataFrame(results)
results_df.to_csv(f'{TBL_DIR}/enrichment_tests.tsv', sep='\t', index=False)
print(f"\nSaved: enrichment_tests.tsv ({len(results_df)} tests)")

# ============================================================
# 6. BGC Cluster-Specific Analysis
# ============================================================
print("\n" + "=" * 70)
print("Step 5: BGC Cluster-Specific GCCGGC Density")
print("=" * 70)

cluster_results = []
for cluster in sorted(bgc_df['bgc_name'].unique()):
    cdf = gene_annot[gene_annot['bgc_name'] == cluster].copy()
    cluster_start = cdf['start'].min()
    cluster_end = cdf['end'].max()
    cluster_len_kb = (cluster_end - cluster_start) / 1000.0

    # Count GCCGGC T1 sites within cluster boundaries
    gccggc_in_cluster = np.sum(
        (gccggc_t1_positions >= cluster_start) & (gccggc_t1_positions <= cluster_end)
    )
    gccggc_density = gccggc_in_cluster / cluster_len_kb if cluster_len_kb > 0 else 0

    # AAGCCCG
    aagcccg_in_cluster = np.sum(
        (aagcccg_t1_positions >= cluster_start) & (aagcccg_t1_positions <= cluster_end)
    )
    aagcccg_density = aagcccg_in_cluster / cluster_len_kb if cluster_len_kb > 0 else 0

    # All 4mC
    all4mc_in_cluster = np.sum(
        (all_4mc_t1_positions >= cluster_start) & (all_4mc_t1_positions <= cluster_end)
    )
    all4mc_density = all4mc_in_cluster / cluster_len_kb if cluster_len_kb > 0 else 0

    cluster_results.append({
        'cluster': cluster,
        'n_genes': len(cdf),
        'start': cluster_start,
        'end': cluster_end,
        'length_kb': cluster_len_kb,
        'region': cdf['region'].iloc[0],  # All BGC genes should be core
        'GCCGGC_sites': gccggc_in_cluster,
        'GCCGGC_density_per_kb': gccggc_density,
        'AAGCCCG_sites': aagcccg_in_cluster,
        'AAGCCCG_density_per_kb': aagcccg_density,
        'All4mC_sites': all4mc_in_cluster,
        'All4mC_density_per_kb': all4mc_density,
    })

    print(f"\n  {cluster.upper()} ({len(cdf)} genes, {cluster_len_kb:.1f} kb, {cluster_start:,}-{cluster_end:,}):")
    print(f"    GCCGGC: {gccggc_in_cluster} sites, {gccggc_density:.3f}/kb")
    print(f"    AAGCCCG: {aagcccg_in_cluster} sites, {aagcccg_density:.3f}/kb")
    print(f"    All 4mC: {all4mc_in_cluster} sites, {all4mc_density:.3f}/kb")

# Genome-wide and core-genome averages
genome_gccggc_density = len(gccggc_t1_positions) / (CHROM_SIZE / 1000)
core_length_bp = ARM_RIGHT_START - ARM_LEFT_END - 1
gccggc_core_count = np.sum(
    (gccggc_t1_positions > ARM_LEFT_END) & (gccggc_t1_positions < ARM_RIGHT_START)
)
core_gccggc_density = gccggc_core_count / (core_length_bp / 1000)

genome_aagcccg_density = len(aagcccg_t1_positions) / (CHROM_SIZE / 1000)
aagcccg_core_count = np.sum(
    (aagcccg_t1_positions > ARM_LEFT_END) & (aagcccg_t1_positions < ARM_RIGHT_START)
)
core_aagcccg_density = aagcccg_core_count / (core_length_bp / 1000)

genome_all4mc_density = len(all_4mc_t1_positions) / (CHROM_SIZE / 1000)
all4mc_core_count = np.sum(
    (all_4mc_t1_positions > ARM_LEFT_END) & (all_4mc_t1_positions < ARM_RIGHT_START)
)
core_all4mc_density = all4mc_core_count / (core_length_bp / 1000)

print(f"\n  Genome-wide averages:")
print(f"    GCCGGC: {genome_gccggc_density:.4f}/kb")
print(f"    AAGCCCG: {genome_aagcccg_density:.4f}/kb")
print(f"    All 4mC: {genome_all4mc_density:.4f}/kb")
print(f"  Core-genome averages:")
print(f"    GCCGGC: {core_gccggc_density:.4f}/kb (n={gccggc_core_count})")
print(f"    AAGCCCG: {core_aagcccg_density:.4f}/kb (n={aagcccg_core_count})")
print(f"    All 4mC: {core_all4mc_density:.4f}/kb (n={all4mc_core_count})")

# Fold vs core average for each cluster
for cr in cluster_results:
    cr['GCCGGC_fold_vs_core'] = cr['GCCGGC_density_per_kb'] / core_gccggc_density if core_gccggc_density > 0 else np.nan
    cr['AAGCCCG_fold_vs_core'] = cr['AAGCCCG_density_per_kb'] / core_aagcccg_density if core_aagcccg_density > 0 else np.nan
    cr['All4mC_fold_vs_core'] = cr['All4mC_density_per_kb'] / core_all4mc_density if core_all4mc_density > 0 else np.nan

print("\n  Per-cluster fold vs core average (GCCGGC):")
for cr in cluster_results:
    print(f"    {cr['cluster'].upper()}: {cr['GCCGGC_fold_vs_core']:.2f}x")

cluster_df = pd.DataFrame(cluster_results)
cluster_df.to_csv(f'{TBL_DIR}/cluster_specific_methylation.tsv', sep='\t', index=False)
print(f"\nSaved: cluster_specific_methylation.tsv ({len(cluster_df)} clusters)")

# ============================================================
# 7. DNA Sequence Motif Control (H22 extension)
# ============================================================
print("\n" + "=" * 70)
print("Step 6: DNA Sequence Motif Control")
print("=" * 70)

# Merge motif counts with BGC status (motif_counts already has 'region', so only merge BGC columns)
motif_merged = motif_counts.merge(
    gene_annot[['gene_id', 'is_BGC', 'bgc_name']],
    on='gene_id', how='inner'
)
print(f"Motif data merged: {len(motif_merged)} genes")

# Compare DNA motif (TGGCCGGC) density in BGC vs non-BGC core genes
core_motif = motif_merged[motif_merged['region'] == 'core'].copy()
bgc_core_motif = core_motif[core_motif['is_BGC']]
non_bgc_core_motif = core_motif[~core_motif['is_BGC']]

dna_vs_methyl_results = []

print("\nDNA motif (TGGCCGGC) density comparison (core genes only):")
for density_col, label in [('extended_density_TGGCCGGC', 'TGGCCGGC_extended'),
                            ('body_density_TGGCCGGC', 'TGGCCGGC_body'),
                            ('promoter_density_TGGCCGGC', 'TGGCCGGC_promoter')]:
    bgc_vals = bgc_core_motif[density_col].values
    non_bgc_vals = non_bgc_core_motif[density_col].values

    bgc_mean = np.mean(bgc_vals)
    non_bgc_mean = np.mean(non_bgc_vals)
    fold = bgc_mean / non_bgc_mean if non_bgc_mean > 0 else np.nan

    stat, pval = mannwhitneyu(bgc_vals, non_bgc_vals, alternative='two-sided')
    print(f"  {label}: BGC={bgc_mean:.4f}, non-BGC={non_bgc_mean:.4f}, "
          f"fold={fold:.3f}, MWU p={pval:.2e}")

    dna_vs_methyl_results.append({
        'comparison': f'DNA_{label}',
        'BGC_core_mean': bgc_mean,
        'BGC_core_median': np.median(bgc_vals),
        'nonBGC_core_mean': non_bgc_mean,
        'nonBGC_core_median': np.median(non_bgc_vals),
        'fold': fold,
        'MWU_stat': stat,
        'MWU_pvalue': pval,
        'n_BGC': len(bgc_vals),
        'n_nonBGC': len(non_bgc_vals),
    })

# Same for AAGCCCG DNA motif
print("\nDNA motif (AAGCCCG) density comparison (core genes only):")
for density_col, label in [('extended_density_AAGCCCG', 'AAGCCCG_extended'),
                            ('body_density_AAGCCCG', 'AAGCCCG_body'),
                            ('promoter_density_AAGCCCG', 'AAGCCCG_promoter')]:
    bgc_vals = bgc_core_motif[density_col].values
    non_bgc_vals = non_bgc_core_motif[density_col].values

    bgc_mean = np.mean(bgc_vals)
    non_bgc_mean = np.mean(non_bgc_vals)
    fold = bgc_mean / non_bgc_mean if non_bgc_mean > 0 else np.nan

    stat, pval = mannwhitneyu(bgc_vals, non_bgc_vals, alternative='two-sided')
    print(f"  {label}: BGC={bgc_mean:.4f}, non-BGC={non_bgc_mean:.4f}, "
          f"fold={fold:.3f}, MWU p={pval:.2e}")

    dna_vs_methyl_results.append({
        'comparison': f'DNA_{label}',
        'BGC_core_mean': bgc_mean,
        'BGC_core_median': np.median(bgc_vals),
        'nonBGC_core_mean': non_bgc_mean,
        'nonBGC_core_median': np.median(non_bgc_vals),
        'fold': fold,
        'MWU_stat': stat,
        'MWU_pvalue': pval,
        'n_BGC': len(bgc_vals),
        'n_nonBGC': len(non_bgc_vals),
    })

# Now compare methylation density (sites per kb within 2kb extended region) -- core only
print("\nMethylation site density comparison (core genes only, 2kb extended):")
for methyl_type in ['GCCGGC', 'AAGCCCG', 'All4mC']:
    density_col = f'density_{methyl_type}'
    bgc_vals = core_genes.loc[core_genes['is_BGC'], density_col].values
    non_bgc_vals = core_genes.loc[~core_genes['is_BGC'], density_col].values

    bgc_mean = np.mean(bgc_vals)
    non_bgc_mean = np.mean(non_bgc_vals)
    fold = bgc_mean / non_bgc_mean if non_bgc_mean > 0 else np.nan

    stat, pval = mannwhitneyu(bgc_vals, non_bgc_vals, alternative='two-sided')
    print(f"  {methyl_type}: BGC={bgc_mean:.4f}, non-BGC={non_bgc_mean:.4f}, "
          f"fold={fold:.3f}, MWU p={pval:.2e}")

    dna_vs_methyl_results.append({
        'comparison': f'Methylation_{methyl_type}_density',
        'BGC_core_mean': bgc_mean,
        'BGC_core_median': np.median(bgc_vals),
        'nonBGC_core_mean': non_bgc_mean,
        'nonBGC_core_median': np.median(non_bgc_vals),
        'fold': fold,
        'MWU_stat': stat,
        'MWU_pvalue': pval,
        'n_BGC': len(bgc_vals),
        'n_nonBGC': len(non_bgc_vals),
    })

# Interpretation: DNA motif fold vs methylation fold
print("\n  INTERPRETATION:")
# Get extended DNA fold and methylation fold for GCCGGC
dna_gccggc = [r for r in dna_vs_methyl_results
               if r['comparison'] == 'DNA_TGGCCGGC_extended'][0]
methyl_gccggc = [r for r in dna_vs_methyl_results
                  if r['comparison'] == 'Methylation_GCCGGC_density'][0]
print(f"    TGGCCGGC DNA motif fold (core): {dna_gccggc['fold']:.3f}")
print(f"    GCCGGC methylation fold (core): {methyl_gccggc['fold']:.3f}")
if dna_gccggc['fold'] > 1.1 and methyl_gccggc['fold'] > 1.1:
    if abs(dna_gccggc['fold'] - methyl_gccggc['fold']) / max(dna_gccggc['fold'], methyl_gccggc['fold']) < 0.2:
        print("    --> Both enriched and similar: SEQUENCE COMPOSITION effect")
    else:
        print("    --> Both enriched but different magnitudes: MIXED effect")
elif dna_gccggc['fold'] <= 1.1 and methyl_gccggc['fold'] > 1.1:
    print("    --> DNA normal but methylation enriched: PREFERENTIAL METHYLATION")
elif dna_gccggc['fold'] > 1.1 and methyl_gccggc['fold'] <= 1.1:
    print("    --> DNA enriched but methylation normal: PREFERENTIAL AVOIDANCE of methylation")
else:
    print("    --> Neither enriched: NO BGC-specific pattern")

dna_vs_methyl_df = pd.DataFrame(dna_vs_methyl_results)
dna_vs_methyl_df.to_csv(f'{TBL_DIR}/DNA_vs_methylation_BGC.tsv', sep='\t', index=False)
print(f"\nSaved: DNA_vs_methylation_BGC.tsv ({len(dna_vs_methyl_df)} comparisons)")

# ============================================================
# 8. Continuous density measure (Wilcoxon rank-sum)
# ============================================================
print("\n" + "=" * 70)
print("Step 7: Continuous Density Comparison (Core Only)")
print("=" * 70)

density_test_results = []
for methyl_type in ['GCCGGC', 'AAGCCCG', 'All4mC']:
    density_col = f'density_{methyl_type}'
    bgc_vals = core_genes.loc[core_genes['is_BGC'], density_col].values
    non_bgc_vals = core_genes.loc[~core_genes['is_BGC'], density_col].values

    stat, pval = mannwhitneyu(bgc_vals, non_bgc_vals, alternative='two-sided')
    # Also compute effect size (rank-biserial correlation)
    n1, n2 = len(bgc_vals), len(non_bgc_vals)
    r = 1 - (2 * stat) / (n1 * n2)

    print(f"\n  {methyl_type} density (sites/kb, 2kb extended region):")
    print(f"    BGC core: mean={np.mean(bgc_vals):.4f}, median={np.median(bgc_vals):.4f}")
    print(f"    non-BGC core: mean={np.mean(non_bgc_vals):.4f}, median={np.median(non_bgc_vals):.4f}")
    print(f"    Wilcoxon rank-sum: U={stat:.0f}, p={pval:.2e}, r={r:.4f}")

    density_test_results.append({
        'motif': methyl_type,
        'BGC_core_n': n1,
        'BGC_core_mean': np.mean(bgc_vals),
        'BGC_core_median': np.median(bgc_vals),
        'nonBGC_core_n': n2,
        'nonBGC_core_mean': np.mean(non_bgc_vals),
        'nonBGC_core_median': np.median(non_bgc_vals),
        'U_stat': stat,
        'p_value': pval,
        'rank_biserial_r': r,
    })

# ============================================================
# 9. VISUALIZATIONS
# ============================================================
print("\n" + "=" * 70)
print("Step 8: Generating Figures")
print("=" * 70)

# --- Figure 1: Chromosome ideogram ---
fig, ax = plt.subplots(1, 1, figsize=(14, 5))

# Draw chromosome
chrom_y = 0.5
chrom_height = 0.12
ax.add_patch(Rectangle((0, chrom_y - chrom_height/2), CHROM_SIZE,
                        chrom_height, facecolor='#E8E8E8', edgecolor='black', linewidth=1))

# Mark arm/core boundaries
ax.axvline(ARM_LEFT_END, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
ax.axvline(ARM_RIGHT_START, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
ax.text(ARM_LEFT_END/2, chrom_y + chrom_height/2 + 0.04, 'Left arm',
        ha='center', fontsize=8, color='gray')
ax.text((ARM_LEFT_END + ARM_RIGHT_START)/2, chrom_y + chrom_height/2 + 0.04, 'Core',
        ha='center', fontsize=8, color='gray')
ax.text((ARM_RIGHT_START + CHROM_SIZE)/2, chrom_y + chrom_height/2 + 0.04, 'Right arm',
        ha='center', fontsize=8, color='gray')

# Plot GCCGGC sites (below chromosome)
for pos in gccggc_t1_positions:
    ax.plot([pos, pos], [chrom_y - chrom_height/2 - 0.02, chrom_y - chrom_height/2 - 0.08],
            color='#2196F3', linewidth=0.3, alpha=0.5)
ax.text(-200000, chrom_y - chrom_height/2 - 0.05, 'GCCGGC\n4mC T1',
        ha='right', va='center', fontsize=7, color='#2196F3')

# Plot BGC clusters (above chromosome)
bgc_colors = {'act': '#E91E63', 'red': '#F44336', 'cda': '#9C27B0', 'cpk': '#FF9800'}
for cluster in sorted(bgc_df['bgc_name'].unique()):
    cdf = gene_annot[gene_annot['bgc_name'] == cluster]
    c_start = cdf['start'].min()
    c_end = cdf['end'].max()
    color = bgc_colors.get(cluster, '#607D8B')
    ax.add_patch(Rectangle((c_start, chrom_y + chrom_height/2 + 0.01),
                            c_end - c_start, 0.06,
                            facecolor=color, edgecolor=color, alpha=0.7))
    ax.text((c_start + c_end)/2, chrom_y + chrom_height/2 + 0.10,
            cluster.upper(), ha='center', fontsize=7, fontweight='bold', color=color)

# AAGCCCG sites (further below)
for pos in aagcccg_t1_positions:
    ax.plot([pos, pos], [chrom_y - chrom_height/2 - 0.10, chrom_y - chrom_height/2 - 0.15],
            color='#4CAF50', linewidth=0.3, alpha=0.6)
ax.text(-200000, chrom_y - chrom_height/2 - 0.125, 'AAGCCCG\n6mA T1',
        ha='right', va='center', fontsize=7, color='#4CAF50')

ax.set_xlim(-500000, CHROM_SIZE + 200000)
ax.set_ylim(0.15, 0.75)
ax.set_xlabel('Genome position (bp)', fontsize=10)
ax.set_title('S. coelicolor M145 Chromosome: BGC Clusters and Methylation Sites',
             fontsize=12, fontweight='bold')
ax.set_yticks([])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)

# Format x-axis in Mb
xticks = np.arange(0, CHROM_SIZE + 1, 1_000_000)
ax.set_xticks(xticks)
ax.set_xticklabels([f'{x/1e6:.0f}' for x in xticks])
ax.set_xlabel('Genome position (Mb)', fontsize=10)

plt.tight_layout()
plt.savefig(f'{FIG_DIR}/chromosome_BGC_GCCGGC_map.pdf', bbox_inches='tight', dpi=300)
plt.savefig(f'{FIG_DIR}/chromosome_BGC_GCCGGC_map.svg', bbox_inches='tight')
plt.close()
print("  Saved: chromosome_BGC_GCCGGC_map.pdf/svg")

# --- Figure 2: Enrichment comparison (All vs Core-only) ---
fig, axes = plt.subplots(1, 3, figsize=(14, 5))

for idx, motif in enumerate(['GCCGGC', 'AAGCCCG', 'All4mC']):
    ax = axes[idx]
    all_res = [r for r in results if r['motif'] == motif and r['subset'] == 'all_genes'][0]
    core_res = [r for r in results if r['motif'] == motif and r['subset'] == 'core_only'][0]

    labels = ['All genes', 'Core only']
    folds = [all_res['fold_enrichment'], core_res['fold_enrichment']]
    pvals = [all_res['pvalue'], core_res['pvalue']]
    colors = ['#42A5F5', '#1565C0']

    bars = ax.bar(labels, folds, color=colors, edgecolor='black', linewidth=0.5)

    # Add significance stars
    for j, (bar, pv) in enumerate(zip(bars, pvals)):
        if pv < 0.001:
            sig = '***'
        elif pv < 0.01:
            sig = '**'
        elif pv < 0.05:
            sig = '*'
        else:
            sig = 'n.s.'
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.03,
                f'{sig}\np={pv:.1e}', ha='center', va='bottom', fontsize=7)

    ax.axhline(1.0, color='red', linestyle='--', linewidth=0.8, label='No enrichment')
    ax.set_ylabel('Fold enrichment (BGC vs non-BGC)')
    ax.set_title(f'{motif}\n(methylation proximity)')
    ax.set_ylim(0, max(folds) * 1.5)

fig.suptitle('H24: BGC Methylation Enrichment -- All Genes vs Core-Only',
             fontsize=13, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(f'{FIG_DIR}/enrichment_comparison.pdf', bbox_inches='tight', dpi=300)
plt.savefig(f'{FIG_DIR}/enrichment_comparison.svg', bbox_inches='tight')
plt.close()
print("  Saved: enrichment_comparison.pdf/svg")

# --- Figure 3: Per-cluster methylation density ---
fig, ax = plt.subplots(figsize=(10, 5))

clusters = [cr['cluster'].upper() for cr in cluster_results]
gccggc_densities = [cr['GCCGGC_density_per_kb'] for cr in cluster_results]
aagcccg_densities = [cr['AAGCCCG_density_per_kb'] for cr in cluster_results]
all4mc_densities = [cr['All4mC_density_per_kb'] for cr in cluster_results]

x = np.arange(len(clusters))
width = 0.25

bars1 = ax.bar(x - width, gccggc_densities, width, label='GCCGGC 4mC',
               color='#2196F3', edgecolor='black', linewidth=0.5)
bars2 = ax.bar(x, aagcccg_densities, width, label='AAGCCCG 6mA',
               color='#4CAF50', edgecolor='black', linewidth=0.5)
bars3 = ax.bar(x + width, all4mc_densities, width, label='All 4mC',
               color='#FF9800', edgecolor='black', linewidth=0.5)

# Add reference lines
ax.axhline(core_gccggc_density, color='#2196F3', linestyle='--', linewidth=0.8, alpha=0.5,
           label=f'Core avg GCCGGC ({core_gccggc_density:.4f}/kb)')
ax.axhline(core_aagcccg_density, color='#4CAF50', linestyle='--', linewidth=0.8, alpha=0.5,
           label=f'Core avg AAGCCCG ({core_aagcccg_density:.4f}/kb)')

ax.set_xlabel('BGC Cluster')
ax.set_ylabel('Methylation site density (sites/kb)')
ax.set_title('Methylation Site Density per BGC Cluster', fontsize=12, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(clusters)
ax.legend(fontsize=8, loc='upper right')

plt.tight_layout()
plt.savefig(f'{FIG_DIR}/cluster_specific_density.pdf', bbox_inches='tight', dpi=300)
plt.savefig(f'{FIG_DIR}/cluster_specific_density.svg', bbox_inches='tight')
plt.close()
print("  Saved: cluster_specific_density.pdf/svg")

# --- Figure 4: Multi-panel comprehensive summary ---
fig = plt.figure(figsize=(16, 12))
gs = gridspec.GridSpec(2, 2, hspace=0.35, wspace=0.3)

# Panel A: BGC geography pie chart
ax_a = fig.add_subplot(gs[0, 0])
bgc_core_n = (bgc_df['region'] == 'core').sum()
bgc_arm_n = (bgc_df['region'] == 'arm').sum()
if bgc_arm_n > 0:
    sizes = [bgc_core_n, bgc_arm_n]
    labels_pie = [f'Core ({bgc_core_n})', f'Arm ({bgc_arm_n})']
    colors_pie = ['#1565C0', '#90CAF9']
else:
    sizes = [bgc_core_n]
    labels_pie = [f'Core ({bgc_core_n})']
    colors_pie = ['#1565C0']
ax_a.pie(sizes, labels=labels_pie, colors=colors_pie, autopct='%1.0f%%',
         startangle=90, textprops={'fontsize': 9})
ax_a.set_title('A. BGC Gene Geography', fontweight='bold', fontsize=11)

# Panel B: All vs Core fold enrichment
ax_b = fig.add_subplot(gs[0, 1])
motifs_display = ['GCCGGC', 'AAGCCCG', 'All4mC']
all_folds = []
core_folds = []
all_pvals = []
core_pvals = []
for motif in motifs_display:
    all_r = [r for r in results if r['motif'] == motif and r['subset'] == 'all_genes'][0]
    core_r = [r for r in results if r['motif'] == motif and r['subset'] == 'core_only'][0]
    all_folds.append(all_r['fold_enrichment'])
    core_folds.append(core_r['fold_enrichment'])
    all_pvals.append(all_r['pvalue'])
    core_pvals.append(core_r['pvalue'])

x = np.arange(len(motifs_display))
width = 0.35
bars1 = ax_b.bar(x - width/2, all_folds, width, label='All genes',
                 color='#42A5F5', edgecolor='black', linewidth=0.5)
bars2 = ax_b.bar(x + width/2, core_folds, width, label='Core only',
                 color='#1565C0', edgecolor='black', linewidth=0.5)

# Significance annotations
for j in range(len(motifs_display)):
    for bar, pv in [(bars1[j], all_pvals[j]), (bars2[j], core_pvals[j])]:
        sig = '***' if pv < 0.001 else ('**' if pv < 0.01 else ('*' if pv < 0.05 else 'n.s.'))
        ax_b.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                  sig, ha='center', va='bottom', fontsize=7)

ax_b.axhline(1.0, color='red', linestyle='--', linewidth=0.8)
ax_b.set_ylabel('Fold enrichment')
ax_b.set_title('B. BGC Enrichment: All vs Core-only', fontweight='bold', fontsize=11)
ax_b.set_xticks(x)
ax_b.set_xticklabels(motifs_display)
ax_b.legend(fontsize=8)
ax_b.set_ylim(0, max(max(all_folds), max(core_folds)) * 1.4)

# Panel C: DNA motif vs Methylation (core-only)
ax_c = fig.add_subplot(gs[1, 0])
dna_fold = dna_gccggc['fold']
methyl_fold = methyl_gccggc['fold']
bars = ax_c.bar(['DNA motif\n(TGGCCGGC)', 'Methylation\n(GCCGGC 4mC)'],
                [dna_fold, methyl_fold],
                color=['#78909C', '#2196F3'], edgecolor='black', linewidth=0.5)
# Add p-values
for bar, pv in zip(bars, [dna_gccggc['MWU_pvalue'], methyl_gccggc['MWU_pvalue']]):
    sig = '***' if pv < 0.001 else ('**' if pv < 0.01 else ('*' if pv < 0.05 else 'n.s.'))
    ax_c.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
              f'{sig}\np={pv:.1e}', ha='center', va='bottom', fontsize=7)
ax_c.axhline(1.0, color='red', linestyle='--', linewidth=0.8)
ax_c.set_ylabel('Fold (BGC vs non-BGC core)')
ax_c.set_title('C. DNA Sequence vs Methylation', fontweight='bold', fontsize=11)
ax_c.set_ylim(0, max(dna_fold, methyl_fold) * 1.6)

# Panel D: Per-cluster GCCGGC density with core average
ax_d = fig.add_subplot(gs[1, 1])
clusters_sorted = sorted(cluster_results, key=lambda x: x['GCCGGC_density_per_kb'], reverse=True)
cluster_names = [cr['cluster'].upper() for cr in clusters_sorted]
cluster_dens = [cr['GCCGGC_density_per_kb'] for cr in clusters_sorted]
cluster_colors = [bgc_colors.get(cr['cluster'], '#607D8B') for cr in clusters_sorted]

bars = ax_d.bar(range(len(cluster_names)), cluster_dens, color=cluster_colors,
                edgecolor='black', linewidth=0.5)
ax_d.axhline(core_gccggc_density, color='blue', linestyle='--', linewidth=1,
             label=f'Core average ({core_gccggc_density:.4f}/kb)')
ax_d.axhline(genome_gccggc_density, color='gray', linestyle=':', linewidth=1,
             label=f'Genome average ({genome_gccggc_density:.4f}/kb)')
ax_d.set_xticks(range(len(cluster_names)))
ax_d.set_xticklabels(cluster_names)
ax_d.set_ylabel('GCCGGC sites/kb')
ax_d.set_title('D. GCCGGC Density per BGC Cluster', fontweight='bold', fontsize=11)
ax_d.legend(fontsize=8)

# Add fold annotations
for i, cr in enumerate(clusters_sorted):
    ax_d.text(i, cr['GCCGGC_density_per_kb'] + 0.005,
              f"{cr['GCCGGC_fold_vs_core']:.1f}x",
              ha='center', va='bottom', fontsize=8)

fig.suptitle('H24: BGC Methylation Enrichment -- Geographic Confound Test',
             fontsize=14, fontweight='bold', y=0.98)
plt.savefig(f'{FIG_DIR}/H24_comprehensive_summary.pdf', bbox_inches='tight', dpi=300)
plt.savefig(f'{FIG_DIR}/H24_comprehensive_summary.svg', bbox_inches='tight')
plt.close()
print("  Saved: H24_comprehensive_summary.pdf/svg")

# ============================================================
# 10. Final Summary
# ============================================================
print("\n" + "=" * 70)
print("FINAL SUMMARY")
print("=" * 70)

# Key comparison
for motif in ['GCCGGC', 'AAGCCCG', 'All4mC']:
    all_res = [r for r in results if r['motif'] == motif and r['subset'] == 'all_genes'][0]
    core_res = [r for r in results if r['motif'] == motif and r['subset'] == 'core_only'][0]
    print(f"\n  {motif}:")
    print(f"    All genes: fold={all_res['fold_enrichment']:.3f}, p={all_res['pvalue']:.2e}")
    print(f"    Core only: fold={core_res['fold_enrichment']:.3f}, p={core_res['pvalue']:.2e}")
    ratio = core_res['fold_enrichment'] / all_res['fold_enrichment'] if all_res['fold_enrichment'] > 0 else np.nan
    print(f"    Ratio (core/all): {ratio:.3f}")
    if core_res['pvalue'] < 0.05:
        if core_res['fold_enrichment'] > 1.2:
            print(f"    --> GENUINE enrichment survives stratification")
        else:
            print(f"    --> Marginal enrichment in core")
    else:
        print(f"    --> Geographic confound CONFIRMED (like H19)")

# DNA vs methylation interpretation
print(f"\n  DNA MOTIF CONTROL:")
print(f"    TGGCCGGC DNA motif fold (core): {dna_gccggc['fold']:.3f}, p={dna_gccggc['MWU_pvalue']:.2e}")
print(f"    GCCGGC methylation fold (core): {methyl_gccggc['fold']:.3f}, p={methyl_gccggc['MWU_pvalue']:.2e}")

# Overall verdict
print(f"\n  HYPOTHESIS VERDICT:")
core_gccggc_res = [r for r in results if r['motif'] == 'GCCGGC' and r['subset'] == 'core_only'][0]
if core_gccggc_res['pvalue'] < 0.05 and core_gccggc_res['fold_enrichment'] > 1.2:
    print(f"    H24 REJECTED: BGC methylation enrichment is GENUINE (not geographic confound)")
    print(f"    Core-only fold = {core_gccggc_res['fold_enrichment']:.3f} >> 1.0")
elif core_gccggc_res['pvalue'] >= 0.05:
    print(f"    H24 SUPPORTED: Geographic confound confirmed (like H19's Simpson's paradox)")
    print(f"    Core-only fold = {core_gccggc_res['fold_enrichment']:.3f} ~ 1.0, p > 0.05")
else:
    print(f"    H24 PARTIAL: Reduced but still detectable enrichment")
    print(f"    Core-only fold = {core_gccggc_res['fold_enrichment']:.3f}")

print("\n" + "=" * 70)
print("Analysis complete.")
print("=" * 70)
