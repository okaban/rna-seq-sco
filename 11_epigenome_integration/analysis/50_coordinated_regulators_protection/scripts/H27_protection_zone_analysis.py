#!/usr/bin/env python3
"""
H27: Protection Zone Characteristics of the 57 Coordinated Regulators

Compares methylation protection zone characteristics between 57 coordinated
regulatory genes (from H6/H8) and ~998 non-coordinated regulatory genes.
Tests whether coordinated regulators have shallower protection zones
(methylation sites closer to TSS).
"""

import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Configuration
# ============================================================
BASE = Path("/Users/okaban/bioinfo/rna-seq")
ANALYSIS = BASE / "11_epigenome_integration/analysis/50_coordinated_regulators_protection"
FIG_DIR = ANALYSIS / "figures"
TBL_DIR = ANALYSIS / "tables"

CHROM_LEN = 8_667_507
ARM_LEFT = 1_500_000
ARM_RIGHT = 7_167_508

BIN_SIZE = 200  # bp
FLANK = 5000    # bp from TSS
N_BOOTSTRAP = 10_000
SEED = 42

np.random.seed(SEED)

# ============================================================
# 1. Load and classify regulatory genes
# ============================================================
print("=" * 70)
print("Step 1: Load and classify regulatory genes")
print("=" * 70)

coord_df = pd.read_csv(
    BASE / "11_epigenome_integration/analysis/29_genomewide_TF_screen/tables/coordinated_regulatory_genes.tsv",
    sep='\t'
)
all_reg_df = pd.read_csv(
    BASE / "11_epigenome_integration/analysis/29_genomewide_TF_screen/tables/all_regulatory_genes.tsv",
    sep='\t'
)

coord_loci = set(coord_df['locus_tag'])
all_reg_df['is_coordinated'] = all_reg_df['locus_tag'].isin(coord_loci)

# Calculate TSS
def get_tss(row):
    if row['strand'] == '+':
        return row['start']
    else:
        return row['end']

all_reg_df['tss'] = all_reg_df.apply(get_tss, axis=1)
coord_df['tss'] = coord_df.apply(get_tss, axis=1)

# Assign region
def assign_region(pos):
    if pos <= ARM_LEFT or pos >= ARM_RIGHT:
        return 'arm'
    return 'core'

all_reg_df['region'] = all_reg_df['tss'].apply(assign_region)

# Classify coordination types from the coordinated genes
# Use coordination_T2 and coordination_T3 columns
coord_types = []
for _, row in coord_df.iterrows():
    types = []
    for tp in ['coordination_T2', 'coordination_T3']:
        ct = row.get(tp, '')
        if pd.notna(ct) and ct not in ['ambiguous', 'methyl_change_no_expr_change', 'no_methyl_change']:
            types.append(ct)
    # Take the most specific coordination
    if 'concordant_derepression' in types:
        coord_types.append('concordant_derepression')  # lost + up
    elif 'concordant_repression' in types:
        coord_types.append('concordant_repression')    # gained + down
    elif 'discordant_gain_up' in types:
        coord_types.append('discordant_gain_up')       # gained + up
    elif 'discordant_loss_down' in types:
        coord_types.append('discordant_loss_down')     # lost + down
    else:
        coord_types.append('other')

coord_df['coordination_primary'] = coord_types

n_coord = all_reg_df['is_coordinated'].sum()
n_noncoord = (~all_reg_df['is_coordinated']).sum()
print(f"  Coordinated: {n_coord}")
print(f"  Non-coordinated: {n_noncoord}")
print(f"  Total: {len(all_reg_df)}")
print(f"\n  Coordination type breakdown:")
for ct, cnt in coord_df['coordination_primary'].value_counts().items():
    print(f"    {ct}: {cnt}")

# ============================================================
# 2. Load methylation data
# ============================================================
print("\n" + "=" * 70)
print("Step 2: Load methylation data")
print("=" * 70)

methyl_df = pd.read_csv(
    BASE / "11_epigenome_integration/analysis/01_integration/high_confidence_sites_weighted.csv"
)
print(f"  Total methylation records: {len(methyl_df)}")
print(f"  Timepoints: {methyl_df['timepoint'].unique()}")
print(f"  Mod types: {methyl_df['mod_type'].unique()}")

# Get unique sites per timepoint
methyl_sites_all = methyl_df[['position', 'timepoint']].drop_duplicates()

# All unique positions across timepoints
all_positions = methyl_df['position'].unique()
print(f"  Unique positions (all timepoints): {len(all_positions)}")

# Per timepoint
for tp in ['T1', 'T2', 'T3']:
    n = methyl_df[methyl_df['timepoint'] == tp]['position'].nunique()
    print(f"    {tp}: {n} unique positions")

# ============================================================
# 3. Per-gene nearest methylation site distance
# ============================================================
print("\n" + "=" * 70)
print("Step 3: Per-gene nearest methylation site distance")
print("=" * 70)

sorted_positions = np.sort(all_positions)

def find_nearest_distance(tss, positions):
    """Find distance to nearest methylation site."""
    idx = np.searchsorted(positions, tss)
    candidates = []
    if idx > 0:
        candidates.append(abs(tss - positions[idx - 1]))
    if idx < len(positions):
        candidates.append(abs(tss - positions[idx]))
    return min(candidates) if candidates else np.nan

def count_sites_within(tss, positions, window=2000):
    """Count methylation sites within window of TSS."""
    left = np.searchsorted(positions, tss - window, side='left')
    right = np.searchsorted(positions, tss + window, side='right')
    return right - left

# Calculate for all regulatory genes
all_reg_df['nearest_methyl_dist'] = all_reg_df['tss'].apply(
    lambda x: find_nearest_distance(x, sorted_positions)
)
all_reg_df['sites_within_2kb'] = all_reg_df['tss'].apply(
    lambda x: count_sites_within(x, sorted_positions, 2000)
)

coord_mask = all_reg_df['is_coordinated']
dist_coord = all_reg_df.loc[coord_mask, 'nearest_methyl_dist'].dropna()
dist_noncoord = all_reg_df.loc[~coord_mask, 'nearest_methyl_dist'].dropna()

stat_dist, p_dist = stats.mannwhitneyu(dist_coord, dist_noncoord, alternative='two-sided')
print(f"\n  Nearest methylation site distance:")
print(f"    Coordinated:     median = {dist_coord.median():.0f} bp, mean = {dist_coord.mean():.0f} bp")
print(f"    Non-coordinated: median = {dist_noncoord.median():.0f} bp, mean = {dist_noncoord.mean():.0f} bp")
print(f"    Mann-Whitney U = {stat_dist:.0f}, p = {p_dist:.2e}")
print(f"    Ratio (coord/noncoord): {dist_coord.median()/dist_noncoord.median():.3f}")

sites_coord = all_reg_df.loc[coord_mask, 'sites_within_2kb']
sites_noncoord = all_reg_df.loc[~coord_mask, 'sites_within_2kb']
stat_sites, p_sites = stats.mannwhitneyu(sites_coord, sites_noncoord, alternative='two-sided')
print(f"\n  Sites within 2 kb of TSS:")
print(f"    Coordinated:     median = {sites_coord.median():.1f}, mean = {sites_coord.mean():.2f}")
print(f"    Non-coordinated: median = {sites_noncoord.median():.1f}, mean = {sites_noncoord.mean():.2f}")
print(f"    Mann-Whitney U = {stat_sites:.0f}, p = {p_sites:.2e}")

# ============================================================
# 4. TSS-centered spatial profile (stratified)
# ============================================================
print("\n" + "=" * 70)
print("Step 4: TSS-centered spatial methylation profile")
print("=" * 70)

bins = np.arange(-FLANK, FLANK + BIN_SIZE, BIN_SIZE)
bin_centers = (bins[:-1] + bins[1:]) / 2

def compute_density_profile(tss_list, positions, bins, strand_list=None):
    """Compute methylation density profile around TSS."""
    counts = np.zeros(len(bins) - 1)
    for i, tss in enumerate(tss_list):
        # For - strand genes, flip the coordinate
        if strand_list is not None and strand_list[i] == '-':
            rel_pos = tss - positions  # flip for minus strand
        else:
            rel_pos = positions - tss
        # Filter to window
        mask = (rel_pos >= -FLANK) & (rel_pos <= FLANK)
        if mask.any():
            hist, _ = np.histogram(rel_pos[mask], bins=bins)
            counts += hist
    n_genes = len(tss_list)
    # density per kb per gene
    density = counts / (BIN_SIZE / 1000) / n_genes
    return density

def bootstrap_density_profile(tss_list, strand_list, positions, bins, n_iter=1000):
    """Bootstrap confidence intervals for density profile."""
    n = len(tss_list)
    boot_densities = np.zeros((n_iter, len(bins) - 1))
    for b in range(n_iter):
        idx = np.random.choice(n, n, replace=True)
        boot_tss = [tss_list[i] for i in idx]
        boot_strand = [strand_list[i] for i in idx]
        boot_densities[b] = compute_density_profile(boot_tss, positions, bins, boot_strand)
    ci_low = np.percentile(boot_densities, 2.5, axis=0)
    ci_high = np.percentile(boot_densities, 97.5, axis=0)
    return ci_low, ci_high

# Prepare data
coord_genes = all_reg_df[all_reg_df['is_coordinated']]
noncoord_genes = all_reg_df[~all_reg_df['is_coordinated']]

coord_tss = coord_genes['tss'].values.tolist()
coord_strand = coord_genes['strand'].values.tolist()
noncoord_tss = noncoord_genes['tss'].values.tolist()
noncoord_strand = noncoord_genes['strand'].values.tolist()
all_tss = all_reg_df['tss'].values.tolist()
all_strand = all_reg_df['strand'].values.tolist()

# Compute profiles
print("  Computing density profiles...")
density_coord = compute_density_profile(coord_tss, sorted_positions, bins, coord_strand)
density_noncoord = compute_density_profile(noncoord_tss, sorted_positions, bins, noncoord_strand)
density_all = compute_density_profile(all_tss, sorted_positions, bins, all_strand)

# Bootstrap CIs (1000 iterations for speed)
print("  Bootstrapping confidence intervals (1000 iterations)...")
ci_low_coord, ci_high_coord = bootstrap_density_profile(
    coord_tss, coord_strand, sorted_positions, bins, n_iter=1000
)
ci_low_noncoord, ci_high_noncoord = bootstrap_density_profile(
    noncoord_tss, noncoord_strand, sorted_positions, bins, n_iter=1000
)

# Per-bin permutation tests (1000 iterations)
print("  Per-bin permutation tests...")
all_tss_arr = np.array(all_tss)
all_strand_arr = np.array(all_strand)
n_coord_total = len(coord_tss)
n_perm = 1000
perm_diffs = np.zeros((n_perm, len(bin_centers)))

for p in range(n_perm):
    perm_idx = np.random.permutation(len(all_tss_arr))
    perm_coord_idx = perm_idx[:n_coord_total]
    perm_noncoord_idx = perm_idx[n_coord_total:]

    d_c = compute_density_profile(
        all_tss_arr[perm_coord_idx].tolist(), sorted_positions, bins,
        all_strand_arr[perm_coord_idx].tolist()
    )
    d_nc = compute_density_profile(
        all_tss_arr[perm_noncoord_idx].tolist(), sorted_positions, bins,
        all_strand_arr[perm_noncoord_idx].tolist()
    )
    perm_diffs[p] = d_c - d_nc

obs_diff = density_coord - density_noncoord
perm_pvals = np.array([
    np.mean(np.abs(perm_diffs[:, i]) >= np.abs(obs_diff[i]))
    for i in range(len(bin_centers))
])

# Find significant bins
sig_bins = perm_pvals < 0.05
print(f"  Significant bins (p < 0.05): {sig_bins.sum()} / {len(bin_centers)}")
if sig_bins.any():
    sig_centers = bin_centers[sig_bins]
    print(f"    Positions: {sig_centers.min():.0f} to {sig_centers.max():.0f} bp")

# ============================================================
# 5. Protection zone metrics comparison
# ============================================================
print("\n" + "=" * 70)
print("Step 5: Protection zone metrics comparison")
print("=" * 70)

def estimate_protection_zone(density, bin_centers, flank_range=3000):
    """Estimate protection zone width and depth."""
    # Baseline = mean density in flanking regions (>3kb from TSS)
    flank_mask = np.abs(bin_centers) > flank_range
    if flank_mask.sum() == 0:
        flank_mask = np.abs(bin_centers) > 2000
    baseline = density[flank_mask].mean()

    if baseline == 0:
        return {'width': 0, 'depth': 0, 'center': 0, 'baseline': 0, 'min_density': 0}

    # Find region where density < baseline (protection zone)
    below = density < baseline

    # Find contiguous region around TSS (bin_center = 0)
    center_idx = np.argmin(np.abs(bin_centers))

    # Expand left
    left_idx = center_idx
    while left_idx > 0 and below[left_idx - 1]:
        left_idx -= 1

    # Expand right
    right_idx = center_idx
    while right_idx < len(bin_centers) - 1 and below[right_idx + 1]:
        right_idx += 1

    if not below[center_idx]:
        # No protection at TSS
        return {
            'width': 0,
            'depth': 0,
            'center': 0,
            'baseline': baseline,
            'min_density': density[center_idx],
            'start': 0,
            'end': 0
        }

    width = bin_centers[right_idx] - bin_centers[left_idx] + BIN_SIZE
    min_density = density[left_idx:right_idx+1].min()
    min_idx = left_idx + np.argmin(density[left_idx:right_idx+1])
    depth = (baseline - min_density) / baseline  # fractional depletion

    return {
        'width': width,
        'depth': depth,
        'center': bin_centers[min_idx],
        'baseline': baseline,
        'min_density': min_density,
        'start': bin_centers[left_idx] - BIN_SIZE/2,
        'end': bin_centers[right_idx] + BIN_SIZE/2
    }

pz_coord = estimate_protection_zone(density_coord, bin_centers)
pz_noncoord = estimate_protection_zone(density_noncoord, bin_centers)
pz_all = estimate_protection_zone(density_all, bin_centers)

print(f"\n  Protection zone metrics:")
print(f"  {'Metric':<30} {'Coordinated':>15} {'Non-coordinated':>18} {'All':>12}")
print(f"  {'-'*75}")
print(f"  {'Width (bp)':<30} {pz_coord['width']:>15.0f} {pz_noncoord['width']:>18.0f} {pz_all['width']:>12.0f}")
print(f"  {'Depth (fractional)':<30} {pz_coord['depth']:>15.3f} {pz_noncoord['depth']:>18.3f} {pz_all['depth']:>12.3f}")
print(f"  {'Min density position (bp)':<30} {pz_coord['center']:>15.0f} {pz_noncoord['center']:>18.0f} {pz_all['center']:>12.0f}")
print(f"  {'Baseline density':<30} {pz_coord['baseline']:>15.3f} {pz_noncoord['baseline']:>18.3f} {pz_all['baseline']:>12.3f}")
print(f"  {'Min density':<30} {pz_coord['min_density']:>15.3f} {pz_noncoord['min_density']:>18.3f} {pz_all['min_density']:>12.3f}")
if pz_coord['width'] > 0:
    print(f"  {'Zone start (bp)':<30} {pz_coord['start']:>15.0f} {pz_noncoord.get('start',0):>18.0f} {pz_all.get('start',0):>12.0f}")
    print(f"  {'Zone end (bp)':<30} {pz_coord['end']:>15.0f} {pz_noncoord.get('end',0):>18.0f} {pz_all.get('end',0):>12.0f}")

# ============================================================
# 6. Coordination type sub-analysis
# ============================================================
print("\n" + "=" * 70)
print("Step 6: Coordination type analysis")
print("=" * 70)

# Merge coordination info
coord_merged = coord_df[['locus_tag', 'coordination_primary']].copy()
all_reg_df = all_reg_df.merge(coord_merged, on='locus_tag', how='left')
all_reg_df['coordination_primary'] = all_reg_df['coordination_primary'].fillna('non_coordinated')

# Per-type analysis
type_results = []
for ctype in coord_df['coordination_primary'].unique():
    subset = all_reg_df[all_reg_df['coordination_primary'] == ctype]
    n = len(subset)

    subset_tss = subset['tss'].values.tolist()
    subset_strand = subset['strand'].values.tolist()

    if n >= 3:
        dens = compute_density_profile(subset_tss, sorted_positions, bins, subset_strand)
        pz = estimate_protection_zone(dens, bin_centers)

        dist_vals = subset['nearest_methyl_dist'].dropna()
        sites_vals = subset['sites_within_2kb']

        type_results.append({
            'coordination_type': ctype,
            'n_genes': n,
            'median_nearest_dist': dist_vals.median(),
            'mean_nearest_dist': dist_vals.mean(),
            'median_sites_2kb': sites_vals.median(),
            'mean_sites_2kb': sites_vals.mean(),
            'protection_zone_width': pz['width'],
            'protection_zone_depth': pz['depth'],
            'pz_center': pz['center'],
            'baseline_density': pz['baseline'],
            'min_density': pz['min_density']
        })

        print(f"\n  {ctype} (n={n}):")
        print(f"    Nearest dist: median={dist_vals.median():.0f} bp, mean={dist_vals.mean():.0f} bp")
        print(f"    Sites <=2kb: median={sites_vals.median():.1f}, mean={sites_vals.mean():.2f}")
        print(f"    Protection zone: width={pz['width']:.0f} bp, depth={pz['depth']:.3f}")

# Add non-coordinated as reference
noncoord_sub = all_reg_df[~all_reg_df['is_coordinated']]
type_results.append({
    'coordination_type': 'non_coordinated',
    'n_genes': len(noncoord_sub),
    'median_nearest_dist': noncoord_sub['nearest_methyl_dist'].dropna().median(),
    'mean_nearest_dist': noncoord_sub['nearest_methyl_dist'].dropna().mean(),
    'median_sites_2kb': noncoord_sub['sites_within_2kb'].median(),
    'mean_sites_2kb': noncoord_sub['sites_within_2kb'].mean(),
    'protection_zone_width': pz_noncoord['width'],
    'protection_zone_depth': pz_noncoord['depth'],
    'pz_center': pz_noncoord['center'],
    'baseline_density': pz_noncoord['baseline'],
    'min_density': pz_noncoord['min_density']
})

type_results_df = pd.DataFrame(type_results)

# ============================================================
# 7. Methylation dynamics at coordinated gene promoters
# ============================================================
print("\n" + "=" * 70)
print("Step 7: Methylation dynamics (T1 -> T2)")
print("=" * 70)

# Per-timepoint positions
t1_positions = np.sort(methyl_df[methyl_df['timepoint'] == 'T1']['position'].unique())
t2_positions = np.sort(methyl_df[methyl_df['timepoint'] == 'T2']['position'].unique())

all_reg_df['sites_2kb_T1'] = all_reg_df['tss'].apply(
    lambda x: count_sites_within(x, t1_positions, 2000)
)
all_reg_df['sites_2kb_T2'] = all_reg_df['tss'].apply(
    lambda x: count_sites_within(x, t2_positions, 2000)
)
all_reg_df['delta_sites'] = all_reg_df['sites_2kb_T2'] - all_reg_df['sites_2kb_T1']

# Get LFC from DESeq2
deseq = pd.read_csv(
    BASE / "04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_2_vs_1.tsv",
    sep='\t'
)
deseq = deseq.rename(columns={'gene_id': 'locus_tag', 'log2FoldChange': 'lfc_T2vsT1', 'padj': 'padj_T2vsT1'})

all_reg_df = all_reg_df.merge(deseq[['locus_tag', 'lfc_T2vsT1', 'padj_T2vsT1']], on='locus_tag', how='left')

# Compare delta_sites between groups
delta_coord = all_reg_df.loc[all_reg_df['is_coordinated'], 'delta_sites']
delta_noncoord = all_reg_df.loc[~all_reg_df['is_coordinated'], 'delta_sites']
stat_delta, p_delta = stats.mannwhitneyu(delta_coord, delta_noncoord, alternative='two-sided')

print(f"\n  Change in methylation site count (T2 - T1, TSS +/-2kb):")
print(f"    Coordinated:     median = {delta_coord.median():.1f}, mean = {delta_coord.mean():.2f}")
print(f"    Non-coordinated: median = {delta_noncoord.median():.1f}, mean = {delta_noncoord.mean():.2f}")
print(f"    Mann-Whitney p = {p_delta:.2e}")

# Correlation: delta_sites vs LFC for coordinated genes
coord_with_lfc = all_reg_df[all_reg_df['is_coordinated'] & all_reg_df['lfc_T2vsT1'].notna()].copy()
if len(coord_with_lfc) > 3:
    rho_delta, p_rho_delta = stats.spearmanr(coord_with_lfc['delta_sites'], coord_with_lfc['lfc_T2vsT1'])
    print(f"\n  Correlation (delta_sites vs LFC_T2vsT1) for 57 coordinated:")
    print(f"    Spearman rho = {rho_delta:.3f}, p = {p_rho_delta:.2e}")

# Same for non-coordinated
noncoord_with_lfc = all_reg_df[~all_reg_df['is_coordinated'] & all_reg_df['lfc_T2vsT1'].notna()].copy()
if len(noncoord_with_lfc) > 3:
    rho_nc, p_nc = stats.spearmanr(noncoord_with_lfc['delta_sites'], noncoord_with_lfc['lfc_T2vsT1'])
    print(f"  Correlation for non-coordinated:")
    print(f"    Spearman rho = {rho_nc:.3f}, p = {p_nc:.2e}")

# ============================================================
# 8. Geographic distribution
# ============================================================
print("\n" + "=" * 70)
print("Step 8: Geographic distribution")
print("=" * 70)

coord_geo = all_reg_df[all_reg_df['is_coordinated']]['region'].value_counts()
noncoord_geo = all_reg_df[~all_reg_df['is_coordinated']]['region'].value_counts()

n_coord_arm = coord_geo.get('arm', 0)
n_coord_core = coord_geo.get('core', 0)
n_noncoord_arm = noncoord_geo.get('arm', 0)
n_noncoord_core = noncoord_geo.get('core', 0)

contingency = [[n_coord_arm, n_coord_core],
               [n_noncoord_arm, n_noncoord_core]]
or_geo, p_geo = stats.fisher_exact(contingency)

print(f"\n  Geographic distribution:")
print(f"    Coordinated:     arm={n_coord_arm} ({100*n_coord_arm/n_coord:.1f}%), core={n_coord_core} ({100*n_coord_core/n_coord:.1f}%)")
print(f"    Non-coordinated: arm={n_noncoord_arm} ({100*n_noncoord_arm/n_noncoord:.1f}%), core={n_noncoord_core} ({100*n_noncoord_core/n_noncoord:.1f}%)")
print(f"    Fisher's exact: OR = {or_geo:.3f}, p = {p_geo:.3f}")

# ============================================================
# 9. Bootstrap confidence intervals
# ============================================================
print("\n" + "=" * 70)
print("Step 9: Bootstrap CIs for mean nearest-site distance difference")
print("=" * 70)

coord_dists = dist_coord.values
noncoord_dists = dist_noncoord.values
obs_mean_diff = coord_dists.mean() - noncoord_dists.mean()

boot_diffs = np.zeros(N_BOOTSTRAP)
for b in range(N_BOOTSTRAP):
    bc = np.random.choice(coord_dists, len(coord_dists), replace=True)
    bnc = np.random.choice(noncoord_dists, len(noncoord_dists), replace=True)
    boot_diffs[b] = bc.mean() - bnc.mean()

ci_95 = np.percentile(boot_diffs, [2.5, 97.5])
print(f"\n  Observed mean difference (coord - noncoord): {obs_mean_diff:.1f} bp")
print(f"  Bootstrap 95% CI: [{ci_95[0]:.1f}, {ci_95[1]:.1f}] bp")
print(f"  CI includes zero: {ci_95[0] <= 0 <= ci_95[1]}")

# Bootstrap p-value
boot_p = np.mean(np.abs(boot_diffs - boot_diffs.mean()) >= np.abs(obs_mean_diff - boot_diffs.mean()))
print(f"  Bootstrap two-tailed p-value: {boot_p:.4f}")

# ============================================================
# 10. Save tables
# ============================================================
print("\n" + "=" * 70)
print("Step 10: Saving tables")
print("=" * 70)

# Gene-level metrics
gene_metrics = all_reg_df[[
    'locus_tag', 'gene_name', 'old_locus_tag', 'product', 'tf_family',
    'start', 'end', 'strand', 'tss', 'region', 'is_coordinated',
    'coordination_primary', 'nearest_methyl_dist', 'sites_within_2kb',
    'sites_2kb_T1', 'sites_2kb_T2', 'delta_sites', 'lfc_T2vsT1', 'padj_T2vsT1'
]].copy()
gene_metrics.to_csv(TBL_DIR / 'gene_level_metrics.tsv', sep='\t', index=False)
print(f"  Saved gene_level_metrics.tsv ({len(gene_metrics)} genes)")

# Group comparison
group_comp = pd.DataFrame([
    {
        'comparison': 'Nearest methylation site distance',
        'coordinated_median': dist_coord.median(),
        'coordinated_mean': dist_coord.mean(),
        'coordinated_sd': dist_coord.std(),
        'noncoordinated_median': dist_noncoord.median(),
        'noncoordinated_mean': dist_noncoord.mean(),
        'noncoordinated_sd': dist_noncoord.std(),
        'stat': stat_dist,
        'p_value': p_dist,
        'test': 'Mann-Whitney U',
        'bootstrap_mean_diff': obs_mean_diff,
        'bootstrap_ci_low': ci_95[0],
        'bootstrap_ci_high': ci_95[1]
    },
    {
        'comparison': 'Sites within 2kb of TSS',
        'coordinated_median': sites_coord.median(),
        'coordinated_mean': sites_coord.mean(),
        'coordinated_sd': sites_coord.std(),
        'noncoordinated_median': sites_noncoord.median(),
        'noncoordinated_mean': sites_noncoord.mean(),
        'noncoordinated_sd': sites_noncoord.std(),
        'stat': stat_sites,
        'p_value': p_sites,
        'test': 'Mann-Whitney U',
        'bootstrap_mean_diff': np.nan,
        'bootstrap_ci_low': np.nan,
        'bootstrap_ci_high': np.nan
    },
    {
        'comparison': 'Delta methylation sites T2-T1 (2kb)',
        'coordinated_median': delta_coord.median(),
        'coordinated_mean': delta_coord.mean(),
        'coordinated_sd': delta_coord.std(),
        'noncoordinated_median': delta_noncoord.median(),
        'noncoordinated_mean': delta_noncoord.mean(),
        'noncoordinated_sd': delta_noncoord.std(),
        'stat': stat_delta,
        'p_value': p_delta,
        'test': 'Mann-Whitney U',
        'bootstrap_mean_diff': np.nan,
        'bootstrap_ci_low': np.nan,
        'bootstrap_ci_high': np.nan
    },
    {
        'comparison': 'Protection zone width',
        'coordinated_median': pz_coord['width'],
        'coordinated_mean': pz_coord['width'],
        'coordinated_sd': np.nan,
        'noncoordinated_median': pz_noncoord['width'],
        'noncoordinated_mean': pz_noncoord['width'],
        'noncoordinated_sd': np.nan,
        'stat': np.nan,
        'p_value': np.nan,
        'test': 'Profile-based estimate',
        'bootstrap_mean_diff': np.nan,
        'bootstrap_ci_low': np.nan,
        'bootstrap_ci_high': np.nan
    },
    {
        'comparison': 'Protection zone depth',
        'coordinated_median': pz_coord['depth'],
        'coordinated_mean': pz_coord['depth'],
        'coordinated_sd': np.nan,
        'noncoordinated_median': pz_noncoord['depth'],
        'noncoordinated_mean': pz_noncoord['depth'],
        'noncoordinated_sd': np.nan,
        'stat': np.nan,
        'p_value': np.nan,
        'test': 'Profile-based estimate',
        'bootstrap_mean_diff': np.nan,
        'bootstrap_ci_low': np.nan,
        'bootstrap_ci_high': np.nan
    },
    {
        'comparison': 'Geographic (arm fraction)',
        'coordinated_median': n_coord_arm / n_coord,
        'coordinated_mean': n_coord_arm / n_coord,
        'coordinated_sd': np.nan,
        'noncoordinated_median': n_noncoord_arm / n_noncoord,
        'noncoordinated_mean': n_noncoord_arm / n_noncoord,
        'noncoordinated_sd': np.nan,
        'stat': or_geo,
        'p_value': p_geo,
        'test': "Fisher's exact (OR)",
        'bootstrap_mean_diff': np.nan,
        'bootstrap_ci_low': np.nan,
        'bootstrap_ci_high': np.nan
    }
])
group_comp.to_csv(TBL_DIR / 'group_comparison.tsv', sep='\t', index=False)
print(f"  Saved group_comparison.tsv")

# Coordination type analysis
type_results_df.to_csv(TBL_DIR / 'coordination_type_analysis.tsv', sep='\t', index=False)
print(f"  Saved coordination_type_analysis.tsv")

# Methylation dynamics
methyl_dynamics = all_reg_df[[
    'locus_tag', 'is_coordinated', 'coordination_primary', 'region',
    'sites_2kb_T1', 'sites_2kb_T2', 'delta_sites', 'lfc_T2vsT1'
]].copy()
methyl_dynamics.to_csv(TBL_DIR / 'methylation_dynamics.tsv', sep='\t', index=False)
print(f"  Saved methylation_dynamics.tsv")

# Spatial profile data
profile_data = pd.DataFrame({
    'bin_center_bp': bin_centers,
    'density_coordinated': density_coord,
    'density_noncoordinated': density_noncoord,
    'density_all': density_all,
    'ci_low_coordinated': ci_low_coord,
    'ci_high_coordinated': ci_high_coord,
    'ci_low_noncoordinated': ci_low_noncoord,
    'ci_high_noncoordinated': ci_high_noncoord,
    'obs_diff': obs_diff,
    'perm_pval': perm_pvals
})
profile_data.to_csv(TBL_DIR / 'spatial_profile_data.tsv', sep='\t', index=False)
print(f"  Saved spatial_profile_data.tsv")

# ============================================================
# 11. Visualization
# ============================================================
print("\n" + "=" * 70)
print("Step 11: Generating figures")
print("=" * 70)

# Color scheme
C_COORD = '#E63946'       # red for coordinated
C_NONCOORD = '#457B9D'    # steel blue for non-coordinated
C_ALL = '#2D3436'         # dark for all
C_COORD_FILL = '#E6394630'
C_NONCOORD_FILL = '#457B9D30'

# --- Panel A: TSS-centered methylation profile ---
fig, ax = plt.subplots(figsize=(10, 6))
ax.fill_between(bin_centers, ci_low_coord, ci_high_coord, alpha=0.15, color=C_COORD)
ax.fill_between(bin_centers, ci_low_noncoord, ci_high_noncoord, alpha=0.15, color=C_NONCOORD)
ax.plot(bin_centers, density_coord, color=C_COORD, linewidth=2, label=f'Coordinated (n={n_coord})')
ax.plot(bin_centers, density_noncoord, color=C_NONCOORD, linewidth=2, label=f'Non-coordinated (n={n_noncoord})')
ax.plot(bin_centers, density_all, color=C_ALL, linewidth=1, linestyle='--', alpha=0.5, label=f'All regulatory (n={len(all_reg_df)})')

# Mark significant bins
if sig_bins.any():
    y_mark = ax.get_ylim()[0] + 0.01 * (ax.get_ylim()[1] - ax.get_ylim()[0])
    sig_x = bin_centers[sig_bins]
    ax.scatter(sig_x, np.full_like(sig_x, y_mark), marker='|', color='red', s=30, alpha=0.5, zorder=5)

ax.axvline(0, color='gray', linestyle=':', alpha=0.5)
ax.set_xlabel('Distance from TSS (bp)', fontsize=12)
ax.set_ylabel('Methylation sites per kb per gene', fontsize=12)
ax.set_title('H27: TSS-Centered Methylation Profile\nCoordinated vs Non-coordinated Regulatory Genes', fontsize=13)
ax.legend(fontsize=10, loc='upper right')
ax.set_xlim(-FLANK, FLANK)

# Annotate protection zones
if pz_coord['width'] > 0:
    ax.axhline(pz_coord['baseline'], color=C_COORD, linestyle=':', alpha=0.3)
if pz_noncoord['width'] > 0:
    ax.axhline(pz_noncoord['baseline'], color=C_NONCOORD, linestyle=':', alpha=0.3)

plt.tight_layout()
fig.savefig(FIG_DIR / 'protection_zone_comparison.pdf', dpi=300, bbox_inches='tight')
fig.savefig(FIG_DIR / 'protection_zone_comparison.svg', bbox_inches='tight')
plt.close()
print("  Saved protection_zone_comparison.pdf/svg")

# --- Panel B: Nearest site distance distribution ---
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Violin/box plot for distance
data_violin = [dist_coord.values, dist_noncoord.values]
parts = axes[0].violinplot(data_violin, positions=[1, 2], showmeans=True, showmedians=True)
for i, pc in enumerate(parts['bodies']):
    pc.set_facecolor([C_COORD, C_NONCOORD][i])
    pc.set_alpha(0.4)
parts['cmeans'].set_color('black')
parts['cmedians'].set_color('red')

axes[0].set_xticks([1, 2])
axes[0].set_xticklabels(['Coordinated\n(n=57)', 'Non-coordinated\n(n=998)'])
axes[0].set_ylabel('Nearest methylation site distance (bp)', fontsize=11)
axes[0].set_title(f'Nearest Site Distance\nMann-Whitney p = {p_dist:.2e}', fontsize=12)

# Overlay box plots
bp = axes[0].boxplot(data_violin, positions=[1, 2], widths=0.15, patch_artist=True,
                     medianprops=dict(color='red', linewidth=2),
                     boxprops=dict(facecolor='white', alpha=0.7))

# Sites within 2kb
data_sites = [sites_coord.values, sites_noncoord.values]
parts2 = axes[1].violinplot(data_sites, positions=[1, 2], showmeans=True, showmedians=True)
for i, pc in enumerate(parts2['bodies']):
    pc.set_facecolor([C_COORD, C_NONCOORD][i])
    pc.set_alpha(0.4)
parts2['cmeans'].set_color('black')
parts2['cmedians'].set_color('red')

axes[1].set_xticks([1, 2])
axes[1].set_xticklabels(['Coordinated\n(n=57)', 'Non-coordinated\n(n=998)'])
axes[1].set_ylabel('Methylation sites within 2 kb of TSS', fontsize=11)
axes[1].set_title(f'Methylation Site Density Near TSS\nMann-Whitney p = {p_sites:.2e}', fontsize=12)

bp2 = axes[1].boxplot(data_sites, positions=[1, 2], widths=0.15, patch_artist=True,
                      medianprops=dict(color='red', linewidth=2),
                      boxprops=dict(facecolor='white', alpha=0.7))

plt.tight_layout()
fig.savefig(FIG_DIR / 'nearest_site_distance.pdf', dpi=300, bbox_inches='tight')
fig.savefig(FIG_DIR / 'nearest_site_distance.svg', bbox_inches='tight')
plt.close()
print("  Saved nearest_site_distance.pdf/svg")

# --- Panel C: Methylation change vs expression change ---
fig, ax = plt.subplots(figsize=(8, 7))

# Plot non-coordinated in background
mask_nc = ~all_reg_df['is_coordinated'] & all_reg_df['lfc_T2vsT1'].notna()
ax.scatter(all_reg_df.loc[mask_nc, 'delta_sites'], all_reg_df.loc[mask_nc, 'lfc_T2vsT1'],
           alpha=0.1, s=15, color=C_NONCOORD, label=f'Non-coordinated (n={mask_nc.sum()})')

# Plot coordinated genes
mask_c = all_reg_df['is_coordinated'] & all_reg_df['lfc_T2vsT1'].notna()
colors_map = {
    'concordant_derepression': '#2ECC71',    # green - lost+up
    'concordant_repression': '#E74C3C',       # red - gained+down
    'discordant_gain_up': '#F39C12',          # orange - gained+up
    'discordant_loss_down': '#9B59B6',        # purple - lost+down
    'other': '#95A5A6'                        # gray
}
for ctype, color in colors_map.items():
    mask_ct = all_reg_df['coordination_primary'] == ctype
    if mask_ct.any():
        ax.scatter(all_reg_df.loc[mask_ct, 'delta_sites'], all_reg_df.loc[mask_ct, 'lfc_T2vsT1'],
                   s=60, color=color, edgecolors='black', linewidth=0.5, zorder=5,
                   label=f'{ctype} (n={mask_ct.sum()})')

ax.axhline(0, color='gray', linestyle=':', alpha=0.5)
ax.axvline(0, color='gray', linestyle=':', alpha=0.5)
ax.set_xlabel('Change in methylation sites (T2 - T1, TSS +/- 2kb)', fontsize=12)
ax.set_ylabel('Expression change (log2FC T2 vs T1)', fontsize=12)

# Add correlation annotation for coordinated
if len(coord_with_lfc) > 3:
    ax.annotate(f'Coordinated: rho={rho_delta:.3f}, p={p_rho_delta:.2e}',
                xy=(0.02, 0.98), xycoords='axes fraction', va='top', fontsize=10,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.5))

ax.set_title('H27: Methylation Change vs Expression Change\nRegulatory Genes (TSS +/- 2kb)', fontsize=13)
ax.legend(fontsize=8, loc='lower right', ncol=2)
plt.tight_layout()
fig.savefig(FIG_DIR / 'methylation_expression_scatter.pdf', dpi=300, bbox_inches='tight')
fig.savefig(FIG_DIR / 'methylation_expression_scatter.svg', bbox_inches='tight')
plt.close()
print("  Saved methylation_expression_scatter.pdf/svg")

# --- Panel D: Geographic distribution ---
fig, axes = plt.subplots(1, 2, figsize=(10, 5))

# Bar chart
categories = ['Coordinated\n(n=57)', 'Non-coordinated\n(n=998)']
arm_pcts = [100 * n_coord_arm / n_coord, 100 * n_noncoord_arm / n_noncoord]
core_pcts = [100 * n_coord_core / n_coord, 100 * n_noncoord_core / n_noncoord]

x = np.arange(2)
width = 0.35
bars1 = axes[0].bar(x, arm_pcts, width, label='Arm', color='#E17055')
bars2 = axes[0].bar(x, core_pcts, width, bottom=arm_pcts, label='Core', color='#74B9FF')
axes[0].set_xticks(x)
axes[0].set_xticklabels(categories)
axes[0].set_ylabel('Percentage (%)')
axes[0].set_title(f"Geographic Distribution\nFisher's exact OR={or_geo:.2f}, p={p_geo:.3f}", fontsize=12)
axes[0].legend()
axes[0].set_ylim(0, 110)

# Add labels
for i, (a, c) in enumerate(zip(arm_pcts, core_pcts)):
    axes[0].text(i, a/2, f'{a:.1f}%', ha='center', va='center', fontsize=10, fontweight='bold')
    axes[0].text(i, a + c/2, f'{c:.1f}%', ha='center', va='center', fontsize=10, fontweight='bold')

# Chromosome map
axes[1].set_xlim(0, CHROM_LEN)
axes[1].set_ylim(-1, 3)

# Draw chromosome
axes[1].barh(1, CHROM_LEN, 0.3, color='#dfe6e9', edgecolor='black', linewidth=0.5)
# Mark arms
axes[1].barh(1, ARM_LEFT, 0.3, color='#fab1a0', alpha=0.5)
axes[1].barh(1, CHROM_LEN - ARM_RIGHT, 0.3, left=ARM_RIGHT, color='#fab1a0', alpha=0.5)

# Plot gene positions
coord_positions = all_reg_df[all_reg_df['is_coordinated']]['tss'].values
noncoord_positions = all_reg_df[~all_reg_df['is_coordinated']]['tss'].values

axes[1].scatter(noncoord_positions, np.full_like(noncoord_positions, 0.2, dtype=float),
                s=3, color=C_NONCOORD, alpha=0.3, label='Non-coordinated')
axes[1].scatter(coord_positions, np.full_like(coord_positions, 1.8, dtype=float),
                s=15, color=C_COORD, marker='v', zorder=5, label='Coordinated')

axes[1].set_xlabel('Chromosome position (bp)', fontsize=11)
axes[1].set_title('Chromosomal Distribution', fontsize=12)
axes[1].legend(fontsize=8, loc='upper center')
axes[1].set_yticks([])

plt.tight_layout()
fig.savefig(FIG_DIR / 'geographic_distribution.pdf', dpi=300, bbox_inches='tight')
fig.savefig(FIG_DIR / 'geographic_distribution.svg', bbox_inches='tight')
plt.close()
print("  Saved geographic_distribution.pdf/svg")

# --- Comprehensive summary (4 panels) ---
fig = plt.figure(figsize=(16, 14))
gs = gridspec.GridSpec(2, 2, hspace=0.35, wspace=0.3)

# Panel A: Profile
ax_a = fig.add_subplot(gs[0, 0])
ax_a.fill_between(bin_centers, ci_low_coord, ci_high_coord, alpha=0.15, color=C_COORD)
ax_a.fill_between(bin_centers, ci_low_noncoord, ci_high_noncoord, alpha=0.15, color=C_NONCOORD)
ax_a.plot(bin_centers, density_coord, color=C_COORD, linewidth=2, label=f'Coordinated (n={n_coord})')
ax_a.plot(bin_centers, density_noncoord, color=C_NONCOORD, linewidth=2, label=f'Non-coordinated (n={n_noncoord})')
if sig_bins.any():
    y_mark = ax_a.get_ylim()[0] + 0.01 * (ax_a.get_ylim()[1] - ax_a.get_ylim()[0])
    ax_a.scatter(bin_centers[sig_bins], np.full(sig_bins.sum(), ax_a.get_ylim()[0]),
                 marker='|', color='red', s=20, alpha=0.5)
ax_a.axvline(0, color='gray', linestyle=':', alpha=0.5)
ax_a.set_xlabel('Distance from TSS (bp)')
ax_a.set_ylabel('Sites per kb per gene')
ax_a.set_title('A. TSS-Centered Methylation Profile')
ax_a.legend(fontsize=8)
ax_a.set_xlim(-FLANK, FLANK)

# Panel B: Distance distribution
ax_b = fig.add_subplot(gs[0, 1])
parts_b = ax_b.violinplot(data_violin, positions=[1, 2], showmeans=True, showmedians=True)
for i, pc in enumerate(parts_b['bodies']):
    pc.set_facecolor([C_COORD, C_NONCOORD][i])
    pc.set_alpha(0.4)
parts_b['cmeans'].set_color('black')
parts_b['cmedians'].set_color('red')
ax_b.boxplot(data_violin, positions=[1, 2], widths=0.15, patch_artist=True,
             medianprops=dict(color='red', linewidth=2),
             boxprops=dict(facecolor='white', alpha=0.7))
ax_b.set_xticks([1, 2])
ax_b.set_xticklabels(['Coordinated', 'Non-coord.'])
ax_b.set_ylabel('Nearest methylation site (bp)')
ax_b.set_title(f'B. Nearest Site Distance\np = {p_dist:.2e}')

# Panel C: Scatter
ax_c = fig.add_subplot(gs[1, 0])
ax_c.scatter(all_reg_df.loc[mask_nc, 'delta_sites'], all_reg_df.loc[mask_nc, 'lfc_T2vsT1'],
             alpha=0.08, s=12, color=C_NONCOORD)
for ctype, color in colors_map.items():
    mask_ct = all_reg_df['coordination_primary'] == ctype
    if mask_ct.any():
        ax_c.scatter(all_reg_df.loc[mask_ct, 'delta_sites'], all_reg_df.loc[mask_ct, 'lfc_T2vsT1'],
                     s=50, color=color, edgecolors='black', linewidth=0.5, zorder=5,
                     label=ctype.replace('_', ' '))
ax_c.axhline(0, color='gray', linestyle=':', alpha=0.5)
ax_c.axvline(0, color='gray', linestyle=':', alpha=0.5)
ax_c.set_xlabel('Delta methylation sites (T2-T1)')
ax_c.set_ylabel('Expression LFC (T2 vs T1)')
if len(coord_with_lfc) > 3:
    ax_c.set_title(f'C. Methylation vs Expression Change\nrho={rho_delta:.3f}, p={p_rho_delta:.2e}')
else:
    ax_c.set_title('C. Methylation vs Expression Change')
ax_c.legend(fontsize=7, loc='lower right')

# Panel D: Geographic
ax_d = fig.add_subplot(gs[1, 1])
bars1 = ax_d.bar(x, arm_pcts, width, label='Arm', color='#E17055')
bars2 = ax_d.bar(x, core_pcts, width, bottom=arm_pcts, label='Core', color='#74B9FF')
ax_d.set_xticks(x)
ax_d.set_xticklabels(['Coordinated', 'Non-coord.'])
ax_d.set_ylabel('Percentage (%)')
ax_d.set_title(f"D. Geographic Distribution\nOR={or_geo:.2f}, p={p_geo:.3f}")
ax_d.legend()
ax_d.set_ylim(0, 110)
for i, (a, c) in enumerate(zip(arm_pcts, core_pcts)):
    ax_d.text(i, a/2, f'{a:.1f}%', ha='center', va='center', fontsize=10, fontweight='bold')
    ax_d.text(i, a + c/2, f'{c:.1f}%', ha='center', va='center', fontsize=10, fontweight='bold')

fig.suptitle('H27: Protection Zone Characteristics of 57 Coordinated Regulators',
             fontsize=15, fontweight='bold', y=0.98)
plt.savefig(FIG_DIR / 'H27_comprehensive_summary.pdf', dpi=300, bbox_inches='tight')
plt.savefig(FIG_DIR / 'H27_comprehensive_summary.svg', bbox_inches='tight')
plt.close()
print("  Saved H27_comprehensive_summary.pdf/svg")

# ============================================================
# 12. Final Summary
# ============================================================
print("\n" + "=" * 70)
print("FINAL SUMMARY")
print("=" * 70)
print(f"""
H27: Protection Zone Characteristics of 57 Coordinated Regulators
===================================================================

1. NEAREST METHYLATION SITE DISTANCE:
   Coordinated:     median = {dist_coord.median():.0f} bp, mean = {dist_coord.mean():.0f} bp
   Non-coordinated: median = {dist_noncoord.median():.0f} bp, mean = {dist_noncoord.mean():.0f} bp
   Mann-Whitney p = {p_dist:.2e}
   Bootstrap 95% CI of mean difference: [{ci_95[0]:.0f}, {ci_95[1]:.0f}] bp

2. SITES WITHIN 2 KB OF TSS:
   Coordinated:     median = {sites_coord.median():.1f}, mean = {sites_coord.mean():.2f}
   Non-coordinated: median = {sites_noncoord.median():.1f}, mean = {sites_noncoord.mean():.2f}
   Mann-Whitney p = {p_sites:.2e}

3. PROTECTION ZONE COMPARISON:
   Coordinated:     width = {pz_coord['width']:.0f} bp, depth = {pz_coord['depth']:.3f}
   Non-coordinated: width = {pz_noncoord['width']:.0f} bp, depth = {pz_noncoord['depth']:.3f}

4. METHYLATION DYNAMICS (T1 -> T2):
   Coordinated delta:     median = {delta_coord.median():.1f}, mean = {delta_coord.mean():.2f}
   Non-coordinated delta: median = {delta_noncoord.median():.1f}, mean = {delta_noncoord.mean():.2f}
   Mann-Whitney p = {p_delta:.2e}
   Correlation (delta vs LFC, coordinated): rho = {rho_delta:.3f}, p = {p_rho_delta:.2e}

5. GEOGRAPHIC DISTRIBUTION:
   Coordinated arm fraction:     {100*n_coord_arm/n_coord:.1f}%
   Non-coordinated arm fraction: {100*n_noncoord_arm/n_noncoord:.1f}%
   Fisher's exact OR = {or_geo:.3f}, p = {p_geo:.3f}

6. SIGNIFICANT BINS IN PROFILE:
   {sig_bins.sum()} / {len(bin_centers)} bins at p < 0.05

7. HYPOTHESIS ASSESSMENT:
   The hypothesis that coordinated regulators have SHALLOWER protection zones
   (sites closer to TSS) needs evaluation based on the data above.
""")

print("Analysis complete.")
