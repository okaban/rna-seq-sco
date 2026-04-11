#!/usr/bin/env python3
"""
H25: Methylation Spatial Gradient Around Regulatory Gene TSS
============================================================
Hypothesis: Methylation site density decreases approaching TSS, with steeper
gradient and wider "protection zone" for regulatory genes vs non-regulatory genes.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch
from scipy import stats
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Configuration
# ============================================================
BASE = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/48_TSS_methylation_gradient")
FIG_DIR = BASE / "figures"
TBL_DIR = BASE / "tables"
CHROM_LEN = 8_667_507
ARM_LEFT = 1_500_000
ARM_RIGHT = 7_167_508
WINDOW = 5000       # ±5 kb around TSS
BIN_SIZE = 200      # 200 bp bins
N_BOOTSTRAP = 1000
N_PERM = 10000
SEED = 42
np.random.seed(SEED)

# ============================================================
# 1. Load and merge TSS data
# ============================================================
print("=" * 60)
print("STEP 1: Loading TSS data")
print("=" * 60)

# Comprehensive TSS table (already has merged TSS positions)
tss_df = pd.read_csv("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/18_tss_analyses/comprehensive_tss_table.csv")
print(f"Comprehensive TSS table: {len(tss_df)} genes")
print(f"Columns: {list(tss_df.columns)}")

# Gene annotation for biotype filtering
gene_ann = pd.read_csv("/Users/okaban/bioinfo/rna-seq/05_annotation/analysis/05_annotation_260128_v1/tables/gene_annotation_basic.tsv", sep='\t')
protein_coding = set(gene_ann[gene_ann['gene_biotype'] == 'protein_coding']['gene_id'].values)
print(f"Protein-coding genes: {len(protein_coding)}")

# Merge: keep only protein-coding genes with TSS
tss_df = tss_df[tss_df['gene_id'].isin(protein_coding)].copy()
print(f"Protein-coding genes with TSS: {len(tss_df)}")

# Filter to experimentally determined TSSs only (Jeong et al., 2016 dRNA-seq)
# GFF_annotation entries use the annotated gene start as a TSS proxy, which is
# not an experimentally validated TSS. Only Jeong2016_dRNA-seq entries have
# genuine TSS positions measured by differential RNA-seq.
tss_df = tss_df[tss_df['tss_source'] == 'Jeong2016_dRNA-seq'].copy()
print(f"Experimentally validated TSSs (Jeong2016 dRNA-seq only): {len(tss_df)}")

# ============================================================
# 2. Load regulatory gene list
# ============================================================
print("\n" + "=" * 60)
print("STEP 2: Loading regulatory genes")
print("=" * 60)

reg_df = pd.read_csv("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/29_genomewide_TF_screen/tables/all_regulatory_genes.tsv", sep='\t')
reg_genes = set(reg_df['locus_tag'].values)
print(f"Regulatory genes: {len(reg_genes)}")

# Subcategories
sigma_factors = set(reg_df[reg_df['tf_family'] == 'Sigma factor']['locus_tag'])
tcs_rr = set(reg_df[reg_df['tf_family'] == 'Response regulator']['locus_tag'])
sarp = set(reg_df[reg_df['tf_family'] == 'SARP']['locus_tag'])
# Other TFs: everything except sigma, TCS (SK+RR), SARP
sensor_kinase = set(reg_df[reg_df['tf_family'] == 'Sensor kinase']['locus_tag'])
other_tf = reg_genes - sigma_factors - tcs_rr - sarp - sensor_kinase

print(f"  Sigma factors: {len(sigma_factors)}")
print(f"  TCS response regulators: {len(tcs_rr)}")
print(f"  SARP-family: {len(sarp)}")
print(f"  Sensor kinases: {len(sensor_kinase)}")
print(f"  Other TFs: {len(other_tf)}")

# Tag genes in TSS table
tss_df['is_regulatory'] = tss_df['gene_id'].isin(reg_genes)
tss_df['reg_subtype'] = 'non_regulatory'
tss_df.loc[tss_df['gene_id'].isin(sigma_factors), 'reg_subtype'] = 'sigma_factor'
tss_df.loc[tss_df['gene_id'].isin(tcs_rr), 'reg_subtype'] = 'TCS_response_reg'
tss_df.loc[tss_df['gene_id'].isin(sarp), 'reg_subtype'] = 'SARP'
tss_df.loc[tss_df['gene_id'].isin(sensor_kinase), 'reg_subtype'] = 'sensor_kinase'
tss_df.loc[tss_df['gene_id'].isin(other_tf), 'reg_subtype'] = 'other_TF'

# Assign region (core/arm)
tss_df['region'] = 'core'
tss_df.loc[(tss_df['tss'] <= ARM_LEFT) | (tss_df['tss'] >= ARM_RIGHT), 'region'] = 'arm'

n_reg = tss_df['is_regulatory'].sum()
n_nonreg = (~tss_df['is_regulatory']).sum()
print(f"\nProtein-coding genes with TSS: regulatory={n_reg}, non-regulatory={n_nonreg}")
print(f"  Core: {(tss_df['region']=='core').sum()}, Arm: {(tss_df['region']=='arm').sum()}")

# ============================================================
# 3. Load methylation site positions
# ============================================================
print("\n" + "=" * 60)
print("STEP 3: Loading methylation sites")
print("=" * 60)

# All 4mC/6mA sites
all_sites = pd.read_csv("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/01_integration/high_confidence_sites_weighted.csv")
print(f"Total methylation site records: {len(all_sites)}")

# GCCGGC 4mC sites
gccggc_df = pd.read_csv("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/37_defense_island_GCCGGC/tables/GCCGGC_sites_by_timepoint.tsv", sep='\t')

# AAGCCCG 6mA sites - extract unique positions per timepoint
aagcccg_map = pd.read_csv("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/36_AAGCCCG_distribution/tables/AAGCCCG_site_gene_mapping.tsv", sep='\t')
aagcccg_unique = aagcccg_map[['position', 'timepoint']].drop_duplicates()

# Organize by motif type and timepoint
def get_positions(df, timepoint):
    """Extract unique positions for a given timepoint."""
    return np.sort(df[df['timepoint'] == timepoint]['position'].unique())

# GCCGGC sites
gccggc_T1 = get_positions(gccggc_df, 'T1')
gccggc_T2 = get_positions(gccggc_df, 'T2')

# AAGCCCG sites
aagcccg_T1 = get_positions(aagcccg_unique, 'T1')
aagcccg_T2 = get_positions(aagcccg_unique, 'T2')

# All 4mC sites
all_4mC_T1 = np.sort(all_sites[(all_sites['mod_type'] == '4mC') & (all_sites['timepoint'] == 'T1')]['position'].unique())
all_4mC_T2 = np.sort(all_sites[(all_sites['mod_type'] == '4mC') & (all_sites['timepoint'] == 'T2')]['position'].unique())

# All 6mA sites
all_6mA_T1 = np.sort(all_sites[(all_sites['mod_type'] == '6mA') & (all_sites['timepoint'] == 'T1')]['position'].unique())
all_6mA_T2 = np.sort(all_sites[(all_sites['mod_type'] == '6mA') & (all_sites['timepoint'] == 'T2')]['position'].unique())

# All sites combined
all_methyl_T1 = np.sort(all_sites[all_sites['timepoint'] == 'T1']['position'].unique())
all_methyl_T2 = np.sort(all_sites[all_sites['timepoint'] == 'T2']['position'].unique())

print(f"GCCGGC 4mC: T1={len(gccggc_T1)}, T2={len(gccggc_T2)}")
print(f"AAGCCCG 6mA: T1={len(aagcccg_T1)}, T2={len(aagcccg_T2)}")
print(f"All 4mC: T1={len(all_4mC_T1)}, T2={len(all_4mC_T2)}")
print(f"All 6mA: T1={len(all_6mA_T1)}, T2={len(all_6mA_T2)}")
print(f"All methylation: T1={len(all_methyl_T1)}, T2={len(all_methyl_T2)}")

# ============================================================
# 4. Core analysis: Calculate spatial profiles
# ============================================================
print("\n" + "=" * 60)
print("STEP 4: Calculating spatial profiles")
print("=" * 60)

bin_edges = np.arange(-WINDOW, WINDOW + BIN_SIZE, BIN_SIZE)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
n_bins = len(bin_centers)

def compute_tss_profile(tss_positions, tss_strands, methyl_positions, bin_edges):
    """
    Compute methylation density profile around TSS.

    For each TSS, find methylation sites within ±WINDOW and compute
    distances relative to TSS (negative = upstream, positive = downstream/gene body).
    For minus-strand genes, distances are flipped.

    Returns: per-gene counts array (n_genes x n_bins), and density array (per bin, per gene, per kb)
    """
    n_bins = len(bin_edges) - 1
    n_genes = len(tss_positions)
    per_gene_counts = np.zeros((n_genes, n_bins), dtype=int)

    # Use sorted methylation positions for efficient binary search
    methyl_sorted = np.sort(methyl_positions)

    for i, (tss, strand) in enumerate(zip(tss_positions, tss_strands)):
        # Find methylation sites within window
        left = np.searchsorted(methyl_sorted, tss - WINDOW, side='left')
        right = np.searchsorted(methyl_sorted, tss + WINDOW, side='right')
        nearby = methyl_sorted[left:right]

        if len(nearby) == 0:
            continue

        # Calculate distances (oriented: negative=upstream, positive=downstream)
        distances = nearby - tss
        if strand == '-':
            distances = -distances  # Flip for minus strand

        # Bin the distances
        bin_idx = np.digitize(distances, bin_edges) - 1
        valid = (bin_idx >= 0) & (bin_idx < n_bins)
        for idx in bin_idx[valid]:
            per_gene_counts[i, idx] += 1

    # Density: sites per kb per gene
    bin_width_kb = BIN_SIZE / 1000
    density = per_gene_counts.sum(axis=0) / (n_genes * bin_width_kb) if n_genes > 0 else np.zeros(n_bins)

    return per_gene_counts, density


def bootstrap_ci(per_gene_counts, n_genes, n_bootstrap=N_BOOTSTRAP, alpha=0.05):
    """Bootstrap 95% CI for density profile."""
    bin_width_kb = BIN_SIZE / 1000
    n_bins = per_gene_counts.shape[1]
    boot_densities = np.zeros((n_bootstrap, n_bins))

    for b in range(n_bootstrap):
        idx = np.random.randint(0, n_genes, size=n_genes)
        sampled = per_gene_counts[idx]
        boot_densities[b] = sampled.sum(axis=0) / (n_genes * bin_width_kb)

    ci_low = np.percentile(boot_densities, 100 * alpha / 2, axis=0)
    ci_high = np.percentile(boot_densities, 100 * (1 - alpha / 2), axis=0)
    return ci_low, ci_high


# Prepare TSS data
tss_positions = tss_df['tss'].values
tss_strands = tss_df['strand'].values
tss_is_reg = tss_df['is_regulatory'].values
tss_region = tss_df['region'].values
tss_subtype = tss_df['reg_subtype'].values

reg_mask = tss_is_reg
nonreg_mask = ~tss_is_reg

# Compute profiles for all methylation site types
results = {}
site_collections = {
    'All_methylation': {'T1': all_methyl_T1, 'T2': all_methyl_T2},
    'All_4mC': {'T1': all_4mC_T1, 'T2': all_4mC_T2},
    'All_6mA': {'T1': all_6mA_T1, 'T2': all_6mA_T2},
    'GCCGGC_4mC': {'T1': gccggc_T1, 'T2': gccggc_T2},
    'AAGCCCG_6mA': {'T1': aagcccg_T1, 'T2': aagcccg_T2},
}

gene_categories = {
    'all': np.ones(len(tss_df), dtype=bool),
    'regulatory': reg_mask,
    'non_regulatory': nonreg_mask,
}

print("Computing spatial profiles...")
for site_name, timepoints in site_collections.items():
    for tp, positions in timepoints.items():
        for cat_name, cat_mask in gene_categories.items():
            key = (site_name, tp, cat_name)
            tss_pos = tss_positions[cat_mask]
            tss_str = tss_strands[cat_mask]
            n_g = len(tss_pos)

            per_gene, density = compute_tss_profile(tss_pos, tss_str, positions, bin_edges)
            ci_low, ci_high = bootstrap_ci(per_gene, n_g)

            results[key] = {
                'density': density,
                'ci_low': ci_low,
                'ci_high': ci_high,
                'per_gene_counts': per_gene,
                'n_genes': n_g,
                'n_sites': len(positions),
            }
            print(f"  {site_name} {tp} {cat_name}: {n_g} genes, {len(positions)} sites, "
                  f"mean density={density.mean():.3f} sites/kb/gene")

# ============================================================
# 5. Core/arm stratification
# ============================================================
print("\n" + "=" * 60)
print("STEP 5: Core/arm stratification")
print("=" * 60)

region_categories = {}
for region_name in ['core', 'arm']:
    for reg_status in ['regulatory', 'non_regulatory']:
        if reg_status == 'regulatory':
            mask = reg_mask & (tss_region == region_name)
        else:
            mask = nonreg_mask & (tss_region == region_name)
        region_categories[f"{region_name}_{reg_status}"] = mask

# Compute for All methylation T1 only
region_results = {}
positions = all_methyl_T1
for cat_name, cat_mask in region_categories.items():
    tss_pos = tss_positions[cat_mask]
    tss_str = tss_strands[cat_mask]
    n_g = len(tss_pos)
    if n_g < 10:
        print(f"  Skipping {cat_name}: only {n_g} genes")
        continue
    per_gene, density = compute_tss_profile(tss_pos, tss_str, positions, bin_edges)
    ci_low, ci_high = bootstrap_ci(per_gene, n_g)
    region_results[cat_name] = {
        'density': density, 'ci_low': ci_low, 'ci_high': ci_high,
        'n_genes': n_g, 'per_gene_counts': per_gene,
    }
    print(f"  {cat_name}: {n_g} genes, mean density={density.mean():.3f}")

# ============================================================
# 6. Sub-category analysis
# ============================================================
print("\n" + "=" * 60)
print("STEP 6: Regulatory sub-category profiles")
print("=" * 60)

subtype_categories = {
    'sigma_factor': tss_subtype == 'sigma_factor',
    'TCS_response_reg': tss_subtype == 'TCS_response_reg',
    'SARP': tss_subtype == 'SARP',
    'sensor_kinase': tss_subtype == 'sensor_kinase',
    'other_TF': tss_subtype == 'other_TF',
}

subtype_results = {}
positions = all_methyl_T1
for cat_name, cat_mask in subtype_categories.items():
    tss_pos = tss_positions[cat_mask]
    tss_str = tss_strands[cat_mask]
    n_g = len(tss_pos)
    if n_g < 5:
        print(f"  Skipping {cat_name}: only {n_g} genes")
        continue
    per_gene, density = compute_tss_profile(tss_pos, tss_str, positions, bin_edges)
    ci_low, ci_high = bootstrap_ci(per_gene, n_g)
    subtype_results[cat_name] = {
        'density': density, 'ci_low': ci_low, 'ci_high': ci_high,
        'n_genes': n_g, 'per_gene_counts': per_gene,
    }
    print(f"  {cat_name}: {n_g} genes, mean density={density.mean():.3f}")

# ============================================================
# 7. Statistical tests
# ============================================================
print("\n" + "=" * 60)
print("STEP 7: Permutation tests (regulatory vs non-regulatory)")
print("=" * 60)

perm_results_list = []

for site_name in site_collections.keys():
    tp = 'T1'
    reg_key = (site_name, tp, 'regulatory')
    nonreg_key = (site_name, tp, 'non_regulatory')

    if reg_key not in results or nonreg_key not in results:
        continue

    reg_density = results[reg_key]['density']
    nonreg_density = results[nonreg_key]['density']

    # Combine per-gene counts for permutation
    reg_counts = results[reg_key]['per_gene_counts']
    nonreg_counts = results[nonreg_key]['per_gene_counts']
    n_reg = reg_counts.shape[0]
    n_nonreg = nonreg_counts.shape[0]
    combined = np.vstack([reg_counts, nonreg_counts])
    n_total = combined.shape[0]

    bin_width_kb = BIN_SIZE / 1000

    # Observed difference per bin
    obs_diff = reg_density - nonreg_density

    # Permutation test
    perm_diffs = np.zeros((N_PERM, n_bins))
    for p in range(N_PERM):
        perm_idx = np.random.permutation(n_total)
        perm_reg = combined[perm_idx[:n_reg]]
        perm_nonreg = combined[perm_idx[n_reg:]]
        perm_reg_density = perm_reg.sum(axis=0) / (n_reg * bin_width_kb)
        perm_nonreg_density = perm_nonreg.sum(axis=0) / (n_nonreg * bin_width_kb)
        perm_diffs[p] = perm_reg_density - perm_nonreg_density

    # Two-sided p-values
    p_values = np.zeros(n_bins)
    for b in range(n_bins):
        p_values[b] = np.mean(np.abs(perm_diffs[:, b]) >= np.abs(obs_diff[b]))
    p_values = np.maximum(p_values, 1.0 / N_PERM)  # floor at 1/N_PERM

    # BH correction
    from statsmodels.stats.multitest import multipletests
    _, p_adj, _, _ = multipletests(p_values, method='fdr_bh')

    for b in range(n_bins):
        perm_results_list.append({
            'site_type': site_name,
            'timepoint': tp,
            'bin_center': bin_centers[b],
            'reg_density': reg_density[b],
            'nonreg_density': nonreg_density[b],
            'density_diff': obs_diff[b],
            'ratio': reg_density[b] / nonreg_density[b] if nonreg_density[b] > 0 else np.nan,
            'p_value': p_values[b],
            'p_adj': p_adj[b],
            'significant': p_adj[b] < 0.05,
        })

    n_sig = np.sum(p_adj < 0.05)
    print(f"  {site_name} T1: {n_sig}/{n_bins} bins significant (FDR<0.05)")

perm_df = pd.DataFrame(perm_results_list)

# ============================================================
# 8. Protection zone metrics
# ============================================================
print("\n" + "=" * 60)
print("STEP 8: Protection zone analysis")
print("=" * 60)

protection_metrics = []

for site_name in site_collections.keys():
    tp = 'T1'
    reg_key = (site_name, tp, 'regulatory')
    nonreg_key = (site_name, tp, 'non_regulatory')

    reg_density = results[reg_key]['density']
    nonreg_density = results[nonreg_key]['density']

    # Find where regulatory density is below non-regulatory
    below = reg_density < nonreg_density

    # Find the TSS-centered protection zone
    # TSS is at bin index where bin_center = 0 (or closest)
    tss_bin = np.argmin(np.abs(bin_centers))

    # Scan outward from TSS to find protection zone boundaries
    # Upstream (negative direction)
    upstream_boundary = -WINDOW
    for b in range(tss_bin, -1, -1):
        if not below[b]:
            upstream_boundary = bin_centers[b]
            break

    # Downstream (positive direction)
    downstream_boundary = WINDOW
    for b in range(tss_bin, n_bins):
        if not below[b]:
            downstream_boundary = bin_centers[b]
            break

    # Find the deepest depletion (minimum ratio near TSS)
    # Focus on ±2kb around TSS
    near_tss = (np.abs(bin_centers) <= 2000)
    ratio_near = np.where(nonreg_density[near_tss] > 0,
                          reg_density[near_tss] / nonreg_density[near_tss],
                          np.nan)
    min_ratio = np.nanmin(ratio_near) if np.any(~np.isnan(ratio_near)) else np.nan
    min_ratio_pos = bin_centers[near_tss][np.nanargmin(ratio_near)] if np.any(~np.isnan(ratio_near)) else np.nan

    # AUC difference (regulatory - non-regulatory) in ±2kb window
    near_mask = np.abs(bin_centers) <= 2000
    auc_diff = np.sum((reg_density[near_mask] - nonreg_density[near_mask]) * BIN_SIZE / 1000)

    # Protection zone: contiguous region around TSS where reg < nonreg
    # Identify the longest contiguous stretch of below=True that includes the TSS bin
    pz_start = tss_bin
    pz_end = tss_bin
    for b in range(tss_bin, -1, -1):
        if below[b]:
            pz_start = b
        else:
            break
    for b in range(tss_bin, n_bins):
        if below[b]:
            pz_end = b
        else:
            break

    pz_width = bin_centers[pz_end] - bin_centers[pz_start] + BIN_SIZE if pz_start != pz_end else 0
    pz_center = (bin_centers[pz_start] + bin_centers[pz_end]) / 2

    metrics = {
        'site_type': site_name,
        'n_reg_genes': results[reg_key]['n_genes'],
        'n_nonreg_genes': results[nonreg_key]['n_genes'],
        'n_sites': results[reg_key]['n_sites'],
        'protection_zone_start_bp': bin_centers[pz_start] if pz_width > 0 else np.nan,
        'protection_zone_end_bp': bin_centers[pz_end] if pz_width > 0 else np.nan,
        'protection_zone_width_bp': pz_width,
        'protection_zone_center_bp': pz_center if pz_width > 0 else np.nan,
        'min_ratio_reg_vs_nonreg': min_ratio,
        'min_ratio_position_bp': min_ratio_pos,
        'AUC_diff_2kb': auc_diff,
        'mean_reg_density_1kb': reg_density[(np.abs(bin_centers) <= 500)].mean(),
        'mean_nonreg_density_1kb': nonreg_density[(np.abs(bin_centers) <= 500)].mean(),
        'mean_reg_density_5kb': reg_density.mean(),
        'mean_nonreg_density_5kb': nonreg_density.mean(),
    }
    protection_metrics.append(metrics)

    print(f"\n  {site_name}:")
    print(f"    Protection zone: {metrics['protection_zone_start_bp']:.0f} to {metrics['protection_zone_end_bp']:.0f} bp ({metrics['protection_zone_width_bp']:.0f} bp)")
    print(f"    Min ratio (reg/nonreg): {min_ratio:.3f} at {min_ratio_pos:.0f} bp")
    print(f"    AUC difference (±2kb): {auc_diff:.4f}")
    print(f"    Mean density ±500bp: reg={metrics['mean_reg_density_1kb']:.4f}, nonreg={metrics['mean_nonreg_density_1kb']:.4f}")

protection_df = pd.DataFrame(protection_metrics)

# ============================================================
# 9. Save tables
# ============================================================
print("\n" + "=" * 60)
print("STEP 9: Saving tables")
print("=" * 60)

# Spatial profile data
profile_rows = []
for (site_name, tp, cat_name), res in results.items():
    for b in range(n_bins):
        profile_rows.append({
            'site_type': site_name,
            'timepoint': tp,
            'gene_category': cat_name,
            'bin_center_bp': bin_centers[b],
            'density_per_kb_per_gene': res['density'][b],
            'ci_low': res['ci_low'][b],
            'ci_high': res['ci_high'][b],
            'n_genes': res['n_genes'],
            'n_sites': res['n_sites'],
        })
profile_df = pd.DataFrame(profile_rows)
profile_df.to_csv(TBL_DIR / "spatial_profile_data.tsv", sep='\t', index=False)
print(f"Saved spatial_profile_data.tsv: {len(profile_df)} rows")

# Protection zone metrics
protection_df.to_csv(TBL_DIR / "protection_zone_metrics.tsv", sep='\t', index=False)
print(f"Saved protection_zone_metrics.tsv: {len(protection_df)} rows")

# Permutation test results
perm_df.to_csv(TBL_DIR / "permutation_test_results.tsv", sep='\t', index=False)
print(f"Saved permutation_test_results.tsv: {len(perm_df)} rows")

# Subcategory profiles
sub_rows = []
for cat_name, res in subtype_results.items():
    for b in range(n_bins):
        sub_rows.append({
            'reg_subtype': cat_name,
            'bin_center_bp': bin_centers[b],
            'density_per_kb_per_gene': res['density'][b],
            'ci_low': res['ci_low'][b],
            'ci_high': res['ci_high'][b],
            'n_genes': res['n_genes'],
        })
# Add non-regulatory as reference
nonreg_res = results[('All_methylation', 'T1', 'non_regulatory')]
for b in range(n_bins):
    sub_rows.append({
        'reg_subtype': 'non_regulatory',
        'bin_center_bp': bin_centers[b],
        'density_per_kb_per_gene': nonreg_res['density'][b],
        'ci_low': nonreg_res['ci_low'][b],
        'ci_high': nonreg_res['ci_high'][b],
        'n_genes': nonreg_res['n_genes'],
    })
sub_df = pd.DataFrame(sub_rows)
sub_df.to_csv(TBL_DIR / "subcategory_profiles.tsv", sep='\t', index=False)
print(f"Saved subcategory_profiles.tsv: {len(sub_df)} rows")

# Region profiles
reg_rows = []
for cat_name, res in region_results.items():
    for b in range(n_bins):
        reg_rows.append({
            'region_category': cat_name,
            'bin_center_bp': bin_centers[b],
            'density_per_kb_per_gene': res['density'][b],
            'ci_low': res['ci_low'][b],
            'ci_high': res['ci_high'][b],
            'n_genes': res['n_genes'],
        })
region_profile_df = pd.DataFrame(reg_rows)
region_profile_df.to_csv(TBL_DIR / "region_stratified_profiles.tsv", sep='\t', index=False)
print(f"Saved region_stratified_profiles.tsv: {len(region_profile_df)} rows")

# ============================================================
# 10. Visualization
# ============================================================
print("\n" + "=" * 60)
print("STEP 10: Creating figures")
print("=" * 60)

# Color scheme
COL_REG = '#D32F2F'       # Red for regulatory
COL_NONREG = '#1565C0'    # Blue for non-regulatory
COL_ALL = '#616161'        # Gray for all genes
COL_REG_LIGHT = '#EF9A9A'
COL_NONREG_LIGHT = '#90CAF9'

# Subtype colors
SUBTYPE_COLORS = {
    'sigma_factor': '#E65100',
    'TCS_response_reg': '#AD1457',
    'SARP': '#6A1B9A',
    'sensor_kinase': '#00695C',
    'other_TF': '#F57F17',
    'non_regulatory': '#1565C0',
}

SUBTYPE_LABELS = {
    'sigma_factor': 'Sigma factors',
    'TCS_response_reg': 'TCS response regulators',
    'SARP': 'SARP family',
    'sensor_kinase': 'Sensor kinases',
    'other_TF': 'Other TFs',
    'non_regulatory': 'Non-regulatory',
}

# ---------- Panel A: Main TSS-centered profile ----------
def plot_profile(ax, site_name, tp, title, show_legend=True):
    """Plot regulatory vs non-regulatory density profile."""
    reg_res = results[(site_name, tp, 'regulatory')]
    nonreg_res = results[(site_name, tp, 'non_regulatory')]
    all_res = results[(site_name, tp, 'all')]

    # CI ribbons
    ax.fill_between(bin_centers, reg_res['ci_low'], reg_res['ci_high'],
                    color=COL_REG, alpha=0.15, linewidth=0)
    ax.fill_between(bin_centers, nonreg_res['ci_low'], nonreg_res['ci_high'],
                    color=COL_NONREG, alpha=0.15, linewidth=0)

    # Lines
    ax.plot(bin_centers, nonreg_res['density'], color=COL_NONREG, linewidth=1.5,
            label=f"Non-regulatory (n={nonreg_res['n_genes']})")
    ax.plot(bin_centers, reg_res['density'], color=COL_REG, linewidth=1.5,
            label=f"Regulatory (n={reg_res['n_genes']})")

    # TSS line
    ax.axvline(0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)

    ax.set_xlabel('Distance from TSS (bp)', fontsize=10)
    ax.set_ylabel('Methylation density\n(sites/kb/gene)', fontsize=10)
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.tick_params(labelsize=9)

    if show_legend:
        ax.legend(fontsize=8, loc='upper right', framealpha=0.9)

    # Add upstream/downstream annotations
    ax.text(-4500, ax.get_ylim()[1] * 0.95, 'Upstream', fontsize=8,
            ha='center', va='top', color='gray', style='italic')
    ax.text(4500, ax.get_ylim()[1] * 0.95, 'Gene body', fontsize=8,
            ha='center', va='top', color='gray', style='italic')


# Figure 1: Main profile (All methylation, T1)
fig1, ax1 = plt.subplots(figsize=(8, 5))
plot_profile(ax1, 'All_methylation', 'T1',
             'Methylation density around TSS (All sites, T1)')
fig1.tight_layout()
fig1.savefig(FIG_DIR / "TSS_methylation_profile.pdf", dpi=300, bbox_inches='tight')
fig1.savefig(FIG_DIR / "TSS_methylation_profile.svg", dpi=300, bbox_inches='tight')
print("Saved TSS_methylation_profile.pdf/svg")

# ---------- Panel B: Motif-specific profiles ----------
fig2, axes2 = plt.subplots(2, 2, figsize=(12, 9))
motifs_to_plot = [
    ('All_4mC', 'T1', 'All 4mC sites (T1)'),
    ('All_6mA', 'T1', 'All 6mA sites (T1)'),
    ('GCCGGC_4mC', 'T1', 'GCCGGC 4mC sites (T1)'),
    ('AAGCCCG_6mA', 'T1', 'AAGCCCG 6mA sites (T1)'),
]
for ax, (sn, tp, title) in zip(axes2.flat, motifs_to_plot):
    plot_profile(ax, sn, tp, title)

fig2.suptitle('Motif-specific methylation profiles around TSS', fontsize=13, fontweight='bold', y=1.01)
fig2.tight_layout()
fig2.savefig(FIG_DIR / "motif_specific_profiles.pdf", dpi=300, bbox_inches='tight')
fig2.savefig(FIG_DIR / "motif_specific_profiles.svg", dpi=300, bbox_inches='tight')
print("Saved motif_specific_profiles.pdf/svg")

# ---------- Panel C: Heatmap of per-bin ratio ----------
fig3, ax3 = plt.subplots(figsize=(10, 4))

# Build ratio matrix (site_type x bin)
site_names_ordered = ['All_methylation', 'All_4mC', 'All_6mA', 'GCCGGC_4mC', 'AAGCCCG_6mA']
site_labels = ['All methylation', 'All 4mC', 'All 6mA', 'GCCGGC 4mC', 'AAGCCCG 6mA']
ratio_matrix = np.zeros((len(site_names_ordered), n_bins))

for i, sn in enumerate(site_names_ordered):
    reg_d = results[(sn, 'T1', 'regulatory')]['density']
    nonreg_d = results[(sn, 'T1', 'non_regulatory')]['density']
    ratio_matrix[i] = np.where(nonreg_d > 0, reg_d / nonreg_d, np.nan)

# Clip ratio for visualization
ratio_clipped = np.clip(ratio_matrix, 0.5, 1.5)

im = ax3.imshow(ratio_clipped, aspect='auto', cmap='RdBu_r',
                vmin=0.5, vmax=1.5,
                extent=[bin_centers[0], bin_centers[-1], len(site_names_ordered) - 0.5, -0.5])
ax3.set_yticks(range(len(site_labels)))
ax3.set_yticklabels(site_labels, fontsize=10)
ax3.set_xlabel('Distance from TSS (bp)', fontsize=10)
ax3.set_title('Regulatory / Non-regulatory density ratio', fontsize=11, fontweight='bold')
ax3.axvline(0, color='black', linestyle='--', linewidth=0.8, alpha=0.7)

# Add significance markers
for i, sn in enumerate(site_names_ordered):
    sub = perm_df[(perm_df['site_type'] == sn) & (perm_df['timepoint'] == 'T1')]
    sig_bins = sub[sub['significant']]['bin_center'].values
    for sb in sig_bins:
        ax3.plot(sb, i, 'k.', markersize=2, alpha=0.5)

cbar = plt.colorbar(im, ax=ax3, label='Density ratio (reg/nonreg)', shrink=0.8)
cbar.ax.axhline(1.0, color='black', linewidth=0.5)

fig3.tight_layout()
fig3.savefig(FIG_DIR / "protection_zone_heatmap.pdf", dpi=300, bbox_inches='tight')
fig3.savefig(FIG_DIR / "protection_zone_heatmap.svg", dpi=300, bbox_inches='tight')
print("Saved protection_zone_heatmap.pdf/svg")

# ---------- Temporal comparison (T1 vs T2) ----------
fig_temp, axes_temp = plt.subplots(1, 2, figsize=(12, 5))

for ax, (cat, color, label) in zip(axes_temp,
    [('regulatory', COL_REG, 'Regulatory'), ('non_regulatory', COL_NONREG, 'Non-regulatory')]):
    t1_res = results[('All_methylation', 'T1', cat)]
    t2_res = results[('All_methylation', 'T2', cat)]

    ax.fill_between(bin_centers, t1_res['ci_low'], t1_res['ci_high'], alpha=0.15, color=color)
    ax.plot(bin_centers, t1_res['density'], color=color, linewidth=1.5, label='T1 (24h)')
    ax.plot(bin_centers, t2_res['density'], color=color, linewidth=1.5, linestyle='--',
            alpha=0.7, label='T2 (36h)')
    ax.axvline(0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.set_xlabel('Distance from TSS (bp)', fontsize=10)
    ax.set_ylabel('Methylation density\n(sites/kb/gene)', fontsize=10)
    ax.set_title(f'{label} genes', fontsize=11, fontweight='bold')
    ax.legend(fontsize=9)

fig_temp.suptitle('Temporal comparison: T1 vs T2', fontsize=13, fontweight='bold', y=1.01)
fig_temp.tight_layout()
fig_temp.savefig(FIG_DIR / "temporal_comparison.pdf", dpi=300, bbox_inches='tight')
fig_temp.savefig(FIG_DIR / "temporal_comparison.svg", dpi=300, bbox_inches='tight')
print("Saved temporal_comparison.pdf/svg")

# ---------- Core/arm stratification ----------
fig_region, axes_region = plt.subplots(1, 2, figsize=(12, 5))

for ax, region in zip(axes_region, ['core', 'arm']):
    reg_key = f"{region}_regulatory"
    nonreg_key = f"{region}_non_regulatory"
    if reg_key in region_results and nonreg_key in region_results:
        r_reg = region_results[reg_key]
        r_nonreg = region_results[nonreg_key]

        ax.fill_between(bin_centers, r_reg['ci_low'], r_reg['ci_high'],
                        color=COL_REG, alpha=0.15)
        ax.fill_between(bin_centers, r_nonreg['ci_low'], r_nonreg['ci_high'],
                        color=COL_NONREG, alpha=0.15)
        ax.plot(bin_centers, r_nonreg['density'], color=COL_NONREG, linewidth=1.5,
                label=f"Non-regulatory (n={r_nonreg['n_genes']})")
        ax.plot(bin_centers, r_reg['density'], color=COL_REG, linewidth=1.5,
                label=f"Regulatory (n={r_reg['n_genes']})")

    ax.axvline(0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.set_xlabel('Distance from TSS (bp)', fontsize=10)
    ax.set_ylabel('Methylation density\n(sites/kb/gene)', fontsize=10)
    ax.set_title(f'{region.capitalize()} region', fontsize=11, fontweight='bold')
    ax.legend(fontsize=8)

fig_region.suptitle('Core vs Arm: TSS methylation profiles (All sites, T1)', fontsize=13, fontweight='bold', y=1.01)
fig_region.tight_layout()
fig_region.savefig(FIG_DIR / "core_arm_profiles.pdf", dpi=300, bbox_inches='tight')
fig_region.savefig(FIG_DIR / "core_arm_profiles.svg", dpi=300, bbox_inches='tight')
print("Saved core_arm_profiles.pdf/svg")

# ---------- Sub-category profiles ----------
fig_sub, ax_sub = plt.subplots(figsize=(9, 6))

for cat_name in ['non_regulatory', 'sigma_factor', 'TCS_response_reg', 'SARP', 'sensor_kinase', 'other_TF']:
    if cat_name == 'non_regulatory':
        res = results[('All_methylation', 'T1', 'non_regulatory')]
    elif cat_name in subtype_results:
        res = subtype_results[cat_name]
    else:
        continue

    color = SUBTYPE_COLORS.get(cat_name, 'gray')
    label = SUBTYPE_LABELS.get(cat_name, cat_name)
    lw = 2.0 if cat_name == 'non_regulatory' else 1.3
    ls = '-' if cat_name == 'non_regulatory' else '-'
    alpha = 1.0 if cat_name == 'non_regulatory' else 0.8

    ax_sub.plot(bin_centers, res['density'], color=color, linewidth=lw,
                label=f"{label} (n={res['n_genes']})", alpha=alpha)

ax_sub.axvline(0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
ax_sub.set_xlabel('Distance from TSS (bp)', fontsize=10)
ax_sub.set_ylabel('Methylation density (sites/kb/gene)', fontsize=10)
ax_sub.set_title('Regulatory sub-type methylation profiles (All sites, T1)', fontsize=11, fontweight='bold')
ax_sub.legend(fontsize=8, loc='upper right')

fig_sub.tight_layout()
fig_sub.savefig(FIG_DIR / "subcategory_profiles.pdf", dpi=300, bbox_inches='tight')
fig_sub.savefig(FIG_DIR / "subcategory_profiles.svg", dpi=300, bbox_inches='tight')
print("Saved subcategory_profiles.pdf/svg")

# ---------- Comprehensive multi-panel figure ----------
print("\nCreating comprehensive multi-panel figure...")
fig_main = plt.figure(figsize=(16, 18))
gs = gridspec.GridSpec(4, 2, figure=fig_main, hspace=0.35, wspace=0.3)

# Panel A: Main profile (All methylation, T1) - large
ax_A = fig_main.add_subplot(gs[0, :])
reg_res = results[('All_methylation', 'T1', 'regulatory')]
nonreg_res = results[('All_methylation', 'T1', 'non_regulatory')]
ax_A.fill_between(bin_centers, reg_res['ci_low'], reg_res['ci_high'],
                  color=COL_REG, alpha=0.15, linewidth=0)
ax_A.fill_between(bin_centers, nonreg_res['ci_low'], nonreg_res['ci_high'],
                  color=COL_NONREG, alpha=0.15, linewidth=0)
ax_A.plot(bin_centers, nonreg_res['density'], color=COL_NONREG, linewidth=2,
          label=f"Non-regulatory (n={nonreg_res['n_genes']})")
ax_A.plot(bin_centers, reg_res['density'], color=COL_REG, linewidth=2,
          label=f"Regulatory (n={reg_res['n_genes']})")
ax_A.axvline(0, color='black', linestyle='--', linewidth=1, alpha=0.5)
ax_A.set_xlabel('Distance from TSS (bp)', fontsize=11)
ax_A.set_ylabel('Methylation density\n(sites/kb/gene)', fontsize=11)
ax_A.set_title('A. TSS-centered methylation density profile (All sites, T1)',
               fontsize=12, fontweight='bold', loc='left')
ax_A.legend(fontsize=10, loc='upper right', framealpha=0.9)
ax_A.text(-4500, ax_A.get_ylim()[1] * 0.95, 'Upstream', fontsize=9,
          ha='center', va='top', color='gray', style='italic')
ax_A.text(4500, ax_A.get_ylim()[1] * 0.95, 'Gene body', fontsize=9,
          ha='center', va='top', color='gray', style='italic')

# Panels B1-B4: Motif-specific profiles
for idx, (sn, tp, title) in enumerate(motifs_to_plot):
    row = 1 + idx // 2
    col = idx % 2
    ax_B = fig_main.add_subplot(gs[row, col])

    reg_res = results[(sn, tp, 'regulatory')]
    nonreg_res = results[(sn, tp, 'non_regulatory')]

    ax_B.fill_between(bin_centers, reg_res['ci_low'], reg_res['ci_high'],
                      color=COL_REG, alpha=0.15, linewidth=0)
    ax_B.fill_between(bin_centers, nonreg_res['ci_low'], nonreg_res['ci_high'],
                      color=COL_NONREG, alpha=0.15, linewidth=0)
    ax_B.plot(bin_centers, nonreg_res['density'], color=COL_NONREG, linewidth=1.5,
              label=f"Non-reg (n={nonreg_res['n_genes']})")
    ax_B.plot(bin_centers, reg_res['density'], color=COL_REG, linewidth=1.5,
              label=f"Regulatory (n={reg_res['n_genes']})")
    ax_B.axvline(0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    ax_B.set_xlabel('Distance from TSS (bp)', fontsize=10)
    ax_B.set_ylabel('Density (sites/kb/gene)', fontsize=10)
    panel_label = chr(ord('B') + idx)
    ax_B.set_title(f'{panel_label}. {title}', fontsize=11, fontweight='bold', loc='left')
    ax_B.legend(fontsize=8, loc='upper right')

# Panel F: Heatmap
ax_F = fig_main.add_subplot(gs[3, 0])
ratio_clipped2 = np.clip(ratio_matrix, 0.5, 1.5)
im2 = ax_F.imshow(ratio_clipped2, aspect='auto', cmap='RdBu_r',
                   vmin=0.5, vmax=1.5,
                   extent=[bin_centers[0], bin_centers[-1], len(site_names_ordered) - 0.5, -0.5])
ax_F.set_yticks(range(len(site_labels)))
ax_F.set_yticklabels(site_labels, fontsize=9)
ax_F.set_xlabel('Distance from TSS (bp)', fontsize=10)
ax_F.set_title('F. Density ratio heatmap (reg/non-reg)', fontsize=11, fontweight='bold', loc='left')
ax_F.axvline(0, color='black', linestyle='--', linewidth=0.8, alpha=0.7)
plt.colorbar(im2, ax=ax_F, label='Ratio', shrink=0.7)

# Panel G: Protection zone width comparison
ax_G = fig_main.add_subplot(gs[3, 1])
pz_widths = protection_df['protection_zone_width_bp'].values
pz_labels = [sl for sl in site_labels]
colors_bar = ['#424242', '#1565C0', '#D32F2F', '#00695C', '#F57F17']
bars = ax_G.barh(range(len(pz_labels)), pz_widths, color=colors_bar, edgecolor='black', linewidth=0.5)
ax_G.set_yticks(range(len(pz_labels)))
ax_G.set_yticklabels(pz_labels, fontsize=9)
ax_G.set_xlabel('Protection zone width (bp)', fontsize=10)
ax_G.set_title('G. Protection zone width by motif', fontsize=11, fontweight='bold', loc='left')
for i, w in enumerate(pz_widths):
    ax_G.text(w + 20, i, f'{w:.0f} bp', va='center', fontsize=9)

fig_main.savefig(FIG_DIR / "H25_comprehensive_summary.pdf", dpi=300, bbox_inches='tight')
fig_main.savefig(FIG_DIR / "H25_comprehensive_summary.svg", dpi=300, bbox_inches='tight')
print("Saved H25_comprehensive_summary.pdf/svg")

plt.close('all')

# ============================================================
# 11. Summary statistics for report
# ============================================================
print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)

# Key findings
for _, row in protection_df.iterrows():
    site = row['site_type']
    pzw = row['protection_zone_width_bp']
    mr = row['min_ratio_reg_vs_nonreg']
    mrp = row['min_ratio_position_bp']
    reg_d = row['mean_reg_density_1kb']
    nonreg_d = row['mean_nonreg_density_1kb']
    print(f"\n{site}:")
    print(f"  Protection zone width: {pzw:.0f} bp")
    print(f"  Min ratio (reg/nonreg): {mr:.3f} at position {mrp:.0f} bp")
    print(f"  Mean density ±500bp: regulatory={reg_d:.4f}, non-regulatory={nonreg_d:.4f}")
    if nonreg_d > 0:
        print(f"  Density reduction at TSS: {(1-reg_d/nonreg_d)*100:.1f}%")

# Check for universal TSS effect
all_res_T1 = results[('All_methylation', 'T1', 'all')]
tss_bin_idx = np.argmin(np.abs(bin_centers))
flanking_bins = (np.abs(bin_centers) > 3000)
tss_bins = (np.abs(bin_centers) <= 500)
flank_density = all_res_T1['density'][flanking_bins].mean()
tss_density = all_res_T1['density'][tss_bins].mean()
print(f"\nUniversal TSS effect (all genes):")
print(f"  Flanking density (>3kb): {flank_density:.4f}")
print(f"  TSS density (±500bp): {tss_density:.4f}")
if flank_density > 0:
    print(f"  TSS depletion: {(1-tss_density/flank_density)*100:.1f}%")

# Temporal comparison
for cat in ['regulatory', 'non_regulatory']:
    t1 = results[('All_methylation', 'T1', cat)]['density'][tss_bins].mean()
    t2 = results[('All_methylation', 'T2', cat)]['density'][tss_bins].mean()
    print(f"\nTemporal change in TSS density ({cat}):")
    print(f"  T1: {t1:.4f}, T2: {t2:.4f}")
    if t1 > 0:
        print(f"  Change: {(t2/t1 - 1)*100:.1f}%")

print("\n" + "=" * 60)
print("H25 analysis complete!")
print("=" * 60)
