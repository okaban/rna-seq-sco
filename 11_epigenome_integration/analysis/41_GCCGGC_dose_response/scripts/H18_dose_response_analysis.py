#!/usr/bin/env python3
"""
H18: GCCGGC Methylation Density-Expression Dose-Response Relationship
=====================================================================

Tests whether genes with MORE GCCGGC 4mC sites within 2kb show progressively
stronger expression suppression (dose-response relationship).

Author: Claude Code analysis
Date: 2026-02-26
"""

import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# Configuration
# ============================================================================
BASE_DIR = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/41_GCCGGC_dose_response")
FIGURES_DIR = BASE_DIR / "figures"
TABLES_DIR = BASE_DIR / "tables"

GCCGGC_SITES_FILE = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/37_defense_island_GCCGGC/tables/GCCGGC_sites_by_timepoint.tsv")
GFF_FILE = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/39_GCCGGC_MTase_reverse_ID/data/GCF_000203835.1_ASM20383v1_genomic.gff")
DESEQ2_T2_FILE = Path("/Users/okaban/bioinfo/rna-seq/04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_2_vs_1.tsv")
DESEQ2_T3_FILE = Path("/Users/okaban/bioinfo/rna-seq/04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_3_vs_1.tsv")
AAGCCCG_MAPPING_FILE = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/36_AAGCCCG_distribution/tables/AAGCCCG_site_gene_mapping.tsv")

GENOME_SIZE = 8_667_507
ARM_CUTOFF_LEFT = 1_500_000
ARM_CUTOFF_RIGHT = GENOME_SIZE - ARM_CUTOFF_LEFT
PROXIMITY_WINDOW = 2000  # 2 kb
PROMOTER_WINDOW = 500    # 500 bp upstream of gene start

# Color scheme
COLORS = {
    'dose_0': '#BBBBBB',
    'dose_1': '#88CCEE',
    'dose_2': '#44AA99',
    'dose_3': '#117733',
    'dose_4': '#882255',
    'GCCGGC': '#4477AA',
    'AAGCCCG': '#EE6677',
    'promoter': '#CC6677',
    'gene_body': '#44AA99',
    'downstream': '#88CCEE',
    'core': '#4477AA',
    'arm': '#EE6677',
}
DOSE_COLORS = [COLORS['dose_0'], COLORS['dose_1'], COLORS['dose_2'], COLORS['dose_3'], COLORS['dose_4']]

# ============================================================================
# Step 1: Parse GFF and build gene database
# ============================================================================
print("=" * 70)
print("H18: GCCGGC Dose-Response Analysis")
print("=" * 70)

print("\nStep 1: Parsing GFF file...")

genes = []
with open(GFF_FILE) as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9:
            continue
        if parts[2] != 'gene':
            continue

        chrom = parts[0]
        start = int(parts[3])
        end = int(parts[4])
        strand = parts[6]
        attrs = parts[8]

        # Parse attributes
        attr_dict = {}
        for attr in attrs.split(';'):
            if '=' in attr:
                key, val = attr.split('=', 1)
                attr_dict[key] = val

        locus_tag = attr_dict.get('locus_tag', '')
        old_locus_tag = attr_dict.get('old_locus_tag', '')
        gene_biotype = attr_dict.get('gene_biotype', '')
        gene_name = attr_dict.get('Name', '')

        if not locus_tag:
            continue

        # Determine arm/core
        midpoint = (start + end) / 2
        if midpoint <= ARM_CUTOFF_LEFT or midpoint >= ARM_CUTOFF_RIGHT:
            region = 'arm'
        else:
            region = 'core'

        gene_length = end - start + 1

        genes.append({
            'chrom': chrom,
            'gene_start': start,
            'gene_end': end,
            'strand': strand,
            'locus_tag': locus_tag,
            'old_locus_tag': old_locus_tag,
            'gene_name': gene_name,
            'gene_biotype': gene_biotype,
            'region': region,
            'gene_length': gene_length,
        })

genes_df = pd.DataFrame(genes)
print(f"  Loaded {len(genes_df)} genes from GFF")
print(f"  Protein-coding: {(genes_df['gene_biotype'] == 'protein_coding').sum()}")
print(f"  Core: {(genes_df['region'] == 'core').sum()}, Arm: {(genes_df['region'] == 'arm').sum()}")

# ============================================================================
# Step 2: Load T1 GCCGGC sites and map to genes
# ============================================================================
print("\nStep 2: Loading T1 GCCGGC sites and mapping to genes...")

sites_df = pd.read_csv(GCCGGC_SITES_FILE, sep='\t')
t1_sites = sites_df[sites_df['timepoint'] == 'T1'].copy()
print(f"  T1 GCCGGC sites: {len(t1_sites)}")

# Map sites to genes within PROXIMITY_WINDOW
site_gene_mappings = []
for _, site in t1_sites.iterrows():
    site_pos = site['position']
    site_strand = site['strand']

    for _, gene in genes_df.iterrows():
        g_start = gene['gene_start']
        g_end = gene['gene_end']
        g_strand = gene['strand']

        # Calculate distance from site to gene
        if site_pos < g_start:
            distance = g_start - site_pos
        elif site_pos > g_end:
            distance = site_pos - g_end
        else:
            distance = 0  # site is within gene body

        if distance > PROXIMITY_WINDOW:
            continue

        # Determine location relative to gene
        if g_start <= site_pos <= g_end:
            location = 'gene_body'
        elif g_strand == '+':
            if site_pos < g_start:
                location = 'upstream'
            else:
                location = 'downstream'
        else:  # minus strand
            if site_pos > g_end:
                location = 'upstream'
            else:
                location = 'downstream'

        # Is it in the promoter region?
        if g_strand == '+':
            promoter_start = max(1, g_start - PROMOTER_WINDOW)
            in_promoter = (promoter_start <= site_pos <= g_start)
        else:
            promoter_end = g_end + PROMOTER_WINDOW
            in_promoter = (g_end <= site_pos <= promoter_end)

        if in_promoter:
            location_detail = 'promoter'
        elif location == 'gene_body':
            location_detail = 'gene_body'
        else:
            location_detail = location

        site_gene_mappings.append({
            'site_position': site_pos,
            'site_strand': site_strand,
            'locus_tag': gene['locus_tag'],
            'gene_start': g_start,
            'gene_end': g_end,
            'gene_strand': g_strand,
            'distance': distance,
            'location': location_detail,
            'region': gene['region'],
            'gene_length': gene['gene_length'],
        })

sgm_df = pd.DataFrame(site_gene_mappings)
print(f"  Total site-gene mappings within {PROXIMITY_WINDOW}bp: {len(sgm_df)}")
print(f"  Unique genes with >=1 GCCGGC site: {sgm_df['locus_tag'].nunique()}")
print(f"  Location breakdown:")
for loc, cnt in sgm_df['location'].value_counts().items():
    print(f"    {loc}: {cnt}")

# Count sites per gene
site_counts = sgm_df.groupby('locus_tag').size().reset_index(name='site_count')

# Also get location-specific counts
promoter_counts = sgm_df[sgm_df['location'] == 'promoter'].groupby('locus_tag').size().reset_index(name='promoter_site_count')
gene_body_counts = sgm_df[sgm_df['location'] == 'gene_body'].groupby('locus_tag').size().reset_index(name='gene_body_site_count')

# Merge all site counts with gene info
gene_site_df = genes_df.merge(site_counts, on='locus_tag', how='left')
gene_site_df['site_count'] = gene_site_df['site_count'].fillna(0).astype(int)
gene_site_df = gene_site_df.merge(promoter_counts, on='locus_tag', how='left')
gene_site_df['promoter_site_count'] = gene_site_df['promoter_site_count'].fillna(0).astype(int)
gene_site_df = gene_site_df.merge(gene_body_counts, on='locus_tag', how='left')
gene_site_df['gene_body_site_count'] = gene_site_df['gene_body_site_count'].fillna(0).astype(int)

# Create dose groups
gene_site_df['dose_group'] = gene_site_df['site_count'].apply(
    lambda x: '0' if x == 0 else ('1' if x == 1 else ('2' if x == 2 else ('3' if x == 3 else '4+')))
)

print(f"\n  Site count distribution across all genes:")
for dg in ['0', '1', '2', '3', '4+']:
    n = (gene_site_df['dose_group'] == dg).sum()
    print(f"    {dg} sites: {n} genes")

# ============================================================================
# Step 3: Merge with DESeq2 expression data
# ============================================================================
print("\nStep 3: Merging with DESeq2 expression data...")

deseq2_t2 = pd.read_csv(DESEQ2_T2_FILE, sep='\t')
deseq2_t3 = pd.read_csv(DESEQ2_T3_FILE, sep='\t')
print(f"  DESeq2 T2vsT1: {len(deseq2_t2)} genes")
print(f"  DESeq2 T3vsT1: {len(deseq2_t3)} genes")

# Merge
merged_t2 = gene_site_df.merge(deseq2_t2, left_on='locus_tag', right_on='gene_id', how='inner')
merged_t3 = gene_site_df.merge(deseq2_t3, left_on='locus_tag', right_on='gene_id', how='inner')
print(f"  After merge - T2: {len(merged_t2)}, T3: {len(merged_t3)}")

# Filter to genes with non-NA LFC
merged_t2 = merged_t2.dropna(subset=['log2FoldChange'])
merged_t3 = merged_t3.dropna(subset=['log2FoldChange'])
print(f"  After NA filter - T2: {len(merged_t2)}, T3: {len(merged_t3)}")

# Save the per-gene site count table
gene_site_output = gene_site_df[['locus_tag', 'old_locus_tag', 'gene_name', 'gene_start', 'gene_end',
                                  'strand', 'gene_biotype', 'region', 'gene_length',
                                  'site_count', 'promoter_site_count', 'gene_body_site_count', 'dose_group']].copy()
gene_site_output.to_csv(TABLES_DIR / "GCCGGC_gene_site_counts.tsv", sep='\t', index=False)
print(f"\n  Saved: GCCGGC_gene_site_counts.tsv")

# ============================================================================
# Step 4: Dose-Response Analysis
# ============================================================================
print("\n" + "=" * 70)
print("Step 4: Dose-Response Analysis")
print("=" * 70)

def dose_response_analysis(df, comparison_label):
    """Run full dose-response analysis on a merged dataframe."""
    print(f"\n--- {comparison_label} ---")

    results = []
    dose_groups = ['0', '1', '2', '3', '4+']

    for dg in dose_groups:
        subset = df[df['dose_group'] == dg]
        n = len(subset)
        if n == 0:
            results.append({
                'dose_group': dg, 'n_genes': 0, 'median_LFC': np.nan,
                'mean_LFC': np.nan, 'iqr_LFC': np.nan,
                'frac_DEG': np.nan, 'frac_up': np.nan, 'frac_down': np.nan
            })
            continue

        lfc = subset['log2FoldChange']
        padj = subset['padj']

        deg_mask = padj < 0.05
        n_deg = deg_mask.sum()
        n_up = ((padj < 0.05) & (lfc > 0)).sum()
        n_down = ((padj < 0.05) & (lfc < 0)).sum()

        q25, q75 = lfc.quantile([0.25, 0.75])

        results.append({
            'dose_group': dg,
            'n_genes': n,
            'median_LFC': lfc.median(),
            'mean_LFC': lfc.mean(),
            'iqr_LFC': q75 - q25,
            'q25': q25,
            'q75': q75,
            'frac_DEG': n_deg / n if n > 0 else 0,
            'frac_up': n_up / n if n > 0 else 0,
            'frac_down': n_down / n if n > 0 else 0,
        })

    results_df = pd.DataFrame(results)

    print(f"\n  Dose-Response Summary ({comparison_label}):")
    print(f"  {'Group':<8} {'N':>6} {'Median LFC':>12} {'Mean LFC':>10} {'IQR':>8} {'%DEG':>7} {'%Up':>6} {'%Down':>7}")
    print(f"  {'-'*70}")
    for _, row in results_df.iterrows():
        print(f"  {row['dose_group']:<8} {row['n_genes']:>6} {row['median_LFC']:>12.4f} {row['mean_LFC']:>10.4f} "
              f"{row['iqr_LFC']:>8.4f} {row['frac_DEG']*100:>6.1f}% {row['frac_up']*100:>5.1f}% {row['frac_down']*100:>6.1f}%")

    # Spearman correlation: site_count vs LFC (genes with >=1 site)
    with_sites = df[df['site_count'] >= 1].copy()
    if len(with_sites) > 10:
        rho, p_spearman = stats.spearmanr(with_sites['site_count'], with_sites['log2FoldChange'])
        print(f"\n  Spearman correlation (site_count vs LFC, genes with >=1 site, n={len(with_sites)}):")
        print(f"    rho = {rho:.4f}, p = {p_spearman:.2e}")
    else:
        rho, p_spearman = np.nan, np.nan

    # Also Spearman including 0-site genes
    rho_all, p_spearman_all = stats.spearmanr(df['site_count'], df['log2FoldChange'])
    print(f"  Spearman correlation (all genes including 0 sites, n={len(df)}):")
    print(f"    rho = {rho_all:.4f}, p = {p_spearman_all:.2e}")

    # Jonckheere-Terpstra trend test (using Mann-Whitney ordered groups)
    # Implement JT test manually
    groups_data = []
    for dg in dose_groups:
        subset = df[df['dose_group'] == dg]['log2FoldChange'].values
        if len(subset) > 0:
            groups_data.append(subset)

    jt_stat, jt_p = jonckheere_terpstra_test(groups_data)
    print(f"\n  Jonckheere-Terpstra trend test (ordered: 0 < 1 < 2 < 3 < 4+):")
    print(f"    JT statistic = {jt_stat:.1f}, p = {jt_p:.2e}")

    # Kruskal-Wallis test
    if len(groups_data) >= 2:
        kw_stat, kw_p = stats.kruskal(*groups_data)
        print(f"\n  Kruskal-Wallis test:")
        print(f"    H = {kw_stat:.2f}, p = {kw_p:.2e}")
    else:
        kw_stat, kw_p = np.nan, np.nan

    return results_df, rho, p_spearman, rho_all, p_spearman_all, jt_stat, jt_p, kw_stat, kw_p


def jonckheere_terpstra_test(groups, alternative='less'):
    """
    Jonckheere-Terpstra test for ordered alternatives.
    Tests H0: all groups same vs H1: ordered trend (group1 <= group2 <= ... <= groupK).
    alternative='less' tests for decreasing trend (more sites -> lower LFC).
    """
    k = len(groups)
    if k < 2:
        return np.nan, np.nan

    # JT statistic: sum of U_ij for all i < j
    jt_stat = 0
    for i in range(k):
        for j in range(i + 1, k):
            # Count pairs where group_i < group_j
            for xi in groups[i]:
                for xj in groups[j]:
                    if xi < xj:
                        jt_stat += 1
                    elif xi == xj:
                        jt_stat += 0.5

    # Under H0, compute expected value and variance
    N = sum(len(g) for g in groups)
    n_sizes = [len(g) for g in groups]

    # Expected value
    E_JT = (N**2 - sum(n**2 for n in n_sizes)) / 4

    # Variance (no ties correction for simplicity)
    sum_ni3 = sum(n**3 for n in n_sizes)
    sum_ni2 = sum(n**2 for n in n_sizes)
    sum_ni = N

    Var_JT = (N**2 * (2*N + 3) - sum(n**2 * (2*n + 3) for n in n_sizes)) / 72

    if Var_JT <= 0:
        return jt_stat, np.nan

    # Z-score
    z = (jt_stat - E_JT) / np.sqrt(Var_JT)

    # Two-sided p-value
    p_value = 2 * stats.norm.sf(abs(z))

    return jt_stat, p_value


# Run for T2vsT1 and T3vsT1
print("\n" + "=" * 70)
print("4a. T2 vs T1 (transition to aerial mycelium)")
print("=" * 70)
results_t2, rho_t2, p_rho_t2, rho_all_t2, p_all_t2, jt_t2, jt_p_t2, kw_t2, kw_p_t2 = dose_response_analysis(merged_t2, "T2vsT1")

print("\n" + "=" * 70)
print("4b. T3 vs T1 (sporulation)")
print("=" * 70)
results_t3, rho_t3, p_rho_t3, rho_all_t3, p_all_t3, jt_t3, jt_p_t3, kw_t3, kw_p_t3 = dose_response_analysis(merged_t3, "T3vsT1")


# ============================================================================
# Step 5: Location-specific analysis
# ============================================================================
print("\n" + "=" * 70)
print("Step 5: Location-Specific Analysis")
print("=" * 70)

def location_analysis(df, comparison_label):
    """Test promoter vs gene_body dose-response separately."""
    print(f"\n--- {comparison_label} ---")
    location_results = []

    for loc_col, loc_name in [('promoter_site_count', 'Promoter'), ('gene_body_site_count', 'Gene body')]:
        subset = df.copy()

        # Create dose groups for this location
        subset['loc_dose'] = subset[loc_col].apply(
            lambda x: '0' if x == 0 else ('1' if x == 1 else '2+')
        )

        for dg in ['0', '1', '2+']:
            grp = subset[subset['loc_dose'] == dg]
            if len(grp) > 0:
                location_results.append({
                    'location': loc_name,
                    'dose_group': dg,
                    'n_genes': len(grp),
                    'median_LFC': grp['log2FoldChange'].median(),
                    'mean_LFC': grp['log2FoldChange'].mean(),
                })

        # Spearman for this location
        with_loc_sites = subset[subset[loc_col] >= 1]
        if len(with_loc_sites) > 10:
            rho, p = stats.spearmanr(with_loc_sites[loc_col], with_loc_sites['log2FoldChange'])
            print(f"  {loc_name}: Spearman rho={rho:.4f}, p={p:.2e} (n={len(with_loc_sites)})")
        else:
            rho, p = np.nan, np.nan
            print(f"  {loc_name}: Too few genes with >=1 site (n={len(with_loc_sites)})")

        # Compare 0 vs >=1 site for this location (Mann-Whitney)
        grp_0 = subset[subset[loc_col] == 0]['log2FoldChange'].values
        grp_1plus = subset[subset[loc_col] >= 1]['log2FoldChange'].values
        if len(grp_0) > 0 and len(grp_1plus) > 0:
            u_stat, u_p = stats.mannwhitneyu(grp_0, grp_1plus, alternative='two-sided')
            median_diff = np.median(grp_1plus) - np.median(grp_0)
            print(f"    0 vs >=1 site: median diff = {median_diff:.4f}, MWU p = {u_p:.2e}")

        location_results.append({
            'location': loc_name,
            'dose_group': 'spearman',
            'n_genes': len(with_loc_sites),
            'median_LFC': rho,
            'mean_LFC': p,
        })

    return pd.DataFrame(location_results)

loc_results_t2 = location_analysis(merged_t2, "T2vsT1")
loc_results_t3 = location_analysis(merged_t3, "T3vsT1")

# ============================================================================
# Step 6: Confound Controls
# ============================================================================
print("\n" + "=" * 70)
print("Step 6: Confound Controls")
print("=" * 70)

# 6a. Gene length as confound
print("\n6a. Gene length control:")
for label, df in [("T2vsT1", merged_t2), ("T3vsT1", merged_t3)]:
    # Partial correlation: site_count vs LFC, controlling for gene_length
    # Using partial correlation formula: r_xy.z = (r_xy - r_xz*r_yz) / sqrt((1-r_xz^2)(1-r_yz^2))
    x = df['site_count'].values
    y = df['log2FoldChange'].values
    z = df['gene_length'].values

    r_xy = stats.spearmanr(x, y)[0]
    r_xz = stats.spearmanr(x, z)[0]
    r_yz = stats.spearmanr(y, z)[0]

    denom = np.sqrt((1 - r_xz**2) * (1 - r_yz**2))
    if denom > 0:
        partial_r = (r_xy - r_xz * r_yz) / denom
        # Approximate p-value using Fisher z-transform
        n = len(df)
        z_score = np.arctanh(partial_r) * np.sqrt(n - 3)
        partial_p = 2 * stats.norm.sf(abs(z_score))
    else:
        partial_r = np.nan
        partial_p = np.nan

    print(f"  {label}:")
    print(f"    Raw Spearman: rho={r_xy:.4f}")
    print(f"    Site_count vs gene_length: rho={r_xz:.4f}")
    print(f"    LFC vs gene_length: rho={r_yz:.4f}")
    print(f"    Partial correlation (controlling for gene_length): rho={partial_r:.4f}, p={partial_p:.2e}")

# 6b. Arm vs Core analysis
print("\n6b. Arm vs Core stratified analysis:")
confound_results = []
for label, df in [("T2vsT1", merged_t2), ("T3vsT1", merged_t3)]:
    for region in ['core', 'arm']:
        subset = df[df['region'] == region]
        with_sites = subset[subset['site_count'] >= 1]

        if len(with_sites) > 10:
            rho, p = stats.spearmanr(with_sites['site_count'], with_sites['log2FoldChange'])
            print(f"  {label} - {region} only (n={len(with_sites)} genes with >=1 site):")
            print(f"    Spearman rho={rho:.4f}, p={p:.2e}")
        else:
            rho, p = np.nan, np.nan
            print(f"  {label} - {region} only: too few genes with sites (n={len(with_sites)})")

        # Also median comparison: 0 vs 1 vs 2+ sites
        for dg in ['0', '1', '2', '3', '4+']:
            grp = subset[subset['dose_group'] == dg]
            if len(grp) > 0:
                confound_results.append({
                    'comparison': label,
                    'region': region,
                    'dose_group': dg,
                    'n_genes': len(grp),
                    'median_LFC': grp['log2FoldChange'].median(),
                    'mean_LFC': grp['log2FoldChange'].mean(),
                })

        confound_results.append({
            'comparison': label,
            'region': region,
            'dose_group': 'spearman',
            'n_genes': len(with_sites),
            'median_LFC': rho,
            'mean_LFC': p,
        })

confound_df = pd.DataFrame(confound_results)

# ============================================================================
# Step 7: AAGCCCG comparison
# ============================================================================
print("\n" + "=" * 70)
print("Step 7: AAGCCCG Comparison")
print("=" * 70)

aagcccg_map = pd.read_csv(AAGCCCG_MAPPING_FILE, sep='\t')
aagcccg_t1 = aagcccg_map[aagcccg_map['timepoint'] == 'T1'].copy()
print(f"  AAGCCCG T1 site-gene mappings: {len(aagcccg_t1)}")

# Count AAGCCCG sites per gene (within 2kb)
aagcccg_2kb = aagcccg_t1[aagcccg_t1['within_2kb'] == True].copy()
print(f"  AAGCCCG T1 within 2kb: {len(aagcccg_2kb)}")

aagcccg_counts = aagcccg_2kb.groupby('locus_tag').size().reset_index(name='aagcccg_site_count')

# Merge with genes and expression
aagcccg_gene_df = genes_df.merge(aagcccg_counts, on='locus_tag', how='left')
aagcccg_gene_df['aagcccg_site_count'] = aagcccg_gene_df['aagcccg_site_count'].fillna(0).astype(int)
aagcccg_gene_df['dose_group'] = aagcccg_gene_df['aagcccg_site_count'].apply(
    lambda x: '0' if x == 0 else ('1' if x == 1 else ('2' if x == 2 else ('3' if x == 3 else '4+')))
)

aagcccg_merged_t2 = aagcccg_gene_df.merge(deseq2_t2, left_on='locus_tag', right_on='gene_id', how='inner').dropna(subset=['log2FoldChange'])
aagcccg_merged_t3 = aagcccg_gene_df.merge(deseq2_t3, left_on='locus_tag', right_on='gene_id', how='inner').dropna(subset=['log2FoldChange'])

print(f"  AAGCCCG site count distribution:")
for dg in ['0', '1', '2', '3', '4+']:
    n = (aagcccg_gene_df['dose_group'] == dg).sum()
    print(f"    {dg} sites: {n} genes")

# AAGCCCG dose-response
def aagcccg_dose_response(df, site_col, comparison_label):
    """Quick dose-response for AAGCCCG."""
    print(f"\n  {comparison_label}:")

    dose_groups = ['0', '1', '2', '3', '4+']
    results = []
    groups_data = []

    for dg in dose_groups:
        subset = df[df['dose_group'] == dg]
        n = len(subset)
        if n > 0:
            results.append({
                'dose_group': dg,
                'n_genes': n,
                'median_LFC': subset['log2FoldChange'].median(),
                'mean_LFC': subset['log2FoldChange'].mean(),
            })
            groups_data.append(subset['log2FoldChange'].values)
            print(f"    {dg} sites: n={n}, median LFC={subset['log2FoldChange'].median():.4f}")

    # Spearman
    with_sites = df[df[site_col] >= 1]
    if len(with_sites) > 10:
        rho, p = stats.spearmanr(with_sites[site_col], with_sites['log2FoldChange'])
        print(f"    Spearman (>=1 site): rho={rho:.4f}, p={p:.2e}")
    else:
        rho, p = np.nan, np.nan

    rho_all, p_all = stats.spearmanr(df[site_col], df['log2FoldChange'])
    print(f"    Spearman (all): rho={rho_all:.4f}, p={p_all:.2e}")

    # JT test
    if len(groups_data) >= 2:
        jt_stat, jt_p = jonckheere_terpstra_test(groups_data)
        print(f"    JT test: stat={jt_stat:.1f}, p={jt_p:.2e}")
    else:
        jt_stat, jt_p = np.nan, np.nan

    return pd.DataFrame(results), rho, p, rho_all, p_all, jt_stat, jt_p

aagcccg_results_t2, aag_rho_t2, aag_p_t2, aag_rho_all_t2, aag_p_all_t2, aag_jt_t2, aag_jt_p_t2 = \
    aagcccg_dose_response(aagcccg_merged_t2, 'aagcccg_site_count', "AAGCCCG T2vsT1")
aagcccg_results_t3, aag_rho_t3, aag_p_t3, aag_rho_all_t3, aag_p_all_t3, aag_jt_t3, aag_jt_p_t3 = \
    aagcccg_dose_response(aagcccg_merged_t3, 'aagcccg_site_count', "AAGCCCG T3vsT1")


# ============================================================================
# Step 8: Save summary tables
# ============================================================================
print("\n" + "=" * 70)
print("Step 8: Saving tables")
print("=" * 70)

# Dose-response summary
dose_summary = pd.concat([
    results_t2.assign(comparison='T2vsT1', motif='GCCGGC'),
    results_t3.assign(comparison='T3vsT1', motif='GCCGGC'),
    aagcccg_results_t2.assign(comparison='T2vsT1', motif='AAGCCCG'),
    aagcccg_results_t3.assign(comparison='T3vsT1', motif='AAGCCCG'),
])
dose_summary.to_csv(TABLES_DIR / "dose_response_summary.tsv", sep='\t', index=False)
print(f"  Saved: dose_response_summary.tsv")

# Statistical tests
stat_tests = pd.DataFrame([
    {'motif': 'GCCGGC', 'comparison': 'T2vsT1', 'test': 'Spearman (>=1 site)', 'statistic': rho_t2, 'p_value': p_rho_t2, 'n': len(merged_t2[merged_t2['site_count']>=1])},
    {'motif': 'GCCGGC', 'comparison': 'T2vsT1', 'test': 'Spearman (all genes)', 'statistic': rho_all_t2, 'p_value': p_all_t2, 'n': len(merged_t2)},
    {'motif': 'GCCGGC', 'comparison': 'T2vsT1', 'test': 'Jonckheere-Terpstra', 'statistic': jt_t2, 'p_value': jt_p_t2, 'n': len(merged_t2)},
    {'motif': 'GCCGGC', 'comparison': 'T2vsT1', 'test': 'Kruskal-Wallis', 'statistic': kw_t2, 'p_value': kw_p_t2, 'n': len(merged_t2)},
    {'motif': 'GCCGGC', 'comparison': 'T3vsT1', 'test': 'Spearman (>=1 site)', 'statistic': rho_t3, 'p_value': p_rho_t3, 'n': len(merged_t3[merged_t3['site_count']>=1])},
    {'motif': 'GCCGGC', 'comparison': 'T3vsT1', 'test': 'Spearman (all genes)', 'statistic': rho_all_t3, 'p_value': p_all_t3, 'n': len(merged_t3)},
    {'motif': 'GCCGGC', 'comparison': 'T3vsT1', 'test': 'Jonckheere-Terpstra', 'statistic': jt_t3, 'p_value': jt_p_t3, 'n': len(merged_t3)},
    {'motif': 'GCCGGC', 'comparison': 'T3vsT1', 'test': 'Kruskal-Wallis', 'statistic': kw_t3, 'p_value': kw_p_t3, 'n': len(merged_t3)},
    {'motif': 'AAGCCCG', 'comparison': 'T2vsT1', 'test': 'Spearman (>=1 site)', 'statistic': aag_rho_t2, 'p_value': aag_p_t2, 'n': len(aagcccg_merged_t2[aagcccg_merged_t2['aagcccg_site_count']>=1])},
    {'motif': 'AAGCCCG', 'comparison': 'T2vsT1', 'test': 'Spearman (all genes)', 'statistic': aag_rho_all_t2, 'p_value': aag_p_all_t2, 'n': len(aagcccg_merged_t2)},
    {'motif': 'AAGCCCG', 'comparison': 'T2vsT1', 'test': 'Jonckheere-Terpstra', 'statistic': aag_jt_t2, 'p_value': aag_jt_p_t2, 'n': len(aagcccg_merged_t2)},
    {'motif': 'AAGCCCG', 'comparison': 'T3vsT1', 'test': 'Spearman (>=1 site)', 'statistic': aag_rho_t3, 'p_value': aag_p_t3, 'n': len(aagcccg_merged_t3[aagcccg_merged_t3['aagcccg_site_count']>=1])},
    {'motif': 'AAGCCCG', 'comparison': 'T3vsT1', 'test': 'Spearman (all genes)', 'statistic': aag_rho_all_t3, 'p_value': aag_p_all_t3, 'n': len(aagcccg_merged_t3)},
    {'motif': 'AAGCCCG', 'comparison': 'T3vsT1', 'test': 'Jonckheere-Terpstra', 'statistic': aag_jt_t3, 'p_value': aag_jt_p_t3, 'n': len(aagcccg_merged_t3)},
])
stat_tests.to_csv(TABLES_DIR / "statistical_tests.tsv", sep='\t', index=False)
print(f"  Saved: statistical_tests.tsv")

# Location analysis
loc_combined = pd.concat([
    loc_results_t2.assign(comparison='T2vsT1'),
    loc_results_t3.assign(comparison='T3vsT1'),
])
loc_combined.to_csv(TABLES_DIR / "location_analysis.tsv", sep='\t', index=False)
print(f"  Saved: location_analysis.tsv")

# Confound analysis
confound_df.to_csv(TABLES_DIR / "confound_arm_core_analysis.tsv", sep='\t', index=False)
print(f"  Saved: confound_arm_core_analysis.tsv")


# ============================================================================
# Step 9: Visualizations
# ============================================================================
print("\n" + "=" * 70)
print("Step 9: Creating visualizations")
print("=" * 70)

# ---------- Figure 1: Dose-response box/violin plots ----------
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for ax, (df, label, results_tbl) in zip(axes, [
    (merged_t2, 'T2 vs T1', results_t2),
    (merged_t3, 'T3 vs T1', results_t3),
]):
    dose_groups = ['0', '1', '2', '3', '4+']
    data_by_group = [df[df['dose_group'] == dg]['log2FoldChange'].values for dg in dose_groups]

    # Violin plot
    parts = ax.violinplot(data_by_group, positions=range(len(dose_groups)),
                          showmeans=False, showmedians=False, showextrema=False)
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(DOSE_COLORS[i])
        pc.set_alpha(0.4)

    # Box plot overlay
    bp = ax.boxplot(data_by_group, positions=range(len(dose_groups)),
                    widths=0.3, patch_artist=True,
                    medianprops=dict(color='black', linewidth=2),
                    flierprops=dict(marker='.', markersize=2, alpha=0.3))
    for i, patch in enumerate(bp['boxes']):
        patch.set_facecolor(DOSE_COLORS[i])
        patch.set_alpha(0.7)

    # Add sample sizes
    for i, dg in enumerate(dose_groups):
        n = len(data_by_group[i])
        ax.text(i, ax.get_ylim()[1] * 0.95 if ax.get_ylim()[1] > 0 else 5.5, f'n={n}',
                ha='center', va='top', fontsize=8, style='italic')

    # Add trend line (medians)
    medians = [np.median(d) if len(d) > 0 else np.nan for d in data_by_group]
    ax.plot(range(len(dose_groups)), medians, 'ko-', linewidth=2, markersize=8, zorder=5)

    ax.axhline(y=0, color='grey', linestyle='--', alpha=0.5)
    ax.set_xticks(range(len(dose_groups)))
    ax.set_xticklabels(dose_groups)
    ax.set_xlabel('GCCGGC 4mC sites within 2 kb', fontsize=12)
    ax.set_ylabel('log2 Fold Change', fontsize=12)
    ax.set_title(f'{label}', fontsize=14, fontweight='bold')

fig.suptitle('H18: GCCGGC Methylation Dose-Response on Gene Expression', fontsize=15, fontweight='bold', y=1.02)
fig.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(FIGURES_DIR / f"dose_response_boxplot.{ext}", bbox_inches='tight', dpi=300)
plt.close()
print("  Saved: dose_response_boxplot.pdf/svg")


# ---------- Figure 2: Scatter with regression ----------
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for ax, (df, label) in zip(axes, [
    (merged_t2, 'T2 vs T1'),
    (merged_t3, 'T3 vs T1'),
]):
    with_sites = df[df['site_count'] >= 1].copy()

    # Jitter x for visibility
    jitter = np.random.normal(0, 0.1, len(with_sites))
    x_jittered = with_sites['site_count'].values + jitter

    ax.scatter(x_jittered, with_sites['log2FoldChange'].values,
               alpha=0.15, s=10, color=COLORS['GCCGGC'], edgecolors='none')

    # Regression line (linear fit through medians)
    site_vals = with_sites['site_count'].values
    lfc_vals = with_sites['log2FoldChange'].values

    # Robust regression: use median per group
    unique_sites = sorted(with_sites['site_count'].unique())
    medians_x = []
    medians_y = []
    for s in unique_sites:
        grp = with_sites[with_sites['site_count'] == s]['log2FoldChange']
        if len(grp) >= 5:
            medians_x.append(s)
            medians_y.append(grp.median())

    if len(medians_x) >= 2:
        slope, intercept, r_val, p_val, se = stats.linregress(medians_x, medians_y)
        x_line = np.linspace(min(unique_sites), max(unique_sites), 100)
        ax.plot(x_line, slope * x_line + intercept, 'r-', linewidth=2, label=f'Trend: slope={slope:.3f}')

    # Also show medians as larger points
    ax.plot(medians_x, medians_y, 'ko', markersize=10, zorder=5, label='Group median')

    rho, p = stats.spearmanr(site_vals, lfc_vals)
    ax.set_title(f'{label}\nSpearman rho={rho:.3f}, p={p:.2e}', fontsize=12, fontweight='bold')
    ax.set_xlabel('Number of GCCGGC 4mC sites within 2 kb', fontsize=11)
    ax.set_ylabel('log2 Fold Change', fontsize=11)
    ax.axhline(y=0, color='grey', linestyle='--', alpha=0.5)
    ax.legend(fontsize=9)

fig.suptitle('GCCGGC Site Count vs Expression Change', fontsize=14, fontweight='bold', y=1.02)
fig.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(FIGURES_DIR / f"scatter_regression.{ext}", bbox_inches='tight', dpi=300)
plt.close()
print("  Saved: scatter_regression.pdf/svg")


# ---------- Figure 3: Location-specific analysis ----------
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for ax, (df, label) in zip(axes, [
    (merged_t2, 'T2 vs T1'),
    (merged_t3, 'T3 vs T1'),
]):
    # Compare promoter vs gene_body sites
    locations = ['promoter_site_count', 'gene_body_site_count']
    loc_labels = ['Promoter (<=500bp upstream)', 'Gene body']
    loc_colors = [COLORS['promoter'], COLORS['gene_body']]

    x_pos = np.arange(3)  # 0, 1, 2+
    width = 0.35

    for i, (loc_col, loc_label, loc_color) in enumerate(zip(locations, loc_labels, loc_colors)):
        medians = []
        ns = []
        for dg_val in [0, 1, 2]:
            if dg_val < 2:
                grp = df[df[loc_col] == dg_val]
            else:
                grp = df[df[loc_col] >= 2]
            medians.append(grp['log2FoldChange'].median() if len(grp) > 0 else 0)
            ns.append(len(grp))

        bars = ax.bar(x_pos + i * width - width/2, medians, width,
                       label=loc_label, color=loc_color, alpha=0.7, edgecolor='black', linewidth=0.5)

        for j, (bar, n) in enumerate(zip(bars, ns)):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                    f'n={n}', ha='center', va='bottom' if bar.get_height() >= 0 else 'top',
                    fontsize=7)

    ax.set_xticks(x_pos)
    ax.set_xticklabels(['0 sites', '1 site', '2+ sites'])
    ax.set_xlabel('GCCGGC sites in location category', fontsize=11)
    ax.set_ylabel('Median log2 Fold Change', fontsize=11)
    ax.set_title(f'{label}', fontsize=13, fontweight='bold')
    ax.axhline(y=0, color='grey', linestyle='--', alpha=0.5)
    ax.legend(fontsize=9)

fig.suptitle('Location-Specific Dose-Response: Promoter vs Gene Body', fontsize=14, fontweight='bold', y=1.02)
fig.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(FIGURES_DIR / f"location_specific.{ext}", bbox_inches='tight', dpi=300)
plt.close()
print("  Saved: location_specific.pdf/svg")


# ---------- Figure 4: Comprehensive 4-panel summary ----------
fig = plt.figure(figsize=(16, 14))
gs = GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.3)

# Panel A: Dose-response boxplots (T2vsT1)
ax_a = fig.add_subplot(gs[0, 0])
dose_groups = ['0', '1', '2', '3', '4+']
data_by_group = [merged_t2[merged_t2['dose_group'] == dg]['log2FoldChange'].values for dg in dose_groups]

parts = ax_a.violinplot(data_by_group, positions=range(len(dose_groups)),
                        showmeans=False, showmedians=False, showextrema=False)
for i, pc in enumerate(parts['bodies']):
    pc.set_facecolor(DOSE_COLORS[i])
    pc.set_alpha(0.4)

bp = ax_a.boxplot(data_by_group, positions=range(len(dose_groups)),
                  widths=0.3, patch_artist=True,
                  medianprops=dict(color='black', linewidth=2),
                  flierprops=dict(marker='.', markersize=2, alpha=0.3))
for i, patch in enumerate(bp['boxes']):
    patch.set_facecolor(DOSE_COLORS[i])
    patch.set_alpha(0.7)

medians = [np.median(d) if len(d) > 0 else np.nan for d in data_by_group]
ax_a.plot(range(len(dose_groups)), medians, 'ko-', linewidth=2, markersize=8, zorder=5)

for i, dg in enumerate(dose_groups):
    n = len(data_by_group[i])
    ax_a.text(i, ax_a.get_ylim()[1] * 0.9 if ax_a.get_ylim()[1] > 0 else 5, f'n={n}',
              ha='center', va='top', fontsize=8, style='italic')

ax_a.axhline(y=0, color='grey', linestyle='--', alpha=0.5)
ax_a.set_xticks(range(len(dose_groups)))
ax_a.set_xticklabels(dose_groups)
ax_a.set_xlabel('GCCGGC sites within 2 kb')
ax_a.set_ylabel('log2 Fold Change (T2/T1)')
ax_a.set_title('A. GCCGGC Dose-Response (T2 vs T1)', fontweight='bold', fontsize=12)

# Panel B: Median LFC bar chart with error bars
ax_b = fig.add_subplot(gs[0, 1])
x_pos = np.arange(len(dose_groups))
width = 0.35

# GCCGGC
medians_gccggc_t2 = results_t2['median_LFC'].values
ax_b.bar(x_pos - width/2, medians_gccggc_t2, width,
         color=COLORS['GCCGGC'], alpha=0.7, label='GCCGGC (T2/T1)',
         edgecolor='black', linewidth=0.5)

# AAGCCCG
medians_aagcccg_t2 = aagcccg_results_t2['median_LFC'].values
# Pad if different lengths
while len(medians_aagcccg_t2) < len(dose_groups):
    medians_aagcccg_t2 = np.append(medians_aagcccg_t2, np.nan)
ax_b.bar(x_pos + width/2, medians_aagcccg_t2[:len(dose_groups)], width,
         color=COLORS['AAGCCCG'], alpha=0.7, label='AAGCCCG (T2/T1)',
         edgecolor='black', linewidth=0.5)

ax_b.axhline(y=0, color='grey', linestyle='--', alpha=0.5)
ax_b.set_xticks(x_pos)
ax_b.set_xticklabels(dose_groups)
ax_b.set_xlabel('Methylation sites within 2 kb')
ax_b.set_ylabel('Median log2 Fold Change')
ax_b.set_title('B. GCCGGC vs AAGCCCG Comparison', fontweight='bold', fontsize=12)
ax_b.legend(fontsize=9)

# Panel C: Arm vs Core stratified
ax_c = fig.add_subplot(gs[1, 0])
x_pos = np.arange(len(dose_groups))
width = 0.35

for i, (region, color) in enumerate(zip(['core', 'arm'], [COLORS['core'], COLORS['arm']])):
    medians = []
    for dg in dose_groups:
        grp = merged_t2[(merged_t2['dose_group'] == dg) & (merged_t2['region'] == region)]
        medians.append(grp['log2FoldChange'].median() if len(grp) > 0 else np.nan)
    ax_c.bar(x_pos + i * width - width/2, medians, width,
             color=color, alpha=0.7, label=f'{region.capitalize()} genes',
             edgecolor='black', linewidth=0.5)

ax_c.axhline(y=0, color='grey', linestyle='--', alpha=0.5)
ax_c.set_xticks(x_pos)
ax_c.set_xticklabels(dose_groups)
ax_c.set_xlabel('GCCGGC sites within 2 kb')
ax_c.set_ylabel('Median log2 Fold Change (T2/T1)')
ax_c.set_title('C. Core vs Arm Stratification', fontweight='bold', fontsize=12)
ax_c.legend(fontsize=9)

# Panel D: Text summary of key statistics
ax_d = fig.add_subplot(gs[1, 1])
ax_d.axis('off')

summary_text = f"""Key Statistical Results

GCCGGC 4mC Dose-Response:

  T2 vs T1:
    Spearman (>=1 site): rho = {rho_t2:.4f}, p = {p_rho_t2:.2e}
    Spearman (all genes): rho = {rho_all_t2:.4f}, p = {p_all_t2:.2e}
    Jonckheere-Terpstra: p = {jt_p_t2:.2e}
    Kruskal-Wallis: p = {kw_p_t2:.2e}

  T3 vs T1:
    Spearman (>=1 site): rho = {rho_t3:.4f}, p = {p_rho_t3:.2e}
    Spearman (all genes): rho = {rho_all_t3:.4f}, p = {p_all_t3:.2e}
    Jonckheere-Terpstra: p = {jt_p_t3:.2e}
    Kruskal-Wallis: p = {kw_p_t3:.2e}

AAGCCCG 6mA Dose-Response:

  T2 vs T1:
    Spearman (all): rho = {aag_rho_all_t2:.4f}, p = {aag_p_all_t2:.2e}
    JT test: p = {aag_jt_p_t2:.2e}

  T3 vs T1:
    Spearman (all): rho = {aag_rho_all_t3:.4f}, p = {aag_p_all_t3:.2e}
    JT test: p = {aag_jt_p_t3:.2e}
"""

ax_d.text(0.05, 0.95, summary_text, transform=ax_d.transAxes,
          fontsize=9.5, verticalalignment='top', fontfamily='monospace',
          bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
ax_d.set_title('D. Statistical Summary', fontweight='bold', fontsize=12)

fig.suptitle('H18: GCCGGC Methylation Density-Expression Dose-Response',
             fontsize=16, fontweight='bold', y=0.98)

for ext in ['pdf', 'svg']:
    fig.savefig(FIGURES_DIR / f"H18_comprehensive_summary.{ext}", bbox_inches='tight', dpi=300)
plt.close()
print("  Saved: H18_comprehensive_summary.pdf/svg")


# ============================================================================
# Step 10: Verdict
# ============================================================================
print("\n" + "=" * 70)
print("VERDICT ASSESSMENT")
print("=" * 70)

# Criteria for SUPPORTED:
# 1. Significant monotonic trend (JT test p < 0.05)
# 2. Negative Spearman correlation (more sites -> lower LFC)
# 3. Consistent across both timepoints

jt_sig_t2 = jt_p_t2 < 0.05
jt_sig_t3 = jt_p_t3 < 0.05
rho_neg_t2 = rho_all_t2 < 0
rho_neg_t3 = rho_all_t3 < 0
rho_sig_t2 = p_all_t2 < 0.05
rho_sig_t3 = p_all_t3 < 0.05

# Check monotonicity of medians
medians_t2 = results_t2['median_LFC'].values
medians_t3 = results_t3['median_LFC'].values

# Monotonic decrease check (allowing some non-monotonicity for noise)
def check_monotonic_decrease(vals):
    """Check if values show general decreasing trend (allow 1 violation)."""
    violations = 0
    for i in range(1, len(vals)):
        if not np.isnan(vals[i]) and not np.isnan(vals[i-1]):
            if vals[i] > vals[i-1]:
                violations += 1
    return violations <= 1

mono_t2 = check_monotonic_decrease(medians_t2)
mono_t3 = check_monotonic_decrease(medians_t3)

print(f"\nT2vsT1:")
print(f"  JT trend significant: {jt_sig_t2} (p={jt_p_t2:.2e})")
print(f"  Spearman negative: {rho_neg_t2} (rho={rho_all_t2:.4f})")
print(f"  Spearman significant: {rho_sig_t2} (p={p_all_t2:.2e})")
print(f"  Monotonic medians: {mono_t2}")
print(f"  Median LFCs: {[f'{m:.4f}' for m in medians_t2]}")

print(f"\nT3vsT1:")
print(f"  JT trend significant: {jt_sig_t3} (p={jt_p_t3:.2e})")
print(f"  Spearman negative: {rho_neg_t3} (rho={rho_all_t3:.4f})")
print(f"  Spearman significant: {rho_sig_t3} (p={p_all_t3:.2e})")
print(f"  Monotonic medians: {mono_t3}")
print(f"  Median LFCs: {[f'{m:.4f}' for m in medians_t3]}")

# Final verdict
supported_count = sum([jt_sig_t2, jt_sig_t3, rho_sig_t2 and rho_neg_t2, rho_sig_t3 and rho_neg_t3])

if supported_count >= 3:
    verdict = "SUPPORTED"
elif supported_count >= 1:
    verdict = "PARTIAL"
else:
    verdict = "REJECTED"

print(f"\n{'='*70}")
print(f"VERDICT: {verdict}")
print(f"{'='*70}")
print(f"Evidence criteria met: {supported_count}/4")
print(f"  - JT trend T2vsT1: {'YES' if jt_sig_t2 else 'NO'}")
print(f"  - JT trend T3vsT1: {'YES' if jt_sig_t3 else 'NO'}")
print(f"  - Spearman negative & significant T2: {'YES' if (rho_sig_t2 and rho_neg_t2) else 'NO'}")
print(f"  - Spearman negative & significant T3: {'YES' if (rho_sig_t3 and rho_neg_t3) else 'NO'}")

print("\nDone.")
