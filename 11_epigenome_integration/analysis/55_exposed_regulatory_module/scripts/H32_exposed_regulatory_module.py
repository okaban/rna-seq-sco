#!/usr/bin/env python3
"""
H32: Exposed Regulatory Module Analysis
========================================
Tests whether the 57 "exposed" transcription factors in S. coelicolor M145
form an interconnected self-regulatory methylation-responsive module, or are
independently scattered across the genome.

Steps:
  1. Load data
  2. Genomic clustering analysis (permutation test, cluster detection)
  3. Co-expression analysis (Spearman correlations, hierarchical clustering)
  4. Coordination type coherence
  5. TCS pair analysis (7 pairs from H8)
  6. Operon context analysis
  7. Functional module detection
  8. Expression dynamics visualization
  9-10. Figures and tables output
"""

import os
import sys
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import seaborn as sns

# ============================================================
# Configuration
# ============================================================
BASE = "/Users/okaban/bioinfo/rna-seq"
ANALYSIS_DIR = f"{BASE}/11_epigenome_integration/analysis/55_exposed_regulatory_module"
FIG_DIR = f"{ANALYSIS_DIR}/figures"
TAB_DIR = f"{ANALYSIS_DIR}/tables"

CHROM_LEN = 8_667_507
ARM_BOUNDARY = 1_500_000  # 1.5 Mb from each end

# Input files
EXPOSED_TSV = f"{BASE}/11_epigenome_integration/analysis/51_exposed_regulators_characteristics/tables/exposed_regulators_full_table.tsv"
ALL_GENES_TSV = f"{BASE}/11_epigenome_integration/analysis/52_shielded_exposed_boundary/tables/all_genes_features.tsv"
COORD_GENES_TSV = f"{BASE}/11_epigenome_integration/analysis/29_genomewide_TF_screen/tables/coordinated_regulatory_genes.tsv"
NORM_COUNTS = f"{BASE}/04_deseq2/analysis/04_deseq2_260128_v1/results/normalized_counts_M145.tsv"
DESEQ_T2 = f"{BASE}/04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_2_vs_1.tsv"
DESEQ_T3 = f"{BASE}/04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_3_vs_1.tsv"
TCS_PAIRS = f"{BASE}/11_epigenome_integration/analysis/31_coordinated_regulators_characterization/tables/TCS_pair_candidates.tsv"

# Ensure output dirs exist
for d in [FIG_DIR, TAB_DIR]:
    os.makedirs(d, exist_ok=True)

np.random.seed(42)

# ============================================================
# STEP 1: Load data
# ============================================================
print("=" * 70)
print("STEP 1: Loading data")
print("=" * 70)

exposed = pd.read_csv(EXPOSED_TSV, sep='\t')
all_genes = pd.read_csv(ALL_GENES_TSV, sep='\t')
coord_genes = pd.read_csv(COORD_GENES_TSV, sep='\t')
norm_counts = pd.read_csv(NORM_COUNTS, sep='\t')
deseq_t2 = pd.read_csv(DESEQ_T2, sep='\t')
deseq_t3 = pd.read_csv(DESEQ_T3, sep='\t')
tcs_pairs = pd.read_csv(TCS_PAIRS, sep='\t')

# Remove empty rows in tcs_pairs
tcs_pairs = tcs_pairs.dropna(subset=['sensor_kinase']).copy()

print(f"Exposed regulators: {len(exposed)}")
print(f"All regulatory genes: {len(all_genes)}")
print(f"Coordinated regulatory genes: {len(coord_genes)}")
print(f"Normalized counts: {norm_counts.shape}")
print(f"TCS pairs: {len(tcs_pairs)}")

# Sample columns in normalized counts
sample_cols = [c for c in norm_counts.columns if c != 'gene_id']
print(f"Sample columns: {sample_cols}")

# Build position lookup from exposed table
exposed['midpoint'] = (exposed['start'] + exposed['end']) / 2
exposed['region_calc'] = exposed['midpoint'].apply(
    lambda x: 'left_arm' if x < ARM_BOUNDARY else ('right_arm' if x > (CHROM_LEN - ARM_BOUNDARY) else 'core')
)

# All regulatory gene positions
all_genes['midpoint'] = (all_genes['start'] + all_genes['end']) / 2

# Match exposed locus_tags
exposed_loci = set(exposed['locus_tag'].values)
print(f"\nExposed locus_tags (first 10): {sorted(exposed_loci)[:10]}")

# ============================================================
# STEP 2: Genomic clustering analysis
# ============================================================
print("\n" + "=" * 70)
print("STEP 2: Genomic clustering analysis")
print("=" * 70)

# Get midpoints for exposed TFs
exp_positions = np.sort(exposed['midpoint'].values)

# Nearest-neighbor distances for exposed TFs
def nearest_neighbor_distances(positions):
    """Calculate nearest-neighbor distance for each position."""
    pos = np.sort(positions)
    if len(pos) < 2:
        return np.array([])
    nn_dists = []
    for i, p in enumerate(pos):
        dists = np.abs(pos - p)
        dists[i] = np.inf
        nn_dists.append(np.min(dists))
    return np.array(nn_dists)

exp_nn_dists = nearest_neighbor_distances(exp_positions)
exp_median_nn = np.median(exp_nn_dists)
exp_mean_nn = np.mean(exp_nn_dists)

print(f"Exposed TFs: median NN distance = {exp_median_nn:,.0f} bp")
print(f"Exposed TFs: mean NN distance = {exp_mean_nn:,.0f} bp")

# Permutation test: sample 57 from all 1,055 regulatory genes
all_reg_positions = all_genes['midpoint'].values
n_perm = 1000
n_exposed = len(exposed)
perm_medians = []
perm_means = []

for i in range(n_perm):
    idx = np.random.choice(len(all_reg_positions), size=n_exposed, replace=False)
    sample_pos = all_reg_positions[idx]
    nn = nearest_neighbor_distances(sample_pos)
    perm_medians.append(np.median(nn))
    perm_means.append(np.mean(nn))

perm_medians = np.array(perm_medians)
perm_means = np.array(perm_means)

# One-sided test: are exposed TFs closer than random?
p_closer = np.mean(perm_medians <= exp_median_nn)
p_farther = np.mean(perm_medians >= exp_median_nn)

print(f"\nPermutation test (n={n_perm}):")
print(f"  Random median NN: {np.median(perm_medians):,.0f} bp (IQR: {np.percentile(perm_medians, 25):,.0f} - {np.percentile(perm_medians, 75):,.0f})")
print(f"  Exposed median NN: {exp_median_nn:,.0f} bp")
print(f"  P(random <= observed): {p_closer:.4f}")
print(f"  P(random >= observed): {p_farther:.4f}")

# Effect size: z-score
z_score = (exp_median_nn - np.mean(perm_medians)) / np.std(perm_medians)
print(f"  Z-score: {z_score:.3f}")

# Identify genomic clusters: exposed TFs within 20kb of each other
CLUSTER_DIST = 20_000  # 20 kb
exp_sorted = exposed.sort_values('midpoint').reset_index(drop=True)

# Single-linkage clustering on midpoints
cluster_id = 0
cluster_labels = [0] * len(exp_sorted)
cluster_labels[0] = cluster_id
for i in range(1, len(exp_sorted)):
    if exp_sorted.loc[i, 'midpoint'] - exp_sorted.loc[i-1, 'midpoint'] <= CLUSTER_DIST:
        cluster_labels[i] = cluster_id
    else:
        cluster_id += 1
        cluster_labels[i] = cluster_id

exp_sorted['genomic_cluster'] = cluster_labels

# Identify true clusters (size >= 2)
cluster_counts = exp_sorted['genomic_cluster'].value_counts()
multi_clusters = cluster_counts[cluster_counts >= 2].index.tolist()

print(f"\n20kb clusters found: {len(multi_clusters)}")
cluster_table_rows = []
for cl in sorted(multi_clusters):
    members = exp_sorted[exp_sorted['genomic_cluster'] == cl]
    span = members['midpoint'].max() - members['midpoint'].min()
    loci = ', '.join(members['locus_tag'].values)
    old_loci = ', '.join(members['old_locus_tag'].fillna('').values)
    families = ', '.join(members['tf_family'].values)
    region = members['region_calc'].mode().values[0] if len(members) > 0 else ''
    print(f"  Cluster {cl}: {len(members)} genes, span={span:,.0f} bp, region={region}")
    print(f"    Loci: {loci}")
    print(f"    Families: {families}")
    cluster_table_rows.append({
        'cluster_id': cl,
        'n_members': len(members),
        'span_bp': span,
        'region': region,
        'locus_tags': loci,
        'old_locus_tags': old_loci,
        'tf_families': families,
        'start': members['start'].min(),
        'end': members['end'].max(),
        'member_midpoints': '; '.join([f"{m:.0f}" for m in members['midpoint'].values])
    })

cluster_df = pd.DataFrame(cluster_table_rows)

# Count singletons
n_singletons = len(cluster_counts[cluster_counts == 1])
n_clustered = len(exp_sorted[exp_sorted['genomic_cluster'].isin(multi_clusters)])
print(f"\nClustered genes: {n_clustered}/{len(exposed)} ({100*n_clustered/len(exposed):.1f}%)")
print(f"Singletons: {n_singletons}")

# ============================================================
# STEP 3: Co-expression analysis
# ============================================================
print("\n" + "=" * 70)
print("STEP 3: Co-expression analysis")
print("=" * 70)

# Extract expression for all regulatory genes
# Use gene_id column from norm_counts matching locus_tag
all_reg_loci = set(all_genes['locus_tag'].values)
exposed_loci_list = list(exposed['locus_tag'].values)

# Get expression matrix for exposed genes
expr_exposed = norm_counts[norm_counts['gene_id'].isin(exposed_loci)].copy()
expr_exposed = expr_exposed.set_index('gene_id')
# Keep only the 9 sample columns
expr_exposed = expr_exposed[sample_cols]

print(f"Exposed genes with expression data: {len(expr_exposed)}/{len(exposed)}")

# Get expression for shielded genes (all regulatory minus exposed)
shielded_loci = all_reg_loci - exposed_loci
expr_shielded = norm_counts[norm_counts['gene_id'].isin(shielded_loci)].copy()
expr_shielded = expr_shielded.set_index('gene_id')
expr_shielded = expr_shielded[sample_cols]

print(f"Shielded genes with expression data: {len(expr_shielded)}")

# Calculate pairwise Spearman correlations for exposed
def pairwise_spearman(df):
    """Calculate pairwise Spearman correlations."""
    genes = df.index.tolist()
    n = len(genes)
    corr_mat = np.zeros((n, n))
    pval_mat = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            if i == j:
                corr_mat[i, j] = 1.0
                pval_mat[i, j] = 0.0
            else:
                rho, pv = stats.spearmanr(df.iloc[i].values, df.iloc[j].values)
                corr_mat[i, j] = rho
                corr_mat[j, i] = rho
                pval_mat[i, j] = pv
                pval_mat[j, i] = pv
    return pd.DataFrame(corr_mat, index=genes, columns=genes), pd.DataFrame(pval_mat, index=genes, columns=genes)

corr_exp, pval_exp = pairwise_spearman(expr_exposed)
print(f"\nExposed 57x57 correlation matrix computed")

# Compare correlation distributions: exposed-exposed, exposed-shielded, shielded-shielded
# Extract upper triangle values for exposed-exposed
mask_upper = np.triu_indices(len(corr_exp), k=1)
ee_corrs = corr_exp.values[mask_upper]

# For exposed-shielded: sample 200 shielded genes for computational tractability
if len(expr_shielded) > 200:
    shielded_sample_idx = np.random.choice(len(expr_shielded), 200, replace=False)
    expr_shielded_sample = expr_shielded.iloc[shielded_sample_idx]
else:
    expr_shielded_sample = expr_shielded

# exposed-shielded correlations
es_corrs = []
for exp_gene in expr_exposed.index:
    for sh_gene in expr_shielded_sample.index:
        rho, _ = stats.spearmanr(expr_exposed.loc[exp_gene].values, expr_shielded_sample.loc[sh_gene].values)
        es_corrs.append(rho)
es_corrs = np.array(es_corrs)

# shielded-shielded correlations (sample pairs for speed)
ss_genes = expr_shielded_sample.index.tolist()
ss_mask = np.triu_indices(len(ss_genes), k=1)
ss_corr_mat = np.zeros((len(ss_genes), len(ss_genes)))
for i in range(len(ss_genes)):
    for j in range(i+1, len(ss_genes)):
        rho, _ = stats.spearmanr(expr_shielded_sample.iloc[i].values, expr_shielded_sample.iloc[j].values)
        ss_corr_mat[i, j] = rho
ss_corrs = ss_corr_mat[ss_mask]

print(f"\nCorrelation distributions:")
print(f"  Exposed-Exposed (n={len(ee_corrs)}): median={np.median(ee_corrs):.4f}, mean={np.mean(ee_corrs):.4f}")
print(f"  Exposed-Shielded (n={len(es_corrs)}): median={np.median(es_corrs):.4f}, mean={np.mean(es_corrs):.4f}")
print(f"  Shielded-Shielded (n={len(ss_corrs)}): median={np.median(ss_corrs):.4f}, mean={np.mean(ss_corrs):.4f}")

# Statistical tests
# Mann-Whitney U: E-E vs E-S
u_ee_es, p_ee_es = stats.mannwhitneyu(ee_corrs, es_corrs, alternative='greater')
r_ee_es = u_ee_es / (len(ee_corrs) * len(es_corrs))  # rank-biserial effect size
# E-E vs S-S
u_ee_ss, p_ee_ss = stats.mannwhitneyu(ee_corrs, ss_corrs, alternative='greater')
r_ee_ss = u_ee_ss / (len(ee_corrs) * len(ss_corrs))

print(f"\n  E-E vs E-S: U={u_ee_es:.0f}, p={p_ee_es:.4e}, rank-biserial r={r_ee_es:.4f}")
print(f"  E-E vs S-S: U={u_ee_ss:.0f}, p={p_ee_ss:.4e}, rank-biserial r={r_ee_ss:.4f}")

# Hierarchical clustering
linkage_mat = linkage(pdist(expr_exposed.values, metric='correlation'), method='ward')

# Co-expression modules: correlation > 0.7 and BH-adjusted p < 0.05
# Flatten upper triangle p-values and apply BH correction
pvals_flat = pval_exp.values[mask_upper]
corrs_flat = corr_exp.values[mask_upper]
reject, padj, _, _ = multipletests(pvals_flat, method='fdr_bh')

# Module detection: find connected components of significant high-correlation pairs
sig_high_corr = (corrs_flat > 0.7) & (padj < 0.05)
gene_names = corr_exp.index.tolist()
n_genes = len(gene_names)

# Build adjacency from significant pairs
from collections import defaultdict

adj = defaultdict(set)
pair_idx = 0
for i in range(n_genes):
    for j in range(i+1, n_genes):
        if sig_high_corr[pair_idx]:
            adj[gene_names[i]].add(gene_names[j])
            adj[gene_names[j]].add(gene_names[i])
        pair_idx += 1

# Find connected components
visited = set()
coexpr_modules = []

def bfs(start, adj, visited):
    queue = [start]
    component = set()
    while queue:
        node = queue.pop(0)
        if node in visited:
            continue
        visited.add(node)
        component.add(node)
        for neighbor in adj.get(node, []):
            if neighbor not in visited:
                queue.append(neighbor)
    return component

for gene in gene_names:
    if gene not in visited and gene in adj:
        component = bfs(gene, adj, visited)
        if len(component) >= 2:
            coexpr_modules.append(sorted(component))

print(f"\nCo-expression modules (rho>0.7, padj<0.05): {len(coexpr_modules)}")
coexpr_module_rows = []
for idx, mod in enumerate(coexpr_modules):
    # Get mean correlation within module
    mod_idx = [gene_names.index(g) for g in mod]
    within_corrs = []
    for i in range(len(mod_idx)):
        for j in range(i+1, len(mod_idx)):
            within_corrs.append(corr_exp.values[mod_idx[i], mod_idx[j]])
    mean_rho = np.mean(within_corrs) if within_corrs else 0

    # Get TF families and coordination types
    mod_info = exposed[exposed['locus_tag'].isin(mod)]
    families = ', '.join(sorted(mod_info['tf_family'].unique()))
    coord_t2 = ', '.join(sorted(mod_info['coordination_T2'].dropna().unique())) if 'coordination_T2' in mod_info.columns else 'N/A'
    coord_t3 = ', '.join(sorted(mod_info['coordination_T3'].dropna().unique())) if 'coordination_T3' in mod_info.columns else 'N/A'
    regions = ', '.join(sorted(mod_info['region_calc'].unique())) if 'region_calc' in mod_info.columns else ', '.join(sorted(mod_info['region'].dropna().unique())) if 'region' in mod_info.columns else 'N/A'

    print(f"  Module {idx+1}: {len(mod)} genes, mean rho={mean_rho:.3f}")
    print(f"    Genes: {', '.join(mod)}")
    print(f"    Families: {families}")

    coexpr_module_rows.append({
        'module_id': idx+1,
        'n_members': len(mod),
        'mean_rho': round(mean_rho, 4),
        'locus_tags': ', '.join(mod),
        'old_locus_tags': ', '.join(mod_info.set_index('locus_tag').loc[mod, 'old_locus_tag'].fillna('').values),
        'tf_families': families,
        'coordination_T2': coord_t2,
        'coordination_T3': coord_t3,
        'regions': regions
    })

coexpr_modules_df = pd.DataFrame(coexpr_module_rows)

# Also count significant pairs
n_sig_pairs = np.sum(sig_high_corr)
n_total_pairs = len(corrs_flat)
print(f"\nSignificant high-correlation pairs: {n_sig_pairs}/{n_total_pairs} ({100*n_sig_pairs/n_total_pairs:.1f}%)")

# ============================================================
# STEP 4: Coordination type coherence
# ============================================================
print("\n" + "=" * 70)
print("STEP 4: Coordination type coherence")
print("=" * 70)

# Coordination type columns not available in updated dataset — skip coherence analysis
exp_coord = exposed.copy()
for col in ['coordination_T2_final', 'coordination_T3_final']:
    exp_coord[col] = 'N/A'

print("  (Coordination type data not available in Jeong2016-filtered dataset; skipped)")

print("\nCoordination type distribution (T2): N/A")
print("Coordination type distribution (T3): N/A")
ct_t2 = pd.Series({'N/A': len(exp_coord)})
ct_t3 = pd.Series({'N/A': len(exp_coord)})

# Save coordination type distribution
coord_dist_rows = []
for t_label, ct_series in [('T2', ct_t2), ('T3', ct_t3)]:
    for ct, n in ct_series.items():
        coord_dist_rows.append({
            'timepoint': t_label,
            'coordination_type': ct,
            'n_genes': n,
            'pct': round(100 * n / len(exp_coord), 1)
        })
coord_dist_df = pd.DataFrame(coord_dist_rows)

# Add genomic cluster membership to exp_coord
exp_coord_with_cluster = exp_coord.merge(
    exp_sorted[['locus_tag', 'genomic_cluster']],
    on='locus_tag', how='left'
)
# Mark whether gene is in a multi-member cluster
exp_coord_with_cluster['in_cluster'] = exp_coord_with_cluster['genomic_cluster'].isin(multi_clusters)

# Chi-square: coordination_T3 vs cluster membership
# Create contingency table
ct_types_t3 = exp_coord_with_cluster['coordination_T3_final'].unique()
contingency_data = []
for ct in ct_types_t3:
    mask = exp_coord_with_cluster['coordination_T3_final'] == ct
    n_in = mask & exp_coord_with_cluster['in_cluster']
    n_out = mask & ~exp_coord_with_cluster['in_cluster']
    contingency_data.append([n_in.sum(), n_out.sum()])

contingency = np.array(contingency_data)
# Filter out rows with all zeros
nonzero_rows = contingency.sum(axis=1) > 0
contingency_filtered = contingency[nonzero_rows]
ct_types_filtered = [ct_types_t3[i] for i in range(len(ct_types_t3)) if nonzero_rows[i]]

if contingency_filtered.shape[0] >= 2 and contingency_filtered.shape[1] >= 2:
    chi2_ct_cluster, p_ct_cluster, dof_ct_cluster, expected_ct = stats.chi2_contingency(contingency_filtered)
    cramers_v = np.sqrt(chi2_ct_cluster / (contingency_filtered.sum() * (min(contingency_filtered.shape) - 1)))
    print(f"\nChi-square: coordination_T3 vs genomic cluster membership")
    print(f"  chi2={chi2_ct_cluster:.3f}, p={p_ct_cluster:.4f}, dof={dof_ct_cluster}, Cramer's V={cramers_v:.4f}")
else:
    chi2_ct_cluster, p_ct_cluster, dof_ct_cluster, cramers_v = np.nan, np.nan, np.nan, np.nan
    print(f"\nChi-square: insufficient data for contingency test")

# Do TFs with same coordination type co-express more?
coord_types_main = ['discordant_gain_up', 'concordant_derepression', 'concordant_repression',
                    'methyl_change_no_expr_change', 'ambiguous']

within_type_corrs = []
between_type_corrs = []

for i in range(len(gene_names)):
    for j in range(i+1, len(gene_names)):
        gi = gene_names[i]
        gj = gene_names[j]
        ci = exp_coord.loc[exp_coord['locus_tag'] == gi, 'coordination_T3_final'].values
        cj = exp_coord.loc[exp_coord['locus_tag'] == gj, 'coordination_T3_final'].values
        if len(ci) > 0 and len(cj) > 0:
            if ci[0] == cj[0]:
                within_type_corrs.append(corr_exp.values[i, j])
            else:
                between_type_corrs.append(corr_exp.values[i, j])

within_type_corrs = np.array(within_type_corrs)
between_type_corrs = np.array(between_type_corrs)

if len(within_type_corrs) > 0 and len(between_type_corrs) > 0:
    u_wt, p_wt = stats.mannwhitneyu(within_type_corrs, between_type_corrs, alternative='greater')
    r_wt = u_wt / (len(within_type_corrs) * len(between_type_corrs))
    print(f"\nWithin-type vs between-type co-expression (T3):")
    print(f"  Within-type: n={len(within_type_corrs)}, median={np.median(within_type_corrs):.4f}")
    print(f"  Between-type: n={len(between_type_corrs)}, median={np.median(between_type_corrs):.4f}")
    print(f"  U={u_wt:.0f}, p={p_wt:.4e}, rank-biserial r={r_wt:.4f}")
else:
    u_wt, p_wt, r_wt = np.nan, np.nan, np.nan
    print("\nInsufficient data for within/between-type comparison")

# ============================================================
# STEP 5: TCS pair analysis
# ============================================================
print("\n" + "=" * 70)
print("STEP 5: TCS pair analysis (7 pairs from H8)")
print("=" * 70)

tcs_analysis_rows = []
for _, row in tcs_pairs.iterrows():
    sk = row['sensor_kinase']
    rr = row['response_regulator']
    sk_old = row['sk_old_locus']
    rr_old = row['rr_old_locus']

    sk_exposed = sk in exposed_loci
    rr_exposed = rr in exposed_loci

    # Expression correlation
    sk_expr = norm_counts[norm_counts['gene_id'] == sk]
    rr_expr = norm_counts[norm_counts['gene_id'] == rr]

    if len(sk_expr) > 0 and len(rr_expr) > 0:
        sk_vals = sk_expr[sample_cols].values.flatten()
        rr_vals = rr_expr[sample_cols].values.flatten()
        rho_pair, p_pair = stats.spearmanr(sk_vals, rr_vals)
    else:
        rho_pair, p_pair = np.nan, np.nan

    # Which partner has methylation?
    sk_coord = coord_genes[coord_genes['locus_tag'] == sk]
    rr_coord = coord_genes[coord_genes['locus_tag'] == rr]

    sk_has_methyl = sk_coord['has_methylation'].values[0] if len(sk_coord) > 0 else False
    rr_has_methyl = rr_coord['has_methylation'].values[0] if len(rr_coord) > 0 else False

    methyl_partner = 'both' if (sk_has_methyl and rr_has_methyl) else \
                     ('sensor_kinase' if sk_has_methyl else \
                     ('response_regulator' if rr_has_methyl else 'neither'))

    print(f"\n  Pair: {sk_old} (SK) - {rr_old} (RR)")
    print(f"    SK exposed: {sk_exposed}, RR exposed: {rr_exposed}")
    print(f"    Methylated partner: {methyl_partner}")
    print(f"    Expression correlation: rho={rho_pair:.3f}, p={p_pair:.4e}" if not np.isnan(rho_pair) else "    No expression data")
    print(f"    SK coordination T2: {row['sk_coordination_T2']}, T3: {row['sk_coordination_T3']}")
    print(f"    RR coordination T2: {row['rr_coordination_T2']}, T3: {row['rr_coordination_T3']}")

    tcs_analysis_rows.append({
        'sensor_kinase': sk,
        'sk_old_locus': sk_old,
        'response_regulator': rr,
        'rr_old_locus': rr_old,
        'sk_exposed': sk_exposed,
        'rr_exposed': rr_exposed,
        'both_exposed': sk_exposed and rr_exposed,
        'methylated_partner': methyl_partner,
        'expression_rho': round(rho_pair, 4) if not np.isnan(rho_pair) else np.nan,
        'expression_pval': p_pair if not np.isnan(p_pair) else np.nan,
        'sk_coordination_T2': row['sk_coordination_T2'],
        'sk_coordination_T3': row['sk_coordination_T3'],
        'rr_coordination_T2': row['rr_coordination_T2'],
        'rr_coordination_T3': row['rr_coordination_T3'],
        'distance_bp': row['distance_bp'],
        'both_in_57_H8': row['both_in_57']
    })

tcs_df = pd.DataFrame(tcs_analysis_rows)

# Summary statistics
n_sk_methyl = (tcs_df['methylated_partner'] == 'sensor_kinase').sum()
n_rr_methyl = (tcs_df['methylated_partner'] == 'response_regulator').sum()
n_both_methyl = (tcs_df['methylated_partner'] == 'both').sum()
print(f"\nMethylation on sensor kinase: {n_sk_methyl}")
print(f"Methylation on response regulator: {n_rr_methyl}")
print(f"Methylation on both: {n_both_methyl}")
print(f"Mean TCS pair correlation: {tcs_df['expression_rho'].mean():.3f}")

# Is methylated partner consistently one type? Binomial test
n_asym = n_sk_methyl + n_rr_methyl
if n_asym > 0:
    p_binom = stats.binomtest(max(n_sk_methyl, n_rr_methyl), n_asym, 0.5).pvalue
    print(f"Binomial test for asymmetry: p={p_binom:.4f}")
else:
    p_binom = np.nan

# ============================================================
# STEP 6: Operon context analysis
# ============================================================
print("\n" + "=" * 70)
print("STEP 6: Operon context analysis")
print("=" * 70)

# Operon proxy: intergenic distance < 150bp and same strand
# Check all pairwise exposed TFs sorted by position
exp_op = exp_sorted.copy()
operon_pairs = []
for i in range(len(exp_op) - 1):
    for j in range(i+1, min(i+5, len(exp_op))):  # Check nearby genes
        g1 = exp_op.iloc[i]
        g2 = exp_op.iloc[j]

        # Same strand?
        if g1['strand'] == g2['strand']:
            # Intergenic distance
            if g1['strand'] == '+':
                intergenic = g2['start'] - g1['end']
            else:
                intergenic = g1['start'] - g2['end']

            # Use absolute intergenic as proxy
            intergenic = abs(g2['start'] - g1['end'])

            if 0 <= intergenic <= 150:
                operon_pairs.append({
                    'gene1': g1['locus_tag'],
                    'gene1_old': g1['old_locus_tag'],
                    'gene2': g2['locus_tag'],
                    'gene2_old': g2['old_locus_tag'],
                    'strand': g1['strand'],
                    'intergenic_dist': intergenic,
                    'gene1_family': g1['tf_family'],
                    'gene2_family': g2['tf_family']
                })

print(f"Exposed-exposed operonic pairs (same strand, <150bp): {len(operon_pairs)}")
for pair in operon_pairs:
    print(f"  {pair['gene1_old']} ({pair['gene1_family']}) - {pair['gene2_old']} ({pair['gene2_family']}), intergenic={pair['intergenic_dist']}bp, strand={pair['strand']}")

# Permutation test: how many operonic pairs expected by random sampling 57 from all_genes?
n_perm_operon = 1000
all_genes_sorted = all_genes.sort_values('start').reset_index(drop=True)

random_operon_counts = []
for _ in range(n_perm_operon):
    idx = np.random.choice(len(all_genes_sorted), size=n_exposed, replace=False)
    sample = all_genes_sorted.iloc[sorted(idx)].reset_index(drop=True)
    count = 0
    for i in range(len(sample) - 1):
        for j in range(i+1, min(i+5, len(sample))):
            g1 = sample.iloc[i]
            g2 = sample.iloc[j]
            if g1['strand'] == g2['strand']:
                intergenic = abs(g2['start'] - g1['end'])
                if 0 <= intergenic <= 150:
                    count += 1
    random_operon_counts.append(count)

random_operon_counts = np.array(random_operon_counts)
p_operon = np.mean(random_operon_counts >= len(operon_pairs))
print(f"\nPermutation test: observed={len(operon_pairs)}, random median={np.median(random_operon_counts):.1f}, p={p_operon:.4f}")

# ============================================================
# STEP 7: Functional module detection
# ============================================================
print("\n" + "=" * 70)
print("STEP 7: Functional module detection")
print("=" * 70)

# Modules: groups of >=3 exposed TFs that are:
# (a) Within 50kb of each other, AND/OR
# (b) Have expression correlation > 0.6, AND
# (c) Share coordination type

# Step 7a: 50kb proximity groups
PROX_DIST = 50_000
cluster_labels_50k = [0] * len(exp_sorted)
cl_id = 0
cluster_labels_50k[0] = cl_id
for i in range(1, len(exp_sorted)):
    if exp_sorted.iloc[i]['midpoint'] - exp_sorted.iloc[i-1]['midpoint'] <= PROX_DIST:
        cluster_labels_50k[i] = cl_id
    else:
        cl_id += 1
        cluster_labels_50k[i] = cl_id
exp_sorted['prox_cluster_50k'] = cluster_labels_50k

# Step 7b: expression correlation network (rho > 0.6)
adj_06 = defaultdict(set)
pair_idx = 0
for i in range(len(gene_names)):
    for j in range(i+1, len(gene_names)):
        if corr_exp.values[i, j] > 0.6:
            adj_06[gene_names[i]].add(gene_names[j])
            adj_06[gene_names[j]].add(gene_names[i])
        pair_idx += 1

# Step 7c: Combine - for each gene, check proximity group + coexpression + coordination
# Strategy: Start with genomic proximity clusters of >=2, expand with coexpression,
# then filter by shared coordination type

# Get coordination type lookup
coord_lookup = {}
for _, row in exp_coord.iterrows():
    coord_lookup[row['locus_tag']] = {
        'T2': row['coordination_T2_final'],
        'T3': row['coordination_T3_final']
    }

# Build functional modules
functional_modules = []

# Approach: enumerate candidate modules from 50kb clusters, coexpression clusters, or their union
# Then score by multi-evidence consistency

# Get 50kb clusters with >= 2 members
prox_50k_counts = pd.Series(cluster_labels_50k).value_counts()
prox_50k_multi = prox_50k_counts[prox_50k_counts >= 2].index.tolist()

candidate_groups = []

# From proximity
for cl in prox_50k_multi:
    members = exp_sorted[exp_sorted['prox_cluster_50k'] == cl]['locus_tag'].tolist()
    candidate_groups.append(('proximity_50k', set(members)))

# From coexpression modules
visited_coexpr = set()
for gene in gene_names:
    if gene not in visited_coexpr and gene in adj_06:
        component = bfs(gene, adj_06, set())
        if len(component) >= 2:
            candidate_groups.append(('coexpression_06', component))
            visited_coexpr.update(component)

# Merge overlapping candidates
all_candidate_members = set()
for _, members in candidate_groups:
    all_candidate_members.update(members)

# For each candidate group, check coordination type coherence
module_rows = []
module_id = 0
for source, members in candidate_groups:
    members_list = sorted(members)
    # Get coordination types
    coord_types_t3 = [coord_lookup.get(g, {}).get('T3', 'unknown') for g in members_list]
    coord_types_t2 = [coord_lookup.get(g, {}).get('T2', 'unknown') for g in members_list]

    # Most common coordination type
    from collections import Counter
    ct_counter_t3 = Counter(coord_types_t3)
    most_common_t3, most_common_count_t3 = ct_counter_t3.most_common(1)[0]
    coherence_t3 = most_common_count_t3 / len(members_list)

    ct_counter_t2 = Counter(coord_types_t2)
    most_common_t2, most_common_count_t2 = ct_counter_t2.most_common(1)[0]
    coherence_t2 = most_common_count_t2 / len(members_list)

    # Check proximity
    member_positions = exp_sorted[exp_sorted['locus_tag'].isin(members)]['midpoint'].values
    if len(member_positions) >= 2:
        max_span = member_positions.max() - member_positions.min()
    else:
        max_span = 0
    is_proximal = max_span <= 50_000

    # Check coexpression
    member_idx = [gene_names.index(g) for g in members_list if g in gene_names]
    if len(member_idx) >= 2:
        within_corrs_mod = []
        for i in range(len(member_idx)):
            for j in range(i+1, len(member_idx)):
                within_corrs_mod.append(corr_exp.values[member_idx[i], member_idx[j]])
        mean_corr = np.mean(within_corrs_mod)
        is_coexpressed = mean_corr > 0.6
    else:
        mean_corr = 0
        is_coexpressed = False

    # Score: proximity (1pt) + coexpression (1pt) + coordination coherence T3 (0-1pt)
    score = (1 if is_proximal else 0) + (1 if is_coexpressed else 0) + coherence_t3

    # Check if this qualifies as a module (>= 3 members with shared features)
    # Or >= 2 with high score
    if len(members_list) >= 2:
        module_id += 1

        info = exposed[exposed['locus_tag'].isin(members)]
        families = ', '.join(sorted(info['tf_family'].unique()))
        regions = ', '.join(sorted(info['region_calc'].unique()))

        module_rows.append({
            'module_id': module_id,
            'source': source,
            'n_members': len(members_list),
            'locus_tags': ', '.join(members_list),
            'old_locus_tags': ', '.join(info.set_index('locus_tag').reindex(members_list)['old_locus_tag'].fillna('').values),
            'tf_families': families,
            'span_bp': max_span,
            'is_proximal': is_proximal,
            'mean_coexpression': round(mean_corr, 4),
            'is_coexpressed': is_coexpressed,
            'dominant_coordination_T3': most_common_t3,
            'coordination_coherence_T3': round(coherence_t3, 3),
            'dominant_coordination_T2': most_common_t2,
            'coordination_coherence_T2': round(coherence_t2, 3),
            'evidence_score': round(score, 3),
            'region': regions
        })

module_df = pd.DataFrame(module_rows)
if len(module_df) > 0:
    module_df = module_df.sort_values('evidence_score', ascending=False).reset_index(drop=True)
    module_df['module_id'] = range(1, len(module_df) + 1)

print(f"Candidate functional modules detected: {len(module_df)}")
if len(module_df) > 0:
    for _, row in module_df.iterrows():
        print(f"  Module {row['module_id']}: {row['n_members']} genes, score={row['evidence_score']:.2f}")
        print(f"    Source: {row['source']}, proximal: {row['is_proximal']}, coexpressed: {row['is_coexpressed']}")
        print(f"    Dominant T3 type: {row['dominant_coordination_T3']} ({row['coordination_coherence_T3']:.0%})")
        print(f"    Genes: {row['locus_tags']}")

# High-confidence modules (score >= 2 and >= 3 members)
hc_modules = module_df[(module_df['evidence_score'] >= 2.0) & (module_df['n_members'] >= 3)] if len(module_df) > 0 else pd.DataFrame()
print(f"\nHigh-confidence modules (score>=2, >=3 members): {len(hc_modules)}")

# ============================================================
# STEP 8: Expression dynamics (z-scored heatmap data)
# ============================================================
print("\n" + "=" * 70)
print("STEP 8: Expression dynamics preparation")
print("=" * 70)

# Z-score the expression data
expr_z = expr_exposed.copy()
expr_z = expr_z.apply(lambda x: (x - x.mean()) / x.std() if x.std() > 0 else 0, axis=1)
print(f"Z-scored expression matrix: {expr_z.shape}")

# Annotation vectors
annot = exposed.set_index('locus_tag').loc[expr_exposed.index]
coord_type_vec = []
for lt in expr_exposed.index:
    ct = exp_coord.loc[exp_coord['locus_tag'] == lt, 'coordination_T3_final'].values
    coord_type_vec.append(ct[0] if len(ct) > 0 else 'unknown')

# ============================================================
# STEP 9: Generate Figures
# ============================================================
print("\n" + "=" * 70)
print("STEP 9: Generating figures")
print("=" * 70)

# Color palettes
coord_colors = {
    'discordant_gain_up': '#E74C3C',      # red
    'concordant_derepression': '#3498DB',   # blue
    'concordant_repression': '#2ECC71',     # green
    'methyl_change_no_expr_change': '#95A5A6',  # gray
    'ambiguous': '#F39C12',                 # orange
    'discordant_loss_down': '#9B59B6',      # purple
    'unknown': '#BDC3C7'
}

family_colors = {
    'TetR': '#E74C3C',
    'Sigma factor': '#3498DB',
    'Sensor kinase': '#2ECC71',
    'Response regulator': '#9B59B6',
    'HTH (other)': '#F39C12',
    'LysR': '#1ABC9C',
    'GntR': '#E67E22',
    'MarR': '#8E44AD',
    'ArsR': '#2980B9',
    'MerR': '#27AE60',
    'WhiB': '#D35400',
    'Xre': '#C0392B',
    'AraC': '#16A085',
    'IclR': '#7F8C8D',
    'LuxR': '#34495E',
}

region_colors = {
    'left_arm': '#E74C3C',
    'core': '#3498DB',
    'right_arm': '#2ECC71'
}

# ===== Figure 1: Genomic distribution =====
print("  Creating genomic_distribution figure...")
fig, axes = plt.subplots(3, 1, figsize=(14, 10), gridspec_kw={'height_ratios': [3, 1, 1]})

# Top panel: chromosome map
ax = axes[0]
# Draw chromosome
ax.barh(0, CHROM_LEN, height=0.3, color='#ECF0F1', edgecolor='#95A5A6', linewidth=1)
# Mark arm boundaries
for boundary in [ARM_BOUNDARY, CHROM_LEN - ARM_BOUNDARY]:
    ax.axvline(boundary, color='#7F8C8D', linestyle='--', linewidth=0.8, alpha=0.5)

# Label regions
ax.text(ARM_BOUNDARY / 2, 0.5, 'Left arm', ha='center', fontsize=10, color='#95A5A6')
ax.text(CHROM_LEN / 2, 0.5, 'Core', ha='center', fontsize=10, color='#95A5A6')
ax.text(CHROM_LEN - ARM_BOUNDARY / 2, 0.5, 'Right arm', ha='center', fontsize=10, color='#95A5A6')

# Plot exposed TFs
for _, gene in exp_sorted.iterrows():
    color = coord_colors.get(gene.get('coordination_T3', 'unknown'), '#BDC3C7')
    ax.plot(gene['midpoint'], 0, 'v', color=color, markersize=8, alpha=0.8)

# Highlight clusters
for _, cl in cluster_df.iterrows():
    rect = Rectangle((cl['start'] - 5000, -0.25), cl['end'] - cl['start'] + 10000, 0.5,
                     linewidth=1.5, edgecolor='red', facecolor='red', alpha=0.15)
    ax.add_patch(rect)

ax.set_xlim(-100000, CHROM_LEN + 100000)
ax.set_ylim(-0.6, 0.8)
ax.set_yticks([])
ax.set_xlabel('Chromosomal position (bp)', fontsize=11)
ax.set_title(f'Genomic Distribution of 57 Exposed Regulatory Genes\n(Clusters within 20kb highlighted in red)', fontsize=12)

# Format x-axis as Mb
ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x/1e6:.1f}'))
ax.set_xlabel('Chromosomal position (Mb)', fontsize=11)

# Legend
legend_elements = [Line2D([0], [0], marker='v', color='w', markerfacecolor=v, label=k, markersize=8)
                   for k, v in coord_colors.items() if k != 'unknown']
ax.legend(handles=legend_elements, loc='upper right', fontsize=7, ncol=2, title='Coordination T3')

# Middle panel: density
ax = axes[1]
positions_mb = exp_sorted['midpoint'].values / 1e6
ax.hist(positions_mb, bins=50, color='#3498DB', alpha=0.7, edgecolor='white')
ax.axvline(ARM_BOUNDARY / 1e6, color='#7F8C8D', linestyle='--', linewidth=0.8)
ax.axvline((CHROM_LEN - ARM_BOUNDARY) / 1e6, color='#7F8C8D', linestyle='--', linewidth=0.8)
ax.set_ylabel('Count', fontsize=10)
ax.set_xlabel('Position (Mb)', fontsize=10)
ax.set_title('Exposed TF density along chromosome', fontsize=11)

# Bottom panel: nearest-neighbor distance histogram comparison
ax = axes[2]
ax.hist(perm_medians / 1e3, bins=30, color='#95A5A6', alpha=0.7, label=f'Random (n={n_perm})', edgecolor='white')
ax.axvline(exp_median_nn / 1e3, color='#E74C3C', linewidth=2, label=f'Observed: {exp_median_nn/1e3:.1f} kb')
ax.set_xlabel('Median nearest-neighbor distance (kb)', fontsize=10)
ax.set_ylabel('Count', fontsize=10)
ax.set_title(f'Permutation test: genomic clustering (p={p_farther:.3f} for dispersal, p={p_closer:.3f} for clustering)', fontsize=11)
ax.legend(fontsize=9)

plt.tight_layout()
for fmt in ['pdf', 'svg']:
    fig.savefig(f"{FIG_DIR}/genomic_distribution.{fmt}", dpi=300, bbox_inches='tight')
plt.close()

# ===== Figure 2: Co-expression heatmap =====
print("  Creating coexpression_heatmap figure...")
fig = plt.figure(figsize=(14, 12))

# Cluster the correlation matrix
row_linkage = linkage(pdist(corr_exp.values, metric='euclidean'), method='ward')
col_linkage = row_linkage  # symmetric

g = sns.clustermap(corr_exp, row_linkage=row_linkage, col_linkage=col_linkage,
                   cmap='RdBu_r', center=0, vmin=-1, vmax=1,
                   figsize=(14, 12),
                   xticklabels=True, yticklabels=True,
                   dendrogram_ratio=0.15,
                   cbar_kws={'label': 'Spearman rho'})
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), fontsize=5, rotation=90)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=5)
g.fig.suptitle('Co-expression Matrix: 57 Exposed Regulatory Genes\n(Spearman correlation, Ward linkage)',
               fontsize=13, y=1.02)

for fmt in ['pdf', 'svg']:
    g.savefig(f"{FIG_DIR}/coexpression_heatmap.{fmt}", dpi=300, bbox_inches='tight')
plt.close('all')

# ===== Figure 3: Coordination type analysis =====
print("  Creating coordination_type_analysis figure...")
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Panel A: T2 and T3 coordination type pie charts
for idx, (t_label, ct_series) in enumerate([('T2', ct_t2), ('T3', ct_t3)]):
    ax = axes[0, idx]
    labels = ct_series.index.tolist()
    sizes = ct_series.values
    colors = [coord_colors.get(l, '#BDC3C7') for l in labels]
    wedges, texts, autotexts = ax.pie(sizes, labels=None, autopct='%1.0f%%',
                                       colors=colors, startangle=90, pctdistance=0.8,
                                       textprops={'fontsize': 8})
    ax.set_title(f'Coordination Type Distribution ({t_label})', fontsize=11)
    ax.legend(labels, loc='center left', bbox_to_anchor=(0.95, 0.5), fontsize=7)

# Panel C: Within-type vs between-type co-expression boxplot
ax = axes[1, 0]
bp_data = [within_type_corrs, between_type_corrs]
bp_labels = [f'Within-type\n(n={len(within_type_corrs)})', f'Between-type\n(n={len(between_type_corrs)})']
bp = ax.boxplot(bp_data, labels=bp_labels, patch_artist=True, widths=0.5)
bp['boxes'][0].set_facecolor('#3498DB')
bp['boxes'][1].set_facecolor('#95A5A6')
for element in ['whiskers', 'caps', 'medians']:
    for line in bp[element]:
        line.set_color('black')
ax.set_ylabel('Spearman rho', fontsize=10)
ax.set_title('Co-expression by Coordination Type', fontsize=11)
if not np.isnan(p_wt):
    sig_marker = '***' if p_wt < 0.001 else ('**' if p_wt < 0.01 else ('*' if p_wt < 0.05 else 'ns'))
    ax.text(1.5, max(np.max(within_type_corrs), np.max(between_type_corrs)) + 0.05,
            f'p={p_wt:.3e}\n{sig_marker}', ha='center', fontsize=9)

# Panel D: Coordination type by region
ax = axes[1, 1]
exp_coord_with_cluster['region_calc'] = exp_coord_with_cluster['locus_tag'].map(
    exposed.set_index('locus_tag')['region_calc'])
ct_region = pd.crosstab(exp_coord_with_cluster['coordination_T3_final'],
                        exp_coord_with_cluster['region_calc'])
ct_region_pct = ct_region.div(ct_region.sum(axis=1), axis=0) * 100
ct_region_pct.plot(kind='barh', stacked=True, ax=ax,
                   color=[region_colors.get(c, '#BDC3C7') for c in ct_region_pct.columns])
ax.set_xlabel('Percentage', fontsize=10)
ax.set_title('Coordination Type by Genomic Region (T3)', fontsize=11)
ax.legend(title='Region', fontsize=8, loc='lower right')

plt.tight_layout()
for fmt in ['pdf', 'svg']:
    fig.savefig(f"{FIG_DIR}/coordination_type_analysis.{fmt}", dpi=300, bbox_inches='tight')
plt.close()

# ===== Figure 4: TCS pairs =====
print("  Creating TCS_pairs figure...")
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Panel A: Expression correlation per pair
ax = axes[0]
pair_labels = [f"{row['sk_old_locus']}\n{row['rr_old_locus']}" for _, row in tcs_df.iterrows()]
rhos = tcs_df['expression_rho'].values
colors_tcs = ['#E74C3C' if m == 'sensor_kinase' else '#3498DB' if m == 'response_regulator' else '#95A5A6'
              for m in tcs_df['methylated_partner']]
bars = ax.bar(range(len(rhos)), rhos, color=colors_tcs, edgecolor='white', width=0.6)
ax.set_xticks(range(len(rhos)))
ax.set_xticklabels(pair_labels, fontsize=7, rotation=45, ha='right')
ax.set_ylabel('Spearman rho', fontsize=10)
ax.set_title('Expression Correlation within TCS Pairs', fontsize=11)
ax.axhline(0, color='gray', linewidth=0.5)
ax.axhline(0.6, color='red', linewidth=0.5, linestyle='--', alpha=0.5, label='rho=0.6')

# Legend for methylated partner
legend_elements = [
    mpatches.Patch(facecolor='#E74C3C', label='Methylated: SK'),
    mpatches.Patch(facecolor='#3498DB', label='Methylated: RR'),
    mpatches.Patch(facecolor='#95A5A6', label='Methylated: both/neither'),
]
ax.legend(handles=legend_elements, fontsize=8)

# Panel B: Coordination asymmetry
ax = axes[1]
asym_data = []
for _, row in tcs_df.iterrows():
    asym_data.append({
        'pair': f"{row['sk_old_locus']}/{row['rr_old_locus']}",
        'SK_T2': row['sk_coordination_T2'],
        'SK_T3': row['sk_coordination_T3'],
        'RR_T2': row['rr_coordination_T2'],
        'RR_T3': row['rr_coordination_T3']
    })

# Show as text table on the plot
ax.axis('off')
cell_text = []
for row in asym_data:
    cell_text.append([row['pair'], row['SK_T2'][:15], row['SK_T3'][:15],
                      row['RR_T2'][:15], row['RR_T3'][:15]])
table = ax.table(cellText=cell_text,
                 colLabels=['TCS Pair', 'SK coord T2', 'SK coord T3', 'RR coord T2', 'RR coord T3'],
                 loc='center', cellLoc='center')
table.auto_set_font_size(False)
table.set_fontsize(7)
table.scale(1, 1.5)
ax.set_title('Coordination Types per TCS Pair', fontsize=11)

plt.tight_layout()
for fmt in ['pdf', 'svg']:
    fig.savefig(f"{FIG_DIR}/TCS_pairs.{fmt}", dpi=300, bbox_inches='tight')
plt.close()

# ===== Figure 5: Expression heatmap =====
print("  Creating expression_heatmap figure...")

# Prepare annotation data
tf_fam_vec = annot['tf_family'].values
region_vec = annot['region_calc'].values

# Create color maps for annotations
tf_fam_colors_map = {fam: family_colors.get(fam, '#BDC3C7') for fam in set(tf_fam_vec)}
region_colors_map = region_colors

row_colors_df = pd.DataFrame({
    'TF family': [tf_fam_colors_map.get(f, '#BDC3C7') for f in tf_fam_vec],
    'Coordination T3': [coord_colors.get(c, '#BDC3C7') for c in coord_type_vec],
    'Region': [region_colors_map.get(r, '#BDC3C7') for r in region_vec]
}, index=expr_z.index)

g2 = sns.clustermap(expr_z, row_linkage=linkage_mat, col_cluster=False,
                    cmap='RdBu_r', center=0, vmin=-2.5, vmax=2.5,
                    figsize=(10, 16),
                    row_colors=row_colors_df,
                    xticklabels=True, yticklabels=True,
                    dendrogram_ratio=0.1,
                    cbar_kws={'label': 'Z-score'})
g2.ax_heatmap.set_xticklabels(g2.ax_heatmap.get_xticklabels(), fontsize=8, rotation=45, ha='right')
g2.ax_heatmap.set_yticklabels(g2.ax_heatmap.get_yticklabels(), fontsize=5)
g2.fig.suptitle('Expression Dynamics: 57 Exposed Regulatory Genes\n(Z-scored normalized counts, 9 samples)',
                fontsize=13, y=1.01)

# Add annotation legends manually below the heatmap
# We'll add small text annotations
for fmt in ['pdf', 'svg']:
    g2.savefig(f"{FIG_DIR}/expression_heatmap.{fmt}", dpi=300, bbox_inches='tight')
plt.close('all')

# ===== Figure 6: Comprehensive summary =====
print("  Creating H32_comprehensive_summary figure...")
fig = plt.figure(figsize=(18, 14))
gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.3)

# Panel A: Genomic distribution (simplified)
ax = fig.add_subplot(gs[0, :])
ax.barh(0, CHROM_LEN, height=0.3, color='#ECF0F1', edgecolor='#95A5A6', linewidth=1)
for boundary in [ARM_BOUNDARY, CHROM_LEN - ARM_BOUNDARY]:
    ax.axvline(boundary, color='#7F8C8D', linestyle='--', linewidth=0.8, alpha=0.5)
for _, gene in exp_sorted.iterrows():
    color = coord_colors.get(gene.get('coordination_T3', 'unknown'), '#BDC3C7')
    ax.plot(gene['midpoint'], 0, 'v', color=color, markersize=6, alpha=0.8)
for _, cl in cluster_df.iterrows():
    rect = Rectangle((cl['start'] - 5000, -0.25), cl['end'] - cl['start'] + 10000, 0.5,
                     linewidth=1.5, edgecolor='red', facecolor='red', alpha=0.15)
    ax.add_patch(rect)
ax.set_xlim(-100000, CHROM_LEN + 100000)
ax.set_ylim(-0.5, 0.6)
ax.set_yticks([])
ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x/1e6:.1f}'))
ax.set_xlabel('Position (Mb)', fontsize=9)
ax.set_title('A. Genomic Distribution of 57 Exposed Regulators', fontsize=11, fontweight='bold')

# Panel B: Permutation test
ax = fig.add_subplot(gs[1, 0])
ax.hist(perm_medians / 1e3, bins=25, color='#95A5A6', alpha=0.7, edgecolor='white')
ax.axvline(exp_median_nn / 1e3, color='#E74C3C', linewidth=2)
ax.set_xlabel('Median NN distance (kb)', fontsize=9)
ax.set_ylabel('Count', fontsize=9)
ax.set_title('B. Permutation Test\n(Genomic Clustering)', fontsize=10, fontweight='bold')
ax.text(0.95, 0.95, f'z={z_score:.2f}\np={min(p_closer, p_farther):.3f}',
        transform=ax.transAxes, ha='right', va='top', fontsize=9,
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

# Panel C: Correlation distributions
ax = fig.add_subplot(gs[1, 1])
ax.boxplot([ee_corrs, es_corrs, ss_corrs],
           labels=['E-E', 'E-S', 'S-S'],
           patch_artist=True, widths=0.5,
           boxprops=dict(facecolor='#3498DB', alpha=0.7))
ax.set_ylabel('Spearman rho', fontsize=9)
ax.set_title('C. Co-expression Comparison', fontsize=10, fontweight='bold')
sig_marker_ee = '***' if p_ee_es < 0.001 else ('**' if p_ee_es < 0.01 else ('*' if p_ee_es < 0.05 else 'ns'))
ax.text(0.95, 0.95, f'E-E vs E-S: {sig_marker_ee}\np={p_ee_es:.2e}',
        transform=ax.transAxes, ha='right', va='top', fontsize=8,
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

# Panel D: Coordination type distribution
ax = fig.add_subplot(gs[1, 2])
ct_t3_main = ct_t3.head(5)
bars = ax.barh(range(len(ct_t3_main)), ct_t3_main.values,
               color=[coord_colors.get(c, '#BDC3C7') for c in ct_t3_main.index])
ax.set_yticks(range(len(ct_t3_main)))
ax.set_yticklabels([c[:20] for c in ct_t3_main.index], fontsize=7)
ax.set_xlabel('Count', fontsize=9)
ax.set_title('D. Coordination Type (T3)', fontsize=10, fontweight='bold')

# Panel E: TCS pair correlations
ax = fig.add_subplot(gs[2, 0])
bars = ax.bar(range(len(tcs_df)), tcs_df['expression_rho'].values,
              color=[colors_tcs[i] for i in range(len(tcs_df))], width=0.6)
ax.set_xticks(range(len(tcs_df)))
ax.set_xticklabels([f"P{i+1}" for i in range(len(tcs_df))], fontsize=8)
ax.set_ylabel('Spearman rho', fontsize=9)
ax.axhline(0, color='gray', linewidth=0.5)
ax.set_title('E. TCS Pair Correlations', fontsize=10, fontweight='bold')

# Panel F: Module evidence summary
ax = fig.add_subplot(gs[2, 1:])
summary_text = (
    f"Module Detection Summary\n"
    f"{'='*40}\n\n"
    f"Genomic clustering (20kb): {len(multi_clusters)} clusters, {n_clustered} genes\n"
    f"Permutation test: z={z_score:.2f}, p={min(p_closer, p_farther):.3f}\n"
    f"Co-expression modules (rho>0.7): {len(coexpr_modules)}\n"
    f"Significant high-corr pairs: {n_sig_pairs}/{n_total_pairs}\n"
    f"Operonic pairs: {len(operon_pairs)} (expected: {np.median(random_operon_counts):.1f})\n"
    f"Functional modules detected: {len(module_df)}\n"
    f"  High-confidence (score>=2, >=3 genes): {len(hc_modules)}\n\n"
    f"E-E median rho: {np.median(ee_corrs):.3f}\n"
    f"E-S median rho: {np.median(es_corrs):.3f}\n"
    f"S-S median rho: {np.median(ss_corrs):.3f}\n"
    f"Within-type vs between-type: p={p_wt:.3e}\n\n"
    f"TCS pairs: {len(tcs_df)}, mean rho={tcs_df['expression_rho'].mean():.3f}\n"
    f"SK methylated: {n_sk_methyl}, RR methylated: {n_rr_methyl}"
)
ax.text(0.05, 0.95, summary_text, transform=ax.transAxes, fontsize=9,
        verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))
ax.axis('off')
ax.set_title('F. Summary Statistics', fontsize=10, fontweight='bold')

fig.suptitle('H32: Exposed Regulatory Module Analysis', fontsize=14, fontweight='bold', y=0.98)

for fmt in ['pdf', 'svg']:
    fig.savefig(f"{FIG_DIR}/H32_comprehensive_summary.{fmt}", dpi=300, bbox_inches='tight')
plt.close()

print("  All figures saved.")

# ============================================================
# STEP 10: Save output tables
# ============================================================
print("\n" + "=" * 70)
print("STEP 10: Saving output tables")
print("=" * 70)

# Table 1: Genomic clusters
cluster_df.to_csv(f"{TAB_DIR}/genomic_clusters.tsv", sep='\t', index=False)
print(f"  genomic_clusters.tsv: {len(cluster_df)} clusters")

# Table 2: Co-expression matrix
corr_exp.to_csv(f"{TAB_DIR}/coexpression_matrix.tsv", sep='\t')
print(f"  coexpression_matrix.tsv: {corr_exp.shape}")

# Table 3: Co-expression modules
coexpr_modules_df.to_csv(f"{TAB_DIR}/coexpression_modules.tsv", sep='\t', index=False)
print(f"  coexpression_modules.tsv: {len(coexpr_modules_df)} modules")

# Table 4: TCS pairs analysis
tcs_df.to_csv(f"{TAB_DIR}/TCS_pairs_analysis.tsv", sep='\t', index=False)
print(f"  TCS_pairs_analysis.tsv: {len(tcs_df)} pairs")

# Table 5: Coordination type distribution
coord_dist_df.to_csv(f"{TAB_DIR}/coordination_type_distribution.tsv", sep='\t', index=False)
print(f"  coordination_type_distribution.tsv: {len(coord_dist_df)} entries")

# Table 6: Statistical tests
stat_tests = []
stat_tests.append({
    'test': 'Genomic clustering permutation',
    'description': 'Exposed TFs closer than random 57/1055?',
    'test_statistic': f'z={z_score:.3f}',
    'observed': f'{exp_median_nn:.0f} bp',
    'expected': f'{np.mean(perm_medians):.0f} bp',
    'p_value': min(p_closer, p_farther),
    'direction': 'clustered' if p_closer < 0.05 else ('dispersed' if p_farther < 0.05 else 'not significant'),
    'effect_size': f'z={z_score:.3f}',
    'correction': 'permutation (n=1000)'
})
stat_tests.append({
    'test': 'E-E vs E-S co-expression',
    'description': 'Exposed-exposed correlations higher than exposed-shielded?',
    'test_statistic': f'U={u_ee_es:.0f}',
    'observed': f'median rho={np.median(ee_corrs):.4f}',
    'expected': f'median rho={np.median(es_corrs):.4f}',
    'p_value': p_ee_es,
    'direction': 'higher' if p_ee_es < 0.05 else 'not significant',
    'effect_size': f'rank-biserial r={r_ee_es:.4f}',
    'correction': 'none (one-sided MWU)'
})
stat_tests.append({
    'test': 'E-E vs S-S co-expression',
    'description': 'Exposed-exposed correlations higher than shielded-shielded?',
    'test_statistic': f'U={u_ee_ss:.0f}',
    'observed': f'median rho={np.median(ee_corrs):.4f}',
    'expected': f'median rho={np.median(ss_corrs):.4f}',
    'p_value': p_ee_ss,
    'direction': 'higher' if p_ee_ss < 0.05 else 'not significant',
    'effect_size': f'rank-biserial r={r_ee_ss:.4f}',
    'correction': 'none (one-sided MWU)'
})
stat_tests.append({
    'test': 'Coordination type vs cluster membership (chi2)',
    'description': 'Association between coordination_T3 and genomic cluster',
    'test_statistic': f'chi2={chi2_ct_cluster:.3f}' if not np.isnan(chi2_ct_cluster) else 'N/A',
    'observed': 'contingency table',
    'expected': 'independence',
    'p_value': p_ct_cluster if not np.isnan(p_ct_cluster) else np.nan,
    'direction': 'associated' if (not np.isnan(p_ct_cluster) and p_ct_cluster < 0.05) else 'not significant',
    'effect_size': f"Cramer's V={cramers_v:.4f}" if not np.isnan(cramers_v) else 'N/A',
    'correction': 'none'
})
stat_tests.append({
    'test': 'Within-type vs between-type co-expression',
    'description': 'Same coordination type = higher co-expression?',
    'test_statistic': f'U={u_wt:.0f}' if not np.isnan(u_wt) else 'N/A',
    'observed': f'within median={np.median(within_type_corrs):.4f}' if len(within_type_corrs) > 0 else 'N/A',
    'expected': f'between median={np.median(between_type_corrs):.4f}' if len(between_type_corrs) > 0 else 'N/A',
    'p_value': p_wt if not np.isnan(p_wt) else np.nan,
    'direction': 'higher' if (not np.isnan(p_wt) and p_wt < 0.05) else 'not significant',
    'effect_size': f'rank-biserial r={r_wt:.4f}' if not np.isnan(r_wt) else 'N/A',
    'correction': 'none (one-sided MWU)'
})
stat_tests.append({
    'test': 'Operon pair enrichment',
    'description': 'More exposed-exposed operonic pairs than random?',
    'test_statistic': f'observed={len(operon_pairs)}',
    'observed': f'{len(operon_pairs)} pairs',
    'expected': f'median={np.median(random_operon_counts):.1f}',
    'p_value': p_operon,
    'direction': 'enriched' if p_operon < 0.05 else 'not significant',
    'effect_size': f'fold={len(operon_pairs)/max(np.mean(random_operon_counts), 0.001):.2f}',
    'correction': 'permutation (n=1000)'
})
if not np.isnan(p_binom):
    stat_tests.append({
        'test': 'TCS methylation asymmetry (binomial)',
        'description': 'Is methylation consistently on SK or RR?',
        'test_statistic': f'SK={n_sk_methyl}, RR={n_rr_methyl}',
        'observed': f'{max(n_sk_methyl, n_rr_methyl)}/{n_asym}',
        'expected': '0.5',
        'p_value': p_binom,
        'direction': 'asymmetric' if p_binom < 0.05 else 'not significant',
        'effect_size': f'SK:RR ratio={n_sk_methyl}:{n_rr_methyl}',
        'correction': 'exact binomial'
    })

stat_df = pd.DataFrame(stat_tests)
stat_df.to_csv(f"{TAB_DIR}/statistical_tests.tsv", sep='\t', index=False)
print(f"  statistical_tests.tsv: {len(stat_df)} tests")

# Table 7: Module summary
module_df.to_csv(f"{TAB_DIR}/module_summary.tsv", sep='\t', index=False)
print(f"  module_summary.tsv: {len(module_df)} modules")

# ============================================================
# Final Summary
# ============================================================
print("\n" + "=" * 70)
print("ANALYSIS COMPLETE - H32 Summary")
print("=" * 70)

# Determine verdict
evidence_for = 0
evidence_against = 0

# Genomic clustering
if min(p_closer, p_farther) < 0.05:
    if p_closer < 0.05:
        evidence_for += 1
        print("[+] Genomic clustering: exposed TFs are CLOSER than random")
    else:
        evidence_against += 1
        print("[-] Genomic clustering: exposed TFs are MORE DISPERSED than random")
else:
    print("[=] Genomic clustering: NOT significantly different from random")

# Co-expression
if p_ee_es < 0.05:
    evidence_for += 1
    print("[+] Co-expression: E-E correlations HIGHER than E-S")
else:
    evidence_against += 1
    print("[-] Co-expression: E-E correlations NOT higher than E-S")

if p_ee_ss < 0.05:
    evidence_for += 1
    print("[+] Co-expression: E-E correlations HIGHER than S-S")
else:
    evidence_against += 1
    print("[-] Co-expression: E-E correlations NOT higher than S-S")

# Within-type coherence
if not np.isnan(p_wt) and p_wt < 0.05:
    evidence_for += 1
    print("[+] Coordination coherence: same-type TFs show HIGHER co-expression")
else:
    evidence_against += 1
    print("[-] Coordination coherence: same-type TFs do NOT show higher co-expression")

# Co-expression modules
if len(coexpr_modules) > 0:
    evidence_for += 1
    print(f"[+] Co-expression modules: {len(coexpr_modules)} modules detected")
else:
    evidence_against += 1
    print("[-] Co-expression modules: NONE detected")

# Operonic pairs
if p_operon < 0.05 and len(operon_pairs) > np.mean(random_operon_counts):
    evidence_for += 1
    print("[+] Operon context: MORE operonic pairs than expected")
else:
    print("[=] Operon context: NOT enriched for operonic pairs")

# High-confidence modules
if len(hc_modules) > 0:
    evidence_for += 1
    print(f"[+] High-confidence modules: {len(hc_modules)} detected")
else:
    evidence_against += 1
    print("[-] High-confidence modules: NONE detected")

print(f"\nEvidence FOR module: {evidence_for}")
print(f"Evidence AGAINST module: {evidence_against}")

if evidence_for >= 4:
    verdict = "SUPPORTED"
elif evidence_for >= 2:
    verdict = "PARTIAL"
else:
    verdict = "REJECTED"

print(f"\nVERDICT: {verdict}")
print(f"\nAll outputs saved to: {ANALYSIS_DIR}")
