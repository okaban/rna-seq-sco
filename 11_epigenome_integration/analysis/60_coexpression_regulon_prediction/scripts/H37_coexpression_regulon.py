#!/usr/bin/env python3
"""
H37: Genome-wide co-expression regulon prediction for 62 exposed TFs
of Streptomyces coelicolor M145.

VECTORIZED VERSION - uses rank-based approach for fast Spearman computation.
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import spearmanr, fisher_exact, mannwhitneyu, rankdata
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib_venn import venn2
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# PATHS
# ============================================================
BASE = '/Users/okaban/bioinfo/rna-seq'
DESEQ_DIR = f'{BASE}/04_deseq2/analysis/04_deseq2_260128_v1/results'
OUT_DIR = f'{BASE}/11_epigenome_integration/analysis/60_coexpression_regulon_prediction'
TABLES = f'{OUT_DIR}/tables'
FIGURES = f'{OUT_DIR}/figures'

COUNTS_FILE = f'{DESEQ_DIR}/normalized_counts_M145.tsv'
EXPOSED_FILE = f'{BASE}/11_epigenome_integration/analysis/51_exposed_regulators_characteristics/tables/exposed_regulators_full_table.tsv'
MODULES_FILE = f'{BASE}/11_epigenome_integration/analysis/55_exposed_regulatory_module/tables/coexpression_modules.tsv'
ANNOT_FILE = f'{BASE}/05_annotation/analysis/05_annotation_260128_v1/tables/gene_annotation_basic.tsv'
ALL_GENES_FILE = f'{BASE}/11_epigenome_integration/analysis/52_shielded_exposed_boundary/tables/all_genes_features.tsv'

DEG_FILES = {
    'T2vsT1': f'{DESEQ_DIR}/DESeq2_M145_2_vs_1.tsv',
    'T3vsT1': f'{DESEQ_DIR}/DESeq2_M145_3_vs_1.tsv',
    'T3vsT2': f'{DESEQ_DIR}/DESeq2_M145_3_vs_2.tsv',
}

np.random.seed(42)

# ============================================================
# HELPER: Vectorized Spearman correlation with p-values
# ============================================================
def fast_spearman_matrix(X, y):
    """
    Compute Spearman correlation between each row of X (n_genes x n_samples) 
    and vector y (n_samples,). Returns rho and p-value arrays.
    Uses rank transform + Pearson for speed.
    """
    n = X.shape[1]
    # Rank each row of X
    X_ranked = np.apply_along_axis(rankdata, 1, X)
    y_ranked = rankdata(y)
    
    # Center
    X_centered = X_ranked - X_ranked.mean(axis=1, keepdims=True)
    y_centered = y_ranked - y_ranked.mean()
    
    # Pearson on ranks = Spearman
    num = X_centered @ y_centered
    denom_x = np.sqrt((X_centered ** 2).sum(axis=1))
    denom_y = np.sqrt((y_centered ** 2).sum())
    
    rho = num / (denom_x * denom_y + 1e-30)
    rho = np.clip(rho, -1, 1)
    
    # P-value using t-distribution approximation
    t_stat = rho * np.sqrt((n - 2) / (1 - rho**2 + 1e-30))
    p_values = 2 * stats.t.sf(np.abs(t_stat), df=n-2)
    
    return rho, p_values


def fast_spearman_single(x, y):
    """Single pair Spearman with p-value."""
    r, p = spearmanr(x, y)
    return r, p

# ============================================================
# LOAD DATA
# ============================================================
print("=" * 70)
print("H37: Co-expression Regulon Prediction for 62 Exposed TFs")
print("=" * 70)

print("\n[Loading data]")
counts = pd.read_csv(COUNTS_FILE, sep='\t', index_col=0)
print(f"  Normalized counts: {counts.shape[0]} genes x {counts.shape[1]} samples")

exposed = pd.read_csv(EXPOSED_FILE, sep='\t')
print(f"  Exposed TFs: {len(exposed)} genes")

modules = pd.read_csv(MODULES_FILE, sep='\t')
print(f"  Co-expression modules: {len(modules)} modules")

annot = pd.read_csv(ANNOT_FILE, sep='\t')
print(f"  Gene annotations: {len(annot)} genes")

all_genes = pd.read_csv(ALL_GENES_FILE, sep='\t')
print(f"  All regulatory genes features: {len(all_genes)} genes")

degs = {}
for name, path in DEG_FILES.items():
    df = pd.read_csv(path, sep='\t', index_col=0)
    degs[name] = df
    n_sig = (df['padj'] < 0.05).sum()
    print(f"  DEGs {name}: {n_sig} significant (padj < 0.05)")

# ============================================================
# PARSE MODULE MEMBERSHIP
# ============================================================
activation_tags = []
repression_tags = []
module_membership = {}

for _, row in modules.iterrows():
    mid = row['module_id']
    tags = [t.strip() for t in row['locus_tags'].split(',')]
    for t in tags:
        module_membership[t] = mid
    if mid in [1, 2, 3]:
        activation_tags.extend(tags)
    elif mid == 4:
        repression_tags.extend(tags)

exposed_tags = set(exposed['locus_tag'].values)
assigned_tags = set(activation_tags + repression_tags)
unassigned = exposed_tags - assigned_tags
if unassigned:
    print(f"\n  Unassigned exposed TFs (not in modules 1-4): {len(unassigned)}")
    for t in sorted(unassigned):
        print(f"    {t}")

print(f"\n  Activation bloc (modules 1-3): {len(activation_tags)} TFs")
print(f"  Repression bloc (module 4): {len(repression_tags)} TFs")

activation_tags = [t for t in activation_tags if t in counts.index]
repression_tags = [t for t in repression_tags if t in counts.index]
all_exposed_in_counts = [t for t in exposed_tags if t in counts.index]
print(f"  Activation bloc in counts: {len(activation_tags)}")
print(f"  Repression bloc in counts: {len(repression_tags)}")
print(f"  Total exposed in counts: {len(all_exposed_in_counts)}")

# ============================================================
# STEP 1: Genome-wide Spearman correlations (VECTORIZED)
# ============================================================
print("\n" + "=" * 70)
print("STEP 1: Genome-wide Spearman correlations (vectorized)")
print("=" * 70)

# Filter expressed genes
gene_sums = counts.sum(axis=1)
# Need at least some variation (std > 0)
gene_stds = counts.std(axis=1)
expressed_mask = (gene_sums > 0) & (gene_stds > 0)
expressed_genes = counts.index[expressed_mask].tolist()
counts_expr = counts.loc[expressed_genes]
n_genes = len(expressed_genes)
print(f"  Expressed genes with variation: {n_genes}")

# Z-score normalize for eigengene
counts_z = counts_expr.subtract(counts_expr.mean(axis=1), axis=0).divide(counts_expr.std(axis=1), axis=0)

# Compute bloc eigengenes
act_in_z = [t for t in activation_tags if t in counts_z.index]
rep_in_z = [t for t in repression_tags if t in counts_z.index]
act_eigengene = counts_z.loc[act_in_z].median(axis=0).values
rep_eigengene = counts_z.loc[rep_in_z].median(axis=0).values

print(f"  Activation eigengene from {len(act_in_z)} genes")
print(f"  Repression eigengene from {len(rep_in_z)} genes")

# Convert counts to numpy matrix for vectorized computation
X = counts_expr.values  # (n_genes, n_samples)

# Eigengene correlations (vectorized)
print("  Computing eigengene correlations (vectorized)...")
rho_act, pval_act = fast_spearman_matrix(X, act_eigengene)
rho_rep, pval_rep = fast_spearman_matrix(X, rep_eigengene)

results_eigen = pd.DataFrame(index=expressed_genes)
results_eigen['rho_activation'] = rho_act
results_eigen['pval_activation'] = pval_act
results_eigen['rho_repression'] = rho_rep
results_eigen['pval_repression'] = pval_rep
results_eigen['padj_activation'] = np.minimum(pval_act * n_genes, 1.0)
results_eigen['padj_repression'] = np.minimum(pval_rep * n_genes, 1.0)

print(f"  Eigengene correlation summary:")
print(f"    Activation: median rho = {np.median(rho_act):.3f}, "
      f"range [{np.min(rho_act):.3f}, {np.max(rho_act):.3f}]")
print(f"    Repression: median rho = {np.median(rho_rep):.3f}, "
      f"range [{np.min(rho_rep):.3f}, {np.max(rho_rep):.3f}]")

# Per-TF correlations (vectorized)
print("  Computing per-TF correlations (vectorized)...")
per_tf_rhos = {}
per_tf_pvals = {}

for tf in all_exposed_in_counts:
    tf_vals = counts_expr.loc[tf].values
    r, p = fast_spearman_matrix(X, tf_vals)
    per_tf_rhos[tf] = r
    per_tf_pvals[tf] = p

print(f"  Per-TF correlations computed for {len(per_tf_rhos)} TFs")

# ============================================================
# STEP 2: Define co-expressed gene sets
# ============================================================
print("\n" + "=" * 70)
print("STEP 2: Define co-expressed gene sets")
print("=" * 70)

thresholds = [0.8, 0.7]
gene_sets = {}

for thresh in thresholds:
    act_pos = set(results_eigen.index[
        (results_eigen['rho_activation'] > thresh) &
        (results_eigen['padj_activation'] < 0.05)
    ])
    act_neg = set(results_eigen.index[
        (results_eigen['rho_activation'] < -thresh) &
        (results_eigen['padj_activation'] < 0.05)
    ])
    rep_pos = set(results_eigen.index[
        (results_eigen['rho_repression'] > thresh) &
        (results_eigen['padj_repression'] < 0.05)
    ])
    rep_neg = set(results_eigen.index[
        (results_eigen['rho_repression'] < -thresh) &
        (results_eigen['padj_repression'] < 0.05)
    ])

    gene_sets[thresh] = {
        'act_pos': act_pos, 'act_neg': act_neg,
        'rep_pos': rep_pos, 'rep_neg': rep_neg
    }

    act_all = act_pos | act_neg
    rep_all = rep_pos | rep_neg
    both = act_all & rep_all
    neither = set(expressed_genes) - act_all - rep_all

    print(f"\n  Threshold |rho| > {thresh} (Bonferroni p < 0.05):")
    print(f"    Activation-correlated (pos): {len(act_pos)}")
    print(f"    Activation-correlated (neg): {len(act_neg)}")
    print(f"    Repression-correlated (pos): {len(rep_pos)}")
    print(f"    Repression-correlated (neg): {len(rep_neg)}")
    print(f"    Both blocs (any): {len(both)}")
    print(f"    Neither: {len(neither)}")

primary_thresh = 0.8
act_regulon_pos = gene_sets[primary_thresh]['act_pos']
act_regulon_neg = gene_sets[primary_thresh]['act_neg']
rep_regulon_pos = gene_sets[primary_thresh]['rep_pos']
rep_regulon_neg = gene_sets[primary_thresh]['rep_neg']

# ============================================================
# STEP 3: Functional enrichment
# ============================================================
print("\n" + "=" * 70)
print("STEP 3: Functional enrichment (keyword-based)")
print("=" * 70)

annot_dict = {}
for _, row in annot.iterrows():
    annot_dict[row['gene_id']] = str(row.get('product', ''))

results_eigen['product'] = [annot_dict.get(g, '') for g in results_eigen.index]

keywords = [
    'regulat', 'transport', 'hydrolase', 'oxidoreductase', 'transferase',
    'synthase', 'ribosom', 'hypothetical', 'ABC', 'cytochrome', 'sigma',
    'kinase', 'secreted', 'membrane', 'peptidase', 'dehydrogenase',
    'reductase', 'lyase', 'polyketide', 'peptide synthetase', 'permease',
    'acyl-CoA', 'methyltransferase', 'helix-turn-helix', 'response regulator',
    'sensor', 'two-component'
]

gene_set_names_all = {
    'act_pos_0.8': act_regulon_pos,
    'act_neg_0.8': act_regulon_neg,
    'rep_pos_0.8': rep_regulon_pos,
    'rep_neg_0.8': rep_regulon_neg,
    'act_pos_0.7': gene_sets[0.7]['act_pos'],
    'act_neg_0.7': gene_sets[0.7]['act_neg'],
    'rep_pos_0.7': gene_sets[0.7]['rep_pos'],
    'rep_neg_0.7': gene_sets[0.7]['rep_neg'],
}

all_products = {g: annot_dict.get(g, '').lower() for g in expressed_genes}
background = set(expressed_genes)

enrichment_results = []
for set_name, gene_set in gene_set_names_all.items():
    if len(gene_set) == 0:
        continue
    for kw in keywords:
        kw_lower = kw.lower()
        in_set_with_kw = sum(1 for g in gene_set if kw_lower in all_products.get(g, ''))
        in_set_without_kw = len(gene_set) - in_set_with_kw
        bg = background - gene_set
        in_bg_with_kw = sum(1 for g in bg if kw_lower in all_products.get(g, ''))
        in_bg_without_kw = len(bg) - in_bg_with_kw
        if in_set_with_kw == 0 and in_bg_with_kw == 0:
            continue
        table = [[in_set_with_kw, in_set_without_kw], [in_bg_with_kw, in_bg_without_kw]]
        odds_ratio, pval = fisher_exact(table, alternative='two-sided')
        total_kw = in_set_with_kw + in_bg_with_kw
        expected = total_kw * len(gene_set) / len(background) if len(background) > 0 else 0
        fold = in_set_with_kw / expected if expected > 0 else np.inf
        enrichment_results.append({
            'gene_set': set_name, 'keyword': kw,
            'observed': in_set_with_kw, 'expected': round(expected, 1),
            'fold_enrichment': round(fold, 2),
            'odds_ratio': round(odds_ratio, 3) if not np.isinf(odds_ratio) else 'inf',
            'pvalue': pval, 'set_size': len(gene_set), 'bg_with_keyword': in_bg_with_kw,
        })

enrich_df = pd.DataFrame(enrichment_results)
if len(enrich_df) > 0:
    enrich_df['padj'] = np.nan
    for sn in enrich_df['gene_set'].unique():
        mask = enrich_df['gene_set'] == sn
        pvals = enrich_df.loc[mask, 'pvalue'].values
        _, padj, _, _ = multipletests(pvals, method='fdr_bh')
        enrich_df.loc[mask, 'padj'] = padj

    sig_enrich = enrich_df[enrich_df['padj'] < 0.1].sort_values(['gene_set', 'padj'])
    print(f"\n  Significant enrichments (FDR < 0.1): {len(sig_enrich)}")
    for _, row in sig_enrich.iterrows():
        direction = 'ENRICHED' if row['fold_enrichment'] > 1 else 'DEPLETED'
        print(f"    {row['gene_set']} | {row['keyword']}: fold={row['fold_enrichment']}, "
              f"obs={row['observed']}, exp={row['expected']}, padj={row['padj']:.4f} [{direction}]")

# ============================================================
# STEP 4: Known Streptomyces pathway overlap
# ============================================================
print("\n" + "=" * 70)
print("STEP 4: Known Streptomyces pathway overlap")
print("=" * 70)

pathway_keywords = {
    'bld_bald': ['bald', 'bldA', 'bldB', 'bldC', 'bldD', 'bldG', 'bldH', 'bldK', 'bldM', 'bldN'],
    'whi_white': ['white', 'whiA', 'whiB', 'whiD', 'whiE', 'whiG', 'whiH', 'whiI', 'whiJ'],
    'ram': ['ram', 'SapB', 'amfS'],
    'sap': ['sap'],
    'chp_chaplins': ['chaplin', 'chpA', 'chpB', 'chpC', 'chpD', 'chpE', 'chpF', 'chpG', 'chpH'],
    'rdl_rodlins': ['rodlin', 'rdlA', 'rdlB'],
    'sigma_factors': ['sigma factor', 'sigma-70', 'ECF sigma', 'anti-sigma'],
    'secondary_metabolism': ['actinorhodin', 'undecylprodigiosin', 'prodiginine', 'CDA',
                              'calcium-dependent antibiotic', 'coelimycin', 'methylenomycin',
                              'polyketide synthase', 'nonribosomal peptide', 'NRPS',
                              'type I PKS', 'type II PKS', 'type III PKS'],
    'two_component': ['two-component', 'sensor histidine kinase', 'response regulator',
                       'sensor kinase', 'histidine kinase'],
    'ABC_transport': ['ABC transporter', 'ATP-binding cassette'],
}

annot_name_dict = {}
for _, row in annot.iterrows():
    annot_name_dict[row['gene_id']] = str(row.get('gene_name', ''))

pathway_genes = {}
for pathway, kw_list in pathway_keywords.items():
    matching = set()
    for gene in expressed_genes:
        prod = all_products.get(gene, '')
        gname = annot_name_dict.get(gene, '').lower()
        for kw in kw_list:
            if kw.lower() in prod or kw.lower() in gname:
                matching.add(gene)
                break
    pathway_genes[pathway] = matching

pathway_overlap_results = []
gene_set_test = {
    'act_pos_0.8': act_regulon_pos,
    'act_neg_0.8': act_regulon_neg,
    'rep_pos_0.8': rep_regulon_pos,
    'rep_neg_0.8': rep_regulon_neg,
}

for pathway, pw_genes in pathway_genes.items():
    if len(pw_genes) == 0:
        continue
    for set_name, gene_set in gene_set_test.items():
        overlap = pw_genes & gene_set
        a = len(overlap)
        b = len(gene_set - pw_genes)
        c = len(pw_genes - gene_set)
        d = len(background - pw_genes - gene_set)
        table = [[a, b], [c, d]]
        odds_ratio, pval = fisher_exact(table, alternative='two-sided')
        expected = len(pw_genes) * len(gene_set) / len(background) if len(background) > 0 else 0
        fold = a / expected if expected > 0 else np.inf
        overlap_genes_str = ', '.join(sorted(overlap)[:10])
        if len(overlap) > 10:
            overlap_genes_str += f' ... (+{len(overlap)-10} more)'
        pathway_overlap_results.append({
            'pathway': pathway, 'gene_set': set_name,
            'pathway_size': len(pw_genes), 'set_size': len(gene_set),
            'overlap': a, 'expected': round(expected, 1),
            'fold_enrichment': round(fold, 2),
            'odds_ratio': round(odds_ratio, 3) if not np.isinf(odds_ratio) else 'inf',
            'pvalue': pval, 'overlap_genes': overlap_genes_str,
        })

pathway_df = pd.DataFrame(pathway_overlap_results)
if len(pathway_df) > 0:
    _, padj, _, _ = multipletests(pathway_df['pvalue'].values, method='fdr_bh')
    pathway_df['padj'] = padj
    sig_pw = pathway_df[pathway_df['pvalue'] < 0.05].sort_values('pvalue')
    print(f"\n  Pathway overlaps (p < 0.05): {len(sig_pw)}")
    for _, row in sig_pw.iterrows():
        print(f"    {row['pathway']} x {row['gene_set']}: overlap={row['overlap']}/{row['pathway_size']}, "
              f"fold={row['fold_enrichment']}, p={row['pvalue']:.4f}")

# ============================================================
# STEP 5: Overlap with DEGs
# ============================================================
print("\n" + "=" * 70)
print("STEP 5: Overlap with DEGs")
print("=" * 70)

deg_overlap_results = []
for transition, deg_df in degs.items():
    sig_up = set(deg_df.index[(deg_df['padj'] < 0.05) & (deg_df['log2FoldChange'] > 0)])
    sig_down = set(deg_df.index[(deg_df['padj'] < 0.05) & (deg_df['log2FoldChange'] < 0)])
    sig_all = sig_up | sig_down
    for set_name, gene_set in gene_set_test.items():
        for deg_dir, deg_set in [('up', sig_up), ('down', sig_down), ('all', sig_all)]:
            overlap = gene_set & deg_set
            a = len(overlap)
            b = len(gene_set - deg_set)
            c = len(deg_set - gene_set)
            d = len(background - gene_set - deg_set)
            table = [[a, b], [c, d]]
            odds_ratio, pval = fisher_exact(table, alternative='two-sided')
            expected = len(deg_set) * len(gene_set) / len(background) if len(background) > 0 else 0
            fold = a / expected if expected > 0 else np.inf
            deg_overlap_results.append({
                'transition': transition, 'deg_direction': deg_dir,
                'gene_set': set_name, 'set_size': len(gene_set),
                'deg_size': len(deg_set), 'overlap': a,
                'expected': round(expected, 1), 'fold_enrichment': round(fold, 2),
                'odds_ratio': round(odds_ratio, 3) if not np.isinf(odds_ratio) else 'inf',
                'pvalue': pval,
            })

deg_overlap_df = pd.DataFrame(deg_overlap_results)
if len(deg_overlap_df) > 0:
    _, padj, _, _ = multipletests(deg_overlap_df['pvalue'].values, method='fdr_bh')
    deg_overlap_df['padj'] = padj
    sig_deg = deg_overlap_df[deg_overlap_df['padj'] < 0.05].sort_values('padj')
    print(f"\n  Significant DEG overlaps (FDR < 0.05): {len(sig_deg)}")
    for _, row in sig_deg.iterrows():
        print(f"    {row['gene_set']} x {row['transition']}_{row['deg_direction']}: "
              f"overlap={row['overlap']}, fold={row['fold_enrichment']}, padj={row['padj']:.2e}")

# ============================================================
# STEP 6: Per-TF regulon size estimation
# ============================================================
print("\n" + "=" * 70)
print("STEP 6: Per-TF regulon size estimation")
print("=" * 70)

regulon_sizes = []
for tf in all_exposed_in_counts:
    rhos = per_tf_rhos[tf]
    pvals = per_tf_pvals[tf]
    padj_vals = np.minimum(pvals * n_genes, 1.0)

    n_pos_08 = int(np.sum((rhos > 0.8) & (padj_vals < 0.05)))
    n_neg_08 = int(np.sum((rhos < -0.8) & (padj_vals < 0.05)))
    n_pos_07 = int(np.sum((rhos > 0.7) & (padj_vals < 0.05)))
    n_neg_07 = int(np.sum((rhos < -0.7) & (padj_vals < 0.05)))

    # Top 5 positively correlated targets (excluding self)
    gene_idx = {g: i for i, g in enumerate(expressed_genes)}
    sorted_idx = np.argsort(-rhos)
    top5 = []
    for j in sorted_idx:
        if expressed_genes[j] != tf and rhos[j] > 0.5:
            top5.append((expressed_genes[j], round(float(rhos[j]), 3)))
        if len(top5) >= 5:
            break

    bloc = 'activation' if tf in activation_tags else ('repression' if tf in repression_tags else 'unassigned')

    tf_info = exposed[exposed['locus_tag'] == tf]
    gene_name = str(tf_info['gene_name'].values[0]) if len(tf_info) > 0 and 'gene_name' in tf_info.columns else ''
    product = str(tf_info['product'].values[0]) if len(tf_info) > 0 else ''
    tf_family = str(tf_info['tf_family'].values[0]) if len(tf_info) > 0 else ''
    old_locus = str(tf_info['old_locus_tag'].values[0]) if len(tf_info) > 0 else ''

    regulon_sizes.append({
        'locus_tag': tf, 'old_locus_tag': old_locus,
        'gene_name': gene_name, 'product': product,
        'tf_family': tf_family, 'bloc': bloc,
        'n_correlated_0.8': n_pos_08, 'n_anticorrelated_0.8': n_neg_08,
        'n_correlated_0.7': n_pos_07, 'n_anticorrelated_0.7': n_neg_07,
        'total_regulon_0.8': n_pos_08 + n_neg_08,
        'total_regulon_0.7': n_pos_07 + n_neg_07,
        'top_targets': '; '.join([f"{g}({r})" for g, r in top5]),
    })

regulon_df = pd.DataFrame(regulon_sizes)

act_sizes = regulon_df[regulon_df['bloc'] == 'activation']['n_correlated_0.8'].values
rep_sizes = regulon_df[regulon_df['bloc'] == 'repression']['n_correlated_0.8'].values

print(f"\n  Activation bloc TFs: n={len(act_sizes)}")
print(f"    Regulon size (rho>0.8): median={np.median(act_sizes):.0f}, "
      f"mean={np.mean(act_sizes):.1f}, range=[{np.min(act_sizes)}, {np.max(act_sizes)}]")
print(f"  Repression bloc TFs: n={len(rep_sizes)}")
print(f"    Regulon size (rho>0.8): median={np.median(rep_sizes):.0f}, "
      f"mean={np.mean(rep_sizes):.1f}, range=[{np.min(rep_sizes)}, {np.max(rep_sizes)}]")

stat_mwu, pval_mwu = mannwhitneyu(act_sizes, rep_sizes, alternative='two-sided')
r_mwu = abs(stat_mwu - len(act_sizes)*len(rep_sizes)/2) / (len(act_sizes)*len(rep_sizes)/2)
print(f"  MWU (act vs rep): U={stat_mwu:.0f}, p={pval_mwu:.4f}, r={r_mwu:.3f}")

print(f"\n  Regulon >= 30 genes: {(regulon_df['n_correlated_0.8'] >= 30).sum()}/{len(regulon_df)}")
print(f"  Regulon >= 100 genes: {(regulon_df['n_correlated_0.8'] >= 100).sum()}/{len(regulon_df)}")
print(f"  Regulon == 0 genes: {(regulon_df['n_correlated_0.8'] == 0).sum()}/{len(regulon_df)}")

# ============================================================
# STEP 7: Cross-bloc exclusivity
# ============================================================
print("\n" + "=" * 70)
print("STEP 7: Cross-bloc exclusivity")
print("=" * 70)

r_cross, p_cross = spearmanr(results_eigen['rho_activation'].values,
                              results_eigen['rho_repression'].values)

for thresh in [0.8, 0.7]:
    gs = gene_sets[thresh]
    set_a = gs['act_pos']
    set_r = gs['rep_pos']
    intersection = set_a & set_r
    union = set_a | set_r
    jaccard = len(intersection) / len(union) if len(union) > 0 else 0

    print(f"\n  Threshold |rho| > {thresh}:")
    print(f"    Act-positive: {len(set_a)}, Rep-positive: {len(set_r)}")
    print(f"    Intersection: {len(intersection)}")
    print(f"    Jaccard index: {jaccard:.4f}")

    a = len(set_a & set_r)
    b = len(set_a - set_r)
    c = len(set_r - set_a)
    d = len(background - set_a - set_r)
    or_val, pval = fisher_exact([[a, b], [c, d]])
    print(f"    Fisher's exact: OR={or_val:.4f}, p={pval:.2e}")

    if len(set_a) > 1:
        act_pos_rep_rhos = results_eigen.loc[list(set_a), 'rho_repression'].values
        stat_w, pval_w = stats.wilcoxon(act_pos_rep_rhos)
        print(f"    Act-pos genes: median rho_repression = {np.median(act_pos_rep_rhos):.3f}, Wilcoxon p = {pval_w:.2e}")

    if len(set_r) > 1:
        rep_pos_act_rhos = results_eigen.loc[list(set_r), 'rho_activation'].values
        stat_w, pval_w = stats.wilcoxon(rep_pos_act_rhos)
        print(f"    Rep-pos genes: median rho_activation = {np.median(rep_pos_act_rhos):.3f}, Wilcoxon p = {pval_w:.2e}")

print(f"\n  Global cross-eigengene correlation: rho={r_cross:.3f}, p={p_cross:.2e}")

# ============================================================
# STEP 8: Statistical controls (permutation tests) - VECTORIZED
# ============================================================
print("\n" + "=" * 70)
print("STEP 8: Statistical controls (permutation tests)")
print("=" * 70)

all_reg_tags = all_genes['locus_tag'].values
all_reg_in_counts = [t for t in all_reg_tags if t in counts_expr.index]
print(f"  All regulatory genes in counts: {len(all_reg_in_counts)}")

n_perm = 1000
real_act_size_08 = len(act_regulon_pos)
real_rep_size_08 = len(rep_regulon_pos)

print(f"  Real activation regulon size (rho > 0.8): {real_act_size_08}")
print(f"  Real repression regulon size (rho > 0.8): {real_rep_size_08}")

# Pre-rank the count matrix once
X_ranked = np.apply_along_axis(rankdata, 1, X)
X_centered = X_ranked - X_ranked.mean(axis=1, keepdims=True)
X_denom = np.sqrt((X_centered ** 2).sum(axis=1))

# Pre-compute Z-scores for all regulatory genes for eigengene sampling
reg_indices = [expressed_genes.index(t) for t in all_reg_in_counts if t in expressed_genes]

print(f"  Running {n_perm} permutations (vectorized)...")

perm_act_sizes = np.zeros(n_perm)
perm_rep_sizes = np.zeros(n_perm)

n_samples = X.shape[1]

for perm_i in range(n_perm):
    if perm_i % 200 == 0:
        print(f"    Permutation {perm_i}/{n_perm}...")

    # Random activation bloc
    rand_act_idx = np.random.choice(reg_indices, size=min(len(activation_tags), len(reg_indices)), replace=False)
    rand_act_z = counts_z.iloc[rand_act_idx].values
    rand_act_eigen = np.median(rand_act_z, axis=0)
    
    # Rank and correlate
    y_ranked = rankdata(rand_act_eigen)
    y_centered = y_ranked - y_ranked.mean()
    y_denom = np.sqrt((y_centered ** 2).sum())
    
    rho_perm = (X_centered @ y_centered) / (X_denom * y_denom + 1e-30)
    rho_perm = np.clip(rho_perm, -1, 1)
    t_stat = rho_perm * np.sqrt((n_samples - 2) / (1 - rho_perm**2 + 1e-30))
    p_perm = 2 * stats.t.sf(np.abs(t_stat), df=n_samples-2)
    padj_perm = np.minimum(p_perm * n_genes, 1.0)
    
    perm_act_sizes[perm_i] = np.sum((rho_perm > 0.8) & (padj_perm < 0.05))

    # Random repression bloc
    rand_rep_idx = np.random.choice(reg_indices, size=min(len(repression_tags), len(reg_indices)), replace=False)
    rand_rep_z = counts_z.iloc[rand_rep_idx].values
    rand_rep_eigen = np.median(rand_rep_z, axis=0)
    
    y_ranked = rankdata(rand_rep_eigen)
    y_centered = y_ranked - y_ranked.mean()
    y_denom = np.sqrt((y_centered ** 2).sum())
    
    rho_perm = (X_centered @ y_centered) / (X_denom * y_denom + 1e-30)
    rho_perm = np.clip(rho_perm, -1, 1)
    t_stat = rho_perm * np.sqrt((n_samples - 2) / (1 - rho_perm**2 + 1e-30))
    p_perm = 2 * stats.t.sf(np.abs(t_stat), df=n_samples-2)
    padj_perm = np.minimum(p_perm * n_genes, 1.0)
    
    perm_rep_sizes[perm_i] = np.sum((rho_perm > 0.8) & (padj_perm < 0.05))

perm_p_act = (np.sum(perm_act_sizes >= real_act_size_08) + 1) / (n_perm + 1)
perm_p_rep = (np.sum(perm_rep_sizes >= real_rep_size_08) + 1) / (n_perm + 1)

act_z = (real_act_size_08 - np.mean(perm_act_sizes)) / (np.std(perm_act_sizes) + 1e-10)
rep_z = (real_rep_size_08 - np.mean(perm_rep_sizes)) / (np.std(perm_rep_sizes) + 1e-10)

print(f"\n  Permutation results ({n_perm} iterations):")
print(f"    Activation bloc:")
print(f"      Real: {real_act_size_08}, Random median: {np.median(perm_act_sizes):.0f}, "
      f"mean: {np.mean(perm_act_sizes):.1f}, 95th: {np.percentile(perm_act_sizes, 95):.0f}")
print(f"      p = {perm_p_act:.4f}, Z = {act_z:.2f}")
print(f"    Repression bloc:")
print(f"      Real: {real_rep_size_08}, Random median: {np.median(perm_rep_sizes):.0f}, "
      f"mean: {np.mean(perm_rep_sizes):.1f}, 95th: {np.percentile(perm_rep_sizes, 95):.0f}")
print(f"      p = {perm_p_rep:.4f}, Z = {rep_z:.2f}")

# ============================================================
# SAVE TABLES
# ============================================================
print("\n" + "=" * 70)
print("SAVING TABLES")
print("=" * 70)

regulon_df.to_csv(f'{TABLES}/per_TF_regulon.tsv', sep='\t', index=False)
print(f"  Saved: tables/per_TF_regulon.tsv ({len(regulon_df)} rows)")

results_eigen['in_act_regulon_0.8'] = results_eigen.index.isin(act_regulon_pos)
results_eigen['in_rep_regulon_0.8'] = results_eigen.index.isin(rep_regulon_pos)
results_eigen['in_act_regulon_0.7'] = results_eigen.index.isin(gene_sets[0.7]['act_pos'])
results_eigen['in_rep_regulon_0.7'] = results_eigen.index.isin(gene_sets[0.7]['rep_pos'])

annot_old = {}
annot_gname = {}
for _, row in annot.iterrows():
    annot_old[row['gene_id']] = row.get('old_locus_tag', '')
    annot_gname[row['gene_id']] = row.get('gene_name', '')

results_eigen['old_locus_tag'] = [annot_old.get(g, '') for g in results_eigen.index]
results_eigen['gene_name_annot'] = [annot_gname.get(g, '') for g in results_eigen.index]

col_order = ['old_locus_tag', 'gene_name_annot', 'product',
             'rho_activation', 'pval_activation', 'padj_activation',
             'rho_repression', 'pval_repression', 'padj_repression',
             'in_act_regulon_0.8', 'in_rep_regulon_0.8',
             'in_act_regulon_0.7', 'in_rep_regulon_0.7']
results_eigen_save = results_eigen[col_order].copy()
results_eigen_save.index.name = 'gene_id'
results_eigen_save.to_csv(f'{TABLES}/bloc_eigengene_correlations.tsv', sep='\t')
print(f"  Saved: tables/bloc_eigengene_correlations.tsv ({len(results_eigen_save)} rows)")

if len(enrich_df) > 0:
    enrich_df.to_csv(f'{TABLES}/functional_enrichment.tsv', sep='\t', index=False)
    print(f"  Saved: tables/functional_enrichment.tsv ({len(enrich_df)} rows)")

if len(pathway_df) > 0:
    pathway_df.to_csv(f'{TABLES}/pathway_overlap.tsv', sep='\t', index=False)
    print(f"  Saved: tables/pathway_overlap.tsv ({len(pathway_df)} rows)")

# Regulon size statistics
size_stats_rows = []
for thresh_col, thresh_label in [('n_correlated_0.8', 'regulon_pos_0.8'),
                                  ('n_anticorrelated_0.8', 'regulon_neg_0.8'),
                                  ('total_regulon_0.8', 'regulon_total_0.8'),
                                  ('n_correlated_0.7', 'regulon_pos_0.7'),
                                  ('total_regulon_0.7', 'regulon_total_0.7')]:
    av = regulon_df[regulon_df['bloc'] == 'activation'][thresh_col].values
    rv = regulon_df[regulon_df['bloc'] == 'repression'][thresh_col].values
    if len(av) > 0 and len(rv) > 0:
        s, p = mannwhitneyu(av, rv, alternative='two-sided')
        n_t = len(av) * len(rv)
        re = abs(s - n_t/2) / (n_t/2)
    else:
        s, p, re = np.nan, np.nan, np.nan
    size_stats_rows.append({
        'metric': thresh_label,
        'activation_bloc': f"median={np.median(av):.0f}, mean={np.mean(av):.1f}",
        'repression_bloc': f"median={np.median(rv):.0f}, mean={np.mean(rv):.1f}",
        'test_statistic': f"U={s:.0f}", 'p_value': p, 'effect_size': f"r={re:.3f}",
    })

pd.DataFrame(size_stats_rows).to_csv(f'{TABLES}/regulon_size_statistics.tsv', sep='\t', index=False)
print(f"  Saved: tables/regulon_size_statistics.tsv")

# Statistical tests
stat_tests = []
for thresh in [0.8, 0.7]:
    gs = gene_sets[thresh]
    set_a = gs['act_pos']
    set_r = gs['rep_pos']
    intersection = set_a & set_r
    union = set_a | set_r
    jac = len(intersection) / len(union) if len(union) > 0 else 0
    a_n = len(set_a & set_r)
    b_n = len(set_a - set_r)
    c_n = len(set_r - set_a)
    d_n = len(background - set_a - set_r)
    or_v, pv = fisher_exact([[a_n, b_n], [c_n, d_n]])
    stat_tests.append({
        'test': f'cross_bloc_exclusivity_rho{thresh}',
        'description': f'Fisher exact: act_pos vs rep_pos overlap at rho>{thresh}',
        'statistic': f'OR={or_v:.4f}', 'p_value': pv,
        'effect_size': f'Jaccard={jac:.4f}',
        'n1': len(set_a), 'n2': len(set_r),
        'interpretation': 'exclusive' if or_v < 0.5 else ('overlap' if or_v > 2 else 'independent'),
    })

stat_tests.append({
    'test': 'eigengene_cross_correlation',
    'description': 'Spearman: rho_activation vs rho_repression genome-wide',
    'statistic': f'rho={r_cross:.4f}', 'p_value': p_cross,
    'effect_size': f'rho={r_cross:.4f}', 'n1': n_genes, 'n2': n_genes,
    'interpretation': 'anti-correlated' if r_cross < -0.5 else ('correlated' if r_cross > 0.5 else 'weak'),
})

stat_tests.append({
    'test': 'regulon_size_bloc_comparison',
    'description': 'MWU: regulon size activation vs repression (rho>0.8)',
    'statistic': f'U={stat_mwu:.0f}', 'p_value': pval_mwu,
    'effect_size': f'act_median={np.median(act_sizes):.0f}, rep_median={np.median(rep_sizes):.0f}',
    'n1': len(act_sizes), 'n2': len(rep_sizes),
    'interpretation': 'different' if pval_mwu < 0.05 else 'not_different',
})

stat_tests.append({
    'test': 'permutation_activation_bloc',
    'description': f'Permutation ({n_perm}x): random {len(activation_tags)} reg genes',
    'statistic': f'real={real_act_size_08}, perm_mean={np.mean(perm_act_sizes):.1f}',
    'p_value': perm_p_act, 'effect_size': f'Z={act_z:.2f}',
    'n1': len(activation_tags), 'n2': n_perm,
    'interpretation': 'significant' if perm_p_act < 0.05 else 'not_significant',
})

stat_tests.append({
    'test': 'permutation_repression_bloc',
    'description': f'Permutation ({n_perm}x): random {len(repression_tags)} reg genes',
    'statistic': f'real={real_rep_size_08}, perm_mean={np.mean(perm_rep_sizes):.1f}',
    'p_value': perm_p_rep, 'effect_size': f'Z={rep_z:.2f}',
    'n1': len(repression_tags), 'n2': n_perm,
    'interpretation': 'significant' if perm_p_rep < 0.05 else 'not_significant',
})

pd.DataFrame(stat_tests).to_csv(f'{TABLES}/statistical_tests.tsv', sep='\t', index=False)
print(f"  Saved: tables/statistical_tests.tsv")

if len(deg_overlap_df) > 0:
    deg_overlap_df.to_csv(f'{TABLES}/deg_overlap.tsv', sep='\t', index=False)
    print(f"  Saved: tables/deg_overlap.tsv")

# ============================================================
# FIGURES
# ============================================================
print("\n" + "=" * 70)
print("GENERATING FIGURES")
print("=" * 70)

ACT_COLOR = '#2166AC'
REP_COLOR = '#B2182B'
BOTH_COLOR = '#7B3294'
NEUTRAL_COLOR = '#999999'

# Figure 1: Venn diagram
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for ax_i, thresh in enumerate([0.8, 0.7]):
    ax = axes[ax_i]
    gs = gene_sets[thresh]
    sa = gs['act_pos']; sr = gs['rep_pos']
    oa = len(sa - sr); orr = len(sr - sa); bn = len(sa & sr)
    if oa + orr + bn > 0:
        v = venn2(subsets=(oa, orr, bn),
                  set_labels=('Activation\nregulons', 'Repression\nregulons'), ax=ax)
        for pid, col in [('10', ACT_COLOR), ('01', REP_COLOR), ('11', BOTH_COLOR)]:
            if v.get_patch_by_id(pid):
                v.get_patch_by_id(pid).set_color(col)
                v.get_patch_by_id(pid).set_alpha(0.5)
    jac = len(sa & sr) / len(sa | sr) if len(sa | sr) > 0 else 0
    ax.set_title(f'|rho| > {thresh}\nJaccard = {jac:.3f}', fontsize=12, fontweight='bold')
fig.suptitle('Activation vs Repression Bloc Regulons', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
fig.savefig(f'{FIGURES}/bloc_regulon_venn.pdf', bbox_inches='tight', dpi=300)
fig.savefig(f'{FIGURES}/bloc_regulon_venn.svg', bbox_inches='tight', dpi=300)
plt.close(fig)
print("  Saved: figures/bloc_regulon_venn.pdf/svg")

# Figure 2: Correlation distributions
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
ax = axes[0]
ax.hist(results_eigen['rho_activation'], bins=100, color=ACT_COLOR, alpha=0.7, edgecolor='none')
ax.axvline(0.8, color='red', ls='--', lw=1.5, label='rho=0.8')
ax.axvline(-0.8, color='red', ls='--', lw=1.5)
ax.axvline(0.7, color='orange', ls='--', lw=1, label='rho=0.7')
ax.axvline(-0.7, color='orange', ls='--', lw=1)
ax.set_xlabel('Spearman rho'); ax.set_ylabel('Number of genes')
ax.set_title('Activation eigengene', fontsize=12, fontweight='bold')
ax.legend(fontsize=9)

ax = axes[1]
ax.hist(results_eigen['rho_repression'], bins=100, color=REP_COLOR, alpha=0.7, edgecolor='none')
ax.axvline(0.8, color='red', ls='--', lw=1.5); ax.axvline(-0.8, color='red', ls='--', lw=1.5)
ax.axvline(0.7, color='orange', ls='--', lw=1); ax.axvline(-0.7, color='orange', ls='--', lw=1)
ax.set_xlabel('Spearman rho'); ax.set_ylabel('Number of genes')
ax.set_title('Repression eigengene', fontsize=12, fontweight='bold')

ax = axes[2]
act_mask = results_eigen['in_act_regulon_0.8']
rep_mask = results_eigen['in_rep_regulon_0.8']
both_mask = act_mask & rep_mask; act_only = act_mask & ~rep_mask; rep_only = rep_mask & ~act_mask
ax.scatter(results_eigen['rho_activation'], results_eigen['rho_repression'], s=1, alpha=0.2, c=NEUTRAL_COLOR)
ax.scatter(results_eigen.loc[act_only, 'rho_activation'], results_eigen.loc[act_only, 'rho_repression'],
           s=8, alpha=0.6, c=ACT_COLOR, label=f'Act ({act_only.sum()})')
ax.scatter(results_eigen.loc[rep_only, 'rho_activation'], results_eigen.loc[rep_only, 'rho_repression'],
           s=8, alpha=0.6, c=REP_COLOR, label=f'Rep ({rep_only.sum()})')
if both_mask.sum() > 0:
    ax.scatter(results_eigen.loc[both_mask, 'rho_activation'], results_eigen.loc[both_mask, 'rho_repression'],
               s=12, alpha=0.8, c=BOTH_COLOR, label=f'Both ({both_mask.sum()})')
ax.axhline(0, color='gray', ls='-', lw=0.5); ax.axvline(0, color='gray', ls='-', lw=0.5)
ax.set_xlabel('rho (activation)'); ax.set_ylabel('rho (repression)')
ax.set_title(f'Cross-bloc (r={r_cross:.3f})', fontsize=12, fontweight='bold')
ax.legend(fontsize=8, loc='upper left')
plt.tight_layout()
fig.savefig(f'{FIGURES}/eigengene_correlation_distribution.pdf', bbox_inches='tight', dpi=300)
fig.savefig(f'{FIGURES}/eigengene_correlation_distribution.svg', bbox_inches='tight', dpi=300)
plt.close(fig)
print("  Saved: figures/eigengene_correlation_distribution.pdf/svg")

# Figure 3: Functional enrichment heatmap
if len(enrich_df) > 0:
    focus_sets = ['act_pos_0.8', 'act_neg_0.8', 'rep_pos_0.8', 'rep_neg_0.8']
    focus_enrich = enrich_df[enrich_df['gene_set'].isin(focus_sets)].copy()
    sig_keywords = focus_enrich[focus_enrich['pvalue'] < 0.1]['keyword'].unique()
    if len(sig_keywords) < 5:
        sig_keywords = focus_enrich.groupby('keyword')['pvalue'].min().sort_values().head(15).index.values
    if len(sig_keywords) > 0:
        fig, ax = plt.subplots(figsize=(10, max(6, len(sig_keywords) * 0.4)))
        pivot_data = focus_enrich[focus_enrich['keyword'].isin(sig_keywords)].pivot_table(
            index='keyword', columns='gene_set', values='fold_enrichment', aggfunc='first')
        pivot_log2 = np.log2(pivot_data.replace(0, np.nan).replace(np.inf, np.nan))
        pivot_pval = focus_enrich[focus_enrich['keyword'].isin(sig_keywords)].pivot_table(
            index='keyword', columns='gene_set', values='pvalue', aggfunc='first')
        im = ax.imshow(pivot_log2.values, cmap='RdBu_r', aspect='auto', vmin=-2, vmax=2)
        ax.set_xticks(range(len(pivot_log2.columns)))
        ax.set_xticklabels([c.replace('_0.8', '') for c in pivot_log2.columns], rotation=45, ha='right', fontsize=10)
        ax.set_yticks(range(len(pivot_log2.index)))
        ax.set_yticklabels(pivot_log2.index, fontsize=9)
        for i in range(len(pivot_log2.index)):
            for j in range(len(pivot_log2.columns)):
                kw = pivot_log2.index[i]; gsn = pivot_log2.columns[j]
                if gsn in pivot_pval.columns and kw in pivot_pval.index:
                    pv = pivot_pval.loc[kw, gsn]
                    if pd.notna(pv):
                        m = '***' if pv < 0.001 else ('**' if pv < 0.01 else ('*' if pv < 0.05 else ''))
                        ax.text(j, i, m, ha='center', va='center', fontsize=8, fontweight='bold')
        plt.colorbar(im, ax=ax, label='log2(fold enrichment)', shrink=0.8)
        ax.set_title('Functional keyword enrichment in co-expression regulons', fontsize=12, fontweight='bold')
        plt.tight_layout()
        fig.savefig(f'{FIGURES}/functional_enrichment.pdf', bbox_inches='tight', dpi=300)
        fig.savefig(f'{FIGURES}/functional_enrichment.svg', bbox_inches='tight', dpi=300)
        plt.close(fig)
        print("  Saved: figures/functional_enrichment.pdf/svg")

# Figure 4: Regulon size distribution
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
data_bp = [act_sizes, rep_sizes]
ax = axes[0]
bp = ax.boxplot(data_bp, labels=['Activation', 'Repression'], patch_artist=True, widths=0.6)
bp['boxes'][0].set_facecolor(ACT_COLOR); bp['boxes'][0].set_alpha(0.5)
bp['boxes'][1].set_facecolor(REP_COLOR); bp['boxes'][1].set_alpha(0.5)
np.random.seed(789)
for i, d in enumerate(data_bp):
    x = np.random.normal(i+1, 0.05, size=len(d))
    ax.scatter(x, d, s=20, alpha=0.6, c=ACT_COLOR if i==0 else REP_COLOR, zorder=5)
ax.axhline(30, color='gray', ls='--', lw=1, label='Typical regulon (~30)')
ax.set_ylabel('Regulon size (rho > 0.8)')
ax.set_title(f'Per-TF regulon size\nMWU p={pval_mwu:.3f}', fontweight='bold')
ax.legend(fontsize=9)

ax = axes[1]
bins = np.arange(0, max(regulon_df['n_correlated_0.8'].max() + 50, 300), 20)
ax.hist(act_sizes, bins=bins, alpha=0.5, color=ACT_COLOR, label=f"Activation (n={len(act_sizes)})")
ax.hist(rep_sizes, bins=bins, alpha=0.5, color=REP_COLOR, label=f"Repression (n={len(rep_sizes)})")
ax.axvline(30, color='gray', ls='--', lw=1)
ax.set_xlabel('Regulon size'); ax.set_ylabel('Number of TFs')
ax.set_title('Regulon size distribution', fontweight='bold')
ax.legend(fontsize=9)
plt.tight_layout()
fig.savefig(f'{FIGURES}/regulon_size_distribution.pdf', bbox_inches='tight', dpi=300)
fig.savefig(f'{FIGURES}/regulon_size_distribution.svg', bbox_inches='tight', dpi=300)
plt.close(fig)
print("  Saved: figures/regulon_size_distribution.pdf/svg")

# Figure 5: Comprehensive multi-panel summary
fig = plt.figure(figsize=(20, 24))
gs_main = gridspec.GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.3)

# A: Eigengene dynamics
ax = fig.add_subplot(gs_main[0, 0])
x_positions = [0, 0, 0, 1, 1, 1, 2, 2, 2]
np.random.seed(123)
for i in range(9):
    tp = x_positions[i]; jit = np.random.uniform(-0.1, 0.1)
    ax.scatter(tp+jit, act_eigengene[i], c=ACT_COLOR, s=40, zorder=5, alpha=0.7)
    ax.scatter(tp+jit, rep_eigengene[i], c=REP_COLOR, s=40, zorder=5, alpha=0.7, marker='s')
for tp in [0, 1, 2]:
    m = [x == tp for x in x_positions]
    am = np.mean([act_eigengene[i] for i in range(9) if m[i]])
    rm = np.mean([rep_eigengene[i] for i in range(9) if m[i]])
    ax.scatter(tp, am, c=ACT_COLOR, s=100, marker='D', edgecolor='black', zorder=10,
               label='Activation' if tp==0 else '')
    ax.scatter(tp, rm, c=REP_COLOR, s=100, marker='D', edgecolor='black', zorder=10,
               label='Repression' if tp==0 else '')
ax.set_xticks([0,1,2]); ax.set_xticklabels(['T1','T2','T3'])
ax.set_ylabel('Eigengene (median Z)'); ax.set_title('A. Bloc eigengene dynamics', fontweight='bold')
ax.legend(fontsize=8)

# B: Distributions
ax = fig.add_subplot(gs_main[0, 1])
ax.hist(rho_act, bins=80, color=ACT_COLOR, alpha=0.5, label='Activation', density=True)
ax.hist(rho_rep, bins=80, color=REP_COLOR, alpha=0.5, label='Repression', density=True)
ax.axvline(0.8, color='black', ls='--', lw=1); ax.axvline(-0.8, color='black', ls='--', lw=1)
ax.set_xlabel('Spearman rho'); ax.set_ylabel('Density')
ax.set_title('B. Correlation distributions', fontweight='bold'); ax.legend(fontsize=8)

# C: Cross-correlation
ax = fig.add_subplot(gs_main[0, 2])
ax.scatter(rho_act, rho_rep, s=1, alpha=0.2, c=NEUTRAL_COLOR)
ax.scatter(results_eigen.loc[act_only, 'rho_activation'], results_eigen.loc[act_only, 'rho_repression'],
           s=5, alpha=0.5, c=ACT_COLOR, label=f'Act ({act_only.sum()})')
ax.scatter(results_eigen.loc[rep_only, 'rho_activation'], results_eigen.loc[rep_only, 'rho_repression'],
           s=5, alpha=0.5, c=REP_COLOR, label=f'Rep ({rep_only.sum()})')
ax.axhline(0, color='gray', ls='-', lw=0.3); ax.axvline(0, color='gray', ls='-', lw=0.3)
ax.set_xlabel('rho (activation)'); ax.set_ylabel('rho (repression)')
ax.set_title(f'C. Cross-bloc (r={r_cross:.3f})', fontweight='bold')
ax.legend(fontsize=7, loc='upper left')

# D: Venn
ax = fig.add_subplot(gs_main[1, 0])
sa08 = gene_sets[0.8]['act_pos']; sr08 = gene_sets[0.8]['rep_pos']
oa_n = len(sa08-sr08); or_n = len(sr08-sa08); b_n = len(sa08&sr08)
if oa_n + or_n + b_n > 0:
    v = venn2(subsets=(oa_n, or_n, b_n), set_labels=('Activation','Repression'), ax=ax)
    for pid, col in [('10',ACT_COLOR),('01',REP_COLOR),('11',BOTH_COLOR)]:
        if v.get_patch_by_id(pid):
            v.get_patch_by_id(pid).set_color(col); v.get_patch_by_id(pid).set_alpha(0.5)
jac08 = len(sa08&sr08)/len(sa08|sr08) if len(sa08|sr08)>0 else 0
ax.set_title(f'D. Regulon overlap (rho>0.8)\nJaccard={jac08:.3f}', fontweight='bold')

# E: Regulon size
ax = fig.add_subplot(gs_main[1, 1])
bp = ax.boxplot(data_bp, labels=['Activation','Repression'], patch_artist=True, widths=0.5)
bp['boxes'][0].set_facecolor(ACT_COLOR); bp['boxes'][0].set_alpha(0.5)
bp['boxes'][1].set_facecolor(REP_COLOR); bp['boxes'][1].set_alpha(0.5)
np.random.seed(456)
for i, d in enumerate(data_bp):
    x = np.random.normal(i+1, 0.04, size=len(d))
    ax.scatter(x, d, s=15, alpha=0.5, c=ACT_COLOR if i==0 else REP_COLOR, zorder=5)
ax.axhline(30, color='gray', ls='--', lw=1)
ax.set_ylabel('Regulon size'); ax.set_title(f'E. Per-TF regulon\nMWU p={pval_mwu:.3f}', fontweight='bold')

# F: Permutation
ax = fig.add_subplot(gs_main[1, 2])
ax.hist(perm_act_sizes, bins=30, color=ACT_COLOR, alpha=0.5, label='Perm act', density=True)
ax.hist(perm_rep_sizes, bins=30, color=REP_COLOR, alpha=0.5, label='Perm rep', density=True)
ax.axvline(real_act_size_08, color=ACT_COLOR, ls='-', lw=2, label=f'Real act ({real_act_size_08})')
ax.axvline(real_rep_size_08, color=REP_COLOR, ls='-', lw=2, label=f'Real rep ({real_rep_size_08})')
ax.set_xlabel('Regulon size (rho>0.8)'); ax.set_ylabel('Density')
ax.set_title(f'F. Permutation test\np_act={perm_p_act:.3f}, p_rep={perm_p_rep:.3f}', fontweight='bold')
ax.legend(fontsize=7)

# G: DEG overlap
ax = fig.add_subplot(gs_main[2, 0])
if len(deg_overlap_df) > 0:
    pd2 = deg_overlap_df[
        (deg_overlap_df['gene_set'].isin(['act_pos_0.8','rep_pos_0.8'])) &
        (deg_overlap_df['deg_direction'].isin(['up','down']))
    ].copy()
    if len(pd2) > 0:
        cats = []; vals = []; cols = []
        for _, row in pd2.iterrows():
            cats.append(f"{row['gene_set'][:3]}_{row['transition']}_{row['deg_direction']}")
            fe = row['fold_enrichment']
            vals.append(np.log2(fe) if fe > 0 and not np.isinf(fe) else 0)
            cols.append(ACT_COLOR if 'act' in row['gene_set'] else REP_COLOR)
        ax.bar(range(len(cats)), vals, color=cols, alpha=0.7)
        ax.set_xticks(range(len(cats))); ax.set_xticklabels(cats, rotation=45, ha='right', fontsize=6)
        ax.axhline(0, color='gray', ls='-', lw=0.5); ax.set_ylabel('log2(fold enrichment)')
        for i, (_, row) in enumerate(pd2.iterrows()):
            m = '***' if row['padj']<0.001 else ('**' if row['padj']<0.01 else ('*' if row['padj']<0.05 else ''))
            ax.text(i, vals[i]+0.05, m, ha='center', fontsize=7)
ax.set_title('G. DEG overlap', fontweight='bold')

# H: Functional enrichment
ax = fig.add_subplot(gs_main[2, 1])
if len(enrich_df) > 0:
    focus = enrich_df[enrich_df['gene_set'].isin(['act_pos_0.8','rep_pos_0.8'])].copy()
    if len(focus) > 0:
        top_kw = focus.groupby('keyword')['pvalue'].min().sort_values().head(10).index
        for gsn, color in [('act_pos_0.8',ACT_COLOR),('rep_pos_0.8',REP_COLOR)]:
            sub = focus[(focus['gene_set']==gsn) & (focus['keyword'].isin(top_kw))].set_index('keyword')
            if len(sub) > 0:
                l2 = np.log2(sub.reindex(top_kw)['fold_enrichment'].replace(0, np.nan))
                offset = 0.15 if 'act' in gsn else -0.15
                ax.barh([y+offset for y in range(len(top_kw))], l2.values, height=0.3,
                        color=color, alpha=0.6, label='Act' if 'act' in gsn else 'Rep')
        ax.set_yticks(range(len(top_kw))); ax.set_yticklabels(top_kw, fontsize=8)
        ax.axvline(0, color='gray', ls='-', lw=0.5); ax.set_xlabel('log2(fold enrichment)')
        ax.legend(fontsize=8)
ax.set_title('H. Top enrichments', fontweight='bold')

# I: Summary text
ax = fig.add_subplot(gs_main[2, 2]); ax.axis('off')
txt = (
    f"SUMMARY STATISTICS\n{'='*35}\n\n"
    f"Genes tested: {n_genes:,}\n"
    f"Activation bloc: {len(activation_tags)} TFs\n"
    f"Repression bloc: {len(repression_tags)} TFs\n\n"
    f"Regulon (rho>0.8):\n"
    f"  Act-corr: {len(act_regulon_pos)}\n"
    f"  Rep-corr: {len(rep_regulon_pos)}\n"
    f"  Overlap: {len(act_regulon_pos&rep_regulon_pos)}\n"
    f"  Jaccard: {jac08:.3f}\n\n"
    f"Regulon (rho>0.7):\n"
    f"  Act-corr: {len(gene_sets[0.7]['act_pos'])}\n"
    f"  Rep-corr: {len(gene_sets[0.7]['rep_pos'])}\n\n"
    f"Cross-eigengene r: {r_cross:.3f}\n"
    f"Perm p(act): {perm_p_act:.3f}\n"
    f"Perm p(rep): {perm_p_rep:.3f}\n\n"
    f"Per-TF (rho>0.8):\n"
    f"  Act median: {np.median(act_sizes):.0f}\n"
    f"  Rep median: {np.median(rep_sizes):.0f}\n"
    f"  MWU p: {pval_mwu:.3f}"
)
ax.text(0.05, 0.95, txt, transform=ax.transAxes, fontsize=10, verticalalignment='top',
        fontfamily='monospace', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
ax.set_title('I. Key statistics', fontweight='bold')

fig.suptitle('H37: Co-expression Regulon Prediction for 62 Exposed TFs\nS. coelicolor M145',
             fontsize=16, fontweight='bold', y=1.01)
plt.tight_layout()
fig.savefig(f'{FIGURES}/H37_comprehensive_summary.pdf', bbox_inches='tight', dpi=300)
fig.savefig(f'{FIGURES}/H37_comprehensive_summary.svg', bbox_inches='tight', dpi=300)
plt.close(fig)
print("  Saved: figures/H37_comprehensive_summary.pdf/svg")

# ============================================================
# FINAL SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)

print(f"\n  Key findings:")
print(f"  1. Activation regulon (rho>0.8): {len(act_regulon_pos)} genes ({len(act_regulon_pos)/n_genes*100:.1f}%)")
print(f"  2. Repression regulon (rho>0.8): {len(rep_regulon_pos)} genes ({len(rep_regulon_pos)/n_genes*100:.1f}%)")
print(f"  3. Overlap: {len(act_regulon_pos&rep_regulon_pos)} (Jaccard={jac08:.3f})")
print(f"  4. Anti-correlated negatively: {len(act_regulon_neg)} act-neg, {len(rep_regulon_neg)} rep-neg")
print(f"  5. Cross-eigengene: r={r_cross:.3f}")
print(f"  6. Permutation: act p={perm_p_act:.4f} (Z={act_z:.1f}), rep p={perm_p_rep:.4f} (Z={rep_z:.1f})")
print(f"  7. Per-TF: act median={np.median(act_sizes):.0f}, rep median={np.median(rep_sizes):.0f}")

print(f"\n  Output: {OUT_DIR}")
print("\n  Done!")
