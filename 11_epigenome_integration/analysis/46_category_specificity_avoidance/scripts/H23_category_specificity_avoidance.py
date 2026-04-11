#!/usr/bin/env python3
"""
H23: Gene Category Specificity of Methylation Avoidance
========================================================
Tests whether methylation avoidance extends beyond regulatory genes to other
functionally important gene categories.

Key question: Is avoidance universal (all important genes) or regulatory-specific?
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import fisher_exact, mannwhitneyu
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import TwoSlopeNorm
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Configuration
# ============================================================
BASE = '/Users/okaban/bioinfo/rna-seq'
ANALYSIS_DIR = f'{BASE}/11_epigenome_integration/analysis/46_category_specificity_avoidance'
FIG_DIR = f'{ANALYSIS_DIR}/figures'
TBL_DIR = f'{ANALYSIS_DIR}/tables'

CHROM_SIZE = 8_667_507
ARM_LEFT_END = 1_500_000
ARM_RIGHT_START = 7_167_508
PROXIMITY_KB = 2  # 2 kb proximity window

# ============================================================
# 1. Load data
# ============================================================
print("=" * 70)
print("H23: Gene Category Specificity of Methylation Avoidance")
print("=" * 70)

# Gene annotations
gene_annot = pd.read_csv(
    f'{BASE}/05_annotation/analysis/05_annotation_260128_v1/tables/gene_annotation_basic.tsv',
    sep='\t'
)
print(f"Gene annotations: {len(gene_annot)} genes")

# DESeq2 results (for expression quintiles)
deseq = pd.read_csv(
    f'{BASE}/04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_2_vs_1.tsv',
    sep='\t'
)
print(f"DESeq2 results: {len(deseq)} genes")

# Regulatory genes (1,055)
reg_genes = pd.read_csv(
    f'{BASE}/11_epigenome_integration/analysis/29_genomewide_TF_screen/tables/all_regulatory_genes.tsv',
    sep='\t'
)
reg_locus_tags = set(reg_genes['locus_tag'].values)
print(f"Regulatory genes: {len(reg_locus_tags)}")

# BGC master
bgc_master = pd.read_csv(
    f'{BASE}/05_annotation/analysis/05_annotation_260128_v1/tables/gene_master_with_BGC.tsv',
    sep='\t'
)
bgc_genes = set(bgc_master.loc[bgc_master['bgc_name'].notna(), 'gene_id'].values)
print(f"BGC-assigned genes: {len(bgc_genes)}")

# GCCGGC sites - T1
gccggc_all = pd.read_csv(
    f'{BASE}/11_epigenome_integration/analysis/37_defense_island_GCCGGC/tables/GCCGGC_sites_by_timepoint.tsv',
    sep='\t'
)
gccggc_t1 = gccggc_all[gccggc_all['timepoint'] == 'T1'].copy()
gccggc_t1_positions = gccggc_t1['position'].unique()
print(f"GCCGGC 4mC T1 sites: {len(gccggc_t1_positions)}")

# AAGCCCG sites - T1 (from gene mapping, get unique positions)
aagcccg_map = pd.read_csv(
    f'{BASE}/11_epigenome_integration/analysis/36_AAGCCCG_distribution/tables/AAGCCCG_site_gene_mapping.tsv',
    sep='\t'
)
aagcccg_t1_map = aagcccg_map[aagcccg_map['timepoint'] == 'T1']
aagcccg_t1_positions = aagcccg_t1_map['position'].unique()
print(f"AAGCCCG 6mA T1 sites (from mapping): {len(aagcccg_t1_positions)}")

# All 4mC sites - T1 (from census)
all_4mc = pd.read_csv(
    f'{BASE}/11_epigenome_integration/analysis/23_expanded_motif_search/4mC_final_census.csv'
)
all_4mc_t1 = all_4mc[all_4mc['timepoint'] == 'T1']
all_4mc_t1_positions = all_4mc_t1['position'].unique()
print(f"All 4mC T1 sites: {len(all_4mc_t1_positions)}")

# ============================================================
# 2. Classify genes into functional categories
# ============================================================
print("\n" + "=" * 70)
print("Step 1: Gene classification")
print("=" * 70)

def classify_gene(row, reg_set, bgc_set):
    """Classify a gene into functional category based on product annotation."""
    gene_id = row['gene_id']
    product = str(row.get('product', '')).lower()

    # Priority order matters - a gene gets assigned to FIRST matching category
    # Regulatory (from curated list)
    if gene_id in reg_set:
        return 'Regulatory'

    # Translation
    if any(kw in product for kw in ['ribosom', 'trna', 'translation factor',
                                     'elongation factor', 'trna synthetase',
                                     'trna ligase', 'aminoacyl-trna']):
        return 'Translation'

    # DNA replication/repair
    if any(kw in product for kw in ['dna polymer', 'helicase', 'gyrase',
                                     'dna ligase', 'recombin', 'topoisomerase',
                                     'dna repair', 'dnaa', 'ssb ', 'primase']):
        return 'DNA_replication_repair'

    # Cell division
    if any(kw in product for kw in ['ftsz', 'division', 'septum', 'cell wall',
                                     'peptidoglycan', 'mur ', 'ftsk', 'ftsw',
                                     'ftsq', 'ftsi', 'penicillin-binding',
                                     'murein']):
        return 'Cell_division'

    # Energy/Respiration
    if any(kw in product for kw in ['cytochrome', 'atp synthase', 'nadh',
                                     'nadph', 'succinate dehydrogenase',
                                     'ubiquinol', 'ferredoxin', 'electron transfer']):
        return 'Energy_respiration'

    # Transport
    if any(kw in product for kw in ['transporter', 'permease', 'mfs ',
                                     'abc transport', 'efflux', 'porin',
                                     'antiporter', 'symporter']):
        return 'Transport'

    # Secondary metabolism (from BGC annotation)
    if gene_id in bgc_set:
        return 'Secondary_metabolism'

    # Hypothetical
    if 'hypothetical' in product:
        return 'Hypothetical'

    return 'Other'


gene_annot['category'] = gene_annot.apply(
    lambda row: classify_gene(row, reg_locus_tags, bgc_genes), axis=1
)

# Assign genomic region
gene_annot['midpoint'] = (gene_annot['start'] + gene_annot['end']) / 2
gene_annot['region'] = gene_annot['midpoint'].apply(
    lambda x: 'arm' if x <= ARM_LEFT_END or x >= ARM_RIGHT_START else 'core'
)

cat_counts = gene_annot['category'].value_counts()
print("\nGene category counts:")
for cat, cnt in cat_counts.items():
    print(f"  {cat}: {cnt}")
print(f"  TOTAL: {len(gene_annot)}")

# ============================================================
# 3. Calculate methylation proximity for each category
# ============================================================
print("\n" + "=" * 70)
print("Step 2: Methylation proximity analysis by category")
print("=" * 70)

def gene_has_proximal_site(gene_start, gene_end, site_positions, window_bp=2000):
    """Check if a gene has any methylation site within window_bp of its start or end."""
    # Check if any site is within window_bp of gene boundaries
    for pos in site_positions:
        if (gene_start - window_bp) <= pos <= (gene_end + window_bp):
            return True
    return False


def calculate_proximity_vectorized(gene_df, site_positions, window_bp=2000):
    """Vectorized check: for each gene, is there a methylation site within window_bp?"""
    sites = np.array(sorted(site_positions))
    proximal = np.zeros(len(gene_df), dtype=bool)

    for i, (_, row) in enumerate(gene_df.iterrows()):
        g_start = row['start']
        g_end = row['end']
        # Binary search for efficiency
        left = np.searchsorted(sites, g_start - window_bp)
        right = np.searchsorted(sites, g_end + window_bp, side='right')
        if left < right:
            proximal[i] = True

    return proximal


# Calculate proximity for all three methylation types
print("\nCalculating proximity (within 2 kb)...")
window = PROXIMITY_KB * 1000

gene_annot['prox_GCCGGC'] = calculate_proximity_vectorized(gene_annot, gccggc_t1_positions, window)
gene_annot['prox_AAGCCCG'] = calculate_proximity_vectorized(gene_annot, aagcccg_t1_positions, window)
gene_annot['prox_All4mC'] = calculate_proximity_vectorized(gene_annot, all_4mc_t1_positions, window)

for methyl_type in ['GCCGGC', 'AAGCCCG', 'All4mC']:
    col = f'prox_{methyl_type}'
    n_prox = gene_annot[col].sum()
    print(f"  {methyl_type}: {n_prox}/{len(gene_annot)} genes proximal ({n_prox/len(gene_annot)*100:.1f}%)")


def compute_category_enrichment(gene_df, prox_col, methyl_label):
    """Compute Fisher enrichment for each gene category vs genome-wide rate."""
    total_genes = len(gene_df)
    total_proximal = gene_df[prox_col].sum()
    genome_rate = total_proximal / total_genes

    results = []
    categories = sorted(gene_df['category'].unique())

    for cat in categories:
        cat_mask = gene_df['category'] == cat
        n_genes = cat_mask.sum()
        n_prox = gene_df.loc[cat_mask, prox_col].sum()

        # Build 2x2 table: [[cat_prox, cat_not_prox], [other_prox, other_not_prox]]
        a = n_prox
        b = n_genes - n_prox
        c = total_proximal - n_prox
        d = (total_genes - n_genes) - c

        table = np.array([[a, b], [c, d]])
        OR, p_val = fisher_exact(table)

        # 95% CI for OR (Woolf method)
        if a > 0 and b > 0 and c > 0 and d > 0:
            log_OR = np.log(OR)
            se_log_OR = np.sqrt(1/a + 1/b + 1/c + 1/d)
            CI_low = np.exp(log_OR - 1.96 * se_log_OR)
            CI_high = np.exp(log_OR + 1.96 * se_log_OR)
        else:
            CI_low = 0
            CI_high = np.inf

        frac_prox = n_prox / n_genes if n_genes > 0 else 0
        fold = frac_prox / genome_rate if genome_rate > 0 else np.nan
        expected = genome_rate * n_genes

        results.append({
            'methyl_type': methyl_label,
            'category': cat,
            'n_genes': n_genes,
            'n_proximal': int(n_prox),
            'expected_proximal': round(expected, 1),
            'frac_proximal': round(frac_prox, 4),
            'genome_rate': round(genome_rate, 4),
            'fold_enrichment': round(fold, 3),
            'odds_ratio': round(OR, 3),
            'OR_CI_low': round(CI_low, 3),
            'OR_CI_high': round(CI_high, 3),
            'p_value': p_val,
            'direction': 'depleted' if OR < 1 else 'enriched'
        })

    return pd.DataFrame(results)


# Run for all three methylation types
all_enrichment = []
for methyl_type, prox_col in [('GCCGGC_4mC_T1', 'prox_GCCGGC'),
                                ('AAGCCCG_6mA_T1', 'prox_AAGCCCG'),
                                ('All_4mC_T1', 'prox_All4mC')]:
    df = compute_category_enrichment(gene_annot, prox_col, methyl_type)
    all_enrichment.append(df)

enrichment_df = pd.concat(all_enrichment, ignore_index=True)

# Bonferroni correction (n categories × 3 methylation types)
n_tests = len(enrichment_df)
enrichment_df['p_bonferroni'] = np.minimum(enrichment_df['p_value'] * n_tests, 1.0)

print("\n--- Category Enrichment Summary ---")
print(f"{'Category':<25} {'Methyl Type':<18} {'n':>5} {'Obs':>5} {'Exp':>6} {'Fold':>6} {'OR':>6} {'p':>10} {'p_bonf':>10}")
print("-" * 105)
for _, row in enrichment_df.sort_values(['category', 'methyl_type']).iterrows():
    p_str = f"{row['p_value']:.2e}" if row['p_value'] < 0.05 else f"{row['p_value']:.3f}"
    pb_str = f"{row['p_bonferroni']:.2e}" if row['p_bonferroni'] < 0.05 else f"{row['p_bonferroni']:.3f}"
    print(f"{row['category']:<25} {row['methyl_type']:<18} {row['n_genes']:>5} {row['n_proximal']:>5} {row['expected_proximal']:>6} {row['fold_enrichment']:>6.2f} {row['odds_ratio']:>6.2f} {p_str:>10} {pb_str:>10}")

# Save
enrichment_df.to_csv(f'{TBL_DIR}/category_enrichment.tsv', sep='\t', index=False)
print(f"\nSaved: {TBL_DIR}/category_enrichment.tsv")

# ============================================================
# 4. Expression quintile analysis
# ============================================================
print("\n" + "=" * 70)
print("Step 3: Expression quintile analysis")
print("=" * 70)

# Merge expression data
gene_with_expr = gene_annot.merge(deseq[['gene_id', 'baseMean']], on='gene_id', how='left')
# Only genes with expression data
gene_expr = gene_with_expr.dropna(subset=['baseMean']).copy()
gene_expr = gene_expr[gene_expr['baseMean'] > 0].copy()

# Quintiles
gene_expr['expr_quintile'] = pd.qcut(gene_expr['baseMean'], 5, labels=['Q1_lowest', 'Q2', 'Q3', 'Q4', 'Q5_highest'])

print(f"Genes with expression data: {len(gene_expr)}")
print(f"\nQuintile baseMean ranges:")
for q in ['Q1_lowest', 'Q2', 'Q3', 'Q4', 'Q5_highest']:
    qdata = gene_expr[gene_expr['expr_quintile'] == q]['baseMean']
    print(f"  {q}: {qdata.min():.1f} - {qdata.max():.1f} (n={len(qdata)})")

# Quintile enrichment
quintile_results = []
for methyl_type, prox_col in [('GCCGGC_4mC_T1', 'prox_GCCGGC'),
                                ('AAGCCCG_6mA_T1', 'prox_AAGCCCG'),
                                ('All_4mC_T1', 'prox_All4mC')]:
    total = len(gene_expr)
    total_prox = gene_expr[prox_col].sum()
    genome_rate = total_prox / total

    for q in ['Q1_lowest', 'Q2', 'Q3', 'Q4', 'Q5_highest']:
        q_mask = gene_expr['expr_quintile'] == q
        n_genes = q_mask.sum()
        n_prox = gene_expr.loc[q_mask, prox_col].sum()

        a = int(n_prox)
        b = int(n_genes - n_prox)
        c = int(total_prox - n_prox)
        d = int((total - n_genes) - c)

        table = np.array([[a, b], [c, d]])
        OR, p_val = fisher_exact(table)

        if a > 0 and b > 0 and c > 0 and d > 0:
            log_OR = np.log(OR)
            se_log_OR = np.sqrt(1/a + 1/b + 1/c + 1/d)
            CI_low = np.exp(log_OR - 1.96 * se_log_OR)
            CI_high = np.exp(log_OR + 1.96 * se_log_OR)
        else:
            CI_low = 0
            CI_high = np.inf

        frac_prox = n_prox / n_genes if n_genes > 0 else 0
        fold = frac_prox / genome_rate if genome_rate > 0 else np.nan

        quintile_results.append({
            'methyl_type': methyl_type,
            'quintile': q,
            'n_genes': n_genes,
            'n_proximal': int(n_prox),
            'expected': round(genome_rate * n_genes, 1),
            'frac_proximal': round(frac_prox, 4),
            'genome_rate': round(genome_rate, 4),
            'fold_enrichment': round(fold, 3),
            'odds_ratio': round(OR, 3),
            'OR_CI_low': round(CI_low, 3),
            'OR_CI_high': round(CI_high, 3),
            'p_value': p_val
        })

quintile_df = pd.DataFrame(quintile_results)
n_qtests = len(quintile_df)
quintile_df['p_bonferroni'] = np.minimum(quintile_df['p_value'] * n_qtests, 1.0)

print("\n--- Expression Quintile Enrichment ---")
print(f"{'Methyl Type':<18} {'Quintile':<12} {'n':>5} {'Obs':>5} {'Exp':>6} {'Fold':>6} {'OR':>6} {'p':>10} {'p_bonf':>10}")
print("-" * 95)
for _, row in quintile_df.iterrows():
    p_str = f"{row['p_value']:.2e}" if row['p_value'] < 0.05 else f"{row['p_value']:.3f}"
    pb_str = f"{row['p_bonferroni']:.2e}" if row['p_bonferroni'] < 0.05 else f"{row['p_bonferroni']:.3f}"
    print(f"{row['methyl_type']:<18} {row['quintile']:<12} {row['n_genes']:>5} {row['n_proximal']:>5} {row['expected']:>6} {row['fold_enrichment']:>6.2f} {row['odds_ratio']:>6.2f} {p_str:>10} {pb_str:>10}")

# Jonckheere-Terpstra trend test (manual implementation)
def jonckheere_terpstra_test(groups, values):
    """
    Jonckheere-Terpstra test for ordered alternatives.
    groups: ordered group labels (list)
    values: dict of group_label -> array of values

    Returns JT statistic and approximate z-score p-value.
    """
    k = len(groups)
    U = 0
    n_total = sum(len(values[g]) for g in groups)
    ns = [len(values[g]) for g in groups]

    for i in range(k - 1):
        for j in range(i + 1, k):
            xi = values[groups[i]]
            xj = values[groups[j]]
            for a in xi:
                for b in xj:
                    if b > a:
                        U += 1
                    elif b == a:
                        U += 0.5

    # Expected value under H0
    N = n_total
    E_U = (N**2 - sum(n**2 for n in ns)) / 4

    # Variance (no ties simplification)
    Var_U = (N**2 * (2*N + 3) - sum(n**2 * (2*n + 3) for n in ns)) / 72

    if Var_U > 0:
        z = (U - E_U) / np.sqrt(Var_U)
        p = 2 * (1 - stats.norm.cdf(abs(z)))
    else:
        z = 0
        p = 1.0

    return U, z, p


print("\n--- Jonckheere-Terpstra Trend Tests ---")
quintile_order = ['Q1_lowest', 'Q2', 'Q3', 'Q4', 'Q5_highest']
jt_results = []

for methyl_type, prox_col in [('GCCGGC_4mC_T1', 'prox_GCCGGC'),
                                ('AAGCCCG_6mA_T1', 'prox_AAGCCCG'),
                                ('All_4mC_T1', 'prox_All4mC')]:
    # Use proximity as 0/1 values for each quintile
    groups_data = {}
    for q in quintile_order:
        q_mask = gene_expr['expr_quintile'] == q
        groups_data[q] = gene_expr.loc[q_mask, prox_col].astype(int).values

    U_stat, z_stat, p_val = jonckheere_terpstra_test(quintile_order, groups_data)

    # Direction: if z < 0, higher expression -> lower proximity
    direction = "decreasing" if z_stat < 0 else "increasing"

    jt_results.append({
        'methyl_type': methyl_type,
        'JT_statistic': round(U_stat, 1),
        'z_score': round(z_stat, 3),
        'p_value': p_val,
        'direction': direction
    })

    print(f"  {methyl_type}: JT U={U_stat:.0f}, z={z_stat:.3f}, p={p_val:.2e}, trend={direction}")

jt_df = pd.DataFrame(jt_results)
jt_df['p_bonferroni'] = np.minimum(jt_df['p_value'] * 3, 1.0)

# Combine and save
quintile_df_out = quintile_df.copy()
quintile_df.to_csv(f'{TBL_DIR}/expression_quintile_analysis.tsv', sep='\t', index=False)
jt_df.to_csv(f'{TBL_DIR}/jt_trend_tests.tsv', sep='\t', index=False)
print(f"\nSaved: {TBL_DIR}/expression_quintile_analysis.tsv")
print(f"Saved: {TBL_DIR}/jt_trend_tests.tsv")

# ============================================================
# 5. Ribosomal protein deep dive
# ============================================================
print("\n" + "=" * 70)
print("Step 4: Ribosomal protein deep dive")
print("=" * 70)

# Extract ribosomal proteins
ribo_mask = gene_annot['product'].str.lower().str.contains('ribosomal protein', na=False)
ribo_genes = gene_annot[ribo_mask].copy()
print(f"Ribosomal protein genes: {len(ribo_genes)}")

# Subgroups
ribo_50S = ribo_genes[ribo_genes['product'].str.contains('50S', na=False)]
ribo_30S = ribo_genes[ribo_genes['product'].str.contains('30S', na=False)]
print(f"  50S subunit: {len(ribo_50S)}")
print(f"  30S subunit: {len(ribo_30S)}")

# Avoidance metrics for ribosomal proteins
ribo_results = []
total_genes = len(gene_annot)

for methyl_type, prox_col in [('GCCGGC_4mC_T1', 'prox_GCCGGC'),
                                ('AAGCCCG_6mA_T1', 'prox_AAGCCCG'),
                                ('All_4mC_T1', 'prox_All4mC')]:
    total_prox = gene_annot[prox_col].sum()
    genome_rate = total_prox / total_genes

    for group_name, group_df in [('All_ribosomal', ribo_genes),
                                  ('50S_subunit', ribo_50S),
                                  ('30S_subunit', ribo_30S),
                                  ('Regulatory', gene_annot[gene_annot['category'] == 'Regulatory'])]:
        n = len(group_df)
        n_prox = group_df[prox_col].sum()

        a = int(n_prox)
        b = int(n - n_prox)
        c = int(total_prox - n_prox)
        d = int((total_genes - n) - c)

        table = np.array([[a, b], [c, d]])
        OR, p_val = fisher_exact(table)

        if a > 0 and b > 0 and c > 0 and d > 0:
            log_OR = np.log(OR)
            se_log_OR = np.sqrt(1/a + 1/b + 1/c + 1/d)
            CI_low = np.exp(log_OR - 1.96 * se_log_OR)
            CI_high = np.exp(log_OR + 1.96 * se_log_OR)
        else:
            CI_low = 0
            CI_high = np.inf

        frac = n_prox / n if n > 0 else 0
        fold = frac / genome_rate if genome_rate > 0 else np.nan

        ribo_results.append({
            'methyl_type': methyl_type,
            'group': group_name,
            'n_genes': n,
            'n_proximal': int(n_prox),
            'frac_proximal': round(frac, 4),
            'genome_rate': round(genome_rate, 4),
            'fold_enrichment': round(fold, 3),
            'odds_ratio': round(OR, 3),
            'OR_CI_low': round(CI_low, 3),
            'OR_CI_high': round(CI_high, 3),
            'p_value': p_val
        })

ribo_df = pd.DataFrame(ribo_results)
n_ribo_tests = len(ribo_df)
ribo_df['p_bonferroni'] = np.minimum(ribo_df['p_value'] * n_ribo_tests, 1.0)

print("\n--- Ribosomal Protein vs Regulatory Comparison ---")
print(f"{'Methyl Type':<18} {'Group':<20} {'n':>5} {'Obs':>4} {'Frac':>6} {'Fold':>6} {'OR':>6} {'p':>10} {'p_bonf':>10}")
print("-" * 100)
for _, row in ribo_df.iterrows():
    p_str = f"{row['p_value']:.2e}" if row['p_value'] < 0.05 else f"{row['p_value']:.3f}"
    pb_str = f"{row['p_bonferroni']:.2e}" if row['p_bonferroni'] < 0.05 else f"{row['p_bonferroni']:.3f}"
    print(f"{row['methyl_type']:<18} {row['group']:<20} {row['n_genes']:>5} {row['n_proximal']:>4} {row['frac_proximal']:>6.3f} {row['fold_enrichment']:>6.2f} {row['odds_ratio']:>6.2f} {p_str:>10} {pb_str:>10}")

ribo_df.to_csv(f'{TBL_DIR}/ribosomal_protein_analysis.tsv', sep='\t', index=False)
print(f"\nSaved: {TBL_DIR}/ribosomal_protein_analysis.tsv")

# ============================================================
# 6. Core/arm stratification
# ============================================================
print("\n" + "=" * 70)
print("Step 5: Core/arm stratified analysis")
print("=" * 70)

stratified_results = []

for region in ['core', 'arm']:
    region_genes = gene_annot[gene_annot['region'] == region].copy()
    total_r = len(region_genes)

    print(f"\n--- {region.upper()} region ({total_r} genes) ---")

    for methyl_type, prox_col in [('GCCGGC_4mC_T1', 'prox_GCCGGC'),
                                    ('AAGCCCG_6mA_T1', 'prox_AAGCCCG'),
                                    ('All_4mC_T1', 'prox_All4mC')]:
        total_prox_r = region_genes[prox_col].sum()
        genome_rate_r = total_prox_r / total_r if total_r > 0 else 0

        for cat in sorted(region_genes['category'].unique()):
            cat_mask = region_genes['category'] == cat
            n = cat_mask.sum()
            n_prox = region_genes.loc[cat_mask, prox_col].sum()

            a = int(n_prox)
            b = int(n - n_prox)
            c = int(total_prox_r - n_prox)
            d = int((total_r - n) - c)

            if d < 0:
                d = 0

            table = np.array([[a, b], [c, d]])
            OR, p_val = fisher_exact(table)

            if a > 0 and b > 0 and c > 0 and d > 0:
                log_OR = np.log(OR)
                se_log_OR = np.sqrt(1/a + 1/b + 1/c + 1/d)
                CI_low = np.exp(log_OR - 1.96 * se_log_OR)
                CI_high = np.exp(log_OR + 1.96 * se_log_OR)
            else:
                CI_low = 0
                CI_high = np.inf

            frac = n_prox / n if n > 0 else 0
            fold = frac / genome_rate_r if genome_rate_r > 0 else np.nan

            stratified_results.append({
                'region': region,
                'methyl_type': methyl_type,
                'category': cat,
                'n_genes': n,
                'n_proximal': int(n_prox),
                'frac_proximal': round(frac, 4),
                'region_rate': round(genome_rate_r, 4),
                'fold_enrichment': round(fold, 3),
                'odds_ratio': round(OR, 3),
                'OR_CI_low': round(CI_low, 3),
                'OR_CI_high': round(CI_high, 3),
                'p_value': p_val
            })

stratified_df = pd.DataFrame(stratified_results)
n_strat_tests = len(stratified_df)
stratified_df['p_bonferroni'] = np.minimum(stratified_df['p_value'] * n_strat_tests, 1.0)

# Print key categories
for region in ['core', 'arm']:
    print(f"\n--- {region.upper()} region, key categories ---")
    subset = stratified_df[(stratified_df['region'] == region)]
    for cat in ['Regulatory', 'Translation', 'DNA_replication_repair', 'Energy_respiration', 'Hypothetical']:
        cat_rows = subset[subset['category'] == cat]
        if len(cat_rows) > 0:
            for _, row in cat_rows.iterrows():
                p_str = f"{row['p_value']:.2e}" if row['p_value'] < 0.05 else f"{row['p_value']:.3f}"
                print(f"  {cat:<25} {row['methyl_type']:<18} n={row['n_genes']:>5} fold={row['fold_enrichment']:>5.2f} OR={row['odds_ratio']:>5.2f} p={p_str}")

stratified_df.to_csv(f'{TBL_DIR}/stratified_by_region.tsv', sep='\t', index=False)
print(f"\nSaved: {TBL_DIR}/stratified_by_region.tsv")

# ============================================================
# 7. Visualization
# ============================================================
print("\n" + "=" * 70)
print("Step 6: Generating figures")
print("=" * 70)

# Color scheme
category_colors = {
    'Regulatory': '#E63946',
    'Translation': '#457B9D',
    'DNA_replication_repair': '#2A9D8F',
    'Cell_division': '#E9C46A',
    'Energy_respiration': '#F4A261',
    'Transport': '#264653',
    'Secondary_metabolism': '#9B2226',
    'Hypothetical': '#ADB5BD',
    'Other': '#6C757D'
}

category_labels = {
    'Regulatory': 'Regulatory (1,055)',
    'Translation': 'Translation',
    'DNA_replication_repair': 'DNA replication/repair',
    'Cell_division': 'Cell division',
    'Energy_respiration': 'Energy/respiration',
    'Transport': 'Transport',
    'Secondary_metabolism': 'Secondary metabolism',
    'Hypothetical': 'Hypothetical',
    'Other': 'Other'
}

# ---- Figure 1: Heatmap ----
fig, ax = plt.subplots(figsize=(8, 6))

# Pivot for heatmap
heatmap_data = enrichment_df.pivot(index='category', columns='methyl_type', values='fold_enrichment')
pval_data = enrichment_df.pivot(index='category', columns='methyl_type', values='p_bonferroni')

# Order categories by mean fold enrichment
cat_order = heatmap_data.mean(axis=1).sort_values().index.tolist()
methyl_order = ['GCCGGC_4mC_T1', 'AAGCCCG_6mA_T1', 'All_4mC_T1']
heatmap_data = heatmap_data.loc[cat_order, methyl_order]
pval_data = pval_data.loc[cat_order, methyl_order]

# Custom colormap: blue (depletion) - white (1.0) - red (enrichment)
vmin = min(0.3, heatmap_data.min().min())
vmax = max(1.7, heatmap_data.max().max())
norm = TwoSlopeNorm(vmin=vmin, vcenter=1.0, vmax=vmax)

im = ax.imshow(heatmap_data.values, cmap='RdBu_r', norm=norm, aspect='auto')
cbar = plt.colorbar(im, ax=ax, label='Fold enrichment (vs genome-wide)')

# Add text annotations
for i in range(len(cat_order)):
    for j in range(len(methyl_order)):
        val = heatmap_data.values[i, j]
        pval = pval_data.values[i, j]
        sig = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else ''
        color = 'white' if abs(val - 1.0) > 0.3 else 'black'
        ax.text(j, i, f'{val:.2f}{sig}', ha='center', va='center', fontsize=9, color=color, fontweight='bold')

# Labels
y_labels = [category_labels.get(c, c) for c in cat_order]
ax.set_yticks(range(len(cat_order)))
ax.set_yticklabels(y_labels, fontsize=10)
ax.set_xticks(range(len(methyl_order)))
ax.set_xticklabels(['GCCGGC 4mC', 'AAGCCCG 6mA', 'All 4mC'], fontsize=10, rotation=30, ha='right')
ax.set_title('H23: Gene Category × Methylation Proximity\nFold Enrichment (2 kb window, T1)', fontsize=12, fontweight='bold')

plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f'{FIG_DIR}/category_avoidance_heatmap.{ext}', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: category_avoidance_heatmap.pdf/svg")

# ---- Figure 2: Forest plot ----
fig, axes = plt.subplots(1, 3, figsize=(15, 7), sharey=True)

for ax_idx, (methyl_type, methyl_label) in enumerate([
    ('GCCGGC_4mC_T1', 'GCCGGC 4mC'),
    ('AAGCCCG_6mA_T1', 'AAGCCCG 6mA'),
    ('All_4mC_T1', 'All 4mC')
]):
    ax = axes[ax_idx]
    subset = enrichment_df[enrichment_df['methyl_type'] == methyl_type].copy()
    subset = subset.sort_values('odds_ratio')

    y_positions = range(len(subset))

    for i, (_, row) in enumerate(subset.iterrows()):
        color = category_colors.get(row['category'], '#999999')
        ci_low = max(row['OR_CI_low'], 0.01)
        ci_high = min(row['OR_CI_high'], 10)

        # Horizontal line for CI
        ax.plot([ci_low, ci_high], [i, i], color=color, linewidth=2, zorder=2)

        # Point for OR
        marker = 's' if row['p_bonferroni'] < 0.05 else 'o'
        ms = 10 if row['p_bonferroni'] < 0.05 else 7
        ax.plot(row['odds_ratio'], i, marker=marker, color=color, markersize=ms,
                markeredgecolor='black', markeredgewidth=0.5, zorder=3)

    # Reference line at OR=1
    ax.axvline(x=1, color='gray', linestyle='--', linewidth=1, zorder=1)

    ax.set_xscale('log')
    ax.set_xlim(0.1, 10)
    ax.set_xlabel('Odds Ratio (log scale)', fontsize=10)
    ax.set_title(methyl_label, fontsize=11, fontweight='bold')

    if ax_idx == 0:
        y_labels = [category_labels.get(row['category'], row['category']) for _, row in subset.iterrows()]
        ax.set_yticks(list(y_positions))
        ax.set_yticklabels(y_labels, fontsize=9)

    ax.grid(axis='x', alpha=0.3)

fig.suptitle('H23: Forest Plot — Odds Ratio by Gene Category\n(squares = Bonferroni p < 0.05, circles = non-significant)',
             fontsize=13, fontweight='bold', y=1.02)
plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f'{FIG_DIR}/forest_plot_categories.{ext}', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: forest_plot_categories.pdf/svg")

# ---- Figure 3: Expression quintile plot ----
fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharey=True)

quintile_labels = {'Q1_lowest': 'Q1\n(lowest)', 'Q2': 'Q2', 'Q3': 'Q3', 'Q4': 'Q4', 'Q5_highest': 'Q5\n(highest)'}

for ax_idx, (methyl_type, methyl_label) in enumerate([
    ('GCCGGC_4mC_T1', 'GCCGGC 4mC'),
    ('AAGCCCG_6mA_T1', 'AAGCCCG 6mA'),
    ('All_4mC_T1', 'All 4mC')
]):
    ax = axes[ax_idx]
    subset = quintile_df[quintile_df['methyl_type'] == methyl_type]

    fracs = [subset[subset['quintile'] == q]['frac_proximal'].values[0] for q in quintile_order]
    genome_rate = subset['genome_rate'].values[0]
    p_vals = [subset[subset['quintile'] == q]['p_bonferroni'].values[0] for q in quintile_order]

    colors = ['#E63946' if p < 0.05 else '#457B9D' for p in p_vals]
    bars = ax.bar(range(5), fracs, color=colors, edgecolor='black', linewidth=0.5, alpha=0.8)

    # Genome-wide rate reference line
    ax.axhline(y=genome_rate, color='gray', linestyle='--', linewidth=1.5, label=f'Genome-wide: {genome_rate:.3f}')

    # Significance stars
    for i, p in enumerate(p_vals):
        if p < 0.001:
            ax.text(i, fracs[i] + 0.005, '***', ha='center', fontsize=10, fontweight='bold')
        elif p < 0.01:
            ax.text(i, fracs[i] + 0.005, '**', ha='center', fontsize=10, fontweight='bold')
        elif p < 0.05:
            ax.text(i, fracs[i] + 0.005, '*', ha='center', fontsize=10, fontweight='bold')

    ax.set_xticks(range(5))
    ax.set_xticklabels([quintile_labels[q] for q in quintile_order], fontsize=9)
    ax.set_xlabel('Expression quintile (baseMean)', fontsize=10)
    ax.set_title(methyl_label, fontsize=11, fontweight='bold')
    ax.legend(fontsize=8, loc='upper right')

    if ax_idx == 0:
        ax.set_ylabel('Fraction proximal (within 2 kb)', fontsize=10)

    # Add JT test result
    jt_row = jt_df[jt_df['methyl_type'] == methyl_type].iloc[0]
    ax.text(0.02, 0.98, f"JT z={jt_row['z_score']:.2f}\np={jt_row['p_value']:.2e}",
            transform=ax.transAxes, fontsize=8, va='top', ha='left',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

fig.suptitle('H23: Expression Level vs Methylation Proximity', fontsize=13, fontweight='bold')
plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f'{FIG_DIR}/expression_quintile_avoidance.{ext}', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: expression_quintile_avoidance.pdf/svg")

# ---- Figure 4: Comprehensive multi-panel summary ----
fig = plt.figure(figsize=(18, 14))
gs = gridspec.GridSpec(2, 3, hspace=0.35, wspace=0.30)

# Panel A: Heatmap (category × motif)
ax_a = fig.add_subplot(gs[0, 0])
im_a = ax_a.imshow(heatmap_data.values, cmap='RdBu_r', norm=norm, aspect='auto')
cbar_a = plt.colorbar(im_a, ax=ax_a, shrink=0.7)
cbar_a.set_label('Fold enrichment', fontsize=8)
for i in range(heatmap_data.shape[0]):
    for j in range(heatmap_data.shape[1]):
        val = heatmap_data.values[i, j]
        pval = pval_data.values[i, j]
        sig = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else ''
        color = 'white' if abs(val - 1.0) > 0.3 else 'black'
        ax_a.text(j, i, f'{val:.2f}{sig}', ha='center', va='center', fontsize=7, color=color, fontweight='bold')

y_labels_a = [category_labels.get(c, c) for c in cat_order]
ax_a.set_yticks(range(len(cat_order)))
ax_a.set_yticklabels(y_labels_a, fontsize=8)
ax_a.set_xticks(range(3))
ax_a.set_xticklabels(['GCCGGC', 'AAGCCCG', 'All 4mC'], fontsize=8, rotation=30, ha='right')
ax_a.set_title('A. Category × Motif Heatmap', fontsize=10, fontweight='bold')

# Panel B: Forest plot (GCCGGC only)
ax_b = fig.add_subplot(gs[0, 1])
subset_b = enrichment_df[enrichment_df['methyl_type'] == 'GCCGGC_4mC_T1'].sort_values('odds_ratio')
for i, (_, row) in enumerate(subset_b.iterrows()):
    color = category_colors.get(row['category'], '#999')
    ci_low = max(row['OR_CI_low'], 0.01)
    ci_high = min(row['OR_CI_high'], 10)
    ax_b.plot([ci_low, ci_high], [i, i], color=color, linewidth=2.5, zorder=2)
    marker = 's' if row['p_bonferroni'] < 0.05 else 'o'
    ms = 9 if row['p_bonferroni'] < 0.05 else 6
    ax_b.plot(row['odds_ratio'], i, marker=marker, color=color, markersize=ms,
              markeredgecolor='black', markeredgewidth=0.5, zorder=3)
ax_b.axvline(x=1, color='gray', linestyle='--', linewidth=1)
ax_b.set_xscale('log')
ax_b.set_xlim(0.1, 10)
ax_b.set_xlabel('Odds Ratio', fontsize=9)
y_labels_b = [category_labels.get(row['category'], row['category']) for _, row in subset_b.iterrows()]
ax_b.set_yticks(range(len(subset_b)))
ax_b.set_yticklabels(y_labels_b, fontsize=8)
ax_b.set_title('B. Forest Plot (GCCGGC 4mC)', fontsize=10, fontweight='bold')
ax_b.grid(axis='x', alpha=0.3)

# Panel C: Expression quintile (GCCGGC)
ax_c = fig.add_subplot(gs[0, 2])
subset_c = quintile_df[quintile_df['methyl_type'] == 'GCCGGC_4mC_T1']
fracs_c = [subset_c[subset_c['quintile'] == q]['frac_proximal'].values[0] for q in quintile_order]
genome_rate_c = subset_c['genome_rate'].values[0]
p_vals_c = [subset_c[subset_c['quintile'] == q]['p_bonferroni'].values[0] for q in quintile_order]
colors_c = ['#E63946' if p < 0.05 else '#457B9D' for p in p_vals_c]
ax_c.bar(range(5), fracs_c, color=colors_c, edgecolor='black', linewidth=0.5, alpha=0.8)
ax_c.axhline(y=genome_rate_c, color='gray', linestyle='--', linewidth=1.5)
for i, p in enumerate(p_vals_c):
    if p < 0.05:
        sig = '***' if p < 0.001 else '**' if p < 0.01 else '*'
        ax_c.text(i, fracs_c[i] + 0.005, sig, ha='center', fontsize=9, fontweight='bold')
ax_c.set_xticks(range(5))
ax_c.set_xticklabels(['Q1\nlow', 'Q2', 'Q3', 'Q4', 'Q5\nhigh'], fontsize=8)
ax_c.set_ylabel('Fraction proximal', fontsize=9)
ax_c.set_title('C. Expression Quintile (GCCGGC)', fontsize=10, fontweight='bold')
jt_c = jt_df[jt_df['methyl_type'] == 'GCCGGC_4mC_T1'].iloc[0]
ax_c.text(0.02, 0.98, f"JT z={jt_c['z_score']:.2f}, p={jt_c['p_value']:.2e}",
          transform=ax_c.transAxes, fontsize=7, va='top',
          bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

# Panel D: Ribosomal protein comparison
ax_d = fig.add_subplot(gs[1, 0])
ribo_gccggc = ribo_df[ribo_df['methyl_type'] == 'GCCGGC_4mC_T1']
groups_d = ribo_gccggc['group'].values
folds_d = ribo_gccggc['fold_enrichment'].values
colors_d = ['#457B9D', '#2A9D8F', '#E9C46A', '#E63946']
bars_d = ax_d.barh(range(len(groups_d)), folds_d, color=colors_d, edgecolor='black', linewidth=0.5, alpha=0.8)
ax_d.axvline(x=1, color='gray', linestyle='--', linewidth=1)
ax_d.set_yticks(range(len(groups_d)))
ax_d.set_yticklabels(['All ribosomal', '50S subunit', '30S subunit', 'Regulatory'], fontsize=9)
ax_d.set_xlabel('Fold enrichment', fontsize=9)
ax_d.set_title('D. Ribosomal vs Regulatory (GCCGGC)', fontsize=10, fontweight='bold')
for i, (_, row) in enumerate(ribo_gccggc.iterrows()):
    p = row['p_bonferroni']
    sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
    ax_d.text(folds_d[i] + 0.02, i, f'{folds_d[i]:.2f} ({sig})', va='center', fontsize=8)

# Panel E: Core vs Arm stratification for key categories
ax_e = fig.add_subplot(gs[1, 1])
key_cats = ['Regulatory', 'Translation', 'DNA_replication_repair', 'Energy_respiration', 'Hypothetical']
strat_gccggc = stratified_df[(stratified_df['methyl_type'] == 'GCCGGC_4mC_T1') &
                              (stratified_df['category'].isin(key_cats))]

bar_width = 0.35
cat_positions = np.arange(len(key_cats))

for region_idx, (region, color) in enumerate([('core', '#457B9D'), ('arm', '#E9C46A')]):
    region_data = strat_gccggc[strat_gccggc['region'] == region]
    folds_e = []
    for cat in key_cats:
        cat_data = region_data[region_data['category'] == cat]
        if len(cat_data) > 0:
            folds_e.append(cat_data['fold_enrichment'].values[0])
        else:
            folds_e.append(np.nan)

    offset = (region_idx - 0.5) * bar_width
    bars = ax_e.barh(cat_positions + offset, folds_e, bar_width,
                      color=color, edgecolor='black', linewidth=0.5, alpha=0.8,
                      label=region.capitalize())

ax_e.axvline(x=1, color='gray', linestyle='--', linewidth=1)
ax_e.set_yticks(cat_positions)
ax_e.set_yticklabels([category_labels.get(c, c) for c in key_cats], fontsize=8)
ax_e.set_xlabel('Fold enrichment', fontsize=9)
ax_e.set_title('E. Core vs Arm (GCCGGC)', fontsize=10, fontweight='bold')
ax_e.legend(fontsize=8)

# Panel F: Summary text
ax_f = fig.add_subplot(gs[1, 2])
ax_f.axis('off')

# Identify significant depleted categories
sig_depleted = enrichment_df[(enrichment_df['p_bonferroni'] < 0.05) & (enrichment_df['direction'] == 'depleted')]
sig_enriched = enrichment_df[(enrichment_df['p_bonferroni'] < 0.05) & (enrichment_df['direction'] == 'enriched')]

# Count unique depleted categories across all motifs
depleted_cats = sig_depleted['category'].unique()
enriched_cats = sig_enriched['category'].unique()

summary_text = (
    "F. KEY FINDINGS\n\n"
    f"Significant depletion (avoidance):\n"
)
for cat in depleted_cats:
    cat_data = sig_depleted[sig_depleted['category'] == cat]
    n_motifs = len(cat_data)
    min_fold = cat_data['fold_enrichment'].min()
    summary_text += f"  - {category_labels.get(cat, cat)}: {n_motifs}/3 motifs, min fold={min_fold:.2f}\n"

summary_text += f"\nSignificant enrichment:\n"
for cat in enriched_cats:
    cat_data = sig_enriched[sig_enriched['category'] == cat]
    n_motifs = len(cat_data)
    max_fold = cat_data['fold_enrichment'].max()
    summary_text += f"  + {category_labels.get(cat, cat)}: {n_motifs}/3 motifs, max fold={max_fold:.2f}\n"

# JT trend summary
summary_text += f"\nExpression quintile trend:\n"
for _, jt_row in jt_df.iterrows():
    sig = "significant" if jt_row['p_bonferroni'] < 0.05 else "non-significant"
    summary_text += f"  {jt_row['methyl_type']}: {jt_row['direction']} ({sig})\n"

# Interpretation
summary_text += f"\nINTERPRETATION:\n"
if len(depleted_cats) == 1 and 'Regulatory' in depleted_cats:
    summary_text += "Regulatory-SPECIFIC avoidance\n"
    summary_text += "(unique protective pressure on regulatory DNA)"
elif len(depleted_cats) > 1 and 'Regulatory' in depleted_cats:
    reg_fold = enrichment_df[(enrichment_df['category'] == 'Regulatory')]['fold_enrichment'].mean()
    other_folds = enrichment_df[(enrichment_df['category'].isin(depleted_cats)) &
                                 (enrichment_df['category'] != 'Regulatory')]['fold_enrichment'].mean()
    if reg_fold < other_folds:
        summary_text += f"Regulatory strongest (mean fold {reg_fold:.2f}),\n"
        summary_text += f"but not unique (other depleted: mean fold {other_folds:.2f})\n"
        summary_text += "-> Graded avoidance with regulatory peak"
    else:
        summary_text += "Universal avoidance across important gene categories"

ax_f.text(0.05, 0.95, summary_text, transform=ax_f.transAxes, fontsize=8,
          verticalalignment='top', fontfamily='monospace',
          bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

fig.suptitle('H23: Gene Category Specificity of Methylation Avoidance — Comprehensive Summary',
             fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f'{FIG_DIR}/H23_comprehensive_summary.{ext}', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: H23_comprehensive_summary.pdf/svg")

# ============================================================
# 8. Final summary statistics
# ============================================================
print("\n" + "=" * 70)
print("FINAL SUMMARY")
print("=" * 70)

print("\n1. Categories with significant DEPLETION (avoidance, Bonferroni p < 0.05):")
for _, row in sig_depleted.sort_values(['category', 'methyl_type']).iterrows():
    print(f"   {row['category']:<25} {row['methyl_type']:<18} fold={row['fold_enrichment']:.2f} OR={row['odds_ratio']:.2f} [{row['OR_CI_low']:.2f}-{row['OR_CI_high']:.2f}] p_bonf={row['p_bonferroni']:.2e}")

print("\n2. Categories with significant ENRICHMENT (Bonferroni p < 0.05):")
for _, row in sig_enriched.sort_values(['category', 'methyl_type']).iterrows():
    print(f"   {row['category']:<25} {row['methyl_type']:<18} fold={row['fold_enrichment']:.2f} OR={row['odds_ratio']:.2f} [{row['OR_CI_low']:.2f}-{row['OR_CI_high']:.2f}] p_bonf={row['p_bonferroni']:.2e}")

print("\n3. Expression quintile trend:")
for _, jt_row in jt_df.iterrows():
    print(f"   {jt_row['methyl_type']}: JT z={jt_row['z_score']:.3f}, p={jt_row['p_value']:.2e} (Bonf: {jt_row['p_bonferroni']:.2e}), direction={jt_row['direction']}")

print("\n4. Ribosomal protein comparison (GCCGGC):")
for _, row in ribo_df[ribo_df['methyl_type'] == 'GCCGGC_4mC_T1'].iterrows():
    print(f"   {row['group']:<20} n={row['n_genes']:>4} fold={row['fold_enrichment']:.2f} OR={row['odds_ratio']:.2f} p_bonf={row['p_bonferroni']:.2e}")

print("\n5. Geographic robustness (core/arm):")
for cat in ['Regulatory', 'Translation', 'Hypothetical']:
    print(f"   {cat}:")
    for region in ['core', 'arm']:
        subset = stratified_df[(stratified_df['category'] == cat) &
                                (stratified_df['region'] == region) &
                                (stratified_df['methyl_type'] == 'GCCGGC_4mC_T1')]
        if len(subset) > 0:
            row = subset.iloc[0]
            print(f"     {region}: fold={row['fold_enrichment']:.2f} OR={row['odds_ratio']:.2f} p={row['p_value']:.2e}")

print("\n" + "=" * 70)
print("Analysis complete. Output files in:")
print(f"  Figures: {FIG_DIR}/")
print(f"  Tables:  {TBL_DIR}/")
print("=" * 70)
