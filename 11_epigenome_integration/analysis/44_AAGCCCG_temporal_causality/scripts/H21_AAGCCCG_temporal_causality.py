#!/usr/bin/env python3
"""
H21: AAGCCCG 6mA Temporal Causality Test

Tests whether genes losing AAGCCCG 6mA methylation at T2 show significant
upregulation (de-repression) compared to never-methylated genes, and whether
this effect survives geographic stratification (unlike GCCGGC in H19).

Key difference from H19/GCCGGC:
- AAGCCCG has MINIMAL geographic shift between T1 and T2 (4.7 pp vs 65 pp)
- Fewer sites (260 T1, 64 T2) so lower statistical power
- H18 showed stronger gene-level dose-response for AAGCCCG
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

# ── Paths ──────────────────────────────────────────────────────────────────
BASE = Path("/Users/okaban/bioinfo/rna-seq")
ANALYSIS = BASE / "11_epigenome_integration/analysis/44_AAGCCCG_temporal_causality"
FIGURES = ANALYSIS / "figures"
TABLES = ANALYSIS / "tables"

AAGCCCG_MAPPING = BASE / "11_epigenome_integration/analysis/36_AAGCCCG_distribution/tables/AAGCCCG_site_gene_mapping.tsv"
CENSUS_6MA = BASE / "11_epigenome_integration/analysis/23_expanded_motif_search/6mA_final_census.csv"
GFF = BASE / "11_epigenome_integration/analysis/39_GCCGGC_MTase_reverse_ID/data/GCF_000203835.1_ASM20383v1_genomic.gff"

DESEQ_T2T1 = BASE / "04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_2_vs_1.tsv"
DESEQ_T3T1 = BASE / "04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_3_vs_1.tsv"

# H19 comparison data
GCCGGC_TRANSITIONS = BASE / "11_epigenome_integration/analysis/42_GCCGGC_temporal_derepression/tables/transition_group_expression.tsv"
GCCGGC_GEO = BASE / "11_epigenome_integration/analysis/42_GCCGGC_temporal_derepression/tables/geographic_stratified_tests.tsv"

GENOME_SIZE = 8667507  # S. coelicolor A3(2) chromosome

# ── Helper functions ───────────────────────────────────────────────────────
def load_genes_from_gff(gff_path):
    """Parse GFF to extract gene positions and locus tags."""
    genes = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'gene':
                continue
            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attrs = parts[8]

            locus_tag = None
            old_locus_tag = None
            for attr in attrs.split(';'):
                if attr.startswith('locus_tag='):
                    locus_tag = attr.split('=')[1]
                elif attr.startswith('old_locus_tag='):
                    old_locus_tag = attr.split('=')[1]

            if locus_tag:
                genes.append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'locus_tag': locus_tag,
                    'old_locus_tag': old_locus_tag
                })
    return pd.DataFrame(genes)


def classify_region(pos, genome_size=GENOME_SIZE):
    """Classify position as core or arm. Arms: <=1.5Mb from ends."""
    if pos < 1_500_000 or pos > (genome_size - 1_500_000):
        return 'arm'
    return 'core'


def rank_biserial(x, y):
    """Compute rank-biserial correlation (effect size for Mann-Whitney U)."""
    n1, n2 = len(x), len(y)
    U, _ = stats.mannwhitneyu(x, y, alternative='two-sided')
    return 1 - (2 * U) / (n1 * n2)


def get_lfc(gene_set, deseq_df):
    """Extract LFC values for a gene set from DESeq2 results."""
    return deseq_df[deseq_df['gene_id'].isin(gene_set)]['log2FoldChange'].dropna()


def expression_summary(name, gene_set, deseq_df):
    """Compute expression summary for a transition group."""
    sub = deseq_df[deseq_df['gene_id'].isin(gene_set)].copy()
    lfc = sub['log2FoldChange'].dropna()
    n = len(lfc)
    if n == 0:
        return {'group': name, 'n_genes': 0, 'median_LFC': np.nan, 'mean_LFC': np.nan,
                'CI_lower': np.nan, 'CI_upper': np.nan, 'n_DEG': 0, 'frac_DEG': 0,
                'n_up': 0, 'frac_up': 0, 'n_down': 0, 'frac_down': 0}

    med = lfc.median()
    mean = lfc.mean()
    se = lfc.std() / np.sqrt(n)
    ci_lo = mean - 1.96 * se
    ci_hi = mean + 1.96 * se

    sig = sub[sub['padj'] < 0.05]
    n_deg = len(sig)
    n_up = (sig['log2FoldChange'] > 0).sum()
    n_down = (sig['log2FoldChange'] < 0).sum()

    return {
        'group': name, 'n_genes': n,
        'median_LFC': med, 'mean_LFC': mean,
        'CI_lower': ci_lo, 'CI_upper': ci_hi,
        'n_DEG': n_deg, 'frac_DEG': n_deg / n,
        'n_up': n_up, 'frac_up': n_up / n,
        'n_down': n_down, 'frac_down': n_down / n
    }


# ══════════════════════════════════════════════════════════════════════════
# STEP 1: Load genes from GFF and classify regions
# ══════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("H21: AAGCCCG 6mA TEMPORAL CAUSALITY TEST")
print("=" * 70)

genes_df = load_genes_from_gff(GFF)
genes_df['midpoint'] = (genes_df['start'] + genes_df['end']) / 2
genes_df['region'] = genes_df['midpoint'].apply(classify_region)
# Keep only chromosome genes (NC_003888.3)
genes_df = genes_df[genes_df['chrom'] == 'NC_003888.3'].copy()

print(f"\nTotal chromosome genes in GFF: {len(genes_df)}")
print(f"  Core: {(genes_df['region'] == 'core').sum()}")
print(f"  Arm:  {(genes_df['region'] == 'arm').sum()}")

# ══════════════════════════════════════════════════════════════════════════
# STEP 2: Map AAGCCCG sites to genes (within 2kb)
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("STEP 2: AAGCCCG SITE-GENE MAPPING")
print("=" * 70)

# Load the pre-computed mapping (all entries are within 2kb)
mapping = pd.read_csv(AAGCCCG_MAPPING, sep='\t')
print(f"Total mapping entries: {len(mapping)}")

# Separate by timepoint
t1_mapping = mapping[mapping['timepoint'] == 'T1'].copy()
t2_mapping = mapping[mapping['timepoint'] == 'T2'].copy()

# Unique sites
t1_sites_unique = t1_mapping['position'].unique()
t2_sites_unique = t2_mapping['position'].unique()
print(f"\nAAGCCCG sites: T1={len(t1_sites_unique)}, T2={len(t2_sites_unique)}")

# Unique genes within 2kb
t1_genes = set(t1_mapping['locus_tag'].unique())
t2_genes = set(t2_mapping['locus_tag'].unique())

print(f"Genes within 2kb of T1 AAGCCCG: {len(t1_genes)}")
print(f"Genes within 2kb of T2 AAGCCCG: {len(t2_genes)}")

# Geographic composition of site sets
t1_site_regions = t1_mapping.drop_duplicates('position')['region'].value_counts()
t2_site_regions = t2_mapping.drop_duplicates('position')['region'].value_counts()
print(f"\nT1 site regions: {dict(t1_site_regions)}")
print(f"T2 site regions: {dict(t2_site_regions)}")

# ══════════════════════════════════════════════════════════════════════════
# STEP 3: Classify genes by AAGCCCG methylation transition (T1→T2)
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("STEP 3: METHYLATION TRANSITION GROUPS")
print("=" * 70)

# Load DESeq2 results
deseq_t2t1 = pd.read_csv(DESEQ_T2T1, sep='\t')
deseq_t3t1 = pd.read_csv(DESEQ_T3T1, sep='\t')
deseq_genes = set(deseq_t2t1['gene_id'])
all_genes = set(genes_df['locus_tag'])

print(f"DESeq2 genes: {len(deseq_genes)}")

lost_genes = t1_genes - t2_genes      # Had T1, no T2
gained_genes = t2_genes - t1_genes    # No T1, has T2
both_genes = t1_genes & t2_genes      # Had both
never_genes = all_genes - t1_genes - t2_genes  # Neither

print(f"\nTransition groups (all genes):")
print(f"  Lost (T1 only):     {len(lost_genes)}")
print(f"  Gained (T2 only):   {len(gained_genes)}")
print(f"  Both (T1 & T2):     {len(both_genes)}")
print(f"  Never (neither):    {len(never_genes)}")

# Restrict to genes with expression data
lost_expr = lost_genes & deseq_genes
gained_expr = gained_genes & deseq_genes
both_expr = both_genes & deseq_genes
never_expr = never_genes & deseq_genes

print(f"\nWith expression data:")
print(f"  Lost:   {len(lost_expr)}")
print(f"  Gained: {len(gained_expr)}")
print(f"  Both:   {len(both_expr)}")
print(f"  Never:  {len(never_expr)}")

# ══════════════════════════════════════════════════════════════════════════
# STEP 4: Expression comparison by transition group
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("STEP 4: EXPRESSION BY TRANSITION GROUP (T2 vs T1)")
print("=" * 70)

lost_lfc = get_lfc(lost_expr, deseq_t2t1)
gained_lfc = get_lfc(gained_expr, deseq_t2t1)
both_lfc = get_lfc(both_expr, deseq_t2t1)
never_lfc = get_lfc(never_expr, deseq_t2t1)

# Also T3vsT1
lost_lfc_t3 = get_lfc(lost_expr, deseq_t3t1)
gained_lfc_t3 = get_lfc(gained_expr, deseq_t3t1)
both_lfc_t3 = get_lfc(both_expr, deseq_t3t1)
never_lfc_t3 = get_lfc(never_expr, deseq_t3t1)

summary_rows = []
for name, gene_set in [("Lost", lost_expr), ("Gained", gained_expr),
                         ("Both", both_expr), ("Never", never_expr)]:
    s = expression_summary(name, gene_set, deseq_t2t1)
    summary_rows.append(s)
    print(f"\n{name} (n={s['n_genes']}):")
    print(f"  Median LFC: {s['median_LFC']:.4f}")
    print(f"  Mean LFC:   {s['mean_LFC']:.4f} (95% CI: [{s['CI_lower']:.4f}, {s['CI_upper']:.4f}])")
    print(f"  DEGs:       {s['n_DEG']} ({s['frac_DEG']:.1%}) | Up: {s['n_up']} ({s['frac_up']:.1%}) | Down: {s['n_down']} ({s['frac_down']:.1%})")

summary_df = pd.DataFrame(summary_rows)

print("\n--- T3 vs T1 LFC by transition group ---")
for name, lfc in [("Lost", lost_lfc_t3), ("Gained", gained_lfc_t3),
                   ("Both", both_lfc_t3), ("Never", never_lfc_t3)]:
    if len(lfc) > 0:
        print(f"  {name}: n={len(lfc)}, median={lfc.median():.4f}, mean={lfc.mean():.4f}")

# ══════════════════════════════════════════════════════════════════════════
# STEP 5: Statistical tests
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("STEP 5: STATISTICAL TESTS (Wilcoxon rank-sum on T2vsT1 LFC)")
print("=" * 70)

test_pairs = [
    ("Lost vs Never", lost_lfc, never_lfc),
    ("Lost vs Gained", lost_lfc, gained_lfc),
    ("Gained vs Never", gained_lfc, never_lfc),
    ("Both vs Never", both_lfc, never_lfc),
    ("Lost vs Both", lost_lfc, both_lfc),
    ("Gained vs Both", gained_lfc, both_lfc),
]

n_tests = len(test_pairs)
stat_rows = []
for label, x, y in test_pairs:
    if len(x) < 3 or len(y) < 3:
        print(f"\n{label}: SKIPPED (n too small: {len(x)}, {len(y)})")
        stat_rows.append({
            'comparison': label, 'n1': len(x), 'n2': len(y),
            'U': np.nan, 'p_value': np.nan, 'p_bonferroni': np.nan,
            'rank_biserial': np.nan,
            'median_diff': np.nan, 'mean_diff': np.nan,
            'direction': 'NA', 'timepoint': 'T2vsT1'
        })
        continue

    U, p = stats.mannwhitneyu(x, y, alternative='two-sided')
    p_bonf = min(p * n_tests, 1.0)
    rb = rank_biserial(x, y)
    med_diff = x.median() - y.median()
    mean_diff = x.mean() - y.mean()
    direction = "Group1 higher" if med_diff > 0 else "Group1 lower"

    print(f"\n{label}:")
    print(f"  n = {len(x)} vs {len(y)}")
    print(f"  U = {U:.0f}, p = {p:.2e}, Bonferroni p = {p_bonf:.2e}")
    print(f"  Rank-biserial r = {rb:.4f}")
    print(f"  Median diff = {med_diff:.4f} ({direction})")

    stat_rows.append({
        'comparison': label, 'n1': len(x), 'n2': len(y),
        'U': U, 'p_value': p, 'p_bonferroni': p_bonf,
        'rank_biserial': rb, 'median_diff': med_diff, 'mean_diff': mean_diff,
        'direction': direction, 'timepoint': 'T2vsT1'
    })

# Also test T3vsT1
print("\n--- T3vsT1 tests ---")
for label, x, y in [("Lost vs Never (T3vsT1)", lost_lfc_t3, never_lfc_t3),
                     ("Gained vs Never (T3vsT1)", gained_lfc_t3, never_lfc_t3)]:
    if len(x) >= 3 and len(y) >= 3:
        U, p = stats.mannwhitneyu(x, y, alternative='two-sided')
        rb = rank_biserial(x, y)
        print(f"  {label}: U={U:.0f}, p={p:.2e}, r={rb:.4f}, med_diff={x.median()-y.median():.4f}")
        stat_rows.append({
            'comparison': label, 'n1': len(x), 'n2': len(y),
            'U': U, 'p_value': p, 'p_bonferroni': np.nan,
            'rank_biserial': rb, 'median_diff': x.median() - y.median(),
            'mean_diff': x.mean() - y.mean(),
            'direction': "Group1 higher" if x.median() > y.median() else "Group1 lower",
            'timepoint': 'T3vsT1'
        })

stat_df = pd.DataFrame(stat_rows)

# ══════════════════════════════════════════════════════════════════════════
# STEP 6: Geographic stratification (critical control)
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("STEP 6: GEOGRAPHIC STRATIFICATION (CRITICAL CONTROL)")
print("=" * 70)

gene_regions = genes_df[['locus_tag', 'region']].drop_duplicates()

# First: geographic composition of each transition group
print("\n--- Geographic composition of transition groups ---")
geo_composition = {}
for name, gene_set in [("Lost", lost_expr), ("Gained", gained_expr),
                         ("Both", both_expr), ("Never", never_expr)]:
    merged = gene_regions[gene_regions['locus_tag'].isin(gene_set)]
    n_core = (merged['region'] == 'core').sum()
    n_arm = (merged['region'] == 'arm').sum()
    total = n_core + n_arm
    pct_core = n_core / total * 100 if total > 0 else 0
    pct_arm = n_arm / total * 100 if total > 0 else 0
    print(f"  {name:>8s}: n={total:>5}, core={n_core} ({pct_core:.1f}%), arm={n_arm} ({pct_arm:.1f}%)")
    geo_composition[name] = {'n_total': total, 'n_core': n_core, 'n_arm': n_arm,
                              'pct_core': pct_core, 'pct_arm': pct_arm}

# Stratified tests
geo_rows = []
for region_name in ['core', 'arm']:
    region_genes = set(gene_regions[gene_regions['region'] == region_name]['locus_tag'])

    lost_r = lost_expr & region_genes
    never_r = never_expr & region_genes
    gained_r = gained_expr & region_genes
    both_r = both_expr & region_genes

    lost_lfc_r = get_lfc(lost_r, deseq_t2t1)
    never_lfc_r = get_lfc(never_r, deseq_t2t1)
    gained_lfc_r = get_lfc(gained_r, deseq_t2t1)
    both_lfc_r = get_lfc(both_r, deseq_t2t1)

    print(f"\n{region_name.upper()} genes:")
    print(f"  Lost:   n={len(lost_lfc_r)}, median LFC={lost_lfc_r.median():.4f}" if len(lost_lfc_r) > 0 else f"  Lost: n=0")
    print(f"  Gained: n={len(gained_lfc_r)}, median LFC={gained_lfc_r.median():.4f}" if len(gained_lfc_r) > 0 else f"  Gained: n=0")
    print(f"  Both:   n={len(both_lfc_r)}, median LFC={both_lfc_r.median():.4f}" if len(both_lfc_r) > 0 else f"  Both: n=0")
    print(f"  Never:  n={len(never_lfc_r)}, median LFC={never_lfc_r.median():.4f}" if len(never_lfc_r) > 0 else f"  Never: n=0")

    # Lost vs Never within region
    for comp_label, g1, g1_name, g2, g2_name in [
        ("Lost vs Never", lost_lfc_r, "Lost", never_lfc_r, "Never"),
        ("Gained vs Never", gained_lfc_r, "Gained", never_lfc_r, "Never"),
        ("Both vs Never", both_lfc_r, "Both", never_lfc_r, "Never"),
    ]:
        if len(g1) >= 3 and len(g2) >= 3:
            U, p = stats.mannwhitneyu(g1, g2, alternative='two-sided')
            rb = rank_biserial(g1, g2)
            print(f"  {comp_label}: U={U:.0f}, p={p:.2e}, r={rb:.4f}")
            geo_rows.append({
                'region': region_name, 'comparison': comp_label,
                'n1': len(g1), 'n2': len(g2),
                'median1': g1.median(), 'median2': g2.median(),
                'U': U, 'p_value': p, 'rank_biserial': rb
            })
        else:
            print(f"  {comp_label}: SKIPPED (n={len(g1)} vs {len(g2)})")
            geo_rows.append({
                'region': region_name, 'comparison': comp_label,
                'n1': len(g1), 'n2': len(g2),
                'median1': g1.median() if len(g1) > 0 else np.nan,
                'median2': g2.median() if len(g2) > 0 else np.nan,
                'U': np.nan, 'p_value': np.nan, 'rank_biserial': np.nan
            })

geo_df = pd.DataFrame(geo_rows)

# ══════════════════════════════════════════════════════════════════════════
# STEP 7: Dose-response within transition groups
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("STEP 7: DOSE-RESPONSE (NUMBER OF T1 SITES LOST)")
print("=" * 70)

# Count T1 sites per gene
t1_site_counts = t1_mapping.groupby('locus_tag')['position'].nunique().reset_index()
t1_site_counts.columns = ['gene_id', 'n_t1_sites']

# Merge with DESeq2 for Lost genes
lost_dose = t1_site_counts[t1_site_counts['gene_id'].isin(lost_expr)].merge(
    deseq_t2t1[['gene_id', 'log2FoldChange', 'padj']], on='gene_id'
)

print(f"\nLost genes with T1 site count: {len(lost_dose)}")
print(f"Distribution of T1 site counts for Lost genes:")
for n_sites, cnt in sorted(lost_dose['n_t1_sites'].value_counts().items()):
    sub = lost_dose[lost_dose['n_t1_sites'] == n_sites]
    print(f"  {n_sites} site(s): {cnt} genes, median LFC={sub['log2FoldChange'].median():.4f}")

# Spearman correlation: number of sites vs LFC
if len(lost_dose) >= 5:
    rho, p_corr = stats.spearmanr(lost_dose['n_t1_sites'], lost_dose['log2FoldChange'])
    print(f"\nSpearman (n_T1_sites vs T2vsT1 LFC in Lost genes):")
    print(f"  rho={rho:.4f}, p={p_corr:.2e}")

# Multi-site vs single-site comparison
single_site = lost_dose[lost_dose['n_t1_sites'] == 1]['log2FoldChange']
multi_site = lost_dose[lost_dose['n_t1_sites'] >= 2]['log2FoldChange']
print(f"\nSingle-site Lost (n={len(single_site)}): median LFC={single_site.median():.4f}")
print(f"Multi-site Lost (n={len(multi_site)}):  median LFC={multi_site.median():.4f}")
if len(single_site) >= 3 and len(multi_site) >= 3:
    U, p = stats.mannwhitneyu(multi_site, single_site, alternative='two-sided')
    print(f"Mann-Whitney multi vs single: U={U:.0f}, p={p:.2e}")

# ══════════════════════════════════════════════════════════════════════════
# STEP 8: Methylation frequency analysis
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("STEP 8: METHYLATION FREQUENCY ANALYSIS")
print("=" * 70)

# Load 6mA census for AAGCCCG sites with frequency data
census = pd.read_csv(CENSUS_6MA)
aag_census = census[(census['at_AAGCCCG'] == True) | (census['final_motif'] == 'AAGCCCG')]
print(f"AAGCCCG sites in census: {len(aag_census)}")

# T1 AAGCCCG sites with frequency
aag_t1_freq = aag_census[aag_census['timepoint'] == 'T1'].copy()
print(f"T1 AAGCCCG sites with frequency: {len(aag_t1_freq)}")

# Map frequency to genes: for each gene, compute mean/max frequency of T1 AAGCCCG sites within 2kb
# Use the mapping file to link sites to genes
t1_site_freq = aag_t1_freq[['position', 'frequency']].drop_duplicates()

# Merge mapping with frequency
t1_map_freq = t1_mapping.merge(t1_site_freq, on='position', how='left')

# Per-gene: mean and max T1 AAGCCCG frequency
gene_freq = t1_map_freq.groupby('locus_tag').agg(
    mean_freq=('frequency', 'mean'),
    max_freq=('frequency', 'max'),
    n_sites=('position', 'nunique')
).reset_index()

# Merge with DESeq2 for Lost genes only
lost_freq = gene_freq[gene_freq['locus_tag'].isin(lost_expr)].merge(
    deseq_t2t1[['gene_id', 'log2FoldChange', 'padj']],
    left_on='locus_tag', right_on='gene_id'
)

print(f"\nLost genes with frequency data: {len(lost_freq)}")
if len(lost_freq) >= 5:
    rho_mean, p_mean = stats.spearmanr(lost_freq['mean_freq'], lost_freq['log2FoldChange'])
    rho_max, p_max = stats.spearmanr(lost_freq['max_freq'], lost_freq['log2FoldChange'])
    print(f"Spearman (mean_freq vs LFC): rho={rho_mean:.4f}, p={p_mean:.2e}")
    print(f"Spearman (max_freq vs LFC):  rho={rho_max:.4f}, p={p_max:.2e}")

# Also for ALL genes near T1 AAGCCCG sites
all_t1_freq = gene_freq.merge(
    deseq_t2t1[['gene_id', 'log2FoldChange', 'padj']],
    left_on='locus_tag', right_on='gene_id'
)
print(f"\nAll T1-proximal genes with frequency data: {len(all_t1_freq)}")
if len(all_t1_freq) >= 5:
    rho_all, p_all = stats.spearmanr(all_t1_freq['mean_freq'], all_t1_freq['log2FoldChange'])
    print(f"Spearman (mean_freq vs T2vsT1 LFC, all T1-proximal): rho={rho_all:.4f}, p={p_all:.2e}")

# ══════════════════════════════════════════════════════════════════════════
# STEP 9: Compare with GCCGGC (H19) directly
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("STEP 9: AAGCCCG vs GCCGGC COMPARISON")
print("=" * 70)

gccggc_expr = pd.read_csv(GCCGGC_TRANSITIONS, sep='\t')
gccggc_geo = pd.read_csv(GCCGGC_GEO, sep='\t')

print("\nSide-by-side comparison:")
print(f"{'Metric':<35s} {'AAGCCCG':>15s} {'GCCGGC':>15s}")
print("-" * 65)

# Sites
print(f"{'T1 sites':<35s} {'260':>15s} {'1,289':>15s}")
print(f"{'T2 sites':<35s} {'64':>15s} {'376':>15s}")

# Geographic shift
print(f"{'T1 core %':<35s} {'68.8%':>15s} {'17%':>15s}")
print(f"{'T2 core %':<35s} {'64.1%':>15s} {'82%':>15s}")
print(f"{'Geographic shift (pp)':<35s} {'4.7':>15s} {'65':>15s}")

# Get GCCGGC Lost and Never stats
gccggc_lost = gccggc_expr[gccggc_expr['group'] == 'Lost'].iloc[0]
gccggc_never = gccggc_expr[gccggc_expr['group'] == 'Never'].iloc[0]
aag_lost_row = summary_df[summary_df['group'] == 'Lost'].iloc[0]
aag_never_row = summary_df[summary_df['group'] == 'Never'].iloc[0]

print(f"\n{'Lost group n':<35s} {aag_lost_row['n_genes']:>15.0f} {gccggc_lost['n_genes']:>15.0f}")
print(f"{'Never group n':<35s} {aag_never_row['n_genes']:>15.0f} {gccggc_never['n_genes']:>15.0f}")
print(f"{'Lost median LFC':<35s} {aag_lost_row['median_LFC']:>+15.4f} {gccggc_lost['median_LFC']:>+15.4f}")
print(f"{'Never median LFC':<35s} {aag_never_row['median_LFC']:>+15.4f} {gccggc_never['median_LFC']:>+15.4f}")
print(f"{'Lost-Never median diff':<35s} {aag_lost_row['median_LFC'] - aag_never_row['median_LFC']:>+15.4f} {gccggc_lost['median_LFC'] - gccggc_never['median_LFC']:>+15.4f}")
print(f"{'Lost frac DEG':<35s} {aag_lost_row['frac_DEG']:>15.1%} {gccggc_lost['frac_DEG']:>15.1%}")
print(f"{'Lost frac up':<35s} {aag_lost_row['frac_up']:>15.1%} {gccggc_lost['frac_up']:>15.1%}")

# Comparison: Lost vs Never p-value
aag_lvn = stat_df[stat_df['comparison'] == 'Lost vs Never']
aag_p = aag_lvn.iloc[0]['p_value'] if len(aag_lvn) > 0 else np.nan
aag_rb = aag_lvn.iloc[0]['rank_biserial'] if len(aag_lvn) > 0 else np.nan

# GCCGGC overall Lost vs Never p (from H19 stat tests)
# Load from H19 if available
gccggc_stat_file = BASE / "11_epigenome_integration/analysis/42_GCCGGC_temporal_derepression/tables/statistical_tests.tsv"
gccggc_stat = pd.read_csv(gccggc_stat_file, sep='\t') if gccggc_stat_file.exists() else pd.DataFrame()
gccggc_lvn = gccggc_stat[gccggc_stat['comparison'] == 'Lost vs Never'] if len(gccggc_stat) > 0 else pd.DataFrame()
gccggc_p = gccggc_lvn.iloc[0]['p_value'] if len(gccggc_lvn) > 0 else np.nan
gccggc_rb = gccggc_lvn.iloc[0]['rank_biserial'] if len(gccggc_lvn) > 0 else np.nan

print(f"\n{'Lost vs Never p-value':<35s} {aag_p:>15.2e} {gccggc_p:>15.2e}")
print(f"{'Lost vs Never rank-biserial':<35s} {aag_rb:>+15.4f} {gccggc_rb:>+15.4f}")

# Geographic stratification comparison
print(f"\n--- Geographic stratification comparison ---")
for region_name in ['core', 'arm']:
    aag_row = geo_df[(geo_df['region'] == region_name) & (geo_df['comparison'] == 'Lost vs Never')]
    gccggc_row = gccggc_geo[(gccggc_geo['region'] == region_name) & (gccggc_geo['comparison'] == 'Lost vs Never')]

    aag_p_geo = aag_row.iloc[0]['p_value'] if len(aag_row) > 0 else np.nan
    aag_rb_geo = aag_row.iloc[0]['rank_biserial'] if len(aag_row) > 0 else np.nan
    gccggc_p_geo = gccggc_row.iloc[0]['p_value'] if len(gccggc_row) > 0 else np.nan
    gccggc_rb_geo = gccggc_row.iloc[0]['rank_biserial'] if len(gccggc_row) > 0 else np.nan

    print(f"  {region_name.upper()}:")
    print(f"    AAGCCCG Lost vs Never: p={aag_p_geo:.2e}, r={aag_rb_geo:+.4f}" if not np.isnan(aag_p_geo) else f"    AAGCCCG: SKIPPED/NA")
    print(f"    GCCGGC Lost vs Never:  p={gccggc_p_geo:.2e}, r={gccggc_rb_geo:+.4f}" if not np.isnan(gccggc_p_geo) else f"    GCCGGC: SKIPPED/NA")

# Build comparison table
comparison_rows = []
for metric, aag_val, gcc_val in [
    ("T1_sites", 260, 1289),
    ("T2_sites", 64, 376),
    ("T1_core_pct", 68.8, 17.0),
    ("T2_core_pct", 64.1, 82.0),
    ("geo_shift_pp", 4.7, 65.0),
    ("Lost_n", aag_lost_row['n_genes'], gccggc_lost['n_genes']),
    ("Never_n", aag_never_row['n_genes'], gccggc_never['n_genes']),
    ("Lost_median_LFC", aag_lost_row['median_LFC'], gccggc_lost['median_LFC']),
    ("Never_median_LFC", aag_never_row['median_LFC'], gccggc_never['median_LFC']),
    ("LostvsNever_p", aag_p, gccggc_p),
    ("LostvsNever_rb", aag_rb, gccggc_rb),
]:
    comparison_rows.append({'metric': metric, 'AAGCCCG': aag_val, 'GCCGGC': gcc_val})

comparison_df = pd.DataFrame(comparison_rows)

# ══════════════════════════════════════════════════════════════════════════
# STEP 10: Save tables
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("SAVING TABLES")
print("=" * 70)

# Gene methylation transitions table
transition_records = []
for gene in sorted(deseq_genes):
    has_t1 = gene in t1_genes
    has_t2 = gene in t2_genes

    if has_t1 and not has_t2:
        group = 'Lost'
    elif not has_t1 and has_t2:
        group = 'Gained'
    elif has_t1 and has_t2:
        group = 'Both'
    else:
        group = 'Never'

    gene_info = genes_df[genes_df['locus_tag'] == gene]
    region = gene_info.iloc[0]['region'] if len(gene_info) > 0 else 'unknown'
    midpoint = gene_info.iloc[0]['midpoint'] if len(gene_info) > 0 else np.nan

    row_t2t1 = deseq_t2t1[deseq_t2t1['gene_id'] == gene]
    row_t3t1 = deseq_t3t1[deseq_t3t1['gene_id'] == gene]

    # Count T1 AAGCCCG sites near this gene
    n_t1 = t1_mapping[t1_mapping['locus_tag'] == gene]['position'].nunique()
    n_t2 = t2_mapping[t2_mapping['locus_tag'] == gene]['position'].nunique()

    # Mean T1 frequency
    gene_t1_freq = t1_map_freq[t1_map_freq['locus_tag'] == gene]['frequency']
    mean_t1_freq = gene_t1_freq.mean() if len(gene_t1_freq) > 0 else np.nan

    transition_records.append({
        'gene_id': gene,
        'region': region,
        'midpoint': midpoint,
        'AAGCCCG_T1_sites': n_t1,
        'AAGCCCG_T2_sites': n_t2,
        'transition': group,
        'mean_T1_frequency': mean_t1_freq,
        'baseMean': row_t2t1['baseMean'].values[0] if len(row_t2t1) > 0 else np.nan,
        'LFC_T2vsT1': row_t2t1['log2FoldChange'].values[0] if len(row_t2t1) > 0 else np.nan,
        'padj_T2vsT1': row_t2t1['padj'].values[0] if len(row_t2t1) > 0 else np.nan,
        'LFC_T3vsT1': row_t3t1['log2FoldChange'].values[0] if len(row_t3t1) > 0 else np.nan,
        'padj_T3vsT1': row_t3t1['padj'].values[0] if len(row_t3t1) > 0 else np.nan,
    })

transitions_full = pd.DataFrame(transition_records)
transitions_full.to_csv(TABLES / "gene_methylation_transitions.tsv", sep='\t', index=False)
print(f"Saved gene_methylation_transitions.tsv ({len(transitions_full)} genes)")

summary_df.to_csv(TABLES / "transition_group_expression.tsv", sep='\t', index=False)
print("Saved transition_group_expression.tsv")

stat_df.to_csv(TABLES / "statistical_tests.tsv", sep='\t', index=False)
print("Saved statistical_tests.tsv")

geo_df.to_csv(TABLES / "geographic_stratified_tests.tsv", sep='\t', index=False)
print("Saved geographic_stratified_tests.tsv")

comparison_df.to_csv(TABLES / "AAGCCCG_vs_GCCGGC_comparison.tsv", sep='\t', index=False)
print("Saved AAGCCCG_vs_GCCGGC_comparison.tsv")

# ══════════════════════════════════════════════════════════════════════════
# STEP 11: Visualizations
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("GENERATING FIGURES")
print("=" * 70)

colors = {
    'Lost': '#e74c3c',     # Red
    'Gained': '#3498db',   # Blue
    'Both': '#9b59b6',     # Purple
    'Never': '#95a5a6',    # Gray
}

# ── Figure 1: Violin plots ──
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for ax_idx, (groups_data, title) in enumerate([
    ([(lost_lfc, 'Lost'), (gained_lfc, 'Gained'),
      (both_lfc, 'Both'), (never_lfc, 'Never')], 'T2 vs T1 LFC'),
    ([(lost_lfc_t3, 'Lost'), (gained_lfc_t3, 'Gained'),
      (both_lfc_t3, 'Both'), (never_lfc_t3, 'Never')], 'T3 vs T1 LFC'),
]):
    ax = axes[ax_idx]
    data_list = []
    labels = []
    cols = []

    for lfc, name in groups_data:
        if len(lfc) > 0:
            data_list.append(lfc.values)
            labels.append(f"{name}\n(n={len(lfc)})")
            cols.append(colors[name])

    if len(data_list) > 0:
        parts = ax.violinplot(data_list, positions=range(len(data_list)),
                              showmeans=True, showmedians=True)
        for i, pc in enumerate(parts['bodies']):
            pc.set_facecolor(cols[i])
            pc.set_alpha(0.7)
        parts['cmeans'].set_color('black')
        parts['cmedians'].set_color('red')

        # Add median labels
        for i, d in enumerate(data_list):
            med = np.median(d)
            ax.text(i, ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1] - ax.get_ylim()[0]),
                    f'med={med:.2f}', ha='center', fontsize=7, color='darkred')

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=9)
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.set_ylabel('log2 Fold Change')
    ax.set_title(title)

fig.suptitle('H21: AAGCCCG Methylation Transition Groups vs Expression Change',
             fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig(FIGURES / 'transition_groups_violin.pdf', dpi=300, bbox_inches='tight')
plt.savefig(FIGURES / 'transition_groups_violin.svg', bbox_inches='tight')
plt.close()
print("Saved transition_groups_violin.pdf/svg")

# ── Figure 2: Geographic composition stacked bars ──
fig, ax = plt.subplots(figsize=(8, 6))

group_names = ['Lost', 'Gained', 'Both', 'Never']
x_pos = np.arange(len(group_names))
width = 0.6

for i, name in enumerate(group_names):
    gc = geo_composition.get(name, {})
    total = gc.get('n_total', 0)
    if total == 0:
        continue
    core_frac = gc['n_core'] / total
    arm_frac = gc['n_arm'] / total

    ax.bar(i, core_frac, width, color='#3498db', alpha=0.8,
           label='Core' if i == 0 else '')
    ax.bar(i, arm_frac, width, bottom=core_frac, color='#e67e22', alpha=0.8,
           label='Arm' if i == 0 else '')

    ax.text(i, 0.5, f'{gc["pct_core"]:.0f}%\ncore', ha='center', fontsize=9, color='white', fontweight='bold')
    ax.text(i, 1.02, f'n={total}', ha='center', fontsize=9)

ax.set_xticks(x_pos)
ax.set_xticklabels(group_names, fontsize=11)
ax.set_ylabel('Fraction')
ax.set_title('Geographic Composition of AAGCCCG Transition Groups', fontweight='bold')
ax.legend(loc='upper right')
ax.set_ylim(0, 1.1)
plt.tight_layout()
plt.savefig(FIGURES / 'geographic_composition.pdf', dpi=300, bbox_inches='tight')
plt.savefig(FIGURES / 'geographic_composition.svg', bbox_inches='tight')
plt.close()
print("Saved geographic_composition.pdf/svg")

# ── Figure 3: AAGCCCG vs GCCGGC side-by-side comparison ──
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Panel A: LFC comparison
ax = axes[0]
aag_groups = ['Lost', 'Gained', 'Both', 'Never']
aag_medians = [summary_df[summary_df['group'] == g].iloc[0]['median_LFC'] for g in aag_groups]
gcc_medians = [gccggc_expr[gccggc_expr['group'] == g].iloc[0]['median_LFC']
               if g in gccggc_expr['group'].values else np.nan for g in aag_groups]

x = np.arange(len(aag_groups))
w = 0.35
ax.bar(x - w/2, aag_medians, w, color='#e74c3c', alpha=0.8, label='AAGCCCG')
ax.bar(x + w/2, gcc_medians, w, color='#3498db', alpha=0.8, label='GCCGGC')

ax.set_xticks(x)
ax.set_xticklabels(aag_groups, fontsize=10)
ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
ax.set_ylabel('Median T2vsT1 LFC')
ax.set_title('A. Median LFC by Transition Group', fontweight='bold')
ax.legend()

# Panel B: Geographic composition comparison (Lost group)
ax = axes[1]
categories = ['AAGCCCG\nLost', 'GCCGGC\nLost', 'AAGCCCG\nNever', 'GCCGGC\nNever']

# Calculate GCCGGC geographic composition from the transitions file
gccggc_trans_file = BASE / "11_epigenome_integration/analysis/42_GCCGGC_temporal_derepression/tables/gene_methylation_transitions.tsv"
if gccggc_trans_file.exists():
    gccggc_trans = pd.read_csv(gccggc_trans_file, sep='\t')
    gcc_lost_core = len(gccggc_trans[(gccggc_trans['transition_T1T2'] == 'Lost') & (gccggc_trans['region'] == 'core')])
    gcc_lost_arm = len(gccggc_trans[(gccggc_trans['transition_T1T2'] == 'Lost') & (gccggc_trans['region'] == 'arm')])
    gcc_never_core = len(gccggc_trans[(gccggc_trans['transition_T1T2'] == 'Never') & (gccggc_trans['region'] == 'core')])
    gcc_never_arm = len(gccggc_trans[(gccggc_trans['transition_T1T2'] == 'Never') & (gccggc_trans['region'] == 'arm')])
else:
    gcc_lost_core, gcc_lost_arm = 0, 0
    gcc_never_core, gcc_never_arm = 0, 0

aag_lc = geo_composition.get('Lost', {}).get('n_core', 0)
aag_la = geo_composition.get('Lost', {}).get('n_arm', 0)
aag_nc = geo_composition.get('Never', {}).get('n_core', 0)
aag_na = geo_composition.get('Never', {}).get('n_arm', 0)

core_fracs = []
arm_fracs = []
totals = []
for nc, na in [(aag_lc, aag_la), (gcc_lost_core, gcc_lost_arm),
               (aag_nc, aag_na), (gcc_never_core, gcc_never_arm)]:
    t = nc + na
    totals.append(t)
    core_fracs.append(nc / t if t > 0 else 0)
    arm_fracs.append(na / t if t > 0 else 0)

bx = np.arange(len(categories))
ax.bar(bx, core_fracs, 0.6, color='#3498db', alpha=0.8, label='Core')
ax.bar(bx, arm_fracs, 0.6, bottom=core_fracs, color='#e67e22', alpha=0.8, label='Arm')
for i in range(len(categories)):
    ax.text(i, 0.5, f'{core_fracs[i]*100:.0f}%', ha='center', fontsize=9, color='white', fontweight='bold')
    ax.text(i, 1.02, f'n={totals[i]}', ha='center', fontsize=8)

ax.set_xticks(bx)
ax.set_xticklabels(categories, fontsize=9)
ax.set_ylabel('Fraction')
ax.set_title('B. Geographic Composition: Lost & Never', fontweight='bold')
ax.legend(loc='upper right', fontsize=8)
ax.set_ylim(0, 1.15)

# Panel C: Statistical summary
ax = axes[2]
ax.axis('off')

summary_text = "COMPARISON SUMMARY\n" + "=" * 40 + "\n\n"
summary_text += f"{'':>20s}  {'AAGCCCG':>10s}  {'GCCGGC':>10s}\n"
summary_text += "-" * 45 + "\n"
summary_text += f"{'T1 sites':>20s}  {'260':>10s}  {'1,289':>10s}\n"
summary_text += f"{'T2 sites':>20s}  {'64':>10s}  {'376':>10s}\n"
summary_text += f"{'Geo shift (pp)':>20s}  {'4.7':>10s}  {'65':>10s}\n"
summary_text += f"{'Lost n':>20s}  {aag_lost_row['n_genes']:>10.0f}  {gccggc_lost['n_genes']:>10.0f}\n"
summary_text += f"{'Lost med LFC':>20s}  {aag_lost_row['median_LFC']:>+10.3f}  {gccggc_lost['median_LFC']:>+10.3f}\n"
summary_text += f"{'Never med LFC':>20s}  {aag_never_row['median_LFC']:>+10.3f}  {gccggc_never['median_LFC']:>+10.3f}\n"
summary_text += f"{'Lost-Never diff':>20s}  {aag_lost_row['median_LFC']-aag_never_row['median_LFC']:>+10.3f}  {gccggc_lost['median_LFC']-gccggc_never['median_LFC']:>+10.3f}\n"
summary_text += f"{'Lost vs Never p':>20s}  {aag_p:>10.2e}  {gccggc_p:>10.2e}\n"
summary_text += f"{'Effect size r':>20s}  {aag_rb:>+10.4f}  {gccggc_rb:>+10.4f}\n\n"

# Add geographic results
summary_text += "GEOGRAPHIC STRATIFICATION\n" + "-" * 40 + "\n"
for region_name in ['core', 'arm']:
    aag_row = geo_df[(geo_df['region'] == region_name) & (geo_df['comparison'] == 'Lost vs Never')]
    gccggc_row = gccggc_geo[(gccggc_geo['region'] == region_name) & (gccggc_geo['comparison'] == 'Lost vs Never')]

    aag_p_r = aag_row.iloc[0]['p_value'] if len(aag_row) > 0 else np.nan
    gcc_p_r = gccggc_row.iloc[0]['p_value'] if len(gccggc_row) > 0 else np.nan

    summary_text += f"{region_name.upper():>10s} Lost vs Never:\n"
    summary_text += f"  AAGCCCG p={aag_p_r:.2e}\n" if not np.isnan(aag_p_r) else f"  AAGCCCG: N/A\n"
    summary_text += f"  GCCGGC  p={gcc_p_r:.2e}\n" if not np.isnan(gcc_p_r) else f"  GCCGGC: N/A\n"

ax.text(0.05, 0.95, summary_text, transform=ax.transAxes, fontsize=8,
        verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.5))
ax.set_title('C. Summary Statistics', fontweight='bold')

fig.suptitle('H21: AAGCCCG vs GCCGGC Temporal Causality Comparison',
             fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(FIGURES / 'AAGCCCG_vs_GCCGGC_comparison.pdf', dpi=300, bbox_inches='tight')
plt.savefig(FIGURES / 'AAGCCCG_vs_GCCGGC_comparison.svg', bbox_inches='tight')
plt.close()
print("Saved AAGCCCG_vs_GCCGGC_comparison.pdf/svg")

# ── Figure 4: Methylation frequency vs LFC scatter ──
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Panel A: All T1-proximal genes
ax = axes[0]
if len(all_t1_freq) > 0:
    scatter = ax.scatter(all_t1_freq['mean_freq'], all_t1_freq['log2FoldChange'],
                         alpha=0.4, s=20, c='#9b59b6', edgecolors='none')
    if len(all_t1_freq) >= 5:
        z = np.polyfit(all_t1_freq['mean_freq'], all_t1_freq['log2FoldChange'], 1)
        p_line = np.poly1d(z)
        x_range = np.linspace(all_t1_freq['mean_freq'].min(), all_t1_freq['mean_freq'].max(), 100)
        ax.plot(x_range, p_line(x_range), 'k--', alpha=0.7)
        rho, p_corr = stats.spearmanr(all_t1_freq['mean_freq'], all_t1_freq['log2FoldChange'])
        ax.text(0.05, 0.95, f'rho={rho:.3f}, p={p_corr:.2e}\nn={len(all_t1_freq)}',
                transform=ax.transAxes, fontsize=10, va='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
ax.set_xlabel('Mean T1 AAGCCCG Methylation Frequency (%)')
ax.set_ylabel('T2 vs T1 log2 Fold Change')
ax.set_title('A. All T1-proximal genes', fontweight='bold')

# Panel B: Lost genes only
ax = axes[1]
if len(lost_freq) > 0:
    ax.scatter(lost_freq['mean_freq'], lost_freq['log2FoldChange'],
               alpha=0.5, s=25, c='#e74c3c', edgecolors='none')
    if len(lost_freq) >= 5:
        z = np.polyfit(lost_freq['mean_freq'], lost_freq['log2FoldChange'], 1)
        p_line = np.poly1d(z)
        x_range = np.linspace(lost_freq['mean_freq'].min(), lost_freq['mean_freq'].max(), 100)
        ax.plot(x_range, p_line(x_range), 'k--', alpha=0.7)
        rho, p_corr = stats.spearmanr(lost_freq['mean_freq'], lost_freq['log2FoldChange'])
        ax.text(0.05, 0.95, f'rho={rho:.3f}, p={p_corr:.2e}\nn={len(lost_freq)}',
                transform=ax.transAxes, fontsize=10, va='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
ax.set_xlabel('Mean T1 AAGCCCG Methylation Frequency (%)')
ax.set_ylabel('T2 vs T1 log2 Fold Change')
ax.set_title('B. Lost genes only', fontweight='bold')

fig.suptitle('H21: T1 AAGCCCG Methylation Frequency vs Expression Change',
             fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig(FIGURES / 'frequency_vs_LFC.pdf', dpi=300, bbox_inches='tight')
plt.savefig(FIGURES / 'frequency_vs_LFC.svg', bbox_inches='tight')
plt.close()
print("Saved frequency_vs_LFC.pdf/svg")

# ── Figure 5: Comprehensive 4-panel summary ──
fig = plt.figure(figsize=(16, 14))
gs = GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.3)

# Panel A: Violin plot T2vsT1
ax = fig.add_subplot(gs[0, 0])
data_list = []
lab_list = []
col_list = []
for lfc, name in [(lost_lfc, 'Lost'), (gained_lfc, 'Gained'),
                    (both_lfc, 'Both'), (never_lfc, 'Never')]:
    if len(lfc) > 0:
        data_list.append(lfc.values)
        lab_list.append(f'{name}\n(n={len(lfc)})')
        col_list.append(colors[name])

if len(data_list) > 0:
    parts = ax.violinplot(data_list, positions=range(len(data_list)),
                          showmeans=True, showmedians=True)
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(col_list[i])
        pc.set_alpha(0.7)
    parts['cmeans'].set_color('black')
    parts['cmedians'].set_color('darkred')

    for i, d in enumerate(data_list):
        med = np.median(d)
        ax.text(i, ax.get_ylim()[1] * 0.88, f'med={med:.2f}', ha='center', fontsize=8, color='darkred')

ax.set_xticks(range(len(lab_list)))
ax.set_xticklabels(lab_list, fontsize=9)
ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
ax.set_ylabel('log2 Fold Change (T2/T1)')
ax.set_title('A. Expression by AAGCCCG Transition', fontweight='bold')

# Panel B: Geographic composition
ax = fig.add_subplot(gs[0, 1])
for i, name in enumerate(group_names):
    gc = geo_composition.get(name, {})
    total = gc.get('n_total', 0)
    if total == 0:
        continue
    core_frac = gc['n_core'] / total
    arm_frac = gc['n_arm'] / total
    ax.bar(i, core_frac, 0.6, color='#3498db', alpha=0.8,
           label='Core' if i == 0 else '')
    ax.bar(i, arm_frac, 0.6, bottom=core_frac, color='#e67e22', alpha=0.8,
           label='Arm' if i == 0 else '')
    ax.text(i, 0.5, f'{gc["pct_core"]:.0f}%', ha='center', fontsize=10, color='white', fontweight='bold')
    ax.text(i, 1.02, f'n={total}', ha='center', fontsize=8)

ax.set_xticks(x_pos)
ax.set_xticklabels(group_names, fontsize=10)
ax.set_ylabel('Fraction')
ax.set_title('B. Geographic Composition', fontweight='bold')
ax.legend(loc='upper right', fontsize=8)
ax.set_ylim(0, 1.15)

# Panel C: Frequency vs LFC for Lost genes
ax = fig.add_subplot(gs[1, 0])
if len(lost_freq) > 0:
    ax.scatter(lost_freq['mean_freq'], lost_freq['log2FoldChange'],
               alpha=0.5, s=25, c='#e74c3c', edgecolors='none')
    if len(lost_freq) >= 5:
        z = np.polyfit(lost_freq['mean_freq'], lost_freq['log2FoldChange'], 1)
        p_line = np.poly1d(z)
        x_range = np.linspace(lost_freq['mean_freq'].min(), lost_freq['mean_freq'].max(), 100)
        ax.plot(x_range, p_line(x_range), 'k--', alpha=0.7)
        rho, p_corr = stats.spearmanr(lost_freq['mean_freq'], lost_freq['log2FoldChange'])
        ax.text(0.05, 0.95, f'rho={rho:.3f}, p={p_corr:.2e}',
                transform=ax.transAxes, fontsize=10, va='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
ax.set_xlabel('Mean T1 AAGCCCG Frequency (%)')
ax.set_ylabel('T2 vs T1 LFC')
ax.set_title('C. Methylation Frequency vs De-repression (Lost)', fontweight='bold')

# Panel D: Statistical summary
ax = fig.add_subplot(gs[1, 1])
ax.axis('off')

summary_text = "STATISTICAL SUMMARY\n" + "=" * 42 + "\n\n"
summary_text += "PRIMARY TEST (T2vsT1 LFC):\n"
lvn = stat_df[(stat_df['comparison'] == 'Lost vs Never') & (stat_df['timepoint'] == 'T2vsT1')]
if len(lvn) > 0 and not np.isnan(lvn.iloc[0]['p_value']):
    r = lvn.iloc[0]
    sig = "***" if r['p_bonferroni'] < 0.001 else "**" if r['p_bonferroni'] < 0.01 else "*" if r['p_bonferroni'] < 0.05 else "ns"
    summary_text += f"  Lost vs Never: p={r['p_value']:.2e}\n"
    summary_text += f"  Bonf. p={r['p_bonferroni']:.2e} {sig}\n"
    summary_text += f"  r={r['rank_biserial']:.4f}\n"
    summary_text += f"  median diff={r['median_diff']:.4f}\n\n"

summary_text += "SECONDARY TESTS:\n"
for comp in ['Gained vs Never', 'Both vs Never']:
    row = stat_df[(stat_df['comparison'] == comp) & (stat_df['timepoint'] == 'T2vsT1')]
    if len(row) > 0 and not np.isnan(row.iloc[0]['p_value']):
        r = row.iloc[0]
        sig = "*" if r['p_bonferroni'] < 0.05 else "ns"
        summary_text += f"  {comp}: p={r['p_value']:.2e} {sig}\n"
        summary_text += f"    r={r['rank_biserial']:.4f}\n"

summary_text += "\nGEOGRAPHIC CONTROLS:\n" + "-" * 30 + "\n"
for _, row in geo_df[geo_df['comparison'] == 'Lost vs Never'].iterrows():
    if np.isnan(row['p_value']):
        summary_text += f"  {row['region']:>5s}: SKIPPED (n={row['n1']:.0f} vs {row['n2']:.0f})\n"
    else:
        sig = "*" if row['p_value'] < 0.05 else "ns"
        summary_text += f"  {row['region']:>5s}: p={row['p_value']:.2e} {sig}\n"
        summary_text += f"        r={row['rank_biserial']:.4f}\n"

ax.text(0.05, 0.95, summary_text, transform=ax.transAxes, fontsize=8,
        verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
ax.set_title('D. Test Results', fontweight='bold')

fig.suptitle('H21: AAGCCCG 6mA Temporal Causality Test', fontsize=15, fontweight='bold', y=0.98)
plt.savefig(FIGURES / 'H21_comprehensive_summary.pdf', dpi=300, bbox_inches='tight')
plt.savefig(FIGURES / 'H21_comprehensive_summary.svg', bbox_inches='tight')
plt.close()
print("Saved H21_comprehensive_summary.pdf/svg")

# ══════════════════════════════════════════════════════════════════════════
# FINAL VERDICT
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("FINAL VERDICT")
print("=" * 70)

print(f"\nAAGCCCG T1->T2 transition groups:")
print(f"  Lost:   n={len(lost_expr):>5}, median T2vsT1 LFC = {lost_lfc.median():>+.4f}")
print(f"  Gained: n={len(gained_expr):>5}, median T2vsT1 LFC = {gained_lfc.median():>+.4f}")
if len(both_lfc) > 0:
    print(f"  Both:   n={len(both_expr):>5}, median T2vsT1 LFC = {both_lfc.median():>+.4f}")
print(f"  Never:  n={len(never_expr):>5}, median T2vsT1 LFC = {never_lfc.median():>+.4f}")

# Check criteria
lvn_row = stat_df[(stat_df['comparison'] == 'Lost vs Never') & (stat_df['timepoint'] == 'T2vsT1')]
if len(lvn_row) > 0 and not np.isnan(lvn_row.iloc[0]['p_value']):
    lost_med = lost_lfc.median()
    never_med = never_lfc.median()
    p_bonf = lvn_row.iloc[0]['p_bonferroni']
    rb = lvn_row.iloc[0]['rank_biserial']

    crit_derepression = lost_med > never_med  # Lost higher LFC = de-repression
    crit_significant = p_bonf < 0.05
    crit_direction = rb > 0  # Lost stochastically larger (check sign convention)

    # Check geographic stratification
    geo_core = geo_df[(geo_df['region'] == 'core') & (geo_df['comparison'] == 'Lost vs Never')]
    geo_arm = geo_df[(geo_df['region'] == 'arm') & (geo_df['comparison'] == 'Lost vs Never')]

    core_survives = False
    arm_survives = False
    if len(geo_core) > 0 and not np.isnan(geo_core.iloc[0]['p_value']):
        if geo_core.iloc[0]['p_value'] < 0.05:
            core_survives = True
    if len(geo_arm) > 0 and not np.isnan(geo_arm.iloc[0]['p_value']):
        if geo_arm.iloc[0]['p_value'] < 0.05:
            arm_survives = True

    geo_survives = core_survives or arm_survives

    print(f"\nCriteria check:")
    print(f"  De-repression (Lost > Never): {crit_derepression} (diff={lost_med - never_med:+.4f})")
    print(f"  Significant (Bonf p<0.05):    {crit_significant} (p_bonf={p_bonf:.2e})")
    print(f"  Geographic stratification:    core={'survives' if core_survives else 'fails'}, arm={'survives' if arm_survives else 'fails'}")

    if crit_derepression and crit_significant and geo_survives:
        verdict = "SUPPORTED"
        core_p_str = "p={:.2e}".format(geo_core.iloc[0]['p_value']) if core_survives else "NS"
        arm_p_str = "p={:.2e}".format(geo_arm.iloc[0]['p_value']) if arm_survives else "NS"
        evidence = (f"AAGCCCG loss -> de-repression is significant (Lost vs Never median LFC diff={lost_med - never_med:+.4f}, "
                    f"p_bonf={p_bonf:.2e}, r={rb:+.4f}) and survives geographic stratification "
                    f"(core: {core_p_str}, arm: {arm_p_str}). "
                    f"Unlike GCCGGC (H19), the effect is not driven by geographic confounding.")
    elif crit_derepression and crit_significant and not geo_survives:
        verdict = "PARTIAL"
        core_p_str = "{:.2e}".format(geo_core.iloc[0]['p_value']) if (len(geo_core) > 0 and not np.isnan(geo_core.iloc[0]['p_value'])) else "NA"
        arm_p_str = "{:.2e}".format(geo_arm.iloc[0]['p_value']) if (len(geo_arm) > 0 and not np.isnan(geo_arm.iloc[0]['p_value'])) else "NA"
        evidence = (f"AAGCCCG loss shows significant de-repression overall (p_bonf={p_bonf:.2e}, r={rb:+.4f}), "
                    f"but does not survive geographic stratification "
                    f"(core p={core_p_str}, arm p={arm_p_str}). "
                    f"May indicate residual geographic confounding or insufficient power within strata.")
    elif crit_derepression and not crit_significant:
        verdict = "PARTIAL"
        evidence = (f"AAGCCCG loss shows directionally correct de-repression (Lost median={lost_med:+.4f} vs "
                    f"Never median={never_med:+.4f}, diff={lost_med-never_med:+.4f}) but not significant after "
                    f"Bonferroni correction (p_bonf={p_bonf:.2e}). Limited power due to small sample size "
                    f"({len(lost_lfc)} Lost genes) may explain non-significance.")
    elif not crit_derepression:
        verdict = "REJECTED"
        evidence = (f"AAGCCCG loss does NOT show de-repression. Lost median LFC={lost_med:+.4f} vs "
                    f"Never median LFC={never_med:+.4f} (diff={lost_med-never_med:+.4f}).")
    else:
        verdict = "REJECTED"
        evidence = "No evidence for de-repression effect."
else:
    verdict = "REJECTED"
    evidence = "Could not perform primary statistical test (insufficient data)."

print(f"\n{'='*70}")
print(f"VERDICT: {verdict}")
print(f"EVIDENCE: {evidence}")
print(f"{'='*70}")

# Save verdict
with open(TABLES / "verdict.txt", "w") as f:
    f.write(f"VERDICT: {verdict}\n")
    f.write(f"EVIDENCE: {evidence}\n")
print(f"\nSaved verdict.txt")

print("\nDone.")
