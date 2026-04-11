#!/usr/bin/env python3
"""
H19: Temporal De-repression Analysis
Tests whether genes losing GCCGGC methylation at T2 show upregulation,
and genes gaining methylation show suppression.
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
ANALYSIS = BASE / "11_epigenome_integration/analysis/42_GCCGGC_temporal_derepression"
FIGURES = ANALYSIS / "figures"
TABLES = ANALYSIS / "tables"

GCCGGC_SITES = BASE / "11_epigenome_integration/analysis/37_defense_island_GCCGGC/tables/GCCGGC_sites_by_timepoint.tsv"
AAGCCCG_MAPPING = BASE / "11_epigenome_integration/analysis/36_AAGCCCG_distribution/tables/AAGCCCG_site_gene_mapping.tsv"
GFF = BASE / "11_epigenome_integration/analysis/39_GCCGGC_MTase_reverse_ID/data/GCF_000203835.1_ASM20383v1_genomic.gff"

DESEQ_T2T1 = BASE / "04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_2_vs_1.tsv"
DESEQ_T3T1 = BASE / "04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_3_vs_1.tsv"
DESEQ_T3T2 = BASE / "04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_3_vs_2.tsv"

GENOME_SIZE = 8667507  # S. coelicolor A3(2) chromosome

# ── Load GFF genes ─────────────────────────────────────────────────────────
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
    """Classify position as core or arm based on S. coelicolor coordinates."""
    # Arms: 0-1.5Mb and 6.4Mb-8.67Mb; Core: 1.5Mb-6.4Mb
    if pos < 1_500_000 or pos > 6_400_000:
        return 'arm'
    return 'core'


def map_sites_to_genes(sites_df, genes_df, distance=2000):
    """Map methylation sites to genes within specified distance."""
    site_gene_pairs = []

    for _, site in sites_df.iterrows():
        site_pos = site['position']
        site_chrom = site['chrom'] if 'chrom' in site.index else 'NC_003888.3'

        for _, gene in genes_df.iterrows():
            if gene['chrom'] != site_chrom:
                continue

            # Calculate distance from site to gene
            if site_pos < gene['start']:
                dist = gene['start'] - site_pos
            elif site_pos > gene['end']:
                dist = site_pos - gene['end']
            else:
                dist = 0  # Inside gene

            if dist <= distance:
                site_gene_pairs.append({
                    'site_position': site_pos,
                    'locus_tag': gene['locus_tag'],
                    'distance': dist,
                    'timepoint': site['timepoint'],
                    'frequency': site.get('frequency', np.nan)
                })

    return pd.DataFrame(site_gene_pairs)


def rank_biserial(x, y):
    """Compute rank-biserial correlation (effect size for Wilcoxon)."""
    n1, n2 = len(x), len(y)
    U, _ = stats.mannwhitneyu(x, y, alternative='two-sided')
    return 1 - (2 * U) / (n1 * n2)


# ══════════════════════════════════════════════════════════════════════════
# STEP 1: Load and classify GCCGGC sites
# ══════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("H19: GCCGGC TEMPORAL DE-REPRESSION ANALYSIS")
print("=" * 70)

# Load GCCGGC sites
sites = pd.read_csv(GCCGGC_SITES, sep='\t')
sites = sites[sites['timepoint'].isin(['T1', 'T2', 'T3'])]

t1_sites = sites[sites['timepoint'] == 'T1'].copy()
t2_sites = sites[sites['timepoint'] == 'T2'].copy()
t3_sites = sites[sites['timepoint'] == 'T3'].copy()

print(f"\nGCCGGC sites: T1={len(t1_sites)}, T2={len(t2_sites)}, T3={len(t3_sites)}")

# Load genes
genes_df = load_genes_from_gff(GFF)
genes_df['region'] = genes_df.apply(
    lambda r: classify_region((r['start'] + r['end']) / 2), axis=1
)
print(f"Total genes in GFF: {len(genes_df)}")
print(f"  Core: {(genes_df['region'] == 'core').sum()}, Arm: {(genes_df['region'] == 'arm').sum()}")

# ══════════════════════════════════════════════════════════════════════════
# STEP 2: Map sites to genes (within 2kb)
# ══════════════════════════════════════════════════════════════════════════
print("\nMapping GCCGGC sites to genes within 2kb...")

t1_mapping = map_sites_to_genes(t1_sites, genes_df, distance=2000)
t2_mapping = map_sites_to_genes(t2_sites, genes_df, distance=2000)
t3_mapping = map_sites_to_genes(t3_sites, genes_df, distance=2000)

t1_genes = set(t1_mapping['locus_tag'].unique())
t2_genes = set(t2_mapping['locus_tag'].unique())
t3_genes = set(t3_mapping['locus_tag'].unique())
all_genes = set(genes_df['locus_tag'])

print(f"Genes near T1 GCCGGC sites: {len(t1_genes)}")
print(f"Genes near T2 GCCGGC sites: {len(t2_genes)}")
print(f"Genes near T3 GCCGGC sites: {len(t3_genes)}")

# ══════════════════════════════════════════════════════════════════════════
# STEP 3: Classify T1→T2 transition groups
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("T1→T2 METHYLATION TRANSITION GROUPS")
print("=" * 70)

lost_genes = t1_genes - t2_genes    # Had T1, no T2
gained_genes = t2_genes - t1_genes  # No T1, has T2
both_genes = t1_genes & t2_genes    # Had both
never_genes = all_genes - t1_genes - t2_genes  # Neither

print(f"Lost (T1→no T2):     {len(lost_genes)} genes")
print(f"Gained (no T1→T2):   {len(gained_genes)} genes")
print(f"Both (T1 & T2):      {len(both_genes)} genes")
print(f"Never (no T1, no T2):{len(never_genes)} genes")

# ══════════════════════════════════════════════════════════════════════════
# STEP 4: Load DESeq2 results
# ══════════════════════════════════════════════════════════════════════════
deseq_t2t1 = pd.read_csv(DESEQ_T2T1, sep='\t')
deseq_t3t1 = pd.read_csv(DESEQ_T3T1, sep='\t')
deseq_t3t2 = pd.read_csv(DESEQ_T3T2, sep='\t')

# DESeq2 gene set
deseq_genes = set(deseq_t2t1['gene_id'])
print(f"\nDESeq2 genes: {len(deseq_genes)}")

# Filter transition groups to genes with DESeq2 data
lost_expr = lost_genes & deseq_genes
gained_expr = gained_genes & deseq_genes
both_expr = both_genes & deseq_genes
never_expr = never_genes & deseq_genes

print(f"\nWith expression data:")
print(f"  Lost:   {len(lost_expr)}")
print(f"  Gained: {len(gained_expr)}")
print(f"  Both:   {len(both_expr)}")
print(f"  Never:  {len(never_expr)}")

# Extract LFC for each group
def get_lfc(gene_set, deseq_df):
    return deseq_df[deseq_df['gene_id'].isin(gene_set)]['log2FoldChange'].dropna()

lost_lfc_t2t1 = get_lfc(lost_expr, deseq_t2t1)
gained_lfc_t2t1 = get_lfc(gained_expr, deseq_t2t1)
both_lfc_t2t1 = get_lfc(both_expr, deseq_t2t1)
never_lfc_t2t1 = get_lfc(never_expr, deseq_t2t1)

lost_lfc_t3t1 = get_lfc(lost_expr, deseq_t3t1)
gained_lfc_t3t1 = get_lfc(gained_expr, deseq_t3t1)
both_lfc_t3t1 = get_lfc(both_expr, deseq_t3t1)
never_lfc_t3t1 = get_lfc(never_expr, deseq_t3t1)

# ══════════════════════════════════════════════════════════════════════════
# STEP 5: Expression summary by transition group
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("EXPRESSION CHANGES BY TRANSITION GROUP (T2 vs T1)")
print("=" * 70)

def expression_summary(name, lfc, deseq_df, gene_set):
    """Compute expression summary for a transition group."""
    sub = deseq_df[deseq_df['gene_id'].isin(gene_set)].copy()
    n = len(sub)
    med = lfc.median()
    mean = lfc.mean()
    se = lfc.std() / np.sqrt(n) if n > 0 else np.nan
    ci_lo = mean - 1.96 * se
    ci_hi = mean + 1.96 * se

    sig = sub[sub['padj'] < 0.05]
    n_deg = len(sig)
    n_up = (sig['log2FoldChange'] > 0).sum()
    n_down = (sig['log2FoldChange'] < 0).sum()
    frac_deg = n_deg / n if n > 0 else 0
    frac_up = n_up / n if n > 0 else 0
    frac_down = n_down / n if n > 0 else 0

    print(f"\n{name} (n={n}):")
    print(f"  Median LFC: {med:.4f}")
    print(f"  Mean LFC:   {mean:.4f} (95% CI: [{ci_lo:.4f}, {ci_hi:.4f}])")
    print(f"  DEGs:       {n_deg} ({frac_deg:.1%}) | Up: {n_up} ({frac_up:.1%}) | Down: {n_down} ({frac_down:.1%})")

    return {
        'group': name, 'n_genes': n,
        'median_LFC': med, 'mean_LFC': mean,
        'CI_lower': ci_lo, 'CI_upper': ci_hi,
        'n_DEG': n_deg, 'frac_DEG': frac_deg,
        'n_up': n_up, 'frac_up': frac_up,
        'n_down': n_down, 'frac_down': frac_down
    }

summary_rows = []
summary_rows.append(expression_summary("Lost", lost_lfc_t2t1, deseq_t2t1, lost_expr))
summary_rows.append(expression_summary("Gained", gained_lfc_t2t1, deseq_t2t1, gained_expr))
summary_rows.append(expression_summary("Both", both_lfc_t2t1, deseq_t2t1, both_expr))
summary_rows.append(expression_summary("Never", never_lfc_t2t1, deseq_t2t1, never_expr))

# T3vsT1
print("\n" + "-" * 50)
print("T3 vs T1 LFC by T1→T2 transition group:")
for name, lfc in [("Lost", lost_lfc_t3t1), ("Gained", gained_lfc_t3t1),
                   ("Both", both_lfc_t3t1), ("Never", never_lfc_t3t1)]:
    print(f"  {name}: median={lfc.median():.4f}, mean={lfc.mean():.4f}")

summary_df = pd.DataFrame(summary_rows)

# ══════════════════════════════════════════════════════════════════════════
# STEP 6: Statistical tests
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("STATISTICAL TESTS (Wilcoxon rank-sum on T2vsT1 LFC)")
print("=" * 70)

test_pairs = [
    ("Lost vs Never", lost_lfc_t2t1, never_lfc_t2t1),
    ("Lost vs Gained", lost_lfc_t2t1, gained_lfc_t2t1),
    ("Gained vs Never", gained_lfc_t2t1, never_lfc_t2t1),
    ("Both vs Never", both_lfc_t2t1, never_lfc_t2t1),
    ("Lost vs Both", lost_lfc_t2t1, both_lfc_t2t1),
    ("Gained vs Both", gained_lfc_t2t1, both_lfc_t2t1),
]

stat_rows = []
for label, x, y in test_pairs:
    if len(x) < 3 or len(y) < 3:
        print(f"\n{label}: SKIPPED (n too small: {len(x)}, {len(y)})")
        stat_rows.append({
            'comparison': label, 'n1': len(x), 'n2': len(y),
            'U': np.nan, 'p_value': np.nan, 'p_bonferroni': np.nan,
            'rank_biserial': np.nan,
            'median_diff': x.median() - y.median() if len(x) > 0 and len(y) > 0 else np.nan,
            'direction': 'NA'
        })
        continue

    U, p = stats.mannwhitneyu(x, y, alternative='two-sided')
    p_bonf = min(p * 6, 1.0)
    rb = rank_biserial(x, y)
    med_diff = x.median() - y.median()
    direction = "Group1 higher" if med_diff > 0 else "Group1 lower"

    print(f"\n{label}:")
    print(f"  n = {len(x)} vs {len(y)}")
    print(f"  U = {U:.0f}, p = {p:.2e}, Bonferroni p = {p_bonf:.2e}")
    print(f"  Rank-biserial r = {rb:.4f}")
    print(f"  Median diff = {med_diff:.4f} ({direction})")

    stat_rows.append({
        'comparison': label, 'n1': len(x), 'n2': len(y),
        'U': U, 'p_value': p, 'p_bonferroni': p_bonf,
        'rank_biserial': rb, 'median_diff': med_diff,
        'direction': direction
    })

stat_df = pd.DataFrame(stat_rows)

# ══════════════════════════════════════════════════════════════════════════
# STEP 7: Confound controls
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("CONFOUND CONTROLS")
print("=" * 70)

# 7a. Geographic stratification
print("\n--- Geographic stratification ---")
gene_regions = genes_df[['locus_tag', 'region']].drop_duplicates()

geo_rows = []
for region_name in ['core', 'arm']:
    region_genes = set(gene_regions[gene_regions['region'] == region_name]['locus_tag'])

    lost_r = lost_expr & region_genes
    never_r = never_expr & region_genes
    gained_r = gained_expr & region_genes

    lost_lfc_r = get_lfc(lost_r, deseq_t2t1)
    never_lfc_r = get_lfc(never_r, deseq_t2t1)
    gained_lfc_r = get_lfc(gained_r, deseq_t2t1)

    print(f"\n{region_name.upper()} genes:")
    print(f"  Lost: n={len(lost_lfc_r)}, median LFC={lost_lfc_r.median():.4f}" if len(lost_lfc_r) > 0 else f"  Lost: n=0")
    print(f"  Never: n={len(never_lfc_r)}, median LFC={never_lfc_r.median():.4f}" if len(never_lfc_r) > 0 else f"  Never: n=0")
    print(f"  Gained: n={len(gained_lfc_r)}, median LFC={gained_lfc_r.median():.4f}" if len(gained_lfc_r) > 0 else f"  Gained: n=0")

    # Lost vs Never within region
    if len(lost_lfc_r) >= 3 and len(never_lfc_r) >= 3:
        U, p = stats.mannwhitneyu(lost_lfc_r, never_lfc_r, alternative='two-sided')
        rb = rank_biserial(lost_lfc_r, never_lfc_r)
        print(f"  Lost vs Never: U={U:.0f}, p={p:.2e}, r={rb:.4f}")
        geo_rows.append({
            'region': region_name, 'comparison': 'Lost vs Never',
            'n1': len(lost_lfc_r), 'n2': len(never_lfc_r),
            'median1': lost_lfc_r.median(), 'median2': never_lfc_r.median(),
            'U': U, 'p_value': p, 'rank_biserial': rb
        })
    else:
        print(f"  Lost vs Never: SKIPPED (n too small)")
        geo_rows.append({
            'region': region_name, 'comparison': 'Lost vs Never',
            'n1': len(lost_lfc_r), 'n2': len(never_lfc_r),
            'median1': lost_lfc_r.median() if len(lost_lfc_r) > 0 else np.nan,
            'median2': never_lfc_r.median() if len(never_lfc_r) > 0 else np.nan,
            'U': np.nan, 'p_value': np.nan, 'rank_biserial': np.nan
        })

    # Gained vs Never within region
    if len(gained_lfc_r) >= 3 and len(never_lfc_r) >= 3:
        U, p = stats.mannwhitneyu(gained_lfc_r, never_lfc_r, alternative='two-sided')
        rb = rank_biserial(gained_lfc_r, never_lfc_r)
        print(f"  Gained vs Never: U={U:.0f}, p={p:.2e}, r={rb:.4f}")
        geo_rows.append({
            'region': region_name, 'comparison': 'Gained vs Never',
            'n1': len(gained_lfc_r), 'n2': len(never_lfc_r),
            'median1': gained_lfc_r.median(), 'median2': never_lfc_r.median(),
            'U': U, 'p_value': p, 'rank_biserial': rb
        })
    else:
        print(f"  Gained vs Never: SKIPPED (n too small)")
        geo_rows.append({
            'region': region_name, 'comparison': 'Gained vs Never',
            'n1': len(gained_lfc_r), 'n2': len(never_lfc_r),
            'median1': gained_lfc_r.median() if len(gained_lfc_r) > 0 else np.nan,
            'median2': never_lfc_r.median() if len(never_lfc_r) > 0 else np.nan,
            'U': np.nan, 'p_value': np.nan, 'rank_biserial': np.nan
        })

geo_df = pd.DataFrame(geo_rows)

# 7b. Baseline expression quintile control
print("\n--- Baseline expression quintile control ---")
deseq_t2t1_full = deseq_t2t1.copy()
deseq_t2t1_full['baseMean_quintile'] = pd.qcut(deseq_t2t1_full['baseMean'], 5, labels=[1,2,3,4,5])

quintile_rows = []
for q in [1, 2, 3, 4, 5]:
    q_genes = set(deseq_t2t1_full[deseq_t2t1_full['baseMean_quintile'] == q]['gene_id'])
    lost_q = lost_expr & q_genes
    never_q = never_expr & q_genes

    lost_lfc_q = get_lfc(lost_q, deseq_t2t1)
    never_lfc_q = get_lfc(never_q, deseq_t2t1)

    row = {'quintile': q, 'n_lost': len(lost_lfc_q), 'n_never': len(never_lfc_q)}

    if len(lost_lfc_q) >= 3 and len(never_lfc_q) >= 3:
        U, p = stats.mannwhitneyu(lost_lfc_q, never_lfc_q, alternative='two-sided')
        rb = rank_biserial(lost_lfc_q, never_lfc_q)
        row.update({
            'lost_median': lost_lfc_q.median(), 'never_median': never_lfc_q.median(),
            'U': U, 'p_value': p, 'rank_biserial': rb
        })
        print(f"  Q{q}: Lost n={len(lost_lfc_q)} (med={lost_lfc_q.median():.3f}) vs Never n={len(never_lfc_q)} (med={never_lfc_q.median():.3f}), p={p:.2e}, r={rb:.3f}")
    else:
        row.update({
            'lost_median': lost_lfc_q.median() if len(lost_lfc_q) > 0 else np.nan,
            'never_median': never_lfc_q.median() if len(never_lfc_q) > 0 else np.nan,
            'U': np.nan, 'p_value': np.nan, 'rank_biserial': np.nan
        })
        print(f"  Q{q}: Lost n={len(lost_lfc_q)} vs Never n={len(never_lfc_q)} -- SKIPPED")

    quintile_rows.append(row)

quintile_df = pd.DataFrame(quintile_rows)

# ══════════════════════════════════════════════════════════════════════════
# STEP 8: T2→T3 transition analysis
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("T2→T3 METHYLATION TRANSITION ANALYSIS")
print("=" * 70)

lost_t2t3 = t2_genes - t3_genes
gained_t2t3 = t3_genes - t2_genes
both_t2t3 = t2_genes & t3_genes
never_t2t3 = all_genes - t2_genes - t3_genes

print(f"Lost (T2→no T3):     {len(lost_t2t3)} genes")
print(f"Gained (no T2→T3):   {len(gained_t2t3)} genes")
print(f"Both (T2 & T3):      {len(both_t2t3)} genes")
print(f"Never (no T2, no T3):{len(never_t2t3)} genes")

# Expression changes T3vsT2
lost_t2t3_expr = lost_t2t3 & deseq_genes
gained_t2t3_expr = gained_t2t3 & deseq_genes
never_t2t3_expr = never_t2t3 & deseq_genes

lost_lfc_t3t2 = get_lfc(lost_t2t3_expr, deseq_t3t2)
gained_lfc_t3t2 = get_lfc(gained_t2t3_expr, deseq_t3t2)
never_lfc_t3t2 = get_lfc(never_t2t3_expr, deseq_t3t2)

print(f"\nT3vsT2 LFC by T2→T3 transition:")
for name, lfc in [("Lost", lost_lfc_t3t2), ("Gained", gained_lfc_t3t2), ("Never", never_lfc_t3t2)]:
    if len(lfc) > 0:
        print(f"  {name}: n={len(lfc)}, median={lfc.median():.4f}, mean={lfc.mean():.4f}")

if len(lost_lfc_t3t2) >= 3 and len(never_lfc_t3t2) >= 3:
    U, p = stats.mannwhitneyu(lost_lfc_t3t2, never_lfc_t3t2, alternative='two-sided')
    rb = rank_biserial(lost_lfc_t3t2, never_lfc_t3t2)
    print(f"\n  Lost vs Never (T3vsT2): U={U:.0f}, p={p:.2e}, r={rb:.4f}")

# ══════════════════════════════════════════════════════════════════════════
# STEP 9: Reciprocal AAGCCCG analysis
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("RECIPROCAL AAGCCCG ANALYSIS")
print("=" * 70)

aag_mapping = pd.read_csv(AAGCCCG_MAPPING, sep='\t')
aag_t1 = aag_mapping[aag_mapping['timepoint'] == 'T1']
aag_t2 = aag_mapping[aag_mapping['timepoint'] == 'T2']

# Keep only genes within 2kb
aag_t1_2kb = aag_t1[aag_t1['within_2kb'] == True]
aag_t2_2kb = aag_t2[aag_t2['within_2kb'] == True]

aag_t1_genes = set(aag_t1_2kb['locus_tag'].unique())
aag_t2_genes = set(aag_t2_2kb['locus_tag'].unique())

aag_lost = aag_t1_genes - aag_t2_genes
aag_gained = aag_t2_genes - aag_t1_genes
aag_both = aag_t1_genes & aag_t2_genes
aag_never = deseq_genes - aag_t1_genes - aag_t2_genes

print(f"AAGCCCG sites: T1={len(aag_t1_2kb['position'].unique())}, T2={len(aag_t2_2kb['position'].unique())}")
print(f"Lost (T1→no T2):     {len(aag_lost)} genes")
print(f"Gained (no T1→T2):   {len(aag_gained)} genes")
print(f"Both (T1 & T2):      {len(aag_both)} genes")
print(f"Never:               {len(aag_never)} genes")

aag_lost_expr = aag_lost & deseq_genes
aag_gained_expr = aag_gained & deseq_genes
aag_never_expr = aag_never & deseq_genes

aag_lost_lfc = get_lfc(aag_lost_expr, deseq_t2t1)
aag_gained_lfc = get_lfc(aag_gained_expr, deseq_t2t1)
aag_never_lfc = get_lfc(aag_never_expr, deseq_t2t1)

print(f"\nAAGCCCG T2vsT1 LFC by transition:")
for name, lfc in [("Lost", aag_lost_lfc), ("Gained", aag_gained_lfc), ("Never", aag_never_lfc)]:
    if len(lfc) > 0:
        print(f"  {name}: n={len(lfc)}, median={lfc.median():.4f}, mean={lfc.mean():.4f}")

aag_stat_rows = []
if len(aag_lost_lfc) >= 3 and len(aag_never_lfc) >= 3:
    U, p = stats.mannwhitneyu(aag_lost_lfc, aag_never_lfc, alternative='two-sided')
    rb = rank_biserial(aag_lost_lfc, aag_never_lfc)
    print(f"\n  AAGCCCG Lost vs Never: U={U:.0f}, p={p:.2e}, r={rb:.4f}")
    aag_stat_rows.append({'comparison': 'AAGCCCG Lost vs Never', 'n1': len(aag_lost_lfc), 'n2': len(aag_never_lfc),
                          'U': U, 'p_value': p, 'rank_biserial': rb,
                          'median_diff': aag_lost_lfc.median() - aag_never_lfc.median()})

if len(aag_gained_lfc) >= 3 and len(aag_never_lfc) >= 3:
    U, p = stats.mannwhitneyu(aag_gained_lfc, aag_never_lfc, alternative='two-sided')
    rb = rank_biserial(aag_gained_lfc, aag_never_lfc)
    print(f"  AAGCCCG Gained vs Never: U={U:.0f}, p={p:.2e}, r={rb:.4f}")
    aag_stat_rows.append({'comparison': 'AAGCCCG Gained vs Never', 'n1': len(aag_gained_lfc), 'n2': len(aag_never_lfc),
                          'U': U, 'p_value': p, 'rank_biserial': rb,
                          'median_diff': aag_gained_lfc.median() - aag_never_lfc.median()})

if len(aag_lost_lfc) >= 3 and len(aag_gained_lfc) >= 3:
    U, p = stats.mannwhitneyu(aag_lost_lfc, aag_gained_lfc, alternative='two-sided')
    rb = rank_biserial(aag_lost_lfc, aag_gained_lfc)
    print(f"  AAGCCCG Lost vs Gained: U={U:.0f}, p={p:.2e}, r={rb:.4f}")
    aag_stat_rows.append({'comparison': 'AAGCCCG Lost vs Gained', 'n1': len(aag_lost_lfc), 'n2': len(aag_gained_lfc),
                          'U': U, 'p_value': p, 'rank_biserial': rb,
                          'median_diff': aag_lost_lfc.median() - aag_gained_lfc.median()})

# ══════════════════════════════════════════════════════════════════════════
# STEP 10: Save tables
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("SAVING TABLES")
print("=" * 70)

# Gene methylation transitions table
transition_records = []
for gene in sorted(all_genes):
    if gene not in deseq_genes:
        continue

    has_t1 = gene in t1_genes
    has_t2 = gene in t2_genes
    has_t3 = gene in t3_genes

    if has_t1 and not has_t2:
        group_t1t2 = 'Lost'
    elif not has_t1 and has_t2:
        group_t1t2 = 'Gained'
    elif has_t1 and has_t2:
        group_t1t2 = 'Both'
    else:
        group_t1t2 = 'Never'

    if has_t2 and not has_t3:
        group_t2t3 = 'Lost'
    elif not has_t2 and has_t3:
        group_t2t3 = 'Gained'
    elif has_t2 and has_t3:
        group_t2t3 = 'Both'
    else:
        group_t2t3 = 'Never'

    gene_info = genes_df[genes_df['locus_tag'] == gene].iloc[0] if gene in genes_df['locus_tag'].values else None
    region = gene_info['region'] if gene_info is not None else 'unknown'

    row_t2t1 = deseq_t2t1[deseq_t2t1['gene_id'] == gene]
    row_t3t1 = deseq_t3t1[deseq_t3t1['gene_id'] == gene]
    row_t3t2 = deseq_t3t2[deseq_t3t2['gene_id'] == gene]

    lfc_t2t1_val = row_t2t1['log2FoldChange'].values[0] if len(row_t2t1) > 0 else np.nan
    padj_t2t1_val = row_t2t1['padj'].values[0] if len(row_t2t1) > 0 else np.nan
    lfc_t3t1_val = row_t3t1['log2FoldChange'].values[0] if len(row_t3t1) > 0 else np.nan
    padj_t3t1_val = row_t3t1['padj'].values[0] if len(row_t3t1) > 0 else np.nan
    lfc_t3t2_val = row_t3t2['log2FoldChange'].values[0] if len(row_t3t2) > 0 else np.nan
    padj_t3t2_val = row_t3t2['padj'].values[0] if len(row_t3t2) > 0 else np.nan
    basemean = row_t2t1['baseMean'].values[0] if len(row_t2t1) > 0 else np.nan

    transition_records.append({
        'gene_id': gene,
        'region': region,
        'GCCGGC_T1': has_t1,
        'GCCGGC_T2': has_t2,
        'GCCGGC_T3': has_t3,
        'transition_T1T2': group_t1t2,
        'transition_T2T3': group_t2t3,
        'baseMean': basemean,
        'LFC_T2vsT1': lfc_t2t1_val,
        'padj_T2vsT1': padj_t2t1_val,
        'LFC_T3vsT1': lfc_t3t1_val,
        'padj_T3vsT1': padj_t3t1_val,
        'LFC_T3vsT2': lfc_t3t2_val,
        'padj_T3vsT2': padj_t3t2_val,
    })

transitions_df = pd.DataFrame(transition_records)
transitions_df.to_csv(TABLES / "gene_methylation_transitions.tsv", sep='\t', index=False)
print(f"Saved gene_methylation_transitions.tsv ({len(transitions_df)} genes)")

# Summary table
summary_df.to_csv(TABLES / "transition_group_expression.tsv", sep='\t', index=False)
print("Saved transition_group_expression.tsv")

# Statistical tests
stat_df.to_csv(TABLES / "statistical_tests.tsv", sep='\t', index=False)
print("Saved statistical_tests.tsv")

# Geographic stratified tests
geo_df.to_csv(TABLES / "geographic_stratified_tests.tsv", sep='\t', index=False)
print("Saved geographic_stratified_tests.tsv")

# Quintile tests
quintile_df.to_csv(TABLES / "basemean_quintile_tests.tsv", sep='\t', index=False)
print("Saved basemean_quintile_tests.tsv")

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

# ── Figure 1: Violin plot of LFC by transition group ──
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for ax_idx, (comparison, title) in enumerate([
    ([(lost_lfc_t2t1, 'Lost'), (gained_lfc_t2t1, 'Gained'),
      (both_lfc_t2t1, 'Both'), (never_lfc_t2t1, 'Never')], 'T2 vs T1 LFC'),
    ([(lost_lfc_t3t1, 'Lost'), (gained_lfc_t3t1, 'Gained'),
      (both_lfc_t3t1, 'Both'), (never_lfc_t3t1, 'Never')], 'T3 vs T1 LFC'),
]):
    ax = axes[ax_idx]
    data_list = []
    labels = []
    cols = []

    for lfc, name in comparison:
        if len(lfc) > 0:
            data_list.append(lfc.values)
            labels.append(f"{name}\n(n={len(lfc)})")
            cols.append(colors[name])

    parts = ax.violinplot(data_list, positions=range(len(data_list)),
                          showmeans=True, showmedians=True)

    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(cols[i])
        pc.set_alpha(0.7)

    parts['cmeans'].set_color('black')
    parts['cmedians'].set_color('red')

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=9)
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.set_ylabel('log2 Fold Change')
    ax.set_title(title)

fig.suptitle('H19: GCCGGC Methylation Transition Groups vs Expression Change',
             fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig(FIGURES / 'transition_groups_violin.pdf', dpi=300, bbox_inches='tight')
plt.savefig(FIGURES / 'transition_groups_violin.svg', bbox_inches='tight')
plt.close()
print("Saved transition_groups_violin.pdf/svg")

# ── Figure 2: Fraction up/down bar chart ──
fig, ax = plt.subplots(figsize=(10, 6))

group_names = ['Lost', 'Gained', 'Both', 'Never']
group_data = {
    'Lost': (lost_expr, deseq_t2t1),
    'Gained': (gained_expr, deseq_t2t1),
    'Both': (both_expr, deseq_t2t1),
    'Never': (never_expr, deseq_t2t1),
}

x_pos = np.arange(len(group_names))
width = 0.6

for i, name in enumerate(group_names):
    gene_set, deseq = group_data[name]
    sub = deseq[deseq['gene_id'].isin(gene_set)]
    n = len(sub)
    if n == 0:
        continue
    sig = sub[sub['padj'] < 0.05]
    n_up = (sig['log2FoldChange'] > 0).sum()
    n_down = (sig['log2FoldChange'] < 0).sum()
    n_ns = n - n_up - n_down

    ax.bar(i, n_up / n, width, color='#e74c3c', alpha=0.8,
           label='Up (padj<0.05)' if i == 0 else '')
    ax.bar(i, n_ns / n, width, bottom=n_up / n, color='#bdc3c7', alpha=0.8,
           label='Not significant' if i == 0 else '')
    ax.bar(i, n_down / n, width, bottom=(n_up + n_ns) / n, color='#3498db', alpha=0.8,
           label='Down (padj<0.05)' if i == 0 else '')

    ax.text(i, 1.02, f'n={n}', ha='center', fontsize=9)

ax.set_xticks(x_pos)
ax.set_xticklabels(group_names, fontsize=11)
ax.set_ylabel('Fraction of genes')
ax.set_title('Fraction Up/Down by GCCGGC Methylation Transition (T2 vs T1)', fontweight='bold')
ax.legend(loc='upper right')
ax.set_ylim(0, 1.1)
plt.tight_layout()
plt.savefig(FIGURES / 'fraction_up_down.pdf', dpi=300, bbox_inches='tight')
plt.savefig(FIGURES / 'fraction_up_down.svg', bbox_inches='tight')
plt.close()
print("Saved fraction_up_down.pdf/svg")

# ── Figure 3: Gene-level scatter (T1 methylation count vs LFC) ──
fig, ax = plt.subplots(figsize=(8, 6))

# Count T1 sites per gene
t1_site_counts = t1_mapping.groupby('locus_tag')['site_position'].nunique().reset_index()
t1_site_counts.columns = ['gene_id', 'n_t1_sites']
t1_scatter = t1_site_counts.merge(deseq_t2t1[['gene_id', 'log2FoldChange', 'padj']], on='gene_id')

ax.scatter(t1_scatter['n_t1_sites'], t1_scatter['log2FoldChange'],
           alpha=0.4, s=20, c='#e74c3c', edgecolors='none')

# Trend line
if len(t1_scatter) > 5:
    from numpy.polynomial import polynomial as P
    z = np.polyfit(t1_scatter['n_t1_sites'], t1_scatter['log2FoldChange'], 1)
    p_line = np.poly1d(z)
    x_range = np.linspace(t1_scatter['n_t1_sites'].min(), t1_scatter['n_t1_sites'].max(), 100)
    ax.plot(x_range, p_line(x_range), 'k--', alpha=0.7)

    rho, p_corr = stats.spearmanr(t1_scatter['n_t1_sites'], t1_scatter['log2FoldChange'])
    ax.text(0.05, 0.95, f'Spearman rho={rho:.3f}, p={p_corr:.2e}',
            transform=ax.transAxes, fontsize=10, va='top')

ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
ax.set_xlabel('Number of T1 GCCGGC sites within 2kb')
ax.set_ylabel('T2 vs T1 log2 Fold Change')
ax.set_title('T1 GCCGGC Methylation Density vs Expression Change at T2', fontweight='bold')
plt.tight_layout()
plt.savefig(FIGURES / 'frequency_vs_LFC_scatter.pdf', dpi=300, bbox_inches='tight')
plt.savefig(FIGURES / 'frequency_vs_LFC_scatter.svg', bbox_inches='tight')
plt.close()
print("Saved frequency_vs_LFC_scatter.pdf/svg")

# ── Figure 4: Geographic context ──
fig, axes = plt.subplots(2, 1, figsize=(16, 8), sharex=True)

# Panel A: T1 and T2 GCCGGC site positions
ax = axes[0]
for _, site in t1_sites.iterrows():
    ax.axvline(x=site['position'] / 1e6, color='#e74c3c', alpha=0.15, linewidth=0.5)
for _, site in t2_sites.iterrows():
    ax.axvline(x=site['position'] / 1e6, color='#3498db', alpha=0.3, linewidth=0.8)

ax.axvspan(0, 1.5, color='lightyellow', alpha=0.3, label='Left arm')
ax.axvspan(1.5, 6.4, color='lightblue', alpha=0.1, label='Core')
ax.axvspan(6.4, GENOME_SIZE / 1e6, color='lightyellow', alpha=0.3, label='Right arm')

ax.set_ylabel('GCCGGC Sites')
ax.set_title('GCCGGC Site Distribution by Timepoint', fontweight='bold')
t1_patch = mpatches.Patch(color='#e74c3c', alpha=0.5, label=f'T1 (n={len(t1_sites)})')
t2_patch = mpatches.Patch(color='#3498db', alpha=0.5, label=f'T2 (n={len(t2_sites)})')
ax.legend(handles=[t1_patch, t2_patch], loc='upper right')
ax.set_yticks([])

# Panel B: Transition groups on chromosome
ax = axes[1]
for gene in lost_expr:
    g = genes_df[genes_df['locus_tag'] == gene]
    if len(g) > 0:
        pos = (g.iloc[0]['start'] + g.iloc[0]['end']) / 2 / 1e6
        ax.scatter(pos, 0.7, c='#e74c3c', s=3, alpha=0.3)

for gene in gained_expr:
    g = genes_df[genes_df['locus_tag'] == gene]
    if len(g) > 0:
        pos = (g.iloc[0]['start'] + g.iloc[0]['end']) / 2 / 1e6
        ax.scatter(pos, 0.3, c='#3498db', s=3, alpha=0.3)

ax.axvspan(0, 1.5, color='lightyellow', alpha=0.3)
ax.axvspan(1.5, 6.4, color='lightblue', alpha=0.1)
ax.axvspan(6.4, GENOME_SIZE / 1e6, color='lightyellow', alpha=0.3)

ax.set_yticks([0.3, 0.7])
ax.set_yticklabels(['Gained', 'Lost'], fontsize=10)
ax.set_xlabel('Chromosome position (Mb)')
ax.set_title('Genes by Methylation Transition Group', fontweight='bold')
ax.set_ylim(0, 1)
ax.set_xlim(0, GENOME_SIZE / 1e6)

plt.tight_layout()
plt.savefig(FIGURES / 'geographic_context.pdf', dpi=300, bbox_inches='tight')
plt.savefig(FIGURES / 'geographic_context.svg', bbox_inches='tight')
plt.close()
print("Saved geographic_context.pdf/svg")

# ── Figure 5: Comprehensive 4-panel summary ──
fig = plt.figure(figsize=(16, 14))
gs = GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.3)

# Panel A: Violin plot T2vsT1
ax = fig.add_subplot(gs[0, 0])
data_list = [lost_lfc_t2t1.values, gained_lfc_t2t1.values]
if len(both_lfc_t2t1) > 0:
    data_list.append(both_lfc_t2t1.values)
data_list.append(never_lfc_t2t1.values)

lab_list = [f'Lost\n(n={len(lost_lfc_t2t1)})', f'Gained\n(n={len(gained_lfc_t2t1)})']
col_list = [colors['Lost'], colors['Gained']]
if len(both_lfc_t2t1) > 0:
    lab_list.append(f'Both\n(n={len(both_lfc_t2t1)})')
    col_list.append(colors['Both'])
lab_list.append(f'Never\n(n={len(never_lfc_t2t1)})')
col_list.append(colors['Never'])

parts = ax.violinplot(data_list, positions=range(len(data_list)),
                      showmeans=True, showmedians=True)
for i, pc in enumerate(parts['bodies']):
    pc.set_facecolor(col_list[i])
    pc.set_alpha(0.7)
parts['cmeans'].set_color('black')
parts['cmedians'].set_color('darkred')

ax.set_xticks(range(len(lab_list)))
ax.set_xticklabels(lab_list, fontsize=9)
ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
ax.set_ylabel('log2 Fold Change (T2/T1)')
ax.set_title('A. Expression by GCCGGC Transition', fontweight='bold')

# Add medians as text
for i, d in enumerate(data_list):
    med = np.median(d)
    ax.text(i, ax.get_ylim()[1] * 0.9, f'med={med:.2f}', ha='center', fontsize=8, color='darkred')

# Panel B: Fraction up/down
ax = fig.add_subplot(gs[0, 1])
for i, name in enumerate(group_names):
    gene_set, deseq = group_data[name]
    sub = deseq[deseq['gene_id'].isin(gene_set)]
    n = len(sub)
    if n == 0:
        continue
    sig = sub[sub['padj'] < 0.05]
    n_up = (sig['log2FoldChange'] > 0).sum()
    n_down = (sig['log2FoldChange'] < 0).sum()
    n_ns = n - n_up - n_down

    ax.bar(i, n_up / n, width, color='#e74c3c', alpha=0.8,
           label='Up' if i == 0 else '')
    ax.bar(i, n_ns / n, width, bottom=n_up / n, color='#bdc3c7', alpha=0.8,
           label='NS' if i == 0 else '')
    ax.bar(i, n_down / n, width, bottom=(n_up + n_ns) / n, color='#3498db', alpha=0.8,
           label='Down' if i == 0 else '')

    ax.text(i, -0.05, f'n={n}', ha='center', fontsize=8)

ax.set_xticks(x_pos)
ax.set_xticklabels(group_names, fontsize=10)
ax.set_ylabel('Fraction')
ax.set_title('B. DEG Fractions by Transition', fontweight='bold')
ax.legend(loc='upper right', fontsize=8)
ax.set_ylim(-0.08, 1.05)

# Panel C: Geographic distribution
ax = fig.add_subplot(gs[1, 0])
for _, site in t1_sites.iterrows():
    ax.axvline(x=site['position'] / 1e6, ymin=0.5, ymax=1, color='#e74c3c', alpha=0.1, linewidth=0.5)
for _, site in t2_sites.iterrows():
    ax.axvline(x=site['position'] / 1e6, ymin=0, ymax=0.5, color='#3498db', alpha=0.2, linewidth=0.8)

ax.axvspan(0, 1.5, color='lightyellow', alpha=0.3)
ax.axvspan(1.5, 6.4, color='lightblue', alpha=0.1)
ax.axvspan(6.4, GENOME_SIZE / 1e6, color='lightyellow', alpha=0.3)

ax.axhline(y=0.5, color='black', linewidth=0.5)
ax.set_yticks([0.25, 0.75])
ax.set_yticklabels(['T2 sites', 'T1 sites'], fontsize=9)
ax.set_xlabel('Chromosome (Mb)')
ax.set_title('C. GCCGGC Site Redistribution', fontweight='bold')
ax.set_xlim(0, GENOME_SIZE / 1e6)

# Panel D: Statistical summary table
ax = fig.add_subplot(gs[1, 1])
ax.axis('off')

# Build summary text
summary_text = "STATISTICAL SUMMARY\n" + "=" * 40 + "\n\n"
for _, row in stat_df.iterrows():
    if pd.isna(row['p_value']):
        continue
    sig_marker = "***" if row['p_bonferroni'] < 0.001 else "**" if row['p_bonferroni'] < 0.01 else "*" if row['p_bonferroni'] < 0.05 else "ns"
    summary_text += f"{row['comparison']}:\n"
    summary_text += f"  p={row['p_value']:.2e} (Bonf. {row['p_bonferroni']:.2e}) {sig_marker}\n"
    summary_text += f"  r={row['rank_biserial']:.3f}, median diff={row['median_diff']:.3f}\n\n"

summary_text += "\nGEOGRAPHIC CONTROLS\n" + "-" * 30 + "\n"
for _, row in geo_df.iterrows():
    if pd.isna(row['p_value']):
        summary_text += f"  {row['region']}: {row['comparison']} - SKIPPED\n"
    else:
        summary_text += f"  {row['region']}: {row['comparison']}\n"
        summary_text += f"    p={row['p_value']:.2e}, r={row['rank_biserial']:.3f}\n"

ax.text(0.05, 0.95, summary_text, transform=ax.transAxes, fontsize=8,
        verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
ax.set_title('D. Test Results', fontweight='bold')

fig.suptitle('H19: GCCGGC Temporal De-repression Analysis', fontsize=15, fontweight='bold', y=0.98)
plt.savefig(FIGURES / 'H19_comprehensive_summary.pdf', dpi=300, bbox_inches='tight')
plt.savefig(FIGURES / 'H19_comprehensive_summary.svg', bbox_inches='tight')
plt.close()
print("Saved H19_comprehensive_summary.pdf/svg")

# ══════════════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("FINAL SUMMARY")
print("=" * 70)

print(f"\nT1→T2 GCCGGC transition groups:")
print(f"  Lost:   n={len(lost_expr):>5}, median T2vsT1 LFC = {lost_lfc_t2t1.median():>+.4f}")
print(f"  Gained: n={len(gained_expr):>5}, median T2vsT1 LFC = {gained_lfc_t2t1.median():>+.4f}")
if len(both_lfc_t2t1) > 0:
    print(f"  Both:   n={len(both_expr):>5}, median T2vsT1 LFC = {both_lfc_t2t1.median():>+.4f}")
print(f"  Never:  n={len(never_expr):>5}, median T2vsT1 LFC = {never_lfc_t2t1.median():>+.4f}")

print(f"\nKey statistical tests:")
for _, row in stat_df.iterrows():
    if pd.isna(row['p_value']):
        continue
    sig = "***" if row['p_bonferroni'] < 0.001 else "**" if row['p_bonferroni'] < 0.01 else "*" if row['p_bonferroni'] < 0.05 else "ns"
    print(f"  {row['comparison']:>20s}: p_bonf={row['p_bonferroni']:.2e} {sig}, r={row['rank_biserial']:+.3f}")

# Verdict
print("\n" + "=" * 70)
lost_med = lost_lfc_t2t1.median()
gained_med = gained_lfc_t2t1.median()
never_med = never_lfc_t2t1.median()

# Check main hypothesis criteria
lost_vs_never_row = stat_df[stat_df['comparison'] == 'Lost vs Never'].iloc[0]
lost_vs_gained_row = stat_df[stat_df['comparison'] == 'Lost vs Gained'].iloc[0]
gained_vs_never_row = stat_df[stat_df['comparison'] == 'Gained vs Never'].iloc[0]

# Criteria for SUPPORTED:
# 1. Lost genes have higher LFC than Never (de-repression)
# 2. Statistical significance after Bonferroni
# 3. Gained genes have lower LFC than Never (suppression)

crit1 = lost_med > never_med  # Lost > Never
crit2 = lost_vs_never_row['p_bonferroni'] < 0.05  # Significant
crit3 = gained_med < never_med  # Gained < Never (suppression)

# Check geographic confound
geo_test = geo_df[(geo_df['comparison'] == 'Lost vs Never') & (~geo_df['p_value'].isna())]
geo_passes = False
for _, gr in geo_test.iterrows():
    if gr['p_value'] < 0.05:
        geo_passes = True

if crit1 and crit2 and crit3 and geo_passes:
    verdict = "SUPPORTED"
    evidence = (f"Lost genes show significant upregulation vs Never (median LFC diff={lost_med - never_med:+.3f}, "
                f"p_bonf={lost_vs_never_row['p_bonferroni']:.2e}). "
                f"Gained genes show suppression (median LFC={gained_med:+.3f}). "
                f"Effect survives geographic stratification.")
elif crit1 and crit2:
    if geo_passes:
        verdict = "PARTIAL"
        evidence = (f"Lost genes show significant upregulation vs Never (p_bonf={lost_vs_never_row['p_bonferroni']:.2e}), "
                    f"but gained genes do not show expected suppression (median LFC={gained_med:+.3f} vs Never={never_med:+.3f}).")
    else:
        verdict = "PARTIAL"
        evidence = (f"Lost genes show significant upregulation vs Never in overall test "
                    f"(p_bonf={lost_vs_never_row['p_bonferroni']:.2e}), "
                    f"but effect does not survive geographic stratification, suggesting possible confounding.")
elif crit1 and not crit2:
    verdict = "PARTIAL"
    evidence = (f"Lost genes show directionally correct upregulation (median LFC={lost_med:+.3f} vs Never={never_med:+.3f}), "
                f"but does not reach significance after Bonferroni correction (p_bonf={lost_vs_never_row['p_bonferroni']:.2e}).")
else:
    verdict = "REJECTED"
    evidence = (f"Lost genes do NOT show expected upregulation. "
                f"Median LFC: Lost={lost_med:+.3f}, Never={never_med:+.3f}.")

print(f"VERDICT: {verdict}")
print(f"EVIDENCE: {evidence}")
print("=" * 70)

# Save verdict for report
with open(TABLES / "verdict.txt", "w") as f:
    f.write(f"VERDICT: {verdict}\n")
    f.write(f"EVIDENCE: {evidence}\n")

print("\nDone.")
