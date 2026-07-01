#!/usr/bin/env python3
"""
H28: Biological Characteristics of 57 Exposed Regulators
=========================================================
Comprehensive analysis of why 57 regulatory genes lack the 1,200 bp protection zone
that shields the other 998 regulators from methylation influence.

Core hypothesis: Exposed regulators have low T1 expression (RNAP not occupying promoter).
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import mannwhitneyu, chi2_contingency, fisher_exact, kruskal, spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from Bio import SeqIO
import warnings

warnings.filterwarnings('ignore')


def rank_biserial_r(U, n1, n2):
    """Rank-biserial correlation r from Mann-Whitney U.
    r = 1 - (2*U)/(n1*n2).  Ranges from -1 to +1.
    """
    return 1.0 - (2.0 * U) / (n1 * n2)


# ============================================================
# Configuration
# ============================================================
BASE = '/Users/okaban/bioinfo/rna-seq'
ANALYSIS_DIR = f'{BASE}/11_epigenome_integration/analysis/51_exposed_regulators_characteristics'
FIG_DIR = f'{ANALYSIS_DIR}/figures'
TAB_DIR = f'{ANALYSIS_DIR}/tables'

COORD_REG = f'{BASE}/11_epigenome_integration/analysis/29_genomewide_TF_screen/tables/coordinated_regulatory_genes.tsv'
ALL_REG = f'{BASE}/11_epigenome_integration/analysis/29_genomewide_TF_screen/tables/all_regulatory_genes.tsv'
DESEQ_T2 = f'{BASE}/04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_2_vs_1.tsv'
DESEQ_T3 = f'{BASE}/04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_3_vs_1.tsv'
NORM_COUNTS = f'{BASE}/04_deseq2/analysis/04_deseq2_260128_v1/results/normalized_counts_M145.tsv'
GENE_ANNOT = f'{BASE}/05_annotation/analysis/05_annotation_260128_v1/tables/gene_annotation_basic.tsv'
H27_METRICS = f'{BASE}/11_epigenome_integration/analysis/50_coordinated_regulators_protection/tables/gene_level_metrics.tsv'
FIMO_GFF = f'{BASE}/13_TF_binding-site/analysis/01_master_TF_list_260206_v1/intermediate/fimo_results/fimo.gff'
GENOME_FASTA = '/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/GCF_000203835.1_ASM20383v1_genomic.fna'

CHR_LEN = 8_667_507
ARM_LEFT_END = 1_500_000
ARM_RIGHT_START = 7_167_508

plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

COLOR_EXPOSED = '#E63946'
COLOR_SHIELDED = '#457B9D'


def region_label(pos):
    if pos <= ARM_LEFT_END or pos >= ARM_RIGHT_START:
        return 'arm'
    return 'core'


def sig_stars(p):
    if p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    return 'ns'


def compare_groups(exp_vals, shi_vals, label):
    """Run Wilcoxon rank-sum and print results."""
    exp_v = exp_vals.dropna()
    shi_v = shi_vals.dropna()
    U, p = mannwhitneyu(exp_v, shi_v, alternative='two-sided')
    r = rank_biserial_r(U, len(exp_v), len(shi_v))
    direction = 'exposed < shielded' if exp_v.median() < shi_v.median() else 'exposed > shielded'
    print(f"\n{label}:")
    print(f"  Exposed (n={len(exp_v)}):  median={exp_v.median():.2f}, mean={exp_v.mean():.2f}")
    print(f"  Shielded (n={len(shi_v)}): median={shi_v.median():.2f}, mean={shi_v.mean():.2f}")
    print(f"  Wilcoxon U={U:.0f}, p={p:.2e}, r={r:.3f} ({direction})")
    return {'feature': label, 'exposed_n': len(exp_v), 'exposed_median': exp_v.median(),
            'exposed_mean': exp_v.mean(), 'shielded_n': len(shi_v),
            'shielded_median': shi_v.median(), 'shielded_mean': shi_v.mean(),
            'U': U, 'p_value': p, 'effect_size_r': r, 'test': 'Wilcoxon',
            'direction': direction, 'significant': sig_stars(p)}


# ============================================================
# 1. Load Data
# ============================================================
print("=" * 70)
print("H28: Biological Characteristics of 57 Exposed Regulators")
print("=" * 70)

coord = pd.read_csv(COORD_REG, sep='\t')
exposed_loci = set(coord['locus_tag'].tolist())
print(f"\nExposed regulators: {len(exposed_loci)}")

all_reg = pd.read_csv(ALL_REG, sep='\t')
all_reg_loci = set(all_reg['locus_tag'].tolist())
shielded_loci = all_reg_loci - exposed_loci
print(f"Total regulators: {len(all_reg_loci)}")
print(f"Shielded regulators: {len(shielded_loci)}")

deseq_t2 = pd.read_csv(DESEQ_T2, sep='\t')
deseq_t2.columns = ['gene_id', 'baseMean', 'log2FC_T2', 'lfcSE_T2', 'pvalue_T2', 'padj_T2']

deseq_t3 = pd.read_csv(DESEQ_T3, sep='\t')
deseq_t3.columns = ['gene_id', 'baseMean_T3', 'log2FC_T3', 'lfcSE_T3', 'pvalue_T3', 'padj_T3']

norm_counts = pd.read_csv(NORM_COUNTS, sep='\t')
t1_cols = [c for c in norm_counts.columns if c.startswith('M145_1_')]
t2_cols = [c for c in norm_counts.columns if c.startswith('M145_2_')]
t3_cols = [c for c in norm_counts.columns if c.startswith('M145_3_')]
norm_counts['T1_mean'] = norm_counts[t1_cols].mean(axis=1)
norm_counts['T2_mean'] = norm_counts[t2_cols].mean(axis=1)
norm_counts['T3_mean'] = norm_counts[t3_cols].mean(axis=1)

annot = pd.read_csv(GENE_ANNOT, sep='\t')
h27 = pd.read_csv(H27_METRICS, sep='\t')

# ============================================================
# Build master table
# ============================================================
master = all_reg[['locus_tag', 'gene_name', 'old_locus_tag', 'product', 'tf_family',
                   'start', 'end', 'strand']].copy()
master['is_exposed'] = master['locus_tag'].isin(exposed_loci)
master['gene_length'] = master['end'] - master['start'] + 1

master = master.merge(deseq_t2[['gene_id', 'baseMean', 'log2FC_T2', 'padj_T2']],
                       left_on='locus_tag', right_on='gene_id', how='left').drop(columns=['gene_id'])
master = master.merge(deseq_t3[['gene_id', 'log2FC_T3', 'padj_T3']],
                       left_on='locus_tag', right_on='gene_id', how='left').drop(columns=['gene_id'])
master = master.merge(norm_counts[['gene_id', 'T1_mean', 'T2_mean', 'T3_mean']],
                       left_on='locus_tag', right_on='gene_id', how='left').drop(columns=['gene_id'])

h27_cols = ['locus_tag', 'nearest_methyl_dist', 'sites_within_2kb', 'sites_2kb_T1',
            'sites_2kb_T2', 'delta_sites']
master = master.merge(h27[h27_cols], on='locus_tag', how='left')

master['midpoint'] = (master['start'] + master['end']) / 2
master['region'] = master['midpoint'].apply(region_label)

coord_types = coord[['locus_tag', 'coordination_T2', 'coordination_T3']].copy()
master = master.merge(coord_types, on='locus_tag', how='left')

# Convenience subsets
exposed_df = master[master['is_exposed']]
shielded_df = master[~master['is_exposed']]

# Collect all stats rows
all_stats = []

# ============================================================
# STEP 1: Expression Level Comparison (PRIMARY HYPOTHESIS TEST)
# ============================================================
print("\n" + "=" * 70)
print("STEP 1: Expression Level Comparison (PRIMARY HYPOTHESIS TEST)")
print("=" * 70)

all_stats.append(compare_groups(exposed_df['baseMean'], shielded_df['baseMean'], 'baseMean'))
all_stats.append(compare_groups(exposed_df['T1_mean'], shielded_df['T1_mean'], 'T1_expression'))
all_stats.append(compare_groups(exposed_df['T2_mean'], shielded_df['T2_mean'], 'T2_expression'))
all_stats.append(compare_groups(exposed_df['T3_mean'], shielded_df['T3_mean'], 'T3_expression'))
all_stats.append(compare_groups(exposed_df['log2FC_T2'], shielded_df['log2FC_T2'], 'log2FC_T2vsT1'))
all_stats.append(compare_groups(exposed_df['log2FC_T3'], shielded_df['log2FC_T3'], 'log2FC_T3vsT1'))
all_stats.append(compare_groups(exposed_df['log2FC_T2'].dropna().abs(), shielded_df['log2FC_T2'].dropna().abs(), '|log2FC_T2vsT1|'))
all_stats.append(compare_groups(exposed_df['log2FC_T3'].dropna().abs(), shielded_df['log2FC_T3'].dropna().abs(), '|log2FC_T3vsT1|'))

# ============================================================
# STEP 2: TF Family Composition
# ============================================================
print("\n" + "=" * 70)
print("STEP 2: TF Family Composition")
print("=" * 70)

family_exposed = exposed_df['tf_family'].value_counts()
family_shielded = shielded_df['tf_family'].value_counts()
all_families = sorted(set(master['tf_family'].unique()))

fisher_results = []
for fam in all_families:
    n_exp_fam = family_exposed.get(fam, 0)
    n_shi_fam = family_shielded.get(fam, 0)
    n_exp_other = 57 - n_exp_fam
    n_shi_other = 998 - n_shi_fam
    table_2x2 = [[n_exp_fam, n_shi_fam], [n_exp_other, n_shi_other]]
    odds, pval = fisher_exact(table_2x2)
    fisher_results.append({
        'tf_family': fam, 'n_exposed': int(n_exp_fam), 'n_shielded': int(n_shi_fam),
        'n_total': int(n_exp_fam + n_shi_fam),
        'pct_exposed': n_exp_fam / (n_exp_fam + n_shi_fam) * 100 if (n_exp_fam + n_shi_fam) > 0 else 0,
        'expected_pct': 57 / 1055 * 100,
        'odds_ratio': odds, 'fisher_p': pval,
        'enrichment': 'enriched' if odds > 1 else 'depleted' if odds < 1 else 'neutral'
    })

fisher_df = pd.DataFrame(fisher_results).sort_values('fisher_p')
print(f"\nExpected exposed %: {57/1055*100:.1f}%")
print(fisher_df[['tf_family', 'n_exposed', 'n_total', 'pct_exposed', 'odds_ratio', 'fisher_p', 'enrichment']].to_string(index=False))

fisher_df.to_csv(f'{TAB_DIR}/TF_family_distribution.tsv', sep='\t', index=False)
print(f"\nSaved: TF_family_distribution.tsv")

# ============================================================
# STEP 3: Gene Structural Features
# ============================================================
print("\n" + "=" * 70)
print("STEP 3: Gene Structural Features")
print("=" * 70)

# Gene length
all_stats.append(compare_groups(exposed_df['gene_length'], shielded_df['gene_length'], 'gene_length'))

# GC content of promoter region
print("\nComputing promoter GC content from genome FASTA...")
genome = SeqIO.to_dict(SeqIO.parse(GENOME_FASTA, 'fasta'))
chr_seq = str(genome['NC_003888.3'].seq)

def get_promoter_gc(row, seq, window=300):
    if row['strand'] == '+':
        tss = int(row['start'])
        prom_start = max(0, tss - window)
        prom_end = tss
    else:
        tss = int(row['end'])
        prom_start = tss
        prom_end = min(len(seq), tss + window)
    prom_seq = seq[prom_start:prom_end].upper()
    if len(prom_seq) == 0:
        return np.nan
    return (prom_seq.count('G') + prom_seq.count('C')) / len(prom_seq)

master['gc_promoter'] = master.apply(lambda r: get_promoter_gc(r, chr_seq), axis=1)
all_stats.append(compare_groups(master.loc[master['is_exposed'], 'gc_promoter'],
                                 master.loc[~master['is_exposed'], 'gc_promoter'], 'gc_promoter'))

# Intergenic distance
all_genes_sorted = annot.sort_values('start').reset_index(drop=True)
all_genes_sorted['prev_end'] = all_genes_sorted['end'].shift(1)
all_genes_sorted['intergenic_dist'] = all_genes_sorted['start'] - all_genes_sorted['prev_end']
all_genes_sorted.loc[all_genes_sorted['intergenic_dist'] < 0, 'intergenic_dist'] = 0
intergenic_map = dict(zip(all_genes_sorted['gene_id'], all_genes_sorted['intergenic_dist']))
master['intergenic_dist'] = master['locus_tag'].map(intergenic_map)
all_stats.append(compare_groups(master.loc[master['is_exposed'], 'intergenic_dist'],
                                 master.loc[~master['is_exposed'], 'intergenic_dist'], 'intergenic_dist'))

# FIMO TF binding site predictions
print("\nLoading FIMO TF binding sites...")
fimo_sites = []
with open(FIMO_GFF) as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) >= 5:
            fimo_sites.append({'start': int(parts[3]), 'end': int(parts[4])})
fimo_df = pd.DataFrame(fimo_sites)
print(f"  Loaded {len(fimo_df)} FIMO sites")

def count_fimo_in_promoter(row, fimo, window_up=500, window_down=100):
    if row['strand'] == '+':
        tss = int(row['start'])
        region_start = tss - window_up
        region_end = tss + window_down
    else:
        tss = int(row['end'])
        region_start = tss - window_down
        region_end = tss + window_up
    return len(fimo[(fimo['start'] >= region_start) & (fimo['end'] <= region_end)])

master['n_fimo_hits'] = master.apply(lambda r: count_fimo_in_promoter(r, fimo_df), axis=1)
all_stats.append(compare_groups(master.loc[master['is_exposed'], 'n_fimo_hits'],
                                 master.loc[~master['is_exposed'], 'n_fimo_hits'], 'n_fimo_hits'))

# ============================================================
# STEP 4: Geographic Distribution
# ============================================================
print("\n" + "=" * 70)
print("STEP 4: Geographic Distribution")
print("=" * 70)

exp_region = master[master['is_exposed']]['region'].value_counts()
shi_region = master[~master['is_exposed']]['region'].value_counts()
n_exp_arm = exp_region.get('arm', 0)
n_exp_core = exp_region.get('core', 0)
n_shi_arm = shi_region.get('arm', 0)
n_shi_core = shi_region.get('core', 0)

geo_table = [[n_exp_arm, n_exp_core], [n_shi_arm, n_shi_core]]
geo_or, geo_p = fisher_exact(geo_table)

print(f"\nExposed: {n_exp_arm} arm ({n_exp_arm/57*100:.1f}%), {n_exp_core} core ({n_exp_core/57*100:.1f}%)")
print(f"Shielded: {n_shi_arm} arm ({n_shi_arm/998*100:.1f}%), {n_shi_core} core ({n_shi_core/998*100:.1f}%)")
print(f"Fisher OR={geo_or:.2f}, p={geo_p:.3f}")

all_stats.append({'feature': 'geographic (arm%)', 'exposed_n': 57, 'exposed_median': n_exp_arm / 57 * 100,
                  'exposed_mean': np.nan, 'shielded_n': 998, 'shielded_median': n_shi_arm / 998 * 100,
                  'shielded_mean': np.nan, 'U': geo_or, 'p_value': geo_p, 'effect_size_r': np.nan,
                  'test': 'Fisher', 'direction': f'OR={geo_or:.2f}', 'significant': sig_stars(geo_p)})

# ============================================================
# STEP 5: Methylation Environment
# ============================================================
print("\n" + "=" * 70)
print("STEP 5: Methylation Environment")
print("=" * 70)

all_stats.append(compare_groups(master.loc[master['is_exposed'], 'sites_within_2kb'],
                                 master.loc[~master['is_exposed'], 'sites_within_2kb'], 'sites_within_2kb'))
all_stats.append(compare_groups(master.loc[master['is_exposed'], 'nearest_methyl_dist'],
                                 master.loc[~master['is_exposed'], 'nearest_methyl_dist'], 'nearest_methyl_dist'))

master['methyl_density_2kb'] = master['sites_within_2kb'] / 4.0
all_stats.append(compare_groups(master.loc[master['is_exposed'], 'methyl_density_2kb'],
                                 master.loc[~master['is_exposed'], 'methyl_density_2kb'], 'methyl_density_2kb'))

# ============================================================
# STEP 6: Operon Context
# ============================================================
print("\n" + "=" * 70)
print("STEP 6: Operon Context")
print("=" * 70)

all_genes_sorted2 = annot.sort_values('start').reset_index(drop=True)
operon_id = 0
operon_map = {}
for i, row in all_genes_sorted2.iterrows():
    if i == 0:
        operon_map[row['gene_id']] = operon_id
        continue
    prev = all_genes_sorted2.loc[i - 1]
    if (row['strand'] == prev['strand'] and
        0 <= row['start'] - prev['end'] <= 200):
        operon_map[row['gene_id']] = operon_id
    else:
        operon_id += 1
        operon_map[row['gene_id']] = operon_id

operon_sizes = pd.Series(operon_map).value_counts().to_dict()
gene_operon_size = {gene: operon_sizes[op] for gene, op in operon_map.items()}
master['operon_size'] = master['locus_tag'].map(gene_operon_size)
master['is_monocistronic'] = master['operon_size'] == 1

n_exp_mono = master.loc[master['is_exposed'], 'is_monocistronic'].sum()
n_exp_poly = 57 - n_exp_mono
n_shi_mono = master.loc[~master['is_exposed'], 'is_monocistronic'].sum()
n_shi_poly = 998 - n_shi_mono

operon_tab = [[n_exp_mono, n_exp_poly], [n_shi_mono, n_shi_poly]]
operon_or, operon_p = fisher_exact(operon_tab)

print(f"\nExposed: {n_exp_mono} monocistronic ({n_exp_mono/57*100:.1f}%), {n_exp_poly} polycistronic ({n_exp_poly/57*100:.1f}%)")
print(f"Shielded: {n_shi_mono} monocistronic ({n_shi_mono/998*100:.1f}%), {n_shi_poly} polycistronic ({n_shi_poly/998*100:.1f}%)")
print(f"Fisher OR={operon_or:.2f}, p={operon_p:.3f}")

all_stats.append({'feature': 'monocistronic (%)', 'exposed_n': 57, 'exposed_median': n_exp_mono / 57 * 100,
                  'exposed_mean': np.nan, 'shielded_n': 998, 'shielded_median': n_shi_mono / 998 * 100,
                  'shielded_mean': np.nan, 'U': operon_or, 'p_value': operon_p, 'effect_size_r': np.nan,
                  'test': 'Fisher', 'direction': f'OR={operon_or:.2f}', 'significant': sig_stars(operon_p)})

all_stats.append(compare_groups(master.loc[master['is_exposed'], 'operon_size'],
                                 master.loc[~master['is_exposed'], 'operon_size'], 'operon_size'))

# ============================================================
# STEP 7: Temporal Expression Pattern
# ============================================================
print("\n" + "=" * 70)
print("STEP 7: Temporal Expression Pattern")
print("=" * 70)

def classify_temporal_pattern(row):
    lfc2 = row.get('log2FC_T2')
    lfc3 = row.get('log2FC_T3')
    padj2 = row.get('padj_T2')
    padj3 = row.get('padj_T3')
    if pd.isna(lfc2) or pd.isna(lfc3):
        return 'no_data'
    sig2 = not pd.isna(padj2) and padj2 < 0.05
    sig3 = not pd.isna(padj3) and padj3 < 0.05
    up2 = sig2 and lfc2 > 1
    down2 = sig2 and lfc2 < -1
    up3 = sig3 and lfc3 > 1
    down3 = sig3 and lfc3 < -1
    if not up2 and not down2 and not up3 and not down3:
        return 'constitutive'
    elif up2 and up3:
        return 'T2T3_up'
    elif down2 and down3:
        return 'T2T3_down'
    elif up2 and not up3 and not down3:
        return 'T2_up'
    elif down2 and not up3 and not down3:
        return 'T2_down'
    elif up3 and not up2 and not down2:
        return 'T3_up'
    elif down3 and not up2 and not down2:
        return 'T3_down'
    elif up2 and down3:
        return 'T2_up_T3_down'
    elif down2 and up3:
        return 'T2_down_T3_up'
    return 'other_dynamic'

master['temporal_pattern'] = master.apply(classify_temporal_pattern, axis=1)

exp_temp = master[master['is_exposed']]['temporal_pattern'].value_counts()
shi_temp = master[~master['is_exposed']]['temporal_pattern'].value_counts()

pattern_order = ['constitutive', 'T2_up', 'T3_up', 'T2T3_up', 'T2_down', 'T3_down',
                 'T2T3_down', 'T2_up_T3_down', 'T2_down_T3_up', 'other_dynamic', 'no_data']

temp_table = pd.DataFrame(index=pattern_order)
temp_table['n_exposed'] = exp_temp.reindex(pattern_order, fill_value=0)
temp_table['n_shielded'] = shi_temp.reindex(pattern_order, fill_value=0)
temp_table['pct_exposed'] = (temp_table['n_exposed'] / 57 * 100).round(1)
temp_table['pct_shielded'] = (temp_table['n_shielded'] / 998 * 100).round(1)

n_exp_const = temp_table.loc['constitutive', 'n_exposed']
n_exp_nodata = temp_table.loc['no_data', 'n_exposed'] if 'no_data' in temp_table.index else 0
n_exp_dynamic = 57 - n_exp_const - n_exp_nodata
n_shi_const = temp_table.loc['constitutive', 'n_shielded']
n_shi_nodata = temp_table.loc['no_data', 'n_shielded'] if 'no_data' in temp_table.index else 0
n_shi_dynamic = 998 - n_shi_const - n_shi_nodata

const_table = [[n_exp_dynamic, n_exp_const], [n_shi_dynamic, n_shi_const]]
const_or, const_p = fisher_exact(const_table)

print(f"\nTemporal expression patterns:")
print(temp_table.to_string())

denom_exp = 57 - n_exp_nodata
denom_shi = 998 - n_shi_nodata
print(f"\nConstitutive vs Dynamic (excluding no_data):")
print(f"  Exposed: {n_exp_dynamic} dynamic ({n_exp_dynamic/denom_exp*100:.1f}%), {n_exp_const} constitutive ({n_exp_const/denom_exp*100:.1f}%)")
print(f"  Shielded: {n_shi_dynamic} dynamic ({n_shi_dynamic/denom_shi*100:.1f}%), {n_shi_const} constitutive ({n_shi_const/denom_shi*100:.1f}%)")
print(f"  Fisher OR={const_or:.2f}, p={const_p:.2e}")

all_stats.append({'feature': 'dynamic_expression (%)', 'exposed_n': denom_exp,
                  'exposed_median': n_exp_dynamic / denom_exp * 100,
                  'exposed_mean': np.nan, 'shielded_n': denom_shi,
                  'shielded_median': n_shi_dynamic / denom_shi * 100,
                  'shielded_mean': np.nan, 'U': const_or, 'p_value': const_p,
                  'effect_size_r': np.nan, 'test': 'Fisher',
                  'direction': f'OR={const_or:.2f}', 'significant': sig_stars(const_p)})

temp_table.index.name = 'temporal_pattern'
temp_table.to_csv(f'{TAB_DIR}/temporal_pattern_distribution.tsv', sep='\t')
print(f"\nSaved: temporal_pattern_distribution.tsv")

# ============================================================
# STEP 8: Statistics Summary
# ============================================================
print("\n" + "=" * 70)
print("STEP 8: Feature Comparison Statistics Summary")
print("=" * 70)

stats_df = pd.DataFrame(all_stats)
stats_df.to_csv(f'{TAB_DIR}/feature_comparison_statistics.tsv', sep='\t', index=False)
print("\n" + stats_df[['feature', 'exposed_median', 'shielded_median', 'p_value', 'effect_size_r', 'significant']].to_string(index=False))

# ============================================================
# Save Full Table of 57 Exposed Regulators
# ============================================================
exposed_full = master[master['is_exposed']].copy()
save_cols = ['locus_tag', 'gene_name', 'old_locus_tag', 'product', 'tf_family',
             'start', 'end', 'strand', 'region', 'baseMean', 'T1_mean', 'T2_mean', 'T3_mean',
             'log2FC_T2', 'log2FC_T3', 'padj_T2', 'padj_T3',
             'gene_length', 'gc_promoter', 'intergenic_dist', 'n_fimo_hits',
             'nearest_methyl_dist', 'sites_within_2kb', 'methyl_density_2kb',
             'operon_size', 'is_monocistronic', 'temporal_pattern',
             'coordination_T2', 'coordination_T3']
save_cols = [c for c in save_cols if c in exposed_full.columns]
exposed_full[save_cols].to_csv(f'{TAB_DIR}/exposed_regulators_full_table.tsv', sep='\t', index=False)
print(f"\nSaved: exposed_regulators_full_table.tsv ({len(exposed_full)} genes)")

# ============================================================
# T1 expression quartile analysis
# ============================================================
print("\n" + "=" * 70)
print("ADDITIONAL: T1 Expression Quartile Analysis")
print("=" * 70)

master_t1 = master[master['T1_mean'].notna()].copy()
master_t1['T1_quartile'] = pd.qcut(master_t1['T1_mean'], q=4, labels=['Q1 (lowest)', 'Q2', 'Q3', 'Q4 (highest)'])
q_table = master_t1.groupby('T1_quartile').agg(
    n_total=('is_exposed', 'count'),
    n_exposed=('is_exposed', 'sum')
).reset_index()
q_table['pct_exposed'] = (q_table['n_exposed'] / q_table['n_total'] * 100).round(1)
print("\nT1 expression quartile breakdown:")
print(q_table.to_string(index=False))

rho_trend, p_trend = spearmanr(master_t1['T1_quartile'].cat.codes, master_t1['is_exposed'].astype(int))
print(f"\nTrend (Spearman): rho={rho_trend:.3f}, p={p_trend:.2e}")

# ============================================================
# FIGURES
# ============================================================
print("\n" + "=" * 70)
print("GENERATING FIGURES")
print("=" * 70)


def make_violin_box(ax, data_list, labels, title, ylabel, p_val):
    """Helper: violin+box plot for 2 groups."""
    vp = ax.violinplot(data_list, positions=[1, 2], showmedians=True, showextrema=False)
    vp['bodies'][0].set_facecolor(COLOR_EXPOSED)
    vp['bodies'][0].set_alpha(0.7)
    vp['bodies'][1].set_facecolor(COLOR_SHIELDED)
    vp['bodies'][1].set_alpha(0.7)
    vp['cmedians'].set_color('black')
    bp = ax.boxplot(data_list, positions=[1, 2], widths=0.15,
                    showfliers=False, patch_artist=True, zorder=3)
    bp['boxes'][0].set_facecolor(COLOR_EXPOSED)
    bp['boxes'][0].set_alpha(0.9)
    bp['boxes'][1].set_facecolor(COLOR_SHIELDED)
    bp['boxes'][1].set_alpha(0.9)
    for element in ['whiskers', 'caps']:
        for item in bp[element]:
            item.set_color('black')
    for m in bp['medians']:
        m.set_color('white')
    ax.set_xticks([1, 2])
    ax.set_xticklabels(labels)
    ax.set_title(f'{title}\np={p_val:.1e} ({sig_stars(p_val)})', fontweight='bold')
    ax.set_ylabel(ylabel)


# -- Get values for plots --
exp_bm = exposed_df['baseMean'].dropna().values
shi_bm = shielded_df['baseMean'].dropna().values
exp_t1v = exposed_df['T1_mean'].dropna().values
shi_t1v = shielded_df['T1_mean'].dropna().values
exp_abs2 = exposed_df['log2FC_T2'].dropna().abs().values
shi_abs2 = shielded_df['log2FC_T2'].dropna().abs().values
exp_abs3 = exposed_df['log2FC_T3'].dropna().abs().values
shi_abs3 = shielded_df['log2FC_T3'].dropna().abs().values

# Get p-values from stats_df
p_bm = stats_df.loc[stats_df['feature'] == 'baseMean', 'p_value'].values[0]
p_t1 = stats_df.loc[stats_df['feature'] == 'T1_expression', 'p_value'].values[0]
p_abs2 = stats_df.loc[stats_df['feature'] == '|log2FC_T2vsT1|', 'p_value'].values[0]
p_abs3 = stats_df.loc[stats_df['feature'] == '|log2FC_T3vsT1|', 'p_value'].values[0]

labels = ['Exposed\n(n=57)', 'Shielded\n(n=998)']

# Panel A: Expression comparison
fig, axes = plt.subplots(1, 4, figsize=(16, 5))
make_violin_box(axes[0], [exp_bm, shi_bm], labels, 'baseMean\n(overall expression)', 'baseMean', p_bm)
make_violin_box(axes[1], [exp_t1v, shi_t1v], labels, 'T1 expression\n(CORE PREDICTION)', 'Normalized counts (T1)', p_t1)
make_violin_box(axes[2], [exp_abs2, shi_abs2], labels, '|log2FC| T2vsT1\n(change magnitude)', '|log2FC|', p_abs2)
make_violin_box(axes[3], [exp_abs3, shi_abs3], labels, '|log2FC| T3vsT1\n(change magnitude)', '|log2FC|', p_abs3)
fig.suptitle('H28: Expression Comparison - Exposed vs Shielded Regulators', fontweight='bold', fontsize=13, y=1.02)
plt.tight_layout()
fig.savefig(f'{FIG_DIR}/expression_comparison.pdf')
fig.savefig(f'{FIG_DIR}/expression_comparison.svg')
print("Saved: expression_comparison.pdf/svg")
plt.close()

# Panel B: TF family composition
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
family_plot = fisher_df.sort_values('n_total', ascending=True)
y_pos = np.arange(len(family_plot))
bar_h = 0.35
ax1.barh(y_pos - bar_h / 2, family_plot['n_exposed'], bar_h,
         color=COLOR_EXPOSED, label='Exposed (57)', edgecolor='white')
ax1.barh(y_pos + bar_h / 2, family_plot['n_shielded'] / 998 * 57, bar_h,
         color=COLOR_SHIELDED, alpha=0.5, label='Shielded (scaled to n=57)', edgecolor='white')
ax1.set_yticks(y_pos)
ax1.set_yticklabels(family_plot['tf_family'], fontsize=9)
ax1.set_xlabel('Count (shielded scaled to n=57)')
ax1.set_title('TF Family Counts', fontweight='bold')
ax1.legend(fontsize=8, loc='lower right')

colors_bar = [COLOR_EXPOSED if r > 1 else COLOR_SHIELDED for r in family_plot['odds_ratio']]
ax2.barh(y_pos, family_plot['pct_exposed'], color=colors_bar, edgecolor='white')
ax2.axvline(57 / 1055 * 100, color='black', linestyle='--', linewidth=1, label=f'Expected ({57/1055*100:.1f}%)')
ax2.set_yticks(y_pos)
ax2.set_yticklabels(family_plot['tf_family'], fontsize=9)
ax2.set_xlabel('% Exposed (of family total)')
ax2.set_title('% Exposed per Family', fontweight='bold')
ax2.legend(fontsize=8)

for i, (_, row) in enumerate(family_plot.iterrows()):
    if row['fisher_p'] < 0.05:
        ax2.text(row['pct_exposed'] + 0.5, i, '*', fontsize=12, fontweight='bold', va='center')

fig.suptitle('H28: TF Family Composition - Exposed vs Shielded', fontweight='bold', fontsize=13, y=1.02)
plt.tight_layout()
fig.savefig(f'{FIG_DIR}/TF_family_composition.pdf')
fig.savefig(f'{FIG_DIR}/TF_family_composition.svg')
print("Saved: TF_family_composition.pdf/svg")
plt.close()

# Panel C: Gene features comparison
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
feature_pairs = [
    (axes[0, 0], 'gene_length', 'Gene Length (bp)'),
    (axes[0, 1], 'gc_promoter', 'Promoter GC Content'),
    (axes[0, 2], 'intergenic_dist', 'Intergenic Distance (bp)'),
    (axes[1, 0], 'n_fimo_hits', 'FIMO TF BS Hits (promoter)'),
    (axes[1, 1], 'sites_within_2kb', 'Methylation Sites in 2kb'),
    (axes[1, 2], 'operon_size', 'Operon Size'),
]
for ax, col, title in feature_pairs:
    ev = master.loc[master['is_exposed'], col].dropna().values
    sv = master.loc[~master['is_exposed'], col].dropna().values
    pv = stats_df.loc[stats_df['feature'] == col, 'p_value'].values
    pv = pv[0] if len(pv) > 0 else 1.0
    make_violin_box(ax, [ev, sv], labels, title, col, pv)

fig.suptitle('H28: Gene Structural Features - Exposed vs Shielded', fontweight='bold', fontsize=13, y=1.02)
plt.tight_layout()
fig.savefig(f'{FIG_DIR}/gene_features_comparison.pdf')
fig.savefig(f'{FIG_DIR}/gene_features_comparison.svg')
print("Saved: gene_features_comparison.pdf/svg")
plt.close()

# Panel D: Chromosome map
fig, ax = plt.subplots(figsize=(14, 4))
shi_pos = master[~master['is_exposed']]['midpoint']
exp_pos = master[master['is_exposed']]['midpoint']
ax.scatter(shi_pos / 1e6, np.zeros(len(shi_pos)), c=COLOR_SHIELDED, alpha=0.3, s=8, label='Shielded (998)')
ax.scatter(exp_pos / 1e6, np.ones(len(exp_pos)) * 0.5, c=COLOR_EXPOSED, s=30, zorder=5,
           label='Exposed (57)', edgecolors='black', linewidth=0.5)
ax.axvline(ARM_LEFT_END / 1e6, color='gray', linestyle='--', alpha=0.5)
ax.axvline(ARM_RIGHT_START / 1e6, color='gray', linestyle='--', alpha=0.5)
ax.fill_betweenx([-0.5, 1.2], 0, ARM_LEFT_END / 1e6, color='lightyellow', alpha=0.3)
ax.fill_betweenx([-0.5, 1.2], ARM_RIGHT_START / 1e6, CHR_LEN / 1e6, color='lightyellow', alpha=0.3)
ax.text(ARM_LEFT_END / 2 / 1e6, 1.1, 'Left Arm', ha='center', fontsize=9, color='gray')
ax.text((ARM_LEFT_END + ARM_RIGHT_START) / 2 / 1e6, 1.1, 'Core', ha='center', fontsize=9, color='gray')
ax.text((ARM_RIGHT_START + CHR_LEN) / 2 / 1e6, 1.1, 'Right Arm', ha='center', fontsize=9, color='gray')
ax.set_xlim(-0.1, CHR_LEN / 1e6 + 0.1)
ax.set_ylim(-0.5, 1.3)
ax.set_xlabel('Chromosome Position (Mb)')
ax.set_yticks([0, 0.5])
ax.set_yticklabels(['Shielded', 'Exposed'])
ax.set_title('H28: Chromosome Map - 57 Exposed Regulators', fontweight='bold')
ax.legend(loc='upper right', fontsize=9)
plt.tight_layout()
fig.savefig(f'{FIG_DIR}/chromosome_map.pdf')
fig.savefig(f'{FIG_DIR}/chromosome_map.svg')
print("Saved: chromosome_map.pdf/svg")
plt.close()

# Panel E: Nearest methylation distance vs baseMean scatter
fig, ax = plt.subplots(figsize=(10, 8))
has_data = master['nearest_methyl_dist'].notna() & master['baseMean'].notna()
plot_data = master[has_data]
shi_p = plot_data[~plot_data['is_exposed']]
exp_p = plot_data[plot_data['is_exposed']]
ax.scatter(shi_p['nearest_methyl_dist'], shi_p['baseMean'],
           c=COLOR_SHIELDED, alpha=0.3, s=15, label='Shielded (998)', zorder=2)
ax.scatter(exp_p['nearest_methyl_dist'], exp_p['baseMean'],
           c=COLOR_EXPOSED, s=40, zorder=5, label='Exposed (57)', edgecolors='black', linewidth=0.5)

exp_nmd_med = master.loc[master['is_exposed'], 'nearest_methyl_dist'].dropna().median()
shi_nmd_med = master.loc[~master['is_exposed'], 'nearest_methyl_dist'].dropna().median()
ax.axvline(exp_nmd_med, color=COLOR_EXPOSED, linestyle='--', alpha=0.7, linewidth=1.5,
           label=f'Exposed median = {exp_nmd_med:.0f} bp')
ax.axvline(shi_nmd_med, color=COLOR_SHIELDED, linestyle='--', alpha=0.7, linewidth=1.5,
           label=f'Shielded median = {shi_nmd_med:.0f} bp')
ax.set_xlabel('Nearest Methylation Distance to TSS (bp)')
ax.set_ylabel('baseMean (expression)')
ax.set_yscale('log')
ax.set_title('H28: Methylation Proximity vs Expression Level\n(All 1,055 Regulatory Genes)', fontweight='bold')
ax.legend(fontsize=9, loc='upper right')
plt.tight_layout()
fig.savefig(f'{FIG_DIR}/methylation_distance_vs_expression.pdf')
fig.savefig(f'{FIG_DIR}/methylation_distance_vs_expression.svg')
print("Saved: methylation_distance_vs_expression.pdf/svg")
plt.close()

# Comprehensive Summary Figure
fig = plt.figure(figsize=(20, 16))
gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.3)

ax_a = fig.add_subplot(gs[0, 0])
make_violin_box(ax_a, [exp_t1v, shi_t1v], ['Exposed (57)', 'Shielded (998)'],
                'A. T1 Expression (Core Prediction)', 'T1 Normalized Counts', p_t1)

ax_b = fig.add_subplot(gs[0, 1])
make_violin_box(ax_b, [exp_bm, shi_bm], ['Exposed (57)', 'Shielded (998)'],
                'B. baseMean (Overall Expression)', 'baseMean', p_bm)

ax_c = fig.add_subplot(gs[0, 2])
make_violin_box(ax_c, [exp_abs2, shi_abs2], ['Exposed (57)', 'Shielded (998)'],
                'C. |log2FC| T2vsT1 (Change Magnitude)', '|log2FC|', p_abs2)

# D: scatter
ax_d = fig.add_subplot(gs[1, 0:2])
ax_d.scatter(shi_p['nearest_methyl_dist'], shi_p['baseMean'],
             c=COLOR_SHIELDED, alpha=0.3, s=12, label='Shielded', zorder=2)
ax_d.scatter(exp_p['nearest_methyl_dist'], exp_p['baseMean'],
             c=COLOR_EXPOSED, s=30, zorder=5, label='Exposed', edgecolors='black', linewidth=0.5)
ax_d.axvline(exp_nmd_med, color=COLOR_EXPOSED, linestyle='--', alpha=0.7, linewidth=1.5)
ax_d.axvline(shi_nmd_med, color=COLOR_SHIELDED, linestyle='--', alpha=0.7, linewidth=1.5)
ax_d.set_xlabel('Nearest Methylation Distance to TSS (bp)')
ax_d.set_ylabel('baseMean (log scale)')
ax_d.set_yscale('log')
ax_d.set_title('D. Methylation Proximity vs Expression', fontweight='bold')
ax_d.legend(fontsize=8)

# E: TF family bar
ax_e = fig.add_subplot(gs[1, 2])
top_fam = fisher_df.nlargest(10, 'n_total')
y_e = np.arange(len(top_fam))
colors_e = [COLOR_EXPOSED if r > 1 else COLOR_SHIELDED for r in top_fam['odds_ratio']]
ax_e.barh(y_e, top_fam['pct_exposed'], color=colors_e, edgecolor='white')
ax_e.axvline(57 / 1055 * 100, color='black', linestyle='--', linewidth=1)
ax_e.set_yticks(y_e)
ax_e.set_yticklabels(top_fam['tf_family'], fontsize=8)
ax_e.set_xlabel('% Exposed')
ax_e.set_title('E. Top TF Families: % Exposed', fontweight='bold')

# F: chromosome map
ax_f = fig.add_subplot(gs[2, :])
ax_f.scatter(shi_pos / 1e6, np.zeros(len(shi_pos)), c=COLOR_SHIELDED, alpha=0.3, s=8, label='Shielded')
ax_f.scatter(exp_pos / 1e6, np.ones(len(exp_pos)) * 0.5, c=COLOR_EXPOSED, s=25, zorder=5,
             label='Exposed', edgecolors='black', linewidth=0.5)
ax_f.axvline(ARM_LEFT_END / 1e6, color='gray', linestyle='--', alpha=0.5)
ax_f.axvline(ARM_RIGHT_START / 1e6, color='gray', linestyle='--', alpha=0.5)
ax_f.fill_betweenx([-0.5, 1.2], 0, ARM_LEFT_END / 1e6, color='lightyellow', alpha=0.3)
ax_f.fill_betweenx([-0.5, 1.2], ARM_RIGHT_START / 1e6, CHR_LEN / 1e6, color='lightyellow', alpha=0.3)
ax_f.set_xlim(-0.1, CHR_LEN / 1e6 + 0.1)
ax_f.set_ylim(-0.5, 1.2)
ax_f.set_xlabel('Chromosome Position (Mb)')
ax_f.set_yticks([0, 0.5])
ax_f.set_yticklabels(['Shielded', 'Exposed'])
ax_f.set_title('F. Chromosome Map', fontweight='bold')
ax_f.legend(fontsize=8, loc='upper right')

fig.suptitle('H28: Comprehensive Characterization of 57 Exposed Regulators',
             fontweight='bold', fontsize=14, y=1.01)
fig.savefig(f'{FIG_DIR}/H28_comprehensive_summary.pdf')
fig.savefig(f'{FIG_DIR}/H28_comprehensive_summary.svg')
print("Saved: H28_comprehensive_summary.pdf/svg")
plt.close()

# ============================================================
# FINAL SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("H28 RESULTS SUMMARY")
print("=" * 70)

n_sig = stats_df[stats_df['significant'] != 'ns'].shape[0]
print(f"\nTotal features tested: {len(stats_df)}")
print(f"Significant (p < 0.05): {n_sig}")

print("\n--- PRIMARY HYPOTHESIS (T1 expression) ---")
row_t1 = stats_df[stats_df['feature'] == 'T1_expression'].iloc[0]
if row_t1['p_value'] < 0.05:
    direction = "LOWER" if row_t1['exposed_median'] < row_t1['shielded_median'] else "HIGHER"
    print(f"  Exposed T1 median={row_t1['exposed_median']:.1f}, Shielded T1 median={row_t1['shielded_median']:.1f}")
    print(f"  p={row_t1['p_value']:.2e}, r={row_t1['effect_size_r']:.3f}")
    if direction == "LOWER":
        print(f"  >>> HYPOTHESIS SUPPORTED: Low T1 expression explains absent protection zone")
    else:
        print(f"  >>> HYPOTHESIS CONTRADICTED: Higher T1 expression, not lower")
else:
    print(f"  Exposed T1 median={row_t1['exposed_median']:.1f}, Shielded T1 median={row_t1['shielded_median']:.1f}")
    print(f"  p={row_t1['p_value']:.2e}")
    print(f"  >>> HYPOTHESIS NOT SUPPORTED")

print("\n--- KEY SIGNIFICANT FEATURES ---")
for _, row in stats_df[stats_df['significant'] != 'ns'].iterrows():
    r_str = f", r={row['effect_size_r']:.3f}" if not pd.isna(row['effect_size_r']) else ""
    print(f"  {row['feature']}: {row['direction']} (p={row['p_value']:.1e}{r_str})")

print("\n--- NON-SIGNIFICANT FEATURES ---")
for _, row in stats_df[stats_df['significant'] == 'ns'].iterrows():
    print(f"  {row['feature']}: p={row['p_value']:.2e}")

print(f"\n{'='*70}")
print("All outputs saved to:")
print(f"  Tables: {TAB_DIR}/")
print(f"  Figures: {FIG_DIR}/")
print(f"{'='*70}")
