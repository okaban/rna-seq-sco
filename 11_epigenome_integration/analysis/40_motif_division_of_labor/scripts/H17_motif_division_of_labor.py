#!/usr/bin/env python3
"""
H17: Motif-Specific "Division of Labor" in Functional Gene Targeting

Tests whether GCCGGC-proximal and AAGCCCG-proximal genes show distinct
expression dynamics: GCCGGC targets (Transport/SecMet enriched) undergo
T2/T3 developmental upregulation, while AAGCCCG targets (Stress/Defense
enriched) maintain constitutive or T1-biased expression.

The motif-specific targeting should have measurable transcriptomic consequences.
"""

import pandas as pd
import numpy as np
from scipy import stats
from collections import OrderedDict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
import warnings
warnings.filterwarnings('ignore')
import os
import urllib.parse

# ============================================================
# Configuration
# ============================================================
BASE_DIR = "/Users/okaban/bioinfo/rna-seq"
EPIGENOME_DIR = f"{BASE_DIR}/11_epigenome_integration"
OUT_DIR = f"{EPIGENOME_DIR}/analysis/40_motif_division_of_labor"
FIG_DIR = f"{OUT_DIR}/figures"
TBL_DIR = f"{OUT_DIR}/tables"

# Input files
AAGCCCG_MAPPING = f"{EPIGENOME_DIR}/analysis/36_AAGCCCG_distribution/tables/AAGCCCG_site_gene_mapping.tsv"
GCCGGC_SITES = f"{EPIGENOME_DIR}/analysis/37_defense_island_GCCGGC/tables/GCCGGC_sites_by_timepoint.tsv"
GFF_FILE = f"{EPIGENOME_DIR}/analysis/39_GCCGGC_MTase_reverse_ID/data/GCF_000203835.1_ASM20383v1_genomic.gff"
H15_ENRICHMENT = f"{EPIGENOME_DIR}/analysis/38_cross_motif_regulatory_avoidance/tables/cross_motif_functional_enrichment.tsv"
DESEQ2_DIR = f"{BASE_DIR}/04_deseq2/analysis/04_deseq2_260128_v1/results"
DESEQ2_T2vsT1 = f"{DESEQ2_DIR}/DESeq2_M145_2_vs_1.tsv"
DESEQ2_T3vsT1 = f"{DESEQ2_DIR}/DESeq2_M145_3_vs_1.tsv"
NORM_COUNTS = f"{DESEQ2_DIR}/normalized_counts_M145.tsv"

WINDOW = 2000  # 2kb proximity window

np.random.seed(42)

# Color palette
COLORS = {
    'GCCGGC': '#2166AC',      # deep blue
    'AAGCCCG': '#B2182B',     # deep red
    'Dual': '#762A83',        # purple
    'Background': '#969696',  # grey
}

CATEGORY_COLORS = {
    'Regulatory/TF': '#E41A1C',
    'Transport': '#377EB8',
    'Secondary metabolism': '#4DAF4A',
    'Stress/Defense': '#FF7F00',
    'Primary metabolism': '#984EA3',
    'DNA/RNA metabolism': '#A65628',
    'Translation': '#F781BF',
    'Membrane/Cell wall': '#999999',
    'Hypothetical': '#CCCCCC',
    'Other': '#666666',
}

# ============================================================
# 1. Load GFF and build gene annotation database
# ============================================================
print("=" * 70)
print("H17: Motif-Specific Division of Labor Analysis")
print("=" * 70)

print("\n[1] Loading gene annotations from GFF...")

genes = []
with open(GFF_FILE) as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue
        if fields[2] != 'gene':
            continue
        chrom = fields[0]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]
        attrs = fields[8]
        attr_dict = {}
        for attr in attrs.split(';'):
            if '=' in attr:
                k, v = attr.split('=', 1)
                attr_dict[k] = v
        locus_tag = attr_dict.get('locus_tag', attr_dict.get('Name', ''))
        old_locus_tag = attr_dict.get('old_locus_tag', '')
        gene_name = attr_dict.get('gene', '')
        gene_biotype = attr_dict.get('gene_biotype', '')
        genes.append({
            'chrom': chrom,
            'start': start,
            'end': end,
            'strand': strand,
            'locus_tag': locus_tag,
            'old_locus_tag': old_locus_tag,
            'gene_name': gene_name,
            'gene_biotype': gene_biotype
        })

genes_df = pd.DataFrame(genes)
print(f"  Total genes loaded: {len(genes_df)}")

# Load CDS product annotations
cds_info = {}
with open(GFF_FILE) as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue
        if fields[2] != 'CDS':
            continue
        attrs = fields[8]
        attr_dict = {}
        for attr in attrs.split(';'):
            if '=' in attr:
                k, v = attr.split('=', 1)
                attr_dict[k] = v
        lt = attr_dict.get('locus_tag', '')
        product = attr_dict.get('product', 'hypothetical protein')
        # URL-decode product names
        product = urllib.parse.unquote(product)
        if lt and lt not in cds_info:
            cds_info[lt] = product
genes_df['product'] = genes_df['locus_tag'].map(cds_info).fillna('unknown')
protein_coding = genes_df[genes_df['gene_biotype'] == 'protein_coding'].copy()
print(f"  Protein-coding genes: {len(protein_coding)}")

# ============================================================
# 2. Functional category classification (same as H15)
# ============================================================
def classify_product(product):
    """Classify gene product into functional categories (identical to H13/H15)."""
    product_lower = product.lower()
    if any(w in product_lower for w in ['transcriptional regulator', 'transcription factor',
                                          'sigma factor', 'response regulator', 'sensor kinase',
                                          'two-component', 'anti-sigma', 'repressor',
                                          'dna-binding', 'regulatory']):
        return 'Regulatory/TF'
    if any(w in product_lower for w in ['transporter', 'permease', 'abc transport', 'mfs ',
                                          'efflux', 'porin', 'channel']):
        return 'Transport'
    if any(w in product_lower for w in ['polyketide', 'non-ribosomal', 'nrps', 'synthase',
                                          'ketoacyl', 'acyl carrier', 'thioesterase',
                                          'actinorhodin', 'undecylprodigiosin',
                                          'antibiotic biosynthesis']):
        return 'Secondary metabolism'
    if any(w in product_lower for w in ['stress', 'heat shock', 'cold shock', 'chaperon',
                                          'protease', 'peptidase', 'dnaj', 'dnak', 'grpe',
                                          'catalase', 'superoxide', 'thioredoxin',
                                          'glutaredoxin', 'universal stress']):
        return 'Stress/Defense'
    if any(w in product_lower for w in ['methyltransferase', 'helicase', 'dnase', 'rnase',
                                          'topoisomerase', 'gyrase', 'recombinase',
                                          'restriction', 'ligase', 'polymerase',
                                          'nuclease', 'dna repair']):
        return 'DNA/RNA metabolism'
    if any(w in product_lower for w in ['dehydrogenase', 'oxidoreductase', 'kinase',
                                          'transferase', 'hydrolase', 'lyase', 'isomerase',
                                          'reductase', 'oxidase', 'synthetase',
                                          'phosphatase', 'esterase']):
        return 'Primary metabolism'
    if any(w in product_lower for w in ['membrane', 'cell wall', 'peptidoglycan',
                                          'lipopolysaccharide', 'penicillin-binding',
                                          'murein']):
        return 'Membrane/Cell wall'
    if any(w in product_lower for w in ['ribosom', 'trna', 'rrna', 'translation',
                                          'elongation factor', 'initiation factor']):
        return 'Translation'
    if 'hypothetical' in product_lower or product_lower == 'unknown':
        return 'Hypothetical'
    return 'Other'

protein_coding = protein_coding.copy()
protein_coding['category'] = protein_coding['product'].apply(classify_product)

cat_counts = protein_coding['category'].value_counts()
print("\n  Genome-wide functional category distribution:")
for cat, n in cat_counts.items():
    print(f"    {cat}: {n} ({100*n/len(protein_coding):.1f}%)")

# ============================================================
# 3. Build motif-proximal gene sets
# ============================================================
print("\n[2] Building motif-proximal gene sets...")

# --- AAGCCCG-proximal genes (from pre-computed mapping) ---
aagcccg_mapping = pd.read_csv(AAGCCCG_MAPPING, sep='\t')
# Filter to T1 and within 2kb
aagcccg_t1 = aagcccg_mapping[(aagcccg_mapping['timepoint'] == 'T1') &
                               (aagcccg_mapping['within_2kb'] == True)]
aagcccg_genes = set(aagcccg_t1['locus_tag'].unique())
print(f"  AAGCCCG-proximal genes (T1, within 2kb): {len(aagcccg_genes)}")
print(f"    From {aagcccg_t1['position'].nunique()} unique AAGCCCG sites")

# --- GCCGGC-proximal genes (map sites to nearby genes) ---
gccggc_sites = pd.read_csv(GCCGGC_SITES, sep='\t')
gccggc_t1 = gccggc_sites[gccggc_sites['timepoint'] == 'T1'].copy()
print(f"  GCCGGC T1 sites: {len(gccggc_t1)}")

# Map GCCGGC sites to nearby genes within 2kb
gccggc_gene_set = set()
gccggc_site_gene_pairs = []

for _, site in gccggc_t1.iterrows():
    site_pos = site['position']
    # Find genes within 2kb window
    nearby = protein_coding[
        (protein_coding['chrom'] == 'NC_003888.3') &
        (
            # Gene overlaps the 2kb window around the site
            (protein_coding['start'] <= site_pos + WINDOW) &
            (protein_coding['end'] >= site_pos - WINDOW)
        )
    ]
    # Also check distances for genes that don't overlap but are nearby
    for _, gene in protein_coding[protein_coding['chrom'] == 'NC_003888.3'].iterrows():
        dist = min(abs(gene['start'] - site_pos), abs(gene['end'] - site_pos))
        if site_pos >= gene['start'] and site_pos <= gene['end']:
            dist = 0
        if dist <= WINDOW:
            gccggc_gene_set.add(gene['locus_tag'])
            gccggc_site_gene_pairs.append({
                'site_position': site_pos,
                'locus_tag': gene['locus_tag'],
                'distance': dist
            })

gccggc_genes = gccggc_gene_set
print(f"  GCCGGC-proximal genes (T1, within 2kb): {len(gccggc_genes)}")

# Overlap
overlap_genes = aagcccg_genes & gccggc_genes
aagcccg_only = aagcccg_genes - gccggc_genes
gccggc_only = gccggc_genes - aagcccg_genes
all_targeted = aagcccg_genes | gccggc_genes

# Background: all protein-coding genes NOT proximal to either motif
background_genes = set(protein_coding['locus_tag']) - all_targeted

print(f"\n  Overlap analysis:")
print(f"    AAGCCCG-only: {len(aagcccg_only)}")
print(f"    GCCGGC-only: {len(gccggc_only)}")
print(f"    Dual-targeted (both motifs): {len(overlap_genes)}")
print(f"    Background (neither): {len(background_genes)}")
jaccard = len(overlap_genes) / len(all_targeted) if len(all_targeted) > 0 else 0
print(f"    Jaccard index: {jaccard:.3f}")

# ============================================================
# 4. Load DESeq2 expression data
# ============================================================
print("\n[3] Loading DESeq2 expression data...")

deseq2_t2 = pd.read_csv(DESEQ2_T2vsT1, sep='\t')
deseq2_t3 = pd.read_csv(DESEQ2_T3vsT1, sep='\t')
norm_counts = pd.read_csv(NORM_COUNTS, sep='\t')

print(f"  T2vsT1: {len(deseq2_t2)} genes")
print(f"  T3vsT1: {len(deseq2_t3)} genes")
print(f"  Normalized counts: {len(norm_counts)} genes x {norm_counts.shape[1]-1} samples")

# Merge T2 and T3 data
expr = deseq2_t2[['gene_id', 'baseMean', 'log2FoldChange', 'padj']].rename(
    columns={'log2FoldChange': 'LFC_T2vsT1', 'padj': 'padj_T2vsT1'})
expr = expr.merge(
    deseq2_t3[['gene_id', 'log2FoldChange', 'padj']].rename(
        columns={'log2FoldChange': 'LFC_T3vsT1', 'padj': 'padj_T3vsT1'}),
    on='gene_id', how='outer'
)

# Add category information
gene_cat = protein_coding[['locus_tag', 'category', 'product']].rename(columns={'locus_tag': 'gene_id'})
expr = expr.merge(gene_cat, on='gene_id', how='left')

# Assign motif proximity groups
def assign_group(gene_id):
    if gene_id in overlap_genes:
        return 'Dual'
    elif gene_id in aagcccg_only:
        return 'AAGCCCG-only'
    elif gene_id in gccggc_only:
        return 'GCCGGC-only'
    else:
        return 'Background'

expr['motif_group'] = expr['gene_id'].apply(assign_group)

# Also create broader groups for primary comparisons
def assign_broad_group(gene_id):
    if gene_id in aagcccg_genes:
        return 'AAGCCCG-proximal'
    elif gene_id in gccggc_genes:
        return 'GCCGGC-proximal'
    else:
        return 'Background'

# For non-overlapping comparison
def assign_exclusive_group(gene_id):
    if gene_id in overlap_genes:
        return 'Dual-targeted'
    elif gene_id in aagcccg_only:
        return 'AAGCCCG-only'
    elif gene_id in gccggc_only:
        return 'GCCGGC-only'
    else:
        return 'Background'

expr['exclusive_group'] = expr['gene_id'].apply(assign_exclusive_group)

# Filter to genes with valid expression data
expr_valid = expr.dropna(subset=['LFC_T2vsT1', 'LFC_T3vsT1']).copy()
print(f"\n  Genes with valid LFC in both comparisons: {len(expr_valid)}")

group_counts = expr_valid['exclusive_group'].value_counts()
for g, n in group_counts.items():
    print(f"    {g}: {n}")

# ============================================================
# 5. Overall LFC distribution comparison
# ============================================================
print("\n[4] Comparing overall LFC distributions...")

summary_stats = []
for group in ['AAGCCCG-only', 'GCCGGC-only', 'Dual-targeted', 'Background']:
    sub = expr_valid[expr_valid['exclusive_group'] == group]
    if len(sub) == 0:
        continue

    # DEG rates (padj < 0.05 and |LFC| > 1)
    deg_t2 = sub[(sub['padj_T2vsT1'] < 0.05) & (sub['LFC_T2vsT1'].abs() > 1)]
    deg_t3 = sub[(sub['padj_T3vsT1'] < 0.05) & (sub['LFC_T3vsT1'].abs() > 1)]

    # Direction
    up_t2 = len(sub[(sub['padj_T2vsT1'] < 0.05) & (sub['LFC_T2vsT1'] > 1)])
    down_t2 = len(sub[(sub['padj_T2vsT1'] < 0.05) & (sub['LFC_T2vsT1'] < -1)])
    up_t3 = len(sub[(sub['padj_T3vsT1'] < 0.05) & (sub['LFC_T3vsT1'] > 1)])
    down_t3 = len(sub[(sub['padj_T3vsT1'] < 0.05) & (sub['LFC_T3vsT1'] < -1)])

    stats_row = {
        'group': group,
        'n_genes': len(sub),
        'median_LFC_T2vsT1': sub['LFC_T2vsT1'].median(),
        'mean_LFC_T2vsT1': sub['LFC_T2vsT1'].mean(),
        'std_LFC_T2vsT1': sub['LFC_T2vsT1'].std(),
        'median_LFC_T3vsT1': sub['LFC_T3vsT1'].median(),
        'mean_LFC_T3vsT1': sub['LFC_T3vsT1'].mean(),
        'std_LFC_T3vsT1': sub['LFC_T3vsT1'].std(),
        'n_DEG_T2vsT1': len(deg_t2),
        'DEG_rate_T2vsT1': len(deg_t2) / len(sub),
        'n_DEG_T3vsT1': len(deg_t3),
        'DEG_rate_T3vsT1': len(deg_t3) / len(sub),
        'n_up_T2': up_t2,
        'n_down_T2': down_t2,
        'up_ratio_T2': up_t2 / max(up_t2 + down_t2, 1),
        'n_up_T3': up_t3,
        'n_down_T3': down_t3,
        'up_ratio_T3': up_t3 / max(up_t3 + down_t3, 1),
    }
    summary_stats.append(stats_row)

    print(f"\n  {group} (n={len(sub)}):")
    print(f"    T2vsT1: median LFC={sub['LFC_T2vsT1'].median():.3f}, "
          f"mean={sub['LFC_T2vsT1'].mean():.3f}, "
          f"DEG rate={100*len(deg_t2)/len(sub):.1f}% "
          f"(up={up_t2}, down={down_t2})")
    print(f"    T3vsT1: median LFC={sub['LFC_T3vsT1'].median():.3f}, "
          f"mean={sub['LFC_T3vsT1'].mean():.3f}, "
          f"DEG rate={100*len(deg_t3)/len(sub):.1f}% "
          f"(up={up_t3}, down={down_t3})")

summary_df = pd.DataFrame(summary_stats)

# ============================================================
# 6. Statistical tests between motif groups
# ============================================================
print("\n[5] Statistical tests between motif groups...")

stat_tests = []
groups = ['AAGCCCG-only', 'GCCGGC-only', 'Dual-targeted']
comparisons = [
    ('AAGCCCG-only', 'GCCGGC-only'),
    ('AAGCCCG-only', 'Background'),
    ('GCCGGC-only', 'Background'),
    ('Dual-targeted', 'Background'),
    ('AAGCCCG-only', 'Dual-targeted'),
    ('GCCGGC-only', 'Dual-targeted'),
]

for g1, g2 in comparisons:
    sub1 = expr_valid[expr_valid['exclusive_group'] == g1]
    sub2 = expr_valid[expr_valid['exclusive_group'] == g2]
    if len(sub1) < 5 or len(sub2) < 5:
        continue

    for metric in ['LFC_T2vsT1', 'LFC_T3vsT1']:
        u_stat, p_val = stats.mannwhitneyu(sub1[metric], sub2[metric], alternative='two-sided')
        # Effect size (rank-biserial correlation)
        n1, n2 = len(sub1), len(sub2)
        r = 1 - (2*u_stat)/(n1*n2)

        stat_tests.append({
            'comparison': f"{g1} vs {g2}",
            'metric': metric,
            'n1': n1,
            'n2': n2,
            'median1': sub1[metric].median(),
            'median2': sub2[metric].median(),
            'U_statistic': u_stat,
            'p_value': p_val,
            'rank_biserial_r': r,
        })

        sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else "n.s."
        print(f"  {g1} vs {g2} [{metric}]: U={u_stat:.0f}, p={p_val:.4e}, r={r:.3f} {sig}")

stat_tests_df = pd.DataFrame(stat_tests)

# Bonferroni correction
n_tests = len(stat_tests_df)
stat_tests_df['p_bonferroni'] = (stat_tests_df['p_value'] * n_tests).clip(upper=1.0)
stat_tests_df['significant_bonferroni'] = stat_tests_df['p_bonferroni'] < 0.05

# ============================================================
# 7. Functional category breakdown by motif proximity
# ============================================================
print("\n[6] Functional category breakdown by motif proximity...")

category_results = []
categories = ['Regulatory/TF', 'Transport', 'Secondary metabolism', 'Stress/Defense',
              'Primary metabolism', 'DNA/RNA metabolism', 'Translation', 'Membrane/Cell wall',
              'Hypothetical', 'Other']

for cat in categories:
    for group in ['AAGCCCG-only', 'GCCGGC-only', 'Dual-targeted', 'Background']:
        sub = expr_valid[(expr_valid['exclusive_group'] == group) & (expr_valid['category'] == cat)]
        if len(sub) == 0:
            category_results.append({
                'category': cat,
                'motif_group': group,
                'n_genes': 0,
                'median_LFC_T2vsT1': np.nan,
                'median_LFC_T3vsT1': np.nan,
                'DEG_rate_T2': np.nan,
                'DEG_rate_T3': np.nan,
            })
            continue

        deg_t2 = sub[(sub['padj_T2vsT1'] < 0.05) & (sub['LFC_T2vsT1'].abs() > 1)]
        deg_t3 = sub[(sub['padj_T3vsT1'] < 0.05) & (sub['LFC_T3vsT1'].abs() > 1)]

        category_results.append({
            'category': cat,
            'motif_group': group,
            'n_genes': len(sub),
            'median_LFC_T2vsT1': sub['LFC_T2vsT1'].median(),
            'mean_LFC_T2vsT1': sub['LFC_T2vsT1'].mean(),
            'median_LFC_T3vsT1': sub['LFC_T3vsT1'].median(),
            'mean_LFC_T3vsT1': sub['LFC_T3vsT1'].mean(),
            'DEG_rate_T2': len(deg_t2) / len(sub) if len(sub) > 0 else np.nan,
            'DEG_rate_T3': len(deg_t3) / len(sub) if len(sub) > 0 else np.nan,
        })

cat_df = pd.DataFrame(category_results)

# Specific hypothesis tests from H15:
# 1. Transport genes near GCCGGC vs all Transport genes
# 2. Stress/Defense genes near AAGCCCG vs all Stress/Defense genes
print("\n  Specific functional category tests:")

specific_tests = []
test_pairs = [
    ('Transport', 'GCCGGC-only', 'Background',
     'Transport genes near GCCGGC vs Background Transport genes'),
    ('Stress/Defense', 'AAGCCCG-only', 'Background',
     'Stress/Defense genes near AAGCCCG vs Background Stress/Defense genes'),
    ('Secondary metabolism', 'GCCGGC-only', 'Background',
     'SecMet genes near GCCGGC vs Background SecMet genes'),
    ('Primary metabolism', 'GCCGGC-only', 'Background',
     'Primary metab genes near GCCGGC vs Background Primary metab genes'),
    ('DNA/RNA metabolism', 'GCCGGC-only', 'Background',
     'DNA/RNA metab genes near GCCGGC vs Background DNA/RNA metab genes'),
]

for cat, g1, g2, desc in test_pairs:
    sub1 = expr_valid[(expr_valid['exclusive_group'] == g1) & (expr_valid['category'] == cat)]
    sub2 = expr_valid[(expr_valid['exclusive_group'] == g2) & (expr_valid['category'] == cat)]

    for metric in ['LFC_T2vsT1', 'LFC_T3vsT1']:
        if len(sub1) >= 3 and len(sub2) >= 3:
            u_stat, p_val = stats.mannwhitneyu(sub1[metric], sub2[metric], alternative='two-sided')
            n1, n2 = len(sub1), len(sub2)
            r = 1 - (2*u_stat)/(n1*n2)
            specific_tests.append({
                'test': desc,
                'metric': metric,
                'n_motif_proximal': n1,
                'n_background': n2,
                'median_proximal': sub1[metric].median(),
                'median_background': sub2[metric].median(),
                'U_statistic': u_stat,
                'p_value': p_val,
                'rank_biserial_r': r,
            })
            sig = "*" if p_val < 0.05 else "n.s."
            print(f"    {desc} [{metric}]: n={n1} vs {n2}, "
                  f"median={sub1[metric].median():.3f} vs {sub2[metric].median():.3f}, "
                  f"p={p_val:.4e} {sig}")
        else:
            print(f"    {desc} [{metric}]: insufficient data (n={len(sub1)} vs {len(sub2)})")

specific_tests_df = pd.DataFrame(specific_tests)

# ============================================================
# 8. Temporal expression classification
# ============================================================
print("\n[7] Temporal expression pattern classification...")

def classify_temporal(row):
    """Classify temporal expression pattern based on T2vsT1 and T3vsT1 LFC."""
    lfc2 = row['LFC_T2vsT1']
    lfc3 = row['LFC_T3vsT1']
    padj2 = row['padj_T2vsT1']
    padj3 = row['padj_T3vsT1']

    sig2 = (padj2 < 0.05) and (abs(lfc2) > 1) if pd.notna(padj2) else False
    sig3 = (padj3 < 0.05) and (abs(lfc3) > 1) if pd.notna(padj3) else False

    if not sig2 and not sig3:
        return 'Constitutive'
    elif sig2 and not sig3:
        if lfc2 > 0:
            return 'T2-up (transient)'
        else:
            return 'T2-down (transient)'
    elif not sig2 and sig3:
        if lfc3 > 0:
            return 'T3-up (late)'
        else:
            return 'T3-down (late)'
    elif sig2 and sig3:
        if lfc2 > 0 and lfc3 > 0:
            if lfc3 > lfc2:
                return 'Monotonic increasing'
            else:
                return 'Sustained up'
        elif lfc2 < 0 and lfc3 < 0:
            if lfc3 < lfc2:
                return 'Monotonic decreasing'
            else:
                return 'Sustained down'
        elif lfc2 > 0 and lfc3 < 0:
            return 'Peak at T2'
        elif lfc2 < 0 and lfc3 > 0:
            return 'Trough at T2'
    return 'Other'

expr_valid = expr_valid.copy()
expr_valid['temporal_class'] = expr_valid.apply(classify_temporal, axis=1)

print("\n  Temporal classification distribution by motif group:")
temporal_summary = pd.crosstab(expr_valid['exclusive_group'], expr_valid['temporal_class'])
print(temporal_summary.to_string())

# Percentage version
temporal_pct = pd.crosstab(expr_valid['exclusive_group'], expr_valid['temporal_class'], normalize='index') * 100
print("\n  Percentages:")
print(temporal_pct.round(1).to_string())

# Chi-squared test for temporal class distribution
print("\n  Chi-squared tests for temporal class distribution:")
for g1, g2 in [('AAGCCCG-only', 'GCCGGC-only'), ('AAGCCCG-only', 'Background'), ('GCCGGC-only', 'Background')]:
    sub = expr_valid[expr_valid['exclusive_group'].isin([g1, g2])]
    ct = pd.crosstab(sub['exclusive_group'], sub['temporal_class'])
    if ct.shape[0] == 2 and ct.shape[1] > 1:
        chi2, p, dof, expected = stats.chi2_contingency(ct)
        sig = "*" if p < 0.05 else "n.s."
        print(f"    {g1} vs {g2}: chi2={chi2:.2f}, dof={dof}, p={p:.4e} {sig}")

# ============================================================
# 9. Dual-targeted gene analysis
# ============================================================
print("\n[8] Dual-targeted gene analysis...")

dual_genes = expr_valid[expr_valid['exclusive_group'] == 'Dual-targeted'].copy()
print(f"  Dual-targeted genes with expression data: {len(dual_genes)}")

if len(dual_genes) > 0:
    print(f"\n  Top dual-targeted genes by absolute LFC_T3vsT1:")
    top_dual = dual_genes.nlargest(min(20, len(dual_genes)), 'LFC_T3vsT1', keep='first')
    for _, row in top_dual.iterrows():
        print(f"    {row['gene_id']}: LFC_T2={row['LFC_T2vsT1']:.2f}, "
              f"LFC_T3={row['LFC_T3vsT1']:.2f}, {row['category']}, {str(row.get('product',''))[:50]}")

# Save dual-targeted genes
dual_out = dual_genes[['gene_id', 'baseMean', 'LFC_T2vsT1', 'padj_T2vsT1',
                        'LFC_T3vsT1', 'padj_T3vsT1', 'category', 'product',
                        'temporal_class']].copy()
dual_out.to_csv(f"{TBL_DIR}/dual_targeted_genes.tsv", sep='\t', index=False)
print(f"  Saved: {TBL_DIR}/dual_targeted_genes.tsv")

# ============================================================
# 10. Save summary tables
# ============================================================
print("\n[9] Saving output tables...")

# Expression summary stats
summary_df.to_csv(f"{TBL_DIR}/expression_summary_stats.tsv", sep='\t', index=False)
print(f"  Saved: expression_summary_stats.tsv")

# Gene set comparison (full table with all genes)
gene_set_out = expr_valid[['gene_id', 'baseMean', 'LFC_T2vsT1', 'padj_T2vsT1',
                            'LFC_T3vsT1', 'padj_T3vsT1', 'category', 'exclusive_group',
                            'temporal_class']].copy()
gene_set_out.to_csv(f"{TBL_DIR}/gene_set_comparison.tsv", sep='\t', index=False)
print(f"  Saved: gene_set_comparison.tsv")

# Functional category by motif
cat_df.to_csv(f"{TBL_DIR}/functional_category_by_motif.tsv", sep='\t', index=False)
print(f"  Saved: functional_category_by_motif.tsv")

# Statistical tests
if len(stat_tests_df) > 0:
    stat_tests_df.to_csv(f"{TBL_DIR}/statistical_tests.tsv", sep='\t', index=False)
    print(f"  Saved: statistical_tests.tsv")

if len(specific_tests_df) > 0:
    specific_tests_df.to_csv(f"{TBL_DIR}/specific_category_tests.tsv", sep='\t', index=False)
    print(f"  Saved: specific_category_tests.tsv")

# ============================================================
# 11. FIGURE 1: Violin plots of LFC distributions
# ============================================================
print("\n[10] Generating figures...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

groups_to_plot = ['AAGCCCG-only', 'GCCGGC-only', 'Dual-targeted', 'Background']
group_colors = [COLORS['AAGCCCG'], COLORS['GCCGGC'], COLORS['Dual'], COLORS['Background']]

for ax_idx, (metric, title) in enumerate([('LFC_T2vsT1', 'T2 vs T1'), ('LFC_T3vsT1', 'T3 vs T1')]):
    ax = axes[ax_idx]
    data_lists = []
    labels = []
    for g in groups_to_plot:
        sub = expr_valid[expr_valid['exclusive_group'] == g][metric].dropna()
        if len(sub) > 0:
            data_lists.append(sub.values)
            labels.append(f"{g}\n(n={len(sub)})")
        else:
            data_lists.append(np.array([0]))
            labels.append(f"{g}\n(n=0)")

    parts = ax.violinplot(data_lists, positions=range(len(groups_to_plot)),
                          showmedians=True, showextrema=False)

    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(group_colors[i])
        pc.set_alpha(0.6)

    parts['cmedians'].set_color('black')

    # Add box plots inside
    bp = ax.boxplot(data_lists, positions=range(len(groups_to_plot)),
                    widths=0.15, showfliers=False, patch_artist=True,
                    medianprops=dict(color='black', linewidth=2),
                    boxprops=dict(linewidth=0.5),
                    whiskerprops=dict(linewidth=0.5),
                    capprops=dict(linewidth=0.5))
    for i, patch in enumerate(bp['boxes']):
        patch.set_facecolor(group_colors[i])
        patch.set_alpha(0.8)

    ax.set_xticks(range(len(groups_to_plot)))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel('log2 Fold Change', fontsize=11)
    ax.set_title(f'LFC Distribution: {title}', fontsize=12, fontweight='bold')
    ax.axhline(y=0, color='grey', linestyle='--', alpha=0.5)
    ax.set_ylim(-8, 8)

    # Add significance annotations
    # Get p-values for AAGCCCG vs GCCGGC comparison
    test_row = stat_tests_df[(stat_tests_df['comparison'] == 'AAGCCCG-only vs GCCGGC-only') &
                              (stat_tests_df['metric'] == metric)]
    if len(test_row) > 0:
        p = test_row.iloc[0]['p_value']
        if p < 0.001:
            sig_text = f'p={p:.2e}'
        elif p < 0.05:
            sig_text = f'p={p:.3f}'
        else:
            sig_text = f'p={p:.3f} n.s.'

        y_max = max(np.percentile(data_lists[0], 95), np.percentile(data_lists[1], 95)) + 1
        ax.plot([0, 0, 1, 1], [y_max, y_max+0.3, y_max+0.3, y_max], 'k-', linewidth=0.8)
        ax.text(0.5, y_max+0.4, sig_text, ha='center', fontsize=8)

plt.tight_layout()
plt.savefig(f"{FIG_DIR}/violin_LFC_by_motif.pdf", dpi=300, bbox_inches='tight')
plt.savefig(f"{FIG_DIR}/violin_LFC_by_motif.svg", dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: violin_LFC_by_motif.pdf/svg")

# ============================================================
# 12. FIGURE 2: Heatmap of functional category x motif proximity
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(16, 8))

for ax_idx, (metric, title) in enumerate([('median_LFC_T2vsT1', 'Median LFC T2vsT1'),
                                            ('median_LFC_T3vsT1', 'Median LFC T3vsT1')]):
    ax = axes[ax_idx]
    heatmap_data = cat_df.pivot_table(index='category', columns='motif_group',
                                       values=metric, aggfunc='first')
    # Reorder
    cat_order = [c for c in categories if c in heatmap_data.index]
    group_order = [g for g in groups_to_plot if g in heatmap_data.columns]
    heatmap_data = heatmap_data.reindex(index=cat_order, columns=group_order)

    # Also get gene counts
    count_data = cat_df.pivot_table(index='category', columns='motif_group',
                                     values='n_genes', aggfunc='first')
    count_data = count_data.reindex(index=cat_order, columns=group_order)

    im = ax.imshow(heatmap_data.values, cmap='RdBu_r', aspect='auto',
                   vmin=-1.5, vmax=1.5)

    ax.set_xticks(range(len(group_order)))
    ax.set_xticklabels(group_order, fontsize=9, rotation=45, ha='right')
    ax.set_yticks(range(len(cat_order)))
    ax.set_yticklabels(cat_order, fontsize=9)

    # Add text annotations
    for i in range(len(cat_order)):
        for j in range(len(group_order)):
            val = heatmap_data.values[i, j]
            n = count_data.values[i, j]
            if pd.notna(val) and pd.notna(n):
                color = 'white' if abs(val) > 0.8 else 'black'
                ax.text(j, i, f'{val:.2f}\n(n={int(n)})', ha='center', va='center',
                       fontsize=7, color=color)

    ax.set_title(title, fontsize=12, fontweight='bold')
    plt.colorbar(im, ax=ax, shrink=0.8, label='Median LFC')

plt.suptitle('Functional Category Expression by Motif Proximity', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(f"{FIG_DIR}/heatmap_category_motif.pdf", dpi=300, bbox_inches='tight')
plt.savefig(f"{FIG_DIR}/heatmap_category_motif.svg", dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: heatmap_category_motif.pdf/svg")

# ============================================================
# 13. FIGURE 3: Temporal class distribution
# ============================================================
fig, ax = plt.subplots(figsize=(14, 7))

temporal_classes = ['Constitutive', 'T2-up (transient)', 'T2-down (transient)',
                    'T3-up (late)', 'T3-down (late)', 'Monotonic increasing',
                    'Monotonic decreasing', 'Sustained up', 'Sustained down',
                    'Peak at T2', 'Trough at T2', 'Other']

temporal_colors = plt.cm.Set3(np.linspace(0, 1, len(temporal_classes)))

groups_for_temp = ['AAGCCCG-only', 'GCCGGC-only', 'Dual-targeted', 'Background']
x = np.arange(len(groups_for_temp))
width = 0.8

bottom = np.zeros(len(groups_for_temp))
for i, tc in enumerate(temporal_classes):
    vals = []
    for g in groups_for_temp:
        sub = expr_valid[expr_valid['exclusive_group'] == g]
        total = len(sub)
        n_class = len(sub[sub['temporal_class'] == tc])
        vals.append(100 * n_class / total if total > 0 else 0)
    vals = np.array(vals)
    ax.bar(x, vals, width, bottom=bottom, label=tc, color=temporal_colors[i], edgecolor='white', linewidth=0.5)
    bottom += vals

ax.set_xticks(x)
ax.set_xticklabels([f"{g}\n(n={len(expr_valid[expr_valid['exclusive_group']==g])})"
                    for g in groups_for_temp], fontsize=10)
ax.set_ylabel('Percentage of genes (%)', fontsize=11)
ax.set_title('Temporal Expression Pattern Distribution by Motif Proximity', fontsize=13, fontweight='bold')
ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8, ncol=1)
ax.set_ylim(0, 105)

plt.tight_layout()
plt.savefig(f"{FIG_DIR}/temporal_class_distribution.pdf", dpi=300, bbox_inches='tight')
plt.savefig(f"{FIG_DIR}/temporal_class_distribution.svg", dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: temporal_class_distribution.pdf/svg")

# ============================================================
# 14. FIGURE 4: Comprehensive 4-panel summary
# ============================================================
fig = plt.figure(figsize=(18, 16))
gs = gridspec.GridSpec(2, 2, hspace=0.35, wspace=0.35)

# Panel A: LFC comparison (AAGCCCG vs GCCGGC) for T3vsT1
ax_a = fig.add_subplot(gs[0, 0])
for group, color, marker, label in [
    ('AAGCCCG-only', COLORS['AAGCCCG'], 'o', 'AAGCCCG-only'),
    ('GCCGGC-only', COLORS['GCCGGC'], 's', 'GCCGGC-only'),
    ('Dual-targeted', COLORS['Dual'], 'D', 'Dual-targeted'),
]:
    sub = expr_valid[expr_valid['exclusive_group'] == group]
    ax_a.scatter(sub['LFC_T2vsT1'], sub['LFC_T3vsT1'],
                c=color, alpha=0.3, s=15, marker=marker, label=f"{label} (n={len(sub)})")

ax_a.axhline(0, color='grey', linestyle='--', alpha=0.5)
ax_a.axvline(0, color='grey', linestyle='--', alpha=0.5)
ax_a.plot([-8, 8], [-8, 8], 'k--', alpha=0.3)
ax_a.set_xlabel('LFC T2vsT1', fontsize=11)
ax_a.set_ylabel('LFC T3vsT1', fontsize=11)
ax_a.set_title('A) Expression Dynamics by Motif Proximity', fontsize=12, fontweight='bold')
ax_a.legend(fontsize=8, loc='upper left')
ax_a.set_xlim(-8, 8)
ax_a.set_ylim(-8, 8)

# Panel B: DEG rate comparison bar chart
ax_b = fig.add_subplot(gs[0, 1])
bar_groups = ['AAGCCCG-only', 'GCCGGC-only', 'Dual-targeted', 'Background']
bar_colors = [COLORS['AAGCCCG'], COLORS['GCCGGC'], COLORS['Dual'], COLORS['Background']]
x_pos = np.arange(len(bar_groups))
bar_width = 0.35

deg_rates_t2 = []
deg_rates_t3 = []
for g in bar_groups:
    row = summary_df[summary_df['group'] == g]
    if len(row) > 0:
        deg_rates_t2.append(100 * row.iloc[0]['DEG_rate_T2vsT1'])
        deg_rates_t3.append(100 * row.iloc[0]['DEG_rate_T3vsT1'])
    else:
        deg_rates_t2.append(0)
        deg_rates_t3.append(0)

bars1 = ax_b.bar(x_pos - bar_width/2, deg_rates_t2, bar_width,
                  label='T2vsT1', color=[c for c in bar_colors], alpha=0.6, edgecolor='black', linewidth=0.5)
bars2 = ax_b.bar(x_pos + bar_width/2, deg_rates_t3, bar_width,
                  label='T3vsT1', color=[c for c in bar_colors], alpha=1.0, edgecolor='black', linewidth=0.5)

# Add value labels
for bar, val in zip(bars1, deg_rates_t2):
    ax_b.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
             f'{val:.1f}%', ha='center', fontsize=7)
for bar, val in zip(bars2, deg_rates_t3):
    ax_b.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
             f'{val:.1f}%', ha='center', fontsize=7)

ax_b.set_xticks(x_pos)
ax_b.set_xticklabels([g.replace('-', '\n') for g in bar_groups], fontsize=9)
ax_b.set_ylabel('DEG Rate (%)', fontsize=11)
ax_b.set_title('B) DEG Rate by Motif Proximity', fontsize=12, fontweight='bold')
ax_b.legend(fontsize=9)

# Panel C: Functional category enrichment in each motif group (from H15)
ax_c = fig.add_subplot(gs[1, 0])
h15_data = pd.read_csv(H15_ENRICHMENT, sep='\t')
h15_t1 = h15_data[h15_data['motif_group'].str.contains('T1')]

# Plot fold enrichment for key categories
key_cats = ['Transport', 'Secondary metabolism', 'Stress/Defense', 'Primary metabolism',
            'DNA/RNA metabolism', 'Regulatory/TF', 'Hypothetical']
motif_labels = {'GCCGGC_4mC_T1': 'GCCGGC 4mC', 'CCGG_other_4mC_T1': 'CCGG other 4mC'}

# Add AAGCCCG if available
for mg in h15_t1['motif_group'].unique():
    if 'AAGCCCG' in mg and 'T1' in mg:
        motif_labels[mg] = mg.replace('_T1', '').replace('_', ' ')

y_pos = np.arange(len(key_cats))
bar_h = 0.25
plotted_motifs = []

for i, (mg, label) in enumerate(motif_labels.items()):
    sub = h15_t1[h15_t1['motif_group'] == mg]
    if len(sub) == 0:
        continue
    folds = []
    for cat in key_cats:
        row = sub[sub['category'] == cat]
        folds.append(row.iloc[0]['fold_enrichment'] if len(row) > 0 else 1.0)
    color = COLORS['GCCGGC'] if 'GCCGGC' in mg else (COLORS['AAGCCCG'] if 'AAGCCCG' in mg else '#4DAF4A')
    ax_c.barh(y_pos + i*bar_h - 0.25, folds, bar_h, label=label, color=color, alpha=0.7, edgecolor='black', linewidth=0.5)
    plotted_motifs.append(mg)

ax_c.axvline(1.0, color='grey', linestyle='--', alpha=0.5)
ax_c.set_yticks(y_pos)
ax_c.set_yticklabels(key_cats, fontsize=9)
ax_c.set_xlabel('Fold Enrichment (from H15)', fontsize=11)
ax_c.set_title('C) Functional Category Enrichment by Motif (H15)', fontsize=12, fontweight='bold')
ax_c.legend(fontsize=8, loc='upper right')

# Panel D: Summary statistics table
ax_d = fig.add_subplot(gs[1, 1])
ax_d.axis('off')

# Create summary text
text_lines = [
    'SUMMARY OF MOTIF-SPECIFIC EXPRESSION PATTERNS',
    '=' * 50,
    '',
]

for _, row in summary_df.iterrows():
    text_lines.append(f"{row['group']} (n={int(row['n_genes'])})")
    text_lines.append(f"  T2vsT1: median LFC = {row['median_LFC_T2vsT1']:.3f}")
    text_lines.append(f"    DEG rate = {100*row['DEG_rate_T2vsT1']:.1f}%  (up={int(row['n_up_T2'])}, down={int(row['n_down_T2'])})")
    text_lines.append(f"  T3vsT1: median LFC = {row['median_LFC_T3vsT1']:.3f}")
    text_lines.append(f"    DEG rate = {100*row['DEG_rate_T3vsT1']:.1f}%  (up={int(row['n_up_T3'])}, down={int(row['n_down_T3'])})")
    text_lines.append('')

# Add key test results
text_lines.append('KEY STATISTICAL TESTS')
text_lines.append('-' * 40)
key_test = stat_tests_df[stat_tests_df['comparison'] == 'AAGCCCG-only vs GCCGGC-only']
for _, row in key_test.iterrows():
    sig = "***" if row['p_bonferroni'] < 0.001 else "**" if row['p_bonferroni'] < 0.01 else "*" if row['p_bonferroni'] < 0.05 else "n.s."
    text_lines.append(f"AAGCCCG vs GCCGGC [{row['metric']}]:")
    text_lines.append(f"  p={row['p_value']:.4e}, r={row['rank_biserial_r']:.3f} {sig}")

text_lines.append('')
text_lines.append(f"Overlap: {len(overlap_genes)} dual-targeted genes")
text_lines.append(f"Jaccard index: {jaccard:.3f}")

ax_d.text(0.05, 0.95, '\n'.join(text_lines), transform=ax_d.transAxes,
         fontsize=8, verticalalignment='top', fontfamily='monospace',
         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
ax_d.set_title('D) Summary Statistics', fontsize=12, fontweight='bold')

plt.suptitle('H17: Motif-Specific Division of Labor in Gene Targeting',
            fontsize=15, fontweight='bold', y=1.01)
plt.savefig(f"{FIG_DIR}/H17_comprehensive_summary.pdf", dpi=300, bbox_inches='tight')
plt.savefig(f"{FIG_DIR}/H17_comprehensive_summary.svg", dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: H17_comprehensive_summary.pdf/svg")

# ============================================================
# 15. Final assessment
# ============================================================
print("\n" + "=" * 70)
print("HYPOTHESIS ASSESSMENT")
print("=" * 70)

# Determine verdict
aagcccg_stats = summary_df[summary_df['group'] == 'AAGCCCG-only'].iloc[0] if len(summary_df[summary_df['group'] == 'AAGCCCG-only']) > 0 else None
gccggc_stats = summary_df[summary_df['group'] == 'GCCGGC-only'].iloc[0] if len(summary_df[summary_df['group'] == 'GCCGGC-only']) > 0 else None

# Key test: AAGCCCG vs GCCGGC LFC difference
key_test_t2 = stat_tests_df[(stat_tests_df['comparison'] == 'AAGCCCG-only vs GCCGGC-only') &
                             (stat_tests_df['metric'] == 'LFC_T2vsT1')]
key_test_t3 = stat_tests_df[(stat_tests_df['comparison'] == 'AAGCCCG-only vs GCCGGC-only') &
                             (stat_tests_df['metric'] == 'LFC_T3vsT1')]

print("\n  Hypothesis: GCCGGC targets undergo T2/T3 upregulation,")
print("             AAGCCCG targets maintain constitutive/T1-biased expression.")
print()

if aagcccg_stats is not None and gccggc_stats is not None:
    print(f"  GCCGGC-only median LFC: T2vsT1={gccggc_stats['median_LFC_T2vsT1']:.3f}, T3vsT1={gccggc_stats['median_LFC_T3vsT1']:.3f}")
    print(f"  AAGCCCG-only median LFC: T2vsT1={aagcccg_stats['median_LFC_T2vsT1']:.3f}, T3vsT1={aagcccg_stats['median_LFC_T3vsT1']:.3f}")
    print()

    # Check if GCCGGC shows T2/T3 upregulation trend
    gccggc_upward = gccggc_stats['median_LFC_T3vsT1'] > 0
    # Check if AAGCCCG is constitutive/T1-biased
    aagcccg_constitutive = abs(aagcccg_stats['median_LFC_T3vsT1']) < abs(gccggc_stats['median_LFC_T3vsT1'])

    # Check statistical significance
    if len(key_test_t2) > 0:
        p_t2 = key_test_t2.iloc[0]['p_value']
        p_t3 = key_test_t3.iloc[0]['p_value'] if len(key_test_t3) > 0 else 1.0

        sig_t2 = p_t2 < 0.05
        sig_t3 = p_t3 < 0.05

        print(f"  T2vsT1 AAGCCCG vs GCCGGC: p={p_t2:.4e} ({'significant' if sig_t2 else 'not significant'})")
        print(f"  T3vsT1 AAGCCCG vs GCCGGC: p={p_t3:.4e} ({'significant' if sig_t3 else 'not significant'})")
    else:
        sig_t2 = sig_t3 = False

    # DEG rate comparison
    deg_diff_t2 = abs(gccggc_stats['DEG_rate_T2vsT1'] - aagcccg_stats['DEG_rate_T2vsT1'])
    deg_diff_t3 = abs(gccggc_stats['DEG_rate_T3vsT1'] - aagcccg_stats['DEG_rate_T3vsT1'])
    print(f"\n  DEG rate difference T2: {100*deg_diff_t2:.1f}%")
    print(f"  DEG rate difference T3: {100*deg_diff_t3:.1f}%")

    # Determine verdict
    evidence_for = 0
    evidence_against = 0

    if gccggc_upward:
        evidence_for += 1
        print("\n  [+] GCCGGC-proximal genes show positive median LFC at T3")
    else:
        evidence_against += 1
        print("\n  [-] GCCGGC-proximal genes do NOT show positive median LFC at T3")

    if aagcccg_constitutive:
        evidence_for += 1
        print("  [+] AAGCCCG-proximal genes show smaller magnitude LFC than GCCGGC")
    else:
        evidence_against += 1
        print("  [-] AAGCCCG-proximal genes do NOT show smaller magnitude LFC than GCCGGC")

    if sig_t2 or sig_t3:
        evidence_for += 1
        print("  [+] Significant difference in LFC between motif groups detected")
    else:
        evidence_against += 1
        print("  [-] No significant LFC difference between motif groups")

    # Temporal class distribution differences
    aagcccg_constitutive_rate = temporal_pct.loc['AAGCCCG-only', 'Constitutive'] if 'AAGCCCG-only' in temporal_pct.index and 'Constitutive' in temporal_pct.columns else 0
    gccggc_constitutive_rate = temporal_pct.loc['GCCGGC-only', 'Constitutive'] if 'GCCGGC-only' in temporal_pct.index and 'Constitutive' in temporal_pct.columns else 0

    if aagcccg_constitutive_rate > gccggc_constitutive_rate:
        evidence_for += 1
        print(f"  [+] AAGCCCG genes have higher constitutive rate ({aagcccg_constitutive_rate:.1f}% vs {gccggc_constitutive_rate:.1f}%)")
    else:
        evidence_against += 1
        print(f"  [-] AAGCCCG genes do NOT have higher constitutive rate ({aagcccg_constitutive_rate:.1f}% vs {gccggc_constitutive_rate:.1f}%)")

print(f"\n  Evidence FOR: {evidence_for}/4")
print(f"  Evidence AGAINST: {evidence_against}/4")

if evidence_for >= 3 and (sig_t2 or sig_t3):
    verdict = "SUPPORTED"
elif evidence_for >= 2:
    verdict = "PARTIAL"
else:
    verdict = "REJECTED"

print(f"\n  VERDICT: {verdict}")

print("\n" + "=" * 70)
print("Analysis complete.")
print("=" * 70)
