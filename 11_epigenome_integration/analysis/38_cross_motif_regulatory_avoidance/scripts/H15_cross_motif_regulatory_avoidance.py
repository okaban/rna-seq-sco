#!/usr/bin/env python3
"""
H15: Cross-Motif Regulatory Avoidance Analysis

Tests whether the "regulatory gene avoidance" pattern found for AAGCCCG (6mA)
sites in H13 is universal across all methylation motifs in S. coelicolor M145.

Key question: Do GCCGGC/CCGG 4mC sites show the same functional category
distribution as AAGCCCG -- depleted near regulatory/TF genes and hypothetical
proteins?
"""

import pandas as pd
import numpy as np
from scipy import stats
from collections import defaultdict, OrderedDict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch, FancyBboxPatch
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Configuration
# ============================================================
BASE_DIR = "/Users/okaban/bioinfo/rna-seq"
EPIGENOME_DIR = f"{BASE_DIR}/11_epigenome_integration"

CENSUS_4mC = f"{EPIGENOME_DIR}/analysis/23_expanded_motif_search/4mC_final_census.csv"
CENSUS_6mA = f"{EPIGENOME_DIR}/analysis/23_expanded_motif_search/6mA_final_census.csv"
GFF_FILE = "/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/genomic.gff"
H13_ENRICHMENT = f"{EPIGENOME_DIR}/analysis/36_AAGCCCG_distribution/tables/functional_enrichment.tsv"

OUT_DIR = f"{EPIGENOME_DIR}/analysis/38_cross_motif_regulatory_avoidance"
FIG_DIR = f"{OUT_DIR}/figures"
TBL_DIR = f"{OUT_DIR}/tables"

GENOME_SIZE = 8667507
WINDOW = 2000  # 2kb proximity window

np.random.seed(42)

# ============================================================
# 1. Load census data and define motif groups
# ============================================================
print("=" * 70)
print("H15: Cross-Motif Regulatory Avoidance Analysis")
print("=" * 70)

census_4mC = pd.read_csv(CENSUS_4mC)
census_6mA = pd.read_csv(CENSUS_6mA)

print(f"\n4mC census: {len(census_4mC)} sites")
print(f"6mA census: {len(census_6mA)} sites")

# Define motif groups for analysis
# GCCGGC-family: TGGCCGGC + GGCCGG (these are the main GCCGGC 4mC motifs)
# CCGG (other): non-GCCGGC CCGG sites
# AAGCCCG 4mC: AAGCCCG sites in the 4mC census
# AAGCCCG 6mA: AAGCCCG sites in the 6mA census

# GCCGGC-family 4mC
gccggc_mask = census_4mC['final_motif'].isin(['TGGCCGGC', 'GGCCGG'])
ccgg_other_mask = census_4mC['final_motif'] == 'CCGG (other)'
aagcccg_4mC_mask = census_4mC['final_motif'] == 'AAGCCCG'

# AAGCCCG 6mA
aagcccg_6mA_mask = census_6mA['final_motif'] == 'AAGCCCG'

# Build motif group dictionary: {group_name: DataFrame}
motif_groups = OrderedDict()

# GCCGGC family (TGGCCGGC + GGCCGG) from 4mC
for tp in ['T1', 'T2', 'T3']:
    sub = census_4mC[gccggc_mask & (census_4mC['timepoint'] == tp)].copy()
    if len(sub) > 0:
        motif_groups[f'GCCGGC_4mC_{tp}'] = sub

# CCGG (other) from 4mC
for tp in ['T1', 'T2', 'T3']:
    sub = census_4mC[ccgg_other_mask & (census_4mC['timepoint'] == tp)].copy()
    if len(sub) > 0:
        motif_groups[f'CCGG_other_4mC_{tp}'] = sub

# AAGCCCG from 4mC
for tp in ['T1', 'T2', 'T3']:
    sub = census_4mC[aagcccg_4mC_mask & (census_4mC['timepoint'] == tp)].copy()
    if len(sub) > 0:
        motif_groups[f'AAGCCCG_4mC_{tp}'] = sub

# AAGCCCG from 6mA
for tp in ['T1', 'T2', 'T3']:
    sub = census_6mA[aagcccg_6mA_mask & (census_6mA['timepoint'] == tp)].copy()
    if len(sub) > 0:
        motif_groups[f'AAGCCCG_6mA_{tp}'] = sub

# All 4mC (any motif)
for tp in ['T1', 'T2', 'T3']:
    sub = census_4mC[census_4mC['timepoint'] == tp].copy()
    if len(sub) > 0:
        motif_groups[f'All_4mC_{tp}'] = sub

print("\nMotif groups:")
for name, df in motif_groups.items():
    print(f"  {name}: {len(df)} sites")

# ============================================================
# 2. Load GFF gene annotations (same method as H13)
# ============================================================
print("\n" + "=" * 70)
print("Loading gene annotations from GFF...")
print("=" * 70)

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
print(f"Total genes loaded: {len(genes_df)}")

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
        if lt and lt not in cds_info:
            cds_info[lt] = product

genes_df['product'] = genes_df['locus_tag'].map(cds_info).fillna('unknown')
protein_coding = genes_df[genes_df['gene_biotype'] == 'protein_coding'].copy()
print(f"Protein-coding genes: {len(protein_coding)}")

# ============================================================
# 3. Functional category classification (SAME as H13)
# ============================================================

def classify_product(product):
    """Classify gene product into functional categories (identical to H13)."""
    product_lower = product.lower()

    # Regulatory
    if any(w in product_lower for w in ['transcriptional regulator', 'transcription factor',
                                          'sigma factor', 'response regulator', 'sensor kinase',
                                          'two-component', 'anti-sigma', 'repressor',
                                          'dna-binding', 'regulatory']):
        return 'Regulatory/TF'

    # Transporters
    if any(w in product_lower for w in ['transporter', 'permease', 'abc transport', 'mfs ',
                                          'efflux', 'porin', 'channel']):
        return 'Transport'

    # Biosynthesis / secondary metabolism
    if any(w in product_lower for w in ['polyketide', 'non-ribosomal', 'nrps', 'synthase',
                                          'ketoacyl', 'acyl carrier', 'thioesterase',
                                          'actinorhodin', 'undecylprodigiosin',
                                          'antibiotic biosynthesis']):
        return 'Secondary metabolism'

    # Stress/defense
    if any(w in product_lower for w in ['stress', 'heat shock', 'cold shock', 'chaperon',
                                          'protease', 'peptidase', 'dnaj', 'dnak', 'grpe',
                                          'catalase', 'superoxide', 'thioredoxin',
                                          'glutaredoxin', 'universal stress']):
        return 'Stress/Defense'

    # DNA/RNA metabolism
    if any(w in product_lower for w in ['methyltransferase', 'helicase', 'dnase', 'rnase',
                                          'topoisomerase', 'gyrase', 'recombinase',
                                          'restriction', 'ligase', 'polymerase',
                                          'nuclease', 'dna repair']):
        return 'DNA/RNA metabolism'

    # Metabolism (primary)
    if any(w in product_lower for w in ['dehydrogenase', 'oxidoreductase', 'kinase',
                                          'transferase', 'hydrolase', 'lyase', 'isomerase',
                                          'reductase', 'oxidase', 'synthetase',
                                          'phosphatase', 'esterase']):
        return 'Primary metabolism'

    # Membrane/cell wall
    if any(w in product_lower for w in ['membrane', 'cell wall', 'peptidoglycan',
                                          'lipopolysaccharide', 'penicillin-binding',
                                          'murein']):
        return 'Membrane/Cell wall'

    # Ribosome/translation
    if any(w in product_lower for w in ['ribosom', 'trna', 'rrna', 'translation',
                                          'elongation factor', 'initiation factor']):
        return 'Translation'

    # Hypothetical
    if 'hypothetical' in product_lower or product_lower == 'unknown':
        return 'Hypothetical'

    return 'Other'

# Classify ALL protein-coding genes for genome-wide expected proportions
protein_coding = protein_coding.copy()
protein_coding['functional_category'] = protein_coding['product'].apply(classify_product)
genome_cat_counts = protein_coding['functional_category'].value_counts()
genome_cat_frac = genome_cat_counts / len(protein_coding)

CATEGORIES = ['Regulatory/TF', 'Hypothetical', 'Primary metabolism', 'Other',
              'Transport', 'DNA/RNA metabolism', 'Stress/Defense',
              'Secondary metabolism', 'Translation', 'Membrane/Cell wall']

print(f"\nGenome-wide functional category distribution ({len(protein_coding)} protein-coding genes):")
for cat in CATEGORIES:
    cnt = genome_cat_counts.get(cat, 0)
    frac = genome_cat_frac.get(cat, 0)
    print(f"  {cat}: {cnt} ({frac*100:.1f}%)")

# ============================================================
# 4. Map methylation sites to proximal genes
# ============================================================
print("\n" + "=" * 70)
print("Mapping methylation sites to proximal genes (2kb window)...")
print("=" * 70)

sorted_genes = protein_coding.sort_values('start').reset_index(drop=True)

def find_proximal_genes(positions, strands, genes_sorted, window=2000):
    """
    For each methylation site, find ALL genes within 'window' bp.
    Returns DataFrame of site-gene pairs.
    """
    gene_starts = genes_sorted['start'].values
    gene_ends = genes_sorted['end'].values
    results = []

    for pos, strand in zip(positions, strands):
        # Find genes within window
        nearby = genes_sorted[
            (genes_sorted['start'] - window <= pos) &
            (genes_sorted['end'] + window >= pos)
        ]

        if len(nearby) == 0:
            # Find absolute nearest gene
            mid_points = (gene_starts + gene_ends) / 2
            dists = np.abs(mid_points - pos)
            nearest_idx = np.argmin(dists)
            nearest_gene = genes_sorted.iloc[nearest_idx]
            dist = min(abs(pos - nearest_gene['start']), abs(pos - nearest_gene['end']))
            results.append({
                'position': pos,
                'site_strand': strand,
                'locus_tag': nearest_gene['locus_tag'],
                'gene_start': nearest_gene['start'],
                'gene_end': nearest_gene['end'],
                'gene_strand': nearest_gene['strand'],
                'product': nearest_gene['product'],
                'functional_category': nearest_gene['functional_category'],
                'distance': dist,
                'within_2kb': dist <= window
            })
            continue

        for _, gene in nearby.iterrows():
            g_start, g_end = gene['start'], gene['end']
            if g_start <= pos <= g_end:
                dist = 0
            elif pos < g_start:
                dist = g_start - pos
            else:
                dist = pos - g_end

            results.append({
                'position': pos,
                'site_strand': strand,
                'locus_tag': gene['locus_tag'],
                'gene_start': gene['start'],
                'gene_end': gene['end'],
                'gene_strand': gene['strand'],
                'product': gene['product'],
                'functional_category': gene['functional_category'],
                'distance': dist,
                'within_2kb': dist <= window
            })

    return pd.DataFrame(results) if results else pd.DataFrame()


# Map all motif groups
motif_mappings = {}
for name, df in motif_groups.items():
    if len(df) > 0:
        mapping = find_proximal_genes(
            df['position'].values,
            df['strand'].values,
            sorted_genes,
            window=WINDOW
        )
        if len(mapping) > 0:
            motif_mappings[name] = mapping
            # Get nearest gene per site for category counting
            nearest = mapping.sort_values('distance').groupby('position').first().reset_index()
            print(f"  {name}: {len(df)} sites -> {len(nearest)} unique site-gene pairs (nearest)")
        else:
            print(f"  {name}: {len(df)} sites -> 0 mappings")

# ============================================================
# 5. Fisher's exact test: functional enrichment per motif
# ============================================================
print("\n" + "=" * 70)
print("Fisher's exact test for functional enrichment per motif group...")
print("=" * 70)

all_enrichment_results = []

for group_name, mapping in motif_mappings.items():
    nearest = mapping.sort_values('distance').groupby('position').first().reset_index()
    n_sites = len(nearest)
    site_cat_counts = nearest['functional_category'].value_counts()

    for cat in CATEGORIES:
        obs_in = site_cat_counts.get(cat, 0)
        obs_out = n_sites - obs_in
        genome_in = genome_cat_counts.get(cat, 0)
        genome_out = len(protein_coding) - genome_in

        # Fisher's exact test
        table = [[obs_in, obs_out], [genome_in, genome_out]]
        odds_ratio, p_val = stats.fisher_exact(table)

        expected = n_sites * (genome_in / len(protein_coding))
        fold_enrich = obs_in / expected if expected > 0 else 0

        # Confidence interval for odds ratio (Woolf's method)
        # log(OR) +/- 1.96 * sqrt(1/a + 1/b + 1/c + 1/d) where table = [[a,b],[c,d]]
        a, b, c, d = obs_in, obs_out, genome_in, genome_out
        if a > 0 and b > 0 and c > 0 and d > 0:
            log_or = np.log(odds_ratio)
            se = np.sqrt(1/a + 1/b + 1/c + 1/d)
            ci_low = np.exp(log_or - 1.96 * se)
            ci_high = np.exp(log_or + 1.96 * se)
        else:
            ci_low = 0
            ci_high = np.inf

        all_enrichment_results.append({
            'motif_group': group_name,
            'category': cat,
            'n_sites': n_sites,
            'observed': obs_in,
            'expected': round(expected, 1),
            'fold_enrichment': round(fold_enrich, 2),
            'odds_ratio': round(odds_ratio, 3),
            'OR_CI_low': round(ci_low, 3),
            'OR_CI_high': round(ci_high, 3),
            'p_value': p_val,
            'direction': 'enriched' if fold_enrich > 1 else 'depleted'
        })

enrichment_df = pd.DataFrame(all_enrichment_results)

# Bonferroni correction: n_categories (10) per motif group
enrichment_df['p_bonferroni'] = enrichment_df.groupby('motif_group')['p_value'].transform(
    lambda x: np.minimum(x * len(CATEGORIES), 1.0)
)

# Print summary table
print(f"\n{'Motif Group':<25} {'Category':<20} {'Obs':>4} {'Exp':>6} {'Fold':>5} {'OR':>6} {'p_bonf':>10} {'Sig':>4}")
print("-" * 90)
for _, row in enrichment_df[enrichment_df['p_bonferroni'] < 0.1].sort_values(['motif_group', 'p_bonferroni']).iterrows():
    sig = '***' if row['p_bonferroni'] < 0.001 else ('**' if row['p_bonferroni'] < 0.01 else ('*' if row['p_bonferroni'] < 0.05 else '.'))
    print(f"{row['motif_group']:<25} {row['category']:<20} {row['observed']:>4} "
          f"{row['expected']:>6} {row['fold_enrichment']:>5.2f} {row['odds_ratio']:>6.3f} "
          f"{row['p_bonferroni']:>10.4f} {sig:>4}")

# ============================================================
# 6. Cross-motif comparison: focus on key categories
# ============================================================
print("\n" + "=" * 70)
print("Cross-motif comparison summary...")
print("=" * 70)

# Focus on T1 timepoint for comparison (largest sample sizes)
focus_groups = ['AAGCCCG_6mA_T1', 'GCCGGC_4mC_T1', 'CCGG_other_4mC_T1',
                'AAGCCCG_4mC_T1', 'All_4mC_T1']
focus_cats = ['Regulatory/TF', 'Hypothetical', 'Primary metabolism',
              'Stress/Defense', 'Secondary metabolism']

print(f"\nT1 comparison (key categories):")
print(f"{'Category':<22}", end="")
for grp in focus_groups:
    label = grp.replace('_T1', '').replace('_', ' ')
    print(f"  {label:>18}", end="")
print()
print("-" * (22 + 20 * len(focus_groups)))

for cat in focus_cats:
    print(f"{cat:<22}", end="")
    for grp in focus_groups:
        row = enrichment_df[(enrichment_df['motif_group'] == grp) &
                            (enrichment_df['category'] == cat)]
        if len(row) > 0:
            fe = row['fold_enrichment'].values[0]
            p = row['p_bonferroni'].values[0]
            sig = '*' if p < 0.05 else ''
            print(f"  {fe:>6.2f}{sig:<11}", end="")
        else:
            print(f"  {'N/A':>18}", end="")
    print()

# ============================================================
# 7. Regulatory family breakdown by motif
# ============================================================
print("\n" + "=" * 70)
print("Regulatory family breakdown by motif...")
print("=" * 70)

# Define regulatory subfamilies
def classify_reg_family(product):
    """Classify regulatory gene into specific family."""
    product_lower = product.lower()
    families = OrderedDict([
        ('TetR', ['tetr']),
        ('MerR', ['merr']),
        ('LysR', ['lysr']),
        ('AraC', ['arac']),
        ('GntR', ['gntr']),
        ('MarR', ['marr']),
        ('LacI', ['laci']),
        ('IclR', ['iclr']),
        ('ArsR', ['arsr']),
        ('WhiB', ['whib', 'wbl']),
        ('Sigma factor', ['sigma factor']),
        ('Response regulator', ['response regulator']),
        ('Sensor kinase', ['sensor kinase']),
        ('Two-component', ['two-component']),
        ('XRE', ['xre']),
        ('SARP', ['sarp', 'streptomyces antibiotic regulatory']),
    ])
    for family, keywords in families.items():
        if any(kw in product_lower for kw in keywords):
            return family
    if classify_product(product) == 'Regulatory/TF':
        return 'Other_regulatory'
    return None

# Get regulatory genes from genome
protein_coding_copy = protein_coding.copy()
protein_coding_copy['reg_family'] = protein_coding_copy['product'].apply(classify_reg_family)
all_reg_genes = protein_coding_copy[protein_coding_copy['reg_family'].notna()].copy()

print(f"\nGenome-wide regulatory gene families:")
reg_fam_counts = all_reg_genes['reg_family'].value_counts()
for fam, cnt in reg_fam_counts.items():
    print(f"  {fam}: {cnt}")

# For each motif group (T1 focus), count regulatory family breakdown
reg_family_results = []

for group_name in ['AAGCCCG_6mA_T1', 'GCCGGC_4mC_T1', 'CCGG_other_4mC_T1', 'All_4mC_T1']:
    if group_name not in motif_mappings:
        continue
    mapping = motif_mappings[group_name]
    within_2kb = mapping[mapping['within_2kb']].copy()

    # Classify reg families for proximal genes
    within_2kb_products = within_2kb.copy()
    within_2kb_products['reg_family'] = within_2kb_products['product'].apply(classify_reg_family)
    reg_proximal = within_2kb_products[within_2kb_products['reg_family'].notna()]

    n_total_sites = mapping['position'].nunique()
    n_near_reg = reg_proximal['position'].nunique()

    print(f"\n{group_name}: {n_near_reg}/{n_total_sites} sites near regulatory genes")

    fam_counts = reg_proximal.groupby('reg_family')['position'].nunique()
    for fam in reg_fam_counts.index:
        obs = fam_counts.get(fam, 0)
        genome_fam_count = reg_fam_counts[fam]
        # Expected: proportion of this family among all genes * n_sites
        # But more accurately: among all proximal genes, what fraction is this family?
        reg_family_results.append({
            'motif_group': group_name,
            'reg_family': fam,
            'sites_near': obs,
            'genome_total': genome_fam_count,
            'total_sites': n_total_sites
        })

        if obs > 0:
            print(f"  {fam}: {obs} sites proximal")

reg_family_df = pd.DataFrame(reg_family_results)

# ============================================================
# 8. Hypothetical protein targeting analysis
# ============================================================
print("\n" + "=" * 70)
print("Hypothetical protein targeting: known-function preference test...")
print("=" * 70)

hyp_analysis_results = []

for group_name in ['AAGCCCG_6mA_T1', 'AAGCCCG_6mA_T2', 'GCCGGC_4mC_T1', 'GCCGGC_4mC_T2',
                    'CCGG_other_4mC_T1', 'All_4mC_T1', 'All_4mC_T2']:
    if group_name not in motif_mappings:
        continue
    mapping = motif_mappings[group_name]
    nearest = mapping.sort_values('distance').groupby('position').first().reset_index()

    n_total = len(nearest)
    n_hyp = (nearest['functional_category'] == 'Hypothetical').sum()
    n_known = n_total - n_hyp

    genome_hyp = genome_cat_counts.get('Hypothetical', 0)
    genome_total = len(protein_coding)
    genome_known = genome_total - genome_hyp

    # Fisher exact: test if methylation sites avoid hypothetical proteins
    table = [[n_hyp, n_known], [genome_hyp, genome_known]]
    or_val, p_val = stats.fisher_exact(table, alternative='less')  # test for depletion

    expected_hyp = n_total * (genome_hyp / genome_total)
    fold = n_hyp / expected_hyp if expected_hyp > 0 else 0

    hyp_analysis_results.append({
        'motif_group': group_name,
        'n_sites': n_total,
        'n_hypothetical': n_hyp,
        'n_known': n_known,
        'pct_hypothetical': round(n_hyp / n_total * 100, 1),
        'expected_hypothetical': round(expected_hyp, 1),
        'fold_enrichment': round(fold, 2),
        'p_depletion': p_val
    })

    print(f"  {group_name}: {n_hyp}/{n_total} hypothetical ({n_hyp/n_total*100:.1f}%) "
          f"vs genome {genome_hyp/genome_total*100:.1f}%, fold={fold:.2f}, "
          f"p(depletion)={p_val:.4f}")

hyp_df = pd.DataFrame(hyp_analysis_results)

# ============================================================
# 9. Load H13 results for comparison
# ============================================================
print("\n" + "=" * 70)
print("Loading H13 AAGCCCG enrichment results for comparison...")
print("=" * 70)

h13_enrichment = pd.read_csv(H13_ENRICHMENT, sep='\t')
h13_enrichment['motif_group'] = 'H13_AAGCCCG_6mA'
print(f"H13 results loaded: {len(h13_enrichment)} rows")

# Validate consistency: our AAGCCCG 6mA T1 should match H13 T1
our_aag_t1 = enrichment_df[enrichment_df['motif_group'] == 'AAGCCCG_6mA_T1']
h13_t1 = h13_enrichment[h13_enrichment['timepoint'] == 'T1']

print("\nValidation: Our AAGCCCG 6mA T1 vs H13 T1 (Regulatory/TF):")
our_reg = our_aag_t1[our_aag_t1['category'] == 'Regulatory/TF']
h13_reg = h13_t1[h13_t1['category'] == 'Regulatory/TF']
if len(our_reg) > 0 and len(h13_reg) > 0:
    print(f"  Our: fold={our_reg['fold_enrichment'].values[0]}, p_bonf={our_reg['p_bonferroni'].values[0]:.4f}")
    print(f"  H13: fold={h13_reg['fold_enrichment'].values[0]}, p_bonf={h13_reg['p_bonferroni'].values[0]:.4f}")

# ============================================================
# 10. Save tables
# ============================================================
print("\n" + "=" * 70)
print("Saving tables...")
print("=" * 70)

# Main cross-motif enrichment table
enrichment_df.to_csv(f"{TBL_DIR}/cross_motif_functional_enrichment.tsv", sep='\t', index=False)
print(f"Saved: cross_motif_functional_enrichment.tsv ({len(enrichment_df)} rows)")

# Summary comparison table
summary_rows = []
for grp in focus_groups:
    sub = enrichment_df[enrichment_df['motif_group'] == grp]
    if len(sub) == 0:
        continue
    n_sites = sub['n_sites'].values[0]
    for cat in CATEGORIES:
        row = sub[sub['category'] == cat]
        if len(row) > 0:
            summary_rows.append({
                'motif_group': grp,
                'n_sites': n_sites,
                'category': cat,
                'fold_enrichment': row['fold_enrichment'].values[0],
                'odds_ratio': row['odds_ratio'].values[0],
                'p_bonferroni': row['p_bonferroni'].values[0],
                'direction': row['direction'].values[0],
                'significant': row['p_bonferroni'].values[0] < 0.05
            })

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(f"{TBL_DIR}/motif_comparison_summary.tsv", sep='\t', index=False)
print(f"Saved: motif_comparison_summary.tsv ({len(summary_df)} rows)")

# Regulatory family table
reg_family_df.to_csv(f"{TBL_DIR}/regulatory_family_by_motif.tsv", sep='\t', index=False)
print(f"Saved: regulatory_family_by_motif.tsv ({len(reg_family_df)} rows)")

# Hypothetical targeting analysis
hyp_df.to_csv(f"{TBL_DIR}/hypothetical_targeting_analysis.tsv", sep='\t', index=False)
print(f"Saved: hypothetical_targeting_analysis.tsv")

# ============================================================
# 11. Generate Figures
# ============================================================
print("\n" + "=" * 70)
print("Generating Figures...")
print("=" * 70)

# Color scheme
COLORS = {
    'AAGCCCG_6mA': '#2166ac',   # blue
    'GCCGGC_4mC': '#b2182b',    # red
    'CCGG_other_4mC': '#4dac26', # green
    'AAGCCCG_4mC': '#762a83',   # purple
    'All_4mC': '#636363',       # gray
}

def get_color(group_name):
    for key, color in COLORS.items():
        if key in group_name:
            return color
    return '#636363'

# ---------- Figure 1 (Panel A): Side-by-side bar chart ----------
fig, ax = plt.subplots(figsize=(14, 6))

# T1 comparison: AAGCCCG 6mA, GCCGGC 4mC, All 4mC
bar_groups = ['AAGCCCG_6mA_T1', 'GCCGGC_4mC_T1', 'All_4mC_T1']
bar_labels = ['AAGCCCG 6mA (T1)', 'GCCGGC 4mC (T1)', 'All 4mC (T1)']
bar_colors = [COLORS['AAGCCCG_6mA'], COLORS['GCCGGC_4mC'], COLORS['All_4mC']]

x = np.arange(len(CATEGORIES))
width = 0.25
offsets = [-width, 0, width]

for i, (grp, label, color) in enumerate(zip(bar_groups, bar_labels, bar_colors)):
    sub = enrichment_df[enrichment_df['motif_group'] == grp]
    folds = [sub[sub['category'] == cat]['fold_enrichment'].values[0] if len(sub[sub['category'] == cat]) > 0 else 0
             for cat in CATEGORIES]
    pvals = [sub[sub['category'] == cat]['p_bonferroni'].values[0] if len(sub[sub['category'] == cat]) > 0 else 1
             for cat in CATEGORIES]

    bars = ax.bar(x + offsets[i], folds, width, label=label, color=color, alpha=0.75,
                  edgecolor='black', linewidth=0.4)

    # Add significance stars
    for j, (fold, p) in enumerate(zip(folds, pvals)):
        if p < 0.05:
            star = '***' if p < 0.001 else ('**' if p < 0.01 else '*')
            ax.text(x[j] + offsets[i], fold + 0.03, star, ha='center', va='bottom',
                    fontsize=8, fontweight='bold', color=color)

ax.axhline(1, color='gray', linestyle='--', linewidth=0.8, label='Expected (genome)')
ax.set_xticks(x)
ax.set_xticklabels(CATEGORIES, rotation=35, ha='right', fontsize=9)
ax.set_ylabel('Fold Enrichment (vs genome)', fontsize=11)
ax.set_title('A. Functional category enrichment by motif type (T1)', fontsize=13, fontweight='bold')
ax.legend(fontsize=9, loc='upper right')
ax.set_ylim(0, max(enrichment_df[enrichment_df['motif_group'].isin(bar_groups)]['fold_enrichment'].max() * 1.3, 2.5))

plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f"{FIG_DIR}/panel_A_bar_chart.{ext}", dpi=150, bbox_inches='tight')
plt.close()
print("Saved: panel_A_bar_chart.pdf/svg")

# ---------- Figure 2 (Panel B): Heatmap of fold enrichment ----------
fig, ax = plt.subplots(figsize=(12, 7))

# Build matrix: categories x (motif x timepoint)
heatmap_groups = ['AAGCCCG_6mA_T1', 'AAGCCCG_6mA_T2',
                  'GCCGGC_4mC_T1', 'GCCGGC_4mC_T2',
                  'CCGG_other_4mC_T1',
                  'AAGCCCG_4mC_T1',
                  'All_4mC_T1', 'All_4mC_T2']
heatmap_labels = ['AAGCCCG\n6mA T1', 'AAGCCCG\n6mA T2',
                  'GCCGGC\n4mC T1', 'GCCGGC\n4mC T2',
                  'CCGG other\n4mC T1',
                  'AAGCCCG\n4mC T1',
                  'All 4mC\nT1', 'All 4mC\nT2']

# Filter to groups that actually exist
existing_groups = [g for g in heatmap_groups if g in enrichment_df['motif_group'].unique()]
existing_labels = [heatmap_labels[heatmap_groups.index(g)] for g in existing_groups]

matrix = np.zeros((len(CATEGORIES), len(existing_groups)))
sig_matrix = np.zeros_like(matrix, dtype=bool)

for j, grp in enumerate(existing_groups):
    sub = enrichment_df[enrichment_df['motif_group'] == grp]
    for i, cat in enumerate(CATEGORIES):
        row = sub[sub['category'] == cat]
        if len(row) > 0:
            fe = row['fold_enrichment'].values[0]
            matrix[i, j] = np.log2(max(fe, 0.01))
            sig_matrix[i, j] = row['p_bonferroni'].values[0] < 0.05

im = ax.imshow(matrix, cmap='RdBu_r', aspect='auto', vmin=-2, vmax=2)
ax.set_xticks(range(len(existing_labels)))
ax.set_xticklabels(existing_labels, fontsize=9)
ax.set_yticks(range(len(CATEGORIES)))
ax.set_yticklabels(CATEGORIES, fontsize=10)

# Add fold values and significance markers
for i in range(len(CATEGORIES)):
    for j in range(len(existing_groups)):
        grp = existing_groups[j]
        sub = enrichment_df[(enrichment_df['motif_group'] == grp) &
                            (enrichment_df['category'] == CATEGORIES[i])]
        if len(sub) > 0:
            fe = sub['fold_enrichment'].values[0]
            p = sub['p_bonferroni'].values[0]
            sig = '*' if p < 0.05 else ''
            text_color = 'white' if abs(matrix[i, j]) > 1.2 else 'black'
            ax.text(j, i, f"{fe:.2f}{sig}", ha='center', va='center',
                    fontsize=8, color=text_color, fontweight='bold' if sig else 'normal')

# Add sample size info
for j, grp in enumerate(existing_groups):
    sub = enrichment_df[enrichment_df['motif_group'] == grp]
    if len(sub) > 0:
        n = sub['n_sites'].values[0]
        ax.text(j, -0.7, f"n={n}", ha='center', va='center', fontsize=8, fontstyle='italic')

cbar = plt.colorbar(im, ax=ax, shrink=0.7, pad=0.02)
cbar.set_label('log2(Fold Enrichment)', fontsize=10)

ax.set_title('B. Cross-motif functional enrichment heatmap\n(* = Bonferroni p < 0.05)',
             fontsize=13, fontweight='bold')

plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f"{FIG_DIR}/panel_B_heatmap.{ext}", dpi=150, bbox_inches='tight')
plt.close()
print("Saved: panel_B_heatmap.pdf/svg")

# ---------- Figure 3 (Panel C): Forest plot for Regulatory/TF ----------
fig, ax = plt.subplots(figsize=(10, 6))

forest_groups = [g for g in existing_groups]
forest_labels_raw = [existing_labels[existing_groups.index(g)] for g in forest_groups]
# Clean labels for forest plot
forest_labels = [l.replace('\n', ' ') for l in forest_labels_raw]

y_positions = np.arange(len(forest_groups))[::-1]

for i, grp in enumerate(forest_groups):
    sub = enrichment_df[(enrichment_df['motif_group'] == grp) &
                        (enrichment_df['category'] == 'Regulatory/TF')]
    if len(sub) > 0:
        or_val = sub['odds_ratio'].values[0]
        ci_low = sub['OR_CI_low'].values[0]
        ci_high = sub['OR_CI_high'].values[0]
        p_bonf = sub['p_bonferroni'].values[0]
        n = sub['n_sites'].values[0]

        color = get_color(grp)
        marker_size = max(4, min(12, np.sqrt(n) * 0.5))

        ax.plot(or_val, y_positions[i], 'o', color=color, markersize=marker_size,
                zorder=3, markeredgecolor='black', markeredgewidth=0.5)

        # Cap CI for display
        ci_high_disp = min(ci_high, 5.0)
        ci_low_disp = max(ci_low, 0.05)

        ax.plot([ci_low_disp, ci_high_disp], [y_positions[i], y_positions[i]],
                '-', color=color, linewidth=2, alpha=0.7)

        if ci_high > 5.0:
            ax.annotate('', xy=(5.0, y_positions[i]), xytext=(4.8, y_positions[i]),
                        arrowprops=dict(arrowstyle='->', color=color, lw=1.5))

        # significance label
        sig_label = f"  OR={or_val:.2f}, p={p_bonf:.3f}" + (" *" if p_bonf < 0.05 else "")
        ax.text(max(ci_high_disp, or_val) + 0.05, y_positions[i], sig_label,
                va='center', fontsize=8, color=color)

ax.axvline(1, color='gray', linestyle='--', linewidth=1)
ax.set_yticks(y_positions)
ax.set_yticklabels(forest_labels, fontsize=10)
ax.set_xlabel('Odds Ratio (Regulatory/TF)', fontsize=11)
ax.set_title('C. Regulatory/TF depletion across motifs\n(Forest plot: OR < 1 = depleted)',
             fontsize=13, fontweight='bold')
ax.set_xlim(0, min(5.0, ax.get_xlim()[1] + 0.5))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f"{FIG_DIR}/panel_C_forest_plot.{ext}", dpi=150, bbox_inches='tight')
plt.close()
print("Saved: panel_C_forest_plot.pdf/svg")

# ---------- Figure 4 (Panel D): Summary diagram ----------
fig, ax = plt.subplots(figsize=(12, 7))
ax.set_xlim(0, 12)
ax.set_ylim(0, 8)
ax.axis('off')

# Title
ax.text(6, 7.5, 'D. Cross-Motif Regulatory Avoidance Summary', fontsize=14,
        fontweight='bold', ha='center', va='center')

# Define summary data
# For each category, check if depleted/enriched across all motif groups at T1
universal_avoided = []
universal_enriched = []
motif_specific = []
no_pattern = []

for cat in CATEGORIES:
    t1_groups = [g for g in existing_groups if 'T1' in g]
    folds = []
    sigs = []
    for grp in t1_groups:
        sub = enrichment_df[(enrichment_df['motif_group'] == grp) &
                            (enrichment_df['category'] == cat)]
        if len(sub) > 0:
            folds.append(sub['fold_enrichment'].values[0])
            sigs.append(sub['p_bonferroni'].values[0] < 0.05)

    if len(folds) > 0:
        mean_fold = np.mean(folds)
        n_depleted = sum(1 for f in folds if f < 0.8)
        n_enriched = sum(1 for f in folds if f > 1.2)
        any_sig_depleted = any(s and f < 1 for s, f in zip(sigs, folds))
        any_sig_enriched = any(s and f > 1 for s, f in zip(sigs, folds))

        if n_depleted >= len(folds) * 0.7 and mean_fold < 0.8:
            universal_avoided.append((cat, mean_fold, any_sig_depleted))
        elif n_enriched >= len(folds) * 0.7 and mean_fold > 1.2:
            universal_enriched.append((cat, mean_fold, any_sig_enriched))
        elif any_sig_depleted or any_sig_enriched:
            motif_specific.append(cat)
        else:
            no_pattern.append(cat)

# Draw boxes
y_pos = 6.5

# Universal avoided
ax.add_patch(FancyBboxPatch((0.3, 4.5), 3.2, 2.2, boxstyle="round,pad=0.1",
             facecolor='#fee0d2', edgecolor='#cb181d', linewidth=2))
ax.text(1.9, 6.4, 'Universally Avoided', fontsize=11, fontweight='bold',
        ha='center', va='center', color='#cb181d')
for i, (cat, fold, sig) in enumerate(universal_avoided):
    marker = ' *' if sig else ''
    ax.text(1.9, 6.0 - i * 0.4, f"{cat} ({fold:.2f}x){marker}",
            fontsize=9, ha='center', va='center')

# Universal enriched
ax.add_patch(FancyBboxPatch((4.1, 4.5), 3.2, 2.2, boxstyle="round,pad=0.1",
             facecolor='#deebf7', edgecolor='#2171b5', linewidth=2))
ax.text(5.7, 6.4, 'Universally Enriched', fontsize=11, fontweight='bold',
        ha='center', va='center', color='#2171b5')
for i, (cat, fold, sig) in enumerate(universal_enriched):
    marker = ' *' if sig else ''
    ax.text(5.7, 6.0 - i * 0.4, f"{cat} ({fold:.2f}x){marker}",
            fontsize=9, ha='center', va='center')

# Motif-specific
ax.add_patch(FancyBboxPatch((7.9, 4.5), 3.8, 2.2, boxstyle="round,pad=0.1",
             facecolor='#f0f0f0', edgecolor='#636363', linewidth=2))
ax.text(9.8, 6.4, 'Motif-Specific / No Pattern', fontsize=11, fontweight='bold',
        ha='center', va='center', color='#636363')
all_other = motif_specific + no_pattern
for i, cat in enumerate(all_other[:6]):
    ax.text(9.8, 6.0 - i * 0.35, cat, fontsize=8, ha='center', va='center')

# Bottom legend
ax.text(6, 4.0, '* = Bonferroni p < 0.05 in at least one motif group',
        fontsize=9, ha='center', va='center', fontstyle='italic')
ax.text(6, 3.5, f'Analysis based on {len(existing_groups)} motif-timepoint groups, '
        f'{len(protein_coding)} protein-coding genes, 2kb proximity window',
        fontsize=8, ha='center', va='center', color='gray')

# Key finding summary
key_text = ("Key Finding: Regulatory/TF depletion is AAGCCCG-specific (6mA),\n"
            "NOT universal across 4mC motifs.\n"
            "GCCGGC 4mC sites show neutral or enriched association with regulatory genes.")

# Check if this is actually what we found
aag_reg = enrichment_df[(enrichment_df['motif_group'] == 'AAGCCCG_6mA_T1') &
                         (enrichment_df['category'] == 'Regulatory/TF')]
gcc_reg = enrichment_df[(enrichment_df['motif_group'] == 'GCCGGC_4mC_T1') &
                         (enrichment_df['category'] == 'Regulatory/TF')]

if len(aag_reg) > 0 and len(gcc_reg) > 0:
    aag_fold = aag_reg['fold_enrichment'].values[0]
    aag_p = aag_reg['p_bonferroni'].values[0]
    gcc_fold = gcc_reg['fold_enrichment'].values[0]
    gcc_p = gcc_reg['p_bonferroni'].values[0]

    if aag_fold < 0.8 and gcc_fold >= 0.8:
        # AAGCCCG depleted, GCCGGC not -> motif-specific
        verdict = "MOTIF-SPECIFIC (not universal)"
    elif aag_fold < 0.8 and gcc_fold < 0.8:
        # Both depleted -> universal
        verdict = "UNIVERSAL pattern"
    else:
        verdict = "COMPLEX pattern"

    ax.add_patch(FancyBboxPatch((1.5, 0.5), 9, 2.5, boxstyle="round,pad=0.2",
                 facecolor='#fffde7', edgecolor='#f57f17', linewidth=2))
    ax.text(6, 2.7, 'VERDICT', fontsize=12, fontweight='bold', ha='center', color='#e65100')
    ax.text(6, 2.1, f'AAGCCCG 6mA T1: Regulatory/TF fold={aag_fold:.2f}, p_bonf={aag_p:.3f}',
            fontsize=9, ha='center')
    ax.text(6, 1.7, f'GCCGGC 4mC T1: Regulatory/TF fold={gcc_fold:.2f}, p_bonf={gcc_p:.3f}',
            fontsize=9, ha='center')
    ax.text(6, 1.1, verdict, fontsize=12, fontweight='bold', ha='center', color='#bf360c')

plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f"{FIG_DIR}/panel_D_summary.{ext}", dpi=150, bbox_inches='tight')
plt.close()
print("Saved: panel_D_summary.pdf/svg")

# ---------- Figure 5: Comprehensive 4-panel figure ----------
fig = plt.figure(figsize=(18, 16))
gs = gridspec.GridSpec(2, 2, hspace=0.35, wspace=0.3)

# Panel A: Bar chart (T1 comparison)
ax = fig.add_subplot(gs[0, 0])
bar_groups_fig = ['AAGCCCG_6mA_T1', 'GCCGGC_4mC_T1', 'All_4mC_T1']
bar_labels_fig = ['AAGCCCG 6mA', 'GCCGGC 4mC', 'All 4mC']
bar_colors_fig = [COLORS['AAGCCCG_6mA'], COLORS['GCCGGC_4mC'], COLORS['All_4mC']]

x = np.arange(len(CATEGORIES))
width = 0.25
offsets = [-width, 0, width]

for i, (grp, label, color) in enumerate(zip(bar_groups_fig, bar_labels_fig, bar_colors_fig)):
    sub = enrichment_df[enrichment_df['motif_group'] == grp]
    folds = [sub[sub['category'] == cat]['fold_enrichment'].values[0]
             if len(sub[sub['category'] == cat]) > 0 else 0 for cat in CATEGORIES]
    pvals = [sub[sub['category'] == cat]['p_bonferroni'].values[0]
             if len(sub[sub['category'] == cat]) > 0 else 1 for cat in CATEGORIES]
    bars = ax.bar(x + offsets[i], folds, width, label=label, color=color, alpha=0.75,
                  edgecolor='black', linewidth=0.3)
    for j, (fold, p) in enumerate(zip(folds, pvals)):
        if p < 0.05:
            star = '***' if p < 0.001 else ('**' if p < 0.01 else '*')
            ax.text(x[j] + offsets[i], fold + 0.02, star, ha='center', va='bottom',
                    fontsize=7, fontweight='bold', color=color)

ax.axhline(1, color='gray', linestyle='--', linewidth=0.8)
ax.set_xticks(x)
ax.set_xticklabels(CATEGORIES, rotation=40, ha='right', fontsize=7)
ax.set_ylabel('Fold Enrichment', fontsize=10)
ax.set_title('A. Fold enrichment by category (T1)', fontsize=11, fontweight='bold')
ax.legend(fontsize=7, loc='upper right')

# Panel B: Heatmap
ax = fig.add_subplot(gs[0, 1])
# Recreate heatmap data
hm_groups = [g for g in existing_groups]
hm_labels = [existing_labels[existing_groups.index(g)] for g in hm_groups]

mat = np.zeros((len(CATEGORIES), len(hm_groups)))
for j, grp in enumerate(hm_groups):
    sub = enrichment_df[enrichment_df['motif_group'] == grp]
    for i, cat in enumerate(CATEGORIES):
        row = sub[sub['category'] == cat]
        if len(row) > 0:
            mat[i, j] = np.log2(max(row['fold_enrichment'].values[0], 0.01))

im = ax.imshow(mat, cmap='RdBu_r', aspect='auto', vmin=-2, vmax=2)
ax.set_xticks(range(len(hm_labels)))
ax.set_xticklabels(hm_labels, fontsize=7)
ax.set_yticks(range(len(CATEGORIES)))
ax.set_yticklabels(CATEGORIES, fontsize=8)

for i in range(len(CATEGORIES)):
    for j in range(len(hm_groups)):
        sub = enrichment_df[(enrichment_df['motif_group'] == hm_groups[j]) &
                            (enrichment_df['category'] == CATEGORIES[i])]
        if len(sub) > 0:
            fe = sub['fold_enrichment'].values[0]
            p = sub['p_bonferroni'].values[0]
            sig = '*' if p < 0.05 else ''
            tc = 'white' if abs(mat[i, j]) > 1.2 else 'black'
            ax.text(j, i, f"{fe:.1f}{sig}", ha='center', va='center',
                    fontsize=6, color=tc, fontweight='bold' if sig else 'normal')

cbar = plt.colorbar(im, ax=ax, shrink=0.7, pad=0.02)
cbar.set_label('log2(FE)', fontsize=8)
ax.set_title('B. Heatmap: all motifs x categories', fontsize=11, fontweight='bold')

# Panel C: Forest plot (Regulatory/TF)
ax = fig.add_subplot(gs[1, 0])
y_pos_fig = np.arange(len(forest_groups))[::-1]

for i, grp in enumerate(forest_groups):
    sub = enrichment_df[(enrichment_df['motif_group'] == grp) &
                        (enrichment_df['category'] == 'Regulatory/TF')]
    if len(sub) > 0:
        or_val = sub['odds_ratio'].values[0]
        ci_low = sub['OR_CI_low'].values[0]
        ci_high = sub['OR_CI_high'].values[0]
        p_bonf = sub['p_bonferroni'].values[0]
        n = sub['n_sites'].values[0]

        color = get_color(grp)
        ms = max(4, min(10, np.sqrt(n) * 0.4))

        ax.plot(or_val, y_pos_fig[i], 'o', color=color, markersize=ms,
                zorder=3, markeredgecolor='black', markeredgewidth=0.5)

        ci_h = min(ci_high, 4.0)
        ci_l = max(ci_low, 0.05)
        ax.plot([ci_l, ci_h], [y_pos_fig[i], y_pos_fig[i]],
                '-', color=color, linewidth=2, alpha=0.7)

        sig_l = f" p={p_bonf:.3f}" + (" *" if p_bonf < 0.05 else "")
        ax.text(max(ci_h, or_val) + 0.05, y_pos_fig[i], sig_l,
                va='center', fontsize=7, color=color)

ax.axvline(1, color='gray', linestyle='--', linewidth=1)
ax.set_yticks(y_pos_fig)
ax.set_yticklabels(forest_labels, fontsize=8)
ax.set_xlabel('Odds Ratio (Regulatory/TF)', fontsize=10)
ax.set_title('C. Regulatory/TF OR by motif', fontsize=11, fontweight='bold')
ax.set_xlim(0, min(4.0, ax.get_xlim()[1]))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel D: Hypothetical depletion comparison
ax = fig.add_subplot(gs[1, 1])
hyp_groups = ['AAGCCCG_6mA_T1', 'GCCGGC_4mC_T1', 'CCGG_other_4mC_T1', 'All_4mC_T1']
hyp_labels_fig = ['AAGCCCG\n6mA', 'GCCGGC\n4mC', 'CCGG other\n4mC', 'All\n4mC']
hyp_colors_fig = [get_color(g) for g in hyp_groups]

hyp_folds = []
hyp_pvals = []
genome_hyp_frac = genome_cat_frac.get('Hypothetical', 0) * 100

for grp in hyp_groups:
    sub = enrichment_df[(enrichment_df['motif_group'] == grp) &
                        (enrichment_df['category'] == 'Hypothetical')]
    if len(sub) > 0:
        hyp_folds.append(sub['fold_enrichment'].values[0])
        hyp_pvals.append(sub['p_bonferroni'].values[0])
    else:
        hyp_folds.append(0)
        hyp_pvals.append(1)

bars = ax.bar(range(len(hyp_groups)), hyp_folds, color=hyp_colors_fig, alpha=0.75,
              edgecolor='black', linewidth=0.5)
ax.axhline(1, color='gray', linestyle='--', linewidth=0.8, label='Expected')

for i, (fold, p) in enumerate(zip(hyp_folds, hyp_pvals)):
    if p < 0.05:
        ax.text(i, fold + 0.03, '*', ha='center', fontsize=12, fontweight='bold')
    ax.text(i, fold - 0.05, f'{fold:.2f}', ha='center', va='top', fontsize=8,
            color='white' if fold > 0.5 else 'black')

ax.set_xticks(range(len(hyp_groups)))
ax.set_xticklabels(hyp_labels_fig, fontsize=8)
ax.set_ylabel('Fold Enrichment', fontsize=10)
ax.set_title('D. Hypothetical protein depletion (T1)', fontsize=11, fontweight='bold')
ax.legend(fontsize=8)

plt.suptitle('H15: Cross-Motif Regulatory Avoidance Analysis\nS. coelicolor M145',
             fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f"{FIG_DIR}/H15_comprehensive_summary.{ext}", dpi=150, bbox_inches='tight')
plt.close()
print("Saved: H15_comprehensive_summary.pdf/svg")

# ============================================================
# 12. Final statistical summary
# ============================================================
print("\n" + "=" * 70)
print("FINAL SUMMARY")
print("=" * 70)

# Key comparisons
print("\n--- Regulatory/TF depletion test ---")
for grp in ['AAGCCCG_6mA_T1', 'AAGCCCG_6mA_T2', 'GCCGGC_4mC_T1', 'GCCGGC_4mC_T2',
            'CCGG_other_4mC_T1', 'AAGCCCG_4mC_T1', 'All_4mC_T1', 'All_4mC_T2']:
    sub = enrichment_df[(enrichment_df['motif_group'] == grp) &
                        (enrichment_df['category'] == 'Regulatory/TF')]
    if len(sub) > 0:
        row = sub.iloc[0]
        sig = "SIG" if row['p_bonferroni'] < 0.05 else "ns"
        print(f"  {grp:<25} fold={row['fold_enrichment']:.2f}, OR={row['odds_ratio']:.3f}, "
              f"p_bonf={row['p_bonferroni']:.4f} [{sig}]")

print("\n--- Hypothetical depletion test ---")
for grp in ['AAGCCCG_6mA_T1', 'AAGCCCG_6mA_T2', 'GCCGGC_4mC_T1', 'GCCGGC_4mC_T2',
            'CCGG_other_4mC_T1', 'AAGCCCG_4mC_T1', 'All_4mC_T1', 'All_4mC_T2']:
    sub = enrichment_df[(enrichment_df['motif_group'] == grp) &
                        (enrichment_df['category'] == 'Hypothetical')]
    if len(sub) > 0:
        row = sub.iloc[0]
        sig = "SIG" if row['p_bonferroni'] < 0.05 else "ns"
        print(f"  {grp:<25} fold={row['fold_enrichment']:.2f}, OR={row['odds_ratio']:.3f}, "
              f"p_bonf={row['p_bonferroni']:.4f} [{sig}]")

print("\n--- Stress/Defense enrichment test ---")
for grp in ['AAGCCCG_6mA_T1', 'GCCGGC_4mC_T1', 'All_4mC_T1']:
    sub = enrichment_df[(enrichment_df['motif_group'] == grp) &
                        (enrichment_df['category'] == 'Stress/Defense')]
    if len(sub) > 0:
        row = sub.iloc[0]
        sig = "SIG" if row['p_bonferroni'] < 0.05 else "ns"
        print(f"  {grp:<25} fold={row['fold_enrichment']:.2f}, OR={row['odds_ratio']:.3f}, "
              f"p_bonf={row['p_bonferroni']:.4f} [{sig}]")

# Determine overall verdict
print("\n" + "=" * 70)
aag_reg_fold = enrichment_df[(enrichment_df['motif_group'] == 'AAGCCCG_6mA_T1') &
                              (enrichment_df['category'] == 'Regulatory/TF')]['fold_enrichment'].values[0]
aag_reg_p = enrichment_df[(enrichment_df['motif_group'] == 'AAGCCCG_6mA_T1') &
                            (enrichment_df['category'] == 'Regulatory/TF')]['p_bonferroni'].values[0]
gcc_reg_fold = enrichment_df[(enrichment_df['motif_group'] == 'GCCGGC_4mC_T1') &
                              (enrichment_df['category'] == 'Regulatory/TF')]['fold_enrichment'].values[0]
gcc_reg_p = enrichment_df[(enrichment_df['motif_group'] == 'GCCGGC_4mC_T1') &
                            (enrichment_df['category'] == 'Regulatory/TF')]['p_bonferroni'].values[0]

aag_hyp_fold = enrichment_df[(enrichment_df['motif_group'] == 'AAGCCCG_6mA_T1') &
                              (enrichment_df['category'] == 'Hypothetical')]['fold_enrichment'].values[0]
gcc_hyp_fold = enrichment_df[(enrichment_df['motif_group'] == 'GCCGGC_4mC_T1') &
                              (enrichment_df['category'] == 'Hypothetical')]['fold_enrichment'].values[0]

print("VERDICT DETERMINATION:")
print(f"  AAGCCCG 6mA T1: Reg/TF fold={aag_reg_fold:.2f} (p_bonf={aag_reg_p:.4f}), "
      f"Hyp fold={aag_hyp_fold:.2f}")
print(f"  GCCGGC 4mC T1:  Reg/TF fold={gcc_reg_fold:.2f} (p_bonf={gcc_reg_p:.4f}), "
      f"Hyp fold={gcc_hyp_fold:.2f}")

if gcc_reg_fold < 0.8 and gcc_reg_p < 0.05:
    if aag_reg_fold < 0.8 and aag_reg_p < 0.05:
        verdict_text = "SUPPORTED - Universal regulatory avoidance across both 6mA and 4mC motifs"
    else:
        verdict_text = "PARTIAL - GCCGGC 4mC shows avoidance but AAGCCCG 6mA replication failed"
elif aag_reg_fold < 0.8 and aag_reg_p < 0.05 and gcc_reg_fold >= 0.8:
    verdict_text = "REJECTED - Regulatory avoidance is AAGCCCG-specific, NOT universal"
elif aag_reg_fold < 0.8 and gcc_reg_fold < 0.8 and gcc_reg_p >= 0.05:
    verdict_text = "PARTIAL - Trend toward universal avoidance but GCCGGC not significant"
else:
    verdict_text = "REJECTED - No clear universal pattern"

print(f"\n  H15 VERDICT: {verdict_text}")
print("=" * 70)

print("\nAnalysis complete.")
print(f"Output directory: {OUT_DIR}")
print(f"Figures: {FIG_DIR}/")
print(f"Tables: {TBL_DIR}/")
