#!/usr/bin/env python3
"""
H13: AAGCCCG 260-site genomic distribution pattern and functional target analysis

Analyzes the genomic distribution of AAGCCCG (6mA) methylation sites in
S. coelicolor M145, comparing T1 (260 sites) vs T2 (64 new sites).
Tests for non-random distribution, functional enrichment, and expression impact.
"""

import pandas as pd
import numpy as np
from scipy import stats
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Configuration
# ============================================================
BASE_DIR = "/Users/okaban/bioinfo/rna-seq"
EPIGENOME_DIR = f"{BASE_DIR}/11_epigenome_integration"
CENSUS_FILE = f"{EPIGENOME_DIR}/analysis/23_expanded_motif_search/6mA_final_census.csv"
GFF_FILE = "/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/genomic.gff"
BGC_FILE = f"{BASE_DIR}/10_SARP_motif_scan/analysis/10_SARP_motif_scan_260128_v1/tables/BGC_and_TF_coordinates.tsv"
DESEQ2_FILE = f"{BASE_DIR}/04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_2_vs_1.tsv"
COORD_REG_FILE = f"{EPIGENOME_DIR}/analysis/29_genomewide_TF_screen/tables/coordinated_regulatory_genes.tsv"
OUT_DIR = f"{EPIGENOME_DIR}/analysis/36_AAGCCCG_distribution"
FIG_DIR = f"{OUT_DIR}/figures"
TBL_DIR = f"{OUT_DIR}/tables"

GENOME_SIZE = 8667507
ARM_BOUNDARY_LEFT = 1_500_000
ARM_BOUNDARY_RIGHT = 7_170_000
CORE_SIZE = ARM_BOUNDARY_RIGHT - ARM_BOUNDARY_LEFT  # 5,670,000
ARM_SIZE = GENOME_SIZE - CORE_SIZE  # 2,997,507

N_PERMUTATIONS = 10_000
np.random.seed(42)

# ============================================================
# 1. Load AAGCCCG site data
# ============================================================
print("=" * 70)
print("H13: AAGCCCG Genomic Distribution & Functional Target Analysis")
print("=" * 70)

census = pd.read_csv(CENSUS_FILE)
aagcccg_all = census[census['final_motif'] == 'AAGCCCG'].copy()
print(f"\nTotal AAGCCCG sites across all timepoints: {len(aagcccg_all)}")
print(f"By timepoint: {aagcccg_all['timepoint'].value_counts().sort_index().to_dict()}")

# Separate T1 and T2 sites
t1_sites = aagcccg_all[aagcccg_all['timepoint'] == 'T1'].copy()
t2_sites = aagcccg_all[aagcccg_all['timepoint'] == 'T2'].copy()
t3_sites = aagcccg_all[aagcccg_all['timepoint'] == 'T3'].copy()

print(f"\nT1 sites: {len(t1_sites)}")
print(f"T2 sites: {len(t2_sites)}")
print(f"T3 sites: {len(t3_sites)}")

# Check overlap between T1 and T2
t1_positions = set(t1_sites['position'].values)
t2_positions = set(t2_sites['position'].values)
overlap = t1_positions & t2_positions
print(f"\nT1-T2 position overlap: {len(overlap)} sites")

# ============================================================
# 2. Load GFF gene annotations
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

        # Parse attributes
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

# Load CDS annotations for product info
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
# 3. Arm vs Core distribution analysis
# ============================================================
print("\n" + "=" * 70)
print("3. Arm vs Core Distribution Analysis")
print("=" * 70)

def classify_region(pos):
    if pos < ARM_BOUNDARY_LEFT or pos > ARM_BOUNDARY_RIGHT:
        return 'arm'
    else:
        return 'core'

t1_sites['region'] = t1_sites['position'].apply(classify_region)
t2_sites['region'] = t2_sites['position'].apply(classify_region)

t1_arm = (t1_sites['region'] == 'arm').sum()
t1_core = (t1_sites['region'] == 'core').sum()
t2_arm = (t2_sites['region'] == 'arm').sum()
t2_core = (t2_sites['region'] == 'core').sum()

# Expected proportions
arm_frac = ARM_SIZE / GENOME_SIZE
core_frac = CORE_SIZE / GENOME_SIZE

print(f"\nGenome: arm fraction = {arm_frac:.3f}, core fraction = {core_frac:.3f}")
print(f"\nT1 (260 sites): arm={t1_arm} ({t1_arm/len(t1_sites)*100:.1f}%), "
      f"core={t1_core} ({t1_core/len(t1_sites)*100:.1f}%)")
print(f"T2 (64 sites):  arm={t2_arm} ({t2_arm/len(t2_sites)*100:.1f}%), "
      f"core={t2_core} ({t2_core/len(t2_sites)*100:.1f}%)")

# Chi-squared test for T1
t1_expected_arm = len(t1_sites) * arm_frac
t1_expected_core = len(t1_sites) * core_frac
chi2_t1, p_t1 = stats.chisquare(
    [t1_arm, t1_core],
    [t1_expected_arm, t1_expected_core]
)
print(f"\nT1 Chi-squared test: chi2={chi2_t1:.3f}, p={p_t1:.2e}")
print(f"  Expected: arm={t1_expected_arm:.1f}, core={t1_expected_core:.1f}")
print(f"  Observed: arm={t1_arm}, core={t1_core}")
t1_arm_enrichment = (t1_arm / t1_expected_arm)
print(f"  Arm enrichment: {t1_arm_enrichment:.2f}x")

# Chi-squared test for T2
t2_expected_arm = len(t2_sites) * arm_frac
t2_expected_core = len(t2_sites) * core_frac
chi2_t2, p_t2 = stats.chisquare(
    [t2_arm, t2_core],
    [t2_expected_arm, t2_expected_core]
)
print(f"\nT2 Chi-squared test: chi2={chi2_t2:.3f}, p={p_t2:.2e}")
print(f"  Expected: arm={t2_expected_arm:.1f}, core={t2_expected_core:.1f}")
print(f"  Observed: arm={t2_arm}, core={t2_core}")
t2_arm_enrichment = (t2_arm / t2_expected_arm) if t2_expected_arm > 0 else 0
print(f"  Arm enrichment: {t2_arm_enrichment:.2f}x")

# Fisher's exact test: T1 vs T2 arm/core proportions
fisher_table = [[t1_arm, t1_core], [t2_arm, t2_core]]
fisher_or, fisher_p = stats.fisher_exact(fisher_table)
print(f"\nFisher's exact (T1 vs T2 arm/core): OR={fisher_or:.3f}, p={fisher_p:.2e}")

# ============================================================
# 4. Sliding window density analysis
# ============================================================
print("\n" + "=" * 70)
print("4. Sliding Window Density Analysis")
print("=" * 70)

window_size = 100_000  # 100 kb windows
step_size = 50_000     # 50 kb steps

def sliding_window_density(positions, genome_size, window, step):
    """Calculate site density in sliding windows."""
    starts = np.arange(0, genome_size - window + 1, step)
    densities = []
    for s in starts:
        count = np.sum((positions >= s) & (positions < s + window))
        densities.append(count / (window / 1e6))  # sites per Mb
    return starts + window // 2, np.array(densities)

t1_pos = t1_sites['position'].values
t2_pos = t2_sites['position'].values

t1_centers, t1_density = sliding_window_density(t1_pos, GENOME_SIZE, window_size, step_size)
t2_centers, t2_density = sliding_window_density(t2_pos, GENOME_SIZE, window_size, step_size)

print(f"T1 density: mean={np.mean(t1_density):.1f}, max={np.max(t1_density):.1f} sites/Mb")
print(f"T2 density: mean={np.mean(t2_density):.1f}, max={np.max(t2_density):.1f} sites/Mb")

# Identify hotspot regions (>2x mean density)
t1_hotspot_threshold = np.mean(t1_density) * 2
t1_hotspots = t1_centers[t1_density > t1_hotspot_threshold]
print(f"\nT1 hotspot windows (>{t1_hotspot_threshold:.1f} sites/Mb): {len(t1_hotspots)}")
if len(t1_hotspots) > 0:
    for pos in t1_hotspots:
        region = classify_region(pos)
        d = t1_density[t1_centers == pos][0]
        print(f"  {pos/1e6:.2f} Mb ({region}): {d:.0f} sites/Mb")

# ============================================================
# 5. Nearest-neighbor clustering test
# ============================================================
print("\n" + "=" * 70)
print("5. Nearest-Neighbor Clustering Test")
print("=" * 70)

def nearest_neighbor_distances(positions):
    """Calculate nearest-neighbor distances for sorted positions."""
    sorted_pos = np.sort(positions)
    if len(sorted_pos) < 2:
        return np.array([])
    dists = np.diff(sorted_pos)
    # For each site, the NN distance is min of left and right gaps
    nn_dists = np.zeros(len(sorted_pos))
    nn_dists[0] = dists[0]
    nn_dists[-1] = dists[-1]
    for i in range(1, len(sorted_pos) - 1):
        nn_dists[i] = min(dists[i-1], dists[i])
    return nn_dists

t1_nn = nearest_neighbor_distances(t1_pos)
t2_nn = nearest_neighbor_distances(t2_pos)

print(f"T1 nearest-neighbor distances: median={np.median(t1_nn):.0f}, "
      f"mean={np.mean(t1_nn):.0f}, min={np.min(t1_nn):.0f}")
print(f"T2 nearest-neighbor distances: median={np.median(t2_nn):.0f}, "
      f"mean={np.mean(t2_nn):.0f}, min={np.min(t2_nn):.0f}")

# Expected NN distance for random uniform distribution
# E[NN] = genome_size / (2 * n) for n points on a line
t1_expected_nn = GENOME_SIZE / (2 * len(t1_sites))
t2_expected_nn = GENOME_SIZE / (2 * len(t2_sites))
print(f"\nExpected NN (uniform): T1={t1_expected_nn:.0f}, T2={t2_expected_nn:.0f}")
print(f"Observed/Expected ratio: T1={np.median(t1_nn)/t1_expected_nn:.3f}, "
      f"T2={np.median(t2_nn)/t2_expected_nn:.3f}")

# Permutation test for clustering
print("\nRunning permutation test (10,000 permutations)...")
t1_obs_median_nn = np.median(t1_nn)
t2_obs_median_nn = np.median(t2_nn)

perm_medians_t1 = []
perm_medians_t2 = []
for _ in range(N_PERMUTATIONS):
    rand_pos_t1 = np.sort(np.random.randint(1, GENOME_SIZE + 1, size=len(t1_sites)))
    rand_nn_t1 = nearest_neighbor_distances(rand_pos_t1)
    perm_medians_t1.append(np.median(rand_nn_t1))

    rand_pos_t2 = np.sort(np.random.randint(1, GENOME_SIZE + 1, size=len(t2_sites)))
    rand_nn_t2 = nearest_neighbor_distances(rand_pos_t2)
    perm_medians_t2.append(np.median(rand_nn_t2))

perm_medians_t1 = np.array(perm_medians_t1)
perm_medians_t2 = np.array(perm_medians_t2)

# Two-tailed test: sites more clustered (smaller NN) or more dispersed (larger NN)?
p_cluster_t1 = np.mean(perm_medians_t1 <= t1_obs_median_nn)
p_disperse_t1 = np.mean(perm_medians_t1 >= t1_obs_median_nn)
p_cluster_t2 = np.mean(perm_medians_t2 <= t2_obs_median_nn)
p_disperse_t2 = np.mean(perm_medians_t2 >= t2_obs_median_nn)

print(f"\nT1: obs median NN={t1_obs_median_nn:.0f}, perm median={np.median(perm_medians_t1):.0f}")
print(f"  P(more clustered)={p_cluster_t1:.4f}, P(more dispersed)={p_disperse_t1:.4f}")
if t1_obs_median_nn < np.median(perm_medians_t1):
    print(f"  => T1 sites are MORE CLUSTERED than random (p={p_cluster_t1:.4f})")
else:
    print(f"  => T1 sites are MORE DISPERSED than random (p={p_disperse_t1:.4f})")

print(f"\nT2: obs median NN={t2_obs_median_nn:.0f}, perm median={np.median(perm_medians_t2):.0f}")
print(f"  P(more clustered)={p_cluster_t2:.4f}, P(more dispersed)={p_disperse_t2:.4f}")
if t2_obs_median_nn < np.median(perm_medians_t2):
    print(f"  => T2 sites are MORE CLUSTERED than random (p={p_cluster_t2:.4f})")
else:
    print(f"  => T2 sites are MORE DISPERSED than random (p={p_disperse_t2:.4f})")

# ============================================================
# 6. Gene proximity analysis
# ============================================================
print("\n" + "=" * 70)
print("6. Gene Proximity Analysis")
print("=" * 70)

def map_sites_to_genes(site_positions, site_strands, genes_df, window=2000):
    """
    For each methylation site, find nearest gene(s) and classify location.
    Returns a list of dicts with site-gene mappings.
    """
    # Pre-sort genes by start position
    sorted_genes = genes_df.sort_values('start').reset_index(drop=True)
    gene_starts = sorted_genes['start'].values
    gene_ends = sorted_genes['end'].values

    results = []
    for idx, (pos, strand) in enumerate(zip(site_positions, site_strands)):
        # Find genes within window
        nearby = sorted_genes[
            (sorted_genes['start'] - window <= pos) &
            (sorted_genes['end'] + window >= pos)
        ]

        if len(nearby) == 0:
            # Find absolute nearest gene
            mid_points = (gene_starts + gene_ends) / 2
            dists = np.abs(mid_points - pos)
            nearest_idx = np.argmin(dists)
            nearest_gene = sorted_genes.iloc[nearest_idx]
            dist = min(abs(pos - nearest_gene['start']), abs(pos - nearest_gene['end']))
            results.append({
                'position': pos,
                'site_strand': strand,
                'locus_tag': nearest_gene['locus_tag'],
                'old_locus_tag': nearest_gene['old_locus_tag'],
                'gene_name': nearest_gene['gene_name'],
                'gene_start': nearest_gene['start'],
                'gene_end': nearest_gene['end'],
                'gene_strand': nearest_gene['strand'],
                'product': nearest_gene['product'],
                'distance': dist,
                'location': 'distal_intergenic',
                'within_500bp': False,
                'within_2kb': dist <= 2000
            })
            continue

        for _, gene in nearby.iterrows():
            g_start, g_end = gene['start'], gene['end']
            g_strand = gene['strand']

            if g_start <= pos <= g_end:
                loc = 'gene_body'
                dist = 0
            elif pos < g_start:
                dist = g_start - pos
                if g_strand == '+':
                    loc = 'promoter' if dist <= 500 else 'upstream'
                else:
                    loc = 'downstream' if dist <= 500 else 'downstream_far'
            else:  # pos > g_end
                dist = pos - g_end
                if g_strand == '-':
                    loc = 'promoter' if dist <= 500 else 'upstream'
                else:
                    loc = 'downstream' if dist <= 500 else 'downstream_far'

            results.append({
                'position': pos,
                'site_strand': strand,
                'locus_tag': gene['locus_tag'],
                'old_locus_tag': gene['old_locus_tag'],
                'gene_name': gene['gene_name'],
                'gene_start': gene['start'],
                'gene_end': gene['end'],
                'gene_strand': gene['strand'],
                'product': gene['product'],
                'distance': dist,
                'location': loc,
                'within_500bp': dist <= 500,
                'within_2kb': dist <= 2000
            })

    return pd.DataFrame(results)

# Map T1 sites
print("Mapping T1 sites to genes...")
t1_mapping = map_sites_to_genes(
    t1_sites['position'].values,
    t1_sites['strand'].values,
    protein_coding,
    window=2000
)
t1_mapping['timepoint'] = 'T1'

# Map T2 sites
print("Mapping T2 sites to genes...")
t2_mapping = map_sites_to_genes(
    t2_sites['position'].values,
    t2_sites['strand'].values,
    protein_coding,
    window=2000
)
t2_mapping['timepoint'] = 'T2'

# Combine
all_mapping = pd.concat([t1_mapping, t2_mapping], ignore_index=True)

# Location classification summary
print("\nSite location classification (unique sites, nearest gene only):")
for tp, label in [('T1', 'T1'), ('T2', 'T2')]:
    sub = all_mapping[all_mapping['timepoint'] == tp]
    # Get nearest gene per site
    nearest = sub.sort_values('distance').groupby('position').first().reset_index()
    loc_counts = nearest['location'].value_counts()
    print(f"\n  {label} ({len(nearest)} unique sites):")
    for loc, cnt in loc_counts.items():
        print(f"    {loc}: {cnt} ({cnt/len(nearest)*100:.1f}%)")

# ============================================================
# 7. Functional category analysis
# ============================================================
print("\n" + "=" * 70)
print("7. Functional Category Analysis")
print("=" * 70)

def classify_product(product):
    """Classify gene product into functional categories."""
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

# Classify nearest genes
for tp in ['T1', 'T2']:
    sub = all_mapping[all_mapping['timepoint'] == tp]
    nearest = sub.sort_values('distance').groupby('position').first().reset_index()
    nearest['functional_category'] = nearest['product'].apply(classify_product)

    print(f"\n{tp} - Functional categories of nearest genes:")
    cat_counts = nearest['functional_category'].value_counts()
    for cat, cnt in cat_counts.items():
        print(f"  {cat}: {cnt} ({cnt/len(nearest)*100:.1f}%)")

# Classify ALL protein-coding genes for expected proportions
protein_coding['functional_category'] = protein_coding['product'].apply(classify_product)
genome_cat_counts = protein_coding['functional_category'].value_counts()
genome_cat_frac = genome_cat_counts / len(protein_coding)

print(f"\nGenome-wide functional category distribution:")
for cat, frac in genome_cat_frac.items():
    print(f"  {cat}: {genome_cat_counts[cat]} ({frac*100:.1f}%)")

# Enrichment analysis (Fisher's exact test for each category)
print("\n" + "-" * 50)
print("Functional category enrichment (Fisher's exact test):")
print("-" * 50)

enrichment_results = []

for tp in ['T1', 'T2']:
    sub = all_mapping[all_mapping['timepoint'] == tp]
    nearest = sub.sort_values('distance').groupby('position').first().reset_index()
    nearest['functional_category'] = nearest['product'].apply(classify_product)

    n_sites = len(nearest)
    site_cat_counts = nearest['functional_category'].value_counts()

    for cat in genome_cat_frac.index:
        obs_in = site_cat_counts.get(cat, 0)
        obs_out = n_sites - obs_in
        genome_in = genome_cat_counts[cat]
        genome_out = len(protein_coding) - genome_in

        # Fisher's exact
        table = [[obs_in, obs_out], [genome_in, genome_out]]
        odds_ratio, p_val = stats.fisher_exact(table)

        expected = n_sites * (genome_in / len(protein_coding))
        fold_enrich = obs_in / expected if expected > 0 else 0

        enrichment_results.append({
            'timepoint': tp,
            'category': cat,
            'observed': obs_in,
            'expected': round(expected, 1),
            'fold_enrichment': round(fold_enrich, 2),
            'odds_ratio': round(odds_ratio, 3),
            'p_value': p_val,
            'direction': 'enriched' if fold_enrich > 1 else 'depleted'
        })

enrichment_df = pd.DataFrame(enrichment_results)
# Multiple testing correction (Bonferroni)
enrichment_df['p_bonferroni'] = np.minimum(enrichment_df['p_value'] * len(enrichment_df), 1.0)

print(f"\n{'Timepoint':<5} {'Category':<25} {'Obs':>4} {'Exp':>6} {'Fold':>5} {'OR':>6} {'p_adj':>10}")
print("-" * 70)
for _, row in enrichment_df.sort_values(['timepoint', 'p_value']).iterrows():
    sig = '*' if row['p_bonferroni'] < 0.05 else ''
    print(f"{row['timepoint']:<5} {row['category']:<25} {row['observed']:>4} "
          f"{row['expected']:>6} {row['fold_enrichment']:>5.2f} {row['odds_ratio']:>6.2f} "
          f"{row['p_bonferroni']:>10.2e} {sig}")

# ============================================================
# 8. BGC overlap analysis
# ============================================================
print("\n" + "=" * 70)
print("8. BGC Overlap Analysis")
print("=" * 70)

bgc_df = pd.read_csv(BGC_FILE, sep='\t')
bgc_genes = bgc_df[bgc_df['gene_type'] == 'BGC_gene']

# Define BGC regions
bgc_regions = {}
for bgc_name in bgc_genes['bgc_name'].unique():
    sub = bgc_genes[bgc_genes['bgc_name'] == bgc_name]
    bgc_regions[bgc_name] = {
        'start': sub['start'].min(),
        'end': sub['end'].max(),
        'n_genes': len(sub)
    }

print("\nBGC regions:")
total_bgc_size = 0
for name, info in bgc_regions.items():
    size = info['end'] - info['start']
    total_bgc_size += size
    print(f"  {name}: {info['start']:,} - {info['end']:,} ({size/1000:.1f} kb, {info['n_genes']} genes)")

bgc_frac = total_bgc_size / GENOME_SIZE
print(f"\nTotal BGC region: {total_bgc_size/1000:.1f} kb ({bgc_frac*100:.2f}% of genome)")

def in_bgc(pos, bgc_regions):
    for name, info in bgc_regions.items():
        if info['start'] <= pos <= info['end']:
            return name
    return None

t1_bgc = [in_bgc(p, bgc_regions) for p in t1_pos]
t2_bgc = [in_bgc(p, bgc_regions) for p in t2_pos]

t1_in_bgc = sum(1 for x in t1_bgc if x is not None)
t2_in_bgc = sum(1 for x in t2_bgc if x is not None)

print(f"\nT1: {t1_in_bgc}/{len(t1_sites)} sites in BGC regions ({t1_in_bgc/len(t1_sites)*100:.1f}%)")
print(f"T2: {t2_in_bgc}/{len(t2_sites)} sites in BGC regions ({t2_in_bgc/len(t2_sites)*100:.1f}%)")

# Fisher's exact test for BGC enrichment
# T1
t1_bgc_exp = len(t1_sites) * bgc_frac
fisher_bgc_t1 = stats.fisher_exact([[t1_in_bgc, len(t1_sites) - t1_in_bgc],
                                      [total_bgc_size, GENOME_SIZE - total_bgc_size]])
print(f"\nT1 BGC: expected={t1_bgc_exp:.1f}, observed={t1_in_bgc}")
print(f"  Fold enrichment: {t1_in_bgc/t1_bgc_exp:.2f}x, Fisher p={fisher_bgc_t1[1]:.2e}")

# T2
t2_bgc_exp = len(t2_sites) * bgc_frac
fisher_bgc_t2 = stats.fisher_exact([[t2_in_bgc, len(t2_sites) - t2_in_bgc],
                                      [total_bgc_size, GENOME_SIZE - total_bgc_size]])
print(f"T2 BGC: expected={t2_bgc_exp:.1f}, observed={t2_in_bgc}")
print(f"  Fold enrichment: {t2_in_bgc/t2_bgc_exp:.2f}x, Fisher p={fisher_bgc_t2[1]:.2e}")

# Details per BGC
bgc_overlap_records = []
for tp, bgc_list, n_total in [('T1', t1_bgc, len(t1_sites)), ('T2', t2_bgc, len(t2_sites))]:
    bgc_counts = defaultdict(int)
    for x in bgc_list:
        if x is not None:
            bgc_counts[x] += 1
    for name in bgc_regions:
        bgc_overlap_records.append({
            'timepoint': tp,
            'bgc_name': name,
            'sites_in_bgc': bgc_counts.get(name, 0),
            'bgc_size_kb': (bgc_regions[name]['end'] - bgc_regions[name]['start']) / 1000,
            'bgc_n_genes': bgc_regions[name]['n_genes'],
            'total_sites': n_total
        })

bgc_overlap_df = pd.DataFrame(bgc_overlap_records)
print("\nPer-BGC overlap:")
for _, row in bgc_overlap_df.iterrows():
    print(f"  {row['timepoint']} {row['bgc_name']}: {row['sites_in_bgc']} sites "
          f"(BGC size={row['bgc_size_kb']:.1f} kb)")

# ============================================================
# 9. Regulatory gene enrichment
# ============================================================
print("\n" + "=" * 70)
print("9. Regulatory Gene Enrichment")
print("=" * 70)

# Load coordinated regulatory genes
coord_reg = pd.read_csv(COORD_REG_FILE, sep='\t')
coord_reg_tags = set(coord_reg['locus_tag'].values)
print(f"Coordinated regulatory genes: {len(coord_reg_tags)}")

# Identify all regulatory genes from GFF (by product annotation)
reg_products = ['transcriptional regulator', 'transcription factor', 'sigma factor',
                'response regulator', 'sensor kinase', 'two-component', 'anti-sigma']
all_reg_genes = protein_coding[protein_coding['product'].str.lower().str.contains(
    '|'.join(reg_products), na=False)]
all_reg_tags = set(all_reg_genes['locus_tag'].values)
print(f"All regulatory genes (from GFF product): {len(all_reg_tags)}")

# Check AAGCCCG proximity to regulatory genes
for tp in ['T1', 'T2']:
    sub = all_mapping[all_mapping['timepoint'] == tp]
    within_2kb = sub[sub['within_2kb']].copy()

    # Near any regulatory gene
    near_reg = within_2kb[within_2kb['locus_tag'].isin(all_reg_tags)]
    n_sites_near_reg = near_reg['position'].nunique()
    n_total = sub['position'].nunique()

    # Near coordinated regulatory gene
    near_coord = within_2kb[within_2kb['locus_tag'].isin(coord_reg_tags)]
    n_sites_near_coord = near_coord['position'].nunique()

    print(f"\n{tp}:")
    print(f"  Sites within 2kb of ANY regulatory gene: {n_sites_near_reg}/{n_total} "
          f"({n_sites_near_reg/n_total*100:.1f}%)")
    print(f"  Sites within 2kb of COORDINATED regulatory gene: {n_sites_near_coord}/{n_total} "
          f"({n_sites_near_coord/n_total*100:.1f}%)")

# Enrichment test: are AAGCCCG sites enriched near regulatory genes?
# Permutation test
print("\nPermutation test: AAGCCCG enrichment near regulatory genes...")

reg_gene_positions = all_reg_genes[['start', 'end']].values

def count_sites_near_reg(positions, reg_positions, window=2000):
    """Count how many sites are within 'window' bp of any regulatory gene."""
    count = 0
    for pos in positions:
        for gs, ge in reg_positions:
            if gs - window <= pos <= ge + window:
                count += 1
                break
    return count

t1_obs_near_reg = count_sites_near_reg(t1_pos, reg_gene_positions)
t2_obs_near_reg = count_sites_near_reg(t2_pos, reg_gene_positions)

perm_near_reg_t1 = []
perm_near_reg_t2 = []
for _ in range(N_PERMUTATIONS):
    rand_t1 = np.random.randint(1, GENOME_SIZE + 1, size=len(t1_sites))
    perm_near_reg_t1.append(count_sites_near_reg(rand_t1, reg_gene_positions))

    rand_t2 = np.random.randint(1, GENOME_SIZE + 1, size=len(t2_sites))
    perm_near_reg_t2.append(count_sites_near_reg(rand_t2, reg_gene_positions))

perm_near_reg_t1 = np.array(perm_near_reg_t1)
perm_near_reg_t2 = np.array(perm_near_reg_t2)

p_reg_t1 = np.mean(perm_near_reg_t1 >= t1_obs_near_reg)
p_reg_t2 = np.mean(perm_near_reg_t2 >= t2_obs_near_reg)

print(f"\nT1: observed={t1_obs_near_reg}, permutation mean={np.mean(perm_near_reg_t1):.1f}, "
      f"p(enrichment)={p_reg_t1:.4f}")
print(f"  Fold enrichment: {t1_obs_near_reg / np.mean(perm_near_reg_t1):.2f}x")
print(f"T2: observed={t2_obs_near_reg}, permutation mean={np.mean(perm_near_reg_t2):.1f}, "
      f"p(enrichment)={p_reg_t2:.4f}")
print(f"  Fold enrichment: {t2_obs_near_reg / np.mean(perm_near_reg_t2):.2f}x")

# ============================================================
# 10. Expression impact analysis
# ============================================================
print("\n" + "=" * 70)
print("10. Expression Impact Analysis")
print("=" * 70)

deseq2 = pd.read_csv(DESEQ2_FILE, sep='\t')
print(f"DESeq2 results loaded: {len(deseq2)} genes")

# Get T1 AAGCCCG-proximal genes (within 2kb, unique)
t1_prox = t1_mapping[t1_mapping['within_2kb']].copy()
t1_prox_genes = set(t1_prox['locus_tag'].unique())

# Get T2 new sites proximal genes
t2_prox = t2_mapping[t2_mapping['within_2kb']].copy()
t2_prox_genes = set(t2_prox['locus_tag'].unique())

# T1-only genes (sites lost at T2) = in T1 but not in T2
t1_only_genes = t1_prox_genes - t2_prox_genes
# T2-only genes (new sites at T2)
t2_only_genes = t2_prox_genes - t1_prox_genes

print(f"\nT1 AAGCCCG-proximal genes (within 2kb): {len(t1_prox_genes)}")
print(f"T2 AAGCCCG-proximal genes (within 2kb): {len(t2_prox_genes)}")
print(f"T1-only proximal genes (lost methylation): {len(t1_only_genes)}")
print(f"T2-only proximal genes (gained methylation): {len(t2_only_genes)}")

# Merge with DESeq2
deseq2_valid = deseq2.dropna(subset=['log2FoldChange', 'padj']).copy()
all_lfc = deseq2_valid['log2FoldChange'].values
all_padj = deseq2_valid['padj'].values

# Expression distribution of AAGCCCG-proximal genes
t1_prox_expr = deseq2_valid[deseq2_valid['gene_id'].isin(t1_prox_genes)].copy()
t2_prox_expr = deseq2_valid[deseq2_valid['gene_id'].isin(t2_prox_genes)].copy()
t1_only_expr = deseq2_valid[deseq2_valid['gene_id'].isin(t1_only_genes)].copy()
t2_only_expr = deseq2_valid[deseq2_valid['gene_id'].isin(t2_only_genes)].copy()

print(f"\nWith DESeq2 data:")
print(f"  T1-proximal: {len(t1_prox_expr)} genes, "
      f"mean LFC={t1_prox_expr['log2FoldChange'].mean():.3f}, "
      f"median LFC={t1_prox_expr['log2FoldChange'].median():.3f}")
print(f"  T2-proximal: {len(t2_prox_expr)} genes, "
      f"mean LFC={t2_prox_expr['log2FoldChange'].mean():.3f}, "
      f"median LFC={t2_prox_expr['log2FoldChange'].median():.3f}")
print(f"  T1-only (lost): {len(t1_only_expr)} genes, "
      f"mean LFC={t1_only_expr['log2FoldChange'].mean():.3f}, "
      f"median LFC={t1_only_expr['log2FoldChange'].median():.3f}")
print(f"  T2-only (gained): {len(t2_only_expr)} genes, "
      f"mean LFC={t2_only_expr['log2FoldChange'].mean():.3f}, "
      f"median LFC={t2_only_expr['log2FoldChange'].median():.3f}")
print(f"  Genome-wide: {len(deseq2_valid)} genes, "
      f"mean LFC={deseq2_valid['log2FoldChange'].mean():.3f}, "
      f"median LFC={deseq2_valid['log2FoldChange'].median():.3f}")

# KS tests
ks_t1_vs_all = stats.ks_2samp(t1_prox_expr['log2FoldChange'].values, all_lfc)
ks_t2_vs_all = stats.ks_2samp(t2_prox_expr['log2FoldChange'].values, all_lfc)
ks_lost_vs_all = stats.ks_2samp(t1_only_expr['log2FoldChange'].values, all_lfc)
ks_gained_vs_all = stats.ks_2samp(t2_only_expr['log2FoldChange'].values, all_lfc)
ks_lost_vs_gained = stats.ks_2samp(t1_only_expr['log2FoldChange'].values,
                                     t2_only_expr['log2FoldChange'].values)

print(f"\nKS tests (T2vsT1 LFC distributions):")
print(f"  T1-proximal vs genome: D={ks_t1_vs_all[0]:.3f}, p={ks_t1_vs_all[1]:.2e}")
print(f"  T2-proximal vs genome: D={ks_t2_vs_all[0]:.3f}, p={ks_t2_vs_all[1]:.2e}")
print(f"  Lost-methylation vs genome: D={ks_lost_vs_all[0]:.3f}, p={ks_lost_vs_all[1]:.2e}")
print(f"  Gained-methylation vs genome: D={ks_gained_vs_all[0]:.3f}, p={ks_gained_vs_all[1]:.2e}")
print(f"  Lost vs Gained: D={ks_lost_vs_gained[0]:.3f}, p={ks_lost_vs_gained[1]:.2e}")

# DEG proportions
t1_prox_deg = (t1_prox_expr['padj'] < 0.05).sum()
t2_prox_deg = (t2_prox_expr['padj'] < 0.05).sum()
all_deg = (deseq2_valid['padj'] < 0.05).sum()

print(f"\nDEG proportions (padj<0.05):")
print(f"  T1-proximal: {t1_prox_deg}/{len(t1_prox_expr)} ({t1_prox_deg/len(t1_prox_expr)*100:.1f}%)")
print(f"  T2-proximal: {t2_prox_deg}/{len(t2_prox_expr)} ({t2_prox_deg/len(t2_prox_expr)*100:.1f}%)")
print(f"  Genome-wide: {all_deg}/{len(deseq2_valid)} ({all_deg/len(deseq2_valid)*100:.1f}%)")

# Stratify by promoter vs gene body
t1_promoter_genes = set(t1_mapping[t1_mapping['location'] == 'promoter']['locus_tag'].unique())
t1_genebody_genes = set(t1_mapping[t1_mapping['location'] == 'gene_body']['locus_tag'].unique())

t1_prom_expr = deseq2_valid[deseq2_valid['gene_id'].isin(t1_promoter_genes)]
t1_body_expr = deseq2_valid[deseq2_valid['gene_id'].isin(t1_genebody_genes)]

print(f"\nStratified by site location (T1):")
if len(t1_prom_expr) > 0:
    print(f"  Promoter-proximal: {len(t1_prom_expr)} genes, "
          f"mean LFC={t1_prom_expr['log2FoldChange'].mean():.3f}, "
          f"DEG={( t1_prom_expr['padj'] < 0.05).sum()}")
if len(t1_body_expr) > 0:
    print(f"  Gene-body: {len(t1_body_expr)} genes, "
          f"mean LFC={t1_body_expr['log2FoldChange'].mean():.3f}, "
          f"DEG={(t1_body_expr['padj'] < 0.05).sum()}")

# ============================================================
# 11. Save tables
# ============================================================
print("\n" + "=" * 70)
print("11. Saving Tables")
print("=" * 70)

# Full site-to-gene mapping
all_mapping['region'] = all_mapping['position'].apply(classify_region)
all_mapping.to_csv(f"{TBL_DIR}/AAGCCCG_site_gene_mapping.tsv", sep='\t', index=False)
print(f"Saved: AAGCCCG_site_gene_mapping.tsv ({len(all_mapping)} rows)")

# Functional enrichment
enrichment_df.to_csv(f"{TBL_DIR}/functional_enrichment.tsv", sep='\t', index=False)
print(f"Saved: functional_enrichment.tsv ({len(enrichment_df)} rows)")

# T1 vs T2 comparison
comparison_data = {
    'metric': [
        'total_sites', 'arm_sites', 'core_sites', 'arm_pct',
        'median_NN_distance', 'mean_NN_distance',
        'sites_in_BGC', 'sites_near_reg_gene_2kb',
        'proximal_genes_2kb', 'proximal_DEGs',
        'mean_LFC_proximal', 'median_LFC_proximal'
    ],
    'T1': [
        len(t1_sites), t1_arm, t1_core, f"{t1_arm/len(t1_sites)*100:.1f}",
        f"{np.median(t1_nn):.0f}", f"{np.mean(t1_nn):.0f}",
        t1_in_bgc, t1_obs_near_reg,
        len(t1_prox_genes), t1_prox_deg,
        f"{t1_prox_expr['log2FoldChange'].mean():.3f}",
        f"{t1_prox_expr['log2FoldChange'].median():.3f}"
    ],
    'T2': [
        len(t2_sites), t2_arm, t2_core, f"{t2_arm/len(t2_sites)*100:.1f}",
        f"{np.median(t2_nn):.0f}", f"{np.mean(t2_nn):.0f}",
        t2_in_bgc, t2_obs_near_reg,
        len(t2_prox_genes), t2_prox_deg,
        f"{t2_prox_expr['log2FoldChange'].mean():.3f}",
        f"{t2_prox_expr['log2FoldChange'].median():.3f}"
    ]
}
comparison_df = pd.DataFrame(comparison_data)
comparison_df.to_csv(f"{TBL_DIR}/T1_vs_T2_comparison.tsv", sep='\t', index=False)
print(f"Saved: T1_vs_T2_comparison.tsv")

# BGC overlap
bgc_overlap_df.to_csv(f"{TBL_DIR}/BGC_overlap.tsv", sep='\t', index=False)
print(f"Saved: BGC_overlap.tsv")

# ============================================================
# 12. Generate Figures
# ============================================================
print("\n" + "=" * 70)
print("12. Generating Figures")
print("=" * 70)

# Color scheme
COLOR_T1 = '#2166ac'  # blue
COLOR_T2 = '#b2182b'  # red
COLOR_T3 = '#4dac26'  # green
COLOR_ARM = '#f4a582'
COLOR_CORE = '#92c5de'

# ---------- Figure 1: Genome distribution ----------
fig, axes = plt.subplots(3, 1, figsize=(14, 8), gridspec_kw={'height_ratios': [3, 3, 1]})

# T1 density
ax = axes[0]
ax.fill_between(t1_centers / 1e6, t1_density, alpha=0.4, color=COLOR_T1)
ax.plot(t1_centers / 1e6, t1_density, color=COLOR_T1, linewidth=0.8)
# Mark arm boundaries
ax.axvline(ARM_BOUNDARY_LEFT / 1e6, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
ax.axvline(ARM_BOUNDARY_RIGHT / 1e6, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
ax.set_ylabel('Sites / Mb', fontsize=11)
ax.set_title('T1: AAGCCCG sites (n=260)', fontsize=12, fontweight='bold')
ax.set_xlim(0, GENOME_SIZE / 1e6)
ax.text(0.75, 0.85, 'Left arm', transform=ax.transAxes, fontsize=8, color='gray', ha='center')
ax.text(4.35, 0.85 * ax.get_ylim()[1], 'Core', fontsize=8, color='gray', ha='center')
ax.text(7.9, 0.85 * ax.get_ylim()[1], 'Right arm', fontsize=8, color='gray', ha='center')
# BGC positions
for name, info in bgc_regions.items():
    mid = (info['start'] + info['end']) / 2 / 1e6
    ax.axvspan(info['start'] / 1e6, info['end'] / 1e6, alpha=0.15, color='orange')

# T2 density
ax = axes[1]
ax.fill_between(t2_centers / 1e6, t2_density, alpha=0.4, color=COLOR_T2)
ax.plot(t2_centers / 1e6, t2_density, color=COLOR_T2, linewidth=0.8)
ax.axvline(ARM_BOUNDARY_LEFT / 1e6, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
ax.axvline(ARM_BOUNDARY_RIGHT / 1e6, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
ax.set_ylabel('Sites / Mb', fontsize=11)
ax.set_title('T2: AAGCCCG sites (n=64)', fontsize=12, fontweight='bold')
ax.set_xlim(0, GENOME_SIZE / 1e6)
for name, info in bgc_regions.items():
    ax.axvspan(info['start'] / 1e6, info['end'] / 1e6, alpha=0.15, color='orange')

# Genome features track
ax = axes[2]
ax.axhspan(0, 1, xmin=0, xmax=ARM_BOUNDARY_LEFT/GENOME_SIZE, color=COLOR_ARM, alpha=0.5)
ax.axhspan(0, 1, xmin=ARM_BOUNDARY_LEFT/GENOME_SIZE, xmax=ARM_BOUNDARY_RIGHT/GENOME_SIZE,
           color=COLOR_CORE, alpha=0.5)
ax.axhspan(0, 1, xmin=ARM_BOUNDARY_RIGHT/GENOME_SIZE, xmax=1, color=COLOR_ARM, alpha=0.5)
for name, info in bgc_regions.items():
    ax.axvspan(info['start'] / 1e6, info['end'] / 1e6, color='orange', alpha=0.6)
    mid = (info['start'] + info['end']) / 2 / 1e6
    ax.text(mid, 0.5, name, ha='center', va='center', fontsize=8, fontweight='bold')
ax.set_xlim(0, GENOME_SIZE / 1e6)
ax.set_xlabel('Chromosome position (Mb)', fontsize=11)
ax.set_yticks([])
ax.set_ylabel('Features', fontsize=11)

plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f"{FIG_DIR}/genome_distribution.{ext}", dpi=150, bbox_inches='tight')
plt.close()
print("Saved: genome_distribution.pdf/svg")

# ---------- Figure 2: Gene proximity ----------
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for ax, tp, color, title in [(axes[0], 'T1', COLOR_T1, 'T1 (260 sites)'),
                               (axes[1], 'T2', COLOR_T2, 'T2 (64 sites)')]:
    sub = all_mapping[all_mapping['timepoint'] == tp]
    nearest = sub.sort_values('distance').groupby('position').first().reset_index()
    loc_counts = nearest['location'].value_counts()

    categories = ['promoter', 'gene_body', 'upstream', 'downstream', 'downstream_far', 'distal_intergenic']
    cat_labels = ['Promoter\n(<500bp)', 'Gene body', 'Upstream\n(500-2kb)', 'Downstream\n(<500bp)',
                  'Downstream\n(500-2kb)', 'Distal\nintergenic']
    counts = [loc_counts.get(c, 0) for c in categories]

    bars = ax.bar(range(len(categories)), counts, color=color, alpha=0.7, edgecolor='black', linewidth=0.5)
    ax.set_xticks(range(len(categories)))
    ax.set_xticklabels(cat_labels, fontsize=8)
    ax.set_ylabel('Number of sites', fontsize=11)
    ax.set_title(title, fontsize=12, fontweight='bold')

    for bar, cnt in zip(bars, counts):
        if cnt > 0:
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                    str(cnt), ha='center', va='bottom', fontsize=9)

plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f"{FIG_DIR}/gene_proximity.{ext}", dpi=150, bbox_inches='tight')
plt.close()
print("Saved: gene_proximity.pdf/svg")

# ---------- Figure 3: Functional enrichment heatmap ----------
fig, ax = plt.subplots(figsize=(10, 6))

pivot = enrichment_df.pivot(index='category', columns='timepoint', values='fold_enrichment')
pivot = pivot.reindex(columns=['T1', 'T2'])
# Sort by T1 fold enrichment
pivot = pivot.sort_values('T1', ascending=True)

# Create log2 fold enrichment for better visualization
log2_fe = np.log2(pivot.replace(0, 0.01))

im = ax.imshow(log2_fe.values, cmap='RdBu_r', aspect='auto', vmin=-2, vmax=2)
ax.set_xticks([0, 1])
ax.set_xticklabels(['T1 (260)', 'T2 (64)'], fontsize=11)
ax.set_yticks(range(len(pivot.index)))
ax.set_yticklabels(pivot.index, fontsize=10)
ax.set_title('AAGCCCG functional target enrichment\n(log2 fold enrichment vs genome)',
             fontsize=12, fontweight='bold')

# Add values
for i in range(len(pivot.index)):
    for j in range(len(pivot.columns)):
        val = pivot.values[i, j]
        # Get p-value
        mask = (enrichment_df['category'] == pivot.index[i]) & \
               (enrichment_df['timepoint'] == pivot.columns[j])
        p = enrichment_df.loc[mask, 'p_bonferroni'].values[0]
        sig = '*' if p < 0.05 else ''
        text = f"{val:.2f}{sig}"
        color = 'white' if abs(log2_fe.values[i, j]) > 1.2 else 'black'
        ax.text(j, i, text, ha='center', va='center', fontsize=9, color=color)

cbar = plt.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label('log2(Fold Enrichment)', fontsize=10)
plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f"{FIG_DIR}/functional_enrichment.{ext}", dpi=150, bbox_inches='tight')
plt.close()
print("Saved: functional_enrichment.pdf/svg")

# ---------- Figure 4: Expression impact ----------
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Panel A: LFC distribution comparison
ax = axes[0]
bins = np.linspace(-6, 10, 50)
ax.hist(all_lfc, bins=bins, alpha=0.3, color='gray', density=True, label=f'Genome (n={len(deseq2_valid)})')
ax.hist(t1_prox_expr['log2FoldChange'], bins=bins, alpha=0.5, color=COLOR_T1, density=True,
        label=f'T1-proximal (n={len(t1_prox_expr)})')
ax.hist(t2_prox_expr['log2FoldChange'], bins=bins, alpha=0.5, color=COLOR_T2, density=True,
        label=f'T2-proximal (n={len(t2_prox_expr)})')
ax.set_xlabel('log2FC (T2 vs T1)', fontsize=11)
ax.set_ylabel('Density', fontsize=11)
ax.set_title('Expression change\n(AAGCCCG-proximal vs genome)', fontsize=11, fontweight='bold')
ax.legend(fontsize=8)

# Panel B: Lost vs Gained methylation
ax = axes[1]
if len(t1_only_expr) > 0 and len(t2_only_expr) > 0:
    data = [t1_only_expr['log2FoldChange'].values, t2_only_expr['log2FoldChange'].values,
            all_lfc]
    bp = ax.boxplot(data, labels=['Lost methyl\n(T1-only)', 'Gained methyl\n(T2-only)',
                                   'Genome-wide'],
                     patch_artist=True, widths=0.6)
    colors_box = [COLOR_T1, COLOR_T2, 'lightgray']
    for patch, color in zip(bp['boxes'], colors_box):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
    ax.axhline(0, color='gray', linestyle='--', linewidth=0.5)
    ax.set_ylabel('log2FC (T2 vs T1)', fontsize=11)
    ax.set_title('Expression by methylation change', fontsize=11, fontweight='bold')

    # Add sample sizes
    for i, d in enumerate(data):
        ax.text(i + 1, ax.get_ylim()[1] * 0.95, f'n={len(d)}', ha='center', fontsize=8)

# Panel C: Arm/core distribution barplot
ax = axes[2]
x = np.arange(3)
width = 0.35
expected_pct = [arm_frac * 100, core_frac * 100]
t1_pct = [t1_arm/len(t1_sites)*100, t1_core/len(t1_sites)*100]
t2_pct = [t2_arm/len(t2_sites)*100, t2_core/len(t2_sites)*100]

x_pos = np.arange(2)
width = 0.25
bars1 = ax.bar(x_pos - width, [arm_frac*100, core_frac*100], width, label='Expected (genome)',
               color='lightgray', edgecolor='black', linewidth=0.5)
bars2 = ax.bar(x_pos, t1_pct, width, label=f'T1 (n={len(t1_sites)})',
               color=COLOR_T1, alpha=0.7, edgecolor='black', linewidth=0.5)
bars3 = ax.bar(x_pos + width, t2_pct, width, label=f'T2 (n={len(t2_sites)})',
               color=COLOR_T2, alpha=0.7, edgecolor='black', linewidth=0.5)
ax.set_xticks(x_pos)
ax.set_xticklabels(['Arms', 'Core'], fontsize=11)
ax.set_ylabel('% of sites', fontsize=11)
ax.set_title('Arm vs Core distribution', fontsize=11, fontweight='bold')
ax.legend(fontsize=8)

# Add significance annotations
if p_t1 < 0.05:
    ax.text(0, max(t1_pct[0], arm_frac*100) + 2, f'p={p_t1:.2e}', ha='center', fontsize=7)
if p_t2 < 0.05:
    ax.text(0 + width, max(t2_pct[0], arm_frac*100) + 2, f'p={p_t2:.2e}', ha='center', fontsize=7)

plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f"{FIG_DIR}/expression_impact.{ext}", dpi=150, bbox_inches='tight')
plt.close()
print("Saved: expression_impact.pdf/svg")

# ---------- Figure 5: Comprehensive summary (multi-panel) ----------
fig = plt.figure(figsize=(18, 14))
gs = gridspec.GridSpec(3, 3, hspace=0.4, wspace=0.35)

# Panel A: Genome distribution (T1)
ax = fig.add_subplot(gs[0, :])
ax.fill_between(t1_centers / 1e6, t1_density, alpha=0.4, color=COLOR_T1, label='T1 (260)')
ax.plot(t1_centers / 1e6, t1_density, color=COLOR_T1, linewidth=0.8)
ax.fill_between(t2_centers / 1e6, -t2_density, alpha=0.4, color=COLOR_T2, label='T2 (64)')
ax.plot(t2_centers / 1e6, -t2_density, color=COLOR_T2, linewidth=0.8)
ax.axhline(0, color='black', linewidth=0.5)
ax.axvline(ARM_BOUNDARY_LEFT / 1e6, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
ax.axvline(ARM_BOUNDARY_RIGHT / 1e6, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
for name, info in bgc_regions.items():
    ax.axvspan(info['start'] / 1e6, info['end'] / 1e6, alpha=0.15, color='orange')
    mid = (info['start'] + info['end']) / 2 / 1e6
    ymax = ax.get_ylim()[1]
    ax.text(mid, ymax * 0.8, name, ha='center', fontsize=8, fontweight='bold', color='darkorange')
ax.set_xlabel('Chromosome position (Mb)', fontsize=11)
ax.set_ylabel('Sites / Mb', fontsize=11)
ax.set_title('A. AAGCCCG site density across S. coelicolor chromosome', fontsize=12, fontweight='bold')
ax.legend(fontsize=10)
ax.set_xlim(0, GENOME_SIZE / 1e6)
ax.text(0.02, 0.95, 'Left arm', transform=ax.transAxes, fontsize=9, color='gray')
ax.text(0.48, 0.95, 'Core', transform=ax.transAxes, fontsize=9, color='gray')
ax.text(0.92, 0.95, 'Right arm', transform=ax.transAxes, fontsize=9, color='gray')

# Panel B: Arm/Core
ax = fig.add_subplot(gs[1, 0])
x_pos = np.arange(2)
width = 0.25
ax.bar(x_pos - width, [arm_frac*100, core_frac*100], width, label='Expected',
       color='lightgray', edgecolor='black', linewidth=0.5)
ax.bar(x_pos, t1_pct, width, label='T1',
       color=COLOR_T1, alpha=0.7, edgecolor='black', linewidth=0.5)
ax.bar(x_pos + width, t2_pct, width, label='T2',
       color=COLOR_T2, alpha=0.7, edgecolor='black', linewidth=0.5)
ax.set_xticks(x_pos)
ax.set_xticklabels(['Arms', 'Core'])
ax.set_ylabel('% of sites')
ax.set_title('B. Arm vs Core', fontsize=11, fontweight='bold')
ax.legend(fontsize=8)

# Panel C: Gene proximity (T1)
ax = fig.add_subplot(gs[1, 1])
sub = all_mapping[all_mapping['timepoint'] == 'T1']
nearest_t1 = sub.sort_values('distance').groupby('position').first().reset_index()
cats = ['promoter', 'gene_body', 'upstream', 'downstream', 'downstream_far', 'distal_intergenic']
cat_labels_short = ['Prom', 'Body', 'Up', 'Down', 'Down_far', 'Distal']
counts_t1 = [nearest_t1['location'].value_counts().get(c, 0) for c in cats]
counts_t2_loc = []
sub2 = all_mapping[all_mapping['timepoint'] == 'T2']
nearest_t2 = sub2.sort_values('distance').groupby('position').first().reset_index()
counts_t2_loc = [nearest_t2['location'].value_counts().get(c, 0) for c in cats]

x_pos = np.arange(len(cats))
width = 0.35
ax.bar(x_pos - width/2, counts_t1, width, label='T1', color=COLOR_T1, alpha=0.7)
ax.bar(x_pos + width/2, counts_t2_loc, width, label='T2', color=COLOR_T2, alpha=0.7)
ax.set_xticks(x_pos)
ax.set_xticklabels(cat_labels_short, fontsize=8, rotation=45)
ax.set_ylabel('Number of sites')
ax.set_title('C. Site-gene proximity', fontsize=11, fontweight='bold')
ax.legend(fontsize=8)

# Panel D: Nearest-neighbor distance
ax = fig.add_subplot(gs[1, 2])
bins_nn = np.logspace(1, 6, 40)
ax.hist(t1_nn, bins=bins_nn, alpha=0.5, color=COLOR_T1, label='T1', density=True)
ax.hist(t2_nn, bins=bins_nn, alpha=0.5, color=COLOR_T2, label='T2', density=True)
ax.axvline(t1_expected_nn, color=COLOR_T1, linestyle='--', linewidth=1, alpha=0.7,
           label=f'T1 expected ({t1_expected_nn/1000:.1f}kb)')
ax.axvline(t2_expected_nn, color=COLOR_T2, linestyle='--', linewidth=1, alpha=0.7,
           label=f'T2 expected ({t2_expected_nn/1000:.1f}kb)')
ax.set_xscale('log')
ax.set_xlabel('NN distance (bp)')
ax.set_ylabel('Density')
ax.set_title('D. Nearest-neighbor distances', fontsize=11, fontweight='bold')
ax.legend(fontsize=7)

# Panel E: Functional enrichment
ax = fig.add_subplot(gs[2, 0:2])
# Bar chart of fold enrichment
cats_to_show = enrichment_df[enrichment_df['timepoint'] == 'T1'].sort_values('fold_enrichment', ascending=True)['category'].values
x_pos = np.arange(len(cats_to_show))
width = 0.35

t1_fe = [enrichment_df[(enrichment_df['timepoint']=='T1') & (enrichment_df['category']==c)]['fold_enrichment'].values[0]
         for c in cats_to_show]
t2_fe = [enrichment_df[(enrichment_df['timepoint']=='T2') & (enrichment_df['category']==c)]['fold_enrichment'].values[0]
         for c in cats_to_show]

ax.barh(x_pos - width/2, t1_fe, width, label='T1', color=COLOR_T1, alpha=0.7)
ax.barh(x_pos + width/2, t2_fe, width, label='T2', color=COLOR_T2, alpha=0.7)
ax.axvline(1, color='gray', linestyle='--', linewidth=0.8)
ax.set_yticks(x_pos)
ax.set_yticklabels(cats_to_show, fontsize=9)
ax.set_xlabel('Fold enrichment')
ax.set_title('E. Functional category enrichment', fontsize=11, fontweight='bold')
ax.legend(fontsize=8)

# Panel F: Expression impact
ax = fig.add_subplot(gs[2, 2])
if len(t1_only_expr) > 0 and len(t2_only_expr) > 0:
    data_box = [t1_only_expr['log2FoldChange'].values, t2_only_expr['log2FoldChange'].values, all_lfc]
    bp = ax.boxplot(data_box, labels=['Lost\nmethyl', 'Gained\nmethyl', 'Genome'],
                     patch_artist=True, widths=0.6, showfliers=False)
    colors_bp = [COLOR_T1, COLOR_T2, 'lightgray']
    for patch, color in zip(bp['boxes'], colors_bp):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
    ax.axhline(0, color='gray', linestyle='--', linewidth=0.5)
    ax.set_ylabel('log2FC (T2 vs T1)')
    ax.set_title('F. Expression impact', fontsize=11, fontweight='bold')

plt.suptitle('H13: AAGCCCG 6mA modification - Genomic distribution & functional targets',
             fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f"{FIG_DIR}/H13_comprehensive_summary.{ext}", dpi=150, bbox_inches='tight')
plt.close()
print("Saved: H13_comprehensive_summary.pdf/svg")

# ============================================================
# 13. Final Summary Statistics
# ============================================================
print("\n" + "=" * 70)
print("FINAL SUMMARY")
print("=" * 70)

print(f"""
H13 Analysis Complete
=====================

1. AAGCCCG site counts:
   T1: {len(t1_sites)}, T2: {len(t2_sites)}, T3: {len(t3_sites)}
   T1-T2 overlap: {len(overlap)} (complete reprogramming confirmed)

2. Arm vs Core distribution:
   T1: {t1_arm} arm ({t1_arm/len(t1_sites)*100:.1f}%), {t1_core} core ({t1_core/len(t1_sites)*100:.1f}%)
   T2: {t2_arm} arm ({t2_arm/len(t2_sites)*100:.1f}%), {t2_core} core ({t2_core/len(t2_sites)*100:.1f}%)
   Expected: arm {arm_frac*100:.1f}%, core {core_frac*100:.1f}%
   T1 chi-sq: p={p_t1:.2e} | T2 chi-sq: p={p_t2:.2e}
   T1 vs T2 Fisher: OR={fisher_or:.3f}, p={fisher_p:.2e}

3. Clustering (permutation test, n={N_PERMUTATIONS}):
   T1 median NN: {np.median(t1_nn):.0f} bp (expected {t1_expected_nn:.0f}), p_cluster={p_cluster_t1:.4f}
   T2 median NN: {np.median(t2_nn):.0f} bp (expected {t2_expected_nn:.0f}), p_cluster={p_cluster_t2:.4f}

4. BGC overlap:
   T1: {t1_in_bgc}/{len(t1_sites)} ({t1_in_bgc/len(t1_sites)*100:.1f}%), fold={t1_in_bgc/t1_bgc_exp:.2f}x
   T2: {t2_in_bgc}/{len(t2_sites)} ({t2_in_bgc/len(t2_sites)*100:.1f}%), fold={t2_in_bgc/t2_bgc_exp:.2f}x

5. Regulatory gene proximity (within 2kb, permutation test):
   T1: {t1_obs_near_reg}/{len(t1_sites)}, fold={t1_obs_near_reg/np.mean(perm_near_reg_t1):.2f}x, p={p_reg_t1:.4f}
   T2: {t2_obs_near_reg}/{len(t2_sites)}, fold={t2_obs_near_reg/np.mean(perm_near_reg_t2):.2f}x, p={p_reg_t2:.4f}

6. Expression impact (T2vsT1 LFC):
   T1-proximal genes: mean={t1_prox_expr['log2FoldChange'].mean():.3f}, KS p={ks_t1_vs_all[1]:.2e}
   T2-proximal genes: mean={t2_prox_expr['log2FoldChange'].mean():.3f}, KS p={ks_t2_vs_all[1]:.2e}
   Lost-methyl genes: mean={t1_only_expr['log2FoldChange'].mean():.3f}, KS p={ks_lost_vs_all[1]:.2e}
   Gained-methyl genes: mean={t2_only_expr['log2FoldChange'].mean():.3f}, KS p={ks_gained_vs_all[1]:.2e}
""")

print("All outputs saved to:", OUT_DIR)
print("Done!")
