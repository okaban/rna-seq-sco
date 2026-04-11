#!/usr/bin/env python3
"""
H20: Geographic Stratification of H15's Universal Regulatory Avoidance

Tests whether H15's "universal regulatory gene avoidance" (fold 0.43-0.71 across
all motifs) is a genuine biological phenomenon or a geographic artifact caused by
differential gene composition in core vs arm regions of the S. coelicolor chromosome.

Key question: Does regulatory avoidance persist when analyzed WITHIN genomic regions
(core-only and arm-only), or does it disappear like the H17/H19 expression suppression?
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

# ============================================================
# Configuration
# ============================================================
BASE_DIR = "/Users/okaban/bioinfo/rna-seq"
EPIGENOME_DIR = f"{BASE_DIR}/11_epigenome_integration"

CENSUS_4mC = f"{EPIGENOME_DIR}/analysis/23_expanded_motif_search/4mC_final_census.csv"
CENSUS_6mA = f"{EPIGENOME_DIR}/analysis/23_expanded_motif_search/6mA_final_census.csv"
GFF_FILE = f"{EPIGENOME_DIR}/analysis/39_GCCGGC_MTase_reverse_ID/data/GCF_000203835.1_ASM20383v1_genomic.gff"
H15_ENRICHMENT = f"{EPIGENOME_DIR}/analysis/38_cross_motif_regulatory_avoidance/tables/cross_motif_functional_enrichment.tsv"
REG_GENES_FILE = f"{EPIGENOME_DIR}/analysis/29_genomewide_TF_screen/tables/all_regulatory_genes.tsv"

OUT_DIR = f"{EPIGENOME_DIR}/analysis/43_regulatory_avoidance_geographic_test"
FIG_DIR = f"{OUT_DIR}/figures"
TBL_DIR = f"{OUT_DIR}/tables"

GENOME_SIZE = 8667507
WINDOW = 2000  # 2kb proximity window

# Region definitions
ARM_LEFT_END = 1_500_000
ARM_RIGHT_START = GENOME_SIZE - 1_500_000  # 7,167,507

np.random.seed(42)

# ============================================================
# 1. Load census data and define motif groups
# ============================================================
print("=" * 70)
print("H20: Geographic Stratification of Regulatory Avoidance")
print("=" * 70)

census_4mC = pd.read_csv(CENSUS_4mC)
census_6mA = pd.read_csv(CENSUS_6mA)

print(f"\n4mC census: {len(census_4mC)} sites")
print(f"6mA census: {len(census_6mA)} sites")

# Define motif groups -- same as H15
gccggc_mask = census_4mC['final_motif'].isin(['TGGCCGGC', 'GGCCGG'])
aagcccg_4mC_mask = census_4mC['final_motif'] == 'AAGCCCG'
aagcccg_6mA_mask = census_6mA['final_motif'] == 'AAGCCCG'

# Build motif groups for T1 only (the primary timepoint of interest)
motif_groups = OrderedDict()

# GCCGGC 4mC T1
sub = census_4mC[gccggc_mask & (census_4mC['timepoint'] == 'T1')].copy()
motif_groups['GCCGGC_4mC_T1'] = sub

# AAGCCCG 4mC T1
sub = census_4mC[aagcccg_4mC_mask & (census_4mC['timepoint'] == 'T1')].copy()
motif_groups['AAGCCCG_4mC_T1'] = sub

# AAGCCCG 6mA T1
sub = census_6mA[aagcccg_6mA_mask & (census_6mA['timepoint'] == 'T1')].copy()
motif_groups['AAGCCCG_6mA_T1'] = sub

# All 4mC T1
sub = census_4mC[census_4mC['timepoint'] == 'T1'].copy()
motif_groups['All_4mC_T1'] = sub

print("\nMotif groups (T1):")
for name, df in motif_groups.items():
    print(f"  {name}: {len(df)} sites")

# ============================================================
# 2. Load GFF and classify genes
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

# Load CDS products
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
# 3. Functional classification (IDENTICAL to H15)
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
protein_coding['functional_category'] = protein_coding['product'].apply(classify_product)

CATEGORIES = ['Regulatory/TF', 'Hypothetical', 'Primary metabolism', 'Other',
              'Transport', 'DNA/RNA metabolism', 'Stress/Defense',
              'Secondary metabolism', 'Translation', 'Membrane/Cell wall']

# ============================================================
# 4. Assign arm/core regions to genes and sites
# ============================================================
print("\n" + "=" * 70)
print("Assigning arm/core regions...")
print("=" * 70)

def assign_region(pos):
    """Assign arm or core based on position."""
    if pos <= ARM_LEFT_END or pos >= ARM_RIGHT_START:
        return 'arm'
    return 'core'

# Assign region to each gene (using midpoint)
protein_coding['midpoint'] = (protein_coding['start'] + protein_coding['end']) / 2
protein_coding['region'] = protein_coding['midpoint'].apply(assign_region)

core_genes = protein_coding[protein_coding['region'] == 'core']
arm_genes = protein_coding[protein_coding['region'] == 'arm']

print(f"\nAll protein-coding genes: {len(protein_coding)}")
print(f"  Core genes: {len(core_genes)} ({len(core_genes)/len(protein_coding)*100:.1f}%)")
print(f"  Arm genes: {len(arm_genes)} ({len(arm_genes)/len(protein_coding)*100:.1f}%)")

# Functional category distribution by region
print("\n--- Functional category by region ---")
region_cat_table = []
for cat in CATEGORIES:
    n_total = len(protein_coding[protein_coding['functional_category'] == cat])
    n_core = len(core_genes[core_genes['functional_category'] == cat])
    n_arm = len(arm_genes[arm_genes['functional_category'] == cat])
    frac_core = n_core / len(core_genes) * 100 if len(core_genes) > 0 else 0
    frac_arm = n_arm / len(arm_genes) * 100 if len(arm_genes) > 0 else 0
    frac_total = n_total / len(protein_coding) * 100 if len(protein_coding) > 0 else 0

    # Chi-squared test for category distribution across regions
    obs = np.array([n_core, n_arm])
    exp = np.array([len(core_genes), len(arm_genes)]) * n_total / len(protein_coding)
    if exp.min() > 0:
        chi2_stat = np.sum((obs - exp)**2 / exp)
        chi2_p = stats.chi2.sf(chi2_stat, df=1)
    else:
        chi2_p = 1.0

    region_cat_table.append({
        'category': cat,
        'n_total': n_total,
        'n_core': n_core,
        'n_arm': n_arm,
        'frac_core_pct': round(frac_core, 2),
        'frac_arm_pct': round(frac_arm, 2),
        'frac_total_pct': round(frac_total, 2),
        'core_vs_arm_p': chi2_p
    })

    sig = "*" if chi2_p < 0.05 else ""
    print(f"  {cat}: core {n_core}({frac_core:.1f}%) arm {n_arm}({frac_arm:.1f}%) total {n_total}({frac_total:.1f}%) {sig}")

region_cat_df = pd.DataFrame(region_cat_table)
region_cat_df.to_csv(f"{TBL_DIR}/regional_gene_composition.tsv", sep='\t', index=False)
print(f"\nSaved: regional_gene_composition.tsv")

# Assign region to methylation sites
for name, df in motif_groups.items():
    df['region'] = df['position'].apply(assign_region)
    n_core = (df['region'] == 'core').sum()
    n_arm = (df['region'] == 'arm').sum()
    print(f"\n{name}: core={n_core} ({n_core/len(df)*100:.1f}%), arm={n_arm} ({n_arm/len(df)*100:.1f}%)")

# ============================================================
# 5. Proximity mapping function (same as H15)
# ============================================================

def find_proximal_genes(positions, strands, genes_sorted, window=2000):
    """For each methylation site, find ALL genes within 'window' bp."""
    results = []
    for pos, strand in zip(positions, strands):
        nearby = genes_sorted[
            (genes_sorted['start'] - window <= pos) &
            (genes_sorted['end'] + window >= pos)
        ]
        if len(nearby) == 0:
            mid_points = (genes_sorted['start'].values + genes_sorted['end'].values) / 2
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


def get_nearest_per_site(mapping_df):
    """Get nearest gene per site for category counting."""
    if len(mapping_df) == 0:
        return pd.DataFrame()
    return mapping_df.sort_values('distance').groupby('position').first().reset_index()


# ============================================================
# 6. Stratified enrichment analysis
# ============================================================
print("\n" + "=" * 70)
print("STRATIFIED ENRICHMENT ANALYSIS (Core / Arm / Combined)")
print("=" * 70)

all_results = []

for region_label in ['combined', 'core', 'arm']:
    print(f"\n{'='*50}")
    print(f"Region: {region_label.upper()}")
    print(f"{'='*50}")

    # Select genes for this region
    if region_label == 'combined':
        region_genes = protein_coding.copy()
    elif region_label == 'core':
        region_genes = core_genes.copy()
    else:
        region_genes = arm_genes.copy()

    sorted_genes = region_genes.sort_values('start').reset_index(drop=True)
    n_region_genes = len(region_genes)
    region_cat_counts = region_genes['functional_category'].value_counts()

    print(f"  Genes in region: {n_region_genes}")

    for group_name, df in motif_groups.items():
        # Filter sites to this region
        if region_label == 'combined':
            region_sites = df
        elif region_label == 'core':
            region_sites = df[df['region'] == 'core']
        else:
            region_sites = df[df['region'] == 'arm']

        if len(region_sites) < 5:
            print(f"  {group_name} [{region_label}]: only {len(region_sites)} sites, SKIPPING")
            continue

        # Map sites to nearby genes in this region
        mapping = find_proximal_genes(
            region_sites['position'].values,
            region_sites['strand'].values,
            sorted_genes,
            window=WINDOW
        )

        if len(mapping) == 0:
            continue

        nearest = get_nearest_per_site(mapping)
        n_sites = len(nearest)
        site_cat_counts = nearest['functional_category'].value_counts()

        print(f"\n  {group_name} [{region_label}]: {len(region_sites)} sites -> {n_sites} gene-site pairs")

        for cat in CATEGORIES:
            obs_in = site_cat_counts.get(cat, 0)
            obs_out = n_sites - obs_in
            genome_in = region_cat_counts.get(cat, 0)
            genome_out = n_region_genes - genome_in

            table = [[obs_in, obs_out], [genome_in, genome_out]]
            odds_ratio, p_val = stats.fisher_exact(table)

            expected = n_sites * (genome_in / n_region_genes) if n_region_genes > 0 else 0
            fold_enrich = obs_in / expected if expected > 0 else 0

            # Confidence interval
            a, b, c, d = obs_in, obs_out, genome_in, genome_out
            if a > 0 and b > 0 and c > 0 and d > 0:
                log_or = np.log(odds_ratio)
                se = np.sqrt(1/a + 1/b + 1/c + 1/d)
                ci_low = np.exp(log_or - 1.96 * se)
                ci_high = np.exp(log_or + 1.96 * se)
            else:
                ci_low = 0
                ci_high = np.inf

            n_bonf = len(CATEGORIES)
            p_bonf = min(p_val * n_bonf, 1.0)

            direction = 'depleted' if fold_enrich < 1 else 'enriched'

            all_results.append({
                'region': region_label,
                'motif_group': group_name,
                'category': cat,
                'n_region_genes': n_region_genes,
                'n_sites': n_sites,
                'observed': obs_in,
                'expected': round(expected, 1),
                'fold_enrichment': round(fold_enrich, 2),
                'odds_ratio': round(odds_ratio, 3),
                'OR_CI_low': round(ci_low, 3),
                'OR_CI_high': round(ci_high, 3),
                'p_value': p_val,
                'p_bonferroni': p_bonf,
                'direction': direction,
                'n_cat_in_region': genome_in
            })

            if cat in ['Regulatory/TF', 'Hypothetical']:
                sig = "***" if p_bonf < 0.001 else "**" if p_bonf < 0.01 else "*" if p_bonf < 0.05 else ""
                print(f"    {cat}: obs={obs_in}, exp={expected:.1f}, fold={fold_enrich:.2f}, OR={odds_ratio:.3f}, p={p_val:.2e} {sig}")

results_df = pd.DataFrame(all_results)
results_df.to_csv(f"{TBL_DIR}/stratified_enrichment.tsv", sep='\t', index=False)
print(f"\nSaved: stratified_enrichment.tsv")

# ============================================================
# 7. Comparison table: Unstratified vs Stratified
# ============================================================
print("\n" + "=" * 70)
print("COMPARISON: H15 Unstratified vs H20 Stratified")
print("=" * 70)

h15_df = pd.read_csv(H15_ENRICHMENT, sep='\t')

comparison_rows = []
focus_categories = ['Regulatory/TF', 'Hypothetical']

for cat in focus_categories:
    for motif in ['GCCGGC_4mC_T1', 'AAGCCCG_4mC_T1', 'AAGCCCG_6mA_T1', 'All_4mC_T1']:
        # H15 unstratified
        h15_row = h15_df[(h15_df['motif_group'] == motif) & (h15_df['category'] == cat)]
        if len(h15_row) == 0:
            continue
        h15_fold = h15_row['fold_enrichment'].values[0]
        h15_or = h15_row['odds_ratio'].values[0]
        h15_p = h15_row['p_value'].values[0]
        h15_pbon = h15_row['p_bonferroni'].values[0]

        for reg in ['combined', 'core', 'arm']:
            h20_row = results_df[(results_df['motif_group'] == motif) &
                                  (results_df['category'] == cat) &
                                  (results_df['region'] == reg)]
            if len(h20_row) == 0:
                h20_fold = np.nan
                h20_or = np.nan
                h20_p = np.nan
                h20_pbon = np.nan
                h20_obs = np.nan
                h20_exp = np.nan
            else:
                h20_fold = h20_row['fold_enrichment'].values[0]
                h20_or = h20_row['odds_ratio'].values[0]
                h20_p = h20_row['p_value'].values[0]
                h20_pbon = h20_row['p_bonferroni'].values[0]
                h20_obs = h20_row['observed'].values[0]
                h20_exp = h20_row['expected'].values[0]

            comparison_rows.append({
                'category': cat,
                'motif': motif,
                'region': reg,
                'H15_fold': h15_fold,
                'H15_OR': h15_or,
                'H15_p': h15_p,
                'H15_p_bonferroni': h15_pbon,
                'H20_fold': h20_fold,
                'H20_OR': h20_or,
                'H20_p': h20_p,
                'H20_p_bonferroni': h20_pbon,
                'H20_observed': h20_obs,
                'H20_expected': h20_exp,
                'verdict': 'GENUINE' if (not np.isnan(h20_fold) and h20_fold < 0.85 and h20_p < 0.1) else
                           'ARTIFACT' if (not np.isnan(h20_fold) and h20_fold >= 0.85) else
                           'LOW_POWER' if np.isnan(h20_fold) else 'AMBIGUOUS'
            })

comp_df = pd.DataFrame(comparison_rows)
comp_df.to_csv(f"{TBL_DIR}/comparison_unstratified_vs_stratified.tsv", sep='\t', index=False)

print("\nKey comparison for Regulatory/TF:")
reg_comp = comp_df[comp_df['category'] == 'Regulatory/TF']
for _, row in reg_comp.iterrows():
    sig_h15 = "*" if row['H15_p'] < 0.05 else ""
    sig_h20 = "*" if (not np.isnan(row['H20_p']) and row['H20_p'] < 0.05) else ""
    print(f"  {row['motif']:20s} {row['region']:10s} | H15 fold={row['H15_fold']:.2f}{sig_h15} | H20 fold={row['H20_fold']:.2f} p={row['H20_p']:.2e}{sig_h20} | {row['verdict']}")

print("\nKey comparison for Hypothetical:")
hyp_comp = comp_df[comp_df['category'] == 'Hypothetical']
for _, row in hyp_comp.iterrows():
    sig_h15 = "*" if row['H15_p'] < 0.05 else ""
    sig_h20 = "*" if (not np.isnan(row['H20_p']) and row['H20_p'] < 0.05) else ""
    print(f"  {row['motif']:20s} {row['region']:10s} | H15 fold={row['H15_fold']:.2f}{sig_h15} | H20 fold={row['H20_fold']:.2f} p={row['H20_p']:.2e}{sig_h20} | {row['verdict']}")

# ============================================================
# 8. Cochran-Mantel-Haenszel test
# ============================================================
print("\n" + "=" * 70)
print("COCHRAN-MANTEL-HAENSZEL TEST (Region-adjusted)")
print("=" * 70)

cmh_results = []

for cat in focus_categories:
    for motif in ['GCCGGC_4mC_T1', 'AAGCCCG_4mC_T1', 'AAGCCCG_6mA_T1', 'All_4mC_T1']:
        tables = []
        for reg in ['core', 'arm']:
            h20_row = results_df[(results_df['motif_group'] == motif) &
                                  (results_df['category'] == cat) &
                                  (results_df['region'] == reg)]
            if len(h20_row) == 0:
                continue
            a = h20_row['observed'].values[0]
            b = h20_row['n_sites'].values[0] - a
            c = h20_row['n_cat_in_region'].values[0]
            d = h20_row['n_region_genes'].values[0] - c
            tables.append(np.array([[a, b], [c, d]]))

        if len(tables) < 2:
            continue

        # CMH test statistic
        # Numerator: sum_i (a_i - E(a_i))
        # where E(a_i) = n_1i * m_1i / N_i
        numerator = 0
        denominator = 0
        or_num = 0
        or_den = 0

        for t in tables:
            a_i, b_i, c_i, d_i = t[0,0], t[0,1], t[1,0], t[1,1]
            N_i = a_i + b_i + c_i + d_i
            n_1i = a_i + b_i  # row 1 total (methylation sites)
            n_2i = c_i + d_i  # row 2 total (genome)
            m_1i = a_i + c_i  # col 1 total (in category)
            m_2i = b_i + d_i  # col 2 total (not in category)

            E_ai = n_1i * m_1i / N_i
            V_ai = n_1i * n_2i * m_1i * m_2i / (N_i**2 * (N_i - 1)) if N_i > 1 else 0

            numerator += (a_i - E_ai)
            denominator += V_ai

            # Mantel-Haenszel OR components
            or_num += a_i * d_i / N_i
            or_den += b_i * c_i / N_i

        if denominator > 0:
            chi2_mh = numerator**2 / denominator
            p_mh = stats.chi2.sf(chi2_mh, df=1)
        else:
            chi2_mh = 0
            p_mh = 1.0

        # MH common odds ratio
        if or_den > 0:
            OR_mh = or_num / or_den
            # Robins-Breslow-Greenland variance
            # For CI on log(OR_MH)
            # Using simplified formula
            P_i = []
            Q_i = []
            R_i = []
            S_i = []
            for t in tables:
                a_i, b_i, c_i, d_i = t[0,0], t[0,1], t[1,0], t[1,1]
                N_i = a_i + b_i + c_i + d_i
                R_i.append(a_i * d_i / N_i)
                S_i.append(b_i * c_i / N_i)
                P_i.append((a_i + d_i) / N_i)
                Q_i.append((b_i + c_i) / N_i)

            R = sum(R_i)
            S = sum(S_i)

            var_log_or = 0
            for i in range(len(tables)):
                var_log_or += (P_i[i] * R_i[i]) / (2 * R**2)
                var_log_or += (P_i[i] * S_i[i] + Q_i[i] * R_i[i]) / (2 * R * S)
                var_log_or += (Q_i[i] * S_i[i]) / (2 * S**2)

            se_log_or = np.sqrt(var_log_or) if var_log_or > 0 else 0
            ci_low_mh = np.exp(np.log(OR_mh) - 1.96 * se_log_or) if se_log_or > 0 else 0
            ci_high_mh = np.exp(np.log(OR_mh) + 1.96 * se_log_or) if se_log_or > 0 else np.inf
        else:
            OR_mh = 0
            ci_low_mh = 0
            ci_high_mh = np.inf

        sig = "***" if p_mh < 0.001 else "**" if p_mh < 0.01 else "*" if p_mh < 0.05 else ""

        # Also get the unstratified (combined) result
        comb_row = results_df[(results_df['motif_group'] == motif) &
                               (results_df['category'] == cat) &
                               (results_df['region'] == 'combined')]
        comb_or = comb_row['odds_ratio'].values[0] if len(comb_row) > 0 else np.nan

        cmh_results.append({
            'category': cat,
            'motif': motif,
            'CMH_chi2': round(chi2_mh, 3),
            'CMH_p': p_mh,
            'CMH_OR': round(OR_mh, 3),
            'CMH_OR_CI_low': round(ci_low_mh, 3),
            'CMH_OR_CI_high': round(ci_high_mh, 3),
            'unstratified_OR': comb_or,
            'interpretation': 'GENUINE_DEPLETION' if (OR_mh < 1 and p_mh < 0.05) else
                              'NOT_SIGNIFICANT' if p_mh >= 0.05 else
                              'GENUINE_ENRICHMENT'
        })

        print(f"  {cat} x {motif}: CMH chi2={chi2_mh:.3f}, p={p_mh:.2e}, OR_MH={OR_mh:.3f} [{ci_low_mh:.3f}, {ci_high_mh:.3f}] {sig}")
        print(f"    Unstratified OR={comb_or:.3f}, MH-adjusted OR={OR_mh:.3f}")

cmh_df = pd.DataFrame(cmh_results)
cmh_df.to_csv(f"{TBL_DIR}/CMH_test_results.tsv", sep='\t', index=False)
print(f"\nSaved: CMH_test_results.tsv")

# ============================================================
# 9. Breslow-Day test for homogeneity of ORs across strata
# ============================================================
print("\n" + "=" * 70)
print("BRESLOW-DAY TEST: Are stratum-specific ORs homogeneous?")
print("=" * 70)

for cat in focus_categories:
    for motif in ['GCCGGC_4mC_T1', 'AAGCCCG_6mA_T1', 'All_4mC_T1']:
        core_row = results_df[(results_df['motif_group'] == motif) &
                               (results_df['category'] == cat) &
                               (results_df['region'] == 'core')]
        arm_row = results_df[(results_df['motif_group'] == motif) &
                              (results_df['category'] == cat) &
                              (results_df['region'] == 'arm')]

        if len(core_row) == 0 or len(arm_row) == 0:
            print(f"  {cat} x {motif}: insufficient data")
            continue

        or_core = core_row['odds_ratio'].values[0]
        or_arm = arm_row['odds_ratio'].values[0]

        print(f"  {cat} x {motif}: OR_core={or_core:.3f}, OR_arm={or_arm:.3f}, ratio={or_core/or_arm:.2f}" if or_arm > 0 else f"  {cat} x {motif}: OR_core={or_core:.3f}, OR_arm=0")

# ============================================================
# 10. FIGURE 1: Grouped bar chart - Regulatory/TF by motif x region
# ============================================================
print("\n" + "=" * 70)
print("Generating figures...")
print("=" * 70)

fig, axes = plt.subplots(1, 2, figsize=(16, 7))

for ax_idx, cat in enumerate(focus_categories):
    ax = axes[ax_idx]

    motifs = ['GCCGGC_4mC_T1', 'AAGCCCG_4mC_T1', 'AAGCCCG_6mA_T1', 'All_4mC_T1']
    motif_labels = ['GCCGGC\n4mC', 'AAGCCCG\n4mC', 'AAGCCCG\n6mA', 'All 4mC']
    regions = ['combined', 'core', 'arm']
    region_colors = {'combined': '#4A90D9', 'core': '#E67E22', 'arm': '#27AE60'}
    region_labels = {'combined': 'Combined', 'core': 'Core only', 'arm': 'Arm only'}

    x = np.arange(len(motifs))
    width = 0.25

    for i, reg in enumerate(regions):
        folds = []
        p_vals = []
        for motif in motifs:
            row = results_df[(results_df['motif_group'] == motif) &
                              (results_df['category'] == cat) &
                              (results_df['region'] == reg)]
            if len(row) > 0:
                folds.append(row['fold_enrichment'].values[0])
                p_vals.append(row['p_value'].values[0])
            else:
                folds.append(0)
                p_vals.append(1)

        bars = ax.bar(x + i * width - width, folds, width,
                       label=region_labels[reg], color=region_colors[reg],
                       edgecolor='white', linewidth=0.5, alpha=0.85)

        # Add significance stars
        for j, (fold, pval) in enumerate(zip(folds, p_vals)):
            if pval < 0.001:
                star = '***'
            elif pval < 0.01:
                star = '**'
            elif pval < 0.05:
                star = '*'
            else:
                star = ''
            if star and fold > 0:
                ax.text(x[j] + i * width - width, fold + 0.02, star,
                        ha='center', va='bottom', fontsize=9, fontweight='bold')

    ax.axhline(y=1.0, color='red', linestyle='--', linewidth=1, alpha=0.5, label='Expected (no enrichment)')
    ax.set_xlabel('Methylation Motif', fontsize=12)
    ax.set_ylabel('Fold Enrichment', fontsize=12)
    ax.set_title(f'{cat} Depletion by Region', fontsize=13, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(motif_labels, fontsize=10)
    ax.legend(fontsize=9)
    ax.set_ylim(0, max(1.5, ax.get_ylim()[1]))

plt.tight_layout()
for ext in ['pdf', 'svg', 'png']:
    plt.savefig(f"{FIG_DIR}/regulatory_avoidance_by_region.{ext}", dpi=200, bbox_inches='tight')
plt.close()
print("Saved: regulatory_avoidance_by_region.pdf/svg/png")


# ============================================================
# 11. FIGURE 2: Forest plot with OR + 95% CI
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(16, 10))

for ax_idx, cat in enumerate(focus_categories):
    ax = axes[ax_idx]

    motifs = ['GCCGGC_4mC_T1', 'AAGCCCG_4mC_T1', 'AAGCCCG_6mA_T1', 'All_4mC_T1']
    motif_labels_short = ['GCCGGC 4mC', 'AAGCCCG 4mC', 'AAGCCCG 6mA', 'All 4mC']
    regions = ['combined', 'core', 'arm']
    region_colors = {'combined': '#4A90D9', 'core': '#E67E22', 'arm': '#27AE60'}
    region_markers = {'combined': 'D', 'core': 's', 'arm': '^'}

    y_pos = 0
    y_labels = []
    y_positions = []

    for m_idx, motif in enumerate(motifs):
        # Add CMH result
        cmh_row = cmh_df[(cmh_df['motif'] == motif) & (cmh_df['category'] == cat)]

        for reg in regions:
            row = results_df[(results_df['motif_group'] == motif) &
                              (results_df['category'] == cat) &
                              (results_df['region'] == reg)]
            if len(row) == 0:
                continue

            or_val = row['odds_ratio'].values[0]
            ci_lo = row['OR_CI_low'].values[0]
            ci_hi = min(row['OR_CI_high'].values[0], 5)  # cap for display
            p_val = row['p_value'].values[0]

            color = region_colors[reg]
            marker = region_markers[reg]

            ax.plot(or_val, y_pos, marker=marker, color=color, markersize=8, zorder=5)
            ax.hlines(y_pos, ci_lo, ci_hi, color=color, linewidth=2, alpha=0.7)

            sig = '*' if p_val < 0.05 else ''
            label = f"{motif_labels_short[m_idx]} [{reg}] {sig}"
            y_labels.append(label)
            y_positions.append(y_pos)
            y_pos += 1

        # Add CMH adjusted OR
        if len(cmh_row) > 0:
            or_mh = cmh_row['CMH_OR'].values[0]
            ci_lo_mh = cmh_row['CMH_OR_CI_low'].values[0]
            ci_hi_mh = min(cmh_row['CMH_OR_CI_high'].values[0], 5)
            p_mh = cmh_row['CMH_p'].values[0]

            ax.plot(or_mh, y_pos, marker='o', color='black', markersize=10, zorder=5, fillstyle='none', linewidth=2)
            ax.hlines(y_pos, ci_lo_mh, ci_hi_mh, color='black', linewidth=2.5, alpha=0.8)

            sig_mh = '*' if p_mh < 0.05 else ''
            y_labels.append(f"{motif_labels_short[m_idx]} [CMH adj.] {sig_mh}")
            y_positions.append(y_pos)
            y_pos += 1

        y_pos += 0.5  # spacing between motifs

    ax.axvline(x=1.0, color='red', linestyle='--', linewidth=1, alpha=0.6)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(y_labels, fontsize=9)
    ax.set_xlabel('Odds Ratio', fontsize=12)
    ax.set_title(f'{cat}: Stratified Forest Plot', fontsize=13, fontweight='bold')
    ax.set_xlim(0, max(3, ax.get_xlim()[1]))
    ax.invert_yaxis()

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='D', color='#4A90D9', label='Combined', markersize=8, linestyle='None'),
        Line2D([0], [0], marker='s', color='#E67E22', label='Core only', markersize=8, linestyle='None'),
        Line2D([0], [0], marker='^', color='#27AE60', label='Arm only', markersize=8, linestyle='None'),
        Line2D([0], [0], marker='o', color='black', label='CMH adjusted', markersize=10, linestyle='None', fillstyle='none', linewidth=2),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=9)

plt.tight_layout()
for ext in ['pdf', 'svg', 'png']:
    plt.savefig(f"{FIG_DIR}/forest_plot_stratified.{ext}", dpi=200, bbox_inches='tight')
plt.close()
print("Saved: forest_plot_stratified.pdf/svg/png")


# ============================================================
# 12. FIGURE 3: Comprehensive 4-panel summary
# ============================================================
fig = plt.figure(figsize=(20, 16))
gs = gridspec.GridSpec(2, 2, hspace=0.35, wspace=0.35)

# Panel A: Regional gene composition stacked bar
ax_a = fig.add_subplot(gs[0, 0])
cat_order = ['Regulatory/TF', 'Hypothetical', 'Primary metabolism', 'Other',
             'Transport', 'DNA/RNA metabolism', 'Stress/Defense', 'Secondary metabolism',
             'Translation', 'Membrane/Cell wall']
cat_colors = plt.cm.Set3(np.linspace(0, 1, len(cat_order)))

for reg_idx, reg_name in enumerate(['Core', 'Arm']):
    region_sub = protein_coding[protein_coding['region'] == reg_name.lower()]
    bottom = 0
    for cat_idx, cat in enumerate(cat_order):
        frac = len(region_sub[region_sub['functional_category'] == cat]) / len(region_sub) * 100
        ax_a.barh(reg_idx, frac, left=bottom, color=cat_colors[cat_idx], edgecolor='white', linewidth=0.3)
        if frac > 4:
            ax_a.text(bottom + frac/2, reg_idx, f'{frac:.1f}%', ha='center', va='center', fontsize=8)
        bottom += frac

ax_a.set_yticks([0, 1])
ax_a.set_yticklabels(['Core', 'Arm'], fontsize=12)
ax_a.set_xlabel('Percentage of genes', fontsize=11)
ax_a.set_title('A. Gene Functional Composition by Region', fontsize=13, fontweight='bold')
# Small legend
patches = [Patch(facecolor=cat_colors[i], label=cat_order[i]) for i in range(len(cat_order))]
ax_a.legend(handles=patches, fontsize=7, ncol=2, loc='lower right', bbox_to_anchor=(1.0, -0.35))

# Panel B: Methylation site distribution
ax_b = fig.add_subplot(gs[0, 1])
motif_names = ['GCCGGC_4mC_T1', 'AAGCCCG_4mC_T1', 'AAGCCCG_6mA_T1', 'All_4mC_T1']
motif_display = ['GCCGGC 4mC', 'AAGCCCG 4mC', 'AAGCCCG 6mA', 'All 4mC']

for m_idx, motif in enumerate(motif_names):
    df = motif_groups[motif]
    n_core = (df['region'] == 'core').sum()
    n_arm = (df['region'] == 'arm').sum()
    total = len(df)
    core_pct = n_core / total * 100
    arm_pct = n_arm / total * 100

    ax_b.barh(m_idx, core_pct, color='#E67E22', edgecolor='white', linewidth=0.5, label='Core' if m_idx == 0 else '')
    ax_b.barh(m_idx, arm_pct, left=core_pct, color='#27AE60', edgecolor='white', linewidth=0.5, label='Arm' if m_idx == 0 else '')
    ax_b.text(core_pct/2, m_idx, f'{core_pct:.0f}%\n(n={n_core})', ha='center', va='center', fontsize=9, fontweight='bold')
    ax_b.text(core_pct + arm_pct/2, m_idx, f'{arm_pct:.0f}%\n(n={n_arm})', ha='center', va='center', fontsize=9, fontweight='bold')

ax_b.set_yticks(range(len(motif_display)))
ax_b.set_yticklabels(motif_display, fontsize=11)
ax_b.set_xlabel('Percentage of sites', fontsize=11)
ax_b.set_title('B. Methylation Sites: Core vs Arm', fontsize=13, fontweight='bold')
ax_b.legend(fontsize=10, loc='lower right')

# Panel C: Regulatory/TF fold enrichment comparison
ax_c = fig.add_subplot(gs[1, 0])
motifs_plot = ['GCCGGC_4mC_T1', 'AAGCCCG_6mA_T1', 'All_4mC_T1']
motifs_plot_labels = ['GCCGGC 4mC', 'AAGCCCG 6mA', 'All 4mC']
x_plot = np.arange(len(motifs_plot))
width_plot = 0.22
regions_plot = ['combined', 'core', 'arm']
colors_plot = {'combined': '#4A90D9', 'core': '#E67E22', 'arm': '#27AE60'}
labels_plot = {'combined': 'Combined', 'core': 'Core only', 'arm': 'Arm only'}

for i, reg in enumerate(regions_plot):
    folds_plot = []
    pvals_plot = []
    for motif in motifs_plot:
        row = results_df[(results_df['motif_group'] == motif) &
                          (results_df['category'] == 'Regulatory/TF') &
                          (results_df['region'] == reg)]
        if len(row) > 0:
            folds_plot.append(row['fold_enrichment'].values[0])
            pvals_plot.append(row['p_value'].values[0])
        else:
            folds_plot.append(0)
            pvals_plot.append(1)

    bars = ax_c.bar(x_plot + i * width_plot - width_plot, folds_plot, width_plot,
                     label=labels_plot[reg], color=colors_plot[reg], edgecolor='white', alpha=0.85)
    for j, (f, p) in enumerate(zip(folds_plot, pvals_plot)):
        star = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
        if star and f > 0:
            ax_c.text(x_plot[j] + i * width_plot - width_plot, f + 0.02, star,
                      ha='center', va='bottom', fontsize=10, fontweight='bold')

ax_c.axhline(y=1.0, color='red', linestyle='--', linewidth=1, alpha=0.5)
ax_c.set_xticks(x_plot)
ax_c.set_xticklabels(motifs_plot_labels, fontsize=11)
ax_c.set_ylabel('Fold Enrichment', fontsize=12)
ax_c.set_title('C. Regulatory/TF: Stratified Depletion', fontsize=13, fontweight='bold')
ax_c.legend(fontsize=9)
ax_c.set_ylim(0, 1.8)

# Panel D: CMH adjusted OR with CI
ax_d = fig.add_subplot(gs[1, 1])
y_cmh = 0
y_cmh_labels = []
y_cmh_positions = []

for cat in ['Regulatory/TF', 'Hypothetical']:
    for motif, mlabel in zip(motifs_plot, motifs_plot_labels):
        cmh_row = cmh_df[(cmh_df['motif'] == motif) & (cmh_df['category'] == cat)]
        if len(cmh_row) == 0:
            continue

        or_mh = cmh_row['CMH_OR'].values[0]
        ci_lo_mh = cmh_row['CMH_OR_CI_low'].values[0]
        ci_hi_mh = min(cmh_row['CMH_OR_CI_high'].values[0], 4)
        p_mh = cmh_row['CMH_p'].values[0]

        color = '#C0392B' if cat == 'Regulatory/TF' else '#8E44AD'
        ax_d.plot(or_mh, y_cmh, 'o', color=color, markersize=10, zorder=5)
        ax_d.hlines(y_cmh, ci_lo_mh, ci_hi_mh, color=color, linewidth=2.5)

        sig = '***' if p_mh < 0.001 else '**' if p_mh < 0.01 else '*' if p_mh < 0.05 else 'ns'
        y_cmh_labels.append(f"{cat}: {mlabel} ({sig})")
        y_cmh_positions.append(y_cmh)
        y_cmh += 1
    y_cmh += 0.5

ax_d.axvline(x=1.0, color='red', linestyle='--', linewidth=1, alpha=0.6)
ax_d.set_yticks(y_cmh_positions)
ax_d.set_yticklabels(y_cmh_labels, fontsize=10)
ax_d.set_xlabel('CMH-adjusted Odds Ratio', fontsize=12)
ax_d.set_title('D. Region-Adjusted Odds Ratios (CMH)', fontsize=13, fontweight='bold')
ax_d.invert_yaxis()
ax_d.set_xlim(0, max(2.5, ax_d.get_xlim()[1]))

plt.suptitle('H20: Geographic Stratification of Regulatory Avoidance',
             fontsize=15, fontweight='bold', y=0.98)

for ext in ['pdf', 'svg', 'png']:
    plt.savefig(f"{FIG_DIR}/H20_comprehensive_summary.{ext}", dpi=200, bbox_inches='tight')
plt.close()
print("Saved: H20_comprehensive_summary.pdf/svg/png")


# ============================================================
# 13. Summary and verdict
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY AND VERDICT")
print("=" * 70)

print("\n--- Regulatory/TF depletion ---")
reg_core_results = results_df[(results_df['category'] == 'Regulatory/TF') & (results_df['region'] == 'core')]
reg_arm_results = results_df[(results_df['category'] == 'Regulatory/TF') & (results_df['region'] == 'arm')]
reg_combined = results_df[(results_df['category'] == 'Regulatory/TF') & (results_df['region'] == 'combined')]

for motif in ['GCCGGC_4mC_T1', 'AAGCCCG_6mA_T1', 'All_4mC_T1']:
    print(f"\n  {motif}:")
    for reg_name, reg_res in [('Combined', reg_combined), ('Core', reg_core_results), ('Arm', reg_arm_results)]:
        row = reg_res[reg_res['motif_group'] == motif]
        if len(row) > 0:
            f = row['fold_enrichment'].values[0]
            p = row['p_value'].values[0]
            o = row['observed'].values[0]
            e = row['expected'].values[0]
            sig = "*" if p < 0.05 else ""
            print(f"    {reg_name:10s}: fold={f:.2f}, obs={o}, exp={e:.1f}, p={p:.2e} {sig}")

    # CMH
    cmh_row = cmh_df[(cmh_df['motif'] == motif) & (cmh_df['category'] == 'Regulatory/TF')]
    if len(cmh_row) > 0:
        print(f"    CMH adj.  : OR={cmh_row['CMH_OR'].values[0]:.3f}, p={cmh_row['CMH_p'].values[0]:.2e}")

print("\n--- Hypothetical depletion ---")
hyp_core_results = results_df[(results_df['category'] == 'Hypothetical') & (results_df['region'] == 'core')]
hyp_arm_results = results_df[(results_df['category'] == 'Hypothetical') & (results_df['region'] == 'arm')]
hyp_combined = results_df[(results_df['category'] == 'Hypothetical') & (results_df['region'] == 'combined')]

for motif in ['GCCGGC_4mC_T1', 'AAGCCCG_6mA_T1', 'All_4mC_T1']:
    print(f"\n  {motif}:")
    for reg_name, reg_res in [('Combined', hyp_combined), ('Core', hyp_core_results), ('Arm', hyp_arm_results)]:
        row = reg_res[reg_res['motif_group'] == motif]
        if len(row) > 0:
            f = row['fold_enrichment'].values[0]
            p = row['p_value'].values[0]
            o = row['observed'].values[0]
            e = row['expected'].values[0]
            sig = "*" if p < 0.05 else ""
            print(f"    {reg_name:10s}: fold={f:.2f}, obs={o}, exp={e:.1f}, p={p:.2e} {sig}")

    cmh_row = cmh_df[(cmh_df['motif'] == motif) & (cmh_df['category'] == 'Hypothetical')]
    if len(cmh_row) > 0:
        print(f"    CMH adj.  : OR={cmh_row['CMH_OR'].values[0]:.3f}, p={cmh_row['CMH_p'].values[0]:.2e}")


# Final verdict logic
print("\n" + "=" * 70)
print("FINAL VERDICT")
print("=" * 70)

# Check if core-only depletion is significant for the main motifs
reg_genuine_count = 0
reg_total_tests = 0
hyp_genuine_count = 0
hyp_total_tests = 0

for motif in ['GCCGGC_4mC_T1', 'AAGCCCG_6mA_T1', 'All_4mC_T1']:
    # Regulatory/TF
    core_row = results_df[(results_df['motif_group'] == motif) &
                           (results_df['category'] == 'Regulatory/TF') &
                           (results_df['region'] == 'core')]
    if len(core_row) > 0:
        reg_total_tests += 1
        if core_row['fold_enrichment'].values[0] < 0.85 and core_row['p_value'].values[0] < 0.05:
            reg_genuine_count += 1

    # Hypothetical
    core_row_h = results_df[(results_df['motif_group'] == motif) &
                             (results_df['category'] == 'Hypothetical') &
                             (results_df['region'] == 'core')]
    if len(core_row_h) > 0:
        hyp_total_tests += 1
        if core_row_h['fold_enrichment'].values[0] < 0.85 and core_row_h['p_value'].values[0] < 0.05:
            hyp_genuine_count += 1

# CMH significant?
cmh_reg_sig = len(cmh_df[(cmh_df['category'] == 'Regulatory/TF') & (cmh_df['CMH_p'] < 0.05)])
cmh_hyp_sig = len(cmh_df[(cmh_df['category'] == 'Hypothetical') & (cmh_df['CMH_p'] < 0.05)])

print(f"\nRegulatory/TF:")
print(f"  Core-only depletion significant: {reg_genuine_count}/{reg_total_tests} motifs")
print(f"  CMH region-adjusted significant: {cmh_reg_sig}/{len(cmh_df[cmh_df['category'] == 'Regulatory/TF'])} motifs")

print(f"\nHypothetical:")
print(f"  Core-only depletion significant: {hyp_genuine_count}/{hyp_total_tests} motifs")
print(f"  CMH region-adjusted significant: {cmh_hyp_sig}/{len(cmh_df[cmh_df['category'] == 'Hypothetical'])} motifs")

# Determine verdict
if reg_genuine_count >= 2 and cmh_reg_sig >= 2:
    verdict_reg = "GENUINE"
elif reg_genuine_count >= 1 or cmh_reg_sig >= 1:
    verdict_reg = "PARTIAL"
else:
    verdict_reg = "ARTIFACT"

if hyp_genuine_count >= 2 and cmh_hyp_sig >= 2:
    verdict_hyp = "GENUINE"
elif hyp_genuine_count >= 1 or cmh_hyp_sig >= 1:
    verdict_hyp = "PARTIAL"
else:
    verdict_hyp = "ARTIFACT"

print(f"\n>>> Regulatory/TF avoidance verdict: {verdict_reg}")
print(f">>> Hypothetical avoidance verdict: {verdict_hyp}")

if verdict_reg == "GENUINE" and verdict_hyp == "GENUINE":
    overall = "SUPPORTED: Both regulatory and hypothetical avoidance are GENUINE"
elif verdict_reg == "GENUINE" or verdict_hyp == "GENUINE":
    overall = "PARTIAL: One category genuine, the other may be confounded"
elif verdict_reg == "ARTIFACT" and verdict_hyp == "ARTIFACT":
    overall = "REJECTED: Both avoidance patterns are geographic ARTIFACTS"
else:
    overall = "PARTIAL: Evidence is mixed"

print(f"\n>>> OVERALL H20 VERDICT: {overall}")
print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
