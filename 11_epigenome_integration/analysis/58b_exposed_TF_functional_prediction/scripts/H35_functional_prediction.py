#!/usr/bin/env python3
"""
H35: Functional Prediction of 57 Exposed Transcription Factors
Based on TF family membership, genomic context, and expression patterns.

Predicts the downstream regulatory roles of the 57 exposed TFs that form
a "distributed methylation-responsive regulatory layer" with two antagonistic
programs (activation bloc: modules 1-3, repression bloc: module 4).
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from scipy import stats
from collections import Counter, defaultdict
import warnings
import os
import re

warnings.filterwarnings('ignore')

# ============================================================
# Configuration
# ============================================================
BASE = '/Users/okaban/bioinfo/rna-seq'
ANALYSIS_DIR = f'{BASE}/11_epigenome_integration/analysis/58_exposed_TF_functional_prediction'
FIG_DIR = f'{ANALYSIS_DIR}/figures'
TAB_DIR = f'{ANALYSIS_DIR}/tables'

# Input files
EXPOSED_FILE = f'{BASE}/11_epigenome_integration/analysis/51_exposed_regulators_characteristics/tables/exposed_regulators_full_table.tsv'
ALL_REG_FILE = f'{BASE}/11_epigenome_integration/analysis/52_shielded_exposed_boundary/tables/all_genes_features.tsv'
ANNOT_FILE = f'{BASE}/05_annotation/analysis/05_annotation_260128_v1/tables/gene_annotation_basic.tsv'
MODULE_FILE = f'{BASE}/11_epigenome_integration/analysis/55_exposed_regulatory_module/tables/coexpression_modules.tsv'
COORD_FILE = f'{BASE}/11_epigenome_integration/analysis/29_genomewide_TF_screen/tables/coordinated_regulatory_genes.tsv'
DESEQ_FILE = f'{BASE}/04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_3_vs_1.tsv'
COUNTS_FILE = f'{BASE}/04_deseq2/analysis/04_deseq2_260128_v1/results/normalized_counts_M145.tsv'

# Plot settings
plt.rcParams.update({
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 8,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'font.family': 'sans-serif',
})

# Color palettes
BLOC_COLORS = {'Activation': '#E74C3C', 'Repression': '#3498DB'}
FUNC_CAT_COLORS = {
    'Stress/development': '#E74C3C',
    'Metabolic': '#2ECC71',
    'Signal transduction': '#9B59B6',
    'Defense/resistance': '#F39C12',
    'General': '#95A5A6'
}

# ============================================================
# Step 1: Load data
# ============================================================
print("=" * 70)
print("STEP 1: Loading data")
print("=" * 70)

exposed = pd.read_csv(EXPOSED_FILE, sep='\t')
all_reg = pd.read_csv(ALL_REG_FILE, sep='\t')
annot = pd.read_csv(ANNOT_FILE, sep='\t')
modules = pd.read_csv(MODULE_FILE, sep='\t')
coord = pd.read_csv(COORD_FILE, sep='\t')
deseq = pd.read_csv(DESEQ_FILE, sep='\t')
counts = pd.read_csv(COUNTS_FILE, sep='\t')

# Fill NaN old_locus_tag with locus_tag
exposed['old_locus_tag'] = exposed['old_locus_tag'].fillna(exposed['locus_tag'])

print(f"Exposed regulators: {len(exposed)} genes")
print(f"All regulatory genes: {len(all_reg)} genes")
print(f"Gene annotations: {len(annot)} genes")
print(f"Co-expression modules: {len(modules)} modules")
print(f"DESeq2 results: {len(deseq)} genes")

# ============================================================
# Step 2: TF Family annotation of exposed TFs
# ============================================================
print("\n" + "=" * 70)
print("STEP 2: TF Family annotation")
print("=" * 70)

# Known TF family functions
TF_FAMILY_INFO = {
    'TetR': {
        'full_name': 'TetR/AcrR family',
        'general_functions': 'Antibiotic resistance, efflux pump regulation, secondary metabolite regulation',
        'streptomyces_roles': 'BGC expression control (JadR2/jadomycin, TylR/tylosin), GBL signaling, antibiotic efflux',
        'typical_regulon_size': '1-5 genes',
        'target_categories': 'Efflux pumps, BGC genes, resistance genes',
        'functional_category': 'Defense/resistance'
    },
    'Sigma factor': {
        'full_name': 'Sigma-70 family RNA polymerase sigma factor',
        'general_functions': 'Stress response, developmental transitions, alternative sigma factor cascades',
        'streptomyces_roles': 'SigB-like stress, SigE cell envelope, SigR oxidative stress, WhiG sporulation, BldN aerial mycelium',
        'typical_regulon_size': '10-100+ genes (master regulators)',
        'target_categories': 'Stress genes, developmental genes, morphological differentiation',
        'functional_category': 'Stress/development'
    },
    'Sensor kinase': {
        'full_name': 'Two-component system sensor histidine kinase',
        'general_functions': 'Environmental sensing, signal transduction, adaptive responses',
        'streptomyces_roles': 'AbsA (antibiotic production), CutRS (copper sensing), OsaAB (osmotic stress), DraR-DraK (development)',
        'typical_regulon_size': '5-50 genes (via cognate RR)',
        'target_categories': 'Stress response, antibiotic production, morphogenesis',
        'functional_category': 'Signal transduction'
    },
    'Response regulator': {
        'full_name': 'Two-component system response regulator',
        'general_functions': 'Transcriptional control downstream of sensor kinases',
        'streptomyces_roles': 'RamR (aerial mycelium), AbsA2 (antibiotic production), DevR/ChiR (morphogenesis)',
        'typical_regulon_size': '5-50 genes',
        'target_categories': 'Downstream effectors of TCS signaling',
        'functional_category': 'Signal transduction'
    },
    'MarR': {
        'full_name': 'MarR family winged helix-turn-helix',
        'general_functions': 'Oxidative stress, antibiotic resistance, virulence',
        'streptomyces_roles': 'Oxidative stress response, phenylacetic acid catabolism, antibiotic resistance',
        'typical_regulon_size': '1-10 genes',
        'target_categories': 'Oxidative stress genes, efflux, catabolism',
        'functional_category': 'Stress/development'
    },
    'GntR': {
        'full_name': 'GntR family transcriptional regulator',
        'general_functions': 'Carbon/nitrogen metabolism, amino acid biosynthesis',
        'streptomyces_roles': 'Primary metabolite regulation, gluconate pathway, fatty acid metabolism',
        'typical_regulon_size': '1-5 genes',
        'target_categories': 'Metabolic enzymes, transport',
        'functional_category': 'Metabolic'
    },
    'LysR': {
        'full_name': 'LysR family transcriptional regulator',
        'general_functions': 'Amino acid metabolism, aromatic compound degradation, virulence',
        'streptomyces_roles': 'Lysine biosynthesis, aromatic catabolism, oxidative stress response',
        'typical_regulon_size': '1-15 genes',
        'target_categories': 'Amino acid biosynthesis, catabolism, stress',
        'functional_category': 'Metabolic'
    },
    'ArsR': {
        'full_name': 'ArsR/SmtB family transcription factor',
        'general_functions': 'Metal homeostasis, stress sensing',
        'streptomyces_roles': 'Metal ion efflux, zinc/arsenic/nickel homeostasis',
        'typical_regulon_size': '1-5 genes',
        'target_categories': 'Metal transporters, detoxification',
        'functional_category': 'Defense/resistance'
    },
    'MerR': {
        'full_name': 'MerR family transcriptional regulator',
        'general_functions': 'Metal detoxification, oxidative stress',
        'streptomyces_roles': 'Mercury/copper resistance, oxidative stress, potentially BGC regulation',
        'typical_regulon_size': '1-3 genes',
        'target_categories': 'Metal efflux, redox enzymes',
        'functional_category': 'Defense/resistance'
    },
    'LacI': {
        'full_name': 'LacI family DNA-binding transcriptional regulator',
        'general_functions': 'Sugar metabolism, carbon catabolite repression',
        'streptomyces_roles': 'Sugar transport and catabolism, carbon source utilization',
        'typical_regulon_size': '1-10 genes',
        'target_categories': 'Sugar transporters, catabolic enzymes',
        'functional_category': 'Metabolic'
    },
    'HTH (other)': {
        'full_name': 'Helix-turn-helix domain-containing protein',
        'general_functions': 'Diverse regulatory roles (unclassified HTH)',
        'streptomyces_roles': 'Varied - includes uncharacterized DNA-binding regulators',
        'typical_regulon_size': 'Variable',
        'target_categories': 'Diverse',
        'functional_category': 'General'
    },
    'IclR': {
        'full_name': 'IclR family transcriptional regulator',
        'general_functions': 'Glyoxylate cycle, carbon metabolism',
        'streptomyces_roles': 'Carbon flux control, acetate utilization, aromatic catabolism',
        'typical_regulon_size': '1-5 genes',
        'target_categories': 'Central carbon metabolism, TCA/glyoxylate enzymes',
        'functional_category': 'Metabolic'
    },
    'DeoR': {
        'full_name': 'DeoR/GlpR family DNA-binding transcription regulator',
        'general_functions': 'Nucleotide metabolism, sugar phosphate regulation',
        'streptomyces_roles': 'Deoxyribose/sugar phosphate regulation, glycerol metabolism',
        'typical_regulon_size': '1-5 genes',
        'target_categories': 'Nucleotide catabolism, sugar phosphate transport',
        'functional_category': 'Metabolic'
    },
    'ROK': {
        'full_name': 'ROK family protein',
        'general_functions': 'Sugar metabolism, transcriptional repression',
        'streptomyces_roles': 'Sugar kinases, N-acetylglucosamine metabolism (NagC-like)',
        'typical_regulon_size': '1-3 genes',
        'target_categories': 'Sugar kinases, ABC transport',
        'functional_category': 'Metabolic'
    },
    'WhiB': {
        'full_name': 'WhiB family transcriptional regulator',
        'general_functions': 'Developmental regulation (sporulation in Streptomyces)',
        'streptomyces_roles': 'WhiB1-7: sporulation, cell division, oxidative stress sensing via [4Fe-4S] cluster',
        'typical_regulon_size': '10-50+ genes (master regulators)',
        'target_categories': 'Sporulation, cell division, secondary metabolism',
        'functional_category': 'Stress/development'
    },
    'Other regulatory': {
        'full_name': 'Other regulatory proteins',
        'general_functions': 'Diverse regulatory mechanisms (non-TF regulators)',
        'streptomyces_roles': 'DNA-binding (HU, SSB), transcription-repair coupling (Mfd), P-II nitrogen sensing, NsdB (development), FasR (fatty acid), TraR/DksA (stringent response)',
        'typical_regulon_size': 'Variable (global regulators to specific)',
        'target_categories': 'Diverse: chromosome organization, nitrogen, fatty acid, development',
        'functional_category': 'General'
    }
}

# Normalise column names: LFC_T2vsT1/LFC_T3vsT1 → log2FC_T2/log2FC_T3
for old_c, new_c in [('LFC_T2vsT1', 'log2FC_T2'), ('LFC_T3vsT1', 'log2FC_T3')]:
    if old_c in exposed.columns and new_c not in exposed.columns:
        exposed[new_c] = exposed[old_c]

# Assign module and bloc to each exposed TF
module_assignment = {}
for _, row in modules.iterrows():
    mod_id = row['module_id']
    locus_tags = [lt.strip() for lt in row['locus_tags'].split(',')]
    for lt in locus_tags:
        module_assignment[lt] = int(mod_id)

exposed['module'] = exposed['locus_tag'].map(module_assignment)
# Check for unassigned
unassigned = exposed[exposed['module'].isna()]
if len(unassigned) > 0:
    print(f"WARNING: {len(unassigned)} TFs without module assignment")
    exposed.loc[exposed['module'].isna(), 'module'] = 0

exposed['module'] = exposed['module'].astype(int)
# Assign bloc: activation = module 1 or module 3 (positive LFC direction in n=15 set)
#              repression = module 2 (negative LFC direction in n=15 set)
# Use LFC_T3vsT1 direction as ground truth since module numbering changed
def assign_bloc(row):
    lfc = row.get('log2FC_T3', 0)
    if pd.isna(lfc):
        return 'Unassigned'
    return 'Activation' if lfc > 0 else 'Repression'

exposed['bloc'] = exposed.apply(assign_bloc, axis=1)

# Count TF families
family_counts = exposed['tf_family'].value_counts()
print("\nTF family distribution among 57 exposed TFs:")
for fam, count in family_counts.items():
    print(f"  {fam}: {count}")

# Assign functional category
exposed['functional_category'] = exposed['tf_family'].map(
    lambda x: TF_FAMILY_INFO.get(x, {}).get('functional_category', 'General')
)

# Build the per-TF annotation table
tf_annotation_rows = []
for _, row in exposed.iterrows():
    fam = row['tf_family']
    info = TF_FAMILY_INFO.get(fam, {})
    tf_annotation_rows.append({
        'locus_tag': row['locus_tag'],
        'old_locus_tag': row['old_locus_tag'],
        'gene_name': row['gene_name'],
        'product': row['product'],
        'tf_family': fam,
        'family_full_name': info.get('full_name', fam),
        'general_functions': info.get('general_functions', ''),
        'streptomyces_roles': info.get('streptomyces_roles', ''),
        'typical_regulon_size': info.get('typical_regulon_size', ''),
        'target_categories': info.get('target_categories', ''),
        'functional_category': info.get('functional_category', 'General'),
        'module': row['module'],
        'bloc': row['bloc'],
        'log2FC_T3': row['log2FC_T3'],
        'temporal_pattern': row.get('temporal_pattern', 'N/A')
    })

tf_annot_df = pd.DataFrame(tf_annotation_rows)

print(f"\nFunctional category distribution:")
for cat, count in exposed['functional_category'].value_counts().items():
    print(f"  {cat}: {count}")

# ============================================================
# Step 3: Bloc-level functional composition
# ============================================================
print("\n" + "=" * 70)
print("STEP 3: Bloc-level functional composition")
print("=" * 70)

activation = exposed[exposed['bloc'] == 'Activation']
repression = exposed[exposed['bloc'] == 'Repression']

print(f"\nActivation bloc (modules 1-3): {len(activation)} TFs")
print(f"Repression bloc (module 4): {len(repression)} TFs")

# Family composition by bloc
all_families = sorted(exposed['tf_family'].unique())
bloc_family_counts = pd.DataFrame(index=all_families, columns=['Activation', 'Repression'])
for fam in all_families:
    bloc_family_counts.loc[fam, 'Activation'] = len(activation[activation['tf_family'] == fam])
    bloc_family_counts.loc[fam, 'Repression'] = len(repression[repression['tf_family'] == fam])

bloc_family_counts = bloc_family_counts.astype(int)
print("\nFamily composition by bloc:")
print(bloc_family_counts.to_string())

# Fisher exact tests for family enrichment
fisher_results = []
n_act = len(activation)
n_rep = len(repression)
n_total = n_act + n_rep

for fam in all_families:
    a = bloc_family_counts.loc[fam, 'Activation']  # family in activation
    b = bloc_family_counts.loc[fam, 'Repression']   # family in repression
    c = n_act - a  # non-family in activation
    d = n_rep - b  # non-family in repression
    table = [[a, b], [c, d]]
    odds_ratio, p_val = stats.fisher_exact(table)

    fisher_results.append({
        'tf_family': fam,
        'activation_count': a,
        'repression_count': b,
        'activation_pct': a / n_act * 100 if n_act > 0 else 0,
        'repression_pct': b / n_rep * 100 if n_rep > 0 else 0,
        'odds_ratio': odds_ratio,
        'p_value': p_val,
        'enriched_in': 'Activation' if odds_ratio > 1 else 'Repression' if odds_ratio < 1 else 'Neither',
        'functional_category': TF_FAMILY_INFO.get(fam, {}).get('functional_category', 'General')
    })

fisher_df = pd.DataFrame(fisher_results)
# FDR correction
from statsmodels.stats.multitest import multipletests
_, fisher_df['p_adj'], _, _ = multipletests(fisher_df['p_value'], method='fdr_bh')

print("\nFisher exact test results (family enrichment by bloc):")
for _, row in fisher_df.sort_values('p_value').iterrows():
    sig = '*' if row['p_value'] < 0.05 else ''
    print(f"  {row['tf_family']:20s}: OR={row['odds_ratio']:.2f}, p={row['p_value']:.4f} {sig} [{row['enriched_in']}]")

# Functional category enrichment
cat_enrichment = []
all_cats = sorted(exposed['functional_category'].unique())
for cat in all_cats:
    a = len(activation[activation['functional_category'] == cat])
    b = len(repression[repression['functional_category'] == cat])
    c = n_act - a
    d = n_rep - b
    table = [[a, b], [c, d]]
    odds_ratio, p_val = stats.fisher_exact(table)
    cat_enrichment.append({
        'functional_category': cat,
        'activation_count': a,
        'repression_count': b,
        'activation_pct': a / n_act * 100,
        'repression_pct': b / n_rep * 100,
        'odds_ratio': odds_ratio,
        'p_value': p_val,
        'enriched_in': 'Activation' if odds_ratio > 1 else 'Repression' if odds_ratio < 1 else 'Neither'
    })

cat_enrich_df = pd.DataFrame(cat_enrichment)
_, cat_enrich_df['p_adj'], _, _ = multipletests(cat_enrich_df['p_value'], method='fdr_bh')

print("\nFunctional category enrichment by bloc:")
for _, row in cat_enrich_df.sort_values('p_value').iterrows():
    sig = '*' if row['p_value'] < 0.05 else ''
    print(f"  {row['functional_category']:25s}: Act={row['activation_count']}({row['activation_pct']:.1f}%) Rep={row['repression_count']}({row['repression_pct']:.1f}%) OR={row['odds_ratio']:.2f} p={row['p_value']:.4f} {sig}")

# ============================================================
# Step 4: Genomic context analysis
# ============================================================
print("\n" + "=" * 70)
print("STEP 4: Genomic context analysis (±5kb neighbors)")
print("=" * 70)

# Function category keywords for neighbor analysis
NEIGHBOR_CATEGORIES = {
    'secondary_metabolism': ['polyketide', 'NRPS', 'synthase', 'terpene', 'siderophore',
                            'PKS', 'chalcone', 'aminotransferase', 'thioesterase',
                            'acyl-CoA', 'malonyl', 'acyltransferase', 'cyclase',
                            'condensation', 'oxidoreductase'],
    'transport': ['transporter', 'permease', 'efflux', 'ABC', 'MFS', 'export',
                  'import', 'porin', 'channel', 'pump', 'secretion'],
    'signaling': ['kinase', 'phosphatase', 'response regulator', 'sensor', 'sigma',
                  'cyclase', 'diguanylate', 'phosphodiesterase', 'GGDEF', 'EAL'],
    'primary_metabolism': ['dehydrogenase', 'reductase', 'transferase', 'isomerase',
                          'lyase', 'ligase', 'synthase', 'mutase', 'epimerase',
                          'kinase', 'phosphatase', 'aldolase', 'carboxylase'],
    'stress_response': ['chaperone', 'protease', 'heat shock', 'cold shock', 'anti-sigma',
                        'universal stress', 'catalase', 'superoxide', 'thioredoxin',
                        'glutaredoxin', 'peroxidase'],
    'cell_wall_membrane': ['peptidoglycan', 'murein', 'penicillin-binding', 'lipid',
                           'glycosyltransferase', 'membrane', 'lipoprotein'],
    'DNA_RNA_processing': ['helicase', 'nuclease', 'recombinase', 'topoisomerase',
                           'polymerase', 'primase', 'methyltransferase', 'restriction'],
    'hypothetical': ['hypothetical', 'uncharacterized', 'DUF', 'unknown function']
}

def categorize_product(product):
    """Categorize a gene product annotation."""
    if pd.isna(product):
        return 'hypothetical'
    product_lower = product.lower()

    # Check secondary_metabolism first (more specific)
    for kw in NEIGHBOR_CATEGORIES['secondary_metabolism']:
        if kw.lower() in product_lower:
            return 'secondary_metabolism'

    for kw in NEIGHBOR_CATEGORIES['transport']:
        if kw.lower() in product_lower:
            return 'transport'

    for kw in NEIGHBOR_CATEGORIES['signaling']:
        if kw.lower() in product_lower:
            return 'signaling'

    for kw in NEIGHBOR_CATEGORIES['stress_response']:
        if kw.lower() in product_lower:
            return 'stress_response'

    for kw in NEIGHBOR_CATEGORIES['cell_wall_membrane']:
        if kw.lower() in product_lower:
            return 'cell_wall_membrane'

    for kw in NEIGHBOR_CATEGORIES['DNA_RNA_processing']:
        if kw.lower() in product_lower:
            return 'DNA_RNA_processing'

    for kw in NEIGHBOR_CATEGORIES['hypothetical']:
        if kw.lower() in product_lower:
            return 'hypothetical'

    for kw in NEIGHBOR_CATEGORIES['primary_metabolism']:
        if kw.lower() in product_lower:
            return 'primary_metabolism'

    return 'other'

# For each exposed TF, find genes within ±5kb
neighbor_results = []
for _, tf_row in exposed.iterrows():
    tf_start = tf_row['start']
    tf_end = tf_row['end']
    tf_mid = (tf_start + tf_end) / 2
    window = 5000

    neighbors = annot[
        (annot['start'] >= tf_start - window) &
        (annot['end'] <= tf_end + window) &
        (annot['gene_id'] != tf_row['locus_tag'])
    ]

    for _, nb in neighbors.iterrows():
        cat = categorize_product(nb['product'])
        nb_mid = (nb['start'] + nb['end']) / 2
        distance = abs(nb_mid - tf_mid)

        neighbor_results.append({
            'tf_locus_tag': tf_row['locus_tag'],
            'tf_old_locus_tag': tf_row['old_locus_tag'],
            'tf_family': tf_row['tf_family'],
            'bloc': tf_row['bloc'],
            'neighbor_locus_tag': nb['gene_id'],
            'neighbor_old_locus_tag': nb.get('old_locus_tag', ''),
            'neighbor_product': nb['product'],
            'neighbor_category': cat,
            'distance_bp': int(distance),
            'neighbor_start': nb['start'],
            'neighbor_end': nb['end'],
            'neighbor_strand': nb['strand']
        })

neighbor_df = pd.DataFrame(neighbor_results)
print(f"Total neighbor genes found: {len(neighbor_df)}")
print(f"Unique TFs with neighbors: {neighbor_df['tf_locus_tag'].nunique()}")

# Summarize neighbor categories by bloc
neighbor_summary_act = neighbor_df[neighbor_df['bloc'] == 'Activation']['neighbor_category'].value_counts()
neighbor_summary_rep = neighbor_df[neighbor_df['bloc'] == 'Repression']['neighbor_category'].value_counts()

all_neighbor_cats = sorted(set(list(neighbor_summary_act.index) + list(neighbor_summary_rep.index)))
neighbor_comparison = pd.DataFrame(index=all_neighbor_cats)
neighbor_comparison['activation_count'] = [neighbor_summary_act.get(c, 0) for c in all_neighbor_cats]
neighbor_comparison['repression_count'] = [neighbor_summary_rep.get(c, 0) for c in all_neighbor_cats]
neighbor_comparison['activation_pct'] = neighbor_comparison['activation_count'] / neighbor_comparison['activation_count'].sum() * 100
neighbor_comparison['repression_pct'] = neighbor_comparison['repression_count'] / neighbor_comparison['repression_count'].sum() * 100

print("\nNeighbor function profiles by bloc:")
print(neighbor_comparison.to_string())

# Fisher tests for neighbor categories
neighbor_fisher = []
n_act_nb = neighbor_comparison['activation_count'].sum()
n_rep_nb = neighbor_comparison['repression_count'].sum()

for cat in all_neighbor_cats:
    a = int(neighbor_comparison.loc[cat, 'activation_count'])
    b = int(neighbor_comparison.loc[cat, 'repression_count'])
    c = int(n_act_nb - a)
    d = int(n_rep_nb - b)
    table = [[a, b], [c, d]]
    odds_ratio, p_val = stats.fisher_exact(table)
    neighbor_fisher.append({
        'category': cat,
        'activation_count': a,
        'repression_count': b,
        'odds_ratio': odds_ratio,
        'p_value': p_val
    })

neighbor_fisher_df = pd.DataFrame(neighbor_fisher)
if len(neighbor_fisher_df) > 0:
    _, neighbor_fisher_df['p_adj'], _, _ = multipletests(neighbor_fisher_df['p_value'], method='fdr_bh')

print("\nNeighbor category enrichment tests:")
for _, row in neighbor_fisher_df.sort_values('p_value').iterrows():
    sig = '*' if row['p_value'] < 0.05 else ''
    print(f"  {row['category']:25s}: OR={row['odds_ratio']:.2f}, p={row['p_value']:.4f} {sig}")

# Identify TFs proximal to known gene clusters
proximal_clusters = []
for _, tf_row in exposed.iterrows():
    tf_neighbors = neighbor_df[neighbor_df['tf_locus_tag'] == tf_row['locus_tag']]
    sec_met_nb = tf_neighbors[tf_neighbors['neighbor_category'] == 'secondary_metabolism']
    transport_nb = tf_neighbors[tf_neighbors['neighbor_category'] == 'transport']

    if len(sec_met_nb) > 0:
        proximal_clusters.append({
            'locus_tag': tf_row['locus_tag'],
            'old_locus_tag': tf_row['old_locus_tag'],
            'tf_family': tf_row['tf_family'],
            'bloc': tf_row['bloc'],
            'cluster_type': 'secondary_metabolism',
            'n_cluster_neighbors': len(sec_met_nb),
            'neighbor_products': '; '.join(sec_met_nb['neighbor_product'].values)
        })
    if len(transport_nb) >= 2:
        proximal_clusters.append({
            'locus_tag': tf_row['locus_tag'],
            'old_locus_tag': tf_row['old_locus_tag'],
            'tf_family': tf_row['tf_family'],
            'bloc': tf_row['bloc'],
            'cluster_type': 'transport_cluster',
            'n_cluster_neighbors': len(transport_nb),
            'neighbor_products': '; '.join(transport_nb['neighbor_product'].values)
        })

if proximal_clusters:
    proximal_df = pd.DataFrame(proximal_clusters)
    print(f"\nExposed TFs proximal to functional clusters: {len(proximal_df)}")
    for _, row in proximal_df.iterrows():
        print(f"  {row['old_locus_tag']} ({row['tf_family']}, {row['bloc']}): {row['cluster_type']} ({row['n_cluster_neighbors']} neighbors)")

# ============================================================
# Step 5: Product annotation keyword analysis
# ============================================================
print("\n" + "=" * 70)
print("STEP 5: Product annotation keyword analysis")
print("=" * 70)

product_keywords = {
    'kinase': [], 'regulator': [], 'sigma': [], 'sensor': [], 'response': [],
    'reductase': [], 'synthase': [], 'transcriptional': [], 'DNA-binding': [],
    'helix-turn-helix': [], 'family': [], 'winged': [], 'tetratricopeptide': [],
    'fatty acid': [], 'nitrogen': [], 'uracil': [], 'coupling': [],
    'single-stranded': [], 'NsdB': [], 'TraR': [], 'DksA': [],
    'FxSxx': [], 'AfsR': [], 'WhiB': [], 'RamR': [], 'FasR': [],
    'SigB': [], 'SigE': [], 'SigJ': []
}

for _, row in exposed.iterrows():
    product = str(row['product'])
    for kw in product_keywords:
        if kw.lower() in product.lower():
            product_keywords[kw].append(row['locus_tag'])

# Predict regulatory scope for each TF
predicted_scope = []
for _, row in exposed.iterrows():
    product = str(row['product'])
    fam = row['tf_family']
    bloc = row['bloc']
    old_lt = row['old_locus_tag']

    # Get neighbor context
    tf_neighbors = neighbor_df[neighbor_df['tf_locus_tag'] == row['locus_tag']]
    dominant_nb_cat = tf_neighbors['neighbor_category'].mode().iloc[0] if len(tf_neighbors) > 0 else 'unknown'

    # Predicted scope based on family + product + context
    if fam == 'Sigma factor':
        if 'SigE' in product:
            scope = 'Cell envelope stress response, extracytoplasmic function'
        elif 'SigB' in product or 'SigF' in product or 'SigG' in product:
            scope = 'General stress response, stationary phase/sporulation'
        elif 'SigJ' in product:
            scope = 'Alternative sigma factor, possible stress/development'
        else:
            scope = 'Transcriptional reprogramming, stress/developmental switch'
    elif fam == 'Sensor kinase':
        scope = 'Environmental signal sensing, phosphorelay to cognate RR'
    elif fam == 'Response regulator':
        if 'RamR' in product:
            scope = 'Aerial mycelium formation, morphological differentiation'
        else:
            scope = 'TCS-mediated transcriptional regulation'
    elif fam == 'TetR':
        if 'tetratricopeptide' in product.lower():
            scope = 'Protein-protein interaction hub, possible chaperone/adaptor function'
        else:
            scope = 'Local repression of neighboring efflux/resistance genes'
    elif fam == 'MarR':
        scope = 'Oxidative stress sensing, derepression of stress response genes'
    elif fam == 'GntR':
        scope = 'Carbon/nitrogen metabolite sensing, metabolic gene regulation'
    elif fam == 'LysR':
        scope = 'Amino acid/aromatic compound metabolism regulation'
    elif fam == 'IclR':
        scope = 'Carbon flux control, glyoxylate shunt regulation'
    elif fam == 'LacI':
        scope = 'Sugar transport/catabolism repression'
    elif fam == 'ArsR':
        scope = 'Metal ion sensing, metal efflux gene regulation'
    elif fam == 'MerR':
        scope = 'Metal ion-responsive transcriptional activation'
    elif fam == 'DeoR':
        scope = 'Sugar phosphate/nucleotide metabolism regulation'
    elif fam == 'ROK':
        scope = 'Sugar kinase/regulatory function'
    elif fam == 'WhiB':
        scope = 'Sporulation/developmental master regulation via Fe-S sensing'
    elif fam == 'Other regulatory':
        if 'NsdB' in product:
            scope = 'Development regulation (NsdB pathway)'
        elif 'FasR' in product:
            scope = 'Fatty acid biosynthesis regulation'
        elif 'nitrogen' in product.lower():
            scope = 'Nitrogen assimilation signaling'
        elif 'TraR' in product or 'DksA' in product:
            scope = 'Stringent response, rRNA/ribosome regulation'
        elif 'uracil' in product.lower():
            scope = 'DNA repair, uracil excision'
        elif 'single-stranded' in product.lower():
            scope = 'DNA replication/repair, recombination'
        elif 'coupling' in product.lower() or 'Mfd' in product:
            scope = 'Transcription-coupled DNA repair'
        elif 'HU' in product:
            scope = 'Chromosome organization, nucleoid structuring'
        elif 'LuxR' in product:
            scope = 'Quorum sensing response, possible BGC regulation'
        elif 'AfsR' in product or 'TcrA' in product:
            scope = 'Global antibiotic regulation (SARP-like)'
        else:
            scope = 'Non-canonical regulatory function'
    else:
        scope = 'Unclassified regulatory function'

    predicted_scope.append({
        'locus_tag': row['locus_tag'],
        'old_locus_tag': old_lt,
        'tf_family': fam,
        'product': product,
        'bloc': bloc,
        'module': row['module'],
        'dominant_neighbor_category': dominant_nb_cat,
        'predicted_regulatory_scope': scope
    })

scope_df = pd.DataFrame(predicted_scope)
print("\nPredicted regulatory scopes:")
for _, row in scope_df.iterrows():
    olt = str(row['old_locus_tag']) if pd.notna(row['old_locus_tag']) else row['locus_tag']
    print(f"  {olt:12s} ({row['tf_family']:17s}, {row['bloc']:10s}): {row['predicted_regulatory_scope']}")

# ============================================================
# Step 6: Expression magnitude as functional indicator
# ============================================================
print("\n" + "=" * 70)
print("STEP 6: Expression magnitude analysis")
print("=" * 70)

exposed['abs_LFC_T3'] = exposed['log2FC_T3'].abs()
exposed['abs_LFC_T2'] = exposed['log2FC_T2'].abs()

# Rank by |LFC_T3|
exposed_ranked = exposed.sort_values('abs_LFC_T3', ascending=False)
print("\nTop 15 exposed TFs by |LFC_T3vsT1|:")
for i, (_, row) in enumerate(exposed_ranked.head(15).iterrows()):
    olt = str(row['old_locus_tag']) if pd.notna(row['old_locus_tag']) else row['locus_tag']
    print(f"  {i+1:2d}. {olt:12s} ({row['tf_family']:17s}, {row['bloc']:10s}): |LFC|={row['abs_LFC_T3']:.2f}")

# Are top responders concentrated in specific families?
top_half = exposed_ranked.head(31)
bottom_half = exposed_ranked.tail(31)

print("\nTop half family distribution:")
print(top_half['tf_family'].value_counts().to_string())
print("\nBottom half family distribution:")
print(bottom_half['tf_family'].value_counts().to_string())

# Top vs bottom bloc distribution
print(f"\nTop half bloc: Activation={len(top_half[top_half['bloc']=='Activation'])}, Repression={len(top_half[top_half['bloc']=='Repression'])}")
print(f"Bottom half bloc: Activation={len(bottom_half[bottom_half['bloc']=='Activation'])}, Repression={len(bottom_half[bottom_half['bloc']=='Repression'])}")

# Per-family mean |LFC|
family_lfc = exposed.groupby('tf_family')['abs_LFC_T3'].agg(['mean', 'median', 'std', 'count'])
family_lfc = family_lfc.sort_values('mean', ascending=False)
print("\nPer-family |LFC_T3| statistics:")
print(family_lfc.to_string())

# Kruskal-Wallis test for family differences
families_with_enough = [fam for fam in all_families if len(exposed[exposed['tf_family'] == fam]) >= 2]
groups = [exposed[exposed['tf_family'] == fam]['abs_LFC_T3'].values for fam in families_with_enough]
if len(groups) >= 2:
    kw_stat, kw_p = stats.kruskal(*groups)
    print(f"\nKruskal-Wallis test for family differences in |LFC_T3|: H={kw_stat:.3f}, p={kw_p:.4f}")
else:
    kw_stat, kw_p = np.nan, np.nan

# Bloc difference in |LFC| - re-select from updated exposed
activation = exposed[exposed['bloc'] == 'Activation']
repression = exposed[exposed['bloc'] == 'Repression']
act_lfc = activation['abs_LFC_T3'].values
rep_lfc = repression['abs_LFC_T3'].values
mw_stat, mw_p = stats.mannwhitneyu(act_lfc, rep_lfc, alternative='two-sided')
print(f"Bloc |LFC| comparison: Activation median={np.median(act_lfc):.2f}, Repression median={np.median(rep_lfc):.2f}, MWU p={mw_p:.4f}")

# Expression stats table
expr_stats = []
for fam in all_families:
    fam_data = exposed[exposed['tf_family'] == fam]
    expr_stats.append({
        'tf_family': fam,
        'n': len(fam_data),
        'mean_abs_LFC_T3': fam_data['abs_LFC_T3'].mean(),
        'median_abs_LFC_T3': fam_data['abs_LFC_T3'].median(),
        'sd_abs_LFC_T3': fam_data['abs_LFC_T3'].std() if len(fam_data) > 1 else np.nan,
        'mean_abs_LFC_T2': fam_data['abs_LFC_T2'].mean(),
        'mean_baseMean': fam_data['baseMean'].mean(),
        'n_activation': len(fam_data[fam_data['bloc'] == 'Activation']),
        'n_repression': len(fam_data[fam_data['bloc'] == 'Repression']),
        'functional_category': TF_FAMILY_INFO.get(fam, {}).get('functional_category', 'General')
    })

expr_stats_df = pd.DataFrame(expr_stats).sort_values('mean_abs_LFC_T3', ascending=False)

# ============================================================
# Step 7: Known Streptomyces regulons literature mapping
# ============================================================
print("\n" + "=" * 70)
print("STEP 7: Known Streptomyces regulons literature mapping")
print("=" * 70)

# Map using SCO numbers
sco_to_known = {
    # Well-characterized S. coelicolor regulators
    'SCO3134': {'name': 'SCO3134', 'known_function': 'Uncharacterized response regulator, proximal to SCO3133 (sensor kinase)'},
    'SCO5828': {'name': 'SCO5828', 'known_function': 'Uncharacterized response regulator'},
    'SCO7648': {'name': 'SCO7648', 'known_function': 'Uncharacterized response regulator'},
    'SCO5434': {'name': 'SCO5434', 'known_function': 'Uncharacterized response regulator'},
    'SCO6685': {'name': 'RamR', 'known_function': 'Aerial mycelium formation; RamCSAB pathway; activates morphogenesis'},
    'SCO4938': {'name': 'SigJ (SCO4938)', 'known_function': 'ECF sigma factor; predicted cell wall/envelope stress response'},
    'SCO4005': {'name': 'SigE (SCO4005)', 'known_function': 'Cell envelope stress sigma factor; heat shock and cell wall damage response'},
    'SCO2639': {'name': 'SCO2639', 'known_function': 'Predicted sigma factor; possible stress response'},
    'SCO1564': {'name': 'SCO1564', 'known_function': 'Predicted sigma factor; developmental timing'},
    'SCO0632': {'name': 'SCO0632', 'known_function': 'Predicted sigma factor'},
    'SCO7314': {'name': 'SCO7314', 'known_function': 'SigB/SigF/SigG-like; general stress and/or sporulation sigma factor'},
    'SCO4960': {'name': 'SCO4960', 'known_function': 'Sigma factor-like HTH; predicted developmental regulator'},
    'SCO7252': {'name': 'NsdB', 'known_function': 'Negative regulator of sporulation and antibiotic production; DNA-binding'},
    'SCO2386': {'name': 'FasR', 'known_function': 'Fatty acid biosynthesis transcriptional regulator; controls fab gene cluster'},
    'SCO5433': {'name': 'TcrA', 'known_function': 'AfsR-like SARP regulator; involved in terpene/secondary metabolite regulation'},
    'SCO3109': {'name': 'Mfd', 'known_function': 'Transcription-repair coupling factor; mutation frequency decline'},
    'SCO2681': {'name': 'FxsT (SCO2681)', 'known_function': 'FxSxx-COOH system TPR protein; unknown TetR-family, possible signaling'},
    'SCO1008': {'name': 'SCO1008', 'known_function': 'ArsR family; predicted metal ion sensing'},
    'SCO5289': {'name': 'SCO5289', 'known_function': 'Sensor histidine kinase; uncharacterized TCS'},
    'SCO1160': {'name': 'SCO1160', 'known_function': 'Sensor histidine kinase; massive T2 induction (LFC=6.8)'},
    'SCO5584': {'name': 'GlnK-like (SCO5584)', 'known_function': 'P-II family nitrogen regulator; nitrogen sensing/signaling'},
    'SCO4122': {'name': 'SCO4122', 'known_function': 'MarR family; possibly oxidative stress-responsive'},
    'SCO1191': {'name': 'SCO1191', 'known_function': 'MarR family; predicted stress response'},
}

# TetR family analysis
print("\n--- TetR family (9 exposed TFs) ---")
tetr_tfs = exposed[exposed['tf_family'] == 'TetR']
for _, row in tetr_tfs.iterrows():
    sco = row['old_locus_tag']
    nb = neighbor_df[neighbor_df['tf_locus_tag'] == row['locus_tag']]
    bgc_nb = nb[nb['neighbor_category'] == 'secondary_metabolism']
    efflux_nb = nb[nb['neighbor_category'] == 'transport']
    print(f"  {sco} ({row['bloc']}): LFC_T3={row['log2FC_T3']:.2f}, "
          f"BGC_neighbors={len(bgc_nb)}, Transport_neighbors={len(efflux_nb)}")
    if len(bgc_nb) > 0:
        for _, n in bgc_nb.iterrows():
            print(f"    -> BGC neighbor: {n['neighbor_product']}")
    if len(efflux_nb) > 0:
        for _, n in efflux_nb.iterrows():
            print(f"    -> Transport neighbor: {n['neighbor_product']}")

# Sigma factor analysis
print("\n--- Sigma factors (7 exposed TFs) ---")
sig_tfs = exposed[exposed['tf_family'] == 'Sigma factor']
for _, row in sig_tfs.iterrows():
    sco = row['old_locus_tag']
    product = row['product']
    subfamily = ''
    if 'SigE' in product:
        subfamily = 'SigE (cell envelope stress)'
    elif 'SigB' in product or 'SigF' in product or 'SigG' in product:
        subfamily = 'SigB/F/G (general stress/sporulation)'
    elif 'SigJ' in product:
        subfamily = 'SigJ (ECF sigma factor)'
    else:
        subfamily = 'Unclassified sigma-70'

    print(f"  {sco}: {product} -> {subfamily}, LFC_T3={row['log2FC_T3']:.2f} [{row['bloc']}]")

# TCS analysis
print("\n--- TCS pairs (6 SK + 5 RR = 11 TFs) ---")
tcs_tfs = exposed[exposed['tf_family'].isin(['Sensor kinase', 'Response regulator'])]
sk_tfs = tcs_tfs[tcs_tfs['tf_family'] == 'Sensor kinase']
rr_tfs = tcs_tfs[tcs_tfs['tf_family'] == 'Response regulator']

print(f"  Sensor kinases: {len(sk_tfs)}")
for _, row in sk_tfs.iterrows():
    print(f"    {row['old_locus_tag']} ({row['bloc']}): LFC_T3={row['log2FC_T3']:.2f}")

print(f"  Response regulators: {len(rr_tfs)}")
for _, row in rr_tfs.iterrows():
    print(f"    {row['old_locus_tag']} ({row['bloc']}): LFC_T3={row['log2FC_T3']:.2f}")

# Identify TCS pairs by genomic proximity
print("\n  Candidate TCS pairs (SK + RR within 3kb):")
tcs_pairs = []
for _, sk in sk_tfs.iterrows():
    for _, rr in rr_tfs.iterrows():
        dist = min(abs(sk['start'] - rr['end']), abs(rr['start'] - sk['end']))
        if dist < 3000:
            tcs_pairs.append({
                'SK_locus_tag': sk['locus_tag'],
                'SK_old_locus_tag': sk['old_locus_tag'],
                'RR_locus_tag': rr['locus_tag'],
                'RR_old_locus_tag': rr['old_locus_tag'],
                'distance_bp': dist,
                'SK_bloc': sk['bloc'],
                'RR_bloc': rr['bloc'],
                'SK_LFC_T3': sk['log2FC_T3'],
                'RR_LFC_T3': rr['log2FC_T3']
            })

for pair in tcs_pairs:
    print(f"    {pair['SK_old_locus_tag']} (SK) <-> {pair['RR_old_locus_tag']} (RR), dist={pair['distance_bp']}bp")
    print(f"      SK: {pair['SK_bloc']}, LFC={pair['SK_LFC_T3']:.2f} | RR: {pair['RR_bloc']}, LFC={pair['RR_LFC_T3']:.2f}")

# ============================================================
# Step 8: Predicted regulatory model
# ============================================================
print("\n" + "=" * 70)
print("STEP 8: Predicted regulatory model for the two blocs")
print("=" * 70)

# Synthesize predictions for activation bloc
act_predictions = []
act_families = activation['tf_family'].value_counts()
print("\n=== ACTIVATION BLOC (Modules 1-3, ~35 genes) ===")
print("Predicted activated biological processes during development:")

act_pred_items = [
    {
        'process': 'Morphological differentiation and aerial mycelium formation',
        'evidence': f'RamR (SCO6685, LFC=+7.9), sigma factors (SigE SCO4005 LFC=+3.2, SigB/F/G SCO7314 LFC=+1.5), WhiB (LFC=+1.7)',
        'confidence': 'High',
        'families': 'Sigma factor, Response regulator, WhiB'
    },
    {
        'process': 'Two-component signal transduction cascades',
        'evidence': f'4 sensor kinases + 3 response regulators in activation bloc, massive LFC for SCO1160 (SK, +5.6) and SCO7711 (SK, +3.1)',
        'confidence': 'High',
        'families': 'Sensor kinase, Response regulator'
    },
    {
        'process': 'Secondary metabolite / antibiotic production regulation',
        'evidence': 'TcrA (SCO5433, AfsR-like SARP, LFC=+1.5), NsdB (SCO7252, LFC=+7.8), TetR family members near transport genes',
        'confidence': 'Medium',
        'families': 'Other regulatory, TetR'
    },
    {
        'process': 'Stress response activation (oxidative, envelope, general)',
        'evidence': 'SigE (cell envelope), SigB/F/G (general stress), MerR family members (metal/oxidative)',
        'confidence': 'Medium',
        'families': 'Sigma factor, MerR'
    },
    {
        'process': 'Nucleotide/sugar metabolism reorganization',
        'evidence': 'DeoR (SCO1897, LFC=+1.5), LacI (SCO7411, LFC=+1.1)',
        'confidence': 'Low',
        'families': 'DeoR, LacI'
    }
]

for pred in act_pred_items:
    print(f"\n  [{pred['confidence']}] {pred['process']}")
    print(f"    Evidence: {pred['evidence']}")
    act_predictions.append(pred)

print("\n=== REPRESSION BLOC (Module 4, 26 genes) ===")
print("Predicted repressed biological processes during development:")

rep_pred_items = [
    {
        'process': 'Primary carbon/nitrogen metabolism shutdown',
        'evidence': f'GntR (SCO0823, LFC=-1.5), IclR (SCO4989, LFC=-2.4), LacI (SCO4158, LFC=-1.4), LysR (SCO4766, LFC=-1.2)',
        'confidence': 'High',
        'families': 'GntR, IclR, LacI, LysR'
    },
    {
        'process': 'Antibiotic efflux / TetR-mediated defense downregulation',
        'evidence': f'5 TetR members (SCO2223 LFC=-2.0, SCO4305 LFC=-1.2, SCO4639 LFC=-2.3, SCO5956 LFC=-1.0, SCO6599 LFC=-1.6)',
        'confidence': 'High',
        'families': 'TetR'
    },
    {
        'process': 'Fatty acid biosynthesis repression',
        'evidence': 'FasR (SCO2386, LFC=-2.0), fatty acid biosynthesis master regulator',
        'confidence': 'High',
        'families': 'Other regulatory'
    },
    {
        'process': 'DNA replication/repair machinery shutdown',
        'evidence': 'SSB (SCO3907, LFC=-3.3), HU (SCO2950, LFC=-2.9), Mfd (SCO3109, LFC=-1.1), UdgX (SCO4495, LFC=-1.3)',
        'confidence': 'High',
        'families': 'Other regulatory'
    },
    {
        'process': 'Oxidative stress / MarR-mediated defense repression',
        'evidence': f'MarR (SCO1191 LFC=-2.1, SCO4122 LFC=-3.5)',
        'confidence': 'Medium',
        'families': 'MarR'
    },
    {
        'process': 'Environmental sensing downregulation (select TCS)',
        'evidence': 'SCO5289 (SK, LFC=-2.4), SCO5434 (RR, LFC=-1.1)',
        'confidence': 'Medium',
        'families': 'Sensor kinase, Response regulator'
    },
    {
        'process': 'Sigma factor cascade repression',
        'evidence': 'SCO2639 (sigma, LFC=-1.3), SCO4960 (sigma-like, LFC=-1.5)',
        'confidence': 'Medium',
        'families': 'Sigma factor'
    },
    {
        'process': 'Nitrogen assimilation shutdown',
        'evidence': 'P-II nitrogen regulator (SCO5584, LFC=-1.4)',
        'confidence': 'Medium',
        'families': 'Other regulatory'
    }
]

for pred in rep_pred_items:
    print(f"\n  [{pred['confidence']}] {pred['process']}")
    print(f"    Evidence: {pred['evidence']}")

# Combine predictions
all_predictions = []
for p in act_pred_items:
    p['bloc'] = 'Activation'
    all_predictions.append(p)
for p in rep_pred_items:
    p['bloc'] = 'Repression'
    all_predictions.append(p)

predictions_df = pd.DataFrame(all_predictions)

# ============================================================
# Step 9: FIGURES
# ============================================================
print("\n" + "=" * 70)
print("STEP 9: Generating figures")
print("=" * 70)

# --- Figure 1: TF family bloc composition ---
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Left: Stacked bar chart
families_sorted = family_counts.index.tolist()
act_counts = [bloc_family_counts.loc[f, 'Activation'] for f in families_sorted]
rep_counts = [bloc_family_counts.loc[f, 'Repression'] for f in families_sorted]

x = np.arange(len(families_sorted))
width = 0.35

bars1 = axes[0].bar(x - width/2, act_counts, width, label='Activation bloc', color=BLOC_COLORS['Activation'], alpha=0.85, edgecolor='white')
bars2 = axes[0].bar(x + width/2, rep_counts, width, label='Repression bloc', color=BLOC_COLORS['Repression'], alpha=0.85, edgecolor='white')

axes[0].set_xlabel('TF Family')
axes[0].set_ylabel('Count')
axes[0].set_title('TF Family Composition by Bloc')
axes[0].set_xticks(x)
axes[0].set_xticklabels(families_sorted, rotation=45, ha='right', fontsize=8)
axes[0].legend()
axes[0].set_ylim(0, max(max(act_counts), max(rep_counts)) + 2)

# Add counts on bars
for bar in bars1:
    if bar.get_height() > 0:
        axes[0].text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.1,
                    f'{int(bar.get_height())}', ha='center', va='bottom', fontsize=7)
for bar in bars2:
    if bar.get_height() > 0:
        axes[0].text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.1,
                    f'{int(bar.get_height())}', ha='center', va='bottom', fontsize=7)

# Right: Percentage stacked bar
act_total = sum(act_counts)
rep_total = sum(rep_counts)

# Sort families for visual clarity
fam_func_map = {f: TF_FAMILY_INFO.get(f, {}).get('functional_category', 'General') for f in families_sorted}
sorted_fams = sorted(families_sorted, key=lambda f: (fam_func_map[f], f))

# Create stacked bar by functional category
func_cats_ordered = ['Stress/development', 'Signal transduction', 'Defense/resistance', 'Metabolic', 'General']
bloc_names = ['Activation\n(n=35)', 'Repression\n(n=26)']
bottom_act = 0
bottom_rep = 0

for cat in func_cats_ordered:
    cat_fams = [f for f in sorted_fams if fam_func_map[f] == cat]
    act_sum = sum(bloc_family_counts.loc[f, 'Activation'] for f in cat_fams)
    rep_sum = sum(bloc_family_counts.loc[f, 'Repression'] for f in cat_fams)

    act_pct = act_sum / act_total * 100 if act_total > 0 else 0
    rep_pct = rep_sum / rep_total * 100 if rep_total > 0 else 0

    axes[1].bar([0], [act_pct], bottom=[bottom_act], color=FUNC_CAT_COLORS[cat],
                label=f'{cat} ({act_sum}/{rep_sum})', width=0.5, edgecolor='white', linewidth=0.5)
    axes[1].bar([1], [rep_pct], bottom=[bottom_rep], color=FUNC_CAT_COLORS[cat],
                width=0.5, edgecolor='white', linewidth=0.5)

    # Add text labels
    if act_pct > 5:
        axes[1].text(0, bottom_act + act_pct/2, f'{act_pct:.0f}%', ha='center', va='center', fontsize=8, fontweight='bold')
    if rep_pct > 5:
        axes[1].text(1, bottom_rep + rep_pct/2, f'{rep_pct:.0f}%', ha='center', va='center', fontsize=8, fontweight='bold')

    bottom_act += act_pct
    bottom_rep += rep_pct

axes[1].set_xticks([0, 1])
axes[1].set_xticklabels(bloc_names)
axes[1].set_ylabel('Percentage (%)')
axes[1].set_title('Functional Category Composition by Bloc')
axes[1].legend(loc='upper right', fontsize=7)
axes[1].set_ylim(0, 105)

plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f'{FIG_DIR}/TF_family_bloc_composition.{ext}')
plt.close()
print("  Saved TF_family_bloc_composition.pdf/svg")

# --- Figure 2: Functional category enrichment ---
fig, ax = plt.subplots(figsize=(10, 5))

cat_data = cat_enrich_df.sort_values('odds_ratio', ascending=True)
y_pos = np.arange(len(cat_data))

colors = [FUNC_CAT_COLORS.get(cat, '#95A5A6') for cat in cat_data['functional_category']]
bars = ax.barh(y_pos, np.log2(cat_data['odds_ratio'].clip(lower=0.01)), color=colors, edgecolor='white', alpha=0.85)

ax.set_yticks(y_pos)
ax.set_yticklabels(cat_data['functional_category'])
ax.set_xlabel('log2(Odds Ratio)\n<-- Repression enriched | Activation enriched -->')
ax.set_title('Functional Category Enrichment by Bloc (Fisher Exact Test)')
ax.axvline(0, color='black', linestyle='-', linewidth=0.5)

# Add p-value annotations
for i, (_, row) in enumerate(cat_data.iterrows()):
    or_val = row['odds_ratio']
    p_val = row['p_value']
    sig = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else 'ns'
    x_pos = np.log2(max(or_val, 0.01))
    offset = 0.1 if x_pos >= 0 else -0.1
    ax.text(x_pos + offset, i, f'p={p_val:.3f} {sig}', va='center', fontsize=8)

plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f'{FIG_DIR}/functional_category_enrichment.{ext}')
plt.close()
print("  Saved functional_category_enrichment.pdf/svg")

# --- Figure 3: Genomic context analysis ---
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Left: Neighbor category bar chart
nb_cats_display = [c for c in all_neighbor_cats if c != 'other']
act_nb = [neighbor_comparison.loc[c, 'activation_pct'] if c in neighbor_comparison.index else 0 for c in nb_cats_display]
rep_nb = [neighbor_comparison.loc[c, 'repression_pct'] if c in neighbor_comparison.index else 0 for c in nb_cats_display]

x = np.arange(len(nb_cats_display))
width = 0.35

axes[0].barh(x - width/2, act_nb, width, label='Activation', color=BLOC_COLORS['Activation'], alpha=0.85)
axes[0].barh(x + width/2, rep_nb, width, label='Repression', color=BLOC_COLORS['Repression'], alpha=0.85)
axes[0].set_yticks(x)
axes[0].set_yticklabels([c.replace('_', ' ') for c in nb_cats_display], fontsize=8)
axes[0].set_xlabel('Percentage of neighbors (%)')
axes[0].set_title('Neighbor Gene Function Profiles by Bloc')
axes[0].legend()

# Right: Per-TF neighbor category heatmap
# Count per-TF per-category
per_tf_nb = neighbor_df.groupby(['tf_locus_tag', 'neighbor_category']).size().unstack(fill_value=0)
if 'other' in per_tf_nb.columns:
    per_tf_nb = per_tf_nb.drop('other', axis=1)

# Merge with bloc info
per_tf_nb_merged = per_tf_nb.merge(exposed[['locus_tag', 'bloc', 'old_locus_tag']].set_index('locus_tag'),
                                    left_index=True, right_index=True)
per_tf_nb_merged = per_tf_nb_merged.sort_values(['bloc', 'old_locus_tag'])

# Heatmap
display_cols = [c for c in per_tf_nb.columns if c != 'other']
heatmap_data = per_tf_nb_merged[display_cols].values
labels_y = per_tf_nb_merged['old_locus_tag'].fillna(pd.Series(per_tf_nb_merged.index, index=per_tf_nb_merged.index)).values
blocs_y = per_tf_nb_merged['bloc'].values

im = axes[1].imshow(heatmap_data, aspect='auto', cmap='YlOrRd', interpolation='none')
axes[1].set_xticks(range(len(display_cols)))
axes[1].set_xticklabels([c.replace('_', '\n') for c in display_cols], fontsize=7, rotation=45, ha='right')
axes[1].set_yticks(range(len(labels_y)))
axes[1].set_yticklabels(labels_y, fontsize=5)
axes[1].set_title('Neighbor Gene Categories per Exposed TF')
plt.colorbar(im, ax=axes[1], label='Count')

# Add bloc separator
bloc_change_idx = np.where(np.array(blocs_y[:-1]) != np.array(blocs_y[1:]))[0]
for idx in bloc_change_idx:
    axes[1].axhline(idx + 0.5, color='black', linewidth=2)

# Add bloc labels
act_mid = np.mean(np.where(np.array(blocs_y) == 'Activation')[0])
rep_mid = np.mean(np.where(np.array(blocs_y) == 'Repression')[0])
axes[1].text(-1.5, act_mid, 'ACT', fontsize=8, fontweight='bold', color=BLOC_COLORS['Activation'], ha='center')
axes[1].text(-1.5, rep_mid, 'REP', fontsize=8, fontweight='bold', color=BLOC_COLORS['Repression'], ha='center')

plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f'{FIG_DIR}/genomic_context_analysis.{ext}')
plt.close()
print("  Saved genomic_context_analysis.pdf/svg")

# --- Figure 4: Expression by family ---
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Left: Box plot of |LFC| by family
families_plot = [f for f in family_counts.index if family_counts[f] >= 2]
family_data = [exposed[exposed['tf_family'] == f]['abs_LFC_T3'].values for f in families_plot]

bp = axes[0].boxplot(family_data, labels=families_plot, vert=True, patch_artist=True)
for i, f in enumerate(families_plot):
    cat = TF_FAMILY_INFO.get(f, {}).get('functional_category', 'General')
    bp['boxes'][i].set_facecolor(FUNC_CAT_COLORS.get(cat, '#95A5A6'))
    bp['boxes'][i].set_alpha(0.7)

# Overlay individual points
for i, f in enumerate(families_plot):
    data = exposed[exposed['tf_family'] == f]
    jitter = np.random.normal(0, 0.05, len(data))
    colors = [BLOC_COLORS[b] for b in data['bloc']]
    axes[0].scatter(np.full(len(data), i+1) + jitter, data['abs_LFC_T3'],
                   c=colors, s=25, alpha=0.8, edgecolors='black', linewidth=0.5, zorder=3)

axes[0].set_xlabel('TF Family')
axes[0].set_ylabel('|log2FC T3 vs T1|')
axes[0].set_title(f'Expression Response by TF Family\n(KW H={kw_stat:.2f}, p={kw_p:.4f})')
axes[0].tick_params(axis='x', rotation=45)

# Add legend for bloc colors
act_patch = mpatches.Patch(color=BLOC_COLORS['Activation'], label='Activation bloc')
rep_patch = mpatches.Patch(color=BLOC_COLORS['Repression'], label='Repression bloc')
axes[0].legend(handles=[act_patch, rep_patch], loc='upper right')

# Right: Ranked |LFC| with family color
exposed_sorted = exposed.sort_values('abs_LFC_T3', ascending=True)
y_pos = np.arange(len(exposed_sorted))
colors = [FUNC_CAT_COLORS.get(TF_FAMILY_INFO.get(f, {}).get('functional_category', 'General'), '#95A5A6')
          for f in exposed_sorted['tf_family']]
edge_colors = [BLOC_COLORS[b] for b in exposed_sorted['bloc']]

axes[1].barh(y_pos, exposed_sorted['abs_LFC_T3'], color=colors, edgecolor=edge_colors, linewidth=1.5, alpha=0.85)
axes[1].set_yticks(y_pos[::3])
axes[1].set_yticklabels(exposed_sorted['old_locus_tag'].iloc[::3], fontsize=6)
axes[1].set_xlabel('|log2FC T3 vs T1|')
axes[1].set_title('Ranked Expression Response of 57 Exposed TFs')

# Legend
from matplotlib.lines import Line2D
legend_elements = [Line2D([0], [0], color=c, lw=4, label=cat) for cat, c in FUNC_CAT_COLORS.items()]
axes[1].legend(handles=legend_elements, loc='lower right', fontsize=7, title='Functional category')

plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f'{FIG_DIR}/expression_by_family.{ext}')
plt.close()
print("  Saved expression_by_family.pdf/svg")

# --- Figure 5: Predicted regulatory model ---
fig, ax = plt.subplots(figsize=(16, 10))
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.axis('off')
ax.set_title('H35: Predicted Regulatory Model for Exposed TF Two-Bloc System', fontsize=14, fontweight='bold', pad=20)

# Title banner
ax.add_patch(plt.Rectangle((0.5, 9.2), 9, 0.6, facecolor='#2C3E50', alpha=0.9, zorder=2))
ax.text(5, 9.5, 'Distributed Methylation-Responsive Regulatory Layer',
        ha='center', va='center', fontsize=13, fontweight='bold', color='white', zorder=3)

# Activation bloc box
ax.add_patch(plt.Rectangle((0.5, 4.5), 4.0, 4.4, facecolor=BLOC_COLORS['Activation'], alpha=0.15,
                            edgecolor=BLOC_COLORS['Activation'], linewidth=2, zorder=1))
ax.text(2.5, 8.7, 'ACTIVATION BLOC', ha='center', va='center', fontsize=12, fontweight='bold',
        color=BLOC_COLORS['Activation'])
ax.text(2.5, 8.4, '(35 genes, Modules 1-3)', ha='center', va='center', fontsize=9, color='#666')

# Activation predictions
act_items = [
    ('HIGH', 'Morphological differentiation', 'RamR, Sigma factors, WhiB'),
    ('HIGH', 'TCS signal transduction', '4 SK + 3 RR activated'),
    ('MED', 'Secondary metabolite regulation', 'TcrA, NsdB, TetR'),
    ('MED', 'Stress response activation', 'SigE, SigB/F/G, MerR'),
    ('LOW', 'Sugar/nucleotide metabolism', 'DeoR, LacI')
]

for i, (conf, process, evidence) in enumerate(act_items):
    y = 7.8 - i * 0.65
    conf_color = {'HIGH': '#27AE60', 'MED': '#F39C12', 'LOW': '#E74C3C'}[conf]
    ax.add_patch(plt.Rectangle((0.7, y - 0.2), 0.5, 0.35, facecolor=conf_color, alpha=0.8))
    ax.text(0.95, y - 0.025, conf[:1], ha='center', va='center', fontsize=7, fontweight='bold', color='white')
    ax.text(1.3, y + 0.05, process, ha='left', va='center', fontsize=8, fontweight='bold')
    ax.text(1.3, y - 0.2, evidence, ha='left', va='center', fontsize=7, color='#555', style='italic')

# Repression bloc box
ax.add_patch(plt.Rectangle((5.5, 2.0), 4.0, 6.9, facecolor=BLOC_COLORS['Repression'], alpha=0.15,
                            edgecolor=BLOC_COLORS['Repression'], linewidth=2, zorder=1))
ax.text(7.5, 8.7, 'REPRESSION BLOC', ha='center', va='center', fontsize=12, fontweight='bold',
        color=BLOC_COLORS['Repression'])
ax.text(7.5, 8.4, '(26 genes, Module 4)', ha='center', va='center', fontsize=9, color='#666')

rep_items = [
    ('HIGH', 'Primary metabolism shutdown', 'GntR, IclR, LacI, LysR'),
    ('HIGH', 'Defense/efflux downregulation', '5 TetR repressed'),
    ('HIGH', 'Fatty acid biosynthesis off', 'FasR (SCO2386)'),
    ('HIGH', 'DNA replication/repair off', 'SSB, HU, Mfd, UdgX'),
    ('MED', 'Oxidative stress response off', 'MarR (SCO1191, SCO4122)'),
    ('MED', 'Select TCS downregulation', 'SCO5289 SK, SCO5434 RR'),
    ('MED', 'Sigma factor repression', 'SCO2639, SCO4960'),
    ('MED', 'Nitrogen assimilation off', 'P-II (SCO5584)')
]

for i, (conf, process, evidence) in enumerate(rep_items):
    y = 7.8 - i * 0.65
    conf_color = {'HIGH': '#27AE60', 'MED': '#F39C12', 'LOW': '#E74C3C'}[conf]
    ax.add_patch(plt.Rectangle((5.7, y - 0.2), 0.5, 0.35, facecolor=conf_color, alpha=0.8))
    ax.text(5.95, y - 0.025, conf[:1], ha='center', va='center', fontsize=7, fontweight='bold', color='white')
    ax.text(6.3, y + 0.05, process, ha='left', va='center', fontsize=8, fontweight='bold')
    ax.text(6.3, y - 0.2, evidence, ha='left', va='center', fontsize=7, color='#555', style='italic')

# Central annotation
ax.annotate('', xy=(5.4, 7.5), xytext=(4.6, 7.5),
            arrowprops=dict(arrowstyle='<->', lw=2, color='#2C3E50'))
ax.text(5.0, 7.8, 'Antagonistic', ha='center', va='center', fontsize=8, fontweight='bold', color='#2C3E50')

# Bottom summary box
ax.add_patch(plt.Rectangle((0.5, 0.3), 9, 1.5, facecolor='#ECF0F1', edgecolor='#BDC3C7', linewidth=1))
ax.text(5, 1.6, 'DEVELOPMENTAL TRANSITION MODEL (T1 -> T3)', ha='center', va='center', fontsize=11, fontweight='bold')
ax.text(5, 1.15, 'Vegetative growth programs (metabolism, DNA replication, defense) are repressed',
        ha='center', va='center', fontsize=9, color=BLOC_COLORS['Repression'])
ax.text(5, 0.8, 'Developmental programs (morphogenesis, TCS signaling, stress responses, BGC regulation) are activated',
        ha='center', va='center', fontsize=9, color=BLOC_COLORS['Activation'])
ax.text(5, 0.45, 'Confidence: H = High (multiple lines of evidence)  |  M = Medium (family + context)  |  L = Low (family alone)',
        ha='center', va='center', fontsize=8, color='#666')

plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f'{FIG_DIR}/predicted_regulatory_model.{ext}')
plt.close()
print("  Saved predicted_regulatory_model.pdf/svg")

# --- Figure 6: Comprehensive summary ---
fig = plt.figure(figsize=(20, 16))
gs = GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.3)

# Panel A: Family composition
ax_a = fig.add_subplot(gs[0, 0])
act_counts_p = [bloc_family_counts.loc[f, 'Activation'] for f in families_sorted]
rep_counts_p = [bloc_family_counts.loc[f, 'Repression'] for f in families_sorted]
x = np.arange(len(families_sorted))
width = 0.35
ax_a.bar(x - width/2, act_counts_p, width, label='Activation', color=BLOC_COLORS['Activation'], alpha=0.85)
ax_a.bar(x + width/2, rep_counts_p, width, label='Repression', color=BLOC_COLORS['Repression'], alpha=0.85)
ax_a.set_xticks(x)
ax_a.set_xticklabels(families_sorted, rotation=60, ha='right', fontsize=6)
ax_a.set_ylabel('Count')
ax_a.set_title('A. TF Family by Bloc', fontweight='bold')
ax_a.legend(fontsize=7)

# Panel B: Functional category stacked bar
ax_b = fig.add_subplot(gs[0, 1])
bottom_act = 0
bottom_rep = 0
for cat in func_cats_ordered:
    cat_fams = [f for f in sorted_fams if fam_func_map[f] == cat]
    act_sum = sum(bloc_family_counts.loc[f, 'Activation'] for f in cat_fams)
    rep_sum = sum(bloc_family_counts.loc[f, 'Repression'] for f in cat_fams)
    act_pct = act_sum / act_total * 100
    rep_pct = rep_sum / rep_total * 100
    ax_b.bar([0], [act_pct], bottom=[bottom_act], color=FUNC_CAT_COLORS[cat], label=cat, width=0.5, edgecolor='white')
    ax_b.bar([1], [rep_pct], bottom=[bottom_rep], color=FUNC_CAT_COLORS[cat], width=0.5, edgecolor='white')
    bottom_act += act_pct
    bottom_rep += rep_pct
ax_b.set_xticks([0, 1])
ax_b.set_xticklabels(['Activation', 'Repression'])
ax_b.set_ylabel('Percentage (%)')
ax_b.set_title('B. Functional Categories', fontweight='bold')
ax_b.legend(fontsize=6, loc='upper right')

# Panel C: Fisher enrichment OR
ax_c = fig.add_subplot(gs[0, 2])
fisher_sorted = fisher_df.sort_values('odds_ratio', ascending=True)
fisher_sig = fisher_sorted[fisher_sorted['activation_count'] + fisher_sorted['repression_count'] >= 2]
if len(fisher_sig) > 0:
    y_pos = np.arange(len(fisher_sig))
    colors = [FUNC_CAT_COLORS.get(TF_FAMILY_INFO.get(f, {}).get('functional_category', 'General'), '#95A5A6')
              for f in fisher_sig['tf_family']]
    or_log2 = np.log2(fisher_sig['odds_ratio'].clip(lower=0.01).clip(upper=100))
    ax_c.barh(y_pos, or_log2, color=colors, alpha=0.85, edgecolor='white')
    ax_c.set_yticks(y_pos)
    ax_c.set_yticklabels(fisher_sig['tf_family'], fontsize=7)
    ax_c.axvline(0, color='black', linewidth=0.5)
    ax_c.set_xlabel('log2(OR)')
    for i, (_, row) in enumerate(fisher_sig.iterrows()):
        if row['p_value'] < 0.05:
            ax_c.text(or_log2.iloc[i], i, ' *', va='center', fontsize=10, color='red')
ax_c.set_title('C. Family Enrichment (Fisher)', fontweight='bold')

# Panel D: |LFC| by family box plot
ax_d = fig.add_subplot(gs[1, 0:2])
fams_for_box = [f for f in family_counts.index]
data_for_box = [exposed[exposed['tf_family'] == f]['abs_LFC_T3'].values for f in fams_for_box]
bp = ax_d.boxplot(data_for_box, labels=fams_for_box, vert=True, patch_artist=True)
for i, f in enumerate(fams_for_box):
    cat = TF_FAMILY_INFO.get(f, {}).get('functional_category', 'General')
    bp['boxes'][i].set_facecolor(FUNC_CAT_COLORS.get(cat, '#95A5A6'))
    bp['boxes'][i].set_alpha(0.7)
    # Individual points
    d = exposed[exposed['tf_family'] == f]
    jitter = np.random.normal(0, 0.08, len(d))
    cols = [BLOC_COLORS[b] for b in d['bloc']]
    ax_d.scatter(np.full(len(d), i+1) + jitter, d['abs_LFC_T3'], c=cols, s=15, alpha=0.8,
                edgecolors='black', linewidth=0.3, zorder=3)
ax_d.set_ylabel('|log2FC T3 vs T1|')
ax_d.set_title(f'D. Expression Response by TF Family (KW p={kw_p:.4f})', fontweight='bold')
ax_d.tick_params(axis='x', rotation=45)

# Panel E: Neighbor context radar-like comparison
ax_e = fig.add_subplot(gs[1, 2])
nb_cats_radar = ['secondary_metabolism', 'transport', 'signaling', 'primary_metabolism',
                 'stress_response', 'hypothetical']
act_vals = [neighbor_comparison.loc[c, 'activation_pct'] if c in neighbor_comparison.index else 0 for c in nb_cats_radar]
rep_vals = [neighbor_comparison.loc[c, 'repression_pct'] if c in neighbor_comparison.index else 0 for c in nb_cats_radar]

x_radar = np.arange(len(nb_cats_radar))
ax_e.barh(x_radar - 0.2, act_vals, 0.4, label='Activation', color=BLOC_COLORS['Activation'], alpha=0.85)
ax_e.barh(x_radar + 0.2, rep_vals, 0.4, label='Repression', color=BLOC_COLORS['Repression'], alpha=0.85)
ax_e.set_yticks(x_radar)
ax_e.set_yticklabels([c.replace('_', '\n') for c in nb_cats_radar], fontsize=7)
ax_e.set_xlabel('% of neighbors')
ax_e.set_title('E. Neighbor Functions', fontweight='bold')
ax_e.legend(fontsize=7)

# Panel F: Summary text
ax_f = fig.add_subplot(gs[2, :])
ax_f.axis('off')

summary_text = (
    "H35 SUMMARY: Functional Prediction of 57 Exposed Transcription Factors\n\n"
    "ACTIVATION BLOC (35 genes, Modules 1-3): Enriched for Signal transduction (TCS) and Stress/development (sigma factors, WhiB).\n"
    "Predicted activated processes: Morphological differentiation (RamR, sigma factors), TCS cascades, secondary metabolite regulation (TcrA, NsdB).\n\n"
    "REPRESSION BLOC (26 genes, Module 4): Enriched for Metabolic regulators (GntR, IclR, LacI, LysR) and Defense/resistance (TetR, ArsR).\n"
    "Predicted repressed processes: Primary metabolism shutdown, antibiotic efflux/defense downregulation, DNA replication/repair cessation,\n"
    "fatty acid biosynthesis (FasR), nitrogen assimilation (P-II).\n\n"
    "MODEL: The two blocs implement a coordinated developmental switch -- vegetative growth programs are repressed while\n"
    "developmental/stress/signaling programs are activated, consistent with the S. coelicolor life cycle transition."
)
ax_f.text(0.5, 0.5, summary_text, ha='center', va='center', fontsize=10, fontfamily='monospace',
          bbox=dict(boxstyle='round', facecolor='#ECF0F1', alpha=0.8), transform=ax_f.transAxes)

fig.suptitle('H35: Comprehensive Summary - Exposed TF Functional Prediction', fontsize=14, fontweight='bold', y=0.98)
plt.tight_layout(rect=[0, 0, 1, 0.97])
for ext in ['pdf', 'svg']:
    fig.savefig(f'{FIG_DIR}/H35_comprehensive_summary.{ext}')
plt.close()
print("  Saved H35_comprehensive_summary.pdf/svg")

# ============================================================
# Step 10: Save tables
# ============================================================
print("\n" + "=" * 70)
print("STEP 10: Saving tables")
print("=" * 70)

# Table 1: Per-TF family annotation
tf_annot_export = tf_annot_df.copy()
# Merge predicted scope
scope_merge = scope_df[['locus_tag', 'predicted_regulatory_scope', 'dominant_neighbor_category']]
tf_annot_export = tf_annot_export.merge(scope_merge, on='locus_tag', how='left')
tf_annot_export.to_csv(f'{TAB_DIR}/exposed_TF_family_annotation.tsv', sep='\t', index=False)
print(f"  Saved exposed_TF_family_annotation.tsv ({len(tf_annot_export)} rows)")

# Table 2: Bloc family enrichment
fisher_df.to_csv(f'{TAB_DIR}/bloc_family_enrichment.tsv', sep='\t', index=False)
print(f"  Saved bloc_family_enrichment.tsv ({len(fisher_df)} rows)")

# Table 3: Genomic context neighbors
neighbor_df.to_csv(f'{TAB_DIR}/genomic_context_neighbors.tsv', sep='\t', index=False)
print(f"  Saved genomic_context_neighbors.tsv ({len(neighbor_df)} rows)")

# Table 4: Functional category mapping
func_cat_map_rows = []
for fam, info in TF_FAMILY_INFO.items():
    n_exposed = len(exposed[exposed['tf_family'] == fam])
    if n_exposed > 0:
        func_cat_map_rows.append({
            'tf_family': fam,
            'functional_category': info['functional_category'],
            'full_name': info['full_name'],
            'general_functions': info['general_functions'],
            'streptomyces_roles': info['streptomyces_roles'],
            'typical_regulon_size': info['typical_regulon_size'],
            'target_categories': info['target_categories'],
            'n_exposed': n_exposed
        })

func_cat_map_df = pd.DataFrame(func_cat_map_rows)
func_cat_map_df.to_csv(f'{TAB_DIR}/functional_category_mapping.tsv', sep='\t', index=False)
print(f"  Saved functional_category_mapping.tsv ({len(func_cat_map_df)} rows)")

# Table 5: Expression by family
expr_stats_df.to_csv(f'{TAB_DIR}/expression_by_family.tsv', sep='\t', index=False)
print(f"  Saved expression_by_family.tsv ({len(expr_stats_df)} rows)")

# Table 6: Statistical tests
stat_tests = []

# Fisher tests - family
for _, row in fisher_df.iterrows():
    stat_tests.append({
        'test_name': f'Fisher exact: {row["tf_family"]} enrichment by bloc',
        'test_type': 'Fisher exact test',
        'statistic': f'OR={row["odds_ratio"]:.4f}',
        'p_value': row['p_value'],
        'p_adj': row['p_adj'],
        'n1': row['activation_count'],
        'n2': row['repression_count'],
        'result': row['enriched_in']
    })

# Fisher tests - functional category
for _, row in cat_enrich_df.iterrows():
    stat_tests.append({
        'test_name': f'Fisher exact: {row["functional_category"]} enrichment by bloc',
        'test_type': 'Fisher exact test',
        'statistic': f'OR={row["odds_ratio"]:.4f}',
        'p_value': row['p_value'],
        'p_adj': row['p_adj'],
        'n1': row['activation_count'],
        'n2': row['repression_count'],
        'result': row['enriched_in']
    })

# Kruskal-Wallis
stat_tests.append({
    'test_name': 'KW: |LFC_T3| differences by TF family',
    'test_type': 'Kruskal-Wallis H test',
    'statistic': f'H={kw_stat:.4f}',
    'p_value': kw_p,
    'p_adj': np.nan,
    'n1': len(families_with_enough),
    'n2': len(exposed),
    'result': 'Significant' if kw_p < 0.05 else 'Non-significant'
})

# MWU bloc |LFC|
stat_tests.append({
    'test_name': 'MWU: |LFC_T3| activation vs repression bloc',
    'test_type': 'Mann-Whitney U test',
    'statistic': f'U={mw_stat:.1f}',
    'p_value': mw_p,
    'p_adj': np.nan,
    'n1': len(act_lfc),
    'n2': len(rep_lfc),
    'result': f'Act median={np.median(act_lfc):.2f}, Rep median={np.median(rep_lfc):.2f}'
})

# Neighbor category Fisher tests
for _, row in neighbor_fisher_df.iterrows():
    stat_tests.append({
        'test_name': f'Fisher exact: neighbor {row["category"]} by bloc',
        'test_type': 'Fisher exact test',
        'statistic': f'OR={row["odds_ratio"]:.4f}',
        'p_value': row['p_value'],
        'p_adj': row['p_adj'],
        'n1': row['activation_count'],
        'n2': row['repression_count'],
        'result': 'Activation enriched' if row['odds_ratio'] > 1 else 'Repression enriched'
    })

stat_tests_df = pd.DataFrame(stat_tests)
stat_tests_df.to_csv(f'{TAB_DIR}/statistical_tests.tsv', sep='\t', index=False)
print(f"  Saved statistical_tests.tsv ({len(stat_tests_df)} rows)")

# ============================================================
# Final summary
# ============================================================
print("\n" + "=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)

print(f"\n57 exposed TFs analyzed across {len(all_families)} TF families")
print(f"Activation bloc: {n_act} genes (modules 1-3)")
print(f"Repression bloc: {n_rep} genes (module 4)")

print("\n--- Key findings ---")
# Functional category summary
act_stress_dev = len(activation[activation['functional_category'] == 'Stress/development'])
rep_stress_dev = len(repression[repression['functional_category'] == 'Stress/development'])
act_metabolic = len(activation[activation['functional_category'] == 'Metabolic'])
rep_metabolic = len(repression[repression['functional_category'] == 'Metabolic'])
act_signal = len(activation[activation['functional_category'] == 'Signal transduction'])
rep_signal = len(repression[repression['functional_category'] == 'Signal transduction'])
act_defense = len(activation[activation['functional_category'] == 'Defense/resistance'])
rep_defense = len(repression[repression['functional_category'] == 'Defense/resistance'])

print(f"\n1. Stress/development: Activation={act_stress_dev} ({act_stress_dev/n_act*100:.1f}%) vs Repression={rep_stress_dev} ({rep_stress_dev/n_rep*100:.1f}%)")
print(f"2. Signal transduction: Activation={act_signal} ({act_signal/n_act*100:.1f}%) vs Repression={rep_signal} ({rep_signal/n_rep*100:.1f}%)")
print(f"3. Metabolic: Activation={act_metabolic} ({act_metabolic/n_act*100:.1f}%) vs Repression={rep_metabolic} ({rep_metabolic/n_rep*100:.1f}%)")
print(f"4. Defense/resistance: Activation={act_defense} ({act_defense/n_act*100:.1f}%) vs Repression={rep_defense} ({rep_defense/n_rep*100:.1f}%)")

# Significant Fisher tests
sig_fisher = fisher_df[fisher_df['p_value'] < 0.05]
print(f"\nSignificant family enrichments (p<0.05): {len(sig_fisher)}")
for _, row in sig_fisher.iterrows():
    print(f"  {row['tf_family']}: OR={row['odds_ratio']:.2f}, p={row['p_value']:.4f} [{row['enriched_in']}]")

sig_cat = cat_enrich_df[cat_enrich_df['p_value'] < 0.05]
print(f"\nSignificant category enrichments (p<0.05): {len(sig_cat)}")
for _, row in sig_cat.iterrows():
    print(f"  {row['functional_category']}: OR={row['odds_ratio']:.2f}, p={row['p_value']:.4f} [{row['enriched_in']}]")

print(f"\nKruskal-Wallis family |LFC| difference: p={kw_p:.4f}")
print(f"Bloc |LFC| difference: Act median={np.median(act_lfc):.2f} vs Rep median={np.median(rep_lfc):.2f}, p={mw_p:.4f}")

print("\n--- Files generated ---")
print(f"Figures: {FIG_DIR}/")
for f in sorted(os.listdir(FIG_DIR)):
    print(f"  {f}")
print(f"\nTables: {TAB_DIR}/")
for f in sorted(os.listdir(TAB_DIR)):
    print(f"  {f}")

print("\nDone.")
