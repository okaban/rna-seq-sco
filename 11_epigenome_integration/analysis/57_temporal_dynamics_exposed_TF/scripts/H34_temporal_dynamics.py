#!/usr/bin/env python3
"""
H34: Temporal Dynamics and Hierarchical Structure of 62 Exposed Transcription Factors
=====================================================================================
Analyzes phase separation between activation and repression blocs,
temporal ordering of TCS pairs, and methylation-expression timing correlations.

Author: Claude Code (H34 analysis)
Date: 2026-02-27
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyArrowPatch
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import pdist
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# Configuration
# ============================================================================
BASE = '/Users/okaban/bioinfo/rna-seq'
ANALYSIS_DIR = f'{BASE}/11_epigenome_integration/analysis/57_temporal_dynamics_exposed_TF'
FIG_DIR = f'{ANALYSIS_DIR}/figures'
TBL_DIR = f'{ANALYSIS_DIR}/tables'

# Input files
EXPOSED_FILE = f'{BASE}/11_epigenome_integration/analysis/51_exposed_regulators_characteristics/tables/exposed_regulators_full_table.tsv'
ALL_GENES_FILE = f'{BASE}/11_epigenome_integration/analysis/52_shielded_exposed_boundary/tables/all_genes_features.tsv'
DESEQ_T2vsT1 = f'{BASE}/04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_2_vs_1.tsv'
DESEQ_T3vsT1 = f'{BASE}/04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_3_vs_1.tsv'
DESEQ_T3vsT2 = f'{BASE}/04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_3_vs_2.tsv'
NORM_COUNTS = f'{BASE}/04_deseq2/analysis/04_deseq2_260128_v1/results/normalized_counts_M145.tsv'
COORDINATED_FILE = f'{BASE}/11_epigenome_integration/analysis/29_genomewide_TF_screen/tables/coordinated_regulatory_genes.tsv'
COEXPR_MODULES = f'{BASE}/11_epigenome_integration/analysis/55_exposed_regulatory_module/tables/coexpression_modules.tsv'
TCS_PAIRS_FILE = f'{BASE}/11_epigenome_integration/analysis/55_exposed_regulatory_module/tables/TCS_pairs_analysis.tsv'
AAGCCCG_FILE = f'{BASE}/11_epigenome_integration/analysis/36_AAGCCCG_distribution/tables/AAGCCCG_site_gene_mapping.tsv'
GCCGGC_FILE = f'{BASE}/11_epigenome_integration/analysis/37_defense_island_GCCGGC/tables/GCCGGC_sites_by_timepoint.tsv'

# Normalized counts sample columns
T1_COLS = ['M145_1_1', 'M145_1_2', 'M145_1_3']
T2_COLS = ['M145_2_1', 'M145_2_3', 'M145_2_4']
T3_COLS = ['M145_3_2', 'M145_3_3', 'M145_3_4']

# Plot style
plt.rcParams.update({
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

# Color palette
COLORS = {
    'early': '#e74c3c',
    'late': '#3498db',
    'gradual': '#9b59b6',
    'activation': '#2ecc71',
    'repression': '#e67e22',
    'up': '#c0392b',
    'down': '#2980b9',
    'exposed': '#e74c3c',
    'shielded': '#95a5a6',
    'module1': '#e74c3c',
    'module2': '#f39c12',
    'module3': '#2ecc71',
    'module4': '#3498db',
}

# ============================================================================
# Step 1: Load data
# ============================================================================
print("=" * 70)
print("STEP 1: Loading data")
print("=" * 70)

exposed_df = pd.read_csv(EXPOSED_FILE, sep='\t')
print(f"  Exposed regulators: {len(exposed_df)} genes")

all_genes_df = pd.read_csv(ALL_GENES_FILE, sep='\t')
print(f"  All regulatory genes: {len(all_genes_df)} genes")

deseq_t2t1 = pd.read_csv(DESEQ_T2vsT1, sep='\t')
deseq_t3t1 = pd.read_csv(DESEQ_T3vsT1, sep='\t')
deseq_t3t2 = pd.read_csv(DESEQ_T3vsT2, sep='\t')
print(f"  DESeq2 T2vsT1: {len(deseq_t2t1)} genes")
print(f"  DESeq2 T3vsT1: {len(deseq_t3t1)} genes")
print(f"  DESeq2 T3vsT2: {len(deseq_t3t2)} genes")

norm_counts = pd.read_csv(NORM_COUNTS, sep='\t')
print(f"  Normalized counts: {len(norm_counts)} genes x {norm_counts.shape[1]-1} samples")

coordinated_df = pd.read_csv(COORDINATED_FILE, sep='\t')
print(f"  Coordinated regulatory genes: {len(coordinated_df)} genes")

coexpr_modules = pd.read_csv(COEXPR_MODULES, sep='\t')
print(f"  Co-expression modules: {len(coexpr_modules)} modules")

tcs_pairs = pd.read_csv(TCS_PAIRS_FILE, sep='\t')
print(f"  TCS pairs: {len(tcs_pairs)} pairs")

aagcccg_sites = pd.read_csv(AAGCCCG_FILE, sep='\t')
print(f"  AAGCCCG sites: {len(aagcccg_sites)} site-gene mappings")

gccggc_sites = pd.read_csv(GCCGGC_FILE, sep='\t')
print(f"  GCCGGC sites: {len(gccggc_sites)} site entries")

# ============================================================================
# Parse module memberships
# ============================================================================
module_membership = {}
for _, row in coexpr_modules.iterrows():
    mod_id = row['module_id']
    tags = [t.strip() for t in str(row['locus_tags']).split(',')]
    for tag in tags:
        module_membership[tag] = mod_id

# Merge DESeq2 results for exposed TFs
exposed_tags = set(exposed_df['locus_tag'])
print(f"\n  Exposed locus_tags: {len(exposed_tags)}")

# Build master table for all 62 exposed TFs
master = exposed_df[['locus_tag', 'gene_name', 'old_locus_tag', 'product', 'tf_family',
                      'start', 'end', 'strand']].copy()

# Add DESeq2 LFC values
deseq_t2t1_map = deseq_t2t1.set_index('gene_id')['log2FoldChange'].to_dict()
deseq_t3t1_map = deseq_t3t1.set_index('gene_id')['log2FoldChange'].to_dict()
deseq_t3t2_map = deseq_t3t2.set_index('gene_id')['log2FoldChange'].to_dict()

padj_t2t1_map = deseq_t2t1.set_index('gene_id')['padj'].to_dict()
padj_t3t1_map = deseq_t3t1.set_index('gene_id')['padj'].to_dict()
padj_t3t2_map = deseq_t3t2.set_index('gene_id')['padj'].to_dict()

master['LFC_T2vsT1'] = master['locus_tag'].map(deseq_t2t1_map)
master['LFC_T3vsT1'] = master['locus_tag'].map(deseq_t3t1_map)
master['LFC_T3vsT2'] = master['locus_tag'].map(deseq_t3t2_map)
master['padj_T2vsT1'] = master['locus_tag'].map(padj_t2t1_map)
master['padj_T3vsT1'] = master['locus_tag'].map(padj_t3t1_map)
master['padj_T3vsT2'] = master['locus_tag'].map(padj_t3t2_map)

# Add module membership
master['module'] = master['locus_tag'].map(module_membership)
master['bloc'] = master['module'].apply(lambda x: 'activation' if x in [1, 2, 3] else ('repression' if x == 4 else 'unassigned'))

# Add coordination type from coordinated_df (if available)
if 'coordination_T2' in coordinated_df.columns:
    coord_t2_map = coordinated_df.set_index('locus_tag')['coordination_T2'].to_dict()
    master['coordination_T2'] = master['locus_tag'].map(coord_t2_map)
else:
    master['coordination_T2'] = np.nan
if 'coordination_T3' in coordinated_df.columns:
    coord_t3_map = coordinated_df.set_index('locus_tag')['coordination_T3'].to_dict()
    master['coordination_T3'] = master['locus_tag'].map(coord_t3_map)
else:
    master['coordination_T3'] = np.nan

# Add normalized counts for trajectory
norm_counts_idx = norm_counts.set_index('gene_id')
for col in T1_COLS + T2_COLS + T3_COLS:
    master[col] = master['locus_tag'].map(norm_counts_idx[col].to_dict())

# Compute mean expression per timepoint
master['T1_mean'] = master[T1_COLS].mean(axis=1)
master['T2_mean'] = master[T2_COLS].mean(axis=1)
master['T3_mean'] = master[T3_COLS].mean(axis=1)

# Compute z-scored expression
for idx in master.index:
    vals = [master.loc[idx, 'T1_mean'], master.loc[idx, 'T2_mean'], master.loc[idx, 'T3_mean']]
    mu = np.mean(vals)
    sd = np.std(vals)
    if sd > 0:
        master.loc[idx, 'T1_z'] = (vals[0] - mu) / sd
        master.loc[idx, 'T2_z'] = (vals[1] - mu) / sd
        master.loc[idx, 'T3_z'] = (vals[2] - mu) / sd
    else:
        master.loc[idx, 'T1_z'] = 0
        master.loc[idx, 'T2_z'] = 0
        master.loc[idx, 'T3_z'] = 0

print(f"\n  Master table: {len(master)} exposed TFs with {master.columns.tolist()[:10]}...")
print(f"  Module distribution: {master['module'].value_counts().to_dict()}")
print(f"  Bloc distribution: {master['bloc'].value_counts().to_dict()}")

# ============================================================================
# Step 2: Classify exposed TFs by temporal response
# ============================================================================
print("\n" + "=" * 70)
print("STEP 2: Temporal classification of exposed TFs")
print("=" * 70)

# Early response: |LFC_T2vsT1|
master['early_magnitude'] = master['LFC_T2vsT1'].abs()
# Late response: |LFC_T3vsT2|
master['late_magnitude'] = master['LFC_T3vsT2'].abs()

# Phase ratio: LFC_T2vsT1 / LFC_T3vsT1
# Only for genes with meaningful total change (|LFC_T3vsT1| > 0.5)
LFC_THRESHOLD = 0.5
master['total_change'] = master['LFC_T3vsT1'].abs()
master['has_meaningful_change'] = master['total_change'] > LFC_THRESHOLD

def compute_phase_ratio(row):
    """Phase ratio = fraction of total change occurring in early transition."""
    lfc_t2t1 = row['LFC_T2vsT1']
    lfc_t3t1 = row['LFC_T3vsT1']
    if abs(lfc_t3t1) <= LFC_THRESHOLD:
        return np.nan
    # Phase ratio: how much of the T1->T3 change happened by T2
    ratio = lfc_t2t1 / lfc_t3t1
    # Clamp to [0, 1+] -- can be >1 if overshoot or <0 if reversal
    return ratio

master['phase_ratio'] = master.apply(compute_phase_ratio, axis=1)

# Temporal classification
def classify_temporal(row):
    pr = row['phase_ratio']
    if pd.isna(pr):
        return 'non-responder'
    elif pr > 0.6:
        return 'early'
    elif pr < 0.4:
        return 'late'
    else:
        return 'gradual'

master['temporal_class'] = master.apply(classify_temporal, axis=1)

# Direction
master['direction'] = master['LFC_T3vsT1'].apply(lambda x: 'up' if x > 0 else 'down')

# Combined classification
master['temporal_direction'] = master['temporal_class'] + '_' + master['direction']

print("\n  Temporal classification (n=62 exposed TFs):")
print(f"  Has meaningful change (|LFC_T3vsT1| > {LFC_THRESHOLD}): {master['has_meaningful_change'].sum()}")
print(f"  Phase ratio stats: mean={master['phase_ratio'].mean():.3f}, "
      f"median={master['phase_ratio'].median():.3f}, "
      f"std={master['phase_ratio'].std():.3f}")

class_counts = master['temporal_class'].value_counts()
print(f"\n  Classification counts:")
for cls, n in class_counts.items():
    print(f"    {cls}: {n}")

print(f"\n  2x3 classification (direction x temporal):")
cross = pd.crosstab(master['direction'], master['temporal_class'])
print(cross.to_string())

# ============================================================================
# Step 3: Bloc-level temporal analysis
# ============================================================================
print("\n" + "=" * 70)
print("STEP 3: Bloc-level temporal analysis")
print("=" * 70)

activation = master[master['bloc'] == 'activation']
repression = master[master['bloc'] == 'repression']

print(f"  Activation bloc: {len(activation)} genes (modules 1-3)")
print(f"  Repression bloc: {len(repression)} genes (module 4)")

# Phase ratio comparison
act_pr = activation['phase_ratio'].dropna()
rep_pr = repression['phase_ratio'].dropna()

print(f"\n  Activation phase ratio: mean={act_pr.mean():.3f}, median={act_pr.median():.3f} (n={len(act_pr)})")
print(f"  Repression phase ratio: mean={rep_pr.mean():.3f}, median={rep_pr.median():.3f} (n={len(rep_pr)})")

# Wilcoxon rank-sum test
if len(act_pr) > 0 and len(rep_pr) > 0:
    stat_w, p_w = stats.mannwhitneyu(act_pr, rep_pr, alternative='two-sided')
    # Effect size (rank-biserial correlation)
    n1, n2 = len(act_pr), len(rep_pr)
    r_rbs = 1 - (2 * stat_w) / (n1 * n2)
    print(f"  Mann-Whitney U test: U={stat_w:.1f}, p={p_w:.4e}, r_rbs={r_rbs:.3f}")
else:
    stat_w, p_w, r_rbs = np.nan, np.nan, np.nan
    print("  Cannot perform test: insufficient data")

# Mean temporal trajectory per bloc
bloc_traj = pd.DataFrame({
    'timepoint': ['T1', 'T2', 'T3'],
    'activation_mean_z': [activation['T1_z'].mean(), activation['T2_z'].mean(), activation['T3_z'].mean()],
    'activation_std_z': [activation['T1_z'].std(), activation['T2_z'].std(), activation['T3_z'].std()],
    'repression_mean_z': [repression['T1_z'].mean(), repression['T2_z'].mean(), repression['T3_z'].mean()],
    'repression_std_z': [repression['T1_z'].std(), repression['T2_z'].std(), repression['T3_z'].std()],
})
print(f"\n  Bloc trajectories (z-scored):")
print(bloc_traj.to_string(index=False))

# Temporal classification within blocs
print(f"\n  Temporal class x Bloc:")
print(pd.crosstab(master[master['bloc'].isin(['activation', 'repression'])]['bloc'],
                  master[master['bloc'].isin(['activation', 'repression'])]['temporal_class']).to_string())

# ============================================================================
# Step 4: Coordination type temporal patterns
# ============================================================================
print("\n" + "=" * 70)
print("STEP 4: Coordination type temporal patterns")
print("=" * 70)

# Use T3 coordination as primary (most change by T3)
coord_types = ['discordant_gain_up', 'concordant_derepression', 'concordant_repression', 'discordant_loss_down']
coord_results = []

for ct in coord_types:
    # Check both T2 and T3 coordination
    mask_t3 = master['coordination_T3'] == ct
    subset = master[mask_t3]
    pr_vals = subset['phase_ratio'].dropna()
    coord_results.append({
        'coordination_type': ct,
        'n_genes': len(subset),
        'n_with_phase_ratio': len(pr_vals),
        'mean_phase_ratio': pr_vals.mean() if len(pr_vals) > 0 else np.nan,
        'median_phase_ratio': pr_vals.median() if len(pr_vals) > 0 else np.nan,
        'std_phase_ratio': pr_vals.std() if len(pr_vals) > 0 else np.nan,
        'mean_LFC_T2vsT1': subset['LFC_T2vsT1'].mean(),
        'mean_LFC_T3vsT2': subset['LFC_T3vsT2'].mean(),
        'dominant_temporal_class': subset['temporal_class'].mode().iloc[0] if len(subset) > 0 else 'NA',
    })
    print(f"  {ct}: n={len(subset)}, phase_ratio mean={pr_vals.mean():.3f}" if len(pr_vals) > 0 else f"  {ct}: n={len(subset)}, no phase ratio data")

coord_results_df = pd.DataFrame(coord_results)

# Kruskal-Wallis test across coordination types (T3)
groups_for_kw = []
group_labels = []
for ct in coord_types:
    mask = master['coordination_T3'] == ct
    vals = master.loc[mask, 'phase_ratio'].dropna()
    if len(vals) >= 2:
        groups_for_kw.append(vals.values)
        group_labels.append(ct)

if len(groups_for_kw) >= 2:
    kw_stat, kw_p = stats.kruskal(*groups_for_kw)
    print(f"\n  Kruskal-Wallis test (T3 coordination types): H={kw_stat:.3f}, p={kw_p:.4e}")
else:
    kw_stat, kw_p = np.nan, np.nan
    print(f"\n  Kruskal-Wallis: insufficient groups (need >= 2, have {len(groups_for_kw)})")

# Also check T2 coordination
print("\n  T2 coordination type distribution:")
for ct in coord_types + ['methyl_change_no_expr_change', 'ambiguous']:
    mask = master['coordination_T2'] == ct
    n = mask.sum()
    if n > 0:
        pr_vals = master.loc[mask, 'phase_ratio'].dropna()
        if len(pr_vals) > 0:
            print(f"    {ct}: n={n}, mean_PR={pr_vals.mean():.3f}")
        else:
            print(f"    {ct}: n={n}, no PR data")

# ============================================================================
# Step 5: Methylation temporal dynamics at exposed TF promoters
# ============================================================================
print("\n" + "=" * 70)
print("STEP 5: Methylation temporal dynamics at exposed TF promoters")
print("=" * 70)

# Use coordinated_df which has per-timepoint methylation counts
# It has: total_methyl_T1, total_methyl_T2, total_methyl_T3, methyl_change
coord_methyl = coordinated_df.set_index('locus_tag')[
    ['total_methyl_T1', 'total_methyl_T2', 'total_methyl_T3', 'methyl_change',
     '6mA_T1_count', '6mA_T2_count', '6mA_T3_count',
     '4mC_T1_count', '4mC_T2_count', '4mC_T3_count']
].to_dict('index')

# For each exposed TF, determine methylation temporal transitions
methyl_timing = []
for _, row in master.iterrows():
    tag = row['locus_tag']
    if tag in coord_methyl:
        m = coord_methyl[tag]
        t1 = m['total_methyl_T1']
        t2 = m['total_methyl_T2']
        t3 = m['total_methyl_T3']

        # Early transition (T1->T2)
        if t1 == 0 and t2 > 0:
            early_methyl = 'gained'
        elif t1 > 0 and t2 == 0:
            early_methyl = 'lost'
        elif t1 > 0 and t2 > 0:
            early_methyl = 'maintained'
        else:
            early_methyl = 'absent'

        # Late transition (T2->T3)
        if t2 == 0 and t3 > 0:
            late_methyl = 'gained'
        elif t2 > 0 and t3 == 0:
            late_methyl = 'lost'
        elif t2 > 0 and t3 > 0:
            late_methyl = 'maintained'
        else:
            late_methyl = 'absent'

        methyl_timing.append({
            'locus_tag': tag,
            'methyl_T1': t1, 'methyl_T2': t2, 'methyl_T3': t3,
            '6mA_T1': m['6mA_T1_count'], '6mA_T2': m['6mA_T2_count'], '6mA_T3': m['6mA_T3_count'],
            '4mC_T1': m['4mC_T1_count'], '4mC_T2': m['4mC_T2_count'], '4mC_T3': m['4mC_T3_count'],
            'methyl_change_str': m['methyl_change'],
            'early_methyl_transition': early_methyl,
            'late_methyl_transition': late_methyl,
            'phase_ratio': row['phase_ratio'],
            'temporal_class': row['temporal_class'],
            'LFC_T2vsT1': row['LFC_T2vsT1'],
            'LFC_T3vsT2': row['LFC_T3vsT2'],
            'direction': row['direction'],
            'bloc': row['bloc'],
        })

methyl_timing_df = pd.DataFrame(methyl_timing)
print(f"  Methylation timing data for {len(methyl_timing_df)} exposed TFs")

# Early methylation transition distribution
print(f"\n  Early (T1->T2) methylation transitions:")
print(methyl_timing_df['early_methyl_transition'].value_counts().to_string())
print(f"\n  Late (T2->T3) methylation transitions:")
print(methyl_timing_df['late_methyl_transition'].value_counts().to_string())

# Cross-tabulate methylation timing with expression timing
# For genes with both methylation changes and meaningful expression changes
methyl_has_change = methyl_timing_df[methyl_timing_df['early_methyl_transition'].isin(['gained', 'lost']) |
                                      methyl_timing_df['late_methyl_transition'].isin(['gained', 'lost'])]
print(f"\n  TFs with methylation changes: {len(methyl_has_change)}")

# Correlate: early expression change with early methylation change
# Encode: gained=+1, lost=-1, maintained/absent=0
def encode_methyl_change(x):
    if x == 'gained': return 1
    elif x == 'lost': return -1
    return 0

methyl_timing_df['early_methyl_code'] = methyl_timing_df['early_methyl_transition'].apply(encode_methyl_change)
methyl_timing_df['late_methyl_code'] = methyl_timing_df['late_methyl_transition'].apply(encode_methyl_change)

# Test: do early expression changes correlate with early methylation changes?
valid_methyl = methyl_timing_df.dropna(subset=['phase_ratio'])
if len(valid_methyl) > 5:
    rho_early, p_early = stats.spearmanr(valid_methyl['early_methyl_code'], valid_methyl['LFC_T2vsT1'])
    rho_late, p_late = stats.spearmanr(valid_methyl['late_methyl_code'], valid_methyl['LFC_T3vsT2'])
    print(f"\n  Spearman: early methylation code vs LFC_T2vsT1: rho={rho_early:.3f}, p={p_early:.4e}")
    print(f"  Spearman: late methylation code vs LFC_T3vsT2: rho={rho_late:.3f}, p={p_late:.4e}")
else:
    rho_early, p_early, rho_late, p_late = np.nan, np.nan, np.nan, np.nan

# Also check: AAGCCCG sites near exposed TFs
# Get TSS for each exposed TF
exposed_tss = {}
for _, row in master.iterrows():
    tag = row['locus_tag']
    if row['strand'] == '+':
        tss = row['start']
    else:
        tss = row['end']
    exposed_tss[tag] = tss

# AAGCCCG sites within 500bp of exposed TF TSS
aagcccg_near_exposed = []
for _, site_row in aagcccg_sites.iterrows():
    pos = site_row['position']
    for tag, tss in exposed_tss.items():
        if abs(pos - tss) <= 500:
            aagcccg_near_exposed.append({
                'locus_tag': tag,
                'site_position': pos,
                'distance_to_TSS': pos - tss,
                'timepoint': site_row['timepoint'],
                'site_strand': site_row['site_strand'],
            })

aagcccg_near_df = pd.DataFrame(aagcccg_near_exposed) if aagcccg_near_exposed else pd.DataFrame()
print(f"\n  AAGCCCG sites within 500bp of exposed TF TSS: {len(aagcccg_near_df)}")
if len(aagcccg_near_df) > 0:
    print(f"    Unique exposed TFs with nearby AAGCCCG: {aagcccg_near_df['locus_tag'].nunique()}")
    print(f"    By timepoint: {aagcccg_near_df['timepoint'].value_counts().to_dict()}")

# GCCGGC sites within 500bp of exposed TF TSS
gccggc_near_exposed = []
for _, site_row in gccggc_sites.iterrows():
    pos = site_row['position']
    for tag, tss in exposed_tss.items():
        if abs(pos - tss) <= 500:
            gccggc_near_exposed.append({
                'locus_tag': tag,
                'site_position': pos,
                'distance_to_TSS': pos - tss,
                'timepoint': site_row['timepoint'],
                'frequency': site_row['frequency'],
            })

gccggc_near_df = pd.DataFrame(gccggc_near_exposed) if gccggc_near_exposed else pd.DataFrame()
print(f"  GCCGGC sites within 500bp of exposed TF TSS: {len(gccggc_near_df)}")
if len(gccggc_near_df) > 0:
    print(f"    Unique exposed TFs with nearby GCCGGC: {gccggc_near_df['locus_tag'].nunique()}")
    print(f"    By timepoint: {gccggc_near_df['timepoint'].value_counts().to_dict()}")

# ============================================================================
# Step 6: TCS pair temporal ordering
# ============================================================================
print("\n" + "=" * 70)
print("STEP 6: TCS pair temporal ordering")
print("=" * 70)

tcs_temporal = []
for _, pair in tcs_pairs.iterrows():
    sk = pair['sensor_kinase']
    rr = pair['response_regulator']
    sk_exposed = pair['sk_exposed']
    rr_exposed = pair['rr_exposed']

    # Get LFC for both partners
    sk_lfc_t2t1 = deseq_t2t1_map.get(sk, np.nan)
    sk_lfc_t3t1 = deseq_t3t1_map.get(sk, np.nan)
    sk_lfc_t3t2 = deseq_t3t2_map.get(sk, np.nan)
    rr_lfc_t2t1 = deseq_t2t1_map.get(rr, np.nan)
    rr_lfc_t3t1 = deseq_t3t1_map.get(rr, np.nan)
    rr_lfc_t3t2 = deseq_t3t2_map.get(rr, np.nan)

    # Phase ratios
    sk_pr = sk_lfc_t2t1 / sk_lfc_t3t1 if abs(sk_lfc_t3t1) > LFC_THRESHOLD else np.nan
    rr_pr = rr_lfc_t2t1 / rr_lfc_t3t1 if abs(rr_lfc_t3t1) > LFC_THRESHOLD else np.nan

    # Exposed partner
    if sk_exposed:
        exposed_partner = 'SK'
        exposed_tag = sk
        shielded_tag = rr
        exposed_pr = sk_pr
        shielded_pr = rr_pr
    else:
        exposed_partner = 'RR'
        exposed_tag = rr
        shielded_tag = sk
        exposed_pr = rr_pr
        shielded_pr = sk_pr

    temporal_lag = exposed_pr - shielded_pr if not (pd.isna(exposed_pr) or pd.isna(shielded_pr)) else np.nan

    # Get normalized counts for trajectory
    sk_counts = norm_counts_idx.loc[sk] if sk in norm_counts_idx.index else pd.Series(dtype=float)
    rr_counts = norm_counts_idx.loc[rr] if rr in norm_counts_idx.index else pd.Series(dtype=float)

    tcs_temporal.append({
        'sensor_kinase': sk,
        'sk_old_locus': pair['sk_old_locus'],
        'response_regulator': rr,
        'rr_old_locus': pair['rr_old_locus'],
        'exposed_partner': exposed_partner,
        'sk_LFC_T2vsT1': sk_lfc_t2t1,
        'sk_LFC_T3vsT1': sk_lfc_t3t1,
        'sk_LFC_T3vsT2': sk_lfc_t3t2,
        'rr_LFC_T2vsT1': rr_lfc_t2t1,
        'rr_LFC_T3vsT1': rr_lfc_t3t1,
        'rr_LFC_T3vsT2': rr_lfc_t3t2,
        'sk_phase_ratio': sk_pr,
        'rr_phase_ratio': rr_pr,
        'exposed_phase_ratio': exposed_pr,
        'shielded_phase_ratio': shielded_pr,
        'temporal_lag': temporal_lag,
        'expression_rho': pair['expression_rho'],
        'sk_coordination_T3': pair['sk_coordination_T3'],
        'rr_coordination_T3': pair['rr_coordination_T3'],
    })

tcs_temporal_df = pd.DataFrame(tcs_temporal)
print(f"  TCS pairs analyzed: {len(tcs_temporal_df)}")

# Test: does the exposed SK respond earlier than the exposed RR?
sk_exposed_pairs = tcs_temporal_df[tcs_temporal_df['exposed_partner'] == 'SK']
rr_exposed_pairs = tcs_temporal_df[tcs_temporal_df['exposed_partner'] == 'RR']
print(f"  SK-exposed pairs: {len(sk_exposed_pairs)}")
print(f"  RR-exposed pairs: {len(rr_exposed_pairs)}")

# Compare exposed phase ratios between SK-exposed and RR-exposed
sk_exp_pr = sk_exposed_pairs['exposed_phase_ratio'].dropna()
rr_exp_pr = rr_exposed_pairs['exposed_phase_ratio'].dropna()
print(f"\n  SK-exposed phase ratios: {sk_exp_pr.values}")
print(f"  RR-exposed phase ratios: {rr_exp_pr.values}")
if len(sk_exp_pr) >= 2 and len(rr_exp_pr) >= 2:
    u_tcs, p_tcs = stats.mannwhitneyu(sk_exp_pr, rr_exp_pr, alternative='two-sided')
    print(f"  Mann-Whitney: U={u_tcs:.1f}, p={p_tcs:.4f}")
else:
    u_tcs, p_tcs = np.nan, np.nan
    print(f"  Insufficient data for test")

# Temporal lag: exposed vs shielded
temporal_lags = tcs_temporal_df['temporal_lag'].dropna()
print(f"\n  Temporal lag (exposed_PR - shielded_PR):")
print(f"    mean={temporal_lags.mean():.3f}, median={temporal_lags.median():.3f}")
if len(temporal_lags) >= 3:
    t_lag, p_lag = stats.wilcoxon(temporal_lags)
    print(f"    Wilcoxon signed-rank: W={t_lag:.1f}, p={p_lag:.4f}")
else:
    t_lag, p_lag = np.nan, np.nan

for _, row in tcs_temporal_df.iterrows():
    print(f"\n  Pair: {row['sk_old_locus']}/{row['rr_old_locus']} "
          f"(exposed={row['exposed_partner']})")
    print(f"    SK LFC: T2vsT1={row['sk_LFC_T2vsT1']:.3f}, T3vsT1={row['sk_LFC_T3vsT1']:.3f}, PR={row['sk_phase_ratio']:.3f}" if not pd.isna(row['sk_phase_ratio']) else f"    SK LFC: T2vsT1={row['sk_LFC_T2vsT1']:.3f}, T3vsT1={row['sk_LFC_T3vsT1']:.3f}, PR=NA")
    print(f"    RR LFC: T2vsT1={row['rr_LFC_T2vsT1']:.3f}, T3vsT1={row['rr_LFC_T3vsT1']:.3f}, PR={row['rr_phase_ratio']:.3f}" if not pd.isna(row['rr_phase_ratio']) else f"    RR LFC: T2vsT1={row['rr_LFC_T2vsT1']:.3f}, T3vsT1={row['rr_LFC_T3vsT1']:.3f}, PR=NA")

# ============================================================================
# Step 7: Hierarchical structure detection
# ============================================================================
print("\n" + "=" * 70)
print("STEP 7: Hierarchical structure detection")
print("=" * 70)

# Cluster TFs by temporal trajectories (z-scored T1, T2, T3)
z_matrix = master[['T1_z', 'T2_z', 'T3_z']].values
# Use Euclidean distance for clustering
dist_mat = pdist(z_matrix, metric='euclidean')
linkage_mat = linkage(dist_mat, method='ward')

# Cut into 4 clusters (to match module structure)
master['trajectory_cluster'] = fcluster(linkage_mat, t=4, criterion='maxclust')

print(f"  Trajectory clusters (Ward linkage, k=4):")
for cl in sorted(master['trajectory_cluster'].unique()):
    subset = master[master['trajectory_cluster'] == cl]
    print(f"    Cluster {cl}: n={len(subset)}, "
          f"mean T1_z={subset['T1_z'].mean():.2f}, T2_z={subset['T2_z'].mean():.2f}, T3_z={subset['T3_z'].mean():.2f}")

# Identify pioneers (earliest responders) and followers (latest responders)
# Use the absolute change at T2 relative to total change as earliness metric
master['earliness_score'] = master['early_magnitude'] / (master['early_magnitude'] + master['late_magnitude'] + 1e-10)

# Top 10 pioneers (highest earliness) and top 10 followers (lowest earliness)
pioneers = master.nlargest(10, 'earliness_score')
followers = master.nsmallest(10, 'earliness_score')

print(f"\n  Top 10 PIONEER TFs (earliest responders):")
for _, row in pioneers.iterrows():
    old_lt = str(row['old_locus_tag']) if not pd.isna(row.get('old_locus_tag', float('nan'))) else row['locus_tag']
    pr_str = f"{row['phase_ratio']:.3f}" if not pd.isna(row['phase_ratio']) else "NA"
    print(f"    {old_lt:10s} ({str(row['tf_family']):20s}) "
          f"earliness={row['earliness_score']:.3f}, PR={pr_str}")

print(f"\n  Top 10 FOLLOWER TFs (latest responders):")
for _, row in followers.iterrows():
    old_lt = str(row['old_locus_tag']) if not pd.isna(row.get('old_locus_tag', float('nan'))) else row['locus_tag']
    pr_str = f"{row['phase_ratio']:.3f}" if not pd.isna(row['phase_ratio']) else "NA"
    print(f"    {old_lt:10s} ({str(row['tf_family']):20s}) "
          f"earliness={row['earliness_score']:.3f}, PR={pr_str}")

# Are pioneer TFs co-expressed with follower TFs?
pioneer_tags = set(pioneers['locus_tag'])
follower_tags = set(followers['locus_tag'])

# Check if pioneers and followers are in the same or different modules
pioneer_modules = pioneers['module'].value_counts()
follower_modules = followers['module'].value_counts()
print(f"\n  Pioneer module distribution: {pioneer_modules.to_dict()}")
print(f"  Follower module distribution: {follower_modules.to_dict()}")

# Granger-like analysis: do early-responding TFs' T2 levels predict late-responding TFs' T3 changes?
# Cross-correlation with 1-timepoint lag
early_responders = master[master['temporal_class'] == 'early']
late_responders = master[master['temporal_class'] == 'late']

print(f"\n  Granger-like analysis:")
print(f"    Early responders: {len(early_responders)}")
print(f"    Late responders: {len(late_responders)}")

if len(early_responders) >= 3 and len(late_responders) >= 3:
    # Mean T2 expression of early responders
    early_mean_T2 = early_responders['T2_z'].mean()
    # Correlate individual late responders' T3 change (T3_z - T2_z) with mean of early responders' T2 level
    # This is a cross-sectional proxy since we have only 3 timepoints
    late_responders_t3_change = late_responders['T3_z'] - late_responders['T2_z']
    print(f"    Early responders mean T2_z: {early_mean_T2:.3f}")
    print(f"    Late responders T3 change (T3_z - T2_z): mean={late_responders_t3_change.mean():.3f}")

    # Better approach: for each late responder, correlate its T3 change with the
    # ensemble mean T2 expression of early responders
    # Since this is n=1 predictor vs n outcomes, compute across replicate variation

    # Alternative: Rank-based comparison - do late responders whose cognate early-responders
    # had higher T2 expression show more T3 change?
    # We can correlate each late responder's LFC_T3vsT2 with the mean LFC_T2vsT1 of
    # early responders in the same module
    print(f"    Note: with only 3 timepoints, true Granger causality cannot be assessed.")
    print(f"    Using cross-correlation proxy instead.")

    # Cross-correlation: for each pair (early_i, late_j),
    # correlate early_i's T2 response with late_j's T3 response
    pairs_data = []
    for _, e_row in early_responders.iterrows():
        for _, l_row in late_responders.iterrows():
            pairs_data.append({
                'early_tag': e_row['locus_tag'],
                'late_tag': l_row['locus_tag'],
                'early_T2_response': e_row['LFC_T2vsT1'],
                'late_T3_response': l_row['LFC_T3vsT2'],
                'same_module': e_row['module'] == l_row['module'],
            })
    pairs_df = pd.DataFrame(pairs_data)
    if len(pairs_df) > 0:
        rho_cross, p_cross = stats.spearmanr(pairs_df['early_T2_response'], pairs_df['late_T3_response'])
        print(f"    Cross-correlation (all pairs): rho={rho_cross:.3f}, p={p_cross:.4e}, n_pairs={len(pairs_df)}")

        # Same module only
        same_mod = pairs_df[pairs_df['same_module']]
        if len(same_mod) > 5:
            rho_sm, p_sm = stats.spearmanr(same_mod['early_T2_response'], same_mod['late_T3_response'])
            print(f"    Cross-correlation (same module): rho={rho_sm:.3f}, p={p_sm:.4e}, n_pairs={len(same_mod)}")
        else:
            rho_sm, p_sm = np.nan, np.nan
else:
    rho_cross, p_cross = np.nan, np.nan
    rho_sm, p_sm = np.nan, np.nan
    print("    Insufficient data for cross-correlation analysis")

# ============================================================================
# Step 8: Comparison with shielded regulatory genes
# ============================================================================
print("\n" + "=" * 70)
print("STEP 8: Comparison with shielded regulatory genes")
print("=" * 70)

# Build temporal metrics for all 1,017 regulatory genes
shielded_df = all_genes_df[all_genes_df['is_exposed'] == 0].copy()
exposed_all = all_genes_df[all_genes_df['is_exposed'] == 1].copy()
print(f"  Shielded regulatory genes: {len(shielded_df)}")
print(f"  Exposed regulatory genes (from all_genes): {len(exposed_all)}")

# Add DESeq2 LFC and padj
for df_subset in [shielded_df, exposed_all]:
    df_subset['LFC_T2vsT1_deseq'] = df_subset['locus_tag'].map(deseq_t2t1_map)
    df_subset['LFC_T3vsT1_deseq'] = df_subset['locus_tag'].map(deseq_t3t1_map)
    df_subset['LFC_T3vsT2_deseq'] = df_subset['locus_tag'].map(deseq_t3t2_map)

# Phase ratio for shielded
def compute_phase_ratio_series(df):
    pr = []
    for _, row in df.iterrows():
        lfc_t2t1 = row.get('LFC_T2vsT1_deseq', row.get('LFC_T2vsT1', np.nan))
        lfc_t3t1 = row.get('LFC_T3vsT1_deseq', row.get('LFC_T3vsT1', np.nan))
        if pd.isna(lfc_t2t1) or pd.isna(lfc_t3t1) or abs(lfc_t3t1) <= LFC_THRESHOLD:
            pr.append(np.nan)
        else:
            pr.append(lfc_t2t1 / lfc_t3t1)
    return pr

shielded_df = shielded_df.copy()
shielded_df['phase_ratio'] = compute_phase_ratio_series(shielded_df)
shielded_df['temporal_class'] = shielded_df['phase_ratio'].apply(
    lambda pr: 'early' if pr > 0.6 else ('late' if pr < 0.4 else ('gradual' if not pd.isna(pr) else 'non-responder'))
)
shielded_df['has_meaningful_change'] = shielded_df.apply(
    lambda row: abs(row.get('LFC_T3vsT1_deseq', row.get('LFC_T3vsT1', 0))) > LFC_THRESHOLD, axis=1
)

# exposed from master
exposed_pr = master['phase_ratio'].dropna()
shielded_pr = shielded_df['phase_ratio'].dropna()

print(f"\n  Exposed phase ratio: mean={exposed_pr.mean():.3f}, median={exposed_pr.median():.3f}, n={len(exposed_pr)}")
print(f"  Shielded phase ratio: mean={shielded_pr.mean():.3f}, median={shielded_pr.median():.3f}, n={len(shielded_pr)}")

# Mann-Whitney test
if len(exposed_pr) > 0 and len(shielded_pr) > 0:
    u_es, p_es = stats.mannwhitneyu(exposed_pr, shielded_pr, alternative='two-sided')
    n1, n2 = len(exposed_pr), len(shielded_pr)
    r_es = 1 - (2 * u_es) / (n1 * n2)
    print(f"  Mann-Whitney U: U={u_es:.1f}, p={p_es:.4e}, r={r_es:.3f}")
else:
    u_es, p_es, r_es = np.nan, np.nan, np.nan

# Temporal class distribution comparison
exp_class = master['temporal_class'].value_counts()
shi_class = shielded_df['temporal_class'].value_counts()
print(f"\n  Temporal class distribution:")
print(f"  {'Class':<15} {'Exposed':>10} {'Shielded':>10} {'Exp%':>10} {'Shi%':>10}")
for cls in ['early', 'late', 'gradual', 'non-responder']:
    e_n = exp_class.get(cls, 0)
    s_n = shi_class.get(cls, 0)
    e_pct = e_n / len(master) * 100
    s_pct = s_n / len(shielded_df) * 100
    print(f"  {cls:<15} {e_n:>10} {s_n:>10} {e_pct:>9.1f}% {s_pct:>9.1f}%")

# Chi-square test for temporal class distribution (exposed vs shielded)
# Combine into contingency table (early/late/gradual only - exclude non-responders)
classes = ['early', 'late', 'gradual']
exp_counts = [exp_class.get(c, 0) for c in classes]
shi_counts = [shi_class.get(c, 0) for c in classes]
contingency = np.array([exp_counts, shi_counts])
if contingency.sum() > 0:
    chi2, p_chi2, dof, expected = stats.chi2_contingency(contingency)
    print(f"\n  Chi-square test (temporal class distribution): chi2={chi2:.3f}, p={p_chi2:.4e}, dof={dof}")
else:
    chi2, p_chi2 = np.nan, np.nan

# Are exposed TFs biased toward early or late?
exp_early_frac = exp_class.get('early', 0) / (exp_class.get('early', 0) + exp_class.get('late', 0) + exp_class.get('gradual', 0) + 1e-10)
shi_early_frac = shi_class.get('early', 0) / (shi_class.get('early', 0) + shi_class.get('late', 0) + shi_class.get('gradual', 0) + 1e-10)
print(f"\n  Early fraction: exposed={exp_early_frac:.3f}, shielded={shi_early_frac:.3f}")

# ============================================================================
# Collect all statistical tests
# ============================================================================
stat_tests = []
stat_tests.append({'test_name': 'Bloc phase ratio comparison (Mann-Whitney U)',
                    'group1': 'Activation (modules 1-3)', 'group2': 'Repression (module 4)',
                    'statistic': stat_w, 'p_value': p_w, 'effect_size': r_rbs,
                    'n1': len(act_pr), 'n2': len(rep_pr), 'interpretation': 'Activation vs repression bloc temporal ordering'})
stat_tests.append({'test_name': 'Coordination type phase ratio (Kruskal-Wallis)',
                    'group1': 'Coordination types (T3)', 'group2': 'N/A',
                    'statistic': kw_stat, 'p_value': kw_p, 'effect_size': np.nan,
                    'n1': sum(len(g) for g in groups_for_kw), 'n2': len(groups_for_kw),
                    'interpretation': 'Differences in temporal response across coordination types'})
stat_tests.append({'test_name': 'Early methylation-expression correlation (Spearman)',
                    'group1': 'Early methyl code', 'group2': 'LFC_T2vsT1',
                    'statistic': rho_early, 'p_value': p_early, 'effect_size': rho_early,
                    'n1': len(valid_methyl), 'n2': len(valid_methyl),
                    'interpretation': 'Do early methylation changes correlate with early expression changes?'})
stat_tests.append({'test_name': 'Late methylation-expression correlation (Spearman)',
                    'group1': 'Late methyl code', 'group2': 'LFC_T3vsT2',
                    'statistic': rho_late, 'p_value': p_late, 'effect_size': rho_late,
                    'n1': len(valid_methyl), 'n2': len(valid_methyl),
                    'interpretation': 'Do late methylation changes correlate with late expression changes?'})
stat_tests.append({'test_name': 'TCS temporal lag (Wilcoxon signed-rank)',
                    'group1': 'Exposed PR', 'group2': 'Shielded PR',
                    'statistic': t_lag, 'p_value': p_lag, 'effect_size': np.nan,
                    'n1': len(temporal_lags), 'n2': len(temporal_lags),
                    'interpretation': 'Does exposed partner respond before/after shielded?'})
stat_tests.append({'test_name': 'SK-exposed vs RR-exposed phase ratio (Mann-Whitney)',
                    'group1': 'SK-exposed pairs', 'group2': 'RR-exposed pairs',
                    'statistic': u_tcs, 'p_value': p_tcs, 'effect_size': np.nan,
                    'n1': len(sk_exp_pr), 'n2': len(rr_exp_pr),
                    'interpretation': 'Do SK-exposed pairs respond earlier than RR-exposed?'})
stat_tests.append({'test_name': 'Exposed vs shielded phase ratio (Mann-Whitney)',
                    'group1': 'Exposed TFs', 'group2': 'Shielded regulatory genes',
                    'statistic': u_es, 'p_value': p_es, 'effect_size': r_es,
                    'n1': len(exposed_pr), 'n2': len(shielded_pr),
                    'interpretation': 'Do exposed TFs show different temporal distribution than shielded?'})
stat_tests.append({'test_name': 'Temporal class distribution (Chi-square)',
                    'group1': 'Exposed', 'group2': 'Shielded',
                    'statistic': chi2, 'p_value': p_chi2, 'effect_size': np.nan,
                    'n1': sum(exp_counts), 'n2': sum(shi_counts),
                    'interpretation': 'Is temporal class distribution different for exposed vs shielded?'})

# Granger cross-correlation
if not pd.isna(rho_cross):
    stat_tests.append({'test_name': 'Cross-correlation early T2 vs late T3 (Spearman)',
                        'group1': 'Early responders T2', 'group2': 'Late responders T3',
                        'statistic': rho_cross, 'p_value': p_cross, 'effect_size': rho_cross,
                        'n1': len(early_responders), 'n2': len(late_responders),
                        'interpretation': 'Granger-like: do early responders T2 predict late responders T3?'})

stat_tests_df = pd.DataFrame(stat_tests)

# ============================================================================
# Step 9: Figures
# ============================================================================
print("\n" + "=" * 70)
print("STEP 9: Generating figures")
print("=" * 70)

# Figure 1: Temporal classification scatter
fig, ax = plt.subplots(figsize=(10, 8))
for tc in ['early', 'late', 'gradual', 'non-responder']:
    for dirn in ['up', 'down']:
        mask = (master['temporal_class'] == tc) & (master['direction'] == dirn)
        subset = master[mask]
        marker = '^' if dirn == 'up' else 'v'
        color = COLORS.get(tc, '#999999')
        alpha = 0.3 if tc == 'non-responder' else 0.8
        ax.scatter(subset['LFC_T2vsT1'], subset['LFC_T3vsT2'],
                  c=color, marker=marker, s=80, alpha=alpha, edgecolors='k', linewidths=0.5,
                  label=f'{tc} ({dirn})')
        # Add labels for notable genes
        for _, row in subset.iterrows():
            if abs(row['LFC_T2vsT1']) > 1.5 or abs(row['LFC_T3vsT2']) > 1.5:
                label = row['old_locus_tag'] if pd.notna(row['old_locus_tag']) and row['old_locus_tag'] != '' else row['locus_tag']
                ax.annotate(label, (row['LFC_T2vsT1'], row['LFC_T3vsT2']),
                           fontsize=7, ha='left', va='bottom', alpha=0.7)

# Add quadrant lines and phase ratio boundaries
ax.axhline(0, color='gray', lw=0.5, ls='--')
ax.axvline(0, color='gray', lw=0.5, ls='--')
ax.set_xlabel('LFC (T2 vs T1) — Early response')
ax.set_ylabel('LFC (T3 vs T2) — Late response')
ax.set_title('H34: Temporal Classification of 62 Exposed TFs')
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f'{FIG_DIR}/temporal_classification.{ext}')
plt.close()
print("  [1/7] temporal_classification.pdf/svg")

# Figure 2: Bloc temporal trajectories
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
timepoints = [1, 2, 3]
tp_labels = ['T1\n(Exponential)', 'T2\n(Transition)', 'T3\n(Stationary)']

# Panel A: Individual gene traces + mean
ax = axes[0]
for _, row in activation.iterrows():
    ax.plot(timepoints, [row['T1_z'], row['T2_z'], row['T3_z']],
            color=COLORS['activation'], alpha=0.15, lw=0.8)
for _, row in repression.iterrows():
    ax.plot(timepoints, [row['T1_z'], row['T2_z'], row['T3_z']],
            color=COLORS['repression'], alpha=0.15, lw=0.8)

# Mean trajectories
act_means = [activation['T1_z'].mean(), activation['T2_z'].mean(), activation['T3_z'].mean()]
act_stds = [activation['T1_z'].std(), activation['T2_z'].std(), activation['T3_z'].std()]
rep_means = [repression['T1_z'].mean(), repression['T2_z'].mean(), repression['T3_z'].mean()]
rep_stds = [repression['T1_z'].std(), repression['T2_z'].std(), repression['T3_z'].std()]

ax.plot(timepoints, act_means, 'o-', color=COLORS['activation'], lw=3, ms=10,
        label=f'Activation (n={len(activation)})', zorder=5)
ax.fill_between(timepoints,
                [m-s for m, s in zip(act_means, act_stds)],
                [m+s for m, s in zip(act_means, act_stds)],
                color=COLORS['activation'], alpha=0.2)
ax.plot(timepoints, rep_means, 's-', color=COLORS['repression'], lw=3, ms=10,
        label=f'Repression (n={len(repression)})', zorder=5)
ax.fill_between(timepoints,
                [m-s for m, s in zip(rep_means, rep_stds)],
                [m+s for m, s in zip(rep_means, rep_stds)],
                color=COLORS['repression'], alpha=0.2)

ax.set_xticks(timepoints)
ax.set_xticklabels(tp_labels)
ax.set_ylabel('Expression (z-score)')
ax.set_title('A. Bloc temporal trajectories')
ax.legend()
ax.axhline(0, color='gray', lw=0.5, ls='--')

# Panel B: Module-level
ax = axes[1]
module_colors = {1: COLORS['module1'], 2: COLORS['module2'],
                 3: COLORS['module3'], 4: COLORS['module4']}
for mod_id in [1, 2, 3, 4]:
    mod_genes = master[master['module'] == mod_id]
    means = [mod_genes['T1_z'].mean(), mod_genes['T2_z'].mean(), mod_genes['T3_z'].mean()]
    stds = [mod_genes['T1_z'].std(), mod_genes['T2_z'].std(), mod_genes['T3_z'].std()]
    ax.plot(timepoints, means, 'o-', color=module_colors[mod_id], lw=2.5, ms=8,
            label=f'Module {mod_id} (n={len(mod_genes)})')
    ax.fill_between(timepoints,
                    [m-s for m, s in zip(means, stds)],
                    [m+s for m, s in zip(means, stds)],
                    color=module_colors[mod_id], alpha=0.15)

ax.set_xticks(timepoints)
ax.set_xticklabels(tp_labels)
ax.set_ylabel('Expression (z-score)')
ax.set_title('B. Module-level temporal trajectories')
ax.legend()
ax.axhline(0, color='gray', lw=0.5, ls='--')

plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f'{FIG_DIR}/bloc_temporal_trajectories.{ext}')
plt.close()
print("  [2/7] bloc_temporal_trajectories.pdf/svg")

# Figure 3: Phase ratio distribution
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Panel A: Exposed vs shielded
ax = axes[0]
bins = np.linspace(-0.5, 2.0, 26)
ax.hist(exposed_pr, bins=bins, color=COLORS['exposed'], alpha=0.7,
        label=f'Exposed (n={len(exposed_pr)})', density=True, edgecolor='k', linewidth=0.5)
ax.hist(shielded_pr, bins=bins, color=COLORS['shielded'], alpha=0.5,
        label=f'Shielded (n={len(shielded_pr)})', density=True, edgecolor='k', linewidth=0.5)
ax.axvline(0.4, color='gray', ls='--', lw=1, label='Late/Gradual boundary')
ax.axvline(0.6, color='gray', ls=':', lw=1, label='Gradual/Early boundary')
ax.set_xlabel('Phase ratio (LFC_T2vsT1 / LFC_T3vsT1)')
ax.set_ylabel('Density')
ax.set_title(f'A. Phase ratio: Exposed vs Shielded\n(MW p={p_es:.3e})')
ax.legend(fontsize=8)

# Panel B: Activation vs repression
ax = axes[1]
ax.hist(act_pr, bins=bins, color=COLORS['activation'], alpha=0.7,
        label=f'Activation (n={len(act_pr)})', density=True, edgecolor='k', linewidth=0.5)
ax.hist(rep_pr, bins=bins, color=COLORS['repression'], alpha=0.5,
        label=f'Repression (n={len(rep_pr)})', density=True, edgecolor='k', linewidth=0.5)
ax.axvline(0.4, color='gray', ls='--', lw=1)
ax.axvline(0.6, color='gray', ls=':', lw=1)
ax.set_xlabel('Phase ratio (LFC_T2vsT1 / LFC_T3vsT1)')
ax.set_ylabel('Density')
ax.set_title(f'B. Phase ratio: Activation vs Repression\n(MW p={p_w:.3e})')
ax.legend(fontsize=8)

plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f'{FIG_DIR}/phase_ratio_distribution.{ext}')
plt.close()
print("  [3/7] phase_ratio_distribution.pdf/svg")

# Figure 4: Coordination type temporal (box plots)
fig, ax = plt.subplots(figsize=(10, 6))
plot_data = []
plot_labels = []
plot_colors = ['#e74c3c', '#2ecc71', '#3498db', '#f39c12', '#9b59b6', '#95a5a6']

# All coordination types (T3) that appear in exposed TFs
coord_types_all = master['coordination_T3'].value_counts()
ct_order = coord_types_all.index.tolist()

bp_data = []
bp_labels = []
for ct in ct_order:
    vals = master.loc[master['coordination_T3'] == ct, 'phase_ratio'].dropna()
    if len(vals) >= 1:
        bp_data.append(vals.values)
        bp_labels.append(f'{ct}\n(n={len(vals)})')

if len(bp_data) >= 1:
    bp = ax.boxplot(bp_data, labels=bp_labels, patch_artist=True, widths=0.6)
    for i, patch in enumerate(bp['boxes']):
        patch.set_facecolor(plot_colors[i % len(plot_colors)])
        patch.set_alpha(0.7)
else:
    ax.text(0.5, 0.5, 'Insufficient data', ha='center', va='center',
            transform=ax.transAxes)

ax.axhline(0.4, color='gray', ls='--', lw=1, alpha=0.5)
ax.axhline(0.6, color='gray', ls=':', lw=1, alpha=0.5)
ax.set_ylabel('Phase ratio')
ax.set_title(f'H34: Phase ratio by coordination type (T3)\n(Kruskal-Wallis p={kw_p:.3e})')
plt.xticks(rotation=30, ha='right')
plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f'{FIG_DIR}/coordination_type_temporal.{ext}')
plt.close()
print("  [4/7] coordination_type_temporal.pdf/svg")

# Figure 5: TCS temporal ordering
fig, axes = plt.subplots(2, 4, figsize=(18, 8))
axes = axes.flatten()

for i, (_, pair) in enumerate(tcs_temporal_df.iterrows()):
    if i >= 7:
        break
    ax = axes[i]
    sk = pair['sensor_kinase']
    rr = pair['response_regulator']

    # Get normalized counts
    if sk in norm_counts_idx.index:
        sk_vals = [norm_counts_idx.loc[sk, T1_COLS].mean(),
                   norm_counts_idx.loc[sk, T2_COLS].mean(),
                   norm_counts_idx.loc[sk, T3_COLS].mean()]
    else:
        sk_vals = [np.nan, np.nan, np.nan]

    if rr in norm_counts_idx.index:
        rr_vals = [norm_counts_idx.loc[rr, T1_COLS].mean(),
                   norm_counts_idx.loc[rr, T2_COLS].mean(),
                   norm_counts_idx.loc[rr, T3_COLS].mean()]
    else:
        rr_vals = [np.nan, np.nan, np.nan]

    # Z-score
    for vals in [sk_vals, rr_vals]:
        mu = np.nanmean(vals)
        sd = np.nanstd(vals)
        if sd > 0:
            for j in range(3):
                vals[j] = (vals[j] - mu) / sd
        else:
            for j in range(3):
                vals[j] = 0

    exposed_color = COLORS['exposed']
    shielded_color = COLORS['shielded']

    if pair['exposed_partner'] == 'SK':
        ax.plot(timepoints, sk_vals, 'o-', color=exposed_color, lw=2.5, ms=8,
                label=f'SK {pair["sk_old_locus"]} (exposed)')
        ax.plot(timepoints, rr_vals, 's--', color=shielded_color, lw=2, ms=8,
                label=f'RR {pair["rr_old_locus"]} (shielded)')
    else:
        ax.plot(timepoints, sk_vals, 'o--', color=shielded_color, lw=2, ms=8,
                label=f'SK {pair["sk_old_locus"]} (shielded)')
        ax.plot(timepoints, rr_vals, 's-', color=exposed_color, lw=2.5, ms=8,
                label=f'RR {pair["rr_old_locus"]} (exposed)')

    ax.set_xticks(timepoints)
    ax.set_xticklabels(['T1', 'T2', 'T3'])
    ax.set_ylabel('Expression (z)')
    lag_str = f'lag={pair["temporal_lag"]:.2f}' if not pd.isna(pair['temporal_lag']) else 'lag=NA'
    ax.set_title(f'{pair["sk_old_locus"]}/{pair["rr_old_locus"]}\n({lag_str})', fontsize=9)
    ax.legend(fontsize=7)
    ax.axhline(0, color='gray', lw=0.5, ls='--')

# Hide unused subplot
axes[7].axis('off')

fig.suptitle('H34: TCS Pair Temporal Ordering (Exposed vs Shielded Partner)', fontsize=13, y=1.02)
plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f'{FIG_DIR}/TCS_temporal_ordering.{ext}')
plt.close()
print("  [5/7] TCS_temporal_ordering.pdf/svg")

# Figure 6: Methylation-expression timing
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Panel A: Early methylation transition vs early expression change
ax = axes[0]
jitter = np.random.normal(0, 0.05, len(methyl_timing_df))
categories = ['absent', 'lost', 'maintained', 'gained']
cat_pos = {c: i for i, c in enumerate(categories)}
methyl_timing_df['early_methyl_pos'] = methyl_timing_df['early_methyl_transition'].map(cat_pos) + jitter[:len(methyl_timing_df)]

for cat in categories:
    mask = methyl_timing_df['early_methyl_transition'] == cat
    subset = methyl_timing_df[mask]
    ax.scatter(subset['early_methyl_pos'], subset['LFC_T2vsT1'],
              c=COLORS.get(cat, '#999999'), s=50, alpha=0.7, edgecolors='k', linewidths=0.5,
              label=f'{cat} (n={len(subset)})')
ax.set_xticks(range(len(categories)))
ax.set_xticklabels(categories)
ax.set_xlabel('Methylation transition (T1→T2)')
ax.set_ylabel('LFC (T2 vs T1)')
ax.set_title(f'A. Early methylation vs expression\n(rho={rho_early:.3f}, p={p_early:.3e})')
ax.axhline(0, color='gray', lw=0.5, ls='--')
ax.legend(fontsize=8)

# Panel B: Late methylation transition vs late expression change
ax = axes[1]
jitter2 = np.random.normal(0, 0.05, len(methyl_timing_df))
methyl_timing_df['late_methyl_pos'] = methyl_timing_df['late_methyl_transition'].map(cat_pos) + jitter2[:len(methyl_timing_df)]

for cat in categories:
    mask = methyl_timing_df['late_methyl_transition'] == cat
    subset = methyl_timing_df[mask]
    ax.scatter(subset['late_methyl_pos'], subset['LFC_T3vsT2'],
              c=COLORS.get(cat, '#999999'), s=50, alpha=0.7, edgecolors='k', linewidths=0.5,
              label=f'{cat} (n={len(subset)})')
ax.set_xticks(range(len(categories)))
ax.set_xticklabels(categories)
ax.set_xlabel('Methylation transition (T2→T3)')
ax.set_ylabel('LFC (T3 vs T2)')
ax.set_title(f'B. Late methylation vs expression\n(rho={rho_late:.3f}, p={p_late:.3e})')
ax.axhline(0, color='gray', lw=0.5, ls='--')
ax.legend(fontsize=8)

# Panel C: Methylation total change timeline
ax = axes[2]
# Show total methylation per timepoint for exposed TFs
methyl_totals = methyl_timing_df[['methyl_T1', 'methyl_T2', 'methyl_T3']].sum()
ax.bar([1, 2, 3], [methyl_totals['methyl_T1'], methyl_totals['methyl_T2'], methyl_totals['methyl_T3']],
       color=['#3498db', '#e74c3c', '#2ecc71'], alpha=0.7, edgecolor='k')
ax.set_xticks([1, 2, 3])
ax.set_xticklabels(['T1', 'T2', 'T3'])
ax.set_ylabel('Total methylation sites at exposed TF promoters')
ax.set_title('C. Methylation site dynamics')

# Add 6mA / 4mC breakdown
for i, tp in enumerate(['T1', 'T2', 'T3']):
    m6a = methyl_timing_df[f'6mA_{tp}'].sum()
    m4c = methyl_timing_df[f'4mC_{tp}'].sum()
    ax.text(i+1, methyl_totals.iloc[i] + 0.5, f'6mA:{int(m6a)}\n4mC:{int(m4c)}',
           ha='center', fontsize=8)

plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f'{FIG_DIR}/methylation_expression_timing.{ext}')
plt.close()
print("  [6/7] methylation_expression_timing.pdf/svg")

# Figure 7: Comprehensive summary (6 panels)
fig = plt.figure(figsize=(20, 14))
gs = gridspec.GridSpec(3, 3, hspace=0.4, wspace=0.35)

# Panel A: Temporal classification scatter (compact)
ax = fig.add_subplot(gs[0, 0])
for tc in ['early', 'late', 'gradual']:
    mask = master['temporal_class'] == tc
    subset = master[mask]
    ax.scatter(subset['LFC_T2vsT1'], subset['LFC_T3vsT2'],
              c=COLORS.get(tc, '#999999'), s=40, alpha=0.7, edgecolors='k', linewidths=0.3,
              label=f'{tc} (n={len(subset)})')
ax.axhline(0, color='gray', lw=0.5, ls='--')
ax.axvline(0, color='gray', lw=0.5, ls='--')
ax.set_xlabel('LFC (T2 vs T1)')
ax.set_ylabel('LFC (T3 vs T2)')
ax.set_title('A. Temporal classification')
ax.legend(fontsize=7)

# Panel B: Bloc trajectories
ax = fig.add_subplot(gs[0, 1])
ax.plot(timepoints, act_means, 'o-', color=COLORS['activation'], lw=2.5, ms=8,
        label=f'Activation (n={len(activation)})')
ax.fill_between(timepoints,
                [m-s for m, s in zip(act_means, act_stds)],
                [m+s for m, s in zip(act_means, act_stds)],
                color=COLORS['activation'], alpha=0.2)
ax.plot(timepoints, rep_means, 's-', color=COLORS['repression'], lw=2.5, ms=8,
        label=f'Repression (n={len(repression)})')
ax.fill_between(timepoints,
                [m-s for m, s in zip(rep_means, rep_stds)],
                [m+s for m, s in zip(rep_means, rep_stds)],
                color=COLORS['repression'], alpha=0.2)
ax.set_xticks(timepoints)
ax.set_xticklabels(['T1', 'T2', 'T3'])
ax.set_ylabel('Expression (z-score)')
ax.set_title(f'B. Bloc trajectories (p={p_w:.3e})')
ax.legend(fontsize=7)
ax.axhline(0, color='gray', lw=0.5, ls='--')

# Panel C: Phase ratio exposed vs shielded
ax = fig.add_subplot(gs[0, 2])
bp = ax.boxplot([exposed_pr.values, shielded_pr.values],
                labels=[f'Exposed\n(n={len(exposed_pr)})', f'Shielded\n(n={len(shielded_pr)})'],
                patch_artist=True, widths=0.5)
bp['boxes'][0].set_facecolor(COLORS['exposed'])
bp['boxes'][0].set_alpha(0.7)
bp['boxes'][1].set_facecolor(COLORS['shielded'])
bp['boxes'][1].set_alpha(0.7)
ax.set_ylabel('Phase ratio')
ax.set_title(f'C. Exposed vs Shielded (p={p_es:.3e})')
ax.axhline(0.4, color='gray', ls='--', lw=0.8, alpha=0.5)
ax.axhline(0.6, color='gray', ls=':', lw=0.8, alpha=0.5)

# Panel D: Coordination type temporal
ax = fig.add_subplot(gs[1, 0:2])
if len(bp_data) >= 1:
    bp2 = ax.boxplot(bp_data, labels=[l.replace('\n', ' ') for l in bp_labels],
                     patch_artist=True, widths=0.5)
    for i, patch in enumerate(bp2['boxes']):
        patch.set_facecolor(plot_colors[i % len(plot_colors)])
        patch.set_alpha(0.7)
else:
    ax.text(0.5, 0.5, 'Insufficient data', ha='center', va='center',
            transform=ax.transAxes)
ax.axhline(0.4, color='gray', ls='--', lw=0.8, alpha=0.5)
ax.axhline(0.6, color='gray', ls=':', lw=0.8, alpha=0.5)
ax.set_ylabel('Phase ratio')
ax.set_title(f'D. Phase ratio by coordination type (KW p={kw_p:.3e})')
plt.setp(ax.get_xticklabels(), rotation=20, ha='right', fontsize=8)

# Panel E: Hierarchical clustering dendrogram
ax = fig.add_subplot(gs[1, 2])
# Color by module
module_colors_map = {1: COLORS['module1'], 2: COLORS['module2'],
                     3: COLORS['module3'], 4: COLORS['module4']}
leaf_colors = {}
for i, (_, row) in enumerate(master.iterrows()):
    mod = row['module']
    leaf_colors[i] = module_colors_map.get(mod, '#999999')

from scipy.cluster.hierarchy import set_link_color_palette
set_link_color_palette(['#333333'])
dend = dendrogram(linkage_mat, ax=ax, labels=master['old_locus_tag'].values,
                  leaf_rotation=90, leaf_font_size=5, color_threshold=0,
                  above_threshold_color='#333333')
ax.set_title('E. Temporal trajectory clustering')
ax.set_ylabel('Ward distance')

# Panel F: Pioneer vs follower table
ax = fig.add_subplot(gs[2, :])
ax.axis('off')

# Create summary text table
summary_text = "F. SUMMARY: Temporal Dynamics of 62 Exposed TFs\n\n"
summary_text += f"Phase ratio distribution: mean={exposed_pr.mean():.3f}, median={exposed_pr.median():.3f}\n"
summary_text += f"  Early responders: {exp_class.get('early', 0)} ({exp_class.get('early', 0)/len(master)*100:.1f}%)\n"
summary_text += f"  Late responders: {exp_class.get('late', 0)} ({exp_class.get('late', 0)/len(master)*100:.1f}%)\n"
summary_text += f"  Gradual: {exp_class.get('gradual', 0)} ({exp_class.get('gradual', 0)/len(master)*100:.1f}%)\n"
summary_text += f"  Non-responder: {exp_class.get('non-responder', 0)} ({exp_class.get('non-responder', 0)/len(master)*100:.1f}%)\n\n"
summary_text += f"Activation bloc (mod 1-3, n={len(activation)}): mean PR={act_pr.mean():.3f}\n"
summary_text += f"Repression bloc (mod 4, n={len(repression)}): mean PR={rep_pr.mean():.3f}\n"
summary_text += f"Bloc comparison: MW p={p_w:.3e}, r={r_rbs:.3f}\n\n"
summary_text += f"Exposed vs Shielded: MW p={p_es:.3e}, r={r_es:.3f}\n"
summary_text += f"  Exposed early fraction: {exp_early_frac:.3f}, Shielded: {shi_early_frac:.3f}\n"

ax.text(0.02, 0.95, summary_text, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

fig.suptitle('H34: Temporal Dynamics and Hierarchical Structure of 62 Exposed TFs',
             fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f'{FIG_DIR}/H34_comprehensive_summary.{ext}')
plt.close()
print("  [7/7] H34_comprehensive_summary.pdf/svg")

# ============================================================================
# Step 10: Save tables
# ============================================================================
print("\n" + "=" * 70)
print("STEP 10: Saving tables")
print("=" * 70)

# Table 1: Temporal classification
table1 = master[['locus_tag', 'old_locus_tag', 'gene_name', 'product', 'tf_family',
                  'module', 'bloc', 'coordination_T2', 'coordination_T3',
                  'LFC_T2vsT1', 'LFC_T3vsT1', 'LFC_T3vsT2',
                  'padj_T2vsT1', 'padj_T3vsT1', 'padj_T3vsT2',
                  'early_magnitude', 'late_magnitude', 'phase_ratio',
                  'temporal_class', 'direction', 'temporal_direction',
                  'earliness_score', 'trajectory_cluster',
                  'T1_mean', 'T2_mean', 'T3_mean',
                  'T1_z', 'T2_z', 'T3_z']].copy()
table1.to_csv(f'{TBL_DIR}/temporal_classification.tsv', sep='\t', index=False)
print(f"  [1/7] temporal_classification.tsv ({len(table1)} rows)")

# Table 2: Bloc comparison
bloc_comp = pd.DataFrame({
    'metric': ['n_genes', 'n_with_phase_ratio', 'mean_phase_ratio', 'median_phase_ratio',
               'std_phase_ratio', 'n_early', 'n_late', 'n_gradual', 'n_nonresponder',
               'mean_T1_z', 'mean_T2_z', 'mean_T3_z',
               'mean_LFC_T2vsT1', 'mean_LFC_T3vsT2', 'mean_LFC_T3vsT1'],
    'activation': [
        len(activation), len(act_pr), act_pr.mean(), act_pr.median(), act_pr.std(),
        (activation['temporal_class'] == 'early').sum(),
        (activation['temporal_class'] == 'late').sum(),
        (activation['temporal_class'] == 'gradual').sum(),
        (activation['temporal_class'] == 'non-responder').sum(),
        activation['T1_z'].mean(), activation['T2_z'].mean(), activation['T3_z'].mean(),
        activation['LFC_T2vsT1'].mean(), activation['LFC_T3vsT2'].mean(), activation['LFC_T3vsT1'].mean(),
    ],
    'repression': [
        len(repression), len(rep_pr), rep_pr.mean(), rep_pr.median(), rep_pr.std(),
        (repression['temporal_class'] == 'early').sum(),
        (repression['temporal_class'] == 'late').sum(),
        (repression['temporal_class'] == 'gradual').sum(),
        (repression['temporal_class'] == 'non-responder').sum(),
        repression['T1_z'].mean(), repression['T2_z'].mean(), repression['T3_z'].mean(),
        repression['LFC_T2vsT1'].mean(), repression['LFC_T3vsT2'].mean(), repression['LFC_T3vsT1'].mean(),
    ],
})
bloc_comp['MW_p_value'] = ''
bloc_comp.loc[bloc_comp['metric'] == 'mean_phase_ratio', 'MW_p_value'] = f'{p_w:.4e}'
bloc_comp.to_csv(f'{TBL_DIR}/bloc_comparison.tsv', sep='\t', index=False)
print(f"  [2/7] bloc_comparison.tsv")

# Table 3: Coordination type temporal
coord_results_df.to_csv(f'{TBL_DIR}/coordination_type_temporal.tsv', sep='\t', index=False)
print(f"  [3/7] coordination_type_temporal.tsv ({len(coord_results_df)} rows)")

# Table 4: TCS temporal analysis
tcs_temporal_df.to_csv(f'{TBL_DIR}/TCS_temporal_analysis.tsv', sep='\t', index=False)
print(f"  [4/7] TCS_temporal_analysis.tsv ({len(tcs_temporal_df)} rows)")

# Table 5: Methylation timing
methyl_timing_df.to_csv(f'{TBL_DIR}/methylation_timing.tsv', sep='\t', index=False)
print(f"  [5/7] methylation_timing.tsv ({len(methyl_timing_df)} rows)")

# Table 6: Statistical tests
stat_tests_df.to_csv(f'{TBL_DIR}/statistical_tests.tsv', sep='\t', index=False)
print(f"  [6/7] statistical_tests.tsv ({len(stat_tests_df)} rows)")

# Table 7: Exposed vs shielded temporal comparison
exp_vs_shi = pd.DataFrame({
    'metric': ['n_total', 'n_meaningful_change', 'n_with_phase_ratio',
               'mean_phase_ratio', 'median_phase_ratio', 'std_phase_ratio',
               'n_early', 'n_late', 'n_gradual', 'n_nonresponder',
               'pct_early', 'pct_late', 'pct_gradual',
               'early_fraction_of_responders'],
    'exposed': [
        len(master), master['has_meaningful_change'].sum(), len(exposed_pr),
        exposed_pr.mean(), exposed_pr.median(), exposed_pr.std(),
        exp_class.get('early', 0), exp_class.get('late', 0),
        exp_class.get('gradual', 0), exp_class.get('non-responder', 0),
        exp_class.get('early', 0)/len(master)*100,
        exp_class.get('late', 0)/len(master)*100,
        exp_class.get('gradual', 0)/len(master)*100,
        exp_early_frac,
    ],
    'shielded': [
        len(shielded_df), shielded_df['has_meaningful_change'].sum(), len(shielded_pr),
        shielded_pr.mean(), shielded_pr.median(), shielded_pr.std(),
        shi_class.get('early', 0), shi_class.get('late', 0),
        shi_class.get('gradual', 0), shi_class.get('non-responder', 0),
        shi_class.get('early', 0)/len(shielded_df)*100,
        shi_class.get('late', 0)/len(shielded_df)*100,
        shi_class.get('gradual', 0)/len(shielded_df)*100,
        shi_early_frac,
    ],
})
exp_vs_shi['test_p_value'] = ''
exp_vs_shi.loc[exp_vs_shi['metric'] == 'mean_phase_ratio', 'test_p_value'] = f'{p_es:.4e}'
exp_vs_shi.loc[exp_vs_shi['metric'] == 'pct_early', 'test_p_value'] = f'chi2 p={p_chi2:.4e}'
exp_vs_shi.to_csv(f'{TBL_DIR}/exposed_vs_shielded_temporal.tsv', sep='\t', index=False)
print(f"  [7/7] exposed_vs_shielded_temporal.tsv")

# ============================================================================
# Final summary
# ============================================================================
print("\n" + "=" * 70)
print("FINAL SUMMARY")
print("=" * 70)

print(f"""
H34 TEMPORAL DYNAMICS OF 62 EXPOSED TFs - KEY FINDINGS
=======================================================

1. TEMPORAL CLASSIFICATION:
   - Early responders (PR>0.6): {exp_class.get('early', 0)} ({exp_class.get('early', 0)/len(master)*100:.1f}%)
   - Late responders (PR<0.4): {exp_class.get('late', 0)} ({exp_class.get('late', 0)/len(master)*100:.1f}%)
   - Gradual (0.4-0.6): {exp_class.get('gradual', 0)} ({exp_class.get('gradual', 0)/len(master)*100:.1f}%)
   - Non-responder (|LFC|<0.5): {exp_class.get('non-responder', 0)} ({exp_class.get('non-responder', 0)/len(master)*100:.1f}%)
   - Mean phase ratio: {exposed_pr.mean():.3f} (median: {exposed_pr.median():.3f})

2. BLOC PHASE SEPARATION:
   - Activation bloc (n={len(activation)}): mean PR = {act_pr.mean():.3f}
   - Repression bloc (n={len(repression)}): mean PR = {rep_pr.mean():.3f}
   - Mann-Whitney: p = {p_w:.3e}, r = {r_rbs:.3f}
   - {'SIGNIFICANT' if p_w < 0.05 else 'NOT significant'} phase separation between blocs

3. COORDINATION TYPE TEMPORAL:
   - Kruskal-Wallis: p = {kw_p:.3e}
   - {'SIGNIFICANT' if kw_p < 0.05 else 'NOT significant'} differences across coordination types

4. METHYLATION-EXPRESSION TIMING:
   - Early methyl vs early expr: rho = {rho_early:.3f}, p = {p_early:.3e}
   - Late methyl vs late expr: rho = {rho_late:.3f}, p = {p_late:.3e}

5. TCS PAIR TEMPORAL ORDERING:
   - SK-exposed pairs: n={len(sk_exposed_pairs)}, RR-exposed: n={len(rr_exposed_pairs)}
   - Temporal lag (exposed-shielded): mean={temporal_lags.mean():.3f}

6. EXPOSED vs SHIELDED:
   - Phase ratio: exposed median={exposed_pr.median():.3f} vs shielded median={shielded_pr.median():.3f}
   - Mann-Whitney: p = {p_es:.3e}, r = {r_es:.3f}
   - Chi-square (class distribution): p = {p_chi2:.3e}

Files saved:
  Figures: {FIG_DIR}/ (7 figure sets, PDF+SVG)
  Tables: {TBL_DIR}/ (7 tables)
""")

print("H34 analysis complete.")
