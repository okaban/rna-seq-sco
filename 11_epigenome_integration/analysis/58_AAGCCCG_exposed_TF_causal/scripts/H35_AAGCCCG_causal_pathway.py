#!/usr/bin/env python3
"""
H35: AAGCCCG Exposed TF Causal Pathway Analysis

Tests whether restricting the AAGCCCG de-repression analysis to the 62 exposed
TFs (which lack protection zones) reveals a causal effect that was invisible
in the genome-wide test (H21: p=0.91).

Key insight: H21 tested all ~800 genes with AAGCCCG loss. But if only the 62
exposed TFs are susceptible (because they lack protection zones), the signal
would be diluted below detectability.
"""

import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# Configuration
# ============================================================================
BASE = Path('/Users/okaban/bioinfo/rna-seq')
ANALYSIS_DIR = BASE / '11_epigenome_integration/analysis/58_AAGCCCG_exposed_TF_causal'
FIG_DIR = ANALYSIS_DIR / 'figures'
TAB_DIR = ANALYSIS_DIR / 'tables'

# Input files
EXPOSED_TF_FILE = BASE / '11_epigenome_integration/analysis/51_exposed_regulators_characteristics/tables/exposed_regulators_full_table.tsv'
ALL_GENES_FILE = BASE / '11_epigenome_integration/analysis/52_shielded_exposed_boundary/tables/all_genes_features.tsv'
AAGCCCG_FILE = BASE / '11_epigenome_integration/analysis/36_AAGCCCG_distribution/tables/AAGCCCG_site_gene_mapping.tsv'
DESEQ_T2_FILE = BASE / '04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_2_vs_1.tsv'
DESEQ_T3_FILE = BASE / '04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_3_vs_1.tsv'
COUNTS_FILE = BASE / '04_deseq2/analysis/04_deseq2_260128_v1/results/normalized_counts_M145.tsv'
COORDINATED_FILE = BASE / '11_epigenome_integration/analysis/29_genomewide_TF_screen/tables/coordinated_regulatory_genes.tsv'
ANNOTATION_FILE = BASE / '05_annotation/analysis/05_annotation_260128_v1/tables/gene_annotation_basic.tsv'
GCCGGC_FILE = BASE / '11_epigenome_integration/analysis/37_defense_island_GCCGGC/tables/GCCGGC_sites_by_timepoint.tsv'

np.random.seed(42)
N_BOOTSTRAP = 10000

# Matplotlib style
plt.rcParams.update({
    'font.size': 10, 'axes.titlesize': 12, 'axes.labelsize': 11,
    'figure.dpi': 150, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
    'font.family': 'sans-serif',
})


# ============================================================================
# Step 1: Load data
# ============================================================================
print("=" * 80)
print("STEP 1: Loading data")
print("=" * 80)

exposed_tf = pd.read_csv(EXPOSED_TF_FILE, sep='\t')
all_genes = pd.read_csv(ALL_GENES_FILE, sep='\t')
aagcccg_map = pd.read_csv(AAGCCCG_FILE, sep='\t')
deseq_t2 = pd.read_csv(DESEQ_T2_FILE, sep='\t')
deseq_t3 = pd.read_csv(DESEQ_T3_FILE, sep='\t')
counts = pd.read_csv(COUNTS_FILE, sep='\t')
coordinated = pd.read_csv(COORDINATED_FILE, sep='\t')
annotation = pd.read_csv(ANNOTATION_FILE, sep='\t')
gccggc_sites = pd.read_csv(GCCGGC_FILE, sep='\t')

print(f"Exposed TFs: {len(exposed_tf)}")
print(f"All regulatory genes: {len(all_genes)}")
print(f"AAGCCCG site-gene mappings: {len(aagcccg_map)}")
print(f"GCCGGC sites: {len(gccggc_sites)}")
print(f"DESeq2 T2vsT1: {len(deseq_t2)} genes")
print(f"Normalized counts: {counts.shape}")

# Merge DESeq2 results
deseq_t2 = deseq_t2.rename(columns={'gene_id': 'locus_tag', 'log2FoldChange': 'LFC_T2vsT1', 'padj': 'padj_T2'})
deseq_t3 = deseq_t3.rename(columns={'gene_id': 'locus_tag', 'log2FoldChange': 'LFC_T3vsT1', 'padj': 'padj_T3'})

# ============================================================================
# Step 2: Map AAGCCCG and GCCGGC sites to exposed TFs
# ============================================================================
print("\n" + "=" * 80)
print("STEP 2: Mapping methylation sites to exposed TFs (TSS ±500bp)")
print("=" * 80)

# Get unique AAGCCCG positions per timepoint
aagcccg_t1_positions = set(aagcccg_map[aagcccg_map['timepoint'] == 'T1']['position'].unique())
aagcccg_t2_positions = set(aagcccg_map[aagcccg_map['timepoint'] == 'T2']['position'].unique())
print(f"AAGCCCG T1 unique positions: {len(aagcccg_t1_positions)}")
print(f"AAGCCCG T2 unique positions: {len(aagcccg_t2_positions)}")
print(f"Overlap: {len(aagcccg_t1_positions & aagcccg_t2_positions)}")

# Get unique GCCGGC positions per timepoint
gccggc_t1_positions = set(gccggc_sites[gccggc_sites['timepoint'] == 'T1']['position'].unique())
gccggc_t2_positions = set(gccggc_sites[gccggc_sites['timepoint'] == 'T2']['position'].unique())
gccggc_t3_positions = set(gccggc_sites[gccggc_sites['timepoint'] == 'T3']['position'].unique())
print(f"\nGCCGGC T1: {len(gccggc_t1_positions)}, T2: {len(gccggc_t2_positions)}, T3: {len(gccggc_t3_positions)}")

# Define TSS for exposed TFs
# For + strand: TSS = start; for - strand: TSS = end
exposed_tf['tss'] = exposed_tf.apply(
    lambda r: r['start'] if r['strand'] == '+' else r['end'], axis=1
)

# Function to find methylation sites near a TSS
def find_sites_near_tss(tss, positions_set, window=500):
    """Find positions within ±window of TSS"""
    return {p for p in positions_set if abs(p - tss) <= window}

# Map AAGCCCG sites to each exposed TF
results = []
for _, tf in exposed_tf.iterrows():
    tss = tf['tss']
    lt = tf['locus_tag']

    # AAGCCCG
    aag_t1 = find_sites_near_tss(tss, aagcccg_t1_positions)
    aag_t2 = find_sites_near_tss(tss, aagcccg_t2_positions)

    # GCCGGC
    gcc_t1 = find_sites_near_tss(tss, gccggc_t1_positions)
    gcc_t2 = find_sites_near_tss(tss, gccggc_t2_positions)
    gcc_t3 = find_sites_near_tss(tss, gccggc_t3_positions)

    # Classify AAGCCCG status
    n_aag_t1 = len(aag_t1)
    n_aag_t2 = len(aag_t2)
    # Since T1 and T2 positions don't overlap:
    # Lost = present at T1, absent at T2 (all T1 sites are lost)
    # Gained = present at T2, absent at T1 (all T2 sites are gained)
    # Never = no site at either timepoint
    if n_aag_t1 > 0 and n_aag_t2 == 0:
        aag_status = 'Lost'
    elif n_aag_t1 == 0 and n_aag_t2 > 0:
        aag_status = 'Gained'
    elif n_aag_t1 > 0 and n_aag_t2 > 0:
        aag_status = 'Both'  # T1 + T2 sites (different positions)
    else:
        aag_status = 'Never'

    # Classify GCCGGC status (T1 vs T2)
    n_gcc_t1 = len(gcc_t1)
    n_gcc_t2 = len(gcc_t2)
    gcc_t1_only = gcc_t1 - gcc_t2
    gcc_t2_only = gcc_t2 - gcc_t1
    gcc_both = gcc_t1 & gcc_t2

    if n_gcc_t1 > 0 and n_gcc_t2 == 0:
        gcc_status = 'Lost'
    elif n_gcc_t1 == 0 and n_gcc_t2 > 0:
        gcc_status = 'Gained'
    elif n_gcc_t1 > 0 and n_gcc_t2 > 0:
        gcc_status = 'Both'
    else:
        gcc_status = 'Never'

    results.append({
        'locus_tag': lt,
        'gene_name': tf.get('gene_name', ''),
        'old_locus_tag': tf.get('old_locus_tag', ''),
        'product': tf.get('product', ''),
        'tf_family': tf.get('tf_family', ''),
        'tss': tss,
        'strand': tf['strand'],
        'region': tf.get('region', ''),
        'n_AAGCCCG_T1': n_aag_t1,
        'n_AAGCCCG_T2': n_aag_t2,
        'n_AAGCCCG_lost': n_aag_t1,  # all T1 sites are lost
        'n_AAGCCCG_gained': n_aag_t2,  # all T2 sites are gained
        'AAGCCCG_status': aag_status,
        'AAGCCCG_T1_positions': ';'.join(map(str, sorted(aag_t1))) if aag_t1 else '',
        'AAGCCCG_T2_positions': ';'.join(map(str, sorted(aag_t2))) if aag_t2 else '',
        'n_GCCGGC_T1': n_gcc_t1,
        'n_GCCGGC_T2': n_gcc_t2,
        'n_GCCGGC_T3': len(gcc_t3),
        'n_GCCGGC_lost': len(gcc_t1_only),
        'n_GCCGGC_maintained': len(gcc_both),
        'GCCGGC_status': gcc_status,
        'LFC_T2vsT1': tf.get('log2FC_T2', np.nan),
        'LFC_T3vsT1': tf.get('log2FC_T3', np.nan),
        'padj_T2': tf.get('padj_T2', np.nan),
        'padj_T3': tf.get('padj_T3', np.nan),
        'baseMean': tf.get('baseMean', np.nan),
        'coordination_T2': tf.get('coordination_T2', ''),
        'coordination_T3': tf.get('coordination_T3', ''),
        'temporal_pattern': tf.get('temporal_pattern', ''),
    })

tf_methyl = pd.DataFrame(results)

# Summary
print(f"\nAAGCCCG status distribution (62 exposed TFs, ±500bp of TSS):")
aag_counts = tf_methyl['AAGCCCG_status'].value_counts()
for status, count in aag_counts.items():
    print(f"  {status}: {count} ({count/len(tf_methyl)*100:.1f}%)")

print(f"\nGCCGGC status distribution:")
gcc_counts = tf_methyl['GCCGGC_status'].value_counts()
for status, count in gcc_counts.items():
    print(f"  {status}: {count} ({count/len(tf_methyl)*100:.1f}%)")

# ============================================================================
# Step 3: Exposed-restricted de-repression test (KEY TEST)
# ============================================================================
print("\n" + "=" * 80)
print("STEP 3: Exposed-restricted de-repression test (KEY TEST)")
print("=" * 80)

# Split exposed TFs by AAGCCCG status
lost_tfs = tf_methyl[tf_methyl['AAGCCCG_status'] == 'Lost']
never_tfs = tf_methyl[tf_methyl['AAGCCCG_status'] == 'Never']
gained_tfs = tf_methyl[tf_methyl['AAGCCCG_status'] == 'Gained']
both_tfs = tf_methyl[tf_methyl['AAGCCCG_status'] == 'Both']

print(f"\nGroup sizes:")
print(f"  AAGCCCG-Lost: {len(lost_tfs)}")
print(f"  AAGCCCG-Never: {len(never_tfs)}")
print(f"  AAGCCCG-Gained: {len(gained_tfs)}")
print(f"  AAGCCCG-Both: {len(both_tfs)}")

stat_tests = []

# LFC_T2vsT1: Lost vs Never
if len(lost_tfs) >= 2 and len(never_tfs) >= 2:
    lfc_lost = lost_tfs['LFC_T2vsT1'].dropna()
    lfc_never = never_tfs['LFC_T2vsT1'].dropna()

    u_stat, p_val = stats.mannwhitneyu(lfc_lost, lfc_never, alternative='two-sided')
    # Effect size: rank-biserial correlation
    n1, n2 = len(lfc_lost), len(lfc_never)
    r_rb = 1 - (2 * u_stat) / (n1 * n2)

    # Cohen's d
    pooled_std = np.sqrt(((n1-1)*lfc_lost.std()**2 + (n2-1)*lfc_never.std()**2) / (n1+n2-2))
    cohens_d = (lfc_lost.mean() - lfc_never.mean()) / pooled_std if pooled_std > 0 else 0

    print(f"\n--- LFC_T2vsT1: Lost vs Never ---")
    print(f"  Lost (n={n1}):  median={lfc_lost.median():.3f}, mean={lfc_lost.mean():.3f}")
    print(f"  Never (n={n2}): median={lfc_never.median():.3f}, mean={lfc_never.mean():.3f}")
    print(f"  MWU U={u_stat:.1f}, p={p_val:.4e}")
    print(f"  Rank-biserial r={r_rb:.3f}")
    print(f"  Cohen's d={cohens_d:.3f}")

    stat_tests.append({
        'test': 'LFC_T2vsT1_Lost_vs_Never',
        'group1': 'AAGCCCG-Lost', 'n1': n1,
        'group2': 'AAGCCCG-Never', 'n2': n2,
        'group1_median': lfc_lost.median(), 'group1_mean': lfc_lost.mean(),
        'group2_median': lfc_never.median(), 'group2_mean': lfc_never.mean(),
        'statistic': u_stat, 'p_value': p_val,
        'effect_size_r': r_rb, 'cohens_d': cohens_d,
        'test_type': 'Mann-Whitney U', 'context': 'Exposed TFs only'
    })

# |LFC_T2vsT1|: Lost vs Never
if len(lost_tfs) >= 2 and len(never_tfs) >= 2:
    abs_lost = lfc_lost.abs()
    abs_never = lfc_never.abs()

    u_stat2, p_val2 = stats.mannwhitneyu(abs_lost, abs_never, alternative='two-sided')
    r_rb2 = 1 - (2 * u_stat2) / (len(abs_lost) * len(abs_never))

    print(f"\n--- |LFC_T2vsT1|: Lost vs Never ---")
    print(f"  Lost (n={len(abs_lost)}):  median={abs_lost.median():.3f}, mean={abs_lost.mean():.3f}")
    print(f"  Never (n={len(abs_never)}): median={abs_never.median():.3f}, mean={abs_never.mean():.3f}")
    print(f"  MWU U={u_stat2:.1f}, p={p_val2:.4e}")
    print(f"  Rank-biserial r={r_rb2:.3f}")

    stat_tests.append({
        'test': '|LFC_T2vsT1|_Lost_vs_Never',
        'group1': 'AAGCCCG-Lost', 'n1': len(abs_lost),
        'group2': 'AAGCCCG-Never', 'n2': len(abs_never),
        'group1_median': abs_lost.median(), 'group1_mean': abs_lost.mean(),
        'group2_median': abs_never.median(), 'group2_mean': abs_never.mean(),
        'statistic': u_stat2, 'p_value': p_val2,
        'effect_size_r': r_rb2, 'cohens_d': np.nan,
        'test_type': 'Mann-Whitney U', 'context': 'Exposed TFs only (absolute LFC)'
    })

# One-sided test: Lost > Never (de-repression prediction = upregulation)
if len(lost_tfs) >= 2 and len(never_tfs) >= 2:
    u_stat_one, p_val_one = stats.mannwhitneyu(lfc_lost, lfc_never, alternative='greater')
    print(f"\n--- LFC_T2vsT1: Lost > Never (one-sided, de-repression) ---")
    print(f"  MWU p={p_val_one:.4e}")

    stat_tests.append({
        'test': 'LFC_T2vsT1_Lost_gt_Never_onesided',
        'group1': 'AAGCCCG-Lost', 'n1': n1,
        'group2': 'AAGCCCG-Never', 'n2': n2,
        'group1_median': lfc_lost.median(), 'group1_mean': lfc_lost.mean(),
        'group2_median': lfc_never.median(), 'group2_mean': lfc_never.mean(),
        'statistic': u_stat_one, 'p_value': p_val_one,
        'effect_size_r': r_rb, 'cohens_d': cohens_d,
        'test_type': 'Mann-Whitney U (one-sided)', 'context': 'De-repression prediction'
    })

# Lost vs Gained (if sample size permits)
if len(gained_tfs) >= 2 and len(lost_tfs) >= 2:
    lfc_gained = gained_tfs['LFC_T2vsT1'].dropna()
    u_lg, p_lg = stats.mannwhitneyu(lfc_lost, lfc_gained, alternative='two-sided')
    r_lg = 1 - (2 * u_lg) / (len(lfc_lost) * len(lfc_gained))

    print(f"\n--- LFC_T2vsT1: Lost vs Gained ---")
    print(f"  Lost (n={len(lfc_lost)}):    median={lfc_lost.median():.3f}")
    print(f"  Gained (n={len(lfc_gained)}): median={lfc_gained.median():.3f}")
    print(f"  MWU p={p_lg:.4e}, r={r_lg:.3f}")

    stat_tests.append({
        'test': 'LFC_T2vsT1_Lost_vs_Gained',
        'group1': 'AAGCCCG-Lost', 'n1': len(lfc_lost),
        'group2': 'AAGCCCG-Gained', 'n2': len(lfc_gained),
        'group1_median': lfc_lost.median(), 'group1_mean': lfc_lost.mean(),
        'group2_median': lfc_gained.median(), 'group2_mean': lfc_gained.mean(),
        'statistic': u_lg, 'p_value': p_lg,
        'effect_size_r': r_lg, 'cohens_d': np.nan,
        'test_type': 'Mann-Whitney U', 'context': 'Exposed TFs only'
    })

# Kruskal-Wallis across all groups with n>=2
groups_kw = []
labels_kw = []
for status, grp in [('Lost', lost_tfs), ('Never', never_tfs), ('Gained', gained_tfs)]:
    vals = grp['LFC_T2vsT1'].dropna()
    if len(vals) >= 2:
        groups_kw.append(vals.values)
        labels_kw.append(f"{status} (n={len(vals)})")
if len(groups_kw) >= 2:
    kw_stat, kw_p = stats.kruskal(*groups_kw)
    print(f"\n--- Kruskal-Wallis across groups ---")
    print(f"  Groups: {', '.join(labels_kw)}")
    print(f"  H={kw_stat:.3f}, p={kw_p:.4e}")

    stat_tests.append({
        'test': 'KW_LFC_T2vsT1_all_groups',
        'group1': 'All AAGCCCG groups', 'n1': sum(len(g) for g in groups_kw),
        'group2': '', 'n2': 0,
        'group1_median': np.nan, 'group1_mean': np.nan,
        'group2_median': np.nan, 'group2_mean': np.nan,
        'statistic': kw_stat, 'p_value': kw_p,
        'effect_size_r': np.nan, 'cohens_d': np.nan,
        'test_type': 'Kruskal-Wallis', 'context': f'Groups: {", ".join(labels_kw)}'
    })


# ============================================================================
# Step 4: Same test for GCCGGC 4mC
# ============================================================================
print("\n" + "=" * 80)
print("STEP 4: GCCGGC 4mC de-repression test (comparison)")
print("=" * 80)

gcc_lost = tf_methyl[tf_methyl['GCCGGC_status'] == 'Lost']
gcc_never = tf_methyl[tf_methyl['GCCGGC_status'] == 'Never']
gcc_gained = tf_methyl[tf_methyl['GCCGGC_status'] == 'Gained']
gcc_both = tf_methyl[tf_methyl['GCCGGC_status'] == 'Both']

print(f"GCCGGC group sizes:")
print(f"  Lost: {len(gcc_lost)}")
print(f"  Never: {len(gcc_never)}")
print(f"  Gained: {len(gcc_gained)}")
print(f"  Both: {len(gcc_both)}")

if len(gcc_lost) >= 2 and len(gcc_never) >= 2:
    gcc_lfc_lost = gcc_lost['LFC_T2vsT1'].dropna()
    gcc_lfc_never = gcc_never['LFC_T2vsT1'].dropna()

    u_gcc, p_gcc = stats.mannwhitneyu(gcc_lfc_lost, gcc_lfc_never, alternative='two-sided')
    n1g, n2g = len(gcc_lfc_lost), len(gcc_lfc_never)
    r_gcc = 1 - (2 * u_gcc) / (n1g * n2g)
    pooled_g = np.sqrt(((n1g-1)*gcc_lfc_lost.std()**2 + (n2g-1)*gcc_lfc_never.std()**2) / (n1g+n2g-2))
    d_gcc = (gcc_lfc_lost.mean() - gcc_lfc_never.mean()) / pooled_g if pooled_g > 0 else 0

    print(f"\n--- LFC_T2vsT1: GCCGGC Lost vs Never ---")
    print(f"  Lost (n={n1g}):  median={gcc_lfc_lost.median():.3f}, mean={gcc_lfc_lost.mean():.3f}")
    print(f"  Never (n={n2g}): median={gcc_lfc_never.median():.3f}, mean={gcc_lfc_never.mean():.3f}")
    print(f"  MWU p={p_gcc:.4e}, r={r_gcc:.3f}, d={d_gcc:.3f}")

    stat_tests.append({
        'test': 'LFC_T2vsT1_GCCGGC_Lost_vs_Never',
        'group1': 'GCCGGC-Lost', 'n1': n1g,
        'group2': 'GCCGGC-Never', 'n2': n2g,
        'group1_median': gcc_lfc_lost.median(), 'group1_mean': gcc_lfc_lost.mean(),
        'group2_median': gcc_lfc_never.median(), 'group2_mean': gcc_lfc_never.mean(),
        'statistic': u_gcc, 'p_value': p_gcc,
        'effect_size_r': r_gcc, 'cohens_d': d_gcc,
        'test_type': 'Mann-Whitney U', 'context': 'GCCGGC comparison, Exposed TFs only'
    })
else:
    print("  Insufficient group sizes for GCCGGC test")
    # Still report what we have
    for status, grp in [('Lost', gcc_lost), ('Never', gcc_never), ('Gained', gcc_gained), ('Both', gcc_both)]:
        if len(grp) > 0:
            vals = grp['LFC_T2vsT1'].dropna()
            print(f"  {status} (n={len(grp)}): median LFC={vals.median():.3f}" if len(vals)>0 else f"  {status}: n={len(grp)}")


# ============================================================================
# Step 5: Dose-response within exposed TFs
# ============================================================================
print("\n" + "=" * 80)
print("STEP 5: Dose-response (n_AAGCCCG_sites_lost vs LFC_T2vsT1)")
print("=" * 80)

# Among all exposed TFs: n_AAGCCCG_lost vs LFC
dose_data = tf_methyl[['locus_tag', 'n_AAGCCCG_lost', 'LFC_T2vsT1']].dropna()
print(f"n_AAGCCCG_lost distribution:")
print(dose_data['n_AAGCCCG_lost'].value_counts().sort_index().to_string())

rho_dose, p_dose = stats.spearmanr(dose_data['n_AAGCCCG_lost'], dose_data['LFC_T2vsT1'])
print(f"\nSpearman: rho={rho_dose:.3f}, p={p_dose:.4e}")

stat_tests.append({
    'test': 'Dose_response_n_AAGCCCG_lost_vs_LFC',
    'group1': 'All exposed TFs', 'n1': len(dose_data),
    'group2': '', 'n2': 0,
    'group1_median': np.nan, 'group1_mean': np.nan,
    'group2_median': np.nan, 'group2_mean': np.nan,
    'statistic': rho_dose, 'p_value': p_dose,
    'effect_size_r': rho_dose, 'cohens_d': np.nan,
    'test_type': 'Spearman correlation', 'context': 'n_AAGCCCG_lost vs LFC_T2vsT1'
})

# Among TFs with at least 1 lost site
dose_with_sites = dose_data[dose_data['n_AAGCCCG_lost'] > 0]
if len(dose_with_sites) >= 5:
    rho_dose2, p_dose2 = stats.spearmanr(dose_with_sites['n_AAGCCCG_lost'], dose_with_sites['LFC_T2vsT1'])
    print(f"\nAmong TFs with >=1 lost site (n={len(dose_with_sites)}):")
    print(f"  Spearman: rho={rho_dose2:.3f}, p={p_dose2:.4e}")

    stat_tests.append({
        'test': 'Dose_response_with_sites_only',
        'group1': 'Exposed TFs with AAGCCCG', 'n1': len(dose_with_sites),
        'group2': '', 'n2': 0,
        'group1_median': np.nan, 'group1_mean': np.nan,
        'group2_median': np.nan, 'group2_mean': np.nan,
        'statistic': rho_dose2, 'p_value': p_dose2,
        'effect_size_r': rho_dose2, 'cohens_d': np.nan,
        'test_type': 'Spearman correlation', 'context': 'TFs with >=1 AAGCCCG lost only'
    })

# Also test for T3vsT1
dose_t3 = tf_methyl[['locus_tag', 'n_AAGCCCG_lost', 'LFC_T3vsT1']].dropna()
rho_dose_t3, p_dose_t3 = stats.spearmanr(dose_t3['n_AAGCCCG_lost'], dose_t3['LFC_T3vsT1'])
print(f"\nDose-response T3vsT1: rho={rho_dose_t3:.3f}, p={p_dose_t3:.4e}")


# ============================================================================
# Step 6: SC_RS17645 expression correlation
# ============================================================================
print("\n" + "=" * 80)
print("STEP 6: SC_RS17645 expression correlation with exposed vs shielded TFs")
print("=" * 80)

# Get SC_RS17645 expression across 9 samples
mtase_row = counts[counts['gene_id'] == 'SC_RS17645']
sample_cols = [c for c in counts.columns if c.startswith('M145_')]
mtase_expr = mtase_row[sample_cols].values.flatten()
print(f"SC_RS17645 expression: {mtase_expr}")

# For each exposed TF, compute Spearman correlation with SC_RS17645
exposed_lts = set(exposed_tf['locus_tag'])
shielded_genes = all_genes[all_genes['is_exposed'] == 0]
shielded_lts = set(shielded_genes['locus_tag'])

corr_results = []
for lt in exposed_lts | shielded_lts:
    gene_row = counts[counts['gene_id'] == lt]
    if len(gene_row) == 0:
        continue
    gene_expr = gene_row[sample_cols].values.flatten()
    if np.std(gene_expr) == 0:
        continue
    rho, p = stats.spearmanr(mtase_expr, gene_expr)
    is_exp = lt in exposed_lts
    corr_results.append({
        'locus_tag': lt,
        'is_exposed': is_exp,
        'spearman_rho': rho,
        'p_value': p,
        'mtase_correlation_direction': 'negative' if rho < 0 else 'positive'
    })

corr_df = pd.DataFrame(corr_results)
corr_exposed = corr_df[corr_df['is_exposed'] == True]
corr_shielded = corr_df[corr_df['is_exposed'] == False]

print(f"\nSC_RS17645 correlations:")
print(f"  With exposed TFs (n={len(corr_exposed)}): median rho={corr_exposed['spearman_rho'].median():.3f}")
print(f"  With shielded TFs (n={len(corr_shielded)}): median rho={corr_shielded['spearman_rho'].median():.3f}")

# Are exposed TF correlations more negative than shielded?
u_corr, p_corr = stats.mannwhitneyu(
    corr_exposed['spearman_rho'].dropna(),
    corr_shielded['spearman_rho'].dropna(),
    alternative='less'  # exposed should be more negative
)
r_corr = 1 - (2 * u_corr) / (len(corr_exposed) * len(corr_shielded))

print(f"  MWU (exposed < shielded): p={p_corr:.4e}, r={r_corr:.3f}")

# Two-sided test too
u_corr2, p_corr2 = stats.mannwhitneyu(
    corr_exposed['spearman_rho'].dropna(),
    corr_shielded['spearman_rho'].dropna(),
    alternative='two-sided'
)

stat_tests.append({
    'test': 'SC_RS17645_corr_exposed_vs_shielded',
    'group1': 'Exposed TFs', 'n1': len(corr_exposed),
    'group2': 'Shielded TFs', 'n2': len(corr_shielded),
    'group1_median': corr_exposed['spearman_rho'].median(),
    'group1_mean': corr_exposed['spearman_rho'].mean(),
    'group2_median': corr_shielded['spearman_rho'].median(),
    'group2_mean': corr_shielded['spearman_rho'].mean(),
    'statistic': u_corr2, 'p_value': p_corr2,
    'effect_size_r': r_corr, 'cohens_d': np.nan,
    'test_type': 'Mann-Whitney U (two-sided)', 'context': 'Correlation with SC_RS17645 MTase'
})

# Fraction of negative correlations
frac_neg_exposed = (corr_exposed['spearman_rho'] < 0).mean()
frac_neg_shielded = (corr_shielded['spearman_rho'] < 0).mean()
print(f"  Fraction negative: exposed={frac_neg_exposed:.3f}, shielded={frac_neg_shielded:.3f}")

# Among exposed TFs, split by AAGCCCG status
corr_with_methyl = corr_exposed.merge(tf_methyl[['locus_tag', 'AAGCCCG_status']], on='locus_tag')
for status in ['Lost', 'Never', 'Gained']:
    sub = corr_with_methyl[corr_with_methyl['AAGCCCG_status'] == status]
    if len(sub) > 0:
        print(f"  Exposed AAGCCCG-{status} (n={len(sub)}): median rho={sub['spearman_rho'].median():.3f}")


# ============================================================================
# Step 7: Temporal specificity — T1→T2 vs T2→T3
# ============================================================================
print("\n" + "=" * 80)
print("STEP 7: Temporal specificity (T1→T2 vs T2→T3)")
print("=" * 80)

# Compute T3vsT2 as LFC_T3vsT1 - LFC_T2vsT1
tf_methyl['LFC_T3vsT2'] = tf_methyl['LFC_T3vsT1'] - tf_methyl['LFC_T2vsT1']

# AAGCCCG Lost vs Never: effect in T2vsT1 vs T3vsT2
if len(lost_tfs) >= 2 and len(never_tfs) >= 2:
    # Recalculate T3vsT2 for lost and never
    lost_t3t2 = tf_methyl[tf_methyl['AAGCCCG_status'] == 'Lost']['LFC_T3vsT2'].dropna()
    never_t3t2 = tf_methyl[tf_methyl['AAGCCCG_status'] == 'Never']['LFC_T3vsT2'].dropna()

    u_t3t2, p_t3t2 = stats.mannwhitneyu(lost_t3t2, never_t3t2, alternative='two-sided')
    r_t3t2 = 1 - (2 * u_t3t2) / (len(lost_t3t2) * len(never_t3t2))

    print(f"\nT1→T2 transition (when SC_RS17645 drops):")
    print(f"  Lost vs Never median LFC diff: {lfc_lost.median() - lfc_never.median():.3f}")
    print(f"  p={p_val:.4e}, r={r_rb:.3f}")

    print(f"\nT2→T3 transition (after SC_RS17645 drop):")
    print(f"  Lost (n={len(lost_t3t2)}): median LFC_T3vsT2={lost_t3t2.median():.3f}")
    print(f"  Never (n={len(never_t3t2)}): median LFC_T3vsT2={never_t3t2.median():.3f}")
    print(f"  p={p_t3t2:.4e}, r={r_t3t2:.3f}")

    stat_tests.append({
        'test': 'LFC_T3vsT2_Lost_vs_Never',
        'group1': 'AAGCCCG-Lost', 'n1': len(lost_t3t2),
        'group2': 'AAGCCCG-Never', 'n2': len(never_t3t2),
        'group1_median': lost_t3t2.median(), 'group1_mean': lost_t3t2.mean(),
        'group2_median': never_t3t2.median(), 'group2_mean': never_t3t2.mean(),
        'statistic': u_t3t2, 'p_value': p_t3t2,
        'effect_size_r': r_t3t2, 'cohens_d': np.nan,
        'test_type': 'Mann-Whitney U', 'context': 'Temporal specificity T2→T3'
    })

    # Also T3vsT1
    lost_t3t1 = tf_methyl[tf_methyl['AAGCCCG_status'] == 'Lost']['LFC_T3vsT1'].dropna()
    never_t3t1 = tf_methyl[tf_methyl['AAGCCCG_status'] == 'Never']['LFC_T3vsT1'].dropna()
    u_t3t1, p_t3t1 = stats.mannwhitneyu(lost_t3t1, never_t3t1, alternative='two-sided')
    r_t3t1 = 1 - (2 * u_t3t1) / (len(lost_t3t1) * len(never_t3t1))

    print(f"\nT3vsT1 (cumulative):")
    print(f"  Lost: median={lost_t3t1.median():.3f}")
    print(f"  Never: median={never_t3t1.median():.3f}")
    print(f"  p={p_t3t1:.4e}, r={r_t3t1:.3f}")

    stat_tests.append({
        'test': 'LFC_T3vsT1_Lost_vs_Never',
        'group1': 'AAGCCCG-Lost', 'n1': len(lost_t3t1),
        'group2': 'AAGCCCG-Never', 'n2': len(never_t3t1),
        'group1_median': lost_t3t1.median(), 'group1_mean': lost_t3t1.mean(),
        'group2_median': never_t3t1.median(), 'group2_mean': never_t3t1.mean(),
        'statistic': u_t3t1, 'p_value': p_t3t1,
        'effect_size_r': r_t3t1, 'cohens_d': np.nan,
        'test_type': 'Mann-Whitney U', 'context': 'Cumulative T3vsT1'
    })


# ============================================================================
# Step 8: Bootstrap confidence intervals
# ============================================================================
print("\n" + "=" * 80)
print("STEP 8: Bootstrap confidence intervals (10,000 iterations)")
print("=" * 80)

def bootstrap_median_diff(group1, group2, n_boot=N_BOOTSTRAP):
    """Bootstrap CI for median difference (group1 - group2)."""
    diffs = []
    for _ in range(n_boot):
        s1 = np.random.choice(group1, size=len(group1), replace=True)
        s2 = np.random.choice(group2, size=len(group2), replace=True)
        diffs.append(np.median(s1) - np.median(s2))
    diffs = np.array(diffs)
    ci_low = np.percentile(diffs, 2.5)
    ci_high = np.percentile(diffs, 97.5)
    return np.median(diffs), ci_low, ci_high, diffs

def bootstrap_mean_diff(group1, group2, n_boot=N_BOOTSTRAP):
    """Bootstrap CI for mean difference."""
    diffs = []
    for _ in range(n_boot):
        s1 = np.random.choice(group1, size=len(group1), replace=True)
        s2 = np.random.choice(group2, size=len(group2), replace=True)
        diffs.append(np.mean(s1) - np.mean(s2))
    diffs = np.array(diffs)
    return np.median(diffs), np.percentile(diffs, 2.5), np.percentile(diffs, 97.5), diffs

bootstrap_results = {}

if len(lfc_lost) >= 2 and len(lfc_never) >= 2:
    # Median diff
    med_diff, ci_lo, ci_hi, boot_diffs = bootstrap_median_diff(lfc_lost.values, lfc_never.values)
    print(f"\nLFC_T2vsT1 Lost vs Never (median diff):")
    print(f"  Bootstrap median diff: {med_diff:.3f}")
    print(f"  95% CI: [{ci_lo:.3f}, {ci_hi:.3f}]")
    print(f"  CI includes 0: {ci_lo <= 0 <= ci_hi}")
    bootstrap_results['median_diff_T2vsT1'] = (med_diff, ci_lo, ci_hi)

    # Mean diff
    mean_diff, ci_lo_m, ci_hi_m, _ = bootstrap_mean_diff(lfc_lost.values, lfc_never.values)
    print(f"\nLFC_T2vsT1 Lost vs Never (mean diff):")
    print(f"  Bootstrap mean diff: {mean_diff:.3f}")
    print(f"  95% CI: [{ci_lo_m:.3f}, {ci_hi_m:.3f}]")
    bootstrap_results['mean_diff_T2vsT1'] = (mean_diff, ci_lo_m, ci_hi_m)

    # Bootstrap for |LFC| diff
    med_abs, ci_lo_a, ci_hi_a, _ = bootstrap_median_diff(lfc_lost.abs().values, lfc_never.abs().values)
    print(f"\n|LFC_T2vsT1| Lost vs Never (median diff):")
    print(f"  Bootstrap: {med_abs:.3f}, 95% CI: [{ci_lo_a:.3f}, {ci_hi_a:.3f}]")
    bootstrap_results['median_abs_diff_T2vsT1'] = (med_abs, ci_lo_a, ci_hi_a)

    # Bootstrap p-value for median diff (fraction of bootstraps <= 0)
    p_boot = np.mean(boot_diffs <= 0)
    print(f"  Bootstrap p (fraction diffs <= 0): {p_boot:.4f}")

# Temporal: T3vsT2
if len(lost_t3t2) >= 2 and len(never_t3t2) >= 2:
    med_t3t2, ci_lo_t3, ci_hi_t3, _ = bootstrap_median_diff(lost_t3t2.values, never_t3t2.values)
    print(f"\nLFC_T3vsT2 Lost vs Never (median diff):")
    print(f"  Bootstrap: {med_t3t2:.3f}, 95% CI: [{ci_lo_t3:.3f}, {ci_hi_t3:.3f}]")
    bootstrap_results['median_diff_T3vsT2'] = (med_t3t2, ci_lo_t3, ci_hi_t3)


# ============================================================================
# Step 9: Comparison with H21 genome-wide result
# ============================================================================
print("\n" + "=" * 80)
print("STEP 9: Comparison with H21 genome-wide result")
print("=" * 80)

# Get all genes with DESeq2 data
deseq_t2_full = deseq_t2[['locus_tag', 'LFC_T2vsT1', 'padj_T2']].copy()

# Get unique AAGCCCG sites per gene (within 500bp of TSS)
# We need to get the TSS for ALL genes
annot = annotation[['gene_id', 'start', 'end', 'strand']].copy()
annot = annot.rename(columns={'gene_id': 'locus_tag'})
annot['tss'] = annot.apply(lambda r: r['start'] if r['strand'] == '+' else r['end'], axis=1)

# For each gene, find AAGCCCG sites within 500bp of TSS
gene_aag_status = {}
for _, row in annot.iterrows():
    lt = row['locus_tag']
    tss = row['tss']
    t1_nearby = find_sites_near_tss(tss, aagcccg_t1_positions, window=500)
    t2_nearby = find_sites_near_tss(tss, aagcccg_t2_positions, window=500)
    if len(t1_nearby) > 0 and len(t2_nearby) == 0:
        gene_aag_status[lt] = 'Lost'
    elif len(t1_nearby) == 0 and len(t2_nearby) > 0:
        gene_aag_status[lt] = 'Gained'
    elif len(t1_nearby) > 0 and len(t2_nearby) > 0:
        gene_aag_status[lt] = 'Both'
    else:
        gene_aag_status[lt] = 'Never'

aag_status_df = pd.DataFrame([
    {'locus_tag': lt, 'AAGCCCG_status_genomewide': status}
    for lt, status in gene_aag_status.items()
])

# Merge with DESeq2
gw_data = deseq_t2_full.merge(aag_status_df, on='locus_tag', how='inner')

# Get regulatory gene sets
all_reg_lts = set(all_genes['locus_tag'])
exposed_reg_lts = set(all_genes[all_genes['is_exposed'] == 1]['locus_tag'])
shielded_reg_lts = set(all_genes[all_genes['is_exposed'] == 0]['locus_tag'])

h21_comparisons = []

for label, gene_set in [
    ('All genes (genome-wide)', None),
    ('All regulatory (1,017)', all_reg_lts),
    ('Shielded regulatory (955)', shielded_reg_lts),
    ('Exposed regulatory (62)', exposed_reg_lts)
]:
    if gene_set is not None:
        subset = gw_data[gw_data['locus_tag'].isin(gene_set)]
    else:
        subset = gw_data

    lost_g = subset[subset['AAGCCCG_status_genomewide'] == 'Lost']['LFC_T2vsT1'].dropna()
    never_g = subset[subset['AAGCCCG_status_genomewide'] == 'Never']['LFC_T2vsT1'].dropna()

    if len(lost_g) >= 2 and len(never_g) >= 2:
        u, p = stats.mannwhitneyu(lost_g, never_g, alternative='two-sided')
        r = 1 - (2 * u) / (len(lost_g) * len(never_g))
        pooled = np.sqrt(((len(lost_g)-1)*lost_g.std()**2 + (len(never_g)-1)*never_g.std()**2) / (len(lost_g)+len(never_g)-2))
        d = (lost_g.mean() - never_g.mean()) / pooled if pooled > 0 else 0

        print(f"\n{label}:")
        print(f"  Lost (n={len(lost_g)}): median={lost_g.median():.3f}, mean={lost_g.mean():.3f}")
        print(f"  Never (n={len(never_g)}): median={never_g.median():.3f}, mean={never_g.mean():.3f}")
        print(f"  MWU p={p:.4e}, r={r:.3f}, d={d:.3f}")

        h21_comparisons.append({
            'gene_set': label,
            'n_lost': len(lost_g), 'n_never': len(never_g),
            'median_lost': lost_g.median(), 'mean_lost': lost_g.mean(),
            'median_never': never_g.median(), 'mean_never': never_g.mean(),
            'median_diff': lost_g.median() - never_g.median(),
            'MWU_p': p, 'rank_biserial_r': r, 'cohens_d': d
        })

        stat_tests.append({
            'test': f'H21_comparison_{label[:20]}',
            'group1': 'AAGCCCG-Lost', 'n1': len(lost_g),
            'group2': 'AAGCCCG-Never', 'n2': len(never_g),
            'group1_median': lost_g.median(), 'group1_mean': lost_g.mean(),
            'group2_median': never_g.median(), 'group2_mean': never_g.mean(),
            'statistic': u, 'p_value': p,
            'effect_size_r': r, 'cohens_d': d,
            'test_type': 'Mann-Whitney U', 'context': f'H21 comparison: {label}'
        })
    else:
        print(f"\n{label}: insufficient data (Lost={len(lost_g)}, Never={len(never_g)})")
        h21_comparisons.append({
            'gene_set': label,
            'n_lost': len(lost_g), 'n_never': len(never_g),
            'median_lost': lost_g.median() if len(lost_g) > 0 else np.nan,
            'mean_lost': lost_g.mean() if len(lost_g) > 0 else np.nan,
            'median_never': never_g.median() if len(never_g) > 0 else np.nan,
            'mean_never': never_g.mean() if len(never_g) > 0 else np.nan,
            'median_diff': np.nan, 'MWU_p': np.nan,
            'rank_biserial_r': np.nan, 'cohens_d': np.nan
        })

h21_comp_df = pd.DataFrame(h21_comparisons)


# ============================================================================
# Step 10: Figures
# ============================================================================
print("\n" + "=" * 80)
print("STEP 10: Generating figures")
print("=" * 80)

# Color palette
COLORS = {
    'Lost': '#e74c3c',
    'Never': '#95a5a6',
    'Gained': '#3498db',
    'Both': '#f39c12',
    'Maintained': '#2ecc71',
    'exposed': '#e74c3c',
    'shielded': '#3498db',
}

# ---------- Figure 1: Exposed vs genome-wide de-repression ----------
fig, axes = plt.subplots(1, 4, figsize=(16, 5))

for idx, row in h21_comp_df.iterrows():
    ax = axes[idx]
    label = row['gene_set']
    short_label = label.split('(')[0].strip()

    if gene_set is not None:
        subset = gw_data[gw_data['locus_tag'].isin(
            all_reg_lts if 'All reg' in label else
            shielded_reg_lts if 'Shielded' in label else
            exposed_reg_lts if 'Exposed' in label else set()
        )]
    else:
        subset = gw_data

    # Get the data for this panel
    if idx == 0:
        subset = gw_data
    elif idx == 1:
        subset = gw_data[gw_data['locus_tag'].isin(all_reg_lts)]
    elif idx == 2:
        subset = gw_data[gw_data['locus_tag'].isin(shielded_reg_lts)]
    else:
        subset = gw_data[gw_data['locus_tag'].isin(exposed_reg_lts)]

    lost_vals = subset[subset['AAGCCCG_status_genomewide'] == 'Lost']['LFC_T2vsT1'].dropna()
    never_vals = subset[subset['AAGCCCG_status_genomewide'] == 'Never']['LFC_T2vsT1'].dropna()

    box_data = []
    box_labels = []
    box_colors = []
    if len(lost_vals) > 0:
        box_data.append(lost_vals.values)
        box_labels.append(f'Lost\n(n={len(lost_vals)})')
        box_colors.append(COLORS['Lost'])
    if len(never_vals) > 0:
        box_data.append(never_vals.values)
        box_labels.append(f'Never\n(n={len(never_vals)})')
        box_colors.append(COLORS['Never'])

    if len(box_data) >= 2:
        bp = ax.boxplot(box_data, positions=[1, 2], widths=0.6, patch_artist=True,
                       showfliers=True, flierprops={'markersize': 3, 'alpha': 0.3})
        for patch, color in zip(bp['boxes'], box_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)
        for med in bp['medians']:
            med.set_color('black')
            med.set_linewidth(2)

        ax.set_xticks([1, 2])
        ax.set_xticklabels(box_labels, fontsize=9)

        p_str = f"p={row['MWU_p']:.2e}" if not np.isnan(row['MWU_p']) else "N/A"
        r_str = f"r={row['rank_biserial_r']:.3f}" if not np.isnan(row['rank_biserial_r']) else ""
        ax.set_title(f"{short_label}\n{p_str}, {r_str}", fontsize=10)
    else:
        ax.text(0.5, 0.5, 'Insufficient data', ha='center', va='center', transform=ax.transAxes)
        ax.set_title(short_label, fontsize=10)

    ax.set_ylabel('LFC T2 vs T1' if idx == 0 else '')
    ax.axhline(y=0, color='black', linestyle='--', alpha=0.3, linewidth=0.8)
    ax.set_ylim(-6, 6)

fig.suptitle('H35: AAGCCCG De-repression Effect — Genome-wide vs Exposed TFs', fontsize=13, y=1.02)
plt.tight_layout()
for fmt in ['pdf', 'svg']:
    fig.savefig(FIG_DIR / f'exposed_vs_genomewide_derepression.{fmt}', bbox_inches='tight')
plt.close()
print("  Saved: exposed_vs_genomewide_derepression.pdf/svg")


# ---------- Figure 2: AAGCCCG Lost vs Never vs Gained LFC (exposed only) ----------
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Panel A: Box/violin for LFC_T2vsT1
ax = axes[0]
groups_to_plot = []
labels_to_plot = []
colors_to_plot = []
for status, color in [('Lost', COLORS['Lost']), ('Never', COLORS['Never']), ('Gained', COLORS['Gained'])]:
    vals = tf_methyl[tf_methyl['AAGCCCG_status'] == status]['LFC_T2vsT1'].dropna()
    if len(vals) >= 1:
        groups_to_plot.append(vals.values)
        labels_to_plot.append(f'{status}\n(n={len(vals)})')
        colors_to_plot.append(color)

positions = list(range(1, len(groups_to_plot) + 1))

# Violin plot
parts = ax.violinplot(groups_to_plot, positions=positions, showmeans=False, showmedians=False, showextrema=False)
for pc, color in zip(parts['bodies'], colors_to_plot):
    pc.set_facecolor(color)
    pc.set_alpha(0.3)

# Overlay box plot
bp = ax.boxplot(groups_to_plot, positions=positions, widths=0.3, patch_artist=True,
               showfliers=False)
for patch, color in zip(bp['boxes'], colors_to_plot):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
for med in bp['medians']:
    med.set_color('black')
    med.set_linewidth(2)

# Scatter individual points
for i, (vals, pos) in enumerate(zip(groups_to_plot, positions)):
    jitter = np.random.normal(0, 0.05, len(vals))
    ax.scatter(pos + jitter, vals, c=colors_to_plot[i], s=25, alpha=0.6, zorder=3, edgecolors='white', linewidth=0.5)

ax.set_xticks(positions)
ax.set_xticklabels(labels_to_plot)
ax.set_ylabel('log2FC (T2 vs T1)')
ax.axhline(y=0, color='black', linestyle='--', alpha=0.3)
ax.set_title('A) LFC distribution by AAGCCCG status\n(Exposed TFs only)')

# Panel B: LFC_T3vsT1
ax = axes[1]
groups_t3 = []
labels_t3 = []
colors_t3 = []
for status, color in [('Lost', COLORS['Lost']), ('Never', COLORS['Never']), ('Gained', COLORS['Gained'])]:
    vals = tf_methyl[tf_methyl['AAGCCCG_status'] == status]['LFC_T3vsT1'].dropna()
    if len(vals) >= 1:
        groups_t3.append(vals.values)
        labels_t3.append(f'{status}\n(n={len(vals)})')
        colors_t3.append(color)

positions_t3 = list(range(1, len(groups_t3) + 1))
parts_t3 = ax.violinplot(groups_t3, positions=positions_t3, showmeans=False, showmedians=False, showextrema=False)
for pc, color in zip(parts_t3['bodies'], colors_t3):
    pc.set_facecolor(color)
    pc.set_alpha(0.3)
bp_t3 = ax.boxplot(groups_t3, positions=positions_t3, widths=0.3, patch_artist=True, showfliers=False)
for patch, color in zip(bp_t3['boxes'], colors_t3):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
for med in bp_t3['medians']:
    med.set_color('black')
    med.set_linewidth(2)
for i, (vals, pos) in enumerate(zip(groups_t3, positions_t3)):
    jitter = np.random.normal(0, 0.05, len(vals))
    ax.scatter(pos + jitter, vals, c=colors_t3[i], s=25, alpha=0.6, zorder=3, edgecolors='white', linewidth=0.5)

ax.set_xticks(positions_t3)
ax.set_xticklabels(labels_t3)
ax.set_ylabel('log2FC (T3 vs T1)')
ax.axhline(y=0, color='black', linestyle='--', alpha=0.3)
ax.set_title('B) LFC distribution by AAGCCCG status\n(T3 vs T1, Exposed TFs only)')

plt.tight_layout()
for fmt in ['pdf', 'svg']:
    fig.savefig(FIG_DIR / f'AAGCCCG_lost_vs_never_LFC.{fmt}', bbox_inches='tight')
plt.close()
print("  Saved: AAGCCCG_lost_vs_never_LFC.pdf/svg")


# ---------- Figure 3: Dose-response ----------
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Panel A: n_AAGCCCG_lost vs LFC_T2vsT1 (all exposed TFs)
ax = axes[0]
x = dose_data['n_AAGCCCG_lost'].values
y = dose_data['LFC_T2vsT1'].values
colors = [COLORS['Lost'] if xi > 0 else COLORS['Never'] for xi in x]
ax.scatter(x, y, c=colors, s=40, alpha=0.7, edgecolors='white', linewidth=0.5)

# Add trend line
if len(x) > 2:
    z = np.polyfit(x, y, 1)
    x_line = np.linspace(x.min(), x.max(), 100)
    ax.plot(x_line, np.polyval(z, x_line), 'k--', alpha=0.5, linewidth=1.5)

ax.set_xlabel('Number of AAGCCCG sites lost (±500bp TSS)')
ax.set_ylabel('log2FC (T2 vs T1)')
ax.axhline(y=0, color='grey', linestyle=':', alpha=0.5)
ax.set_title(f'A) Dose-response (all 62 exposed TFs)\nrho={rho_dose:.3f}, p={p_dose:.2e}')

# Panel B: Among TFs with sites
ax = axes[1]
if len(dose_with_sites) >= 3:
    x2 = dose_with_sites['n_AAGCCCG_lost'].values
    y2 = dose_with_sites['LFC_T2vsT1'].values
    ax.scatter(x2, y2, c=COLORS['Lost'], s=50, alpha=0.7, edgecolors='white', linewidth=0.5)

    if len(x2) > 2:
        z2 = np.polyfit(x2, y2, 1)
        x_line2 = np.linspace(x2.min(), x2.max(), 100)
        ax.plot(x_line2, np.polyval(z2, x_line2), 'k--', alpha=0.5, linewidth=1.5)

    ax.set_title(f'B) Dose-response (TFs with AAGCCCG only)\nrho={rho_dose2:.3f}, p={p_dose2:.2e}')
else:
    ax.text(0.5, 0.5, 'Insufficient data', ha='center', va='center', transform=ax.transAxes)
    ax.set_title('B) TFs with AAGCCCG sites only')

ax.set_xlabel('Number of AAGCCCG sites lost (±500bp TSS)')
ax.set_ylabel('log2FC (T2 vs T1)')
ax.axhline(y=0, color='grey', linestyle=':', alpha=0.5)

plt.tight_layout()
for fmt in ['pdf', 'svg']:
    fig.savefig(FIG_DIR / f'dose_response.{fmt}', bbox_inches='tight')
plt.close()
print("  Saved: dose_response.pdf/svg")


# ---------- Figure 4: SC_RS17645 correlation ----------
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Panel A: Distribution of correlations
ax = axes[0]
ax.hist(corr_exposed['spearman_rho'].dropna(), bins=20, alpha=0.7, color=COLORS['exposed'],
        label=f'Exposed (n={len(corr_exposed)})', density=True, edgecolor='white')
ax.hist(corr_shielded['spearman_rho'].dropna(), bins=30, alpha=0.4, color=COLORS['shielded'],
        label=f'Shielded (n={len(corr_shielded)})', density=True, edgecolor='white')
ax.axvline(corr_exposed['spearman_rho'].median(), color=COLORS['exposed'], linestyle='--', linewidth=2)
ax.axvline(corr_shielded['spearman_rho'].median(), color=COLORS['shielded'], linestyle='--', linewidth=2)
ax.set_xlabel('Spearman rho with SC_RS17645')
ax.set_ylabel('Density')
ax.legend(fontsize=9)
ax.set_title(f'A) MTase-TF expression correlation\np={p_corr2:.2e}')

# Panel B: By AAGCCCG status within exposed
ax = axes[1]
for status, color in [('Lost', COLORS['Lost']), ('Never', COLORS['Never']), ('Gained', COLORS['Gained'])]:
    sub = corr_with_methyl[corr_with_methyl['AAGCCCG_status'] == status]
    if len(sub) > 0:
        ax.hist(sub['spearman_rho'].dropna(), bins=12, alpha=0.6, color=color,
                label=f'{status} (n={len(sub)})', density=True, edgecolor='white')
ax.set_xlabel('Spearman rho with SC_RS17645')
ax.set_ylabel('Density')
ax.legend(fontsize=9)
ax.set_title('B) By AAGCCCG status (exposed only)')

# Panel C: Scatter - top examples
ax = axes[2]
# Show a few example exposed TFs
top_neg = corr_exposed.nsmallest(3, 'spearman_rho')
top_pos = corr_exposed.nlargest(3, 'spearman_rho')

for _, row in pd.concat([top_neg, top_pos]).iterrows():
    lt = row['locus_tag']
    gene_row = counts[counts['gene_id'] == lt]
    if len(gene_row) == 0:
        continue
    gene_expr = gene_row[sample_cols].values.flatten()
    color = COLORS['Lost'] if row['spearman_rho'] < 0 else COLORS['Gained']
    ax.scatter(mtase_expr, gene_expr, s=20, alpha=0.7, label=f"{lt} (rho={row['spearman_rho']:.2f})")

ax.set_xlabel('SC_RS17645 expression')
ax.set_ylabel('TF expression')
ax.legend(fontsize=7, loc='upper left')
ax.set_title('C) Example exposed TF correlations')

plt.tight_layout()
for fmt in ['pdf', 'svg']:
    fig.savefig(FIG_DIR / f'SC_RS17645_correlation.{fmt}', bbox_inches='tight')
plt.close()
print("  Saved: SC_RS17645_correlation.pdf/svg")


# ---------- Figure 5: Temporal specificity ----------
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Panel A: T2vsT1 effect
ax = axes[0]
for status, color in [('Lost', COLORS['Lost']), ('Never', COLORS['Never'])]:
    vals = tf_methyl[tf_methyl['AAGCCCG_status'] == status]['LFC_T2vsT1'].dropna()
    if len(vals) > 0:
        ax.hist(vals, bins=12, alpha=0.6, color=color, label=f'{status} (n={len(vals)})',
                density=True, edgecolor='white')
ax.axvline(x=0, color='black', linestyle='--', alpha=0.3)
ax.set_xlabel('log2FC (T2 vs T1)')
ax.set_ylabel('Density')
ax.legend(fontsize=9)
ax.set_title(f'A) T1→T2 (MTase drops)\np={p_val:.2e}')

# Panel B: T3vsT2 effect
ax = axes[1]
for status, color in [('Lost', COLORS['Lost']), ('Never', COLORS['Never'])]:
    vals = tf_methyl[tf_methyl['AAGCCCG_status'] == status]['LFC_T3vsT2'].dropna()
    if len(vals) > 0:
        ax.hist(vals, bins=12, alpha=0.6, color=color, label=f'{status} (n={len(vals)})',
                density=True, edgecolor='white')
ax.axvline(x=0, color='black', linestyle='--', alpha=0.3)
ax.set_xlabel('log2FC (T3 vs T2)')
ax.set_ylabel('Density')
ax.legend(fontsize=9)
ax.set_title(f'B) T2→T3 (after MTase drop)\np={p_t3t2:.2e}')

# Panel C: Effect size comparison
ax = axes[2]
transitions = ['T2 vs T1', 'T3 vs T2', 'T3 vs T1']
effects = []
cis = []

for trans, vals_lost, vals_never in [
    ('T2 vs T1', lfc_lost, lfc_never),
    ('T3 vs T2', lost_t3t2, never_t3t2),
    ('T3 vs T1', lost_t3t1, never_t3t1)
]:
    med_diff_val = vals_lost.median() - vals_never.median()
    effects.append(med_diff_val)
    # Quick bootstrap for CI
    boot_d = []
    for _ in range(N_BOOTSTRAP):
        s1 = np.random.choice(vals_lost.values, len(vals_lost), replace=True)
        s2 = np.random.choice(vals_never.values, len(vals_never), replace=True)
        boot_d.append(np.median(s1) - np.median(s2))
    cis.append((np.percentile(boot_d, 2.5), np.percentile(boot_d, 97.5)))

ci_low = [c[0] for c in cis]
ci_high = [c[1] for c in cis]
yerr = [np.array(effects) - np.array(ci_low), np.array(ci_high) - np.array(effects)]

ax.barh(transitions, effects, xerr=yerr, color=['#e74c3c', '#3498db', '#2ecc71'],
        alpha=0.7, edgecolor='black', linewidth=0.5, capsize=5)
ax.axvline(x=0, color='black', linestyle='--', alpha=0.5)
ax.set_xlabel('Median LFC difference (Lost - Never)')
ax.set_title('C) Effect size by transition\n(with 95% bootstrap CI)')

plt.tight_layout()
for fmt in ['pdf', 'svg']:
    fig.savefig(FIG_DIR / f'temporal_specificity.{fmt}', bbox_inches='tight')
plt.close()
print("  Saved: temporal_specificity.pdf/svg")


# ---------- Figure 6: Comprehensive summary ----------
fig = plt.figure(figsize=(20, 16))
gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.35)

# Panel A: H21 comparison bar chart
ax = fig.add_subplot(gs[0, 0])
gene_sets_short = ['All genes', 'All reg.', 'Shielded', 'Exposed']
p_values = h21_comp_df['MWU_p'].values
r_values = h21_comp_df['rank_biserial_r'].values
med_diffs = h21_comp_df['median_diff'].values

bar_colors = ['#bdc3c7', '#95a5a6', '#3498db', '#e74c3c']
bars = ax.bar(gene_sets_short, [-np.log10(p) if not np.isnan(p) and p > 0 else 0 for p in p_values],
              color=bar_colors, alpha=0.7, edgecolor='black', linewidth=0.5)
ax.axhline(y=-np.log10(0.05), color='red', linestyle='--', alpha=0.5, label='p=0.05')
ax.set_ylabel('-log10(p)')
ax.set_title('A) Signal strength by gene set')
ax.legend(fontsize=8)

# Panel B: Box plot Lost vs Never (exposed)
ax = fig.add_subplot(gs[0, 1])
if len(groups_to_plot) >= 2:
    bp = ax.boxplot(groups_to_plot[:2], positions=[1, 2], widths=0.5, patch_artist=True, showfliers=True,
                   flierprops={'markersize': 3, 'alpha': 0.4})
    for patch, color in zip(bp['boxes'], [COLORS['Lost'], COLORS['Never']]):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
    for med in bp['medians']:
        med.set_color('black')
        med.set_linewidth(2)
    for i, (vals, pos) in enumerate(zip(groups_to_plot[:2], [1, 2])):
        jitter = np.random.normal(0, 0.06, len(vals))
        ax.scatter(pos + jitter, vals, c=[COLORS['Lost'], COLORS['Never']][i], s=20, alpha=0.5, zorder=3)
    ax.set_xticks([1, 2])
    ax.set_xticklabels([labels_to_plot[0], labels_to_plot[1]])
ax.axhline(y=0, color='black', linestyle='--', alpha=0.3)
ax.set_ylabel('log2FC (T2 vs T1)')
ax.set_title(f'B) KEY TEST: Lost vs Never\np={p_val:.2e}, r={r_rb:.3f}')

# Panel C: Dose-response
ax = fig.add_subplot(gs[0, 2])
x = dose_data['n_AAGCCCG_lost'].values
y = dose_data['LFC_T2vsT1'].values
colors_dose = [COLORS['Lost'] if xi > 0 else COLORS['Never'] for xi in x]
ax.scatter(x, y, c=colors_dose, s=30, alpha=0.6, edgecolors='white', linewidth=0.5)
if len(x) > 2:
    z = np.polyfit(x, y, 1)
    ax.plot(np.linspace(x.min(), x.max(), 100), np.polyval(z, np.linspace(x.min(), x.max(), 100)),
            'k--', alpha=0.5, linewidth=1.5)
ax.set_xlabel('n AAGCCCG sites lost')
ax.set_ylabel('log2FC (T2 vs T1)')
ax.set_title(f'C) Dose-response\nrho={rho_dose:.3f}, p={p_dose:.2e}')

# Panel D: MTase correlation exposed vs shielded
ax = fig.add_subplot(gs[1, 0])
ax.hist(corr_shielded['spearman_rho'].dropna(), bins=30, alpha=0.4, color=COLORS['shielded'],
        label=f'Shielded (n={len(corr_shielded)})', density=True, edgecolor='white')
ax.hist(corr_exposed['spearman_rho'].dropna(), bins=15, alpha=0.7, color=COLORS['exposed'],
        label=f'Exposed (n={len(corr_exposed)})', density=True, edgecolor='white')
ax.axvline(corr_exposed['spearman_rho'].median(), color=COLORS['exposed'], linestyle='--', linewidth=2)
ax.axvline(corr_shielded['spearman_rho'].median(), color=COLORS['shielded'], linestyle='--', linewidth=2)
ax.set_xlabel('Spearman rho with SC_RS17645')
ax.legend(fontsize=8)
ax.set_title(f'D) SC_RS17645 correlation\np={p_corr2:.2e}')

# Panel E: Temporal specificity
ax = fig.add_subplot(gs[1, 1])
ax.barh(transitions, effects, xerr=yerr, color=['#e74c3c', '#3498db', '#2ecc71'],
        alpha=0.7, edgecolor='black', linewidth=0.5, capsize=5)
ax.axvline(x=0, color='black', linestyle='--', alpha=0.5)
ax.set_xlabel('Median LFC diff (Lost - Never)')
ax.set_title('E) Temporal specificity\n(95% bootstrap CI)')

# Panel F: Bootstrap distribution
ax = fig.add_subplot(gs[1, 2])
if 'median_diff_T2vsT1' in bootstrap_results:
    bd = bootstrap_results['median_diff_T2vsT1']
    # Regenerate bootstrap samples for histogram
    boot_hist = []
    for _ in range(N_BOOTSTRAP):
        s1 = np.random.choice(lfc_lost.values, len(lfc_lost), replace=True)
        s2 = np.random.choice(lfc_never.values, len(lfc_never), replace=True)
        boot_hist.append(np.median(s1) - np.median(s2))
    ax.hist(boot_hist, bins=50, alpha=0.7, color='#9b59b6', edgecolor='white', density=True)
    ax.axvline(x=0, color='red', linestyle='--', linewidth=2, label='Null (0)')
    ax.axvline(x=bd[0], color='black', linewidth=2, label=f'Observed: {bd[0]:.3f}')
    ax.axvline(x=bd[1], color='black', linestyle=':', linewidth=1)
    ax.axvline(x=bd[2], color='black', linestyle=':', linewidth=1)
    ax.set_xlabel('Median LFC diff (Lost - Never)')
    ax.set_ylabel('Density')
    ax.legend(fontsize=8)
    ax.set_title(f'F) Bootstrap (n={N_BOOTSTRAP})\n95% CI: [{bd[1]:.2f}, {bd[2]:.2f}]')

# Panel G: AAGCCCG site counts per exposed TF
ax = fig.add_subplot(gs[2, 0])
site_counts = tf_methyl['n_AAGCCCG_T1'].value_counts().sort_index()
ax.bar(site_counts.index, site_counts.values, color='#e74c3c', alpha=0.7, edgecolor='black', linewidth=0.5)
ax.set_xlabel('Number of AAGCCCG T1 sites (±500bp TSS)')
ax.set_ylabel('Number of exposed TFs')
ax.set_title(f'G) AAGCCCG site distribution\n(62 exposed TFs)')

# Panel H: GCCGGC comparison
ax = fig.add_subplot(gs[2, 1])
gcc_groups = []
gcc_labels_p = []
gcc_colors_p = []
for status, color in [('Lost', COLORS['Lost']), ('Never', COLORS['Never']), ('Gained', COLORS['Gained']), ('Both', COLORS['Both'])]:
    vals = tf_methyl[tf_methyl['GCCGGC_status'] == status]['LFC_T2vsT1'].dropna()
    if len(vals) >= 1:
        gcc_groups.append(vals.values)
        gcc_labels_p.append(f'{status}\n(n={len(vals)})')
        gcc_colors_p.append(color)

if len(gcc_groups) >= 2:
    bp_g = ax.boxplot(gcc_groups, positions=list(range(1, len(gcc_groups)+1)), widths=0.5,
                      patch_artist=True, showfliers=True, flierprops={'markersize': 3, 'alpha': 0.4})
    for patch, color in zip(bp_g['boxes'], gcc_colors_p):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
    for med in bp_g['medians']:
        med.set_color('black')
        med.set_linewidth(2)
    ax.set_xticks(list(range(1, len(gcc_groups)+1)))
    ax.set_xticklabels(gcc_labels_p, fontsize=9)
ax.axhline(y=0, color='black', linestyle='--', alpha=0.3)
ax.set_ylabel('log2FC (T2 vs T1)')
ax.set_title('H) GCCGGC 4mC comparison\n(Exposed TFs)')

# Panel I: Causal pathway model diagram
ax = fig.add_subplot(gs[2, 2])
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.axis('off')

# Draw the pathway
boxes = [
    (1, 8, 'SC_RS17645\n(HsdM MTase)', '#3498db'),
    (1, 5.5, 'AAGCCCG 6mA\nat exposed TF TSS', '#e74c3c'),
    (1, 3, '62 Exposed TFs\n(no protection zone)', '#f39c12'),
    (6, 8, 'T1: High expression\n(LFC=0)', '#3498db'),
    (6, 5.5, 'T1: Methylated (260 sites)\nT2: Lost (64 remain)', '#e74c3c'),
    (6, 3, 'T1→T2: Expression change\n(63% early responders)', '#f39c12'),
]

for x, y, text, color in boxes:
    bbox = FancyBboxPatch((x-0.8, y-0.6), 3.5, 1.2, boxstyle="round,pad=0.1",
                          facecolor=color, alpha=0.2, edgecolor=color, linewidth=1.5)
    ax.add_patch(bbox)
    ax.text(x+0.95, y, text, ha='center', va='center', fontsize=7, fontweight='bold')

# Arrows
arrow_style = dict(arrowstyle='->', color='black', lw=1.5)
ax.annotate('', xy=(1, 6.3), xytext=(1, 7.2), arrowprops=arrow_style)
ax.annotate('', xy=(1, 3.8), xytext=(1, 4.7), arrowprops=arrow_style)
ax.annotate('', xy=(5, 8), xytext=(3.5, 8), arrowprops=dict(arrowstyle='->', color='grey', lw=1))
ax.annotate('', xy=(5, 5.5), xytext=(3.5, 5.5), arrowprops=dict(arrowstyle='->', color='grey', lw=1))
ax.annotate('', xy=(5, 3), xytext=(3.5, 3), arrowprops=dict(arrowstyle='->', color='grey', lw=1))

ax.text(5, 1.5, f"H35 KEY RESULT:\nLost vs Never p={p_val:.2e}\nr={r_rb:.3f}",
        ha='center', va='center', fontsize=9, fontweight='bold',
        bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.3))
ax.set_title('I) Causal pathway model', fontsize=10)

fig.suptitle('H35: AAGCCCG Causal Pathway — SC_RS17645 → Exposed TFs', fontsize=14, y=0.98)
for fmt in ['pdf', 'svg']:
    fig.savefig(FIG_DIR / f'H35_comprehensive_summary.{fmt}', bbox_inches='tight')
plt.close()
print("  Saved: H35_comprehensive_summary.pdf/svg")


# ============================================================================
# Step 11: Save tables
# ============================================================================
print("\n" + "=" * 80)
print("STEP 11: Saving tables")
print("=" * 80)

# Table 1: Per-TF methylation status
tf_methyl.to_csv(TAB_DIR / 'exposed_TF_methylation_status.tsv', sep='\t', index=False)
print(f"  Saved: exposed_TF_methylation_status.tsv ({len(tf_methyl)} rows)")

# Table 2: De-repression test results
derep_rows = []
# AAGCCCG
for status in ['Lost', 'Never', 'Gained', 'Both']:
    sub = tf_methyl[tf_methyl['AAGCCCG_status'] == status]
    if len(sub) > 0:
        vals = sub['LFC_T2vsT1'].dropna()
        vals_t3 = sub['LFC_T3vsT1'].dropna()
        derep_rows.append({
            'motif': 'AAGCCCG', 'group': status, 'n': len(sub),
            'median_LFC_T2vsT1': vals.median() if len(vals) > 0 else np.nan,
            'mean_LFC_T2vsT1': vals.mean() if len(vals) > 0 else np.nan,
            'std_LFC_T2vsT1': vals.std() if len(vals) > 0 else np.nan,
            'median_LFC_T3vsT1': vals_t3.median() if len(vals_t3) > 0 else np.nan,
        })
# GCCGGC
for status in ['Lost', 'Never', 'Gained', 'Both']:
    sub = tf_methyl[tf_methyl['GCCGGC_status'] == status]
    if len(sub) > 0:
        vals = sub['LFC_T2vsT1'].dropna()
        vals_t3 = sub['LFC_T3vsT1'].dropna()
        derep_rows.append({
            'motif': 'GCCGGC', 'group': status, 'n': len(sub),
            'median_LFC_T2vsT1': vals.median() if len(vals) > 0 else np.nan,
            'mean_LFC_T2vsT1': vals.mean() if len(vals) > 0 else np.nan,
            'std_LFC_T2vsT1': vals.std() if len(vals) > 0 else np.nan,
            'median_LFC_T3vsT1': vals_t3.median() if len(vals_t3) > 0 else np.nan,
        })
pd.DataFrame(derep_rows).to_csv(TAB_DIR / 'derepression_test.tsv', sep='\t', index=False)
print(f"  Saved: derepression_test.tsv")

# Table 3: SC_RS17645 correlations
corr_save = corr_exposed.merge(tf_methyl[['locus_tag', 'AAGCCCG_status', 'n_AAGCCCG_lost',
                                           'LFC_T2vsT1', 'product', 'tf_family']], on='locus_tag', how='left')
corr_save.to_csv(TAB_DIR / 'SC_RS17645_correlations.tsv', sep='\t', index=False)
print(f"  Saved: SC_RS17645_correlations.tsv ({len(corr_save)} rows)")

# Table 4: All statistical tests
stat_df = pd.DataFrame(stat_tests)
stat_df.to_csv(TAB_DIR / 'statistical_tests.tsv', sep='\t', index=False)
print(f"  Saved: statistical_tests.tsv ({len(stat_df)} rows)")

# Table 5: H21 comparison
h21_comp_df.to_csv(TAB_DIR / 'H21_comparison.tsv', sep='\t', index=False)
print(f"  Saved: H21_comparison.tsv")


# ============================================================================
# Step 12: Summary
# ============================================================================
print("\n" + "=" * 80)
print("SUMMARY: H35 AAGCCCG Exposed TF Causal Pathway")
print("=" * 80)

print(f"\n1. AAGCCCG site mapping (±500bp TSS, 62 exposed TFs):")
for status in ['Lost', 'Never', 'Gained', 'Both']:
    n = len(tf_methyl[tf_methyl['AAGCCCG_status'] == status])
    print(f"   {status}: {n} ({n/62*100:.1f}%)")

print(f"\n2. KEY TEST: Exposed-restricted de-repression")
print(f"   Lost vs Never (LFC T2vsT1): p={p_val:.4e}, r={r_rb:.3f}, d={cohens_d:.3f}")
print(f"   Lost median: {lfc_lost.median():.3f}, Never median: {lfc_never.median():.3f}")
print(f"   Diff: {lfc_lost.median() - lfc_never.median():.3f}")
if 'median_diff_T2vsT1' in bootstrap_results:
    bd = bootstrap_results['median_diff_T2vsT1']
    print(f"   Bootstrap 95% CI for median diff: [{bd[1]:.3f}, {bd[2]:.3f}]")

print(f"\n3. H21 comparison (genome-wide → exposed-restricted):")
for _, row in h21_comp_df.iterrows():
    p_str = f"p={row['MWU_p']:.2e}" if not np.isnan(row['MWU_p']) else "p=N/A"
    r_str = f"r={row['rank_biserial_r']:.3f}" if not np.isnan(row['rank_biserial_r']) else "r=N/A"
    print(f"   {row['gene_set']}: {p_str}, {r_str} (n_lost={row['n_lost']}, n_never={row['n_never']})")

print(f"\n4. Dose-response: rho={rho_dose:.3f}, p={p_dose:.2e}")
print(f"\n5. SC_RS17645 correlation:")
print(f"   Exposed median rho: {corr_exposed['spearman_rho'].median():.3f}")
print(f"   Shielded median rho: {corr_shielded['spearman_rho'].median():.3f}")
print(f"   Difference p: {p_corr2:.4e}")

print(f"\n6. Temporal specificity:")
print(f"   T2vsT1 effect: median diff={lfc_lost.median() - lfc_never.median():.3f}, p={p_val:.2e}")
if len(lost_t3t2) >= 2 and len(never_t3t2) >= 2:
    print(f"   T3vsT2 effect: median diff={lost_t3t2.median() - never_t3t2.median():.3f}, p={p_t3t2:.2e}")

# Determine overall verdict
print("\n" + "=" * 80)
is_significant = p_val < 0.05
is_directional = lfc_lost.median() > lfc_never.median()  # de-repression = upregulation
is_exposed_specific = True  # check if exposed p < genome-wide p

gw_p = h21_comp_df[h21_comp_df['gene_set'].str.startswith('All genes')]['MWU_p'].values[0]
exp_p = h21_comp_df[h21_comp_df['gene_set'].str.startswith('Exposed')]['MWU_p'].values[0]

print(f"\nVERDICT:")
print(f"  Significant (p<0.05): {is_significant}")
print(f"  Direction (Lost > Never = de-repression): {is_directional}")
print(f"  Exposed-specific (exposed p < genome-wide p): {exp_p < gw_p if not np.isnan(exp_p) and not np.isnan(gw_p) else 'N/A'}")
if is_significant and is_directional:
    verdict = "SUPPORTED: AAGCCCG loss causes de-repression specifically in exposed TFs"
elif is_significant and not is_directional:
    verdict = "PARTIALLY SUPPORTED: Significant but unexpected direction"
elif not is_significant and abs(r_rb) > 0.2:
    verdict = "SUGGESTIVE: Non-significant but moderate effect size"
else:
    verdict = "NOT SUPPORTED: No evidence for exposed-specific AAGCCCG causal effect"
print(f"  {verdict}")
print("=" * 80)
