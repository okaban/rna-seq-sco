#!/usr/bin/env python3
"""
H33: Exposed TF Genomic Neighborhood Transcriptional Impact Analysis

Tests whether genes in the genomic neighborhood of 57 exposed TFs show
expression changes correlated with the TF's methylation coordination type.

Author: Claude Code
Date: 2026-02-27
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

# ============================================================
# Configuration
# ============================================================
BASE = Path("/Users/okaban/bioinfo/rna-seq")
ANALYSIS_DIR = BASE / "11_epigenome_integration/analysis/56_exposed_TF_neighborhood"
FIG_DIR = ANALYSIS_DIR / "figures"
TAB_DIR = ANALYSIS_DIR / "tables"

CHROM_LEN = 8_667_507  # NC_003888.3 linear chromosome

WINDOWS = [10_000, 20_000, 50_000]  # ±10kb, ±20kb, ±50kb

# Input files
GENE_ANNOT = BASE / "05_annotation/analysis/05_annotation_260128_v1/tables/gene_annotation_basic.tsv"
EXPOSED_TF = BASE / "11_epigenome_integration/analysis/51_exposed_regulators_characteristics/tables/exposed_regulators_full_table.tsv"
ALL_REG = BASE / "11_epigenome_integration/analysis/52_shielded_exposed_boundary/tables/all_genes_features.tsv"
DESEQ_T2T1 = BASE / "04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_2_vs_1.tsv"
DESEQ_T3T1 = BASE / "04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_3_vs_1.tsv"
DESEQ_T3T2 = BASE / "04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_3_vs_2.tsv"
COORD_REG = BASE / "11_epigenome_integration/analysis/29_genomewide_TF_screen/tables/coordinated_regulatory_genes.tsv"
COEXPR_MOD = BASE / "11_epigenome_integration/analysis/55_exposed_regulatory_module/tables/coexpression_modules.tsv"

# Plotting style
plt.rcParams.update({
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'font.family': 'sans-serif',
})

N_PERM = 1000
np.random.seed(42)


def save_fig(fig, name):
    """Save figure as both PDF and SVG."""
    fig.savefig(FIG_DIR / f"{name}.pdf")
    fig.savefig(FIG_DIR / f"{name}.svg")
    plt.close(fig)
    print(f"  Saved: {name}.pdf/.svg")


# ============================================================
# Step 1: Load data
# ============================================================
print("=" * 70)
print("STEP 1: Loading data")
print("=" * 70)

# Gene annotation (all 8,275 genes)
annot = pd.read_csv(GENE_ANNOT, sep='\t')
print(f"  Gene annotation: {len(annot)} genes")
# Compute midpoint for distance calculations
annot['midpoint'] = (annot['start'] + annot['end']) / 2

# Exposed TFs (57)
exposed = pd.read_csv(EXPOSED_TF, sep='\t')
exposed_ids = set(exposed['locus_tag'].tolist())
print(f"  Exposed TFs: {len(exposed)}")

# All regulatory genes (1,017)
all_reg = pd.read_csv(ALL_REG, sep='\t')
all_reg_ids = set(all_reg['locus_tag'].tolist())
shielded_reg = all_reg[all_reg['is_exposed'] == 0].copy()
shielded_ids = set(shielded_reg['locus_tag'].tolist())
print(f"  All regulatory genes: {len(all_reg)} (shielded: {len(shielded_ids)}, exposed: {sum(all_reg['is_exposed'])})")

# DESeq2 results
deseq_t2t1 = pd.read_csv(DESEQ_T2T1, sep='\t')
deseq_t3t1 = pd.read_csv(DESEQ_T3T1, sep='\t')
deseq_t3t2 = pd.read_csv(DESEQ_T3T2, sep='\t')
print(f"  DESeq2 T2vsT1: {len(deseq_t2t1)} genes")
print(f"  DESeq2 T3vsT1: {len(deseq_t3t1)} genes")
print(f"  DESeq2 T3vsT2: {len(deseq_t3t2)} genes")

# Coordinated regulatory genes
coord_reg = pd.read_csv(COORD_REG, sep='\t')
print(f"  Coordinated regulatory genes: {len(coord_reg)}")

# Co-expression modules
coexpr = pd.read_csv(COEXPR_MOD, sep='\t')
print(f"  Co-expression modules: {len(coexpr)}")

# Build module assignment dict
module_assignment = {}
for _, row in coexpr.iterrows():
    mod_id = row['module_id']
    tags = [t.strip() for t in row['locus_tags'].split(',')]
    for t in tags:
        module_assignment[t] = mod_id
print(f"  Module assignments: {len(module_assignment)} TFs across {coexpr['module_id'].nunique()} modules")

# Merge DESeq2 data into a single reference
# Use gene_id from DESeq2, which matches gene_id in annotation
expr_data = annot[['gene_id', 'start', 'end', 'strand', 'midpoint', 'product']].copy()

# Merge T2vsT1
t2t1_cols = deseq_t2t1[['gene_id', 'baseMean', 'log2FoldChange', 'padj']].rename(
    columns={'log2FoldChange': 'LFC_T2vsT1', 'padj': 'padj_T2vsT1', 'baseMean': 'baseMean'})
expr_data = expr_data.merge(t2t1_cols, on='gene_id', how='left')

# Merge T3vsT1
t3t1_cols = deseq_t3t1[['gene_id', 'log2FoldChange', 'padj']].rename(
    columns={'log2FoldChange': 'LFC_T3vsT1', 'padj': 'padj_T3vsT1'})
expr_data = expr_data.merge(t3t1_cols, on='gene_id', how='left')

# Merge T3vsT2
t3t2_cols = deseq_t3t2[['gene_id', 'log2FoldChange', 'padj']].rename(
    columns={'log2FoldChange': 'LFC_T3vsT2', 'padj': 'padj_T3vsT2'})
expr_data = expr_data.merge(t3t2_cols, on='gene_id', how='left')

# Compute absolute LFC values
expr_data['absLFC_T2vsT1'] = expr_data['LFC_T2vsT1'].abs()
expr_data['absLFC_T3vsT1'] = expr_data['LFC_T3vsT1'].abs()

# DEG flags (padj < 0.05)
expr_data['DEG_T2vsT1'] = expr_data['padj_T2vsT1'] < 0.05
expr_data['DEG_T3vsT1'] = expr_data['padj_T3vsT1'] < 0.05

print(f"  Expression data merged: {len(expr_data)} genes with expression info")
print(f"    With LFC T2vsT1: {expr_data['LFC_T2vsT1'].notna().sum()}")
print(f"    With LFC T3vsT1: {expr_data['LFC_T3vsT1'].notna().sum()}")

# Map locus_tag (exposed TF uses locus_tag) to gene_id (annotation uses gene_id)
# They should be the same format (SC_RS...)
# Verify
exposed_in_annot = set(exposed['locus_tag']) & set(annot['gene_id'])
print(f"  Exposed TFs found in annotation: {len(exposed_in_annot)}/57")

# Build exposed_coords from available columns (Jeong2016-filtered set lacks some old columns)
_col_map = {'log2FC_T2': 'LFC_T2vsT1', 'log2FC_T3': 'LFC_T3vsT1'}
for old_c, new_c in _col_map.items():
    if old_c not in exposed.columns and new_c in exposed.columns:
        exposed[old_c] = exposed[new_c]
for c in ['coordination_T2', 'coordination_T3', 'padj_T2', 'padj_T3', 'temporal_pattern']:
    if c not in exposed.columns:
        exposed[c] = np.nan
exposed_use_cols = ['locus_tag', 'start', 'end', 'strand', 'coordination_T2', 'coordination_T3',
                    'log2FC_T2', 'log2FC_T3', 'padj_T2', 'padj_T3', 'baseMean',
                    'product', 'tf_family', 'old_locus_tag', 'temporal_pattern']
exposed_coords = exposed[[c for c in exposed_use_cols if c in exposed.columns]].copy()
for c in exposed_use_cols:
    if c not in exposed_coords.columns:
        exposed_coords[c] = np.nan
exposed_coords['midpoint'] = (exposed_coords['start'] + exposed_coords['end']) / 2

# Get shielded TF coordinates from all_reg
shielded_coords = all_reg[all_reg['is_exposed'] == 0][
    ['locus_tag', 'start', 'end', 'strand', 'LFC_T2vsT1', 'LFC_T3vsT1']].copy()
shielded_coords['midpoint'] = (shielded_coords['start'] + shielded_coords['end']) / 2

# ============================================================
# Step 2: Define neighborhoods
# ============================================================
print("\n" + "=" * 70)
print("STEP 2: Defining neighborhoods")
print("=" * 70)


def get_neighbors(tf_id, tf_start, tf_end, window, expr_df, exclude_ids):
    """Get all genes within ±window of the TF, excluding certain gene IDs."""
    tf_mid = (tf_start + tf_end) / 2
    # Calculate distance from TF midpoint to gene midpoint
    mask = (
        (expr_df['midpoint'] >= tf_mid - window) &
        (expr_df['midpoint'] <= tf_mid + window) &
        (~expr_df['gene_id'].isin(exclude_ids))
    )
    # Handle chromosome edges
    if tf_mid - window < 0:
        mask = mask & (expr_df['midpoint'] >= 0)
    if tf_mid + window > CHROM_LEN:
        mask = mask & (expr_df['midpoint'] <= CHROM_LEN)

    neighbors = expr_df[mask].copy()
    neighbors['distance'] = (neighbors['midpoint'] - tf_mid).abs()
    neighbors['tf_id'] = tf_id
    return neighbors


# Build neighborhood for all windows
all_neighborhoods = {}
for w in WINDOWS:
    neighborhood_records = []
    for _, tf in exposed_coords.iterrows():
        tf_id = tf['locus_tag']
        nbrs = get_neighbors(tf_id, tf['start'], tf['end'], w, expr_data, exposed_ids)
        if len(nbrs) > 0:
            nbrs['window'] = w
            neighborhood_records.append(nbrs)
    if neighborhood_records:
        all_neighborhoods[w] = pd.concat(neighborhood_records, ignore_index=True)
    print(f"  Window ±{w//1000}kb: {len(all_neighborhoods[w])} neighbor-TF pairs, "
          f"mean {len(all_neighborhoods[w])/57:.1f} neighbors/TF")

# Also build shielded TF neighborhoods (for comparison)
shielded_neighborhoods = {}
for w in WINDOWS:
    neighborhood_records = []
    for _, tf in shielded_coords.iterrows():
        tf_id = tf['locus_tag']
        nbrs = get_neighbors(tf_id, tf['start'], tf['end'], w, expr_data, all_reg_ids)
        if len(nbrs) > 0:
            nbrs['window'] = w
            neighborhood_records.append(nbrs)
    if neighborhood_records:
        shielded_neighborhoods[w] = pd.concat(neighborhood_records, ignore_index=True)
    print(f"  Shielded ±{w//1000}kb: {len(shielded_neighborhoods[w])} neighbor-TF pairs, "
          f"mean {len(shielded_neighborhoods[w])/len(shielded_coords):.1f} neighbors/TF")

# For exposed neighborhoods, EXCLUDE other regulatory genes (exposed already excluded)
# The task says exclude exposed TFs only, but let's keep non-exposed regulatory genes
# as they may be legitimate targets

# ============================================================
# Step 3: Neighborhood expression analysis
# ============================================================
print("\n" + "=" * 70)
print("STEP 3: Neighborhood expression analysis")
print("=" * 70)

# Background: all non-TF genes
background_ids = set(expr_data['gene_id']) - all_reg_ids
background = expr_data[expr_data['gene_id'].isin(background_ids)].copy()
bg_absLFC_T2 = background['absLFC_T2vsT1'].dropna()
bg_absLFC_T3 = background['absLFC_T3vsT1'].dropna()
bg_deg_T2 = background['DEG_T2vsT1'].dropna().mean()
bg_deg_T3 = background['DEG_T3vsT1'].dropna().mean()
print(f"  Background (non-TF): {len(background)} genes, "
      f"mean |LFC| T2={bg_absLFC_T2.mean():.3f}, T3={bg_absLFC_T3.mean():.3f}")
print(f"  Background DEG rate: T2={bg_deg_T2:.3f}, T3={bg_deg_T3:.3f}")

stat_records = []

for w in WINDOWS:
    print(f"\n  --- Window ±{w//1000}kb ---")
    exp_nb = all_neighborhoods[w]
    shi_nb = shielded_neighborhoods[w]

    # Per-TF mean |LFC|
    exp_per_tf_t2 = exp_nb.groupby('tf_id')['absLFC_T2vsT1'].mean().dropna()
    exp_per_tf_t3 = exp_nb.groupby('tf_id')['absLFC_T3vsT1'].mean().dropna()
    shi_per_tf_t2 = shi_nb.groupby('tf_id')['absLFC_T2vsT1'].mean().dropna()
    shi_per_tf_t3 = shi_nb.groupby('tf_id')['absLFC_T3vsT1'].mean().dropna()

    # Per-TF DEG rate
    def deg_rate(grp, col):
        valid = grp[col].dropna()
        if len(valid) == 0:
            return np.nan
        return valid.mean()

    exp_deg_t2 = exp_nb.groupby('tf_id').apply(lambda g: deg_rate(g, 'DEG_T2vsT1')).dropna()
    exp_deg_t3 = exp_nb.groupby('tf_id').apply(lambda g: deg_rate(g, 'DEG_T3vsT1')).dropna()
    shi_deg_t2 = shi_nb.groupby('tf_id').apply(lambda g: deg_rate(g, 'DEG_T2vsT1')).dropna()
    shi_deg_t3 = shi_nb.groupby('tf_id').apply(lambda g: deg_rate(g, 'DEG_T3vsT1')).dropna()

    # Wilcoxon tests: exposed vs shielded
    stat_exp_shi_t2 = stats.mannwhitneyu(exp_per_tf_t2, shi_per_tf_t2, alternative='two-sided')
    stat_exp_shi_t3 = stats.mannwhitneyu(exp_per_tf_t3, shi_per_tf_t3, alternative='two-sided')

    # Wilcoxon tests: exposed vs background (use per-gene background, compare to pooled exposed neighbor |LFC|)
    # For a fair test, compare per-TF mean to background mean via permutation-like approach
    # Actually: compare the distribution of per-TF means to a single background value is not valid
    # Better: compare the pooled neighbor |LFC| to background |LFC|
    exp_pool_t2 = exp_nb['absLFC_T2vsT1'].dropna()
    exp_pool_t3 = exp_nb['absLFC_T3vsT1'].dropna()
    stat_exp_bg_t2 = stats.mannwhitneyu(exp_pool_t2, bg_absLFC_T2, alternative='two-sided')
    stat_exp_bg_t3 = stats.mannwhitneyu(exp_pool_t3, bg_absLFC_T3, alternative='two-sided')

    # DEG rate comparison
    stat_deg_exp_shi_t2 = stats.mannwhitneyu(exp_deg_t2, shi_deg_t2, alternative='two-sided')
    stat_deg_exp_shi_t3 = stats.mannwhitneyu(exp_deg_t3, shi_deg_t3, alternative='two-sided')

    print(f"    Exposed neighborhood mean |LFC| T2={exp_per_tf_t2.mean():.4f}, T3={exp_per_tf_t3.mean():.4f}")
    print(f"    Shielded neighborhood mean |LFC| T2={shi_per_tf_t2.mean():.4f}, T3={shi_per_tf_t3.mean():.4f}")
    print(f"    Exposed vs Shielded |LFC|: T2 p={stat_exp_shi_t2.pvalue:.4e}, T3 p={stat_exp_shi_t3.pvalue:.4e}")
    print(f"    Exposed vs Background |LFC|: T2 p={stat_exp_bg_t2.pvalue:.4e}, T3 p={stat_exp_bg_t3.pvalue:.4e}")
    print(f"    Exposed DEG rate T2={exp_deg_t2.mean():.4f}, T3={exp_deg_t3.mean():.4f}")
    print(f"    Shielded DEG rate T2={shi_deg_t2.mean():.4f}, T3={shi_deg_t3.mean():.4f}")
    print(f"    DEG rate exposed vs shielded: T2 p={stat_deg_exp_shi_t2.pvalue:.4e}, T3 p={stat_deg_exp_shi_t3.pvalue:.4e}")

    # Effect sizes (rank-biserial r = 1 - 2U/(n1*n2))
    n1, n2 = len(exp_per_tf_t2), len(shi_per_tf_t2)
    r_t2 = 1 - 2 * stat_exp_shi_t2.statistic / (n1 * n2) if n1 * n2 > 0 else np.nan
    n1, n2 = len(exp_per_tf_t3), len(shi_per_tf_t3)
    r_t3 = 1 - 2 * stat_exp_shi_t3.statistic / (n1 * n2) if n1 * n2 > 0 else np.nan

    for comp, metric, timepoint, stat_val, pval, eff in [
        ('Exposed_vs_Shielded', 'mean_absLFC', 'T2vsT1', stat_exp_shi_t2.statistic, stat_exp_shi_t2.pvalue, r_t2),
        ('Exposed_vs_Shielded', 'mean_absLFC', 'T3vsT1', stat_exp_shi_t3.statistic, stat_exp_shi_t3.pvalue, r_t3),
        ('Exposed_vs_Background', 'pooled_absLFC', 'T2vsT1', stat_exp_bg_t2.statistic, stat_exp_bg_t2.pvalue, np.nan),
        ('Exposed_vs_Background', 'pooled_absLFC', 'T3vsT1', stat_exp_bg_t3.statistic, stat_exp_bg_t3.pvalue, np.nan),
        ('Exposed_vs_Shielded', 'DEG_rate', 'T2vsT1', stat_deg_exp_shi_t2.statistic, stat_deg_exp_shi_t2.pvalue, np.nan),
        ('Exposed_vs_Shielded', 'DEG_rate', 'T3vsT1', stat_deg_exp_shi_t3.statistic, stat_deg_exp_shi_t3.pvalue, np.nan),
    ]:
        stat_records.append({
            'step': 'Step3_neighborhood_expression',
            'window_kb': w // 1000,
            'comparison': comp,
            'metric': metric,
            'timepoint': timepoint,
            'statistic': stat_val,
            'p_value': pval,
            'effect_size_r': eff,
            'exposed_mean': exp_per_tf_t2.mean() if timepoint == 'T2vsT1' and 'absLFC' in metric else (exp_per_tf_t3.mean() if timepoint == 'T3vsT1' and 'absLFC' in metric else (exp_deg_t2.mean() if timepoint == 'T2vsT1' else exp_deg_t3.mean())),
            'control_mean': shi_per_tf_t2.mean() if timepoint == 'T2vsT1' and comp == 'Exposed_vs_Shielded' and 'absLFC' in metric else (shi_per_tf_t3.mean() if timepoint == 'T3vsT1' and comp == 'Exposed_vs_Shielded' and 'absLFC' in metric else (bg_absLFC_T2.mean() if timepoint == 'T2vsT1' and 'absLFC' in metric else (bg_absLFC_T3.mean() if 'absLFC' in metric else np.nan))),
            'n_exposed': len(exp_per_tf_t2),
            'n_control': len(shi_per_tf_t2) if 'Shielded' in comp else len(bg_absLFC_T2),
        })


# ============================================================
# Step 4: Directional concordance
# ============================================================
print("\n" + "=" * 70)
print("STEP 4: Directional concordance analysis")
print("=" * 70)

# Use ±20kb as the primary window
w_primary = 20_000
exp_nb_20 = all_neighborhoods[w_primary]
shi_nb_20 = shielded_neighborhoods[w_primary]

# For each exposed TF, determine its direction at T3
concordance_records = []

for _, tf in exposed_coords.iterrows():
    tf_id = tf['locus_tag']
    tf_lfc_t3 = tf['log2FC_T3']
    tf_padj_t3 = tf['padj_T3']
    coord_t3 = tf['coordination_T3']

    # Determine TF direction at T3
    if pd.isna(tf_lfc_t3):
        tf_dir = 'unknown'
    elif tf_lfc_t3 > 0:
        tf_dir = 'up'
    else:
        tf_dir = 'down'

    # Get neighbors for this TF
    nbrs = exp_nb_20[exp_nb_20['tf_id'] == tf_id].copy()
    nbrs_with_lfc = nbrs.dropna(subset=['LFC_T3vsT1'])

    if len(nbrs_with_lfc) == 0:
        concordance_records.append({
            'tf_id': tf_id, 'tf_dir_T3': tf_dir, 'coordination_T3': coord_t3,
            'n_neighbors': 0, 'n_concordant': 0, 'concordance_rate': np.nan,
            'binom_p': np.nan
        })
        continue

    # Count concordant (same direction as TF)
    if tf_dir == 'up':
        n_concordant = (nbrs_with_lfc['LFC_T3vsT1'] > 0).sum()
    elif tf_dir == 'down':
        n_concordant = (nbrs_with_lfc['LFC_T3vsT1'] < 0).sum()
    else:
        n_concordant = np.nan

    n_total = len(nbrs_with_lfc)
    conc_rate = n_concordant / n_total if n_total > 0 and not pd.isna(n_concordant) else np.nan

    # Binomial test: is concordance > 50%?
    if not pd.isna(n_concordant) and n_total > 0:
        binom_p = stats.binomtest(int(n_concordant), n_total, 0.5, alternative='greater').pvalue
    else:
        binom_p = np.nan

    concordance_records.append({
        'tf_id': tf_id, 'tf_dir_T3': tf_dir, 'coordination_T3': coord_t3,
        'n_neighbors': n_total, 'n_concordant': int(n_concordant) if not pd.isna(n_concordant) else 0,
        'concordance_rate': conc_rate, 'binom_p': binom_p
    })

concordance_df = pd.DataFrame(concordance_records)

# Count significant concordances
sig_conc = concordance_df[(concordance_df['binom_p'] < 0.05) & (concordance_df['concordance_rate'] > 0.5)]
print(f"  TFs with significant directional concordance (p<0.05): {len(sig_conc)}/57")
print(f"  Overall mean concordance rate: {concordance_df['concordance_rate'].mean():.4f}")

# Group by coordination type
print("\n  Concordance by coordination type (T3):")
coord_types = concordance_df.groupby('coordination_T3').agg(
    n_TFs=('tf_id', 'count'),
    mean_concordance=('concordance_rate', 'mean'),
    median_concordance=('concordance_rate', 'median'),
    n_sig=('binom_p', lambda x: (x < 0.05).sum())
).reset_index()
for _, row in coord_types.iterrows():
    print(f"    {row['coordination_T3']}: n={row['n_TFs']}, "
          f"mean conc={row['mean_concordance']:.3f}, median={row['median_concordance']:.3f}, "
          f"n_sig={row['n_sig']}")

# Compare concordance: exposed vs shielded TFs
# Compute concordance for shielded TFs at T3
shielded_concordance = []
for _, tf in shielded_coords.iterrows():
    tf_id = tf['locus_tag']
    tf_lfc_t3 = tf['LFC_T3vsT1']
    if pd.isna(tf_lfc_t3):
        continue
    tf_dir = 'up' if tf_lfc_t3 > 0 else 'down'

    nbrs = shi_nb_20[shi_nb_20['tf_id'] == tf_id].copy()
    nbrs_with_lfc = nbrs.dropna(subset=['LFC_T3vsT1'])
    if len(nbrs_with_lfc) == 0:
        continue

    if tf_dir == 'up':
        n_concordant = (nbrs_with_lfc['LFC_T3vsT1'] > 0).sum()
    else:
        n_concordant = (nbrs_with_lfc['LFC_T3vsT1'] < 0).sum()

    n_total = len(nbrs_with_lfc)
    conc_rate = n_concordant / n_total
    shielded_concordance.append(conc_rate)

exposed_conc_rates = concordance_df['concordance_rate'].dropna()
shielded_conc_rates = pd.Series(shielded_concordance).dropna()

stat_conc_exp_shi = stats.mannwhitneyu(exposed_conc_rates, shielded_conc_rates, alternative='two-sided')
print(f"\n  Concordance comparison (±20kb):")
print(f"    Exposed mean: {exposed_conc_rates.mean():.4f} (n={len(exposed_conc_rates)})")
print(f"    Shielded mean: {shielded_conc_rates.mean():.4f} (n={len(shielded_conc_rates)})")
print(f"    Mann-Whitney U p={stat_conc_exp_shi.pvalue:.4e}")

stat_records.append({
    'step': 'Step4_directional_concordance',
    'window_kb': 20,
    'comparison': 'Exposed_vs_Shielded_concordance',
    'metric': 'concordance_rate_T3',
    'timepoint': 'T3vsT1',
    'statistic': stat_conc_exp_shi.statistic,
    'p_value': stat_conc_exp_shi.pvalue,
    'effect_size_r': 1 - 2 * stat_conc_exp_shi.statistic / (len(exposed_conc_rates) * len(shielded_conc_rates)),
    'exposed_mean': exposed_conc_rates.mean(),
    'control_mean': shielded_conc_rates.mean(),
    'n_exposed': len(exposed_conc_rates),
    'n_control': len(shielded_conc_rates),
})

# Concordance by coordination type: for "up" TFs - are neighbors more up?
# For "down" TFs - are neighbors more down?
print("\n  Coordination type directional analysis:")
for coord_type in concordance_df['coordination_T3'].unique():
    sub = concordance_df[concordance_df['coordination_T3'] == coord_type]
    mean_c = sub['concordance_rate'].mean()
    n_up = (sub['tf_dir_T3'] == 'up').sum()
    n_down = (sub['tf_dir_T3'] == 'down').sum()
    print(f"    {coord_type}: n_up={n_up}, n_down={n_down}, mean_concordance={mean_c:.3f}")

# ============================================================
# Step 5: Activation vs Repression bloc neighborhoods
# ============================================================
print("\n" + "=" * 70)
print("STEP 5: Activation vs Repression bloc neighborhoods")
print("=" * 70)

# Module 1-3 = activation bloc, Module 4 = repression bloc
activation_tfs = set()
repression_tfs = set()
for tf_id, mod in module_assignment.items():
    if mod in [1, 2, 3]:
        activation_tfs.add(tf_id)
    elif mod == 4:
        repression_tfs.add(tf_id)

# Only consider exposed TFs
activation_exposed = activation_tfs & exposed_ids
repression_exposed = repression_tfs & exposed_ids
unassigned = exposed_ids - activation_tfs - repression_tfs
print(f"  Activation bloc (M1-3): {len(activation_exposed)} exposed TFs")
print(f"  Repression bloc (M4): {len(repression_exposed)} exposed TFs")
print(f"  Unassigned: {len(unassigned)} exposed TFs")

# Get neighbors for each bloc using ±20kb
act_nb = exp_nb_20[exp_nb_20['tf_id'].isin(activation_exposed)].copy()
rep_nb = exp_nb_20[exp_nb_20['tf_id'].isin(repression_exposed)].copy()

print(f"  Activation bloc neighbors: {len(act_nb)}")
print(f"  Repression bloc neighbors: {len(rep_nb)}")

# Keyword enrichment in product annotation
keywords = ['secondary metabolite', 'transport', 'signal', 'kinase', 'synthase',
            'oxidoreductase', 'transferase', 'hydrolase', 'dehydrogenase',
            'ABC', 'MFS', 'permease', 'polyketide', 'peptide synthetase',
            'cytochrome', 'regulator', 'hypothetical', 'membrane', 'secreted',
            'ribosom', 'protease', 'reductase', 'methyltransferase']

bloc_comparison_records = []

for kw in keywords:
    act_count = act_nb['product'].str.contains(kw, case=False, na=False).sum()
    rep_count = rep_nb['product'].str.contains(kw, case=False, na=False).sum()
    act_frac = act_count / len(act_nb) if len(act_nb) > 0 else 0
    rep_frac = rep_count / len(rep_nb) if len(rep_nb) > 0 else 0

    # Fisher's exact test
    act_yes = act_count
    act_no = len(act_nb) - act_count
    rep_yes = rep_count
    rep_no = len(rep_nb) - rep_count
    if act_yes + rep_yes > 0:
        odds, fisher_p = stats.fisher_exact([[act_yes, act_no], [rep_yes, rep_no]])
    else:
        odds, fisher_p = np.nan, np.nan

    bloc_comparison_records.append({
        'keyword': kw,
        'activation_count': act_count,
        'activation_fraction': act_frac,
        'repression_count': rep_count,
        'repression_fraction': rep_frac,
        'odds_ratio': odds,
        'fisher_p': fisher_p
    })

bloc_df = pd.DataFrame(bloc_comparison_records)
sig_kw = bloc_df[bloc_df['fisher_p'] < 0.05]
print(f"\n  Significant keyword differences (p<0.05): {len(sig_kw)}")
for _, row in sig_kw.iterrows():
    print(f"    {row['keyword']}: act={row['activation_fraction']:.3f}, rep={row['repression_fraction']:.3f}, "
          f"OR={row['odds_ratio']:.2f}, p={row['fisher_p']:.4f}")

# Mean baseMean of neighbors
act_basemean = act_nb['baseMean'].dropna()
rep_basemean = rep_nb['baseMean'].dropna()
stat_basemean = stats.mannwhitneyu(act_basemean, rep_basemean, alternative='two-sided')
print(f"\n  BaseMean comparison:")
print(f"    Activation neighbors: median={act_basemean.median():.1f} (n={len(act_basemean)})")
print(f"    Repression neighbors: median={rep_basemean.median():.1f} (n={len(rep_basemean)})")
print(f"    Mann-Whitney p={stat_basemean.pvalue:.4e}")

# DEG rate comparison
act_deg_t3 = act_nb['DEG_T3vsT1'].dropna().mean()
rep_deg_t3 = rep_nb['DEG_T3vsT1'].dropna().mean()
# Per-TF DEG rate
act_per_tf_deg = act_nb.groupby('tf_id')['DEG_T3vsT1'].mean().dropna()
rep_per_tf_deg = rep_nb.groupby('tf_id')['DEG_T3vsT1'].mean().dropna()
if len(act_per_tf_deg) > 0 and len(rep_per_tf_deg) > 0:
    stat_deg_bloc = stats.mannwhitneyu(act_per_tf_deg, rep_per_tf_deg, alternative='two-sided')
else:
    stat_deg_bloc = type('obj', (object,), {'pvalue': np.nan, 'statistic': np.nan})()
print(f"\n  DEG rate (T3vsT1):")
print(f"    Activation: {act_deg_t3:.4f}")
print(f"    Repression: {rep_deg_t3:.4f}")
print(f"    Per-TF comparison p={stat_deg_bloc.pvalue:.4e}")

# Mean |LFC| comparison
act_lfc = act_nb['absLFC_T3vsT1'].dropna()
rep_lfc = rep_nb['absLFC_T3vsT1'].dropna()
stat_lfc_bloc = stats.mannwhitneyu(act_lfc, rep_lfc, alternative='two-sided')
print(f"\n  Mean |LFC| T3vsT1:")
print(f"    Activation: {act_lfc.mean():.4f} (median={act_lfc.median():.4f})")
print(f"    Repression: {rep_lfc.mean():.4f} (median={rep_lfc.median():.4f})")
print(f"    Mann-Whitney p={stat_lfc_bloc.pvalue:.4e}")

# Direction of neighbor change
act_up = (act_nb['LFC_T3vsT1'] > 0).sum()
act_down = (act_nb['LFC_T3vsT1'] < 0).sum()
rep_up = (rep_nb['LFC_T3vsT1'] > 0).sum()
rep_down = (rep_nb['LFC_T3vsT1'] < 0).sum()
print(f"\n  Neighbor direction (T3vsT1):")
print(f"    Activation: {act_up} up, {act_down} down ({act_up/(act_up+act_down)*100:.1f}% up)")
print(f"    Repression: {rep_up} up, {rep_down} down ({rep_up/(rep_up+rep_down)*100:.1f}% up)")
try:
    chi2_dir, p_dir = stats.chi2_contingency([[act_up, act_down], [rep_up, rep_down]])[:2]
    print(f"    Chi-squared test: chi2={chi2_dir:.2f}, p={p_dir:.4e}")
except ValueError:
    chi2_dir, p_dir = np.nan, np.nan
    print("    Chi-squared test: skipped (zero cell in contingency table)")

for metric, pval in [('baseMean', stat_basemean.pvalue),
                      ('DEG_rate_T3', stat_deg_bloc.pvalue),
                      ('absLFC_T3', stat_lfc_bloc.pvalue),
                      ('direction_chi2', p_dir)]:
    stat_records.append({
        'step': 'Step5_bloc_comparison',
        'window_kb': 20,
        'comparison': 'Activation_vs_Repression',
        'metric': metric,
        'timepoint': 'T3vsT1',
        'statistic': np.nan,
        'p_value': pval,
        'effect_size_r': np.nan,
        'exposed_mean': np.nan,
        'control_mean': np.nan,
        'n_exposed': len(activation_exposed),
        'n_control': len(repression_exposed),
    })


# ============================================================
# Step 6: BGC proximity
# ============================================================
print("\n" + "=" * 70)
print("STEP 6: BGC proximity analysis")
print("=" * 70)

# BGC-related keywords
bgc_keywords = ['polyketide', 'non-ribosomal peptide', 'NRPS', 'PKS',
                 'terpene', 'siderophore', 'lantipeptide', 'lanthipeptide',
                 'melanin', 'ectoine', 'butyrolactone', 'bacteriocin',
                 'desferrioxamine', 'actinorhodin', 'undecylprodigiosin',
                 'prodiginine', 'coelimycin', 'calcium-dependent antibiotic',
                 'CDA', 'geosmin', 'methylisoborneol', 'hopanoid',
                 'chalcone', 'type I PKS', 'type II PKS', 'type III PKS']

# Check exposed TF products for BGC-related terms
print("  Checking exposed TF products for BGC-related terms:")
bgc_tfs = []
for _, tf in exposed_coords.iterrows():
    product = str(tf['product']).lower()
    for kw in bgc_keywords:
        if kw.lower() in product:
            bgc_tfs.append((tf['locus_tag'], tf['old_locus_tag'], tf['product'], kw))
            break

if bgc_tfs:
    for tf_id, old_tag, prod, kw in bgc_tfs:
        print(f"    {tf_id} ({old_tag}): {prod} [matched: {kw}]")
else:
    print("    No exposed TFs have BGC-related product annotations")

# Check neighbors for BGC-related products
print("\n  Checking ±50kb neighbors for BGC-related products:")
bgc_neighbor_count = 0
for _, tf in exposed_coords.iterrows():
    tf_id = tf['locus_tag']
    nbrs = all_neighborhoods[50_000]
    tf_nbrs = nbrs[nbrs['tf_id'] == tf_id]
    for _, nbr in tf_nbrs.iterrows():
        product = str(nbr['product']).lower()
        for kw in bgc_keywords:
            if kw.lower() in product:
                bgc_neighbor_count += 1
                break

print(f"    BGC-related neighbors within ±50kb: {bgc_neighbor_count}")

# More broadly, check for secondary metabolite pathway genes
sec_met_keywords = ['synthase', 'synthetase', 'cyclase', 'epimerase',
                     'acyltransferase', 'thioesterase', 'condensation',
                     'phosphopantetheinyl', 'PPTase']
sec_met_neighbors = 0
for _, tf in exposed_coords.iterrows():
    tf_id = tf['locus_tag']
    nbrs = all_neighborhoods[50_000]
    tf_nbrs = nbrs[nbrs['tf_id'] == tf_id]
    for _, nbr in tf_nbrs.iterrows():
        product = str(nbr['product']).lower()
        for kw in sec_met_keywords:
            if kw.lower() in product:
                sec_met_neighbors += 1
                break
print(f"    Secondary metabolism-related neighbors: {sec_met_neighbors}")

# Check each TF specifically
print("\n  Exposed TFs potentially near BGC clusters (based on neighbor products):")
bgc_proximal_tfs = []
for _, tf in exposed_coords.iterrows():
    tf_id = tf['locus_tag']
    nbrs = all_neighborhoods[20_000]
    tf_nbrs = nbrs[nbrs['tf_id'] == tf_id]
    bgc_nbr_count = 0
    bgc_products = []
    for _, nbr in tf_nbrs.iterrows():
        product = str(nbr['product']).lower()
        for kw in bgc_keywords + sec_met_keywords:
            if kw.lower() in product:
                bgc_nbr_count += 1
                bgc_products.append(nbr['product'])
                break
    if bgc_nbr_count >= 2:  # At least 2 BGC-related neighbors
        bgc_proximal_tfs.append({
            'tf_id': tf_id,
            'old_locus_tag': tf['old_locus_tag'],
            'product': tf['product'],
            'n_bgc_neighbors': bgc_nbr_count,
            'bgc_products': '; '.join(bgc_products[:5])
        })

if bgc_proximal_tfs:
    for tf_info in bgc_proximal_tfs:
        print(f"    {tf_info['tf_id']} ({tf_info['old_locus_tag']}): "
              f"{tf_info['n_bgc_neighbors']} BGC neighbors - {tf_info['bgc_products'][:80]}")
else:
    print("    No exposed TFs with >=2 BGC-related neighbors within ±20kb")


# ============================================================
# Step 7: Distance decay
# ============================================================
print("\n" + "=" * 70)
print("STEP 7: Distance decay analysis")
print("=" * 70)

# Use ±50kb window to capture full distance range
exp_nb_50 = all_neighborhoods[50_000]
shi_nb_50 = shielded_neighborhoods[50_000]

# Bin distances
distance_bins = np.arange(0, 51_000, 5_000)  # 0-5kb, 5-10kb, ..., 45-50kb
bin_labels = [f"{b//1000}-{(b+5000)//1000}kb" for b in distance_bins[:-1]]

# Exposed: per-bin mean |LFC|
exp_nb_50_valid = exp_nb_50.dropna(subset=['absLFC_T3vsT1']).copy()
exp_nb_50_valid['dist_bin'] = pd.cut(exp_nb_50_valid['distance'], bins=distance_bins, labels=bin_labels, right=False)
exp_dist_decay = exp_nb_50_valid.groupby('dist_bin', observed=True).agg(
    mean_absLFC_T3=('absLFC_T3vsT1', 'mean'),
    median_absLFC_T3=('absLFC_T3vsT1', 'median'),
    mean_absLFC_T2=('absLFC_T2vsT1', 'mean'),
    n_genes=('gene_id', 'count')
).reset_index()

# Shielded: per-bin mean |LFC|
shi_nb_50_valid = shi_nb_50.dropna(subset=['absLFC_T3vsT1']).copy()
shi_nb_50_valid['dist_bin'] = pd.cut(shi_nb_50_valid['distance'], bins=distance_bins, labels=bin_labels, right=False)
shi_dist_decay = shi_nb_50_valid.groupby('dist_bin', observed=True).agg(
    mean_absLFC_T3=('absLFC_T3vsT1', 'mean'),
    median_absLFC_T3=('absLFC_T3vsT1', 'median'),
    mean_absLFC_T2=('absLFC_T2vsT1', 'mean'),
    n_genes=('gene_id', 'count')
).reset_index()

print("  Distance decay (exposed TF neighborhoods, T3vsT1):")
for _, row in exp_dist_decay.iterrows():
    print(f"    {row['dist_bin']}: mean |LFC|={row['mean_absLFC_T3']:.4f}, n={row['n_genes']}")

# Per-TF Spearman correlation: distance vs |LFC|
per_tf_corr = []
for tf_id in exposed_coords['locus_tag']:
    tf_nbrs = exp_nb_50_valid[exp_nb_50_valid['tf_id'] == tf_id]
    if len(tf_nbrs) >= 5:
        rho, p = stats.spearmanr(tf_nbrs['distance'], tf_nbrs['absLFC_T3vsT1'])
        per_tf_corr.append({'tf_id': tf_id, 'rho': rho, 'p': p, 'n': len(tf_nbrs)})

per_tf_corr_df = pd.DataFrame(per_tf_corr)
mean_rho = per_tf_corr_df['rho'].mean()
# One-sample t-test: is mean rho different from 0?
t_stat, t_p = stats.ttest_1samp(per_tf_corr_df['rho'].dropna(), 0)
print(f"\n  Per-TF distance-|LFC| correlation:")
print(f"    Mean Spearman rho: {mean_rho:.4f}")
print(f"    t-test vs 0: t={t_stat:.3f}, p={t_p:.4e}")
print(f"    TFs with negative rho (closer = higher |LFC|): {(per_tf_corr_df['rho'] < 0).sum()}/{len(per_tf_corr_df)}")

# Same for shielded
shi_per_tf_corr = []
for tf_id in shielded_coords['locus_tag']:
    tf_nbrs = shi_nb_50_valid[shi_nb_50_valid['tf_id'] == tf_id]
    if len(tf_nbrs) >= 5:
        rho, p = stats.spearmanr(tf_nbrs['distance'], tf_nbrs['absLFC_T3vsT1'])
        shi_per_tf_corr.append({'tf_id': tf_id, 'rho': rho, 'p': p, 'n': len(tf_nbrs)})

shi_per_tf_corr_df = pd.DataFrame(shi_per_tf_corr)
shi_mean_rho = shi_per_tf_corr_df['rho'].mean()
shi_t, shi_p = stats.ttest_1samp(shi_per_tf_corr_df['rho'].dropna(), 0)
print(f"\n  Shielded TF distance-|LFC| correlation:")
print(f"    Mean Spearman rho: {shi_mean_rho:.4f}")
print(f"    t-test vs 0: t={shi_t:.3f}, p={shi_p:.4e}")

# Compare exposed vs shielded rho
stat_rho = stats.mannwhitneyu(per_tf_corr_df['rho'].dropna(), shi_per_tf_corr_df['rho'].dropna(), alternative='two-sided')
print(f"    Exposed vs Shielded rho: MWU p={stat_rho.pvalue:.4e}")

stat_records.append({
    'step': 'Step7_distance_decay',
    'window_kb': 50,
    'comparison': 'Exposed_distance_rho_vs_0',
    'metric': 'Spearman_rho_distance_absLFC',
    'timepoint': 'T3vsT1',
    'statistic': t_stat,
    'p_value': t_p,
    'effect_size_r': mean_rho,
    'exposed_mean': mean_rho,
    'control_mean': 0,
    'n_exposed': len(per_tf_corr_df),
    'n_control': np.nan,
})
stat_records.append({
    'step': 'Step7_distance_decay',
    'window_kb': 50,
    'comparison': 'Exposed_vs_Shielded_distance_rho',
    'metric': 'Spearman_rho_distance_absLFC',
    'timepoint': 'T3vsT1',
    'statistic': stat_rho.statistic,
    'p_value': stat_rho.pvalue,
    'effect_size_r': np.nan,
    'exposed_mean': mean_rho,
    'control_mean': shi_mean_rho,
    'n_exposed': len(per_tf_corr_df),
    'n_control': len(shi_per_tf_corr_df),
})


# ============================================================
# Step 8: Permutation test
# ============================================================
print("\n" + "=" * 70)
print("STEP 8: Permutation test (1,000 iterations)")
print("=" * 70)

# For the main finding: neighborhood mean |LFC| at ±20kb
# Observe: mean per-TF |LFC| across 57 exposed TFs
exp_per_tf_t3_20 = exp_nb_20.groupby('tf_id')['absLFC_T3vsT1'].mean().dropna()
observed_mean_absLFC = exp_per_tf_t3_20.mean()

# Also observe concordance rate
observed_mean_concordance = concordance_df['concordance_rate'].dropna().mean()

all_gene_ids = expr_data['gene_id'].tolist()
all_gene_coords = expr_data[['gene_id', 'start', 'end', 'midpoint']].copy()

perm_absLFC = []
perm_concordance = []
perm_deg_rate = []

for perm_i in range(N_PERM):
    # Random sample of 57 genes
    sample_ids = np.random.choice(all_gene_ids, size=57, replace=False)
    sample_coords = all_gene_coords[all_gene_coords['gene_id'].isin(sample_ids)]

    perm_records = []
    for _, g in sample_coords.iterrows():
        g_id = g['gene_id']
        nbrs = get_neighbors(g_id, g['start'], g['end'], w_primary, expr_data, set(sample_ids))
        if len(nbrs) > 0:
            perm_records.append(nbrs)

    if not perm_records:
        perm_absLFC.append(np.nan)
        perm_concordance.append(np.nan)
        perm_deg_rate.append(np.nan)
        continue

    perm_nb = pd.concat(perm_records, ignore_index=True)

    # Mean |LFC|
    per_gene_lfc = perm_nb.groupby('tf_id')['absLFC_T3vsT1'].mean().dropna()
    perm_absLFC.append(per_gene_lfc.mean())

    # Concordance rate
    conc_rates = []
    for g_id in sample_ids:
        g_row = expr_data[expr_data['gene_id'] == g_id]
        if len(g_row) == 0:
            continue
        g_lfc = g_row['LFC_T3vsT1'].values[0]
        if pd.isna(g_lfc):
            continue
        g_dir = 'up' if g_lfc > 0 else 'down'
        g_nbrs = perm_nb[perm_nb['tf_id'] == g_id]
        g_nbrs_valid = g_nbrs.dropna(subset=['LFC_T3vsT1'])
        if len(g_nbrs_valid) == 0:
            continue
        if g_dir == 'up':
            n_conc = (g_nbrs_valid['LFC_T3vsT1'] > 0).sum()
        else:
            n_conc = (g_nbrs_valid['LFC_T3vsT1'] < 0).sum()
        conc_rates.append(n_conc / len(g_nbrs_valid))
    perm_concordance.append(np.mean(conc_rates) if conc_rates else np.nan)

    # DEG rate
    per_gene_deg = perm_nb.groupby('tf_id')['DEG_T3vsT1'].mean().dropna()
    perm_deg_rate.append(per_gene_deg.mean())

    if (perm_i + 1) % 200 == 0:
        print(f"  Completed {perm_i + 1}/{N_PERM} permutations")

perm_absLFC = np.array(perm_absLFC)
perm_concordance = np.array(perm_concordance)
perm_deg_rate = np.array(perm_deg_rate)

# Empirical p-values
p_absLFC = (np.nansum(perm_absLFC >= observed_mean_absLFC) + 1) / (np.sum(~np.isnan(perm_absLFC)) + 1)
p_concordance = (np.nansum(perm_concordance >= observed_mean_concordance) + 1) / (np.sum(~np.isnan(perm_concordance)) + 1)

observed_deg_rate = exp_nb_20.groupby('tf_id')['DEG_T3vsT1'].mean().dropna().mean()
p_deg = (np.nansum(perm_deg_rate >= observed_deg_rate) + 1) / (np.sum(~np.isnan(perm_deg_rate)) + 1)

print(f"\n  Permutation results (±20kb, T3vsT1):")
print(f"    Mean |LFC|: observed={observed_mean_absLFC:.4f}, perm mean={np.nanmean(perm_absLFC):.4f}, "
      f"perm SD={np.nanstd(perm_absLFC):.4f}, empirical p={p_absLFC:.4f}")
print(f"    Concordance: observed={observed_mean_concordance:.4f}, perm mean={np.nanmean(perm_concordance):.4f}, "
      f"perm SD={np.nanstd(perm_concordance):.4f}, empirical p={p_concordance:.4f}")
print(f"    DEG rate: observed={observed_deg_rate:.4f}, perm mean={np.nanmean(perm_deg_rate):.4f}, "
      f"perm SD={np.nanstd(perm_deg_rate):.4f}, empirical p={p_deg:.4f}")

z_absLFC = (observed_mean_absLFC - np.nanmean(perm_absLFC)) / np.nanstd(perm_absLFC)
z_concordance = (observed_mean_concordance - np.nanmean(perm_concordance)) / np.nanstd(perm_concordance)

print(f"    Z-score |LFC|: {z_absLFC:.3f}")
print(f"    Z-score concordance: {z_concordance:.3f}")

stat_records.append({
    'step': 'Step8_permutation',
    'window_kb': 20,
    'comparison': 'Observed_vs_Permutation',
    'metric': 'mean_absLFC_T3',
    'timepoint': 'T3vsT1',
    'statistic': z_absLFC,
    'p_value': p_absLFC,
    'effect_size_r': np.nan,
    'exposed_mean': observed_mean_absLFC,
    'control_mean': np.nanmean(perm_absLFC),
    'n_exposed': 57,
    'n_control': N_PERM,
})
stat_records.append({
    'step': 'Step8_permutation',
    'window_kb': 20,
    'comparison': 'Observed_vs_Permutation',
    'metric': 'concordance_rate',
    'timepoint': 'T3vsT1',
    'statistic': z_concordance,
    'p_value': p_concordance,
    'effect_size_r': np.nan,
    'exposed_mean': observed_mean_concordance,
    'control_mean': np.nanmean(perm_concordance),
    'n_exposed': 57,
    'n_control': N_PERM,
})

# ============================================================
# Step 9: Figures
# ============================================================
print("\n" + "=" * 70)
print("STEP 9: Generating figures")
print("=" * 70)

# ------ Figure 1: Neighborhood LFC comparison ------
fig, axes = plt.subplots(2, 3, figsize=(14, 9))
fig.suptitle('H33: Neighborhood Expression Magnitude — Exposed vs Shielded vs Background', fontsize=14, fontweight='bold')

for i, w in enumerate(WINDOWS):
    # T2vsT1
    ax = axes[0, i]
    exp_vals = all_neighborhoods[w].groupby('tf_id')['absLFC_T2vsT1'].mean().dropna()
    shi_vals = shielded_neighborhoods[w].groupby('tf_id')['absLFC_T2vsT1'].mean().dropna()

    # Sample background to 200 random genes for visualization
    bg_sample = bg_absLFC_T2.sample(min(200, len(bg_absLFC_T2)), random_state=42)

    data = [exp_vals.values, shi_vals.values, bg_sample.values]
    parts = ax.violinplot(data, positions=[1, 2, 3], showmeans=True, showmedians=True)
    for pc, color in zip(parts['bodies'], ['#e74c3c', '#3498db', '#95a5a6']):
        pc.set_facecolor(color)
        pc.set_alpha(0.7)
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(['Exposed\n(n=57)', f'Shielded\n(n={len(shi_vals)})', f'Background\n(n={len(bg_sample)})'])
    ax.set_ylabel('Mean |log2FC| per TF')
    ax.set_title(f'T2vsT1 (±{w//1000}kb)')

    # T3vsT1
    ax = axes[1, i]
    exp_vals = all_neighborhoods[w].groupby('tf_id')['absLFC_T3vsT1'].mean().dropna()
    shi_vals = shielded_neighborhoods[w].groupby('tf_id')['absLFC_T3vsT1'].mean().dropna()
    bg_sample = bg_absLFC_T3.sample(min(200, len(bg_absLFC_T3)), random_state=42)

    data = [exp_vals.values, shi_vals.values, bg_sample.values]
    parts = ax.violinplot(data, positions=[1, 2, 3], showmeans=True, showmedians=True)
    for pc, color in zip(parts['bodies'], ['#e74c3c', '#3498db', '#95a5a6']):
        pc.set_facecolor(color)
        pc.set_alpha(0.7)
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(['Exposed\n(n=57)', f'Shielded\n(n={len(shi_vals)})', f'Background\n(n={len(bg_sample)})'])
    ax.set_ylabel('Mean |log2FC| per TF')
    ax.set_title(f'T3vsT1 (±{w//1000}kb)')

plt.tight_layout()
save_fig(fig, 'neighborhood_LFC_comparison')

# ------ Figure 2: Directional concordance ------
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle('H33: Directional Concordance of TF Neighborhoods', fontsize=14, fontweight='bold')

# Panel A: By coordination type
ax = axes[0]
coord_order = concordance_df.groupby('coordination_T3')['concordance_rate'].mean().sort_values(ascending=False).index
conc_by_type = concordance_df.groupby('coordination_T3')['concordance_rate'].agg(['mean', 'sem', 'count']).reindex(coord_order)
colors = ['#e74c3c' if m > 0.5 else '#3498db' for m in conc_by_type['mean']]
bars = ax.bar(range(len(conc_by_type)), conc_by_type['mean'], yerr=conc_by_type['sem'],
              color=colors, alpha=0.8, capsize=3)
ax.axhline(0.5, color='gray', linestyle='--', linewidth=1, label='Random (50%)')
ax.set_xticks(range(len(conc_by_type)))
ax.set_xticklabels([f"{ct}\n(n={int(n)})" for ct, n in zip(conc_by_type.index, conc_by_type['count'])],
                    rotation=30, ha='right', fontsize=8)
ax.set_ylabel('Concordance rate')
ax.set_title('By coordination type (T3)')
ax.legend()

# Panel B: Exposed vs Shielded
ax = axes[1]
exp_mean = exposed_conc_rates.mean()
shi_mean = shielded_conc_rates.mean()
exp_sem = exposed_conc_rates.sem()
shi_sem = shielded_conc_rates.sem()
ax.bar([0, 1], [exp_mean, shi_mean], yerr=[exp_sem, shi_sem],
       color=['#e74c3c', '#3498db'], alpha=0.8, capsize=5)
ax.axhline(0.5, color='gray', linestyle='--', linewidth=1)
ax.set_xticks([0, 1])
ax.set_xticklabels([f'Exposed\n(n={len(exposed_conc_rates)})',
                     f'Shielded\n(n={len(shielded_conc_rates)})'])
ax.set_ylabel('Mean concordance rate')
ax.set_title(f'Exposed vs Shielded (p={stat_conc_exp_shi.pvalue:.3e})')

plt.tight_layout()
save_fig(fig, 'directional_concordance')

# ------ Figure 3: Activation vs Repression neighborhoods ------
if len(rep_lfc) > 0 and len(act_lfc) > 0:
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    fig.suptitle('H33: Activation Bloc vs Repression Bloc Neighbor Characteristics', fontsize=14, fontweight='bold')

    # Panel A: |LFC| distribution
    ax = axes[0]
    data = [act_lfc.values, rep_lfc.values]
    parts = ax.violinplot(data, positions=[1, 2], showmeans=True, showmedians=True)
    for pc, color in zip(parts['bodies'], ['#e67e22', '#8e44ad']):
        pc.set_facecolor(color)
        pc.set_alpha(0.7)
    ax.set_xticks([1, 2])
    ax.set_xticklabels([f'Activation\n(M1-3, {len(activation_exposed)} TFs)',
                         f'Repression\n(M4, {len(repression_exposed)} TFs)'])
    ax.set_ylabel('|log2FC| T3vsT1')
    ax.set_title(f'Neighbor |LFC| (p={stat_lfc_bloc.pvalue:.3e})')

    # Panel B: DEG rate
    ax = axes[1]
    ax.bar([0, 1], [act_deg_t3, rep_deg_t3],
           color=['#e67e22', '#8e44ad'], alpha=0.8)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Activation', 'Repression'])
    ax.set_ylabel('DEG rate (T3vsT1, padj<0.05)')
    ax.set_title(f'Neighbor DEG rate (p={stat_deg_bloc.pvalue:.3e})')

    # Panel C: Direction pie charts
    ax = axes[2]
    sizes_act = [act_up, act_down]
    sizes_rep = [rep_up, rep_down]
    colors_dir = ['#27ae60', '#c0392b']
    ax.pie(sizes_act, labels=['Up', 'Down'], colors=colors_dir, autopct='%1.1f%%',
           startangle=90, center=(0.25, 0.5), radius=0.35)
    ax.pie(sizes_rep, labels=['Up', 'Down'], colors=colors_dir, autopct='%1.1f%%',
           startangle=90, center=(0.75, 0.5), radius=0.35)
    ax.text(0.25, 0.95, 'Activation', ha='center', fontsize=11, fontweight='bold', transform=ax.transAxes)
    ax.text(0.75, 0.95, 'Repression', ha='center', fontsize=11, fontweight='bold', transform=ax.transAxes)
    ax.set_title(f'Neighbor direction T3vsT1 (chi2 p={p_dir:.3e})')
    ax.set_xlim(-0.15, 1.15)
    ax.set_ylim(-0.05, 1.05)
    ax.set_aspect('equal')

    plt.tight_layout()
    save_fig(fig, 'activation_vs_repression_neighborhoods')
else:
    print("  Bloc comparison figure skipped (repression bloc has no neighbors with n=15)")

# ------ Figure 4: Distance decay ------
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle('H33: Distance Decay of Expression Change Near TFs', fontsize=14, fontweight='bold')

# Panel A: T3vsT1
ax = axes[0]
x_vals = np.arange(len(exp_dist_decay))
ax.plot(x_vals, exp_dist_decay['mean_absLFC_T3'].values, 'o-', color='#e74c3c',
        label=f'Exposed (n=57)', linewidth=2, markersize=6)
ax.plot(x_vals, shi_dist_decay['mean_absLFC_T3'].values, 's-', color='#3498db',
        label=f'Shielded (n={len(shielded_coords)})', linewidth=2, markersize=6)
ax.axhline(bg_absLFC_T3.mean(), color='gray', linestyle='--', label='Background')
ax.set_xticks(x_vals)
ax.set_xticklabels(bin_labels, rotation=45, ha='right', fontsize=8)
ax.set_xlabel('Distance from TF')
ax.set_ylabel('Mean |log2FC| T3vsT1')
ax.set_title('T3vsT1')
ax.legend(fontsize=9)

# Panel B: T2vsT1
ax = axes[1]
ax.plot(x_vals, exp_dist_decay['mean_absLFC_T2'].values, 'o-', color='#e74c3c',
        label='Exposed', linewidth=2, markersize=6)
ax.plot(x_vals, shi_dist_decay['mean_absLFC_T2'].values, 's-', color='#3498db',
        label='Shielded', linewidth=2, markersize=6)
ax.axhline(bg_absLFC_T2.mean(), color='gray', linestyle='--', label='Background')
ax.set_xticks(x_vals)
ax.set_xticklabels(bin_labels, rotation=45, ha='right', fontsize=8)
ax.set_xlabel('Distance from TF')
ax.set_ylabel('Mean |log2FC| T2vsT1')
ax.set_title('T2vsT1')
ax.legend(fontsize=9)

plt.tight_layout()
save_fig(fig, 'distance_decay')

# ------ Figure 5: Chromosome neighborhood map ------
fig, ax = plt.subplots(figsize=(16, 8))
fig.suptitle('H33: Chromosome Map of 57 Exposed TF Neighborhoods', fontsize=14, fontweight='bold')

# Draw chromosome as horizontal bar
chrom_y = 0.5
ax.barh(chrom_y, CHROM_LEN, height=0.08, color='#bdc3c7', edgecolor='black', linewidth=0.5)

# Add TF positions
for idx, (_, tf) in enumerate(exposed_coords.iterrows()):
    tf_mid = tf['midpoint']
    tf_lfc = tf['log2FC_T3']
    mod_id = module_assignment.get(tf['locus_tag'], 0)

    # Color by module
    mod_colors = {1: '#e74c3c', 2: '#e67e22', 3: '#f1c40f', 4: '#8e44ad', 0: '#95a5a6'}
    color = mod_colors.get(mod_id, '#95a5a6')

    # Y position: alternate above/below
    y_offset = 0.15 + (idx % 4) * 0.12 if idx % 2 == 0 else -(0.15 + (idx % 4) * 0.12)
    y_pos = chrom_y + y_offset

    # Draw line from chromosome to marker
    ax.plot([tf_mid, tf_mid], [chrom_y, y_pos], color=color, linewidth=0.5, alpha=0.5)

    # Arrow direction based on LFC
    marker = '^' if tf_lfc > 0 else 'v'
    size = min(max(abs(tf_lfc) * 20, 15), 80)
    ax.scatter(tf_mid, y_pos, marker=marker, s=size, c=color, edgecolors='black', linewidth=0.5, zorder=5)

    # Neighborhood shading (±20kb)
    rect = plt.Rectangle((tf_mid - 20_000, chrom_y - 0.035), 40_000, 0.07,
                          alpha=0.15, color=color, linewidth=0)
    ax.add_patch(rect)

# Legend
legend_elements = [
    mpatches.Patch(facecolor='#e74c3c', label='Module 1'),
    mpatches.Patch(facecolor='#e67e22', label='Module 2'),
    mpatches.Patch(facecolor='#f1c40f', label='Module 3'),
    mpatches.Patch(facecolor='#8e44ad', label='Module 4'),
    mpatches.Patch(facecolor='#95a5a6', label='Unassigned'),
    plt.Line2D([0], [0], marker='^', color='w', markerfacecolor='gray', markersize=8, label='Up at T3'),
    plt.Line2D([0], [0], marker='v', color='w', markerfacecolor='gray', markersize=8, label='Down at T3'),
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=9, ncol=2)

ax.set_xlim(-100_000, CHROM_LEN + 100_000)
ax.set_ylim(-0.5, 1.5)
ax.set_xlabel('Chromosome position (bp)')
ax.set_yticks([])

# Add scale bar
ax.plot([0, 1_000_000], [-0.35, -0.35], 'k-', linewidth=2)
ax.text(500_000, -0.42, '1 Mb', ha='center', fontsize=9)

plt.tight_layout()
save_fig(fig, 'chromosome_neighborhood_map')

# ------ Figure 6: Comprehensive summary ------
fig = plt.figure(figsize=(18, 22))
gs = GridSpec(4, 3, figure=fig, hspace=0.35, wspace=0.3)
fig.suptitle('H33: Comprehensive Summary — Exposed TF Neighborhood Impact', fontsize=16, fontweight='bold', y=0.98)

# Panel A: Neighborhood |LFC| comparison (±20kb, T3vsT1)
ax = fig.add_subplot(gs[0, 0])
exp_vals_20_t3 = all_neighborhoods[20_000].groupby('tf_id')['absLFC_T3vsT1'].mean().dropna()
shi_vals_20_t3 = shielded_neighborhoods[20_000].groupby('tf_id')['absLFC_T3vsT1'].mean().dropna()
bp = ax.boxplot([exp_vals_20_t3.values, shi_vals_20_t3.values, bg_absLFC_T3.sample(200, random_state=42).values],
                labels=['Exposed', 'Shielded', 'Background'],
                patch_artist=True, widths=0.6)
colors_box = ['#e74c3c', '#3498db', '#95a5a6']
for patch, color in zip(bp['boxes'], colors_box):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
ax.set_ylabel('Mean |log2FC| T3vsT1')
ax.set_title('A. Neighborhood |LFC| (±20kb)', fontweight='bold')

# Panel B: Concordance by coordination type
ax = fig.add_subplot(gs[0, 1])
conc_by_type_reset = conc_by_type.reset_index()
colors_coord = ['#e74c3c' if m > 0.55 else '#f39c12' if m > 0.5 else '#3498db'
                for m in conc_by_type_reset['mean']]
bars = ax.barh(range(len(conc_by_type_reset)), conc_by_type_reset['mean'],
               xerr=conc_by_type_reset['sem'], color=colors_coord, alpha=0.8, capsize=3)
ax.axvline(0.5, color='gray', linestyle='--', linewidth=1)
ax.set_yticks(range(len(conc_by_type_reset)))
ax.set_yticklabels(conc_by_type_reset['coordination_T3'], fontsize=8)
ax.set_xlabel('Concordance rate')
ax.set_title('B. Directional concordance', fontweight='bold')

# Panel C: Permutation distribution
ax = fig.add_subplot(gs[0, 2])
ax.hist(perm_absLFC[~np.isnan(perm_absLFC)], bins=30, color='#95a5a6', alpha=0.7, edgecolor='black', linewidth=0.5)
ax.axvline(observed_mean_absLFC, color='#e74c3c', linewidth=2, label=f'Observed ({observed_mean_absLFC:.3f})')
ax.axvline(np.nanmean(perm_absLFC), color='gray', linestyle='--', label=f'Perm mean ({np.nanmean(perm_absLFC):.3f})')
ax.set_xlabel('Mean |log2FC| T3vsT1')
ax.set_ylabel('Count')
ax.set_title(f'C. Permutation test (p={p_absLFC:.3f})', fontweight='bold')
ax.legend(fontsize=8)

# Panel D: Distance decay
ax = fig.add_subplot(gs[1, 0])
x_vals = np.arange(len(exp_dist_decay))
ax.plot(x_vals, exp_dist_decay['mean_absLFC_T3'].values, 'o-', color='#e74c3c', label='Exposed', linewidth=2)
ax.plot(x_vals, shi_dist_decay['mean_absLFC_T3'].values, 's-', color='#3498db', label='Shielded', linewidth=2)
ax.axhline(bg_absLFC_T3.mean(), color='gray', linestyle='--')
ax.set_xticks(x_vals[::2])
ax.set_xticklabels([bin_labels[i] for i in range(0, len(bin_labels), 2)], rotation=45, ha='right', fontsize=7)
ax.set_ylabel('Mean |log2FC| T3vsT1')
ax.set_title('D. Distance decay', fontweight='bold')
ax.legend(fontsize=8)

# Panel E: Bloc comparison |LFC|
ax = fig.add_subplot(gs[1, 1])
bp2 = ax.boxplot([act_lfc.values, rep_lfc.values], labels=['Activation\n(M1-3)', 'Repression\n(M4)'],
                  patch_artist=True, widths=0.5)
bp2['boxes'][0].set_facecolor('#e67e22')
bp2['boxes'][1].set_facecolor('#8e44ad')
for b in bp2['boxes']:
    b.set_alpha(0.7)
ax.set_ylabel('|log2FC| T3vsT1')
ax.set_title(f'E. Bloc neighbors (p={stat_lfc_bloc.pvalue:.3e})', fontweight='bold')

# Panel F: Bloc neighbor direction
ax = fig.add_subplot(gs[1, 2])
x = np.arange(2)
w_bar = 0.35
act_total = act_up + act_down
rep_total = rep_up + rep_down
ax.bar(x - w_bar/2, [act_up/act_total*100, rep_up/rep_total*100], w_bar,
       color='#27ae60', alpha=0.8, label='Up')
ax.bar(x + w_bar/2, [act_down/act_total*100, rep_down/rep_total*100], w_bar,
       color='#c0392b', alpha=0.8, label='Down')
ax.set_xticks(x)
ax.set_xticklabels(['Activation', 'Repression'])
ax.set_ylabel('% neighbors')
ax.set_title(f'F. Neighbor direction (chi2 p={p_dir:.3e})', fontweight='bold')
ax.legend(fontsize=8)

# Panel G: Per-TF concordance histogram
ax = fig.add_subplot(gs[2, 0])
ax.hist(concordance_df['concordance_rate'].dropna(), bins=20, color='#e74c3c', alpha=0.7,
        edgecolor='black', linewidth=0.5, label='Exposed')
ax.hist(shielded_conc_rates, bins=20, color='#3498db', alpha=0.5,
        edgecolor='black', linewidth=0.5, label='Shielded')
ax.axvline(0.5, color='gray', linestyle='--')
ax.set_xlabel('Concordance rate')
ax.set_ylabel('Count')
ax.set_title('G. Concordance distribution', fontweight='bold')
ax.legend(fontsize=8)

# Panel H: Keyword enrichment
ax = fig.add_subplot(gs[2, 1:])
top_kw = bloc_df.nsmallest(12, 'fisher_p')
x = np.arange(len(top_kw))
w_bar = 0.35
ax.barh(x - w_bar/2, top_kw['activation_fraction'] * 100, w_bar, color='#e67e22', alpha=0.8, label='Activation')
ax.barh(x + w_bar/2, top_kw['repression_fraction'] * 100, w_bar, color='#8e44ad', alpha=0.8, label='Repression')
ax.set_yticks(x)
ax.set_yticklabels(top_kw['keyword'], fontsize=8)
ax.set_xlabel('% neighbors with keyword')
ax.set_title('H. Product keyword enrichment (top 12)', fontweight='bold')
ax.legend(fontsize=8)
# Mark significant
for i, (_, row) in enumerate(top_kw.iterrows()):
    if row['fisher_p'] < 0.05:
        ax.text(max(row['activation_fraction'], row['repression_fraction']) * 100 + 0.5,
                i, '*', fontsize=14, color='red', va='center')

# Panel I: Chromosome map (simplified)
ax = fig.add_subplot(gs[3, :])
ax.barh(0.5, CHROM_LEN, height=0.15, color='#ecf0f1', edgecolor='black', linewidth=0.5)
for _, tf in exposed_coords.iterrows():
    tf_mid = tf['midpoint']
    mod_id = module_assignment.get(tf['locus_tag'], 0)
    mod_colors = {1: '#e74c3c', 2: '#e67e22', 3: '#f1c40f', 4: '#8e44ad', 0: '#95a5a6'}
    color = mod_colors.get(mod_id, '#95a5a6')
    marker = '^' if tf['log2FC_T3'] > 0 else 'v'
    y = 0.85 if tf['log2FC_T3'] > 0 else 0.15
    ax.scatter(tf_mid, y, marker=marker, s=40, c=color, edgecolors='black', linewidth=0.3, zorder=5)
    ax.plot([tf_mid, tf_mid], [0.425, y], color=color, linewidth=0.3, alpha=0.5)
ax.set_xlim(-100_000, CHROM_LEN + 100_000)
ax.set_ylim(-0.1, 1.1)
ax.set_xlabel('Chromosome position (bp)')
ax.set_yticks([])
ax.set_title('I. 57 Exposed TFs on NC_003888.3 (color = module, direction = T3 LFC)', fontweight='bold')

save_fig(fig, 'H33_comprehensive_summary')


# ============================================================
# Step 10: Tables
# ============================================================
print("\n" + "=" * 70)
print("STEP 10: Saving tables")
print("=" * 70)

# Table 1: Per-TF neighborhood stats
per_tf_records = []
for _, tf in exposed_coords.iterrows():
    tf_id = tf['locus_tag']
    rec = {
        'tf_id': tf_id,
        'old_locus_tag': tf['old_locus_tag'],
        'product': tf['product'],
        'tf_family': tf['tf_family'],
        'module': module_assignment.get(tf_id, np.nan),
        'coordination_T3': tf['coordination_T3'],
        'LFC_T2': tf['log2FC_T2'],
        'LFC_T3': tf['log2FC_T3'],
        'tf_start': tf['start'],
        'tf_end': tf['end'],
        'tf_strand': tf['strand'],
    }
    for w in WINDOWS:
        nb = all_neighborhoods[w]
        tf_nb = nb[nb['tf_id'] == tf_id]
        rec[f'n_neighbors_{w//1000}kb'] = len(tf_nb)
        rec[f'mean_absLFC_T2_{w//1000}kb'] = tf_nb['absLFC_T2vsT1'].mean()
        rec[f'mean_absLFC_T3_{w//1000}kb'] = tf_nb['absLFC_T3vsT1'].mean()
        valid = tf_nb['DEG_T3vsT1'].dropna()
        rec[f'DEG_rate_T3_{w//1000}kb'] = valid.mean() if len(valid) > 0 else np.nan

    # Concordance (from concordance_df)
    conc_row = concordance_df[concordance_df['tf_id'] == tf_id]
    if len(conc_row) > 0:
        rec['concordance_rate'] = conc_row.iloc[0]['concordance_rate']
        rec['concordance_binom_p'] = conc_row.iloc[0]['binom_p']
    else:
        rec['concordance_rate'] = np.nan
        rec['concordance_binom_p'] = np.nan

    # Distance correlation
    corr_row = per_tf_corr_df[per_tf_corr_df['tf_id'] == tf_id]
    if len(corr_row) > 0:
        rec['distance_rho'] = corr_row.iloc[0]['rho']
        rec['distance_rho_p'] = corr_row.iloc[0]['p']
    else:
        rec['distance_rho'] = np.nan
        rec['distance_rho_p'] = np.nan

    per_tf_records.append(rec)

per_tf_df = pd.DataFrame(per_tf_records)
per_tf_df.to_csv(TAB_DIR / 'per_TF_neighborhood.tsv', sep='\t', index=False)
print(f"  Saved: per_TF_neighborhood.tsv ({len(per_tf_df)} TFs)")

# Table 2: Neighborhood genes (all pairs)
all_nb_genes = []
for w in [20_000]:  # Primary window
    nb = all_neighborhoods[w]
    for _, row in nb.iterrows():
        all_nb_genes.append({
            'tf_id': row['tf_id'],
            'neighbor_gene_id': row['gene_id'],
            'distance': row['distance'],
            'neighbor_product': row['product'],
            'LFC_T2vsT1': row['LFC_T2vsT1'],
            'LFC_T3vsT1': row['LFC_T3vsT1'],
            'padj_T2vsT1': row['padj_T2vsT1'],
            'padj_T3vsT1': row['padj_T3vsT1'],
            'baseMean': row['baseMean'],
            'DEG_T2vsT1': row['DEG_T2vsT1'],
            'DEG_T3vsT1': row['DEG_T3vsT1'],
            'window_kb': w // 1000,
        })

nb_genes_df = pd.DataFrame(all_nb_genes)
nb_genes_df.to_csv(TAB_DIR / 'neighborhood_genes.tsv', sep='\t', index=False)
print(f"  Saved: neighborhood_genes.tsv ({len(nb_genes_df)} pairs)")

# Table 3: Bloc comparison
bloc_df.to_csv(TAB_DIR / 'bloc_comparison.tsv', sep='\t', index=False)
print(f"  Saved: bloc_comparison.tsv ({len(bloc_df)} keywords)")

# Table 4: Statistical tests
stat_df = pd.DataFrame(stat_records)
stat_df.to_csv(TAB_DIR / 'statistical_tests.tsv', sep='\t', index=False)
print(f"  Saved: statistical_tests.tsv ({len(stat_df)} tests)")

# Table 5: Permutation results
perm_results = pd.DataFrame({
    'metric': ['mean_absLFC_T3', 'concordance_rate', 'DEG_rate_T3'],
    'observed': [observed_mean_absLFC, observed_mean_concordance, observed_deg_rate],
    'perm_mean': [np.nanmean(perm_absLFC), np.nanmean(perm_concordance), np.nanmean(perm_deg_rate)],
    'perm_sd': [np.nanstd(perm_absLFC), np.nanstd(perm_concordance), np.nanstd(perm_deg_rate)],
    'z_score': [z_absLFC, z_concordance,
                (observed_deg_rate - np.nanmean(perm_deg_rate)) / np.nanstd(perm_deg_rate)],
    'empirical_p': [p_absLFC, p_concordance, p_deg],
    'n_permutations': [N_PERM, N_PERM, N_PERM],
    'window_kb': [20, 20, 20],
})
perm_results.to_csv(TAB_DIR / 'permutation_results.tsv', sep='\t', index=False)
print(f"  Saved: permutation_results.tsv")


# ============================================================
# Final summary
# ============================================================
print("\n" + "=" * 70)
print("H33 ANALYSIS COMPLETE")
print("=" * 70)

print(f"\nKey findings:")
print(f"  1. Neighborhood |LFC| (±20kb, T3vsT1):")
print(f"     Exposed: {exp_per_tf_t3_20.mean():.4f}, Shielded: {shi_vals_20_t3.mean():.4f}")
w20_stat = [r for r in stat_records if r['step'] == 'Step3_neighborhood_expression'
            and r['window_kb'] == 20 and r['metric'] == 'mean_absLFC'
            and r['timepoint'] == 'T3vsT1' and r['comparison'] == 'Exposed_vs_Shielded']
if w20_stat:
    print(f"     p={w20_stat[0]['p_value']:.4e}, r={w20_stat[0]['effect_size_r']:.3f}")

print(f"  2. Directional concordance:")
print(f"     Exposed: {exposed_conc_rates.mean():.4f}, Shielded: {shielded_conc_rates.mean():.4f}")
print(f"     p={stat_conc_exp_shi.pvalue:.4e}")

print(f"  3. Permutation test (mean |LFC|):")
print(f"     Z={z_absLFC:.3f}, empirical p={p_absLFC:.4f}")

print(f"  4. Distance decay:")
print(f"     Exposed mean rho={mean_rho:.4f} (p={t_p:.4e})")
print(f"     Shielded mean rho={shi_mean_rho:.4f}")

print(f"  5. Activation vs Repression:")
print(f"     |LFC|: act={act_lfc.mean():.4f} vs rep={rep_lfc.mean():.4f} (p={stat_lfc_bloc.pvalue:.4e})")
print(f"     Direction: chi2 p={p_dir:.4e}")

print(f"\nOutput files:")
print(f"  Figures: {FIG_DIR}")
print(f"  Tables: {TAB_DIR}")
print(f"\nDone!")
