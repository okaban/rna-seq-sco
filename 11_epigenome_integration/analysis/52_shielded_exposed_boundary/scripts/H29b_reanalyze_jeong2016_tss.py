#!/usr/bin/env python3
"""
H29b: Re-analysis of Shielded/Exposed boundary using Jeong2016 experimental TSSs only.

Motivation:
  The original H29 analysis used GFF annotation gene-start positions as TSS proxies
  for all 1,055 regulatory genes. However, only 450 of these have experimentally
  validated TSSs from Jeong et al. (2016) dRNA-seq.  Using unannotated gene starts
  as TSSs introduces noise in nearest_methyl_distance, inflating the AUC and sample
  size in Fig 2c.

Changes vs H29:
  - Filter to 449 regulatory genes with Jeong2016 dRNA-seq TSS
  - Recompute nearest_methyl_distance from the correct experimental TSS
  - Re-run ROC / Mann-Whitney analysis
  - Overwrite all_genes_features.tsv (old version archived)

Outputs (overwrite):
  tables/all_genes_features.tsv        — 449-gene feature matrix
  tables/ROC_analysis.tsv              — updated ROC results
"""

import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve
import warnings
warnings.filterwarnings('ignore')

BASE = Path('/Users/okaban/bioinfo/rna-seq')
EPIGENOME = BASE / '11_epigenome_integration' / 'analysis'
TBL_DIR = EPIGENOME / '52_shielded_exposed_boundary' / 'tables'

# ── Archive old tables ────────────────────────────────────────────────────────
import shutil, datetime
tag = '260312'
for fname in ['all_genes_features.tsv', 'ROC_analysis.tsv']:
    src = TBL_DIR / fname
    if src.exists():
        dst = TBL_DIR / f'{src.stem}_ARCHIVED_{tag}.tsv'
        if not dst.exists():
            shutil.copy(src, dst)
            print(f'Archived: {dst.name}')

# ── 1. Load inputs ────────────────────────────────────────────────────────────
print('Loading data...')

# Jeong2016 experimental TSS (protein-coding regulatory genes only)
tss_src = pd.read_csv(EPIGENOME / '18_tss_analyses' / 'comprehensive_tss_table.csv')
jeong_tss = (tss_src[tss_src['tss_source'] == 'Jeong2016_dRNA-seq']
             [['gene_id', 'tss', 'strand']]
             .rename(columns={'gene_id': 'locus_tag', 'tss': 'exp_tss'}))
print(f'Jeong2016 TSS entries: {len(jeong_tss)}')

# Previous feature matrix (for is_exposed, LFC, baseMean, etc.)
old_feat = pd.read_csv(TBL_DIR / f'all_genes_features_ARCHIVED_{tag}.tsv', sep='\t')
print(f'Old feature matrix: {len(old_feat)} genes '
      f'(exposed={old_feat["is_exposed"].sum()}, shielded={(old_feat["is_exposed"]==0).sum()})')

# HC methylation sites (union across all timepoints — same as original analysis)
methyl_all = pd.read_csv(
    EPIGENOME / '01_integration' / 'high_confidence_sites_weighted.csv')
methyl_pos = np.sort(methyl_all['position'].unique())
print(f'HC methylation sites: {len(methyl_pos)} unique positions')

# ── 2. Build Jeong2016-filtered feature matrix ────────────────────────────────
print('\nBuilding filtered feature matrix...')

df = old_feat.merge(jeong_tss[['locus_tag', 'exp_tss']], on='locus_tag', how='inner')
print(f'After Jeong2016 filter: {len(df)} genes '
      f'(exposed={df["is_exposed"].sum()}, shielded={(df["is_exposed"]==0).sum()})')

# Replace TSS with experimental value
df['tss_old'] = df['tss']
df['tss'] = df['exp_tss']
df.drop(columns=['exp_tss'], inplace=True)

# ── 3. Recompute nearest_methyl_distance using experimental TSS ───────────────
print('Recomputing nearest_methyl_distance...')

def nearest_methyl_distance(tss, sorted_positions):
    idx = np.searchsorted(sorted_positions, tss)
    candidates = []
    if idx > 0:
        candidates.append(abs(tss - sorted_positions[idx - 1]))
    if idx < len(sorted_positions):
        candidates.append(abs(tss - sorted_positions[idx]))
    return min(candidates) if candidates else np.nan

def count_sites_within(tss, sorted_positions, window=2000):
    lo = np.searchsorted(sorted_positions, tss - window, side='left')
    hi = np.searchsorted(sorted_positions, tss + window, side='right')
    return hi - lo

df['nearest_methyl_distance'] = df['tss'].apply(
    lambda t: nearest_methyl_distance(t, methyl_pos))
df['n_methyl_sites_2kb'] = df['tss'].apply(
    lambda t: count_sites_within(t, methyl_pos))

print(f'  Exposed median nearest dist:  '
      f'{df[df["is_exposed"]==1]["nearest_methyl_distance"].median():.0f} bp')
print(f'  Shielded median nearest dist: '
      f'{df[df["is_exposed"]==0]["nearest_methyl_distance"].median():.0f} bp')

# ── 4. Drop rows missing essential features ───────────────────────────────────
essential = ['baseMean', 'nearest_methyl_distance', 'n_methyl_sites_2kb', 'gene_length']
before = len(df)
df = df.dropna(subset=essential).copy()
print(f'After dropping NaN rows: {before} → {len(df)} '
      f'(exposed={df["is_exposed"].sum()}, shielded={(df["is_exposed"]==0).sum()})')

# ── 5. ROC analysis ───────────────────────────────────────────────────────────
print('\nRunning ROC analysis...')

y = df['is_exposed'].values
features = {
    'nearest_methyl_distance': ('lower_exposed', df['nearest_methyl_distance'].values),
    'n_methyl_sites_2kb':       ('higher_exposed', df['n_methyl_sites_2kb'].values),
    'baseMean':                  ('higher_exposed', df['baseMean'].values),
    'log2_baseMean_p1':          ('higher_exposed', df['log2_baseMean_p1'].values),
    'gene_length':               ('higher_exposed', df['gene_length'].values),
}

roc_rows = []
for feat_name, (direction, vals) in features.items():
    mask = ~np.isnan(vals)
    yv, xv = y[mask], vals[mask]
    if len(np.unique(yv)) < 2:
        continue
    auc = roc_auc_score(yv, xv)
    if auc < 0.5:
        auc = roc_auc_score(yv, -xv)
        direction = 'lower_exposed' if direction == 'higher_exposed' else 'higher_exposed'

    fpr, tpr, thresholds = roc_curve(yv, xv if direction == 'higher_exposed' else -xv)
    youden = tpr - fpr
    best_idx = np.argmax(youden)
    best_thresh = (thresholds[best_idx]
                   if direction == 'higher_exposed'
                   else -thresholds[best_idx])

    roc_rows.append({
        'feature': feat_name,
        'direction': direction,
        'AUC': round(auc, 4),
        'optimal_threshold': round(best_thresh, 1),
        'n_exposed': int(yv.sum()),
        'n_shielded': int((yv == 0).sum()),
    })
    print(f'  {feat_name}: AUC={auc:.4f}, threshold={best_thresh:.1f} ({direction})')

roc_df = pd.DataFrame(roc_rows).sort_values('AUC', ascending=False)

# Mann-Whitney for nearest_methyl_distance
exposed_dist  = df[df['is_exposed'] == 1]['nearest_methyl_distance']
shielded_dist = df[df['is_exposed'] == 0]['nearest_methyl_distance']
U, p_mwu = stats.mannwhitneyu(exposed_dist, shielded_dist, alternative='less')
print(f'\nMann-Whitney U (exposed < shielded): U={U:.0f}, p={p_mwu:.2e}')
print(f'Exposed   median = {exposed_dist.median():.0f} bp (n={len(exposed_dist)})')
print(f'Shielded  median = {shielded_dist.median():.0f} bp (n={len(shielded_dist)})')

# Optimal threshold for nearest_methyl_distance from ROC
nmd_row = roc_df[roc_df['feature'] == 'nearest_methyl_distance'].iloc[0]
opt_thresh = nmd_row['optimal_threshold']
auc_nmd    = nmd_row['AUC']
print(f'\nnearest_methyl_distance: AUC={auc_nmd:.4f}, Youden threshold={opt_thresh:.0f} bp')

# ── 6. Save ───────────────────────────────────────────────────────────────────
save_cols = [c for c in [
    'locus_tag', 'gene_name', 'old_locus_tag', 'product', 'tf_family',
    'start', 'end', 'strand', 'tss', 'region', 'is_exposed',
    'baseMean', 'log2_baseMean_p1', 'nearest_methyl_distance',
    'n_methyl_sites_2kb', 'gene_length', 'LFC_T2vsT1', 'LFC_T3vsT1', 'n_FIMO_hits',
] if c in df.columns]

df[save_cols].to_csv(TBL_DIR / 'all_genes_features.tsv', sep='\t', index=False)
print(f'\nSaved all_genes_features.tsv: {len(df)} rows')

roc_df.to_csv(TBL_DIR / 'ROC_analysis.tsv', sep='\t', index=False)
print(f'Saved ROC_analysis.tsv: {len(roc_df)} rows')

print('\n=== Summary ===')
print(f'Total regulatory genes (Jeong2016 TSS): {len(df)}')
print(f'  Exposed:  {df["is_exposed"].sum()}')
print(f'  Shielded: {(df["is_exposed"]==0).sum()}')
print(f'AUC (nearest_methyl_distance): {auc_nmd:.4f}')
print(f'Youden threshold: {opt_thresh:.0f} bp')
print(f'Mann-Whitney p: {p_mwu:.2e}')
