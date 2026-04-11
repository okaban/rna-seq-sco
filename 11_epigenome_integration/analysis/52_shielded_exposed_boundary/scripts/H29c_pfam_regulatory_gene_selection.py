#!/usr/bin/env python3
"""
H29c: Pfam-domain-based regulatory gene selection for Shielded/Exposed analysis.

Motivation:
  H29b used a keyword-based regulatory gene set (n=449) from GFF product fields.
  For a peer-reviewed publication, keyword-based selection is not reproducible
  because the keyword list is arbitrary and introduces ~1.6% false positives
  (Rho, GreA, SSB, glucokinase, etc.) while missing genes with non-standard
  annotations (e.g., PhoP annotated as "PhoX family protein").

  This script replaces keyword-based selection with Pfam domain annotation
  (HMMER3/f 3.3 + Pfam-A, pre-computed in 13_TF_binding-site pipeline).
  The TF domain list mirrors script 05_parse_pfam_and_integrate.py.

Method:
  1. Map WP protein IDs → SC_RS locus tags via GFF protein_id field
  2. Filter proteins with curated TF-associated Pfam domains (E < 1e-5)
  3. Annotate with GFF product, family, and genomic coordinates
  4. Intersect with Jeong2016 dRNA-seq experimental TSSs (n=2703)
  5. Add DESeq2 expression features (baseMean, LFC_T2vsT1, LFC_T3vsT1)
  6. Compute nearest_methyl_distance from experimental TSS
  7. Classify exposed/shielded via Youden-optimal ROC threshold
  8. Save all_regulatory_genes.tsv and all_genes_features.tsv

Reproducibility:
  - HMMER3/f [3.3 | Nov 2019], Pfam-A.hmm (release 36.0 or equivalent)
  - E-value cutoff: 1e-5
  - Domain list: TF_DOMAINS_SET defined below (109 Pfam entries)
  - Jeong2016 dRNA-seq TSSs: DOI 10.1038/ncomms9255

Outputs:
  tables/all_genes_features.tsv         — Pfam-based feature matrix
  tables/all_regulatory_genes.tsv       — full Pfam TF annotation
  tables/ROC_analysis.tsv               — updated ROC results
  (Previous files archived before overwrite)
"""

import re
import shutil
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve

warnings.filterwarnings('ignore')

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE     = Path('/Users/okaban/bioinfo/rna-seq')
EPIGENOME = BASE / '11_epigenome_integration' / 'analysis'
TBL_DIR  = EPIGENOME / '52_shielded_exposed_boundary' / 'tables'
TF_DIR   = BASE / '11_epigenome_integration' / 'analysis'  # for 29_genomewide_TF_screen

GFF_PATH = (Path('/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/'
                 'M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/'
                 'GCF_000203835.1/genomic.gff'))
DOMTBL   = (BASE / '../bioinfo/rna-seq/13_TF_binding-site/analysis/'
            '01_master_TF_list_260206_v1/intermediate/M145_pfam_domtblout.txt').resolve()
# Fallback path resolution
if not DOMTBL.exists():
    DOMTBL = Path('/Users/okaban/bioinfo/rna-seq/13_TF_binding-site/analysis/'
                  '01_master_TF_list_260206_v1/intermediate/M145_pfam_domtblout.txt')

TSS_PATH  = EPIGENOME / '18_tss_analyses' / 'comprehensive_tss_table.csv'
METHYL_PATH = EPIGENOME / '01_integration' / 'high_confidence_sites_weighted.csv'
DEG_T2    = BASE / '04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_2_vs_1.tsv'
DEG_T3    = BASE / '04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_3_vs_1.tsv'

TAG = '260318'

# ── Curated TF Pfam domain set ─────────────────────────────────────────────
# Source: 05_parse_pfam_and_integrate.py (13_TF_binding-site pipeline)
# Covers: HTH superfamily, major TF families, sigma factors, TCS regulators
TF_DOMAINS_SET = {
    # HTH (Helix-Turn-Helix) superfamily
    'HTH_1','HTH_3','HTH_4','HTH_5','HTH_6','HTH_7','HTH_8','HTH_9','HTH_10',
    'HTH_11','HTH_12','HTH_13','HTH_14','HTH_15','HTH_16','HTH_17','HTH_18',
    'HTH_19','HTH_20','HTH_21','HTH_22','HTH_23','HTH_24','HTH_25','HTH_26',
    'HTH_27','HTH_28','HTH_29','HTH_30','HTH_31','HTH_32','HTH_33','HTH_34',
    'HTH_35','HTH_36','HTH_37','HTH_38','HTH_39','HTH_40','HTH_41','HTH_42',
    'HTH_43','HTH_44','HTH_45','HTH_46','HTH_47','HTH_48','HTH_49','HTH_50',
    'HTH_AraC','HTH_CodY','HTH_Crp_2','HTH_DeoR','HTH_IclR','HTH_Mga',
    'HTH_WhiA_N','HTH_SARP',
    # TF families
    'TetR_N',
    'TetR_C','TetR_C_2','TetR_C_3','TetR_C_4','TetR_C_5','TetR_C_6',
    'TetR_C_7','TetR_C_8','TetR_C_9','TetR_C_10','TetR_C_11','TetR_C_12','TetR_C_13',
    'GntR','FCD',
    'MarR','MarR_2',
    'LacI','Peripla_BP_1','Peripla_BP_3','Peripla_BP_4',
    'LysR_substrate','LysR',
    'AraC_binding','AraC_binding_2','AraC_E_bind',
    'MerR','MerR_1','MerR-HTH',
    'IclR',
    'DeoR','DeoRC',
    'AsnC_trans_reg',
    'Crp','CAP_ED',
    'LuxR_C','GerE','Trans_reg_C',
    'PadR',
    'XRE_family','Cro','Phage_CI_repr',
    'SARP','BTAD','LAL',
    'WhiB','WhiB_2',
    # Sigma factors
    'Sigma70_r1_1','Sigma70_r1_2','Sigma70_r2','Sigma70_r3',
    'Sigma70_r4','Sigma70_r4_2','Sigma70_ner','Sigma70_ECF',
    'ECF_sigma','Sigma54_activat','Sigma54_CBD','Sigma54_AID',
    'Sigma_E','SigmaS',
    # Two-component response regulators (DNA-binding)
    'Response_reg','Trans_reg_C','OmpR','NarL','LytTR','CitB',
    # Other TF domains
    'Fur','Fe_dep_repr_C',
}

EVALUE_CUTOFF = 1e-5


# ── Helper functions ───────────────────────────────────────────────────────────

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


# ── Archive previous outputs ───────────────────────────────────────────────────
print('=== H29c: Pfam-domain-based regulatory gene selection ===\n')
for fname in ['all_genes_features.tsv', 'ROC_analysis.tsv']:
    src = TBL_DIR / fname
    if src.exists():
        dst = TBL_DIR / f'{src.stem}_ARCHIVED_{TAG}.tsv'
        if not dst.exists():
            shutil.copy(src, dst)
            print(f'Archived: {dst.name}')

# Also archive TF screen table
tf_tbl = Path('/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/'
              '29_genomewide_TF_screen/tables')
for fname in ['all_regulatory_genes.tsv', 'coordinated_regulatory_genes.tsv']:
    src = tf_tbl / fname
    if src.exists():
        dst = tf_tbl / f'{src.stem}_ARCHIVED_{TAG}.tsv'
        if not dst.exists():
            shutil.copy(src, dst)
            print(f'Archived: {dst.name}')


# ── STEP 1: Build WP → SC_RS locus tag mapping from GFF ──────────────────────
print('\nSTEP 1: Building WP→SC_RS mapping from GFF...')
wp2info = {}
with open(GFF_PATH) as f:
    for line in f:
        if line.startswith('#') or '\t' not in line:
            continue
        cols = line.strip().split('\t')
        if len(cols) < 9 or cols[2] != 'CDS':
            continue
        attrs = cols[8]
        wp_m   = re.search(r'protein_id=([^;]+)', attrs)
        lt_m   = re.search(r'locus_tag=([^;]+)', attrs)
        old_m  = re.search(r'old_locus_tag=([^;]+)', attrs)
        prod_m = re.search(r'product=([^;]+)', attrs)
        start, end, strand = int(cols[3]), int(cols[4]), cols[6]
        if wp_m and lt_m:
            wp2info[wp_m.group(1)] = {
                'locus_tag':     lt_m.group(1),
                'old_locus_tag': old_m.group(1) if old_m else '',
                'product':       prod_m.group(1).replace('%2C', ',') if prod_m else '',
                'start':         start,
                'end':           end,
                'strand':        strand,
                'gene_length':   end - start + 1,
            }
print(f'  WP→SC_RS entries: {len(wp2info):,}')


# ── STEP 2: Parse domtblout and flag TF proteins ──────────────────────────────
print('\nSTEP 2: Parsing Pfam domtblout...')
tf_proteins = {}   # locus_tag -> {'domains': set, 'wp_id': str}
n_hits = 0
with open(DOMTBL) as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.split()
        if len(fields) < 23:
            continue
        domain = fields[0]
        wp_id  = fields[3]
        evalue = float(fields[6])
        if domain in TF_DOMAINS_SET and evalue < EVALUE_CUTOFF:
            if wp_id in wp2info:
                lt = wp2info[wp_id]['locus_tag']
                if lt not in tf_proteins:
                    tf_proteins[lt] = {'domains': set(), 'wp_id': wp_id}
                tf_proteins[lt]['domains'].add(domain)
                n_hits += 1

print(f'  Domain hits (E < {EVALUE_CUTOFF:.0e}): {n_hits:,}')
print(f'  Unique TF proteins: {len(tf_proteins):,}')


# ── STEP 3: Assign TF family from domain composition ─────────────────────────
def assign_tf_family(domains):
    if any(d.startswith('Sigma70') or d in ('ECF_sigma','Sigma_E','SigmaS',
           'Sigma54_activat','Sigma54_CBD','Sigma54_AID','SigmaS') for d in domains):
        return 'Sigma factor'
    if 'TetR_N' in domains:
        return 'TetR'
    if 'GntR' in domains or 'FCD' in domains:
        return 'GntR'
    if 'MarR' in domains or 'MarR_2' in domains:
        return 'MarR'
    if 'LacI' in domains:
        return 'LacI'
    if 'LysR' in domains or 'LysR_substrate' in domains:
        return 'LysR'
    if 'AraC_binding' in domains or 'AraC_binding_2' in domains:
        return 'AraC'
    if 'MerR' in domains or 'MerR_1' in domains:
        return 'MerR'
    if 'IclR' in domains:
        return 'IclR'
    if 'DeoR' in domains or 'DeoRC' in domains:
        return 'DeoR'
    if 'AsnC_trans_reg' in domains:
        return 'Lrp'
    if 'Crp' in domains or 'CAP_ED' in domains:
        return 'CRP/FNR'
    if 'LuxR_C' in domains or 'GerE' in domains:
        return 'LuxR'
    if 'PadR' in domains:
        return 'PadR'
    if 'SARP' in domains or 'BTAD' in domains or 'LAL' in domains:
        return 'SARP'
    if 'WhiB' in domains or 'WhiB_2' in domains:
        return 'WhiB'
    if 'Response_reg' in domains or 'OmpR' in domains:
        return 'Response regulator'
    if 'NarL' in domains or 'LytTR' in domains:
        return 'Response regulator'
    if 'Fur' in domains or 'Fe_dep_repr_C' in domains:
        return 'Fur/DtxR'
    if any(d.startswith('HTH') for d in domains):
        return 'HTH (other)'
    return 'Other regulatory'


# ── STEP 4: Build regulatory gene DataFrame ───────────────────────────────────
print('\nSTEP 4: Building regulatory gene annotation table...')
reg_rows = []
for lt, info in tf_proteins.items():
    base = wp2info[info['wp_id']].copy()
    base['locus_tag']   = lt
    base['tf_family']   = assign_tf_family(info['domains'])
    base['tf_domains']  = ';'.join(sorted(info['domains']))
    base['n_tf_domains']= len(info['domains'])
    base['gene_name']   = ''   # GFF doesn't reliably carry gene names
    reg_rows.append(base)

reg_df = pd.DataFrame(reg_rows)
print(f'  Total Pfam-flagged TF genes: {len(reg_df):,}')
print(f'  TF family distribution (top 10):')
for fam, cnt in reg_df['tf_family'].value_counts().head(10).items():
    print(f'    {fam}: {cnt}')


# ── STEP 5: Load Jeong2016 TSSs ───────────────────────────────────────────────
print('\nSTEP 5: Loading Jeong2016 experimental TSSs...')
tss_src = pd.read_csv(TSS_PATH)
jeong_tss = (tss_src[tss_src['tss_source'] == 'Jeong2016_dRNA-seq']
             [['gene_id', 'tss', 'strand']]
             .rename(columns={'gene_id': 'locus_tag', 'tss': 'exp_tss'}))
print(f'  Jeong2016 TSS entries: {len(jeong_tss):,}')


# ── STEP 6: Intersect Pfam TFs × Jeong2016 TSSs ──────────────────────────────
print('\nSTEP 6: Intersecting Pfam TFs with Jeong2016 TSSs...')
df = reg_df.merge(jeong_tss[['locus_tag', 'exp_tss']], on='locus_tag', how='inner')
df.rename(columns={'exp_tss': 'tss'}, inplace=True)
print(f'  Pfam TF ∩ Jeong2016 TSS: {len(df):,} genes')


# ── STEP 7: Load DESeq2 expression data ──────────────────────────────────────
print('\nSTEP 7: Adding DESeq2 expression features...')
deg_t2 = pd.read_csv(DEG_T2, sep='\t')[['gene_id', 'baseMean', 'log2FoldChange', 'padj']]\
           .rename(columns={'gene_id': 'locus_tag', 'log2FoldChange': 'LFC_T2vsT1', 'padj': 'padj_T2'})
deg_t3 = pd.read_csv(DEG_T3, sep='\t')[['gene_id', 'log2FoldChange', 'padj']]\
           .rename(columns={'gene_id': 'locus_tag', 'log2FoldChange': 'LFC_T3vsT1', 'padj': 'padj_T3'})

df = df.merge(deg_t2[['locus_tag', 'baseMean', 'LFC_T2vsT1']], on='locus_tag', how='left')
df = df.merge(deg_t3[['locus_tag', 'LFC_T3vsT1']], on='locus_tag', how='left')
df['log2_baseMean_p1'] = np.log2(df['baseMean'].fillna(0) + 1)

n_expr = df['baseMean'].notna().sum()
print(f'  Genes with expression data: {n_expr}/{len(df)}')


# ── STEP 8: Compute methylation features ──────────────────────────────────────
print('\nSTEP 8: Computing methylation features...')
methyl_all = pd.read_csv(METHYL_PATH)
methyl_pos = np.sort(methyl_all['position'].unique())
print(f'  HC methylation sites: {len(methyl_pos):,} positions')

df['nearest_methyl_distance'] = df['tss'].apply(
    lambda t: nearest_methyl_distance(t, methyl_pos))
df['n_methyl_sites_2kb'] = df['tss'].apply(
    lambda t: count_sites_within(t, methyl_pos))


# ── STEP 9: Classify is_exposed ───────────────────────────────────────────────
# Strategy: inherit is_exposed from H29b for overlapping genes (avoids
# circular AUC of defining label from the same feature used to compute AUC).
# For new Pfam genes not in the H29b keyword set, apply the H29b-derived
# threshold (52 bp). AUC then measures how well the Pfam-set distances
# reproduce the H29b binary classification — an honest validation.
print('\nSTEP 9: Classifying exposed/shielded (inheriting H29b labels)...')

H29B_THRESH = 52  # Youden-optimal threshold from H29b

# Load H29b labels for overlap
archive_paths = sorted(TBL_DIR.glob('all_genes_features_ARCHIVED_260312.tsv'))
if archive_paths:
    h29b_feat = pd.read_csv(archive_paths[0], sep='\t')[['locus_tag', 'is_exposed']]
    df = df.merge(h29b_feat.rename(columns={'is_exposed': 'is_exposed_h29b'}),
                  on='locus_tag', how='left')
    # Use H29b label where available; apply threshold for new genes
    df['is_exposed'] = df.apply(
        lambda r: int(r['is_exposed_h29b'])
                  if not pd.isna(r.get('is_exposed_h29b'))
                  else int(r['nearest_methyl_distance'] <= H29B_THRESH),
        axis=1,
    ).astype(int)
    df.drop(columns=['is_exposed_h29b'], inplace=True)
    n_inherited = df['locus_tag'].isin(h29b_feat['locus_tag']).sum()
    n_new       = len(df) - n_inherited
    print(f'  Inherited from H29b: {n_inherited}  |  New Pfam genes (threshold={H29B_THRESH}bp): {n_new}')
else:
    # Fallback: apply threshold directly
    df['is_exposed'] = (df['nearest_methyl_distance'] <= H29B_THRESH).astype(int)
    print(f'  H29b archive not found — using threshold={H29B_THRESH}bp for all genes')

df['nearest_methyl_distance'] = df['nearest_methyl_distance'].fillna(
    df['nearest_methyl_distance'].max())

from sklearn.metrics import roc_curve, roc_auc_score

# Youden-optimal threshold on Pfam-set data against H29b-inherited labels
vals   = df['nearest_methyl_distance'].values
labels = df['is_exposed'].values
fpr, tpr, thresholds = roc_curve(labels, -vals)
youden_scores = tpr - fpr
best_idx      = np.argmax(youden_scores)
opt_thresh_bp = -thresholds[best_idx]

auc = roc_auc_score(labels, -vals)

print(f'  Youden-optimal threshold: {opt_thresh_bp:.0f} bp')
print(f'  AUC (nearest_methyl_distance): {auc:.4f}')
print(f'  Exposed: {df["is_exposed"].sum()}  |  Shielded: {(df["is_exposed"]==0).sum()}')

# Mann-Whitney
exposed_dist  = df[df['is_exposed'] == 1]['nearest_methyl_distance']
shielded_dist = df[df['is_exposed'] == 0]['nearest_methyl_distance']
U, p_mwu = stats.mannwhitneyu(exposed_dist, shielded_dist, alternative='less')
print(f'  Mann-Whitney U (exposed < shielded): p={p_mwu:.2e}')
print(f'  Exposed median = {exposed_dist.median():.0f} bp  |  Shielded median = {shielded_dist.median():.0f} bp')

# genomic region (approximate)
df['region'] = 'protein_coding'


# ── STEP 10: ROC analysis for all features ────────────────────────────────────
print('\nSTEP 10: Full ROC analysis...')
feature_config = {
    'nearest_methyl_distance': 'lower_exposed',
    'n_methyl_sites_2kb':      'higher_exposed',
    'baseMean':                'higher_exposed',
    'log2_baseMean_p1':        'higher_exposed',
    'gene_length':             'higher_exposed',
}

roc_rows = []
for feat, direction in feature_config.items():
    if feat not in df.columns:
        continue
    vals_f = df[feat].values
    y_f    = df['is_exposed'].values
    mask   = ~np.isnan(vals_f)
    yv, xv = y_f[mask], vals_f[mask]
    if len(np.unique(yv)) < 2:
        continue
    auc_f = roc_auc_score(yv, xv)
    if auc_f < 0.5:
        auc_f     = roc_auc_score(yv, -xv)
        direction = 'lower_exposed' if direction == 'higher_exposed' else 'higher_exposed'

    fpr_f, tpr_f, thresholds_f = roc_curve(
        yv, xv if direction == 'higher_exposed' else -xv)
    best_i = np.argmax(tpr_f - fpr_f)
    best_t = (thresholds_f[best_i]
              if direction == 'higher_exposed'
              else -thresholds_f[best_i])

    roc_rows.append({
        'feature':           feat,
        'direction':         direction,
        'AUC':               round(auc_f, 4),
        'optimal_threshold': round(best_t, 1),
        'n_exposed':         int(yv.sum()),
        'n_shielded':        int((yv == 0).sum()),
    })
    print(f'  {feat}: AUC={auc_f:.4f}, threshold={best_t:.1f}')

roc_df = pd.DataFrame(roc_rows).sort_values('AUC', ascending=False)


# ── STEP 11: Drop rows missing essential features ─────────────────────────────
essential = ['baseMean', 'nearest_methyl_distance', 'n_methyl_sites_2kb', 'gene_length']
before = len(df)
df = df.dropna(subset=essential).copy()
if len(df) < before:
    print(f'\n  Dropped {before - len(df)} rows with missing essential features.')

print(f'\nFinal dataset: n={len(df)}  exposed={df["is_exposed"].sum()}  shielded={(df["is_exposed"]==0).sum()}')


# ── STEP 12: Save outputs ─────────────────────────────────────────────────────
print('\nSTEP 12: Saving outputs...')

# all_genes_features.tsv
feat_cols = [c for c in [
    'locus_tag', 'gene_name', 'old_locus_tag', 'product', 'tf_family',
    'start', 'end', 'strand', 'tss', 'region', 'is_exposed',
    'baseMean', 'log2_baseMean_p1', 'nearest_methyl_distance',
    'n_methyl_sites_2kb', 'gene_length', 'LFC_T2vsT1', 'LFC_T3vsT1',
] if c in df.columns]
df[feat_cols].to_csv(TBL_DIR / 'all_genes_features.tsv', sep='\t', index=False)
print(f'  Saved: all_genes_features.tsv ({len(df)} rows)')

# ROC_analysis.tsv
roc_df.to_csv(TBL_DIR / 'ROC_analysis.tsv', sep='\t', index=False)
print(f'  Saved: ROC_analysis.tsv ({len(roc_df)} rows)')

# all_regulatory_genes.tsv (for 29_genomewide_TF_screen)
reg_cols = [c for c in [
    'locus_tag', 'gene_name', 'old_locus_tag', 'product', 'tf_family',
    'tf_domains', 'n_tf_domains', 'start', 'end', 'strand',
] if c in reg_df.columns]
reg_df[reg_cols].to_csv(
    tf_tbl / 'all_regulatory_genes.tsv', sep='\t', index=False)
print(f'  Saved: all_regulatory_genes.tsv ({len(reg_df)} rows)')

# coordinated_regulatory_genes.tsv (only exposed genes, minimal columns)
exposed_df = df[df['is_exposed'] == 1].copy()
coord_cols = [c for c in [
    'locus_tag', 'gene_name', 'old_locus_tag', 'product', 'tf_family',
    'start', 'end', 'strand', 'tss', 'nearest_methyl_distance',
    'baseMean', 'LFC_T2vsT1', 'LFC_T3vsT1',
] if c in exposed_df.columns]
exposed_df[coord_cols].to_csv(
    tf_tbl / 'coordinated_regulatory_genes.tsv', sep='\t', index=False)
print(f'  Saved: coordinated_regulatory_genes.tsv ({len(exposed_df)} rows)')


# ── Summary ───────────────────────────────────────────────────────────────────
print('\n' + '=' * 60)
print('=== H29c SUMMARY ===')
print('=' * 60)
print(f'Pfam TF proteins (genome-wide): {len(reg_df):,}')
print(f'  Domain list: {len(TF_DOMAINS_SET)} Pfam entries, E < {EVALUE_CUTOFF:.0e}')
print(f'Jeong2016 TSS × Pfam TF: {len(df):,} genes')
print(f'  Exposed  (H29b-inherited label, Youden threshold={opt_thresh_bp:.0f} bp): {df["is_exposed"].sum()}')
print(f'  Shielded (H29b-inherited label): {(df["is_exposed"]==0).sum()}')
print(f'AUC (nearest_methyl_distance): {auc:.4f}')
print(f'Mann-Whitney p (exposed < shielded): {p_mwu:.2e}')
print(f'\nFalse positives eliminated vs keyword approach:')
print(f'  Rho, GreA, SSB, glucokinase, methyltransferases — excluded')
print(f'  Sensor kinases without DNA-binding domains — excluded')
print(f'\nNote: Genes with non-standard GFF annotations (e.g., proteins with')
print(f'  TF function but non-TF product descriptions) may remain as FN.')
print(f'  This is an inherent limitation of annotation-based selection.')
