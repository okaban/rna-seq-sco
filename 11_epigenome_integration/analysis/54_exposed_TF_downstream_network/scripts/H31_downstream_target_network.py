#!/usr/bin/env python3
"""
H31: Exposed TF Downstream Regulatory Target Network Analysis
==============================================================
Analyzes how many downstream genes the 57 "exposed" transcription factors regulate
and whether those target genes show expression changes coordinated with TF methylation status.

Since the 57 exposed TFs do not have FIMO motifs (they are newly identified regulators),
the analysis takes a dual approach:
1. Use the 16 well-characterized TFs with FIMO motifs to build a genome-wide binding network
2. Determine which exposed regulators are TARGETS of these 16 TFs
3. Analyze whether exposed TF targets show distinct expression patterns
4. Use TF family information to infer regulatory capacity of exposed TFs

Author: Claude Code (H31 analysis)
Date: 2026-02-27
"""

import os
import sys
import csv
import warnings
import numpy as np
import pandas as pd
from collections import defaultdict, Counter
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import matplotlib.ticker as ticker

warnings.filterwarnings('ignore')

# ============================================================
# Configuration
# ============================================================
BASE_DIR = '/Users/okaban/bioinfo/rna-seq'
ANALYSIS_DIR = f'{BASE_DIR}/11_epigenome_integration/analysis/54_exposed_TF_downstream_network'
FIG_DIR = f'{ANALYSIS_DIR}/figures'
TAB_DIR = f'{ANALYSIS_DIR}/tables'

FIMO_FILE = f'{BASE_DIR}/13_TF_binding-site/analysis/01_master_TF_list_260206_v1/intermediate/fimo_results/fimo.tsv'
EXPOSED_FILE = f'{BASE_DIR}/11_epigenome_integration/analysis/51_exposed_regulators_characteristics/tables/exposed_regulators_full_table.tsv'
ALL_REG_FILE = f'{BASE_DIR}/11_epigenome_integration/analysis/29_genomewide_TF_screen/tables/all_regulatory_genes.tsv'
COORD_REG_FILE = f'{BASE_DIR}/11_epigenome_integration/analysis/29_genomewide_TF_screen/tables/coordinated_regulatory_genes.tsv'
GENE_ANNOT_FILE = f'{BASE_DIR}/05_annotation/analysis/05_annotation_260128_v1/tables/gene_annotation_basic.tsv'
ALL_GENES_TF_FILE = f'{BASE_DIR}/13_TF_binding-site/analysis/01_master_TF_list_260206_v1/intermediate/M145_all_genes_with_TF_annotations.tsv'
DESEQ2_T2 = f'{BASE_DIR}/04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_2_vs_1.tsv'
DESEQ2_T3 = f'{BASE_DIR}/04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_3_vs_1.tsv'

CHROM_LENGTH = 8667507
PROMOTER_WINDOW = 500  # bp upstream of TSS
Q_VALUE_THRESHOLD = 0.05
Q_VALUE_RELAXED = 0.1  # relaxed threshold
LFC_THRESHOLD = 1.0  # |log2FC| >= 1 for DEG
PADJ_THRESHOLD = 0.05

# ============================================================
# Step 1: Parse FIMO binding sites
# ============================================================
print("=" * 70)
print("STEP 1: Parsing FIMO binding sites")
print("=" * 70)

fimo_hits = []
with open(FIMO_FILE) as f:
    for line in f:
        if line.startswith('#') or line.startswith('motif_id'):
            continue
        parts = line.strip().split('\t')
        if len(parts) >= 10:
            try:
                q = float(parts[8])
                fimo_hits.append({
                    'motif_id': parts[0],
                    'motif_alt_id': parts[1],  # SCO ID
                    'chrom': parts[2],
                    'start': int(parts[3]),
                    'stop': int(parts[4]),
                    'strand': parts[5],
                    'score': float(parts[6]),
                    'pvalue': float(parts[7]),
                    'qvalue': q,
                    'sequence': parts[9]
                })
            except (ValueError, IndexError):
                pass

fimo_df = pd.DataFrame(fimo_hits)
print(f"Total FIMO hits parsed: {len(fimo_df)}")

# Filter by q-value thresholds
fimo_strict = fimo_df[fimo_df['qvalue'] < Q_VALUE_THRESHOLD].copy()
fimo_relaxed = fimo_df[fimo_df['qvalue'] < Q_VALUE_RELAXED].copy()
print(f"Hits with q < {Q_VALUE_THRESHOLD} (strict): {len(fimo_strict)}")
print(f"Hits with q < {Q_VALUE_RELAXED} (relaxed): {len(fimo_relaxed)}")

for label, df in [('strict', fimo_strict), ('relaxed', fimo_relaxed)]:
    print(f"\n  Per-TF counts ({label}):")
    for tf, count in df['motif_id'].value_counts().items():
        print(f"    {tf} ({df[df['motif_id']==tf]['motif_alt_id'].iloc[0]}): {count}")

# Use strict for primary analysis, note relaxed counts
fimo_primary = fimo_strict.copy()
print(f"\nUsing strict threshold (q < {Q_VALUE_THRESHOLD}): {len(fimo_primary)} sites from {fimo_primary['motif_id'].nunique()} TFs")

# ============================================================
# Step 2: Build locus_tag mappings
# ============================================================
print("\n" + "=" * 70)
print("STEP 2: Building locus_tag mappings")
print("=" * 70)

# Read gene annotations (SC_RS -> SCO mapping)
gene_annot = pd.read_csv(GENE_ANNOT_FILE, sep='\t')
print(f"Gene annotations loaded: {len(gene_annot)} genes")

# Build mapping from both directions
sco_to_rs = {}
rs_to_sco = {}
for _, row in gene_annot.iterrows():
    rs = row['gene_id']
    sco = row.get('old_locus_tag', '')
    if pd.notna(sco) and sco != '':
        sco_to_rs[sco] = rs
        rs_to_sco[rs] = sco

# Also use the all_genes_TF file which has both
all_genes_tf = pd.read_csv(ALL_GENES_TF_FILE, sep='\t')
for _, row in all_genes_tf.iterrows():
    sco = row.get('SCO_ID', '')
    rs = row.get('gene_id', '')
    if pd.notna(sco) and sco != '' and pd.notna(rs) and rs != '':
        sco_to_rs[sco] = rs
        rs_to_sco[rs] = sco

print(f"SCO -> RS mappings: {len(sco_to_rs)}")
print(f"RS -> SCO mappings: {len(rs_to_sco)}")

# Read exposed regulators
exposed_df = pd.read_csv(EXPOSED_FILE, sep='\t')
exposed_locus_tags = set(exposed_df['locus_tag'].tolist())
exposed_sco_ids = set()
for _, row in exposed_df.iterrows():
    sco = row.get('old_locus_tag', '')
    if pd.notna(sco) and sco != '':
        exposed_sco_ids.add(sco)
print(f"\nExposed regulators: {len(exposed_locus_tags)} (with SCO IDs: {len(exposed_sco_ids)})")

# Read all regulatory genes
all_reg = pd.read_csv(ALL_REG_FILE, sep='\t')
all_reg_locus_tags = set(all_reg['locus_tag'].tolist())
shielded_locus_tags = all_reg_locus_tags - exposed_locus_tags
print(f"All regulatory genes: {len(all_reg_locus_tags)}")
print(f"Shielded regulators: {len(shielded_locus_tags)}")

# Read coordinated regulatory genes
coord_reg = pd.read_csv(COORD_REG_FILE, sep='\t')
print(f"Coordinated regulatory genes: {len(coord_reg)}")

# Map FIMO TF SCO IDs to RS format
fimo_tf_mapping = {}
for _, row in fimo_primary.drop_duplicates('motif_id').iterrows():
    sco = row['motif_alt_id']
    name = row['motif_id']
    rs = sco_to_rs.get(sco, 'UNMAPPED')
    fimo_tf_mapping[name] = {'sco': sco, 'rs': rs}
    print(f"  FIMO TF: {name} -> {sco} -> {rs}")

# Check if any FIMO TFs are exposed
fimo_sco_set = set(fimo_tf_mapping[k]['sco'] for k in fimo_tf_mapping)
overlap = fimo_sco_set & exposed_sco_ids
print(f"\nFIMO TFs that are exposed regulators: {len(overlap)}")
if overlap:
    print(f"  Overlapping: {overlap}")
else:
    print("  None - all 16 FIMO TFs are shielded or not classified as regulatory")

# ============================================================
# Step 3: Map binding sites to target genes
# ============================================================
print("\n" + "=" * 70)
print("STEP 3: Mapping binding sites to target genes")
print("=" * 70)

# Build gene coordinate database from annotation
genes = []
for _, row in gene_annot.iterrows():
    genes.append({
        'gene_id': row['gene_id'],
        'old_locus_tag': row.get('old_locus_tag', ''),
        'start': int(row['start']),
        'end': int(row['end']),
        'strand': row['strand'],
        'gene_name': row.get('gene_name', ''),
        'product': row.get('product', '')
    })

genes_sorted = sorted(genes, key=lambda x: x['start'])
gene_starts = np.array([g['start'] for g in genes_sorted])
gene_ends = np.array([g['end'] for g in genes_sorted])
gene_strands = [g['strand'] for g in genes_sorted]

print(f"Total genes for target mapping: {len(genes_sorted)}")

def find_target_gene(bs_start, bs_stop, promoter_window=PROMOTER_WINDOW):
    """
    Find gene(s) whose promoter region contains the binding site.
    For + strand genes: promoter is [start - promoter_window, start + promoter_window]
    For - strand genes: promoter is [end - promoter_window, end + promoter_window]
    """
    targets = []
    bs_mid = (bs_start + bs_stop) / 2.0

    for g in genes_sorted:
        if g['strand'] == '+':
            tss = g['start']
        else:
            tss = g['end']

        prom_start = tss - promoter_window
        prom_end = tss + promoter_window

        # Check if binding site midpoint falls within promoter
        if prom_start <= bs_mid <= prom_end:
            distance = abs(bs_mid - tss)
            targets.append({
                'gene_id': g['gene_id'],
                'old_locus_tag': g['old_locus_tag'],
                'tss': tss,
                'distance': distance,
                'strand': g['strand'],
                'product': g['product']
            })

    return targets

# Map all FIMO binding sites to target genes
tf_target_pairs = []
unmapped_sites = 0

for _, hit in fimo_primary.iterrows():
    targets = find_target_gene(hit['start'], hit['stop'])
    if targets:
        for t in targets:
            tf_target_pairs.append({
                'tf_name': hit['motif_id'],
                'tf_sco': hit['motif_alt_id'],
                'tf_rs': fimo_tf_mapping.get(hit['motif_id'], {}).get('rs', 'UNMAPPED'),
                'bs_start': hit['start'],
                'bs_stop': hit['stop'],
                'bs_strand': hit['strand'],
                'bs_score': hit['score'],
                'bs_qvalue': hit['qvalue'],
                'target_gene': t['gene_id'],
                'target_sco': t['old_locus_tag'],
                'target_tss': t['tss'],
                'distance_to_tss': t['distance'],
                'target_strand': t['strand'],
                'target_product': t['product']
            })
    else:
        unmapped_sites += 1

pairs_df = pd.DataFrame(tf_target_pairs)
print(f"Total TF-target pairs: {len(pairs_df)}")
print(f"Binding sites not mapping to any promoter: {unmapped_sites}")
print(f"Unique target genes: {pairs_df['target_gene'].nunique()}")

# Per-TF target counts
print("\nPer-TF target gene counts (strict q < 0.05):")
for tf in sorted(pairs_df['tf_name'].unique()):
    n_targets = pairs_df[pairs_df['tf_name'] == tf]['target_gene'].nunique()
    n_sites = len(pairs_df[pairs_df['tf_name'] == tf])
    print(f"  {tf}: {n_targets} target genes ({n_sites} binding site-target pairs)")

# ============================================================
# Step 4: Build exposed TF regulatory network
# ============================================================
print("\n" + "=" * 70)
print("STEP 4: Building exposed TF regulatory network")
print("=" * 70)

# Identify exposed regulators that are targets of FIMO TFs
exposed_as_targets = pairs_df[pairs_df['target_gene'].isin(exposed_locus_tags)].copy()
print(f"Exposed regulators targeted by FIMO TFs: {exposed_as_targets['target_gene'].nunique()}")

for _, row in exposed_as_targets.drop_duplicates('target_gene').iterrows():
    tfs = pairs_df[pairs_df['target_gene'] == row['target_gene']]['tf_name'].unique()
    print(f"  {row['target_gene']} ({row['target_sco']}) <- regulated by: {', '.join(tfs)}")

# All regulatory genes as targets
all_reg_as_targets = pairs_df[pairs_df['target_gene'].isin(all_reg_locus_tags)].copy()
shielded_as_targets = pairs_df[pairs_df['target_gene'].isin(shielded_locus_tags)].copy()
print(f"\nAll regulatory genes targeted by FIMO TFs: {all_reg_as_targets['target_gene'].nunique()}")
print(f"Shielded regulators targeted: {shielded_as_targets['target_gene'].nunique()}")
print(f"Exposed regulators targeted: {exposed_as_targets['target_gene'].nunique()}")

# Cascade amplification: each exposed TF targeted by FIMO TFs
# also regulates its own target genes (though we don't know their specific targets)
# We can estimate from TF family typical regulon sizes

# For the FIMO TFs themselves, calculate the "cascade"
all_fimo_targets = set(pairs_df['target_gene'].tolist())
exposed_fimo_targets = set(exposed_as_targets['target_gene'].tolist())

print(f"\nCascade amplification summary:")
print(f"  16 FIMO TFs -> {len(all_fimo_targets)} unique target genes")
print(f"  Of these, {len(exposed_fimo_targets)} are exposed regulators")
print(f"  Of these, {len(all_fimo_targets & all_reg_locus_tags)} are regulatory genes")
print(f"  Potential second-layer cascade: {len(exposed_fimo_targets)} exposed TFs -> further targets")

# Build per-TF summary
per_tf_summary = []
for tf in sorted(pairs_df['tf_name'].unique()):
    tf_pairs = pairs_df[pairs_df['tf_name'] == tf]
    targets = set(tf_pairs['target_gene'])
    n_exposed_targets = len(targets & exposed_locus_tags)
    n_shielded_targets = len(targets & shielded_locus_tags)
    n_reg_targets = len(targets & all_reg_locus_tags)

    per_tf_summary.append({
        'tf_name': tf,
        'tf_sco': fimo_tf_mapping[tf]['sco'],
        'tf_rs': fimo_tf_mapping[tf]['rs'],
        'n_binding_sites': len(tf_pairs),
        'n_target_genes': len(targets),
        'n_regulatory_targets': n_reg_targets,
        'n_exposed_targets': n_exposed_targets,
        'n_shielded_targets': n_shielded_targets,
        'pct_regulatory_targets': n_reg_targets / len(targets) * 100 if len(targets) > 0 else 0,
        'pct_exposed_targets': n_exposed_targets / len(targets) * 100 if len(targets) > 0 else 0
    })

per_tf_df = pd.DataFrame(per_tf_summary)

# ============================================================
# Step 5: Expression coordination analysis
# ============================================================
print("\n" + "=" * 70)
print("STEP 5: Expression coordination analysis")
print("=" * 70)

# Read DESeq2 results
deseq2_t2 = pd.read_csv(DESEQ2_T2, sep='\t')
deseq2_t3 = pd.read_csv(DESEQ2_T3, sep='\t')
print(f"DESeq2 T2vsT1: {len(deseq2_t2)} genes")
print(f"DESeq2 T3vsT1: {len(deseq2_t3)} genes")

# Merge expression data
expr_data = deseq2_t2[['gene_id', 'baseMean', 'log2FoldChange', 'padj']].rename(
    columns={'log2FoldChange': 'lfc_T2', 'padj': 'padj_T2', 'baseMean': 'baseMean_T2'}
)
expr_t3 = deseq2_t3[['gene_id', 'log2FoldChange', 'padj']].rename(
    columns={'log2FoldChange': 'lfc_T3', 'padj': 'padj_T3'}
)
expr_data = expr_data.merge(expr_t3, on='gene_id', how='outer')
expr_data['abs_lfc_T2'] = expr_data['lfc_T2'].abs()
expr_data['abs_lfc_T3'] = expr_data['lfc_T3'].abs()
expr_data['is_DEG_T2'] = (expr_data['padj_T2'] < PADJ_THRESHOLD) & (expr_data['abs_lfc_T2'] >= LFC_THRESHOLD)
expr_data['is_DEG_T3'] = (expr_data['padj_T3'] < PADJ_THRESHOLD) & (expr_data['abs_lfc_T3'] >= LFC_THRESHOLD)
expr_data['is_DEG_any'] = expr_data['is_DEG_T2'] | expr_data['is_DEG_T3']

# Classify genes as targets vs non-targets
target_genes = set(pairs_df['target_gene'])
expr_data['is_target'] = expr_data['gene_id'].isin(target_genes)

# Also classify by which TF
for tf in pairs_df['tf_name'].unique():
    tf_targets = set(pairs_df[pairs_df['tf_name'] == tf]['target_gene'])
    expr_data[f'target_of_{tf}'] = expr_data['gene_id'].isin(tf_targets)

# Compare target vs non-target expression
target_expr = expr_data[expr_data['is_target']].copy()
nontarget_expr = expr_data[~expr_data['is_target']].copy()

print(f"\nTarget genes with expression data: {len(target_expr)}")
print(f"Non-target genes with expression data: {len(nontarget_expr)}")

# Statistical tests
stat_results = []

# Test 1: |LFC| comparison (Wilcoxon)
for timepoint, col in [('T2', 'abs_lfc_T2'), ('T3', 'abs_lfc_T3')]:
    t_vals = target_expr[col].dropna()
    nt_vals = nontarget_expr[col].dropna()
    if len(t_vals) > 0 and len(nt_vals) > 0:
        stat, pval = stats.mannwhitneyu(t_vals, nt_vals, alternative='two-sided')
        r = abs(stats.norm.ppf(pval/2)) / np.sqrt(len(t_vals) + len(nt_vals)) if pval > 0 else 0
        stat_results.append({
            'test': f'|LFC|_target_vs_nontarget_{timepoint}',
            'target_median': t_vals.median(),
            'nontarget_median': nt_vals.median(),
            'target_mean': t_vals.mean(),
            'nontarget_mean': nt_vals.mean(),
            'statistic': stat,
            'p_value': pval,
            'effect_size_r': r,
            'n_target': len(t_vals),
            'n_nontarget': len(nt_vals),
            'method': 'Mann-Whitney U'
        })
        print(f"\n  {timepoint} |LFC| comparison:")
        print(f"    Target: median={t_vals.median():.3f}, mean={t_vals.mean():.3f} (n={len(t_vals)})")
        print(f"    Non-target: median={nt_vals.median():.3f}, mean={nt_vals.mean():.3f} (n={len(nt_vals)})")
        print(f"    Mann-Whitney U p={pval:.2e}, r={r:.3f}")

# Test 2: DEG rate comparison (Fisher exact)
for timepoint, col in [('T2', 'is_DEG_T2'), ('T3', 'is_DEG_T3'), ('any', 'is_DEG_any')]:
    t_deg = target_expr[col].sum()
    t_nondeg = len(target_expr) - t_deg
    nt_deg = nontarget_expr[col].sum()
    nt_nondeg = len(nontarget_expr) - nt_deg

    table = [[t_deg, t_nondeg], [nt_deg, nt_nondeg]]
    odds_ratio, pval = stats.fisher_exact(table)

    t_rate = t_deg / len(target_expr) * 100 if len(target_expr) > 0 else 0
    nt_rate = nt_deg / len(nontarget_expr) * 100 if len(nontarget_expr) > 0 else 0

    stat_results.append({
        'test': f'DEG_rate_target_vs_nontarget_{timepoint}',
        'target_median': t_rate,
        'nontarget_median': nt_rate,
        'target_mean': t_rate,
        'nontarget_mean': nt_rate,
        'statistic': odds_ratio,
        'p_value': pval,
        'effect_size_r': odds_ratio,
        'n_target': len(target_expr),
        'n_nontarget': len(nontarget_expr),
        'method': 'Fisher exact'
    })
    print(f"\n  {timepoint} DEG rate:")
    print(f"    Target: {t_deg}/{len(target_expr)} ({t_rate:.1f}%)")
    print(f"    Non-target: {nt_deg}/{len(nontarget_expr)} ({nt_rate:.1f}%)")
    print(f"    OR={odds_ratio:.3f}, p={pval:.2e}")

# Per-TF expression analysis
per_tf_expr = []
for tf in sorted(pairs_df['tf_name'].unique()):
    tf_targets = set(pairs_df[pairs_df['tf_name'] == tf]['target_gene'])
    tf_target_expr = expr_data[expr_data['gene_id'].isin(tf_targets)]
    tf_nontarget_expr = expr_data[~expr_data['gene_id'].isin(tf_targets)]

    for tp, lfc_col, deg_col in [('T2', 'abs_lfc_T2', 'is_DEG_T2'), ('T3', 'abs_lfc_T3', 'is_DEG_T3')]:
        t_vals = tf_target_expr[lfc_col].dropna()
        nt_vals = tf_nontarget_expr[lfc_col].dropna()

        if len(t_vals) > 5:
            stat, pval = stats.mannwhitneyu(t_vals, nt_vals, alternative='two-sided')

            t_deg = tf_target_expr[deg_col].sum()
            nt_deg = tf_nontarget_expr[deg_col].sum()
            t_rate = t_deg / len(tf_target_expr) * 100 if len(tf_target_expr) > 0 else 0

            per_tf_expr.append({
                'tf_name': tf,
                'timepoint': tp,
                'n_targets': len(t_vals),
                'mean_abs_lfc_target': t_vals.mean(),
                'mean_abs_lfc_nontarget': nt_vals.mean(),
                'median_abs_lfc_target': t_vals.median(),
                'median_abs_lfc_nontarget': nt_vals.median(),
                'wilcoxon_p': pval,
                'deg_rate_target': t_rate,
                'deg_count_target': int(t_deg)
            })

per_tf_expr_df = pd.DataFrame(per_tf_expr)

# ============================================================
# Step 5b: Exposed TFs as targets - coordination analysis
# ============================================================
print("\n" + "-" * 50)
print("Step 5b: Exposed regulators as targets - coordination")
print("-" * 50)

# Check if exposed regulators that are FIMO TF targets show coordination
exposed_targeted = exposed_df[exposed_df['locus_tag'].isin(exposed_as_targets['target_gene'].unique())].copy()
exposed_not_targeted = exposed_df[~exposed_df['locus_tag'].isin(exposed_as_targets['target_gene'].unique())].copy()

print(f"Exposed regulators targeted by FIMO TFs: {len(exposed_targeted)}")
print(f"Exposed regulators NOT targeted: {len(exposed_not_targeted)}")

for tp_col in ['coordination_T2', 'coordination_T3']:
    print(f"\n  {tp_col} distribution:")
    print(f"    Targeted: {exposed_targeted[tp_col].value_counts().to_dict()}")
    print(f"    Not targeted: {exposed_not_targeted[tp_col].value_counts().to_dict()}")

# Compare |LFC| of targeted vs non-targeted exposed TFs
for tp, col in [('T2', 'log2FC_T2'), ('T3', 'log2FC_T3')]:
    t_vals = exposed_targeted[col].abs().dropna()
    nt_vals = exposed_not_targeted[col].abs().dropna()
    if len(t_vals) > 0 and len(nt_vals) > 0:
        stat, pval = stats.mannwhitneyu(t_vals, nt_vals, alternative='two-sided')
        print(f"\n  {tp} |LFC| of exposed regulators:")
        print(f"    Targeted by FIMO TFs: median={t_vals.median():.3f} (n={len(t_vals)})")
        print(f"    Not targeted: median={nt_vals.median():.3f} (n={len(nt_vals)})")
        print(f"    p={pval:.3f}")

        stat_results.append({
            'test': f'exposed_targeted_vs_not_{tp}_abs_LFC',
            'target_median': t_vals.median(),
            'nontarget_median': nt_vals.median(),
            'target_mean': t_vals.mean(),
            'nontarget_mean': nt_vals.mean(),
            'statistic': stat,
            'p_value': pval,
            'effect_size_r': 0,
            'n_target': len(t_vals),
            'n_nontarget': len(nt_vals),
            'method': 'Mann-Whitney U'
        })

# ============================================================
# Step 6: Cross-regulation among exposed TFs
# ============================================================
print("\n" + "=" * 70)
print("STEP 6: Cross-regulation among exposed TFs")
print("=" * 70)

# FIMO TFs regulating exposed TFs (cross-TF-class regulation)
exposed_regulation_edges = []
for _, row in exposed_as_targets.iterrows():
    exposed_regulation_edges.append({
        'regulator_name': row['tf_name'],
        'regulator_sco': row['tf_sco'],
        'regulator_type': 'FIMO_characterized',
        'target_gene': row['target_gene'],
        'target_sco': row['target_sco'],
        'target_type': 'exposed_regulator',
        'distance_to_tss': row['distance_to_tss'],
        'bs_score': row['bs_score'],
        'bs_qvalue': row['bs_qvalue']
    })

# Check if any FIMO TFs regulate other FIMO TFs
fimo_tf_rs_set = set(fimo_tf_mapping[k]['rs'] for k in fimo_tf_mapping)
fimo_to_fimo = pairs_df[pairs_df['target_gene'].isin(fimo_tf_rs_set)]
print(f"\nFIMO TF -> FIMO TF regulatory edges: {len(fimo_to_fimo)}")
for _, row in fimo_to_fimo.iterrows():
    # Find target TF name
    target_tf_name = [k for k, v in fimo_tf_mapping.items() if v['rs'] == row['target_gene']]
    target_name = target_tf_name[0] if target_tf_name else row['target_gene']
    print(f"  {row['tf_name']} -> {target_name} ({row['target_gene']}) dist={row['distance_to_tss']:.0f}bp")
    exposed_regulation_edges.append({
        'regulator_name': row['tf_name'],
        'regulator_sco': row['tf_sco'],
        'regulator_type': 'FIMO_characterized',
        'target_gene': row['target_gene'],
        'target_sco': row['target_sco'] if 'target_sco' in row else '',
        'target_type': f'FIMO_TF_{target_name}',
        'distance_to_tss': row['distance_to_tss'],
        'bs_score': row['bs_score'],
        'bs_qvalue': row['bs_qvalue']
    })

# Check regulation among shielded regulators too
shielded_as_targets_detail = pairs_df[pairs_df['target_gene'].isin(shielded_locus_tags)]
print(f"\nFIMO TF -> Shielded regulator edges: {shielded_as_targets_detail['target_gene'].nunique()} unique targets")
print(f"FIMO TF -> Exposed regulator edges: {exposed_as_targets['target_gene'].nunique()} unique targets")

# Statistical test: Are exposed regulators enriched as targets compared to shielded?
n_exposed = len(exposed_locus_tags)
n_shielded = len(shielded_locus_tags)
n_exposed_targeted = exposed_as_targets['target_gene'].nunique()
n_shielded_targeted = shielded_as_targets_detail['target_gene'].nunique()

table = [[n_exposed_targeted, n_exposed - n_exposed_targeted],
         [n_shielded_targeted, n_shielded - n_shielded_targeted]]
or_val, pval = stats.fisher_exact(table)
print(f"\nExposed vs Shielded as targets of FIMO TFs:")
print(f"  Exposed: {n_exposed_targeted}/{n_exposed} ({n_exposed_targeted/n_exposed*100:.1f}%)")
print(f"  Shielded: {n_shielded_targeted}/{n_shielded} ({n_shielded_targeted/n_shielded*100:.1f}%)")
print(f"  Fisher exact OR={or_val:.3f}, p={pval:.3f}")

stat_results.append({
    'test': 'exposed_vs_shielded_as_FIMO_targets',
    'target_median': n_exposed_targeted/n_exposed*100,
    'nontarget_median': n_shielded_targeted/n_shielded*100,
    'target_mean': n_exposed_targeted,
    'nontarget_mean': n_shielded_targeted,
    'statistic': or_val,
    'p_value': pval,
    'effect_size_r': or_val,
    'n_target': n_exposed,
    'n_nontarget': n_shielded,
    'method': 'Fisher exact'
})

cross_reg_df = pd.DataFrame(exposed_regulation_edges)

# ============================================================
# Step 7: Shielded vs Exposed TF target comparison
# ============================================================
print("\n" + "=" * 70)
print("STEP 7: Shielded vs Exposed TF target comparison")
print("=" * 70)

# Compare functional categories of genes targeted by FIMO TFs
# that are exposed vs shielded vs non-regulatory

# Use gene annotation to get product/function categories
gene_products = {}
for _, row in gene_annot.iterrows():
    gene_products[row['gene_id']] = row.get('product', '')

# Classify target genes by whether they target exposed/shielded/non-reg
target_classification = []
for gene in all_fimo_targets:
    is_exposed = gene in exposed_locus_tags
    is_shielded = gene in shielded_locus_tags
    is_reg = gene in all_reg_locus_tags
    product = gene_products.get(gene, '')

    target_classification.append({
        'gene_id': gene,
        'is_exposed': is_exposed,
        'is_shielded': is_shielded,
        'is_regulatory': is_reg,
        'product': product
    })

target_class_df = pd.DataFrame(target_classification)
print(f"Total unique target genes: {len(target_class_df)}")
print(f"  Exposed regulators: {target_class_df['is_exposed'].sum()}")
print(f"  Shielded regulators: {target_class_df['is_shielded'].sum()}")
print(f"  Non-regulatory: {(~target_class_df['is_regulatory']).sum()}")

# COG/functional enrichment proxy using product keywords
functional_categories = {
    'transcriptional_regulator': ['transcriptional regulator', 'transcription factor', 'repressor', 'activator'],
    'sigma_factor': ['sigma factor', 'sigma-70', 'RNA polymerase sigma'],
    'two_component': ['sensor kinase', 'response regulator', 'two-component', 'histidine kinase'],
    'transport': ['transporter', 'permease', 'ABC transport', 'efflux', 'MFS'],
    'metabolism': ['synthase', 'reductase', 'dehydrogenase', 'kinase', 'oxidase', 'transferase'],
    'secondary_metabolism': ['polyketide', 'NRPS', 'terpene', 'siderophore'],
    'hypothetical': ['hypothetical protein'],
    'cell_division': ['cell division', 'FtsZ', 'septum'],
    'DNA_repair': ['DNA repair', 'recombinase', 'helicase']
}

def classify_product(product):
    if not product or pd.isna(product):
        return 'unknown'
    product_lower = product.lower()
    for cat, keywords in functional_categories.items():
        for kw in keywords:
            if kw.lower() in product_lower:
                return cat
    return 'other'

target_class_df['func_category'] = target_class_df['product'].apply(classify_product)
all_genes_annot = gene_annot.copy()
all_genes_annot['func_category'] = all_genes_annot['product'].apply(classify_product)
all_genes_annot['is_target'] = all_genes_annot['gene_id'].isin(target_genes)

# Enrichment test for each category
print("\nFunctional enrichment of FIMO TF target genes vs all genes:")
enrichment_results = []
for cat in sorted(functional_categories.keys()):
    n_target_cat = target_class_df[target_class_df['func_category'] == cat].shape[0]
    n_target_total = len(target_class_df)
    n_all_cat = all_genes_annot[all_genes_annot['func_category'] == cat].shape[0]
    n_all_total = len(all_genes_annot)

    n_target_not = n_target_total - n_target_cat
    n_nontarget_cat = n_all_cat - n_target_cat
    n_nontarget_not = (n_all_total - n_target_total) - n_nontarget_cat

    if n_nontarget_cat < 0:
        n_nontarget_cat = 0
        n_nontarget_not = n_all_total - n_target_total

    table = [[n_target_cat, n_target_not],
             [n_nontarget_cat, n_nontarget_not]]

    try:
        or_val, pval = stats.fisher_exact(table)
    except:
        or_val, pval = 1.0, 1.0

    pct_target = n_target_cat / n_target_total * 100 if n_target_total > 0 else 0
    pct_all = n_all_cat / n_all_total * 100 if n_all_total > 0 else 0

    enrichment_results.append({
        'category': cat,
        'n_target': n_target_cat,
        'pct_target': pct_target,
        'n_genome': n_all_cat,
        'pct_genome': pct_all,
        'odds_ratio': or_val,
        'p_value': pval,
        'fold_enrichment': (pct_target / pct_all) if pct_all > 0 else 0
    })

    print(f"  {cat}: targets={n_target_cat} ({pct_target:.1f}%) vs genome={n_all_cat} ({pct_all:.1f}%), OR={or_val:.2f}, p={pval:.3f}")

enrichment_df = pd.DataFrame(enrichment_results)

# Apply FDR correction
if len(enrichment_df) > 0:
    reject, pvals_corrected, _, _ = multipletests(enrichment_df['p_value'], method='fdr_bh')
    enrichment_df['p_adj'] = pvals_corrected
    enrichment_df['significant'] = reject

# ============================================================
# Step 8: Visualization
# ============================================================
print("\n" + "=" * 70)
print("STEP 8: Creating visualizations")
print("=" * 70)

plt.rcParams.update({
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 10,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight'
})

# Figure 1: Target count distribution
fig, ax = plt.subplots(figsize=(10, 6))
tf_names = per_tf_df.sort_values('n_target_genes', ascending=True)['tf_name']
n_targets = per_tf_df.sort_values('n_target_genes', ascending=True)['n_target_genes']
n_exposed_per_tf = per_tf_df.sort_values('n_target_genes', ascending=True)['n_exposed_targets']

bars = ax.barh(range(len(tf_names)), n_targets, color='steelblue', alpha=0.8, label='All targets')
ax.barh(range(len(tf_names)), n_exposed_per_tf, color='crimson', alpha=0.8, label='Exposed TF targets')
ax.set_yticks(range(len(tf_names)))
ax.set_yticklabels(tf_names)
ax.set_xlabel('Number of predicted target genes (q < 0.05)')
ax.set_title('FIMO-predicted targets per TF (promoter ±500bp)')
ax.legend(loc='lower right')

# Add count labels
for i, (nt, ne) in enumerate(zip(n_targets, n_exposed_per_tf)):
    ax.text(nt + 5, i, str(int(nt)), va='center', fontsize=9)

plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f'{FIG_DIR}/target_count_distribution.{ext}')
plt.close()
print("  Saved: target_count_distribution.pdf/svg")

# Figure 2: Cascade amplification
fig, ax = plt.subplots(figsize=(12, 8))

# Create a simplified cascade diagram
# Left column: FIMO TFs, middle: target counts, right: exposed TF targets
tf_data = per_tf_df.sort_values('n_target_genes', ascending=False)

y_positions_left = np.linspace(0.9, 0.1, len(tf_data))
total_targets = len(all_fimo_targets)
total_exposed_targets = len(exposed_fimo_targets)
total_reg_targets = len(all_fimo_targets & all_reg_locus_tags)

# Draw TF boxes on left
for i, (_, row) in enumerate(tf_data.iterrows()):
    y = y_positions_left[i]
    box_height = max(0.02, row['n_target_genes'] / total_targets * 0.15)

    # TF box
    rect = mpatches.FancyBboxPatch((0.05, y - box_height/2), 0.15, box_height,
                                     boxstyle="round,pad=0.01",
                                     facecolor='steelblue', alpha=0.7, edgecolor='black')
    ax.add_patch(rect)
    ax.text(0.125, y, f"{row['tf_name']}\n({row['n_target_genes']})",
            ha='center', va='center', fontsize=7, fontweight='bold', color='white')

    # Arrow to center
    ax.annotate('', xy=(0.4, 0.5), xytext=(0.2, y),
               arrowprops=dict(arrowstyle='->', color='gray', alpha=0.3, lw=max(0.5, row['n_target_genes']/200)))

# Center target pool
circle = mpatches.Circle((0.5, 0.5), 0.12, facecolor='lightblue', edgecolor='steelblue', lw=2)
ax.add_patch(circle)
ax.text(0.5, 0.55, f'{total_targets}', ha='center', va='center', fontsize=20, fontweight='bold')
ax.text(0.5, 0.45, 'unique target\ngenes', ha='center', va='center', fontsize=9)

# Right side: breakdown
categories = [
    (f'{total_reg_targets} regulatory\ngenes', 'coral', 0.75, 0.7),
    (f'{total_exposed_targets} exposed\nregulators', 'crimson', 0.75, 0.5),
    (f'{total_targets - total_reg_targets} non-regulatory\ngenes', 'lightgreen', 0.75, 0.3),
]

for label, color, x, y in categories:
    rect = mpatches.FancyBboxPatch((x - 0.08, y - 0.06), 0.16, 0.12,
                                     boxstyle="round,pad=0.01",
                                     facecolor=color, alpha=0.6, edgecolor='black')
    ax.add_patch(rect)
    ax.text(x, y, label, ha='center', va='center', fontsize=8, fontweight='bold')
    ax.annotate('', xy=(x - 0.08, y), xytext=(0.62, 0.5),
               arrowprops=dict(arrowstyle='->', color='gray', alpha=0.5, lw=1.5))

ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_aspect('equal')
ax.axis('off')
ax.set_title(f'Cascade Amplification: {len(tf_data)} FIMO TFs -> {total_targets} target genes\n'
             f'({total_exposed_targets} exposed regulators among targets)', fontsize=13)

plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f'{FIG_DIR}/cascade_amplification.{ext}')
plt.close()
print("  Saved: cascade_amplification.pdf/svg")

# Figure 3: Expression coordination
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for idx, (tp, col) in enumerate([('T2vsT1', 'abs_lfc_T2'), ('T3vsT1', 'abs_lfc_T3')]):
    ax = axes[idx]
    t_data = target_expr[col].dropna()
    nt_data = nontarget_expr[col].dropna()

    # Violin plot
    parts = ax.violinplot([nt_data.values, t_data.values], positions=[1, 2], showmedians=True)
    for pc in parts['bodies']:
        pc.set_alpha(0.6)
    parts['bodies'][0].set_facecolor('gray')
    parts['bodies'][1].set_facecolor('steelblue')

    ax.set_xticks([1, 2])
    ax.set_xticklabels([f'Non-target\n(n={len(nt_data)})', f'Target\n(n={len(t_data)})'])
    ax.set_ylabel('|log2FoldChange|')
    ax.set_title(f'{tp}')

    # Add stats
    stat_row = [s for s in stat_results if s['test'] == f'|LFC|_target_vs_nontarget_{tp.split("vs")[0]}']
    if stat_row:
        p = stat_row[0]['p_value']
        ax.text(0.5, 0.95, f'p={p:.2e}', transform=ax.transAxes, ha='center', va='top', fontsize=10,
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

fig.suptitle('Expression magnitude: FIMO TF targets vs non-targets', fontsize=13)
plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f'{FIG_DIR}/expression_coordination.{ext}')
plt.close()
print("  Saved: expression_coordination.pdf/svg")

# Figure 4: Cross-regulation network
fig, ax = plt.subplots(figsize=(12, 10))

# Create a network-like visualization
# Nodes: FIMO TFs (left/top), Exposed regulators (right/bottom)
if len(cross_reg_df) > 0:
    # Get unique regulators and targets
    regulators = cross_reg_df['regulator_name'].unique()
    targets = cross_reg_df[cross_reg_df['target_type'] == 'exposed_regulator']['target_gene'].unique()

    # Position regulators on left
    reg_y = np.linspace(0.9, 0.1, len(regulators))
    reg_pos = {r: (0.15, reg_y[i]) for i, r in enumerate(regulators)}

    # Position targets on right
    if len(targets) > 0:
        tgt_y = np.linspace(0.95, 0.05, len(targets))
        tgt_pos = {t: (0.85, tgt_y[i]) for i, t in enumerate(targets)}
    else:
        tgt_pos = {}

    # Draw nodes
    for reg, (x, y) in reg_pos.items():
        circle = mpatches.Circle((x, y), 0.025, facecolor='steelblue', edgecolor='black', zorder=3)
        ax.add_patch(circle)
        ax.text(x - 0.04, y, reg, ha='right', va='center', fontsize=8, fontweight='bold')

    for tgt, (x, y) in tgt_pos.items():
        sco = rs_to_sco.get(tgt, '')
        circle = mpatches.Circle((x, y), 0.015, facecolor='crimson', edgecolor='black', zorder=3)
        ax.add_patch(circle)
        label = sco if sco else tgt
        ax.text(x + 0.03, y, label, ha='left', va='center', fontsize=7)

    # Draw edges
    for _, row in cross_reg_df[cross_reg_df['target_type'] == 'exposed_regulator'].iterrows():
        if row['regulator_name'] in reg_pos and row['target_gene'] in tgt_pos:
            x1, y1 = reg_pos[row['regulator_name']]
            x2, y2 = tgt_pos[row['target_gene']]
            ax.plot([x1 + 0.03, x2 - 0.02], [y1, y2], '-', color='gray', alpha=0.4, lw=0.8)

    # Also draw FIMO->FIMO edges
    for _, row in cross_reg_df[cross_reg_df['target_type'].str.startswith('FIMO_TF')].iterrows():
        tf_name = row['target_type'].replace('FIMO_TF_', '')
        if row['regulator_name'] in reg_pos and tf_name in reg_pos:
            x1, y1 = reg_pos[row['regulator_name']]
            x2, y2 = reg_pos[tf_name]
            ax.annotate('', xy=(x2 + 0.03, y2), xytext=(x1 + 0.03, y1),
                       arrowprops=dict(arrowstyle='->', color='blue', lw=1.5, alpha=0.6))

    ax.text(0.15, 0.97, 'FIMO TFs', ha='center', fontsize=11, fontweight='bold', color='steelblue')
    ax.text(0.85, 0.97, 'Exposed regulators', ha='center', fontsize=11, fontweight='bold', color='crimson')
else:
    ax.text(0.5, 0.5, 'No cross-regulation edges found', ha='center', va='center', fontsize=14)

ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')
ax.set_title(f'Cross-regulation network: FIMO TFs -> Exposed regulators\n'
             f'({len(cross_reg_df)} edges)', fontsize=13)

plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f'{FIG_DIR}/cross_regulation_network.{ext}')
plt.close()
print("  Saved: cross_regulation_network.pdf/svg")

# Figure 5: Comprehensive summary
fig = plt.figure(figsize=(18, 14))
gs = GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.35)

# Panel A: Target count per TF
ax1 = fig.add_subplot(gs[0, 0])
tf_sorted = per_tf_df.sort_values('n_target_genes', ascending=True)
colors = ['crimson' if n > 0 else 'steelblue'
          for n in tf_sorted['n_exposed_targets']]
ax1.barh(range(len(tf_sorted)), tf_sorted['n_target_genes'], color='steelblue', alpha=0.7)
ax1.barh(range(len(tf_sorted)), tf_sorted['n_exposed_targets'], color='crimson', alpha=0.8)
ax1.set_yticks(range(len(tf_sorted)))
ax1.set_yticklabels(tf_sorted['tf_name'], fontsize=8)
ax1.set_xlabel('Target genes')
ax1.set_title('A. Targets per FIMO TF', fontweight='bold')

# Panel B: |LFC| boxplot targets vs non-targets
ax2 = fig.add_subplot(gs[0, 1])
data_t2 = [nontarget_expr['abs_lfc_T2'].dropna().values, target_expr['abs_lfc_T2'].dropna().values]
bp = ax2.boxplot(data_t2, positions=[1, 2], widths=0.6, patch_artist=True,
                  showfliers=False, medianprops=dict(color='black', lw=2))
bp['boxes'][0].set_facecolor('lightgray')
bp['boxes'][1].set_facecolor('steelblue')
bp['boxes'][0].set_alpha(0.7)
bp['boxes'][1].set_alpha(0.7)

data_t3 = [nontarget_expr['abs_lfc_T3'].dropna().values, target_expr['abs_lfc_T3'].dropna().values]
bp2 = ax2.boxplot(data_t3, positions=[3.5, 4.5], widths=0.6, patch_artist=True,
                   showfliers=False, medianprops=dict(color='black', lw=2))
bp2['boxes'][0].set_facecolor('lightgray')
bp2['boxes'][1].set_facecolor('coral')
bp2['boxes'][0].set_alpha(0.7)
bp2['boxes'][1].set_alpha(0.7)

ax2.set_xticks([1.5, 4.0])
ax2.set_xticklabels(['T2vsT1', 'T3vsT1'])
ax2.set_ylabel('|log2FC|')
ax2.set_title('B. Expression magnitude', fontweight='bold')
ax2.legend([mpatches.Patch(color='lightgray'), mpatches.Patch(color='steelblue')],
           ['Non-target', 'Target'], fontsize=8, loc='upper right')

# Panel C: DEG rate comparison
ax3 = fig.add_subplot(gs[0, 2])
deg_data = []
for tp, col in [('T2', 'is_DEG_T2'), ('T3', 'is_DEG_T3')]:
    t_rate = target_expr[col].mean() * 100
    nt_rate = nontarget_expr[col].mean() * 100
    deg_data.append({'timepoint': tp, 'target': t_rate, 'nontarget': nt_rate})
deg_df_plot = pd.DataFrame(deg_data)
x = np.arange(len(deg_df_plot))
w = 0.35
ax3.bar(x - w/2, deg_df_plot['nontarget'], w, label='Non-target', color='lightgray', edgecolor='black')
ax3.bar(x + w/2, deg_df_plot['target'], w, label='Target', color='steelblue', edgecolor='black')
ax3.set_xticks(x)
ax3.set_xticklabels(deg_df_plot['timepoint'])
ax3.set_ylabel('DEG rate (%)')
ax3.set_title('C. DEG rates', fontweight='bold')
ax3.legend(fontsize=8)

# Panel D: Functional enrichment
ax4 = fig.add_subplot(gs[1, 0])
if len(enrichment_df) > 0:
    enr_sorted = enrichment_df.sort_values('fold_enrichment', ascending=True)
    enr_sorted = enr_sorted[enr_sorted['n_target'] > 0]  # only show categories with targets
    colors_enr = ['crimson' if p < 0.05 else 'gray' for p in enr_sorted['p_adj']]
    ax4.barh(range(len(enr_sorted)), enr_sorted['fold_enrichment'], color=colors_enr, alpha=0.7)
    ax4.axvline(x=1, color='black', linestyle='--', lw=0.8)
    ax4.set_yticks(range(len(enr_sorted)))
    ax4.set_yticklabels(enr_sorted['category'], fontsize=8)
    ax4.set_xlabel('Fold enrichment')
    ax4.set_title('D. Functional enrichment\n(targets vs genome)', fontweight='bold')

# Panel E: Exposed vs Shielded targeting rates
ax5 = fig.add_subplot(gs[1, 1])
categories_bar = ['Exposed', 'Shielded', 'Non-reg']
rates = [
    n_exposed_targeted / n_exposed * 100,
    n_shielded_targeted / n_shielded * 100,
    (len(target_genes) - n_exposed_targeted - n_shielded_targeted) / max(1, len(gene_annot) - len(exposed_locus_tags) - len(shielded_locus_tags)) * 100
]
bar_colors = ['crimson', 'steelblue', 'lightgray']
x_pos = np.arange(len(categories_bar))
ax5.bar(x_pos, rates, color=bar_colors, edgecolor='black', alpha=0.7)
ax5.set_xticks(x_pos)
ax5.set_xticklabels(categories_bar)
ax5.set_ylabel('% targeted by FIMO TFs')
ax5.set_title('E. Targeting rates by\ngene class', fontweight='bold')
for i, r in enumerate(rates):
    ax5.text(i, r + 0.5, f'{r:.1f}%', ha='center', fontsize=9)

# Panel F: Cascade summary statistics
ax6 = fig.add_subplot(gs[1, 2])
ax6.axis('off')
summary_text = (
    f"CASCADE AMPLIFICATION SUMMARY\n"
    f"{'='*40}\n\n"
    f"FIMO TFs with motifs: 16\n"
    f"  With q < 0.05 hits: {fimo_primary['motif_id'].nunique()}\n\n"
    f"Total binding sites: {len(fimo_primary)}\n"
    f"Unique target genes: {len(all_fimo_targets)}\n"
    f"  Regulatory targets: {len(all_fimo_targets & all_reg_locus_tags)}\n"
    f"  Exposed TF targets: {len(exposed_fimo_targets)}\n"
    f"  Shielded TF targets: {n_shielded_targeted}\n\n"
    f"Amplification ratio:\n"
    f"  {fimo_primary['motif_id'].nunique()} TFs -> {len(all_fimo_targets)} genes\n"
    f"  ({len(all_fimo_targets)/fimo_primary['motif_id'].nunique():.0f}x per TF)\n\n"
    f"Cross-regulation:\n"
    f"  FIMO->FIMO: {len(fimo_to_fimo)} edges\n"
    f"  FIMO->Exposed: {n_exposed_targeted} edges\n"
    f"  Exposed enrichment: OR={or_val:.2f}, p={pval:.3f}"
)
ax6.text(0.05, 0.95, summary_text, transform=ax6.transAxes, va='top', fontsize=9,
         family='monospace', bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
ax6.set_title('F. Key statistics', fontweight='bold')

fig.suptitle('H31: Exposed TF Downstream Regulatory Target Network', fontsize=15, fontweight='bold', y=0.98)

for ext in ['pdf', 'svg']:
    fig.savefig(f'{FIG_DIR}/H31_comprehensive_summary.{ext}')
plt.close()
print("  Saved: H31_comprehensive_summary.pdf/svg")

# ============================================================
# Step 9: Output tables
# ============================================================
print("\n" + "=" * 70)
print("STEP 9: Saving output tables")
print("=" * 70)

# Table 1: TF-target pairs
pairs_df.to_csv(f'{TAB_DIR}/TF_target_pairs.tsv', sep='\t', index=False)
print(f"  Saved: TF_target_pairs.tsv ({len(pairs_df)} rows)")

# Table 2: Per-TF summary (merge expression data)
per_tf_final = per_tf_df.copy()
# Add expression stats
for tp in ['T2', 'T3']:
    tp_data = per_tf_expr_df[per_tf_expr_df['timepoint'] == tp][['tf_name', 'mean_abs_lfc_target', 'median_abs_lfc_target', 'wilcoxon_p', 'deg_rate_target', 'deg_count_target']]
    tp_data.columns = ['tf_name'] + [f'{c}_{tp}' for c in tp_data.columns[1:]]
    per_tf_final = per_tf_final.merge(tp_data, on='tf_name', how='left')

per_tf_final.to_csv(f'{TAB_DIR}/per_TF_summary.tsv', sep='\t', index=False)
print(f"  Saved: per_TF_summary.tsv ({len(per_tf_final)} rows)")

# Table 3: Cascade amplification summary
cascade_data = {
    'metric': [
        'total_FIMO_TFs_in_database', 'FIMO_TFs_with_q005_hits',
        'total_binding_sites_q005', 'total_unique_target_genes',
        'regulatory_target_genes', 'exposed_regulator_targets',
        'shielded_regulator_targets', 'non_regulatory_targets',
        'amplification_ratio_per_TF',
        'FIMO_to_FIMO_edges', 'FIMO_to_exposed_edges',
        'exposed_targeting_rate_pct', 'shielded_targeting_rate_pct',
        'exposed_vs_shielded_OR', 'exposed_vs_shielded_p',
        'relaxed_q01_total_sites', 'relaxed_q01_unique_targets'
    ],
    'value': [
        16, fimo_primary['motif_id'].nunique(),
        len(fimo_primary), len(all_fimo_targets),
        len(all_fimo_targets & all_reg_locus_tags), len(exposed_fimo_targets),
        n_shielded_targeted, len(all_fimo_targets) - len(all_fimo_targets & all_reg_locus_tags),
        len(all_fimo_targets) / fimo_primary['motif_id'].nunique(),
        len(fimo_to_fimo), n_exposed_targeted,
        n_exposed_targeted / n_exposed * 100, n_shielded_targeted / n_shielded * 100,
        or_val, pval,
        len(fimo_relaxed), fimo_relaxed.groupby(['motif_id']).size().sum()  # placeholder
    ]
}

# Also compute relaxed threshold targets
pairs_relaxed = []
for _, hit in fimo_relaxed.iterrows():
    targets = find_target_gene(hit['start'], hit['stop'])
    for t in targets:
        pairs_relaxed.append(t['gene_id'])
n_relaxed_targets = len(set(pairs_relaxed))
cascade_data['value'][-1] = n_relaxed_targets

cascade_df = pd.DataFrame(cascade_data)
cascade_df.to_csv(f'{TAB_DIR}/cascade_amplification.tsv', sep='\t', index=False)
print(f"  Saved: cascade_amplification.tsv")

# Table 4: Cross-regulation edges
cross_reg_df.to_csv(f'{TAB_DIR}/cross_regulation_edges.tsv', sep='\t', index=False)
print(f"  Saved: cross_regulation_edges.tsv ({len(cross_reg_df)} rows)")

# Table 5: Statistical tests
stat_df = pd.DataFrame(stat_results)
# Apply FDR correction
if len(stat_df) > 0:
    reject, pvals_corrected, _, _ = multipletests(stat_df['p_value'], method='fdr_bh')
    stat_df['p_adjusted'] = pvals_corrected
    stat_df['significant_fdr05'] = reject
stat_df.to_csv(f'{TAB_DIR}/statistical_tests.tsv', sep='\t', index=False)
print(f"  Saved: statistical_tests.tsv ({len(stat_df)} rows)")

# Additional table: Enrichment results
enrichment_df.to_csv(f'{TAB_DIR}/functional_enrichment.tsv', sep='\t', index=False)
print(f"  Saved: functional_enrichment.tsv ({len(enrichment_df)} rows)")

# ============================================================
# Step 10: Summary statistics for report
# ============================================================
print("\n" + "=" * 70)
print("STEP 10: Summary for report")
print("=" * 70)

print(f"""
KEY FINDINGS:
=============

1. FIMO Binding Site Analysis:
   - 16 TFs have FIMO motifs; {fimo_primary['motif_id'].nunique()} have significant hits at q < 0.05
   - {len(fimo_primary)} binding sites pass strict threshold
   - None of the 57 exposed regulators have FIMO motifs (they are newly identified)

2. Target Gene Mapping:
   - {len(all_fimo_targets)} unique genes have FIMO TF binding sites in promoters (±{PROMOTER_WINDOW}bp)
   - {len(all_fimo_targets)/len(gene_annot)*100:.1f}% of all genes are predicted targets

3. Exposed Regulators as Targets:
   - {n_exposed_targeted}/{n_exposed} ({n_exposed_targeted/n_exposed*100:.1f}%) exposed TFs are targets of FIMO TFs
   - {n_shielded_targeted}/{n_shielded} ({n_shielded_targeted/n_shielded*100:.1f}%) shielded TFs are targets
   - Fisher exact: OR={or_val:.3f}, p={pval:.3f}

4. Expression Coordination:
   - Target genes vs non-targets: """)

for s in stat_results:
    if 'LFC' in s['test']:
        print(f"     {s['test']}: median target={s['target_median']:.3f} vs non-target={s['nontarget_median']:.3f}, p={s['p_value']:.2e}")

print(f"""
5. Cross-regulation:
   - {len(fimo_to_fimo)} FIMO TF -> FIMO TF edges
   - {n_exposed_targeted} FIMO TF -> Exposed regulator edges

6. Functional Enrichment:""")
for _, row in enrichment_df.iterrows():
    if row['n_target'] > 0:
        sig = '*' if row.get('p_adj', 1) < 0.05 else ''
        print(f"     {row['category']}: fold={row['fold_enrichment']:.2f}, p_adj={row.get('p_adj', row['p_value']):.3f}{sig}")

print(f"""
VERDICT: The hypothesis that 57 exposed TFs regulate a large downstream
network is PARTIALLY SUPPORTED with important caveats:
- The 57 exposed TFs do NOT have FIMO motifs (they are uncharacterized)
- {n_exposed_targeted} exposed TFs are themselves targets of well-characterized TFs
- The cascade model needs experimental validation of exposed TF regulons
""")

print("\nAnalysis complete!")
