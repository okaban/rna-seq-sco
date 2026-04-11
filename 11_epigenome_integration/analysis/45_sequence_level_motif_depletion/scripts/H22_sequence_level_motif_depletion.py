#!/usr/bin/env python3
"""
H22: Sequence-Level Motif Depletion Near Regulatory Genes

Tests whether R-M recognition motifs (TGGCCGGC, AAGCCCG) are depleted
at the DNA SEQUENCE level near regulatory genes, independent of methylation status.

This distinguishes:
  1. Sequence-level avoidance (evolutionary counter-selection)
  2. Protein occupancy model (TF/RNAP blocking MTase access)
"""

import re
import numpy as np
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from scipy import stats
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Configuration
# ============================================================
BASE_DIR = Path('/Users/okaban/bioinfo/rna-seq')
ANALYSIS_DIR = BASE_DIR / '11_epigenome_integration/analysis/45_sequence_level_motif_depletion'
FIGURES_DIR = ANALYSIS_DIR / 'figures'
TABLES_DIR = ANALYSIS_DIR / 'tables'

GENOME_FASTA = Path('/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/GCF_000203835.1_ASM20383v1_genomic.fna')
GENE_ANNOTATION = BASE_DIR / '05_annotation/analysis/05_annotation_260128_v1/tables/gene_annotation_basic.tsv'
REGULATORY_GENES = BASE_DIR / '11_epigenome_integration/analysis/29_genomewide_TF_screen/tables/all_regulatory_genes.tsv'

# H15 data for methylation-level comparison
H15_ENRICHMENT = BASE_DIR / '11_epigenome_integration/analysis/38_cross_motif_regulatory_avoidance/tables/cross_motif_functional_enrichment.tsv'
# H20 data for stratified comparison
H20_STRATIFIED = BASE_DIR / '11_epigenome_integration/analysis/43_regulatory_avoidance_geographic_test/tables/stratified_enrichment.tsv'

# Genomic regions
CHROMOSOME_LENGTH = 8_667_507
ARM_BOUNDARY = 1_500_000  # Arms: <=1.5Mb or >=7,167,508

# Motifs to scan (forward and reverse complement)
MOTIFS = {
    'TGGCCGGC': {
        'forward': 'TGGCCGGC',
        'revcomp': 'GCCGGCCA',
        'methylation_motif': 'GCCGGC',
        'label': 'GCCGGC (4mC recognition)'
    },
    'AAGCCCG': {
        'forward': 'AAGCCCG',
        'revcomp': 'CGGGCTT',
        'methylation_motif': 'AAGCCCG',
        'label': 'AAGCCCG (6mA recognition)'
    }
}

# Regions to analyze
PROMOTER_UPSTREAM = 500
PROMOTER_DOWNSTREAM = 100
EXTENDED_UPSTREAM = 2000

N_PERMUTATIONS = 10_000

print("=" * 70)
print("H22: Sequence-Level Motif Depletion Near Regulatory Genes")
print("=" * 70)

# ============================================================
# Step 1: Load genome sequence
# ============================================================
print("\n[Step 1] Loading genome sequence...")
genome_records = {rec.id: rec for rec in SeqIO.parse(str(GENOME_FASTA), 'fasta')}
genome_seq = str(genome_records['NC_003888.3'].seq).upper()
genome_len = len(genome_seq)
print(f"  Genome length: {genome_len:,} bp (NC_003888.3)")

# ============================================================
# Step 2: Scan genome for motif occurrences
# ============================================================
print("\n[Step 2] Scanning genome for motif occurrences...")

def find_motif_positions(sequence, motif_fwd, motif_rev):
    """Find all positions of a motif on both strands."""
    positions = []
    # Forward strand
    for match in re.finditer(motif_fwd, sequence):
        positions.append((match.start(), '+'))
    # Reverse strand (complement motif on forward sequence)
    for match in re.finditer(motif_rev, sequence):
        positions.append((match.start(), '-'))
    return sorted(positions, key=lambda x: x[0])

motif_positions = {}
for motif_name, motif_info in MOTIFS.items():
    positions = find_motif_positions(genome_seq, motif_info['forward'], motif_info['revcomp'])
    motif_positions[motif_name] = positions
    n_fwd = sum(1 for _, s in positions if s == '+')
    n_rev = sum(1 for _, s in positions if s == '-')
    density = len(positions) / (genome_len / 1000)
    print(f"  {motif_name}: {len(positions)} total ({n_fwd} fwd + {n_rev} rev), density = {density:.3f}/kb")

# ============================================================
# Step 3: Load gene annotations
# ============================================================
print("\n[Step 3] Loading gene annotations...")
genes_df = pd.read_csv(GENE_ANNOTATION, sep='\t')
# Filter to protein-coding genes on main chromosome
genes_df = genes_df[(genes_df['gene_biotype'] == 'protein_coding') &
                     (genes_df['contig'] == 'NC_003888.3')].copy()
print(f"  Total protein-coding genes on NC_003888.3: {len(genes_df)}")

# Load regulatory genes
reg_df = pd.read_csv(REGULATORY_GENES, sep='\t')
regulatory_loci = set(reg_df['locus_tag'].values)
print(f"  Regulatory genes: {len(regulatory_loci)}")

# Mark regulatory genes
genes_df['is_regulatory'] = genes_df['gene_id'].isin(regulatory_loci)
n_reg = genes_df['is_regulatory'].sum()
n_nonreg = (~genes_df['is_regulatory']).sum()
print(f"  Matched regulatory: {n_reg}, Non-regulatory: {n_nonreg}")

# Assign core/arm
def assign_region(start, end):
    midpoint = (start + end) / 2
    if midpoint <= ARM_BOUNDARY or midpoint >= (CHROMOSOME_LENGTH - ARM_BOUNDARY):
        return 'arm'
    return 'core'

genes_df['region'] = genes_df.apply(lambda r: assign_region(r['start'], r['end']), axis=1)
print(f"  Core genes: {(genes_df['region'] == 'core').sum()}, Arm genes: {(genes_df['region'] == 'arm').sum()}")

# ============================================================
# Step 4: Count motifs per gene in promoter, body, and extended regions
# ============================================================
print("\n[Step 4] Counting motifs per gene in each region...")

def count_motifs_in_region(positions_list, region_start, region_end):
    """Count motif occurrences within a genomic region."""
    # Clamp to valid range
    region_start = max(0, region_start)
    region_end = min(genome_len, region_end)
    count = 0
    for pos, strand in positions_list:
        if region_start <= pos < region_end:
            count += 1
    return count

# Pre-sort positions for binary search
from bisect import bisect_left, bisect_right

def count_motifs_fast(sorted_pos_array, region_start, region_end):
    """Fast motif counting using binary search on position arrays."""
    region_start = max(0, region_start)
    region_end = min(genome_len, region_end)
    left = bisect_left(sorted_pos_array, region_start)
    right = bisect_right(sorted_pos_array, region_end - 1)
    return right - left

# Create position arrays for fast lookup
motif_pos_arrays = {}
for motif_name, positions in motif_positions.items():
    pos_array = np.array([p[0] for p in positions])
    motif_pos_arrays[motif_name] = pos_array

results = []
for idx, gene in genes_df.iterrows():
    gene_start = gene['start']
    gene_end = gene['end']
    strand = gene['strand']
    gene_length = gene_end - gene_start

    for motif_name in MOTIFS:
        pos_array = motif_pos_arrays[motif_name]

        # Define regions based on strand
        if strand == '+':
            prom_start = gene_start - PROMOTER_UPSTREAM
            prom_end = gene_start + PROMOTER_DOWNSTREAM
            ext_start = gene_start - EXTENDED_UPSTREAM
            ext_end = gene_end
        else:
            prom_start = gene_end - PROMOTER_DOWNSTREAM
            prom_end = gene_end + PROMOTER_UPSTREAM
            ext_start = gene_start
            ext_end = gene_end + EXTENDED_UPSTREAM

        # Count motifs
        prom_count = count_motifs_fast(pos_array, prom_start, prom_end)
        body_count = count_motifs_fast(pos_array, gene_start, gene_end)
        ext_count = count_motifs_fast(pos_array, ext_start, ext_end)

        # Region lengths
        prom_len = min(prom_end, genome_len) - max(prom_start, 0)
        body_len = gene_length
        ext_len = min(ext_end, genome_len) - max(ext_start, 0)

        results.append({
            'gene_id': gene['gene_id'],
            'gene_name': gene.get('gene_name', ''),
            'start': gene_start,
            'end': gene_end,
            'strand': strand,
            'gene_length': gene_length,
            'is_regulatory': gene['is_regulatory'],
            'region': gene['region'],
            'product': gene.get('product', ''),
            'motif': motif_name,
            'promoter_count': prom_count,
            'promoter_length': prom_len,
            'body_count': body_count,
            'body_length': body_len,
            'extended_count': ext_count,
            'extended_length': ext_len,
            'promoter_density': prom_count / (prom_len / 1000) if prom_len > 0 else 0,
            'body_density': body_count / (body_len / 1000) if body_len > 0 else 0,
            'extended_density': ext_count / (ext_len / 1000) if ext_len > 0 else 0,
        })

motif_df = pd.DataFrame(results)
print(f"  Computed motif counts for {len(motif_df)} gene-motif pairs")

# ============================================================
# Step 5: Compare regulatory vs non-regulatory genes
# ============================================================
print("\n[Step 5] Comparing regulatory vs non-regulatory genes...")

stat_results = []

for motif_name in MOTIFS:
    mdf = motif_df[motif_df['motif'] == motif_name]
    reg = mdf[mdf['is_regulatory']]
    nonreg = mdf[~mdf['is_regulatory']]

    for zone_name, count_col, len_col, density_col in [
        ('promoter', 'promoter_count', 'promoter_length', 'promoter_density'),
        ('gene_body', 'body_count', 'body_length', 'body_density'),
        ('extended_2kb', 'extended_count', 'extended_length', 'extended_density'),
    ]:
        # Mean densities
        reg_density = reg[density_col].mean()
        nonreg_density = nonreg[density_col].mean()
        fold = reg_density / nonreg_density if nonreg_density > 0 else np.nan

        # Total counts and lengths for overall density
        reg_total_count = reg[count_col].sum()
        reg_total_len = reg[len_col].sum()
        nonreg_total_count = nonreg[count_col].sum()
        nonreg_total_len = nonreg[len_col].sum()

        reg_overall_density = reg_total_count / (reg_total_len / 1000) if reg_total_len > 0 else 0
        nonreg_overall_density = nonreg_total_count / (nonreg_total_len / 1000) if nonreg_total_len > 0 else 0
        overall_fold = reg_overall_density / nonreg_overall_density if nonreg_overall_density > 0 else np.nan

        # Wilcoxon rank-sum test on per-gene counts
        try:
            wilcox_stat, wilcox_p = stats.mannwhitneyu(
                reg[count_col].values, nonreg[count_col].values, alternative='two-sided'
            )
        except:
            wilcox_stat, wilcox_p = np.nan, np.nan

        # Fisher exact test: genes WITH motif vs WITHOUT motif
        reg_with = (reg[count_col] > 0).sum()
        reg_without = (reg[count_col] == 0).sum()
        nonreg_with = (nonreg[count_col] > 0).sum()
        nonreg_without = (nonreg[count_col] == 0).sum()

        fisher_table = [[reg_with, reg_without], [nonreg_with, nonreg_without]]
        try:
            fisher_or, fisher_p = stats.fisher_exact(fisher_table)
        except:
            fisher_or, fisher_p = np.nan, np.nan

        # Proportion with motif
        reg_prop = reg_with / len(reg) if len(reg) > 0 else 0
        nonreg_prop = nonreg_with / len(nonreg) if len(nonreg) > 0 else 0

        stat_results.append({
            'motif': motif_name,
            'zone': zone_name,
            'region': 'combined',
            'n_reg': len(reg),
            'n_nonreg': len(nonreg),
            'reg_mean_density': reg_density,
            'nonreg_mean_density': nonreg_density,
            'mean_density_fold': fold,
            'reg_overall_density': reg_overall_density,
            'nonreg_overall_density': nonreg_overall_density,
            'overall_density_fold': overall_fold,
            'reg_with_motif': reg_with,
            'reg_without_motif': reg_without,
            'nonreg_with_motif': nonreg_with,
            'nonreg_without_motif': nonreg_without,
            'reg_prop_with': reg_prop,
            'nonreg_prop_with': nonreg_prop,
            'fisher_OR': fisher_or,
            'fisher_p': fisher_p,
            'wilcoxon_U': wilcox_stat,
            'wilcoxon_p': wilcox_p,
        })

        print(f"  {motif_name} | {zone_name:15s} | reg_density={reg_density:.4f}/kb, nonreg={nonreg_density:.4f}/kb, "
              f"fold={fold:.3f}, Fisher_p={fisher_p:.2e}, Wilcox_p={wilcox_p:.2e}")

# ============================================================
# Step 6: Core/arm stratification
# ============================================================
print("\n[Step 6] Core/arm stratified analysis...")

for region_label in ['core', 'arm']:
    for motif_name in MOTIFS:
        mdf = motif_df[(motif_df['motif'] == motif_name) & (motif_df['region'] == region_label)]
        reg = mdf[mdf['is_regulatory']]
        nonreg = mdf[~mdf['is_regulatory']]

        for zone_name, count_col, len_col, density_col in [
            ('promoter', 'promoter_count', 'promoter_length', 'promoter_density'),
            ('gene_body', 'body_count', 'body_length', 'body_density'),
            ('extended_2kb', 'extended_count', 'extended_length', 'extended_density'),
        ]:
            reg_density = reg[density_col].mean()
            nonreg_density = nonreg[density_col].mean()
            fold = reg_density / nonreg_density if nonreg_density > 0 else np.nan

            reg_total_count = reg[count_col].sum()
            reg_total_len = reg[len_col].sum()
            nonreg_total_count = nonreg[count_col].sum()
            nonreg_total_len = nonreg[len_col].sum()

            reg_overall_density = reg_total_count / (reg_total_len / 1000) if reg_total_len > 0 else 0
            nonreg_overall_density = nonreg_total_count / (nonreg_total_len / 1000) if nonreg_total_len > 0 else 0
            overall_fold = reg_overall_density / nonreg_overall_density if nonreg_overall_density > 0 else np.nan

            try:
                wilcox_stat, wilcox_p = stats.mannwhitneyu(
                    reg[count_col].values, nonreg[count_col].values, alternative='two-sided'
                )
            except:
                wilcox_stat, wilcox_p = np.nan, np.nan

            reg_with = (reg[count_col] > 0).sum()
            reg_without = (reg[count_col] == 0).sum()
            nonreg_with = (nonreg[count_col] > 0).sum()
            nonreg_without = (nonreg[count_col] == 0).sum()

            try:
                fisher_or, fisher_p = stats.fisher_exact([[reg_with, reg_without], [nonreg_with, nonreg_without]])
            except:
                fisher_or, fisher_p = np.nan, np.nan

            reg_prop = reg_with / len(reg) if len(reg) > 0 else 0
            nonreg_prop = nonreg_with / len(nonreg) if len(nonreg) > 0 else 0

            stat_results.append({
                'motif': motif_name,
                'zone': zone_name,
                'region': region_label,
                'n_reg': len(reg),
                'n_nonreg': len(nonreg),
                'reg_mean_density': reg_density,
                'nonreg_mean_density': nonreg_density,
                'mean_density_fold': fold,
                'reg_overall_density': reg_overall_density,
                'nonreg_overall_density': nonreg_overall_density,
                'overall_density_fold': overall_fold,
                'reg_with_motif': reg_with,
                'reg_without_motif': reg_without,
                'nonreg_with_motif': nonreg_with,
                'nonreg_without_motif': nonreg_without,
                'reg_prop_with': reg_prop,
                'nonreg_prop_with': nonreg_prop,
                'fisher_OR': fisher_or,
                'fisher_p': fisher_p,
                'wilcoxon_U': wilcox_stat,
                'wilcoxon_p': wilcox_p,
            })

    print(f"  {region_label}: computed all motif x zone statistics")

stats_df = pd.DataFrame(stat_results)

# ============================================================
# Step 7: Compare DNA motif depletion vs methylation site depletion (H15/H20)
# ============================================================
print("\n[Step 7] Comparing DNA motif depletion vs methylation site depletion...")

# Load H15 data
h15_df = pd.read_csv(H15_ENRICHMENT, sep='\t')
# Load H20 data
h20_df = pd.read_csv(H20_STRATIFIED, sep='\t')

comparison_rows = []

# H15 methylation fold enrichment for Regulatory/TF (T1 data = highest power)
for motif_name, h15_motif_group, h20_motif_group in [
    ('TGGCCGGC', 'GCCGGC_4mC_T1', 'GCCGGC_4mC_T1'),
    ('AAGCCCG', 'AAGCCCG_6mA_T1', 'AAGCCCG_6mA_T1'),
]:
    # H15 unstratified methylation fold
    h15_row = h15_df[(h15_df['motif_group'] == h15_motif_group) & (h15_df['category'] == 'Regulatory/TF')]
    if len(h15_row) > 0:
        methyl_fold = h15_row['fold_enrichment'].values[0]
        methyl_p = h15_row['p_value'].values[0]
    else:
        methyl_fold, methyl_p = np.nan, np.nan

    # DNA motif fold (extended region, combined)
    dna_row = stats_df[(stats_df['motif'] == motif_name) &
                        (stats_df['zone'] == 'extended_2kb') &
                        (stats_df['region'] == 'combined')]
    if len(dna_row) > 0:
        dna_fold = dna_row['overall_density_fold'].values[0]
        dna_fisher_p = dna_row['fisher_p'].values[0]
        dna_wilcox_p = dna_row['wilcoxon_p'].values[0]
    else:
        dna_fold, dna_fisher_p, dna_wilcox_p = np.nan, np.nan, np.nan

    # Ratio: DNA fold / Methylation fold
    # If both are < 1 (depleted), ratio > 1 means DNA more depleted than methylation
    ratio = dna_fold / methyl_fold if methyl_fold > 0 else np.nan

    # Contribution: if DNA fold = methylation fold, 100% sequence-level
    # if DNA fold = 1.0 but methylation fold < 1, 100% protein occupancy
    if methyl_fold < 1:
        sequence_contribution_pct = ((1.0 - dna_fold) / (1.0 - methyl_fold)) * 100 if (1.0 - methyl_fold) != 0 else np.nan
    else:
        sequence_contribution_pct = np.nan

    comparison_rows.append({
        'motif': motif_name,
        'region': 'combined',
        'methylation_fold': methyl_fold,
        'methylation_p': methyl_p,
        'dna_motif_fold': dna_fold,
        'dna_fisher_p': dna_fisher_p,
        'dna_wilcoxon_p': dna_wilcox_p,
        'fold_ratio_dna_over_methyl': ratio,
        'sequence_contribution_pct': sequence_contribution_pct,
    })

    print(f"  {motif_name:10s} | Methyl fold={methyl_fold:.3f} (p={methyl_p:.2e}), "
          f"DNA fold={dna_fold:.3f} (Fisher p={dna_fisher_p:.2e}), "
          f"Ratio={ratio:.3f}, Seq contribution={sequence_contribution_pct:.1f}%")

    # Stratified (core/arm) from H20
    for region_label in ['core', 'arm']:
        h20_row = h20_df[(h20_df['motif_group'] == h20_motif_group) &
                          (h20_df['category'] == 'Regulatory/TF') &
                          (h20_df['region'] == region_label)]
        if len(h20_row) > 0:
            methyl_fold_strat = h20_row['fold_enrichment'].values[0]
            methyl_p_strat = h20_row['p_value'].values[0]
        else:
            methyl_fold_strat, methyl_p_strat = np.nan, np.nan

        dna_row_strat = stats_df[(stats_df['motif'] == motif_name) &
                                  (stats_df['zone'] == 'extended_2kb') &
                                  (stats_df['region'] == region_label)]
        if len(dna_row_strat) > 0:
            dna_fold_strat = dna_row_strat['overall_density_fold'].values[0]
            dna_fisher_p_strat = dna_row_strat['fisher_p'].values[0]
            dna_wilcox_p_strat = dna_row_strat['wilcoxon_p'].values[0]
        else:
            dna_fold_strat, dna_fisher_p_strat, dna_wilcox_p_strat = np.nan, np.nan, np.nan

        ratio_strat = dna_fold_strat / methyl_fold_strat if methyl_fold_strat > 0 else np.nan
        if methyl_fold_strat < 1:
            seq_contrib_strat = ((1.0 - dna_fold_strat) / (1.0 - methyl_fold_strat)) * 100 if (1.0 - methyl_fold_strat) != 0 else np.nan
        else:
            seq_contrib_strat = np.nan

        comparison_rows.append({
            'motif': motif_name,
            'region': region_label,
            'methylation_fold': methyl_fold_strat,
            'methylation_p': methyl_p_strat,
            'dna_motif_fold': dna_fold_strat,
            'dna_fisher_p': dna_fisher_p_strat,
            'dna_wilcoxon_p': dna_wilcox_p_strat,
            'fold_ratio_dna_over_methyl': ratio_strat,
            'sequence_contribution_pct': seq_contrib_strat,
        })

comparison_df = pd.DataFrame(comparison_rows)

# ============================================================
# Step 8: Permutation test
# ============================================================
print(f"\n[Step 8] Permutation test ({N_PERMUTATIONS:,} iterations)...")

np.random.seed(42)
n_regulatory = genes_df['is_regulatory'].sum()
all_gene_ids = genes_df['gene_id'].values

perm_results = []

for motif_name in MOTIFS:
    mdf = motif_df[motif_df['motif'] == motif_name].set_index('gene_id')

    # Observed: mean extended density for regulatory genes
    obs_density = mdf.loc[mdf['is_regulatory'], 'extended_density'].mean()
    obs_prop = (mdf.loc[mdf['is_regulatory'], 'extended_count'] > 0).mean()
    obs_prom_density = mdf.loc[mdf['is_regulatory'], 'promoter_density'].mean()

    # Null distribution
    null_densities = np.zeros(N_PERMUTATIONS)
    null_props = np.zeros(N_PERMUTATIONS)
    null_prom_densities = np.zeros(N_PERMUTATIONS)

    # Pre-extract arrays for speed
    all_ext_density = mdf['extended_density'].values
    all_ext_has = (mdf['extended_count'] > 0).astype(int).values
    all_prom_density = mdf['promoter_density'].values
    n_genes = len(mdf)

    for i in range(N_PERMUTATIONS):
        sample_idx = np.random.choice(n_genes, size=n_regulatory, replace=False)
        null_densities[i] = all_ext_density[sample_idx].mean()
        null_props[i] = all_ext_has[sample_idx].mean()
        null_prom_densities[i] = all_prom_density[sample_idx].mean()

    # P-value: proportion of null <= observed (one-sided, testing depletion)
    p_density = (null_densities <= obs_density).sum() / N_PERMUTATIONS
    p_prop = (null_props <= obs_prop).sum() / N_PERMUTATIONS
    p_prom = (null_prom_densities <= obs_prom_density).sum() / N_PERMUTATIONS

    null_mean = null_densities.mean()
    null_std = null_densities.std()
    z_score = (obs_density - null_mean) / null_std if null_std > 0 else np.nan

    perm_results.append({
        'motif': motif_name,
        'zone': 'extended_2kb',
        'observed_density': obs_density,
        'null_mean': null_mean,
        'null_std': null_std,
        'z_score': z_score,
        'p_value': p_density,
        'fold_vs_null': obs_density / null_mean if null_mean > 0 else np.nan,
    })
    perm_results.append({
        'motif': motif_name,
        'zone': 'promoter',
        'observed_density': obs_prom_density,
        'null_mean': null_prom_densities.mean(),
        'null_std': null_prom_densities.std(),
        'z_score': (obs_prom_density - null_prom_densities.mean()) / null_prom_densities.std() if null_prom_densities.std() > 0 else np.nan,
        'p_value': p_prom,
        'fold_vs_null': obs_prom_density / null_prom_densities.mean() if null_prom_densities.mean() > 0 else np.nan,
    })

    print(f"  {motif_name} extended: obs={obs_density:.4f}, null={null_mean:.4f}+/-{null_std:.4f}, "
          f"z={z_score:.2f}, p={p_density:.4f}")
    print(f"  {motif_name} promoter: obs={obs_prom_density:.4f}, null={null_prom_densities.mean():.4f}, "
          f"p={p_prom:.4f}")

perm_df = pd.DataFrame(perm_results)

# ============================================================
# Save tables
# ============================================================
print("\n[Saving tables]...")

# Gene motif counts (wide format) - merge instead of pivot to preserve all genes
motif_tg = motif_df[motif_df['motif'] == 'TGGCCGGC'][['gene_id', 'promoter_count', 'body_count', 'extended_count',
    'promoter_density', 'body_density', 'extended_density']].rename(
    columns={c: f'{c}_TGGCCGGC' for c in ['promoter_count', 'body_count', 'extended_count',
             'promoter_density', 'body_density', 'extended_density']})
motif_aa = motif_df[motif_df['motif'] == 'AAGCCCG'][['gene_id', 'promoter_count', 'body_count', 'extended_count',
    'promoter_density', 'body_density', 'extended_density']].rename(
    columns={c: f'{c}_AAGCCCG' for c in ['promoter_count', 'body_count', 'extended_count',
             'promoter_density', 'body_density', 'extended_density']})
gene_base = motif_df[motif_df['motif'] == 'TGGCCGGC'][['gene_id', 'gene_name', 'start', 'end', 'strand',
    'gene_length', 'is_regulatory', 'region', 'product']].copy()
gene_motif_wide = gene_base.merge(motif_tg, on='gene_id').merge(motif_aa, on='gene_id')
gene_motif_wide.to_csv(TABLES_DIR / 'gene_motif_counts.tsv', sep='\t', index=False)
print(f"  gene_motif_counts.tsv: {len(gene_motif_wide)} genes")

# Depletion statistics
stats_df.to_csv(TABLES_DIR / 'depletion_statistics.tsv', sep='\t', index=False)
print(f"  depletion_statistics.tsv: {len(stats_df)} rows")

# Comparison with methylation
comparison_df.to_csv(TABLES_DIR / 'sequence_vs_methylation_comparison.tsv', sep='\t', index=False)
print(f"  sequence_vs_methylation_comparison.tsv: {len(comparison_df)} rows")

# Stratified by region (subset of stats_df)
stratified = stats_df[stats_df['region'] != 'combined'].copy()
stratified.to_csv(TABLES_DIR / 'stratified_by_region.tsv', sep='\t', index=False)
print(f"  stratified_by_region.tsv: {len(stratified)} rows")

# Permutation results
perm_df.to_csv(TABLES_DIR / 'permutation_test_results.tsv', sep='\t', index=False)
print(f"  permutation_test_results.tsv: {len(perm_df)} rows")

# ============================================================
# Figure 1: Motif density comparison
# ============================================================
print("\n[Figure 1] Motif density comparison...")

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for ax_idx, motif_name in enumerate(MOTIFS):
    ax = axes[ax_idx]
    combined = stats_df[(stats_df['motif'] == motif_name) & (stats_df['region'] == 'combined')]

    zones = ['promoter', 'gene_body', 'extended_2kb']
    zone_labels = ['Promoter\n(-500 to +100)', 'Gene Body', 'Extended\n(-2kb to end)']

    reg_densities = [combined[combined['zone'] == z]['reg_mean_density'].values[0] for z in zones]
    nonreg_densities = [combined[combined['zone'] == z]['nonreg_mean_density'].values[0] for z in zones]
    folds = [combined[combined['zone'] == z]['mean_density_fold'].values[0] for z in zones]
    fisher_ps = [combined[combined['zone'] == z]['fisher_p'].values[0] for z in zones]
    wilcox_ps = [combined[combined['zone'] == z]['wilcoxon_p'].values[0] for z in zones]

    x = np.arange(len(zones))
    width = 0.35

    bars1 = ax.bar(x - width/2, reg_densities, width, label='Regulatory', color='#e74c3c', alpha=0.8)
    bars2 = ax.bar(x + width/2, nonreg_densities, width, label='Non-regulatory', color='#3498db', alpha=0.8)

    # Significance stars
    for i, (fp, wp) in enumerate(zip(fisher_ps, wilcox_ps)):
        min_p = min(fp, wp)
        max_y = max(reg_densities[i], nonreg_densities[i])
        if min_p < 0.001:
            ax.text(i, max_y * 1.05, '***', ha='center', fontsize=12, fontweight='bold')
        elif min_p < 0.01:
            ax.text(i, max_y * 1.05, '**', ha='center', fontsize=12, fontweight='bold')
        elif min_p < 0.05:
            ax.text(i, max_y * 1.05, '*', ha='center', fontsize=12, fontweight='bold')
        else:
            ax.text(i, max_y * 1.05, 'ns', ha='center', fontsize=10, color='gray')

    # Fold annotations
    for i, f in enumerate(folds):
        ax.text(i, -0.03 * ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 0,
                f'fold={f:.2f}', ha='center', fontsize=8, color='gray', va='top')

    ax.set_xlabel('')
    ax.set_ylabel('Motif density (per kb)')
    ax.set_title(f'{MOTIFS[motif_name]["label"]}')
    ax.set_xticks(x)
    ax.set_xticklabels(zone_labels, fontsize=9)
    ax.legend(fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

plt.suptitle('H22: DNA Motif Density Near Regulatory vs Non-Regulatory Genes', fontsize=13, fontweight='bold')
plt.tight_layout()
for fmt in ['pdf', 'svg', 'png']:
    plt.savefig(FIGURES_DIR / f'motif_density_comparison.{fmt}', dpi=200, bbox_inches='tight')
plt.close()
print("  Saved motif_density_comparison.pdf/svg/png")

# ============================================================
# Figure 2: Sequence vs methylation depletion
# ============================================================
print("\n[Figure 2] Sequence vs methylation depletion...")

fig, axes = plt.subplots(1, 2, figsize=(10, 5))

for ax_idx, motif_name in enumerate(MOTIFS):
    ax = axes[ax_idx]
    rows = comparison_df[comparison_df['motif'] == motif_name]

    regions = ['combined', 'core', 'arm']
    region_labels = ['Combined', 'Core only', 'Arm only']

    methyl_folds = []
    dna_folds = []
    for r in regions:
        row = rows[rows['region'] == r]
        if len(row) > 0:
            methyl_folds.append(row['methylation_fold'].values[0])
            dna_folds.append(row['dna_motif_fold'].values[0])
        else:
            methyl_folds.append(np.nan)
            dna_folds.append(np.nan)

    x = np.arange(len(regions))
    width = 0.35

    bars1 = ax.bar(x - width/2, dna_folds, width, label='DNA motif', color='#2ecc71', alpha=0.8, edgecolor='black', linewidth=0.5)
    bars2 = ax.bar(x + width/2, methyl_folds, width, label='Methylation site', color='#9b59b6', alpha=0.8, edgecolor='black', linewidth=0.5)

    ax.axhline(y=1.0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.set_ylabel('Fold enrichment near regulatory genes')
    ax.set_title(f'{MOTIFS[motif_name]["label"]}')
    ax.set_xticks(x)
    ax.set_xticklabels(region_labels)
    ax.legend(fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylim(0, 1.5)

    # Add fold values
    for i, (df, mf) in enumerate(zip(dna_folds, methyl_folds)):
        if not np.isnan(df):
            ax.text(x[i] - width/2, df + 0.02, f'{df:.2f}', ha='center', fontsize=8, fontweight='bold')
        if not np.isnan(mf):
            ax.text(x[i] + width/2, mf + 0.02, f'{mf:.2f}', ha='center', fontsize=8, fontweight='bold')

plt.suptitle('H22: DNA Sequence vs Methylation Depletion Near Regulatory Genes', fontsize=13, fontweight='bold')
plt.tight_layout()
for fmt in ['pdf', 'svg', 'png']:
    plt.savefig(FIGURES_DIR / f'sequence_vs_methylation_depletion.{fmt}', dpi=200, bbox_inches='tight')
plt.close()
print("  Saved sequence_vs_methylation_depletion.pdf/svg/png")

# ============================================================
# Figure 3: H22 Comprehensive Summary (multi-panel)
# ============================================================
print("\n[Figure 3] Comprehensive summary...")

fig = plt.figure(figsize=(16, 14))
gs = gridspec.GridSpec(3, 2, hspace=0.4, wspace=0.3)

# Panel A: Genome-wide motif counts
ax_a = fig.add_subplot(gs[0, 0])
motif_labels = []
motif_counts = []
for mn in MOTIFS:
    motif_labels.append(mn)
    motif_counts.append(len(motif_positions[mn]))
ax_a.barh(motif_labels, motif_counts, color=['#e74c3c', '#3498db'], alpha=0.8)
ax_a.set_xlabel('Total motif occurrences (both strands)')
ax_a.set_title('A. Genome-wide Motif Counts', fontweight='bold')
for i, v in enumerate(motif_counts):
    density = v / (genome_len / 1000)
    ax_a.text(v + 10, i, f'{v:,} ({density:.2f}/kb)', va='center', fontsize=9)
ax_a.spines['top'].set_visible(False)
ax_a.spines['right'].set_visible(False)

# Panel B: Fold enrichment by zone (extended region, both motifs)
ax_b = fig.add_subplot(gs[0, 1])
for ax_idx, motif_name in enumerate(MOTIFS):
    combined = stats_df[(stats_df['motif'] == motif_name) & (stats_df['region'] == 'combined')]
    zones = ['promoter', 'gene_body', 'extended_2kb']
    folds = [combined[combined['zone'] == z]['overall_density_fold'].values[0] for z in zones]
    x_pos = np.arange(len(zones))
    offset = -0.2 + ax_idx * 0.4
    color = '#e74c3c' if ax_idx == 0 else '#3498db'
    bars = ax_b.bar(x_pos + offset, folds, 0.35, label=motif_name, color=color, alpha=0.8)
    for i, f in enumerate(folds):
        ax_b.text(x_pos[i] + offset, f + 0.01, f'{f:.2f}', ha='center', fontsize=8)

ax_b.axhline(y=1.0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
ax_b.set_xticks(np.arange(3))
ax_b.set_xticklabels(['Promoter', 'Gene body', 'Extended 2kb'])
ax_b.set_ylabel('Fold (regulatory / non-regulatory)')
ax_b.set_title('B. DNA Motif Density Fold by Gene Zone', fontweight='bold')
ax_b.legend(fontsize=9)
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)

# Panel C: Core/arm stratification
ax_c = fig.add_subplot(gs[1, 0])
strat_data = []
for motif_name in MOTIFS:
    for region in ['core', 'arm']:
        row = stats_df[(stats_df['motif'] == motif_name) &
                        (stats_df['zone'] == 'extended_2kb') &
                        (stats_df['region'] == region)]
        if len(row) > 0:
            strat_data.append({
                'motif': motif_name,
                'region': region,
                'fold': row['overall_density_fold'].values[0],
                'fisher_p': row['fisher_p'].values[0],
            })

strat_plot = pd.DataFrame(strat_data)
x_labels = []
folds_plot = []
colors_plot = []
for _, r in strat_plot.iterrows():
    x_labels.append(f'{r["motif"]}\n{r["region"]}')
    folds_plot.append(r['fold'])
    colors_plot.append('#e74c3c' if 'TGGCCGGC' in r['motif'] else '#3498db')

bars = ax_c.bar(range(len(x_labels)), folds_plot, color=colors_plot, alpha=0.8, edgecolor='black', linewidth=0.5)
ax_c.axhline(y=1.0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
ax_c.set_xticks(range(len(x_labels)))
ax_c.set_xticklabels(x_labels, fontsize=8)
ax_c.set_ylabel('Fold (reg / non-reg)')
ax_c.set_title('C. Geographic Stratification (Extended 2kb)', fontweight='bold')
for i, (f, r) in enumerate(zip(folds_plot, strat_plot.itertuples())):
    sig = '***' if r.fisher_p < 0.001 else '**' if r.fisher_p < 0.01 else '*' if r.fisher_p < 0.05 else 'ns'
    ax_c.text(i, f + 0.01, f'{f:.2f}\n{sig}', ha='center', fontsize=8)
ax_c.spines['top'].set_visible(False)
ax_c.spines['right'].set_visible(False)

# Panel D: DNA vs Methylation comparison
ax_d = fig.add_subplot(gs[1, 1])
comp_combined = comparison_df[comparison_df['region'] == 'combined']
x_motifs = comp_combined['motif'].values
x_pos = np.arange(len(x_motifs))
dna_bars = comp_combined['dna_motif_fold'].values
met_bars = comp_combined['methylation_fold'].values
width = 0.35
bars1 = ax_d.bar(x_pos - width/2, dna_bars, width, label='DNA motif fold', color='#2ecc71', alpha=0.8)
bars2 = ax_d.bar(x_pos + width/2, met_bars, width, label='Methylation fold (H15)', color='#9b59b6', alpha=0.8)
ax_d.axhline(y=1.0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
ax_d.set_xticks(x_pos)
ax_d.set_xticklabels(x_motifs)
ax_d.set_ylabel('Fold enrichment near regulatory genes')
ax_d.set_title('D. Sequence vs Methylation Depletion', fontweight='bold')
ax_d.legend(fontsize=9)
for i, (d, m) in enumerate(zip(dna_bars, met_bars)):
    ax_d.text(x_pos[i] - width/2, d + 0.01, f'{d:.2f}', ha='center', fontsize=9, fontweight='bold')
    ax_d.text(x_pos[i] + width/2, m + 0.01, f'{m:.2f}', ha='center', fontsize=9, fontweight='bold')
ax_d.spines['top'].set_visible(False)
ax_d.spines['right'].set_visible(False)
ax_d.set_ylim(0, 1.5)

# Panel E: Permutation test results
ax_e = fig.add_subplot(gs[2, 0])
# Plot null distribution for TGGCCGGC extended
mdf_perm = motif_df[motif_df['motif'] == 'TGGCCGGC'].set_index('gene_id')
all_ext_density_tg = mdf_perm['extended_density'].values
obs_tg = mdf_perm.loc[mdf_perm['is_regulatory'], 'extended_density'].mean()
null_tg = np.array([all_ext_density_tg[np.random.choice(len(all_ext_density_tg), size=n_regulatory, replace=False)].mean()
                     for _ in range(N_PERMUTATIONS)])

ax_e.hist(null_tg, bins=50, color='gray', alpha=0.7, density=True, label='Null distribution')
ax_e.axvline(x=obs_tg, color='red', linewidth=2, linestyle='--', label=f'Observed: {obs_tg:.4f}')
ax_e.set_xlabel('Mean motif density (per kb)')
ax_e.set_ylabel('Density')
ax_e.set_title('E. Permutation Test: TGGCCGGC Extended Region', fontweight='bold')
ax_e.legend(fontsize=9)
ax_e.spines['top'].set_visible(False)
ax_e.spines['right'].set_visible(False)

# Panel F: Permutation test for AAGCCCG
ax_f = fig.add_subplot(gs[2, 1])
mdf_perm2 = motif_df[motif_df['motif'] == 'AAGCCCG'].set_index('gene_id')
all_ext_density_aa = mdf_perm2['extended_density'].values
obs_aa = mdf_perm2.loc[mdf_perm2['is_regulatory'], 'extended_density'].mean()
null_aa = np.array([all_ext_density_aa[np.random.choice(len(all_ext_density_aa), size=n_regulatory, replace=False)].mean()
                     for _ in range(N_PERMUTATIONS)])

ax_f.hist(null_aa, bins=50, color='gray', alpha=0.7, density=True, label='Null distribution')
ax_f.axvline(x=obs_aa, color='red', linewidth=2, linestyle='--', label=f'Observed: {obs_aa:.4f}')
ax_f.set_xlabel('Mean motif density (per kb)')
ax_f.set_ylabel('Density')
ax_f.set_title('F. Permutation Test: AAGCCCG Extended Region', fontweight='bold')
ax_f.legend(fontsize=9)
ax_f.spines['top'].set_visible(False)
ax_f.spines['right'].set_visible(False)

plt.suptitle('H22: Sequence-Level Motif Depletion Near Regulatory Genes\nComprehensive Summary',
             fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()
for fmt in ['pdf', 'svg', 'png']:
    plt.savefig(FIGURES_DIR / f'H22_comprehensive_summary.{fmt}', dpi=200, bbox_inches='tight')
plt.close()
print("  Saved H22_comprehensive_summary.pdf/svg/png")

# ============================================================
# Summary
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print("\n--- DNA Motif Density Fold (regulatory / non-regulatory) ---")
for motif_name in MOTIFS:
    for zone in ['promoter', 'gene_body', 'extended_2kb']:
        row = stats_df[(stats_df['motif'] == motif_name) &
                        (stats_df['zone'] == zone) &
                        (stats_df['region'] == 'combined')]
        if len(row) > 0:
            r = row.iloc[0]
            sig = 'YES' if min(r['fisher_p'], r['wilcoxon_p']) < 0.05 else 'NO'
            print(f"  {motif_name:10s} {zone:15s}: fold={r['mean_density_fold']:.3f}, "
                  f"overall_fold={r['overall_density_fold']:.3f}, "
                  f"Fisher_p={r['fisher_p']:.2e}, Wilcox_p={r['wilcoxon_p']:.2e}, Significant={sig}")

print("\n--- DNA vs Methylation Comparison ---")
for _, r in comparison_df[comparison_df['region'] == 'combined'].iterrows():
    print(f"  {r['motif']:10s}: DNA_fold={r['dna_motif_fold']:.3f}, Methyl_fold={r['methylation_fold']:.3f}, "
          f"Seq_contribution={r['sequence_contribution_pct']:.1f}%")

print("\n--- Permutation Test ---")
for _, r in perm_df.iterrows():
    print(f"  {r['motif']:10s} {r['zone']:15s}: obs={r['observed_density']:.4f}, "
          f"null={r['null_mean']:.4f}, z={r['z_score']:.2f}, p={r['p_value']:.4f}")

print("\n--- Interpretation ---")
for motif_name in MOTIFS:
    comp_row = comparison_df[(comparison_df['motif'] == motif_name) & (comparison_df['region'] == 'combined')]
    if len(comp_row) > 0:
        r = comp_row.iloc[0]
        dna_f = r['dna_motif_fold']
        met_f = r['methylation_fold']
        seq_pct = r['sequence_contribution_pct']

        if dna_f < 0.85 and met_f < 0.85:
            if abs(dna_f - met_f) < 0.15:
                model = "Primarily SEQUENCE-LEVEL counter-selection"
            elif dna_f > met_f:
                model = "BOTH mechanisms: sequence-level + protein occupancy"
            else:
                model = "Primarily SEQUENCE-LEVEL counter-selection (DNA more depleted)"
        elif dna_f >= 0.85 and met_f < 0.85:
            model = "Primarily PROTEIN OCCUPANCY model"
        else:
            model = "No significant depletion"

        print(f"  {motif_name}: {model}")
        print(f"    DNA fold={dna_f:.3f}, Methylation fold={met_f:.3f}, Sequence contribution={seq_pct:.1f}%")

print("\nDone!")
