#!/usr/bin/env python3
"""
H36: Evolutionary Conservation Analysis of 62 Exposed Transcription Factors
============================================================================

Compares evolutionary conservation proxies between 62 "exposed" regulators
(lacking methylation protection zones) and 993 "shielded" regulators in
Streptomyces coelicolor A3(2) M145.

Approaches:
- Sequence composition metrics (GC3, Nc, gene length)
- Annotation quality (gene_name, old_locus_tag, hypothetical status)
- Chromosomal position (core vs arm, distance to oriC)
- Expression level as conservation proxy
- Integrated composite conservation score
"""

import os
import sys
import warnings
import math
from collections import Counter, defaultdict
from itertools import product as itertools_product

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import mannwhitneyu, fisher_exact, wilcoxon, kruskal, spearmanr
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import StandardScaler
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.ticker as mticker

warnings.filterwarnings('ignore')

# ============================================================
# Configuration
# ============================================================
BASE_DIR = "/Users/okaban/bioinfo/rna-seq"
ANALYSIS_DIR = f"{BASE_DIR}/11_epigenome_integration/analysis/59_exposed_TF_conservation"
FIGURES_DIR = f"{ANALYSIS_DIR}/figures"
TABLES_DIR = f"{ANALYSIS_DIR}/tables"

GENE_ANNOTATION = f"{BASE_DIR}/05_annotation/analysis/05_annotation_260128_v1/tables/gene_annotation_basic.tsv"
EXPOSED_TABLE = f"{BASE_DIR}/11_epigenome_integration/analysis/51_exposed_regulators_characteristics/tables/exposed_regulators_full_table.tsv"
ALL_GENES_TABLE = f"{BASE_DIR}/11_epigenome_integration/analysis/52_shielded_exposed_boundary/tables/all_genes_features.tsv"
GENOME_FASTA = f"{BASE_DIR}/11_epigenome_integration/data/NC_003888.3.fna"

# Chromosome parameters
CHROM_LENGTH = 8_667_507
ORIC_POSITION = 4_272_693  # oriC near center of linear chromosome
ARM_BOUNDARY = 1_500_000   # core/arm boundary at 1.5 Mb from each end
LEFT_ARM_END = ARM_BOUNDARY
RIGHT_ARM_START = CHROM_LENGTH - ARM_BOUNDARY

# Standard codon table
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# Group codons by amino acid
AA_TO_CODONS = defaultdict(list)
for codon, aa in CODON_TABLE.items():
    if aa != '*':
        AA_TO_CODONS[aa].append(codon)

# Amino acids with synonymous codons (exclude Met, Trp)
SYNONYMOUS_AAS = {aa: codons for aa, codons in AA_TO_CODONS.items() if len(codons) > 1}

# Rare codon definition for S. coelicolor (GC-rich organism)
# AT-rich codons are rare in high-GC organisms
RARE_CODONS_HIGH_GC = {
    'TTA', 'TCA', 'ACA', 'ATA', 'AGA', 'GGA',  # A-ending
    'TTT', 'TCT', 'ACT', 'ATT', 'AGT', 'GGT',  # T-ending
    'AAA', 'AAT', 'TAT', 'GAT',  # AA-rich
}

# ============================================================
# Helper functions
# ============================================================
def load_genome(fasta_path):
    """Load genome sequence from FASTA file."""
    print(f"Loading genome from {fasta_path}...")
    seq = []
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith('>'):
                seq.append(line.strip().upper())
    genome = ''.join(seq)
    print(f"  Genome length: {len(genome):,} bp, GC: {(genome.count('G')+genome.count('C'))/len(genome)*100:.1f}%")
    return genome


def reverse_complement(seq):
    """Return reverse complement of DNA sequence."""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq))


def extract_gene_sequence(genome, start, end, strand):
    """Extract gene sequence (1-based coordinates)."""
    seq = genome[start-1:end]
    if strand == '-':
        seq = reverse_complement(seq)
    return seq


def calculate_gc_content(seq):
    """Calculate GC content of a sequence."""
    if len(seq) == 0:
        return np.nan
    gc = seq.count('G') + seq.count('C')
    return gc / len(seq)


def calculate_gc3(cds):
    """Calculate GC content at third codon position."""
    if len(cds) < 3:
        return np.nan
    third_positions = []
    for i in range(0, len(cds) - 2, 3):
        codon = cds[i:i+3]
        if len(codon) == 3 and 'N' not in codon:
            third_positions.append(codon[2])
    if not third_positions:
        return np.nan
    gc3 = sum(1 for b in third_positions if b in ('G', 'C')) / len(third_positions)
    return gc3


def calculate_nc(cds):
    """
    Calculate effective number of codons (Nc) using Wright (1990) formula.
    Nc = 2 + 9/F2 + 1/F3 + 5/F4 + 3/F6
    where Fn is the average homozygosity for amino acids with n-fold degeneracy.
    """
    if len(cds) < 3:
        return np.nan

    # Count codons
    codon_counts = Counter()
    for i in range(0, len(cds) - 2, 3):
        codon = cds[i:i+3]
        if len(codon) == 3 and 'N' not in codon and CODON_TABLE.get(codon) != '*':
            codon_counts[codon] += 1

    if not codon_counts:
        return np.nan

    # Group by degeneracy class
    degeneracy_classes = {2: [], 3: [], 4: [], 6: []}

    for aa, codons in SYNONYMOUS_AAS.items():
        n_fold = len(codons)
        # Get counts for this amino acid's codons
        counts = [codon_counts.get(c, 0) for c in codons]
        total = sum(counts)
        if total <= 1:
            continue

        # Calculate F (homozygosity) = sum(pi^2) adjusted
        # F = (n * sum(pi^2) - 1) / (n - 1)
        pi_sq = sum((c/total)**2 for c in counts)
        F = (total * pi_sq - 1) / (total - 1)

        if n_fold in degeneracy_classes:
            degeneracy_classes[n_fold].append(F)
        elif n_fold == 6:
            degeneracy_classes[6].append(F)

    # Handle Ile (3-fold) separately if in 3-fold class
    # Ser and Arg are 6-fold, Leu is 6-fold
    # Standard: 2-fold (9 aa), 3-fold (1 aa: Ile), 4-fold (5 aa), 6-fold (3 aa: Ser, Leu, Arg)

    # Average F for each degeneracy class
    Nc = 2  # Met + Trp (no synonyms)

    for n_fold, F_values in degeneracy_classes.items():
        if F_values:
            avg_F = np.mean(F_values)
            if avg_F > 0:
                if n_fold == 2:
                    Nc += 9 / avg_F
                elif n_fold == 3:
                    Nc += 1 / avg_F
                elif n_fold == 4:
                    Nc += 5 / avg_F
                elif n_fold == 6:
                    Nc += 3 / avg_F
            else:
                # F=0 means maximum bias
                if n_fold == 2:
                    Nc += 9
                elif n_fold == 3:
                    Nc += 1
                elif n_fold == 4:
                    Nc += 5
                elif n_fold == 6:
                    Nc += 3
        else:
            # Use maximum (no bias) if no data
            if n_fold == 2:
                Nc += 9
            elif n_fold == 3:
                Nc += 1
            elif n_fold == 4:
                Nc += 5
            elif n_fold == 6:
                Nc += 3

    # Nc ranges from 20 (extreme bias) to 61 (no bias)
    return min(max(Nc, 20), 61)


def calculate_rare_codon_freq(cds):
    """Calculate fraction of rare codons (AT-rich codons rare in high-GC organism)."""
    if len(cds) < 3:
        return np.nan

    total_codons = 0
    rare_count = 0
    for i in range(0, len(cds) - 2, 3):
        codon = cds[i:i+3]
        if len(codon) == 3 and 'N' not in codon and CODON_TABLE.get(codon) != '*':
            total_codons += 1
            if codon in RARE_CODONS_HIGH_GC:
                rare_count += 1

    if total_codons == 0:
        return np.nan
    return rare_count / total_codons


def classify_region(start, end):
    """Classify gene position as arm or core."""
    midpoint = (start + end) / 2
    if midpoint <= LEFT_ARM_END or midpoint >= RIGHT_ARM_START:
        return 'arm'
    else:
        return 'core'


def distance_to_oriC(start, end):
    """Calculate distance from gene midpoint to oriC."""
    midpoint = (start + end) / 2
    return abs(midpoint - ORIC_POSITION)


def effect_size_r(U, n1, n2):
    """Calculate effect size r from Mann-Whitney U statistic."""
    z = (U - n1 * n2 / 2) / math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
    return abs(z) / math.sqrt(n1 + n2)


def cohens_d(group1, group2):
    """Calculate Cohen's d effect size."""
    n1, n2 = len(group1), len(group2)
    var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
    pooled_std = math.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
    if pooled_std == 0:
        return 0
    return (np.mean(group1) - np.mean(group2)) / pooled_std


# ============================================================
# MAIN ANALYSIS
# ============================================================
print("=" * 70)
print("H36: Evolutionary Conservation of Exposed TFs")
print("=" * 70)

# ----------------------------------------------------------
# Step 1: Load data
# ----------------------------------------------------------
print("\n--- Step 1: Loading data ---")

# Gene annotation
gene_annot = pd.read_csv(GENE_ANNOTATION, sep='\t')
print(f"Gene annotation: {len(gene_annot)} genes")

# All regulatory genes
all_reg = pd.read_csv(ALL_GENES_TABLE, sep='\t')
print(f"All regulatory genes: {len(all_reg)} genes")
print(f"  Exposed: {all_reg['is_exposed'].sum()}")
print(f"  Shielded: {(~all_reg['is_exposed'].astype(bool)).sum()}")

# Exposed regulators
exposed = pd.read_csv(EXPOSED_TABLE, sep='\t')
exposed_loci = set(exposed['locus_tag'])
print(f"Exposed regulators (detailed table): {len(exposed)}")

# Merge with gene annotation to get protein_id
all_reg = all_reg.merge(
    gene_annot[['gene_id', 'protein_id', 'gene_biotype']].rename(columns={'gene_id': 'locus_tag'}),
    on='locus_tag', how='left'
)

# Also ensure old_locus_tag and gene_name from gene_annot are used where needed
annot_lookup = gene_annot.set_index('gene_id')[['old_locus_tag', 'gene_name', 'product', 'protein_id']].to_dict('index')

# Create group labels
all_reg['group'] = all_reg['is_exposed'].apply(lambda x: 'Exposed' if x else 'Shielded')
all_reg['is_exposed_bool'] = all_reg['is_exposed'].astype(bool)

print(f"\nExposed: {all_reg['group'].value_counts().get('Exposed', 0)}")
print(f"Shielded: {all_reg['group'].value_counts().get('Shielded', 0)}")

# ----------------------------------------------------------
# Step 2: Survey available Streptomyces genomes
# ----------------------------------------------------------
print("\n--- Step 2: Surveying available Streptomyces genomes ---")

genomes_dir = f"{BASE_DIR}/11_epigenome_integration/data/ncbi_genomes_all_species/extracted/ncbi_dataset/data"
gcf_dirs = [d for d in os.listdir(genomes_dir) if d.startswith('GCF_')]
print(f"Total genome assemblies available: {len(gcf_dirs)}")

# Check assembly report for species names
import json
species_info = {}
report_file = os.path.join(genomes_dir, 'assembly_data_report.jsonl')
if os.path.exists(report_file):
    with open(report_file) as f:
        for line in f:
            try:
                d = json.loads(line)
                acc = d.get('accession', '')
                org = d.get('organism', {}).get('organismName', 'Unknown')
                species_info[acc] = org
            except:
                pass
    print(f"Species info loaded for {len(species_info)} assemblies")
    # Show some species
    sample = list(species_info.items())[:10]
    for acc, sp in sample:
        print(f"  {acc}: {sp}")

# ----------------------------------------------------------
# Step 3: Gene-level conservation using annotation proxies
# ----------------------------------------------------------
print("\n--- Step 3: Annotation-based conservation proxies ---")

# Approach C: Gene name and product conservation
# has_gene_name: genes with assigned gene names are better characterized
# has_old_locus_tag: present in original S. coelicolor annotation (SCO numbers)
# is_hypothetical: hypothetical proteins are often lineage-specific

all_reg['has_gene_name'] = all_reg['gene_name'].notna() & (all_reg['gene_name'] != '')
all_reg['has_old_locus_tag'] = all_reg['old_locus_tag'].notna() & (all_reg['old_locus_tag'] != '')
all_reg['is_hypothetical'] = all_reg['product'].str.contains('hypothetical', case=False, na=False)
all_reg['has_named_product'] = ~all_reg['is_hypothetical']

# Also check protein_id presence
all_reg['has_protein_id'] = all_reg['protein_id'].notna() & (all_reg['protein_id'] != '')

for col in ['has_gene_name', 'has_old_locus_tag', 'has_named_product', 'has_protein_id']:
    exp_frac = all_reg.loc[all_reg['group']=='Exposed', col].mean()
    shi_frac = all_reg.loc[all_reg['group']=='Shielded', col].mean()
    print(f"  {col}: Exposed={exp_frac:.3f}, Shielded={shi_frac:.3f}")

# ----------------------------------------------------------
# Step 4: Sequence feature analysis
# ----------------------------------------------------------
print("\n--- Step 4: Extracting and analyzing gene sequences ---")

genome = load_genome(GENOME_FASTA)
genome_gc = calculate_gc_content(genome)
print(f"Genome GC content: {genome_gc:.4f}")

# Calculate genome-wide GC3 for reference
# Use all coding genes from annotation
all_cds_gc3 = []
for _, row in gene_annot.iterrows():
    if row.get('gene_biotype') == 'protein_coding':
        try:
            seq = extract_gene_sequence(genome, int(row['start']), int(row['end']), row['strand'])
            gc3 = calculate_gc3(seq)
            if not np.isnan(gc3):
                all_cds_gc3.append(gc3)
        except:
            pass

genome_avg_gc3 = np.mean(all_cds_gc3) if all_cds_gc3 else 0.9  # high GC3 expected
print(f"Genome average GC3 (from {len(all_cds_gc3)} coding genes): {genome_avg_gc3:.4f}")

# Calculate metrics for each regulatory gene
metrics = []
for _, row in all_reg.iterrows():
    locus = row['locus_tag']
    start = int(row['start'])
    end = int(row['end'])
    strand = row['strand']

    # Extract sequence
    try:
        cds = extract_gene_sequence(genome, start, end, strand)
    except:
        cds = ''

    # GC content (gene_length already in all_reg)
    gc_content = calculate_gc_content(cds) if cds else np.nan

    # GC3
    gc3 = calculate_gc3(cds) if cds else np.nan

    # GC3 deviation from genome average
    gc3_deviation = abs(gc3 - genome_avg_gc3) if not np.isnan(gc3) else np.nan

    # Effective number of codons
    nc = calculate_nc(cds) if cds else np.nan

    # Rare codon frequency
    rare_freq = calculate_rare_codon_freq(cds) if cds else np.nan

    # Distance to oriC
    dist_oriC = distance_to_oriC(start, end)

    # Region classification
    region_calc = classify_region(start, end)

    metrics.append({
        'locus_tag': locus,
        'gc_content': gc_content,
        'gc3': gc3,
        'gc3_deviation': gc3_deviation,
        'nc': nc,
        'rare_codon_freq': rare_freq,
        'dist_oriC': dist_oriC,
        'region_calc': region_calc,
    })

metrics_df = pd.DataFrame(metrics)
all_reg = all_reg.merge(metrics_df, on='locus_tag', how='left')

print(f"\nSequence metrics calculated for {len(metrics_df)} genes")
print(f"  Median gene length: Exposed={all_reg.loc[all_reg['group']=='Exposed', 'gene_length'].median():.0f}, "
      f"Shielded={all_reg.loc[all_reg['group']=='Shielded', 'gene_length'].median():.0f}")
print(f"  Median GC3: Exposed={all_reg.loc[all_reg['group']=='Exposed', 'gc3'].median():.4f}, "
      f"Shielded={all_reg.loc[all_reg['group']=='Shielded', 'gc3'].median():.4f}")
print(f"  Median Nc: Exposed={all_reg.loc[all_reg['group']=='Exposed', 'nc'].median():.1f}, "
      f"Shielded={all_reg.loc[all_reg['group']=='Shielded', 'nc'].median():.1f}")

# ----------------------------------------------------------
# Step 5: Functional constraint indicators (Fisher exact tests)
# ----------------------------------------------------------
print("\n--- Step 5: Functional constraint indicators ---")

annotation_results = []

for col, label in [('has_gene_name', 'Has gene name'),
                    ('has_old_locus_tag', 'Has SCO locus tag'),
                    ('has_named_product', 'Named product (non-hypothetical)'),
                    ('has_protein_id', 'Has protein_id (WP_)')]:

    exp = all_reg[all_reg['group'] == 'Exposed']
    shi = all_reg[all_reg['group'] == 'Shielded']

    a = exp[col].sum()  # exposed + feature
    b = len(exp) - a    # exposed - feature
    c = shi[col].sum()  # shielded + feature
    d = len(shi) - c    # shielded - feature

    table = [[a, b], [c, d]]
    odds_ratio, p_val = fisher_exact(table)

    exp_pct = a / len(exp) * 100
    shi_pct = c / len(shi) * 100

    print(f"  {label}: Exposed={exp_pct:.1f}%, Shielded={shi_pct:.1f}%, OR={odds_ratio:.3f}, p={p_val:.4e}")

    annotation_results.append({
        'feature': label,
        'exposed_count': int(a),
        'exposed_total': len(exp),
        'exposed_pct': exp_pct,
        'shielded_count': int(c),
        'shielded_total': len(shi),
        'shielded_pct': shi_pct,
        'odds_ratio': odds_ratio,
        'p_value': p_val,
        'test': 'Fisher exact'
    })

annot_df = pd.DataFrame(annotation_results)
annot_df.to_csv(f"{TABLES_DIR}/annotation_quality.tsv", sep='\t', index=False)
print(f"\nSaved: {TABLES_DIR}/annotation_quality.tsv")

# ----------------------------------------------------------
# Step 6: Chromosomal position as conservation indicator
# ----------------------------------------------------------
print("\n--- Step 6: Chromosomal position analysis ---")

# Core vs arm distribution
exp_core = all_reg.loc[all_reg['group']=='Exposed', 'region_calc'].value_counts().get('core', 0)
exp_arm = all_reg.loc[all_reg['group']=='Exposed', 'region_calc'].value_counts().get('arm', 0)
shi_core = all_reg.loc[all_reg['group']=='Shielded', 'region_calc'].value_counts().get('core', 0)
shi_arm = all_reg.loc[all_reg['group']=='Shielded', 'region_calc'].value_counts().get('arm', 0)

table_region = [[exp_core, exp_arm], [shi_core, shi_arm]]
or_region, p_region = fisher_exact(table_region)

print(f"  Core/Arm distribution:")
print(f"    Exposed: core={exp_core} ({exp_core/(exp_core+exp_arm)*100:.1f}%), arm={exp_arm} ({exp_arm/(exp_core+exp_arm)*100:.1f}%)")
print(f"    Shielded: core={shi_core} ({shi_core/(shi_core+shi_arm)*100:.1f}%), arm={shi_arm} ({shi_arm/(shi_core+shi_arm)*100:.1f}%)")
print(f"    Fisher OR={or_region:.3f}, p={p_region:.4e}")

# Distance to oriC comparison
exp_dist = all_reg.loc[all_reg['group']=='Exposed', 'dist_oriC'].dropna()
shi_dist = all_reg.loc[all_reg['group']=='Shielded', 'dist_oriC'].dropna()
U_dist, p_dist = mannwhitneyu(exp_dist, shi_dist, alternative='two-sided')
r_dist = effect_size_r(U_dist, len(exp_dist), len(shi_dist))

print(f"\n  Distance to oriC:")
print(f"    Exposed median: {exp_dist.median()/1e6:.2f} Mb")
print(f"    Shielded median: {shi_dist.median()/1e6:.2f} Mb")
print(f"    Mann-Whitney p={p_dist:.4e}, r={r_dist:.3f}")

# ----------------------------------------------------------
# Step 7: Expression level as conservation proxy
# ----------------------------------------------------------
print("\n--- Step 7: Expression level analysis ---")

exp_basemean = all_reg.loc[all_reg['group']=='Exposed', 'baseMean'].dropna()
shi_basemean = all_reg.loc[all_reg['group']=='Shielded', 'baseMean'].dropna()

U_bm, p_bm = mannwhitneyu(exp_basemean, shi_basemean, alternative='two-sided')
r_bm = effect_size_r(U_bm, len(exp_basemean), len(shi_basemean))

print(f"  baseMean: Exposed median={exp_basemean.median():.1f}, Shielded median={shi_basemean.median():.1f}")
print(f"  Mann-Whitney p={p_bm:.4e}, r={r_bm:.3f}")

# Expression breadth: use LFC_T2 and LFC_T3 to assess constitutive vs dynamic
# Already established: 0% exposed are constitutive vs ~38% shielded
# Use absolute LFC as proxy for expression variability
for lfc_col, label in [('LFC_T2vsT1', 'T2vsT1'), ('LFC_T3vsT1', 'T3vsT1')]:
    exp_lfc = all_reg.loc[all_reg['group']=='Exposed', lfc_col].dropna().abs()
    shi_lfc = all_reg.loc[all_reg['group']=='Shielded', lfc_col].dropna().abs()
    U_l, p_l = mannwhitneyu(exp_lfc, shi_lfc, alternative='two-sided')
    r_l = effect_size_r(U_l, len(exp_lfc), len(shi_lfc))
    print(f"  |{label}|: Exposed median={exp_lfc.median():.3f}, Shielded median={shi_lfc.median():.3f}, p={p_l:.4e}, r={r_l:.3f}")

# Constitutive gene fraction
# Define constitutive: |LFC_T2| < 0.5 AND |LFC_T3| < 0.5
all_reg['is_constitutive'] = (all_reg['LFC_T2vsT1'].abs() < 0.5) & (all_reg['LFC_T3vsT1'].abs() < 0.5)
exp_const = all_reg.loc[all_reg['group']=='Exposed', 'is_constitutive'].sum()
shi_const = all_reg.loc[all_reg['group']=='Shielded', 'is_constitutive'].sum()
exp_n = (all_reg['group']=='Exposed').sum()
shi_n = (all_reg['group']=='Shielded').sum()

print(f"\n  Constitutive genes (|LFC|<0.5 both timepoints):")
print(f"    Exposed: {exp_const}/{exp_n} ({exp_const/exp_n*100:.1f}%)")
print(f"    Shielded: {shi_const}/{shi_n} ({shi_const/shi_n*100:.1f}%)")

# ----------------------------------------------------------
# Step 8: Integrated composite conservation score
# ----------------------------------------------------------
print("\n--- Step 8: Computing composite conservation score ---")

# Metrics for composite score:
# 1. gene_length (longer = more conserved, typically) → z-score, positive = more conserved
# 2. has_gene_name (yes=1 = conserved) → binary
# 3. gc3_deviation (lower = more conserved) → z-score, NEGATE so lower = more conserved
# 4. nc (lower = more codon bias = more conserved for highly expressed) → z-score, NEGATE
# 5. baseMean (higher = more conserved) → z-score, positive
# 6. has_named_product (yes=1 = more conserved) → binary
# 7. has_old_locus_tag (yes=1 = original annotation = more conserved) → binary
# 8. rare_codon_freq (lower = more conserved in high-GC) → z-score, NEGATE

score_cols = []

# Continuous metrics: z-score normalize
for col, negate, label in [
    ('gene_length', False, 'gene_length_z'),
    ('gc3_deviation', True, 'gc3_deviation_z'),
    ('nc', True, 'nc_z'),
    ('baseMean', False, 'baseMean_z'),
    ('rare_codon_freq', True, 'rare_codon_z'),
]:
    vals = all_reg[col].copy()
    valid = vals.dropna()
    if len(valid) > 0:
        z = (vals - valid.mean()) / valid.std()
        if negate:
            z = -z
        all_reg[label] = z
        score_cols.append(label)
        print(f"  {label}: mean(z) Exposed={all_reg.loc[all_reg['group']=='Exposed', label].mean():.3f}, "
              f"Shielded={all_reg.loc[all_reg['group']=='Shielded', label].mean():.3f}")

# Binary metrics: convert to 0/1 and z-score
for col, label in [
    ('has_gene_name', 'has_gene_name_z'),
    ('has_named_product', 'has_named_product_z'),
    ('has_old_locus_tag', 'has_old_locus_tag_z'),
]:
    vals = all_reg[col].astype(float)
    z = (vals - vals.mean()) / vals.std()
    all_reg[label] = z
    score_cols.append(label)

# Compute composite score (mean of z-scores)
all_reg['composite_conservation'] = all_reg[score_cols].mean(axis=1)

exp_score = all_reg.loc[all_reg['group']=='Exposed', 'composite_conservation'].dropna()
shi_score = all_reg.loc[all_reg['group']=='Shielded', 'composite_conservation'].dropna()

U_comp, p_comp = mannwhitneyu(exp_score, shi_score, alternative='two-sided')
r_comp = effect_size_r(U_comp, len(exp_score), len(shi_score))
d_comp = cohens_d(exp_score.values, shi_score.values)

print(f"\n  Composite conservation score:")
print(f"    Exposed: mean={exp_score.mean():.4f}, median={exp_score.median():.4f}")
print(f"    Shielded: mean={shi_score.mean():.4f}, median={shi_score.median():.4f}")
print(f"    Mann-Whitney U={U_comp:.0f}, p={p_comp:.4e}")
print(f"    Effect size r={r_comp:.3f}, Cohen's d={d_comp:.3f}")

# ROC analysis: can composite score predict exposed status?
from sklearn.metrics import roc_auc_score
valid = all_reg[['composite_conservation', 'is_exposed_bool']].dropna()
if len(valid) > 0:
    # Note: lower conservation -> more likely exposed, so negate
    y_true = valid['is_exposed_bool'].astype(int).values
    y_score = -valid['composite_conservation'].values  # negate: less conserved = more likely exposed

    try:
        auc_score = roc_auc_score(y_true, y_score)
        fpr, tpr, thresholds = roc_curve(y_true, y_score)
        print(f"    ROC AUC = {auc_score:.3f}")
    except:
        auc_score = np.nan
        fpr, tpr = [0, 1], [0, 1]
        print("    ROC AUC: could not compute")
else:
    auc_score = np.nan
    fpr, tpr = [0, 1], [0, 1]

# ----------------------------------------------------------
# Compile all statistical tests
# ----------------------------------------------------------
print("\n--- Compiling statistical tests ---")

stat_tests = []

# Continuous metrics: Wilcoxon rank-sum (Mann-Whitney U)
for col, label in [
    ('gene_length', 'Gene length'),
    ('gc_content', 'GC content'),
    ('gc3', 'GC3 (third codon position)'),
    ('gc3_deviation', 'GC3 deviation from genome average'),
    ('nc', 'Effective number of codons (Nc)'),
    ('rare_codon_freq', 'Rare codon frequency'),
    ('dist_oriC', 'Distance to oriC'),
    ('baseMean', 'Base mean expression'),
]:
    exp_vals = all_reg.loc[all_reg['group']=='Exposed', col].dropna()
    shi_vals = all_reg.loc[all_reg['group']=='Shielded', col].dropna()

    if len(exp_vals) > 0 and len(shi_vals) > 0:
        U, p = mannwhitneyu(exp_vals, shi_vals, alternative='two-sided')
        r = effect_size_r(U, len(exp_vals), len(shi_vals))
        d = cohens_d(exp_vals.values, shi_vals.values)

        stat_tests.append({
            'metric': label,
            'test': 'Mann-Whitney U',
            'exposed_n': len(exp_vals),
            'exposed_median': exp_vals.median(),
            'exposed_mean': exp_vals.mean(),
            'shielded_n': len(shi_vals),
            'shielded_median': shi_vals.median(),
            'shielded_mean': shi_vals.mean(),
            'U_statistic': U,
            'p_value': p,
            'effect_size_r': r,
            'cohens_d': d,
            'direction': 'Exposed > Shielded' if exp_vals.median() > shi_vals.median() else 'Exposed < Shielded'
        })

# Fisher exact tests (from annotation_results)
for ar in annotation_results:
    stat_tests.append({
        'metric': ar['feature'],
        'test': 'Fisher exact',
        'exposed_n': ar['exposed_total'],
        'exposed_median': ar['exposed_pct'],
        'exposed_mean': ar['exposed_pct'],
        'shielded_n': ar['shielded_total'],
        'shielded_median': ar['shielded_pct'],
        'shielded_mean': ar['shielded_pct'],
        'U_statistic': ar['odds_ratio'],
        'p_value': ar['p_value'],
        'effect_size_r': np.nan,
        'cohens_d': np.nan,
        'direction': f"OR={ar['odds_ratio']:.3f}"
    })

# Core/arm test
stat_tests.append({
    'metric': 'Core enrichment',
    'test': 'Fisher exact',
    'exposed_n': exp_core + exp_arm,
    'exposed_median': exp_core / (exp_core + exp_arm) * 100,
    'exposed_mean': exp_core / (exp_core + exp_arm) * 100,
    'shielded_n': shi_core + shi_arm,
    'shielded_median': shi_core / (shi_core + shi_arm) * 100,
    'shielded_mean': shi_core / (shi_core + shi_arm) * 100,
    'U_statistic': or_region,
    'p_value': p_region,
    'effect_size_r': np.nan,
    'cohens_d': np.nan,
    'direction': f"OR={or_region:.3f}"
})

# Composite score test
stat_tests.append({
    'metric': 'Composite conservation score',
    'test': 'Mann-Whitney U',
    'exposed_n': len(exp_score),
    'exposed_median': exp_score.median(),
    'exposed_mean': exp_score.mean(),
    'shielded_n': len(shi_score),
    'shielded_median': shi_score.median(),
    'shielded_mean': shi_score.mean(),
    'U_statistic': U_comp,
    'p_value': p_comp,
    'effect_size_r': r_comp,
    'cohens_d': d_comp,
    'direction': 'Exposed < Shielded' if exp_score.median() < shi_score.median() else 'Exposed > Shielded'
})

# |LFC| tests
for lfc_col, label in [('LFC_T2vsT1', '|LFC T2vsT1|'), ('LFC_T3vsT1', '|LFC T3vsT1|')]:
    exp_vals = all_reg.loc[all_reg['group']=='Exposed', lfc_col].dropna().abs()
    shi_vals = all_reg.loc[all_reg['group']=='Shielded', lfc_col].dropna().abs()
    U, p = mannwhitneyu(exp_vals, shi_vals, alternative='two-sided')
    r = effect_size_r(U, len(exp_vals), len(shi_vals))
    d = cohens_d(exp_vals.values, shi_vals.values)
    stat_tests.append({
        'metric': label,
        'test': 'Mann-Whitney U',
        'exposed_n': len(exp_vals),
        'exposed_median': exp_vals.median(),
        'exposed_mean': exp_vals.mean(),
        'shielded_n': len(shi_vals),
        'shielded_median': shi_vals.median(),
        'shielded_mean': shi_vals.mean(),
        'U_statistic': U,
        'p_value': p,
        'effect_size_r': r,
        'cohens_d': d,
        'direction': 'Exposed > Shielded' if exp_vals.median() > shi_vals.median() else 'Exposed < Shielded'
    })

# Constitutive test
table_const = [[exp_const, exp_n - exp_const], [shi_const, shi_n - shi_const]]
or_const, p_const = fisher_exact(table_const)
stat_tests.append({
    'metric': 'Constitutive expression (|LFC|<0.5)',
    'test': 'Fisher exact',
    'exposed_n': exp_n,
    'exposed_median': exp_const / exp_n * 100,
    'exposed_mean': exp_const / exp_n * 100,
    'shielded_n': shi_n,
    'shielded_median': shi_const / shi_n * 100,
    'shielded_mean': shi_const / shi_n * 100,
    'U_statistic': or_const,
    'p_value': p_const,
    'effect_size_r': np.nan,
    'cohens_d': np.nan,
    'direction': f"OR={or_const:.3f}"
})

stat_df = pd.DataFrame(stat_tests)
stat_df.to_csv(f"{TABLES_DIR}/statistical_tests.tsv", sep='\t', index=False)
print(f"Saved: {TABLES_DIR}/statistical_tests.tsv")
print(f"\nSignificant tests (p < 0.05):")
for _, row in stat_df[stat_df['p_value'] < 0.05].iterrows():
    print(f"  {row['metric']}: p={row['p_value']:.4e}, {row['direction']}")

# ----------------------------------------------------------
# Save per-gene conservation table
# ----------------------------------------------------------
print("\n--- Saving per-gene conservation table ---")

output_cols = ['locus_tag', 'gene_name', 'old_locus_tag', 'product', 'tf_family',
               'start', 'end', 'strand', 'region', 'group',
               'gene_length', 'gc_content', 'gc3', 'gc3_deviation', 'nc', 'rare_codon_freq',
               'dist_oriC', 'region_calc', 'baseMean',
               'has_gene_name', 'has_old_locus_tag', 'has_named_product', 'is_hypothetical',
               'has_protein_id', 'is_constitutive',
               'LFC_T2vsT1', 'LFC_T3vsT1',
               'composite_conservation']

# Only include columns that exist
output_cols = [c for c in output_cols if c in all_reg.columns]
all_reg[output_cols].to_csv(f"{TABLES_DIR}/per_gene_conservation.tsv", sep='\t', index=False)
print(f"Saved: {TABLES_DIR}/per_gene_conservation.tsv")

# Codon usage table
codon_cols = ['locus_tag', 'gene_name', 'old_locus_tag', 'group',
              'gene_length', 'gc_content', 'gc3', 'gc3_deviation', 'nc', 'rare_codon_freq']
codon_cols = [c for c in codon_cols if c in all_reg.columns]
all_reg[codon_cols].to_csv(f"{TABLES_DIR}/codon_usage.tsv", sep='\t', index=False)
print(f"Saved: {TABLES_DIR}/codon_usage.tsv")

# Exposed vs shielded comparison table
comparison_data = []
for col in ['gene_length', 'gc_content', 'gc3', 'gc3_deviation', 'nc', 'rare_codon_freq',
            'dist_oriC', 'baseMean', 'composite_conservation']:
    if col in all_reg.columns:
        exp_vals = all_reg.loc[all_reg['group']=='Exposed', col].dropna()
        shi_vals = all_reg.loc[all_reg['group']=='Shielded', col].dropna()
        comparison_data.append({
            'metric': col,
            'exposed_n': len(exp_vals),
            'exposed_mean': exp_vals.mean(),
            'exposed_median': exp_vals.median(),
            'exposed_std': exp_vals.std(),
            'shielded_n': len(shi_vals),
            'shielded_mean': shi_vals.mean(),
            'shielded_median': shi_vals.median(),
            'shielded_std': shi_vals.std(),
        })

comp_df = pd.DataFrame(comparison_data)
comp_df.to_csv(f"{TABLES_DIR}/exposed_vs_shielded_comparison.tsv", sep='\t', index=False)
print(f"Saved: {TABLES_DIR}/exposed_vs_shielded_comparison.tsv")

# ============================================================
# FIGURES
# ============================================================
print("\n" + "=" * 70)
print("GENERATING FIGURES")
print("=" * 70)

# Color palette
COLOR_EXPOSED = '#E74C3C'
COLOR_SHIELDED = '#3498DB'
COLOR_PALETTE = [COLOR_SHIELDED, COLOR_EXPOSED]

def save_fig(fig, name):
    """Save figure in both PDF and SVG."""
    fig.savefig(f"{FIGURES_DIR}/{name}.pdf", dpi=300, bbox_inches='tight')
    fig.savefig(f"{FIGURES_DIR}/{name}.svg", dpi=300, bbox_inches='tight')
    print(f"  Saved: {FIGURES_DIR}/{name}.pdf/svg")
    plt.close(fig)


# ----------------------------------------------------------
# Figure 1: Conservation metrics comparison (multi-panel box plots)
# ----------------------------------------------------------
print("\n--- Figure 1: Conservation metrics comparison ---")

fig, axes = plt.subplots(2, 4, figsize=(20, 10))
fig.suptitle('Conservation Proxy Metrics: Exposed vs Shielded Regulators', fontsize=14, fontweight='bold', y=0.98)

plot_metrics = [
    ('gene_length', 'Gene Length (bp)', False),
    ('gc_content', 'GC Content', False),
    ('gc3', 'GC3 (3rd codon position)', False),
    ('gc3_deviation', 'GC3 Deviation from Genome', False),
    ('nc', 'Effective Number of Codons (Nc)', False),
    ('rare_codon_freq', 'Rare Codon Frequency', False),
    ('baseMean', 'Base Mean Expression', True),
    ('dist_oriC', 'Distance to oriC (bp)', False),
]

for idx, (col, title, log_scale) in enumerate(plot_metrics):
    ax = axes[idx // 4, idx % 4]

    data_shi = all_reg.loc[all_reg['group']=='Shielded', col].dropna()
    data_exp = all_reg.loc[all_reg['group']=='Exposed', col].dropna()

    bp = ax.boxplot([data_shi, data_exp],
                    labels=['Shielded\n(n={})'.format(len(data_shi)),
                            'Exposed\n(n={})'.format(len(data_exp))],
                    patch_artist=True, widths=0.6,
                    medianprops=dict(color='black', linewidth=2))

    bp['boxes'][0].set_facecolor(COLOR_SHIELDED)
    bp['boxes'][0].set_alpha(0.6)
    bp['boxes'][1].set_facecolor(COLOR_EXPOSED)
    bp['boxes'][1].set_alpha(0.6)

    if log_scale and data_shi.min() > 0 and data_exp.min() > 0:
        ax.set_yscale('log')

    # Add p-value
    if len(data_exp) > 0 and len(data_shi) > 0:
        U, p = mannwhitneyu(data_exp, data_shi, alternative='two-sided')
        r = effect_size_r(U, len(data_exp), len(data_shi))
        p_str = f"p={p:.2e}" if p < 0.001 else f"p={p:.4f}"
        ax.set_title(f'{title}\n{p_str}, r={r:.3f}', fontsize=10)
    else:
        ax.set_title(title, fontsize=10)

    ax.tick_params(labelsize=9)

plt.tight_layout(rect=[0, 0, 1, 0.95])
save_fig(fig, 'conservation_metrics_comparison')

# ----------------------------------------------------------
# Figure 2: Annotation quality bar charts
# ----------------------------------------------------------
print("\n--- Figure 2: Annotation quality ---")

fig, axes = plt.subplots(1, 4, figsize=(16, 5))
fig.suptitle('Annotation Quality as Conservation Proxy: Exposed vs Shielded', fontsize=13, fontweight='bold')

for idx, (col, title) in enumerate([
    ('has_gene_name', 'Has Gene Name'),
    ('has_old_locus_tag', 'Has SCO Locus Tag'),
    ('has_named_product', 'Named Product\n(non-hypothetical)'),
    ('has_protein_id', 'Has Protein ID (WP_)'),
]):
    ax = axes[idx]

    exp_pct = all_reg.loc[all_reg['group']=='Exposed', col].mean() * 100
    shi_pct = all_reg.loc[all_reg['group']=='Shielded', col].mean() * 100

    bars = ax.bar(['Shielded', 'Exposed'], [shi_pct, exp_pct],
                  color=[COLOR_SHIELDED, COLOR_EXPOSED], alpha=0.7, edgecolor='black')

    # Get Fisher test result
    match = annot_df[annot_df['feature'] == [r for r in annot_df['feature'] if col.replace('has_', 'Has ').replace('_', ' ') in r.lower() or title.split('\n')[0].lower() in r.lower()][0]] if False else None

    # Find matching result
    for _, ar in annot_df.iterrows():
        if col.replace('has_', '').replace('_', ' ') in ar['feature'].lower() or \
           title.split('\n')[0].lower().replace('has ', '') in ar['feature'].lower():
            p_str = f"p={ar['p_value']:.2e}" if ar['p_value'] < 0.001 else f"p={ar['p_value']:.4f}"
            ax.set_title(f"{title}\nOR={ar['odds_ratio']:.2f}, {p_str}", fontsize=10)
            break
    else:
        ax.set_title(title, fontsize=10)

    ax.set_ylabel('Percentage (%)')
    ax.set_ylim(0, 105)

    for bar, pct in zip(bars, [shi_pct, exp_pct]):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f'{pct:.1f}%', ha='center', va='bottom', fontsize=10)

plt.tight_layout(rect=[0, 0, 1, 0.92])
save_fig(fig, 'annotation_quality')

# ----------------------------------------------------------
# Figure 3: Chromosomal position distribution
# ----------------------------------------------------------
print("\n--- Figure 3: Chromosomal position ---")

fig, axes = plt.subplots(1, 3, figsize=(18, 5))
fig.suptitle('Chromosomal Position: Exposed vs Shielded Regulators', fontsize=13, fontweight='bold')

# Panel A: Distribution along chromosome
ax = axes[0]
shi_pos = all_reg.loc[all_reg['group']=='Shielded', 'start'] / 1e6
exp_pos = all_reg.loc[all_reg['group']=='Exposed', 'start'] / 1e6

ax.hist(shi_pos, bins=50, alpha=0.5, color=COLOR_SHIELDED, label=f'Shielded (n={len(shi_pos)})', density=True)
ax.hist(exp_pos, bins=20, alpha=0.7, color=COLOR_EXPOSED, label=f'Exposed (n={len(exp_pos)})', density=True)
ax.axvline(LEFT_ARM_END/1e6, color='gray', linestyle='--', alpha=0.7, label='Arm/Core boundary')
ax.axvline(RIGHT_ARM_START/1e6, color='gray', linestyle='--', alpha=0.7)
ax.axvline(ORIC_POSITION/1e6, color='green', linestyle=':', alpha=0.7, label='oriC')
ax.set_xlabel('Chromosomal Position (Mb)')
ax.set_ylabel('Density')
ax.set_title('Position Distribution')
ax.legend(fontsize=8)

# Panel B: Core vs arm bar chart
ax = axes[1]
exp_total = exp_core + exp_arm
shi_total = shi_core + shi_arm
data = pd.DataFrame({
    'Group': ['Shielded', 'Shielded', 'Exposed', 'Exposed'],
    'Region': ['Core', 'Arm', 'Core', 'Arm'],
    'Fraction': [shi_core/shi_total*100, shi_arm/shi_total*100,
                 exp_core/exp_total*100, exp_arm/exp_total*100],
    'Count': [shi_core, shi_arm, exp_core, exp_arm]
})

x = np.arange(2)
width = 0.35
core_vals = [shi_core/shi_total*100, exp_core/exp_total*100]
arm_vals = [shi_arm/shi_total*100, exp_arm/exp_total*100]
bars1 = ax.bar(x - width/2, core_vals, width, label='Core', color='#2ECC71', alpha=0.7, edgecolor='black')
bars2 = ax.bar(x + width/2, arm_vals, width, label='Arm', color='#F39C12', alpha=0.7, edgecolor='black')

ax.set_ylabel('Percentage (%)')
ax.set_title(f'Core vs Arm\nFisher OR={or_region:.2f}, p={p_region:.4f}')
ax.set_xticks(x)
ax.set_xticklabels(['Shielded', 'Exposed'])
ax.legend()
ax.set_ylim(0, 100)

for bars, vals, counts in [(bars1, core_vals, [shi_core, exp_core]),
                            (bars2, arm_vals, [shi_arm, exp_arm])]:
    for bar, val, cnt in zip(bars, vals, counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f'{val:.0f}%\n(n={cnt})', ha='center', va='bottom', fontsize=8)

# Panel C: Distance to oriC
ax = axes[2]
bp = ax.boxplot([shi_dist/1e6, exp_dist/1e6],
                labels=[f'Shielded\n(n={len(shi_dist)})', f'Exposed\n(n={len(exp_dist)})'],
                patch_artist=True, widths=0.6,
                medianprops=dict(color='black', linewidth=2))
bp['boxes'][0].set_facecolor(COLOR_SHIELDED)
bp['boxes'][0].set_alpha(0.6)
bp['boxes'][1].set_facecolor(COLOR_EXPOSED)
bp['boxes'][1].set_alpha(0.6)
ax.set_ylabel('Distance to oriC (Mb)')
p_str = f"p={p_dist:.2e}" if p_dist < 0.001 else f"p={p_dist:.4f}"
ax.set_title(f'Distance to oriC\n{p_str}, r={r_dist:.3f}')

plt.tight_layout(rect=[0, 0, 1, 0.92])
save_fig(fig, 'chromosomal_position')

# ----------------------------------------------------------
# Figure 4: Composite conservation score
# ----------------------------------------------------------
print("\n--- Figure 4: Composite conservation score ---")

fig, axes = plt.subplots(1, 3, figsize=(18, 5))
fig.suptitle('Composite Conservation Score: Exposed vs Shielded', fontsize=13, fontweight='bold')

# Panel A: Distribution
ax = axes[0]
bins = np.linspace(all_reg['composite_conservation'].dropna().min(),
                    all_reg['composite_conservation'].dropna().max(), 30)
ax.hist(shi_score, bins=bins, alpha=0.5, color=COLOR_SHIELDED, label=f'Shielded (n={len(shi_score)})', density=True)
ax.hist(exp_score, bins=bins, alpha=0.7, color=COLOR_EXPOSED, label=f'Exposed (n={len(exp_score)})', density=True)
ax.axvline(shi_score.median(), color=COLOR_SHIELDED, linestyle='--', linewidth=2)
ax.axvline(exp_score.median(), color=COLOR_EXPOSED, linestyle='--', linewidth=2)
p_str = f"p={p_comp:.2e}" if p_comp < 0.001 else f"p={p_comp:.4f}"
ax.set_title(f'Score Distribution\nMWU {p_str}, r={r_comp:.3f}')
ax.set_xlabel('Composite Conservation Score')
ax.set_ylabel('Density')
ax.legend(fontsize=9)

# Panel B: Box plot
ax = axes[1]
bp = ax.boxplot([shi_score, exp_score],
                labels=[f'Shielded\n(n={len(shi_score)})', f'Exposed\n(n={len(exp_score)})'],
                patch_artist=True, widths=0.6,
                medianprops=dict(color='black', linewidth=2))
bp['boxes'][0].set_facecolor(COLOR_SHIELDED)
bp['boxes'][0].set_alpha(0.6)
bp['boxes'][1].set_facecolor(COLOR_EXPOSED)
bp['boxes'][1].set_alpha(0.6)
ax.set_ylabel('Composite Conservation Score')
ax.set_title(f"Cohen's d = {d_comp:.3f}")

# Panel C: ROC curve
ax = axes[2]
if not np.isnan(auc_score):
    ax.plot(fpr, tpr, color=COLOR_EXPOSED, linewidth=2, label=f'AUC = {auc_score:.3f}')
    ax.plot([0, 1], [0, 1], color='gray', linestyle='--', alpha=0.7)
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title(f'ROC: Conservation Predicts Exposed\nAUC = {auc_score:.3f}')
    ax.legend(fontsize=11)
else:
    ax.text(0.5, 0.5, 'ROC not available', ha='center', va='center')

plt.tight_layout(rect=[0, 0, 1, 0.92])
save_fig(fig, 'composite_conservation_score')

# ----------------------------------------------------------
# Figure 5: Conservation by TF family
# ----------------------------------------------------------
print("\n--- Figure 5: Conservation by TF family ---")

# Get TF families with enough genes
family_counts = all_reg['tf_family'].value_counts()
major_families = family_counts[family_counts >= 10].index.tolist()

fig, axes = plt.subplots(1, 2, figsize=(16, 6))
fig.suptitle('Conservation Score by TF Family', fontsize=13, fontweight='bold')

# Panel A: Box plot of conservation by family
ax = axes[0]
family_data = []
family_labels = []
family_colors = []
for fam in sorted(major_families):
    fam_data = all_reg.loc[all_reg['tf_family']==fam, 'composite_conservation'].dropna()
    if len(fam_data) >= 5:
        family_data.append(fam_data.values)
        n_exp = all_reg.loc[(all_reg['tf_family']==fam) & (all_reg['group']=='Exposed')].shape[0]
        family_labels.append(f'{fam}\n(n={len(fam_data)}, E={n_exp})')
        family_colors.append(COLOR_EXPOSED if n_exp > 0 else COLOR_SHIELDED)

if family_data:
    bp = ax.boxplot(family_data, labels=family_labels, patch_artist=True,
                    medianprops=dict(color='black', linewidth=1.5))
    for patch, color in zip(bp['boxes'], family_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.5)

ax.set_ylabel('Composite Conservation Score')
ax.set_title('Conservation by TF Family (major families)')
ax.tick_params(axis='x', rotation=45, labelsize=8)
ax.axhline(0, color='gray', linestyle='--', alpha=0.5)

# Panel B: Exposed fraction vs conservation score by family
ax = axes[1]
fam_stats = []
for fam in sorted(major_families):
    fam_df = all_reg[all_reg['tf_family']==fam]
    if len(fam_df) >= 5:
        n_exp = (fam_df['group']=='Exposed').sum()
        exp_frac = n_exp / len(fam_df) * 100
        mean_cons = fam_df['composite_conservation'].mean()
        fam_stats.append({'family': fam, 'exposed_frac': exp_frac,
                         'mean_conservation': mean_cons, 'n': len(fam_df)})

if fam_stats:
    fam_stats_df = pd.DataFrame(fam_stats)
    sizes = fam_stats_df['n'] * 5
    ax.scatter(fam_stats_df['mean_conservation'], fam_stats_df['exposed_frac'],
              s=sizes, alpha=0.7, c=fam_stats_df['exposed_frac'], cmap='RdYlBu_r',
              edgecolors='black', linewidth=0.5)

    for _, row in fam_stats_df.iterrows():
        ax.annotate(row['family'], (row['mean_conservation'], row['exposed_frac']),
                   fontsize=7, ha='center', va='bottom')

    ax.set_xlabel('Mean Composite Conservation Score')
    ax.set_ylabel('% Exposed in Family')
    ax.set_title('Exposed Fraction vs Conservation by Family')

    # Correlation
    if len(fam_stats_df) > 3:
        rho, p_rho = spearmanr(fam_stats_df['mean_conservation'], fam_stats_df['exposed_frac'])
        ax.annotate(f'Spearman rho={rho:.3f}, p={p_rho:.3f}', xy=(0.05, 0.95),
                   xycoords='axes fraction', fontsize=9, va='top')

plt.tight_layout(rect=[0, 0, 1, 0.93])
save_fig(fig, 'conservation_by_TF_family')

# ----------------------------------------------------------
# Figure 6: Comprehensive summary (multi-panel)
# ----------------------------------------------------------
print("\n--- Figure 6: Comprehensive summary ---")

fig = plt.figure(figsize=(24, 18))
gs = GridSpec(3, 4, figure=fig, hspace=0.35, wspace=0.3)
fig.suptitle('H36: Evolutionary Conservation of 62 Exposed Transcription Factors\n'
             'Comprehensive Analysis Summary', fontsize=15, fontweight='bold', y=0.99)

# Panel A: Gene length
ax = fig.add_subplot(gs[0, 0])
data_s = all_reg.loc[all_reg['group']=='Shielded', 'gene_length'].dropna()
data_e = all_reg.loc[all_reg['group']=='Exposed', 'gene_length'].dropna()
bp = ax.boxplot([data_s, data_e], labels=['Shielded', 'Exposed'],
                patch_artist=True, widths=0.6, medianprops=dict(color='black', linewidth=2))
bp['boxes'][0].set_facecolor(COLOR_SHIELDED); bp['boxes'][0].set_alpha(0.6)
bp['boxes'][1].set_facecolor(COLOR_EXPOSED); bp['boxes'][1].set_alpha(0.6)
U, p = mannwhitneyu(data_e, data_s, alternative='two-sided')
r = effect_size_r(U, len(data_e), len(data_s))
ax.set_title(f'A. Gene Length\np={p:.2e}, r={r:.3f}', fontsize=10)
ax.set_ylabel('bp')

# Panel B: GC3
ax = fig.add_subplot(gs[0, 1])
data_s = all_reg.loc[all_reg['group']=='Shielded', 'gc3'].dropna()
data_e = all_reg.loc[all_reg['group']=='Exposed', 'gc3'].dropna()
bp = ax.boxplot([data_s, data_e], labels=['Shielded', 'Exposed'],
                patch_artist=True, widths=0.6, medianprops=dict(color='black', linewidth=2))
bp['boxes'][0].set_facecolor(COLOR_SHIELDED); bp['boxes'][0].set_alpha(0.6)
bp['boxes'][1].set_facecolor(COLOR_EXPOSED); bp['boxes'][1].set_alpha(0.6)
U, p = mannwhitneyu(data_e, data_s, alternative='two-sided')
r = effect_size_r(U, len(data_e), len(data_s))
ax.set_title(f'B. GC3\np={p:.2e}, r={r:.3f}', fontsize=10)
ax.set_ylabel('GC at 3rd codon position')

# Panel C: Nc
ax = fig.add_subplot(gs[0, 2])
data_s = all_reg.loc[all_reg['group']=='Shielded', 'nc'].dropna()
data_e = all_reg.loc[all_reg['group']=='Exposed', 'nc'].dropna()
bp = ax.boxplot([data_s, data_e], labels=['Shielded', 'Exposed'],
                patch_artist=True, widths=0.6, medianprops=dict(color='black', linewidth=2))
bp['boxes'][0].set_facecolor(COLOR_SHIELDED); bp['boxes'][0].set_alpha(0.6)
bp['boxes'][1].set_facecolor(COLOR_EXPOSED); bp['boxes'][1].set_alpha(0.6)
U, p = mannwhitneyu(data_e, data_s, alternative='two-sided')
r = effect_size_r(U, len(data_e), len(data_s))
ax.set_title(f'C. Effective Nc\np={p:.2e}, r={r:.3f}', fontsize=10)
ax.set_ylabel('Nc')

# Panel D: baseMean
ax = fig.add_subplot(gs[0, 3])
bp = ax.boxplot([shi_basemean, exp_basemean], labels=['Shielded', 'Exposed'],
                patch_artist=True, widths=0.6, medianprops=dict(color='black', linewidth=2))
bp['boxes'][0].set_facecolor(COLOR_SHIELDED); bp['boxes'][0].set_alpha(0.6)
bp['boxes'][1].set_facecolor(COLOR_EXPOSED); bp['boxes'][1].set_alpha(0.6)
ax.set_yscale('log')
p_str = f"p={p_bm:.2e}" if p_bm < 0.001 else f"p={p_bm:.4f}"
ax.set_title(f'D. Expression Level\n{p_str}, r={r_bm:.3f}', fontsize=10)
ax.set_ylabel('baseMean (log scale)')

# Panel E: Annotation quality grouped bar chart
ax = fig.add_subplot(gs[1, 0])
features = ['Gene\nName', 'SCO\nTag', 'Named\nProduct', 'Protein\nID']
exp_pcts = []
shi_pcts = []
for col in ['has_gene_name', 'has_old_locus_tag', 'has_named_product', 'has_protein_id']:
    exp_pcts.append(all_reg.loc[all_reg['group']=='Exposed', col].mean() * 100)
    shi_pcts.append(all_reg.loc[all_reg['group']=='Shielded', col].mean() * 100)

x = np.arange(len(features))
width = 0.35
ax.bar(x - width/2, shi_pcts, width, label='Shielded', color=COLOR_SHIELDED, alpha=0.7, edgecolor='black')
ax.bar(x + width/2, exp_pcts, width, label='Exposed', color=COLOR_EXPOSED, alpha=0.7, edgecolor='black')
ax.set_ylabel('Percentage (%)')
ax.set_title('E. Annotation Quality', fontsize=10)
ax.set_xticks(x)
ax.set_xticklabels(features, fontsize=8)
ax.legend(fontsize=8)
ax.set_ylim(0, 105)

# Panel F: Core/Arm distribution
ax = fig.add_subplot(gs[1, 1])
x = np.arange(2)
width = 0.35
core_vals = [shi_core/(shi_core+shi_arm)*100, exp_core/(exp_core+exp_arm)*100]
arm_vals = [shi_arm/(shi_core+shi_arm)*100, exp_arm/(exp_core+exp_arm)*100]
ax.bar(x - width/2, core_vals, width, label='Core', color='#2ECC71', alpha=0.7, edgecolor='black')
ax.bar(x + width/2, arm_vals, width, label='Arm', color='#F39C12', alpha=0.7, edgecolor='black')
ax.set_ylabel('%')
ax.set_title(f'F. Core vs Arm\nFisher p={p_region:.4f}', fontsize=10)
ax.set_xticks(x)
ax.set_xticklabels(['Shielded', 'Exposed'])
ax.legend(fontsize=8)
ax.set_ylim(0, 100)

# Panel G: Composite score distribution
ax = fig.add_subplot(gs[1, 2])
bins_c = np.linspace(all_reg['composite_conservation'].dropna().min(),
                      all_reg['composite_conservation'].dropna().max(), 25)
ax.hist(shi_score, bins=bins_c, alpha=0.5, color=COLOR_SHIELDED, label='Shielded', density=True)
ax.hist(exp_score, bins=bins_c, alpha=0.7, color=COLOR_EXPOSED, label='Exposed', density=True)
ax.axvline(shi_score.median(), color=COLOR_SHIELDED, linestyle='--', linewidth=2)
ax.axvline(exp_score.median(), color=COLOR_EXPOSED, linestyle='--', linewidth=2)
p_str = f"p={p_comp:.2e}" if p_comp < 0.001 else f"p={p_comp:.4f}"
ax.set_title(f'G. Composite Score\n{p_str}, d={d_comp:.3f}', fontsize=10)
ax.set_xlabel('Score')
ax.set_ylabel('Density')
ax.legend(fontsize=8)

# Panel H: ROC
ax = fig.add_subplot(gs[1, 3])
if not np.isnan(auc_score):
    ax.plot(fpr, tpr, color=COLOR_EXPOSED, linewidth=2)
    ax.plot([0, 1], [0, 1], color='gray', linestyle='--', alpha=0.7)
    ax.set_xlabel('FPR')
    ax.set_ylabel('TPR')
    ax.set_title(f'H. ROC Curve\nAUC = {auc_score:.3f}', fontsize=10)
# aspect ratio left to auto for better layout

# Panel I: Rare codon frequency
ax = fig.add_subplot(gs[2, 0])
data_s = all_reg.loc[all_reg['group']=='Shielded', 'rare_codon_freq'].dropna()
data_e = all_reg.loc[all_reg['group']=='Exposed', 'rare_codon_freq'].dropna()
bp = ax.boxplot([data_s, data_e], labels=['Shielded', 'Exposed'],
                patch_artist=True, widths=0.6, medianprops=dict(color='black', linewidth=2))
bp['boxes'][0].set_facecolor(COLOR_SHIELDED); bp['boxes'][0].set_alpha(0.6)
bp['boxes'][1].set_facecolor(COLOR_EXPOSED); bp['boxes'][1].set_alpha(0.6)
U, p = mannwhitneyu(data_e, data_s, alternative='two-sided')
r = effect_size_r(U, len(data_e), len(data_s))
ax.set_title(f'I. Rare Codon Frequency\np={p:.2e}, r={r:.3f}', fontsize=10)
ax.set_ylabel('Fraction rare codons')

# Panel J: |LFC| comparison
ax = fig.add_subplot(gs[2, 1])
data_s_t2 = all_reg.loc[all_reg['group']=='Shielded', 'LFC_T2vsT1'].dropna().abs()
data_e_t2 = all_reg.loc[all_reg['group']=='Exposed', 'LFC_T2vsT1'].dropna().abs()
data_s_t3 = all_reg.loc[all_reg['group']=='Shielded', 'LFC_T3vsT1'].dropna().abs()
data_e_t3 = all_reg.loc[all_reg['group']=='Exposed', 'LFC_T3vsT1'].dropna().abs()

positions = [1, 2, 4, 5]
bp = ax.boxplot([data_s_t2, data_e_t2, data_s_t3, data_e_t3],
                positions=positions, patch_artist=True, widths=0.6,
                medianprops=dict(color='black', linewidth=2))
colors = [COLOR_SHIELDED, COLOR_EXPOSED, COLOR_SHIELDED, COLOR_EXPOSED]
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color); patch.set_alpha(0.6)

ax.set_xticks([1.5, 4.5])
ax.set_xticklabels(['T2 vs T1', 'T3 vs T1'])
ax.set_ylabel('|log2 Fold Change|')
ax.set_title('J. Expression Variability\n(Higher = Less Conserved)', fontsize=10)

# Panel K: Chromosomal distribution density
ax = fig.add_subplot(gs[2, 2])
# KDE plot
from scipy.stats import gaussian_kde
pos_s = all_reg.loc[all_reg['group']=='Shielded', 'start'].values / 1e6
pos_e = all_reg.loc[all_reg['group']=='Exposed', 'start'].values / 1e6

kde_s = gaussian_kde(pos_s, bw_method=0.1)
kde_e = gaussian_kde(pos_e, bw_method=0.2)
x_range = np.linspace(0, CHROM_LENGTH/1e6, 200)
ax.fill_between(x_range, kde_s(x_range), alpha=0.4, color=COLOR_SHIELDED, label='Shielded')
ax.plot(x_range, kde_e(x_range), color=COLOR_EXPOSED, linewidth=2, label='Exposed')
ax.axvline(LEFT_ARM_END/1e6, color='gray', linestyle='--', alpha=0.5)
ax.axvline(RIGHT_ARM_START/1e6, color='gray', linestyle='--', alpha=0.5)
ax.axvline(ORIC_POSITION/1e6, color='green', linestyle=':', alpha=0.7, label='oriC')
ax.set_xlabel('Position (Mb)')
ax.set_ylabel('Density')
ax.set_title('K. Chromosomal Density', fontsize=10)
ax.legend(fontsize=7)

# Panel L: Summary text
ax = fig.add_subplot(gs[2, 3])
ax.axis('off')

# Count significant results
sig_tests = stat_df[stat_df['p_value'] < 0.05]
nonsig_tests = stat_df[stat_df['p_value'] >= 0.05]

summary_lines = [
    "L. SUMMARY OF FINDINGS",
    "",
    f"Regulatory genes: 1,017",
    f"  Exposed: 62 | Shielded: 955",
    "",
    f"Significant tests: {len(sig_tests)}/{len(stat_df)}",
    "",
    "KEY RESULTS:",
]

for _, row in sig_tests.iterrows():
    p_str = f"p={row['p_value']:.1e}"
    metric_short = row['metric']
    if len(metric_short) > 30:
        metric_short = metric_short[:28] + ".."
    summary_lines.append(f"  {metric_short}")
    summary_lines.append(f"    {p_str}, {row['direction']}")

summary_lines.append("")
summary_lines.append(f"Composite: d={d_comp:.3f}")
if not np.isnan(auc_score):
    summary_lines.append(f"ROC AUC: {auc_score:.3f}")

summary_lines.append("")
summary_lines.append("CONCLUSION:")
if p_comp < 0.05:
    if exp_score.median() < shi_score.median():
        summary_lines.append("Exposed TFs show LOWER")
        summary_lines.append("conservation proxies")
    else:
        summary_lines.append("Exposed TFs show HIGHER")
        summary_lines.append("conservation proxies")
else:
    summary_lines.append("No significant difference")
    summary_lines.append("in composite conservation")

summary_text = "\n".join(summary_lines)

ax.text(0.05, 0.95, summary_text, transform=ax.transAxes,
        fontsize=7.5, verticalalignment='top', fontfamily='monospace',
        linespacing=1.3,
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

plt.tight_layout(rect=[0, 0, 1, 0.95])
save_fig(fig, 'H36_comprehensive_summary')

# ============================================================
# FINAL SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("H36 ANALYSIS COMPLETE")
print("=" * 70)

print(f"\nFiles generated:")
print(f"  Tables:")
for f in sorted(os.listdir(TABLES_DIR)):
    print(f"    {f}")
print(f"  Figures:")
for f in sorted(os.listdir(FIGURES_DIR)):
    print(f"    {f}")

print(f"\n--- Summary Statistics ---")
print(f"Total regulatory genes analyzed: {len(all_reg)}")
print(f"  Exposed: {(all_reg['group']=='Exposed').sum()}")
print(f"  Shielded: {(all_reg['group']=='Shielded').sum()}")

print(f"\nComposite Conservation Score:")
print(f"  Exposed:  mean={exp_score.mean():.4f}, median={exp_score.median():.4f}")
print(f"  Shielded: mean={shi_score.mean():.4f}, median={shi_score.median():.4f}")
print(f"  Mann-Whitney p={p_comp:.4e}, Cohen's d={d_comp:.3f}")
if not np.isnan(auc_score):
    print(f"  ROC AUC = {auc_score:.3f}")

print(f"\nSignificant differences (p<0.05):")
for _, row in sig_tests.iterrows():
    print(f"  {row['metric']}: p={row['p_value']:.4e}, {row['direction']}")

print(f"\nNon-significant (p>=0.05):")
for _, row in nonsig_tests.iterrows():
    print(f"  {row['metric']}: p={row['p_value']:.4f}")

print("\nDone.")
