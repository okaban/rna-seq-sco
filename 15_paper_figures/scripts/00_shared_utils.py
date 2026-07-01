"""
Shared utilities for paper figure generation.
S. coelicolor A3(2) M145 Nanopore methylome analysis.
"""

import sys
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import SeqIO

# ── Base paths ──────────────────────────────────────────────────────────────
# Session-agnostic resolution: prefer first existing path
def _resolve_base():
    import os
    candidates = [
        '/sessions/busy-eloquent-bell/mnt/rna-seq',
        '/sessions/zen-modest-davinci/mnt/rna-seq',
    ]
    for c in candidates:
        if os.path.isdir(c):
            return Path(c)
    # fallback: use the directory two levels above this file (scripts/ → 15_paper_figures/ → rna-seq/)
    return Path(__file__).resolve().parents[2]

BASE = _resolve_base()
EPIGENOME = BASE / '11_epigenome_integration' / 'analysis'
METHYL = BASE
FIG_DIR = BASE / '15_paper_figures' / 'figures' / 'main'
FIG_SUP_DIR = BASE / '15_paper_figures' / 'figures' / 'supplementary'
TABLE_SUP_DIR = BASE / '15_paper_figures' / 'tables' / 'supplementary'

# ── Colors ──────────────────────────────────────────────────────────────────
COL_4mC = '#C26B6B'  # muted rose (Tol-like)
COL_6mA = '#4477AA'  # muted blue (Tol)
COL_BOTH = '#7E57C2'
COL_GRAY = '#BBBBBB'  # neutral grey
COL_DARK = '#37474F'
COL_GREEN = '#43A047'
COL_ORANGE = '#FB8C00'

# Gatekeeper model colors
COL_ACTIVATION = '#43A047'    # Green — activation bloc
COL_REPRESSION = '#FB8C00'    # Orange — repression bloc
COL_EXPOSED = '#7E57C2'       # Purple — exposed TFs
COL_SHIELDED = '#B0BEC5'      # Gray — shielded TFs
COL_ARTIFACT = '#FFB74D'      # Light orange — Simpson's artifact
COL_SIGNAL = '#1565C0'        # Blue — genuine signal

# ── Timepoint labels ────────────────────────────────────────────────────
TP_LABELS = ['T1 (12 h)', 'T2 (24 h)', 'T3 (50 h)']
TP_LABELS_NL = ['T1\n(12 h)', 'T2\n(24 h)', 'T3\n(50 h)']

MOTIF_COLORS = {
    'CCGG': '#E53935',
    'GGCCGG': '#C62828',
    'TGGCCGGC': '#B71C1C',
    'AAGCCCG': '#7E57C2',
    'GATC': '#1565C0',
    'CCGKCA': '#0D47A1',
    'unassigned': '#B0BEC5',
    'other': '#78909C',
}

# ── Style ───────────────────────────────────────────────────────────────────
STYLE = {
    'font.family': 'Arial',
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
    'axes.linewidth': 1.0,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'pdf.fonttype': 42,      # TrueType (editable in Illustrator)
    'ps.fonttype': 42,
    'svg.fonttype': 'none',  # Text as text, not paths
}


def apply_style():
    """Apply publication-quality matplotlib style."""
    plt.rcParams.update(STYLE)
    warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')


def add_panel_label(ax, label, x=-0.12, y=1.08, fontsize=14):
    """Add bold panel label (a, b, c, ...) to axes."""
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=fontsize, fontweight='bold', va='top', ha='left')


def save_figure(fig, path, formats=('pdf', 'svg')):
    """Save figure in multiple formats."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        out = path.with_suffix(f'.{fmt}')
        fig.savefig(out, format=fmt, dpi=300, bbox_inches='tight')
        print(f'  Saved: {out}')
    plt.close(fig)


# ── Data loaders ────────────────────────────────────────────────────────────

def load_methylation_census():
    """Load 4mC and 6mA census data.

    Returns:
        df_4mc, df_6ma: DataFrames with columns:
            chrom, position, strand, mod_type, timepoint, frequency,
            sequence, center_base, at_AAGCCCG, final_motif
    """
    path_4mc = EPIGENOME / '23_expanded_motif_search' / '4mC_final_census.csv'
    path_6ma = EPIGENOME / '23_expanded_motif_search' / '6mA_final_census.csv'
    df_4mc = pd.read_csv(path_4mc)
    df_6ma = pd.read_csv(path_6ma)
    return df_4mc, df_6ma


def load_methylation_hc_all():
    """Load all HC methylation sites with timepoint info (including overlaps).

    Returns:
        DataFrame with columns:
            chrom, position, strand, mod_type, timepoint,
            n_reps, total_coverage, weighted_mod_freq, unweighted_mod_freq
    """
    path = EPIGENOME / '01_integration' / 'high_confidence_sites_weighted.csv'
    df = pd.read_csv(path)
    print(f'  HC all sites: {len(df)} rows '
          f'(4mC: {(df.mod_type=="4mC").sum()}, 6mA: {(df.mod_type=="6mA").sum()})')
    return df


def load_methylation_unique_positions():
    """Load unique methylation positions (collapsed across timepoints).

    Returns:
        df_unique: DataFrame with columns: position, mod_type, final_motif, strand
    """
    df_4mc, df_6ma = load_methylation_census()
    df_all = pd.concat([df_4mc, df_6ma], ignore_index=True)
    df_unique = (df_all.drop_duplicates(subset=['position', 'mod_type'])
                 [['chrom', 'position', 'strand', 'mod_type', 'final_motif']]
                 .sort_values('position')
                 .reset_index(drop=True))
    return df_unique


def load_tss_jeong2016():
    """Load TSS data, filtered to Jeong2016 dRNA-seq experimental TSS only.

    Returns:
        DataFrame with 2,703 genes (Jeong2016 experimental TSS)
    """
    path = EPIGENOME / '18_tss_analyses' / 'comprehensive_tss_table.csv'
    df = pd.read_csv(path)
    df_j = df[df['tss_source'] == 'Jeong2016_dRNA-seq'].copy()
    print(f'  Jeong2016 TSS loaded: {len(df_j)} genes')
    return df_j


def load_reference_genome():
    """Load reference genome as SeqRecord.

    Returns:
        SeqRecord (NC_003888.3, 8,667,507 bp)
    """
    path = METHYL / 'data' / 'ref.fa'
    record = SeqIO.read(path, 'fasta')
    print(f'  Genome: {record.id}, {len(record.seq):,} bp')
    return record


def load_reference_gbk():
    """Load reference GenBank file for gene annotations.

    Returns:
        SeqRecord with features
    """
    path = METHYL / 'data' / 'ref.gbk'
    record = SeqIO.read(path, 'genbank')
    print(f'  GenBank: {record.id}, {len(record.features)} features')
    return record


def parse_meme_pwm(meme_file):
    """Parse MEME meme.txt file and extract position weight matrices.

    Args:
        meme_file: Path to meme.txt

    Returns:
        dict: {motif_name: {'pwm': pd.DataFrame, 'nsites': int, 'evalue': str, 'width': int}}
    """
    meme_file = Path(meme_file)
    text = meme_file.read_text()

    motifs = {}
    # Find all motif blocks
    pattern = (r'Motif (\S+) MEME-(\d+) position-specific probability matrix\n'
               r'-+\n'
               r'letter-probability matrix: alength= 4 w= (\d+) nsites= (\d+) E= (\S+)\s*\n'
               r'((?:\s*[\d.]+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s*\n)+)')
    for m in re.finditer(pattern, text):
        name = m.group(1)
        rank = int(m.group(2))
        width = int(m.group(3))
        nsites = int(m.group(4))
        evalue = m.group(5)
        rows = []
        for line in m.group(6).strip().split('\n'):
            vals = [float(x) for x in line.split()]
            rows.append(vals)
        pwm = pd.DataFrame(rows, columns=['A', 'C', 'G', 'T'])
        motifs[f'MEME-{rank}'] = {
            'name': name,
            'pwm': pwm,
            'nsites': nsites,
            'evalue': evalue,
            'width': width,
        }

    print(f'  Parsed {len(motifs)} motifs from {meme_file.name}')
    return motifs


def load_tf_binding_sites():
    """Load TF binding site master table.

    Returns:
        DataFrame with BS_source-based tier classification added.
    """
    path = (BASE / '13_TF_binding-site' / 'analysis' /
            '01_master_TF_list_260206_v1' /
            'master_TF_binding_sites_M145_ARCHIVED_260224.tsv')
    df = pd.read_csv(path, sep='\t')

    # Tier classification
    tier1_sources = ['RegPrecise', 'ZorroAranda2022_Curated_Strong', 'FIMO_curated_motif']
    df['tier'] = df['BS_source'].apply(lambda x: 1 if x in tier1_sources else 2)

    n1 = (df['tier'] == 1).sum()
    n2 = (df['tier'] == 2).sum()
    print(f'  TF BS loaded: {len(df)} total (Tier 1: {n1}, Tier 2: {n2})')
    return df


def load_mtase_expression():
    """Load MTase genes with expression data.

    Returns:
        DataFrame with 22 MTase genes and DESeq2 expression data.
    """
    path = EPIGENOME / '11_rm_system_identification' / 'mtase_genes_with_expression.csv'
    df = pd.read_csv(path)
    print(f'  MTase genes: {len(df)}')
    return df


def load_rebase_conservation():
    """Load REBASE motif conservation matrix.

    Returns:
        DataFrame with conservation rates per motif.
    """
    path = EPIGENOME / '21_genuswide_motif_conservation' / 'rebase_motif_conservation_matrix.csv'
    df = pd.read_csv(path)
    return df


def load_5mc_4mc_data():
    """Load 5mC vs 4mC comparison data at CCGG sites.

    Returns:
        df_detailed, df_summary: detailed per-sample and genome-wide summary
    """
    base = EPIGENOME / '24_5mC_vs_4mC_CCGG' / 'tables'
    df_det = pd.read_csv(base / 'CCGG_4mC_vs_5mC_detailed.csv')
    df_sum = pd.read_csv(base / 'genome_wide_4mC_vs_5mC_summary.csv')
    return df_det, df_sum


def load_motif_summary():
    """Load comprehensive motif summary.

    Returns:
        DataFrame with motif, mod, sites, m145_oe, rebase, etc.
    """
    path = EPIGENOME / '23_expanded_motif_search' / 'comprehensive_motif_summary.csv'
    df = pd.read_csv(path)
    return df


def load_integrated_expression():
    """Load integrated methylation-expression table.

    Returns:
        DataFrame with 8,083 genes.
    """
    path = (EPIGENOME / 'archive' / 'v1_weighted_minreps2' /
            '01_integration' / 'integrated_methyl_expression.csv')
    df = pd.read_csv(path)
    return df


# ── Gatekeeper model data loaders ──────────────────────────────────────────

def load_spatial_profile():
    """Load TSS methylation gradient spatial profile data (analysis/48).

    Returns:
        DataFrame with columns: site_type, timepoint, gene_category,
        bin_center_bp, density_per_kb_per_gene, ci_low, ci_high, n_genes, n_sites
    """
    path = EPIGENOME / '48_TSS_methylation_gradient' / 'tables' / 'spatial_profile_data.tsv'
    df = pd.read_csv(path, sep='\t')
    print(f'  Spatial profile: {len(df)} rows')
    return df


def load_all_genes_features():
    """Load shielded/exposed classification for all regulatory genes (analysis/52).

    Returns:
        DataFrame with 1,017 regulatory genes and features:
        locus_tag, gene_name, old_locus_tag, product, tf_family, start, end, strand,
        tss, region, is_exposed, baseMean, nearest_methyl_distance,
        n_methyl_sites_2kb, gene_length, LFC_T2vsT1, LFC_T3vsT1, n_FIMO_hits
    """
    path = EPIGENOME / '52_shielded_exposed_boundary' / 'tables' / 'all_genes_features.tsv'
    df = pd.read_csv(path, sep='\t')
    print(f'  All genes features: {len(df)} genes ({df["is_exposed"].sum()} exposed)')
    return df


def load_exposed_regulators():
    """Load 62 exposed TF full annotation (analysis/51).

    Returns:
        DataFrame with 62 exposed TFs and comprehensive annotation.
    """
    path = EPIGENOME / '51_exposed_regulators_characteristics' / 'tables' / 'exposed_regulators_full_table.tsv'
    df = pd.read_csv(path, sep='\t')
    print(f'  Exposed regulators: {len(df)} TFs')
    return df


def load_coexpression_matrix():
    """Load 62×62 co-expression matrix (analysis/55).

    Returns:
        DataFrame (62×62) with locus_tags as both index and columns.
    """
    path = EPIGENOME / '55_exposed_regulatory_module' / 'tables' / 'coexpression_matrix.tsv'
    df = pd.read_csv(path, sep='\t', index_col=0)
    print(f'  Co-expression matrix: {df.shape}')
    return df


def load_module_summary():
    """Load regulatory module summary (analysis/55).

    Returns:
        DataFrame with module IDs, members, evidence scores.
    """
    path = EPIGENOME / '55_exposed_regulatory_module' / 'tables' / 'module_summary.tsv'
    df = pd.read_csv(path, sep='\t')
    print(f'  Modules: {len(df)}')
    return df


def load_temporal_classification():
    """Load temporal dynamics classification for 62 exposed TFs (analysis/57).

    Returns:
        DataFrame with bloc (activation/repression/unassigned), phase_ratio,
        temporal_class, T1/T2/T3 z-scores, LFC values.
    """
    path = EPIGENOME / '57_temporal_dynamics_exposed_TF' / 'tables' / 'temporal_classification.tsv'
    df = pd.read_csv(path, sep='\t')
    print(f'  Temporal classification: {len(df)} TFs '
          f'(act={sum(df["bloc"]=="activation")}, rep={sum(df["bloc"]=="repression")})')
    return df


def load_bloc_family_enrichment():
    """Load TF family enrichment by bloc (analysis/58).

    Returns:
        DataFrame with tf_family, activation_count, repression_count,
        odds_ratio, p_value, enriched_in, functional_category, p_adj.
    """
    path = EPIGENOME / '58_exposed_TF_functional_prediction' / 'tables' / 'bloc_family_enrichment.tsv'
    df = pd.read_csv(path, sep='\t')
    df = df.dropna(subset=['tf_family'])
    print(f'  Bloc family enrichment: {len(df)} families')
    return df


def load_bloc_comparison():
    """Load bloc-level summary statistics (analysis/57).

    Returns:
        DataFrame with activation/repression aggregate metrics.
    """
    path = EPIGENOME / '57_temporal_dynamics_exposed_TF' / 'tables' / 'bloc_comparison.tsv'
    df = pd.read_csv(path, sep='\t')
    return df


def load_geographic_redistribution():
    """Load GCCGGC geographic redistribution data (analysis/37).

    Returns:
        df_summary: core/arm distribution by timepoint
        df_sites: individual GCCGGC site positions by timepoint
    """
    base = EPIGENOME / '37_defense_island_GCCGGC' / 'tables'
    df_summary = pd.read_csv(base / 'geographic_distribution_summary.tsv', sep='\t')
    df_sites = pd.read_csv(base / 'GCCGGC_sites_by_timepoint.tsv', sep='\t')
    print(f'  GCCGGC sites: T1={len(df_sites[df_sites.timepoint=="T1"])}, '
          f'T2={len(df_sites[df_sites.timepoint=="T2"])}, '
          f'T3={len(df_sites[df_sites.timepoint=="T3"])}')
    return df_summary, df_sites


def load_gene_methylation_transitions(motif='GCCGGC'):
    """Load gene-level methylation transition data (analysis/42 or 44).

    Args:
        motif: 'GCCGGC' (analysis/42) or 'AAGCCCG' (analysis/44)

    Returns:
        df_transitions: per-gene methylation transition and expression data
        df_geographic: geographic stratified test results
    """
    if motif == 'GCCGGC':
        base = EPIGENOME / '42_GCCGGC_temporal_derepression' / 'tables'
    elif motif == 'AAGCCCG':
        base = EPIGENOME / '44_AAGCCCG_temporal_causality' / 'tables'
    else:
        raise ValueError(f'Unknown motif: {motif}')
    df_trans = pd.read_csv(base / 'gene_methylation_transitions.tsv', sep='\t')
    df_geo = pd.read_csv(base / 'geographic_stratified_tests.tsv', sep='\t')
    print(f'  {motif} transitions: {len(df_trans)} genes')
    return df_trans, df_geo


def load_roc_analysis(analysis='boundary'):
    """Load ROC analysis results (analysis/52 or 53).

    Args:
        analysis: 'boundary' (analysis/52) or 'sequence' (analysis/53)

    Returns:
        DataFrame with feature, AUC, CI, threshold, etc.
    """
    if analysis == 'boundary':
        path = EPIGENOME / '52_shielded_exposed_boundary' / 'tables' / 'ROC_analysis.tsv'
    elif analysis == 'sequence':
        path = EPIGENOME / '53_TSS_sequence_determinants' / 'tables' / 'ROC_analysis.tsv'
    else:
        raise ValueError(f'Unknown analysis: {analysis}')
    df = pd.read_csv(path, sep='\t')
    print(f'  ROC ({analysis}): {len(df)} features')
    return df


def load_expression_quintile():
    """Load expression quintile analysis (analysis/52).

    Returns:
        DataFrame with quintile, n_genes, n_exposed, frac_exposed, etc.
    """
    path = EPIGENOME / '52_shielded_exposed_boundary' / 'tables' / 'expression_quintile.tsv'
    df = pd.read_csv(path, sep='\t')
    return df


def load_tfbs_spatial_profile():
    """Load TF BS methylation spatial profile (analysis/49).

    Returns:
        DataFrame with methylation density around TF binding sites.
    """
    path = EPIGENOME / '49_TF_BS_methylation_protection' / 'tables' / 'spatial_profile_data.tsv'
    df = pd.read_csv(path, sep='\t')
    print(f'  TFBS spatial profile: {len(df)} rows')
    return df


def load_tcs_pairs():
    """Load TCS pair analysis (analysis/55).

    Returns:
        DataFrame with 7 TCS pairs showing exposed/shielded asymmetry.
    """
    path = EPIGENOME / '55_exposed_regulatory_module' / 'tables' / 'TCS_pairs_analysis.tsv'
    df = pd.read_csv(path, sep='\t')
    print(f'  TCS pairs: {len(df)}')
    return df


def load_neighborhood_results():
    """Load neighborhood effect analysis (analysis/56).

    Returns:
        df_permutation: permutation test results
        df_stats: statistical test results
    """
    base = EPIGENOME / '56_exposed_TF_neighborhood' / 'tables'
    df_perm = pd.read_csv(base / 'permutation_results.tsv', sep='\t')
    df_stats = pd.read_csv(base / 'statistical_tests.tsv', sep='\t')
    return df_perm, df_stats


def load_exposed_methylation_status():
    """Load AAGCCCG/GCCGGC methylation status for 62 exposed TFs (analysis/58_causal).

    Returns:
        DataFrame with per-TF AAGCCCG and GCCGGC methylation counts.
    """
    path = EPIGENOME / '58_AAGCCCG_exposed_TF_causal' / 'tables' / 'exposed_TF_methylation_status.tsv'
    df = pd.read_csv(path, sep='\t')
    print(f'  Exposed methylation status: {len(df)} TFs')
    return df


# ── Utility functions ───────────────────────────────────────────────────────

def mm_to_inch(mm):
    """Convert millimeters to inches."""
    return mm / 25.4


def format_pvalue(p):
    """Format p-value for display."""
    if pd.isna(p):
        return 'N/A'
    if p < 0.001:
        return f'{p:.1e}'
    if p < 0.01:
        return f'{p:.3f}'
    if p < 0.05:
        return f'{p:.2f}'
    return f'{p:.2f}'


if __name__ == '__main__':
    apply_style()
    print('=== Testing shared utilities ===')
    print()

    print('Loading methylation census...')
    df_4mc, df_6ma = load_methylation_census()
    print(f'  4mC: {len(df_4mc)} rows, 6mA: {len(df_6ma)} rows')

    print('Loading unique positions...')
    df_uniq = load_methylation_unique_positions()
    print(f'  Unique: {len(df_uniq)} positions')

    print('Loading Jeong2016 TSS...')
    df_tss = load_tss_jeong2016()

    print('Loading reference genome...')
    ref = load_reference_genome()

    print('Loading MEME 4mC...')
    m4 = parse_meme_pwm(EPIGENOME / 'archive' / 'v1_weighted_minreps2' /
                         '07_motif_analysis' / 'meme_4mC' / 'meme.txt')
    for k, v in m4.items():
        print(f'  {k}: {v["name"]} (w={v["width"]}, n={v["nsites"]}, E={v["evalue"]})')

    print('Loading MEME 6mA...')
    m6 = parse_meme_pwm(EPIGENOME / 'archive' / 'v1_weighted_minreps2' /
                         '07_motif_analysis' / 'meme_6mA' / 'meme.txt')
    for k, v in m6.items():
        print(f'  {k}: {v["name"]} (w={v["width"]}, n={v["nsites"]}, E={v["evalue"]})')

    print()
    print('=== All utilities OK ===')
