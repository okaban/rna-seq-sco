#!/usr/bin/env python3
"""
Null model simulations for methylation pattern validation.
Tests:
1. -10 box methylation depletion: Is it explained by sequence composition alone?
2. SARP zone protection: Is BGC methylation absence significant vs random expectation?
"""

import pandas as pd
import numpy as np
from pathlib import Path
from Bio import SeqIO
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# Paths
BASE = Path('/Users/okaban/bioinfo/rna-seq')
GFF_PATH = BASE / 'reference/GCF_000203835.1/genomic.gff'
GENOME_FASTA = Path('/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/GCF_000203835.1_ASM20383v1_genomic.fna')
METHYL_SITES = BASE / '11_epigenome_integration/analysis/01_integration/high_confidence_sites_weighted.csv'
TSS_TABLE = BASE / '11_epigenome_integration/analysis/18_tss_analyses/comprehensive_tss_table.csv'
BGC_DEF = BASE / '05_annotation/analysis/05_annotation_260128_v1/tables/BGC_definition_manual.tsv'
OUT_DIR = BASE / '11_epigenome_integration/analysis/18_tss_analyses'
OUT_DIR.mkdir(parents=True, exist_ok=True)

N_SIMULATIONS = 10000

def load_genome():
    for rec in SeqIO.parse(GENOME_FASTA, 'fasta'):
        if 'NC_003888' in rec.id:
            return str(rec.seq).upper()
    return None

def load_tss_data():
    tss = pd.read_csv(TSS_TABLE)
    tss = tss.dropna(subset=['tss'])
    tss['tss'] = tss['tss'].astype(int)
    return tss

def load_methylation():
    return pd.read_csv(METHYL_SITES)

def count_motifs_in_region(genome_seq, start, end):
    """Count CCGG and AAGCCCG motifs in a region."""
    if start < 0:
        start = 0
    if end > len(genome_seq):
        end = len(genome_seq)
    region = genome_seq[start:end]
    ccgg = region.count('CCGG') + region.count('CCGG'[::-1].translate(str.maketrans('ACGT','TGCA')))
    aagcccg = region.count('AAGCCCG') + region.count('CCCGCTT')
    return ccgg, aagcccg

def get_gc_content(genome_seq, start, end):
    if start < 0:
        start = 0
    if end > len(genome_seq):
        end = len(genome_seq)
    region = genome_seq[start:end]
    if len(region) == 0:
        return 0
    gc = sum(1 for b in region if b in 'GC')
    return gc / len(region)


# =============================================================
# TEST 1: -10 box methylation depletion vs sequence composition
# =============================================================
def test_minus10_depletion(genome_seq, tss_df, methyl_df):
    print("=" * 60)
    print("TEST 1: -10 box methylation depletion")
    print("=" * 60)

    # Define regions relative to TSS
    regions = {
        'minus10': (-20, 0),       # -10 box region
        'minus35': (-45, -25),     # -35 box region
        'control': (-100, -60),    # control region
        'upstream': (-200, -100),  # far upstream
    }

    # For each gene with experimental TSS, extract region sequences
    exp_tss = tss_df[tss_df['tss_source'] == 'Jeong2016_dRNA-seq'].copy()
    print(f"Using {len(exp_tss)} genes with experimental TSS")

    # Count actual motifs in each region (sequence composition)
    region_motif_counts = {r: {'ccgg': 0, 'aagcccg': 0, 'total_bp': 0} for r in regions}
    region_gc = {r: [] for r in regions}

    for _, row in exp_tss.iterrows():
        tss_pos = int(row['tss'])
        strand = row['strand']
        for rname, (rel_start, rel_end) in regions.items():
            if strand == '+':
                abs_start = tss_pos + rel_start
                abs_end = tss_pos + rel_end
            else:
                abs_start = tss_pos - rel_end
                abs_end = tss_pos - rel_start
            ccgg, aagcccg = count_motifs_in_region(genome_seq, abs_start, abs_end)
            region_motif_counts[rname]['ccgg'] += ccgg
            region_motif_counts[rname]['aagcccg'] += aagcccg
            region_motif_counts[rname]['total_bp'] += (rel_end - rel_start)
            region_gc[rname].append(get_gc_content(genome_seq, abs_start, abs_end))

    print("\n--- Motif density by region (per kb) ---")
    print(f"{'Region':<15} {'CCGG/kb':>10} {'AAGCCCG/kb':>12} {'GC%':>8}")
    for rname in regions:
        total_kb = region_motif_counts[rname]['total_bp'] / 1000
        ccgg_density = region_motif_counts[rname]['ccgg'] / total_kb
        aagcccg_density = region_motif_counts[rname]['aagcccg'] / total_kb
        gc_mean = np.mean(region_gc[rname]) * 100
        print(f"{rname:<15} {ccgg_density:>10.2f} {aagcccg_density:>12.2f} {gc_mean:>7.1f}%")

    # Count actual methylation sites per region
    region_methyl_counts = {r: {'6mA': 0, '4mC': 0} for r in regions}

    for _, row in exp_tss.iterrows():
        tss_pos = int(row['tss'])
        strand = row['strand']
        gene_id = row['gene_id']
        for rname, (rel_start, rel_end) in regions.items():
            if strand == '+':
                abs_start = tss_pos + rel_start
                abs_end = tss_pos + rel_end
            else:
                abs_start = tss_pos - rel_end
                abs_end = tss_pos - rel_start
            sites = methyl_df[
                (methyl_df['position'] >= abs_start) &
                (methyl_df['position'] < abs_end)
            ]
            for mt in ['6mA', '4mC']:
                region_methyl_counts[rname][mt] += len(sites[sites['mod_type'] == mt])

    print("\n--- Actual methylation site counts ---")
    print(f"{'Region':<15} {'6mA':>8} {'4mC':>8} {'Total':>8} {'Density/kb':>12}")
    for rname in regions:
        total = region_methyl_counts[rname]['6mA'] + region_methyl_counts[rname]['4mC']
        total_kb = region_motif_counts[rname]['total_bp'] / 1000
        density = total / total_kb if total_kb > 0 else 0
        print(f"{rname:<15} {region_methyl_counts[rname]['6mA']:>8} {region_methyl_counts[rname]['4mC']:>8} {total:>8} {density:>12.3f}")

    # Simulation: given the motif density in -10 region, what's the expected methylation?
    # Null model: methylation probability per motif is constant across regions
    print(f"\n--- Null model simulation (n={N_SIMULATIONS}) ---")

    # Calculate genome-wide methylation rate per motif
    # For CCGG: what fraction of genomic CCGGs are methylated?
    total_ccgg_genome = genome_seq.count('CCGG')
    total_4mC = len(methyl_df[methyl_df['mod_type'] == '4mC'].drop_duplicates(['position']))
    # Rough: unique 4mC positions / total CCGG sites
    methyl_rate_ccgg = total_4mC / (total_ccgg_genome * 2)  # *2 for both strands approx

    # For each region, expected methylation = motif_count * methyl_rate
    print(f"\nGenome-wide: {total_ccgg_genome} CCGG sites, {total_4mC} unique 4mC positions")
    print(f"Rough methylation rate per CCGG: {methyl_rate_ccgg:.4f}")

    # More precise: simulate by randomly placing methylation sites
    # keeping the total number fixed, but distributing based on motif counts
    control_density = (region_methyl_counts['control']['6mA'] + region_methyl_counts['control']['4mC'])
    control_bp = region_motif_counts['control']['total_bp']
    minus10_density_obs = (region_methyl_counts['minus10']['6mA'] + region_methyl_counts['minus10']['4mC'])
    minus10_bp = region_motif_counts['minus10']['total_bp']

    # Motif-based null: expected methylation proportional to motif count
    control_motifs = region_motif_counts['control']['ccgg'] + region_motif_counts['control']['aagcccg']
    minus10_motifs = region_motif_counts['minus10']['ccgg'] + region_motif_counts['minus10']['aagcccg']

    total_methyl_both = control_density + minus10_density_obs
    if control_motifs + minus10_motifs > 0:
        expected_ratio = minus10_motifs / (control_motifs + minus10_motifs)
    else:
        expected_ratio = 0.5

    # Binomial simulation
    simulated_minus10 = np.random.binomial(
        total_methyl_both, expected_ratio, N_SIMULATIONS
    )
    p_value_motif = np.mean(simulated_minus10 <= minus10_density_obs)

    print(f"\nMotif-based null model:")
    print(f"  Control region motifs (CCGG+AAGCCCG): {control_motifs}")
    print(f"  -10 region motifs (CCGG+AAGCCCG): {minus10_motifs}")
    print(f"  Expected ratio for -10: {expected_ratio:.4f}")
    print(f"  Observed -10 methylation: {minus10_density_obs}")
    print(f"  Expected -10 methylation (mean): {np.mean(simulated_minus10):.1f}")
    print(f"  P-value (observed ≤ simulated): {p_value_motif:.6f}")

    # Also do length-based null
    expected_ratio_len = minus10_bp / (control_bp + minus10_bp)
    simulated_minus10_len = np.random.binomial(
        total_methyl_both, expected_ratio_len, N_SIMULATIONS
    )
    p_value_len = np.mean(simulated_minus10_len <= minus10_density_obs)

    print(f"\nLength-based null model:")
    print(f"  Expected ratio for -10 (by bp): {expected_ratio_len:.4f}")
    print(f"  Expected -10 methylation (mean): {np.mean(simulated_minus10_len):.1f}")
    print(f"  P-value (observed ≤ simulated): {p_value_len:.6f}")

    # Key question: does motif-based null EXPLAIN the depletion?
    if p_value_motif > 0.05:
        verdict_motif = "EXPLAINED by sequence composition"
    else:
        verdict_motif = "NOT fully explained by sequence composition"

    if p_value_len < 0.05:
        verdict_len = "Significant depletion vs length-based null"
    else:
        verdict_len = "Not significant vs length-based null"

    print(f"\n--- VERDICT ---")
    print(f"  Length-based null: {verdict_len} (p={p_value_len:.6f})")
    print(f"  Motif-based null: {verdict_motif} (p={p_value_motif:.6f})")

    # Figure
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    # Panel A: motif density by region
    region_names = list(regions.keys())
    ccgg_densities = [region_motif_counts[r]['ccgg'] / (region_motif_counts[r]['total_bp']/1000) for r in region_names]
    aagcccg_densities = [region_motif_counts[r]['aagcccg'] / (region_motif_counts[r]['total_bp']/1000) for r in region_names]

    x = np.arange(len(region_names))
    w = 0.35
    axes[0].bar(x - w/2, ccgg_densities, w, label='CCGG', color='#4ECDC4')
    axes[0].bar(x + w/2, aagcccg_densities, w, label='AAGCCCG', color='#FF6B6B')
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(['-10 box', '-35 box', 'Control', 'Upstream'], fontsize=9)
    axes[0].set_ylabel('Motif density (per kb)')
    axes[0].set_title('A. Target motif density by region')
    axes[0].legend(fontsize=8)

    # Panel B: observed vs expected methylation
    obs_densities = []
    exp_densities_motif = []
    for rname in region_names:
        total_kb = region_motif_counts[rname]['total_bp'] / 1000
        obs = (region_methyl_counts[rname]['6mA'] + region_methyl_counts[rname]['4mC']) / total_kb
        obs_densities.append(obs)
        motif_total = region_motif_counts[rname]['ccgg'] + region_motif_counts[rname]['aagcccg']
        # Expected based on control methylation rate per motif
        if control_motifs > 0:
            rate = control_density / control_motifs
            exp = (motif_total * rate) / total_kb
        else:
            exp = 0
        exp_densities_motif.append(exp)

    axes[1].bar(x - w/2, obs_densities, w, label='Observed', color='#2C3E50')
    axes[1].bar(x + w/2, exp_densities_motif, w, label='Expected\n(motif-based)', color='#95A5A6', alpha=0.7)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(['-10 box', '-35 box', 'Control', 'Upstream'], fontsize=9)
    axes[1].set_ylabel('Methylation density (sites/kb)')
    axes[1].set_title('B. Observed vs motif-based expected')
    axes[1].legend(fontsize=8)

    # Panel C: simulation distribution
    axes[2].hist(simulated_minus10, bins=30, color='#95A5A6', alpha=0.7, label='Motif-based null')
    axes[2].axvline(minus10_density_obs, color='#E74C3C', lw=2, label=f'Observed ({minus10_density_obs})')
    axes[2].axvline(np.mean(simulated_minus10), color='#3498DB', lw=1.5, ls='--',
                     label=f'Expected ({np.mean(simulated_minus10):.1f})')
    axes[2].set_xlabel('Methylation sites in -10 region')
    axes[2].set_ylabel('Simulation count')
    axes[2].set_title(f'C. Null model (n={N_SIMULATIONS})\np={p_value_motif:.4f}')
    axes[2].legend(fontsize=8)

    plt.tight_layout()
    fig.savefig(OUT_DIR / 'null_model_minus10_box.png', dpi=150, bbox_inches='tight')
    fig.savefig(OUT_DIR / 'null_model_minus10_box.pdf', bbox_inches='tight')
    plt.close()
    print(f"Saved: null_model_minus10_box.png/pdf")

    return {
        'p_motif': p_value_motif,
        'p_length': p_value_len,
        'observed': minus10_density_obs,
        'expected_motif': np.mean(simulated_minus10),
        'expected_length': np.mean(simulated_minus10_len),
        'minus10_motifs': minus10_motifs,
        'control_motifs': control_motifs,
        'region_motif_counts': region_motif_counts,
        'region_methyl_counts': region_methyl_counts,
    }


# =============================================================
# TEST 2: SARP zone protection
# =============================================================
def test_sarp_zone_protection(tss_df, methyl_df):
    print("\n" + "=" * 60)
    print("TEST 2: SARP zone protection in BGC genes")
    print("=" * 60)

    # Load BGC definitions
    bgc = pd.read_csv(BGC_DEF, sep='\t')
    bgc_genes = set(bgc['gene_id'].dropna().values)
    print(f"BGC genes: {len(bgc_genes)}")

    # Get genes with experimental TSS
    exp_tss = tss_df[tss_df['tss_source'] == 'Jeong2016_dRNA-seq'].copy()

    # For each gene, check if SARP zone (-60 to -20 from TSS) has methylation
    def has_sarp_zone_methyl(row, methyl_sites):
        tss_pos = int(row['tss'])
        strand = row['strand']
        if strand == '+':
            zone_start = tss_pos - 60
            zone_end = tss_pos - 20
        else:
            zone_start = tss_pos + 20
            zone_end = tss_pos + 60
        sites = methyl_sites[
            (methyl_sites['position'] >= zone_start) &
            (methyl_sites['position'] < zone_end)
        ]
        return len(sites) > 0

    # All unique methylation positions
    methyl_unique = methyl_df.drop_duplicates(['position', 'mod_type'])

    # Check all genes with exp TSS
    results = []
    for _, row in exp_tss.iterrows():
        lt = row['gene_id']
        is_bgc = lt in bgc_genes
        has_methyl = has_sarp_zone_methyl(row, methyl_unique)
        results.append({'gene_id': lt, 'is_bgc': is_bgc, 'has_sarp_methyl': has_methyl})

    res_df = pd.DataFrame(results)

    # Contingency table
    bgc_with_methyl = res_df[(res_df['is_bgc']) & (res_df['has_sarp_methyl'])].shape[0]
    bgc_without_methyl = res_df[(res_df['is_bgc']) & (~res_df['has_sarp_methyl'])].shape[0]
    nonbgc_with_methyl = res_df[(~res_df['is_bgc']) & (res_df['has_sarp_methyl'])].shape[0]
    nonbgc_without_methyl = res_df[(~res_df['is_bgc']) & (~res_df['has_sarp_methyl'])].shape[0]

    total_bgc = bgc_with_methyl + bgc_without_methyl
    total_nonbgc = nonbgc_with_methyl + nonbgc_without_methyl

    print(f"\n--- Contingency Table (exp TSS genes only) ---")
    print(f"{'':>20} {'SARP methyl+':>15} {'SARP methyl-':>15} {'Total':>10}")
    print(f"{'BGC':>20} {bgc_with_methyl:>15} {bgc_without_methyl:>15} {total_bgc:>10}")
    print(f"{'Non-BGC':>20} {nonbgc_with_methyl:>15} {nonbgc_without_methyl:>15} {total_nonbgc:>10}")

    # Fisher's exact test
    table = [[bgc_with_methyl, bgc_without_methyl],
             [nonbgc_with_methyl, nonbgc_without_methyl]]
    odds_ratio, p_fisher = stats.fisher_exact(table, alternative='less')
    print(f"\nFisher's exact test (BGC depleted?):")
    print(f"  Odds ratio: {odds_ratio:.4f}")
    print(f"  P-value (one-sided, less): {p_fisher:.6f}")

    # Genome-wide methylation rate in SARP zone
    genome_rate = (nonbgc_with_methyl + bgc_with_methyl) / (total_bgc + total_nonbgc)
    expected_bgc_methyl = genome_rate * total_bgc
    print(f"\nGenome-wide SARP zone methylation rate: {genome_rate:.4f} ({genome_rate*100:.1f}%)")
    print(f"Expected BGC genes with SARP methyl: {expected_bgc_methyl:.1f}")
    print(f"Observed BGC genes with SARP methyl: {bgc_with_methyl}")

    # Binomial test
    p_binom = stats.binomtest(bgc_with_methyl, total_bgc, genome_rate, alternative='less').pvalue
    print(f"Binomial test p-value: {p_binom:.6f}")

    # Simulation: randomly sample total_bgc genes from all exp TSS genes
    # and count how many have SARP zone methylation
    all_methyl_status = res_df['has_sarp_methyl'].values
    simulated_counts = np.zeros(N_SIMULATIONS)
    for i in range(N_SIMULATIONS):
        sampled = np.random.choice(all_methyl_status, size=total_bgc, replace=False)
        simulated_counts[i] = np.sum(sampled)

    p_sim = np.mean(simulated_counts <= bgc_with_methyl)
    print(f"\nPermutation test (n={N_SIMULATIONS}):")
    print(f"  Mean simulated: {np.mean(simulated_counts):.1f}")
    print(f"  Observed: {bgc_with_methyl}")
    print(f"  P-value: {p_sim:.6f}")

    # Now also test with ALL genes (not just exp TSS)
    print("\n--- Extended analysis: all genes ---")
    results_all = []
    for _, row in tss_df.dropna(subset=['tss']).iterrows():
        lt = row['gene_id']
        is_bgc = lt in bgc_genes
        has_methyl = has_sarp_zone_methyl(row, methyl_unique)
        results_all.append({'gene_id': lt, 'is_bgc': is_bgc, 'has_sarp_methyl': has_methyl})

    res_all_df = pd.DataFrame(results_all)
    bgc_m_all = res_all_df[(res_all_df['is_bgc']) & (res_all_df['has_sarp_methyl'])].shape[0]
    bgc_nm_all = res_all_df[(res_all_df['is_bgc']) & (~res_all_df['has_sarp_methyl'])].shape[0]
    nonbgc_m_all = res_all_df[(~res_all_df['is_bgc']) & (res_all_df['has_sarp_methyl'])].shape[0]
    nonbgc_nm_all = res_all_df[(~res_all_df['is_bgc']) & (~res_all_df['has_sarp_methyl'])].shape[0]

    total_bgc_all = bgc_m_all + bgc_nm_all
    total_nonbgc_all = nonbgc_m_all + nonbgc_nm_all
    genome_rate_all = (nonbgc_m_all + bgc_m_all) / (total_bgc_all + total_nonbgc_all)

    print(f"{'':>20} {'SARP methyl+':>15} {'SARP methyl-':>15} {'Total':>10}")
    print(f"{'BGC':>20} {bgc_m_all:>15} {bgc_nm_all:>15} {total_bgc_all:>10}")
    print(f"{'Non-BGC':>20} {nonbgc_m_all:>15} {nonbgc_nm_all:>15} {total_nonbgc_all:>10}")

    table_all = [[bgc_m_all, bgc_nm_all], [nonbgc_m_all, nonbgc_nm_all]]
    or_all, p_all = stats.fisher_exact(table_all, alternative='less')
    print(f"Fisher p-value (all genes): {p_all:.6f}, OR={or_all:.4f}")
    print(f"Genome-wide rate: {genome_rate_all:.4f} ({genome_rate_all*100:.1f}%)")
    print(f"Expected BGC: {genome_rate_all * total_bgc_all:.1f}, Observed: {bgc_m_all}")

    # Figure
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

    # Panel A: Observed vs expected
    categories = ['BGC genes', 'Non-BGC genes']
    observed = [bgc_with_methyl / total_bgc * 100 if total_bgc > 0 else 0,
                nonbgc_with_methyl / total_nonbgc * 100 if total_nonbgc > 0 else 0]
    expected = [genome_rate * 100, genome_rate * 100]

    x = np.arange(len(categories))
    axes[0].bar(x - 0.2, observed, 0.35, label='Observed', color=['#E74C3C', '#3498DB'])
    axes[0].bar(x + 0.2, expected, 0.35, label='Expected\n(genome-wide rate)',
                color='#95A5A6', alpha=0.7)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(categories)
    axes[0].set_ylabel('% genes with SARP zone methylation')
    axes[0].set_title(f'A. SARP zone methylation rate\n(Fisher p={p_fisher:.4f})')
    axes[0].legend(fontsize=8)
    for i, (o, e) in enumerate(zip(observed, expected)):
        axes[0].text(i - 0.2, o + 0.3, f'{o:.1f}%', ha='center', fontsize=9)

    # Panel B: Permutation distribution
    axes[1].hist(simulated_counts, bins=30, color='#95A5A6', alpha=0.7)
    axes[1].axvline(bgc_with_methyl, color='#E74C3C', lw=2,
                     label=f'Observed BGC ({bgc_with_methyl})')
    axes[1].axvline(np.mean(simulated_counts), color='#3498DB', lw=1.5, ls='--',
                     label=f'Expected ({np.mean(simulated_counts):.1f})')
    axes[1].set_xlabel('# BGC genes with SARP zone methylation')
    axes[1].set_ylabel('Simulation count')
    axes[1].set_title(f'B. Permutation test (n={N_SIMULATIONS})\np={p_sim:.4f}')
    axes[1].legend(fontsize=8)

    plt.tight_layout()
    fig.savefig(OUT_DIR / 'null_model_sarp_zone.png', dpi=150, bbox_inches='tight')
    fig.savefig(OUT_DIR / 'null_model_sarp_zone.pdf', bbox_inches='tight')
    plt.close()
    print(f"Saved: null_model_sarp_zone.png/pdf")

    return {
        'p_fisher': p_fisher,
        'p_binom': p_binom,
        'p_sim': p_sim,
        'bgc_with_methyl': bgc_with_methyl,
        'total_bgc': total_bgc,
        'genome_rate': genome_rate,
        'expected': expected_bgc_methyl,
    }


def main():
    print("Loading data...")
    genome_seq = load_genome()
    print(f"Genome length: {len(genome_seq):,} bp")
    tss_df = load_tss_data()
    print(f"Genes with TSS: {len(tss_df)}")
    methyl_df = load_methylation()
    print(f"Methylation sites: {len(methyl_df)}")

    # Test 1
    result1 = test_minus10_depletion(genome_seq, tss_df, methyl_df)

    # Test 2
    result2 = test_sarp_zone_protection(tss_df, methyl_df)

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"\n-10 box depletion:")
    print(f"  Length-based null p = {result1['p_length']:.6f}")
    print(f"  Motif-based null p = {result1['p_motif']:.6f}")
    if result1['p_motif'] > 0.05:
        print(f"  → Depletion is EXPLAINED by low motif density (sequence composition)")
    else:
        print(f"  → Depletion is NOT fully explained by motif density; additional mechanism")

    print(f"\nSARP zone BGC protection:")
    print(f"  Fisher p = {result2['p_fisher']:.6f}")
    print(f"  Permutation p = {result2['p_sim']:.6f}")
    print(f"  Observed/Expected = {result2['bgc_with_methyl']}/{result2['expected']:.1f}")
    if result2['p_sim'] < 0.05:
        print(f"  → BGC SARP zones are significantly depleted of methylation")
    else:
        print(f"  → BGC SARP zone methylation absence is NOT significant (random expectation)")

    # Save results
    summary = pd.DataFrame([
        {'test': 'minus10_length_null', 'p_value': result1['p_length'],
         'observed': result1['observed'], 'expected': result1['expected_length']},
        {'test': 'minus10_motif_null', 'p_value': result1['p_motif'],
         'observed': result1['observed'], 'expected': result1['expected_motif']},
        {'test': 'sarp_zone_fisher', 'p_value': result2['p_fisher'],
         'observed': result2['bgc_with_methyl'], 'expected': result2['expected']},
        {'test': 'sarp_zone_permutation', 'p_value': result2['p_sim'],
         'observed': result2['bgc_with_methyl'], 'expected': result2['expected']},
    ])
    summary.to_csv(OUT_DIR / 'null_model_results.csv', index=False)
    print(f"\nSaved: null_model_results.csv")


if __name__ == '__main__':
    main()
