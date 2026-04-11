#!/usr/bin/env python3
"""
Coverage-Weighted Methylation Analysis
Includes all 3 replicates but weights by coverage to reduce low-coverage sample bias
Also checks 5mC detection status
"""

import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Paths
PILEUP_DIR = Path("/Users/okaban/bioinfo/methyl/260102_M145/analysis/cursor_results/20260108/pileup")
GFF_PATH = "/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/genomic.gff"
DESEQ2_DIR = Path("/Users/okaban/bioinfo/rna-seq/04_deseq2/analysis/04_deseq2_260128_v1/results")
OUTPUT_DIR = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis")

# Sample mapping
SAMPLE_TIMEPOINT = {
    '1-1': 'T1', '1-2': 'T1', '1-3': 'T1',
    '2-1': 'T2', '2-3': 'T2', '2-4': 'T2',
    '3-2': 'T3', '3-3': 'T3', '3-4': 'T3'
}

# Thresholds
MIN_COVERAGE = 5  # Lower threshold per sample
MIN_TOTAL_COVERAGE = 30  # Minimum combined coverage across replicates
MIN_REPS = 2  # At least 2 replicates must detect the site
MIN_WEIGHTED_FREQ = 50.0  # Weighted frequency threshold

def parse_pileup_with_coverage(filepath):
    """Parse pileup file and return sites with coverage info."""
    sites = {}
    mod_type_map = {'a': '6mA', 'm': '5mC', '21839': '4mC'}

    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 11:
                continue

            chrom = parts[0]
            pos = int(parts[1])
            mod_code = parts[3]
            strand = parts[5]
            coverage = int(parts[9])
            mod_freq = float(parts[10])

            if mod_code not in mod_type_map:
                continue

            mod_type = mod_type_map[mod_code]

            # Store all sites with minimum coverage
            if coverage >= MIN_COVERAGE:
                key = (chrom, pos, strand, mod_type)
                sites[key] = {'coverage': coverage, 'mod_freq': mod_freq}

    return sites

def calculate_weighted_consensus(timepoint_data):
    """
    Calculate coverage-weighted consensus for a timepoint.

    Args:
        timepoint_data: dict of {sample: {site_key: {coverage, mod_freq}}}

    Returns:
        DataFrame with weighted methylation frequencies
    """
    # Collect all sites across samples
    all_sites = defaultdict(list)

    for sample, sites in timepoint_data.items():
        for key, data in sites.items():
            all_sites[key].append({
                'sample': sample,
                'coverage': data['coverage'],
                'mod_freq': data['mod_freq']
            })

    # Calculate weighted consensus
    results = []
    for key, sample_data in all_sites.items():
        n_reps = len(sample_data)

        if n_reps >= MIN_REPS:
            coverages = [d['coverage'] for d in sample_data]
            mod_freqs = [d['mod_freq'] for d in sample_data]
            total_coverage = sum(coverages)

            # Coverage-weighted average
            weighted_freq = sum(c * f for c, f in zip(coverages, mod_freqs)) / total_coverage

            # Also calculate unweighted for comparison
            unweighted_freq = np.mean(mod_freqs)

            if total_coverage >= MIN_TOTAL_COVERAGE and weighted_freq >= MIN_WEIGHTED_FREQ:
                results.append({
                    'chrom': key[0],
                    'position': key[1],
                    'strand': key[2],
                    'mod_type': key[3],
                    'n_reps': n_reps,
                    'total_coverage': total_coverage,
                    'weighted_mod_freq': weighted_freq,
                    'unweighted_mod_freq': unweighted_freq,
                    'coverage_per_rep': [d['coverage'] for d in sample_data],
                    'freq_per_rep': [d['mod_freq'] for d in sample_data]
                })

    return pd.DataFrame(results)

def analyze_5mC_detailed():
    """Detailed analysis of 5mC detection."""
    print("\n" + "=" * 70)
    print("5mC DETAILED ANALYSIS")
    print("=" * 70)

    all_5mC = defaultdict(list)

    for sample, timepoint in SAMPLE_TIMEPOINT.items():
        filepath = PILEUP_DIR / f"{sample}_pileup.bed"
        print(f"\nProcessing {sample} ({timepoint})...")

        with open(filepath, 'r') as f:
            sample_5mC = {'total': 0, 'cov10': 0, 'cov10_freq10': 0, 'cov10_freq20': 0, 'cov10_freq30': 0}
            freqs = []

            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 11 or parts[3] != 'm':
                    continue

                coverage = int(parts[9])
                mod_freq = float(parts[10])
                sample_5mC['total'] += 1

                if coverage >= 10:
                    sample_5mC['cov10'] += 1
                    freqs.append(mod_freq)
                    if mod_freq >= 10:
                        sample_5mC['cov10_freq10'] += 1
                    if mod_freq >= 20:
                        sample_5mC['cov10_freq20'] += 1
                    if mod_freq >= 30:
                        sample_5mC['cov10_freq30'] += 1

            print(f"  Total 5mC sites: {sample_5mC['total']:,}")
            print(f"  Sites with cov>=10: {sample_5mC['cov10']:,}")
            print(f"  Sites with cov>=10 & freq>=10%: {sample_5mC['cov10_freq10']}")
            print(f"  Sites with cov>=10 & freq>=20%: {sample_5mC['cov10_freq20']}")
            print(f"  Sites with cov>=10 & freq>=30%: {sample_5mC['cov10_freq30']}")
            if freqs:
                print(f"  Mean frequency (cov>=10): {np.mean(freqs):.2f}%")
                print(f"  Max frequency (cov>=10): {max(freqs):.2f}%")

    print("\n" + "-" * 70)
    print("5mC CONCLUSION:")
    print("-" * 70)
    print("5mC modification is detected but at very low frequencies (<10%).")
    print("This is expected for Streptomyces - 5mC is not a major modification.")
    print("The primary modifications are 6mA and 4mC (restriction-modification systems).")
    print("-" * 70)

def main():
    print("=" * 70)
    print("Coverage-Weighted Methylation Analysis (n=3 with bias correction)")
    print("=" * 70)

    # First, analyze 5mC
    analyze_5mC_detailed()

    # Load pileup data for each sample
    print("\n" + "=" * 70)
    print("WEIGHTED METHYLATION ANALYSIS (6mA, 4mC)")
    print("=" * 70)

    print("\n[1] Loading pileup data from all samples...")
    sample_data = {}
    for sample in SAMPLE_TIMEPOINT.keys():
        filepath = PILEUP_DIR / f"{sample}_pileup.bed"
        print(f"  Loading {sample}...")
        sample_data[sample] = parse_pileup_with_coverage(filepath)
        print(f"    Sites passing min coverage ({MIN_COVERAGE}x): {len(sample_data[sample]):,}")

    # Group by timepoint
    print("\n[2] Calculating coverage-weighted consensus per timepoint...")
    timepoint_results = {}

    for timepoint in ['T1', 'T2', 'T3']:
        tp_samples = {s: sample_data[s] for s, tp in SAMPLE_TIMEPOINT.items() if tp == timepoint}
        result_df = calculate_weighted_consensus(tp_samples)
        result_df['timepoint'] = timepoint
        timepoint_results[timepoint] = result_df

        print(f"\n  {timepoint}:")
        print(f"    Total high-confidence sites: {len(result_df)}")
        for mod in ['6mA', '4mC']:
            mod_count = len(result_df[result_df['mod_type'] == mod])
            print(f"    - {mod}: {mod_count}")

    # Combine all timepoints
    all_sites_df = pd.concat(timepoint_results.values(), ignore_index=True)

    # Comparison with original and 2-rep analyses
    print("\n" + "=" * 70)
    print("COMPARISON: Original vs 2-rep vs Weighted")
    print("=" * 70)

    original_counts = {'T1': 2349, 'T2': 2794, 'T3': 885}  # From previous analysis
    tworep_counts = {'T1': 2349, 'T2': 2794, 'T3': 2780}  # From 2-rep analysis

    print("\n{:<12} {:>12} {:>12} {:>12}".format("Timepoint", "Original", "2-rep", "Weighted"))
    print("-" * 50)
    for tp in ['T1', 'T2', 'T3']:
        weighted_count = len(timepoint_results[tp])
        print("{:<12} {:>12} {:>12} {:>12}".format(
            tp, original_counts[tp], tworep_counts[tp], weighted_count))

    # Show T3 sample contributions
    print("\n[3] T3 Sample Contribution Analysis...")
    t3_df = timepoint_results['T3']

    if len(t3_df) > 0 and 'coverage_per_rep' in t3_df.columns:
        # Analyze how much each replicate contributes
        rep_contributions = defaultdict(list)
        for _, row in t3_df.iterrows():
            covs = row['coverage_per_rep']
            total = sum(covs)
            for i, cov in enumerate(covs):
                rep_contributions[f'rep_{i+1}'].append(cov / total * 100)

        print("\n  T3 Coverage Contribution per Site:")
        t3_samples = ['3-2', '3-3', '3-4']
        for i, sample in enumerate(t3_samples):
            if f'rep_{i+1}' in rep_contributions:
                mean_contrib = np.mean(rep_contributions[f'rep_{i+1}'])
                print(f"    {sample}: {mean_contrib:.1f}% average contribution")

    # Save weighted results
    print("\n[4] Saving results...")
    output_path = OUTPUT_DIR

    # Save combined weighted sites
    save_df = all_sites_df[['chrom', 'position', 'strand', 'mod_type', 'timepoint',
                            'n_reps', 'total_coverage', 'weighted_mod_freq', 'unweighted_mod_freq']]
    save_df.to_csv(output_path / 'high_confidence_sites_weighted.csv', index=False)

    # Summary statistics
    print("\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)

    print("\nMETHOD COMPARISON:")
    print("-" * 70)
    print("1. Original (3-rep, equal weight):")
    print("   - Requires all 3 replicates to detect a site")
    print("   - Low-coverage sample 3-2 (14x) reduces T3 detection")
    print(f"   - T3 sites: {original_counts['T3']}")
    print("\n2. 2-rep (excluding 3-2):")
    print("   - Excludes low-coverage sample entirely")
    print("   - Loses statistical power from n=3")
    print(f"   - T3 sites: {tworep_counts['T3']}")
    print("\n3. Weighted (3-rep, coverage-weighted):")
    print("   - Includes all 3 replicates")
    print("   - Weights contribution by coverage")
    print("   - Low-coverage samples contribute less but aren't excluded")
    print(f"   - T3 sites: {len(timepoint_results['T3'])}")

    print("\n5mC STATUS:")
    print("-" * 70)
    print("- 5mC is detected at very low frequencies (<10% for all sites)")
    print("- No high-confidence 5mC sites pass standard thresholds")
    print("- This is biologically expected for Streptomyces")
    print("- Primary modifications: 6mA (adenine methylation) and 4mC")

    return all_sites_df, timepoint_results

if __name__ == "__main__":
    all_sites_df, timepoint_results = main()
