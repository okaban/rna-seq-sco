#!/usr/bin/env python3
"""
T3 Re-analysis with 2 replicates (excluding low-coverage sample 3-2)
Uses samples 3-3 and 3-4 only for T3 methylation calling
"""

import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Paths
PILEUP_DIR = "/Users/okaban/bioinfo/methyl/260102_M145/analysis/cursor_results/20260108/pileup"
ORIGINAL_SITES = "/Users/okaban/bioinfo/methyl/260102_M145/analysis/cursor_results/20260108/high_confidence_sites.csv"
OUTPUT_DIR = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis"

# Thresholds (matching original analysis)
MIN_COVERAGE = 10
MIN_MOD_FREQ = 50.0

def parse_pileup(filepath):
    """Parse pileup BED file to extract methylation sites."""
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

            if coverage >= MIN_COVERAGE and mod_freq >= MIN_MOD_FREQ:
                key = (chrom, pos, strand, mod_type)
                sites[key] = {'coverage': coverage, 'mod_freq': mod_freq}

    return sites

def find_consensus_sites(sample_sites_list, min_reps=2):
    """Find sites present in at least min_reps samples."""
    site_counts = defaultdict(list)

    for sites in sample_sites_list:
        for key, data in sites.items():
            site_counts[key].append(data)

    consensus = []
    for key, data_list in site_counts.items():
        if len(data_list) >= min_reps:
            coverages = [d['coverage'] for d in data_list]
            mod_freqs = [d['mod_freq'] for d in data_list]
            consensus.append({
                'chrom': key[0],
                'position': key[1],
                'strand': key[2],
                'mod_type': key[3],
                'mean_coverage': np.mean(coverages),
                'std_coverage': np.std(coverages),
                'mean_mod_freq': np.mean(mod_freqs),
                'std_mod_freq': np.std(mod_freqs),
                'min_mod_freq': min(mod_freqs),
                'max_mod_freq': max(mod_freqs),
                'n_reps': len(data_list)
            })

    return pd.DataFrame(consensus)

def main():
    print("=" * 60)
    print("T3 Re-analysis: 2 Replicates (3-3 and 3-4 only)")
    print("=" * 60)

    # Load original data for comparison
    print("\n[1] Loading original 3-replicate results...")
    original_df = pd.read_csv(ORIGINAL_SITES)
    original_T3 = original_df[original_df['timepoint'] == 'T3']
    print(f"Original T3 sites (3 reps): {len(original_T3)}")
    print(f"  - 6mA: {len(original_T3[original_T3['mod_type'] == '6mA'])}")
    print(f"  - 4mC: {len(original_T3[original_T3['mod_type'] == '4mC'])}")

    # Parse T3 pileup files (excluding 3-2)
    print("\n[2] Parsing T3 pileup files (3-3 and 3-4 only)...")
    t3_samples = ['3-3', '3-4']
    t3_sites = []

    for sample in t3_samples:
        filepath = Path(PILEUP_DIR) / f"{sample}_pileup.bed"
        print(f"  Processing {sample}...")
        sites = parse_pileup(filepath)
        print(f"    Sites passing threshold: {len(sites)}")
        t3_sites.append(sites)

    # Find consensus sites (present in both replicates)
    print("\n[3] Finding consensus sites (2/2 replicates)...")
    new_T3_df = find_consensus_sites(t3_sites, min_reps=2)
    new_T3_df['timepoint'] = 'T3'

    print(f"New T3 sites (2 reps): {len(new_T3_df)}")
    print(f"  - 6mA: {len(new_T3_df[new_T3_df['mod_type'] == '6mA'])}")
    print(f"  - 4mC: {len(new_T3_df[new_T3_df['mod_type'] == '4mC'])}")

    # Comparison
    print("\n" + "=" * 60)
    print("COMPARISON: Original (3 reps) vs New (2 reps)")
    print("=" * 60)

    comparison = {
        'Metric': ['Total sites', '6mA sites', '4mC sites', 'Mean coverage', 'Mean mod_freq'],
        'Original (3 reps)': [
            len(original_T3),
            len(original_T3[original_T3['mod_type'] == '6mA']),
            len(original_T3[original_T3['mod_type'] == '4mC']),
            f"{original_T3['mean_coverage'].mean():.1f}",
            f"{original_T3['mean_mod_freq'].mean():.1f}%"
        ],
        'New (2 reps)': [
            len(new_T3_df),
            len(new_T3_df[new_T3_df['mod_type'] == '6mA']),
            len(new_T3_df[new_T3_df['mod_type'] == '4mC']),
            f"{new_T3_df['mean_coverage'].mean():.1f}" if len(new_T3_df) > 0 else "N/A",
            f"{new_T3_df['mean_mod_freq'].mean():.1f}%" if len(new_T3_df) > 0 else "N/A"
        ],
        'Change': [
            f"+{len(new_T3_df) - len(original_T3)} ({(len(new_T3_df)/len(original_T3)-1)*100:+.1f}%)" if len(original_T3) > 0 else "N/A",
            f"+{len(new_T3_df[new_T3_df['mod_type'] == '6mA']) - len(original_T3[original_T3['mod_type'] == '6mA'])}",
            f"+{len(new_T3_df[new_T3_df['mod_type'] == '4mC']) - len(original_T3[original_T3['mod_type'] == '4mC'])}",
            "-",
            "-"
        ]
    }

    comp_df = pd.DataFrame(comparison)
    print(comp_df.to_string(index=False))

    # Save new T3 sites
    output_path = Path(OUTPUT_DIR)
    new_T3_df.to_csv(output_path / 'T3_2rep_methylation_sites.csv', index=False)

    # Create combined dataset (T1, T2 original + T3 new)
    print("\n[4] Creating combined dataset with new T3...")
    t1_t2 = original_df[original_df['timepoint'].isin(['T1', 'T2'])]
    combined_df = pd.concat([t1_t2, new_T3_df], ignore_index=True)
    combined_df.to_csv(output_path / 'high_confidence_sites_T3_2rep.csv', index=False)

    print(f"\nCombined dataset saved:")
    print(f"  - T1: {len(combined_df[combined_df['timepoint'] == 'T1'])} sites")
    print(f"  - T2: {len(combined_df[combined_df['timepoint'] == 'T2'])} sites")
    print(f"  - T3 (2 reps): {len(combined_df[combined_df['timepoint'] == 'T3'])} sites")

    # Summary statistics
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    increase_total = len(new_T3_df) - len(original_T3)
    increase_6mA = len(new_T3_df[new_T3_df['mod_type'] == '6mA']) - len(original_T3[original_T3['mod_type'] == '6mA'])
    increase_4mC = len(new_T3_df[new_T3_df['mod_type'] == '4mC']) - len(original_T3[original_T3['mod_type'] == '4mC'])

    print(f"Removing low-coverage sample 3-2 results in:")
    print(f"  - Total T3 sites: {len(original_T3)} -> {len(new_T3_df)} ({increase_total:+d}, {increase_total/len(original_T3)*100:+.1f}%)")
    print(f"  - 6mA sites: {increase_6mA:+d}")
    print(f"  - 4mC sites: {increase_4mC:+d}")

    return new_T3_df, combined_df

if __name__ == "__main__":
    new_T3_df, combined_df = main()
