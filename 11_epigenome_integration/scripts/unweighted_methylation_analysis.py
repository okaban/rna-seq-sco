#!/usr/bin/env python3
"""
Unweighted Methylation Analysis (v2)
Simple consensus without coverage weighting.
Outputs both MIN_REPS=2 and MIN_REPS=3 for comparison.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# Paths
PILEUP_DIR = Path("/Users/okaban/bioinfo/methyl/260102_M145/analysis/cursor_results/20260108/pileup")
OUTPUT_DIR = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/01_integration")

# Sample mapping
SAMPLE_TIMEPOINT = {
    '1-1': 'T1', '1-2': 'T1', '1-3': 'T1',
    '2-1': 'T2', '2-3': 'T2', '2-4': 'T2',
    '3-2': 'T3', '3-3': 'T3', '3-4': 'T3'
}

# Thresholds
MIN_COVERAGE = 5       # Minimum per-sample coverage
MIN_TOTAL_COVERAGE = 30  # Minimum combined coverage across replicates
MIN_MOD_FREQ = 50.0    # Minimum average frequency (%)


def parse_pileup(filepath):
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

            if coverage >= MIN_COVERAGE:
                key = (chrom, pos, strand, mod_type)
                sites[key] = {'coverage': coverage, 'mod_freq': mod_freq}
    return sites


def calculate_consensus(timepoint_data, min_reps):
    """
    Simple unweighted consensus for a timepoint.
    mod_freq = simple mean of replicates that detected the site.
    """
    all_sites = defaultdict(list)
    for sample, sites in timepoint_data.items():
        for key, data in sites.items():
            all_sites[key].append(data)

    results = []
    for key, sample_data in all_sites.items():
        n_reps = len(sample_data)
        if n_reps >= min_reps:
            coverages = [d['coverage'] for d in sample_data]
            mod_freqs = [d['mod_freq'] for d in sample_data]
            total_coverage = sum(coverages)
            mean_freq = np.mean(mod_freqs)

            if total_coverage >= MIN_TOTAL_COVERAGE and mean_freq >= MIN_MOD_FREQ:
                results.append({
                    'chrom': key[0],
                    'position': key[1],
                    'strand': key[2],
                    'mod_type': key[3],
                    'n_reps': n_reps,
                    'total_coverage': total_coverage,
                    'mean_mod_freq': round(mean_freq, 2),
                    'min_freq': round(min(mod_freqs), 2),
                    'max_freq': round(max(mod_freqs), 2),
                })
    return pd.DataFrame(results)


def main():
    print("=" * 70)
    print("Unweighted Methylation Analysis v2")
    print("  - Simple mean of replicates (no coverage weighting)")
    print(f"  - MIN_COVERAGE={MIN_COVERAGE}, MIN_TOTAL_COVERAGE={MIN_TOTAL_COVERAGE}")
    print(f"  - MIN_MOD_FREQ={MIN_MOD_FREQ}%")
    print("=" * 70)

    # Load pileup data
    print("\n[1] Loading pileup data...")
    sample_data = {}
    for sample in SAMPLE_TIMEPOINT:
        filepath = PILEUP_DIR / f"{sample}_pileup.bed"
        print(f"  Loading {sample}...", end="", flush=True)
        sample_data[sample] = parse_pileup(filepath)
        print(f"  {len(sample_data[sample]):,} sites (cov>={MIN_COVERAGE})")

    # Calculate consensus for both MIN_REPS=2 and MIN_REPS=3
    for min_reps in [2, 3]:
        print(f"\n{'=' * 70}")
        print(f"[2] Consensus with MIN_REPS={min_reps}")
        print("=" * 70)

        all_dfs = []
        for tp in ['T1', 'T2', 'T3']:
            tp_samples = {s: sample_data[s] for s, t in SAMPLE_TIMEPOINT.items() if t == tp}
            df = calculate_consensus(tp_samples, min_reps)
            df['timepoint'] = tp
            all_dfs.append(df)

            n_6mA = len(df[df['mod_type'] == '6mA'])
            n_4mC = len(df[df['mod_type'] == '4mC'])
            print(f"  {tp}: {len(df):>5} sites  (6mA={n_6mA}, 4mC={n_4mC})")

        combined = pd.concat(all_dfs, ignore_index=True)
        total = len(combined)
        unique_pos = combined.groupby(['chrom', 'position', 'strand', 'mod_type']).ngroups
        print(f"  Total: {total} site-timepoint pairs, {unique_pos} unique positions")

        # Save
        out_cols = ['chrom', 'position', 'strand', 'mod_type', 'timepoint',
                    'n_reps', 'total_coverage', 'mean_mod_freq', 'min_freq', 'max_freq']
        out_path = OUTPUT_DIR / f"high_confidence_sites_unweighted_minreps{min_reps}.csv"
        combined[out_cols].to_csv(out_path, index=False)
        print(f"  Saved: {out_path.name}")

    # T3 sample coverage comparison
    print(f"\n{'=' * 70}")
    print("[3] T3 Sample Coverage Summary")
    print("=" * 70)
    for sample in ['3-2', '3-3', '3-4']:
        sites = sample_data[sample]
        covs = [d['coverage'] for d in sites.values()]
        print(f"  {sample}: median coverage = {np.median(covs):.0f}x, "
              f"mean = {np.mean(covs):.1f}x, "
              f"sites = {len(sites):,}")

    print(f"\n{'=' * 70}")
    print("DONE")
    print("=" * 70)


if __name__ == "__main__":
    main()
