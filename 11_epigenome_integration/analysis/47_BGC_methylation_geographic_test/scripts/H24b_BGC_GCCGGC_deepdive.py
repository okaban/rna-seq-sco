"""Deep-dive: is GCCGGC 4mC enrichment at BGCs more than a sequence-count artifact?

Tests, using per-site occupancy (methylation frequency) rather than site counts:
  1. T1 per-site occupancy of GCCGGC sites: BGC-internal vs non-BGC (genome-wide & core-only)
  2. Temporal erasure: # BGC-internal GCCGGC sites and mean occupancy at T1/T2/T3
  3. Site density BGC vs non-BGC (context)
"""
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import mannwhitneyu

BASE = Path(__file__).resolve().parents[4]
SITES = BASE / '11_epigenome_integration/analysis/37_defense_island_GCCGGC/tables/GCCGGC_sites_by_timepoint.tsv'
BGC = BASE / '11_epigenome_integration/analysis/72_antismash_genomewide/bgc_summary.tsv'

sites = pd.read_csv(SITES, sep='\t')
bgc = pd.read_csv(BGC, sep='\t')
print(f"GCCGGC site-rows: {len(sites)} | BGC regions: {len(bgc)}")
print("site columns:", list(sites.columns))

# mark BGC-internal sites by position interval
intervals = list(zip(bgc['start'], bgc['end']))
def in_bgc(pos):
    return any(s <= pos <= e for s, e in intervals)
sites['in_bgc'] = sites['position'].apply(in_bgc)

for tp in ['T1', 'T2', 'T3']:
    s = sites[sites['timepoint'] == tp]
    b = s[s['in_bgc']]
    n = s[~s['in_bgc']]
    print(f"\n=== {tp} ===")
    print(f"  GCCGGC sites: {len(b)} BGC-internal / {len(n)} non-BGC")
    if len(b) >= 3 and len(n) >= 3:
        # per-site occupancy (frequency)
        U, p = mannwhitneyu(b['frequency'], n['frequency'], alternative='two-sided')
        print(f"  per-site occupancy (frequency %): BGC median={b['frequency'].median():.1f} "
              f"vs non-BGC median={n['frequency'].median():.1f}  (MWU p={p:.3g})")
        # core-only (control geography)
        bc = b[b['region'] == 'core']; nc = n[n['region'] == 'core']
        if len(bc) >= 3 and len(nc) >= 3:
            Uc, pc = mannwhitneyu(bc['frequency'], nc['frequency'], alternative='two-sided')
            print(f"  core-only occupancy: BGC median={bc['frequency'].median():.1f} (n={len(bc)}) "
                  f"vs non-BGC median={nc['frequency'].median():.1f} (n={len(nc)})  (MWU p={pc:.3g})")

# temporal erasure of BGC-internal GCCGGC
print("\n=== Temporal: BGC-internal GCCGGC site count & mean occupancy ===")
for tp in ['T1', 'T2', 'T3']:
    b = sites[(sites['timepoint'] == tp) & (sites['in_bgc'])]
    allt = sites[sites['timepoint'] == tp]
    print(f"  {tp}: {len(b)} BGC sites (mean occ {b['frequency'].mean():.1f}%) | "
          f"genome {len(allt)} sites (mean occ {allt['frequency'].mean():.1f}%)")
