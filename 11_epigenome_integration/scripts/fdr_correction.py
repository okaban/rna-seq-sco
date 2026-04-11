#!/usr/bin/env python3
"""
Apply Benjamini-Hochberg FDR correction to all correlation and enrichment tests.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

BASE = Path('/Users/okaban/bioinfo/rna-seq')
TSS_DIR = BASE / '11_epigenome_integration/analysis/18_tss_analyses'

def correct_distance_correlations():
    """FDR correction for C1 distance-stratified correlations."""
    print("=" * 60)
    print("1. C1 Distance-stratified correlations")
    print("=" * 60)

    df = pd.read_csv(TSS_DIR / 'C1_distance_stratified_correlation.csv')
    print(f"Total tests: {len(df)}")

    # Apply BH correction
    reject, pvals_corrected, _, _ = multipletests(df['p_value'], method='fdr_bh')
    df['p_adj'] = pvals_corrected
    df['significant_fdr'] = reject

    # Show significant results
    sig = df[df['significant_fdr']]
    print(f"Significant after FDR: {len(sig)}")
    if len(sig) > 0:
        print(sig[['mod_type', 'comparison', 'distance_bin', 'spearman_r', 'p_value', 'p_adj']].to_string(index=False))
    else:
        print("  No tests survive FDR correction")

    # Show borderline (p_adj < 0.1)
    borderline = df[(df['p_adj'] < 0.1) & (~df['significant_fdr'])]
    if len(borderline) > 0:
        print(f"\nBorderline (FDR < 0.1): {len(borderline)}")
        print(borderline[['mod_type', 'comparison', 'distance_bin', 'spearman_r', 'p_value', 'p_adj']].to_string(index=False))

    df.to_csv(TSS_DIR / 'C1_distance_correlation_fdr.csv', index=False)
    return df


def correct_sarp_zone():
    """FDR for SARP zone correlations."""
    print("\n" + "=" * 60)
    print("2. SARP zone correlations")
    print("=" * 60)

    df = pd.read_csv(TSS_DIR / 'sarp_zone_correlation.csv')
    print(f"Total tests: {len(df)}")

    if len(df) > 1:
        reject, pvals_corrected, _, _ = multipletests(df['p_value'], method='fdr_bh')
        df['p_adj'] = pvals_corrected
        df['significant_fdr'] = reject
    else:
        df['p_adj'] = df['p_value']
        df['significant_fdr'] = df['p_value'] < 0.05

    sig = df[df['significant_fdr']]
    print(f"Significant after FDR: {len(sig)}")
    print(df[['mod_type', 'comparison', 'spearman_r', 'p_value', 'p_adj']].to_string(index=False))
    return df


def correct_gained_lost_motifs():
    """FDR for gained/lost motif Fisher tests."""
    print("\n" + "=" * 60)
    print("3. Gained/Lost motif Fisher tests")
    print("=" * 60)

    df = pd.read_csv(TSS_DIR / 'temporal_gained_lost_motifs.csv')

    # Fisher tests were done per comparison × mod_type × motif
    # Check if p-values exist
    if 'fisher_p' not in df.columns:
        # Need to reconstruct from the data
        print("No fisher_p column found. Computing from gained/lost percentages...")
        # The script computed Fisher's tests and printed them but may not have saved p-values
        # Let's recompute
        from scipy.stats import fisher_exact

        results = []
        for comp in df['comparison'].unique():
            for mod in df['mod_type'].unique():
                for motif in df['motif'].unique():
                    if motif == 'CCGG_rc' or motif == 'CCCGCTT':
                        continue
                    gained = df[(df['comparison'] == comp) & (df['mod_type'] == mod) &
                               (df['category'] == 'gained') & (df['motif'] == motif)]
                    lost = df[(df['comparison'] == comp) & (df['mod_type'] == mod) &
                             (df['category'] == 'lost') & (df['motif'] == motif)]
                    if len(gained) == 0 or len(lost) == 0:
                        continue
                    g_count = int(gained['count'].values[0])
                    g_n = int(gained['n_sites'].values[0])
                    l_count = int(lost['count'].values[0])
                    l_n = int(lost['n_sites'].values[0])

                    table = [[g_count, g_n - g_count],
                             [l_count, l_n - l_count]]
                    try:
                        odds, p = fisher_exact(table)
                        results.append({
                            'comparison': comp, 'mod_type': mod, 'motif': motif,
                            'gained_pct': g_count / g_n * 100 if g_n > 0 else 0,
                            'lost_pct': l_count / l_n * 100 if l_n > 0 else 0,
                            'gained_n': g_n, 'lost_n': l_n,
                            'odds_ratio': odds, 'fisher_p': p
                        })
                    except:
                        pass

        fisher_df = pd.DataFrame(results)
    else:
        fisher_df = df[df['fisher_p'].notna()].copy()

    if len(fisher_df) == 0:
        print("No Fisher test results to correct")
        return None

    print(f"Total Fisher tests: {len(fisher_df)}")

    reject, pvals_corrected, _, _ = multipletests(fisher_df['fisher_p'], method='fdr_bh')
    fisher_df['p_adj'] = pvals_corrected
    fisher_df['significant_fdr'] = reject

    sig = fisher_df[fisher_df['significant_fdr']].sort_values('p_adj')
    print(f"Significant after FDR: {len(sig)}")
    if len(sig) > 0:
        cols = ['comparison', 'mod_type', 'motif', 'gained_pct', 'lost_pct', 'fisher_p', 'p_adj']
        available_cols = [c for c in cols if c in sig.columns]
        print(sig[available_cols].to_string(index=False))

    fisher_df.to_csv(TSS_DIR / 'temporal_gained_lost_motifs_fdr.csv', index=False)
    return fisher_df


def correct_promoter_genebody():
    """FDR for C2 promoter vs gene body correlations."""
    print("\n" + "=" * 60)
    print("4. C2 Promoter vs Gene body correlations")
    print("=" * 60)

    df = pd.read_csv(TSS_DIR / 'C2_promoter_vs_genebody_correlation.csv')
    print(f"Total tests: {len(df)}")

    reject, pvals_corrected, _, _ = multipletests(df['p_value'], method='fdr_bh')
    df['p_adj'] = pvals_corrected
    df['significant_fdr'] = reject

    sig = df[df['significant_fdr']]
    print(f"Significant after FDR: {len(sig)}")
    if len(sig) > 0:
        print(sig[['mod_type', 'comparison', 'region', 'spearman_r', 'p_value', 'p_adj']].to_string(index=False))

    df.to_csv(TSS_DIR / 'C2_promoter_genebody_fdr.csv', index=False)
    return df


def combined_correction():
    """All tests combined for overall FDR."""
    print("\n" + "=" * 60)
    print("5. COMBINED FDR (all tests pooled)")
    print("=" * 60)

    all_pvals = []
    all_labels = []

    # C1
    c1 = pd.read_csv(TSS_DIR / 'C1_distance_stratified_correlation.csv')
    for _, row in c1.iterrows():
        all_pvals.append(row['p_value'])
        all_labels.append(f"C1:{row['mod_type']}:{row['comparison']}:{row['distance_bin']}")

    # C2
    c2 = pd.read_csv(TSS_DIR / 'C2_promoter_vs_genebody_correlation.csv')
    for _, row in c2.iterrows():
        all_pvals.append(row['p_value'])
        all_labels.append(f"C2:{row['mod_type']}:{row['comparison']}:{row['region']}")

    # SARP zone
    sarp = pd.read_csv(TSS_DIR / 'sarp_zone_correlation.csv')
    for _, row in sarp.iterrows():
        all_pvals.append(row['p_value'])
        all_labels.append(f"SARP:{row['mod_type']}:{row['comparison']}")

    print(f"Total correlation tests: {len(all_pvals)}")

    reject, pvals_corrected, _, _ = multipletests(all_pvals, method='fdr_bh')

    combined = pd.DataFrame({
        'test': all_labels,
        'p_value': all_pvals,
        'p_adj': pvals_corrected,
        'significant': reject
    }).sort_values('p_adj')

    sig = combined[combined['significant']]
    print(f"Significant after combined FDR: {len(sig)}")
    if len(sig) > 0:
        print(sig.head(20).to_string(index=False))
    else:
        print("  No tests survive combined FDR correction")

    # Show top 10
    print(f"\nTop 10 by adjusted p-value:")
    print(combined.head(10).to_string(index=False))

    combined.to_csv(TSS_DIR / 'all_correlations_fdr.csv', index=False)
    return combined


def main():
    c1 = correct_distance_correlations()
    sarp = correct_sarp_zone()
    fisher = correct_gained_lost_motifs()
    c2 = correct_promoter_genebody()
    combined = combined_correction()

    print("\n" + "=" * 60)
    print("FDR CORRECTION SUMMARY")
    print("=" * 60)

    print(f"\nC1 distance correlations: {c1['significant_fdr'].sum()}/{len(c1)} survive FDR")
    if c2 is not None:
        print(f"C2 promoter/genebody: {c2['significant_fdr'].sum()}/{len(c2)} survive FDR")
    if fisher is not None:
        print(f"Gained/Lost Fisher tests: {fisher['significant_fdr'].sum()}/{len(fisher)} survive FDR")
    print(f"Combined (all correlations): {combined['significant'].sum()}/{len(combined)} survive FDR")


if __name__ == '__main__':
    main()
