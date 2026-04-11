#!/usr/bin/env python3
"""
Supplementary Tables ST1-ST6
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import importlib
_utils = importlib.import_module('00_shared_utils')
for _attr in dir(_utils):
    if not _attr.startswith('_'):
        globals()[_attr] = getattr(_utils, _attr)

import numpy as np
import pandas as pd

OUT = TABLE_SUP_DIR


def st1_all_methylation_sites():
    """ST1: Complete methylation site catalog."""
    print('ST1: All methylation sites...')
    df_4mc, df_6ma = load_methylation_census()
    df_all = pd.concat([df_4mc, df_6ma], ignore_index=True)

    # Select and order columns
    cols = ['chrom', 'position', 'strand', 'mod_type', 'timepoint',
            'frequency', 'final_motif', 'sequence']
    df_out = df_all[cols].sort_values(['mod_type', 'position', 'timepoint'])

    out_path = OUT / 'ST1_methylation_sites.tsv'
    df_out.to_csv(out_path, sep='\t', index=False)
    print(f'  {len(df_out):,} rows → {out_path.name}')


def st2_motif_summary():
    """ST2: Motif attribution summary with REBASE conservation."""
    print('ST2: Motif summary + REBASE...')
    df_motif = load_motif_summary()
    df_rebase = load_rebase_conservation()

    # Merge REBASE info onto motif summary
    rebase_this = df_rebase[df_rebase['source'] == 'This study'][
        ['motif', 'rate_pct', 'ci_low', 'ci_high', 'n_match', 'n_total']
    ].rename(columns={
        'rate_pct': 'rebase_conservation_pct',
        'ci_low': 'rebase_ci_low',
        'ci_high': 'rebase_ci_high',
        'n_match': 'rebase_n_match',
        'n_total': 'rebase_n_total',
    })

    df_out = df_motif.merge(rebase_this, on='motif', how='left')

    out_path = OUT / 'ST2_motif_summary.tsv'
    df_out.to_csv(out_path, sep='\t', index=False)
    print(f'  {len(df_out)} motifs → {out_path.name}')


def st3_mtase_genes():
    """ST3: MTase gene catalog with expression."""
    print('ST3: MTase genes...')
    df = load_mtase_expression()

    out_path = OUT / 'ST3_mtase_genes.tsv'
    df.to_csv(out_path, sep='\t', index=False)
    print(f'  {len(df)} genes → {out_path.name}')


def st4_bs_methylation_overlap():
    """ST4: TF binding site × methylation site overlaps (Tier 1, Jeong2016 TSS)."""
    print('ST4: BS-methylation overlaps...')

    # Import BS resolution from Figure 4 script
    fig4 = importlib.import_module('04_figure4_tf_bs_depletion')
    bs_resolved = fig4.resolve_bs_jeong2016()

    df_unique = load_methylation_unique_positions()
    meth_pos = df_unique['position'].values
    meth_mod = df_unique['mod_type'].values
    meth_motif = df_unique['final_motif'].values

    overlaps = []
    for _, bs in bs_resolved.iterrows():
        mask = (meth_pos >= bs['abs_start']) & (meth_pos <= bs['abs_end'])
        idx = np.where(mask)[0]
        for i in idx:
            overlaps.append({
                'TF_name': bs['TF_name'],
                'BS_source': bs['BS_source'],
                'BS_abs_start': bs['abs_start'],
                'BS_abs_end': bs['abs_end'],
                'meth_position': meth_pos[i],
                'meth_mod_type': meth_mod[i],
                'meth_motif': meth_motif[i],
            })

    df_out = pd.DataFrame(overlaps)
    out_path = OUT / 'ST4_BS_methylation_overlaps.tsv'
    df_out.to_csv(out_path, sep='\t', index=False)
    print(f'  {len(df_out)} overlaps → {out_path.name}')


def st5_temporal_classification():
    """ST5: Temporal variation classification per unique position."""
    print('ST5: Temporal classification...')
    df_4mc, df_6ma = load_methylation_census()
    df_all = pd.concat([df_4mc, df_6ma], ignore_index=True)

    # For each unique (position, mod_type), determine timepoints present
    grouped = df_all.groupby(['position', 'mod_type']).agg(
        timepoints=('timepoint', lambda x: ','.join(sorted(x.unique()))),
        n_timepoints=('timepoint', 'nunique'),
        mean_frequency=('frequency', 'mean'),
        final_motif=('final_motif', 'first'),
        strand=('strand', 'first'),
    ).reset_index()

    # Classification based on timepoint of detection
    # Note: each HC site is detected at exactly one timepoint in this dataset
    def classify(row):
        tps = set(row['timepoints'].split(','))
        if 'T1' in tps and len(tps) == 1:
            return 'T1-specific'
        elif 'T2' in tps and len(tps) == 1:
            return 'T2-specific'
        elif 'T3' in tps and len(tps) == 1:
            return 'T3-specific'
        elif tps == {'T1', 'T2', 'T3'}:
            return 'All timepoints'
        else:
            return 'Multi-timepoint'

    grouped['temporal_class'] = grouped.apply(classify, axis=1)

    out_path = OUT / 'ST5_temporal_classification.tsv'
    grouped.to_csv(out_path, sep='\t', index=False)
    print(f'  {len(grouped):,} positions → {out_path.name}')
    print(f'  Classes: {grouped["temporal_class"].value_counts().to_dict()}')


def st6_tss_promoter_methylation():
    """ST6: Jeong2016 TSS genes + promoter methylation status."""
    print('ST6: TSS + promoter methylation...')
    df_tss = load_tss_jeong2016()
    df_unique = load_methylation_unique_positions()

    pos_4mc = np.sort(df_unique[df_unique['mod_type'] == '4mC']['position'].values)
    pos_6ma = np.sort(df_unique[df_unique['mod_type'] == '6mA']['position'].values)

    # For each gene, count methylation sites in promoter (-300 to TSS)
    n_4mc_list = []
    n_6ma_list = []
    for _, row in df_tss.iterrows():
        tss = int(row['tss'])
        strand = row['strand']
        if strand == '+':
            prom_start = max(0, tss - 300)
            prom_end = tss
        else:
            prom_start = tss + 1
            prom_end = tss + 301

        # Use searchsorted for fast counting
        n_4mc = (np.searchsorted(pos_4mc, prom_end, side='right') -
                 np.searchsorted(pos_4mc, prom_start, side='left'))
        n_6ma = (np.searchsorted(pos_6ma, prom_end, side='right') -
                 np.searchsorted(pos_6ma, prom_start, side='left'))
        n_4mc_list.append(n_4mc)
        n_6ma_list.append(n_6ma)

    df_tss = df_tss.copy()
    df_tss['promoter_4mC_sites'] = n_4mc_list
    df_tss['promoter_6mA_sites'] = n_6ma_list
    df_tss['promoter_methylated'] = (
        (df_tss['promoter_4mC_sites'] > 0) |
        (df_tss['promoter_6mA_sites'] > 0)
    )

    out_path = OUT / 'ST6_TSS_promoter_methylation.tsv'
    df_tss.to_csv(out_path, sep='\t', index=False)
    n_meth = df_tss['promoter_methylated'].sum()
    print(f'  {len(df_tss)} genes, {n_meth} with promoter methylation → {out_path.name}')


def main():
    print('=== Supplementary Tables ===')
    OUT.mkdir(parents=True, exist_ok=True)

    st1_all_methylation_sites()
    st2_motif_summary()
    st3_mtase_genes()
    st4_bs_methylation_overlap()
    st5_temporal_classification()
    st6_tss_promoter_methylation()

    print('\n=== All supplementary tables complete ===')


if __name__ == '__main__':
    main()
