#!/usr/bin/env python3
"""
Reviewer A5 fix: Promote ST5 (exposed TFs full annotation) to a main-text
table. Output:

  1. tables/supplementary/ST5_exposed_TFs_full_n57.tsv
       - Full 57-row annotation, all columns (Online Resource)
  2. tables/main/Table1_exposed_TFs_main.tsv
       - Curated main-text table: 57 rows × 9 most informative columns
  3. tables/main/Table1_exposed_TFs_main.html
       - Same as TSV, formatted for inline LaTeX inclusion
"""

from pathlib import Path
import importlib
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
_utils = importlib.import_module('00_shared_utils')
for _attr in dir(_utils):
    if not _attr.startswith('_'):
        globals()[_attr] = getattr(_utils, _attr)


SRC = (EPIGENOME / '51_exposed_regulators_characteristics' /
       'tables' / 'exposed_regulators_full_table.tsv')
TEMPORAL = (EPIGENOME / '57_temporal_dynamics_exposed_TF' /
            'tables' / 'temporal_classification.tsv')

OUT_TABLE_SUP = TABLE_SUP_DIR / 'ST5_exposed_TFs_full_n57.tsv'
TABLE_MAIN_DIR = BASE / '15_paper_figures' / 'tables' / 'main'
OUT_TABLE_MAIN = TABLE_MAIN_DIR / 'Table1_exposed_TFs_main.tsv'
OUT_TABLE_HTML = TABLE_MAIN_DIR / 'Table1_exposed_TFs_main.html'


def main():
    print('=== Table 1 / ST5: 57 exposed TFs ===')
    df = pd.read_csv(SRC, sep='\t')
    print(f'  Loaded: {len(df)} TFs')
    assert len(df) == 57, f'Expected 57 exposed TFs, got {len(df)}'

    # Try to merge with temporal classification (bloc: act/rep)
    if TEMPORAL.exists():
        try:
            tc = pd.read_csv(TEMPORAL, sep='\t')
            tc_cols = [c for c in ['locus_tag', 'bloc', 'temporal_class']
                       if c in tc.columns]
            df = df.merge(tc[tc_cols], on='locus_tag', how='left')
            print(f'  Merged temporal_classification: '
                  f'bloc fields = {tc_cols[1:]}')
        except Exception as exc:
            print(f'  Skip temporal merge: {exc}')

    # Save full ST5 (Online Resource)
    OUT_TABLE_SUP.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_TABLE_SUP, sep='\t', index=False)
    print(f'  ST5 (full): {OUT_TABLE_SUP}  ({len(df)} rows)')

    # Build curated main-text Table 1
    main_cols = ['locus_tag', 'old_locus_tag', 'gene_name', 'tf_family',
                 'region', 'nearest_methyl_dist', 'baseMean',
                 'log2FC_T2', 'log2FC_T3']
    if 'bloc' in df.columns:
        main_cols.insert(4, 'bloc')

    sub = df[main_cols].copy()

    # Round numeric columns
    for c in ('baseMean',):
        if c in sub.columns:
            sub[c] = sub[c].round(1)
    for c in ('log2FC_T2', 'log2FC_T3'):
        if c in sub.columns:
            sub[c] = sub[c].round(2)

    # Sort: arm first, then by descending baseMean
    region_rank = {'arm': 0, 'core': 1}
    sub['__r'] = sub['region'].map(region_rank).fillna(2)
    sub = sub.sort_values(['__r', 'baseMean'],
                          ascending=[True, False]).drop(columns='__r')

    OUT_TABLE_MAIN.parent.mkdir(parents=True, exist_ok=True)
    sub.to_csv(OUT_TABLE_MAIN, sep='\t', index=False)
    print(f'  Table 1 (main): {OUT_TABLE_MAIN}  ({len(sub)} rows × '
          f'{len(sub.columns)} cols)')

    # HTML version (LaTeX-friendly inclusion)
    rename = {
        'locus_tag': 'Locus tag',
        'old_locus_tag': 'SCO',
        'gene_name': 'Gene',
        'tf_family': 'TF family',
        'bloc': 'Bloc',
        'region': 'Region',
        'nearest_methyl_dist': 'Nearest methyl. (bp)',
        'baseMean': 'baseMean',
        'log2FC_T2': 'LFC T2 vs T1',
        'log2FC_T3': 'LFC T3 vs T1',
    }
    sub_h = sub.rename(columns=rename)
    sub_h.to_html(OUT_TABLE_HTML, index=False, na_rep='—',
                  border=0, classes='exposed-tf-main')
    print(f'  Table 1 (HTML): {OUT_TABLE_HTML}')

    # Summary breakdown
    print('\n  Summary:')
    if 'region' in df.columns:
        print(f'    region: {df["region"].value_counts().to_dict()}')
    if 'bloc' in df.columns:
        print(f'    bloc:   {df["bloc"].value_counts(dropna=False).to_dict()}')
    if 'tf_family' in df.columns:
        top_fam = df['tf_family'].value_counts().head(5).to_dict()
        print(f'    top TF families: {top_fam}')
    print('=== Done ===')


if __name__ == '__main__':
    main()
