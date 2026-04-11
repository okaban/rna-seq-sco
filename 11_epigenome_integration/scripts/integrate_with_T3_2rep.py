#!/usr/bin/env python3
"""
Re-run integration analysis with T3 2-replicate data
Compare results with original 3-replicate analysis
"""

import pandas as pd
import numpy as np
from pathlib import Path
import re
from collections import defaultdict
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Configuration
PROMOTER_UPSTREAM = 300
PROMOTER_DOWNSTREAM = 50

# Paths
GFF_PATH = "/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/genomic.gff"
METHYL_ORIGINAL = "/Users/okaban/bioinfo/methyl/260102_M145/analysis/cursor_results/20260108/high_confidence_sites.csv"
METHYL_NEW = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/high_confidence_sites_T3_2rep.csv"
DESEQ2_DIR = "/Users/okaban/bioinfo/rna-seq/04_deseq2/analysis/04_deseq2_260128_v1/results"
OUTPUT_DIR = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis"

def parse_gff(gff_path):
    genes = []
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'gene':
                continue
            locus_match = re.search(r'locus_tag=([^;]+)', parts[8])
            if locus_match:
                genes.append({
                    'gene_id': locus_match.group(1),
                    'chrom': parts[0],
                    'start': int(parts[3]),
                    'end': int(parts[4]),
                    'strand': parts[6],
                    'tss': int(parts[3]) if parts[6] == '+' else int(parts[4])
                })
    return pd.DataFrame(genes)

def load_deseq2_results(deseq2_dir):
    results = {}
    for file_suffix, label in [('M145_2_vs_1', 'T2_vs_T1'), ('M145_3_vs_1', 'T3_vs_T1'), ('M145_3_vs_2', 'T3_vs_T2')]:
        path = Path(deseq2_dir) / f"DESeq2_{file_suffix}.tsv"
        if path.exists():
            results[label] = pd.read_csv(path, sep='\t')
    return results

def assign_methyl_to_promoters(genes_df, methyl_df, upstream=300, downstream=50):
    promoter_methyl = defaultdict(list)
    for _, gene in genes_df.iterrows():
        tss = gene['tss']
        strand = gene['strand']
        if strand == '+':
            prom_start, prom_end = tss - upstream, tss + downstream
        else:
            prom_start, prom_end = tss - downstream, tss + upstream
        mask = (methyl_df['chrom'] == gene['chrom']) & (methyl_df['position'] >= prom_start) & (methyl_df['position'] <= prom_end)
        for _, site in methyl_df[mask].iterrows():
            promoter_methyl[gene['gene_id']].append({
                'mod_type': site['mod_type'],
                'timepoint': site['timepoint'],
                'mean_mod_freq': site['mean_mod_freq']
            })
    return promoter_methyl

def create_integrated_table(genes_df, promoter_methyl, deseq2_results):
    records = []
    for gene_id in genes_df['gene_id'].unique():
        record = {'gene_id': gene_id}
        for comparison, df in deseq2_results.items():
            gene_data = df[df['gene_id'] == gene_id]
            if len(gene_data) > 0:
                record[f'log2FC_{comparison}'] = gene_data.iloc[0]['log2FoldChange']
                record[f'padj_{comparison}'] = gene_data.iloc[0]['padj']
        for mod_type in ['6mA', '4mC']:
            for tp in ['T1', 'T2', 'T3']:
                sites = [s for s in promoter_methyl.get(gene_id, []) if s['mod_type'] == mod_type and s['timepoint'] == tp]
                record[f'{mod_type}_{tp}_count'] = len(sites)
                record[f'{mod_type}_{tp}_mean_freq'] = np.mean([s['mean_mod_freq'] for s in sites]) if sites else 0
        records.append(record)
    return pd.DataFrame(records)

def calculate_methyl_change(df):
    for mod in ['6mA', '4mC']:
        df[f'{mod}_change_T2_vs_T1'] = df[f'{mod}_T2_mean_freq'] - df[f'{mod}_T1_mean_freq']
        df[f'{mod}_change_T3_vs_T1'] = df[f'{mod}_T3_mean_freq'] - df[f'{mod}_T1_mean_freq']
        df[f'{mod}_change_T3_vs_T2'] = df[f'{mod}_T3_mean_freq'] - df[f'{mod}_T2_mean_freq']
    return df

def correlation_analysis(df):
    results = []
    for mod in ['6mA', '4mC']:
        for comp in ['T2_vs_T1', 'T3_vs_T1', 'T3_vs_T2']:
            methyl_col, expr_col = f'{mod}_change_{comp}', f'log2FC_{comp}'
            if methyl_col in df.columns and expr_col in df.columns:
                mask = df[methyl_col].abs() > 0
                subset = df[mask].dropna(subset=[methyl_col, expr_col])
                if len(subset) >= 10:
                    corr, pval = stats.spearmanr(subset[methyl_col], subset[expr_col])
                    results.append({'mod_type': mod, 'comparison': comp, 'n_genes': len(subset), 'spearman_r': corr, 'p_value': pval})
    return pd.DataFrame(results)

def main():
    print("=" * 70)
    print("Comparison: Original (3 rep) vs T3-2rep Integration Analysis")
    print("=" * 70)

    genes_df = parse_gff(GFF_PATH)
    deseq2_results = load_deseq2_results(DESEQ2_DIR)

    # Original analysis
    print("\n[1] Original analysis (T3 with 3 replicates)...")
    methyl_orig = pd.read_csv(METHYL_ORIGINAL)
    prom_orig = assign_methyl_to_promoters(genes_df, methyl_orig, PROMOTER_UPSTREAM, PROMOTER_DOWNSTREAM)
    int_orig = create_integrated_table(genes_df, prom_orig, deseq2_results)
    int_orig = calculate_methyl_change(int_orig)
    corr_orig = correlation_analysis(int_orig)

    # New analysis with T3 2-rep
    print("[2] New analysis (T3 with 2 replicates, excluding 3-2)...")
    methyl_new = pd.read_csv(METHYL_NEW)
    prom_new = assign_methyl_to_promoters(genes_df, methyl_new, PROMOTER_UPSTREAM, PROMOTER_DOWNSTREAM)
    int_new = create_integrated_table(genes_df, prom_new, deseq2_results)
    int_new = calculate_methyl_change(int_new)
    corr_new = correlation_analysis(int_new)

    # Comparison
    print("\n" + "=" * 70)
    print("PROMOTER METHYLATION COMPARISON")
    print("=" * 70)

    for mod in ['6mA', '4mC']:
        print(f"\n{mod}:")
        for tp in ['T1', 'T2', 'T3']:
            orig_count = (int_orig[f'{mod}_{tp}_count'] > 0).sum()
            new_count = (int_new[f'{mod}_{tp}_count'] > 0).sum()
            diff = new_count - orig_count
            print(f"  {tp}: {orig_count} -> {new_count} genes ({diff:+d})")

    print("\n" + "=" * 70)
    print("CORRELATION ANALYSIS COMPARISON")
    print("=" * 70)

    print("\n{:<8} {:<12} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}".format(
        "ModType", "Comparison", "n(orig)", "r(orig)", "p(orig)", "n(new)", "r(new)", "p(new)"))
    print("-" * 80)

    for mod in ['6mA', '4mC']:
        for comp in ['T2_vs_T1', 'T3_vs_T1', 'T3_vs_T2']:
            orig_row = corr_orig[(corr_orig['mod_type'] == mod) & (corr_orig['comparison'] == comp)]
            new_row = corr_new[(corr_new['mod_type'] == mod) & (corr_new['comparison'] == comp)]

            if len(orig_row) > 0 and len(new_row) > 0:
                print("{:<8} {:<12} {:>10} {:>10.3f} {:>10.4f} {:>10} {:>10.3f} {:>10.4f}".format(
                    mod, comp,
                    orig_row.iloc[0]['n_genes'], orig_row.iloc[0]['spearman_r'], orig_row.iloc[0]['p_value'],
                    new_row.iloc[0]['n_genes'], new_row.iloc[0]['spearman_r'], new_row.iloc[0]['p_value']))

    # Save results
    output_path = Path(OUTPUT_DIR)
    int_new.to_csv(output_path / 'integrated_methyl_expression_T3_2rep.csv', index=False)
    corr_new.to_csv(output_path / 'correlation_analysis_T3_2rep.csv', index=False)

    # Create comparison summary
    print("\n" + "=" * 70)
    print("KEY FINDINGS")
    print("=" * 70)

    # T3 vs T1 comparison improvement
    for mod in ['6mA', '4mC']:
        orig = corr_orig[(corr_orig['mod_type'] == mod) & (corr_orig['comparison'] == 'T3_vs_T1')]
        new = corr_new[(corr_new['mod_type'] == mod) & (corr_new['comparison'] == 'T3_vs_T1')]
        if len(orig) > 0 and len(new) > 0:
            print(f"\n{mod} T3_vs_T1:")
            print(f"  Genes analyzed: {orig.iloc[0]['n_genes']} -> {new.iloc[0]['n_genes']} ({new.iloc[0]['n_genes'] - orig.iloc[0]['n_genes']:+d})")
            print(f"  Correlation: {orig.iloc[0]['spearman_r']:.3f} -> {new.iloc[0]['spearman_r']:.3f}")
            print(f"  P-value: {orig.iloc[0]['p_value']:.4f} -> {new.iloc[0]['p_value']:.4f}")

    print("\n\nResults saved to:", output_path)

if __name__ == "__main__":
    main()
