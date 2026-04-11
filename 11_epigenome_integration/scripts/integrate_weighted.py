#!/usr/bin/env python3
"""
Integration analysis with coverage-weighted methylation data
"""

import pandas as pd
import numpy as np
from pathlib import Path
import re
from collections import defaultdict
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Paths
GFF_PATH = "/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/genomic.gff"
METHYL_WEIGHTED = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/high_confidence_sites_weighted.csv"
DESEQ2_DIR = Path("/Users/okaban/bioinfo/rna-seq/04_deseq2/analysis/04_deseq2_260128_v1/results")
OUTPUT_DIR = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis")

PROMOTER_UPSTREAM = 300
PROMOTER_DOWNSTREAM = 50

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
        path = deseq2_dir / f"DESeq2_{file_suffix}.tsv"
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
                'mean_mod_freq': site['weighted_mod_freq']
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
    print("Integration Analysis with Coverage-Weighted Methylation Data")
    print("=" * 70)

    # Load data
    genes_df = parse_gff(GFF_PATH)
    methyl_df = pd.read_csv(METHYL_WEIGHTED)
    deseq2_results = load_deseq2_results(DESEQ2_DIR)

    print(f"\nLoaded {len(genes_df)} genes")
    print(f"Loaded {len(methyl_df)} weighted methylation sites")

    # Assign to promoters
    promoter_methyl = assign_methyl_to_promoters(genes_df, methyl_df, PROMOTER_UPSTREAM, PROMOTER_DOWNSTREAM)
    print(f"Genes with promoter methylation: {len(promoter_methyl)}")

    # Create integrated table
    int_df = create_integrated_table(genes_df, promoter_methyl, deseq2_results)
    int_df = calculate_methyl_change(int_df)

    # Correlation analysis
    corr_df = correlation_analysis(int_df)

    print("\n" + "=" * 70)
    print("CORRELATION RESULTS (Weighted)")
    print("=" * 70)
    print("\n{:<8} {:<12} {:>10} {:>12} {:>12}".format("Mod", "Comparison", "n_genes", "Spearman r", "p-value"))
    print("-" * 60)
    for _, row in corr_df.iterrows():
        sig = "***" if row['p_value'] < 0.001 else "**" if row['p_value'] < 0.01 else "*" if row['p_value'] < 0.05 else ""
        print("{:<8} {:<12} {:>10} {:>12.3f} {:>12.4f} {}".format(
            row['mod_type'], row['comparison'], row['n_genes'], row['spearman_r'], row['p_value'], sig))

    # Compare all three methods
    print("\n" + "=" * 70)
    print("COMPARISON: All Methods")
    print("=" * 70)

    # Previous results
    original_corr = {
        ('6mA', 'T2_vs_T1'): (361, 0.169, 0.0012),
        ('6mA', 'T3_vs_T1'): (330, -0.050, 0.3673),
        ('6mA', 'T3_vs_T2'): (343, -0.166, 0.0021),
        ('4mC', 'T2_vs_T1'): (474, 0.172, 0.0002),
        ('4mC', 'T3_vs_T1'): (346, -0.041, 0.4446),
        ('4mC', 'T3_vs_T2'): (396, -0.093, 0.0632)
    }

    tworep_corr = {
        ('6mA', 'T2_vs_T1'): (361, 0.169, 0.0012),
        ('6mA', 'T3_vs_T1'): (449, 0.064, 0.1730),
        ('6mA', 'T3_vs_T2'): (458, -0.059, 0.2056),
        ('4mC', 'T2_vs_T1'): (474, 0.172, 0.0002),
        ('4mC', 'T3_vs_T1'): (430, 0.083, 0.0857),
        ('4mC', 'T3_vs_T2'): (453, -0.003, 0.9522)
    }

    print("\n{:<8} {:<12} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}".format(
        "Mod", "Comparison", "r(orig)", "r(2rep)", "r(wgt)", "p(orig)", "p(2rep)", "p(wgt)"))
    print("-" * 80)

    for _, row in corr_df.iterrows():
        mod, comp = row['mod_type'], row['comparison']
        orig = original_corr.get((mod, comp), (0, 0, 1))
        two = tworep_corr.get((mod, comp), (0, 0, 1))
        print("{:<8} {:<12} {:>8.3f} {:>8.3f} {:>8.3f} {:>8.4f} {:>8.4f} {:>8.4f}".format(
            mod, comp, orig[1], two[1], row['spearman_r'], orig[2], two[2], row['p_value']))

    # Promoter summary
    print("\n" + "=" * 70)
    print("PROMOTER METHYLATION SUMMARY (Weighted)")
    print("=" * 70)

    for mod in ['6mA', '4mC']:
        print(f"\n{mod}:")
        for tp in ['T1', 'T2', 'T3']:
            count = (int_df[f'{mod}_{tp}_count'] > 0).sum()
            print(f"  {tp}: {count} genes")

    # Save results
    int_df.to_csv(OUTPUT_DIR / 'integrated_methyl_expression_weighted.csv', index=False)
    corr_df.to_csv(OUTPUT_DIR / 'correlation_analysis_weighted.csv', index=False)

    print(f"\nResults saved to: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()
