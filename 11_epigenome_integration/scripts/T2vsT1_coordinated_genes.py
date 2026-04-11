#!/usr/bin/env python3
"""
T2 vs T1 Coordinated Methylation-Expression Analysis (Original method)
Focus on genes with both methylation and expression changes
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
METHYL_PATH = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/01_integration/high_confidence_sites_weighted.csv"
DESEQ2_PATH = "/Users/okaban/bioinfo/rna-seq/04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_2_vs_1.tsv"
OUTPUT_DIR = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/01_integration")

PROMOTER_UPSTREAM = 300
PROMOTER_DOWNSTREAM = 50

def parse_gff(gff_path):
    """Parse GFF to get gene info including product annotations."""
    genes = {}
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            if parts[2] == 'gene':
                locus_match = re.search(r'locus_tag=([^;]+)', parts[8])
                old_match = re.search(r'old_locus_tag=([^;]+)', parts[8])
                if locus_match:
                    gene_id = locus_match.group(1)
                    genes[gene_id] = {
                        'chrom': parts[0],
                        'start': int(parts[3]),
                        'end': int(parts[4]),
                        'strand': parts[6],
                        'tss': int(parts[3]) if parts[6] == '+' else int(parts[4]),
                        'old_locus_tag': old_match.group(1) if old_match else '',
                        'product': ''
                    }

            elif parts[2] == 'CDS':
                locus_match = re.search(r'locus_tag=([^;]+)', parts[8])
                product_match = re.search(r'product=([^;]+)', parts[8])
                if locus_match and product_match:
                    gene_id = locus_match.group(1)
                    if gene_id in genes:
                        genes[gene_id]['product'] = product_match.group(1)

    return genes

def main():
    print("=" * 80)
    print("T2 vs T1 Coordinated Methylation-Expression Analysis")
    print("Method: Original (3-rep, equal weight)")
    print("=" * 80)

    # Load gene annotations
    print("\n[1] Loading gene annotations...")
    genes = parse_gff(GFF_PATH)
    print(f"  Loaded {len(genes)} genes")

    # Load methylation data (T1 and T2 only)
    print("\n[2] Loading methylation data...")
    methyl_df = pd.read_csv(METHYL_PATH)
    methyl_t1 = methyl_df[methyl_df['timepoint'] == 'T1']
    methyl_t2 = methyl_df[methyl_df['timepoint'] == 'T2']
    print(f"  T1 sites: {len(methyl_t1)}")
    print(f"  T2 sites: {len(methyl_t2)}")

    # Load DESeq2 results
    print("\n[3] Loading DESeq2 results (T2 vs T1)...")
    deseq2_df = pd.read_csv(DESEQ2_PATH, sep='\t')
    print(f"  Total genes: {len(deseq2_df)}")

    # Assign methylation to promoters
    print("\n[4] Assigning methylation sites to promoter regions...")

    gene_methyl = {}
    for gene_id, gene_info in genes.items():
        tss = gene_info['tss']
        strand = gene_info['strand']
        chrom = gene_info['chrom']

        if strand == '+':
            prom_start, prom_end = tss - PROMOTER_UPSTREAM, tss + PROMOTER_DOWNSTREAM
        else:
            prom_start, prom_end = tss - PROMOTER_DOWNSTREAM, tss + PROMOTER_UPSTREAM

        # T1 methylation
        t1_sites = methyl_t1[(methyl_t1['chrom'] == chrom) &
                            (methyl_t1['position'] >= prom_start) &
                            (methyl_t1['position'] <= prom_end)]

        # T2 methylation
        t2_sites = methyl_t2[(methyl_t2['chrom'] == chrom) &
                            (methyl_t2['position'] >= prom_start) &
                            (methyl_t2['position'] <= prom_end)]

        gene_methyl[gene_id] = {
            '6mA_T1': t1_sites[t1_sites['mod_type'] == '6mA']['weighted_mod_freq'].tolist(),
            '6mA_T2': t2_sites[t2_sites['mod_type'] == '6mA']['weighted_mod_freq'].tolist(),
            '4mC_T1': t1_sites[t1_sites['mod_type'] == '4mC']['weighted_mod_freq'].tolist(),
            '4mC_T2': t2_sites[t2_sites['mod_type'] == '4mC']['weighted_mod_freq'].tolist(),
        }

    # Find genes with coordinated changes
    print("\n[5] Identifying genes with coordinated methylation-expression changes...")

    coordinated_genes = []

    for gene_id, gene_info in genes.items():
        # Get expression data
        expr_data = deseq2_df[deseq2_df['gene_id'] == gene_id]
        if len(expr_data) == 0:
            continue

        log2fc = expr_data.iloc[0]['log2FoldChange']
        padj = expr_data.iloc[0]['padj']
        basemean = expr_data.iloc[0]['baseMean']

        if pd.isna(padj):
            continue

        # Get methylation data
        methyl = gene_methyl[gene_id]

        for mod_type in ['6mA', '4mC']:
            t1_freqs = methyl[f'{mod_type}_T1']
            t2_freqs = methyl[f'{mod_type}_T2']

            # Check if there's methylation in at least one timepoint
            if len(t1_freqs) == 0 and len(t2_freqs) == 0:
                continue

            # Calculate methylation change
            t1_mean = np.mean(t1_freqs) if t1_freqs else 0
            t2_mean = np.mean(t2_freqs) if t2_freqs else 0
            methyl_change = t2_mean - t1_mean

            # Count sites
            t1_count = len(t1_freqs)
            t2_count = len(t2_freqs)

            # Determine methylation change category
            if t1_count == 0 and t2_count > 0:
                methyl_category = "Gained"
            elif t1_count > 0 and t2_count == 0:
                methyl_category = "Lost"
            elif methyl_change > 10:
                methyl_category = "Increased"
            elif methyl_change < -10:
                methyl_category = "Decreased"
            else:
                methyl_category = "Stable"

            # Determine expression change category
            if padj < 0.05:
                if log2fc > 1:
                    expr_category = "Up"
                elif log2fc < -1:
                    expr_category = "Down"
                else:
                    expr_category = "Mild"
            else:
                expr_category = "NS"

            # Record if there's any methylation AND significant expression change
            if (t1_count > 0 or t2_count > 0) and padj < 0.05 and abs(log2fc) > 0.5:
                coordinated_genes.append({
                    'gene_id': gene_id,
                    'old_locus_tag': gene_info['old_locus_tag'],
                    'product': gene_info['product'],
                    'mod_type': mod_type,
                    'T1_sites': t1_count,
                    'T2_sites': t2_count,
                    'T1_mean_freq': round(t1_mean, 1),
                    'T2_mean_freq': round(t2_mean, 1),
                    'methyl_change': round(methyl_change, 1),
                    'methyl_category': methyl_category,
                    'log2FC': round(log2fc, 2),
                    'padj': padj,
                    'baseMean': round(basemean, 1),
                    'expr_category': expr_category,
                    'coordination': f"{methyl_category}_{expr_category}"
                })

    result_df = pd.DataFrame(coordinated_genes)

    # Summary
    print("\n" + "=" * 80)
    print("RESULTS SUMMARY")
    print("=" * 80)

    print(f"\nTotal gene-modification pairs with coordinated changes: {len(result_df)}")

    if len(result_df) > 0:
        print("\n--- By Modification Type ---")
        for mod in ['6mA', '4mC']:
            mod_df = result_df[result_df['mod_type'] == mod]
            print(f"\n{mod}: {len(mod_df)} genes")
            if len(mod_df) > 0:
                print(f"  Expression Up: {len(mod_df[mod_df['expr_category'] == 'Up'])}")
                print(f"  Expression Down: {len(mod_df[mod_df['expr_category'] == 'Down'])}")

        print("\n--- By Coordination Pattern ---")
        coord_counts = result_df['coordination'].value_counts()
        for pattern, count in coord_counts.items():
            print(f"  {pattern}: {count}")

        # Top genes by absolute log2FC
        print("\n" + "=" * 80)
        print("TOP COORDINATED GENES (by |log2FC|)")
        print("=" * 80)

        top_df = result_df.sort_values('log2FC', key=abs, ascending=False).head(30)

        print("\n{:<12} {:<10} {:<6} {:>5} {:>5} {:>8} {:>8} {:>8}".format(
            "Gene", "Old_Locus", "Mod", "T1#", "T2#", "ΔMethyl", "log2FC", "Category"))
        print("-" * 80)

        for _, row in top_df.iterrows():
            print("{:<12} {:<10} {:<6} {:>5} {:>5} {:>8.1f} {:>8.2f} {:>8}".format(
                row['gene_id'][:12],
                str(row['old_locus_tag'])[:10] if row['old_locus_tag'] else '-',
                row['mod_type'],
                row['T1_sites'],
                row['T2_sites'],
                row['methyl_change'],
                row['log2FC'],
                row['coordination']))

        # Positive correlation examples (methyl up + expr up, or methyl down + expr down)
        print("\n" + "=" * 80)
        print("POSITIVE CORRELATION EXAMPLES (Methylation ↔ Expression same direction)")
        print("=" * 80)

        positive_corr = result_df[
            ((result_df['methyl_category'].isin(['Gained', 'Increased'])) & (result_df['expr_category'] == 'Up')) |
            ((result_df['methyl_category'].isin(['Lost', 'Decreased'])) & (result_df['expr_category'] == 'Down'))
        ].sort_values('log2FC', key=abs, ascending=False)

        print(f"\nFound {len(positive_corr)} genes with positive correlation")

        if len(positive_corr) > 0:
            print("\n{:<12} {:<8} {:<6} {:>6} {:>6} {:>7} {:>7} {:<30}".format(
                "Gene", "OldLocus", "Mod", "T1→T2", "ΔMeth", "log2FC", "padj", "Product"))
            print("-" * 100)

            for _, row in positive_corr.head(20).iterrows():
                sites_change = f"{row['T1_sites']}→{row['T2_sites']}"
                product = str(row['product'])[:30] if row['product'] else '-'
                print("{:<12} {:<8} {:<6} {:>6} {:>7.1f} {:>7.2f} {:>7.1e} {:<30}".format(
                    row['gene_id'][:12],
                    str(row['old_locus_tag'])[:8] if row['old_locus_tag'] else '-',
                    row['mod_type'],
                    sites_change,
                    row['methyl_change'],
                    row['log2FC'],
                    row['padj'],
                    product))

        # Negative correlation examples
        print("\n" + "=" * 80)
        print("NEGATIVE CORRELATION EXAMPLES (Methylation ↔ Expression opposite direction)")
        print("=" * 80)

        negative_corr = result_df[
            ((result_df['methyl_category'].isin(['Gained', 'Increased'])) & (result_df['expr_category'] == 'Down')) |
            ((result_df['methyl_category'].isin(['Lost', 'Decreased'])) & (result_df['expr_category'] == 'Up'))
        ].sort_values('log2FC', key=abs, ascending=False)

        print(f"\nFound {len(negative_corr)} genes with negative correlation")

        if len(negative_corr) > 0:
            print("\n{:<12} {:<8} {:<6} {:>6} {:>6} {:>7} {:>7} {:<30}".format(
                "Gene", "OldLocus", "Mod", "T1→T2", "ΔMeth", "log2FC", "padj", "Product"))
            print("-" * 100)

            for _, row in negative_corr.head(20).iterrows():
                sites_change = f"{row['T1_sites']}→{row['T2_sites']}"
                product = str(row['product'])[:30] if row['product'] else '-'
                print("{:<12} {:<8} {:<6} {:>6} {:>7.1f} {:>7.2f} {:>7.1e} {:<30}".format(
                    row['gene_id'][:12],
                    str(row['old_locus_tag'])[:8] if row['old_locus_tag'] else '-',
                    row['mod_type'],
                    sites_change,
                    row['methyl_change'],
                    row['log2FC'],
                    row['padj'],
                    product))

    # Save results
    print("\n[6] Saving results...")
    result_df.to_csv(OUTPUT_DIR / 'T2vsT1_coordinated_genes.csv', index=False)

    if len(positive_corr) > 0:
        positive_corr.to_csv(OUTPUT_DIR / 'T2vsT1_positive_correlation.csv', index=False)
    if len(negative_corr) > 0:
        negative_corr.to_csv(OUTPUT_DIR / 'T2vsT1_negative_correlation.csv', index=False)

    print(f"\nResults saved to: {OUTPUT_DIR}")

    return result_df, positive_corr, negative_corr

if __name__ == "__main__":
    result_df, positive_corr, negative_corr = main()
