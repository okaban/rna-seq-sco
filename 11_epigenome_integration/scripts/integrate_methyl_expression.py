#!/usr/bin/env python3
"""
Epigenome-Transcriptome Integration Analysis
Streptomyces coelicolor A3(2) M145

Integrates methylation data (6mA, 4mC) with gene expression changes.
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
PROMOTER_UPSTREAM = 300  # bp upstream of TSS
PROMOTER_DOWNSTREAM = 50  # bp downstream of TSS
THREADS = 10  # Use 10 of 14 cores

# Paths
GFF_PATH = "/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/genomic.gff"
METHYL_PATH = "/Users/okaban/bioinfo/methyl/260102_M145/analysis/cursor_results/20260108/high_confidence_sites.csv"
DESEQ2_DIR = "/Users/okaban/bioinfo/rna-seq/04_deseq2/analysis/04_deseq2_260128_v1/results"
OUTPUT_DIR = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis"

def parse_gff(gff_path):
    """Parse GFF file to extract gene coordinates."""
    genes = []
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            if parts[2] == 'gene':
                chrom = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attrs = parts[8]
                # Extract locus_tag
                locus_match = re.search(r'locus_tag=([^;]+)', attrs)
                if locus_match:
                    locus_tag = locus_match.group(1)
                    # Get old_locus_tag if available
                    old_match = re.search(r'old_locus_tag=([^;]+)', attrs)
                    old_locus = old_match.group(1) if old_match else None
                    genes.append({
                        'gene_id': locus_tag,
                        'old_locus_tag': old_locus,
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'tss': start if strand == '+' else end
                    })
    return pd.DataFrame(genes)

def load_methylation_data(methyl_path):
    """Load high-confidence methylation sites."""
    df = pd.read_csv(methyl_path)
    print(f"Loaded {len(df)} methylation sites")
    print(f"Modification types: {df['mod_type'].value_counts().to_dict()}")
    print(f"Timepoints: {df['timepoint'].value_counts().to_dict()}")
    return df

def load_deseq2_results(deseq2_dir):
    """Load DESeq2 differential expression results."""
    results = {}
    comparisons = [
        ('M145_2_vs_1', 'T2_vs_T1'),
        ('M145_3_vs_1', 'T3_vs_T1'),
        ('M145_3_vs_2', 'T3_vs_T2')
    ]
    for file_suffix, label in comparisons:
        path = Path(deseq2_dir) / f"DESeq2_{file_suffix}.tsv"
        if path.exists():
            df = pd.read_csv(path, sep='\t')
            results[label] = df
            print(f"Loaded {label}: {len(df)} genes")
    return results

def assign_methyl_to_promoters(genes_df, methyl_df, upstream=300, downstream=50):
    """Assign methylation sites to gene promoter regions."""
    promoter_methyl = defaultdict(list)

    for _, gene in genes_df.iterrows():
        gene_id = gene['gene_id']
        tss = gene['tss']
        strand = gene['strand']
        chrom = gene['chrom']

        if strand == '+':
            prom_start = tss - upstream
            prom_end = tss + downstream
        else:
            prom_start = tss - downstream
            prom_end = tss + upstream

        # Find methylation sites in promoter region
        mask = (
            (methyl_df['chrom'] == chrom) &
            (methyl_df['position'] >= prom_start) &
            (methyl_df['position'] <= prom_end)
        )
        sites = methyl_df[mask]

        if len(sites) > 0:
            for _, site in sites.iterrows():
                promoter_methyl[gene_id].append({
                    'position': site['position'],
                    'mod_type': site['mod_type'],
                    'timepoint': site['timepoint'],
                    'mean_mod_freq': site['mean_mod_freq'],
                    'distance_to_tss': site['position'] - tss if strand == '+' else tss - site['position']
                })

    return promoter_methyl

def create_integrated_table(genes_df, promoter_methyl, deseq2_results):
    """Create integrated table with methylation and expression data."""
    records = []

    for gene_id in genes_df['gene_id'].unique():
        gene_info = genes_df[genes_df['gene_id'] == gene_id].iloc[0]

        record = {
            'gene_id': gene_id,
            'old_locus_tag': gene_info['old_locus_tag'],
            'chrom': gene_info['chrom'],
            'start': gene_info['start'],
            'end': gene_info['end'],
            'strand': gene_info['strand']
        }

        # Add DESeq2 results
        for comparison, df in deseq2_results.items():
            gene_data = df[df['gene_id'] == gene_id]
            if len(gene_data) > 0:
                record[f'log2FC_{comparison}'] = gene_data.iloc[0]['log2FoldChange']
                record[f'padj_{comparison}'] = gene_data.iloc[0]['padj']
                record[f'baseMean_{comparison}'] = gene_data.iloc[0]['baseMean']

        # Add methylation data
        if gene_id in promoter_methyl:
            methyl_sites = promoter_methyl[gene_id]

            for mod_type in ['6mA', '4mC', '5mC']:
                for tp in ['T1', 'T2', 'T3']:
                    sites = [s for s in methyl_sites if s['mod_type'] == mod_type and s['timepoint'] == tp]
                    record[f'{mod_type}_{tp}_count'] = len(sites)
                    if sites:
                        record[f'{mod_type}_{tp}_mean_freq'] = np.mean([s['mean_mod_freq'] for s in sites])
                    else:
                        record[f'{mod_type}_{tp}_mean_freq'] = 0
        else:
            for mod_type in ['6mA', '4mC', '5mC']:
                for tp in ['T1', 'T2', 'T3']:
                    record[f'{mod_type}_{tp}_count'] = 0
                    record[f'{mod_type}_{tp}_mean_freq'] = 0

        records.append(record)

    return pd.DataFrame(records)

def calculate_methyl_change(integrated_df):
    """Calculate methylation changes between timepoints."""
    for mod_type in ['6mA', '4mC']:
        # T2 vs T1
        integrated_df[f'{mod_type}_change_T2_vs_T1'] = (
            integrated_df[f'{mod_type}_T2_mean_freq'] - integrated_df[f'{mod_type}_T1_mean_freq']
        )
        # T3 vs T1
        integrated_df[f'{mod_type}_change_T3_vs_T1'] = (
            integrated_df[f'{mod_type}_T3_mean_freq'] - integrated_df[f'{mod_type}_T1_mean_freq']
        )
        # T3 vs T2
        integrated_df[f'{mod_type}_change_T3_vs_T2'] = (
            integrated_df[f'{mod_type}_T3_mean_freq'] - integrated_df[f'{mod_type}_T2_mean_freq']
        )
    return integrated_df

def correlation_analysis(integrated_df):
    """Perform correlation analysis between methylation and expression."""
    results = []

    for mod_type in ['6mA', '4mC']:
        for comparison in ['T2_vs_T1', 'T3_vs_T1', 'T3_vs_T2']:
            methyl_col = f'{mod_type}_change_{comparison}'
            expr_col = f'log2FC_{comparison}'

            if methyl_col in integrated_df.columns and expr_col in integrated_df.columns:
                # Filter genes with methylation changes
                mask = integrated_df[methyl_col].abs() > 0
                subset = integrated_df[mask].dropna(subset=[methyl_col, expr_col])

                if len(subset) >= 10:
                    corr, pval = stats.spearmanr(subset[methyl_col], subset[expr_col])
                    results.append({
                        'mod_type': mod_type,
                        'comparison': comparison,
                        'n_genes': len(subset),
                        'spearman_r': corr,
                        'p_value': pval
                    })

    return pd.DataFrame(results)

def identify_genes_with_methyl_expression_changes(integrated_df, log2fc_threshold=1.0, padj_threshold=0.05):
    """Identify genes with both methylation and expression changes."""
    significant_genes = []

    for comparison in ['T2_vs_T1', 'T3_vs_T1', 'T3_vs_T2']:
        expr_col = f'log2FC_{comparison}'
        padj_col = f'padj_{comparison}'

        if expr_col not in integrated_df.columns:
            continue

        for mod_type in ['6mA', '4mC']:
            methyl_col = f'{mod_type}_change_{comparison}'
            if methyl_col not in integrated_df.columns:
                continue

            # Significant DE genes with methylation changes
            mask = (
                (integrated_df[padj_col] < padj_threshold) &
                (integrated_df[expr_col].abs() > log2fc_threshold) &
                (integrated_df[methyl_col].abs() > 10)  # >10% methylation change
            )

            subset = integrated_df[mask].copy()
            subset['comparison'] = comparison
            subset['mod_type'] = mod_type
            significant_genes.append(subset[['gene_id', 'old_locus_tag', expr_col, padj_col, methyl_col, 'comparison', 'mod_type']])

    if significant_genes:
        return pd.concat(significant_genes, ignore_index=True)
    return pd.DataFrame()

def generate_summary_stats(integrated_df, methyl_df):
    """Generate summary statistics."""
    summary = {
        'Total genes': len(integrated_df),
        'Genes with promoter methylation': (integrated_df[[c for c in integrated_df.columns if '_count' in c]].sum(axis=1) > 0).sum(),
    }

    for mod_type in ['6mA', '4mC']:
        for tp in ['T1', 'T2', 'T3']:
            col = f'{mod_type}_{tp}_count'
            if col in integrated_df.columns:
                summary[f'Genes with {mod_type} in {tp}'] = (integrated_df[col] > 0).sum()

    return summary

def main():
    print("=" * 60)
    print("Epigenome-Transcriptome Integration Analysis")
    print("Streptomyces coelicolor A3(2) M145")
    print("=" * 60)

    # Create output directory
    output_path = Path(OUTPUT_DIR)
    output_path.mkdir(parents=True, exist_ok=True)

    # Load data
    print("\n[1] Loading gene annotations...")
    genes_df = parse_gff(GFF_PATH)
    print(f"Loaded {len(genes_df)} genes")

    print("\n[2] Loading methylation data...")
    methyl_df = load_methylation_data(METHYL_PATH)

    print("\n[3] Loading DESeq2 results...")
    deseq2_results = load_deseq2_results(DESEQ2_DIR)

    print("\n[4] Assigning methylation sites to promoter regions...")
    print(f"Promoter region: -{PROMOTER_UPSTREAM} to +{PROMOTER_DOWNSTREAM} bp from TSS")
    promoter_methyl = assign_methyl_to_promoters(genes_df, methyl_df, PROMOTER_UPSTREAM, PROMOTER_DOWNSTREAM)
    print(f"Genes with promoter methylation: {len(promoter_methyl)}")

    print("\n[5] Creating integrated table...")
    integrated_df = create_integrated_table(genes_df, promoter_methyl, deseq2_results)
    integrated_df = calculate_methyl_change(integrated_df)

    print("\n[6] Performing correlation analysis...")
    correlation_df = correlation_analysis(integrated_df)
    print("\nCorrelation Results (Methylation change vs Expression change):")
    print(correlation_df.to_string(index=False))

    print("\n[7] Identifying genes with coordinated changes...")
    significant_df = identify_genes_with_methyl_expression_changes(integrated_df)
    if len(significant_df) > 0:
        print(f"Found {len(significant_df)} gene-modification combinations with significant changes")

    print("\n[8] Generating summary statistics...")
    summary = generate_summary_stats(integrated_df, methyl_df)
    for key, value in summary.items():
        print(f"  {key}: {value}")

    # Save outputs
    print("\n[9] Saving results...")
    integrated_df.to_csv(output_path / 'integrated_methyl_expression.csv', index=False)
    correlation_df.to_csv(output_path / 'correlation_analysis.csv', index=False)
    if len(significant_df) > 0:
        significant_df.to_csv(output_path / 'significant_coordinated_changes.csv', index=False)

    # Save summary
    with open(output_path / 'analysis_summary.txt', 'w') as f:
        f.write("Epigenome-Transcriptome Integration Analysis Summary\n")
        f.write("=" * 60 + "\n\n")
        f.write("Data Overview:\n")
        for key, value in summary.items():
            f.write(f"  {key}: {value}\n")
        f.write("\nCorrelation Analysis:\n")
        f.write(correlation_df.to_string(index=False))
        f.write("\n")

    print(f"\nResults saved to: {output_path}")
    print("Done!")

    return integrated_df, correlation_df, significant_df

if __name__ == "__main__":
    integrated_df, correlation_df, significant_df = main()
