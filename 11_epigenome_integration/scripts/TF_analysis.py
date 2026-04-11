#!/usr/bin/env python3
"""
Transcription Factor Analysis in Methylation-Expression Correlation
Including SARP family search
"""

import pandas as pd
from pathlib import Path

DATA_PATH = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis")

# Known SARP family genes in S. coelicolor
SARP_GENES = {
    'SCO5085': 'actII-ORF4 (Act cluster activator)',
    'SCO5877': 'redD (Red cluster activator)',
    'SCO3217': 'cdaR (CDA cluster activator)',
    'SCO6269': 'cpkO (Cpk cluster regulator)',
    'SCO4425': 'afsR (global regulator)',
    'SCO4423': 'afsR2',
    'SCO6992': 'nagE2 (SARP-like)',
}

# TF family keywords
TF_KEYWORDS = [
    'transcription', 'regulator', 'activator', 'repressor',
    'DNA-binding', 'HTH', 'helix-turn-helix', 'sigma',
    'TetR', 'LuxR', 'MarR', 'AraC', 'LysR', 'GntR', 'WhiB',
    'response regulator', 'two-component', 'MerR', 'ArsR',
    'SARP', 'antibiotic regulatory'
]

def load_data():
    coord_df = pd.read_csv(DATA_PATH / 'T2vsT1_coordinated_genes.csv')
    pos_df = pd.read_csv(DATA_PATH / 'T2vsT1_positive_correlation.csv')
    neg_df = pd.read_csv(DATA_PATH / 'T2vsT1_negative_correlation.csv')
    return coord_df, pos_df, neg_df

def is_transcription_factor(product):
    if pd.isna(product):
        return False
    product_lower = product.lower()
    return any(kw.lower() in product_lower for kw in TF_KEYWORDS)

def find_sarp_genes(df):
    """Find known SARP family genes."""
    sarp_matches = []
    for old_locus, description in SARP_GENES.items():
        matches = df[df['old_locus_tag'] == old_locus]
        if len(matches) > 0:
            for _, row in matches.iterrows():
                sarp_matches.append({
                    'gene_id': row['gene_id'],
                    'old_locus_tag': old_locus,
                    'sarp_description': description,
                    'mod_type': row['mod_type'],
                    'methyl_change': row['methyl_change'],
                    'log2FC': row['log2FC'],
                    'coordination': row['coordination']
                })
    return pd.DataFrame(sarp_matches)

def find_all_tfs(df, correlation_type):
    """Find all transcription factors."""
    tfs = []
    for _, row in df.iterrows():
        product = row.get('product', '')
        if is_transcription_factor(product):
            tfs.append({
                'gene_id': row['gene_id'],
                'old_locus_tag': row['old_locus_tag'],
                'product': product,
                'mod_type': row['mod_type'],
                'T1_sites': row['T1_sites'],
                'T2_sites': row['T2_sites'],
                'methyl_change': row['methyl_change'],
                'methyl_category': row['methyl_category'],
                'log2FC': row['log2FC'],
                'padj': row['padj'],
                'coordination': row['coordination'],
                'correlation_type': correlation_type
            })
    return pd.DataFrame(tfs)

def main():
    print("=" * 80)
    print("Transcription Factor Analysis (Including SARP Family)")
    print("=" * 80)

    coord_df, pos_df, neg_df = load_data()

    # Search for SARP genes
    print("\n[1] Searching for known SARP family genes...")
    sarp_in_coord = find_sarp_genes(coord_df)
    if len(sarp_in_coord) > 0:
        print(f"\nFound {len(sarp_in_coord)} SARP genes with coordinated changes:")
        for _, row in sarp_in_coord.iterrows():
            print(f"  {row['old_locus_tag']}: {row['sarp_description']}")
            print(f"    Methylation: {row['methyl_change']:+.1f}%, Expression: log2FC={row['log2FC']:+.2f}")
    else:
        print("\nNo known SARP genes found with coordinated methylation-expression changes.")
        print("(Note: actII-ORF4, redD, cdaR, cpkO may not have promoter methylation)")

    # Find all TFs
    print("\n[2] Finding all transcription factors...")

    pos_tfs = find_all_tfs(pos_df, 'Positive')
    neg_tfs = find_all_tfs(neg_df, 'Negative')
    all_tfs = pd.concat([pos_tfs, neg_tfs], ignore_index=True)

    print(f"\nTotal transcription factors found: {len(all_tfs)}")
    print(f"  - Positive correlation: {len(pos_tfs)}")
    print(f"  - Negative correlation: {len(neg_tfs)}")

    # Positive correlation TFs
    print("\n" + "=" * 80)
    print("POSITIVE CORRELATION TRANSCRIPTION FACTORS")
    print("(Methylation and expression change in same direction)")
    print("=" * 80)

    if len(pos_tfs) > 0:
        pos_tfs_sorted = pos_tfs.sort_values('log2FC', key=abs, ascending=False)
        print("\n{:<12} {:<10} {:>6} {:>8} {:>8} {:<35}".format(
            "Gene", "OldLocus", "Mod", "ΔMethyl", "log2FC", "Product"))
        print("-" * 90)
        for _, row in pos_tfs_sorted.iterrows():
            product = str(row['product'])[:35] if row['product'] else '-'
            print("{:<12} {:<10} {:>6} {:>+8.1f} {:>+8.2f} {:<35}".format(
                row['gene_id'][:12],
                str(row['old_locus_tag'])[:10] if pd.notna(row['old_locus_tag']) else '-',
                row['mod_type'],
                row['methyl_change'],
                row['log2FC'],
                product))

    # Negative correlation TFs
    print("\n" + "=" * 80)
    print("NEGATIVE CORRELATION TRANSCRIPTION FACTORS")
    print("(Methylation and expression change in opposite direction)")
    print("=" * 80)

    if len(neg_tfs) > 0:
        neg_tfs_sorted = neg_tfs.sort_values('log2FC', key=abs, ascending=False)
        print("\n{:<12} {:<10} {:>6} {:>8} {:>8} {:<35}".format(
            "Gene", "OldLocus", "Mod", "ΔMethyl", "log2FC", "Product"))
        print("-" * 90)
        for _, row in neg_tfs_sorted.iterrows():
            product = str(row['product'])[:35] if row['product'] else '-'
            print("{:<12} {:<10} {:>6} {:>+8.1f} {:>+8.2f} {:<35}".format(
                row['gene_id'][:12],
                str(row['old_locus_tag'])[:10] if pd.notna(row['old_locus_tag']) else '-',
                row['mod_type'],
                row['methyl_change'],
                row['log2FC'],
                product))

    # Key findings
    print("\n" + "=" * 80)
    print("KEY TRANSCRIPTION FACTORS")
    print("=" * 80)

    key_tfs = all_tfs[all_tfs['log2FC'].abs() > 2].sort_values('log2FC', key=abs, ascending=False)

    print(f"\n{len(key_tfs)} TFs with |log2FC| > 2:")
    for _, row in key_tfs.iterrows():
        direction = "↑" if row['log2FC'] > 0 else "↓"
        methyl_dir = "↑" if row['methyl_change'] > 0 else "↓" if row['methyl_change'] < 0 else "="
        print(f"\n  {row['old_locus_tag']} ({row['gene_id']})")
        print(f"    {row['product']}")
        print(f"    Methylation: {methyl_dir} ({row['methyl_change']:+.1f}%), Expression: {direction} (log2FC={row['log2FC']:+.2f})")
        print(f"    Pattern: {row['coordination']} ({row['correlation_type']} correlation)")

    # Save results
    all_tfs.to_csv(DATA_PATH / 'TF_coordinated_changes.csv', index=False)
    print(f"\n\nResults saved to: {DATA_PATH / 'TF_coordinated_changes.csv'}")

if __name__ == "__main__":
    main()
