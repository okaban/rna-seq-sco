#!/usr/bin/env python3
"""
BGC (Biosynthetic Gene Cluster) Analysis for Methylation-Expression Correlation
Focus on Act, Red, CDA, Cpk and other known BGCs in S. coelicolor
"""

import pandas as pd
import numpy as np
from pathlib import Path

# Paths
DATA_PATH = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis")
OUTPUT_PATH = DATA_PATH

# Known BGC regions in S. coelicolor A3(2) M145
# Based on antiSMASH and literature annotations
# Format: (name, start_SCO, end_SCO, description)
KNOWN_BGCS = {
    'Act': {
        'description': 'Actinorhodin (type II PKS, blue pigment)',
        'old_locus_range': ('SCO5071', 'SCO5092'),
        'genes': ['SCO5071', 'SCO5072', 'SCO5073', 'SCO5074', 'SCO5075', 'SCO5076',
                 'SCO5077', 'SCO5078', 'SCO5079', 'SCO5080', 'SCO5081', 'SCO5082',
                 'SCO5083', 'SCO5084', 'SCO5085', 'SCO5086', 'SCO5087', 'SCO5088',
                 'SCO5089', 'SCO5090', 'SCO5091', 'SCO5092'],
        'key_genes': {
            'SCO5085': 'actII-ORF4 (pathway-specific activator)',
            'SCO5087': 'actVA-ORF5 (tailoring)',
            'SCO5090': 'actVI-ORFA (ketoreductase)',
        }
    },
    'Red': {
        'description': 'Undecylprodigiosin (prodiginine, red pigment)',
        'old_locus_range': ('SCO5877', 'SCO5898'),
        'genes': ['SCO5877', 'SCO5878', 'SCO5879', 'SCO5880', 'SCO5881', 'SCO5882',
                 'SCO5883', 'SCO5884', 'SCO5885', 'SCO5886', 'SCO5887', 'SCO5888',
                 'SCO5889', 'SCO5890', 'SCO5891', 'SCO5892', 'SCO5893', 'SCO5894',
                 'SCO5895', 'SCO5896', 'SCO5897', 'SCO5898'],
        'key_genes': {
            'SCO5881': 'redZ (response regulator)',
            'SCO5877': 'redD (pathway-specific activator)',
        }
    },
    'CDA': {
        'description': 'Calcium-Dependent Antibiotic (lipopeptide)',
        'old_locus_range': ('SCO3210', 'SCO3249'),
        'genes': [f'SCO{i}' for i in range(3210, 3250)],
        'key_genes': {
            'SCO3217': 'cdaR (pathway-specific activator)',
            'SCO3230': 'cdaPS1 (NRPS module)',
        }
    },
    'Cpk': {
        'description': 'Cryptic Type I PKS (coelimycin P1)',
        'old_locus_range': ('SCO6269', 'SCO6288'),
        'genes': [f'SCO{i}' for i in range(6269, 6289)],
        'key_genes': {
            'SCO6269': 'cpkO (TetR regulator)',
            'SCO6273': 'cpkA (PKS)',
        }
    },
    'Germicidin': {
        'description': 'Germicidin (type III PKS)',
        'old_locus_range': ('SCO7221', 'SCO7222'),
        'genes': ['SCO7221', 'SCO7222'],
        'key_genes': {}
    },
    'Hopanoid': {
        'description': 'Hopanoid biosynthesis',
        'old_locus_range': ('SCO6759', 'SCO6771'),
        'genes': [f'SCO{i}' for i in range(6759, 6772)],
        'key_genes': {}
    },
    'Desferrioxamine': {
        'description': 'Desferrioxamine (siderophore)',
        'old_locus_range': ('SCO2782', 'SCO2785'),
        'genes': ['SCO2782', 'SCO2783', 'SCO2784', 'SCO2785'],
        'key_genes': {
            'SCO2782': 'desA (decarboxylase)',
        }
    },
    'Coelichelin': {
        'description': 'Coelichelin (siderophore, NRPS)',
        'old_locus_range': ('SCO0489', 'SCO0499'),
        'genes': [f'SCO0{i}' for i in range(489, 500)],
        'key_genes': {}
    },
    'Geosmin': {
        'description': 'Geosmin (terpene, earthy smell)',
        'old_locus_range': ('SCO6073', 'SCO6073'),
        'genes': ['SCO6073'],
        'key_genes': {}
    },
    'Albaflavenone': {
        'description': 'Albaflavenone (sesquiterpene)',
        'old_locus_range': ('SCO5222', 'SCO5223'),
        'genes': ['SCO5222', 'SCO5223'],
        'key_genes': {}
    },
    'SapB': {
        'description': 'SapB (lanthipeptide, morphogenetic)',
        'old_locus_range': ('SCO6681', 'SCO6685'),
        'genes': ['SCO6681', 'SCO6682', 'SCO6683', 'SCO6684', 'SCO6685'],
        'key_genes': {
            'SCO6682': 'ramS (SapB precursor)',
            'SCO6685': 'ramR (response regulator)',
        }
    },
    'Melanin': {
        'description': 'Melanin pigment',
        'old_locus_range': ('SCO2700', 'SCO2701'),
        'genes': ['SCO2700', 'SCO2701'],
        'key_genes': {}
    },
    'Eicosapentaenoic_acid': {
        'description': 'PUFA synthase cluster',
        'old_locus_range': ('SCO0124', 'SCO0129'),
        'genes': [f'SCO0{i}' for i in range(124, 130)],
        'key_genes': {}
    },
}

def load_coordinated_genes():
    """Load the coordinated genes data."""
    coord_df = pd.read_csv(DATA_PATH / 'T2vsT1_coordinated_genes.csv')
    pos_df = pd.read_csv(DATA_PATH / 'T2vsT1_positive_correlation.csv')
    neg_df = pd.read_csv(DATA_PATH / 'T2vsT1_negative_correlation.csv')
    return coord_df, pos_df, neg_df

def find_bgc_genes(coord_df, bgc_name, bgc_info):
    """Find genes from a BGC in the coordinated genes list."""
    bgc_genes = bgc_info['genes']

    # Match by old_locus_tag
    matches = coord_df[coord_df['old_locus_tag'].isin(bgc_genes)]
    return matches

def analyze_all_bgcs(coord_df, pos_df, neg_df):
    """Analyze all known BGCs."""
    results = []
    detailed_results = []

    for bgc_name, bgc_info in KNOWN_BGCS.items():
        # Find in all coordinated genes
        all_matches = find_bgc_genes(coord_df, bgc_name, bgc_info)
        pos_matches = find_bgc_genes(pos_df, bgc_name, bgc_info)
        neg_matches = find_bgc_genes(neg_df, bgc_name, bgc_info)

        results.append({
            'BGC': bgc_name,
            'Description': bgc_info['description'],
            'Total_genes_in_cluster': len(bgc_info['genes']),
            'Coordinated_genes': len(all_matches),
            'Positive_correlation': len(pos_matches),
            'Negative_correlation': len(neg_matches),
        })

        # Store detailed info for each match
        for _, row in all_matches.iterrows():
            old_locus = row['old_locus_tag']
            is_key = old_locus in bgc_info.get('key_genes', {})
            key_desc = bgc_info.get('key_genes', {}).get(old_locus, '')

            detailed_results.append({
                'BGC': bgc_name,
                'gene_id': row['gene_id'],
                'old_locus_tag': old_locus,
                'mod_type': row['mod_type'],
                'T1_sites': row['T1_sites'],
                'T2_sites': row['T2_sites'],
                'methyl_change': row['methyl_change'],
                'methyl_category': row['methyl_category'],
                'log2FC': row['log2FC'],
                'padj': row['padj'],
                'expr_category': row['expr_category'],
                'coordination': row['coordination'],
                'is_key_gene': is_key,
                'key_gene_function': key_desc,
                'product': row['product'] if 'product' in row else ''
            })

    return pd.DataFrame(results), pd.DataFrame(detailed_results)

def main():
    print("=" * 80)
    print("BGC Analysis: Methylation-Expression Correlation in Biosynthetic Gene Clusters")
    print("=" * 80)

    # Load data
    print("\n[1] Loading coordinated genes data...")
    coord_df, pos_df, neg_df = load_coordinated_genes()
    print(f"  Total coordinated genes: {len(coord_df)}")

    # Analyze BGCs
    print("\n[2] Analyzing known BGCs...")
    summary_df, detailed_df = analyze_all_bgcs(coord_df, pos_df, neg_df)

    # Print summary
    print("\n" + "=" * 80)
    print("BGC SUMMARY")
    print("=" * 80)

    print("\n{:<20} {:<12} {:>8} {:>8} {:>8}".format(
        "BGC", "Type", "Total", "Pos", "Neg"))
    print("-" * 60)

    for _, row in summary_df.sort_values('Coordinated_genes', ascending=False).iterrows():
        if row['Coordinated_genes'] > 0:
            print("{:<20} {:<12} {:>8} {:>8} {:>8}".format(
                row['BGC'],
                row['Description'][:12],
                row['Coordinated_genes'],
                row['Positive_correlation'],
                row['Negative_correlation']))

    # Detailed analysis for Act and Red
    print("\n" + "=" * 80)
    print("DETAILED BGC ANALYSIS")
    print("=" * 80)

    for bgc in ['Act', 'Red', 'CDA', 'Cpk', 'SapB', 'Desferrioxamine']:
        bgc_detail = detailed_df[detailed_df['BGC'] == bgc]
        if len(bgc_detail) > 0:
            print(f"\n### {bgc} ({KNOWN_BGCS[bgc]['description']}) ###")
            print(f"Genes with coordinated methylation-expression changes: {len(bgc_detail)}")
            print()

            for _, row in bgc_detail.iterrows():
                key_marker = " [KEY]" if row['is_key_gene'] else ""
                print(f"  {row['old_locus_tag']} ({row['gene_id']}){key_marker}")
                print(f"    Modification: {row['mod_type']}, Sites: {row['T1_sites']}→{row['T2_sites']}")
                print(f"    Methylation: {row['methyl_change']:+.1f}% ({row['methyl_category']})")
                print(f"    Expression: log2FC={row['log2FC']:+.2f}, padj={row['padj']:.2e} ({row['expr_category']})")
                print(f"    Coordination: {row['coordination']}")
                if row['is_key_gene']:
                    print(f"    Function: {row['key_gene_function']}")
                if row['product']:
                    print(f"    Product: {row['product']}")
                print()
        else:
            print(f"\n### {bgc} ({KNOWN_BGCS[bgc]['description']}) ###")
            print("No genes with coordinated methylation-expression changes found.")

    # Save results
    print("\n[3] Saving results...")
    summary_df.to_csv(OUTPUT_PATH / 'BGC_summary.csv', index=False)
    if len(detailed_df) > 0:
        detailed_df.to_csv(OUTPUT_PATH / 'BGC_detailed.csv', index=False)

    print(f"\nResults saved to: {OUTPUT_PATH}")

    # Return for further analysis
    return summary_df, detailed_df

if __name__ == "__main__":
    summary_df, detailed_df = main()
