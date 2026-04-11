#!/usr/bin/env python3
"""
Test Fang et al. 2022 model: Are transcription factor promoters enriched
for methylation changes compared to genome average?
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

BASE = Path('/Users/okaban/bioinfo/rna-seq')
GENE_MASTER = BASE / '05_annotation/analysis/05_annotation_260128_v1/tables/gene_master_DESeq2.tsv'
INTEGRATED = BASE / '11_epigenome_integration/analysis/01_integration/integrated_methyl_expression_weighted.csv'
TSS_TABLE = BASE / '11_epigenome_integration/analysis/18_tss_analyses/comprehensive_tss_table.csv'
METHYL_SITES = BASE / '11_epigenome_integration/analysis/01_integration/high_confidence_sites_weighted.csv'
BGC_REG = BASE / '05_annotation/analysis/05_annotation_260128_v1/tables/gene_master_with_BGC_regulators.tsv'
OUT_DIR = BASE / '11_epigenome_integration/analysis/18_tss_analyses'


def identify_tf_genes(gene_master):
    """Identify transcription factor / regulator genes from annotation."""
    tf_keywords = [
        'transcriptional regulator', 'transcription factor',
        'DNA-binding', 'response regulator', 'two-component',
        'sigma factor', 'MarR', 'TetR', 'LysR', 'GntR', 'AraC',
        'LacI', 'IclR', 'MerR', 'ArsR', 'WhiB', 'SARP',
        'regulatory protein', 'repressor', 'activator',
    ]
    tf_mask = gene_master['product'].fillna('').str.contains(
        '|'.join(tf_keywords), case=False, regex=True
    )
    return set(gene_master[tf_mask]['gene_id'].values)


def load_bgc_regulators():
    """Load BGC-associated regulator genes."""
    try:
        df = pd.read_csv(BGC_REG, sep='\t')
        # Check for regulator column
        for col in df.columns:
            if 'regulator' in col.lower() or 'is_regulator' in col.lower():
                return set(df[df[col].fillna(False).astype(bool)]['gene_id'].values)
        return set()
    except:
        return set()


def test_tf_methylation_enrichment():
    print("Loading data...")

    # Load gene master
    gene_master = pd.read_csv(GENE_MASTER, sep='\t')
    print(f"Total genes: {len(gene_master)}")

    # Identify TFs
    tf_genes = identify_tf_genes(gene_master)
    print(f"Transcription factors / regulators: {len(tf_genes)}")

    # Load integrated methylation-expression data
    integrated = pd.read_csv(INTEGRATED)
    print(f"Integrated genes: {len(integrated)}")

    # Load TSS table
    tss = pd.read_csv(TSS_TABLE)
    tss = tss.dropna(subset=['tss'])

    # Load methylation sites
    methyl = pd.read_csv(METHYL_SITES)

    # For each gene, count promoter methylation sites (-300 to +50 from TSS)
    print("\nCounting promoter methylation per gene...")
    methyl_unique = methyl.drop_duplicates(['position', 'mod_type', 'timepoint'])

    gene_methyl = {}
    for _, row in tss.iterrows():
        lt = row['gene_id']
        tss_pos = int(row['tss'])
        strand = row['strand']
        if strand == '+':
            prom_start = tss_pos - 300
            prom_end = tss_pos + 50
        else:
            prom_start = tss_pos - 50
            prom_end = tss_pos + 300

        sites = methyl_unique[
            (methyl_unique['position'] >= prom_start) &
            (methyl_unique['position'] < prom_end)
        ]
        gene_methyl[lt] = {
            'total_sites': len(sites),
            '6mA_sites': len(sites[sites['mod_type'] == '6mA']),
            '4mC_sites': len(sites[sites['mod_type'] == '4mC']),
            'is_tf': lt in tf_genes,
        }

    methyl_df = pd.DataFrame.from_dict(gene_methyl, orient='index')
    methyl_df.index.name = 'locus_tag'

    # Test 1: Do TFs have more promoter methylation sites?
    print("\n" + "=" * 60)
    print("TEST 1: TF promoter methylation density")
    print("=" * 60)

    tf_sites = methyl_df[methyl_df['is_tf']]['total_sites']
    nontf_sites = methyl_df[~methyl_df['is_tf']]['total_sites']

    print(f"\nTF genes (n={len(tf_sites)}):")
    print(f"  Mean promoter methyl sites: {tf_sites.mean():.2f}")
    print(f"  Median: {tf_sites.median():.1f}")
    print(f"  % with ≥1 site: {(tf_sites > 0).mean()*100:.1f}%")

    print(f"\nNon-TF genes (n={len(nontf_sites)}):")
    print(f"  Mean promoter methyl sites: {nontf_sites.mean():.2f}")
    print(f"  Median: {nontf_sites.median():.1f}")
    print(f"  % with ≥1 site: {(nontf_sites > 0).mean()*100:.1f}%")

    # Mann-Whitney U test
    u_stat, p_mw = stats.mannwhitneyu(tf_sites, nontf_sites, alternative='greater')
    print(f"\nMann-Whitney U (TF > non-TF): p = {p_mw:.6f}")

    # Fisher's test: TF with methylation vs without
    tf_meth = (tf_sites > 0).sum()
    tf_nometh = (tf_sites == 0).sum()
    nontf_meth = (nontf_sites > 0).sum()
    nontf_nometh = (nontf_sites == 0).sum()

    table = [[tf_meth, tf_nometh], [nontf_meth, nontf_nometh]]
    odds, p_fisher = stats.fisher_exact(table, alternative='greater')
    print(f"Fisher's (TF methylated rate > non-TF): OR={odds:.3f}, p={p_fisher:.6f}")

    # By mod type
    for mod in ['6mA', '4mC']:
        col = f'{mod}_sites'
        tf_m = methyl_df[methyl_df['is_tf']][col]
        nontf_m = methyl_df[~methyl_df['is_tf']][col]
        u, p = stats.mannwhitneyu(tf_m, nontf_m, alternative='greater')
        print(f"\n{mod}: TF mean={tf_m.mean():.2f}, non-TF mean={nontf_m.mean():.2f}, MW p={p:.6f}")

    # Test 2: Among DEGs, are TFs enriched in methylation-associated expression changes?
    print("\n" + "=" * 60)
    print("TEST 2: TF enrichment in methylation-correlated DEGs")
    print("=" * 60)

    # Find genes where methylation change and expression change are correlated
    # Use integrated data
    comparisons = ['T2_vs_T1', 'T3_vs_T1', 'T3_vs_T2']
    mod_types = ['6mA', '4mC']

    results = []
    for comp in comparisons:
        padj_col = f'padj_{comp}'
        lfc_col = f'log2FC_{comp}'

        if padj_col not in integrated.columns:
            continue

        for mod in mod_types:
            freq_col = f'{mod}_T1_mean_freq'
            if freq_col not in integrated.columns:
                continue

            # Genes with both methylation and significant expression change
            sub = integrated[
                (integrated[padj_col] < 0.05) &
                (integrated[freq_col].notna()) &
                (integrated[freq_col] > 0)
            ].copy()

            if len(sub) < 10:
                continue

            sub['is_tf'] = sub['gene_id'].isin(tf_genes)
            n_tf = sub['is_tf'].sum()
            n_total = len(sub)

            # Background TF rate
            all_degs = integrated[integrated[padj_col] < 0.05]
            all_degs_tf = all_degs['gene_id'].isin(tf_genes).sum()
            bg_rate = all_degs_tf / len(all_degs) if len(all_degs) > 0 else 0

            # Binomial test
            p_binom = stats.binomtest(n_tf, n_total, bg_rate, alternative='greater').pvalue

            results.append({
                'comparison': comp, 'mod_type': mod,
                'n_methyl_degs': n_total, 'n_tf': n_tf,
                'tf_rate': n_tf / n_total * 100,
                'bg_tf_rate': bg_rate * 100,
                'p_binom': p_binom,
            })
            print(f"{comp} {mod}: {n_tf}/{n_total} TFs ({n_tf/n_total*100:.1f}%) vs bg {bg_rate*100:.1f}% (p={p_binom:.4f})")

    # Test 3: Specific check — do methylated TFs control BGC expression?
    print("\n" + "=" * 60)
    print("TEST 3: Methylated TFs with BGC relevance")
    print("=" * 60)

    # Find TFs with promoter methylation that are also DEGs
    tf_with_methyl = methyl_df[(methyl_df['is_tf']) & (methyl_df['total_sites'] > 0)]
    tf_loci = set(tf_with_methyl.index)

    # Cross with expression data
    tf_expr = integrated[integrated['gene_id'].isin(tf_loci)].copy()

    # Check known BGC regulators
    known_regulators = {
        'SC_RS25685': 'actII-orf4 (Act pathway-specific activator)',
        'SC_RS26285': 'redD (Red pathway-specific activator)',
        'SC_RS26295': 'redZ (Red response regulator)',
        'SC_RS18490': 'cdaR (CDA SARP activator)',
        'SC_RS31380': 'cpkO (Cpk SARP activator)',
        'SC_RS25680': 'absA2 (response regulator)',
        'SC_RS25530': 'afsR (global SARP)',
        'SC_RS25540': 'afsS (sigma-like)',
    }

    print(f"\nMethylated TFs: {len(tf_loci)}")
    print(f"\nKnown BGC regulators with promoter methylation:")
    for locus, name in known_regulators.items():
        if locus in tf_loci:
            sites = tf_with_methyl.loc[locus]
            print(f"  {locus} ({name}): {int(sites['total_sites'])} sites "
                  f"(6mA:{int(sites['6mA_sites'])}, 4mC:{int(sites['4mC_sites'])})")

            # Expression changes
            expr_row = tf_expr[tf_expr['gene_id'] == locus]
            if len(expr_row) > 0:
                for comp in comparisons:
                    lfc_col = f'log2FC_{comp}'
                    padj_col = f'padj_{comp}'
                    if lfc_col in expr_row.columns:
                        lfc = expr_row[lfc_col].values[0]
                        padj = expr_row[padj_col].values[0]
                        sig = '*' if padj < 0.05 else 'ns'
                        print(f"    {comp}: log2FC={lfc:.2f} ({sig})")
        else:
            # Check if it exists in methyl_df at all
            if locus in methyl_df.index:
                print(f"  {locus} ({name}): NO promoter methylation")

    # Test 4: Temporal dynamics — do TFs gain/lose methylation differently?
    print("\n" + "=" * 60)
    print("TEST 4: TF methylation temporal dynamics")
    print("=" * 60)

    # Count per-timepoint promoter methylation for TFs vs non-TFs
    for tp in ['T1', 'T2', 'T3']:
        tp_methyl = methyl[methyl['timepoint'] == tp].drop_duplicates(['position', 'mod_type'])
        tf_counts = []
        nontf_counts = []

        for _, row in tss.iterrows():
            lt = row['gene_id']
            tss_pos = int(row['tss'])
            strand = row['strand']
            if strand == '+':
                prom_start = tss_pos - 300
                prom_end = tss_pos + 50
            else:
                prom_start = tss_pos - 50
                prom_end = tss_pos + 300

            n_sites = len(tp_methyl[
                (tp_methyl['position'] >= prom_start) &
                (tp_methyl['position'] < prom_end)
            ])

            if lt in tf_genes:
                tf_counts.append(n_sites)
            else:
                nontf_counts.append(n_sites)

        tf_mean = np.mean(tf_counts)
        nontf_mean = np.mean(nontf_counts)
        u, p = stats.mannwhitneyu(tf_counts, nontf_counts, alternative='greater')
        print(f"{tp}: TF mean={tf_mean:.2f}, non-TF mean={nontf_mean:.2f}, "
              f"ratio={tf_mean/nontf_mean:.2f}, MW p={p:.4f}")

    # Figure
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    # Panel A: TF vs non-TF methylation distribution
    bins = np.arange(0, 20, 1)
    axes[0].hist(nontf_sites, bins=bins, alpha=0.6, label=f'Non-TF (n={len(nontf_sites)})',
                  color='#3498DB', density=True)
    axes[0].hist(tf_sites, bins=bins, alpha=0.6, label=f'TF (n={len(tf_sites)})',
                  color='#E74C3C', density=True)
    axes[0].set_xlabel('Promoter methylation sites')
    axes[0].set_ylabel('Density')
    axes[0].set_title(f'A. Promoter methylation density\n(MW p={p_mw:.4f})')
    axes[0].legend(fontsize=8)

    # Panel B: % methylated by category
    categories = ['TF', 'Non-TF']
    pct_methylated = [
        (tf_sites > 0).mean() * 100,
        (nontf_sites > 0).mean() * 100
    ]
    mean_sites = [tf_sites.mean(), nontf_sites.mean()]

    x = np.arange(2)
    bars = axes[1].bar(x, pct_methylated, color=['#E74C3C', '#3498DB'], alpha=0.8)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(categories)
    axes[1].set_ylabel('% genes with promoter methylation')
    axes[1].set_title(f'B. % genes with ≥1 promoter methyl site\n(Fisher p={p_fisher:.4f})')
    for i, (pct, bar) in enumerate(zip(pct_methylated, bars)):
        axes[1].text(bar.get_x() + bar.get_width()/2, pct + 0.5,
                     f'{pct:.1f}%', ha='center', fontsize=10)

    # Panel C: By modification type
    mod_data = {}
    for mod in ['6mA', '4mC']:
        col = f'{mod}_sites'
        tf_pct = (methyl_df[methyl_df['is_tf']][col] > 0).mean() * 100
        nontf_pct = (methyl_df[~methyl_df['is_tf']][col] > 0).mean() * 100
        mod_data[mod] = {'TF': tf_pct, 'Non-TF': nontf_pct}

    x = np.arange(2)
    w = 0.35
    axes[2].bar(x - w/2, [mod_data['6mA']['TF'], mod_data['6mA']['Non-TF']],
                w, label='6mA', color='#E74C3C', alpha=0.8)
    axes[2].bar(x + w/2, [mod_data['4mC']['TF'], mod_data['4mC']['Non-TF']],
                w, label='4mC', color='#3498DB', alpha=0.8)
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(categories)
    axes[2].set_ylabel('% genes with promoter methylation')
    axes[2].set_title('C. By modification type')
    axes[2].legend(fontsize=8)

    plt.tight_layout()
    fig.savefig(OUT_DIR / 'tf_methylation_enrichment.png', dpi=150, bbox_inches='tight')
    fig.savefig(OUT_DIR / 'tf_methylation_enrichment.pdf', bbox_inches='tight')
    plt.close()
    print(f"\nSaved: tf_methylation_enrichment.png/pdf")

    # Save summary
    summary = pd.DataFrame(results) if results else pd.DataFrame()
    summary.to_csv(OUT_DIR / 'tf_methylation_enrichment_summary.csv', index=False)

    # Save methylated TF list
    tf_detail = methyl_df[methyl_df['is_tf']].copy()
    tf_detail = tf_detail[tf_detail['total_sites'] > 0].sort_values('total_sites', ascending=False)
    tf_detail.to_csv(OUT_DIR / 'methylated_tf_list.csv')
    print(f"Saved: methylated_tf_list.csv ({len(tf_detail)} TFs with promoter methylation)")


if __name__ == '__main__':
    test_tf_methylation_enrichment()
