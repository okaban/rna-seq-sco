#!/usr/bin/env python3
"""
DMG Selectivity Deep Analysis
==============================
5つの追加解析で、メチル化変動の機能選択性を多角的に検証する。

解析1: DEG∩DMG（協調変動遺伝子）のKEGG/GO/COGエンリッチメント
解析2: COGカテゴリ別メチル化-発現相関
解析3: プロモーター vs 遺伝子本体メチル化の機能比較
解析4: V-Defense MTase群・GO isomerase群の発現変動検証
解析5: TF/BGC遺伝子に限定したDMG頻度検定
"""

import os
import warnings
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
from pathlib import Path

warnings.filterwarnings('ignore')

# ============================================================
# Configuration
# ============================================================
PROJECT_ROOT = Path("/Users/okaban/bioinfo/rna-seq")
ANALYSIS_DIR = PROJECT_ROOT / "14_DMG_functional_enrichment" / "analysis" / "14_DMG_selectivity_260206_v1"
ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
(ANALYSIS_DIR / "figures").mkdir(exist_ok=True)
(ANALYSIS_DIR / "tables").mkdir(exist_ok=True)

# Input files
INTEGRATED_METHYL = PROJECT_ROOT / "11_epigenome_integration" / "analysis" / "01_integration" / "integrated_methyl_expression_weighted.csv"
GENE_MASTER = PROJECT_ROOT / "05_annotation" / "analysis" / "05_annotation_260128_v1" / "tables" / "gene_master_DESeq2.tsv"
GENE_ANNOT = PROJECT_ROOT / "05_annotation" / "analysis" / "05_annotation_260128_v1" / "tables" / "gene_annotation_basic.tsv"
COG_CLASS = PROJECT_ROOT / "12_supplementary_figures" / "analysis" / "12_supplementary_260202_v1" / "tables" / "gene_COG_classification.tsv"
BGC_DEF = PROJECT_ROOT / "05_annotation" / "analysis" / "05_annotation_260128_v1" / "tables" / "BGC_definition_manual.tsv"
GENE_MASTER_BGC_REG = PROJECT_ROOT / "05_annotation" / "analysis" / "05_annotation_260128_v1" / "tables" / "gene_master_with_BGC_regulators.tsv"
COORDINATED_T2T1 = PROJECT_ROOT / "11_epigenome_integration" / "analysis" / "01_integration" / "T2vsT1_coordinated_genes.csv"
COORDINATED_T3T1 = PROJECT_ROOT / "11_epigenome_integration" / "analysis" / "16_t3_coordinated" / "T3vsT1_coordinated_genes.csv"
COORDINATED_T3T2 = PROJECT_ROOT / "11_epigenome_integration" / "analysis" / "16_t3_coordinated" / "T3vsT2_coordinated_genes.csv"
HIGH_CONF_SITES = PROJECT_ROOT / "11_epigenome_integration" / "analysis" / "01_integration" / "high_confidence_sites_weighted.csv"
TSS_TABLE = PROJECT_ROOT / "11_epigenome_integration" / "analysis" / "18_tss_analyses" / "comprehensive_tss_table.csv"
KEGG_CACHE = PROJECT_ROOT / "14_DMG_functional_enrichment" / "analysis" / "14_DMG_enrichment_260206_v1" / "tables" / "kegg_cache.tsv"

DMG_CHANGE_THRESHOLD = 10.0
DEG_PADJ = 0.05
DEG_LOG2FC = 1.0
MIN_GENES = 3

plt.rcParams.update({
    'font.size': 10, 'axes.titlesize': 12, 'axes.labelsize': 10,
    'figure.dpi': 150, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
    'font.family': 'sans-serif',
})

COG_NAMES = {
    'C': 'Energy production', 'D': 'Cell division', 'E': 'Amino acid metabolism',
    'F': 'Nucleotide metabolism', 'G': 'Carbohydrate metabolism', 'H': 'Coenzyme metabolism',
    'I': 'Lipid metabolism', 'J': 'Translation', 'K': 'Transcription',
    'L': 'Replication/Repair', 'M': 'Cell wall', 'N': 'Motility',
    'O': 'Chaperones/PTM', 'P': 'Inorganic ion transport', 'Q': 'Secondary metabolism',
    'R': 'General function', 'S': 'Unknown', 'T': 'Signal transduction',
    'U': 'Secretion', 'V': 'Defense',
}


def save_figure(fig, name):
    fig.savefig(ANALYSIS_DIR / "figures" / f"{name}.pdf", format='pdf')
    fig.savefig(ANALYSIS_DIR / "figures" / f"{name}.svg", format='svg')
    plt.close(fig)
    print(f"  Saved: {name}.pdf / .svg")


def run_ora(gene_list, background, term_to_genes, term_names, min_genes=MIN_GENES):
    gene_set = set(gene_list)
    bg_set = set(background)
    n_total = len(bg_set)
    n_list = len(gene_set & bg_set)
    results = []
    for term_id, term_genes_raw in term_to_genes.items():
        term_genes = term_genes_raw & bg_set
        n_term = len(term_genes)
        if n_term < min_genes:
            continue
        overlap = gene_set & term_genes
        n_overlap = len(overlap)
        if n_overlap < 1:
            continue
        a, b = n_overlap, n_list - n_overlap
        c, d = n_term - n_overlap, n_total - n_list - n_term + n_overlap
        odds_ratio, p_value = stats.fisher_exact([[a, b], [c, d]], alternative='greater')
        gene_ratio = n_overlap / n_list if n_list > 0 else 0
        bg_ratio = n_term / n_total if n_total > 0 else 0
        fold_enrichment = gene_ratio / bg_ratio if bg_ratio > 0 else 0
        results.append({
            'term_id': term_id, 'term_name': term_names.get(term_id, term_id),
            'count': n_overlap, 'total_in_term': n_term,
            'gene_ratio': f"{n_overlap}/{n_list}", 'fold_enrichment': round(fold_enrichment, 3),
            'pvalue': p_value, 'genes': ';'.join(sorted(overlap))
        })
    if not results:
        return pd.DataFrame()
    result_df = pd.DataFrame(results)
    if len(result_df) > 0:
        _, padj, _, _ = multipletests(result_df['pvalue'].values, method='fdr_bh')
        result_df['padj'] = padj
    result_df = result_df.sort_values('pvalue').reset_index(drop=True)
    return result_df


# ============================================================
# Load common data
# ============================================================
def load_all_data():
    print("=" * 60)
    print("Loading data...")
    print("=" * 60)

    methyl = pd.read_csv(INTEGRATED_METHYL)
    gene_master = pd.read_csv(GENE_MASTER, sep='\t')
    gene_annot = pd.read_csv(GENE_ANNOT, sep='\t')
    cog = pd.read_csv(COG_CLASS, sep='\t')
    cog['COG_letter'] = cog['COG_category'].str.extract(r'^([A-Z])')

    # KEGG mapping
    kegg_cache = pd.read_csv(KEGG_CACHE, sep='\t')
    pathway_to_genes = defaultdict(set)
    pathway_names = {}
    for _, row in kegg_cache.iterrows():
        pathway_to_genes[row['pathway_id']].add(row['gene'])
        pathway_names[row['pathway_id']] = row['pathway_name']
    id_map = gene_annot[['gene_id', 'old_locus_tag']].dropna().set_index('gene_id')['old_locus_tag'].to_dict()
    kegg_genes = set(kegg_cache['gene'].unique())
    all_sco_bg = set(id_map.values()) & kegg_genes

    # GO mapping
    go_mapping = gene_annot[['gene_id', 'ontology_term']].dropna(subset=['ontology_term'])
    go_to_genes = defaultdict(set)
    go_bg = set()
    for _, row in go_mapping.iterrows():
        gene = row['gene_id']
        terms = [t.strip() for t in str(row['ontology_term']).split(',') if t.strip().startswith('GO:')]
        if terms:
            go_bg.add(gene)
            for t in terms:
                go_to_genes[t].add(gene)
    go_names = {t: t for t in go_to_genes}

    # COG mapping
    cog_to_genes = defaultdict(set)
    for _, row in cog.dropna(subset=['COG_letter']).iterrows():
        cog_to_genes[row['COG_letter']].add(row['gene_id'])
    cog_bg = set(cog['gene_id'].unique())

    print(f"  Methylation data: {len(methyl)} genes")
    print(f"  KEGG background: {len(all_sco_bg)} genes")
    print(f"  GO background: {len(go_bg)} genes")
    print(f"  COG background: {len(cog_bg)} genes")

    return {
        'methyl': methyl, 'gene_master': gene_master, 'gene_annot': gene_annot,
        'cog': cog, 'id_map': id_map,
        'pathway_to_genes': pathway_to_genes, 'pathway_names': pathway_names, 'all_sco_bg': all_sco_bg,
        'go_to_genes': go_to_genes, 'go_names': go_names, 'go_bg': go_bg,
        'cog_to_genes': cog_to_genes, 'cog_bg': cog_bg, 'COG_NAMES': COG_NAMES,
    }


# ============================================================
# Analysis 1: DEG∩DMG Coordinated Gene Enrichment
# ============================================================
def analysis1_coordinated_enrichment(data):
    print("\n" + "=" * 60)
    print("Analysis 1: DEG∩DMG Coordinated Gene Enrichment")
    print("=" * 60)

    coord_files = {
        'T2vsT1': COORDINATED_T2T1,
        'T3vsT1': COORDINATED_T3T1,
        'T3vsT2': COORDINATED_T3T2,
    }

    all_results = {}
    summary_rows = []

    for comp, fpath in coord_files.items():
        if not fpath.exists():
            print(f"  {comp}: file not found, skipping")
            continue

        df = pd.read_csv(fpath)
        print(f"\n  {comp}: {len(df)} entries loaded")

        # Filter to significant coordinated genes (DEG AND DMG)
        if 'padj' in df.columns:
            sig = df[(df['padj'] < DEG_PADJ) & (df['methyl_change'].abs() > DMG_CHANGE_THRESHOLD)].copy()
        else:
            sig = df[df['methyl_change'].abs() > DMG_CHANGE_THRESHOLD].copy()

        if 'correlation_type' not in sig.columns and 'coordination' in sig.columns:
            # Derive correlation_type from coordination
            sig['correlation_type'] = sig['coordination'].apply(
                lambda x: 'Positive' if x in ['Gained_Up', 'Lost_Down'] else
                          'Negative' if x in ['Gained_Down', 'Lost_Up'] else 'Other'
            )

        # Get unique genes per group
        groups = {
            'all': sig['gene_id'].unique().tolist(),
            'positive': sig[sig['correlation_type'] == 'Positive']['gene_id'].unique().tolist(),
            'negative': sig[sig['correlation_type'] == 'Negative']['gene_id'].unique().tolist(),
        }

        if 'methyl_direction' in sig.columns:
            groups['gained'] = sig[sig['methyl_direction'] == 'Gained']['gene_id'].unique().tolist()
            groups['lost'] = sig[sig['methyl_direction'] == 'Lost']['gene_id'].unique().tolist()
        elif 'methyl_category' in sig.columns:
            groups['gained'] = sig[sig['methyl_category'] == 'Gained']['gene_id'].unique().tolist()
            groups['lost'] = sig[sig['methyl_category'] == 'Lost']['gene_id'].unique().tolist()

        for group_name, genes in groups.items():
            if len(genes) < MIN_GENES:
                continue

            # KEGG
            sco_genes = [data['id_map'][g] for g in genes if g in data['id_map'] and data['id_map'][g] in data['all_sco_bg']]
            if len(sco_genes) >= MIN_GENES:
                kegg_res = run_ora(sco_genes, data['all_sco_bg'], data['pathway_to_genes'], data['pathway_names'])
                key = f"coordinated_{comp}_{group_name}_KEGG"
                all_results[key] = kegg_res
                n_sig = len(kegg_res[kegg_res['padj'] < 0.05]) if not kegg_res.empty else 0
                print(f"    {comp} {group_name} KEGG: {len(sco_genes)} genes -> {n_sig} sig pathways")
                if not kegg_res.empty:
                    kegg_res.to_csv(ANALYSIS_DIR / "tables" / f"A1_KEGG_{comp}_{group_name}.tsv", sep='\t', index=False)

            # GO
            genes_go = [g for g in genes if g in data['go_bg']]
            if len(genes_go) >= MIN_GENES:
                go_res = run_ora(genes_go, data['go_bg'], data['go_to_genes'], data['go_names'])
                key = f"coordinated_{comp}_{group_name}_GO"
                all_results[key] = go_res
                n_sig = len(go_res[go_res['padj'] < 0.05]) if not go_res.empty else 0
                print(f"    {comp} {group_name} GO: {len(genes_go)} genes -> {n_sig} sig terms")
                if not go_res.empty:
                    go_res.to_csv(ANALYSIS_DIR / "tables" / f"A1_GO_{comp}_{group_name}.tsv", sep='\t', index=False)

            # COG
            genes_cog = [g for g in genes if g in data['cog_bg']]
            if len(genes_cog) >= MIN_GENES:
                cog_res = run_ora(genes_cog, data['cog_bg'], data['cog_to_genes'], COG_NAMES, min_genes=1)
                key = f"coordinated_{comp}_{group_name}_COG"
                all_results[key] = cog_res
                n_sig = len(cog_res[cog_res['padj'] < 0.05]) if not cog_res.empty else 0
                print(f"    {comp} {group_name} COG: {len(genes_cog)} genes -> {n_sig} sig categories")
                if not cog_res.empty:
                    cog_res.to_csv(ANALYSIS_DIR / "tables" / f"A1_COG_{comp}_{group_name}.tsv", sep='\t', index=False)

            summary_rows.append({
                'comparison': comp, 'group': group_name, 'n_genes': len(genes),
                'n_KEGG_sig': len(all_results.get(f"coordinated_{comp}_{group_name}_KEGG", pd.DataFrame()).query('padj < 0.05')) if f"coordinated_{comp}_{group_name}_KEGG" in all_results and not all_results[f"coordinated_{comp}_{group_name}_KEGG"].empty else 0,
                'n_GO_sig': len(all_results.get(f"coordinated_{comp}_{group_name}_GO", pd.DataFrame()).query('padj < 0.05')) if f"coordinated_{comp}_{group_name}_GO" in all_results and not all_results[f"coordinated_{comp}_{group_name}_GO"].empty else 0,
                'n_COG_sig': len(all_results.get(f"coordinated_{comp}_{group_name}_COG", pd.DataFrame()).query('padj < 0.05')) if f"coordinated_{comp}_{group_name}_COG" in all_results and not all_results[f"coordinated_{comp}_{group_name}_COG"].empty else 0,
            })

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(ANALYSIS_DIR / "tables" / "A1_coordinated_enrichment_summary.tsv", sep='\t', index=False)
    print("\n  Summary:")
    print(summary_df.to_string(index=False))

    return all_results, summary_df


# ============================================================
# Analysis 2: COG Category-Specific Methylation-Expression Correlation
# ============================================================
def analysis2_category_correlation(data):
    print("\n" + "=" * 60)
    print("Analysis 2: COG Category-Specific Correlation")
    print("=" * 60)

    methyl = data['methyl']
    cog = data['cog']

    # Merge methylation data with COG
    merged = methyl.merge(cog[['gene_id', 'COG_letter']].dropna(), on='gene_id', how='inner')

    comparisons = {
        'T2vsT1': {'methyl_cols': ['4mC_change_T2_vs_T1', '6mA_change_T2_vs_T1'],
                    'expr_col': 'log2FC_T2_vs_T1', 'padj_col': 'padj_T2_vs_T1'},
        'T3vsT1': {'methyl_cols': ['4mC_change_T3_vs_T1', '6mA_change_T3_vs_T1'],
                    'expr_col': 'log2FC_T3_vs_T1', 'padj_col': 'padj_T3_vs_T1'},
        'T3vsT2': {'methyl_cols': ['4mC_change_T3_vs_T2', '6mA_change_T3_vs_T2'],
                    'expr_col': 'log2FC_T3_vs_T2', 'padj_col': 'padj_T3_vs_T2'},
    }

    all_results = []
    categories = sorted([c for c in merged['COG_letter'].unique() if c != 'S'])

    for comp, cols in comparisons.items():
        for mod_col in cols['methyl_cols']:
            mod_type = '4mC' if '4mC' in mod_col else '6mA'
            for cat in categories:
                subset = merged[merged['COG_letter'] == cat].dropna(subset=[mod_col, cols['expr_col']])
                # Filter to genes with some methylation change
                subset = subset[subset[mod_col].abs() > 0]
                n = len(subset)
                if n < 10:
                    all_results.append({
                        'comparison': comp, 'mod_type': mod_type, 'COG': cat,
                        'COG_name': COG_NAMES.get(cat, ''), 'n_genes': n,
                        'spearman_r': np.nan, 'pvalue': np.nan
                    })
                    continue

                r, p = stats.spearmanr(subset[mod_col], subset[cols['expr_col']])
                all_results.append({
                    'comparison': comp, 'mod_type': mod_type, 'COG': cat,
                    'COG_name': COG_NAMES.get(cat, ''), 'n_genes': n,
                    'spearman_r': round(r, 4), 'pvalue': p
                })

        # Also compute whole-genome correlation for reference
        for mod_col in cols['methyl_cols']:
            mod_type = '4mC' if '4mC' in mod_col else '6mA'
            subset_all = merged.dropna(subset=[mod_col, cols['expr_col']])
            subset_all = subset_all[subset_all[mod_col].abs() > 0]
            if len(subset_all) >= 10:
                r, p = stats.spearmanr(subset_all[mod_col], subset_all[cols['expr_col']])
                all_results.append({
                    'comparison': comp, 'mod_type': mod_type, 'COG': 'ALL',
                    'COG_name': 'All genes', 'n_genes': len(subset_all),
                    'spearman_r': round(r, 4), 'pvalue': p
                })

    result_df = pd.DataFrame(all_results)
    result_df.to_csv(ANALYSIS_DIR / "tables" / "A2_category_correlation.tsv", sep='\t', index=False)

    # Heatmap: 4mC correlations across COG categories and comparisons
    for mod_type in ['4mC', '6mA']:
        subset = result_df[(result_df['mod_type'] == mod_type) & (result_df['COG'] != 'ALL')]
        if subset.empty:
            continue

        pivot = subset.pivot_table(index='COG', columns='comparison', values='spearman_r')
        pivot_p = subset.pivot_table(index='COG', columns='comparison', values='pvalue')

        # Add COG names
        pivot.index = [f"{c} - {COG_NAMES.get(c, '')}" for c in pivot.index]

        fig, ax = plt.subplots(figsize=(8, max(6, len(pivot) * 0.35)))
        mask = pivot.isna()
        sns.heatmap(pivot, annot=True, fmt='.3f', cmap='RdBu_r', center=0,
                    vmin=-0.4, vmax=0.4, mask=mask, ax=ax,
                    linewidths=0.5, cbar_kws={'label': 'Spearman r'})

        # Mark significant cells
        for i, cat in enumerate(pivot.index):
            for j, comp in enumerate(pivot.columns):
                orig_cat = cat.split(' - ')[0]
                p_row = pivot_p[(pivot_p.index == orig_cat)]
                if not p_row.empty and comp in p_row.columns:
                    p_val = p_row[comp].values[0]
                    if not np.isnan(p_val) and p_val < 0.05:
                        ax.text(j + 0.5, i + 0.85, '*', ha='center', va='center',
                                fontsize=14, fontweight='bold', color='black')

        ax.set_title(f'{mod_type} Methylation-Expression Correlation by COG Category\n(* p < 0.05)')
        fig.tight_layout()
        save_figure(fig, f'A2_correlation_heatmap_{mod_type}')

    # Print key findings
    sig_results = result_df[(result_df['pvalue'] < 0.05) & (result_df['COG'] != 'ALL')]
    print(f"\n  Significant category-specific correlations (p<0.05): {len(sig_results)}")
    if not sig_results.empty:
        print(sig_results[['comparison', 'mod_type', 'COG', 'COG_name', 'n_genes', 'spearman_r', 'pvalue']].to_string(index=False))

    return result_df


# ============================================================
# Analysis 3: Promoter vs Gene Body Methylation Function
# ============================================================
def analysis3_promoter_vs_body(data):
    print("\n" + "=" * 60)
    print("Analysis 3: Promoter vs Gene Body Methylation")
    print("=" * 60)

    if not TSS_TABLE.exists() or not HIGH_CONF_SITES.exists():
        print("  TSS table or sites file not found, skipping")
        return None

    tss_df = pd.read_csv(TSS_TABLE)
    sites = pd.read_csv(HIGH_CONF_SITES)

    # Get TSS positions per gene
    tss_info = tss_df[['gene_id', 'chrom', 'tss', 'strand']].dropna(subset=['tss']).drop_duplicates(subset=['gene_id'])
    tss_info['tss'] = tss_info['tss'].astype(int)

    # Assign each methylation site to nearest gene's TSS distance
    gene_site_distances = []
    gene_annot = data['gene_annot']
    gene_coords = gene_annot[['gene_id', 'contig', 'start', 'end', 'strand']].dropna()

    # For each site, find the gene it belongs to and compute TSS distance
    # Use a simpler approach: merge sites with gene info
    for _, gene_row in tss_info.iterrows():
        gene_id = gene_row['gene_id']
        chrom = gene_row['chrom']
        tss_pos = gene_row['tss']
        gene_strand = gene_row['strand']

        gene_sites = sites[sites['chrom'] == chrom].copy()
        if gene_sites.empty:
            continue

        # Get gene boundaries
        gene_info = gene_coords[gene_coords['gene_id'] == gene_id]
        if gene_info.empty:
            continue
        g_start = gene_info.iloc[0]['start']
        g_end = gene_info.iloc[0]['end']

        # Sites within extended gene region (TSS-500 to gene end + 100)
        if gene_strand == '+':
            region_start = tss_pos - 500
            region_end = g_end + 100
        else:
            region_start = g_start - 100
            region_end = tss_pos + 500

        nearby = gene_sites[(gene_sites['position'] >= region_start) & (gene_sites['position'] <= region_end)]
        for _, site in nearby.iterrows():
            if gene_strand == '+':
                dist = site['position'] - tss_pos
            else:
                dist = tss_pos - site['position']

            gene_site_distances.append({
                'gene_id': gene_id,
                'site_position': site['position'],
                'mod_type': site['mod_type'],
                'timepoint': site['timepoint'],
                'tss_distance': dist,
                'region': 'promoter' if dist < 0 else ('tss_proximal' if dist <= 200 else 'gene_body')
            })

    if not gene_site_distances:
        print("  No site-gene assignments found")
        return None

    dist_df = pd.DataFrame(gene_site_distances)
    print(f"  Total site-gene assignments: {len(dist_df)}")
    print(f"  Promoter: {(dist_df['region'] == 'promoter').sum()}, TSS proximal: {(dist_df['region'] == 'tss_proximal').sum()}, Gene body: {(dist_df['region'] == 'gene_body').sum()}")

    # Get unique genes per region
    promoter_genes = dist_df[dist_df['region'] == 'promoter']['gene_id'].unique().tolist()
    body_genes = dist_df[dist_df['region'] == 'gene_body']['gene_id'].unique().tolist()
    proximal_genes = dist_df[dist_df['region'] == 'tss_proximal']['gene_id'].unique().tolist()

    print(f"  Promoter genes: {len(promoter_genes)}, Gene body genes: {len(body_genes)}, TSS proximal genes: {len(proximal_genes)}")

    # COG enrichment for each region
    region_results = {}
    for region_name, genes in [('promoter', promoter_genes), ('tss_proximal', proximal_genes), ('gene_body', body_genes)]:
        genes_in_bg = [g for g in genes if g in data['cog_bg']]
        if len(genes_in_bg) < MIN_GENES:
            continue
        cog_res = run_ora(genes_in_bg, data['cog_bg'], data['cog_to_genes'], COG_NAMES, min_genes=1)
        region_results[region_name] = cog_res
        n_sig = len(cog_res[cog_res['padj'] < 0.05]) if not cog_res.empty else 0
        print(f"  {region_name} COG: {len(genes_in_bg)} genes -> {n_sig} sig categories")
        if not cog_res.empty:
            cog_res.to_csv(ANALYSIS_DIR / "tables" / f"A3_COG_{region_name}.tsv", sep='\t', index=False)

    # KEGG enrichment
    for region_name, genes in [('promoter', promoter_genes), ('tss_proximal', proximal_genes), ('gene_body', body_genes)]:
        sco_genes = [data['id_map'][g] for g in genes if g in data['id_map'] and data['id_map'][g] in data['all_sco_bg']]
        if len(sco_genes) < MIN_GENES:
            continue
        kegg_res = run_ora(sco_genes, data['all_sco_bg'], data['pathway_to_genes'], data['pathway_names'])
        n_sig = len(kegg_res[kegg_res['padj'] < 0.05]) if not kegg_res.empty else 0
        print(f"  {region_name} KEGG: {len(sco_genes)} genes -> {n_sig} sig pathways")
        if not kegg_res.empty:
            kegg_res.to_csv(ANALYSIS_DIR / "tables" / f"A3_KEGG_{region_name}.tsv", sep='\t', index=False)

    # Plot COG distribution comparison
    cog_copy = data['cog'].dropna(subset=['COG_letter'])
    all_cog_genes = set(cog_copy['gene_id'].unique())
    categories = sorted([c for c in COG_NAMES.keys() if c != 'S'])

    fig, ax = plt.subplots(figsize=(14, 6))
    prom_set = set(promoter_genes) & all_cog_genes
    body_set = set(body_genes) & all_cog_genes
    prox_set = set(proximal_genes) & all_cog_genes

    prom_fracs, body_fracs, prox_fracs, bg_fracs = [], [], [], []
    for cat in categories:
        cat_genes = set(cog_copy[cog_copy['COG_letter'] == cat]['gene_id'])
        prom_fracs.append(len(prom_set & cat_genes) / len(prom_set) * 100 if prom_set else 0)
        body_fracs.append(len(body_set & cat_genes) / len(body_set) * 100 if body_set else 0)
        prox_fracs.append(len(prox_set & cat_genes) / len(prox_set) * 100 if prox_set else 0)
        bg_fracs.append(len(all_cog_genes & cat_genes) / len(all_cog_genes) * 100)

    x = np.arange(len(categories))
    w = 0.2
    ax.bar(x - 1.5 * w, bg_fracs, w, label='Background', color='#999999', alpha=0.6)
    ax.bar(x - 0.5 * w, prom_fracs, w, label='Promoter (<TSS)', color='#d73027', alpha=0.8)
    ax.bar(x + 0.5 * w, prox_fracs, w, label='TSS proximal (0-200bp)', color='#fc8d59', alpha=0.8)
    ax.bar(x + 1.5 * w, body_fracs, w, label='Gene body (>200bp)', color='#4575b4', alpha=0.8)

    ax.set_xticks(x)
    ax.set_xticklabels([f"{c}\n{COG_NAMES[c][:12]}" for c in categories], fontsize=7, rotation=45, ha='right')
    ax.set_ylabel('Proportion (%)')
    ax.set_title('COG Distribution: Promoter vs TSS Proximal vs Gene Body Methylation')
    ax.legend(fontsize=8)
    fig.tight_layout()
    save_figure(fig, 'A3_COG_promoter_vs_body')

    return region_results


# ============================================================
# Analysis 4: MTase / Isomerase Expression Verification
# ============================================================
def analysis4_expression_verification(data):
    print("\n" + "=" * 60)
    print("Analysis 4: MTase / Isomerase Expression Verification")
    print("=" * 60)

    gm = data['gene_master']

    # V-Defense MTase genes (from previous analysis)
    defense_genes = [
        'SC_RS05450', 'SC_RS06105', 'SC_RS06595', 'SC_RS09775', 'SC_RS11875',
        'SC_RS12500', 'SC_RS12890', 'SC_RS13205', 'SC_RS13305', 'SC_RS15340',
        'SC_RS19420', 'SC_RS38275', 'SC_RS38280', 'SC_RS39465', 'SC_RS40100'
    ]

    # GO isomerase genes
    isomerase_genes = [
        'SC_RS04810', 'SC_RS04820', 'SC_RS11045', 'SC_RS15835',
        'SC_RS15855', 'SC_RS18320', 'SC_RS18855', 'SC_RS35070'
    ]

    def check_genes(gene_list, label):
        subset = gm[gm['gene_id'].isin(gene_list)].copy()
        results = []
        for _, row in subset.iterrows():
            r = {'gene_id': row['gene_id']}
            if 'old_locus_tag' in row:
                r['old_locus_tag'] = row['old_locus_tag']
            if 'product' in row:
                r['product'] = row['product']

            for comp, lfc_col, padj_col in [
                ('T2vsT1', 'log2FoldChange_2_vs_1', 'padj_2_vs_1'),
                ('T3vsT1', 'log2FoldChange_3_vs_1', 'padj_3_vs_1'),
                ('T3vsT2', 'log2FoldChange_3_vs_2', 'padj_3_vs_2'),
            ]:
                lfc = row.get(lfc_col, np.nan)
                padj = row.get(padj_col, np.nan)
                is_deg = (not pd.isna(padj)) and (padj < DEG_PADJ) and (not pd.isna(lfc)) and (abs(lfc) > DEG_LOG2FC)
                r[f'{comp}_log2FC'] = round(lfc, 3) if not pd.isna(lfc) else np.nan
                r[f'{comp}_padj'] = padj
                r[f'{comp}_DEG'] = is_deg
            results.append(r)

        result_df = pd.DataFrame(results)
        result_df.to_csv(ANALYSIS_DIR / "tables" / f"A4_{label}_expression.tsv", sep='\t', index=False)

        # Summary
        for comp in ['T2vsT1', 'T3vsT1', 'T3vsT2']:
            n_deg = result_df[f'{comp}_DEG'].sum()
            n_total = len(result_df)
            print(f"    {label} {comp}: {n_deg}/{n_total} are DEGs")

        return result_df

    print("\n  --- V-Defense MTase Genes ---")
    defense_df = check_genes(defense_genes, 'defense_MTase')

    print("\n  --- GO Isomerase Genes ---")
    iso_df = check_genes(isomerase_genes, 'GO_isomerase')

    # Also check methylation data
    methyl = data['methyl']
    print("\n  --- Methylation changes for Defense MTase genes ---")
    defense_methyl = methyl[methyl['gene_id'].isin(defense_genes)][
        ['gene_id', '4mC_change_T2_vs_T1', '4mC_change_T3_vs_T1', '4mC_change_T3_vs_T2',
         '6mA_change_T2_vs_T1', '6mA_change_T3_vs_T1', '6mA_change_T3_vs_T2']
    ]
    defense_methyl.to_csv(ANALYSIS_DIR / "tables" / "A4_defense_MTase_methylation.tsv", sep='\t', index=False)
    print(defense_methyl.to_string(index=False))

    print("\n  --- Methylation changes for GO Isomerase genes ---")
    iso_methyl = methyl[methyl['gene_id'].isin(isomerase_genes)][
        ['gene_id', '4mC_change_T2_vs_T1', '4mC_change_T3_vs_T1', '4mC_change_T3_vs_T2',
         '6mA_change_T2_vs_T1', '6mA_change_T3_vs_T1', '6mA_change_T3_vs_T2']
    ]
    iso_methyl.to_csv(ANALYSIS_DIR / "tables" / "A4_GO_isomerase_methylation.tsv", sep='\t', index=False)
    print(iso_methyl.to_string(index=False))

    return defense_df, iso_df


# ============================================================
# Analysis 5: TF / BGC DMG Frequency Test
# ============================================================
def analysis5_tf_bgc_dmg_frequency(data):
    print("\n" + "=" * 60)
    print("Analysis 5: TF / BGC DMG Frequency Test")
    print("=" * 60)

    methyl = data['methyl']

    # Load TF/regulator info
    if GENE_MASTER_BGC_REG.exists():
        reg_df = pd.read_csv(GENE_MASTER_BGC_REG, sep='\t')
        tf_genes = set(reg_df[reg_df['is_regulator'] == True]['gene_id'].tolist())
        print(f"  TF/Regulator genes: {len(tf_genes)}")
    else:
        # Fallback: identify TFs from product name
        gm = data['gene_master']
        tf_mask = gm['product'].str.contains('regulator|transcription|DNA-binding', case=False, na=False)
        tf_genes = set(gm[tf_mask]['gene_id'].tolist())
        print(f"  TF genes (from product name): {len(tf_genes)}")

    # Load BGC gene info
    if BGC_DEF.exists():
        bgc_df = pd.read_csv(BGC_DEF, sep='\t')
        # Map BGC gene_ids through gene_master to get current IDs
        gm = data['gene_master']
        # BGC_definition uses gene_id which may be old annotation
        # Cross-reference via old_locus_tag
        bgc_locus = set(bgc_df['old_locus_tag'].dropna().tolist())
        bgc_genes_via_locus = set(gm[gm['old_locus_tag'].isin(bgc_locus)]['gene_id'].tolist())
        # Also try direct gene_id match
        bgc_genes_direct = set(bgc_df['gene_id'].dropna().tolist()) & set(methyl['gene_id'].tolist())
        bgc_genes = bgc_genes_via_locus | bgc_genes_direct
        print(f"  BGC genes: {len(bgc_genes)}")
    else:
        bgc_genes = set()
        print("  BGC definition file not found")

    # Compute DMG status for all genes
    comparisons = {
        'T2vsT1': ('6mA_change_T2_vs_T1', '4mC_change_T2_vs_T1'),
        'T3vsT1': ('6mA_change_T3_vs_T1', '4mC_change_T3_vs_T1'),
        'T3vsT2': ('6mA_change_T3_vs_T2', '4mC_change_T3_vs_T2'),
    }

    all_genes = set(methyl['gene_id'].tolist())
    results = []

    for comp, (col_6mA, col_4mC) in comparisons.items():
        c6 = methyl[col_6mA].fillna(0)
        c4 = methyl[col_4mC].fillna(0)
        dmg_mask = (c6.abs() > DMG_CHANGE_THRESHOLD) | (c4.abs() > DMG_CHANGE_THRESHOLD)
        dmg_genes = set(methyl[dmg_mask]['gene_id'].tolist())
        non_dmg_genes = all_genes - dmg_genes

        for class_name, class_genes in [('TF', tf_genes), ('BGC', bgc_genes)]:
            class_in_all = class_genes & all_genes
            non_class = all_genes - class_in_all

            if not class_in_all:
                continue

            # 2x2: DMG vs non-DMG × class vs non-class
            a = len(class_in_all & dmg_genes)
            b = len(non_class & dmg_genes)
            c = len(class_in_all & non_dmg_genes)
            d = len(non_class & non_dmg_genes)

            odds_ratio, p_value = stats.fisher_exact([[a, b], [c, d]])

            dmg_rate_class = a / len(class_in_all) * 100 if class_in_all else 0
            dmg_rate_non = b / len(non_class) * 100 if non_class else 0

            results.append({
                'comparison': comp, 'gene_class': class_name,
                'n_class': len(class_in_all), 'n_DMG_in_class': a,
                'DMG_rate_class': round(dmg_rate_class, 1),
                'n_non_class': len(non_class), 'n_DMG_in_non_class': b,
                'DMG_rate_non_class': round(dmg_rate_non, 1),
                'odds_ratio': round(odds_ratio, 3), 'pvalue': p_value,
                'significant': p_value < 0.05
            })

            print(f"  {comp} {class_name}: DMG rate = {dmg_rate_class:.1f}% ({a}/{len(class_in_all)}) vs non-{class_name}: {dmg_rate_non:.1f}% ({b}/{len(non_class)}), OR={odds_ratio:.3f}, p={p_value:.4f}")

    result_df = pd.DataFrame(results)
    result_df.to_csv(ANALYSIS_DIR / "tables" / "A5_TF_BGC_DMG_frequency.tsv", sep='\t', index=False)

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for idx, class_name in enumerate(['TF', 'BGC']):
        ax = axes[idx]
        subset = result_df[result_df['gene_class'] == class_name]
        if subset.empty:
            continue

        x = np.arange(len(subset))
        w = 0.35
        ax.bar(x - w / 2, subset['DMG_rate_class'], w, label=class_name, color='#d73027', alpha=0.8)
        ax.bar(x + w / 2, subset['DMG_rate_non_class'], w, label=f'non-{class_name}', color='#4575b4', alpha=0.8)

        ax.set_xticks(x)
        ax.set_xticklabels(subset['comparison'])
        ax.set_ylabel('DMG Rate (%)')
        ax.set_title(f'DMG Frequency: {class_name} vs non-{class_name}')
        ax.legend()

        for i, row in enumerate(subset.itertuples()):
            sig = '*' if row.pvalue < 0.05 else 'ns'
            y_max = max(row.DMG_rate_class, row.DMG_rate_non_class)
            ax.text(i, y_max + 1, f'OR={row.odds_ratio:.2f}\n{sig}', ha='center', fontsize=8)

    fig.tight_layout()
    save_figure(fig, 'A5_TF_BGC_DMG_frequency')

    return result_df


# ============================================================
# Main
# ============================================================
def main():
    print("=" * 60)
    print("DMG Selectivity Deep Analysis (5 analyses)")
    print(f"Output: {ANALYSIS_DIR}")
    print("=" * 60)

    data = load_all_data()

    # Analysis 1
    a1_results, a1_summary = analysis1_coordinated_enrichment(data)

    # Analysis 2
    a2_results = analysis2_category_correlation(data)

    # Analysis 3
    a3_results = analysis3_promoter_vs_body(data)

    # Analysis 4
    a4_defense, a4_iso = analysis4_expression_verification(data)

    # Analysis 5
    a5_results = analysis5_tf_bgc_dmg_frequency(data)

    print("\n" + "=" * 60)
    print("ALL ANALYSES COMPLETE")
    print(f"Output: {ANALYSIS_DIR}")
    print("=" * 60)


if __name__ == '__main__':
    main()
