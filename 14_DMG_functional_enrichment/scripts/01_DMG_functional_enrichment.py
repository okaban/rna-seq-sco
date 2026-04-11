#!/usr/bin/env python3
"""
DMG Functional Enrichment Analysis
===================================
DMG (Differentially Methylated Genes) に対する機能エンリッチメント解析
- KEGG pathway enrichment (via KEGG REST API)
- GO term enrichment
- COG category enrichment
- DEG enrichment との比較
- モチーフ別 (CCGG / AAGCCCG) 機能分類

Usage:
    conda activate rna-seq
    python 01_DMG_functional_enrichment.py
"""

import os
import sys
import time
import warnings
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
import requests
from collections import defaultdict
from pathlib import Path

warnings.filterwarnings('ignore')

# ============================================================
# Configuration
# ============================================================
PROJECT_ROOT = Path("/Users/okaban/bioinfo/rna-seq")
ANALYSIS_DIR = PROJECT_ROOT / "14_DMG_functional_enrichment" / "analysis" / "14_DMG_enrichment_260206_v1"
ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
(ANALYSIS_DIR / "figures").mkdir(exist_ok=True)
(ANALYSIS_DIR / "tables").mkdir(exist_ok=True)

# Input files
INTEGRATED_METHYL = PROJECT_ROOT / "11_epigenome_integration" / "analysis" / "01_integration" / "integrated_methyl_expression_weighted.csv"
GENE_MASTER = PROJECT_ROOT / "05_annotation" / "analysis" / "05_annotation_260128_v1" / "tables" / "gene_master_DESeq2.tsv"
GENE_ANNOT = PROJECT_ROOT / "05_annotation" / "analysis" / "05_annotation_260128_v1" / "tables" / "gene_annotation_basic.tsv"
COG_CLASS = PROJECT_ROOT / "12_supplementary_figures" / "analysis" / "12_supplementary_260202_v1" / "tables" / "gene_COG_classification.tsv"
MOTIF_CCGG = PROJECT_ROOT / "11_epigenome_integration" / "analysis" / "18_tss_analyses" / "motif_CCGG_spatial_detail.csv"
MOTIF_AAGCCCG = PROJECT_ROOT / "11_epigenome_integration" / "analysis" / "18_tss_analyses" / "motif_AAGCCCG_spatial_detail.csv"

# Existing DEG KEGG enrichment for comparison
DEG_KEGG_DIR = PROJECT_ROOT / "12_supplementary_figures" / "analysis" / "12_supplementary_260202_v1" / "tables"

# Thresholds
DMG_CHANGE_THRESHOLD = 10.0  # 10% methylation change (same as existing)
DEG_PADJ_THRESHOLD = 0.05
DEG_LOG2FC_THRESHOLD = 1.0
ENRICHMENT_PADJ_THRESHOLD = 0.05
MIN_GENES_PER_TERM = 3

# Comparisons
COMPARISONS = ['T2vsT1', 'T3vsT1', 'T3vsT2']

# Plot settings
plt.rcParams.update({
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'font.family': 'sans-serif',
})


# ============================================================
# 1. Data Loading
# ============================================================
def load_data():
    """Load all required data files."""
    print("=" * 60)
    print("Loading data...")
    print("=" * 60)

    methyl = pd.read_csv(INTEGRATED_METHYL)
    print(f"  Integrated methylation data: {len(methyl)} genes")

    gene_master = pd.read_csv(GENE_MASTER, sep='\t')
    print(f"  Gene master: {len(gene_master)} genes")

    gene_annot = pd.read_csv(GENE_ANNOT, sep='\t')
    print(f"  Gene annotation: {len(gene_annot)} genes")

    cog = pd.read_csv(COG_CLASS, sep='\t')
    print(f"  COG classification: {len(cog)} genes")

    return methyl, gene_master, gene_annot, cog


# ============================================================
# 2. DMG Extraction
# ============================================================
def extract_dmgs(methyl_df):
    """Extract DMGs stratified by comparison, modification type, and direction."""
    print("\n" + "=" * 60)
    print("Extracting DMGs...")
    print("=" * 60)

    col_map = {
        'T2vsT1': {
            '6mA': '6mA_change_T2_vs_T1',
            '4mC': '4mC_change_T2_vs_T1'
        },
        'T3vsT1': {
            '6mA': '6mA_change_T3_vs_T1',
            '4mC': '4mC_change_T3_vs_T1'
        },
        'T3vsT2': {
            '6mA': '6mA_change_T3_vs_T2',
            '4mC': '4mC_change_T3_vs_T2'
        }
    }

    dmg_sets = {}
    summary_rows = []

    for comp in COMPARISONS:
        for mod_type in ['6mA', '4mC']:
            col = col_map[comp][mod_type]
            change = methyl_df[col].fillna(0)

            hyper = methyl_df.loc[change > DMG_CHANGE_THRESHOLD, 'gene_id'].tolist()
            hypo = methyl_df.loc[change < -DMG_CHANGE_THRESHOLD, 'gene_id'].tolist()
            all_dmg = list(set(hyper + hypo))

            dmg_sets[f"{comp}_{mod_type}_hyper"] = hyper
            dmg_sets[f"{comp}_{mod_type}_hypo"] = hypo
            dmg_sets[f"{comp}_{mod_type}_all"] = all_dmg

            summary_rows.append({
                'comparison': comp,
                'mod_type': mod_type,
                'hyper': len(hyper),
                'hypo': len(hypo),
                'total_DMG': len(all_dmg)
            })

        # Combined (either 6mA or 4mC)
        col_6mA = col_map[comp]['6mA']
        col_4mC = col_map[comp]['4mC']
        change_6mA = methyl_df[col_6mA].fillna(0)
        change_4mC = methyl_df[col_4mC].fillna(0)

        hyper_combined = methyl_df.loc[
            (change_6mA > DMG_CHANGE_THRESHOLD) | (change_4mC > DMG_CHANGE_THRESHOLD),
            'gene_id'
        ].tolist()
        hypo_combined = methyl_df.loc[
            (change_6mA < -DMG_CHANGE_THRESHOLD) | (change_4mC < -DMG_CHANGE_THRESHOLD),
            'gene_id'
        ].tolist()
        all_combined = methyl_df.loc[
            (change_6mA.abs() > DMG_CHANGE_THRESHOLD) | (change_4mC.abs() > DMG_CHANGE_THRESHOLD),
            'gene_id'
        ].tolist()

        dmg_sets[f"{comp}_combined_hyper"] = hyper_combined
        dmg_sets[f"{comp}_combined_hypo"] = hypo_combined
        dmg_sets[f"{comp}_combined_all"] = all_combined

        summary_rows.append({
            'comparison': comp,
            'mod_type': 'combined',
            'hyper': len(hyper_combined),
            'hypo': len(hypo_combined),
            'total_DMG': len(all_combined)
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(ANALYSIS_DIR / "tables" / "DMG_summary.tsv", sep='\t', index=False)
    print("\nDMG Summary:")
    print(summary_df.to_string(index=False))

    return dmg_sets, summary_df


# ============================================================
# 3. DEG Extraction
# ============================================================
def extract_degs(gene_master):
    """Extract DEGs for comparison with DMG enrichment."""
    print("\n" + "=" * 60)
    print("Extracting DEGs...")
    print("=" * 60)

    deg_col_map = {
        'T2vsT1': ('log2FoldChange_2_vs_1', 'padj_2_vs_1'),
        'T3vsT1': ('log2FoldChange_3_vs_1', 'padj_3_vs_1'),
        'T3vsT2': ('log2FoldChange_3_vs_2', 'padj_3_vs_2'),
    }

    deg_sets = {}
    for comp, (lfc_col, padj_col) in deg_col_map.items():
        mask_up = (gene_master[padj_col] < DEG_PADJ_THRESHOLD) & (gene_master[lfc_col] > DEG_LOG2FC_THRESHOLD)
        mask_down = (gene_master[padj_col] < DEG_PADJ_THRESHOLD) & (gene_master[lfc_col] < -DEG_LOG2FC_THRESHOLD)

        deg_sets[f"{comp}_up"] = gene_master.loc[mask_up, 'gene_id'].tolist()
        deg_sets[f"{comp}_down"] = gene_master.loc[mask_down, 'gene_id'].tolist()
        deg_sets[f"{comp}_all"] = gene_master.loc[mask_up | mask_down, 'gene_id'].tolist()

        print(f"  {comp}: Up={len(deg_sets[f'{comp}_up'])}, Down={len(deg_sets[f'{comp}_down'])}, Total={len(deg_sets[f'{comp}_all'])}")

    return deg_sets


# ============================================================
# 4. KEGG Pathway Enrichment
# ============================================================
def fetch_kegg_data():
    """Fetch KEGG pathway data from REST API."""
    print("\n  Fetching KEGG pathway data from REST API...")

    cache_file = ANALYSIS_DIR / "tables" / "kegg_cache.tsv"
    if cache_file.exists():
        print("  Using cached KEGG data.")
        kegg_df = pd.read_csv(cache_file, sep='\t')
        gene_to_pathways = defaultdict(set)
        pathway_names = {}
        for _, row in kegg_df.iterrows():
            gene_to_pathways[row['gene']].add(row['pathway_id'])
            pathway_names[row['pathway_id']] = row['pathway_name']
        return gene_to_pathways, pathway_names

    # Fetch pathway list
    url_pathways = "https://rest.kegg.jp/list/pathway/sco"
    resp = requests.get(url_pathways)
    resp.raise_for_status()
    pathway_names = {}
    for line in resp.text.strip().split('\n'):
        parts = line.split('\t')
        pw_id = parts[0].replace('path:', '')
        pw_name = parts[1].split(' - ')[0].strip()
        pathway_names[pw_id] = pw_name

    print(f"  Fetched {len(pathway_names)} KEGG pathways")
    time.sleep(1)

    # Fetch gene-pathway links
    url_links = "https://rest.kegg.jp/link/sco/pathway"
    resp = requests.get(url_links)
    resp.raise_for_status()
    gene_to_pathways = defaultdict(set)
    cache_rows = []
    for line in resp.text.strip().split('\n'):
        parts = line.split('\t')
        pw_id = parts[0].replace('path:', '')
        gene_sco = parts[1].replace('sco:', '')
        gene_to_pathways[gene_sco].add(pw_id)
        cache_rows.append({'pathway_id': pw_id, 'pathway_name': pathway_names.get(pw_id, ''), 'gene': gene_sco})

    print(f"  Fetched gene-pathway links for {len(gene_to_pathways)} genes")

    # Cache
    pd.DataFrame(cache_rows).to_csv(cache_file, sep='\t', index=False)

    return gene_to_pathways, pathway_names


def run_ora(gene_list, background, term_to_genes, term_names, min_genes=MIN_GENES_PER_TERM):
    """
    Run Over-Representation Analysis with Fisher's exact test.
    Returns DataFrame with enrichment results.
    """
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

        # 2x2 contingency table
        a = n_overlap
        b = n_list - n_overlap
        c = n_term - n_overlap
        d = n_total - n_list - n_term + n_overlap

        odds_ratio, p_value = stats.fisher_exact([[a, b], [c, d]], alternative='greater')

        gene_ratio = n_overlap / n_list if n_list > 0 else 0
        bg_ratio = n_term / n_total if n_total > 0 else 0
        fold_enrichment = gene_ratio / bg_ratio if bg_ratio > 0 else 0

        results.append({
            'term_id': term_id,
            'term_name': term_names.get(term_id, term_id),
            'count': n_overlap,
            'total_in_term': n_term,
            'gene_ratio': f"{n_overlap}/{n_list}",
            'fold_enrichment': round(fold_enrichment, 3),
            'pvalue': p_value,
            'genes': ';'.join(sorted(overlap))
        })

    if not results:
        return pd.DataFrame()

    result_df = pd.DataFrame(results)
    _, padj, _, _ = multipletests(result_df['pvalue'].values, method='fdr_bh')
    result_df['padj'] = padj
    result_df = result_df.sort_values('pvalue').reset_index(drop=True)

    return result_df


def run_kegg_enrichment(dmg_sets, deg_sets, gene_annot):
    """Run KEGG pathway enrichment for DMGs and DEGs."""
    print("\n" + "=" * 60)
    print("Running KEGG pathway enrichment...")
    print("=" * 60)

    gene_to_pathways, pathway_names = fetch_kegg_data()

    # Build pathway_to_genes (inverted)
    pathway_to_genes = defaultdict(set)
    for gene, pws in gene_to_pathways.items():
        for pw in pws:
            pathway_to_genes[pw].add(gene)

    # Map gene_id -> old_locus_tag (SCO####) for KEGG
    id_map = gene_annot[['gene_id', 'old_locus_tag']].dropna().set_index('gene_id')['old_locus_tag'].to_dict()
    # Background: all genes with SCO IDs that are in KEGG
    all_sco = set(id_map.values()) & set(gene_to_pathways.keys())
    print(f"  Background genes with KEGG annotation: {len(all_sco)}")

    def convert_to_sco(gene_list):
        return [id_map[g] for g in gene_list if g in id_map and id_map[g] in all_sco]

    kegg_results = {}

    # DMG enrichment
    for comp in COMPARISONS:
        for mod_type in ['6mA', '4mC', 'combined']:
            for direction in ['hyper', 'hypo', 'all']:
                key = f"{comp}_{mod_type}_{direction}"
                genes = dmg_sets.get(key, [])
                sco_genes = convert_to_sco(genes)
                if len(sco_genes) < MIN_GENES_PER_TERM:
                    continue

                result = run_ora(sco_genes, all_sco, pathway_to_genes, pathway_names)
                kegg_results[f"DMG_{key}"] = result
                if not result.empty:
                    sig = result[result['padj'] < ENRICHMENT_PADJ_THRESHOLD]
                    print(f"  DMG {key}: {len(sco_genes)} genes -> {len(sig)} significant pathways (padj<0.05)")

    # DEG enrichment (for comparison)
    for comp in COMPARISONS:
        for direction in ['up', 'down', 'all']:
            key = f"{comp}_{direction}"
            genes = deg_sets.get(key, [])
            sco_genes = convert_to_sco(genes)
            if len(sco_genes) < MIN_GENES_PER_TERM:
                continue

            result = run_ora(sco_genes, all_sco, pathway_to_genes, pathway_names)
            kegg_results[f"DEG_{key}"] = result
            if not result.empty:
                sig = result[result['padj'] < ENRICHMENT_PADJ_THRESHOLD]
                print(f"  DEG {key}: {len(sco_genes)} genes -> {len(sig)} significant pathways (padj<0.05)")

    # Save all results
    for name, df in kegg_results.items():
        if not df.empty:
            df.to_csv(ANALYSIS_DIR / "tables" / f"KEGG_{name}.tsv", sep='\t', index=False)

    return kegg_results


# ============================================================
# 5. GO Term Enrichment
# ============================================================
def run_go_enrichment(dmg_sets, deg_sets, gene_annot):
    """Run GO term enrichment for DMGs."""
    print("\n" + "=" * 60)
    print("Running GO term enrichment...")
    print("=" * 60)

    # Build gene-to-GO mapping from gene_annotation_basic.tsv
    go_mapping = gene_annot[['gene_id', 'ontology_term']].dropna(subset=['ontology_term'])
    go_mapping = go_mapping[go_mapping['ontology_term'] != '']

    gene_to_go = {}
    go_to_genes = defaultdict(set)
    for _, row in go_mapping.iterrows():
        gene = row['gene_id']
        terms = [t.strip() for t in str(row['ontology_term']).split(',') if t.strip().startswith('GO:')]
        gene_to_go[gene] = set(terms)
        for t in terms:
            go_to_genes[t].add(gene)

    background = set(gene_to_go.keys())
    print(f"  Genes with GO annotation: {len(background)}")
    print(f"  GO terms: {len(go_to_genes)}")

    go_names = {t: t for t in go_to_genes}  # Use GO ID as name (no obo file)

    go_results = {}

    # DMG enrichment
    for comp in COMPARISONS:
        for mod_type in ['6mA', '4mC', 'combined']:
            for direction in ['hyper', 'hypo', 'all']:
                key = f"{comp}_{mod_type}_{direction}"
                genes = dmg_sets.get(key, [])
                genes_in_bg = [g for g in genes if g in background]
                if len(genes_in_bg) < MIN_GENES_PER_TERM:
                    continue

                result = run_ora(genes_in_bg, background, go_to_genes, go_names)
                go_results[f"DMG_{key}"] = result
                if not result.empty:
                    sig = result[result['padj'] < ENRICHMENT_PADJ_THRESHOLD]
                    print(f"  DMG {key}: {len(genes_in_bg)} genes -> {len(sig)} significant GO terms (padj<0.05)")

    # DEG enrichment
    for comp in COMPARISONS:
        for direction in ['up', 'down', 'all']:
            key = f"{comp}_{direction}"
            genes = deg_sets.get(key, [])
            genes_in_bg = [g for g in genes if g in background]
            if len(genes_in_bg) < MIN_GENES_PER_TERM:
                continue

            result = run_ora(genes_in_bg, background, go_to_genes, go_names)
            go_results[f"DEG_{key}"] = result
            if not result.empty:
                sig = result[result['padj'] < ENRICHMENT_PADJ_THRESHOLD]
                print(f"  DEG {key}: {len(genes_in_bg)} genes -> {len(sig)} significant GO terms (padj<0.05)")

    for name, df in go_results.items():
        if not df.empty:
            df.to_csv(ANALYSIS_DIR / "tables" / f"GO_{name}.tsv", sep='\t', index=False)

    return go_results


# ============================================================
# 6. COG Category Enrichment
# ============================================================
COG_CATEGORY_NAMES = {
    'C': 'Energy production',
    'D': 'Cell division',
    'E': 'Amino acid metabolism',
    'F': 'Nucleotide metabolism',
    'G': 'Carbohydrate metabolism',
    'H': 'Coenzyme metabolism',
    'I': 'Lipid metabolism',
    'J': 'Translation',
    'K': 'Transcription',
    'L': 'Replication/Repair',
    'M': 'Cell wall',
    'N': 'Motility',
    'O': 'Chaperones/PTM',
    'P': 'Inorganic ion transport',
    'Q': 'Secondary metabolism',
    'R': 'General function',
    'S': 'Unknown',
    'T': 'Signal transduction',
    'U': 'Secretion',
    'V': 'Defense',
}


def run_cog_enrichment(dmg_sets, deg_sets, cog_df):
    """Run COG category enrichment for DMGs."""
    print("\n" + "=" * 60)
    print("Running COG category enrichment...")
    print("=" * 60)

    # Parse COG categories (extract single letter code)
    cog_df = cog_df.copy()
    cog_df['COG_letter'] = cog_df['COG_category'].str.extract(r'^([A-Z])')
    cog_df = cog_df.dropna(subset=['COG_letter'])

    # Build mapping
    cog_to_genes = defaultdict(set)
    for _, row in cog_df.iterrows():
        cog_to_genes[row['COG_letter']].add(row['gene_id'])

    background = set(cog_df['gene_id'].unique())
    print(f"  Genes with COG annotation: {len(background)}")

    cog_results = {}

    # DMG enrichment
    for comp in COMPARISONS:
        for mod_type in ['6mA', '4mC', 'combined']:
            for direction in ['hyper', 'hypo', 'all']:
                key = f"{comp}_{mod_type}_{direction}"
                genes = dmg_sets.get(key, [])
                genes_in_bg = [g for g in genes if g in background]
                if len(genes_in_bg) < MIN_GENES_PER_TERM:
                    continue

                result = run_ora(genes_in_bg, background, cog_to_genes, COG_CATEGORY_NAMES, min_genes=1)
                cog_results[f"DMG_{key}"] = result
                if not result.empty:
                    sig = result[result['padj'] < ENRICHMENT_PADJ_THRESHOLD]
                    print(f"  DMG {key}: {len(genes_in_bg)} genes -> {len(sig)} significant COG categories (padj<0.05)")

    # DEG enrichment
    for comp in COMPARISONS:
        for direction in ['up', 'down', 'all']:
            key = f"{comp}_{direction}"
            genes = deg_sets.get(key, [])
            genes_in_bg = [g for g in genes if g in background]
            if len(genes_in_bg) < MIN_GENES_PER_TERM:
                continue

            result = run_ora(genes_in_bg, background, cog_to_genes, COG_CATEGORY_NAMES, min_genes=1)
            cog_results[f"DEG_{key}"] = result
            if not result.empty:
                sig = result[result['padj'] < ENRICHMENT_PADJ_THRESHOLD]
                print(f"  DEG {key}: {len(genes_in_bg)} genes -> {len(sig)} significant COG categories (padj<0.05)")

    for name, df in cog_results.items():
        if not df.empty:
            df.to_csv(ANALYSIS_DIR / "tables" / f"COG_{name}.tsv", sep='\t', index=False)

    return cog_results


# ============================================================
# 7. Motif-Specific Functional Classification
# ============================================================
def run_motif_classification(cog_df, gene_annot):
    """Classify methylated genes by motif (CCGG vs AAGCCCG) and analyze COG/KEGG."""
    print("\n" + "=" * 60)
    print("Running motif-specific functional classification...")
    print("=" * 60)

    motif_sets = {}

    # CCGG motif genes (methylated)
    if MOTIF_CCGG.exists():
        ccgg = pd.read_csv(MOTIF_CCGG)
        ccgg_methylated = ccgg[ccgg['is_methylated'] == True]['gene_id'].unique().tolist()
        ccgg_all = ccgg['gene_id'].unique().tolist()
        motif_sets['CCGG_methylated'] = ccgg_methylated
        motif_sets['CCGG_all'] = ccgg_all
        print(f"  CCGG: {len(ccgg_methylated)} methylated genes / {len(ccgg_all)} total genes with motif")

    # AAGCCCG motif genes (methylated)
    if MOTIF_AAGCCCG.exists():
        aagcccg = pd.read_csv(MOTIF_AAGCCCG)
        aagcccg_methylated = aagcccg[aagcccg['is_methylated'] == True]['gene_id'].unique().tolist()
        aagcccg_all = aagcccg['gene_id'].unique().tolist()
        motif_sets['AAGCCCG_methylated'] = aagcccg_methylated
        motif_sets['AAGCCCG_all'] = aagcccg_all
        print(f"  AAGCCCG: {len(aagcccg_methylated)} methylated genes / {len(aagcccg_all)} total genes with motif")

    # COG enrichment for motif genes
    cog_copy = cog_df.copy()
    cog_copy['COG_letter'] = cog_copy['COG_category'].str.extract(r'^([A-Z])')
    cog_copy = cog_copy.dropna(subset=['COG_letter'])

    cog_to_genes = defaultdict(set)
    for _, row in cog_copy.iterrows():
        cog_to_genes[row['COG_letter']].add(row['gene_id'])

    background = set(cog_copy['gene_id'].unique())

    motif_cog_results = {}
    for motif_name, genes in motif_sets.items():
        genes_in_bg = [g for g in genes if g in background]
        if len(genes_in_bg) < MIN_GENES_PER_TERM:
            continue
        result = run_ora(genes_in_bg, background, cog_to_genes, COG_CATEGORY_NAMES, min_genes=1)
        motif_cog_results[motif_name] = result
        if not result.empty:
            sig = result[result['padj'] < ENRICHMENT_PADJ_THRESHOLD]
            print(f"  Motif {motif_name}: {len(genes_in_bg)} genes -> {len(sig)} significant COG categories")

    for name, df in motif_cog_results.items():
        if not df.empty:
            df.to_csv(ANALYSIS_DIR / "tables" / f"COG_motif_{name}.tsv", sep='\t', index=False)

    # KEGG enrichment for motif genes
    gene_to_pathways, pathway_names = fetch_kegg_data()
    pathway_to_genes = defaultdict(set)
    for gene, pws in gene_to_pathways.items():
        for pw in pws:
            pathway_to_genes[pw].add(gene)

    id_map = gene_annot[['gene_id', 'old_locus_tag']].dropna().set_index('gene_id')['old_locus_tag'].to_dict()
    all_sco = set(id_map.values()) & set(gene_to_pathways.keys())

    motif_kegg_results = {}
    for motif_name, genes in motif_sets.items():
        sco_genes = [id_map[g] for g in genes if g in id_map and id_map[g] in all_sco]
        if len(sco_genes) < MIN_GENES_PER_TERM:
            continue
        result = run_ora(sco_genes, all_sco, pathway_to_genes, pathway_names)
        motif_kegg_results[motif_name] = result
        if not result.empty:
            sig = result[result['padj'] < ENRICHMENT_PADJ_THRESHOLD]
            print(f"  Motif {motif_name} KEGG: {len(sco_genes)} genes -> {len(sig)} significant pathways")

    for name, df in motif_kegg_results.items():
        if not df.empty:
            df.to_csv(ANALYSIS_DIR / "tables" / f"KEGG_motif_{name}.tsv", sep='\t', index=False)

    return motif_sets, motif_cog_results, motif_kegg_results


# ============================================================
# 8. Figure Generation
# ============================================================
def save_figure(fig, name):
    """Save figure as PDF and SVG."""
    fig.savefig(ANALYSIS_DIR / "figures" / f"{name}.pdf", format='pdf')
    fig.savefig(ANALYSIS_DIR / "figures" / f"{name}.svg", format='svg')
    plt.close(fig)
    print(f"  Saved: {name}.pdf / .svg")


def plot_dmg_summary(dmg_summary):
    """Plot DMG counts by comparison, modification type, and direction."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))

    for idx, comp in enumerate(COMPARISONS):
        ax = axes[idx]
        data = dmg_summary[dmg_summary['comparison'] == comp].copy()
        data = data[data['mod_type'] != 'combined']

        x = np.arange(len(data))
        width = 0.35

        ax.bar(x - width / 2, data['hyper'], width, label='Hyper', color='#d73027', alpha=0.8)
        ax.bar(x + width / 2, data['hypo'], width, label='Hypo', color='#4575b4', alpha=0.8)

        ax.set_xticks(x)
        ax.set_xticklabels(data['mod_type'])
        ax.set_title(comp)
        ax.set_ylabel('Number of DMGs')
        ax.legend()

        for i, row in enumerate(data.itertuples()):
            ax.text(i - width / 2, row.hyper + 2, str(row.hyper), ha='center', va='bottom', fontsize=8)
            ax.text(i + width / 2, row.hypo + 2, str(row.hypo), ha='center', va='bottom', fontsize=8)

    fig.suptitle('DMG Counts by Comparison, Modification Type, and Direction', fontweight='bold')
    fig.tight_layout()
    save_figure(fig, 'DMG_summary_counts')


def plot_kegg_enrichment_top(kegg_results, top_n=15):
    """Plot top KEGG enrichment results for DMGs (combined, all directions)."""
    for comp in COMPARISONS:
        key = f"DMG_{comp}_combined_all"
        if key not in kegg_results or kegg_results[key].empty:
            continue

        df = kegg_results[key].head(top_n).copy()
        df = df.sort_values('fold_enrichment', ascending=True)

        fig, ax = plt.subplots(figsize=(10, max(4, len(df) * 0.35)))
        colors = ['#d73027' if p < 0.05 else '#fee08b' if p < 0.1 else '#cccccc' for p in df['padj']]
        bars = ax.barh(range(len(df)), df['fold_enrichment'], color=colors, edgecolor='gray', linewidth=0.5)

        ax.set_yticks(range(len(df)))
        ax.set_yticklabels(df['term_name'])
        ax.set_xlabel('Fold Enrichment')
        ax.set_title(f'KEGG Pathway Enrichment - DMG {comp}')

        for i, (_, row) in enumerate(df.iterrows()):
            ax.text(row['fold_enrichment'] + 0.05, i, f"n={row['count']}", va='center', fontsize=8)

        # Legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#d73027', label='padj < 0.05'),
            Patch(facecolor='#fee08b', label='padj < 0.1'),
            Patch(facecolor='#cccccc', label='padj >= 0.1'),
        ]
        ax.legend(handles=legend_elements, loc='lower right', fontsize=8)

        fig.tight_layout()
        save_figure(fig, f'KEGG_DMG_{comp}_top{top_n}')


def plot_cog_distribution_comparison(cog_results, dmg_sets, deg_sets, cog_df):
    """Plot COG distribution comparison between DMGs and DEGs."""
    cog_copy = cog_df.copy()
    cog_copy['COG_letter'] = cog_copy['COG_category'].str.extract(r'^([A-Z])')
    cog_copy = cog_copy.dropna(subset=['COG_letter'])

    all_cog_genes = set(cog_copy['gene_id'].unique())

    for comp in COMPARISONS:
        fig, ax = plt.subplots(figsize=(14, 6))

        categories = sorted(COG_CATEGORY_NAMES.keys())
        categories = [c for c in categories if c != 'S']  # Exclude Unknown for cleaner plot

        dmg_genes = set(dmg_sets.get(f"{comp}_combined_all", [])) & all_cog_genes
        deg_genes = set(deg_sets.get(f"{comp}_all", [])) & all_cog_genes
        bg_genes = all_cog_genes

        dmg_fracs = []
        deg_fracs = []
        bg_fracs = []

        for cat in categories:
            cat_genes = set(cog_copy[cog_copy['COG_letter'] == cat]['gene_id'])
            dmg_fracs.append(len(dmg_genes & cat_genes) / len(dmg_genes) * 100 if dmg_genes else 0)
            deg_fracs.append(len(deg_genes & cat_genes) / len(deg_genes) * 100 if deg_genes else 0)
            bg_fracs.append(len(bg_genes & cat_genes) / len(bg_genes) * 100 if bg_genes else 0)

        x = np.arange(len(categories))
        width = 0.25

        ax.bar(x - width, bg_fracs, width, label='Background', color='#999999', alpha=0.6)
        ax.bar(x, deg_fracs, width, label='DEGs', color='#4575b4', alpha=0.8)
        ax.bar(x + width, dmg_fracs, width, label='DMGs', color='#d73027', alpha=0.8)

        ax.set_xticks(x)
        ax.set_xticklabels([f"{c}\n{COG_CATEGORY_NAMES[c][:12]}" for c in categories], fontsize=7, rotation=45, ha='right')
        ax.set_ylabel('Proportion (%)')
        ax.set_title(f'COG Distribution: DMGs vs DEGs vs Background ({comp})')
        ax.legend()

        fig.tight_layout()
        save_figure(fig, f'COG_distribution_{comp}')


def plot_motif_cog_comparison(motif_cog_results, motif_sets, cog_df):
    """Plot COG distribution comparison between CCGG and AAGCCCG methylated genes."""
    cog_copy = cog_df.copy()
    cog_copy['COG_letter'] = cog_copy['COG_category'].str.extract(r'^([A-Z])')
    cog_copy = cog_copy.dropna(subset=['COG_letter'])

    all_cog_genes = set(cog_copy['gene_id'].unique())

    if 'CCGG_methylated' not in motif_sets or 'AAGCCCG_methylated' not in motif_sets:
        print("  Skipping motif COG comparison (data not available)")
        return

    ccgg_genes = set(motif_sets['CCGG_methylated']) & all_cog_genes
    aagcccg_genes = set(motif_sets['AAGCCCG_methylated']) & all_cog_genes

    if not ccgg_genes or not aagcccg_genes:
        print("  Skipping motif COG comparison (no genes in COG)")
        return

    categories = sorted(COG_CATEGORY_NAMES.keys())
    categories = [c for c in categories if c != 'S']

    fig, ax = plt.subplots(figsize=(14, 6))

    ccgg_fracs = []
    aagcccg_fracs = []
    bg_fracs = []

    for cat in categories:
        cat_genes = set(cog_copy[cog_copy['COG_letter'] == cat]['gene_id'])
        ccgg_fracs.append(len(ccgg_genes & cat_genes) / len(ccgg_genes) * 100)
        aagcccg_fracs.append(len(aagcccg_genes & cat_genes) / len(aagcccg_genes) * 100)
        bg_fracs.append(len(all_cog_genes & cat_genes) / len(all_cog_genes) * 100)

    x = np.arange(len(categories))
    width = 0.25

    ax.bar(x - width, bg_fracs, width, label='Background', color='#999999', alpha=0.6)
    ax.bar(x, ccgg_fracs, width, label='CCGG (4mC)', color='#4575b4', alpha=0.8)
    ax.bar(x + width, aagcccg_fracs, width, label='AAGCCCG (6mA)', color='#d73027', alpha=0.8)

    ax.set_xticks(x)
    ax.set_xticklabels([f"{c}\n{COG_CATEGORY_NAMES[c][:12]}" for c in categories], fontsize=7, rotation=45, ha='right')
    ax.set_ylabel('Proportion (%)')
    ax.set_title('COG Distribution: CCGG vs AAGCCCG Methylated Genes')
    ax.legend()

    fig.tight_layout()
    save_figure(fig, 'COG_motif_CCGG_vs_AAGCCCG')


def plot_kegg_dmg_vs_deg_heatmap(kegg_results):
    """Create heatmap comparing DMG vs DEG KEGG enrichment across comparisons."""
    # Collect all significant pathways
    all_pathways = set()
    for key, df in kegg_results.items():
        if df.empty:
            continue
        sig = df[df['padj'] < 0.1]
        all_pathways.update(sig['term_name'].tolist())

    if not all_pathways:
        print("  No significant KEGG pathways for heatmap")
        return

    # Build matrix: rows=pathways, cols=conditions
    conditions_dmg = [f"DMG_{comp}_combined_all" for comp in COMPARISONS]
    conditions_deg = [f"DEG_{comp}_all" for comp in COMPARISONS]
    all_conditions = conditions_dmg + conditions_deg
    col_labels = [f"DMG {c}" for c in COMPARISONS] + [f"DEG {c}" for c in COMPARISONS]

    pathways_sorted = sorted(all_pathways)
    matrix = np.full((len(pathways_sorted), len(all_conditions)), np.nan)

    for j, cond in enumerate(all_conditions):
        if cond in kegg_results and not kegg_results[cond].empty:
            df = kegg_results[cond]
            for i, pw in enumerate(pathways_sorted):
                row = df[df['term_name'] == pw]
                if not row.empty:
                    pval = row.iloc[0]['padj']
                    matrix[i, j] = -np.log10(pval) if pval > 0 else 10

    # Filter to rows with at least one non-NaN
    mask = ~np.all(np.isnan(matrix), axis=1)
    matrix = matrix[mask]
    pathways_sorted = [p for p, m in zip(pathways_sorted, mask) if m]

    if len(pathways_sorted) == 0:
        print("  No data for heatmap")
        return

    # Limit to top 30 pathways by max -log10(padj)
    if len(pathways_sorted) > 30:
        max_vals = np.nanmax(matrix, axis=1)
        top_idx = np.argsort(max_vals)[-30:]
        matrix = matrix[top_idx]
        pathways_sorted = [pathways_sorted[i] for i in top_idx]

    fig, ax = plt.subplots(figsize=(10, max(6, len(pathways_sorted) * 0.35)))
    im = ax.imshow(matrix, aspect='auto', cmap='YlOrRd', interpolation='nearest')

    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=45, ha='right', fontsize=9)
    ax.set_yticks(range(len(pathways_sorted)))
    ax.set_yticklabels(pathways_sorted, fontsize=8)

    # Add significance markers
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if not np.isnan(matrix[i, j]):
                if matrix[i, j] > -np.log10(0.05):
                    ax.text(j, i, '*', ha='center', va='center', fontsize=10, fontweight='bold')

    plt.colorbar(im, ax=ax, label='-log10(padj)', shrink=0.8)
    ax.set_title('KEGG Pathway Enrichment: DMGs vs DEGs', fontweight='bold')
    fig.tight_layout()
    save_figure(fig, 'KEGG_DMG_vs_DEG_heatmap')


def plot_kegg_enrichment_by_direction(kegg_results):
    """Plot KEGG enrichment comparison: hyper vs hypo DMGs."""
    for comp in COMPARISONS:
        key_hyper = f"DMG_{comp}_combined_hyper"
        key_hypo = f"DMG_{comp}_combined_hypo"

        df_hyper = kegg_results.get(key_hyper, pd.DataFrame())
        df_hypo = kegg_results.get(key_hypo, pd.DataFrame())

        if df_hyper.empty and df_hypo.empty:
            continue

        # Collect top pathways from both
        top_pathways = set()
        if not df_hyper.empty:
            top_pathways.update(df_hyper.head(10)['term_name'].tolist())
        if not df_hypo.empty:
            top_pathways.update(df_hypo.head(10)['term_name'].tolist())

        if not top_pathways:
            continue

        pathways = sorted(top_pathways)

        hyper_fe = []
        hypo_fe = []
        for pw in pathways:
            h = df_hyper[df_hyper['term_name'] == pw] if not df_hyper.empty else pd.DataFrame()
            l = df_hypo[df_hypo['term_name'] == pw] if not df_hypo.empty else pd.DataFrame()
            hyper_fe.append(h.iloc[0]['fold_enrichment'] if not h.empty else 0)
            hypo_fe.append(-l.iloc[0]['fold_enrichment'] if not l.empty else 0)

        fig, ax = plt.subplots(figsize=(10, max(5, len(pathways) * 0.35)))

        y = np.arange(len(pathways))
        ax.barh(y, hyper_fe, height=0.4, label='Hyper DMG', color='#d73027', alpha=0.8)
        ax.barh(y, hypo_fe, height=0.4, label='Hypo DMG', color='#4575b4', alpha=0.8)

        ax.set_yticks(y)
        ax.set_yticklabels(pathways, fontsize=8)
        ax.set_xlabel('Fold Enrichment (← Hypo | Hyper →)')
        ax.axvline(0, color='black', linewidth=0.5)
        ax.set_title(f'KEGG Enrichment: Hyper vs Hypo DMGs ({comp})')
        ax.legend(loc='lower right')

        fig.tight_layout()
        save_figure(fig, f'KEGG_hyper_vs_hypo_{comp}')


# ============================================================
# 9. Summary Statistics
# ============================================================
def generate_summary(dmg_sets, deg_sets, kegg_results, go_results, cog_results, motif_sets):
    """Generate comprehensive summary table."""
    print("\n" + "=" * 60)
    print("Generating summary statistics...")
    print("=" * 60)

    rows = []
    for comp in COMPARISONS:
        row = {'comparison': comp}

        # DMG counts
        row['DMG_6mA_hyper'] = len(dmg_sets.get(f"{comp}_6mA_hyper", []))
        row['DMG_6mA_hypo'] = len(dmg_sets.get(f"{comp}_6mA_hypo", []))
        row['DMG_4mC_hyper'] = len(dmg_sets.get(f"{comp}_4mC_hyper", []))
        row['DMG_4mC_hypo'] = len(dmg_sets.get(f"{comp}_4mC_hypo", []))
        row['DMG_combined'] = len(dmg_sets.get(f"{comp}_combined_all", []))

        # DEG counts
        row['DEG_up'] = len(deg_sets.get(f"{comp}_up", []))
        row['DEG_down'] = len(deg_sets.get(f"{comp}_down", []))

        # Significant KEGG pathways
        key_dmg = f"DMG_{comp}_combined_all"
        key_deg = f"DEG_{comp}_all"
        if key_dmg in kegg_results and not kegg_results[key_dmg].empty:
            row['KEGG_DMG_sig'] = len(kegg_results[key_dmg][kegg_results[key_dmg]['padj'] < 0.05])
        else:
            row['KEGG_DMG_sig'] = 0
        if key_deg in kegg_results and not kegg_results[key_deg].empty:
            row['KEGG_DEG_sig'] = len(kegg_results[key_deg][kegg_results[key_deg]['padj'] < 0.05])
        else:
            row['KEGG_DEG_sig'] = 0

        # Significant GO terms
        key_dmg_go = f"DMG_{comp}_combined_all"
        key_deg_go = f"DEG_{comp}_all"
        if key_dmg_go in go_results and not go_results[key_dmg_go].empty:
            row['GO_DMG_sig'] = len(go_results[key_dmg_go][go_results[key_dmg_go]['padj'] < 0.05])
        else:
            row['GO_DMG_sig'] = 0
        if key_deg_go in go_results and not go_results[key_deg_go].empty:
            row['GO_DEG_sig'] = len(go_results[key_deg_go][go_results[key_deg_go]['padj'] < 0.05])
        else:
            row['GO_DEG_sig'] = 0

        # Significant COG categories
        key_dmg_cog = f"DMG_{comp}_combined_all"
        key_deg_cog = f"DEG_{comp}_all"
        if key_dmg_cog in cog_results and not cog_results[key_dmg_cog].empty:
            row['COG_DMG_sig'] = len(cog_results[key_dmg_cog][cog_results[key_dmg_cog]['padj'] < 0.05])
        else:
            row['COG_DMG_sig'] = 0
        if key_deg_cog in cog_results and not cog_results[key_deg_cog].empty:
            row['COG_DEG_sig'] = len(cog_results[key_deg_cog][cog_results[key_deg_cog]['padj'] < 0.05])
        else:
            row['COG_DEG_sig'] = 0

        rows.append(row)

    summary_df = pd.DataFrame(rows)
    summary_df.to_csv(ANALYSIS_DIR / "tables" / "enrichment_summary.tsv", sep='\t', index=False)
    print("\nEnrichment Summary:")
    print(summary_df.to_string(index=False))

    return summary_df


# ============================================================
# Main
# ============================================================
def main():
    print("=" * 60)
    print("DMG Functional Enrichment Analysis")
    print(f"Output: {ANALYSIS_DIR}")
    print("=" * 60)

    # 1. Load data
    methyl, gene_master, gene_annot, cog = load_data()

    # 2. Extract DMGs
    dmg_sets, dmg_summary = extract_dmgs(methyl)

    # 3. Extract DEGs
    deg_sets = extract_degs(gene_master)

    # 4. KEGG enrichment
    kegg_results = run_kegg_enrichment(dmg_sets, deg_sets, gene_annot)

    # 5. GO enrichment
    go_results = run_go_enrichment(dmg_sets, deg_sets, gene_annot)

    # 6. COG enrichment
    cog_results = run_cog_enrichment(dmg_sets, deg_sets, cog)

    # 7. Motif-specific classification
    motif_sets, motif_cog_results, motif_kegg_results = run_motif_classification(cog, gene_annot)

    # 8. Figures
    print("\n" + "=" * 60)
    print("Generating figures...")
    print("=" * 60)
    plot_dmg_summary(dmg_summary)
    plot_kegg_enrichment_top(kegg_results)
    plot_cog_distribution_comparison(cog_results, dmg_sets, deg_sets, cog)
    plot_motif_cog_comparison(motif_cog_results, motif_sets, cog)
    plot_kegg_dmg_vs_deg_heatmap(kegg_results)
    plot_kegg_enrichment_by_direction(kegg_results)

    # 9. Summary
    summary = generate_summary(dmg_sets, deg_sets, kegg_results, go_results, cog_results, motif_sets)

    print("\n" + "=" * 60)
    print("DONE. All results saved to:")
    print(f"  {ANALYSIS_DIR}")
    print("=" * 60)


if __name__ == '__main__':
    main()
