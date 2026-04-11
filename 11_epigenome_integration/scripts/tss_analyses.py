#!/usr/bin/env python3
"""
TSS-based Epigenome-Transcriptome Integration Analyses
Streptomyces coelicolor A3(2) M145

Analyses:
  A1: Import experimental TSS from Jeong et al. 2016 (Nature Commun.)
  B1: Metagene profile of methylation around TSS
  C1: Distance-stratified correlation (methylation vs expression)
  C2: Promoter vs gene-body methylation effect comparison
  D1: Sigma factor binding site (-35/-10) overlap with methylation
"""

import pandas as pd
import numpy as np
from pathlib import Path
import re
from collections import defaultdict
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
import matplotlib.ticker as ticker
import warnings
warnings.filterwarnings('ignore')

# ── Paths ──
GFF_PATH = "/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/genomic.gff"
JEONG_XLSX = "/Users/okaban/bioinfo/rna-seq/reference/Jeong_Nat-Commun_2016/41467_2016_BFncomms11605_MOESM1058_ESM.xlsx"
METHYL_SITES = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/01_integration/high_confidence_sites_weighted.csv"
INTEGRATED_CSV = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/01_integration/integrated_methyl_expression_weighted.csv"
GENOME_FASTA = "/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/GCF_000203835.1_ASM20383v1_genomic.fna"
OUTPUT_DIR = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/18_tss_analyses")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Plot style ──
plt.rcParams.update({
    'font.size': 10,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'font.family': 'sans-serif',
})


# ============================================================
# A1: Import experimental TSS from Jeong et al. 2016
# ============================================================

def parse_gff(gff_path):
    """Parse GFF to extract gene coords with SCO -> SC_RS mapping."""
    genes = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'gene':
                continue
            attrs = parts[8]
            m_new = re.search(r'locus_tag=([^;\s]+)', attrs)
            m_old = re.search(r'old_locus_tag=([^;\s]+)', attrs)
            m_name = re.search(r'Name=([^;\s]+)', attrs)
            m_gene = re.search(r';gene=([^;\s]+)', attrs)
            if m_new:
                genes.append({
                    'gene_id': m_new.group(1),
                    'old_locus_tag': m_old.group(1) if m_old else None,
                    'gene_name': m_gene.group(1) if m_gene else None,
                    'chrom': parts[0],
                    'start': int(parts[3]),
                    'end': int(parts[4]),
                    'strand': parts[6],
                    'gff_tss': int(parts[3]) if parts[6] == '+' else int(parts[4]),
                })
    return pd.DataFrame(genes)


def load_jeong_tss(xlsx_path):
    """Load experimental TSS from Jeong et al. 2016 Supplementary Data 1."""
    import openpyxl
    wb = openpyxl.load_workbook(xlsx_path, read_only=True)
    ws = wb['Sheet1']
    rows = list(ws.iter_rows(min_row=3, values_only=True))  # skip description + header
    records = []
    for r in rows:
        if r[0] is None:
            continue
        records.append({
            'tss_id': r[0],
            'tss_position': int(r[1]),
            'strand': r[2],
            'abundance': int(r[3]) if r[3] else 0,
            'category': r[4],
            'sco_gene': r[5],
            'reference': r[6] if r[6] else '',
        })
    wb.close()
    return pd.DataFrame(records)


def run_a1():
    """A1: Build comprehensive TSS table merging experimental + GFF-based."""
    print("=" * 60)
    print("A1: Importing experimental TSS (Jeong et al. 2016)")
    print("=" * 60)

    genes_df = parse_gff(GFF_PATH)
    print(f"  GFF genes loaded: {len(genes_df)}")

    jeong_df = load_jeong_tss(JEONG_XLSX)
    print(f"  Jeong TSSs loaded: {len(jeong_df)}")
    print(f"  Categories: {jeong_df['category'].value_counts().to_dict()}")

    # Filter primary TSSs only (category P) for main analysis
    primary_tss = jeong_df[jeong_df['category'] == 'P'].copy()
    print(f"  Primary TSSs (P): {len(primary_tss)}")

    # Build SCO -> SC_RS mapping
    sco_to_scrs = {}
    for _, g in genes_df.iterrows():
        if g['old_locus_tag']:
            sco_to_scrs[g['old_locus_tag']] = g['gene_id']

    # Map Jeong gene IDs to SC_RS
    primary_tss['gene_id'] = primary_tss['sco_gene'].map(sco_to_scrs)
    mapped = primary_tss.dropna(subset=['gene_id'])
    print(f"  Primary TSSs mapped to SC_RS: {len(mapped)} / {len(primary_tss)}")

    # For genes with multiple primary TSSs, keep the most abundant
    mapped_best = mapped.sort_values('abundance', ascending=False).drop_duplicates('gene_id', keep='first')
    print(f"  Unique genes with experimental TSS: {len(mapped_best)}")

    # Merge into genes_df
    tss_map = mapped_best.set_index('gene_id')['tss_position'].to_dict()
    genes_df['experimental_tss'] = genes_df['gene_id'].map(tss_map)
    genes_df['tss_source'] = np.where(genes_df['experimental_tss'].notna(), 'Jeong2016_dRNA-seq', 'GFF_annotation')
    genes_df['tss'] = genes_df['experimental_tss'].fillna(genes_df['gff_tss']).astype(int)

    # Evaluate GFF vs experimental TSS offset
    has_both = genes_df.dropna(subset=['experimental_tss']).copy()
    has_both['tss_offset'] = has_both['experimental_tss'].astype(int) - has_both['gff_tss'].astype(int)
    # For minus strand, offset sense is reversed
    has_both.loc[has_both['strand'] == '-', 'tss_offset'] *= -1

    print(f"\n  TSS offset (experimental - GFF, strand-corrected):")
    print(f"    Median: {has_both['tss_offset'].median():.0f} bp")
    print(f"    Mean:   {has_both['tss_offset'].mean():.1f} bp")
    print(f"    Std:    {has_both['tss_offset'].std():.1f} bp")
    print(f"    Within ±50 bp:  {(has_both['tss_offset'].abs() <= 50).sum()} ({(has_both['tss_offset'].abs() <= 50).mean()*100:.1f}%)")
    print(f"    Within ±100 bp: {(has_both['tss_offset'].abs() <= 100).sum()} ({(has_both['tss_offset'].abs() <= 100).mean()*100:.1f}%)")

    # Save comprehensive TSS table
    genes_df.to_csv(OUTPUT_DIR / 'comprehensive_tss_table.csv', index=False)

    # Save offset histogram
    fig, ax = plt.subplots(figsize=(8, 4))
    offsets = has_both['tss_offset'].clip(-500, 500)
    ax.hist(offsets, bins=100, color='steelblue', edgecolor='none', alpha=0.8)
    ax.axvline(0, color='red', ls='--', lw=1, label='Perfect match')
    ax.axvline(offsets.median(), color='orange', ls='--', lw=1, label=f'Median={offsets.median():.0f} bp')
    ax.set_xlabel('TSS offset (experimental − GFF annotation, bp)')
    ax.set_ylabel('Number of genes')
    ax.set_title(f'GFF vs Experimental TSS Offset (n={len(has_both)})')
    ax.legend()
    fig.savefig(OUTPUT_DIR / 'A1_tss_offset_histogram.png')
    plt.close(fig)

    # Save all Jeong TSSs (including secondary, internal, antisense)
    jeong_df['gene_id'] = jeong_df['sco_gene'].map(sco_to_scrs)
    jeong_df.to_csv(OUTPUT_DIR / 'jeong2016_all_tss.csv', index=False)

    # Summary stats
    summary = {
        'total_jeong_tss': len(jeong_df),
        'primary_tss': len(primary_tss),
        'mapped_to_scrs': len(mapped),
        'unique_genes_with_exp_tss': len(mapped_best),
        'total_gff_genes': len(genes_df),
        'genes_with_exp_tss': (genes_df['tss_source'] == 'Jeong2016_dRNA-seq').sum(),
        'genes_gff_only': (genes_df['tss_source'] == 'GFF_annotation').sum(),
        'median_offset_bp': has_both['tss_offset'].median(),
        'mean_offset_bp': has_both['tss_offset'].mean(),
    }

    print(f"\n  Summary:")
    print(f"    Genes with experimental TSS: {summary['genes_with_exp_tss']} ({summary['genes_with_exp_tss']/summary['total_gff_genes']*100:.1f}%)")
    print(f"    Genes with GFF-only TSS:     {summary['genes_gff_only']} ({summary['genes_gff_only']/summary['total_gff_genes']*100:.1f}%)")

    return genes_df, jeong_df, summary


# ============================================================
# B1: Metagene methylation profile around TSS
# ============================================================

def run_b1(genes_df):
    """B1: Metagene profile of 6mA/4mC density around TSS."""
    print("\n" + "=" * 60)
    print("B1: Metagene methylation profile around TSS")
    print("=" * 60)

    methyl_df = pd.read_csv(METHYL_SITES)
    print(f"  Methylation sites loaded: {len(methyl_df)}")

    window = 1000  # ±1000 bp from TSS
    bin_size = 20  # 20 bp bins

    bins = np.arange(-window, window + bin_size, bin_size)
    bin_centers = (bins[:-1] + bins[1:]) / 2

    # Pre-index methylation sites by position for fast lookup
    methyl_positions = {}
    for tp in ['T1', 'T2', 'T3']:
        for mod in ['6mA', '4mC']:
            subset = methyl_df[(methyl_df['timepoint'] == tp) & (methyl_df['mod_type'] == mod)]
            methyl_positions[(tp, mod)] = set(subset['position'].values)

    results = {(tp, mod): np.zeros(len(bin_centers)) for tp in ['T1', 'T2', 'T3'] for mod in ['6mA', '4mC']}
    gene_count = 0

    for _, gene in genes_df.iterrows():
        tss = gene['tss']
        strand = gene['strand']
        gene_count += 1

        for tp in ['T1', 'T2', 'T3']:
            for mod in ['6mA', '4mC']:
                positions = methyl_positions[(tp, mod)]
                for pos in range(tss - window, tss + window + 1):
                    if pos in positions:
                        if strand == '+':
                            rel_pos = pos - tss
                        else:
                            rel_pos = tss - pos
                        bin_idx = np.searchsorted(bins, rel_pos, side='right') - 1
                        if 0 <= bin_idx < len(bin_centers):
                            results[(tp, mod)][bin_idx] += 1

    # Normalize per gene per kb
    for key in results:
        results[key] = results[key] / gene_count / (bin_size / 1000)

    # ── Figure 1: Combined metagene plot (6mA + 4mC, all timepoints) ──
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    colors_tp = {'T1': '#2166ac', 'T2': '#d6604d', 'T3': '#4daf4a'}
    tp_labels = {'T1': 'T1 (early)', 'T2': 'T2 (mid)', 'T3': 'T3 (late)'}

    for mod, ax in zip(['6mA', '4mC'], axes):
        for tp in ['T1', 'T2', 'T3']:
            # Smooth with rolling average
            y = pd.Series(results[(tp, mod)]).rolling(3, center=True, min_periods=1).mean().values
            ax.plot(bin_centers, y, color=colors_tp[tp], lw=1.5, label=tp_labels[tp])
        ax.axvline(0, color='grey', ls='--', lw=0.8, alpha=0.6)
        ax.axvspan(-50, -10, alpha=0.08, color='purple', label='-10 box region')
        ax.axvspan(-60, -35, alpha=0.05, color='blue')
        ax.set_ylabel(f'{mod} density\n(sites / gene / kb)')
        ax.set_title(f'{mod} methylation around TSS')
        ax.legend(loc='upper right', fontsize=8)

    axes[1].set_xlabel('Distance from TSS (bp)')
    axes[0].text(-window * 0.95, axes[0].get_ylim()[1] * 0.9, 'Upstream', fontsize=9, ha='left', color='grey')
    axes[0].text(window * 0.95, axes[0].get_ylim()[1] * 0.9, 'Downstream', fontsize=9, ha='right', color='grey')

    fig.suptitle(f'Metagene methylation profile (n = {gene_count} genes)', fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / 'B1_metagene_methylation_profile.png')
    fig.savefig(OUTPUT_DIR / 'B1_metagene_methylation_profile.pdf')
    plt.close(fig)

    # ── Figure 2: Zoomed promoter region (−500 to +200) ──
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    zoom_mask = (bin_centers >= -500) & (bin_centers <= 200)

    for mod, ax in zip(['6mA', '4mC'], axes):
        for tp in ['T1', 'T2', 'T3']:
            y = pd.Series(results[(tp, mod)]).rolling(3, center=True, min_periods=1).mean().values
            ax.plot(bin_centers[zoom_mask], y[zoom_mask], color=colors_tp[tp], lw=1.8, label=tp_labels[tp])
        ax.axvline(0, color='green', ls='-', lw=1.2, alpha=0.7, label='TSS')
        ax.axvspan(-45, -30, alpha=0.15, color='purple', label='-35 box')
        ax.axvspan(-15, -5, alpha=0.15, color='orange', label='-10 box')
        ax.set_ylabel(f'{mod} density\n(sites / gene / kb)')
        ax.legend(loc='upper right', fontsize=8)

    axes[1].set_xlabel('Distance from TSS (bp)')
    fig.suptitle('Methylation in promoter region (zoomed)', fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / 'B1_metagene_promoter_zoom.png')
    fig.savefig(OUTPUT_DIR / 'B1_metagene_promoter_zoom.pdf')
    plt.close(fig)

    # ── Figure 3: Separate by TSS source (experimental vs GFF) ──
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

    for source, ls, label_sfx in [('Jeong2016_dRNA-seq', '-', 'exp.'), ('GFF_annotation', '--', 'GFF')]:
        sub_genes = genes_df[genes_df['tss_source'] == source]
        n_genes_src = len(sub_genes)
        src_results = {mod: np.zeros(len(bin_centers)) for mod in ['6mA', '4mC']}

        for _, gene in sub_genes.iterrows():
            tss = gene['tss']
            strand = gene['strand']
            for mod in ['6mA', '4mC']:
                all_pos = set()
                for tp in ['T1', 'T2', 'T3']:
                    all_pos.update(methyl_positions[(tp, mod)])
                for pos in range(tss - window, tss + window + 1):
                    if pos in all_pos:
                        rel_pos = (pos - tss) if strand == '+' else (tss - pos)
                        bin_idx = np.searchsorted(bins, rel_pos, side='right') - 1
                        if 0 <= bin_idx < len(bin_centers):
                            src_results[mod][bin_idx] += 1

        for mod_idx, mod in enumerate(['6mA', '4mC']):
            y = pd.Series(src_results[mod] / n_genes_src / (bin_size / 1000)).rolling(3, center=True, min_periods=1).mean().values
            axes[mod_idx].plot(bin_centers, y, ls=ls, lw=1.5,
                               label=f'{label_sfx} (n={n_genes_src})',
                               color='steelblue' if source == 'Jeong2016_dRNA-seq' else 'coral')

    for ax, mod in zip(axes, ['6mA', '4mC']):
        ax.axvline(0, color='grey', ls='--', lw=0.8)
        ax.set_ylabel(f'{mod} density\n(sites / gene / kb)')
        ax.set_title(f'{mod}: Experimental TSS vs GFF-based TSS')
        ax.legend(fontsize=8)
    axes[1].set_xlabel('Distance from TSS (bp)')
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / 'B1_metagene_by_tss_source.png')
    plt.close(fig)

    # Save numerical data
    metagene_data = pd.DataFrame({'bin_center': bin_centers})
    for tp in ['T1', 'T2', 'T3']:
        for mod in ['6mA', '4mC']:
            metagene_data[f'{mod}_{tp}'] = results[(tp, mod)]
    metagene_data.to_csv(OUTPUT_DIR / 'B1_metagene_data.csv', index=False)

    print(f"  Figures saved to {OUTPUT_DIR}/B1_*")
    return results, bin_centers


# ============================================================
# C1: Distance-stratified correlation analysis
# ============================================================

def run_c1(genes_df):
    """C1: Stratify methylation by TSS distance and correlate with expression."""
    print("\n" + "=" * 60)
    print("C1: Distance-stratified methylation-expression correlation")
    print("=" * 60)

    methyl_df = pd.read_csv(METHYL_SITES)
    integrated_df = pd.read_csv(INTEGRATED_CSV)

    # Distance bins: upstream promoter → around TSS → gene body
    distance_bins = [
        (-300, -200, 'distal_upstream'),
        (-200, -100, 'mid_upstream'),
        (-100, -50,  'proximal_upstream'),
        (-50,    0,  'core_promoter'),
        (0,    +50,  'TSS_proximal'),
        (+50,  +200, 'early_gene_body'),
        (+200, +500, 'mid_gene_body'),
        (+500, +1000, 'distal_gene_body'),
    ]

    # For each gene, assign each methylation site to a distance bin
    gene_tss = genes_df.set_index('gene_id')[['tss', 'strand', 'start', 'end']].to_dict('index')

    records = []
    for _, site in methyl_df.iterrows():
        pos = site['position']
        tp = site['timepoint']
        mod = site['mod_type']
        freq = site['weighted_mod_freq']

        # Find nearest gene and compute distance to TSS
        for gid, ginfo in gene_tss.items():
            tss = ginfo['tss']
            strand = ginfo['strand']
            g_start = ginfo['start']
            g_end = ginfo['end']

            # Check if site is within gene region ± 500 bp
            if strand == '+':
                rel_pos = pos - tss
                if -300 <= rel_pos <= 1000:
                    records.append({
                        'gene_id': gid, 'position': pos, 'rel_pos': rel_pos,
                        'mod_type': mod, 'timepoint': tp, 'freq': freq
                    })
            else:
                rel_pos = tss - pos
                if -300 <= rel_pos <= 1000:
                    records.append({
                        'gene_id': gid, 'position': pos, 'rel_pos': rel_pos,
                        'mod_type': mod, 'timepoint': tp, 'freq': freq
                    })

    site_gene_df = pd.DataFrame(records)
    print(f"  Site-gene assignments: {len(site_gene_df)}")

    # Assign distance bin
    def assign_bin(rel_pos):
        for lo, hi, name in distance_bins:
            if lo <= rel_pos < hi:
                return name
        return None

    site_gene_df['distance_bin'] = site_gene_df['rel_pos'].apply(assign_bin)
    site_gene_df = site_gene_df.dropna(subset=['distance_bin'])

    # Aggregate: per gene, per bin, compute mean freq change across timepoints
    # We need T1 vs T2 freq changes per bin
    pivot = site_gene_df.pivot_table(
        index=['gene_id', 'distance_bin', 'mod_type'],
        columns='timepoint',
        values='freq',
        aggfunc='mean'
    ).reset_index()

    comparisons = [('T2', 'T1', 'T2_vs_T1'), ('T3', 'T1', 'T3_vs_T1'), ('T3', 'T2', 'T3_vs_T2')]
    for tp_a, tp_b, label in comparisons:
        if tp_a in pivot.columns and tp_b in pivot.columns:
            pivot[f'methyl_change_{label}'] = pivot[tp_a] - pivot[tp_b]

    # Merge with expression data
    expr_cols = ['gene_id', 'log2FC_T2_vs_T1', 'padj_T2_vs_T1',
                 'log2FC_T3_vs_T1', 'padj_T3_vs_T1',
                 'log2FC_T3_vs_T2', 'padj_T3_vs_T2']
    available_cols = [c for c in expr_cols if c in integrated_df.columns]
    merged = pivot.merge(integrated_df[available_cols], on='gene_id', how='inner')

    # Correlation per distance bin
    corr_results = []
    for mod in ['6mA', '4mC']:
        for label in ['T2_vs_T1', 'T3_vs_T1', 'T3_vs_T2']:
            mc = f'methyl_change_{label}'
            ec = f'log2FC_{label}'
            if mc not in merged.columns or ec not in merged.columns:
                continue
            for _, (lo, hi, bname) in enumerate(distance_bins):
                sub = merged[(merged['distance_bin'] == bname) & (merged['mod_type'] == mod)].dropna(subset=[mc, ec])
                if len(sub) >= 10:
                    r, p = stats.spearmanr(sub[mc], sub[ec])
                    corr_results.append({
                        'mod_type': mod, 'comparison': label,
                        'distance_bin': bname, 'bin_range': f'{lo} to {hi}',
                        'n_genes': len(sub), 'spearman_r': r, 'p_value': p,
                    })

    corr_df = pd.DataFrame(corr_results)
    corr_df.to_csv(OUTPUT_DIR / 'C1_distance_stratified_correlation.csv', index=False)
    print(f"\n  Correlation results:")
    for _, row in corr_df[corr_df['comparison'] == 'T2_vs_T1'].iterrows():
        sig = '***' if row['p_value'] < 0.001 else '**' if row['p_value'] < 0.01 else '*' if row['p_value'] < 0.05 else ''
        print(f"    {row['mod_type']} {row['distance_bin']:>20s}: r={row['spearman_r']:+.3f} p={row['p_value']:.1e} n={row['n_genes']} {sig}")

    # ── Figure: Heatmap of correlations by distance bin ──
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    bin_order = [b[2] for b in distance_bins]

    for idx, mod in enumerate(['6mA', '4mC']):
        ax = axes[idx]
        sub = corr_df[corr_df['mod_type'] == mod]
        if len(sub) == 0:
            continue
        pivot_r = sub.pivot_table(index='distance_bin', columns='comparison', values='spearman_r')
        pivot_p = sub.pivot_table(index='distance_bin', columns='comparison', values='p_value')
        # Reorder
        comp_order = ['T2_vs_T1', 'T3_vs_T1', 'T3_vs_T2']
        existing_bins = [b for b in bin_order if b in pivot_r.index]
        existing_comps = [c for c in comp_order if c in pivot_r.columns]
        pivot_r = pivot_r.reindex(index=existing_bins, columns=existing_comps)
        pivot_p = pivot_p.reindex(index=existing_bins, columns=existing_comps)

        im = ax.imshow(pivot_r.values, cmap='RdBu_r', vmin=-0.3, vmax=0.3, aspect='auto')
        ax.set_xticks(range(len(existing_comps)))
        ax.set_xticklabels(existing_comps, rotation=45, ha='right')
        ax.set_yticks(range(len(existing_bins)))
        ax.set_yticklabels(existing_bins)

        # Annotate with r values and significance
        for i in range(len(existing_bins)):
            for j in range(len(existing_comps)):
                r_val = pivot_r.values[i, j]
                p_val = pivot_p.values[i, j]
                if np.isnan(r_val):
                    continue
                sig = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else ''
                ax.text(j, i, f'{r_val:.2f}{sig}', ha='center', va='center', fontsize=8,
                        color='white' if abs(r_val) > 0.15 else 'black')

        ax.set_title(f'{mod}')
        plt.colorbar(im, ax=ax, label='Spearman r', shrink=0.8)

    fig.suptitle('Methylation-Expression Correlation by Distance from TSS', fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / 'C1_distance_correlation_heatmap.png')
    fig.savefig(OUTPUT_DIR / 'C1_distance_correlation_heatmap.pdf')
    plt.close(fig)

    print(f"  Figures saved to {OUTPUT_DIR}/C1_*")
    return corr_df


# ============================================================
# C2: Promoter vs gene-body methylation effect
# ============================================================

def run_c2(genes_df):
    """C2: Compare promoter vs gene-body methylation impact on expression."""
    print("\n" + "=" * 60)
    print("C2: Promoter vs gene-body methylation effect comparison")
    print("=" * 60)

    methyl_df = pd.read_csv(METHYL_SITES)
    integrated_df = pd.read_csv(INTEGRATED_CSV)

    gene_tss = genes_df.set_index('gene_id')[['tss', 'strand', 'start', 'end']].to_dict('index')

    # Classify each site as promoter (-300 to +50 from TSS) or gene_body (+51 to gene end)
    site_records = []
    for _, site in methyl_df.iterrows():
        pos = site['position']
        tp = site['timepoint']
        mod = site['mod_type']
        freq = site['weighted_mod_freq']

        for gid, ginfo in gene_tss.items():
            tss = ginfo['tss']
            strand = ginfo['strand']
            g_start, g_end = ginfo['start'], ginfo['end']

            if strand == '+':
                rel_pos = pos - tss
                in_gene = g_start <= pos <= g_end
            else:
                rel_pos = tss - pos
                in_gene = g_start <= pos <= g_end

            if -300 <= rel_pos <= 50:
                region = 'promoter'
            elif rel_pos > 50 and in_gene:
                region = 'gene_body'
            else:
                continue

            site_records.append({
                'gene_id': gid, 'region': region, 'mod_type': mod,
                'timepoint': tp, 'freq': freq
            })

    site_df = pd.DataFrame(site_records)
    print(f"  Classified sites: {len(site_df)}")
    print(f"  Promoter: {(site_df['region']=='promoter').sum()}, Gene body: {(site_df['region']=='gene_body').sum()}")

    # Aggregate per gene per region per timepoint
    agg = site_df.groupby(['gene_id', 'region', 'mod_type', 'timepoint']).agg(
        mean_freq=('freq', 'mean'),
        n_sites=('freq', 'count')
    ).reset_index()

    # Pivot timepoints
    pivot = agg.pivot_table(
        index=['gene_id', 'region', 'mod_type'],
        columns='timepoint', values='mean_freq', aggfunc='mean'
    ).reset_index()

    for tp_a, tp_b, label in [('T2', 'T1', 'T2_vs_T1'), ('T3', 'T1', 'T3_vs_T1'), ('T3', 'T2', 'T3_vs_T2')]:
        if tp_a in pivot.columns and tp_b in pivot.columns:
            pivot[f'methyl_change_{label}'] = pivot[tp_a] - pivot[tp_b]

    # Merge with expression
    expr_cols = ['gene_id', 'log2FC_T2_vs_T1', 'padj_T2_vs_T1',
                 'log2FC_T3_vs_T1', 'padj_T3_vs_T1',
                 'log2FC_T3_vs_T2', 'padj_T3_vs_T2']
    available_cols = [c for c in expr_cols if c in integrated_df.columns]
    merged = pivot.merge(integrated_df[available_cols], on='gene_id', how='inner')

    # Correlation per region
    corr_results = []
    for mod in ['6mA', '4mC']:
        for label in ['T2_vs_T1', 'T3_vs_T1', 'T3_vs_T2']:
            mc = f'methyl_change_{label}'
            ec = f'log2FC_{label}'
            if mc not in merged.columns or ec not in merged.columns:
                continue
            for region in ['promoter', 'gene_body']:
                sub = merged[(merged['region'] == region) & (merged['mod_type'] == mod)].dropna(subset=[mc, ec])
                if len(sub) >= 10:
                    r, p = stats.spearmanr(sub[mc], sub[ec])
                    corr_results.append({
                        'mod_type': mod, 'comparison': label,
                        'region': region, 'n_genes': len(sub),
                        'spearman_r': r, 'p_value': p,
                    })

    corr_df = pd.DataFrame(corr_results)
    corr_df.to_csv(OUTPUT_DIR / 'C2_promoter_vs_genebody_correlation.csv', index=False)

    print(f"\n  Correlation results (Promoter vs Gene body):")
    for _, row in corr_df.iterrows():
        sig = '***' if row['p_value'] < 0.001 else '**' if row['p_value'] < 0.01 else '*' if row['p_value'] < 0.05 else ''
        print(f"    {row['mod_type']} {row['region']:>10s} {row['comparison']}: r={row['spearman_r']:+.3f} p={row['p_value']:.1e} n={row['n_genes']} {sig}")

    # ── Figure: Side-by-side scatter plots ──
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    comp_labels = ['T2_vs_T1', 'T3_vs_T1', 'T3_vs_T2']

    for row_idx, mod in enumerate(['6mA', '4mC']):
        for col_idx, label in enumerate(comp_labels):
            ax = axes[row_idx, col_idx]
            mc = f'methyl_change_{label}'
            ec = f'log2FC_{label}'
            if mc not in merged.columns or ec not in merged.columns:
                continue

            for region, color, marker in [('promoter', '#d6604d', 'o'), ('gene_body', '#2166ac', 's')]:
                sub = merged[(merged['region'] == region) & (merged['mod_type'] == mod)].dropna(subset=[mc, ec])
                if len(sub) > 0:
                    ax.scatter(sub[mc], sub[ec], c=color, marker=marker, s=15, alpha=0.4, label=region, edgecolors='none')

            ax.axhline(0, color='grey', ls='--', lw=0.5)
            ax.axvline(0, color='grey', ls='--', lw=0.5)
            ax.set_xlabel(f'Δ methylation freq ({label})')
            ax.set_ylabel(f'log2FC expression ({label})')
            ax.set_title(f'{mod} — {label}')
            if row_idx == 0 and col_idx == 0:
                ax.legend(fontsize=8, markerscale=2)

    fig.suptitle('Promoter vs Gene-body Methylation Effect on Expression', fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / 'C2_promoter_vs_genebody_scatter.png')
    fig.savefig(OUTPUT_DIR / 'C2_promoter_vs_genebody_scatter.pdf')
    plt.close(fig)

    # ── Figure: Bar chart comparing r values ──
    fig, ax = plt.subplots(figsize=(10, 5))
    bar_data = corr_df[corr_df['comparison'] == 'T2_vs_T1'].copy()
    if len(bar_data) > 0:
        x = np.arange(len(bar_data))
        labels = [f"{r['mod_type']}\n{r['region']}" for _, r in bar_data.iterrows()]
        colors = ['#d6604d' if r['region'] == 'promoter' else '#2166ac' for _, r in bar_data.iterrows()]
        bars = ax.bar(x, bar_data['spearman_r'], color=colors, edgecolor='white', width=0.6)
        for i, (_, r) in enumerate(bar_data.iterrows()):
            sig = '***' if r['p_value'] < 0.001 else '**' if r['p_value'] < 0.01 else '*' if r['p_value'] < 0.05 else 'ns'
            ax.text(i, r['spearman_r'] + 0.01 * np.sign(r['spearman_r']),
                    f"r={r['spearman_r']:.3f}\n{sig}\nn={r['n_genes']}",
                    ha='center', va='bottom' if r['spearman_r'] >= 0 else 'top', fontsize=8)
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.axhline(0, color='grey', lw=0.5)
        ax.set_ylabel('Spearman r')
        ax.set_title('Methylation-Expression Correlation: Promoter vs Gene Body (T2 vs T1)')
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / 'C2_promoter_vs_genebody_barplot.png')
    fig.savefig(OUTPUT_DIR / 'C2_promoter_vs_genebody_barplot.pdf')
    plt.close(fig)

    print(f"  Figures saved to {OUTPUT_DIR}/C2_*")
    return corr_df


# ============================================================
# D1: Sigma factor binding site overlap with methylation
# ============================================================

def load_genome_sequence(fasta_path):
    """Load genome sequence from FASTA."""
    seq_parts = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                continue
            seq_parts.append(line.strip())
    return ''.join(seq_parts)


def find_sigma_motifs(sequence, genes_df):
    """Find -35 and -10 box motifs in promoter regions."""
    # σ70 (HrdB) consensus for Streptomyces:
    # -35 box: TTGAC[A/G] (relaxed: [TC]TGAC[ACGT])
    # -10 box: TA[ACGT][ACGT][ACGT]T (relaxed Pribnow: TA.{3}T)
    # For Streptomyces, σHrdB is the primary sigma factor

    sigma70_35 = re.compile(r'[TC]TGAC[ACGT]')
    sigma70_10 = re.compile(r'TA[ACGT]{3}T')

    # Also look for σWhiG (sporulation sigma factor) which is GC-rich
    # σBldN motif: CGC pattern

    motif_records = []
    genome_len = len(sequence)

    for _, gene in genes_df.iterrows():
        tss = gene['tss']
        strand = gene['strand']
        gid = gene['gene_id']

        # Extract promoter region (-100 to +10 from TSS)
        if strand == '+':
            prom_start = max(0, tss - 101)  # 0-indexed
            prom_end = min(genome_len, tss + 10)
            prom_seq = sequence[prom_start:prom_end]
        else:
            prom_start = max(0, tss - 11)
            prom_end = min(genome_len, tss + 100)
            prom_seq = sequence[prom_start:prom_end]
            # Reverse complement
            comp = str.maketrans('ACGTacgt', 'TGCAtgca')
            prom_seq = prom_seq.translate(comp)[::-1]

        prom_seq_upper = prom_seq.upper()

        # Search for -35 box (expected around position 55-70 in our 111-bp window = -45 to -30 from TSS)
        for m in sigma70_35.finditer(prom_seq_upper):
            motif_start = m.start()
            # Convert to TSS-relative position
            rel_pos = motif_start - 100  # since prom_start is tss-101
            if -60 <= rel_pos <= -25:
                motif_records.append({
                    'gene_id': gid, 'motif_type': '-35_box', 'motif_seq': m.group(),
                    'rel_pos': rel_pos, 'strand': strand,
                })

        # Search for -10 box (expected around position 85-100 = -15 to 0 from TSS)
        for m in sigma70_10.finditer(prom_seq_upper):
            motif_start = m.start()
            rel_pos = motif_start - 100
            if -25 <= rel_pos <= 0:
                motif_records.append({
                    'gene_id': gid, 'motif_type': '-10_box', 'motif_seq': m.group(),
                    'rel_pos': rel_pos, 'strand': strand,
                })

    return pd.DataFrame(motif_records)


def run_d1(genes_df):
    """D1: σ factor binding site overlap with methylation sites."""
    print("\n" + "=" * 60)
    print("D1: Sigma factor binding site – methylation overlap")
    print("=" * 60)

    methyl_df = pd.read_csv(METHYL_SITES)

    # Load genome
    print("  Loading genome sequence...")
    genome_seq = load_genome_sequence(GENOME_FASTA)
    print(f"  Genome length: {len(genome_seq):,} bp")

    # Find sigma motifs
    print("  Searching for σ70 (HrdB) -35/-10 box motifs...")
    motif_df = find_sigma_motifs(genome_seq, genes_df)
    print(f"  Found {len(motif_df)} motif instances")
    print(f"    -35 box: {(motif_df['motif_type']=='-35_box').sum()}")
    print(f"    -10 box: {(motif_df['motif_type']=='-10_box').sum()}")

    # Genes with motifs
    genes_with_35 = set(motif_df[motif_df['motif_type'] == '-35_box']['gene_id'])
    genes_with_10 = set(motif_df[motif_df['motif_type'] == '-10_box']['gene_id'])
    print(f"  Genes with -35 box: {len(genes_with_35)}")
    print(f"  Genes with -10 box: {len(genes_with_10)}")
    print(f"  Genes with both:    {len(genes_with_35 & genes_with_10)}")

    # Map methylation sites to TSS-relative positions
    gene_tss = genes_df.set_index('gene_id')[['tss', 'strand']].to_dict('index')

    methyl_rel = []
    for _, site in methyl_df.iterrows():
        pos = site['position']
        for gid, ginfo in gene_tss.items():
            tss = ginfo['tss']
            strand = ginfo['strand']
            if strand == '+':
                rel_pos = pos - tss
            else:
                rel_pos = tss - pos
            if -100 <= rel_pos <= 10:
                methyl_rel.append({
                    'gene_id': gid, 'position': pos, 'rel_pos': rel_pos,
                    'mod_type': site['mod_type'], 'timepoint': site['timepoint'],
                    'freq': site['weighted_mod_freq'],
                })

    methyl_rel_df = pd.DataFrame(methyl_rel)
    print(f"\n  Methylation sites in promoter core (-100 to +10): {len(methyl_rel_df)}")

    # Check overlap: methylation within ±3 bp of motif boundaries
    overlap_records = []
    for _, motif in motif_df.iterrows():
        gid = motif['gene_id']
        motif_rel_start = motif['rel_pos']
        motif_len = len(motif['motif_seq'])
        motif_rel_end = motif_rel_start + motif_len

        nearby_methyl = methyl_rel_df[
            (methyl_rel_df['gene_id'] == gid) &
            (methyl_rel_df['rel_pos'] >= motif_rel_start - 3) &
            (methyl_rel_df['rel_pos'] <= motif_rel_end + 3)
        ]

        for _, m in nearby_methyl.iterrows():
            overlap_records.append({
                'gene_id': gid,
                'motif_type': motif['motif_type'],
                'motif_seq': motif['motif_seq'],
                'motif_rel_pos': motif_rel_start,
                'methyl_rel_pos': m['rel_pos'],
                'mod_type': m['mod_type'],
                'timepoint': m['timepoint'],
                'methyl_freq': m['freq'],
                'overlap_type': 'within' if motif_rel_start <= m['rel_pos'] <= motif_rel_end else 'flanking',
            })

    overlap_df = pd.DataFrame(overlap_records)
    overlap_df.to_csv(OUTPUT_DIR / 'D1_sigma_methylation_overlap.csv', index=False)

    print(f"  Methylation-motif overlaps: {len(overlap_df)}")
    if len(overlap_df) > 0:
        print(f"    Within motif: {(overlap_df['overlap_type']=='within').sum()}")
        print(f"    Flanking:     {(overlap_df['overlap_type']=='flanking').sum()}")
        print(f"    By motif type:")
        for mt in ['-35_box', '-10_box']:
            n = (overlap_df['motif_type'] == mt).sum()
            print(f"      {mt}: {n}")

    # ── Statistical test: is methylation enriched at σ binding sites? ──
    # Compare observed methylation density at motifs vs random promoter regions
    # Count methylation sites in -35 box region (-45 to -25) vs control (-100 to -60)

    methyl_at_35 = methyl_rel_df[(methyl_rel_df['rel_pos'] >= -45) & (methyl_rel_df['rel_pos'] <= -25)]
    methyl_at_10 = methyl_rel_df[(methyl_rel_df['rel_pos'] >= -20) & (methyl_rel_df['rel_pos'] <= 0)]
    methyl_ctrl  = methyl_rel_df[(methyl_rel_df['rel_pos'] >= -100) & (methyl_rel_df['rel_pos'] <= -60)]

    n_genes = len(genes_df)
    region_sizes = {'-35_region': 21, '-10_region': 21, 'control': 41}

    density_35 = len(methyl_at_35) / (n_genes * region_sizes['-35_region'] / 1000)
    density_10 = len(methyl_at_10) / (n_genes * region_sizes['-10_region'] / 1000)
    density_ctrl = len(methyl_ctrl) / (n_genes * region_sizes['control'] / 1000)

    print(f"\n  Methylation density (sites per gene per kb):")
    print(f"    -35 region (-45 to -25): {density_35:.2f}  (n={len(methyl_at_35)})")
    print(f"    -10 region (-20 to  0):  {density_10:.2f}  (n={len(methyl_at_10)})")
    print(f"    Control    (-100 to -60): {density_ctrl:.2f}  (n={len(methyl_ctrl)})")

    # Chi-squared test: -35 region vs control
    from scipy.stats import chi2_contingency
    obs_35 = len(methyl_at_35)
    obs_ctrl = len(methyl_ctrl)
    exp_ratio = region_sizes['-35_region'] / region_sizes['control']

    table_35 = [[obs_35, n_genes * region_sizes['-35_region'] - obs_35],
                [obs_ctrl, n_genes * region_sizes['control'] - obs_ctrl]]
    if obs_35 > 0 and obs_ctrl > 0:
        chi2_35, p_35, _, _ = chi2_contingency(table_35)
        print(f"\n  χ² test (-35 vs control): χ²={chi2_35:.2f}, p={p_35:.2e}")

    obs_10 = len(methyl_at_10)
    table_10 = [[obs_10, n_genes * region_sizes['-10_region'] - obs_10],
                [obs_ctrl, n_genes * region_sizes['control'] - obs_ctrl]]
    if obs_10 > 0 and obs_ctrl > 0:
        chi2_10, p_10, _, _ = chi2_contingency(table_10)
        print(f"  χ² test (-10 vs control): χ²={chi2_10:.2f}, p={p_10:.2e}")

    # ── Figure: methylation density across promoter with motif regions highlighted ──
    fig, ax = plt.subplots(figsize=(10, 5))
    positions = methyl_rel_df['rel_pos'].values
    ax.hist(positions, bins=np.arange(-100, 12, 2), color='steelblue', edgecolor='none', alpha=0.7, density=False)
    ax.axvspan(-45, -25, alpha=0.2, color='purple', label='-35 box region')
    ax.axvspan(-20, 0, alpha=0.2, color='orange', label='-10 box region')
    ax.axvline(0, color='green', ls='-', lw=1.5, label='TSS')
    ax.set_xlabel('Distance from TSS (bp)')
    ax.set_ylabel('Number of methylation sites')
    ax.set_title('Methylation site distribution in core promoter region')
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / 'D1_promoter_methylation_distribution.png')
    fig.savefig(OUTPUT_DIR / 'D1_promoter_methylation_distribution.pdf')
    plt.close(fig)

    # ── Figure: Overlap by mod_type and motif_type ──
    if len(overlap_df) > 0:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        for idx, mod in enumerate(['6mA', '4mC']):
            ax = axes[idx]
            sub = overlap_df[overlap_df['mod_type'] == mod]
            if len(sub) > 0:
                counts = sub.groupby(['motif_type', 'overlap_type']).size().unstack(fill_value=0)
                counts.plot(kind='bar', ax=ax, color=['#d6604d', '#2166ac'])
                ax.set_title(f'{mod} at σ70 motifs')
                ax.set_ylabel('Count')
                ax.set_xlabel('')
                ax.legend(title='Overlap type')
            ax.tick_params(axis='x', rotation=0)
        fig.suptitle('Methylation sites overlapping σ70 binding motifs', fontsize=13, y=1.02)
        fig.tight_layout()
        fig.savefig(OUTPUT_DIR / 'D1_sigma_overlap_barplot.png')
        fig.savefig(OUTPUT_DIR / 'D1_sigma_overlap_barplot.pdf')
        plt.close(fig)

    # Save motif data
    motif_df.to_csv(OUTPUT_DIR / 'D1_sigma_motifs.csv', index=False)

    # Enrichment summary
    enrichment = {
        'n_methyl_at_35': len(methyl_at_35),
        'n_methyl_at_10': len(methyl_at_10),
        'n_methyl_ctrl': len(methyl_ctrl),
        'density_35': density_35,
        'density_10': density_10,
        'density_ctrl': density_ctrl,
    }

    print(f"\n  Figures saved to {OUTPUT_DIR}/D1_*")
    return overlap_df, motif_df, enrichment


# ============================================================
# Main
# ============================================================

def main():
    print("TSS-based Epigenome-Transcriptome Integration Analyses")
    print("S. coelicolor A3(2) M145")
    print("=" * 60)

    # A1
    genes_df, jeong_df, a1_summary = run_a1()

    # B1
    b1_results, b1_bins = run_b1(genes_df)

    # C1
    c1_corr = run_c1(genes_df)

    # C2
    c2_corr = run_c2(genes_df)

    # D1
    d1_overlap, d1_motifs, d1_enrichment = run_d1(genes_df)

    print("\n" + "=" * 60)
    print("All analyses complete.")
    print(f"Output directory: {OUTPUT_DIR}")
    print("=" * 60)


if __name__ == '__main__':
    main()
