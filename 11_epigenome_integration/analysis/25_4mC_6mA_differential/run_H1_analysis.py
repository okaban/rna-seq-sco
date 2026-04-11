#!/usr/bin/env python3
"""
H1 Analysis: Compare genomic distribution patterns of 4mC vs 6mA methylation sites.
S. coelicolor A3(2) M145 RNA-seq + Epigenome Integration
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
from scipy import stats
from collections import defaultdict
import re
import sys
import warnings
warnings.filterwarnings('ignore')

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/25_4mC_6mA_differential"
FIG_DIR = f"{BASE}/figures"
TBL_DIR = f"{BASE}/tables"

SITES_FILE = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/01_integration/high_confidence_sites_weighted.csv"
GFF_FILE = "/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/genomic.gff"
INTEG_FILE = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/01_integration/integrated_methyl_expression_weighted.csv"
DESEQ_FILE = "/Users/okaban/bioinfo/rna-seq/04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_2_vs_1.tsv"

# ---------------------------------------------------------------------------
# Matplotlib settings
# ---------------------------------------------------------------------------
plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
print("=" * 70)
print("H1: 4mC vs 6mA Genomic Distribution Analysis")
print("=" * 70)

sites = pd.read_csv(SITES_FILE)
integ = pd.read_csv(INTEG_FILE)
deseq = pd.read_csv(DESEQ_FILE, sep='\t')

print(f"\nLoaded {len(sites)} methylation sites ({sites['mod_type'].value_counts().to_dict()})")
print(f"Loaded {len(integ)} integrated gene records")
print(f"Loaded {len(deseq)} DESeq2 results")

# ---------------------------------------------------------------------------
# 2. Parse GFF -> gene coordinates
# ---------------------------------------------------------------------------
genes = []
with open(GFF_FILE) as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9:
            continue
        if parts[2] != 'gene':
            continue
        chrom = parts[0]
        start = int(parts[3])
        end = int(parts[4])
        strand = parts[6]
        attrs = parts[8]
        # Extract locus_tag
        m = re.search(r'locus_tag=([^;]+)', attrs)
        locus_tag = m.group(1) if m else None
        if locus_tag:
            genes.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': strand,
                'gene_id': locus_tag
            })

genes_df = pd.DataFrame(genes)
print(f"Parsed {len(genes_df)} genes from GFF")

# ---------------------------------------------------------------------------
# 3. Genomic Region Classification
# ---------------------------------------------------------------------------
PROMOTER_UP = 300   # upstream of TSS
PROMOTER_DOWN = 50  # downstream of TSS

def classify_site(pos, chrom, genes_chrom):
    """Classify a site into genomic region relative to nearest gene."""
    best_region = 'Intergenic'
    best_dist = float('inf')
    best_gene = None

    for _, g in genes_chrom.iterrows():
        gstart, gend, gstrand = g['start'], g['end'], g['strand']
        gene_len = gend - gstart + 1

        # Determine TSS
        tss = gstart if gstrand == '+' else gend

        # Check promoter
        if gstrand == '+':
            prom_start = tss - PROMOTER_UP
            prom_end = tss + PROMOTER_DOWN
        else:
            prom_start = tss - PROMOTER_DOWN
            prom_end = tss + PROMOTER_UP

        if prom_start <= pos <= prom_end:
            dist = abs(pos - tss)
            if dist < best_dist:
                best_dist = dist
                best_region = 'Promoter'
                best_gene = g['gene_id']
            continue

        # Check gene body
        if gstart <= pos <= gend:
            # Determine relative position within gene
            if gstrand == '+':
                rel_pos = (pos - gstart) / gene_len
            else:
                rel_pos = (gend - pos) / gene_len

            if rel_pos < 0.25:
                region = "5' gene body"
            elif rel_pos < 0.75:
                region = "Mid gene body"
            else:
                region = "3' gene body"

            dist = 0  # inside gene
            if dist < best_dist or (dist == best_dist and best_region == 'Intergenic'):
                best_dist = dist
                best_region = region
                best_gene = g['gene_id']

    return best_region, best_gene


# Build sorted gene index per chromosome for faster lookup
# Use interval-based approach for efficiency
def classify_sites_vectorized(sites_df, genes_df):
    """Classify all methylation sites using vectorized operations where possible."""
    results = []
    
    for chrom in sites_df['chrom'].unique():
        sites_chrom = sites_df[sites_df['chrom'] == chrom].copy()
        genes_chrom = genes_df[genes_df['chrom'] == chrom].copy().sort_values('start')
        
        if len(genes_chrom) == 0:
            sites_chrom['region'] = 'Intergenic'
            sites_chrom['nearest_gene'] = None
            results.append(sites_chrom)
            continue
        
        positions = sites_chrom['position'].values
        g_starts = genes_chrom['start'].values
        g_ends = genes_chrom['end'].values
        g_strands = genes_chrom['strand'].values
        g_ids = genes_chrom['gene_id'].values
        
        regions = []
        nearest_genes = []
        
        for pos in positions:
            # Find nearby genes using binary search
            idx = np.searchsorted(g_starts, pos)
            # Check genes in a window around the position
            check_range = range(max(0, idx - 5), min(len(g_starts), idx + 5))
            
            best_region = 'Intergenic'
            best_gene = None
            best_priority = 99  # lower = higher priority
            best_dist = float('inf')
            
            for i in check_range:
                gs, ge, gstr, gid = g_starts[i], g_ends[i], g_strands[i], g_ids[i]
                gene_len = ge - gs + 1
                tss = gs if gstr == '+' else ge
                
                # Check promoter
                if gstr == '+':
                    prom_s = tss - PROMOTER_UP
                    prom_e = tss + PROMOTER_DOWN
                else:
                    prom_s = tss - PROMOTER_DOWN
                    prom_e = tss + PROMOTER_UP
                
                if prom_s <= pos <= prom_e:
                    dist = abs(pos - tss)
                    if 0 < best_priority or dist < best_dist:
                        best_priority = 0
                        best_dist = dist
                        best_region = 'Promoter'
                        best_gene = gid
                    continue
                
                # Check gene body
                if gs <= pos <= ge:
                    if gstr == '+':
                        rel = (pos - gs) / gene_len
                    else:
                        rel = (ge - pos) / gene_len
                    
                    if rel < 0.25:
                        region = "5' gene body"
                    elif rel < 0.75:
                        region = "Mid gene body"
                    else:
                        region = "3' gene body"
                    
                    if best_priority > 1:
                        best_priority = 1
                        best_region = region
                        best_gene = gid
            
            regions.append(best_region)
            nearest_genes.append(best_gene)
        
        sites_chrom['region'] = regions
        sites_chrom['nearest_gene'] = nearest_genes
        results.append(sites_chrom)
    
    return pd.concat(results, ignore_index=True)

print("\nClassifying methylation sites into genomic regions...")
sites_classified = classify_sites_vectorized(sites, genes_df)
print("Done.")

# ---------------------------------------------------------------------------
# 4. Summary table: region distribution by mod_type
# ---------------------------------------------------------------------------
region_order = ['Promoter', "5' gene body", 'Mid gene body', "3' gene body", 'Intergenic']

# Overall distribution
dist_table = sites_classified.groupby(['mod_type', 'region']).size().unstack(fill_value=0)
dist_table = dist_table.reindex(columns=region_order, fill_value=0)
dist_pct = dist_table.div(dist_table.sum(axis=1), axis=0) * 100

print("\n--- Region Distribution (%) ---")
print(dist_pct.round(1).to_string())

# By timepoint
dist_tp = sites_classified.groupby(['mod_type', 'timepoint', 'region']).size().unstack(fill_value=0)
dist_tp = dist_tp.reindex(columns=region_order, fill_value=0)
dist_tp_pct = dist_tp.div(dist_tp.sum(axis=1), axis=0) * 100

print("\n--- Region Distribution by Timepoint (%) ---")
print(dist_tp_pct.round(1).to_string())

# Save tables
dist_table.to_csv(f"{TBL_DIR}/region_distribution_counts.tsv", sep='\t')
dist_pct.round(2).to_csv(f"{TBL_DIR}/region_distribution_pct.tsv", sep='\t')
dist_tp_pct.round(2).to_csv(f"{TBL_DIR}/region_distribution_by_timepoint_pct.tsv", sep='\t')

# ---------------------------------------------------------------------------
# 5. Stacked Bar Plot — Overall 4mC vs 6mA
# ---------------------------------------------------------------------------
colors = {
    'Promoter': '#e74c3c',
    "5' gene body": '#f39c12',
    'Mid gene body': '#2ecc71',
    "3' gene body": '#3498db',
    'Intergenic': '#95a5a6'
}

fig, axes = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'width_ratios': [1, 2]})

# Panel A: Overall comparison
ax = axes[0]
mod_types = ['4mC', '6mA']
bottoms = {mt: 0 for mt in mod_types}
x_pos = np.arange(len(mod_types))

for region in region_order:
    heights = [dist_pct.loc[mt, region] if mt in dist_pct.index else 0 for mt in mod_types]
    bars = ax.bar(x_pos, heights, bottom=[bottoms[mt] for mt in mod_types],
                  color=colors[region], label=region, edgecolor='white', linewidth=0.5)
    # Add percentage labels for regions > 5%
    for i, (h, mt) in enumerate(zip(heights, mod_types)):
        if h > 5:
            ax.text(x_pos[i], bottoms[mt] + h/2, f'{h:.1f}%',
                    ha='center', va='center', fontsize=8, fontweight='bold', color='white')
    for mt, h in zip(mod_types, heights):
        bottoms[mt] += h

ax.set_xticks(x_pos)
ax.set_xticklabels(mod_types, fontsize=12)
ax.set_ylabel('Percentage of sites (%)')
ax.set_title('A) Overall Distribution')
ax.set_ylim(0, 100)
ax.legend(loc='upper right', fontsize=7, framealpha=0.9)

# Panel B: By timepoint
ax = axes[1]
groups = [('4mC', 'T1'), ('4mC', 'T2'), ('4mC', 'T3'),
          ('6mA', 'T1'), ('6mA', 'T2'), ('6mA', 'T3')]
labels = ['4mC\nT1', '4mC\nT2', '4mC\nT3', '6mA\nT1', '6mA\nT2', '6mA\nT3']
x_pos = np.arange(len(groups))
bottoms = [0] * len(groups)

for region in region_order:
    heights = []
    for mt, tp in groups:
        if (mt, tp) in dist_tp_pct.index:
            heights.append(dist_tp_pct.loc[(mt, tp), region])
        else:
            heights.append(0)
    ax.bar(x_pos, heights, bottom=bottoms, color=colors[region], label=region,
           edgecolor='white', linewidth=0.5)
    for i, (h, b) in enumerate(zip(heights, bottoms)):
        if h > 7:
            ax.text(x_pos[i], b + h/2, f'{h:.0f}%',
                    ha='center', va='center', fontsize=7, fontweight='bold', color='white')
    bottoms = [b + h for b, h in zip(bottoms, heights)]

ax.set_xticks(x_pos)
ax.set_xticklabels(labels, fontsize=9)
ax.set_ylabel('Percentage of sites (%)')
ax.set_title('B) Distribution by Timepoint')
ax.set_ylim(0, 100)

# Add separator line between 4mC and 6mA groups
ax.axvline(2.5, color='black', linestyle='--', linewidth=0.8, alpha=0.5)

plt.suptitle('Genomic Distribution of 4mC vs 6mA Methylation Sites', fontsize=13, fontweight='bold', y=1.02)
plt.tight_layout()
fig.savefig(f"{FIG_DIR}/H1_region_distribution_stacked_bar.pdf")
fig.savefig(f"{FIG_DIR}/H1_region_distribution_stacked_bar.svg")
plt.close()
print("\nSaved stacked bar plot.")

# ---------------------------------------------------------------------------
# 6. Region-specific Correlation: methylation change vs expression change
# ---------------------------------------------------------------------------
print("\n--- Region-specific Correlation Analysis ---")

# Merge integrated data with DESeq2
deseq_sub = deseq[['gene_id', 'log2FoldChange', 'padj']].rename(
    columns={'log2FoldChange': 'log2FC_deseq', 'padj': 'padj_deseq'}
)
integ_merged = integ.merge(deseq_sub, on='gene_id', how='left')

# Also classify each gene's methylation sites as promoter or gene body
# Use the site classifications
sites_gene_region = sites_classified[sites_classified['nearest_gene'].notna()].copy()

# For each gene, determine if it has promoter and/or gene body methylation
gene_region_types = sites_gene_region.groupby(['nearest_gene', 'mod_type', 'region']).size().reset_index(name='count')

# Create flags for promoter vs gene body
gene_region_types['is_promoter'] = gene_region_types['region'] == 'Promoter'
gene_region_types['is_genebody'] = gene_region_types['region'].isin(["5' gene body", "Mid gene body", "3' gene body"])

# Build correlation dataset
corr_results = []

for mod_type in ['4mC', '6mA']:
    change_col = f'{mod_type.lower().replace("m","m")}_change_T2_vs_T1'
    # The column names use exact mod_type prefix
    if mod_type == '6mA':
        change_col = '6mA_change_T2_vs_T1'
    else:
        change_col = '4mC_change_T2_vs_T1'
    
    # Get genes with this mod type in promoter
    prom_genes = set(gene_region_types[
        (gene_region_types['mod_type'] == mod_type) & 
        (gene_region_types['is_promoter'])
    ]['nearest_gene'].values)
    
    # Get genes with this mod type in gene body
    gb_genes = set(gene_region_types[
        (gene_region_types['mod_type'] == mod_type) & 
        (gene_region_types['is_genebody'])
    ]['nearest_gene'].values)
    
    for region_name, gene_set in [('Promoter', prom_genes), ('Gene Body', gb_genes)]:
        subset = integ_merged[integ_merged['gene_id'].isin(gene_set)].copy()
        subset = subset.dropna(subset=[change_col, 'log2FC_deseq'])
        
        # Encode methylation change: positive = gained, negative = lost, 0 = stable
        methyl_change = subset[change_col].values
        # Convert to direction: +1, -1, 0
        methyl_dir = np.sign(methyl_change)
        expr_change = subset['log2FC_deseq'].values
        
        if len(subset) >= 10:
            rho, pval = stats.spearmanr(methyl_change, expr_change)
            rho_dir, pval_dir = stats.spearmanr(methyl_dir, expr_change)
        else:
            rho, pval = np.nan, np.nan
            rho_dir, pval_dir = np.nan, np.nan
        
        corr_results.append({
            'mod_type': mod_type,
            'region': region_name,
            'n_genes': len(subset),
            'spearman_rho_continuous': round(rho, 4) if not np.isnan(rho) else np.nan,
            'pval_continuous': pval,
            'spearman_rho_direction': round(rho_dir, 4) if not np.isnan(rho_dir) else np.nan,
            'pval_direction': pval_dir,
        })
        
        sig = '*' if pval < 0.05 else 'ns'
        print(f"  {mod_type} {region_name}: n={len(subset)}, rho={rho:.4f}, p={pval:.2e} [{sig}]")

corr_df = pd.DataFrame(corr_results)
corr_df.to_csv(f"{TBL_DIR}/region_specific_correlation.tsv", sep='\t', index=False)

# ---------------------------------------------------------------------------
# 7. Correlation scatter plots
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(10, 9))

plot_idx = 0
for mod_type in ['4mC', '6mA']:
    change_col = f'{mod_type}_change_T2_vs_T1'
    
    prom_genes = set(gene_region_types[
        (gene_region_types['mod_type'] == mod_type) & 
        (gene_region_types['is_promoter'])
    ]['nearest_gene'].values)
    
    gb_genes = set(gene_region_types[
        (gene_region_types['mod_type'] == mod_type) & 
        (gene_region_types['is_genebody'])
    ]['nearest_gene'].values)
    
    for j, (region_name, gene_set) in enumerate([('Promoter', prom_genes), ('Gene Body', gb_genes)]):
        ax = axes[plot_idx // 2][plot_idx % 2]
        plot_idx += 1
        
        subset = integ_merged[integ_merged['gene_id'].isin(gene_set)].copy()
        subset = subset.dropna(subset=[change_col, 'log2FC_deseq'])
        
        if len(subset) < 5:
            ax.text(0.5, 0.5, 'Insufficient data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'{mod_type} - {region_name}')
            continue
        
        x = subset[change_col].values
        y = subset['log2FC_deseq'].values
        rho, pval = stats.spearmanr(x, y)
        
        ax.scatter(x, y, alpha=0.4, s=15, c='#2c3e50', edgecolors='none')
        
        # Add trend line
        z = np.polyfit(x, y, 1)
        p = np.poly1d(z)
        x_line = np.linspace(x.min(), x.max(), 100)
        ax.plot(x_line, p(x_line), 'r-', linewidth=1.5, alpha=0.7)
        
        sig_str = f'p={pval:.2e}' if pval < 0.001 else f'p={pval:.4f}'
        ax.set_title(f'{mod_type} - {region_name}\n(rho={rho:.3f}, {sig_str}, n={len(subset)})')
        ax.set_xlabel(f'{mod_type} methylation change (T2 vs T1)')
        ax.set_ylabel('log2FC expression (T2 vs T1)')
        ax.axhline(0, color='grey', linewidth=0.5, linestyle='--')
        ax.axvline(0, color='grey', linewidth=0.5, linestyle='--')

plt.suptitle('Methylation Change vs Expression Change by Region', fontsize=13, fontweight='bold')
plt.tight_layout()
fig.savefig(f"{FIG_DIR}/H1_methylation_expression_correlation.pdf")
fig.savefig(f"{FIG_DIR}/H1_methylation_expression_correlation.svg")
plt.close()
print("\nSaved correlation scatter plots.")

# ---------------------------------------------------------------------------
# 8. Fisher's Exact Test: Promoter enrichment 4mC vs 6mA
# ---------------------------------------------------------------------------
print("\n--- Statistical Test: Promoter Enrichment ---")

# Contingency table: rows = mod_type (4mC, 6mA), cols = (Promoter, Non-promoter)
sites_4mC = sites_classified[sites_classified['mod_type'] == '4mC']
sites_6mA = sites_classified[sites_classified['mod_type'] == '6mA']

n_4mC_prom = (sites_4mC['region'] == 'Promoter').sum()
n_4mC_other = len(sites_4mC) - n_4mC_prom
n_6mA_prom = (sites_6mA['region'] == 'Promoter').sum()
n_6mA_other = len(sites_6mA) - n_6mA_prom

contingency = np.array([[n_4mC_prom, n_4mC_other],
                         [n_6mA_prom, n_6mA_other]])

print(f"\nContingency Table (Promoter vs Non-Promoter):")
print(f"  {'':>10s}  {'Promoter':>10s}  {'Non-Prom':>10s}  {'Total':>10s}  {'% Prom':>8s}")
print(f"  {'4mC':>10s}  {n_4mC_prom:>10d}  {n_4mC_other:>10d}  {len(sites_4mC):>10d}  {100*n_4mC_prom/len(sites_4mC):>7.1f}%")
print(f"  {'6mA':>10s}  {n_6mA_prom:>10d}  {n_6mA_other:>10d}  {len(sites_6mA):>10d}  {100*n_6mA_prom/len(sites_6mA):>7.1f}%")

# Fisher's exact test
odds_ratio, fisher_p = stats.fisher_exact(contingency)
print(f"\n  Fisher's exact test: OR = {odds_ratio:.3f}, p = {fisher_p:.2e}")

# Chi-squared test (for completeness)
chi2, chi2_p, dof, expected = stats.chi2_contingency(contingency)
print(f"  Chi-squared test:   chi2 = {chi2:.2f}, p = {chi2_p:.2e}, dof = {dof}")

# Also test each region
print("\n--- Chi-squared Test: Full Region Distribution ---")
full_contingency = dist_table.reindex(['4mC', '6mA']).values
chi2_full, p_full, dof_full, exp_full = stats.chi2_contingency(full_contingency)
print(f"  Chi-squared = {chi2_full:.2f}, p = {p_full:.2e}, dof = {dof_full}")

# Save Fisher test results
fisher_results = pd.DataFrame({
    'test': ['Fisher_exact_promoter', 'Chi2_promoter', 'Chi2_all_regions'],
    'statistic': [odds_ratio, chi2, chi2_full],
    'p_value': [fisher_p, chi2_p, p_full],
    'note': [
        f'4mC promoter: {n_4mC_prom}/{len(sites_4mC)} ({100*n_4mC_prom/len(sites_4mC):.1f}%); 6mA promoter: {n_6mA_prom}/{len(sites_6mA)} ({100*n_6mA_prom/len(sites_6mA):.1f}%)',
        f'dof={dof}',
        f'dof={dof_full}, regions: {region_order}'
    ]
})
fisher_results.to_csv(f"{TBL_DIR}/statistical_tests.tsv", sep='\t', index=False)

# ---------------------------------------------------------------------------
# 9. Supplementary figure: region distribution heatmap
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 4))

# Heatmap of % distribution
heatmap_data = dist_tp_pct.reset_index()
heatmap_data['label'] = heatmap_data['mod_type'] + ' ' + heatmap_data['timepoint']
heatmap_pivot = heatmap_data.set_index('label')[region_order]

sns.heatmap(heatmap_pivot, annot=True, fmt='.1f', cmap='YlOrRd', ax=ax,
            linewidths=0.5, cbar_kws={'label': '% of sites'})
ax.set_title('Genomic Region Distribution (%) by Modification Type and Timepoint')
ax.set_ylabel('')
plt.tight_layout()
fig.savefig(f"{FIG_DIR}/H1_region_distribution_heatmap.pdf")
fig.savefig(f"{FIG_DIR}/H1_region_distribution_heatmap.svg")
plt.close()
print("\nSaved heatmap figure.")

# ---------------------------------------------------------------------------
# 10. Summary
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("KEY FINDINGS SUMMARY")
print("=" * 70)

# Biggest difference between 4mC and 6mA
for region in region_order:
    pct_4mC = dist_pct.loc['4mC', region] if '4mC' in dist_pct.index else 0
    pct_6mA = dist_pct.loc['6mA', region] if '6mA' in dist_pct.index else 0
    diff = pct_4mC - pct_6mA
    if abs(diff) > 1:
        direction = "enriched in 4mC" if diff > 0 else "enriched in 6mA"
        print(f"  {region}: {pct_4mC:.1f}% (4mC) vs {pct_6mA:.1f}% (6mA) -> {direction} by {abs(diff):.1f}pp")

print(f"\n  Promoter enrichment: 4mC={100*n_4mC_prom/len(sites_4mC):.1f}% vs 6mA={100*n_6mA_prom/len(sites_6mA):.1f}%")
print(f"    Fisher's OR={odds_ratio:.3f}, p={fisher_p:.2e}")

print(f"\n  Full region distribution (chi2): p={p_full:.2e}")

print(f"\n  Correlation (methylation change vs expression, T2 vs T1):")
for _, row in corr_df.iterrows():
    sig = '***' if row['pval_continuous'] < 0.001 else ('**' if row['pval_continuous'] < 0.01 else ('*' if row['pval_continuous'] < 0.05 else 'ns'))
    print(f"    {row['mod_type']} {row['region']}: rho={row['spearman_rho_continuous']:.4f}, p={row['pval_continuous']:.2e} [{sig}], n={row['n_genes']}")

print("\n" + "=" * 70)
print("Output files:")
print(f"  Figures: {FIG_DIR}/")
print(f"  Tables:  {TBL_DIR}/")
print("=" * 70)
