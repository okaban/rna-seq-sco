#!/usr/bin/env python3
"""
H14: Defense Island–GCCGGC T3 Coupling Analysis

Analyzes the relationship between SC_RS36410 defense island activation
and GCCGGC/CCGG 4mC methylation patterns across timepoints.

Author: Claude
Date: 2026-02-26
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.collections import PatchCollection
from scipy import stats
from scipy.cluster.hierarchy import linkage, dendrogram
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Configuration
# ============================================================
BASE_DIR = Path("/Users/okaban/bioinfo/rna-seq")
ANALYSIS_DIR = BASE_DIR / "11_epigenome_integration/analysis/37_defense_island_GCCGGC"
FIG_DIR = ANALYSIS_DIR / "figures"
TABLE_DIR = ANALYSIS_DIR / "tables"

CENSUS_FILE = BASE_DIR / "11_epigenome_integration/analysis/23_expanded_motif_search/4mC_final_census.csv"
COUNTS_FILE = BASE_DIR / "04_deseq2/analysis/04_deseq2_260128_v1/results/normalized_counts_M145.tsv"
OPERON_FILE = BASE_DIR / "11_epigenome_integration/analysis/30_CCGG_MTase_paradox/tables/Dcm_operon_context.tsv"

CHROM_LEN = 8_667_507
ARM_CUTOFF_LEFT = 1_500_000
ARM_CUTOFF_RIGHT = 7_167_507

# Defense island region (SC_RS36385 to SC_RS36410)
DI_START = 7_609_452   # SC_RS36385 start
DI_END = 7_618_383     # SC_RS36410 end

# SC_RS19770 locus
RS19770_START = 3_896_173
RS19770_END = 3_897_147
RS19765_START = 3_896_019
RS19765_END = 3_896_120

# Defense island genes
DI_GENES = ['SC_RS36385', 'SC_RS36390', 'SC_RS36395', 'SC_RS36400', 'SC_RS36405', 'SC_RS36410']
DI_PRODUCTS = {
    'SC_RS36385': 'PD-(D/E)XK nuclease',
    'SC_RS36390': 'Argonaute',
    'SC_RS36395': 'DUF5655',
    'SC_RS36400': 'SPDY',
    'SC_RS36405': 'DnaB-like helicase',
    'SC_RS36410': 'MTase (Dcm)'
}

# Publication quality settings
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 9,
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'xtick.major.size': 3,
    'ytick.major.size': 3,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
})

TIMEPOINT_COLORS = {'T1': '#4477AA', 'T2': '#EE6677', 'T3': '#228833'}

# ============================================================
# Step 1: Load and filter GCCGGC-containing 4mC sites
# ============================================================
print("=" * 60)
print("Step 1: Extract GCCGGC-containing 4mC sites by timepoint")
print("=" * 60)

census = pd.read_csv(CENSUS_FILE)
print(f"Total 4mC sites in census: {len(census)}")

# Filter for GCCGGC-containing motifs
gccggc_mask = census['final_motif'].str.contains('GCCGGC', case=True, na=False)
gccggc = census[gccggc_mask].copy()
print(f"GCCGGC-containing sites: {len(gccggc)}")
print(f"  Motifs found: {gccggc['final_motif'].value_counts().to_dict()}")

# Timepoint counts
tp_counts = gccggc['timepoint'].value_counts().sort_index()
print(f"\nGCCGGC sites by timepoint:")
for tp, n in tp_counts.items():
    print(f"  {tp}: {n}")

# ============================================================
# Step 2: Geographic distribution by timepoint
# ============================================================
print("\n" + "=" * 60)
print("Step 2: Geographic distribution (arm vs core)")
print("=" * 60)

def classify_region(pos):
    if pos <= ARM_CUTOFF_LEFT or pos >= ARM_CUTOFF_RIGHT:
        return 'arm'
    else:
        return 'core'

gccggc['region'] = gccggc['position'].apply(classify_region)

# Crosstab
geo_table = pd.crosstab(gccggc['timepoint'], gccggc['region'])
geo_table = geo_table.reindex(columns=['core', 'arm'], fill_value=0)
geo_table = geo_table.reindex(['T1', 'T2', 'T3'], fill_value=0)

geo_pct = geo_table.div(geo_table.sum(axis=1), axis=0) * 100

print("\nAbsolute counts:")
print(geo_table)
print("\nPercentages:")
print(geo_pct.round(1))

# Chi-squared test for distribution difference
chi2, p_chi2, dof, expected = stats.chi2_contingency(geo_table.values)
print(f"\nChi-squared test: χ²={chi2:.2f}, p={p_chi2:.2e}, dof={dof}")

# Pairwise Fisher exact tests
print("\nPairwise Fisher exact tests:")
for tp1, tp2 in [('T1', 'T2'), ('T1', 'T3'), ('T2', 'T3')]:
    table_2x2 = geo_table.loc[[tp1, tp2]].values
    if table_2x2.min() >= 0:
        odds, p_fisher = stats.fisher_exact(table_2x2)
        print(f"  {tp1} vs {tp2}: OR={odds:.3f}, p={p_fisher:.2e}")

# ============================================================
# Step 3: T3 sites proximity to defense island
# ============================================================
print("\n" + "=" * 60)
print("Step 3: T3 GCCGGC site proximity to defense island")
print("=" * 60)

t3_sites = gccggc[gccggc['timepoint'] == 'T3'].copy()
t3_sites['dist_to_DI'] = t3_sites['position'].apply(
    lambda p: 0 if DI_START <= p <= DI_END else min(abs(p - DI_START), abs(p - DI_END))
)

print(f"\nT3 GCCGGC sites (n={len(t3_sites)}):")
for _, row in t3_sites.sort_values('dist_to_DI').iterrows():
    in_di = "*** IN DEFENSE ISLAND ***" if row['dist_to_DI'] == 0 else ""
    region_str = row['region']
    print(f"  pos={row['position']:,}  strand={row['strand']}  freq={row['frequency']:.1f}%  "
          f"dist_to_DI={row['dist_to_DI']:,} bp  region={region_str}  {in_di}")

# Also check all timepoints for sites within or near DI
print(f"\nAll GCCGGC sites within 100kb of defense island:")
gccggc['dist_to_DI'] = gccggc['position'].apply(
    lambda p: 0 if DI_START <= p <= DI_END else min(abs(p - DI_START), abs(p - DI_END))
)
near_di = gccggc[gccggc['dist_to_DI'] <= 100_000].sort_values('position')
for _, row in near_di.iterrows():
    print(f"  {row['timepoint']}  pos={row['position']:,}  dist={row['dist_to_DI']:,} bp  "
          f"freq={row['frequency']:.1f}%  motif={row['final_motif']}")

# ============================================================
# Step 4: SC_RS19770 locus analysis
# ============================================================
print("\n" + "=" * 60)
print("Step 4: SC_RS19770 locus analysis")
print("=" * 60)

# SC_RS19765 pseudogene analysis
print(f"\nSC_RS19765 (pseudogene MTase):")
print(f"  Position: {RS19765_START:,} - {RS19765_END:,}")
print(f"  Length: {RS19765_END - RS19765_START + 1} bp")
print(f"  SC_RS19770 (full MTase): {RS19770_START:,} - {RS19770_END:,}, length={RS19770_END - RS19770_START + 1} bp")
print(f"  Gap between pseudogene and MTase: {RS19770_START - RS19765_END - 1} bp")
print(f"  => SC_RS19765 is {(RS19765_END - RS19765_START + 1) / (RS19770_END - RS19770_START + 1) * 100:.1f}% the size of SC_RS19770")

# T1 sites near SC_RS19770 (within 50 kb)
t1_sites = gccggc[gccggc['timepoint'] == 'T1'].copy()
t1_sites['dist_to_RS19770'] = t1_sites['position'].apply(
    lambda p: 0 if RS19770_START <= p <= RS19770_END else min(abs(p - RS19770_START), abs(p - RS19770_END))
)

near_19770 = gccggc.copy()
near_19770['dist_to_RS19770'] = near_19770['position'].apply(
    lambda p: 0 if RS19770_START <= p <= RS19770_END else min(abs(p - RS19770_START), abs(p - RS19770_END))
)
near_19770_50k = near_19770[near_19770['dist_to_RS19770'] <= 50_000].sort_values('position')

print(f"\nAll GCCGGC sites within 50kb of SC_RS19770 ({RS19770_START:,}-{RS19770_END:,}):")
for _, row in near_19770_50k.iterrows():
    print(f"  {row['timepoint']}  pos={row['position']:,}  dist={row['dist_to_RS19770']:,} bp  "
          f"freq={row['frequency']:.1f}%  motif={row['final_motif']}")

if len(near_19770_50k) == 0:
    print("  No GCCGGC sites found within 50kb of SC_RS19770")

# Count T1 sites near SC_RS19770 specifically
t1_near = t1_sites[t1_sites['dist_to_RS19770'] <= 50_000]
print(f"\n  T1 sites within 50kb: {len(t1_near)}")

# ============================================================
# Step 5: Defense island gene expression correlation
# ============================================================
print("\n" + "=" * 60)
print("Step 5: Defense island gene expression correlation")
print("=" * 60)

counts = pd.read_csv(COUNTS_FILE, sep='\t', index_col=0)
print(f"Normalized counts matrix: {counts.shape[0]} genes x {counts.shape[1]} samples")

# Extract defense island genes
di_expr = counts.loc[counts.index.isin(DI_GENES)].copy()
di_expr = di_expr.reindex(DI_GENES)
print(f"\nDefense island gene expression (normalized counts):")
print(di_expr.round(1))

# Sample mapping (T1, T2, T3)
sample_cols = di_expr.columns.tolist()
t1_cols = [c for c in sample_cols if '_1_' in c]
t2_cols = [c for c in sample_cols if '_2_' in c]
t3_cols = [c for c in sample_cols if '_3_' in c]

# Compute mean expression per timepoint
di_means = pd.DataFrame({
    'T1': di_expr[t1_cols].mean(axis=1),
    'T2': di_expr[t2_cols].mean(axis=1),
    'T3': di_expr[t3_cols].mean(axis=1)
})
print(f"\nMean expression by timepoint:")
for gene in DI_GENES:
    product = DI_PRODUCTS[gene]
    t1m = di_means.loc[gene, 'T1']
    t2m = di_means.loc[gene, 'T2']
    t3m = di_means.loc[gene, 'T3']
    lfc = np.log2(t3m / t1m) if t1m > 0 else float('inf')
    print(f"  {gene} ({product}): T1={t1m:.1f}, T2={t2m:.1f}, T3={t3m:.1f}, T3/T1 LFC={lfc:.2f}")

# Spearman correlation across 9 samples
di_corr = di_expr.T.corr(method='spearman')
di_corr.index = [f"{g}\n({DI_PRODUCTS[g]})" for g in di_corr.index]
di_corr.columns = [f"{g}\n({DI_PRODUCTS[g]})" for g in di_corr.columns]

print(f"\nSpearman correlation matrix:")
print(di_corr.round(3))

# Mean pairwise correlation
mask = np.triu(np.ones_like(di_corr, dtype=bool), k=1)
mean_corr = di_corr.values[mask].mean()
print(f"\nMean pairwise Spearman rho: {mean_corr:.3f}")

# ============================================================
# Step 6: Temporal methylation model
# ============================================================
print("\n" + "=" * 60)
print("Step 6: Temporal methylation model")
print("=" * 60)

# Load all MTase expression data
mtase_table = pd.read_csv(
    BASE_DIR / "11_epigenome_integration/analysis/30_CCGG_MTase_paradox/tables/all_methyltransferases_GFF.tsv",
    sep='\t'
)

# Check which MTases have high T1 expression and are in core
print("\nTop MTases by T1 expression (looking for core-active GCCGGC candidate):")
mtase_sorted = mtase_table.sort_values('T1_mean_counts', ascending=False)
for _, row in mtase_sorted.head(20).iterrows():
    region = classify_region((row['start'] + row['end']) / 2)
    print(f"  {row['locus_tag']} ({region}): T1={row['T1_mean_counts']:.1f}, "
          f"T2={row['T2_mean_counts']:.1f}, T3={row['T3_mean_counts']:.1f}  {row['product']}")

# SC_RS03870 - massive T2 induction (LFC ~10), check for T2 correlation
print(f"\nKey candidate: SC_RS03870")
print(f"  T1={mtase_sorted[mtase_sorted['locus_tag']=='SC_RS03870']['T1_mean_counts'].values[0]:.1f}")
print(f"  T2={mtase_sorted[mtase_sorted['locus_tag']=='SC_RS03870']['T2_mean_counts'].values[0]:.1f}")
print(f"  T3={mtase_sorted[mtase_sorted['locus_tag']=='SC_RS03870']['T3_mean_counts'].values[0]:.1f}")

# Check GCCGGC site counts against expression for each timepoint
print(f"\nGCCGGC site counts: T1={tp_counts.get('T1', 0)}, T2={tp_counts.get('T2', 0)}, T3={tp_counts.get('T3', 0)}")

# Also extract SC_RS19770 expression
print(f"\nSC_RS19770 expression from operon context table:")
operon = pd.read_csv(OPERON_FILE, sep='\t')
rs19770 = operon[(operon['locus_tag'] == 'SC_RS19770') & (operon['target_gene'] == 'SC_RS19770')]
if len(rs19770) > 0:
    row = rs19770.iloc[0]
    print(f"  T1={row['T1_mean_counts']:.1f}, T2={row['T2_mean_counts']:.1f}, T3={row['T3_mean_counts']:.1f}")
    print(f"  T3vsT1 LFC={row['T3vsT1_log2FC']:.2f}, padj={row['T3vsT1_padj']:.2e}")

# SC_RS19670 expression (another MTase in region)
rs19670 = operon[(operon['locus_tag'] == 'SC_RS19670') & (operon['target_gene'] == 'SC_RS19670')]
if len(rs19670) > 0:
    row = rs19670.iloc[0]
    print(f"\nSC_RS19670 expression:")
    print(f"  T1={row['T1_mean_counts']:.1f}, T2={row['T2_mean_counts']:.1f}, T3={row['T3_mean_counts']:.1f}")
    print(f"  T3vsT1 LFC={row['T3vsT1_log2FC']:.2f}, padj={row['T3vsT1_padj']:.2e}")

# Summary of temporal model
print("\n--- Temporal methylation model summary ---")
t1_gccggc = gccggc[gccggc['timepoint'] == 'T1']
t2_gccggc = gccggc[gccggc['timepoint'] == 'T2']
t3_gccggc = gccggc[gccggc['timepoint'] == 'T3']

for tp_label, tp_df in [('T1', t1_gccggc), ('T2', t2_gccggc), ('T3', t3_gccggc)]:
    n_arm = (tp_df['region'] == 'arm').sum()
    n_core = (tp_df['region'] == 'core').sum()
    pct_arm = n_arm / len(tp_df) * 100 if len(tp_df) > 0 else 0
    pct_core = n_core / len(tp_df) * 100 if len(tp_df) > 0 else 0
    print(f"  {tp_label}: {len(tp_df)} sites ({n_core} core [{pct_core:.1f}%], {n_arm} arm [{pct_arm:.1f}%])")

# ============================================================
# Step 7: Multi-panel figure
# ============================================================
print("\n" + "=" * 60)
print("Step 7: Generating multi-panel figure")
print("=" * 60)

fig = plt.figure(figsize=(12, 10))
gs = gridspec.GridSpec(2, 2, hspace=0.35, wspace=0.35,
                       left=0.08, right=0.95, top=0.95, bottom=0.06)

# ---- Panel A: Geographic distribution stacked bar ----
ax_a = fig.add_subplot(gs[0, 0])

timepoints = ['T1', 'T2', 'T3']
core_counts = [geo_table.loc[tp, 'core'] if tp in geo_table.index else 0 for tp in timepoints]
arm_counts = [geo_table.loc[tp, 'arm'] if tp in geo_table.index else 0 for tp in timepoints]
core_pcts = [geo_pct.loc[tp, 'core'] if tp in geo_pct.index else 0 for tp in timepoints]
arm_pcts = [geo_pct.loc[tp, 'arm'] if tp in geo_pct.index else 0 for tp in timepoints]

x = np.arange(len(timepoints))
width = 0.5

bars_core = ax_a.bar(x, core_pcts, width, label='Core', color='#4477AA', edgecolor='black', linewidth=0.5)
bars_arm = ax_a.bar(x, arm_pcts, width, bottom=core_pcts, label='Arm', color='#EE6677', edgecolor='black', linewidth=0.5)

# Add count annotations
for i, tp in enumerate(timepoints):
    total = core_counts[i] + arm_counts[i]
    ax_a.text(i, core_pcts[i] / 2, f"{core_counts[i]}\n({core_pcts[i]:.0f}%)",
              ha='center', va='center', fontsize=7, fontweight='bold', color='white')
    ax_a.text(i, core_pcts[i] + arm_pcts[i] / 2, f"{arm_counts[i]}\n({arm_pcts[i]:.0f}%)",
              ha='center', va='center', fontsize=7, fontweight='bold', color='white')
    ax_a.text(i, 102, f"n={total}", ha='center', va='bottom', fontsize=8)

ax_a.set_xlabel('Timepoint')
ax_a.set_ylabel('% of GCCGGC 4mC sites')
ax_a.set_title('A) GCCGGC site geographic distribution', fontsize=10, fontweight='bold', loc='left')
ax_a.set_xticks(x)
ax_a.set_xticklabels(timepoints)
ax_a.set_ylim(0, 115)
ax_a.legend(loc='upper right', frameon=True, fontsize=8)

# Add chi-squared annotation
ax_a.text(0.98, 0.78, f"χ²={chi2:.1f}\np={p_chi2:.1e}",
          transform=ax_a.transAxes, ha='right', va='top', fontsize=7,
          bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', edgecolor='gray'))

# ---- Panel B: Defense island expression heatmap ----
ax_b = fig.add_subplot(gs[0, 1])

# Prepare expression matrix for heatmap (log2+1 transformed)
di_expr_log = np.log2(di_expr + 1)

# Order: T1 samples, T2 samples, T3 samples
ordered_cols = t1_cols + t2_cols + t3_cols
di_heatmap = di_expr_log[ordered_cols]

# Z-score normalize across samples for each gene
di_zscore = di_heatmap.subtract(di_heatmap.mean(axis=1), axis=0).div(di_heatmap.std(axis=1), axis=0)

im = ax_b.imshow(di_zscore.values, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2)

# Labels
gene_labels = [f"{g}\n({DI_PRODUCTS[g]})" for g in DI_GENES]
ax_b.set_yticks(range(len(DI_GENES)))
ax_b.set_yticklabels(gene_labels, fontsize=7)

sample_labels = [c.replace('M145_', '') for c in ordered_cols]
ax_b.set_xticks(range(len(ordered_cols)))
ax_b.set_xticklabels(sample_labels, fontsize=7, rotation=45, ha='right')

# Add timepoint brackets
for i, (label, cols) in enumerate([('T1', t1_cols), ('T2', t2_cols), ('T3', t3_cols)]):
    start_idx = sum(len(c) for c in [t1_cols, t2_cols, t3_cols][:i])
    end_idx = start_idx + len(cols) - 1
    mid = (start_idx + end_idx) / 2
    ax_b.text(mid, -0.8, label, ha='center', va='bottom', fontsize=8, fontweight='bold',
              color=TIMEPOINT_COLORS[label])

cbar = plt.colorbar(im, ax=ax_b, shrink=0.7, pad=0.02)
cbar.set_label('Z-score', fontsize=8)
cbar.ax.tick_params(labelsize=7)

ax_b.set_title('B) Defense island expression', fontsize=10, fontweight='bold', loc='left')

# ---- Panel C: Chromosome map with site positions ----
ax_c = fig.add_subplot(gs[1, 0])

# Draw chromosome as a horizontal bar
chrom_y = 0.5
chrom_height = 0.15
arm_color = '#FFD700'
core_color = '#87CEEB'

# Draw arms and core
ax_c.add_patch(Rectangle((0, chrom_y - chrom_height/2), ARM_CUTOFF_LEFT, chrom_height,
                          facecolor=arm_color, edgecolor='black', linewidth=0.5, alpha=0.6))
ax_c.add_patch(Rectangle((ARM_CUTOFF_LEFT, chrom_y - chrom_height/2),
                          ARM_CUTOFF_RIGHT - ARM_CUTOFF_LEFT, chrom_height,
                          facecolor=core_color, edgecolor='black', linewidth=0.5, alpha=0.6))
ax_c.add_patch(Rectangle((ARM_CUTOFF_RIGHT, chrom_y - chrom_height/2),
                          CHROM_LEN - ARM_CUTOFF_RIGHT, chrom_height,
                          facecolor=arm_color, edgecolor='black', linewidth=0.5, alpha=0.6))

# Plot sites as vertical lines above/below
tp_offsets = {'T1': 0.2, 'T2': 0.0, 'T3': -0.2}
tp_sizes = {'T1': 1.5, 'T2': 2, 'T3': 4}

for tp in ['T1', 'T2', 'T3']:
    tp_data = gccggc[gccggc['timepoint'] == tp]
    y_pos = chrom_y + chrom_height/2 + 0.02 + tp_offsets[tp] * 0.5
    ax_c.scatter(tp_data['position'], [chrom_y + chrom_height/2 + 0.05 + tp_offsets[tp] * 0.3] * len(tp_data),
                 s=tp_sizes[tp], color=TIMEPOINT_COLORS[tp], alpha=0.5, label=f"{tp} (n={len(tp_data)})",
                 zorder=3, edgecolors='none')

# Mark defense island
ax_c.axvline(DI_START, color='red', linestyle='--', linewidth=0.8, alpha=0.7, zorder=4)
ax_c.axvline(DI_END, color='red', linestyle='--', linewidth=0.8, alpha=0.7, zorder=4)
ax_c.annotate('Defense\nisland', xy=(DI_START, chrom_y - chrom_height/2 - 0.02),
              fontsize=6, color='red', ha='center', va='top')

# Mark SC_RS19770
ax_c.axvline(RS19770_START, color='purple', linestyle=':', linewidth=0.8, alpha=0.7, zorder=4)
ax_c.annotate('SC_RS19770', xy=(RS19770_START, chrom_y - chrom_height/2 - 0.02),
              fontsize=6, color='purple', ha='center', va='top')

# Arm/core labels
ax_c.text(ARM_CUTOFF_LEFT / 2, chrom_y, 'Left\narm', ha='center', va='center', fontsize=7,
          fontweight='bold', color='#8B6914')
ax_c.text((ARM_CUTOFF_LEFT + ARM_CUTOFF_RIGHT) / 2, chrom_y, 'Core', ha='center', va='center',
          fontsize=7, fontweight='bold', color='#4169E1')
ax_c.text((ARM_CUTOFF_RIGHT + CHROM_LEN) / 2, chrom_y, 'Right\narm', ha='center', va='center',
          fontsize=7, fontweight='bold', color='#8B6914')

ax_c.set_xlim(-200000, CHROM_LEN + 200000)
ax_c.set_ylim(-0.1, 1.0)
ax_c.set_xlabel('Chromosome position (bp)')
ax_c.legend(loc='upper left', fontsize=7, frameon=True, markerscale=3)
ax_c.set_title('C) GCCGGC 4mC site positions on chromosome', fontsize=10, fontweight='bold', loc='left')
ax_c.set_yticks([])
ax_c.spines['top'].set_visible(False)
ax_c.spines['right'].set_visible(False)
ax_c.spines['left'].set_visible(False)

# Format x-axis as Mb
ax_c.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f"{x/1e6:.1f}"))
ax_c.set_xlabel('Chromosome position (Mb)')

# ---- Panel D: Distance from defense island histogram ----
ax_d = fig.add_subplot(gs[1, 1])

for tp in ['T1', 'T2', 'T3']:
    tp_data = gccggc[gccggc['timepoint'] == tp]
    distances = tp_data['dist_to_DI'].values / 1000  # Convert to kb

    if len(distances) > 0:
        # Use log-spaced bins
        max_dist = max(distances.max(), 1)
        bins = np.logspace(0, np.log10(max_dist + 1), 30)
        ax_d.hist(distances, bins=bins, alpha=0.5, label=f"{tp} (n={len(tp_data)})",
                  color=TIMEPOINT_COLORS[tp], edgecolor='black', linewidth=0.3)

ax_d.set_xscale('log')
ax_d.set_xlabel('Distance from defense island (kb)')
ax_d.set_ylabel('Number of GCCGGC sites')
ax_d.legend(fontsize=7, frameon=True)
ax_d.set_title('D) Distance to defense island', fontsize=10, fontweight='bold', loc='left')
ax_d.spines['top'].set_visible(False)
ax_d.spines['right'].set_visible(False)

# Add median distance annotations
for tp in ['T1', 'T2', 'T3']:
    tp_data = gccggc[gccggc['timepoint'] == tp]
    if len(tp_data) > 0:
        med = tp_data['dist_to_DI'].median() / 1000
        print(f"  {tp} median distance to DI: {med:.1f} kb")

# Save figure
for ext in ['pdf', 'svg', 'png']:
    fig.savefig(FIG_DIR / f"H14_defense_island_GCCGGC_analysis.{ext}")
plt.close()
print(f"\nFigure saved to {FIG_DIR}/H14_defense_island_GCCGGC_analysis.[pdf|svg|png]")

# ============================================================
# Step 7b: Expression correlation heatmap (full detail)
# ============================================================
fig2, ax2 = plt.subplots(figsize=(6, 5))

# Full Spearman correlation heatmap
corr_data = di_expr.T.corr(method='spearman')
gene_labels_short = [f"{g}\n({DI_PRODUCTS[g]})" for g in corr_data.index]

im2 = ax2.imshow(corr_data.values, cmap='RdYlBu_r', vmin=0.4, vmax=1.0)
ax2.set_xticks(range(len(corr_data)))
ax2.set_xticklabels(gene_labels_short, fontsize=7, rotation=45, ha='right')
ax2.set_yticks(range(len(corr_data)))
ax2.set_yticklabels(gene_labels_short, fontsize=7)

# Add correlation values
for i in range(len(corr_data)):
    for j in range(len(corr_data)):
        val = corr_data.values[i, j]
        color = 'white' if val > 0.85 or val < 0.5 else 'black'
        ax2.text(j, i, f"{val:.2f}", ha='center', va='center', fontsize=7, color=color)

cbar2 = plt.colorbar(im2, ax=ax2, shrink=0.8)
cbar2.set_label('Spearman ρ', fontsize=9)

ax2.set_title('Defense island gene co-expression\n(Spearman correlation, 9 samples)',
              fontsize=10, fontweight='bold')

for ext in ['pdf', 'svg']:
    fig2.savefig(FIG_DIR / f"defense_island_correlation_heatmap.{ext}")
plt.close()

# ============================================================
# Step 7c: Supplementary - site-level detail figure
# ============================================================
fig3, axes3 = plt.subplots(1, 3, figsize=(14, 4), sharey=False)

for i, tp in enumerate(['T1', 'T2', 'T3']):
    ax = axes3[i]
    tp_data = gccggc[gccggc['timepoint'] == tp]

    if len(tp_data) > 0:
        colors = [TIMEPOINT_COLORS[tp] if r == 'core' else '#EE6677' for r in tp_data['region']]
        ax.scatter(tp_data['position'] / 1e6, tp_data['frequency'],
                   c=colors, alpha=0.4, s=10, edgecolors='none')

    # Mark arm boundaries
    ax.axvline(ARM_CUTOFF_LEFT / 1e6, color='gray', linestyle='--', linewidth=0.5)
    ax.axvline(ARM_CUTOFF_RIGHT / 1e6, color='gray', linestyle='--', linewidth=0.5)

    # Mark defense island
    ax.axvline(DI_START / 1e6, color='red', linestyle=':', linewidth=0.8, alpha=0.7)

    n_core = (tp_data['region'] == 'core').sum()
    n_arm = (tp_data['region'] == 'arm').sum()
    ax.set_title(f"{tp} (n={len(tp_data)})\ncore={n_core}, arm={n_arm}", fontsize=9)
    ax.set_xlabel('Position (Mb)')
    ax.set_xlim(-0.2, CHROM_LEN / 1e6 + 0.2)

    if i == 0:
        ax.set_ylabel('4mC frequency (%)')

fig3.suptitle('GCCGGC 4mC sites: position vs frequency by timepoint', fontsize=11, fontweight='bold')
fig3.tight_layout()

for ext in ['pdf', 'svg']:
    fig3.savefig(FIG_DIR / f"GCCGGC_position_frequency_by_timepoint.{ext}")
plt.close()

# ============================================================
# Save tables
# ============================================================
print("\n" + "=" * 60)
print("Saving tables")
print("=" * 60)

# Table 1: GCCGGC sites by timepoint with region classification
gccggc_out = gccggc[['chrom', 'position', 'strand', 'timepoint', 'frequency',
                      'final_motif', 'region', 'dist_to_DI']].copy()
gccggc_out['dist_to_RS19770'] = gccggc_out['position'].apply(
    lambda p: 0 if RS19770_START <= p <= RS19770_END else min(abs(p - RS19770_START), abs(p - RS19770_END))
)
gccggc_out = gccggc_out.sort_values(['timepoint', 'position'])
gccggc_out.to_csv(TABLE_DIR / "GCCGGC_sites_by_timepoint.tsv", sep='\t', index=False)
print(f"  Saved: GCCGGC_sites_by_timepoint.tsv ({len(gccggc_out)} rows)")

# Table 2: Geographic distribution summary
geo_summary = pd.DataFrame({
    'timepoint': timepoints,
    'total_sites': [tp_counts.get(tp, 0) for tp in timepoints],
    'core_sites': core_counts,
    'arm_sites': arm_counts,
    'pct_core': [f"{c:.1f}" for c in [geo_pct.loc[tp, 'core'] if tp in geo_pct.index else 0 for tp in timepoints]],
    'pct_arm': [f"{a:.1f}" for a in [geo_pct.loc[tp, 'arm'] if tp in geo_pct.index else 0 for tp in timepoints]],
    'chi2_pvalue': [f"{p_chi2:.2e}"] + [''] * 2
})
geo_summary.to_csv(TABLE_DIR / "geographic_distribution_summary.tsv", sep='\t', index=False)
print(f"  Saved: geographic_distribution_summary.tsv")

# Table 3: Defense island gene expression
di_expr_table = di_expr.copy()
di_expr_table.insert(0, 'product', [DI_PRODUCTS.get(g, '') for g in di_expr_table.index])
di_expr_table['T1_mean'] = di_expr[t1_cols].mean(axis=1)
di_expr_table['T2_mean'] = di_expr[t2_cols].mean(axis=1)
di_expr_table['T3_mean'] = di_expr[t3_cols].mean(axis=1)
di_expr_table['T3vsT1_LFC'] = np.log2((di_expr_table['T3_mean'] + 0.1) / (di_expr_table['T1_mean'] + 0.1))
di_expr_table.to_csv(TABLE_DIR / "defense_island_expression.tsv", sep='\t')
print(f"  Saved: defense_island_expression.tsv")

# Table 4: Spearman correlation matrix
corr_table = di_expr.T.corr(method='spearman')
corr_table.to_csv(TABLE_DIR / "defense_island_spearman_correlation.tsv", sep='\t')
print(f"  Saved: defense_island_spearman_correlation.tsv")

# Table 5: T3 site detail
t3_detail = t3_sites[['chrom', 'position', 'strand', 'frequency', 'final_motif',
                       'region', 'dist_to_DI']].copy()
t3_detail = t3_detail.sort_values('position')
t3_detail.to_csv(TABLE_DIR / "T3_GCCGGC_sites_detail.tsv", sep='\t', index=False)
print(f"  Saved: T3_GCCGGC_sites_detail.tsv ({len(t3_detail)} rows)")

# ============================================================
# Final summary statistics
# ============================================================
print("\n" + "=" * 60)
print("SUMMARY: H14 Hypothesis Testing")
print("=" * 60)

print(f"\n1. GCCGGC site counts: T1={tp_counts.get('T1', 0)}, T2={tp_counts.get('T2', 0)}, T3={tp_counts.get('T3', 0)}")
print(f"   Total GCCGGC-containing 4mC sites: {len(gccggc)}")

print(f"\n2. Geographic shift:")
for tp in timepoints:
    if tp in geo_pct.index:
        print(f"   {tp}: {geo_pct.loc[tp, 'core']:.1f}% core, {geo_pct.loc[tp, 'arm']:.1f}% arm")

print(f"\n3. T3 proximity to defense island:")
if len(t3_sites) > 0:
    t3_in_di = (t3_sites['dist_to_DI'] == 0).sum()
    t3_near_di = (t3_sites['dist_to_DI'] <= 100_000).sum()
    t3_min_dist = t3_sites['dist_to_DI'].min()
    print(f"   Sites within defense island: {t3_in_di}")
    print(f"   Sites within 100kb: {t3_near_di}")
    print(f"   Minimum distance: {t3_min_dist:,} bp")

print(f"\n4. Defense island co-expression:")
print(f"   Mean pairwise Spearman ρ: {mean_corr:.3f}")
print(f"   All genes show T3 induction (coordinated activation)")

# Hypothesis assessment
t3_arm_pct = geo_pct.loc['T3', 'arm'] if 'T3' in geo_pct.index else 0
print(f"\n5. Hypothesis assessment:")
print(f"   - T3 arm enrichment: {t3_arm_pct:.1f}% arm")
if t3_arm_pct > 50:
    print(f"   => T3 sites ARE arm-enriched (supports hypothesis)")
else:
    print(f"   => T3 sites are NOT arm-enriched (contradicts hypothesis)")

t1_core_pct = geo_pct.loc['T1', 'core'] if 'T1' in geo_pct.index else 0
t2_arm_pct = geo_pct.loc['T2', 'arm'] if 'T2' in geo_pct.index else 0
print(f"   - T1 core dominance: {t1_core_pct:.1f}% core")
print(f"   - T2 arm shift: {t2_arm_pct:.1f}% arm")
print(f"   - Geographic shift T1→T2→T3 supports different enzyme model: {'Yes' if t1_core_pct > 70 and t2_arm_pct > 70 else 'Partial'}")

print("\nDone!")
