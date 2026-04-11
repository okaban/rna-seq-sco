#!/usr/bin/env python3
"""
H8 Analysis: Geographic and functional characterization of 62 coordinated regulatory genes.

Tests the hypothesis that the 62 methylation-expression coordinated regulators are:
1. Enriched at chromosomal arms (vs core)
2. Enriched in stress response / secondary metabolism functional categories
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# CONFIGURATION
# ============================================================
CHROM_LEN = 8_667_507
CORE_START = 1_500_000
CORE_END = 6_500_000

BASE_DIR = '/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/31_coordinated_regulators_characterization'
FIG_DIR = f'{BASE_DIR}/figures'
TBL_DIR = f'{BASE_DIR}/tables'

COORD_FILE = '/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/29_genomewide_TF_screen/tables/coordinated_regulatory_genes.tsv'
ALL_REG_FILE = '/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/29_genomewide_TF_screen/tables/all_regulatory_genes.tsv'
COG_FILE = '/Users/okaban/bioinfo/rna-seq/12_supplementary_figures/analysis/12_supplementary_260202_v1/tables/gene_COG_classification.tsv'
GENE_MASTER_FILE = '/Users/okaban/bioinfo/rna-seq/05_annotation/analysis/05_annotation_260128_v1/tables/gene_master_with_BGC.tsv'
GENE_ANNOT_FILE = '/Users/okaban/bioinfo/rna-seq/05_annotation/analysis/05_annotation_260128_v1/tables/gene_annotation_basic.tsv'

# Consistent style
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

# ============================================================
# DATA LOADING
# ============================================================
print("=" * 70)
print("H8 ANALYSIS: Coordinated Regulatory Genes Characterization")
print("=" * 70)

coord_df = pd.read_csv(COORD_FILE, sep='\t')
all_reg_df = pd.read_csv(ALL_REG_FILE, sep='\t')
cog_df = pd.read_csv(COG_FILE, sep='\t')
gene_master = pd.read_csv(GENE_MASTER_FILE, sep='\t', low_memory=False)
gene_annot = pd.read_csv(GENE_ANNOT_FILE, sep='\t')

print(f"Coordinated regulatory genes: {len(coord_df)}")
print(f"All regulatory genes: {len(all_reg_df)}")
print(f"Genome-wide genes (COG file): {len(cog_df)}")
print(f"Genome-wide genes (master): {len(gene_master)}")

# ============================================================
# HELPER FUNCTIONS
# ============================================================
def classify_region(pos):
    """Classify genomic position as core or arm."""
    if pos < CORE_START or pos > CORE_END:
        return 'arm'
    return 'core'

def midpoint(row):
    return (row['start'] + row['end']) / 2

def fisher_test_arm_enrichment(test_arm, test_core, ref_arm, ref_core):
    """Fisher's exact test for arm enrichment."""
    table = [[test_arm, test_core], [ref_arm, ref_core]]
    odds_ratio, p_value = stats.fisher_exact(table, alternative='greater')
    return odds_ratio, p_value

# ============================================================
# STEP 1: GEOGRAPHIC DISTRIBUTION
# ============================================================
print("\n" + "=" * 70)
print("STEP 1: GEOGRAPHIC DISTRIBUTION")
print("=" * 70)

# Classify 62 coordinated genes
coord_df['midpoint'] = coord_df.apply(midpoint, axis=1)
coord_df['region'] = coord_df['midpoint'].apply(classify_region)

# Classify all 1,055 regulatory genes
all_reg_df['midpoint'] = all_reg_df.apply(midpoint, axis=1)
all_reg_df['region'] = all_reg_df['midpoint'].apply(classify_region)

# Classify all genome-wide genes
gene_annot['midpoint'] = (gene_annot['start'] + gene_annot['end']) / 2
gene_annot['region'] = gene_annot['midpoint'].apply(classify_region)

# Counts
coord_arm = (coord_df['region'] == 'arm').sum()
coord_core = (coord_df['region'] == 'core').sum()
all_reg_arm = (all_reg_df['region'] == 'arm').sum()
all_reg_core = (all_reg_df['region'] == 'core').sum()
genome_arm = (gene_annot['region'] == 'arm').sum()
genome_core = (gene_annot['region'] == 'core').sum()

print(f"\n62 Coordinated regulators: {coord_arm} arm, {coord_core} core ({coord_arm/len(coord_df)*100:.1f}% arm)")
print(f"1,055 All regulators: {all_reg_arm} arm, {all_reg_core} core ({all_reg_arm/len(all_reg_df)*100:.1f}% arm)")
print(f"Genome-wide: {genome_arm} arm, {genome_core} core ({genome_arm/len(gene_annot)*100:.1f}% arm)")

# Arm fraction expected based on chromosome length
arm_fraction = (CORE_START + (CHROM_LEN - CORE_END)) / CHROM_LEN
print(f"\nExpected arm fraction by length: {arm_fraction*100:.1f}%")

# Fisher's exact tests
or1, p1 = fisher_test_arm_enrichment(coord_arm, coord_core, genome_arm, genome_core)
or2, p2 = fisher_test_arm_enrichment(coord_arm, coord_core, all_reg_arm, all_reg_core)
or3, p3 = fisher_test_arm_enrichment(all_reg_arm, all_reg_core, genome_arm, genome_core)

print(f"\nFisher's exact (one-sided, arm enrichment):")
print(f"  62 coordinated vs genome-wide: OR={or1:.3f}, p={p1:.4e}")
print(f"  62 coordinated vs 1,055 regulators: OR={or2:.3f}, p={p2:.4e}")
print(f"  1,055 regulators vs genome-wide: OR={or3:.3f}, p={or3:.4e}")

# Geographic distribution by coordination type
print("\n--- Geographic distribution by coordination type ---")
coord_types_t2 = coord_df.groupby('coordination_T2')['region'].value_counts().unstack(fill_value=0)
coord_types_t3 = coord_df.groupby('coordination_T3')['region'].value_counts().unstack(fill_value=0)
print("\nT2 coordination types by region:")
print(coord_types_t2)
print("\nT3 coordination types by region:")
print(coord_types_t3)

# Classify combined coordination pattern
def classify_combined(row):
    patterns = []
    for tp in ['coordination_T2', 'coordination_T3']:
        c = row[tp]
        if 'concordant_derepression' in c:
            patterns.append('derepression')
        elif 'concordant_repression' in c:
            patterns.append('repression')
        elif 'discordant' in c:
            patterns.append('discordant')
        else:
            patterns.append('other')
    return '_'.join(patterns)

coord_df['combined_pattern'] = coord_df.apply(classify_combined, axis=1)

# Save geographic distribution table
geo_stats = pd.DataFrame({
    'group': ['62 Coordinated', '1,055 All regulators', 'Genome-wide'],
    'n_arm': [coord_arm, all_reg_arm, genome_arm],
    'n_core': [coord_core, all_reg_core, genome_core],
    'n_total': [len(coord_df), len(all_reg_df), len(gene_annot)],
    'pct_arm': [coord_arm/len(coord_df)*100, all_reg_arm/len(all_reg_df)*100, genome_arm/len(gene_annot)*100],
})
geo_stats['comparison'] = ['vs genome: OR={:.3f}, p={:.2e}'.format(or1, p1),
                            'vs genome: OR={:.3f}, p={:.2e}'.format(or3, p3),
                            'reference']
geo_stats.to_csv(f'{TBL_DIR}/geographic_distribution_stats.tsv', sep='\t', index=False)
print(f"\nSaved: {TBL_DIR}/geographic_distribution_stats.tsv")

# ============================================================
# FIGURE 1: Geographic Distribution (Chromosome Ideogram)
# ============================================================
print("\n--- Generating Figure 1: Geographic Distribution ---")

fig, axes = plt.subplots(3, 1, figsize=(14, 8), gridspec_kw={'height_ratios': [2, 1, 1.5]})

# Panel A: Chromosome ideogram with gene positions
ax = axes[0]
ax.set_xlim(0, CHROM_LEN)
ax.set_ylim(-1.5, 3.5)

# Draw chromosome bar
chrom_y = 1.5
chrom_h = 0.4
# Core region
ax.add_patch(FancyBboxPatch((CORE_START, chrom_y - chrom_h/2), CORE_END - CORE_START, chrom_h,
                              boxstyle="round,pad=0", facecolor='#E0E0E0', edgecolor='black', linewidth=1.5))
# Left arm
ax.add_patch(FancyBboxPatch((0, chrom_y - chrom_h/2), CORE_START, chrom_h,
                              boxstyle="round,pad=0", facecolor='#FFCCCC', edgecolor='black', linewidth=1.5))
# Right arm
ax.add_patch(FancyBboxPatch((CORE_END, chrom_y - chrom_h/2), CHROM_LEN - CORE_END, chrom_h,
                              boxstyle="round,pad=0", facecolor='#FFCCCC', edgecolor='black', linewidth=1.5))

# Labels
ax.text(CORE_START/2, chrom_y, 'Left arm', ha='center', va='center', fontsize=9, fontweight='bold', color='#990000')
ax.text((CORE_START+CORE_END)/2, chrom_y, 'Core region (1.5-6.5 Mb)', ha='center', va='center', fontsize=9, fontweight='bold', color='#333333')
ax.text((CORE_END+CHROM_LEN)/2, chrom_y, 'Right arm', ha='center', va='center', fontsize=9, fontweight='bold', color='#990000')

# Plot gene positions
# Concordant derepression
concordant_derep = coord_df[coord_df['coordination_T2'].str.contains('concordant_derepression') |
                             coord_df['coordination_T3'].str.contains('concordant_derepression')]
# Concordant repression
concordant_repr = coord_df[coord_df['coordination_T2'].str.contains('concordant_repression') |
                            coord_df['coordination_T3'].str.contains('concordant_repression')]
# Discordant
discordant = coord_df[~coord_df.index.isin(concordant_derep.index) & ~coord_df.index.isin(concordant_repr.index)]

# Mark positions on chromosome
for _, row in concordant_derep.iterrows():
    ax.plot([row['midpoint'], row['midpoint']], [chrom_y + chrom_h/2, chrom_y + 0.8], color='#2196F3', linewidth=1.2, alpha=0.8)
    ax.plot(row['midpoint'], chrom_y + 0.8, 'v', color='#2196F3', markersize=6, alpha=0.8)

for _, row in concordant_repr.iterrows():
    ax.plot([row['midpoint'], row['midpoint']], [chrom_y - chrom_h/2, chrom_y - 0.8], color='#F44336', linewidth=1.2, alpha=0.8)
    ax.plot(row['midpoint'], chrom_y - 0.8, '^', color='#F44336', markersize=6, alpha=0.8)

for _, row in discordant.iterrows():
    ax.plot([row['midpoint'], row['midpoint']], [chrom_y + chrom_h/2, chrom_y + 0.5], color='#9E9E9E', linewidth=0.8, alpha=0.6)
    ax.plot(row['midpoint'], chrom_y + 0.5, 'o', color='#9E9E9E', markersize=4, alpha=0.6)

# Top candidate labels
top_candidates = ['SC_RS10435', 'SC_RS31385', 'SC_RS35525']
for lt in top_candidates:
    row = coord_df[coord_df['locus_tag'] == lt].iloc[0]
    ax.annotate(lt, xy=(row['midpoint'], chrom_y + 0.8), xytext=(row['midpoint'], chrom_y + 1.5),
                fontsize=7, fontweight='bold', color='#1565C0',
                arrowprops=dict(arrowstyle='->', color='#1565C0', lw=0.8),
                ha='center', va='bottom')

# Mb scale
for mb in range(0, 9):
    ax.axvline(mb * 1e6, color='gray', linewidth=0.3, linestyle=':', alpha=0.5, ymin=0.1, ymax=0.9)
    ax.text(mb * 1e6, -1.3, f'{mb} Mb', ha='center', fontsize=8, color='gray')

ax.legend([plt.Line2D([0], [0], marker='v', color='#2196F3', linestyle='none', markersize=8),
           plt.Line2D([0], [0], marker='^', color='#F44336', linestyle='none', markersize=8),
           plt.Line2D([0], [0], marker='o', color='#9E9E9E', linestyle='none', markersize=6)],
          ['Concordant derepression', 'Concordant repression', 'Discordant/Other'],
          loc='upper right', fontsize=8, framealpha=0.9)

ax.set_title('A. Chromosomal distribution of 62 coordinated regulatory genes', fontsize=12, fontweight='bold', loc='left')
ax.set_yticks([])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)

# Panel B: Density comparison
ax2 = axes[1]
bins = np.linspace(0, CHROM_LEN, 50)
bin_centers = (bins[:-1] + bins[1:]) / 2

# All regulatory genes density
hist_all, _ = np.histogram(all_reg_df['midpoint'], bins=bins)
hist_coord, _ = np.histogram(coord_df['midpoint'], bins=bins)
hist_genome, _ = np.histogram(gene_annot['midpoint'], bins=bins)

# Normalize to density
hist_all_norm = hist_all / hist_all.sum() if hist_all.sum() > 0 else hist_all
hist_coord_norm = hist_coord / hist_coord.sum() if hist_coord.sum() > 0 else hist_coord
hist_genome_norm = hist_genome / hist_genome.sum() if hist_genome.sum() > 0 else hist_genome

ax2.fill_between(bin_centers, hist_genome_norm, alpha=0.2, color='gray', label=f'Genome-wide (n={len(gene_annot):,})')
ax2.plot(bin_centers, hist_all_norm, color='#FF9800', linewidth=1.5, label=f'All 1,055 regulators')
ax2.plot(bin_centers, hist_coord_norm, color='#E91E63', linewidth=2.5, label=f'62 coordinated')

ax2.axvspan(0, CORE_START, alpha=0.05, color='red')
ax2.axvspan(CORE_END, CHROM_LEN, alpha=0.05, color='red')
ax2.set_ylabel('Density', fontsize=10)
ax2.set_title('B. Gene density distribution along chromosome', fontsize=12, fontweight='bold', loc='left')
ax2.legend(fontsize=8, loc='upper right')
ax2.set_xlim(0, CHROM_LEN)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

# Panel C: Bar chart comparison
ax3 = axes[2]
categories = ['62 Coordinated\nregulators', '1,055 All\nregulators', 'Genome-wide\n(7,825 genes)']
arm_pcts = [coord_arm/len(coord_df)*100, all_reg_arm/len(all_reg_df)*100, genome_arm/len(gene_annot)*100]
core_pcts = [100-p for p in arm_pcts]

x = np.arange(len(categories))
width = 0.35
bars1 = ax3.bar(x - width/2, arm_pcts, width, label='Chromosomal arms', color='#EF5350', edgecolor='black', linewidth=0.5)
bars2 = ax3.bar(x + width/2, core_pcts, width, label='Core region', color='#42A5F5', edgecolor='black', linewidth=0.5)

# Add percentage labels
for bar, pct in zip(bars1, arm_pcts):
    ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1, f'{pct:.1f}%', ha='center', fontsize=9, fontweight='bold')
for bar, pct in zip(bars2, core_pcts):
    ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1, f'{pct:.1f}%', ha='center', fontsize=9, fontweight='bold')

# Expected line
ax3.axhline(arm_fraction*100, color='red', linestyle='--', linewidth=1, alpha=0.7, label=f'Expected arm % by length ({arm_fraction*100:.1f}%)')

ax3.set_ylabel('Percentage (%)', fontsize=10)
ax3.set_xticks(x)
ax3.set_xticklabels(categories)
ax3.legend(fontsize=8, loc='upper right')

# Add Fisher test p-values
ax3.annotate(f'Fisher p={p1:.2e}', xy=(0, arm_pcts[0]), xytext=(0.5, arm_pcts[0]+8),
             fontsize=8, ha='center', color='#B71C1C',
             arrowprops=dict(arrowstyle='->', color='#B71C1C', lw=0.8))

ax3.set_title('C. Arm vs core distribution comparison', fontsize=12, fontweight='bold', loc='left')
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.set_ylim(0, 100)

plt.tight_layout()
for fmt in ['pdf', 'svg']:
    fig.savefig(f'{FIG_DIR}/geographic_distribution.{fmt}')
plt.close()
print(f"Saved: {FIG_DIR}/geographic_distribution.pdf/svg")


# ============================================================
# STEP 2: COORDINATION PATTERN ANALYSIS
# ============================================================
print("\n" + "=" * 70)
print("STEP 2: COORDINATION PATTERN ANALYSIS")
print("=" * 70)

# Classify genes by coordination types at T2 and T3
coord_type_categories = [
    'concordant_derepression',
    'concordant_repression',
    'discordant_gain_up',
    'discordant_loss_down',
    'methyl_change_no_expr_change',
    'ambiguous'
]

# T2 counts by TF family
print("\n--- Coordination at T2 ---")
t2_summary = coord_df.groupby(['coordination_T2', 'tf_family']).size().unstack(fill_value=0)
print(t2_summary)

print("\n--- Coordination at T3 ---")
t3_summary = coord_df.groupby(['coordination_T3', 'tf_family']).size().unstack(fill_value=0)
print(t3_summary)

# Identify genes with concordant_derepression at BOTH T2 and T3
both_derep = coord_df[(coord_df['coordination_T2'] == 'concordant_derepression') &
                       (coord_df['coordination_T3'] == 'concordant_derepression')]
print(f"\n*** Concordant derepression at BOTH T2 and T3: {len(both_derep)} genes ***")
for _, r in both_derep.iterrows():
    print(f"  {r['locus_tag']} ({r['old_locus_tag']}) - {r['product']} [{r['tf_family']}]")
    print(f"    log2FC T2={r['log2FC_T2vsT1_deseq']:.2f}, T3={r['log2FC_T3vsT1_deseq']:.2f}")

# Both repression
both_repr = coord_df[(coord_df['coordination_T2'] == 'concordant_repression') &
                      (coord_df['coordination_T3'] == 'concordant_repression')]
print(f"\n*** Concordant repression at BOTH T2 and T3: {len(both_repr)} genes ***")
for _, r in both_repr.iterrows():
    print(f"  {r['locus_tag']} ({r['old_locus_tag']}) - {r['product']} [{r['tf_family']}]")
    print(f"    log2FC T2={r['log2FC_T2vsT1_deseq']:.2f}, T3={r['log2FC_T3vsT1_deseq']:.2f}")

# Coordination type summary by family
family_coord_t2 = coord_df.groupby('tf_family')['coordination_T2'].value_counts().unstack(fill_value=0)
family_coord_t3 = coord_df.groupby('tf_family')['coordination_T3'].value_counts().unstack(fill_value=0)

# Save coordination summary table
coord_summary = coord_df[['locus_tag', 'gene_name', 'old_locus_tag', 'product', 'tf_family',
                            'start', 'end', 'region', 'coordination_T2', 'coordination_T3',
                            'log2FC_T2vsT1_deseq', 'log2FC_T3vsT1_deseq',
                            'methyl_change']].copy()
coord_summary['both_T2T3_derepression'] = ((coord_df['coordination_T2'] == 'concordant_derepression') &
                                             (coord_df['coordination_T3'] == 'concordant_derepression'))
coord_summary['both_T2T3_repression'] = ((coord_df['coordination_T2'] == 'concordant_repression') &
                                           (coord_df['coordination_T3'] == 'concordant_repression'))
coord_summary.to_csv(f'{TBL_DIR}/coordination_type_summary.tsv', sep='\t', index=False)

# Family coordination summary
family_summary_rows = []
for fam in coord_df['tf_family'].unique():
    fam_df = coord_df[coord_df['tf_family'] == fam]
    n = len(fam_df)
    n_derep_t2 = (fam_df['coordination_T2'] == 'concordant_derepression').sum()
    n_derep_t3 = (fam_df['coordination_T3'] == 'concordant_derepression').sum()
    n_repr_t2 = (fam_df['coordination_T2'] == 'concordant_repression').sum()
    n_repr_t3 = (fam_df['coordination_T3'] == 'concordant_repression').sum()
    n_disc_t2 = fam_df['coordination_T2'].str.contains('discordant').sum()
    n_disc_t3 = fam_df['coordination_T3'].str.contains('discordant').sum()
    family_summary_rows.append({
        'tf_family': fam,
        'n_genes': n,
        'concordant_derep_T2': n_derep_t2,
        'concordant_derep_T3': n_derep_t3,
        'concordant_repr_T2': n_repr_t2,
        'concordant_repr_T3': n_repr_t3,
        'discordant_T2': n_disc_t2,
        'discordant_T3': n_disc_t3,
    })
family_summary = pd.DataFrame(family_summary_rows).sort_values('n_genes', ascending=False)
family_summary.to_csv(f'{TBL_DIR}/coordination_by_family.tsv', sep='\t', index=False)
print(f"\nSaved: {TBL_DIR}/coordination_by_family.tsv")

# ============================================================
# FIGURE 2: Coordination type classification
# ============================================================
print("\n--- Generating Figure 2: Coordination Types ---")

fig, axes = plt.subplots(1, 3, figsize=(16, 6))

# Panel A: T2 coordination types
ax = axes[0]
t2_counts = coord_df['coordination_T2'].value_counts()
colors_t2 = {
    'concordant_derepression': '#2196F3',
    'concordant_repression': '#F44336',
    'discordant_gain_up': '#FF9800',
    'discordant_loss_down': '#9C27B0',
    'methyl_change_no_expr_change': '#607D8B',
    'ambiguous': '#BDBDBD',
}
bars = ax.barh(range(len(t2_counts)), t2_counts.values,
               color=[colors_t2.get(k, '#999') for k in t2_counts.index])
ax.set_yticks(range(len(t2_counts)))
ax.set_yticklabels([s.replace('_', ' ') for s in t2_counts.index], fontsize=8)
for i, (v, k) in enumerate(zip(t2_counts.values, t2_counts.index)):
    ax.text(v + 0.3, i, str(v), va='center', fontsize=9, fontweight='bold')
ax.set_xlabel('Number of genes')
ax.set_title('A. T2 coordination types', fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel B: T3 coordination types
ax = axes[1]
t3_counts = coord_df['coordination_T3'].value_counts()
bars = ax.barh(range(len(t3_counts)), t3_counts.values,
               color=[colors_t2.get(k, '#999') for k in t3_counts.index])
ax.set_yticks(range(len(t3_counts)))
ax.set_yticklabels([s.replace('_', ' ') for s in t3_counts.index], fontsize=8)
for i, (v, k) in enumerate(zip(t3_counts.values, t3_counts.index)):
    ax.text(v + 0.3, i, str(v), va='center', fontsize=9, fontweight='bold')
ax.set_xlabel('Number of genes')
ax.set_title('B. T3 coordination types', fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel C: TF family distribution
ax = axes[2]
family_counts = coord_df['tf_family'].value_counts()
colors_family = sns.color_palette('Set2', n_colors=len(family_counts))
wedges, texts, autotexts = ax.pie(family_counts.values, labels=None, autopct='%1.0f%%',
                                    colors=colors_family, startangle=90, pctdistance=0.85)
for t in autotexts:
    t.set_fontsize(8)
ax.legend(family_counts.index, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=8)
ax.set_title('C. TF family distribution', fontweight='bold')

plt.tight_layout()
for fmt in ['pdf', 'svg']:
    fig.savefig(f'{FIG_DIR}/coordination_types.{fmt}')
plt.close()
print(f"Saved: {FIG_DIR}/coordination_types.pdf/svg")


# ============================================================
# STEP 3: FUNCTIONAL ANNOTATION (COG) ENRICHMENT
# ============================================================
print("\n" + "=" * 70)
print("STEP 3: FUNCTIONAL ANNOTATION (COG) ENRICHMENT")
print("=" * 70)

# Merge COG with gene sets
coord_cog = coord_df.merge(cog_df[['gene_id', 'COG_category']], left_on='locus_tag', right_on='gene_id', how='left')
all_reg_cog = all_reg_df.merge(cog_df[['gene_id', 'COG_category']], left_on='locus_tag', right_on='gene_id', how='left')

print(f"\nCOG annotation coverage:")
print(f"  62 coordinated: {coord_cog['COG_category'].notna().sum()}/{len(coord_cog)} ({coord_cog['COG_category'].notna().sum()/len(coord_cog)*100:.1f}%)")
print(f"  1,055 regulators: {all_reg_cog['COG_category'].notna().sum()}/{len(all_reg_cog)} ({all_reg_cog['COG_category'].notna().sum()/len(all_reg_cog)*100:.1f}%)")
print(f"  Genome-wide: {cog_df['COG_category'].notna().sum()}/{len(cog_df)} ({cog_df['COG_category'].notna().sum()/len(cog_df)*100:.1f}%)")

# COG category descriptions
cog_descriptions = {
    'C': 'C - Energy production',
    'D': 'D - Cell cycle',
    'E': 'E - Amino acid metabolism',
    'F': 'F - Nucleotide metabolism',
    'G': 'G - Carbohydrate metabolism',
    'H': 'H - Coenzyme metabolism',
    'I': 'I - Lipid metabolism',
    'J': 'J - Translation',
    'K': 'K - Transcription',
    'L': 'L - Replication/Repair',
    'M': 'M - Cell wall/membrane',
    'N': 'N - Cell motility',
    'O': 'O - Post-translational mod.',
    'P': 'P - Inorganic ion transport',
    'Q': 'Q - Secondary metabolism',
    'R': 'R - General function',
    'S': 'S - Unknown',
    'T': 'T - Signal transduction',
    'U': 'U - Intracellular trafficking',
    'V': 'V - Defense mechanisms',
    'W': 'W - Extracellular structures',
    'X': 'X - Mobilome',
}

# Extract single-letter COG
def get_cog_letter(cog_str):
    if pd.isna(cog_str) or cog_str == '':
        return 'na'
    return cog_str[0]

coord_cog['cog_letter'] = coord_cog['COG_category'].apply(get_cog_letter)
all_reg_cog['cog_letter'] = all_reg_cog['COG_category'].apply(get_cog_letter)
cog_df['cog_letter'] = cog_df['COG_category'].apply(get_cog_letter)

# Count distribution
coord_cog_counts = coord_cog['cog_letter'].value_counts()
all_reg_cog_counts = all_reg_cog['cog_letter'].value_counts()
genome_cog_counts = cog_df['cog_letter'].value_counts()

# Fisher's exact enrichment for each COG category
enrichment_rows = []
all_cog_letters = sorted(set(coord_cog_counts.index) | set(all_reg_cog_counts.index) | set(genome_cog_counts.index))
for letter in all_cog_letters:
    if letter == 'na':
        continue
    # 62 coordinated vs genome
    a = coord_cog_counts.get(letter, 0)
    b = len(coord_cog) - a
    c = genome_cog_counts.get(letter, 0)
    d = len(cog_df) - c
    if a + c > 0:
        or_val, p_val = stats.fisher_exact([[a, b], [c, d]])
        enrichment_rows.append({
            'COG_category': cog_descriptions.get(letter, letter),
            'coordinated_62_count': a,
            'coordinated_62_pct': a/len(coord_cog)*100,
            'all_reg_1055_count': all_reg_cog_counts.get(letter, 0),
            'all_reg_1055_pct': all_reg_cog_counts.get(letter, 0)/len(all_reg_cog)*100,
            'genome_count': c,
            'genome_pct': c/len(cog_df)*100,
            'odds_ratio_vs_genome': or_val,
            'pvalue_vs_genome': p_val,
        })

enrichment_df = pd.DataFrame(enrichment_rows).sort_values('pvalue_vs_genome')
enrichment_df.to_csv(f'{TBL_DIR}/COG_enrichment_results.tsv', sep='\t', index=False)
print(f"\nSaved: {TBL_DIR}/COG_enrichment_results.tsv")

# Print top enriched
print("\nTop COG categories (62 coordinated vs genome):")
for _, r in enrichment_df.head(10).iterrows():
    sig = '*' if r['pvalue_vs_genome'] < 0.05 else ''
    print(f"  {r['COG_category']}: {r['coordinated_62_count']}/{len(coord_cog)} ({r['coordinated_62_pct']:.1f}%) "
          f"vs genome {r['genome_count']}/{len(cog_df)} ({r['genome_pct']:.1f}%) "
          f"OR={r['odds_ratio_vs_genome']:.2f} p={r['pvalue_vs_genome']:.3e}{sig}")


# ============================================================
# FIGURE 3: COG enrichment comparison
# ============================================================
print("\n--- Generating Figure 3: COG Enrichment ---")

fig, axes = plt.subplots(2, 1, figsize=(14, 10))

# Panel A: COG distribution comparison
ax = axes[0]
relevant_cogs = sorted([l for l in all_cog_letters if l != 'na' and l in cog_descriptions])

coord_pcts = [coord_cog_counts.get(l, 0)/len(coord_cog)*100 for l in relevant_cogs]
allreg_pcts = [all_reg_cog_counts.get(l, 0)/len(all_reg_cog)*100 for l in relevant_cogs]
genome_pcts = [genome_cog_counts.get(l, 0)/len(cog_df)*100 for l in relevant_cogs]

x = np.arange(len(relevant_cogs))
width = 0.25
ax.bar(x - width, coord_pcts, width, label='62 Coordinated', color='#E91E63', edgecolor='black', linewidth=0.5)
ax.bar(x, allreg_pcts, width, label='1,055 Regulators', color='#FF9800', edgecolor='black', linewidth=0.5)
ax.bar(x + width, genome_pcts, width, label='Genome-wide', color='#607D8B', edgecolor='black', linewidth=0.5)

# Mark significant enrichments
for i, l in enumerate(relevant_cogs):
    row = enrichment_df[enrichment_df['COG_category'].str.startswith(l)]
    if len(row) > 0 and row.iloc[0]['pvalue_vs_genome'] < 0.05:
        ax.text(i - width, coord_pcts[i] + 0.5, '*', fontsize=14, fontweight='bold', color='red', ha='center')

ax.set_xticks(x)
ax.set_xticklabels([f'{l}\n{cog_descriptions.get(l, "").split(" - ")[1][:15]}' if l in cog_descriptions else l
                     for l in relevant_cogs], fontsize=7, rotation=45, ha='right')
ax.set_ylabel('Percentage (%)')
ax.legend(fontsize=9)
ax.set_title('A. COG category distribution: 62 coordinated vs 1,055 regulators vs genome', fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel B: Enrichment/depletion heatmap
ax2 = axes[1]
# Calculate log2 fold enrichment vs genome
enrichment_data = {}
for l in relevant_cogs:
    coord_frac = coord_cog_counts.get(l, 0) / len(coord_cog) if len(coord_cog) > 0 else 0
    allreg_frac = all_reg_cog_counts.get(l, 0) / len(all_reg_cog) if len(all_reg_cog) > 0 else 0
    genome_frac = genome_cog_counts.get(l, 0) / len(cog_df) if len(cog_df) > 0 else 0

    if genome_frac > 0:
        enrichment_data[l] = {
            '62 Coordinated': np.log2((coord_frac + 0.001) / (genome_frac + 0.001)),
            '1,055 Regulators': np.log2((allreg_frac + 0.001) / (genome_frac + 0.001)),
        }
    else:
        enrichment_data[l] = {'62 Coordinated': 0, '1,055 Regulators': 0}

enrichment_heatmap = pd.DataFrame(enrichment_data).T
# Filter out categories with very few genes
mask = enrichment_heatmap.index.isin([l for l in relevant_cogs if genome_cog_counts.get(l, 0) >= 20])
enrichment_heatmap = enrichment_heatmap[mask]

sns.heatmap(enrichment_heatmap.T, cmap='RdBu_r', center=0, annot=True, fmt='.1f',
            linewidths=0.5, ax=ax2, cbar_kws={'label': 'log2(fold enrichment vs genome)'},
            xticklabels=[f'{l} - {cog_descriptions.get(l, "").split(" - ")[1][:20]}' if l in cog_descriptions else l
                         for l in enrichment_heatmap.index])
ax2.set_title('B. log2 fold enrichment vs genome-wide background', fontweight='bold')
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45, ha='right', fontsize=8)

plt.tight_layout()
for fmt in ['pdf', 'svg']:
    fig.savefig(f'{FIG_DIR}/COG_enrichment_comparison.{fmt}')
plt.close()
print(f"Saved: {FIG_DIR}/COG_enrichment_comparison.pdf/svg")


# ============================================================
# BGC PROXIMITY ANALYSIS
# ============================================================
print("\n--- BGC Proximity Analysis ---")

# Extract BGC boundaries from gene_master
bgc_genes = gene_master[gene_master['bgc_name'].notna() & (gene_master['bgc_name'] != '')]
if len(bgc_genes) > 0:
    bgc_boundaries = bgc_genes.groupby('bgc_name').agg(
        bgc_start=('start', 'min'),
        bgc_end=('end', 'max'),
        n_genes=('gene_id', 'count')
    ).reset_index()
    print(f"\nBGC clusters found: {len(bgc_boundaries)}")
    print(bgc_boundaries.to_string(index=False))

    # Check 62 genes for BGC proximity (within 20 kb)
    PROXIMITY_KB = 20
    bgc_proximity_rows = []
    for _, gene in coord_df.iterrows():
        gene_mid = gene['midpoint']
        for _, bgc in bgc_boundaries.iterrows():
            dist_to_start = abs(gene_mid - bgc['bgc_start'])
            dist_to_end = abs(gene_mid - bgc['bgc_end'])
            min_dist = min(dist_to_start, dist_to_end)
            # Check if inside or near
            if bgc['bgc_start'] <= gene_mid <= bgc['bgc_end']:
                min_dist = 0
            if min_dist <= PROXIMITY_KB * 1000:
                bgc_proximity_rows.append({
                    'locus_tag': gene['locus_tag'],
                    'old_locus_tag': gene['old_locus_tag'],
                    'product': gene['product'],
                    'tf_family': gene['tf_family'],
                    'gene_midpoint': int(gene_mid),
                    'bgc_name': bgc['bgc_name'],
                    'bgc_start': bgc['bgc_start'],
                    'bgc_end': bgc['bgc_end'],
                    'distance_bp': int(min_dist),
                    'coordination_T2': gene['coordination_T2'],
                    'coordination_T3': gene['coordination_T3'],
                })

    bgc_proximity = pd.DataFrame(bgc_proximity_rows)
    if len(bgc_proximity) > 0:
        print(f"\nGenes within {PROXIMITY_KB} kb of BGC boundaries:")
        for _, r in bgc_proximity.iterrows():
            print(f"  {r['locus_tag']} ({r['old_locus_tag']}) - {r['distance_bp']/1000:.1f} kb from {r['bgc_name']}")
    else:
        print(f"\nNo genes within {PROXIMITY_KB} kb of BGC boundaries")
else:
    print("No BGC annotations found in gene_master")
    bgc_proximity = pd.DataFrame()
    bgc_boundaries = pd.DataFrame()


# ============================================================
# STEP 4: NEIGHBORHOOD ANALYSIS FOR TOP CANDIDATES
# ============================================================
print("\n" + "=" * 70)
print("STEP 4: NEIGHBORHOOD ANALYSIS (TOP CANDIDATES)")
print("=" * 70)

top_candidates = ['SC_RS10435', 'SC_RS31385', 'SC_RS35525']
WINDOW = 10_000  # +/- 10 kb

fig, axes = plt.subplots(len(top_candidates), 1, figsize=(16, 5*len(top_candidates)))

for idx, candidate in enumerate(top_candidates):
    ax = axes[idx]

    # Get candidate info
    cand_info = coord_df[coord_df['locus_tag'] == candidate].iloc[0]
    cand_mid = cand_info['midpoint']
    cand_start = cand_info['start']
    cand_end = cand_info['end']

    # Find neighboring genes
    window_start = max(0, cand_start - WINDOW)
    window_end = min(CHROM_LEN, cand_end + WINDOW)

    neighbors = gene_annot[(gene_annot['start'] >= window_start) & (gene_annot['end'] <= window_end)].copy()

    # Merge with expression data
    neighbors = neighbors.merge(
        gene_master[['gene_id', 'baseMean', 'log2FoldChange_2_vs_1', 'padj_2_vs_1',
                      'log2FoldChange_3_vs_1', 'padj_3_vs_1', 'bgc_name', 'bgc_role']],
        on='gene_id', how='left'
    )

    # Also merge with COG
    neighbors = neighbors.merge(cog_df[['gene_id', 'COG_category']], on='gene_id', how='left')

    print(f"\n{candidate} ({cand_info['old_locus_tag']}) - {cand_info['product']}")
    print(f"  Position: {cand_start:,}-{cand_end:,} ({cand_info['region']})")
    print(f"  Coordination: T2={cand_info['coordination_T2']}, T3={cand_info['coordination_T3']}")
    print(f"  log2FC: T2={cand_info['log2FC_T2vsT1_deseq']:.2f}, T3={cand_info['log2FC_T3vsT1_deseq']:.2f}")
    print(f"  Neighbors in +/-10 kb window: {len(neighbors)}")

    for _, n in neighbors.iterrows():
        deg_mark = ''
        if pd.notna(n.get('padj_2_vs_1')) and n.get('padj_2_vs_1', 1) < 0.01:
            deg_mark += ' DEG_T2'
        if pd.notna(n.get('padj_3_vs_1')) and n.get('padj_3_vs_1', 1) < 0.01:
            deg_mark += ' DEG_T3'
        bgc_mark = f" [{n['bgc_name']}/{n['bgc_role']}]" if pd.notna(n.get('bgc_name')) and n.get('bgc_name') else ''
        gene_label = n['gene_name'] if pd.notna(n.get('gene_name')) and n.get('gene_name') else n['gene_id']
        product_str = str(n['product'])[:50] if pd.notna(n.get('product')) else 'unknown'
        print(f"    {gene_label} ({n['gene_id']}) {n['start']:,}-{n['end']:,} {n['strand']} "
              f"{product_str}{bgc_mark}{deg_mark}")

    # Draw neighborhood diagram
    ax.set_xlim(window_start, window_end)
    ax.set_ylim(-2, 3)

    # Draw genes as arrows
    for _, n in neighbors.iterrows():
        gene_start = n['start']
        gene_end = n['end']
        gene_len = gene_end - gene_start
        is_candidate = (n['gene_id'] == candidate)

        # Color coding
        if is_candidate:
            color = '#E91E63'
            alpha = 1.0
            lw = 2
        elif pd.notna(n.get('bgc_name')) and n.get('bgc_name'):
            color = '#4CAF50'
            alpha = 0.8
            lw = 1
        else:
            # Color by COG
            cog = n.get('COG_category', '')
            if pd.notna(cog) and cog.startswith('K'):
                color = '#FF9800'
            elif pd.notna(cog) and cog.startswith('T'):
                color = '#2196F3'
            elif pd.notna(cog) and cog.startswith('Q'):
                color = '#9C27B0'
            else:
                color = '#B0BEC5'
            alpha = 0.7
            lw = 0.8

        y_pos = 0.5
        h = 0.6

        if n['strand'] == '+':
            arrow = mpatches.FancyArrow(gene_start, y_pos, gene_len, 0,
                                         width=h, head_width=h*1.2, head_length=min(gene_len*0.15, 300),
                                         fc=color, ec='black', linewidth=lw, alpha=alpha)
        else:
            arrow = mpatches.FancyArrow(gene_end, y_pos, -gene_len, 0,
                                         width=h, head_width=h*1.2, head_length=min(gene_len*0.15, 300),
                                         fc=color, ec='black', linewidth=lw, alpha=alpha)
        ax.add_patch(arrow)

        # Gene label
        label = n['gene_name'] if pd.notna(n.get('gene_name')) and n.get('gene_name') else n['gene_id'].replace('SC_RS', '')
        ax.text((gene_start + gene_end)/2, 1.4, label, ha='center', va='bottom',
                fontsize=7, fontweight='bold' if is_candidate else 'normal',
                color='#C2185B' if is_candidate else 'black', rotation=30)

    # Draw scale bar
    scale_bp = 2000
    ax.plot([window_start + 500, window_start + 500 + scale_bp], [-1.3, -1.3], 'k-', linewidth=2)
    ax.text(window_start + 500 + scale_bp/2, -1.6, f'{scale_bp/1000:.0f} kb', ha='center', fontsize=8)

    ax.set_title(f'{chr(65+idx)}. {candidate} ({cand_info["old_locus_tag"]}) - {cand_info["product"]} '
                 f'[{cand_info["tf_family"]}]\n'
                 f'Concordant derepression T2+T3 | log2FC: T2={cand_info["log2FC_T2vsT1_deseq"]:.2f}, T3={cand_info["log2FC_T3vsT1_deseq"]:.2f}',
                 fontsize=10, fontweight='bold', loc='left')
    ax.set_xlabel(f'Chromosome position (bp)')
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

plt.tight_layout()
for fmt in ['pdf', 'svg']:
    fig.savefig(f'{FIG_DIR}/top_candidate_neighborhoods.{fmt}')
plt.close()
print(f"\nSaved: {FIG_DIR}/top_candidate_neighborhoods.pdf/svg")


# ============================================================
# STEP 5: TWO-COMPONENT SYSTEM (TCS) ANALYSIS
# ============================================================
print("\n" + "=" * 70)
print("STEP 5: TWO-COMPONENT SYSTEM (TCS) ANALYSIS")
print("=" * 70)

# Find sensor kinases and response regulators among 62 genes
sensor_kinases = coord_df[coord_df['tf_family'] == 'Sensor kinase'].copy()
response_regulators = coord_df[coord_df['tf_family'] == 'Response regulator'].copy()

print(f"\nSensor kinases: {len(sensor_kinases)}")
for _, sk in sensor_kinases.iterrows():
    print(f"  {sk['locus_tag']} ({sk['old_locus_tag']}) pos={sk['start']:,}-{sk['end']:,} {sk['strand']}")

print(f"\nResponse regulators: {len(response_regulators)}")
for _, rr in response_regulators.iterrows():
    print(f"  {rr['locus_tag']} ({rr['old_locus_tag']}) pos={rr['start']:,}-{rr['end']:,} {rr['strand']}")

# Check for cognate pairs (adjacent on genome, within 2 kb of each other)
TCS_PAIR_DIST = 5000
tcs_pairs = []

# Also look for SK/RR pairs in the full regulatory list (not just the 62)
all_sk = all_reg_df[all_reg_df['tf_family'] == 'Sensor kinase']
all_rr = all_reg_df[all_reg_df['tf_family'] == 'Response regulator']

for _, sk in sensor_kinases.iterrows():
    for _, rr in all_rr.iterrows():
        dist = min(abs(sk['start'] - rr['end']), abs(sk['end'] - rr['start']))
        if dist <= TCS_PAIR_DIST:
            both_coordinated = rr['locus_tag'] in coord_df['locus_tag'].values
            tcs_pairs.append({
                'sensor_kinase': sk['locus_tag'],
                'sk_old_locus': sk['old_locus_tag'],
                'sk_start': sk['start'],
                'sk_end': sk['end'],
                'sk_strand': sk['strand'],
                'sk_coordination_T2': sk['coordination_T2'],
                'sk_coordination_T3': sk['coordination_T3'],
                'response_regulator': rr['locus_tag'],
                'rr_old_locus': rr['old_locus_tag'],
                'rr_start': rr['start'],
                'rr_end': rr['end'],
                'rr_strand': rr['strand'],
                'rr_coordination_T2': rr['coordination_T2'],
                'rr_coordination_T3': rr['coordination_T3'],
                'distance_bp': dist,
                'both_in_62': both_coordinated,
            })

for _, rr in response_regulators.iterrows():
    for _, sk in all_sk.iterrows():
        dist = min(abs(sk['start'] - rr['end']), abs(sk['end'] - rr['start']))
        if dist <= TCS_PAIR_DIST:
            both_coordinated = sk['locus_tag'] in coord_df['locus_tag'].values
            # Avoid duplicates
            existing = [(p['sensor_kinase'], p['response_regulator']) for p in tcs_pairs]
            if (sk['locus_tag'], rr['locus_tag']) not in existing:
                tcs_pairs.append({
                    'sensor_kinase': sk['locus_tag'],
                    'sk_old_locus': sk['old_locus_tag'],
                    'sk_start': sk['start'],
                    'sk_end': sk['end'],
                    'sk_strand': sk['strand'],
                    'sk_coordination_T2': sk['coordination_T2'],
                    'sk_coordination_T3': sk['coordination_T3'],
                    'response_regulator': rr['locus_tag'],
                    'rr_old_locus': rr['old_locus_tag'],
                    'rr_start': rr['start'],
                    'rr_end': rr['end'],
                    'rr_strand': rr['strand'],
                    'rr_coordination_T2': rr['coordination_T2'],
                    'rr_coordination_T3': rr['coordination_T3'],
                    'distance_bp': dist,
                    'both_in_62': both_coordinated,
                })

tcs_pairs_df = pd.DataFrame(tcs_pairs)
if len(tcs_pairs_df) > 0:
    tcs_pairs_df = tcs_pairs_df.sort_values('distance_bp')
    tcs_pairs_df.to_csv(f'{TBL_DIR}/TCS_pair_candidates.tsv', sep='\t', index=False)
    print(f"\nFound {len(tcs_pairs_df)} potential TCS pairs (within {TCS_PAIR_DIST/1000:.0f} kb):")
    for _, p in tcs_pairs_df.iterrows():
        both_mark = " *** BOTH COORDINATED ***" if p['both_in_62'] else ""
        print(f"  SK: {p['sensor_kinase']} ({p['sk_old_locus']}) -- "
              f"RR: {p['response_regulator']} ({p['rr_old_locus']}) "
              f"dist={p['distance_bp']:,} bp{both_mark}")
        print(f"    SK coord: T2={p['sk_coordination_T2']}, T3={p['sk_coordination_T3']}")
        print(f"    RR coord: T2={p['rr_coordination_T2']}, T3={p['rr_coordination_T3']}")
else:
    print("\nNo TCS pairs found within distance threshold")
    tcs_pairs_df = pd.DataFrame()

# Save TCS candidates (even if empty)
tcs_pairs_df.to_csv(f'{TBL_DIR}/TCS_pair_candidates.tsv', sep='\t', index=False)
print(f"\nSaved: {TBL_DIR}/TCS_pair_candidates.tsv")


# ============================================================
# ADDITIONAL: Arm enrichment by coordination TYPE
# ============================================================
print("\n" + "=" * 70)
print("ADDITIONAL: Arm Enrichment by Coordination Type")
print("=" * 70)

# Genes with ANY concordant derepression (T2 or T3)
derep_genes = coord_df[(coord_df['coordination_T2'] == 'concordant_derepression') |
                        (coord_df['coordination_T3'] == 'concordant_derepression')]
repr_genes = coord_df[(coord_df['coordination_T2'] == 'concordant_repression') |
                       (coord_df['coordination_T3'] == 'concordant_repression')]
disc_genes = coord_df[~coord_df.index.isin(derep_genes.index) & ~coord_df.index.isin(repr_genes.index)]

for label, subset in [('Concordant derepression', derep_genes),
                        ('Concordant repression', repr_genes),
                        ('Discordant/Other', disc_genes)]:
    n_arm = (subset['region'] == 'arm').sum()
    n_core = (subset['region'] == 'core').sum()
    n_total = len(subset)
    if n_total > 0:
        pct_arm = n_arm / n_total * 100
        or_val, p_val = fisher_test_arm_enrichment(n_arm, n_core, genome_arm, genome_core)
        print(f"  {label}: {n_arm}/{n_total} arm ({pct_arm:.1f}%), OR={or_val:.2f}, p={p_val:.3e}")


# ============================================================
# SUMMARY FIGURE: Comprehensive Overview
# ============================================================
print("\n--- Generating Summary Figure ---")

fig = plt.figure(figsize=(18, 14))
gs = fig.add_gridspec(3, 3, hspace=0.4, wspace=0.35)

# Panel 1: Chromosome map (wide, top)
ax1 = fig.add_subplot(gs[0, :])
ax1.set_xlim(0, CHROM_LEN)
ax1.set_ylim(-0.5, 2.5)

# Draw chromosome
chrom_y = 1.0
chrom_h = 0.3
ax1.add_patch(FancyBboxPatch((CORE_START, chrom_y - chrom_h/2), CORE_END - CORE_START, chrom_h,
                              boxstyle="round,pad=0", facecolor='#F5F5F5', edgecolor='black', linewidth=1))
ax1.add_patch(FancyBboxPatch((0, chrom_y - chrom_h/2), CORE_START, chrom_h,
                              boxstyle="round,pad=0", facecolor='#FFEBEE', edgecolor='black', linewidth=1))
ax1.add_patch(FancyBboxPatch((CORE_END, chrom_y - chrom_h/2), CHROM_LEN - CORE_END, chrom_h,
                              boxstyle="round,pad=0", facecolor='#FFEBEE', edgecolor='black', linewidth=1))

# Gene marks colored by coordination
for _, row in coord_df.iterrows():
    if row['coordination_T2'] == 'concordant_derepression' or row['coordination_T3'] == 'concordant_derepression':
        color, marker = '#1565C0', 'v'
        y_offset = 0.6
    elif row['coordination_T2'] == 'concordant_repression' or row['coordination_T3'] == 'concordant_repression':
        color, marker = '#C62828', '^'
        y_offset = -0.6
    else:
        color, marker = '#757575', 'o'
        y_offset = 0.4
    ax1.plot(row['midpoint'], chrom_y + y_offset, marker, color=color, markersize=5, alpha=0.8)

for mb in range(0, 9):
    ax1.text(mb * 1e6, -0.3, f'{mb}', ha='center', fontsize=8, color='gray')
ax1.text(CHROM_LEN/2, -0.5, 'Mb', ha='center', fontsize=9, color='gray')

ax1.set_title('62 Coordinated Regulatory Genes on S. coelicolor Chromosome', fontsize=13, fontweight='bold')
ax1.set_yticks([])
for spine in ax1.spines.values():
    spine.set_visible(False)

# Panel 2: Arm vs Core bar
ax2 = fig.add_subplot(gs[1, 0])
groups = ['62\nCoordinated', '1,055\nRegulators', 'Genome']
arm_vals = [coord_arm/len(coord_df)*100, all_reg_arm/len(all_reg_df)*100, genome_arm/len(gene_annot)*100]
core_vals = [100-v for v in arm_vals]
ax2.bar(groups, arm_vals, color='#EF5350', label='Arms', edgecolor='black', linewidth=0.5)
ax2.bar(groups, core_vals, bottom=arm_vals, color='#42A5F5', label='Core', edgecolor='black', linewidth=0.5)
ax2.axhline(arm_fraction*100, color='red', linestyle='--', linewidth=1, alpha=0.7)
ax2.set_ylabel('Percentage (%)')
ax2.legend(fontsize=8)
ax2.set_title('Arm vs Core', fontweight='bold', fontsize=10)
ax2.text(0, arm_vals[0]+2, f'p={p1:.1e}', ha='center', fontsize=7, color='red')

# Panel 3: Coordination types
ax3 = fig.add_subplot(gs[1, 1])
# Simplified: count unique patterns
pattern_summary = {
    'Derepression\n(any T)': len(derep_genes),
    'Repression\n(any T)': len(repr_genes),
    'Discordant\n/Other': len(disc_genes),
}
colors_p = ['#1565C0', '#C62828', '#757575']
bars = ax3.bar(pattern_summary.keys(), pattern_summary.values(), color=colors_p, edgecolor='black', linewidth=0.5)
for bar, val in zip(bars, pattern_summary.values()):
    ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5, str(val), ha='center', fontweight='bold')
ax3.set_ylabel('Number of genes')
ax3.set_title('Coordination Patterns', fontweight='bold', fontsize=10)
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)

# Panel 4: TF family
ax4 = fig.add_subplot(gs[1, 2])
top_fams = family_counts.head(8)
colors_f = sns.color_palette('Set2', n_colors=len(top_fams))
bars = ax4.barh(range(len(top_fams)), top_fams.values, color=colors_f, edgecolor='black', linewidth=0.5)
ax4.set_yticks(range(len(top_fams)))
ax4.set_yticklabels(top_fams.index, fontsize=8)
for i, v in enumerate(top_fams.values):
    ax4.text(v + 0.2, i, str(v), va='center', fontsize=9)
ax4.set_xlabel('Count')
ax4.set_title('TF Family Distribution', fontweight='bold', fontsize=10)
ax4.spines['top'].set_visible(False)
ax4.spines['right'].set_visible(False)
ax4.invert_yaxis()

# Panel 5: Top COG enrichments
ax5 = fig.add_subplot(gs[2, 0:2])
# Show top 8 most enriched/depleted COGs
enrichment_show = enrichment_df[enrichment_df['coordinated_62_pct'] > 0].head(10)
x_pos = np.arange(len(enrichment_show))
ax5.bar(x_pos - 0.15, enrichment_show['coordinated_62_pct'], 0.3, label='62 Coordinated', color='#E91E63')
ax5.bar(x_pos + 0.15, enrichment_show['genome_pct'], 0.3, label='Genome', color='#607D8B')
ax5.set_xticks(x_pos)
ax5.set_xticklabels([c.split(' - ')[1][:18] if ' - ' in c else c for c in enrichment_show['COG_category']],
                     rotation=45, ha='right', fontsize=8)
ax5.set_ylabel('Percentage (%)')
ax5.legend(fontsize=8)
ax5.set_title('COG Category Comparison', fontweight='bold', fontsize=10)
ax5.spines['top'].set_visible(False)
ax5.spines['right'].set_visible(False)

# Panel 6: TCS pairs summary
ax6 = fig.add_subplot(gs[2, 2])
ax6.axis('off')
tcs_text = "TCS Pair Analysis\n" + "=" * 25 + "\n\n"
tcs_text += f"Sensor kinases: {len(sensor_kinases)}\n"
tcs_text += f"Response regulators: {len(response_regulators)}\n\n"
if len(tcs_pairs_df) > 0:
    both_coordinated = tcs_pairs_df[tcs_pairs_df['both_in_62']]
    tcs_text += f"Cognate pairs found: {len(tcs_pairs_df)}\n"
    tcs_text += f"Both coordinated: {len(both_coordinated)}\n\n"
    for _, p in tcs_pairs_df.head(3).iterrows():
        tcs_text += f"SK: {p['sensor_kinase']}\nRR: {p['response_regulator']}\n"
        tcs_text += f"Dist: {p['distance_bp']:,} bp\n"
        tcs_text += f"Both in 62: {p['both_in_62']}\n\n"
else:
    tcs_text += "No cognate pairs found\nwithin 5 kb"

ax6.text(0.05, 0.95, tcs_text, transform=ax6.transAxes, fontsize=8, verticalalignment='top',
         fontfamily='monospace', bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

fig.suptitle('H8: Geographic and Functional Characterization of 62 Coordinated Regulatory Genes',
             fontsize=14, fontweight='bold', y=0.98)

plt.tight_layout(rect=[0, 0, 1, 0.96])
for fmt in ['pdf', 'svg']:
    fig.savefig(f'{FIG_DIR}/H8_summary_overview.{fmt}')
plt.close()
print(f"Saved: {FIG_DIR}/H8_summary_overview.pdf/svg")


# ============================================================
# FINAL SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("H8 ANALYSIS SUMMARY")
print("=" * 70)

print(f"""
GEOGRAPHIC DISTRIBUTION:
  62 coordinated: {coord_arm}/{len(coord_df)} arm ({coord_arm/len(coord_df)*100:.1f}%)
  1,055 regulators: {all_reg_arm}/{len(all_reg_df)} arm ({all_reg_arm/len(all_reg_df)*100:.1f}%)
  Genome-wide: {genome_arm}/{len(gene_annot)} arm ({genome_arm/len(gene_annot)*100:.1f}%)
  Expected by length: {arm_fraction*100:.1f}%

  Fisher's exact (62 vs genome): OR={or1:.3f}, p={p1:.2e}
  Fisher's exact (62 vs 1,055): OR={or2:.3f}, p={p2:.2e}

COORDINATION PATTERNS:
  Concordant derepression (any T): {len(derep_genes)} genes
  Concordant repression (any T): {len(repr_genes)} genes
  Discordant/Other: {len(disc_genes)} genes

  Both T2+T3 derepression: {len(both_derep)} genes
  Both T2+T3 repression: {len(both_repr)} genes

TOP TF FAMILIES:
{family_counts.head(5).to_string()}

TCS PAIRS: {len(tcs_pairs_df)} candidates found

H8 HYPOTHESIS VERDICT:
  Geographic enrichment at arms: {'SUPPORTED' if p1 < 0.05 else 'NOT SUPPORTED'} (p={p1:.2e})
  Functional enrichment: See COG analysis above
""")

print("\nAll outputs saved to:")
print(f"  Figures: {FIG_DIR}/")
print(f"  Tables: {TBL_DIR}/")
print("\nDone!")
