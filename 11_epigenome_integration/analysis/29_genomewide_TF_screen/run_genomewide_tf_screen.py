#!/usr/bin/env python3
"""
H6: Genome-wide regulatory gene methylation-expression screening
Screens all ~865 regulatory genes from GFF against methylation and expression data.
"""

import pandas as pd
import numpy as np
import re
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns

# ========== Configuration ==========
BASE_DIR = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/29_genomewide_TF_screen"
FIG_DIR = os.path.join(BASE_DIR, "figures")
TAB_DIR = os.path.join(BASE_DIR, "tables")

GFF_PATH = "/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/genomic.gff"
METHYL_EXPR_PATH = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/01_integration/integrated_methyl_expression_weighted.csv"
HC_SITES_PATH = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/01_integration/high_confidence_sites_weighted.csv"
DEG_T2_PATH = "/Users/okaban/bioinfo/rna-seq/04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_2_vs_1.tsv"
DEG_T3_PATH = "/Users/okaban/bioinfo/rna-seq/04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_3_vs_1.tsv"
LIT_TF_PATH = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/12_grn_tf_methylation/literature_tf_master.csv"

os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TAB_DIR, exist_ok=True)

# ========== 1. Parse GFF for all regulatory genes ==========
print("=" * 60)
print("Step 1: Parsing GFF for regulatory genes...")
print("=" * 60)

# Keywords to identify regulatory genes from product field
REGULATORY_KEYWORDS = [
    r'regulator', r'transcription', r'DNA-binding', r'DNA binding',
    r'response regulator', r'sensor kinase', r'sensor histidine kinase',
    r'sigma factor', r'sigma-70', r'ECF sigma',
    r'repressor', r'activator',
    r'MarR', r'TetR', r'GntR', r'LysR', r'AraC', r'LacI', r'IclR',
    r'MerR', r'WhiB', r'Lrp', r'XRE', r'SARP', r'AfsR', r'ArsR',
    r'helix-turn-helix', r'HTH',
    r'two-component', r'anti-sigma',
    r'Crp/Fnr', r'CRP/FNR',
    r'DeoR', r'PadR', r'RpiR', r'ROK',
    r'winged helix',
]

# Compile regex pattern (case-insensitive)
keyword_pattern = re.compile('|'.join(REGULATORY_KEYWORDS), re.IGNORECASE)

# TF family classification patterns (order matters for specificity)
FAMILY_PATTERNS = [
    (r'TetR', re.compile(r'TetR', re.IGNORECASE)),
    (r'MarR', re.compile(r'MarR', re.IGNORECASE)),
    (r'GntR', re.compile(r'GntR', re.IGNORECASE)),
    (r'LysR', re.compile(r'LysR', re.IGNORECASE)),
    (r'AraC', re.compile(r'AraC', re.IGNORECASE)),
    (r'LacI', re.compile(r'LacI', re.IGNORECASE)),
    (r'IclR', re.compile(r'IclR', re.IGNORECASE)),
    (r'MerR', re.compile(r'MerR', re.IGNORECASE)),
    (r'WhiB', re.compile(r'WhiB|Wbl', re.IGNORECASE)),
    (r'Lrp', re.compile(r'Lrp|feast.famine', re.IGNORECASE)),
    (r'XRE', re.compile(r'XRE', re.IGNORECASE)),
    (r'SARP', re.compile(r'SARP', re.IGNORECASE)),
    (r'ArsR', re.compile(r'ArsR', re.IGNORECASE)),
    (r'DeoR', re.compile(r'DeoR', re.IGNORECASE)),
    (r'PadR', re.compile(r'PadR', re.IGNORECASE)),
    (r'ROK', re.compile(r'ROK', re.IGNORECASE)),
    (r'Sigma factor', re.compile(r'sigma factor|sigma-70|ECF sigma|σ', re.IGNORECASE)),
    (r'Sensor kinase', re.compile(r'sensor.*(kinase|histidine)|histidine kinase', re.IGNORECASE)),
    (r'Response regulator', re.compile(r'response regulator', re.IGNORECASE)),
    (r'Two-component', re.compile(r'two-component', re.IGNORECASE)),
    (r'Anti-sigma', re.compile(r'anti-sigma', re.IGNORECASE)),
    (r'CRP/FNR', re.compile(r'Crp|FNR|CRP', re.IGNORECASE)),
    (r'HTH (other)', re.compile(r'helix-turn-helix|HTH|winged helix', re.IGNORECASE)),
    (r'Repressor', re.compile(r'repressor', re.IGNORECASE)),
    (r'Activator', re.compile(r'activator', re.IGNORECASE)),
    (r'Other regulatory', re.compile(r'regulator|transcription|DNA.binding', re.IGNORECASE)),
]

def classify_tf_family(product):
    """Classify a product description into a TF family."""
    for family_name, pattern in FAMILY_PATTERNS:
        if pattern.search(product):
            return family_name
    return 'Other regulatory'

def parse_gff_attributes(attr_str):
    """Parse GFF9 attribute string into a dict."""
    attrs = {}
    for item in attr_str.split(';'):
        if '=' in item:
            key, val = item.split('=', 1)
            attrs[key] = val.replace('%20', ' ').replace('%2C', ',').replace('%3B', ';').replace('%3D', '=').replace('%25', '%').replace('%09', '\t').replace('%0A', '\n').replace('%26', '&').replace('%23', '#')
    return attrs

# Parse GFF
regulatory_genes = []
gene_names = {}  # locus_tag -> gene name
old_locus_tags = {}  # locus_tag -> old_locus_tag (SCO number)

with open(GFF_PATH, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue

        feature_type = fields[2]
        attrs = parse_gff_attributes(fields[8])

        # Collect gene names and old_locus_tag from gene features
        if feature_type == 'gene':
            lt = attrs.get('locus_tag', '')
            if 'Name' in attrs and attrs['Name'] != lt:
                gene_names[lt] = attrs['Name']
            if 'old_locus_tag' in attrs:
                old_locus_tags[lt] = attrs['old_locus_tag']

        # Only process CDS features
        if feature_type != 'CDS':
            continue

        product = attrs.get('product', '')
        locus_tag = attrs.get('locus_tag', '')

        if not locus_tag:
            continue

        # Check if product matches regulatory keywords
        if keyword_pattern.search(product):
            family = classify_tf_family(product)
            gene_name = gene_names.get(locus_tag, '')
            sco_number = old_locus_tags.get(locus_tag, '')

            regulatory_genes.append({
                'locus_tag': locus_tag,
                'gene_name': gene_name,
                'old_locus_tag': sco_number,
                'product': product,
                'tf_family': family,
                'start': int(fields[3]),
                'end': int(fields[4]),
                'strand': fields[6],
            })

reg_df = pd.DataFrame(regulatory_genes)
# Remove duplicates (some CDS may appear multiple times for pseudo etc.)
reg_df = reg_df.drop_duplicates(subset='locus_tag', keep='first')
reg_df = reg_df.sort_values('start').reset_index(drop=True)

print(f"Total regulatory genes extracted: {len(reg_df)}")
print(f"\nFamily breakdown:")
family_counts = reg_df['tf_family'].value_counts()
for fam, cnt in family_counts.items():
    print(f"  {fam}: {cnt}")

# ========== 2. Load methylation-expression data ==========
print("\n" + "=" * 60)
print("Step 2: Loading methylation-expression data...")
print("=" * 60)

methyl_expr = pd.read_csv(METHYL_EXPR_PATH)
print(f"Methylation-expression table: {len(methyl_expr)} genes")

# Load DESeq2 data
deg_t2 = pd.read_csv(DEG_T2_PATH, sep='\t')
deg_t3 = pd.read_csv(DEG_T3_PATH, sep='\t')
print(f"DESeq2 T2vsT1: {len(deg_t2)} genes")
print(f"DESeq2 T3vsT1: {len(deg_t3)} genes")

# Load HC sites
hc_sites = pd.read_csv(HC_SITES_PATH)
print(f"High-confidence sites: {len(hc_sites)} sites")

# Load literature TF list
lit_tf = pd.read_csv(LIT_TF_PATH)
lit_tf_locus_tags = set(lit_tf['locus_tag'].dropna().values)
print(f"Literature TF list: {len(lit_tf_locus_tags)} TFs")

# ========== 3. Cross-validation ==========
print("\n" + "=" * 60)
print("Step 3: Cross-validating locus_tags...")
print("=" * 60)

reg_locus_tags = set(reg_df['locus_tag'].values)
methyl_locus_tags = set(methyl_expr['gene_id'].values)

matched = reg_locus_tags & methyl_locus_tags
unmatched = reg_locus_tags - methyl_locus_tags
print(f"Regulatory genes in methyl-expr data: {len(matched)}/{len(reg_locus_tags)}")
if unmatched:
    print(f"Unmatched locus_tags ({len(unmatched)}): {sorted(list(unmatched))[:10]}...")

# Check literature TFs in reg_df
lit_in_gff = lit_tf_locus_tags & reg_locus_tags
lit_not_in_gff = lit_tf_locus_tags - reg_locus_tags
print(f"\nLiterature TFs found in GFF regulatory set: {len(lit_in_gff)}/{len(lit_tf_locus_tags)}")
if lit_not_in_gff:
    missing_names = lit_tf[lit_tf['locus_tag'].isin(lit_not_in_gff)][['name','locus_tag','product']].to_string(index=False)
    print(f"Literature TFs NOT in GFF regulatory set (non-canonical regulators):")
    print(missing_names)

# ========== 4. Integrate data ==========
print("\n" + "=" * 60)
print("Step 4: Integrating methylation, expression, and regulatory data...")
print("=" * 60)

# Merge regulatory genes with methylation-expression data
df = reg_df.merge(methyl_expr, left_on='locus_tag', right_on='gene_id', how='left')

# Add DESeq2 data explicitly
deg_t2_sub = deg_t2[['gene_id', 'log2FoldChange', 'padj']].rename(
    columns={'log2FoldChange': 'log2FC_T2vsT1_deseq', 'padj': 'padj_T2vsT1_deseq'})
deg_t3_sub = deg_t3[['gene_id', 'log2FoldChange', 'padj']].rename(
    columns={'log2FoldChange': 'log2FC_T3vsT1_deseq', 'padj': 'padj_T3vsT1_deseq'})

df = df.merge(deg_t2_sub, left_on='locus_tag', right_on='gene_id', how='left', suffixes=('', '_t2'))
df = df.merge(deg_t3_sub, left_on='locus_tag', right_on='gene_id', how='left', suffixes=('', '_t3'))

# Clean up duplicate gene_id columns
for col in ['gene_id', 'gene_id_t2', 'gene_id_t3']:
    if col in df.columns:
        df.drop(col, axis=1, inplace=True)

# Calculate total methylation sites per gene per timepoint
df['total_methyl_T1'] = df['6mA_T1_count'].fillna(0) + df['4mC_T1_count'].fillna(0)
df['total_methyl_T2'] = df['6mA_T2_count'].fillna(0) + df['4mC_T2_count'].fillna(0)
df['total_methyl_T3'] = df['6mA_T3_count'].fillna(0) + df['4mC_T3_count'].fillna(0)
df['has_methylation'] = (df['total_methyl_T1'] + df['total_methyl_T2'] + df['total_methyl_T3']) > 0

# Methylation change status
def methyl_status(row):
    """Determine methylation change: gained, lost, stable, or none."""
    t1 = row['total_methyl_T1']
    t2 = row['total_methyl_T2']
    t3 = row['total_methyl_T3']
    if t1 == 0 and t2 == 0 and t3 == 0:
        return 'none'
    if t1 == t2 == t3:
        return 'constitutive'
    # Check T2 vs T1
    changes = []
    if t2 > t1:
        changes.append('gained_T2')
    elif t2 < t1:
        changes.append('lost_T2')
    if t3 > t1:
        changes.append('gained_T3')
    elif t3 < t1:
        changes.append('lost_T3')
    if not changes:
        return 'constitutive'
    return ','.join(changes)

df['methyl_change'] = df.apply(methyl_status, axis=1)

# DEG status
DEG_LFC_CUTOFF = 1.0
DEG_PADJ_CUTOFF = 0.05

df['DEG_T2vsT1'] = (df['padj_T2vsT1_deseq'].fillna(1) < DEG_PADJ_CUTOFF) & \
                    (df['log2FC_T2vsT1_deseq'].fillna(0).abs() >= DEG_LFC_CUTOFF)
df['DEG_T3vsT1'] = (df['padj_T3vsT1_deseq'].fillna(1) < DEG_PADJ_CUTOFF) & \
                    (df['log2FC_T3vsT1_deseq'].fillna(0).abs() >= DEG_LFC_CUTOFF)
df['is_DEG'] = df['DEG_T2vsT1'] | df['DEG_T3vsT1']

# Direction of expression change
def expr_direction(lfc, padj, cutoff_lfc=1.0, cutoff_padj=0.05):
    if pd.isna(padj) or padj >= cutoff_padj or pd.isna(lfc):
        return 'ns'
    if lfc >= cutoff_lfc:
        return 'up'
    elif lfc <= -cutoff_lfc:
        return 'down'
    return 'ns'

df['expr_dir_T2'] = df.apply(lambda r: expr_direction(r.get('log2FC_T2vsT1_deseq'), r.get('padj_T2vsT1_deseq')), axis=1)
df['expr_dir_T3'] = df.apply(lambda r: expr_direction(r.get('log2FC_T3vsT1_deseq'), r.get('padj_T3vsT1_deseq')), axis=1)

# Coordination assessment
def coordination(methyl_ch, expr_dir, comparison):
    """
    Assess methylation-expression coordination.
    Gained methylation + decreased expression = concordant (repression)
    Lost methylation + increased expression = concordant (de-repression)
    Opposite patterns = discordant
    """
    if methyl_ch == 'none' or methyl_ch == 'constitutive':
        return 'no_methyl_change'
    if expr_dir == 'ns':
        return 'methyl_change_no_expr_change'

    gained = f'gained_{comparison}' in methyl_ch
    lost = f'lost_{comparison}' in methyl_ch

    if gained and expr_dir == 'down':
        return 'concordant_repression'
    elif gained and expr_dir == 'up':
        return 'discordant_gain_up'
    elif lost and expr_dir == 'up':
        return 'concordant_derepression'
    elif lost and expr_dir == 'down':
        return 'discordant_loss_down'
    return 'ambiguous'

df['coordination_T2'] = df.apply(lambda r: coordination(r['methyl_change'], r['expr_dir_T2'], 'T2'), axis=1)
df['coordination_T3'] = df.apply(lambda r: coordination(r['methyl_change'], r['expr_dir_T3'], 'T3'), axis=1)

# Mark if in literature list
df['in_literature_37'] = df['locus_tag'].isin(lit_tf_locus_tags)

# ========== 5. Summary statistics ==========
print("\n" + "=" * 60)
print("Step 5: Summary statistics")
print("=" * 60)

n_total = len(df)
n_methylated = df['has_methylation'].sum()
n_deg = df['is_DEG'].sum()
n_both = ((df['has_methylation']) & (df['is_DEG'])).sum()
n_in_lit = df['in_literature_37'].sum()
n_not_lit = n_total - n_in_lit

print(f"Total regulatory genes: {n_total}")
print(f"With methylation: {n_methylated} ({100*n_methylated/n_total:.1f}%)")
print(f"DEGs (|LFC|>=1, padj<0.05): {n_deg} ({100*n_deg/n_total:.1f}%)")
print(f"Both methylated + DEG: {n_both}")
print(f"In literature 37 list: {n_in_lit}")
print(f"NOT in literature 37 list: {n_not_lit}")

# Methylated but not in literature list
methyl_not_lit = df[(df['has_methylation']) & (~df['in_literature_37'])]
print(f"\nMethylated regulatory genes NOT in literature 37: {len(methyl_not_lit)}")

# Coordination
coord_cols = ['coordination_T2', 'coordination_T3']
for col in coord_cols:
    print(f"\n{col}:")
    vc = df[col].value_counts()
    for k, v in vc.items():
        print(f"  {k}: {v}")

# Coordinated genes
coordinated_types = ['concordant_repression', 'concordant_derepression', 'discordant_gain_up', 'discordant_loss_down']
coord_mask = (df['coordination_T2'].isin(coordinated_types)) | (df['coordination_T3'].isin(coordinated_types))
coordinated_df = df[coord_mask].copy()
print(f"\nTotal genes with methylation-expression coordination: {len(coordinated_df)}")

# ========== 6. Save tables ==========
print("\n" + "=" * 60)
print("Step 6: Saving tables...")
print("=" * 60)

# Table 1: All regulatory genes
cols_to_save = [
    'locus_tag', 'gene_name', 'old_locus_tag', 'product', 'tf_family',
    'start', 'end', 'strand',
    '6mA_T1_count', '6mA_T2_count', '6mA_T3_count',
    '4mC_T1_count', '4mC_T2_count', '4mC_T3_count',
    'total_methyl_T1', 'total_methyl_T2', 'total_methyl_T3',
    'has_methylation', 'methyl_change',
    'log2FC_T2vsT1_deseq', 'padj_T2vsT1_deseq',
    'log2FC_T3vsT1_deseq', 'padj_T3vsT1_deseq',
    'DEG_T2vsT1', 'DEG_T3vsT1', 'is_DEG',
    'expr_dir_T2', 'expr_dir_T3',
    'coordination_T2', 'coordination_T3',
    'in_literature_37'
]
df[cols_to_save].to_csv(os.path.join(TAB_DIR, 'all_regulatory_genes.tsv'), sep='\t', index=False)
print(f"Saved: all_regulatory_genes.tsv ({len(df)} genes)")

# Table 2: Family-level methylation summary
family_summary = df.groupby('tf_family').agg(
    n_genes=('locus_tag', 'count'),
    n_methylated=('has_methylation', 'sum'),
    n_DEG=('is_DEG', 'sum'),
    n_both_methyl_DEG=('locus_tag', lambda x: ((df.loc[x.index, 'has_methylation']) & (df.loc[x.index, 'is_DEG'])).sum()),
    n_in_lit37=('in_literature_37', 'sum'),
).reset_index()
family_summary['methylation_rate'] = (family_summary['n_methylated'] / family_summary['n_genes'] * 100).round(1)
family_summary['DEG_rate'] = (family_summary['n_DEG'] / family_summary['n_genes'] * 100).round(1)
family_summary = family_summary.sort_values('n_genes', ascending=False)
family_summary.to_csv(os.path.join(TAB_DIR, 'regulatory_methylation_summary.tsv'), sep='\t', index=False)
print(f"Saved: regulatory_methylation_summary.tsv ({len(family_summary)} families)")

# Table 3: Coordinated regulatory genes
if len(coordinated_df) > 0:
    coordinated_df[cols_to_save].to_csv(os.path.join(TAB_DIR, 'coordinated_regulatory_genes.tsv'), sep='\t', index=False)
    print(f"Saved: coordinated_regulatory_genes.tsv ({len(coordinated_df)} genes)")
else:
    # Save empty file with header
    pd.DataFrame(columns=cols_to_save).to_csv(os.path.join(TAB_DIR, 'coordinated_regulatory_genes.tsv'), sep='\t', index=False)
    print("Saved: coordinated_regulatory_genes.tsv (0 genes - no coordination found)")

# Table 4: Missing from literature list but methylated/DEG
missing_lit = df[(~df['in_literature_37']) & ((df['has_methylation']) | (df['is_DEG']))].copy()
missing_lit = missing_lit.sort_values(['has_methylation', 'is_DEG'], ascending=[False, False])
missing_lit[cols_to_save].to_csv(os.path.join(TAB_DIR, 'missing_from_literature_list.tsv'), sep='\t', index=False)
print(f"Saved: missing_from_literature_list.tsv ({len(missing_lit)} genes)")

# ========== 7. Notable regulatory genes check ==========
print("\n" + "=" * 60)
print("Step 7: Notable regulatory genes check...")
print("=" * 60)

# Check for well-known TFs that should be in a comprehensive list
notable_names = [
    'whiA', 'whiB', 'whiG', 'whiH', 'whiI', 'whiJ',
    'glnR', 'phoP', 'phoR', 'scbR', 'ramR', 'ramC',
    'abrC1', 'abrC2', 'abrC3', 'argR', 'cutR',
    'relA', 'osaB', 'mtrA', 'mtrB',
    'matR', 'rapA1', 'rapA2',
]

# Search by gene name in GFF
notable_found = []
for nm in notable_names:
    matches = df[df['gene_name'].str.lower() == nm.lower()]
    if len(matches) > 0:
        row = matches.iloc[0]
        notable_found.append({
            'gene_name': nm,
            'locus_tag': row['locus_tag'],
            'old_locus_tag': row['old_locus_tag'],
            'tf_family': row['tf_family'],
            'has_methylation': row['has_methylation'],
            'is_DEG': row['is_DEG'],
            'in_literature_37': row['in_literature_37'],
        })
    else:
        # Try searching in old_locus_tag or by partial name match
        pass

if notable_found:
    notable_df = pd.DataFrame(notable_found)
    print("Notable regulatory genes found in GFF extract:")
    print(notable_df.to_string(index=False))

    # Those not in literature 37 that are noteworthy
    notable_not_lit = notable_df[~notable_df['in_literature_37']]
    if len(notable_not_lit) > 0:
        print(f"\nNotable TFs NOT in literature 37 list:")
        print(notable_not_lit.to_string(index=False))
else:
    print("No notable genes found by gene name (may need SCO lookup)")

# ========== 8. Figures ==========
print("\n" + "=" * 60)
print("Step 8: Generating figures...")
print("=" * 60)

plt.rcParams.update({
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

# --- Figure 1: TF family methylation rate bar chart ---
fig1, ax1 = plt.subplots(figsize=(12, 7))

# Filter families with >= 3 genes for meaningful comparison
fam_plot = family_summary[family_summary['n_genes'] >= 3].copy()
fam_plot = fam_plot.sort_values('methylation_rate', ascending=True)

colors = plt.cm.RdYlBu_r(np.linspace(0.2, 0.8, len(fam_plot)))
bars = ax1.barh(range(len(fam_plot)), fam_plot['methylation_rate'].values, color=colors)

ax1.set_yticks(range(len(fam_plot)))
ax1.set_yticklabels([f"{r['tf_family']} (n={r['n_genes']})" for _, r in fam_plot.iterrows()], fontsize=9)
ax1.set_xlabel('Methylation Rate (%)')
ax1.set_title('Methylation Rate by TF Family (families with n >= 3 genes)')

# Add value labels
for i, (_, row) in enumerate(fam_plot.iterrows()):
    ax1.text(row['methylation_rate'] + 0.5, i,
             f"{row['methylation_rate']:.1f}% ({int(row['n_methylated'])}/{int(row['n_genes'])})",
             va='center', fontsize=8)

ax1.set_xlim(0, max(fam_plot['methylation_rate'].values) * 1.3 if len(fam_plot) > 0 else 100)
plt.tight_layout()

for ext in ['pdf', 'svg']:
    fig1.savefig(os.path.join(FIG_DIR, f'tf_family_methylation_rate.{ext}'), format=ext)
print("Saved: tf_family_methylation_rate.pdf/svg")
plt.close(fig1)

# --- Figure 2: Literature 37 vs all regulatory genes comparison ---
fig2, axes2 = plt.subplots(1, 3, figsize=(14, 5))

# 2a: Methylation rate comparison
lit_df = df[df['in_literature_37']]
nonlit_df = df[~df['in_literature_37']]

categories = ['Literature 37', f'Other Regulatory\n(n={len(nonlit_df)})', f'All Regulatory\n(n={len(df)})']
methyl_rates = [
    100 * lit_df['has_methylation'].sum() / len(lit_df) if len(lit_df) > 0 else 0,
    100 * nonlit_df['has_methylation'].sum() / len(nonlit_df) if len(nonlit_df) > 0 else 0,
    100 * df['has_methylation'].sum() / len(df) if len(df) > 0 else 0,
]

bar_colors = ['#e74c3c', '#3498db', '#2ecc71']
bars2a = axes2[0].bar(categories, methyl_rates, color=bar_colors, edgecolor='black', linewidth=0.5)
axes2[0].set_ylabel('Methylation Rate (%)')
axes2[0].set_title('Methylation Rate')
for bar, rate in zip(bars2a, methyl_rates):
    axes2[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                  f'{rate:.1f}%', ha='center', va='bottom', fontsize=9)

# 2b: DEG rate comparison
deg_rates = [
    100 * lit_df['is_DEG'].sum() / len(lit_df) if len(lit_df) > 0 else 0,
    100 * nonlit_df['is_DEG'].sum() / len(nonlit_df) if len(nonlit_df) > 0 else 0,
    100 * df['is_DEG'].sum() / len(df) if len(df) > 0 else 0,
]
bars2b = axes2[1].bar(categories, deg_rates, color=bar_colors, edgecolor='black', linewidth=0.5)
axes2[1].set_ylabel('DEG Rate (%)')
axes2[1].set_title('DEG Rate (|LFC|>=1, padj<0.05)')
for bar, rate in zip(bars2b, deg_rates):
    axes2[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                  f'{rate:.1f}%', ha='center', va='bottom', fontsize=9)

# 2c: Both methylated + DEG
both_rates = [
    100 * ((lit_df['has_methylation']) & (lit_df['is_DEG'])).sum() / len(lit_df) if len(lit_df) > 0 else 0,
    100 * ((nonlit_df['has_methylation']) & (nonlit_df['is_DEG'])).sum() / len(nonlit_df) if len(nonlit_df) > 0 else 0,
    100 * ((df['has_methylation']) & (df['is_DEG'])).sum() / len(df) if len(df) > 0 else 0,
]
bars2c = axes2[2].bar(categories, both_rates, color=bar_colors, edgecolor='black', linewidth=0.5)
axes2[2].set_ylabel('Rate (%)')
axes2[2].set_title('Methylated + DEG')
for bar, rate in zip(bars2c, both_rates):
    axes2[2].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                  f'{rate:.1f}%', ha='center', va='bottom', fontsize=9)

fig2.suptitle('Literature 37 TFs vs Genome-wide Regulatory Genes', fontsize=13, y=1.02)
plt.tight_layout()

for ext in ['pdf', 'svg']:
    fig2.savefig(os.path.join(FIG_DIR, f'literature_vs_genomewide_comparison.{ext}'), format=ext)
print("Saved: literature_vs_genomewide_comparison.pdf/svg")
plt.close(fig2)

# --- Figure 3: Coordination heatmap (if any coordinated genes exist) ---
# Show all methylated regulatory genes with their expression and coordination status
methyl_genes = df[df['has_methylation']].copy()
if len(methyl_genes) > 0:
    # Prepare heatmap data
    heatmap_data = methyl_genes[['locus_tag', 'gene_name', 'tf_family',
                                  'total_methyl_T1', 'total_methyl_T2', 'total_methyl_T3',
                                  'log2FC_T2vsT1_deseq', 'log2FC_T3vsT1_deseq',
                                  'coordination_T2', 'coordination_T3',
                                  'in_literature_37']].copy()

    # Create display label
    heatmap_data['label'] = heatmap_data.apply(
        lambda r: f"{r['gene_name']} ({r['locus_tag']})" if r['gene_name'] else r['locus_tag'], axis=1)

    # Sort by family then methylation count
    heatmap_data['total_methyl'] = heatmap_data['total_methyl_T1'] + heatmap_data['total_methyl_T2'] + heatmap_data['total_methyl_T3']
    heatmap_data = heatmap_data.sort_values(['tf_family', 'total_methyl'], ascending=[True, False])

    n_genes_hm = len(heatmap_data)
    fig_height = max(6, n_genes_hm * 0.3)
    fig3, axes3 = plt.subplots(1, 2, figsize=(12, fig_height),
                                gridspec_kw={'width_ratios': [1, 1]})

    # Left panel: methylation counts
    methyl_matrix = heatmap_data[['total_methyl_T1', 'total_methyl_T2', 'total_methyl_T3']].values
    im1 = axes3[0].imshow(methyl_matrix, aspect='auto', cmap='YlOrRd', interpolation='nearest')
    axes3[0].set_yticks(range(n_genes_hm))
    axes3[0].set_yticklabels(heatmap_data['label'].values, fontsize=7)
    axes3[0].set_xticks([0, 1, 2])
    axes3[0].set_xticklabels(['T1', 'T2', 'T3'])
    axes3[0].set_title('Methylation Sites (count)')
    plt.colorbar(im1, ax=axes3[0], shrink=0.5)

    # Add text annotations
    for i in range(n_genes_hm):
        for j in range(3):
            val = methyl_matrix[i, j]
            if val > 0:
                axes3[0].text(j, i, str(int(val)), ha='center', va='center', fontsize=6, color='black')

    # Right panel: expression log2FC
    expr_matrix = heatmap_data[['log2FC_T2vsT1_deseq', 'log2FC_T3vsT1_deseq']].fillna(0).values
    vmax = max(abs(expr_matrix.min()), abs(expr_matrix.max()), 2)
    im2 = axes3[1].imshow(expr_matrix, aspect='auto', cmap='RdBu_r',
                           interpolation='nearest', vmin=-vmax, vmax=vmax)
    axes3[1].set_yticks(range(n_genes_hm))
    axes3[1].set_yticklabels(heatmap_data['label'].values, fontsize=7)
    axes3[1].set_xticks([0, 1])
    axes3[1].set_xticklabels(['T2vsT1', 'T3vsT1'])
    axes3[1].set_title('Expression log2FC')
    plt.colorbar(im2, ax=axes3[1], shrink=0.5)

    # Add coordination markers
    for i in range(n_genes_hm):
        for j, col in enumerate(['coordination_T2', 'coordination_T3']):
            coord = heatmap_data.iloc[i][col]
            if 'concordant' in coord:
                axes3[1].text(j, i, '*', ha='center', va='center', fontsize=10,
                             color='gold', fontweight='bold')
            elif 'discordant' in coord:
                axes3[1].text(j, i, 'x', ha='center', va='center', fontsize=8,
                             color='white', fontweight='bold')

    # Mark literature TFs
    for i in range(n_genes_hm):
        if heatmap_data.iloc[i]['in_literature_37']:
            axes3[0].text(-0.7, i, '\u25cf', ha='center', va='center', fontsize=8, color='red',
                         transform=axes3[0].get_yaxis_transform())

    fig3.suptitle('Methylated Regulatory Genes: Methylation & Expression', fontsize=12, y=1.01)
    plt.tight_layout()

    for ext in ['pdf', 'svg']:
        fig3.savefig(os.path.join(FIG_DIR, f'methylated_regulatory_heatmap.{ext}'), format=ext)
    print("Saved: methylated_regulatory_heatmap.pdf/svg")
    plt.close(fig3)
else:
    print("No methylated regulatory genes found - skipping heatmap")

# --- Figure 4: Methylation by type (6mA vs 4mC) across TF families ---
fig4, ax4 = plt.subplots(figsize=(12, 7))

# For each family, count genes with 6mA and 4mC
fam_methyl_type = []
for fam in fam_plot['tf_family'].values:
    fam_df = df[df['tf_family'] == fam]
    n = len(fam_df)
    n_6mA = ((fam_df['6mA_T1_count'].fillna(0) + fam_df['6mA_T2_count'].fillna(0) + fam_df['6mA_T3_count'].fillna(0)) > 0).sum()
    n_4mC = ((fam_df['4mC_T1_count'].fillna(0) + fam_df['4mC_T2_count'].fillna(0) + fam_df['4mC_T3_count'].fillna(0)) > 0).sum()
    fam_methyl_type.append({'tf_family': fam, 'n_genes': n, '6mA_genes': n_6mA, '4mC_genes': n_4mC})

fam_mt = pd.DataFrame(fam_methyl_type)
fam_mt = fam_mt.sort_values('n_genes', ascending=True)

y_pos = range(len(fam_mt))
bar_height = 0.35
ax4.barh([y - bar_height/2 for y in y_pos],
         100 * fam_mt['6mA_genes'] / fam_mt['n_genes'],
         height=bar_height, color='#e74c3c', label='6mA', edgecolor='black', linewidth=0.3)
ax4.barh([y + bar_height/2 for y in y_pos],
         100 * fam_mt['4mC_genes'] / fam_mt['n_genes'],
         height=bar_height, color='#3498db', label='4mC', edgecolor='black', linewidth=0.3)

ax4.set_yticks(list(y_pos))
ax4.set_yticklabels([f"{r['tf_family']} (n={r['n_genes']})" for _, r in fam_mt.iterrows()], fontsize=9)
ax4.set_xlabel('Genes with Modification (%)')
ax4.set_title('6mA vs 4mC by TF Family')
ax4.legend()
plt.tight_layout()

for ext in ['pdf', 'svg']:
    fig4.savefig(os.path.join(FIG_DIR, f'tf_family_6mA_vs_4mC.{ext}'), format=ext)
print("Saved: tf_family_6mA_vs_4mC.pdf/svg")
plt.close(fig4)

# ========== 9. Final summary ==========
print("\n" + "=" * 60)
print("FINAL SUMMARY")
print("=" * 60)

print(f"\nTotal regulatory genes from GFF: {n_total}")
print(f"TF families identified: {len(family_summary)}")
print(f"Genes with methylation: {n_methylated} ({100*n_methylated/n_total:.1f}%)")
print(f"  - 6mA only: {((df[['6mA_T1_count','6mA_T2_count','6mA_T3_count']].fillna(0).sum(axis=1) > 0) & (df[['4mC_T1_count','4mC_T2_count','4mC_T3_count']].fillna(0).sum(axis=1) == 0)).sum()}")
print(f"  - 4mC only: {((df[['6mA_T1_count','6mA_T2_count','6mA_T3_count']].fillna(0).sum(axis=1) == 0) & (df[['4mC_T1_count','4mC_T2_count','4mC_T3_count']].fillna(0).sum(axis=1) > 0)).sum()}")
print(f"  - Both: {((df[['6mA_T1_count','6mA_T2_count','6mA_T3_count']].fillna(0).sum(axis=1) > 0) & (df[['4mC_T1_count','4mC_T2_count','4mC_T3_count']].fillna(0).sum(axis=1) > 0)).sum()}")

print(f"\nDEGs among regulatory genes: {n_deg} ({100*n_deg/n_total:.1f}%)")
print(f"Methylated + DEG: {n_both}")
print(f"Coordinated (concordant or discordant): {len(coordinated_df)}")

# Methylation change breakdown
methyl_change_vc = df[df['has_methylation']]['methyl_change'].value_counts()
print(f"\nMethylation change breakdown (methylated genes only):")
for k, v in methyl_change_vc.items():
    print(f"  {k}: {v}")

# Show all methylated regulatory genes not in lit 37
print(f"\n--- Methylated regulatory genes NOT in literature 37 ({len(methyl_not_lit)}) ---")
if len(methyl_not_lit) > 0:
    display_cols = ['locus_tag', 'gene_name', 'old_locus_tag', 'product', 'tf_family',
                    'total_methyl_T1', 'total_methyl_T2', 'total_methyl_T3',
                    'methyl_change', 'is_DEG', 'expr_dir_T2', 'expr_dir_T3',
                    'coordination_T2', 'coordination_T3']
    print(methyl_not_lit[display_cols].to_string(index=False))

# Show all coordinated genes
if len(coordinated_df) > 0:
    print(f"\n--- Coordinated regulatory genes ({len(coordinated_df)}) ---")
    print(coordinated_df[display_cols].to_string(index=False))

print("\n\nDone! All outputs saved to:")
print(f"  Tables: {TAB_DIR}")
print(f"  Figures: {FIG_DIR}")
