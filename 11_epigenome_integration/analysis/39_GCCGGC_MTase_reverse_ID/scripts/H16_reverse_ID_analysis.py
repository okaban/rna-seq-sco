#!/usr/bin/env python3
"""
H16: T1 GCCGGC MTase Expression-Correlated Reverse Identification
Compute Spearman correlation between MTase expression dynamics and GCCGGC site counts,
identify top candidates, analyze genomic neighborhoods, and cross-validate.
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import warnings
warnings.filterwarnings('ignore')

# Paths
BASE = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/39_GCCGGC_MTase_reverse_ID"
CATALOG = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/33_GCCGGC_RM_identification/tables/M145_MTase_complete_catalog.tsv"
H10_CANDIDATES = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/33_GCCGGC_RM_identification/tables/H10_ranked_candidates.tsv"
BLAST_RESULTS = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/35_GCCGGC_MTase_BLAST/tables/blast_results_Svi27968I_full.txt"
GCCGGC_SITES = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/37_defense_island_GCCGGC/tables/GCCGGC_sites_by_timepoint.tsv"
GFF_FILE = f"{BASE}/data/GCF_000203835.1_ASM20383v1_genomic.gff"

# GCCGGC site counts per timepoint
SITE_COUNTS = [1289, 407, 21]

print("=" * 80)
print("H16: GCCGGC MTase Expression-Correlated Reverse Identification")
print("=" * 80)

# =============================================================================
# 1. Load MTase catalog
# =============================================================================
print("\n--- Step 1: Load MTase catalog ---")
df = pd.read_csv(CATALOG, sep='\t')
print(f"Loaded {len(df)} MTases from catalog")
print(f"Columns: {list(df.columns)}")

# Show category distribution
print(f"\nCategory distribution:")
for cat, cnt in df['category'].value_counts().items():
    print(f"  {cat}: {cnt}")

# =============================================================================
# 2. Compute expression-site correlation (Spearman)
# =============================================================================
print("\n--- Step 2: Compute expression-site Spearman correlation ---")

def compute_spearman(row):
    """Compute Spearman correlation between expression and site counts."""
    expr = [row['T1_mean'], row['T2_mean'], row['T3_mean']]
    # Skip if any NaN
    if any(pd.isna(expr)):
        return np.nan, np.nan
    # If all expression values identical, rho is undefined
    if len(set(expr)) == 1:
        return 0.0, 1.0
    rho, pval = stats.spearmanr(expr, SITE_COUNTS)
    return rho, pval

results = []
for _, row in df.iterrows():
    rho, pval = compute_spearman(row)
    results.append({'rho': rho, 'pval': pval})

corr_df = pd.DataFrame(results)
df['spearman_rho'] = corr_df['rho']
df['spearman_pval'] = corr_df['pval']

# Sort by rho descending
df_sorted = df.sort_values('spearman_rho', ascending=False).reset_index(drop=True)

print(f"\nTop 20 MTases by Spearman rho (expression vs GCCGGC site counts):")
print(f"{'Rank':<5} {'Locus':<14} {'rho':<8} {'T1_expr':<10} {'T2_expr':<10} {'T3_expr':<10} {'Category':<35} {'Product':<50}")
print("-" * 160)
for i, row in df_sorted.head(20).iterrows():
    print(f"{i+1:<5} {row['locus_tag']:<14} {row['spearman_rho']:<8.3f} {row['T1_mean']:<10.1f} {row['T2_mean']:<10.1f} {row['T3_mean']:<10.1f} {str(row['category']):<35} {str(row['product'])[:50]}")

# =============================================================================
# 3. Filter top candidates: rho > 0.8 AND T1 > 50
# =============================================================================
print("\n--- Step 3: Filter candidates (rho > 0.8 AND T1 > 50) ---")
candidates = df_sorted[(df_sorted['spearman_rho'] > 0.8) & (df_sorted['T1_mean'] > 50)].copy()
print(f"Found {len(candidates)} candidates meeting criteria")

# Also show a relaxed filter: rho >= 0.5 AND T1 > 50
relaxed = df_sorted[(df_sorted['spearman_rho'] >= 0.5) & (df_sorted['T1_mean'] > 50)].copy()
print(f"Relaxed filter (rho >= 0.5 AND T1 > 50): {len(relaxed)} candidates")

print(f"\n--- Top candidates (rho > 0.8, T1 > 50) ---")
for _, row in candidates.iterrows():
    print(f"  {row['locus_tag']}: rho={row['spearman_rho']:.3f}, T1={row['T1_mean']:.1f}, T2={row['T2_mean']:.1f}, T3={row['T3_mean']:.1f}")
    print(f"    Product: {row['product']}")
    print(f"    Category: {row['category']}")
    if 'protein_id' in row and pd.notna(row['protein_id']):
        print(f"    Protein: {row['protein_id']}")
    print()

# =============================================================================
# 4. Domain classification: Exclude non-DNA MTases
# =============================================================================
print("\n--- Step 4: Domain classification ---")

# Define non-DNA MTase categories to flag
non_dna_categories = [
    'RNA MTase', 'Protein MTase', 'Small molecule MTase',
    'Histidine MTase', 'DNA repair MTase', 'O-methyltransferase',
    'Terpenoid MTase', 'FkbM family MTase', 'FxLD system MTase'
]

# DNA-relevant categories
dna_categories = [
    'DNA MTase (unspecified)', 'DNA cytosine MTase (Dcm-like, 5mC)',
    'N-6 adenine DNA MTase', 'BREX system MTase'
]

# Unclassified but potentially DNA-related
unclassified = [
    'SAM-dependent MTase (unclassified)', 'MTase (general)',
    'Type 11 MTase'
]

def classify_candidate(row):
    cat = row['category']
    prod = str(row['product']).lower()

    if cat in non_dna_categories:
        return 'NON-DNA', f'Category: {cat}'
    elif cat in dna_categories:
        if 'n-6' in prod or 'adenine' in prod:
            return 'DNA-N6A', f'N-6 adenine specific (not cytosine)'
        elif 'brex' in prod:
            return 'DNA-BREX', f'BREX system (typically m6A)'
        elif 'cytosine' in prod:
            return 'DNA-C', f'DNA cytosine methyltransferase'
        else:
            return 'DNA-UNSPECIFIED', f'DNA MTase, specificity unknown'
    elif cat in unclassified:
        # Check product for clues
        if 'rna' in prod or 'rrna' in prod or 'trna' in prod:
            return 'NON-DNA', f'RNA MTase (from product)'
        elif 'protein' in prod or 'peptide' in prod:
            return 'NON-DNA', f'Protein MTase (from product)'
        elif 'dna' in prod:
            if 'cytosine' in prod:
                return 'DNA-C', f'DNA cytosine MTase (from product)'
            elif 'adenine' in prod or 'n-6' in prod:
                return 'DNA-N6A', f'DNA adenine MTase (from product)'
            else:
                return 'DNA-UNSPECIFIED', f'DNA MTase (from product)'
        else:
            return 'UNCLASSIFIED', f'SAM-dependent, substrate unknown'
    else:
        return 'OTHER', f'Category: {cat}'

# Apply classification to all candidates
for idx, row in candidates.iterrows():
    cls, note = classify_candidate(row)
    candidates.loc[idx, 'domain_class'] = cls
    candidates.loc[idx, 'domain_note'] = note

print("\nCandidate domain classification:")
for _, row in candidates.iterrows():
    print(f"  {row['locus_tag']}: {row['domain_class']} - {row['domain_note']}")

# Prioritize: DNA-C > DNA-UNSPECIFIED > UNCLASSIFIED > DNA-BREX/N6A > NON-DNA
priority_map = {
    'DNA-C': 1, 'DNA-UNSPECIFIED': 2, 'UNCLASSIFIED': 3,
    'DNA-BREX': 4, 'DNA-N6A': 5, 'NON-DNA': 6, 'OTHER': 7
}
candidates['domain_priority'] = candidates['domain_class'].map(priority_map)

# =============================================================================
# 5. Genomic neighborhood analysis (+/- 10 kb)
# =============================================================================
print("\n--- Step 5: Genomic neighborhood analysis ---")

# Parse GFF to get all CDS features
gff_genes = []
with open(GFF_FILE) as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9:
            continue
        if parts[2] != 'CDS':
            continue
        attrs = {}
        for item in parts[8].split(';'):
            if '=' in item:
                k, v = item.split('=', 1)
                attrs[k] = v
        locus = attrs.get('locus_tag', '')
        product = attrs.get('product', '')
        protein_id = attrs.get('protein_id', '')
        gff_genes.append({
            'chrom': parts[0],
            'start': int(parts[3]),
            'end': int(parts[4]),
            'strand': parts[6],
            'locus_tag': locus,
            'product': product,
            'protein_id': protein_id
        })

gff_df = pd.DataFrame(gff_genes)
print(f"Loaded {len(gff_df)} CDS features from GFF")

# R-M related keywords
rm_keywords = [
    'restriction', 'endonuclease', 'nuclease', 'specificity',
    'methyltransferase', 'methylase', 'hsd', 'modification',
    'defense', 'argonaute', 'abortive', 'BREX', 'pgl',
    'cas', 'crispr', 'toxin', 'antitoxin', 'restriction-modification'
]

def find_neighbors(locus_tag, window=10000):
    """Find genes within +/- window bp of a given locus_tag."""
    target = gff_df[gff_df['locus_tag'] == locus_tag]
    if target.empty:
        # Try matching without exact match - partial
        target = gff_df[gff_df['locus_tag'].str.contains(locus_tag.replace('SC_RS', ''), na=False)]
    if target.empty:
        return pd.DataFrame(), None, None

    t = target.iloc[0]
    center = (t['start'] + t['end']) // 2
    neighbors = gff_df[
        (gff_df['chrom'] == t['chrom']) &
        (gff_df['start'] >= t['start'] - window) &
        (gff_df['end'] <= t['end'] + window) &
        (gff_df['locus_tag'] != locus_tag)
    ].copy()

    # Flag R-M related neighbors
    def is_rm_related(product):
        prod_lower = str(product).lower()
        return any(kw in prod_lower for kw in rm_keywords)

    neighbors['rm_related'] = neighbors['product'].apply(is_rm_related)
    return neighbors, t['start'], t['end']

# Analyze neighborhood for all candidates with rho > 0.8 AND T1 > 50
# Also include some key candidates for comprehensive analysis
analysis_targets = list(candidates['locus_tag'].values)
# Add relaxed candidates that are DNA-related
for _, row in relaxed.iterrows():
    cls, _ = classify_candidate(row)
    if cls.startswith('DNA') and row['locus_tag'] not in analysis_targets:
        analysis_targets.append(row['locus_tag'])

print(f"\nAnalyzing neighborhoods for {len(analysis_targets)} targets:")

neighbor_results = []
for lt in analysis_targets:
    neighbors, start, end = find_neighbors(lt)
    if neighbors.empty:
        print(f"\n  {lt}: No CDS neighbors found within 10kb (or locus not in GFF)")
        continue

    rm_neighbors = neighbors[neighbors['rm_related']]
    print(f"\n  {lt} ({start}-{end}):")
    print(f"    Total neighbors: {len(neighbors)}, R-M related: {len(rm_neighbors)}")

    for _, n in rm_neighbors.iterrows():
        print(f"    ** {n['locus_tag']}: {n['product']} ({n['start']}-{n['end']} {n['strand']})")
        neighbor_results.append({
            'target_locus': lt,
            'neighbor_locus': n['locus_tag'],
            'neighbor_product': n['product'],
            'neighbor_start': n['start'],
            'neighbor_end': n['end'],
            'neighbor_strand': n['strand'],
            'neighbor_protein_id': n['protein_id'],
            'rm_related': True
        })

    # Also show all neighbors for context
    if len(rm_neighbors) == 0:
        print(f"    No R-M related genes within 10kb")
        for _, n in neighbors.iterrows():
            print(f"      {n['locus_tag']}: {n['product'][:60]}")

# =============================================================================
# 6. Expression dynamics comparison
# =============================================================================
print("\n--- Step 6: Expression dynamics comparison ---")

# Reference: SC_RS17645 (AAGCCCG, T1=245, 260 sites -> ratio ~1.06 sites/count)
ref_locus = 'SC_RS17645'
ref_row = df[df['locus_tag'] == ref_locus].iloc[0]
ref_ratio = 260 / ref_row['T1_mean']  # AAGCCCG sites / T1 expression
print(f"Reference: {ref_locus} (AAGCCCG)")
print(f"  T1={ref_row['T1_mean']:.1f}, T2={ref_row['T2_mean']:.1f}, T3={ref_row['T3_mean']:.1f}")
print(f"  Sites=260, Ratio (sites/expr) = {ref_ratio:.2f}")

print(f"\nExpression dynamics for top candidates:")
print(f"{'Locus':<14} {'T1_expr':<10} {'T2_expr':<10} {'T3_expr':<10} {'T1_ratio':<12} {'T2_ratio':<12} {'T3_ratio':<12} {'Category'}")
print("-" * 110)

# For all MTases with rho >= 0.5 and T1 > 50 (candidates worth considering)
for _, row in relaxed.iterrows():
    t1_ratio = SITE_COUNTS[0] / row['T1_mean'] if row['T1_mean'] > 0 else float('inf')
    t2_ratio = SITE_COUNTS[1] / row['T2_mean'] if row['T2_mean'] > 0 else float('inf')
    t3_ratio = SITE_COUNTS[2] / row['T3_mean'] if row['T3_mean'] > 0 else float('inf')
    print(f"{row['locus_tag']:<14} {row['T1_mean']:<10.1f} {row['T2_mean']:<10.1f} {row['T3_mean']:<10.1f} {t1_ratio:<12.2f} {t2_ratio:<12.2f} {t3_ratio:<12.2f} {row['category']}")

# =============================================================================
# 7. Cross-validation with existing data
# =============================================================================
print("\n--- Step 7: Cross-validation ---")

# 7a. Check H10 ranked candidates
print("\n7a. Cross-reference with H10 ranked candidates:")
h10 = pd.read_csv(H10_CANDIDATES, sep='\t')
h10_loci = set(h10['locus_tag'].values)
for _, row in candidates.iterrows():
    lt = row['locus_tag']
    if lt in h10_loci:
        h10_row = h10[h10['locus_tag'] == lt].iloc[0]
        print(f"  {lt}: FOUND in H10 (score={h10_row['evidence_score']}, rank=H10)")
    else:
        print(f"  {lt}: NOT in H10 candidates")

# 7b. Check BLAST hits
print("\n7b. Cross-reference with H12 BLAST hits:")
blast_hits = ['WP_254613693.1', 'WP_011031229.1', 'WP_011028829.1', 'WP_011030310.1', 'WP_161270170.1']
blast_loci = {
    'WP_254613693.1': 'SC_RS19770',
    'WP_011031229.1': 'SC_RS36410',
    'WP_011028829.1': 'SC_RS18645',  # approximate
    'WP_011030310.1': 'SC_RS30050',
    'WP_161270170.1': 'SC_RS25950'
}
for _, row in candidates.iterrows():
    lt = row['locus_tag']
    pid = row.get('protein_id', '')
    if pid in blast_hits:
        print(f"  {lt} ({pid}): FOUND in BLAST hits")
    else:
        print(f"  {lt}: NOT in BLAST hits")

# =============================================================================
# 8. Comprehensive ranking table
# =============================================================================
print("\n--- Step 8: Comprehensive candidate ranking ---")

# Create comprehensive scoring for all rho >= 0.5 candidates
comprehensive = relaxed.copy()

# Score components:
# 1. Spearman rho (0-3 points)
# 2. T1 expression adequacy (0-3 points, based on sites/expr ratio)
# 3. Domain classification (0-3 points)
# 4. Genomic neighborhood (0-2 points, bonus for R-M neighbors)
# 5. Cross-validation (0-2 points, bonus for H10/BLAST hits)

def score_candidate(row):
    score = 0
    details = []

    # 1. Spearman rho
    rho = row['spearman_rho']
    if rho == 1.0:
        score += 3
        details.append(f"rho=1.0 (+3)")
    elif rho >= 0.8:
        score += 2
        details.append(f"rho={rho:.2f} (+2)")
    elif rho >= 0.5:
        score += 1
        details.append(f"rho={rho:.2f} (+1)")

    # 2. T1 expression (sites/expr ratio close to reference ~1.06)
    t1 = row['T1_mean']
    ratio = SITE_COUNTS[0] / t1 if t1 > 0 else float('inf')
    if 0.5 <= ratio <= 10:  # reasonable range
        score += 3
        details.append(f"T1={t1:.0f}, ratio={ratio:.1f} (+3)")
    elif 10 < ratio <= 30:
        score += 2
        details.append(f"T1={t1:.0f}, ratio={ratio:.1f} (+2)")
    elif t1 > 50:
        score += 1
        details.append(f"T1={t1:.0f}, ratio={ratio:.1f} (+1)")
    else:
        details.append(f"T1={t1:.0f}, too low (+0)")

    # 3. Domain classification
    cls, note = classify_candidate(row)
    if cls == 'DNA-C':
        score += 3
        details.append(f"DNA cytosine MTase (+3)")
    elif cls == 'DNA-UNSPECIFIED':
        score += 2
        details.append(f"DNA MTase unspecified (+2)")
    elif cls == 'UNCLASSIFIED':
        score += 1
        details.append(f"Unclassified SAM MTase (+1)")
    elif cls == 'DNA-BREX':
        score += 1
        details.append(f"BREX (typically m6A, +1)")
    elif cls == 'DNA-N6A':
        score += 0
        details.append(f"N6-adenine specific (+0)")
    else:
        score += 0
        details.append(f"Non-DNA MTase (+0)")

    # 4. R-M neighbor bonus
    lt = row['locus_tag']
    has_rm_neighbor = any(r['target_locus'] == lt for r in neighbor_results)
    if has_rm_neighbor:
        score += 2
        details.append(f"R-M neighbor (+2)")

    # 5. Cross-validation
    if lt in h10_loci:
        score += 1
        details.append(f"H10 hit (+1)")
    pid = row.get('protein_id', '')
    if pid in blast_hits:
        score += 1
        details.append(f"BLAST hit (+1)")

    return score, '; '.join(details)

scores = []
for idx, row in comprehensive.iterrows():
    s, d = score_candidate(row)
    scores.append({'composite_score': s, 'scoring_details': d})

score_df = pd.DataFrame(scores, index=comprehensive.index)
comprehensive = pd.concat([comprehensive, score_df], axis=1)
comprehensive = comprehensive.sort_values('composite_score', ascending=False)

print(f"\nTop 15 candidates by composite score:")
print(f"{'Rank':<5} {'Locus':<14} {'Score':<7} {'rho':<8} {'T1':<10} {'T2':<10} {'T3':<10} {'Category':<35}")
print("-" * 120)
for rank, (_, row) in enumerate(comprehensive.head(15).iterrows(), 1):
    print(f"{rank:<5} {row['locus_tag']:<14} {row['composite_score']:<7} {row['spearman_rho']:<8.3f} {row['T1_mean']:<10.1f} {row['T2_mean']:<10.1f} {row['T3_mean']:<10.1f} {str(row['category']):<35}")
    print(f"      Details: {row['scoring_details']}")

# =============================================================================
# 9. Save output tables
# =============================================================================
print("\n--- Step 9: Save output tables ---")

# 9a. Full correlation table
out_cols = ['locus_tag', 'start', 'end', 'strand', 'product', 'protein_id',
            'T1_mean', 'T2_mean', 'T3_mean', 'category', 'spearman_rho', 'spearman_pval']
df_sorted[out_cols].to_csv(f"{BASE}/tables/ranked_MTase_correlations.tsv", sep='\t', index=False)
print(f"Saved ranked_MTase_correlations.tsv ({len(df_sorted)} MTases)")

# 9b. Top candidates detail
detail_cols = list(comprehensive.columns)
comprehensive.to_csv(f"{BASE}/tables/top_candidates_detail.tsv", sep='\t', index=False)
print(f"Saved top_candidates_detail.tsv ({len(comprehensive)} candidates)")

# 9c. Genomic neighbors
if neighbor_results:
    pd.DataFrame(neighbor_results).to_csv(f"{BASE}/tables/genomic_neighbors.tsv", sep='\t', index=False)
    print(f"Saved genomic_neighbors.tsv ({len(neighbor_results)} R-M related neighbors)")
else:
    print("No R-M related neighbors found")

# =============================================================================
# 10. Generate figures
# =============================================================================
print("\n--- Step 10: Generate figures ---")

# Figure 1: Spearman rho distribution with candidate highlighting
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# 1a: Scatter: T1 expression vs Spearman rho
ax = axes[0]
colors = []
for _, row in df_sorted.iterrows():
    cls, _ = classify_candidate(row)
    if cls == 'DNA-C':
        colors.append('red')
    elif cls.startswith('DNA'):
        colors.append('orange')
    elif cls == 'UNCLASSIFIED':
        colors.append('steelblue')
    elif cls == 'NON-DNA':
        colors.append('gray')
    else:
        colors.append('lightgray')

ax.scatter(df_sorted['T1_mean'], df_sorted['spearman_rho'],
           c=colors, alpha=0.7, edgecolors='black', linewidths=0.3, s=40)

# Highlight top candidates
for _, row in comprehensive.head(5).iterrows():
    ax.annotate(row['locus_tag'],
                (row['T1_mean'], row['spearman_rho']),
                fontsize=7, fontweight='bold',
                xytext=(5, 5), textcoords='offset points')

ax.axhline(y=0.8, color='red', linestyle='--', alpha=0.5, label='rho=0.8')
ax.axhline(y=0.5, color='orange', linestyle='--', alpha=0.5, label='rho=0.5')
ax.axvline(x=50, color='green', linestyle='--', alpha=0.5, label='T1=50')
ax.set_xlabel('T1 Expression (normalized counts)', fontsize=11)
ax.set_ylabel('Spearman rho\n(expression vs GCCGGC sites)', fontsize=11)
ax.set_title('MTase Expression-Site Correlation', fontsize=12, fontweight='bold')
ax.set_xscale('log')

# Legend
legend_elements = [
    mpatches.Patch(color='red', label='DNA cytosine MTase'),
    mpatches.Patch(color='orange', label='Other DNA MTase'),
    mpatches.Patch(color='steelblue', label='Unclassified SAM MTase'),
    mpatches.Patch(color='gray', label='Non-DNA MTase'),
]
ax.legend(handles=legend_elements, fontsize=8, loc='lower left')

# 1b: Histogram of Spearman rho
ax = axes[1]
ax.hist(df_sorted['spearman_rho'].dropna(), bins=20, color='steelblue', alpha=0.7, edgecolor='black')
ax.axvline(x=0.8, color='red', linestyle='--', label='rho=0.8 threshold')
ax.axvline(x=0.5, color='orange', linestyle='--', label='rho=0.5 threshold')
ax.set_xlabel('Spearman rho', fontsize=11)
ax.set_ylabel('Number of MTases', fontsize=11)
ax.set_title('Distribution of Expression-Site Correlations', fontsize=12, fontweight='bold')
ax.legend(fontsize=9)

# Count and annotate
n_high = len(df_sorted[df_sorted['spearman_rho'] > 0.8])
n_med = len(df_sorted[(df_sorted['spearman_rho'] >= 0.5) & (df_sorted['spearman_rho'] <= 0.8)])
ax.text(0.95, 0.95, f'rho > 0.8: {n_high}\n0.5 <= rho <= 0.8: {n_med}',
        transform=ax.transAxes, ha='right', va='top', fontsize=9,
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

plt.tight_layout()
plt.savefig(f"{BASE}/figures/correlation_scatter.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved correlation_scatter.png")

# Figure 2: Expression timelines for top candidates
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

timepoints = ['T1', 'T2', 'T3']
tp_x = [1, 2, 3]

# 2a: GCCGGC site counts
ax = axes[0, 0]
ax.plot(tp_x, SITE_COUNTS, 'ko-', linewidth=2, markersize=8, label='GCCGGC sites')
ax.set_xticks(tp_x)
ax.set_xticklabels(timepoints)
ax.set_ylabel('Number of 4mC sites')
ax.set_title('GCCGGC 4mC Site Counts', fontweight='bold')
ax.legend()
for i, (x, y) in enumerate(zip(tp_x, SITE_COUNTS)):
    ax.annotate(str(y), (x, y), fontsize=10, fontweight='bold',
                xytext=(0, 10), textcoords='offset points', ha='center')

# 2b: Top candidates expression
ax = axes[0, 1]
top_n = min(8, len(comprehensive))
cmap = plt.cm.tab10
for rank, (_, row) in enumerate(comprehensive.head(top_n).iterrows()):
    expr = [row['T1_mean'], row['T2_mean'], row['T3_mean']]
    label = f"{row['locus_tag']} (rho={row['spearman_rho']:.2f})"
    ax.plot(tp_x, expr, 'o-', color=cmap(rank), linewidth=1.5, markersize=6, label=label)

ax.set_xticks(tp_x)
ax.set_xticklabels(timepoints)
ax.set_ylabel('Expression (normalized counts)')
ax.set_title('Top Candidate Expression Dynamics', fontweight='bold')
ax.legend(fontsize=7, loc='upper right')

# 2c: Normalized expression comparison (scaled to T1=100%)
ax = axes[1, 0]
# Normalize site counts
site_norm = [s / SITE_COUNTS[0] * 100 for s in SITE_COUNTS]
ax.plot(tp_x, site_norm, 'k--', linewidth=3, markersize=10, label='GCCGGC sites', marker='s')

for rank, (_, row) in enumerate(comprehensive.head(5).iterrows()):
    expr = [row['T1_mean'], row['T2_mean'], row['T3_mean']]
    if expr[0] > 0:
        expr_norm = [e / expr[0] * 100 for e in expr]
    else:
        continue
    label = f"{row['locus_tag']}"
    ax.plot(tp_x, expr_norm, 'o-', color=cmap(rank), linewidth=1.5, markersize=6, label=label)

ax.set_xticks(tp_x)
ax.set_xticklabels(timepoints)
ax.set_ylabel('Relative level (T1 = 100%)')
ax.set_title('Normalized Expression vs Site Counts', fontweight='bold')
ax.legend(fontsize=8)
ax.axhline(y=100, color='gray', linestyle=':', alpha=0.3)

# 2d: Sites-per-expression ratio
ax = axes[1, 1]
ref_ratios = [260 / ref_row['T1_mean'], 64 / ref_row['T2_mean'],
              # AAGCCCG T3 reference - use a reasonable value
              0]  # AAGCCCG T3 is largely lost
# Actually compute from known data
aagcccg_sites = [260, 64, 0]  # approximate
for rank, (_, row) in enumerate(comprehensive.head(5).iterrows()):
    ratios = []
    for i, tp in enumerate(['T1_mean', 'T2_mean', 'T3_mean']):
        if row[tp] > 0:
            ratios.append(SITE_COUNTS[i] / row[tp])
        else:
            ratios.append(0)
    label = f"{row['locus_tag']}"
    ax.plot(tp_x, ratios, 'o-', color=cmap(rank), linewidth=1.5, markersize=6, label=label)

# Add AAGCCCG reference
aagcccg_ratios = []
for i, tp in enumerate(['T1_mean', 'T2_mean', 'T3_mean']):
    if ref_row[tp] > 0 and i < len(aagcccg_sites):
        aagcccg_ratios.append(aagcccg_sites[i] / ref_row[tp])
    else:
        aagcccg_ratios.append(0)
ax.plot(tp_x, aagcccg_ratios, 'k--', linewidth=2, marker='s', markersize=8, label='SC_RS17645 (AAGCCCG ref)')

ax.set_xticks(tp_x)
ax.set_xticklabels(timepoints)
ax.set_ylabel('Sites / Expression')
ax.set_title('Site-to-Expression Ratio', fontweight='bold')
ax.legend(fontsize=7)

plt.tight_layout()
plt.savefig(f"{BASE}/figures/expression_timelines.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved expression_timelines.png")

# Figure 3: Genomic neighborhood diagram for top candidates
fig, axes = plt.subplots(len(analysis_targets[:6]), 1, figsize=(16, 3 * min(6, len(analysis_targets))))

if len(analysis_targets) == 1:
    axes = [axes]

for ax_idx, lt in enumerate(analysis_targets[:6]):
    ax = axes[ax_idx] if len(analysis_targets[:6]) > 1 else axes[0]
    neighbors, start, end = find_neighbors(lt, window=10000)

    if neighbors.empty:
        ax.text(0.5, 0.5, f'{lt}: No data', ha='center', va='center', transform=ax.transAxes)
        ax.set_title(f'{lt} Genomic Neighborhood (±10 kb)')
        continue

    # Get target info from catalog
    target_row = df[df['locus_tag'] == lt]
    if not target_row.empty:
        t_start = target_row.iloc[0]['start']
        t_end = target_row.iloc[0]['end']
    else:
        t_start = start
        t_end = end

    center = (t_start + t_end) // 2

    # Draw genes
    all_genes = pd.concat([
        neighbors,
        pd.DataFrame([{'locus_tag': lt, 'start': t_start, 'end': t_end,
                       'strand': '+', 'product': 'TARGET', 'rm_related': False}])
    ])

    for _, gene in all_genes.iterrows():
        g_start = gene['start']
        g_end = gene['end']
        g_strand = gene['strand']
        is_target = (gene['locus_tag'] == lt)
        is_rm = gene.get('rm_related', False)

        if is_target:
            color = 'red'
            alpha = 1.0
        elif is_rm:
            color = 'orange'
            alpha = 0.9
        else:
            color = 'lightblue'
            alpha = 0.6

        y = 0.5 if g_strand == '+' else -0.5
        width = g_end - g_start
        rect = mpatches.FancyBboxPatch((g_start, y - 0.3), width, 0.6,
                                        boxstyle="round,pad=0",
                                        facecolor=color, edgecolor='black',
                                        linewidth=0.5, alpha=alpha)
        ax.add_patch(rect)

        # Label important genes
        if is_target or is_rm:
            label_text = gene['locus_tag']
            ax.text(g_start + width/2, y + 0.5 if g_strand == '+' else y - 0.5,
                   label_text, fontsize=6, ha='center', va='center', fontweight='bold')

    ax.set_xlim(center - 11000, center + 11000)
    ax.set_ylim(-1.5, 1.5)
    ax.axhline(y=0, color='black', linewidth=0.5)

    # Get expression info
    cat_row = df[df['locus_tag'] == lt]
    if not cat_row.empty:
        t1 = cat_row.iloc[0]['T1_mean']
        rho = cat_row.iloc[0]['spearman_rho'] if 'spearman_rho' in cat_row.columns else 'N/A'
        ax.set_title(f'{lt} (T1={t1:.0f}, rho={rho}) - Genomic Neighborhood ±10 kb', fontweight='bold', fontsize=10)
    else:
        ax.set_title(f'{lt} Genomic Neighborhood ±10 kb', fontweight='bold', fontsize=10)

    ax.set_xlabel('Genome position (bp)')
    ax.set_yticks([])

plt.tight_layout()
plt.savefig(f"{BASE}/figures/genomic_neighborhoods.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved genomic_neighborhoods.png")

# =============================================================================
# 11. Final Summary
# =============================================================================
print("\n" + "=" * 80)
print("FINAL SUMMARY")
print("=" * 80)

print(f"\nGCCGGC 4mC site counts: T1={SITE_COUNTS[0]}, T2={SITE_COUNTS[1]}, T3={SITE_COUNTS[2]}")
print(f"Total MTases analyzed: {len(df)}")
print(f"MTases with rho > 0.8: {len(df_sorted[df_sorted['spearman_rho'] > 0.8])}")
print(f"  of which T1 > 50: {len(candidates)}")

# Determine the strongest candidate
top = comprehensive.head(1).iloc[0] if len(comprehensive) > 0 else None
if top is not None:
    print(f"\n*** TOP CANDIDATE ***")
    print(f"  {top['locus_tag']}: {top['product']}")
    print(f"  Category: {top['category']}")
    print(f"  Spearman rho: {top['spearman_rho']:.3f}")
    print(f"  T1={top['T1_mean']:.1f}, T2={top['T2_mean']:.1f}, T3={top['T3_mean']:.1f}")
    print(f"  Composite score: {top['composite_score']}")
    print(f"  Scoring: {top['scoring_details']}")

# List all candidates with score >= 5
print(f"\n--- All candidates with composite score >= 5 ---")
for _, row in comprehensive[comprehensive['composite_score'] >= 5].iterrows():
    print(f"  {row['locus_tag']}: score={row['composite_score']}, rho={row['spearman_rho']:.3f}, T1={row['T1_mean']:.1f}, cat={row['category']}")

# Key insight about required T1 expression
print(f"\n--- Expression level analysis ---")
print(f"If the responsible MTase has ~1 site/count (like AAGCCCG reference),")
print(f"  T1 expression should be ~1289 to explain 1289 sites")
print(f"  T2 expression should be ~407 to explain 407 sites")
print(f"  T3 expression should be ~21 to explain 21 sites")
print(f"\nMTases with T1 > 500 AND rho=1.0:")
high_expr_corr = df_sorted[(df_sorted['T1_mean'] > 500) & (df_sorted['spearman_rho'] == 1.0)]
for _, row in high_expr_corr.iterrows():
    cls, _ = classify_candidate(row)
    print(f"  {row['locus_tag']}: T1={row['T1_mean']:.1f}, T2={row['T2_mean']:.1f}, T3={row['T3_mean']:.1f}, cat={row['category']}, class={cls}")

print(f"\nMTases with T1 > 1000 AND rho=1.0:")
very_high = df_sorted[(df_sorted['T1_mean'] > 1000) & (df_sorted['spearman_rho'] == 1.0)]
for _, row in very_high.iterrows():
    cls, _ = classify_candidate(row)
    print(f"  {row['locus_tag']}: T1={row['T1_mean']:.1f}, T2={row['T2_mean']:.1f}, T3={row['T3_mean']:.1f}, cat={row['category']}, class={cls}")

print("\n--- Analysis complete ---")
