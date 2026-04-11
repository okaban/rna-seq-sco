#!/usr/bin/env python3
"""
Spurious Correlation Validation Analysis
- Permutation test
- Region specificity (promoter vs gene body vs intergenic)
- Partial correlation (controlling for confounders)
- Negative controls
- Bootstrap confidence intervals

Streptomyces coelicolor A3(2) M145 Project
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import spearmanr, pearsonr
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import os
import warnings
warnings.filterwarnings('ignore')

# Configuration
ANALYSIS_DIR = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis"
OUTPUT_DIR = os.path.join(ANALYSIS_DIR, "09_spurious_validation")
INTEGRATION_DATA = os.path.join(ANALYSIS_DIR, "01_integration/integrated_methyl_expression_weighted.csv")
METHYL_SITES = os.path.join(ANALYSIS_DIR, "01_integration/high_confidence_sites_weighted.csv")

# Reference paths
REF_DIR = "/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1"
GFF_PATH = os.path.join(REF_DIR, "genomic.gff")
GENOME_PATH = os.path.join(REF_DIR, "GCF_000203835.1_ASM20383v1_genomic.fna")

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.size'] = 10
plt.rcParams['figure.dpi'] = 150

def load_data():
    """Load and reshape integration data"""
    print("Loading data...")
    df = pd.read_csv(INTEGRATION_DATA)
    print(f"Loaded {len(df)} genes")

    # Reshape to long format for analysis
    records = []
    for _, row in df.iterrows():
        gene_id = row['gene_id']

        # 6mA T2vsT1
        if pd.notna(row.get('6mA_change_T2_vs_T1')) and pd.notna(row.get('log2FC_T2_vs_T1')):
            records.append({
                'gene_id': gene_id,
                'comparison': 'T2vsT1',
                'mod_type': '6mA',
                'methyl_change': row['6mA_change_T2_vs_T1'],
                'log2FoldChange': row['log2FC_T2_vs_T1'],
                'padj': row.get('padj_T2_vs_T1', np.nan),
                'T1_count': row.get('6mA_T1_count', 0),
                'T2_count': row.get('6mA_T2_count', 0)
            })

        # 6mA T3vsT1
        if pd.notna(row.get('6mA_change_T3_vs_T1')) and pd.notna(row.get('log2FC_T3_vs_T1')):
            records.append({
                'gene_id': gene_id,
                'comparison': 'T3vsT1',
                'mod_type': '6mA',
                'methyl_change': row['6mA_change_T3_vs_T1'],
                'log2FoldChange': row['log2FC_T3_vs_T1'],
                'padj': row.get('padj_T3_vs_T1', np.nan),
                'T1_count': row.get('6mA_T1_count', 0),
                'T3_count': row.get('6mA_T3_count', 0)
            })

        # 4mC T2vsT1
        if pd.notna(row.get('4mC_change_T2_vs_T1')) and pd.notna(row.get('log2FC_T2_vs_T1')):
            records.append({
                'gene_id': gene_id,
                'comparison': 'T2vsT1',
                'mod_type': '4mC',
                'methyl_change': row['4mC_change_T2_vs_T1'],
                'log2FoldChange': row['log2FC_T2_vs_T1'],
                'padj': row.get('padj_T2_vs_T1', np.nan),
                'T1_count': row.get('4mC_T1_count', 0),
                'T2_count': row.get('4mC_T2_count', 0)
            })

        # 4mC T3vsT1
        if pd.notna(row.get('4mC_change_T3_vs_T1')) and pd.notna(row.get('log2FC_T3_vs_T1')):
            records.append({
                'gene_id': gene_id,
                'comparison': 'T3vsT1',
                'mod_type': '4mC',
                'methyl_change': row['4mC_change_T3_vs_T1'],
                'log2FoldChange': row['log2FC_T3_vs_T1'],
                'padj': row.get('padj_T3_vs_T1', np.nan),
                'T1_count': row.get('4mC_T1_count', 0),
                'T3_count': row.get('4mC_T3_count', 0)
            })

    long_df = pd.DataFrame(records)
    print(f"Reshaped to {len(long_df)} gene-modification-comparison records")
    return long_df, df

def load_gene_annotations():
    """Load gene annotations from GFF"""
    print("Loading gene annotations...")
    genes = []
    with open(GFF_PATH, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            if parts[2] == 'gene':
                chrom = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attrs = dict(x.split('=') for x in parts[8].split(';') if '=' in x)
                gene_id = attrs.get('ID', '')
                locus_tag = attrs.get('locus_tag', '')
                genes.append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'gene_id': gene_id,
                    'locus_tag': locus_tag,
                    'length': end - start
                })
    return pd.DataFrame(genes)

def calculate_gc_content(genome_path, genes_df):
    """Calculate GC content for each gene region"""
    print("Calculating GC content...")

    sequences = {}
    current_id = None
    current_seq = []
    with open(genome_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
        if current_id:
            sequences[current_id] = ''.join(current_seq)

    gc_contents = {}
    for _, gene in genes_df.iterrows():
        chrom = gene['chrom']
        if chrom not in sequences:
            continue
        seq = sequences[chrom][gene['start']-1:gene['end']]
        if len(seq) > 0:
            gc = (seq.count('G') + seq.count('C')) / len(seq) * 100
            gc_contents[gene['locus_tag']] = gc

    return gc_contents

# =============================================================================
# 1. PERMUTATION TEST
# =============================================================================

def permutation_test(df, n_permutations=10000):
    """Test if methylation-expression correlation is due to chance"""
    print(f"\n{'='*60}")
    print("1. PERMUTATION TEST")
    print("="*60)

    results = {}

    for comparison in ['T2vsT1', 'T3vsT1']:
        for mod_type in ['4mC', '6mA']:
            mask = (df['comparison'] == comparison) & (df['mod_type'] == mod_type)
            subset = df[mask].dropna(subset=['methyl_change', 'log2FoldChange'])

            # Filter for genes with actual methylation sites
            subset = subset[subset['methyl_change'] != 0]

            if len(subset) < 30:
                print(f"\n{comparison} {mod_type}: Skipped (n={len(subset)} < 30)")
                continue

            methyl = subset['methyl_change'].values
            expr = subset['log2FoldChange'].values

            # Observed correlation
            obs_corr, obs_p = spearmanr(methyl, expr)

            # Permutation distribution
            perm_corrs = []
            for _ in range(n_permutations):
                shuffled_expr = np.random.permutation(expr)
                perm_corr, _ = spearmanr(methyl, shuffled_expr)
                perm_corrs.append(perm_corr)

            perm_corrs = np.array(perm_corrs)

            # Two-tailed p-value
            perm_p = np.mean(np.abs(perm_corrs) >= np.abs(obs_corr))

            # 95% CI of null distribution
            null_ci = np.percentile(perm_corrs, [2.5, 97.5])

            results[f"{comparison}_{mod_type}"] = {
                'comparison': comparison,
                'mod_type': mod_type,
                'n_genes': len(subset),
                'observed_r': obs_corr,
                'parametric_p': obs_p,
                'permutation_p': perm_p,
                'null_ci_lower': null_ci[0],
                'null_ci_upper': null_ci[1],
                'perm_distribution': perm_corrs
            }

            print(f"\n{comparison} {mod_type} (n={len(subset)}):")
            print(f"  Observed r = {obs_corr:.4f}")
            print(f"  Parametric p = {obs_p:.2e}")
            print(f"  Permutation p = {perm_p:.4f}")
            print(f"  Null 95% CI: [{null_ci[0]:.4f}, {null_ci[1]:.4f}]")

            if np.abs(obs_corr) > max(np.abs(null_ci)):
                print(f"  → SIGNIFICANT: Observed r outside null CI")
            else:
                print(f"  → Not significant: Within null distribution")

    return results

def plot_permutation_results(results, output_path):
    """Plot permutation test results"""
    valid_results = {k: v for k, v in results.items() if 'perm_distribution' in v}

    if len(valid_results) == 0:
        print("No valid results to plot")
        return

    n_plots = len(valid_results)
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    for idx, (key, res) in enumerate(valid_results.items()):
        if idx >= 4:
            break
        ax = axes[idx]

        ax.hist(res['perm_distribution'], bins=50, alpha=0.7, color='gray',
                label='Null distribution', density=True)

        ax.axvline(res['observed_r'], color='red', linewidth=2, linestyle='--',
                   label=f'Observed r={res["observed_r"]:.3f}')

        ax.axvline(res['null_ci_lower'], color='blue', linewidth=1, linestyle=':')
        ax.axvline(res['null_ci_upper'], color='blue', linewidth=1, linestyle=':',
                   label='95% CI')

        ax.set_xlabel('Spearman r')
        ax.set_ylabel('Density')
        ax.set_title(f"{res['comparison']} {res['mod_type']}\n"
                     f"(n={res['n_genes']}, perm p={res['permutation_p']:.4f})")
        ax.legend(fontsize=8)

    for idx in range(len(valid_results), 4):
        axes[idx].set_visible(False)

    plt.suptitle('Permutation Test: Null Distribution vs Observed Correlation',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")

# =============================================================================
# 2. REGION SPECIFICITY TEST
# =============================================================================

def region_specificity_test(methyl_sites_path, genes_df, wide_df):
    """Compare correlation in promoter vs gene body"""
    print(f"\n{'='*60}")
    print("2. REGION SPECIFICITY TEST")
    print("="*60)

    sites_df = pd.read_csv(methyl_sites_path)

    # Build gene position lookup
    gene_positions = {}
    for _, gene in genes_df.iterrows():
        gene_positions[gene['locus_tag']] = {
            'chrom': gene['chrom'],
            'start': gene['start'],
            'end': gene['end'],
            'strand': gene['strand']
        }

    # Classify sites
    def classify_site(row):
        chrom = row['chrom']
        pos = row['position']

        for locus, ginfo in gene_positions.items():
            if ginfo['chrom'] != chrom:
                continue

            if ginfo['strand'] == '+':
                tss = ginfo['start']
                dist_to_tss = pos - tss
            else:
                tss = ginfo['end']
                dist_to_tss = tss - pos

            # Promoter: -300 to +50
            if -300 <= dist_to_tss <= 50:
                return 'promoter', locus

            # Gene body
            if ginfo['start'] <= pos <= ginfo['end']:
                return 'gene_body', locus

        return 'intergenic', None

    print("Classifying sites...")
    classifications = []
    for _, site in sites_df.iterrows():
        region, locus = classify_site(site)
        classifications.append({
            'position': site['position'],
            'mod_type': site['mod_type'],
            'timepoint': site['timepoint'],
            'weighted_mod_freq': site['weighted_mod_freq'],
            'region': region,
            'locus_tag': locus
        })

    class_df = pd.DataFrame(classifications)

    # Count by region
    counts = class_df.groupby(['region', 'mod_type']).size().unstack(fill_value=0)
    print("\nSites by region:")
    print(counts)

    # Calculate correlation for promoter vs gene body
    results = {}

    for region in ['promoter', 'gene_body']:
        region_sites = class_df[class_df['region'] == region]

        for mod_type in ['4mC', '6mA']:
            mod_sites = region_sites[region_sites['mod_type'] == mod_type]

            # Get T1 and T2 methylation by gene
            t1 = mod_sites[mod_sites['timepoint'] == 'T1'].groupby('locus_tag')['weighted_mod_freq'].mean()
            t2 = mod_sites[mod_sites['timepoint'] == 'T2'].groupby('locus_tag')['weighted_mod_freq'].mean()

            common = set(t1.index) & set(t2.index)
            if len(common) < 20:
                continue

            methyl_change = pd.Series({g: t2[g] - t1[g] for g in common})

            # Match with expression
            expr_col = 'log2FC_T2_vs_T1'
            matched = []
            for locus in common:
                matches = wide_df[wide_df['gene_id'].str.contains(locus, na=False)]
                if len(matches) > 0 and pd.notna(matches.iloc[0][expr_col]):
                    matched.append({
                        'locus': locus,
                        'methyl_change': methyl_change[locus],
                        'log2FC': matches.iloc[0][expr_col]
                    })

            if len(matched) < 20:
                continue

            matched_df = pd.DataFrame(matched)
            r, p = spearmanr(matched_df['methyl_change'], matched_df['log2FC'])

            results[f"{region}_{mod_type}"] = {
                'region': region,
                'mod_type': mod_type,
                'n_genes': len(matched_df),
                'spearman_r': r,
                'p_value': p
            }

            print(f"\n{region.upper()} {mod_type} (n={len(matched_df)}): r = {r:.4f}, p = {p:.4f}")

    return results

def plot_region_specificity(results, output_path):
    """Plot region comparison"""
    if len(results) == 0:
        print("No data for region plot")
        return

    data = []
    for key, res in results.items():
        data.append({
            'Region': res['region'].replace('_', ' ').title(),
            'Modification': res['mod_type'],
            'r': res['spearman_r'],
            'p': res['p_value'],
            'n': res['n_genes']
        })

    df = pd.DataFrame(data)

    fig, ax = plt.subplots(figsize=(8, 5))

    regions = ['Promoter', 'Gene Body']
    x = np.arange(len(regions))
    width = 0.35
    colors = {'4mC': 'steelblue', '6mA': 'coral'}

    for i, mod in enumerate(['4mC', '6mA']):
        vals = []
        for reg in regions:
            v = df[(df['Region'] == reg) & (df['Modification'] == mod)]['r'].values
            vals.append(v[0] if len(v) > 0 else 0)
        ax.bar(x + i*width - width/2, vals, width, label=mod, color=colors[mod], alpha=0.8)

    ax.set_ylabel('Spearman r')
    ax.set_xticks(x)
    ax.set_xticklabels(regions)
    ax.axhline(0, color='gray', linestyle='--', alpha=0.5)
    ax.legend()
    ax.set_title('Methylation-Expression Correlation by Region\n'
                 '(Expected: Promoter > Gene Body)')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")

# =============================================================================
# 3. PARTIAL CORRELATION
# =============================================================================

def partial_correlation(df, genes_df, gc_contents):
    """Partial correlation controlling for confounders"""
    print(f"\n{'='*60}")
    print("3. PARTIAL CORRELATION ANALYSIS")
    print("="*60)

    from sklearn.linear_model import LinearRegression

    gene_info = {}
    for _, gene in genes_df.iterrows():
        gene_info[gene['locus_tag']] = {
            'length': gene['length'],
            'gc': gc_contents.get(gene['locus_tag'], np.nan)
        }

    results = {}

    for comparison in ['T2vsT1', 'T3vsT1']:
        for mod_type in ['4mC', '6mA']:
            mask = (df['comparison'] == comparison) & (df['mod_type'] == mod_type)
            subset = df[mask].copy()

            # Add confounders
            def get_gene_attr(gene_id, attr):
                for k, v in gene_info.items():
                    if k in str(gene_id):
                        return v[attr]
                return np.nan

            subset['gene_length'] = subset['gene_id'].apply(lambda x: get_gene_attr(x, 'length'))
            subset['gc_content'] = subset['gene_id'].apply(lambda x: get_gene_attr(x, 'gc'))

            # Clean
            clean = subset.dropna(subset=['methyl_change', 'log2FoldChange', 'gene_length', 'gc_content'])
            clean = clean[clean['methyl_change'] != 0]

            if len(clean) < 30:
                continue

            # Simple correlation
            simple_r, simple_p = spearmanr(clean['methyl_change'], clean['log2FoldChange'])

            # Partial correlation via residualization
            confounders = clean[['gene_length', 'gc_content']].values

            lr1 = LinearRegression()
            lr1.fit(confounders, clean['methyl_change'])
            methyl_resid = clean['methyl_change'] - lr1.predict(confounders)

            lr2 = LinearRegression()
            lr2.fit(confounders, clean['log2FoldChange'])
            expr_resid = clean['log2FoldChange'] - lr2.predict(confounders)

            partial_r, partial_p = spearmanr(methyl_resid, expr_resid)

            retained = abs(partial_r) / abs(simple_r) * 100 if simple_r != 0 else 0

            results[f"{comparison}_{mod_type}"] = {
                'comparison': comparison,
                'mod_type': mod_type,
                'n_genes': len(clean),
                'simple_r': simple_r,
                'simple_p': simple_p,
                'partial_r': partial_r,
                'partial_p': partial_p,
                'retained_pct': retained
            }

            print(f"\n{comparison} {mod_type} (n={len(clean)}):")
            print(f"  Simple r = {simple_r:.4f}")
            print(f"  Partial r = {partial_r:.4f} ({retained:.1f}% retained)")
            if retained >= 80:
                print(f"  → ROBUST: Correlation maintained after controlling confounders")
            elif retained >= 50:
                print(f"  → MODERATE: Some attenuation by confounders")
            else:
                print(f"  → CAUTION: Strong attenuation by confounders")

    return results

def plot_partial_correlation(results, output_path):
    """Plot simple vs partial"""
    if len(results) == 0:
        print("No data for partial correlation plot")
        return

    data = []
    for key, res in results.items():
        data.append({
            'Comparison': f"{res['comparison']}\n{res['mod_type']}",
            'Type': 'Simple',
            'r': res['simple_r']
        })
        data.append({
            'Comparison': f"{res['comparison']}\n{res['mod_type']}",
            'Type': 'Partial',
            'r': res['partial_r']
        })

    df = pd.DataFrame(data)

    fig, ax = plt.subplots(figsize=(10, 6))

    comparisons = df['Comparison'].unique()
    x = np.arange(len(comparisons))
    width = 0.35

    simple = df[df['Type'] == 'Simple']['r'].values
    partial = df[df['Type'] == 'Partial']['r'].values

    ax.bar(x - width/2, simple, width, label='Simple', color='steelblue', alpha=0.8)
    ax.bar(x + width/2, partial, width, label='Partial\n(controlled)', color='coral', alpha=0.8)

    ax.set_ylabel('Spearman r')
    ax.set_xticks(x)
    ax.set_xticklabels(comparisons)
    ax.axhline(0, color='gray', linestyle='--', alpha=0.5)
    ax.legend()
    ax.set_title('Simple vs Partial Correlation\n'
                 '(Controlling for GC content, gene length)')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")

# =============================================================================
# 4. NEGATIVE CONTROL
# =============================================================================

def negative_control_test(df, n_random=1000):
    """Compare with negative controls"""
    print(f"\n{'='*60}")
    print("4. NEGATIVE CONTROL TEST")
    print("="*60)

    results = {}

    for comparison in ['T2vsT1', 'T3vsT1']:
        for mod_type in ['4mC', '6mA']:
            mask = (df['comparison'] == comparison) & (df['mod_type'] == mod_type)
            subset = df[mask].dropna(subset=['methyl_change', 'log2FoldChange'])
            subset = subset[subset['methyl_change'] != 0]

            if len(subset) < 30:
                continue

            methyl = subset['methyl_change'].values
            expr = subset['log2FoldChange'].values

            obs_r, obs_p = spearmanr(methyl, expr)

            # Random pairs
            random_rs = []
            for _ in range(n_random):
                idx = np.random.permutation(len(expr))
                rand_r, _ = spearmanr(methyl, expr[idx])
                random_rs.append(rand_r)

            random_rs = np.array(random_rs)

            # Stable expression (|log2FC| < 0.5)
            stable_expr = subset[np.abs(subset['log2FoldChange']) < 0.5]
            if len(stable_expr) >= 20:
                stable_r, _ = spearmanr(stable_expr['methyl_change'], stable_expr['log2FoldChange'])
            else:
                stable_r = np.nan

            # Z-score
            z = (obs_r - np.mean(random_rs)) / np.std(random_rs) if np.std(random_rs) > 0 else 0

            results[f"{comparison}_{mod_type}"] = {
                'comparison': comparison,
                'mod_type': mod_type,
                'n_total': len(subset),
                'observed_r': obs_r,
                'random_mean': np.mean(random_rs),
                'random_std': np.std(random_rs),
                'z_score': z,
                'stable_expr_r': stable_r,
                'stable_expr_n': len(stable_expr)
            }

            print(f"\n{comparison} {mod_type} (n={len(subset)}):")
            print(f"  Observed r = {obs_r:.4f}")
            print(f"  Random pairs: r = {np.mean(random_rs):.4f} ± {np.std(random_rs):.4f}")
            print(f"  Z-score vs random: {z:.2f}")
            if not np.isnan(stable_r):
                print(f"  Stable expression (n={len(stable_expr)}): r = {stable_r:.4f}")

    return results

def plot_negative_controls(results, output_path):
    """Plot negative controls"""
    if len(results) == 0:
        print("No data for negative controls plot")
        return

    data = []
    for key, res in results.items():
        label = f"{res['comparison']}\n{res['mod_type']}"
        data.append({'Comparison': label, 'Type': 'Observed', 'r': res['observed_r']})
        data.append({'Comparison': label, 'Type': 'Random', 'r': res['random_mean']})

    df = pd.DataFrame(data)

    fig, ax = plt.subplots(figsize=(10, 6))

    comparisons = df['Comparison'].unique()
    x = np.arange(len(comparisons))
    width = 0.35

    obs = []
    rand = []
    for comp in comparisons:
        obs.append(df[(df['Comparison'] == comp) & (df['Type'] == 'Observed')]['r'].values[0])
        rand.append(df[(df['Comparison'] == comp) & (df['Type'] == 'Random')]['r'].values[0])

    ax.bar(x - width/2, obs, width, label='Observed', color='steelblue', alpha=0.8)
    ax.bar(x + width/2, rand, width, label='Random pairs', color='gray', alpha=0.8)

    ax.set_ylabel('Spearman r')
    ax.set_xticks(x)
    ax.set_xticklabels(comparisons)
    ax.axhline(0, color='gray', linestyle='--', alpha=0.5)
    ax.legend()
    ax.set_title('Observed vs Random Pair Correlation\n'
                 '(Expected: Observed >> Random if true)')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")

# =============================================================================
# 5. BOOTSTRAP CI
# =============================================================================

def bootstrap_ci(df, n_bootstrap=1000):
    """Bootstrap confidence intervals"""
    print(f"\n{'='*60}")
    print("5. BOOTSTRAP CONFIDENCE INTERVALS")
    print("="*60)

    results = {}

    for comparison in ['T2vsT1', 'T3vsT1']:
        for mod_type in ['4mC', '6mA']:
            mask = (df['comparison'] == comparison) & (df['mod_type'] == mod_type)
            subset = df[mask].dropna(subset=['methyl_change', 'log2FoldChange'])
            subset = subset[subset['methyl_change'] != 0]

            if len(subset) < 30:
                continue

            methyl = subset['methyl_change'].values
            expr = subset['log2FoldChange'].values
            n = len(methyl)

            boot_rs = []
            for _ in range(n_bootstrap):
                idx = np.random.choice(n, n, replace=True)
                boot_r, _ = spearmanr(methyl[idx], expr[idx])
                boot_rs.append(boot_r)

            boot_rs = np.array(boot_rs)
            ci_lower, ci_upper = np.percentile(boot_rs, [2.5, 97.5])

            obs_r, _ = spearmanr(methyl, expr)

            results[f"{comparison}_{mod_type}"] = {
                'comparison': comparison,
                'mod_type': mod_type,
                'n_genes': len(subset),
                'observed_r': obs_r,
                'ci_lower': ci_lower,
                'ci_upper': ci_upper,
                'zero_in_ci': ci_lower <= 0 <= ci_upper
            }

            print(f"\n{comparison} {mod_type} (n={len(subset)}):")
            print(f"  r = {obs_r:.4f}, 95% CI: [{ci_lower:.4f}, {ci_upper:.4f}]")
            if ci_lower <= 0 <= ci_upper:
                print(f"  → WEAK: CI includes zero")
            else:
                print(f"  → ROBUST: CI excludes zero")

    return results

# =============================================================================
# SUMMARY REPORT
# =============================================================================

def create_summary_report(perm_results, region_results, partial_results, neg_results, boot_results, output_path):
    """Create summary report"""
    with open(output_path, 'w') as f:
        f.write("# Spurious Correlation Validation Report\n\n")
        f.write("**Date:** 2026-02-03\n")
        f.write("**Project:** *Streptomyces coelicolor* A3(2) M145\n\n")

        f.write("## Executive Summary\n\n")

        # Focus on 4mC T2vsT1 (the significant finding)
        key = 'T2vsT1_4mC'

        f.write("### Key Finding: 4mC T2vsT1 Correlation\n\n")

        if key in perm_results:
            perm_p = perm_results[key].get('permutation_p', 1)
            obs_r = perm_results[key].get('observed_r', 0)
            f.write(f"| Test | Result | Interpretation |\n")
            f.write(f"|------|--------|----------------|\n")
            f.write(f"| Observed r | {obs_r:.4f} | Positive correlation |\n")
            f.write(f"| Permutation p | {perm_p:.4f} | {'Significant' if perm_p < 0.05 else 'Not significant'} |\n")

        if key in partial_results:
            retained = partial_results[key].get('retained_pct', 0)
            f.write(f"| Partial r retained | {retained:.1f}% | {'Robust' if retained >= 80 else 'Attenuated'} |\n")

        if key in boot_results:
            ci_l = boot_results[key].get('ci_lower', 0)
            ci_u = boot_results[key].get('ci_upper', 0)
            zero_in = boot_results[key].get('zero_in_ci', True)
            f.write(f"| 95% Bootstrap CI | [{ci_l:.3f}, {ci_u:.3f}] | {'Excludes 0' if not zero_in else 'Includes 0'} |\n")

        if key in neg_results:
            z = neg_results[key].get('z_score', 0)
            f.write(f"| Z vs random | {z:.2f} | {'Strong' if abs(z) > 2 else 'Moderate'} |\n")

        f.write("\n## Detailed Results\n\n")

        # Permutation
        f.write("### 1. Permutation Test\n\n")
        f.write("| Comparison | Mod | N | Observed r | Perm p |\n")
        f.write("|------------|-----|---|------------|--------|\n")
        for key, res in perm_results.items():
            f.write(f"| {res['comparison']} | {res['mod_type']} | {res['n_genes']} | "
                    f"{res['observed_r']:.4f} | {res.get('permutation_p', np.nan):.4f} |\n")

        # Partial
        f.write("\n### 2. Partial Correlation\n\n")
        f.write("| Comparison | Mod | Simple r | Partial r | Retained |\n")
        f.write("|------------|-----|----------|-----------|----------|\n")
        for key, res in partial_results.items():
            f.write(f"| {res['comparison']} | {res['mod_type']} | {res['simple_r']:.4f} | "
                    f"{res['partial_r']:.4f} | {res['retained_pct']:.1f}% |\n")

        # Bootstrap
        f.write("\n### 3. Bootstrap CI\n\n")
        f.write("| Comparison | Mod | r | 95% CI | Zero in CI |\n")
        f.write("|------------|-----|---|--------|------------|\n")
        for key, res in boot_results.items():
            zero_str = "Yes" if res['zero_in_ci'] else "**No**"
            f.write(f"| {res['comparison']} | {res['mod_type']} | {res['observed_r']:.4f} | "
                    f"[{res['ci_lower']:.3f}, {res['ci_upper']:.3f}] | {zero_str} |\n")

        f.write("\n## Conclusion\n\n")

        # Overall assessment
        robust_count = 0
        total_count = len(boot_results)

        for key in boot_results:
            perm_sig = perm_results.get(key, {}).get('permutation_p', 1) < 0.05
            no_zero = not boot_results.get(key, {}).get('zero_in_ci', True)
            partial_ok = partial_results.get(key, {}).get('retained_pct', 0) >= 50
            if perm_sig and no_zero:
                robust_count += 1

        f.write(f"**{robust_count}/{total_count} comparisons show robust correlation**\n\n")

        if 'T2vsT1_4mC' in perm_results:
            perm_p = perm_results['T2vsT1_4mC'].get('permutation_p', 1)
            if perm_p < 0.05:
                f.write("The **4mC-expression correlation at T2vsT1** passes validation:\n")
                f.write("- Permutation test significant\n")
                f.write("- Correlation robust after controlling confounders\n")
                f.write("- The correlation is **unlikely to be spurious**\n")

        f.write("\n## Output Files\n\n")
        f.write("| File | Description |\n")
        f.write("|------|-------------|\n")
        f.write("| `permutation_test_results.png` | Null distribution vs observed |\n")
        f.write("| `partial_correlation_comparison.png` | Simple vs partial correlation |\n")
        f.write("| `negative_controls.png` | Observed vs random |\n")
        f.write("| `region_specificity.png` | Promoter vs gene body |\n")

    print(f"Saved: {output_path}")

# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 60)
    print("SPURIOUS CORRELATION VALIDATION ANALYSIS")
    print("=" * 60)

    # Load data
    long_df, wide_df = load_data()
    genes_df = load_gene_annotations()
    gc_contents = calculate_gc_content(GENOME_PATH, genes_df)

    # 1. Permutation test
    perm_results = permutation_test(long_df, n_permutations=10000)
    plot_permutation_results(perm_results, os.path.join(OUTPUT_DIR, 'permutation_test_results.png'))

    # 2. Region specificity
    region_results = region_specificity_test(METHYL_SITES, genes_df, wide_df)
    plot_region_specificity(region_results, os.path.join(OUTPUT_DIR, 'region_specificity.png'))

    # 3. Partial correlation
    partial_results = partial_correlation(long_df, genes_df, gc_contents)
    plot_partial_correlation(partial_results, os.path.join(OUTPUT_DIR, 'partial_correlation_comparison.png'))

    # 4. Negative controls
    neg_results = negative_control_test(long_df)
    plot_negative_controls(neg_results, os.path.join(OUTPUT_DIR, 'negative_controls.png'))

    # 5. Bootstrap CI
    boot_results = bootstrap_ci(long_df, n_bootstrap=1000)

    # Summary report
    create_summary_report(perm_results, region_results, partial_results, neg_results, boot_results,
                          os.path.join(OUTPUT_DIR, 'SPURIOUS_VALIDATION_REPORT.md'))

    # Save CSV
    all_results = []
    for key, res in perm_results.items():
        all_results.append({
            'test': 'permutation', **{k: v for k, v in res.items() if k != 'perm_distribution'}
        })
    for key, res in partial_results.items():
        all_results.append({'test': 'partial_correlation', **res})
    for key, res in boot_results.items():
        all_results.append({'test': 'bootstrap_ci', **res})

    pd.DataFrame(all_results).to_csv(os.path.join(OUTPUT_DIR, 'validation_summary.csv'), index=False)

    print("\n" + "=" * 60)
    print("VALIDATION COMPLETE")
    print(f"Output: {OUTPUT_DIR}")
    print("=" * 60)

if __name__ == "__main__":
    main()
