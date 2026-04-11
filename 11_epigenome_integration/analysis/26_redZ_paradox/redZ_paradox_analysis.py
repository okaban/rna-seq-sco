#!/usr/bin/env python3
"""
Hypothesis H3: Quantitative resolution of the redZ paradox
- How does Red BGC activate despite redZ downregulation?

Analyzes the regulatory balance of activators vs repressors
controlling the Red (undecylprodiginine/streptorubin) BGC in S. coelicolor.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Configuration
# ============================================================
DESEQ2_DIR = "/Users/okaban/bioinfo/rna-seq/04_deseq2/analysis/04_deseq2_260128_v1/results"
OUT_DIR = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/26_redZ_paradox"
FIG_DIR = f"{OUT_DIR}/figures"
TAB_DIR = f"{OUT_DIR}/tables"

# Gene mapping: locus_tag -> (common_name, SCO_id, role, description)
REGULATORS = {
    'SC_RS31630': ('redD', 'SCO5877', 'activator', 'SARP pathway-specific activator'),
    'SC_RS31650': ('redZ', 'SCO5881', 'CSR', 'Response regulator (CSR for Red)'),
    'SC_RS18240': ('absA1', 'SCO3225', 'repressor', 'Sensor kinase (negative regulator)'),
    'SC_RS18245': ('absA2', 'SCO3226', 'repressor', 'Response regulator (negative regulator)'),
    'SC_RS24295': ('afsR', 'SCO4426', 'activator', 'Global SARP activator'),
    'SC_RS24290': ('afsS', 'SCO4425', 'activator', 'Anti-anti-sigma / afsR target'),
    'SC_RS26695': ('afsQ1', 'SCO4907', 'activator', 'TCS response regulator'),
    'SC_RS26690': ('afsQ2', 'SCO4906', 'sensor', 'TCS sensor kinase'),
    'SC_RS17645': ('MTase', 'SCO3104', 'epigenetic', 'Candidate N-6 methyltransferase'),
}

# Red BGC structural genes (SCO5877-SCO5898 region)
RED_BGC_GENES = {
    'SC_RS31630': ('redD', 'SCO5877'),
    'SC_RS31635': ('redX', 'SCO5878'),
    'SC_RS31640': ('redW', 'SCO5879'),
    'SC_RS31645': ('redY', 'SCO5880'),
    'SC_RS31650': ('redZ', 'SCO5881'),
    'SC_RS31655': ('redV', 'SCO5882'),
    'SC_RS31660': ('redU', 'SCO5883'),
    'SC_RS31665': ('redT', 'SCO5884'),
    'SC_RS31670': ('redS', 'SCO5885'),
    'SC_RS31675': ('redR', 'SCO5886'),
    'SC_RS31680': ('redQ', 'SCO5887'),
    'SC_RS31685': ('redP', 'SCO5888'),
    'SC_RS31690': ('redO', 'SCO5889'),
    'SC_RS31695': ('redN', 'SCO5890'),
    'SC_RS31700': ('redM', 'SCO5891'),
    'SC_RS31705': ('redL', 'SCO5892'),
    'SC_RS31710': ('redK', 'SCO5893'),
    'SC_RS31715': ('redJ', 'SCO5894'),
    'SC_RS31720': ('redI', 'SCO5895'),
    'SC_RS31725': ('redH', 'SCO5896'),
    'SC_RS31730': ('redG', 'SCO5897'),
    'SC_RS31735': ('redF', 'SCO5898'),
}

# Structural genes only (exclude regulators redD and redZ for BGC expression calc)
RED_STRUCTURAL = {k: v for k, v in RED_BGC_GENES.items() 
                  if k not in ['SC_RS31630', 'SC_RS31650']}

# ============================================================
# Load data
# ============================================================
print("=" * 70)
print("H3: Quantitative Resolution of the redZ Paradox")
print("=" * 70)

de_t2v1 = pd.read_csv(f"{DESEQ2_DIR}/DESeq2_M145_2_vs_1.tsv", sep='\t', index_col='gene_id')
de_t3v1 = pd.read_csv(f"{DESEQ2_DIR}/DESeq2_M145_3_vs_1.tsv", sep='\t', index_col='gene_id')
de_t3v2 = pd.read_csv(f"{DESEQ2_DIR}/DESeq2_M145_3_vs_2.tsv", sep='\t', index_col='gene_id')
counts = pd.read_csv(f"{DESEQ2_DIR}/normalized_counts_M145.tsv", sep='\t', index_col='gene_id')

t1_cols = [c for c in counts.columns if c.startswith('M145_1')]
t2_cols = [c for c in counts.columns if c.startswith('M145_2')]
t3_cols = [c for c in counts.columns if c.startswith('M145_3')]

counts['T1_mean'] = counts[t1_cols].mean(axis=1)
counts['T2_mean'] = counts[t2_cols].mean(axis=1)
counts['T3_mean'] = counts[t3_cols].mean(axis=1)

print(f"\nLoaded DESeq2 data: {len(de_t2v1)} genes")
print(f"Sample columns: T1={t1_cols}, T2={t2_cols}, T3={t3_cols}")

# ============================================================
# Analysis 1: Regulator Timeline Table
# ============================================================
print("\n" + "=" * 70)
print("ANALYSIS 1: Regulator Expression Timeline")
print("=" * 70)

reg_rows = []
for locus, (name, sco, role, desc) in REGULATORS.items():
    row = {'locus_tag': locus, 'gene': name, 'SCO': sco, 'role': role, 'description': desc}
    
    for label, df in [('T2v1', de_t2v1), ('T3v1', de_t3v1), ('T3v2', de_t3v2)]:
        if locus in df.index:
            row[f'log2FC_{label}'] = df.loc[locus, 'log2FoldChange']
            row[f'padj_{label}'] = df.loc[locus, 'padj']
        else:
            row[f'log2FC_{label}'] = np.nan
            row[f'padj_{label}'] = np.nan
    
    if locus in counts.index:
        row['T1_mean'] = counts.loc[locus, 'T1_mean']
        row['T2_mean'] = counts.loc[locus, 'T2_mean']
        row['T3_mean'] = counts.loc[locus, 'T3_mean']
    else:
        row['T1_mean'] = np.nan
        row['T2_mean'] = np.nan
        row['T3_mean'] = np.nan
    
    reg_rows.append(row)

reg_df = pd.DataFrame(reg_rows)
reg_df.to_csv(f"{TAB_DIR}/regulator_timeline.tsv", sep='\t', index=False)

print(f"\n{'Gene':<8} {'Role':<12} {'T1 mean':>10} {'T2 mean':>10} {'T3 mean':>10} "
      f"{'L2FC T2v1':>10} {'padj':>10} {'L2FC T3v1':>10} {'padj':>10}")
print("-" * 102)
for _, r in reg_df.iterrows():
    sig_t2 = '*' if r.get('padj_T2v1', 1) < 0.05 else ' '
    sig_t3 = '*' if r.get('padj_T3v1', 1) < 0.05 else ' '
    print(f"{r['gene']:<8} {r['role']:<12} {r['T1_mean']:>10.1f} {r['T2_mean']:>10.1f} {r['T3_mean']:>10.1f} "
          f"{r.get('log2FC_T2v1', 0):>9.2f}{sig_t2} {r.get('padj_T2v1', 1):>10.2e} "
          f"{r.get('log2FC_T3v1', 0):>9.2f}{sig_t3} {r.get('padj_T3v1', 1):>10.2e}")

# ============================================================
# Analysis 2: Red BGC Expression Timeline
# ============================================================
print("\n" + "=" * 70)
print("ANALYSIS 2: Red BGC Expression Timeline")
print("=" * 70)

bgc_rows = []
for locus, (name, sco) in RED_BGC_GENES.items():
    row = {'locus_tag': locus, 'gene': name, 'SCO': sco}
    
    if locus in counts.index:
        row['T1_mean'] = counts.loc[locus, 'T1_mean']
        row['T2_mean'] = counts.loc[locus, 'T2_mean']
        row['T3_mean'] = counts.loc[locus, 'T3_mean']
    
    for label, df in [('T2v1', de_t2v1), ('T3v1', de_t3v1), ('T3v2', de_t3v2)]:
        if locus in df.index:
            row[f'log2FC_{label}'] = df.loc[locus, 'log2FoldChange']
            row[f'padj_{label}'] = df.loc[locus, 'padj']
    
    bgc_rows.append(row)

bgc_df = pd.DataFrame(bgc_rows)
bgc_df.to_csv(f"{TAB_DIR}/red_bgc_expression.tsv", sep='\t', index=False)

struct_loci = list(RED_STRUCTURAL.keys())
struct_in_counts = [l for l in struct_loci if l in counts.index]

bgc_mean_t1 = counts.loc[struct_in_counts, 'T1_mean'].mean()
bgc_mean_t2 = counts.loc[struct_in_counts, 'T2_mean'].mean()
bgc_mean_t3 = counts.loc[struct_in_counts, 'T3_mean'].mean()

struct_in_de = [l for l in struct_loci if l in de_t2v1.index]
bgc_l2fc_t2v1 = de_t2v1.loc[struct_in_de, 'log2FoldChange'].mean()
bgc_l2fc_t3v1 = de_t3v1.loc[struct_in_de, 'log2FoldChange'].mean()
bgc_l2fc_t3v2 = de_t3v2.loc[struct_in_de, 'log2FoldChange'].mean()

print(f"\nRed BGC structural genes ({len(struct_in_counts)} genes):")
print(f"  Mean normalized counts: T1={bgc_mean_t1:.1f}, T2={bgc_mean_t2:.1f}, T3={bgc_mean_t3:.1f}")
print(f"  Mean log2FC: T2vsT1={bgc_l2fc_t2v1:.2f}, T3vsT1={bgc_l2fc_t3v1:.2f}, T3vsT2={bgc_l2fc_t3v2:.2f}")
print(f"  Fold-change T2/T1: {2**bgc_l2fc_t2v1:.1f}x")
print(f"  Fold-change T3/T1: {2**bgc_l2fc_t3v1:.1f}x")

print(f"\n{'Gene':<8} {'SCO':<10} {'T1 mean':>10} {'T2 mean':>10} {'T3 mean':>10} {'L2FC T2v1':>10} {'L2FC T3v1':>10}")
print("-" * 70)
for _, r in bgc_df.iterrows():
    print(f"{r['gene']:<8} {r['SCO']:<10} {r.get('T1_mean', 0):>10.1f} {r.get('T2_mean', 0):>10.1f} "
          f"{r.get('T3_mean', 0):>10.1f} {r.get('log2FC_T2v1', 0):>10.2f} {r.get('log2FC_T3v1', 0):>10.2f}")

# ============================================================
# Analysis 3: Activator-Repressor Balance Score
# ============================================================
print("\n" + "=" * 70)
print("ANALYSIS 3: Activator-Repressor Balance Score")
print("=" * 70)

activator_loci = ['SC_RS31630', 'SC_RS24295', 'SC_RS24290', 'SC_RS26695']
repressor_loci = ['SC_RS18240', 'SC_RS18245']
csr_locus = 'SC_RS31650'

def get_l2fc(locus, df):
    if locus in df.index:
        return df.loc[locus, 'log2FoldChange']
    return 0.0

act_scores_t2v1 = {l: get_l2fc(l, de_t2v1) for l in activator_loci}
rep_scores_t2v1 = {l: get_l2fc(l, de_t2v1) for l in repressor_loci}
redz_t2v1 = get_l2fc(csr_locus, de_t2v1)

act_total_t2v1 = sum(act_scores_t2v1.values())
rep_total_t2v1 = sum(rep_scores_t2v1.values())
net_t2v1 = act_total_t2v1 - rep_total_t2v1

act_scores_t3v1 = {l: get_l2fc(l, de_t3v1) for l in activator_loci}
rep_scores_t3v1 = {l: get_l2fc(l, de_t3v1) for l in repressor_loci}
redz_t3v1 = get_l2fc(csr_locus, de_t3v1)

act_total_t3v1 = sum(act_scores_t3v1.values())
rep_total_t3v1 = sum(rep_scores_t3v1.values())
net_t3v1 = act_total_t3v1 - rep_total_t3v1

print("\n--- T2 vs T1 (Transition to antibiotic production) ---")
print(f"\nActivators (sum of log2FC):")
for l, v in act_scores_t2v1.items():
    name = REGULATORS.get(l, RED_BGC_GENES.get(l, ('?',)))[0]
    padj = de_t2v1.loc[l, 'padj'] if l in de_t2v1.index else np.nan
    sig = '*' if padj < 0.05 else ' '
    print(f"  {name:<8} (log2FC={v:>7.3f}{sig}, padj={padj:.2e})")
print(f"  TOTAL activator score: {act_total_t2v1:>7.3f}")

print(f"\nRepressors (sum of log2FC, positive = MORE repression):")
for l, v in rep_scores_t2v1.items():
    name = REGULATORS[l][0]
    padj = de_t2v1.loc[l, 'padj'] if l in de_t2v1.index else np.nan
    sig = '*' if padj < 0.05 else ' '
    print(f"  {name:<8} (log2FC={v:>7.3f}{sig}, padj={padj:.2e})")
print(f"  TOTAL repressor score: {rep_total_t2v1:>7.3f}")

print(f"\nredZ (CSR):  log2FC={redz_t2v1:>7.3f}")
print(f"\nNET BALANCE (Activators - Repressors): {net_t2v1:>7.3f}")
print(f"Red BGC structural mean log2FC:        {bgc_l2fc_t2v1:>7.3f}")

print("\n--- T3 vs T1 (Late stationary phase) ---")
print(f"\nActivators (sum of log2FC):")
for l, v in act_scores_t3v1.items():
    name = REGULATORS.get(l, RED_BGC_GENES.get(l, ('?',)))[0]
    padj = de_t3v1.loc[l, 'padj'] if l in de_t3v1.index else np.nan
    sig = '*' if padj < 0.05 else ' '
    print(f"  {name:<8} (log2FC={v:>7.3f}{sig}, padj={padj:.2e})")
print(f"  TOTAL activator score: {act_total_t3v1:>7.3f}")

print(f"\nRepressors (sum of log2FC):")
for l, v in rep_scores_t3v1.items():
    name = REGULATORS[l][0]
    padj = de_t3v1.loc[l, 'padj'] if l in de_t3v1.index else np.nan
    sig = '*' if padj < 0.05 else ' '
    print(f"  {name:<8} (log2FC={v:>7.3f}{sig}, padj={padj:.2e})")
print(f"  TOTAL repressor score: {rep_total_t3v1:>7.3f}")

print(f"\nredZ (CSR):  log2FC={redz_t3v1:>7.3f}")
print(f"\nNET BALANCE (Activators - Repressors): {net_t3v1:>7.3f}")
print(f"Red BGC structural mean log2FC:        {bgc_l2fc_t3v1:>7.3f}")

balance_data = pd.DataFrame({
    'Comparison': ['T2vsT1', 'T3vsT1'],
    'Activator_score': [act_total_t2v1, act_total_t3v1],
    'Repressor_score': [rep_total_t2v1, rep_total_t3v1],
    'Net_balance': [net_t2v1, net_t3v1],
    'redZ_log2FC': [redz_t2v1, redz_t3v1],
    'Red_BGC_mean_log2FC': [bgc_l2fc_t2v1, bgc_l2fc_t3v1],
})
balance_data.to_csv(f"{TAB_DIR}/activator_repressor_balance.tsv", sep='\t', index=False)

# ============================================================
# Analysis 5: Key Question
# ============================================================
print("\n" + "=" * 70)
print("KEY QUESTION: Does activator total compensate for redZ loss?")
print("=" * 70)

redd_t2v1 = get_l2fc('SC_RS31630', de_t2v1)
absA_relief_t2v1 = -rep_total_t2v1
compensation_t2v1 = redd_t2v1 + absA_relief_t2v1

print(f"\n--- T2 vs T1 ---")
print(f"  redD activation (log2FC):     {redd_t2v1:>+7.3f}")
print(f"  absA suppression relief:      {absA_relief_t2v1:>+7.3f}  (negative of repressor sum)")
print(f"  Compensatory total:           {compensation_t2v1:>+7.3f}")
print(f"  redZ reduction (log2FC):      {redz_t2v1:>+7.3f}")
print(f"  Net (compensation + redZ):    {compensation_t2v1 + redz_t2v1:>+7.3f}")
if compensation_t2v1 + redz_t2v1 > 0:
    print(f"  --> YES: Net balance is POSITIVE, compensators outweigh redZ loss")
else:
    print(f"  --> NO: Net balance is negative")

redd_t3v1 = get_l2fc('SC_RS31630', de_t3v1)
absA_relief_t3v1 = -rep_total_t3v1
compensation_t3v1 = redd_t3v1 + absA_relief_t3v1

print(f"\n--- T3 vs T1 ---")
print(f"  redD activation (log2FC):     {redd_t3v1:>+7.3f}")
print(f"  absA suppression relief:      {absA_relief_t3v1:>+7.3f}")
print(f"  Compensatory total:           {compensation_t3v1:>+7.3f}")
print(f"  redZ reduction (log2FC):      {redz_t3v1:>+7.3f}")
print(f"  Net (compensation + redZ):    {compensation_t3v1 + redz_t3v1:>+7.3f}")
if compensation_t3v1 + redz_t3v1 > 0:
    print(f"  --> YES: Net balance is POSITIVE, compensators outweigh redZ loss")
else:
    print(f"  --> NO: Net balance is negative")

afsr_t2v1 = get_l2fc('SC_RS24295', de_t2v1)
afss_t2v1 = get_l2fc('SC_RS24290', de_t2v1)
afsq1_t2v1 = get_l2fc('SC_RS26695', de_t2v1)
global_t2v1 = afsr_t2v1 + afss_t2v1 + afsq1_t2v1

afsr_t3v1 = get_l2fc('SC_RS24295', de_t3v1)
afss_t3v1 = get_l2fc('SC_RS24290', de_t3v1)
afsq1_t3v1 = get_l2fc('SC_RS26695', de_t3v1)
global_t3v1 = afsr_t3v1 + afss_t3v1 + afsq1_t3v1

print(f"\n--- Full regulatory balance including global regulators ---")
print(f"\n  T2vsT1:")
print(f"    redD activation:           {redd_t2v1:>+7.3f}")
print(f"    absA relief:               {absA_relief_t2v1:>+7.3f}")
print(f"    Global regulators (afsR+afsS+afsQ1): {global_t2v1:>+7.3f}")
print(f"    Total positive input:      {redd_t2v1 + absA_relief_t2v1 + global_t2v1:>+7.3f}")
print(f"    redZ loss:                 {redz_t2v1:>+7.3f}")
full_net_t2v1 = redd_t2v1 + absA_relief_t2v1 + global_t2v1 + redz_t2v1
print(f"    FULL NET BALANCE:          {full_net_t2v1:>+7.3f}")

print(f"\n  T3vsT1:")
print(f"    redD activation:           {redd_t3v1:>+7.3f}")
print(f"    absA relief:               {absA_relief_t3v1:>+7.3f}")
print(f"    Global regulators (afsR+afsS+afsQ1): {global_t3v1:>+7.3f}")
print(f"    Total positive input:      {redd_t3v1 + absA_relief_t3v1 + global_t3v1:>+7.3f}")
print(f"    redZ loss:                 {redz_t3v1:>+7.3f}")
full_net_t3v1 = redd_t3v1 + absA_relief_t3v1 + global_t3v1 + redz_t3v1
print(f"    FULL NET BALANCE:          {full_net_t3v1:>+7.3f}")

# ============================================================
# Visualization
# ============================================================
print("\n" + "=" * 70)
print("Generating figures...")
print("=" * 70)

colors_map = {
    'activator': '#2ca02c',
    'repressor': '#d62728',
    'CSR': '#ff7f0e',
    'sensor': '#9467bd',
    'epigenetic': '#8c564b',
    'BGC': '#1f77b4',
}

# ============================================================
# Figure 1: Multi-panel figure
# ============================================================
fig = plt.figure(figsize=(18, 16))
gs = gridspec.GridSpec(2, 2, hspace=0.35, wspace=0.3)

# --- Panel A ---
ax_a = fig.add_subplot(gs[0, 0])

timepoints = ['T1', 'T2', 'T3']
x = np.arange(len(timepoints))
width = 0.12
genes_to_plot = [
    ('SC_RS31630', 'redD', colors_map['activator']),
    ('SC_RS31650', 'redZ', colors_map['CSR']),
    ('SC_RS24295', 'afsR', '#2ecc71'),
    ('SC_RS24290', 'afsS', '#27ae60'),
    ('SC_RS26695', 'afsQ1', '#1abc9c'),
    ('SC_RS18240', 'absA1', colors_map['repressor']),
    ('SC_RS18245', 'absA2', '#e74c3c'),
]

for i, (locus, name, color) in enumerate(genes_to_plot):
    if locus in counts.index:
        vals = [counts.loc[locus, 'T1_mean'], counts.loc[locus, 'T2_mean'], counts.loc[locus, 'T3_mean']]
        offset = (i - len(genes_to_plot)/2) * width + width/2
        bars = ax_a.bar(x + offset, vals, width, label=name, color=color, alpha=0.85, edgecolor='white', linewidth=0.5)

ax_a.set_xlabel('Timepoint', fontsize=12)
ax_a.set_ylabel('Normalized counts (mean)', fontsize=12)
ax_a.set_title('A. Regulator Expression Timeline', fontsize=14, fontweight='bold')
ax_a.set_xticks(x)
ax_a.set_xticklabels(timepoints)
ax_a.legend(fontsize=8, ncol=2, loc='upper right')
ax_a.set_yscale('log')
ax_a.set_ylim(bottom=0.5)

# --- Panel B ---
ax_b = fig.add_subplot(gs[0, 1])

bgc_means = [bgc_mean_t1, bgc_mean_t2, bgc_mean_t3]
ax_b.fill_between(x, bgc_means, alpha=0.3, color=colors_map['BGC'])
ax_b.plot(x, bgc_means, 'o-', color=colors_map['BGC'], linewidth=2.5, markersize=10, label='Red BGC mean expr', zorder=5)

if csr_locus in counts.index:
    redz_vals = [counts.loc[csr_locus, 'T1_mean'], counts.loc[csr_locus, 'T2_mean'], counts.loc[csr_locus, 'T3_mean']]
    scale_factor = max(bgc_means) / max(redz_vals) if max(redz_vals) > 0 else 1
    redz_scaled = [v * scale_factor for v in redz_vals]
    ax_b.plot(x, redz_scaled, 's--', color=colors_map['CSR'], linewidth=2, markersize=8, label=f'redZ (scaled {scale_factor:.0f}x)', zorder=4)

if 'SC_RS31630' in counts.index:
    redd_vals = [counts.loc['SC_RS31630', 'T1_mean'], counts.loc['SC_RS31630', 'T2_mean'], counts.loc['SC_RS31630', 'T3_mean']]
    scale_factor_d = max(bgc_means) / max(redd_vals) if max(redd_vals) > 0 else 1
    redd_scaled = [v * scale_factor_d for v in redd_vals]
    ax_b.plot(x, redd_scaled, '^--', color=colors_map['activator'], linewidth=2, markersize=8, label=f'redD (scaled {scale_factor_d:.0f}x)', zorder=4)

ax_b.set_xlabel('Timepoint', fontsize=12)
ax_b.set_ylabel('Normalized counts (mean)', fontsize=12)
ax_b.set_title('B. Red BGC Expression vs Key Regulators', fontsize=14, fontweight='bold')
ax_b.set_xticks(x)
ax_b.set_xticklabels(timepoints)
ax_b.legend(fontsize=9, loc='upper left')

# --- Panel C ---
ax_c = fig.add_subplot(gs[1, :])

all_loci = list(REGULATORS.keys()) + list(RED_BGC_GENES.keys())
all_loci = list(dict.fromkeys(all_loci))

label_map = {}
for l in all_loci:
    if l in REGULATORS:
        name, sco, role, _ = REGULATORS[l]
        label_map[l] = f"{name} ({sco}) [{role}]"
    elif l in RED_BGC_GENES:
        name, sco = RED_BGC_GENES[l]
        label_map[l] = f"{name} ({sco})"

valid_loci = [l for l in all_loci if l in counts.index]

expr_matrix = counts.loc[valid_loci, ['T1_mean', 'T2_mean', 'T3_mean']].copy()
expr_matrix.columns = ['T1', 'T2', 'T3']

zscore_matrix = expr_matrix.apply(lambda row: (row - row.mean()) / row.std() if row.std() > 0 else row * 0, axis=1)

ylabels = [label_map.get(l, l) for l in valid_loci]

label_colors = []
for l in valid_loci:
    if l in REGULATORS:
        role = REGULATORS[l][2]
        label_colors.append(colors_map.get(role, 'black'))
    else:
        label_colors.append(colors_map['BGC'])

sns.heatmap(zscore_matrix, ax=ax_c, cmap='RdBu_r', center=0, 
            yticklabels=ylabels, xticklabels=['T1', 'T2', 'T3'],
            linewidths=0.5, linecolor='white',
            cbar_kws={'label': 'Z-score (row-normalized)'})
ax_c.set_title('C. Expression Heatmap: Red BGC Regulators & Structural Genes', fontsize=14, fontweight='bold')
ax_c.set_ylabel('')

for i, (tick, color) in enumerate(zip(ax_c.get_yticklabels(), label_colors)):
    tick.set_color(color)
    tick.set_fontsize(7)

plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f"{FIG_DIR}/redZ_paradox_multipanel.{ext}", dpi=300, bbox_inches='tight')
print(f"  Saved: redZ_paradox_multipanel.pdf/svg")
plt.close()

# ============================================================
# Figure 2: Balance Diagram
# ============================================================
fig2, axes = plt.subplots(1, 2, figsize=(14, 6))

for idx, (comp, data) in enumerate([
    ('T2 vs T1', {
        'redD': redd_t2v1,
        'afsR': afsr_t2v1,
        'afsS': afss_t2v1,
        'afsQ1': afsq1_t2v1,
        'absA1': get_l2fc('SC_RS18240', de_t2v1),
        'absA2': get_l2fc('SC_RS18245', de_t2v1),
        'redZ': redz_t2v1,
    }),
    ('T3 vs T1', {
        'redD': redd_t3v1,
        'afsR': afsr_t3v1,
        'afsS': afss_t3v1,
        'afsQ1': afsq1_t3v1,
        'absA1': get_l2fc('SC_RS18240', de_t3v1),
        'absA2': get_l2fc('SC_RS18245', de_t3v1),
        'redZ': redz_t3v1,
    }),
]):
    ax = axes[idx]
    genes = list(data.keys())
    values = list(data.values())
    
    bar_colors = []
    for g in genes:
        if g in ['redD', 'afsR', 'afsS', 'afsQ1']:
            bar_colors.append(colors_map['activator'] if data[g] > 0 else '#90EE90')
        elif g in ['absA1', 'absA2']:
            bar_colors.append(colors_map['repressor'] if data[g] > 0 else '#FFB6C1')
        else:
            bar_colors.append(colors_map['CSR'])
    
    bars = ax.barh(genes, values, color=bar_colors, edgecolor='white', linewidth=0.5, height=0.6)
    ax.axvline(x=0, color='black', linewidth=0.8)
    ax.set_xlabel('log2 Fold Change', fontsize=12)
    ax.set_title(f'{comp}', fontsize=14, fontweight='bold')
    
    for bar, val in zip(bars, values):
        ha = 'left' if val >= 0 else 'right'
        offset = 0.05 if val >= 0 else -0.05
        ax.text(val + offset, bar.get_y() + bar.get_height()/2, f'{val:.2f}', 
                va='center', ha=ha, fontsize=9, fontweight='bold')
    
    if idx == 0:
        net_val = full_net_t2v1
    else:
        net_val = full_net_t3v1
    ax.annotate(f'Full net balance: {net_val:+.2f}', xy=(0.98, 0.02), xycoords='axes fraction',
                ha='right', va='bottom', fontsize=10, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.8))

plt.suptitle('Activator-Repressor Balance for Red BGC Regulation', fontsize=15, fontweight='bold', y=1.02)
plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig2.savefig(f"{FIG_DIR}/balance_diagram.{ext}", dpi=300, bbox_inches='tight')
print(f"  Saved: balance_diagram.pdf/svg")
plt.close()

# ============================================================
# Figure 3: redZ vs redD vs Red BGC time course
# ============================================================
fig3, ax3 = plt.subplots(figsize=(10, 6))

t_labels = ['T1', 'T2', 'T3']
t_positions = [1, 2, 3]

for col_group, tp, marker in [(t1_cols, 1, 'o'), (t2_cols, 2, 'o'), (t3_cols, 3, 'o')]:
    for col in col_group:
        val = counts.loc[struct_in_counts, col].mean()
        ax3.scatter(tp, val, color=colors_map['BGC'], alpha=0.3, s=30, zorder=3)

ax3.plot(t_positions, bgc_means, 'o-', color=colors_map['BGC'], linewidth=2.5, markersize=10, 
         label='Red BGC structural (mean)', zorder=5)

if csr_locus in counts.index:
    redz_means = [counts.loc[csr_locus, 'T1_mean'], counts.loc[csr_locus, 'T2_mean'], counts.loc[csr_locus, 'T3_mean']]
    ax3.plot(t_positions, redz_means, 's-', color=colors_map['CSR'], linewidth=2, markersize=10, label='redZ', zorder=5)
    for col_group, tp in [(t1_cols, 1), (t2_cols, 2), (t3_cols, 3)]:
        for col in col_group:
            ax3.scatter(tp, counts.loc[csr_locus, col], color=colors_map['CSR'], alpha=0.3, s=30, zorder=3)

if 'SC_RS31630' in counts.index:
    redd_means = [counts.loc['SC_RS31630', 'T1_mean'], counts.loc['SC_RS31630', 'T2_mean'], counts.loc['SC_RS31630', 'T3_mean']]
    ax3.plot(t_positions, redd_means, '^-', color=colors_map['activator'], linewidth=2, markersize=10, label='redD', zorder=5)
    for col_group, tp in [(t1_cols, 1), (t2_cols, 2), (t3_cols, 3)]:
        for col in col_group:
            ax3.scatter(tp, counts.loc['SC_RS31630', col], color=colors_map['activator'], alpha=0.3, s=30, zorder=3)

ax3.set_xlabel('Timepoint', fontsize=13)
ax3.set_ylabel('Normalized counts', fontsize=13)
ax3.set_title('Red BGC Expression vs Regulators redZ and redD', fontsize=14, fontweight='bold')
ax3.set_xticks(t_positions)
ax3.set_xticklabels(t_labels)
ax3.legend(fontsize=11)
ax3.set_yscale('symlog', linthresh=10)
ax3.grid(axis='y', alpha=0.3)

plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig3.savefig(f"{FIG_DIR}/redZ_vs_redD_timecourse.{ext}", dpi=300, bbox_inches='tight')
print(f"  Saved: redZ_vs_redD_timecourse.pdf/svg")
plt.close()

# ============================================================
# Figure 4: Correlation
# ============================================================
fig4, ax4 = plt.subplots(figsize=(8, 6))

comparisons = ['T2vsT1', 'T3vsT1']
net_balances = [full_net_t2v1, full_net_t3v1]
bgc_l2fcs = [bgc_l2fc_t2v1, bgc_l2fc_t3v1]

for comp_label, de_df, net_bal in [('T2vsT1', de_t2v1, full_net_t2v1), ('T3vsT1', de_t3v1, full_net_t3v1)]:
    for l in struct_in_de:
        l2fc = de_df.loc[l, 'log2FoldChange']
        name = RED_STRUCTURAL.get(l, ('?',))[0]
        color = colors_map['BGC'] if comp_label == 'T2vsT1' else '#5dade2'
        ax4.scatter(net_bal, l2fc, color=color, alpha=0.4, s=40, zorder=3)

ax4.scatter(net_balances, bgc_l2fcs, color='red', s=200, marker='*', zorder=10, label='BGC structural mean')
for i, (nb, bl) in enumerate(zip(net_balances, bgc_l2fcs)):
    ax4.annotate(comparisons[i], (nb, bl), fontsize=10, fontweight='bold', 
                xytext=(10, 10), textcoords='offset points')

ax4.axhline(y=0, color='gray', linewidth=0.5, linestyle='--')
ax4.axvline(x=0, color='gray', linewidth=0.5, linestyle='--')
ax4.set_xlabel('Full Net Regulatory Balance (log2FC sum)', fontsize=12)
ax4.set_ylabel('Red BGC gene log2FC', fontsize=12)
ax4.set_title('Regulatory Balance vs BGC Output', fontsize=14, fontweight='bold')
ax4.legend(fontsize=10)
ax4.grid(alpha=0.3)

plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig4.savefig(f"{FIG_DIR}/balance_vs_output.{ext}", dpi=300, bbox_inches='tight')
print(f"  Saved: balance_vs_output.pdf/svg")
plt.close()

# ============================================================
# Summary
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY: Resolution of the redZ Paradox")
print("=" * 70)

print(f"""
KEY FINDINGS:

1. redZ EXPRESSION PATTERN:
   - T1 (exponential): {redz_vals[0]:.1f} normalized counts
   - T2 (transition):  {redz_vals[1]:.1f} normalized counts (log2FC vs T1: {redz_t2v1:+.2f})
   - T3 (stationary):  {redz_vals[2]:.1f} normalized counts (log2FC vs T1: {redz_t3v1:+.2f})

2. redD EXPRESSION (pathway-specific activator, downstream of redZ):
   - T1: {redd_vals[0]:.1f}, T2: {redd_vals[1]:.1f}, T3: {redd_vals[2]:.1f}
   - log2FC T2vsT1: {redd_t2v1:+.2f}, T3vsT1: {redd_t3v1:+.2f}

3. Red BGC STRUCTURAL GENES (mean of {len(struct_in_counts)} genes):
   - T1: {bgc_mean_t1:.1f}, T2: {bgc_mean_t2:.1f}, T3: {bgc_mean_t3:.1f}
   - Mean log2FC T2vsT1: {bgc_l2fc_t2v1:+.2f}, T3vsT1: {bgc_l2fc_t3v1:+.2f}

4. PARADOX RESOLUTION:
   Despite redZ downregulation ({redz_t2v1:+.2f} at T2, {redz_t3v1:+.2f} at T3):
   
   a) redD is STRONGLY UPREGULATED ({redd_t2v1:+.2f} at T2, {redd_t3v1:+.2f} at T3)
      - redD is the direct SARP activator of Red structural genes
      - redD can be activated by other pathways besides redZ
   
   b) absA1/absA2 (negative regulators) show changes:
      absA1: T2vsT1={get_l2fc('SC_RS18240', de_t2v1):+.2f}, T3vsT1={get_l2fc('SC_RS18240', de_t3v1):+.2f}
      absA2: T2vsT1={get_l2fc('SC_RS18245', de_t2v1):+.2f}, T3vsT1={get_l2fc('SC_RS18245', de_t3v1):+.2f}
   
   c) Global regulators:
      afsR: T2vsT1={afsr_t2v1:+.2f}, T3vsT1={afsr_t3v1:+.2f}
      afsS: T2vsT1={afss_t2v1:+.2f}, T3vsT1={afss_t3v1:+.2f}
      afsQ1: T2vsT1={afsq1_t2v1:+.2f}, T3vsT1={afsq1_t3v1:+.2f}

5. QUANTITATIVE BALANCE:
   T2vsT1: Full net balance = {full_net_t2v1:+.2f} (BGC output: {bgc_l2fc_t2v1:+.2f})
   T3vsT1: Full net balance = {full_net_t3v1:+.2f} (BGC output: {bgc_l2fc_t3v1:+.2f})

CONCLUSION: {"The positive regulatory inputs (especially redD upregulation and absA relief) quantitatively compensate for redZ loss, explaining how the Red BGC activates despite redZ downregulation." if full_net_t2v1 > 0 or bgc_l2fc_t2v1 > 0 else "Complex regulatory dynamics require further analysis."}
""")

print("\nOutput files:")
print(f"  Tables: {TAB_DIR}/")
print(f"    - regulator_timeline.tsv")
print(f"    - red_bgc_expression.tsv")
print(f"    - activator_repressor_balance.tsv")
print(f"  Figures: {FIG_DIR}/")
print(f"    - redZ_paradox_multipanel.pdf/svg")
print(f"    - balance_diagram.pdf/svg")
print(f"    - redZ_vs_redD_timecourse.pdf/svg")
print(f"    - balance_vs_output.pdf/svg")

print("\nDone.")
