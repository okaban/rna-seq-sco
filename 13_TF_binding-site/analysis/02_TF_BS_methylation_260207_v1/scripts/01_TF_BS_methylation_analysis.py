#!/usr/bin/env python3
"""
01: TF Binding Site Methylation Analysis
==========================================
Overlay methylation site positions (4mC/6mA) onto TF binding site coordinates,
assess spatial enrichment/depletion, temporal dynamics, and correlations with
target gene expression changes.

Data handling:
  - FIMO: absolute genomic coordinates (120 BS)
  - RegPrecise: relative positions → converted to absolute using target gene TSS
  - ZorroAranda Curated Strong: TF-target pairs → use target gene promoter region
  - ZorroAranda MEME: BS_sequence → exact match genome search
"""

import re
import gzip
from pathlib import Path

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# ── Paths ──
PROJECT = Path('/Users/okaban/bioinfo/rna-seq')
BS_DIR = PROJECT / '13_TF_binding-site' / 'analysis' / '01_master_TF_list_260206_v1'
EPIGENOME = PROJECT / '11_epigenome_integration'
OUT_DIR = PROJECT / '13_TF_binding-site' / 'analysis' / '02_TF_BS_methylation_260207_v1'
FIG_DIR = OUT_DIR / 'figures'
TBL_DIR = OUT_DIR / 'tables'
FIG_DIR.mkdir(parents=True, exist_ok=True)
TBL_DIR.mkdir(parents=True, exist_ok=True)

GENOME_SIZE = 8_667_507  # NC_003888.3
GENOME_FASTA = (PROJECT / 'reference' / 'GCF_000203835.1_ASM20383v1_genomic.fna')

# ── Style ──
plt.rcParams.update({
    'font.family': 'Arial', 'font.size': 10, 'axes.titlesize': 12,
    'axes.labelsize': 11, 'xtick.labelsize': 9, 'ytick.labelsize': 9,
    'legend.fontsize': 9, 'figure.dpi': 300, 'savefig.dpi': 300,
    'savefig.bbox': 'tight', 'axes.linewidth': 1.0,
    'axes.spines.top': False, 'axes.spines.right': False,
})
COL_4mC = '#E53935'
COL_6mA = '#1565C0'


def count_in_range(sorted_arr, lo, hi):
    """Binary search count of elements in [lo, hi]."""
    return int(np.searchsorted(sorted_arr, hi, 'right') - np.searchsorted(sorted_arr, lo, 'left'))


# ============================================================
# 1. Load and prepare data
# ============================================================
print("=" * 70)
print("TF BINDING SITE METHYLATION ANALYSIS")
print("=" * 70)

# -- Binding sites (raw) --
print("\n1. Loading binding site data...")
bs_raw = pd.read_csv(BS_DIR / 'master_TF_binding_sites_M145.tsv', sep='\t')

# -- Gene info for coordinate conversion --
gene_info = pd.read_csv(BS_DIR / 'intermediate' / 'M145_gene_basic_info.tsv', sep='\t')
# Build lookup: gene_id -> (start, end, strand)
# Also old_locus_tag (SCO_ID) -> gene_id mapping
gene_coord = {}
sco_to_geneid = {}
for _, g in gene_info.iterrows():
    gid = g.get('gene_id', '')
    sco = g.get('old_locus_tag', g.get('SCO_ID', ''))
    start = g.get('start', np.nan)
    end = g.get('end', np.nan)
    strand = g.get('strand', '+')
    if pd.notna(start) and pd.notna(end):
        gene_coord[gid] = (int(start), int(end), strand)
        if sco:
            sco_to_geneid[sco] = gid
            gene_coord[sco] = (int(start), int(end), strand)

# -- TSS table --
tss_df = pd.read_csv(EPIGENOME / 'analysis' / '18_tss_analyses' / 'comprehensive_tss_table.csv')
tss_map = dict(zip(tss_df['gene_id'], tss_df['tss']))
for _, t in tss_df.iterrows():
    olt = t.get('old_locus_tag', '')
    if pd.notna(olt) and olt:
        tss_map[olt] = t['tss']

# ============================================================
# 2. Resolve absolute coordinates for all BS sources
# ============================================================
print("\n2. Resolving absolute genomic coordinates...")

resolved_bs = []

for _, row in bs_raw.iterrows():
    src = row['BS_source']
    tf_name = row['TF_name']
    tf_gid = row['TF_gene_id']
    tg_sco = row.get('TG_SCO_ID', '')
    tg_gid = row.get('TG_gene_id', '')
    bs_seq = row.get('BS_sequence', '')

    rec = {
        'TF_name': tf_name, 'TF_gene_id': tf_gid,
        'TG_SCO_ID': tg_sco if pd.notna(tg_sco) else '',
        'TG_gene_id': tg_gid if pd.notna(tg_gid) else '',
        'BS_source': src,
        'BS_sequence': bs_seq if pd.notna(bs_seq) else '',
        'abs_start': np.nan, 'abs_end': np.nan,
    }

    if src == 'FIMO_curated_motif' and pd.notna(row['BS_start']):
        rec['abs_start'] = int(row['BS_start'])
        rec['abs_end'] = int(row['BS_end'])

    elif src == 'RegPrecise' and pd.notna(row['BS_start']):
        # BS_start = relative position (e.g. -147) from target gene
        rel_pos = int(row['BS_start'])
        tg_key = tg_gid if tg_gid else tg_sco
        tss_val = tss_map.get(tg_key) or tss_map.get(tg_sco)
        if tss_val and pd.notna(tss_val):
            tss_val = int(tss_val)
            coord = gene_coord.get(tg_key) or gene_coord.get(tg_sco)
            if coord:
                strand = coord[2]
                seq_len = len(bs_seq) if pd.notna(bs_seq) and bs_seq else 18
                if strand == '+':
                    abs_s = tss_val + rel_pos
                    abs_e = abs_s + seq_len - 1
                else:
                    abs_e = tss_val - rel_pos
                    abs_s = abs_e - seq_len + 1
                rec['abs_start'] = max(1, abs_s)
                rec['abs_end'] = abs_e

    elif src == 'ZorroAranda2022_Curated_Strong':
        # No coordinates, use target gene promoter region (-300 to +50 of TSS)
        tg_key = tg_gid if tg_gid else tg_sco
        tss_val = tss_map.get(tg_key) or tss_map.get(tg_sco)
        coord = gene_coord.get(tg_key) or gene_coord.get(tg_sco)
        if tss_val and pd.notna(tss_val) and coord:
            tss_val = int(tss_val)
            strand = coord[2]
            if strand == '+':
                rec['abs_start'] = max(1, tss_val - 300)
                rec['abs_end'] = tss_val + 50
            else:
                rec['abs_start'] = max(1, tss_val - 50)
                rec['abs_end'] = tss_val + 300

    # ZorroAranda2022_MEME: skip (too many, no coordinates, low quality)

    resolved_bs.append(rec)

bs_df = pd.DataFrame(resolved_bs)
bs_df = bs_df.dropna(subset=['abs_start', 'abs_end'])
bs_df['abs_start'] = bs_df['abs_start'].astype(int)
bs_df['abs_end'] = bs_df['abs_end'].astype(int)
bs_df['BS_center'] = ((bs_df['abs_start'] + bs_df['abs_end']) / 2).astype(int)
bs_df['BS_width'] = bs_df['abs_end'] - bs_df['abs_start'] + 1

# Tier assignment
tier1_sources = {'FIMO_curated_motif', 'RegPrecise', 'ZorroAranda2022_Curated_Strong'}
bs_df['tier'] = bs_df['BS_source'].apply(lambda x: 'Tier1' if x in tier1_sources else 'Tier2')

print(f"   Resolved BS with coordinates: {len(bs_df)}")
for src, cnt in bs_df['BS_source'].value_counts().items():
    print(f"     {src}: {cnt}")
print(f"   Tier 1: {(bs_df['tier']=='Tier1').sum()}")
print(f"   Tier 2: {(bs_df['tier']=='Tier2').sum()}")

# ============================================================
# 3. Load methylation data
# ============================================================
print("\n3. Loading methylation data...")
census_4mC = pd.read_csv(EPIGENOME / 'analysis' / '23_expanded_motif_search' / '4mC_final_census.csv')
census_6mA = pd.read_csv(EPIGENOME / 'analysis' / '23_expanded_motif_search' / '6mA_final_census.csv')
methyl_all = pd.concat([census_4mC, census_6mA], ignore_index=True)

methyl_unique = methyl_all.drop_duplicates(subset=['chrom', 'position', 'mod_type'])[
    ['chrom', 'position', 'strand', 'mod_type', 'final_motif']].copy()
print(f"   Unique methylation positions: {len(methyl_unique)} (4mC: {(methyl_unique.mod_type=='4mC').sum()}, 6mA: {(methyl_unique.mod_type=='6mA').sum()})")

# Temporal pivot
methyl_temporal = methyl_all.pivot_table(
    index=['chrom', 'position', 'mod_type', 'final_motif'],
    columns='timepoint', values='frequency', aggfunc='first'
).reset_index()
methyl_temporal.columns.name = None
for tp in ['T1', 'T2', 'T3']:
    if tp not in methyl_temporal.columns:
        methyl_temporal[tp] = np.nan

# Sorted position arrays
chr_methyl = methyl_unique[methyl_unique['chrom'] == 'NC_003888.3']
pos_4mC = np.sort(chr_methyl[chr_methyl['mod_type'] == '4mC']['position'].values)
pos_6mA = np.sort(chr_methyl[chr_methyl['mod_type'] == '6mA']['position'].values)
pos_all = np.sort(chr_methyl['position'].values)

# Expression data
print("\n4. Loading expression data...")
expr_df = pd.read_csv(EPIGENOME / 'analysis' / '01_integration' / 'integrated_methyl_expression.csv')
print(f"   Genes: {len(expr_df)}")

# TF master
tf_master = pd.read_csv(BS_DIR / 'master_TF_list_M145_with_BS.tsv', sep='\t')


# ============================================================
# ANALYSIS 1: Spatial Overlap
# ============================================================
print("\n" + "=" * 70)
print("ANALYSIS 1: SPATIAL OVERLAP")
print("=" * 70)

overlap_results = []
for tier_label in ['Tier1', 'All']:
    subset = bs_df if tier_label == 'All' else bs_df[bs_df['tier'] == 'Tier1']
    n_bs = len(subset)
    if n_bs == 0:
        continue

    for win_label, ext in [('direct', 0), ('±50bp', 50), ('±200bp', 200)]:
        n4, n6, na, bs_hit = 0, 0, 0, 0
        total_w = 0
        for _, r in subset.iterrows():
            s, e = r['abs_start'] - ext, r['abs_end'] + ext
            total_w += (e - s + 1)
            c4 = count_in_range(pos_4mC, s, e)
            c6 = count_in_range(pos_6mA, s, e)
            ca = count_in_range(pos_all, s, e)
            n4 += c4; n6 += c6; na += ca
            if ca > 0:
                bs_hit += 1

        e4 = len(pos_4mC) * total_w / GENOME_SIZE
        e6 = len(pos_6mA) * total_w / GENOME_SIZE
        ea = len(pos_all) * total_w / GENOME_SIZE
        f4 = n4 / e4 if e4 > 0 else 0
        f6 = n6 / e6 if e6 > 0 else 0
        fa = na / ea if ea > 0 else 0
        p4 = stats.poisson.sf(n4 - 1, e4) if n4 > e4 else stats.poisson.cdf(n4, e4)
        p6 = stats.poisson.sf(n6 - 1, e6) if n6 > e6 else stats.poisson.cdf(n6, e6)
        pa = stats.poisson.sf(na - 1, ea) if na > ea else stats.poisson.cdf(na, ea)

        overlap_results.append({
            'tier': tier_label, 'window': win_label, 'n_BS': n_bs,
            'total_width_bp': total_w,
            'BS_with_methyl': bs_hit,
            'BS_with_methyl_pct': round(bs_hit / n_bs * 100, 1),
            'n_4mC_obs': n4, 'n_4mC_exp': round(e4, 1), 'fold_4mC': round(f4, 3), 'p_4mC': p4,
            'n_6mA_obs': n6, 'n_6mA_exp': round(e6, 1), 'fold_6mA': round(f6, 3), 'p_6mA': p6,
            'n_any_obs': na, 'n_any_exp': round(ea, 1), 'fold_any': round(fa, 3), 'p_any': pa,
        })
        print(f"  {tier_label} {win_label}: {n_bs} BS, {bs_hit} ({bs_hit/n_bs*100:.1f}%) with methyl | "
              f"4mC {n4}(exp {e4:.1f}, fold={f4:.2f}, p={p4:.3g}) | "
              f"6mA {n6}(exp {e6:.1f}, fold={f6:.2f}, p={p6:.3g})")

pd.DataFrame(overlap_results).to_csv(TBL_DIR / 'T1_BS_methylation_overlap_summary.tsv', sep='\t', index=False)
print(f"\n   Saved: T1_BS_methylation_overlap_summary.tsv")


# ============================================================
# ANALYSIS 2: Metagene Profile
# ============================================================
print("\n" + "=" * 70)
print("ANALYSIS 2: METAGENE PROFILE")
print("=" * 70)

WINDOW = 2000
BIN = 100
bins = np.arange(-WINDOW, WINDOW + BIN, BIN)
bin_c = (bins[:-1] + bins[1:]) / 2

bs_tier1 = bs_df[bs_df['tier'] == 'Tier1']
centers = bs_tier1['BS_center'].values
n_c = len(centers)
print(f"   Tier 1 BS centers: {n_c}")

h4 = np.zeros(len(bin_c))
h6 = np.zeros(len(bin_c))
for c in centers:
    for pos_arr, hist in [(pos_4mC, h4), (pos_6mA, h6)]:
        lo_idx = np.searchsorted(pos_arr, c - WINDOW, 'left')
        hi_idx = np.searchsorted(pos_arr, c + WINDOW, 'right')
        for p in pos_arr[lo_idx:hi_idx]:
            d = p - c
            idx = int((d + WINDOW) // BIN)
            if 0 <= idx < len(hist):
                hist[idx] += 1

d4 = h4 / n_c / (BIN / 1000)
d6 = h6 / n_c / (BIN / 1000)
exp4 = len(pos_4mC) / (GENOME_SIZE / 1000)
exp6 = len(pos_6mA) / (GENOME_SIZE / 1000)

fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
axes[0].bar(bin_c, d4, width=BIN * 0.9, color=COL_4mC, alpha=0.7, label='4mC')
axes[0].axhline(exp4, color='gray', ls='--', lw=1, label=f'Genome avg ({exp4:.3f}/kb)')
axes[0].axvline(0, color='black', ls='-', lw=0.8, alpha=0.5)
axes[0].set_ylabel('4mC sites / BS / kb')
axes[0].set_title('A. 4mC Density Around TF Binding Site Centers (Tier 1)', fontweight='bold', loc='left')
axes[0].legend(fontsize=8)

axes[1].bar(bin_c, d6, width=BIN * 0.9, color=COL_6mA, alpha=0.7, label='6mA')
axes[1].axhline(exp6, color='gray', ls='--', lw=1, label=f'Genome avg ({exp6:.3f}/kb)')
axes[1].axvline(0, color='black', ls='-', lw=0.8, alpha=0.5)
axes[1].set_ylabel('6mA sites / BS / kb')
axes[1].set_xlabel('Distance from BS center (bp)')
axes[1].set_title('B. 6mA Density Around TF Binding Site Centers (Tier 1)', fontweight='bold', loc='left')
axes[1].legend(fontsize=8)

plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(FIG_DIR / f'F1_metagene_methylation_around_BS.{ext}')
plt.close()
print(f"   Saved: F1_metagene_methylation_around_BS")


# ============================================================
# ANALYSIS 3: Temporal Dynamics at BS
# ============================================================
print("\n" + "=" * 70)
print("ANALYSIS 3: TEMPORAL DYNAMICS")
print("=" * 70)

overlap_sites = []
for _, bsr in bs_tier1.iterrows():
    s, e = bsr['abs_start'] - 50, bsr['abs_end'] + 50
    mask = (methyl_temporal['chrom'] == 'NC_003888.3') & \
           (methyl_temporal['position'] >= s) & (methyl_temporal['position'] <= e)
    for _, m in methyl_temporal[mask].iterrows():
        direct = bsr['abs_start'] <= m['position'] <= bsr['abs_end']
        overlap_sites.append({
            'TF_name': bsr['TF_name'], 'TF_gene_id': bsr['TF_gene_id'],
            'TG_gene_id': bsr['TG_gene_id'], 'BS_source': bsr['BS_source'],
            'abs_start': bsr['abs_start'], 'abs_end': bsr['abs_end'],
            'methyl_pos': m['position'], 'mod_type': m['mod_type'],
            'final_motif': m['final_motif'], 'direct_overlap': direct,
            'T1_freq': m.get('T1', np.nan), 'T2_freq': m.get('T2', np.nan),
            'T3_freq': m.get('T3', np.nan),
        })

ov_df = pd.DataFrame(overlap_sites)
print(f"   Methylation sites at Tier 1 BS (±50bp): {len(ov_df)}")

if len(ov_df) > 0:
    print(f"     Direct overlap: {ov_df['direct_overlap'].sum()}")
    print(f"     4mC: {(ov_df.mod_type=='4mC').sum()}, 6mA: {(ov_df.mod_type=='6mA').sum()}")

    def classify(r):
        t1, t2, t3 = r['T1_freq'], r['T2_freq'], r['T3_freq']
        h1 = pd.notna(t1) and t1 > 0
        h2 = pd.notna(t2) and t2 > 0
        h3 = pd.notna(t3) and t3 > 0
        if h1 and not h2 and not h3: return 'Lost_T2T3'
        if h1 and h2 and not h3: return 'Lost_T3'
        if not h1 and (h2 or h3): return 'Gained'
        if h1 and h2 and h3:
            v = [x for x in [t1, t2, t3] if pd.notna(x)]
            return 'Dynamic' if max(v) - min(v) > 10 else 'Stable'
        return 'Dynamic' if (h1 or h2 or h3) else 'Other'

    ov_df['dynamics'] = ov_df.apply(classify, axis=1)
    for cat, cnt in ov_df['dynamics'].value_counts().items():
        print(f"     {cat}: {cnt}")

    ov_df.to_csv(TBL_DIR / 'T2_BS_methylation_temporal_dynamics.tsv', sep='\t', index=False)

    # Heatmap
    plot_df = ov_df.copy()
    plot_df['label'] = plot_df.apply(
        lambda r: f"{r['TF_name']}|{r['mod_type']}|{int(r['methyl_pos'])}", axis=1)
    heat = plot_df[['label', 'T1_freq', 'T2_freq', 'T3_freq']].set_index('label')
    heat.columns = ['T1', 'T2', 'T3']
    heat = heat.dropna(how='all')
    # Remove duplicate labels, keep first
    heat = heat[~heat.index.duplicated(keep='first')]

    if len(heat) > 0:
        fig_h, ax_h = plt.subplots(figsize=(5, max(3, len(heat) * 0.35)))
        sns.heatmap(heat, cmap='YlOrRd', annot=True, fmt='.0f', linewidths=0.5,
                    cbar_kws={'label': 'Methylation freq (%)'}, ax=ax_h, vmin=0, vmax=100)
        ax_h.set_title('Methylation at TF Binding Sites (Tier 1, ±50bp)', fontweight='bold', fontsize=10)
        ax_h.set_xlabel('Timepoint'); ax_h.set_ylabel('')
        plt.tight_layout()
        for ext in ['pdf', 'svg']:
            fig_h.savefig(FIG_DIR / f'F2_BS_methylation_temporal_heatmap.{ext}')
        plt.close()
        print(f"   Saved: F2_BS_methylation_temporal_heatmap")
else:
    ov_df = pd.DataFrame()
    print("   No overlap found")


# ============================================================
# ANALYSIS 4: BS Methylation vs Target Expression
# ============================================================
print("\n" + "=" * 70)
print("ANALYSIS 4: BS METHYLATION vs TARGET EXPRESSION")
print("=" * 70)

if len(ov_df) > 0:
    from matplotlib.lines import Line2D

    has_tg = ov_df[ov_df['TG_gene_id'].notna() & (ov_df['TG_gene_id'] != '')].copy()
    print(f"   Overlap sites with target gene: {len(has_tg)}")

    expr_cols = ['gene_id', 'log2FC_T2_vs_T1', 'padj_T2_vs_T1',
                 'log2FC_T3_vs_T1', 'padj_T3_vs_T1',
                 'log2FC_T3_vs_T2', 'padj_T3_vs_T2']

    if len(has_tg) > 0:
        merged = has_tg.merge(expr_df[expr_cols], left_on='TG_gene_id', right_on='gene_id', how='left')
        merged['delta_T2vsT1'] = merged['T2_freq'] - merged['T1_freq']
        merged['delta_T3vsT1'] = merged['T3_freq'] - merged['T1_freq']
        merged.to_csv(TBL_DIR / 'T3_BS_methylation_expression_pairs.tsv', sep='\t', index=False)
        print(f"   Saved: T3_BS_methylation_expression_pairs.tsv")

    # --- Approach: Lost/Gained/No-methylation BS vs target gene expression ---
    # Build per-BS classification: Lost, Gained, No methylation
    # For each BS with a target gene, classify methylation state and get target expression
    bs_with_tg = bs_df[bs_df['TG_gene_id'].notna() & (bs_df['TG_gene_id'] != '')].copy()
    print(f"   BS with target gene info: {len(bs_with_tg)}")

    bs_expr_list = []
    for _, bsr in bs_with_tg.iterrows():
        s, e = bsr['abs_start'] - 50, bsr['abs_end'] + 50
        mask = (methyl_temporal['chrom'] == 'NC_003888.3') & \
               (methyl_temporal['position'] >= s) & (methyl_temporal['position'] <= e)
        matched = methyl_temporal[mask]

        # Target expression
        tg = bsr['TG_gene_id']
        tg_expr = expr_df[expr_df['gene_id'] == tg]
        if len(tg_expr) == 0:
            continue
        lfc_t2 = tg_expr.iloc[0].get('log2FC_T2_vs_T1', np.nan)
        lfc_t3 = tg_expr.iloc[0].get('log2FC_T3_vs_T1', np.nan)

        if len(matched) == 0:
            bs_expr_list.append({
                'methyl_class': 'No methylation', 'log2FC_T2_vs_T1': lfc_t2,
                'log2FC_T3_vs_T1': lfc_t3, 'TF_name': bsr['TF_name'],
                'TG_gene_id': tg, 'n_methyl': 0,
            })
        else:
            # Classify: at least one Lost? at least one Gained?
            cats = []
            for _, m in matched.iterrows():
                t1h = pd.notna(m.get('T1')) and m.get('T1', 0) > 0
                t2h = pd.notna(m.get('T2')) and m.get('T2', 0) > 0
                t3h = pd.notna(m.get('T3')) and m.get('T3', 0) > 0
                if t1h and not t2h and not t3h:
                    cats.append('Lost')
                elif not t1h and (t2h or t3h):
                    cats.append('Gained')
                else:
                    cats.append('Other')
            # Dominant pattern
            if 'Lost' in cats and 'Gained' not in cats:
                mc = 'Lost at BS'
            elif 'Gained' in cats and 'Lost' not in cats:
                mc = 'Gained at BS'
            else:
                mc = 'Mixed'
            bs_expr_list.append({
                'methyl_class': mc, 'log2FC_T2_vs_T1': lfc_t2,
                'log2FC_T3_vs_T1': lfc_t3, 'TF_name': bsr['TF_name'],
                'TG_gene_id': tg, 'n_methyl': len(matched),
            })

    bse_df = pd.DataFrame(bs_expr_list)
    print(f"\n   BS-target pairs by methylation class:")
    for c, n in bse_df['methyl_class'].value_counts().items():
        print(f"     {c}: {n}")

    # Figure: boxplot of target log2FC by methylation class
    fig_s, axes_s = plt.subplots(1, 2, figsize=(13, 5.5))
    class_order = ['No methylation', 'Lost at BS', 'Gained at BS', 'Mixed']
    class_colors = {'No methylation': '#BDBDBD', 'Lost at BS': '#E53935',
                    'Gained at BS': '#1565C0', 'Mixed': '#7E57C2'}

    for ax, yv, lab in [
        (axes_s[0], 'log2FC_T2_vs_T1', 'T2 vs T1'),
        (axes_s[1], 'log2FC_T3_vs_T1', 'T3 vs T1'),
    ]:
        present = [c for c in class_order if c in bse_df['methyl_class'].values]
        d = bse_df[bse_df['methyl_class'].isin(present)].dropna(subset=[yv])
        if len(d) < 3:
            ax.text(0.5, 0.5, 'Insufficient data', ha='center', va='center',
                    transform=ax.transAxes, fontsize=12, color='gray')
            ax.set_title(lab, fontweight='bold', loc='left')
            continue

        bp = ax.boxplot(
            [d[d['methyl_class'] == c][yv].values for c in present],
            labels=[f"{c}\n(n={len(d[d['methyl_class']==c])})" for c in present],
            patch_artist=True, widths=0.6, showfliers=True,
            flierprops={'markersize': 3, 'alpha': 0.4},
        )
        for patch, c in zip(bp['boxes'], present):
            patch.set_facecolor(class_colors.get(c, 'gray'))
            patch.set_alpha(0.7)

        ax.axhline(0, color='gray', ls='--', lw=0.5, alpha=0.5)
        ax.set_ylabel('Target gene log2FC')
        ax.set_title(f'{lab}', fontweight='bold', loc='left', fontsize=11)

        # Stats: compare Lost/Gained vs No-methyl if n≥5
        nometh = d[d['methyl_class'] == 'No methylation'][yv].dropna()
        for mc in ['Lost at BS', 'Gained at BS']:
            grp = d[d['methyl_class'] == mc][yv].dropna()
            if len(grp) >= 5 and len(nometh) >= 5:
                u, p = stats.mannwhitneyu(grp, nometh, alternative='two-sided')
                med_diff = grp.median() - nometh.median()
                print(f"   {lab} | {mc} vs No-methyl: U={u}, p={p:.3g}, Δmedian={med_diff:.3f}")

    plt.suptitle('Target Gene Expression by BS Methylation Status', fontweight='bold', fontsize=12, y=1.02)
    plt.tight_layout()
    for ext in ['pdf', 'svg']:
        fig_s.savefig(FIG_DIR / f'F3_BS_methylation_vs_expression.{ext}')
    plt.close()
    print(f"   Saved: F3_BS_methylation_vs_expression")
else:
    print("   Skipped")


# ============================================================
# ANALYSIS 5: Motif & TF Family
# ============================================================
print("\n" + "=" * 70)
print("ANALYSIS 5: MOTIF & TF FAMILY")
print("=" * 70)

# Per-BS methylation flag (all BS, ±50bp)
bs_flags = []
for _, r in bs_df.iterrows():
    s, e = r['abs_start'] - 50, r['abs_end'] + 50
    bs_flags.append({
        'TF_name': r['TF_name'], 'TF_gene_id': r['TF_gene_id'],
        'BS_source': r['BS_source'], 'tier': r['tier'],
        'has_methyl': count_in_range(pos_all, s, e) > 0,
        'n_4mC': count_in_range(pos_4mC, s, e),
        'n_6mA': count_in_range(pos_6mA, s, e),
    })
bsf = pd.DataFrame(bs_flags)

# TF family
tf_fam = dict(zip(tf_master['gene_id'], tf_master['TF_family']))
bsf['TF_family'] = bsf['TF_gene_id'].map(tf_fam).fillna('Unknown')

# Motif composition at BS
if len(ov_df) > 0:
    motif_bs = ov_df['final_motif'].value_counts()
    motif_genome = methyl_unique['final_motif'].value_counts()
    mc = pd.DataFrame({'at_BS': motif_bs, 'genome': motif_genome}).fillna(0).astype(int)
    mc['at_BS_pct'] = (mc['at_BS'] / mc['at_BS'].sum() * 100).round(1)
    mc['genome_pct'] = (mc['genome'] / mc['genome'].sum() * 100).round(1)
    print("\n   Motif composition at BS vs genome:")
    print(mc.to_string())

# TF family methylation rate
fam = bsf.groupby('TF_family').agg(
    n_BS=('has_methyl', 'count'), n_hit=('has_methyl', 'sum')
).reset_index()
fam['rate'] = (fam['n_hit'] / fam['n_BS'] * 100).round(1)
fam = fam.sort_values('n_BS', ascending=False)
print("\n   TF family methylation rate (all BS, ±50bp):")
print(fam.to_string(index=False))

fam.to_csv(TBL_DIR / 'T4_motif_TF_family_summary.tsv', sep='\t', index=False)

# Figure
plot_f = fam[fam['n_BS'] >= 5].sort_values('rate', ascending=True)
if len(plot_f) > 0:
    fig_f, ax_f = plt.subplots(figsize=(7, max(3, len(plot_f) * 0.4)))
    ax_f.barh(range(len(plot_f)), plot_f['rate'], color='#546E7A', edgecolor='black', linewidth=0.5)
    ax_f.set_yticks(range(len(plot_f)))
    ax_f.set_yticklabels(plot_f['TF_family'])
    ax_f.set_xlabel('BS with methylation (±50bp) (%)')
    ax_f.set_title('TF Family: BS Methylation Rate (≥5 BS)', fontweight='bold', fontsize=10)
    for i, (_, rr) in enumerate(plot_f.iterrows()):
        ax_f.text(rr['rate'] + 0.5, i, f"{int(rr['n_hit'])}/{int(rr['n_BS'])}", va='center', fontsize=8)
    avg = bsf['has_methyl'].mean() * 100
    ax_f.axvline(avg, color='red', ls='--', lw=1, label=f'Average ({avg:.1f}%)')
    ax_f.legend(fontsize=8)
    plt.tight_layout()
    for ext in ['pdf', 'svg']:
        fig_f.savefig(FIG_DIR / f'F4_TF_family_BS_methylation_rate.{ext}')
    plt.close()
    print(f"   Saved: F4_TF_family_BS_methylation_rate")


# ============================================================
# ANALYSIS 6: Key TF Highlights
# ============================================================
print("\n" + "=" * 70)
print("ANALYSIS 6: KEY TF HIGHLIGHTS")
print("=" * 70)

key_tfs = ['redZ', 'actII-ORF4', 'bldD', 'afsR', 'sigR', 'sigE', 'sigB',
           'phoP', 'dasR', 'glnR', 'Rex', 'DraR', 'AbrC3']

key_results = []
for tf_name in key_tfs:
    tf_bs = bs_df[bs_df['TF_name'].str.lower() == tf_name.lower()]
    if len(tf_bs) == 0:
        print(f"   {tf_name}: no BS data")
        continue

    n_total = len(tf_bs)
    n_tier1 = (tf_bs['tier'] == 'Tier1').sum()
    n_methyl_bs = 0
    methyl_details = []

    for _, bsr in tf_bs.iterrows():
        s, e = bsr['abs_start'] - 50, bsr['abs_end'] + 50
        mask = (methyl_temporal['chrom'] == 'NC_003888.3') & \
               (methyl_temporal['position'] >= s) & (methyl_temporal['position'] <= e)
        matched = methyl_temporal[mask]
        if len(matched) > 0:
            n_methyl_bs += 1
            for _, m in matched.iterrows():
                methyl_details.append({
                    'TF_name': tf_name, 'BS_start': bsr['abs_start'], 'BS_end': bsr['abs_end'],
                    'BS_source': bsr['BS_source'], 'TG_gene_id': bsr.get('TG_gene_id', ''),
                    'methyl_pos': m['position'], 'mod_type': m['mod_type'],
                    'motif': m['final_motif'],
                    'T1': m.get('T1', np.nan), 'T2': m.get('T2', np.nan), 'T3': m.get('T3', np.nan),
                })

    tf_gids = tf_bs['TF_gene_id'].unique()
    for gid in tf_gids:
        tf_e = expr_df[expr_df['gene_id'] == gid]
        if len(tf_e) > 0:
            re_ = tf_e.iloc[0]
            key_results.append({
                'TF_name': tf_name, 'TF_gene_id': gid,
                'n_BS': n_total, 'n_BS_tier1': n_tier1,
                'n_BS_with_methyl': n_methyl_bs, 'n_methyl_sites': len(methyl_details),
                'TF_log2FC_T2vsT1': re_.get('log2FC_T2_vs_T1', np.nan),
                'TF_padj_T2vsT1': re_.get('padj_T2_vs_T1', np.nan),
                'TF_log2FC_T3vsT1': re_.get('log2FC_T3_vs_T1', np.nan),
                'TF_padj_T3vsT1': re_.get('padj_T3_vs_T1', np.nan),
            })

    if methyl_details:
        print(f"\n   {tf_name} ({n_total} BS, {n_tier1} Tier1): {n_methyl_bs} BS with methyl, {len(methyl_details)} sites")
        for md in methyl_details[:5]:
            ts = ' '.join([f"T{i}={md[f'T{i}']:.0f}%" if pd.notna(md[f'T{i}']) else f"T{i}=N/A" for i in [1,2,3]])
            print(f"     pos={md['methyl_pos']}, {md['mod_type']}, {md['motif']}, {ts}")
    else:
        print(f"   {tf_name}: {n_total} BS, no methyl overlap")

kr = pd.DataFrame(key_results)
if len(kr) > 0:
    kr.to_csv(TBL_DIR / 'T5_key_TF_BS_methylation_detail.tsv', sep='\t', index=False)
    print(f"\n   Saved: T5_key_TF_BS_methylation_detail.tsv")

# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"  Tables: {len(list(TBL_DIR.glob('*.tsv')))} TSV files")
print(f"  Figures: {len(list(FIG_DIR.glob('*.pdf')))} PDF + SVG pairs")
print("Done.")
