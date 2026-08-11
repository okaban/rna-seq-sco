#!/usr/bin/env python3
"""
75_tss_downstream_protection: TSS上流・下流GCCGGC保護ゾーン解析
=================================================================

目的:
    GCCGGC (4mC) メチル化フリーゾーンがTSS上流だけでなく
    下流 (+1〜+293bp) にも存在するかを定量し、
    Shielded/Exposed遺伝子間で比較する。

パネル:
    A: TSS中心化メタプロファイル (-500〜+500bp, 25bpビン)
       - Shielded遺伝子 vs Exposed遺伝子
       - 95% CI バンド付き
    B: 上流密度 vs 下流密度の散布図 (293bp範囲)
       - 遺伝子ごとの点、色分け

データ:
    - all_genes_features_unified_n57.tsv: 1,055 regulatory genes, TSS, is_exposed
    - 4mC_final_census.csv: GCCGGC (TGGCCGGC/GGCCGG) 4mCサイト T1
    - comprehensive_tss_table.csv: 全遺伝子TSS (Jeong2016 + GFF)

実施日: 2026-05-10
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
from pathlib import Path

# ── パス設定 ────────────────────────────────────────────────────────────────
BASE      = Path('/Users/okaban/bioinfo/rna-seq')
EPIGENOME = BASE / '11_epigenome_integration' / 'analysis'
ANALYSIS  = EPIGENOME / '75_tss_downstream_protection'
FIG_DIR   = ANALYSIS / 'figures'
TAB_DIR   = ANALYSIS / 'tables'
FIG_SUP   = BASE / '15_paper_figures' / 'figures' / 'supp'

FIG_DIR.mkdir(parents=True, exist_ok=True)
TAB_DIR.mkdir(parents=True, exist_ok=True)
FIG_SUP.mkdir(parents=True, exist_ok=True)

# ── パラメータ ────────────────────────────────────────────────────────────
BIN_SIZE   = 25       # bp
HALF_WIN   = 500      # TSS ± 500bp
DOWN_BOUND = 293      # downstream boundary (bp) — 既存の293bp閾値に対応
SIGMA_BINS = 1.5      # Gaussian smoothing (bins)

# カラー (shared_utils に合わせる)
COL_SHIELDED = '#B0BEC5'   # gray
COL_EXPOSED  = '#7E57C2'   # purple
COL_ALPHA    = 0.20

STYLE = {
    'font.family': 'Arial',
    'font.size': 10,
    'axes.titlesize': 11,
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'axes.linewidth': 1.0,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'pdf.fonttype': 42,
    'svg.fonttype': 'none',
}
plt.rcParams.update(STYLE)

# ── 1. データ読み込み ────────────────────────────────────────────────────────
print("=" * 65)
print("75_tss_downstream_protection: TSS上流・下流GCCGGC保護ゾーン解析")
print("=" * 65)

# (a) 規制遺伝子フィーチャー表 (is_exposed含む)
feat_path = EPIGENOME / '52_shielded_exposed_boundary' / 'tables' / 'all_genes_features_unified_n57.tsv'
df_feat = pd.read_csv(feat_path, sep='\t')
print(f"\n[1a] all_genes_features_unified: {len(df_feat)} genes")
print(f"     Exposed: {df_feat['is_exposed'].sum()}, Shielded: {(df_feat['is_exposed']==0).sum()}")

# (b) GCCGGC 4mC サイト (T1)
census_path = EPIGENOME / '23_expanded_motif_search' / '4mC_final_census.csv'
df_census = pd.read_csv(census_path)
gccggc_t1 = (df_census[
    (df_census['final_motif'].isin(['TGGCCGGC', 'GGCCGG'])) &
    (df_census['timepoint'] == 'T1')
].drop_duplicates(subset=['position'])).copy()
gccggc_pos = gccggc_t1['position'].values
print(f"\n[1b] GCCGGC 4mC sites (T1, unique positions): {len(gccggc_pos)}")

# (c) 補助: 全遺伝子TSS (Jeong2016 優先)
tss_path = EPIGENOME / '18_tss_analyses' / 'comprehensive_tss_table.csv'
df_tss_all = pd.read_csv(tss_path)
# Jeong2016 優先でマージ
jeong = df_tss_all[df_tss_all['tss_source'] == 'Jeong2016_dRNA-seq'][['gene_id','tss','strand']].rename(columns={'gene_id':'locus_tag'})
print(f"[1c] Jeong2016 experimental TSS: {len(jeong)} genes")

# df_feat に Jeong TSS をマージ（元のTSSで補完）
df_genes = df_feat.copy()
df_genes = df_genes.merge(jeong.rename(columns={'tss':'jeong_tss'}), on='locus_tag', how='left', suffixes=('','_j'))
# Jeong優先: 欠損なら元のTSSを使用
df_genes['tss_use'] = df_genes['jeong_tss'].combine_first(df_genes['tss'])
n_jeong = df_genes['jeong_tss'].notna().sum()
print(f"     Regulatory genes with Jeong TSS: {n_jeong} / {len(df_genes)}")
print(f"     Using annotation TSS for remaining: {len(df_genes) - n_jeong}")

# strand列が重複した場合の処理
if 'strand_j' in df_genes.columns:
    df_genes = df_genes.drop(columns=['strand_j'])

# ── 2. 各遺伝子のGCCGGC相対位置を計算 ──────────────────────────────────────
print("\n[2] Computing per-gene GCCGGC relative positions...")

# ビン定義
bin_edges   = np.arange(-HALF_WIN, HALF_WIN + BIN_SIZE, BIN_SIZE)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
n_bins      = len(bin_centers)

# 各遺伝子のビン別サイト数を蓄積
def compute_gene_profiles(df_subset, gccggc_positions, window=HALF_WIN, bin_size=BIN_SIZE):
    """Returns matrix (n_genes, n_bins) of GCCGGC site counts."""
    edges = np.arange(-window, window + bin_size, bin_size)
    profiles = []
    for _, row in df_subset.iterrows():
        tss    = row['tss_use']
        strand = row['strand']
        if pd.isna(tss):
            profiles.append(np.full(len(edges)-1, np.nan))
            continue
        # 近傍のサイトのみ処理 (速度向上)
        mask = (gccggc_positions >= tss - window - 50) & (gccggc_positions <= tss + window + 50)
        nearby = gccggc_positions[mask]
        if strand == '+':
            rel = nearby - tss
        else:
            rel = tss - nearby
        counts, _ = np.histogram(rel, bins=edges)
        profiles.append(counts.astype(float))
    return np.array(profiles)

# Shielded / Exposed 分割
mask_sh = df_genes['is_exposed'] == 0
mask_ex = df_genes['is_exposed'] == 1

df_sh = df_genes[mask_sh].reset_index(drop=True)
df_ex = df_genes[mask_ex].reset_index(drop=True)

print(f"  Shielded: {len(df_sh)} genes, Exposed: {len(df_ex)} genes")

prof_sh = compute_gene_profiles(df_sh, gccggc_pos)
prof_ex = compute_gene_profiles(df_ex, gccggc_pos)

print(f"  Profile matrices: Shielded={prof_sh.shape}, Exposed={prof_ex.shape}")

# ── 3. 密度計算 (sites/kb/gene) ──────────────────────────────────────────────
bin_kb = BIN_SIZE / 1000.0

def profile_stats(prof_mat):
    """Mean, CI95 (bootstrap-free: SEM×1.96) density per kb."""
    n_valid  = (~np.isnan(prof_mat)).sum(axis=0)
    mean_d   = np.nanmean(prof_mat, axis=0) / bin_kb
    sem_d    = np.nanstd(prof_mat, axis=0, ddof=1) / np.sqrt(np.maximum(n_valid, 1)) / bin_kb
    ci_low   = mean_d - 1.96 * sem_d
    ci_high  = mean_d + 1.96 * sem_d
    return mean_d, ci_low, ci_high

def gaussian_smooth(arr, sigma=SIGMA_BINS):
    from scipy.ndimage import gaussian_filter1d
    return gaussian_filter1d(arr, sigma=sigma)

mean_sh, ci_lo_sh, ci_hi_sh = profile_stats(prof_sh)
mean_ex, ci_lo_ex, ci_hi_ex = profile_stats(prof_ex)

# スムージング
sm_sh    = gaussian_smooth(mean_sh)
sm_ex    = gaussian_smooth(mean_ex)
sm_lo_sh = gaussian_smooth(ci_lo_sh)
sm_hi_sh = gaussian_smooth(ci_hi_sh)
sm_lo_ex = gaussian_smooth(ci_lo_ex)
sm_hi_ex = gaussian_smooth(ci_hi_ex)

print(f"\n[3] Density stats:")
# 上流・下流の平均密度 (全遺伝子プール)
up_mask   = bin_centers < 0
down_mask = (bin_centers >= 0) & (bin_centers <= DOWN_BOUND)
print(f"  Shielded — upstream mean: {mean_sh[up_mask].mean():.4f}  downstream mean: {mean_sh[down_mask].mean():.4f}")
print(f"  Exposed  — upstream mean: {mean_ex[up_mask].mean():.4f}  downstream mean: {mean_ex[down_mask].mean():.4f}")

# ── 4. 遺伝子ごとの上流・下流密度 (Panel B用) ──────────────────────────────
def per_gene_density(prof_mat, bin_centers, lo, hi):
    """Per-gene density in region [lo, hi) bp from TSS."""
    mask = (bin_centers >= lo) & (bin_centers < hi)
    region_bp = (hi - lo)  # bp
    counts = np.nansum(prof_mat[:, mask], axis=1)
    return counts / (region_bp / 1000.0)   # sites/kb

up_sh  = per_gene_density(prof_sh, bin_centers, -HALF_WIN, 0)
dn_sh  = per_gene_density(prof_sh, bin_centers, 0, DOWN_BOUND)
up_ex  = per_gene_density(prof_ex, bin_centers, -HALF_WIN, 0)
dn_ex  = per_gene_density(prof_ex, bin_centers, 0, DOWN_BOUND)

# ── 5. 統計検定 ──────────────────────────────────────────────────────────────
print("\n[4] Statistical tests (Shielded vs Exposed):")

def mw_test(a, b, label):
    stat, p = stats.mannwhitneyu(a, b, alternative='two-sided')
    print(f"  {label}: U={stat:.0f}, p={p:.3e} (n_sh={len(a)}, n_ex={len(b)})")
    return p

p_up = mw_test(up_sh, up_ex, "Upstream density  [-500, 0)")
p_dn = mw_test(dn_sh, dn_ex, f"Downstream density [0, +{DOWN_BOUND})")

# 上流 vs 下流 比率 (Shielded)
ratio_sh = (up_sh + 1e-6) / (dn_sh + 1e-6)
ratio_ex = (up_ex + 1e-6) / (dn_ex + 1e-6)
p_ratio  = mw_test(ratio_sh, ratio_ex, "Upstream/Downstream ratio")

# テーブル保存
df_sh_stats = pd.DataFrame({
    'locus_tag': df_sh['locus_tag'],
    'category': 'Shielded',
    'upstream_density': up_sh,
    'downstream_density': dn_sh,
})
df_ex_stats = pd.DataFrame({
    'locus_tag': df_ex['locus_tag'],
    'category': 'Exposed',
    'upstream_density': up_ex,
    'downstream_density': dn_ex,
})
df_stats = pd.concat([df_sh_stats, df_ex_stats], ignore_index=True)
df_stats.to_csv(TAB_DIR / 'per_gene_upstream_downstream_density.tsv', sep='\t', index=False)
print(f"\n  Saved: per_gene_upstream_downstream_density.tsv")

# ビン別プロファイルテーブル
df_profile = pd.DataFrame({
    'bin_center_bp': bin_centers,
    'shielded_mean': mean_sh,
    'shielded_ci_lo': ci_lo_sh,
    'shielded_ci_hi': ci_hi_sh,
    'exposed_mean': mean_ex,
    'exposed_ci_lo': ci_lo_ex,
    'exposed_ci_hi': ci_hi_ex,
})
df_profile.to_csv(TAB_DIR / 'tss_metaprofile_binned.tsv', sep='\t', index=False)
print(f"  Saved: tss_metaprofile_binned.tsv")

# ── 6. 図生成 ────────────────────────────────────────────────────────────────
print("\n[5] Generating figures...")

fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
fig.subplots_adjust(wspace=0.38, left=0.09, right=0.97, top=0.90, bottom=0.13)

# ─── Panel A: メタプロファイル ─────────────────────────────────────────────
ax = axes[0]
x = bin_centers

# Shielded (gray)
ax.fill_between(x, sm_lo_sh, sm_hi_sh, color=COL_SHIELDED, alpha=COL_ALPHA)
ax.plot(x, sm_sh, color=COL_SHIELDED, linewidth=2.0, label=f'Shielded (n={len(df_sh)})')

# Exposed (purple)
ax.fill_between(x, sm_lo_ex, sm_hi_ex, color=COL_EXPOSED, alpha=COL_ALPHA)
ax.plot(x, sm_ex, color=COL_EXPOSED, linewidth=2.0, label=f'Exposed (n={len(df_ex)})')

# TSS線
ax.axvline(0, color='black', linewidth=1.0, linestyle='--', alpha=0.6)
# +293bp境界線
ax.axvline(DOWN_BOUND, color='#E53935', linewidth=1.0, linestyle=':', alpha=0.8, label=f'+{DOWN_BOUND} bp boundary')
# 上流境界 (-500)
ax.axvspan(-HALF_WIN, 0, color='#1565C0', alpha=0.04, zorder=0)
# 下流領域 (0〜+293)
ax.axvspan(0, DOWN_BOUND, color='#E53935', alpha=0.04, zorder=0)

ax.set_xlabel('Distance from TSS (bp)', fontsize=10)
ax.set_ylabel('GCCGGC 4mC density\n(sites kb⁻¹ gene⁻¹)', fontsize=10)
ax.set_title('GCCGGC protection at TSS\n(T1, regulatory genes)', fontsize=10, fontweight='bold')
ax.set_xlim(-HALF_WIN, HALF_WIN)
ax.set_xticks([-500, -293, 0, 293, 500])
ax.set_xticklabels(['-500', '-293', 'TSS', '+293', '+500'], fontsize=8)

# 領域ラベル
ymax = ax.get_ylim()[1]
ax.text(-HALF_WIN/2, ymax * 0.92, 'Upstream', ha='center', fontsize=8, color='#1565C0', alpha=0.8)
ax.text(DOWN_BOUND/2, ymax * 0.92, 'Gene body\n(proximal)', ha='center', fontsize=7.5, color='#E53935', alpha=0.8)

ax.legend(fontsize=8, loc='upper right', framealpha=0.8)
ax.text(-0.12, 1.06, 'a', transform=ax.transAxes, fontsize=14, fontweight='bold', va='top')

# p値アノテーション
p_str_dn = f'p={p_dn:.2e}' if p_dn >= 1e-10 else 'p<1e-10'
ax.text(0.62, 0.10, f'Downstream\n{p_str_dn}',
        transform=ax.transAxes, fontsize=7.5, color='#E53935',
        ha='center', va='bottom',
        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7, edgecolor='#E53935', linewidth=0.5))

# ─── Panel B: 上流密度 vs 下流密度 散布図 ──────────────────────────────────
ax2 = axes[1]

jitter = np.random.RandomState(42).uniform(-0.005, 0.005, len(df_sh))
ax2.scatter(up_sh + jitter, dn_sh, color=COL_SHIELDED, alpha=0.55, s=22,
            label=f'Shielded (n={len(df_sh)})', linewidths=0, zorder=3)
ax2.scatter(up_ex, dn_ex, color=COL_EXPOSED, alpha=0.85, s=40,
            label=f'Exposed (n={len(df_ex)})', edgecolors='white', linewidths=0.5, zorder=5)

# 対角線 (upstream = downstream)
lim_max = max(np.nanmax(up_sh), np.nanmax(up_ex), np.nanmax(dn_sh), np.nanmax(dn_ex)) * 1.05
ax2.plot([0, lim_max], [0, lim_max], 'k--', linewidth=0.8, alpha=0.3, zorder=1)

# 各群の重心
ax2.scatter(np.median(up_sh), np.median(dn_sh), color=COL_SHIELDED, s=90,
            marker='D', edgecolors='black', linewidths=1.0, zorder=6)
ax2.scatter(np.median(up_ex), np.median(dn_ex), color=COL_EXPOSED, s=90,
            marker='D', edgecolors='black', linewidths=1.0, zorder=6)

ax2.set_xlabel(f'Upstream GCCGGC density\n(sites kb⁻¹, −500 to 0 bp)', fontsize=9.5)
ax2.set_ylabel(f'Downstream GCCGGC density\n(sites kb⁻¹, 0 to +{DOWN_BOUND} bp)', fontsize=9.5)
ax2.set_title('Upstream vs Downstream protection\n(per gene, T1)', fontsize=10, fontweight='bold')
ax2.set_xlim(-0.02, lim_max)
ax2.set_ylim(-0.02, lim_max)
ax2.legend(fontsize=8, loc='upper left', framealpha=0.8)
ax2.text(-0.12, 1.06, 'b', transform=ax2.transAxes, fontsize=14, fontweight='bold', va='top')

# 統計注釈
p_str_up = f'p={p_up:.2e}' if p_up >= 1e-10 else 'p<1e-10'
p_str_dn2 = f'p={p_dn:.2e}' if p_dn >= 1e-10 else 'p<1e-10'
textstr = (f'Upstream: {p_str_up}\nDownstream: {p_str_dn2}')
ax2.text(0.97, 0.04, textstr, transform=ax2.transAxes,
         fontsize=7.5, ha='right', va='bottom',
         bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8,
                   edgecolor='gray', linewidth=0.5))

# ── 保存 ─────────────────────────────────────────────────────────────────────
out_stem = FIG_DIR / 'tss_upstream_downstream_profile'
for fmt in ('png', 'pdf', 'svg'):
    fig.savefig(out_stem.with_suffix(f'.{fmt}'), format=fmt,
                dpi=300, bbox_inches='tight')
    print(f"  Saved: {out_stem.with_suffix(f'.{fmt}')}")

plt.close(fig)

# Supp figureとしてコピー
import shutil
supp_path = FIG_SUP / 'SuppFig_protection_downstream.png'
shutil.copy(out_stem.with_suffix('.png'), supp_path)
print(f"  Copied to: {supp_path}")

# ── 7. サマリー出力 ────────────────────────────────────────────────────────
print("\n" + "=" * 65)
print("SUMMARY")
print("=" * 65)

up_ratio_sh = np.median(up_sh)
dn_ratio_sh = np.median(dn_sh)
up_ratio_ex = np.median(up_ex)
dn_ratio_ex = np.median(dn_ex)

print(f"\nShielded genes (n={len(df_sh)}):")
print(f"  Median upstream density  : {up_ratio_sh:.4f} sites/kb")
print(f"  Median downstream density: {dn_ratio_sh:.4f} sites/kb")
print(f"  Ratio (up/down)          : {up_ratio_sh/(dn_ratio_sh+1e-9):.2f}")

print(f"\nExposed genes (n={len(df_ex)}):")
print(f"  Median upstream density  : {up_ratio_ex:.4f} sites/kb")
print(f"  Median downstream density: {dn_ratio_ex:.4f} sites/kb")
print(f"  Ratio (up/down)          : {up_ratio_ex/(dn_ratio_ex+1e-9):.2f}")

print(f"\nInterpretation:")
# 上流vs下流を比べる
sh_asym = up_ratio_sh / (dn_ratio_sh + 1e-9)
ex_asym = up_ratio_ex / (dn_ratio_ex + 1e-9)

if sh_asym > 1.5:
    print("  → Shielded: 上流保護が下流より顕著 (上流優位モデル)")
else:
    print("  → Shielded: 上流・下流ともに保護 (TSS周辺対称モデル)")

if p_dn < 0.05:
    print(f"  → Downstream depletion is SIGNIFICANT (p={p_dn:.3e})")
    print("     解釈: 保護ゾーンはTSS下流(+293bp)にも実在する")
else:
    print(f"  → Downstream depletion is NOT significant (p={p_dn:.3e})")
    print("     解釈: 保護は主にTSS上流に限定される")

print("\nDone.")
