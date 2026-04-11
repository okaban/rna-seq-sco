"""
64_timepoint_TSS: タイムポイント別TSS周辺メチル化密度解析
==========================================================
対応論点: C-2

方針:
  - T1/T2/T3 それぞれのタイムポイントについて、TSS周辺±5kbのメチル化密度プロファイルを計算
  - 4mCと6mAを別々に解析
  - 制御遺伝子 vs 非制御遺伝子での比較
  - 保護ゾーン（TSS直上流の低密度領域）の幅・深さのタイムポイント依存性
  - 既存の18_tss_analysesデータ（T1/T2のみ）にT3を追加して比較

実施日: 2026-04-11
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path

# === パス設定 ===
ANALYSIS_DIR = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/64_timepoint_TSS")
DATA_PATH = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/01_integration/high_confidence_sites_weighted.csv")
TSS_PATH = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/18_tss_analyses/jeong2016_all_tss.csv")
GENE_ANNOT = Path("/Users/okaban/bioinfo/rna-seq/05_annotation/analysis/05_annotation_260128_v1/tables/gene_annotation_basic.tsv")
GENE_SET = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/40_motif_division_of_labor/tables/gene_set_comparison.tsv")

OUT_TABLES = ANALYSIS_DIR / "tables"
OUT_FIGS   = ANALYSIS_DIR / "figures"

TIMEPOINTS = ['T1', 'T2', 'T3']
MOD_TYPES = ['4mC', '6mA']
WINDOW = 5000  # ±5kb
BIN_SIZE = 200  # 200bp bins
COLORS = {'T1': '#4C72B0', 'T2': '#DD8452', 'T3': '#55A868'}

# カテゴリ定義
REG_CATEGORIES = {'Regulatory/TF'}

def load_data():
    """データ読み込み"""
    methyl = pd.read_csv(DATA_PATH)
    tss = pd.read_csv(TSS_PATH)
    annot = pd.read_csv(GENE_ANNOT, sep='\t')
    gene_set = pd.read_csv(GENE_SET, sep='\t')
    return methyl, tss, annot, gene_set

def compute_metagene_profile(methyl_df, tss_df, timepoint, mod_type,
                              window=5000, bin_size=200):
    """
    TSS中心のメタジーンプロファイルを計算
    各TSSについて±window内のメチル化サイト密度をbin_size単位で集計
    """
    subset = methyl_df[(methyl_df['timepoint'] == timepoint) &
                       (methyl_df['mod_type'] == mod_type)].copy()

    bins = np.arange(-window, window + bin_size, bin_size)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    n_bins = len(bin_centers)

    # TSSごとにカウント
    all_counts = np.zeros(n_bins)
    n_tss = 0

    for _, tss_row in tss_df.iterrows():
        tss_pos = tss_row['tss_position']
        strand = tss_row['strand']

        # TSS相対座標に変換
        rel_pos = subset['position'] - tss_pos
        if strand == '-':
            rel_pos = -rel_pos  # マイナス鎖は反転

        # ウィンドウ内フィルタ
        in_window = (rel_pos >= -window) & (rel_pos < window)
        rel_in_window = rel_pos[in_window].values

        if len(rel_in_window) > 0:
            counts, _ = np.histogram(rel_in_window, bins=bins)
            all_counts += counts
        n_tss += 1

    # bins/kb/TSS に正規化
    bin_kb = bin_size / 1000
    density = all_counts / (n_tss * bin_kb)

    return pd.DataFrame({
        'bin_center': bin_centers,
        'density': density,
        'total_sites': all_counts,
        'n_tss': n_tss,
    })

def find_protection_zone(profile_df, upstream_range=(-500, 0)):
    """
    保護ゾーンの指標を計算
    TSS直上流の密度が下流より低い領域を保護ゾーンとみなす
    """
    upstream = profile_df[
        (profile_df['bin_center'] >= upstream_range[0]) &
        (profile_df['bin_center'] < upstream_range[1])
    ]['density'].mean()

    downstream = profile_df[
        (profile_df['bin_center'] >= 0) &
        (profile_df['bin_center'] < 500)
    ]['density'].mean()

    far_upstream = profile_df[
        (profile_df['bin_center'] >= -2000) &
        (profile_df['bin_center'] < -500)
    ]['density'].mean()

    ratio = upstream / far_upstream if far_upstream > 0 else np.nan
    return {
        'upstream_density': upstream,
        'downstream_density': downstream,
        'far_upstream_density': far_upstream,
        'protection_ratio': ratio,  # <1 = protected
    }

def plot_timepoint_profiles(profiles_dict, mod_type, gene_cat, out_path):
    """
    T1/T2/T3のメタジーンプロファイルを重ねてプロット
    profiles_dict: {tp: DataFrame}
    """
    fig, ax = plt.subplots(figsize=(10, 5))

    for tp in TIMEPOINTS:
        if tp not in profiles_dict:
            continue
        df = profiles_dict[tp]
        x = df['bin_center'] / 1000  # kb
        y = df['density']
        ax.plot(x, y, color=COLORS[tp], linewidth=2, label=tp)
        ax.fill_between(x, y, alpha=0.15, color=COLORS[tp])

    ax.axvline(0, color='red', linestyle='--', linewidth=1, label='TSS')
    ax.set_xlabel('Distance from TSS (kb)', fontsize=12)
    ax.set_ylabel('Methylation density (sites/kb/TSS)', fontsize=12)
    ax.set_title(f'{mod_type} methylation profile around TSS\n({gene_cat} genes)', fontsize=13)
    ax.legend()
    ax.set_xlim(-WINDOW / 1000, WINDOW / 1000)

    plt.tight_layout()
    plt.savefig(out_path, bbox_inches='tight')
    plt.savefig(str(out_path).replace('.pdf', '.svg'), format='svg', bbox_inches='tight')
    plt.close()
    print(f"  保存: {out_path}")

def plot_protection_metrics_barplot(metrics_df, out_path):
    """保護ゾーン指標のバープロット"""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for col, (ax, metric, ylabel) in enumerate(zip(
        axes,
        ['protection_ratio', 'upstream_density'],
        ['Protection ratio (upstream/far-upstream)', 'Upstream density (sites/kb/TSS)']
    )):
        pivot = metrics_df.pivot(index='timepoint', columns='mod_type', values=metric)
        pivot = pivot.reindex(TIMEPOINTS)
        x = np.arange(len(TIMEPOINTS))
        width = 0.35

        for i, mod in enumerate(['4mC', '6mA']):
            if mod in pivot.columns:
                vals = pivot[mod].values
                bars = ax.bar(x + i * width, vals, width, label=mod,
                              color=['#d62728', '#1f77b4'][i], alpha=0.75)
                for bar, val in zip(bars, vals):
                    if not np.isnan(val):
                        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.005,
                                f'{val:.3f}', ha='center', va='bottom', fontsize=8)

        ax.set_xticks(x + width / 2)
        ax.set_xticklabels(TIMEPOINTS)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(f'Protection zone: {metric}', fontsize=11)
        ax.legend()
        if metric == 'protection_ratio':
            ax.axhline(1.0, color='gray', linestyle='--', linewidth=0.8)

    plt.tight_layout()
    plt.savefig(out_path, bbox_inches='tight')
    plt.savefig(str(out_path).replace('.pdf', '.svg'), format='svg', bbox_inches='tight')
    plt.close()
    print(f"  保存: {out_path}")

# ============================================================
# メイン
# ============================================================
if __name__ == '__main__':
    print("=== 64_timepoint_TSS ===\n")

    # データ読み込み
    print("1. データ読み込み...")
    methyl, tss, annot, gene_set = load_data()

    # 制御遺伝子リスト
    reg_gene_ids = set(gene_set[gene_set['category'].isin(REG_CATEGORIES)]['gene_id'].astype(str))
    print(f"   TSS数: {len(tss)}")
    print(f"   制御遺伝子数: {len(reg_gene_ids)}")
    print(f"   メチル化サイト: {len(methyl)}")

    # TSS-gene_id マッピング（Primary TSSのみ）
    tss_primary = tss[tss['category'] == 'P'].copy()
    tss_primary['is_regulatory'] = tss_primary['gene_id'].astype(str).isin(reg_gene_ids)
    print(f"   Primary TSS数: {len(tss_primary)}")
    print(f"   制御遺伝子TSS: {tss_primary['is_regulatory'].sum()}")

    # 遺伝子カテゴリ別TSS
    tss_sets = {
        'All': tss_primary,
        'Regulatory': tss_primary[tss_primary['is_regulatory']],
        'Non-regulatory': tss_primary[~tss_primary['is_regulatory']],
    }

    # メタジーンプロファイル計算
    print("\n2. メタジーンプロファイル計算...")
    all_profiles = {}  # {(mod_type, gene_cat, tp): DataFrame}
    metrics_rows = []

    for mod in MOD_TYPES:
        for gene_cat, tss_subset in tss_sets.items():
            profiles = {}
            print(f"\n  [{mod}] {gene_cat} (n_TSS={len(tss_subset)})...")
            for tp in TIMEPOINTS:
                print(f"   {tp}...", end=' ', flush=True)
                prof = compute_metagene_profile(methyl, tss_subset, tp, mod,
                                                window=WINDOW, bin_size=BIN_SIZE)
                profiles[tp] = prof
                all_profiles[(mod, gene_cat, tp)] = prof

                # 保護ゾーン指標
                pz = find_protection_zone(prof)
                metrics_rows.append({
                    'mod_type': mod,
                    'gene_cat': gene_cat,
                    'timepoint': tp,
                    **pz,
                })
                print(f"ratio={pz['protection_ratio']:.3f}")

            # プロファイルプロット
            plot_timepoint_profiles(
                profiles, mod, gene_cat,
                OUT_FIGS / f"C2_profile_{mod}_{gene_cat.replace(' ', '_')}.pdf"
            )

    # 保護ゾーン指標テーブル
    metrics_df = pd.DataFrame(metrics_rows)
    metrics_df.to_csv(OUT_TABLES / "C2_protection_zone_metrics.tsv", sep='\t', index=False, float_format='%.5f')
    print(f"\n   保存: {OUT_TABLES}/C2_protection_zone_metrics.tsv")

    # 保護ゾーン指標の要約表示
    print("\n3. 保護ゾーン指標サマリー:")
    print(metrics_df[metrics_df['gene_cat'] == 'All'][
        ['mod_type', 'timepoint', 'protection_ratio', 'upstream_density', 'far_upstream_density']
    ].to_string(index=False))

    # バープロット
    metrics_all = metrics_df[metrics_df['gene_cat'] == 'All'].copy()
    plot_protection_metrics_barplot(
        metrics_all,
        OUT_FIGS / "C2_protection_zone_barplot_all.pdf"
    )

    # 制御 vs 非制御の比較プロット
    for mod in MOD_TYPES:
        fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
        for i, tp in enumerate(TIMEPOINTS):
            ax = axes[i]
            for gene_cat, color in [('Regulatory', '#d62728'), ('Non-regulatory', '#1f77b4')]:
                prof = all_profiles.get((mod, gene_cat, tp))
                if prof is None:
                    continue
                x = prof['bin_center'] / 1000
                y = prof['density']
                ax.plot(x, y, color=color, linewidth=1.8, label=gene_cat)
                ax.fill_between(x, y, alpha=0.12, color=color)
            ax.axvline(0, color='red', linestyle='--', linewidth=0.8)
            ax.set_title(tp, fontsize=12)
            ax.set_xlabel('Distance from TSS (kb)', fontsize=10)
            if i == 0:
                ax.set_ylabel('Methylation density (sites/kb/TSS)', fontsize=10)
            ax.legend(fontsize=8)
            ax.set_xlim(-WINDOW / 1000, WINDOW / 1000)
        fig.suptitle(f'{mod} methylation profile: Regulatory vs Non-regulatory', fontsize=13)
        plt.tight_layout()
        out = OUT_FIGS / f"C2_reg_vs_nonreg_{mod}_panel.pdf"
        plt.savefig(out, bbox_inches='tight')
        plt.savefig(str(out).replace('.pdf', '.svg'), format='svg', bbox_inches='tight')
        plt.close()
        print(f"  保存: {out}")

    # プロファイルデータをテーブルに保存
    dfs = []
    for (mod, gene_cat, tp), prof in all_profiles.items():
        prof = prof.copy()
        prof['mod_type'] = mod
        prof['gene_cat'] = gene_cat
        prof['timepoint'] = tp
        dfs.append(prof)
    all_prof_df = pd.concat(dfs, ignore_index=True)
    all_prof_df.to_csv(OUT_TABLES / "C2_all_profiles.tsv", sep='\t', index=False, float_format='%.6f')
    print(f"   保存: {OUT_TABLES}/C2_all_profiles.tsv")

    print("\n=== 完了 ===")
