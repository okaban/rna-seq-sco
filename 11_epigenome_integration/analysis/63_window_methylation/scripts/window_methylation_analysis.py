"""
63_window_methylation: ゲノムワイドウィンドウ単位メチル化密度解析
==================================================================
対応論点: B-1

方針:
  - ゲノムを一定幅（50kb）のウィンドウに分割
  - 各ウィンドウのメチル化サイト密度（sites/kb）をT1/T2/T3で集計
  - 4mCと6mAを別々に可視化
  - タイムポイント間の変化（T2-T1, T3-T2, T3-T1）を可視化
  - コア/アーム領域の注記（oriC位置）

実施日: 2026-04-11
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

# === パス設定 ===
ANALYSIS_DIR = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/63_window_methylation")
DATA_PATH = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/01_integration/high_confidence_sites_weighted.csv")

OUT_TABLES = ANALYSIS_DIR / "tables"
OUT_FIGS   = ANALYSIS_DIR / "figures"

# S. coelicolor A3(2) M145 ゲノムパラメータ
GENOME_SIZE = 8667507  # bp
CHROM = "NC_003888.3"
ORI_C = 4000000       # oriCおよびコア/アーム境界の概算（arm: <2Mb or >6Mb）
ARM_BOUNDARY_LEFT = 2_000_000
ARM_BOUNDARY_RIGHT = 6_000_000

WINDOW_SIZE = 50_000  # 50 kb

TIMEPOINTS = ['T1', 'T2', 'T3']
MOD_TYPES = ['4mC', '6mA']

COLORS = {
    'T1': '#4C72B0',
    'T2': '#DD8452',
    'T3': '#55A868',
}

def make_windows(genome_size, window_size):
    """ゲノムウィンドウ定義"""
    starts = np.arange(0, genome_size, window_size)
    ends = np.minimum(starts + window_size, genome_size)
    centers = (starts + ends) / 2
    return pd.DataFrame({'start': starts, 'end': ends, 'center': centers})

def count_sites_per_window(df_mod, windows, timepoint, mod_type):
    """指定タイムポイント・モディフィケーションのウィンドウ別サイト数"""
    subset = df_mod[(df_mod['timepoint'] == timepoint) & (df_mod['mod_type'] == mod_type)].copy()
    counts = []
    for _, w in windows.iterrows():
        n = ((subset['position'] >= w['start']) & (subset['position'] < w['end'])).sum()
        window_kb = (w['end'] - w['start']) / 1000
        counts.append(n / window_kb)  # sites per kb
    return np.array(counts)

def mean_freq_per_window(df_mod, windows, timepoint, mod_type):
    """指定タイムポイント・モディフィケーションのウィンドウ別平均メチル化頻度"""
    subset = df_mod[(df_mod['timepoint'] == timepoint) & (df_mod['mod_type'] == mod_type)].copy()
    means = []
    for _, w in windows.iterrows():
        mask = (subset['position'] >= w['start']) & (subset['position'] < w['end'])
        vals = subset.loc[mask, 'weighted_mod_freq']
        means.append(vals.mean() if len(vals) > 0 else np.nan)
    return np.array(means)

def add_arm_core_shading(ax, arm_left, arm_right, alpha=0.08):
    """コア/アーム領域のシェーディング"""
    ax.axvspan(0, arm_left / 1e6, color='gray', alpha=alpha, label='Arm region')
    ax.axvspan(arm_right / 1e6, GENOME_SIZE / 1e6, color='gray', alpha=alpha)

def plot_density_comparison(windows, density_dict, mod_type, out_path):
    """T1/T2/T3のサイト密度を重ねてプロット"""
    fig, axes = plt.subplots(2, 1, figsize=(18, 8), sharex=True)
    x = windows['center'].values / 1e6  # Mb

    # 上段: 密度
    ax = axes[0]
    for tp in TIMEPOINTS:
        y = density_dict[tp]
        ax.fill_between(x, y, alpha=0.4, color=COLORS[tp], label=tp)
        ax.plot(x, y, linewidth=0.5, color=COLORS[tp])
    add_arm_core_shading(ax, ARM_BOUNDARY_LEFT, ARM_BOUNDARY_RIGHT)
    ax.set_ylabel('Sites per kb', fontsize=11)
    ax.set_title(f'{mod_type} methylation site density (50kb windows)', fontsize=13)
    ax.legend(loc='upper right')
    ax.set_xlim(0, GENOME_SIZE / 1e6)

    # 下段: T3-T1 差分
    ax2 = axes[1]
    delta = density_dict['T3'] - density_dict['T1']
    colors_bar = ['#d62728' if v > 0 else '#1f77b4' for v in delta]
    ax2.bar(x, delta, width=WINDOW_SIZE / 1e6, color=colors_bar, alpha=0.7)
    add_arm_core_shading(ax2, ARM_BOUNDARY_LEFT, ARM_BOUNDARY_RIGHT)
    ax2.axhline(0, color='black', linewidth=0.8)
    ax2.set_ylabel('ΔDensity (T3−T1)', fontsize=11)
    ax2.set_xlabel('Genomic position (Mb)', fontsize=11)
    ax2.set_title(f'{mod_type} density change T3−T1', fontsize=12)

    plt.tight_layout()
    plt.savefig(out_path, bbox_inches='tight')
    plt.savefig(str(out_path).replace('.pdf', '.svg'), format='svg', bbox_inches='tight')
    plt.close()
    print(f"  保存: {out_path}")

def plot_freq_heatmap(windows, freq_dict, mod_type, out_path):
    """T1/T2/T3の平均メチル化頻度をヒートマップ表示"""
    data = np.vstack([freq_dict['T1'], freq_dict['T2'], freq_dict['T3']])
    x = windows['center'].values / 1e6

    fig, ax = plt.subplots(figsize=(18, 3))
    im = ax.imshow(data, aspect='auto', cmap='YlOrRd',
                   extent=[0, GENOME_SIZE / 1e6, -0.5, 2.5],
                   vmin=0, vmax=100, origin='lower')
    ax.set_yticks([0, 1, 2])
    ax.set_yticklabels(['T1', 'T2', 'T3'])
    ax.set_xlabel('Genomic position (Mb)', fontsize=11)
    ax.set_title(f'{mod_type} mean methylation frequency (%) – 50kb windows', fontsize=12)
    plt.colorbar(im, ax=ax, label='Mean freq (%)', shrink=0.8)

    # Arm shading
    for ymin, ymax in [(-0.5, 2.5)]:
        ax.axvspan(0, ARM_BOUNDARY_LEFT / 1e6, color='blue', alpha=0.05)
        ax.axvspan(ARM_BOUNDARY_RIGHT / 1e6, GENOME_SIZE / 1e6, color='blue', alpha=0.05)

    plt.tight_layout()
    plt.savefig(out_path, bbox_inches='tight')
    plt.savefig(str(out_path).replace('.pdf', '.svg'), format='svg', bbox_inches='tight')
    plt.close()
    print(f"  保存: {out_path}")

def plot_timepoint_comparison_panel(windows, density_dict_4mC, density_dict_6mA, out_path):
    """4mC/6mA並列パネル、T1→T2→T3変化"""
    fig, axes = plt.subplots(2, 3, figsize=(22, 8), sharex=True, sharey='row')
    x = windows['center'].values / 1e6
    comparisons = [('T1', 'T2'), ('T2', 'T3'), ('T1', 'T3')]

    for col, (tp_a, tp_b) in enumerate(comparisons):
        for row, (mod_type, d_dict) in enumerate([('4mC', density_dict_4mC), ('6mA', density_dict_6mA)]):
            ax = axes[row][col]
            delta = d_dict[tp_b] - d_dict[tp_a]
            colors_bar = ['#d62728' if v > 0 else '#1f77b4' for v in delta]
            ax.bar(x, delta, width=WINDOW_SIZE / 1e6 * 0.9, color=colors_bar, alpha=0.75)
            ax.axhline(0, color='black', linewidth=0.8)
            add_arm_core_shading(ax, ARM_BOUNDARY_LEFT, ARM_BOUNDARY_RIGHT)
            ax.set_title(f'{mod_type}: {tp_b}−{tp_a}', fontsize=11)
            if col == 0:
                ax.set_ylabel('Δ sites/kb', fontsize=10)
            if row == 1:
                ax.set_xlabel('Position (Mb)', fontsize=10)

    plt.suptitle('Methylation density changes across timepoints (50kb windows)', fontsize=13, y=1.01)
    plt.tight_layout()
    plt.savefig(out_path, bbox_inches='tight')
    plt.savefig(str(out_path).replace('.pdf', '.svg'), format='svg', bbox_inches='tight')
    plt.close()
    print(f"  保存: {out_path}")

# ============================================================
# メイン
# ============================================================
if __name__ == '__main__':
    print("=== 63_window_methylation ===\n")

    # データ読み込み
    print("1. データ読み込み...")
    df = pd.read_csv(DATA_PATH)
    print(f"   総サイト数: {len(df)}")
    print(f"   タイムポイント分布:\n{df.groupby(['timepoint','mod_type']).size().to_string()}")

    # ウィンドウ作成
    print(f"\n2. ウィンドウ設定 (size={WINDOW_SIZE//1000}kb)...")
    windows = make_windows(GENOME_SIZE, WINDOW_SIZE)
    print(f"   ウィンドウ数: {len(windows)}")

    # 密度計算
    print("\n3. ウィンドウ別密度計算...")
    density = {}
    freq = {}
    for mod in MOD_TYPES:
        density[mod] = {}
        freq[mod] = {}
        for tp in TIMEPOINTS:
            print(f"   {mod} {tp}...")
            density[mod][tp] = count_sites_per_window(df, windows, tp, mod)
            freq[mod][tp] = mean_freq_per_window(df, windows, tp, mod)

    # テーブル出力
    print("\n4. テーブル出力...")
    for mod in MOD_TYPES:
        out_df = windows.copy()
        for tp in TIMEPOINTS:
            out_df[f'density_{tp}'] = density[mod][tp]
            out_df[f'mean_freq_{tp}'] = freq[mod][tp]
        out_df['delta_T3_T1'] = out_df['density_T3'] - out_df['density_T1']
        out_path = OUT_TABLES / f"B1_window_methylation_{mod}.tsv"
        out_df.to_csv(out_path, sep='\t', index=False, float_format='%.4f')
        print(f"   保存: {out_path}")

    # 統計サマリー
    print("\n5. 統計サマリー...")
    for mod in MOD_TYPES:
        print(f"\n  [{mod}]")
        total_sites = {}
        for tp in TIMEPOINTS:
            n = len(df[(df['mod_type'] == mod) & (df['timepoint'] == tp)])
            total_sites[tp] = n
            print(f"   {tp}: {n} sites total, mean density={np.nanmean(density[mod][tp]):.3f} sites/kb")
        # 最高密度ウィンドウ
        for tp in TIMEPOINTS:
            idx = np.argmax(density[mod][tp])
            w = windows.iloc[idx]
            print(f"   {tp} 最高密度ウィンドウ: {w['start']/1e6:.2f}-{w['end']/1e6:.2f}Mb "
                  f"({density[mod][tp][idx]:.2f} sites/kb)")

    # 図出力
    print("\n6. 図出力...")
    for mod in MOD_TYPES:
        plot_density_comparison(
            windows, density[mod], mod,
            OUT_FIGS / f"B1_density_comparison_{mod}.pdf"
        )
        plot_freq_heatmap(
            windows, freq[mod], mod,
            OUT_FIGS / f"B1_freq_heatmap_{mod}.pdf"
        )

    plot_timepoint_comparison_panel(
        windows, density['4mC'], density['6mA'],
        OUT_FIGS / "B1_timepoint_changes_panel.pdf"
    )

    # コア vs アーム 集計
    print("\n7. コア/アーム別集計...")
    core_mask = (windows['center'] >= ARM_BOUNDARY_LEFT) & (windows['center'] <= ARM_BOUNDARY_RIGHT)
    arm_mask = ~core_mask

    summary_rows = []
    for mod in MOD_TYPES:
        for tp in TIMEPOINTS:
            d = density[mod][tp]
            summary_rows.append({
                'mod_type': mod, 'timepoint': tp, 'region': 'core',
                'mean_density': np.nanmean(d[core_mask]),
                'total_windows': core_mask.sum(),
            })
            summary_rows.append({
                'mod_type': mod, 'timepoint': tp, 'region': 'arm',
                'mean_density': np.nanmean(d[arm_mask]),
                'total_windows': arm_mask.sum(),
            })
    summary_df = pd.DataFrame(summary_rows)
    print(summary_df.to_string(index=False))
    summary_df.to_csv(OUT_TABLES / "B1_core_arm_density_summary.tsv", sep='\t', index=False, float_format='%.4f')

    print("\n=== 完了 ===")
