"""
67_spatial_integration: DEG vs DNAメチル化サイトの空間的統合解析
================================================================
目的: DEGとメチル化サイトの位置的関係を定量化し、
     「メチル化は発現の直接制御因子ではない」という主張を空間的データで補強する。

解析1: DEG vs non-DEG の最近傍メチル化距離比較
解析2: 上昇DEG vs 低下DEG のプロモーターメチル化密度比較
解析3: 時系列発現変化とプロモーターメチル化の相関

実施日: 2026-04-19
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# scikit-posthocsが無い場合はDunn's testを手動実装
try:
    import scikit_posthocs as sp
    HAS_SP = True
except ImportError:
    HAS_SP = False

# ============================================================
# パス設定
# ============================================================
BASE_DIR   = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/67_spatial_integration")
OUT_TABLES = BASE_DIR / "tables"
OUT_FIGS   = BASE_DIR / "figures"
OUT_TABLES.mkdir(parents=True, exist_ok=True)
OUT_FIGS.mkdir(parents=True, exist_ok=True)

METHYL_PATH = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/01_integration/high_confidence_sites_weighted.csv")
TSS_PATH    = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/18_tss_analyses/comprehensive_tss_table.csv")
DEG_DIR     = Path("/Users/okaban/bioinfo/rna-seq/04_deseq2/analysis/04_deseq2_260128_v1/results")
COUNTS_PATH = DEG_DIR / "normalized_counts_M145.tsv"

COMPARISONS = {
    "T2vsT1": "DESeq2_M145_2_vs_1.tsv",
    "T3vsT1": "DESeq2_M145_3_vs_1.tsv",
    "T3vsT2": "DESeq2_M145_3_vs_2.tsv",
}

PADJ_CUTOFF = 0.05
PROMOTER_WINDOW = 500  # TSS ±500bp

# カラーパレット
COLORS = {
    "non-DEG":  "#AAAAAA",
    "up":       "#D62728",
    "down":     "#1F77B4",
    "4mC":      "#E377C2",
    "6mA":      "#17BECF",
    "T1":       "#4C72B0",
    "T2":       "#DD8452",
    "T3":       "#55A868",
    "mono_inc": "#D62728",
    "mono_dec": "#1F77B4",
    "other":    "#AAAAAA",
}

# ============================================================
# データ読み込み
# ============================================================
def load_data():
    print("[INFO] データ読み込み中...")
    methyl = pd.read_csv(METHYL_PATH)
    tss    = pd.read_csv(TSS_PATH)
    counts = pd.read_csv(COUNTS_PATH, sep='\t')

    degs = {}
    for label, fname in COMPARISONS.items():
        df = pd.read_csv(DEG_DIR / fname, sep='\t')
        degs[label] = df

    print(f"  メチル化サイト: {len(methyl):,} 行")
    print(f"  TSS: {len(tss):,} 遺伝子")
    print(f"  DEG (T2vsT1): {len(degs['T2vsT1']):,} 遺伝子")
    return methyl, tss, counts, degs


def get_deg_status(df_deg, padj_cutoff=PADJ_CUTOFF):
    """DEGをup/down/non-DEGに分類した Series を返す"""
    status = {}
    for _, row in df_deg.iterrows():
        g = row['gene_id']
        if pd.notna(row['padj']) and row['padj'] < padj_cutoff:
            status[g] = 'up' if row['log2FoldChange'] > 0 else 'down'
        else:
            status[g] = 'non-DEG'
    return status


# ============================================================
# 統計ヘルパー
# ============================================================
def rank_biserial_r(x, y):
    """Wilcoxon rank-biserial correlation (effect size)"""
    n1, n2 = len(x), len(y)
    u_stat, _ = stats.mannwhitneyu(x, y, alternative='two-sided')
    r = 1 - (2 * u_stat) / (n1 * n2)
    return r


def dunn_test_manual(groups: dict):
    """簡易Dunn's検定 (Bonferroni補正) — scikit_posthocs不要版"""
    import itertools
    all_vals = np.concatenate(list(groups.values()))
    n_total  = len(all_vals)
    ranks    = stats.rankdata(all_vals)

    start = 0
    group_ranks = {}
    for k, v in groups.items():
        group_ranks[k] = ranks[start:start + len(v)]
        start += len(v)

    pairs = list(itertools.combinations(groups.keys(), 2))
    results = []
    for g1, g2 in pairs:
        n1 = len(groups[g1])
        n2 = len(groups[g2])
        mean_r1 = group_ranks[g1].mean()
        mean_r2 = group_ranks[g2].mean()
        se = np.sqrt((n_total * (n_total + 1) / 12) * (1 / n1 + 1 / n2))
        z = (mean_r1 - mean_r2) / se
        p = 2 * stats.norm.sf(abs(z))
        results.append({"group1": g1, "group2": g2, "z": z, "p_raw": p,
                         "p_bonf": min(p * len(pairs), 1.0)})
    return pd.DataFrame(results)


# ============================================================
# 解析1: DEG vs non-DEG 最近傍メチル化距離比較
# ============================================================
def analysis1_nearest_distance(methyl, tss, degs):
    print("\n[解析1] DEG vs non-DEG 最近傍メチル化距離比較")

    # 4mCサイトのみ (主解析) と6mAサイト (補助) を分離
    methyl_4mC = methyl[methyl['mod_type'] == '4mC'].copy()
    methyl_6mA = methyl[methyl['mod_type'] == '6mA'].copy()

    # 全タイムポイント合算した位置セット（重複除去）
    sites_4mC = methyl_4mC[['chrom', 'position']].drop_duplicates()
    sites_6mA = methyl_6mA[['chrom', 'position']].drop_duplicates()

    # TSS座標を取得
    tss_df = tss[['gene_id', 'chrom', 'tss']].dropna(subset=['tss']).copy()
    tss_df['tss'] = tss_df['tss'].astype(int)

    results_all = []

    for comp_label, deg_df in degs.items():
        print(f"  比較: {comp_label}")
        status_map = get_deg_status(deg_df)

        for mod_label, sites in [('4mC', sites_4mC), ('6mA', sites_6mA)]:
            # 染色体ごとに最近傍計算
            dists = {}
            for chrom, grp in tss_df.groupby('chrom'):
                s = sites[sites['chrom'] == chrom]['position'].values
                if len(s) == 0:
                    for _, row in grp.iterrows():
                        dists[row['gene_id']] = np.nan
                    continue
                s_sorted = np.sort(s)
                for _, row in grp.iterrows():
                    tss_pos = row['tss']
                    idx = np.searchsorted(s_sorted, tss_pos)
                    candidates = []
                    if idx > 0:
                        candidates.append(abs(tss_pos - s_sorted[idx - 1]))
                    if idx < len(s_sorted):
                        candidates.append(abs(tss_pos - s_sorted[idx]))
                    dists[row['gene_id']] = min(candidates) if candidates else np.nan

            tss_df['nn_dist'] = tss_df['gene_id'].map(dists)
            tss_df['deg_status'] = tss_df['gene_id'].map(status_map).fillna('non-DEG')

            for status in ['up', 'down', 'non-DEG']:
                vals = tss_df[tss_df['deg_status'] == status]['nn_dist'].dropna().values
                results_all.append({
                    'comparison': comp_label,
                    'mod_type': mod_label,
                    'deg_status': status,
                    'n': len(vals),
                    'median_dist': np.median(vals) if len(vals) > 0 else np.nan,
                    'mean_dist': np.mean(vals) if len(vals) > 0 else np.nan,
                })

    results_df = pd.DataFrame(results_all)

    # 統計検定: DEG(up+down合算) vs non-DEG per comparison per mod_type
    stat_rows = []
    for comp_label, deg_df in degs.items():
        status_map = get_deg_status(deg_df)
        for mod_label, sites in [('4mC', sites_4mC), ('6mA', sites_6mA)]:
            dists = {}
            for chrom, grp in tss_df.groupby('chrom'):
                s = sites[sites['chrom'] == chrom]['position'].values
                if len(s) == 0:
                    for _, row in grp.iterrows():
                        dists[row['gene_id']] = np.nan
                    continue
                s_sorted = np.sort(s)
                for _, row in grp.iterrows():
                    tss_pos = row['tss']
                    idx = np.searchsorted(s_sorted, tss_pos)
                    candidates = []
                    if idx > 0:
                        candidates.append(abs(tss_pos - s_sorted[idx - 1]))
                    if idx < len(s_sorted):
                        candidates.append(abs(tss_pos - s_sorted[idx]))
                    dists[row['gene_id']] = min(candidates) if candidates else np.nan

            tss_df2 = tss_df.copy()
            tss_df2['nn_dist'] = tss_df2['gene_id'].map(dists)
            tss_df2['deg_status'] = tss_df2['gene_id'].map(status_map).fillna('non-DEG')

            deg_vals    = tss_df2[tss_df2['deg_status'].isin(['up', 'down'])]['nn_dist'].dropna().values
            nondeg_vals = tss_df2[tss_df2['deg_status'] == 'non-DEG']['nn_dist'].dropna().values
            up_vals     = tss_df2[tss_df2['deg_status'] == 'up']['nn_dist'].dropna().values
            down_vals   = tss_df2[tss_df2['deg_status'] == 'down']['nn_dist'].dropna().values

            if len(deg_vals) > 1 and len(nondeg_vals) > 1:
                stat, pval = stats.mannwhitneyu(deg_vals, nondeg_vals, alternative='two-sided')
                r = rank_biserial_r(deg_vals, nondeg_vals)
                stat_rows.append({
                    'comparison': comp_label, 'mod_type': mod_label,
                    'test': 'Wilcoxon(DEG vs nonDEG)',
                    'n_group1': len(deg_vals), 'n_group2': len(nondeg_vals),
                    'median_g1': np.median(deg_vals), 'median_g2': np.median(nondeg_vals),
                    'U_stat': stat, 'p_value': pval, 'rank_biserial_r': r,
                })

    stat_df = pd.DataFrame(stat_rows)
    stat_df.to_csv(OUT_TABLES / "A1_nearest_dist_stats.tsv", sep='\t', index=False)
    results_df.to_csv(OUT_TABLES / "A1_nearest_dist_summary.tsv", sep='\t', index=False)
    print(f"  -> テーブル保存完了")

    # ========== 可視化 ==========
    # 再計算して可視化用データを収集
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle("Analysis 1: Nearest Methylation Site Distance\nDEG vs non-DEG", fontsize=14, fontweight='bold')

    for col_idx, (comp_label, deg_df) in enumerate(degs.items()):
        status_map = get_deg_status(deg_df)
        for row_idx, (mod_label, sites) in enumerate([('4mC', sites_4mC), ('6mA', sites_6mA)]):
            ax = axes[row_idx, col_idx]

            dists = {}
            for chrom, grp in tss_df.groupby('chrom'):
                s = sites[sites['chrom'] == chrom]['position'].values
                if len(s) == 0:
                    for _, row in grp.iterrows():
                        dists[row['gene_id']] = np.nan
                    continue
                s_sorted = np.sort(s)
                for _, row in grp.iterrows():
                    tss_pos = row['tss']
                    idx = np.searchsorted(s_sorted, tss_pos)
                    candidates = []
                    if idx > 0:
                        candidates.append(abs(tss_pos - s_sorted[idx - 1]))
                    if idx < len(s_sorted):
                        candidates.append(abs(tss_pos - s_sorted[idx]))
                    dists[row['gene_id']] = min(candidates) if candidates else np.nan

            tmp = tss_df.copy()
            tmp['nn_dist'] = tmp['gene_id'].map(dists)
            tmp['deg_status'] = tmp['gene_id'].map(status_map).fillna('non-DEG')

            # violin plot
            groups = ['non-DEG', 'up', 'down']
            plot_data = [tmp[tmp['deg_status'] == g]['nn_dist'].dropna().values for g in groups]
            plot_data_log = [np.log10(d + 1) for d in plot_data]

            parts = ax.violinplot(plot_data_log, positions=[0, 1, 2], showmedians=True, showextrema=False)
            col_list = [COLORS['non-DEG'], COLORS['up'], COLORS['down']]
            for pc, c in zip(parts['bodies'], col_list):
                pc.set_facecolor(c)
                pc.set_alpha(0.7)
            parts['cmedians'].set_color('black')
            parts['cmedians'].set_linewidth(2)

            ax.set_xticks([0, 1, 2])
            ax.set_xticklabels([f"non-DEG\n(n={len(plot_data[0])})",
                                 f"up\n(n={len(plot_data[1])})",
                                 f"down\n(n={len(plot_data[2])})"])
            ax.set_ylabel("log10(Distance + 1) [bp]")
            ax.set_title(f"{comp_label} | {mod_label}")

            # p値を追加
            if comp_label in stat_df['comparison'].values:
                row_stat = stat_df[(stat_df['comparison'] == comp_label) &
                                   (stat_df['mod_type'] == mod_label)]
                if len(row_stat) > 0:
                    pv = row_stat.iloc[0]['p_value']
                    r  = row_stat.iloc[0]['rank_biserial_r']
                    pstr = f"p={pv:.3f}" if pv >= 0.001 else f"p={pv:.2e}"
                    ax.text(0.5, 0.97, f"Wilcoxon: {pstr}\nr={r:.3f}",
                            transform=ax.transAxes, ha='center', va='top',
                            fontsize=8, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()
    fig.savefig(OUT_FIGS / "A1_nearest_dist_violin.png", dpi=300, bbox_inches='tight')
    plt.close()

    # CDF plot (T2vsT1, 4mC のみ)
    fig2, axes2 = plt.subplots(1, 3, figsize=(15, 5))
    fig2.suptitle("Analysis 1: Cumulative Distribution of Nearest Methylation Distance", fontsize=13, fontweight='bold')

    for col_idx, (comp_label, deg_df) in enumerate(degs.items()):
        ax = axes2[col_idx]
        status_map = get_deg_status(deg_df)

        for mod_label, sites in [('4mC', sites_4mC)]:
            dists = {}
            for chrom, grp in tss_df.groupby('chrom'):
                s = sites[sites['chrom'] == chrom]['position'].values
                if len(s) == 0:
                    for _, row in grp.iterrows():
                        dists[row['gene_id']] = np.nan
                    continue
                s_sorted = np.sort(s)
                for _, row in grp.iterrows():
                    tss_pos = row['tss']
                    idx = np.searchsorted(s_sorted, tss_pos)
                    candidates = []
                    if idx > 0:
                        candidates.append(abs(tss_pos - s_sorted[idx - 1]))
                    if idx < len(s_sorted):
                        candidates.append(abs(tss_pos - s_sorted[idx]))
                    dists[row['gene_id']] = min(candidates) if candidates else np.nan

            tmp = tss_df.copy()
            tmp['nn_dist'] = tmp['gene_id'].map(dists)
            tmp['deg_status'] = tmp['gene_id'].map(status_map).fillna('non-DEG')

        for status, color in [('non-DEG', COLORS['non-DEG']), ('up', COLORS['up']), ('down', COLORS['down'])]:
            vals = np.sort(tmp[tmp['deg_status'] == status]['nn_dist'].dropna().values)
            cdf  = np.arange(1, len(vals) + 1) / len(vals)
            ax.plot(vals, cdf, color=color, label=f"{status} (n={len(vals)})", linewidth=2)

        ax.set_xscale('log')
        ax.set_xlabel("Distance to nearest 4mC site [bp]")
        ax.set_ylabel("Cumulative fraction")
        ax.set_title(f"{comp_label}")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig2.savefig(OUT_FIGS / "A1_nearest_dist_cdf.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("  -> 図保存完了")

    return stat_df


# ============================================================
# 解析2: up DEG vs down DEG vs non-DEG プロモーターメチル化密度
# ============================================================
def analysis2_promoter_density(methyl, tss, degs):
    print("\n[解析2] プロモーターメチル化密度 (up vs down vs non-DEG)")

    # プロモーター領域 (TSS ±500bp) のメチル化サイト数をカウント
    tss_df = tss[['gene_id', 'chrom', 'tss', 'strand']].dropna(subset=['tss']).copy()
    tss_df['tss'] = tss_df['tss'].astype(int)

    stat_rows = []
    all_density = []

    for comp_label, deg_df in degs.items():
        print(f"  比較: {comp_label}")
        status_map = get_deg_status(deg_df)

        for mod_label in ['4mC', '6mA']:
            sites = methyl[methyl['mod_type'] == mod_label][['chrom', 'position']].drop_duplicates()

            # 各遺伝子のプロモーター密度
            densities = {}
            for chrom, grp in tss_df.groupby('chrom'):
                s = sites[sites['chrom'] == chrom]['position'].values
                for _, row in grp.iterrows():
                    tss_pos = row['tss']
                    lo = tss_pos - PROMOTER_WINDOW
                    hi = tss_pos + PROMOTER_WINDOW
                    n_sites = np.sum((s >= lo) & (s <= hi))
                    # density: sites per kb
                    density = n_sites / (2 * PROMOTER_WINDOW / 1000)
                    densities[row['gene_id']] = density

            tmp = tss_df.copy()
            tmp['density'] = tmp['gene_id'].map(densities)
            tmp['deg_status'] = tmp['gene_id'].map(status_map).fillna('non-DEG')
            tmp['comparison'] = comp_label
            tmp['mod_type']   = mod_label
            all_density.append(tmp)

            # Kruskal-Wallis
            groups_data = {
                'up':      tmp[tmp['deg_status'] == 'up']['density'].dropna().values,
                'down':    tmp[tmp['deg_status'] == 'down']['density'].dropna().values,
                'non-DEG': tmp[tmp['deg_status'] == 'non-DEG']['density'].dropna().values,
            }
            if all(len(v) > 0 for v in groups_data.values()):
                h_stat, kw_pval = stats.kruskal(*groups_data.values())

                stat_rows.append({
                    'comparison': comp_label, 'mod_type': mod_label,
                    'test': 'Kruskal-Wallis',
                    'H_stat': h_stat, 'p_value': kw_pval,
                    'n_up': len(groups_data['up']),
                    'n_down': len(groups_data['down']),
                    'n_nonDEG': len(groups_data['non-DEG']),
                    'median_up': np.median(groups_data['up']),
                    'median_down': np.median(groups_data['down']),
                    'median_nonDEG': np.median(groups_data['non-DEG']),
                })

                # Dunn's test (post-hoc)
                dunn = dunn_test_manual(groups_data)
                dunn['comparison'] = comp_label
                dunn['mod_type'] = mod_label
                stat_rows_dunn = dunn.to_dict('records')

                print(f"    {mod_label}: KW H={h_stat:.3f}, p={kw_pval:.4f}")
                for r in stat_rows_dunn:
                    print(f"      Dunn {r['group1']} vs {r['group2']}: z={r['z']:.3f}, p_bonf={r['p_bonf']:.4f}")

    density_df = pd.concat(all_density, ignore_index=True)
    stat_df    = pd.DataFrame([r for r in stat_rows if 'H_stat' in r])
    stat_df.to_csv(OUT_TABLES / "A2_promoter_density_stats.tsv", sep='\t', index=False)

    # Dunn's test テーブルを別途保存
    dunn_rows = []
    for comp_label, deg_df in degs.items():
        status_map = get_deg_status(deg_df)
        for mod_label in ['4mC', '6mA']:
            sites = methyl[methyl['mod_type'] == mod_label][['chrom', 'position']].drop_duplicates()
            densities = {}
            for chrom, grp in tss_df.groupby('chrom'):
                s = sites[sites['chrom'] == chrom]['position'].values
                for _, row in grp.iterrows():
                    tss_pos = row['tss']
                    lo = tss_pos - PROMOTER_WINDOW
                    hi = tss_pos + PROMOTER_WINDOW
                    n_sites = np.sum((s >= lo) & (s <= hi))
                    densities[row['gene_id']] = n_sites / (2 * PROMOTER_WINDOW / 1000)
            tmp = tss_df.copy()
            tmp['density'] = tmp['gene_id'].map(densities)
            tmp['deg_status'] = tmp['gene_id'].map(status_map).fillna('non-DEG')
            groups_data = {
                'up':      tmp[tmp['deg_status'] == 'up']['density'].dropna().values,
                'down':    tmp[tmp['deg_status'] == 'down']['density'].dropna().values,
                'non-DEG': tmp[tmp['deg_status'] == 'non-DEG']['density'].dropna().values,
            }
            if all(len(v) > 0 for v in groups_data.values()):
                dunn = dunn_test_manual(groups_data)
                dunn['comparison'] = comp_label
                dunn['mod_type'] = mod_label
                dunn_rows.append(dunn)
    if dunn_rows:
        pd.concat(dunn_rows, ignore_index=True).to_csv(
            OUT_TABLES / "A2_dunn_test.tsv", sep='\t', index=False)

    # ========== 可視化 ==========
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle(f"Analysis 2: Promoter Methylation Density (TSS ±{PROMOTER_WINDOW}bp)\n"
                  "up-DEG vs down-DEG vs non-DEG", fontsize=14, fontweight='bold')

    for col_idx, comp_label in enumerate(COMPARISONS.keys()):
        deg_df = degs[comp_label]
        status_map = get_deg_status(deg_df)

        for row_idx, mod_label in enumerate(['4mC', '6mA']):
            ax = axes[row_idx, col_idx]
            sites = methyl[methyl['mod_type'] == mod_label][['chrom', 'position']].drop_duplicates()

            densities = {}
            for chrom, grp in tss_df.groupby('chrom'):
                s = sites[sites['chrom'] == chrom]['position'].values
                for _, row in grp.iterrows():
                    tss_pos = row['tss']
                    lo = tss_pos - PROMOTER_WINDOW
                    hi = tss_pos + PROMOTER_WINDOW
                    densities[row['gene_id']] = np.sum((s >= lo) & (s <= hi)) / (2 * PROMOTER_WINDOW / 1000)

            tmp = tss_df.copy()
            tmp['density'] = tmp['gene_id'].map(densities)
            tmp['deg_status'] = tmp['gene_id'].map(status_map).fillna('non-DEG')

            groups = ['non-DEG', 'up', 'down']
            plot_data = [tmp[tmp['deg_status'] == g]['density'].dropna().values + 1e-6
                         for g in groups]

            bp = ax.boxplot(plot_data, positions=[0, 1, 2], patch_artist=True,
                            showfliers=False, medianprops=dict(color='black', linewidth=2))
            for patch, c in zip(bp['boxes'], [COLORS['non-DEG'], COLORS['up'], COLORS['down']]):
                patch.set_facecolor(c)
                patch.set_alpha(0.7)

            ax.set_xticks([0, 1, 2])
            ax.set_xticklabels([f"non-DEG\n(n={len(plot_data[0])})",
                                 f"up\n(n={len(plot_data[1])})",
                                 f"down\n(n={len(plot_data[2])})"])
            ax.set_ylabel("Methylation density [sites/kb]")
            ax.set_title(f"{comp_label} | {mod_label}")

            # KW p値を表示
            kw_row = stat_df[(stat_df['comparison'] == comp_label) & (stat_df['mod_type'] == mod_label)]
            if len(kw_row) > 0:
                pv = kw_row.iloc[0]['p_value']
                pstr = f"KW p={pv:.3f}" if pv >= 0.001 else f"KW p={pv:.2e}"
                ax.text(0.5, 0.97, pstr, transform=ax.transAxes, ha='center', va='top',
                        fontsize=9, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()
    fig.savefig(OUT_FIGS / "A2_promoter_density_boxplot.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("  -> テーブル・図保存完了")

    return stat_df


# ============================================================
# 解析3: 時系列発現変化とプロモーターメチル化
# ============================================================
def analysis3_timeseries(methyl, tss, counts):
    print("\n[解析3] 時系列発現変化とプロモーターメチル化")

    # サンプル列を整理
    T1_cols = [c for c in counts.columns if 'M145_1' in c]
    T2_cols = [c for c in counts.columns if 'M145_2' in c]
    T3_cols = [c for c in counts.columns if 'M145_3' in c]

    counts = counts.set_index('gene_id')
    mean_expr = pd.DataFrame({
        'T1': counts[T1_cols].mean(axis=1),
        'T2': counts[T2_cols].mean(axis=1),
        'T3': counts[T3_cols].mean(axis=1),
    })
    mean_expr = mean_expr[mean_expr.max(axis=1) > 1]  # 低発現除去

    # 単調増加/減少の分類
    def monotonic_class(row):
        if row['T1'] < row['T2'] < row['T3']:
            return 'mono_inc'
        elif row['T1'] > row['T2'] > row['T3']:
            return 'mono_dec'
        else:
            return 'other'

    mean_expr['trend'] = mean_expr.apply(monotonic_class, axis=1)
    trend_counts = mean_expr['trend'].value_counts()
    print(f"  単調増加: {trend_counts.get('mono_inc', 0)}, "
          f"単調減少: {trend_counts.get('mono_dec', 0)}, "
          f"その他: {trend_counts.get('other', 0)}")

    tss_df = tss[['gene_id', 'chrom', 'tss']].dropna(subset=['tss']).copy()
    tss_df['tss'] = tss_df['tss'].astype(int)
    tss_df['trend'] = tss_df['gene_id'].map(mean_expr['trend']).fillna('other')

    # タイムポイント別プロモーター密度
    density_rows = []
    for mod_label in ['4mC', '6mA']:
        for tp in ['T1', 'T2', 'T3']:
            sites = methyl[(methyl['mod_type'] == mod_label) &
                           (methyl['timepoint'] == tp)][['chrom', 'position']].drop_duplicates()

            for chrom, grp in tss_df.groupby('chrom'):
                s = sites[sites['chrom'] == chrom]['position'].values
                for _, row in grp.iterrows():
                    tss_pos = row['tss']
                    lo = tss_pos - PROMOTER_WINDOW
                    hi = tss_pos + PROMOTER_WINDOW
                    n_sites = np.sum((s >= lo) & (s <= hi))
                    density_rows.append({
                        'gene_id': row['gene_id'],
                        'trend': row['trend'],
                        'mod_type': mod_label,
                        'timepoint': tp,
                        'density': n_sites / (2 * PROMOTER_WINDOW / 1000),
                    })

    density_df = pd.DataFrame(density_rows)

    # 群ごとの時系列平均
    summary = density_df.groupby(['trend', 'mod_type', 'timepoint'])['density'].agg(
        ['mean', 'sem', 'count']).reset_index()
    summary.to_csv(OUT_TABLES / "A3_timeseries_methylation.tsv", sep='\t', index=False)

    # 単調群の遺伝子リストを保存
    trend_df = mean_expr[mean_expr['trend'].isin(['mono_inc', 'mono_dec'])][['T1', 'T2', 'T3', 'trend']].reset_index()
    trend_df.to_csv(OUT_TABLES / "A3_monotonic_genes.tsv", sep='\t', index=False)

    # ========== 可視化 ==========
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("Analysis 3: Promoter Methylation Dynamics\n(Monotonically Increasing vs Decreasing Genes)",
                 fontsize=13, fontweight='bold')

    for col_idx, mod_label in enumerate(['4mC', '6mA']):
        ax = axes[col_idx]
        sub = summary[summary['mod_type'] == mod_label]

        for trend, color, ls in [('mono_inc', COLORS['mono_inc'], '-'),
                                   ('mono_dec', COLORS['mono_dec'], '--'),
                                   ('other',    COLORS['other'],    ':')]:
            d = sub[sub['trend'] == trend].sort_values('timepoint')
            if len(d) == 0:
                continue
            n = d['count'].iloc[0]
            label = {'mono_inc': f'Mono. increase (n={n})',
                     'mono_dec': f'Mono. decrease (n={n})',
                     'other':    f'Other (n={n})'}[trend]
            ax.errorbar(d['timepoint'], d['mean'], yerr=d['sem'],
                        color=color, linestyle=ls, marker='o', linewidth=2,
                        capsize=4, label=label)

        ax.set_xlabel("Timepoint")
        ax.set_ylabel("Methylation density [sites/kb]")
        ax.set_title(f"{mod_label} methylation")
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig(OUT_FIGS / "A3_timeseries_methylation.png", dpi=300, bbox_inches='tight')
    plt.close()

    # 発現変化のヒートマップ (mono_inc / mono_dec 遺伝子の例)
    fig2, axes2 = plt.subplots(1, 2, figsize=(10, 5))
    fig2.suptitle("Analysis 3: Expression Trajectory of Monotonic Gene Groups", fontsize=12, fontweight='bold')

    for col_idx, (trend, color) in enumerate([('mono_inc', COLORS['mono_inc']),
                                               ('mono_dec', COLORS['mono_dec'])]):
        ax = axes2[col_idx]
        gene_ids = mean_expr[mean_expr['trend'] == trend].index.tolist()
        if len(gene_ids) == 0:
            ax.text(0.5, 0.5, 'No genes', ha='center', va='center')
            continue

        # 正規化して描画 (最大値で正規化)
        expr_sub = mean_expr.loc[gene_ids, ['T1', 'T2', 'T3']]
        expr_sub_norm = expr_sub.div(expr_sub.max(axis=1) + 1e-6, axis=0)

        # 各遺伝子をlightに、平均を強調
        for _, row in expr_sub_norm.iterrows():
            ax.plot(['T1', 'T2', 'T3'], row.values, color=color, alpha=0.05, linewidth=0.5)
        mean_norm = expr_sub_norm.mean()
        ax.plot(['T1', 'T2', 'T3'], mean_norm.values, color=color, linewidth=3,
                marker='o', label='Mean', zorder=5)

        title = "Monotonically Increasing" if trend == 'mono_inc' else "Monotonically Decreasing"
        ax.set_title(f"{title}\n(n={len(gene_ids)})")
        ax.set_ylabel("Normalized expression")
        ax.set_xlabel("Timepoint")
        ax.legend()
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig2.savefig(OUT_FIGS / "A3_monotonic_expression.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("  -> テーブル・図保存完了")

    return summary


# ============================================================
# サマリーテーブル生成
# ============================================================
def generate_summary_table(stat1, stat2, summary3):
    print("\n[サマリー] 統合サマリーテーブル生成")

    rows = []

    # Analysis 1
    for _, r in stat1.iterrows():
        rows.append({
            'Analysis': 'A1_nearest_dist',
            'Comparison': r['comparison'],
            'ModType': r['mod_type'],
            'Test': r['test'],
            'Statistic': f"U={r['U_stat']:.1f}",
            'P_value': r['p_value'],
            'EffectSize': f"r={r['rank_biserial_r']:.3f}",
            'Group1_median': f"{r['median_g1']:.1f}bp",
            'Group2_median': f"{r['median_g2']:.1f}bp",
            'Interpretation': "NS" if r['p_value'] >= 0.05 else "Significant",
        })

    # Analysis 2
    for _, r in stat2.iterrows():
        rows.append({
            'Analysis': 'A2_promoter_density',
            'Comparison': r['comparison'],
            'ModType': r['mod_type'],
            'Test': r['test'],
            'Statistic': f"H={r['H_stat']:.3f}",
            'P_value': r['p_value'],
            'EffectSize': '',
            'Group1_median': f"up={r['median_up']:.3f}",
            'Group2_median': f"down={r['median_down']:.3f}",
            'Interpretation': "NS" if r['p_value'] >= 0.05 else "Significant",
        })

    df = pd.DataFrame(rows)
    df.to_csv(OUT_TABLES / "spatial_integration_summary.tsv", sep='\t', index=False)
    print(f"  -> {OUT_TABLES / 'spatial_integration_summary.tsv'}")
    return df


# ============================================================
# Main
# ============================================================
def main():
    print("=" * 60)
    print("67_spatial_integration: DEG vs 메チル化サイト空間統合解析")
    print("=" * 60)

    methyl, tss, counts, degs = load_data()

    stat1   = analysis1_nearest_distance(methyl, tss, degs)
    stat2   = analysis2_promoter_density(methyl, tss, degs)
    summary3 = analysis3_timeseries(methyl, tss, counts)

    summary_df = generate_summary_table(stat1, stat2, summary3)

    print("\n" + "=" * 60)
    print("[完了] 全解析終了")
    print(f"  テーブル: {OUT_TABLES}")
    print(f"  図:       {OUT_FIGS}")
    print("=" * 60)

    print("\n--- Summary Table ---")
    print(summary_df[['Analysis', 'Comparison', 'ModType', 'P_value', 'Interpretation']].to_string(index=False))


if __name__ == "__main__":
    main()
