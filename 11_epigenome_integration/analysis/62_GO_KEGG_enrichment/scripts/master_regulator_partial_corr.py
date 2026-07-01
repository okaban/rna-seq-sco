"""
B-5: マスターレギュレーター偏相関解析
=========================================
メチル化近傍遺伝子の発現変化がBldA/BldD等のマスターレギュレーターで
説明されるかを偏相関で検証する。

マスターレギュレーター:
  BldA (SCO4009) = SC_RS22155  hybrid sensor kinase/response regulator
  BldD (SCO3778) = SC_RS21010  threonine--tRNA ligase (old annotation: bldD)
  BldM (SCO3063) = SC_RS17430  response regulator TF
  WhiA (SCO3579) = SC_RS20030  wblA TF
  WhiB (SCO3034) = SC_RS17285  WhiB family TF

実施日: 2026-04-18
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pingouin as pg
from scipy import stats
from statsmodels.stats.multitest import multipletests
from pathlib import Path

# === パス設定 ===
ANALYSIS_DIR = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/62_GO_KEGG_enrichment")
DATA_DIR     = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis")
GENE_ANNOT   = Path("/Users/okaban/bioinfo/rna-seq/05_annotation/analysis/05_annotation_260128_v1/tables/gene_annotation_basic.tsv")
TSS_TABLE    = DATA_DIR / "18_tss_analyses/comprehensive_tss_table.csv"
METH_SITES   = DATA_DIR / "07_motif_analysis/methylation_site_sequences.csv"
DESEQ_T2T1   = Path("/Users/okaban/bioinfo/rna-seq/04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_2_vs_1.tsv")
OUT_TABLES   = ANALYSIS_DIR / "tables"
OUT_FIGS     = ANALYSIS_DIR / "figures"

# マスターレギュレーター gene_id マッピング（old_locus_tag → gene_id）
MASTER_REGS = {
    'BldA (SCO4009)': 'SC_RS22155',
    'BldD (SCO3778)': 'SC_RS21010',
    'BldM (SCO3063)': 'SC_RS17430',
    'WhiA (SCO3579)': 'SC_RS20030',
    'WhiB (SCO3034)': 'SC_RS17285',
}

WINDOW_SIZE = 2000  # ±2kb
NORM_COUNTS = Path("/Users/okaban/bioinfo/rna-seq/04_deseq2/analysis/04_deseq2_260128_v1/results/normalized_counts_M145.tsv")
# サンプル列名: M145_1_1/2/3 = T1, M145_2_1/3/4 = T2, M145_3_2/3/4 = T3
T1_COLS = ['M145_1_1', 'M145_1_2', 'M145_1_3']
T2_COLS = ['M145_2_1', 'M145_2_3', 'M145_2_4']


def load_data():
    annot  = pd.read_csv(GENE_ANNOT, sep='\t')
    tss_df = pd.read_csv(TSS_TABLE)
    meth   = pd.read_csv(METH_SITES)
    deseq  = pd.read_csv(DESEQ_T2T1, sep='\t')
    counts = pd.read_csv(NORM_COUNTS, sep='\t')
    return annot, tss_df, meth, deseq, counts


def get_gccggc_t1_sites(meth: pd.DataFrame) -> pd.DataFrame:
    """T1 GCCGGCモチーフ（4mC）サイト座標を抽出"""
    t1_4mc = meth[(meth['timepoint'] == 'T1') & (meth['mod_type'] == '4mC')].copy()
    # GCCGGC または GGCCGG（reverse complement）を含むもの
    mask = (t1_4mc['sequence'].str.contains('GCCGGC', regex=False) |
            t1_4mc['sequence'].str.contains('GGCCGG', regex=False))
    sites = t1_4mc[mask][['chrom', 'position']].copy()
    print(f"  GCCGGC T1 4mC sites: {len(sites)}")
    return sites


def assign_proximity(annot: pd.DataFrame, tss_df: pd.DataFrame,
                     sites: pd.DataFrame, window: int = 2000) -> pd.Series:
    """
    各遺伝子にGCCGGCサイト±2kb内に存在するか（1/0）を割り当て
    遺伝子のTSS座標を基準にする
    """
    # TSS + start/end を結合
    gene_coords = tss_df[['gene_id', 'chrom', 'start', 'end', 'tss']].copy()

    site_positions = sites['position'].values

    # ±2kb ウィンドウ内サイト数を計算
    proximal = {}
    for _, gene in gene_coords.iterrows():
        tss = float(gene['tss'])
        lo  = tss - window
        hi  = tss + window
        n_sites = ((site_positions >= lo) & (site_positions <= hi)).sum()
        proximal[gene['gene_id']] = int(n_sites > 0)

    return pd.Series(proximal, name='gccggc_proximal')


def residualize(y: np.ndarray, covariates: np.ndarray) -> np.ndarray:
    """
    OLSでy ~ covariatesを回帰し残差を返す
    （偏相関の手動計算: partial_corr(X,Y|Z) = pearsonr(resid(X|Z), resid(Y|Z))）
    """
    if covariates.ndim == 1:
        covariates = covariates[:, np.newaxis]
    # intercept追加
    X = np.column_stack([np.ones(len(y)), covariates])
    # OLS: beta = (X'X)^-1 X'y
    try:
        beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
        return y - X @ beta
    except Exception:
        return y


def partial_corr_manual(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple:
    """
    手動偏相関: pearsonr(resid(x|z), resid(y|z))
    Returns (r, p)
    """
    resid_x = residualize(x.astype(float), z.astype(float))
    resid_y = residualize(y.astype(float), z.astype(float))
    return stats.pearsonr(resid_x, resid_y)


def run_incremental_partial_corr(df: pd.DataFrame, target: str, covariate: str,
                                 control_vars: list) -> pd.DataFrame:
    """
    コントロール変数を1つずつ追加しながら偏相関を計算
    OLS残差ベースの手動偏相関を使用（rank deficiency 問題を回避）
    """
    rows = []
    sub_all = df[[target, covariate] + control_vars].dropna()

    x_all = sub_all[target].values.astype(float)
    y_all = sub_all[covariate].values.astype(float)
    z_all = sub_all[control_vars].values.astype(float)
    n_all = len(sub_all)

    # 単純相関
    r0, p0 = stats.pearsonr(x_all, y_all)
    rows.append({'controlled_for': '(none)', 'r': round(r0, 4), 'p': round(p0, 8), 'n': n_all})

    # 1変数ずつ
    for i, var in enumerate(control_vars):
        sub = df[[target, covariate, var]].dropna()
        x = sub[target].values.astype(float)
        y = sub[covariate].values.astype(float)
        z = sub[var].values.astype(float)
        try:
            r, p = partial_corr_manual(x, y, z)
            rows.append({'controlled_for': var, 'r': round(float(r), 4), 'p': round(float(p), 8), 'n': len(sub)})
        except Exception as e:
            rows.append({'controlled_for': var, 'r': np.nan, 'p': np.nan, 'n': len(sub)})

    # 全変数まとめて
    try:
        r, p = partial_corr_manual(x_all, y_all, z_all)
        rows.append({'controlled_for': 'all_MRs', 'r': round(float(r), 4), 'p': round(float(p), 8), 'n': n_all})
    except Exception as e:
        rows.append({'controlled_for': 'all_MRs', 'r': np.nan, 'p': np.nan, 'n': n_all})

    return pd.DataFrame(rows)


def plot_comparison(result_df: pd.DataFrame, out_path: Path):
    """単純相関 vs 偏相関の比較プロット"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))

    # r の変化
    labels = result_df['controlled_for'].tolist()
    r_vals = result_df['r'].tolist()
    colors = ['#2196F3' if i == 0 else '#FF5722' if i == len(labels)-1 else '#9E9E9E'
              for i in range(len(labels))]

    bars = ax1.barh(range(len(labels)), r_vals, color=colors, edgecolor='white', height=0.6)
    ax1.set_yticks(range(len(labels)))
    ax1.set_yticklabels(labels, fontsize=9)
    ax1.set_xlabel('Pearson r', fontsize=10)
    ax1.set_title('Correlation: GCCGGC proximity ~ log2FC(T2/T1)', fontsize=10)
    ax1.axvline(0, color='black', linewidth=0.8)
    ax1.axvline(r_vals[0], color='#2196F3', linewidth=1.2, linestyle='--', alpha=0.6)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # -log10(p) の変化
    p_vals = result_df['p'].values
    p_vals = np.where(p_vals == 0, 1e-300, p_vals)
    log_p = -np.log10(p_vals)
    ax2.barh(range(len(labels)), log_p, color=colors, edgecolor='white', height=0.6)
    ax2.set_yticks(range(len(labels)))
    ax2.set_yticklabels(labels, fontsize=9)
    ax2.set_xlabel('-log10(p-value)', fontsize=10)
    ax2.set_title('Significance', fontsize=10)
    ax2.axvline(-np.log10(0.05), color='red', linewidth=1, linestyle='--', alpha=0.7, label='p=0.05')
    ax2.legend(fontsize=8)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(str(out_path), format='pdf', bbox_inches='tight')
    plt.savefig(str(out_path).replace('.pdf', '.svg'), format='svg', bbox_inches='tight')
    plt.savefig(str(out_path).replace('.pdf', '.png'), format='png', bbox_inches='tight', dpi=150)
    plt.close()
    print(f"  図保存: {out_path.name}")


# ============================================================
# メイン
# ============================================================
if __name__ == '__main__':
    print("=== B-5: マスターレギュレーター偏相関解析 ===\n")

    # 1. データ読み込み
    print("1. データ読み込み...")
    annot, tss_df, meth, deseq, counts = load_data()
    print(f"   DESeq2 T2vsT1: {len(deseq)} genes")
    print(f"   正規化カウント: {len(counts)} genes × {counts.shape[1]-1} samples")

    # 2. マスターレギュレーター発現変化の確認
    print("\n2. マスターレギュレーター T1→T2 log2FC:")
    mr_log2fc = {}
    mr_display = []
    for name, gid in MASTER_REGS.items():
        row = deseq[deseq['gene_id'] == gid]
        if len(row) > 0:
            lfc = float(row['log2FoldChange'].iloc[0])
            pval = row['padj'].iloc[0]
            mr_log2fc[name] = lfc
            mr_display.append({'regulator': name, 'gene_id': gid,
                                'log2FC': round(lfc, 3), 'padj': pval})
            print(f"   {name}: log2FC={lfc:.3f}, padj={pval:.2e}")
        else:
            print(f"   {name} ({gid}): DESeq2結果なし")

    # 3. GCCGGC近傍遺伝子の割り当て
    print(f"\n3. GCCGGC±{WINDOW_SIZE}bp近傍遺伝子の割り当て...")
    sites = get_gccggc_t1_sites(meth)
    proximity = assign_proximity(annot, tss_df, sites, window=WINDOW_SIZE)
    n_proximal = proximity.sum()
    n_total = len(proximity)
    print(f"   近傍遺伝子数: {n_proximal} / {n_total} ({100*n_proximal/n_total:.1f}%)")

    # 4. 解析用データフレーム構築
    print("\n4. 解析データフレーム構築...")
    df = deseq[['gene_id', 'log2FoldChange']].copy()
    df.columns = ['gene_id', 'log2FC_T2T1']
    df['log2FC_T2T1'] = pd.to_numeric(df['log2FC_T2T1'], errors='coerce')
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=['log2FC_T2T1'])

    # 近傍フラグを追加
    df['gccggc_proximal'] = df['gene_id'].map(proximity).fillna(0).astype(int)

    # --- 共発現スコア: 各遺伝子とMRのPearson r（全9サンプル）---
    # 各遺伝子の発現プロファイル（log1p変換して相関）
    all_sample_cols = T1_COLS + T2_COLS
    # T2-T1の追加サンプル列を確認
    available_cols = [c for c in all_sample_cols if c in counts.columns]
    print(f"   使用サンプル列: {available_cols}")

    counts_sub = counts.set_index('gene_id')[available_cols].copy()
    counts_log = np.log1p(counts_sub)

    # MR共発現スコアを各遺伝子について計算
    mr_cols = []
    for name, gid in MASTER_REGS.items():
        col = f"coexpr_{name.split('(')[0].strip().replace(' ', '_')}"
        mr_cols.append(col)
        if gid in counts_log.index:
            mr_expr = counts_log.loc[gid].values
            def corr_with_mr(gene_expr):
                if gene_expr.std() < 1e-10 or mr_expr.std() < 1e-10:
                    return 0.0
                r, _ = stats.pearsonr(gene_expr, mr_expr)
                return r if not np.isnan(r) else 0.0
            coexpr_vals = counts_log.apply(corr_with_mr, axis=1)
            df[col] = df['gene_id'].map(coexpr_vals).fillna(0)
            print(f"   {name} 共発現スコア: mean={df[col].mean():.3f}, std={df[col].std():.3f}")
        else:
            df[col] = 0.0
            print(f"   {name} ({gid}): カウントデータなし")

    print(f"\n   解析対象遺伝子数: {len(df)}")
    print(f"   GCCGGC近傍: {df['gccggc_proximal'].sum()} 遺伝子")

    # 5. 単純相関 vs 偏相関（コントロール変数を順次追加）
    print("\n5. 偏相関解析...")
    result_df = run_incremental_partial_corr(
        df, target='gccggc_proximal', covariate='log2FC_T2T1', control_vars=mr_cols
    )
    result_df.columns = ['controlled_for', 'r', 'p', 'n']
    result_df['delta_r'] = result_df['r'] - result_df.loc[0, 'r']
    result_df['-log10p'] = -np.log10(result_df['p'].clip(1e-300))

    print("\n  結果:")
    print(result_df.to_string(index=False))

    # 6. 追加解析: ポイントバイシリアル相関（連続×2値）
    print("\n6. ポイントバイシリアル相関 (point-biserial)...")
    proximal_fc   = df[df['gccggc_proximal'] == 1]['log2FC_T2T1'].dropna()
    distal_fc     = df[df['gccggc_proximal'] == 0]['log2FC_T2T1'].dropna()
    r_pb, p_pb    = stats.pointbiserialr(df['gccggc_proximal'], df['log2FC_T2T1'].fillna(0))
    t_stat, p_t   = stats.ttest_ind(proximal_fc, distal_fc)
    print(f"   近傍遺伝子 (n={len(proximal_fc)}): mean log2FC = {proximal_fc.mean():.3f} ± {proximal_fc.std():.3f}")
    print(f"   遠位遺伝子 (n={len(distal_fc)}): mean log2FC = {distal_fc.mean():.3f} ± {distal_fc.std():.3f}")
    print(f"   Point-biserial r = {r_pb:.4f}, p = {p_pb:.4e}")
    print(f"   t-test: t={t_stat:.3f}, p={p_t:.4e}")

    # 7. 保存
    out_pcorr = OUT_TABLES / "B5_MR_partial_correlation.tsv"
    result_df.to_csv(out_pcorr, sep='\t', index=False)
    print(f"\n  保存: {out_pcorr.name}")

    # MR summary
    mr_df = pd.DataFrame(mr_display)
    mr_out = OUT_TABLES / "B5_master_regulator_log2FC.tsv"
    mr_df.to_csv(mr_out, sep='\t', index=False)
    print(f"  保存: {mr_out.name}")

    # 8. 可視化
    print("\n7. 可視化...")
    plot_comparison(result_df, OUT_FIGS / "B5_partial_correlation.pdf")

    # 9. サマリー
    r_simple  = float(result_df.loc[0, 'r'])
    p_simple  = float(result_df.loc[0, 'p'])
    r_partial = float(result_df[result_df['controlled_for'] == 'all_MRs']['r'].iloc[0])
    p_partial = float(result_df[result_df['controlled_for'] == 'all_MRs']['p'].iloc[0])

    print("\n" + "="*60)
    print("解析結果サマリー")
    print("="*60)
    print(f"  単純相関 (GCCGGC proximity ~ log2FC T2vsT1):")
    print(f"    r = {r_simple:.4f}, p = {p_simple:.4e}")
    print(f"  偏相関 (全MRコントロール後):")
    print(f"    r = {r_partial:.4f}, p = {p_partial:.4e}")
    delta = r_partial - r_simple
    pct   = 100 * abs(delta) / max(abs(r_simple), 1e-9)
    print(f"  Δr = {delta:+.4f} ({pct:.1f}% 変化)")
    if abs(r_partial) >= abs(r_simple) * 0.7:
        print(f"  解釈: MRコントロール後もrが維持 → メチル化の独立した寄与あり")
    else:
        print(f"  解釈: MRコントロール後にrが大きく低下 → MR依存的な発現変化")

    print("\n  各MRコントロール後のr:")
    for _, row in result_df.iterrows():
        print(f"    {row['controlled_for']:30s}: r={row['r']:.4f}, p={row['p']:.4e}")

    print("\n=== 完了 ===")
