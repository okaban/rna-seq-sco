"""
70_methylation_dynamics: A-5 メチル化レート dynamics T1/T2/T3
=============================================================
GCCGGC 4mC と AAGCCCG 6mA の per-site 時系列変動（dynamic vs stable）分類

実施日: 2026-04-20
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from scipy import stats
from pathlib import Path
import re

# === パス設定 ===
ANALYSIS_DIR = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/70_methylation_dynamics")
DATA_PATH = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/01_integration/high_confidence_sites_weighted.csv")
GENOME_PATH = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/data/NC_003888.3.fna")

OUT_TABLES = ANALYSIS_DIR / "tables"
OUT_FIGS   = ANALYSIS_DIR / "figures"

# 染色体 core/arm 境界 (プロジェクト既定値)
ARM_LEFT_END     = 1_500_000
ARM_RIGHT_START  = 7_167_508

# Dynamic 分類閾値
DELTA_THRESHOLDS = [0.2, 0.3, 0.4]
MAIN_THRESHOLD   = 0.3

# カラーパレット
# Unified palette "Rich & Calm" (2026-07): timepoints use the grayscale ramp,
# categories use muted tones threaded from the shared palette.
COLORS = {
    'dynamic': '#A64B44',   # unified rich muted red (changing)
    'stable':  '#9AA7B0',   # unified neutral grey (stable)
    'T1':      '#AAB3BB',   # grayscale ramp — early / light
    'T2':      '#6B7783',   # grayscale ramp — mid
    'T3':      '#293039',   # grayscale ramp — late / dark
}

PATTERN_COLORS = {
    'Monotone increase': '#3E7256',   # deep muted green
    'Monotone decrease': '#A64B44',   # rich muted red
    'Peak (T2 max)':     '#C0803A',   # calm amber
    'Valley (T2 min)':   '#6E5495',   # muted purple
    'Other':             '#9AA7B0',   # neutral grey
}


# ============================================================
# 1. ゲノム読み込み & モチーフ割り当て
# ============================================================

def load_genome(fasta_path: Path) -> str:
    """NC_003888.3 の配列を文字列として返す (1-based → 0-based index)"""
    seq_parts = []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith('>'):
                continue
            seq_parts.append(line.strip().upper())
    return ''.join(seq_parts)


def assign_motif(genome: str, chrom: str, position: int, strand: str, mod_type: str,
                 flank: int = 15) -> str:
    """
    参照配列からモチーフを判定。
    Returns: 'GCCGGC', 'AAGCCCG', or 'other'
    """
    pos0 = position - 1  # 0-based
    start = max(0, pos0 - flank)
    end   = min(len(genome), pos0 + flank + 1)
    ctx   = genome[start:end]

    # サイト中心位置を文脈内相対座標で確認
    center_in_ctx = pos0 - start

    if mod_type == '4mC':
        # GCCGGC: センターC が motif の index 1 or 4 (0-based) にある
        # + 鎖: G-C-C-G-G-C  → メチル化C は index 1 (2番目) or index 4 (内部GGCのC)
        # 実際の PacBio 解析では内部 C (index 1, 2) が対象
        # シンプルに: contextにGCCGGCがあり、中心位置がGCCGGC内のCかどうか確認
        for m in re.finditer('GCCGGC', ctx):
            motif_start = m.start()
            # GCCGGC の C は index 1, 2 (GCCGGC: G=0,C=1,C=2,G=3,G=4,C=5)
            c_positions_in_motif = [motif_start + 1, motif_start + 2, motif_start + 5]
            if center_in_ctx in c_positions_in_motif:
                return 'GCCGGC'
        # reverse complement: GCCGGC は palindromic
        for m in re.finditer('GCCGGC', ctx):
            motif_start = m.start()
            c_positions_in_motif = [motif_start + 1, motif_start + 2, motif_start + 5]
            if center_in_ctx in c_positions_in_motif:
                return 'GCCGGC'
        # AAGCCCG 内の C: AAGCCCG → A=0,A=1,G=2,C=3,C=4,C=5,G=6
        # 4mC might also occur in AAGCCCG context
        for m in re.finditer('AAGCCCG', ctx):
            motif_start = m.start()
            c_positions = [motif_start + 3, motif_start + 4, motif_start + 5]
            if center_in_ctx in c_positions:
                return 'AAGCCCG'
        return 'other'

    elif mod_type == '6mA':
        # AAGCCCG: A が index 0 or 1
        for m in re.finditer('AAGCCCG', ctx):
            motif_start = m.start()
            a_positions = [motif_start, motif_start + 1]
            if center_in_ctx in a_positions:
                return 'AAGCCCG'
        # reverse complement of AAGCCCG = CGGGCTT → 6mA on A
        for m in re.finditer('CGGGCTT', ctx):
            motif_start = m.start()
            # A on reverse: CGGGCTT の T(6) は元鎖の A
            t_positions = [motif_start + 5, motif_start + 6]
            if center_in_ctx in t_positions:
                return 'AAGCCCG'
        return 'other'

    return 'other'


def assign_motifs_fast(df: pd.DataFrame, genome: str) -> pd.Series:
    """
    ベクトル化実装: 各サイトのモチーフを返す。
    ゲノム文字列から ±15bp コンテキストを取り出し、正規表現でモチーフ確認。
    """
    motifs = []
    genome_len = len(genome)

    for _, row in df.iterrows():
        pos0   = int(row['position']) - 1
        mod    = row['mod_type']
        strand = row['strand']
        flank  = 15
        start  = max(0, pos0 - flank)
        end    = min(genome_len, pos0 + flank + 1)
        ctx    = genome[start:end]
        center = pos0 - start

        motif = 'other'

        if mod == '4mC':
            for m in re.finditer('GCCGGC', ctx):
                ms = m.start()
                if center in (ms + 1, ms + 2, ms + 5):
                    motif = 'GCCGGC'
                    break
            if motif == 'other':
                for m in re.finditer('AAGCCCG', ctx):
                    ms = m.start()
                    if center in (ms + 3, ms + 4, ms + 5):
                        motif = 'AAGCCCG'
                        break

        elif mod == '6mA':
            for m in re.finditer('AAGCCCG', ctx):
                ms = m.start()
                if center in (ms, ms + 1):
                    motif = 'AAGCCCG'
                    break
            if motif == 'other':
                for m in re.finditer('CGGGCTT', ctx):
                    ms = m.start()
                    if center in (ms + 5, ms + 6):
                        motif = 'AAGCCCG'
                        break

        motifs.append(motif)

    return pd.Series(motifs, index=df.index)


# ============================================================
# 2. データ読み込み & Wide テーブル作成
# ============================================================

def build_wide_table(df: pd.DataFrame, mod_type: str) -> pd.DataFrame:
    """
    指定 mod_type のサイトについて T1/T2/T3 を横持ちに変換。
    3 timepoint すべてにデータがあるサイトのみ残す。
    """
    sub = df[df['mod_type'] == mod_type].copy()
    wide = sub.pivot_table(
        index=['chrom', 'position', 'strand', 'motif'],
        columns='timepoint',
        values='weighted_mod_freq',
        aggfunc='first',
    ).reset_index()
    wide.columns.name = None
    # 3 timepoint 完全なサイトのみ
    wide = wide.dropna(subset=['T1', 'T2', 'T3'])
    wide = wide.rename(columns={'T1': 'T1_freq', 'T2': 'T2_freq', 'T3': 'T3_freq'})
    return wide


# ============================================================
# 3. Dynamic/Stable 分類
# ============================================================

def classify_dynamics(wide: pd.DataFrame, threshold: float) -> pd.DataFrame:
    df = wide.copy()
    df['delta'] = df[['T1_freq', 'T2_freq', 'T3_freq']].max(axis=1) - \
                  df[['T1_freq', 'T2_freq', 'T3_freq']].min(axis=1)
    df['category'] = np.where(df['delta'] > threshold, 'dynamic', 'stable')
    return df


def classify_pattern(row) -> str:
    """Dynamic sites のパターン分類"""
    t1, t2, t3 = row['T1_freq'], row['T2_freq'], row['T3_freq']
    if t1 < t2 < t3:
        return 'Monotone increase'
    elif t1 > t2 > t3:
        return 'Monotone decrease'
    elif t2 >= t1 and t2 >= t3 and t2 > max(t1, t3):
        return 'Peak (T2 max)'
    elif t2 <= t1 and t2 <= t3 and t2 < min(t1, t3):
        return 'Valley (T2 min)'
    else:
        return 'Other'


def assign_region(position: float) -> str:
    if position <= ARM_LEFT_END or position >= ARM_RIGHT_START:
        return 'arm'
    return 'core'


# ============================================================
# 4. 統計検定
# ============================================================

def core_arm_chisq(df: pd.DataFrame, motif_name: str) -> dict:
    """Dynamic vs Stable × Core vs Arm のカイ二乗検定"""
    df = df.copy()
    df['region'] = df['position'].apply(assign_region)

    ct = pd.crosstab(df['category'], df['region'])
    # ensure both categories present
    for cat in ('dynamic', 'stable'):
        if cat not in ct.index:
            ct.loc[cat] = 0
    for reg in ('core', 'arm'):
        if reg not in ct.columns:
            ct[reg] = 0
    ct = ct.loc[['dynamic', 'stable'], ['core', 'arm']]

    # Use Fisher's exact test when expected frequencies are too small
    try:
        chi2, p, dof, expected = stats.chi2_contingency(ct)
        if (expected < 5).any():
            raise ValueError("Small expected frequencies")
    except ValueError:
        # Fall back to Fisher's exact test (2x2 table)
        odds, p = stats.fisher_exact(ct.values)
        chi2 = float('nan')
        dof = 1
    total = len(df)
    n_dyn = (df['category'] == 'dynamic').sum()
    n_stable = (df['category'] == 'stable').sum()

    dyn_core_pct = (ct.loc['dynamic', 'core'] / n_dyn * 100) if n_dyn > 0 else 0
    dyn_arm_pct  = (ct.loc['dynamic', 'arm']  / n_dyn * 100) if n_dyn > 0 else 0
    stab_core_pct = (ct.loc['stable', 'core'] / n_stable * 100) if n_stable > 0 else 0
    stab_arm_pct  = (ct.loc['stable', 'arm']  / n_stable * 100) if n_stable > 0 else 0

    return {
        'motif': motif_name,
        'n_total': total,
        'n_dynamic': n_dyn,
        'n_stable':  n_stable,
        'dynamic_pct': n_dyn / total * 100,
        'dynamic_core': ct.loc['dynamic', 'core'],
        'dynamic_arm':  ct.loc['dynamic', 'arm'],
        'dynamic_core_pct': dyn_core_pct,
        'dynamic_arm_pct':  dyn_arm_pct,
        'stable_core': ct.loc['stable', 'core'],
        'stable_arm':  ct.loc['stable', 'arm'],
        'stable_core_pct': stab_core_pct,
        'stable_arm_pct':  stab_arm_pct,
        'chi2': chi2,
        'p_value': p,
        'dof': dof,
        'crosstab': ct,
    }


# ============================================================
# 5. 可視化
# ============================================================

def plot_violin(gccggc_df: pd.DataFrame, aagcccg_df: pd.DataFrame):
    """
    Figure 1: T1/T2/T3 の per-site メチル化頻度分布 (violin) +
    delta distribution (histogram)
    全サイトが nearly-dynamic なので dynamic/stable 分類より
    T1→T2→T3 の時系列変化を可視化する
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Methylation frequency dynamics: per-site trajectories',
                 fontsize=13, fontweight='bold')

    tp_cols   = ['T1_freq', 'T2_freq', 'T3_freq']
    tp_labels = ['T1', 'T2', 'T3']

    # Upper row: violin per timepoint
    for ax, df, motif_label in [
        (axes[0, 0], gccggc_df, 'GCCGGC 4mC\n(T1→T2→T3 methylation frequency)'),
        (axes[0, 1], aagcccg_df, 'AAGCCCG 6mA\n(T1→T2→T3 methylation frequency)'),
    ]:
        data_by_tp = [df[c].dropna().values for c in tp_cols]
        # skip empty
        valid = [(d, tp) for d, tp in zip(data_by_tp, tp_labels) if len(d) > 0]

        positions = list(range(1, len(valid) + 1))
        vp = ax.violinplot([d for d, _ in valid], positions=positions,
                            showmedians=True, showextrema=True, widths=0.6)

        for i, (body, (d, tp)) in enumerate(zip(vp['bodies'], valid)):
            body.set_facecolor(COLORS[tp])
            body.set_alpha(0.7)
        vp['cmedians'].set_color('black')
        vp['cmedians'].set_linewidth(2)

        ax.set_xticks(positions)
        ax.set_xticklabels([tp for _, tp in valid])
        ax.set_ylabel('Methylation frequency (%)', fontsize=10)
        ax.set_title(motif_label, fontsize=10)
        ax.set_ylim(0, 105)
        n_all = len(df)
        ax.text(0.98, 0.98, f'n={n_all} sites\n(all 3 TPs present)',
                transform=ax.transAxes, ha='right', va='top', fontsize=8,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # Lower row: delta distribution histogram with threshold lines
    for ax, df, motif_label in [
        (axes[1, 0], gccggc_df, 'GCCGGC 4mC — delta distribution'),
        (axes[1, 1], aagcccg_df, 'AAGCCCG 6mA — delta distribution'),
    ]:
        deltas = df['delta'].values
        ax.hist(deltas, bins=40, color='#6B7783', alpha=0.8, edgecolor='white')
        for thr, col, ls in [(0.2, '#C0803A', '--'), (0.3, '#A64B44', '-'), (0.4, '#6E5495', ':')]:
            n_above = (deltas > thr).sum()
            pct = n_above / len(deltas) * 100
            ax.axvline(thr, color=col, linewidth=1.5, linestyle=ls,
                       label=f'>{thr}: {n_above} ({pct:.0f}%)')
        ax.set_xlabel('Delta (max_freq − min_freq) %', fontsize=10)
        ax.set_ylabel('Number of sites', fontsize=10)
        ax.set_title(motif_label, fontsize=10)
        ax.legend(fontsize=8, title='Dynamic threshold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.tight_layout()
    out = OUT_FIGS / 'A5_dynamics_violin.png'
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  Saved: {out}')


def plot_chromosomal_distribution(gccggc_df: pd.DataFrame, aagcccg_df: pd.DataFrame,
                                   genome_len: int = 8_667_507):
    """
    Figure 2: 染色体上のメチル化サイト分布密度
    dynamic=赤, stable=灰色, arm/core 境界を縦線で示す
    """
    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
    fig.suptitle('Chromosomal distribution of methylation sites\n(dynamic=red, stable=grey)',
                 fontsize=12, fontweight='bold')

    bin_size = 100_000  # 100 kb bins
    bins = np.arange(0, genome_len + bin_size, bin_size)

    for ax, df, motif_label in [
        (axes[0], gccggc_df, 'GCCGGC (4mC)'),
        (axes[1], aagcccg_df, 'AAGCCCG (6mA)'),
    ]:
        for cat, color, alpha in [('stable', COLORS['stable'], 0.7),
                                   ('dynamic', COLORS['dynamic'], 0.9)]:
            sub = df[df['category'] == cat]
            ax.hist(sub['position'], bins=bins, color=color, alpha=alpha,
                    label=f'{cat} (n={len(sub)})')

        # arm/core 境界
        for boundary in [ARM_LEFT_END, ARM_RIGHT_START]:
            ax.axvline(boundary, color='navy', linewidth=1.2, linestyle='--', alpha=0.7)

        ax.text(ARM_LEFT_END / 2, ax.get_ylim()[1] * 0.9, 'arm',
                ha='center', fontsize=9, color='navy')
        ax.text((ARM_LEFT_END + ARM_RIGHT_START) / 2, ax.get_ylim()[1] * 0.9, 'core',
                ha='center', fontsize=9, color='navy')
        ax.text((ARM_RIGHT_START + genome_len) / 2, ax.get_ylim()[1] * 0.9, 'arm',
                ha='center', fontsize=9, color='navy')

        ax.set_ylabel('Site count per 100 kb', fontsize=9)
        ax.set_title(motif_label, fontsize=10)
        ax.legend(fontsize=8)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    axes[1].set_xlabel('Chromosomal position (Mb)', fontsize=10)
    axes[1].xaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda x, _: f'{x/1e6:.1f}')
    )
    plt.tight_layout()
    out = OUT_FIGS / 'A5_chromosomal_distribution.png'
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  Saved: {out}')


def plot_pattern_classification(gccggc_dyn: pd.DataFrame, aagcccg_dyn: pd.DataFrame):
    """
    Figure 3: Dynamic サイトのパターン分類 (bar chart + example trajectories)
    """
    fig = plt.figure(figsize=(14, 8))
    gs = GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.35)

    ax_bar_gccggc  = fig.add_subplot(gs[0, 0])
    ax_bar_aagcccg = fig.add_subplot(gs[1, 0])
    ax_traj1 = fig.add_subplot(gs[0, 1])
    ax_traj2 = fig.add_subplot(gs[0, 2])
    ax_traj3 = fig.add_subplot(gs[1, 1])
    ax_traj4 = fig.add_subplot(gs[1, 2])

    pattern_order = ['Monotone increase', 'Monotone decrease',
                     'Peak (T2 max)', 'Valley (T2 min)', 'Other']

    for ax, df, motif_label in [
        (ax_bar_gccggc, gccggc_dyn, 'GCCGGC (4mC)'),
        (ax_bar_aagcccg, aagcccg_dyn, 'AAGCCCG (6mA)'),
    ]:
        counts = df['pattern'].value_counts()
        y_vals = [counts.get(p, 0) for p in pattern_order]
        colors  = [PATTERN_COLORS[p] for p in pattern_order]
        bars = ax.barh(pattern_order, y_vals, color=colors, edgecolor='white')
        for bar, v in zip(bars, y_vals):
            pct = v / sum(y_vals) * 100 if sum(y_vals) > 0 else 0
            ax.text(bar.get_width() + 0.3, bar.get_y() + bar.get_height() / 2,
                    f'{v} ({pct:.0f}%)', va='center', fontsize=8)
        ax.set_title(motif_label, fontsize=10)
        ax.set_xlabel('Count', fontsize=9)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # Example trajectories (1 representative site per pattern × 2 motifs)
    tp_labels = ['T1', 'T2', 'T3']
    tp_cols   = ['T1_freq', 'T2_freq', 'T3_freq']
    traj_axes = [ax_traj1, ax_traj2, ax_traj3, ax_traj4]

    patterns_to_show = ['Monotone increase', 'Monotone decrease',
                        'Peak (T2 max)', 'Valley (T2 min)']

    for ax, pat in zip(traj_axes, patterns_to_show):
        # GCCGGC example
        g_sub = gccggc_dyn[gccggc_dyn['pattern'] == pat]
        a_sub = aagcccg_dyn[aagcccg_dyn['pattern'] == pat]

        plotted = False
        for df_sub, label, linestyle in [(g_sub, 'GCCGGC', '-'),
                                          (a_sub, 'AAGCCCG', '--')]:
            if len(df_sub) > 0:
                # pick site with largest delta as representative
                rep = df_sub.nlargest(1, 'delta').iloc[0]
                vals = [rep[c] for c in tp_cols]
                ax.plot(tp_labels, vals, marker='o', linestyle=linestyle,
                        color=PATTERN_COLORS[pat], label=label, linewidth=1.8)
                plotted = True

        ax.set_title(pat, fontsize=9, color=PATTERN_COLORS[pat], fontweight='bold')
        ax.set_ylabel('Methylation freq (%)', fontsize=8)
        ax.set_ylim(0, 105)
        ax.legend(fontsize=7)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if not plotted:
            ax.text(0.5, 0.5, 'No data', transform=ax.transAxes,
                    ha='center', va='center', fontsize=10, color='grey')

    fig.suptitle('Dynamic site pattern classification', fontsize=12, fontweight='bold')
    out = OUT_FIGS / 'A5_pattern_classification.png'
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  Saved: {out}')


# ============================================================
# 6. サマリーテーブル
# ============================================================

def build_summary(results: dict, thresholds: list) -> pd.DataFrame:
    rows = []
    for thr in thresholds:
        for motif in ('GCCGGC', 'AAGCCCG'):
            df = results[thr][motif]
            n_dyn   = (df['category'] == 'dynamic').sum()
            n_stable = (df['category'] == 'stable').sum()
            n_total  = len(df)
            rows.append({
                'threshold': thr,
                'motif':     motif,
                'n_total':   n_total,
                'n_dynamic': n_dyn,
                'n_stable':  n_stable,
                'pct_dynamic': round(n_dyn / n_total * 100, 1) if n_total else 0,
            })
    return pd.DataFrame(rows)


# ============================================================
# Main
# ============================================================

def main():
    print("=== A-5 Methylation Dynamics Analysis ===\n")

    # 1. Load genome
    print("Loading reference genome...")
    genome = load_genome(GENOME_PATH)
    print(f"  Genome length: {len(genome):,} bp")

    # 2. Load methylation data
    print("Loading methylation data...")
    df = pd.read_csv(DATA_PATH)
    print(f"  Total records: {len(df):,}")

    # 3. Assign motifs
    print("Assigning motifs from sequence context...")
    # Use unique site positions (not per-timepoint) to speed up
    unique_sites = df[['chrom', 'position', 'strand', 'mod_type']].drop_duplicates()
    unique_sites = unique_sites.copy()
    unique_sites['motif'] = assign_motifs_fast(unique_sites, genome)
    motif_map = unique_sites.set_index(['chrom', 'position', 'strand', 'mod_type'])['motif']
    df['motif'] = df.set_index(['chrom', 'position', 'strand', 'mod_type']).index.map(motif_map).values

    # Motif assignment summary
    motif_counts = unique_sites.groupby(['mod_type', 'motif']).size()
    print("\n  Motif assignment (unique sites):")
    print(motif_counts.to_string())

    # 4. Build wide tables (all 3 TPs present)
    print("\nBuilding wide tables (sites with all 3 timepoints)...")
    wide_4mC = build_wide_table(df, '4mC')
    wide_6mA = build_wide_table(df, '6mA')
    print(f"  4mC sites (all TPs): {len(wide_4mC):,}")
    print(f"  6mA sites (all TPs): {len(wide_6mA):,}")

    # Filter by motif
    gccggc_4mC = wide_4mC[wide_4mC['motif'] == 'GCCGGC'].copy()
    aagcccg_6mA = wide_6mA[wide_6mA['motif'] == 'AAGCCCG'].copy()
    print(f"\n  GCCGGC 4mC (all TPs): {len(gccggc_4mC):,}")
    print(f"  AAGCCCG 6mA (all TPs): {len(aagcccg_6mA):,}")

    # 5. Dynamic/Stable classification at all thresholds
    print("\nClassifying dynamic vs stable sites...")
    results = {}
    for thr in DELTA_THRESHOLDS:
        results[thr] = {}
        results[thr]['GCCGGC']  = classify_dynamics(gccggc_4mC,  thr)
        results[thr]['AAGCCCG'] = classify_dynamics(aagcccg_6mA, thr)

    # Main threshold results
    gccggc_main  = results[MAIN_THRESHOLD]['GCCGGC']
    aagcccg_main = results[MAIN_THRESHOLD]['AAGCCCG']

    # 6. Pattern classification for dynamic sites
    print("Classifying patterns for dynamic sites...")
    for motif, df_main in [('GCCGGC', gccggc_main), ('AAGCCCG', aagcccg_main)]:
        dyn_mask = df_main['category'] == 'dynamic'
        df_main.loc[dyn_mask, 'pattern'] = df_main[dyn_mask].apply(classify_pattern, axis=1)
        df_main.loc[~dyn_mask, 'pattern'] = 'stable'

    # 7. Save per-site tables
    print("\nSaving per-site tables...")

    # GCCGGC
    gccggc_out = gccggc_main[
        ['chrom', 'position', 'strand', 'T1_freq', 'T2_freq', 'T3_freq', 'delta', 'category', 'pattern']
    ].sort_values('delta', ascending=False)
    gccggc_out['region'] = gccggc_out['position'].apply(assign_region)
    gccggc_out.to_csv(OUT_TABLES / 'A5_GCCGGC_site_dynamics.tsv', sep='\t', index=False)
    print(f"  Saved: A5_GCCGGC_site_dynamics.tsv ({len(gccggc_out)} sites)")

    # AAGCCCG
    aagcccg_out = aagcccg_main[
        ['chrom', 'position', 'strand', 'T1_freq', 'T2_freq', 'T3_freq', 'delta', 'category', 'pattern']
    ].sort_values('delta', ascending=False)
    aagcccg_out['region'] = aagcccg_out['position'].apply(assign_region)
    aagcccg_out.to_csv(OUT_TABLES / 'A5_AAGCCCG_site_dynamics.tsv', sep='\t', index=False)
    print(f"  Saved: A5_AAGCCCG_site_dynamics.tsv ({len(aagcccg_out)} sites)")

    # 8. Core/Arm chi-square tests
    print("\nCore/Arm chi-square tests...")
    gccggc_stats  = core_arm_chisq(gccggc_main,  'GCCGGC')
    aagcccg_stats = core_arm_chisq(aagcccg_main, 'AAGCCCG')

    for s in [gccggc_stats, aagcccg_stats]:
        print(f"\n  {s['motif']}:")
        print(f"    n_total={s['n_total']}, dynamic={s['n_dynamic']} ({s['dynamic_pct']:.1f}%)")
        print(f"    Dynamic: core={s['dynamic_core']} ({s['dynamic_core_pct']:.1f}%), arm={s['dynamic_arm']} ({s['dynamic_arm_pct']:.1f}%)")
        print(f"    Stable:  core={s['stable_core']} ({s['stable_core_pct']:.1f}%), arm={s['stable_arm']} ({s['stable_arm_pct']:.1f}%)")
        print(f"    chi2={s['chi2']:.3f}, p={s['p_value']:.4f} (dof={s['dof']})")

    # 9. Summary table (sensitivity analysis)
    summary = build_summary(results, DELTA_THRESHOLDS)
    summary.to_csv(OUT_TABLES / 'A5_dynamics_summary.tsv', sep='\t', index=False)
    print(f"\n  Saved: A5_dynamics_summary.tsv")

    # 10. Pattern statistics
    pattern_stats = []
    for motif, df_m in [('GCCGGC', gccggc_main), ('AAGCCCG', aagcccg_main)]:
        dyn = df_m[df_m['category'] == 'dynamic']
        for pat, cnt in dyn['pattern'].value_counts().items():
            pattern_stats.append({
                'motif': motif,
                'pattern': pat,
                'count': cnt,
                'pct_of_dynamic': round(cnt / len(dyn) * 100, 1) if len(dyn) > 0 else 0,
            })
    pd.DataFrame(pattern_stats).to_csv(
        OUT_TABLES / 'A5_pattern_statistics.tsv', sep='\t', index=False
    )

    # 11. Plots
    print("\nGenerating figures...")
    plot_violin(gccggc_main, aagcccg_main)

    gccggc_dyn  = gccggc_main[gccggc_main['category'] == 'dynamic']
    aagcccg_dyn = aagcccg_main[aagcccg_main['category'] == 'dynamic']
    plot_chromosomal_distribution(gccggc_main, aagcccg_main)
    plot_pattern_classification(gccggc_dyn, aagcccg_dyn)

    # 12. Report output
    print("\n" + "=" * 60)
    print("=== RESULTS SUMMARY ===")
    print("=" * 60)
    print(f"\nThreshold: delta > {MAIN_THRESHOLD} (sensitivity: {DELTA_THRESHOLDS})")

    for motif, df_m, stats in [
        ('GCCGGC 4mC', gccggc_main, gccggc_stats),
        ('AAGCCCG 6mA', aagcccg_main, aagcccg_stats),
    ]:
        n_dyn   = (df_m['category'] == 'dynamic').sum()
        n_total = len(df_m)
        print(f"\n{motif}:")
        print(f"  Total sites (all 3 TPs): {n_total}")
        print(f"  Dynamic: {n_dyn} ({n_dyn/n_total*100:.1f}%)")
        print(f"  Stable:  {n_total - n_dyn} ({(n_total-n_dyn)/n_total*100:.1f}%)")
        print(f"  Core/Arm chi2={stats['chi2']:.3f}, p={stats['p_value']:.4f}")
        dyn = df_m[df_m['category'] == 'dynamic']
        if len(dyn) > 0:
            print(f"  Pattern breakdown:")
            for pat, cnt in dyn['pattern'].value_counts().items():
                print(f"    {pat}: {cnt} ({cnt/len(dyn)*100:.0f}%)")

    print("\nDone.")
    return {
        'gccggc': gccggc_main,
        'aagcccg': aagcccg_main,
        'gccggc_stats': gccggc_stats,
        'aagcccg_stats': aagcccg_stats,
        'summary': summary,
    }


if __name__ == '__main__':
    main()
