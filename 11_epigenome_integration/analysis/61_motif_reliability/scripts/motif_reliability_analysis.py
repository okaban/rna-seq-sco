"""
61_motif_reliability: モチーフ解析信頼性強化
=============================================
対応論点:
  A-1: GCCGGCのどのC、AAGCCCGのどのAがメチル化されているか
  A-2: ゲノム全体モチーフ総数 vs メチル化率（閾値50%）
  A-3: Unassignedカテゴリの定義・内訳

実施日: 2026-04-11
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
from collections import Counter
import re

# === パス設定 ===
ANALYSIS_DIR = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/61_motif_reliability")
DATA_DIR = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis")
GENOME_FASTA = Path("/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/GCF_000203835.1_ASM20383v1_genomic.fna")
SITE_SEQ_CSV = DATA_DIR / "07_motif_analysis/methylation_site_sequences.csv"
HIGHCONF_CSV = DATA_DIR / "01_integration/high_confidence_sites_weighted.csv"

OUT_TABLES = ANALYSIS_DIR / "tables"
OUT_FIGS   = ANALYSIS_DIR / "figures"

# ============================================================
# ゲノム読み込み
# ============================================================
def load_genome(fasta_path):
    genome = {}
    with open(fasta_path) as f:
        chrom = None
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if chrom:
                    genome[chrom] = "".join(seq).upper()
                chrom = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if chrom:
            genome[chrom] = "".join(seq).upper()
    return genome

# ============================================================
# モチーフ総数カウント（ゲノム全体）
# ============================================================
def count_motif_in_genome(genome, motif):
    """正鎖＋負鎖両方でモチーフをカウント"""
    rc_map = str.maketrans("ACGT", "TGCA")
    rc_motif = motif.translate(rc_map)[::-1]
    total = 0
    positions = []  # (chrom, pos, strand)
    for chrom, seq in genome.items():
        for m in re.finditer(motif, seq):
            total += 1
            positions.append((chrom, m.start(), '+'))
        if rc_motif != motif:  # パリンドロームでなければ負鎖もカウント
            for m in re.finditer(rc_motif, seq):
                total += 1
                positions.append((chrom, m.start(), '-'))
    return total, positions

# ============================================================
# A-3: Unassigned定義の明確化
# ============================================================
def analyze_unassigned(site_df):
    """6mAサイトのモチーフ帰属を分析"""
    m6a = site_df[site_df['mod_type'] == '6mA'].copy()

    results = []
    for tp in ['T1', 'T2', 'T3']:
        df = m6a[m6a['timepoint'] == tp]
        total = len(df)
        if total == 0:
            continue

        # モチーフ検索（±15bp window内）
        has_aagcccg_fwd = df['sequence'].str.contains('AAGCCCG', na=False)
        has_aagcccg_rev = df['sequence'].str.contains('CGGGCTT', na=False)  # RC of AAGCCCG
        has_aagcccg     = has_aagcccg_fwd | has_aagcccg_rev
        has_gatc        = df['sequence'].str.contains('GATC', na=False)

        only_aagcccg = has_aagcccg & ~has_gatc
        only_gatc    = has_gatc & ~has_aagcccg
        both         = has_aagcccg & has_gatc
        unassigned   = ~has_aagcccg & ~has_gatc

        results.append({
            'timepoint': tp,
            'total': total,
            'AAGCCCG_only': only_aagcccg.sum(),
            'GATC_only': only_gatc.sum(),
            'both': both.sum(),
            'unassigned': unassigned.sum(),
            'AAGCCCG_pct': round((only_aagcccg.sum() + both.sum()) / total * 100, 1),
            'GATC_pct': round((only_gatc.sum() + both.sum()) / total * 100, 1),
            'unassigned_pct': round(unassigned.sum() / total * 100, 1),
        })

    return pd.DataFrame(results)

def analyze_unassigned_sequences(site_df):
    """Unassignedサイトのde novoモチーフ探索（±5bp k-mer）"""
    m6a = site_df[(site_df['mod_type'] == '6mA') & (site_df['timepoint'] == 'T1')].copy()
    has_aagcccg = m6a['sequence'].str.contains('AAGCCCG', na=False) | m6a['sequence'].str.contains('CGGGCTT', na=False)
    has_gatc    = m6a['sequence'].str.contains('GATC', na=False)
    unassigned  = m6a[~has_aagcccg & ~has_gatc]

    # センター±5bp = positions 10-20 (0-indexed) in 31-mer
    center_seqs = unassigned['sequence'].str[10:21]  # 11-mer centered
    kmer_counts = Counter()
    for seq in center_seqs:
        if len(seq) >= 7:
            for i in range(len(seq) - 6):
                kmer_counts[seq[i:i+7]] += 1

    top_kmers = pd.DataFrame(kmer_counts.most_common(20), columns=['kmer', 'count'])
    top_kmers['pct'] = (top_kmers['count'] / len(unassigned) * 100).round(1)
    return top_kmers, len(unassigned)

# ============================================================
# A-1: メチル化塩基の特定（GCCGGC内のC、AAGCCCG内のA）
# ============================================================
def find_methylated_position_in_motif(site_df, motif, mod_type):
    """
    motif内でメチル化された塩基が何番目かを特定する。

    31-mer中心(pos=15)がメチル化塩基。
    motifがseq中のどこに現れるか検索し、
    center(15)がmotif内の何番目かを計算する。
    """
    df = site_df[(site_df['mod_type'] == mod_type) & (site_df['timepoint'] == 'T1')].copy()
    rc_map = str.maketrans("ACGT", "TGCA")
    rc_motif = motif.translate(rc_map)[::-1]

    positions_in_motif = []

    for _, row in df.iterrows():
        seq = row['sequence']  # 31-mer
        center = 15  # 0-indexed, center of 31-mer = the modified base

        # 正鎖側でmotif探索
        for m in re.finditer(motif, seq):
            rel_pos = center - m.start()
            if 0 <= rel_pos < len(motif):
                positions_in_motif.append({
                    'strand_in_seq': '+',
                    'motif': motif,
                    'pos_in_motif': rel_pos,
                    'base_at_pos': motif[rel_pos],
                })

        # 逆鎖モチーフを探索（RC of motif）
        if rc_motif != motif:
            for m in re.finditer(rc_motif, seq):
                # RC motif内でのrelative position
                rel_pos = center - m.start()
                if 0 <= rel_pos < len(rc_motif):
                    # motif上の対応位置（RC変換）
                    pos_in_original = len(motif) - 1 - rel_pos
                    positions_in_motif.append({
                        'strand_in_seq': '-',
                        'motif': motif,
                        'pos_in_motif': pos_in_original,
                        'base_at_pos': motif[pos_in_original],
                    })

    return pd.DataFrame(positions_in_motif)

# ============================================================
# A-2: モチーフ充足率
# ============================================================
def compute_motif_occupancy(genome, site_seq_df, motifs_config):
    """
    各モチーフについて:
      - ゲノム全体の出現サイト数（正鎖＋負鎖）
      - 各タイムポイントでそのモチーフに一致するメチル化サイト数
      - 充足率(%)

    motifs_config: [(motif_seq, mod_type, label), ...]
    """
    rc_map = str.maketrans("ACGT", "TGCA")

    rows = []
    for motif, mod_type, label in motifs_config:
        rc_motif = motif.translate(rc_map)[::-1]
        genome_count, _ = count_motif_in_genome(genome, motif)

        for tp in ['T1', 'T2', 'T3']:
            df_tp = site_seq_df[
                (site_seq_df['mod_type'] == mod_type) &
                (site_seq_df['timepoint'] == tp)
            ]
            # モチーフに一致するサイトだけカウント
            has_motif = (
                df_tp['sequence'].str.contains(motif, na=False) |
                (df_tp['sequence'].str.contains(rc_motif, na=False) if rc_motif != motif else False)
            )
            meth_count = has_motif.sum()

            rows.append({
                'motif': motif,
                'mod_type': mod_type,
                'label': label,
                'timepoint': tp,
                'genome_sites_total': genome_count,
                'methylated_sites': meth_count,
                'occupancy_pct': round(meth_count / genome_count * 100, 2) if genome_count > 0 else 0,
            })

    return pd.DataFrame(rows)

# ============================================================
# 図の作成
# ============================================================
def plot_unassigned_breakdown(unassigned_df, out_path):
    """A-3: タイムポイント別モチーフ帰属のスタックドバー"""
    fig, ax = plt.subplots(figsize=(7, 4))

    tps = unassigned_df['timepoint'].tolist()
    x = np.arange(len(tps))
    width = 0.5

    colors = {'AAGCCCG_only': '#2196F3', 'GATC_only': '#FF9800', 'both': '#9C27B0', 'unassigned': '#9E9E9E'}
    labels = {'AAGCCCG_only': 'AAGCCCG', 'GATC_only': 'GATC', 'both': 'Both', 'unassigned': 'Unassigned'}

    bottom = np.zeros(len(tps))
    for key, color in colors.items():
        vals = unassigned_df[key].values
        pct  = vals / unassigned_df['total'].values * 100
        bars = ax.bar(x, pct, width, bottom=bottom, color=color, label=labels[key], edgecolor='white', linewidth=0.5)
        # パーセント注釈
        for i, (p, b) in enumerate(zip(pct, bottom)):
            if p >= 3:
                ax.text(i, b + p/2, f'{p:.0f}%', ha='center', va='center', fontsize=9, color='white', fontweight='bold')
        bottom += pct

    ax.set_xticks(x)
    ax.set_xticklabels([f'{tp}\n(n={n:,})' for tp, n in zip(tps, unassigned_df['total'])])
    ax.set_ylabel('Sites (%)')
    ax.set_title('6mA Site Motif Assignment by Timepoint\n(A-3: Unassigned Category Definition)')
    ax.legend(loc='lower right', framealpha=0.9)
    ax.set_ylim(0, 105)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(str(out_path).replace('.pdf', '.pdf'), format='pdf')
    plt.savefig(str(out_path).replace('.pdf', '.svg'), format='svg')
    plt.close()

def plot_position_in_motif(pos_df, motif, mod_type, out_path):
    """A-1: モチーフ内の修飾塩基位置のバーチャート"""
    if len(pos_df) == 0:
        print(f"  警告: {motif}の位置データなし")
        return

    counts = pos_df['pos_in_motif'].value_counts().sort_index()
    total = len(pos_df)

    fig, ax = plt.subplots(figsize=(max(5, len(motif)*0.8 + 2), 4))

    # モチーフの各位置を可視化
    mod_base = 'A' if mod_type == '6mA' else 'C'
    colors = []
    for i, base in enumerate(motif):
        if base == mod_base and i in counts.index:
            colors.append('#E53935')  # 修飾塩基候補
        elif base == mod_base:
            colors.append('#FFCDD2')  # 修飾可能だが未検出
        else:
            colors.append('#90CAF9')  # 非修飾塩基

    x = np.arange(len(motif))
    bars_data = [counts.get(i, 0) for i in range(len(motif))]
    pct_data  = [c / total * 100 for c in bars_data]

    bars = ax.bar(x, pct_data, color=colors, edgecolor='gray', linewidth=0.5, width=0.6)

    for i, (p, c) in enumerate(zip(pct_data, bars_data)):
        if p > 0:
            ax.text(i, p + 0.5, f'{p:.1f}%\n(n={c})', ha='center', va='bottom', fontsize=8)

    # 軸設定
    ax.set_xticks(x)
    ax.set_xticklabels([f'pos{i}\n{base}' for i, base in enumerate(motif)], fontsize=10)
    ax.set_ylabel('Sites (%)')
    ax.set_title(f'{mod_type} Methylated Position within {motif} Motif\n(n={total} sites, T1)')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    legend_patches = [
        mpatches.Patch(color='#E53935', label=f'Modified {mod_base} (detected)'),
        mpatches.Patch(color='#90CAF9', label='Other bases'),
    ]
    ax.legend(handles=legend_patches, loc='upper right')

    plt.tight_layout()
    pdf_path = str(out_path)
    plt.savefig(pdf_path, format='pdf')
    plt.savefig(pdf_path.replace('.pdf', '.svg'), format='svg')
    plt.close()
    print(f"  保存: {pdf_path}")

def plot_motif_occupancy(occ_df, out_path):
    """A-2: モチーフ別充足率（ゲノム全体に対するメチル化割合）"""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    motifs = occ_df['label'].unique()
    colors = ['#1565C0', '#E53935', '#2E7D32', '#6A1B9A']

    for ax, tp in zip(axes, ['T1', 'T2']):
        df_tp = occ_df[occ_df['timepoint'] == tp]
        x = np.arange(len(motifs))
        for i, (label, color) in enumerate(zip(motifs, colors)):
            row = df_tp[df_tp['label'] == label]
            if len(row) == 0:
                continue
            occ = row['occupancy_pct'].values[0]
            meth = row['methylated_sites'].values[0]
            total = row['genome_sites_total'].values[0]
            ax.bar(i, occ, color=color, edgecolor='gray', linewidth=0.5, width=0.5)
            ax.text(i, occ + 0.05, f'{occ:.2f}%\n({meth:,}/{total:,})',
                   ha='center', va='bottom', fontsize=8)

        ax.set_xticks(x)
        ax.set_xticklabels(motifs, fontsize=9)
        ax.set_ylabel('Methylated / Genome sites (%)')
        ax.set_title(f'{tp}: Motif Occupancy (freq ≥ 50%)')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.suptitle('A-2: Methylation Occupancy by Motif\n(methylated sites / total genomic motif occurrences)',
                 fontsize=11, y=1.01)
    plt.tight_layout()
    plt.savefig(str(out_path), format='pdf')
    plt.savefig(str(out_path).replace('.pdf', '.svg'), format='svg')
    plt.close()


# ============================================================
# メイン実行
# ============================================================
if __name__ == '__main__':
    print("=== 61_motif_reliability_analysis ===\n")

    # データ読み込み
    print("1. データ読み込み...")
    site_df = pd.read_csv(SITE_SEQ_CSV)
    highconf_df = pd.read_csv(HIGHCONF_CSV)
    print(f"   site_seq: {len(site_df)} rows, highconf: {len(highconf_df)} rows")

    print("2. ゲノム読み込み...")
    genome = load_genome(GENOME_FASTA)
    genome_len = sum(len(v) for v in genome.values())
    print(f"   ゲノム長: {genome_len:,} bp")

    # --- A-3: Unassigned定義 ---
    print("\n--- A-3: Unassigned定義 ---")
    unassigned_df = analyze_unassigned(site_df)
    print(unassigned_df.to_string(index=False))
    unassigned_df.to_csv(OUT_TABLES / "A3_6mA_motif_assignment_by_timepoint.tsv", sep='\t', index=False)

    top_kmers, n_unasgn = analyze_unassigned_sequences(site_df)
    print(f"\n   Unassigned T1サイト数: {n_unasgn}")
    print("   Top 7-mers in Unassigned sites:")
    print(top_kmers.head(10).to_string(index=False))
    top_kmers.to_csv(OUT_TABLES / "A3_unassigned_top_kmers.tsv", sep='\t', index=False)

    plot_unassigned_breakdown(unassigned_df, OUT_FIGS / "A3_6mA_motif_assignment.pdf")
    print("   → figures/A3_6mA_motif_assignment.pdf")

    # --- A-1: メチル化塩基位置特定 ---
    print("\n--- A-1: メチル化塩基位置 ---")

    # AAGCCCG内のA
    print("  [6mA / AAGCCCG]")
    pos6mA = find_methylated_position_in_motif(site_df, 'AAGCCCG', '6mA')
    if len(pos6mA) > 0:
        print(f"   マッチ数: {len(pos6mA)}")
        print(f"   pos_in_motif分布:\n{pos6mA['pos_in_motif'].value_counts().sort_index()}")
        pos6mA.to_csv(OUT_TABLES / "A1_AAGCCCG_methylated_position.tsv", sep='\t', index=False)
        plot_position_in_motif(pos6mA, 'AAGCCCG', '6mA', OUT_FIGS / "A1_AAGCCCG_methylated_position.pdf")

    # GATC内のA
    print("  [6mA / GATC]")
    posGATC = find_methylated_position_in_motif(site_df, 'GATC', '6mA')
    if len(posGATC) > 0:
        print(f"   マッチ数: {len(posGATC)}")
        print(f"   pos_in_motif分布:\n{posGATC['pos_in_motif'].value_counts().sort_index()}")
        posGATC.to_csv(OUT_TABLES / "A1_GATC_methylated_position.tsv", sep='\t', index=False)
        plot_position_in_motif(posGATC, 'GATC', '6mA', OUT_FIGS / "A1_GATC_methylated_position.pdf")

    # GCCGGC内のC（4mC）
    print("  [4mC / GCCGGC]")
    pos4mC = find_methylated_position_in_motif(site_df, 'GCCGGC', '4mC')
    if len(pos4mC) > 0:
        print(f"   マッチ数: {len(pos4mC)}")
        print(f"   pos_in_motif分布:\n{pos4mC['pos_in_motif'].value_counts().sort_index()}")
        pos4mC.to_csv(OUT_TABLES / "A1_GCCGGC_methylated_position.tsv", sep='\t', index=False)
        plot_position_in_motif(pos4mC, 'GCCGGC', '4mC', OUT_FIGS / "A1_GCCGGC_methylated_position.pdf")

    # --- A-2: モチーフ充足率 ---
    print("\n--- A-2: モチーフ充足率 ---")
    motifs_config = [
        ('AAGCCCG', '6mA', 'AAGCCCG\n(6mA)'),
        ('GATC',    '6mA', 'GATC\n(6mA)'),
        ('GCCGGC',  '4mC', 'GCCGGC\n(4mC)'),
    ]
    occ_df = compute_motif_occupancy(genome, site_df, motifs_config)
    print(occ_df[['label','timepoint','genome_sites_total','methylated_sites','occupancy_pct']].to_string(index=False))
    occ_df.to_csv(OUT_TABLES / "A2_motif_occupancy.tsv", sep='\t', index=False)
    plot_motif_occupancy(occ_df, OUT_FIGS / "A2_motif_occupancy.pdf")
    print("   → figures/A2_motif_occupancy.pdf")

    print("\n=== 完了 ===")
    print(f"Tables: {OUT_TABLES}")
    print(f"Figures: {OUT_FIGS}")
