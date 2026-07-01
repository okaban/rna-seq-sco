"""
62_GO_KEGG_enrichment: 位置層別KEGG/COGエンリッチメント
=========================================================
対応論点: E-1b (位置層別)

分類カテゴリ:
  - promoter      : TSS上流1〜500bp
  - 5UTR_approx   : TSS下流0〜100bp（5'-UTR近似）
  - CDS_internal  : CDS内（上記以外のgene body）
  - intergenic    : 遺伝子間領域

実施日: 2026-04-18
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from pathlib import Path
import json
import sys
import re

# === パス設定 ===
ANALYSIS_DIR = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/62_GO_KEGG_enrichment")
DATA_DIR     = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis")
GENE_ANNOT   = Path("/Users/okaban/bioinfo/rna-seq/05_annotation/analysis/05_annotation_260128_v1/tables/gene_annotation_basic.tsv")
TSS_TABLE    = DATA_DIR / "18_tss_analyses/comprehensive_tss_table.csv"
METH_SITES   = DATA_DIR / "07_motif_analysis/methylation_site_sequences.csv"
COG_FILE     = Path("/Users/okaban/bioinfo/rna-seq/12_supplementary_figures/analysis/12_supplementary_260202_v1/tables/gene_COG_classification.tsv")
KEGG_CACHE   = ANALYSIS_DIR / "tables/kegg_sco_pathways_cache.json"

OUT_TABLES = ANALYSIS_DIR / "tables"
OUT_FIGS   = ANALYSIS_DIR / "figures"


# ============================================================
# 1. データ読み込み
# ============================================================
def load_data():
    annot = pd.read_csv(GENE_ANNOT, sep='\t')
    tss_df = pd.read_csv(TSS_TABLE)
    meth = pd.read_csv(METH_SITES)
    cog = pd.read_csv(COG_FILE, sep='\t')

    # 必要列のみ
    tss_df = tss_df[['gene_id','chrom','start','end','strand','tss']].copy()
    tss_df.columns = ['gene_id','chrom','gene_start','gene_end','strand','tss']

    return annot, tss_df, meth, cog


# ============================================================
# 2. メチル化サイトにモチーフ割り当て + 位置カテゴリ割り当て
# ============================================================
MOTIF_PATTERNS = {
    'GCCGGC': re.compile(r'GCCGGC'),
    'GGCCGG': re.compile(r'GGCCGG'),   # GCCGGC reverse complement
    'AAGCCCG': re.compile(r'AAGCCCG'),
    'CGGGCTT': re.compile(r'CGGGCTT'),  # AAGCCCG reverse complement
}

def assign_motif(seq):
    """配列コンテキストからモチーフを割り当て（GCCGGC優先）"""
    if not isinstance(seq, str):
        return None
    if MOTIF_PATTERNS['GCCGGC'].search(seq) or MOTIF_PATTERNS['GGCCGG'].search(seq):
        return 'GCCGGC'
    if MOTIF_PATTERNS['AAGCCCG'].search(seq) or MOTIF_PATTERNS['CGGGCTT'].search(seq):
        return 'AAGCCCG'
    return None


def assign_position_category(site_pos, gene_start, gene_end, strand, tss):
    """
    メチル化サイト位置 → カテゴリ割り当て

    Parameters
    ----------
    site_pos : int   ゲノム座標（0-based）
    gene_start, gene_end : int  遺伝子の開始/終了（0-based, GFF準拠）
    strand : str  '+' or '-'
    tss : int    TSS座標（0-based）

    Returns
    -------
    str: 'promoter' | '5UTR_approx' | 'CDS_internal' | 'intergenic'
    """
    p = site_pos
    if strand == '+':
        dist_from_tss = p - tss            # 正: downstream, 負: upstream
    else:
        dist_from_tss = tss - p            # 正: downstream, 負: upstream（逆鎖）

    # 遺伝子内かどうか
    in_gene = (gene_start <= p <= gene_end)

    if not in_gene:
        # プロモーター: TSS上流1〜500bp
        if -500 <= dist_from_tss < -1:
            return 'promoter'
        return 'intergenic'

    # 遺伝子内
    if 0 <= dist_from_tss <= 100:
        return '5UTR_approx'
    elif dist_from_tss < 0:
        # 遺伝子内だが上流（TSSより先）→ intergenic扱い
        return 'intergenic'
    else:
        return 'CDS_internal'


def classify_sites(meth: pd.DataFrame, tss_df: pd.DataFrame):
    """
    各メチル化サイトにモチーフとゲノム位置カテゴリを割り当て
    最近傍遺伝子（TSS距離で判断）を割り当てる
    """
    # T1のみ対象
    m = meth[meth['timepoint'] == 'T1'].copy()
    m['motif'] = m['sequence'].apply(assign_motif)
    m = m[m['motif'].notna()].reset_index(drop=True)
    print(f"  T1 GCCGGC sites: {(m['motif']=='GCCGGC').sum()}")
    print(f"  T1 AAGCCCG sites: {(m['motif']=='AAGCCCG').sum()}")

    rows = []
    for _, site in m.iterrows():
        chrom = site['chrom']
        pos   = int(site['position'])
        motif = site['motif']

        # 同染色体遺伝子のみ
        candidates = tss_df[tss_df['chrom'] == chrom].copy()
        if len(candidates) == 0:
            rows.append({'chrom': chrom, 'position': pos, 'motif': motif,
                         'category': 'intergenic', 'gene_id': None})
            continue

        # 全遺伝子との距離を計算してプロモーター or CDS内の遺伝子を探す
        assigned = False
        for _, gene in candidates.iterrows():
            gs = int(gene['gene_start'])
            ge = int(gene['gene_end'])
            st = gene['strand']
            tss = int(gene['tss'])
            gid = gene['gene_id']

            in_gene = (gs <= pos <= ge)
            if st == '+':
                dist = pos - tss
            else:
                dist = tss - pos

            # プロモーターまたはCDS内に分類
            if in_gene or (-500 <= dist < 0):
                cat = assign_position_category(pos, gs, ge, st, tss)
                rows.append({'chrom': chrom, 'position': pos, 'motif': motif,
                             'category': cat, 'gene_id': gid,
                             'dist_from_tss': dist})
                assigned = True
                # 最初のマッチを使用（複数遺伝子重複は先勝ち）
                # より良い実装: 最も近いものを選ぶ
                break

        if not assigned:
            rows.append({'chrom': chrom, 'position': pos, 'motif': motif,
                         'category': 'intergenic', 'gene_id': None,
                         'dist_from_tss': None})

    classified = pd.DataFrame(rows)
    print(f"\n  位置カテゴリ分布:")
    for motif_name, grp in classified.groupby('motif'):
        print(f"  [{motif_name}]")
        print(grp['category'].value_counts().to_string())
    return classified


def classify_sites_vectorized(meth: pd.DataFrame, tss_df: pd.DataFrame):
    """
    ベクトル化バージョン: 各サイトに最近傍遺伝子を割り当てて位置カテゴリを付与
    """
    m = meth[meth['timepoint'] == 'T1'].copy()
    m['motif'] = m['sequence'].apply(assign_motif)
    m = m[m['motif'].notna()].reset_index(drop=True)
    print(f"  T1 GCCGGC sites: {(m['motif']=='GCCGGC').sum()}")
    print(f"  T1 AAGCCCG sites: {(m['motif']=='AAGCCCG').sum()}")

    results = []
    chrom = 'NC_003888.3'
    genes_sorted = tss_df[tss_df['chrom'] == chrom].copy()
    genes_sorted = genes_sorted.sort_values('gene_start').reset_index(drop=True)

    for _, site in m.iterrows():
        pos   = int(site['position'])
        motif = site['motif']

        # 1. gene bodyに含まれるか確認
        in_gene_mask = (genes_sorted['gene_start'] <= pos) & (pos <= genes_sorted['gene_end'])
        in_gene_rows = genes_sorted[in_gene_mask]

        if len(in_gene_rows) > 0:
            # 複数ある場合は最も近いTSSの遺伝子を選ぶ
            g = in_gene_rows.iloc[0]
            tss = int(g['tss'])
            dist = pos - tss if g['strand'] == '+' else tss - pos
            cat = assign_position_category(
                pos, int(g['gene_start']), int(g['gene_end']), g['strand'], tss
            )
            results.append({
                'position': pos, 'motif': motif, 'category': cat,
                'gene_id': g['gene_id'], 'dist_from_tss': dist
            })
            continue

        # 2. プロモーター領域か確認（TSS上流1-500bp）
        prom_candidates = []
        for _, gene in genes_sorted.iterrows():
            tss = int(gene['tss'])
            if gene['strand'] == '+':
                dist = pos - tss
            else:
                dist = tss - pos
            if -500 <= dist <= -1:
                prom_candidates.append((abs(dist), gene, dist))

        if prom_candidates:
            prom_candidates.sort(key=lambda x: x[0])
            _, gene, dist = prom_candidates[0]
            results.append({
                'position': pos, 'motif': motif, 'category': 'promoter',
                'gene_id': gene['gene_id'], 'dist_from_tss': dist
            })
            continue

        # 3. intergenic
        results.append({
            'position': pos, 'motif': motif, 'category': 'intergenic',
            'gene_id': None, 'dist_from_tss': None
        })

    classified = pd.DataFrame(results)
    print(f"\n  位置カテゴリ分布:")
    for motif_name, grp in classified.groupby('motif'):
        print(f"  [{motif_name}]")
        print(grp['category'].value_counts().to_string())
    return classified


# ============================================================
# 3. KEGG エンリッチメント
# ============================================================
def load_kegg_cache(cache_path):
    with open(cache_path) as f:
        return json.load(f)


def map_kegg_to_locus(pw2genes, annot_df):
    locus2gene = {}
    for _, row in annot_df.iterrows():
        if pd.notna(row.get('old_locus_tag')) and row['old_locus_tag']:
            sco = str(row['old_locus_tag']).strip()
            locus2gene[f"sco:{sco.lower()}"] = row['gene_id']
            locus2gene[sco] = row['gene_id']

    pw2gene_ids = {}
    for pid, info in pw2genes.items():
        mapped = []
        for kg in info.get('genes', []):
            if kg.lower() in locus2gene:
                mapped.append(locus2gene[kg.lower()])
        if mapped:
            pw2gene_ids[pid] = {'name': info['name'], 'gene_ids': mapped}
    return pw2gene_ids


def run_kegg_enrichment(query_genes, background_genes, pw2gene_ids, min_genes=2):
    query_set = set(query_genes) & set(background_genes)
    bg_set = set(background_genes)
    if len(query_set) == 0:
        return pd.DataFrame()
    N = len(bg_set)
    k = len(query_set)

    rows = []
    for pid, info in pw2gene_ids.items():
        pw_genes = set(info['gene_ids']) & bg_set
        if len(pw_genes) < min_genes:
            continue
        pw_in_query = pw_genes & query_set
        if len(pw_in_query) == 0:
            continue

        a = len(pw_in_query)
        b = k - a
        c = len(pw_genes) - a
        d = N - a - b - c
        if d < 0:
            continue

        odds, pval = fisher_exact([[a, b], [c, d]], alternative='greater')
        rows.append({
            'pathway_id': pid,
            'pathway_name': info['name'],
            'n_query': a,
            'n_pathway': len(pw_genes),
            'n_query_total': k,
            'odds_ratio': odds,
            'enrichment': (a / k) / (len(pw_genes) / N) if k > 0 else 0,
            'pvalue': pval,
        })

    if not rows:
        return pd.DataFrame()
    result = pd.DataFrame(rows)
    if len(result) > 1:
        _, padj, _, _ = multipletests(result['pvalue'], method='fdr_bh')
        result['padj'] = padj
    else:
        result['padj'] = result['pvalue']
    return result.sort_values('pvalue')


# ============================================================
# 4. COGエンリッチメント
# ============================================================
def run_cog_enrichment(query_genes, background_genes, cog_df, min_genes=3):
    """COGカテゴリ別のFisher's exact test"""
    query_set = set(query_genes) & set(background_genes)
    bg_set = set(background_genes)
    if len(query_set) == 0:
        return pd.DataFrame()

    # COGカテゴリを単一文字に正規化
    cog_map = {}
    for _, row in cog_df.iterrows():
        gid = row['gene_id']
        cat_str = str(row['COG_category'])
        # "X - Description" 形式 → カテゴリ文字と説明を分離
        parts = cat_str.split(' - ', 1)
        cat_code = parts[0].strip()
        cat_desc = parts[1].strip() if len(parts) > 1 else cat_str
        cog_map[gid] = (cat_code, cat_desc)

    # カテゴリ別の遺伝子セット作成
    cat2genes = {}
    cat2desc  = {}
    for gid, (code, desc) in cog_map.items():
        if gid not in bg_set:
            continue
        # 複数カテゴリ（例: "JK"）に対応
        for c in code:
            cat2genes.setdefault(c, set()).add(gid)
            cat2desc[c] = desc if len(code) == 1 else f"Multiple ({code})"

    N = len(bg_set)
    k = len(query_set)

    rows = []
    for cat, cat_genes in cat2genes.items():
        cat_in_bg = cat_genes & bg_set
        if len(cat_in_bg) < min_genes:
            continue
        cat_in_query = cat_genes & query_set
        if len(cat_in_query) == 0:
            continue

        a = len(cat_in_query)
        b = k - a
        c_n = len(cat_in_bg) - a
        d = N - a - b - c_n
        if d < 0:
            continue

        odds, pval = fisher_exact([[a, b], [c_n, d]], alternative='greater')
        rows.append({
            'COG_category': cat,
            'COG_description': cat2desc.get(cat, cat),
            'n_query': a,
            'n_category': len(cat_in_bg),
            'n_query_total': k,
            'odds_ratio': odds,
            'enrichment': (a / k) / (len(cat_in_bg) / N) if k > 0 else 0,
            'pvalue': pval,
        })

    if not rows:
        return pd.DataFrame()
    result = pd.DataFrame(rows)
    if len(result) > 1:
        _, padj, _, _ = multipletests(result['pvalue'], method='fdr_bh')
        result['padj'] = padj
    else:
        result['padj'] = result['pvalue']
    return result.sort_values('pvalue')


# ============================================================
# 5. サマリーテーブル作成
# ============================================================
def make_summary_table(classified, all_genes, pw2gene_ids, cog_df, top_n=5):
    """
    カテゴリ × モチーフ のサマリーテーブルを作成
    """
    categories = ['promoter', '5UTR_approx', 'CDS_internal', 'intergenic']
    motifs = ['GCCGGC', 'AAGCCCG']
    motif_labels = {'GCCGGC': 'GCCGGC 4mC', 'AAGCCCG': 'AAGCCCG 6mA'}

    kegg_rows = []
    cog_rows  = []

    for motif in motifs:
        label = motif_labels[motif]
        for cat in categories:
            subset = classified[(classified['motif'] == motif) & (classified['category'] == cat)]
            genes  = [g for g in subset['gene_id'].dropna().unique().tolist()]
            n_sites = len(subset)

            if len(genes) < 2 or cat == 'intergenic':
                kegg_rows.append({'motif': label, 'category': cat,
                                  'n_sites': n_sites, 'n_genes': len(genes),
                                  'pathway': '(遺伝子数不足)', 'OR': '-', 'padj': '-'})
                cog_rows.append({'motif': label, 'category': cat,
                                 'n_sites': n_sites, 'n_genes': len(genes),
                                 'COG_category': '(遺伝子数不足)', 'COG_description': '-',
                                 'OR': '-', 'padj': '-'})
                continue

            # KEGG
            kegg_res = run_kegg_enrichment(genes, all_genes, pw2gene_ids)
            if len(kegg_res) == 0:
                kegg_rows.append({'motif': label, 'category': cat,
                                  'n_sites': n_sites, 'n_genes': len(genes),
                                  'pathway': '(シグナルなし)', 'OR': '-', 'padj': '-'})
            else:
                for i, row in kegg_res.head(top_n).iterrows():
                    kegg_rows.append({
                        'motif': label, 'category': cat,
                        'n_sites': n_sites if i == kegg_res.index[0] else '',
                        'n_genes': len(genes) if i == kegg_res.index[0] else '',
                        'pathway': row['pathway_name'],
                        'OR': round(row['odds_ratio'], 3),
                        'padj': round(row['padj'], 4),
                    })

            # COG
            cog_res = run_cog_enrichment(genes, all_genes, cog_df)
            if len(cog_res) == 0:
                cog_rows.append({'motif': label, 'category': cat,
                                 'n_sites': n_sites, 'n_genes': len(genes),
                                 'COG_category': '(シグナルなし)', 'COG_description': '-',
                                 'OR': '-', 'padj': '-'})
            else:
                for i, row in cog_res.head(top_n).iterrows():
                    cog_rows.append({
                        'motif': label, 'category': cat,
                        'n_sites': n_sites if i == cog_res.index[0] else '',
                        'n_genes': len(genes) if i == cog_res.index[0] else '',
                        'COG_category': row['COG_category'],
                        'COG_description': row['COG_description'],
                        'OR': round(row['odds_ratio'], 3),
                        'padj': round(row['padj'], 4),
                    })

    return pd.DataFrame(kegg_rows), pd.DataFrame(cog_rows)


# ============================================================
# メイン
# ============================================================
if __name__ == '__main__':
    print("=== 位置層別 KEGG/COG エンリッチメント ===\n")

    # 1. データ読み込み
    print("1. データ読み込み...")
    annot, tss_df, meth, cog = load_data()
    print(f"   アノテーション: {len(annot)} genes")
    print(f"   TSS table: {len(tss_df)} genes")
    print(f"   Meth sites: {len(meth)} (全タイムポイント)")
    print(f"   COG: {len(cog)} genes")

    # 2. サイト分類
    print("\n2. メチル化サイトの位置カテゴリ分類...")
    classified = classify_sites_vectorized(meth, tss_df)
    classified.to_csv(OUT_TABLES / "F1_classified_meth_sites.tsv", sep='\t', index=False)
    print(f"   保存: F1_classified_meth_sites.tsv")

    # 3. KEGG キャッシュ読み込み
    print("\n3. KEGGキャッシュ読み込み...")
    pw2genes_raw = load_kegg_cache(KEGG_CACHE)
    pw2gene_ids  = map_kegg_to_locus(pw2genes_raw, annot)
    all_genes    = annot['gene_id'].tolist()
    print(f"   KEGGパスウェイ: {len(pw2gene_ids)}")
    print(f"   バックグラウンド: {len(all_genes)} genes")

    # 4. サマリーテーブル作成
    print("\n4. カテゴリ別エンリッチメント計算...")
    kegg_summary, cog_summary = make_summary_table(classified, all_genes, pw2gene_ids, cog)

    # 保存
    kegg_out = OUT_TABLES / "F1_stratified_KEGG_by_position_v2.tsv"
    cog_out  = OUT_TABLES / "F1_stratified_COG_by_position.tsv"
    kegg_summary.to_csv(kegg_out, sep='\t', index=False)
    cog_summary.to_csv(cog_out, sep='\t', index=False)
    print(f"   保存: {kegg_out.name}")
    print(f"   保存: {cog_out.name}")

    # 5. 詳細結果の表示
    print("\n" + "="*60)
    print("KEGG エンリッチメント サマリー")
    print("="*60)
    print(kegg_summary.to_string(index=False))

    print("\n" + "="*60)
    print("COG エンリッチメント サマリー")
    print("="*60)
    print(cog_summary.to_string(index=False))

    # 6. 有意な結果のハイライト
    print("\n" + "="*60)
    print("有意な結果 (padj < 0.1)")
    print("="*60)
    for label, df in [("KEGG", kegg_summary), ("COG", cog_summary)]:
        sig = df[df['padj'].apply(lambda x: isinstance(x, (int, float)) and x < 0.1)]
        print(f"\n[{label}] {len(sig)} 有意な組み合わせ:")
        if len(sig) > 0:
            print(sig.to_string(index=False))

    print("\n=== 完了 ===")
