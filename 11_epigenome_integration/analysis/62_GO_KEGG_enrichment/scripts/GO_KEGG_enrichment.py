"""
62_GO_KEGG_enrichment: メチル化遺伝子のGO/KEGG解析
===================================================
対応論点: E-1

方針:
  - GOエンリッチメント: gene_annotation_basic.tsvのGO termを使用、Fisher's exact test + Benjamini-Hochberg補正
  - KEGG: KEGG REST API で S. coelicolor M145 (sco) のパスウェイ取得 → gene-pathway mapping
  - 対象遺伝子群: GCCGGC近傍、AAGCCCG近傍、dual-targeted（H17データから）

実施日: 2026-04-11
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from pathlib import Path
import requests
import time
import json

# === パス設定 ===
ANALYSIS_DIR = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/62_GO_KEGG_enrichment")
DATA_DIR = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis")
GENE_ANNOT = Path("/Users/okaban/bioinfo/rna-seq/05_annotation/analysis/05_annotation_260128_v1/tables/gene_annotation_basic.tsv")
GENE_SET   = DATA_DIR / "40_motif_division_of_labor/tables/gene_set_comparison.tsv"

OUT_TABLES = ANALYSIS_DIR / "tables"
OUT_FIGS   = ANALYSIS_DIR / "figures"

# ============================================================
# GO Enrichment
# ============================================================
def load_go_annotations(annot_path):
    """gene_annotation_basic.tsvからGO term辞書を構築"""
    df = pd.read_csv(annot_path, sep='\t')
    gene2go = {}
    go2genes = {}
    for _, row in df.iterrows():
        gene = str(row['gene_id'])  # plain Python str to avoid pandas dtype issues
        if pd.isna(row['ontology_term']) or row['ontology_term'] == '':
            continue
        terms = [t.strip() for t in str(row['ontology_term']).split(',') if t.strip().startswith('GO:')]
        gene2go[gene] = terms
        for t in terms:
            go2genes.setdefault(t, set()).add(gene)
    # go2genes の値を frozenset に変換（すべて plain str）
    go2genes = {k: frozenset(v) for k, v in go2genes.items()}
    return gene2go, go2genes, df

def run_go_enrichment(query_genes, background_genes, go2genes, min_genes=3, fdr_threshold=0.05):
    """Fisher's exact test によるGO enrichment"""
    query_set = set(str(g) for g in query_genes)
    bg_set = set(str(g) for g in background_genes)

    # queryはbgに含まれることを確認
    query_in_bg = query_set & bg_set
    bg_only = bg_set - query_in_bg

    N = len(bg_set)           # background total
    k = len(query_in_bg)      # query in background

    rows = []
    for term, term_genes in go2genes.items():
        term_in_bg = term_genes & bg_set
        if len(term_in_bg) < min_genes:
            continue
        term_in_query = term_genes & query_in_bg

        # 2×2 contingency table
        a = len(term_in_query)          # query & term
        b = len(query_in_bg) - a        # query & not-term
        c = len(term_in_bg) - a         # not-query & term
        d = N - a - b - c               # not-query & not-term

        if a == 0:
            continue

        odds, pval = fisher_exact([[a, b], [c, d]], alternative='greater')
        rows.append({
            'GO_term': term,
            'n_query': a,
            'n_term_in_bg': len(term_in_bg),
            'n_query_total': k,
            'n_bg_total': N,
            'enrichment': (a/k) / (len(term_in_bg)/N) if k > 0 else 0,
            'odds_ratio': odds,
            'pvalue': pval,
            'gene_list': ','.join(sorted(term_in_query)),
        })

    if not rows:
        return pd.DataFrame()
    result = pd.DataFrame(rows)
    _, padj, _, _ = multipletests(result['pvalue'], method='fdr_bh')
    result['padj'] = padj
    result = result.sort_values('pvalue')
    return result[result['padj'] < fdr_threshold] if fdr_threshold else result.sort_values('pvalue')

# ============================================================
# KEGG via REST API
# ============================================================
def fetch_kegg_pathway_genes(organism='sco', cache_path=None):
    """KEGG REST APIでpathway-gene mappingを取得"""
    if cache_path and Path(cache_path).exists():
        with open(cache_path) as f:
            return json.load(f)

    print("  KEGG APIから pathway一覧取得中...")
    # パスウェイ一覧取得
    url = f"http://rest.kegg.jp/list/pathway/{organism}"
    try:
        r = requests.get(url, timeout=30)
        if r.status_code != 200:
            print(f"  KEGG API エラー: {r.status_code}")
            return {}
    except Exception as e:
        print(f"  KEGG API 接続失敗: {e}")
        return {}

    pathways = {}
    for line in r.text.strip().split('\n'):
        if '\t' in line:
            pid, name = line.split('\t', 1)
            pathways[pid] = name

    print(f"  取得パスウェイ数: {len(pathways)}")

    # 各パスウェイの遺伝子リスト取得
    pw2genes = {}
    for i, (pid, pname) in enumerate(pathways.items()):
        if i % 20 == 0:
            print(f"    {i}/{len(pathways)} 処理中...")
        try:
            url_g = f"http://rest.kegg.jp/get/{pid}"
            r2 = requests.get(url_g, timeout=30)
            if r2.status_code == 200:
                genes = []
                in_gene = False
                for line in r2.text.split('\n'):
                    if line.startswith('GENE'):
                        in_gene = True
                    elif in_gene and line.startswith(' '):
                        parts = line.strip().split()
                        if parts:
                            genes.append(f"{organism}:{parts[0]}")
                    elif in_gene and not line.startswith(' '):
                        break
                if genes:
                    pw2genes[pid] = {'name': pname, 'genes': genes}
        except:
            pass
        time.sleep(0.1)

    if cache_path:
        with open(cache_path, 'w') as f:
            json.dump(pw2genes, f)

    return pw2genes

def map_kegg_to_locus(pw2genes, annot_df):
    """KEGG遺伝子ID (sco:SCO0001形式) → locus_tag → gene_id マッピング"""
    # old_locus_tag (SCO番号) でマッピング
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

def run_kegg_enrichment(query_genes, background_genes, pw2gene_ids, min_genes=2, fdr_threshold=0.2):
    """KEGG pathway enrichment"""
    query_set = set(query_genes) & set(background_genes)
    bg_set = set(background_genes)
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

        odds, pval = fisher_exact([[a, b], [c, d]], alternative='greater')
        rows.append({
            'pathway_id': pid,
            'pathway_name': info['name'],
            'n_query': a,
            'n_pathway': len(pw_genes),
            'enrichment': (a/k) / (len(pw_genes)/N) if k > 0 else 0,
            'odds_ratio': odds,
            'pvalue': pval,
            'gene_list': ','.join(sorted(pw_in_query)),
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
# 可視化
# ============================================================
def plot_enrichment_dotplot(result_df, title, out_path, top_n=15, term_col='GO_term', name_col=None):
    """エンリッチメント結果のドットプロット"""
    if len(result_df) == 0:
        print(f"  警告: {title} - 有意なtermなし")
        return

    df = result_df.head(top_n).copy()
    df['-log10(padj)'] = -np.log10(df['padj'].clip(1e-10))

    # term labelの作成
    if name_col and name_col in df.columns:
        df['label'] = df[name_col].str[:50]
    else:
        df['label'] = df[term_col]

    fig, ax = plt.subplots(figsize=(9, max(4, len(df)*0.35 + 1.5)))

    scatter = ax.scatter(
        df['enrichment'], range(len(df)),
        s=df['n_query'] * 8 + 20,
        c=df['-log10(padj)'],
        cmap='YlOrRd', vmin=0, vmax=max(5, df['-log10(padj)'].max()),
        edgecolors='gray', linewidths=0.5, zorder=3
    )

    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df['label'], fontsize=8)
    ax.set_xlabel('Fold Enrichment', fontsize=10)
    ax.set_title(title, fontsize=11)
    ax.axvline(1, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.grid(axis='x', alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    cbar = plt.colorbar(scatter, ax=ax, shrink=0.5)
    cbar.set_label('-log10(FDR)', fontsize=9)

    plt.tight_layout()
    plt.savefig(str(out_path), format='pdf')
    plt.savefig(str(out_path).replace('.pdf', '.svg'), format='svg')
    plt.close()
    print(f"  保存: {out_path}")

# ============================================================
# メイン
# ============================================================
if __name__ == '__main__':
    print("=== 62_GO_KEGG_enrichment ===\n")

    # --- データ読み込み ---
    print("1. データ読み込み...")
    gene2go, go2genes, _ = load_go_annotations(GENE_ANNOT)
    annot_df = pd.read_csv(GENE_ANNOT, sep='\t')
    gene_sets = pd.read_csv(GENE_SET, sep='\t')
    print(f"   遺伝子アノテーション: {len(annot_df)} genes, GO term保有: {len(gene2go)}")
    print(f"   遺伝子セット: {len(gene_sets)} genes")
    print(f"   exclusive_group分布:\n{gene_sets['exclusive_group'].value_counts()}")

    # --- 遺伝子群の定義 ---
    print("\n2. 遺伝子群設定...")
    all_genes = gene_sets['gene_id'].tolist()
    gccggc_genes = gene_sets[gene_sets['exclusive_group'].str.contains('GCCGGC', na=False)]['gene_id'].tolist()
    aagcccg_genes = gene_sets[gene_sets['exclusive_group'].str.contains('AAGCCCG', na=False)]['gene_id'].tolist()
    dual_genes = gene_sets[gene_sets['exclusive_group'] == 'Dual-targeted']['gene_id'].tolist()

    print(f"   GCCGGC-only: {len(gccggc_genes)}")
    print(f"   AAGCCCG-only: {len(aagcccg_genes)}")
    print(f"   Dual: {len(dual_genes)}")
    print(f"   Background: {len(all_genes)}")

    # --- GO Enrichment ---
    print("\n3. GO Enrichment (Fisher's exact test + BH)...")

    analyses = [
        ('GCCGGC-proximal (4mC)', gccggc_genes),
        ('AAGCCCG-proximal (6mA)', aagcccg_genes),
        ('Dual-targeted', dual_genes),
    ]

    go_all_results = {}
    for name, query in analyses:
        print(f"\n  [{name}] n={len(query)}")
        result = run_go_enrichment(query, all_genes, go2genes, min_genes=3, fdr_threshold=None)
        if len(result) > 0:
            sig = result[result['padj'] < 0.05]
            print(f"   有意GO terms (FDR<0.05): {len(sig)}")
            if len(sig) > 0:
                print(f"   Top 5:")
                print(sig[['GO_term','n_query','n_term_in_bg','enrichment','pvalue','padj']].head(5).to_string(index=False))
        go_all_results[name] = result
        safe_name = name.replace(' ', '_').replace('/', '_').replace('(', '').replace(')', '')
        result.to_csv(OUT_TABLES / f"E1_GO_enrichment_{safe_name}.tsv", sep='\t', index=False)
        plot_enrichment_dotplot(
            result[result['padj'] < 0.1].head(20),
            f'GO Enrichment: {name}',
            OUT_FIGS / f"E1_GO_enrichment_{safe_name}.pdf",
        )

    # --- KEGG Enrichment ---
    print("\n4. KEGG Pathway Enrichment...")
    cache_file = OUT_TABLES / "kegg_sco_pathways_cache.json"
    pw2genes = fetch_kegg_pathway_genes('sco', cache_path=str(cache_file))

    if pw2genes:
        print(f"   取得パスウェイ数: {len(pw2genes)}")
        pw2gene_ids = map_kegg_to_locus(pw2genes, annot_df)
        print(f"   マッピング済みパスウェイ: {len(pw2gene_ids)}")

        kegg_results = {}
        for name, query in analyses:
            print(f"\n  [{name}]")
            result = run_kegg_enrichment(query, all_genes, pw2gene_ids, min_genes=2, fdr_threshold=None)
            if len(result) > 0:
                sig = result[result['padj'] < 0.2]
                print(f"   有意KEGG pathways (FDR<0.2): {len(sig)}")
                if len(sig) > 0:
                    print(f"   Top 5:")
                    print(sig[['pathway_id','pathway_name','n_query','enrichment','pvalue','padj']].head(5).to_string(index=False))
            kegg_results[name] = result
            safe_name = name.replace(' ', '_').replace('/', '_').replace('(', '').replace(')', '')
            result.to_csv(OUT_TABLES / f"E1_KEGG_enrichment_{safe_name}.tsv", sep='\t', index=False)
            if len(result) > 0:
                plot_enrichment_dotplot(
                    result.head(15),
                    f'KEGG Enrichment: {name}',
                    OUT_FIGS / f"E1_KEGG_enrichment_{safe_name}.pdf",
                    term_col='pathway_id',
                    name_col='pathway_name',
                )
    else:
        print("  KEGG APIが利用不可のため、COGベースの代替解析を使用します")
        # COGカテゴリを使用（既存データから）
        cog_data = pd.read_csv(DATA_DIR / "46_category_specificity_avoidance/tables/category_enrichment.tsv", sep='\t')
        for name, _ in analyses:
            motif_key = 'GCCGGC_4mC_T1' if 'GCCGGC' in name else 'AAGCCCG_6mA_T1'
            subset = cog_data[cog_data['methyl_type'] == motif_key].copy()
            subset = subset[subset['p_bonferroni'] < 0.1].sort_values('fold_enrichment', ascending=False)
            print(f"\n  [{name}] COGエンリッチメント（既存H23データ）:")
            print(subset[['category','n_proximal','fold_enrichment','p_bonferroni']].to_string(index=False))

    print("\n=== 完了 ===")
