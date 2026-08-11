# A群（即修正5件）図の再生成レポート

- 日付: 2026-05-04
- 対象論文: S. coelicolor A3(2) M145 RNA-seq + epigenome integration
- スコープ: 17件のレビューコメントのうち、A群（即修正可能な視覚化）5件
- 作業ディレクトリ: `15_paper_figures/`

---

## 結論サマリ

5件すべて完了。出力は `15_paper_figures/figures/main/` に PDF + SVG。
旧版PDFは `15_paper_figures/figures/main/archive/preA-group_260504/` に退避。

| # | 項目 | 入力スクリプト | 出力ファイル | サンプルサイズ |
|---|---|---|---|---|
| A1 | Fig 1c デュアルスケール化 | `01_figure1_landscape.py` (修正) | `Figure1_methylation_landscape.{pdf,svg}` | 4mC 6,080 / 6mA 6,349 sites |
| A2/A4 | Fig 8 (Unassigned 6mA削除 + 品質向上) | `29_figure8_positional_stratification.py` (新規) | `Figure8_positional_stratification.{pdf,svg}` | GCCGGC 1,659 / AAGCCCG 1,078 sites (T1) |
| A3 | Fig 7 KEGG bubble 品質向上 | `28_figure7_kegg_bubble.py` (新規) | `Figure7_kegg_bubble.{pdf,svg}` | 263 enrichment行, 上位12 pathway |
| A5 | Fig 4 ST5 格上げ + 数値更新 | `27_figure4_shielded_exposed.py` (修正) + `30_table_exposed_TFs_main.py` (新規) | `Figure4_shielded_exposed.{pdf,svg}` + `Table1_exposed_TFs_main.{tsv,html}` + `ST5_exposed_TFs_full_n57.tsv` | 1,055 reg. genes (1,051 TSSあり、57 exposed / 994 shielded) |

---

## 重要な数値変更: n=57/955 → n=57/998

### 検証

`exposed_regulators_full_table.tsv` (51_exposed_regulators_characteristics) と
`all_regulatory_genes.tsv` (29_genomewide_TF_screen) のクロスチェックで確定:

- **総制御遺伝子: 1,055** (`all_regulatory_genes.tsv`)
- **Exposed: 57** (`exposed_regulators_full_table.tsv`)
- **Shielded: 998** (= 1,055 − 57)

旧版 (260312 archive: 1017/57/955) と異なる。
論文本文の `n=57` は `n=57`、`n=955` は `n=998` に置換が必要。

### 実装

新規統合データを作成:
- `11_epigenome_integration/analysis/52_shielded_exposed_boundary/tables/all_genes_features_unified_n57.tsv`
- 1,055行、`is_exposed` フラグを最新の57遺伝子で再付与
- 既存の260312アーカイブ (1,017行) + 38遺伝子分の `nearest_methyl_distance` を新規計算
- 4遺伝子はTSSが取れず ROC 解析からは除外 → 図では 1,051遺伝子 (57 exposed / 994 shielded)

---

## 各項目の詳細

### A1: Fig 1c デュアルスケール

**問題**: CDS が約80%を占めるため、Promoter / 5'UTR / Intergenic の差が線形バーでは判別困難。

**実装** (`01_figure1_landscape.py` の `panel_c_genomic_distribution`):
- メインパネルは線形 % のまま据え置き (CDS 含む全カテゴリ表示)
- インセット軸を新設 (`mpl_toolkits.axes_grid1.inset_locator.inset_axes`)
- インセットは Promoter / 5'UTR / Intergenic 限定で 0–15% レンジに自動拡大
- 有意マーカー (▲▼ + *) もインセット側に重複表示

### A2 + A4: Fig 8 全面再構成 (Unassigned 6mA 削除)

**問題**: 旧 Fig 8a の Unassigned 6mA カテゴリは生物学的モチーフ非対応のため削除を要望。
全体としてもクオリティ不足。

**実装** (`29_figure8_positional_stratification.py` 新規):
- Panel A: 水平スタックバー、GCCGGC (4mC) と AAGCCCG (m4C/6mA) のみ表示
- Panel A: Promoter / 5'UTR / CDS internal / Intergenic の4カテゴリ
- Panel A: 各セグメントに % ラベル、合計サイト数を Y 軸に明記
- Panel B: 位置層別 KEGG 富化のバブルヒートマップ (motif × category × pathway)
- Panel B: バブルサイズ = −log10(FDR)、色 = log2(OR)、★ = FDR < 0.05、太枠 = FDR < 0.10
- データソース: `62_GO_KEGG_enrichment/tables/F1_classified_meth_sites.tsv` および `F1_stratified_KEGG_by_position_v2.tsv`

### A3: Fig 7 KEGG bubble 品質向上

**問題**: 旧 Fig 7 (`Figure7_enrichment.pdf`) は単軸に近く、3モチーフ比較が不明瞭。

**実装** (`28_figure7_kegg_bubble.py` 新規):
- 入力: `E1_KEGG_enrichment_GCCGGC-proximal_4mC.tsv` 他3ファイル
- レイアウト: pathway × motif の2D bubble
- バブル色: GCCGGC (赤) / AAGCCCG (紫) / Dual-targeted (青)
- バブルサイズ: pathway 内 query 遺伝子数
- 縁: 太線 = FDR < 0.10、細線 = それ以外
- ★ overlay: FDR < 0.05
- 軸: pathway 名は短縮ルール適用 (`Streptomyces coelicolor` → `S. coelicolor`、 `biosynthesis` → `biosynth.`)
- 表示 pathway: 12個 (best padj 順)

主要発見 (FDR < 0.05):
- sco00975 Biosynthesis of various siderophores (padj_min = 9.9e-04)
- sco02024 Quorum sensing (padj_min = 0.024)
- sco00790 Folate biosynthesis (padj_min = 0.029)
- sco00190 Oxidative phosphorylation (padj_min = 0.039)
- sco00550 Peptidoglycan biosynthesis (padj_min = 0.041)

### A5: Fig 4 数値更新 + ST5 本文表化

**問題**: Fig 4 は 260312 アーカイブの古い数値 (n=57/955) を使用していた。
さらに ST5 (57 exposed TFs) は補助表のままで本文には載っていない。

**実装**:
- `27_figure4_shielded_exposed.py` を `all_genes_features_unified_n57.tsv` を読むよう変更
- 新規 `30_table_exposed_TFs_main.py` で:
  - `ST5_exposed_TFs_full_n57.tsv` (Online Resource 用、29カラム × 57行)
  - `Table1_exposed_TFs_main.tsv` (本文用、10カラム × 57行、arm先頭、baseMean降順)
  - `Table1_exposed_TFs_main.html` (LaTeX 取り込み用)
- 出力先: `15_paper_figures/tables/main/Table1_exposed_TFs_main.*`

57 exposed TF の構成:
- region: arm 26 / core 31
- bloc (temporal): activation 56 / unassigned 1 (repression なし — 全 exposed が活性化型)
- 主要TFファミリー: HTH (other) 10、Other regulatory 10、TetR 9、Sigma factor 6、Sensor kinase 6

---

## 残課題

- `15_paper_figures/scripts/16_new_supplementary_tables.py` のコメント `ST5: 57 exposed TFs` を `57 exposed TFs` に統一
- 本文中の `n=57` `n=955` `n=1,017` 表記の検索 → `n=57` `n=998` `n=1,055` に置換
- `paper_figures_generation.py` 等の旧スクリプトに残る古い数値の点検
- 旧 `Figure7_enrichment.pdf` を引用していた本文の差し替え (`Figure7_kegg_bubble.pdf` に変更)
- 旧 `Figure8_positional_stratification.pdf` (Unassigned 6mA入り) を引用していた本文 / キャプションの差し替え

---

## 次のアクション

A群は完了。次に着手すべきは:

1. B群 #6,#7,#8 (統計記述追補) — 並行可
2. C群 #11,#12 (既存データ再構成) — 1日
3. D群 #16,#17 (Figure 5/6 再設計) — 3日
