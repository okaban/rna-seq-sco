# 12_supplementary_figures

## 概要

補足図表の生成ステップ。DEGベン図、GO/KEGGエンリッチメント解析、COG分類、カバレッジトラック等の補足的な可視化・解析を行う。

## ディレクトリ構造

```
12_supplementary_figures/
├── analysis/
│   └── 12_supplementary_260202_v1/     # 補足図表結果
│       ├── figures/                     # 各種補足図
│       ├── tables/                      # 補足テーブル
│       └── logs/                        # 実行ログ
├── data/                                # 入力データ
├── scripts/
│   ├── 01_venn_diagram_DEGs.R           # DEGベン図
│   ├── 02_GO_KEGG_enrichment.R          # GO/KEGGエンリッチメント
│   ├── 02b_KEGG_enrichment.R            # KEGG追加解析
│   ├── 03_COG_classification.R          # COG機能分類
│   ├── 04_coverage_tracks.R             # カバレッジトラック
│   └── 05_report_tables_visualization.R # レポート用テーブル
```

## 主要な出力

- `analysis/12_supplementary_260202_v1/figures/` - ベン図、エンリッチメントプロット、COG分類図等
- `analysis/12_supplementary_260202_v1/tables/` - GO/KEGGエンリッチメント結果テーブル
