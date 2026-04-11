# 14_DMG_functional_enrichment

## 概要

差次的メチル化遺伝子（DMG）の機能エンリッチメント解析ステップ。DMGの機能カテゴリ偏りを統計的に評価し、メチル化の選択性（selectivity）を深掘り解析する。

## ディレクトリ構造

```
14_DMG_functional_enrichment/
├── analysis/
│   ├── 14_DMG_enrichment_260206_v1/       # 機能エンリッチメント解析
│   │   ├── figures/                        # エンリッチメント可視化
│   │   └── tables/                         # エンリッチメント結果テーブル
│   └── 14_DMG_selectivity_260206_v1/      # メチル化選択性解析
│       ├── figures/                        # 選択性可視化
│       └── tables/                         # 選択性統計テーブル
├── data/                                   # 入力データ
├── scripts/
│   ├── 01_DMG_functional_enrichment.py    # エンリッチメント解析
│   └── 02_DMG_selectivity_deep_analysis.py # 選択性深掘り解析
```

## 主要な出力

- `analysis/14_DMG_enrichment_260206_v1/tables/` - GO/KEGG/COGエンリッチメント結果
- `analysis/14_DMG_enrichment_260206_v1/figures/` - エンリッチメントプロット
- `analysis/14_DMG_selectivity_260206_v1/tables/` - メチル化選択性の統計結果
- `analysis/14_DMG_selectivity_260206_v1/figures/` - 選択性パターン図
