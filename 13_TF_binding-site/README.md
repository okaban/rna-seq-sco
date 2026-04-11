# 13_TF_binding-site

## 概要

全転写因子（TF）のマスターリスト作成と、TF結合サイトにおけるDNAメチル化解析ステップ。文献・データベース情報を統合してTFリストを構築し、結合サイトのメチル化状態を評価する。

## ディレクトリ構造

```
13_TF_binding-site/
├── analysis/
│   ├── 01_master_TF_list_260206_v1/      # TFマスターリスト作成
│   │   ├── data/                          # 取得データ
│   │   ├── intermediate/                  # 中間ファイル
│   │   ├── literature/                    # 文献情報
│   │   ├── reports/                       # 解析レポート
│   │   └── scripts/                       # 実行スクリプト
│   └── 02_TF_BS_methylation_260207_v1/   # TF結合サイトメチル化解析
│       ├── figures/                       # メチル化可視化
│       ├── tables/                        # メチル化統計テーブル
│       ├── reports/                       # 解析レポート
│       └── scripts/                       # 実行スクリプト
└── 260206_make-list-TF.md                # TFリスト作成プロンプト
```

## 主要な出力

- `analysis/01_master_TF_list_260206_v1/master_TF_list_M145*.tsv` - TFマスターリスト
- `analysis/01_master_TF_list_260206_v1/master_TF_binding_sites_M145*.tsv` - TF結合サイトリスト
- `analysis/02_TF_BS_methylation_260207_v1/tables/` - 結合サイトメチル化統計
- `analysis/02_TF_BS_methylation_260207_v1/figures/` - メチル化パターン図
