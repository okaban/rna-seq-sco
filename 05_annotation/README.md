# 05_annotation

## 概要

DEGに対する機能アノテーション付与ステップ。*S. coelicolor* A3(2)の遺伝子にGO, KEGG, COGなどの機能情報を紐付ける。

## ディレクトリ構造

```
05_annotation/
├── analysis/
│   └── 05_annotation_260128_v1/ # アノテーション結果
│       ├── tables/               # アノテーションテーブル
│       ├── scripts/              # 実行時生成スクリプト
│       └── logs/                 # 実行ログ
├── scripts/
│   └── run_annotation_M145.py   # アノテーションスクリプト（Python）
└── 05_annotation_M145_prompt.md # 実行プロンプト
```

## 主要な出力

- `analysis/05_annotation_260128_v1/tables/` - 遺伝子機能アノテーションテーブル
- `analysis/05_annotation_260128_v1/annotation_report_M145.md` - アノテーションレポート
