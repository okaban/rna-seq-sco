# 03_quantification

## 概要

BAMファイルから遺伝子ごとのリードカウントを定量するステップ。featureCountsなどを用いてカウントマトリクスを生成する。

## ディレクトリ構造

```
03_quantification/
├── analysis/
│   └── 03_quant_260128_v1/      # 定量結果
│       ├── counts/               # カウントマトリクス
│       ├── figures/              # 定量QC図
│       ├── summary/              # 定量サマリー統計
│       └── logs/                 # 実行ログ
├── scripts/
│   └── environment_rnaseq_quant_260128.yml
└── 03_quant_M145_prompt.md       # 実行プロンプト
```

## 主要な出力

- `analysis/03_quant_260128_v1/counts/` - 遺伝子カウントマトリクス（DESeq2入力用）
- `analysis/03_quant_260128_v1/summary/` - アサインメント率サマリー
- `analysis/03_quant_260128_v1/quant_report_M145.md` - 定量レポート
