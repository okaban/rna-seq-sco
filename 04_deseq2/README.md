# 04_deseq2

## 概要

DESeq2を用いた差次的発現遺伝子（DEG）解析ステップ。3タイムポイント間（M145_1 vs M145_2 vs M145_3）のペアワイズ比較を実施し、正規化・統計検定・可視化を行う。

## ディレクトリ構造

```
04_deseq2/
├── analysis/
│   └── 04_deseq2_260128_v1/     # DESeq2解析結果
│       ├── results/              # DEGリスト（TSV）
│       ├── figures/              # PCA, MA plot, Volcano plot等
│       ├── rds/                  # RDS中間ファイル
│       └── logs/                 # 実行ログ
├── scripts/
│   ├── run_deseq2_M145.R        # DESeq2メインスクリプト
│   ├── run_deseq2_M145_figures.R # 図表生成スクリプト
│   └── environment_rnaseq_deseq2_260128.yml
└── 04_deseq2_M145_prompt.md     # 実行プロンプト
```

## 主要な出力

- `analysis/04_deseq2_260128_v1/results/` - DEGテーブル（log2FC, padj等）
- `analysis/04_deseq2_260128_v1/figures/` - PCA, MAプロット, ボルケーノプロット
- `analysis/04_deseq2_260128_v1/rds/` - DESeqDataSetオブジェクト（再利用可能）
- `analysis/04_deseq2_260128_v1/deseq2_report_M145.md` - DESeq2解析レポート
