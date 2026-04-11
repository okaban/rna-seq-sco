# 01_qc

## 概要

RNA-seqリードの品質管理ステップ。FastQC（生リード）、fastp（トリミング・フィルタリング）、MultiQC（統合レポート）を実行する。

## ディレクトリ構造

```
01_qc/
├── analysis/
│   └── 01_qc_260127_v1/       # QC実行結果
│       ├── raw_fastqc/         # 生リードのFastQC結果
│       ├── fastp/              # トリミング後のリード・レポート
│       ├── trimmed_fastqc/     # トリミング後のFastQC結果
│       ├── multiqc/            # MultiQC統合レポート
│       ├── figures/            # QC可視化図
│       ├── logs/               # 実行ログ
│       └── sample_table.tsv    # サンプル情報テーブル
├── scripts/
│   └── environment_rnaseq_260127.yml  # Conda環境定義
└── 01_qc_M145_check_prompt.md  # QC検証用プロンプト
```

## 主要な出力

- `analysis/01_qc_260127_v1/multiqc/` - 全サンプルの統合QCレポート
- `analysis/01_qc_260127_v1/fastp/` - トリミング済みリード（FASTQ）
- `analysis/01_qc_260127_v1/qc_report_M145.md` - QCサマリーレポート
- `analysis/01_qc_260127_v1/sample_table.tsv` - サンプル管理テーブル（M145_1/2/3, 各n=3）
