# 02_alignment

## 概要

トリミング済みリードを *S. coelicolor* A3(2) 参照ゲノムにアラインメントするステップ。HISAT2およびSTARのインデックスを作成し、マッピングを実行する。

## ディレクトリ構造

```
02_alignment/
├── analysis/
│   └── 02_alignment_260127_v1/  # アラインメント結果
│       ├── bam/                  # ソート済みBAMファイル
│       ├── fastq_decompressed/   # 解凍済みFASTQ
│       ├── summary/              # マッピング率サマリー
│       └── logs/                 # 実行ログ
├── hisat2_index/                 # HISAT2ゲノムインデックス
├── star_index/                   # STARゲノムインデックス
├── scripts/
│   ├── environment_rnaseq_alignment_260127.yml
│   └── ungz.sh                   # FASTQ解凍スクリプト
└── 02_alignment_M145_prompt.md   # 実行プロンプト
```

## 主要な出力

- `analysis/02_alignment_260127_v1/bam/` - ソート済みBAMファイル（全9サンプル）
- `analysis/02_alignment_260127_v1/summary/` - マッピング率統計
- `analysis/02_alignment_260127_v1/alignment_report_M145.md` - アラインメントレポート
