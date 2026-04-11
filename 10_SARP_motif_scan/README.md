# 10_SARP_motif_scan

## 概要

SARPファミリーTFの結合モチーフスキャンステップ。SARP結合モチーフのPWMを構築し、FIMOを用いてプロモーター領域の結合サイト予測を行う。

## ディレクトリ構造

```
10_SARP_motif_scan/
├── analysis/
│   └── 10_SARP_motif_scan_260128_v1/  # モチーフスキャン結果
│       ├── pwm/                        # Position Weight Matrices
│       ├── promoters/                  # プロモーター配列
│       ├── fimo_raw/                   # FIMO生出力
│       ├── figures/                    # モチーフ可視化
│       ├── tables/                     # 結合サイト予測結果
│       └── logs/                       # 実行ログ
├── scripts/
│   ├── run_SARP_motif_scan_M145.R     # モチーフスキャン本体
│   ├── run_SARP_fimo_scan_M145.R      # FIMOスキャン
│   └── run_SARP_aggregate_M145.R      # 結果集約
└── 10_SARP_motif_scan_M145_prompt.md  # 実行プロンプト
```

## 主要な出力

- `analysis/.../fimo_raw/` - FIMO結合サイト予測の生データ
- `analysis/.../tables/` - 統合結合サイトテーブル
- `analysis/.../figures/` - モチーフロゴ・分布図
- `analysis/.../SARP_motif_scan_report_M145.md` - モチーフスキャンレポート
