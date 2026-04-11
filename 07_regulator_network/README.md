# 07_regulator_network

## 概要

転写調節因子ネットワーク解析ステップ。DEGに含まれる転写因子・シグマ因子・二成分制御系等の制御関係を整理し、調節ネットワークを構築する。

## ディレクトリ構造

```
07_regulator_network/
├── analysis/
│   └── 07_regulator_network_260128_v1/  # ネットワーク解析結果
│       ├── figures/                      # ネットワーク可視化図
│       ├── tables/                       # 調節因子リスト・エッジリスト
│       └── logs/                         # 実行ログ
├── scripts/
│   └── run_regulator_network_M145.R     # ネットワーク構築スクリプト
└── 07_regulator_network_M145_prompt.md  # 実行プロンプト
```

## 主要な出力

- `analysis/07_regulator_network_260128_v1/figures/` - 調節ネットワーク図
- `analysis/07_regulator_network_260128_v1/tables/` - 調節因子・ターゲット遺伝子テーブル
- `analysis/07_regulator_network_260128_v1/regulator_network_report_M145.md` - ネットワーク解析レポート
