# 06_BGC_dynamics

## 概要

生合成遺伝子クラスター（BGC）の発現動態解析ステップ。antiSMASH予測BGCの時系列発現変動を可視化・解析し、二次代謝産物生合成の活性化タイミングを評価する。

## ディレクトリ構造

```
06_BGC_dynamics/
├── analysis/
│   └── 06_BGC_dynamics_260128_v1/  # BGC動態解析結果
│       ├── figures/                 # BGC発現ヒートマップ等
│       ├── tables/                  # BGCごとの発現統計
│       └── logs/                    # 実行ログ
├── scripts/
│   └── run_BGC_dynamics_M145.R     # BGC動態解析スクリプト
└── 06_BGC_dynamics_M145_prompt.md  # 実行プロンプト
```

## 主要な出力

- `analysis/06_BGC_dynamics_260128_v1/figures/` - BGC発現パターン図
- `analysis/06_BGC_dynamics_260128_v1/tables/` - BGC発現量テーブル
- `analysis/06_BGC_dynamics_260128_v1/BGC_dynamics_report_M145.md` - BGC動態レポート
