# 09_SARP_integration

## 概要

SARPファミリー転写因子の統合解析ステップ。HMMプロファイル（PF00486, PF00931, PF03704, PF13424）を用いたSARPドメイン検索とRNA-seq発現データの統合を行う。

## ディレクトリ構造

```
09_SARP_integration/
├── analysis/
│   ├── 09_SARP_integration_260128_v1/  # v1解析結果
│   └── 09_SARP_integration_260128_v2/  # v2解析結果（改良版）
│       ├── figures/                     # SARP発現パターン図
│       ├── tables/                      # SARPリスト・発現テーブル
│       └── logs/                        # 実行ログ
├── scripts/
│   ├── run_SARP_integration_M145.R     # v1スクリプト
│   ├── run_SARP_integration_v2_M145.R  # v2スクリプト
│   ├── SARP_domains_combined.hmm*      # 統合HMMプロファイル
│   └── PF*.hmm                         # 個別Pfamドメイン
└── 09_SARP_integration_M145_prompt*.md # 実行プロンプト
```

## 主要な出力

- `analysis/09_SARP_integration_260128_v2/tables/` - SARP同定結果・発現量テーブル
- `analysis/09_SARP_integration_260128_v2/figures/` - SARP発現動態図
- `analysis/09_SARP_integration_260128_v2/SARP_report_M145_v2.md` - SARPレポート
