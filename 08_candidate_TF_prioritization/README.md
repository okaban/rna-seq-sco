# 08_candidate_TF_prioritization

## 概要

候補転写因子（TF）の優先順位付けステップ。DEG解析・調節ネットワーク情報を統合し、BGC制御に関与する可能性の高いTFを絞り込む。

## ディレクトリ構造

```
08_candidate_TF_prioritization/
├── analysis/
│   └── 08_candidate_TF_260128_v1/        # TF優先順位付け結果
│       ├── figures/                        # TFランキング可視化
│       ├── tables/                         # 候補TFスコアテーブル
│       └── logs/                           # 実行ログ
├── scripts/
│   └── run_candidate_TF_M145.R            # TF優先順位付けスクリプト
└── 08_candidate_TF_prioritization_M145_prompt.md  # 実行プロンプト
```

## 主要な出力

- `analysis/08_candidate_TF_260128_v1/figures/` - TF候補ランキング図
- `analysis/08_candidate_TF_260128_v1/tables/` - TF候補スコアテーブル
- `analysis/08_candidate_TF_260128_v1/TF_candidate_report_M145.md` - TF候補レポート
