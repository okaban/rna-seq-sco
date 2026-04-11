# 07_regulator_network_M145_prompt.md

# 07_REGULATOR_NETWORK プロンプト（M145 RNA-seq, BGC と制御因子の共変動解析）

あなたは、*Streptomyces coelicolor* A3(2) M145 の RNA-seq データを用いて、
主要4 BGC（act, red, cda, cpk）の発現と、ゲノム全体の regulator 遺伝子（is_regulator=TRUE）との共変動を解析し、
「どの regulator がどの BGC の ON/OFF に関与していそうか」を定量的に示すバイオインフォマティクスエージェントです。

このステップの目的:

1. 877 個の regulator 遺伝子と、4 BGC スコア（act/red/cda/cpk）の **相関構造** を解析する。
2. 各 BGC と高相関な regulator 候補を **ランキングと図**として提示する。
3. 既知の制御ネットワーク（ScbR/ScbR2, cpkO/kasO, actII-orf4, redD/redZ, absA1/absA2, AfsQ1/Q2, DraR, DasR など）との整合性を簡単にコメントする。
4. 最後に、結果を要約・解釈した Markdown レポート `regulator_network_report_M145.md` を必ず出力する。

> **文献背景**: ScbR/ScbR2 が gamma-butyrolactone シグナル系を介して cpk を正に・act/red を間接的に負に制御する報告や、cpkO/kasO が coelimycin 活性化と他抗生物質への前駆体配分に関与する報告など、複数 BGC 間のクロストークが知られている。

---

## 0. 前提

以下のステップが完了している前提で動いてください:

- `COMMON_PROMPT_M145.md`
- `01_qc_M145_prompt.md`
- `02_alignment_M145_prompt.md`
- `03_quant_M145_prompt.md`
- `04_deseq2_M145_prompt.md`
- `05_annotation_M145_prompt.md`（regulator フラグ付きマスターテーブル）
- `06_BGC_dynamics_M145_prompt.md`（BGC スコアと BGC ダイナミクス）

特に存在すると仮定するファイル:

- `04_deseq2/analysis/04_deseq2_260128_v1/results/normalized_counts_M145.tsv`
- `04_deseq2/analysis/04_deseq2_260128_v1/rds/rld.rds`（rlog オブジェクト）
- `05_annotation/analysis/05_annotation_260128_v1/tables/gene_master_with_BGC_regulators.tsv`
- `06_BGC_dynamics/analysis/06_BGC_dynamics_260128_v1/tables/BGC_sample_scores.tsv`
- `06_BGC_dynamics/analysis/06_BGC_dynamics_260128_v1/BGC_dynamics_report_M145.md`

---

## 1. パスと変数定義

```bash
# Regulator ネットワーク解析ステップ用ベースディレクトリ
REG_ROOT="/Users/okaban/bioinfo/rna-seq/07_regulator_network"

# 既存ステップのルート
DESEQ_ROOT="/Users/okaban/bioinfo/rna-seq/04_deseq2"
ANNOT_ROOT="/Users/okaban/bioinfo/rna-seq/05_annotation"
BGC_ROOT="/Users/okaban/bioinfo/rna-seq/06_BGC_dynamics"

# 最新 run（必要に応じて更新）
DESEQ_RUN_DIR="${DESEQ_ROOT}/analysis/04_deseq2_260128_v1"
ANNOT_RUN_DIR="${ANNOT_ROOT}/analysis/05_annotation_260128_v1"
BGC_RUN_DIR="${BGC_ROOT}/analysis/06_BGC_dynamics_260128_v1"

# 入力ファイル
NORMALIZED_COUNTS_FILE="${DESEQ_RUN_DIR}/results/normalized_counts_M145.tsv"
RLD_FILE="${DESEQ_RUN_DIR}/rds/rld.rds"
GENE_MASTER_FILE="${ANNOT_RUN_DIR}/tables/gene_master_with_BGC_regulators.tsv"
BGC_SAMPLE_SCORES_FILE="${BGC_RUN_DIR}/tables/BGC_sample_scores.tsv"

# 出力ルート
REG_OUT_ROOT="${REG_ROOT}/analysis"
REG_SCRIPT_DIR="${REG_ROOT}/scripts"

# run ID
REG_RUN_DATE=$(date +%y%m%d)         # 例: 260128
REG_RUN_ID="07_regulator_network_${REG_RUN_DATE}_v1"
REG_RUN_DIR="${REG_OUT_ROOT}/${REG_RUN_ID}"
```

標準ディレクトリ構造:

```text
${REG_RUN_DIR}/
├── tables/        # 相関行列、ランキングなど
├── figures/       # 散布図・ヒートマップなど
├── logs/          # ログ
└── regulator_network_report_M145.md
```

---

## 2. 環境・ログセットアップ

Conda 環境は引き続き `rna-seq` を使用し、R (DESeq2, tidyverse, ggplot2, pheatmap 等) が利用可能とします。

```bash
conda activate rna-seq

mkdir -p "${REG_RUN_DIR}/tables" \
         "${REG_RUN_DIR}/figures" \
         "${REG_RUN_DIR}/logs" \
         "${REG_SCRIPT_DIR}"

LOG_FILE="${REG_RUN_DIR}/logs/regulator_network_pipeline.log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "=== 07_regulator_network started at $(date) ==="
echo "REG_RUN_DIR: ${REG_RUN_DIR}"
echo "GENE_MASTER_FILE: ${GENE_MASTER_FILE}"
echo "BGC_SAMPLE_SCORES_FILE: ${BGC_SAMPLE_SCORES_FILE}"
```

必要に応じて環境 export:

```bash
conda env export > "${REG_SCRIPT_DIR}/environment_rnaseq_regnet_${REG_RUN_DATE}.yml"
```

---

## 3. 入力データの読み込みと準備（R スクリプト）

`${REG_SCRIPT_DIR}` に `run_regulator_network_M145.R` を作成し、以下を実装してください（ここではやるべき処理の仕様だけ書きます）。

### 3.1 rlog / 正規化カウントとサンプル情報

R スクリプト内で:

- `rld.rds` を読み込み、rld オブジェクトを得る。
  - `assay(rld)` は gene × sample の rlog 変換行列。
- rlog の `colnames`（サンプル名）から、condition（M145_1, M145_2, M145_3）と replicate を再構築するか、
  既存の `quant_sample_table.tsv` があれば再利用して `colData(rld)` と整合させる。
- `BGC_sample_scores.tsv` を読み込み、各サンプルについて:
  - `act_score`, `red_score`, `cda_score`, `cpk_score`
  が得られるようにする。

### 3.2 Regulator 遺伝子リスト

`gene_master_with_BGC_regulators.tsv` を読み込み、`is_regulator == TRUE` の遺伝子を抽出。

少なくとも以下の列を保持:

- `gene_id`（SC_RSxxxxx 形式）
- `old_locus_tag`（SCOxxxx 形式）
- `gene_name`
- `product`
- `is_regulator`
- 可能なら既知の regulator 種類（SARP, two-component, sigma, GBL receptor など）

> **ID 体系の注意**: マスターテーブルの列名は `old_locus_tag`（旧 SCO 表記）です。`locus_tag` ではないので注意してください。

抽出された regulator 遺伝子（約 877 個）について、`assay(rld)` から対応する発現行列（regulator × sample）をサブセットする。

---

## 4. BGC スコアと Regulator 発現の共変動解析

### 4.1 Regulator × BGC 相関行列の計算

各 regulator 遺伝子について、次を計算してください:

- サンプルごとの rlog 発現ベクトル（長さ 9）
- BGC スコアベクトル（act/red/cda/cpk, 各長さ 9）

それぞれについて、Pearson 相関係数を求める:

- `corr(regulator, act_score)`
- `corr(regulator, red_score)`
- `corr(regulator, cda_score)`
- `corr(regulator, cpk_score)`

結果を `tables/regulator_BGC_correlation.tsv` に保存。

形式例:

```text
gene_id	old_locus_tag	gene_name	product	corr_act	corr_red	corr_cda	corr_cpk
SC_RS33660	SCO6282	cpkO	...	0.95	0.80	0.60	0.99
...
```

相関の絶対値（|corr|）も別列で持ってよいです。

### 4.2 高相関 regulator のランキング

各 BGC ごとに:

- |corr| で降順ソートし、上位 N（例: 上位30）regulator を抽出。
- `tables/regulator_top_hits_by_BGC.tsv` として保存。

形式例:

```text
bgc_name	gene_id	old_locus_tag	gene_name	product	corr	sign	corr_rank	is_known
act	SC_RSxxxxx	SCO5800	...	...	0.93	+	1	TRUE
act	SC_RSyyyyy	SCOyyyy	...	...	0.91	+	2	FALSE
...
cpk	SC_RS33660	SCO6282	cpkO	...	0.99	+	1	TRUE
...
```

- `sign` は `+`（正の相関）か `-`（負の相関）。
- `is_known` は、文献で BGC 制御に関与すると知られている regulator なら `TRUE`、それ以外は `FALSE`。

known list（手動定義してよい）:

| regulator_name | gene_id | old_locus_tag | 主な制御対象 |
|---|---|---|---|
| actII-orf4 | SC_RS27570 | SCO5082 | act CSR (SARP) |
| redD | SC_RS31630 | SCO5877 | red CSR (SARP) |
| redZ | SC_RS31650 | SCO5881 | red (response regulator) |
| absA1 | SC_RS18240 | SCO3225 | cda CSR (sensor kinase), pleiotropic |
| absA2 | SC_RS18245 | SCO3226 | cda CSR (response regulator), pleiotropic |
| cpkO/kasO | SC_RS33660 | SCO6282 | cpk CSR (SARP-like) |
| scbR | SC_RS33575 | SCO6265 | gamma-butyrolactone receptor |
| scbR2 | SC_RS33680 | SCO6286 | gamma-butyrolactone receptor (cpk 内) |
| atrA | SC_RS25560 | SCO4677 | TetR family, actII-orf4 上位 |

> **注**: GFF アノテーションには `cdaR` という gene_name は存在しません。cda クラスター内の regulator としては absA1 (SC_RS18240) / absA2 (SC_RS18245) および SCO3217 (SC_RS18200, SARP family auto-detected) が該当します。`cpkN` も GFF には明示的な gene_name としては存在しないため、cpk クラスター内の regulator は cpkO/kasO (SC_RS33660)、scbR2 (SC_RS33680)、SC_RS33650 (SCO6280, SARP auto-detected)、SC_RS33690 (SCO6288, SARP auto-detected) を SC_RS ID で参照してください。

---

## 5. 代表的 Regulator のタイムコースと BGC スコアの図示

### 5.1 BGC スコア vs Regulator 発現の散布図

各 BGC について、代表的な regulator（既知＋新規候補）を 2〜3 個選び、
サンプルごとの BGC スコアと rlog 発現の散布図を描きます。

- x軸: BGC スコア（例: cpk_score）
- y軸: regulator rlog 発現（例: cpkO）
- ポイント: サンプル（条件で色分け）

出力例:

```text
${REG_RUN_DIR}/figures/scatter_cpk_score_vs_cpkO_M145.pdf
${REG_RUN_DIR}/figures/scatter_cpk_score_vs_scbR2_M145.pdf
...
```

### 5.2 タイムコースライン（Regulator + BGC スコア）

各 BGC に対して、代表的 regulator 1〜2 個について:

- 条件ごとの平均発現（M145_1, 2, 3）を計算し、
- BGC スコアと同じ図に二軸 or 正規化して重ねて描く。

例:

- 上図: cpk_score（M145_1→2→3）
- 下図: cpkO/kasO (SC_RS33660), scbR2 (SC_RS33680) の rlog 平均値（同じ x 軸）

出力例:

```text
${REG_RUN_DIR}/figures/timecourse_cpk_score_with_cpkO_scbR2_M145.pdf
${REG_RUN_DIR}/figures/timecourse_act_score_with_actII_orf4_M145.pdf
${REG_RUN_DIR}/figures/timecourse_red_score_with_redD_redZ_M145.pdf
${REG_RUN_DIR}/figures/timecourse_cda_score_with_absA1_absA2_M145.pdf
```

これにより:

- 「cpkO/kasO の立ち上がり → cpk スコアの立ち上がり」の関係
- 「ScbR/ScbR2 のタイミングと act/red/cda/cpk の切り替わり」

などを視覚的に確認できます。

---

## 6. 既知ネットワークとの比較（定性的なチェック）

`regulator_top_hits_by_BGC.tsv` と文献情報を用いて、少なくとも以下を確認し、後のレポートでコメントする準備をしてください:

- cpk スコアと cpkO/kasO (SC_RS33660) の相関は高いか（期待通り）。
- act スコアと actII-orf4 (SC_RS27570) の相関、red スコアと redD (SC_RS31630) / redZ (SC_RS31650) の相関、cda スコアと absA1 (SC_RS18240) / absA2 (SC_RS18245) の相関が高いか。
- ScbR (SC_RS33575) / ScbR2 (SC_RS33680) が cpk スコアと正の/負の相関を持ち、act/red/cda とは逆方向の傾向を持つか。

これら既知のパターンが再現できているなら、「解析パイプラインの妥当性確認」として重要。

また、新規候補についても:

- act スコアに高相関だが、これまで BGC 制御に直接関与すると報告されていない TF
- 複数 BGC スコア（例えば act と cda）に同時に高相関な global regulator

などを、候補としてフラグを立てておくとよいです。

---

## 7. テーブル出力まとめ

少なくとも以下の TSV を出力してください:

- `tables/regulator_BGC_correlation.tsv`
  - 全 regulator × 4 BGC の相関係数。
- `tables/regulator_top_hits_by_BGC.tsv`
  - 各 BGC ごとに上位 N（例: 30）regulator のランキング。
- （任意）`tables/regulator_BGC_known_vs_novel.tsv`
  - `is_known` フラグ付きで、「既知ネットワークを再現している regulator」と「新規候補」を整理したもの。

---

## 8. レポート `regulator_network_report_M145.md`

このステップの最後に、必ず `${REG_RUN_DIR}/regulator_network_report_M145.md` を作成し、結果を要約・解釈してください。

含めるべき内容:

### 8.1 このステップの目的

「主要4 BGC の発現と regulator 群の共変動から、BGC 制御ネットワークの候補リンクを抽出した」。

### 8.2 解析の概要

- 使用した入力（rlog, BGC_sample_scores, gene_master_with_BGC_regulators）
- 相関解析の方法（Pearson, n=9サンプル）
- 上位ヒットの抽出基準（例: |corr|>0.9 などあれば）

### 8.3 主な結果（既知ネットワークの再現）

- act スコアと actII-orf4 (SC_RS27570) の強い正相関の有無
- red スコアと redD (SC_RS31630) / redZ (SC_RS31650) の相関の有無
- cda スコアと absA1 (SC_RS18240) / absA2 (SC_RS18245) の相関の有無
- cpk スコアと cpkO/kasO (SC_RS33660) の相関の有無
- ScbR/ScbR2 や AfsQ1/Q2, DraR/S, DasR の BGC スコアとの相関パターンが、先行研究の知見（cpk 過剰で ACT/RED/CDA が抑制されるなど）とどう対応するか。

### 8.4 主な結果（新規候補）

- 各 BGC について、新規に見えてきた高相関 regulator（gene_id, product）を数例挙げ、その発現パターンを短くコメント。
- 複数 BGC と同時に相関する global regulator 候補。

### 8.5 生物学的解釈と今後の実験候補

- どの regulator をノックアウト／過剰発現すれば、act/red/cda/cpk の生産に影響を与えそうか。
- 特に coelimycin（cpk）と他抗生物質のトレードオフを調整する制御ノード（例: cpkO/kasO, ScbR2）について簡潔に触れる。

### 8.6 このステップで見えた重要なポイントを 3〜5 行でまとめる

例:

- 「BGC スコアと既知 CSR（actII-orf4, redD, absA1/absA2, cpkO/kasO）の相関パターンは先行知見とよく整合し、本データセットの転写ネットワークが教科書的挙動を再現している」
- 「一方で、act と cda を同時に制御していそうな新規 TF 候補や、cpk と逆相関する global regulator 候補が複数見出され、二次代謝フラックス制御の新しいノードとなりうる」
- 「これらの候補は、今後の遺伝学的解析や代謝工学設計の出発点として有望である」。

---

## 9. 成果物チェックリスト

07_regulator_network ステップ完了の目安:

- [ ] `${REG_RUN_DIR}/logs/regulator_network_pipeline.log`
- [ ] `${REG_RUN_DIR}/tables/regulator_BGC_correlation.tsv`
- [ ] `${REG_RUN_DIR}/tables/regulator_top_hits_by_BGC.tsv`
- [ ] （任意）`${REG_RUN_DIR}/tables/regulator_BGC_known_vs_novel.tsv`
- [ ] `${REG_RUN_DIR}/figures/scatter_*_M145.pdf`（BGCスコア vs regulator 散布図）
- [ ] `${REG_RUN_DIR}/figures/timecourse_*_M145.pdf`（BGCスコア＋regulator タイムコース）
- [ ] `${REG_RUN_DIR}/regulator_network_report_M145.md`
