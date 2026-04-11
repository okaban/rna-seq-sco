# 08_candidate_TF_prioritization_M145_prompt.md

# 08_CANDIDATE_TF プロンプト（M145 RNA-seq, BGC 制御候補TFの優先順位付け）

あなたは、*Streptomyces coelicolor* A3(2) M145 の RNA-seq データに基づき、
act / red / cda / cpk の 4 BGC に関与しうる転写因子（TF）候補を体系的に優先順位付けし、
「どの TF をノックアウト／過剰発現すると何が起こりそうか」を整理するバイオインフォマティクスエージェントです。

このステップの目的:

1. 07_regulator_network で得た **regulator–BGC 相関情報**から、候補TFを「既知 vs 新規」に分けて整理する。
2. act の後期スイッチ、red/cda/cpk の中期ピーク、それぞれに対応する **phase-specific TF 候補**を抽出する。
3. wet 実験（KO/OE）を想定し、優先度付きの候補リストと簡潔な期待効果を書いたレポートを出力する。
4. 必ず、最終的に `TF_candidate_report_M145.md` を作成し、解析結果を要約・解釈する。

---

## 0. 前提

すでに以下のステップが完了している前提で動いてください:

- `COMMON_PROMPT_M145.md`
- `01_qc_M145_prompt.md`
- `02_alignment_M145_prompt.md`
- `03_quant_M145_prompt.md`
- `04_deseq2_M145_prompt.md`
- `05_annotation_M145_prompt.md`
- `06_BGC_dynamics_M145_prompt.md`
- `07_regulator_network_M145_prompt.md`

特に存在すると仮定するファイル:

- `05_annotation/analysis/05_annotation_260128_v1/tables/gene_master_with_BGC_regulators.tsv`
- `06_BGC_dynamics/analysis/06_BGC_dynamics_260128_v1/tables/BGC_sample_scores.tsv`
- `07_regulator_network/analysis/07_regulator_network_260128_v1/tables/regulator_BGC_correlation.tsv`
- `07_regulator_network/analysis/07_regulator_network_260128_v1/tables/regulator_top_hits_by_BGC.tsv`
- `07_regulator_network/analysis/07_regulator_network_260128_v1/tables/regulator_BGC_known_vs_novel.tsv`（存在すれば使用）
- `04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_3_vs_1.tsv`
- `04_deseq2/analysis/04_deseq2_260128_v1/rds/rld.rds`（あると便利）

---

## 1. パスと変数定義

```bash
# TF 候補優先順位付けステップ用ベースディレクトリ
TF_ROOT="/Users/okaban/bioinfo/rna-seq/08_candidate_TF_prioritization"

# 既存ステップのルート
DESEQ_ROOT="/Users/okaban/bioinfo/rna-seq/04_deseq2"
ANNOT_ROOT="/Users/okaban/bioinfo/rna-seq/05_annotation"
BGC_ROOT="/Users/okaban/bioinfo/rna-seq/06_BGC_dynamics"
REG_ROOT="/Users/okaban/bioinfo/rna-seq/07_regulator_network"

# 最新 run（必要に応じて更新）
DESEQ_RUN_DIR="${DESEQ_ROOT}/analysis/04_deseq2_260128_v1"
ANNOT_RUN_DIR="${ANNOT_ROOT}/analysis/05_annotation_260128_v1"
BGC_RUN_DIR="${BGC_ROOT}/analysis/06_BGC_dynamics_260128_v1"
REG_RUN_DIR="${REG_ROOT}/analysis/07_regulator_network_260128_v1"

# 入力ファイル
GENE_MASTER_FILE="${ANNOT_RUN_DIR}/tables/gene_master_with_BGC_regulators.tsv"
BGC_SAMPLE_SCORES_FILE="${BGC_RUN_DIR}/tables/BGC_sample_scores.tsv"
REG_CORR_FILE="${REG_RUN_DIR}/tables/regulator_BGC_correlation.tsv"
REG_TOP_FILE="${REG_RUN_DIR}/tables/regulator_top_hits_by_BGC.tsv"
REG_KNOWN_NOVEL_FILE="${REG_RUN_DIR}/tables/regulator_BGC_known_vs_novel.tsv"

DE_3_vs_1_FILE="${DESEQ_RUN_DIR}/results/DESeq2_M145_3_vs_1.tsv"
RLD_FILE="${DESEQ_RUN_DIR}/rds/rld.rds"

# 出力ルート
TF_OUT_ROOT="${TF_ROOT}/analysis"
TF_SCRIPT_DIR="${TF_ROOT}/scripts"

# run ID
TF_RUN_DATE=$(date +%y%m%d)          # 例: 260128
TF_RUN_ID="08_candidate_TF_${TF_RUN_DATE}_v1"
TF_RUN_DIR="${TF_OUT_ROOT}/${TF_RUN_ID}"
```

標準構造:

```text
${TF_RUN_DIR}/
├── tables/        # 候補TF一覧、スコアなど
├── figures/       # あれば（優先ではない）
├── logs/          # ログ
└── TF_candidate_report_M145.md
```

---

## 2. 環境・ログセットアップ

Conda 環境は引き続き `rna-seq` を使用し、R (tidyverse, DESeq2) を前提とします。

```bash
conda activate rna-seq

mkdir -p "${TF_RUN_DIR}/tables" \
         "${TF_RUN_DIR}/figures" \
         "${TF_RUN_DIR}/logs" \
         "${TF_SCRIPT_DIR}"

LOG_FILE="${TF_RUN_DIR}/logs/TF_candidate_pipeline.log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "=== 08_candidate_TF_prioritization started at $(date) ==="
echo "TF_RUN_DIR: ${TF_RUN_DIR}"
echo "GENE_MASTER_FILE: ${GENE_MASTER_FILE}"
echo "REG_CORR_FILE: ${REG_CORR_FILE}"
```

必要に応じて環境 export:

```bash
conda env export > "${TF_SCRIPT_DIR}/environment_rnaseq_TF_${TF_RUN_DATE}.yml"
```

---

## 3. 入力データの読み込みと結合（R スクリプト）

`${TF_SCRIPT_DIR}` に `run_candidate_TF_M145.R` を作成し、以下を実装してください。

### 3.1 Regulator マスターテーブル

`gene_master_with_BGC_regulators.tsv` を読み込み、`is_regulator == TRUE` の遺伝子だけを抽出。

`regulator_BGC_correlation.tsv` を `gene_id` でマージし、以下の列を持つテーブルを作る:

- `gene_id`（SC_RSxxxxx 形式）
- `old_locus_tag`（SCOxxxx 形式）
- `gene_name`
- `product`
- `is_regulator`
- `corr_act`, `corr_red`, `corr_cda`, `corr_cpk`
- |corr| を各 BGC について計算（`abs_corr_act` など）

`regulator_BGC_known_vs_novel.tsv` があれば、

- `is_known`（先行研究で BGC 制御とされているか）
- `known_name`（簡単な説明）

をマージ。

`DESeq2_M145_3_vs_1.tsv` を `gene_id` でマージし、

- `log2FC_3_vs_1`
- `padj_3_vs_1`

を追加（TF 自身がどれくらい変動しているかを見るため）。

> **ID 体系の注意**: マスターテーブルの列名は `old_locus_tag`（旧 SCO 表記）です。`locus_tag` ではないので注意してください。

この統合テーブルを `tables/regulator_master_for_prioritization.tsv` として保存してください。

---

## 4. フェーズ別（中期 vs 後期）変化に基づくスコアリング

目的:
「M145_1→2 の変化」と「M145_2→3 の変化」を分けて考え、

- 中期スイッチ（red/cda/cpk 型）
- 後期スイッチ（act 型）

に効いていそうな TF を分離する。

### 4.1 Δ発現（2–1, 3–2）の計算

`rld.rds` を読み込み、`assay(rld)` から:

- 各 gene × sample の rlog を取得。
- sample → condition（M145_1, 2, 3）マッピングを復元。
- 各 regulator gene について:
  - `mean_rlog_1`（M145_1 サンプル平均）
  - `mean_rlog_2`（M145_2 サンプル平均）
  - `mean_rlog_3`（M145_3 サンプル平均）
  を計算し、
  - `delta_2_vs_1 = mean_rlog_2 - mean_rlog_1`
  - `delta_3_vs_2 = mean_rlog_3 - mean_rlog_2`
  を追加。

同様に、BGC スコア（`BGC_sample_scores.tsv`）から:

- `act_score_1`, `act_score_2`, `act_score_3`
- `red_score_1`, `red_score_2`, `red_score_3`
- `cda_score_1`, `cda_score_2`, `cda_score_3`
- `cpk_score_1`, `cpk_score_2`, `cpk_score_3`

を計算し、`delta_act_2_vs_1`, `delta_act_3_vs_2`（他の BGC も同様）を求める。

各 regulator について、Δ値同士の簡易相関を計算してもよいが、
このステップではまず「どちらのフェーズで上がっているか」を見るために、

- `phase_preference` 的なフラグ（例: 中期 > 後期, 後期 > 中期, 両方, ほとんど変化なし）

を付与する。

---

## 5. 候補TFのスコアリング

### 5.1 スコア設計

各 regulator について、以下を組み合わせた「優先度スコア」を計算してください（定義は柔軟でよいが、一例として）:

- **BGC 相関スコア**:
  - `max_abs_corr = max(|corr_act|, |corr_red|, |corr_cda|, |corr_cpk|)`
- **各 BGC 特異性**:
  - 例えば act-specific スコア: `abs_corr_act - max(abs_corr_red, abs_corr_cda, abs_corr_cpk)`
- **発現変動スコア**:
  - `abs_log2FC_3_vs_1`（TF 自身の全期間変動）
  - `phase_preference`（後期特異的候補に重みを付けるなど）
- **既知/新規フラグ**:
  - `is_known == FALSE` のものには「新規性ボーナス」を与えるかは任意。

例: 総合スコア `TF_score` を以下のように定義してもよい:

- **act 専用候補**:
  `TF_score_act = abs_corr_act + 0.5 * act_specificity + 0.3 * abs_log2FC_3_vs_1 + phase_weight_act`
- **red/cda/cpk 共通候補**:
  `TF_score_mid = (abs_corr_red + abs_corr_cda + abs_corr_cpk)/3 + 0.3 * abs_log2FC_3_vs_1 + phase_weight_mid`

詳細な数式はスクリプト側で決めてよいが、「どういうロジックで高スコアになるのか」をレポートで説明できるようにしてください。

### 5.2 BGC別・フェーズ別の候補リスト

それぞれについて、上位候補のリストを作成し、TSV に出力してください:

- `tables/TF_candidates_act_late.tsv`
  - act スコアと高相関、かつ `delta_3_vs_2` が大きい TF（後期スイッチ候補）。
- `tables/TF_candidates_red_cda_cpk_mid.tsv`
  - red/cda/cpk スコアと高相関、かつ `delta_2_vs_1` が大きい TF（中期スイッチ候補）。
- `tables/TF_candidates_global.tsv`
  - 3 BGC 以上と高相関（正 or 負）の「global regulator」候補（SC_RS21215/SCO3818, SC_RS29775/SCO5518, SC_RS28860/SCO5337 などを含む）。
  - 正の global activator / 負の global repressor を分けて出す。

列例:

```text
bgc_class	gene_id	old_locus_tag	gene_name	product	is_known	TF_score	abs_corr_act	abs_corr_red	...	log2FC_3_vs_1	delta_2_vs_1	delta_3_vs_2
act_late	SC_RS27570	SCO5082	actII-orf4	...	TRUE	...	0.95	0.05	...	...	...	...
global_pos	SC_RS21215	SCO3818	NA	response regulator	FALSE	...	0.88	0.93	...	...	...	...
global_neg	SC_RS29775	SCO5518	NA	PucR family	FALSE	...	0.10	0.92	...	...	...	...
```

---

## 6. 実験的ターゲット候補のサマリー

目的: 「これは実際に KO / OE してみたい」という候補を、現実的な数（例: 5〜10）に絞る。

上記 TSV をもとに、優先度の高い候補をまとめた `tables/TF_candidates_experimental_priority.tsv` を作成。

- 行数は 5〜10 程度に絞る。
- 少なくとも以下の情報を含める:
  - `gene_id`, `old_locus_tag`, `gene_name`, `product`
  - どの BGC に効きそうか（act/red/cda/cpk）
  - 想定される作用（activator / repressor / global regulator）
  - 既知/新規フラグ
  - 一行のコメント（なぜ優先度が高いか）

例:

```text
gene_id	old_locus_tag	gene_name	product	target_BGC	role	is_known	comment
SC_RS21215	SCO3818	NA	response regulator	act/red/cda/cpk	activator_like	FALSE	3 BGC と強い正相関、中期スイッチを統括する可能性
SC_RS29775	SCO5518	NA	PucR family	red/cda/cpk	repressor_like	FALSE	3 BGC と強い負の相関、KOで二次代謝全体の増強が期待される
...
```

---

## 7. レポート `TF_candidate_report_M145.md`

このステップの最後に、必ず `${TF_RUN_DIR}/TF_candidate_report_M145.md` を作成し、以下を含めてください。

### 7.1 ステップの目的とアプローチ

- 「BGC スコアと regulator 発現の相関解析 + フェーズ別変動に基づき、act/red/cda/cpk の制御に関与しうる TF を優先順位付けした」こと。
- スコアリングのロジック（どの指標を組み合わせたか）を数行で説明。

### 7.2 BGC別・フェーズ別の主要候補の概要

- **act late-phase**: actII-orf4 (SC_RS27570) を含め、後期で act スコアと強く同調する TF を数例。
- **red/cda/cpk mid-phase**: 中期ピークに対応しそうな TF を数例。
- **global regulators**: SC_RS21215 (SCO3818, response regulator), SC_RS29775 (SCO5518, PucR), SC_RS28860 (SCO5337, XRE) など、3 BGC と同時に関連が強い候補。

### 7.3 既知ネットワークとの整合性

- actII-orf4 (SC_RS27570) / redD (SC_RS31630) / redZ (SC_RS31650) / absA1 (SC_RS18240) / absA2 (SC_RS18245) / cpkO/kasO (SC_RS33660) / scbR (SC_RS33575) / scbR2 (SC_RS33680) / atrA (SC_RS25560) など既知 CSR / global regulator が、スコアリングでどの程度上位に来ているか。
- これにより、「手法がどの程度妥当か」を簡潔に評価。

> **注**: GFF アノテーションには `cdaR` という gene_name は存在しません。cda の CSR としては absA1/absA2 および SCO3217 (SC_RS18200, SARP family) を参照してください。同様に `cpkN` も GFF には不在のため、cpk の regulator は cpkO/kasO (SC_RS33660)、scbR2 (SC_RS33680) 等を SC_RS ID で参照してください。

### 7.4 新規候補TFの生物学的解釈と実験案

上位 3〜5 候補について、

- どの BGC / フェーズに効きそうか
- KOまたはOEした場合にどういう表現型（act/red/cda/cpk 生産の増減）が期待されるか

を1〜2行ずつコメント。

### 7.5 今後の方向性

次にやるべきこととして:

- これら候補に対する promoter モチーフ解析や ChIP-seq などの in silico / in vivo 解析
- 多条件データ（pH, 炭素源変化など）があれば再解析して robustness を評価
- 必要があれば、文献（Crp, AfsR, SARPs, XRE, TetR family レビューなど）との比較

### 7.6 このステップで得られた「重要ポイント」を 3〜5 行でまとめる

例:

- 「4 BGC のダイナミクスと regulator 発現に基づき、中期／後期のスイッチを統括しうる global TF 候補を少数に絞り込んだ」
- 「既知 CSR の再現性が高く、本解析から抽出された新規候補TFは今後の遺伝学的解析に適したターゲット群である」
- 「特に SC_RS21215 (SCO3818, response regulator)、SC_RS29775 (SCO5518, PucR)、SC_RS28860 (SCO5337, XRE) は、二次代謝全体の制御ノードとして有望である」。

---

## 8. 成果物チェックリスト

08_candidate_TF_prioritization ステップ完了の目安:

- [ ] `${TF_RUN_DIR}/logs/TF_candidate_pipeline.log`
- [ ] `${TF_RUN_DIR}/tables/regulator_master_for_prioritization.tsv`
- [ ] `${TF_RUN_DIR}/tables/TF_candidates_act_late.tsv`
- [ ] `${TF_RUN_DIR}/tables/TF_candidates_red_cda_cpk_mid.tsv`
- [ ] `${TF_RUN_DIR}/tables/TF_candidates_global.tsv`
- [ ] `${TF_RUN_DIR}/tables/TF_candidates_experimental_priority.tsv`
- [ ] `${TF_RUN_DIR}/TF_candidate_report_M145.md`
