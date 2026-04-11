# 09_SARP_integration_M145_prompt.md

# 09_SARP_INTEGRATION プロンプト（M145 RNA-seq, SARPファミリー統合解析）

あなたは、*Streptomyces coelicolor* A3(2) M145 の RNA-seq 解析において、
SARPファミリー転写因子（Streptomyces Antibiotic Regulatory Proteins）を **網羅的に同定・分類**し、
既に構築済みの BGC・regulator ネットワークに統合するバイオインフォマティクスエージェントです。

このステップの目的:

1. M145 ゲノム中の **全 SARPファミリータンパク質** を同定し、small / medium / large / SARP-LAL の4サブタイプに分類する。
2. 05〜08で構築したマスターテーブル・相関結果に `is_SARP` と `SARP_subtype` を付与し、
   BGC 発現・TFネットワークの中での SARP の位置付けを明確にする。
3. act/red/cda/cpk を含む既知 BGC と、SARP との関係（cluster-situated vs global）を整理し、
   SARPにフォーカスしたサマリレポート `SARP_report_M145.md` を作成する。

最後に、解析結果の要約と解釈を含む **Markdown レポート** を必ず出力してください。

---

## 0. 前提

以下のステップが完了している前提で動いてください:

- `COMMON_PROMPT_M145.md`
- `01_qc_M145_prompt.md`
- `02_alignment_M145_prompt.md`
- `03_quant_M145_prompt.md`
- `04_deseq2_M145_prompt.md`
- `05_annotation_M145_prompt.md`
- `06_BGC_dynamics_M145_prompt.md`
- `07_regulator_network_M145_prompt.md`
- `08_candidate_TF_prioritization_M145_prompt.md`

特に存在すると仮定するファイル:

- `05_annotation/analysis/05_annotation_260128_v1/tables/gene_master_with_BGC_regulators.tsv`
- `05_annotation/analysis/05_annotation_260128_v1/tables/BGC_definition_manual.tsv`
- `06_BGC_dynamics/analysis/06_BGC_dynamics_260128_v1/tables/BGC_sample_scores.tsv`
- `07_regulator_network/analysis/07_regulator_network_260128_v1/tables/regulator_BGC_correlation.tsv`
- `08_candidate_TF_prioritization/analysis/08_candidate_TF_260128_v1/tables/regulator_master_for_prioritization.tsv`
- M145 リファレンスプロテオーム:
  - `/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/protein.faa`

---

## 1. パスと変数定義

```bash
# SARP 統合解析ステップ用ベースディレクトリ
SARP_ROOT="/Users/okaban/bioinfo/rna-seq/09_SARP_integration"

# 既存ステップのルート
ANNOT_ROOT="/Users/okaban/bioinfo/rna-seq/05_annotation"
BGC_ROOT="/Users/okaban/bioinfo/rna-seq/06_BGC_dynamics"
REG_ROOT="/Users/okaban/bioinfo/rna-seq/07_regulator_network"
TF_ROOT="/Users/okaban/bioinfo/rna-seq/08_candidate_TF_prioritization"

REF_DIR="/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1"
REF_PROT="${REF_DIR}/protein.faa"

# 最新 run（必要に応じて更新）
ANNOT_RUN_DIR="${ANNOT_ROOT}/analysis/05_annotation_260128_v1"
BGC_RUN_DIR="${BGC_ROOT}/analysis/06_BGC_dynamics_260128_v1"
REG_RUN_DIR="${REG_ROOT}/analysis/07_regulator_network_260128_v1"
TF_RUN_DIR="${TF_ROOT}/analysis/08_candidate_TF_260128_v1"

# 入力ファイル
GENE_MASTER_FILE="${ANNOT_RUN_DIR}/tables/gene_master_with_BGC_regulators.tsv"
BGC_DEF_MANUAL_FILE="${ANNOT_RUN_DIR}/tables/BGC_definition_manual.tsv"
BGC_SAMPLE_SCORES_FILE="${BGC_RUN_DIR}/tables/BGC_sample_scores.tsv"
REG_CORR_FILE="${REG_RUN_DIR}/tables/regulator_BGC_correlation.tsv"
TF_REG_MASTER_FILE="${TF_RUN_DIR}/tables/regulator_master_for_prioritization.tsv"

# 出力ルート
SARP_OUT_ROOT="${SARP_ROOT}/analysis"
SARP_SCRIPT_DIR="${SARP_ROOT}/scripts"

# run ID
SARP_RUN_DATE=$(date +%y%m%d)        # 例: 260128
SARP_RUN_ID="09_SARP_integration_${SARP_RUN_DATE}_v1"
SARP_RUN_DIR="${SARP_OUT_ROOT}/${SARP_RUN_ID}"
```

標準構造:

```text
${SARP_RUN_DIR}/
├── tables/        # SARPリスト、サマリ
├── figures/       # SARP×BGCヒートマップなど（任意）
├── logs/          # ログ
└── SARP_report_M145.md
```

---

## 2. 環境・ログセットアップ

Conda 環境は `rna-seq` を使用し、R (tidyverse) と HMMER または InterPro/Pfam 検索用ツールが利用可能とします。

```bash
conda activate rna-seq

mkdir -p "${SARP_RUN_DIR}/tables" \
         "${SARP_RUN_DIR}/figures" \
         "${SARP_RUN_DIR}/logs" \
         "${SARP_SCRIPT_DIR}"

LOG_FILE="${SARP_RUN_DIR}/logs/SARP_integration_pipeline.log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "=== 09_SARP_integration started at $(date) ==="
echo "SARP_RUN_DIR: ${SARP_RUN_DIR}"
echo "REF_PROT: ${REF_PROT}"
echo "GENE_MASTER_FILE: ${GENE_MASTER_FILE}"
```

環境情報:

```bash
conda env export > "${SARP_SCRIPT_DIR}/environment_rnaseq_SARP_${SARP_RUN_DATE}.yml"
```

---

## 3. SARPドメインプロファイルの準備

SARPs に共通するのは、「N末端 OmpR型 DNA-binding domain（DBD）＋ C末端 bacterial transcriptional activation domain（BTAD）」から成る SARP ドメインです。

Pfam などで SARP 関連として使われている代表的プロファイルの例:

- **PF00486**（Trans_reg_C / HTH系 DBD、SARPのN末端領域）
- **PF03704**（BTAD = Bacterial Transcriptional Activator Domain）
- **PF00931**（NB-ARC、medium/large SARP で付加）
- **PF13424**（TPR_12、large SARP で付加）

ここでは簡略化のために:

- DBD: PF00486
- BTAD: PF03704

を最低限検出し、「少なくとも DBD + BTAD を持つもの」を SARP 候補とみなします。
可能なら NB-ARC（PF00931）と TPR（PF13424）も見て、サブタイプ分類に使います。

### 3.1 HMMER 用データ（任意）

手元に Pfam HMM（PF00486, PF03704, PF00931, PF13424）がある場合は `hmmsearch` を使ってもよいです。

ない場合は、事前に別ステップで取得したものを `SARP_SCRIPT_DIR` に配置し、それを参照してください。

---

## 4. M145 全 ORF から SARPファミリー候補を抽出

`SARP_SCRIPT_DIR` に、SARP検出用のシェルまたは Python/R スクリプト（例: `detect_SARP_domains_M145.sh`）を作成し、以下を行います。

### 4.1 ドメイン検出

`REF_PROT`（protein.faa）に対して:

- PF00486, PF03704, PF00931, PF13424 の HMM を `hmmsearch` 等で走らせる。
- 各ヒットについて:
  - タンパク質ID（RefSeq の protein ID）
  - ドメインID（PF番号）
  - 開始/終了位置、E-value
  を抽出し、`tables/SARP_domain_hits_raw.tsv` として保存。

- DBD + BTAD の組み合わせ検出:
  - 同一タンパク質に PF00486 と PF03704 両方が存在するものを SARP 候補とする。
  - NB-ARC (PF00931) と TPR_12 (PF13424) の有無を記録しておく（サブタイプ分類用）。

- 結果を `tables/SARP_candidates_domain_based.tsv` にまとめる。

列例:

```text
protein_id	gene_id	old_locus_tag	has_DBD	has_BTAD	has_NB_ARC	has_TPR
WP_003xxxxxx	SC_RS27585	SCO5085	TRUE	TRUE	FALSE	FALSE
...
```

`gene_id` / `old_locus_tag` は、既に持っている GFF/マスターテーブルとの対応を使って補完してください。

---

## 5. SARPサブタイプ分類と gene_master への統合

### 5.1 SARP サブタイプ定義

Yan & Xia 2024 などのレビューに従い、サイズとドメイン構成に基づいて SARP を4つに分類します:

- **small SARP**:
  - 長さ: おおよそ 250–350 aa
  - ドメイン: DBD + BTAD のみ（NB-ARC, TPR なし）
  - 例: RedD / SC_RS31630 / SCO5877 (red), SCO5085 / SC_RS27585 (act)

- **medium SARP**:
  - 長さ: おおよそ 550–650 aa
  - ドメイン: DBD + BTAD + NB-ARC
  - 例: SCO3217 / SC_RS18200 (cda), SCO6280 / SC_RS33650 (cpk), SCO6288 / SC_RS33690 (cpk)

- **large SARP**:
  - 長さ: 900–1100 aa 前後
  - ドメイン: DBD + BTAD + NB-ARC + TPR
  - 例: AfsR / SC_RS24295 / SCO4426

- **SARP-LAL**:
  - 長さ: 大型
  - ドメイン: N末端に SARPドメイン（DBD + BTAD）、C末端側に LAL様領域（ATP/GTP結合ドメイン）
  - 他種の SanG, PolR, PimR などと同様のアーキテクチャ。

実装上は:
- `protein.faa` から各タンパク質の長さを取得。
- ドメイン組み合わせ + 長さでサブタイプを決める（閾値はざっくりで良いが、後で手動で気になるものをチェック可能に）。
- 結果を `tables/SARP_list_M145.tsv` として出力:

```text
gene_id	old_locus_tag	gene_name	product	protein_id	length	has_DBD	has_BTAD	has_NB_ARC	has_TPR	SARP_subtype
SC_RS27585	SCO5085	NA	AfsR/SARP family transcriptional regulator	WP_003...	~300	TRUE	TRUE	FALSE	FALSE	small
...
```

> **ID 体系の注意**: マスターテーブルの遺伝子 ID 列は `gene_id`（SC_RSxxxxx 形式）、旧 SCO 表記の列は `old_locus_tag` です。`locus_tag` ではないので注意してください。

> **アノテーション上の注意点**:
> - actII-orf4（SC_RS27570 / SCO5082）は GFF 上で「TetR family transcriptional regulator」と注釈されており、product フィールドには "SARP" の文字が**含まれません**。文献上は SARP として扱われますが、ドメイン検索でDBD + BTADが検出されるかどうかで判定してください。
> - cpkO/kasO（SC_RS33660 / SCO6282）は GFF 上で「SDR family oxidoreductase」と注釈されており、SARP とは記載されていません。文献上は SARP/regulator として議論されることがありますが、ドメイン検出の結果に従ってください。
> - `cdaR` および `cpkN` という gene_name は GFF アノテーションに**存在しません**。cda クラスター内の SARP は SCO3217（SC_RS18200, AfsR/SARP family）、cpk クラスター内の SARP は SCO6280（SC_RS33650）および SCO6288（SC_RS33690）を参照してください。

### 5.2 gene_master への統合

`GENE_MASTER_FILE`（`gene_master_with_BGC_regulators.tsv`）に対して:

- `is_SARP`（TRUE/FALSE）
- `SARP_subtype`（small/medium/large/SARP-LAL/NA）

の2列を追加し、`tables/gene_master_with_SARP.tsv` として保存。

また、`TF_REG_MASTER_FILE`（`regulator_master_for_prioritization.tsv`）が存在する場合は、同様に `is_SARP` / `SARP_subtype` を追加した `tables/regulator_master_with_SARP.tsv` も作成してください。

---

## 6. SARP × BGC / ダイナミクスのサマリ

SARPにフォーカスしたビューを作成します。

### 6.1 SARP-BGC 関連付け

`SARP_list_M145.tsv` と `BGC_definition_manual.tsv` を用いて:

各 SARP 遺伝子が:
- どの BGC クラスター内に座っているか（cluster-situated SARP / CSR）
- BGC の外にあるか（global/pleiotropic SARP 候補）

を判定。

結果を `tables/SARP_BGC_mapping.tsv` として保存:

```text
gene_id	old_locus_tag	gene_name	product	SARP_subtype	bgc_name	role_in_BGC
SC_RS27585	SCO5085	NA	AfsR/SARP family TR	small	act	CSR
SC_RS31630	SCO5877	NA	AfsR/SARP family TR	small	red	CSR
SC_RS18200	SCO3217	NA	AfsR/SARP family TR	medium	cda	CSR
SC_RS33650	SCO6280	NA	AfsR/SARP family TR	medium	cpk	CSR
SC_RS33690	SCO6288	NA	AfsR/SARP family TR	medium	cpk	CSR
SC_RS24295	SCO4426	afsR	transcriptional regulator AfsR	large	none	global
...
```

> **参考**: GFF アノテーションから product フィールドに "AfsR/SARP family" を含む遺伝子は以下の 9 個が確認されています:
> SC_RS06435 (SCO0898), SC_RS18200 (SCO3217, cda), SC_RS24295 (SCO4426, afsR),
> SC_RS25495 (SCO4663), SC_RS27585 (SCO5085, act), SC_RS29350 (SCO5433, tcrA),
> SC_RS31630 (SCO5877, red), SC_RS33650 (SCO6280, cpk), SC_RS33690 (SCO6288, cpk)。
> HMMER ドメイン検索によりこれ以外にも SARP 候補が見つかる可能性があります（特に actII-orf4 など product フィールドが "TetR" と注釈されているもの）。

### 6.2 SARP 発現・BGC 相関

`regulator_BGC_correlation.tsv` と `gene_master_with_SARP.tsv` を使って:

全 SARP について:
- `log2FC_3_vs_1`, `padj_3_vs_1`
- `corr_act`, `corr_red`, `corr_cda`, `corr_cpk`
- `delta_2_vs_1`, `delta_3_vs_2`（存在すれば）

をまとめた `tables/SARP_BGC_activity_summary.tsv` を作成。

ここで特に注目するのは:

- **CSR**: act クラスター SARP（SC_RS27585 / SCO5085）、redD（SC_RS31630 / SCO5877）、cda SARP（SC_RS18200 / SCO3217）、cpk SARP（SC_RS33650 / SCO6280, SC_RS33690 / SCO6288）の挙動が、それぞれの BGCスコアとどう整合しているか。
- **actII-orf4**（SC_RS27570 / SCO5082）: ドメイン検索で SARP と判定された場合は、act CSR としてここに含める。
- **Global SARP**: AfsR（SC_RS24295 / SCO4426）の発現と BGCスコアの関係（多くの BGC と同調するか、特定 BGCに偏っているか）。

### 6.3 ヒートマップ（任意）

- 行: SARP 遺伝子
- 列: `corr_act`, `corr_red`, `corr_cda`, `corr_cpk` または `log2FC_3_vs_1`, `delta_2_vs_1`, `delta_3_vs_2`
- 出力:
  - `${SARP_RUN_DIR}/figures/SARP_BGC_correlation_heatmap_M145.pdf`
  - `${SARP_RUN_DIR}/figures/SARP_expression_heatmap_M145.pdf`

---

## 7. SARP候補の位置付けと 08 結果との統合

`TF_candidates_experimental_priority.tsv`（08での9候補）と SARP 情報を突き合わせます。

- 9候補の中に含まれる SARP があるかをチェック（例: actII-orf4 がドメイン検索で SARP と判定された場合など）。
- `SARP_BGC_activity_summary.tsv` から、08で拾えていないが面白い SARP（例: BGC外にありつつ BGCスコアと高相関な SARP）を抽出。
- 実験優先候補に SARP を追加する場合は、`tables/TF_candidates_experimental_priority_with_SARP.tsv` のような形で、SARP拡張版の優先度リストを作成してもよいです。

---

## 8. SARP レポート `SARP_report_M145.md`

このステップの最後に、必ず `${SARP_RUN_DIR}/SARP_report_M145.md` を作成し、以下を含めてください。

### 8.1 目的とアプローチ

「M145ゲノムから全SARPをドメイン構成に基づいて同定・分類し、既存の BGC / TF ネットワークに統合した」こと。

### 8.2 SARP の数とサブタイプの分布

- small / medium / large / SARP-LAL がそれぞれ何個あったか。
- そのうち何個が BGC 内 CSR、何個が BGC外 global/pleiotropic か。

### 8.3 主要 BGC に関連する SARP の挙動

- **act**: actII-orf4（SC_RS27570, ドメイン検索で判定）と SCO5085（SC_RS27585, act クラスター SARP）の発現・相関パターン。
- **red**: redD（SC_RS31630 / SCO5877, small SARP）
- **cda**: SCO3217（SC_RS18200, medium SARP）
- **cpk**: SCO6280（SC_RS33650）＋ SCO6288（SC_RS33690）
- これらが 06/07/08 の結果とどう整合するか（中期 vs 後期、BGCスコアとの相関）。

### 8.4 global SARP（AfsR など）の位置付け

- AfsR（SC_RS24295 / SCO4426）のサブタイプ（large SARP, NB-ARC/TPRを持つ）と、Act/Red/CDA への既知の正の効果。
- 今回の M145 データで AfsR の発現・BGC相関がどう見えるか（強い相関か、弱いが存在するか）。

### 8.5 新規に注目される SARP 候補

- BGC 外にありつつ、複数 BGC のスコアと高相関な SARP、あるいは特定 BGC と高相関な SARP。
- それらが small/medium/large/SARP-LAL のどれか、またどの BGC クラス（PKS/NRPS など）に紐づくか。

### 8.6 実験的示唆

act/red/cda/cpk の制御設計において、

- CSR SARP（act クラスター SARP, redD, cda SARP, cpk SARP）
- global SARP（AfsR）

をどう組み合わせて KO/OE するべきか、1〜2段階のモデルとして簡潔に書く。

silent / cryptic BGC の活性化に使えそうな SARP（SARP-LAL型や BGC外 small SARP）の候補も触れておく。

### 8.7 重要ポイント（3〜5行）

---

## 9. 成果物チェックリスト

09_SARP_integration ステップ完了の目安:

- [ ] `${SARP_RUN_DIR}/logs/SARP_integration_pipeline.log`
- [ ] `${SARP_RUN_DIR}/tables/SARP_domain_hits_raw.tsv`
- [ ] `${SARP_RUN_DIR}/tables/SARP_candidates_domain_based.tsv`
- [ ] `${SARP_RUN_DIR}/tables/SARP_list_M145.tsv`
- [ ] `${SARP_RUN_DIR}/tables/gene_master_with_SARP.tsv`
- [ ] `${SARP_RUN_DIR}/tables/regulator_master_with_SARP.tsv`
- [ ] `${SARP_RUN_DIR}/tables/SARP_BGC_mapping.tsv`
- [ ] `${SARP_RUN_DIR}/tables/SARP_BGC_activity_summary.tsv`
- [ ] （任意）`${SARP_RUN_DIR}/figures/SARP_BGC_correlation_heatmap_M145.pdf`
- [ ] `${SARP_RUN_DIR}/SARP_report_M145.md`
