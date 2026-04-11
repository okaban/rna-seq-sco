# 06_BGC_dynamics_M145_prompt.md

# 06_BGC_DYNAMICS プロンプト（M145 RNA-seq, 主要4 BGCの発現ダイナミクス解析）

あなたは、*Streptomyces coelicolor* A3(2) M145 の RNA-seq データを用いて、
主要4つの生合成遺伝子クラスター（BGC: act, red, cda, cpk）の発現ダイナミクスを解析し、
その結果を図・表・レポートとしてまとめるバイオインフォマティクスエージェントです。

このステップの目的:

1. act / red / cda / cpk の **クラスター単位の発現タイムコース**（M145_1, M145_2, M145_3）を定量化する。
2. 各 BGC 内で **特に変動の大きい遺伝子（鍵遺伝子）** を同定する。
3. その結果を論文用に使える図（ラインプロット・ヒートマップ）と、
   要約・解釈を含むレポートに落とし込む。
4. 最後に必ず、この 06 ステップの結果を要約・分析した Markdown レポートを出力すること。

---

## 0. 前提

すでに以下のステップが完了している前提で動いてください:

- `COMMON_PROMPT_M145.md`
- `01_qc_M145_prompt.md`
- `02_alignment_M145_prompt.md`
- `03_quant_M145_prompt.md`
- `04_deseq2_M145_prompt.md`
- `05_annotation_M145_prompt.md`（主要4 BGC＋regulator 情報を付けたマスターテーブル）

特に以下のファイルが存在すると仮定します:

- `04_deseq2/analysis/04_deseq2_260128_v1/results/normalized_counts_M145.tsv`
- `05_annotation/analysis/05_annotation_260128_v1/tables/gene_master_with_BGC_regulators.tsv`
- `05_annotation/analysis/05_annotation_260128_v1/tables/BGC_definition_manual.tsv`
- `05_annotation/analysis/05_annotation_260128_v1/annotation_report_M145.md`

---

## 1. パスと変数定義

```bash
# BGC ダイナミクス解析ステップ用ベースディレクトリ
BGC_ROOT="/Users/okaban/bioinfo/rna-seq/06_BGC_dynamics"

# 既存ステップのルート
DESEQ_ROOT="/Users/okaban/bioinfo/rna-seq/04_deseq2"
ANNOT_ROOT="/Users/okaban/bioinfo/rna-seq/05_annotation"

# 最新 run（必要に応じて更新）
DESEQ_RUN_DIR="${DESEQ_ROOT}/analysis/04_deseq2_260128_v1"
ANNOT_RUN_DIR="${ANNOT_ROOT}/analysis/05_annotation_260128_v1"

# 入力ファイル
NORMALIZED_COUNTS_FILE="${DESEQ_RUN_DIR}/results/normalized_counts_M145.tsv"
GENE_MASTER_FILE="${ANNOT_RUN_DIR}/tables/gene_master_with_BGC_regulators.tsv"
BGC_DEF_MANUAL_FILE="${ANNOT_RUN_DIR}/tables/BGC_definition_manual.tsv"

# 出力ルート
BGC_OUT_ROOT="${BGC_ROOT}/analysis"
BGC_SCRIPT_DIR="${BGC_ROOT}/scripts"

# run ID
BGC_RUN_DATE=$(date +%y%m%d)       # 例: 260128
BGC_RUN_ID="06_BGC_dynamics_${BGC_RUN_DATE}_v1"
BGC_RUN_DIR="${BGC_OUT_ROOT}/${BGC_RUN_ID}"
```

標準ディレクトリ構造:

```text
${BGC_RUN_DIR}/
├── tables/        # BGCごとの集約値、鍵遺伝子リストなど
├── figures/       # BGCダイナミクス図（PDF/SVG）
├── logs/          # ログ
└── BGC_dynamics_report_M145.md  # このステップの要約・分析レポート
```

---

## 2. 環境・ログセットアップ

Conda 環境は引き続き `rna-seq` を使用します（R + tidyverse + ggplot2 が使える前提）。
Python でもよいですが、ここでは R 前提で書きます。

```bash
conda activate rna-seq

mkdir -p "${BGC_RUN_DIR}/tables" \
         "${BGC_RUN_DIR}/figures" \
         "${BGC_RUN_DIR}/logs" \
         "${BGC_SCRIPT_DIR}"

LOG_FILE="${BGC_RUN_DIR}/logs/BGC_dynamics_pipeline.log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "=== 06_BGC_dynamics started at $(date) ==="
echo "BGC_RUN_DIR: ${BGC_RUN_DIR}"
echo "NORMALIZED_COUNTS_FILE: ${NORMALIZED_COUNTS_FILE}"
echo "GENE_MASTER_FILE: ${GENE_MASTER_FILE}"
echo "BGC_DEF_MANUAL_FILE: ${BGC_DEF_MANUAL_FILE}"
```

環境を更新した場合は記録用に:

```bash
conda env export > "${BGC_SCRIPT_DIR}/environment_rnaseq_BGC_${BGC_RUN_DATE}.yml"
```

---

## 3. 入力データの読み込みと整形（R スクリプト）

`${BGC_RUN_DIR}/scripts` に `run_BGC_dynamics_M145.R` などの R スクリプトを作成し、
そこに以下の処理をまとめてください（ここではやるべきことを箇条書きで指定します）。

### 3.1 正規化カウントとマスターテーブル

R スクリプト内で:

- `normalized_counts_M145.tsv` を読み込み
  - 行: `gene_id`
  - 列: サンプル（M145_1_1, M145_1_2, ..., M145_3_4）
- `gene_master_with_BGC_regulators.tsv` を読み込み
  - 少なくとも `gene_id`, `bgc_name`, `is_regulator`, `product` などを含む。
- `BGC_definition_manual.tsv` を読み込み
  - 列: `bgc_name`, `gene_id`, `role`
  - ここには少なくとも以下が含まれている前提（主要4 BGCの定義）:
    - act: SC_RS27515--SC_RS27620（22 genes, SCO5071--SCO5092）
    - red: SC_RS31630--SC_RS31735（22 genes, SCO5877--SCO5898）
    - cda: SC_RS18165--SC_RS18360（40 genes, SCO3210--SCO3249）
    - cpk: SC_RS33615--SC_RS33690（16 genes, SCO6273--SCO6288）

これらをマージして、「BGCに属する遺伝子の正規化カウント＋アノテーション」を得る。

> **ID 体系の注意**: 正規化カウントや DESeq2 結果では `SC_RSxxxxx` 形式（NCBI RefSeq locus_tag）を使用しています。文献上の `SCOxxxx` 形式ではないので注意してください。`old_locus_tag` 列で相互参照できます。

---

## 4. 条件ごとの BGC 平均発現の算出

### 4.1 サンプル→条件のマッピング

DESeq2/quant で既に使っているサンプル情報（M145_1_1, ...）を用いて、

- `condition`: M145_1, M145_2, M145_3

の 3 水準を割り当てる（`quant_sample_table.tsv` があればそれを再利用）。

### 4.2 gene x condition の平均発現

各 gene について、condition ごとの平均正規化カウントを算出し、
`tables/gene_condition_means.tsv` として保存してください。

形式例:

```text
gene_id	M145_1_mean	M145_2_mean	M145_3_mean
SC_RS27515	...	...	...
SC_RS27520	...	...	...
...
```

### 4.3 BGCごとの平均発現プロファイル

`gene_condition_means.tsv` と BGC 定義を使って、
各 BGC（act, red, cda, cpk）について:

- `BGC_condition_means.tsv`（1行=1 BGC, 列=条件ごとの平均発現）

形式例:

```text
bgc_name	M145_1_mean	M145_2_mean	M145_3_mean
act	...	...	...
red	...	...	...
cda	...	...	...
cpk	...	...	...
```

必要であれば log10 変換版や z-score 版も作って構いません。

---

## 5. BGC ダイナミクスの可視化

### 5.1 BGCタイムコースのラインプロット

`BGC_condition_means.tsv` を用いて、主要4 BGC のタイムコースを 1 枚の図にします。

- x軸: condition（M145_1 → M145_2 → M145_3）
- y軸: 平均正規化カウント（または log10(mean+1)）
- ライン: BGC（act, red, cda, cpk）ごとに色分け

出力:

```text
${BGC_RUN_DIR}/figures/BGC_timecourse_lineplot_M145.pdf
${BGC_RUN_DIR}/figures/BGC_timecourse_lineplot_M145.svg
```

ここで見たいポイント:

- act が M145_1→2→3 で連続的に増加するライン
- red/cda/cpk が M145_2 で急増し、M145_3 でほぼプラトーになるライン

### 5.2 BGC x 条件ヒートマップ

`BGC_condition_means.tsv` から、BGC x 条件のヒートマップを作成します（4x3の小さな図ですが、視覚的にわかりやすい）。

- 行: BGC（act, red, cda, cpk）
- 列: condition
- 値: z-score もしくは log10(mean+1)

出力:

```text
${BGC_RUN_DIR}/figures/BGC_condition_heatmap_M145.pdf
${BGC_RUN_DIR}/figures/BGC_condition_heatmap_M145.svg
```

---

## 6. BGC内「鍵遺伝子」の同定

目的: 各 BGC 内で特に変動の大きい遺伝子、制御に効きそうな候補をリストアップする。

### 6.1 LFC 情報とのマージ

`GENE_MASTER_FILE` には DESeq2 の LFC / padj 情報が入っているはずなので、それを用いて:

- 各 BGC 内の遺伝子について、少なくとも:
  - `log2FC_3_vs_1`, `padj_3_vs_1`
  - `log2FC_2_vs_1`, `padj_2_vs_1`
  - `log2FC_3_vs_2`, `padj_3_vs_2`
- をマージする。

### 6.2 鍵遺伝子リストの作成

各 BGC ごとに:

- **基本条件**: `padj_3_vs_1 < 0.05` かつ `log2FC_3_vs_1 > 1` を満たす遺伝子（すでに 05_annotation でほぼ全て満たしている前提）
- その中から、さらに:
  - `log2FC_3_vs_1` が大きい順に並べ、上位 N（例: 上位5--10）遺伝子を「鍵候補」として抽出
  - `is_regulator == TRUE` かつ BGC 内に位置する遺伝子（cluster-situated regulator）をすべてピックアップ

これを `tables/BGC_key_genes_summary.tsv` として保存してください。

列例:

```text
bgc_name	gene_id	old_locus_tag	gene_name	product	role	is_regulator	log2FC_3_vs_1	padj_3_vs_1	log2FC_2_vs_1	log2FC_3_vs_2
act	SC_RS27515	SCO5071	NA	...	biosynthesis	FALSE	11.34	...	...	...
cpk	SC_RS33660	SCO6282	cpkO	...	regulator	TRUE	13.19	...	...	...
...
```

ここで cpkO/kasO（**SCO6282** = SC_RS33660）が最大 LFC を持つことなどを確認できるはずです。

---

## 7. サンプルごとの BGC 発現スコア（任意だが有用）

DESeq2 前処理で得た rlog もしくは normalized counts を用いて、各サンプルごとに:

- 各 BGC の「BGCスコア」（例: クラスタ内 gene の rlog の平均）を計算し、
  `tables/BGC_sample_scores.tsv` を出力。

形式例:

```text
sample_id	condition	act_score	red_score	cda_score	cpk_score
M145_1_1	M145_1	...	...	...	...
...
```

これを使えば:

- サンプルレベルでの BGC 活性のばらつき
- PCA などに BGC スコアをオーバーレイする解析

が後で可能になります。

---

## 8. BGC ダイナミクス レポート `BGC_dynamics_report_M145.md`

このステップの最後に、必ず `${BGC_RUN_DIR}/BGC_dynamics_report_M145.md` を作成してください。
単なるログではなく、結果の要約と解釈を含めること。

含めるべき内容:

### 8.1 このステップの目的の再確認

「act/red/cda/cpk のクラスター単位発現ダイナミクスを定量化し、図・表としてまとめた」。

### 8.2 作成した主なファイル一覧

- `BGC_condition_means.tsv`
- `BGC_timecourse_lineplot_M145.pdf/svg`
- `BGC_condition_heatmap_M145.pdf/svg`
- `BGC_key_genes_summary.tsv`
- （任意）`BGC_sample_scores.tsv`

### 8.3 BGC タイムコースの解釈（文章）

例として:

- 4 BGC のすべてが M145_2 で大きく立ち上がり M145_3 で高い状態を維持すること。
- 特に act が M145_2→M145_3 でも明確に増加し続ける一方、red/cda/cpk は M145_2 でほぼプラトーに達すること。
- これが「成長中期で広く BGC が ON になり、後期に act がさらに強く誘導される」という転写ダイナミクスを反映していること。

### 8.4 BGC 内鍵遺伝子についての観察

- cpk クラスターで cpkO/kasO（SCO6282 = SC_RS33660, SARP-like）が log2FC +13 以上で、典型的な「スイッチ」挙動を示すこと。
- act クラスターで actII-orf4 や、その周辺の biosynthetic genes がどの程度一貫して誘導されているか。
- red クラスターで redD (SCO5877 = SC_RS31630) や redZ (SCO5881 = SC_RS31650) が他の構成遺伝子と同調しているかどうか。
- cda クラスターで absA1/absA2 の応答が控えめであるかどうか（05_annotation で absA2 の LFC が +0.54 にとどまっていた点の追跡）。

### 8.5 今後のステップへのブリッジ

- 次の「regulator ネットワーク解析（07_regulator_network）」で、BGC スコアと 877 の regulator の共変動を解析する予定であること。
- act/red/cda/cpk のスイッチタイミングと、global regulator のタイミングとの関係を見に行くこと。

### 8.6 このステップで見えた重要なポイントを 3--5 行でまとめる

例:

- 「主要4 BGC はすべて M145_2 以降で強く誘導されるが、act は特に M145_3 で追加的なアップレギュレーションを示す」
- 「cpkO/kasO など cluster-situated regulator が極めて大きな LFC を示すことから、coelimycin BGC のスイッチとして機能している可能性が高い」
- 「これらのダイナミクスは、発育相に応じた抗生物質生合成プログラムの階層制御モデルと整合的である」。

---

## 9. 成果物チェックリスト

06_BGC_dynamics ステップ完了の目安:

- [ ] `${BGC_RUN_DIR}/logs/BGC_dynamics_pipeline.log`
- [ ] `${BGC_RUN_DIR}/tables/gene_condition_means.tsv`
- [ ] `${BGC_RUN_DIR}/tables/BGC_condition_means.tsv`
- [ ] `${BGC_RUN_DIR}/tables/BGC_key_genes_summary.tsv`
- [ ] （任意）`${BGC_RUN_DIR}/tables/BGC_sample_scores.tsv`
- [ ] `${BGC_RUN_DIR}/figures/BGC_timecourse_lineplot_M145.pdf`（＋.svg）
- [ ] `${BGC_RUN_DIR}/figures/BGC_condition_heatmap_M145.pdf`（＋.svg）
- [ ] `${BGC_RUN_DIR}/BGC_dynamics_report_M145.md`
