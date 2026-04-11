# 05_annotation_M145_prompt.md

# 05_ANNOTATION プロンプト（M145 RNA-seq, 機能・BGC・regulator 統合アノテーション）

あなたは、*Streptomyces coelicolor* A3(2) M145 の RNA-seq 解析において、
DESeq2 までの結果を統合的に解釈するための **マスターテーブル（遺伝子アノテーション＋DE 情報）** を構築し、
その要約・分析レポートを出力するバイオインフォマティクスエージェントです。

このステップの目的は、後続で行う:

1. DEGs リストからの機能・経路レベル解釈（GO/KEGG/COG enrichment など）
2. BGC・クラスター単位の解析（主要 4 BGC: act・red・cda・cpk）
3. 転写制御・ネットワークレベルの解析（regulator と BGC の関係）

の **共通土台となるアノテーション付きテーブル** を作ることです。

必ず、解析が終わったタイミングで「05_annotation の結果を要約・分析した Markdown レポート」を出力してください（詳細は最後のセクション）。

---

## 0. 前提

すでに以下が完了している前提で動いてください:

- `COMMON_PROMPT_M145.md`
- `01_qc_M145_prompt.md`
- `02_alignment_M145_prompt.md`（HISAT2）
- `03_quant_M145_prompt.md`（featureCounts）
- `04_deseq2_M145_prompt.md`（DESeq2）

とくに以下のファイルが存在する:

- `03_quantification/analysis/03_quant_260128_v1/counts/featureCounts_M145.txt`
- `04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_2_vs_1.tsv`
- `04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_3_vs_1.tsv`
- `04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_3_vs_2.tsv`
- `04_deseq2/analysis/04_deseq2_260128_v1/results/normalized_counts_M145.tsv`

---

## 1. パスと変数定義

以下をデフォルトとして使用してください（実際の run ID に合わせて適宜書き換える）。

```bash
# アノテーションステップ用のベースディレクトリ
ANNOT_ROOT="/Users/okaban/bioinfo/rna-seq/05_annotation"

# 既存ステップ
QUANT_ROOT="/Users/okaban/bioinfo/rna-seq/03_quantification"
DESEQ_ROOT="/Users/okaban/bioinfo/rna-seq/04_deseq2"

# 最新 run（必要に応じて更新）
QUANT_RUN_DIR="${QUANT_ROOT}/analysis/03_quant_260128_v1"
DESEQ_RUN_DIR="${DESEQ_ROOT}/analysis/04_deseq2_260128_v1"

# 入力ファイル
FEATURECOUNTS_FILE="${QUANT_RUN_DIR}/counts/featureCounts_M145.txt"
NORMALIZED_COUNTS_FILE="${DESEQ_RUN_DIR}/results/normalized_counts_M145.tsv"
DE_2_vs_1_FILE="${DESEQ_RUN_DIR}/results/DESeq2_M145_2_vs_1.tsv"
DE_3_vs_1_FILE="${DESEQ_RUN_DIR}/results/DESeq2_M145_3_vs_1.tsv"
DE_3_vs_2_FILE="${DESEQ_RUN_DIR}/results/DESeq2_M145_3_vs_2.tsv"

# リファレンスアノテーション
REF_DIR="/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1"
REF_GFF="${REF_DIR}/genomic.gff"
REF_GTF="${REF_DIR}/genomic.gtf"
REF_GBFF="${REF_DIR}/genomic.gbff"
REF_PROT="${REF_DIR}/protein.faa"

# 出力ルート
ANNOT_OUT_ROOT="${ANNOT_ROOT}/analysis"
ANNOT_SCRIPT_DIR="${ANNOT_ROOT}/scripts"

# run ID
ANNOT_RUN_DATE=$(date +%y%m%d)        # 例: 260128
ANNOT_RUN_ID="05_annotation_${ANNOT_RUN_DATE}_v1"
ANNOT_RUN_DIR="${ANNOT_OUT_ROOT}/${ANNOT_RUN_ID}"
```

`ANNOT_RUN_DIR` の標準構造:

```text
${ANNOT_RUN_DIR}/
├── tables/        # マスターテーブル、サブテーブル
├── logs/          # ログ
├── scripts/       # R / Python スクリプト
└── annotation_report_M145.md  # このステップの要約・分析レポート
```

---

## 2. 環境セットアップとログ

Conda 環境は引き続き `rna-seq` を使用し、R / Python どちらを使ってもよいが、一貫性のため:

- 軽いパース・整形: Python (pandas) or R (tidyverse)
- Bioconductor 系パッケージが必要なら R 内で

とする。

```bash
conda activate rna-seq

mkdir -p "${ANNOT_RUN_DIR}/tables" \
         "${ANNOT_RUN_DIR}/logs" \
         "${ANNOT_RUN_DIR}/scripts"

LOG_FILE="${ANNOT_RUN_DIR}/logs/annotation_pipeline.log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "=== 05_annotation started at $(date) ==="
echo "ANNOT_RUN_DIR: ${ANNOT_RUN_DIR}"
echo "FEATURECOUNTS_FILE: ${FEATURECOUNTS_FILE}"
echo "NORMALIZED_COUNTS_FILE: ${NORMALIZED_COUNTS_FILE}"
```

必要に応じて環境 export:

```bash
conda env export > "${ANNOT_SCRIPT_DIR}/environment_rnaseq_annotation_${ANNOT_RUN_DATE}.yml"
```

---

## 3. gene ベースの基本アノテーションテーブル作成

### 3.1 GFF からの情報抽出

`REF_GFF` から、少なくとも以下の列を持つ gene-level テーブルを作成し、`tables/gene_annotation_basic.tsv` として保存:

- `gene_id`（DESeq2 の gene_id と対応するもの = GFF の `locus_tag`、例: `SC_RS27510`）
- `old_locus_tag`（旧 Sanger 表記、例: `SCO5070`。GFF の `old_locus_tag` 属性から取得）
- `gene_name`（GFF の `gene=` 属性。例: `actII`, `actVA`, `redW` など。ない場合は NA）
- `product`（機能注釈。**注意: GFF では `product` は CDS レコードにのみ存在し、gene レコードには含まれない。** gene の `locus_tag` と CDS の `locus_tag` をキーにして結合する）
- `start`, `end`, `strand`, `contig`
- `gene_biotype`（`protein_coding`, `rRNA`, `tRNA`, `pseudogene` など）

可能なら:

- `Ontology_term`（GO アノテーション。CDS レコードの `Ontology_term` 属性から取得可能な遺伝子がある）
- `protein_id`（RefSeq protein accession、例: `WP_011026778.1`）

> **ID 体系の注意**: 本プロジェクトの DESeq2 結果では `SC_RSxxxxx`（NCBI RefSeq locus_tag）形式を使用しています。文献では旧 Sanger 表記 `SCOxxxx` が使われることが多いため、`old_locus_tag` 列を必ず含めて両体系を相互参照できるようにしてください。

### 3.2 DESeq2 / counts 情報のマージ

`gene_annotation_basic.tsv` に対して:

- `DESeq2_M145_2_vs_1.tsv`
- `DESeq2_M145_3_vs_1.tsv`
- `DESeq2_M145_3_vs_2.tsv`
- `normalized_counts_M145.tsv`（各サンプルの正規化カウント）

を `gene_id` でマージし、統合マトリクスを `tables/gene_master_DESeq2.tsv` として保存。

推奨する列の構成（例）:

- **基本情報**: `gene_id`, `old_locus_tag`, `gene_name`, `product`, `contig`, `start`, `end`, `strand`
- **DESeq2（各コントラストごと）**: `log2FC_2_vs_1`, `padj_2_vs_1`, `log2FC_3_vs_1`, `padj_3_vs_1`, `log2FC_3_vs_2`, `padj_3_vs_2`, `baseMean` など
- **正規化カウント（列が多いので必要に応じて）**: `norm_M145_1_1`, `norm_M145_1_2`, ..., `norm_M145_3_4`

---

## 4. 主要 4 BGC（act / red / cda / cpk）の情報統合

*S. coelicolor* A3(2) には代表的な 4 つの二次代謝産物 BGC があり、本プロジェクトではこれら **全て** を必須の解析対象とします:

| BGC | 産物 | コア範囲 (SCO) | RefSeq locus_tag 範囲 | 遺伝子数 |
|-----|------|---------------|----------------------|---------|
| **act** | actinorhodin | SCO5071--SCO5092 | SC_RS27515--SC_RS27620 | 22 |
| **red** | undecylprodigiosin | SCO5877--SCO5898 | SC_RS31630--SC_RS31735 | 22 |
| **cda** | calcium-dependent antibiotic | SCO3210--SCO3249 | SC_RS18165--SC_RS18360 | 40 |
| **cpk** | coelimycin P1 (CPK) | SCO6273--SCO6288 | SC_RS33615--SC_RS33690 | 16 |

このステップで、少なくとも以下を明示的に行ってください。

### 4.1 リテラチャーベースの BGC 定義

文献や既存アノテーションに基づき、上記 4 BGC を gene_id レベルで定義します。

- **actinorhodin（act クラスター）**
  - コア範囲: SCO5071--SCO5092 → **SC_RS27515--SC_RS27620**（22 遺伝子）
  - 拡張領域: SCO5070 (SC_RS27510) を含めて 23 遺伝子としてもよい
  - 主要遺伝子: actII-orf4 = SC_RS27570 (SCO5082, SARP regulator), actVA = SC_RS27560 (SCO5080), actVB = SC_RS27620 (SCO5092)

- **undecylprodigiosin（red クラスター）**
  - 範囲: SCO5877--SCO5898 → **SC_RS31630--SC_RS31735**（22 遺伝子）
  - 主要遺伝子: redD = SC_RS31630 (SCO5877, SARP regulator), redW = SC_RS31640 (SCO5879), redM = SC_RS31700 (SCO5891)

- **calcium-dependent antibiotic（cda クラスター）**
  - 範囲: SCO3210--SCO3249 → **SC_RS18165--SC_RS18360**（40 遺伝子）
  - 主要遺伝子: absA1 = SC_RS18240 (SCO3225, two-component sensor kinase), absA2 = SC_RS18245 (SCO3226, response regulator), hppD = SC_RS18260 (SCO3229), asnO = SC_RS18295 (SCO3236)
  - 注: absA1/absA2 は cda クラスター内に位置するが、act・red を含む複数 BGC の発現を制御する広域レギュレーターとしても知られる

- **coelimycin P1（cpk / CPK クラスター）**
  - 範囲: SCO6273--SCO6288 → **SC_RS33615--SC_RS33690**（16 遺伝子）
  - 注: SC_RS33620 (SCO6274) は pseudogene として注釈されている
  - 主要遺伝子: cpkO/kasO = SC_RS33660 (SCO6282, SARP-like regulator), scbR2 = SC_RS33680 (SCO6286, gamma-butyrolactone receptor)
  - 注: scbR (SCO6265, SC_RS33575) および scbA/scbB/scbC (SCO6264--6267) は cpk クラスター境界の直近上流に位置し、gamma-butyrolactone シグナル系を介して cpk を制御する

これら 4 BGC の全遺伝子を `tables/BGC_definition_manual.tsv` にまとめてください。

形式例:

```text
bgc_name	gene_id	old_locus_tag	gene_name	role
act	SC_RS27515	SCO5071	NA	biosynthesis
act	SC_RS27520	SCO5072	NA	biosynthesis
act	SC_RS27570	SCO5082	actII	regulator
...
red	SC_RS31630	SCO5877	redD	regulator
red	SC_RS31640	SCO5879	redW	biosynthesis
...
cda	SC_RS18165	SCO3210	NA	biosynthesis
cda	SC_RS18240	SCO3225	absA1	regulator
cda	SC_RS18245	SCO3226	absA2	regulator
...
cpk	SC_RS33615	SCO6273	NA	biosynthesis
cpk	SC_RS33660	SCO6282	cpkO	regulator
cpk	SC_RS33680	SCO6286	scbR2	regulator
...
```

`role` は `biosynthesis` / `regulator` / `transport` / `resistance` / `unknown` など簡単な分類で構いません。GFF の `product` 注釈や文献情報をもとに判断してください。

### 4.2 gene_master と BGC 情報のマージ

`gene_master_DESeq2.tsv` に、`bgc_name` と `bgc_role` 列を追加してマージし、`tables/gene_master_with_BGC.tsv` を作成:

- BGC 非該当遺伝子は `bgc_name = NA`
- 主要 4 BGC の遺伝子には `bgc_name = act` / `red` / `cda` / `cpk` が入る

これにより、DESeq2 情報と BGC アノテーションが 1 テーブルで見られるようにします。

### 4.3 主要 4 BGC の DE サマリテーブル

主要 4 BGC（act・red・cda・cpk）の遺伝子を抽出し、`tables/BGC_major4_DE_summary.tsv` を作成してください。

列例:

- `bgc_name`
- `gene_id`
- `old_locus_tag`
- `gene_name`
- `product`
- `role`
- `log2FC_3_vs_1`, `padj_3_vs_1`
- `log2FC_2_vs_1`, `padj_2_vs_1`
- `log2FC_3_vs_2`, `padj_3_vs_2`
- （必要なら normalized counts）

これが「各 BGC の生産に寄与しうる候補遺伝子」を見に行くための基礎テーブルになります。

---

## 5. その他 BGC / regulator のタグ付け（あとで使う前処理）

ここでは「タグ付け」まででよく、解析（enrichment やネットワーク）は後続ステップで実行します。

### 5.1 BGC 全般

- 主要 4 BGC（act, red, cda, cpk）は Section 4 で必ず `BGC_definition_manual.tsv` に含まれている。
- それ以外の BGC（hopene, albaflavenone, geosmin, desferrioxamine など）について、antiSMASH 等の結果があれば、その出力から追加の BGC 名と gene_id リストを抽出し、`tables/BGC_definition_all.tsv` として保存。
- `gene_master_with_BGC.tsv` の `bgc_name` をもう少し包括的に埋める。
- antiSMASH 結果がない場合は、主要 4 BGC の手動定義で進め、`BGC_definition_all.tsv` は「可能なら」の扱いとする。

### 5.2 転写因子・regulator タグ

- `product` 注釈（product 列）から、SARP, LuxR, AraC, two-component system, sigma factor など regulator らしき遺伝子に `is_regulator = TRUE` フラグを付ける。
- 特に文献で知られている以下の regulator は手動でも構わないので `regulator_name` / `regulator_type` を追加:

| regulator_name | old_locus_tag | gene_id | regulator_type |
|---|---|---|---|
| actII-orf4 | SCO5082 | SC_RS27570 | SARP (act cluster-situated) |
| redD | SCO5877 | SC_RS31630 | SARP (red cluster-situated) |
| redZ | SCO5881 | SC_RS31650 | response regulator |
| atrA | SCO4677 | SC_RS25560 | TetR family |
| scbR | SCO6265 | SC_RS33575 | gamma-butyrolactone receptor |
| scbR2 | SCO6286 | SC_RS33680 | gamma-butyrolactone receptor |
| cpkO / kasO | SCO6282 | SC_RS33660 | SARP-like (cpk cluster-situated) |
| absA1 | SCO3225 | SC_RS18240 | two-component sensor kinase (cda cluster-situated) |
| absA2 | SCO3226 | SC_RS18245 | response regulator (cda cluster-situated) |

- 結果を `tables/gene_master_with_BGC_regulators.tsv` として保存してください。

---

## 6. 05_annotation の要約・分析レポート `annotation_report_M145.md`

このステップの最後に、必ず `${ANNOT_RUN_DIR}/annotation_report_M145.md` を作成してください。
レポートには、単なるログではなく「結果の要約・分析」を含めてください。

最低限含める内容:

### 6.1 このステップの目的の再確認

「DESeq2 結果とゲノムアノテーション、BGC/regulator 情報を統合したマスターテーブルを作成し、次の enrichment 解析や BGC 解析の土台を整えること」。

### 6.2 作成したテーブルの一覧と簡単な説明

- `gene_annotation_basic.tsv`: gene の基本情報
- `gene_master_DESeq2.tsv`: gene x DESeq2 x normalized counts
- `gene_master_with_BGC.tsv`: 上記に BGC 情報を付与
- `gene_master_with_BGC_regulators.tsv`: さらに regulator タグを付与
- `BGC_definition_manual.tsv`: 主要 4 BGC（act/red/cda/cpk）の手動定義
- `BGC_major4_DE_summary.tsv`: 主要 4 BGC の DE サマリ

### 6.3 主要 4 BGC についての簡単な観察・分析

各クラスターについて:

- act クラスター内の何割の遺伝子が M145_3 vs M145_1 で有意 up か（padj < 0.05 & log2FC > 1）
- red クラスターでも同様に、どの条件で誘導されているか
- cda クラスター遺伝子の発現変動パターン（act/red と同調しているか、異なるか）
- cpk クラスター遺伝子の発現変動パターン
- 各クラスタ内で特に LFC の大きい regulator や transport 系遺伝子があるかどうか

### 6.4 今後のステップへのブリッジ

- 機能・経路 enrichment 解析では、このマスターテーブルを使って GO/KEGG/COG の enrichment を行う予定であること
- BGC 単位解析では、`BGC_major4_DE_summary.tsv` と全 BGC 定義を使って、タイムポイントごとの BGC 発現プロファイルを作る予定であること
- 転写制御解析では、`is_regulator` フラグ付き遺伝子と BGC 発現との関連を見に行く予定であること

### 6.5 「このステップで見えてきた重要なポイント」を 3--5 行程度でまとめる

例:

- act クラスター遺伝子の大半が後期（M145_3）で強く誘導されていること
- red クラスターは act とは異なるタイミングまたは振る舞いを示していること
- cda / cpk クラスターの発現タイミングが act・red と同調しているか、独立したパターンを示しているか
- 一部の regulator（actII-orf4, redD, absA1/absA2, cpkO/kasO, scbR2 など）が BGC 発現と同調して変動しており、今後の候補として注目に値すること

---

## 7. 成果物チェックリスト

05_annotation ステップが完了したとみなす条件:

- [ ] `${ANNOT_RUN_DIR}/logs/annotation_pipeline.log`
- [ ] `${ANNOT_RUN_DIR}/tables/gene_annotation_basic.tsv`
- [ ] `${ANNOT_RUN_DIR}/tables/gene_master_DESeq2.tsv`
- [ ] `${ANNOT_RUN_DIR}/tables/gene_master_with_BGC.tsv`
- [ ] `${ANNOT_RUN_DIR}/tables/gene_master_with_BGC_regulators.tsv`
- [ ] `${ANNOT_RUN_DIR}/tables/BGC_definition_manual.tsv`（主要 4 BGC: act/red/cda/cpk）
- [ ] `${ANNOT_RUN_DIR}/tables/BGC_major4_DE_summary.tsv`
- [ ] （可能なら）`${ANNOT_RUN_DIR}/tables/BGC_definition_all.tsv`
- [ ] `${ANNOT_RUN_DIR}/annotation_report_M145.md`
