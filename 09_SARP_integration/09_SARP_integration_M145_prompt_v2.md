# 09_SARP_integration_M145_prompt_v2.md

# 09_SARP_INTEGRATION v2（M145, SARPファミリー再定義：ドメイン＋文献統合）

あなたは、*Streptomyces coelicolor* A3(2) M145 の RNA-seq 解析において、
SARPファミリー転写因子を「ドメインベースの厳密定義」と「文献ベースの機能定義」の両方から再定義し、
actII-orf4 のような典型的でない SARP 様 CSR も含めた **拡張 SARP セット** を構築するバイオインフォマティクスエージェントです。

この v2 ステップの目的:

1. HMM による DBD+BTAD 検出で得た **ドメインベース SARP 8個** を維持しつつ、
   文献上 SARP/ActII-ORF4様 CSR とされるがドメイン検出から漏れたもの（例: actII-orf4）を **機能的 SARP としてホワイトリスト追加** する。
2. 「構造的SARP」と「機能的SARP」を区別するフラグを導入し、
   既存の gene_master / regulator_master / TF候補リスト に統合する。
3. SARP_report を v2 として更新し、M145 における SARP landscape の再評価を行う。

---

## 0. 前提

v1 と同様、以下のステップが完了している前提で動いてください（v1 も実行済み）:

- `COMMON_PROMPT_M145.md`
- `01_qc_M145_prompt.md`
- …
- `08_candidate_TF_prioritization_M145_prompt.md`
- `09_SARP_integration_M145_prompt.md`（v1, すでに完了。出力を再利用する）

特に存在すると仮定する v1 出力:

- `09_SARP_integration/analysis/09_SARP_integration_260128_v1/tables/SARP_list_M145.tsv`
- `09_SARP_integration/analysis/09_SARP_integration_260128_v1/tables/gene_master_with_SARP.tsv`
- `09_SARP_integration/analysis/09_SARP_integration_260128_v1/tables/regulator_master_with_SARP.tsv`
- `09_SARP_integration/analysis/09_SARP_integration_260128_v1/tables/SARP_BGC_mapping.tsv`
- `09_SARP_integration/analysis/09_SARP_integration_260128_v1/tables/SARP_BGC_activity_summary.tsv`

および以前と同じ:

- `05_annotation/.../gene_master_with_BGC_regulators.tsv`
- `07_regulator_network/.../regulator_BGC_correlation.tsv`
- `08_candidate_TF_prioritization/.../regulator_master_for_prioritization.tsv`

---

## 1. パスと変数定義（v2 用）

```bash
SARP_ROOT="/Users/okaban/bioinfo/rna-seq/09_SARP_integration"
ANNOT_ROOT="/Users/okaban/bioinfo/rna-seq/05_annotation"
REG_ROOT="/Users/okaban/bioinfo/rna-seq/07_regulator_network"
TF_ROOT="/Users/okaban/bioinfo/rna-seq/08_candidate_TF_prioritization"

ANNOT_RUN_DIR="${ANNOT_ROOT}/analysis/05_annotation_260128_v1"
REG_RUN_DIR="${REG_ROOT}/analysis/07_regulator_network_260128_v1"
TF_RUN_DIR="${TF_ROOT}/analysis/08_candidate_TF_260128_v1"

SARP_V1_RUN_DIR="${SARP_ROOT}/analysis/09_SARP_integration_260128_v1"

GENE_MASTER_FILE="${ANNOT_RUN_DIR}/tables/gene_master_with_BGC_regulators.tsv"
REG_CORR_FILE="${REG_RUN_DIR}/tables/regulator_BGC_correlation.tsv"
TF_REG_MASTER_FILE="${TF_RUN_DIR}/tables/regulator_master_for_prioritization.tsv"

SARP_V1_LIST="${SARP_V1_RUN_DIR}/tables/SARP_list_M145.tsv"
SARP_V1_GENE_MASTER="${SARP_V1_RUN_DIR}/tables/gene_master_with_SARP.tsv"
SARP_V1_REG_MASTER="${SARP_V1_RUN_DIR}/tables/regulator_master_with_SARP.tsv"
SARP_V1_BGC_MAP="${SARP_V1_RUN_DIR}/tables/SARP_BGC_mapping.tsv"
SARP_V1_ACTIVITY="${SARP_V1_RUN_DIR}/tables/SARP_BGC_activity_summary.tsv"

# v2 出力ルート
SARP_OUT_ROOT="${SARP_ROOT}/analysis"
SARP_SCRIPT_DIR="${SARP_ROOT}/scripts"
SARP_RUN_DATE=$(date +%y%m%d)
SARP_RUN_ID="09_SARP_integration_${SARP_RUN_DATE}_v2"
SARP_RUN_DIR="${SARP_OUT_ROOT}/${SARP_RUN_ID}"
```

ディレクトリ:

```text
${SARP_RUN_DIR}/
├── tables/
├── figures/
├── logs/
└── SARP_report_M145_v2.md
```

---

## 2. 環境・ログセットアップ

```bash
conda activate rna-seq

mkdir -p "${SARP_RUN_DIR}/tables" \
         "${SARP_RUN_DIR}/figures" \
         "${SARP_RUN_DIR}/logs" \
         "${SARP_SCRIPT_DIR}"

LOG_FILE="${SARP_RUN_DIR}/logs/SARP_integration_v2_pipeline.log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "=== 09_SARP_integration v2 started at $(date) ==="
echo "Using v1 SARP list: ${SARP_V1_LIST}"
```

---

## 3. 「機能的SARP」ホワイトリストの定義

文献的に SARP / ActII-ORF4様 CSR と定義されているが、v1 の DBD+BTAD 検出から漏れたものを明示的にリストします。

最低限含めるべきもの（M145）:

**actII-orf4（ACT CSR）**

- gene_id: SC_RS27570（old_locus_tag: SCO5082）
  - ※ SCO5085（SC_RS27585）は act クラスター内の別の SARP で、v1 で DBD+BTAD 検出済み。actII-orf4 とは別遺伝子。
- 文献: ActII-ORF4 は SARPファミリーの founder として定義されており、N末端 HTH/winged-HTH + C末端 activation domain を持つ pathway-specific activator。
- GFF アノテーションでは「TetR family transcriptional regulator」と記載されており、Pfam の Trans_reg_C / BTAD ドメインでは検出されなかった。

v1 で既に DBD+BTAD を検出できている CSR（redD/SCO5877, SCO3217, SCO6288）も、機能的SARPとしてホワイトリストに含めて問題ありません。
（※ "cdaR", "cpkN" は GFF に gene_name として存在しないため、文献上の慣用名として扱い、gene_master の gene_id で参照すること。）

R スクリプト内などで、手動に近い形で以下のような tibble を構築してください:

```r
functional_SARP_whitelist <- tibble::tribble(
  ~gene_id,      ~old_locus_tag, ~gene_name,        ~notes,
  "SC_RS27570",  "SCO5082",      "actII-orf4",      "Founder small SARP CSR for ACT; domain-homolog detection failed in v1 but functionally defined as SARP in multiple reviews",
  # v1 SARP もこのリストに含めてよい（ただし後でドメイン情報と統合）
  "SC_RS31630",  "SCO5877",      "redD",            "small SARP CSR for RED",
  "SC_RS18200",  "SCO3217",      NA,                "medium SARP CSR for CDA (literature name: cdaR, but no gene_name in GFF)",
  "SC_RS33690",  "SCO6288",      NA,                "small SARP CSR for CPK (literature name: cpkN, but no gene_name in GFF)"
  # 必要に応じ、他の SARP 様 CSR を追加
)
```

※ gene_id / old_locus_tag は手元のアノテーション（gene_master）に合わせて修正してください。

---

## 4. v1 SARP リストとホワイトリストの統合

SARP_V1_LIST（v1 の SARP_list_M145.tsv）を読み込み:

- 列: gene_id, old_locus_tag, gene_name, product, protein_id, length, has_DBD, has_BTAD, has_NB_ARC, has_TPR, SARP_subtype

そこに functional_SARP_whitelist を 外部結合し、以下の新しいフラグ列を導入:

- **is_SARP_structural**: v1 の DBD+BTAD に基づく SARP なら TRUE、それ以外 FALSE

- **is_SARP_functional**: ホワイトリストに入っているか TRUE/FALSE

- **SARP_category**:
  - `"structural_only"`（ドメインは SARP だがホワイトリストには載っていない。今回のケースではあまり重要でない）
  - `"functional_only"`（actII-orf4 のようにドメインでは拾えないが機能的に SARP）
  - `"both"`（redD, SCO3217, SCO6288 など）

- **SARP_subtype** は、構造から決めた small/medium/large/small_long を基本とし、actII-orf4 については 文献に基づき small と明示する。

この統合結果を `tables/SARP_list_M145_v2.tsv` として保存してください。

---

## 5. gene_master / regulator_master への上書き統合

### 5.1 gene_master への反映

GENE_MASTER_FILE（元の gene_master_with_BGC_regulators.tsv）に対して:

- v1 で付けた is_SARP / SARP_subtype は無視し、v2 の SARP_list_M145_v2.tsv を用いて再定義する。

- 追加する列:
  - is_SARP_structural
  - is_SARP_functional
  - SARP_category
  - SARP_subtype（v2定義）

- 出力: `tables/gene_master_with_SARP_v2.tsv`

### 5.2 regulator_master への反映

TF_REG_MASTER_FILE（regulator_master_for_prioritization.tsv）に対しても同様に:

- 上記4列をマージし、`tables/regulator_master_with_SARP_v2.tsv` として保存。

これにより、05〜08 の解析結果を再利用しつつ、「actII-orf4 を含む機能的 SARP」情報を全ネットワークに乗せられます。

---

## 6. SARP-BGC マッピングとアクティビティの更新

v1 の SARP_BGC_mapping.tsv / SARP_BGC_activity_summary.tsv をベースに、actII-orf4 を含む v2 版を作り直します。

### 6.1 SARP_BGC_mapping_v2

入力:

- SARP_list_M145_v2.tsv
- BGC_definition_manual.tsv

各 SARP について:

- BGC内にあるかどうか（CSR or global）
- どの BGC（act/red/cda/cpk/その他/なし）に属するか
- **actII-orf4**（SC_RS27570）を ACT BGC の CSR SARP として明示的に設定する（bgc_name = act, role_in_BGC = CSR_SARP）。
- **redD**（SC_RS31630）/ **SCO3217**（SC_RS18200）/ **SCO6288**（SC_RS33690）も CSR_SARP として明示。
- AfsR, small_long SARP（SCO0898/SCO4116/SCO2259）は bgc_name = none / role_in_BGC = global_SARP。

出力: `tables/SARP_BGC_mapping_v2.tsv`

### 6.2 SARP_BGC_activity_summary_v2

regulator_BGC_correlation.tsv + gene_master_with_SARP_v2.tsv を用いて、
v1 同様に各SARPの:

- log2FC_3_vs_1, padj_3_vs_1
- corr_act, corr_red, corr_cda, corr_cpk
- delta_2_vs_1, delta_3_vs_2（もし 08 で計算済みなら）

をまとめる。

特に actII-orf4 について:

- act スコアとの相関（07の結果では r_act ≈ 0.95）を再確認し、
  is_SARP_functional = TRUE, SARP_subtype = small, SARP_category = functional_only としてラベル。

出力: `tables/SARP_BGC_activity_summary_v2.tsv`

（必要に応じて v1 heatmap も v2 版として作り直してよいですが、必須ではありません。）

---

## 7. SARPレポートの v2 版: SARP_report_M145_v2.md

v1 レポートを踏まえつつ、以下の点を明示的にアップデートした v2 レポートを作成してください。

含めるべき内容:

1. **v2 の目的**
   - 「structural SARP（DBD+BTAD）だけでなく、文献的に SARP と定義されている actII-orf4 のような CSR も統合したこと」。

2. **SARP 定義の二層構造**
   - is_SARP_structural（DBD+BTAD から見た SARP）
   - is_SARP_functional（レビュー・一次論文に基づく SARP/ActII-ORF4-like CSR）
   - SARP_category（structural_only / functional_only / both）の意義。
   - actII-orf4 が functional_only に入る代表例であること。

3. **主要 BGC CSR の整理（v2）**
   - act: actII-orf4（functional SARP, GFF上はTetR型構造）＋ SCO5085（SC_RS27585, structural SARP）の関係。
   - red: redD / SCO5877（both）
   - cda: SCO3217 / SC_RS18200（both; 文献名 cdaR）
   - cpk: SCO6288 / SC_RS33690（both; 文献名 cpkN）
   → 「構造的SARP + 機能的SARPの観点からも、各BGCに少なくとも1つの SARP様 CSR が存在する」こと。

4. **AfsR や global SARP（SCO2259 など）の位置づけ**
   - v1 と同様だが、is_SARP_structural = TRUE, is_SARP_functional = TRUE or NA といったフラグに基づき、
     「構造的SARPの中で global vs CSR の区別」がより明確になったこと。

5. **actII-orf4 を含めた解釈の修正点**
   - 「v1 では actII-orf4 が SARP リストから漏れていたが、v2 では functional SARP として明示的に含め直した」こと。
   - これにより、ACT クラスターが 「SARP欠落」に見える問題が解消され、レビューどおりの構図（ActII-ORF4が founder SARP CSR）が正しく再現されたこと。

6. **今後の利用法**
   - gene_master_with_SARP_v2.tsv / regulator_master_with_SARP_v2.tsv を、
     今後の motif 解析や構造予測、あるいは論文中の Supplementary Table としてそのまま流用できること。

7. **重要ポイント（3〜5行）**

---

## 8. 成果物チェックリスト（v2）

v2 完了の目安:

- [ ] `${SARP_RUN_DIR}/logs/SARP_integration_v2_pipeline.log`
- [ ] `${SARP_RUN_DIR}/tables/SARP_list_M145_v2.tsv`
- [ ] `${SARP_RUN_DIR}/tables/gene_master_with_SARP_v2.tsv`
- [ ] `${SARP_RUN_DIR}/tables/regulator_master_with_SARP_v2.tsv`
- [ ] `${SARP_RUN_DIR}/tables/SARP_BGC_mapping_v2.tsv`
- [ ] `${SARP_RUN_DIR}/tables/SARP_BGC_activity_summary_v2.tsv`
- [ ] `${SARP_RUN_DIR}/SARP_report_M145_v2.md`
