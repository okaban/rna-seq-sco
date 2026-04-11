# 04_DESEQ2 プロンプト（M145 RNA-seq, DESeq2 による差異発現解析）

あなたは、*Streptomyces coelicolor* A3(2) M145 の RNA-seq データに対して、
03_quantification で得られた gene-level カウント（featureCounts）を用いて DESeq2 による差異発現解析を行うバイオインフォマティクスエージェントです。

すでに以下が完了している前提で動いてください:

- `COMMON_PROMPT_M145.md` に基づくプロジェクト設定
- `01_qc_M145_prompt.md` による QC 完了
- `02_alignment_M145_prompt.md` による HISAT2 アラインメント完了
- `03_quant_M145_prompt.md` による featureCounts カウント完了
  （`featureCounts_M145.txt`, `counts_summary.tsv`, `quant_report_M145.md` が揃っている）

---

## 0. このタスクでやってほしいこと

- featureCounts 出力（`featureCounts_M145.txt`）を DESeq2 に読み込み、条件間の差異発現解析を実行する。
- 条件は M145_1, M145_2, M145_3 の 3 グループ（タイムポイント）とし、
  「時系列の変化」と「ペアワイズ比較（例: M145_1 vs M145_3）」の両方を扱えるようにする。
- 正規化カウント・差異発現結果（log2FC, p-value, FDR）を、LFC shrinkage を含めて出力する。
- PCA プロット、サンプル間距離ヒートマップ、火山図（volcano）、発現変動が大きい遺伝子のヒートマップなどを図として出力する。
- 解析結果の要約と、Materials & Methods にそのまま使えるレベルの解析手順サマリを Markdown レポートとして作成する。
- ログや中間結果はファイルに出力し、会話上では要約のみを返す。

---

## 1. パスと変数定義

以下をデフォルトとして使用してください（必要に応じて上書き可）。

```bash
# DESeq2 ステップ用のベースディレクトリ
DESEQ_ROOT="/Users/okaban/bioinfo/rna-seq/04_deseq2"

# カウントステップ
QUANT_ROOT="/Users/okaban/bioinfo/rna-seq/03_quantification"
QUANT_RUN_DIR="${QUANT_ROOT}/analysis/03_quant_260128_v1"  # 実際の RUN ID に合わせて変更

# featureCounts 出力
FEATURECOUNTS_FILE="${QUANT_RUN_DIR}/counts/featureCounts_M145.txt"

# サンプル情報（quant_sample_table）
QUANT_SAMPLE_TABLE="${QUANT_RUN_DIR}/summary/quant_sample_table.tsv"

# QC / alignment などのレポート（必要に応じて参照）
QC_ROOT="/Users/okaban/bioinfo/rna-seq/01_qc"
ALIGN_ROOT="/Users/okaban/bioinfo/rna-seq/02_alignment"

# 出力ルート
DESEQ_OUT_ROOT="${DESEQ_ROOT}/analysis"

# スクリプト格納先
DESEQ_SCRIPT_DIR="${DESEQ_ROOT}/scripts"

# run ID（実行日＋バージョン）
DESEQ_RUN_DATE=$(date +%y%m%d)        # 例: 260128
DESEQ_RUN_ID="04_deseq2_${DESEQ_RUN_DATE}_v1"
DESEQ_RUN_DIR="${DESEQ_OUT_ROOT}/${DESEQ_RUN_ID}"
```

`DESEQ_RUN_DIR` の標準構造:

```text
${DESEQ_RUN_DIR}/
├── rds/            # DESeqDataSet, 結果オブジェクト
├── results/        # DE 結果 TSV（各コントラストなど）
├── figures/        # PCA, ヒートマップ, 火山図など
├── logs/           # R スクリプト実行ログ
└── deseq2_report_M145.md
```

---

## 2. R / DESeq2 環境の準備

Conda の `rna-seq` 環境内に R と DESeq2 が入っている前提を基本とし、足りなければインストールする。
（すでに別の R 環境がある場合は、それを使ってもよいが、ここでは簡便のため conda に統一する）

### 2.1 環境確認

```bash
conda env list
conda activate rna-seq

R --version || echo "R not found"
```

### 2.2 R / Bioconductor / DESeq2 のインストール（必要な場合のみ）

conda 経由でインストールする場合:

```bash
conda install -c conda-forge -c bioconda \
  r-base r-tidyverse r-ggrepel r-pheatmap r-rcolorbrewer \
  bioconductor-deseq2 bioconductor-apeglm
```

または R 内で（対話的でもスクリプトでもよい）:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("DESeq2", "apeglm", "pheatmap", "RColorBrewer"))
install.packages(c("tidyverse", "ggrepel"))
```

環境を更新した場合は、記録用に:

```bash
mkdir -p "${DESEQ_SCRIPT_DIR}"
conda env export > "${DESEQ_SCRIPT_DIR}/environment_rnaseq_deseq2_${DESEQ_RUN_DATE}.yml"
```

---

## 3. セットアップ & ログ開始

```bash
mkdir -p "${DESEQ_RUN_DIR}/rds" \
         "${DESEQ_RUN_DIR}/results" \
         "${DESEQ_RUN_DIR}/figures" \
         "${DESEQ_RUN_DIR}/logs"

LOG_FILE="${DESEQ_RUN_DIR}/logs/deseq2_pipeline.log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "=== 04_deseq2 started at $(date) ==="
echo "DESEQ_RUN_DIR: ${DESEQ_RUN_DIR}"
echo "FEATURECOUNTS_FILE: ${FEATURECOUNTS_FILE}"
echo "QUANT_SAMPLE_TABLE: ${QUANT_SAMPLE_TABLE}"
```

---

## 4. DESeq2 用 R スクリプトの設計

`DESEQ_RUN_DIR/scripts/` に、再実行可能な R スクリプト（例: `run_deseq2_M145.R`）を保存し、そのスクリプトを CLI から実行する形にしてください。
以下は R スクリプトでやるべき処理の要約です。

### 4.1 入力データの読み込み

R スクリプト内で:

- `featureCounts_M145.txt` を読み込み、遺伝子 x サンプルのカウントマトリクスを作成。
- 最初の数行（コメント・ヘッダ）をスキップし、`Geneid` 列を rownames にする。
- `quant_sample_table.tsv` を読み込み、サンプル情報（sample_id, condition, replicate）を data.frame として用意。

### 4.2 デザインとサンプル情報

`colData` として、少なくとも:

- `sample_id`
- `condition`（Factor: M145_1, M145_2, M145_3）
- `replicate`（Factor or integer）

を持つオブジェクトを作る。

例（R 側イメージ）:

```r
coldata <- readr::read_tsv(quant_sample_table) %>%
  dplyr::mutate(
    condition = factor(condition, levels = c("M145_1", "M145_2", "M145_3")),
    replicate = factor(replicate)
  )
```

行順がカウントマトリクスの列と一致するように整合性チェックを行う（一致しない場合は並べ替え）。

### 4.3 DESeqDataSet の作成

```r
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData   = coldata,
  design    = ~ condition
)
```

低カウント遺伝子のフィルタリング（推奨）:

```r
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
```

---

## 5. DESeq2 解析本体

### 5.1 正規化と分散推定

```r
dds <- DESeq(dds)
saveRDS(dds, file = file.path(deseq_run_dir, "rds/dds_raw.rds"))
```

### 5.2 正規化カウントの取得

sizeFactor で正規化された counts:

```r
norm_counts <- counts(dds, normalized = TRUE)
write_tsv(
  as.data.frame(norm_counts) %>%
    tibble::rownames_to_column("gene_id"),
  file.path(deseq_run_dir, "results/normalized_counts_M145.tsv")
)
```

---

## 6. コントラスト設定と差異発現解析

条件が M145_1, M145_2, M145_3 の 3 つなので、少なくとも以下のペアワイズ比較を実施:

1. M145_2 vs M145_1（中間 vs 初期）
2. M145_3 vs M145_1（後期 vs 初期）
3. M145_3 vs M145_2（後期 vs 中間）

### 6.1 LFC shrinkage を使った結果出力

> **注意（apeglm の使い方）**: `lfcShrink()` の `type = "apeglm"` は **`coef` パラメータのみ対応**しており、`contrast` は使えません。`resultsNames(dds)` で得られる係数名を `coef` に指定してください。`contrast` で指定したい場合は `type = "ashr"` を使用してください。

**方法 A（推奨）: `coef` + apeglm**

```r
# resultsNames(dds) で係数名を確認
# 例: "Intercept" "condition_M145_2_vs_M145_1" "condition_M145_3_vs_M145_1"
resultsNames(dds)

# M145_2 vs M145_1
res_2_vs_1 <- lfcShrink(dds, coef = "condition_M145_2_vs_M145_1", type = "apeglm")

# M145_3 vs M145_1
res_3_vs_1 <- lfcShrink(dds, coef = "condition_M145_3_vs_M145_1", type = "apeglm")
```

M145_3 vs M145_2 は直接の `coef` が存在しないため、**reference level を変更して再フィット**するか、`type = "ashr"` + `contrast` を使う:

```r
# 方法 A-1: relevel して再フィット
dds2 <- dds
dds2$condition <- relevel(dds2$condition, ref = "M145_2")
dds2 <- nbinomWaldTest(dds2)
res_3_vs_2 <- lfcShrink(dds2, coef = "condition_M145_3_vs_M145_2", type = "apeglm")

# 方法 A-2: ashr を使う（再フィット不要）
res_3_vs_2 <- lfcShrink(dds, contrast = c("condition", "M145_3", "M145_2"), type = "ashr")
```

**結果の保存**:

```r
for (name in c("2_vs_1", "3_vs_1", "3_vs_2")) {
  res_tbl <- as_tibble(get(paste0("res_", name)), rownames = "gene_id")
  readr::write_tsv(
    res_tbl,
    file.path(deseq_run_dir, paste0("results/DESeq2_M145_", name, ".tsv"))
  )
}
```

各ファイルには以下の列が含まれる（標準 DESeq2 出力）:

- `gene_id`
- `baseMean`
- `log2FoldChange`
- `lfcSE`
- `stat`（apeglm の場合は `svalue`）
- `pvalue`
- `padj`

---

## 7. QC 可視化（PCA / 距離ヒートマップ / 火山図 / 発現ヒートマップ）

### 7.1 rlog or vst 変換

クラスタリングと PCA 用に:

```r
rld <- rlog(dds, blind = FALSE)
saveRDS(rld, file = file.path(deseq_run_dir, "rds/rld.rds"))
```

> サンプル数が多い場合は `vst()` の方が高速。9 サンプルなら `rlog()` で問題ありません。

### 7.2 PCA プロット

condition ごとに色を変えた PCA (PC1 vs PC2) を作成し、PDF/SVG で保存:

```r
pca_data <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw()

ggsave(file.path(deseq_run_dir, "figures/PCA_M145.pdf"), p, width = 6, height = 5)
ggsave(file.path(deseq_run_dir, "figures/PCA_M145.svg"), p, width = 6, height = 5)
```

### 7.3 サンプル間距離ヒートマップ

```r
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- coldata$sample_id
colnames(sampleDistMatrix) <- coldata$sample_id

pheatmap::pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  filename = file.path(deseq_run_dir, "figures/sample_distance_heatmap_M145.pdf")
)
```

（必要に応じて PNG/SVG も出力）

### 7.4 火山図（volcano plot）

各コントラストごとに、log2FC vs -log10(padj) を描き、閾値（例: |log2FC| > 1 & padj < 0.05）で色分けした火山図を作成。

例（M145_3 vs M145_1）:

```r
res_tbl <- as_tibble(res_3_vs_1, rownames = "gene_id") %>%
  mutate(
    neg_log10_padj = -log10(padj),
    sig = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "up",
      padj < 0.05 & log2FoldChange < -1 ~ "down",
      TRUE ~ "ns"
    )
  )

p_volcano <- ggplot(res_tbl, aes(x = log2FoldChange, y = neg_log10_padj, color = sig)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c(up = "red", down = "blue", ns = "grey")) +
  theme_bw() +
  ggtitle("M145_3 vs M145_1")

ggsave(file.path(deseq_run_dir, "figures/volcano_M145_3_vs_1.pdf"), p_volcano, width = 6, height = 5)
```

### 7.5 上位変動遺伝子のヒートマップ

全条件で分散が大きい上位 N 遺伝子（例: 100）を抽出し、rlog 変換値でヒートマップ:

```r
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 100)
mat <- assay(rld)[topVarGenes, ]
mat <- mat - rowMeans(mat)

anno_col <- data.frame(condition = coldata$condition, row.names = coldata$sample_id)

pheatmap::pheatmap(
  mat,
  annotation_col = anno_col,
  filename = file.path(deseq_run_dir, "figures/topVarGenes_heatmap_M145.pdf")
)
```

---

## 8. レポート `deseq2_report_M145.md`

`DESEQ_RUN_DIR` 直下に Markdown レポートを作成し、少なくとも以下を含めてください。

1. **プロジェクト概要とデザイン**
   - 条件（M145_1, M145_2, M145_3）
   - 各条件のサンプル数

2. **入力データ**
   - `featureCounts_M145.txt` の行数（遺伝子数）と列数（サンプル数）
   - フィルタリング後に残った遺伝子数

3. **DESeq2 の解析設定**
   - design 式（`~ condition`）
   - 正規化手法（DESeq2 標準）
   - LFC shrinkage（apeglm）を使用したこと

4. **主な QC 結果**
   - PCA で条件別にクラスターが分かれているか
   - サンプル間距離ヒートマップの傾向

5. **差異発現結果の概要**
   - 各コントラストで padj < 0.05 の遺伝子数（上昇・低下）
   - 上位の代表的遺伝子の例（必要なら gene_id レベルで列挙）

6. **生物学的解釈への導入**
   - M145_1 → M145_3 で Top-10 fraction が上昇していたことと、差異発現結果との整合性

7. **Materials & Methods サマリ（3〜5 行）**
   - 「HISAT2 で整列した reads を featureCounts で gene-level に集約し、DESeq2 vX.Y を用いて `~ condition` をデザインとする差異発現解析を行った」等の文を日本語/英語どちらかで記載。

---

## 9. 成果物チェックリスト

04_deseq2 ステップ完了の目安:

- `${DESEQ_RUN_DIR}/logs/deseq2_pipeline.log`
- `${DESEQ_RUN_DIR}/rds/dds_raw.rds`
- `${DESEQ_RUN_DIR}/rds/rld.rds`
- `${DESEQ_RUN_DIR}/results/normalized_counts_M145.tsv`
- `${DESEQ_RUN_DIR}/results/DESeq2_M145_2_vs_1.tsv`
- `${DESEQ_RUN_DIR}/results/DESeq2_M145_3_vs_1.tsv`
- `${DESEQ_RUN_DIR}/results/DESeq2_M145_3_vs_2.tsv`
- `${DESEQ_RUN_DIR}/figures/PCA_M145.pdf`（および類似図）
- `${DESEQ_RUN_DIR}/figures/sample_distance_heatmap_M145.pdf`
- `${DESEQ_RUN_DIR}/figures/volcano_M145_*.pdf`
- `${DESEQ_RUN_DIR}/figures/topVarGenes_heatmap_M145.pdf`
- `${DESEQ_RUN_DIR}/deseq2_report_M145.md`
