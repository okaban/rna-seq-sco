# 02_ALIGNMENT プロンプト（M145 RNA-seq, STAR アラインメント）

あなたは、*Streptomyces coelicolor* A3(2) M145 の RNA-seq データに対して、STAR を用いたリファレンスゲノムへのアラインメント（02_alignment ステップ）を実行するバイオインフォマティクスエージェントです。
すでに `COMMON_PROMPT_M145.md` と `01_qc_M145_prompt.md` に基づいて QC が完了している前提で動いてください。

---

## 0. このタスクでやってほしいこと

- QC 済み fastq（RAW または trimmed）の情報をもとに、STAR で M145 リファレンスゲノムにアラインメントを実施する。
- リファレンス（GCF_000203835.1）から STAR 用インデックスを作成する。
- 各サンプルについて、ソート済み BAM（および必要ならインデックス .bai）を作成する。
- アラインメントに関する QC 指標（マッピング率など）を集計し、Markdown レポートにまとめる。
- 実行ログをファイルに集約し、後続解析（featureCounts, FADU, DESeq2）にそのまま渡せる状態にする。
- ClaudeCode のトークン消費を抑えるため、詳細ログはファイルに書き出し、チャットには要約だけ返すこと。

---

## 1. 前提とパス・変数定義

### 1.1 共通前提

- QC ステップはすでに実行済みであり、直近の run ディレクトリ（例:
  `/Users/okaban/bioinfo/rna-seq/01_qc/analysis/01_qc_260127_v1`）に `sample_table.tsv` と `qc_report_M145.md` が存在する。
- `sample_table.tsv` から、どの fastq（RAW/trimmed）を使うかを決定できる（必要に応じて `qc_report_M145.md` を参照）。

### 1.2 パスと変数

以下をデフォルトとします（必要に応じて冒頭で上書き可能）:

```bash
# アラインメントステップ用のベースディレクトリ
ALIGN_ROOT="/Users/okaban/bioinfo/rna-seq/02_alignment"

# QC の最新 run ディレクトリ（必要に応じて更新）
QC_ROOT="/Users/okaban/bioinfo/rna-seq/01_qc"
QC_RUN_DIR="${QC_ROOT}/analysis/01_qc_260127_v1"  # ← 実際の run ID に合わせて更新

# QC で作成されたサンプル表
QC_SAMPLE_TABLE="${QC_RUN_DIR}/sample_table.tsv"

# QC のレポート（RAW vs trimmed の判断に使用）
QC_REPORT="${QC_RUN_DIR}/qc_report_M145.md"

# リファレンスディレクトリ（GCF_000203835.1）
REF_DIR="/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1"

# リファレンスファイル（STAR インデックス作成に使用）
REF_FASTA="${REF_DIR}/GCF_000203835.1_ASM20383v1_genomic.fna"
REF_GTF="${REF_DIR}/genomic.gtf"   # または genomic.gff（ツールに応じて）

# 出力ルート
ALIGN_OUT_ROOT="${ALIGN_ROOT}/analysis"

# スクリプト格納先
ALIGN_SCRIPT_DIR="${ALIGN_ROOT}/scripts"

# STAR インデックス出力先
STAR_INDEX_DIR="${ALIGN_ROOT}/star_index"

# run ID（実行日＋バージョン）
ALIGN_RUN_DATE=$(date +%y%m%d)        # 例: 260127
ALIGN_RUN_ID="02_alignment_${ALIGN_RUN_DATE}_v1"
ALIGN_RUN_DIR="${ALIGN_OUT_ROOT}/${ALIGN_RUN_ID}"

# 利用スレッド数（MacBook Pro のコア数に応じて調整）
NTHREADS=6
```

`ALIGN_RUN_DIR` の標準構造:

```text
${ALIGN_RUN_DIR}/
├── star_index/       # STAR genome index（必要ならこの run にコピー or リンク）
├── bam/              # ソート済み BAM と .bai
├── logs/             # STAR の Log.out / Log.final.out など
├── summary/          # マッピング統計の集約 TSV / Markdown
└── pipeline.log      # 全体ログ
```

---

## 2. Conda 環境の確認とツール準備

RNA-seq 用 Conda 環境名は `rna-seq` を想定。

- `rna-seq` 環境が存在するか確認。
- **存在しない場合**:
  01_qc のタイミングで作っているはずだが、存在しない場合は新規作成。
- **存在する場合**:
  `conda activate rna-seq`。
- STAR, samtools が使用可能か確認。
- 足りなければインストール。

```bash
# 環境一覧
conda env list

# 必要なら（まだなければ）環境作成
conda create -n rna-seq -c bioconda -c conda-forge \
  python=3.9 \
  star \
  samtools

# アクティベート
conda activate rna-seq

# ツール確認
command -v STAR     >/dev/null || echo "STAR missing"
command -v samtools >/dev/null || echo "samtools missing"

# 足りなければ追加インストール
conda install -c bioconda -c conda-forge star samtools

# 環境 export（任意）
mkdir -p "${ALIGN_SCRIPT_DIR}"
conda env export > "${ALIGN_SCRIPT_DIR}/environment_rnaseq_alignment_${ALIGN_RUN_DATE}.yml"
```

---

## 3. セットアップ & ログ開始

```bash
# 出力ディレクトリ作成
mkdir -p "${ALIGN_RUN_DIR}/bam" \
         "${ALIGN_RUN_DIR}/logs" \
         "${ALIGN_RUN_DIR}/summary"

# STAR インデックスディレクトリ
STAR_INDEX_DIR="${ALIGN_RUN_DIR}/star_index"

mkdir -p "${STAR_INDEX_DIR}"

# ログ設定
LOG_FILE="${ALIGN_RUN_DIR}/pipeline.log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "=== 02_alignment started at $(date) ==="
echo "ALIGN_ROOT: ${ALIGN_ROOT}"
echo "ALIGN_RUN_DIR: ${ALIGN_RUN_DIR}"
echo "QC_RUN_DIR: ${QC_RUN_DIR}"
echo "REF_FASTA: ${REF_FASTA}"
echo "REF_GTF: ${REF_GTF}"
conda info --envs | grep '*' || true
uname -a || true
```

---

## 4. STAR インデックス作成

すでに同じ FASTA/注釈に対する STAR インデックスがある場合は再利用してよいですが、ここでは「この run 用にインデックスを作る」手順を明示します。

```bash
STAR \
  --runMode genomeGenerate \
  --genomeDir "${STAR_INDEX_DIR}" \
  --genomeFastaFiles "${REF_FASTA}" \
  --sjdbGTFfile "${REF_GTF}" \
  --runThreadN "${NTHREADS}"
```

メモリ不足が起きる場合は、`--genomeSAsparseD` の調整などを検討（必要になったときで良い）。

---

## 5. 使用する fastq の決定（RAW vs trimmed）

`QC_SAMPLE_TABLE` と `QC_REPORT` を読み、以下の方針で fastq を選択:

- **QC レポートで「トリミング不要」とされている場合**:
  生データ（RAW fastq）を使用。
- **QC レポートで「トリミングを実施し、trimmed fastq を使用する」とされている場合**:
  `QC_RUN_DIR/fastp/` 内の `*_trimmed_R1.fastq.gz`, `*_trimmed_R2.fastq.gz` を使用。

このロジックに基づき、`ALIGN_RUN_DIR/summary/` に、実際にアライメントに使用する fastq の一覧をまとめた `alignment_input_table.tsv` を作成してください:

```text
sample_id	condition	replicate	R1_fastq	R2_fastq	source
M145_1_1	M145_1	1	/path/to/..._R1.fastq.gz	/path/to/..._R2.fastq.gz	RAW
M145_1_2	M145_1	2	/path/to/..._trimmed_R1.fastq.gz	/path/to/..._trimmed_R2.fastq.gz	TRIMMED
...
```

---

## 6. STAR によるアラインメント実行

`alignment_input_table.tsv` をループし、各サンプルについて STAR を実行します。
以下は 1 サンプルの例です:

```bash
STAR \
  --genomeDir "${STAR_INDEX_DIR}" \
  --readFilesIn <R1_fastq> <R2_fastq> \
  --readFilesCommand zcat \
  --runThreadN "${NTHREADS}" \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix "${ALIGN_RUN_DIR}/bam/<sample_id>." \
  --quantMode TranscriptomeSAM \
  --outSAMattributes NH HI AS nM MD \
  --outFilterMultimapNmax 1 \
  --outFilterMismatchNmax 10
```

推奨設定（細かい調整は今後でも良い）:

- `--outSAMtype BAM SortedByCoordinate`
  → ソート済み BAM を直接出力。
- `--outFilterMultimapNmax 1`
  → ユニークマッピングのみを残す（バクテリアでよく使う設定）。
- `--outFilterMismatchNmax 10`
  → 許容ミスマッチ数（read 長や品質に応じて調整可）。
- `--quantMode TranscriptomeSAM`
  → 後で transcriptome ベースの定量に使いたい場合。

STAR が生成する主なファイル（サンプル `<sample_id>`）:

- `${ALIGN_RUN_DIR}/bam/<sample_id>.Aligned.sortedByCoord.out.bam`
- `${ALIGN_RUN_DIR}/bam/<sample_id>.Log.out`
- `${ALIGN_RUN_DIR}/bam/<sample_id>.Log.final.out`
- `${ALIGN_RUN_DIR}/bam/<sample_id>.SJ.out.tab`
- `${ALIGN_RUN_DIR}/bam/<sample_id>.Aligned.toTranscriptome.out.bam`（quantMode 使用時）

実行後、`samtools index` で BAM にインデックスをつけておくと便利です:

```bash
samtools index "${ALIGN_RUN_DIR}/bam/<sample_id>.Aligned.sortedByCoord.out.bam"
```

全サンプル分の BAM と BAI を作成してください。

---

## 7. アラインメント QC の集計

STAR の `Log.final.out` は、マッピング率などの QC 情報を含んでいます。
各サンプルの `Log.final.out` から主要指標（例: 総リード数、ユニークマッピング率、多重マッピング率など）を抽出し、`ALIGN_RUN_DIR/summary/alignment_stats.tsv` を作成してください。

例:

```text
sample_id	total_reads	uniquely_mapped	uniquely_mapped_percent	multi_mapped_percent	unmapped_percent
M145_1_1	12345678	12000000	97.2	1.1	1.7
...
```

---

## 8. アラインメントレポート `alignment_report_M145.md`

`ALIGN_RUN_DIR` 直下に Markdown レポートを作成してください。
推奨構成:

- プロジェクト概要
- 使用したリファレンス（FASTA / GTF のパス）
- 使用した fastq（RAW/trimmed）の一覧（`alignment_input_table.tsv` の要約）
- STAR の実行パラメータ（代表例）
- アラインメント QC 結果の概要（`alignment_stats.tsv` の要約）
- M145 特有の注意点（高GC, オペロンなど）とアラインメントへの影響
- 次ステップ（featureCounts/FADU → DESeq2）への接続メモ
- Materials & Methods にそのまま使える 3〜5 行のアラインメント手順サマリ（英語/日本語どちらでも）

---

## 9. 出力物のチェックリスト

このステップが完了したと判断するには、少なくとも以下が揃っている必要があります:

- `ALIGN_RUN_DIR/pipeline.log`
- `ALIGN_RUN_DIR/star_index/`（STAR インデックス）
- `ALIGN_RUN_DIR/bam/<sample_id>.Aligned.sortedByCoord.out.bam`（全9サンプル分）
- `ALIGN_RUN_DIR/bam/<sample_id>.Aligned.sortedByCoord.out.bam.bai`（インデックス）
- `ALIGN_RUN_DIR/bam/<sample_id>.Log.final.out`（全サンプル分）
- `ALIGN_RUN_DIR/summary/alignment_input_table.tsv`
- `ALIGN_RUN_DIR/summary/alignment_stats.tsv`
- `ALIGN_RUN_DIR/alignment_report_M145.md`
