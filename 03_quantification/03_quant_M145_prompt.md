# 03_QUANTIFICATION プロンプト（M145 RNA-seq, featureCounts/FADU によるカウント）

あなたは、*Streptomyces coelicolor* A3(2) M145 の RNA-seq データに対して、
「HISAT2 で整列済み BAM から遺伝子ごとのカウントを取得する」ステップ（03_quantification）を実行するバイオインフォマティクスエージェントです。

すでに以下が完了している前提で動いてください:

- `COMMON_PROMPT_M145.md` に基づくプロジェクト設定
- `01_qc_M145_prompt.md` による QC 完了
- `02_alignment` ステップにより、HISAT2 アラインメント済み BAM と `alignment_report_M145.md` が揃っている

---

## 0. このタスクでやってほしいこと

- 02_alignment で得られた 9 サンプルの sorted BAM（HISAT2）を用いて、**遺伝子ごとの read count** を取得する。
- メインは **featureCounts** を用いた gene-level カウント、オプションとして **FADU**（バクテリア用）によるカウントも試せるようにする。
- 得られたカウントを DESeq2 にそのまま渡せる形式（行＝遺伝子、列＝サンプル）にまとめる。
- カウント結果の簡単な QC とサマリ（総カウント分布、サンプルごとの library size など）を Markdown レポートにまとめる。
- 実行ログをファイルに集約し、ClaudeCode のトークン消費は最小限に抑える（詳細はログ・TSV・図に、チャットにはサマリのみ）。

---

## 1. パスと変数定義

以下をデフォルトとして使用してください（必要に応じて上書き可）。

```bash
# カウントステップ用のベースディレクトリ
QUANT_ROOT="/Users/okaban/bioinfo/rna-seq/03_quantification"

# アラインメントステップのベース
ALIGN_ROOT="/Users/okaban/bioinfo/rna-seq/02_alignment"

# 直近の alignment run ディレクトリ（必要に応じて更新）
ALIGN_RUN_DIR="${ALIGN_ROOT}/analysis/02_alignment_260127_v1"  # 実際の RUN ID に合わせて変更

# BAM ディレクトリ
BAM_DIR="${ALIGN_RUN_DIR}/bam"

# サンプル一覧（alignment_input_table.tsv を再利用）
ALIGN_SAMPLE_TABLE="${ALIGN_RUN_DIR}/summary/alignment_input_table.tsv"

# リファレンスディレクトリ（GCF_000203835.1）
REF_DIR="/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1"

# アノテーション（featureCounts 用）
REF_GTF="${REF_DIR}/genomic.gtf"   # 必要に応じて genomic.gff に変更可

# 出力ルート
QUANT_OUT_ROOT="${QUANT_ROOT}/analysis"

# スクリプト格納先
QUANT_SCRIPT_DIR="${QUANT_ROOT}/scripts"

# run ID（実行日＋バージョン）
QUANT_RUN_DATE=$(date +%y%m%d)          # 例: 260128
QUANT_RUN_ID="03_quant_${QUANT_RUN_DATE}_v1"
QUANT_RUN_DIR="${QUANT_OUT_ROOT}/${QUANT_RUN_ID}"

# 利用スレッド数（MacBook Pro のコア数に応じて調整）
NTHREADS=6
```

`QUANT_RUN_DIR` の標準構造:

```text
${QUANT_RUN_DIR}/
├── counts/         # featureCounts / FADU のカウント結果
├── logs/           # featureCounts / FADU のログ
├── summary/        # カウント QC・サマリ
├── figures/        # 補足図（ライブラリサイズなど）
└── pipeline.log    # 全体ログ
```

---

## 2. Conda 環境とツール準備

RNA-seq 用 Conda 環境名は `rna-seq` を想定します。

- `conda env list` で `rna-seq` 環境の有無を確認。
- 既に存在し、featureCounts（subread）、samtools が使える場合はそのまま使用。
- 足りない場合のみインストール。

```bash
# 環境一覧
conda env list

# アクティベート
conda activate rna-seq

# ツール存在チェック
command -v featureCounts >/dev/null || echo "featureCounts missing"
command -v samtools       >/dev/null || echo "samtools missing"

# 足りないツールのインストール
conda install -c bioconda -c conda-forge subread samtools

# FADU を試す場合（任意）：conda パッケージがなければ、後で GitHub から配置する（ここでは必須にしない）
# 例: git clone https://github.com/IGS/FADU.git （手動セットアップ前提）
```

環境を更新した場合は、再現性のために:

```bash
mkdir -p "${QUANT_SCRIPT_DIR}"
conda env export > "${QUANT_SCRIPT_DIR}/environment_rnaseq_quant_${QUANT_RUN_DATE}.yml"
```

---

## 3. セットアップ & ログ開始

```bash
# 出力ディレクトリ作成
mkdir -p "${QUANT_RUN_DIR}/counts" \
         "${QUANT_RUN_DIR}/logs" \
         "${QUANT_RUN_DIR}/summary" \
         "${QUANT_RUN_DIR}/figures"

# ログ設定
LOG_FILE="${QUANT_RUN_DIR}/pipeline.log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "=== 03_quantification started at $(date) ==="
echo "QUANT_ROOT: ${QUANT_ROOT}"
echo "QUANT_RUN_DIR: ${QUANT_RUN_DIR}"
echo "ALIGN_RUN_DIR: ${ALIGN_RUN_DIR}"
echo "BAM_DIR: ${BAM_DIR}"
echo "REF_GTF: ${REF_GTF}"
conda info --envs | grep '*' || true
uname -a || true
```

---

## 4. 対象 BAM とサンプル情報の確認

### 4.1 BAM ファイル一覧

`BAM_DIR` 内の BAM ファイル（`*.Aligned.sortedByCoord.out.bam`）を列挙し、サンプル ID と紐付けてください。

例: `M145_1_1.Aligned.sortedByCoord.out.bam` → sample_id: `M145_1_1`

### 4.2 サンプル表の作成

`ALIGN_SAMPLE_TABLE` と BAM ファイルを照合し、最終的なカウント対象を `QUANT_RUN_DIR/summary/quant_sample_table.tsv` にまとめてください。

推奨形式:

```text
sample_id	condition	replicate	bam_path
M145_1_1	M145_1	1	/path/to/M145_1_1.Aligned.sortedByCoord.out.bam
...
```

不足・不一致があれば、ログとサマリに明記し、どう修正するべきかコメントしてください。

---

## 5. featureCounts によるカウント（メイン）

### 5.1 コマンドの基本形

9 サンプルの BAM を一度に与えて、遺伝子ごとのカウントマトリクスを作成します。

```bash
featureCounts \
  -T "${NTHREADS}" \
  -p --countReadPairs \
  -B \
  -C \
  -s 0 \
  -t gene \
  -g gene_id \
  -a "${REF_GTF}" \
  -o "${QUANT_RUN_DIR}/counts/featureCounts_M145.txt" \
  ${BAM_DIR}/*.bam \
  2> "${QUANT_RUN_DIR}/logs/featureCounts_M145.log"
```

オプションの意味（paired-end 原核生物想定）:

- `-T`: スレッド数
- `-p --countReadPairs`: ペアエンドとしてフラグメント単位でカウント（subread ≥2.0 では `--countReadPairs` が必須）
- `-B`: 正しくペアリングされたリードのみカウント
- `-C`: chimeric フラグメントを除外
- `-s 0`: ストランド非特異（unstranded）。ライブラリがストランド特異的な場合は `-s 1`（forward）または `-s 2`（reverse）に変更すること
- `-t gene`: GTF の feature type として `gene` を使用。本 GTF は `gene`（8,276 件）、`CDS`（8,200 件）、`exon`（87 件＝rRNA/ncRNA のみ）を含むため、原核生物の全遺伝子カウントには **`-t gene` が適切**
- `-g gene_id`: `gene_id` 属性ごとに集約（例: `SC_RS02065`）。`locus_tag` でも同値

> **注意**: `-s`（strandedness）の値が不正確だとカウントが大幅に変わります。
> 実行前に `featureCounts` の `-s 0`, `-s 1`, `-s 2` を 1 サンプルで試し、assign されるリード数が最大になる設定を選択してください。
> 判断がつかない場合は、`salmon` や `infer_experiment.py`（RSeQC）で strandedness を推定するのも有効です。

### 5.2 出力形式

`featureCounts_M145.txt` はヘッダ付きのタブ区切りテキストで、行が遺伝子、列がサンプルになっています。
このファイルを後で R/DESeq2 に入力する予定です（04_deseq2 ステップ）。

---

## 6. （任意）FADU によるバクテリア特化カウント

原核生物特有のオペロン構造や短い遺伝子の見落としを補うため、FADU を使うと、HTSeq/featureCounts ではカウントがつかない小遺伝子を検出しやすくなります。

- GitHub: https://github.com/IGS/FADU

使う場合:

1. BAM（HISAT2）と GFF を入力として FADU を実行し、
2. `QUANT_RUN_DIR/counts/FADU_M145.tsv` を別に出力、
3. 将来、featureCounts との比較に使えるようにする。

※ 今回のプロンプトでは **必須ではなく、余力があれば導入する項目** に留めてください。

---

## 7. カウント結果の QC とサマリ

### 7.1 ライブラリサイズとカウント分布

`featureCounts_M145.txt` を読み込み、以下の集計を行って `QUANT_RUN_DIR/summary/counts_summary.tsv` に書き出してください:

各サンプルの:

- 総カウント数（total read count）
- 非ゼロ遺伝子数
- 上位 10 遺伝子のカウント合計割合 など

例:

```text
sample_id	total_counts	n_genes_nonzero	top10_genes_fraction
M145_1_1	xxx	yyyy	0.35
...
```

### 7.2 図の作成（任意だが推奨）

DESeq2 前の QC として、`QUANT_RUN_DIR/figures/` に以下のような PDF/SVG を保存:

- 各サンプルのライブラリサイズ（total_counts）の棒グラフ
- 各サンプルの log10(count+1) の分布（箱ひげ図 or density）

ファイル名例:

```text
${QUANT_RUN_DIR}/figures/library_size_barplot.pdf
${QUANT_RUN_DIR}/figures/library_size_barplot.svg
${QUANT_RUN_DIR}/figures/count_distribution_boxplot.pdf
```

---

## 8. カウントレポート `quant_report_M145.md`

`QUANT_RUN_DIR` 直下に Markdown レポートを作成し、以下を含めてください:

1. プロジェクト概要
2. 入力 BAM の概要（サンプル数、条件、HISAT2 で整列済みであること）
3. 使用したアノテーション（GTF/GFF のパスとバージョン）
4. featureCounts の実行コマンドと主要オプション
5. 各サンプルのライブラリサイズと非ゼロ遺伝子数の概要（`counts_summary.tsv` の要約）
6. M145（高GC・オペロン）に特有の注意点と、featureCounts/FADU の位置づけ
7. 次ステップ（DESeq2）への接続メモ
8. Materials & Methods に流用できる 3〜5 行のサマリ

---

## 9. 成果物チェックリスト

このステップが完了したと判断するには、少なくとも以下が揃っている必要があります:

- `${QUANT_RUN_DIR}/pipeline.log`
- `${QUANT_RUN_DIR}/counts/featureCounts_M145.txt`
- （任意）`${QUANT_RUN_DIR}/counts/FADU_M145.tsv`
- `${QUANT_RUN_DIR}/summary/quant_sample_table.tsv`
- `${QUANT_RUN_DIR}/summary/counts_summary.tsv`
- `${QUANT_RUN_DIR}/quant_report_M145.md`
- （任意）`${QUANT_RUN_DIR}/figures/*.pdf`, `*.svg`
