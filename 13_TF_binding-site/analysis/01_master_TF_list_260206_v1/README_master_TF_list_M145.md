# README: master_TF_list_M145.tsv

**作成日**: 2026-02-06
**プロジェクト**: *Streptomyces coelicolor* A3(2) M145 全転写因子マスターテーブル
**解析ディレクトリ**: `13_TF_binding-site/analysis/01_master_TF_list_260206_v1/`

---

## 概要

M145株ゲノム内の全転写因子（TF）、σ因子、およびレギュレーター候補を網羅したマスターテーブル。4つのデータソースを統合して作成。

---

## 使用したリファレンス

- **RefSeq アセンブリ**: GCF_000203835.1 (ASM20383v1)
- **アノテーション**: NCBI RefSeq GCF_000203835.1-RS_2026_01_06
- **リファレンスファイル**:
  - `GCF_000203835.1_ASM20383v1_genomic.fna` （ゲノム配列）
  - `genomic.gff` （遺伝子アノテーション）
  - `protein.faa` （タンパク質配列）

---

## データソース

### 1. Pfam/HMMERドメインスキャン
- **ツール**: HMMER 3.4 (Aug 2023) `hmmscan`
- **データベース**: Pfam-A.hmm（InterPro/Pfam FTP、27,481 HMMs）
  - URL: `https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz`
- **パラメータ**: `--cpu 6 -E 1e-5 --domE 1e-5 --noali`
- **TF関連ドメイン**: HTH系、TetR、GntR、MarR、LysR、AraC、MerR、LuxR、XRE、SARP、WhiB、Sigma70、ECF、REC+Trans_reg_C等
- **結果**: 14,823ドメインヒット、835タンパク質にTF関連ドメイン検出

### 2. Zorro-Aranda et al. 2022
- **論文**: Zorro-Aranda, A., Martínez-Antonio, A., & Freyre-González, J. A. (2022). Curation, inference, and assessment of a globally reconstructed gene regulatory network for *Streptomyces coelicolor* A3(2). *Scientific Reports*, 12, 2704. https://doi.org/10.1038/s41598-022-06658-x
- **データ**: Supplementary Tables (MOESM2) から104個のレギュレーターを抽出
- **ダウンロード**: 自動ダウンロード成功
  - `Zorro-Aranda_2022_MOESM2_SupplementaryTables.xlsx`（618 KB）
  - URL: `https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-022-06658-x/MediaObjects/41598_2022_6658_MOESM2_ESM.xlsx`
- **保存先**: `literature/Zorro-Aranda_2022_MOESM2_SupplementaryTables.xlsx`

### 3. Castro-Melchor et al. 2010
- **論文**: Castro-Melchor, M., Charaniya, S., Karypis, G., Takano, E., & Hu, W.-S. (2010). Genome-wide inference of regulatory networks in *Streptomyces coelicolor*. *BMC Genomics*, 11, 578. https://doi.org/10.1186/1471-2164-11-578
- **データ**: Additional file 3 (遺伝子-オペロンマッピング) + Additional file 4 (692 regulator cistrons) から1,080レギュレーター候補遺伝子を抽出
- **ダウンロード**: 自動ダウンロード成功（16 supplementary files）
  - URL base: `https://static-content.springer.com/esm/art%3A10.1186%2F1471-2164-11-578/MediaObjects/`
- **注意**: 692 cistrons は「少なくとも1つのレギュレーター遺伝子を含むオペロン単位」であり、cistron内の全遺伝子（1,080個）がマークされている。そのため、オペロン同居遺伝子（TF自体ではないもの）も含まれる。`confidence` カラムで区別可能。

### 4. σ因子レビュー
- **主要文献**: Sun, D., Liu, C., Zhu, J., & Liu, W. (2017). Connecting Metabolic Pathways: Sigma Factors in *Streptomyces* spp. *Frontiers in Microbiology*, 8, 2546. https://doi.org/10.3389/fmicb.2017.02546
- **ダウンロード**: PDFの自動ダウンロード成功
  - 保存先: `literature/Sun_2017_sigma_factors_Streptomyces_FrontMicrobiol_8_2546.pdf`
- **データ**: 文献から21個の named σ因子をキュレーション + GFF product列から59個の追加σ因子候補を同定（合計80個）
- **S. coelicolor A3(2) σ因子数**: 64個（Table 1: Group 1: 1, Group 2: 3, Group 3: 10, Group 4/ECF: 50）。本リストの80個はSCP1/SCP2プラスミド遺伝子と2026年RefSeq再注釈による追加を含む。

---

## カラム説明

| カラム名 | 説明 |
|---------|------|
| SCO_ID | 旧locus_tag（SCO####形式）。文献引用に使用 |
| gene_id | 現行RefSeq gene_id（SC_RS#####形式）。解析の主キー |
| gene_name | 遺伝子名（例: sigR, bldD, tetR） |
| contig | chromosome / SCP1 / SCP2 |
| start | 開始位置 |
| end | 終了位置 |
| strand | + / - |
| pfam_TF_flag | Pfam TF関連ドメインの有無（yes/no） |
| tf_related_domains | 検出されたTF関連Pfamドメイン（;区切り） |
| ZorroAranda_regulator_type | Zorro-Aranda 2022の分類（global/sigma/TF/TCS/NA） |
| CastroMelchor_regulator_flag | Castro-Melchor 2010のregulator cistronメンバー（1/0） |
| sigma_flag | σ因子フラグ（yes/no） |
| sigma_group | σ因子グループ（sigma70_group1/group2/group3/ECF/NA） |
| TF_family | TFファミリー分類（TetR, GntR, SARP, etc.） |
| source_flags | データソースフラグ（;区切り） |
| confidence | 信頼度レベル（high/medium/low） |
| product | NCBI RefSeqのproduct名 |

---

## 信頼度レベルの定義

| レベル | 基準 | 推奨用途 |
|--------|------|---------|
| **high** | σ因子、またはPfam TFドメイン + Zorro-Aranda文献、またはZorro-Aranda + Castro-Melchor 両方に登場 | 確実なTF/レギュレーター。全解析に使用可 |
| **medium** | Pfam TFドメインのみ、Zorro-Arandaのみ、Pfam + Castro-Melchor、またはCastro-Melchor + product名がregulator/repressor等 | TF候補として使用。結合サイト解析の主要対象 |
| **low** | Castro-Melchor cistronメンバーのみ（Pfam TFドメインなし、product名に手がかりなし） | レギュレーターcistron内のオペロン同居遺伝子の可能性。フィルタリング推奨 |

**推奨**: `confidence` が `high` または `medium` のエントリ（1,262個）を標準的な解析セットとして使用。

---

## QC サマリ

### 基本統計

| 項目 | 数値 |
|------|------|
| 全CDS数 | 7,996 |
| TF/レギュレーター候補（全体） | 1,705 |
| pfam_TF_flag == yes | 828 |
| ZorroAranda_regulator_type != NA | 104 |
| CastroMelchor_regulator_flag == 1 | 1,053 |
| sigma_flag == yes | 80 |

### 信頼度分布

| Confidence | 数 |
|-----------|-----|
| high | 155 |
| medium | 1,107 |
| low | 443 |
| **合計** | **1,705** |

### σ因子グループ分布

| Group | 数 |
|-------|-----|
| ECF (Group 4) | 45 |
| sigma70_unclassified | 21 |
| sigma70_group3 (alternative) | 10 |
| sigma70_group2 | 3 |
| sigma70_group1 (primary) | 1 |
| **合計** | **80** |

### TFファミリー分布（high + medium）

| TF Family | 数 |
|-----------|-----|
| HTH_other | 180 |
| TetR | 151 |
| Sigma | 135 |
| LuxR | 100 |
| MarR | 69 |
| GntR | 57 |
| LacI | 50 |
| Unclassified_TF | 46 |
| LysR | 41 |
| MerR | 34 |
| TCS_response_regulator | 24 |
| AraC | 16 |
| AsnC | 14 |
| Unknown_TF | 14 |
| SARP | 13 |
| Other_TF | 10 |
| PadR | 9 |
| DeoR | 7 |
| Lsr2 | 1 |

---

## スクリプト一覧

| スクリプト | 内容 |
|-----------|------|
| `scripts/01_extract_gene_info.py` | GFFから全遺伝子基本情報を抽出 |
| `scripts/02_parse_castro_melchor.py` | Castro-Melchor 2010 のcistron→SCO_IDマッピング |
| `scripts/03_create_sigma_factor_list.py` | σ因子リスト作成（文献 + GFFアノテーション） |
| `scripts/04_standardize_zorro_aranda.py` | Zorro-Aranda 2022 データの標準化 |
| `scripts/05_parse_pfam_and_integrate.py` | Pfamパース + 全データ統合 |
| `scripts/06_add_confidence_and_finalize.py` | 信頼度レベル追加 + 最終出力 |

---

## 中間ファイル

| ファイル | 内容 |
|---------|------|
| `intermediate/M145_gene_basic_info.tsv` | 全8,083遺伝子の基本情報 |
| `intermediate/M145_pfam_domtblout.txt` | HMMER domtblout結果（14,823ヒット） |
| `intermediate/ZorroAranda2022_regulators_standardized.tsv` | Zorro-Aranda標準化済みレギュレーター |
| `intermediate/CastroMelchor2010_regulators_dedup.tsv` | Castro-Melchor重複除去済みレギュレーター |
| `intermediate/SigmaFactors_M145_from_lit.tsv` | σ因子リスト |
| `intermediate/M145_all_genes_with_TF_annotations.tsv` | 全遺伝子のTFアノテーション付き完全テーブル |

---

## 参考文献

1. Zorro-Aranda, A., Martínez-Antonio, A., & Freyre-González, J. A. (2022). Curation, inference, and assessment of a globally reconstructed gene regulatory network for *Streptomyces coelicolor* A3(2). *Scientific Reports*, 12, 2704. https://doi.org/10.1038/s41598-022-06658-x

2. Castro-Melchor, M., Charaniya, S., Karypis, G., Takano, E., & Hu, W.-S. (2010). Genome-wide inference of regulatory networks in *Streptomyces coelicolor*. *BMC Genomics*, 11, 578. https://doi.org/10.1186/1471-2164-11-578

3. Sun, D., Liu, C., Zhu, J., & Liu, W. (2017). Connecting Metabolic Pathways: Sigma Factors in *Streptomyces* spp. *Frontiers in Microbiology*, 8, 2546. https://doi.org/10.3389/fmicb.2017.02546

4. Bentley, S. D., et al. (2002). Complete genome sequence of the model actinomycete *Streptomyces coelicolor* A3(2). *Nature*, 417, 141–147.

---

*最終更新: 2026-02-06*
