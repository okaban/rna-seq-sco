# 解析レポート: *Streptomyces coelicolor* A3(2) M145 転写因子マスターリスト & Binding Site統合

**作成日**: 2026-02-06
**解析ディレクトリ**: `13_TF_binding-site/analysis/01_master_TF_list_260206_v1/`
**実行環境**: conda `rna-seq` (Python 3.11.14, pandas 3.0.0, HMMER 3.4, MEME Suite 5.5.9)

---

## 1. 解析の目的

*Streptomyces coelicolor* A3(2) M145株ゲノムに存在する全転写因子（TF）、σ因子、および転写制御因子候補を網羅的にリストアップし、さらに各TFの結合サイト（binding site）の有無・配列・ゲノム上の位置を紐付けた統合データベースを構築する。

本解析は2つの段階からなる：

- **Phase 1**: 複数データソースを統合したTF候補マスターリストの作成
- **Phase 2**: 実験的・計算的binding site情報の収集と統合

---

## 2. 使用したリファレンス

| 項目 | 値 |
|------|-----|
| RefSeqアセンブリ | GCF_000203835.1 (ASM20383v1) |
| アノテーション | NCBI RefSeq GCF_000203835.1-RS_2026_01_06 |
| ゲノム構成 | chromosome (NC_003888.3, 8.67 Mb) + SCP1 (NC_003903.1, 356 kb) + SCP2 (NC_003904.1, 31 kb) |
| GC含量 | ~72% |
| 全遺伝子数 | 8,083 (protein_coding: 7,996) |

---

## 3. Phase 1: TF候補マスターリスト作成

### 3.1 データソース

#### Source 1: Pfam/HMMERドメインスキャン

- **ツール**: HMMER 3.4 `hmmscan`
- **データベース**: Pfam-A.hmm (27,481 HMMs, InterPro/Pfam FTP)
- **パラメータ**: `--cpu 6 -E 1e-5 --domE 1e-5 --noali`
- **結果**: 14,823ドメインヒット → **835タンパク質**にTF関連ドメイン検出

TF関連ドメインとして以下を使用: HTH系 (HTH_1, HTH_3, HTH_5, HTH_6, HTH_7, HTH_8, HTH_11, HTH_12, HTH_13, HTH_17, HTH_18, HTH_19, HTH_20, HTH_21, HTH_22, HTH_23, HTH_24, HTH_25, HTH_26, HTH_27, HTH_28, HTH_29, HTH_31, HTH_32, HTH_34, HTH_36, HTH_38, HTH_40, HTH_42, HTH_43, HTH_45, HTH_46, HTH_47, HTH_AraC), TetR_N, GntR, MarR, MarR_2, LysR_substrate, AraC_binding, AraC_binding_2, MerR, MerR_1, MerR-HTH, LuxR_C_like, Trans_reg_C (with REC), XRE_family, BTAD, SARP, WhiB, Sigma70_r2, Sigma70_r3, Sigma70_r4, Sigma70_r4_2, Sigma70_ECF, SigmaE_PrsW, PadR, DeoR, AsnC_trans_reg, Lsr2, IclR, FUR, CodY, CRP

#### Source 2: Zorro-Aranda et al. 2022

- **論文**: *Scientific Reports*, 12, 2704
- **データ**: MOESM2 Supplementary Tables → **104レギュレーター**
  - 内訳: TF 56, global 32, sigma 9, TCS 7

#### Source 3: Castro-Melchor et al. 2010

- **論文**: *BMC Genomics*, 11, 578
- **データ**: Additional File 3 (遺伝子→オペロンマッピング) + Additional File 4 (692 regulator cistrons)
- **結果**: 692 cistrons → **1,080ユニークレギュレーター関連遺伝子**
- **注意**: cistron（オペロン）単位でのマーキングのため、TF自体ではないオペロン同居遺伝子も含む

#### Source 4: σ因子レビュー

- **主要文献**: Sun et al. 2017, *Front. Microbiol.*, 8, 2546
- **データ**: 文献から21個のnamed σ因子 + GFF product列から59個の追加σ因子候補
- **結果**: 合計 **80 σ因子**
  - Group 1 (primary): 1, Group 2: 3, Group 3: 10, ECF (Group 4): 45, Unclassified: 21

### 3.2 統合方法

1. GFFから全遺伝子の基本情報 (gene_id, SCO_ID, contig, start, end, strand) を抽出
2. hmmscanで全protein_codingに対しPfam-Aスキャンを実行し、TF関連ドメイン保有遺伝子をフラグ付け
3. ドメイン構成に基づくTFファミリー分類（優先度ベースのルールで割り当て）
4. Zorro-Aranda、Castro-Melchor、σ因子リストをgene_idベースでleft join
5. いずれかのソースでTF/レギュレーターとして同定された遺伝子のみ抽出
6. 信頼度レベル（high/medium/low）を付与

### 3.3 信頼度レベル

| レベル | 基準 | 数 |
|--------|------|-----|
| **high** | σ因子、またはPfamドメイン + Zorro-Aranda文献、またはZorro-Aranda + Castro-Melchor 両方 | 155 |
| **medium** | Pfamドメインのみ、Zorro-Arandaのみ、Pfam + Castro-Melchor、またはCastro-Melchor + product名がregulator等 | 1,107 |
| **low** | Castro-Melchor cistronメンバーのみ（Pfam TFドメインなし、product名に手がかりなし） | 443 |

### 3.4 結果サマリ

| 項目 | 数値 |
|------|------|
| **TF/レギュレーター候補合計** | **1,705** |
| pfam_TF_flag == yes | 828 |
| ZorroAranda登録あり | 104 |
| CastroMelchor_regulator_flag == 1 | 1,053 |
| sigma_flag == yes | 80 |
| 推奨解析セット (high + medium) | 1,262 |

#### ゲノムロケーション分布

| Contig | TF候補数 |
|--------|---------|
| chromosome | 1,684 |
| SCP1 | 20 |
| SCP2 | 1 |

#### ソースフラグ分布

| ソース組み合わせ | 数 |
|-----------------|-----|
| CastroMelchor only | 734 |
| pfam only | 612 |
| pfam + CastroMelchor | 189 |
| pfam + ZorroAranda + CastroMelchor | 60 |
| pfam + ZorroAranda | 38 |
| sigma only | 36 |
| sigma + pfam | 14 |
| sigma + pfam + ZorroAranda | 7 |
| sigma + pfam + CastroMelchor | 6 |
| ZorroAranda only | 5 |
| sigma + pfam + ZorroAranda + CastroMelchor | 2 |
| CastroMelchor + ZorroAranda | 1 |
| sigma + CastroMelchor | 1 |

#### TFファミリー分布 (上位15)

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

#### High confidence TFのファミリー分布

| TF Family | 数 |
|-----------|-----|
| Sigma | 79 |
| Unclassified_TF | 38 |
| TCS_response_regulator | 12 |
| TetR | 9 |
| GntR | 5 |
| LuxR | 4 |
| LysR | 3 |
| MarR | 2 |
| AraC | 1 |
| DeoR | 1 |
| HTH_other | 1 |

---

## 4. Phase 2: Binding Site情報の統合

### 4.1 データソース

#### BS Source 1: ユーザーキュレーションモチーフ + FIMO

文献から専門家がキュレーションした16個のTF binding motifコンセンサス配列をIUPACコードで記述し、MEME形式に変換後、FIMOでゲノムワイドスキャンを実施した。

**キュレーションされたモチーフ一覧:**

| TF名 | SCO_ID | コンセンサス配列 | モチーフ長 | 構造 | エビデンス |
|------|--------|----------------|-----------|------|----------|
| GlnR | SCO4159 | GTNACNNNNNNGANAC | 16 bp | single_box | DNase I footprint; EMSA |
| PhoP | SCO7637 | GTTCACCNNNNGTTCACC | 18 bp | direct_repeat | ChIP-chip; DNase I footprint; EMSA |
| DasR | SCO5231 | ANTGGTCTAGACCANT | 16 bp | palindrome | ChIP-seq; EMSA |
| AfsQ1 | SCO4907 | GTNACNNNNNNGTNAC | 16 bp | direct_repeat | DNase I footprint |
| DraR | SCO3063 | AMAAWYMAKCA | 11 bp | consensus | DNase I footprint; EMSA |
| Crp | SCO3571 | TGTGANNNNNNTCACA | 16 bp | palindrome | ChIP-chip; homology |
| ScbR | SCO6265 | TTTGGNNNNNNNNNNNNNNNNNNCCAAA | 28 bp | inverted_repeat | DNase I footprint; MEME |
| AbrC3 | SCO4596 | GAASGSGRMS | 10 bp | consensus | EMSA; ChIP-qPCR |
| BldD | SCO1489 | AGTGANNNNNNTCACT | 16 bp | palindrome | ChIP-chip; DNase I footprint; EMSA |
| AdpA | SCO2792 | TGGCSNGWWY | 10 bp | consensus | EMSA; DNase I footprint; lacZ reporter |
| HrdB | SCO5820 | TTGACANNNNNNNNNNNNNNNNNTAGAAT | 28 bp | promoter -35/-10 | promoter compilation; IVT |
| SigR | SCO5216 | GGAACNNNNNNNNNNNNNNNNCGTT | 24 bp | promoter -35/-10 | IVT; ChIP-chip |
| SigE | SCO5147 | GGAACNNNNNNNNNNNNNNNNNGTT | 24 bp | promoter -35/-10 | ChIP-seq; IVT |
| SigB | SCO0600 | GNNTNNNNNNNNNNNNNNNNNGGGTAC | 26 bp | promoter -35/-10 | consensus-directed search |
| ArgR | SCO1576 | TGAATAANNNNNNTTATTCA | 20 bp | palindrome | EMSA; DNase I footprint |
| NdgR | SCO5552 | BTYCANYWSYBGGAC | 15 bp | palindrome | ChIP-seq; EMSA |

**FIMO実行パラメータ:**
- 閾値: p-value < 1e-4
- バックグラウンド頻度: A=0.14, C=0.36, G=0.36, T=0.14 (GC ~72%)
- 結果: 58,460ヒット（全モチーフ合計）
- 統合時フィルタ: q-value < 0.01 → **120ヒット**

**FIMO ヒット内訳 (q < 0.01):**

| TF | ヒット数 |
|----|---------|
| SigB | 68 |
| SigE | 35 |
| SigR | 16 |
| PhoP | 1 |

他のモチーフ（GlnR, DasR, BldD等）はq < 0.01では検出されなかった。これはGC-richゲノムにおいて短いまたは縮退度の高いモチーフの偽陽性が多く、厳しいq-value閾値でフィルタリングされたことによる。p-value閾値のみでは多数のヒットがあり、モチーフ探索（p < 1e-4での全58,460ヒット）はFIMO結果ディレクトリに保存されている。

#### BS Source 2: Zorro-Aranda 2022 MEME推定binding site

- **データ**: MOESM3 Supplementary File 2, MEME_BSシート
- **内容**: MEME de novoモチーフ発見によるTF-target gene間の推定binding site
- **レコード数**: 27,919 (72 TFs)
- **統合時処理**: TF-TG pair毎にp-value最小のbest hitを1件保持 → **23,907レコード**

#### BS Source 3: RegPrecise

- **データ**: RegPrecise database (比較ゲノミクスベースのregulon予測)
- **レコード数**: **211レコード** (29 regulons / 31 unique TFs)
- **エビデンス**: computational_comparative

#### BS Source 4: Zorro-Aranda 2022 キュレーション済み相互作用

- **データ**: MOESM2 Table 1（文献キュレーション済みTF→target相互作用）
- **全レコード**: 9,714 interactions (101 TFs)
- **Strong evidenceのみ抽出**: **459レコード** (59 TFs)
- **実験手法**: EMSA, DNase I footprint, ChIP-chip, ChIP-seq, S1 mapping, IVT, RT-PCR等

### 4.2 統合方法

1. 各ソースのbinding site情報を共通スキーマに変換（TF_name, TF_SCO_ID, TF_gene_id, TG_SCO_ID, TG_gene_id, BS_contig, BS_start, BS_end, BS_strand, BS_sequence, BS_score, BS_pvalue, BS_qvalue, BS_source, BS_evidence）
2. 全レコードを結合して `master_TF_binding_sites_M145.tsv` を作成
3. gene_id（およびSCO_ID）ベースでTFマスターリストに集約（has_binding_site, BS_count, BS_sources, BS_evidence_types, motif_consensus）

### 4.3 結果サマリ

#### Binding Siteマスターテーブル

| 項目 | 数値 |
|------|------|
| **総レコード数** | **24,697** |
| ユニークTF数 | 97 |
| TG情報ありレコード | 24,218 |
| ユニークTF-TGペア | 24,218 |

#### ソース別レコード数

| Source | Records | Unique TFs | Evidence Type |
|--------|---------|-----------|---------------|
| ZorroAranda2022_MEME | 23,907 | 72 | computational (MEME) |
| ZorroAranda2022_Curated_Strong | 459 | 59 | experimental (various) |
| RegPrecise | 211 | 31 | computational (comparative genomics) |
| FIMO_curated_motif | 120 | 4 | computational (FIMO) |

#### TFレベル集約結果

| 項目 | 数値 |
|------|------|
| **BS情報ありTF数** | **82 / 1,705 (4.8%)** |
| モチーフコンセンサスありTF数 | 15 |

#### 信頼度レベル別BS保有率

| Confidence | BS保有数 / 全数 | 割合 |
|-----------|----------------|------|
| high | 56 / 155 | 36.1% |
| medium | 24 / 1,107 | 2.2% |
| low | 2 / 443 | 0.5% |

高信頼度TFの36%にbinding site情報が存在する一方、中・低信頼度では情報量が限定的であることは、well-studiedなTFほど文献や公開DBに情報が蓄積されていることと整合する。

#### モチーフコンセンサスが定義された15 TF

| SCO_ID | Gene名 | TFファミリー | モチーフ配列 | BS数 | ソース |
|--------|--------|------------|------------|------|--------|
| SCO0600 | (sigB) | Sigma | GNNTNNNNNNNNNNNNNNNNNGGGTAC | 512 | FIMO; Curated; MEME |
| SCO1489 | bldD | Unclassified_TF | AGTGANNNNNNTCACT | 143 | Curated; MEME |
| SCO1576 | (argR) | Unclassified_TF | TGAATAANNNNNNTTATTCA | 95 | RegPrecise; Curated; MEME |
| SCO2792 | adpA | Unclassified_TF | TGGCSNGWWY | 1,526 | Curated; MEME |
| SCO3063 | (draR) | TCS_response_regulator | AMAAWYMAKCA | 44 | Curated; MEME |
| SCO3571 | (crp) | Unclassified_TF | TGTGANNNNNNTCACA | 1,088 | Curated; MEME |
| SCO4159 | glnR | GntR | GTNACNNNNNNGANAC | 285 | Curated; MEME |
| SCO4596 | (abrC3) | TCS_response_regulator | GAASGSGRMS | 16 | Curated; MEME |
| SCO4907 | (afsQ1) | TCS_response_regulator | GTNACNNNNNNGTNAC | 28 | Curated; MEME |
| SCO5147 | (sigE) | Sigma | GGAACNNNNNNNNNNNNNNNNNGTT | 429 | FIMO; Curated; MEME |
| SCO5216 | sigR | Sigma | GGAACNNNNNNNNNNNNNNNNCGTT | 273 | FIMO; Curated; MEME |
| SCO5231 | dasR | GntR | ANTGGTCTAGACCANT | 375 | Curated; MEME |
| SCO5552 | ndgR | Unclassified_TF | BTYCANYWSYBGGAC | 27 | Curated; MEME |
| SCO6265 | scbR | TetR | TTTGGNNNNNNNNNNNNNNNNNNCCAAA | 1,508 | Curated; MEME |
| SCO7637 | (phoP) | TCS_response_regulator | GTTCACCNNNNGTTCACC | 1 | FIMO |

---

## 5. 出力ファイル

### 5.1 主要出力

| ファイル | 行数 | カラム数 | 内容 |
|---------|------|---------|------|
| `master_TF_list_M145.tsv` | 1,705 | 17 | TF候補マスターリスト |
| `master_TF_list_M145_with_BS.tsv` | 1,705 | 22 | ↑にBS情報5カラム追加 |
| `master_TF_binding_sites_M145.tsv` | 24,697 | 15 | Binding Site詳細テーブル |

### 5.2 カラム定義

#### master_TF_list_M145.tsv (17カラム)

| カラム | 説明 |
|-------|------|
| SCO_ID | 旧locus_tag (SCO####) |
| gene_id | 現行RefSeq gene_id (SC_RS#####) |
| gene_name | 遺伝子名 |
| contig | chromosome / SCP1 / SCP2 |
| start | 開始位置 |
| end | 終了位置 |
| strand | + / - |
| pfam_TF_flag | Pfam TFドメインの有無 (yes/no) |
| tf_related_domains | 検出Pfamドメイン (;区切り) |
| ZorroAranda_regulator_type | Zorro-Aranda 2022分類 |
| CastroMelchor_regulator_flag | Castro-Melchor 2010 cistronメンバー (1/0) |
| sigma_flag | σ因子フラグ (yes/no) |
| sigma_group | σ因子グループ |
| TF_family | TFファミリー分類 |
| source_flags | データソース (;区切り) |
| confidence | 信頼度 (high/medium/low) |
| product | NCBI RefSeq product名 |

#### 追加カラム (master_TF_list_M145_with_BS.tsv)

| カラム | 説明 |
|-------|------|
| has_binding_site | Binding site情報の有無 (yes/no) |
| BS_count | Binding siteレコード数 |
| BS_sources | データソース (;区切り) |
| BS_evidence_types | エビデンスタイプ (;区切り) |
| motif_consensus | キュレーション済みモチーフコンセンサス配列 |

#### master_TF_binding_sites_M145.tsv (15カラム)

| カラム | 説明 |
|-------|------|
| TF_name | TF名 |
| TF_SCO_ID | TFのSCO_ID |
| TF_gene_id | TFのgene_id |
| TG_SCO_ID | ターゲット遺伝子のSCO_ID |
| TG_gene_id | ターゲット遺伝子のgene_id |
| BS_contig | Binding siteのcontig |
| BS_start | 開始位置 |
| BS_end | 終了位置 |
| BS_strand | 鎖 (+/-) |
| BS_sequence | Binding site配列 |
| BS_score | スコア |
| BS_pvalue | p-value |
| BS_qvalue | q-value |
| BS_source | データソース |
| BS_evidence | エビデンスタイプ |

### 5.3 中間ファイル

| ファイル | 内容 |
|---------|------|
| `intermediate/M145_gene_basic_info.tsv` | 全8,083遺伝子の基本情報 |
| `intermediate/M145_pfam_domtblout.txt` | HMMER domtblout (14,823ヒット) |
| `intermediate/ZorroAranda2022_regulators_standardized.tsv` | Zorro-Aranda標準化レギュレーター |
| `intermediate/CastroMelchor2010_regulators_dedup.tsv` | Castro-Melchor重複除去レギュレーター |
| `intermediate/SigmaFactors_M145_from_lit.tsv` | σ因子リスト |
| `intermediate/M145_all_genes_with_TF_annotations.tsv` | 全遺伝子TFアノテーション完全テーブル |
| `intermediate/curated_TF_motifs_all_sources.tsv` | キュレーション済みモチーフ16件 |
| `intermediate/curated_TF_motifs.meme` | MEME形式モチーフファイル |
| `intermediate/fimo_results/` | FIMOゲノムワイドスキャン結果 |
| `intermediate/fimo_binding_sites_processed.tsv` | FIMO処理済み結果 |
| `intermediate/ZorroAranda2022_MEME_binding_sites.tsv` | Zorro-Aranda MEMEサイト |
| `intermediate/ZorroAranda2022_Inferred_BS_pairs.tsv` | Zorro-Aranda推定ペア |
| `intermediate/ZorroAranda2022_curated_interactions.tsv` | Zorro-Arandaキュレーション相互作用 |
| `intermediate/RegPrecise_binding_sites.tsv` | RegPreciseサイト |
| `intermediate/Literature_known_TF_motifs.tsv` | 文献既知モチーフ |

---

## 6. 処理パイプライン

全処理は `scripts/` 以下の9本のPythonスクリプトで再現可能である。

| # | スクリプト | 内容 | 入力 | 出力 |
|---|-----------|------|------|------|
| 01 | `01_extract_gene_info.py` | GFFから全遺伝子基本情報を抽出 | genomic.gff | M145_gene_basic_info.tsv |
| 02 | `02_parse_castro_melchor.py` | Castro-Melchor cistron→SCO_IDマッピング | AddFile3, AddFile4 | CastroMelchor2010_regulators_dedup.tsv |
| 03 | `03_create_sigma_factor_list.py` | σ因子リスト作成 | 文献 + GFF | SigmaFactors_M145_from_lit.tsv |
| 04 | `04_standardize_zorro_aranda.py` | Zorro-Aranda データ標準化 | MOESM2 | ZorroAranda2022_regulators_standardized.tsv |
| 05 | `05_parse_pfam_and_integrate.py` | Pfamパース + 全ソース統合 | domtblout + 各中間ファイル | master_TF_list_M145.tsv |
| 06 | `06_add_confidence_and_finalize.py` | 信頼度レベル追加 | master_TF_list_M145.tsv | master_TF_list_M145.tsv (更新) |
| 07 | `07_extract_zorro_binding_sites.py` | Zorro-Aranda BS情報抽出 | MOESM2, MOESM3 | MEME_BS + curated interactions TSV |
| 08 | `08_create_motif_table_and_fimo.py` | モチーフテーブル + FIMOスキャン | キュレーションモチーフ + ゲノム | curated_TF_motifs.meme + fimo_results/ |
| 09 | `09_integrate_binding_sites.py` | 全BS情報統合 | 全中間ファイル | master_TF_binding_sites_M145.tsv, master_TF_list_M145_with_BS.tsv |

---

## 7. 考察と留意点

### 7.1 TF候補数について

本解析で同定された1,705候補は、*S. coelicolor* ゲノムの約21%に相当する。この数は一見多いが、その主因はCastro-Melchor 2010のcistron（オペロン）単位でのレギュレーターフラグ付けにある。734遺伝子がCastro-Melchorソースのみで抽出されており、その多くはレギュレーターオペロン内の同居遺伝子（レギュレーター自体ではない）である可能性が高い。

実用上は:
- **高信頼度 (155)**: 確実なTF/σ因子。全解析に使用可
- **中信頼度 (1,107)**: Pfamドメインベースの候補。ドメイン構成を確認の上使用
- **高+中 (1,262)**: 標準的な解析セット
- **低信頼度 (443)**: フィルタリング推奨

### 7.2 Binding Site情報の網羅性

82/1,705 TF (4.8%) にbinding site情報が存在するが、高信頼度TFに限ると56/155 (36.1%) と高い割合を示す。これはwell-characterizedなTFほど文献に蓄積されている結果を反映している。

Binding site情報が無い1,623 TFのうち、多くはPfamドメインの存在からTFと推定されるが、直接的な結合サイト実験が未実施のものである。今後、ChIP-seq等の網羅的解析によりカバレッジ向上が期待される。

### 7.3 FIMO結果の解釈

FIMOゲノムワイドスキャンではp < 1e-4で58,460ヒットが得られたが、q < 0.01の厳しいフィルタリング後は120ヒット（4 TFのみ）に減少した。*S. coelicolor* の高GC含量（~72%）はGCリッチなモチーフの偽陽性を増大させるため、厳格な多重検定補正が必要である。

SigB (68), SigE (35), SigR (16) のECF/alternativeシグマ因子プロモーター配列が大半を占めるのは、これらのモチーフが長く（24-26 bp）特異的であるため、GC-richゲノムでも統計的に有意なヒットが得られることによる。一方、短い・縮退度の高いモチーフ（DraR: 11 bp, AbrC3: 10 bp等）はq < 0.01では検出されなかった。

### 7.4 SCO_IDに関する注意

文献ではSCO_ID（旧locus_tag）が広く使用されているが、本解析ではgene_id（SC_RS#####形式）を主キーとして使用している。SCO_IDとgene_idの対応は `intermediate/M145_gene_basic_info.tsv` で確認可能。ユーザーからの重要な修正として:

- **PhoP**: SCO7637 (旧データベースのSCO4230は誤り)
- **SigE**: SCO5147 (旧データベースのSCO3356は誤り)

---

## 8. 参考文献

1. Zorro-Aranda, A., Martínez-Antonio, A., & Freyre-González, J. A. (2022). Curation, inference, and assessment of a globally reconstructed gene regulatory network for *Streptomyces coelicolor* A3(2). *Scientific Reports*, 12, 2704.
2. Castro-Melchor, M., Charaniya, S., Karypis, G., Takano, E., & Hu, W.-S. (2010). Genome-wide inference of regulatory networks in *Streptomyces coelicolor*. *BMC Genomics*, 11, 578.
3. Sun, D., Liu, C., Zhu, J., & Liu, W. (2017). Connecting Metabolic Pathways: Sigma Factors in *Streptomyces* spp. *Frontiers in Microbiology*, 8, 2546.
4. Bentley, S. D., et al. (2002). Complete genome sequence of the model actinomycete *Streptomyces coelicolor* A3(2). *Nature*, 417, 141-147.
5. Tiffert, Y., et al. (2008). The *Streptomyces coelicolor* GlnR regulon: identification of new GlnR targets and evidence for a central role of GlnR in nitrogen metabolism. *Journal of Bacteriology*, 190(15), 5451-5461.
6. Allenby, N. E. E., et al. (2012). Genome wide transcriptional analysis of the phosphate starvation stimulon of *Bacillus subtilis*. *Journal of Bacteriology*, 194(15), 4088-4095. (PhoP motif reference)
7. Colson, S., et al. (2007). The pleiotropic regulator DasR links N-acetylglucosamine utilization, morphogenesis, and antibiotic production in *Streptomyces coelicolor*. *Molecular Microbiology*, 63(5), 1345-1360.
8. Schumacher, M. A., et al. (2017). The crystal structure of the BldD-c-di-GMP complex reveals a novel regulatory role of c-di-GMP in *Streptomyces* morphological development. *Nucleic Acids Research*, 45(2), 1127-1139.
9. Kim, M. S., et al. (2012). Conservation of thiol-oxidative stress responses regulated by SigR orthologues in actinomycetes. *Molecular Microbiology*, 85(2), 326-344.
10. Novichkov, P. S., et al. (2013). RegPrecise 3.0 — A resource for genome-scale exploration of transcriptional regulation in bacteria. *BMC Genomics*, 14, 745.

---

*最終更新: 2026-02-06*
