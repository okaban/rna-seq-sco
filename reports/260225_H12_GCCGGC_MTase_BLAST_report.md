# H12: M.Svi27968I BLAST相同性検索によるGCCGGC N4-C MTase同定

**日付**: 2026-02-25
**解析ディレクトリ**: `11_epigenome_integration/analysis/35_GCCGGC_MTase_BLAST/`
**判定**: **H10候補覆し — SC_RS19770が最有力（E=0.007）、SC_RS24685は相同性なし**

---

## 1. 背景と仮説

H10でGCCGGC認識N4-C MTaseの候補をアノテーション+発現パターンで絞り込み、SC_RS24685（score 9/10）を最有力候補とした。しかしH10は配列相同性に基づく検証が未実施であった。

**H12仮説**: M.Svi27968I（*S. violascens* ATCC 27968、REBASE唯一の*Streptomyces* GCCGGC Nm4C産生酵素、GenBank: QPA00148.1）をM145プロテオームにBLAST検索すると、SC_RS24685が最高スコアでヒットする。

## 2. 方法

1. NCBI datasets CLIでCP029377（*S. violascens* ATCC 27968）プロテオームをダウンロード
2. QPA00148.1（M.Svi27968I, 429 aa）をクエリ配列として抽出
3. M145プロテオーム（7,872配列）でBLASTデータベース構築（makeblastdb）
4. blastp 2.17.0+で検索

## 3. 結果

### 3.1 BLAST Top Hits

| Rank | Protein ID | Locus Tag | Product | E-value | Identity | Similarity | Coverage |
|------|-----------|-----------|---------|---------|----------|------------|----------|
| **1** | **WP_254613693.1** | **SC_RS19770** | DNA cytosine methyltransferase | **0.007** | **29%** | **41%** | 109/429 aa |
| **2** | **WP_011031229.1** | **SC_RS36410** | DNA cytosine methyltransferase | **0.047** | **26%** | **42%** | 118/429 aa |
| 3 | WP_011028829.1 | SC_RS18680 | class II DAHP synthase | 0.076 | 29% | 42% | (非MTase) |
| 4 | WP_011030310.1 | SC_RS30175 | 16S rRNA methyltransferase RsmD | 0.97 | 43% | 55% | (RNA MTase) |
| 5 | WP_161270170.1 | SC_RS25950 | class I SAM-dependent MTase | 1.5 | 47% | 60% | 30/429 aa |

### 3.2 H10候補との比較

| H10 Rank | Locus Tag | H10 Score | BLAST E-value | BLAST Rank |
|----------|-----------|-----------|---------------|------------|
| **1** | **SC_RS24685** | **9/10** | **ヒットなし (E>10)** | **圏外** |
| 2 | SC_RS28835 | 9/10 | ヒットなし | 圏外 |
| 74 (最下位) | **SC_RS19770** | **3/10** | **0.007** | **1位** |
| 75 (最下位) | **SC_RS36410** | **3/10** | **0.047** | **2位** |

### 3.3 SC_RS19770 / SC_RS36410の発現プロファイル

| Locus Tag | T1 (norm count) | T2 | T3 | LFC T2vsT1 | LFC T3vsT1 | 備考 |
|-----------|-----------------|----|----|------------|------------|------|
| SC_RS19770 | 11.0 | 4.6 | 62.1 | -1.27 | +2.49 | Pseudogene近傍 |
| SC_RS36410 | 1.3 | 4.2 | 64.9 | +1.66 | +5.59 | **Defense island** (H7) |

### 3.4 発現パラドックスの解釈

**問題**: T1で1,516 GCCGGC 4mCサイトが存在するのに、両候補MTaseのT1発現はきわめて低い（11, 1.3）。

**解釈（H7との統合）**:
- H7で判明した事実: CCGGサイトはタイムポイント間で**完全非重複（Jaccard=0.000）**、T1=core (77%), T2=arms (89%)
- **T1 core-enriched GCCGGC sites (1,516個)の責任酵素は依然不明** — SC_RS19770/SC_RS36410では説明不可
- **SC_RS36410のT3急誘導（LFC=+5.59）**は、defense island co-inductionと一致（H7）
- SC_RS36410は**T3 arm-enriched sites (30個)のみ**を生成する可能性
- T1 core sites → **未知の第3酵素**（M145プロテオーム中にM.Svi27968I相同性なし）
  - 可能性: (a) 非常に遠い相同性でBLAST検出不能、(b) 収束進化で同じモチーフを認識する非相同酵素、(c) T1培養前に既にメチル化されており維持のみ必要

### 3.5 重要な発見: アノテーションの誤解

M.Svi27968I自体のNCBIアノテーションが「DNA (cytosine-**5**-)-methyltransferase」であるにもかかわらず、REBASEではN**4**-methylcytosine産生と確認されている。これは、SC_RS19770/SC_RS36410の「DNA cytosine methyltransferase（Dcm-like）」アノテーションが5mC特異的とは限らないことを示す。

## 4. 結論

| 項目 | 結果 |
|------|------|
| H12仮説（SC_RS24685が最高ヒット） | **棄却** |
| M.Svi27968I最近縁M145タンパク | **SC_RS19770** (E=0.007, 29% identity) |
| 2番目 | **SC_RS36410** (E=0.047, 26% identity) — defense island |
| H10 top候補SC_RS24685 | BLAST圏外（M.Svi27968Iとの相同性なし） |
| T1 core GCCGGC sites (1,516) 責任酵素 | **未同定のまま**（既知N4-C MTaseとの相同性でM145内に該当なし） |
| SC_RS36410の役割 | T3 defense island活性化に伴うarm-enriched sites生成（少数） |

## 5. 出力ファイル

### Tables
| ファイル | 内容 |
|---------|------|
| `blast_results_Svi27968I_full.txt` | blastp全結果（5ヒット） |

### Data
| ファイル | 内容 |
|---------|------|
| `M_Svi27968I.fasta` | M.Svi27968Iタンパク配列 (QPA00148.1, 429 aa) |
| `M_PfrJS2V.fasta` | M.PfrJS2Vタンパク配列 |
| `M145_protein.faa` | M145全プロテオーム (7,872 seq) |
| `m145_prot.*` | M145 BLASTデータベース |
| `sviolascens_*` | S. violascensプロテオーム・GFF・BLASTデータベース |

---

*Analysis directory: `11_epigenome_integration/analysis/35_GCCGGC_MTase_BLAST/`*
*Generated: 2026-02-25*
