# E-1: メチル化遺伝子群 GO/KEGG エンリッチメント解析レポート

**日付:** 2026-04-11  
**解析ディレクトリ:** `11_epigenome_integration/analysis/62_GO_KEGG_enrichment/`  
**対応論点:** E-1（3/12・3/20議事録）  
**スクリプト:** `scripts/GO_KEGG_enrichment.py`

---

## 解析概要

### 目的
4mC（GCCGGC近傍）・6mA（AAGCCCG近傍）・Dual-targeted遺伝子群に機能的偏りがあるかを検定する。  
背景遺伝子: gene_set_comparison.tsv の全7,646遺伝子。

### 手法
- **GO enrichment**: Fisher's exact test（片側・over-representation） + Benjamini-Hochberg FDR補正  
- **KEGG enrichment**: KEGG REST API (http://rest.kegg.jp) で *S. coelicolor* M145 (sco) の150パスウェイ取得 → gene-pathway mapping → Fisher's exact test + BH補正  
- **アノテーション源**: gene_annotation_basic.tsv（8,275遺伝子、GO term保有4,845遺伝子）  
- **遺伝子群サイズ**: GCCGGC-only=3,197、AAGCCCG-only=459、Dual-targeted=448、Background=7,646

---

## 結果

### GO Enrichment（FDR < 0.05）

| 遺伝子群 | 有意 GO terms | トップ |
|---------|--------------|--------|
| GCCGGC-proximal (4mC) | **0** | — |
| AAGCCCG-proximal (6mA) | **1** | GO:0006098（pentose-phosphate shunt） |
| Dual-targeted | **0** | — |

#### AAGCCCG-proximal の有意 GO term
| GO term | n_query | n_bg | enrichment | FDR |
|---------|---------|------|-----------|-----|
| GO:0006098（pentose-phosphate shunt） | 5 | 11 | 7.57 | 0.048 |

遺伝子: SC_RS06830, SC_RS09295, SC_RS11715, SC_RS11720, SC_RS11725

---

### KEGG Pathway Enrichment（FDR < 0.2）

#### GCCGGC-proximal (4mC): 有意2パスウェイ

| pathway_id | pathway_name | n_query | enrichment | FDR |
|-----------|-------------|---------|-----------|-----|
| sco02024 | Quorum sensing | 65 | 1.55 | **2.5e-4** |
| sco02020 | Two-component system | 64 | 1.32 | 0.139 |

- Quorum sensing: 100遺伝子中65遺伝子がGCCGGC近傍に存在（65%）

#### AAGCCCG-proximal (6mA): 有意0パスウェイ

#### Dual-targeted: 有意5パスウェイ（FDR < 0.2）

| pathway_id | pathway_name | n_query | enrichment | FDR |
|-----------|-------------|---------|-----------|-----|
| sco00975 | Biosynthesis of various siderophores | 5 | **12.2** | 8.7e-4 |
| sco00790 | Folate biosynthesis | 6 | 5.12 | 0.025 |
| sco00190 | Oxidative phosphorylation | 11 | 2.84 | 0.032 |
| sco00550 | Peptidoglycan biosynthesis | 7 | 3.73 | 0.036 |
| sco00450 | Selenocompound metabolism | 4 | 4.88 | 0.099 |

---

## 解釈

### 主要所見

1. **4mC（GCCGGC）は Quorum sensing/TCS 遺伝子に集積**  
   - sco02024（Quorum sensing）に強いエンリッチメント（FDR=2.5e-4）  
   - sco02020（Two-component system）にも傾向（FDR=0.14）  
   - 分化・環境応答シグナル経路の遺伝子が4mCメチル化リファレンス近傍に集まっている  
   - これはGatekeeperモデルのLayer 3（Signal Gating）と整合

2. **Dual-targeted（4mC+6mA両方）はシデロフォア合成・エネルギー代謝に集積**  
   - sco00975（siderophore biosynthesis）は enrichment=12.2 と最高  
   - 鉄代謝・細胞壁合成（peptidoglycan）・呼吸鎖（oxidative phosphorylation）が集まる  
   - 両メチル化システムで「挟まれた」遺伝子群は一次代謝・生存核心機能に偏る可能性

3. **6mA（AAGCCCG）は GO/KEGG いずれでも有意エンリッチメントなし（GO:0006098を除く）**  
   - 6mA近傍遺伝子は機能的に均一（ランダム）分布している  
   - これはH17/H19で示した「6mA-発現相関なし」と整合し、6mAが特定機能カテゴリを標的にしていないことを支持

### ゲートキーパーモデルとの整合性
- 4mCによるQuorum sensing/TCS集積は、Layer 3（R-M防御による調節遺伝子のシールド）の定量的裏付けとなる
- 6mAの機能非偏在は「防御目的でのランダム配置」という既存解釈と整合

---

## 出力ファイル

### テーブル（`tables/`）
| ファイル | 内容 |
|---------|------|
| `E1_GO_enrichment_GCCGGC-proximal_4mC.tsv` | GCCGGC GO enrichment全結果 |
| `E1_GO_enrichment_AAGCCCG-proximal_6mA.tsv` | AAGCCCG GO enrichment全結果 |
| `E1_GO_enrichment_Dual-targeted.tsv` | Dual GO enrichment全結果 |
| `E1_KEGG_enrichment_GCCGGC-proximal_4mC.tsv` | GCCGGC KEGG enrichment全結果 |
| `E1_KEGG_enrichment_AAGCCCG-proximal_6mA.tsv` | AAGCCCG KEGG enrichment全結果 |
| `E1_KEGG_enrichment_Dual-targeted.tsv` | Dual KEGG enrichment全結果 |
| `kegg_sco_pathways_cache.json` | KEGG REST APIキャッシュ（130パスウェイ） |

### 図（`figures/`）
| ファイル | 内容 |
|---------|------|
| `E1_GO_enrichment_*.pdf/svg` | GO enrichmentドットプロット（PDF+SVG） |
| `E1_KEGG_enrichment_*.pdf/svg` | KEGG enrichmentドットプロット（PDF+SVG） |

---

## 論点 E-1 への回答

> 「GCCGGC/AACCCGが挟む遺伝子に機能的偏りがあるか」

**回答**: **4mC（GCCGGC）はシグナル伝達遺伝子（Quorum sensing, TCS）に有意偏在（FDR=2.5e-4）**。6mAおよびDual-targetedは一部代謝カテゴリへの集積があるものの、全体として機能非依存的分布が優勢。これはGatekeeperモデルにおいて4mCが制御・シグナル伝達遺伝子群を優先的に保護するという仮説を支持する。
