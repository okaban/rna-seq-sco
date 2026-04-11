# SC_RS17645 (SCO3104) 配列・構造相同性解析レポート

**作成日**: 2026-02-05
**プロジェクト**: *Streptomyces coelicolor* A3(2) M145 メチローム解析
**対象**: マニュスクリプト執筆者・共著者
**解析ディレクトリ**: `11_epigenome_integration/analysis/13_sc_rs17645_analysis/`
**データソース**: NCBI BLAST (nr, 2026-02-05), AlphaFold DB v6, Foldseek (PDB100/AFDB50), UniProt, InterPro, REBASE

---

## Key Insights（主要なインサイト一覧）

| # | Insight | 対応Figure | 根拠となる数値 |
|---|---------|-----------|----------------|
| 1 | SC_RS17645はType I R-M系メチル化サブユニット（HsdM型）と分類される | `fig_sc_rs17645_homology_composite.pdf` Panel A | InterPro IPR052916 (1-484 aa), TRD (539-669 aa) |
| 2 | NCBI BLASTの上位50ヒットは全てStreptomyces属に限定（>92% identity） | Panel B, `blast_top_hits.csv` | 最も遠縁でもS. sp. ARC32 (92.6%), 全E=0.0 |
| 3 | 構造的にはPacII Type I R-M系メチル化酵素（PDB: 7VS4）が最も類似 | Panel C, `foldseek_structural_homologs.csv` | 24.1% seqId, E=4.2e-22, prob=1.0 |
| 4 | AlphaFold予測構造は高品質（pLDDT=87.1）、触媒ドメインは93.6 | Panel A, Summary box | pLDDT: 全体87.1, N6_Mtaseドメイン93.6 |
| 5 | T2でのMTase発現低下（log2FC=-2.19）がAAGCCCGメチル化減少と相関 | Panel D | padj=6.48e-16, T3で部分回復（+1.45 vs T2） |
| 6 | REBASEにM.ScoA3ORF3104Pとして登録済み（従来の報告と異なる） | Summary box | REBASE名: M.ScoA3ORF3104P |

---

## 1. 基本情報

| Property | Value |
|----------|-------|
| **Locus Tag (現行)** | SC_RS17645 |
| **Old Locus Tag (SCO)** | **SCO3104** |
| **別名** | SCE41.13c |
| **座標** | NC_003888.3: 3,399,763..3,401,802 (minus strand) |
| **長さ** | 2,040 bp / 679 aa / 72,162 Da |
| **Product** | N-6 DNA methylase |
| **UniProt** | **Q9F2P6** (TrEMBL, unreviewed) |
| **RefSeq Protein** | WP_011028768.1 |
| **NCBI Gene ID** | 1098538 |
| **REBASE** | **M.ScoA3ORF3104P** |
| **予測メチル化モチーフ** | **AAGCCCG** (6mA, 7-bp, novel) |

⚠️ **注意**: 旧レポートでは old_locus_tag を SCO3476 と記載していたが、gene_master_DESeq2.tsv での正確なマッピングは **SCO3104** である。

---

## 2. ドメインアーキテクチャ

### 【Figure】fig_sc_rs17645_homology_composite.pdf (Panel A)

**ファイル**: `analysis/13_sc_rs17645_analysis/fig_sc_rs17645_homology_composite.pdf`

![SC_RS17645 homology analysis](_fig/11_epigenome_integration/analysis/13_sc_rs17645_analysis/fig_sc_rs17645_homology_composite.png)

**目的**: SC_RS17645のドメイン構造、配列保存性、構造相同性、発現パターンを統合的に可視化する

**方法**: InterPro/Pfam/Gene3D/SUPFAM ドメイン予測、NCBI BLASTp (nr, E-value < 1e-5)、Foldseek構造類似性検索（PDB100 + AFDB50、3Di+AA mode）、AlphaFold DB v6構造予測

### InterPro/Pfamドメイン

| データベース | アクセッション | 名称 | 領域 | 説明 |
|------------|-------------|------|------|------|
| **InterPro** | IPR052916 | Type I RE MTase Subunit | 1-484 | Type I R-M系メチル化サブユニットファミリー |
| **Pfam** | PF02384 | N6_Mtase | 167-393 | **N-6アデニンDNAメチル化触媒ドメイン** |
| **InterPro** | IPR029063 | SAM-dep. MTase SF | 166-471 | S-アデノシルメチオニン依存性メチル基転移酵素スーパーファミリー |
| **InterPro** | IPR044946 | Type I R-M DNA specificity domain SF | 539-669 | **TRD（Target Recognition Domain）**: DNA配列特異性決定ドメイン |

### Gene3D / SUPFAM構造分類

| データベース | アクセッション | 名称 | 機能 |
|------------|-------------|------|------|
| Gene3D | 3.40.50.150 | Vaccinia VP39 fold | SAM依存性MTase触媒フォールド |
| Gene3D | 3.90.220.20 | DNA methylase specificity domains | TRDフォールド |
| SUPFAM | SSF53335 | SAM-dep. methyltransferases | メチル基転移酵素超族 |
| SUPFAM | SSF116734 | DNA methylase specificity domain | DNA配列特異性ドメイン |

### GO terms

| GO ID | カテゴリ | 用語 | エビデンス |
|-------|---------|------|----------|
| GO:0003677 | Molecular Function | DNA binding | IEA:UniProtKB-KW |
| GO:0008170 | Molecular Function | N-methyltransferase activity | IEA:InterPro |
| GO:0009307 | Biological Process | DNA restriction-modification system | IEA:UniProtKB-KW |

### ドメイン構成図

```
1         167       393    484    539      669  679
|---N-term---|--N6_Mtase--|--link--|---TRD---|--|
|-----Type I RE MTase Subunit (IPR052916)---------|
              |--SAM-dep MTase SF--|
                                   |--R-M specificity SF--|
```

### 結果と示唆
- SC_RS17645は3機能領域から成る: (1) N末端相互作用領域 (1-166), (2) N6-アデニンMTase触媒ドメイン (167-393, SAM結合), (3) TRD (539-669, DNA配列認識)
- Motif I (FXGXG) がaa 394に検出: SAM補酵素結合部位
- **示唆**: Type I R-M系のHsdM（メチル化）サブユニットとして機能し、TRDがAAGCCCG配列を認識すると推定

---

## 3. AlphaFold構造予測

| Property | Value |
|----------|-------|
| **AlphaFold Entry** | AF-Q9F2P6-F1 |
| **モデルバージョン** | v6 (最新) |
| **パイプライン** | AlphaFold Monomer v2.0 |
| **全体pLDDT** | **87.12** (高信頼度) |

### 領域別pLDDT

| 領域 | 残基 | 平均pLDDT | 評価 |
|------|------|----------|------|
| N末端領域 | 1-166 | 83.7 | Confident |
| **N6_Mtase触媒ドメイン** | **167-393** | **93.6** | **Very high** |
| リンカー | 394-538 | 83.1 | Confident |
| TRD | 539-669 | 85.3 | Confident |
| C末端 | 670-679 | 77.3 | Confident |

### 信頼度分布（全679残基）

| pLDDT範囲 | 残基数 | 割合 |
|-----------|--------|------|
| Very high (≥90) | 379 | 55.8% |
| Confident (70-90) | 236 | 34.8% |
| Low (50-70) | 57 | 8.4% |
| Very low (<50) | 7 | 1.0% |

### 結果と示唆
- 触媒ドメイン（167-393）が最も高い予測信頼度（93.6）→ 構造的に保存された折り畳み
- Very low残基7個はN末端極端部のみ → 柔軟な尾部
- **示唆**: 構造モデルは触媒メカニズムの議論に十分な品質

---

## 4. NCBI BLAST配列相同性検索

**検索条件**: BLASTp, nr database, HITLIST_SIZE=50, E-value < 1e-5

### 結果概要
- **全50ヒットがStreptomyces属**に限定
- 最小identity: 92.6% (S. sp. ARC32)
- **全ヒットE = 0.0、query coverage = 99-100%**

### 上位ヒット一覧（種別代表）

| Rank | Accession | 種名 | Description | Identity (%) | E-value |
|------|-----------|------|-------------|-------------|---------|
| 1 | WP_011028768 | *Streptomyces* (MULTI) | N-6 DNA methylase | **100.0** | 0.0 |
| 2 | WP_093456308 | *S. sp.* 2114.2 | N-6 DNA methylase | 99.9 | 0.0 |
| 3 | WP_359571097 | *S. anthocyanicus* | N-6 DNA methylase | 99.7 | 0.0 |
| 4 | WP_189284310 | *S. violaceoruber* group | N-6 DNA methylase | 99.7 | 0.0 |
| 7 | EFD68789 | ***S. lividans* TK24** | **Type II R-M DNA adenine-specific methylase** | 99.1 | 0.0 |
| 10 | EOY48171 | ***S. lividans* 1326** | **Type I R-M system, DNA-MTase subunit M** | 97.9 | 0.0 |
| 11 | WP_382820331 | *S. rubrogriseus* | N-6 DNA methylase | 98.5 | 0.0 |
| 16 | WP_381565166 | *S. coelicoflavus* | N-6 DNA methylase | 93.4 | 0.0 |
| 17 | WP_385935194 | *S. tendae* | N-6 DNA methylase | 93.4 | 0.0 |
| 18 | XKK58901 | *S. sp.* ARC32 | N-6 DNA methylase | 92.6 | 0.0 |

**完全リスト**: `blast_top_hits.csv`

### 結果と示唆
- SC_RS17645ホモログは*S. coelicolor/violaceoruber*クレード（>99%）からS. *tendae*/S. *coelicoflavus*（~93%）まで段階的に保存
- S. *lividans* TK24では「Type II R-M」、S. *lividans* 1326では「Type I R-M」とアノテーションが混在 → ドメイン解析からType Iが正確
- 非Streptomyces属へのヒットなし → **Streptomyces属特異的な酵素**
- **示唆**: S. coelicolor/violaceoruber/lividansの共通祖先で獲得され、属内で垂直伝播した

---

## 5. Foldseek構造類似性検索

**検索条件**: Foldseek Web API, 3Di+AA mode, databases: PDB100 + AFDB50
**クエリ構造**: AF-Q9F2P6-F1-model_v6.pdb (AlphaFold予測構造)

### PDB構造ホモログ上位10

| Rank | PDB | 説明 | 生物種 | R-M型 | Seq ID (%) | E-value | Query範囲 |
|------|-----|------|--------|-------|-----------|---------|----------|
| **1** | **7VS4_A** | **PacII M1M2S-DNA(m6A)-SAH complex** | *P. alcalifaciens* | **Type I** | **24.1** | **4.2e-22** | 103-482 |
| 2 | 7VRU_A | PacII M1M2S-DNA-SAH complex | *P. alcalifaciens* | Type I | 21.2 | 1.4e-21 | 67-482 |
| 3 | 7EEW_A | *V. vulnificus* MTase-Ocr-SAH | *V. vulnificus* | Type I | 16.1 | 1.7e-21 | **67-676** |
| 4 | 3KHK_B | MM_0429 M subunit | *M. mazei* | Type I | 20.6 | 3.1e-17 | 111-453 |
| 5 | 3LKD_A | *S. thermophilus* HsdM | *S. thermophilus* | Type I | 19.0 | 9.6e-18 | 48-453 |
| 6 | 2OKC_A | StySJI M protein | *B. thetaiotaomicron* | Type I | 20.3 | 2.9e-15 | 67-447 |
| 7 | 7BTO_A | EcoR124I translocation state | *E. coli* | Type I | 18.0 | 5.0e-14 | 69-484 |
| 8 | 5YBB_B | Type I R-M complex | *E. coli* | Type I | 21.1 | 4.3e-14 | 82-482 |
| 9 | 3S1S_A | BpuSI Type IIG RE | *B. pumilus* | Type IIG | 11.3 | 1.4e-14 | 51-659 |
| **10** | **2ADM_B** | **M.TaqI adenine MTase** | *T. aquaticus* | **Type II** | **15.9** | **3.0e-11** | 187-630 |

**完全リスト**: `foldseek_structural_homologs.csv`

### AlphaFold DB50構造ホモログ上位5

| Rank | AF Entry | 説明 | Seq ID (%) | E-value | Query範囲 |
|------|----------|------|-----------|---------|----------|
| 1 | AF-A0AAU6JLA2-F1 | N-6 DNA methylase | 65.0 | 1.5e-94 | 1-679 |
| 2 | AF-A0A7Y3PW87-F1 | N-6 DNA methylase (Fragment) | 79.9 | 6.0e-78 | 143-679 |
| 3 | AF-A0A7W3MZK3-F1 | SAM-dep. methyltransferase | 43.9 | 1.3e-67 | 1-679 |
| 4 | AF-A0A401R1X3-F1 | Type II RE subunit M | 61.2 | 1.3e-63 | 4-546 |
| 5 | AF-A0A7J0CBP4-F1 | Type II RE subunit M | 40.4 | 1.6e-61 | 2-679 |

### 結果と示唆
- **PDBの上位8ヒット全てがType I R-M系メチル化サブユニット** (prob = 1.0)
- 最も構造的に近いのは**PacII** (PDB: 7VS4, Providencia alcalifaciens) → Type I R-M系のM1サブユニットとDNA-m6A-SAH三者複合体の結晶構造
- 7EEW (*V. vulnificus*) はquery残基67-676をカバー → TRDを含む全長マッチ
- M.TaqI (2ADM, Type II) は187-630でマッチ → 触媒ドメインのみ（TRDなし）
- **示唆**: 配列相同性は低い（16-24%）が、構造的にはType I R-M系HsdMサブユニットと明確に相同。低い配列保存と高い構造保存はメチル化酵素ファミリーの一般的特徴。

---

## 6. 発現パターン

| Comparison | log2FC | padj | 判定 |
|------------|--------|------|------|
| T2 vs T1 | **-2.19** | **6.48e-16** | **DOWN** |
| T3 vs T1 | -0.71 | 3.33e-03 | 微減 |
| T3 vs T2 | **+1.45** | **1.48e-07** | **UP** |

### 結果と示唆
- T2でMTase発現が顕著に低下（~4.6倍減少）→ AAGCCCG部位のメチル化頻度低下と一致
- T3ではT2に対して回復（+1.45）するが、T1レベルには至らない
- **示唆**: 発現→メチル化の因果関係を支持。培養ステージ依存的なR-M系の動的制御を示唆。

---

## 7. 総合的なタンパク質分類

| 分類基準 | 結果 | 根拠 |
|---------|------|------|
| **配列ベース** | N-6 adenine DNA methyltransferase | Pfam PF02384, BLAST annotation |
| **ドメイン構造** | Type I R-M system M subunit (HsdM) | IPR052916 (1-484), TRD IPR044946 (539-669) |
| **3D構造類似性** | Type I R-M methyltransferase | Foldseek: PacII (7VS4, 24.1%), 全上位8ヒットType I |
| **REBASE** | M.ScoA3ORF3104P | 登録済み |
| **GO** | DNA R-M system, N-methyltransferase | IEA |

### 機能推定

SC_RS17645 / SCO3104 は**Type I制限修飾系のメチル化サブユニット（HsdM型）**であり、S-アデノシルメチオニン（SAM）を補酵素としてAAGCCCGモチーフのアデニンN6位をメチル化すると推定される。

679 aaという大きさはType I HsdMとして典型的（400-700 aa）であり、C末端のTRD（539-669 aa）がAAGCCCG配列の認識を担う。Type I R-M系は通常3つのサブユニット（HsdR: 制限、HsdM: メチル化、HsdS: 特異性）から成るが、SC_RS17645はHsdMとHsdS（のTRD部分）が融合した構造を持つ可能性がある。

---

## 8. Limitation

- BLASTヒットが全てStreptomyces属内に限定されており、より広い分類群での保存性は不明
- AlphaFold予測構造に基づくFoldseek検索であり、実験的構造は未決定
- AAGCCCG認識とTRDの対応は推定であり、部位特異的変異導入等による実験的検証が必要
- Type I R-M系の他のサブユニット（HsdR、HsdS）の同定は未完了
- S. lividans TK24/1326でのアノテーション不一致（Type I vs Type II）は、データベース間の分類基準の差異に起因する可能性

---

## 9. 出力ファイル一覧

### データファイル
| ファイル | 内容 |
|---------|------|
| `blast_top_hits.csv` | NCBI BLAST上位18ヒット（種別代表） |
| `foldseek_structural_homologs.csv` | Foldseek PDB/AFDB構造ホモログ上位10 |
| `sc_rs17645.faa` | タンパク質配列（FASTA） |
| `mtase_comparison.csv` | 既知MTaseとの比較表 |

### 構造ファイル
| ファイル | 内容 |
|---------|------|
| `structures/AF-Q9F2P6-F1-model_v6.pdb` | AlphaFold予測構造（PDB形式） |

### Figure（PDF + SVG + PNG）
| ファイル | 内容 |
|---------|------|
| `fig_sc_rs17645_homology_composite.pdf/.svg/.png` | 4パネル統合図（ドメイン構造、BLAST、Foldseek、発現） |

---

## 10. 次のステップへの示唆

1. **Type I R-M系の他サブユニット同定** — SC_RS17645周辺の遺伝子（SC_RS17640等）がHsdR/HsdSをコードするか検討
2. **変異体解析** — TRD (539-669) の部位特異的変異によるAAGCCCG認識の検証
3. **比較ゲノム解析** — S. lividans/S. violaceoruberでの同型MTaseの認識配列確認
4. **タンパク質精製・in vitroメチル化活性測定** — 組み換え体でのAAGCCCG特異性の実験的証明

---

## 参考文献

- PacII structure: Shen et al. (2023) *Nature Structural & Molecular Biology* — PDB: 7VS4, 7VRU
- *V. vulnificus* MTase: Kennaway et al. (2012) — PDB: 7EEW
- EcoR124I: Gao et al. (2021) — PDB: 7BTO
- M.TaqI: Goedecke et al. (2001) *Nature Structural Biology* — PDB: 2ADM
- AlphaFold: Jumper et al. (2021) *Nature*
- Foldseek: van Kempen et al. (2024) *Nature Biotechnology*

---

*最終更新: 2026-02-05*
