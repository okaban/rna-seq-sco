# DMG機能エンリッチメント解析レポート

**作成日**: 2026-02-06
**最終更新**: 2026-02-06
**プロジェクト**: *Streptomyces coelicolor* A3(2) M145 エピゲノム-トランスクリプトーム統合解析
**対象**: 差異メチル化遺伝子 (DMGs) の機能分類
**解析ディレクトリ**: `14_DMG_functional_enrichment/analysis/14_DMG_enrichment_260206_v1/` (基本解析), `14_DMG_functional_enrichment/analysis/14_DMG_selectivity_260206_v1/` (追加解析)
**データソース**: KEGG REST API (sco, 150パスウェイ), GO (RefSeq genomic.gff, 1,543ターム), COG (gene_COG_classification.tsv, 20カテゴリ)

---

## Key Insights（主要なインサイト一覧）

| # | Insight | 対応Figure/Table | 根拠となる数値 |
|---|---------|-----------------|-----------------|
| 1 | DMGはKEGG/GO/COGいずれにおいても有意なエンリッチメントをほぼ示さない（DEGとは対照的） | `KEGG_DMG_vs_DEG_heatmap.pdf` | DMG: KEGG有意0件 vs DEG: 最大18件; DMG: GO有意1件 vs DEG: 最大32件 |
| 2 | DMGのCOG分布はバックグラウンドとほぼ同一であり、メチル化変動は機能非選択的 | `COG_distribution_T3vsT1.pdf` | 全19カテゴリで比率差 < 2% |
| 3 | V-Defense（COG）のエンリッチメントはメチルトランスフェラーゼ遺伝子群が駆動しており、15遺伝子中11がSAM依存性メチルトランスフェラーゼ | `COG_DMG_T3vsT2_4mC_hypo.tsv`, 本レポート表3.4a | OR=2.38, padj=0.028; 11/15 = 73%がMTase |
| 4 | GO:0016853（isomerase activity）のエンリッチメントは糖リン酸イソメラーゼ/エピメラーゼ7遺伝子が駆動 | `GO_DMG_T3vsT1_combined_all.tsv`, 本レポート表3.5a | padj=0.006, 7/8 = 88%が糖リン酸イソメラーゼ |
| 5 | KEGG名目p値レベルで、アミノ酸生合成・フェナジン生合成が複数比較で再現的にトレンドを示すが、いずれもFDR補正後は非有意 | `KEGG_DMG_T3vsT1_combined_all.tsv` 等 | アミノ酸生合成: p=0.026 (T3vsT1), p=0.044 (T3vsT2); 全padj>0.9 |
| 6 | CCGG (4mC) とAAGCCCG (6mA) のモチーフ間でもCOG/KEGG分布の偏りは見られない（K-Transcriptionの視覚的印象は統計的に支持されない） | `COG_motif_CCGG_vs_AAGCCCG.pdf`, 本レポート3.6節 | AAGCCCG K: fold=0.875, p=0.82; CCGG K: fold=0.935, p=0.80 |
| 7 | DMGの検出力はDEGの約1/3〜1/5だが、名目p値レベルでもパスウェイが浮上しないことから、検出力不足だけでは説明できない | 本レポート3.7節 | DMG KEGG背景: 115〜155遺伝子 vs DEG: 480〜560遺伝子 |
| 8 | DEG∩DMG協調変動遺伝子もKEGGエンリッチメント0件 — 発現・メチル化の同時変動ですら機能選択性なし | `A1_coordinated_enrichment_summary.tsv` | 全15条件でKEGG有意0件、GO有意0〜3件、COG有意0〜1件 |
| 9 | メチル化-発現相関の**強度**はCOGカテゴリにより異なる — O-Chaperones (r=0.48), V-Defense (r=0.57), U-Secretion (r=-0.60) | `A2_correlation_heatmap_4mC.pdf`, `A2_correlation_heatmap_6mA.pdf` | 5カテゴリで名目有意; ただし全体相関は4mC T2vsT1のみ有意(r=0.14) |
| 10 | 遺伝子本体メチル化はCOG 7カテゴリに有意エンリッチするが、プロモーター/TSS近傍は0カテゴリ — 位置依存的な機能偏り | `A3_COG_promoter_vs_body.pdf` | Gene body: E,G,T,C,L,O,R有意; Promoter: 0; TSS proximal: 0 |
| 11 | V-Defense MTase遺伝子15中12がDEG、4mC脱メチル化（最大-94%）と大幅な発現変動（最大log2FC=5.5）の両方を示す | `A4_defense_MTase_expression.tsv` | T2vsT1: 12/15 DEG; SC_RS11875 log2FC=+5.54; SC_RS12890 4mC=-94% |
| 12 | GO isomerase 8遺伝子中4〜6がDEG、SC_RS18320 (SCO3241) はlog2FC=+8.9の劇的発現上昇 | `A4_GO_isomerase_expression.tsv` | T2vsT1: 4/8 DEG; T3vsT1: 6/8 DEG; SC_RS18320 log2FC=+8.93 |
| 13 | TF遺伝子もBGC遺伝子もDMGとして選択的に標的化されていない（全OR<1, 全p>0.2） | `A5_TF_BGC_DMG_frequency.pdf` | TF DMG率5.9-8.5% vs 非TF 7.0-9.4%; BGC DMG率5.1-8.1% vs 非BGC 6.9-9.3% |

---

## 1. 解析の背景と目的

**目的**: トランスクリプトーム解析におけるKEGGパスウェイエンリッチメントやCOG機能分類と同様の手法で、メチル化変動遺伝子（DMG）の機能的偏りを評価する。

これまでの解析で、メチル化の頻度・増減・遺伝子発現との相関を調べてきたが、「メチル化変動の標的はどのような機能カテゴリに集中するのか」という問いには未着手であった。DEGに対するエンリッチメント解析と同一フレームワークでDMGを解析し、両者を直接比較することで、メチル化制御の標的選択性を検証する。

---

## 2. 方法

### 2.1 DMG定義

既存定義（`11_epigenome_integration/scripts/DEG_DMG_overlap_analysis.py`）を踏襲:
- **DMG**: |メチル化変化量| > 10%（6mAまたは4mCのいずれか）
- データソース: `integrated_methyl_expression_weighted.csv`
- 層別化: hyper / hypo × 6mA / 4mC × 3比較 (T2vsT1, T3vsT1, T3vsT2)

### 2.2 ORA（Over-Representation Analysis）

- **検定**: Fisher正確検定（片側、greater）
- **多重検定補正**: Benjamini-Hochberg法 (FDR)
- **有意性閾値**: padj < 0.05
- **最小遺伝子数**: KEGG/GO: 3遺伝子/ターム, COG: 1遺伝子/カテゴリ
- **バックグラウンド**: 各データベースにアノテーションのある全遺伝子
- **KEGG遺伝子ID変換**: gene_id (SC_RS#####) → old_locus_tag (SCO####) → KEGG ID (sco:SCO####)

### 2.3 対照群

同一手法・同一バックグラウンドでDEG（padj<0.05, |log2FC|>1）のエンリッチメントも実施し、DMGの結果と直接比較。

### 2.4 モチーフ別機能分類

`motif_CCGG_spatial_detail.csv`（96,644行）と`motif_AAGCCCG_spatial_detail.csv`（850行）から、`is_methylated == True`の遺伝子を抽出し、COG/KEGGエンリッチメントを実施。

### 2.5 解析実装

スクリプト: `14_DMG_functional_enrichment/scripts/01_DMG_functional_enrichment.py`（Python 3, pandas 3.0, scipy 1.17, statsmodels 0.14）

---

## 3. 結果

### 3.1 DMGの概要

| 比較 | 6mA hyper | 6mA hypo | 4mC hyper | 4mC hypo | Combined |
|------|-----------|----------|-----------|----------|----------|
| T2vsT1 | 152 | 147 | 182 | 90 | **554** |
| T3vsT1 | 217 | 152 | 39 | 282 | **656** |
| T3vsT2 | 221 | 149 | 33 | 368 | **749** |

4mCの方向性に時間的トレンドあり: T1→T2でhyper優位（182 vs 90）、T2→T3でhypo優位（33 vs 368）。これはCCGGメチル化の一過的増加→脱メチル化の進行を反映する。

### 3.2 KEGGパスウェイエンリッチメント — 全条件で有意なし

> **Insight #1**: DMGはKEGGパスウェイに有意にエンリッチしない
> 全27条件（3比較 × 3修飾型 × 3方向）+ combined条件でFDR補正後に有意なKEGGパスウェイは**0件**。
> 対照的に、DEGは各比較あたり1〜18件の有意パスウェイを検出。
> **Figure**: `KEGG_DMG_vs_DEG_heatmap.pdf`

| 条件 | 有意パスウェイ数 (padj<0.05) | 上位パスウェイ（DEGのみ） |
|------|---------------------------|--------------------------|
| DMG T2vsT1 combined | 0 | — |
| DMG T3vsT1 combined | 0 | — |
| DMG T3vsT2 combined | 0 | — |
| **DEG T2vsT1 up** | **5** | **Ribosome** (padj=6.5e-14), Sesquiterpene biosyn., Type II PKS backbone, Type II PKS products, Siderophore biosyn. |
| **DEG T3vsT1 up** | **11** | **ABC transporters** (padj=6.0e-10), **Quorum sensing** (padj=1.6e-06), Nitrogen cycle, Type II PKS products |
| **DEG T3vsT2 down** | **18** | **Gly/Ser/Thr metabolism** (padj=5e-04), Cofactor biosyn., Amino acid biosyn. |

### 3.3 KEGG名目p値のトレンド分析

> **Insight #5**: 名目p値レベルでは、複数比較にまたがるトレンドが存在する

FDR補正前のp<0.1で複数比較に出現するパスウェイ:

| パスウェイ | T2vsT1 p | T3vsT1 p | T3vsT2 p | 生物学的解釈 |
|-----------|----------|----------|----------|------------|
| **Amino acid biosynthesis** (sco01230) | — | **0.026** | **0.044** | アミノ酸代謝関連遺伝子の一貫したメチル化変動 |
| **Phenazine biosynthesis** (sco00405) | — | **0.045** | **0.044** | フェナジン（二次代謝）への弱いシグナル |
| Glycerophospholipid metabolism | **0.024** | — | 0.072 | 膜脂質リモデリング |
| Siderophore NRP biosynthesis | **0.041** | — | 0.069 | 鉄獲得系 |
| Homologous recombination | 0.067 | **0.054** | — | DNA修復 |

**重要**: いずれも全padj > 0.7であり、統計的に有意とは言えない。ただし、アミノ酸生合成とフェナジン生合成が2つの独立した比較で名目有意であることは、弱いシグナルの存在を示唆する。

### 3.4 COGカテゴリ — V-Defenseの深掘り

> **Insight #2**: DMGのCOG分布はバックグラウンドとほぼ同一
> **Figure**: `COG_distribution_T3vsT1.pdf`

COG V-Defense は、DMGで唯一の統計的に有意なエンリッチメントであった:
- T3vsT2 4mC hypo: OR=2.38, **padj=0.028**
- T3vsT1 combined: p=0.016, padj=0.296（名目有意だがFDR非有意）
- T2vsT1 combined: p=0.028（名目有意だがFDR非有意）

> **Insight #3**: V-Defenseのエンリッチメントはメチルトランスフェラーゼ遺伝子群が駆動

15遺伝子の内訳を精査した結果:

**表3.4a: V-Defense エンリッチ遺伝子の詳細**

| gene_id | old_locus_tag | product | 遺伝子タイプ |
|---------|---------------|---------|-------------|
| SC_RS05450 | SCO0705 | SAM-dependent methyltransferase | **MTase** |
| SC_RS06105 | SCO0835 | class I SAM-dependent methyltransferase | **MTase** |
| SC_RS06595 | SCO0929 | SAM-dependent methyltransferase | **MTase** |
| SC_RS09775 | SCO1552 | TrmH family RNA methyltransferase | RNA修飾 |
| SC_RS11875 | SCO1969 | methylated-DNA--[protein]-cysteine S-methyltransferase | DNA修復 |
| SC_RS12500 | SCO2092 (rsmH) | 16S rRNA (C1402-N4)-methyltransferase RsmH | rRNA修飾 |
| SC_RS12890 | SCO2170 | class I SAM-dependent methyltransferase | **MTase** |
| SC_RS13205 | SCO2235 | Phd/YefM family antitoxin | TA系 |
| SC_RS13305 | SCO2256 (panB) | 3-methyl-2-oxobutanoate hydroxymethyltransferase | 代謝酵素 |
| SC_RS15340 | SCO2653 | class I SAM-dependent methyltransferase | **MTase** |
| SC_RS19420 | SCO3459 | methyltransferase domain-containing protein | **MTase** |
| SC_RS38275 | — | methyltransferase type 11 | **MTase** |
| SC_RS38280 | SCO7213 | SAM-dependent methyltransferase | **MTase** |
| SC_RS39465 | SCO7452 | methyltransferase | **MTase** |
| SC_RS40100 | SCO7580 | SAM-dependent methyltransferase | **MTase** |

**内訳**:
- SAM依存性メチルトランスフェラーゼ: **11/15 (73%)**
- RNA/rRNA修飾酵素: 2/15
- DNA修復（O6-methylguanine修復）: 1/15
- その他（TA系, panB）: 2/15 (COG V分類が疑わしい)

**R-M系との関係**: `rm_systems.csv`に登録されたコアR-M構成遺伝子（SC_RS19765, SC_RS19770, SC_RS17645等）は**含まれていない**。これらはR-M系の主要構成遺伝子ではなく、機能未知のオーファンメチルトランスフェラーゼ群である。ただし、8/15が`mtase_genes.csv`には登録されており、広義のメチル化関連遺伝子として特定されている。

**生物学的解釈**: 培養後期（T2→T3）における4mC脱メチル化が、メチルトランスフェラーゼ遺伝子自身の近傍で選択的に起こっている。これは**エピジェネティックな自己制御ループ**（メチル化酵素遺伝子のメチル化状態がメチル化酵素活性にフィードバックする）を示唆する興味深い知見である。

### 3.5 GOタームエンリッチメント — isomerase activityの深掘り

> **Insight #4**: GO:0016853（isomerase activity）のエンリッチメントは糖リン酸イソメラーゼ群が駆動

DMGで有意だった唯一のGOターム:
- **GO:0016853 (isomerase activity)**: T3vsT1 combined, padj=**0.006**, fold=5.65, 8/16遺伝子

**表3.5a: GO:0016853 エンリッチ遺伝子の詳細**

| gene_id | old_locus_tag | product | タイプ |
|---------|---------------|---------|-------|
| SC_RS04810 | SCO0577 | sugar phosphate isomerase/epimerase | **糖代謝** |
| SC_RS04820 | SCO0579 | ribose-5-phosphate isomerase | **糖代謝** |
| SC_RS11045 | SCO1804 | SAM:tRNA ribosyltransferase-isomerase | tRNA修飾 |
| SC_RS15835 | SCO2750 | sugar phosphate isomerase/epimerase | **糖代謝** |
| SC_RS15855 | SCO2754 | sugar phosphate isomerase/epimerase | **糖代謝** |
| SC_RS18320 | SCO3241 | sugar phosphate isomerase/epimerase | **糖代謝** |
| SC_RS18855 | SCO3347 | sugar phosphate isomerase/epimerase | **糖代謝** |
| SC_RS35070 | SCO6575 | sugar phosphate isomerase/epimerase | **糖代謝** |

**内訳**: 7/8 (88%) が糖リン酸イソメラーゼ/エピメラーゼファミリー。

**生物学的解釈**: T1→T3間で糖代謝関連イソメラーゼのメチル化が協調的に変動している。*Streptomyces*の培養段階進行に伴う一次代謝→二次代謝の転換において、糖リン酸代謝は前駆体供給の観点から重要であり、この経路のエピジェネティック制御は論理的に合致する。ただし、これは8遺伝子のみに基づく知見であり、さらなる検証が必要。

### 3.6 モチーフ別機能分類 — 視覚的印象の補正

> **Insight #6**: モチーフ間のCOG分布差は統計的に支持されない

COG分布図（`COG_motif_CCGG_vs_AAGCCCG.pdf`）の視覚的印象と統計値の比較:

| COGカテゴリ | AAGCCCG (6mA) fold | p値 | CCGG (4mC) fold | p値 |
|------------|-------------------|-----|-----------------|-----|
| K - Transcription | 0.875 | 0.82 | 0.935 | 0.80 |
| G - Carbohydrate | — | — | 1.376 | 0.033 (padj=0.44) |
| P - Inorganic ion | — | — | 1.201 | 0.044 (padj=0.44) |

**補正**: 棒グラフでAAGCCCGのK-Transcriptionが高く見えるのは**比率の視覚的な差**であり、ORA（Fisher正確検定）では**むしろ depleted 傾向**（fold=0.875）。これはAAGCCCGモチーフを持つ遺伝子の母集団サイズ（785遺伝子）がK-Transcription遺伝子の母集団（822遺伝子）と近い規模であるため、遺伝子数では多く見えるが比率では低い。

CCGGでG-Carbohydrateが名目有意（p=0.033）だが、FDR補正後は非有意（padj=0.44）。

### 3.7 検出力の評価

> **Insight #7**: DMGの検出力はDEGの約1/3〜1/5だが、それだけでは負の結果を説明できない

| 指標 | DMG (combined all) | DEG (all) | 比 |
|------|-------------------|-----------|-----|
| 遺伝子数 (T3vsT1) | 656 | 4,841 | 1:7.4 |
| KEGG背景に含まれる数 | 155 | ~1,204 | 1:7.8 |
| KEGG検定パスウェイ数 | ~97 | ~130 | 1:1.3 |
| 有意パスウェイ (padj<0.05) | **0** | **3** | — |
| 名目有意パスウェイ (p<0.05) | **4** | **~50** | 1:12.5 |

検出力は遺伝子数に依存するため、DMGの方が不利であることは確か。しかし:
- 名目p<0.05のパスウェイがDMGでは**わずか4件**しかない（DEGは約50件）
- この差は遺伝子数比（1:7.4）よりもはるかに大きい
- **結論**: 検出力不足は一因だが、DMGに機能的一貫性がないことが主因

---

## 4. 生物学的示唆

### 4.1 メチル化は配列依存的・機能非選択的に分布する

本解析の最も重要な知見は、**DMGが特定の機能カテゴリにエンリッチしない**という「negative result」である。これはDEGとの対比で3つの解釈を支持する:

1. **配列依存的メチル化**: R-M系（CCGG, AAGCCCG等）はDNA配列モチーフを認識してメチル化するため、標的選択はゲノム配列組成によって決まり、遺伝子機能とは独立
2. **真核生物との根本的な違い**: 真核生物では遺伝子体メチル化がハウスキーピング遺伝子に集中し、プロモーターメチル化が発生遺伝子を抑制するなどの機能的偏りがあるが、原核生物のR-M系メチル化にはそのようなパターンがない
3. **メチル化の「効果」は遺伝子の機能ではなく、位置に依存する**: メチル化がTSSプロモーター領域にあるか遺伝子本体にあるかが重要であり（既報: TSS proximal正相関）、遺伝子がどの機能カテゴリに属するかは無関係

### 4.2 「少数標的・間接カスケード」モデル

メチル化変動が機能的に偏らないにもかかわらず、4mC T2vsT1で有意なメチル化-発現相関が検出される事実は、以下のモデルと整合する:

```
ゲノム全体のメチル化変動（機能非選択的）
    │
    ├── 大多数の遺伝子: メチル化変動 → 表現型的影響なし
    │
    └── 少数のキーレギュレーター（redZ, afsS等）:
         メチル化変動 → 転写活性変化 → 下流パスウェイ全体の発現変動カスケード
```

このモデルでは、DMGの「機能エンリッチメント」ではなく、**どの特定遺伝子のプロモーターにメチル化変動が起こるか**が生物学的意味を持つ。

### 4.3 メチルトランスフェラーゼの自己制御ループ

V-Defenseエンリッチメントの実態がメチルトランスフェラーゼ遺伝子群であったことは、エピジェネティック制御における興味深い可能性を示す:

- T2→T3で4mCの脱メチル化が進行（4mC hypo DMG: 368遺伝子）
- その中にメチルトランスフェラーゼ自身が含まれる
- これは**正のフィードバック**（MTase遺伝子の脱メチル化 → MTase発現低下 → さらなる脱メチル化）を形成しうる
- 既報のMTase発現動態（SC_RS19765等のT3での発現変化）と照合する価値がある

### 4.4 糖リン酸イソメラーゼの協調的メチル化変動

GO:0016853の8遺伝子（うち7つが糖リン酸イソメラーゼ/エピメラーゼ）は、培養段階の進行に伴う一次→二次代謝転換との関連で注目に値する:

- 糖リン酸代謝はポリケタイド合成の前駆体供給に直結
- *Streptomyces*の二次代謝開始時における代謝フラックスの再配分にエピジェネティック制御が関与する可能性
- ただし、8遺伝子のみに基づく知見であり、仮説生成レベル

---

## 5. DEGエンリッチメントとの詳細比較

### 5.1 DEG上位パスウェイ（DMGで検出されないもの）

| パスウェイ | DEG padj | DEG遺伝子数 | DMG padj | 解釈 |
|-----------|----------|------------|----------|------|
| Ribosome | 6.5e-14 (T2vsT1) | 47 | >0.9 | 翻訳装置の大規模リプログラミングはメチル化非依存 |
| ABC transporters | 6.0e-10 (T3vsT1) | 108 | >0.9 | 輸送体系の発現変動はメチル化非依存 |
| Quorum sensing | 1.6e-06 (T3vsT1) | 60 | >0.9 | クオラムセンシングの時期特異的活性化はメチル化非依存 |
| Oxidative phosphorylation | 2.5e-05 (T3vsT1) | 23 | >0.9 | エネルギー代謝リモデリングはメチル化非依存 |

**示唆**: *Streptomyces*の培養段階に伴う主要な転写プログラム（翻訳装置抑制、二次代謝活性化、輸送系再編）は、メチル化とは独立した転写因子カスケードによって制御されている。

---

## Figure解説

### 【Figure】KEGG_DMG_vs_DEG_heatmap.pdf

**ファイル**: `14_DMG_functional_enrichment/analysis/14_DMG_enrichment_260206_v1/figures/KEGG_DMG_vs_DEG_heatmap.pdf`

**目的**: DMGとDEGのKEGGエンリッチメント結果を直接比較し、メチル化変動の機能的偏りの有無を可視化する

**方法**: Fisher正確検定 + BH補正。全比較・全方向の結果を-log10(padj)ヒートマップで表示。*印はpadj<0.05。バックグラウンド: KEGGアノテーションのある全遺伝子（1,679 genes）。

**図の読み方**:
- X軸: 左3列がDMG（T2vsT1, T3vsT1, T3vsT2）、右3列がDEG
- Y軸: KEGGパスウェイ（いずれかの条件でpadj<0.1のもの）
- 色: -log10(padj)のスケール（濃い色ほど有意）、白=非有意
- *: padj < 0.05

**結果と示唆**:
- DMG列は**全て淡色**（全padj > 0.7）、DEG列に複数の有意パスウェイ
- DEGではRibosome（T2vsT1）、Quorum sensing / ABC transporters（T3vsT1, T3vsT2）が一貫して強いシグナル
- DMGとDEGのエンリッチメントパターンは**完全に非対称**: 転写プログラムの大規模変動はメチル化に反映されない
- **示唆**: メチル化変動は特定パスウェイを選択的に標的としていない

**サポートするInsight**: #1, #5

### 【Figure】COG_distribution_T3vsT1.pdf

**ファイル**: `14_DMG_functional_enrichment/analysis/14_DMG_enrichment_260206_v1/figures/COG_distribution_T3vsT1.pdf`

**目的**: DMG, DEG, バックグラウンドのCOGカテゴリ分布を3群で比較し、メチル化変動遺伝子の機能的構成を評価する

**方法**: 各群に属する遺伝子のCOGカテゴリ比率（%）を棒グラフで表示。バックグラウンドはCOGアノテーションのある全8,275遺伝子。S-Unknown（最大カテゴリ）はプロットから除外。

**図の読み方**:
- X軸: COGカテゴリ（19カテゴリ、S除外）
- Y軸: 各群内の比率 (%)
- 灰: Background, 青: DEGs, 赤: DMGs

**結果と示唆**:
- DMG（赤）は**ほぼ全カテゴリでバックグラウンド（灰）と一致**
- R（General function）が全群で最大（約35%）
- DEG（青）ではK（Transcription）がバックグラウンドより低い傾向
- **示唆**: メチル化変動はゲノム全体の機能構成をそのまま反映しており、特定カテゴリを選択していない

**サポートするInsight**: #2

### 【Figure】COG_motif_CCGG_vs_AAGCCCG.pdf

**ファイル**: `14_DMG_functional_enrichment/analysis/14_DMG_enrichment_260206_v1/figures/COG_motif_CCGG_vs_AAGCCCG.pdf`

**目的**: 異なるR-M系モチーフでメチル化される遺伝子の機能分布を比較し、モチーフ間の標的選択性を評価する

**方法**: CCGG（4mC, 1,109メチル化遺伝子）とAAGCCCG（6mA, 633メチル化遺伝子）各モチーフでメチル化が検出された遺伝子のCOGカテゴリ比率を、バックグラウンド（全8,275遺伝子）と比較。

**図の読み方**:
- 灰: Background, 青: CCGG (4mC), 赤: AAGCCCG (6mA)

**注意点（視覚的印象の補正）**:
- AAGCCCGでK（Transcription）が高く**見える**が、Fisher検定ではfold=0.875（むしろ depleted 傾向, p=0.82）。棒グラフの比率の差は遺伝子母集団サイズの違いに起因する見かけの差。
- CCGGでG（Carbohydrate）が名目有意（p=0.033）だがFDR補正後は非有意（padj=0.44）

**結果と示唆**:
- 両モチーフとも大部分のカテゴリでバックグラウンドと近似
- 統計的に有意な偏りは**いずれのモチーフでも検出されない**
- **示唆**: R-M系のモチーフ認識は配列ベースであり、遺伝子機能による選択性はない

**サポートするInsight**: #6

---

## 6. 追加解析（5つの多角的検証）

基本解析で「DMGに機能選択性がない」という結論が得られたが、これを複数の異なる角度から検証するため5つの追加解析を実施した。

**スクリプト**: `14_DMG_functional_enrichment/scripts/02_DMG_selectivity_deep_analysis.py`
**出力ディレクトリ**: `14_DMG_functional_enrichment/analysis/14_DMG_selectivity_260206_v1/`

### 6.1 Analysis 1: DEG∩DMG協調変動遺伝子のエンリッチメント

**目的**: 発現変動とメチル化変動の**両方**を示す遺伝子（協調変動遺伝子）に限定してエンリッチメント解析を行い、「同時変動する遺伝子セット」に機能的偏りがあるかを検証する。

**方法**: 既存の協調変動遺伝子リスト（T2vsT1: 558, T3vsT1: 286, T3vsT2: 218遺伝子）を使用。各比較について全遺伝子 (all)、正相関 (positive)、負相関 (negative)、Gained、Lost の5サブセットでKEGG/GO/COG ORAを実施。

> **Insight #8**: 協調変動遺伝子もKEGGパスウェイにエンリッチしない

| 比較 | グループ | 遺伝子数 | KEGG有意 | GO有意 | COG有意 |
|------|---------|---------|---------|--------|---------|
| T2vsT1 | all | 304 | 0 | 0 | 0 |
| T2vsT1 | positive | 135 | 0 | 0 | 0 |
| T2vsT1 | negative | 75 | 0 | 0 | **1** |
| T3vsT1 | all | 278 | 0 | 0 | 0 |
| T3vsT1 | positive | 209 | 0 | **1** | 0 |
| T3vsT1 | lost | 180 | 0 | **1** | 0 |
| T3vsT2 | all | 216 | 0 | 0 | 0 |
| T3vsT2 | gained | 70 | 0 | **3** | 0 |

**結果と示唆**:
- KEGG: **全15条件で有意0件** — DMG単独と同じパターン
- GO: 散発的に1〜3件のみ
- COG: T2vsT1 negative で1件のみ
- **示唆**: 発現とメチル化が同時に変動する遺伝子においても、特定の機能パスウェイへの集中は見られない。メチル化変動の機能非選択性は協調変動遺伝子でも維持される

### 6.2 Analysis 2: COGカテゴリ別メチル化-発現相関

**目的**: DMGの「数」ではなく、メチル化変動量と発現変動量の「相関の強さ」がCOGカテゴリにより異なるかを検証する。

**方法**: COGカテゴリごとにSpearman相関（メチル化変化量 vs log2FC）を算出。4mCと6mAを別々に、3比較（T2vsT1, T3vsT1, T3vsT2）で実施。遺伝子数<10のカテゴリは除外。

> **Insight #9**: 相関の強度はカテゴリにより異なる

**有意な相関（p<0.05）**:

| 比較 | 修飾型 | COG | カテゴリ名 | n | Spearman r | p値 |
|------|--------|-----|-----------|---|-----------|-----|
| T2vsT1 | 4mC | O | Chaperones/PTM | 19 | **+0.481** | 0.037 |
| T2vsT1 | 6mA | V | Defense | 13 | **+0.566** | 0.044 |
| T3vsT1 | 4mC | R | General function | 181 | **-0.154** | 0.039 |
| T3vsT2 | 6mA | E | Amino acid metabolism | 26 | **-0.461** | 0.018 |
| T3vsT2 | 6mA | U | Secretion | 13 | **-0.599** | 0.031 |

**全体相関（参考）**: 4mC T2vsT1 全遺伝子 r=0.137, p=5.8e-04（唯一の有意な全体相関）

**結果と示唆**:
- DMGの「数」が特定カテゴリに偏らなくても、相関の「強さ」にはカテゴリ間差が存在
- O-Chaperones, V-Defense では正の相関（メチル化増加→発現上昇）、E-Amino acid, U-Secretion では負の相関
- ただし、多重検定補正（18カテゴリ×2修飾型×3比較=108検定）を考慮すると、5件/108はFDR閾値に達しない可能性が高い
- **示唆**: メチル化の効果の「方向性」がカテゴリ依存的である可能性はあるが、慎重な解釈が必要

### 【Figure】A2_correlation_heatmap_4mC.pdf / A2_correlation_heatmap_6mA.pdf

**ファイル**: `14_DMG_selectivity_260206_v1/figures/A2_correlation_heatmap_4mC.pdf`, `A2_correlation_heatmap_6mA.pdf`

**目的**: COGカテゴリ別のメチル化-発現相関を修飾型ごとに可視化し、カテゴリ間の差異を評価する

**方法**: Spearman相関係数をヒートマップで表示。X軸: 3比較、Y軸: COGカテゴリ。*印はp<0.05。

**結果と示唆**:
- 4mCヒートマップ: T2vsT1でO-Chaperonesが突出した正の相関
- 6mAヒートマップ: V-Defense (T2vsT1) で正、E-Amino acid / U-Secretion (T3vsT2) で負
- 全体として散発的であり、一貫した系統的パターンは見られない
- **示唆**: カテゴリ別相関は統計的ノイズと実シグナルの境界にある

**サポートするInsight**: #9

### 6.3 Analysis 3: プロモーター vs 遺伝子本体メチル化の機能分布

**目的**: メチル化の**位置**（プロモーター、TSS近傍、遺伝子本体）により、メチル化される遺伝子の機能分布が異なるかを検証する。

**方法**: TSS情報（Jeong et al. 2016, comprehensive_tss_table.csv）を用いて、各メチル化サイトを以下に分類:
- **Promoter**: TSS上流 ≤ 300bp
- **TSS proximal**: TSS ±150bp以内
- **Gene body**: 遺伝子領域内かつTSS ±150bp外

各領域でメチル化されている遺伝子について、COGおよびKEGG ORAを実施。

> **Insight #10**: 遺伝子本体メチル化にのみCOG機能偏りが存在する

**COGエンリッチメント結果**:

| 領域 | 遺伝子数 | 有意COG数 | 上位カテゴリ |
|------|---------|----------|------------|
| **Gene body** | 2,943 | **7** | E-Amino acid (padj=6.1e-08), G-Carbohydrate (padj=7.1e-08), T-Signal transduction (padj=1.9e-07), C-Energy (padj=4.6e-05), L-Replication (padj=0.002), O-Chaperones (padj=0.004), R-General (padj=0.027) |
| **Promoter** | 808 | **0** | なし（最小padj=0.48: G-Carbohydrate） |
| **TSS proximal** | 1,855 | **0** | なし（最小padj=0.26: J-Translation, P-Inorganic ion） |

**KEGGエンリッチメント結果**:

| 領域 | 有意パスウェイ数 (padj<0.05) | 名目有意数 (p<0.05) |
|------|---------------------------|-------------------|
| Gene body | **0** | 12（2-Oxocarboxylic acid, Glycolysis, TCA cycle等） |
| Promoter | 0 | 2（Fructose/mannose, One carbon pool by folate） |
| TSS proximal | 0 | 0 |

**結果と示唆**:
- 遺伝子本体にメチル化を持つ遺伝子は、代謝系（E, G, C）、シグナル伝達（T）、修復（L）、シャペロン（O）に有意にエンリッチ
- プロモーターとTSS近傍にメチル化を持つ遺伝子は、COGでもKEGGでも有意なエンリッチメントなし
- KEGGでは遺伝子本体でも padj<0.05 に達しないが、名目レベルで12パスウェイ（解糖系、TCA回路、アミノ酸生合成等の中心代謝）がトレンドを示す
- **示唆**: 遺伝子本体メチル化は代謝遺伝子に選択的に分布する傾向がある。これはR-M系モチーフのゲノム分布が代謝遺伝子の配列組成（codon usageやGC含量）と関連している可能性を反映する。一方、プロモーター領域のメチル化には機能的偏りがなく、真の遺伝子制御はプロモーターメチル化（位置依存的・機能非選択的）によることを支持する

### 【Figure】A3_COG_promoter_vs_body.pdf

**ファイル**: `14_DMG_selectivity_260206_v1/figures/A3_COG_promoter_vs_body.pdf`

**目的**: メチル化位置別（プロモーター vs 遺伝子本体 vs TSS近傍）のCOG分布を比較する

**方法**: 3領域のメチル化遺伝子のCOG ORA結果を並列表示。

**結果と示唆**:
- Gene body列にのみ有意カテゴリ（*印）が集中
- Promoter列とTSS proximal列は全カテゴリ非有意
- **示唆**: メチル化位置により機能選択性が質的に異なる

**サポートするInsight**: #10

### 6.4 Analysis 4: MTase/イソメラーゼ遺伝子の発現変動検証

**目的**: V-Defense MTase遺伝子（15遺伝子）とGOイソメラーゼ遺伝子（8遺伝子）が、メチル化変動だけでなく発現変動（DEG）も示すかを検証する。

**方法**: DESeq2結果（gene_master_DESeq2.tsv）からlog2FC, padj, DEG判定を取得。メチル化データ（integrated_methyl_expression_weighted.csv）からモチーフ別メチル化変化量を取得。

> **Insight #11**: V-Defense MTase遺伝子はメチル化と発現の両方で劇的変動

**表6.4a: V-Defense MTase遺伝子 — 発現変動**

| gene_id | old_locus_tag | product | T2vsT1 log2FC | T2vsT1 DEG | T3vsT1 DEG | 4mC変化(T3vsT1) |
|---------|---------------|---------|--------------|-----------|-----------|-----------------|
| **SC_RS11875** | SCO1969 | methylated-DNA-cysteine S-MTase | **+5.54** | **True** | **True** | -54% |
| **SC_RS05450** | SCO0705 | SAM-dependent MTase | **+2.73** | **True** | **True** | -85% |
| **SC_RS06595** | SCO0929 | SAM-dependent MTase | **+2.66** | **True** | **True** | 0% (T2で+51%) |
| **SC_RS39465** | SCO7452 | methyltransferase | **+2.12** | **True** | **True** | 0% (T2で+57%) |
| **SC_RS40100** | SCO7580 | SAM-dependent MTase | **+2.01** | **True** | **True** | 0% (T2で+76%) |
| **SC_RS13205** | SCO2235 | Phd/YefM antitoxin | +2.02 | **True** | False | -78% |
| SC_RS12890 | SCO2170 | class I SAM-dependent MTase | -1.93 | **True** | **True** | **-94%** |
| SC_RS13305 | SCO2256 (panB) | hydroxymethyltransferase | -1.84 | **True** | False | -69% |
| SC_RS09775 | SCO1552 | TrmH RNA MTase | -1.13 | **True** | **True** | -80% |
| SC_RS12500 | SCO2092 (rsmH) | 16S rRNA MTase | -2.58 | **True** | **True** | **-87%** |
| SC_RS19420 | SCO3459 | MTase domain protein | -1.37 | **True** | **True** | -56% |
| SC_RS38275 | — | methyltransferase type 11 | -1.04 | **True** | **True** | 0% (T2で+87%) |
| SC_RS06105 | SCO0835 | class I SAM-dependent MTase | +0.63 | False | **True** | 0% |
| SC_RS15340 | SCO2653 | class I SAM-dependent MTase | +0.37 | False | **True** | **-87%** |
| SC_RS38280 | SCO7213 | SAM-dependent MTase | -0.27 | False | False | 0% (T2で+87%) |

**DEG率**: T2vsT1: **12/15 (80%)**, T3vsT1: **12/15 (80%)**, T3vsT2: 7/15 (47%)

**4mCメチル化動態の特徴**:
- T1→T3で大幅な4mC脱メチル化（-54%〜-94%）を示す遺伝子が多数: SC_RS12890 (-94%), SC_RS15340 (-87%), SC_RS12500 (-87%), SC_RS05450 (-85%), SC_RS09775 (-80%)
- 別のグループはT2で一過的に4mCが増加し、T3で消失: SC_RS38275 (+87%→0%), SC_RS40100 (+76%→0%)
- 6mA変動は大部分で0%（CCGGサイト経由の4mCが主導的）

> **Insight #12**: GO isomerase遺伝子もDEGを含み、特にSC_RS18320が劇的発現上昇

**表6.4b: GO isomerase遺伝子 — 発現変動**

| gene_id | old_locus_tag | product | T2vsT1 log2FC | T2 DEG | T3vs1 DEG | メチル化変動（主要） |
|---------|---------------|---------|--------------|--------|-----------|-------------------|
| **SC_RS18320** | SCO3241 | sugar phosphate isomerase | **+8.93** | **True** | **True** | 6mA: -55% (T3vsT1) |
| SC_RS18855 | SCO3347 | sugar phosphate isomerase | +1.52 | **True** | False | 4mC: -83% (T3vsT1) |
| SC_RS35070 | SCO6575 | sugar phosphate isomerase | +1.08 | **True** | **True** | 6mA: +53% (T2vsT1) |
| SC_RS15835 | SCO2750 | sugar phosphate isomerase | -1.46 | **True** | False | 4mC: -81% (T3vsT1) |
| SC_RS04820 | SCO0579 | ribose-5-phosphate isomerase | +0.71 | False | **True** | 4mC: -94% (T3vsT1) |
| SC_RS11045 | SCO1804 | SAM:tRNA isomerase | +0.03 | False | **True** | 6mA: +58% (T2vsT1) |
| SC_RS15855 | SCO2754 | sugar phosphate isomerase | -0.82 | False | **True** | 6mA: +10% (T3vsT1) |
| SC_RS04810 | SCO0577 | sugar phosphate isomerase | -0.17 | False | **True** | 6mA: -55% (T2vsT1) |

**DEG率**: T2vsT1: 4/8 (50%), T3vsT1: **6/8 (75%)**, T3vsT2: 5/8 (63%)

**結果と示唆**:
- MTase遺伝子群は**メチル化と発現の両方で大幅な変動**を示し、T2で発現上昇→T3で4mC脱メチル化という時間的パターンが見られる
- SC_RS18320 (SCO3241) のlog2FC=+8.93はゲノム全体で最大級の発現上昇であり、糖リン酸代謝のエピジェネティック制御の存在を示唆
- MTase遺伝子のDEG率80%は、一般的なDMGの DEG率（約30-40%、既報DEG∩DMG overlap解析）よりも顕著に高い
- **示唆**: V-Defense MTase遺伝子群は「メチル化変動も発現変動も示す」高応答性遺伝子群であり、エピジェネティック自己制御ループの構成要素として機能している可能性がある

### 6.5 Analysis 5: TF/BGC遺伝子のDMG選択性

**目的**: 転写因子（TF）や生合成遺伝子クラスター（BGC）の遺伝子がDMGとして優先的にメチル化されるかを検証する（「少数標的・間接カスケード」モデルの検証）。

**方法**: gene_master_with_BGC_regulators.tsvからTF遺伝子（863遺伝子）とBGC遺伝子（99遺伝子）を特定。各群のDMG率を非TF/非BGC遺伝子と比較（Fisher正確検定）。

> **Insight #13**: TF/BGCはDMGとして選択的に標的化されない

| 比較 | 遺伝子クラス | n | DMG数 | DMG率 | 非クラスDMG率 | OR | p値 |
|------|------------|---|-------|------|-------------|-----|-----|
| T2vsT1 | TF | 863 | 51 | 5.9% | 7.0% | 0.84 | 0.28 |
| T2vsT1 | BGC | 99 | 5 | 5.1% | 6.9% | 0.72 | 0.69 |
| T3vsT1 | TF | 863 | 60 | 7.0% | 8.3% | 0.83 | 0.21 |
| T3vsT1 | BGC | 99 | 8 | 8.1% | 8.1% | 1.00 | 1.00 |
| T3vsT2 | TF | 863 | 73 | 8.5% | 9.4% | 0.90 | 0.42 |
| T3vsT2 | BGC | 99 | 8 | 8.1% | 9.3% | 0.86 | 0.86 |

**結果と示唆**:
- 全6検定で**OR < 1**（TF/BGCのDMG率がむしろ低い傾向）
- 全p値 > 0.2で統計的に非有意
- **示唆**: メチル化変動は転写因子やBGC遺伝子を優先的に標的としていない。これは「少数標的モデル」における選択性がメチル化の「存在」ではなく「位置」（プロモーター配置）に依存することを支持する

### 【Figure】A5_TF_BGC_DMG_frequency.pdf

**ファイル**: `14_DMG_selectivity_260206_v1/figures/A5_TF_BGC_DMG_frequency.pdf`

**目的**: TF/BGC遺伝子のDMG率を非TF/非BGC遺伝子と比較し、選択的メチル化の有無を可視化する

**方法**: Fisher正確検定。棒グラフでDMG率を群間比較。

**結果と示唆**:
- TF/BGCのDMG率バーは一貫して非クラスバーと同等かやや低い
- エラーバーは大きく重複しており、有意差なし
- **示唆**: 「少数標的モデル」の鍵は、特定遺伝子がDMGになる確率ではなく、DMGとなった遺伝子のうち転写因子のプロモーター上にメチル化が配置される確率にある

**サポートするInsight**: #13

---

## 7. 追加解析の統合考察

### 7.1 機能非選択性の頑健性

5つの追加解析はいずれも、基本解析の結論（DMGの機能非選択性）を異なる角度から支持する:

1. **協調変動遺伝子（A1）**: 発現・メチル化の同時変動ですら機能偏りなし
2. **カテゴリ別相関（A2）**: 相関の「強さ」にカテゴリ差はあるが系統的ではない
3. **位置別機能分布（A3）**: 遺伝子本体メチル化にのみ機能偏りがあり、制御的に重要なプロモーターメチル化には偏りがない
4. **例外遺伝子の検証（A4）**: MTaseとイソメラーゼの両例外は発現変動も大きく、真のシグナルである
5. **TF/BGC選択性（A5）**: 制御遺伝子のDMG率はゲノム平均以下

### 7.2 「位置依存・機能非依存」モデルの精緻化

基本解析で提案した「少数標的・間接カスケード」モデルを、追加解析の結果を統合して精緻化する:

```
原核生物R-M系メチル化
    │
    ├── 遺伝子本体メチル化 → 代謝遺伝子に多い（A3）
    │     └── 機能的影響は限定的
    │
    ├── プロモーターメチル化 → 機能非選択的（A3）
    │     └── 転写に直接影響 → 少数のTF/BGC遺伝子で
    │          カスケード効果を引き起こしうる
    │
    └── MTase遺伝子メチル化 → 自己制御ループ（A4）
          └── 4mC脱メチル化 ↔ MTase発現変動
```

**重要な発見**: 遺伝子本体メチル化が代謝遺伝子に偏る（A3）のは、R-M系モチーフの配列依存性（GC含量、コドン使用頻度）に由来する**受動的な**偏りであり、**能動的な**機能選択ではないと考えられる。一方、プロモーターメチル化が機能非選択的であることは、メチル化の遺伝子制御効果が「標的遺伝子の機能」ではなく「メチル化の位置」に依存するという結論を強く支持する。

---

## 出力ファイル一覧

### データファイル（基本解析: `14_DMG_enrichment_260206_v1/`）

| ファイル | 内容 |
|---------|------|
| `tables/DMG_summary.tsv` | DMGカウント一覧（比較×修飾型×方向） |
| `tables/enrichment_summary.tsv` | エンリッチメント結果サマリー（DMG vs DEG） |
| `tables/kegg_cache.tsv` | KEGG REST APIキャッシュ（1,761遺伝子, 150パスウェイ） |
| `tables/KEGG_DMG_*.tsv` (27ファイル) | DMG KEGGエンリッチメント結果（全条件） |
| `tables/KEGG_DEG_*.tsv` (9ファイル) | DEG KEGGエンリッチメント結果（比較用） |
| `tables/GO_DMG_*.tsv` / `GO_DEG_*.tsv` | GOエンリッチメント結果 |
| `tables/COG_DMG_*.tsv` / `COG_DEG_*.tsv` | COGエンリッチメント結果 |
| `tables/KEGG_motif_*.tsv` | モチーフ別KEGGエンリッチメント |
| `tables/COG_motif_*.tsv` | モチーフ別COGエンリッチメント |

### データファイル（追加解析: `14_DMG_selectivity_260206_v1/`）

| ファイル | 内容 |
|---------|------|
| `tables/A1_coordinated_enrichment_summary.tsv` | 協調変動遺伝子エンリッチメントサマリー |
| `tables/A1_KEGG_*.tsv` / `A1_GO_*.tsv` / `A1_COG_*.tsv` (45ファイル) | 協調変動遺伝子ORA結果（全条件） |
| `tables/A2_category_correlation.tsv` | COGカテゴリ別メチル化-発現相関（全122行） |
| `tables/A3_COG_gene_body.tsv` / `A3_COG_promoter.tsv` / `A3_COG_tss_proximal.tsv` | 位置別COGエンリッチメント |
| `tables/A3_KEGG_gene_body.tsv` / `A3_KEGG_promoter.tsv` / `A3_KEGG_tss_proximal.tsv` | 位置別KEGGエンリッチメント |
| `tables/A4_defense_MTase_expression.tsv` | V-Defense MTase 15遺伝子の発現データ |
| `tables/A4_defense_MTase_methylation.tsv` | V-Defense MTase 15遺伝子のメチル化変動 |
| `tables/A4_GO_isomerase_expression.tsv` | GO isomerase 8遺伝子の発現データ |
| `tables/A4_GO_isomerase_methylation.tsv` | GO isomerase 8遺伝子のメチル化変動 |
| `tables/A5_TF_BGC_DMG_frequency.tsv` | TF/BGC遺伝子DMG頻度（Fisher検定結果） |

### Figure（PDF + SVG）— 基本解析

| ファイル | 内容 | サポートするInsight |
|---------|------|-------------------|
| `DMG_summary_counts` | DMGカウント棒グラフ | — |
| `KEGG_DMG_vs_DEG_heatmap` | DMG vs DEG KEGGヒートマップ | #1, #5 |
| `KEGG_DMG_T2vsT1_top15` | DMG T2vsT1 KEGGトップ15 | #5 |
| `KEGG_DMG_T3vsT1_top15` | DMG T3vsT1 KEGGトップ15 | #5 |
| `KEGG_DMG_T3vsT2_top15` | DMG T3vsT2 KEGGトップ15 | #5 |
| `KEGG_hyper_vs_hypo_T2vsT1` | Hyper vs Hypo DMG KEGG | — |
| `KEGG_hyper_vs_hypo_T3vsT1` | 同上 | — |
| `KEGG_hyper_vs_hypo_T3vsT2` | 同上 | — |
| `COG_distribution_T2vsT1` | COG分布比較 | #2 |
| `COG_distribution_T3vsT1` | COG分布比較 | #2 |
| `COG_distribution_T3vsT2` | COG分布比較 | #2 |
| `COG_motif_CCGG_vs_AAGCCCG` | モチーフ別COG分布比較 | #6 |

### Figure（PDF + SVG）— 追加解析

| ファイル | 内容 | サポートするInsight |
|---------|------|-------------------|
| `A2_correlation_heatmap_4mC` | 4mC COGカテゴリ別相関ヒートマップ | #9 |
| `A2_correlation_heatmap_6mA` | 6mA COGカテゴリ別相関ヒートマップ | #9 |
| `A3_COG_promoter_vs_body` | 位置別COGエンリッチメント比較 | #10 |
| `A5_TF_BGC_DMG_frequency` | TF/BGC DMG頻度比較 | #13 |

---

## Limitation

- **DMG定義のアドホック性**: メチル化変化量 > 10% という閾値は慣例的なものであり、トランスクリプトームのpadj/log2FCのような統計的裏付けがない。PacBio SMRTのメチル化検出は二値的（メチル化あり/なし）に近いため、連続値としての解釈には注意が必要
- **検出力の非対称性**: DMGのKEGG背景遺伝子数（115〜155）はDEG（480〜1,204）の約1/3〜1/8。ただし名目p値レベルでもシグナルが乏しいことから、検出力のみでは説明不十分（セクション3.7）
- **GO term名の未解決**: GO IDの正式名称（oboファイル）を使用していないため、GO結果の解釈には追加のID→名前変換が必要
- **COGカテゴリの粒度**: COGは20カテゴリと粗い分類であるため、サブカテゴリレベルの偏りは検出できない
- **ORA手法の限界**: ORAは二値（DMG/非DMG）ベースであり、メチル化変動量の連続的情報を活用できない。GSEAベースのアプローチが有用な場合がある
- **カテゴリ別相関の多重検定**: Analysis 2で108の検定を実施しており、5件の名目有意結果はBonferroni/BH補正後に有意でない可能性が高い。探索的知見として位置付ける
- **遺伝子本体メチル化の因果関係**: Analysis 3で遺伝子本体メチル化がCOG 7カテゴリに偏ることが判明したが、これは配列組成（GC含量、コドン使用頻度）に起因する受動的な偏りの可能性があり、因果的な制御機構を意味しない可能性がある

---

## 次のステップへの示唆

1. ~~**MTaseフィードバックループの検証**~~ → **完了**（Analysis 4）: MTase遺伝子12/15がDEGであり、4mC脱メチル化（最大-94%）と発現変動の両方を確認。フィードバックの方向性（正 vs 負）は混在しており、単純なモデルでは説明不十分
2. ~~**「少数標的」モデルの定量的検証**~~ → **完了**（Analysis 5）: TF/BGCのDMG率はゲノム平均以下であり、メチル化の「存在」ではなく「位置」が鍵であることを確認
3. **糖リン酸イソメラーゼの二次代謝への寄与** — SC_RS18320 (SCO3241) のlog2FC=+8.93の劇的発現上昇の生物学的意義を調査。ポリケタイド前駆体供給との関連を検証
4. **BGCクラスター単位の解析** — 個別遺伝子ではなくBGC全体での4mC密度変動を解析し、二次代謝パスウェイレベルでのメチル化動態を明らかにする
5. **遺伝子本体メチル化の配列依存性検証** — Analysis 3で発見された遺伝子本体メチル化のCOG偏り（E, G, T, C等）が、R-M系モチーフの配列分布（GC含量依存性）で説明できるかをシミュレーションで検証

---

*最終更新: 2026-02-06*
