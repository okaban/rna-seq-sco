# DEGs機能エンリッチメント解析レポート

**解析内容**: ベン図、GO/KEGGエンリッチメント、COG分類、BGCカバレッジトラック

**日付**: 2026-02-02
**プロジェクト**: M145 RNA-seq
**対象生物**: *Streptomyces coelicolor* A3(2) M145
**解析ディレクトリ**: `12_supplementary_figures/analysis/12_supplementary_260202_v1/`

---

## Key Insights（主要な発見）

本解析から得られた重要なインサイトと、それを支持するFigureを以下にまとめる。

**Summary Figure**: `report_key_insights_summary.pdf` - 全インサイトの俯瞰図

| # | インサイト | 支持するFigure | 関連セクション |
|---|-----------|---------------|---------------|
| 1 | **T2が転写リプログラミングのピーク**: リボソーム抑制はT2で最大（3.34倍）、T3では緩和（2.76倍） | `report_FE_change_patterns.pdf`, `KEGG_bubble_T2vsT1_down.pdf` | §2 |
| 2 | **Prodigiosin生合成の早期誘導**: T2で既に強く誘導（2.83倍）、T3では維持（1.83倍） | `report_FE_change_patterns.pdf`, `report_BGC_expression_heatmap.pdf` | §2, §5 |
| 3 | **約1,000遺伝子がコア応答遺伝子**: 3比較すべてで変動する遺伝子群が存在 | `venn_all_DEGs.pdf` | §1 |
| 4 | **成長-二次代謝トレードオフ**: 翻訳カテゴリUp:Down = 1:6、二次代謝は上昇のみ | `report_COG_net_change.pdf`, `COG_barplot_T3vsT1.pdf` | §4 |
| 5 | **BGCクラスター間の誘導タイミング差異**: red（早期）→ cda（中期）→ act/cpk（後期） | `report_BGC_expression_heatmap.pdf`, `coverage_tracks_all_BGC.pdf` | §5 |
| 6 | **栄養枯渇への適応応答**: ABC輸送体108遺伝子、シデロフォア合成の発現上昇 | `report_pathway_pattern_summary.pdf`, `KEGG_bubble_T3vsT1_up.pdf` | §2, §3 |

---

## 解析の背景と目的

トランスクリプトーム解析の事例ファイル（`00_memo/rna-seq_workflow_260202.md`）と本プロジェクトを比較し、不足していた以下の4つの解析を追加実施した：

1. **ベン図（Venn Diagram）**: 3比較間のDEGs重複解析
2. **GO/KEGGエンリッチメント解析**: 機能・パスウェイの統計的濃縮
3. **COG分類解析**: 機能カテゴリ別の遺伝子分布
4. **ゲノムカバレッジトラック**: BGC領域のリードカバレッジ可視化

---

## 1. ベン図解析：DEGsの経時的パターン

### 主要な発見

| カテゴリ | 遺伝子数 | 解釈 |
|---------|---------|------|
| **3比較すべてで変動** | 1,041 | コア応答遺伝子（培養全期間で持続的に変動） |
| **T2vsT1 & T3vsT1のみ** | 1,793 | T1からの累積的変化（培養初期からの一貫した変動） |
| **T3vsT1 & T3vsT2のみ** | 1,486 | 後期特異的変化（T3で特に顕著に変動） |
| **T3vsT1のみ** | 521 | 長期培養特異的遺伝子 |

### インサイト

> **Insight #3**: 約1,000遺伝子がコア応答遺伝子として培養全期間を通じて変動し続けている。これは*S. coelicolor*の成長段階依存的な代謝シフトの根幹をなす遺伝子群である。
>
> **Figure**: `venn_all_DEGs.pdf`

- 後期（T3）に向けて変動遺伝子数が増加（T3vsT1で4,841 DEGs、最大）
- Up-regulated とDown-regulated で重複パターンが異なる（`venn_up_DEGs.pdf`, `venn_down_DEGs.pdf`）

### 出力ファイル
| ファイル | 内容 |
|---------|------|
| `figures/venn_all_DEGs.pdf` | 全DEGs（3比較の重複） |
| `figures/venn_up_DEGs.pdf` | Up-regulated DEGs |
| `figures/venn_down_DEGs.pdf` | Down-regulated DEGs |
| `tables/venn_DEGs_summary.tsv` | 重複カテゴリ別統計 |

---

## 2. KEGG パスウェイエンリッチメント解析

### 培養後期（T3 vs T1）での主要な変化

#### Up-regulated パスウェイ（発現上昇）

| パスウェイ | 遺伝子数 | Fold Enrichment | p-value | 意義 |
|-----------|---------|-----------------|---------|------|
| **ABC transporters** | 108 | 1.66 | 1.0e-10 | 栄養枯渇に対する輸送能力の増強 |
| **Quorum sensing** | 60 | 1.77 | 8.5e-8 | 細胞密度依存的な調節の活性化 |
| **Nitrogen cycle** | 12 | 2.74 | 1.9e-5 | 窒素代謝の調整 |
| **Type II polyketide products** | 10 | 2.70 | 1.4e-4 | ポリケタイド系抗生物質生合成 |
| **Type II polyketide backbone** | 6 | 2.97 | 1.4e-3 | PKS骨格合成の誘導 |
| **Prodigiosin biosynthesis** | 16 | 1.83 | 3.2e-3 | 赤色色素（undecylprodigiosin）生合成 |

#### Down-regulated パスウェイ（発現低下）

| パスウェイ | 遺伝子数 | Fold Enrichment | p-value | 意義 |
|-----------|---------|-----------------|---------|------|
| **Ribosome** | 55 | 2.76 | 1.1e-25 | タンパク質合成能力の低下（成長減速） |
| **Metabolic pathways** | 408 | 1.50 | 4.8e-23 | 一次代謝全般の抑制 |
| **Biosynthesis of cofactors** | 91 | 2.17 | 1.8e-16 | 補因子合成の調整 |
| **Biosynthesis of amino acids** | 79 | 2.09 | 3.1e-13 | アミノ酸合成の抑制 |
| **Nucleotide metabolism** | 35 | 2.79 | 2.0e-11 | DNA/RNA合成の減速 |
| **Pyrimidine metabolism** | 29 | 3.06 | 2.4e-11 | ピリミジン合成の大幅抑制 |

### 培養中期（T2 vs T1）での主要な変化

#### Up-regulated パスウェイ（発現上昇）

| パスウェイ | 遺伝子数 | Fold Enrichment | p-value | 意義 |
|-----------|---------|-----------------|---------|------|
| **Prodigiosin biosynthesis** | 18 | **2.83** | 1.8e-6 | 赤色色素生合成の早期誘導 |
| **Biosynthesis of secondary metabolites** | 170 | 1.35 | 4.0e-6 | 二次代謝全般の活性化開始 |
| **Sulfur metabolism** | 12 | 2.58 | 3.9e-4 | 硫黄代謝の活性化 |
| **Pyruvate metabolism** | 25 | 1.82 | 7.4e-4 | 中心炭素代謝の調整 |
| **Siderophore biosynthesis** | 6 | 3.06 | 3.8e-3 | 鉄獲得システムの活性化 |
| **Sesquiterpenoid biosynthesis** | 4 | 4.08 | 3.6e-3 | テルペノイド生合成 |

#### Down-regulated パスウェイ（発現低下）

| パスウェイ | 遺伝子数 | Fold Enrichment | p-value | 意義 |
|-----------|---------|-----------------|---------|------|
| **Ribosome** | 59 | **3.34** | 8.0e-35 | タンパク質合成の顕著な抑制 |
| **Biosynthesis of cofactors** | 61 | 1.64 | 1.5e-5 | 補因子合成の調整 |
| **Pyrimidine metabolism** | 20 | 2.39 | 3.2e-5 | ピリミジン合成の抑制 |
| **Nucleotide metabolism** | 24 | 2.16 | 4.6e-5 | ヌクレオチド合成の抑制 |
| **Thiamine metabolism** | 9 | 3.05 | 4.5e-4 | チアミン合成の抑制 |
| **DNA replication** | 11 | 2.43 | 1.7e-3 | DNA複製の減速 |

### インサイト

> **Insight #1**: T2が転写リプログラミングのピークである。リボソーム遺伝子の抑制はT2で最大（3.34倍濃縮）であり、T3（2.76倍）では緩和される。これは移行期（T2）で最も劇的な代謝シフトが起きていることを示す。
>
> **Figure**: `KEGG_bubble_T2vsT1_down.pdf` vs `KEGG_bubble_T3vsT1_down.pdf`

> **Insight #2**: Prodigiosin生合成はT2で既に強く誘導されている（2.83倍）。T3では維持されるが濃縮度は低下（1.83倍）。redクラスターは最も早期に活性化するBGCである。
>
> **Figure**: `KEGG_bubble_T2vsT1_up.pdf`, `coverage_track_red.pdf`

> **Insight #6**: 栄養枯渇への適応応答として、ABC輸送体（108遺伝子）、シデロフォア合成、クオラムセンシング関連遺伝子が発現上昇している。
>
> **Figure**: `KEGG_bubble_T3vsT1_up.pdf`

### T3 vs T1 と T2 vs T1 の詳細比較

#### 比較表のピックアップ基準

以下の基準でパスウェイを網羅的に分類した：

1. **両方で同方向に有意**: T2 vs T1 と T3 vs T1 の両方で p < 0.05、かつ同じ方向（Up/Down）
2. **方向反転**: T2とT3で発現変動の方向が逆転（Up→Down または Down→Up）
3. **時期特異的**: 一方の比較でのみ有意（p < 0.05）
4. **Fold Enrichment変化**: 時間経過での濃縮度の増減を「↑」「↓」「→」で表記

**凡例**: FE = Fold Enrichment、n.s. = not significant (p ≥ 0.05)

---

#### パスウェイ変動パターンの全体像

**Figure**: `report_pathway_pattern_summary.pdf` - 63パスウェイの変動パターン分類

##### 同方向パターン（Up/Down分類が適用可能）

| パターン | Up | Down | 計 |
|---------|-----|------|-----|
| 両方で同方向に有意 | 3 | 14 | 17 |
| T2のみ有意 | 7 | 0 | 7 |
| T3のみ有意 | 9 | 28 | 37 |
| **小計** | 19 | 42 | **61** |

##### 方向反転パターン（T2とT3で変動方向が逆転）

| パターン | パスウェイ数 |
|---------|-------------|
| T2 Up → T3 Down | 2 |
| T2 Down → T3 Up | 0 |
| **小計** | **2** |

**総計: 63パスウェイ**

> **注**: T2でのみDown-regulatedのパスウェイは0件。これはT2で抑制が始まったパスウェイはすべてT3でも持続的に抑制されていることを示す。

---

#### 方向反転パスウェイ（T2 Up → T3 Down）

**Figure**: `report_FE_change_patterns.pdf` - T2→T3のFold Enrichment変化を矢印で可視化

**最も重要な変動パターン**: T2で一過性に誘導された後、T3で抑制に転じる

| パスウェイ | T2 vs T1 (Up) | T3 vs T1 (Down) | 解釈 |
|-----------|---------------|-----------------|------|
| | 遺伝子数 (FE) p値 | 遺伝子数 (FE) p値 | |
| **Metabolic pathways** | 302 (1.16) 1.1e-3 | 408 (1.50) **4.8e-23** | 一次代謝は一時的に活性化後、後期で大幅抑制 |
| **Biosynthesis of secondary metabolites** | 170 (1.35) **4.0e-6** | 202 (1.53) 1.4e-12 | 二次代謝全般はT2で誘導開始、T3では特定BGCに分化 |

> **Insight**: T2は代謝の「移行期」であり、一次代謝と二次代謝の両方が一時的に活性化する。T3では一次代謝が抑制され、特定の二次代謝経路（Type II PKS等）のみが活性を維持する。

---

#### Up-regulated パスウェイの比較

##### 両比較で有意なUp-regulated（3パスウェイ）

| パスウェイ | T2 vs T1 | T3 vs T1 | 変化 | 解釈 |
|-----------|----------|----------|------|------|
| | 遺伝子数 (FE) p値 | 遺伝子数 (FE) p値 | | |
| **Prodigiosin biosynthesis** | 18 (**2.83**) 1.8e-6 | 16 (1.83) 3.2e-3 | ↓ | T2で最大誘導、T3で維持 |
| **Siderophore biosynthesis** | 6 (**3.06**) 3.8e-3 | 7 (2.60) 2.7e-3 | ↓ | 鉄獲得は早期から活性化 |
| **Starch and sucrose metabolism** | 22 (1.66) 6.0e-3 | 30 (1.65) 7.4e-4 | → | 糖代謝は持続的に活性 |

##### T2でのみUp-regulated（早期一過性応答、7パスウェイ）

| パスウェイ | T2 vs T1 | T3 vs T1 | 解釈 |
|-----------|----------|----------|------|
| | 遺伝子数 (FE) p値 | | |
| **Sulfur metabolism** | 12 (2.58) 3.9e-4 | n.s. | 硫黄代謝の一過性活性化 |
| **Pyruvate metabolism** | 25 (1.82) 7.4e-4 | n.s. | 中心炭素代謝の一時的調整 |
| **Valine, leucine, isoleucine biosynthesis** | 10 (2.40) 2.6e-3 | n.s. | 分岐鎖アミノ酸の早期合成 |
| **Sulfur cycle** | 6 (3.06) 3.8e-3 | n.s. | 硫黄循環の一過性活性化 |
| **Sesquiterpenoid biosynthesis** | 4 (**4.08**) 3.6e-3 | n.s. | テルペノイド生合成の早期活性 |
| **Pantothenate and CoA biosynthesis** | 13 (1.90) 9.3e-3 | n.s. | 補酵素A合成の一過性活性 |
| **beta-Alanine metabolism** | 10 (2.04) 1.2e-2 | n.s. | β-アラニン代謝の一過性活性 |

##### T3でのみUp-regulated（後期特異的応答、9パスウェイ）

| パスウェイ | T2 vs T1 | T3 vs T1 | 解釈 |
|-----------|----------|----------|------|
| | | 遺伝子数 (FE) p値 | |
| **ABC transporters** | n.s. | 108 (1.66) **1.0e-10** | 栄養枯渇への適応輸送 |
| **Quorum sensing** | n.s. | 60 (1.77) **8.5e-8** | 細胞密度依存的調節 |
| **Nitrogen cycle** | n.s. | 12 (2.74) 1.9e-5 | 窒素代謝の後期活性化 |
| **Type II PKS products** | n.s. | 10 (2.70) 1.4e-4 | actinorhodin等の後期誘導 |
| **Type II PKS backbone** | n.s. | 6 (2.97) 1.4e-3 | PKS骨格合成 |
| **Ascorbate and aldarate metabolism** | n.s. | 11 (2.18) 1.9e-3 | 後期代謝調整 |
| **Tetracycline biosynthesis** | n.s. | 7 (2.60) 2.7e-3 | 後期二次代謝 |
| **Degradation of flavonoids** | n.s. | 8 (2.38) 3.6e-3 | フラボノイド分解 |
| **Other glycan degradation** | n.s. | 9 (2.06) 9.4e-3 | 多糖分解の後期活性 |

---

#### Down-regulated パスウェイの比較

##### 両比較で有意なDown-regulated（14パスウェイ）

**Figure**: `report_both_down_FE_comparison.pdf` - 14パスウェイのFE比較ドットプロット

| パスウェイ | T2 vs T1 | T3 vs T1 | 変化 | 解釈 |
|-----------|----------|----------|------|------|
| | 遺伝子数 (FE) p値 | 遺伝子数 (FE) p値 | | |
| **Ribosome** | 59 (**3.34**) 8.0e-35 | 55 (2.76) 1.1e-25 | ↓ | T2で最大抑制、T3で緩和 |
| **Biosynthesis of cofactors** | 61 (1.64) 1.5e-5 | 91 (**2.17**) 1.8e-16 | ↑ | 抑制が持続・強化 |
| **Pyrimidine metabolism** | 20 (2.39) 3.2e-5 | 29 (**3.06**) 2.4e-11 | ↑ | 抑制が持続・強化 |
| **Nucleotide metabolism** | 24 (2.16) 4.6e-5 | 35 (**2.79**) 2.0e-11 | ↑ | 抑制が持続・強化 |
| **Homologous recombination** | 13 (2.49) 4.5e-4 | 16 (**2.72**) 1.2e-5 | ↑ | DNA修復の持続的抑制 |
| **Thiamine metabolism** | 9 (**3.05**) 4.5e-4 | 10 (3.01) 1.5e-4 | → | 持続的な抑制 |
| **Aminoacyl-tRNA biosynthesis** | 15 (0.70) 1.4e-3 | 24 (0.99) 2.1e-9 | ↑ | tRNA合成の持続的抑制 |
| **DNA replication** | 11 (2.43) 1.7e-3 | 15 (**2.93**) 5.0e-6 | ↑ | 抑制が持続・強化 |
| **Amino sugar/nucleotide sugar** | 22 (1.73) 3.7e-3 | 27 (**1.89**) 2.1e-4 | ↑ | 持続的な抑制 |
| **Biosynthesis of nucleotide sugars** | 18 (1.76) 6.7e-3 | 23 (**2.00**) 2.0e-4 | ↑ | 抑制が持続・強化 |
| **Ubiquinone biosynthesis** | 7 (2.57) 8.1e-3 | 8 (**2.61**) 3.2e-3 | → | 持続的な抑制 |
| **Alanine, aspartate, glutamate** | 15 (1.84) 8.3e-3 | 19 (**2.06**) 4.3e-4 | ↑ | 抑制が持続・強化 |
| **Purine metabolism** | 25 (1.55) 1.1e-2 | 40 (**2.20**) 3.2e-8 | ↑ | 抑制が持続・強化 |
| **Mismatch repair** | 9 (2.09) 1.5e-2 | 13 (**2.68**) 1.1e-4 | ↑ | 抑制が持続・強化 |

##### T2でのみDown-regulated（早期一過性抑制）

**該当パスウェイ: 0件**

> T2で抑制が始まったパスウェイはすべてT3でも持続している。これは一次代謝・成長関連経路の抑制が不可逆的であることを示唆する。

##### T3でのみDown-regulated（後期特異的抑制、28パスウェイ）

上位12パスウェイを示す（全28パスウェイはtables/KEGG_enrichment_all.tsvを参照）:

| パスウェイ | T2 vs T1 | T3 vs T1 | 解釈 |
|-----------|----------|----------|------|
| | | 遺伝子数 (FE) p値 | |
| **Biosynthesis of amino acids** | n.s. | 79 (2.09) **3.1e-13** | アミノ酸合成の後期抑制 |
| **Glycine, serine, threonine** | n.s. | 30 (2.67) **3.3e-9** | アミノ酸代謝の後期抑制 |
| **Porphyrin metabolism** | n.s. | 28 (2.61) 2.4e-8 | ポルフィリン合成の後期抑制 |
| **Oxidative phosphorylation** | n.s. | 36 (2.10) 7.5e-7 | エネルギー代謝の後期抑制 |
| **Peptidoglycan biosynthesis** | n.s. | 23 (2.73) 1.3e-7 | 細胞壁合成停止 |
| **Carbon metabolism** | n.s. | 62 (1.64) 9.1e-6 | 炭素代謝の後期抑制 |
| **Base excision repair** | n.s. | 15 (2.79) 1.3e-5 | DNA修復の後期抑制 |
| **RNA degradation** | n.s. | 11 (2.87) 4.9e-5 | RNA分解の後期抑制 |
| **One carbon pool by folate** | n.s. | 15 (2.55) 7.1e-5 | 葉酸代謝の後期抑制 |
| **Carbon fixation** | n.s. | 14 (2.49) 1.8e-4 | 炭素固定の後期抑制 |
| **Pentose phosphate pathway** | n.s. | 21 (2.05) 2.4e-4 | PPP経路の後期抑制 |
| **Glycolysis / Gluconeogenesis** | n.s. | 26 (1.75) 1.1e-3 | 解糖系の後期抑制 |

---

#### 比較から得られる生物学的インサイト

**Figures**: `report_pathway_pattern_summary.pdf`, `report_FE_change_patterns.pdf`, `report_both_down_FE_comparison.pdf`

| パターン | 観察 | 解釈 |
|---------|------|------|
| **方向反転** | Metabolic pathways: T2 Up → T3 Down | **T2は代謝の移行期**。一次・二次代謝が一時的に共活性化 |
| **方向反転** | Secondary metabolites: T2 Up → T3 Down | T3では特定BGCのみ活性、全般的な二次代謝は抑制 |
| **T2 > T3 抑制** | Ribosome: FE 3.34 → 2.76 | **T2が転写リプログラミングのピーク** |
| **T2 > T3 誘導** | Prodigiosin: FE 2.83 → 1.83 | **red BGCは早期誘導型** |
| **T3のみUp** | ABC/Quorum sensing/Type II PKS | **栄養枯渇応答・後期BGCは後期に顕在化** |
| **T3のみDown (28)** | アミノ酸合成、エネルギー代謝等 | **一次代謝の抑制は後期で大幅拡大** |
| **T2のみDown = 0** | 全T2 Downが T3でも持続 | **成長抑制は不可逆的** |

### T3 vs T2 を実施しなかった理由

1. **生物学的解釈の焦点**: T1を基準とした比較が最も解釈しやすい
2. **冗長性の回避**: T3 vs T2 は T3 vs T1 と T2 vs T1 の差分として理解可能
3. **解析自体は実施済み**: DESeq2結果は`04_deseq2/`に存在、必要に応じて追加可能

### 出力ファイル
| ファイル | 内容 |
|---------|------|
| `figures/KEGG_bubble_T3vsT1_up.pdf` | T3vsT1 Up バブルチャート |
| `figures/KEGG_bubble_T3vsT1_down.pdf` | T3vsT1 Down バブルチャート |
| `figures/KEGG_bubble_T2vsT1_up.pdf` | T2vsT1 Up バブルチャート |
| `figures/KEGG_bubble_T2vsT1_down.pdf` | T2vsT1 Down バブルチャート |
| `figures/report_pathway_pattern_summary.pdf` | **パスウェイ変動パターン全体像（論文用）** |
| `figures/report_FE_change_patterns.pdf` | **FE変化矢印図（論文用）** |
| `figures/report_both_down_FE_comparison.pdf` | **持続的Down-regulatedパスウェイ比較（論文用）** |
| `tables/KEGG_enrichment_all.tsv` | 全KEGG結果（82 pathways） |

---

## 3. GO エンリッチメント解析

### 培養後期（T3 vs T1）の特徴

#### Up-regulated（24 GO terms）

| GO term | 遺伝子数 | Fold Enrichment | 機能 |
|---------|---------|-----------------|------|
| GO:0008643 | 58 | 2.33 | carbohydrate transport |
| GO:0055085 | 153 | 1.44 | transmembrane transport |
| GO:0042626 | 122 | 1.50 | ATPase-coupled transport |
| GO:0016491 | 129 | 1.47 | oxidoreductase activity |
| GO:0004553 | 32 | 1.97 | hydrolase activity (glycosyl) |

#### Down-regulated（53 GO terms）

| GO term | 遺伝子数 | Fold Enrichment | 機能 |
|---------|---------|-----------------|------|
| GO:0006412 | 60 | **3.10** | translation |
| GO:0003735 | 54 | **3.11** | structural constituent of ribosome |
| GO:0005840 | 30 | 2.85 | ribosome |
| GO:0003723 | 29 | 2.62 | RNA binding |
| GO:0009252 | 21 | 2.46 | peptidoglycan biosynthesis |
| GO:0006260 | 17 | 2.39 | DNA replication |

### インサイト

> **翻訳関連遺伝子の劇的な発現低下**（3.1倍濃縮）は、KEGGのRibosomeパスウェイの結果と一致し、細胞が増殖モードから生存モードへ移行していることを裏付ける。
>
> **Figure**: `GO_bubble_T3vsT1_down.pdf`, `GO_bar_T3vsT1_down.pdf`

### 出力ファイル
| ファイル | 内容 |
|---------|------|
| `figures/GO_bubble_T3vsT1_up.pdf` | T3vsT1 Up バブルチャート |
| `figures/GO_bubble_T3vsT1_down.pdf` | T3vsT1 Down バブルチャート |
| `figures/GO_bar_T3vsT1_up.pdf` | T3vsT1 Up バーチャート |
| `figures/GO_bar_T3vsT1_down.pdf` | T3vsT1 Down バーチャート |
| `tables/GO_enrichment_all.tsv` | 全GO結果（115 terms） |

---

## 4. COG機能分類解析

**注**: NCBI COG 2024には26の機能カテゴリが定義されているが、本データセットで検出されたのは20カテゴリ（A, B, W, X, Y, Zは真核生物特異的または可動遺伝因子関連のため細菌では稀）。T3 vs T1では19カテゴリ、T2 vs T1では20カテゴリ（Nが追加）が検出された。

### T3 vs T1 での機能カテゴリ別変動（全19カテゴリ）

| カテゴリ | Up | Down | Net | 解釈 |
|---------|-----|------|-----|------|
| **R - General function** | 1064 | 645 | **+419** | 多機能遺伝子の活性化 |
| **S - Unknown** | 443 | 319 | +124 | 未知機能遺伝子の活性化 |
| **P - Inorganic ion transport** | 227 | 102 | **+125** | イオン輸送の活性化 |
| **C - Energy production** | 205 | 120 | +85 | 代謝リモデリング |
| **H - Coenzyme metabolism** | 76 | 26 | +50 | 補酵素合成の活性化 |
| **G - Carbohydrate metabolism** | 81 | 49 | +32 | 糖代謝の活性化 |
| **I - Lipid metabolism** | 32 | 14 | +18 | 脂質代謝の活性化 |
| **Q - Secondary metabolism** | 12 | 0 | **+12** | 二次代謝の一方向的誘導 |
| **V - Defense** | 53 | 41 | +12 | 防御機構の活性化 |
| **M - Cell envelope** | 36 | 31 | +5 | 細胞壁維持（ほぼ均衡） |
| **O - Protein modification** | 61 | 56 | +5 | タンパク質修飾（ほぼ均衡） |
| **T - Signal transduction** | 100 | 97 | +3 | シグナル伝達（ほぼ均衡） |
| **U - Secretion** | 42 | 46 | -4 | 分泌系（ほぼ均衡） |
| **D - Cell division** | 1 | 10 | -9 | 細胞分裂の停止 |
| **F - Nucleotide metabolism** | 7 | 17 | -10 | ヌクレオチド合成の抑制 |
| **E - Amino acid metabolism** | 95 | 107 | -12 | アミノ酸代謝の抑制傾向 |
| **L - Replication/Repair** | 43 | 59 | -16 | DNA複製の抑制 |
| **K - Transcription** | 168 | 218 | -50 | 転写調節の変動 |
| **J - Translation** | 20 | 118 | **-98** | リボソーム合成の大幅抑制 |

### T2 vs T1 での機能カテゴリ別変動（全20カテゴリ）

| カテゴリ | Up | Down | Net | 解釈 |
|---------|-----|------|-----|------|
| **R - General function** | 807 | 587 | **+220** | 多機能遺伝子の活性化 |
| **S - Unknown** | 443 | 280 | +44 | 未知機能遺伝子の活性化 |
| **C - Energy production** | 137 | 109 | +28 | エネルギー代謝の活性化 |
| **H - Coenzyme metabolism** | 46 | 26 | +20 | 補酵素合成の活性化 |
| **O - Protein modification** | 51 | 36 | +15 | タンパク質修飾の活性化 |
| **V - Defense** | 45 | 33 | +12 | 防御機構の活性化 |
| **Q - Secondary metabolism** | 9 | 0 | **+9** | 二次代謝の一方向的誘導 |
| **U - Secretion** | 43 | 34 | +9 | 分泌系の活性化 |
| **I - Lipid metabolism** | 19 | 12 | +7 | 脂質代謝の活性化 |
| **T - Signal transduction** | 79 | 76 | +3 | シグナル伝達（ほぼ均衡） |
| **E - Amino acid metabolism** | 83 | 81 | +2 | アミノ酸代謝（ほぼ均衡） |
| **F - Nucleotide metabolism** | 8 | 9 | -1 | ヌクレオチド代謝（ほぼ均衡） |
| **N - Cell motility** | 0 | 1 | -1 | 細胞運動性の抑制 |
| **M - Cell envelope** | 22 | 25 | -3 | 細胞壁（ほぼ均衡） |
| **D - Cell division** | 2 | 8 | -6 | 細胞分裂の抑制傾向 |
| **G - Carbohydrate metabolism** | 42 | 53 | -11 | 糖代謝の抑制傾向 |
| **P - Inorganic ion transport** | 109 | 124 | -15 | イオン輸送の抑制傾向 |
| **L - Replication/Repair** | 26 | 46 | -20 | DNA複製の抑制 |
| **K - Transcription** | 134 | 205 | -71 | 転写調節の大幅変動 |
| **J - Translation** | 17 | 100 | **-83** | リボソーム合成の大幅抑制 |

### T2 vs T1 → T3 vs T1 の変化傾向

| カテゴリ | T2vsT1 Net | T3vsT1 Net | 変化 | 解釈 |
|---------|------------|------------|------|------|
| **P - Inorganic ion transport** | -15 | +125 | ↑↑ | 後期に急活性化 |
| **G - Carbohydrate metabolism** | -11 | +32 | ↑ | 後期に活性化に転換 |
| **R - General function** | +220 | +419 | ↑ | 持続的活性化・強化 |
| **J - Translation** | -83 | -98 | ↓ | 抑制が持続・強化 |
| **Q - Secondary metabolism** | +9 | +12 | → | 一貫した活性化 |

### インサイト

> **Insight #4**: 成長-二次代謝トレードオフの明確な証拠。翻訳カテゴリ（J）はUp:Down = 1:6という顕著な抑制を示す一方、二次代謝カテゴリ（Q）は発現上昇のみ（Down = 0）で一方向的な誘導を示す。T2からT3にかけて、イオン輸送（P）と糖代謝（G）が抑制から活性化へ反転する代謝リモデリングが観察される。
>
> **Figure**: `report_COG_net_change.pdf`（T3vsT1）, `report_COG_net_change_T2vsT1.pdf`（T2vsT1）

- **細胞分裂（D）の抑制**（Up:Down = 1:10）は増殖停止を裏付ける
- **イオン輸送（P）の大幅活性化**（Net +125）は栄養枯渇への適応を反映
- **T2→T3でイオン輸送（P）が-15から+125へ反転**：後期の栄養取り込み戦略の変化

### 出力ファイル
| ファイル | 内容 |
|---------|------|
| `figures/COG_barplot_T3vsT1.pdf` | T3vsT1 機能カテゴリ別バーチャート |
| `figures/COG_barplot_T2vsT1.pdf` | T2vsT1 機能カテゴリ別バーチャート |
| `figures/COG_scatter_up_vs_down.pdf` | Up/Down比較散布図 |
| `figures/report_COG_net_change.pdf` | **COG Net変化水平棒グラフ T3vsT1（論文用・全19カテゴリ）** |
| `figures/report_COG_net_change_T2vsT1.pdf` | **COG Net変化水平棒グラフ T2vsT1（論文用・全20カテゴリ）** |
| `tables/COG_classification_summary.tsv` | カテゴリ別集計（全比較） |
| `tables/report_COG_T3vsT1.tsv` | T3vsT1 全19カテゴリ詳細 |
| `tables/report_COG_T2vsT1.tsv` | T2vsT1 全20カテゴリ詳細 |
| `tables/gene_COG_classification.tsv` | 全遺伝子のCOGカテゴリ |

---

## 5. BGCカバレッジトラック解析

### 4つの主要BGCクラスターの発現動態

| BGC | 領域 | T1 | T2 | T3 | パターン |
|-----|------|-----|-----|-----|----------|
| **act (Actinorhodin)** | 5,513,809-5,535,091 | 低 | 中 | 高 | 後期誘導型 |
| **red (Undecylprodigiosin)** | 6,432,812-6,464,206 | 低 | **高** | 高 | 中期から誘導 |
| **cda (CDA)** | 3,519,449-3,602,320 | 低 | 中 | 中 | 中期誘導・維持 |
| **cpk (Coelimycin P1)** | 6,900,898-6,948,414 | 低 | 中 | 高 | 後期誘導型 |

### インサイト

> **Insight #5**: BGCクラスター間で誘導タイミングに明確な差異がある。red（早期：T2から高発現）→ cda（中期誘導・維持）→ act/cpk（後期誘導型）の順序で活性化する。これは各BGCの制御機構の違い（レギュレーター、シグナル応答）を反映している。
>
> **Figure**: `report_BGC_expression_heatmap.pdf`（発現動態ヒートマップ）, `coverage_tracks_all_BGC.pdf`

- **すべてのBGCがT1で低発現**：増殖期には二次代謝が抑制されていることを確認（`coverage_track_*.pdf`）
- **red クラスターの早期誘導**はKEGGのProdigiosin biosynthesisの結果（T2で2.83倍）と一致

### 出力ファイル
| ファイル | 内容 |
|---------|------|
| `figures/coverage_track_act.pdf` | Actinorhodinクラスター |
| `figures/coverage_track_red.pdf` | Undecylprodigiosinクラスター |
| `figures/coverage_track_cda.pdf` | CDAクラスター |
| `figures/coverage_track_cpk.pdf` | Coelimycin P1クラスター |
| `figures/coverage_tracks_all_BGC.pdf` | 4BGC統合図 |
| `figures/report_BGC_expression_heatmap.pdf` | **BGC発現動態ヒートマップ（論文用）** |

---

## 総合的結論

### *S. coelicolor* M145 の培養段階依存的なトランスクリプトーム変動モデル

```
T1 (初期)           T2 (中期)            T3 (後期)
    │                   │                    │
    ▼                   ▼                    ▼
┌──────────┐       ┌──────────┐        ┌──────────┐
│  増殖期  │  →→→  │  移行期  │  →→→   │  定常期  │
│          │       │          │        │          │
│・リボソーム│      │・red誘導  │       │・act/cpk最大│
│  高発現  │       │・リボソーム│       │・リボソーム↓│
│・一次代謝│       │  最大抑制 │       │・二次代謝最大│
│・DNA複製 │       │・QS活性化 │       │・分裂停止   │
└──────────┘       └──────────┘        └──────────┘
     ↑                  ↑                    ↑
  venn図で         KEGGで最大の          カバレッジで
  T1基準確認       変動を確認            BGC発現確認
```

---

## 次のステップへの示唆

1. **BGC特異的な転写因子の同定**
   - 各BGCの誘導タイミングを制御するレギュレーターの特定
   - 関連解析: `07_regulator_network/`, `09_SARP_integration/`

2. **メチローム統合解析との連携**
   - `11_epigenome_integration/` の結果と組み合わせ、エピジェネティック制御の解明

3. **代謝フラックス解析**
   - 発現変動と実際の代謝産物量の対応確認

---

## 使用スクリプト

| スクリプト | 機能 |
|-----------|------|
| `01_venn_diagram_DEGs.R` | ベン図作成 |
| `02_GO_KEGG_enrichment_simple.R` | GOエンリッチメント |
| `02b_KEGG_enrichment.R` | KEGGエンリッチメント |
| `03_COG_classification.R` | COG機能分類 |
| `04_coverage_tracks.R` | カバレッジトラック |
| `05_report_tables_visualization.R` | **論文用グラフ作成（PDF/SVG）** |

---

## 論文用Figure一覧（Publication-quality）

共同研究先への共有および論文投稿用に作成した高品質グラフ（PDF/SVG形式）。

| Figure | ファイル名 | 対応するインサイト/主張 |
|--------|-----------|------------------------|
| **Fig. S1** | `report_key_insights_summary.pdf` | 全インサイトの俯瞰図 |
| **Fig. S2** | `report_pathway_pattern_summary.pdf` | 63パスウェイの変動パターン分類（Insight #6: 栄養枯渇応答） |
| **Fig. S3** | `report_FE_change_patterns.pdf` | FE変化の矢印図（Insight #1: T2ピーク, #2: Prodigiosin早期誘導） |
| **Fig. S4** | `report_BGC_expression_heatmap.pdf` | BGC発現動態ヒートマップ（Insight #5: BGC誘導タイミング） |
| **Fig. S5a** | `report_COG_net_change.pdf` | COG機能分類Net変化 T3vsT1（全19カテゴリ、Insight #4: 成長-二次代謝トレードオフ） |
| **Fig. S5b** | `report_COG_net_change_T2vsT1.pdf` | COG機能分類Net変化 T2vsT1（全20カテゴリ、T2時点での代謝状態） |
| **Fig. S6** | `report_both_down_FE_comparison.pdf` | 持続的Down-regulatedパスウェイ（T2で始まった抑制の不可逆性） |

**ファイル形式**: PDF（印刷用）, SVG（編集用）
**保存先**: `12_supplementary_figures/analysis/12_supplementary_260202_v1/figures/`
**COGカテゴリ数**: T3vsT1=19, T2vsT1=20（Nカテゴリのみ差異、NCBI COG 2024標準の26カテゴリのうち本データで検出された数）

---

*レポート作成: Claude Code*
*解析環境: conda rna-seq*
*最終更新: 2026-02-02*
