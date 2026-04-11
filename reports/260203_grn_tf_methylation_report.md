# GRN TF-Methylation Integration Report

**日付:** 2026-02-03
**プロジェクト:** *Streptomyces coelicolor* A3(2) M145 エピゲノム-トランスクリプトーム統合解析
**解析ディレクトリ:** `11_epigenome_integration/analysis/12_grn_tf_methylation/`

---

## Key Insights（主要なインサイト一覧）

| # | Insight | 対応Figure/Table | 根拠となる数値 |
|---|---------|-----------------|---------------|
| 1 | **redZはGRN内で唯一メチル化-発現協調変動を示すTF** | `tf_hierarchy_methylation_overview.png`, `TF_coordinated_changes.csv` | 6mA Lost + log2FC=-2.25, p=1.28e-17 |
| 2 | **SC_RS17645 (N-6 MTase) 発現低下がredZメチル化と連動** | `epigenetic_cascade_diagram.png` | MTase log2FC=-2.19 (T2vsT1) |
| 3 | **redDパラドックス: redZ低下にもかかわらず高発現** | `triad_network_visualization.png` | redZ: -2.25, redD: +4.77 |
| 4 | **Red BGC遺伝子は9/10が上昇発現** | `bgc_expression_summary.png` | Mean log2FC=2.85 |
| 5 | **GRN TFの大半（34/37）は非メチル化** | `tf_methylation_summary.csv` | Tier1: 3/14, Tier2: 0/17, Tier3: 1/6 メチル化 |

---

## 1. 解析の背景と目的

### 1.1 背景

前解析において、M145では以下の知見が得られている：

1. **新規メチル化モチーフ AAGCCCG (6mA)** - REBASEに未登録
2. **SC_RS17645 (N-6 DNA methylase)** - AAGCCCGメチル化の候補酵素、T2で発現低下
3. **4mC-発現相関 (r=0.137)** - T2時点で有意な正の相関

本解析では、**遺伝子制御ネットワーク (GRN)** の階層構造（Global → Pleiotropic → CSR）において、メチル化がどのように転写因子活性に影響し、最終的に二次代謝産物（特にAct, Red）生合成に寄与するかを統合的に解析する。

### 1.2 目的

1. 文献に基づくM145 GRN階層構造の構築（37 TFs）
2. GRN TFのプロモーターメチル化状態解析
3. TF-メチル化-標的遺伝子のtriad関係解析
4. エピジェネティック制御ネットワークの可視化

---

## 2. 方法

### 2.1 GRN階層構造の構築

文献（Zorro-Aranda et al., 2022等）およびデータベース（Abasy Atlas, DBSCR）に基づき、以下の階層でTFを分類：

| Tier | カテゴリ | 遺伝子数 | 代表的TF |
|------|---------|---------|---------|
| 1 | Global regulators | 14 | bldA, bldD, adpA, afsR |
| 2 | Pleiotropic regulators | 8 | absA1/2, afsQ1/2, wblA |
| 2 | Sigma factors | 9 | hrdB, hrdD, sigR, sigF |
| 3 | CSR (Cluster-Situated Regulators) | 6 | actII-ORF4, redD, redZ |

### 2.2 メチル化-発現統合

既存のメチル化データ（`integrated_methyl_expression_weighted.csv`）とGRN TFを照合し、各TFのメチル化状態を分類：
- **Gained in T2**: T2でメチル化増加
- **Lost in T2**: T2でメチル化減少
- **Stable**: 変化なし
- **Unmethylated**: メチル化なし

### 2.3 Triad解析

既知の制御関係（REGULATORY_NETWORK）を用いて、TF→標的遺伝子の関係を追跡し、メチル化-発現-標的の三者関係を解析。

---

## 3. 結果

### 3.1 GRN TFのメチル化状態

**37 TFのうち4 TFのみがメチル化サイトを保持**

| TF | Tier | T1 sites | T2 sites | Status | log2FC | BGC |
|----|------|----------|----------|--------|--------|-----|
| **redZ** | 3 | 2 | 2 | **協調変動** | -2.25 | red |
| bldN | 1 | 2 | 1 | Lost in T2 | +0.50 | - |
| afsR | 1 | 1 | 1 | Stable | -0.26 | - |
| afsS | 1 | 1 | 1 | Stable | -0.03 | - |

> **Insight #1**: redZはGRN内で唯一メチル化-発現協調変動を示すTF
> TF_coordinated_changes.csvにおいて、redZ (SC_RS27300) は6mAメチル化消失と発現低下が協調している唯一のGRN TFである。これは**Red BGCへのエピジェネティック制御の直接的証拠**となる。
> **根拠**: 6mA Lost (56.0 → 0 freq), log2FC=-2.25, padj=1.28e-17

### 3.2 redZ → redD → Red BGC カスケード

```
SC_RS17645 (N-6 MTase)
    ↓ 発現低下 (log2FC=-2.19)
AAGCCCG メチル化減少
    ↓
redZ プロモーター脱メチル化
    ↓
redZ 発現低下 (log2FC=-2.25)
    ↓
redD → Red BGC ???
```

> **Insight #2**: SC_RS17645 (N-6 MTase) 発現低下がredZメチル化と連動
> メチラーゼ活性低下→メチル化レベル低下という分子機構を支持する。SC_RS17645はAAGCCCGモチーフをメチル化する候補酵素である。
> **根拠**: MTase log2FC=-2.19, padj=6.48e-16

### 3.3 redDパラドックス

**問題**: redZはredDの正の制御因子であるが、redZ発現低下にもかかわらずredDは強く上昇発現している。

| 遺伝子 | 役割 | log2FC | 予測 | 実測 |
|--------|-----|--------|------|------|
| redZ | redD活性化因子 | -2.25 | redD↓ | redD↑ |
| absA2 | redD抑制因子 | +6.29 | redD↓ | redD↑ |
| **redD** | Red活性化因子 | **+4.77** | - | **強い上昇** |

> **Insight #3**: redDパラドックス - redZ低下にもかかわらず高発現
> AbsA2は抑制因子であり、その発現上昇はredDを抑制すべきだが、redDは上昇している。これは以下を示唆する：
> 1. 他の活性化因子（AfsR, BldD, papR2等）の寄与
> 2. AbsA2活性の飽和または翻訳後制御
> 3. 制御ネットワークのフィードバック構造
> **根拠**: absA2 log2FC=+6.29, redD log2FC=+4.77

### 3.4 Red BGC発現

> **Insight #4**: Red BGC遺伝子は9/10が上昇発現
> Red BGC領域の10遺伝子を解析した結果、9遺伝子がT2で上昇発現（log2FC > 1）を示した。平均log2FC=2.85であり、**Red経路は強く活性化**されている。
> **根拠**: n_up=9, n_down=0, mean=2.85

---

## 4. Figure解説

### 4.1 tf_hierarchy_methylation_overview.png

**ファイル**: `analysis/12_grn_tf_methylation/tf_promoter_methylation/tf_hierarchy_methylation_overview.png`

**説明**: GRN階層構造における転写因子のメチル化と発現を4パネルで可視化。

**パネル構成**:
- 左上: Tier別メチル化サイト数（T1/T2/T3）
- 右上: メチル化状態の分布（pie chart）
- 左下: TF-BGC相関散布図（Act vs Red）
- 右下: 主要TFの発現変動（bar chart）

**読み取り方**:
- Tier 1（Global）にメチル化TFが集中
- 大半のTFは非メチル化
- actII-ORF4はAct BGCと強い正相関（r=0.95）

**サポートするInsight**: #1, #5

### 4.2 regulatory_cascade_methylation.png

**ファイル**: `analysis/12_grn_tf_methylation/tf_promoter_methylation/regulatory_cascade_methylation.png`

**説明**: GRN制御カスケードにメチル化状態をオーバーレイした図。ノード色がメチル化状態、発現方向は矢印で表示。

**読み取り方**:
- オレンジ: メチル化増加（Gained）
- 青: メチル化減少（Lost）
- 灰色: 非メチル化
- 矢印: 制御関係（→活性化、⊣抑制）

**サポートするInsight**: #1, #2

### 4.3 triad_network_visualization.png

**ファイル**: `analysis/12_grn_tf_methylation/triad_analysis/triad_network_visualization.png`

**説明**: 左パネル: TF-Methylation-Target制御ネットワーク。右パネル: エピジェネティックカスケード詳細図。

**読み取り方**:
- 左: redZのメチル化消失（赤ノード）がネットワーク内で目立つ
- 右: MTase→Methylation→TF→Targetの線形カスケード
- パラドックス（redD上昇）の位置が明示

**サポートするInsight**: #2, #3

### 4.4 epigenetic_cascade_diagram.png

**ファイル**: `analysis/12_grn_tf_methylation/network_visualization/epigenetic_cascade_diagram.png`

**説明**: エピジェネティックカスケードの詳細フローチャート。各ステップに根拠データを付記。

**サポートするInsight**: #2, #3

### 4.5 comprehensive_epigenetic_network.png

**ファイル**: `analysis/12_grn_tf_methylation/network_visualization/comprehensive_epigenetic_network.png`

**説明**: MTase→Global TF→Pleiotropic→CSR→BGCの完全な制御ネットワーク。論文用Figure候補。

**サポートするInsight**: #1, #2, #3, #4, #5

### 4.6 bgc_expression_summary.png

**ファイル**: `analysis/12_grn_tf_methylation/triad_analysis/bgc_expression_summary.png`

**説明**: Act, Red, CDA, Cpk各BGCの遺伝子発現変動を一覧表示。

**サポートするInsight**: #4

---

## 5. データファイル一覧

| ファイル | 形式 | 説明 | サポートするInsight |
|---------|------|------|-------------------|
| `hierarchical_network_tfs.csv` | CSV | 37 TFの階層構造データ | #5 |
| `regulatory_interactions.csv` | CSV | 19の制御相互作用 | #3 |
| `tf_methylation_summary.csv` | CSV | TFメチル化状態サマリー | #1, #5 |
| `bgc_regulator_methylation.csv` | CSV | BGC制御因子のメチル化詳細 | #1 |
| `tf_target_triads.csv` | CSV | 23のTF-標的triad | #3 |
| `key_insights_triad.csv` | CSV | Triad解析からの主要知見 | #1, #2, #3 |

---

## 6. 結論

### 6.1 主要な発見

1. **redZはM145 GRNにおいて唯一のエピジェネティック制御ターゲット**
   - 6mAメチル化消失と発現低下が協調
   - Red BGCへの直接的影響を示唆

2. **SC_RS17645 (N-6 MTase) → AAGCCCG → redZの分子経路**
   - メチラーゼ発現低下→メチル化減少→TF発現変動
   - 新規R-M系によるエピジェネティック制御

3. **制御ネットワークの複雑性**
   - redDパラドックス: 単純な線形モデルでは説明不可
   - 複数の活性化因子・抑制因子の相互作用

### 6.2 二次代謝への示唆

- **Red BGC**: 強く活性化（9/10遺伝子上昇、redD +4.77）
- エピジェネティック制御が**redZ→redD経路**に影響
- しかし、最終的なredD活性化は**他の制御機構**が支配的

### 6.3 限界と今後の課題

1. **redDパラドックスの解明**: 翻訳後修飾、タンパク質相互作用の検討
2. **他のBGC**: Act, CDA, CpkでのMTase-TF経路探索
3. **実験的検証**: SC_RS17645ノックアウト株でのメチル化・発現解析

---

## 7. 次のステップ

| 優先度 | タスク | 目的 |
|--------|--------|------|
| **高** | 論文用Figure作成 | Step B実行 |
| **高** | redDパラドックス解析 | 制御機構の詳細解明 |
| **中** | ChIP-seq統合 | TF結合部位とメチル化の関係 |
| **低** | 他株比較 | M145特異的エピジェネティクス |

---

## 8. 参考文献

- Zorro-Aranda et al. (2022) Nucleic Acids Research - M145 GRN (5386 genes, 9707 interactions)
- Abasy Atlas - https://abasy.ccg.unam.mx/
- DBSCR - Streptomyces coelicolor regulatory database

---

*生成日: 2026-02-03*
*解析スクリプト:*
- `scripts/build_grn_hierarchy.py`
- `scripts/tf_promoter_methylation_analysis.py`
- `scripts/tf_methylation_target_triad.py`
- `scripts/epigenetic_network_visualization.py`
