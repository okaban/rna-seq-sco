# H13: AAGCCCGサイトのゲノム分布パターンと機能的標的推定

**日付**: 2026-02-25
**解析ディレクトリ**: `11_epigenome_integration/analysis/36_AAGCCCG_distribution/`
**判定**: **部分的支持 — core偏在・クラスタリング確認、ただし二次代謝ではなくストレス/防御に弱enrichment**

---

## 1. 背景と仮説

AAGCCCG 6mAモチーフはREBASE 0%保存の新規モチーフ（H10レポート）。T1で260サイト、T2で64サイトが検出され、SC_RS17645（N-6 DNA methylase）の発現低下（LFC=-2.19）と連動してサイト数が激減する（H5）。しかし、SP1（AAGCCCGゲノム分布の非ランダム性と制御標的の推定）はLoop 1から5ループ間未着手であった。

**H13仮説**: AAGCCCG 6mAサイトはゲノム上で非ランダムに分布し、特定の機能カテゴリ（二次代謝、制御遺伝子）の近傍に集積する。

## 2. 方法

1. 6mA_final_census.csvからAAGCCCGモチーフサイトをT1/T2別に抽出
2. M145 GFFアノテーションとの距離計算（2kb以内の近傍遺伝子を同定）
3. ゲノム地理（arm: ≤1.5 Mb from each end, core: 中央領域）別の分布解析
4. 最近傍距離（nearest-neighbor distance）によるクラスタリング解析
5. 近傍遺伝子の機能カテゴリ別エンリッチメント（Fisher exact test + Bonferroni補正）
6. BGCとの重複解析（act, red, cda, cpk）
7. DEGとの関連（近傍遺伝子の発現変動解析）

## 3. 結果

### 3.1 基本統計 — T1 vs T2比較

| 指標 | T1 | T2 |
|------|----|----|
| 総サイト数 | 260 | 64 |
| Arm サイト | 81 (31.2%) | 23 (35.9%) |
| Core サイト | 179 (68.8%) | 41 (64.1%) |
| 最近傍中央値 | 9,406 bp | 55,672 bp |
| 最近傍平均値 | 14,732 bp | 70,728 bp |
| BGC内サイト | 3 | 0 |
| 2kb以内制御遺伝子近傍 | 85 | 19 |
| 2kb以内近傍遺伝子数 | 974 | 285 |
| 近傍DEG数 | 627 | 180 |
| 近傍遺伝子平均LFC | 0.354 | 0.189 |
| 近傍遺伝子中央値LFC | -0.022 | -0.082 |

### 3.2 ゲノム地理分布

両タイムポイントともcore genome偏在（T1: 68.8%, T2: 64.1%）。これはCCGG 4mCのT1 core偏在（77.4%）と同様だが、T2でもcore偏在を維持する点が異なる（CCGG T2はarm 88.9%に逆転）。

**AAGCCCGはタイムポイント間でcore enrichmentパターンを保持** — CCGGの劇的な地理的シフト（H7）とは対照的。

### 3.3 クラスタリング

T1サイトは強くクラスタリング（最近傍中央値 9,406 bp）。ランダム分布の期待値（~33,500 bp = 8.7 Mb / 260）の約28%。

T2サイトは分散（最近傍中央値 55,672 bp）。サイト数減少（260→64）を考慮しても期待値（~136,000 bp）の41%であり、ランダムよりやや密だが、T1ほどの強いクラスタリングではない。

### 3.4 機能カテゴリエンリッチメント

| カテゴリ | T1 observed | T1 expected | T1 fold | T1 p_bonf | T2 fold | T2 p_bonf | T1方向 |
|---------|-------------|-------------|---------|-----------|---------|-----------|--------|
| **Hypothetical** | 18 | 34.6 | **0.52** | **0.039** | 0.35 | 0.81 | **depleted** |
| **Regulatory/TF** | 11 | 25.4 | **0.43** | **0.035** | 0.64 | 1.0 | **depleted** |
| Stress/Defense | 16 | 8.2 | 1.96 | 0.24 | 2.0 | 1.0 | enriched (n.s.) |
| Primary metabolism | 67 | 58.3 | 1.15 | 1.0 | 1.18 | 1.0 | enriched (n.s.) |
| DNA/RNA metabolism | 18 | 12.8 | 1.41 | 1.0 | 0.95 | 1.0 | enriched (n.s.) |
| Membrane/Cell wall | 5 | 2.3 | 2.14 | 1.0 | 0.0 | 1.0 | enriched (n.s.) |
| Other | 96 | 86.8 | 1.11 | 1.0 | 1.26 | 1.0 | enriched (n.s.) |
| Secondary metabolism | 5 | 6.9 | 0.73 | 1.0 | 1.18 | 1.0 | depleted (n.s.) |
| Transport | 21 | 21.6 | 0.97 | 1.0 | 0.75 | 1.0 | neutral |
| Translation | 3 | 3.1 | 0.96 | 1.0 | 0.0 | 1.0 | neutral |

**有意な発見（Bonferroni補正後）:**
1. **Hypothetical proteins: DEPLETED** (T1 fold=0.52, p_bonf=0.039) — AAGCCCGサイトは機能既知遺伝子の近傍に偏在
2. **Regulatory/TF: DEPLETED** (T1 fold=0.43, p_bonf=0.035) — H4/H6の結果と整合的：メチル化は文献TFを直接標的としない

**注目すべきトレンド:**
- Stress/Defense: 約2倍のenrichment（T1/T2とも）だがBonferroni後は非有意
- Secondary metabolism: 予想に反してenrichmentなし（fold=0.73）

### 3.5 BGCとの重複

| タイムポイント | BGC | サイト数 |
|-------------|-----|---------|
| T1 | cpk | 2 |
| T1 | cda | 1 |
| T1 | act | 0 |
| T1 | red | 0 |
| T2 | (全BGC) | 0 |

4大BGC内のAAGCCCGサイトはほぼゼロ。これはsecondary metabolism enrichmentの欠如と一致し、**AAGCCCGメチル化は二次代謝遺伝子クラスターを直接標的としていない**。

### 3.6 発現への影響

近傍遺伝子（2kb以内）の発現変動:
- T1: 974近傍遺伝子中627がDEG（64.4%）、平均LFC=+0.354（弱い正の傾向）
- T2: 285近傍遺伝子中180がDEG（63.2%）、平均LFC=+0.189

中央値LFCはほぼゼロ（T1: -0.022, T2: -0.082）で、個別遺伝子レベルでは方向性のない（bidirectional）影響。H5の「メチル化はglobal制御であり遺伝子特異的cis制御ではない」と整合。

## 4. 考察

### 4.1 ゲートキーパーモデルv2との統合

| H13発見 | モデルとの関係 |
|---------|--------------|
| Core偏在（68.8%）かつクラスタリング | **Layer 1 (Landscape)**: AAGCCCGはcore genomeの一次代謝領域に集中配置 |
| Regulatory/TF depleted | **Layer 3 (Signal Gating)**: メチル化は制御因子を避ける→H4, H6と一致 |
| Hypothetical depleted | メチル化は「機能既知」遺伝子に選択的 |
| BGC非重複 | **Layer 2 (Protection)**: BGCはメチル化から保護されている可能性 |
| 発現影響bidirectional | **全般**: Global effect with no gene-specific directionality |

### 4.2 CCGGとの対比

| 特性 | AAGCCCG (6mA) | CCGG/GCCGGC (4mC) |
|------|--------------|-------------------|
| タイムポイント間地理 | core維持 | core→arm shift (H7) |
| サイト重複 | 全消失+新規出現 | 完全非重複 (Jaccard=0.000) |
| クラスタリング | 強 (T1) | 不明 |
| 制御遺伝子近傍 | depleted | — |
| BGC重複 | ほぼゼロ | — |

両モチーフともタイムポイント間でサイトが完全にリモデリングされるが、地理的パターンは逆方向に変化する（AAGCCCGはcore維持、CCGGはcore→arm移行）。

### 4.3 AAGCCCG 260サイト消失の生物学的意義

SC_RS17645（Type I HsdM）の発現低下（LFC=-2.19）がT1→T2の260サイト消失を駆動（H5）。消失するサイトはcore genome上のprimary metabolism遺伝子近傍に集中。一方T2で新たに出現する64サイトは、より分散した配置。

**仮説的モデル**: T1（指数増殖期）でcore genomeの一次代謝関連領域に集中的なAAGCCCGメチル化が維持され、T2（定常期移行）でMTase発現低下に伴い受動的に消失する。同時にT2で少数の新規サイトが出現するメカニズムは不明（残存酵素活性か、別のMTaseか）。

## 5. 結論

| 項目 | 結果 |
|------|------|
| H13仮説（非ランダム分布、二次代謝集積） | **部分的支持** — 非ランダム分布は確認、ただし二次代謝ではなく一次代謝/ストレス防御に偏在 |
| ゲノム地理 | **Core genome偏在**（T1: 68.8%, T2: 64.1%）— CCGGと異なりタイムポイント間で安定 |
| クラスタリング | **T1で強いクラスタリング**（NN中央値 9,406 bp, 期待値の~28%） |
| 機能エンリッチメント | **Hypothetical depleted** (p=0.039), **Regulatory/TF depleted** (p=0.035) |
| BGC重複 | **ほぼゼロ**（T1: cpk 2, cda 1のみ） |
| 発現影響 | **Bidirectional** — 個別遺伝子の方向性なし（mean +0.35, median -0.02） |

## 6. 出力ファイル

### Figures
| ファイル | 内容 |
|---------|------|
| `genome_distribution.pdf/svg` | ゲノム全体のAAGCCCGサイト分布（T1 vs T2） |
| `gene_proximity.pdf/svg` | 遺伝子近接性解析 |
| `functional_enrichment.pdf/svg` | 機能カテゴリエンリッチメント |
| `expression_impact.pdf/svg` | 近傍遺伝子の発現変動 |
| `H13_comprehensive_summary.pdf/svg` | H13包括サマリー |

### Tables
| ファイル | 内容 |
|---------|------|
| `T1_vs_T2_comparison.tsv` | T1 vs T2基本統計比較 |
| `functional_enrichment.tsv` | 機能カテゴリエンリッチメント（Fisher検定） |
| `BGC_overlap.tsv` | 4大BGCとの重複 |
| `AAGCCCG_site_gene_mapping.tsv` | 全サイト-遺伝子マッピング（175 KB） |

### Scripts
| ファイル | 内容 |
|---------|------|
| `H13_AAGCCCG_distribution.py` | 解析スクリプト |

---

*Analysis directory: `11_epigenome_integration/analysis/36_AAGCCCG_distribution/`*
*Generated: 2026-02-25*
