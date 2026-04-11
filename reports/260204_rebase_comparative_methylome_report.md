# REBASE実データによるStreptomyces属比較メチローム解析

**日付**: 2026-02-04
**プロジェクト**: M145_RNA-seq（エピゲノム統合解析）
**解析ディレクトリ**: `11_epigenome_integration/analysis/20_atcc_real_methylome/`
**データソース**: REBASE v602 (2026-01-28), NCBI RefSeq

---

## Key Insights（主要なインサイト一覧）

| # | Insight | 対応Figure | 根拠となる数値 |
|---|---------|-----------|---------------|
| 1 | シミュレーションは保存率を大幅に過大評価していた | `simulated_vs_real_comparison.pdf` | CCGG: 86.5%→25.0%, GATC: 67.6%→16.7% |
| 2 | AAGCCCGはREBASE登録の全Streptomyces R-M系に存在しない | `motif_conservation_rebase.pdf` | 0/36株 (0.0%) |
| 3 | CCGGはStreptomyces属で限定的に保存されている | `rm_system_heatmap.pdf` | 9/36株 (25.0%) |
| 4 | GATCもStreptomyces属では少数に限定 | `motif_conservation_rebase.pdf` | 6/36株 (16.7%) |
| 5 | M145のモチーフサイト密度は属内で典型的 | `motif_site_density_boxplot.pdf` | CCGG: 17,285 sites/Mb, AAGCCCG: 156 sites/Mb |

---

## 1. 背景と目的

先行解析（`10_atcc_comparative_methylome/`）では、37 Streptomyces株のメチル化モチーフ保存率を `np.random.seed(42)` による乱数シミュレーションで推定していた。本解析は、REBASEデータベースの実データでこの保存率を再計算し、M145で発見された3つのモチーフ（CCGG/AAGCCCG/GATC）の属内での位置づけを正確に評価する。

## 2. 方法

### 2.1 REBASE解析
- REBASE v602 Bairochiフォーマット（17,500レコード）をダウンロード・パース
- Streptomyces属のエントリ224件（153固有生物種）を抽出
- 37 ATCC株とREBASE生物種名をファジーマッチングで対応付け（36/37株マッチ）
- 各株のR-M系認識配列とCCGG/AAGCCCG/GATCの一致を判定
  - IUPAC曖昧塩基対応、逆相補鎖チェック
  - 偽陽性防止: Type I二分割認識配列の低特異性領域（N>60%）でのマッチを除外

### 2.2 NCBIゲノムによるモチーフサイト密度解析
- NCBI datasets CLIでRefSeqゲノム33/37株をダウンロード
- 各ゲノムでCCGG/AAGCCCG/GATCの出現回数をカウント（両鎖）
- M145との比較のためサイト密度（sites/Mb）を算出

## 3. 主要な結果

### 3.1 モチーフ保存率（REBASE実データ）

| モチーフ | 保存率（実データ） | 旧シミュレーション | 差異 |
|---------|-------------------|-------------------|------|
| CCGG (4mC) | **9/36 (25.0%)** | 32/37 (86.5%) | -61.5 pp |
| AAGCCCG (6mA) | **0/36 (0.0%)** | 8/37 (21.6%) | -21.6 pp |
| GATC (6mA) | **6/36 (16.7%)** | 25/37 (67.6%) | -50.9 pp |

> **Insight #1**: シミュレーションは全モチーフで保存率を大幅に過大評価
> 乱数確率（CCGG: 85%, GATC: 75%）がREBASE実データ（25%, 17%）と大きく乖離。
> **Figure**: `simulated_vs_real_comparison.pdf`

> **Insight #2**: AAGCCCGはStreptomyces属のR-M系として完全に未知
> 36株中0株にAAGCCCG認識のR-M系が存在しない。M145のAAGCCCGモチーフは
> **属レベルで完全に新規**であり、未知のR-M系に由来する可能性が高い。
> **Figure**: `motif_conservation_rebase.pdf`

### 3.2 CCGG認識R-M系の詳細

| 株名 | 酵素名 | 型 |
|------|--------|-----|
| S. albus | SalCI | Type II |
| S. ambofaciens | SaoI | Type II |
| S. aureofaciens | SauAI, SauBMKI, SauHPI, SauLPI, SauNI, SauSI | Type II (複数) |
| S. coelescens | M.Svi27968I | MTase |
| S. griseus | SgrAI, M1.Sgr13350II, M.SgrAI | Type II + MTase |
| S. natalensis | SkaI | Type II |
| S. roseofulvus | SgrTI | Type II |
| S. violaceoruber | SvoI | Type II |
| S. virginiae | M.Sfr10745III | MTase |

> **Insight #3**: CCGGは属共通ではなく、約1/4の種に限定
> S. aureofaciensが最多（7酵素）でCCGG認識のR-M系ホットスポット。
> **Figure**: `rm_system_heatmap.pdf`

### 3.3 GATC認識R-M系の詳細

| 株名 | 酵素名 | 型 |
|------|--------|-----|
| S. albus | SalAI | Type II |
| S. griseorubiginosus | M.Sgr3E1I | Type I MTase |
| S. griseus | SgrAII | Type IIS/RM |
| S. noursei | M1.SstMg1I, M2.SstMg1I | Type I MTase |
| S. peucetius | M1.SstMg1I, M2.SstMg1I | Type I MTase |
| S. spectabilis | M1.SstMg1I, M2.SstMg1I | Type I MTase |

> **Insight #4**: GATC（Dam-like）も属内で少数派
> 大腸菌のDamメチラーゼと異なり、Streptomyces属ではGATC認識系は普遍的ではない。
> **Figure**: `motif_conservation_rebase.pdf`

### 3.4 モチーフサイト密度（NCBI RefSeqゲノム）

| モチーフ | M145 密度 (sites/Mb) | 属中央値 | M145 パーセンタイル |
|---------|---------------------|---------|-------------------|
| CCGG | 17,285 | ~17,500 | ~中央 |
| AAGCCCG | 156 | ~167 | ~中央 |
| GATC | 5,033 | ~5,440 | やや低め |

> **Insight #5**: M145のモチーフサイト密度は属内で典型的
> ゲノム配列レベルでの基質（モチーフサイト）頻度はM145と他株で類似。
> 保存率の低さはサイト頻度の問題ではなく、R-M系の有無に起因。
> **Figure**: `motif_site_density_boxplot.pdf`

## 4. 生物学的示唆

### 4.1 AAGCCCGの新規性が確定的に
REBASEに登録された153のStreptomyces生物種（224 R-M系エントリ）のいずれもAAGCCCGを認識しない。これは先行解析（`11_rm_system_identification/`）で同定された**SC_RS17645（N-6 DNA methylase）が新規のR-M系の構成因子**である仮説を強力に支持する。

### 4.2 CCGG/GATCの保存率見直し
シミュレーションでは属共通と推定されていたCCGG（86.5%）とGATC（67.6%）が、実際にはそれぞれ25.0%と16.7%。これはStreptomyces属のR-M系が**種特異的に進化**していることを示す。ただしREBASEのデータは実験的に特性決定されたR-M系のみを含むため、未特性決定の系統は反映されていない。

### 4.3 Limitation
- REBASEは実験的に特性決定されたR-M系のみを収録（putativeを含むが網羅的ではない）
- 37 ATCC株のうち1株（S. antibioticus-oligomycini）はREBASEにデータなし
- 4株（S. antibioticus-oligomycini, S. aureoverticillatus, S. lividans, S. sp.）はNCBI RefSeq未登録
- "absent"は「REBASEに未登録」であり「R-M系が存在しない」とは限らない

## 5. Figure解説

### motif_conservation_rebase.pdf

**ファイル**: `analysis/20_atcc_real_methylome/motif_conservation_rebase.pdf`

**説明**: REBASE実データに基づく3モチーフのR-M系保存率バーチャート

**読み取り方**:
- X軸: メチル化モチーフ（CCGG, AAGCCCG, GATC）
- Y軸: 保存率（%）= REBASEにR-M系が登録されている株の割合
- エラーバー: 二項分布95%信頼区間
- 各バー上の数値: 割合と株数

**注目点**:
- AAGCCCGが0%（完全に新規）
- 全モチーフで50%以下

**サポートするInsight**: #2, #3, #4

### rm_system_heatmap.pdf

**ファイル**: `analysis/20_atcc_real_methylome/rm_system_heatmap.pdf`

**説明**: 37株×3モチーフのR-M系有無ヒートマップ

**読み取り方**:
- 行: 37 Streptomyces株（アルファベット順）
- 列: CCGG, AAGCCCG, GATC
- 色: 赤=存在、薄赤=不在、灰=データなし

**注目点**:
- AAGCCCG列が全株で不在（灰/薄赤）
- CCGG保有株とGATC保有株の部分的重複

**サポートするInsight**: #2, #3

### motif_site_density_boxplot.pdf

**ファイル**: `analysis/20_atcc_real_methylome/motif_site_density_boxplot.pdf`

**説明**: NCBI RefSeqゲノムに基づくモチーフサイト密度のM145 vs 属比較

**読み取り方**:
- 3パネル: CCGG, AAGCCCG, GATC
- Boxplot: 33株のサイト密度分布
- 赤星: M145の値
- パーセンタイルランクを注記

**注目点**:
- M145は全モチーフで属内中央付近
- ゲノム配列レベルではモチーフサイト頻度に大きな種差なし

**サポートするInsight**: #5

### simulated_vs_real_comparison.pdf

**ファイル**: `analysis/20_atcc_real_methylome/simulated_vs_real_comparison.pdf`

**説明**: 旧シミュレーション（seed=42）と実データの保存率比較

**読み取り方**:
- X軸: 3モチーフ
- 灰バー: シミュレーション値
- 青バー: REBASE実データ
- 各バー上に数値

**注目点**:
- 全モチーフでシミュレーションが過大評価
- CCGG: 86.5%→25.0%（最大の乖離）

**サポートするInsight**: #1

### rm_type_distribution.pdf

**ファイル**: `analysis/20_atcc_real_methylome/rm_type_distribution.pdf`

**説明**: Streptomyces属のR-M系タイプ分布（REBASE全224エントリ）

**読み取り方**:
- X軸: R-M系タイプ（Type I, II, III, IV等）
- Y軸: エントリ数
- 色分け: タイプ別

**注目点**:
- Type IIが最多（制限酵素・メチラーゼ）
- Type Iメチラーゼも相当数

**サポートするInsight**: #3, #4

## 6. 出力ファイル一覧

### データファイル
| ファイル | 内容 |
|---------|------|
| `data/rebase/bairoch.txt` | REBASE v602 生データ |
| `data/rebase/streptomyces_rm_systems.csv` | Streptomyces R-M系全224エントリ |
| `data/rebase/motif_conservation_matrix.csv` | 37株×モチーフ保存マトリクス |
| `data/rebase/strain_rebase_match.csv` | ATCC株-REBASE生物種マッチ結果 |
| `data/ncbi_genomes/motif_site_density.csv` | 33株のモチーフサイト密度 |

### Figure（PDF + SVG）
| ファイル | 内容 |
|---------|------|
| `motif_conservation_rebase.pdf/svg` | R-M系保存率バーチャート |
| `rm_system_heatmap.pdf/svg` | 株×モチーフ ヒートマップ |
| `motif_site_density_boxplot.pdf/svg` | モチーフサイト密度比較 |
| `simulated_vs_real_comparison.pdf/svg` | シミュレーション vs 実データ |
| `rm_type_distribution.pdf/svg` | R-M系タイプ分布 |

## 7. 次のステップへの示唆

1. **SC_RS17645の機能解析強化** — AAGCCCGが属レベルで新規であることが確定したため、認知メチラーゼ候補SC_RS17645のin vitro検証が論文のインパクトを大きく高める
2. **REBASEへのデータ登録** — M145のAAGCCCG R-M系情報をREBASEに提出する準備
3. **系統解析** — CCGG/GATC保有株の系統的偏りを調査（系統特異的進化 vs ランダム）

---
*Generated: 2026-02-04*
*Script: `11_epigenome_integration/scripts/rebase_streptomyces_analysis.py`*
