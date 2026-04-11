# AAGCCCG メチル化カスケードによるBGC制御因子のエピゲノム制御解析レポート

**日付:** 2026-02-04
**プロジェクト:** *Streptomyces coelicolor* A3(2) M145 エピゲノム-トランスクリプトーム統合解析
**解析ディレクトリ:** `11_epigenome_integration/analysis/19_bgc_regulator_overview/`
**前提レポート:** `260204_bgc_regulator_methylation_landscape_report.md`

---

## Key Insights（主要なインサイト一覧）

| # | Insight | 対応Figure/Table | 根拠となる数値 |
|---|---------|-----------------|---------------|
| 1 | **afsS/redZのメチル化部位はいずれもAAGCCCGモチーフ上に位置し、プロモーターの機能的領域内にある** | `aagcccg_cascade_analysis.png` Panel A-1, A-2 | afsS: TSS上流172bp, redZ: TSS上流108bp |
| 2 | **SC_RS17645 (MTase) の発現低下（4.6倍）がafsS/redZのメチル化消失と完全に同期する** | `aagcccg_cascade_analysis.png` Panel B | MTase LFC=-2.19 (T2vT1); afsS 4mC: 83%→0%; redZ 6mA: 59%→0% |
| 3 | **T3でMTaseが部分回復してもメチル化は回復しない（受動的脱メチル化モデル）** | `aagcccg_cascade_analysis.png` Panel B | MTase T3: 145 counts (回復); 4mC/6mA: 0%のまま |
| 4 | **ゲノム全体でAAGCCCGの78.8%がメチル化されており、afsS/redZはその中からT2で選択的に脱メチル化された部位** | `aagcccg_cascade_analysis.png` Panel C, `aagcccg_cascade_summary.tsv` | 669/849サイトがメチル化; 163遺伝子で発現と協調変動 |
| 5 | **BGC制御因子のうちAAGCCCGプロモーターメチル化を持つのはafsS/redZのみ** | `aagcccg_cascade_analysis.png` Panel D | 12制御因子中2因子のみがAAGCCCG+メチル化+消失の全条件を満たす |
| 6 | **3段階カスケードモデルの提唱: SC_RS17645↓ → AAGCCCG脱メチル化 → afsS/redZ転写変動 → BGC発現制御** | `aagcccg_cascade_analysis.png` 全パネル | 時間的整合性・モチーフ一致・プロモーター位置の3条件が揃う |

---

## 1. 解析の背景と目的

### 1.1 背景

前解析（`260204_bgc_regulator_methylation_landscape_report.md`）において、27のBGC制御因子のうちafsS（4mC消失）、redZ（6mA消失）、bldN（6mA減少）の3因子のみが発現変動とメチル化変動の両方を示すことを報告した。

しかし、以下の重要な問題が未解決であった：

1. メチル化部位のプロモーター内の**正確な位置**とσ因子認識配列との空間的関係
2. メチル化部位が認識する**モチーフの同一性**（CCGGなのかAAGCCCGなのか）
3. 候補メチルトランスフェラーゼ SC_RS17645 の発現動態とメチル化消失の**時間的整合性**
4. ゲノム全体のAAGCCCGメチル化パターンの中での**afsS/redZの位置づけ**

### 1.2 目的

1. **A-1**: afsS/redZプロモーターのメチル化部位スキーマティック図を作成し、TSS・σ因子モチーフ・AAGCCCGモチーフ・メチル化部位の空間的関係を明示する
2. **A-2**: SC_RS17645（候補MTase）の発現動態とafsS/redZメチル化消失の時間的対応を図示・評価する
3. **A-3**: ゲノム全体のAAGCCCG含有プロモーターのメチル化パターンを解析し、afsS/redZの脱メチル化がゲノムワイドなプログラムの一部か特異的事象かを評価する

---

## 2. 方法

### 2.1 メチル化部位の座標特定

高信頼度メチル化サイトデータ（`/Users/okaban/bioinfo/methyl/260102_M145/analysis/cursor_results/20260108/high_confidence_sites.csv`）から、afsS（SC_RS22980）およびredZ（SC_RS27300）の遺伝子領域内・プロモーター領域のメチル化サイトを特定した。

### 2.2 プロモーター構造情報の統合

以下のデータを統合してスキーマティック図を作成した：
- **TSS座標**: `18_tss_analyses/comprehensive_tss_table.csv`（Jeong et al. 2016 dRNA-seq）
- **σ因子モチーフ**: `18_tss_analyses/D1_sigma_motifs.csv`（-10 box配列と位置）
- **AAGCCCGモチーフ空間分布**: `18_tss_analyses/motif_AAGCCCG_spatial_detail.csv`（is_methylatedフラグ含む）

### 2.3 MTase発現データ

SC_RS17645のDESeq2結果（3比較）および正規化カウント（9サンプル）を抽出し、メチル化頻度の時系列データと並列で図示した。

### 2.4 ゲノムワイドAAGCCCG解析

- `14_aagcccg_promoter_analysis/aagcccg_promoter_counts.csv`: 8,083プロモーターのAAGCCCGモチーフ有無
- `18_tss_analyses/motif_AAGCCCG_spatial_detail.csv`: 849サイトのメチル化状態（is_methylated）
- 上記を統合し、メチル化率・TSS近傍メチル化率・発現協調変動率を算出した

### 2.5 スクリプト

`11_epigenome_integration/scripts/plot_aagcccg_cascade_analysis.py`

---

## 3. 結果

### 3.1 A-1: afsS/redZ プロモーターのメチル化部位の空間的配置

#### afsS (SC_RS22980)

| 要素 | ゲノム座標 | TSSからの距離 | 備考 |
|------|-----------|-------------|------|
| TSS | 4,576,134 | 0 | Jeong2016 dRNA-seq |
| -10 box | ~4,576,146 | -12 bp | 配列: TAGACT |
| -35 box (推定) | ~4,576,163 | -29 bp | -10 boxから-17 bp |
| **4mC部位** | **4,576,306** | **-172 bp** | **AAGCCCGモチーフ上** |
| AAGCCCGモチーフ | 4,576,305 | -171 bp | is_methylated=True |

> **Insight #1（afsS）**: 4mC部位はTSSから172 bp上流のAAGCCCGモチーフ上に位置する。
> この距離はプロモーター領域の典型的な制御要素配置範囲内（-50〜-300 bp）であり、
> メチル化がプロモーター活性に影響し得る空間的条件を満たしている。
> **Figure**: `aagcccg_cascade_analysis.png` Panel A-1

**メチル化頻度の時間変動**:
- T1: 83.19%（coverage 24.0）→ T2: 0% → T3: 0%

#### redZ (SC_RS27300)

| 要素 | ゲノム座標 | TSSからの距離 | 備考 |
|------|-----------|-------------|------|
| TSS | 5,461,950 | 0 | GFFアノテーション |
| -10 box | ~5,461,962 | -12 bp | 配列: TAACGT |
| **6mA部位** | **5,462,058** | **-108 bp** | **AAGCCCGモチーフ上** |
| AAGCCCGモチーフ(1) | 5,462,053 | -103 bp | is_methylated=True |
| AAGCCCGモチーフ(2) | 5,462,367 | -417 bp | is_methylated=False |

> **Insight #1（redZ）**: 6mA部位はTSSから108 bp上流のAAGCCCGモチーフ上に位置する。
> 同一プロモーターに2つのAAGCCCGモチーフが存在するが、メチル化されていたのは
> TSS近位の1つのみであり、位置依存的なメチル化選択性が示唆される。
> **Figure**: `aagcccg_cascade_analysis.png` Panel A-2

**メチル化頻度の時間変動**:
- T1: 58.98%（coverage 46.3）→ T2: 0% → T3: 0%

#### 両プロモーターの共通点

- メチル化部位はいずれも**AAGCCCGモチーフ上**に位置
- TSS上流 **100-200 bp** の制御的に重要な領域に存在
- T1→T2で**完全消失**（段階的低下ではなくバイナリな変化）
- T3でも回復しない（不可逆的な脱メチル化）

### 3.2 A-2: SC_RS17645 (MTase) の発現動態とメチル化消失の時間的対応

#### SC_RS17645 発現データ

| 比較 | LFC | padj | 方向 |
|------|-----|------|------|
| T2 vs T1 | **-2.19** | 6.5e-16 | 強く低下（4.6倍） |
| T3 vs T1 | -0.71 | 3.3e-3 | 軽度低下 |
| T3 vs T2 | **+1.45** | 1.5e-7 | 回復 |

**正規化カウント**: T1: ~245 → T2: ~52 → T3: ~145

#### 時間的整合性の評価

```
タイムポイント   SC_RS17645 (MTase)     afsS 4mC      redZ 6mA
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    T1          245 counts (高発現)     83.2% ✓       59.0% ✓
    T2           52 counts (4.6x低下)    0.0% ✗        0.0% ✗
    T3          145 counts (部分回復)    0.0% ✗        0.0% ✗
```

> **Insight #2**: MTase発現の急激な低下（T1→T2, 4.6倍）と、両プロモーターの
> メチル化完全消失が同一タイムポイント間で発生する。これは「MTase発現低下 →
> DNA複製時のメチル化維持不能 → 受動的脱メチル化」モデルと整合する。
> **Figure**: `aagcccg_cascade_analysis.png` Panel B

> **Insight #3**: T3でMTase発現が部分的に回復（145 counts）しても、一旦消失した
> メチル化は回復しない。これは以下の理由で説明可能：
> (a) 新規合成DNAの半メチル化状態からの維持メチル化には十分なMTase量が必要
> (b) T2での細胞分裂による希釈効果が不可逆的
> (c) MTaseの回復が閾値に達していない
> **Figure**: `aagcccg_cascade_analysis.png` Panel B

### 3.3 A-3: ゲノムワイドAAGCCCGメチル化パターン

#### 全体統計

| 指標 | 値 |
|------|-----|
| 解析プロモーター総数 | 8,083 |
| AAGCCCG含有プロモーター | 416 (5.1%) |
| AAGCCCGサイト総数（TSS ±500bp） | 849 |
| メチル化AAGCCCGサイト | 669 (**78.8%**) |
| 非メチル化AAGCCCGサイト | 180 (21.2%) |
| TSS近傍（±200bp）のメチル化率 | 373/479 (**77.9%**) |
| 発現と協調変動する遺伝子 | 163/416 (**39.2%**) |

> **Insight #4**: AAGCCCGモチーフの約79%がメチル化を受けており、これはAAGCCCGが
> SC_RS17645メチルトランスフェラーゼの真の認識配列であることを強く支持する。
> 21%の非メチル化サイトは、配列コンテキスト依存的なメチル化選択性、または
> クロマチン構造による接近阻害を反映している可能性がある。
> **Figure**: `aagcccg_cascade_analysis.png` Panel C

#### BGC制御因子のAAGCCCGステータス

| 制御因子 | AAGCCCG in promoter | メチル化 | T2で消失 |
|----------|:---:|:---:|:---:|
| **afsS** | **Yes** | **Yes** | **Yes** |
| **redZ** | **Yes** | **Yes** | **Yes** |
| afsR | No | Yes (6mA, 別モチーフ) | No (stable) |
| bldN | No | Yes (6mA, 別モチーフ) | Yes |
| actII-ORF4 | No | No | — |
| redD | No | No | — |
| cdaR | No | No | — |
| cpkO/kasO | No | No | — |
| bldA | No | No | — |
| bldD | No | No | — |
| adpA | No | No | — |
| papR2 | No | No | — |

> **Insight #5**: 12の主要BGC制御因子のうち、AAGCCCGモチーフをプロモーターに持ち、
> かつメチル化を受け、かつT2で消失する全条件を満たすのは**afsSとredZの2因子のみ**。
> 両者はいずれも「外部シグナル → 二次代謝」のシグナル統合ノードに位置する点で
> 共通しており、AAGCCCGメチル化がカスケードの特定の階層を選択的に修飾する
> メカニズムを示唆する。
> **Figure**: `aagcccg_cascade_analysis.png` Panel D

---

## 4. 3段階カスケードモデル

本解析の結果を統合すると、以下のエピゲノム制御カスケードが提唱できる：

```
Step 1: SC_RS17645 (N-6 DNA methylase) 発現低下
        │  T1→T2で4.6倍低下 (LFC = -2.19, padj = 6.5e-16)
        │  原因: 培養段階の移行に伴う上流シグナル（未同定）
        ▼
Step 2: AAGCCCGモチーフでの脱メチル化（受動的メカニズム）
        │  afsS プロモーター: 4mC 83% → 0% (TSS -172bp)
        │  redZ プロモーター: 6mA 59% → 0% (TSS -108bp)
        │  ゲノム全体で163遺伝子に影響する協調的プログラム
        ▼
Step 3: 標的遺伝子の転写変動 → BGC制御への影響
        ├─ afsS 低下 (LFC = -1.39) → Act産生タイミング制御
        │   ※ afsS KOではAct完全消失 (Lee et al. 2002)
        │   → 4mC消失がafsS発現の「計画的低下」を許容
        └─ redZ 低下 (LFC = -2.25 T2vT1) → redDカスケードへの影響
            ※ bldA翻訳制御と組み合わさった多層的制御
            → 6mA消失がredZの転写制御に追加的な制御層を提供
```

### 4.1 モデルの根拠と限界

**支持する証拠:**
1. **モチーフの一致**: 両メチル化部位がAAGCCCGモチーフ上に位置（偶然の一致確率は低い）
2. **時間的整合性**: MTase低下とメチル化消失が同一タイムポイント間で発生
3. **プロモーター内位置**: TSS上流100-200 bpの制御的に重要な領域
4. **AAGCCCGの高メチル化率**: ゲノム全体で78.8%がメチル化 → 真の認識配列
5. **シグナル統合ノードへの限局**: BGC制御因子の中でafsS/redZのみが該当

**限界:**
1. **因果関係は未実証**: 相関データのみであり、メチル化消失が転写変動の原因であることは示していない
2. **サンプルサイズ**: 各タイムポイントn=3の生物学的複製
3. **時間分解能**: 3タイムポイントでは、MTase低下→メチル化消失→転写変動の時間的順序を厳密に区別できない
4. **4mC vs 6mA**: afsSは4mC、redZは6mAと修飾型が異なり、同一のMTaseが両方を担うかは不明
5. **afsRの安定メチル化**: afsRは6mA安定（AAGCCCGではないモチーフ）であり、本カスケードとは独立の制御系

---

## 5. Figure解説

### aagcccg_cascade_analysis.png

**ファイル**: `11_epigenome_integration/analysis/19_bgc_regulator_overview/aagcccg_cascade_analysis.png`

**説明**: AAGCCCGメチル化カスケードの3つの優先度A解析を統合した4パネル構成のFigure。

#### Panel A-1: afsS プロモーター構造

**読み取り方**:
- X軸: TSSからの距離（bp）。TSS = 0、左が上流（プロモーター方向）
- 黒い太線: DNA骨格
- 青い矢印 (TSS): 転写開始点
- オレンジ枠 (-10): σ因子-10 box（配列: TAGACT）
- 青枠 (-35): 推定-35 box
- 黄色破線枠: AAGCCCGモチーフ領域
- ロリポッププロット: T1（赤）/T2（橙）/T3（緑）のメチル化頻度
  - 塗りつぶし丸 = メチル化あり（頻度%表示）
  - 白抜き丸+破線 = メチル化消失（0%）
- 緑矢印: 遺伝子本体
- 右上ボックス: 発現変動（LFC、有意性）

**注目点**: 4mC部位がTSS上流172 bpのAAGCCCGモチーフ上に位置し、T1の83%からT2で完全消失。

**サポートするInsight**: #1, #6

#### Panel A-2: redZ プロモーター構造

**読み取り方**: Panel A-1と同一形式。

**注目点**: 6mA部位がTSS上流108 bpに位置。afsS（172 bp）よりTSSに近い。T1で59%の中程度のメチル化がT2で完全消失。

**サポートするInsight**: #1, #6

#### Panel B: SC_RS17645 (MTase) 発現 vs メチル化消失

**読み取り方**:
- X軸: 3タイムポイント
- 左Y軸（紫）: SC_RS17645の正規化発現量（エラーバー = SD, n=3）
- 右Y軸: メチル化頻度（%）
  - 赤丸破線: afsS 4mC頻度
  - 青三角破線: redZ 6mA頻度
  - × マーク: メチル化消失（0%）

**注目点**: T1→T2でMTase発現が4.6倍低下し、同時に両プロモーターのメチル化が消失。T3でMTaseが部分回復してもメチル化は回復しない。

**サポートするInsight**: #2, #3

#### Panel C: ゲノムワイドAAGCCCGメチル化

**読み取り方**:
- 棒グラフ: メチル化/非メチル化AAGCCCGサイト数
- オレンジ棒内: afsS/redZを含むメチル化サイトグループ

**注目点**: AAGCCCGの78.8%がメチル化 → 高い認識効率。

**サポートするInsight**: #4

#### Panel D: BGC制御因子のAAGCCCGステータス

**読み取り方**:
- Y軸: 12のBGC制御因子
- X軸: 3つの条件（AAGCCCGプロモーター存在、メチル化あり、メチル化消失）
- 色付き丸+チェック = 条件を満たす
- 灰色丸 = 条件を満たさない
- 黄色ハイライト行 = 全3条件を満たす因子

**注目点**: afsS/redZのみが全3条件を満たす。

**サポートするInsight**: #5

---

## 6. 出力ファイル一覧

| ファイル | 場所 | 内容 |
|---------|------|------|
| `aagcccg_cascade_analysis.png` | `19_bgc_regulator_overview/` | メインFigure (300 dpi, 4パネル) |
| `aagcccg_cascade_analysis.pdf` | `19_bgc_regulator_overview/` | PDF版 |
| `aagcccg_cascade_analysis.svg` | `19_bgc_regulator_overview/` | SVG版（編集可能） |
| `aagcccg_cascade_summary.tsv` | `19_bgc_regulator_overview/` | ゲノムワイド統計サマリー |
| `plot_aagcccg_cascade_analysis.py` | `scripts/` | Figure生成スクリプト |

---

## 7. 論文への組み込み提案

### 7.1 Results節（推奨表現）

> "Examination of the exact methylation site positions revealed that the 4mC site in the *afsS* promoter (position −172 relative to TSS) and the 6mA site in the *redZ* promoter (position −108) both reside within AAGCCCG motifs, the novel recognition sequence of the SC_RS17645 methyltransferase identified in this study. Strikingly, SC_RS17645 expression decreased 4.6-fold at T2 (log₂FC = −2.19, *P*_adj = 6.5 × 10⁻¹⁶), coinciding precisely with the complete loss of methylation at both promoters (Fig. X). Genome-wide, 78.8% of AAGCCCG motifs within 500 bp of annotated TSSs were methylated, and 163 of 416 AAGCCCG-containing promoters showed coordinated changes in methylation and gene expression (39.2%). Among 12 major BGC regulatory genes, only *afsS* and *redZ* harbored methylated AAGCCCG motifs in their promoter regions, and both lost methylation at T2."

### 7.2 Discussion節（推奨表現）

> "These findings support a three-step epigenetic cascade model in which (1) decreased expression of the SC_RS17645 methyltransferase during the growth phase transition leads to (2) passive demethylation at AAGCCCG motifs in key regulatory promoters, which in turn (3) modulates the transcription of signal-integration nodes in the BGC regulatory hierarchy. The irreversibility of demethylation despite partial recovery of methyltransferase expression at T3 is consistent with a passive, replication-dependent demethylation mechanism rather than active demethylase activity.
>
> Notably, the two affected regulators—AfsS, a master integrator of nutritional stress and secondary metabolism (Lee et al. 2002; Lian et al. 2008), and RedZ, the bldA-dependent gatekeeper of Red biosynthesis (White et al. 1997)—occupy analogous positions in the regulatory hierarchy as nodes linking global signals to pathway-specific activation. This convergence suggests that AAGCCCG methylation may serve as a selective epigenetic filter acting specifically at signal-integration points, rather than broadly regulating the entire BGC cascade."

### 7.3 論文Figure構成案

本Figureは以下のいずれかで論文に組み込むことを提案する：

- **Figure 3** (メイン): GRN-メチル化統合解析として、既存Figure 3を本Figureで置換
- **Figure 5** (追加メインFigure): カスケードモデルの独立Figure
- **Supplementary Figure**: Panel C/Dをメイン、Panel A/Bをサプリメンタリー

---

## 8. 次のステップへの示唆

### 8.1 直近で可能な解析（既存データ）

1. **163遺伝子の機能エンリッチメント**: AAGCCCGで発現協調変動する163遺伝子のGO/KEGG/COGエンリッチメントにより、AAGCCCGメチル化の標的機能カテゴリを同定
2. **ATCC比較**: 他のStreptomyces種でのAAGCCCGモチーフ保存性を、afsS/redZプロモーター領域に限定して再評価
3. **redZプロモーターの2つ目のAAGCCCG**: TSS -417bpに位置する非メチル化AAGCCCGの配列コンテキストを比較し、メチル化選択性の決定因子を推定

### 8.2 論文準備

1. 本レポートのFigure + 前レポートのFigureを統合し、論文Figureのfinal版を作成
2. Supplementary Tableとして、27制御因子の完全な発現・メチル化データを整理
3. カスケードモデルの概念図（schematic diagram）を作成

---

## 9. 参考文献

- Lee, P. C., Umeyama, T., & Horinouchi, S. (2002). AfsS is a target of AfsR, a transcriptional factor with ATPase activity that globally controls secondary metabolism in *Streptomyces coelicolor* A3(2). *Molecular Microbiology*, 43(6), 1413–1430.
- Lian, W., et al. (2008). Genome-wide transcriptome analysis reveals that a pleiotropic antibiotic regulator, AfsS, modulates nutritional stress response in *Streptomyces coelicolor* A3(2). *BMC Genomics*, 9, 56.
- White, J., Bibb, M. J., & Buttner, M. J. (1997). bldA dependence of undecylprodigiosin production in *Streptomyces coelicolor* A3(2) involves a pathway-specific regulatory gene, redZ. *Journal of Bacteriology*, 179(3), 627–635.
- Jeong, Y., et al. (2016). The dynamic transcriptional and translational landscape of the model antibiotic producer *Streptomyces coelicolor* A3(2). *Nature Communications*, 7, 11605.
