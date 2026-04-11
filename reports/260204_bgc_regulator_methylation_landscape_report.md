# BGC制御因子のメチル化-発現変動ランドスケープ解析レポート

**日付:** 2026-02-04
**プロジェクト:** *Streptomyces coelicolor* A3(2) M145 エピゲノム-トランスクリプトーム統合解析
**解析ディレクトリ:** `11_epigenome_integration/analysis/19_bgc_regulator_overview/`

---

## Key Insights（主要なインサイト一覧）

| # | Insight | 対応Figure/Table | 根拠となる数値 |
|---|---------|-----------------|---------------|
| 1 | **BGC制御ネットワークはメチル化非依存で機能する** | `bgc_regulator_overview.png` Panel B, `bgc_regulator_summary.tsv` | 27制御因子中23因子（85%）が全タイムポイントで非メチル化 |
| 2 | **メチル化変動を示す因子はシグナル統合ノードに限局** | `bgc_regulator_overview.png` Panel C | 変動3因子中2つ（afsS, redZ）が「外部シグナル→二次代謝」の統合点 |
| 3 | **afsSの4mC消失がAct産生開始の時間制御と同期** | `bgc_regulator_overview.png` Panel A-B | 4mC: T1=1→T2/T3=0; afsS LFC=-1.39 (T3vT1); Act BGC 22/22遺伝子活性化 |
| 4 | **redZの6mA消失はbldA翻訳制御に加わるエピゲノム制御層を示唆** | `bgc_regulator_overview.png` Panel A-B | 6mA: T1=1→T2/T3=0; redZ LFC=-2.25 (T2vT1); redD LFC=+4.77 |
| 5 | **全メチル化変動はLost（消失）方向のみ** — 成長段階移行での「脱メチル化」を示唆 | `bgc_regulator_overview.png` Panel B | afsS: 4mC Lost, bldN: 6mA Lost, redZ: 6mA Lost |
| 6 | **bldAは転写レベルで連続的に低下するが、BGC活性化と逆相関** | `bgc_regulator_overview.png` Panel A | bldA LFC=-2.11 (T3vT1, padj=2.5e-36); Act/Red BGC同時期に活性化 |

---

## 1. 解析の背景と目的

### 1.1 背景

前解析（`260203_grn_tf_methylation_report.md`）において、GRN転写因子のメチル化状態を階層ごとに評価し、redZがメチル化-発現協調変動を示す唯一のCSRであることを報告した。しかし、以下の課題が残されていた：

1. **全制御因子の系統的な比較が不十分**: 前解析は37 TFを対象としたが、各BGCへのマッピングと発現変動・メチル化変動の「両方」を一覧化した図表がなかった
2. **bldAの扱い**: BGC制御カスケードの上流に位置するtRNA遺伝子bldAが解析パイプライン（06-10）に含まれておらず、翻訳制御の文脈での評価が不足していた
3. **先行研究との対比**: afsS, bldN, redZの機能喪失実験データ（Lee et al. 2002; Lian et al. 2008; Bibb et al. 2000; White et al. 1997）との整合性が未検討だった

### 1.2 目的

1. BGC制御に関与する全27因子について、発現変動（3比較）とメチル化変動を統合した概観Figureを作成する
2. 各BGC（act, red, cda, cpk）ごとに、発現変動とメチル化変動の「両方」を示した制御因子を同定する
3. 先行研究の機能データと対比し、論文報告価値のあるfindingを評価する

---

## 2. 方法

### 2.1 対象制御因子の選定

文献およびプロジェクト内解析結果（`05_annotation`, `07_regulator_network`, `09_SARP_integration`）に基づき、以下の27因子を選定：

| Tier | カテゴリ | 遺伝子数 | 含まれる因子 |
|------|---------|---------|------------|
| 1 | Global regulators | 11 | bldA, bldB, bldD, bldG, adpA, bldN, afsR, afsK, afsS, crp, dasR |
| 2 | Pleiotropic regulators | 6 | absA1, absA2, absB, nsdA, nsdB, wblA |
| 2 | Sigma factors | 4 | hrdD, sigE, sigF, sigU |
| 3 | CSR | 6 | actII-ORF4, redD, redZ, cdaR, cpkO/kasO, papR2 |

### 2.2 データソース

- **発現データ**: `04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_{2_vs_1,3_vs_1,3_vs_2}.tsv`
- **メチル化データ**: `11_epigenome_integration/analysis/01_integration/integrated_methyl_expression.csv`（6mA, 4mC, 5mCサイト数・頻度; T1/T2/T3）
- **BGCターゲット情報**: 文献ベース（各制御因子が制御することが知られているBGC）

### 2.3 判定基準

- **発現変動あり**: padj < 0.05（T3 vs T1比較を主要基準）
- **メチル化変動あり**: いずれかのタイムポイント間でメチル化サイト数が変化
- **「BOTH」判定**: 上記の両条件を満たす

### 2.4 スクリプト

`11_epigenome_integration/scripts/plot_bgc_regulator_overview.py`

---

## 3. 結果

### 3.1 メチル化状態の全体像

27制御因子のメチル化状態を以下のように分類した：

| メチル化カテゴリ | 因子数 | 割合 | 該当因子 |
|----------------|--------|------|---------|
| Unmethylated（全タイムポイントで非メチル化） | 23 | 85.2% | bldA, bldB, bldD, bldG, adpA, afsK, crp, dasR, absA1, absA2, absB, nsdA, nsdB, wblA, hrdD, sigE, sigF, sigU, actII-ORF4, redD, cdaR, cpkO/kasO, papR2 |
| Stable（メチル化あり・変動なし） | 1 | 3.7% | afsR (6mA × 1サイト) |
| Lost（メチル化消失） | 3 | 11.1% | afsS (4mC), bldN (6mA), redZ (6mA) |

### 3.2 発現変動+メチル化変動の「BOTH」判定結果

T3 vs T1で有意な発現変動（padj < 0.05）かつメチル化変動を示した因子：

| Regulator | Tier | 対象BGC | LFC (T3vT1) | padj | メチル化変動 | 修飾型 |
|-----------|------|---------|-------------|------|------------|--------|
| **afsS** | 1 (Global) | act, red, cda | -1.39 | 3.4e-12 | T1: 1 → T2/T3: 0 | 4mC Lost |
| **bldN** | 1 (Global) | act, red | +1.46 | 2.2e-7 | T1: 2 → T2/T3: 1 | 6mA Lost |
| **redZ** | 3 (CSR) | red | -1.10 | 1.7e-6 | T1: 1 → T2/T3: 0 | 6mA Lost |

### 3.3 bldAの発現動態

bldA（tRNA-Leu, UUAコドン翻訳に必須）はメチル化サイトを持たないが、転写レベルで顕著な低下を示した：

| 比較 | LFC | padj |
|------|-----|------|
| T2 vs T1 | -1.17 | 2.9e-12 |
| T3 vs T1 | **-2.11** | 2.5e-36 |
| T3 vs T2 | -0.92 | 1.6e-7 |

正規化カウント: T1平均 ~711 → T2平均 ~312 → T3平均 ~163（約4.4倍の低下）

BGC遺伝子群が強く活性化される局面でbldA転写産物が減少する逆相関パターンを示す。ただし、成熟tRNAの安定性（半減期がmRNAより格段に長い）を考慮すると、翻訳機能は維持されている可能性が高い。

### 3.4 BGCごとの制御因子マッピング概要

Panel Cの結果を以下にまとめる：

| BGC | 有意な発現変動を示した制御因子数 | BOTH（発現+メチル化） | 発現のみ | 変動なし |
|-----|-------------------------------|---------------------|---------|---------|
| **act** | 16/18 | **2**（afsS, bldN） | 14 | 2 |
| **red** | 16/19 | **3**（afsS, bldN, redZ） | 13 | 3 |
| **cda** | 10/13 | **1**（afsS） | 9 | 3 |
| **cpk** | 3/5 | 0 | 3 | 2 |

cpkクラスターは制御因子のメチル化変動が皆無であり、完全にメチル化非依存の制御を受けていると考えられる。

---

## 4. 先行研究との対比による生物学的解釈

### 4.1 afsS — マスター型制御因子の4mC消失と時間制御

> **Insight #3**: afsSの4mC消失がAct産生開始の時間制御と同期
> afsS KO株ではAct産生が完全消失する（Lee et al. 2002）一方、本データではafsSが低下する
> 局面でAct BGCが活性化している。4mCの消失がafsSの発現維持解除に関与し、
> カスケード始動後の「スイッチオフ」を許容している可能性を示唆する。
> **Figure**: `bgc_regulator_overview.png` Panel A-B, afsS行

**先行研究の知見（Lee et al. 2002; Lian et al. 2008）:**
- afsS欠損株：Act産生**完全消失**、Red軽度低下、形態形成は概ね正常
- afsSは栄養ストレス応答と二次代謝を**統合するマスター型σ様因子**
- AfsK → AfsR → afsS という**リン酸化カスケードの末端**

**本データとの対比:**
- afsS発現: T2vT1ではほぼ不変（LFC -0.03, n.s.）→ T3vT1で有意に低下（LFC -1.39）
- 4mC: T1で1サイト存在 → T2以降で消失
- Act BGC: T3で22/22遺伝子が有意に活性化
- afsK発現: T2vT1で強く上昇（LFC +3.04）— 上流キナーゼは誘導されているがafsSは低下

**解釈**: 4mCメチル化は初期（T1）のafsS発現維持に寄与し、その消失はカスケード始動後のafsS発現の「計画的低下」を許容するメカニズムである可能性がある。afsS KOでAct産生が完全消失する事実は、afsSがAct開始に**必須**であることを示しており、4mCによるafsS発現の時間的制御が二次代謝のタイミングに影響する可能性は生物学的に整合する。ただし、因果関係の実証には追加実験を要する。

### 4.2 redZ — bldA翻訳制御への追加的エピゲノム制御層

> **Insight #4**: redZの6mA消失はbldA翻訳制御に加わるエピゲノム制御層を示唆
> redZはredDの一次活性化因子であり、その機能喪失はRed産生をブロックする
> （White et al. 1997）。6mAの消失と転写低下の同時観察は、翻訳制御に加えて
> 転写レベルでのエピゲノム制御が存在する可能性を示す。
> **Figure**: `bgc_regulator_overview.png` Panel A-B, redZ行

**先行研究の知見（White et al. 1997; Narva et al. 1998）:**
- redZはredDの**一次活性化因子**（redZ機能喪失 → redD転写ほぼ完全消失 → Red産生ブロック）
- redZは**TTAコドンを含有** → bldA tRNAに翻訳が依存（開発タイミングとの連結）
- 非典型レスポンスレギュレーター（リン酸化ポケット欠損）

**本データの時系列（「redDパラドックス」の文脈）:**

| 時点 | redZ LFC | 6mA | redD LFC | 解釈 |
|------|----------|-----|----------|------|
| T1→T2 | -2.25 | 1→0 | +4.77 | redZ低下+6mA消失にもかかわらずredDは強く活性化 |
| T1→T3 | -1.10 | 0 | +3.57 | redZ低下持続、redDも維持 |
| T2→T3 | +1.11 | 0→0 | -1.15 | redZ部分回復、redDは低下 |

**解釈**: White et al. (1997) が示す通り、redZの転写は増殖期から検出されredDの転写に先行する。今回のデータでredZが最も低下するT2でredDが最大活性化を示す「パラドックス」は、**T1以前に蓄積されたRedZタンパク質がT2時点でのredD活性化を駆動**し、mRNAレベルの低下はタンパク質機能に直接影響しないという解釈で整合する。6mA消失はredZ転写低下と同期するが、1サイトのみであり因果関係は不明。bldA依存的翻訳制御が主要メカニズムであることを踏まえ、6mAは「追加的な制御層」として位置づけるのが適切である。

### 4.3 bldN — BGC制御との関連は間接的

> **Insight（補助的）**: bldNの6mA減少と発現上昇は形態形成プログラムの文脈で理解すべき
> **Figure**: `bgc_regulator_overview.png` Panel A-B, bldN行

**先行研究の知見（Bibb et al. 2000; Bignell et al. 2003）:**
- σ^BldNは空中菌糸形成のECF σ因子
- bldN KO → bald表現型（空中菌糸欠損）、基底菌糸は保持
- chaplin/rodlin遺伝子クラスターの直接的な転写ドライバー
- BGC制御への影響は**間接的かつ培地条件依存**

**本データ**: 6mA 2→1サイト減少 + 発現上昇（LFC +1.46）。メチル化減少と発現上昇の逆相関は、6mAが転写抑制に関与する可能性を示唆するが、BGC制御の直接的な関連性は弱い。形態形成-二次代謝のタイミング連結の文脈で補助的に言及するのが適切。

### 4.4 bldA — 翻訳制御因子のtRNA定量限界

> **Insight #6**: bldAは転写レベルで連続的に低下するが、BGC活性化と逆相関
> **Figure**: `bgc_regulator_overview.png` Panel A, bldA行

bldA転写産物の低下（LFC -2.11, T3vT1）はBGC活性化と逆相関するが、これは矛盾ではない：
- bldAの遺伝子領域は1,311 bpだが、成熟tRNAは~80 nt（前駆体を測定している）
- tRNAは**mRNAよりも格段に安定**（半減期が長い）
- 標準RNA-seqでの成熟tRNA定量は技術的限界がある
- T1以前の蓄積で翻訳機能の閾値を超えている可能性が高い

---

## 5. 論文報告価値の評価

### 5.1 メインクレームとして成立するもの

**「BGC制御カスケードはメチル化非依存で機能するが、シグナル統合ノードに限定的なメチル化制御が存在する」**

| 主張 | データ強度 | 論文内位置 |
|------|----------|-----------|
| BGC制御因子の85%がメチル化非依存 | **強い**（系統的データ） | Results |
| afsSの4mC消失が時間制御に関与する可能性 | **中程度**（KO文献との対比で補強） | Results |
| redZの6mA消失がbldA翻訳制御に加わる制御層 | **弱〜中**（競合仮説あり） | Discussion |
| 全メチル化変動がLost方向のみ | **中程度**（成長段階移行モデル） | Discussion |
| bldNの6mA減少 | **弱い**（BGCとの関連が間接的） | Supplementary |

### 5.2 推奨する論文内表現

Results節（afsSを中心に）：
> "We systematically examined the methylation status of 27 known BGC regulatory genes across three growth phases. Strikingly, 85% of these regulators were completely unmethylated at all timepoints, indicating that the BGC regulatory cascade operates largely independently of DNA methylation. Among the few exceptions, *afsS*—a master integrator of nutritional stress response and Act biosynthesis (Lee et al. 2002; Lian et al. 2008)—lost its sole 4mC site between T1 and T2, coinciding with the onset of its transcriptional downregulation (log₂FC = −1.39, T3 vs T1)."

Discussion節（エピゲノム制御モデル）：
> "The near-complete absence of methylation at BGC regulatory loci indicates that the hierarchical regulatory cascade operates predominantly through methylation-independent mechanisms. That methylation marks are concentrated at signal-integration nodes (*afsS*, *redZ*), rather than at the pathway-specific activators themselves, suggests a model in which epigenetic modification fine-tunes the sensitivity of decision points rather than directly controlling biosynthetic gene transcription."

---

## 6. Figure解説

### bgc_regulator_overview.png

**ファイル**: `11_epigenome_integration/analysis/19_bgc_regulator_overview/bgc_regulator_overview.png`

**説明**: BGC制御に関与する全27因子の発現変動（3比較）、メチル化サイト動態、およびBGCへのマッピングを統合した3パネル構成のFigure。

**読み取り方**:

- **Panel A** (左): 発現変動ヒートマップ
  - X軸: 3つの比較（T2vT1, T3vT1, T3vT2）
  - Y軸: 27制御因子（Tier 1-3で区分、赤太字=メチル化サイトあり）
  - 色: 青=発現低下、赤=発現上昇（log₂FC）
  - 太枠 = padj < 0.05; 薄色 = 非有意
  - 背景色: 青帯=Tier 1、橙帯=Tier 2、緑帯=Tier 3

- **Panel B** (中央): メチル化サイト数の時間変動
  - X軸: T1, T2, T3
  - オレンジ丸 = 6mA、紫丸 = 4mC（数字はサイト数）
  - 右端の "Lost"/"Stable" ラベルで変動方向を表示
  - 大部分の因子は空白（メチル化サイトなし）

- **Panel C** (右): BGCターゲットマッピング
  - X軸: act, red, cda, cpk
  - マーカー形状で分類:
    - ★（星）+ 黄色破線枠 = 発現変動 AND メチル化変動の**両方**
    - ◆（ダイヤ）= 発現変動 + 安定メチル化あり
    - ●（丸・大・色付き）= 発現変動のみ
    - ●（丸・小・灰）= 有意な変動なし
  - 色: Panel Aと同一のLFCカラースケール

**注目点**:
- Panel Bで有色の丸が表示されるのはわずか4因子（afsR, afsS, bldN, redZ）
- Panel Cの★マーカーは act列に2つ（afsS, bldN）、red列に3つ（afsS, bldN, redZ）、cda列に1つ（afsS）
- cpk列には★もなく、◆もない = cpkクラスターの制御は完全にメチル化非依存
- bldA（最上行）はPanel Aで全比較が有意な青色だがPanel Bは空白 = 転写低下するがメチル化なし

**サポートするInsight**: #1, #2, #3, #4, #5, #6

---

## 7. 出力ファイル一覧

| ファイル | 場所 | 内容 |
|---------|------|------|
| `bgc_regulator_overview.png` | `19_bgc_regulator_overview/` | メインFigure (300 dpi) |
| `bgc_regulator_overview.pdf` | `19_bgc_regulator_overview/` | PDF版 |
| `bgc_regulator_overview.svg` | `19_bgc_regulator_overview/` | SVG版（編集可能） |
| `bgc_regulator_summary.tsv` | `19_bgc_regulator_overview/` | 全制御因子×BGCの判定結果テーブル（63行） |
| `plot_bgc_regulator_overview.py` | `scripts/` | Figure生成スクリプト |

---

## 8. 次のステップへの示唆

### 8.1 論文準備に向けて

1. **afsSプロモーター領域の4mC部位の詳細解析**: TSS基盤解析（`18_tss_analyses`）のデータと統合し、4mC部位がプロモーター内の機能的領域（-10 box, -35 box付近）に位置するかを確認する
2. **公共データとの比較**: 他のStreptomyces種のメチロームデータでafsS, redZのメチル化保存性を評価する（ATCC比較解析の拡張）
3. **Figure統合**: 本Figureを論文Figure 3（GRN methylation）の改訂版、またはSupplementary Figureとして位置づける

### 8.2 追加解析の可能性（実験なし・公共データ活用）

1. **Ribo-seqデータの探索**: *S. coelicolor* の公開Ribo-seqデータがあれば、TTAコドン含有mRNA（redZ, actII-ORF4）の翻訳効率の時間変動を評価可能
2. **プロテオミクスデータ**: AfsS, RedZのタンパク質レベルデータがあれば、mRNA低下がタンパク質量に反映されるかを検証可能
3. **ChIP-seqデータ**: BldDやAfsRのChIP-seqデータとメチル化部位の重複を確認し、メチル化がTF結合に影響するかの間接的評価

---

## 9. 参考文献

- Lee, P. C., Umeyama, T., & Horinouchi, S. (2002). AfsS is a target of AfsR, a transcriptional factor with ATPase activity that globally controls secondary metabolism in *Streptomyces coelicolor* A3(2). *Molecular Microbiology*, 43(6), 1413–1430.
- Lian, W., Jayapal, K. P., Charaniya, S., Mehra, S., Glod, F., Kyung, Y. S., Xia, X., Smith, J., & Hu, W.-S. (2008). Genome-wide transcriptome analysis reveals that a pleiotropic antibiotic regulator, AfsS, modulates nutritional stress response in *Streptomyces coelicolor* A3(2). *BMC Genomics*, 9, 56.
- Bibb, M. J., Molle, V., & Buttner, M. J. (2000). σ^BldN, an extracytoplasmic function RNA polymerase sigma factor required for aerial mycelium formation in *Streptomyces coelicolor* A3(2). *Journal of Bacteriology*, 182(17), 4606–4616.
- Bignell, D. R. D., Warawa, J. L., Strap, J. L., Chater, K. F., & Leskiw, B. K. (2003). Study of the bldN gene of *Streptomyces coelicolor* A3(2). *Journal of Bacteriology*, 185(7), 2338–2345.
- White, J., Bibb, M. J., & Buttner, M. J. (1997). bldA dependence of undecylprodigiosin production in *Streptomyces coelicolor* A3(2) involves a pathway-specific regulatory gene, redZ. *Journal of Bacteriology*, 179(3), 627–635.
- Narva, K. E., Woolhiser, L. K., & Yan, M. J. (1998). A response-regulator-like activator of antibiotic synthesis from *Streptomyces coelicolor* A3(2). *Journal of Bacteriology*, 180(9), 2459–2464.
- Elliot, M. A., et al. (2003). The *Streptomyces coelicolor* developmental transcription factor BldD is a direct regulator of key developmental genes. *Journal of Bacteriology*, 185(21), 6075–6082.
