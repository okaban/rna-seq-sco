# H6: Genome-wide Regulatory Gene Methylation-Expression Screening

**Created:** 2026-02-24
**Last updated:** 2026-02-24
**Analysis directory:** `11_epigenome_integration/analysis/29_genomewide_TF_screen/`

---

## Key Insights

| # | Finding | Significance |
|---|---------|-------------|
| 1 | GFFから1,055個の制御関連遺伝子を抽出（25ファミリー） | 文献37 TFリストの約28倍の対象を網羅 |
| 2 | 165/1,055 (15.6%) がメチル化を保有 | 文献37リスト (6.7%) より2.3倍高いメチル化率 |
| 3 | 62個の制御因子でメチル化-発現協調変動を検出 | 文献37リストでは0件だったが、ゲノムワイドでは多数存在 |
| 4 | 協調変動の62個は全て37リスト外 | 文献バイアスが深刻 -- 重要TFを見落としていた |
| 5 | Concordant（repression + derepression）= T2: 11件, T3: 17件 | メチル化による転写抑制モデルを支持するケースが存在 |
| 6 | Discordant（gain+up, loss+down）= T2: 18件, T3: 25件 | 単純な抑制モデルだけでは説明不可 |
| 7 | MerR (28.6%), ROK (21.7%), Sensor kinase (20.5%) が高メチル化率 | SARP, Repressor, Anti-sigma, Lrp は0% |
| 8 | ramR (SCO6685) が唯一の文献TFでメチル化+DEG+協調あり | T3でgained methylation + upregulation (discordant) |
| 9 | 6mA主体 87遺伝子, 4mC主体 60遺伝子, 両方 18遺伝子 | 修飾タイプにファミリー間偏りあり（Sigma factor: 6mA優位, MerR: 4mC優位） |
| 10 | 文献37リストはDEG率73.3%だがメチル化率6.7%と低い | 二次代謝CSRは転写制御中心でエピジェネティック制御の標的ではない |

---

## 1. Background and Hypothesis

### 背景
Loop 2のH4解析で、文献ベースの37個のTF（bldA/B/C/D, afsR, actII-ORF4, redD/Z, dasRなど）を再スクリーニングしたところ、3/37のみがメチル化を保持し、全てconstitutiveで協調変動はゼロだった。しかし、S. coelicolor M145のGFFには約865個以上の制御関連遺伝子が存在する。

### 仮説
> GFF全制御因子（~1,055）を対象にすると、37 TFリスト外のTFにメチル化-発現協調が存在する

### 結論
**仮説は支持された。** 62個の制御因子でメチル化-発現協調変動が検出され、全て文献37リスト外であった。

---

## 2. Methods

### 2.1 GFFからの制御遺伝子抽出
NCBI RefSeq GFF (GCF_000203835.1) のCDS featureから、product fieldに以下のキーワードを含む遺伝子を抽出:
- TFファミリー名: TetR, MarR, GntR, LysR, AraC, LacI, IclR, MerR, WhiB, Lrp, XRE, SARP, ArsR, DeoR, PadR, ROK
- 機能語: regulator, transcription, DNA-binding, response regulator, sensor kinase, sigma factor, repressor, activator, helix-turn-helix, two-component, anti-sigma, CRP/FNR

### 2.2 データ統合
- **メチル化データ**: `integrated_methyl_expression_weighted.csv` (6mA/4mC x T1/T2/T3)
- **発現データ**: DESeq2 T2vsT1, T3vsT1 (|log2FC| >= 1, padj < 0.05)
- **文献TFリスト**: `literature_tf_master.csv` (37 TF)

### 2.3 協調変動の判定
- **Concordant repression**: gained methylation + decreased expression
- **Concordant derepression**: lost methylation + increased expression
- **Discordant**: gained+up or lost+down

### 2.4 クロスバリデーション
- GFF抽出1,055 locus_tagのうち1,033 (97.9%) がmethyl-expr dataに存在
- 文献37 TFのうち30がGFF制御遺伝子セットに含まれる（7個は非典型的制御因子: bldB, afsK, afsS, absB, nsdB, papR2, hrdA）

---

## 3. Results

### 3.1 制御遺伝子の全体像

| Metric | Value |
|--------|-------|
| Total regulatory genes | 1,055 |
| TF families | 25 |
| With methylation (any timepoint) | 165 (15.6%) |
| 6mA only | 87 |
| 4mC only | 60 |
| Both 6mA + 4mC | 18 |
| DEGs (any comparison) | 654 (62.0%) |
| Methylated + DEG | 114 (10.8%) |
| In literature 37 list | 30 (of 37; 7 non-canonical) |

### 3.2 文献37 vs ゲノムワイド比較

| Metric | Literature 37 | Other Regulatory (n=1,025) | All Regulatory (n=1,055) |
|--------|:---:|:---:|:---:|
| Methylation rate | **6.7%** | **15.9%** | 15.6% |
| DEG rate | **73.3%** | 61.7% | 62.0% |
| Methylated + DEG | **3.3%** | **11.0%** | 10.8% |

文献37リストは**転写レベルでの発現変動は高い (73.3%)**が、**メチル化保有率は著しく低い (6.7%)**。これはCSR（cluster-situated regulator）やグローバルTFが主にタンパク質-タンパク質相互作用やリン酸化で制御されていることと整合する。一方、ゲノムワイドの制御因子は2.3倍高いメチル化率を示した。

### 3.3 TFファミリー別メチル化率

| Family | n | Methylated | Rate (%) | DEGs | Both |
|--------|:-:|:---:|:---:|:---:|:---:|
| MerR | 28 | 8 | **28.6** | 15 | 3 |
| ROK | 23 | 5 | **21.7** | 14 | 3 |
| Sensor kinase | 83 | 17 | **20.5** | 55 | 12 |
| LacI | 35 | 7 | **20.0** | 20 | 3 |
| Sigma factor | 78 | 15 | **19.2** | 52 | 10 |
| LysR | 42 | 8 | 19.0 | 18 | 6 |
| HTH (other) | 163 | 28 | 17.2 | 97 | 22 |
| Other regulatory | 126 | 21 | 16.7 | 90 | 17 |
| TetR | 178 | 26 | 14.6 | 106 | 15 |
| Response regulator | 81 | 12 | 14.8 | 52 | 10 |
| MarR | 46 | 3 | 6.5 | 23 | 2 |
| GntR | 48 | 5 | 10.4 | 33 | 5 |
| **SARP** | **7** | **0** | **0.0** | **6** | **0** |
| **Repressor** | **7** | **0** | **0.0** | **6** | **0** |
| **Anti-sigma** | **3** | **0** | **0.0** | **2** | **0** |
| **Lrp** | **16** | **0** | **0.0** | **13** | **0** |

MerR (28.6%) が最高メチル化率。SARP, Repressor, Anti-sigma, Lrpはメチル化が完全に欠如。

### 3.4 メチル化変動パターン

| Pattern | Count |
|---------|:-----:|
| Constitutive (T1=T2=T3) | 58 |
| Gained T2 | 22 |
| Gained T3 | 22 |
| Gained T2+T3 | 15 |
| Lost T2 | 13 |
| Lost T3 | 18 |
| Lost T2+T3 | 16 |

変動メチル化 (gained/lost) は107/165 (64.8%) で、constitutiveは35.2%のみ。ゲノムワイドでは動的なメチル化パターンが多数。

### 3.5 メチル化-発現協調変動

#### T2vsT1

| Category | Count |
|----------|:-----:|
| Concordant repression (gained + down) | 5 |
| Concordant derepression (lost + up) | 6 |
| Discordant gain+up | 10 |
| Discordant loss+down | 8 |
| Methyl change, no expr change | 59 |

#### T3vsT1

| Category | Count |
|----------|:-----:|
| Concordant repression (gained + down) | 8 |
| Concordant derepression (lost + up) | 9 |
| Discordant gain+up | 15 |
| Discordant loss+down | 10 |
| Methyl change, no expr change | 42 |

#### 合計: 62個のユニークな協調変動制御因子

Concordant（メチル化による転写抑制/脱抑制モデルと一致）は計28件（T2: 11, T3: 17）、Discordant（逆方向）は43件（T2: 18, T3: 25）。

### 3.6 注目すべき協調変動遺伝子

#### Concordant repression (gained methylation + downregulation)

| Locus tag | SCO | Product | Family | Methyl change | Timepoint |
|-----------|-----|---------|--------|:---:|:---:|
| SC_RS07910 | SCO1191 | MarR family regulator | MarR | gained_T2 | T2 |
| SC_RS06045 | SCO0823 | GntR family regulator | GntR | gained_T3 | T3 |
| SC_RS07000 | SCO1008 | ArsR/SmtB family TF | ArsR | gained_T2+T3 | T3 |
| SC_RS20775 | SCO3729 | HTH domain protein | HTH | gained_T3 | T3 |
| SC_RS21680 | SCO3907 | ssDNA-binding protein | Other | gained_T2 | T2 |
| SC_RS23710 | SCO4305 | TetR/AcrR regulator | TetR | gained_T3 | T3 |
| SC_RS24635 | SCO4495 | UdgX uracil-DNA binding | Other | gained_T2+T3 | **T2+T3** |
| SC_RS25990 | SCO4766 | LysR family regulator | LysR | gained_T3 | T3 |
| SC_RS28630 | SCO5289 | Sensor histidine kinase | Sensor kinase | gained_T2+T3 | T3 |
| SC_RS29355 | SCO5434 | Response regulator | Response reg. | gained_T3 | T3 |
| SC_RS36820 | SCO6924 | HTH domain protein | HTH | gained_T2 | T2 |
| SC_RS40040 | SCO7568 | ROK family protein | ROK | gained_T2 | T2 |

**SC_RS24635 (SCO4495)** はT2, T3の両方でconcordant repressionを示す唯一の遺伝子。UdgX family uracil-DNA binding proteinで、DNA修復/修飾と関連。

#### Concordant derepression (lost methylation + upregulation)

| Locus tag | SCO | Product | Family | Methyl change | Timepoint |
|-----------|-----|---------|--------|:---:|:---:|
| SC_RS10435 | SCO1684 | Tetratricopeptide repeat | TetR | lost_T2+T3 | **T2+T3** |
| SC_RS31385 | SCO5828 | Response regulator | Response reg. | lost_T2+T3 | **T2+T3** |
| SC_RS35525 | SCO6668 | Sensor histidine kinase | Sensor kinase | lost_T2+T3 | **T2+T3** |
| SC_RS09835 | SCO1564 | Sigma-70 factor | Sigma factor | lost_T2 | T2 |
| SC_RS04270 | SCO0471 | HTH domain protein | HTH | lost_T2 | T2 |
| SC_RS22135 | SCO4005 | SigE family sigma | Sigma factor | lost_T2+T3 | T3 |
| SC_RS33075 | SCO6165 | TraR/DksA regulator | Other | lost_T2+T3 | T3 |
| SC_RS40435 | SCO7648 | Response regulator TF | Response reg. | lost_T3 | T3 |
| SC_RS37035 | -- | WhiB family regulator | WhiB | lost_T2 | T2 |

**SC_RS10435 (SCO1684)**, **SC_RS31385 (SCO5828)**, **SC_RS35525 (SCO6668)** の3遺伝子はT2, T3両方でconcordant derepressionを示す。特にSC_RS35525はsensor histidine kinaseで、二成分制御系のシグナル伝達上流に位置し、メチル化消失による脱抑制がシグナルカスケードに影響する可能性がある。

### 3.7 ramR (SCO6685) -- 唯一の文献TFで協調変動

ramR (SC_RS35610) は文献37リスト外だがnotable TFとして確認:
- **メチル化**: T1=1, T2=1, T3=2 (gained_T3)
- **発現**: T2vsT1 up, T3vsT1 up
- **協調**: T3でdiscordant (gained + up)
- ramRはaerialミセリウム形成と抗生物質産生の正の調節因子。メチル化増加にもかかわらず発現が上昇しており、他の転写活性化因子による制御が優位と考えられる。

### 3.8 Notable TFs not in literature 37 list

| Gene | Locus tag | SCO | Family | Methylated | DEG |
|------|-----------|-----|--------|:---:|:---:|
| whiA | SC_RS11780 | SCO1950 | Other regulatory | No | Yes |
| whiH | SC_RS31335 | SCO5819 | Other regulatory | No | Yes |
| whiJ | SC_RS24870 | SCO4543 | Other regulatory | No | Yes |
| glnR | SC_RS22985 | SCO4159 | Response regulator | No | Yes |
| ramR | SC_RS35610 | SCO6685 | Response regulator | **Yes** | Yes |
| mtrA | SC_RS17165 | SCO3013 | Response regulator | No | Yes |
| mtrB | SC_RS17160 | SCO3012 | Sensor kinase | No | Yes |

whiA, whiH, whiJ, glnR, mtrA/Bは全て発現変動遺伝子だがメチル化は保持していない。ramRのみがメチル化+DEGの両方を示す。

---

## 4. Figures

### Figure 1: TF Family Methylation Rate
`figures/tf_family_methylation_rate.pdf`

TFファミリー別のメチル化保有率（n >= 3のファミリーのみ）。MerR (28.6%) が最高、SARP/Repressor/Anti-sigma/Lrp/PadR/Activatorは0%。ファミリーによってメチル化感受性が大きく異なる。

### Figure 2: Literature 37 vs Genome-wide Comparison
`figures/literature_vs_genomewide_comparison.pdf`

3パネル比較: (A) メチル化率 -- 文献37は6.7%と低く、ゲノムワイド15.6%の半分以下。(B) DEG率 -- 文献37は73.3%と高い。(C) Methylated+DEG -- 文献37は3.3%、ゲノムワイドは10.8%で3.3倍。

### Figure 3: Methylated Regulatory Gene Heatmap
`figures/methylated_regulatory_heatmap.pdf`

165個の全メチル化制御遺伝子のヒートマップ。左パネル: メチル化サイト数（T1/T2/T3）、右パネル: 発現log2FC。*=concordant, x=discordant。赤丸=文献37リスト所属。afsR (SC_RS24295) がconstitutive methylation (1,1,1) で唯一の文献TF。

### Figure 4: 6mA vs 4mC by TF Family
`figures/tf_family_6mA_vs_4mC.pdf`

6mAと4mCの分布はファミリー間で偏りがある。Sigma factor: 6mA (15.4%) >> 4mC (2.6%)。MerR: 4mC (21.4%) >> 6mA (10.7%)。Sensor kinase: 4mC (15.7%) > 6mA (7.2%)。修飾タイプの嗜好性がTFファミリーごとに異なる。

---

## 5. Biological Interpretation

### 5.1 文献37 TFリストの限界
文献ベースの37 TFは二次代謝の「古典的」制御因子に偏っており、メチル化研究の対象として不適切だった。これらのTFは主にタンパク質リン酸化（AfsR-AfsK）、リガンド結合（DasR-GlcNAc）、c-di-GMP（BldD）などの翻訳後修飾で制御されており、エピジェネティック修飾の標的ではない。

### 5.2 メチル化による制御の標的
ゲノムワイドスクリーニングにより、**メチル化による制御が優先的に働く制御因子群**が明らかになった:
- **MerR family** (28.6%): 金属応答転写因子。環境ストレス応答のエピジェネティック制御
- **Sensor kinase** (20.5%): 二成分制御系のセンサー部分。シグナル伝達上流のゲーティング
- **Sigma factor** (19.2%): RNAポリメラーゼのプロモーター認識。転写開始レベルでの制御

### 5.3 協調変動の方向性
Concordant（メチル化が抑制的に働くモデル）とDiscordant（逆方向）が拮抗している。これは:
1. メチル化の位置（プロモーター vs gene body）による効果の違い
2. 他の制御因子による転写活性化がメチル化抑制を上回るケース
3. メチル化がactivator結合を阻害 vs repressor結合を阻害する二方向効果

を示唆する。特にSC_RS24635 (concordant T2+T3) やSC_RS10435/SC_RS31385/SC_RS35525 (concordant derepression T2+T3) は**再現性のある協調変動**であり、機能検証の優先候補。

### 5.4 SARP完全非メチル化の意義
7個のSARP family全遺伝子がメチル化ゼロで、6/7がDEG。SARPは二次代謝BGCのcluster-situated regulatorで、転写制御は上位カスケード（AfsR, BldD等）による直接活性化に依存しており、エピジェネティック制御の必要がないことを示す。

---

## 6. Output Files

### Tables (`tables/`)
| File | Rows | Description |
|------|:----:|-------------|
| `all_regulatory_genes.tsv` | 1,055 | 全制御因子の統合テーブル |
| `regulatory_methylation_summary.tsv` | 25 | ファミリー別メチル化率サマリー |
| `coordinated_regulatory_genes.tsv` | 62 | 協調変動を示す制御因子 |
| `missing_from_literature_list.tsv` | 682 | 37リスト外でメチル化またはDEGの制御因子 |

### Figures (`figures/`)
| File | Format | Description |
|------|--------|-------------|
| `tf_family_methylation_rate` | PDF, SVG | TFファミリー別メチル化率 |
| `literature_vs_genomewide_comparison` | PDF, SVG | 文献37 vs ゲノムワイド比較 |
| `methylated_regulatory_heatmap` | PDF, SVG | メチル化制御因子ヒートマップ |
| `tf_family_6mA_vs_4mC` | PDF, SVG | 6mA vs 4mC修飾タイプ比較 |

### Script
| File | Description |
|------|-------------|
| `run_genomewide_tf_screen.py` | 解析スクリプト（GFF解析、データ統合、Figure生成） |

---

## 7. Implications for Next Steps

1. **SC_RS24635 (SCO4495, UdgX)**: T2+T3で再現性のあるconcordant repressionを示す最優先候補。DNA修復との関連で機能的意義が高い。
2. **SC_RS10435/SC_RS31385/SC_RS35525**: T2+T3 concordant derepressionの3遺伝子。二成分制御系センサーキナーゼ (SC_RS35525) は特に重要。
3. **MerR family**: 最高メチル化率28.6%。環境ストレス応答のメチル化制御ネットワークとして体系的解析の価値がある。
4. **文献リスト拡張の必要性**: 37 TFでは不十分。ゲノムワイド1,055遺伝子から機能的に重要な制御因子をエピジェネティック観点で再選定する必要がある。

---

*Analysis performed with conda environment `rna-seq`.*
