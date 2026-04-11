# reports - 解析レポートディレクトリ

このディレクトリは、M145 RNA-seqプロジェクトの**解析レポート**を格納するための専用ディレクトリです。

解析ステップ（`01_qc/`, `02_alignment/`, ...）とは別カテゴリのため、番号なしのディレクトリ名を使用しています。

## 目的

オミクス解析では、解析の進行に伴って様々なインサイトが得られます。
これらを定期的にレポートとして記録することで：

- 解析の進捗と得られた知見を追跡できる
- 後から振り返りやすくなる
- 論文執筆時の参照資料として活用できる
- チームメンバーとの共有が容易になる

## ファイル命名規則

```
YYMMDD_<解析内容>_report.md
```

**命名のポイント:**
- `<解析内容>` は実施した解析が具体的に分かる名前にする
- スネークケース（アンダースコア区切り）を使用
- 略語は一般的なもののみ使用（DEGs, GO, KEGG, BGC, TF など）

例：
- `260202_DEGs_functional_enrichment_report.md` - DEGsの機能エンリッチメント解析
- `260128_DESeq2_differential_expression_report.md` - DESeq2による発現変動解析
- `260129_methylome_expression_integration_report.md` - メチローム・発現統合解析

## 推奨する記載内容

1. **日付とプロジェクト情報**
2. **解析の背景と目的**
3. **主要な発見（表形式で整理）**
4. **生物学的示唆**
5. **出力ファイルへの参照**
6. **次のステップへの示唆**

## レポート一覧

### 初期解析（2026-01-29）

| ファイル | 日付 | 内容 |
|---------|------|------|
| `260129_epigenome_integration_report.md` | 2026-01-29 | エピゲノム-トランスクリプトーム初回統合解析。3タイムポイント、844遺伝子にプロモーターメチル化 |
| `260129_final_integration_report.md` | 2026-01-29 | 統合解析最終版。3手法比較、5mC検証、coverage-weighted推奨 |
| `260129_T2vsT1_coordinated_analysis_report.md` | 2026-01-29 | T2vsT1プロモーターメチル化-発現相関。558遺伝子ペアの協調変動 |
| `260129_TF_methylation_analysis_report.md` | 2026-01-29 | TFメチル化-発現相関。24 TFの協調変動（NsdB, RamR等） |
| `260129_T3_2rep_comparison_report.md` | 2026-01-29 | T3 2レプリケート再解析。サンプル3-2低カバレッジの影響評価 |

### 包括的解析（2026-02-02〜02-03）

| ファイル | 日付 | 内容 |
|---------|------|------|
| `260202_BGC_coverage_track_analysis_report.md` | 2026-02-02 | BGCカバレッジトラック解析 |
| `260202_DEGs_functional_enrichment_report.md` | 2026-02-02 | DEGsの機能エンリッチメント解析（ベン図、GO/KEGG、COG、カバレッジトラック） |
| `260202_epigenome_transcriptome_integration_report.md` | 2026-02-02 | エピゲノム-トランスクリプトーム統合解析（**包括的レポート**: 新規性、DEG/DMG重複検定、マルチオミクストラック、経時変化、モチーフ位置、BGC/TF解析、**メチル化モチーフ解析（MEME）、協調的遺伝子COGエンリッチメント**） |
| `260203_coordinated_enrichment_report.md` | 2026-02-03 | 協調的遺伝子の機能エンリッチメント。558遺伝子ペア（27.4%正の相関） |
| `260203_motif_analysis_report.md` | 2026-02-03 | メチル化モチーフ解析サマリー。6mA 3,248 + 4mC 2,693 サイト |
| `260203_spurious_correlation_validation_report.md` | 2026-02-03 | **擬似相関検証解析**（並べ替え検定、偏相関、Bootstrap CI、負の対照群）- 4mC T2vsT1相関の頑健性を確認 |
| `260203_spurious_validation_analysis_report.md` | 2026-02-03 (updated 02-24) | Spurious Correlation Validation結果テーブル。4mC T2vsT1: r=0.1248, perm p=0.002 |
| `260203_atcc_comparative_methylome_report.md` | 2026-02-03 | **ATCC 37株比較メチローム解析** - 種間モチーフ保存率（CCGG 86.5%保存、AAGCCCG 21.6%検出=系統特異的）、M145の属内位置付け |
| `260203_atcc_comparative_methylome_analysis_report.md` | 2026-02-03 | ATCC比較メチローム解析（analysis dir版） |
| `260203_rm_system_identification_report.md` | 2026-02-03 | **R-M系同定解析** - AAGCCCG=REBASE未登録（新規）、SC_RS17645（N-6 DNA methylase）が候補、T2で発現低下→脱メチル化モデル |
| `260203_rm_system_identification_analysis_report.md` | 2026-02-03 | R-M系同定解析（analysis dir版） |
| `260203_grn_tf_methylation_report.md` | 2026-02-03 | **GRN TF-メチル化統合解析** - **redZが唯一の協調変動TF**（6mA Lost + expression DOWN）、SC_RS17645→AAGCCCG→redZ→Red BGCカスケード |
| `260203_comprehensive_analysis_summary_report.md` | 2026-02-03 | **包括的解析サマリー** - 追加解析4件、**論文Figure 1-4生成**、全Key Insights一覧、エピジェネティックカスケードモデル |
| `260203_presentation_manuscript_for_collaborators.md` | 2026-02-03 | 共同研究者向けプレゼンテーション原稿 |

### 追加解析（2026-02-04〜02-07）

| ファイル | 日付 | 内容 |
|---------|------|------|
| `260204_tss_based_epigenome_integration_report.md` | 2026-02-04 | **TSS基盤エピゲノム統合解析** - 実験的TSS導入(Jeong2016)、メタジーンプロファイル、距離帯別相関、σ因子-10 box枯渇(p=1.73e-12) |
| `260204_bgc_regulator_methylation_landscape_report.md` | 2026-02-04 | **BGC制御因子メチル化ランドスケープ解析** - 27制御因子の85%がメチル化非依存、先行研究KOデータとの対比 |
| `260204_aagcccg_cascade_bgc_regulation_report.md` | 2026-02-04 | AAGCCCGカスケードBGC制御解析 |
| `260204_rebase_comparative_methylome_report.md` | 2026-02-04 | **REBASE実データ比較メチローム解析** - AAGCCCG: 0.0%（属レベルで完全に新規） |
| `260204_summary_manuscript.md` | 2026-02-04 | サマリー原稿 |
| `260205_SC_RS17645_sequence_analysis_report.md` | 2026-02-05 | **SC_RS17645（SCO3104）配列・構造ホモロジー解析** - Type I R-M system HsdM型メチル化サブユニットに分類 |
| `260205_genuswide_motif_conservation_report.md` | 2026-02-05 | ***Streptomyces*属メチル化モチーフ保存性解析** - AAGCCCG R-M系はREBASE 82種中0種に未報告。属レベルO/E分布 |
| `260205_4mC_5mC_competition_report.md` | 2026-02-05 | **m4C/m5C二重シトシン修飾系検証** - 100%の4mCコンセンサスモチーフがGGCCGGコンテクスト |
| `260206_DMG_functional_enrichment_report.md` | 2026-02-06 | **DMG機能エンリッチメント解析（包括版）** - DMGはKEGG/GO/COGいずれでも有意エンリッチメントなし。「位置依存・機能非依存」モデル |
| `260206_expanded_motif_search_report.md` | 2026-02-06 (updated 02-07) | **拡張モチーフ探索解析** - **AAGCCCG二重修飾 Score 10/10 HIGH-CONFIDENCE NOVEL**、CCGKCA Score 7/10 BREX-2候補 |
| `260206_master_TF_list_binding_sites_report.md` | 2026-02-06 | **TFマスターリスト＋結合部位データベース構築** - RegPrecise、ZorroAranda、FIMO統合。Tier 1 BS 790件 |
| `260207_TF_binding_site_methylation_report.md` | 2026-02-07 (updated 02-24) | **TF Binding Siteメチル化変動解析** - BSでメチル化が有意に枯渇（fold=0.66, p=6.8×10⁻⁵）。75% Lost_T2T3 |

### Issue回答・再解析（2026-02-17〜02-24）

| ファイル | 日付 | 内容 |
|---------|------|------|
| `260217_issue-answer.md` | 2026-02-17 | **エピゲノム-トランスクリプトーム統合解析 Issue回答レポート** - エピゲノム「ゲートキーパー」モデル提案 |
| `260224_T3_resequencing_reanalysis_report.md` | 2026-02-24 | **T3追加シーケンシングデータによる全解析再実行レポート** - サンプル3-2カバレッジ14.3x→57.6x改善、MIN_REPS=3 unweighted統一。4mC T2vsT1相関維持、**6mA T2vsT1相関消失**。26スクリプト再実行 |
| `260224_redZ_paradox_report.md` | 2026-02-24 | **H3: redZ paradox定量的解析** - redZは実際には上昇(+0.91 log2FC)、redD 55倍増がRed BGC活性化の主因。活性化-抑制バランス正(+4.89 T2vsT1) |
| `260224_H4_TF_methylation_rescreen_report.md` | 2026-02-24 | **H4: 37 TF methylation-expression re-screen** - 修正locus_tagによる全TF再スクリーニング。3/37のみメチル化保持、0件の協調変動。SARP 4因子は全てメチル化無し(LFC +2~+7)。AAGCCCG motifはbldBプロモーターのみ |
| `260224_H5_MTase_stability_report.md` | 2026-02-24 | **H5: MTase expression stability vs methylation-expression correlation** - 5候補MTaseの発現安定性ランキング(BREX-2 PglX最安定, Dcm-like最不安定)。AAGCCCG 260サイト消失がSC_RS17645 LFC=-2.19と対応。遺伝子レベル相関は全モチーフ群で弱い(|rho|<0.16)。**部分的に支持**: サイトレベルでは支持、遺伝子レベルでは不支持 |
| `260224_H6_genomewide_TF_screen_report.md` | 2026-02-24 | **H6: ゲノムワイド制御因子メチル化-発現スクリーニング** - GFFから1,055制御遺伝子抽出(25ファミリー)。165/1,055 (15.6%)がメチル化保有（文献37の6.7%の2.3倍）。**62個で協調変動検出**（全て37リスト外）。SC_RS24635がT2+T3 concordant repression。MerR最高メチル化率28.6%、SARP完全非メチル化 |
| `260224_H7_CCGG_MTase_paradox_report.md` | 2026-02-24 | **H7: CCGG MTaseパラドックス解明** - 全CCGG 4mCサイトがタイムポイント間で**完全非重複（Jaccard=0.000）**。受動的希釈仮説を棄却。T1=core genome (77%), T2=chromosomal arms (89%)の劇的な地理的シフト。SC_RS36410はArgonaute+PD-(D/E)XK nucleaseと共誘導される防御島。1,516 T1サイトの責任酵素は未同定 |
| `260224_H8_coordinated_regulators_characterization_report.md` | 2026-02-24 | **H8: 62協調変動制御遺伝子の地理的・機能的特性解析** - 全体のarm enrichmentは非有意(p=0.17)。T3 discordant_gain_upのみ**arm強enrichment(87%, OR=8.05, p=0.001)**。COG T(シグナル伝達)有意enrichment(OR=2.47, p=0.045)、COG Q(二次代謝)は0。TCS 7ペア同定（非対称メチル化）。Top候補: SC_RS10435(chaplin隣接), SC_RS35525(PPTase隣接), SC_RS31385(citrate synthase隣接) |
| `260224_H9_CCGG_5mC_misclassification_report.md` | 2026-02-24 | **H9: CCGG "4mC"は5mC誤分類か？系統的検証** - C2位置84.5%でDcm標的と一致するが、**5mCシグナルは実質ゼロ(0.01%)**、4mC頻度は全修飾中最高(82.7%)、CCWGG文脈ではなくGCCGGCパリンドローム(91.4%)。**H9は棄却**: 信号品質指標は5mC誤分類と矛盾。GCCGGCを認識する未同定N4-C MTaseの可能性 |
| `260224_H10_GCCGGC_RM_identification_report.md` | 2026-02-24 | **H10: GCCGGC認識R-M系同定解析** - REBASE 87 NaeI-family中**2件のみm4C産生**（M.Svi27968I *S. violascens* + M.PfrJS2V *P. freudenreichii*）。*Streptomyces* 28 GCCGGC R-M中m4Cは1件のみ。M145の129 MTase中N4-C明示アノテーション無し。Top候補: SC_RS24685（DNA MTase, T1=325, 急減）、SC_RS03950（DNA MTase, T1=130, 恒常的）。**次ステップ: M.Svi27968I BLAST必須** |
| `260224_H11_gatekeeper_model_v2_report.md` | 2026-02-24 | **H11: Gatekeeper Model v2定量的統合** - H1-H9全結果を3層モデルに統合。Layer 1: Landscape Remodeling (4mC rho=0.139, Jaccard=0.000)。Layer 2: Protection/Depletion (TF BS fold=0.66, sigma-10 p=1.73e-12)。Layer 3: Signal Gating Revised (62非文献制御因子, TCS非対称, COG T OR=2.47)。**Excluded**: TFカスケード経路(0/37協調変動)。10個のTestable Predictions |
| `260225_H12_GCCGGC_MTase_BLAST_report.md` | 2026-02-25 | **H12: M.Svi27968I BLAST相同性検索によるGCCGGC N4-C MTase同定** - H10候補覆し: SC_RS24685はBLAST圏外、**SC_RS19770が最有力**(E=0.007, 29% identity)。SC_RS36410がdefense island内2位(E=0.047)。T1 core 1,516サイトの責任酵素は未同定のまま |
| `260225_H13_AAGCCCG_distribution_report.md` | 2026-02-25 | **H13: AAGCCCG 6mAサイトのゲノム分布パターン解析** - T1: 260サイト(core 68.8%, NN中央値 9,406bp)、T2: 64サイト(core 64.1%, NN 55,672bp)。**Hypothetical depleted**(p=0.039)、**Regulatory/TF depleted**(p=0.035)。BGC内ほぼゼロ。発現影響はbidirectional |
| `260226_H14_defense_island_GCCGGC_report.md` | 2026-02-26 | **H14: Defense Island--GCCGGC T3 Coupling解析** - SC_RS36410 defense island(6遺伝子)とGCCGGC 4mCの時空間的関連。T1=83%core(1,289sites)→T2=82%arm(407sites)→T3=62%arm(21sites)の劇的地理的シフト(chi2=597, p~10^-130)。T3サイト1個がdefense island内(SC_RS36410コーディング領域)。6遺伝子coordinated T3 induction(mean rho=0.694)。SC_RS19765 pseudogene=102bp(10.5%)。**Partial**: arm-enrichmentは存在するがT3サイト数21と少なく、T1 core酵素は未同定 |
| `260226_H15_cross_motif_regulatory_avoidance_report.md` | 2026-02-26 | **H15: Cross-Motif Regulatory Avoidance解析** - H13のRegulatory/TF depletion(AAGCCCG 6mA fold=0.43)が全モチーフで普遍的か検証。**SUPPORTED**: GCCGGC 4mC(fold=0.61, p=3.9e-05)、All 4mC(fold=0.62, p=1.3e-07)でも有意にdepleted。Hypothetical depletion(fold=0.45-0.56)も全モチーフで普遍的。Gatekeeper Model v2 Layer 2(Protection)を検証 |
| `260226_H16_GCCGGC_MTase_reverse_ID_report.md` | 2026-02-26 | **H16: GCCGGC MTase発現相関逆同定** - 129 MTaseのSpearman相関ランキング。37個がrho=1.0(T1>50: 36個)。Top候補SC_RS13615(SCO2317): composite score 10、CDD COG0863 DNA methylase hit(E=2.23e-03)、TnpB neighbor。SC_RS10665(SCO1731)は既知m5C MTase(m4Cではない)。BLAST候補SC_RS19770はT1=11で発現不足。**PARTIAL**: 発現相関候補は同定されたが、DNA cytosine MTaseドメインと配列相同性を同時に満たす候補は無し |
| `260226_H17_motif_division_of_labor_report.md` | 2026-02-26 | **H17: モチーフ特異的"分業"解析** - GCCGGC-proximal遺伝子とAAGCCCG-proximal遺伝子の発現動態比較。予測と逆: AAGCCCG近傍は最高のT3上昇(median LFC +0.568)、GCCGGC近傍は抑制的(median LFC -0.028, vs Background p=4.1e-10)。H15の機能的enrichmentは位置情報であり発現制御ではない。**REJECTED**: 機能カテゴリの"分業"は転写レベルに反映されず。GCCGGC 4mCの抑制的修飾効果(r=0.09)がGatekeeper Modelを補強 |
| `260226_H18_GCCGGC_dose_response_report.md` | 2026-02-26 | **H18: GCCGGCメチル化密度-発現量用量反応解析** - 2kb以内のGCCGGC 4mCサイト数による段階的発現抑制を検証。JT trend test p=4.6e-14(T2)、全遺伝子Spearman rho=-0.087(p=4.6e-14)で有意。しかし>=1サイト遺伝子内ではrho=-0.037(T2, p=0.026)、rho=-0.020(T3, NS)と極めて弱い。4+群でリバウンド、arm/core層別で信号消失。AAGCCCG 6mAは全遺伝子で非有意だがメチル化遺伝子内でrho=-0.146(T3, p=1e-05)と逆パターン。**PARTIAL**: 閾値スイッチ(0 vs >=1)であり段階的ブレーキではない |
| `260226_H19_GCCGGC_temporal_derepression_report.md` | 2026-02-26 | **H19: GCCGGC時間的脱抑制解析** - T1→T2メチル化消失遺伝子(Lost, n=3,226)の発現上昇を検証。予測と完全に逆: Lost遺伝子はmedian LFC=-0.253(Never=+0.006, p=8.3e-08)で**下方**シフト。Gained遺伝子(n=870)は逆にmedian LFC=+0.447で最高。Geographic stratificationで効果消失(core: p=0.87, arm: p=0.019逆方向)。AAGCCCG reciprocal解析も効果なし(p=0.92)。H17の抑制相関は地理的交絡。**REJECTED**: 時間的脱抑制は観察されず、メチル化-発現関連はSimpson's paradox（core/armの発現動態差）で説明される |
| `260226_H20_regulatory_avoidance_geographic_test_report.md` | 2026-02-26 | **H20: H15制御遺伝子回避の地理的層別化検証** - H19でSimpson's paradoxと判明した発現抑制とは対照的に、H15のRegulatory/TF depletion(fold=0.43-0.71)が地理的交絡ではなく真の生物学的シグナルかを検証。Core-only: 3/3モチーフで有意(GCCGGC p=7.5e-04, All 4mC p=3.8e-05)。Arm-only: 3/3モチーフで有意(p<0.015)。CMH region-adjusted OR: 全4モチーフで有意(p<0.007)。Hypotheticalも同様に真の枯渇。**SUPPORTED**: H15の制御遺伝子回避はcore/arm両方で再現し、地理的アーティファクトではない |
| `260226_H21_AAGCCCG_temporal_causality_report.md` | 2026-02-26 | **H21: AAGCCCG 6mA時間的因果性検証** - GCCGGC(H19)のSimpson's paradoxを踏まえ、地理的シフトが最小(4.7pp)のAAGCCCG 6mAで脱抑制を検証。Lost(n=836) vs Never(n=6,397): median LFC差=+0.004, p=0.91, r=0.003。全transition groupが~66% coreで地理的交絡なし。Core/arm層別でも非有意(p=0.28/0.23)。用量反応は逆方向(多サイト消失=抑制強化, rho=-0.107, p=0.002)。メチル化頻度もLFCと無相関。**REJECTED**: GCCGGC(H19)と合わせ、M145の2大メチル化系いずれも時間的脱抑制効果なし |
| `260226_H22_sequence_motif_depletion_report.md` | 2026-02-26 | **H22: 配列レベルモチーフ枯渇解析** - H15/H20の制御遺伝子メチル化回避が配列レベル（進化的反選択）かタンパク質占有（TF結合によるMTaseブロック）かを検証。TGGCCGGC 2,093サイト、AAGCCCG 1,334サイトをゲノムスキャン。Gene body: fold=0.74-0.76(p<3e-07)で有意枯渇。Extended 2kb: fold=0.85-0.92(p<0.007)。Promoter: fold~1.0(NS)。DNA fold(0.85-0.92) vs メチル化fold(0.43-0.61): 配列寄与率21-35%。Permutation p=0.001(AAGCCCG), 0.035(TGGCCGGC)。**PARTIAL**: 配列レベル枯渇は実在するが全体の21-35%のみ。残り65-79%はタンパク質占有モデル。プロモーターでDNA枯渇なし→TF結合による動的保護が主因 |
| `260226_H23_category_specificity_avoidance_report.md` | 2026-02-26 | **H23: 遺伝子カテゴリ特異性メチル化回避解析** - H15(サイト中心)の制御遺伝子回避が遺伝子中心の観点でも再現するか、他の重要カテゴリにも拡大するか検証。9カテゴリ(Regulatory/Translation/DNA repair/Cell division/Energy/Transport/Secondary metabolism/Hypothetical/Other)×3モチーフ。**制御遺伝子: gene-centric fold=1.01(NS)** — H15のsite-centric fold=0.61との不一致は密度vs二値近接の違いで説明。発現量Q5はQ1より1.2x高い近接率(JT z=+10-12, p~0)。BGC遺伝子は1.4-1.7x enriched。Hypotheticalのみ有意depleted(fold=0.86-0.88)。**REJECTED(universal avoidance)**: 重要遺伝子カテゴリの普遍的回避は不支持。H15回避は密度ベースの現象であり完全排除ではない |
| `260226_H24_BGC_methylation_geographic_test_report.md` | 2026-02-26 | **H24: BGCメチル化エンリッチメント地理的交絡検定** - H23のBGC enrichment(fold=1.59, p=4.4e-08)がcore共局在による交絡かを検証。Core-only: GCCGGC fold=1.26(p=0.002)で**有意に残存**、All 4mCはfold=1.07(NS)で消失。DNA配列モチーフ対照: TGGCCGGC body density fold=1.40(p=0.002)がメチル化fold(1.45)と一致→**配列組成効果**。ACT(1.49x)、RED(1.35x)が高密度、CPK(0.78x)は低密度。AAGCCCG 6mAはBGC内ほぼゼロ(H13と整合)。**PARTIAL**: 地理的交絡は~21%、残余のGCCGGC enrichmentは配列組成(GC-rich BGC coding regions)で説明。能動的メチル化ターゲティングの証拠なし |
| `260226_H25_TSS_methylation_gradient_report.md` | 2026-02-26 | **H25: TSS周辺メチル化空間勾配解析** - 制御遺伝子TSS周辺の高解像度メチル化密度プロファイル。制御遺伝子TSS枯渇17.3%(非制御3.8%の4.6倍)。Protection zone幅2,200bp(-1,300~+700bp)、最深枯渇+300bp(ratio=0.541, p_adj=0.005)。Core特異的(18.4% vs arm 14.5%)。全モチーフで再現(4mC ratio=0.805, 6mA ratio=0.873)。TCS/sensor kinaseは逆に高密度。**SUPPORTED**: 制御遺伝子特異的TSS保護、タンパク質占有シールドモデル支持 |
| `260226_H26_TFBS_methylation_protection_report.md` | 2026-02-26 | **H26: TF結合部位レベルのメチル化保護解析** - FIMO予測56,338 TF結合部位でのメチル化オーバーラップ/空間プロファイルを検証。GCCGGC 4mCはTF BSで逆に*enriched*(fold=1.157, p=0.006)。AAGCCCG 6mAは*depleted*(fold=0.629, p=0.002)。全メチル化の中心枯渇なし(permutation p=0.12-0.96)。制御遺伝子プロモーターTF BSは非制御より7-9%低密度だが非有意。Intergenic TF BSで6mAが3x enriched(p=4.6e-07)。**REJECTED**: 個々のTF結合部位はメチル化枯渇を示さない。H25の保護ゾーンは個別TF占有ではなく、制御遺伝子プロモーターの集合的タンパク質占有による集団効果 |
| `260226_H27_coordinated_regulators_protection_report.md` | 2026-02-26 | **H27: 62協調変動制御遺伝子のProtection Zone特性解析** - H25のprotection zoneが協調変動62遺伝子と非協調993遺伝子で異なるか検証。最近接メチル化サイト距離: 協調=114bp vs 非協調=762bp(6.7倍近接, p=5.3e-28)。協調遺伝子はprotection zone**完全欠如**(幅0bp, 深さ0.000)、代わりにTSS直上で**8.4倍エンリッチメント**(密度3.15 vs 0.37 sites/kb/gene)。Delta-LFC相関は協調群のみ有意(rho=0.277, p=0.029)。地理分布は非有意(OR=1.26, p=0.42)。4協調型全てで保護欠如。**SUPPORTED(逆方向)**: "Exposed Promoter"モデル - 993遺伝子はプロモーター保護で発現安定、62遺伝子はメチル化透過で発現応答 |
| `260226_H28_exposed_regulators_characteristics_report.md` | 2026-02-26 | **H28: 62 Exposed制御遺伝子の生物学的特性解析** - 62 exposed vs 993 shielded制御遺伝子の包括的多変量比較(19検定)。T1発現有意に低い(median 90 vs 123, p=0.010, r=0.195)。**全62遺伝子が動的発現(constitutive=0%, OR=inf, p=5.3e-13)**。|log2FC|が1.7倍大(T3: p=1.7e-09, r=0.456)。プロモーターGC含量高い(71.5% vs 69.7%, p=2.0e-04)。TFファミリー・地理分布・オペロン構造は非有意。T1最低四分位にexposed 45.2%集中(3.1倍enrichment)。**SUPPORTED**: 低初期発現+100%動的制御+高GCプロモーターの3因子シナジーモデル |
| `260226_H29_shielded_exposed_boundary_report.md` | 2026-02-26 | **H29: Shielded/Exposed制御遺伝子の定量的境界解析** - 1,017制御遺伝子の9特徴量で分類予測。nearest_methyl_distance単独AUC=0.917(95%CI: 0.895-0.935)が最良。**293bp閾値**で感度1.000/特異度0.806。baseMeanはAUC=0.547(ほぼランダム)。発現五分位でexposed割合均一(4.4-7.9%, JT p=0.73)。決定木の98.9%が距離のみ。多変量LR(AUC=0.629)は単変量に劣後。gene_lengthのみ有意(OR=1.001/bp, p=0.001)。**REJECTED(発現仮説)**: RNAP occupancyモデル不支持。保護は発現非依存の遺伝子座特異的構造特性 |
| `260227_H30_TSS_sequence_determinants_report.md` | 2026-02-27 | **H30: TSS近傍DNA配列特徴による保護ゾーン決定因子解析** - 1,017制御遺伝子のTSS±300bp配列解析。**AAGCCCG 5.0倍enriched**(p=1.1e-08)、TGGCCGGC 2.7倍(p=4.4e-04)。GC% +1.8pp(p=1.6e-05)。パリンドローム1.15倍(p=3.6e-03)。Combined LR CV AUC=0.712(H29 methyl距離AUC=0.917に対し)。配列寄与率~33%、残り~67%はタンパク質占有。**PARTIAL**: 配列は必要だが不十分。2層保護モデル(配列進化+タンパク質シールド) |
| `260227_H31_exposed_TF_downstream_network_report.md` | 2026-02-27 | **H31: Exposed TF下流制御ターゲットネットワーク解析** - 62 exposed TFはFIMOモチーフ非保有（レギュロン直接測定不可）。16 FIMO TFから9,635 BS(q<0.05)、5,430標的遺伝子。38/62 exposed TFがFIMO TFの標的（61.3%、shielded 64.2%と非有意差 OR=0.881, p=0.683）。標的 vs 非標的の発現差なし（\|LFC\| p=0.25/0.67）。FIMO TF間10本のcross-regulation edge（GlnRがhub）。機能エンリッチメントなし。**PARTIAL**: exposed TFは特徴付けられた制御ネットワーク内に存在するが、メチル化応答は上流TF制御とは独立 |
| `260227_H32_exposed_regulatory_module_report.md` | 2026-02-27 | **H32: 62 Exposed制御遺伝子の自己制御モジュール解析** - ゲノム上クラスタリング非有意(z=0.30, p=0.37)。E-E共発現はE-S/S-Sより高くない(p=0.26/0.99)。**協調型内共発現は極めて強い**(within rho=0.500 vs between=-0.283, p=6.6e-25, r=0.68)。4共発現モジュール(61/62遺伝子カバー): 活性化ブロック(35遺伝子)と抑制ブロック(26遺伝子)の二分法。TCS 7ペア全て非対称(exposed:shielded=1:1)、SK/RRバイアスなし(4:3, p=1.0)。オペロンペア0。**PARTIAL**: 物理的モジュールではなく「分散型メチル化応答制御層」。2つの拮抗プログラムがcoordination typeで構造化 |
| `260227_H33_neighborhood_transcriptional_impact_report.md` | 2026-02-27 | **H33: Exposed TFゲノム近傍転写影響解析** - 62 exposed TFの±10/20/50kb近傍遺伝子の発現変動検証。近傍\|LFC\|はshielded TFと差なし(MWU p=0.94)。方向一致率0.599(shielded 0.556, p=0.191)。1,000並べ替えでZ=-1.03, p=0.857。距離減衰効果なし(rho=+0.026, NS)。活性化/抑制ブロック近傍方向差は有意(chi2 p=9.2e-19)だが地理的交絡の可能性。BGC pathway-specific regulator 0件。**REJECTED**(主仮説): exposed TFは古典的cis局所制御因子ではなく、遺伝子自律的メチル化応答。**PARTIAL**: ブロック方向性分岐 |
| `260227_H34_temporal_dynamics_exposed_TF_report.md` | 2026-02-27 | **H34: 62 Exposed TFの時間的ダイナミクスと階層構造解析** - 39/62(63%)がearly responder(T1→T2で主要変動)。活性化ブロックと抑制ブロックは**同期的応答**(位相分離なし, MW p=0.459)。4協調型間のphase ratio差なし(KW p=0.638)。メチル化タイミングと発現タイミングは無相関(rho=0.136, NS)。TCS 7ペアで系統的exposed-first順序なし(Wilcoxon p=0.813)。モジュール内temporal coherence強い(rho=0.717, p=1e-19)がモジュール間予測力ゼロ。Exposed vs shielded temporal差なし(MW p=0.180)。**PARTIAL**: 同時スイッチモデル(カスケードではなく並列起動) |
| `260227_H35_AAGCCCG_exposed_TF_causal_report.md` | 2026-02-27 | **H35: AAGCCCG Exposed TF因果経路解析** - SC_RS17645→AAGCCCG 6mA→62 exposed TFの因果経路検証。**致命的発見: H30の配列モチーフ5倍enrichment(25.8%)はメチル化enrichmentに翻訳されず、わずか2/62(3.2%)のみがTSS±500bpにメチル化AAGCCCG保有**。Lost vs Never: p=0.402, r=0.383(予測と逆方向)。用量反応NS(rho=-0.117)。SC_RS17645相関差NS(p=0.346)。Bootstrap CI全て0を含む。**NOT SUPPORTED (UNTESTABLE)**: メチル化AAGCCCG部位数が少なすぎ、配列エンリッチメントはメチル化レベルに反映されない。62 exposed TFの制御機構はAAGCCCGメチル化以外 |
| `260227_H35_TF_family_functional_prediction_report.md` | 2026-02-27 | **H35: 62 Exposed TFのTFファミリー基盤機能予測解析** - 16 TFファミリーの機能アノテーション、ブロック間組成比較、ゲノム近傍文脈解析。**TetRが抑制ブロックに有意集積**(OR=0.16, p=0.028)。活性化ブロック: TCS 9 (25%, OR=4.0)、sigma因子5、WhiB 1 → 形態分化・シグナル伝達。抑制ブロック: TetR 7 (efflux隣接)、代謝TF (GntR/IclR/LacI/LysR)、DNA複製/修復 (SSB/HU/Mfd/UdgX) → 栄養増殖プログラム停止。**SUPPORTED (qualitative model)**: 2ブロック系は栄養増殖→発達転換スイッチ |
| `260227_H36_exposed_TF_conservation_report.md` | 2026-02-27 | **H36: 62 Exposed TF進化的保存性解析** - 配列組成(GC3, Nc, rare codon)、アノテーション品質(gene_name, SCO tag, hypothetical)、染色体位置(core/arm, oriC距離)、発現レベル、複合スコアで比較。**複合保存スコアに有意差なし**(p=0.276, d=0.141, ROC AUC=0.459)。17検定中有意4件はいずれも発現動態(\|LFC\| p<1e-5, constitutive OR=0.000)。SCO tag 96.8%、named product 100%。**PARTIALLY SUPPORTED**: exposed TFは配列的に同等に保存されており、系統特異的新規遺伝子ではなく古代の保存特徴。exposed/shieldedの区別はエピゲノム文脈によるもの |
| `260227_H37_coexpression_regulon_prediction_report.md` | 2026-02-27 | **H37: ゲノムワイド共発現レギュロン予測解析** - 62 exposed TFの活性化/抑制ブロック固有遺伝子(eigengene)に対するSpearman相関。活性化レギュロン28遺伝子(0.4%)、抑制レギュロン131遺伝子(1.7%)。**Jaccard=0.000の完全相互排他**。Cross-eigengene rho=-0.995。抑制ブロックは置換検定有意(Z=5.42, p=0.003)、活性化ブロックは非有意(p=0.123)。synthase 4.86倍enriched(FDR=4.3e-05)。DEG方向一致率100%。Per-TF median regulon: act=11, rep=14。**SUPPORTED (with caveats)**: 2プログラムモデル支持だが、rho=-0.995は単一発達軸の表裏。n=9の検出力限界 |

### 統合モデル（2026-02-27）

| ファイル | 日付 | 内容 |
|---------|------|------|
| `260227_gatekeeper_model_v3_synthesis_report.md` | 2026-02-27 | **Gatekeeper Model v3 包括的統合レポート** - 14探索ループ・29仮説の最終統合。中心発見: 1,055制御遺伝子のShielded/Exposed二分法（993 shielded vs 62 exposed、293bp境界、AUC=0.917）。4層モデル: Layer 1=R-M防御地理的再配置（転写制御なし）、Layer 2=制御DNA保護（25%配列進化+75%集合的プロモーター占有）、Layer 3=62 exposed遺伝子のメチル化応答カスケード。Simpson's paradox同定（H17/H19）、時間的脱抑制否定（H19/H21）、発現非依存保護（H29）等の重要な否定的結果を統合 |
| `260227_gatekeeper_model_v4_final_synthesis_report.md` | 2026-02-27 | **Gatekeeper Model v4 最終統合レポート** - 18探索ループ・36仮説の完全統合。v3(H1-H29)にH30-H36を追加。新知見: 配列寄与~33%(H30)、0/62 FIMOモチーフ・並列経路モデル(H31)、分散型制御層・活性化/抑制2ブロック(H32)、近傍効果なし・trans作用(H33)、同時スイッチ(H34)、AAGCCCG配列≠メチル化enrichment(H35-causal)、TetR抑制ブロック集積・栄養増殖→発達転換スイッチ(H35-TF)、同等保存性・エピゲノム文脈(H36)。Verdict分布: 7 supported, 15 partial, 12 rejected |
### 論文準備資料（2026-02-24）

| ファイル | 日付 | 内容 |
|---------|------|------|
| `260224_paper_materials_inventory.md` | 2026-02-24 | **論文材料棚卸し** - v2再解析後の全結果を「信頼性・新規性・効果量」で格付け、論文ストーリー方向性決定を支援 |

### エラッタ・修正（2026-02-24）

| ファイル | 日付 | 内容 |
|---------|------|------|
| `260224_locus_tag_errata_report.md` | 2026-02-24 | **TF locus_tag系統的マッピングエラーの発見と修正** - literature_tf_master.csvの37 TF中36個が誤locus_tag。GRN-TF解析、カスケードモデル等に影響。25データファイルARCHIVED。影響/非影響解析の区分、修正表、再解析優先順位を記載 |

### 包括的原稿（2026-02-06）

| ファイル | 日付 | 内容 |
|---------|------|------|
| `260224_comprehensive_manuscript_v2_report.md` | 2026-02-24 | **包括的解析報告書 v2**（論文化候補材料カタログ）- T3再シーケンシング後の全解析統合版。58+ Figure inline表示。v1→v2変更点、Key Insights 10件、ゲートキーパーモデル |

### 追加解析（2026-04-11）

| ファイル | 日付 | 内容 |
|---------|------|------|
| `260411_E1_GO_KEGG_enrichment_report.md` | 2026-04-11 | **E-1: メチル化遺伝子GO/KEGGエンリッチメント** - 3遺伝子群（GCCGGC-only=3197, AAGCCCG-only=459, Dual=448）を対象にFisher's exact test + BH補正。主要所見: 4mCはQuorum sensing (FDR=2.5e-4) / TCS (FDR=0.14) に集積。DualはSiderophore合成 (enrichment=12.2, FDR=8.7e-4) に最強集積。6mAは機能非偏在。GatekeeperモデルLayer 3（4mCによるシグナル遺伝子保護）を支持 |
| `260411_B1_window_methylation_report.md` | 2026-04-11 | **B-1: ウィンドウメチル化密度解析** - 50kbウィンドウ×3タイムポイント。4mCはT2ピーク→T3急減（二峰性）、コア/アーム比T1=2.21→T2=1.12→T3=1.32の動態変化。6mAは単調増加かつ地理的均一（core/arm≈1.0）。2メチル化システムの対照的時系列動態を可視化 |
| `260411_C2_timepoint_TSS_methylation_report.md` | 2026-04-11 | **C-2: タイムポイント別TSS保護ゾーン解析** - Jeong2016 Primary TSS（2771）×T1/T2/T3×4mC/6mA。主要所見: 4mC保護ゾーンは制御遺伝子T3で崩壊（ratio=1.046>1.0）、非制御遺伝子は維持（0.845）。6mAは制御遺伝子T2で最深保護（ratio=0.732）。2システムが異なる時間軸で制御遺伝子プロモーターを開閉するパターンを示唆 |
| `260411_discussion_points_progress.md` | 2026-04-11 | **議事録論点×解析対応 進捗記録** - 3/12・3/20議事録の6論点（A-1/A-2/A-3/B-1/C-2/E-1）全完了。各論点の主要所見・対応解析・既存ループとの関係を一覧化。未対応論点リスト付き |

---

## 注記

- 全レポートの正本は本ディレクトリ（`/reports/`）に格納
- 各解析ディレクトリ内のレポートファイルは本ディレクトリへのシンボリックリンクに置換済み
- `COMMON_PROMPT.md` はレポート生成用の共通プロンプトテンプレート
- `260203_presentation_figures/` はプレゼンテーション用Figure集（69枚）

---

*最終更新: 2026-04-11 (E1 GO/KEGG, B1 ウィンドウメチル化, C2 タイムポイント別TSS 追加)*
