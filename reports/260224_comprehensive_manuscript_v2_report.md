# M145 RNA-seq/Methylome統合解析 包括的解析報告書 v2

**作成日**: 2026-02-24
**プロジェクト**: *Streptomyces coelicolor* A3(2) M145 トランスクリプトーム-エピゲノム統合解析
**前版**: 260204_summary_manuscript.md (v1, weighted MIN_REPS=2)
**本版の位置付け**: T3追加シーケンシングデータ（2026-02-24統合）による**全解析再実行版**。論文化に向けた全候補材料の包括的カタログ。

### v1 → v2 主要変更点

| 項目 | v1 (旧) | v2 (本版) |
|------|---------|----------|
| コンセンサス手法 | カバレッジ加重平均 | **単純平均（unweighted）** |
| MIN_REPS | 2 | **3（全3レプリケート一致を要求）** |
| サンプル3-2カバレッジ | 14.3x | **57.6x (+303%)** |
| High-confidence sites | 11,778 | **12,429 (+5.5%)** |
| 4mC T2vsT1相関 | r=0.172 (p=2.3e-5) | **r=0.125 (p=0.002) — 維持** |
| 6mA T2vsT1相関 | r=0.169 (p=2.1e-5) | **r=-0.038 (p=0.35) — 消失** |
| FDR補正（67テスト） | 未実施 | **1/67のみ有意（6mA T3vsT1 gene body, p_adj=0.049）** |

---

## 目次

**I. 序論**
  研究背景と目的 / 実験デザイン / 解析パイプライン概要

**II. Results**
  1. 培養段階依存的トランスクリプトームリプログラミング
  2. BGC発現動態と二段階活性化モデル
  3. メチル化ランドスケープとモチーフの同定
  4. メチル化-発現相関の統合解析
  5. 擬似相関検証とFDR多重検定補正
  6. 新規AAGCCCG R-M系の発見と検証
  7. DMG機能エンリッチメント解析
  8. TF Binding Siteメチル化変動解析
  9. GRN-TFメチル化統合解析とAAGCCCGカスケード
  10. m4C/m5C二重シトシン修飾系

**III. Discussion**
  提案モデル / 先行研究との比較 / Key Insights / 知見の堅牢性

**IV. 研究の限界と今後**

**V. 付録**
  全Figure一覧 / 使用ソフトウェア

---

## 研究背景と目的

### 背景

- **放線菌** (*Streptomyces*属) は抗生物質、免疫抑制剤など有用二次代謝産物の主要な生産者
- *S. coelicolor* A3(2) M145は放線菌のモデル生物として広く研究されている
- 二次代謝の活性化は培養段階に強く依存（成長期→定常期の移行）
- エピジェネティック制御（DNAメチル化）の関与は未解明

### 先行研究

- **Pisciotta et al. (2018)**: SCO1731（m5C MTase）KO株で形態分化と二次代謝に影響
- **Pisciotta et al. (2023)**: BS-seqによりM145の5mCメチロームを報告（3,360サイト、GGCmCGG/GCCmCGモチーフ）
- **Fang et al. (2022)**: *S. roseosporus* L30でSroLm3（4mC MTase）がダプトマイシン生合成を制御
- しかし、6mA/4mCの**アデニン/シトシンメチル化**と転写制御のTSS基準統合解析は未実施

### 目的

1. M145の3培養段階（T1, T2, T3）における全ゲノム的転写変動の解明
2. DNAメチル化（6mA, 4mC）と遺伝子発現の関連をTSS基準で解析
3. 新規エピジェネティック制御機構の同定と**多重検定補正による検証**

---

## 実験デザイン

| 条件 | タイムポイント | 生物学的レプリケート |
|------|--------------|-------------------|
| M145_1 (T1) | 初期（増殖期） | n = 3 |
| M145_2 (T2) | 中期（移行期） | n = 3 |
| M145_3 (T3) | 後期（定常期） | n = 3 |

**合計**: 9サンプル × 2プラットフォーム（Illumina RNA-seq + Nanopore メチローム）

### T3追加シーケンシング（v2新規）

| サンプル | 旧カバレッジ | 新カバレッジ | 改善率 |
|---------|-----------|-----------|--------|
| 3-2 | 14.3x | **57.6x** | +303% |
| 3-4 | 62.9x | **75.4x** | +20% |

---

## 解析パイプライン概要

### トランスクリプトーム解析

```
Raw FASTQ (PE 150bp)
    ↓ FastQC v0.12.1
    ↓ fastp v1.1.0 (アダプター除去)
    ↓ HISAT2 v2.2.1 (--no-spliced-alignment)
    ↓ SAMtools v1.21
    ↓ featureCounts v2.1.1 (-s 0 unstranded)
    ↓ DESeq2 v1.46.0 + apeglm v1.28.0
```

### メチローム解析（v2更新）

```
Nanopore BAM (9サンプル、3-2/3-4は追加データ統合済み)
    ↓ modkit pileup (修飾塩基コール)
    ↓ コンセンサス: 単純平均, MIN_REPS=3, MIN_COVERAGE=5, MIN_MOD_FREQ≥50%
    ↓ 実験的TSS統合 (Jeong et al. 2016 dRNA-seq)
    ↓ TSS基準プロモーター領域定義 (-300〜+50bp)
```

---

---

# Results

---

## 1. 培養段階依存的トランスクリプトームリプログラミング

### 【Fig. 1A】PCAによる培養段階間の転写プロファイル分離

![PCA_M145](../04_deseq2/analysis/04_deseq2_260128_v1/figures/PCA_M145.svg)

**方法**: DESeq2 rlog変換後、plotPCA()。7,646遺伝子、9サンプル。

**結果**:
- PC1（83%）がT1 vs T2/T3を分離 → 増殖期→移行期が最大の転写変動
- PC2（16%）がT2 vs T3を分離 → 二次的だが明確な変化
- レプリケート間の高い再現性

### 【Fig. S2】サンプル間階層的クラスタリング

![sample_distance_heatmap_M145](../04_deseq2/analysis/04_deseq2_260128_v1/figures/sample_distance_heatmap_M145.svg)

**結果**: T2-T3が最も類似、T1が最も孤立。増殖期→移行期の遷移が支配的イベント。

### DEG数サマリー

| 比較 | Up | Down | 合計 | 全遺伝子中割合 |
|------|-----|------|------|--------------|
| T2 vs T1 | 2,563 | 2,697 | 5,260 | **69%** |
| T3 vs T1 | 3,290 | 2,898 | 6,188 | **81%** |
| T3 vs T2 | 2,639 | 2,904 | 5,543 | **73%** |

### 【Fig. S3a-c】Volcanoプロット

![volcano_T2vsT1](../04_deseq2/analysis/04_deseq2_260128_v1/figures/volcano_M145_2_vs_1.svg)

![volcano_T3vsT1](../04_deseq2/analysis/04_deseq2_260128_v1/figures/volcano_M145_3_vs_1.svg)

![volcano_T3vsT2](../04_deseq2/analysis/04_deseq2_260128_v1/figures/volcano_M145_3_vs_2.svg)

---

## 2. BGC発現動態と二段階活性化モデル

### BGC平均発現の経時変化

| BGC | T1 | T2 | T3 | FC(T2/T1) | FC(T3/T1) |
|-----|-----|-----|-----|-----------|-----------|
| act | 65 | 99 | 9,491 | 1.5x | **145x** |
| red | 65 | 682 | 966 | 10.4x | **14.8x** |
| cda | 57 | 3,309 | 3,085 | **58.4x** | 54.4x |
| cpk | 52 | 9,465 | 6,144 | **181.7x** | 118x |

### 【Fig. 1D】BGC経時変化

![BGC_timecourse](../06_BGC_dynamics/analysis/06_BGC_dynamics_260128_v1/figures/BGC_timecourse_lineplot_M145.svg)

**二段階活性化モデル**:
1. **Phase I（T1→T2）**: red, cda, cpkが10〜180倍に急増。actはわずか
2. **Phase II（T2→T3）**: actのみが追加的に96倍増加

### 【Fig. 1E】Act BGCヒートマップ

![act_heatmap](../11_epigenome_integration/analysis/02_publication_figures/act_bgc_expression_heatmap.svg)

- T3 vs T1で全22/22遺伝子が有意上昇（平均LFC = +8.34）
- T2 vs T1で13/22が前段階的誘導済み（平均LFC = +0.83）
- **Act BGCは「Phase II爆発型」活性化**

### 【Fig. 1F】Red BGCヒートマップ

![red_heatmap](../11_epigenome_integration/analysis/02_publication_figures/red_bgc_expression_heatmap.svg)

- 21/22遺伝子がT3 vs T1で有意上昇
- T2 vs T1で既にLFC > 2の遺伝子が多数 → Phase I早期誘導
- プロモーターにメチル化サイトをほぼ持たない
- **Red BGCは「Phase I早期誘導・持続型」**

---

## 3. メチル化ランドスケープとモチーフの同定

### High-Confidence Sites（v2更新）

| メチル化型 | T1 | T2 | T3 | 合計 | ユニーク位置 |
|-----------|------|------|------|------|-----------|
| 6mA | 1,934 | 2,120 | 2,295 | 6,349 | 3,248 |
| 4mC | 1,987 | 2,446 | 1,647 | 6,080 | 2,693 |
| **合計** | **3,921** | **4,566** | **3,942** | **12,429** | **5,941** |

**v1との比較**: 4mC T3が1,077→1,647（+52.9%）と最大改善。3-2カバレッジ改善の直接的効果。

### 【Fig. 2A】モチーフスコアリングヒートマップ

![motif_scoring](../11_epigenome_integration/analysis/23_expanded_motif_search/figures/integration_scoring_heatmap.svg)

### 【Fig. 2B】モチーフ帰属円グラフ

![motif_pie](../11_epigenome_integration/analysis/23_expanded_motif_search/figures/motif_attribution_piechart.svg)

### 同定されたメチル化モチーフ

| モチーフ | 修飾型 | サイト数 | スコア | 分類 | 候補MTase |
|---------|-------|---------|--------|------|----------|
| **AAGCCCG** | 4mC+6mA (DUAL) | 973+360 | **10/10** | HIGH-CONFIDENCE NOVEL | SC_RS17645 |
| CCGG/TGGCCGGC | 4mC | 2,034 (75.9%) | 4/10 | ESTABLISHED | SC_RS19770/SC_RS36410 |
| CCGKCA | 6mA | 153 (4.8%) | 7/10 | CANDIDATE | SC_RS28835/SC_RS35335 (BREX-2) |
| GATC | 6mA | 37 (1.2%) | 4/10 | ESTABLISHED | Unknown |
| GAACCGG | 6mA | 61 (1.9%) | 4/10 | ESTABLISHED/INSUFFICIENT | Unknown orphan |
| CGGCAACC | 6mA | 57 (1.8%) | 3/10 | ESTABLISHED/INSUFFICIENT | Unknown orphan |

### 【Fig. 2C】REBASE属間保存性

![rebase_conservation](../11_epigenome_integration/analysis/23_expanded_motif_search/figures/rebase_refseq_conservation.svg)

**AAGCCCG**: REBASE *Streptomyces* 82種中0種で報告 — **属レベルで完全に新規**

### 【Fig. 2D】O/E比較

![oe_ratio](../11_epigenome_integration/analysis/23_expanded_motif_search/figures/oe_ratio_comparison.svg)

### 【Fig. 2E】モチーフ時系列動態

![motif_temporal](../11_epigenome_integration/analysis/02_publication_figures/motif_temporal_dynamics.svg)

**時系列パターン**:
- 4mCサイト総数: T1=1,987 → T2=2,446（+23%）→ T3=1,647（-17%）
- 6mAサイト総数: T1=1,934 → T2=2,120 → T3=2,295（+19%）
- **4mCはT2でピーク後減少、6mAは単調増加** → 異なる動態パターン

### 【Fig. 2F】タイムポイント別MEMEロゴ

![meme_logos](../11_epigenome_integration/analysis/02_publication_figures/timepoint_meme_logos.svg)

### 【Fig. 2G】モチーフ動態（Fold Change）

![motif_foldchange](../11_epigenome_integration/analysis/02_publication_figures/timepoint_motif_foldchange.svg)

---

## 4. メチル化-発現相関の統合解析

### ★ v2における最重要更新事項

**6mA T2vsT1相関の消失**:
| 比較 | v1 (r, p) | v2 (r, p) | 判定 |
|------|-----------|-----------|------|
| 4mC T2vsT1 | r=0.172, p=2.3e-5 | **r=0.125, p=0.002** | **維持（robust）** |
| 6mA T2vsT1 | r=0.169, p=2.1e-5 | **r=-0.038, p=0.35** | **消失（artifact）** |
| 6mA T3vsT1 | r=-0.176, p=4.3e-6 | **r=-0.093, p=0.016** | 弱化 |
| 4mC T3vsT1 | r=0.061, p=0.15 | r=-0.003, p=0.94 | 非有意（維持） |

### 【Fig. 3A】TSS周辺メタジーンプロファイル

![metagene](../11_epigenome_integration/analysis/18_tss_analyses/B1_metagene_methylation_profile.png)

**方法**: 実験的TSS (Jeong et al. 2016) を基準に±2000bpのメチル化密度を算出

**結果**: σ因子 -10 box領域（TSS -10bp付近）でメチル化が有意に枯渇（p=1.73e-12）

### 【Fig. 3B】TSS距離帯別相関ヒートマップ

![correlation_heatmap](../11_epigenome_integration/analysis/18_tss_analyses/C1_distance_correlation_heatmap.png)

### 【Fig. 3C】6パネル散布図（カテゴリ別統計）

![6panel_correlation](../11_epigenome_integration/analysis/02_publication_figures/Fig1_methylation_expression_correlation_6panel.svg)

**カテゴリ別統計**:
- **4mC T2vsT1**: Gained遺伝子の中央値LFC = +0.569, Lost = -0.417（Mann-Whitney p=1.1e-4, **FDR有意**）
- **6mA T2vsT1**: Gained/Lost間で有意差なし（p=0.027, FDR後非有意）

### 【Fig. 3D】協調変動パターン

![concordance](../11_epigenome_integration/analysis/02_publication_figures/concordance_patterns.svg)

**Concordance Summary**:

| 比較 | 協調 | 不協調 | その他 | 総数 |
|------|------|--------|--------|------|
| T2 vs T1 | 148 (18.3%) | 128 (15.8%) | 532 (65.8%) | 808 |
| T3 vs T1 | 205 (72.4%) | 78 (27.6%) | 0 | 283 |
| T3 vs T2 | 136 (73.9%) | 48 (26.1%) | 0 | 184 |

### FDR多重検定補正（67テスト、v2新規）

67個の独立したSpearman相関テスト（6修飾型×距離帯×比較）にBH-FDR補正を適用。

| テスト | p_raw | p_adj | 有意 |
|--------|-------|-------|------|
| 6mA T3vsT1 gene body | 0.0007 | **0.049** | **Yes** |
| 4mC T3vsT2 mid gene body | 0.003 | 0.094 | No |
| 他65テスト | — | >0.10 | No |

**結論**: FDR補正後、**67テスト中わずか1件のみ有意**。メチル化-発現相関は全般的に弱く、直接的な転写制御よりも間接的・位置依存的な関係。

---

## 5. 擬似相関検証

### 【Fig. 5A】並べ替え検定

![permutation](../11_epigenome_integration/analysis/09_spurious_validation/permutation_test_results.png)

### 検証結果サマリー

| 比較 | 観測r | 並べ替えp | Bootstrap CI | ゼロ含む？ | 判定 |
|------|-------|----------|-------------|----------|------|
| 4mC T2vsT1 | 0.125 | **0.002** | [0.046, 0.194] | No | **Robust** |
| 6mA T2vsT1 | -0.038 | 0.349 | [-0.114, 0.047] | Yes | Non-significant |
| 4mC T3vsT1 | -0.003 | 0.937 | [-0.088, 0.089] | Yes | Non-significant |
| 6mA T3vsT1 | -0.093 | **0.015** | [-0.168, -0.012] | No | Significant |

### 【Fig. 5B】偏相関（GC含量補正）

![partial_correlation](../11_epigenome_integration/analysis/09_spurious_validation/partial_correlation_comparison.png)

**4mC T2vsT1**: 偏相関r=0.139（p=0.0005）→ GC含量補正後に**増強**。真の生物学的シグナル。

---

## 6. 新規AAGCCCG R-M系の発見と検証

### AAGCCCG の特徴

| 属性 | 値 |
|------|-----|
| 修飾型 | **二重修飾（4mC + 6mA）** |
| 4mCサイト | 973（36.3%） |
| 6mAサイト | 360（11.2%） |
| REBASE報告 | **0/82種（属レベル新規）** |
| RefSeq保存 | 低い（系統特異的） |
| 候補MTase | **SC_RS17645**（N-6 DNA methylase, HsdM型） |
| スコア | **10/10 (HIGH-CONFIDENCE NOVEL)** |
| ゲノムO/E | 0.65（回避傾向） |

### 【Fig. 6A】AAGCCCGプロモーター解析

![aagcccg_promoter](../11_epigenome_integration/analysis/14_aagcccg_promoter_analysis/aagcccg_promoter_analysis.png)

**結果**: 協調変動遺伝子のプロモーターでAAGCCCGが有意に濃縮（**OR=12.89, p=1.1e-107**）

### SC_RS17645（SCO3104）の特性

- **配列ホモロジー**: Type I R-M system HsdMサブユニット（N-6 adenine methylase domain）
- **発現変動**: T2 vs T1で log2FC = -2.19（**発現低下**）
- **カスケードモデル**: SC_RS17645発現低下 → AAGCCCG 6mA減少 → redZプロモーター脱メチル化 → redZ発現低下

### 【Fig. S-SC_RS17645】構造ホモロジー

![sc_rs17645_homology](../11_epigenome_integration/analysis/13_sc_rs17645_analysis/fig_sc_rs17645_homology_composite.svg)

---

## 7. DMG機能エンリッチメント解析（v2新規）

### DMG（Differentially Methylated Gene）サマリー

| 比較 | hyper | hypo | 合計 |
|------|-------|------|------|
| T2vsT1 | 336 | 237 | **566** |
| T3vsT1 | 281 | 348 | **619** |
| T3vsT2 | 250 | 411 | **648** |

### 【Fig. 7A】DMGサマリー

![dmg_summary](../14_DMG_functional_enrichment/analysis/14_DMG_enrichment_260206_v1/figures/DMG_summary_counts.svg)

### 【Fig. 7B】KEGG DMG vs DEG ヒートマップ

![kegg_heatmap](../14_DMG_functional_enrichment/analysis/14_DMG_enrichment_260206_v1/figures/KEGG_DMG_vs_DEG_heatmap.svg)

### 【Fig. 7C】COG分布

![cog_t2vst1](../14_DMG_functional_enrichment/analysis/14_DMG_enrichment_260206_v1/figures/COG_distribution_T2vsT1.svg)

### 結果

DMGはKEGG、GO、COGいずれのデータベースにおいても**有意な機能エンリッチメントを示さなかった**。

**「位置依存・機能非依存」モデル**:
- メチル化はゲノム上の特定モチーフ位置に依存して発生
- 遺伝子機能による選択はない
- メチル化変動が特定の生物学的パスウェイを標的としているわけではない
- エピジェネティック制御は**間接的**（特定のTF binding siteを介した経路依存的制御）

### DEG-DMG重複検定

| 比較 | Overlap | OR | p値 | 有意 |
|------|---------|-----|-----|------|
| T2vsT1 | 267 | 1.00 | 1.0 | No |
| **T3vsT1** | **407** | **1.35** | **5.6e-4** | **Yes** |
| T3vsT2 | 301 | 1.16 | 0.069 | No |

### 【Fig. 7D-F】DEG-DMG Venn図

![venn_t2vst1](../11_epigenome_integration/analysis/03_overlap_analysis/venn_DEG_DMG_T2vsT1.png)

![venn_t3vst1](../11_epigenome_integration/analysis/03_overlap_analysis/venn_DEG_DMG_T3vsT1.png)

![venn_t3vst2](../11_epigenome_integration/analysis/03_overlap_analysis/venn_DEG_DMG_T3vsT2.png)

### 【Fig. 7G】4mC相関ヒートマップ（Selectivity解析）

![selectivity_4mC](../14_DMG_functional_enrichment/analysis/14_DMG_selectivity_260206_v1/figures/A2_correlation_heatmap_4mC.svg)

### 【Fig. 7H】COG: プロモーター vs 遺伝子本体

![cog_promoter_body](../14_DMG_functional_enrichment/analysis/14_DMG_selectivity_260206_v1/figures/A3_COG_promoter_vs_body.svg)

### 【Fig. 7I】TF/BGC DMG頻度

![tf_bgc_dmg](../14_DMG_functional_enrichment/analysis/14_DMG_selectivity_260206_v1/figures/A5_TF_BGC_DMG_frequency.svg)

---

## 8. TF Binding Siteメチル化変動解析（v2新規）

### 空間的重複解析

| Window | BS数 | BS+methyl | fold | p値 |
|--------|------|-----------|------|-----|
| Direct overlap | 782 | 63 (8.1%) | **0.659** | **6.8e-5** |
| ±50bp | 782 | 100 (12.8%) | 0.746 | 3.6e-4 |
| ±200bp | 782 | 210 (26.9%) | 0.810 | 2.2e-4 |

**核心的発見**: TF binding siteでメチル化が有意に**枯渇**（fold=0.66, p<0.001）。メチル化はTFの結合を妨害する位置を回避する傾向。

### 【Fig. 8A】BS周辺メチル化メタジーンプロファイル

![bs_metagene](../13_TF_binding-site/analysis/02_TF_BS_methylation_260207_v1/figures/F1_metagene_methylation_around_BS.svg)

### 【Fig. 8B】BS近傍メチル化の時系列ヒートマップ

![bs_temporal](../13_TF_binding-site/analysis/02_TF_BS_methylation_260207_v1/figures/F2_BS_methylation_temporal_heatmap.svg)

**時系列パターン**: BS近傍メチル化サイトの**75.2%がLost_T2T3パターン**（T1で検出→T2/T3で消失）。移行期でのTF binding site周辺の脱メチル化が広範に発生。

### 【Fig. 8C】BSメチル化変動 vs ターゲット発現変動

![bs_expression](../13_TF_binding-site/analysis/02_TF_BS_methylation_260207_v1/figures/F3_BS_methylation_vs_expression.svg)

### 【Fig. 8D】TFファミリー別BSメチル化率

![tf_family](../13_TF_binding-site/analysis/02_TF_BS_methylation_260207_v1/figures/F4_TF_family_BS_methylation_rate.svg)

---

## 9. GRN-TFメチル化統合解析とAAGCCCGカスケード

### TFプロモーターメチル化（37 TFスキャン）

37個の文献既知GRN TFのプロモーターを検査した結果、**わずか5個のみがプロモーターメチル化を保持**。

### 【Fig. 9A】制御カスケードと協調変動

![cascade](../11_epigenome_integration/analysis/12_grn_tf_methylation/tf_promoter_methylation/regulatory_cascade_methylation.png)

### 【Fig. 9B】TF-ターゲットTriadネットワーク

![triad](../11_epigenome_integration/analysis/12_grn_tf_methylation/triad_analysis/triad_network_visualization.png)

### 【Fig. 9C】BGC発現サマリー

![bgc_summary](../11_epigenome_integration/analysis/12_grn_tf_methylation/triad_analysis/bgc_expression_summary.png)

### redZ: 唯一の協調変動TF

| 属性 | 値 |
|------|-----|
| 遺伝子 | SC_RS27300 (redZ/SCO5881) |
| 機能 | CSR (Cluster-Specific Regulator), Red BGC |
| メチル化変動 | **6mA Lost** (T1→T2でプロモーター脱メチル化) |
| 発現変動 | **DOWN** (log2FC = -2.25, T2 vs T1) |
| 協調性 | **Positive** (Lost methylation + Down expression) |

### 【Fig. 9D】エピジェネティック制御ネットワーク

![network](../11_epigenome_integration/analysis/12_grn_tf_methylation/network_visualization/comprehensive_epigenetic_network.png)

### 【Fig. 9E】カスケード模式図

![cascade_diagram](../11_epigenome_integration/analysis/12_grn_tf_methylation/network_visualization/epigenetic_cascade_diagram.png)

### 提案するAAGCCCGカスケードモデル

```
環境シグナル（T2タイムポイント）
    ↓
SC_RS17645 (N-6 MTase) 発現低下 (log2FC = -2.19)
    ↓
AAGCCCG 6mA メチル化の全体的減少
    ↓
redZ プロモーター脱メチル化
    ↓
redZ 発現低下 (log2FC = -2.25)
    ↓
redD + Red BGC（複合制御）
```

**パラドックス**: redZが低下にもかかわらずredDとRed BGCは上昇 → AbsA2抑圧解除または他の活性化因子による補償

---

## 10. m4C/m5C二重シトシン修飾系

### 【Fig. 10A】4mC/5mC Competition

![4mc_5mc](../11_epigenome_integration/analysis/22_4mC_5mC_competition/fig_4mC_5mC_competition.svg)

**発見**: 100%の4mCコンセンサスモチーフがGGCCGGコンテクストに存在 → Dcm-likeパラドックス

---

## 追加Figure: 個別遺伝子可視化

### 【Fig. S-TD】経時トラジェクトリー（全遺伝子概要）

![temporal_all](../11_epigenome_integration/analysis/05_temporal_dynamics/temporal_dynamics_all_genes.png)

### 【Fig. S-TD-RamR】RamR経時トラジェクトリー

![temporal_ramr](../11_epigenome_integration/analysis/05_temporal_dynamics/temporal_dynamics_RamR.png)

### 【Fig. S-TD-NsdB】NsdB経時トラジェクトリー

![temporal_nsdb](../11_epigenome_integration/analysis/05_temporal_dynamics/temporal_dynamics_NsdB.png)

### 【Fig. S-TD-RedOx】Red oxygenase経時トラジェクトリー

![temporal_red](../11_epigenome_integration/analysis/05_temporal_dynamics/temporal_dynamics_Red_oxygenase.png)

### 【Fig. S-TD-ActNmrA】Act NmrA経時トラジェクトリー

![temporal_act](../11_epigenome_integration/analysis/05_temporal_dynamics/temporal_dynamics_Act_NmrA.png)

### 【Fig. S-TD-Cpk】Cpk carboxylase経時トラジェクトリー

![temporal_cpk](../11_epigenome_integration/analysis/05_temporal_dynamics/temporal_dynamics_Cpk_carboxylase.png)

### 【Fig. S-TD-CDA】CDA trpC経時トラジェクトリー

![temporal_cda](../11_epigenome_integration/analysis/05_temporal_dynamics/temporal_dynamics_CDA_trpC.png)

### 【Fig. S-TD-Traj】メチル化-発現トラジェクトリー

![trajectory](../11_epigenome_integration/analysis/05_temporal_dynamics/methylation_expression_trajectories.png)

---

## 追加Figure: マルチオミクストラック

### 【Fig. S-MT-RamR】RamRマルチオミクストラック

![mt_ramr](../11_epigenome_integration/analysis/04_multiomics_tracks/multiomics_track_RamR.svg)

### 【Fig. S-MT-NsdB】NsdBマルチオミクストラック

![mt_nsdb](../11_epigenome_integration/analysis/04_multiomics_tracks/multiomics_track_NsdB.svg)

### 【Fig. S-MT-Red】Red SCO5897マルチオミクストラック

![mt_red](../11_epigenome_integration/analysis/04_multiomics_tracks/multiomics_track_Red_SCO5897.svg)

### 【Fig. S-MT-Act】Act SCO5079マルチオミクストラック

![mt_act](../11_epigenome_integration/analysis/04_multiomics_tracks/multiomics_track_Act_SCO5079.svg)

### 【Fig. S-MT-Cpk】Cpk SCO6284マルチオミクストラック

![mt_cpk](../11_epigenome_integration/analysis/04_multiomics_tracks/multiomics_track_Cpk_SCO6284.svg)

---

## 追加Figure: プロモーター構造模式図

### 【Fig. S-PS】プロモーター構造サマリー

![promoter_summary](../11_epigenome_integration/analysis/06_motif_schematics/promoter_architecture_summary.png)

### 【Fig. S-PS-All】全遺伝子プロモーター模式図

![motif_all](../11_epigenome_integration/analysis/06_motif_schematics/motif_schematic_all_genes.png)

---

## 追加Figure: TSS関連解析

### 【Fig. S-TSS-Null】-10 box帰無モデル

![null_minus10](../11_epigenome_integration/analysis/18_tss_analyses/null_model_minus10_box.png)

### 【Fig. S-TSS-MTase】MTaseメチル化動態

![mtase_dynamics](../11_epigenome_integration/analysis/18_tss_analyses/mtase_methylation_dynamics.png)

### 【Fig. S-TSS-SARP】SARPゾーン空間解析

![sarp_spatial](../11_epigenome_integration/analysis/18_tss_analyses/sarp_methylation_mechanism_schematic.png)

### 【Fig. S-TSS-TF】TFメチル化エンリッチメント

![tf_enrichment](../11_epigenome_integration/analysis/18_tss_analyses/tf_methylation_enrichment.png)

### 【Fig. S-TSS-Temporal】Gained/Lostモチーフバイアス

![temporal_gained_lost](../11_epigenome_integration/analysis/18_tss_analyses/temporal_gained_lost_motif_barplot.png)

---

## 追加Figure: Act BGCエピジェネティック解析

### 【Fig. S-ActBGC】Act vs Red BGCエピジェネティック比較

![act_bgc_epigenetic](../11_epigenome_integration/analysis/15_act_bgc_epigenetic/act_bgc_epigenetic_analysis.png)

**結果**: Act/Red BGC遺伝子はプロモーターにメチル化サイトをほぼ持たない → 直接的エピジェネティック制御なし。制御はredZを介した間接的経路。

---

## 追加Figure: T3協調遺伝子比較

### 【Fig. S-T3Coord】T3協調遺伝子比較

![t3_coordinated](../11_epigenome_integration/analysis/16_t3_coordinated/t3_coordinated_comparison.png)

---

## 追加Figure: 属間比較解析

### 【Fig. S-Genus】属全種メチル化モチーフ保存性

![genus_conservation](../11_epigenome_integration/analysis/21_genuswide_motif_conservation/fig_motif_conservation_composite.svg)

### 【Fig. S-BGCReg】BGC制御因子メチル化ランドスケープ

![bgc_regulator](../11_epigenome_integration/analysis/19_bgc_regulator_overview/bgc_regulator_overview.svg)

### 【Fig. S-Cascade】AAGCCCGカスケード解析

![aagcccg_cascade](../11_epigenome_integration/analysis/19_bgc_regulator_overview/aagcccg_cascade_analysis.svg)

---

## 追加Figure: 論文用総合Figure

### 【Paper Fig. 1】研究概要

![paper_fig1](../11_epigenome_integration/analysis/17_paper_figures/Figure1_overview.svg)

### 【Paper Fig. 2】AAGCCCG R-M系

![paper_fig2](../11_epigenome_integration/analysis/17_paper_figures/Figure2_AAGCCCG_system.svg)

### 【Paper Fig. 3】GRNメチル化

![paper_fig3](../11_epigenome_integration/analysis/17_paper_figures/Figure3_GRN_methylation.svg)

### 【Paper Fig. 4】Act vs Red BGC

![paper_fig4](../11_epigenome_integration/analysis/17_paper_figures/Figure4_Act_vs_Red.svg)

---

---

# Discussion

## 提案するモデル: エピゲノム「ゲートキーパー」モデル

本研究の結果を統合すると、M145におけるDNAメチル化は従来想定されていた「マスタースイッチ」型の直接的転写制御ではなく、以下の**3層構造のゲートキーパーモデル**で記述される。

### 層1: 位置依存的・機能非依存的メチル化（「ランドスケープ層」）

- DMGはKEGG/GO/COGいずれでも有意な機能エンリッチメントを示さない
- メチル化はゲノム上のモチーフ配列（CCGG, AAGCCCG等）の位置に依存
- ゲノム全体に散在する約5,900の高信頼度メチル化位置が存在

### 層2: TF binding site回避（「保護層」）

- TF BSでメチル化が有意に枯渇（fold=0.66, p=6.8e-5）
- 75%のBS近傍メチル化がT1→T2で消失 → 移行期での脱メチル化
- メチル化はTF結合を妨害しうる位置を選択的に「回避」または「解除」

### 層3: 特定カスケードによる間接制御（「シグナル層」）

- SC_RS17645（N-6 MTase）→ AAGCCCG → redZプロモーター → Red BGC
- 37個のGRN TFのうち**わずか1個（redZ）のみが協調変動**
- エピジェネティック制御は汎用的ではなく、特定のカスケードに限定

## Key Insights

### Impact: HIGH — 本研究の新規性

| # | 発見 | インパクト | エビデンスの強さ |
|---|------|----------|----------------|
| 1 | **AAGCCCG二重修飾モチーフ（4mC+6mA）** | REBASE 82種中0種で報告。*Streptomyces*属で初の同一配列二重修飾 | スコア10/10 |
| 2 | **SC_RS17645 = Type I HsdM型 N-6 MTase** | AAGCCCG修飾の候補酵素同定。配列・構造ホモロジーで確認 | BLAST + Pfam |
| 3 | **4mC T2vsT1相関のrobustness** | 擬似相関検証（並べ替え p=0.002, Bootstrap CI exclude 0, 偏相関で増強） | 4重検証 |
| 4 | **6mA T2vsT1相関の消失** | v1で報告されたr=0.169がv2でr=-0.038に。MIN_REPS=3で再現せず | **重要な否定的結果** |
| 5 | **FDR補正後67テスト中1件のみ有意** | メチル化-発現相関は全般的に弱い。直接制御ではなく間接制御 | BH-FDR |
| 6 | **TF BS methylation depletion** | fold=0.66 (p=6.8e-5)。メチル化はTF結合部位を回避 | Tier 1 BSベース |
| 7 | **DMG機能非選択性** | DMGはKEGG/GO/COGいずれでもenrichmentなし。「位置依存・機能非依存」 | 3データベース |
| 8 | **redZ = 唯一の協調変動TF** | 37 GRN TF中唯一。6mA Lost + expression DOWN | Coordinated screen |
| 9 | **CCGKCA = BREX-2候補** | Score 7/10。SC_RS28835/SC_RS35335がcandidate MTase | MEME + REBASE |
| 10 | **エピゲノム「ゲートキーパー」モデル** | 3層構造（ランドスケープ/保護/シグナル）の統合モデル | 全解析統合 |

### Impact: MEDIUM — 知見の堅牢性

| 検証 | 結果 | 意義 |
|------|------|------|
| 並べ替え検定（N=1000） | 4mC T2vsT1: p=0.002 | 偶然では説明不能 |
| 偏相関（GC含量補正） | r: 0.125 → 0.139 | 交絡因子ではない |
| Bootstrap CI (N=1000) | [0.046, 0.194] | ゼロを含まない |
| 負の対照群 | r ≈ 0 | 陽性結果の特異性確認 |
| FDR多重検定 | 67テスト中1件有意 | 過度な多重比較補正でも一部残存 |

---

# 研究の限界と今後

## 研究の限界

1. **相関研究の本質的限界**: メチル化と発現の相関は因果関係を証明しない。MTase KO実験が必要
2. **タイムポイント数**: 3点のみ。連続的な動態追跡には不十分
3. **生物学的レプリケート**: n=3。統計的検出力に限界
4. **STREME/MEMEは旧データで実行**: モチーフ発見自体は Feb 3 実行、ただし下流の統合スコアリングはFeb 24更新データで再評価済み
5. **SC_RS17645の機能未確認**: in silico同定のみ。KOまたはCRISPRによる実験的検証が必要
6. **TF binding siteの網羅性**: Tier 1 BSは790件のみ。RegPrecise/ZorroAranda/FIMOの範囲に限定

## 今後の検証実験（優先順位順）

1. **SC_RS17645 (N-6 MTase) の遺伝学的検証**
   - CRISPRi/KOによるAAGCCCGメチル化への影響
   - redZ/Red BGC発現への影響

2. **AAGCCCG二重修飾の生化学的検証**
   - 精製SC_RS17645による in vitro メチル化アッセイ
   - 修飾配列特異性の確認

3. **redZプロモーターのメチル化感受性**
   - メチル化/非メチル化プロモーターでの in vitro 転写活性比較
   - EMSA (Electrophoretic Mobility Shift Assay) によるTF結合へのメチル化影響

4. **時系列の高解像度化**
   - 5-7点のタイムコースで動態を追跡
   - 特にT1→T2遷移の詳細な分解

---

# 付録

## 使用ソフトウェア

| ソフトウェア | バージョン | 用途 |
|-------------|----------|------|
| FastQC | 0.12.1 | 品質管理 |
| fastp | 1.1.0 | アダプター除去 |
| HISAT2 | 2.2.1 | アライメント |
| SAMtools | 1.21 | BAM操作 |
| featureCounts | 2.1.1 | カウント |
| DESeq2 | 1.46.0 | 差次的発現解析 |
| apeglm | 1.28.0 | LFC shrinkage |
| modkit | — | 修飾塩基コール |
| MEME Suite | — | モチーフ発見 |
| Python | 3.11 | 統合解析・可視化 |
| R | 4.4.2 | 統計解析 |

## 全Figure一覧

### 本文Figure

| # | Figure ID | セクション | 内容 |
|---|-----------|----------|------|
| 1 | Fig. 1A | Sec. 1 | PCA |
| 2 | Fig. S2 | Sec. 1 | サンプル距離ヒートマップ |
| 3 | Fig. S3a-c | Sec. 1 | Volcanoプロット ×3 |
| 4 | Fig. 1D | Sec. 2 | BGC経時変化 |
| 5 | Fig. 1E | Sec. 2 | Act BGCヒートマップ |
| 6 | Fig. 1F | Sec. 2 | Red BGCヒートマップ |
| 7 | Fig. 2A | Sec. 3 | モチーフスコアリング |
| 8 | Fig. 2B | Sec. 3 | モチーフ帰属円グラフ |
| 9 | Fig. 2C | Sec. 3 | REBASE保存性 |
| 10 | Fig. 2D | Sec. 3 | O/E比較 |
| 11 | Fig. 2E | Sec. 3 | モチーフ時系列 |
| 12 | Fig. 2F | Sec. 3 | MEMEロゴ |
| 13 | Fig. 2G | Sec. 3 | Fold Change |
| 14 | Fig. 3A | Sec. 4 | メタジーンプロファイル |
| 15 | Fig. 3B | Sec. 4 | 距離帯別相関 |
| 16 | Fig. 3C | Sec. 4 | 6パネル散布図 |
| 17 | Fig. 3D | Sec. 4 | Concordanceパターン |
| 18 | Fig. 5A | Sec. 5 | 並べ替え検定 |
| 19 | Fig. 5B | Sec. 5 | 偏相関 |
| 20 | Fig. 6A | Sec. 6 | AAGCCCGプロモーター |
| 21 | Fig. 7A-I | Sec. 7 | DMG enrichment ×9 |
| 22 | Fig. 8A-D | Sec. 8 | TF BS methylation ×4 |
| 23 | Fig. 9A-E | Sec. 9 | GRN-TF cascade ×5 |
| 24 | Fig. 10A | Sec. 10 | 4mC/5mC competition |

### Supplementary Figure

| # | Figure ID | 内容 |
|---|-----------|------|
| S1-S8 | Fig. S-TD-* | 経時トラジェクトリー ×8 |
| S9-S13 | Fig. S-MT-* | マルチオミクストラック ×5 |
| S14-S15 | Fig. S-PS-* | プロモーター模式図 ×2 |
| S16-S20 | Fig. S-TSS-* | TSS関連解析 ×5 |
| S21 | Fig. S-ActBGC | Act BGCエピジェネティック |
| S22 | Fig. S-T3Coord | T3協調遺伝子比較 |
| S23 | Fig. S-Genus | 属間保存性 |
| S24-S25 | Fig. S-BGCReg, S-Cascade | BGC制御因子 ×2 |
| S26 | Fig. S-SC_RS17645 | SC_RS17645ホモロジー |
| P1-P4 | Paper Fig. 1-4 | 論文用総合Figure ×4 |

**総Figure数**: 本文24 + Supplementary 30+ + Paper 4 = **58+ Figure**

---

*最終更新: 2026-02-06*
*データバージョン: v2 (unweighted MIN_REPS=3, T3 re-sequenced, 2026-02-24)*
*全Figure: 2026-02-06再生成済み*
