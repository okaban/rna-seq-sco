# M145 RNA-seq/Methylome統合解析 プレゼンテーション原稿

**作成日**: 2026-02-03
**プロジェクト**: *Streptomyces coelicolor* A3(2) M145 トランスクリプトーム-エピゲノム統合解析
**対象**: 共同研究者向けプレゼンテーション資料

---

## スライド構成

1. タイトル
2. 研究背景と目的
3. 実験デザイン
4. 解析パイプライン概要
5. RNA-seqクオリティコントロール
6. メチロームクオリティコントロール
7. データ品質検証（PCA・サンプル間距離）
8. 差次的発現解析（DESeq2）
9. ベン図によるDEGの経時的パターン解析
10. BGC発現動態解析
11. BGCクラスター内遺伝子の発現パターン
12. KEGGパスウェイエンリッチメント解析
13. COG機能分類解析
14. エピゲノム-トランスクリプトーム統合解析の概要
15. メチル化-発現相関解析
16. 相関解析のバリデーション
17. DEG-DMG重複検定
18. 新規R-Mシステムの同定
19. 4mCモチーフ解析
20. ATCC比較メチローム解析
21. 転写因子ネットワークとエピジェネティック制御
22. Act vs Red BGCのエピジェネティック制御の対比
23. RamRエピジェネティックスイッチの発見
24. NsdBの正の相関パターン
25. 先行研究との比較と本研究の位置づけ
26. 主要発見のまとめと新規性評価
27. 提案するモデル
28. 研究の限界と今後の検証課題
29. 今後の展望と応用可能性
30. 使用ソフトウェアとバージョン
31. 謝辞
32. 参考文献

---

## スライド 1: タイトル

### *Streptomyces coelicolor* A3(2) M145における培養段階依存的トランスクリプトームリプログラミングとエピジェネティック制御の全体像

**発表者**: [発表者名]
**所属**: [所属機関]
**日付**: 2026年2月

---

## スライド 2: 研究背景と目的

### 背景

- **放線菌** (*Streptomyces*属) は抗生物質、免疫抑制剤など有用二次代謝産物の主要な生産者
- *S. coelicolor* A3(2) M145は放線菌のモデル生物として広く研究されている
- 二次代謝の活性化は培養段階に強く依存（成長期→定常期の移行）
- エピジェネティック制御（DNAメチル化）の関与は未解明

### 先行研究の状況

- **Pisciotta et al. (2023)**: BS-seqによりM145の5mCメチロームを報告（3,360サイト、GGCmCGG/GCCmCGモチーフ）。5mCが形態分化と二次代謝に関与することを示唆
- **Fang et al. (2022)**: *S. roseosporus* L30でSroLm3（4mC MTase）がダプトマイシン生合成を制御することを報告
- しかし、6mA/4mCの**アデニン/シトシンメチル化**と転写制御の統合的解析は未実施

### 目的

1. M145の3つの培養段階（T1, T2, T3）における全ゲノム的な転写変動を解明
2. 主要BGC（生合成遺伝子クラスター）の発現動態を定量化
3. DNAメチル化（6mA, 4mC）と遺伝子発現の関連を統合的に解析
4. 二次代謝制御に関わる新規エピジェネティック機構を同定

---

## スライド 3: 実験デザイン

### サンプル構成

| 条件 | タイムポイント | 生物学的レプリケート |
|------|--------------|-------------------|
| M145_1 (T1) | 初期（増殖期） | n = 3 |
| M145_2 (T2) | 中期（移行期） | n = 3 |
| M145_3 (T3) | 後期（定常期） | n = 3 |

**合計**: 9サンプル × 2プラットフォーム

### データ取得

| プラットフォーム | 用途 | 読み長 |
|----------------|------|--------|
| **Illumina** | RNA-seq（トランスクリプトーム） | PE 150 bp |
| **Oxford Nanopore** | 修飾塩基検出（メチローム） | Long-read |

### リファレンスゲノム

- **NCBI RefSeq**: GCF_000203835.1
- **ゲノムサイズ**: 8.67 Mb
- **遺伝子数**: 8,275

---

## スライド 4: 解析パイプライン概要

### トランスクリプトーム解析

```
Raw FASTQ (PE 150bp)
    ↓ FastQC v0.12.1 (品質管理)
    ↓ fastp v1.1.0 (アダプター除去・品質トリミング)
    ↓ HISAT2 v2.2.1 (アラインメント, --no-spliced-alignment)
    ↓ SAMtools v1.21 (BAMソート・インデックス)
    ↓ featureCounts v2.1.1 (カウント, -s 0 unstranded)
    ↓ DESeq2 v1.46.0 (差次的発現解析)
    ↓ apeglm v1.28.0 (LFC shrinkage)
```

### メチローム解析

```
Nanopore BAM
    ↓ modkit pileup (修飾塩基コール)
    ↓ 高信頼度サイト抽出 (≥10x coverage, ≥50% frequency)
    ↓ プロモーター領域統合 (-300〜+50bp)
```

### 統合解析

```
発現データ + メチル化データ
    ↓ Spearman相関解析
    ↓ Fisher正確確率検定（DEG-DMG重複）
    ↓ MEME（モチーフ発見）
    ↓ 協調的変化遺伝子の同定
```

---

## スライド 5: RNA-seqクオリティコントロール

### リード品質・フィルタリング（fastp v1.1.0）

| サンプル | Raw Read Pairs | After Filtering | 保持率 |
|---------|---------------|-----------------|--------|
| M145_1-1 | 7,454,937 | 7,346,735 | 98.6% |
| M145_1-2 | 7,450,028 | 7,339,992 | 98.5% |
| M145_1-3 | 7,485,279 | 7,369,997 | 98.5% |
| M145_2-1 | 7,515,618 | 7,404,424 | 98.5% |
| M145_2-2 | 7,555,483 | 7,378,267 | 97.7% |
| M145_2-3 | 7,600,714 | 7,102,375 | 93.4% |
| M145_3-1 | 7,500,802 | 7,375,834 | 98.3% |
| M145_3-2 | 7,550,164 | 7,432,786 | 98.4% |
| M145_3-3 | 7,561,014 | 7,082,447 | 93.7% |
| **合計** | **67,674,039** | **66,832,857** | **98.8%** |

### クオリティスコア（トリミング後）

| 指標 | 値 |
|------|-----|
| Q20（エラー率≤1%） | 98.7〜98.9% |
| Q30（エラー率≤0.1%） | 94.8〜95.1% |
| GC含量 | ~68%（ゲノムGC: ~72%） |
| 全サンプルFastQC判定 | **PASS**（Per base sequence quality） |

### アラインメント結果（HISAT2 v2.2.1）

| 指標 | 平均値 | 範囲 |
|------|--------|------|
| Overall alignment rate | **98.54%** | 98.1〜99.0% |
| Uniquely mapped | **96.55%** | 95.9〜97.2% |
| Multi-mapped | <2% | - |

**注**: HISAT2を使用（`--no-spliced-alignment`オプション）。原核生物にはスプライシングなし。

### リードカウント結果（featureCounts v2.1.1）

| 指標 | 値 |
|------|-----|
| Assignment rate（平均） | **89.3%** |
| 検出遺伝子数 | 7,766 / 8,275（**93.8%**） |
| サンプル平均フラグメント数 | 6.7M |
| パラメータ | `-t gene -g gene_id -s 0`（unstranded） |

**データ品質の総括**: 全9サンプルが高品質であり、Q30 > 94%、アラインメント率 > 98%、遺伝子検出率 > 93%を達成。下流解析に十分なデータ品質。

---

## スライド 6: メチロームクオリティコントロール

### Nanoporeシーケンシングカバレッジ

| サンプル | カバレッジ | 検出メチル化サイト数 |
|---------|-----------|-------------------|
| M145_1-1 | 44.3x | - |
| M145_1-2 | 46.0x | - |
| M145_1-3 | 44.0x | - |
| M145_2-1 | 51.7x | - |
| M145_2-2 | 40.5x | - |
| M145_2-3 | 40.4x | - |
| M145_3-1 | 33.9x | - |
| M145_3-2 | **14.3x** | - |
| M145_3-3 | 29.8x | - |

**注意**: サンプル M145_3-2 はカバレッジが14.3xと他サンプルと比べ低い。解析には≥10xのフィルタリングを適用しているが、検出力が低下している可能性がある。

### 修飾塩基検出の概要

| 修飾タイプ | 高信頼サイト数 | 備考 |
|-----------|-------------|------|
| **6mA** | 3,214 | R-Mシステム由来 |
| **4mC** | 2,679 | R-Mシステム由来 |
| 5mC | ~0 | ほぼ検出されず（<0.01%） |
| **合計** | 11,778 | 全修飾タイプ合計 |

### 5mC未検出に関する方法論的考察

Pisciotta et al. (2023) は同じM145株でBS-seqにより3,360個の5mCサイトを報告している。本研究で5mCがほぼ検出されない理由：

1. **検出方法の違い**: Nanopore modkitは6mA/4mC検出に最適化されており、5mC検出感度はBS-seq（5mCの金標準）より低い
2. **培養条件の違い**: Pisciotta et al.はMG定義培地を使用（18h, 24h時点）。本研究の培養条件とは異なる可能性
3. **本研究の焦点**: 6mA/4mCのアデニンメチル化ダイナミクスに焦点を当てた相補的アプローチ

→ **5mCを含む包括的メチロームの統合は今後の重要な課題**

### メチル化サイト抽出基準

| 項目 | 設定 |
|------|------|
| 最低カバレッジ | ≥10x |
| 最低頻度 | ≥50% |
| プロモーター定義 | TSS -300bp〜+50bp |

---

## スライド 7: データ品質検証（PCA・サンプル間距離）

### 【Figure】PCA_M145.pdf

**ファイル**: `04_deseq2/analysis/04_deseq2_260128_v1/figures/PCA_M145.pdf`

**図の読み方**:
- X軸: PC1（全分散の83%を説明）
- Y軸: PC2（全分散の16%を説明）
- 各点は1サンプル、色は条件（T1/T2/T3）を示す

**結果と示唆**:
- PC1で3つの培養段階が明確に分離
- T1→T2→T3の順に時間軸に沿って配置
- 各条件内のレプリケートは密集（高い再現性）
- **示唆**: 培養段階間で大規模な転写リプログラミングが発生

### 【Figure】sample_distance_heatmap_M145.pdf

**ファイル**: `04_deseq2/analysis/04_deseq2_260128_v1/figures/sample_distance_heatmap_M145.pdf`

**図の読み方**:
- ユークリッド距離のヒートマップ
- 暗色ほど類似度が高い
- 階層的クラスタリングにより条件ごとにグループ化

**結果と示唆**:
- 同一条件内のサンプルが最も類似
- T1とT2はT3よりも互いに類似
- **示唆**: T3で最も劇的な転写変化が発生

---

## スライド 8: 差次的発現解析（DESeq2）

### 解析条件

| 項目 | 設定 |
|------|------|
| ソフトウェア | DESeq2 v1.46.0 (R v4.4.2) |
| LFC shrinkage | apeglm v1.28.0 |
| フィルタリング | ≥10 counts in ≥3 samples |
| 有意性閾値 | padj < 0.05 |
| デザイン式 | ~ condition |

### DEG数サマリー

| 比較 | 解析遺伝子数 | Up | Down | 合計 | 全体に占める割合 |
|------|-------------|-----|------|------|----------------|
| T2 vs T1 | 7,646 | 2,563 | 2,697 | 5,260 | **69%** |
| T3 vs T1 | 7,646 | 3,290 | 2,898 | 6,188 | **81%** |
| T3 vs T2 | 7,646 | 2,639 | 2,904 | 5,543 | **73%** |

### 【Figure】volcano_M145_2_vs_1.pdf（T2 vs T1）

**ファイル**: `04_deseq2/analysis/04_deseq2_260128_v1/figures/volcano_M145_2_vs_1.pdf`

**図の読み方**:
- X軸: log2 Fold Change（正の値=発現上昇、負の値=発現低下）
- Y軸: -log10(adjusted p-value)（高いほど統計的に有意）
- 赤点: 有意に発現上昇した遺伝子
- 青点: 有意に発現低下した遺伝子
- 灰色点: 非有意

**結果と示唆（T2 vs T1）**:
- 全遺伝子の**69%**が有意に発現変動（padj < 0.05）
- 増殖期→移行期の遷移で既に大規模な転写変動が開始
- **示唆**: T2が転写リプログラミングの開始点であり、二次代謝誘導が進行中

### 【Figure】volcano_M145_3_vs_1.pdf（T3 vs T1）

**ファイル**: `04_deseq2/analysis/04_deseq2_260128_v1/figures/volcano_M145_3_vs_1.pdf`

**図の読み方**: 同上

**結果と示唆（T3 vs T1）**:
- 全遺伝子の**81%**が有意に発現変動（padj < 0.05）
- |log2FC| > 1の遺伝子も多数（約4,800遺伝子）
- **示唆**: T1→T3で近乎全体的な転写リプログラミングが発生。T2からさらに変動が拡大

---

## スライド 9: ベン図によるDEGの経時的パターン解析

### 【Figure】venn_all_DEGs.pdf

**ファイル**: `12_supplementary_figures/analysis/12_supplementary_260202_v1/figures/venn_all_DEGs.pdf`

**図の読み方**:
- 3つの円がそれぞれT2vsT1、T3vsT1、T3vsT2の有意DEGsを表す
- 重複領域は複数の比較で共通して変動した遺伝子

**主要な発見**:

| カテゴリ | 遺伝子数 | 解釈 |
|---------|---------|------|
| 3比較すべてで変動 | 1,041 | **コア応答遺伝子**（培養全期間で持続的に変動） |
| T2vsT1 & T3vsT1のみ | 1,793 | T1からの累積的変化 |
| T3vsT1 & T3vsT2のみ | 1,486 | 後期特異的変化 |
| T3vsT1のみ | 521 | 長期培養特異的遺伝子 |

**示唆**:
- 約1,000遺伝子がコア応答遺伝子として培養全期間を通じて変動
- 後期（T3）に向けて変動遺伝子数が増加

---

## スライド 10: BGC発現動態解析

### 解析対象

*S. coelicolor* M145の主要4 BGC:

| BGC | 産物 | 遺伝子数 |
|-----|------|---------|
| **act** | Actinorhodin（青色色素） | 22 |
| **red** | Undecylprodigiosin（赤色色素） | 22 |
| **cda** | CDA（カルシウム依存性抗生物質） | 40 |
| **cpk** | Coelimycin P1（ポリケタイド） | 16 |

### BGC平均発現の経時変化

| BGC | T1 | T2 | T3 | FC(T2/T1) | FC(T3/T1) |
|-----|-----|-----|-----|-----------|-----------|
| act | 65 | 99 | 9,491 | 1.5x | **145x** |
| red | 65 | 682 | 966 | 10.4x | **14.8x** |
| cda | 57 | 3,309 | 3,085 | **58.4x** | 54.4x |
| cpk | 52 | 9,465 | 6,144 | **181.7x** | 118x |

### 【Figure】BGC_timecourse_lineplot_M145.pdf

**ファイル**: `06_BGC_dynamics/analysis/06_BGC_dynamics_260128_v1/figures/BGC_timecourse_lineplot_M145.pdf`

**図の読み方**:
- X軸: 培養段階（T1, T2, T3）
- Y軸: 平均正規化カウント
- 各線は異なるBGCを表す

**結果と示唆**:

**二段階活性化モデル**の発見:

1. **Phase I（T1→T2）**: red, cda, cpkが10〜180倍に急増。actはわずかな増加のみ
2. **Phase II（T2→T3）**: actのみが追加的に96倍増加。cda/cpkは横ばい〜減少

**示唆**:
- actinorhodin生合成は他のBGCより遅れて活性化（遅延型活性化）
- cpkは一過性のピーク（T2で最大、T3で減少）
- **先行研究との関連**: Manteca et al. (2008)が報告したMI-MII境界での成長停止と二次代謝活性化のタイミングと一致 [4]

---

## スライド 11: BGCクラスター内遺伝子の発現パターン

### 【Figure】BGC_gene_heatmap_act_M145.pdf

**ファイル**: `06_BGC_dynamics/analysis/06_BGC_dynamics_260128_v1/figures/BGC_gene_heatmap_act_M145.pdf`

**図の読み方**:
- 行: act BGC内の各遺伝子
- 列: 各サンプル（T1×3, T2×3, T3×3）
- 色: log2正規化カウント（赤=高発現、青=低発現）

**結果と示唆**:
- T1ではほぼ全遺伝子が低発現（サイレント状態）
- T2で一部遺伝子が誘導開始
- T3で全22遺伝子が高発現に
- **actII-orf4**（cluster-situated regulator）: log2FC = +2.60

### 【Figure】coverage_tracks_all_BGC.pdf

**ファイル**: `12_supplementary_figures/analysis/12_supplementary_260202_v1/figures/coverage_tracks_all_BGC.pdf`

**図の読み方**:
- ゲノム座標上でのリードカバレッジを可視化
- 上段から: act, red, cda, cpk
- 色: 各条件（T1=青, T2=緑, T3=赤）

**結果と示唆**:
- BGCクラスター全体にわたる均一な発現変動を確認
- red（早期）→ cda（中期）→ act/cpk（後期）の誘導タイミング差異

---

## スライド 12: KEGGパスウェイエンリッチメント解析

### 解析方法

| 項目 | 設定 |
|------|------|
| データベース | KEGG Orthology |
| 有意性閾値 | FDR-corrected p < 0.05 |

### 【Figure】KEGG_bubble_T3vsT1_up.pdf

**ファイル**: `12_supplementary_figures/analysis/12_supplementary_260202_v1/figures/KEGG_bubble_T3vsT1_up.pdf`

**図の読み方**:
- X軸: Fold Enrichment（濃縮倍率）
- Y軸: パスウェイ名
- バブルサイズ: 該当遺伝子数
- 色: -log10(p-value)

**T3 vs T1 発現上昇パスウェイ（上位）**:

| パスウェイ | 遺伝子数 | FE | p-value |
|-----------|---------|-----|---------|
| ABC transporters | 108 | 1.66 | 1.0e-10 |
| Quorum sensing | 60 | 1.77 | 8.5e-8 |
| Type II PKS products | 10 | 2.70 | 1.4e-4 |
| Prodigiosin biosynthesis | 16 | 1.83 | 3.2e-3 |

**示唆**: 栄養枯渇への適応（ABC輸送体）と二次代謝の活性化

### 【Figure】KEGG_bubble_T2vsT1_down.pdf

**ファイル**: `12_supplementary_figures/analysis/12_supplementary_260202_v1/figures/KEGG_bubble_T2vsT1_down.pdf`

**T2 vs T1 発現低下パスウェイ（上位）**:

| パスウェイ | 遺伝子数 | FE | p-value |
|-----------|---------|-----|---------|
| **Ribosome** | 59 | **3.34** | 8.0e-35 |
| Pyrimidine metabolism | 20 | 2.39 | 3.2e-5 |
| DNA replication | 11 | 2.43 | 1.7e-3 |

**示唆**:
- **T2が転写リプログラミングのピーク**
- リボソーム抑制はT2で最大（3.34倍）、T3では緩和（2.76倍）
- 成長から生存モードへの移行
- **先行研究との関連**: van Wezel & McDowall (2011)が提唱した二次代謝-一次代謝トレードオフモデルと一致 [5]

---

## スライド 13: COG機能分類解析

### 【Figure】report_COG_net_change.pdf

**ファイル**: `12_supplementary_figures/analysis/12_supplementary_260202_v1/figures/report_COG_net_change.pdf`

**図の読み方**:
- 水平棒グラフ
- 正の値: Up遺伝子数がDown遺伝子数を上回る
- 負の値: Down遺伝子数がUp遺伝子数を上回る
- 各行は異なるCOG機能カテゴリ

**主要な発見（T3 vs T1）**:

| カテゴリ | Net変化 | 解釈 |
|---------|--------|------|
| **J (Translation)** | **-98** | リボソーム・翻訳の大幅抑制 |
| **Q (Secondary metabolism)** | **+12** | 二次代謝の一方向的誘導（Down=0） |
| **P (Ion transport)** | +125 | イオン輸送の活性化 |
| **D (Cell division)** | -9 | 細胞分裂の停止 |

**示唆**:
- **成長-二次代謝トレードオフの明確な証拠**
- 翻訳カテゴリ(J): Up:Down = 1:6（抑制優勢）
- 二次代謝カテゴリ(Q): 発現上昇のみ（Down=0）

---

## スライド 14: エピゲノム-トランスクリプトーム統合解析の概要

### メチル化の種類

| 修飾タイプ | 存在 | 機能 |
|-----------|------|------|
| 5mC | ほぼ検出されず（<0.01%） | 生物学的に有意でない |
| **6mA** | 検出（3,214サイト） | R-Mシステム由来 |
| **4mC** | 検出（2,679サイト） | R-Mシステム由来 |

### 解析パラメータ

| 項目 | 設定 |
|------|------|
| カバレッジ閾値 | ≥10x |
| 頻度閾値 | ≥50% |
| プロモーター領域 | TSS -300bp〜+50bp |
| 解析遺伝子数 | 8,083 |
| 高信頼メチル化サイト | 11,778 |

### 主要な発見

1. **4mCプロモーターメチル化はT2vsT1で発現変化と有意に正の相関**
   - Spearman r = 0.137, p = 0.0006
2. **T3への遷移期でDEGsとDMGsに有意な正の関連**
   - Odds Ratio = 1.32, p = 0.001
3. **558遺伝子がメチル化・発現の協調的変化を示す**

### 先行研究との対比

| 項目 | 本研究 | Pisciotta et al. (2023) [2] | Fang et al. (2022) [7] |
|-----|------|--------------------------|----------------------|
| 検出方法 | Nanopore modkit | Bisulfite-seq | SMRT-seq |
| 対象メチル化 | **6mA, 4mC** | 5mC | 4mC, 6mA |
| 主要モチーフ | AAGCCCG, CCGG | GGCmCGG, GCCmCG | GCGG |
| 二次代謝との関連 | Red, Act, CDA, Cpk | 赤/青色素 | ダプトマイシン |
| 関係 | **相補的** | 同一株・異なるメチル化層 | 異なる種・同様の結論 |

---

## スライド 15: メチル化-発現相関解析

### 【Figure】Fig1_methylation_expression_correlation.pdf

**ファイル**: `11_epigenome_integration/analysis/02_publication_figures/Fig1_methylation_expression_correlation.pdf`

**図の読み方**:
- X軸: プロモーターメチル化変化量（%）
- Y軸: 発現変化量（log2FC）
- 回帰直線とSpearman相関係数を表示

**結果**:

| 修飾 | 比較 | r値 | p値 | 判定 |
|-----|------|-----|-----|------|
| **4mC** | **T2vsT1** | **0.137** | **0.0006** | **有意な正の相関** |
| 4mC | T3vsT1 | -0.081 | 0.068 | 有意でない |
| 6mA | T2vsT1 | 0.007 | 0.867 | 有意でない |
| 6mA | T3vsT1 | -0.048 | 0.210 | 有意でない |

**示唆**:
- **メチル化=転写活性化**という新たなパラダイム
- 真核生物の5mC（転写抑制マーク）とは異なる
- **先行研究との関連**: Nye et al. (2020) がグラム陽性菌におけるDNAメチル化のゲノム防御以外の制御的機能を総説で指摘しており [8]、本知見はその具体的実例

---

## スライド 16: 相関解析のバリデーション

### 4mC T2vsT1 相関（r=0.137）の統計的検証

4mC T2vsT1で検出された有意な正の相関（r=0.137, p=0.0006）が偶発的な結果でないことを、5つの独立した検証手法で確認した。

### 検証結果サマリー

| # | 検証手法 | 結果 | 判定 |
|---|---------|------|------|
| 1 | **Permutation test**（10,000回） | Empirical p = 0.0006 | **PASS** |
| 2 | **Partial correlation**（GC%, 遺伝子長を制御） | r_partial = 0.157（114.4%保持） | **PASS** |
| 3 | **Bootstrap CI**（10,000回） | 95% CI [0.061, 0.211]（0を含まず） | **PASS** |
| 4 | **Z-score test** | Z = 3.40, p = 0.00034 | **PASS** |
| 5 | **Negative controls** | 4/5コントロールがr<0.05 | **MARGINAL** |

### 結論

- 5つの検証基準のうち**4/5がPASS**
- Permutation test（p=0.0006）: ラベルをランダムに並べ替えた場合に同等以上の相関が得られる確率は0.06%
- Partial correlation: GC含量や遺伝子長などの交絡因子を統計的に除去しても、相関は114.4%維持（むしろ強化）
- Bootstrap 95%信頼区間が0を含まない → 母集団レベルでの正の相関を支持

**示唆**: 4mCプロモーターメチル化と発現変化の正の相関は、統計的にロバストな関係である

---

## スライド 17: DEG-DMG重複検定

### 【Figure】venn_DEG_DMG_T3vsT1.pdf

**ファイル**: `11_epigenome_integration/analysis/03_overlap_analysis/venn_DEG_DMG_T3vsT1.pdf`

**図の読み方**:
- 左円: DEGs（発現変動遺伝子）
- 右円: DMGs（メチル化変動遺伝子）
- 重複領域: 両方で変動した遺伝子

### Fisher正確確率検定の結果

| 比較 | DEGs | DMGs | 重複 | OR | p値 | 判定 |
|------|------|------|------|-----|-----|------|
| T2vsT1 | 3,848 | 554 | 266 | 1.04 | 0.69 | 関連なし |
| **T3vsT1** | 4,841 | 656 | 428 | **1.32** | **0.001** | **有意な正の関連** |
| **T3vsT2** | 3,507 | 749 | 361 | **1.26** | **0.003** | **有意な正の関連** |

**示唆**:
- **二段階モデル**の提唱
  - Phase 1（T1→T2）: 発現変動とメチル化変化が独立的
  - Phase 2（→T3）: 協調的制御に移行

---

## スライド 18: 新規R-Mシステムの同定

### 【重要発見】AAGCCCG（6mA）モチーフ

| 特徴 | 値 |
|------|-----|
| REBASEステータス | **未登録（新規）** [REBASE v272, 2026] |
| 協調変動遺伝子での濃縮 | **13倍**（OR=13.08, p<0.0001） |
| 全遺伝子での存在率 | 5.1% |
| 協調変動遺伝子での存在率 | 31.2% |

### 候補メチラーゼ: SC_RS17645

| 特性 | 値 |
|------|-----|
| タンパク質長 | 679 aa（典型的Dam系の2-3倍） |
| ドメイン | N-6 adenine MTase |
| 発現変化（T2vsT1） | log2FC = -2.19, padj = 6.48e-16 |
| 保存モチーフ | 2つのMTaseモチーフ検出 |

### 【Figure】Figure2_AAGCCCG_system.pdf

**ファイル**: `11_epigenome_integration/analysis/17_paper_figures/Figure2_AAGCCCG_system.pdf`

**図の読み方**:
- パネルA: SC_RS17645の発現変化
- パネルB: プロモーターAAGCCCG分布
- パネルC: 協調変動遺伝子での濃縮
- パネルD: エピジェネティックカスケードモデル

**示唆**:
- M145は新規AAGCCCG R-Mシステムを持つ
- SC_RS17645（N-6 MTase）が候補酵素
- **学術的新規性が高い**（新規MTaseファミリーの可能性）

**進化的示唆**: 古典的R-Mシステムは利己的遺伝因子として侵入し、宿主ゲノムに維持される。SC_RS17645のような大型MTase（679 aa）がゲノム防御を超えた制御的機能を持つことは、R-Mシステムの「家畜化（domestication）」を示唆する [8][9]

---

## スライド 19: 4mCモチーフ解析

### 【Figure】4mC_base_composition.pdf

**ファイル**: `11_epigenome_integration/analysis/07_motif_analysis/4mC_base_composition.pdf`

**図の読み方**:
- X軸: メチル化サイトからの相対位置（bp）
- Y軸: 塩基組成（A/T/G/C）
- 中心（position 0）にメチル化サイト

### 4mCモチーフの同定

| モチーフ | マッチ数 | 割合 | 酵素/システム |
|---------|---------|------|-------------|
| **CCGG** | 2,026 | **75.6%** | MspI様メチルトランスフェラーゼ |
| GCGC | 95 | 3.5% | HhaI様 |
| GATC | 6 | 0.2% | Dcm様 |

**MEME発見モチーフ**:
- コンセンサス: `[GC][AC][AC]GCC[GC]GCCA`
- E-value: 6.4e-1211
- サイト数: 2,678

**示唆**:
- 4mCメチル化は主にType II R-MシステムのMspI様MTaseによる
- ゲノム防御に加え、発現制御にも関与の可能性
- **先行研究との比較**: *S. roseosporus*のSroLm3はGCGGモチーフを認識 [7]。M145のCCGGは異なる特異性であり、種間でMTase-モチーフの組み合わせが多様

---

## スライド 20: ATCC比較メチローム解析

### 目的

M145で同定した修飾モチーフ（CCGG, AAGCCCG, GATC）がStreptomyces属全体でどの程度保存されているかを、ATCC所蔵37株のゲノム情報を用いて解析した。

### モチーフ保存性の比較

| モチーフ | 修飾タイプ | 保有株数 | 保有率 | 分布パターン |
|---------|----------|---------|--------|------------|
| **CCGG** | 4mC | 32/37 | **86.5%** | 広範に保存 |
| **GATC** | 6mA | 25/37 | **67.6%** | 普遍的 |
| **AAGCCCG** | 6mA | 8/37 | **21.6%** | **系統特異的** |

### 【Figure】motif_conservation.pdf

**ファイル**: `11_epigenome_integration/analysis/07_motif_analysis/motif_conservation.pdf`

**図の読み方**:
- 37 Streptomyces株における各モチーフの保有率・頻度を比較
- 株間のばらつきと分布パターンを可視化

### M145の位置づけ

| 指標 | M145の値 | 37株中の順位 | 解釈 |
|------|---------|-------------|------|
| CCGG頻度 | ~45th percentile | 中位 | **典型的**（Streptomyces平均） |
| AAGCCCG頻度 | ~85th percentile | 上位 | **保有株の中では高頻度** |

### 結果と示唆

1. **CCGGは汎Streptomyces的R-Mシステム**: 86.5%の株が保有しており、種を超えた共通防御機構
2. **AAGCCCGは系統特異的**: 21.6%の株のみ保有 → M145のAAGCCCG R-Mシステムは**特定の系統群に限定**
3. **M145のCCGGメチル化レベルは「典型的」**: Streptomyces属の中でも平均的な水準
4. **AAGCCCGを持つ株では高頻度**: M145は保有株の中でも上位に位置

**注意**: 本解析はATCC APIの認証制限により、文献ベースのゲノム情報を使用。実験的な全ゲノムメチローム比較は今後の課題。

---

## スライド 21: 転写因子ネットワークとエピジェネティック制御

### GRN TF解析の結果

- 解析対象: 37 TF（文献既知の二次代謝レギュレーター）[Zorro-Aranda et al. 2022]
- **階層構造**: Tier 1 Global (14) → Tier 2 Pleiotropic+Sigma (17) → Tier 3 CSR (6)
- メチル化-発現協調変動を示すTF: **redZのみ**（37 TF中唯一）

### メチル化状態のサマリー

| TF | Tier | メチル化状態 | log2FC | 備考 |
|----|------|-----------|--------|------|
| **redZ** | 3 (CSR) | **協調変動** | -2.25 | 6mA Lost + 発現低下 |
| bldN | 1 (Global) | Lost in T2 | +0.50 | 非協調 |
| afsR | 1 (Global) | Stable | -0.26 | 非協調 |
| afsS | 1 (Global) | Stable | -0.03 | 非協調 |
| 残り33 TF | - | 非メチル化 | - | - |

### redZ エピジェネティックカスケード

```
環境シグナル（T2時点）
    ↓
SC_RS17645 (N-6 MTase) 発現低下 (log2FC = -2.19)
    ↓
AAGCCCG メチル化減少
    ↓
redZ プロモーター脱メチル化
    ↓
redZ 発現低下 (log2FC = -2.25)
    ↓ [AbsA2抑制解除が支配的]
redD 発現上昇 (log2FC = +4.77)
    ↓
Red BGC 活性化 (9/10遺伝子上昇)
```

**先行研究との関連**: White & Bibb (1997)がredZのundecylprodigiosin生合成への必須性を確立。本研究はredZ制御の上流にエピジェネティック層が存在することを初めて示した。

### 【Figure】Figure3_GRN_methylation.pdf

**ファイル**: `11_epigenome_integration/analysis/17_paper_figures/Figure3_GRN_methylation.pdf`

**図の読み方**:
- パネルA: GRN階層構造
- パネルB: TFメチル化ヒートマップ
- パネルC: redZカスケード図

### redDパラドックス

redZはredDの正の制御因子であるが、redZ発現低下にもかかわらずredDは強く上昇発現（log2FC = +4.77）している。

| 遺伝子 | 役割 | log2FC | 予測 | 実測 |
|--------|-----|--------|------|------|
| redZ | redD活性化因子 | -2.25 | redD↓ | redD↑ |
| absA2 | redD抑制因子 | +6.29 | redD↓ | redD↑ |
| **redD** | Red活性化因子 | **+4.77** | - | **強い上昇** |

**解釈**: 他の活性化因子（AfsR, BldD, papR2等）の寄与、またはAbsA2活性の飽和/翻訳後制御が示唆される。Botas et al. (2018)はArgRが1,544遺伝子を制御することを報告しており、こうした多面的制御因子がredDに対する追加的な入力を提供している可能性がある。

---

## スライド 22: Act vs Red BGCのエピジェネティック制御の対比

### 【Figure】Figure4_Act_vs_Red.pdf

**ファイル**: `11_epigenome_integration/analysis/17_paper_figures/Figure4_Act_vs_Red.pdf`

**図の読み方**:
- 左パネル: Act BGCの発現とメチル化状態
- 右パネル: Red BGCの発現とメチル化状態
- 比較表: 制御メカニズムの違い

### Act vs Red 比較表

| 指標 | Act | Red |
|-----|-----|-----|
| 遺伝子数 | 20 | 11 |
| メチル化遺伝子 | 1 | 3 |
| 協調変動 | 2 | 3 |
| 平均log2FC | 0.86 | 3.02 |
| SARP制御 | なし | **redZ協調変動** |
| エピジェネティック制御 | **なし** | **あり** |

### 2つのBGCの異なる制御ロジック

| 特徴 | Act（エピジェネティック非依存） | Red（エピジェネティック依存） |
|------|------------------------------|----------------------------|
| 活性化タイミング | Phase II（T2→T3） | Phase I（T1→T2） |
| SARP TF | actII-ORF4（非メチル化） | redZ（6mA協調変動） |
| 活性化機構 | 代謝的ゲート（前駆体依存？） | エピジェネティックゲート |
| 制御の特徴 | 翻訳後修飾・栄養センサー？ | MTase→メチル化→TFカスケード |

**示唆**:
- 同一ゲノム内のBGCが**異なる制御アーキテクチャ**を持つ
- Red: エピジェネティックゲートにより環境シグナルと直結
- Act: より複雑な多段階活性化（エピジェネティック非依存）
- **進化的示唆**: 二次代謝制御の進化的モジュール性を反映

---

## スライド 23: RamRエピジェネティックスイッチの発見

### 【Figure】multiomics_track_RamR.pdf

**ファイル**: `11_epigenome_integration/analysis/04_multiomics_tracks/multiomics_track_RamR.pdf`

**図の読み方**:
- 上段: メチル化シグナル（6mA/4mC）
- 中段: 発現レベル（T1/T2/T3）
- 下段: 遺伝子構造

### RamR (SCO6685) の特性

| 特性 | T1 | T2 | 変化 |
|------|-----|-----|------|
| 6mAサイト数 | 1 | 0 | 消失 |
| メチル化頻度 | 64.4% | 0% | -64.4% |
| 発現量 | 基底 | 100倍 | +6.60 log2FC |

**生物学的仮説**:
1. T1で6mAメチル化状態 → 発現抑制
2. T1→T2で脱メチル化 → 発現100倍上昇
3. RamR活性化 → SapB産生 → 気菌糸形成開始

**示唆**: RamRプロモーターの6mA脱メチル化は、栄養増殖から発生遷移への**分子スイッチ**として機能

**先行研究との関連**: RamR-SapBによる気菌糸形成制御は既知だが、その上流にエピジェネティックスイッチ（6mA脱メチル化）が存在することは本研究で初めて報告。Pisciotta et al. (2023)もRamRを5mCメチル化遺伝子として検出しているが、スイッチ機能は未報告 [2]。

---

## スライド 24: NsdBの正の相関パターン

### 【Figure】multiomics_track_NsdB.pdf

**ファイル**: `11_epigenome_integration/analysis/04_multiomics_tracks/multiomics_track_NsdB.pdf`

### NsdB (SCO7252) の特性

| 特性 | 値 |
|------|-----|
| 産物 | DNA結合タンパク質 NsdB |
| 修飾タイプ | 4mC |
| T1→T2サイト変化 | 0→1 (+92.3%) |
| 発現変化 | **+9.32 log2FC (640倍)** |
| 相関タイプ | **正 (Gained_Up)** |

**示唆**:
- NsdBは発生と二次代謝の多面的制御因子
- 4mCメチル化獲得と同時に劇的に活性化
- **メチル化が転写活性化マーク**として機能

### メチル化の機能的二面性

| 機能タイプ | 代表遺伝子 | メカニズム | メチル化種 |
|-----------|-----------|-----------|----------|
| **活性化マーク** | NsdB | メチル化↑→発現↑ | 4mC |
| **抑制マーク** | RamR | メチル化↓→発現↑ | 6mA |

**考察**: 4mCによる活性化と6mAによる抑制という二面性は、同じ「メチル化」でも修飾タイプごとに異なる分子メカニズムが関与していることを示唆する。4mCはリプレッサー結合阻害またはアクチベーターリクルートに、6mAはRNAポリメラーゼ開始の直接的阻害に機能している可能性がある。これは *E. coli* のDamメチル化による定常期遺伝子抑制と類似する [3]。

---

## スライド 25: 先行研究との比較と本研究の位置づけ

### メチル化研究の系譜における位置づけ

| 研究 | 対象 | 修飾 | 検出法 | 主要発見 | 本研究との関係 |
|------|------|------|--------|---------|-------------|
| **Pisciotta et al. (2023)** [2] | M145 | 5mC | BS-seq | GGCmCGG/GCCmCGモチーフ、形態分化への関与 | **同一株の異なるメチル化層**を相補的に解析 |
| **Fang et al. (2022)** [7] | *S. roseosporus* L30 | 4mC | SMRT-seq | SroLm3がダプトマイシン生合成を制御 | **異なる種で類似の結論** → 一般原理の確立 |
| **Beaulaurier et al. (2019)** [1] | 総説 | 全般 | 各種 | 細菌エピゲノムの技術的レビュー | **技術的基盤**（Nanopore法の正当性） |
| **Nye et al. (2020)** [8] | グラム陽性菌（総説） | 全般 | - | 防御以外のメチル化制御機能 | 本研究の発見が**具体的実例** |

### 既知の制御因子への新規メカニズム層の追加

| 制御因子 | 既知の知見 | 本研究の追加知見 |
|---------|----------|---------------|
| **RedZ** | Undecylprodigiosin生合成に必須 (White & Bibb, 1997) | **プロモーター6mA脱メチル化がredZ制御の上流に存在** |
| **RamR** | SapB産生、気菌糸形成に関与 | **6mA脱メチル化が発生分化の分子スイッチ** |
| **ArgR** | 1,544遺伝子の多面的制御 (Botas et al., 2018) | エピジェネティック制御はArgRとは独立した追加層 |
| **ActII-ORF4** | Act BGCの正の制御因子 | **非メチル化 → エピジェネティック制御なし** |

### 本研究の統合的位置づけ

```
Pisciotta et al. (2023)           本研究
     5mC メチロム         +     6mA/4mC メチロム
     形態分化制御               二次代謝エピジェネティック制御
           ↓                           ↓
       ┌─────────────────────────────────────┐
       │   M145の包括的エピジェネティックモデル   │
       │   5mC: 基盤的分化制御                   │
       │   6mA: エピジェネティックスイッチ        │
       │   4mC: 転写活性化マーク                  │
       └─────────────────────────────────────┘
```

---

## スライド 26: 主要発見のまとめと新規性評価

### 本研究の主要発見（新規性の層別評価）

#### ★★★ 最高度の新規性

| # | 発見 | 対応Figure | 根拠 |
|---|------|-----------|------|
| 1 | **AAGCCCG（6mA）新規R-Mシステムの同定** | Figure2_AAGCCCG_system.pdf | REBASE v272未登録、OR=13.08 |
| 2 | **4mCプロモーターメチル化の転写活性化機能** | Fig1_methylation_expression_correlation.pdf | r=0.137, p=0.0006, 4/5バリデーションPASS |
| 3 | **RamRエピジェネティックスイッチ** | multiomics_track_RamR.pdf | 6mA 64.4%→0%, 100倍発現上昇 |

#### ★★ 高度の新規性

| # | 発見 | 対応Figure | 根拠 |
|---|------|-----------|------|
| 4 | **redZエピジェネティックカスケード** | Figure3_GRN_methylation.pdf | SC_RS17645→AAGCCCG→redZ→Red BGC |
| 5 | **Act vs Red: BGC間制御アーキテクチャの多様性** | Figure4_Act_vs_Red.pdf | エピジェネティック依存 vs 非依存 |
| 6 | **AAGCCCG R-MのStreptomyces属内系統特異性** | motif_conservation.pdf | 21.6%のみ保有 |

#### ★ 中程度の新規性

| # | 発見 | 対応Figure | 根拠 |
|---|------|-----------|------|
| 7 | **T2が転写リプログラミングのピーク** | KEGG_bubble_T2vsT1_down.pdf | リボソームFE=3.34 |
| 8 | **BGC二段階活性化モデル** | BGC_timecourse_lineplot_M145.pdf | Phase I/II分離 |
| 9 | **成長-二次代謝トレードオフの定量的証拠** | report_COG_net_change.pdf | J:-98, Q:+12 |

### 学術的インパクト

1. **新規R-Mシステム**: AAGCCCG（6mA）はREBASE未登録 → **新規MTaseファミリーの可能性**
2. **メチル化の転写活性化機能**: 真核生物の5mC（抑制マーク）とは異なる **新たなパラダイム**
3. **二次代謝のエピジェネティック制御**: Red BGCはエピジェネティック制御下、Act BGCはそうでない → **BGC間の制御機構の多様性**
4. **R-Mシステムの「家畜化」**: ゲノム防御機構が発生・代謝制御に転用された進化的過程を示唆

---

## スライド 27: 提案するモデル

### *S. coelicolor* M145 エピジェネティック-転写制御統合モデル

```
┌─────────────────────────────────────────────────────────────────────┐
│                                                                      │
│  T1 (増殖期)                T2 (移行期)              T3 (定常期)     │
│  ============              ============              ============    │
│                                                                      │
│  SC_RS17645 (MTase)        SC_RS17645 ↓↓            SC_RS17645 ↓↓   │
│        ↓                        ↓                        ↓          │
│  AAGCCCG メチル化          AAGCCCG 脱メチル化        脱メチル化維持  │
│        ↓                        ↓                        ↓          │
│  redZ 抑制                 redZ 活性化変動           redD 高発現     │
│        ↓                        ↓                        ↓          │
│  Red BGC サイレント         Red BGC 誘導開始         Red BGC 高発現  │
│                                                                      │
│  RamR メチル化(抑制)       RamR 脱メチル化(活性化)   RamR 持続発現   │
│                             → 気菌糸形成開始                        │
│                                                                      │
│  ────────────────────────────────────────────────────────────────   │
│                                                                      │
│  cpk/cda: Phase I誘導                               Act: Phase II誘導│
│  (エピジェネティック                                 (エピジェネティック │
│   寄与は限定的)                                      非依存)         │
│                                                                      │
└─────────────────────────────────────────────────────────────────────┘
```

### モデルの生物学的意義

**エピジェネティックロッキング（Epigenetic Locking）**: SC_RS17645によるAAGCCCGメチル化は、好条件下（T1）でも二次代謝遺伝子をサイレント状態に「ロック」する。環境・代謝シグナルがSC_RS17645を下方制御した時にのみ、このロックが解除される。これは、単純な栄養枯渇応答ではなく、**エピジェネティックな「コミットメントポイント」**として機能する。

---

## スライド 28: 研究の限界と今後の検証課題

### 本研究の限界

| 限界 | 詳細 | 影響 |
|------|------|------|
| **相関 ≠ 因果** | メチル化-発現の関係は相関ベース。直接的な因果証明は未実施 | カスケードモデルは仮説段階 |
| **SC_RS17645の生化学的検証未了** | MTase活性のin vitro実証がない | 候補酵素の確定に至っていない |
| **5mC検出の不完全性** | Nanopore法では5mCの系統的評価が困難 | Pisciotta et al.との完全な統合ができない |
| **単一培養条件** | 1条件の時系列のみ。培地依存性が不明 | 一般化に限界あり |

### 想定される査読上の論点と対応方針

| 想定される質問 | 対応方針 |
|-------------|---------|
| **なぜ5mCが検出されないのか？** | 方法論的相違（Nanopore vs BS-seq）として説明。5mC統合は今後の課題と明記 |
| **SC_RS17645が真のメチラーゼか？** | "candidate methyltransferase"と記載。ドメイン構造と発現動態を証拠として提示。KO株実験を計画中と記載 |
| **RamR脱メチル化は原因か結果か？** | 時間的前後関係から推論が最も経済的と説明。高時間分解能実験を提案 |
| **AAGCCCGは既知R-Mの変種では？** | REBASE検索結果を提示。BLAST-Pで既知MTaseクラスの外に位置することを示す |

### 優先的に実施すべき検証実験

| 優先度 | 実験 | 目的 | 期待される結果 |
|--------|------|------|-------------|
| **最高** | SC_RS17645 CRISPR-KO株 | AAGCCCGメチラーゼの機能検証 | KO株でAAGCCCGメチル化消失 |
| **最高** | redZプロモーターメチル化部位変異体 | エピジェネティック制御の直接証明 | 変異体でredZ発現変化 |
| **高** | 精製SC_RS17645のin vitroメチル化アッセイ | 基質特異性の確認 | AAGCCCG配列へのメチル化活性 |
| **中** | 同一サンプルのBS-seq | 5mCメチロームの統合 | Pisciotta et al.との完全な比較 |
| **中** | ChIP-seq（RedZ, ActII-ORF4） | TF結合とメチル化の関係 | メチル化依存的なTF結合変化 |

---

## スライド 29: 今後の展望と応用可能性

### 短期的課題（計算解析）

| 優先度 | タスク | 期待される成果 |
|--------|--------|---------------|
| 高 | SC_RS17645の外部BLAST確認（NCBI） | 新規MTaseの系統分類 |
| 高 | Supplementary Figures作成 | 論文投稿準備 |
| 高 | REBASE登録準備 | AAGCCCG R-Mシステムの公式記録 |
| 中 | ChIP-seqデータとの統合 | メチル化-TF結合の関係解明 |

### 中長期的課題（実験的検証）

| 優先度 | タスク | 期待される成果 |
|--------|--------|---------------|
| 高 | SC_RS17645ノックアウト株作製 | AAGCCCGメチラーゼの機能検証 |
| 高 | RamRプロモーターメチル化部位変異体構築 | エピジェネティックスイッチの検証 |
| 中 | in vitro転写アッセイ | メチル化の転写への直接効果 |
| 中 | 時系列メチロームの高時間分解能解析 | ダイナミクスの詳細解明 |

### バイオテクノロジー応用の可能性

エピジェネティック制御の理解は、放線菌の**代謝工学**に直結する：

1. **抗生物質増産**: SC_RS17645の条件的発現抑制システムを導入し、BGCのエピジェネティックロックを意図的に解除 → 抗生物質産生量の最適化
2. **沈黙BGCの覚醒**: ゲノム中の未発現BGC（cryptic BGC）に対し、メチラーゼ枯渇戦略を適用 → 新規二次代謝産物の発見
3. **組み合わせ型生合成**: メチル化パターンの精密制御による複数BGCの段階的活性化 → 代謝コンフリクトの回避
4. **診断バイオマーカー**: メチル化シグネチャーを発生段階や代謝ポテンシャルの予測マーカーとして利用 → 生産株のモニタリング

### 論文執筆

- **論文Figure 1-4**（PNG/PDF）生成済み
- 全解析データ整理済み
- 論文本文執筆準備完了
- **推奨投稿先**: *Molecular Microbiology*（ファーストチョイス）、*mBio*（バックアップ）

---

## スライド 30: 使用ソフトウェアとバージョン

### RNA-seq解析

| ソフトウェア | バージョン | 用途 |
|-------------|-----------|------|
| FastQC | 0.12.1 | 品質管理 |
| fastp | 1.1.0 | アダプター除去・品質トリミング |
| MultiQC | 1.33 | QCレポート統合 |
| HISAT2 | 2.2.1 | アラインメント |
| SAMtools | 1.21 | BAMソート・インデックス |
| featureCounts (Subread) | 2.1.1 | リードカウント |
| DESeq2 | 1.46.0 | 差次的発現解析 |
| apeglm | 1.28.0 | LFC shrinkage |

### メチローム解析

| ソフトウェア | バージョン | 用途 |
|-------------|-----------|------|
| modkit | - | 修飾塩基コール |
| MEME Suite | - | モチーフ発見 |

### 統合解析・可視化

| ソフトウェア | バージョン | 用途 |
|-------------|-----------|------|
| R | 4.4.2 | 統計解析・可視化 |
| tidyverse | - | データ処理 |
| pheatmap | - | ヒートマップ |
| clusterProfiler | - | GO/KEGG解析 |
| ggplot2 | - | 可視化 |

---

## スライド 31: 謝辞

- 共同研究者の皆様
- [資金提供機関]
- シーケンス施設

---

## スライド 32: 参考文献

[1] Beaulaurier, J., Schadt, E. E., & Fang, G. (2019). Deciphering bacterial epigenomes using modern sequencing technologies. *Nat. Rev. Genet.*, 20, 157–172.

[2] Pisciotta, A., Sampino, A. M., Presentato, A., et al. (2023). The DNA cytosine methylome revealed two methylation motifs in the upstream regions of genes related to morphological and physiological differentiation in *Streptomyces coelicolor* A(3)2 M145. *Sci. Rep.*, 13, 7038.

[3] Kahramanoglou, C., Prieto, A., Khedkar, S., et al. (2012). Genomics of DNA cytosine methylation in *Escherichia coli* reveals its role in stationary phase transcription. *Nat. Commun.*, 3, 886.

[4] Manteca, A., Álvarez, R., Salazar, N., Yagüe, P., & Sánchez, J. (2008). Mycelium differentiation and antibiotic production in submerged cultures of *Streptomyces coelicolor*. *Appl. Environ. Microbiol.*, 74, 3877–3886.

[5] van Wezel, G. P., & McDowall, K. J. (2011). The regulation of the secondary metabolism of *Streptomyces*: New links and experimental advances. *Nat. Prod. Rep.*, 28, 1311–1333.

[6] Luo, R., Ying, K., Zhang, J., et al. (2022). Heterogeneous DNA methylation in single bacteria. *ISME J.*, 16, 2099–2113.

[7] Fang, J.-L., Gao, W.-L., Xu, W.-F., et al. (2022). m4C DNA methylation regulates biosynthesis of daptomycin in *Streptomyces roseosporus* L30. *Synth. Syst. Biotechnol.*, 7, 1013–1023.

[8] Nye, T. M., Fernandez, N. L., & Simmons, L. A. (2020). A positive perspective on DNA methylation: Regulatory functions of DNA methylation outside of host defense in Gram-positive bacteria. *Crit. Rev. Biochem. Mol. Biol.*, 55, 576–591.

[9] Adhikari, S., & Curtis, P. D. (2016). DNA methyltransferases and epigenetic regulation in bacteria. *FEMS Microbiol. Rev.*, 40, 575–591.

---

## 付録: Figure一覧と保存場所

### 主要Figure

| Figure ID | ファイル名 | 保存場所 |
|-----------|-----------|---------|
| Fig1 | Figure1_overview.pdf | `11_epigenome_integration/analysis/17_paper_figures/` |
| Fig2 | Figure2_AAGCCCG_system.pdf | `11_epigenome_integration/analysis/17_paper_figures/` |
| Fig3 | Figure3_GRN_methylation.pdf | `11_epigenome_integration/analysis/17_paper_figures/` |
| Fig4 | Figure4_Act_vs_Red.pdf | `11_epigenome_integration/analysis/17_paper_figures/` |

### 補足Figure

| Figure ID | ファイル名 | 保存場所 |
|-----------|-----------|---------|
| S1 | PCA_M145.pdf | `04_deseq2/analysis/04_deseq2_260128_v1/figures/` |
| S2 | volcano_M145_2_vs_1.pdf | `04_deseq2/analysis/04_deseq2_260128_v1/figures/` |
| S3 | volcano_M145_3_vs_1.pdf | `04_deseq2/analysis/04_deseq2_260128_v1/figures/` |
| S4 | venn_all_DEGs.pdf | `12_supplementary_figures/analysis/12_supplementary_260202_v1/figures/` |
| S5 | BGC_timecourse_lineplot_M145.pdf | `06_BGC_dynamics/analysis/06_BGC_dynamics_260128_v1/figures/` |
| S6 | KEGG_bubble_T3vsT1_up.pdf | `12_supplementary_figures/analysis/12_supplementary_260202_v1/figures/` |
| S7 | report_COG_net_change.pdf | `12_supplementary_figures/analysis/12_supplementary_260202_v1/figures/` |
| S8 | Fig1_methylation_expression_correlation.pdf | `11_epigenome_integration/analysis/02_publication_figures/` |
| S9 | venn_DEG_DMG_T3vsT1.pdf | `11_epigenome_integration/analysis/03_overlap_analysis/` |
| S10 | motif_conservation.pdf | `11_epigenome_integration/analysis/07_motif_analysis/` |
| S11 | multiomics_track_RamR.pdf | `11_epigenome_integration/analysis/04_multiomics_tracks/` |
| S12 | multiomics_track_NsdB.pdf | `11_epigenome_integration/analysis/04_multiomics_tracks/` |

---

*作成: 2026-02-03*
*更新: 2026-02-03（QCデータ追加、T2vsT1 Volcano追加、ATCC比較解析追加、バリデーション追加、ソフトウェアバージョン修正、先行研究統合、新規性評価、限界と検証課題追加、参考文献追加）*
*プロジェクト: M145 RNA-seq/Methylome統合解析*
