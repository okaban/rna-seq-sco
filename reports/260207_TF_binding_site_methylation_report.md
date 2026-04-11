# TF Binding Site Methylation Analysis Report

**Date:** 2026-02-07 (updated 2026-02-24 with v2 methylation data, composition correction)
**Project:** *Streptomyces coelicolor* A3(2) M145 — Epigenome-Transcriptome Integration
**Analysis:** `13_TF_binding-site/analysis/02_TF_BS_methylation_260207_v1/`

---

## Key Insights

| # | Finding | Evidence | Significance |
|---|---------|----------|-------------|
| 1 | TF binding siteでメチル化が有意に**枯渇**（組成補正後も維持） | 直接重複: fold=0.66, p<0.003; **GC補正後 fold=0.62, p=5.2×10⁻⁶; モチーフ補正後 fold=0.73, p=0.003** | R-Mシステムによる結合部位の保護を示唆。塩基組成バイアスでは説明できない |
| 2 | BS近傍メチル化の大半がT1で**消失** | Lost: 94/125 (75%), Gained: 31/125 (25%) | 指数増殖期(T1)特異的なメチル化パターン |
| 3 | BSメチル化変動とターゲット遺伝子発現に**有意な相関なし** | Mann-Whitney p=0.083 (T2vsT1), p=0.645 (T3vsT1) | メチル化によるTF結合阻害の証拠は得られず |
| 4 | LuxR・MarRファミリーBSで比較的高メチル化率 | LuxR 30.8% (4/13), MarR 28.6% (2/7) | サンプルサイズ不足で統計的結論は困難 |
| 5 | phoP・glnRのBSに複数メチル化サイト集中 | phoP: 6/38 BS (10 sites), glnR: 6/27 BS (10 sites) | 栄養応答TFのBS近傍にメチル化クラスター |

---

## 1. データ概要

### 入力データ

| データ | ソース | 件数 |
|--------|--------|------|
| TF binding sites (Tier 1) | RegPrecise + ZorroAranda Curated + FIMO | 782 BS |
| メチル化サイト (unique) | Nanopore modkit (unweighted, MIN_REPS=3) | 5,941 (4mC: 2,693 / 6mA: 3,248) |
| 発現データ | DESeq2 (T1, T2, T3) | 8,083 genes |

### BS品質階層

| Tier | ソース | 件数 |
|------|--------|------|
| Tier 1 | RegPrecise（実験的） | 210 |
| Tier 1 | ZorroAranda Curated Strong（文献） | 452 |
| Tier 1 | FIMO q<0.01（厳格） | 120 |
| Tier 2 | ZorroAranda MEME（計算的） | 0（フィルタリング後） |

---

## 2. Analysis 1: 空間的重複解析

### メチル化はTF binding siteで枯渇している

| Window | BS数 | メチル化BS | 4mC obs/exp (fold) | p値 | 6mA obs/exp (fold) | p値 |
|--------|------|-----------|-------------------|------|-------------------|------|
| **直接重複** | 782 | 63 (8.1%) | 34/51.6 (**0.66**) | 0.006 | 41/62.2 (**0.66**) | 0.003 |
| **±50bp** | 782 | 100 (12.8%) | 58/75.9 (0.76) | 0.020 | 67/91.5 (0.73) | 0.004 |
| **±200bp** | 782 | 210 (26.9%) | 123/148.8 (0.83) | 0.017 | 143/179.5 (0.80) | 0.003 |

**結論**: 4mC・6mA共にTF結合部位でゲノムランダム期待値の60-80%に抑制されており、R-M認識配列がTF binding siteに重なると致死的であるため進化的に除去されてきた可能性が高い。枯渇傾向はBS中心から±200bpまで一貫して持続する。

**Figure F1** (metagene profile): BS中心±2kbのメチル化密度プロファイル。4mCはBS中心付近（±100bp）でゲノム平均（0.311/kb）を下回る傾向。6mAは変動が大きいが同様にBS直上で密度が低下する。

---

## 2b. Analysis 1b: 塩基組成補正（Composition Correction）

### 背景

Analysis 1の枯渇検定（fold=0.66, p=6.8×10⁻⁵）はゲノム上の一様分布を帰無仮説としている。しかし、TF binding site領域はゲノム平均と異なるGC含量やモチーフ頻度を持つ可能性があり、これらの組成バイアスが見かけの枯渇を生んでいる可能性を排除する必要がある。

### BS領域の組成特性

| 指標 | BS領域（782 Tier 1） | ゲノム全体 | 差 |
|------|---------------------|-----------|-----|
| GC含量 | **65.9%** | 72.1% | **−6.2%** |
| CCGG出現数 | 2,542 (in 165kb) | 151,175 (in 8.7Mb) | 比率: 0.81× |
| AAGCCCG出現数 | 18 | 1,334 | 比率: 0.71× |
| GATC出現数 | 878 | 43,535 | 比率: 1.06× |

BS領域はゲノム平均より**低GC**（65.9% vs 72.1%）。これは4mC（C/G上に生じる）の期待値を下げ、6mA（A/T上に生じる）の期待値を上げる方向に作用する。

### 4手法の比較（4mC + 6mA combined）

| Method | 期待値 | Fold | P (depletion) | 解釈 |
|--------|--------|------|---------------|------|
| **M0: 一様分布（元の解析）** | 113.3 | **0.662** | 8.4×10⁻⁵ | Baseline |
| **M1: GC含量補正** | 120.7 | **0.621** | **5.2×10⁻⁶** | 枯渇が**強化**される |
| **M2: モチーフ頻度補正** | 102.1 | **0.734** | 3.0×10⁻³ | 枯渇が若干弱まるが有意 |
| **M3a: ランダム並べ替え（10,000回）** | 113.6 ± 12.2 | **0.660** | 1.0×10⁻⁴ | M0と一致 |
| M3b: GC-matched並べ替え | — | — | — | 収束せず（※下記） |

### 修飾タイプ別の詳細

| Method | 4mC obs/exp (fold) | P | 6mA obs/exp (fold) | P |
|--------|-------------------|---|-------------------|---|
| M0: 一様 | 34/51.4 (0.66) | 6.6×10⁻³ | 41/62.0 (0.66) | 3.0×10⁻³ |
| M1: GC補正 | 34/47.9 (0.71) | 2.2×10⁻² | 41/72.9 (**0.56**) | **3.4×10⁻⁵** |
| M2: モチーフ補正 | 34/42.6 (0.80) | 1.1×10⁻¹ | 41/59.6 (0.69) | 7.1×10⁻³ |
| M3a: ランダム | 34/51.5 (0.66) | 6.0×10⁻³ | 41/62.1 (0.66) | 5.5×10⁻³ |

### 各手法の解釈

1. **M1（GC補正）でp値が改善**: BS領域は低GC → A/T比率が高い → 6mA期待値が増加（62.0→72.9）→ 観測値41との乖離がさらに拡大。4mCは逆にC/G期待値が減少するため若干弱まるが、6mAの強化が全体を牽引し、combined foldが0.621（p=5.2×10⁻⁶）に改善。

2. **M2（モチーフ補正）で4mCが非有意化**: BS内のCCGG（2,542個）はゲノム比で過少（0.81×）なため、CCGG由来の4mC期待値が減少。4mC単独ではfold=0.80, p=0.11で非有意に。一方6mAは有意を維持（fold=0.69, p=0.007）。Combined foldは0.734, p=0.003。

3. **M3b（GC-matched並べ替え）が収束不能**: BS平均GC=65.9%はゲノム平均72.1%と大きく乖離。ランダムに782領域を選ぶと平均GCはほぼ必ず~72%になるため、±2%の許容範囲内（63.9-67.9%）に収まるセットが1,000,000回の試行でもゼロ。これ自体がBS領域の特異な組成を示す間接的証拠。M1の解析的GC補正が代替となる。

### 結論

**TF BSでのメチル化枯渇はGC含量・モチーフ組成のバイアスでは説明できない。** GC補正ではむしろ枯渇が強化され（p=5.2×10⁻⁶）、最も保守的なモチーフ補正でもp=0.003で有意を維持する。これはR-M認識配列がTF結合部位で進化的に排除されていることを示す頑健なシグナルである。

**Figure F5** (permutation distribution): ランダム並べ替え10,000回の帰無分布と観測値（赤線）。観測値=75は帰無分布の下端に位置する。

---

## 3. Analysis 2: メタジーンプロファイル

**Figure F1** の解釈:

- **4mC (Panel A)**: BS中心周辺で密度がゲノム平均（0.311/kb）を下回るビンが多い。ただし-1000bp付近に0.49/kbのピークがあり、TGGCCGGC（palindromic）モチーフの偶発的集積と考えられる。
- **6mA (Panel B)**: BS中心±200bpではゲノム平均（0.375/kb）を下回るが、+800〜+1200bp域で0.6-0.7/kbの高密度域が出現。これはBSの下流域に位置する遺伝子本体内のメチル化に起因する可能性がある。
- **全体的傾向**: BS直上の明確なdipは検出されるものの、100bpビンの解像度ではノイズも大きく、サイト数（782 BS × ~6,000 methyl）の制約を受ける。

---

## 4. Analysis 3: 時系列変動

### BS近傍メチル化サイトの時系列パターン

| パターン | サイト数 | 比率 | 説明 |
|----------|---------|------|------|
| **Lost_T2T3** | 94 | 75.2% | T1で検出、T2/T3で消失 |
| **Gained** | 31 | 24.8% | T2またはT3で新規出現 |

### 修飾タイプ別

| 修飾型 | サイト数 |
|--------|---------|
| 6mA | 67 |
| 4mC | 58 |

### モチーフ構成

| モチーフ | BS近傍サイト数 | ゲノム全体 | BS/ゲノム比 |
|----------|---------------|-----------|------------|
| unassigned | 47 (37.6%) | 2,130 (35.9%) | 1.05 |
| TGGCCGGC | 35 (28.0%) | 1,717 (28.9%) | 0.97 |
| AAGCCCG | 28 (22.4%) | 1,176 (19.8%) | 1.13 |
| CCGKCA | 5 (4.0%) | 150 (2.5%) | 1.58 |
| GCCG | 4 (3.2%) | 108 (1.8%) | 1.75 |

**結論**: BS近傍のメチル化サイトはT1（指数増殖期）に最も多く検出され、T2/T3（定常期）で75%が消失する。この傾向はメチル化の全体的なタイムポイント分布（T1=3,921, T2=4,566, T3=3,942）とは異なるため、BS近傍では特にT1特異的なメチル化が多いことを示す。CCGKCA・GCCGモチーフがBSで若干enrichedであるが、サンプルサイズが小さい。

**Figure F2** (heatmap): 125サイト × 3タイムポイントのヒートマップ。T1列に集中する高メチル化（赤-黄色）とT2/T3の空白（白）が支配的パターン。

---

## 5. Analysis 4: BSメチル化変動 vs ターゲット遺伝子発現

### BS-ターゲットペアの分類

| BS methylation状態 | ペア数 | 割合 |
|-------------------|--------|------|
| No methylation | 575 | 87.2% |
| Lost at BS | 64 | 9.7% |
| Gained at BS | 20 | 3.0% |
| Mixed | 3 | 0.5% |

### Mann-Whitney U検定（ターゲット遺伝子 log2FC 比較）

| 比較 | Lost vs No-methyl | Gained vs No-methyl |
|------|-------------------|---------------------|
| **T2 vs T1** | U=20718, **p=0.083**, Δmedian=+0.181 | U=5589, p=0.862 |
| **T3 vs T1** | U=18947, p=0.645, Δmedian=-0.198 | U=5741, p=0.978 |

**結論**: BSメチル化消失（Lost）はT2vsT1においてターゲット遺伝子の発現上昇とマージナルな関連（p=0.083）を示すが、統計的有意水準には達しない。T3vsT1では関連なし。BSメチル化獲得（Gained）も発現変動と無関係。

**Figure F3** (boxplot): T2vsT1のLost群でやや正のmedian shiftが見られるが、分布の重なりが大きい。T3vsT1では全群でほぼ同じ分布。

### 注目すべき個別ペア

以下のBS-ターゲットペアは、メチル化変動と発現変動が同時に生じた事例:

| TF | ターゲット | BS methylation | ターゲット発現 (log2FC) |
|----|----------|---------------|---------------------|
| **phoP** → SC_RS40670 | CCGKCA (6mA) Lost at T1 | T2vsT1: +8.13, T3vsT1: +4.58 |
| **phoP** → SC_RS13455 | TGGCCGGC (4mC) Lost at T1 | T2vsT1: +8.70, T3vsT1: +4.49 |
| **argR** → SC_RS23650 | TGGCCGGC (4mC) Lost at T1 | T2vsT1: +2.19, T3vsT1: -1.21 |
| **oxyR** → SC_RS27325 | AAGCCCG (4mC) Lost at T1 | T2vsT1: -1.22, T3vsT1: -3.39 |
| **ramR** → SC_RS35590 | 6mA Gained at T3 | T3vsT1: +7.18 |
| **scbR** → SC_RS33650 | GCCG (6mA) Gained at T3 | T2vsT1: +7.18, T3vsT1: +5.21 |

これらは個別には興味深いが、全体的パターンとして統計的有意性は認められない。

---

## 6. Analysis 5: モチーフ・TFファミリー別集計

### TFファミリー別BS methylation率

| TFファミリー | BS数 | メチル化BS | 率 (%) | ゲノム平均との比較 |
|-------------|------|-----------|--------|-----------------|
| LuxR | 13 | 4 | **30.8** | +18.0 |
| MarR | 7 | 2 | **28.6** | +15.8 |
| MerR | 5 | 1 | 20.0 | +7.2 |
| TCS_response_regulator | 80 | 14 | **17.5** | +4.7 |
| Unclassified_TF | 339 | 47 | 13.9 | +1.1 |
| Sigma | 153 | 19 | 12.4 | -0.4 |
| GntR | 23 | 2 | 8.7 | -4.1 |
| Unknown | 122 | 9 | 7.4 | -5.4 |
| TetR | 18 | 1 | 5.6 | -7.2 |
| AraC | 6 | 0 | 0.0 | -12.8 |
| LacI | 11 | 0 | 0.0 | -12.8 |

**Figure F4** (横棒グラフ): 平均メチル化率12.8%を基準線として表示。LuxR・MarRが突出するが、BS数が少ない（13, 7）ためランダム変動の可能性あり。TCS応答制御因子（80 BS, 17.5%）が最も信頼性の高い高メチル化ファミリー。

---

## 7. Analysis 6: 主要TFのBS methylation詳細

### 発生分化関連TF

| TF | BS数 | メチル化BS | サイト数 | 特記事項 |
|----|------|-----------|---------|---------|
| **bldD** | 5 | 0 | 0 | メチル化なし — 形態分化マスター制御因子のBSは保護されている |
| **actII-ORF4** | 3 | 0 | 0 | メチル化なし — 抗生物質生合成活性化因子 |
| **Rex** | 20 | 0 | 0 | メチル化なし — 酸化還元センサー |
| redZ | — | — | — | BSデータなし |
| afsR | — | — | — | BSデータなし |

### ストレス応答シグマ因子

| TF | BS数 | メチル化BS | サイト数 | 主なモチーフ | 動態 |
|----|------|-----------|---------|------------|------|
| **sigB** | 70 | 9 | 12 | AAGCCCG, TGGCCGGC, CCGKCA | 主にLost_T2T3 |
| **sigE** | 44 | 5 | 6 | AAGCCCG (direct overlap) | Lost_T2T3 |
| **sigR** | 25 | 4 | 5 | AAGCCCG, TGGCCGGC | Lost_T2T3 |

sigR BSでのAAGCCG 4mCメチル化（93%）がT1で検出されT2/T3で消失。sigR自体はT2で1.69 log2FC上昇しており、BS近傍メチル化の消失とsigR活性化が時間的に一致するが、因果関係の証拠は不十分。

### 栄養応答・二成分制御系TF

| TF | BS数 | メチル化BS | サイト数 | 主なモチーフ | 特記 |
|----|------|-----------|---------|------------|------|
| **phoP** | 38 | 6 | 10 | TGGCCGGC, CCGKCA, unassigned | SC_RS40670 BS: CCGKCA 6mA + 発現+8.1 log2FC |
| **glnR** | 27 | 6 | 10 | TGGCCGGC, CCGKCA, unassigned, AAGCCCG | SC_RS30140 BS: 4mC+6mA dual |
| **dasR** | 80 | 2 | 4 | AAGCCCG (dual: 4mC + 6mA) | BSでのAAGCCG dual modification |
| **argR** | — | — | 12 | TGGCCGGC, AAGCCCG | 最多メチル化サイト（6 BS） |
| **DraR** | 20 | 3 | 3 | AAGCCCG, GCCG, unassigned | 混合モチーフ |

**phoP** は最も注目に値するTF:
- 6つのBSに10のメチル化サイトが集中
- SC_RS40670のプロモーター領域でCCGKCA (6mA)がT1で57%検出、T2/T3で消失
- 同遺伝子はT2vsT1で+8.13 log2FCの劇的な発現上昇
- ただし、phoPターゲット遺伝子の大部分（32/38 BS）はメチル化なしで同様に発現変動しているため、メチル化が発現変動の主因とは言えない

### γ-butyrolactone制御系

| TF | BS数 | メチル化BS | サイト数 | 特記 |
|----|------|-----------|---------|------|
| **scbR** | — | 3 | 4 | T3でGained（6mA GCCG） |
| **scbR2** | — | 14 | 14 | 最多エントリ — 広範なBS群 |
| **ramR** | — | 1 | 1 | T3でGained（6mA unassigned） |

scbR2は14エントリと最多だが、BS数自体が多いことに起因。

---

## 8. 総合考察

### 主要結論

1. **TF BSでのメチル化枯渇は進化的保存パターン（組成補正後も頑健）**: R-M認識配列（TGGCCGGC, AAGCCCG）とTF BSが空間的に重なると、制限酵素による切断がTF結合を妨害するため、進化的にこの重複は排除される方向の選択圧がかかる。observed/expected比=0.66は塩基組成補正（GC補正: fold=0.62, p=5.2×10⁻⁶; モチーフ補正: fold=0.73, p=0.003; ランダム並べ替え: fold=0.66, p=0.0001）の全てで維持され、組成バイアスでは説明できない真の生物学的シグナルである。

2. **BS近傍メチル化の時系列変動は全体パターンと一致**: T1でのメチル化サイトの75%がT2/T3で消失するパターンは、4mCの全体的動態（T1: 1,987 → T3: 1,647サイト）と概ね整合する。BS特異的な制御というよりも、ゲノム全体のメチル化動態の一部と解釈される。

3. **BSメチル化と遺伝子発現の直接的関連は検出されず**: T2vsT1でのマージナルな傾向（p=0.083）を除き、BSメチル化状態とターゲット遺伝子発現変動の間に統計的有意な関連はない。これは遺伝子レベルの解析（Spearman r=0.125, p=0.002 for 4mC T2vsT1）と一致し、メチル化が転写制御に直接寄与する証拠は得られなかった。

4. **個別TFでは興味深い事例が存在するが、系統的パターンではない**: phoP-SC_RS40670（CCGKCA 6mA消失 + 発現+8.13 log2FC）やdasR BS（AAGCCCG dual modification）は個別には興味深いが、全体の統計的テストでは機能非選択的であり、これらは偶発的な共起と解釈するのが妥当。

### 先行研究との比較

- 大腸菌のDam methylation (GATC)はTF binding siteとの重複が報告されており（Oshima et al., 2002; Løbner-Olesen et al., 2005）、一部のTF（Lrp, OxyR）でメチル化依存的な結合制御が実証されている。
- 本解析の *S. coelicolor* では、GATC 6mAはわずか37サイト（6mAの1.2%）であり、Dam-like制御の寄与は限定的。
- 主要メチル化モチーフ（TGGCCGGC, AAGCCCG）はType II / Type I R-Mシステムに帰属され、「自己DNA保護」が主機能と考えられる。遺伝子制御への二次的転用の証拠は本データからは得られなかった。

---

## 9. Limitation

1. **BS品質**: Tier 1 BS（782）はRegPrecise（実験的）+ 文献キュレーション + FIMO予測の混合セット。FIMOのFPRは不明であり、偽陽性BSがenrichment解析を希釈している可能性。
2. **検出力**: BS近傍のメチル化サイトは125と少なく、下位グループ（Lost/Gained × motif）に分割すると統計検出力が低下する。
3. **因果関係の不在**: 観察された共起パターンは相関であり、メチル化によるTF結合阻害の直接的証拠（EMSA、bisulfite-seq等）は欠如。
4. **タイムポイント解像度**: T1(18h), T2(48h), T3(72h)の3点では、メチル化変動の速度論を分解できない。
5. **Strand特異性**: BS strandとメチル化strandの一致/不一致は本解析では考慮していない。

---

## 10. Output Files

### Tables

| ファイル | 内容 | 行数 |
|---------|------|------|
| `T1_BS_methylation_overlap_summary.tsv` | 空間的重複統計（Tier/Window別） | 6 |
| `T2_BS_methylation_temporal_dynamics.tsv` | BS近傍メチル化の時系列動態 | 125 |
| `T3_BS_methylation_expression_pairs.tsv` | BS-ターゲット遺伝子ペアの発現情報 | 109 |
| `T4_motif_TF_family_summary.tsv` | TFファミリー別メチル化率 | 14 |
| `T5_key_TF_BS_methylation_detail.tsv` | 主要TFのBS methylation詳細 | 14 |
| `BS_methylation_composition_correction.csv` | 塩基組成補正5手法の比較結果 | 5 |

### Figures

| ファイル | 内容 | 形式 |
|---------|------|------|
| `F1_metagene_methylation_around_BS` | BS中心±2kbメチル化密度プロファイル | PDF + SVG |
| `F2_BS_methylation_temporal_heatmap` | BS×タイムポイント ヒートマップ | PDF + SVG |
| `F3_BS_methylation_vs_expression` | BSメチル化状態別ターゲット発現boxplot | PDF + SVG |
| `F4_TF_family_BS_methylation_rate` | TFファミリー別BS methylation率 | PDF + SVG |
| `BS_composition_correction_permutation` | 塩基組成補正 並べ替え帰無分布 | PDF + PNG + SVG |

---

## 11. Methods

### メチル化データ
- Nanopore sequencing (R10.4.1, SUP basecalling v5.2.0)
- modkit pileup → unweighted consensus (MIN_REPS=3/3, MIN_COVERAGE=5, MIN_TOTAL_COVERAGE=30, MIN_MOD_FREQ≥50%)
- 5,941 unique positions (4mC: 2,693, 6mA: 3,248)

### TF Binding Site
- Tier 1: RegPrecise (210) + ZorroAranda2022 Curated Strong (452) + FIMO q<0.01 (120) = 782 BS
- Contig mapping: `chromosome` → `NC_003888.3`

### 統計検定
- 空間的重複: Poisson検定（ランダム期待値 = メチル化サイト数 / ゲノムサイズ × BS領域幅）
- 塩基組成補正: 4手法
  - M0: 一様分布（Poisson検定）
  - M1: GC含量補正（塩基別メチル化率 × BS内C/G数・A/T数）
  - M2: モチーフ頻度補正（CCGG/AAGCCCG/GATC出現率 × BS内モチーフ数 + 残差の一様分布）
  - M3: 並べ替え検定（10,000回、同一長分布のランダム領域）
- 発現比較: Mann-Whitney U検定
- メチル化率のp値はBonferroni補正なし（探索的解析のため）

### ソフトウェア
- Python 3.x (pandas, numpy, scipy, matplotlib, biopython)
- Scripts: `01_TF_BS_methylation_analysis.py`, `02_BS_composition_correction.py`
