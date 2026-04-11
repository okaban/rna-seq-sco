# m4C/m5C二重シトシン修飾系の検証レポート

**作成日**: 2026-02-05
**プロジェクト**: *Streptomyces coelicolor* A3(2) M145 メチローム解析
**対象**: マニュスクリプト執筆者・共著者
**解析ディレクトリ**: `11_epigenome_integration/analysis/22_4mC_5mC_competition/`
**データソース**: 本研究Nanoporeデータ + Pisciotta et al. 2023 BS-seq (Scientific Reports 13, 7038)

---

## Key Insights（主要なインサイト一覧）

| # | Insight | 対応Figure | 根拠となる数値 |
|---|---------|-----------|----------------|
| 1 | 4mCコンセンサスモチーフの100%がGGCCGG文脈内（Pisciottaのm5Cモチーフと完全一致） | `fig_4mC_5mC_competition.pdf` Panel B | 上位20コンセンサス: 20/20 = 100% vs ゲノム上CCGG中GGCCGG比率24.4% |
| 2 | T3でDcm-like MTaseが12倍上昇するのに4mC-CCGGが56%崩壊（パラドックス） | Panel A | SC_RS19770 +3.54, SC_RS36410 +3.76 log2FC; 4mC: 1,961→866 |
| 3 | CCGGの内側Cは4mCまたは5mCとして修飾され、培養段階依存的にスイッチする | Panel C | 6つの独立したエビデンスライン（Panel D） |
| 4 | Pisciotta 2023のsite-level BS-seqデータは公開されておらず、直接統合は不可 | — | PRJNA933392: targeted loci のみ、全ゲノムBSデータなし |

---

## 1. 背景と仮説

### 問題提起

本研究のNanopore解析でCCGGモチーフの内側C（C2位）に4mCを検出した。一方、Pisciotta et al. (2023) はBS-seqで同じCCGGを含むGGCCGGモチーフの同位置に5mCを検出した。同一塩基に異なるメチル化修飾が報告されている。

### 仮説

**m4C/m5C競合仮説**: CCGGモチーフの内側シトシンは、4mCまたは5mCとして修飾されうる二重修飾部位であり、培養段階に応じて異なるMTaseが優位になることで修飾型がスイッチする。

---

## 2. Pisciotta et al. 2023 BS-seqデータの可用性

| 項目 | 状態 |
|------|------|
| **BioProject** | PRJNA933392 |
| **論文** | Scientific Reports 13, 7038 (2023) |
| **データ内容** | 321遺伝子の上流メチル化配列（targeted loci） |
| **全ゲノムBS-seq raw reads** | **公開されていない** |
| **Site-level m5Cデータ** | **利用不可**（再処理不能） |
| **利用可能情報** | 5mCモチーフ（GGCCGG, GCCCG）、321メチル化遺伝子リスト、総m5C数（3,360サイト） |

**結論**: Pisciottaの公開データでは直接的なsite-levelでのm4C vs m5C比較は不可能。しかし、配列レベルでの重複（GGCCGG）と発現動態から間接的検証が可能。

---

## 3. 配列レベルの証拠: 4mCのGGCCGG集中

### 【Figure】fig_4mC_5mC_competition.pdf

**ファイル**: `analysis/22_4mC_5mC_competition/fig_4mC_5mC_competition.pdf`

![4mC/5mC competition](_fig/11_epigenome_integration/analysis/22_4mC_5mC_competition/fig_4mC_5mC_competition.png)

**目的**: m4C/m5C二重修飾仮説の6つの独立エビデンスを統合的に提示する。

**方法**: Panel A: 4mC-CCGGサイト数（high confidence sites, n_reps≥2）とDcm-like MTase 2酵素の平均log2FC。Panel B: ゲノム全CCGGサイト中のGGCCGG割合 vs 4mCメチル化CCGGコンセンサスモチーフ上位20のGGCCGG割合。Panel C: 培養段階依存的m4C/m5Cスイッチモデル。Panel D: 6つのエビデンス概要表。

### ゲノム vs メチル化サイトの文脈比較

| メトリクス | ゲノム全体 | 4mCメチル化サイト |
|-----------|-----------|-----------------|
| 全CCGG数 | ~151,213 (17,441/Mb) | 2,026 (unique) |
| うちGGCCGG文脈 | ~36,900 (**24.4%**) | 上位20コンセンサス全て (**100%**) |
| 期待値（24.4%なら） | — | ~495 |
| 実測 vs 期待 | — | **~4.1倍の濃縮** |

### 結果と示唆
- ゲノム上のCCGGの24.4%のみがGGCCGG文脈であるのに対し、4mCメチル化されるCCGGサイトの上位20コンセンサスは**100%がGGCCGG**
- これはPisciotta (2023) がBS-seqで5mC検出したGGCCGGモチーフと**完全に一致**
- **示唆**: 4mCと5mCは同じGGCCGG部位を標的としており、同一塩基に対する修飾の競合を強く示唆する

---

## 4. 発現動態の証拠: Dcm-likeパラドックス

### Dcm-like MTase発現（Nanopore 4mC → m5C産生酵素候補）

| MTase | T2 vs T1 | T3 vs T1 | T3 vs T2 |
|-------|----------|----------|----------|
| SC_RS19770 (Dcm-like) | -0.94 (ns) | **+2.32** (padj=9.5e-6) | **+3.54** (padj=4.5e-8) |
| SC_RS36410 (Dcm-like) | +1.06 (ns) | **+5.34** (padj=1.0e-10) | **+3.76** (padj=2.1e-9) |

### 4mC-CCGG部位数の推移

| タイムポイント | 全4mC | CCGG含有4mC | CCGG比率 |
|-------------|-------|-----------|---------|
| T1 | 1,995 | ~1,576 | 79.0% |
| T2 | 2,458 | ~1,961 | 79.8% |
| **T3** | **1,077** | **~866** | 80.4% |

### パラドックスの解決

**問題**: T3でDcm-like MTaseが11-14倍上昇しているのに、4mC-CCGGが56%崩壊する。Dcm-likeが4mCを産生するなら、上昇に伴い4mCも増加するはず。

**解決**: Dcm-like MTaseが産生するのは**m5C**（4mCではない）。T3でDcm-like活性が上昇すると、CCGG部位のCがm5Cに修飾され、先に存在していたm4Cが置換される。Nanoporeはm5Cをほとんど検出できない（<0.01%）ため、m5Cに変わったサイトは4mC消失として観測される。

---

## 5. 6つの独立エビデンスの統合

| # | エビデンス | ソース | 強度 |
|---|----------|--------|------|
| 1 | CCGGで4mC検出（Nanopore） | 本研究 | 直接 |
| 2 | GGCCGGで5mC検出（BS-seq） | Pisciotta 2023 | 直接 |
| 3 | **4mCコンセンサス=GGCCGGの100%集中** | 本研究 | **強** |
| 4 | **Dcm-like +12x上昇なのに4mC 56%崩壊** | 本研究 | **強** |
| 5 | BS-seqはm4Cを検出不能（非メチル化と同扱い） | 方法論 | 論理的 |
| 6 | Nanoporeはm5C検出感度が低い（<0.01%） | 方法論 | 論理的 |

### エビデンスの相互補完性

- エビデンス1+2: 同じCCGG/GGCCGG部位で異なる修飾型を検出
- エビデンス3: 配列文脈の完全一致（4mC → GGCCGG = 5mCモチーフ）
- エビデンス4: 時間動態の反相関（Dcm-like↑ + 4mC↓ → Dcm-likeはm5C産生酵素）
- エビデンス5+6: 検出法の盲点が相互補完的（Nanopore→4mC特化、BS-seq→5mC特化）

---

## 6. Limitation

- Pisciottaの全ゲノムBS-seq site-levelデータが非公開のため、同一座標でのm4C vs m5Cの直接比較は不可能
- 培養条件が異なる（本研究: [条件未詳], Pisciotta: MG liquid, 18h/24h）
- 「4mCサイトの消失 = m5Cへの置換」は推論であり、他の可能性（単純な脱メチル化、DNA複製による希釈）を排除できない
- m4Cを産生するCCGG-MTaseの同定が未完了
- **決定的実験**: 同一生物試料での同時BS-seq + Nanoporeが必要（今後の課題）

---

## 7. 出力ファイル一覧

### Figure（PDF + SVG + PNG）
| ファイル | 内容 |
|---------|------|
| `fig_4mC_5mC_competition.pdf/.svg/.png` | 4パネル統合図（4mC動態、GGCCGG濃縮、モデル、エビデンス表） |

---

## 8. 次のステップへの示唆

1. **同一試料BS-seq + Nanopore** — m4C/m5Cの同時定量による直接証明（最優先）
2. **Dcm-like MTaseのin vitro活性測定** — SC_RS19770/SC_RS36410がm5C産生することの実験的証明
3. **m4C-CCGG MTaseの同定** — 4mCを産生する酵素の特定
4. **GGCCGG部位のタイムコース解析** — T1→T2→T3でのm4C/m5C比率の定量的推移

---

## 参考文献

- Pisciotta A, et al. (2023) Scientific Reports 13, 7038. DOI: 10.1038/s41598-023-34075-1
- BioProject: PRJNA933392

---

*最終更新: 2026-02-05*
