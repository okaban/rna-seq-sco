# Figure Legend vs. 実図 整合性チェックレポート

**作成日**: 2026-05-23  
**チェック対象**: S. coelicolor M145 論文  
**Legendソース**: `manuscript/01_results.md` / `manuscript/full_manuscript.md`  
**図PDFソース**: `figures/main/` 内の各 Figure PDF  
**手法**: pdfplumber でPDFテキスト抽出 → 原稿記載値と照合

---

## 凡例（Legend）の所在について

本論文には独立した "Figure Legends" セクションが存在せず、図の説明は `01_results.md` および `full_manuscript.md` の本文中に埋め込まれている（各図への参照と統計値が本文内に記述）。以下では本文記載値を "Legend/原稿" として参照する。

---

## 不一致一覧

### 🔴 Category A — 数値の不一致（修正必須）

| 図 | パネル | 原稿の記載 | 実図の内容 | 状態 |
|---|---|---|---|---|
| **Fig. 1b** | 全体 site counts | "4mC: **3,294** high-confidence sites" (union) | 図内: "Unique positions: 4mC: **2,693**" | 🔴 不一致 |
| **Fig. 1b** | 全体 site counts | "6mA: **2,734** sites" (union) | 図内: "Unique positions: 6mA: **3,248**" | 🔴 不一致（4mCと6mAの大小関係も逆転） |
| **Fig. 2c / Fig. 3b / Fig. 4** | ROC AUC | AUC = **0.917** (95% CI: 0.895–0.935) | Figure2c: **0.908** / Figure3: **0.923** / Figure4: **0.910** | 🔴 全PDFで異なる値（0.917と一致するものなし） |
| **Fig. 2c / Fig. 3b** | Shielded n数 | n = **998** shielded | Figure2c: **994** / Figure2_RM: **994** / Figure4: **994** (panel D: "962 of 998") | 🔴 主要結果図で994（原稿は998） |
| **Fig. 3b / Fig. 4b** | Expression level AUC | baseMean AUC = **0.547** | Figure3: **0.444** / Figure4: **0.543** | 🔴 3つの値が全て異なる |
| **Fig. 5 / Fig. 6** | 早期応答遺伝子数 | "**63%** (39/57 genes)" are early responders | Figure5 panel c: "**36/57** early-dominant" | 🔴 39/57（68%）vs 36/57（63%）不一致 |
| **Fig. 5b** | 活性化bloc/抑制bloc遺伝子数 | "activation bloc (**~35** genes)" / "repression bloc (**~26** genes)" | Figure5 panel b: "Up (n=**34**) / Down (n=**23**)" (計57) | 🔴 Figure6は35/26で一致するがFigure5は34/23 |
| **Fig. 2 panel C** | Mann-Whitney p値 | *p* = **5.3 × 10⁻²⁸** | Figure2_RM_redistribution panel C: p = **1.7e-25** | 🔴 3桁違う |
| **Fig. 6 panel d** | Exposed TF数（Summary欄） | n = **57** exposed / n = **998** shielded | Figure6 summary: "**62** exposed, **993** shielded" | 🔴🔴 最重要：summary図が旧バージョンの数値 |

---

### 🟡 Category B — パネル内容記述の不一致（要確認）

| 図 | パネル | 原稿の記載 | 実図の内容 | 状態 |
|---|---|---|---|---|
| **Fig. 3b** | AUC comparison | "nearest distance AUC = 0.917; expression level AUC = 0.547" | Figure3: distance AUC=0.923, expression AUC=**0.444** | 🟡 上記Aと重複・分母も要確認 |
| **Fig. 2 panel B** | TSS metagene n数 | "7,646 expressed genes" | Figure2: "(n=2,646, all genes)" | 🟡 5,000遺伝子差。定義（全遺伝子 vs 発現遺伝子）の違いか要確認 |
| **Fig. 3c** | Feature AUC bar label | "~33% sequence / ~67% protein occupancy" | Figure3 panel c: "~33% sequence / ~67% protein occupancy" ✅ | 🟢 一致 |
| **Fig. 3d** | GCCGGC at TFBS fold | "GCCGGC fold = **1.157**" | Figure3 panel d: "GCCGGC fold = **1.16**" | 🟢 実質一致（小数点以下3桁 vs 2桁） |
| **Fig. 3e** | JT p値 | "JT *p* = **0.730**" | Figure3 panel e: "JT p = **0.730**" ✅ | 🟢 一致 |
| **Fig. 5c** | Phase ratio threshold | "phase ratio > **0.6**" → 早期応答 | Figure5 panel c: "Phase ratio = **0.6**" ✅ | 🟢 一致 |
| **Fig. 6c** | 同時スイッチ p値 | "Mann-Whitney *p* = **0.459**" | Figure6 panel c: "p = **0.459**" ✅ | 🟢 一致 |
| **Fig. 5e** | TCS pairs | "7/7 identified pairs" | Figure5 panel e: "7/7 identified pairs" ✅ | 🟢 一致 |

---

### ℹ️ Category C — ファイル命名/パネル番号の問題（既知・参照）

これらは `Writing/panel_numbering_check.md` で既に詳細に報告済み。要約のみ記載：

| 問題 | 概要 |
|---|---|
| `Figure5_exposed_TFs.pdf` | 内容は原稿の **Fig. 4**（co-expression）に相当、ファイル名と不一致 |
| `Figure4_shielded_exposed.pdf` | 現行原稿に対応する引用なし（旧バージョン） |
| `01_results.md` 行61 | "Fig. 3c" とあるが内容はphase scatter（Figure5 panel c に相当）→ 誤参照 |
| `full_manuscript.md` 行37 | "Fig. 5a/5b" = Simpson's paradox → `01_results.md` では "Fig. S9a/S9b" と矛盾 |

---

## 図別詳細チェック

### Figure 1 — Methylation Landscape (`Figure1_methylation_landscape.pdf`)

| 項目 | 原稿記載 | 図内容 | 状態 |
|---|---|---|---|
| 総4mCサイト数（union） | 3,294 sites | Unique positions: 2,693 | 🔴 |
| 総6mAサイト数（union） | 2,734 sites | Unique positions: 3,248 | 🔴（大小関係も逆） |
| T1 4mC sites（全修飾） | 1,289 (GCCGGC特異的) | 1,987（全4mC） | ℹ️ 定義の違い（後者はGCCGGC以外も含む） |
| T1 6mA sites | 1,934 | 1,934 | 🟢 |
| MEME motif (4mC) | "GCCGGC" | SAMGCCSGCCA (n=2,678) | 🟢 実質一致 |
| MEME motif (6mA) | "AAGCCCG" | GVSAAGCCCGVC (n=656) | 🟢 実質一致 |
| Genomic feature分布 | "promoter, coding, intergenic...proportion" | Fig 1c: Bonferroni補正あり、▲/▼表示 | 🟢 概念一致（詳細数値は補足） |
| 解析対象timepoints | T1 (12h), T2 (24h), T3 (50h) | T1 (12h), T2 (24h), T3 (50h) | 🟢 |

**注記**: `Figure1a` の landscape内にある per-timepoint サイト数（T1: 4mC=1,987 / 6mA=1,934）は**全4mC修飾**の数であり、原稿の「1,289 at T1」はGCCGGCモチーフ特異的サイト数。この定義の使い分けが原稿内で明示されているか要確認。

---

### Figure 1C — SC_RS17645 domain (`Figure1C_SC_RS17645_domain.pdf`)

| 項目 | 原稿記載 | 図内容 | 状態 |
|---|---|---|---|
| 遺伝子名 | SC_RS17645 (SCO3104) | "SC_RS17645 (SCO3104) · 679 aa" | 🟢 |
| ドメイン | HsdM-type MTase | SAM-dep. MTase superfamily + TRD-like specificity | 🟢 |
| Solo MTase根拠 | 原稿: "no cognate HsdR or HsdS homologue within genome" | 図: "Solo MTase no cognate HsdR or HsdS within 50 kb" | ⚠️ 原稿は「ゲノム全体で見つからない」、図は「50 kb以内で見つからない」と範囲が異なる |

---

### Figure 1D — Layer Independence (`Figure1D_layer_independence.pdf`)

| 項目 | 原稿記載 | 図内容 | 状態 |
|---|---|---|---|
| Jaccard index (T1-T2 GCCGGC) | 0.000 | "Jaccard = 0.000" | ✅ 一致するが文脈が異なる（図はLayer独立性、原稿はT1-T2サイト非重複） |
| p値 | 原稿未記載 | p = 6.4e-26 | ℹ️ 原稿に対応する記述なし |
| Universe (解析遺伝子数) | 未記載 | 762 genes | ℹ️ 原稿に数値なし |

---

### Figure 2 — RM Redistribution (`Figure2_RM_redistribution.pdf`)

| 項目 | 原稿記載 | 図内容 | 状態 |
|---|---|---|---|
| Panel A ラベル | 小文字 "a" | **大文字 "A"** | 🟡 表記不統一 |
| Panel B ラベル | 小文字 "b" | **大文字 "B"** | 🟡 表記不統一 |
| Panel B (TSS metagene) n数 | "7,646 expressed genes" | "(n=2,646, all genes)" | 🟡 定義差？要確認 |
| Panel C: Shielded n | 998 | **994** | 🔴 |
| Panel C: Mann-Whitney p値 | 5.3 × 10⁻²⁸ | **1.7e-25** | 🔴 |
| Panel C: 中央値 Shielded | 762 bp (median) | 747 bp (median) | 🔴 微少差（15 bp） |
| Panel D: Shielded total n | 998 | "962 **of 998**" (←図内に998と記載あり) | ⚠️ Panel Dは998と認識、Panel Cは994と表示（同一PDF内で不統一） |

---

### Figure 2c — Shielded/Exposed (`Figure2c_shielded_exposed.pdf`)

| 項目 | 原稿記載 | 図内容 | 状態 |
|---|---|---|---|
| Shielded n | **998** | **994** | 🔴 |
| Exposed n | 57 | 57 | 🟢 |
| ROC AUC | **0.917** | **0.908** | 🔴 |
| 感度（sensitivity） | 100% | sens = 1.00 | 🟢 |
| 特異度（specificity） | 80.6% | spec = 0.80 | 🟢 実質一致 |
| 閾値 | 293 bp | 293 bp | 🟢 |

---

### Figure 3 — Protection Zone (`Figure3_protection_zone.pdf`)

| 項目 | 原稿記載 | 図内容 | 状態 |
|---|---|---|---|
| Panel b: AUC (nearest distance) | **0.917** | **0.923** | 🔴 |
| Panel b: AUC (expression level) | **0.547** | **0.444** | 🔴 |
| Panel c: Sequence contribution | "~51%" (AUC-based) / "~33%" (motif-based) | "~33% sequence / ~67% protein occupancy" | 🟡 motif-based値を使用（AUC-based値は非表示） |
| Panel d: GCCGGC fold at TFBS | 1.157 | 1.16 | 🟢 |
| Panel e: JT p値 | 0.730 | 0.730 | 🟢 |
| Panel e: baseMean AUC | 0.547 | 0.547 | 🟢 |

**注記**: panel bのAUCが原稿(0.917)・Figure2c(0.908)・Figure3(0.923)・Figure4(0.910)で全て異なる。どの計算が最新かを確定し、全図・原稿を統一する必要がある。

---

### Figure 4 — Shielded/Exposed (旧バージョン) (`Figure4_shielded_exposed.pdf`)

| 項目 | 原稿記載 | 図内容 | 状態 |
|---|---|---|---|
| Panel a: Shielded n | 998 | **994** | 🔴 |
| Panel b: AUC (nearest) | 0.917 | **0.910** | 🔴 |
| Panel b: AUC (expression) | 0.547 | **0.543** | 🔴 |
| Panel c: JT p値 | 0.730 | 0.730 | 🟢 |
| Panel c: baseMean AUC | 0.543 | 0.543 | ⚠️ 図は0.543だが原稿テキストは0.547 |

**注記**: このPDFは現行原稿に対応する本文引用がない（旧バージョン）。`archive/` への移動を推奨。

---

### Figure 5 — Exposed TFs (`Figure5_exposed_TFs.pdf`)

> ⚠️ **注意**: このファイルの内容は原稿の **Fig. 4**（co-expression解析）に相当（`panel_numbering_check.md` 参照）。

| 項目 | 原稿記載 | 図内容 | 状態 |
|---|---|---|---|
| Panel a: 解析遺伝子数 | 57 exposed TFs | "57 exposed TFs" | 🟢 |
| Panel b: Activation遺伝子数 | **~35 genes** | Up (n=**34**) | 🔴 1遺伝子差 |
| Panel b: Repression遺伝子数 | **~26 genes** | Down (n=**23**) | 🔴 3遺伝子差（合計は34+23=57で一致） |
| Panel c: 早期応答遺伝子数 | **39/57 (63%)** | **36/57** early-dominant | 🔴 |
| Panel d: TF family composition | HTH enrichedなど | 図内バーチャートで確認 | 🟢 内容一致 |
| Panel e: TCS pairs | 7/7 | 7/7 | 🟢 |
| Panel e: SK:RR split | 4:3 | 図内: 4SK / 3RR | 🟢 |

---

### Figure 6 — Switch Model (`Figure6_switch_model.pdf`)

| 項目 | 原稿記載 | 図内容 | 状態 |
|---|---|---|---|
| Panel a/b: gene lists | 各遺伝子名は本文列挙 | 図内のgene list | 🟢 主要遺伝子一致（ramR, tcrA, fasR, sigJなど） |
| Panel c: Activation n | ~35 | **35** | 🟢 |
| Panel c: Repression n | ~26 | **26** | 🟢 |
| Panel c: simultaneous p値 | 0.459 | "p = 0.459, simultaneous" | 🟢 |
| Panel d (Summary): Exposed n | **57** | **62** | 🔴🔴 最重要不一致 |
| Panel d (Summary): Shielded n | **998** | **993** | 🔴 |
| Panel d (Summary): AUC | 0.917 | 0.917 | 🟢 |
| Panel d (Summary): threshold | 293 bp | 293 bp | 🟢 |
| Panel d (Summary): rho T1→T2 | 0.136, p=0.299 | "ρ = 0.136, p = 0.30" | 🟢 |
| Panel d (Summary): rho T2→T3 | 0.012, p=0.927 | "ρ = 0.012, p = 0.93" | 🟢 |

---

### Figure 7 — KEGG Bubble (`Figure7_kegg_bubble.pdf`)

| 項目 | 原稿記載 | 図内容 | 状態 |
|---|---|---|---|
| Quorum Sensing (GCCGGC) | 65/100 genes, fold=1.55, FDR=2.5×10⁻⁴ | FDR<0.05表示あり（QS, GCCGGC） | 🟢 概念一致 |
| Siderophore (Dual) | 5/7 genes, fold=12.2, FDR=8.7×10⁻⁴ | FDR<0.05表示あり（Siderophore, Dual） | 🟢 |
| AAGCCCG-proximal: enrichment | "no significant KEGG enrichment" | 図内 AAGCCCG列にFDR<0.05なし | 🟢 |
| n numbers in legend | GCCGGC n=3,197 / AAGCCCG n=459 / Dual n=448 | 図内未記載（凡例なし） | ℹ️ |

---

### Figure 8 — Positional Stratification (`Figure8_positional_stratification.pdf`)

| 項目 | 原稿記載 | 図内容 | 状態 |
|---|---|---|---|
| Panel A: AAGCCCG n | 1,078 sites (figure内) | n=1,078 | 🟢 |
| Panel A: GCCGGC n | 1,659 sites (figure内) | n=1,659 | 🟢 |
| Panel A: AAGCCCG 分布 | 82% CDS | 82% CDS | 🟢 |

---

## 優先度別 修正サマリー

| 優先度 | 問題 | 対象 | 推奨アクション |
|---|---|---|---|
| 🔴🔴 最最高 | **Figure 6 panel d の Exposed n = 62**（原稿は57） | `Figure6_switch_model.pdf` | Summaryテキストを n=57/n=998 に更新して再生成。これは投稿版に残ると即座にリジェクト理由になる |
| 🔴 最高 | **ROC AUC 値が4つの値に分散**（0.908/0.910/0.917/0.923） | Fig2c, Fig3, Fig4, 原稿全体 | 最新の計算結果（authoritative run）で1つの値に統一し、全PDF・全原稿を更新 |
| 🔴 最高 | **総サイト数 4mC 3,294 vs 図 2,693 / 6mA 2,734 vs 図 3,248**（大小関係も逆転） | Fig 1b, 原稿Abstract/Results | 原稿 or 図の定義（union方法）を統一し数値を修正 |
| 🔴 高 | **Shielded n = 998（原稿）vs 994（図）** | Fig2c, Fig2_RM, Fig4 | 正しい数値を確定し全PDFと原稿を統一 |
| 🔴 高 | **Mann-Whitney p値: 5.3×10⁻²⁸（原稿）vs 1.7e-25（図）** | Fig2_RM panel C | 再計算で正しいp値を確定し統一 |
| 🔴 高 | **Expression level AUC: 0.444/0.543/0.547 の3値が混在** | Fig3, Fig4, 原稿 | 正しい値に統一 |
| 🔴 高 | **早期応答遺伝子 39/57（原稿）vs 36/57（図）** | Fig5 panel c | 計算再確認し修正 |
| 🔴 中 | **Bloc sizes: 35/26（原稿・Fig6）vs 34/23（Fig5）** | Fig5 panel b | Fig5のラベルを35/26に修正 or 計算確認 |
| 🟡 中 | **Fig 1C: "50 kb以内"（図）vs "ゲノム全体"（原稿）のHsdR/HsdS非存在** | Fig1C, 原稿 | 範囲の記述を統一 |
| 🟡 中 | **Fig 2 panel B: n=2,646 vs 原稿の7,646 expressed genes** | Fig2_RM panel B | 定義の明確化 |
| 🟡 低 | **Fig 2 A/B パネルラベル大文字 vs 原稿小文字** | Fig2_RM_redistribution.pdf | PDF再生成時に a/b に統一 |
| ℹ️ | **Fig 3c: AUCバー表示が motif-based(33%)のみ（AUC-basedの51%は非表示）** | Fig3 panel c | 注釈追加を検討 |

---

## 一致が確認された重要統計値

| 統計値 | 原稿 | 図 | 状態 |
|---|---|---|---|
| GCCGGC sites T1 | 1,289 | 1,289 (Fig1_GCCGGC) | ✅ |
| GCCGGC sites T2 | 407 | 407 (Fig1_GCCGGC) | ✅ |
| GCCGGC sites T3 | 21 | 21 (Fig1_GCCGGC) | ✅ |
| AAGCCCG sites T1 | 260 | 260 (Fig1_panel) | ✅ |
| AAGCCCG sites T2 | 64 | 64 (Fig1_panel) | ✅ |
| T1 core fraction (GCCGGC) | 83% | 83% (Fig1_panel) | ✅ |
| T2 arm fraction (GCCGGC) | 82% | 82% (Fig1_panel) | ✅ |
| Chi-squared for redistribution | *p* ~ 10⁻¹³⁰ | (原稿のみ記載) | ✅ |
| co-modification OR | "921-fold" enrichment | "Fold-enrichment: 921×" (Fig2_co) | ✅ |
| co-modification Fisher OR | 145,448 (figure_claims_map) | OR = 145,448 (Fig2_co) | ✅ |
| ROC threshold | 293 bp | 293 bp (全図共通) | ✅ |
| Sensitivity at threshold | 100% | 1.00 | ✅ |
| Exposed n | 57 | 57（Fig6 summary以外全図） | ✅ |
| Simultaneous switch p | 0.459 | 0.459 (Fig6) | ✅ |
| rho (T1→T2) methylation-expression | 0.136, p=0.299 | 0.136, p=0.30 (Fig6) | ✅ |
| rho (T2→T3) methylation-expression | 0.012, p=0.927 | 0.012, p=0.93 (Fig6) | ✅ |
| JT p値（expression independence） | 0.730 | 0.730 (Fig3, Fig4) | ✅ |
| TFBS methylation fold | 1.157 | 1.16 (Fig3) | ✅ |
| 7/7 TCS asymmetry | 7/7 | 7/7 (Fig5) | ✅ |
| Jaccard index T1-T2 GCCGGC | 0.000 | 0.000 (Fig1D) | ✅ |
| Figure 1D: p値 | 未記載 | 6.4e-26 | ℹ️ |

---

## 重大不一致の根本原因（推定）

1. **ROC AUC / n数の分散**:  最も可能性が高い原因は、解析パラメータの変更（例: 高信頼度サイトのフィルタリング閾値の変更）が図の一部にのみ反映され、他の図・原稿が古いバージョンのまま残っていること。

2. **Figure 6 panel d の n=62/993**: Figure 6 が以前のバージョン（閾値や分類基準が異なる時期）のまま更新されていない。原稿の n=57 は最終確定値であり、この図だけが旧解析結果を表示している。

3. **総サイト数の逆転（4mC 3,294 vs 6mA 2,734 → 図は逆）**: 「高信頼度サイト」の定義変更か、修飾タイプ（4mC vs 6mA）のラベリングエラーの可能性。Early draftでは4mCと6mAの union 数が現在の原稿値と異なっていた可能性もある。

4. **早期応答遺伝子 39 vs 36**: phase ratio 閾値の適用方法またはデータセット（57遺伝子中の除外基準）の違いと推測。

---

*チェック実施者: Claude (Cowork mode) / 2026-05-23*  
*前回チェック (`panel_numbering_check.md`) との関係: 本レポートは数値整合性に特化。パネル番号問題は `panel_numbering_check.md` を参照のこと。*
