# 図パネル番号整合性チェックレポート

**作成日**: 2026-05-23  
**チェック対象**: S. coelicolor M145 論文 — `01_results.md` / `full_manuscript.md` vs. PDF図ファイル

---

## 1. 各PDFに実際に存在するパネル一覧

| PDFファイル | 実際のパネル | 備考 |
|---|---|---|
| `Figure1_methylation_landscape.pdf` | **a, b, c, d** | ✅ 4パネル確認 |
| `Figure2_RM_redistribution.pdf` | **A, B** (大文字) | 原稿の小文字 a, b に対して大文字表記 |
| `Figure2c_shielded_exposed.pdf` | **単一パネル**（旧 a=距離分布のみ） | 2026-07-01: ROC/AUC panel b を除去（AUC非提示方針）。孤児図（投稿本体は `Figure2_RM_redistribution.pdf` panel C を使用）。is_exposed も 62/989 に更新 |
| `Figure3_protection_zone.pdf` | **a, b, c, d, e** | ✅ 5パネル確認 |
| `Figure4_shielded_exposed.pdf` | **a, b, c** のみ | d, e は**存在しない** |
| `Figure5_exposed_TFs.pdf` | **(a), (b), (c), (d), (e)** | ⚠️ ファイル名は "Figure5" だが内容は原稿の "Fig. 4" に相当 |
| `Figure6_switch_model.pdf` | **a, b, c, d** | ✅ 4パネル確認 (d は未参照) |

---

## 2. 不一致箇所リスト

### 🔴 重大な不一致（修正必須）

---

#### 不一致 #1 — PDFファイル命名とパネル番号の根本的ズレ

**問題の構造:**

| 状況 | 内容 |
|---|---|
| `Figure5_exposed_TFs.pdf` のパネル | (a) co-expression matrix, (b) bloc trajectories, (c) phase scatter, (d) TF family, (e) TCS asymmetry |
| 原稿での対応図番号 | **Fig. 4a–4e** （両manuscript共通） |
| `Figure4_shielded_exposed.pdf` のパネル | (a) promoter methylation exposure, (b) ROC, (c) expression-independence |
| 原稿での対応図番号 | **対応なし** — 原稿中のFig. 4参照はco-expression内容（= Figure5ファイル）を指す |

**影響する参照すべて:**

| ファイル | 行 | 本文の引用 | 実際のパネル所在 | 状態 |
|---|---|---|---|---|
| `01_results.md` | 47 | Fig. 4a, 4b (co-expression) | `Figure5_exposed_TFs.pdf` panels (a),(b) | 🔴 ファイル名不一致 |
| `01_results.md` | 51 | (Fig. 4b) activation bloc | `Figure5_exposed_TFs.pdf` panel (b) | 🔴 ファイル名不一致 |
| `01_results.md` | 53 | (Fig. 4a) repression bloc | `Figure5_exposed_TFs.pdf` panel (a) | 🔴 ファイル名不一致 |
| `01_results.md` | 55 | (Fig. 4e) TCS asymmetry | `Figure5_exposed_TFs.pdf` panel (e) | 🔴 ファイル名不一致 |
| `01_results.md` | 63 | (Fig. 4c) phase scatter | `Figure5_exposed_TFs.pdf` panel (c) | 🔴 ファイル名不一致 |
| `full_manuscript.md` | 69 | (Fig. 4a) co-expression | `Figure5_exposed_TFs.pdf` panel (a) | 🔴 ファイル名不一致 |
| `full_manuscript.md` | 73 | (Fig. 4d, Fig. 6b) | `Figure5_exposed_TFs.pdf` panel (d) | 🔴 ファイル名不一致 |
| `full_manuscript.md` | 75 | (Fig. 4d, Fig. 6a) | `Figure5_exposed_TFs.pdf` panel (d) | 🔴 ファイル名不一致 |
| `full_manuscript.md` | 77 | (Fig. 4e) TCS asymmetry | `Figure5_exposed_TFs.pdf` panel (e) | 🔴 ファイル名不一致 |
| `full_manuscript.md` | 83 | (Fig. 4c) phase scatter | `Figure5_exposed_TFs.pdf` panel (c) | 🔴 ファイル名不一致 |
| `full_manuscript.md` | 85 | (Fig. 4b), (Fig. 4c) | `Figure5_exposed_TFs.pdf` panels (b),(c) | 🔴 ファイル名不一致 |

**推奨修正:**  
`Figure5_exposed_TFs.pdf` → `Figure4_exposed_TFs.pdf` にリネーム（またはFig. 4の参照をすべてFig. 5に変更）。  
`Figure4_shielded_exposed.pdf` は現在の原稿に対応する参照がないため、archiveに移動またはFig. 3の旧バージョンとして整理。

---

#### 不一致 #2 — `01_results.md` 行61: Fig. 3c の内容記述が誤り

| 項目 | 内容 |
|---|---|
| **ファイル** | `01_results.md` |
| **行** | 61 |
| **本文引用** | `Fig. 3c` |
| **本文の記述** | "Figure 3c plots the **two transition magnitudes** directly (x-axis: \|LFC T1→T2\|; y-axis: \|LFC T2→T3\|), with the phase ratio = 0.6 threshold shown as the dashed line" |
| **実際の Fig. 3c** | `Figure3_protection_zone.pdf` panel c = **Feature discriminative power**（AUCバーチャート: Distance 0.923, Sequence 0.589, GC% 0.605, Expression 0.556） |
| **phase scatter の実在** | `Figure5_exposed_TFs.pdf` panel (c) = "Phase scatter (T1→T2 vs T2→T3)" |
| **状態** | 🔴 **不一致** — 記述内容と図パネルが合っていない |

**推奨修正:**  
`01_results.md` 行61の `Fig. 3c` → `Fig. 4c` に変更（phase scatterはFig. 4c相当）。  
cf. `full_manuscript.md` 行83では同内容に `(Fig. 6c)` および `(Fig. 4c)` を使用しており、01_results.mdは修正が遅れている。

---

#### 不一致 #3 — `full_manuscript.md` 行37: Fig. 5a/5b が内容と不一致

| 項目 | 内容 |
|---|---|
| **ファイル** | `full_manuscript.md` |
| **行** | 37 |
| **本文引用** | `(Fig. 5a)` および `(Fig. 5b)` |
| **本文の内容** | Simpson's paradox: GCCGGC地理的confoundの消失 (5a)、AAGCCCG null result (5b) |
| **実際の Fig. 5a** | `Figure5_exposed_TFs.pdf` panel (a) = **Spearman co-expression matrix** (57 exposed TFs) |
| **実際の Fig. 5b** | `Figure5_exposed_TFs.pdf` panel (b) = **Bloc expression trajectories** |
| **01_results.md の同内容** | 行13 では `(Fig. S9a)` `(Fig. S9b)` として補足図参照 |
| **状態** | 🔴 **不一致** — Fig. 5a/5b の内容がファイルと完全不一致、かつ01_results.mdとも矛盾 |

**推奨修正:**  
`full_manuscript.md` 行37の `(Fig. 5a)` → `(Fig. S9a)`、`(Fig. 5b)` → `(Fig. S9b)` に戻すか、  
またはSimpson's paradox用の新しいメインFigureを用意して明確に番号を割り振る。

---

#### 不一致 #4 — `full_manuscript.md` 内でFig. 5が2種類の内容を指す矛盾

| 項目 | 内容 |
|---|---|
| **ファイル** | `full_manuscript.md` |
| **行37** | `(Fig. 5a)` `(Fig. 5b)` = Simpson's paradox |
| **行103** | "Figure 5 presents an integrative **four-panel** summary ... Panel A, B, C, D ... (Fig. 5)" |
| **状態** | 🔴 **自己矛盾** — 同一の "Fig. 5" に、Simpson's paradox (a,b) と integrative summary (A-D) の2通りの内容が割り当てられている |

**推奨修正:**  
行37のFig. 5参照をFig. S9に戻し、行103のFig. 5をintegrative summaryとして明確化する。

---

### 🟡 注意が必要な箇所（要確認）

---

#### 不一致 #5 — `01_results.md` 行51: Fig. 3b "see inset" がSCO1160に使われている

| 項目 | 内容 |
|---|---|
| **ファイル** | `01_results.md` |
| **行** | 51 |
| **本文引用** | `(Fig. 3b; see inset)` |
| **本文の内容** | "*SCO1160* (peak LFC = +6.8)... encoding a sensor kinase whose cognate response regulator is **shielded** (Fig. 3b; see inset)" |
| **実際のFig. 3b** | `Figure3_protection_zone.pdf` panel b = ROC curve (AUC = 0.923) でSCO1160の言及なし |
| **full_manuscript.md の同内容** | 行73: `(Fig. 4d, Fig. 6b)` — TF family compositionとswitch model |
| **状態** | 🟡 **要確認** — ROC figureにSCO1160のinsetがあるかPNG等で視覚確認が必要。なければ誤参照。 |

**推奨修正（insetが存在しない場合）:**  
`(Fig. 3b; see inset)` → `(Fig. 4d; Fig. 6b)` に変更（full_manuscript.mdに合わせる）。

---

#### 不一致 #6 — ROC AUC値の不一致（PDFと本文で異なる）

| 出典 | AUC値 |
|---|---|
| `01_results.md` 行31, 39 | **0.917** (95% CI: 0.895–0.935) |
| `full_manuscript.md` 行61 | **0.917** |
| `Figure2c_shielded_exposed.pdf` panel b | 0.908 |
| `Figure3_protection_zone.pdf` panel b | 0.923 |
| `Figure4_shielded_exposed.pdf` panel b | 0.910 |
| **状態** | 🟡 **要確認** — 3つのPDFすべてが原稿値0.917と異なる |

**推奨修正:**  
最新の計算結果でFig. 3b (ROC) のPDFを再生成し、0.917を示すバージョンをfinalとして使用。または原稿の数値を最新の計算値に揃える。

---

#### 不一致 #7 — `Figure2_RM_redistribution.pdf` のパネルラベルが大文字

| 項目 | 内容 |
|---|---|
| **PDFのパネル表記** | **A** (Geographic redistribution) / **B** (TSS metagene profile) |
| **原稿の表記** | `(Fig. 2a)` / `(Fig. 2b)` (小文字) |
| **状態** | 🟡 **軽微な不一致** — 内容は一致するが表記統一が必要 |

---

### ✅ 一致が確認された箇所

| 参照 | 対応するPDF/パネル | 状態 |
|---|---|---|
| Fig. 1a (landscape) | `Figure1_methylation_landscape.pdf` panel a | ✅ |
| Fig. 1b (site counts) | `Figure1_methylation_landscape.pdf` panel b | ✅ |
| Fig. 1c (genomic distribution) | `Figure1_methylation_landscape.pdf` panel c | ✅ |
| Fig. 1d (MEME motif) | `Figure1_methylation_landscape.pdf` panel d | ✅ |
| Fig. 2a (geographic redistribution) | `Figure2_RM_redistribution.pdf` panel A | ✅ (大文字注意) |
| Fig. 2b (TSS metagene) | `Figure2_RM_redistribution.pdf` panel B | ✅ (大文字注意) |
| Fig. 2c (shielded/exposed) | `Figure2c_shielded_exposed.pdf` | ✅ |
| Fig. 3a (methylation at TSS) | `Figure3_protection_zone.pdf` panel a | ✅ |
| Fig. 3b (ROC) | `Figure3_protection_zone.pdf` panel b | ✅ (AUC値要確認) |
| Fig. 3c (feature AUC) | `Figure3_protection_zone.pdf` panel c | ✅ |
| Fig. 3d (TFBS profile) | `Figure3_protection_zone.pdf` panel d | ✅ |
| Fig. 3e (expression-independence) | `Figure3_protection_zone.pdf` panel e | ✅ |
| Fig. 6a/6b/6c (switch model) | `Figure6_switch_model.pdf` panels a, b, c | ✅ |

---

## 3. 修正優先度サマリー

| 優先度 | 問題 | 推奨アクション |
|---|---|---|
| 🔴 最高 | `Figure5_exposed_TFs.pdf` がFig. 4内容を持つ（ファイル命名ミス） | PDFをFigure4としてリネーム or 全参照をFig. 5に統一 |
| 🔴 最高 | `01_results.md` 行61: Fig. 3c = phase scatter（誤） | → Fig. 4c に修正 |
| 🔴 最高 | `full_manuscript.md` 行37: Fig. 5a/5b = Simpson's paradox（Fig. 5内容と矛盾） | → Fig. S9a/S9b に戻す |
| 🔴 最高 | `full_manuscript.md` 内でFig. 5が2種類の内容を指す自己矛盾 | Fig. 5の用途を一本化 |
| 🟡 中 | `01_results.md` 行51: Fig. 3b "see inset" でSCO1160参照（視覚的inset要確認） | PNG目視確認後、必要なら → (Fig. 4d, Fig. 6b) |
| 🟡 中 | ROC AUC値: 原稿0.917 vs PDFs 0.908/0.910/0.923 | 最新計算結果でPDF再生成 |
| 🟡 低 | `Figure2_RM_redistribution.pdf` パネル大文字A/B | PDF再生成時に小文字a/bに統一 |
| ℹ️ 情報 | `Figure4_shielded_exposed.pdf` が原稿参照に対応しない（旧バージョン） | archiveディレクトリへ移動を推奨 |

---

## 4. 参照整理: 両ファイル間の相違点

| 内容 | `01_results.md` | `full_manuscript.md` |
|---|---|---|
| Simpson's paradox | Fig. S9a, S9b | **Fig. 5a, 5b** 🔴 |
| SCO1160 sensor kinase | Fig. 3b; see inset | **Fig. 4d, Fig. 6b** |
| Phase scatter (timing) | Fig. 3c (行61) / Fig. 4c (行63) | Fig. 6c / Fig. 4c |
| TCS asymmetry | Fig. 4e | Fig. 4e |
| Integrative summary | Fig. 5 (4パネル) | Fig. 5 (矛盾あり) |
| Switch model blocs | 記載なし (Fig. 6なし) | Fig. 6a, 6b, 6c |

