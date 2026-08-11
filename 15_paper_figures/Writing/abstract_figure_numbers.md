# Abstract ↔ 図内数値 照合レポート

**作成日**: 2026-05-23  
**対象**: S. coelicolor M145 論文 — `05_abstract.md` の数値 vs 対応図PDF内数値

---

## 凡例

| 記号 | 意味 |
|---|---|
| ✅ | 一致（Abstract記載値 = 図内値） |
| 🔴 | 不一致（明確な数値差） |
| ⚠️ | 要注意（文脈不一致・図に明示なし） |
| ℹ️ | 図内に未表示（本文・表のみ） |
| 🚫 | **Abstractに存在しない**（ユーザー指定リストの誤り） |

---

## 1. 照合表（Abstract記載全数値）

### 1-A. 主要分類数値

| # | 数値 | Abstract記載 | 対応パネル | 図内表示値 | 状態 | 備考 |
|---|---|---|---|---|---|---|
| 1 | Exposed TFs（数） | **57** | `Figure2_RM_redistribution.pdf` panel C, D | 57 (panel C), n=57 (panel D) | ✅ | 一致 |
| 2 | Shielded TFs（数） | **998** | `Figure2_RM_redistribution.pdf` panel C | **994** (panel C) | 🔴 | panel C は n=994、panel D のみ "n=962 of **998**" と表示。4の差 |
| 3 | Exposed TFs（割合） | **5.4%** | implied (57/1,055) | 明示なし | ⚠️ | 計算値は5.40%で正確。図に% label なし |
| 4 | Shielded TFs（割合） | **94.6%（998/1,055）** | implied | 明示なし | ⚠️ | 計算値は94.60%。図に% label なし |

### 1-B. ROC / 分類境界

| # | 数値 | Abstract記載 | 対応パネル | 図内表示値 | 状態 | 備考 |
|---|---|---|---|---|---|---|
| 5 | AUC（pooled） | **0.917** | `Figure3_protection_zone.pdf` panel b | **0.923** | 🔴 | 既知の不一致（`panel_numbering_check.md` #6参照）。`Figure2c_shielded_exposed.pdf` panel b は **0.908** とさらに乖離 |
| 6 | 距離境界 | **293 bp** | `Figure3_protection_zone.pdf` panel b, `Figure2c_shielded_exposed.pdf` panel b, `Figure2_RM_redistribution.pdf` panel B, C | 293 bp（全PDF一致） | ✅ | 全図で一致 |

### 1-C. 保護ゾーン・地理的混同

| # | 数値 | Abstract記載 | 対応パネル | 図内表示値 | 状態 | 備考 |
|---|---|---|---|---|---|---|
| 7 | 保護ゾーン幅 | **~2,200 bp** | `Figure2_RM_redistribution.pdf` panel B（TSS metagene） | ラベル **なし** | ⚠️ | 図のTSSメタジーンプロファイルに "~2,200 bp" の明示ラベルはなく "293 bp threshold" のみ表示 |
| 8 | p値（unstratified） | **8.3 × 10⁻⁸** | `FigS9_negative_results.pdf` panel a | **8.3 × 10⁻⁸** | ✅ | 完全一致 |
| 9 | p値（core-only） | **0.87** | `FigS9_negative_results.pdf` panel a | **0.87** | ✅ | 完全一致 |
| 10 | Jaccard（T1–T2） | **0.000** | 図なし（Table S1参照） | Table S8 text: "Jaccard=0.000" | ℹ️ | 現行主要図には未表示。`Figure2_RM_redistribution.pdf` caption内にも記載なし |

### 1-D. 発現解析数値

| # | 数値 | Abstract記載 | 対応パネル | 図内表示値 | 状態 | 備考 |
|---|---|---|---|---|---|---|
| 11 | TetR p値 | **p = 0.028** | 図なし（本文のみ） | — | ℹ️ | `Figure5_exposed_TFs.pdf` panel d（TF family bar chart）に TetR 本数は示されるが p 値ラベルなし |
| 12 | T1–T2 rho | **0.136, p = 0.299** | 図なし（本文のみ） | — | ℹ️ | Negative results section テキストに記載あり。対応図パネルは現行稿に存在しない |
| 13 | T2–T3 rho | **0.012, p = 0.927** | 図なし（本文のみ） | — | ℹ️ | 同上 |

### 1-E. ユーザー指定リストの要精査項目

| # | 数値 | 記載場所 | 判定 | 備考 |
|---|---|---|---|---|
| 14 | Fisher OR = 145,448（95% CI: 36,131–585,503） | **Methods** (`02_methods.md`) | 🚫 **Abstract に記載なし** | co-modification enrichment の OR。Methods の計算根拠として記載。図にも未表示 |
| 15 | p < 10⁻³⁰⁰ | **Abstract に記載なし** | 🚫 **Abstract に記載なし** | 原稿内にも検索されない（Results/Methods に相当する記述なし） |

---

## 2. 重大不一致サマリー

### 🔴 要修正（数値が一致しない）

#### 不一致 A — AUC 0.917 vs 図内 0.923 / 0.908

| 出典 | 値 |
|---|---|
| Abstract（`05_abstract.md`） | **0.917** |
| Results（`01_results.md` 行31, 39） | **0.917**（95% CI: 0.895–0.935） |
| `Figure3_protection_zone.pdf` panel b | **0.923** |
| `Figure2c_shielded_exposed.pdf` panel b | **0.908** |
| `Figure4_shielded_exposed.pdf` panel b（旧版） | **0.910** |

**影響**: Abstractの中心的主張の数値が全PDFと食い違う。最終バージョンのPDFで0.917を示す図を確定させる必要がある。

---

#### 不一致 B — Shielded TFs n=998 vs 図内 n=994

| 出典 | 値 |
|---|---|
| Abstract | **998** |
| Results（`01_results.md` 行37） | **998** |
| `Figure2_RM_redistribution.pdf` panel C（箱ひげ図） | **994** |
| `Figure2_RM_redistribution.pdf` panel D（発現変動量） | "962 of **998**"（998は一致） |
| `Figure2c_shielded_exposed.pdf`（OCR解析）| "(n=994)" と読まれる可能性あり |

**影響**: 同一PDFの2パネル間でも不一致。Panel C のn=994は古い集計値（4遺伝子分の処理条件変更の可能性）。

---

### ⚠️ 要確認（明示的ラベルなし）

#### 注意 C — ~2,200 bp が図中に未ラベル

- Abstract と Results 本文に「~2,200 bp methylation-free protection zone」として記載
- `Figure2_RM_redistribution.pdf` panel B（TSS metagene profile）が対応図だが、図中に「2,200 bp」の明示ラベルなし
- 「293 bp threshold」のみラベルあり
- **推奨**: panel B に ±1,100 bp の範囲を示すブレースまたは "~2,200 bp protection zone" のアノテーションを追加

#### 注意 D — Jaccard = 0.000 が主要図に未表示

- 図には未表示。現状は Results 本文記述と Table S8 のみ
- 現行の `Figure2_RM_redistribution.pdf` panel A（geographic redistribution）キャプションへの注記、または panel 内テキストとして追加を推奨

#### 注意 E — S9c の % 値が本文と微差

- `FigS9_negative_results.pdf` panel c: 露出配列あり **25.8%**、遮蔽配列あり **6.3%**
- 本文: **24.6%** (14/57)、**6.4%** (64/994)
- Abstract には含まれないが、Results 本文引用と補足図の間で軽微な不一致あり

---

## 3. ℹ️ Abstract数値のうち対応する図が存在しない項目

以下の数値は Abstract に記載されているが、現行の主要図（Fig. 1–5）では**数値ラベルとして表示されていない**。本文記述のみで根拠が示されている。

| 数値 | Abstract の文脈 | 推奨対応 |
|---|---|---|
| Jaccard = 0.000 | RM系の完全な地理的再分布 | `Figure2_RM_redistribution.pdf` panel A にテキスト注記追加 |
| ~2,200 bp | 保護ゾーンの幅 | `Figure2_RM_redistribution.pdf` panel B にブレース注記追加 |
| p = 0.028（TetR） | repression bloc の TetR 濃縮 | `Figure5_exposed_TFs.pdf` panel d に有意性ラベル追加（探索的である旨の注記付き） |
| rho = 0.136, p = 0.299 | メチル化変化と発現変化の無相関（T1–T2） | Scatter plot 図（現行未作成）を追加 |
| rho = 0.012, p = 0.927 | 同上（T2–T3） | 同上 |

---

## 4. 対応PDFファイル一覧

| 原稿図番号 | 対応PDFファイル（`figures/main/`） | 確認パネル |
|---|---|---|
| Fig. 1a–1d | `Figure1_methylation_landscape.pdf` | a, b, c, d |
| Fig. 2a | `Figure2_RM_redistribution.pdf` | A（大文字） |
| Fig. 2b | `Figure2_RM_redistribution.pdf` | B（TSS metagene） |
| Fig. 2c | `Figure2c_shielded_exposed.pdf` | a（距離分布）, b（ROC） |
| Fig. 3a–3e | `Figure3_protection_zone.pdf` | a, b, c, d, e |
| Fig. 4a–4e | `Figure5_exposed_TFs.pdf`（ファイル名に注意） | (a)–(e) |
| Fig. S9a–S9d | `figures/supplementary/FigS9_negative_results.pdf` | a, b, c, d |

> **注**: `Figure5_exposed_TFs.pdf` のファイル名は "Figure5" だが内容は原稿の **Fig. 4** に対応する（既知のファイル命名ミス、`panel_numbering_check.md` 不一致#1 参照）。

---

*作成: Claude (Cowork) — 2026-05-23*
