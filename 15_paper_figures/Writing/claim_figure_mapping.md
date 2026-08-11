# 主張–図パネル 双方向チェック
*S. coelicolor* M145 メチローム論文  
対象: `manuscript/01_results.md` × `figures/main/`  
実施日: 2026-05-23

---

## 凡例

| 記号 | 意味 |
|---|---|
| ✅ | 対応あり・問題なし |
| ⚠️ | 軽微な不整合（バージョン差・数値ずれ等） |
| 🔴 | 重大な不整合（パネル欠損・孤立パネル・未引用図） |

---

## PDF ↔ 原稿 図番号 対応表

PDFファイルの命名番号と原稿の図番号が **ずれている**点に注意。

| PDFファイル名 | 主要コンテンツ | 原稿図番号（推定） |
|---|---|---|
| `Figure1_methylation_landscape.pdf` | 染色体全体のメチル化マップ, 4mC/6mA時系列, モチーフ, 組成 | **Fig. 1** (a,b,c,d) |
| `Figure2_RM_redistribution.pdf` | Core→Arm 再分布 (A), TSS metagene (B), Shielded/Exposed距離 (C) | **Fig. 2** (a,b,c) |
| `Figure3_protection_zone.pdf` | TSS脱メチル化プロファイル, ROC, 特徴量識別力, TFBS, 発現独立性 | **Fig. 3** (a,b,c,d,e) |
| `Figure5_exposed_TFs.pdf` | 57 TF共発現行列, 時系列軌跡, 応答タイミング, TFファミリー, TCS非対称 | **Fig. 4** (a,b,c,d,e) |
| `Figure6_switch_model.pdf` | Repression/Activationブロックのスイッチモデル | **未引用** |
| `Figure7_kegg_bubble.pdf` | KEGGパスウェイ富化バブルチャート | **未引用** |
| `Figure8_positional_stratification.pdf` | ゲノム位置別KEGG富化 | **未引用** |
| `SuppFig_protection_zone_temporal.pdf` | 時系列保護ゾーン変化 (a,b,c,d) | **Fig. S16** (a,b,c,d) |

---

## Forward チェック（主張 → 対応パネルが存在するか）

### Fig. 1

| 主張（本文より） | 引用 | PDFパネル | 状態 |
|---|---|---|---|
| "Both modification types showed dramatic temporal dynamics... 4mC sites declining from 1,289 at T1 to 407 at T2 and 21 at T3" | (Fig. 1a, 1b) | `Figure1_methylation_landscape.pdf` panels a, b ✓ | ✅ |
| "neither modification type showed preferential targeting of specific genomic features" | (Fig. 1c) | `Figure1_methylation_landscape.pdf` panel c ✓ | ✅ |
| "De novo motif discovery (MEME)... GCCGGC for 4mC and AAGCCCG for 6mA" | (Fig. 1d) | `Figure1_methylation_landscape.pdf` panel d ✓ | ✅ |
| "*SC_RS17645*... downregulated 4.6-fold... only MTase showing significant temporal change" | (Fig. 1, Table S3) | `Figure1_methylation_landscape.pdf` 全体 ✓ | ✅ |
| "5mC signal was 0.01%... modification is 4mC" | (Fig. S1) | 補足図（mainフォルダ外） | ✅（補足） |

### Fig. 2

| 主張 | 引用 | PDFパネル | 状態 |
|---|---|---|---|
| "GCCGGC 4mC system: dramatic geographic redistribution... core-to-arm shift" | (Fig. 2a) | `Figure2_RM_redistribution.pdf` panel A ✓ | ✅ |
| "Regulatory gene TSSs showed pronounced depletion zone ~2,200 bp" | (Fig. 2b) | `Figure2_RM_redistribution.pdf` panel B (TSS metagene) ✓ | ✅ |
| "57 (5.4%) exposed regulators... 998 (94.6%) shielded regulators" | (Fig. 2c) | `Figure2_RM_redistribution.pdf` panel C ✓ | ✅ |

### Fig. 3

| 主張 | 引用 | PDFパネル | 状態 |
|---|---|---|---|
| "Deepest depletion at +300 bp... 54.1% of flanking levels" | (Fig. 3a) | `Figure3_protection_zone.pdf` panel a ✓ | ✅ |
| "ROC: AUC = **0.917**, optimal threshold 293 bp, sens=1.00, spec=0.806" | (Fig. 3b) | `Figure3_protection_zone.pdf` panel b **AUC = 0.923** と表示 | ⚠️ **数値不一致** |
| "Combined logistic regression CV AUC = **0.712 ± 0.049**" | (Fig. 3c) | `Figure3_protection_zone.pdf` panel c は個別特徴量AUC棒グラフのみ（0.923, 0.589, 0.605, 0.556）。**0.712 の複合モデルは図に非表示** | 🔴 **コンテンツ不一致** |
| "GCCGGC sites enriched at TFBS (fold = 1.157, p = 0.006), no depletion" | (Fig. 3d) | `Figure3_protection_zone.pdf` panel d ✓ (fold=1.16 と表示) | ✅ |
| "Expression AUC = 0.547, JT test p = 0.730" | (Fig. 3e) | `Figure3_protection_zone.pdf` panel e ✓ (baseMean AUC=0.547, JT p=0.730) | ✅ |

### Fig. 4

| 主張 | 引用 | PDFパネル | 状態 |
|---|---|---|---|
| "Co-expression analysis revealed two antagonistic blocs" | (Fig. 4a, 4b) | `Figure5_exposed_TFs.pdf` panels a (co-expression matrix), b (temporal trajectories) ✓ | ✅ |
| "No significant temporal separation between blocs (Mann-Whitney p = 0.459)" | (Fig. 4c) | `Figure5_exposed_TFs.pdf` panel c (response timing scatter) ✓ | ✅ |
| *(Fig. 4d は本文で引用なし)* | — | `Figure5_exposed_TFs.pdf` panel d **存在**（TF family by functional category） | 🔴 **孤立パネル** |
| "7/7 identified TCS pairs showed asymmetry (one partner exposed)" | (Fig. 4e) | `Figure5_exposed_TFs.pdf` panel e ✓ (TCS pair asymmetry) | ✅ |

### Fig. 5

| 主張 | 引用 | PDFパネル | 状態 |
|---|---|---|---|
| "Figure 5 presents an integrative four-panel summary... Panel A (chromosome geography), B (protection zone), C (expression heatmap), D (4-layer model)" | (Fig. 5) | **統合済みPDFが存在しない**。`Figure5_exposed_TFs.pdf` はFig. 4相当のコンテンツ。`Figure5_gatekeeper_model.png` は1パネルのみ。 | 🔴 **Fig. 5 PDFが未組立** |

---

## Backward チェック（図パネル → 本文引用があるか）

### `Figure1_methylation_landscape.pdf`

| パネル | 本文引用 | 状態 |
|---|---|---|
| a | "(Fig. 1a, 1b)" | ✅ |
| b | "(Fig. 1b)" | ✅ |
| c | "(Fig. 1c)" | ✅ |
| d | "(Fig. 1d)" | ✅ |

### `Figure2_RM_redistribution.pdf`

| パネル | 本文引用 | 状態 |
|---|---|---|
| A | "(Fig. 2a)" | ✅ |
| B | "(Fig. 2b)" | ✅ |
| C | "(Fig. 2c)" | ✅ |

### `Figure3_protection_zone.pdf`

| パネル | 本文引用 | 状態 |
|---|---|---|
| a | "(Fig. 3a)" | ✅ |
| b | "(Fig. 3b)" | ⚠️ AUC数値差: PDF=0.923, 本文=0.917 |
| c | "(Fig. 3c)" | 🔴 本文はlogistic regression AUC=0.712を示すと記載するが、PDFは個別特徴量AUCのみ |
| d | "(Fig. 3d)" | ✅ |
| e | "(Fig. 3e)" | ✅ |

### `Figure5_exposed_TFs.pdf`（原稿 Fig. 4 相当）

| パネル | 本文引用 | 状態 |
|---|---|---|
| a | "(Fig. 4a, 4b)" | ✅ |
| b | "(Fig. 4a, 4b)" | ✅ |
| c | "(Fig. 4c)" | ✅ |
| **d** | **引用なし** | 🔴 **孤立パネル**（TF family by functional category） |
| e | "(Fig. 4e)" | ✅ |

### `SuppFig_protection_zone_temporal.pdf`（原稿 Fig. S16 相当）

| パネル | 本文引用 | 状態 |
|---|---|---|
| a | "(Fig. S16)" | ✅ |
| b | "(Fig. S16)" | ✅ |
| c | "(Fig. S16c, d)" | ✅ |
| d | "(Fig. S16c, d)" | ✅ |

### 未引用 PDF（完全孤立）

| PDFファイル | パネル | 本文引用 | 内容 | 状態 |
|---|---|---|---|---|
| `Figure6_switch_model.pdf` | a,b,c,d | **なし** | Repression/Activationブロックのネットワーク図 | 🔴 **孤立図** |
| `Figure7_kegg_bubble.pdf` | (なし) | **なし** | KEGG富化バブルチャート | 🔴 **孤立図** |
| `Figure8_positional_stratification.pdf` | A,B | **なし** | ゲノム位置別KEGG | 🔴 **孤立図** |
| `Figure2_co_modification.pdf` | A,B | **なし** | 6mA/4mC共修飾OR分析 | ⚠️ 本文で内容言及あり（「921-fold enrichment」）が図引用なし |
| `Figure2c_shielded_exposed.pdf` | a,b | **なし** | Shielded/Exposed距離分布+ROC（旧バージョン） | ⚠️ Figure2_RM_redistribution.pdf と重複（旧版と推定） |
| `Figure4_shielded_exposed.pdf` | a,b,c | **なし** | 距離バイオリン+ROC+発現分位（AUC=0.910） | ⚠️ Figure3_protection_zone.pdf の旧版と推定（AUC数値差あり） |
| `Figure1C_SC_RS17645_domain.pdf` | (なし) | **なし** | SC_RS17645 ドメイン構造 | ⚠️ 本文で言及（「TRD specificity domain」等）だが図引用なし |
| `Figure1D_layer_independence.pdf` | (なし) | **なし** | Layer 1/2 重複解析（Jaccard=0.000） | ⚠️ 本文で言及（「Jaccard index = 0.000」）だが図引用なし |
| `Figure1_GCCGGC_arm_distribution.pdf` | a,b | **なし** | GCCGGC core/arm棒グラフ（旧版） | ⚠️ Figure2_RM_redistribution.pdf の旧版と推定 |
| `Figure1_panel_GCCGGC_distribution.pdf` | a,b | **なし** | GCCGGC 密度折れ線グラフ（旧版） | ⚠️ 同上 |

---

## 不整合サマリー（優先度別）

### 🔴 要対処（高優先度）

| # | 問題 | 詳細 |
|---|---|---|
| 1 | **Fig. 4d が孤立パネル** | `Figure5_exposed_TFs.pdf` panel d（TF family by functional category）が本文で一度も引用されていない。テキストを追加するか、パネルを削除する必要あり。 |
| 2 | **Fig. 3c のコンテンツ不一致** | 本文は「複合logistic regression CV AUC = 0.712 ± 0.049」を Fig. 3c で示すと記載。しかし実際のパネル c は個別特徴量AUC棒グラフのみで、0.712 という値が図に存在しない。図の更新 or 本文の修正が必要。 |
| 3 | **Fig. 5 の統合PDFが未作成** | 本文最終段落で「Figure 5 presents an integrative four-panel summary (A–D)」と記載しているが、対応する統合PDFが存在しない。`Figure5_gatekeeper_model.png`（推定panel D）のみ存在。Fig. 5を組み立てる必要あり。 |
| 4 | **Figures 6, 7, 8 が孤立** | PDFは完成しているが本文に引用が一切ない。これらを Results に組み込むか、削除するかを決定する必要あり。 |

### ⚠️ 要確認（中優先度）

| # | 問題 | 詳細 |
|---|---|---|
| 5 | **Fig. 3b の AUC 数値差** | 本文 0.917 vs. `Figure3_protection_zone.pdf` panel b の 0.923。どちらが最終版の解析結果か確認が必要。（Figure4_shielded_exposed.pdfは0.910を示しており三者が異なる） |
| 6 | **Figure2_co_modification.pdf が未引用** | 「921-fold co-modification enrichment」は本文で言及されているが "(Fig. Xn)" 形式の引用がない。補足図として引用するか削除するか判断が必要。 |
| 7 | **複数の旧バージョンPDFが残存** | `Figure2c_shielded_exposed.pdf`, `Figure4_shielded_exposed.pdf`, `Figure1_GCCGGC_arm_distribution.pdf`, `Figure1_panel_GCCGGC_distribution.pdf` は最新版のPDFと内容が重複している。`archive/` フォルダへの移動を推奨。 |
| 8 | **Figure1D, Figure1C が未引用** | Jaccard=0.000 の重複解析結果（1D）とSC_RS17645ドメイン図（1C）は本文で数値は述べられているが図引用がない。 |

---

## 補足：本文中の全図パネル引用リスト

```
本文引用 (grep 結果, 頻度降順):
2回: (Fig. S9b), (Fig. S15), (Fig. 3c)
1回: (Fig. S9d), (Fig. S9c), (Fig. S9a), (Fig. S9), (Fig. S7),
     (Fig. S4), (Fig. S3), (Fig. S16), (Fig. S14d)-(S14a),
     (Fig. S1), (Fig. 5), (Fig. 4e), (Fig. 4c), (Fig. 4b),
     (Fig. 4a), (Fig. 3e), (Fig. 3d), (Fig. 3b), (Fig. 3a),
     (Fig. 2c), (Fig. 2b), (Fig. 2a), (Fig. 1d), (Fig. 1c),
     (Fig. 1b)

括弧外での複合引用: "(Fig. 1a, 1b)", "(Fig. S16c, d)", "(Fig. S9a, Fig. ...)",
                    "(Fig. 3b; see inset)", "(Fig. 4a, 4b)"
```

---

*生成: Claude (Cowork mode) · 2026-05-23*
