# カラー・Supplementary Figure チェックレポート

生成日: 2026-05-23  
対象: S. coelicolor M145 methylation 論文 (15_paper_figures)

---

## D: カラー・記号の統一性チェック

### D-1. Shielded / Exposed カラー

**正規定義 (`00_shared_utils.py` lines 38–39):**
- `COL_EXPOSED  = '#7E57C2'`（紫）
- `COL_SHIELDED = '#B0BEC5'`（グレー）

#### スクリプト別カラー使用状況

| スクリプト | 対応図 | COL_EXPOSED | COL_SHIELDED | 統一? |
|---|---|---|---|---|
| `00_shared_utils.py` | (canonical) | `#7E57C2` 紫 | `#B0BEC5` グレー | ✅ 正規 |
| `10_new_figure1_overview.py` (panel_d) | Fig 1 panel D | **`#E74C3C` 赤** | **`#999999` グレー** | ⚠️ **不整合** |
| `11_new_figure2_protection.py` | Fig 2/3 | shared_utils経由 #7E57C2 | shared_utils経由 #B0BEC5 | ✅ |
| `12_new_figure3_exposed_TFs.py` | Fig 3/5 (Exposed TFs) | shared_utils経由 #7E57C2 | shared_utils経由 #B0BEC5 | ✅ |
| `13_new_figure4_negative_results.py` | Fig 4 | shared_utils経由 #7E57C2 | shared_utils経由 #B0BEC5 | ✅ |
| `14_new_figure5_switch_model.py` | Fig 5/6 (switch) | shared_utils経由 #7E57C2 | shared_utils経由 #B0BEC5 | ✅ |
| `02d_figure2c_shielded_exposed_redesign.py` | Fig 2c | `#7E57C2` (独自定義) | `#B0BEC5` (独自定義) | ✅ |
| `run_figure5_n57.py` | Fig 5 (n=57) | `#7E57C2` (独自定義) | `#B0BEC5` (独自定義) | ✅ |
| `05_figure5_gatekeeper_model.py` | Figure5_gatekeeper_model (旧?) | **`#F44336` 赤** | — | ⚠️ **不整合 (旧スクリプト)** |

#### 🔴 問題点 (Shielded/Exposed色)

1. **`10_new_figure1_overview.py` の `panel_d` 関数 (line 284–285)**  
   - `col_exposed = '#E74C3C'`（赤 = Tomato）を使用  
   - `col_shielded = '#999999'`（canonical `#B0BEC5` と異なるグレー）  
   - この関数は `shared_utils` をインポート済みだが、`COL_EXPOSED` を使わず**ローカル上書き**している  
   - → **Figure 1 panel D のみ Exposed が赤、他は紫** → 読者に混乱を与える可能性大

2. **`05_figure5_gatekeeper_model.py` (line 44)**  
   - `COL_EXPOSED = '#F44336'`（赤）  
   - `Figure5_gatekeeper_model.png` が main/ に存在するため現用図に影響あり  
   - → **Figure 5 gatekeeper panel でも Exposed が赤**

---

### D-2. T1 / T2 / T3 タイムポイントカラー

**⚠️ shared_utils に T1/T2/T3 の正規カラー定義なし** — 各スクリプトが独自に定義。

| スクリプト | 対応図 | T1 | T2 | T3 | 整合性 |
|---|---|---|---|---|---|
| `02b_figure2_RM_redistribution.py` | Fig 2 panel A | `#E53935` 赤 | `#FB8C00` 橙 | `#1565C0` 濃青 | — |
| `21_figS_bgc_methylation_track.py` | FigS_BGC track | **`#43A047` 緑** | `#FB8C00` 橙 | `#1E88E5` 青 | ⚠️ T1が緑 |
| `23_figS_bgc_integrated_track.py` | FigS_BGC integrated | **`#43A047` 緑** | `#FB8C00` 橙 | `#1E88E5` 青 | ⚠️ T1が緑 |
| `25_figS15_window_methylation.py` | FigS15 | `#4C72B0` 青 | `#DD8452` 橙 | `#55A868` 緑 | ⚠️ Seabornデフォルト, T1=青 |
| `26_figS16_timepoint_TSS.py` | FigS16 | `#4C72B0` 青 | `#DD8452` 橙 | `#55A868` 緑 | ⚠️ Seabornデフォルト, T1=青 |
| `07_supplementary_figures.py` (旧) | 旧 Supp | `#FFCDD2` 薄ピンク | `#C5CAE9` 薄紫 | `#C8E6C9` 薄緑 | ⚠️ 旧スクリプト |
| `20_figS12_protection_controls.py` | FigS12 | **`COL_EXPOSED`=#7E57C2** 紫 | `#FB8C00` 橙 | — | ⚠️ T1にExposed色を流用 |

#### 🔴 問題点 (T1/T2/T3色)

- T1/T2/T3 に統一カラーが存在しない。T1だけで赤・緑・青・紫 の4種類が混在。
- 特に `20_figS12` で T1 = `COL_EXPOSED`（紫）を流用しているのは意味的に混乱を招く。
- **推奨対応**: `00_shared_utils.py` に `TP_COLORS = {'T1': ..., 'T2': ..., 'T3': ...}` を追加し、全スクリプトで参照する。

---

### D-3. 記号 (Marker) 統一性

| スクリプト | Exposed marker | Shielded marker | 統一? |
|---|---|---|---|
| `12_new_figure3_exposed_TFs.py` | `'D'` ◆ ダイヤモンド | `'o'` ● 丸 | ✅ |
| `run_figure5_n57.py` | `'D'` ◆ ダイヤモンド | `'o'` ● 丸 | ✅ |

→ マーカーを明示しているスクリプト2件は統一されている。他の図では散布図やバープロットのため記号なし。

---

### D まとめ（カラー）

| 要素 | Fig 1 (panel D) | Fig 2/3 | Fig 3 (Exp TFs) | Fig 4 | Fig 5 (switch) | Fig 5 (gatekeeper) | 統一? |
|---|---|---|---|---|---|---|---|
| Exposed色 | **#E74C3C 赤** | #7E57C2 紫 | #7E57C2 紫 | #7E57C2 紫 | #7E57C2 紫 | **#F44336 赤** | ⚠️ Fig1-D・gatekeeper が不整合 |
| Shielded色 | **#999999 灰** | #B0BEC5 灰 | #B0BEC5 灰 | #B0BEC5 灰 | #B0BEC5 灰 | — | ⚠️ Fig1-D が微妙に違う |
| T1色 | (なし) | #E53935 赤 | — | — | — | (なし) | ⚠️ スクリプト間で不統一 |
| T2色 | (なし) | #FB8C00 橙 | — | — | — | (なし) | △ T2橙は複数一致 |
| T3色 | (なし) | #1565C0 濃青 | — | — | — | (なし) | ⚠️ スクリプト間で不統一 |
| Exposed marker | — | — | ◆ 'D' | — | ◆ 'D' | — | ✅ |
| Shielded marker | — | — | ● 'o' | — | ● 'o' | — | ✅ |

---

## E: Supplementary Figure 存在確認

### E-1. 本文中の引用 vs ファイル存在

本文（`full_manuscript.md` / `01_results.md`）から抽出された Supplementary Figure 引用と対応ファイルの存在確認：

| 引用 | 対応ファイル（supplementary/） | 状態 |
|---|---|---|
| Fig. S1 | FigS1_motif_landscape.pdf / .svg | ✅ |
| Fig. S2 | FigS2_simpsons_paradox_detail.pdf / .svg | ✅ |
| Fig. S3 | FigS3_regulatory_avoidance_CMH.pdf / .svg | ✅ |
| Fig. S4 | FigS4_sequence_motif_depletion.pdf / .svg | ✅ |
| Fig. S7 | FigS7_conservation_metrics.pdf / .svg | ✅ |
| Fig. S9 | FigS9_negative_results.pdf / .svg | ✅ |
| Fig. S9a / b / c / d | FigS9_negative_results.pdf（サブパネル） | ✅ (単一ファイル内) |
| Fig. S14a / b / c / d | FigS14_motif_reliability.pdf / .svg | ✅ (単一ファイル内) |
| Fig. S15 | FigS15_window_methylation.pdf / .svg | ✅ |
| Fig. S16 / S16c | FigS16_timepoint_TSS_protection.pdf / .svg | ✅ (S16c はサブパネル) |

### E-2. ファイルは存在するが本文（results/full）に引用なし

| ファイル | 状態 | 備考 |
|---|---|---|
| FigS5_exposed_TF_annotation.pdf | 🟡 引用なし（outline には記載あり） | `00_outline.md` の Fig. S5 に対応 |
| FigS6_TCS_pair_analysis.pdf | 🟡 引用なし（outline には記載あり） | `00_outline.md` の Fig. S6 に対応 |
| FigS8_protection_zone.pdf | 🟡 引用なし（どのmd文書にも Fig. S8 参照なし） | — |
| FigS10_per_motif_temporal_regional.pdf | 🟡 引用なし | — |
| FigS11_bgc_methylation.pdf | 🟡 引用なし | — |
| FigS12_protection_controls.pdf | 🟡 引用なし | — |
| FigS13_comprehensive_motif_analysis.pdf | 🟡 引用なし | — |
| FigS_act/cda/cpk/red_methylation_track.pdf (×4) | 🟡 引用なし | BGCトラック図 |
| FigS_act/cda/cpk/red_integrated_track.pdf (×4) | 🟡 引用なし | BGCトラック図 |

### E-3. 命名コンフリクト

| 問題 | 詳細 |
|---|---|
| ⚠️ FigS8 と FigureS8 の二重存在 | `FigS8_protection_zone.pdf`（旧）と `FigureS8_enrichment_bubble.pdf`（新・PNG付）が両方 supplementary/ に存在。どちらが正式か要確認。 |

### E-4. プレースホルダー確認

```
grep "Figure XX" → 結果なし ✅
grep "Fig. XX"   → 結果なし ✅
```
→ **プレースホルダーは残存していない。**

---

## 総合サマリー

### 優先度 高 🔴

1. **`10_new_figure1_overview.py` panel_d のカラー修正**  
   - `col_exposed = '#E74C3C'` → `COL_EXPOSED`（`#7E57C2`）に変更  
   - `col_shielded = '#999999'` → `COL_SHIELDED`（`#B0BEC5`）に変更  
   - 現状 Figure 1 のみ Exposed が赤で他の全図と異なる

2. **`05_figure5_gatekeeper_model.py` の Exposed カラー修正**  
   - `COL_EXPOSED = '#F44336'` → `#7E57C2` に変更（または shared_utils import 追加）  
   - `Figure5_gatekeeper_model.png` が main/ に存在

3. **FigS8 命名コンフリクトの解消**  
   - `FigS8_protection_zone.pdf` と `FigureS8_enrichment_bubble.pdf` の二重存在を整理  
   - 正式な S8 を決定し、不要な方を archive/ へ移動

### 優先度 中 🟡

4. **T1/T2/T3 カラーを `shared_utils` に正規化**  
   - `00_shared_utils.py` に `TP_COLORS = {'T1': '...', 'T2': '...', 'T3': '...'}` を追加  
   - `20_figS12`、`21/23_figS_bgc`、`25/26_figS15/16` を統一

5. **引用されていない Supp Figs (S5, S6, S8, S10–S13, BGC tracks) の扱いを確認**  
   - Methods やCaption リストで言及しているか確認  
   - 不要ならば archive/ へ移動

### 優先度 低 🟢

6. **`20_figS12` での T1 = `COL_EXPOSED` 使用の意味確認**  
   - T1 にタイムポイントカラーではなく Exposed カラーを意図的に使用しているか要確認

