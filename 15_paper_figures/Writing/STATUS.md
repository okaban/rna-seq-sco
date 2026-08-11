# Figure Generation Status

Last updated: 2026-05-23

---

## n=57 Update (2026-05-23)

Exposed TF数を旧値 n=62 → 正しい値 **n=57**（Shielded: n=998）に修正。
以下の図とデータファイルを再生成・更新した。

### データファイル更新

| ファイル | 変更内容 |
|---|---|
| `52_shielded_exposed_boundary/tables/all_genes_features_unified_n57.tsv` | 正データ（1055 genes, 57 exposed / 998 shielded）— 変更なし（既に正） |
| `52_shielded_exposed_boundary/tables/expression_quintile.tsv` | n=62→n=57 で再生成。各五分位 exposed 合計 = 57 ✓ |
| `55_exposed_regulatory_module/tables/coexpression_matrix.tsv` | 57×57 — 変更なし（既に正） |
| `57_temporal_dynamics_exposed_TF/tables/temporal_classification.tsv` | 57 TFs — 変更なし（既に正） |

### スクリプト修正

| スクリプト | 変更内容 |
|---|---|
| `27_figure4_shielded_exposed.py` | 293bp operating point: spec 0.806→0.802, annotation 0.81→0.80; baseMean AUC 0.547→0.543 |
| `12_new_figure3_exposed_TFs.py` | タイトル「62」→「57」, 出力先 `new_Figure3_exposed_TFs`→`Figure5_exposed_TFs`, all_genes_features を n57 版に |
| `run_figure5_n57.py` | 新規スタンドアローンスクリプト（Figure5 高速生成用） |
| `make_fig4_png.py` | Figure 4 PNG 生成用スタンドアローンスクリプト |

### 生成ファイル (figures/main/)

| ファイル | 生成日時 | 内容 |
|---|---|---|
| `Figure2_RM_redistribution.pdf/png/svg` | 2026-05-09 | パネルC(バイオリン): Shielded n=998 / Exposed n=57; パネルD: |LFC| variability。スクリプトは `all_genes_features_unified_n57.tsv` を参照済み。 |
| `Figure4_shielded_exposed.pdf/png/svg` | **2026-05-23** | パネルA: バイオリン Shielded n=994 / Exposed n=57; パネルB: ROC AUC=0.910, 293bp spec=0.80; パネルC: 発現五分位 n=57 |
| `Figure5_exposed_TFs.pdf/png/svg` | **2026-05-23** | パネルA: 57×57 co-expression heatmap; パネルB: 時系列軌跡 (up n=34/down n=23); パネルC/D/E: 更新済み |

### 数値確認（n=57 データ）

| 項目 | 値 | 論文記載値 | 状態 |
|---|---|---|---|
| Exposed TF | 57 | 57 | ✓ |
| Shielded TF | 998 | 998 | ✓ |
| 総TF | 1,055 | 1,055 | ✓ |
| Exposed 割合 | 5.4% | 5.4% | ✓ |
| 293bp 感度 | 1.000 | 1.00 | ✓ |
| 293bp 特異度 | 0.800 | 0.81* | ✓ (*NaN除外で 994/998 genes 使用) |
| Distance AUC (full ROC) | 0.910 | 0.917† | ✓ (†論文値は cross-val mean=0.9166) |
| Expression AUC | 0.543 | 0.547→0.543 | ✓ 更新済み |
| Quintile 合計 Exposed | 57 | 57 | ✓ |

### 注意事項

- Figure 4 パネルA/B/C は `nearest_methyl_distance` の NaN を除外するため n=1051（Exposed 57 / Shielded 994）を使用。全体の分類は n=1055（57/998）。
- Figure 5 の `temporal_classification.tsv`（現行）は bloc="activation" が 56、"unassigned" が 1。n=62 時代の arc (35 act / 26 rep) と構造が異なる。`direction` 列（up=34/down=23）を用いて heatmap の順序を代替した。
- AUC=0.917（論文）は 10-fold CV の mean AUC（`cross_validation.tsv` 参照）; 全データ ROC は 0.9096。

---

## 残作業

- [ ] Abstract の Exposed TF 数記述を n=57 に確認
- [ ] Fig 2 co-modification (Figure2_co_modification.pdf) — 現状確認要
- [ ] Fig 5 model description の最終調整
