# Presentation Figures README

**作成日**: 2026-02-03
**用途**: NotebookLMへの読み込み用図表一覧

---

## 概要

このディレクトリには、プレゼンテーション原稿（`260203_presentation_manuscript_for_collaborators.md`）で言及した図表が格納されています。

**総ファイル数**: 68 PDFファイル

---

## PDFファイル一覧

### 1. QC・DESeq2解析

| ファイル名 | 内容 |
|-----------|------|
| `PCA_M145.pdf` | PCAプロット（PC1 83%, PC2 16%） |
| `sample_distance_heatmap_M145.pdf` | サンプル間ユークリッド距離ヒートマップ |
| `volcano_M145_3_vs_1.pdf` | Volcanoプロット（T3 vs T1） |
| `volcano_M145_2_vs_1.pdf` | Volcanoプロット（T2 vs T1） |
| `topVarGenes_heatmap_M145.pdf` | 上位100変動遺伝子ヒートマップ |

### 2. ベン図解析

| ファイル名 | 内容 |
|-----------|------|
| `venn_all_DEGs.pdf` | 全DEGsベン図（3比較の重複） |
| `venn_up_DEGs.pdf` | Up-regulated DEGsベン図 |
| `venn_down_DEGs.pdf` | Down-regulated DEGsベン図 |

### 3. BGC発現動態解析

| ファイル名 | 内容 |
|-----------|------|
| `BGC_timecourse_lineplot_M145.pdf` | BGCタイムコース（線形スケール） |
| `BGC_timecourse_lineplot_log10_M145.pdf` | BGCタイムコース（log10スケール） |
| `BGC_condition_heatmap_M145.pdf` | BGC×条件ヒートマップ |
| `BGC_gene_heatmap_act_M145.pdf` | Act BGC内遺伝子ヒートマップ |
| `BGC_gene_heatmap_red_M145.pdf` | Red BGC内遺伝子ヒートマップ |
| `BGC_gene_heatmap_cda_M145.pdf` | CDA BGC内遺伝子ヒートマップ |
| `BGC_gene_heatmap_cpk_M145.pdf` | Cpk BGC内遺伝子ヒートマップ |
| `BGC_sample_scores_boxplot_M145.pdf` | サンプルごとBGCスコア箱ひげ図 |
| `coverage_tracks_all_BGC.pdf` | 4 BGC統合カバレッジトラック |
| `coverage_track_act.pdf` | Act BGCカバレッジトラック |
| `coverage_track_red.pdf` | Red BGCカバレッジトラック |
| `coverage_track_cda.pdf` | CDA BGCカバレッジトラック |
| `coverage_track_cpk.pdf` | Cpk BGCカバレッジトラック |
| `report_BGC_expression_heatmap.pdf` | BGC発現動態ヒートマップ（論文用） |

### 4. KEGG/GO/COGエンリッチメント解析

| ファイル名 | 内容 |
|-----------|------|
| `KEGG_bubble_T3vsT1_up.pdf` | KEGG T3vsT1 Upバブルチャート |
| `KEGG_bubble_T3vsT1_down.pdf` | KEGG T3vsT1 Downバブルチャート |
| `KEGG_bubble_T2vsT1_up.pdf` | KEGG T2vsT1 Upバブルチャート |
| `KEGG_bubble_T2vsT1_down.pdf` | KEGG T2vsT1 Downバブルチャート |
| `GO_bubble_T3vsT1_up.pdf` | GO T3vsT1 Upバブルチャート |
| `GO_bubble_T3vsT1_down.pdf` | GO T3vsT1 Downバブルチャート |
| `GO_bar_T3vsT1_up.pdf` | GO T3vsT1 Upバーチャート |
| `GO_bar_T3vsT1_down.pdf` | GO T3vsT1 Downバーチャート |
| `COG_barplot_T3vsT1.pdf` | COG T3vsT1棒グラフ |
| `COG_barplot_T2vsT1.pdf` | COG T2vsT1棒グラフ |
| `COG_scatter_up_vs_down.pdf` | COG Up/Down散布図 |
| `report_COG_net_change.pdf` | COG Net変化水平棒グラフ（論文用） |
| `report_pathway_pattern_summary.pdf` | パスウェイ変動パターン分類 |
| `report_FE_change_patterns.pdf` | FE変化矢印図 |
| `report_key_insights_summary.pdf` | Key Insights俯瞰図 |

### 5. エピゲノム統合解析

| ファイル名 | 内容 |
|-----------|------|
| `Fig1_methylation_expression_correlation.pdf` | メチル化-発現相関散布図 |
| `Fig2_coordination_patterns.pdf` | 協調パターン分布 |
| `Fig3_volcano_methylation.pdf` | メチル化Volcanoプロット |
| `Fig4_summary_statistics.pdf` | 統計サマリー |
| `Fig5_top_genes_heatmap.pdf` | トップ遺伝子ヒートマップ |
| `venn_DEG_DMG_T2vsT1.pdf` | DEG-DMG重複ベン図（T2vsT1） |
| `venn_DEG_DMG_T3vsT1.pdf` | DEG-DMG重複ベン図（T3vsT1） |
| `venn_DEG_DMG_T3vsT2.pdf` | DEG-DMG重複ベン図（T3vsT2） |
| `contingency_tables.pdf` | 2×2分割表 |
| `coordination_summary.pdf` | 協調遺伝子サマリー |

### 6. R-Mシステム・モチーフ解析

| ファイル名 | 内容 |
|-----------|------|
| `4mC_base_composition.pdf` | 4mCサイト周辺塩基組成 |
| `6mA_base_composition.pdf` | 6mAサイト周辺塩基組成 |
| `motif_conservation.pdf` | モチーフ保存性解析 |
| `motif_heatmap.pdf` | モチーフヒートマップ |
| `m145_comparison_boxplot.pdf` | M145比較箱ひげ図 |

### 7. マルチオミクストラック・経時変化

| ファイル名 | 内容 |
|-----------|------|
| `multiomics_track_RamR.pdf` | RamRマルチオミクストラック |
| `multiomics_track_NsdB.pdf` | NsdBマルチオミクストラック |
| `multiomics_track_Act_SCO5079.pdf` | Act SCO5079トラック |
| `multiomics_track_Red_SCO5897.pdf` | Red SCO5897トラック |
| `multiomics_track_Cpk_SCO6284.pdf` | Cpk SCO6284トラック |
| `multiomics_overview_all_targets.pdf` | 全ターゲット統合 |
| `temporal_dynamics_RamR.pdf` | RamR経時変化 |
| `temporal_dynamics_NsdB.pdf` | NsdB経時変化 |
| `temporal_dynamics_all_genes.pdf` | 全遺伝子経時変化 |
| `methylation_expression_trajectories.pdf` | メチル化-発現軌跡 |

### 8. プロモーター構造

| ファイル名 | 内容 |
|-----------|------|
| `promoter_architecture_summary.pdf` | プロモーター構造サマリー |
| `motif_schematic_all_genes.pdf` | モチーフ模式図 |

### 9. 論文Main Figure

| ファイル名 | 内容 |
|-----------|------|
| `Figure1_overview.pdf` | 研究概要（論文Figure 1） |
| `Figure2_AAGCCCG_system.pdf` | AAGCCCG R-Mシステム（論文Figure 2） |
| `Figure3_GRN_methylation.pdf` | GRNメチル化（論文Figure 3） |
| `Figure4_Act_vs_Red.pdf` | Act vs Red比較（論文Figure 4） |

---

## NotebookLM読み込み用ファイル

1. **原稿**: `260203_presentation_manuscript_for_collaborators.md`
2. **図表**: 本ディレクトリ内の68個のPDFファイル

---

*作成: 2026-02-03*
*更新: PDFのみに整理*
