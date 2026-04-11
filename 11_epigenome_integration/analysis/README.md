# 11_epigenome_integration/analysis - 解析結果ディレクトリ

*Streptomyces coelicolor* A3(2) M145 エピゲノム-トランスクリプトーム統合解析

## ディレクトリ構造

```
analysis/
├── 01_integration/          # 基本統合データ（メチル化-発現統合CSV、相関解析結果）
├── 02_publication_figures/  # 論文用Figure（Fig1-Fig5）
├── 03_overlap_analysis/     # DEG-DMG重複検定（ベン図、分割表）
├── 04_multiomics_tracks/    # ゲノムブラウザ風トラック（RamR, NsdB等）
├── 05_temporal_dynamics/    # 経時変化プロット（T1→T2→T3）
├── 06_motif_schematics/     # プロモーター構造模式図
├── 07_motif_analysis/       # メチル化モチーフ解析（MEME）
├── 08_coordinated_enrichment/ # 協調的変化遺伝子のCOGエンリッチメント
├── 09_spurious_validation/  # 擬似相関検証（Permutation, 偏相関等）
├── 10_atcc_comparative_methylome/ # ATCC 37株比較メチローム解析（旧・シミュレーション）
├── 20_atcc_real_methylome/       # REBASE実データによる比較メチローム解析
│
├── figures -> 02_publication_figures  # 後方互換性シンボリックリンク
├── overlap_analysis -> 03_overlap_analysis
├── multiomics_tracks -> 04_multiomics_tracks
├── temporal_dynamics -> 05_temporal_dynamics
├── motif_schematics -> 06_motif_schematics
├── motif_analysis -> 07_motif_analysis
└── coordinated_enrichment -> 08_coordinated_enrichment
```

## 各ディレクトリの内容

### 01_integration/
メチル化-発現統合の基本データ
- `integrated_methyl_expression_weighted.csv` - 重み付け統合データ（主解析用）
- `correlation_analysis_weighted.csv` - Spearman相関解析結果
- `high_confidence_sites_weighted.csv` - 高信頼度メチル化サイト
- `T2vsT1_coordinated_genes.csv` - 協調的変化遺伝子リスト
- `BGC_detailed.csv`, `TF_coordinated_changes.csv` - BGC/TF解析結果

### 02_publication_figures/
論文用の主要Figure
- `Fig1_methylation_expression_correlation.pdf` - メチル化-発現散布図
- `Fig2_coordination_patterns.pdf` - 協調パターン分布
- `Fig3_volcano_methylation.pdf` - Volcanoプロット
- `Fig4_summary_statistics.pdf` - 統計サマリー
- `Fig5_top_genes_heatmap.pdf` - トップ遺伝子ヒートマップ

### 03_overlap_analysis/
DEG-DMG重複のFisher検定
- `venn_DEG_DMG_*.pdf` - ベン図（T2vsT1, T3vsT1, T3vsT2）
- `contingency_tables.pdf` - 2×2分割表
- `DEG_DMG_overlap_statistics.csv` - 統計結果

### 04_multiomics_tracks/
重要遺伝子座のゲノムブラウザ風可視化
- `multiomics_track_RamR.pdf` - SCO6685（6mA脱メチル化→活性化）
- `multiomics_track_NsdB.pdf` - SCO7252（4mCメチル化獲得→活性化）
- `multiomics_track_Act_*.pdf`, `*_Red_*.pdf`, `*_Cpk_*.pdf` - BGCトラック

### 05_temporal_dynamics/
T1→T2→T3の経時変化
- `temporal_dynamics_*.pdf` - 個別遺伝子の経時プロット
- `methylation_expression_trajectories.pdf` - 2D軌跡プロット

### 06_motif_schematics/
プロモーター構造と-35/-10 box関係
- `promoter_architecture_summary.pdf` - プロモーター構造サマリー
- `motif_schematic_*.pdf` - 個別遺伝子の模式図

### 07_motif_analysis/
MEME de novoモチーフ発見
- `meme_6mA/`, `meme_4mC/` - MEME出力
- `*_base_composition.png` - 塩基組成プロット
- `*_known_motifs.csv` - 既知R-Mモチーフとのマッチング

### 08_coordinated_enrichment/
協調的変化遺伝子の機能エンリッチメント
- `COG_comparison_pos_vs_neg.png` - 正vs負相関のCOG比較
- `coordination_summary.png` - 相関タイプ分布

### 09_spurious_validation/
擬似相関検証（完了）
- Permutation test（並べ替え検定）
- 偏相関分析（交絡因子制御）
- Bootstrap信頼区間
- 負の対照群との比較
- **結果**: 4mC T2vsT1相関は頑健（p=0.0006）

### 10_atcc_comparative_methylome/（旧・シミュレーション）
ATCC Genome Portal 37株比較メチローム解析（**シミュレーションデータ**、`20_atcc_real_methylome/` で置換）
- `comparative_methylome_data.csv` - 38株（M145含む）の比較データ（乱数生成）
- **注意**: このディレクトリの保存率は `np.random.seed(42)` によるシミュレーション値

### 20_atcc_real_methylome/（実データ）
REBASE v602 + NCBI RefSeqによる実データ比較メチローム解析
- `motif_conservation_rebase.pdf/svg` - REBASE R-M系保存率バーチャート
- `rm_system_heatmap.pdf/svg` - 37株×モチーフ ヒートマップ
- `motif_site_density_boxplot.pdf/svg` - モチーフサイト密度 M145 vs 属
- `simulated_vs_real_comparison.pdf/svg` - シミュレーション vs 実データ比較
- `rm_type_distribution.pdf/svg` - R-M系タイプ分布
- **主要発見**:
  - CCGG (4mC): **25.0%** (9/36株) — シミュレーション86.5%から大幅低下
  - AAGCCCG (6mA): **0.0%** (0/36株) — **属レベルで完全に新規**
  - GATC (6mA): **16.7%** (6/36株) — シミュレーション67.6%から大幅低下

### 11_rm_system_identification/ ~ 16_t3_coordinated/
R-M系同定、GRN-TFメチル化統合、SC_RS17645配列解析、AAGCCCGプロモーター解析、Act BGCエピジェネティック制御、T3協調変動

### 17_paper_figures/ ~ 19_bgc_regulator_overview/
論文Figure生成、TSS解析、BGCレギュレーター概観

### 21_genuswide_motif_conservation/ ~ 24_5mC_vs_4mC_CCGG/
属レベルモチーフ保存性、4mC/5mC二重修飾、拡張モチーフ探索、CCGG修飾種判定

### 25_4mC_6mA_differential/ ~ 26_redZ_paradox/
**Loop 1**: 4mC vs 6mA差異的制御（H1）、redZパラドックス解明（H3）

### 27_TF_methylation_rescreen/ ~ 34_gatekeeper_model_v2/
**Loop 2-3**: TF再スクリーニング（H4）、MTase安定性（H5）、ゲノムワイドTFスクリーン（H6）、CCGG MTaseパラドックス（H7）、協調変動レギュレーター特性（H8）、CCGG 5mC誤分類検証（H9）、GCCGGC R-M同定（H10）、Gatekeeperモデルv2（H11）

### 35_GCCGGC_MTase_BLAST/ ~ 42_GCCGGC_temporal_derepression/
**Loop 4-6**: GCCGGC MTase BLAST（H12）、AAGCCCG分布（H13）、Defense Island（H14）、Cross-motif回避（H15）、GCCGGC MTase逆同定（H16）、モチーフ分業（H17）、用量反応（H18）、時間的脱抑制（H19）

### 43_regulatory_avoidance_geographic_test/ ~ 52_shielded_exposed_boundary/
**Loop 7-10**: 地理的交絡検証（H20）、AAGCCCG因果性（H21）、配列レベル枯渇（H22）、カテゴリ特異性（H23）、BGCメチル化（H24）、TSSメチル化勾配（H25）、TFBS保護（H26）、協調レギュレーター保護（H27）、Exposed特性（H28）、Shielded/Exposed境界（H29）

### 53_TSS_sequence_determinants/ ~ 60_coexpression_regulon_prediction/
**Loop 11-14**: TSS配列決定因子（H30）、Exposed TF下流ネットワーク（H31）、自己制御モジュール（H32）、近傍転写影響（H33）、時間的ダイナミクス（H34）、AAGCCCG因果経路（H35-causal）、TFファミリー機能予測（H35-TF/58b）、進化的保存性（H36）、共発現レギュロン（H37）

### archive/
旧バージョンのデータファイル（v1_weighted_minreps2等）

## シンボリックリンクについて

後方互換性のため、古いディレクトリ名へのシンボリックリンクを維持しています。
既存のスクリプトやレポートは引き続き動作します。

**例**: `figures/` → `02_publication_figures/`

## 最終更新

- 2026-03-05: README更新（21-60番台の記述追加、58番重複修正→58b）
- 2026-02-04: REBASE実データ比較メチローム解析追加（20_atcc_real_methylome/）
- 2026-02-03: ディレクトリ構造を番号付きに整理、シンボリックリンク追加
- 2026-02-02: DEG-DMG重複検定、マルチオミクストラック、経時変化プロット追加
