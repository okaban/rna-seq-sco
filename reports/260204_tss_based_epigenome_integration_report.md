# TSS-based Epigenome-Transcriptome Integration Report

**日付**: 2026-02-04
**プロジェクト**: M145_RNA-seq (*Streptomyces coelicolor* A3(2) M145)
**解析ディレクトリ**: `11_epigenome_integration/analysis/18_tss_analyses/`
**スクリプト**: `11_epigenome_integration/scripts/tss_analyses.py`

---

## Key Insights（主要なインサイト一覧）

| # | Insight | 対応Figure | 根拠となる数値 |
|---|---------|-----------|---------------|
| 1 | GFF TSSは実験的TSSより中央値45 bp上流にずれており、プロモーター解析の精度に影響する | `A1_tss_offset_histogram.png` | Median offset = -45 bp, 52.4%が±50 bp以内 |
| 2 | TSS直近に明瞭なメチル化枯渇域（"TSS dip"）が存在し、6mA・4mCともに確認される | `B1_metagene_methylation_profile.png`, `B1_metagene_promoter_zoom.png` | TSS±50 bp領域で密度が周辺の約50%に低下 |
| 3 | 実験的TSSを使用するとメタジーンプロファイルのTSS dipがより鮮明になり、TSSデータの品質が解析解像度を直接左右する | `B1_metagene_by_tss_source.png` | Exp. TSS (n=2,703) vs GFF (n=5,380) |
| 4 | T3で4mCがゲノムワイドに激減し、メチル化の時間的ダイナミクスが顕著 | `B1_metagene_methylation_profile.png` | T3の4mC密度: T1/T2の約1/3 |
| 5 | TSS直下流（+0~+50 bp）のメチル化のみが発現と正の相関を示し、位置依存的な制御を示唆 | `C1_distance_correlation_heatmap.png` | 6mA TSS_proximal: r=+0.29 (T2vsT1); 4mC: r=+0.25* (p=0.041) |
| 6 | 4mC proximal_upstream（-100~-50 bp）はT3vsT1で強い負の相関を示す | `C1_distance_correlation_heatmap.png` | r=-0.40* (p<0.05, T3_vs_T1) |
| 7 | 6mAの発現への影響は後期(T3)で顕在化し、プロモーター・遺伝子本体ともに負の相関 | `C2_promoter_vs_genebody_scatter.png` | Promoter: r=-0.128* (p=0.015); Gene body: r=-0.110*** (p=7.4e-4) |
| 8 | ~~-10 box領域のメチル化枯渇は配列組成（低モチーフ密度）で説明される~~ **[REVISED]** | `null_model_minus10_box.png` | モチーフ帰無モデル p=0.43; -10 boxのCCGG密度=14.6/kb vs control 32.2/kb |
| 9 | AAGCCCG・CCGGモチーフはTSS上流~155 bpに集中し、σ因子結合部位とは空間的に分離 | `motif_tss_position_distribution.png`, `motif_promoter_architecture_map.png` | AAGCCCG中央値: -155 bp, -35 boxとの距離中央値: 106 bp |
| 10 | AAGCCCGはσ70 boxと物理的重複ゼロ。-10 box枯渇の主因はターゲットモチーフの配列非親和性 | `motif_sigma_distance_histogram.png`, `motif_core_promoter_zoom.png` | AAGCCCG-σ overlap: 0件, CCGG-σ overlap: <0.1% |
| 11 | ~~SARP zone BGC保護~~ → ランダム期待値と差なし（Fisher p=0.63） **[REVISED]** | `null_model_sarp_zone.png` | 期待値0.5 vs 観測0; ゲノム全体でもSARP zoneメチル化率はわずか1.6% |
| 12 | メチル化モチーフ組成はタイムポイント間で安定だが、Gained/Lostサイトで有意に異なる | `temporal_key_motifs_barplot.png`, `temporal_gained_lost_motif_barplot.png` | 4mC-CCGG: T1 79.0%→T3 80.4%; Gained vs Lost差 p<1e-4 |
| 13 | T2vsT1とT3vsT2でGained/Lost間のモチーフ濃縮パターンが逆転する | `temporal_gained_lost_motif_barplot.png` | T2vsT1: CCGG→gained; T3vsT2: CCGG→lost, AAGCCCG→gained |
| 14 | AAGCCCG上の4mC獲得は発現低下と連動し、CCGG上の4mC消失は発現上昇と連動する | `temporal_expression_motif_heatmap.png` | 4mC gained×DOWN: AAGCCCG 65-100%; 4mC lost×UP: CCGG 81-84% |
| 15 | TCGA（SARPコアモチーフ）は後期(T3vsT2)の6mA新規獲得サイトに濃縮される | `temporal_gained_lost_motif_barplot.png` | T3vsT2 6mA gained: TCGA 8.1% vs lost 5.0% (p=0.021) |
| 16 | **Dcm-like MTaseがT3で10-40倍発現上昇するにもかかわらずm4C-CCGGが56%崩壊→m4C≠m5C二重修飾系仮説** | `mtase_methylation_dynamics.png` | SC_RS36410 T3vsT2 log2FC=+3.76*; 4mC: 2458→1077 |
| 17 | SC_RS17645の発現動態がAAGCCCGメチル化のGained/Lostパターンを予測する | `mtase_methylation_dynamics.png` | T2vsT1: -2.19* → AAGCCCG lost; T3vsT2: +1.45* → AAGCCCG gained |
| 18 | 転写因子プロモーターのメチル化はゲノム平均と差なし（Fangモデル不支持） | `tf_methylation_enrichment.png` | TF: 15.0% vs non-TF: 17.1% (Fisher p=0.94); BGCレギュレーター8遺伝子全てメチル化なし |
| 19 | FDR補正後、相関解析で唯一生残するのはC2:6mA:T3vsT1:gene_body (r=-0.110, FDR p=0.049) | `all_correlations_fdr.csv` | 66検定中1件のみFDR<0.05; Fisher検定は3/36件生残 |

---

## 1. 解析の背景と目的

本解析では、dRNA-seqで実験的に同定されたTSS（Jeong et al. 2016, Nature Commun.）を導入し、GFFアノテーションベースのTSSを精緻化した上で、TSS基準のエピゲノム解析パイプラインを構築した。

**目的:**
- (A1) 実験的TSSの取り込みとGFF TSSとの比較評価
- (B1) TSS基準のメタジーンメチル化プロファイル（論文Figure品質）
- (C1) TSS距離帯別のメチル化-発現相関マッピング
- (C2) プロモーター vs 遺伝子本体メチル化の効果分離
- (D1) σ70 (HrdB) 結合部位とメチル化の空間的関係

---

## 2. A1: 実験的TSSの取り込み

### 2.1 データソース

- **Jeong et al. 2016** (Nature Commun. 7:11605): S. coelicolorのdRNA-seq解析
- Supplementary Data 1から3,570 TSS（Primary: 2,771, Secondary: 333, Antisense: 256, Intergenic: 131, Internal: 79）
- SCO遺伝子IDをGFFのold_locus_tag経由でSC_RSにマッピング

### 2.2 主要結果

| 指標 | 値 |
|------|-----|
| Jeong TSS総数 | 3,570 |
| Primary TSS | 2,771 |
| SC_RSにマッピング成功 | 2,703 (97.5%) |
| 全GFF遺伝子数 | 8,083 |
| 実験的TSSカバー率 | 33.4% |
| GFF-onlyの遺伝子 | 5,380 (66.6%) |

### 2.3 TSS Offset解析

> **Insight #1**: GFFアノテーションのTSS（= ORF開始位置）は、実験的TSS（= 転写開始点）に比べて中央値45 bp上流にある。
> これはGFFがATG（翻訳開始コドン）を遺伝子開始としてアノテートしているため、5' UTR長（~45 bp）がオフセットとして現れたもの。
> **Figure**: `A1_tss_offset_histogram.png`

| パラメータ | 値 |
|-----------|-----|
| Offset中央値 | -45 bp |
| Offset平均値 | -72.6 bp |
| 標準偏差 | 91.1 bp |
| ±50 bp以内 | 1,416 (52.4%) |
| ±100 bp以内 | 1,978 (73.2%) |

**解釈**: S. coelicolorの5' UTRは中央値~45 bpと比較的短く、バクテリアの典型的範囲（20-100 bp）に収まる。しかし27%以上の遺伝子で100 bp以上のオフセットがあり、GFF-only TSSに基づくプロモーター解析には系統的バイアスが含まれることを示す。

---

## 3. B1: メタジーンメチル化プロファイル

### 3.1 全遺伝子プロファイル（±1000 bp）

> **Insight #2**: TSS直近に明瞭なメチル化枯渇域（"TSS dip"）が存在する。
> 6mA、4mCともにTSS周辺（約-100~+100 bp）で密度が低下し、特に4mCの枯渇が顕著。
> これは転写開始部位における配列選択圧またはクロマチンアクセシビリティの反映と考えられる。
> **Figure**: `B1_metagene_methylation_profile.png`, `B1_metagene_promoter_zoom.png`

**観察された主要パターン:**

1. **TSS dip**: 6mA・4mCともにTSS±50 bp付近で密度が最小値を示す
2. **4mCのT3での劇的減少**: T3の4mC密度はT1/T2の約1/3に低下（Insight #4）
3. **6mAの上流ピーク**: TSS上流約-200~-300 bpに6mA密度のピークが存在
4. **非対称パターン**: メチル化は上流（プロモーター側）で下流（遺伝子本体側）より高密度

### 3.2 TSS source別比較

> **Insight #3**: 実験的TSS使用時にTSS dipがより鮮明になる。
> GFF-only遺伝子では45 bpのオフセットにより信号がblurされ、解析解像度が低下する。
> **Figure**: `B1_metagene_by_tss_source.png`

**重要な含意**: 今後のプロモーター解析では、実験的TSSが利用可能な2,703遺伝子を優先し、GFF-only遺伝子にはオフセット補正を検討すべき。

---

## 4. C1: TSS距離帯別メチル化-発現相関

### 4.1 距離帯定義

| 距離帯 | TSS相対位置 | 機能的意味 |
|--------|-----------|-----------|
| distal_upstream | -300 ~ -200 bp | 遠位プロモーター |
| mid_upstream | -200 ~ -100 bp | 中間プロモーター |
| proximal_upstream | -100 ~ -50 bp | 近位プロモーター |
| core_promoter | -50 ~ 0 bp | コアプロモーター |
| TSS_proximal | 0 ~ +50 bp | TSS直下流/5' UTR |
| early_gene_body | +50 ~ +200 bp | 遺伝子本体前半 |
| mid_gene_body | +200 ~ +500 bp | 遺伝子本体中間 |
| distal_gene_body | +500 ~ +1000 bp | 遺伝子本体後半 |

### 4.2 主要結果

> **Insight #5**: TSS直下流（0~+50 bp）のメチル化変動のみが発現変動と正の相関を示す。
> 6mA: r=+0.29 (p=0.084, T2 vs T1), 4mC: r=+0.25 (p=0.041, T2 vs T1)
> これはプロモーター上流のメチル化が抑制的であるのに対し、TSS直下流のメチル化は転写活性と連動する可能性を示唆する。
> **Figure**: `C1_distance_correlation_heatmap.png`

> **Insight #6**: 4mCのproximal_upstream (-100~-50 bp) はT3 vs T1で強い負の相関（r=-0.40）を示す。
> 後期タイムポイントでの4mC上流メチル化の増加が発現抑制と連動することを示し、4mCの位置特異的な抑制機能を示唆。
> **Figure**: `C1_distance_correlation_heatmap.png`

**距離依存パターンのまとめ（T2 vs T1）:**

| 領域 | 6mA | 4mC | 解釈 |
|------|-----|-----|------|
| Upstream (-300 ~ -50) | 弱い負/正 (ns) | 弱い負 (ns) | シグナル不明瞭 |
| Core promoter (-50 ~ 0) | +0.13 (ns) | +0.01 (ns) | 中立 |
| **TSS proximal (0 ~ +50)** | **+0.29 (trend)** | **+0.25*** | **正の相関** |
| Gene body (+50 ~ +1000) | 弱い負 (ns) | 弱い正/負 (ns) | 効果小 |

---

## 5. C2: プロモーター vs 遺伝子本体メチル化

### 5.1 主要結果

> **Insight #7**: 6mAの発現抑制効果は後期（T3）で顕在化する。
> T2 vs T1では6mAプロモーターメチル化は非有意（r=0.048, ns）であるが、
> T3 vs T1で有意な負の相関に転じる（promoter: r=-0.128, p=0.015; gene body: r=-0.110, p=7.4e-4）。
> これは培養後期での6mAの蓄積とそれに伴う転写抑制を示唆する。
> **Figure**: `C2_promoter_vs_genebody_scatter.png`, `C2_promoter_vs_genebody_barplot.png`

**比較結果一覧（全comparisons）:**

| 修飾 | 領域 | T2 vs T1 | T3 vs T1 | T3 vs T2 |
|------|------|----------|----------|----------|
| 6mA | Promoter | r=+0.048 (ns) | r=-0.128* | r=-0.121* |
| 6mA | Gene body | r=-0.012 (ns) | r=-0.110*** | r=-0.007 (ns) |
| 4mC | Promoter | r=+0.001 (ns) | r=-0.044 (ns) | r=+0.013 (ns) |
| 4mC | Gene body | r=-0.011 (ns) | r=-0.038 (ns) | r=+0.025 (ns) |

**解釈:**
- 6mAはプロモーター・遺伝子本体ともに後期で抑制的に作用
- 4mCはいずれの領域でも有意な全体相関を示さない（ただしC1で距離帯別に見ると位置特異的シグナルあり）
- プロモーター vs 遺伝子本体の効果サイズ差は小さく、6mAの抑制効果は領域非特異的

---

## 6. D1: σ因子結合部位とメチル化の重なり

### 6.1 σ70 (HrdB) モチーフ探索結果

| 指標 | -35 box | -10 box |
|------|---------|---------|
| 検出モチーフ数 | 684 | 1,991 |
| モチーフ保有遺伝子数 | 656 | 1,948 |
| 両方を保有 | 183 |  |
| コンセンサス | [TC]TGAC[N] | TA[N]{3}T |

### 6.2 メチル化密度の位置特異性

> **Insight #8**: -10 box領域はメチル化が有意に枯渇している（p=1.73e-12）。
> 当初の仮説「メチル化がσ因子結合を阻害する」とは逆に、-10 box領域にメチル化が進化的に排除されていることが判明。
> これはσ因子結合に必須の配列にメチル化ターゲットモチーフが存在しにくいという配列レベルの制約を反映する可能性がある。
> **Figure**: `D1_promoter_methylation_distribution.png`

| 領域 | TSS相対位置 | メチル化サイト数 | 密度 (sites/gene/kb) | χ² p-value vs control |
|------|-----------|----------------|---------------------|----------------------|
| -35 region | -45 ~ -25 bp | 199 | 1.17 | 0.149 (ns) |
| **-10 region** | **-20 ~ 0 bp** | **107** | **0.63** | **1.73e-12 (***)** |
| Control | -100 ~ -60 bp | 441 | 1.33 | — |

### 6.3 直接的なモチーフ-メチル化重複

- メチル化-σモチーフ重複: **21件のみ**（うちwithin: 15, flanking: 6）
- 大部分が-10 box（19/21）
- メチル化がσ結合部位を直接修飾するケースは極めてまれ

### 6.4 仮説の更新

当初仮説: 「メチル化がσ因子結合部位に集積し、転写開始を阻害する」
→ **改訂仮説**: 「-10 box領域はメチル化を回避するよう進化的に制約されており、σ因子の結合は配列レベルでprotectされている。メチル化による転写制御は、σ結合部位以外のプロモーター領域（-100~-50 bp付近）で間接的に作用する可能性がある」

---

## Figure解説

### A1_tss_offset_histogram.png

**ファイル**: `analysis/18_tss_analyses/A1_tss_offset_histogram.png`

**説明**: GFFアノテーションベースTSSと実験的TSS（Jeong et al. 2016 dRNA-seq）のオフセット分布

**読み取り方**:
- X軸: TSS offset（実験的TSS − GFF TSS, bp, ストランド補正済み）
- Y軸: 遺伝子数
- 赤破線: 完全一致（offset = 0）
- 橙破線: 中央値（-45 bp）

**注目点**:
- 分布は0付近にピークがあるが、負の方向に長い裾を持つ
- 中央値が-45 bpであることは、GFF TSSが実際のTSSの上流（5' UTR + ORF開始位置）を指すことを反映
- ~300 bp以上離れたケースも存在し、一部はリーダーレス転写やmisannotationの可能性

**サポートするInsight**: #1

### B1_metagene_methylation_profile.png

**ファイル**: `analysis/18_tss_analyses/B1_metagene_methylation_profile.png`

**説明**: 全8,083遺伝子のTSSを基準としたメタジーンメチル化密度プロファイル（±1000 bp）

**読み取り方**:
- X軸: TSSからの距離（bp）。負=上流（プロモーター）、正=下流（遺伝子本体）
- Y軸: メチル化密度（sites / gene / kb）
- 上段: 6mA、下段: 4mC
- 青=T1(early), 赤=T2(mid), 緑=T3(late)
- 紫/青の帯: -10/-35 box付近

**注目点**:
- TSS±100 bp付近の明瞭な枯渇域（"TSS dip"）
- 4mCのT3での劇的減少（緑線が全領域で低い）
- 6mAは上流-200~-300 bpにピーク

**サポートするInsight**: #2, #4

### B1_metagene_promoter_zoom.png

**ファイル**: `analysis/18_tss_analyses/B1_metagene_promoter_zoom.png`

**説明**: プロモーター領域のズーム表示（-500 ~ +200 bp）。-35 box, -10 box, TSSの位置を明示

**読み取り方**:
- 紫帯: -35 box (-45~-30 bp)
- 橙帯: -10 box (-15~-5 bp)
- 緑線: TSS

**注目点**:
- -10 box〜TSS付近で密度が最低値を示す
- -35 box領域では密度がやや回復

**サポートするInsight**: #2, #8

### B1_metagene_by_tss_source.png

**ファイル**: `analysis/18_tss_analyses/B1_metagene_by_tss_source.png`

**説明**: 実験的TSS (n=2,703) vs GFF-only TSS (n=5,380) でのメタジーン比較

**読み取り方**:
- 青実線: 実験的TSS遺伝子、橙破線: GFF-only遺伝子
- タイムポイントは全期間を合算

**注目点**:
- 実験的TSS群でTSS dipがより鮮明
- GFF-only群ではTSS dipがblurされて浅い

**サポートするInsight**: #3

### C1_distance_correlation_heatmap.png

**ファイル**: `analysis/18_tss_analyses/C1_distance_correlation_heatmap.png`

**説明**: 8つのTSS距離帯 × 3比較 × 2修飾タイプのSpearman相関ヒートマップ

**読み取り方**:
- 行: 距離帯（上流→下流の順）
- 列: タイムポイント比較（T2vsT1, T3vsT1, T3vsT2）
- 色: Spearman r（赤=正の相関、青=負の相関）
- アスタリスク: *p<0.05, **p<0.01, ***p<0.001

**注目点**:
- TSS_proximal行のT2vsT1が赤い（正の相関）
- 4mC proximal_upstreamのT3vsT1が濃い青（強い負の相関 r=-0.40）
- 6mAは後期で全般的に負の相関が強まる

**サポートするInsight**: #5, #6

### D1_promoter_methylation_distribution.png

**ファイル**: `analysis/18_tss_analyses/D1_promoter_methylation_distribution.png`

**説明**: コアプロモーター領域（-100 ~ +10 bp）のメチル化サイト分布ヒストグラム

**読み取り方**:
- X軸: TSSからの距離（bp）
- Y軸: メチル化サイト数（全タイムポイント合計）
- 紫帯: -35 box領域
- 橙帯: -10 box領域
- 緑線: TSS

**注目点**:
- -10 box〜TSS領域で明らかにサイト数が減少
- -35 box領域では比較的高いメチル化が維持
- -60 bp付近に小さなピーク

**サポートするInsight**: #8

---

## 6b. 追加解析: メチル化モチーフ-TSS-σ因子結合部位の空間的関係

### 6b.1 モチーフのTSS相対位置分布

| モチーフ | 出現数 | 保有遺伝子数 | TSS相対位置中央値 | 機能的位置 |
|---------|--------|-------------|----------------|----------|
| AAGCCCG (6mA) | 849 | 785 | **-155 bp** | 遠位～中間プロモーター領域 |
| CCGG (4mC) | 96,643 | 8,082 | **-156 bp** | 全域（GCリッチゲノムで遍在） |

> **Insight #9**: AAGCCCG・CCGGともにTSSの上流約150 bpに中央値を持ち、-35/-10 boxとは空間的に分離している。
> AAGCCCGは-35 boxとの距離中央値が106 bp、-10 boxとは139 bp離れており、σ因子結合部位への直接干渉は生じにくい。
> **Figure**: `motif_tss_position_distribution.png`, `motif_promoter_architecture_map.png`

### 6b.2 σ因子結合部位との距離

| モチーフ | 対象box | 距離中央値 | ±5 bp以内 | ±10 bp以内 | 重複（モチーフ長以内） |
|---------|--------|----------|----------|-----------|-------------------|
| AAGCCCG | -35 box | **106 bp** | 0 (0.0%) | 1 (1.4%) | 0 |
| AAGCCCG | -10 box | **139 bp** | 0 (0.0%) | 1 (0.4%) | 0 |
| CCGG | -35 box | **125 bp** | 67 (0.9%) | 196 (2.6%) | 57 |
| CCGG | -10 box | **150 bp** | 169 (0.7%) | 456 (2.0%) | 108 |

> **Insight #10**: AAGCCCGモチーフはσ70の-35/-10 boxと物理的に重複するケースが**ゼロ**。
> CCGGは遍在性ゆえ少数の重複（57件/-35, 108件/-10）が見られるが、全体の0.1%未満。
> メチル化モチーフとσ因子結合部位は配列レベルで「棲み分け」ており、D1の-10 box枯渇の主因は**ターゲットモチーフの配列非親和性**であると結論。
> **Figure**: `motif_sigma_distance_histogram.png`, `motif_core_promoter_zoom.png`

### 6b.3 メチル化状態と位置の関係

| モチーフ | メチル化あり | メチル化なし | 位置差(中央値) | Mann-Whitney U p |
|---------|-----------|-----------|-------------|-----------------|
| AAGCCCG | 669 (78.8%) | 180 (21.2%) | -162 vs -142 bp | 0.657 (ns) |
| CCGG | 1,191 (1.2%) | 95,452 (98.8%) | -181 vs -155 bp | 0.417 (ns) |

AAGCCCGの78.8%が実際にメチル化されている一方、CCGGは1.2%のみ。CCGGは配列モチーフとしてはコアプロモーターにも存在するが、-10 box付近のCCGGは選択的にメチル化を回避している可能性がある。

### 6b.4 コアプロモーター領域のズーム観察

ズーム図（`motif_core_promoter_zoom.png`）から:
- **AAGCCCG**: コアプロモーター(-80~+30 bp)に119件のみ。-35 box領域にも散発的に出現するが密度は低い。-10 box〜TSS付近は**ほぼ空白**
- **CCGG**: -10 box領域(-20~0 bp)で明らかに密度が低下（周辺の約40%に減少）。TSS直後(+25~+30 bp)にピーク

### 6b.5 統合的解釈: プロモーターの「メチル化フリーゾーン」

```
Position:  -500    -300    -200    -100    -50   -35  -10  TSS  +50
              │       │       │       │      │     │    │    │    │
AAGCCCG:     ████████████████████████│      │     │    │    │    │
             [高密度: -500〜-100]     │ [低密度]  [空白]   [空白]
                                     │      │     │    │    │    │
CCGG:        ████████████████████████│██████│████ │[枯渇]│    │████
             [高密度: 遍在]          │      │     │    │    │    │
                                     │      │     │    │    │    │
σ70 boxes:                           │      │ [-35]│[-10]│ TSS │
                                     │      │     │    │    │    │
METHYL-FREE ZONE:                    │      │     ├────┤    │
                                     │      │     │ ←20bp→  │
```

-10 box〜TSS付近に「メチル化フリーゾーン」が存在し、これは以下の2つのメカニズムによる:
1. **配列レベル**: メチル化ターゲットモチーフ（特にAAGCCCG）がこの領域に配列的に存在しにくい
2. **選択的メチル化回避**: CCGG（遍在モチーフ）がこの領域に存在してもメチル化されにくい

---

## 6c. 追加解析: SARP結合領域とメチル化の空間的統合

### 背景

放線菌のSARP（Streptomyces Antibiotic Regulatory Protein）は、-35 box付近（-60~-20 bp）にダイレクトリピート（TCGAGC[G/C], 11 bp周期）で結合し、σHrdB R4が-35 boxに直接接触する代わりにSARPプロトマー上に載ることで転写活性化を行う（Krysenko et al. 2024, Nature Commun.）。従って、**SARP結合領域のメチル化はσ因子ではなくSARP結合を阻害し得る**。

### 6c.1 SARP zone (-60~-20 bp) のメチル化

| 指標 | 値 |
|------|-----|
| SARP zone内メチル化サイト数 | 344 (6mA: 238, 4mC: 106) |
| メチル化を持つ遺伝子数 | 158 |
| うちBGC遺伝子 | 2 (cdaV: 4mC, cdaK: 6mA) |
| AAGCCCG motif in SARP zone | 42 genes (BGC: 0) |
| CCGG motif in SARP zone | 3,957 genes (BGC: 48) |

### 6c.2 SARP motif – メチル化モチーフ共存

| 指標 | 遺伝子数 |
|------|---------|
| SARP motif (TCGAGC等) in zone | 1,401 |
| メチル化ターゲットモチーフ in zone | 3,974 |
| 実際のメチル化サイト in zone | 158 |
| SARP motif ∩ メチル化モチーフ | 583 |
| SARP motif ∩ 実際のメチル化 | **15** |
| Triple overlap (3者共存) | **13** |

> **Insight #11**: SARP結合モチーフと実際のメチル化サイトが同一遺伝子のSARP zoneで共存するケースは**わずか15遺伝子**。配列レベルではSARP motifとメチル化ターゲットが583遺伝子で共存するが、実際にメチル化されるのはその2.6%にすぎない。SARP結合領域はメチル化から機能的に保護されている可能性がある。
> **Figure**: `sarp_methylation_mechanism_schematic.png`

### 6c.3 SARP zone メチル化-発現相関（ゲノムワイド）

| 修飾 | 比較 | Spearman r | p-value | n |
|------|------|-----------|---------|---|
| 6mA | T2 vs T1 | +0.11 | 0.46 (ns) | 47 |
| 6mA | T3 vs T1 | -0.06 | 0.68 (ns) | 47 |
| 4mC | T2 vs T1 | -0.01 | 0.95 (ns) | 31 |

SARP zone単独では有意な相関が検出されない。これは:
1. サンプルサイズの制約（n=31~47）
2. SARP依存的プロモーターとそうでないプロモーターが混在
3. SARP zone内でもメチル化が機能的に影響するサブセットが限定的

### 6c.4 BGC遺伝子の個別解析

SARP zoneにメチル化を持つBGC遺伝子は**2遺伝子のみ**:
- **cdaV** (SC_RS18355): 4mCが-20 bpに全タイムポイントで安定存在（T1: 84.8%, T2: 88.4%, T3: 83.3%）→ CDA BGCのSARP依存的制御に影響する可能性
- **cdaK** (SC_RS18310): 6mAが-59 bpにT2のみ出現（51.2%）→ 一過性メチル化

Act, Red, Cpk BGCのSARP zoneにはメチル化が**見られない**。これは、主要BGCのSARP結合領域がメチル化から保護されていることを意味する。

### 6c.5 統合的解釈

**仮説の検証結果:**

当初仮説: 「SARP zone (-60~-20 bp)のメチル化がSARP結合を阻害し、BGC遺伝子の発現を抑制する」

**結論: 仮説は部分的に支持されるが、主要BGCでは関連が乏しい**

1. **SARP zoneは全般的にメチル化密度が低い** → D1のコアプロモーター枯渇パターンの一部
2. **主要BGC (Act, Red, Cpk)のSARP zoneにメチル化なし** → これらのBGCではメチル化によるSARP阻害は起きていない
3. **CDA BGCのcdaVは例外**: SARP zone内に安定した4mCを持ち、SARP結合への影響が考えられる
4. **ゲノムワイドでは158遺伝子がSARP zoneにメチル化を持つ** → 非BGC遺伝子のSARP依存的制御（AfsR等のグローバルSARP）への影響は今後の検討課題

> **Figure**: `region_correlation_with_sarp_zone.png` — 全9領域の相関を並べて比較。SARP zoneは弱い正の相関（6mA T2vsT1: r=+0.11）を示すが非有意。TSS proximalの正相関（4mC T2vsT1: r=+0.25*）が最も強いシグナル。

---

## 6d. 追加解析: タイムポイント別メチル化モチーフと発現関連モチーフ

### 解析の背景

先行解析ではメチル化モチーフを全タイムポイントプールで評価していた。本解析では以下の4つの問いに回答する:
1. T1/T2/T3ごとのモチーフ含有率は変化するか？
2. Gained/Lostサイト間でモチーフ組成は異なるか？
3. 発現変動方向（UP/DOWN）によってモチーフ組成は変わるか？
4. これらのパターンは比較間（T2vsT1, T3vsT1, T3vsT2）で一貫するか？

### 6d.1 タイムポイント別モチーフ含有率

> **Insight #12**: メチル化モチーフの全体的組成はT1→T2→T3で極めて安定している。
> 主要モチーフの含有率変動は2%以内であり、メチル化酵素のターゲット特異性は培養期間を通じて一定。
> ただしT3で4mCサイト総数が激減（2,458→1,077）するため、"何が残るか"より"何が消えるか"に注目すべき。
> **Figure**: `temporal_key_motifs_barplot.png`, `temporal_motif_heatmap.png`

| モチーフ | 修飾 | T1 (%) | T2 (%) | T3 (%) | 傾向 |
|---------|------|--------|--------|--------|------|
| CCGG | 4mC | 79.0 (n=1,995) | 79.8 (n=2,458) | 80.4 (n=1,077) | 安定 |
| AAGCCCG | 4mC | 35.3 | 34.9 | 34.3 | 安定 |
| CCGG | 6mA | 49.8 (n=1,889) | 49.4 (n=2,102) | 47.7 (n=2,257) | 微減 |
| AAGCCCG | 6mA | 22.9 | 22.1 | 20.9 | 微減 |
| TCGA | 6mA | 5.6 | 5.3 | 6.4 | **T3で微増** |
| GCGC | 4mC | 16.3 | 15.5 | 13.6 | **T3で漸減** |

### 6d.2 Gained/Lostサイト間のモチーフ差異

> **Insight #13**: Gainedサイト（新規獲得）とLostサイト（消失）ではモチーフ組成が有意に異なり、その差異パターンはT2vsT1とT3vsT2で逆転する。
> これは培養フェーズの移行に伴い、メチル化の「ターゲット選択性」が変化することを示唆する。
> **Figure**: `temporal_gained_lost_motif_barplot.png`

**T2 vs T1（指数増殖期→遷移期）:**

| 修飾 | Gainedに濃縮 | Lostに濃縮 |
|------|-------------|-----------|
| 4mC | **CCGG** 77.0% vs 62.2% (p=5.47e-5) | **AAGCCCG** 51.2% vs 38.7% (p=1.89e-3) |
| 6mA | 有意差なし | 有意差なし |

**T3 vs T2（遷移期→定常期）:**

| 修飾 | Gainedに濃縮 | Lostに濃縮 |
|------|-------------|-----------|
| 4mC | **AAGCCCG** 72.7% vs 36.3% (p=3.22e-5) | **CCGG** 78.9% vs 57.6% (p=8.35e-3) |
| 6mA | **TCGA** 8.1% vs 5.0% (p=0.021) | — |

**逆転パターンの解釈:**
- **T2vsT1**: CCGGモチーフ部位が新規にメチル化を獲得（4mCの拡大）。AAGCCCGモチーフ部位のメチル化が消失（6mA関連部位の選択的除去）
- **T3vsT2**: 逆にCCGGモチーフ部位が大規模にメチル化を消失（T3での4mC崩壊の主体）。AAGCCCG部位が新規にメチル化を獲得（少数だが有意）
- この逆転は、CCGGメチル化（Dcm型）とAAGCCCGメチル化（独自モチーフ）が異なるメチルトランスフェラーゼに制御されており、培養フェーズごとに活性バランスが変化することを示唆する

### 6d.3 発現関連メチル化モチーフ

> **Insight #14**: AAGCCCGモチーフ上の4mC獲得は発現低下と連動し（gained×DOWN: AAGCCCG 65-100%）、CCGGモチーフ上の4mC消失は発現上昇と連動する（lost×UP: CCGG 81-84%）。
> モチーフの種類によってメチル化の発現への効果方向が異なることを示す。
> **Figure**: `temporal_expression_motif_heatmap.png`

**4mC Gained × Expression カテゴリ別（T2vsT1）:**

| 発現変動 | CCGG (%) | AAGCCCG (%) | n |
|---------|----------|-------------|---|
| UP | 80.4 | 37.0 | 46 |
| DOWN | **55.0** | **65.0** | 20 |
| nonDEG | 74.5 | 40.0 | 55 |

→ 4mCが獲得されて発現が**低下**した遺伝子ではAAGCCCGが65%と突出して高く、CCGGは55%と低い。

**4mC Lost × Expression カテゴリ別（T3vsT2）:**

| 発現変動 | CCGG (%) | AAGCCCG (%) | n |
|---------|----------|-------------|---|
| UP | **83.7** | 30.4 | 92 |
| DOWN | 63.5 | **50.6** | 85 |
| nonDEG | 77.4 | 42.9 | 84 |

→ 4mCが消失して発現が**上昇**した遺伝子ではCCGGが84%と最も高い。一方、消失しても発現が**低下**した遺伝子ではAAGCCCGが51%と高い（他の因子による発現低下を示唆）。

### 6d.4 TCGA（SARPコアモチーフ）の後期濃縮

> **Insight #15**: TCGA（SARP結合モチーフTCGAGCの4bpコア）はT3vsT2の6mA新規獲得サイトに有意に濃縮される（8.1% vs 5.0%, p=0.021）。
> これは後期（定常期）でSARP関連領域に新規6mAメチル化が獲得されることを示し、二次代謝制御のエピジェネティック調節に寄与する可能性がある。
> **Figure**: `temporal_gained_lost_motif_barplot.png`

### 6d.5 統合モデル: フェーズ依存的メチル化ターゲティング

```
T1→T2（増殖→遷移）:
  4mC: CCGG部位にメチル化を獲得 → 遺伝子サイレンシングには限定的
       AAGCCCG部位のメチル化を消失
  6mA: モチーフ選択性に有意差なし

T2→T3（遷移→定常）:
  4mC: CCGG部位から大規模にメチル化消失 → 発現上昇と連動（4mC抑制の解除）
       少数のAAGCCCG部位に新規メチル化 → 発現低下と連動（標的的サイレンシング）
  6mA: TCGA部位（SARP関連）に新規メチル化 → 二次代謝制御への関与を示唆
```

---

## 6e. 統計的検証: 帰無モデル・FDR補正・先行研究統合

### 6e.1 帰無モデルシミュレーション

#### -10 box メチル化枯渇の再評価

> **Insight #8 [REVISED]**: -10 box領域のメチル化枯渇は、**メチル化ターゲットモチーフの配列非親和性（低密度）で完全に説明される**。
> -10 box (TA[N]₃T)のATリッチ配列はGC含量60.7%と低く、CCGG密度は14.6/kbでcontrol領域(32.2/kb)の半分以下。
> モチーフ密度を考慮した帰無モデルではp=0.43と非有意。
> **Figure**: `null_model_minus10_box.png`

| 領域 | GC% | CCGG/kb | AAGCCCG/kb | 観測メチル化 | モチーフ帰無期待値 |
|------|-----|---------|-----------|------------|----------------|
| -10 box | 60.7% | 14.6 | 0.04 | 28 | 29.4 |
| -35 box | 67.0% | 30.1 | 0.30 | 44 | — |
| Control | 69.1% | 32.2 | 0.20 | 132 | — |
| Upstream | 70.8% | 33.4 | 0.14 | 357 | — |

**結論**: -10 boxのメチル化枯渇には進化的選択圧の仮説は不要。ATリッチ配列にGCリッチなメチル化モチーフが少ないという単純な配列組成の帰結である。

#### SARP zone BGC保護の再評価

> **Insight #11 [REVISED]**: SARP zone (40bp幅)のメチル化密度はゲノム全体でわずか1.6%であり、BGCでの不在はランダム期待値と差がない（Fisher p=0.63, 並べ替えp=0.63）。

| 指標 | 値 |
|------|-----|
| ゲノム全体SARP zoneメチル化率 | 1.6% |
| BGC (exp TSS, n=29) メチル化遺伝子 | 0 |
| ランダム期待値 | 0.5 |
| Fisher p (BGC枯渇?) | 0.633 |
| 全遺伝子 (n=8,083) でのBGCメチル化 | 2/99 (2.0%) vs non-BGC 1.9% |

### 6e.2 FDR多重検定補正

> **Insight #19**: 全相関解析66検定にBenjamini-Hochberg補正を適用した結果、FDR<0.05で生残するのは**1件のみ**: C2:6mA:T3vsT1:gene_body (r=-0.110, FDR p=0.049)。

| 解析カテゴリ | 検定数 | FDR<0.05 | 最強シグナル |
|-------------|--------|----------|------------|
| C1距離帯別相関 | 48 | 0 | 4mC:T3vsT2:mid_gene_body (p_adj=0.094) |
| C2プロモーター/遺伝子本体 | 12 | **1** | 6mA:T3vsT1:gene_body (r=-0.110, p_adj=0.049) |
| SARP zone相関 | 6 | 0 | — |
| **Combined** | **66** | **1** | 同上 |

**Gained/Lost Fisher検定（別枠、36検定）:**

| 比較 | モチーフ | Gained% | Lost% | FDR p |
|------|---------|---------|-------|-------|
| T2vsT1 4mC | **CCGG** | 77.0 | 62.2 | **0.0010** |
| T3vsT2 4mC | **AAGCCCG** | 72.7 | 36.3 | **0.0010** |
| T2vsT1 4mC | AAGCCCG | 38.7 | 51.2 | **0.023** |

→ FDR補正後も3件が有意。**モチーフ動態（Gained/Lost間の組成差）は最も堅牢なシグナル**。

### 6e.3 転写因子メチル化濃縮検定（Fangモデル検証）

> **Insight #18**: 転写因子（TF）プロモーターのメチル化率はゲノム平均と有意差なし。Fang et al. (2022)が*S. roseosporus*で示した「m4CがTFプロモーターを介してBGCを間接制御する」モデルは、*S. coelicolor*では支持されない。

| カテゴリ | n | メチル化率 | MW p (TF > non-TF) |
|---------|---|----------|-------------------|
| TF/レギュレーター | 824 | 15.0% | 0.923 |
| 非TF | 7,259 | 17.1% | — |

**BGCレギュレーター個別確認**: actII-orf4, redD, redZ, cdaR, cpkO, absA2, afsR, afsS — **8遺伝子全てプロモーターメチル化なし**。

### 6e.4 MTase発現-メチル化動態の対応

> **Insight #16**: Dcm-like MTase（SC_RS19770, SC_RS36410）はT3で10〜40倍に発現上昇するが、Nanoporeで検出されるm4C-CCGGは56%減少する。この矛盾は、Dcm-like MTaseがm5C（5-メチルシトシン）を生成し、我々のNanoporeが検出するm4C（N4-メチルシトシン）とは異なる修飾であることを示唆する。
> **Figure**: `mtase_methylation_dynamics.png`

| 酵素 | 機能 | T2vsT1 | T3vsT2 | 予測メチル化動態 | 実観測 | 一致? |
|------|------|--------|--------|---------------|--------|------|
| SC_RS17645 | N-6 MTase | **-2.19*** | **+1.45*** | AAGCCCG: T2消失→T3回復 | Gained/Lost逆転 | **一致** |
| SC_RS19770 | Dcm-like | -0.94 (ns) | **+3.54*** | CCGG m5C増加 | m4C崩壊 | **矛盾** |
| SC_RS36410 | Dcm-like | +1.06 (ns) | **+3.76*** | CCGG m5C増加 | m4C崩壊 | **矛盾** |
| SC_RS10665 | SCO1731 m5C | **-1.16*** | -0.23 (ns) | m5C減少 | Pisciotta整合 | **一致** |

> **Insight #17**: SC_RS17645（N-6 DNA methylase, AAGCCCG候補）の発現動態は、AAGCCCGメチル化のGained/Lostパターンと時間的に一致する。T2で酵素発現が-2.19 log2FC低下→AAGCCCG部位が消失、T3で+1.45 log2FC回復→AAGCCCG部位が獲得。
> **Figure**: `mtase_methylation_dynamics.png`

### 6e.5 先行研究との統合

| 先行研究 | 本データとの整合性 | 根拠 |
|---------|:---:|------|
| Pisciotta 2018 (SCO1731 KO) | **整合** | SCO1731のT2/T3での発現低下はm5C減少と一致 |
| Pisciotta 2023 (BS-seq) | **相補的** | BS-seq(m5C)とNanopore(m4C)は異なる修飾。CCGGの二重修飾系を示唆 |
| Fang 2022 (m4C→TF→BGC) | **不支持** | TFプロモーターのメチル化率はゲノム平均と差なし |
| 5-azacytidine研究群 | **整合** | 脱メチル化→Act活性化の方向はT3での4mC崩壊→発現上昇と一致 |
| González-Cerón 2009 | **整合** | SCO3261/3262の劇的発現低下はR-M系の発育段階依存制御と一致 |

### 6e.6 m4C / m5C 二重シトシンメチル化仮説

本解析の最も重要な発見は、以下のエビデンスラインから導かれる**CCGGモチーフ上のm4C/m5C二重修飾系の存在**である:

1. **Dcm-like MTaseのT3発現爆発（+3.5〜+3.8 log2FC）** vs **m4Cの56%崩壊**: 同一酵素が両方を行うなら矛盾
2. **Pisciotta 2023のBS-seq**: m5CをGGCmCGG/GCCmCGモチーフで検出（CCGGを含む）
3. **我々のNanopore**: m4CをCCGGモチーフで検出
4. **BS-seqはm4Cを検出しない**: Bisulfiteはm5Cを保護するがm4Cは保護しない

**提案モデル:**
```
CCGG部位の二重修飾:
  m5C系: Dcm-like MTase (SC_RS19770, SC_RS36410) → T3で活性化
  m4C系: 未同定酵素 → T3で活性低下（→ m4C崩壊）

同一CCGGが時期によって異なる修飾を受ける可能性
  T1-T2: m4C優位（Nanopore検出）
  T3: m5C優位（BS-seq検出相当）→ m4C酵素の不活性化
```

---

## 7. 出力ファイル一覧

| ファイル | 内容 |
|---------|------|
| `comprehensive_tss_table.csv` | 全8,083遺伝子のTSS情報（実験的TSS + GFF TSS統合） |
| `jeong2016_all_tss.csv` | Jeong et al. 2016の全3,570 TSSデータ |
| `A1_tss_offset_histogram.png` | GFF vs 実験的TSSのオフセット分布 |
| `B1_metagene_methylation_profile.png/pdf` | メタジーンプロファイル（±1000 bp） |
| `B1_metagene_promoter_zoom.png/pdf` | プロモーター領域ズーム（-500~+200 bp） |
| `B1_metagene_by_tss_source.png` | TSS source別メタジーン比較 |
| `B1_metagene_data.csv` | メタジーン数値データ |
| `C1_distance_stratified_correlation.csv` | 距離帯別相関データ |
| `C1_distance_correlation_heatmap.png/pdf` | 距離帯別相関ヒートマップ |
| `C2_promoter_vs_genebody_correlation.csv` | プロモーター/遺伝子本体別相関データ |
| `C2_promoter_vs_genebody_scatter.png/pdf` | 散布図 |
| `C2_promoter_vs_genebody_barplot.png/pdf` | 相関係数比較棒グラフ |
| `D1_sigma_motifs.csv` | σ70モチーフ一覧 |
| `D1_sigma_methylation_overlap.csv` | メチル化-σモチーフ重複データ |
| `D1_promoter_methylation_distribution.png/pdf` | コアプロモーターメチル化分布 |
| `D1_sigma_overlap_barplot.png/pdf` | σモチーフ重複棒グラフ |
| `motif_tss_sigma_spatial.csv` | モチーフ-TSS-σ因子空間データ |
| `motif_tss_position_distribution.png/pdf` | モチーフTSS相対位置分布 |
| `motif_sigma_distance_histogram.png/pdf` | モチーフ-σ距離ヒストグラム |
| `motif_promoter_architecture_map.png/pdf` | プロモーター構造上のモチーフ密度マップ |
| `motif_core_promoter_zoom.png/pdf` | コアプロモーター領域ズーム |
| `motif_methylated_vs_unmethylated_position.png/pdf` | メチル化/非メチル化モチーフ位置比較 |
| `motif_AAGCCCG_spatial_detail.csv` | AAGCCCG詳細空間データ |
| `motif_CCGG_spatial_detail.csv` | CCGG詳細空間データ |
| `sarp_methylation_mechanism_schematic.png/pdf` | SARPメカニズム模式図 |
| `sarp_zone_bgc_vs_nonbgc.png/pdf` | SARP zone BGC/非BGC比較 |
| `sarp_zone_bgc_heatmap.png/pdf` | BGC遺伝子SARP zoneヒートマップ |
| `sarp_zone_bgc_detail.csv` | BGC SARP zoneメチル化詳細 |
| `sarp_zone_correlation.csv` | SARP zone相関データ |
| `region_correlation_with_sarp_zone.png/pdf` | 全領域相関比較（SARP zone含む） |
| `temporal_motif_enrichment.csv` | タイムポイント別モチーフ含有率データ |
| `temporal_gained_lost_motifs.csv` | Gained/Lostサイト別モチーフ含有率 + Fisher検定 |
| `temporal_expression_associated_motifs.csv` | 発現カテゴリ×メチル化変動×モチーフのクロス集計 |
| `temporal_motif_heatmap.png/pdf` | タイムポイント別モチーフ含有率ヒートマップ |
| `temporal_gained_lost_motif_barplot.png/pdf` | Gained vs Lostサイトのモチーフ比較棒グラフ |
| `temporal_expression_motif_heatmap.png/pdf` | 発現変動×メチル化変動別モチーフヒートマップ |
| `temporal_key_motifs_barplot.png/pdf` | 主要モチーフのタイムポイント別含有率棒グラフ |
| `null_model_minus10_box.png/pdf` | -10 box枯渇の帰無モデル検証 |
| `null_model_sarp_zone.png/pdf` | SARP zone BGC保護の帰無モデル検証 |
| `null_model_results.csv` | 帰無モデルp値一覧 |
| `C1_distance_correlation_fdr.csv` | C1相関FDR補正結果 |
| `C2_promoter_genebody_fdr.csv` | C2相関FDR補正結果 |
| `temporal_gained_lost_motifs_fdr.csv` | Gained/Lost Fisher検定FDR補正結果 |
| `all_correlations_fdr.csv` | 全相関解析FDR補正結果（66検定） |
| `tf_methylation_enrichment.png/pdf` | 転写因子メチル化濃縮検定 |
| `tf_methylation_enrichment_summary.csv` | TF濃縮検定サマリ |
| `methylated_tf_list.csv` | メチル化TF一覧（124遺伝子） |
| `mtase_methylation_dynamics.png/pdf` | MTase発現-メチル化動態統合図 |
| `mtase_expression_summary.csv` | MTase発現サマリ |

---

## 8. 生物学的示唆と統合モデル

### 8.1 TSS周辺のメチル化ランドスケープモデル

```
Position:  -300    -200    -100    -50   -35  -10  TSS  +50    +200   +500
                                    │     │    │    │
Methylation: ███████████████████████│     │    │    │████████████████████
density:     [high]  [high]  [med]  │[med]│[low]│[low]│  [med]  [med]
                                    │     │    │    │
                              ┌─────┘     │    │    └────┐
                              │   -35 box │-10 │         │
                              │   (normal)│(dep)│  TSS   │
                              │           │    │  dip   │
                              └───────────┘    └────────┘
```

### 8.2 メチル化の位置依存的機能

1. **-10 box枯渇（進化的制約）**: σHrdBの結合に必須な-10 box配列はメチル化ターゲットを含みにくく進化的に保存
2. **TSS proximal正相関**: TSS直下流（5' UTR）のメチル化は転写活性と正に連動 → 活発に転写される遺伝子のDNAがメチル化されやすい（accessibility effect）
3. **上流プロモーター抑制（後期）**: -100~-50 bpの4mC、および全領域の6mAがT3で抑制的に作用 → 後期での遺伝子サイレンシングに関与

### 8.3 モチーフ依存的メチル化ダイナミクス（新規）

1. **CCGG（Dcm型4mC）**: T2で拡大→T3で大規模消失。消失が発現上昇と連動し、4mCの抑制解除メカニズムを示唆
2. **AAGCCCG（独自モチーフ）**: T2で消失→T3で少数が新規獲得。獲得が発現低下と連動し、標的的サイレンシングに関与する可能性
3. **TCGA（SARP関連モチーフ）**: T3で6mA新規獲得サイトに濃縮。二次代謝遷移期でのSARP関連エピジェネティック制御を示唆
4. **2つのメチル化酵素系の独立制御**: CCGGとAAGCCCGの逆転パターンは、異なるメチルトランスフェラーゼが培養フェーズに応じて異なる活性制御を受けることを示す

### 8.4 先行研究との整合性

- **Pisciotta 2018**: SCO1731(SC_RS10665)のT2/T3での発現低下は、m5Cが発育初期に最高で後期に低下するという報告と完全に一致
- **Pisciotta 2023**: BS-seqで検出されたm5CモチーフとNanoporeのm4Cモチーフの比較から、同一CCGG上の二重修飾系の存在を示唆
- **Fang 2022**: TFプロモーターのm4Cメチル化による間接的BGC制御モデルは、*S. coelicolor*では支持されなかった（TF=15.0% vs non-TF=17.1%, ns）
- **5-azacytidine研究群**: 脱メチル化→二次代謝活性化の方向性は、T3での4mC崩壊→遺伝子発現上昇パターンと一致
- **-10 box枯渇**: 「進化的選択圧」仮説は否定され、ATリッチ配列のモチーフ密度低下で説明可能

### 8.5 知見の堅牢性まとめ

| 知見 | FDR後 | 帰無モデル | 先行研究 | 最終評価 |
|------|:-----:|:---------:|:-------:|:-------:|
| TSS dip（メタジーンプロファイル） | N/A | N/A | 新規 | **堅牢** |
| T3での4mC崩壊 | N/A | N/A | 整合 | **堅牢** |
| Gained/Lostモチーフ逆転 | **3/36生残** | — | — | **堅牢** |
| SC_RS17645↔AAGCCCG対応 | N/A | — | 新規 | **注目** |
| m4C≠m5C二重修飾系仮説 | N/A | — | Pisciotta相補 | **注目** |
| 6mA gene body T3抑制 | **1/66生残** | — | 整合 | **やや堅牢** |
| -10 box枯渇 → 選択圧 | — | **否定(p=0.43)** | — | **棄却** |
| SARP zone BGC保護 | — | **否定(p=0.63)** | — | **棄却** |
| C1距離帯別相関(多数) | 0/48生残 | — | — | **弱い** |
| TFメチル化濃縮(Fangモデル) | — | — | **不支持** | **棄却** |

---

## 9. 論文ストーリー提案と次のステップ

### 9.1 論文に含めるべき内容（追加実験なし）

| 優先度 | 内容 | Figure候補 | 根拠 |
|:------:|------|-----------|------|
| **1** | TSS基準メタジーンメチル化プロファイル | Fig. 1 | *Streptomyces*初、記述的だが新規性あり |
| **2** | MTase発現↔モチーフ動態の対応 | Fig. 2 | SC_RS17645↔AAGCCCG、Dcm↔CCGG矛盾 |
| **3** | m4C/m5C二重修飾系仮説 | Fig. 3 | Dcm発現上昇 vs m4C崩壊の矛盾 |
| **4** | Gained/Lostモチーフ逆転パターン | Fig. 4 | FDR補正後も有意(3件) |
| Supp | FDR補正済み相関解析 | Table S1 | 堅牢性の担保 |
| Supp | 帰無モデル検証 | Fig. S1 | 棄却された仮説の明示 |

### 9.2 論文で主張できること / できないこと

**主張できること:**
1. *S. coelicolor*のTSS周辺にメチル化枯渇域が存在する（記述的事実）
2. 培養フェーズ間で異なるモチーフ（CCGG vs AAGCCCG）のメチル化が逆方向に動態する（FDR有意）
3. SC_RS17645の発現動態がAAGCCCGメチル化パターンと時間的に一致する（相関的証拠）
4. Dcm-like MTaseの発現上昇とm4C崩壊の矛盾がm4C≠m5Cの状況証拠を提供する

**主張できないこと:**
1. メチル化が発現の「原因」であること（相関のみ）
2. SC_RS17645がAAGCCCGを修飾すること（候補の提示のみ）
3. SARP結合のメチル化による阻害（帰無モデルで棄却）
4. -10 box枯渇の生物学的意義（配列組成で説明済み）

### 9.3 将来の実験提案（論文Discussion向け）

1. **SC_RS17645 KO株の作製**: AAGCCCG-6mAの消失を確認し、発現変動を解析（因果関係証明）
2. **Bisulfite-seq × Nanopore二重解析**: 同一サンプルでm5CとM4Cを分離定量（二重修飾系の直接証明）
3. **m4C酵素の同定**: CCGG上のm4Cを生成する酵素のスクリーニング（既知Dcm-like MTaseでないことの確認）

---

*レポート最終更新: 2026-02-04*
*スクリプト: tss_analyses.py, motif_tss_spatial_analysis.py, sarp_methylation_spatial.py, temporal_motif_expression_analysis.py, null_model_simulation.py, fdr_correction.py, tf_methylation_enrichment.py, mtase_methylation_dynamics.py*
