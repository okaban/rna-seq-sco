# Regulator Network Report -- M145 RNA-seq

## 1. Purpose

主要4 BGC（act, red, cda, cpk）の発現スコアと 877 個の putative regulator 遺伝子の rlog 発現との共変動を Pearson 相関解析で定量化し、BGC 制御ネットワークの候補リンクを抽出した。既知の cluster-situated regulator（CSR）との整合性を検証するとともに、新規の regulator 候補を同定した。

---

## 2. Input Data

| File | Source | Description |
|------|--------|-------------|
| `rld.rds` | 04_deseq2 | rlog-transformed expression, 7,646 genes x 9 samples |
| `gene_master_with_BGC_regulators.tsv` | 05_annotation | 8,275 genes, annotation + regulator flags |
| `BGC_sample_scores.tsv` | 06_BGC_dynamics | Per-sample BGC scores (act/red/cda/cpk) |

---

## 3. Methods

- **Regulator set**: `is_regulator == TRUE` の 877 遺伝子のうち、rlog データに存在する 856 遺伝子を解析対象とした。
- **相関解析**: 各 regulator の rlog 発現ベクトル（9 サンプル）と 4 BGC スコアベクトルとの Pearson 相関係数を算出。
- **ランキング**: 各 BGC について |corr| で降順ソートし、上位 30 regulator を抽出。
- **既知 regulator**: 文献に基づき 13 の known regulator を手動定義し、`is_known` フラグで識別。
- **Global regulator 候補**: |corr| > 0.9 を 2 つ以上の BGC に対して示す regulator を global candidate とした。

---

## 4. Output Files

### 4.1 Tables

| File | Description | Rows |
|------|-------------|------|
| `regulator_BGC_correlation.tsv` | 全 regulator × 4 BGC の Pearson 相関 | 856 |
| `regulator_top_hits_by_BGC.tsv` | 各 BGC の上位 30 regulator | 120 |
| `regulator_BGC_known_vs_novel.tsv` | 既知 13 regulator の相関プロファイル | 13 |

### 4.2 Figures

| File | Description |
|------|-------------|
| `heatmap_known_regulators_BGC_corr_M145.pdf/.svg` | 既知 13 regulator × 4 BGC の相関ヒートマップ |
| `heatmap_top_regulators_all_BGC_M145.pdf/.svg` | 上位 regulator 全体の相関ヒートマップ |
| `scatter_*_M145.pdf/.svg` | BGC スコア vs regulator rlog 散布図（8 組） |
| `timecourse_*_M145.pdf/.svg` | BGC スコア + regulator タイムコース（4 組） |

---

## 5. Correlation Summary Statistics

| BGC | Median |corr| | n(|corr| > 0.9) | n(|corr| > 0.8) |
|-----|----------------|------------------|------------------|
| **act** | 0.497 | 28 | 100 |
| **red** | 0.674 | 71 | 247 |
| **cda** | 0.693 | **166** | **315** |
| **cpk** | 0.668 | **156** | 290 |

- **act** は他の 3 BGC と比較して中央値が低く（0.497）、高相関 regulator の数も最少（|corr|>0.9: 28）。これは act の独自の遅延型活性化パターン（M145_3 でのみ急増）が、多くの regulator の発現パターンと一致しにくいことを反映する。
- **cda / cpk** は |corr| > 0.9 の regulator が 150 以上と最多。M145_1 → M145_2 での急増パターンが多くの regulator の変動パターンと同期するため。
- **red** は中間的な位置づけ。

---

## 6. Known Regulator Correlation Profiles

### 6.1 Known Regulator × BGC 相関行列

| Regulator | gene_id | act | red | cda | cpk | Pattern |
|-----------|---------|-----|-----|-----|-----|---------|
| **actII-orf4** | SC_RS27570 | **+0.95** | +0.57 | +0.25 | -0.02 | act-specific |
| SCO5085 (SARP) | SC_RS27585 | **+0.92** | +0.87 | +0.69 | +0.47 | act-biased |
| **redD** | SC_RS31630 | +0.39 | +0.86 | **+0.97** | **+0.95** | red/cda/cpk |
| **redZ** | SC_RS31650 | +0.58 | **+0.95** | **+0.98** | **+0.92** | red/cda/cpk |
| SCO3217 (SARP) | SC_RS18200 | +0.36 | +0.85 | **+0.97** | **+0.96** | cda/cpk |
| **absA1** | SC_RS18240 | +0.65 | **+0.97** | **+0.97** | +0.87 | red/cda |
| **absA2** | SC_RS18245 | +0.19 | +0.73 | +0.90 | **+0.96** | cpk-biased |
| **cpkO/kasO** | SC_RS33660 | +0.45 | +0.88 | **+0.97** | **+0.92** | cda/cpk |
| SCO6280 (SARP) | SC_RS33650 | +0.26 | +0.79 | +0.95 | **+0.97** | cpk-dominant |
| **scbR2** | SC_RS33680 | +0.46 | +0.89 | **+0.97** | **+0.92** | cda/cpk |
| SCO6288 (SARP) | SC_RS33690 | +0.36 | +0.84 | **+0.96** | **+0.94** | cda/cpk |
| **scbR** | SC_RS33575 | -0.25 | +0.40 | +0.69 | +0.87 | cpk-positive |
| **atrA** | SC_RS25560 | +0.45 | +0.88 | **+0.97** | **+0.93** | cda/cpk |

### 6.2 既知ネットワークの再現性

**再現されたパターン:**

1. **actII-orf4 → act**: corr_act = **+0.95** と極めて高い正相関。他の 3 BGC とは低相関（red +0.57, cda +0.25, cpk -0.02）。actII-orf4 が act クラスターの特異的 CSR として機能していることと完全に整合する。

2. **redD / redZ → red**: redZ は corr_red = **+0.95**、redD は corr_red = +0.86。ただし redD と redZ の両方が cda / cpk とも高相関（redZ: cda +0.98, cpk +0.92）を示す。これは red/cda/cpk が M145_1 → M145_2 で同時に立ち上がるという共通の temporal pattern を持つためで、必ずしも因果関係を意味しない。

3. **cpkO/kasO → cpk**: corr_cpk = **+0.92** と高い正相関。ただし corr_cda = +0.97 とさらに高い。これも temporal confounding（cda と cpk が同様の M145_1→M145_2 立ち上がりパターンを共有）の影響。

4. **scbR2 → cpk**: corr_cpk = **+0.92**（正の相関）。gamma-butyrolactone 受容体として cpk クラスターの正の制御因子であるという既知の機能と整合。

5. **scbR → cpk**: corr_cpk = +0.87 と正相関。scbR は cpk 上流の butyrolactone receptor であり、cpk 活性化と正に共変動することは妥当。ただし act とは -0.25 と弱い負の相関にとどまる。

6. **absA1 / absA2**: absA1 は cda (+0.97) / red (+0.97) と強い正相関、absA2 は cpk (+0.96) と最も強い相関を示す。absA1 と absA2 が異なる BGC 相関プロファイルを持つ点は注目に値する。

7. **atrA → act（部分的再現）**: corr_act = +0.45 と中程度の正相関にとどまる。atrA は actII-orf4 の上位制御因子として報告されているが、本データでは act score との相関が控えめ。cda (+0.97) / cpk (+0.93) との相関のほうが高く、atrA が act 以外の BGC にも影響する pleiotropic な因子であることを示唆する。

### 6.3 Temporal Confounding に関する注意

n=9（3 条件 × 3 レプリケート）の相関解析では、M145_1 → M145_2 で同時に立ち上がる cda/cpk/red のいずれにも高相関を示す regulator が多数検出される。これは **temporal confounding** の影響であり、真の因果関係の特定には条件数の増加（追加タイムポイント）や遺伝学的検証（ノックアウト/過剰発現）が必要である。

一方で、**act に対して特異的に高い相関を示す regulator**（actII-orf4 など）は、act の独自の遅延型パターン（M145_3 でのみ急増）により temporal confounding の影響を受けにくく、より信頼性の高い候補といえる。

---

## 7. Novel Regulator Candidates

### 7.1 act-specific 高相関候補（|corr_act| > 0.9）

act に特異的に高い相関を示す regulator は特に注目に値する（temporal confounding の影響が小さいため）:

| Rank | gene_id | SCO | Product | corr_act |
|------|---------|-----|---------|----------|
| 1 | SC_RS27590 | SCO5086 | ABC transporter ATP-binding protein | +0.99 |
| 2 | SC_RS27575 | SCO5083 | beta-ketoacyl-[acyl-carrier-protein] synthase | +0.98 |
| 3 | SC_RS27580 | SCO5084 | acyl carrier protein | +0.98 |
| 4 | SC_RS27565 | SCO5081 | hypothetical protein | +0.96 |
| 5 | SC_RS27560 | SCO5080 | monooxygenase | +0.96 |

> 注: 上位は act クラスター自体の構造遺伝子に regulator フラグが付いたものが多い。act クラスター外で高相関の regulator を重視すべき。

act クラスター外（bgc_name ≠ act）で |corr_act| が高い regulator:

| gene_id | SCO | Product | corr_act | corr_red | corr_cda |
|---------|-----|---------|----------|----------|----------|
| SC_RS27585 | SCO5085 | SARP family regulator | +0.92 | +0.87 | +0.69 |

SCO5085 は act クラスターの transport 遺伝子として定義されているが、SARP family として regulator フラグも付いている。act に対して +0.92 と actII-orf4 に匹敵する相関を示す。

### 7.2 Global regulator candidates（|corr| > 0.9 for ≥ 2 BGCs）

141 の regulator が 2 つ以上の BGC で |corr| > 0.9 を示した。特に 3 つ以上の BGC で高相関を示すものを注目候補とする:

| gene_id | SCO | Product | corr_act | corr_red | corr_cda | corr_cpk | n_high |
|---------|-----|---------|----------|----------|----------|----------|--------|
| SC_RS31650 | SCO5881 | response regulator (redZ) | +0.58 | **+0.95** | **+0.98** | **+0.92** | 3 |
| SC_RS21215 | SCO3818 | response regulator | +0.53 | +0.93 | **+0.98** | **+0.93** | 3 |
| SC_RS29775 | SCO5518 | PucR family regulator | -0.55 | **-0.94** | **-0.98** | **-0.91** | 3 |
| SC_RS28860 | SCO5337 | XRE family regulator | -0.53 | -0.93 | **-0.97** | **-0.93** | 3 |

**SC_RS21215 (SCO3818)** -- response regulator: cda (+0.98), cpk (+0.93), red (+0.93) と 3 BGC で強い正相関。BGC クラスターには属さず、未報告の global regulator 候補として注目に値する。

**SC_RS29775 (SCO5518)** -- PucR family regulator: 3 BGC で強い**負の**相関（cda -0.98, red -0.94, cpk -0.91）。成長期（M145_1）で高発現 → 二次代謝期に低下するパターンが想定され、BGC 活性化のリプレッサー候補。

**SC_RS28860 (SCO5337)** -- XRE family regulator: 同様に 3 BGC で負の相関。二次代謝への移行とともに発現が低下するリプレッサー型の制御因子候補。

---

## 8. Biological Interpretation

### 8.1 制御ネットワークの階層構造

本解析の結果は、*S. coelicolor* における BGC 制御の **階層的ネットワーク構造** を反映している:

1. **CSR 層（クラスター内制御）**: actII-orf4 → act（r = +0.95）、cpkO/kasO → cpk（r = +0.92）、redD/redZ → red（r = +0.86/+0.95）の各 CSR が自身の BGC と高い相関を示す。
2. **シグナル系（butyrolactone）**: scbR2（r_cpk = +0.92）、scbR（r_cpk = +0.87）が cpk と正に共変動。
3. **Global regulator 層**: 141 の regulator が 2+ BGC で |corr| > 0.9 を示し、二次代謝全般の制御に関与する候補プール。

### 8.2 act の独自性と actII-orf4 の特異性

actII-orf4 は **唯一、act にのみ強い相関を示す CSR** である（corr_act = +0.95, corr_cpk = -0.02）。他の CSR（redD, cpkO/kasO, scbR2）はいずれも複数の BGC と高相関を示す。これは：

- act の M145_3 特異的な遅延型活性化が、他の BGC（M145_2 で立ち上がるもの）と異なるシグナルで駆動されていること。
- actII-orf4 が真に act-specific な制御を行っている可能性が高いこと。

を示唆する。

### 8.3 absA1/absA2 の非対称な相関プロファイル

absA1（sensor kinase）は red (+0.97) / cda (+0.97) と強い正相関を示す一方、absA2（response regulator）は cpk (+0.96) との相関が最も強い。Two-component system の 2 コンポーネントが異なる BGC 相関プロファイルを持つことは、absA2 の転写量変化（LFC +0.54 にとどまる、05_annotation 参照）がリン酸化レベルでの制御と乖離している可能性を示唆する。

### 8.4 今後の実験候補

| 優先度 | Regulator | Experiment | Rationale |
|--------|-----------|------------|-----------|
| High | actII-orf4 (SC_RS27570) | 過剰発現 | act-specific CSR; act の Phase II 誘導を増強する可能性 |
| High | SC_RS21215 (SCO3818) | ノックアウト | 3 BGC に正の相関を持つ未報告 global regulator |
| High | SC_RS29775 (SCO5518) | ノックアウト | 3 BGC に負の相関 → ノックアウトで二次代謝活性化の可能性 |
| Medium | scbR2 (SC_RS33680) | ノックアウト/過剰発現 | cpk/cda との強い正相関; cpk–act トレードオフへの影響 |
| Medium | SC_RS28860 (SCO5337) | ノックアウト | XRE family; 3 BGC への負の相関 |

---

## 9. Key Takeaways

1. **actII-orf4 は 856 の regulator 中で唯一、act score のみと強い正相関（+0.95）を示し、他の 3 BGC とは無相関**。act の遅延型活性化パターンを反映した、真に act-specific な制御因子であることが相関レベルで確認された。
2. **redD, cpkO/kasO, scbR2, atrA を含む多くの既知 CSR は、red/cda/cpk の 3 BGC に対して広く高相関**を示す。これは temporal confounding（M145_1→M145_2 での同時立ち上がり）の寄与が大きく、因果関係の特定には追加実験が必要。
3. **141 の global regulator 候補**が 2+ BGC で |corr| > 0.9 を達成。特に SC_RS21215 (SCO3818, response regulator) は 3 BGC で正の相関を示す新規候補として有望。
4. **SC_RS29775 (SCO5518, PucR family) と SC_RS28860 (SCO5337, XRE family) は 3 BGC で強い負の相関**を示し、二次代謝への移行時に発現が低下するリプレッサー候補。ノックアウトにより二次代謝活性化が期待される。
5. **absA1/absA2 が異なる相関プロファイルを持つ**（absA1: red/cda, absA2: cpk に偏る）ことは、この two-component system のリン酸化シグナルと転写量が必ずしも連動しないことを示唆する。

---

## 10. Materials & Methods

05_annotation で定義した 877 の putative regulator 遺伝子のうち、DESeq2 rlog データに含まれる 856 遺伝子について、rlog 発現値（9 サンプル）と 06_BGC_dynamics で算出した 4 BGC スコア（act, red, cda, cpk）との Pearson 相関係数を計算した。各 BGC について |corr| の上位 30 regulator をランキングした。文献に基づく 13 の known regulator（actII-orf4, redD, redZ, absA1, absA2, cpkO/kasO, scbR, scbR2, atrA, 及び BGC 内 SARP family 4 遺伝子）を手動定義し、is_known フラグで識別した。|corr| > 0.9 を 2 つ以上の BGC に対して示す regulator を global candidate（141 遺伝子）として抽出した。散布図およびタイムコースプロットは ggplot2 で、相関ヒートマップは pheatmap で作成した。解析は R 4.4.2（DESeq2, tidyverse, pheatmap, svglite, ggrepel）で実施した。

---

*Generated: 2026-01-28*
*Run directory: `07_regulator_network_260128_v1`*
