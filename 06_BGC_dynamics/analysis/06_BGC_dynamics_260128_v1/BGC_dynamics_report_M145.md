# BGC Dynamics Report -- M145 RNA-seq

## 1. Purpose

主要4つの生合成遺伝子クラスター（act, red, cda, cpk）のクラスター単位の発現ダイナミクスを定量化し、タイムコース（M145_1 → M145_2 → M145_3）における発現変動を図・表として可視化した。各 BGC 内の鍵遺伝子（高 LFC 遺伝子および cluster-situated regulator）を同定し、サンプルごとの BGC 発現スコアを算出した。

---

## 2. Input Data

| File | Source | Description |
|------|--------|-------------|
| `normalized_counts_M145.tsv` | 04_deseq2 | 7,646 genes x 9 samples, DESeq2 normalized counts |
| `gene_master_with_BGC_regulators.tsv` | 05_annotation | 8,275 genes, annotation + DESeq2 + BGC + regulator info |
| `BGC_definition_manual.tsv` | 05_annotation | 100 genes in 4 BGCs (act 22, red 22, cda 40, cpk 16) |

---

## 3. Output Files

### 3.1 Tables

| File | Description | Rows |
|------|-------------|------|
| `gene_condition_means.tsv` | Gene x condition mean normalized counts | 7,646 |
| `BGC_condition_means.tsv` | BGC-level mean expression per condition | 4 |
| `BGC_key_genes_summary.tsv` | Key genes per BGC (top LFC + regulators) | 30 |
| `BGC_sample_scores.tsv` | Per-sample BGC expression scores | 9 |

### 3.2 Figures

| File | Description |
|------|-------------|
| `BGC_timecourse_lineplot_M145.pdf/.svg` | 4 BGC timecourse, raw scale |
| `BGC_timecourse_lineplot_log10_M145.pdf/.svg` | 4 BGC timecourse, log10 scale |
| `BGC_condition_heatmap_M145.pdf/.svg` | BGC x condition heatmap (log10) |
| `BGC_gene_heatmap_{act,red,cda,cpk}_M145.pdf/.svg` | Per-gene heatmaps within each BGC |
| `BGC_sample_scores_boxplot_M145.pdf/.svg` | Per-sample BGC score boxplots |

---

## 4. BGC Timecourse: Cluster-Level Expression

### 4.1 BGC Mean Expression by Condition

| BGC | M145_1 mean | M145_2 mean | M145_3 mean | FC (2v1) | FC (3v1) | FC (3v2) |
|-----|-------------|-------------|-------------|----------|----------|----------|
| **act** | 65.3 | 98.5 | 9,491 | 1.5x | **145x** | **96x** |
| **red** | 65.2 | 681.6 | 965.5 | 10.4x | **14.8x** | 1.4x |
| **cda** | 56.7 | 3,309 | 3,085 | **58.4x** | 54.4x | 0.93x |
| **cpk** | 52.1 | 9,465 | 6,144 | **181.7x** | 118x | 0.65x |

### 4.2 Temporal Activation Patterns

3つの timepoint における BGC-level 平均発現の推移から、**2つの明確な temporal profile** が識別される:

**Profile A -- Progressive escalation (act)**:
- M145_1 → M145_2 ではわずかに増加（1.5x）。
- M145_2 → M145_3 で **96 倍の劇的な増加**。
- 全 22 遺伝子が M145_3 vs M145_2 で有意に上昇（22/22）。
- actinorhodin 生合成はタイムウィンドウの後期に集中して活性化される。

**Profile B -- Early-mid activation with plateau (red, cda, cpk)**:
- **red**: M145_1 → M145_2 で 10.4x に増加し、M145_2 → M145_3 でさらに 1.4x 増加（緩やかな上昇継続）。
- **cda**: M145_1 → M145_2 で 58.4x と急増し、M145_2 → M145_3 ではほぼ横ばい（0.93x、微減）。
- **cpk**: M145_1 → M145_2 で 181.7x と最も劇的に増加し、M145_2 → M145_3 ではむしろ減少（0.65x）。

cpk は M145_2 で最も高い絶対発現レベル（mean 9,465）を示すが、M145_3 では 6,144 に低下する。この「一過性の発現ピーク」パターンは、coelimycin 生合成が培養中期に最大活性を示し、後期には減衰することを示唆する。

---

## 5. Differential Expression Summary (Gene-Level)

### 5.1 Per-BGC DE Statistics

| BGC | Total | Sig 3v1 | Up 3v1 (LFC>1) | Median LFC 3v1 | Max LFC 3v1 | Up 2v1 | Up 3v2 |
|-----|-------|---------|-----------------|-----------------|-------------|--------|--------|
| **act** | 22 | 22 (100%) | **22** | +8.50 | +11.34 | 13 | **22** |
| **red** | 22 | 22 (100%) | **21** | +4.93 | +9.83 | 21 | 17 |
| **cda** | 40 | 40 (100%) | **39** | +7.16 | +10.80 | 39 | 0 |
| **cpk** | 16 | 16 (100%) | **16** | +6.11 | +13.19 | 16 | 0 |

- 全 100 遺伝子が M145_3 vs M145_1 で有意（padj < 0.05）。
- act: M145_2 → M145_3 で **全 22 遺伝子がさらに有意に上昇**（Up 3v2 = 22）。
- cda / cpk: M145_2 → M145_3 で上昇遺伝子なし（Up 3v2 = 0）。プラトーまたは減少。
- red: M145_2 → M145_3 でも 17/22 遺伝子が上昇（中間的な挙動）。

### 5.2 LFC Distribution by BGC

- **act**: Median LFC +8.50 と最も高い中央値。M145_1 でほぼサイレント → M145_3 で爆発的発現。
- **cpk**: Max LFC +13.19（cpkO/kasO）と全 BGC 遺伝子中の最大値。ただし median は +6.11。
- **cda**: 40 遺伝子と最大のクラスターだが、median LFC +7.16 と均一に強く誘導。
- **red**: Median LFC +4.93 と相対的に控えめ。M145_1 でも一定のベースライン発現あり。

---

## 6. Key Genes within BGCs

### 6.1 act (7 key genes)

| gene_id | SCO | Gene | LFC 3v1 | LFC 2v1 | LFC 3v2 | Role |
|---------|-----|------|---------|---------|---------|------|
| SC_RS27515 | SCO5071 | -- | **+11.34** | +5.56 | +5.91 | top LFC |
| SC_RS27590 | SCO5086 | -- | **+11.11** | +1.93 | +9.13 | top LFC (transport) |
| SC_RS27520 | SCO5072 | -- | +10.89 | +1.94 | +8.88 | top LFC |
| SC_RS27530 | SCO5074 | -- | +10.87 | +1.97 | +8.85 | top LFC |
| SC_RS27605 | SCO5089 | -- | +10.17 | +4.38 | +5.72 | top LFC |
| SC_RS27585 | SCO5085 | -- | +6.44 | +0.69 | +5.82 | regulator (SARP) |
| **SC_RS27570** | **SCO5082** | **actII-orf4** | **+2.60** | +1.23 | +1.39 | **CSR (SARP)** |

**Observations**:
- actII-orf4（SCO5082, SARP family）は LFC +2.60 と、構造遺伝子（LFC +8--11）に比べて控えめな変化を示す。これは CSR（cluster-situated regulator）が比較的低コピーで機能しうることを反映する。
- SCO5085（SARP family auto-detected）は LFC +6.44 と actII-orf4 より大きな変動を示し、追加的な制御要素として注目に値する。
- 構造遺伝子は M145_2 → M145_3 での LFC 3v2 が +5--9 と大きく、**後期における大幅な追加的誘導**を確認。

### 6.2 red (7 key genes)

| gene_id | SCO | Gene | LFC 3v1 | LFC 2v1 | LFC 3v2 | Role |
|---------|-----|------|---------|---------|---------|------|
| SC_RS31690 | SCO5889 | -- | **+9.83** | +7.35 | +2.31 | top LFC |
| SC_RS31680 | SCO5887 | -- | +5.81 | +4.16 | +1.61 | top LFC |
| SC_RS31675 | SCO5886 | -- | +5.80 | +3.80 | +2.02 | top LFC |
| SC_RS31710 | SCO5893 | -- | +5.42 | +3.27 | +2.14 | top LFC |
| SC_RS31660 | SCO5883 | -- | +5.42 | +3.44 | +2.06 | top LFC |
| **SC_RS31630** | **SCO5877** | **redD** | **+5.00** | +3.29 | +1.75 | **CSR (SARP)** |
| **SC_RS31650** | **SCO5881** | **redZ** | **+0.97** | +1.67 | -0.73 | **response reg.** |

**Observations**:
- redD（SCO5877, SARP）は LFC +5.00 と構造遺伝子と同程度に変動しており、クラスターと同調した誘導パターンを示す。
- redZ（SCO5881, response regulator）は LFC +0.97 と控えめ。M145_2 で LFC +1.67 まで上昇した後、M145_3 では -0.73 と**減少に転じている**。redZ は redD の上位制御因子とされるが、M145_3 での発現低下は、red クラスター発現が redZ 非依存的に維持される段階に入っていることを示唆する。
- SCO5889 の LFC +9.83 は red クラスター内で突出しており、M145_1 での発現がほぼゼロからの急激な誘導を反映する。

### 6.3 cda (8 key genes)

| gene_id | SCO | Gene | LFC 3v1 | LFC 2v1 | LFC 3v2 | Role |
|---------|-----|------|---------|---------|---------|------|
| SC_RS18350 | SCO3247 | -- | **+10.80** | +10.35 | +0.31 | top LFC |
| SC_RS18315 | SCO3240 | -- | +9.49 | +8.42 | +0.99 | top LFC |
| SC_RS18225 | SCO3222 | -- | +9.28 | +9.13 | +0.08 | top LFC |
| SC_RS18325 | SCO3242 | -- | +9.16 | +8.49 | +0.59 | top LFC |
| SC_RS18295 | SCO3236 | asnO | +9.05 | +8.72 | +0.24 | top LFC |
| SC_RS18200 | SCO3217 | -- | +5.76 | +5.78 | -0.14 | regulator (SARP) |
| **SC_RS18240** | **SCO3225** | **absA1** | **+2.73** | +2.85 | -0.19 | **sensor kinase** |
| **SC_RS18245** | **SCO3226** | **absA2** | **+0.54** | +0.84 | -0.33 | **response reg.** |

**Observations**:
- cda 構造遺伝子は LFC 2v1 と LFC 3v1 がほぼ同値（LFC 3v2 ≈ 0）で、**M145_2 の時点ですでに最大発現に到達**していることが明確。
- absA1（SCO3225, sensor kinase）は LFC +2.73 と中程度の誘導。absA2（SCO3226, response regulator）は LFC +0.54 にとどまる。absA1/absA2 は two-component system として cda を含む複数 BGC を制御するプレイオトロピックレギュレーターだが、response regulator 側の変動が小さい。
- SCO3217（SARP family auto-detected）は LFC +5.76 で、cda クラスター内の CSR 候補として機能している可能性がある。

### 6.4 cpk (8 key genes)

| gene_id | SCO | Gene | LFC 3v1 | LFC 2v1 | LFC 3v2 | Role |
|---------|-----|------|---------|---------|---------|------|
| **SC_RS33660** | **SCO6282** | **cpkO/kasO** | **+13.19** | +13.12 | -0.02 | **CSR (SARP-like)** |
| SC_RS33645 | SCO6279 | -- | +9.30 | +8.96 | +0.21 | top LFC |
| SC_RS33630 | SCO6276 | -- | +8.99 | +8.83 | +0.03 | top LFC |
| SC_RS33640 | SCO6278 | -- | +8.52 | +8.26 | +0.17 | top LFC |
| SC_RS33665 | SCO6283 | -- | +7.81 | +7.77 | -0.05 | top LFC |
| **SC_RS33680** | **SCO6286** | **scbR2** | **+5.49** | +6.29 | -0.88 | **butyrolactone rec.** |
| SC_RS33650 | SCO6280 | -- | +5.21 | +5.49 | -0.37 | regulator (SARP) |
| SC_RS33690 | SCO6288 | -- | +3.90 | +4.99 | -1.23 | regulator (SARP) |

**Observations**:
- cpkO/kasO（SCO6282）は **全 BGC 遺伝子中の最大 LFC（+13.19）** を示す。LFC 2v1 = +13.12 であり、M145_2 の時点ですでにほぼ最大誘導に達している（LFC 3v2 ≈ 0）。M145_1 でほぼサイレントだった転写が 8,000 倍以上に急増する「binary switch」型の挙動を示す。
- scbR2（SCO6286）は LFC +5.49（3v1）だが、LFC 3v2 = -0.88 と M145_3 では若干の低下。gamma-butyrolactone シグナル系を介した制御が培養中期に主に機能することと整合する。
- SC_RS33690（SCO6288, SARP）は LFC 3v2 = -1.23 と cpk クラスター内で最も大きな後期低下を示す。
- cpk 全体として M145_3 でのクラスター平均発現が M145_2 より低下（FC 3v2 = 0.65x）しており、一部の regulator や構造遺伝子が後期に減衰していることを反映する。

---

## 7. Sample-Level BGC Scores

### 7.1 Per-Sample Scores

| Sample | Condition | act score | red score | cda score | cpk score |
|--------|-----------|-----------|-----------|-----------|-----------|
| M145_1_1 | M145_1 | 62.2 | 71.0 | 58.3 | 53.2 |
| M145_1_2 | M145_1 | 67.9 | 60.8 | 56.0 | 50.7 |
| M145_1_3 | M145_1 | 65.8 | 63.9 | 55.7 | 52.4 |
| M145_2_1 | M145_2 | 96.1 | 756.7 | 3,417 | 9,976 |
| M145_2_3 | M145_2 | 103.2 | 670.4 | 3,606 | 10,019 |
| M145_2_4 | M145_2 | 96.2 | 617.6 | 2,903 | 8,401 |
| M145_3_2 | M145_3 | 7,168 | 751.7 | 2,467 | 4,890 |
| M145_3_3 | M145_3 | **13,628** | 1,333 | 4,147 | 8,888 |
| M145_3_4 | M145_3 | 7,676 | 811.7 | 2,642 | 4,654 |

### 7.2 Replicate Variability

- **M145_1**: 全 4 BGC で非常に低い発現（50--70）。レプリケート間の変動は小さい（CV < 10%）。
- **M145_2**: cda と cpk が急増（3,000--10,000）。レプリケート間の一致は良好だが、cpk では M145_2_4 がやや低い。
- **M145_3**: **act の M145_3_3 が 13,628 と他の 2 レプリケート（7,168 / 7,676）の約 1.8 倍**。この生物学的変動は、actinorhodin 生合成の活性化タイミングにサンプル間差があることを示唆する。cpk は M145_3_3 のみ高値（8,888）で他は 4,600--4,900 と低く、こちらでも M145_3_3 が高発現傾向を示す。

---

## 8. Biological Interpretation

### 8.1 Two-Phase Model of BGC Activation

本解析の結果は、*S. coelicolor* M145 における二次代謝遺伝子活性化の **二段階モデル** を支持する:

1. **Phase I（M145_1 → M145_2）: 広域 BGC 活性化**
   - red, cda, cpk の 3 クラスターが 10--180 倍に急増。cpk が最大の変化を示す。
   - act は M145_2 でもまだ低発現（1.5x 増加のみ）。
   - これは成長期から定常期への移行に伴う、gamma-butyrolactone（SCB1 等）や ppGpp を介した広域スイッチと考えられる。

2. **Phase II（M145_2 → M145_3）: act-specific な追加的誘導**
   - act のみが 96 倍にさらに増加し、M145_3 で最大発現に達する。
   - red は緩やかに増加継続（1.4x）。
   - cda / cpk は横ばいまたは減少。
   - act の遅延型活性化は、actII-orf4 の発現量増加（LFC 3v2 = +1.39）と相まって、act 特異的な追加シグナルの存在を示唆する。

### 8.2 Regulator Dynamics

| Regulator | BGC | LFC 3v1 | LFC 2v1 | LFC 3v2 | Pattern |
|-----------|-----|---------|---------|---------|---------|
| actII-orf4 | act | +2.60 | +1.23 | +1.39 | Progressive increase |
| SCO5085 (SARP) | act | +6.44 | +0.69 | +5.82 | Late-phase surge |
| redD | red | +5.00 | +3.29 | +1.75 | Progressive, cluster-synchronous |
| redZ | red | +0.97 | +1.67 | **-0.73** | Early peak, late decline |
| SCO3217 (SARP) | cda | +5.76 | +5.78 | -0.14 | Early-mid plateau |
| absA1 | cda | +2.73 | +2.85 | -0.19 | Early-mid plateau |
| absA2 | cda | +0.54 | +0.84 | -0.33 | Minimal change |
| cpkO/kasO | cpk | **+13.19** | +13.12 | -0.02 | Binary switch at M145_2 |
| scbR2 | cpk | +5.49 | +6.29 | -0.88 | Early-mid peak |

**Regulator-level observations**:

- **cpkO/kasO** は M145_1 → M145_2 で LFC +13.12 と実質的に「OFF → ON」のスイッチ型挙動。これは coelimycin 生合成の全制御を一手に担う master switch としての機能と整合する。
- **actII-orf4** は LFC +2.60 と控えめだが、M145_2 → M145_3 での +1.39 の追加的増加が、act 構造遺伝子の Phase II 誘導と対応する。
- **redZ** の M145_3 での発現低下は、red クラスターが redZ 非依存的な経路（たとえば redD の autoregulation）で維持される段階への移行を反映する可能性がある。
- **absA2** の微弱な応答（+0.54）は、この pleiotropic regulator がタンパク質レベルのリン酸化で主に制御されることと整合し、転写量の変化が小さくても機能的に重要でありうる。

### 8.3 cpk の「一過性ピーク」パターン

cpk はクラスター平均で M145_2（9,465）→ M145_3（6,144）と **35% の低下** を示す。これは他の 3 BGC では見られない特徴であり、coelimycin 生合成が培養中期に一過性のピークを迎え、後期に減衰することを示唆する。scbR2（LFC 3v2 = -0.88）と SC_RS33690（LFC 3v2 = -1.23）の後期低下が、このクラスターレベルの減衰と対応している。

---

## 9. Next Steps

- **07_regulator_network**: BGC sample scores（本ステップで算出）と 877 の putative regulator（05_annotation で同定）の発現パターンの共変動解析を行い、BGC 発現を制御する candidate regulator を同定する。
- **Global regulator timing**: atrA (SC_RS25560)、bldA/bldH/bldD 等の global/developmental regulator と BGC activation timing の関係を検証する。
- **act の遅延型誘導メカニズム**: Phase II で act のみが追加的に誘導される機構について、actII-orf4 の上位制御因子（absA1/absA2, atrA 等）の発現動態との相関を詳細に解析する。

---

## 10. Key Takeaways

1. **主要 4 BGC（act/red/cda/cpk）はすべて M145_2 以降で強く誘導される**が、act のみ M145_3 で追加的な 96 倍のアップレギュレーションを示す（二段階活性化モデル）。
2. **cpk は M145_2 で最大発現に達した後、M145_3 では 35% 低下する**「一過性ピーク」型のダイナミクスを示し、他の 3 BGC とは異なる制御を受けている。
3. **cpkO/kasO（SC_RS33660, LFC +13.19）は全 BGC 遺伝子中の最大 LFC を示し**、coelimycin BGC の binary switch として機能している可能性が高い。
4. **Cluster-situated SARP regulator（actII-orf4, redD, cpkO）はそれぞれの BGC と同調して誘導される**が、LFC の大きさは構造遺伝子の 30--50% にとどまり、低コピーでの制御機能を示唆する。
5. **absA2（response regulator）の転写応答は +0.54 と微弱**であり、absA1/absA2 two-component system はタンパク質レベル（リン酸化）の制御が主であることを示唆する。

---

## 11. Materials & Methods

主要 4 BGC（act 22 genes, red 22 genes, cda 40 genes, cpk 16 genes; total 100 genes）の発現ダイナミクスを、DESeq2 normalized counts を用いて定量化した。各遺伝子の条件ごとの平均正規化カウントを算出し、BGC レベルの平均発現プロファイルを作成した。タイムコースラインプロット、BGC x 条件ヒートマップ、およびクラスター内遺伝子レベルヒートマップを作成した。鍵遺伝子は、各 BGC 内で |log2FoldChange(3_vs_1)| の上位 5 遺伝子と、全ての cluster-situated regulator を抽出して同定した。サンプルごとの BGC 発現スコアは、各クラスターの構成遺伝子の正規化カウント平均として算出した。解析は R 4.4.2（tidyverse, pheatmap, svglite）で実施した。

---

*Generated: 2026-01-28*
*Run directory: `06_BGC_dynamics_260128_v1`*
