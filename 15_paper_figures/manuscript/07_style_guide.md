# Style Guide: M145 Gatekeeper Model Paper

**Purpose**: このファイルは論文執筆時に全セクションで参照し、用語・数値・表記の一貫性を保つためのもの。
各セッション冒頭でこのファイルを読み込むこと。

---

## 1. 用語統一表 (Terminology)

| 使用する表現 | 使用しない表現 | 備考 |
|-------------|--------------|------|
| Gatekeeper model | gatekeeper mechanism, gatekeeper hypothesis | 固有名詞として大文字 |
| shielded regulators | protected genes, methylation-free genes | 小文字、名詞として使う |
| exposed regulators | unprotected genes, methylation-accessible | 小文字 |
| protection zone | exclusion zone, methylation-free zone | |
| R-M system | restriction-modification system | 初出時のみフルスペル |
| GCCGGC (N4-methylcytosine) | GCCGGC 4mC | 初出時のみ "N4-methylcytosine (4mC)" |
| AAGCCCG (N6-methyladenine) | AAGCCCG 6mA | 初出時のみ "N6-methyladenine (6mA)" |
| activation bloc | activation module, upregulated cluster | "bloc" は仏語由来だが生物学で通用 |
| repression bloc | repression module, downregulated cluster | |
| Simpson's paradox | Simpson's effect, ecological fallacy | |
| vegetative-to-developmental transition | growth-to-development switch | |
| linear chromosome | Streptomyces chromosome | core/arm構造を強調する文脈で |
| core region | central region | 1.5–7.17 Mb |
| arm regions | terminal regions | 0–1.5 Mb (left), 7.17–8.67 Mb (right) |
| two-component system (TCS) | two-component signal transduction | 初出時のみフルスペル |

## 2. 生物名・株名

| 表記 | 用途 |
|------|------|
| *Streptomyces coelicolor* A3(2) strain M145 | 初出時 |
| *S. coelicolor* M145 | 2回目以降 |
| M145 | 文脈が明確な場合 |

## 3. Key Statistics（数値照合用）

### Methylation systems
| Statistic | Value | Figure |
|-----------|-------|--------|
| GCCGGC sites T1 | 1,289 (83% core) | Fig 1a |
| GCCGGC sites T2 | 407 (82% arm) | Fig 1a |
| GCCGGC sites T3 | 21 (62% arm) | Fig 1a |
| GCCGGC T1–T2 Jaccard overlap | 0.000 | — |
| GCCGGC geographic chi-squared | 597 (p ~ 10^-130) | — |
| AAGCCCG sites T1 | 260 (69% core) | — |
| AAGCCCG sites T2 | 64 (64% core) | — |
| SC_RS17645 MTase LFC T2vsT1 | -2.19 | — |
| AAGCCCG de-repression r | 0.003 (p = 0.91) | Fig 4b |

### Protection zone (Layer 2)
| Statistic | Value | Figure |
|-----------|-------|--------|
| Total regulatory genes | 1,055 | — |
| Shielded regulators | 993 (94.1%) | Fig 1c |
| Exposed regulators | 62 (5.9%) | Fig 1c |
| Protection zone width | 2,200 bp (-1,300 to +700) | Fig 1b, 2a |
| Deepest depletion | +300 bp (ratio = 0.541) | Fig 1b |
| TSS methylation ratio (exp/shi) | 8.4x | — |
| Distance boundary threshold | 293 bp | Fig 1c, 2b |
| Distance AUC | 0.917 (95% CI: 0.895–0.935) | Fig 2b |
| baseMean AUC | 0.547 | Fig 2c, 2e |
| Expression quintile JT p | 0.730 | Fig 2e |
| Sequence contribution (CV AUC) | 0.712 (5-fold CV) | Fig 2c |
| Sequence contribution estimate | ~33% | Fig 2c |
| Protein occupancy contribution | ~67% | Fig 2c |
| GCCGGC fold at TFBS | 1.157 (counter-intuitive enrichment) | Fig 2d |
| Mann-Whitney (exposed vs shielded) | p = 5.3 x 10^-28 | Fig 1c |

### 62 Exposed TFs (Layer 3)
| Statistic | Value | Figure |
|-----------|-------|--------|
| Activation bloc size | ~35 genes | Fig 1d, 3 |
| Repression bloc size | ~26 genes | Fig 1d, 3 |
| Unassigned | 1 gene | — |
| TetR enrichment in repression | OR = 0.16, p = 0.028 | Fig 3d |
| TCS pairs with asymmetric split | 7/7 (100%) | Fig 3e |
| Early responders (phase ratio > 0.6) | 39/62 (63%) | Fig 3c |
| Bloc phase separation p | 0.459 (simultaneous) | Fig 3c |
| Within-type co-expression rho | 0.500 (p = 6.6 x 10^-25) | Fig 3a |
| Between-type co-expression rho | -0.283 | — |
| Module-internal coherence rho | 0.717 (p = 1 x 10^-19) | — |
| Methylation–expression timing rho | 0.136 (p = 0.299, NS) | — |
| FIMO motifs for exposed TFs | 0/62 | — |
| Neighborhood permutation p | 0.857 | Fig 4d |
| Neighborhood distance decay rho | +0.026 | Fig 4d |
| Methylated AAGCCCG at exposed TSS | 2/62 (3.2%) | Fig 4c |
| AAGCCCG sequence at exposed TSS | ~16/62 (25.8%) | Fig 4c |

### Simpson's Paradox
| Statistic | Value | Figure |
|-----------|-------|--------|
| GCCGGC de-repression unstratified | p = 8.3 x 10^-8 | Fig 4a |
| GCCGGC de-repression core-only | p = 0.87 | Fig 4a |
| BGC enrichment unstratified | fold = 1.66, p = 4.4 x 10^-8 | — |
| BGC enrichment core-only | fold = 1.07 (NS) | — |

### Activation bloc key genes
| Gene | Old locus | LFC (T3vsT1) | Function |
|------|-----------|--------------|----------|
| RamR | SCO6685 | +7.9 | Aerial mycelium |
| NsdB | — | +7.8 | Negative regulator of development |
| SCO1160 | SCO1160 | +5.6 | Sensor kinase |

### Conservation (H36)
| Statistic | Value |
|-----------|-------|
| Exposed SCO locus rate | 96.8% |
| Composite conservation p | 0.276 (NS) |
| Conservation AUC | 0.459 (below random) |

## 4. Figure / Table 番号

### Main Figures
| # | Short title | Script |
|---|-------------|--------|
| Fig 1 | Gatekeeper Model overview | 10_new_figure1_overview.py |
| Fig 2 | Protection zone characterization | 11_new_figure2_protection.py |
| Fig 3 | 62 exposed TFs — developmental switch | 12_new_figure3_exposed_TFs.py |
| Fig 4 | Negative results | 13_new_figure4_negative_results.py |
| Fig 5 | Vegetative-to-developmental switch model | 14_new_figure5_switch_model.py |

### Supplementary Figures (予定)
| # | Content |
|---|---------|
| S1 | Methylation landscape overview (旧Fig1) |
| S2 | 4mC reclassification + AAGCCCG (旧Fig2+3) |
| S3 | Simpson's paradox 6-panel detail |
| S4 | Regulatory avoidance forest plot (CMH) |
| S5 | Sequence-level motif depletion |
| S6 | 62 exposed TF complete annotation |
| S7 | TCS pair detailed analysis |
| S8 | Conservation metrics (H36) |

### Supplementary Tables (予定)
| # | Content |
|---|---------|
| ST1 | All HC methylation sites |
| ST2 | Motif summary + REBASE conservation |
| ST3 | MTase genes (22 genes) |
| ST4 | 1,055 regulatory genes: shielded/exposed classification |
| ST5 | 62 exposed TFs: full annotation |
| ST6 | TCS pair analysis |
| ST7 | Simpson's Paradox statistics |
| ST8 | Hypothesis ledger (H1–H36) |

## 5. 略語一覧 (Abbreviations)

| Abbreviation | Full form | 初出で定義 |
|-------------|-----------|-----------|
| 4mC | N4-methylcytosine | Yes |
| 6mA | N6-methyladenine | Yes |
| R-M | restriction-modification | Yes |
| TSS | transcription start site | Yes |
| TCS | two-component system | Yes |
| TF | transcription factor | Yes |
| BS | binding site | Yes |
| LFC | log2 fold change | Yes |
| AUC | area under the ROC curve | Yes |
| ROC | receiver operating characteristic | Yes |
| DEG | differentially expressed gene | Yes |
| BGC | biosynthetic gene cluster | Yes |
| FIMO | Find Individual Motif Occurrences | Yes |
| CMH | Cochran-Mantel-Haenszel (test) | Yes |
| JT | Jonckheere-Terpstra (test) | Yes |
| MW | Mann-Whitney U (test) | Yes |
| CV | cross-validation | Yes |
| SK | sensor kinase | contextual |
| RR | response regulator | contextual |

## 6. 書式ルール

- p値: p < 0.001 は科学表記 (例: p = 5.3 x 10^-28)、p >= 0.001 は小数 (例: p = 0.028)
- 統計量: U, rho, r は斜体
- 遺伝子名: 斜体 (*ramR*, *nsdB*)、タンパク質名は立体 (RamR, NsdB)
- 種名: 斜体 (*Streptomyces coelicolor*)
- 数値の有効桁: AUC は小数3桁、p値は状況に応じて、LFCは小数1–2桁
- 染色体位置: Mb単位 (例: 1.5 Mb)
- Figure参照: "(Fig. 1a)" 形式、Supplementaryは "(Fig. S1)"
- bp, kb, Mb: 単位は半角スペースの後 (例: 293 bp, 2.2 kb)

## 7. データ参照パス

```
BASE = /Users/okaban/bioinfo/rna-seq
EPIGENOME = BASE/11_epigenome_integration/analysis
METHYL = /Users/okaban/bioinfo/methyl/260102_M145
FIGURES = BASE/15_paper_figures/figures/main
REPORT = BASE/reports/260227_gatekeeper_model_v4_final_synthesis_report.md
```
