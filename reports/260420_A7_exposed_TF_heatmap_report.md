# A-7: Exposed Regulatory Genes Expression Heatmap — Report

**Date**: 2026-04-20  
**Analysis dir**: `11_epigenome_integration/analysis/71_exposed_tf_heatmap/`  
**Figure**: Figure 5(b) candidate

---

## Overview

57個のExposed regulatory genesの発現をVST-normalized counts(log₂)でZ-score変換し、Wardクラスタリングによるheatmapを作成した。

## Results

### T2 DEG数

| Metric | Value |
|--------|-------|
| 解析対象遺伝子数 | **57** (旧カウント62は更新前) |
| T2でup-regulated (padj<0.05, LFC>1) | **23 / 57 (40.4%)** |

### Wardクラスター構造 (k=3)

| Cluster | n | T2 DEGs | 発現パターン | 解釈 |
|---------|---|---------|------------|------|
| C1 | 23 | 0 | T1高 → T2/T3低下 | T1特異的発現 (early downregulated) |
| C2 | 9 | 0 | T3高 | Late responders |
| C3 | 25 | **23** | T2/T3高 | T2活性化クラスター（メチル化応答候補） |

**注目**: 23個のT2 DEGsは全てCluster 3に集中。Gatekeeper Model v4 Layer 3の「T2でのexposed TF協調的upregulation」を直接支持。

### Top T2-upregulated遺伝子 (Cluster 3)

| Gene | SCO | Product | LFC T2vsT1 | padj |
|------|-----|---------|-----------|------|
| — | SCO7252 | DNA-binding protein NsdB | 9.53 | 1.9×10⁻²⁵² |
| ramR | SCO6685 | TCS response regulator RamR | 8.85 | 3.3×10⁻¹¹ |
| — | SCO1160 | Sensor histidine kinase | 7.90 | 6.7×10⁻⁸⁰ |
| — | SCO1564 | σ70 sigma factor | 5.53 | 3.4×10⁻²⁴ |
| — | SCO1227 | HTH transcriptional regulator | 4.86 | 1.1×10⁻³⁹ |

## Output Files

- **Figure**: `11_epigenome_integration/analysis/71_exposed_tf_heatmap/figures/A7_exposed_TF_heatmap.png` (300 dpi)
- **Table**: `11_epigenome_integration/analysis/71_exposed_tf_heatmap/tables/A7_exposed_TF_expression.tsv`

## Methods Note

- 発現データ: DESeq2 size-factor normalized counts (log₂+1変換後z-score)
- VSTカウントは利用不可のためnormalized countsを使用
- 遺伝子リスト: `51_exposed_regulators_characteristics/tables/exposed_regulators_full_table.tsv` (57遺伝子)
- クラスタリング: Ward linkage, Euclidean distance, 列順はT1/T2/T3固定
- T2 DEGアノテーション: 左側赤バーで表示 (padj<0.05 & LFC>1)
