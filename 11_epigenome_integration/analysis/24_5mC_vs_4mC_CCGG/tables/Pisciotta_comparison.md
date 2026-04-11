# Pisciotta et al. (2023) vs 本研究 — 直接比較

## 論文情報

**Pisciotta et al.** "The DNA cytosine methylome revealed two methylation motifs in the upstream regions of genes related to morphological and physiological differentiation in *Streptomyces coelicolor* A(3)2 M145"
*Scientific Reports* 13, 7038 (2023). DOI: 10.1038/s41598-023-34075-1

## 方法論比較

| 項目 | Pisciotta et al. (2023) | 本研究 |
|------|------------------------|--------|
| 手法 | Whole-genome bisulfite sequencing (WGBS) | Oxford Nanopore (R10.4.1) |
| 検出修飾 | 「5mC」（※BS-seqは5mC/4mCを区別不能） | 6mA, 4mC, 5mC（独立検出） |
| Basecallモデル | N/A | dorado SUP v5.2.0 + **4mC_5mC@v1** |
| 5mC/4mC区別 | **不可能** | **可能** |
| 6mA検出 | **不可能** | **可能** |
| 株 | M145 | M145 |
| 培地 | MG液体 | 固体培地 |
| 総メチル化サイト | 3,360 ("5mC") | 12,429 (4mC: 6,080 + 6mA: 6,349) |

## モチーフ対応

### モチーフ1: GGCCGG

| 項目 | Pisciotta | 本研究 |
|------|-----------|--------|
| 報告モチーフ | GGC**m**CGG | GGCCGG / CCGG |
| 修飾タイプ | 「5mC」 | **4mC** |
| 5mC検出 | BS-seq陽性 | Nanopore **ゼロ** (9サンプル全て) |
| 4mC検出 | BS-seqでは区別不能 | HC sites: 998-1,533/サンプル |
| R-M系 | 未同定 | 既知（汎細菌的CCGG R-M） |

**結論**: BS-seqが報告した「5mC」はN4-methylcytosine (4mC)の誤分類

### モチーフ2: AAGCCCG

| 項目 | Pisciotta | 本研究 |
|------|-----------|--------|
| 報告モチーフ | GCC**m**CG / AAGCC**m**CG | **AAGCCCG** |
| 修飾タイプ | 「5mC」 | **4mC + 6mA（二重修飾）** |
| 4mC成分 | BS-seqでは5mCと区別不能 | 973 sites (36.3%) |
| 6mA成分 | **BS-seqでは原理的に検出不能** | 360 sites (11.2%) |
| REBASE | 未確認 | **0/82種 — 属レベルで新規** |
| 候補MTase | SCO1731（Kim et al. 2018） | **SC_RS17645 (SCO3104)** — HsdM型 |

**結論**:
1. BS-seqが報告した「5mC」は4mCの誤分類
2. さらに、BS-seqでは6mA成分が完全に見落とされていた
3. このモチーフはREBASE未登録の**完全に新規なR-Mシステム**

## BS-seqの根本的限界（本研究で実証）

```
BS-seq:
  未メチル化 C → U（変換）
  5mC → C（変換されない）
  4mC → C（変換されない）  ← 5mCと区別不能！
  6mA → 検出不能

Nanopore (4mC_5mC@v1):
  5mC → 独立検出（本データ: ゼロ）
  4mC → 独立検出（本データ: 2,693 unique sites）
  6mA → 独立検出（本データ: 3,248 unique sites）
```

## 本研究からの主張

1. *S. coelicolor* M145のCCGG/GGCCGGサイトにおけるシトシン修飾は**全て4mC**であり、5mCは検出されない
2. Pisciotta et al.のGCCCG/AAGCCCGモチーフは**4mCの誤分類**に加え、**6mA成分が完全に見落とされていた**
3. これはBS-seqが原理的に4mC/5mCを区別できないことに起因する方法論的限界である
4. 細菌メチロームにおいてBS-seqの結果は慎重に解釈すべきであり、Nanopore等の直接検出法による検証が必要
