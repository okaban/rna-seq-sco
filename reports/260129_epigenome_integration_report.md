# Epigenome-Transcriptome Integration Analysis Report

**Project:** *Streptomyces coelicolor* A3(2) M145
**Date:** 2026-01-29
**Analysis Directory:** `/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/`

---

## 1. Executive Summary

This report presents the integrated analysis of epigenome (DNA methylation: 6mA, 4mC) and transcriptome (RNA-seq) data from *Streptomyces coelicolor* A3(2) M145 across three timepoints (T1, T2, T3).

### Key Findings

1. **844 genes** have methylation sites in their promoter regions (-300 to +50 bp from TSS)
2. **Positive correlation** between promoter methylation changes and gene expression changes at T2 vs T1 (Spearman r = 0.17, p < 0.002 for both 6mA and 4mC)
3. **898 gene-modification combinations** show coordinated methylation and expression changes
4. **Data quality issue:** Sample 3-2 (T3) has significantly lower sequencing coverage (14.3x vs 27-52x), affecting T3 methylation detection

---

## 2. Data Overview

### 2.1 Epigenome Data (Nanopore Methylation)

| Metric | Value |
|--------|-------|
| Total high-confidence methylation sites | 6,028 |
| 6mA sites | 2,734 |
| 4mC sites | 3,294 |

**Distribution by Timepoint:**

| Timepoint | 6mA Sites | 4mC Sites | Total |
|-----------|-----------|-----------|-------|
| T1 | 983 | 1,366 | 2,349 |
| T2 | 1,093 | 1,701 | 2,794 |
| T3 | 658 | 227 | 885 |

**Note:** The low number of T3 sites is due to low coverage in sample 3-2 (14.3x vs 27-52x average).

### 2.2 Transcriptome Data (RNA-seq DESeq2)

| Comparison | Upregulated (padj<0.05) | Downregulated (padj<0.05) |
|------------|------------------------|---------------------------|
| T2 vs T1 | 2,563 | 2,697 |
| T3 vs T1 | 3,290 | 2,898 |
| T3 vs T2 | 2,639 | 2,904 |

Total genes analyzed: 7,646

---

## 3. Epigenome Data Quality Assessment

### 3.1 Sequencing Coverage by Sample

| Sample | Timepoint | Avg Coverage | Status |
|--------|-----------|--------------|--------|
| 1-1 | T1 | 37.1x | OK |
| 1-2 | T1 | 34.7x | OK |
| 1-3 | T1 | 51.7x | OK |
| 2-1 | T2 | 41.5x | OK |
| 2-3 | T2 | 37.6x | OK |
| 2-4 | T2 | 36.1x | OK |
| **3-2** | **T3** | **14.3x** | **LOW** |
| 3-3 | T3 | 33.3x | OK |
| 3-4 | T3 | 27.2x | OK |

### 3.2 Quality Assessment Conclusion

- **Pileup files:** All 9 samples processed correctly (1.2-1.3 GB each)
- **Issue:** Sample 3-2 has ~40% of average coverage, likely due to lower sequencing depth
- **Impact:** Reduced sensitivity for T3 methylation detection
- **Recommendation:** Consider re-sequencing sample 3-2 or use 2-replicate analysis for T3

---

## 4. Integration Analysis Results

### 4.1 Promoter Methylation Summary

- **Promoter region defined:** -300 to +50 bp from TSS
- **Total genes in genome:** 8,083
- **Genes with promoter methylation:** 844 (10.4%)

| Modification | T1 | T2 | T3 |
|-------------|-----|-----|-----|
| 6mA | 261 genes | 271 genes | 194 genes |
| 4mC | 347 genes | 401 genes | 52 genes |

### 4.2 Correlation Analysis: Methylation vs Expression

| Modification | Comparison | N genes | Spearman r | P-value | Interpretation |
|-------------|------------|---------|------------|---------|----------------|
| 6mA | T2 vs T1 | 361 | **0.169** | **0.001** | Significant positive |
| 6mA | T3 vs T1 | 330 | -0.050 | 0.367 | Not significant |
| 6mA | T3 vs T2 | 343 | **-0.166** | **0.002** | Significant negative |
| 4mC | T2 vs T1 | 474 | **0.172** | **0.0002** | Significant positive |
| 4mC | T3 vs T1 | 346 | -0.041 | 0.445 | Not significant |
| 4mC | T3 vs T2 | 396 | -0.093 | 0.063 | Marginal |

### 4.3 Key Interpretation

1. **T1 → T2 transition:** Both 6mA and 4mC promoter methylation changes show **positive correlation** with gene expression changes (r ≈ 0.17, p < 0.002)
   - Genes gaining methylation tend to increase expression
   - Genes losing methylation tend to decrease expression

2. **T2 → T3 transition:** Weak negative correlation observed for 6mA (r = -0.17, p = 0.002)
   - This may indicate a regulatory shift or be affected by low T3 coverage

3. **Overall pattern:** DNA methylation (both 6mA and 4mC) appears to have a **positive regulatory role** in early growth phase (T1→T2), which may shift in later phase

---

## 5. Genes with Coordinated Methylation and Expression Changes

Identified **898 gene-modification combinations** where:
- Differential expression: |log2FC| > 1 and padj < 0.05
- Methylation change: >10% change in modification frequency

These candidates are saved in: `significant_coordinated_changes.csv`

---

## 6. Biological Implications

### 6.1 DNA Methylation in Streptomyces

- **6mA (N6-methyladenine):** Common in bacteria, often associated with DNA replication and gene regulation
- **4mC (N4-methylcytosine):** Restriction-modification system component, may have regulatory roles

### 6.2 Temporal Dynamics

The observed pattern suggests:
1. **Active methylation-mediated regulation** during T1→T2 transition (possibly exponential to stationary phase)
2. **Different regulatory dynamics** at T3 (possibly related to secondary metabolism onset)

### 6.3 Limitations

- T3 data quality is compromised (sample 3-2 low coverage)
- 5mC was not detected in high-confidence sites (may require different detection parameters)
- Correlation ≠ causation; functional validation needed

---

## 7. Output Files

| File | Description |
|------|-------------|
| `integrated_methyl_expression.csv` | Full integration table (8,083 genes × methylation + expression) |
| `correlation_analysis.csv` | Spearman correlation results |
| `significant_coordinated_changes.csv` | Genes with coordinated changes |
| `analysis_summary.txt` | Summary statistics |

---

## 8. Methods

### 8.1 Methylation Data Processing
- **Source:** Nanopore sequencing (modkit pileup output)
- **Filtering:** High-confidence sites with ≥3 replicates, ≥10x coverage, ≥50% modification frequency

### 8.2 Transcriptome Data Processing
- **Alignment:** STAR to GCF_000203835.1
- **Quantification:** featureCounts
- **Differential expression:** DESeq2 (padj < 0.05, |log2FC| > 1)

### 8.3 Integration
- **Promoter definition:** -300 to +50 bp from TSS
- **Correlation:** Spearman rank correlation
- **Machine specs:** 14 cores, 48 GB RAM (macOS)

---

## 9. Recommendations

1. **Re-sequence sample 3-2** or perform sensitivity analysis excluding this sample
2. **Functional validation** of top candidate genes with coordinated changes
3. **Motif analysis** around methylated promoters to identify potential regulatory elements
4. **ChIP-seq integration** if available for transcription factors

---

*Report generated by Claude Code integration analysis pipeline*
