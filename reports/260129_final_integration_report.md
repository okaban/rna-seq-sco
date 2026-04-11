# Final Epigenome-Transcriptome Integration Report

**Project:** *Streptomyces coelicolor* A3(2) M145
**Date:** 2026-01-29
**Analysis Directory:** `/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/`

---

## Executive Summary

This report presents the comprehensive integration of epigenome (DNA methylation) and transcriptome (RNA-seq) data, including three different analytical approaches and detailed 5mC verification.

### Key Conclusions

1. **5mC is virtually absent** in *S. coelicolor* - verified from raw data
2. **6mA and 4mC are the primary modifications** (restriction-modification systems)
3. **4mC shows consistent positive correlation** with gene expression at T1→T2
4. **Sample 3-2 quality issue** significantly impacts T3 analysis
5. **Coverage-weighted analysis** is the recommended approach for n=3

---

## 1. 5mC Verification (Detailed)

### 1.1 Raw Data Evidence

**BAM file methylation tags confirmed:**
```
MM:Z:A+a  → 6mA (N6-methyladenine)
MM:Z:C+21839 → 4mC (N4-methylcytosine)
MM:Z:C+m  → 5mC (5-methylcytosine)
```

All three modification types are being detected by the basecaller.

### 1.2 5mC Frequency Distribution

| Sample | Total 5mC sites | cov≥10x | freq≥10% | freq≥20% | Max freq |
|--------|----------------|---------|----------|----------|----------|
| 1-1 | 6,464,535 | 6,202,115 | 10 | 1 | 21.4% |
| 1-2 | 6,460,976 | 6,199,912 | 3 | 0 | 18.2% |
| 1-3 | 6,554,472 | 6,242,521 | 4 | 0 | 16.7% |
| 2-1 | 6,520,944 | 6,248,242 | 2 | 1 | 21.1% |
| 2-3 | 6,511,379 | 6,247,117 | 0 | 0 | 8.3% |
| 2-4 | 6,495,936 | 6,246,165 | 1 | 0 | 10.0% |
| 3-2 | 6,358,598 | 5,131,803 | 60* | 0 | 10.0% |
| 3-3 | 6,483,683 | 6,243,556 | 1 | 0 | 10.0% |
| 3-4 | 6,422,161 | 6,226,873 | 3 | 0 | 10.5% |

*Sample 3-2 shows elevated false-positive rate due to low coverage

### 1.3 5mC Conclusion

**5mC is NOT a significant modification in *S. coelicolor*.**

- Mean modification frequency: <0.01%
- Maximum observed: 21.4% (single site)
- No sites pass high-confidence threshold (≥50%)
- **Biologically expected:** 5mC is rare in bacteria; 6mA and 4mC dominate

---

## 2. Sample Quality Assessment

### 2.1 Sequencing Coverage

| Sample | Timepoint | Avg Coverage | Status |
|--------|-----------|--------------|--------|
| 1-1 | T1 | 37.1x | ✓ Good |
| 1-2 | T1 | 34.7x | ✓ Good |
| 1-3 | T1 | 51.7x | ✓ Excellent |
| 2-1 | T2 | 41.5x | ✓ Good |
| 2-3 | T2 | 37.6x | ✓ Good |
| 2-4 | T2 | 36.1x | ✓ Good |
| **3-2** | **T3** | **14.3x** | **⚠️ Low** |
| 3-3 | T3 | 33.3x | ✓ Good |
| 3-4 | T3 | 27.2x | ✓ Acceptable |

### 2.2 Impact on Analysis

Sample 3-2's low coverage (40% of average) causes:
- Reduced T3 methylation site detection
- Potential false negatives
- Biased correlation estimates

---

## 3. Three Analytical Approaches Compared

### 3.1 Method Descriptions

| Method | Description | T3 Sites | Pros | Cons |
|--------|-------------|----------|------|------|
| **Original (3-rep)** | Require detection in all 3 replicates | 885 | Highest confidence | Low-coverage sample bottleneck |
| **2-rep (exclude 3-2)** | Use only 3-3 and 3-4 for T3 | 2,780 | Avoids bias | Loses replicate |
| **Weighted (3-rep)** | Coverage-weighted average, all samples | 3,334 | Uses all data, reduces bias | Lower per-site confidence |

### 3.2 Site Detection Comparison

| Timepoint | Original | 2-rep | Weighted |
|-----------|----------|-------|----------|
| T1 | 2,349 | 2,349 | 3,884 |
| T2 | 2,794 | 2,794 | 4,560 |
| T3 | 885 | 2,780 | 3,334 |

### 3.3 Correlation Analysis Comparison

| Mod | Comparison | r (Original) | p | r (2-rep) | p | r (Weighted) | p |
|-----|------------|--------------|---|-----------|---|--------------|---|
| 6mA | T2_vs_T1 | **0.169** | **0.001** | **0.169** | **0.001** | 0.007 | 0.867 |
| 6mA | T3_vs_T1 | -0.050 | 0.367 | 0.064 | 0.173 | -0.048 | 0.210 |
| 6mA | T3_vs_T2 | **-0.166** | **0.002** | -0.059 | 0.206 | -0.043 | 0.264 |
| 4mC | T2_vs_T1 | **0.172** | **0.0002** | **0.172** | **0.0002** | **0.137** | **0.0006** |
| 4mC | T3_vs_T1 | -0.041 | 0.445 | 0.083 | 0.086 | -0.081 | 0.068 |
| 4mC | T3_vs_T2 | -0.093 | 0.063 | -0.003 | 0.952 | -0.042 | 0.311 |

---

## 4. Key Biological Findings

### 4.1 Consistent Finding Across All Methods

**4mC promoter methylation shows positive correlation with gene expression at T1→T2**

- Original: r = 0.172, p = 0.0002
- 2-rep: r = 0.172, p = 0.0002
- Weighted: r = 0.137, p = 0.0006

This is a robust finding that persists regardless of analytical method.

### 4.2 6mA Findings

The 6mA correlation at T2_vs_T1 is sensitive to analytical approach:
- Original/2-rep: Significant (r = 0.17, p = 0.001)
- Weighted: Not significant (r = 0.007, p = 0.87)

This suggests the 6mA-expression relationship may be driven by a subset of high-confidence sites.

### 4.3 T3 Interpretation

T3 comparisons show no consistent significant correlations across methods, possibly due to:
1. True biological change in regulatory dynamics
2. Residual sample quality effects
3. Transition to secondary metabolism phase

### 4.4 T3 Sample Contribution (Weighted Method)

| Sample | Avg Coverage Contribution |
|--------|---------------------------|
| 3-2 | 21.6% |
| 3-3 | 45.5% |
| 3-4 | 35.1% |

The weighted method appropriately down-weights sample 3-2's contribution.

---

## 5. Recommendations

### 5.1 For This Dataset

1. **Use weighted method** for comprehensive analysis (includes all data, reduces bias)
2. **Report 4mC T2_vs_T1 correlation** as the primary finding (robust across methods)
3. **Interpret T3 results with caution** due to sample quality issue
4. **Do not report 5mC** - it is not biologically relevant for this organism

### 5.2 For Publications

Report all three methods with appropriate caveats:
- Note sample 3-2 quality issue
- Emphasize findings consistent across methods
- Discuss potential impact of sample quality on T3 results

### 5.3 For Future Experiments

1. **Consider re-sequencing sample 3-2** if T3 is critical
2. **Set minimum coverage threshold** of 20x per sample
3. **Use coverage-weighted analysis** as default for variable-quality datasets

---

## 6. Output Files

| File | Description |
|------|-------------|
| **Core Analysis** | |
| `integrated_methyl_expression.csv` | Original 3-rep integration |
| `integrated_methyl_expression_T3_2rep.csv` | 2-rep T3 integration |
| `integrated_methyl_expression_weighted.csv` | Weighted integration |
| **Correlation Results** | |
| `correlation_analysis.csv` | Original correlations |
| `correlation_analysis_T3_2rep.csv` | 2-rep correlations |
| `correlation_analysis_weighted.csv` | Weighted correlations |
| **Methylation Sites** | |
| `high_confidence_sites_T3_2rep.csv` | Combined with 2-rep T3 |
| `high_confidence_sites_weighted.csv` | Weighted consensus sites |
| `T3_2rep_methylation_sites.csv` | T3 sites (2-rep only) |
| **Reports** | |
| `INTEGRATION_REPORT.md` | Initial report |
| `T3_2rep_comparison_report.md` | 2-rep comparison |
| `FINAL_INTEGRATION_REPORT.md` | This comprehensive report |

---

## 7. Methods Summary

### 7.1 Methylation Analysis
- **Platform:** Oxford Nanopore (direct DNA sequencing)
- **Basecalling:** With methylation detection (6mA, 4mC, 5mC)
- **Processing:** modkit pileup for per-site quantification
- **Thresholds:** ≥10x coverage, ≥50% modification frequency, ≥2 replicates

### 7.2 Transcriptome Analysis
- **Platform:** Illumina RNA-seq
- **Alignment:** STAR to GCF_000203835.1
- **Quantification:** featureCounts
- **Differential expression:** DESeq2 (padj < 0.05, |log2FC| > 1)

### 7.3 Integration
- **Promoter definition:** -300 to +50 bp from TSS
- **Correlation:** Spearman rank correlation
- **Resources:** 14 cores, 48 GB RAM (macOS)

---

## 8. Conclusion

The integration of epigenome and transcriptome data reveals that **4mC promoter methylation is positively correlated with gene expression changes during the T1→T2 transition** in *S. coelicolor* M145. This finding is robust across all three analytical approaches.

The sample 3-2 quality issue necessitates careful interpretation of T3 results, but the coverage-weighted analysis provides a reasonable compromise between including all data and minimizing bias.

5mC was confirmed to be virtually absent in this organism through detailed verification of raw data, consistent with the biology of actinobacteria which primarily use 6mA and 4mC in their restriction-modification systems.

---

*Report generated by Claude Code integration analysis pipeline*
