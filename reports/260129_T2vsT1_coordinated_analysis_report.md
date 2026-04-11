# Promoter DNA Methylation and Gene Expression Correlation Analysis

## *Streptomyces coelicolor* A3(2) M145: T2 vs T1 Comparison

**Date:** 2026-01-29
**Analysis Method:** Original (3-replicate, equal weight)
**Focus:** T2 vs T1 transition only

---

## Abstract

This analysis examines the relationship between promoter DNA methylation (6mA and 4mC) and gene expression changes during the T1 to T2 growth phase transition in *Streptomyces coelicolor* A3(2) M145. We identified 558 gene-modification pairs showing coordinated methylation and expression changes, with a significant positive correlation between promoter methylation dynamics and transcriptional regulation. Notably, multiple biosynthetic gene clusters (BGCs) including Act, Red, CDA, and Cpk show coordinated methylation-expression changes, suggesting epigenetic regulation of secondary metabolism.

---

## 1. Introduction

DNA methylation in bacteria, particularly N6-methyladenine (6mA) and N4-methylcytosine (4mC), plays crucial roles in restriction-modification systems and gene regulation. Unlike eukaryotic 5-methylcytosine (5mC), these modifications are prevalent in actinobacteria and may contribute to transcriptional control.

This study integrates nanopore-based methylation profiling with RNA-seq differential expression analysis to characterize the methylation-expression relationship during growth phase transition.

---

## 2. Materials and Methods

### 2.1 Data Sources

| Data Type | Platform | Samples | Reference |
|-----------|----------|---------|-----------|
| Methylation | Oxford Nanopore | 9 (3×3 timepoints) | GCF_000203835.1 |
| Transcriptome | Illumina RNA-seq | 9 (3×3 timepoints) | GCF_000203835.1 |

### 2.2 Analysis Parameters

- **Promoter region:** -300 to +50 bp from transcription start site (TSS)
- **High-confidence methylation:** ≥10× coverage, ≥50% modification frequency, ≥2/3 replicates
- **Differential expression:** DESeq2, padj < 0.05, |log₂FC| > 0.5
- **Methylation change categories:**
  - Gained: 0 sites in T1, ≥1 site in T2
  - Lost: ≥1 site in T1, 0 sites in T2
  - Increased: >10% frequency increase
  - Decreased: >10% frequency decrease
  - Stable: <10% change

### 2.3 Statistical Analysis

- Spearman rank correlation for methylation-expression relationship
- Mann-Whitney U test for group comparisons
- All analyses performed in Python 3.x with scipy, pandas, numpy

---

## 3. Results

### 3.1 Overview of Coordinated Changes

| Metric | Count |
|--------|-------|
| Total gene-modification pairs analyzed | 558 |
| Genes with 6mA promoter methylation | 230 |
| Genes with 4mC promoter methylation | 328 |

### 3.2 Correlation Analysis

#### 3.2.1 6mA Promoter Methylation

| Statistic | Value |
|-----------|-------|
| Spearman ρ | 0.169 |
| p-value | 0.0012 |
| n (genes) | 361 |

**Interpretation:** Significant positive correlation (p < 0.01) between 6mA methylation change and gene expression change.

#### 3.2.2 4mC Promoter Methylation

| Statistic | Value |
|-----------|-------|
| Spearman ρ | 0.172 |
| p-value | 0.0002 |
| n (genes) | 474 |

**Interpretation:** Highly significant positive correlation (p < 0.001) between 4mC methylation change and gene expression change.

### 3.3 Correlation Direction Summary

| Category | Count | Percentage |
|----------|-------|------------|
| **Positive correlation** | 153 | 27.4% |
| (Methyl ↑ & Expr ↑) or (Methyl ↓ & Expr ↓) | | |
| **Negative correlation** | 89 | 16.0% |
| (Methyl ↑ & Expr ↓) or (Methyl ↓ & Expr ↑) | | |
| **Other patterns** | 316 | 56.6% |
| (Stable methylation or mild expression) | | |

**Key finding:** Positive correlation genes outnumber negative correlation genes by 1.7:1, supporting a model where promoter methylation positively regulates transcription in *S. coelicolor*.

---

## 4. Biosynthetic Gene Cluster (BGC) Analysis

### 4.1 BGC Summary

A total of **14 genes from 7 different BGCs** showed coordinated methylation-expression changes during the T1→T2 transition.

| BGC | Product | Coordinated Genes | Positive Corr. | Negative Corr. |
|-----|---------|-------------------|----------------|----------------|
| **CDA** | Calcium-Dependent Antibiotic | 4 | 0 | 1 |
| **Act** | Actinorhodin (blue pigment) | 2 | 0 | 1 |
| **Desferrioxamine** | Siderophore | 2 | 0 | 0 |
| **Coelichelin** | Siderophore (NRPS) | 2 | 1 | 0 |
| **SapB** | Lanthipeptide (morphogenetic) | 2 | 1 | 1 |
| **Red** | Undecylprodigiosin (red pigment) | 1 | 1 | 0 |
| **Cpk** | Coelimycin P1 (cryptic PKS) | 1 | 1 | 0 |

### 4.2 Act Cluster (Actinorhodin) Analysis

**Actinorhodin** is the blue-pigmented polyketide antibiotic, a hallmark of *S. coelicolor* secondary metabolism.

| Gene | Old Locus | Modification | Δ Methylation | log₂FC | Pattern |
|------|-----------|--------------|---------------|--------|---------|
| SC_RS27535 | **SCO5075** | 4mC | +7.4% (Stable) | **+1.05** | Stable_Up |
| SC_RS27555 | **SCO5079** | 4mC | -70.8% (Lost) | **+1.44** | Lost_Up |

**Interpretation:** Both Act cluster genes show **expression upregulation** at T2. SCO5079 (NmrA/HSCARG family protein) loses 4mC methylation while being upregulated, suggesting possible **methylation-mediated repression release**.

### 4.3 Red Cluster (Undecylprodigiosin) Analysis

**Undecylprodigiosin** is the red-pigmented prodiginine antibiotic.

| Gene | Old Locus | Modification | Δ Methylation | log₂FC | Pattern |
|------|-----------|--------------|---------------|--------|---------|
| SC_RS31730 | **SCO5897** | 4mC | **+71.5% (Gained)** | **+2.60** | **Gained_Up** |

**Interpretation:** SCO5897 (aromatic ring-hydroxylating oxygenase) shows **strong positive correlation** - gains 4mC methylation and is strongly upregulated (6× increase). This supports a model of **methylation-activated transcription** in the Red cluster.

### 4.4 CDA Cluster (Calcium-Dependent Antibiotic) Analysis

**CDA** is a lipopeptide antibiotic produced by NRPS machinery.

| Gene | Old Locus | Function | Modification | Δ Methylation | log₂FC | Pattern |
|------|-----------|----------|--------------|---------------|--------|---------|
| SC_RS18170 | SCO3211 | TrpC | 4mC | -0.6% (Stable) | **+8.09** | Stable_Up |
| SC_RS18280 | SCO3233 | Hydrolase | 4mC | -13.7% (Decreased) | **+6.59** | Decreased_Up |
| SC_RS18320 | SCO3241 | Epimerase | 6mA | +7.6% (Stable) | **+8.93** | Stable_Up |
| SC_RS18355 | SCO3248 | KS (β-ketoacyl synthase) | 4mC | +2.5% (Stable) | **+8.79** | Stable_Up |

**Interpretation:** The CDA cluster shows **massive transcriptional activation** (up to 500× increase for SCO3241). Most genes maintain stable methylation, but SCO3233 shows decreased 4mC with upregulation - a **negative correlation** suggesting possible derepression mechanism.

### 4.5 Cpk Cluster (Coelimycin P1) Analysis

**Coelimycin P1** is the cryptic type I polyketide, recently characterized.

| Gene | Old Locus | Modification | Δ Methylation | log₂FC | Pattern |
|------|-----------|--------------|---------------|--------|---------|
| SC_RS33670 | **SCO6284** | 6mA | **+64.4% (Gained)** | **+7.36** | **Gained_Up** |

**Interpretation:** SCO6284 (acyl-CoA carboxylase) is **strongly upregulated with 6mA gain**. This enzyme provides malonyl-CoA for polyketide biosynthesis, suggesting **coordinated metabolic and epigenetic regulation** of secondary metabolism.

### 4.6 SapB/Ram Cluster Analysis

**SapB** is a lanthipeptide involved in aerial mycelium formation.

| Gene | Old Locus | Function | Modification | Δ Methylation | log₂FC | Pattern |
|------|-----------|----------|--------------|---------------|--------|---------|
| SC_RS35605 | SCO6684 | ABC transporter | 6mA | **+58.2% (Gained)** | **+5.65** | Gained_Up |
| SC_RS35610 | **SCO6685 (ramR)** | Response regulator | 6mA | **-64.4% (Lost)** | **+6.60** | Lost_Up |

**Key Finding:** **RamR (SCO6685)** is a critical pathway-specific activator of the SapB cluster. It shows a **negative correlation pattern** - loses 6mA methylation while being strongly upregulated. This suggests that **6mA methylation may repress ramR**, and demethylation during T1→T2 transition leads to its activation.

### 4.7 Siderophore Clusters Analysis

#### Desferrioxamine Cluster

| Gene | Old Locus | Modification | Δ Methylation | log₂FC | Pattern |
|------|-----------|--------------|---------------|--------|---------|
| SC_RS16000 | SCO2783 | 4mC | -4.2% (Stable) | **+6.41** | Stable_Up |
| SC_RS16005 | SCO2784 | 6mA | -4.0% (Stable) | **+5.93** | Stable_Up |

**Interpretation:** Both desferrioxamine biosynthesis genes are **strongly upregulated** (30-80× increase) with stable methylation, indicating transcriptional activation **independent of methylation changes**.

---

## 5. Discussion

### 5.1 Positive Correlation Dominance

The predominance of positive correlation (153 vs 89 genes) suggests that in *S. coelicolor*, promoter DNA methylation may act as a **positive regulatory mark** rather than a repressive one. This contrasts with the canonical view of DNA methylation as a silencing mechanism in eukaryotes.

### 5.2 BGC-Specific Methylation Patterns

Our analysis reveals **cluster-specific epigenetic regulation**:

1. **Positive correlation (methylation activates):**
   - Red: SCO5897 (4mC gain → upregulation)
   - Cpk: SCO6284 (6mA gain → upregulation)
   - SapB: SCO6684 (6mA gain → upregulation)

2. **Negative correlation (methylation represses):**
   - Act: SCO5079 (4mC loss → upregulation)
   - SapB/RamR: SCO6685 (6mA loss → upregulation)
   - CDA: SCO3233 (4mC decrease → upregulation)

3. **Expression-dominant (stable methylation):**
   - CDA: Most genes (massive upregulation, stable methylation)
   - Desferrioxamine: Both genes (strong upregulation, stable methylation)

### 5.3 RamR as a Key Epigenetically Regulated Target

The discovery that **ramR (SCO6685)** shows coordinated demethylation and activation is particularly significant:

- RamR is the response regulator controlling SapB production
- SapB is essential for aerial mycelium formation
- The T1→T2 transition involves morphological differentiation
- **6mA demethylation of ramR may be a switch for developmental transition**

### 5.4 Secondary Metabolism Activation

The T1→T2 transition shows widespread **BGC activation**:
- 7 of 13 analyzed BGCs have at least one gene with coordinated changes
- All BGC genes with coordinated changes show **upregulation** at T2
- This is consistent with the transition from vegetative growth (T1) to secondary metabolism phase (T2)

### 5.5 Biological Model

Based on our findings, we propose:

```
T1 (Vegetative Growth)
    ↓
    • Specific methylation patterns maintain BGC genes in low/basal expression
    ↓
T1→T2 Transition
    ↓
    • Coordinated methylation changes at key regulatory genes (e.g., ramR)
    • Both methylation gain (activation) and loss (derepression) contribute
    ↓
T2 (Secondary Metabolism)
    ↓
    • BGC genes activated
    • Antibiotic production begins
    • Morphological differentiation initiates
```

---

## 6. Conclusions

1. **Significant positive correlation** exists between promoter methylation (6mA, 4mC) and gene expression during the T1→T2 transition

2. **14 BGC genes from 7 clusters** show coordinated methylation-expression changes

3. **Act and Red clusters** both contain genes with coordinated changes, supporting epigenetic regulation of pigmented antibiotic biosynthesis

4. **RamR (SCO6685)**, a key developmental regulator, shows strong demethylation-activation correlation, suggesting it as a potential **epigenetic switch** for differentiation

5. The **CDA and Cpk clusters** show strong transcriptional activation with varying methylation patterns

6. Both **positive (methylation-activation)** and **negative (methylation-repression)** regulatory mechanisms operate in *S. coelicolor* BGCs

---

## 7. Figures

| Figure | Description | File |
|--------|-------------|------|
| Fig. 1 | Methylation vs expression scatter plot | `Fig1_methylation_expression_correlation.pdf` |
| Fig. 2 | Coordination patterns by category | `Fig2_coordination_patterns.pdf` |
| Fig. 3 | Volcano-style methylation plot | `Fig3_volcano_methylation.pdf` |
| Fig. 4 | Summary statistics panel | `Fig4_summary_statistics.pdf` |
| Fig. 5 | Top genes heatmap | `Fig5_top_genes_heatmap.pdf` |

---

## 8. Data Availability

| File | Description |
|------|-------------|
| `T2vsT1_coordinated_genes.csv` | All 558 gene-modification pairs |
| `T2vsT1_positive_correlation.csv` | 153 positive correlation genes |
| `T2vsT1_negative_correlation.csv` | 89 negative correlation genes |
| `BGC_summary.csv` | BGC analysis summary |
| `BGC_detailed.csv` | Detailed BGC gene information |

---

## 9. References

1. Bentley SD et al. (2002) Complete genome sequence of the model actinomycete *Streptomyces coelicolor* A3(2). Nature 417:141-147.
2. Flärdh K, Buttner MJ (2009) Streptomyces morphogenetics: dissecting differentiation in a filamentous bacterium. Nat Rev Microbiol 7:36-49.
3. van Wezel GP, McDowall KJ (2011) The regulation of the secondary metabolism of Streptomyces: new links and experimental advances. Nat Prod Rep 28:1311-1333.

---

*Report generated by integrated epigenome-transcriptome analysis pipeline*
*Streptomyces coelicolor A3(2) M145 Project*
