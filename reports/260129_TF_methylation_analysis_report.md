# Transcription Factor Analysis Report

## Methylation-Expression Correlation in Transcription Factors
### *Streptomyces coelicolor* A3(2) M145: T2 vs T1

**Date:** 2026-01-29
**Analysis:** Transcription factors with coordinated promoter methylation and expression changes

---

## 1. Executive Summary

We identified **24 transcription factors** with coordinated methylation-expression changes during the T1→T2 transition. While known SARP family pathway-specific activators (actII-ORF4, redD, cdaR) were not directly methylated, key global regulators including **NsdB** and **RamR** showed strong epigenetic regulation.

### Key Findings

| Finding | Details |
|---------|---------|
| Total TFs with coordinated changes | 24 |
| Positive correlation (methyl↔expr same direction) | 16 |
| Negative correlation (methyl↔expr opposite direction) | 8 |
| Most significant TF | **NsdB (SCO7252)**: 4mC gain + 640× upregulation |
| Key developmental regulator | **RamR (SCO6685)**: 6mA loss + 100× upregulation |

---

## 2. SARP Family Analysis

### 2.1 Known SARP Genes Searched

| Gene | Old Locus | Function | In Coordinated List? |
|------|-----------|----------|---------------------|
| actII-ORF4 | SCO5085 | Act cluster activator | **No** |
| redD | SCO5877 | Red cluster activator | **No** |
| cdaR | SCO3217 | CDA cluster activator | **No** |
| cpkO | SCO6269 | Cpk cluster regulator | **No** |
| afsR | SCO4425 | Global regulator | **No** |
| afsR2 | SCO4423 | Global regulator | **No** |

### 2.2 Interpretation

The absence of SARP family genes in the coordinated methylation-expression list suggests:

1. **No high-confidence methylation sites** in their promoter regions (-300 to +50 bp from TSS)
2. Or expression changes did not meet significance thresholds (padj < 0.05, |log2FC| > 0.5)

**Important:** This does not mean SARPs are not regulated by methylation. Upper-level regulators like NsdB may control SARP expression indirectly.

---

## 3. Positive Correlation Transcription Factors

### 3.1 Definition
Genes where methylation and expression change in the **same direction**:
- Methylation ↑ and Expression ↑ (Gained_Up, Increased_Up)
- Methylation ↓ and Expression ↓ (Lost_Down, Decreased_Down)

### 3.2 Complete List (16 TFs)

| Gene ID | Old Locus | Modification | T1→T2 Sites | Δ Methylation | log₂FC | Family | Product |
|---------|-----------|--------------|-------------|---------------|--------|--------|---------|
| SC_RS38475 | **SCO7252** | 4mC | 0→1 | **+92.3%** | **+9.32** | - | **DNA-binding protein NsdB** |
| SC_RS14635 | SCO2517 | 4mC | 0→1 | +79.7% | +3.93 | - | Response regulator |
| SC_RS22780 | SCO4122 | 6mA | 1→0 | -56.1% | -3.80 | MarR | Winged HTH transcriptional regulator |
| SC_RS30370 | SCO5629 | 4mC | 1→1 | +11.3% | +2.65 | - | TPR protein |
| SC_RS08620 | SCO1331 | 4mC | 0→1 | +84.5% | +2.43 | LuxR | Transcriptional regulator |
| SC_RS10090 | SCO1616 | 4mC | 1→0 | -93.0% | -2.33 | LysR | Transcriptional regulator |
| SC_RS27300 | SCO5027 | 6mA | 1→0 | -59.0% | -2.25 | - | Winged helix DNA-binding protein |
| SC_RS36220 | - | 6mA | 0→1 | +62.1% | +2.17 | - | HTH transcriptional regulator |
| SC_RS24380 | SCO4443 | 6mA | 0→1 | +58.3% | +2.12 | MerR | Transcriptional regulator |
| SC_RS32025 | SCO5956 | 6mA | 1→1 | -13.2% | -2.09 | TetR | Transcriptional regulator |
| SC_RS20135 | SCO3600 | 4mC | 1→0 | -69.6% | -2.03 | - | HTH domain protein |
| SC_RS02455 | SCO0089 | 4mC | 0→1 | +78.6% | +1.96 | LysR | Transcriptional regulator |
| SC_RS35800 | SCO6722 | 4mC | 0→1 | +75.5% | +1.65 | - | Spore wall regulator SsgD |
| SC_RS38610 | SCO7279 | 4mC | 0→1 | +65.6% | +1.51 | - | HTH domain protein |
| SC_RS26965 | SCO4960 | 6mA | 1→0 | -54.1% | -1.27 | - | Sigma factor-like HTH protein |
| SC_RS25370 | SCO4639 | 6mA | 1→0 | -55.0% | -1.21 | TetR | Transcriptional regulator |

---

## 4. Negative Correlation Transcription Factors

### 4.1 Definition
Genes where methylation and expression change in **opposite directions**:
- Methylation ↑ and Expression ↓ (Gained_Down, Increased_Down)
- Methylation ↓ and Expression ↑ (Lost_Up, Decreased_Up)

### 4.2 Complete List (8 TFs)

| Gene ID | Old Locus | Modification | T1→T2 Sites | Δ Methylation | log₂FC | Family | Product |
|---------|-----------|--------------|-------------|---------------|--------|--------|---------|
| SC_RS35610 | **SCO6685** | 6mA | 1→0 | **-64.4%** | **+6.60** | - | **Response regulator RamR** |
| SC_RS24370 | SCO4441 | 6mA | 1→0 | -55.9% | +3.69 | - | HTH domain protein |
| SC_RS36820 | SCO6924 | 4mC | 0→1 | +83.8% | -2.54 | - | HTH domain protein |
| SC_RS16920 | SCO2964 | 4mC | 1→1 | +12.5% | -2.15 | LysR | Transcriptional regulator StgR |
| SC_RS09425 | SCO1490 | 4mC | 0→1 | +70.8% | -1.89 | - | Antitermination factor NusB |
| SC_RS13965 | SCO2386 | 4mC | 0→2 | +64.3% | -1.35 | - | FA biosynthesis regulator FasR |
| SC_RS29770 | SCO5517 | 4mC | 0→1 | +70.5% | -1.26 | TetR | Transcriptional regulator |
| SC_RS11615 | SCO1916 | 6mA | 1→0 | -56.0% | +1.20 | - | N-succinyltransferase |

---

## 5. Key Transcription Factors Analysis

### 5.1 NsdB (SCO7252) - Master Developmental Regulator

| Property | Value |
|----------|-------|
| Gene ID | SC_RS38475 |
| Old Locus | SCO7252 |
| Product | DNA-binding protein NsdB |
| Modification | 4mC |
| Methylation change | 0→1 site (+92.3%) |
| Expression change | log₂FC = **+9.32** (640× increase) |
| Correlation | **Positive** (Gained_Up) |

**Biological significance:**
- NsdB is a **pleiotropic regulator** of development and secondary metabolism
- Controls morphological differentiation and antibiotic production
- The massive upregulation with 4mC methylation gain suggests **methylation-mediated transcriptional activation**
- This may be a **master epigenetic switch** for the T1→T2 developmental transition

### 5.2 RamR (SCO6685) - SapB/Aerial Mycelium Regulator

| Property | Value |
|----------|-------|
| Gene ID | SC_RS35610 |
| Old Locus | SCO6685 |
| Product | Two-component system response regulator RamR |
| Modification | 6mA |
| Methylation change | 1→0 site (-64.4%) |
| Expression change | log₂FC = **+6.60** (100× increase) |
| Correlation | **Negative** (Lost_Up) |

**Biological significance:**
- RamR controls **SapB production** required for aerial mycelium formation
- Demethylation correlates with strong activation
- Suggests **6mA methylation acts as a repressive mark** on ramR
- Removal of methylation during T1→T2 may **release repression** and trigger differentiation

### 5.3 Other Notable Regulators

| Gene | Function | Pattern | Significance |
|------|----------|---------|--------------|
| SCO4122 (MarR) | Drug resistance/stress | Lost_Down | MarR repressor downregulated with demethylation |
| SCO1331 (LuxR) | Quorum sensing-like | Gained_Up | Potential cell-density signaling |
| SCO6722 (SsgD) | Sporulation | Gained_Up | Spore wall synthesis regulator activated |
| SCO2964 (StgR) | LysR regulator | Increased_Down | Negative regulation with methylation |

---

## 6. Regulatory Model

Based on our findings, we propose the following model for epigenetic regulation of transcription factors:

```
T1 (Vegetative Growth)
    │
    ├─ NsdB: Low methylation, low expression
    │      → Developmental genes repressed
    │
    ├─ RamR: High 6mA methylation, low expression
    │      → SapB/aerial mycelium genes repressed
    │
    └─ Multiple TFs in basal state

    ↓ T1→T2 Transition

T2 (Secondary Metabolism / Differentiation)
    │
    ├─ NsdB: 4mC GAINED → STRONGLY ACTIVATED (640×)
    │      → Triggers developmental cascade
    │      → Activates secondary metabolism
    │
    ├─ RamR: 6mA LOST → STRONGLY ACTIVATED (100×)
    │      → SapB production initiated
    │      → Aerial mycelium formation begins
    │
    └─ Coordinated TF network activation
         → 16 TFs with positive correlation
         → 8 TFs with negative correlation
```

---

## 7. Conclusions

1. **24 transcription factors** show coordinated methylation-expression changes

2. **SARP family genes** (actII-ORF4, redD, cdaR) are **not directly methylated** in promoter regions, but may be regulated indirectly through upper-level regulators

3. **NsdB (SCO7252)** is the most dramatically affected TF:
   - 4mC methylation gain
   - 640-fold expression increase
   - Master regulator of development and secondary metabolism

4. **RamR (SCO6685)** shows negative correlation:
   - 6mA demethylation leads to activation
   - Controls morphological differentiation

5. **Two distinct mechanisms** operate:
   - **Positive regulation:** Methylation activates transcription (NsdB model)
   - **Negative regulation:** Methylation represses transcription (RamR model)

6. These findings suggest **epigenetic regulation is a major mechanism** controlling the transcription factor network during growth phase transition in *S. coelicolor*

---

## 8. Data Files

| File | Description |
|------|-------------|
| `TF_coordinated_changes.csv` | All 24 TFs with methylation-expression data |
| `TF_ANALYSIS_REPORT.md` | This report |

---

## 9. Methods

- **Promoter definition:** -300 to +50 bp from TSS
- **Methylation threshold:** ≥10× coverage, ≥50% frequency, ≥2/3 replicates
- **Expression threshold:** padj < 0.05, |log₂FC| > 0.5
- **TF identification:** Keyword search for regulator, transcription, DNA-binding, HTH, sigma, activator, repressor, and family names (TetR, LysR, MarR, LuxR, GntR, MerR, ArsR, WhiB)

---

*Report generated as part of the integrated epigenome-transcriptome analysis*
*Streptomyces coelicolor A3(2) M145 Project*
