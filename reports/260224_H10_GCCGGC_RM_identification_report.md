# H10: Identification of the GCCGGC-Recognizing R-M System Responsible for 4mC Modification

**Date:** 2026-02-24
**Project:** S. coelicolor A3(2) M145 Epigenome-Transcriptome Integration
**Analysis:** 33_GCCGGC_RM_identification
**Status:** PARTIALLY SUPPORTED

---

## 1. Background

H9 (Loop 4) confirmed that CCGG 4mC is genuine N4-methylcytosine (not 5mC misclassification). Key findings from H9:
- 91.4% of CCGG 4mC sites are in the **GCCGGC palindromic context** (specifically the TGGCCGGC extended motif)
- **1,706 sites** detected at T1, declining to 486 (T2) and 30 (T3)
- The 5mC signal at these sites is essentially zero (0.01%)
- All three known DNA cytosine MTases in M145 (SC_RS19765/pseudogene, SC_RS19770, SC_RS36410) are **Dcm-like (produce 5mC)** and near-silent at T1

## 2. Hypothesis

**H10:** REBASE contains enzymes recognizing GCCGGC that produce N4-methylcytosine. A homolog exists in the M145 genome, potentially among the 129 MTases catalogued from GFF, or as an unannotated/mis-annotated gene.

## 3. Methods

### 3.1 REBASE Search
Parsed REBASE v602 bairoch format file (249,370 lines) for all entries with recognition sequences containing GCCGGC or GGCCGG. Extracted enzyme name, organism, modification type (m4C/m5C/m6A), and enzyme type classification.

### 3.2 Streptomyces Cross-reference
Cross-referenced with the 224 Streptomyces R-M system entries from `streptomyces_rm_systems.csv`.

### 3.3 M145 MTase Catalog
Used the H7-derived MTase catalog (129 entries from GFF annotations) merged with DESeq2 normalized expression counts across three timepoints (T1=18h, T2=36h, T3=60h).

### 3.4 Candidate Ranking
Scored candidates based on: DNA-specific annotation (+3), non-Dcm-like (+2), T1 expression level (+1-3), correlation with site count decline (+1-2), substrate specificity penalties.

## 4. Results

### 4.1 REBASE GCCGGC/GGCCGG Enzyme Census

| Metric | Count |
|--------|-------|
| Total REBASE entries matching GCCGGC/GGCCGG | 95 |
| GCCGGC exact (NaeI family) | 87 |
| CGCCGGCG (Sse232I family) | 4 |
| GGCCGGCC (FseI family) | 3 |
| YGCCGGCR (CcrNAV family) | 2 |

Among the 87 NaeI-family (GCCGGC) entries:

| Modification type | Count | Fraction |
|------------------|-------|----------|
| Unknown (REase only, no MTase characterized) | 49 | 56.3% |
| 5-methylcytosine (m5C) | 36 | 41.4% |
| **N4-methylcytosine (m4C)** | **2** | **2.3%** |

### 4.2 The Two Known GCCGGC m4C-Producing Enzymes

| Enzyme | Organism | REBASE Accession | Modification | Position |
|--------|----------|-----------------|-------------|----------|
| **M.Svi27968I** | *Streptomyces violascens* ATCC 27968 | RB455182 | Nm4C | Position 2 |
| **M.PfrJS2V** | *Propionibacterium freudenreichii* JS2 | RB177132 | Nm4C | Position 3 |

**Key observations:**
- Both organisms are **Actinobacteria**, phylogenetically close to *S. coelicolor*
- M.Svi27968I is from a *Streptomyces* species, making it the most relevant precedent
- The modification at position 2 of GCCGGC means the first C is methylated: G**C**CGGC -- consistent with our CCGG context where C at position 1 within GC**C**GGC is the modified base
- M.PfrJS2V modifies position 3 (GCC**G**GC -- but this is cytosine, so actually GC**C**GGC = position 3 from 5')
- The genome of *S. violascens* ATCC 27968 is available (GenBank CP029377, 6.92 Mb linear chromosome, 73.08% GC)
- Reference: Mao X., Gao K. (unpublished observations, REBASE)

### 4.3 Streptomyces GCCGGC R-M Landscape

Among 224 Streptomyces R-M entries in REBASE, **28 recognize GCCGGC** from **22 unique species**:

| Category | Count |
|----------|-------|
| REases (no MTase characterized) | 17 |
| MTases producing m5C | 11 |
| **MTases producing m4C** | **1 (M.Svi27968I)** |

GCCGGC is the **2nd most common recognition sequence** in *Streptomyces* R-M systems (after CCGCGG). Multiple species encode this system: *S. aureofaciens* (6 entries from different strains), *S. achromogenes*, *S. albus*, *S. albofaciens*, *S. ghanaensis*, *S. griseus*, *S. karnatakensis*, *S. longwoodensis*, *S. lusitanus*, *S. moderatus*, *S. xanthophaeus*, among others.

**The overwhelming majority produce m5C; only M.Svi27968I produces m4C.** This makes M145's GCCGGC m4C modification rare but precedented within the genus.

### 4.4 M145 MTase Catalog and Classification

The M145 genome encodes **129 MTases** (from GFF annotations). Classification:

| Category | Count |
|----------|-------|
| SAM-dependent MTase (unclassified) | 59 |
| RNA MTase | 18 |
| Protein MTase | 14 |
| Small molecule MTase | 13 |
| **DNA MTase (unspecified)** | **4** |
| MTase (general) | 4 |
| O-methyltransferase | 3 |
| **DNA cytosine MTase (Dcm-like, 5mC)** | **3** |
| DNA repair MTase | 2 |
| BREX system MTase (m6A) | 2 |
| N-6 adenine DNA MTase | 2 |
| Other | 5 |

### 4.5 DNA-Specific MTase Candidates

The 4 "DNA MTase (unspecified)" entries are the primary candidates:

| Locus tag | Product | Position | Strand | T1 expr | T2 expr | T3 expr | T2vsT1 log2FC |
|-----------|---------|----------|--------|---------|---------|---------|---------------|
| **SC_RS03950** | class I SAM-dependent DNA MTase | 428,657-429,403 | + | **130.2** | 160.9 | 102.3 | +0.29 |
| **SC_RS24685** | class I SAM-dependent DNA MTase | 4,924,984-4,925,628 | - | **324.8** | 33.8 | 66.1 | -3.22 |
| SC_RS19670 | DNA-methyltransferase | 3,878,757-3,879,503 | - | 3.8 | 2.9 | 48.0 | -0.21 |
| SC_RS36625 | DNA-methyltransferase | 7,654,146-7,654,901 | + | 2.7 | 3.9 | 46.0 | +0.31 |

**SC_RS19670** and **SC_RS36625** are annotated as generic "DNA-methyltransferase" without cytosine/adenine specificity, making them theoretical candidates. However, both are **near-silent at T1** (< 5 counts), making them unlikely to produce 1,516 modified sites. They dramatically increase at T3 (48 and 46 counts respectively).

**SC_RS24685** and **SC_RS03950** are annotated as "class I SAM-dependent DNA methyltransferase" and are well-expressed at T1.

### 4.6 Top 10 Ranked Candidates

| Rank | Score | Locus tag | Product | T1 | T2 | T3 | Key evidence |
|------|-------|-----------|---------|-----|-----|-----|-------------|
| 1 | 9 | **SC_RS24685** | class I SAM-dependent DNA MTase | 325 | 34 | 66 | DNA MTase, strong T1, dramatic decline |
| 2 | 9 | SC_RS28835 | BREX-2 PglX | 604 | 395 | 382 | DNA MTase, declining, but typically m6A |
| 3 | 8 | **SC_RS12530** | methyltransferase | 733 | 280 | 154 | Strong T1, monotonic decline |
| 4 | 8 | SC_RS12025 | class I SAM-dependent MTase | 231 | 96 | 71 | Strong T1, monotonic decline |
| 5 | 8 | **SC_RS03950** | class I SAM-dependent DNA MTase | 130 | 161 | 102 | DNA MTase annotation, constitutive |
| 6 | 8 | SC_RS07160 | methyltransferase | 758 | 220 | 135 | Strong T1, ~6-fold decline |
| 7 | 8 | SC_RS10665 | class I SAM-dependent MTase | 443 | 195 | 166 | Strong T1, declining |
| 8 | 8 | SC_RS25365 | class I SAM-dependent MTase | 434 | 267 | 146 | Strong T1, declining |
| 9 | 8 | SC_RS12890 | class I SAM-dependent MTase | 1462 | 379 | 157 | Very strong T1, ~9-fold decline |
| 10 | 8 | SC_RS13615 | class I SAM-dependent MTase | 410 | 162 | 93 | Strong T1, declining |

### 4.7 Assessment of Top Candidates

**SC_RS24685 (Rank 1, Score 9):**
- Annotated as "class I SAM-dependent DNA methyltransferase" (WP_003974444.1)
- T1 expression = 325 (adequate for 1,516 sites)
- **Dramatic T2 decline** (log2FC = -3.22, padj = 5.0e-36), partially recovers at T3
- This expression pattern (T1 >> T2 < T3) does not perfectly match the site count decline (T1 >> T2 >> T3), but the T1-to-T2 collapse parallels the 3-fold site reduction
- Located at position 4.92 Mb (core genome region)

**SC_RS03950 (Rank 5, Score 8):**
- Annotated as "class I SAM-dependent DNA methyltransferase" (WP_037666971.1)
- **Constitutive expression** across all timepoints (130, 161, 102)
- If this enzyme continuously methylates GCCGGC, site loss could be due to replication dilution without re-methylation of daughter strands
- Located at position 0.43 Mb

**SC_RS12530 (Rank 3, Score 8):**
- Annotated simply as "methyltransferase" (WP_030869367.1)
- Very strong T1 expression (733) with monotonic decline
- Substrate specificity unknown from annotation alone

### 4.8 Critical Negative Finding

**No M145 MTase has an explicit N4-cytosine annotation.** The only N4-C annotation in the GFF is for:
- SC_RS12500: "16S rRNA (cytosine(1402)-N(4))-methyltransferase RsmH" -- this targets rRNA, not DNA

This means the responsible enzyme is either:
1. A generic "SAM-dependent MTase" or "DNA methyltransferase" whose N4-C specificity is not reflected in the annotation
2. An unannotated or mis-annotated gene
3. A known protein with cryptic DNA N4-C MTase activity

## 5. Discussion

### 5.1 GCCGGC m4C is Rare but Precedented in Streptomyces

Out of 87 GCCGGC-recognizing R-M systems in REBASE, only 2 (2.3%) produce m4C. The overwhelming default is m5C. However, one of these two is **M.Svi27968I from *Streptomyces violascens***, demonstrating that GCCGGC m4C exists within the genus.

### 5.2 The M.Svi27968I Homology Search Imperative

The most direct path to identifying the M145 enzyme would be a **BLAST search of M.Svi27968I protein sequence against the M145 proteome**. The genome of *S. violascens* ATCC 27968 is available (GenBank CP029377), enabling extraction of the M.Svi27968I coding sequence. A homolog among M145's 59 unclassified "SAM-dependent MTases" or 4 generic "DNA MTases" would be the definitive identification.

### 5.3 Expression-Site Count Correlation

The GCCGGC site count shows a dramatic monotonic decline: 1,516 --> 486 --> 30.

Several MTases show expression patterns paralleling this decline:
- SC_RS12890: 1462 --> 379 --> 157 (9.3-fold decline)
- SC_RS07160: 758 --> 220 --> 135 (5.6-fold decline)
- SC_RS12530: 733 --> 280 --> 154 (4.8-fold decline)
- SC_RS24685: 325 --> 34 --> 66 (T1-to-T2 collapse, 4.9-fold)

However, expression correlation alone cannot distinguish DNA MTases from the many metabolic MTases that also decline during stationary phase. The annotation-based filter (DNA MTase + adequate T1 expression + not Dcm-like) narrows to **SC_RS24685** and **SC_RS03950** as the strongest candidates.

### 5.4 Alternative Hypothesis: Constitutive Enzyme + Replication Dilution

If the MTase is constitutively expressed (like SC_RS03950), the site count decline could be explained by:
- T1: Active methylation on freshly replicated DNA
- T2-T3: Growth slowdown, but ongoing cell division without complete re-methylation = passive dilution

This model does not require the MTase expression to decline.

## 6. Key Conclusions

| Finding | Detail |
|---------|--------|
| REBASE GCCGGC m4C precedent | 2 enzymes globally, including M.Svi27968I (*S. violascens*) |
| Streptomyces GCCGGC prevalence | 28/224 entries (12.5%), 2nd most common motif |
| M145 N4-C MTase annotation | **None** -- no gene has explicit N4-cytosine annotation |
| Top candidate by evidence | SC_RS24685 (DNA MTase, T1=325, dramatic decline) |
| Top candidate (constitutive) | SC_RS03950 (DNA MTase, T1=130, stable expression) |
| Required next step | **BLAST M.Svi27968I against M145 proteome** |

## 7. Verdict

**H10: PARTIALLY SUPPORTED**

- **SUPPORTED:** REBASE confirms that GCCGGC m4C exists, including in *Streptomyces* (M.Svi27968I)
- **NOT RESOLVED:** No M145 gene has explicit N4-cytosine MTase annotation. Candidates exist (SC_RS24685, SC_RS03950, SC_RS12530) but cannot be confirmed without sequence homology analysis
- **Next step:** BLAST search of M.Svi27968I protein against M145 proteome to identify the homolog among the 59 unclassified SAM-dependent MTases

## 8. Output Files

### Figures
| File | Description |
|------|-------------|
| `fig1_REBASE_GCCGGC_catalog.pdf/.svg` | REBASE GCCGGC enzyme catalog with modification types |
| `fig2_M145_MTase_heatmap.pdf/.svg` | M145 DNA/unclassified MTase expression heatmap |
| `fig3_candidate_expression_vs_sites.pdf/.svg` | Candidate expression vs GCCGGC site count timeline |
| `fig4_comprehensive_summary.pdf/.svg` | Comprehensive summary with all evidence panels |

### Tables
| File | Description |
|------|-------------|
| `REBASE_GCCGGC_GGCCGG_entries.tsv` | All 95 REBASE entries matching GCCGGC/GGCCGG |
| `REBASE_GCCGGC_NaeI_family.tsv` | 87 NaeI-family (exact GCCGGC) entries |
| `Streptomyces_GCCGGC_RM_systems.tsv` | 28 Streptomyces GCCGGC entries |
| `M145_MTase_complete_catalog.tsv` | Complete M145 MTase catalog with expression |
| `H10_ranked_candidates.tsv` | Ranked candidate list with evidence scores |

### Script
- `scripts/H10_GCCGGC_RM_identification.py`

---

*Analysis directory:* `/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/33_GCCGGC_RM_identification/`
*Report:* `/Users/okaban/bioinfo/rna-seq/reports/260224_H10_GCCGGC_RM_identification_report.md`
