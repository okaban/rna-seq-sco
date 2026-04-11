# H16: T1 GCCGGC MTase Expression-Correlated Reverse Identification

**Date**: 2026-02-26
**Analysis**: `11_epigenome_integration/analysis/39_GCCGGC_MTase_reverse_ID/`
**Hypothesis**: Among M145's MTases, one whose expression dynamics (high T1 -> declining T2/T3) correlates with GCCGGC site counts (1289->407->21), possesses DNA cytosine methyltransferase domain architecture, and has an REase partner nearby.

---

## 1. Background

GCCGGC 4mC sites in *S. coelicolor* M145 show dramatic temporal dynamics: 1,289 at T1, 407 at T2, and only 21 at T3. Previous attempts to identify the responsible N4-cytosine methyltransferase failed:

- **H10 (annotation-based)**: Top candidate SC_RS24685 (SCO4573, DNA MTase, T1=325, score 9/10) had no BLAST similarity to known GCCGGC MTases
- **H12 (BLAST-based)**: SC_RS19770 (E=0.007) and SC_RS36410 (E=0.047) showed sequence similarity to M.Svi27968I (a known GCCGGC m4C MTase from *S. violascens*), but both have very low T1 expression (11 and 1.3, respectively), far insufficient to sustain 1,289 sites

The present analysis takes a **reverse approach**: instead of starting from annotation or sequence homology, we rank all 129 MTases in M145 by how well their expression dynamics correlate with the observed GCCGGC site counts.

## 2. Methods

### 2.1 Spearman Correlation
For each of the 129 MTases in the complete catalog, Spearman rank correlation was computed between the expression vector [T1_mean, T2_mean, T3_mean] and the GCCGGC site count vector [1289, 407, 21].

### 2.2 Filtering
Candidates required:
- Spearman rho > 0.8 (strong positive correlation with site dynamics)
- T1 expression > 50 (sufficient to sustain methylation activity)

### 2.3 Domain Classification
Each candidate was classified using NCBI CDD (Conserved Domain Database) and product annotations:
- **DNA-C**: DNA cytosine methyltransferase
- **DNA-UNSPECIFIED**: DNA MTase, specificity unknown
- **UNCLASSIFIED**: SAM-dependent MTase, substrate unknown
- **NON-DNA**: RNA, protein, small molecule, or other non-DNA MTase

### 2.4 Genomic Neighborhood
For each candidate, genes within +/- 10 kb were examined for R-M system components (restriction endonucleases, specificity subunits, defense-related genes).

### 2.5 Composite Scoring
A 13-point composite score integrated:
- Spearman rho (0-3 points)
- T1 expression adequacy / sites-per-expression ratio (0-3 points)
- Domain classification (0-3 points)
- R-M neighbor presence (0-2 points)
- Cross-validation with H10/H12 (0-2 points)

## 3. Results

### 3.1 Spearman Correlation Distribution

Of 129 MTases, **37 showed rho > 0.8** (all with rho = 1.0 due to the 3-timepoint monotonic decrease). Of these, **36 had T1 expression > 50**. With only 3 timepoints, a monotonically decreasing expression pattern automatically yields rho = 1.0, so the correlation filter alone is insufficient -- domain classification and neighborhood analysis are essential for discrimination.

### 3.2 Top Candidates by Composite Score

| Rank | Locus | SCO | Score | rho | T1 | T2 | T3 | Category | Key Evidence |
|------|-------|-----|-------|-----|-----|-----|-----|----------|-------------|
| 1 | **SC_RS12890** | SCO2170 | **10** | 1.00 | 1462 | 379 | 157 | SAM MTase (unclassified) | Sites/expr=0.88; R-M neighbor (SC_RS12910 MTase); H10 hit |
| 2 | **SC_RS25365** | SCO4638 | **10** | 1.00 | 434 | 267 | 146 | SAM MTase (unclassified) | McrA (type IV REase) 10 kb away; H10 hit |
| 3 | **SC_RS12025** | SCO1897 | **10** | 1.00 | 231 | 96 | 71 | SAM MTase (unclassified) | R-M neighbor (SC_RS11995 MTase + helicase); H10 hit |
| 4 | **SC_RS12530** | SCO2243 | **10** | 1.00 | 733 | 280 | 154 | MTase (general) | Neighbor RsmH MTase; H10 hit |
| 5 | **SC_RS13615** | SCO2317 | **10** | 1.00 | 410 | 162 | 93 | SAM MTase (unclassified) | **CDD: COG0863 YhdJ DNA modification methylase hit (E=2.23e-03)**; neighbor TnpB (IS200/IS605); H10 hit |
| 6 | SC_RS03950 | SCO0381 | 9 | 0.50 | 130 | 161 | 102 | DNA MTase (unspecified) | Annotated DNA MTase; R-M neighbor; H10 hit |
| 7 | SC_RS38280 | SCO7400 | 9 | 1.00 | 101 | 83 | 51 | SAM MTase (unclassified) | Adjacent Type 11 MTase; H10 hit |

### 3.3 Domain Analysis (NCBI CDD)

CDD searches for the top 5 candidates revealed:

| Candidate | Top CDD Hits | DNA MTase Domain? |
|-----------|-------------|-------------------|
| SC_RS12890 (WP_011028174.1) | COG2227 (UbiG), PRK11036 (CmoM), COG2226 (UbiE) | **No** - small molecule/tRNA MTase domains |
| SC_RS25365 (WP_011029786.1) | COG2226 (UbiE), COG2227 (UbiG), COG4106 | **No** - small molecule MTase domains |
| SC_RS12025 (WP_003976818.1) | COG2226 (UbiE), Methyltransf_11/25 | **No** - ubiquinone/menaquinone biosynthesis |
| SC_RS12530 (WP_030869367.1) | PRK11036 (CmoM), COG2227, COG2226 | **No** - tRNA/lipid MTase domains |
| **SC_RS13615** (WP_003976496.1) | COG2226 (UbiE, E=3.15e-35), **COG0863 (YhdJ DNA modification methylase, E=2.23e-03)** | **Weak hit** - minor COG0863 DNA methylase signal |

**Critical finding**: SC_RS13615 is the only top candidate with any CDD evidence for DNA methyltransferase activity (COG0863/YhdJ, though at a marginal E-value).

### 3.4 SCO1731 (SC_RS10665): Known DNA Cytosine MTase

A literature-critical finding: **SC_RS10665 = SCO1731**, previously characterized as a DNA cytosine methyltransferase that modulates actinorhodin production (Pisciotta et al. 2018, Sci Rep). Key data:
- rho = 1.0, T1 = 443, T2 = 195, T3 = 166
- Sites/expression ratio at T1: 1289/443 = 2.91
- **However**: SCO1731 was shown to produce **m5C** (5-methylcytosine), not m4C (N4-methylcytosine)
- The 2023 Pisciotta et al. study (Sci Rep) identified GGC**m5**CGG and GCC**m5**CG motifs by bisulfite sequencing
- Our GCCGGC sites are **4mC** (N4-methylcytosine), a fundamentally different modification

**Conclusion**: SCO1731/SC_RS10665 is the m5C methyltransferase for GGCCGG/GCCG but is **not** the enzyme responsible for the m4C at GCCGGC sites. The two modifications (m5C by SCO1731 and m4C by an unknown enzyme) co-exist at overlapping sequence contexts.

### 3.5 Genomic Neighborhood Analysis

Key findings for top candidates:

| Candidate | Notable Neighbors within 10 kb |
|-----------|-------------------------------|
| **SC_RS12890** | SC_RS12910 (class I SAM MTase, 3 kb away, convergent orientation) |
| **SC_RS25365** | **SC_RS25315/mcrA** (type IV restriction endonuclease McrA, 8 kb away) |
| **SC_RS13615** | SC_RS13590 (class I SAM MTase, 4 kb away); SC_RS13585 (TnpB, IS200/IS605 transposase accessory); SC_RS13610 (prenyltransferase) |
| SC_RS12025 | SC_RS11995 (class I SAM MTase); SC_RS12030 (ATP-dependent helicase) |
| SC_RS12530 | SC_RS12500 (RsmH, 16S rRNA MTase -- not R-M related) |

The **SC_RS25365 -- mcrA** pairing is particularly notable: McrA is a methyl-dependent restriction endonuclease that cleaves methylated DNA. This represents a type IV R-M system architecture.

### 3.6 Expression Dynamics Comparison

Reference: SC_RS17645 (AAGCCCG MTase): T1=245, 260 sites, ratio = 1.06 sites/count

| Candidate | T1 | T2 | T3 | T1 ratio | T2 ratio | T3 ratio | Decline pattern |
|-----------|-----|-----|-----|----------|----------|----------|----------------|
| GCCGGC sites | 1289 | 407 | 21 | - | - | - | 100% -> 31.6% -> 1.6% |
| SC_RS12890 | 1462 | 379 | 157 | 0.88 | 1.07 | 0.13 | 100% -> 25.9% -> 10.7% |
| SC_RS10665 (SCO1731) | 443 | 195 | 166 | 2.91 | 2.09 | 0.13 | 100% -> 44.0% -> 37.5% |
| SC_RS25365 | 434 | 267 | 146 | 2.97 | 1.52 | 0.14 | 100% -> 61.5% -> 33.6% |
| SC_RS13615 | 410 | 162 | 93 | 3.14 | 2.51 | 0.23 | 100% -> 39.5% -> 22.6% |

**SC_RS12890** has the best sites-to-expression ratio (0.88 at T1, close to the AAGCCCG reference of 1.06) and the steepest decline that most closely mirrors the GCCGGC site count decrease. However, its CDD profile shows no DNA MTase domains.

### 3.7 Cross-Validation Summary

| Candidate | In H10? | In H12 BLAST? | CDD DNA MTase? | R-M Neighbor? | Literature? |
|-----------|---------|--------------|----------------|---------------|-------------|
| SC_RS12890 | Yes (score 8) | No | No | SC_RS12910 (MTase) | Unknown |
| SC_RS25365 | Yes (score 8) | No | No | mcrA (REase) | Unknown |
| SC_RS13615 | Yes (score 8) | No | Weak (COG0863) | TnpB | Unknown |
| SC_RS12025 | Yes (score 8) | No | No | MTase + helicase | Unknown |
| SC_RS10665 (SCO1731) | Yes (score 8) | No | m5C producer | No R-M | Known m5C MTase |
| SC_RS19770 | Yes (score 3) | **Yes** (E=0.007) | DNA cytosine MTase | Defense island | T1=11 (too low) |
| SC_RS36410 | Yes (score 3) | **Yes** (E=0.047) | DNA cytosine MTase | Defense island | T1=1.3 (too low) |

## 4. Synthesis: The Two-Enzyme Problem

The analysis reveals a fundamental tension between two lines of evidence:

### Expression-correlated candidates (this analysis)
- SC_RS12890: Best expression fit (T1=1462, ratio=0.88) but **no DNA MTase domain**
- SC_RS13615: Strong expression (T1=410) with **weak CDD DNA methylase hit** (COG0863)
- SC_RS10665/SCO1731: Known cytosine MTase but produces **m5C, not m4C**

### Sequence-homology candidates (H12)
- SC_RS19770: BLAST hit to M.Svi27968I (E=0.007), annotated DNA cytosine MTase, but **T1=11** (rho=-0.5)
- SC_RS36410: BLAST hit (E=0.047), DNA cytosine MTase in defense island, but **T1=1.3** (rho=-1.0)

**No single candidate satisfies all criteria**: high T1 expression, declining dynamics, DNA cytosine MTase domain, and BLAST homology to known GCCGGC MTases.

### Possible Explanations

1. **Post-transcriptional regulation**: The responsible MTase may be translated more efficiently than its mRNA levels suggest, or stabilized at the protein level. SC_RS19770 or SC_RS36410 could potentially maintain high protein levels despite low transcript counts.

2. **Cooperative/redundant system**: Multiple MTases may contribute to GCCGGC methylation, each individually insufficient but collectively responsible.

3. **Cryptic DNA MTase**: One of the "unclassified SAM-dependent MTases" (SC_RS12890, SC_RS13615) may have diverged sufficiently from known DNA MTases to be unrecognizable by current CDD profiles, yet retain DNA cytosine methylation activity. The weak COG0863 hit in SC_RS13615 supports this possibility.

4. **SCO1731 dual activity**: SCO1731 (SC_RS10665) may produce both m5C and m4C at GCCGGC sites. This has not been tested experimentally. The Pisciotta et al. studies used bisulfite sequencing, which detects m5C but not m4C, leaving the question open.

## 5. Revised Candidate Ranking

Based on all evidence, prioritized candidates for experimental validation:

| Priority | Locus | Rationale |
|----------|-------|-----------|
| **1** | **SC_RS13615** (SCO2317) | Composite score 10; rho=1.0; T1=410; **only candidate with CDD DNA methylase signal** (COG0863); TnpB defense neighbor; expression dynamics match |
| **2** | **SC_RS10665** (SCO1731) | Known m5C MTase at GGC**m5**CGG; potential dual m4C activity; rho=1.0; T1=443; literature-validated DNA cytosine MTase |
| **3** | **SC_RS12890** (SCO2170) | Best expression-site ratio (0.88); rho=1.0; T1=1462; paired with SC_RS12910 MTase; no DNA domain signal but highest expression match |
| **4** | **SC_RS19770** | BLAST homology to M.Svi27968I (E=0.007); annotated DNA cytosine MTase; low T1 but possible post-transcriptional regulation |
| **5** | **SC_RS25365** (SCO4638) | Composite score 10; mcrA REase 8 kb away (Type IV R-M); but no DNA domain |

## 6. Key Figures

- `figures/correlation_scatter.png`: T1 expression vs Spearman rho scatter, colored by domain class
- `figures/expression_timelines.png`: Top candidate expression dynamics alongside GCCGGC site counts, normalized comparison, and sites-per-expression ratios
- `figures/genomic_neighborhoods.png`: Genomic context diagrams for top candidates

## 7. Output Files

| File | Description |
|------|-------------|
| `tables/ranked_MTase_correlations.tsv` | All 129 MTases ranked by Spearman rho |
| `tables/top_candidates_detail.tsv` | 57 candidates (rho >= 0.5, T1 > 50) with composite scores |
| `tables/genomic_neighbors.tsv` | 48 R-M related neighbor genes |
| `scripts/H16_reverse_ID_analysis.py` | Complete analysis script |

## 8. Verdict

### **PARTIAL** -- Expression-correlated candidates identified but none conclusively confirmed as the GCCGGC N4-C MTase

**Supporting evidence**:
- 36 MTases show perfect expression-site correlation (rho=1.0) with T1 > 50
- SC_RS13615 (SCO2317) is the strongest novel candidate: composite score 10, CDD DNA methylase domain hit (COG0863, E=2.23e-03), TnpB defense neighbor, declining expression (410 -> 162 -> 93)
- SC_RS10665 (SCO1731) is a literature-validated DNA cytosine MTase with matching expression dynamics, but produces m5C not m4C
- The sites-to-expression ratio analysis shows SC_RS12890 most closely matches the expected MTase activity level

**Unsupporting evidence**:
- No candidate simultaneously satisfies all four criteria: (i) high T1 expression, (ii) declining dynamics, (iii) confirmed DNA cytosine MTase domain, and (iv) BLAST homology to GCCGGC MTases
- CDD searches for the top 5 expression-correlated candidates show small molecule/tRNA MTase domains, not DNA MTase domains (except the weak COG0863 hit in SC_RS13615)
- The BLAST-validated candidates (SC_RS19770, SC_RS36410) have T1 expression far too low to explain 1,289 sites
- With only 3 timepoints, the Spearman correlation analysis cannot distinguish the true GCCGGC MTase from any other gene with monotonically decreasing expression (37 of 129 MTases show rho=1.0)

**Next steps for experimental validation**:
1. Gene knockout of SC_RS13615 and SC_RS10665 followed by SMRT-seq to test GCCGGC 4mC loss
2. Protein-level analysis (Western blot or proteomics) for SC_RS19770 to assess post-transcriptional regulation
3. In vitro methylation assay with purified SC_RS13615 and SC_RS10665 protein on GCCGGC-containing substrates

---

*Analysis: `/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/39_GCCGGC_MTase_reverse_ID/`*
*Script: `scripts/H16_reverse_ID_analysis.py`*
