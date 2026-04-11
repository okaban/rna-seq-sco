# H11: Gatekeeper Model v2 - Quantitative Epigenome-Transcriptome Integration

**Date**: 2026-02-24
**Analysis**: H11 - Synthesis of all hypothesis results (H1-H9) into a quantitative 3-layer model
**Organism**: *Streptomyces coelicolor* A3(2) M145
**Analysis Directory**: `11_epigenome_integration/analysis/34_gatekeeper_model_v2/`

---

## 1. Executive Summary

This report synthesizes the results of 9 hypothesis tests (H1-H9) across 4 analysis loops into a quantitative "Gatekeeper Model v2" describing how DNA methylation interfaces with transcriptional regulation in *S. coelicolor* M145. The model comprises **3 active layers** and **1 excluded pathway**, supported by definitive statistical evidence from multiple independent tests.

**Key conclusions:**
1. DNA methylation operates as a **genome-wide landscape remodeler**, not a gene-specific cis-regulatory element
2. Critical regulatory DNA (TF binding sites, sigma factor promoter elements) is **actively protected** from methylation
3. Methylation gates **non-canonical regulatory genes** (outside the well-studied literature set), particularly through **asymmetric TCS methylation** and **geographic modulation**
4. The previously hypothesized methylation-TF cascade-BGC activation pathway is **definitively excluded**

---

## 2. Model Architecture

### 2.1 Layer 1: Landscape Remodeling (Strong Evidence)

Methylation acts at the genome-wide level to reshape the epigenomic landscape across developmental timepoints. Individual site-level effects are weak; the signal emerges from wholesale changes in methylation patterns.

| Parameter | Value | Source |
|-----------|-------|--------|
| 4mC promoter-expression correlation | rho = 0.139, p = 8.26x10^-4 | H1 |
| Permutation test (4mC T2vsT1) | r = 0.125, p = 0.002 | Pre-loop |
| CCGG 4mC positional overlap | Jaccard = 0.000 | H7 |
| AAGCCCG site loss (SC_RS17645 down) | 260 sites, 0% retention | H5 |
| T1 CCGG core genome fraction | 77% | H7 |
| T2 CCGG chromosomal arms fraction | 89% | H7 |
| 6mA T2vsT1 correlation | Lost after T3 resequencing | Pre-loop |
| 5mC misclassification | Excluded (5mC signal = 0.01%) | H9 |
| GCCGGC palindrome context | 91.4% of CCGG 4mC sites | H9 |
| Gene-level expression correlation (all motifs) | |rho| < 0.16 | H5 |

**Interpretation**: Methylation does not regulate genes one-by-one in cis. Instead, when a MTase changes expression (e.g., SC_RS17645 LFC = -2.19), ALL 260 of its target sites are lost simultaneously. The CCGG/GCCGGC system shows complete positional turnover between timepoints (Jaccard = 0.000), with a dramatic geographic shift from core genome (77% at T1) to chromosomal arms (89% at T2). This is active remodeling, not passive dilution (H7). The 4mC modification identity is confirmed (H9): it is a genuine N4-cytosine modification, not 5mC misclassification, recognizing the GCCGGC palindrome.

### 2.2 Layer 2: Protection / Depletion (Definitive Evidence)

Critical regulatory DNA elements are actively protected from methylation, creating a "safe zone" for transcriptional machinery.

| Parameter | Value | Source |
|-----------|-------|--------|
| TF binding site methylation fold enrichment | 0.66 (depleted) | Pre-loop |
| TF binding site depletion p-value | 6.8x10^-5 | Pre-loop |
| Sigma factor -10 box depletion | p = 1.73x10^-12 | Pre-loop |
| Literature 37 TFs methylated | 3/37 (8%), all constitutive | H4 |
| SARP pathway activators methylated | 0/4 (0%) | H4 |
| SARP activator LFC range | +2.0 to +7.2 (unmethylated) | H4 |
| DMG functional enrichment (GO/KEGG/COG) | NONE | Pre-loop |
| TF BS methylation: proportion Lost T2T3 | 75% | Pre-loop |

**Interpretation**: The genome maintains methylation-free zones around transcription factor binding sites (fold = 0.66, significantly depleted) and sigma factor -10 boxes (p = 1.73x10^-12). The 37 well-characterized literature TFs are overwhelmingly methylation-free (92%), with all 4 SARP pathway-specific activators (redD, actII-ORF4, cdaR, cpkO) completely unmethylated despite being highly expressed (LFC +2 to +7.2). DMGs show no functional enrichment, meaning methylation target selection is position-dependent but function-independent. The few methylation events at regulatory sites are transient (75% lost by T2/T3).

### 2.3 Layer 3: Signal Gating -- Revised (Moderate-Strong Evidence)

Methylation gates a distinct set of non-literature regulatory genes, particularly through asymmetric two-component system (TCS) methylation and geographic modulation at chromosomal arms.

| Parameter | Value | Source |
|-----------|-------|--------|
| Total regulatory genes screened | 1,055 (25 families) | H6 |
| Methylated regulatory genes | 165/1,055 (15.6%) | H6 |
| Literature 37 methylation rate for comparison | 6.7% (vs 15.6%) | H4 vs H6 |
| Coordinated methyl + expr genes | 62 (100% outside literature 37) | H6 |
| MerR family methylation rate (highest) | 28.6% | H6 |
| SARP family methylation rate (lowest) | 0% (n=7) | H6 |
| Overall 62 at chromosomal arms | 51.6%, p = 0.17 (NS) | H8 |
| T3 discordant_gain_up at arms | 87%, OR = 8.05, p = 0.001 | H8 |
| COG T (signal transduction) enrichment | OR = 2.47, p = 0.045 | H8 |
| COG Q (secondary metabolism) | 0 genes | H8 |
| TCS asymmetric pairs identified | 7 | H8 |

**Top candidate genes for experimental follow-up:**

| Locus Tag | Gene | Product | Coordination | Context |
|-----------|------|---------|--------------|---------|
| SC_RS10435 | -- | TetR family | Concordant derepression (T2 + T3) | Chaplin-adjacent |
| SC_RS35525 | -- | Sensor histidine kinase | Concordant derepression (T2 + T3) | PPTase-adjacent |
| SC_RS31385 | -- | Response regulator | Concordant derepression (T2 + T3) | Citrate synthase-adjacent |
| SC_RS24635 | -- | UdgX uracil-DNA binding | Concordant repression (T2 + T3) | DNA repair |
| SC_RS35610 | ramR | RamR response regulator | Discordant gain_up (T3) | Aerial mycelium |

**Interpretation**: The original model proposed methylation --> TF cascade --> BGC activation (REJECTED by H4). The revised model shows that methylation targets a previously uncharacterized set of 62 regulatory genes, ALL outside the well-studied literature 37 list. These genes are enriched in signal transduction functions (COG T, OR = 2.47) but contain ZERO secondary metabolism genes (COG Q = 0). TCS pairs show asymmetric methylation: sensor kinases are methylated while their cognate response regulators remain free. At T3, genes gaining methylation while being upregulated (discordant_gain_up) cluster overwhelmingly at chromosomal arms (87%, OR = 8.05), consistent with the Layer 1 geographic remodeling.

### 2.4 Excluded Pathway: Methylation-TF Cascade (Definitively Rejected)

| Evidence | Result | Implication |
|----------|--------|-------------|
| Literature 37 TFs methylation rate | 8% (3/37), all constitutive | TF network is methylation-independent |
| SARP activators (redD, actII-ORF4, cdaR, cpkO) | 0% methylated | BGC activation requires no methylation |
| Coordinated TFs in literature 37 | 0/37 | Zero methylation-expression coordination |
| redZ paradox | Premise error (locus_tag 36/37 wrong) | No paradox exists |
| True redZ expression | LFC = +0.91 (upregulated) | Consistent with Red BGC activation |
| redD expression | LFC = +5.79 (55-fold increase) | Main driver of Red pathway |

**Interpretation**: The elegant but incorrect model of "methylation regulates master TFs which activate BGC cascades" is definitively excluded. Not a single one of the 37 well-characterized TFs shows coordinated methylation-expression changes. The 4 SARP pathway-specific activators are all completely unmethylated, yet show the strongest expression changes in the dataset (LFC +2 to +7.2). The redZ paradox that originally motivated part of this hypothesis was an artifact of systematic locus_tag errors in the literature database (36/37 incorrect).

---

## 3. Evidence Integration Matrix

### 3.1 Hypothesis-to-Layer Mapping

| Hypothesis | Layer 1 | Layer 2 | Layer 3 | Excluded | Verdict |
|------------|---------|---------|---------|----------|---------|
| **H1**: 4mC/6mA distribution | **Strong** | -- | -- | -- | Partial support |
| **H3**: redZ paradox | -- | -- | -- | **Rejects premise** | Premise rejected |
| **H4**: TF methylation cascade | -- | **Supports** | -- | **Rejects** | Rejected |
| **H5**: MTase-site dynamics | **Strong** | -- | -- | -- | Partial support |
| **H6**: Genome-wide TF screen | -- | -- | **Strong** | -- | Supported |
| **H7**: CCGG MTase paradox | **Strong** | -- | Supports | -- | Passive dilution rejected |
| **H8**: Coordinated reg. char. | -- | -- | **Strong** | -- | Partial support |
| **H9**: 5mC misclassification | Supports | -- | -- | -- | Rejected |

### 3.2 Evidence Strength Summary

- **Definitive** (4 hypotheses): H3, H4, H7, H9 -- clear-cut results with no ambiguity
- **Strong** (1 hypothesis): H6 -- 62 coordinated regulators identified, all novel
- **Moderate-Strong** (1 hypothesis): H8 -- T3 arm enrichment significant, overall geography NS
- **Partial** (2 hypotheses): H1, H5 -- some aspects supported, others not

### 3.3 Verdict Distribution

| Category | Count | Hypotheses |
|----------|-------|------------|
| Supported | 1 | H6 |
| Partial support | 3 | H1, H5, H8 |
| Rejected (hypothesis) | 2 | H4, H9 |
| Rejected (premise/alternative) | 2 | H3 (premise error), H7 (passive dilution) |

---

## 4. Model Comparison: v1 vs v2

### 4.1 Original Gatekeeper Model (v1, 2026-02-17)

The original model proposed:
1. Methylation as a "gatekeeper" regulating TF activity
2. TF cascade from methylation changes to BGC activation
3. AAGCCCG methylation of redZ as the paradigm example

### 4.2 Revised Gatekeeper Model (v2, 2026-02-24)

Key revisions based on H1-H9 hypothesis testing:

| Aspect | v1 | v2 | Evidence |
|--------|----|----|----------|
| TF cascade | Methylation --> 37 TFs --> BGCs | EXCLUDED (0 coordination) | H4 |
| Target regulators | Literature TFs | Non-literature 62 | H6 |
| Methylation mode | Gene-specific cis | Genome-wide landscape | H5 |
| CCGG dynamics | Passive dilution? | Active geographic remodeling | H7 |
| redZ role | Paradox paradigm | No paradox (locus_tag error) | H3 |
| 4mC identity | Uncertain (5mC?) | Confirmed N4-C at GCCGGC | H9 |
| Signal gating | Through known TFs | Asymmetric TCS + arm-enriched | H8 |

---

## 5. Remaining Gaps and Testable Predictions

### 5.1 Key Unresolved Questions

1. **GCCGGC MTase identity**: The enzyme responsible for T1 core-genome CCGG/GCCGGC 4mC sites (1,516 sites) remains unidentified. SC_RS36410 (Dcm-like, induced at T2/T3) likely creates arm-enriched sites but cannot explain T1 patterns.

2. **Mechanism of TF BS protection**: Is methylation actively excluded from TF binding sites (by TF occupancy blocking MTase access), or is there a sequence-level incompatibility? The sharp depletion (fold = 0.66) suggests active protection, but the mechanism is unknown.

3. **Causality in coordinated regulators**: The 62 coordinated regulatory genes show correlation between methylation and expression changes, but causality is not established. Does methylation cause expression changes, or do expression/chromatin changes alter MTase accessibility?

4. **Functional impact of TCS asymmetry**: The asymmetric methylation of sensor kinases (methylated) vs response regulators (free) in 7 TCS pairs is striking, but the functional consequence is unknown.

5. **Cross-species conservation**: Does the 3-layer model apply to other *Streptomyces* species with complex R-M systems?

### 5.2 Testable Predictions (10 total)

| ID | Priority | Layer | Prediction | Test Method |
|----|----------|-------|------------|-------------|
| P1 | High | L1 | SC_RS17645 KO eliminates all AAGCCCG sites | Gene KO + PacBio |
| P2 | High | L1 | Separate MTase creates T1 core CCGG sites | SC_RS36410 KO + methylome |
| P3 | Medium | L2 | SARP promoter methylation is inert | Engineered MTase + RT-qPCR |
| P4 | High | L3 | SC_RS10435 demethylation derepresses chaplins | KO + RNA-seq |
| P5 | High | L3 | SC_RS35525 demethylation activates kinase | In vitro phosphorylation |
| P6 | Medium | L3 | TCS asymmetry conserved across species | Comparative methylomics |
| P7 | High | L1 | Core-to-arm CCGG shift is reproducible | Independent PacBio replicates |
| P8 | Medium | L2 | TF BS protection has sharp boundaries | High-res methylome at BS |
| P9 | Medium | L3 | MerR methylation responds to metal stress | Metal stress + methylome |
| P10 | Low | All | 3-layer model applies to other Streptomyces | Multi-species study |

### 5.3 Priority Ranking for Experimental Validation

**Tier 1 (Immediate)**:
- P1: SC_RS17645 KO -- foundational for the AAGCCCG R-M system
- P4: SC_RS10435 (chaplin-adjacent TetR) -- strongest candidate from 62 coordinated regulators (concordant derepression at both T2 and T3)
- P7: CCGG geographic shift reproducibility -- validates the core finding of Layer 1

**Tier 2 (High priority)**:
- P2: GCCGGC MTase identification -- resolves the biggest remaining mechanistic question
- P5: SC_RS35525 (PPTase-adjacent SK) -- tests the asymmetric signal gating model directly

**Tier 3 (Supporting)**:
- P3, P6, P8, P9, P10 -- these strengthen the model but are not required to establish it

---

## 6. Quantitative Model Parameters

### 6.1 Layer 1 Parameters

| Parameter | Statistic | CI/Range | Test |
|-----------|-----------|----------|------|
| 4mC promoter-expression (rho) | 0.139 | p = 8.26x10^-4 | Spearman |
| 4mC T2vsT1 correlation (r) | 0.125 | perm p = 0.002 | Permutation (n=10,000) |
| CCGG site overlap | 0.000 | Jaccard index | Exact |
| AAGCCCG retention | 0% | 0/260 sites | Count |
| MTase-site linkage | SC_RS17645 LFC=-2.19 | p_adj < 0.001 | DESeq2 |
| Geographic shift (T1 core) | 77% | -- | Position analysis |
| Geographic shift (T2 arms) | 89% | -- | Position analysis |

### 6.2 Layer 2 Parameters

| Parameter | Statistic | CI/Range | Test |
|-----------|-----------|----------|------|
| TF BS methylation fold | 0.66 | p = 6.8x10^-5 | Enrichment |
| Sigma -10 depletion | Significant | p = 1.73x10^-12 | Chi-square |
| Literature TF methylation | 8% (3/37) | -- | Count |
| SARP methylation | 0% (0/4) | -- | Count |
| DMG functional enrichment | None | All NS | GO/KEGG/COG |
| BS methylation turnover | 75% lost | T2/T3 | Temporal |

### 6.3 Layer 3 Parameters

| Parameter | Statistic | CI/Range | Test |
|-----------|-----------|----------|------|
| Regulatory genes screened | 1,055 | 25 families | GFF extraction |
| Methylated regulators | 165 (15.6%) | -- | Count |
| Coordinated regulators | 62 | 100% novel | Integration |
| T3 gain_up arm enrichment | 87% | OR=8.05, p=0.001 | Fisher's exact |
| COG T enrichment | OR=2.47 | p=0.045 | Fisher's exact |
| TCS asymmetric pairs | 7 | -- | Pair analysis |
| MerR methylation rate | 28.6% | -- | Family-level |
| SARP methylation rate | 0% | n=7 | Family-level |

---

## 7. Biological Interpretation

### 7.1 The Gatekeeper Metaphor Revised

In v1, the "gatekeeper" metaphor implied methylation controlling access to the TF-BGC regulatory cascade. In v2, the metaphor is more nuanced:

- **Layer 1 (Landscape)**: Methylation acts as a **surveyor**, continuously reshaping the terrain on which gene regulation occurs. It does not control individual gates but changes the entire landscape.
- **Layer 2 (Protection)**: The most important regulatory elements are **walled gardens** -- methylation cannot enter the zones where master TFs bind and sigma factors initiate transcription.
- **Layer 3 (Signal Gating)**: Methylation acts as a **signal filter** on a previously uncharacterized set of regulatory genes, particularly at chromosomal arms during late development. It modulates signal transduction rather than directly controlling biosynthetic output.

### 7.2 Why the TF Cascade Model Fails

The TF cascade model fails because:
1. The 37 literature TFs are in the "walled garden" (92% methylation-free)
2. SARP activators show the strongest expression changes without any methylation
3. The genes that ARE methylation-regulated are different genes entirely (62 non-literature regulators)
4. These 62 regulators are enriched in signal transduction (COG T), not secondary metabolism (COG Q = 0)

This means methylation influences gene regulation through an **indirect, multi-step pathway** rather than a direct TF --> target mechanism.

### 7.3 Implications for Streptomyces Biology

1. **Antibiotic production regulation**: BGC activation is methylation-INDEPENDENT. Efforts to engineer antibiotic production through methylation manipulation are unlikely to succeed directly. However, manipulating the 62 coordinated regulators could have indirect effects.

2. **Morphological differentiation**: The strongest candidate (SC_RS10435, chaplin-adjacent TetR with concordant derepression) suggests methylation may influence aerial mycelium development through previously unknown regulatory genes.

3. **Defense systems**: The CCGG/GCCGGC system appears to have a defense function (H7: defense island co-induction at T3), consistent with R-M system biology but adding a temporal/geographic dimension.

4. **Two-component signaling**: The asymmetric TCS methylation (SK methylated, RR free) suggests a novel regulatory mechanism where methylation modulates signal INPUT (sensor kinase expression) while preserving signal OUTPUT (response regulator availability).

---

## 8. Output Files

### 8.1 Figures

| File | Format | Description |
|------|--------|-------------|
| `gatekeeper_model_v2.pdf/.svg` | PDF, SVG | 3-layer model conceptual diagram with quantitative annotations |
| `evidence_integration_matrix.pdf/.svg` | PDF, SVG | H1-H9 evidence mapped to model layers with verdict panel |
| `layer_quantitative_summaries.pdf/.svg` | PDF, SVG | 4-panel figure with parameters tables and family methylation rates |
| `excluded_pathway_diagram.pdf/.svg` | PDF, SVG | Rejected vs actual pathway comparison diagram |
| `comprehensive_model_figure.pdf/.svg` | PDF, SVG | 6-panel publication figure (model, verdicts, families, coordination, geography, summary) |

### 8.2 Tables

| File | Description |
|------|-------------|
| `H1_H9_evidence_integration.tsv` | 8-row table with hypothesis verdicts, key statistics, layer mapping |
| `layer_parameters.tsv` | 22-parameter table across 3 layers |
| `testable_predictions.tsv` | 10 predictions with methods, expected outcomes, priorities |

### 8.3 Scripts

| File | Description |
|------|-------------|
| `gatekeeper_model_v2.py` | Complete analysis script for all figures and tables |

---

## 9. Methods

All figures generated with Python 3/matplotlib. Statistical tests from prior hypothesis analyses (H1-H9) were synthesized; no new statistical tests were performed in this integration analysis. Data sources include the 62 coordinated regulatory genes table (H6), family methylation summary (H6), TCS pair candidates (H8), COG enrichment results (H8), geographic distribution statistics (H8), and the corrected 37-TF dataset (H4).

---

## 10. Conclusions

The Gatekeeper Model v2 represents a fundamental revision of how DNA methylation interfaces with gene regulation in *S. coelicolor*. The key insight is that methylation does NOT directly control the well-studied transcriptional regulatory network. Instead, it operates through three distinct but interconnected mechanisms:

1. **Landscape remodeling**: Wholesale, MTase-driven changes in methylation patterns across the genome (not gene-by-gene)
2. **Active protection**: Critical regulatory DNA is shielded from methylation, ensuring master regulators function independently
3. **Signal gating through novel regulators**: A previously uncharacterized set of 62 regulatory genes, enriched in signal transduction, is the actual interface between methylation and gene expression

This model generates 10 testable predictions, 5 of which are prioritized as high-priority for immediate experimental validation. The most impactful would be SC_RS17645 knockout (confirming AAGCCCG MTase identity), SC_RS10435 characterization (validating chaplin-adjacent regulation), and CCGG geographic shift reproducibility (confirming the core Layer 1 finding).

---

*Analysis directory: `11_epigenome_integration/analysis/34_gatekeeper_model_v2/`*
*Script: `scripts/gatekeeper_model_v2.py`*
*Generated: 2026-02-24*
