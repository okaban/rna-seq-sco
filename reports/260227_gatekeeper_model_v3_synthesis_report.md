# Gatekeeper Model v3: Comprehensive Synthesis of M145 Epigenome-Transcriptome Integration

**Date:** 2026-02-27
**Scope:** Integration of 29 hypotheses across 14 exploration loops
**Status:** Final synthesis

---

## Executive Summary

In *Streptomyces coelicolor* A3(2) M145, DNA methylation (4mC/6mA) does **not** directly control transcription genome-wide. Instead, the methylation systems interact with the regulatory network through a "Shielded/Exposed" dichotomy: 998 of 1,055 regulatory genes maintain a ~1.2 kb protection zone that excludes methylation from the TSS (via collective promoter protein occupancy + evolutionary sequence divergence), rendering them methylation-insensitive. The remaining 57 regulatory genes lack this protection (methylation at median 114 bp from TSS), and these alone show methylation-expression coordination. This Shielded/Exposed dichotomy, defined by a sharp 293 bp boundary (AUC = 0.917), resolves the paradox of widespread methylation coexisting with minimal transcriptional impact.

---

## 1. The Central Finding: Shielded vs Exposed Regulatory Promoters

The single most important discovery of the 14-loop exploration is the identification of two structurally distinct classes within the 1,055 regulatory genes of *S. coelicolor* M145 (H27, Loop 13):

| Property | Shielded (n=998, 94.6%) | Exposed (n=57, 5.4%) |
|----------|:------------------------:|:---------------------:|
| Protection zone width | 1,200 bp | **0 bp** |
| Median nearest methylation site | 762 bp | **114 bp** |
| TSS methylation density | 0.37 sites/kb/gene | **3.15 sites/kb/gene** (8.4x enriched) |
| Methylation-expression correlation | rho = 0.063, p = 0.051 (NS) | **rho = 0.277, p = 0.029** |
| Constitutive expression | 38.0% | **0.0%** |
| \|log2FC\| T3vsT1 | 0.93 | **1.53** (1.7x larger) |
| Statistical separation | Mann-Whitney p = **5.3 x 10^-28** | 6.7x distance ratio |

**Interpretation:** The 998 shielded regulators are embedded in constitutively protected promoter environments where collective protein occupancy physically blocks methyltransferase (MTase) access. Because methylation cannot penetrate these promoters, genome-wide methylation changes have no effect on their transcription. The 57 exposed regulators lack this structural protection; methylation sites lie within or immediately adjacent to their promoters (-300 to +100 bp from TSS), placing them at positions where they could directly alter sigma factor recognition, TF binding affinity, or DNA melting at the transcription start site.

**Boundary definition (H29):** The optimal classifier between exposed and shielded is the nearest methylation-to-TSS distance, with a sharp threshold at **293 bp** (AUC = 0.917, sensitivity = 1.000, specificity = 0.806). No other feature -- including expression level (AUC = 0.547), genomic region, TF family, gene length, or number of predicted TF binding sites -- adds predictive power. The protection zone is an expression-independent, locus-specific structural property.

---

## 2. The Four-Layer Model

The Gatekeeper Model v3 organizes the epigenome-transcriptome relationship into four mechanistic layers, each supported by specific hypotheses and statistical evidence.

### Layer 1: R-M Defense Geographic Redistribution

**Key hypotheses:** H7 (CCGG/GCCGGC complete remodeling), H14 (defense island-GCCGGC coupling), H19 (Simpson's paradox)

**Core evidence:**

- GCCGGC 4mC sites undergo **complete positional remodeling** between timepoints (Jaccard similarity = 0.000 for all pairwise comparisons), with a dramatic geographic shift: T1 = 83% core (1,289 sites) --> T2 = 82% arm (407 sites) --> T3 = 62% arm (21 sites). The chi-squared statistic for this geographic redistribution is 597 (p ~ 10^-130).
- The defense island containing SC_RS36410 (Argonaute + PD-(D/E)XK nuclease + MTase) is coordinatedly induced at T3 (6-gene mean Spearman rho = 0.694), and one T3 GCCGGC site falls within the SC_RS36410 coding region (possible self-methylation).
- The GCCGGC R-M system responsible for the T1 core methylation program (1,289 sites) remains **unidentified** despite three independent approaches: annotation-based scoring (H10, SC_RS24685 ranked first but BLAST-negative), BLAST homology (H12, SC_RS19770 E = 0.007 but T1 expression = 11), and expression correlation (H16, SC_RS13615 top composite score but no DNA C-MTase domain).

**What this layer does NOT do:** The geographic redistribution does **not** regulate transcription. H17's apparent suppressive effect of GCCGGC proximity (p = 4.1 x 10^-10) was demonstrated by H19 to be a **Simpson's paradox**: T1 GCCGGC sites are core-biased (83%), core genes decline at T2 as part of normal development, and the apparent methylation-expression association vanishes upon geographic stratification (core-only: p = 0.87; arm-only: p = 0.019, reversed direction). The complete site remodeling with geographic redistribution is consistent with R-M defense biology (phage protection, self-DNA marking) rather than gene expression regulation.

### Layer 2: Regulatory DNA Protection (Three Sub-mechanisms)

**Overview:** All methylation systems in M145 -- regardless of modification type (4mC, 6mA) or recognition sequence (GCCGGC, AAGCCCG, CCGG) -- show significant depletion near regulatory genes and hypothetical protein genes (H15: fold = 0.43-0.71, p < 0.002 for all motifs). This depletion is a **genuine biological signal** surviving geographic stratification (H20: CMH-adjusted p < 0.007 for all motifs), not a Simpson's paradox artifact. The same stratification methodology that exposed H17/H19 as artifacts confirms H15/H20 as genuine.

This regulatory DNA protection operates through three complementary sub-mechanisms:

#### 2a: Evolutionary Counter-Selection (~25% contribution)

**Key hypothesis:** H22 (sequence-level motif depletion)

The R-M recognition motif DNA sequences (TGGCCGGC and AAGCCCG) are physically depleted from regulatory gene neighborhoods at the DNA level:

| Zone | TGGCCGGC fold | p-value | AAGCCCG fold | p-value |
|------|:---:|:---:|:---:|:---:|
| Gene body | **0.74** | 2.98 x 10^-7 | **0.76** | 1.26 x 10^-7 |
| Extended 2 kb | **0.92** | 6.79 x 10^-3 | **0.85** | 1.05 x 10^-3 |
| Promoter | 0.88 (NS) | 0.167 | 1.02 (NS) | 1.00 |

Critically, **gene bodies** show the strongest sequence depletion (fold = 0.74-0.76), while **promoters show no sequence depletion at all** (fold ~ 1.0). Comparing DNA-level depletion (fold = 0.85-0.92 at 2 kb) with methylation-level depletion (fold = 0.43-0.61 from H15), sequence counter-selection explains only **21-35%** of the total methylation avoidance. Permutation tests confirm the depletion is significant (p = 0.001 for AAGCCCG, p = 0.035 for TGGCCGGC at the extended 2 kb region).

This sub-mechanism represents a "hardwired" evolutionary baseline: over evolutionary time, R-M recognition motifs have been counter-selected from regulatory gene bodies (where coding sequence constraints permit), providing basal protection independent of transcriptional state.

#### 2b: Collective Promoter Occupancy Shield (~75% contribution)

**Key hypotheses:** H25 (TSS methylation gradient), H26 (TF BS-level test)

High-resolution spatial profiling (200 bp bins, +/-5 kb from TSS) reveals a **2,200 bp protection zone** around regulatory gene TSS (-1,300 to +700 bp) where methylation density is depleted by **17.3%** compared to flanking regions (vs. 3.8% for non-regulatory genes; a 4.6-fold difference in protection magnitude):

| Metric | Value |
|--------|-------|
| Protection zone width | 2,200 bp (-1,300 to +700 relative to TSS) |
| TSS depletion (regulatory) | 17.3% |
| TSS depletion (non-regulatory) | 3.8% |
| Deepest depletion position | +300 bp (reg/nonreg ratio = 0.541, p_adj = 0.005) |
| Core-specific depletion | 18.4% (core regulatory) vs -1.8% (arm non-regulatory) |

The deepest depletion at +300 bp downstream of TSS corresponds to the position of the RNAP-sigma factor transcription initiation complex, and the asymmetric upstream extension (-1,300 bp) reflects the footprint of multiple upstream-binding TFs at complex regulatory promoters.

Crucially, H26 demonstrated that **individual TF binding sites do NOT show methylation depletion**. GCCGGC 4mC is paradoxically *enriched* at predicted TF BS (fold = 1.157, p = 0.006), and spatial profiles around TF BS are completely flat (permutation p = 0.12-0.96). The protection zone is therefore a **collective, gene-level phenomenon**: no single TF binding event is sufficient to exclude MTase, but the cumulative occupancy of multiple TFs + RNAP + sigma factors at regulatory gene promoters creates an aggregate exclusion zone.

#### 2c: Shielded/Exposed Dichotomy

**Key hypotheses:** H27 (shielded vs exposed), H28 (exposed characteristics), H29 (293 bp boundary)

The protection zone is not a universal feature of all regulatory genes. Rather, it defines a binary dichotomy:

- **998 shielded regulators**: Full protection zone (1,200 bp width, 44.3% depth), methylation kept at median 762 bp from TSS, methylation-expression correlation non-significant.
- **57 exposed regulators**: No protection zone (0 bp width, 0.000 depth), methylation enriched 8.4x at TSS (median 114 bp from TSS), methylation-expression correlation significant (rho = 0.277, p = 0.029).

**What determines exposed status?** H28 found three distinguishing features: (1) lower T1 expression (median 90 vs 123, p = 0.010), (2) 100% dynamic expression with zero constitutive genes (OR = infinity, p = 5.3 x 10^-13), (3) higher promoter GC content (71.5% vs 69.7%, p = 2.0 x 10^-4). However, H29 definitively showed that **expression level does not predict classification** (baseMean AUC = 0.547, essentially random). The 293 bp boundary captures 98.9% of decision tree importance as the sole feature, and exposed regulators are distributed uniformly across all expression quintiles (4.4-7.9% per quintile, JT p = 0.730).

The **causal direction is reversed**: it is not "low expression causes loss of protection" but rather "loss of locus-specific structural protection allows methylation to approach TSS, which then contributes to expression variability." The protection zone is an expression-independent, locus-specific structural property, likely determined by sequence-level features (H22 evolutionary counter-selection), nucleoid-associated protein binding, or DNA topology -- not by transcriptional activity.

### Layer 3: Methylation-Responsive Regulatory Cascade (Limited to 62 Genes)

**Key hypotheses:** H6 (genome-wide regulator screen), H8 (coordinated regulator characterization), H27 (exposed promoter model), H28 (exposed characteristics)

The 57 exposed regulators constitute the sole pathway through which methylation dynamics could influence gene expression. Their characteristics:

- **Discovery (H6):** Genome-wide screen of 1,055 GFF-annotated regulatory genes identified 165 (15.6%) with methylation in their vicinity. Of these, 62 show coordinated methylation-expression changes. All 62 are outside the literature TF list of 37 -- the literature list was biased toward TF cascade apex genes (SARPs, master regulators), which are universally shielded.
- **Coordination types (H8):** discordant_gain_up (22), concordant_derepression (12), concordant_repression (12), discordant_loss_down (16). All four types show identical absence of protection zones and similar TSS-proximal methylation.
- **Functional enrichment (H8):** COG T (signal transduction) is enriched (OR = 2.47, p = 0.045). COG Q (secondary metabolism) = 0. The exposed regulators are regulatory layer components, not biosynthetic genes. Seven TCS cognate pairs were identified, all showing asymmetric methylation (one partner methylated, the other not).
- **Geographic distribution (H8, H27):** No significant arm/core bias overall (OR = 1.26, p = 0.42), except that the T3 discordant_gain_up subgroup shows strong arm enrichment (87%, OR = 8.05, p = 0.001).

**Critical caveat:** While the 57 exposed regulators show methylation-expression coordination and represent plausible targets of methylation-based regulation, the evidence remains correlational. The temporal causality tests (H19, H21) failed to demonstrate that methylation loss causes expression changes for any methylation system. The coordination could reflect shared responses to developmental signals rather than a causal methylation-to-expression pathway.

---

## 3. What Methylation Does NOT Do (Negative Results Equally Important)

The negative findings of this project are as significant as the positive ones:

### 3.1 Does not directly control individual gene transcription (H5, H19, H21)

- **H5:** Gene-level methylation-expression correlation is weak across all motifs (|rho| < 0.16, non-significant). Methylation operates as a global/landscape-level phenomenon, not a gene-specific cis-regulatory element.
- **H19 (GCCGGC):** Genes losing methylation at T2 show no de-repression. Instead, they show significant *downregulation* (median LFC = -0.253 vs Never = +0.006), entirely explained by geographic confounding (core-only: p = 0.87).
- **H21 (AAGCCCG):** The cleanest temporal test (minimal geographic shift, 4.7 pp). Lost vs Never: p = 0.91, r = 0.003. No de-repression at all.

### 3.2 Simpson's paradox: apparent effects are geographic confounds (H19, H24)

- **H17 to H19:** The cross-sectional "suppressive effect" of GCCGGC proximity (p = 4.1 x 10^-10, H17) was unmasked as a Simpson's paradox. T1 sites are 83% core-located; core genes have lower LFC at T2 independent of methylation. Geographic stratification completely eliminates the signal.
- **H24:** BGC methylation enrichment (fold = 1.66, H23) is explained by (a) 100% core co-localization of BGC genes and (b) GC-rich PKS/NRPS coding sequences containing more GCCGGC motifs at the DNA level. No active methylation targeting.

### 3.3 Expression level does not determine protection (H29)

- baseMean AUC = 0.547 (near random) for predicting exposed/shielded status
- No trend across expression quintiles (JT p = 0.730)
- RNAP occupancy model is rejected; protection is a structural property

### 3.4 Individual TF binding does not exclude methylation (H26)

- TF BS spatial profiles are flat (permutation p = 0.12-0.96)
- GCCGGC is enriched at TF BS (fold = 1.157, p = 0.006)
- Protection requires collective promoter-wide occupancy, not individual TF binding

### 3.5 TF cascade apex regulators are not methylation-controlled (H4)

- 37 literature TFs: only 3/37 (8%) have nearby methylation, zero show coordinated changes
- All four SARPs (redD, actII-ORF4, cdaR, cpkO) have zero methylation
- The canonical regulatory cascade operates independently of methylation

---

## 4. Complete Hypothesis Ledger

| ID | Hypothesis | Loop | Verdict | One-line Summary |
|----|-----------|:----:|---------|-----------------|
| H1 | 4mC is promoter-biased, 6mA gene-body-dispersed, causing differential expression correlation | 1 | **Partial** | Distribution differs (p = 1.56 x 10^-14) but 6mA is more promoter-proximal; 4mC promoter correlates with expression (rho = 0.139, p = 8.3 x 10^-4) |
| H2 | MTase expression stability causes 4mC/6mA correlation difference | 1 | Not tested | Deferred to H5 |
| H3 | absA de-repression compensates for redZ to activate Red BGC | 1 | **Premise rejected** | True redZ (SC_RS31650) is upregulated (LFC = +0.91); paradox does not exist; 36/37 literature locus_tag errors discovered |
| H4 | Corrected locus_tags reveal TF methylation-expression coordination | 2 | **Rejected** | 3/37 TFs methylated, all constitutive; 0/37 show coordination; 0/4 SARPs methylated |
| H5 | Stable MTase expression yields stronger gene-level methylation-expression correlation | 2 | **Partial** | Site-level dynamics confirmed (AAGCCCG 260 sites lost with SC_RS17645 downregulation), but gene-level correlations are weak (all |rho| < 0.16) |
| H6 | Genome-wide regulator screen finds coordinated genes outside literature 37 | 3 | **Supported** | 1,055 regulators from GFF; 62 show coordination, all outside literature list; MerR highest methylation rate (28.6%) |
| H7 | CCGG 4mC decrease is passive dilution from Dcm-like MTase expression loss | 3 | **Rejected** | Complete non-overlap (Jaccard = 0.000); T1 = 77% core to T2 = 89% arm; active remodeling, not passive dilution |
| H8 | 57 coordinated regulators are arm-enriched and secondary metabolism-associated | 4 | **Partial** | Overall arm enrichment NS (p = 0.17), but T3 discordant_gain_up = 87% arm (OR = 8.05, p = 0.001); COG T enriched (OR = 2.47); 7 TCS pairs with asymmetric methylation |
| H9 | CCGG "4mC" is Dcm-derived 5mC misclassified by Nanopore | 4 | **Rejected** | 5mC signal = 0.01%; 4mC frequency highest of all modifications (82.7%); 91.4% in GCCGGC palindrome, not CCWGG |
| H10 | REBASE + M145 MTase catalog identifies GCCGGC N4-C MTase | 5 | **Partial** | Only 2/87 NaeI-family produce m4C; SC_RS24685 top by annotation (9/10) but later BLAST-negative (H12) |
| H11 | Gatekeeper Model v2 integrates H1-H9 into 3-layer quantitative model | 5 | **Completed** | 3-layer model built; 10 testable predictions defined; superseded by v3 (this report) |
| H12 | M.Svi27968I BLAST identifies GCCGGC MTase as SC_RS24685 | 6 | **Rejected** | SC_RS19770 is BLAST top hit (E = 0.007, 29% identity); SC_RS24685 is out of range (E > 10); T1 enzyme remains unidentified |
| H13 | AAGCCCG sites are non-randomly distributed, enriched near secondary metabolism | 6 | **Partial** | Non-random distribution confirmed; Regulatory/TF depleted (fold = 0.43, p = 0.035); Hypothetical depleted (fold = 0.52, p = 0.039); secondary metabolism NOT enriched |
| H14 | SC_RS36410 defense island activation produces T3 arm-enriched GCCGGC methylation | 7 | **Partial** | Geographic shift confirmed (chi2 = 597, p ~ 10^-130); 1 T3 site in defense island; 6-gene coordinated induction; but T3 = only 21 sites total |
| H15 | GCCGGC/CCGG 4mC also show regulatory gene avoidance (like AAGCCCG) | 7 | **Supported** | Universal avoidance: AAGCCCG fold = 0.43, GCCGGC fold = 0.61, All 4mC fold = 0.62 (all p < 0.02); Hypothetical also depleted |
| H16 | Expression dynamics reverse-identify GCCGGC MTase | 8 | **Partial** | SC_RS13615 (SCO2317) top composite score; COG0863 DNA methylase domain; but no candidate satisfies all four criteria simultaneously |
| H17 | GCCGGC and AAGCCCG show motif-specific transcriptional "division of labor" | 8 | **Rejected** | GCCGGC proximity associated with suppression (p = 4.1 x 10^-10) but this is later shown to be geographic confound (H19) |
| H18 | GCCGGC site dose-response on nearby gene expression | 9 | **Partial** | Threshold switch (0 vs >=1), not graded dose-response; within methylated genes rho = -0.037 (T2); arm/core stratification eliminates signal |
| H19 | Genes losing GCCGGC methylation at T2 show de-repression (upregulation) | 9 | **Rejected** | Lost genes show downregulation (median LFC = -0.253); Simpson's paradox identified; core-only p = 0.87; AAGCCCG reciprocal p = 0.92 |
| H20 | H15's regulatory avoidance survives geographic stratification | 10 | **Supported** | Core-only and arm-only both significant (p < 0.015); CMH-adjusted p < 0.007 for all motifs; arm ORs often lower (stronger depletion) |
| H21 | AAGCCCG loss causes temporal de-repression (clean test, minimal geographic shift) | 10 | **Rejected** | Lost vs Never: p = 0.91, r = 0.003; no de-repression; no geographic confound either (clean null) |
| H22 | R-M motif sequences are depleted from regulatory gene DNA | 11 | **Partial** | Gene body depleted (fold = 0.74-0.76, p < 3 x 10^-7); promoter NOT depleted (fold ~ 1.0); sequence explains 21-35% of methylation avoidance |
| H23 | Methylation avoidance extends to all essential gene categories | 11 | **Rejected** | Only regulatory genes show avoidance; ribosomal proteins, DNA repair, cell division, energy metabolism show no depletion; BGC genes enriched |
| H24 | BGC methylation enrichment is a geographic confound | 12 | **Partial** | All 4mC enrichment disappears in core-only; GCCGGC enrichment partially remains but explained by DNA sequence composition (GC-rich PKS/NRPS) |
| H25 | Methylation density gradient around regulatory TSS is steeper | 12 | **Supported** | 2,200 bp protection zone; 17.3% depletion (vs 3.8%); deepest at +300 bp; all motifs independently replicate |
| H26 | Individual TF binding sites show methylation depletion | 13 | **Rejected** | Flat spatial profiles (p = 0.12-0.96); GCCGGC enriched at TF BS (fold = 1.16); protection is collective not individual |
| H27 | 57 coordinated regulators have shallower protection zones | 13 | **Supported (opposite)** | Not shallower but **absent**; 8.4x enrichment at TSS; 6.7x closer (p = 5.3 x 10^-28); Shielded/Exposed dichotomy |
| H28 | Exposed regulators have low T1 expression (RNAP not occupying promoter) | 14 | **Supported** | T1 lower (p = 0.010); 100% dynamic (p = 5.3 x 10^-13); higher GC (p = 2.0 x 10^-4); but H29 reverses causal interpretation |
| H29 | baseMean threshold separates shielded from exposed (AUC > 0.7) | 14 | **Rejected** | baseMean AUC = 0.547 (random); nearest_methyl_distance AUC = 0.917 at 293 bp; protection is expression-independent |

**Verdict distribution across 28 testable hypotheses:**

| Verdict | Count | Percentage |
|---------|:-----:|:----------:|
| Supported | 5 | 17.9% |
| Supported (opposite direction) | 1 | 3.6% |
| Partial | 10 | 35.7% |
| Rejected | 10 | 35.7% |
| Premise rejected | 1 | 3.6% |
| Not tested | 1 | 3.6% |

---

## 5. Key Figures for Publication

The following figures represent the strongest candidates for a research publication:

### Main Figure (Multi-panel)

1. **TSS methylation spatial profile** (H25): Regulatory vs non-regulatory density with smoothed CI ribbons showing the 2,200 bp protection zone. The deepest depletion at +300 bp from TSS is the spatial fingerprint of promoter occupancy shielding.
   - Source: `48_TSS_methylation_gradient/figures/H25_comprehensive_summary.pdf`

2. **Shielded vs Exposed protection zone comparison** (H27): Side-by-side density profiles for 998 shielded vs 57 exposed regulators, showing the complete absence of protection and the 8.4x TSS enrichment spike.
   - Source: `50_coordinated_regulators_protection/figures/H27_comprehensive_summary.pdf`

3. **293 bp boundary ROC and distance distributions** (H29): ROC curve showing AUC = 0.917 for nearest_methyl_distance, with violin/box plots of the bimodal distance distribution.
   - Source: `52_shielded_exposed_boundary/figures/H29_comprehensive_summary.pdf`

### Supporting Figures

4. **Simpson's paradox demonstration** (H19): Unstratified vs geographic-stratified LFC distributions for Lost/Gained/Never groups, illustrating how the apparent methylation-expression association dissolves under stratification.
   - Source: `42_GCCGGC_temporal_derepression/figures/H19_comprehensive_summary.pdf`

5. **Genuine regulatory avoidance vs geographic artifact** (H20): Forest plot of stratified odds ratios showing H15's regulatory avoidance survives the same methodology that exposed H17/H19 as artifacts.
   - Source: `43_regulatory_avoidance_geographic_test/figures/H20_comprehensive_summary.pdf`

6. **DNA vs methylation depletion** (H22): Paired comparison showing sequence-level (21-35%) vs protein occupancy (65-79%) contributions to regulatory gene methylation avoidance.
   - Source: `45_sequence_level_motif_depletion/figures/H22_comprehensive_summary.png`

7. **Exposed regulators: 100% dynamic expression** (H28): Temporal pattern distribution showing zero constitutive expression among exposed regulators (OR = infinity, p = 5.3 x 10^-13).
   - Source: `51_exposed_regulators_characteristics/figures/H28_comprehensive_summary.pdf`

8. **GCCGGC geographic redistribution** (H14): Chromosome ideogram showing T1 core (83%) to T2 arm (82%) to T3 arm (57%) shift.
   - Source: `37_defense_island_GCCGGC/figures/H14_defense_island_GCCGGC_analysis.pdf`

---

## 6. Unresolved Questions and Future Directions

### 6.1 Identity of the T1 GCCGGC MTase

The enzyme responsible for the 1,289 GCCGGC 4mC sites at T1 (83% core) remains unidentified despite three independent approaches:

- Annotation-based scoring (H10): SC_RS24685 top-ranked but BLAST-negative
- BLAST homology (H12): SC_RS19770 top hit (E = 0.007) but T1 expression = 11
- Expression correlation (H16): SC_RS13615 (SCO2317) top composite score but no DNA C-MTase domain match

GCCGGC m4C production is globally rare (2/87 NaeI-family enzymes in REBASE) and exceptional within *Streptomyces* (1/28 GCCGGC R-M systems). Experimental approaches (targeted gene KO + SMRT-seq) are required.

### 6.2 Structural determinants of the 293 bp boundary

With expression level rejected as a predictor (AUC = 0.547), what defines the locus-specific protection? Candidates include:

- DNA sequence features beyond R-M motif density (H22's 25% contribution may underestimate structural effects)
- Nucleoid-associated protein (NAP) binding sites (HupA, IHF, Lsr2 homologs)
- DNA supercoiling topology and replication-dependent accessibility
- Long-range chromosomal organization (Hi-C data would be informative)

### 6.3 Downstream targets of the 57 exposed regulators

The 57 exposed regulators are enriched for signal transduction (COG T, OR = 2.47) and include 7 TCS pairs with asymmetric methylation. Their regulatory targets remain uncharacterized. Integration with the FIMO binding site predictions (56,338 hits genome-wide) could map the indirect regulatory cascade from methylation to expression.

### 6.4 Cross-species conservation of the model

Key questions:

- Do other *Streptomyces* species show the same Shielded/Exposed dichotomy?
- Is the ~293 bp boundary conserved across Actinobacteria?
- Does the AAGCCCG system (completely novel in the genus; 0/82 REBASE species) represent a recent horizontal acquisition?

### 6.5 Causal validation

While the Shielded/Exposed model provides a compelling mechanistic framework, the temporal causality tests (H19, H21) failed for both major methylation systems. Key experiments:

- MTase knockout strains to test whether loss of methylation at specific exposed promoters alters expression
- In vitro methylation of exposed promoter sequences to test direct effects on transcription factor binding
- ChIP-seq for nucleoid-associated proteins to map protection zone determinants

---

## 7. Methodological Lessons

### 7.1 Simpson's paradox and geographic stratification

The most impactful methodological lesson is the systematic application of geographic stratification to distinguish genuine signals from geographic confounds in the *S. coelicolor* linear chromosome:

| Signal | Unstratified p-value | Stratified result | Verdict |
|--------|:---:|:---:|:---:|
| GCCGGC suppression (H17) | 4.1 x 10^-10 | Core: p = 0.87 | **ARTIFACT** |
| GCCGGC de-repression (H19) | 8.3 x 10^-8 | Core: p = 0.87 | **ARTIFACT** |
| BGC enrichment All 4mC (H24) | 4.4 x 10^-8 | Core-only: fold = 1.07 (NS) | **ARTIFACT** |
| Regulatory avoidance (H20) | 1.4 x 10^-7 | Core: p = 3.8 x 10^-5; Arm: p = 7.3 x 10^-3 | **GENUINE** |
| AAGCCCG de-repression (H21) | p = 0.91 | Core: p = 0.28; Arm: p = 0.23 | **GENUINE NULL** |

The linear chromosome of *S. coelicolor*, with its core/arm functional differentiation, creates systematic confounds whenever a feature is geographically biased (e.g., methylation sites: 83% core at T1). Any apparent methylation-expression association MUST be validated by geographic stratification before causal interpretation. The Cochran-Mantel-Haenszel test is the recommended approach.

### 7.2 Site-centric vs gene-centric analysis

H15 (site-centric) found strong regulatory avoidance (fold = 0.43-0.61), while H23 (gene-centric) found fold = 1.01 (NS). This discrepancy is not contradictory: site-centric analysis measures methylation density around genes, while gene-centric analysis measures binary proximity (has/lacks a site within 2 kb). Regulatory genes show reduced methylation density but are not completely excluded from methylation (46.9% have at least one GCCGGC site within 2 kb). The protection mechanism is probabilistic density modulation, not all-or-nothing exclusion.

### 7.3 Importance of null results

Of 28 testable hypotheses, 10 were rejected and 1 premise-rejected (39.3%). These null results were not failures but critical advances:

- **H4 rejection** (TF cascade not methylation-controlled) redirected the search from literature TFs to genome-wide regulators, leading to the 57 exposed regulators (H6)
- **H19 rejection** (Simpson's paradox) protected against false positive claims of methylation-mediated transcription control
- **H26 rejection** (individual TF BS not protected) clarified that protection is a collective phenomenon, not attributable to individual binding events
- **H29 rejection** (expression does not predict protection) reversed the causal interpretation from H28, establishing protection as a structural property

---

## Appendix: Integrated Gatekeeper Model v3 Schematic

```
S. coelicolor M145 Epigenome-Transcriptome Model (v3)

LAYER 1: R-M DEFENSE GEOGRAPHIC REDISTRIBUTION
[H7, H14, H19]
- GCCGGC 4mC: T1=1,289 sites (83% core) --> T2=407 (82% arm) --> T3=21 (57% arm)
- AAGCCCG 6mA: T1=260 (69% core) --> T2=64 (64% core) [more stable]
- Complete positional remodeling (Jaccard = 0.000)
- Function: Phage defense, self-DNA marking
- NOT transcriptional regulation (H19 Simpson's paradox)

              |
              v

LAYER 2: REGULATORY DNA PROTECTION
[H15, H20, H22, H25, H26]

  2a. Evolutionary counter-selection (~25%)     2b. Collective promoter occupancy (~75%)
  [H22]                                         [H25, H26]
  - Gene body motif depletion                   - 2,200 bp protection zone at TSS
    (fold = 0.74-0.76)                          - 17.3% depletion (4.6x vs non-reg)
  - Promoter: NO depletion (fold ~ 1.0)         - +300 bp deepest (RNAP position)
  - Hardwired basal protection                  - No individual TF BS depletion (H26)
                                                - Collective multi-protein effect

              |
              v

  2c. SHIELDED / EXPOSED DICHOTOMY [H27, H29]
  ============================================
  |                                            |
  |    SHIELDED (998 genes, 94.6%)             |    EXPOSED (57 genes, 5.4%)
  |    Protection zone: 1,200 bp               |    Protection zone: 0 bp
  |    Nearest methyl: 762 bp                  |    Nearest methyl: 114 bp
  |    Methylation-insensitive                 |    8.4x TSS enrichment
  |    38% constitutive expression             |    0% constitutive expression
  |    --> NO methylation effect               |    --> Methylation-responsive
  |                                            |
  ==============================================
             |                                          |
             v                                          v

EXCLUDED                                   LAYER 3: METHYLATION-RESPONSIVE
(NO PATHWAY)                               REGULATORY CASCADE (57 genes)
                                           [H6, H8, H27, H28]
- Literature TF cascade (H4: 0/37)         - Signal transduction enriched (COG T, p=0.045)
- Direct transcriptional control            - 7 TCS pairs (asymmetric methylation)
  (H5, H19, H21: all null)                 - 100% dynamic expression
- Dose-dependent repression                 - Methylation-expression rho = 0.277
  (H18: threshold only, confounded)         - Causal link: UNPROVEN (H19/H21 null)
                                            - 293 bp boundary (AUC = 0.917)
```

---

*Analysis performed: 2026-02-24 through 2026-02-27 (14 exploration loops)*
*This report integrates findings from 29 hypotheses tested across 52 analysis directories.*
*Scripts and data: `/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/`*
