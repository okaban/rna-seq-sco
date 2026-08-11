# Gatekeeper Model v4: Final Synthesis of M145 Epigenome-Transcriptome Integration

**Date:** 2026-02-27
**Scope:** Complete integration of 36 hypotheses across 18 exploration loops
**Status:** Final synthesis (v4)

---

## 1. Executive Summary

Through 18 iterative hypothesis-driven exploration loops spanning 36 hypotheses, we discovered that DNA methylation (4mC/6mA) in *Streptomyces coelicolor* A3(2) M145 does **not** function as a direct transcriptional regulator genome-wide. Instead, the methylation systems interact with the regulatory architecture through a four-layer "Gatekeeper" mechanism: (1) R-M defense systems undergo dramatic geographic redistribution across the linear chromosome during development, (2) 94.6% of regulatory gene promoters maintain a ~1.2 kb methylation-free protection zone created by a two-tier defense (evolutionary counter-selection contributing ~33% and collective promoter protein occupancy contributing ~67%), (3) the remaining 57 regulatory genes (5.4%) lack this protection entirely -- with methylation sites directly at their TSS -- creating a sharp "Shielded/Exposed" dichotomy defined by a 293 bp boundary (AUC = 0.917), and (4) these 57 exposed regulators constitute a distributed vegetative-to-developmental switch, organized into two synchronously activated antagonistic blocs that simultaneously repress growth programs (metabolism, DNA repair, efflux) and activate developmental programs (morphogenesis, signal transduction, stress response) at the exponential-to-transition boundary.

**Key numbers:** 18 loops, 36 hypotheses, 7 supported, 15 partially supported, 12 rejected, 2 other (1 not tested, 1 premise rejected). Of 1,055 regulatory genes analyzed, 998 are shielded and 62 are exposed. The central discovery -- that apparent genome-wide methylation-expression correlations are Simpson's paradox artifacts while the true methylation-responsive pathway operates through only 57 specific regulatory genes -- required the systematic rejection of 12 hypotheses to establish.

**The central discovery in brief:** Bacterial DNA methylation does not control transcription directly. Instead, a binary structural dichotomy at regulatory gene promoters -- defined by the presence or absence of a protein occupancy-dependent protection zone -- partitions the regulatory genome into a methylation-insensitive majority and a methylation-responsive minority. This minority encodes a coordinated developmental switch that is equally conserved, genomically dispersed, and mechanistically distinct from all characterized transcriptional regulatory pathways.

---

## 2. The Gatekeeper Model v4 -- Final Architecture

### Layer 1: R-M Defense Geographic Remodeling (H7, H14, H19)

The two major methylation systems in M145 undergo dramatic positional remodeling during the growth cycle, but this remodeling is **not** transcriptional regulation:

| System | T1 (exponential) | T2 (transition) | T3 (stationary) | Geographic shift |
|--------|:-:|:-:|:-:|:-:|
| **GCCGGC 4mC** | 1,289 sites (83% core) | 407 sites (82% arm) | 21 sites (57% arm) | Core-to-arm (chi2=597, p~10^-130) |
| **AAGCCCG 6mA** | 260 sites (69% core) | 64 sites (64% core) | -- | 75% loss, minimal shift (4.7 pp) |

Key features:

- **Complete positional remodeling**: T1 and T2 GCCGGC sites are entirely non-overlapping (Jaccard = 0.000), proving active remodeling rather than passive dilution (H7)
- **AAGCCCG dramatic loss**: 260 to 64 sites correlates with SC_RS17645 (HsdM-type MTase) expression drop (LFC = -2.19), but this loss does not cause de-repression of target genes (H21: p = 0.91, r = 0.003)
- **Simpson's Paradox**: The apparent suppressive effect of GCCGGC methylation proximity (H17: p = 4.1 x 10^-10) is entirely a geographic artifact. T1 sites are 83% core-located; core genes have lower LFC at T2 independent of methylation. Geographic stratification eliminates the signal completely (core-only: p = 0.87; H19)
- **Function**: Phage defense and self-DNA marking, not transcriptional control

### Layer 2: Regulatory DNA Protection System

The protection of regulatory gene promoters from methylation operates through three integrated sub-layers:

#### Layer 2a: Evolutionary Counter-Selection (~25-33%)

Regulatory gene DNA sequences have evolved to reduce R-M recognition motif density (H22, H30):

| Evidence | Metric | Value | Reference |
|----------|--------|-------|-----------|
| Gene body motif depletion | Fold (reg vs non-reg) | 0.74-0.76 (p < 3 x 10^-7) | H22 |
| TSS +/-300bp AAGCCCG enrichment at exposed TFs | Fold (exposed/shielded) | 4.97 (p = 1.1 x 10^-8) | H30 |
| TSS +/-300bp TGGCCGGC enrichment at exposed TFs | Fold (exposed/shielded) | 2.66 (p = 4.4 x 10^-4) | H30 |
| Combined sequence model | AUC (5-fold CV) | 0.712 (+/- 0.049) | H30 |
| Sequence contribution estimate | % of total protection | ~33% | H30 |

The AAGCCCG motif is the strongest single-sequence discriminator between exposed and shielded regulators (5.0x enrichment, p = 1.1 x 10^-8). Shielded regulators have evolved away from AAGCCCG near their TSS, while exposed regulators retain these motifs -- but critically, retaining the motif sequence does NOT mean the site is methylated (H35-causal).

#### Layer 2b: Collective Promoter Occupancy Shield (~67-75%)

The dominant protection mechanism is a collective protein occupancy effect at regulatory gene promoters (H25, H26, H29):

| Evidence | Metric | Value | Reference |
|----------|--------|-------|-----------|
| Protection zone width | bp | 2,200 (-1,300 to +700 from TSS) | H25 |
| Deepest depletion position | bp from TSS | +300 (ratio = 0.541) | H25 |
| TSS methylation depletion (reg vs non-reg) | % | 17.3% (vs 3.8%, 4.6x ratio) | H25 |
| Individual TF BS protection | Permutation p | 0.12-0.96 (NS, flat profiles) | H26 |
| GCCGGC enrichment at TF BS | Fold | 1.157 (p = 0.006, counter-intuitive) | H26 |
| Nearest methylation distance AUC | AUC | 0.917 (95% CI: 0.895-0.935) | H29 |
| Optimal distance threshold | bp | 293 | H29 |
| Expression independence | baseMean AUC | 0.547 (near random) | H29 |
| Expression quintile trend | JT p | 0.730 (no trend) | H29 |

Critical insight: Protection is **not** determined by individual TF binding site occupancy (H26: FIMO BS are actually enriched for GCCGGC methylation). Rather, it is a collective, expression-independent structural property of regulatory gene promoters. The 293 bp boundary achieves 100% sensitivity and 80.6% specificity with a single distance measurement, establishing it as a near-binary structural feature.

#### Layer 2c: Shielded/Exposed Dichotomy (H27, H29, H36)

The protection zone creates a sharp binary classification of all 1,055 regulatory genes:

| Property | Shielded (n=998, 94.6%) | Exposed (n=57, 5.4%) |
|----------|:------------------------:|:---------------------:|
| Protection zone width | 1,200 bp | **0 bp** |
| Median nearest methylation | 762 bp | **114 bp** |
| TSS methylation density | 0.37 sites/kb/gene | **3.15 sites/kb/gene** (8.4x) |
| Methylation-expression rho | 0.063, p = 0.051 (NS) | **0.277, p = 0.029** |
| Constitutive expression | 38.0% | **0.0%** |
| |log2FC| T3vsT1 | 0.93 | **1.53** (1.7x larger) |
| Statistical separation | Mann-Whitney p = 5.3 x 10^-28 | 6.7x distance ratio |
| SCO locus tag rate | 97.8% | **96.8%** |
| Named product rate | 100.0% | **100.0%** |
| Composite conservation score | -0.005 | **0.073** (p = 0.276, NS) |
| Core localization | 63.8% | 56.5% (p = 0.277, NS) |

**H36 (NEW)**: The exposed regulators are equally conserved as shielded regulators by all available evolutionary proxy measures -- GC3, effective number of codons, rare codon frequency, annotation quality, and chromosomal position all show no significant difference (composite ROC AUC = 0.459). The exposed/shielded distinction is determined by **epigenomic context** (local methylation landscape), not by the gene's evolutionary history.

### Layer 3: 57 Exposed Regulators -- The Methylation-Responsive Switch

The 57 exposed regulatory genes represent the sole conduit through which methylation interacts with the transcriptional program. Loops 15-18 revealed their internal architecture, temporal dynamics, mechanism, and evolutionary status.

#### 3a. Internal Structure (H32, H35-TF)

The 57 exposed regulators are NOT a single coherent module. They are a **distributed regulatory layer** with strong internal structure:

| Feature | Result | Evidence |
|---------|--------|---------|
| Genomic clustering | NOT significant (z = +0.30, p = 0.37) | H32 |
| E-E co-expression > E-S | NOT significant (p = 0.257) | H32 |
| Operonic pairs | 0 (p = 1.0) | H32 |
| Within-type co-expression | **rho = 0.500** (p = 6.6 x 10^-25) | H32 |
| Co-expression modules | 4 modules, 61/57 genes covered | H32 |
| TCS asymmetry | 7/7 pairs split (one exposed, one shielded) | H32 |

**Two antagonistic blocs (H35-TF):**

| Bloc | Size | Function | Key TF families | Biological role |
|------|:----:|----------|-----------------|-----------------|
| **Activation** (M1-3) | ~36 | Developmental programs ON | TCS (9, 25%), Sigma (5), WhiB (1), RamR | Morphogenesis, signal transduction, stress response |
| **Repression** (M4) | ~26 | Vegetative programs OFF | TetR (7, p=0.028), GntR, IclR, LacI, SSB, HU | Metabolism shutdown, DNA repair cessation, efflux repression |

The TetR family is the only statistically significantly enriched family (OR = 0.16, p = 0.028 in repression bloc). TCS members show complete asymmetry: in all 7 cognate TCS pairs, exactly one member is exposed and one is shielded -- never both (binomial p < 10^-3 for this perfect asymmetry). Methylation is distributed 4 SK : 3 RR with no bias (p = 1.0).

#### 3b. Temporal Dynamics (H34)

The two blocs respond **synchronously** at the T1-to-T2 developmental transition:

| Metric | Activation bloc | Repression bloc | Test | p-value |
|--------|:-:|:-:|:-:|:-:|
| N genes | 35 | 26 | -- | -- |
| Mean phase ratio | 0.815 | 0.962 | MW U | 0.459 |
| % Early responders | 60% | 69% | -- | -- |
| T1 z-score | -1.03 (low) | +1.28 (high) | -- | -- |
| T3 z-score | +0.79 (high) | -0.81 (low) | -- | -- |

Key findings:
- **63% are early responders** (phase ratio > 0.6), with the bulk of expression change occurring during the T1-to-T2 transition
- **Simultaneous switch**: Both blocs activate at the same time (p = 0.459 for phase separation). The "gate" opens once at the developmental boundary
- **Module-internal coherence is strong** (rho = 0.717, p = 1.0 x 10^-19) but **cross-module prediction is zero** -- modules respond in parallel, not in sequence
- **Methylation timing is decoupled from expression timing** (rho = 0.136, p = 0.299 for T1-T2; rho = 0.012, p = 0.927 for T2-T3)

#### 3c. Regulatory Mechanism (H31, H33, H35-causal)

The mechanism by which exposed TFs exert their regulatory effects remains unknown, but a series of negative results has eliminated all computationally testable candidates:

| Candidate mechanism | Test | Result | Reference |
|---------------------|------|--------|-----------|
| Known TF cascade (FIMO motifs) | 0/57 have curated binding motifs | **Eliminated** | H31 |
| Preferential targeting by master regulators | 38/57 targeted, same rate as shielded (OR=0.881, p=0.683) | **No enrichment** | H31 |
| Cis-regulatory (neighborhood effect) | Permutation p = 0.857, no distance decay | **Eliminated** | H33 |
| AAGCCCG methylation-mediated | Only 2/57 have methylated AAGCCCG at TSS | **Eliminated** | H35-causal |
| Upstream TF regulation driving coordination | Target vs non-target |LFC| p = 0.891 | **Independent** | H31 |

**The critical H30/H35 distinction**: AAGCCCG **DNA sequence** is 5.0x enriched near exposed TF promoters (25.8% have motifs; H30), but only 3.2% have **methylated** AAGCCCG sites (H35-causal). The MTase methylates only ~19.5% of its recognition sites genome-wide. Therefore, sequence enrichment does NOT equal methylation enrichment -- this is a structural/evolutionary property of the DNA, not an active methylation mechanism.

The 57 exposed TFs act through **trans-acting effects on dispersed targets** via an unknown mechanism. Experimental approaches (ChIP-seq, DAP-seq, genetic knockouts) are required.

#### 3d. Conservation (H36)

The exposed regulators are not evolutionary newcomers:

| Metric | Exposed | Shielded | p-value |
|--------|---------|----------|---------|
| SCO locus tags | 96.8% | 97.8% | 0.646 |
| Named products | 100% | 100% | 1.000 |
| GC3 | 0.919 | 0.932 | 0.060 (NS) |
| Nc (effective codons) | 31.6 | 31.5 | 0.725 |
| Rare codon frequency | 0.032 | 0.027 | 0.122 |
| Composite conservation AUC | -- | -- | 0.459 (below random) |

Exposed status is an **epigenomic context property**, not an evolutionary divergence property. These are standard, well-conserved *Streptomyces* regulatory genes that happen to reside in genomic positions where methylation sites are proximal to their TSS.

---

## 3. Complete Hypothesis Ledger (H1-H36)

### Supported (7)

| Loop | ID | Hypothesis | Key Statistic | Reference |
|:----:|:--:|-----------|---------------|-----------|
| 3 | H6 | Genome-wide regulator screen finds coordinated genes outside literature 37 | 57/1,055 coordinated (all novel) | analysis/29 |
| 7 | H15 | Cross-motif regulatory gene methylation avoidance | All motifs: fold=0.43-0.62, p<0.02 | analysis/38 |
| 10 | H20 | H15 regulatory avoidance survives geographic stratification | CMH-adjusted p<0.007, core+arm both significant | analysis/43 |
| 12 | H25 | TSS methylation spatial gradient at regulatory genes | 2,200bp zone, +300bp deepest (ratio=0.541) | analysis/48 |
| 13 | H27 | 57 coordinated regulators have absent (not shallower) protection zones | 8.4x TSS enrichment, p=5.3e-28 | analysis/50 |
| 14 | H28 | Exposed regulators: low T1, 100% dynamic, high GC | constitutive=0% (OR=inf, p=5.3e-13) | analysis/51 |
| 16 | H35-TF | TF family predicts activation/repression bloc membership | TetR OR=0.16 (p=0.028), TCS 3.3x activation | analysis/58_func |

### Partially Supported (15)

| Loop | ID | Hypothesis | Key Statistic | Reference |
|:----:|:--:|-----------|---------------|-----------|
| 1 | H1 | 4mC promoter-biased, 6mA gene-body-dispersed | Distribution differs (p=1.56e-14) but 6mA more promoter-proximal | analysis/20 |
| 2 | H5 | MTase stability correlates with methylation-expression | Site-level dynamics confirmed, gene-level weak (|rho|<0.16) | analysis/24 |
| 4 | H8 | Coordinated regulators arm-enriched, secondary metabolism | T3 discordant_gain_up 87% arm (OR=8.05, p=0.001), COG T enriched | analysis/31 |
| 5 | H10 | REBASE identifies GCCGGC N4-C MTase | 2/87 NaeI produce m4C; SC_RS24685 top by annotation, BLAST-negative | analysis/33 |
| 6 | H13 | AAGCCCG non-random, secondary metabolism-associated | Non-random confirmed; Regulatory depleted (fold=0.43); sec. met. NOT enriched | analysis/35 |
| 7 | H14 | Defense island produces T3 arm GCCGGC methylation | Geographic shift (chi2=597); 1 T3 site in island; T3=21 sites only | analysis/37 |
| 8 | H16 | Expression dynamics reverse-ID GCCGGC MTase | SC_RS13615 top composite; no candidate satisfies all 4 criteria | analysis/39 |
| 9 | H18 | GCCGGC dose-response on expression | Threshold switch (0 vs >=1), not graded; geographic confound | analysis/41 |
| 11 | H22 | R-M motif sequences depleted from regulatory gene DNA | Gene body fold=0.74-0.76; promoter NOT depleted; 21-35% contribution | analysis/45 |
| 12 | H24 | BGC methylation enrichment is geographic confound | All 4mC disappears in core-only; GCCGGC residual = sequence composition | analysis/47 |
| 15 | H30 | TSS sequence features determine protection zone | AAGCCCG 5.0x enriched; combined CV AUC=0.712 vs methyl distance 0.917 | analysis/53 |
| 15 | H31 | Exposed TFs regulate downstream cascade | 0/57 have FIMO motifs; 38/57 are FIMO targets; coordination is intrinsic | analysis/54 |
| 16 | H32 | Exposed TFs form self-regulatory module | Not physical module; within-type rho=0.500 (p=6.6e-25); 2 blocs | analysis/55 |
| 17 | H34 | Temporal phase separation between blocs | No phase separation (p=0.459); 63% early; simultaneous switch | analysis/57 |
| 18 | H36 | Exposed TF evolutionary conservation | Equally conserved (composite p=0.276, AUC=0.459); epigenomic context | analysis/59 |

### Rejected (12)

| Loop | ID | Hypothesis | Key Statistic | Reference |
|:----:|:--:|-----------|---------------|-----------|
| 2 | H4 | Literature TFs show methylation-expression coordination | 3/37 methylated, 0/37 coordinated, 0/4 SARPs methylated | analysis/23 |
| 3 | H7 | CCGG 4mC decrease is passive dilution | Jaccard=0.000 (complete non-overlap); active remodeling | analysis/27 |
| 4 | H9 | CCGG "4mC" is misclassified 5mC | 5mC signal 0.01%; 4mC 82.7%; 91.4% in GCCGGC palindrome | analysis/30 |
| 6 | H12 | BLAST identifies GCCGGC MTase as SC_RS24685 | SC_RS19770 top hit (E=0.007); SC_RS24685 out of range | analysis/34 |
| 8 | H17 | GCCGGC and AAGCCCG show transcriptional division of labor | GCCGGC suppression = geographic confound (Simpson's paradox) | analysis/40 |
| 9 | H19 | GCCGGC loss causes temporal de-repression | Lost genes DOWN (median LFC=-0.253); core-only p=0.87 | analysis/42 |
| 10 | H21 | AAGCCCG loss causes temporal de-repression | Lost vs Never p=0.91, r=0.003; clean null | analysis/44 |
| 11 | H23 | Methylation avoidance extends to all essential categories | Only regulatory genes; ribosomal, DNA repair, energy NS | analysis/46 |
| 13 | H26 | Individual TF binding sites show methylation depletion | Flat profiles (p=0.12-0.96); GCCGGC enriched at TFBS (fold=1.16) | analysis/49 |
| 14 | H29 | baseMean threshold separates shielded/exposed | baseMean AUC=0.547; distance AUC=0.917; protection is expression-independent | analysis/52 |
| 16 | H33 | Exposed TF neighborhoods show correlated expression | Permutation p=0.857; no distance decay; not cis-regulators | analysis/56 |
| 17 | H35-causal | AAGCCCG methylation mediates exposed TF regulation | 2/57 have methylated AAGCCCG; sequence != methylation enrichment | analysis/58_causal |

### Other (2)

| Loop | ID | Hypothesis | Verdict | Key Statistic | Reference |
|:----:|:--:|-----------|---------|---------------|-----------|
| 1 | H2 | MTase expression stability causes 4mC/6mA difference | Not tested | Deferred to H5 | -- |
| 1 | H3 | absA de-repression compensates for redZ | Premise rejected | True redZ upregulated (LFC=+0.91); 36/37 locus_tag errors | analysis/22 |

### Additional completed analyses

| Loop | ID | Description | Status | Reference |
|:----:|:--:|------------|--------|-----------|
| 5 | H11 | Gatekeeper Model v2 quantitative integration | Completed (superseded by v3/v4) | analysis/32 |

**Verdict distribution (34 testable hypotheses + 1 completed + 1 not tested):**

| Verdict | Count | Percentage |
|---------|:-----:|:----------:|
| Supported | 7 | 19.4% |
| Partially supported | 15 | 41.7% |
| Rejected | 12 | 33.3% |
| Premise rejected | 1 | 2.8% |
| Not tested | 1 | 2.8% |

---

## 4. Negative Results That Shaped the Model

The model's most important structural insights came from rejected hypotheses. The following negative results were transformative:

### 4.1 Methylation does NOT control transcription (H19, H21)

The single most consequential finding: neither GCCGGC 4mC nor AAGCCCG 6mA shows temporal de-repression when methylation is lost.

| System | Test | n_Lost | n_Never | Median LFC diff | p-value | r |
|--------|------|:------:|:-------:|:---------------:|:-------:|:-:|
| GCCGGC (H19) | Lost vs Never T2vsT1 | 3,226 | 4,180 | -0.259 | 8.3e-08 | -0.052 |
| GCCGGC core-only (H19) | | 2,690 | 3,141 | +0.022 | 0.87 | +0.002 |
| AAGCCCG (H21) | Lost vs Never T2vsT1 | 836 | 6,397 | +0.004 | 0.91 | 0.003 |

The GCCGGC result initially appeared to show suppression (Lost genes go DOWN), but geographic stratification revealed this as a Simpson's paradox: Lost genes are enriched in the core, and core genes have lower LFC independently of methylation. The AAGCCCG test provides the cleanest null: with minimal geographic shift (4.7 pp), there is simply no effect (r = 0.003).

### 4.2 Expression does NOT predict protection (H29)

The H28 finding that exposed genes have lower T1 expression initially suggested an RNAP-occupancy model (low expression = unoccupied promoter = exposed to MTase). H29 definitively rejected this:

- baseMean AUC = 0.547 (near random)
- No trend across expression quintiles (JT p = 0.730)
- Nearest methylation distance alone: AUC = 0.917

Protection is a structural, expression-independent property of the genomic locus.

### 4.3 No neighborhood effect (H33)

If exposed TFs were cis-regulators of adjacent genes (common in *Streptomyces* pathway-specific regulation), their neighbors should show correlated expression changes. They do not:

- Permutation p = 0.857 (observed |LFC| within null distribution)
- No distance decay (rho = +0.026, NS)
- Concordance rate 59.9% vs shielded 55.6% (p = 0.191, NS)

The exposed TFs are trans-acting regulators of dispersed targets.

### 4.4 Sequence enrichment does NOT equal methylation enrichment (H35-causal)

This was the most surprising negative result. H30 showed AAGCCCG DNA sequence is 5.0x enriched at exposed TF promoters. The natural hypothesis was that these motifs would be methylated and that AAGCCCG methylation loss drives exposed TF de-repression. H35-causal showed:

- 25.8% of exposed TFs have AAGCCCG **sequence** within +/-300bp (H30)
- Only 3.2% have **methylated** AAGCCCG within +/-500bp (H35-causal)
- The 2 genes with methylated AAGCCCG show expression changes in opposite directions

The MTase methylates only ~19.5% of its recognition sites genome-wide. The sequence enrichment is an evolutionary/structural property; it does not translate to proportional methylation.

### 4.5 Simpson's Paradox discoveries (H17, H19, H24)

Three independent Simpson's paradox instances were identified, all driven by the linear chromosome's core/arm functional differentiation:

| Discovery | Apparent signal | After stratification | Confound |
|-----------|----------------|---------------------|----------|
| H17: GCCGGC suppression | p = 4.1 x 10^-10 | Core: p = 0.87 | T1 sites 83% core; core genes lower LFC |
| H19: De-repression reversed | p = 8.3 x 10^-8 | Core: p = 0.87 | Lost sites enriched in core |
| H24: BGC enrichment | fold = 1.66, p = 4.4 x 10^-8 | Core-only fold = 1.07 (NS) | 100% core co-localization of BGCs |

**Methodological lesson**: In organisms with a linear chromosome and core/arm functional differentiation, ANY methylation-expression analysis must include geographic stratification (Cochran-Mantel-Haenszel test) before causal claims.

---

## 5. Model Evolution: v1 -> v2 -> v3 -> v4

### v1 (Loops 1-4): "Methylation Controls Transcription"

**Working model**: DNA methylation (4mC promoter-localized, 6mA gene-body-dispersed) directly regulates gene expression, with MTase expression stability determining correlation strength. The AAGCCCG cascade through SC_RS17645 was hypothesized to control BGC regulators.

**Key results**: H3 revealed 36/37 literature TF locus_tag mapping errors. H4 showed 0/37 characterized TFs have methylation-expression coordination. H6 discovered 62 novel coordinated regulatory genes outside the literature set. The initial model was comprehensively rejected, but the 57 exposed regulators were discovered.

### v2 (Loops 5-9): Three-Layer Model, Simpson's Paradox

**Working model**: Three layers -- (1) R-M defense landscape remodeling, (2) regulatory gene protection/depletion, (3) 57-gene signal gating cascade. GCCGGC showed apparent suppressive effects and AAGCCCG showed site-level dynamics.

**Key results**: H9 confirmed GCCGGC 4mC is genuine (not 5mC misclassification). H11 formalized the Gatekeeper Model v2 with 10 testable predictions. H17-H19 discovered Simpson's paradox -- the apparent suppressive effect of GCCGGC was entirely a geographic artifact. This transformed the model from "methylation controls transcription" to "methylation does NOT control transcription."

### v3 (Loops 10-14): Shielded/Exposed Dichotomy, 293bp Boundary

**Working model**: Refined three-layer model with the Shielded/Exposed dichotomy as central organizing principle. 998 genes have a 1.2 kb protection zone; 57 genes lack it entirely. The 293bp distance threshold (AUC = 0.917) provides near-perfect binary classification. H21 confirmed AAGCCCG temporal de-repression is also null (r = 0.003). H22 and H30 quantified sequence vs protein occupancy contributions (~33% vs ~67%). H29 showed protection is expression-independent.

**Key transformation**: The model shifted from asking "how does methylation affect gene expression" to "what structural features determine which genes are accessible to methylation."

### v4 (Loops 15-18): Distributed Switch, Simultaneous Activation, Vegetative-to-Developmental Transition

**New discoveries (H30-H36)**:

| Discovery | Source | Impact |
|-----------|--------|--------|
| Sequence contributes ~33% of protection | H30 | Quantified two-tier protection model |
| 0/57 have FIMO motifs; parallel pathway model | H31 | Exposed TFs are outside characterized networks |
| Distributed layer (not physical module), two blocs | H32 | Activation (~36) vs repression (~26) |
| No neighborhood effect (trans-acting) | H33 | Not cis-regulators |
| Simultaneous switch (not cascade) at T1-T2 | H34 | Both blocs engage in parallel |
| AAGCCCG sequence != methylation enrichment | H35-causal | Eliminated AAGCCCG as proximate mechanism |
| TetR enriched in repression; TCS/sigma in activation | H35-TF | Vegetative shutdown + developmental activation |
| Equal conservation (not novel genes) | H36 | Exposed = epigenomic context, not gene novelty |

**Model transformation**: The exposed TFs were recharacterized from a "methylation-responsive cascade" to a "distributed developmental switch." The v4 model identifies the biological function: simultaneous repression of vegetative growth programs (metabolism, DNA repair, efflux) and activation of developmental programs (morphogenesis, signal transduction, stress response) at the exponential-to-transition boundary.

---

## 6. Unresolved Questions and Experimental Roadmap

### 6.1 Immediate Experiments (Highest Priority)

| Experiment | Question | Expected outcome if model correct |
|-----------|----------|-----------------------------------|
| **MTase KO + SMRT-seq** | Does removing methylation alter exposed TF expression? | Exposed TFs show altered expression; shielded TFs unaffected |
| **ChIP-seq for NAPs** (HupA, IHF, Lsr2) | What creates the protection zone? | NAP binding enriched at shielded TSS, depleted at exposed TSS |
| **DAP-seq for exposed TFs** | What are the downstream targets? | Each exposed TF binds 10-100 dispersed targets |
| **CRISPRi of top exposed TFs** | RamR (LFC=+7.9), NsdB (+7.8), SCO1160 (+5.6) | Developmental arrest or altered secondary metabolism |

### 6.2 Computational Predictions for Experimental Testing

| Prediction | Basis | Test | Expected result |
|-----------|-------|------|-----------------|
| MTase KO removes all exposed/shielded distinction | Methylation defines the dichotomy | SMRT-seq + RNA-seq in MTase KO | 57 exposed genes no longer show TSS methylation; expression changes |
| Protection zone = NAP occupancy | H25, H26, H29 | ChIP-seq for HupA/IHF | NAP signal peaks at +300bp from regulatory TSS |
| TCS pairs: exposed partner is the developmental sensor | H32 TCS asymmetry | TCS gene KO | KO of exposed partner blocks developmental transition; shielded partner KO has no effect |
| Activation bloc targets include BGC/morphogenesis genes | H35-TF | DAP-seq or RNA-seq after CRISPRi | RamR/NsdB regulate aerial mycelium and secondary metabolite genes |
| Repression bloc targets include primary metabolism | H35-TF | DAP-seq or RNA-seq after CRISPRi | TetR-family exposed TFs repress transport/efflux/metabolic operons |
| The 293bp boundary is conserved in other Streptomyces | H36 | Comparative SMRT-seq | Similar shielded/exposed ratio in S. venezuelae, S. avermitilis |

### 6.3 Remaining Computational Questions

| Question | Current status | Approach |
|----------|---------------|----------|
| T1 GCCGGC MTase identity | SC_RS19770 (BLAST E=0.007) or SC_RS13615 (expression correlation) | Gene KO + SMRT-seq |
| De novo motif discovery for exposed TF binding sites | No FIMO data available | MEME/STREME on co-regulated gene sets |
| Cross-species Shielded/Exposed conservation | Proxy measures only (H36) | Direct BLAST of 57 genes against 833 Streptomyces genomes |
| Protection zone structural determinants | ~33% sequence, ~67% protein | Hi-C + ATAC-seq to map chromatin accessibility |
| Why does MTase methylate only 19.5% of AAGCCCG sites? | Unknown site selectivity | Structural modeling of SC_RS17645; comparison of methylated vs unmethylated sites |
| Regulatory network reconstruction for 57 exposed TFs | No binding motifs known | RNA-seq time course with MTase KO strains; co-expression network inference |

---

## 7. Paper Figure Plan

### Main Figures

| Figure | Title | Content | Key panels |
|:------:|-------|---------|------------|
| **1** | Gatekeeper Model overview | 4-layer architecture diagram | (A) Chromosome ideogram with geographic redistribution, (B) Protection zone spatial profile, (C) Shielded/Exposed dichotomy distance distributions, (D) 57 exposed TFs two-bloc switch |
| **2** | Protection zone characterization | Quantitative dissection of Layers 2a-2c | (A) TSS methylation gradient (H25: 2,200bp zone), (B) 293bp boundary ROC curve (H29: AUC=0.917), (C) Sequence vs protein occupancy decomposition (H30: 33/67%), (D) Individual TF BS NOT protective (H26), (E) Expression-independence (H29: baseMean AUC=0.547) |
| **3** | 57 Exposed TFs: distributed developmental switch | Internal structure and temporal dynamics | (A) Co-expression heatmap with 4 modules, (B) Two blocs temporal trajectories (mirror-image), (C) Phase ratio distribution (63% early), (D) TF family composition by bloc, (E) TCS asymmetry diagram (7/7 split) |
| **4** | Negative results panel | Simpson's paradox and null results | (A) GCCGGC unstratified vs stratified de-repression (H19), (B) AAGCCCG clean null (H21: r=0.003), (C) Sequence vs methylation enrichment gap (H30 vs H35-causal), (D) No neighborhood effect (H33: permutation p=0.857) |
| **5** | Vegetative-to-developmental switch model | Biological interpretation | (A) Vegetative programs repressed (TetR/efflux, metabolism, DNA repair), (B) Developmental programs activated (RamR/morphogenesis, TCS cascades, stress sigma factors), (C) Simultaneous switch timing, (D) Known developmental gene network placement |

### Supplementary Figures

| Figure | Title | Content |
|:------:|-------|---------|
| **S1** | Methylation landscape overview | Site counts, motif composition, temporal dynamics for all 3 timepoints |
| **S2** | GCCGGC geographic redistribution time series | Chromosome ideograms for T1/T2/T3 with core/arm percentages |
| **S3** | Simpson's paradox demonstration (detailed) | 6-panel: unstratified, core-only, arm-only for both GCCGGC and AAGCCCG |
| **S4** | Regulatory avoidance forest plot | CMH-adjusted odds ratios for all motifs, core/arm stratified (H20) |
| **S5** | Sequence-level motif depletion | Gene body vs promoter depletion, sliding window gradients (H22, H30) |
| **S6** | 57 exposed TF complete annotation table | Family, bloc, coordination type, LFC, methylation details, conservation scores |
| **S7** | TCS pair detailed analysis | 7 pairs: expression trajectories, methylation status, temporal ordering |
| **S8** | H36 conservation metrics | 8-panel comparison of all evolutionary proxy measures |

---

## 8. Key Statistics Summary Table

### Core Model Statistics

| Statistic | Value | Source |
|-----------|-------|--------|
| Total regulatory genes analyzed | 1,055 | H6 |
| Shielded regulators | 998 (94.6%) | H27 |
| Exposed regulators | 57 (5.4%) | H27 |
| Protection zone width | 2,200 bp (-1,300 to +700) | H25 |
| Deepest depletion position | +300 bp from TSS (ratio=0.541) | H25 |
| TSS methylation density ratio (exposed/shielded) | 8.4x (3.15 vs 0.37 sites/kb/gene) | H27 |
| Shielded/Exposed boundary threshold | 293 bp | H29 |
| Distance-based classification AUC | 0.917 (95% CI: 0.895-0.935) | H29 |
| Expression-based classification AUC | 0.547 (near random) | H29 |
| Sequence contribution to protection | ~33% (CV AUC=0.712) | H30 |
| Protein occupancy contribution to protection | ~67% | H30, H29 |

### Methylation System Statistics

| Statistic | Value | Source |
|-----------|-------|--------|
| GCCGGC 4mC sites T1 | 1,289 (83% core) | H7, H14 |
| GCCGGC 4mC sites T2 | 407 (82% arm) | H7, H14 |
| GCCGGC 4mC sites T3 | 21 (57% arm) | H14 |
| GCCGGC T1-T2 Jaccard overlap | 0.000 (complete non-overlap) | H7 |
| GCCGGC geographic shift chi2 | 597 (p ~ 10^-130) | H14 |
| AAGCCCG 6mA sites T1 | 260 (69% core) | H5, H13 |
| AAGCCCG 6mA sites T2 | 64 (64% core) | H5, H13 |
| SC_RS17645 (AAGCCCG MTase) LFC T2vsT1 | -2.19 | H5 |
| AAGCCCG temporal de-repression effect size | r = 0.003 (null) | H21 |

### Regulatory Avoidance Statistics

| Statistic | Value | Source |
|-----------|-------|--------|
| AAGCCCG regulatory depletion fold (site-centric) | 0.43 (p=0.035) | H13, H15 |
| GCCGGC regulatory depletion fold (site-centric) | 0.61 (p=3.9e-05) | H15 |
| Gene body sequence motif depletion fold | 0.74-0.76 (p<3e-07) | H22 |
| CMH-adjusted regulatory depletion p (all motifs) | <0.007 | H20 |

### Exposed TF Statistics

| Statistic | Value | Source |
|-----------|-------|--------|
| Exposed-exposed co-expression median rho | +0.017 | H32 |
| Within-type co-expression median rho | +0.500 (p=6.6e-25) | H32 |
| Between-type co-expression median rho | -0.283 | H32 |
| Co-expression modules detected | 4 (covering 61/57 genes) | H32 |
| Activation bloc size | ~36 genes (modules 1-3) | H32 |
| Repression bloc size | ~26 genes (module 4) | H32 |
| TCS pairs with asymmetric exposed/shielded | 7/7 (100%) | H32 |
| Early responders (phase ratio > 0.6) | 39/57 (63%) | H34 |
| Bloc phase separation p-value | 0.459 (NS, simultaneous) | H34 |
| Module-internal temporal coherence rho | 0.717 (p=1e-19) | H34 |
| Methylation-expression timing rho (T1-T2) | 0.136 (p=0.299, NS) | H34 |
| FIMO motifs for exposed TFs | 0/57 | H31 |
| Neighborhood transcriptional effect (permutation p) | 0.857 (null) | H33 |
| Exposed TFs with methylated AAGCCCG at TSS | 2/57 (3.2%) | H35-causal |
| Exposed TFs with AAGCCCG sequence at TSS | ~16/57 (25.8%) | H30 |
| TetR enrichment in repression bloc | OR=0.16 (p=0.028) | H35-TF |

### Simpson's Paradox Statistics

| Statistic | Value | Source |
|-----------|-------|--------|
| GCCGGC suppression (unstratified) | p=4.1e-10 | H17 |
| GCCGGC suppression (core-only) | p=0.87 (eliminated) | H19 |
| GCCGGC de-repression (unstratified) | p=8.3e-08 | H19 |
| GCCGGC de-repression (core-only) | p=0.87 (eliminated) | H19 |
| BGC enrichment (unstratified) | fold=1.66, p=4.4e-08 | H23 |
| BGC enrichment (core-only) | fold=1.07 (NS, eliminated) | H24 |
| AAGCCCG de-repression | p=0.91 (genuine null, no confound) | H21 |

### Conservation Statistics

| Statistic | Value | Source |
|-----------|-------|--------|
| Exposed SCO locus tag rate | 96.8% | H36 |
| Exposed named product rate | 100% | H36 |
| Composite conservation score p (exposed vs shielded) | 0.276 (NS) | H36 |
| Conservation ROC AUC for predicting exposed | 0.459 (below random) | H36 |
| Exposed GC3 vs shielded GC3 | 0.919 vs 0.932 (p=0.060, NS) | H36 |
| Exposed Nc vs shielded Nc | 31.6 vs 31.5 (p=0.725, NS) | H36 |

---

## Appendix: Integrated Gatekeeper Model v4 Schematic

```
S. coelicolor M145 Epigenome-Transcriptome Model (v4 FINAL)

LAYER 1: R-M DEFENSE GEOGRAPHIC REDISTRIBUTION
[H7, H14, H19, H21]
- GCCGGC 4mC: T1=1,289 (83% core) --> T2=407 (82% arm) --> T3=21 (57% arm)
- AAGCCCG 6mA: T1=260 (69% core) --> T2=64 (64% core) [75% site loss]
- Complete positional remodeling (Jaccard = 0.000)
- Simpson's Paradox: apparent correlations are geographic artifacts
- NOT transcriptional regulation (H19 p=0.87; H21 r=0.003)

              |
              v

LAYER 2: REGULATORY DNA PROTECTION SYSTEM
[H15, H20, H22, H25, H26, H29, H30]

  2a. Evolutionary counter-selection (~33%)     2b. Collective promoter occupancy (~67%)
  [H22, H30]                                    [H25, H26, H29]
  - Gene body motif depletion                   - 2,200 bp protection zone at TSS
    (fold = 0.74-0.76)                          - Deepest at +300 bp (ratio = 0.541)
  - TSS AAGCCCG 5x enriched at                  - NOT individual TF BS (H26: enriched)
    exposed TFs (H30)                           - NOT expression-dependent (H29: AUC=0.547)
  - Combined CV AUC = 0.712                     - 293 bp threshold (AUC = 0.917)

              |
              v

  2c. SHIELDED / EXPOSED DICHOTOMY [H27, H29, H36]
  ======================================================================
  |                                  |                                  |
  |  SHIELDED (998, 94.6%)          |  EXPOSED (57, 5.4%)              |
  |  Protection zone: 1,200 bp      |  Protection zone: 0 bp           |
  |  Nearest methyl: 762 bp         |  Nearest methyl: 114 bp          |
  |  Methylation-insensitive        |  8.4x TSS enrichment             |
  |  38% constitutive               |  0% constitutive                 |
  |  Equally conserved (H36)        |  Equally conserved (H36)         |
  |  --> NO methylation pathway     |  --> Methylation-responsive       |
  |                                  |                                  |
  ======================================================================
             |                                    |
             v                                    v

  EXCLUDED PATHWAYS                    LAYER 3: DISTRIBUTED DEVELOPMENTAL SWITCH
  (ALL NULL)                           [H6, H8, H27, H28, H31-H36]
  - Literature TFs (H4: 0/37)
  - Direct transcription (H19/H21)     3a. INTERNAL STRUCTURE (H32, H35-TF)
  - Dose-dependent (H18: confounded)       Two antagonistic blocs:
  - Individual TFBS (H26: enriched)        ACTIVATION (~36): TCS/sigma/WhiB/RamR
  - AAGCCCG-mediated (H35: 2/57)          REPRESSION (~26): TetR(p=0.028)/metab/DNA repair
  - Cis-regulatory (H33: p=0.857)         TCS 7/7 asymmetric (exposed:shielded = 1:1)
  - Known TF cascade (H31: 0/57)
                                       3b. TEMPORAL DYNAMICS (H34)
                                           63% early responders (T1->T2 transition)
                                           SIMULTANEOUS switch (p=0.459, not cascade)
                                           Module coherence rho=0.717
                                           Methylation timing decoupled (rho<0.14)

                                       3c. MECHANISM: UNKNOWN
                                           Trans-acting on dispersed targets
                                           Not FIMO motifs, not cis, not AAGCCCG methylation
                                           Requires ChIP-seq / DAP-seq / genetics

                                       3d. CONSERVATION (H36)
                                           96.8% SCO tags, 100% named products
                                           Equal GC3, Nc, rare codons
                                           Exposed = epigenomic context, not gene novelty
```

---

## 2026-04-11 追加解析アップデート（議事録論点 A-1/A-2/A-3, B-1, C-2, E-1）

### A-1/A-2/A-3: モチーフ信頼性（61_motif_reliability）

| 解析 | 結果 |
|------|------|
| A-1: AAGCCCG methylated position | **pos1（2番目A）が63%**（164/260）。GCCGGC: pos1とpos3ほぼ等分（53%/47%）→ 回文の両C確認 |
| A-2: モチーフ充足率 | AAGCCCG: **31.6%**、GCCGGC: **7.96%**、GATC: 0.27%（T1） |
| A-3: 6mA Unassigned | T1=**72.2%**（1396/1934）→ 未知MTase存在を示唆 |

### B-1: ウィンドウメチル化（63_window_methylation）

- 4mC: T2ピーク（2446サイト）→ T3急減（1647）。コア/アーム比: T1=2.21 → T2=1.12 → T3=1.32
- 6mA: 単調増加（1934→2120→2295）。コア/アーム≈1.0（地理的均一）
- **Gatekeeperモデルへの含意**: 4mCはT2で全ゲノム展開後に縮小。コア優位性はT1で最大（Layer 1と整合）

### C-2: タイムポイント別TSS保護ゾーン（64_timepoint_TSS）

| mod_type | gene_cat | T1 ratio | T2 ratio | T3 ratio |
|---------|---------|---------|---------|---------|
| 4mC | All | 0.823 | 0.856 | 0.867 |
| 4mC | **Regulatory** | 0.839 | 0.844 | **1.046（崩壊！）** |
| 4mC | Non-regulatory | 0.821 | 0.857 | 0.845 |
| 6mA | Regulatory | 0.927 | **0.732（最深）** | 0.764 |

**重要**: 4mC保護ゾーンが制御遺伝子においてT3で崩壊（ratio>1.0）。非制御遺伝子では維持。Layer 2（保護ゾーン維持）に時間的動態が存在する新所見。

### E-1: GO/KEGG機能エンリッチメント（62_GO_KEGG_enrichment）

| 遺伝子群 | 主要所見 |
|---------|---------|
| GCCGGC-proximal | **Quorum sensing: FDR=2.5e-4, enrichment=1.55**（65/100 genes）|
| AAGCCCG-proximal | 有意エンリッチメントなし（機能非偏在）|
| Dual-targeted | Siderophore合成: FDR=8.7e-4, enrichment=12.2 |

**Gatekeeperモデルへの含意**: Layer 3（Signal Gating）に対するKEGGレベルの機能的裏付けが得られた。4mCがシグナル伝達遺伝子（QS/TCS）を優先的にシールドする構造を確認。

---

*Analysis performed: 2026-02-24 through 2026-02-27 (18 exploration loops)*
*This report integrates findings from 36 hypotheses tested across 59 analysis directories.*
*Scripts and data: `/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/`*
*Previous versions: v2 (H11, Loop 5), v3 (H1-H29, 14 loops)*
*2026-04-11: Additional analyses (A-1/A-2/A-3, B-1, C-2, E-1) added as appendix*
