# Paper Outline: The Methylation-Exclusion Architecture of Bacterial DNA Regulation

**Working title**: A restriction-modification methylation-exclusion mechanism partitions the regulatory genome into a methylation-insensitive majority and a developmental-switch minority in *Streptomyces coelicolor*

**Target journal**: TBD (format assumes broad-interest journal: Nature Microbiology / PNAS / eLife)

**Writing order**: Results → Methods → Discussion → Introduction → Abstract

---

## Abstract (~250 words)

**Background** (2 sentences): Bacterial DNA methylation is widely assumed to regulate transcription, yet genome-wide evidence for this function remains inconsistent. *Streptomyces coelicolor* M145, with two active methylation systems (GCCGGC 4mC and AAGCCCG 6mA), provides an ideal model to test this assumption using integrated epigenomic and transcriptomic data across developmental stages.

**Results** (4-5 sentences): We show that apparent genome-wide methylation-expression correlations are spurious associations driven by geographic confounding from the linear chromosome's core/arm functional differentiation. Instead, a four-layer methylation-exclusion architecture governs methylation-transcription interactions: (1) R-M systems undergo complete geographic redistribution (Jaccard = 0.000) without transcriptional consequences; (2) 94.6% of regulatory gene promoters (998/1,055) maintain a ~2.2 kb methylation-free protection zone through collective promoter occupancy (~67%) and evolutionary counter-selection (~33%); (3) 57 exposed regulators (5.4%) lack this protection entirely, defined by a 293 bp boundary (AUC = 0.917); (4) these 57 genes constitute a distributed vegetative-to-developmental switch organized into two synchronous antagonistic blocs.

**Conclusion** (2 sentences): DNA methylation does not control transcription directly in *S. coelicolor*. Instead, a structural dichotomy at regulatory gene promoters channels methylation effects through a specific minority of regulators that encode a coordinated developmental transition, revealing an unexpected organizational principle of bacterial regulatory genomes.

---

## Introduction (~1,500 words, 6 paragraphs)

### Para 1: DNA methylation in bacteria — established functions vs open questions
- R-M systems: defense against foreign DNA (well established)
- Orphan MTases and phase variation: gene regulation in pathogens (Dam, CcrM)
- Gap: genome-wide regulatory roles remain poorly tested outside a few model systems
- **Cite**: Blow et al. 2016 (PLoS Genet), Casadesus & Low 2006 (Microbiol Mol Biol Rev), Sanchez-Romero & Casadesus 2020

### Para 2: *Streptomyces* as a model — linear chromosome and complex development
- Linear chromosome with core/arm architecture (8.67 Mb)
- Complex developmental cycle: vegetative growth → aerial mycelium → sporulation
- Secondary metabolism coupled to development
- Known developmental regulators: BldA/BldH cascade, WhiB/WhiG, sigma factors
- No previous integrated epigenome-transcriptome analysis
- **Cite**: Bentley et al. 2002 (Nature), Hopwood 2007, Flardh & Buttner 2009

### Para 3: M145 methylome — two systems, dramatic temporal dynamics
- GCCGGC (N4-methylcytosine, 4mC): 1,289 → 407 → 21 sites across T1/T2/T3
- AAGCCCG (N6-methyladenine, 6mA): 260 → 64 sites
- PacBio SMRT-seq + RNA-seq across 3 developmental timepoints
- Initial expectation: temporal methylation changes regulate developmental gene expression

### Para 4: The challenge of confounded inference in linear chromosomes
- Linear chromosome creates geographic functional bias (core: essential/housekeeping; arms: adaptive/secondary)
- Methylation systems have geographic biases (T1 GCCGGC: 83% core)
- Any correlation between methylation proximity and expression change may be geographic confound
- Simpson's paradox: aggregate correlations can reverse when stratified
- Importance of geographic stratification (CMH test) for causal inference

### Para 5: From descriptive methylomics to mechanistic model — iterative hypothesis testing
- 18 exploration loops, 36 hypotheses: systematic construction of the model
- 12 rejected hypotheses were the most informative (negative results as structural insights)
- The model emerged from the failure of the "methylation controls transcription" paradigm
- Brief overview of the four-layer methylation-exclusion architecture

### Para 6: Scope and aims of this study
- Aim 1: Determine whether DNA methylation directly regulates transcription genome-wide
- Aim 2: Characterize the structural basis of methylation depletion at regulatory gene promoters
- Aim 3: Identify and characterize the methylation-responsive regulatory minority
- Aim 4: Determine the biological function of methylation-responsive regulators in the developmental transition

---

## Results (~4,500 words, organized by Figures 1-6)

### Section 1: R-M systems undergo geographic redistribution without transcriptional consequences (Fig. 2a)

**Para 1**: GCCGGC geographic redistribution
- T1: 1,289 sites (83% core) → T2: 407 sites (82% arm) → T3: 21 sites (62% arm)
- Chi-squared = 597, p ~ 10^-130 for geographic shift
- Jaccard = 0.000: complete positional non-overlap between T1 and T2 → active remodeling, not passive dilution
- Consistent with R-M defense function (phage pressure shifts to arms during development)
- **Data**: Fig. 2a, Table S1

**Para 2**: AAGCCCG 6mA loss
- T1: 260 sites (69% core) → T2: 64 sites (64% core): 75% loss with minimal geographic shift (4.7 pp)
- Correlates with SC_RS17645 (HsdM-type MTase) downregulation (LFC = -2.19)
- **Data**: Fig. 1, Table S2

**Para 3**: Neither system causes transcriptional de-repression
- GCCGGC: Lost vs Never genes LFC comparison → unstratified p = 8.3 × 10^-8 (apparent effect)
- Geographic stratification: core-only p = 0.87 → Simpson's paradox (Fig. 5a)
- AAGCCCG: r = 0.003, p = 0.91 — clean null with no geographic confound (Fig. 5b)
- Conclusion: methylation changes are defense system dynamics, not transcriptional regulation

### Section 2: A 2,200 bp protection zone shields regulatory gene promoters (Fig. 2b, Fig. 3)

**Para 4**: Spatial profile of methylation depletion
- 1,055 regulatory genes vs all genes: methylation density profile ±5 kb from TSS
- Regulatory genes show 2,200 bp zone of depletion (-1,300 to +700 from TSS)
- Deepest depletion at +300 bp downstream of TSS (ratio = 0.541 of flanking density)
- Non-regulatory genes: no depletion
- **Data**: Fig. 2b, Fig. 3a

**Para 5**: Protection is expression-independent
- baseMean AUC = 0.547 (near random) for predicting exposed/shielded status
- No trend across expression quintiles (JT p = 0.730)
- Constitutive expression: 38% of shielded, 0% of exposed — but this is a consequence, not a cause
- **Data**: Fig. 3b, Fig. 3e

**Para 6**: Two-tier protection mechanism
- Sequence contribution: AAGCCCG motif 4.10× enriched (mean count) at exposed TF promoters; combined sequence model CV AUC = 0.712 → ~33% contribution
- Protein occupancy contribution: nearest methylation distance AUC = 0.917 >> sequence AUC → ~67% contribution
- Individual TF binding sites are NOT protective: GCCGGC fold = 1.157 at TFBS (counter-intuitive enrichment)
- Collective promoter occupancy (NAPs, RNAP, general factors) rather than specific TF binding
- **Data**: Fig. 3b, Fig. 3c, Fig. 3d

### Section 3: A sharp 293 bp boundary defines the Shielded/Exposed dichotomy (Fig. 2c)

**Para 7**: Binary classification of regulatory genes
- 998 shielded (94.6%): median nearest methylation = 762 bp
- 57 exposed (5.4%): median nearest methylation = 114 bp
- Mann-Whitney p = 1.7 × 10^-25; 6.7× distance ratio
- 293 bp threshold: AUC = 0.917 (95% CI: 0.895-0.935), 100% sensitivity, 80.6% specificity
- This is a near-binary structural feature, not a graded continuum
- **Data**: Fig. 2c, Fig. 3b

**Para 8**: Exposed regulators are equally conserved
- SCO locus tag rate: 96.8% (exposed) vs 97.8% (shielded)
- GC3, effective number of codons, rare codon frequency: all NS
- Composite conservation AUC = 0.459 (below random)
- Exposed status is epigenomic context, not evolutionary novelty
- **Data**: Fig. S7

### Section 4: 57 exposed regulators form a distributed developmental switch (Fig. 4)

**Para 9**: Internal structure — two antagonistic blocs
- Co-expression analysis reveals 4 modules covering 61/57 genes
- Within-type co-expression rho = 0.500 (p = 6.6 × 10^-25)
- Between-type co-expression rho = -0.283 → mirror-image trajectories
- NOT a physical genomic cluster: z = +0.30, p = 0.37; 0 operonic pairs
- **Data**: Fig. 4a

**Para 10**: Activation bloc (34 genes)
- TCS (9 genes, 26%), sigma factors (5), WhiB (1), RamR
- Key genes: RamR (SCO6685, LFC = +7.9), NsdB (LFC = +7.8), SCO1160 (sensor kinase, LFC = +5.6)
- Function: morphogenesis, signal transduction, stress response → developmental programs ON
- **Data**: Fig. 4d, Fig. 6b

**Para 11**: Repression bloc (23 genes)
- TetR family significantly enriched (OR = 0.16, p = 0.028)
- Also GntR, IclR, LacI, SSB, HU
- Function: metabolism shutdown, DNA repair cessation, efflux repression → vegetative programs OFF
- **Data**: Fig. 4d, Fig. 6a

**Para 12**: TCS asymmetry
- 7/7 cognate TCS pairs: exactly one exposed, one shielded
- Binomial p < 10^-3 for this perfect asymmetry
- Distribution: 4 SK : 3 RR exposed (no bias, p = 1.0)
- Suggests methylation selectively exposes the developmental-sensing component
- **Data**: Fig. 4e

### Section 5: Simultaneous activation at the vegetative-to-developmental transition (Fig. 6)

**Para 13**: Temporal dynamics
- 63% early responders (phase ratio > 0.6): bulk of change at T1→T2 transition
- Both blocs engage simultaneously (phase separation p = 0.459)
- Module-internal coherence: rho = 0.717 (p = 1.0 × 10^-19)
- Modules respond in parallel, not in cascade
- **Data**: Fig. 4b, Fig. 4c, Fig. 6c

**Para 14**: Methylation timing is decoupled from expression timing
- T1-T2: rho = 0.136, p = 0.299 (NS)
- T2-T3: rho = 0.012, p = 0.927 (NS)
- Methylation changes and expression changes are temporally independent
- The "gate" opens once at the developmental boundary; methylation provides the structural context, not the temporal signal

### Section 6: Negative results define the model's boundaries (Fig. 5)

**Para 15**: Simpson's paradox — three instances
- GCCGGC suppression: p = 4.1 × 10^-10 unstratified → p = 0.87 core-only
- GCCGGC de-repression: p = 8.3 × 10^-8 unstratified → p = 0.87 core-only
- BGC enrichment: fold = 1.66, p = 4.4 × 10^-8 → fold = 1.07 (NS) core-only
- AAGCCCG provides internal control: p = 0.91 with no geographic confound
- **Data**: Fig. 5a, Fig. 5b, Fig. S2

**Para 16**: Sequence ≠ Methylation
- AAGCCCG DNA sequence 4.10× enriched (mean motif count) at exposed TF promoters (24.6%, 14/57 have motifs vs 6.4%, 64/994 shielded)
- Only 3.2% have methylated AAGCCCG (MTase methylates ~19.5% of recognition sites)
- Sequence enrichment is evolutionary/structural; methylation enrichment requires active MTase targeting
- **Data**: Fig. 5c

**Para 17**: No neighborhood effect
- Permutation p = 0.857; no distance decay (rho = +0.026)
- Concordance rate: 59.9% (exposed) vs 55.6% (shielded), p = 0.191
- Exposed TFs are trans-acting regulators of dispersed targets, not cis-regulators
- **Data**: Fig. 5d

**Para 18**: Mechanism remains unknown
- 0/57 have curated FIMO binding motifs
- 38/62 are FIMO targets of known TFs (same rate as shielded: OR = 0.881, p = 0.683)
- Not AAGCCCG methylation-mediated (2/57)
- Not cis-regulatory, not cascade, not known TF network
- Experimental approaches required: ChIP-seq, DAP-seq, CRISPRi

---

## Discussion (~2,000 words, 7 paragraphs)

### Para 1: Summary of central findings
- DNA methylation does NOT control transcription directly in *S. coelicolor*
- Instead, a methylation-exclusion architecture channels methylation effects through 57 specific regulatory genes
- The Shielded/Exposed dichotomy (AUC = 0.917) is the central organizing principle
- The 57 exposed regulators encode a coordinated vegetative-to-developmental switch

### Para 2: Reinterpretation of bacterial methylation-transcription paradigm
- Dam/CcrM paradigm: direct regulatory roles at specific loci (phase variation, replication timing)
- Our findings suggest an alternative paradigm for R-M associated methylation: structural partitioning of the regulatory genome
- Not all methylation is regulatory; most methylation in M145 is defense-related
- The 5.4% (57/1,055) exposed regulators represent the true methylation-responsive fraction
- Compare with eukaryotic methylation: CpG islands near promoters are also depleted (analogy with protection zone)

### Para 3: Simpson's paradox as a methodological warning
- Three independent instances of Simpson's paradox in a single organism
- Linear chromosomes with core/arm differentiation are common in Actinobacteria
- Any epigenome-transcriptome correlation study in these organisms MUST include geographic stratification
- Propose CMH test as standard for organisms with non-random chromosome architecture
- Broader relevance: analogous confounds may exist in other bacteria with replication-coupled methylation biases

### Para 4: The protection zone — comparison with known mechanisms
- 2,200 bp zone centered on TSS resembles nucleosome-free regions in eukaryotes
- Collective protein occupancy (~67%) as the dominant mechanism
- Individual TF binding is NOT protective (counter-intuitive GCCGGC enrichment at TFBS)
- NAP proteins (HupA, IHF, Lsr2) as likely candidates for the protection zone
- Predict: ChIP-seq for NAPs will show enrichment at shielded TSS, depletion at exposed TSS

### Para 5: The 57 exposed regulators — a novel regulatory layer
- Not part of any characterized pathway (0/57 FIMO motifs)
- Distinct from BldA/BldH cascade, WhiB/WhiG pathway, known sigma factor cascades
- TCS asymmetry (7/7 split) suggests methylation selectively gates one component of sensory systems
- Parallel, not sequential, activation: the switch operates as a single gate-opening event
- TetR enrichment in repression bloc (p = 0.028): consistent with role in shutting down vegetative functions

### Para 6: Evolutionary implications
- Exposed regulators are equally conserved (H36): not evolutionary newcomers
- Exposed status determined by epigenomic context (local methylation landscape), not by gene identity
- Suggests the methylation-exclusion architecture is a conserved organizational principle in Streptomyces
- Predicts similar Shielded/Exposed dichotomy in S. venezuelae, S. avermitilis, S. griseus
- The ~33% sequence contribution implies co-evolution of motif density and regulatory function

### Para 7: Limitations and future directions
- **Computational limitations**: correlation-based analysis cannot establish causation
- **Key experiments needed**: (1) MTase KO + SMRT-seq to test if removing methylation alters exposed TF expression, (2) ChIP-seq for NAPs to identify protection zone determinants, (3) DAP-seq for exposed TFs to map downstream targets, (4) CRISPRi of RamR/NsdB to test developmental necessity
- **Mechanism unknown**: the proximate mechanism linking methylation exposure to coordinated expression remains unidentified
- **Single strain**: analysis limited to M145; comparative SMRT-seq across Streptomyces species needed
- **Temporal resolution**: 3 timepoints may miss intermediate dynamics

---

## Methods (~2,500 words)

### Bacterial strain and growth conditions
- *S. coelicolor* A3(2) strain M145
- Growth conditions, media, timepoints (T1: exponential, T2: transition, T3: stationary)
- Biological replicates

### DNA methylation detection (SMRT-seq)
- PacBio SMRT sequencing
- Methylation calling pipeline (kinetic modification detection)
- High-confidence site filtering criteria
- Motif identification (GCCGGC, AAGCCCG)
- Timepoint-specific site calling

### RNA-seq and differential expression
- RNA extraction, library preparation
- Alignment and quantification pipeline
- DESeq2 analysis: T2 vs T1, T3 vs T1, T3 vs T2
- Definition of DEGs (adjusted p < 0.05, |LFC| > 1)
- baseMean calculation

### Regulatory gene classification
- Definition of regulatory genes (1,055 genes): TFs, TCS components, sigma factors
- Source databases and annotation criteria
- TSS assignment methodology

### Methylation-TSS distance analysis
- Nearest methylation site distance calculation for each regulatory gene
- Protection zone spatial profile: TSS-centered ±5 kb density calculation with 100 bp bins
- Comparison: regulatory vs non-regulatory genes
- Bootstrap confidence intervals

### Shielded/Exposed classification
- ROC analysis for nearest_methyl_distance threshold
- 293 bp optimal boundary determination (Youden's J statistic)
- AUC calculation with 95% CI (DeLong method or bootstrap)
- Mann-Whitney U test for distance distributions
- baseMean as negative control predictor

### Sequence determinants of protection
- FIMO motif scanning: GCCGGC (TGGCCGGC palindrome) and AAGCCCG within ±300 bp of TSS
- Logistic regression model: sequence features predicting exposed/shielded status
- 5-fold cross-validation for sequence model AUC
- Decomposition: sequence contribution (~33%) vs residual protein occupancy (~67%)

### TF binding site methylation analysis
- FIMO: known TF binding motifs from CollecTF and CIS-BP
- TFBS-centered methylation profile (±2 kb)
- Enrichment fold calculation (observed/random ratio at TFBS)

### Geographic stratification and Simpson's paradox
- Core region: 1.5-7.17 Mb; Arm regions: 0-1.5 Mb (left), 7.17-8.67 Mb (right)
- Cochran-Mantel-Haenszel (CMH) test for geographic stratification
- Gene methylation transitions: "Lost" (methylated at T1, unmethylated at T2) vs "Never" (never methylated)
- LFC comparison: Lost vs Never, both unstratified and stratified by core/arm

### Co-expression analysis and module detection
- 57 × 57 Spearman co-expression matrix across all timepoints/replicates
- Hierarchical clustering with Ward's method
- Module detection (4 modules, 61/57 coverage)
- Within-type vs between-type co-expression rho comparison
- Bloc assignment: activation (modules 1-3) vs repression (module 4)

### Temporal dynamics analysis
- Z-score normalization of expression across timepoints
- Phase ratio: proportion of total expression change occurring at T1→T2
- Early responder threshold: phase ratio > 0.6
- Mann-Whitney U test for phase separation between blocs
- Module-internal temporal coherence (Spearman rho)

### TF family enrichment analysis
- TF family annotation (TetR, TCS, sigma, WhiB, GntR, IclR, LacI, etc.)
- Fisher's exact test for family enrichment in activation vs repression bloc
- TCS pair analysis: exposed/shielded status of cognate SK-RR pairs
- Binomial test for asymmetry

### Neighborhood analysis
- Flanking gene expression change (±10 genes from each exposed TF)
- Permutation test (10,000 permutations): observed vs random mean |LFC|
- Distance decay: Spearman correlation between genomic distance and |LFC|
- Concordance rate: fraction of neighbors with same direction LFC

### Conservation analysis
- SCO locus tag presence/absence (proxy for conservation in S. coelicolor reference)
- GC3 (third codon position GC content)
- Effective number of codons (Nc)
- Rare codon frequency
- Composite conservation score: combined Z-score
- ROC analysis for conservation predicting exposed status

### Statistical analysis
- All analyses performed in Python 3.x with scipy, statsmodels, scikit-learn
- Multiple testing correction: Benjamini-Hochberg where applicable
- p-value reporting: scientific notation for p < 0.001, decimal for p >= 0.001
- Effect sizes reported alongside p-values throughout
- Visualization: matplotlib, seaborn

---

## Figures and Tables

### Main Figures
| # | Title | Script |
|---|-------|--------|
| Fig. 1 | Methylation Landscape overview (Circos, site counts, genomic distribution, logos) | `01_figure1_landscape.py` |
| Fig. 2 | Methylation-exclusion architecture overview | `10_new_figure1_overview.py` |
| Fig. 3 | Protection zone characterization | `11_new_figure2_protection.py` |
| Fig. 4 | 57 exposed TFs — distributed developmental switch | `12_new_figure3_exposed_TFs.py` |
| Fig. 5 | Negative results panel | `13_new_figure4_negative_results.py` |
| Fig. 6 | Vegetative-to-developmental switch model | `14_new_figure5_switch_model.py` |

### Supplementary Figures
| # | Content |
|---|---------|
| Fig. S1 | 4mC reclassification + AAGCCCG novel motif characterization |
| Fig. S2 | Simpson's paradox detailed (6-panel: unstratified/core/arm × GCCGGC/AAGCCCG) |
| Fig. S3 | Regulatory avoidance forest plot (CMH-adjusted) |
| Fig. S4 | Sequence-level motif depletion (gene body + promoter) |
| Fig. S5 | 57 exposed TF complete annotation |
| Fig. S6 | TCS pair detailed analysis |
| Fig. S7 | Conservation metrics (H36) |
| Fig. S14 | R-M recognition motif characterization: methylated position in AAGCCCG (pos1=63%) and GCCGGC (palindromic pos1/pos3), motif occupancy (AAGCCCG=31.6%/GCCGGC=7.96%/GATC=0.27%), 6mA Unassigned fraction (72.2% at T1) — `24_figS14_motif_reliability.py` |
| Fig. S15 | Genome-wide methylation density dynamics (50 kb windows): 4mC T2-peak/T3-decline, 6mA monotonic increase, core/arm ratio dynamics — `25_figS15_window_methylation.py` |
| Fig. S16 | Timepoint-resolved TSS protection zone: 4mC regulatory gene protection collapses at T3 (ratio=1.046), 6mA deepest protection at T2, non-regulatory genes show stable protection — `26_figS16_timepoint_TSS.py` |

### Supplementary Tables
| # | Content |
|---|---------|
| Table S1 | All high-confidence methylation sites (position, mod_type, motif, frequency, timepoint) |
| Table S2 | Motif summary + REBASE conservation |
| Table S3 | MTase genes (22 genes × expression/domain) |
| Table S4 | 1,055 regulatory genes: shielded/exposed classification + all features |
| Table S5 | 57 exposed TFs: full annotation (family, bloc, module, temporal class, LFC, conservation) |
| Table S6 | TCS pair analysis (7 pairs × exposed/shielded status) |
| Table S7 | Simpson's paradox statistics (3 instances × unstratified/stratified) |
| Table S8 | Complete hypothesis ledger (H1-H36 × verdict, key statistic) |

---

## Drafting Notes

### Key narrative arc
1. **Hook**: Bacterial methylation is assumed to regulate transcription, but is this true genome-wide?
2. **Surprise**: Apparent correlations are Simpson's paradox artifacts
3. **Discovery**: Instead of genome-wide regulation, a structural dichotomy partitions the regulatory genome
4. **Mechanism**: The protection zone is created by collective protein occupancy, not individual TF binding
5. **Biology**: The 57 exposed regulators encode a coordinated developmental switch
6. **Impact**: Reframes bacterial epigenetics from "methylation controls genes" to "methylation partitions regulatory architecture"

### Tone
- Confident but measured: strong evidence for structural findings, honest about mechanistic unknowns
- Negative results presented as positive contributions (Simpson's paradox, expression-independence, no neighborhood effect)
- Avoid overclaiming causation from correlational data
- Emphasize the iterative, hypothesis-driven discovery process in Discussion

### Word count targets
| Section | Target | Notes |
|---------|--------|-------|
| Abstract | 250 | Structured: background/results/conclusion |
| Introduction | 1,500 | 6 paragraphs |
| Results | 4,500 | 18 paragraphs, organized by Fig. 1-6 |
| Discussion | 2,000 | 7 paragraphs |
| Methods | 2,500 | 14 subsections |
| **Total** | **~10,750** | Excluding references |

### Cross-reference checklist
- [ ] All statistics in text match `99_style_guide.md` Key Statistics
- [ ] All figure citations match the 6 main + 7 supplementary figure structure
- [ ] Terminology consistent with style guide Section 1
- [ ] All abbreviations defined at first use (style guide Section 5)
- [ ] p-value formatting follows style guide Section 6
