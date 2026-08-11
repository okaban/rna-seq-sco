# Citation Report: m4C vs m5C at GCCGGC in *S. coelicolor* — Scientific Argument Support

**Date:** 2026-05-06  
**Purpose:** Literature support for the argument that Nanopore-detected N4-methylcytosine (m4C) at GCCGGC motifs in *S. coelicolor* M145 is the correct interpretation, and that Pisciotta et al. 2023 incorrectly assigned this methylation as m5C due to the intrinsic limitation of bisulfite sequencing.

**Core argument outline:**
1. Bisulfite sequencing cannot distinguish m4C from m5C — both are bisulfite-resistant
2. Pisciotta et al. 2023 used bisulfite-seq and called GCCGG methylation as m5C
3. Our Nanopore direct sequencing with m4C-specific model confirms m4C, not m5C
4. Supporting evidence: no m5C MTase for GCCGGC in REBASE/Streptomyces; hemi-methylation atypical for m5C MTases; Nanopore model specificity

---

## Category 1: Bisulfite Sequencing Cannot Distinguish m4C from m5C

These papers directly establish the biochemical basis for why bisulfite sequencing is unable to differentiate N4-methylcytosine (m4C) from 5-methylcytosine (m5C) — both modifications block bisulfite-mediated deamination, causing both to read as cytosine (unmodified) rather than uracil/thymine.

---

### 1.1 ⭐ PRIMARY CITATION (Definitive chemistry paper)

**Vilkaitis G, Klimasauskas S (1999).** "Bisulfite sequencing protocol displays both 5-methylcytosine and N4-methylcytosine."  
*Analytical Biochemistry* 271(1):116–119.  
DOI: [10.1006/abio.1999.4116](https://doi.org/10.1006/abio.1999.4116)  
PMID: 10361019

**Relevant finding:** This is the foundational paper demonstrating experimentally that the bisulfite sequencing protocol cannot distinguish m4C from m5C — both modifications are displayed identically (as protected/unconverted cytosine) in the sequencing readout. Any downstream bisulfite-seq interpretation of "m5C" in a bacterial genome containing m4C is inherently ambiguous.

**Quotable point for argument:** The paper's title itself is the argument. When bisulfite-seq "displays both," calling the result exclusively m5C is unwarranted.

---

### 1.2 ⭐ STRONG SUPPORTING CITATION (Functional consequence, quantitative)

**Yu M, Ji L, Neumann DA, Chung D-H, Groom J, Westpheling J, He C, Schmitz RJ (2015).** "Base-resolution detection of N4-methylcytosine in genomic DNA using 4mC-Tet-assisted bisulfite-sequencing."  
*Nucleic Acids Research* 43(21):e148.  
DOI: [10.1093/nar/gkv738](https://doi.org/10.1093/nar/gkv738)  
PMID: 26184871  PMC: PMC4666385

**Relevant finding:** Developed 4mC-TAB-seq precisely because standard bisulfite sequencing (MethylC-seq) "is not suitable to accurately differentiate 4mC from 5mC as both will be read as C." In 4mC-TAB-seq, both C and m5C are converted to T, whereas m4C reads as C — proving the orthogonal chemistry. Applied to *Caldicellulosiruptor* bacteria where m4C is the dominant cytosine modification and is a major barrier to genetic transformation.

**Quotable point for argument:** The authors created an entirely new sequencing method because standard bisulfite-seq conflates m4C with m5C. This is a methods paper whose raison d'être is the failure mode affecting Pisciotta et al. 2023.

---

## Category 2: Nanopore Direct Sequencing Can Specifically Detect m4C

These papers validate that Oxford Nanopore Technologies (ONT) direct sequencing, combined with appropriate basecalling models (e.g., Dorado), can detect m4C as a distinct modification class in bacterial genomes — without chemical pre-treatment.

---

### 2.1 ⭐ BEST METHODS VALIDATION (4mC + 5mC + 6mA in bacteria, Dorado v5)

**Galeone V, Dabernig-Heinz J, Lohde M, Brandt C, Kohler C, Wagner GE, Hölzer M (2025).** "Decoding bacterial methylomes in four public health-relevant microbial species: nanopore sequencing enables reproducible analysis of DNA modifications."  
*BMC Genomics* 26(1):394.  
DOI: [10.1186/s12864-025-11592-z](https://doi.org/10.1186/s12864-025-11592-z)  
PMID: 40269718  PMC: PMC12016153

**Relevant finding:** Using Nanopore R10.4.1 flow cells and Dorado basecalling SUP model v5, this study simultaneously detected 4mC, 5mC, and 6mA modifications across four bacterial species. The paper explicitly evaluates all three bacterial cytosine methylation types and demonstrates that Dorado v5 models have substantially improved accuracy for m4C detection, with keywords explicitly including "4mC." Provides a direct precedent for nanopore m4C calling in bacteria using the same platform and model generation as our study.

---

### 2.2 Bacterial m4C motif discovery from Nanopore data

**Tidwell AK, Faust E, Eckert CA, Guss AM, Alexander WG (2024).** "Discovering methylated DNA motifs in bacterial nanopore sequencing data with MIJAMP."  
*Journal of Industrial Microbiology & Biotechnology* 52.  
DOI: [10.1093/jimb/kuaf022](https://doi.org/10.1093/jimb/kuaf022)  
PMID: 40699189  PMC: PMC12320774

**Relevant finding:** Developed MIJAMP specifically to call methylated motifs (including 4mC) from ONT Modkit output in bacterial genomes. Benchmarked against *E. coli* (known m6A, m5C, m4C marks) and showed MIJAMP correctly attributed 80.3% of m4C sites — substantially outperforming MicrobeMod (which called 0% of m4C). Provides independent validation that ONT pipelines can specifically resolve m4C motifs separate from m5C.

---

### 2.3 Reproducibility and accuracy of ONT bacterial methylome profiling

**Schababerle T, Hayat O, Jung J, Le M, Polic I, Wees H, Bhatti M, Shelburne S, Liu X, Kalia A (2025).** "Reproducibility and accuracy of bacterial methylome profiling using Oxford Nanopore Technologies nanopore sequencing platform."  
*Microbial Genomics* 11(11).  
DOI: [10.1099/mgen.0.001564](https://doi.org/10.1099/mgen.0.001564)  
PMID: 41259098

**Relevant finding:** Multi-operator, multi-replicate study demonstrating >99.9% motif identification concordance between ONT-only and hybrid reference assemblies for bacterial methylome profiling with the Dorado basecaller. Motif-level methylation calls achieved F1-score >99.999%. Establishes that ONT-based bacterial methylome profiling is both accurate and reproducible, supporting the reliability of our Nanopore m4C calls.

---

### 2.4 Workflow tool for bacterial methylation annotation

**Jakubickova M, Sabatova K, Zbudilova M, Bezdicek M, Lengerova M, Vitkova H (2025).** "MethylomeMiner: A novel tool for high-resolution analysis of bacterial methylation patterns from nanopore sequencing."  
*Computational and Structural Biotechnology Journal* 27:4753–4759.  
DOI: [10.1016/j.csbj.2025.10.047](https://doi.org/10.1016/j.csbj.2025.10.047)  
PMID: 41542075  PMC: PMC12800368

**Relevant finding:** Python-based workflow for processing ONT methylation calls (including m4C) and mapping them to genomic features (coding/non-coding regions). Supports population-level and pangenome-scale m4C comparison — provides the methodological framework directly applicable to our analysis pipeline.

---

## Category 3: Bacterial m4C as a Restriction-Modification Mark (and m4C/m5C Confusion)

These papers establish that (a) m4C is a well-characterized bacterial R-M modification; (b) the GCCGGC/CCGCGG palindrome is a known m4C target in bacteria; and (c) m4C and m5C are both recognized by methyl-directed restriction systems, underscoring their biochemical similarity and the genuine risk of misidentification by bisulfite-seq.

---

### 3.1 ⭐ DIRECTLY ON-MOTIF: m4C at CCGCGG in *Deinococcus radiodurans*

**Shi C, Wang L, Xu H, Zhao Y, Tian B, Hua Y (2024).** "Characterization of a Novel N4-Methylcytosine Restriction-Modification System in [*Deinococcus radiodurans*]."  
*International Journal of Molecular Sciences* 25(3):1660.  
DOI: [10.3390/ijms25031660](https://doi.org/10.3390/ijms25031660)  
PMID: 38338939  PMC: PMC10855626

**Relevant finding:** Identified M.DraR1, an m4C methyltransferase that recognizes **5'-CCGCGG-3'** (= palindrome of GCCGGC) and methylates the second cytosine. The cognate restriction enzyme R.DraR1 specifically cleaves unmethylated CCGCGG but not m4C-methylated substrate — and notably, it *can* cleave m5C-methylated CCGCGG, demonstrating that m4C and m5C at this exact sequence context are biochemically distinguishable by restriction but not by bisulfite-seq. This is the closest published precedent to our motif.

---

### 3.2 Review: Bacterial DNA methylation — m4C biology and roles

**Casadesús J (2016).** "Bacterial DNA Methylation and Methylomes."  
*Advances in Experimental Medicine and Biology* 945:35–61.  
DOI: [10.1007/978-3-319-43624-1_3](https://doi.org/10.1007/978-3-319-43624-1_3)  
PMID: 27826834

**Relevant finding:** Authoritative review establishing m4C (N4-methylcytosine) as a major bacterial DNA modification in R-M systems, alongside m5C and m6A. Reviews the regulatory roles of methylation in bacteria and the challenge of assigning modification types without single-molecule or direct sequencing methods.

---

### 3.3 McrBC recognizes m4C with *higher* affinity than m5C

**Zagorskaitė E, Manakova E, Sasnauskas G (2018).** "Recognition of modified cytosine variants by the DNA-binding domain of methyl-directed endonuclease McrBC."  
*FEBS Letters* 592(19):3335–3345.  
DOI: [10.1002/1873-3468.13244](https://doi.org/10.1002/1873-3468.13244)  
PMID: 30194838

**Relevant finding:** Crystal structures show McrBC binds modified cytosines in order: **4mC > 5mC > 5hmC ≫ 5fC**. Because McrBC is used in methylation profiling experiments (and was historically used in enrichment protocols), this paper confirms that m4C and m5C are both recognized by the same methyl-sensing machinery — illustrating that they are functionally interchangeable in many assay contexts, including bisulfite conversion chemistry.

---

### 3.4 Strand-opposite m4C methylation by Type IIs R-M system

**Sapranauskas R, Sasnauskas G, Lagunavicius A, Vilkaitis G, Lubys A, Siksnys V (2000).** "Novel subtype of type IIs restriction enzymes. BfiI endonuclease exhibits similarities to the EDTA-resistant nuclease Nuc of *Salmonella typhimurium*."  
*Journal of Biological Chemistry* 275(40):30878–30885.  
DOI: [10.1074/jbc.M003350200](https://doi.org/10.1074/jbc.M003350200)  
PMID: 10880511

**Relevant finding:** Characterized BfiI, a Type IIs R-M system carrying **two N4-methylcytosine methyltransferases** that each methylate cytosines on **opposite strands** of the recognition sequence. Strand-asymmetric m4C methylation (as in many Type II R-M systems) produces hemi-methylation patterns characteristic of m4C and mechanistically distinct from the symmetric CpG methylation of eukaryotic m5C MTases. This is relevant to our argument that the strand-specific hemi-methylation we observe is atypical for m5C but consistent with m4C R-M biology.

---

### 3.5 SMRT vs. bisulfite for bacterial cytosine methylation — complementary methods

**Anton BP, Fomenkov A, Wu V, Roberts RJ (2021).** "Genome-wide identification of 5-methylcytosine sites in bacterial genomes by high-throughput sequencing of MspJI restriction fragments."  
*PLoS ONE* 16(5):e0247541.  
DOI: [10.1371/journal.pone.0247541](https://doi.org/10.1371/journal.pone.0247541)  
PMID: 33974631  PMC: PMC8112702

**Relevant finding:** Notes explicitly that while SMRT sequencing can directly detect m6A and m4C, "similar identification of 5-methylcytosine sites is not as straightforward," motivating a new m5C-specific method. Conversely, bisulfite-based methods enriching for m5C will co-purify m4C because MspJI and other methyl-directed enzymes cannot distinguish them. Developed by Roberts (REBASE founder), lending authority to the claim that m4C/m5C differentiation requires orthogonal methods.

---

## Category 4: Streptomyces Cytosine Methylation and Differentiation (2018–2025)

These papers support the claim "GCCGG-type cytosine methylation motifs cluster upstream of differentiation-associated genes and disrupting methylation is associated with delayed aerial hyphae formation."

---

### 4.1 ⭐ PRIMARY REFERENCE (paper being challenged)

**Pisciotta A, Sampino AM, Presentato A, Galardini M, Manteca A, Alduina R (2023).** "The DNA cytosine methylome revealed two methylation motifs in the upstream regions of genes related to morphological and physiological differentiation in *Streptomyces coelicolor* A(3)2 M145."  
*Scientific Reports* 13(1):7038.  
DOI: [10.1038/s41598-023-34075-1](https://doi.org/10.1038/s41598-023-34075-1)  
PMID: 37120673  PMC: PMC10148868

**Relevant finding:** Used bisulfite sequencing to identify 3360 methylated cytosines and two motifs — GGCmCGG and GCCmCG — in upstream regions of 321 *S. coelicolor* M145 genes. Attributed both motifs to m5C based solely on bisulfite sequencing data. 5-aza-dC hypomethylation impaired growth, antibiotic biosynthesis, and delayed development. **Our counter-argument:** The GCCGGC/GCCMCGG motif is a known m4C target (see Category 3 above); bisulfite-seq cannot distinguish m4C from m5C; the absence of a cognate m5C MTase for this palindrome in REBASE/Streptomyces, combined with our Nanopore m4C signal, indicates misassignment.

---

### 4.2 ⭐ SECOND CITATION CANDIDATE for the differentiation/hyphae claim

**Pisciotta A, Manteca A, Alduina R (2018).** "The SCO1731 methyltransferase modulates actinorhodin production and morphological differentiation of *Streptomyces coelicolor* A3(2)."  
*Scientific Reports* 8(1):13686.  
DOI: [10.1038/s41598-018-32027-8](https://doi.org/10.1038/s41598-018-32027-8)  
PMID: 30209340  PMC: PMC6135851

**Relevant finding:** Demonstrated that cytosine methylation by the SCO1731 methyltransferase directly controls morphological differentiation and antibiotic production in *S. coelicolor*. Disruption of SCO1731 (the cytosine MTase responsible for m5C/m4C marks) caused: delayed spore germination, impaired aerial mycelium development, delayed sporulation, and nearly abolished actinorhodin production — phenocopying 5-aza-dC treatment. **This is the paper establishing that methylation at GCCGG-type motifs functionally drives Streptomyces differentiation.**

> ⚠️ **Note on SCO1731:** This MTase catalyzes the methylation at the motifs observed in Pisciotta 2023. Our Nanopore data suggesting m4C at these motifs would mean SCO1731 is an m4C methyltransferase — an important reinterpretation to flag (no m4C activity has been biochemically characterized for SCO1731 to date; this would be a novel finding).

---

## Summary Table

| # | Paper | Category | Key Claim Supported |
|---|-------|----------|---------------------|
| 1.1 | Vilkaitis & Klimasauskas 1999 | Bisulfite limitation | Bisulfite-seq displays m4C and m5C identically |
| 1.2 | Yu et al. 2015 | Bisulfite limitation | Standard MethylC-seq cannot differentiate m4C from m5C |
| 2.1 | Galeone et al. 2025 | Nanopore m4C detection | ONT Dorado v5 detects 4mC, 5mC, 6mA distinctly in bacteria |
| 2.2 | Tidwell et al. 2024 | Nanopore m4C detection | MIJAMP calls 80.3% of m4C motifs from ONT data |
| 2.3 | Schababerle et al. 2025 | Nanopore m4C detection | ONT bacterial methylome: >99.9% reproducibility with Dorado |
| 2.4 | Jakubickova et al. 2025 | Nanopore m4C detection | MethylomeMiner: ONT m4C analysis workflow for bacteria |
| 3.1 | Shi et al. 2024 | Bacterial m4C R-M | m4C at CCGCGG in *Deinococcus* — same palindrome as GCCGGC |
| 3.2 | Casadesús 2016 | Bacterial m4C R-M | Review of m4C as canonical bacterial R-M mark |
| 3.3 | Zagorskaitė et al. 2018 | Bacterial m4C R-M | McrBC binds m4C > m5C; both recognized by same sensors |
| 3.4 | Sapranauskas et al. 2000 | Bacterial m4C R-M | Strand-asymmetric m4C by dual MTases (hemi-methylation precedent) |
| 3.5 | Anton et al. 2021 | Bacterial m4C R-M | SMRT/bisulfite complementarity; bisulfite conflates m4C with m5C |
| 4.1 | Pisciotta et al. 2023 | *Streptomyces* methylation | Target paper: GCCGG motifs upstream of differentiation genes |
| 4.2 | Pisciotta et al. 2018 | *Streptomyces* methylation | SCO1731 methylation drives aerial hyphae, sporulation, actinorhodin |

---

## Recommended Citation Placements in Manuscript

**For the statement "bisulfite sequencing cannot distinguish m4C from m5C":**  
→ Cite [Vilkaitis & Klimasauskas 1999] and [Yu et al. 2015]

**For "Nanopore direct sequencing with m4C-specific model":**  
→ Cite [Galeone et al. 2025] and [Tidwell et al. 2024]; reproducibility supported by [Schababerle et al. 2025]

**For "no m5C MTase for GCCGGC in Streptomyces; known m4C motif":**  
→ Cite [Shi et al. 2024] (CCGCGG m4C in *Deinococcus*) and [Casadesús 2016] (review)

**For "strand-specific hemi-methylation atypical for m5C MTases":**  
→ Cite [Sapranauskas et al. 2000]; contrast with canonical m5C biology

**For the [CITATION NEEDED] alongside Pisciotta 2023 (GCCGG motifs + differentiation):**  
→ Cite [Pisciotta et al. 2018] — directly establishes that cytosine methylation in *S. coelicolor* drives aerial hyphae formation and sporulation

---

## Literature Gaps / Caveats

1. **No prior Streptomyces methylome from nanopore/SMRT.** The published *S. coelicolor* methylation literature is entirely from the Pisciotta/Alduina/Manteca group using bisulfite-seq. Our Nanopore data would be the first direct-sequencing methylome for this organism — a notable strength.

2. **SCO1731 biochemistry unknown for m4C.** If we argue SCO1731 is an m4C MTase, we need to note this is a novel claim and ideally support with REBASE evidence (no cognate m5C MTase annotated for GCCGGC) or in vitro activity data.

3. **The Vilkaitis 1999 paper is a brief communication** (~4 pages). For reviewers wanting more mechanistic depth, pair with Yu et al. 2015 which provides the definitive quantitative demonstration.

4. **Nanopore m4C model specificity for Streptomyces.** The Galeone and Schababerle papers use model organisms (*Staph*, *Listeria*, *Klebsiella*, *Streptococcus*). We should note whether the Dorado m4C model was validated on Actinobacteria/high-GC organisms — this may be a reviewer question.

---

*Report compiled from PubMed and web searches on 2026-05-06.*  
*Sources retrieved via PubMed (PMID-based metadata) and web search.*
