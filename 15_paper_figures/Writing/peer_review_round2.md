# Peer Review — Round 2
**Manuscript:** S. coelicolor A3(2) M145 DNA Methylation × RNA-seq Epigenome Integration  
**Review date:** 2026-05-11  
**Target journals reviewed against:** Nature Microbiology / EMBO Journal / Nucleic Acids Research  
**Reviewer stance:** Independent expert reviewer; no prior knowledge of authors' work assumed.

---

## Executive Summary

This manuscript presents a genome-scale analysis of two DNA methylation systems (GCCGGC 4mC and AAGCCCG 6mA) across three developmental stages of *S. coelicolor* M145, integrated with matched RNA-seq. The central methodological contribution — demonstrating that apparent methylation-expression correlations are Simpson's paradox artifacts driven by the linear chromosome's core/arm architecture — is genuinely novel and methodologically important. The four-layer "methylation-exclusion" model and the Shielded/Exposed dichotomy are conceptually original. However, the manuscript carries serious weaknesses that would prevent acceptance at Nature Microbiology in its current form, and require substantial revision even for EMBO Journal. The most critical issues are: (1) a numerical inconsistency between Introduction and the rest of the paper (5.9% vs 5.4%); (2) the unresolved status of 72.2% of 6mA sites that undermine the "two independent methylation systems" claim; (3) the complete absence of genetic perturbation or direct mechanistic validation; (4) Fisher OR numerical inconsistencies across manuscript sections; and (5) residual editorial artefacts (development comment tags) that must be removed before submission.

---

## Major Concerns

### MC1: Internal numerical inconsistency — 5.4% vs 5.9% exposed regulator fraction

The proportion of exposed regulatory genes is stated as **5.4%** (57/1,055) in the Abstract, Results, and Discussion — which is arithmetically correct. However, the Introduction (paragraph 4) states "the remaining **5.9%** of regulatory genes — 57 'exposed regulators'." This discrepancy is not trivial: 5.9% would imply a denominator of ~966, not 1,055. This inconsistency suggests the Introduction was drafted at an earlier stage of the analysis when the regulatory gene set had a different size, and was never reconciled with the finalized 1,055-gene count. A Nature Microbiology reviewer encountering this on a first read will immediately doubt the internal consistency of the entire numerical framework. **This must be corrected throughout before submission.** In the same Introduction paragraph, the phrase "the remaining 5.9% of regulatory genes" should read "the remaining 5.4% (57/1,055)."

### MC2: The "two independent methylation systems" claim is undermined by 72.2% unassigned 6mA sites

The manuscript's Big Claim 1 — that *S. coelicolor* M145 harbors "two independent methylation systems" (GCCGGC 4mC and AAGCCCG 6mA) — is prominently stated in the Abstract and Introduction. However, the Results section reveals that **72.2% of 6mA sites at T1 (1,396/1,934) do not match AAGCCCG, GATC, or their reverse complements within a ±15 bp window**, a fraction increasing to 80.3% at T2 and 82.7% at T3. STREME analysis failed to identify a consensus motif for these unassigned sites (best E = 0.12). This means the 6mA landscape is dominated by modification from at least one, possibly more, **entirely uncharacterized MTase system(s)**, accounting for the majority of observed 6mA sites at every timepoint.

This has two major consequences. First, the claim of "two methylation systems" is misleading — the data reveal at minimum three, with the third being uncharacterized and dominant. Second, the internal control logic relies on AAGCCCG having "minimal geographic bias" to argue that geographic stratification is unnecessary for the 6mA de-repression test. But the 6mA de-repression test uses **all** 6mA sites, the majority of which come from the uncharacterized system. If that system has non-trivial geographic bias, the internal control argument fails. The authors must either (a) reframe the paper to explicitly acknowledge three-or-more methylation systems and restrict claims accordingly, or (b) demonstrate that the unassigned 6mA sites also lack geographic bias and do not confound the null result.

### MC3: SC_RS17645 as "solo MTase" — Type I architecture inconsistency and direct biochemical evidence absent

The figure_claims_map and the Introduction present SC_RS17645 (SCO3104) as a "solo MTase" for AAGCCCG 6mA, and Figure 1C is described as providing "structural evidence" for a solo MTase. However, the Methods section identifies SC_RS17645 as an **HsdM-type** methyltransferase. In canonical Type I R-M systems, HsdM functions as part of a heteromeric complex with HsdR (restriction subunit) and HsdS (specificity/target-recognition subunit) — a solo HsdM is not the default assumption. The manuscript never addresses whether cognate HsdR and HsdS subunits exist in M145. If they do, this is not a solo MTase but part of a Type I system; if they do not, this needs explicit genomic evidence (absence of RM-associated HsdR/HsdS in the genomic neighborhood). 

Furthermore, the attribution of AAGCCCG modification to SC_RS17645 rests solely on: (1) motif co-occurrence with HsdM domain annotation, and (2) correlated expression downregulation (LFC = −2.19 at T2 vs T1). The 4.6-fold downregulation of SC_RS17645 at T2 corresponds to a 75% loss of AAGCCCG sites — this correlation is suggestive but not probative; 20 other MTase genes are expressed in M145 and could contribute. The manuscript correctly labels SC_RS17645 as a "candidate" in the Results section, but the figure_claims_map and Introduction use the term "solo MTase" without qualification — language that will attract reviewer criticism at a Nature Microbiology-level journal. Either provide in vitro methyltransferase activity data (even cell-free assay), or consistently use "candidate" language and remove "solo MTase" from all presentation-level descriptions. Verification by SC_RS17645 knockout and rescue would be the definitive experiment.

### MC4: Fisher OR numerical inconsistency across manuscript sections

The co-modification Fisher's exact odds ratio is reported as **OR = 138,440** (95% CI 34,383–557,000) in the Methods section, but as **OR = 145,448** in the figure_claims_map. The CLAUDE.md project file reports "OR=921×" as the headline co-modification figure, which refers to the hypergeometric fold-enrichment (not the Fisher OR). These three numbers — 138,440, 145,448, and 921 — measure related but distinct quantities and are currently presented in ways that risk conflation. A reviewer who reads the figure legend (presumably showing 145,448) and then the Methods (138,440) will flag this as an error. The correct Fisher OR from `statsmodels` must be verified against a single authoritative run of the analysis, the discrepancy resolved, and a clear distinction drawn between (a) the 921× fold-enrichment over permutation null and (b) the Fisher OR. These should appear in different contexts and never be interchanged.

### MC5: AUC = 0.917 uses pooled all-timepoint methylation — a potential circularity

The Shielded/Exposed classification (AUC = 0.917, threshold = 293 bp) is computed using the nearest methylation distance from each regulatory gene TSS, where "nearest site" is taken across **all timepoints and both modification types combined** (Methods: "across all timepoints and types"). This pooling decision is not trivial. At T3, only 21 GCCGGC sites remain genome-wide — an average of less than one site per 400 kb. If T3 sites alone were used, essentially all genes would appear "shielded" by the 293 bp criterion simply due to sparsity. The high AUC may therefore partly reflect the T1 methylation pattern, where site density is highest. The manuscript should report the AUC computed at each timepoint independently (T1-only, T2-only, T3-only) to determine whether the Shielded/Exposed dichotomy is stable across developmental stages or is driven by the high-density T1 dataset. If AUC degrades substantially at T2 or T3, the biological interpretation changes.

### MC6: No genetic perturbation — the mechanism is unvalidated and the protection zone is inferred

The manuscript makes strong structural claims: that a ~2,200 bp methylation-free protection zone is maintained by "collective protein occupancy at promoters" (NAPs + basal transcription apparatus), and that this protection is structural rather than transcription-dependent. The sole evidence for protein occupancy is the absence of an alternative explanation (expression-independence, individual TF BS non-protection) and a citation to existing literature on NAP distributions. No ChIP-seq, DAP-seq, DNase-seq, ATAC-seq, or any direct occupancy measurement is presented. The argument is: "the gap between AUC = 0.917 and AUC = 0.712 (sequence-only) = 0.205 = protein occupancy." This is a residual attribution — attributing unexplained variance to a favored hypothesis — which does not constitute evidence for that hypothesis.

Similarly, the claim that 57 exposed regulators constitute "a parallel regulatory layer" is purely descriptive: co-expression is demonstrated (within-bloc ρ = 0.717), but the coordinating mechanism is explicitly unknown (0/57 have curated binding motifs; upstream TF targeting fails; cis-regulation excluded). Nature Microbiology will require that at minimum one or two of the top exposed regulators (e.g., *ramR*, *nsdB*) be validated by directed genetics (knockout or CRISPRi) demonstrating that their exposed status is functionally relevant. The Discussion's Limitations section is admirably honest about this gap, but honesty about a gap does not substitute for closing it.

### MC7: Residual editorial artefacts must be removed before submission

The Results section contains multiple `<!-- E14: ... -->` HTML comment tags that are clearly development notes for the authors, not for publication. These include:
- "E14: Candidate for SuppFig (temporal decoupling パネル, ρ=0.136 NS)" 
- "E14: Fig. S9 → 正式 Supplementary Figure XX 番号に更新予定"
- References to figure numbers as "Supplementary Figure XX" (placeholder)
- A Japanese-language author note embedded in Results prose

These tags indicate that figure numbering is still unresolved (some figures are planned to move from main to supplementary), and at least two panel references (Fig. S9a/S9b) in the main Results text need updating. Submission with these tags present would result in immediate editorial rejection. A complete pass to remove all comment tags, resolve all figure numbers, and audit all cross-references is mandatory before submission. The same applies to "Fig. 3e" which is used to label two distinct panels (expression quintile analysis and TCS asymmetry figure) in the Results section.

---

## Minor Concerns

**Causal framing drift between sections.** The Abstract's concluding sentence — "DNA methylation does not directly control transcription in *S. coelicolor*" — is a strong negative causal claim. The evidence is correlational (absence of correlation). The manuscript cannot rule out indirect methylation effects (e.g., methylation-dependent protein recruitment that affects transcription indirectly). The Discussion qualifies this appropriately ("observational, correlational evidence"), but the Abstract statement is stronger than the data support. Recommend revising to: "Our data provide no evidence that DNA methylation directly controls transcription in *S. coelicolor*."

**"36 hypotheses tested across 18 iterative analysis cycles" in Introduction.** This phrasing reads as a research diary entry rather than a scientific statement. Reviewers may interpret it as post-hoc hypothesis generation, which raises multiple-testing concerns. The claim that "12 rejected hypotheses were required to reveal the true pathway" is not an argument for the model's validity — it is a description of the discovery process that is better left to supplementary notes or removed entirely.

**The "Gatekeeper" model name is introduced in the Introduction but absent from the Abstract and Results.** Terminological inconsistency across sections will confuse readers. Either use "Gatekeeper" consistently throughout (Abstract → Introduction → Results → Discussion) or drop it and use "four-layer methylation-exclusion architecture" uniformly.

**TCS 7/7 asymmetry with p = 0.83.** The Discussion devotes substantial prose to the observation that all 7 identified cognate TCS pairs contain exactly one exposed member, followed by a permutation test yielding p = 0.83 (not significant, n = 7 pairs). The paper correctly concludes this is statistically indistinguishable from chance. However, presenting a non-significant observation at this length, with a biological interpretation ("differentially gating input and output arms"), risks being read as data dredging. We recommend either moving this to a supplementary note or condensing it to two sentences: the observation, the permutation p-value, and a statement that this awaits testing with larger TCS datasets.

**TetR enrichment p = 0.028 without multiple testing correction.** The Results and Discussion report TetR-family enrichment in the repression bloc at p = 0.028 (Fisher's exact). Given that multiple TF families were tested for enrichment, a multiple-testing correction (Bonferroni or BH-FDR) should be applied. If TetR enrichment does not survive correction, the claim that the repression bloc shows "statistically significant enrichment for TetR family" (p = 0.028) must be qualified accordingly.

**One sample at 14.3x sequencing coverage.** The Methods report one of nine Nanopore samples at 14.3x coverage ("eight samples at 27–52x; one sample at 14.3x"). For reliable methylation calling with modkit, 14x is at the lower boundary. The paper should report (a) which sample this is, (b) whether excluding it or down-weighting it changes any key result, and (c) what modkit's minimum recommended coverage is for 4mC/6mA calling confidence at this threshold.

**The protection zone "deepest depletion at +300 bp."** The Results state "the deepest depletion occurred at +300 bp downstream of the TSS, where methylation density dropped to 54.1% of flanking region levels." A depletion maximum *inside* the gene body (+300 bp) rather than upstream of the TSS is biologically unusual and warrants mechanistic comment. Is this an artifact of TSS position assignment errors (many TSSs assigned from RefSeq start codons rather than dRNA-seq), or does it reflect a real biological pattern (methylation exclusion from the transcription initiation complex extending into the early transcribed region)? The 75% of TSSs assigned from RefSeq annotation (not dRNA-seq) will systematically overestimate the coding-sequence start relative to the true TSS, potentially placing the "protection zone" incorrectly. The supplementary analysis using only the 257 dRNA-seq-derived TSSs (AUC = 0.923, 318 bp) is reassuring but should be presented as a primary sensitivity analysis in the main text, not buried in supplementary notes.

**OR = ∞ for constitutive expression.** "None of the 57 exposed regulators showed constitutive expression (0%), compared to 38% of shielded regulators (odds ratio = infinity, p = 5.3 × 10^-13)." OR = ∞ arises from complete separation (zero cells in a Fisher 2×2 table). While Fisher's exact p is valid, the OR = ∞ notation should be replaced with a statement of the observed cell counts and a note that the OR is undefined but the result is Fisher-exact significant at p < 10^-12.

**Missing comparator: prior *Streptomyces* epigenomics literature.** The Introduction and Discussion cite Dam, CcrM, and generically "genome-wide surveys" (Blow et al., Beaulaurier et al.) but do not engage with any prior *Streptomyces*-specific methylation or epigenomics work. Is this truly the first genome-wide methylation × expression study in *Streptomyces*? If so, this should be stated explicitly. If not, prior work should be cited and positioned. REBASE entries for *S. coelicolor* M145 should also be cited to contextualize the two identified R-M systems within the predicted RM repertoire.

**Sequence-AUC decomposition framing inconsistency.** The Results report the sequence component of protection as "approximately one-quarter to one-third" (from motif depletion metrics) or "~51%" (from AUC decomposition). The Discussion's two-tier model figure description states "~33% sequence, ~67% protein occupancy" in one place and "~51%/~49%" in another. The acknowledgment that "two metrics measure different things" is appropriate, but the inconsistent framing within a single paper is confusing. Recommend standardizing: present the AUC-decomposition estimate (51%/49%) as the primary metric (since it relates directly to the classifier performance being claimed), with the motif-depletion figure (33%) in a supplementary note with explicit explanation of why they differ.

**Coverage of GATC 6mA (0.27% occupancy).** GATC is a Dam methylase target. Its 0.27% occupancy in M145 is mentioned in passing but never interpreted. Does M145 have a Dam methylase? Is there published evidence? The near-zero GATC occupancy is potentially informative (suggesting M145 lacks Dam) and deserves one sentence of interpretation.

---

## Strengths

**Methodological innovation: Simpson's paradox characterization in bacterial epigenomics.** The identification of three independent Simpson's paradox instances — driven by the linear chromosome's core/arm architecture — is the paper's most significant and generalizable contribution. The methodological recommendation (CMH stratification as standard practice for linear-chromosome organisms) is important and publication-worthy on its own. This analysis is rigorous, reproducible, and directly addresses a gap in the field.

**Transparent negative results reporting.** The "Negative results define the scope of the model" section, and the explicit discussion of the unidentified proximate mechanism for the 57 exposed regulators, represent exemplary scientific transparency. The paper does not oversell its findings and explicitly delineates what is observational versus causal. This is increasingly rare and commendable.

**Quantitative rigor of the Shielded/Exposed dichotomy.** The ROC analysis with bootstrap CI, the Youden's J threshold selection, the multiple negative-control AUC comparisons (expression level AUC = 0.547; evolutionary conservation AUC = 0.459), and the CMH-adjusted ORs for the protection zone all reflect careful statistical practice. The 293 bp threshold and AUC = 0.917 are well-supported by the data presented.

**Multidimensional validation of the protection zone.** The two-tier mechanism (sequence counter-selection + protein occupancy) is examined from multiple angles — expression-independence, TF binding site analysis, CMH stratification, evolutionary conservation — each providing an independent line of evidence against simpler alternative explanations. The finding that individual TF binding sites are not protective (fold = 1.157 enrichment rather than depletion) is a counterintuitive but important negative result that strengthens the collective-occupancy model.

**Internal control design.** Using AAGCCCG 6mA (minimal geographic bias) as an internal control for the geographic-confounding argument, within the same dataset, is an elegant design feature that strengthens the Simpson's paradox demonstration.

**Nine-sample Nanopore + RNA-seq design.** Three biological replicates per timepoint, matched methylome and transcriptome, is a substantial experimental investment for a prokaryotic epigenomics study. The overall sequencing quality (98.54% alignment rate for RNA-seq; 27–52x coverage for Nanopore on 8 of 9 samples) is appropriate.

**Rigorous evolutionary controls for the Shielded/Exposed dichotomy.** The multi-metric evolutionary analysis (SCO locus tag rate, GC3, codon usage, named product rate, composite AUC = 0.459) effectively rules out the alternative hypothesis that exposed regulators are evolutionarily recent acquisitions. This pre-empts an obvious reviewer objection.

---

## Recommended Journal Determination

### Primary recommendation: **EMBO Journal**

**Rationale:** EMBO Journal has a strong track record of publishing genome-scale bacterial regulatory studies with quantitative computational frameworks, including papers that present important negative results alongside structural models. The paper's scope — a single organism, three timepoints, correlational epigenomics — fits EMBO Journal's appetite for systematic genomic analyses that advance conceptual frameworks without requiring complete mechanistic dissection. The methodological contribution (Simpson's paradox in bacterial epigenomics) and the Shielded/Exposed architecture are sufficiently novel for EMBO Journal's scope. The journal's word limits and figure allowances also accommodate the paper's current structure better than Nature Microbiology.

**Conditions for EMBO Journal acceptance:** Resolution of all Major Concerns listed above (especially MC1, MC4, MC7 as immediate editorial issues; MC2 and MC6 as scientific issues that need text revision), plus addressing the minor concern about figure number placeholders. EMBO Journal will also expect the coverage/TSS-sensitivity analyses to be prominently presented. A revised abstract framing the Simpson's paradox contribution as the primary finding (with the structural model as the secondary finding) would better match EMBO Journal's editorial priorities.

### Secondary recommendation (with additional experiments): **Nature Microbiology**

**Rationale:** Nature Microbiology is achievable if the following are added: (a) at minimum a methyltransferase knockout (preferably *SC_RS17645*) demonstrating that loss of AAGCCCG methylation preferentially affects exposed regulators' expression; (b) ChIP-seq for one or two NAPs (ideally HupA and HupS) demonstrating differential occupancy at shielded vs exposed TSSs. These two experiments would transform the paper from a correlational epigenomics study into a mechanistic epigenomics paper, which is what Nature Microbiology requires. The Simpson's paradox analysis and the structural model are genuinely Nature Microbiology-level conceptual contributions; the missing element is experimental causal validation.

**Not recommended: Nucleic Acids Research.** NAR is appropriate for database papers and methodological contributions. This manuscript's central contribution is a biological model of *S. coelicolor* epigenomics with a methodological component, not a database or purely computational tool. NAR would undervalue the biological novelty of the Shielded/Exposed framework.

---

## Priority Action List for Revision

The following items must be addressed before any journal submission, ranked by urgency:

1. **[CRITICAL — MC7]** Remove all `<!-- E14: ... -->` comment tags from the manuscript. Resolve all "Supplementary Figure XX" placeholders with definitive figure numbers. Fix the duplicate "Fig. 3e" panel reference. Conduct a full cross-reference audit (main figures vs supplementary).

2. **[CRITICAL — MC1]** Correct 5.9% → 5.4% in the Introduction. Audit all occurrences of the exposed-regulator percentage across the full manuscript.

3. **[CRITICAL — MC4]** Reconcile the Fisher OR discrepancy (138,440 vs 145,448). Run a single authoritative analysis, report the correct number with CI, and ensure the 921× fold-enrichment and Fisher OR are clearly distinguished in all figures and text.

4. **[HIGH — MC2]** Reframe the "two methylation systems" claim to acknowledge the majority-unassigned 6mA fraction as evidence of a third (uncharacterized) system. Assess geographic bias of unassigned 6mA sites to validate the internal control logic.

5. **[HIGH — MC3]** Clarify the Type I vs solo MTase classification of SC_RS17645. Present genomic neighborhood data (presence/absence of HsdR/HsdS). Replace "solo MTase" with "candidate MTase" in all non-Methods sections. Optionally, add a bioinformatic analysis of the M145 RM repertoire from REBASE.

6. **[HIGH — MC5]** Report AUC = 0.917 computed separately at T1, T2, and T3 (or at minimum T1-only vs all-timepoints pooled) to demonstrate robustness of the protection zone classification across the developmental time course.

7. **[MEDIUM — MC6]** Revise causal language in the Abstract and Introduction to consistently use correlational framing. Replace "DNA methylation does not directly control transcription" with "no evidence for direct methylation control of transcription." Move the dRNA-seq TSS sensitivity analysis (AUC = 0.923) to the main text as a primary result.

8. **[MEDIUM]** Apply multiple-testing correction to TF family enrichment analyses (TetR p = 0.028). Report corrected p-values or explicitly state that no correction was applied and why.

9. **[MEDIUM]** Remove or substantially condense the TCS 7/7 asymmetry discussion given p = 0.83 (not significant). If retained, restrict to factual description only.

10. **[LOW]** Standardize "Gatekeeper" model nomenclature across all sections or remove it entirely. Report which sample has 14.3x Nanopore coverage and demonstrate that excluding it does not alter key results.

---

*Review conducted by: Claude (independent automated reviewer), 2026-05-11. Based on full reading of: Abstract, Introduction, Results, Discussion, Methods (full text), and figure_claims_map.pdf.*
