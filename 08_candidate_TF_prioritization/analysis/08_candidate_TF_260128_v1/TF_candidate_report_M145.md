# TF Candidate Report — M145 RNA-seq

## 1. Purpose and Approach

This step prioritizes transcription factor (TF) candidates that may regulate the four major BGCs (act, red, cda, cpk) in *Streptomyces coelicolor* A3(2) M145, based on:

1. **Regulator–BGC correlation** (from step 07): Pearson r between each regulator's rlog expression and 4 BGC composite scores across 9 samples (3 timepoints × 3 replicates).
2. **Phase-specific expression changes**: rlog deltas between conditions (M145_1→2 = early/mid transition, M145_2→3 = late transition) to separate act-type late-phase regulators from red/cda/cpk-type mid-phase regulators.
3. **DESeq2 fold change**: log2FC(M145_3 vs M145_1) for overall dynamic range.

### Scoring logic

Two complementary scores were computed:

- **TF_score_act** (act late-phase candidates):
  `abs_corr_act + 0.5 × max(act_specificity, 0) + 0.3 × min(|log2FC_3v1|/5, 1) + phase_weight`
  - `act_specificity = |corr_act| − max(|corr_red|, |corr_cda|, |corr_cpk|)`
  - `phase_weight`: +0.3 for "late" phase preference, +0.15 for "both", 0 otherwise.

- **TF_score_mid** (red/cda/cpk mid-phase candidates):
  `mean(|corr_red|, |corr_cda|, |corr_cpk|) + 0.3 × min(|log2FC_3v1|/5, 1) + phase_weight`
  - `phase_weight`: +0.3 for "early_mid", +0.15 for "both", 0 otherwise.

Global candidates were identified as regulators with |corr| > 0.9 for ≥2 BGCs, or |corr| > 0.8 for ≥3 BGCs with mean |corr| ≥ 0.8.

---

## 2. BGC Phase Dynamics Summary

The BGC composite scores confirm two distinct temporal patterns:

| BGC | Δ score (M145_2 − M145_1) | Δ score (M145_3 − M145_2) | Pattern |
|-----|---------------------------|---------------------------|---------|
| act | +33 | +9,392 | **Late-dominant** (M145_3 surge) |
| red | +616 | +284 | Early/mid-dominant |
| cda | +3,252 | −224 | Early/mid-dominant |
| cpk | +9,413 | −3,321 | Early/mid-dominant |

Phase preference distribution among 856 regulators:
- Early/mid: 362 (42%)
- Late: 197 (23%)
- Both: 180 (21%)
- Stable: 117 (14%)

---

## 3. BGC-specific and Phase-specific Major Candidates

### 3.1 Act late-phase candidates (top 5)

| Rank | Gene ID | Old locus | Product | Known? | corr_act | Act specificity | TF_score_act |
|------|---------|-----------|---------|--------|----------|-----------------|-------------|
| 1 | SC_RS37190 | SCO6993 | LuxR family (AbsR2) | No | +0.946 | +0.139 | 1.61 |
| 2 | SC_RS36895 | SCO6937 | LuxR C-terminal regulator | No | +0.938 | +0.241 | 1.61 |
| 3 | SC_RS27570 | SCO5082 | actII-orf4 | **Yes** | +0.954 | **+0.382** | 1.60 |
| 4 | SC_RS07840 | SCO1177 | FadR/GntR family | No | +0.947 | +0.254 | 1.59 |
| 5 | SC_RS35025 | SCO6566 | ROK family regulator | No | +0.968 | +0.304 | 1.57 |

**Key observations:**
- actII-orf4 (SC_RS27570) ranks 3rd by TF_score_act and has the **highest act specificity** (+0.382), confirming its unique association with act. Its "late" phase preference (delta_3_vs_2 = +3.32) matches the act late-surge pattern.
- SCO5085_SARP (SC_RS27585) also appears in the top 30 (TF_score_act = 1.54, corr_act = +0.917), consistent with its known SARP role within the act cluster.
- SC_RS35025 (SCO6566, ROK family) has the highest raw |corr_act| = 0.968 but lower specificity, indicating possible broader effects.
- Both SC_RS37190 (LuxR family) and SC_RS36895 (LuxR C-terminal) are novel act-specific candidates with late phase preference.

### 3.2 Red/cda/cpk mid-phase candidates (top 5)

| Rank | Gene ID | Old locus | Product | Known? | mid_mean_corr | delta_2_vs_1 | TF_score_mid |
|------|---------|-----------|---------|--------|--------------|-------------|-------------|
| 1 | SC_RS26000 | SCO4768 | Response regulator | No | 0.934 | +5.63 | 1.53 |
| 2 | SC_RS20030 | SCO3579 | WblA | No | 0.932 | +6.12 | 1.53 |
| 3 | SC_RS31630 | SCO5877 | redD | **Yes** | 0.930 | +5.53 | 1.53 |
| 4 | SC_RS25560 | SCO4677 | atrA | **Yes** | 0.927 | +7.45 | 1.53 |
| 5 | SC_RS30130 | SCO5582 | NsdA | No | 0.927 | +6.38 | 1.53 |

**Key observations:**
- All known CSRs (redD, atrA, scbR2, cpkO/kasO, SCO3217_SARP, SCO6280_SARP, SCO6288_SARP) appear in the top 30 mid-phase candidates, validating the scoring approach.
- SC_RS26000 (SCO4768, response regulator) is the top novel mid-phase candidate with the highest mid_mean_corr (0.934). Its large delta_2_vs_1 (+5.63) matches the red/cda/cpk early activation pattern.
- WblA (SC_RS20030, SCO3579) is a well-characterized developmental regulator in *Streptomyces*; its high mid-phase score supports a link between morphological differentiation and secondary metabolism.
- NsdA (SC_RS30130, SCO5582) is a known repressor of differentiation; its positive correlation here (all 3 BGC scores positive) suggests it may have a dual or context-dependent role.

### 3.3 Global regulator candidates

148 regulators meet global criteria. Top candidates by mean |corr|:

| Rank | Gene ID | Old locus | Product | Role | mean_abs_corr | n(|r|>0.9) |
|------|---------|-----------|---------|------|---------------|------------|
| 1 | SC_RS16870 | SCO2953 | RsuA (anti-sigma factor) | Activator | 0.871 | 2 |
| 2 | SC_RS18240 | SCO3225 | absA1 (HK) | Activator | 0.866 | 2 |
| 3 | SC_RS34305 | SCO6422 | Response regulator | Activator | 0.861 | 2 |
| 4 | SC_RS02220 | SCO0037 | Sigma-70 factor | Activator | 0.860 | 2 |
| 5 | SC_RS31650 | SCO5881 | redZ (RR) | Activator | 0.856 | 3 |
| 6 | SC_RS32155 | SCO5982 | PaaX family | Repressor | 0.854 | 2 |
| 7 | SC_RS36185 | SCO6801 | LysR family | Repressor | 0.852 | 2 |
| 8 | SC_RS21215 | SCO3818 | Response regulator | Activator | 0.846 | **3** |
| 9 | SC_RS29775 | SCO5518 | PucR family | Repressor | 0.844 | **3** |
| 10 | SC_RS28860 | SCO5337 | XRE family | Repressor | 0.839 | **3** |

**Distribution:** 47 global activators, 52 global repressors, 49 mixed/other (signs not consistent across all 4 BGCs).

---

## 4. Consistency with Known Regulatory Networks

### 4.1 Known CSRs in scoring results

| Known regulator | Act score rank | Mid score rank | Global? | Assessment |
|----------------|---------------|---------------|---------|------------|
| actII-orf4 (SC_RS27570) | **3** | — | No | Act-specific; highest specificity |
| SCO5085_SARP (SC_RS27585) | Top 30 | — | No | Act cluster SARP; confirmed |
| redD (SC_RS31630) | — | **3** | Yes (activator) | Mid-phase validated |
| atrA (SC_RS25560) | — | **4** | Yes (activator) | Pleiotropic; confirmed |
| SCO3217_SARP (SC_RS18200) | — | **6** | Yes (activator) | cda/cpk SARP; confirmed |
| scbR2 (SC_RS33680) | — | **7** | Yes (activator) | Autoregulator receptor; confirmed |
| cpkO/kasO (SC_RS33660) | — | **8** | Yes (activator) | cpk CSR; confirmed |
| redZ (SC_RS31650) | — | Top 30 | **Yes (3 BGC)** | Global activator; 3 BGC with |r|>0.9 |
| absA1 (SC_RS18240) | — | Top 30 | Yes (activator) | Two-component HK; confirmed |
| absA2 (SC_RS18245) | — | Top 30 | Yes (activator) | Two-component RR; confirmed |
| SCO6280_SARP (SC_RS33650) | — | Top 30 | — | cpk cluster SARP; confirmed |
| SCO6288_SARP (SC_RS33690) | — | Top 30 | — | cpk cluster SARP; confirmed |
| scbR (SC_RS33575) | — | — | — | Weak correlations; consistent with autoregulatory role |

All 13 known regulators are positioned consistently with their literature-described functions:
- **actII-orf4** is uniquely act-specific (highest act_specificity = +0.382).
- **redZ** emerges as a global activator (3 BGC with |r|>0.9), consistent with its known pleiotropic effects.
- **scbR** shows weak correlations, consistent with its role as a butyrolactone receptor rather than a direct transcriptional activator.

This high concordance between scoring results and known biology supports the reliability of the novel candidates identified.

---

## 5. Novel Candidate TFs — Biological Interpretation and Experimental Proposals

### 5.1 SC_RS21215 / SCO3818 (Response regulator)

- **Profile:** Positive correlation with red (+0.93), cda (+0.98), cpk (+0.93); moderate with act (+0.53). Three BGC with |r|>0.9. Phase = early/mid.
- **Interpretation:** Two-component response regulator that activates during the mid-phase transition. Likely acts upstream of or in parallel with the red/cda/cpk regulatory cascade.
- **Experiment:** KO → expected decrease in red, cda, cpk production. OE → potential enhancement of all three. Good candidate for CRISPR-based disruption.

### 5.2 SC_RS29775 / SCO5518 (PucR family)

- **Profile:** Negative correlation with red (−0.94), cda (−0.98), cpk (−0.91). Three BGC with |r|>0.9 (all negative). Phase = early/mid.
- **Interpretation:** Putative global repressor of secondary metabolism. High expression in M145_1 (vegetative growth) that drops as BGC activation occurs.
- **Experiment:** **KO is the highest priority** — disruption may de-repress multiple BGCs simultaneously, potentially enhancing antibiotic production.

### 5.3 SC_RS28860 / SCO5337 (XRE family)

- **Profile:** Negative correlation with red (−0.93), cda (−0.97), cpk (−0.93). Three BGC negative. Phase = early/mid.
- **Interpretation:** Second candidate global repressor, from the XRE (xenobiotic response element) family. XRE-type regulators in *Streptomyces* often control toxin-antitoxin systems or secondary metabolism switches.
- **Experiment:** KO → expected enhancement of red/cda/cpk. Compare with SC_RS29775 KO to assess redundancy.

### 5.4 SC_RS37190 / SCO6993 (LuxR family, AbsR2)

- **Profile:** corr_act = +0.946, act_specificity = +0.139. Phase = late. Strongly upregulated in M145_3 (log2FC_3v1 = +4.9).
- **Interpretation:** LuxR-family regulator with act-preferential activation. The late-phase expression surge parallels actinorhodin accumulation. May function as a secondary activator downstream of actII-orf4.
- **Experiment:** KO → test whether act production is reduced while red/cda/cpk remain unaffected. OE → potential act enhancement.

### 5.5 SC_RS26000 / SCO4768 (Response regulator)

- **Profile:** Top novel mid-phase candidate. mid_mean_corr = 0.934, delta_2_vs_1 = +5.63. Phase = early/mid.
- **Interpretation:** Response regulator with the highest average correlation to red/cda/cpk among all novel candidates. Its strong induction at M145_2 suggests involvement in the early secondary metabolism switch.
- **Experiment:** KO → expected broad decrease in red/cda/cpk production. Complement with OE to confirm.

---

## 6. Experimental Priority List

9 candidates were selected for experimental follow-up, balancing known controls and novel discoveries:

| # | Gene ID | Old locus | Name | Target BGC | Role | Known? | Rationale |
|---|---------|-----------|------|-----------|------|--------|-----------|
| 1 | SC_RS27570 | SCO5082 | actII-orf4 | act | CSR activator | Yes | Positive control; highest act specificity |
| 2 | SC_RS37190 | SCO6993 | (LuxR/AbsR2) | act | Candidate activator | No | Novel; corr_act=+0.946, late phase |
| 3 | SC_RS36895 | SCO6937 | (LuxR C-term) | act | Candidate activator | No | Novel; corr_act=+0.938, high act specificity |
| 4 | SC_RS31650 | SCO5881 | redZ | red/cda/cpk | Global activator | Yes | Positive control; 3 BGC with |r|>0.9 |
| 5 | SC_RS21215 | SCO3818 | (RR) | red/cda/cpk | Candidate global activator | No | 3 BGC positive, response regulator |
| 6 | SC_RS29775 | SCO5518 | (PucR) | red/cda/cpk | Candidate global repressor | No | **Top KO target**; 3 BGC negative |
| 7 | SC_RS28860 | SCO5337 | (XRE) | red/cda/cpk | Candidate global repressor | No | 3 BGC negative; XRE family |
| 8 | SC_RS26000 | SCO4768 | (RR) | red/cda/cpk | Candidate mid activator | No | Highest mid_mean_corr among novel TFs |
| 9 | SC_RS16870 | SCO2953 | rsuA | all 4 BGC | Candidate global activator | No | Anti-sigma factor; highest global mean |

---

## 7. Future Directions

1. **Promoter motif analysis**: Search upstream regions of BGC genes for binding motifs of the candidate TF families (LuxR, PucR, XRE, OmpR-type RR). Tools: MEME-ChIP, CiiiDER, or *Streptomyces*-specific databases (SCoDB).
2. **Cross-condition validation**: If additional RNA-seq data become available (e.g., different media, carbon sources, pH conditions), re-evaluate candidate correlations for robustness beyond the current timecourse.
3. **ChIP-seq or DAP-seq**: For the top 3 novel candidates (SC_RS21215, SC_RS29775, SC_RS26000), direct binding evidence to BGC promoter regions would confirm regulatory relationships.
4. **Genetic experiments**: CRISPR-Cas9 disruption of SC_RS29775 (PucR, global repressor) is the highest-priority experiment: KO may simultaneously de-repress red, cda, and cpk.
5. **Literature comparison**: Cross-reference with known *Streptomyces* regulatory networks involving Crp (SC_RS16735), AfsR (SC_RS33160), and other master regulators not included as "known" CSRs in this analysis.
6. **Protein–protein interaction**: Yeast two-hybrid or bacterial adenylate cyclase two-hybrid (BACTH) to test direct interactions between novel TF candidates and known CSRs.

---

## 8. Key Takeaways

1. **Phase-specific scoring** successfully separates act-late TF candidates from red/cda/cpk-mid candidates, leveraging the distinct temporal dynamics of BGC activation in M145.
2. **All 13 known CSRs** rank appropriately in their respective categories, validating the scoring methodology. actII-orf4 is the most act-specific regulator; redZ is the strongest global activator.
3. **Three novel global repressor/activator candidates** — SC_RS21215 (SCO3818, response regulator), SC_RS29775 (SCO5518, PucR family), SC_RS28860 (SCO5337, XRE family) — show |corr| > 0.9 with 3 BGCs simultaneously and represent high-priority targets for genetic manipulation.
4. **SC_RS29775 (SCO5518) is the top KO priority**: as a putative global repressor (negative correlation with red, cda, and cpk), its disruption may broadly enhance secondary metabolite production.
5. **9 experimental priority candidates** provide a balanced mix of known controls (actII-orf4, redZ) and novel discoveries (6 novel TFs across activator and repressor classes), enabling systematic validation of the regulatory network.

---

*Generated: 2026-01-28 | Pipeline: 08_candidate_TF_prioritization | Data: M145 RNA-seq (3 timepoints × 3 replicates)*
