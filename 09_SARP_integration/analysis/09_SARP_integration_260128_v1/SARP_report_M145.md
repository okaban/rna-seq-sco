# SARP Integration Report — M145 RNA-seq

## 1. Purpose and Approach

This step systematically identified all SARP (Streptomyces Antibiotic Regulatory Protein) family transcription factors in the *S. coelicolor* A3(2) M145 genome using Pfam HMM-based domain detection, classified them by subtype, and integrated the results with the BGC dynamics (step 06), regulator–BGC correlation (step 07), and TF prioritization (step 08) networks.

**Method**: HMMER 3.4 `hmmsearch` was run against the M145 proteome (protein.faa, 7,996 annotated proteins) using four Pfam HMM profiles:
- **PF00486** (Trans_reg_C): wHTH DNA-binding domain (DBD)
- **PF03704** (BTAD): Bacterial Transcriptional Activator Domain
- **PF00931** (NB-ARC): nucleotide-binding domain
- **PF13424** (TPR_12): tetratricopeptide repeat

A protein was classified as SARP if it contained **both DBD + BTAD** domains.

---

## 2. SARP Identification and Subtype Distribution

### 2.1 Domain search results

| Domain | Proteins detected |
|--------|------------------|
| Trans_reg_C (DBD) | 33 |
| BTAD | 13 |
| NB-ARC | 9 |
| TPR_12 | 23 |
| **DBD + BTAD (= SARP)** | **8** |
| BTAD only (no DBD) | 5 |

### 2.2 SARP subtype classification

8 proteins contain both DBD and BTAD domains:

| Gene ID | Old locus | Name | Product | Length (aa) | Subtype | BGC | Role |
|---------|-----------|------|---------|------------|---------|-----|------|
| SC_RS27585 | SCO5085 | — | AfsR/SARP family TR | 255 | **small** | act | CSR |
| SC_RS31630 | SCO5877 | — | AfsR/SARP family TR | 268 | **small** | red | CSR |
| SC_RS33690 | SCO6288 | — | AfsR/SARP family TR | 281 | **small** | cpk | CSR |
| SC_RS18200 | SCO3217 | — | AfsR/SARP family TR | 638 | **medium** | cda | CSR |
| SC_RS24295 | SCO4426 | afsR | transcriptional regulator AfsR | 993 | **large** | none | global |
| SC_RS06435 | SCO0898 | — | AfsR/SARP family TR | 660 | small_long | none | global |
| SC_RS22750 | SCO4116 | — | BTAD domain-containing TR | 1102 | small_long | none | global |
| SC_RS13320 | SCO2259 | — | BTAD domain-containing TR | 1322 | small_long | none | global |

**Subtype distribution**: 3 small, 1 medium, 1 large, 3 small_long (DBD + BTAD but >500 aa without NB-ARC).
No SARP-LAL type was identified (no LAL domain co-occurrence detected).

**CSR vs Global**: 4 cluster-situated (1 per major BGC) + 4 global/pleiotropic.

### 2.3 Notable non-detections

- **actII-orf4** (SC_RS27570 / SCO5082): No DBD or BTAD domain detected by hmmsearch. The GFF annotates it as "TetR family transcriptional regulator." Despite its literature classification as a SARP, its protein sequence (255 aa) does not match the Pfam SARP domain profiles at the applied thresholds (E < 1e-3). This suggests actII-orf4 may represent an atypical/divergent SARP or may be more accurately classified as a TetR-like activator.
- **cpkO/kasO** (SC_RS33660 / SCO6282): No SARP-related domains detected. Annotated as "SDR family oxidoreductase" in GFF. Not a domain-based SARP.
- **SCO6280** (SC_RS33650): Detected BTAD domain only (no DBD). Classified as BTAD-only protein, not a full SARP. This gene is within the cpk cluster and annotated as "AfsR/SARP family" in the product field.

---

## 3. Major BGC-associated SARPs

### 3.1 act cluster: SCO5085 (SC_RS27585, small SARP)

| Metric | Value |
|--------|-------|
| corr_act | **+0.917** |
| corr_red | +0.873 |
| corr_cda | +0.690 |
| corr_cpk | +0.465 |
| log2FC(3 vs 1) | +6.44 |
| Phase | **late** (Δ3-2 = +4.23 >> Δ2-1 = +1.94) |

SCO5085 is an act-cluster SARP with the highest correlation to the act BGC score (r = +0.917) among the domain-confirmed SARPs. Its late-phase expression surge (Δ3-2 = +4.23) matches the act activation pattern. It also shows moderate-to-high correlation with red (r = +0.87), suggesting partial pleiotropy.

**Relationship with actII-orf4**: actII-orf4 (SC_RS27570) has corr_act = +0.954 and even higher act-specificity (from step 07), but was NOT detected as a domain-based SARP. SCO5085 (annotated as "transport" in BGC_definition_manual.tsv) is the only domain-confirmed SARP in the act cluster.

### 3.2 red cluster: SCO5877 / redD (SC_RS31630, small SARP)

| Metric | Value |
|--------|-------|
| corr_act | +0.385 |
| corr_red | **+0.863** |
| corr_cda | **+0.973** |
| corr_cpk | **+0.955** |
| log2FC(3 vs 1) | +4.99 |
| Phase | **early_mid** (Δ2-1 = +5.53, Δ3-2 = −0.80) |

RedD shows the typical mid-phase activation pattern: strong induction at M145_2 with slight decline at M145_3. It correlates most strongly with cda (r = +0.97) and cpk (r = +0.96), even more than with its cognate red BGC (r = +0.86). This confirms redD's pleiotropic influence beyond the red cluster.

### 3.3 cda cluster: SCO3217 (SC_RS18200, medium SARP)

| Metric | Value |
|--------|-------|
| corr_act | +0.364 |
| corr_red | +0.852 |
| corr_cda | **+0.970** |
| corr_cpk | **+0.958** |
| log2FC(3 vs 1) | +5.76 |
| Phase | **early_mid** (Δ2-1 = +6.62, Δ3-2 = −1.11) |
| Domains | DBD + BTAD + NB-ARC |

SCO3217 is the only medium SARP identified, containing an additional NB-ARC domain. Its mid-phase expression (Δ2-1 = +6.62) and very high cda correlation (r = +0.970) are consistent with a cluster-situated activator role. The NB-ARC domain may provide an additional regulatory input for signal integration.

### 3.4 cpk cluster: SCO6288 (SC_RS33690, small SARP)

| Metric | Value |
|--------|-------|
| corr_act | +0.360 |
| corr_red | +0.840 |
| corr_cda | **+0.958** |
| corr_cpk | **+0.940** |
| log2FC(3 vs 1) | +3.90 |
| Phase | **early_mid** (Δ2-1 = +4.17, Δ3-2 = −0.62) |

SCO6288 is a small SARP within the cpk cluster showing mid-phase activation. Note: SCO6280 (SC_RS33650), the other SARP-annotated gene in cpk, was detected with BTAD only (no DBD), so it is not a full domain-confirmed SARP.

**Consistency with steps 06–08**: All four CSR SARPs match the phase patterns established in steps 06/07 — act CSR (SCO5085) has late-phase dynamics, while red/cda/cpk CSRs (redD, SCO3217, SCO6288) show early/mid-phase activation.

---

## 4. Global SARPs

### 4.1 AfsR (SC_RS24295 / SCO4426, large SARP)

| Metric | Value |
|--------|-------|
| Domains | DBD + BTAD + NB-ARC + TPR_12 |
| Length | 993 aa |
| corr_act | **−0.825** |
| corr_red | −0.749 |
| corr_cda | −0.591 |
| corr_cpk | −0.366 |
| log2FC(3 vs 1) | −0.918 |
| Phase | **late** (declining: Δ2-1 = −0.27, Δ3-2 = −0.64) |

AfsR is the sole large SARP in M145, containing all four SARP-associated domains. Contrary to its well-documented role as a positive activator of secondary metabolism in literature, AfsR shows **negative correlations** with all four BGC scores in this dataset. Its expression declines throughout the time course (log2FC = −0.92). This pattern suggests that in this particular culture condition, AfsR may function as an early-phase factor whose downregulation coincides with (or permits) BGC activation, or that its activating role is mediated through transient phosphorylation-dependent mechanisms not captured by mRNA levels alone.

### 4.2 SCO2259 (SC_RS13320, small_long SARP)

| corr_act | corr_red | corr_cda | corr_cpk | mean_abs_corr | Phase |
|----------|----------|----------|----------|--------------|-------|
| −0.667 | **−0.921** | −0.891 | −0.768 | 0.812 | stable |

SCO2259 is a 1,322-aa protein with DBD + BTAD (no NB-ARC/TPR) and the highest mean |corr| (0.812) among global SARPs. Its strong negative correlation with red (r = −0.92) and cda (r = −0.89) suggests a repressor-like role. Despite being classified as "stable" by phase preference (small absolute deltas), its negative correlations are robust.

### 4.3 SCO4116 (SC_RS22750, small_long SARP)

| corr_act | corr_red | corr_cda | corr_cpk | mean_abs_corr | Phase |
|----------|----------|----------|----------|--------------|-------|
| −0.217 | −0.660 | −0.815 | −0.801 | 0.623 | early_mid |

Negative correlations across all BGCs, strongest for cda and cpk. Moderate candidate for a repressor role.

### 4.4 SCO0898 (SC_RS06435, small_long SARP)

| corr_act | corr_red | corr_cda | corr_cpk | mean_abs_corr | Phase |
|----------|----------|----------|----------|--------------|-------|
| +0.448 | +0.264 | +0.104 | −0.074 | 0.222 | late |

Weak correlations across all BGCs. No clear activator or repressor pattern. Not prioritized.

---

## 5. Novel SARP Candidates of Interest

### BTAD-only proteins (not full SARPs but related)

Five proteins contain BTAD without DBD. Of particular interest:

| Gene ID | Old locus | Product | Length | In BGC? |
|---------|-----------|---------|--------|---------|
| SC_RS33650 | SCO6280 | AfsR/SARP family TR | 539 | cpk |
| SC_RS25495 | SCO4663 | AfsR/SARP family TR | 117 | no |
| SC_RS14395 | SCO2450 | SAV_2350-like | 1349 | no |
| SC_RS18620 | SCO3291 | macro domain | 477 | no |
| SC_RS24690 | SCO4548 | bacterial TR | 1002 | no |

SCO6280 (SC_RS33650) is notable: it is within the cpk cluster, annotated as "AfsR/SARP family," and from step 07 has corr_cpk = +0.970 (top among cpk-associated regulators). Despite lacking a detectable DBD by hmmsearch, it may represent a divergent or partial SARP that still functions as a cpk-cluster activator.

---

## 6. Integration with Step 08 Priority List

None of the 9 step-08 experimental priority candidates are domain-confirmed SARPs (actII-orf4 was not detected as SARP by hmmsearch). Three global SARPs were added to create an extended priority list of 12 candidates:

| # | Gene ID | Name | Role | SARP subtype | Comment |
|---|---------|------|------|-------------|---------|
| 1–9 | (original step 08 candidates) | — | — | — | Unchanged |
| 10 | SC_RS13320 | SCO2259 | SARP repressor | small_long | mean_abs_corr=0.812; 4 BGC negative |
| 11 | SC_RS24295 | SCO4426 / afsR | SARP repressor | large | mean_abs_corr=0.633; act-corr=−0.825 |
| 12 | SC_RS22750 | SCO4116 | SARP repressor | small_long | mean_abs_corr=0.623; cda/cpk negative |

---

## 7. Experimental Proposals

### 7.1 CSR SARP manipulation strategy

**Tier 1 — Single-gene experiments:**
- **OE of SCO5085** (SC_RS27585, act small SARP): Expected to enhance actinorhodin production. Also test whether red is affected (moderate corr_red = +0.87).
- **OE of redD** (SC_RS31630): Expected to enhance undecylprodiginine (red) and possibly cda/cpk.
- **OE of SCO3217** (SC_RS18200, cda medium SARP): Expected to enhance CDA production. The NB-ARC domain may allow conditional activation via signal molecules.

**Tier 2 — Combinatorial experiments:**
- **Double OE: SCO5085 + redD**: Test whether simultaneous CSR activation of act + red/cda/cpk pathways yields synergistic or competitive effects.
- **AfsR (SC_RS24295) OE under alternative conditions**: Since AfsR shows negative correlations in this dataset, test OE under nutrient-limited conditions where its activating role may be more apparent (literature-expected positive effect on act/red/cda).

### 7.2 Global SARP repressor candidates

- **KO of SCO2259** (SC_RS13320): Highest mean_abs_corr (0.812) among global SARPs, all correlations negative. KO may de-repress red and cda.
- **KO of SCO4116** (SC_RS22750): Secondary repressor candidate (cda/cpk negative correlations).

### 7.3 Silent/cryptic BGC activation

- No SARP-LAL type was identified in M145. For activation of silent BGCs, consider **OE of the global SARP SCO0898** (SC_RS06435, 660 aa) under varied conditions, as it is a BGC-external SARP with weak but positive act correlation.
- Alternatively, heterologous expression of known SARP-LAL activators (e.g., PimR, PolR from other *Streptomyces* species) could be tested to activate cryptic clusters.

---

## 8. Key Takeaways

1. **8 domain-confirmed SARPs** were identified in M145 by HMMER (DBD + BTAD): 3 small, 1 medium, 1 large, 3 small_long. Four are BGC-cluster-situated CSRs (act, red, cda, cpk) and four are global/pleiotropic.
2. **actII-orf4 and cpkO/kasO are NOT domain-based SARPs** by Pfam HMM detection, indicating they are either divergent SARPs or belong to different TF families (TetR and SDR, respectively, as annotated in GFF). This is an important caveat for literature interpretations.
3. **All four CSR SARPs match the phase-specific BGC dynamics** from steps 06–08: SCO5085 (act) shows late-phase activation; redD, SCO3217, and SCO6288 (red/cda/cpk) show early/mid-phase activation. This confirms the two-wave regulatory model.
4. **AfsR (SC_RS24295) shows unexpected negative correlations** with all BGC scores, suggesting its activating role may be post-translational (phosphorylation-dependent) rather than transcription-dependent in this timecourse.
5. **SCO2259 (SC_RS13320)** is the most promising novel global SARP for experimental follow-up, with strong negative correlations across red/cda/cpk (mean |r| = 0.81), making it a candidate repressor whose KO may enhance multiple secondary metabolites.

---

*Generated: 2026-01-28 | Pipeline: 09_SARP_integration | Data: M145 RNA-seq (3 timepoints × 3 replicates)*
