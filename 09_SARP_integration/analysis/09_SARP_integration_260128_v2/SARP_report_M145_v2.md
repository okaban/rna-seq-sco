# SARP Integration Report v2 — M145 RNA-seq

## 1. Purpose and Approach (v2)

### 1.1 v1 recap

In v1, all SARP family transcription factors in *S. coelicolor* A3(2) M145 were identified by Pfam HMM-based domain detection (HMMER 3.4). A protein was classified as SARP if it contained **both DBD (PF00486/Trans_reg_C) and BTAD (PF03704)** domains. This yielded **8 domain-confirmed (structural) SARPs**.

However, **actII-orf4** (SC_RS27570 / SCO5082) — the canonical pathway-specific activator of the actinorhodin (act) BGC and the historical "founder" of the SARP family — was **not detected** by hmmsearch. Its GFF annotation is "TetR family transcriptional regulator ActII," and it lacks detectable Trans_reg_C or BTAD domains at the applied Pfam thresholds.

### 1.2 v2 objective

v2 introduces a **dual-layer SARP definition** to resolve this gap:

1. **Structural SARP** (`is_SARP_structural`): proteins with both DBD + BTAD domains detected by hmmsearch (v1 definition, 8 proteins).
2. **Functional SARP** (`is_SARP_functional`): proteins defined as SARP or ActII-ORF4-like CSRs in the primary literature and review articles, added via a curated whitelist.

The union produces an **extended SARP set of 9 proteins**, with actII-orf4 added as the sole `functional_only` member.

---

## 2. SARP Definition: Dual-Layer Classification

### 2.1 Classification flags

| Flag | Definition |
|------|-----------|
| `is_SARP_structural` | DBD + BTAD detected by hmmsearch (v1 criterion) |
| `is_SARP_functional` | Curated from literature as SARP / ActII-ORF4-like CSR |
| `SARP_category` | `both`, `structural_only`, or `functional_only` |

### 2.2 Category distribution (v2)

| SARP_category | Count | Members |
|---------------|-------|---------|
| **both** | 3 | SCO5877/redD (red CSR), SCO3217 (cda CSR), SCO6288 (cpk CSR) |
| **structural_only** | 5 | SCO5085 (act), SCO4426/AfsR, SCO0898, SCO4116, SCO2259 |
| **functional_only** | 1 | **SCO5082/actII-orf4** (act CSR) |
| **Total** | **9** | |

### 2.3 Rationale for actII-orf4 as functional SARP

ActII-ORF4 (SCO5082 / SC_RS27570) is the historical founder member of the SARP family, originally defined by its N-terminal HTH/winged-HTH DNA-binding motif and C-terminal activation domain. Despite this, Pfam HMM profiles for Trans_reg_C and BTAD do not detect matching domains in its 259-aa sequence. The NCBI GFF annotates it as a "TetR family transcriptional regulator," suggesting its domain architecture has diverged sufficiently from the canonical Pfam SARP models.

By adding it to the functional whitelist with `SARP_subtype = small` (consistent with its size and literature classification), v2 ensures that the ACT cluster is no longer missing its cognate SARP CSR.

---

## 3. Complete SARP List (v2)

| # | Gene ID | Old locus | Name | Length (aa) | Subtype | BGC | Role | Category |
|---|---------|-----------|------|------------|---------|-----|------|----------|
| 1 | SC_RS27585 | SCO5085 | — | 255 | small | act | CSR | structural_only |
| 2 | SC_RS31630 | SCO5877 | redD | 268 | small | red | CSR | both |
| 3 | SC_RS33690 | SCO6288 | — | 281 | small | cpk | CSR | both |
| 4 | SC_RS18200 | SCO3217 | — | 638 | medium | cda | CSR | both |
| 5 | SC_RS24295 | SCO4426 | afsR | 993 | large | none | global | structural_only |
| 6 | SC_RS06435 | SCO0898 | — | 660 | small_long | none | global | structural_only |
| 7 | SC_RS22750 | SCO4116 | — | 1102 | small_long | none | global | structural_only |
| 8 | SC_RS13320 | SCO2259 | — | 1322 | small_long | none | global | structural_only |
| 9 | **SC_RS27570** | **SCO5082** | **actII-orf4** | **259** | **small** | **act** | **CSR** | **functional_only** |

**Subtype distribution**: 4 small (including actII-orf4), 1 medium, 1 large, 3 small_long.

**CSR vs Global**: 5 cluster-situated (2 in act, 1 each in red/cda/cpk) + 4 global/pleiotropic.

---

## 4. Major BGC CSR SARPs (v2)

### 4.1 act cluster: Two SARPs with distinct origins

The act cluster now has **two SARP-category regulators**:

#### actII-orf4 (SC_RS27570 / SCO5082) — functional_only

| Metric | Value |
|--------|-------|
| SARP_category | **functional_only** |
| SARP_subtype | small (literature) |
| corr_act | **+0.954** |
| corr_red | +0.572 |
| corr_cda | +0.250 |
| corr_cpk | −0.025 |
| log2FC(3 vs 1) | +2.60 |
| Phase | **late** (delta_2_vs_1 = −0.84, delta_3_vs_2 = +3.32) |

actII-orf4 has the **highest act-specificity** among all SARPs (r_act = +0.954) with relatively low correlations to other BGCs. Its late-phase dynamics (sharp induction at M145_3) are consistent with its role as the direct transcriptional activator of act biosynthetic genes. Despite being annotated as "TetR family" in GFF, its functional behavior is unambiguously that of an act-specific SARP CSR.

#### SCO5085 (SC_RS27585) — structural_only

| Metric | Value |
|--------|-------|
| SARP_category | **structural_only** |
| SARP_subtype | small (domain-based) |
| corr_act | **+0.917** |
| corr_red | +0.873 |
| corr_cda | +0.690 |
| corr_cpk | +0.465 |
| log2FC(3 vs 1) | +6.44 |
| Phase | **late** (delta_2_vs_1 = +1.94, delta_3_vs_2 = +4.23) |

SCO5085 is a domain-confirmed SARP (DBD + BTAD) located within the act cluster. It has a high act correlation (r = +0.917) but also shows substantial red correlation (r = +0.87), suggesting broader regulatory influence. Its log2FC (+6.44) is much larger than actII-orf4's (+2.60), indicating a more dramatic transcriptional change.

**Relationship between actII-orf4 and SCO5085**: Both are late-phase, act-correlated SARPs within the act cluster, but actII-orf4 shows higher act-specificity (r = 0.954 vs 0.917) and lower pleiotropy (corr_red = 0.57 vs 0.87). actII-orf4 is the well-characterized pathway-specific activator; SCO5085 (annotated as "AfsR/SARP family transcriptional regulator" in GFF; annotated as "transport" in BGC_definition_manual.tsv) may serve an auxiliary or distinct regulatory role.

### 4.2 red cluster: redD / SCO5877 (SC_RS31630) — both

| Metric | Value |
|--------|-------|
| SARP_category | both |
| corr_act | +0.385 |
| corr_red | **+0.863** |
| corr_cda | **+0.973** |
| corr_cpk | **+0.955** |
| log2FC(3 vs 1) | +5.00 |
| Phase | **early_mid** (delta_2_vs_1 = +5.53, delta_3_vs_2 = −0.80) |

RedD (SCO5877) is a structural + functional SARP (category `both`). It correlates strongly with cda (r = +0.97) and cpk (r = +0.96), even exceeding its cognate red correlation (r = +0.86). Its early/mid-phase activation pattern is distinct from the late-phase act CSRs. Literature name: redD.

### 4.3 cda cluster: SCO3217 (SC_RS18200) — both

| Metric | Value |
|--------|-------|
| SARP_category | both |
| Domains | DBD + BTAD + NB-ARC |
| corr_act | +0.364 |
| corr_red | +0.852 |
| corr_cda | **+0.970** |
| corr_cpk | **+0.958** |
| log2FC(3 vs 1) | +5.76 |
| Phase | **early_mid** (delta_2_vs_1 = +6.62, delta_3_vs_2 = −1.11) |

SCO3217 is the sole medium SARP, containing an NB-ARC domain in addition to DBD + BTAD. Its strong cda correlation (r = +0.970) and early/mid-phase expression are consistent with a cluster-situated activator role. Literature name: cdaR (not in GFF as gene_name).

### 4.4 cpk cluster: SCO6288 (SC_RS33690) — both

| Metric | Value |
|--------|-------|
| SARP_category | both |
| corr_act | +0.360 |
| corr_red | +0.840 |
| corr_cda | **+0.958** |
| corr_cpk | **+0.940** |
| log2FC(3 vs 1) | +3.90 |
| Phase | **early_mid** (delta_2_vs_1 = +4.17, delta_3_vs_2 = −0.62) |

SCO6288 is a small SARP with `both` classification. It shows strong cpk correlation (r = +0.940) and early/mid-phase activation. Literature name: cpkN (not in GFF as gene_name).

### 4.5 Phase consistency across CSR SARPs (v2)

| CSR SARP | BGC | Phase | Pattern |
|----------|-----|-------|---------|
| actII-orf4 | act | late | delta_3_vs_2 = +3.32 >> delta_2_vs_1 = −0.84 |
| SCO5085 | act | late | delta_3_vs_2 = +4.23 >> delta_2_vs_1 = +1.94 |
| redD | red | early_mid | delta_2_vs_1 = +5.53, delta_3_vs_2 = −0.80 |
| SCO3217 | cda | early_mid | delta_2_vs_1 = +6.62, delta_3_vs_2 = −1.11 |
| SCO6288 | cpk | early_mid | delta_2_vs_1 = +4.17, delta_3_vs_2 = −0.62 |

The two-wave model from steps 06–08 is confirmed in v2: act CSRs (actII-orf4 and SCO5085) activate in the **late phase** (M145_3), while red/cda/cpk CSRs activate in the **early/mid phase** (M145_2). Including actII-orf4 strengthens this conclusion by providing the canonical act activator alongside SCO5085.

---

## 5. Global SARPs (v2)

All four global SARPs are `structural_only` (domain-detected, not on the functional whitelist). Their behavior is unchanged from v1.

### 5.1 AfsR (SC_RS24295 / SCO4426, large SARP)

| Metric | Value |
|--------|-------|
| Domains | DBD + BTAD + NB-ARC + TPR_12 |
| corr_act | **−0.825** |
| corr_red | −0.749 |
| log2FC(3 vs 1) | −0.918 |
| Phase | late (declining) |

AfsR shows negative correlations with all BGC scores, consistent with v1. Its activating role may be post-translational (phosphorylation-dependent) rather than transcription-level in this time course.

### 5.2 SCO2259 (SC_RS13320, small_long SARP)

| corr_act | corr_red | corr_cda | corr_cpk | mean_abs_corr | Phase |
|----------|----------|----------|----------|--------------|-------|
| −0.667 | **−0.921** | −0.891 | −0.768 | **0.812** | stable |

Highest mean |corr| among global SARPs. Strong negative correlations suggest a repressor-like role.

### 5.3 SCO4116 (SC_RS22750, small_long SARP)

| corr_act | corr_red | corr_cda | corr_cpk | mean_abs_corr | Phase |
|----------|----------|----------|----------|--------------|-------|
| −0.217 | −0.660 | −0.815 | −0.801 | 0.623 | early_mid |

Moderate negative correlations, strongest for cda/cpk.

### 5.4 SCO0898 (SC_RS06435, small_long SARP)

| corr_act | corr_red | corr_cda | corr_cpk | mean_abs_corr | Phase |
|----------|----------|----------|----------|--------------|-------|
| +0.448 | +0.264 | +0.104 | −0.074 | 0.222 | late |

Weak correlations. Not prioritized.

---

## 6. Interpretation Updates from v1 to v2

### 6.1 ACT cluster: "SARP gap" resolved

In v1, the act cluster appeared to lack a canonical SARP CSR because actII-orf4 was not detected by hmmsearch. SCO5085 was the only domain-confirmed SARP in the act region, but its BGC_definition_manual.tsv annotation as "transport" made its CSR status ambiguous.

In v2, actII-orf4 is explicitly included as a `functional_only` SARP. This resolves the apparent gap and restores the literature-consistent picture: **ActII-ORF4 is the founder SARP and the primary pathway-specific activator of actinorhodin biosynthesis**. SCO5085 (structural_only) serves as a secondary domain-confirmed SARP in the same cluster.

### 6.2 All four major BGCs have SARP CSRs

| BGC | CSR SARP(s) | Category |
|-----|-------------|----------|
| act | actII-orf4 (functional_only) + SCO5085 (structural_only) | 2 SARPs |
| red | redD / SCO5877 (both) | 1 SARP |
| cda | SCO3217 (both) | 1 SARP |
| cpk | SCO6288 (both) | 1 SARP |

Every major BGC has at least one SARP-type CSR, consistent with the SARP family's central role in *Streptomyces* secondary metabolism regulation.

### 6.3 Structural vs. functional: Why the distinction matters

The `SARP_category` flag distinguishes:
- **Domain-based confidence**: `structural_only` and `both` SARPs have Pfam domain evidence.
- **Literature-based annotation**: `functional_only` members (actII-orf4) are supported by extensive experimental characterization but lack detectable Pfam SARP domains, possibly due to sequence divergence.

This distinction is important for:
1. **Computational predictions**: Domain-based SARPs can be reliably detected in unannotated genomes. Functional-only SARPs require manual curation.
2. **Regulatory mechanism**: actII-orf4's divergent domain architecture (TetR-like in Pfam) may reflect a distinct DNA-binding mode compared to canonical SARP DBDs.
3. **Downstream analyses**: Motif prediction and structural modeling should treat domain-confirmed and functional-only SARPs differently.

---

## 7. Output Files (v2)

| File | Description | Rows |
|------|------------|------|
| `SARP_list_M145_v2.tsv` | Extended SARP list with dual classification | 9 |
| `gene_master_with_SARP_v2.tsv` | Full gene table with v2 SARP flags | 8,275 |
| `regulator_master_with_SARP_v2.tsv` | Regulator table with v2 SARP flags | 856 |
| `SARP_BGC_mapping_v2.tsv` | SARP-to-BGC assignment with category | 9 |
| `SARP_BGC_activity_summary_v2.tsv` | Expression/correlation summary per SARP | 9 |

New columns added (vs v1): `is_SARP_structural`, `is_SARP_functional`, `SARP_category`, `SARP_subtype_v2`.

---

## 8. Future Use

- `gene_master_with_SARP_v2.tsv` and `regulator_master_with_SARP_v2.tsv` can be directly used as Supplementary Tables in publications.
- The v2 SARP flags enable filtering by structural evidence, functional evidence, or both, supporting different analysis contexts (e.g., genome-wide computational scans vs. targeted experimental validation).
- For motif analysis (e.g., SARP binding site prediction), `structural_only` and `both` SARPs share canonical DBD domains suitable for PWM construction, while actII-orf4 (`functional_only`) requires separate treatment due to its divergent DNA-binding architecture.

---

## 9. Key Takeaways

1. **9 SARPs identified in v2** (8 structural + 1 functional-only): actII-orf4 (SCO5082) was added as a functional SARP to resolve its absence in v1's domain-based detection.
2. **Dual-layer classification** (`is_SARP_structural` / `is_SARP_functional` / `SARP_category`) provides transparent tracking of evidence type, distinguishing domain-confirmed SARPs from literature-curated ones.
3. **All four major BGCs (act, red, cda, cpk) now have explicit SARP CSRs**, with act having two (actII-orf4 + SCO5085). This restores the literature-consistent regulatory architecture.
4. **actII-orf4 has the highest act-specificity** (r_act = +0.954) among all SARPs, confirming its role as the primary act activator despite its TetR-like Pfam annotation. Its late-phase dynamics (delta_3_vs_2 = +3.32) match the act BGC activation pattern.
5. **The two-wave CSR activation model is reinforced**: act CSRs (late phase) vs. red/cda/cpk CSRs (early/mid phase), now supported by both structural and functional SARP evidence.

---

*Generated: 2026-01-28 | Pipeline: 09_SARP_integration v2 | Data: M145 RNA-seq (3 timepoints x 3 replicates)*
