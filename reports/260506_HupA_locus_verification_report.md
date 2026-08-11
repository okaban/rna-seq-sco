# HupA Locus Verification Report — *S. coelicolor* M145
**Date:** 2026-05-06  
**Context:** Manuscript discussion cites Salerno et al. 2009 for HupA (vegetative) and HupS (sporulation-specific). QC flagged a discrepancy between SCO4950 (currently annotated as *narI*) and SCO2950 (annotated as "HU family DNA-binding protein").

---

## Summary

| Question | Answer |
|---|---|
| Correct SCO for HupA? | **SCO2950** (SC_RS16855) |
| SCO2950 matches constitutive/vegetative pattern? | **Yes** — highest expression at T1, monotonically declining |
| Should manuscript cite SCO2950 or SCO4950 for HupA? | **SCO2950** |
| Is this a discrepancy requiring correction? | **Yes — SCO4950 must not be cited as HupA** |

---

## 1. Literature Confirmation: Salerno et al. 2009

Salerno et al. 2009 (*J. Bacteriol.* 191:6489–6498; PMID 19717607) explicitly assigned locus tags in *S. coelicolor* M145:

- **HupA = SCO2950** — constitutive, expressed in vegetative hyphae; conventional HU-type protein similar to *E. coli* HUα/HUβ
- **HupS = SCO5556** — developmentally regulated; two-domain protein; strongly upregulated in sporogenic aerial hyphae in a *whiA/whiG/whiI*-dependent manner; required for spore nucleoid condensation and heat resistance

This assignment is independently corroborated by the STRING database (entry 100226.SCO5556 = hup2/HupS) and by a 2024 follow-up study (Microbial Cell Factories, doi:10.1186/s12934-024-02549-0) that uses the same locus tags.

> **Note on PMID:** The question references PMID 19549799. The confirmed Salerno 2009 HupA/HupS paper is PMID 19717607. If the manuscript cites 19549799, that PMID should be verified — it may refer to a related Salerno et al. paper on *Streptomyces* Dps-like proteins.

---

## 2. Project Annotation Data

From `05_annotation/analysis/05_annotation_260128_v1/tables/gene_annotation_basic.tsv` and `gene_master_DESeq2.tsv`:

| Gene ID | Old Locus | Gene Name | Product (RefSeq NC_003888.3) | Verdict |
|---|---|---|---|---|
| SC_RS16855 | SCO2950 | — | **HU family DNA-binding protein** | ✅ **HupA** |
| SC_RS26915 | SCO4950 | narI | Respiratory nitrate reductase subunit gamma | ❌ Not HupA |
| SC_RS29990 | SCO5556 | — | HU family DNA-binding protein | HupS (per Salerno 2009) |

SCO4950 / SC_RS26915 carries the gene name `narI` and product description "respiratory nitrate reductase subunit gamma" in all current annotation tables. It is unambiguously a metabolic enzyme, not a nucleoid-associated protein.

---

## 3. Expression Data

Normalized counts from DESeq2 (columns: T1-rep1, T1-rep2, T1-rep3 | T2-rep1, T2-rep2, T2-rep3 | T3-rep1, T3-rep2, T3-rep3):

### SCO2950 (SC_RS16855) — HU family DNA-binding protein = **HupA**

| Timepoint | Rep1 | Rep2 | Rep3 | Mean |
|---|---|---|---|---|
| T1 (vegetative) | 21,517 | 22,706 | 22,881 | **22,368** |
| T2 (transition) | 7,140 | 7,073 | 6,485 | **6,899** |
| T3 (sporulation) | 2,473 | 4,251 | 2,071 | **2,932** |

- log₂FC T2 vs T1 = −1.65 (padj = 2.8 × 10⁻¹²) — significant drop at transition
- log₂FC T3 vs T1 = **−2.90** (padj = 2.7 × 10⁻³⁴) — 7.5-fold decrease by sporulation
- **Interpretation:** Highest absolute expression of any nucleoid-associated protein at T1. Monotonically declining. This is fully consistent with Salerno 2009's description of HupA as the dominant vegetative NAP.

### SCO4950 (SC_RS26915) — *narI*, nitrate reductase γ-subunit (**NOT HupA**)

| Timepoint | Rep1 | Rep2 | Rep3 | Mean |
|---|---|---|---|---|
| T1 (vegetative) | 49 | 41 | 49 | **46** |
| T2 (transition) | 1,321 | 1,946 | 1,534 | **1,600** |
| T3 (sporulation) | 149 | 201 | 138 | **163** |

- log₂FC T2 vs T1 = **+5.09** (padj = 3.2 × 10⁻⁹⁸) — strongly induced at transition
- Expression profile: very low at T1, peaks at T2, drops at T3
- **Interpretation:** Classic nitrogen metabolism gene expression pattern, not that of a constitutive NAP. Completely incompatible with HupA biology.

### SCO5556 (SC_RS29990) — HU family DNA-binding protein = **HupS** (per Salerno 2009)

| Timepoint | Rep1 | Rep2 | Rep3 | Mean |
|---|---|---|---|---|
| T1 (vegetative) | 1,876 | 1,947 | 1,906 | **1,910** |
| T2 (transition) | 431 | 721 | 730 | **627** |
| T3 (sporulation) | 335 | 595 | 331 | **420** |

- log₂FC T3 vs T1 = −1.56 (padj = 2.8 × 10⁻⁹) — moderate decline
- **Note on HupS vs SCO1175:** The project's prior HiC/FIRE integration report (260506) identified SCO1175 (SC_RS07825) as the strongly sporulation-induced HupS candidate (log₂FC T3 vs T1 = +5.15). SCO1175 is currently annotated as a pseudogene in RefSeq, and SCO5556 as the second HU-family coding gene. This discrepancy may reflect re-annotation between the 2009 genome and the current RefSeq assembly. **The HupS locus assignment warrants separate verification** (see Section 5).

---

## 4. Verdict: Manuscript Correction Required

The manuscript must **not cite SCO4950 as HupA**. The evidence is unambiguous:

1. Salerno et al. 2009 assigns hupA = **SCO2950**, not SCO4950.
2. SCO4950 is *narI* (nitrate reductase) in all current RefSeq annotation tables — confirmed in this project's own `gene_annotation_basic.tsv` and `gene_master_DESeq2.tsv`.
3. SCO2950 has the correct molecular annotation ("HU family DNA-binding protein"), the correct GO terms (GO:0030261 chromosome condensation, GO:0003677 DNA binding, GO:0030527 structural constituent of chromatin), and the correct expression profile (constitutively high in T1 vegetative hyphae, declining through sporulation).

**The Discussion should be updated to cite SCO2950 (SC_RS16855) as HupA.**

If the original error arose from a stale SCO numbering (e.g., an older annotation where SCO4950 carried the hupA annotation), that should be explicitly noted in a revision comment or author note, as it could affect reproducibility if other groups try to locate the gene.

---

## 5. Open Question: HupS Locus

The Salerno 2009 paper assigns hupS = SCO5556. However, this project's expression data shows SCO5556 has moderate constitutive expression that **declines** during sporulation — inconsistent with sporulation-specific induction. Conversely, SCO1175 (currently annotated as a pseudogene) shows dramatic sporulation induction (T1 mean ~1.3 → T3 mean ~26 normalized counts; log₂FC +5.15, padj = 1.1 × 10⁻⁴).

**Possible explanations:**
- SCO1175 and SCO5556 may have had their annotations swapped between the 2003 genome release (used by Salerno 2009) and the current RefSeq NC_003888.3.
- Alternatively, Salerno's "SCO5556 = hupS" was based on gene identity/synteny rather than directly on the SCO numbers in the published genome, and a re-annotation moved one HU-family gene to a new locus tag.

**Recommendation:** Before finalizing the manuscript, cross-check SCO5556 and SCO1175 protein sequences against the Salerno 2009 Supplementary Figure S1 (HupA/HupS alignment). The sporulation-induced protein (whiA/whiG/whiH-dependent) should be the two-domain HupS; its sequence will be distinct from the single-domain HupA. This sequence-level check will definitively resolve the HupS locus tag.

---

## 6. Recommended Manuscript Changes

1. **Replace** all instances of "SCO4950" cited as HupA with **"SCO2950 (hupA)"**.
2. **Add** the RefSeq locus tag SC_RS16855 for precision.
3. **Verify** HupS locus (SCO5556 vs SCO1175) by protein sequence comparison before submission.
4. **Check PMID** cited for Salerno 2009 — confirm whether 19549799 or 19717607 is the intended reference.

---

## Sources

- Salerno et al. 2009, *J. Bacteriol.* 191:6489–6498, PMID 19717607  
  https://journals.asm.org/doi/10.1128/jb.00709-09
- PubMed entry: https://pubmed.ncbi.nlm.nih.gov/19717607/
- STRING DB SCO5556 (hup2): https://string-db.org/network/100226.SCO5556
- 2024 follow-up: Microbial Cell Factories, https://link.springer.com/article/10.1186/s12934-024-02549-0
- Project annotation: `05_annotation/analysis/05_annotation_260128_v1/tables/gene_annotation_basic.tsv`
- Project DESeq2 master table: `05_annotation/analysis/05_annotation_260128_v1/tables/gene_master_DESeq2.tsv`
