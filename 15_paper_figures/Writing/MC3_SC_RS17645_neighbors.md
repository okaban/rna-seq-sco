# MC3: SC_RS17645 (SCO3104) — Cognate HsdR/HsdS Confirmation

**Reviewer concern:** SC_RS17645をソロMTaseと主張しているが、cognate HsdR/HsdSの有無を確認していない。

**Analysis date:** 2026-05-17  
**Data source:** `11_epigenome_integration/analysis/73_SC_RS17645/` (GCF_000203835.1 GFF3, InterPro, REBASE, tblastn)  
**Genome:** *S. coelicolor* A3(2) M145, NC_003888.3 (AL645882.2 = GenBank)

---

## 1. Gene identity and annotation

| Field | Value | Source |
|-------|-------|--------|
| Locus tag | SC_RS17645 | RefSeq GFF3 |
| Old locus | SCO3104 | GFF3 |
| Coordinates | 3,399,763–3,401,802 (−) | GFF3 |
| Product | N-6 DNA methylase | RefSeq |
| Protein | WP_011028768.1 (679 aa) | NCBI |
| REBASE name | M.ScoA3ORF3104P | REBASE |
| Predicted target | AAGCCCG (6mA, novel 7-mer) | This work |
| InterPro family | **IPR052916** Type I RE MTase Subunit (aa 1–484) | InterPro |
| Pfam catalytic | PF02384 N6_Mtase (aa 167–393) | Pfam |
| TRD domain | **IPR044946** Type I R-M DNA specificity domain (aa 539–669) | InterPro |

**Key domain architecture:** SC_RS17645 contains both a catalytic N6-MTase domain (Pfam PF02384) and a C-terminal Target Recognition Domain (TRD, IPR044946) characteristic of Type I R-M system MTase/specificity subunits (HsdM+HsdS fused). This architecture resembles a *fused* HsdM–HsdS module rather than a free-standing simple MTase.

---

## 2. ±50 kb neighborhood scan for HsdR/HsdS

Scan of GCF_000203835.1 GFF3, genomic region 3,350,000–3,450,000 bp (±50 kb around SC_RS17645 at 3,399,763):

| Position rel. to target | Locus | Product | HsdR/HsdS? |
|------------------------|-------|---------|------------|
| −5 | SC_RS17620 (SCO3099) | Cytochrome P450 family | No |
| −4 | SC_RS17625 (SCO3100) | NTP pyrophosphohydrolase | No |
| −3 | SC_RS17630 (SCO3101) | SurA N-terminal domain | No |
| −2 | SC_RS17635 (SCO3102) | pkaE (Ser/Thr protein kinase) | No |
| −1 | SC_RS17640 (SCO3103) | Hypothetical protein | No |
| **0** | **SC_RS17645 (SCO3104)** | **N-6 DNA methylase** | **TARGET** |
| +1 | SC_RS17650 (SCO3105) | Hypothetical protein | No |
| +2 | SC_RS17655 (SCO3106) | Lipoprotein | No |
| +3 | SC_RS17660 (SCO3107) | **HNH endonuclease family protein** | Endonuclease, not HsdR |
| +4 | SC_RS17665 (SCO3108) | Antitoxin | No |
| +5 | SC_RS17670 (SCO3109) | mfd (transcription-repair coupling factor) | No |

**±50 kb result: No HsdR (Type I restriction enzyme) or HsdS (specificity subunit) annotated.**

---

## 3. Whole-genome HsdR/HsdS scan

Systematic scan of all 7,769 CDS product names in GCF_000203835.1 GFF3 for terms: "HsdR", "HsdS", "HsdM", "Type I restriction", "restriction-modification specificity", "specificity subunit":

**Result: No gene in the M145 genome is annotated as HsdR, HsdS, or a cognate Type I restriction subunit.**

The HNH endonuclease at SC_RS17660 (+3; 642 bp away) is a different structural class (HNH family, not REase), and no function prediction links it to Type I RM-system restriction.

---

## 4. REBASE and conservation context

- REBASE name: **M.ScoA3ORF3104P** — listed as an orphan methyltransferase (no cognate R recorded in REBASE for M145)  
- tblastn across 833 *Streptomyces* genomes: 666/833 (80%) have a homolog (≥30% identity); high-identity homologs (≥99%) restricted to *S. coelicolor* / *S. lividans* / *S. violaceoruber* clade  
- In all surveyed near-identical relatives (S. lividans, S. violaceoruber, S. violaceolatus), the genomic neighborhood context is conserved with no HsdR/HsdS flanking genes

---

## 5. Nuance: Type I MTase subunit domain architecture

SC_RS17645 carries **IPR044946 (Type I R-M DNA specificity domain / TRD)** at its C-terminus (aa 539–669), in addition to the catalytic N6-MTase domain. TRD domains are normally found in:
(a) full HsdS specificity subunits of classical Type I systems (HsdR–HsdM–HsdS triplet)
(b) fused HsdM+HsdS units in Type I systems lacking a separate HsdS gene

The presence of a TRD in SC_RS17645 is thus notable, but:
- No HsdR partner is encoded anywhere in the M145 genome  
- No HsdS-standalone gene is found in the neighborhood or genome-wide  
- Expression profile: SC_RS17645 is significantly downregulated at T2 (log₂FC=−2.19, padj=6.5×10⁻¹⁶) and partially recovers at T3 (+1.45 vs T2), a pattern inconsistent with constitutively expressed RM defense systems but consistent with a regulatory MTase

**Interpretation:** SC_RS17645 appears to be an **ancestrally Type-I-like MTase that has been decoupled from its restriction partner** — either through gene loss or evolutionary divergence to a solo-MTase function. This is not unprecedented; "orphan MTases" in bacterial genomes frequently arise from RM system decay.

---

## 6. Judgment and recommended manuscript modification

**Decision: Maintain "solo MTase" claim, but add domain architecture caveat and GenBank citation.**

### Updated claim language:

Replace:
> "SC_RS17645 (SCO3104) encodes a solo N-6 adenine DNA methyltransferase..."

With:
> "SC_RS17645 (SCO3104) encodes an N-6 adenine DNA methyltransferase (REBASE: M.ScoA3ORF3104P) whose protein architecture includes both a catalytic N6-MTase domain (Pfam PF02384) and a C-terminal Target Recognition Domain (TRD; IPR044946), characteristic of Type I R-M system MTase-specificity subunits. However, no cognate HsdR (restriction endonuclease) or standalone HsdS gene is annotated within ±50 kb or elsewhere in the M145 genome (NCBI GenBank AL645882.2; GFF3 product-name scan), and REBASE lists no cognate restriction partner for this locus. SC_RS17645 is therefore classified as a functional solo MTase—likely representing an ancestrally Type-I-derived enzyme that has been decoupled from its restriction partner through genomic decay or evolutionary repurposing."

### Methods addition (1 sentence in Epigenome Analysis section):

> "The absence of cognate HsdR/HsdS partners for SC_RS17645 was confirmed by (i) systematic GFF3 product-name scanning across all 7,769 CDSs of GCF_000203835.1 (NC_003888.3), (ii) manual inspection of the ±50 kb chromosomal neighborhood (positions 3,350,000–3,450,000), and (iii) REBASE cross-reference (M.ScoA3ORF3104P has no associated restriction entry); no Type I restriction subunit or specificity subunit was identified."

---

## 7. Summary recommendation

| Criterion | Result | Implication |
|-----------|--------|-------------|
| HsdR within ±50 kb | **None** | Main claim holds |
| HsdS within ±50 kb | **None** | Main claim holds |
| HsdR/HsdS anywhere in genome | **None** | Main claim holds |
| Protein domain architecture | TRD domain present | Add caveat in Results/Methods |
| REBASE classification | Orphan MTase | Cite explicitly |
| Nearest endonuclease | HNH (+3; unrelated class) | Not cognate |

**Verdict:** "Solo MTase" framing is **justified and maintainable**. The TRD domain should be acknowledged to pre-empt reviewer concern, framed as ancestral architecture consistent with gene-loss evolution rather than evidence for an active Type I system.

---

## Output files

- `11_epigenome_integration/analysis/73_SC_RS17645/SC_RS17645_annotation_summary.tsv` — full annotation table  
- `11_epigenome_integration/analysis/73_SC_RS17645/gene_neighborhood.tsv` — ±5 gene neighborhood  
- `11_epigenome_integration/analysis/73_SC_RS17645/tblastn_conservation_summary.txt` — cross-species conservation  
