# H7: CCGG MTase Paradox -- Resolution through Zero-Overlap Analysis

**Created:** 2026-02-24
**Last updated:** 2026-02-24
**Analysis directory:** `/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/30_CCGG_MTase_paradox/`

---

## Key Insights

| # | Finding | Significance |
|---|---------|-------------|
| 1 | **ZERO positional overlap** across all timepoints (Jaccard = 0.000 for T1-T2, T1-T3, T2-T3) | Completely rules out passive dilution; rules out maintenance methylation |
| 2 | **T1: 1,516 CCGG 4mC sites (core genome), T2: 486 (chromosomal arms), T3: 30 (scattered)** | Dramatic core-to-arm geographic shift between timepoints |
| 3 | **SC_RS19770/SC_RS36410 and 4 flanking genes co-induced at T3** as a defense island (LFC +1.1 to +5.3) | Dcm-like MTases are part of a larger defense system, not isolated enzymes |
| 4 | **PD-(D/E)XK nuclease** (SC_RS36385) and **Argonaute** (SC_RS36390) flank SC_RS36410 | Suggests R-M + Argonaute-based defense, not classical Dcm methylation |
| 5 | **SC_RS19765** (102 bp) is a truncated cytosine MTase pseudogene fragment, 53 bp upstream of SC_RS19770 | SC_RS19770 region has undergone duplication/rearrangement |
| 6 | Frequency distributions shift downward: T1 mean 82.7% --> T2 77.5% --> T3 67.7% | Each timepoint's de novo methylation is progressively less efficient |
| 7 | All known cytosine MTases are near-silent at T1 (SC_RS19770: 11 counts, SC_RS36410: 1.3 counts) | **No identified enzyme can account for 1,516 CCGG 4mC sites at T1** |
| 8 | T2 sites are median 4.6 kb from nearest T1 site (not adjacent) | Excludes "spreading" or "sliding" from T1 positions |

---

## 1. Background and Hypothesis

### The Paradox

In H5, a paradox was identified regarding CCGG/TGGCCGGC-associated 4mC sites:

- **T1** has 1,516 CCGG-related 4mC sites (76.3% of all 4mC)
- The only candidate Dcm-like MTases (SC_RS19770 and SC_RS36410) are near-silent at T1 (mean 11 and 1.3 normalized counts)
- At **T3**, SC_RS36410 expression surges (LFC = +5.34***) but 4mC sites plummet to just 30

### Hypotheses Tested

| Hypothesis | Description |
|-----------|-------------|
| (a) Unknown MTase | T1 4mC is produced by an unidentified methyltransferase |
| (b) Passive dilution | T1 4mC is semi-conservatively maintained; MTase expression not needed |
| (c) 5mC vs 4mC confusion | SC_RS19770/36410 catalyze 5mC (Dcm-like), and 4mC comes from another enzyme |
| (d) Detection artifact | T1 4mC sites are false positives from replication-associated signals |

---

## 2. Site Overlap Analysis (Verification 1)

### Result: Complete Non-Overlap

All CCGG-related 4mC sites are **completely unique to each timepoint**:

| Comparison | Sites | Jaccard Index |
|-----------|-------|---------------|
| T1 total | 1,516 | -- |
| T2 total | 486 | -- |
| T3 total | 30 | -- |
| T1 intersect T2 | **0** | **0.000** |
| T1 intersect T3 | **0** | **0.000** |
| T2 intersect T3 | **0** | **0.000** |

This zero overlap is consistent across all sub-motifs:

| Motif | T1 | T2 | T3 | T1-T2 overlap | T1-T3 overlap |
|-------|----|----|-----|---------------|---------------|
| TGGCCGGC | 1,289 | 407 | 21 | 0 | 0 |
| CCGG | 101 | 39 | 6 | 0 | 0 |
| GGCCGG | 10 | 6 | 0 | 0 | 0 |
| CCGC | 116 | 34 | 3 | 0 | 0 |

**Interpretation:** This definitively **rules out hypothesis (b) -- passive dilution**. If methylation were semi-conservatively maintained, a subset of T1 sites would persist at T2 with reduced frequency. Instead, every T2 site is at a position never methylated at T1.

### Figure: Venn Diagram

![Venn Diagram](../11_epigenome_integration/analysis/30_CCGG_MTase_paradox/figures/CCGG_site_venn_diagram.pdf)

Three non-overlapping circles for T1 (1,516), T2 (486), and T3 (30). Both the all-motif and TGGCCGGC-only panels show identical zero-overlap pattern.

---

## 3. Frequency Distribution Analysis (Verification 2)

### Frequency Statistics

| Timepoint | n | Mean (%) | Median (%) | SD | % >= 80% | % < 60% |
|-----------|---|---------|-----------|-----|----------|---------|
| T1 | 1,516 | 82.7 | 85.5 | 12.0 | 66.0% | 6.4% |
| T2 | 486 | 77.5 | 79.1 | 11.8 | 47.3% | 10.7% |
| T3 | 30 | 67.7 | 64.9 | 12.7 | 20.0% | 33.3% |

The progressive frequency decrease (T1 > T2 > T3) is statistically significant:
- **T1 vs T2 frequency:** Mann-Whitney U = 337,688, p = 1.06e-18 (T1 > T2)

### Per-Motif Frequency

TGGCCGGC (dominant motif) follows the same pattern: T1 mean 83.7% --> T2 78.4% --> T3 68.3%.

### Figure: Frequency Distributions

![Frequency Distribution](../11_epigenome_integration/analysis/30_CCGG_MTase_paradox/figures/CCGG_frequency_distribution.pdf)

- Upper-left: All CCGG-related 4mC sites, histograms shifted left from T1 to T3
- Upper-right: TGGCCGGC-only, same pattern
- Lower-left: CDF curves showing T3 sites have substantially lower frequencies
- Lower-right: Bar chart showing dramatic decline across all motif subtypes

**Interpretation:** Since there is zero positional overlap, the frequency decline across timepoints does NOT reflect dilution of the same sites. Instead, each timepoint's de novo methylation event is less efficient (fewer sites, lower frequency), consistent with declining enzymatic activity or substrate availability.

---

## 4. Genomic Context Analysis (Verification 4)

### Dramatic Core-to-Arm Geographic Shift

| Timepoint | Core genome (1.5-6.5 Mb) | Chromosomal arms | % Arm |
|-----------|--------------------------|-------------------|-------|
| T1 | 1,174 (77.4%) | 342 (22.6%) | 22.6% |
| T2 | 54 (11.1%) | **432 (88.9%)** | **88.9%** |
| T3 | 9 (30.0%) | 21 (70.0%) | 70.0% |

This is one of the most striking findings: **T1 4mC sites are concentrated in the core genome, but T2 sites shift almost exclusively to the chromosomal arms.**

### Figure: Genome-Wide Density

![Genome Density](../11_epigenome_integration/analysis/30_CCGG_MTase_paradox/figures/CCGG_genome_density.pdf)

- Top panel: Genomic CCGG density (relatively uniform, ~800-1000 per 50 kb)
- T1 panel: Sites distributed across the core genome (1.5-6.5 Mb), mirroring CCGG density
- T2 panel: Sites concentrated at 0-1.5 Mb and 6.5-8.7 Mb (chromosomal arms)
- T3 panel: Only 30 scattered sites remaining

### Intragenic Bias

| Timepoint | Intragenic | Intergenic |
|-----------|-----------|------------|
| T1 | 1,423 (93.9%) | 93 (6.1%) |
| T2 | 449 (92.4%) | 37 (7.6%) |
| T3 | 30 (100.0%) | 0 (0.0%) |

The very high intragenic fraction (~93%) is expected given the dense coding nature of the *S. coelicolor* genome (87% coding).

### T2 Site Proximity to T1 Sites

T2 TGGCCGGC sites are NOT adjacent to former T1 sites:
- Mean distance to nearest T1 site: **8,872 bp**
- Median distance: **4,602 bp**
- Within 100 bp: only 14/407 sites (3.4%)
- Within 1 kb: only 54/407 sites (13.3%)

This rules out a "spreading" or "sliding" mechanism from T1 positions.

---

## 5. Methyltransferase Catalogue (Verification 3)

### All DNA Cytosine MTases in the Genome

| Locus tag | Product | Size | T1 expr | T2 expr | T3 expr | T3vsT1 LFC |
|-----------|---------|------|---------|---------|---------|------------|
| SC_RS19765 | DNA cytosine methyltransferase | 102 bp (pseudogene) | N/A | N/A | N/A | N/A |
| SC_RS19770 | DNA cytosine methyltransferase | 975 bp | 11.0 | 4.6 | 62.1 | +2.32*** |
| SC_RS36410 | DNA cytosine methyltransferase | 1,272 bp | 1.3 | 4.2 | 64.9 | +5.34*** |

Additional non-cytosine DNA MTase candidates:

| Locus tag | Product | T1 expr | T2 expr | T3 expr |
|-----------|---------|---------|---------|---------|
| SC_RS17645 | N-6 DNA methylase (Type I HsdM) | 245.1 | 51.7 | 145.1 |
| SC_RS19670 | DNA-methyltransferase | 3.8 | 2.9 | 48.0 |
| SC_RS36625 | DNA-methyltransferase | 2.7 | 3.9 | 46.0 |
| SC_RS03950 | Class I SAM-dep. DNA MTase | 130.2 | 160.9 | 102.3 |
| SC_RS24685 | Class I SAM-dep. DNA MTase | 324.8 | 33.8 | 66.1 |
| SC_RS28835 | BREX-2 PglX (adenine MTase) | 604.3 | 395.2 | 381.7 |
| SC_RS35335 | BREX-2 PglX (adenine MTase) | 1,666.3 | 1,989.3 | 979.8 |

**Key finding:** No DNA cytosine methyltransferase has sufficient expression at T1 to explain 1,516 4mC sites. The 129 total MTases in the GFF include many non-DNA-modifying enzymes (RNA MTases, metabolic MTases). After filtering for DNA-specific cytosine MTases, only SC_RS19770 and SC_RS36410 remain -- both near-silent at T1.

SC_RS19765 is a 102 bp pseudogene fragment, only 53 bp upstream of SC_RS19770, likely resulting from a duplication/rearrangement event.

---

## 6. Operon Context of Dcm-like MTases (Verification 5)

### SC_RS36410 Defense Island

SC_RS36410 is embedded in a co-induced defense island:

| Locus tag | Product | T1 | T2 | T3 | T3vsT1 LFC |
|-----------|---------|----|----|-----|-------------|
| SC_RS36385 | **PD-(D/E)XK nuclease** | 33.4 | 27.4 | 77.5 | +1.13** |
| SC_RS36390 | **Argonaute/Piwi** | 17.3 | 20.3 | 114.7 | +2.65*** |
| SC_RS36395 | DUF5655 domain | 17.2 | 46.7 | 73.4 | +2.03*** |
| SC_RS36400 | SPDY domain | 8.3 | 6.2 | 57.8 | +2.61*** |
| SC_RS36405 | **DnaB-like helicase** | 4.7 | 6.1 | 59.3 | +3.48*** |
| **SC_RS36410** | **DNA cytosine MTase** | **1.3** | **4.2** | **64.9** | **+5.34***** |

The presence of:
- **PD-(D/E)XK nuclease** (potential restriction endonuclease)
- **Argonaute/Piwi** (RNA/DNA-guided nuclease for phage defense)
- **DnaB-like helicase** (unwinding for defense processing)

strongly suggests this is a **defense island** combining R-M-like methylation with Argonaute-based immunity. The entire region is coordinately induced at T3.

### SC_RS19770 Region

| Locus tag | Product | T1 | T2 | T3 | T3vsT1 LFC |
|-----------|---------|----|----|-----|-------------|
| SC_RS19750 | Type IV secretory system | -- | -- | -- | -- |
| SC_RS19760 | **DnaB-like helicase** | 10.2 | 4.5 | 89.0 | +3.00*** |
| SC_RS19765 | DNA cytosine MTase (pseudogene) | N/A | N/A | N/A | N/A |
| **SC_RS19770** | **DNA cytosine MTase** | **11.0** | **4.6** | **62.1** | **+2.32***** |
| SC_RS19775 | DUF4913 domain | -- | -- | -- | -- |

SC_RS19770 also has a DnaB-like helicase neighbor (SC_RS19760) that is co-induced at T3 with similar kinetics.

### No REase Found

No gene annotated as "restriction endonuclease" was found within 10 kb of either SC_RS19770 or SC_RS36410 by keyword search. However, the PD-(D/E)XK nuclease (SC_RS36385) near SC_RS36410 belongs to a superfamily that includes many restriction endonucleases and could function as a REase.

### Figure: Expression Timeline

![MTase Expression](../11_epigenome_integration/analysis/30_CCGG_MTase_paradox/figures/MTase_expression_timeline_focused.pdf)

- Panel A: All DNA MTases -- SC_RS19770 and SC_RS36410 show the characteristic "silent at T1/T2, induced at T3" pattern
- Panel B: SC_RS36385-SC_RS36410 defense island genes are coordinately induced at T3

---

## 7. Paradox Resolution

### What the data show

1. **Zero positional overlap**: Each timepoint has entirely different CCGG 4mC sites
2. **Geographic shift**: T1 sites are in the core genome; T2 sites are in chromosomal arms
3. **No maintenance**: The pattern is incompatible with semi-conservative maintenance methylation
4. **No identified T1 enzyme**: All cytosine MTases are near-silent at T1
5. **T3 defense induction**: The Dcm-like MTases are part of defense islands induced at T3, but 4mC sites are nearly gone by T3

### Hypothesis evaluation

| Hypothesis | Verdict | Evidence |
|-----------|---------|---------|
| **(b) Passive dilution** | **REJECTED** | Zero overlap definitively rules this out. Passive dilution requires persistence of a subset of sites. |
| **(d) Detection artifact** | **PARTIALLY SUPPORTED** | The stochastic, non-overlapping pattern with high frequencies (>80%) is consistent with systematic false positives, but the clear frequency hierarchy (T1 > T2 > T3) and geographic structure argue against pure randomness. |
| **(a) Unknown MTase** | **POSSIBLE but problematic** | No cytosine MTase has sufficient T1 expression. Would require a completely uncharacterized enzyme or post-transcriptional regulation. |
| **(c) 5mC vs 4mC confusion** | **REMAINS VIABLE** | If Dcm-like MTases produce 5mC (as expected for Dcm family), nanopore basecalling might misclassify some 5mC signals as 4mC in specific sequence contexts. The geographic shift could reflect differential 5mC:4mC misclassification rates in different genomic regions. |

### New hypothesis: Stochastic de novo methylation by transiently active enzyme(s)

The zero-overlap, declining-count, geographic-shift pattern is best explained by:

1. **Transient de novo methylation**: At each timepoint, a methyltransferase (possibly with very brief expression windows not captured by our RNA-seq sampling) methylates a stochastic subset of CCGG sites
2. **No maintenance**: Once established, the methylation is NOT maintained across cell divisions
3. **Geographic preference shifts**: The enzyme or its accessibility changes between timepoints (core at T1, arms at T2)
4. **Declining activity**: The enzyme's overall activity diminishes across the growth curve

Alternatively, the **technical artifact hypothesis** gains credibility from the observation that:
- The pattern mirrors what would be expected from stochastic false positives that vary with sequencing depth and DNA quality
- However, the clear geographic structure (core vs arm shift) argues against a purely random technical artifact

---

## 8. Output Files

### Tables

| File | Description |
|------|-------------|
| `tables/CCGG_site_overlap_between_timepoints.tsv` | T1/T2/T3 site overlap counts |
| `tables/CCGG_site_overlap_per_motif.tsv` | Per-motif (TGGCCGGC, CCGG, GGCCGG, CCGC) overlap and Jaccard indices |
| `tables/CCGG_frequency_distribution.tsv` | Frequency statistics per timepoint |
| `tables/all_methyltransferases_GFF.tsv` | All 129 GFF MTases with expression data |
| `tables/CCGG_genomic_context.tsv` | Genomic context (gene, region) for each CCGG 4mC site |
| `tables/Dcm_operon_context.tsv` | 10 kb flanking genes for SC_RS19770, SC_RS36410, SC_RS19670, SC_RS36625 |

### Figures (PDF + SVG)

| File | Description |
|------|-------------|
| `figures/CCGG_site_venn_diagram` | Venn diagram of T1/T2/T3 CCGG site overlap (zero for all) |
| `figures/CCGG_frequency_distribution` | Frequency histograms, CDF, and per-motif bar chart |
| `figures/CCGG_genome_density` | Genome-wide density showing T1-core to T2-arm shift |
| `figures/MTase_expression_timeline_focused` | DNA MTase expression timeline + SC_RS36410 defense island |
| `figures/MTase_expression_timeline` | All 129 MTases (comprehensive version) |

---

## 9. Implications for the Project

### For the AAGCCCG story
The zero-overlap pattern seen here for CCGG 4mC **parallels the AAGCCCG result from H5** (also zero overlap). This suggests that the complete positional turnover between timepoints is a general property of methylation in this system, not specific to one motif.

### For the "gatekeeper" model
If methylation positions change completely between timepoints, the "gatekeeper" model must be interpreted as: methylation acts as a **transient, position-variable** gate that opens different genomic regions at different developmental stages, rather than a stable lock at fixed positions.

### For the 4mC vs 5mC question
The finding that Dcm-like MTases (which typically produce 5mC) are the only CCGG-targeting enzymes, yet the modifications are called as 4mC, strengthens hypothesis (c). Oxford Nanopore basecallers have known difficulties distinguishing 4mC from 5mC, particularly in GC-rich sequence contexts characteristic of Streptomyces.

### For defense island biology
The discovery that SC_RS36410 is co-induced with an Argonaute and PD-(D/E)XK nuclease reveals a previously unrecognized **composite defense system** in *S. coelicolor*. This system combines:
- **R-M-like activity** (MTase + nuclease)
- **Argonaute-mediated defense** (RNA/DNA-guided cleavage)
- **Helicase activity** (DNA unwinding for substrate access)

This defense island is activated specifically at T3 (stationary phase), consistent with the known upregulation of phage defense during stress.
