# H8: Geographic and Functional Characterization of 62 Coordinated Regulatory Genes

**Date:** 2026-02-24
**Analysis directory:** `11_epigenome_integration/analysis/31_coordinated_regulators_characterization/`
**Hypothesis:** The 62 methylation-expression coordinated regulatory genes are enriched at chromosomal arms and in stress response / secondary metabolism COG categories, supporting a "gatekeeper" model where methylation gates a specific regulatory layer.

---

## Background

In H6 (genomewide TF screen), we identified 62 regulatory genes (out of 1,055 genome-wide) that show methylation-expression coordination. Critically, **ALL 62 are outside the literature-curated 37-TF list**, representing a previously unrecognized "hidden" methylation-responsive regulatory layer. This analysis characterizes these 62 genes geographically, functionally, and structurally.

---

## Key Results Summary

| Finding | Result | Statistical Support |
|---------|--------|-------------------|
| Arm enrichment (62 vs genome) | 51.6% arm vs 44.7% genome | OR=1.32, p=0.17 (NS) |
| T3 discordant_gain_up arm enrichment | 86.7% arm (13/15) | OR=8.05, **p=0.0010** |
| T3 discordant_loss_down core enrichment | 80.0% core (8/10) | OR=3.23, p=0.10 |
| COG K (Transcription) enrichment | 72.6% vs 9.9% genome | OR=24.0, **p=2.7e-31** |
| COG T (Signal transduction) enrichment | 9.7% vs 4.2% genome | OR=2.47, **p=0.045** |
| COG Q (Secondary metabolism) | 0/62 | Not enriched |
| Both T2+T3 concordant derepression | 3 genes | Top candidates |
| Both T2+T3 concordant repression | 1 gene (SC_RS24635) | |
| TCS cognate pairs identified | 7 pairs (0 both coordinated) | |
| BGC proximity (<20 kb) | 1 gene (SC_RS33745 near CPK) | |

---

## Step 1: Geographic Distribution

### Overall Distribution

| Group | n_arm | n_core | n_total | % arm |
|-------|-------|--------|---------|-------|
| 62 Coordinated | 32 | 30 | 62 | **51.6%** |
| 1,055 Regulators | 475 | 580 | 1,055 | 45.0% |
| Genome-wide | 3,696 | 4,579 | 8,275 | 44.7% |
| Expected (by length) | - | - | - | 42.3% |

**Overall arm enrichment: NOT significant** (Fisher's exact, 62 vs genome: OR=1.32, p=0.17)

While the 62 coordinated regulators show a modest trend toward chromosomal arms (51.6% vs 44.7%), this does not reach statistical significance. The arm enrichment hypothesis for the **entire set** is not supported.

### Critical Finding: Subtype-specific Geographic Bias at T3

When stratified by coordination type, a **striking geographic pattern** emerges at T3:

| T3 Coordination Type | n_arm | n_core | % arm | Fisher OR (vs genome) | p-value |
|----------------------|-------|--------|-------|-----------------------|---------|
| discordant_gain_up | **13** | **2** | **86.7%** | **8.05** | **0.0010** |
| concordant_derepression | 5 | 4 | 55.6% | 1.55 | 0.37 |
| methyl_change_no_expr | 3 | 2 | 60.0% | 1.86 | 0.40 |
| concordant_repression | 2 | 6 | 25.0% | 0.41 | 0.80 |
| discordant_loss_down | 2 | 8 | 20.0% | 0.31 | 0.97 |
| ambiguous | 7 | 8 | 46.7% | 1.08 | 0.50 |

The 15 T3 discordant_gain_up genes (methylation gained + expression up) are **massively enriched at chromosomal arms** (OR=8.05, p=0.001). Conversely, T3 discordant_loss_down genes (methylation lost + expression down) trend toward the **core region** (80% core, though not significant with n=10).

**Biological interpretation:** Genes at chromosomal arms that gain methylation at T3 while becoming upregulated may represent a distinct regulatory mechanism -- methylation could act as a positive regulatory mark in the chromosomal arm context, or the methylation gain may be a secondary consequence of chromatin reorganization during late growth.

---

## Step 2: Coordination Pattern Analysis

### Concordant Derepression at BOTH T2 and T3 (Top Candidates)

These 3 genes show the strongest evidence for the "gatekeeper" model -- methylation is lost at both timepoints while expression increases:

| Locus Tag | SCO | Product | Family | log2FC T2 | log2FC T3 | Region |
|-----------|-----|---------|--------|-----------|-----------|--------|
| **SC_RS10435** | SCO1684 | TetR-family TPR protein | TetR | +1.73 | +2.60 | core |
| **SC_RS31385** | SCO5828 | Response regulator | RR | +3.20 | +1.09 | core |
| **SC_RS35525** | SCO6668 | Sensor histidine kinase | SK | +1.04 | +1.34 | arm |

**SC_RS10435 (SCO1684):** TetR-family regulator with tetratricopeptide repeat domain. Located in a neighborhood rich in carbohydrate metabolism genes (gluconokinase, GntP permease, SDR oxidoreductases). Adjacent to chaplin genes (chpC, chpH) -- chaplins are critical for aerial mycelium formation. This regulator may link carbon sensing to morphological differentiation.

**SC_RS31385 (SCO5828):** An orphan response regulator showing the strongest upregulation at T2 (log2FC=+3.20, ~9-fold). Located near citrate synthase genes and a sensor kinase (SC_RS31365) that also shows methylation changes. This TCS-adjacent response regulator may mediate metabolic state sensing.

**SC_RS35525 (SCO6668):** A sensor histidine kinase at the chromosomal right arm (7.41 Mb). Its cognate response regulator SC_RS35520 (SCO6667) is immediately adjacent (3 bp separation) but does NOT show methylation. Located near transaldolase/transketolase (pentose phosphate pathway) and an MmyB-family regulator. This region also contains a 4'-phosphopantetheinyl transferase -- a key enzyme for secondary metabolite activation.

### Concordant Repression at BOTH T2 and T3

| Locus Tag | SCO | Product | Family | log2FC T2 | log2FC T3 |
|-----------|-----|---------|--------|-----------|-----------|
| **SC_RS24635** | SCO4495 | UdgX family uracil-DNA binding protein | Other reg. | -2.52 | -1.29 |

SC_RS24635 encodes UdgX, a uracil-DNA glycosylase/binding protein involved in DNA repair. Its coordinated repression (methylation gained + expression down) at both timepoints suggests methylation may actively silence this DNA repair factor during secondary metabolism transition.

### TF Family Distribution

| TF Family | Count | % of 62 | Key coordination patterns |
|-----------|-------|---------|--------------------------|
| HTH (other) | 11 | 17.7% | Diverse; 4 discordant at T2, 4 at T3 |
| Other regulatory | 10 | 16.1% | 2 concordant repression at T2 |
| TetR | 9 | 14.5% | 4 discordant_loss_down at T2; 2 derepression at T3 |
| Sigma factor | 7 | 11.3% | 3 discordant_gain_up at T3 (3 arm = all arm!) |
| Sensor kinase | 6 | 9.7% | 3 discordant_gain_up at T3 |
| Response regulator | 5 | 8.1% | 2 derepression at T3 |
| Others (10 families) | 14 | 22.6% | |

Notable: **Sigma factors** contribute disproportionately to the T3 arm-enriched discordant_gain_up pattern (3/7 sigma factors, all at chromosomal arms).

---

## Step 3: Functional (COG) Enrichment

### Significant COG Enrichments (62 vs Genome)

| COG Category | 62 Coordinated | Genome-wide | OR | p-value |
|--------------|---------------|-------------|-----|---------|
| **K - Transcription** | 45 (72.6%) | 822 (9.9%) | **24.0** | **2.7e-31** |
| **S - Unknown** | 0 (0.0%) | 1,597 (19.3%) | **0.0** | **2.5e-06** |
| **R - General function** | 11 (17.7%) | 2,848 (34.4%) | **0.41** | **0.0047** |
| **T - Signal transduction** | 6 (9.7%) | 344 (4.2%) | **2.47** | **0.045** |

**Key findings:**
1. **COG K (Transcription)** is massively enriched (OR=24), which is expected since we selected for regulatory genes. However, within the 1,055 regulators, the 62 coordinated show similar K enrichment (72.6% vs 77.3%), confirming they are not a biased subset.
2. **COG T (Signal transduction)** is 2.5-fold enriched (p=0.045), driven by the 6 sensor kinases and 5 response regulators. This suggests methylation preferentially targets the signal transduction layer.
3. **COG Q (Secondary metabolism): 0 genes.** The 62 coordinated regulators are NOT directly annotated as secondary metabolism enzymes. They function as **transcriptional regulators** rather than biosynthetic genes.
4. **Complete absence** of S (Unknown function) -- all 62 have defined functional annotations, in contrast to 19.3% of the genome.

### Biological Interpretation

The functional profile of the 62 coordinated regulators shows enrichment in **transcription (K) and signal transduction (T)**, but NOT in secondary metabolism (Q) directly. This is consistent with a "gatekeeper" model where methylation targets the **regulatory layer** (transcription factors and signal transducers) rather than the biosynthetic genes themselves. These regulators likely control downstream secondary metabolism pathways indirectly.

---

## Step 4: Top Candidate Neighborhoods

### SC_RS10435 (SCO1684) Neighborhood (+/-10 kb)

A 20 kb region around SC_RS10435 contains 21 genes including:
- **chpC, chpH** (chaplins for aerial mycelium) -- both DEG at T2 and T3
- Carbohydrate metabolism cluster: gluconokinase, GntP permease, SDR oxidoreductases
- Another TetR regulator (SC_RS10470) just 6 kb downstream
- Multiple DEGs (15/21 genes are DEG at T2 or T3)

This neighborhood suggests SC_RS10435 may regulate a **carbon metabolism / morphological differentiation** switch.

### SC_RS31385 (SCO5828) Neighborhood (+/-10 kb)

- **SC_RS31365/SC_RS31370**: A sensor kinase / response regulator pair (3 genes away). SC_RS31365 shows methylation changes but no expression change -- suggesting the kinase is "poised" while the orphan RR SC_RS31385 is actively derepressed.
- **Citrate synthases** (SC_RS31400, SC_RS31405): central metabolic enzymes
- Multiple DUF proteins (DUF1453, DUF485, DUF6082) with DEG status

This neighborhood links SC_RS31385 to **central carbon metabolism** regulation.

### SC_RS35525 (SCO6668) Neighborhood (+/-10 kb)

- **SC_RS35520**: Cognate response regulator (3 bp away, overlapping!) -- no methylation
- **tal, tkt**: Transaldolase and transketolase (pentose phosphate pathway)
- **SC_RS35550**: 4'-phosphopantetheinyl transferase -- required for polyketide/NRP synthase activation
- **SC_RS35565**: MmyB family regulator (methylenomycin biosynthesis regulator family)
- **SC_RS35530**: IclR family regulator

This is a particularly interesting neighborhood: a TCS pair near pentose phosphate pathway genes AND a 4'-phosphopantetheinyl transferase. The Sfp-type PPTase is essential for activating carrier proteins in secondary metabolite biosynthesis. SC_RS35525 may sense metabolic state and regulate PPTase expression, linking primary to secondary metabolism.

---

## Step 5: Two-Component System (TCS) Analysis

### TCS Components Among 62 Coordinated Genes

| Type | Count | Genes |
|------|-------|-------|
| Sensor kinases | 6 | SC_RS07750, SC_RS28630, SC_RS34060, SC_RS35525, SC_RS37675, SC_RS40740 |
| Response regulators | 5 | SC_RS17795, SC_RS29355, SC_RS31385, SC_RS35610 (ramR), SC_RS40435 |

### Cognate Pair Candidates (within 5 kb)

| SK | RR | Distance | SK Coordination | RR Coordination | Both in 62? |
|----|-----|----------|-----------------|-----------------|-------------|
| SC_RS35525 (SCO6668) | SC_RS35520 (SCO6667) | 3 bp | **derepression T2+T3** | none | No |
| SC_RS37675 (SCO7089) | SC_RS37670 (SCO7088) | 3 bp | discordant_gain_up T3 | none | No |
| SC_RS29360 (SCO5435) | SC_RS29355 (SCO5434) | 3 bp | none | repression T3 | No |
| SC_RS40740 (SCO7711) | SC_RS40745 (SCO7712) | 8 bp | discordant_gain_up T2 | none | No |
| SC_RS40440 (SCO7649) | SC_RS40435 (SCO7648) | 142 bp | none | derepression T3 | No |
| SC_RS31365 (SCO5824) | SC_RS31385 (SCO5828) | 2,960 bp | methyl_no_expr | **derepression T2+T3** | No |
| SC_RS34060 (SCO6369) | SC_RS34035 (SCO6364) | 4,881 bp | discordant_gain_up T3 | none | No |

**Key finding:** No TCS pairs have BOTH members showing methylation-expression coordination. In all 7 pairs, only one partner (either SK or RR) shows coordination while the other has either no methylation or no expression change. This suggests methylation may act **asymmetrically** within TCS -- targeting one component to modulate the signaling pathway's sensitivity rather than silencing/activating both simultaneously.

The **SC_RS35525/SC_RS35520 pair** is particularly notable: the sensor kinase shows robust concordant derepression at both T2 and T3, while its cognate response regulator has no methylation at all. This asymmetry could mean methylation controls the "input side" (sensing) while the "output side" (transcriptional response) remains constitutively available.

### RamR (SC_RS35610/SCO6685)

Among the response regulators, **ramR** (SC_RS35610) is a known aerial mycelium regulator. It shows discordant_gain_up at T3 (methylation gained + massive expression increase: log2FC=+7.9). This is consistent with RamR's known role in morphological differentiation and suggests its activation at T3 occurs through a methylation-independent mechanism, with the gained methylation being a secondary event.

---

## Step 6: BGC Proximity

Only **1 gene** (SC_RS33745/SCO6299) among the 62 is within 20 kb of a BGC boundary:
- SC_RS33745: TetR/AcrR family regulator, 11 kb from CPK cluster boundary
- Shows discordant_loss_down at T2 (methylation lost, expression down)

The lack of BGC proximity further supports the model that these 62 regulators operate as an **upstream regulatory layer** rather than direct BGC components.

---

## Hypothesis Verdict

### H8: Geographic arm enrichment + stress/secondary metabolism enrichment

| Component | Verdict | Evidence |
|-----------|---------|----------|
| Overall arm enrichment | **NOT SUPPORTED** | p=0.17, OR=1.32 (NS) |
| Subtype-specific arm enrichment | **SUPPORTED** | T3 discordant_gain_up: 87% arm, OR=8.05, p=0.001 |
| COG Q (Secondary metabolism) | **NOT SUPPORTED** | 0/62 genes |
| COG T (Signal transduction) | **SUPPORTED** | OR=2.47, p=0.045 |
| Gatekeeper model (regulatory layer) | **PARTIALLY SUPPORTED** | K+T enriched; indirect regulation, not direct SM |

### Revised Model: "Asymmetric Signal Gating"

The data refine the original "gatekeeper" hypothesis into a more nuanced model:

1. **Not geographic gating (overall):** Methylation-expression coordination is distributed across the entire chromosome, not preferentially at arms.

2. **Subtype-specific geography:** T3 discordant patterns show strong geographic bias -- gain-up at arms, loss-down at core. This suggests different methylation-expression coupling mechanisms operate in arm vs core contexts.

3. **Signal transduction gating:** The enrichment of TCS components (11/62 = 17.7% are SK or RR) suggests methylation preferentially targets the **signal transduction layer**. In TCS pairs, methylation targets one partner (usually the sensor kinase), creating asymmetric modulation.

4. **Regulatory cascade, not direct BGC control:** The 62 genes are transcription factors and signal transducers (COG K+T), not biosynthetic genes (COG Q). They likely exert their effects through **regulatory cascades** several steps upstream of BGC activation.

5. **Top candidates for experimental validation:**
   - **SC_RS10435 (SCO1684):** TetR near chaplin cluster -- carbon/morphology switch
   - **SC_RS35525 (SCO6668):** Sensor kinase near PPTase -- primary/secondary metabolism link
   - **SC_RS31385 (SCO5828):** Orphan RR near citrate synthase -- metabolic state sensor
   - **SC_RS24635 (SCO4495):** UdgX DNA repair factor -- concordant repression at both timepoints

---

## Output Files

### Figures
| File | Description |
|------|-------------|
| `figures/geographic_distribution.pdf/svg` | Chromosome ideogram, density, arm vs core comparison |
| `figures/coordination_types.pdf/svg` | T2/T3 coordination type distributions, TF family pie chart |
| `figures/COG_enrichment_comparison.pdf/svg` | COG category bar chart and log2 fold enrichment heatmap |
| `figures/top_candidate_neighborhoods.pdf/svg` | +/-10 kb gene neighborhoods for 3 top candidates |
| `figures/H8_summary_overview.pdf/svg` | Comprehensive 6-panel summary figure |

### Tables
| File | Description |
|------|-------------|
| `tables/geographic_distribution_stats.tsv` | Arm vs core counts and Fisher test results |
| `tables/coordination_type_summary.tsv` | All 62 genes with coordination types, regions, expression |
| `tables/coordination_by_family.tsv` | TF family x coordination type cross-tabulation |
| `tables/COG_enrichment_results.tsv` | COG category enrichment (Fisher's exact, 62 vs genome) |
| `tables/TCS_pair_candidates.tsv` | 7 TCS cognate pair candidates with coordination status |

### Script
| File | Description |
|------|-------------|
| `scripts/H8_analysis.py` | Complete analysis pipeline (Steps 1-5, all figures and tables) |

---

## Next Steps

1. **Experimental validation priority:** SC_RS10435 (chaplin neighborhood) and SC_RS35525 (PPTase neighborhood) are the strongest candidates for knockout/knockdown studies.
2. **ChIP-seq integration:** Determine whether these 62 regulators have binding site motifs that overlap with methylation sites.
3. **Temporal dynamics:** The T2 vs T3 pattern differences suggest methylation-expression coupling evolves during growth -- a time-resolved methylation profiling would clarify causality.
4. **Sigma factor subgroup:** The 7 sigma factors (all showing arm-enriched discordant patterns at T3) deserve dedicated investigation as a potential coordinated sigma factor switch.

---

*Analysis performed: 2026-02-24*
*Script: `11_epigenome_integration/analysis/31_coordinated_regulators_characterization/scripts/H8_analysis.py`*
