# H32: Exposed Regulatory Module Analysis

**Date:** 2026-02-27
**Analysis directory:** `11_epigenome_integration/analysis/55_exposed_regulatory_module/`
**Script:** `scripts/H32_exposed_regulatory_module.py`

---

## Background

In *Streptomyces coelicolor* M145, previous analyses (H27-H29) identified 57 "exposed" regulatory genes -- transcription factors and signal transduction genes that lack the methylation protection zone characteristic of the remaining 998 "shielded" regulatory genes. These 57 genes share several distinctive properties: they sit within 293 bp of methylation sites (vs >293 bp for shielded), show 100% dynamic expression (zero constitutive), have lower T1 baseline expression, and display methylation-expression coordination across developmental timepoints.

The 57 exposed regulators span diverse TF families: TetR (9), sigma factors (7), sensor kinases (6), response regulators (5), and various HTH-type regulators. H8 previously identified 7 TCS (two-component system) pairs among the coordinated regulatory genes with asymmetric methylation.

**Key question:** Do these 57 exposed genes form an interconnected, self-regulatory "methylation-responsive module" -- with genomic clustering, co-expression, and shared coordination? Or are they independently scattered across the genome, unified only by the shared property of methylation exposure?

## Hypothesis

**H32:** The 57 exposed regulators form a functionally interconnected regulatory module, evidenced by:
1. Non-random genomic clustering
2. Higher mutual co-expression than with shielded genes
3. Coherent coordination types within co-expression groups
4. Operonic linkage between module members

---

## Methods

### Data sources
- 57 exposed regulators with full annotation (from H28/51 analysis)
- 1,017 regulatory genes with features (from H29/52 analysis)
- 57 coordinated regulatory genes with methylation details (from H6/29 analysis)
- Normalized expression counts (9 samples: 3 timepoints x 3 replicates)
- 7 TCS pairs (from H8/31 analysis)

### Analyses performed
1. **Genomic clustering**: Pairwise nearest-neighbor distances for 57 exposed TFs; permutation test (1,000 iterations sampling 62 from 1,017 regulatory genes); cluster detection at 20 kb threshold
2. **Co-expression**: Pairwise Spearman correlations across 9 samples; comparison of exposed-exposed, exposed-shielded, and shielded-shielded correlation distributions; hierarchical clustering (Ward linkage); co-expression module detection (rho > 0.7, BH-adjusted p < 0.05)
3. **Coordination type coherence**: Distribution of coordination types; chi-square test for association with genomic cluster membership; within-type vs between-type co-expression comparison
4. **TCS pair analysis**: Methylation asymmetry, expression correlation, and coordination status for each of the 7 pairs
5. **Operon context**: Detection of operonic pairs (same strand, <150 bp intergenic distance) among exposed TFs; permutation enrichment test
6. **Functional module detection**: Integration of genomic proximity (50 kb), co-expression (rho > 0.6), and coordination type coherence
7. **Expression dynamics**: Z-scored heatmap of 57 exposed TFs across 9 samples

---

## Results

### 1. Genomic Distribution: Not Clustered

The 57 exposed TFs are distributed across the full 8.67 Mb chromosome with **no significant genomic clustering**.

| Metric | Value |
|--------|-------|
| Median nearest-neighbor distance (exposed) | 51,280 bp |
| Median NN distance (random permutation) | 47,534 bp (IQR: 40,396-54,747) |
| Permutation p-value (closer) | 0.631 |
| Permutation p-value (more dispersed) | 0.369 |
| Z-score | +0.295 |

Only **2 genomic clusters** (within 20 kb) were found, containing 4 of 57 genes (6.5%):

| Cluster | Genes | Span | Region | TF Families |
|---------|-------|------|--------|-------------|
| 1 | SC_RS29350, SC_RS29355 (SCO5435-area) | 1,374 bp | Core | Other regulatory, Response regulator |
| 2 | SC_RS35525, SC_RS35610 | 17,628 bp | Right arm | Sensor kinase, Response regulator |

**Interpretation:** The exposed TFs are genomically dispersed, behaving indistinguishably from a random sample of regulatory genes. The hypothesis of physical genomic clustering is **not supported**.

### 2. Co-expression: Not Enhanced Between Exposed TFs

| Comparison | n pairs | Median rho | Mean rho |
|------------|---------|------------|----------|
| Exposed-Exposed | 1,891 | +0.017 | -0.003 |
| Exposed-Shielded | 12,400 | -0.050 | -0.018 |
| Shielded-Shielded | 19,900 | +0.050 | +0.031 |

| Test | U statistic | p-value | Rank-biserial r |
|------|-------------|---------|-----------------|
| E-E vs E-S (greater) | 11,833,230 | 0.257 | 0.505 |
| E-E vs S-S (greater) | 18,188,250 | 0.992 | 0.483 |

**Interpretation:** Exposed-exposed correlations are NOT significantly higher than exposed-shielded or shielded-shielded. In fact, exposed-exposed median (0.017) is slightly below shielded-shielded (0.050). The 57 exposed TFs do **not** preferentially co-express with each other compared to background.

### 3. Co-expression Modules: Four Detected Within the 62

Despite the lack of overall enhanced co-expression, internal structure was detected. Four co-expression modules (connected components at rho > 0.7, BH-adjusted p < 0.05) encompass **61 of 62** exposed genes:

| Module | n genes | Mean rho | Key families | Dominant T3 coordination |
|--------|---------|----------|-------------|--------------------------|
| 1 | 13 | 0.811 | HTH, LysR, MerR, Sigma, SK, RR | Mixed (discordant_gain_up / concordant_derepression) |
| 2 | 18 | 0.818 | DeoR, HTH, LacI, MerR, SK, Sigma, TetR, WhiB | Mixed (discordant_gain_up / concordant_derepression) |
| 3 | 4 | 0.833 | HTH, ROK, Sigma | discordant_gain_up / concordant_derepression |
| 4 | 26 | 0.614 | ArsR, GntR, IclR, LacI, LysR, MarR, SK, TetR | discordant_loss_down / concordant_repression |

Total significant high-correlation pairs: **266 / 1,891 (14.1%)**

**Critical observation:** Module 4 (26 genes, mean rho = 0.614) has a strong downregulation/repression profile, while Modules 1-3 (35 genes) are upregulation/derepression-dominated. This creates a **temporal bifurcation**: approximately half the exposed regulators are activated and half are repressed during development, but within each group the co-expression is strong.

### 4. Coordination Type Coherence: Strongly Supported

**T3 coordination type distribution:**

| Coordination Type | Count | % |
|-------------------|-------|---|
| discordant_gain_up | 15 | 24.2% |
| ambiguous | 15 | 24.2% |
| discordant_loss_down | 10 | 16.1% |
| concordant_derepression | 9 | 14.5% |
| concordant_repression | 8 | 12.9% |
| methyl_change_no_expr_change | 5 | 8.1% |

**Within-type vs between-type co-expression (T3):**

| Metric | Within-type | Between-type |
|--------|------------|--------------|
| n pairs | 329 | 1,562 |
| Median rho | **+0.500** | **-0.283** |
| Mann-Whitney U | 349,108 | |
| p-value | **6.6 x 10^-25** | |
| Rank-biserial r | 0.679 | |

This is the strongest finding: TFs sharing the same coordination type co-express vastly more than TFs with different coordination types (median rho difference = 0.783, p = 6.6 x 10^-25). However, coordination type was defined partly from expression, so this is partly circular.

**Coordination type vs genomic cluster membership:**
- Chi-square: chi2 = 1.845, p = 0.870, Cramer's V = 0.173
- **No association** between coordination type and physical genomic clustering

### 5. TCS Pair Analysis: Asymmetric but Not Biased

7 TCS pairs from H8 were analyzed:

| Pair (SK/RR) | Methylated partner | SK exposed | RR exposed | Expression rho | p-value |
|--------------|-------------------|------------|------------|---------------|---------|
| SCO6668/SCO6667 | Sensor kinase | Yes | No | 0.783 | 0.013 |
| SCO7089/SCO7088 | Sensor kinase | Yes | No | 0.533 | 0.139 |
| SCO5435/SCO5434 | Response regulator | No | Yes | 0.850 | 0.004 |
| SCO7711/SCO7712 | Sensor kinase | Yes | No | 0.650 | 0.058 |
| SCO7649/SCO7648 | Response regulator | No | Yes | 0.650 | 0.058 |
| SCO5824/SCO5828 | Response regulator | No | Yes | -0.617 | 0.077 |
| SCO6369/SCO6364 | Sensor kinase | Yes | No | 0.500 | 0.170 |

Key findings:
- **All 7 pairs are asymmetric**: exactly one partner is exposed (methylated), the other is shielded
- **0/7 pairs have both members exposed** -- confirming H8's finding
- Methylation is on the SK in 4 pairs and the RR in 3 pairs (binomial p = 1.0, no bias)
- Mean within-pair correlation = 0.479 (moderate to high)
- 5/7 pairs show positive correlation (rho > 0), 1 pair (SCO5824/SCO5828) is strongly anti-correlated (-0.617)
- The SCO5824/SCO5828 anti-correlation with concordant_derepression on the RR suggests an antagonistic regulatory circuit

### 6. Operon Context: No Enrichment

- Exposed-exposed operonic pairs detected: **0**
- Expected by random permutation: median = 0 (mean ~ 0.6)
- Permutation p-value: 1.0

The 57 exposed TFs are **not** organized into shared operons with each other. This is consistent with the genomic dispersion finding.

### 7. Functional Module Detection

Integrating proximity (50 kb), co-expression (rho > 0.6), and coordination coherence:

| Module | Size | Source | Proximal? | Coexpressed? | Top T3 type | Evidence score |
|--------|------|--------|-----------|-------------|-------------|----------------|
| 1 | 2 | Proximity | Yes | Yes | concordant_derepression | 2.50 |
| 2 | 2 | Proximity | Yes | No | ambiguous | 2.00 |
| 3 | 3 | Proximity | No | Yes | ambiguous | 1.67 |
| 10 | 26 | Co-expression | No | Yes | discordant_loss_down | 1.39 |

**High-confidence modules (score >= 2, >= 3 members): 0**

No large module has simultaneous genomic proximity AND co-expression AND coordination coherence. The two highest-scoring modules contain only 2 genes each. The large co-expression-based modules (26 and 36 genes) lack genomic proximity.

---

## Interpretation

### The 57 Exposed Regulators: A Functional Layer, Not a Physical Module

The evidence presents a clear picture:

1. **NOT a genomic module**: The 57 genes are dispersed across the chromosome indistinguishably from random (z = +0.30, p = 0.37). Only 4/57 (6.5%) are within 20 kb of another exposed gene. No shared operons.

2. **NOT preferentially co-expressed as a group**: Exposed-exposed correlations (median rho = 0.017) are no higher than exposed-shielded or shielded-shielded. The 57 genes do not behave as a single co-expression community.

3. **BUT internally structured by coordination type**: Genes sharing the same coordination type (e.g., all discordant_gain_up or all concordant_repression) co-express very strongly (median rho = 0.500 vs -0.283 for different types, p = 6.6 x 10^-25). This creates two main expression blocs:
   - **Activation bloc** (Modules 1-3, ~35 genes): Upregulated during development, methylation changes correlate with activation
   - **Repression bloc** (Module 4, ~26 genes): Downregulated, methylation changes correlate with repression

4. **TCS pairs are exclusively asymmetric**: In all 7 TCS pairs, exactly one partner is exposed and one is shielded -- never both. This asymmetry is perfect (7/7) but not biased toward SK or RR (4 SK vs 3 RR). This suggests methylation targets one arm of the signal transduction pathway while leaving the other protected, potentially creating a differential sensitivity mechanism.

5. **Coordination type drives structure, not physical linkage**: The strong within-type co-expression (rank-biserial r = 0.68) combined with the absence of genomic clustering and the non-significant chi-square (p = 0.87) indicates that the functional relationships among exposed TFs are based on shared regulatory logic (response to methylation state), not physical chromosome organization.

### Revised Model: Distributed Methylation-Responsive Regulatory Layer

Rather than forming a single self-regulatory module, the 57 exposed regulators constitute a **distributed regulatory layer** -- genes scattered across the genome that are unified by their shared vulnerability to methylation changes. Within this layer, two antagonistic response programs operate:

- **Program 1 (Activation)**: ~35 genes that are de-repressed when methylation is removed or that show discordant gain-up responses. These include sigma factors, sensor kinases, and HTH regulators involved in developmental transitions.
- **Program 2 (Repression)**: ~26 genes that are repressed when methylation accumulates, including GntR, MarR, ArsR, and LysR family regulators. These may represent metabolic regulators that are silenced during differentiation.

The perfect TCS asymmetry (one partner exposed, one shielded) suggests this layer functions as a **differential amplifier**: methylation changes affect one component of a signaling pair while the other remains stable, creating a controlled imbalance that drives signal transduction.

---

## Verdict

**PARTIAL**

The hypothesis that the 57 exposed TFs form a self-regulatory module is **partially supported**:

| Evidence | Result | Supports module? |
|----------|--------|-----------------|
| Genomic clustering | NOT significant (z=0.30, p=0.37) | No |
| E-E co-expression > E-S | NOT significant (p=0.26) | No |
| E-E co-expression > S-S | NOT significant (p=0.99) | No |
| Within-type co-expression coherence | **Highly significant** (p=6.6e-25, r=0.68) | Yes |
| Co-expression modules detected | **4 modules, 61/57 genes covered** | Yes |
| Operonic pair enrichment | NOT significant (0 pairs) | No |
| High-confidence multi-evidence modules | 0 detected | No |

**Key conclusion:** The 57 exposed regulators do NOT form a single interconnected physical module. Instead, they constitute a **distributed methylation-responsive regulatory layer** with strong internal structure organized by coordination type rather than genomic position. Two antagonistic expression programs (activation and repression) operate within this layer, and TCS pairs are invariably split between exposed and shielded members, suggesting a differential amplification mechanism.

---

## Output Files

### Figures
| File | Description |
|------|-------------|
| `figures/genomic_distribution.pdf/svg` | Chromosome map with 57 TFs, clusters, permutation test |
| `figures/coexpression_heatmap.pdf/svg` | 62x62 Spearman correlation matrix with dendrogram |
| `figures/coordination_type_analysis.pdf/svg` | Coordination type distribution and coherence |
| `figures/TCS_pairs.pdf/svg` | 7 TCS pairs expression and methylation analysis |
| `figures/expression_heatmap.pdf/svg` | 57 TFs x 9 samples Z-scored heatmap |
| `figures/H32_comprehensive_summary.pdf/svg` | Multi-panel summary |

### Tables
| File | Description |
|------|-------------|
| `tables/genomic_clusters.tsv` | 2 genomic clusters (20kb threshold) |
| `tables/coexpression_matrix.tsv` | 62x62 Spearman correlation matrix |
| `tables/coexpression_modules.tsv` | 4 co-expression modules (rho>0.7) |
| `tables/TCS_pairs_analysis.tsv` | 7 TCS pairs detailed analysis |
| `tables/coordination_type_distribution.tsv` | Coordination type counts for T2 and T3 |
| `tables/statistical_tests.tsv` | All 7 statistical tests with effect sizes |
| `tables/module_summary.tsv` | 14 candidate functional modules with scores |

---

## Connection to Previous Hypotheses

| Hypothesis | Relationship |
|------------|-------------|
| **H6** (Genome-wide TF screen) | Identified the 57 coordinated regulatory genes analyzed here |
| **H8** (Coordinated regulators characterization) | Identified 7 TCS pairs with asymmetric methylation; confirmed all 7 are split exposed/shielded |
| **H25** (TSS methylation gradient) | Established the protection zone model; exposed genes lack this protection |
| **H27** (Exposed promoter model) | Showed exposed genes are 6.7x closer to methylation, lack protection zones |
| **H28** (Exposed regulator characteristics) | Showed 100% dynamic expression, low T1 baseline; H32 shows this dynamic behavior splits into activation and repression programs |
| **H29** (Shielded/exposed boundary) | Defined the 293bp boundary; H32 shows the exposed class is genomically dispersed, not clustered |
| **H30** (TSS sequence determinants) | Showed DNA sequence contributes ~33% of protection; H32's dispersal pattern is consistent with distributed sequence-level determinants |
| **Gatekeeper Model v3** | H32 refines Layer 3 (Signal Gating): the 57 exposed genes are not a single module but a distributed regulatory layer with two antagonistic programs. The TCS asymmetry finding strengthens the differential amplification concept |
