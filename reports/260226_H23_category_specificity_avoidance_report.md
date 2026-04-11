# H23: Gene Category Specificity of Methylation Avoidance

**Date**: 2026-02-26
**Analysis directory**: `11_epigenome_integration/analysis/46_category_specificity_avoidance/`
**Script**: `scripts/H23_category_specificity_avoidance.py`

---

## Background

H15 and H20 established that methylation sites (GCCGGC 4mC, AAGCCCG 6mA) are depleted near regulatory genes when analyzed from a **site-centric** perspective (i.e., fewer methylation sites target regulatory gene neighborhoods than expected). H20 confirmed this avoidance is geographically robust (present in both core and arm regions).

**H23 asks the converse question**: is this avoidance specific to regulatory genes, or does it extend to other functionally important gene categories (translation, DNA replication/repair, energy metabolism, etc.)?

- If avoidance is universal across essential categories: purifying selection against R-M damage to any critical gene
- If avoidance is regulatory-specific: something unique about regulatory DNA (e.g., TF binding site sensitivity) requires protection
- If intermediate (regulatory strongest, others moderate): graded avoidance with regulatory peak

## Methodology

### Gene Classification

8,275 genes were classified into 9 functional categories based on product annotation keywords and curated lists:

| Category | n genes | Classification basis |
|----------|---------|---------------------|
| **Regulatory** | 1,055 | Curated list from H6 genome-wide TF screen |
| **Other** | 4,951 | All remaining genes |
| **Hypothetical** | 1,077 | "hypothetical protein" in product |
| **Transport** | 693 | transporter, permease, MFS, ABC, efflux |
| **Translation** | 148 | ribosomal protein, tRNA, translation factor |
| **DNA replication/repair** | 115 | DNA polymerase, helicase, gyrase, recombinase |
| **Energy/respiration** | 109 | cytochrome, ATP synthase, NADH |
| **Secondary metabolism** | 81 | BGC-assigned genes (from gene_master_with_BGC) |
| **Cell division** | 46 | FtsZ, division, septum, cell wall, peptidoglycan |

### Methylation Proximity

For each gene, determined whether any methylation site falls within 2 kb of the gene boundaries (start - 2 kb to end + 2 kb). Three methylation datasets at T1:
- **GCCGGC 4mC**: 1,289 sites (N4-cytosine methylation at GCCGGC palindrome)
- **AAGCCCG 6mA**: 260 sites (N6-adenine methylation)
- **All 4mC**: 1,987 sites (all 4mC modifications from census)

### Statistical Tests
- Fisher's exact test for each category vs genome-wide rate
- Bonferroni correction across all category x motif tests (27 tests)
- Jonckheere-Terpstra trend test for expression quintile analysis
- Core/arm stratification (core: 1.5-7.17 Mb; arms: flanking regions)

## Results

### 1. Category-Level Enrichment/Depletion

**Genome-wide proximity rates**: GCCGGC 46.9%, AAGCCCG 12.0%, All 4mC 59.7%

| Category | GCCGGC fold | GCCGGC p_bonf | AAGCCCG fold | AAGCCCG p_bonf | All 4mC fold | All 4mC p_bonf |
|----------|-------------|---------------|--------------|----------------|--------------|----------------|
| **Hypothetical** | **0.86** | **1.2e-04** | 1.00 | 1.0 | **0.88** | **1.8e-05** |
| **Regulatory** | 1.01 | 1.0 | 0.91 | 1.0 | 0.99 | 1.0 |
| Translation | 1.04 | 1.0 | 1.12 | 1.0 | **1.23** | **0.013** |
| DNA replication/repair | 1.00 | 1.0 | 1.30 | 1.0 | 1.00 | 1.0 |
| Cell division | 1.16 | 1.0 | 1.08 | 1.0 | 1.24 | 1.0 |
| Energy/respiration | 1.29 | 0.13 | 1.14 | 1.0 | 1.17 | 1.0 |
| Transport | 1.08 | 1.0 | 1.01 | 1.0 | 1.06 | 1.0 |
| **Secondary metabolism** | **1.66** | **4.2e-07** | 0.92 | 1.0 | **1.43** | **2.6e-05** |
| Other | 1.00 | 1.0 | 1.01 | 1.0 | 1.00 | 1.0 |

**Key finding**: From the **gene-centric** perspective (what fraction of genes in each category have nearby methylation?), **Regulatory genes show NO avoidance** (fold ~1.0, all p > 0.24). This contrasts sharply with the H15 **site-centric** result (fold = 0.61, p = 3.9e-06).

Only two categories show significant effects:
1. **Hypothetical genes**: Depleted for GCCGGC (fold=0.86) and All 4mC (fold=0.88) -- fewer hypothetical genes near methylation sites than expected
2. **Secondary metabolism genes**: Strongly enriched for GCCGGC (fold=1.66, OR=4.01) and All 4mC (fold=1.43, OR=3.92) -- BGC genes are methylation-rich
3. **Translation genes**: Modestly enriched for All 4mC (fold=1.23, OR=1.91)

### 2. Reconciling H15 (Site-Centric) vs H23 (Gene-Centric)

The apparent contradiction between H15 and H23 is resolved by understanding the two perspectives:

| Approach | Question | Regulatory result | Interpretation |
|----------|----------|-------------------|----------------|
| **H15 site-centric** | "Given a methylation site, is the nearest gene regulatory?" | fold=0.61, p=3.9e-06 | Fewer sites **target** regulatory neighborhoods |
| **H23 gene-centric** | "Given a regulatory gene, is there a nearby methylation site?" | fold=1.01, NS | The **fraction** of regulatory genes with nearby sites is genome-average |

**Explanation**: With 1,289 GCCGGC sites across 8.67 Mb, 46.9% of all genes are within 2 kb of at least one site. At this high density, the binary proximity metric saturates -- nearly half of all genes in every category are "proximal" regardless. The H15 site-centric approach, which counts the number of sites targeting each category, has greater resolution because it captures **methylation density** rather than binary proximity.

This means: regulatory genes are not methylation-free zones, but they are **methylation-sparse** zones. They have fewer sites per gene than expected, even though most of them (~47%) have at least one site somewhere in a 2 kb window.

### 3. Expression Quintile Analysis

Genes divided into quintiles by baseMean (average expression across conditions):

| Quintile | baseMean range | GCCGGC fold | GCCGGC p_bonf | All 4mC fold | All 4mC p_bonf |
|----------|---------------|-------------|---------------|--------------|----------------|
| Q1 (lowest) | 4.7 - 40.8 | **0.82** | **1.1e-12** | **0.83** | **2.4e-19** |
| Q2 | 40.8 - 94.4 | **0.91** | **7.3e-04** | **0.91** | **6.7e-06** |
| Q3 | 94.4 - 229.5 | 0.98 | 1.0 | 0.99 | 1.0 |
| Q4 | 229.6 - 674.2 | **1.08** | **0.015** | **1.06** | **5.4e-03** |
| Q5 (highest) | 674.4 - 1.26M | **1.22** | **2.9e-19** | **1.20** | **1.7e-30** |

**Jonckheere-Terpstra trend tests**:
- GCCGGC: z = +10.0, p ~ 0 (highly significant increasing trend)
- AAGCCCG: z = +2.1, p = 0.037 (nominally significant, Bonf NS)
- All 4mC: z = +12.0, p ~ 0 (highly significant increasing trend)

**Critical result**: The trend is **the opposite of avoidance**. Highly expressed genes (Q5) are 1.2x MORE likely to have nearby methylation sites, while lowly expressed genes (Q1) are 0.82x LESS likely. Methylation proximity is positively correlated with expression level.

### 4. Ribosomal Protein Deep Dive

65 ribosomal protein genes (40 x 50S, 24 x 30S) as a test case for highly conserved, essential genes:

| Group | n | GCCGGC fold | GCCGGC p | All 4mC fold | All 4mC p_bonf |
|-------|---|-------------|----------|--------------|----------------|
| All ribosomal | 65 | 0.85 | 0.32 | 1.29 | 0.058 |
| 50S subunit | 40 | 0.85 | 0.43 | 1.38 | **0.038** |
| 30S subunit | 24 | 0.89 | 0.69 | 1.19 | 1.0 |
| Regulatory | 1,055 | 1.01 | 0.74 | 0.99 | 1.0 |

**No significant avoidance for ribosomal proteins for GCCGGC**. The 50S subunit shows weak **enrichment** for All 4mC (fold=1.38, p=0.038 after Bonferroni). This is consistent with the expression quintile finding: ribosomal protein genes are among the most highly expressed genes and thus fall in Q5 (enriched).

### 5. Core/Arm Stratification

Key categories in GCCGGC analysis:

| Category | Core fold | Core p | Arm fold | Arm p |
|----------|-----------|--------|----------|-------|
| Regulatory | 1.00 | 0.93 | 1.06 | 0.40 |
| Translation | **0.81** | **0.012** | **1.94** | **0.026** |
| Hypothetical | 0.97 | 0.41 | **0.80** | **0.002** |
| Energy/respiration | **1.31** | **0.001** | 1.06 | 0.86 |

Translation genes show opposite patterns by region: depleted in core (fold=0.81) but enriched in arms (fold=1.94). This likely reflects the clustering of ribosomal protein operons (largely in core) versus scattered translation-related genes in arms. Hypothetical gene depletion is arm-specific (fold=0.80, p=0.002), consistent with arm regions having both more hypothetical genes and distinct methylation patterns.

## Interpretation

### The H15/H23 Dichotomy Resolves the Mechanism

The combination of H15 (site-centric avoidance) and H23 (gene-centric non-avoidance) reveals a nuanced picture:

1. **Methylation sites actively avoid regulatory gene neighborhoods** (H15: fewer sites per kb near regulatory genes)
2. **But regulatory genes are not in methylation-free zones** (H23: similar fraction of regulatory genes have at least one nearby site)
3. **Methylation density, not binary proximity, drives the signal**: The avoidance is about reduced **density** of modifications near regulatory genes, not complete exclusion

### Expression-Methylation Positive Correlation

The strongest signal in H23 is the **positive** correlation between expression level and methylation proximity (JT z = +10-12). This is the opposite of what a simple "methylation damages important genes" model would predict. Possible explanations:

1. **Open chromatin accessibility**: Highly expressed genes are in more accessible chromatin, making them more available to methyltransferases
2. **Gene density confound**: Highly expressed genes (e.g., ribosomal operons) tend to cluster in gene-dense core regions where methylation sites are also dense
3. **Positive regulatory role**: Methylation near highly expressed genes may serve a regulatory rather than damaging function

### Verdict

**H23 rejects the "universal avoidance of important genes" model.** The data support a more complex picture:

- **No gene category shows the strong avoidance seen in H15's site-centric analysis** when tested gene-centrically
- **Only Hypothetical genes** show significant depletion (modest, fold ~0.86-0.88)
- **Secondary metabolism (BGC) genes** are significantly ENRICHED for methylation (fold 1.43-1.66)
- **Highly expressed genes** are MORE likely to be methylation-proximal, not less
- The H15 regulatory avoidance reflects reduced **methylation density** per gene, not reduced **probability of any nearby site**

**This refines the Gatekeeper Model**: methylation avoidance operates at the level of site density modulation near regulatory sequences, not through wholesale exclusion of methylation from important gene neighborhoods. The mechanism is more consistent with selective pressure on **specific regulatory elements** (promoters, TF binding sites) rather than broad genomic neighborhoods.

## Output Files

### Figures
| File | Description |
|------|-------------|
| `figures/category_avoidance_heatmap.pdf/svg` | Category x motif fold enrichment heatmap |
| `figures/forest_plot_categories.pdf/svg` | Forest plot with OR and 95% CI for each category |
| `figures/expression_quintile_avoidance.pdf/svg` | Expression level vs methylation proximity |
| `figures/H23_comprehensive_summary.pdf/svg` | Multi-panel comprehensive summary |

### Tables
| File | Description |
|------|-------------|
| `tables/category_enrichment.tsv` | Per-category enrichment statistics (9 categories x 3 motifs) |
| `tables/expression_quintile_analysis.tsv` | Expression quintile analysis (5 quintiles x 3 motifs) |
| `tables/jt_trend_tests.tsv` | Jonckheere-Terpstra trend test results |
| `tables/ribosomal_protein_analysis.tsv` | Ribosomal protein deep dive |
| `tables/stratified_by_region.tsv` | Core/arm stratified results |

### Script
| File | Description |
|------|-------------|
| `scripts/H23_category_specificity_avoidance.py` | Complete analysis pipeline |

## Relationship to Other Hypotheses

| Hypothesis | Finding | H23 Update |
|-----------|---------|------------|
| H15 | Regulatory/TF depleted (site-centric fold=0.61) | Gene-centric fold=1.01 -- avoidance is density-based, not exclusion-based |
| H20 | H15 avoidance robust to core/arm stratification | H23 confirms no gene-centric avoidance in either region |
| H18 | GCCGGC dose-response is threshold (0 vs >=1), not graded | Consistent: binary proximity metric captures threshold; density metric captures graded signal |
| H17 | GCCGGC-proximal genes show suppression vs AAGCCCG activation | Expression quintile shows highly expressed genes are MORE proximal |

## Conclusions

1. **Regulatory-specific avoidance is a density phenomenon**: Regulatory genes have fewer methylation sites per kb of flanking DNA, but similar probability of having at least one site nearby (46-47%)
2. **No universal avoidance of essential genes**: Translation, DNA replication/repair, energy metabolism, and cell division genes show no significant avoidance
3. **Expression level positively correlates with methylation proximity**: Highly expressed genes (Q5) have 1.2x higher proximity than lowly expressed genes (Q1), driven by GCCGGC and All 4mC
4. **Secondary metabolism (BGC) genes are methylation hotspots**: 1.4-1.7x enrichment, consistent with BGCs being targets of epigenetic regulation
5. **Hypothetical genes are the only depleted category** (gene-centrically), possibly because they tend to be shorter, less conserved, and located in methylation-sparse regions

---

*Analysis performed with Python 3, scipy, numpy, pandas, matplotlib.*
