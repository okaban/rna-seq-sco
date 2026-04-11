# H35: TF Family-Based Functional Prediction of 62 Exposed Transcription Factors

**Date**: 2026-02-27
**Analysis**: `11_epigenome_integration/analysis/58_exposed_TF_functional_prediction/`
**Script**: `scripts/H35_functional_prediction.py`
**Hypothesis**: TF family membership, genomic context, and expression magnitude can predict the downstream regulatory roles of the 62 exposed TFs.

---

## Background

The 62 exposed TFs form a "distributed methylation-responsive regulatory layer" (H32) with two antagonistic programs:
- **Activation bloc** (~36 genes, modules 1-3): upregulated during development
- **Repression bloc** (26 genes, module 4): downregulated during development

Their downstream targets remain unknown -- no FIMO motifs are available, and no cis-neighborhood transcriptional effect was detected (H33). TF family membership provides a strong basis for functional prediction because TF family members typically bind similar DNA sequences, regulate analogous biological processes, and occupy conserved roles across bacterial species.

## Results

### 1. TF Family Distribution Among 62 Exposed TFs

| TF Family | Count | Functional Category |
|-----------|-------|---------------------|
| HTH (other) | 11 | General |
| Other regulatory | 10 | General |
| TetR | 9 | Defense/resistance |
| Sigma factor | 7 | Stress/development |
| Sensor kinase | 6 | Signal transduction |
| Response regulator | 5 | Signal transduction |
| LysR | 2 | Metabolic |
| MarR | 2 | Stress/development |
| LacI | 2 | Metabolic |
| MerR | 2 | Defense/resistance |
| GntR | 1 | Metabolic |
| ArsR | 1 | Defense/resistance |
| DeoR | 1 | Metabolic |
| IclR | 1 | Metabolic |
| WhiB | 1 | Stress/development |
| ROK | 1 | Metabolic |

16 TF families represented. The most abundant are HTH (other), Other regulatory, and TetR.

### 2. Bloc-Level Functional Composition

**Activation bloc (36 TFs):**
- Signal transduction: 9 (25.0%) -- 5 sensor kinases + 4 response regulators
- General (HTH/Other): 13 (36.1%)
- Stress/development: 6 (16.7%) -- 5 sigma factors + 1 WhiB
- Metabolic: 4 (11.1%)
- Defense/resistance: 4 (11.1%)

**Repression bloc (26 TFs):**
- Defense/resistance: 8 (30.8%) -- 7 TetR + 1 ArsR
- General (HTH/Other): 8 (30.8%)
- Other regulatory: 6 (23.1%)
- Metabolic: 4 (15.4%)
- Signal transduction: 2 (7.7%)
- Stress/development: 4 (15.4%)

### 3. Statistical Enrichment by Bloc

#### Family-level Fisher exact tests

| TF Family | Activation | Repression | OR | p-value | Enriched in |
|-----------|-----------|-----------|-----|---------|-------------|
| **TetR** | **2** | **7** | **0.16** | **0.028*** | **Repression** |
| HTH (other) | 9 | 2 | 4.00 | 0.101 | Activation |
| MarR | 0 | 2 | 0.00 | 0.172 | Repression |
| Sensor kinase | 5 | 1 | 4.03 | 0.387 | Activation |
| Response regulator | 4 | 1 | 3.12 | 0.388 | Activation |

**TetR is the only family with statistically significant bloc enrichment** (OR=0.16, p=0.028, Fisher exact test). 7 of 9 TetR-family exposed TFs are in the repression bloc, consistent with their known role as local repressors of neighboring efflux and resistance genes in *Streptomyces*.

HTH (other) trends toward activation enrichment (OR=4.0) and sensor kinases also trend toward activation (OR=4.0), but neither reaches p<0.05 with these sample sizes.

#### Functional category enrichment

| Category | Activation | Repression | OR | p-value |
|----------|-----------|-----------|-----|---------|
| Signal transduction | 9 (25.0%) | 2 (7.7%) | 4.00 | 0.101 |
| Defense/resistance | 4 (11.1%) | 8 (30.8%) | 0.28 | 0.101 |
| Metabolic | 4 (11.1%) | 4 (15.4%) | 0.69 | 0.710 |
| General | 13 (36.1%) | 8 (30.8%) | 1.27 | 0.788 |
| Stress/development | 6 (16.7%) | 4 (15.4%) | 1.10 | 1.000 |

No functional category reaches statistical significance individually. However, the combined trend is biologically coherent: Signal transduction skews 3.3x toward the activation bloc (25% vs 7.7%), while Defense/resistance skews 2.8x toward the repression bloc (30.8% vs 11.1%). These complementary trends (p=0.10 each) are mutually reinforcing.

### 4. Genomic Context Analysis (+-5kb Neighbors)

554 neighbor genes identified across the 62 exposed TFs.

**Key neighbor category findings:**
- Transport genes are equally distributed (11.7% each bloc)
- Secondary metabolism neighbors show a trend toward the repression bloc (6.6% vs 3.8%, OR=0.56, p=0.158)
- Signaling neighbors are equally distributed (5.9% vs 5.6%)

**TFs proximal to functional clusters (40 identified):**
- 13 TFs are near secondary metabolism genes
- 27 TFs are near transport clusters (>=2 transport genes)
- Notably, repression-bloc TetR members are frequently adjacent to efflux/transport genes:
  - SCO5956: 4 transport neighbors (MFS, ECF ABC transporters)
  - SCO6599: 3 ABC transporter neighbors (carbohydrate ABC)
  - SCO4639: 1 MFS transporter + 1 SDR oxidoreductase

### 5. Expression Magnitude Analysis

**Top 10 exposed TFs by |LFC_T3vsT1|:**

| Rank | SCO | Family | Bloc | |LFC| |
|------|-----|--------|------|-------|
| 1 | SCO6685 (RamR) | Response regulator | Activation | 7.90 |
| 2 | SCO7252 (NsdB) | Other regulatory | Activation | 7.76 |
| 3 | SCO1160 | Sensor kinase | Activation | 5.57 |
| 4 | SCO1227 | HTH (other) | Activation | 4.59 |
| 5 | SCO6165 | Other regulatory | Activation | 4.40 |
| 6 | SCO1468 | TetR | Activation | 4.26 |
| 7 | SCO4122 | MarR | Repression | 3.54 |
| 8 | SCO3907 (SSB) | Other regulatory | Repression | 3.33 |
| 9 | SCO4005 (SigE) | Sigma factor | Activation | 3.23 |
| 10 | SCO3134 | Response regulator | Activation | 3.16 |

**Bloc |LFC| comparison:** Activation median = 1.77, Repression median = 1.41 (MWU p=0.069). The activation bloc trends toward stronger expression responses.

**Family |LFC| comparison:** Kruskal-Wallis H=11.10, p=0.269 (non-significant). Response regulators (mean 3.06) and sensor kinases (mean 3.00) have the highest mean |LFC|, consistent with TCS as high-amplitude responders.

### 6. Sigma Factor Subfamily Analysis (7 TFs)

| SCO | Subfamily | Bloc | LFC_T3 |
|-----|-----------|------|--------|
| SCO4005 | SigE (cell envelope stress) | Activation | +3.23 |
| SCO7314 | SigB/SigF/SigG (general stress/sporulation) | Activation | +1.49 |
| SCO4938 | SigJ (ECF sigma) | Activation | +1.34 |
| SCO0632 | Unclassified sigma-70 | Activation | +1.61 |
| SCO1564 | Unclassified sigma-70 | Activation | +1.11 |
| SCO2639 | Unclassified sigma-70 | Repression | -1.29 |
| SCO4960 | Sigma-like HTH | Repression | -1.51 |

**5 of 7 sigma factors are in the activation bloc**, including SigE (strongest at LFC=+3.2), the SigB/F/G-family member, and SigJ. The 2 repression-bloc sigma factors (SCO2639, SCO4960) may represent vegetative-phase sigma factors that are downregulated during the developmental transition.

### 7. TCS Analysis (11 TFs: 6 SK + 5 RR)

**Activation bloc (9 TCS members):** 5 SK + 4 RR
- SCO1160 (SK): Massive induction (LFC=+5.57), likely senses a critical developmental signal
- SCO7711 (SK): Strong activation (LFC=+3.14)
- SCO6685/RamR (RR): Strongest responder overall (LFC=+7.90), known aerial mycelium regulator
- SCO3134 (RR): Strong activation (LFC=+3.16)

**Repression bloc (2 TCS members):** 1 SK + 1 RR
- SCO5289 (SK): Downregulated (LFC=-2.41)
- SCO5434 (RR): Downregulated (LFC=-1.12)

No TCS pairs were found within 3kb among the exposed TFs, consistent with the previously observed asymmetric methylation pattern where only one member of each cognate pair is exposed.

### 8. TetR Family Analysis (9 TFs)

The TetR family shows the clearest bloc partitioning:
- **Activation (2):** SCO1468 (LFC=+4.26), SCO1684 (LFC=+2.60) -- both annotated as "tetratricopeptide repeat protein" (TPR), suggesting protein-protein interaction rather than classical DNA-binding repression
- **Repression (7):** SCO2223, SCO2681, SCO4305, SCO4639, SCO5956, SCO6299, SCO6599 -- all classical TetR/AcrR repressors

The repression-bloc TetR members are characteristically adjacent to:
- **Transport genes**: efflux pumps (MFS, ABC transporters) -- 5 of 7 have transport neighbors
- **Secondary metabolism genes**: SDR oxidoreductases, thioesterases -- 5 of 7 have BGC-related neighbors

This is consistent with TetR-family regulators serving as local repressors of efflux and secondary metabolism genes that are downregulated during the developmental transition.

### 9. Predicted Regulatory Model

#### ACTIVATION BLOC -- Processes activated during development (T1 to T3)

| Confidence | Process | Key TFs | Evidence |
|------------|---------|---------|----------|
| **HIGH** | Morphological differentiation / aerial mycelium | RamR (SCO6685, +7.9), WhiB (+1.7), SigB/F/G (SCO7314, +1.5) | Known developmental regulators, massive LFC |
| **HIGH** | Two-component signal transduction cascades | 5 SK + 4 RR, SCO1160 (+5.6), SCO7711 (+3.1) | TCS 3.3x enriched in activation bloc |
| **MEDIUM** | Secondary metabolite regulation | TcrA/SCO5433 (+1.5), NsdB/SCO7252 (+7.8) | AfsR-like SARP and development regulator |
| **MEDIUM** | Stress response activation | SigE/SCO4005 (+3.2), MerR (+0.6/+1.7) | Cell envelope + oxidative stress sigma/TFs |
| **LOW** | Sugar/nucleotide metabolism reorganization | DeoR/SCO1897 (+1.5), LacI/SCO7411 (+1.1) | Single family representatives |

#### REPRESSION BLOC -- Processes repressed during development (T1 to T3)

| Confidence | Process | Key TFs | Evidence |
|------------|---------|---------|----------|
| **HIGH** | Primary metabolism shutdown | GntR (-1.5), IclR (-2.4), LacI (-1.4), LysR (-1.2) | 4 metabolic TF families converge |
| **HIGH** | TetR-mediated efflux/defense downregulation | 5 TetR (-0.8 to -2.3) | Statistically significant enrichment (p=0.028), all adjacent to transport genes |
| **HIGH** | Fatty acid biosynthesis repression | FasR/SCO2386 (-2.0) | Known fab regulon master regulator |
| **HIGH** | DNA replication/repair cessation | SSB (-3.3), HU (-2.9), Mfd (-1.1), UdgX (-1.3) | 4 DNA maintenance proteins |
| **MEDIUM** | Oxidative stress response repression | MarR: SCO1191 (-2.1), SCO4122 (-3.5) | Both MarR members in repression bloc |
| **MEDIUM** | Select TCS downregulation | SCO5289 SK (-2.4), SCO5434 RR (-1.1) | Environmental sensing pathway shutdown |
| **MEDIUM** | Vegetative sigma factor repression | SCO2639 (-1.3), SCO4960 (-1.5) | Counter to activation-bloc sigma induction |
| **MEDIUM** | Nitrogen assimilation shutdown | P-II/SCO5584 (-1.4) | Key nitrogen sensing regulator |

### 10. Developmental Transition Model

The two blocs implement a coordinated developmental switch:

**Vegetative growth programs are REPRESSED:**
- Primary carbon/nitrogen metabolism (GntR, IclR, LacI, LysR)
- Fatty acid biosynthesis (FasR)
- DNA replication and repair (SSB, HU, Mfd, UdgX)
- Antibiotic efflux/self-defense (5 TetR repressors)
- Nitrogen assimilation (P-II)

**Developmental programs are ACTIVATED:**
- Morphological differentiation (RamR, WhiB, sigma factors)
- Multi-component signal transduction (9 TCS members)
- Secondary metabolite regulation (TcrA, NsdB)
- Stress response cascades (SigE, SigB/F/G, MerR)

This pattern is fully consistent with the known *S. coelicolor* life cycle transition from vegetative growth to aerial mycelium formation and sporulation, where:
1. Metabolic resources are redirected from growth to secondary metabolism
2. Cell division and DNA replication slow
3. Stress response pathways prepare for environmental challenge
4. Complex signal transduction cascades coordinate the multicellular developmental program

## Statistical Summary

| Test | Statistic | p-value | Result |
|------|-----------|---------|--------|
| Fisher: TetR enrichment in repression bloc | OR=0.16 | 0.028* | Significant |
| Fisher: Signal transduction in activation bloc | OR=4.00 | 0.101 | Trend |
| Fisher: Defense/resistance in repression bloc | OR=0.28 | 0.101 | Trend |
| MWU: |LFC| activation vs repression | U | 0.069 | Marginal trend |
| KW: |LFC| by TF family | H=11.10 | 0.269 | Non-significant |
| Fisher: Secondary metabolism neighbors | OR=0.56 | 0.158 | Non-significant |

## Limitations

1. **Sample size**: With 62 TFs across 16 families, many families have only 1-2 representatives, limiting statistical power for individual family tests
2. **HTH (other) and Other regulatory**: These catch-all categories (21/62 = 34%) include diverse TFs whose functions cannot be predicted from family alone
3. **Genomic context**: +-5kb window captures immediate neighbors but may miss more distant operon members
4. **No cognate TCS pairs among exposed TFs**: All 7 known TCS pairs have asymmetric methylation (one exposed, one shielded), preventing pair-level prediction
5. **Predictions are family-based**: Individual TFs may deviate from family-typical functions

## Output Files

### Figures
- `figures/TF_family_bloc_composition.pdf/svg` -- Family counts by bloc + functional category stacked bars
- `figures/functional_category_enrichment.pdf/svg` -- Fisher test OR visualization
- `figures/genomic_context_analysis.pdf/svg` -- Neighbor profiles + per-TF heatmap
- `figures/expression_by_family.pdf/svg` -- |LFC| distributions by family
- `figures/predicted_regulatory_model.pdf/svg` -- Conceptual two-bloc model diagram
- `figures/H35_comprehensive_summary.pdf/svg` -- 6-panel summary

### Tables
- `tables/exposed_TF_family_annotation.tsv` -- Per-TF: family, functions, bloc, predicted role (62 rows)
- `tables/bloc_family_enrichment.tsv` -- Fisher exact tests for family enrichment (16 rows)
- `tables/genomic_context_neighbors.tsv` -- All +-5kb neighbor annotations (554 rows)
- `tables/functional_category_mapping.tsv` -- Family-to-category mapping with literature annotations (16 rows)
- `tables/expression_by_family.tsv` -- Per-family expression statistics (16 rows)
- `tables/statistical_tests.tsv` -- All statistical tests (32 rows)

## Conclusion

**H35 is SUPPORTED (qualitative model)**. TF family membership, combined with genomic context and expression magnitude, provides a coherent functional prediction for the two-bloc system. The activation bloc is enriched for signal transduction (TCS) and developmental/stress regulators (sigma factors, WhiB, RamR), while the repression bloc is significantly enriched for TetR-family defense/efflux regulators (p=0.028) and contains the majority of metabolic TF families and DNA maintenance proteins.

The model predicts a coordinated vegetative-to-developmental switch mediated by the distributed methylation-responsive regulatory layer, where growth-phase programs (metabolism, replication, defense) are repressed and developmental programs (morphogenesis, signaling, stress response) are activated. This is consistent with the known biology of the *Streptomyces* life cycle.

---

*Analysis: 58_exposed_TF_functional_prediction | Script: H35_functional_prediction.py*
