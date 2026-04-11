# H26: TF Binding Site-Level Methylation Protection

**Date:** 2026-02-26
**Analysis directory:** `11_epigenome_integration/analysis/49_TF_BS_methylation_protection/`
**Status:** REJECTED

---

## Background

H25 (SUPPORTED) demonstrated a 2,200 bp "protection zone" around regulatory gene TSS where methylation density is depleted by 17.3% (4.6x more than non-regulatory genes). H22 showed that promoter regions have no DNA motif depletion, supporting a protein occupancy model where TF/RNAP binding physically blocks MTase access. This analysis tests whether the occupancy effect is detectable at individual TF binding sites using FIMO predictions.

## Hypothesis

Predicted TF binding sites (56,338 FIMO hits) show methylation depletion in their central region, and the effect is stronger at TF BS located near regulatory gene promoters.

## Methods

### Data
- **TF binding sites**: 56,338 FIMO motif matches (p < 1e-4) from genome-wide scan
- **Methylation sites**: GCCGGC T1 (1,289), AAGCCCG T1 (260), All 4mC T1 (1,987), All 6mA T1 (1,088)
- **Genomic context**: Promoter-regulatory (3,274), promoter-non-regulatory (22,813), intergenic (1,076), intragenic (29,175)

### Statistical methods
- Binomial test for overlap enrichment/depletion (±10 bp)
- Spatial density profiles (50 bp bins, ±2 kb from TF BS center)
- 1,000-iteration permutation tests for central depletion
- Fisher exact tests for context-specific comparisons

## Results

### 1. TF BS do NOT show methylation depletion

| Methylation Type | n Overlap | Expected | Fold | p-value | Direction |
|-----------------|-----------|----------|------|---------|-----------|
| GCCGGC T1 | 373 | 322 | **1.157** | **0.006** | **Enriched** |
| AAGCCCG T1 | 41 | 65 | **0.629** | **0.002** | **Depleted** |
| All 4mC T1 | 503 | 496 | 1.014 | 0.752 | NS |
| All 6mA T1 | 542 | 483 | 1.122 | 0.008 | Enriched |
| All T1 | 1,007 | 974 | 1.034 | 0.286 | NS |

**GCCGGC is significantly enriched** at TF BS (fold=1.16, p=0.006), opposite to the predicted depletion. **AAGCCCG is depleted** (fold=0.63, p=0.002), but this is the only motif showing the predicted pattern.

### 2. Spatial profiles show no central depletion

| Methylation Type | Central Density | Flank Density | Ratio | Permutation p |
|-----------------|----------------|---------------|-------|---------------|
| GCCGGC T1 | 0.150 | 0.149 | 1.011 | 0.744 (NS) |
| All 4mC T1 | 0.226 | 0.229 | 0.984 | 0.120 (NS) |
| All 6mA T1 | 0.229 | 0.223 | 1.028 | 0.960 (NS) |
| All T1 | 0.455 | 0.452 | 1.006 | 0.700 (NS) |

No methylation type shows significant central depletion at TF BS. The spatial profiles are essentially flat.

### 3. Regulatory vs non-regulatory promoter TF BS

| Context | All 4mC Central Density | All T1 Central Density | Reg/NonReg Ratio |
|---------|------------------------|----------------------|-----------------|
| Regulatory promoter | 0.198 | 0.405 | — |
| Non-regulatory promoter | 0.218 | 0.438 | — |
| Ratio | 0.909 | 0.926 | reg < nonreg |

Regulatory promoter TF BS have ~7-9% lower methylation density than non-regulatory, but this is modest and consistent with the aggregate gene-level effect (H25), not a TF BS-specific phenomenon.

### 4. Context-specific patterns

| Context | GCCGGC Fold | All 4mC Fold | All 6mA Fold | Notable |
|---------|------------|-------------|-------------|---------|
| Promoter (reg) | 0.956 (NS) | 0.966 (NS) | 1.028 (NS) | No enrichment or depletion |
| Promoter (nonreg) | 1.064 (NS) | 0.940 (NS) | 1.058 (NS) | No enrichment or depletion |
| Intergenic | 0.323 (NS) | 0.315 (p=0.03) | **3.021 (p=4.6e-07)** | 6mA massively enriched |
| Intragenic | **1.284 (p=4e-04)** | 1.103 (NS) | 1.113 (NS) | GCCGGC enriched |

- **Intergenic TF BS**: 6mA is 3x enriched (p=4.6e-07), while 4mC is 3x depleted. Intergenic regulatory elements show strong motif-specific bias.
- **Intragenic TF BS**: GCCGGC is enriched (fold=1.28, p=4e-04), consistent with the gene body methylation pattern.

## Interpretation

### Why the hypothesis is rejected

The prediction was that individual TF binding sites would show methylation depletion due to protein occupancy blocking MTase access. Instead:

1. **No central depletion at TF BS**: Spatial profiles are flat (permutation p=0.12-0.96). If TF binding physically blocks MTase, we should see a "valley" at the BS center — we do not.

2. **GCCGGC is enriched at TF BS**: The opposite of depletion. This likely reflects the GC-rich nature of both GCCGGC recognition sequences and many TF binding motifs (both are palindromic GC-rich sequences in GC-rich Streptomyces genome).

3. **The H25 protection zone is an aggregate phenomenon**: The 2,200 bp protection zone at regulatory gene TSS is NOT explained by individual TF BS occupancy. Rather, it reflects the cumulative effect of multiple proteins (RNAP, sigma factors, multiple TFs, nucleoid-associated proteins) collectively excluding MTase from regulatory gene promoter regions.

### Reconciliation with H25

- **H25**: Regulatory gene TSS shows 17.3% methylation depletion in a 2,200 bp zone
- **H26**: Individual TF BS show no depletion
- **Resolution**: The protection zone is a **collective, gene-level** phenomenon. Regulatory gene promoters have higher overall protein occupancy (more TFs binding, constitutive RNAP engagement) that creates an aggregate exclusion zone. No single TF binding event is sufficient to exclude MTase, but the combined occupancy at regulatory gene promoters is.

### Unexpected finding: Intergenic 6mA enrichment at TF BS

6mA shows 3x enrichment at intergenic TF BS (p=4.6e-07). This could reflect:
- 6mA preferentially marking regulatory intergenic regions
- A functional connection between 6mA and transcription factor binding sites in regulatory regions
- This warrants further investigation

## Verdict

**REJECTED** — Individual TF binding sites do not show methylation depletion. The H25 protection zone at regulatory gene TSS is a collective, aggregate phenomenon of promoter-wide protein occupancy, not attributable to individual TF binding events. GCCGGC is paradoxically enriched at TF BS. AAGCCCG shows weak depletion but this may reflect different sequence composition preferences rather than occupancy effects.

## Output Files

### Figures
| File | Description |
|------|-------------|
| `figures/TFBS_methylation_profile.pdf/svg` | TF BS-centered methylation density profile |
| `figures/regulatory_vs_nonreg_TFBS.pdf/svg` | Regulatory vs non-regulatory promoter comparison |
| `figures/overlap_by_context.pdf/svg` | Genomic context comparison |
| `figures/H26_comprehensive_summary.pdf/svg` | Multi-panel summary |

### Tables
| File | Description |
|------|-------------|
| `tables/TFBS_methylation_overlap.tsv` | Per-TF-BS overlap data (56,338 entries) |
| `tables/spatial_profile_data.tsv` | Binned density values |
| `tables/context_analysis.tsv` | Genomic context statistics |
| `tables/statistical_tests.tsv` | All test results (17 tests) |

### Script
- `scripts/H26_TFBS_methylation_protection.py`

---

*Analysis performed: 2026-02-26*
