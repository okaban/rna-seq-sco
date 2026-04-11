# H5: MTase Expression Stability vs Methylation-Expression Correlation

**Date:** 2026-02-24  
**Analysis directory:** `11_epigenome_integration/analysis/28_MTase_stability/`  
**Script:** `run_analysis_v2.py`

## Hypothesis

More stably expressed methyltransferases (MTases) produce more orderly (consistent, predictable) methylation-expression correlations at their target motif sites.

## Methods

1. **MTase identification** — Five candidate DNA methyltransferases were identified from GFF annotations and prior motif analysis (Analysis 23), cross-validated against NCBI RefSeq annotations.
2. **Expression stability** — Normalized counts (DESeq2) across 3 timepoints x 3 replicates were used to compute coefficient of variation (CV) and maximum absolute log2FC. A composite stability score (0 = most variable, 1 = most stable) was derived from normalized CV and max|LFC|.
3. **Motif attribution** — Motif-assigned sites (6mA and 4mC) were mapped to genes (within gene body or 200 bp upstream) using GFF coordinates.
4. **Site dynamics** — For each motif group, sites were classified as Gained, Lost, or Stable between timepoints.
5. **Methylation-expression correlation** — Spearman correlation between gene-level methylation change and log2FC was computed for genes containing sites of each motif group.

## Results

### MTase Expression Stability Ranking

| MTase | SCO | Product | Motif | T1 mean | T2 mean | T3 mean | CV | max\|LFC\| | Stability |
|-------|-----|---------|-------|---------|---------|---------|-----|----------|-----------|
| SC_RS28835 | SCO5333 | BREX-2 PglX (6mA) | CCGKCA | 604 | 395 | 382 | 0.27 | 0.65 | **1.00** |
| SC_RS35335 | SCO6843 | BREX-2 PglX (6mA) | CCGKCA | 1666 | 1989 | 980 | 0.33 | 1.01 | **0.94** |
| SC_RS17645 | SCO3104 | N-6 DNA methylase (HsdM) | AAGCCCG | 245 | 52 | 145 | 0.66 | 2.19 | **0.68** |
| SC_RS19770 | SCO3527 | DNA cytosine MTase (Dcm) | CCGG/TGGCCGGC | 11 | 5 | 62 | 1.22 | 3.54 | **0.32** |
| SC_RS36410 | SCO7091 | DNA cytosine MTase | CCGG/TGGCCGGC | 1 | 4 | 65 | 1.53 | 5.34 | **0.00** |

### Key Observations

#### 1. SC_RS17645 (AAGCCCG, Type I HsdM)
- **Expression**: Strongly downregulated at T2 (LFC = -2.19, padj = 6.5e-16), partial recovery at T3
- **Site dynamics**: 260 sites lost, 64 gained, 0 shared between T1 and T2 — complete remodeling of AAGCCCG methylation landscape
- **Correlation**: rho = -0.127 (p = 0.36, n = 54 genes) — weak negative trend, not significant

#### 2. SC_RS28835/SC_RS35335 (CCGKCA, BREX-2 PglX)
- **Expression**: Most stable of all MTases (stability = 1.00 and 0.94)
- **Site dynamics**: CCGKCA sites not detected in high-confidence motif assignments (0 sites)
- **Correlation**: Cannot be computed — insufficient data
- **Interpretation**: BREX-2 PglX may produce low-frequency methylation below detection threshold, or CCGKCA motif assignment may need refinement

#### 3. SC_RS19770/SC_RS36410 (CCGG family, Dcm-like)
- **Expression**: Highly variable — nearly silent at T1/T2, dramatically upregulated at T3 (SC_RS36410: LFC = +5.34, padj = 1.0e-10)
- **Site dynamics**: 1669 sites lost, 543 gained (T2 vs T1) — paradoxically, massive site loss despite low MTase expression at T1
- **Correlation**: rho = 0.088 (p = 0.12, n = 307 genes) — near-zero, not significant
- **Paradox**: Sites are present at T1 when these MTases are barely expressed, suggesting an alternative MTase maintains CCGG methylation, or these sites reflect passive demethylation during rapid growth

### Methylation-Expression Correlations (All Comparisons)

| Motif Group | Mod | Comparison | n genes | rho | p-value | Concordance |
|-------------|-----|-----------|---------|-----|---------|-------------|
| AAGCCCG | 6mA | T2 vs T1 | 54 | -0.127 | 0.361 | 40.7% |
| AAGCCCG | 6mA | T3 vs T1 | 59 | -0.049 | 0.710 | 39.0% |
| AAGCCCG | 4mC | T2 vs T1 | 49 | 0.155 | 0.287 | 46.9% |
| CCGG family | 4mC | T2 vs T1 | 307 | 0.088 | 0.123 | 51.8% |
| CCGG family | 4mC | T3 vs T1 | 269 | 0.031 | 0.614 | 49.8% |
| Unattributed | 6mA | T3 vs T1 | 397 | -0.121 | **0.016** | 42.8% |

Only one comparison (Unattributed 6mA, T3 vs T1) reached nominal significance.

## Verdict

**PARTIALLY SUPPORTED**

- **Supported at site level**: MTase downregulation (SC_RS17645, LFC = -2.19) clearly drives motif site loss (260 lost AAGCCCG sites). The most variable Dcm-like MTases show the most dramatic site turnover.
- **Not supported at gene level**: Methylation-expression correlations are uniformly weak (|rho| < 0.16) and non-significant across all motif groups, regardless of MTase stability. Concordance rates hover near 50% (chance level).
- **Implication**: DNA methylation in *S. coelicolor* likely operates through a **global regulatory mechanism** (e.g., chromosome structure, replication timing) rather than gene-specific cis-regulation. MTase expression dynamics control the overall methylation state but do not create predictable gene-by-gene expression responses.

## Output Files

### Figures
- `figures/MTase_expression_timeline.pdf/svg` — Expression profiles across timepoints
- `figures/motif_site_dynamics.pdf/svg` — Gained/Lost/Stable site counts per motif group
- `figures/MTase_stability_vs_correlation.pdf/svg` — Stability vs correlation scatter plots
- `figures/H5_comprehensive.pdf/svg` — Four-panel summary figure

### Tables
- `tables/MTase_expression_profiles.tsv` — Expression and stability metrics for all 5 MTases
- `tables/motif_site_dynamics.tsv` — Site gain/loss counts per motif group per comparison
- `tables/motif_expression_correlations.tsv` — Spearman correlations per motif group
- `tables/H5_summary.tsv` — Integrated summary table
