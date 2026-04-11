# H3: Quantitative Resolution of the redZ Paradox

**Date:** 2026-02-24  
**Analysis directory:** `11_epigenome_integration/analysis/26_redZ_paradox/`  
**Hypothesis:** How does the Red BGC activate despite redZ downregulation?

## Background

The Red (undecylprodiginine/streptorubin) biosynthetic gene cluster (BGC) in *S. coelicolor* A3(2) M145 is regulated by a hierarchical cascade. The cluster-situated regulator (CSR) **redZ** (SCO5881) is a response regulator that activates **redD** (SCO5877), a SARP-family pathway-specific activator that directly induces structural gene transcription. Additionally, the **absA1/absA2** two-component system acts as a negative regulator, while global regulators **afsR**, **afsS**, and **afsQ1** provide additional regulatory inputs.

The "redZ paradox" refers to the expectation that if redZ is downregulated, the Red BGC should not activate -- yet the data shows robust Red BGC activation across timepoints.

## Key Finding: There Is No Paradox

**Critical observation: redZ is NOT downregulated.** The data shows redZ is actually *upregulated* across all timepoints:

| Gene | T1 (exp) | T2 (trans) | T3 (stat) | log2FC T2vsT1 | log2FC T3vsT1 |
|------|----------|------------|-----------|---------------|---------------|
| redZ | 690.1 | 1,319.5 | 1,377.3 | +0.91* | +0.97* |

redZ expression nearly doubles from T1 to T2/T3 (padj < 0.05 in both comparisons). Therefore, the original premise of "redZ downregulation" is not supported by this dataset.

## Regulatory Landscape

### Regulator Expression Timeline

| Gene | Role | T1 mean | T2 mean | T3 mean | L2FC T2v1 | sig | L2FC T3v1 | sig |
|------|------|---------|---------|---------|-----------|-----|-----------|-----|
| redD | activator | 30.6 | 1,737.5 | 1,008.5 | +5.79 | * | +5.00 | * |
| redZ | CSR | 690.1 | 1,319.5 | 1,377.3 | +0.91 | * | +0.97 | * |
| absA1 | repressor | 265.5 | 1,255.1 | 1,802.2 | +2.20 | * | +2.73 | * |
| absA2 | repressor | 481.1 | 950.2 | 712.4 | +0.95 | * | +0.54 | * |
| afsR | activator | 985.2 | 818.4 | 517.1 | -0.26 | ns | -0.92 | * |
| afsS | activator | 253.5 | 2,155.4 | 2,194.6 | +3.04 | * | +3.07 | * |
| afsQ1 | activator | 1,935.3 | 696.9 | 379.8 | -1.44 | * | -2.33 | * |
| afsQ2 | sensor | 2,120.5 | 871.0 | 422.3 | -1.26 | * | -2.31 | * |
| MTase (SCO3104) | epigenetic | 245.1 | 51.7 | 145.1 | -2.19 | * | -0.71 | * |

### Red BGC Structural Gene Activation

The Red BGC structural genes (20 genes, excluding regulators redD and redZ) show massive upregulation:

- **Mean normalized counts:** T1=35.7, T2=596.9, T3=942.8
- **Mean log2FC:** T2vsT1=+3.85 (14.5-fold), T3vsT1=+5.03 (32.7-fold)
- All structural genes show strong, consistent upregulation

## Activator-Repressor Balance Analysis

### Balance Scores

| Component | T2vsT1 | T3vsT1 |
|-----------|--------|--------|
| Activator sum (redD + afsR + afsS + afsQ1) | +7.13 | +4.82 |
| Repressor sum (absA1 + absA2) | +3.15 | +3.27 |
| Net balance (Act - Rep) | +3.98 | +1.56 |
| redZ log2FC | +0.91 | +0.97 |
| Red BGC structural mean log2FC | +3.85 | +5.03 |

### Detailed Regulatory Accounting (T2vsT1)

| Input | log2FC | Direction |
|-------|--------|-----------|
| redD activation | +5.79 | Strong positive |
| afsS activation | +3.04 | Strong positive |
| redZ increase | +0.91 | Moderate positive |
| afsR change | -0.26 | Negligible |
| afsQ1 decrease | -1.44 | Moderate negative |
| absA1 increase | +2.20 | Counter-productive (more repression) |
| absA2 increase | +0.95 | Counter-productive (more repression) |

### Key Interpretation

1. **redD is the dominant driver** -- with a log2FC of +5.79 (55-fold increase) at T2, redD overwhelms all other regulatory signals. redD's massive upregulation is the primary mechanism activating the Red BGC structural genes.

2. **afsS provides strong support** -- the afsS anti-anti-sigma factor (+3.04 log2FC) is strongly co-upregulated, likely contributing to the activation cascade via afsR-mediated signaling.

3. **absA1 upregulation is paradoxical** -- the negative regulator absA1 is also upregulated (+2.20), which would be expected to *repress* Red production. This may represent a homeostatic feedback mechanism or may indicate that absA1's repressive activity requires phospho-transfer from absA2, which shows only moderate change.

4. **afsQ1/afsQ2 decline** -- the afsQ TCS system is strongly downregulated, removing one positive input. However, this loss is more than compensated by the redD and afsS increases.

5. **Full net balance is positive** at both timepoints (+4.89 at T2, +2.52 at T3), consistent with the observed Red BGC activation.

## Revised Interpretation

The original hypothesis framed this as a "paradox" suggesting redZ downregulation should prevent Red BGC activation. However, **the data shows redZ is upregulated, not downregulated**. The regulatory picture is:

- redZ is modestly upregulated (+0.91 log2FC), maintaining its activating role
- redD responds disproportionately (+5.79 log2FC), suggesting signal amplification in the redZ -> redD cascade
- The massive redD response could also reflect additional inputs from global regulators (afsS, afsR) converging on the redD promoter
- The true regulatory puzzle is why absA1 is simultaneously upregulated -- this may represent negative feedback to prevent runaway Red production

## Output Files

### Tables
- `tables/regulator_timeline.tsv` -- Expression and DESeq2 stats for all regulators
- `tables/red_bgc_expression.tsv` -- Expression data for all 22 Red BGC genes
- `tables/activator_repressor_balance.tsv` -- Quantitative balance scores

### Figures
- `figures/redZ_paradox_multipanel.pdf/svg` -- Multi-panel overview (A: regulator bars, B: BGC vs regulators, C: heatmap)
- `figures/balance_diagram.pdf/svg` -- Horizontal bar chart of regulatory log2FC balance
- `figures/redZ_vs_redD_timecourse.pdf/svg` -- Time course with individual replicates
- `figures/balance_vs_output.pdf/svg` -- Regulatory balance vs BGC output scatter

## Conclusions

1. **No redZ paradox exists in this dataset** -- redZ is upregulated, not downregulated
2. **redD is the dominant activator** of Red BGC structural genes, with 55-fold upregulation at T2
3. **afsS provides a strong co-regulatory signal** (+3.04 log2FC)
4. **The absA system is simultaneously upregulated**, suggesting active negative feedback during Red production
5. **Net regulatory balance is positive**, consistent with the observed 14.5-fold (T2) to 32.7-fold (T3) increase in Red BGC structural gene expression
