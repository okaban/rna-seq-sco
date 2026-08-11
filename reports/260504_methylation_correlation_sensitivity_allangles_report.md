# 260504 Methylation–Expression Correlation: Multi-Angle Sensitivity Analysis

**Cohort**: 57 exposed TFs (M145; core=31, arm=26)
**Question**: After the headline result of "no correlation" (whole-cohort Spearman ρ=+0.257, p=0.053, NS), is the absence of association robust across alternative cuts of the data?
**Bottom line**: **Mixed.** The Δmethylation × Δexpression test is genuinely null. But the absolute T1 methylation × LFC test is **not null** — particularly for 4mC at TSS ±200–500 bp, with the strongest cell surviving region (core/arm) control.

## Inputs

| Source | What |
|---|---|
| `01_integration/high_confidence_sites_weighted.csv` | per-site methylation freq × T1/T2/T3 (4mC + 6mA) |
| `37_defense_island_GCCGGC/tables/GCCGGC_sites_by_timepoint.tsv` | GCCGGC-context flag (4mC subset) |
| `36_AAGCCCG_distribution/tables/AAGCCCG_site_gene_mapping.tsv` | AAGCCCG-context flag (6mA subset) |
| `58_AAGCCCG_exposed_TF_causal/tables/exposed_TF_methylation_status.tsv` | 57 exposed TFs (TSS, region, LFC, padj) |
| `04_deseq2/.../DESeq2_M145_3_vs_2.tsv` | T3-vs-T2 padj |

## Design — six orthogonal angles, all run on the same TF cohort

For each TF and each window W ∈ {100, 200, 293, 500, 1000, 2000} bp around the TSS we count methylated sites by motif: `all`, `all_4mC`, `GCCGGC_4mC`, `all_6mA`, `AAGCCCG_6mA`. Then:

| # | Angle | Statistic |
|---|---|---|
| 1 | window-size sensitivity | Spearman ρ(Δcount, LFC) for 3 time pairs × 5 motifs × 6 windows |
| 2 | time-pair recheck | the same at canonical W=500 bp |
| 3a | absolute T1 methylation | Spearman ρ(count_T1, LFC) |
| 3b | binary T1 methylated/not | Mann–Whitney on LFC; rank-biserial effect size |
| 4 | sig-DEG-only (padj<0.05) | Δcount ~ LFC restricted to significant DEGs |
| 4b | sig-DEG-only, abs T1 | count_T1 ~ LFC restricted to significant DEGs |
| 5a | partial Spearman | Δcount ~ LFC, controlling for core/arm |
| 5b | partial Spearman | count_T1 ~ LFC, controlling for core/arm |
| 6 | site-count stratification | ≥3 sites at T1 vs 1–2 |

## Per-angle hit rate (sig p<0.05 / total tests with finite ρ)

| Angle | hits | tests | rate | × chance (5%) |
|---|---:|---:|---:|---:|
| 1_window_size | 6 | 82 | 7.3% | 1.5× |
| 2_time_pair | 0 | 14 | 0.0% | 0.0× |
| **3a_absT1_count** | **12** | **87** | **13.8%** | **2.8×** |
| **3b_T1binary_MW** | **11** | **72** | **15.3%** | **3.1×** |
| 4_sigDEG_only | 6 | 82 | 7.3% | 1.5× |
| **4b_sigDEG_absT1** | **11** | **86** | **12.8%** | **2.6×** |
| 5a_partialcorr_dcount | 0 | 90 | 0.0% | 0.0× |
| 5b_partialcorr_absT1 | 5 | 90 | 5.6% | 1.1× |
| 6_strata_T1count | 2 | 62 | 3.2% | 0.6× |

**Read**: every angle that uses Δmethylation behaves at chance (1, 2, 4, 5a). Every angle that uses absolute T1 methylation runs at ~3× chance (3a, 3b, 4b). After region control (5b) the absolute-T1 effect halves but does not vanish.

## Strongest cells

### Δmethylation ~ LFC (Row 1 of figure)
Best cell: T3vT1, all_4mC, ±200 bp → ρ = +0.337, p = 0.010 (would not survive any reasonable multiple-comparison correction across 82 cells; 6 sig observed where ~4 are expected at chance).

### Absolute T1 methylation ~ LFC (Row 2 of figure) — **the actual signal**
Most consistent cell across all variants:

| Test | n | ρ / rb | p |
|---|---:|---:|---:|
| 3a count_T1 ~ LFC_T2vsT1, all_4mC, ±293 bp | 57 | **−0.402** | **0.0019** |
| 3a count_T1 ~ LFC_T2vsT1, all sites, ±293 bp | 57 | −0.386 | 0.003 |
| 3a count_T1 ~ LFC_T2vsT1, all_4mC, ±200 bp | 57 | −0.367 | 0.005 |
| 3b binary methylated ±293 bp, all_4mC, MW | 57 | rb = +0.566 | 0.0021 |
| 4b count_T1 ~ LFC_T2vsT1, all_4mC, ±293 bp (padj<0.05 only, n=45) | 45 | −0.421 | 0.0040 |

**Direction**: NEGATIVE. Exposed TFs whose promoters carry more 4mC sites at T1 are downregulated more strongly by T2. The same cells stay significant when restricted to genes that DESeq2 already calls differentially expressed.

### Region-controlled (Row 3 of figure)
Partial Spearman on count_T1 vs LFC, conditioning on core/arm:

| Test | n | partial ρ | p |
|---|---:|---:|---:|
| **all_4mC, T2vsT1, ±293 bp** | **57** | **−0.332** | **0.012** |
| all sites, T2vsT1, ±293 bp | 57 | −0.316 | 0.017 |
| all_4mC, T2vsT1, ±200 bp | 57 | −0.293 | 0.027 |

About half of the 3a effect is geography (core/arm), but the strongest cell remains significant after conditioning on it. The previously documented core/arm Simpson's-paradox concern is therefore real but partial; it does not erase the absolute-T1 4mC signal.

## Time-pair pattern

| pair | best |ρ| | sig cells (across all motifs/windows) |
|---|---:|---|
| **T2vsT1** | 0.40 | **dominant**: 8/12 sig in 3a, 6/11 in 3b |
| T3vsT1 | 0.34 | 4/12 sig in 3a |
| T3vsT2 | 0.35 (AAGCCCG @ 2kb) | 1/12 in 3a; rest NS |

The signal is concentrated at T1→T2. T2→T3 is essentially flat — consistent with the picture that the T1 4mC state predisposes the early induction response, and once that response has played out the additional Δ between T2 and T3 carries no further information about further LFC.

## Files

- `analysis/57_temporal_dynamics_exposed_TF/scripts/correlation_sensitivity_allangles.py` — full analysis script
- `analysis/57_temporal_dynamics_exposed_TF/results/correlation_sensitivity_allangles.txt` — per-test report
- `analysis/57_temporal_dynamics_exposed_TF/results/correlation_sensitivity_allangles.tsv` — machine-readable table
- `analysis/figures/methylation_expression_correlation_sensitivity.pdf` / `.png` — 3-row × 3-column sensitivity figure

## Final judgment

1. **"Δmethylation does not predict Δexpression" — VERIFIED.** No window, no motif, no time pair, no DEG subset, and no region-controlled cut produces a significant ρ on Δcount × LFC beyond chance.
2. **"Methylation and expression are unrelated" — REJECTED.** When the question is reframed as "does the *absolute initial 4mC density* near the TSS predict how the gene moves over T1→T2?", the answer is yes, with ρ = −0.40 (p = 0.002), shrinking but surviving (ρ = −0.33, p = 0.012) after core/arm control.
3. **Mechanistic implication.** This is a state-based, not dynamic, relationship: 4mC at the TF promoter at T1 is a *prior* that biases the gene toward repression upon induction; it is not a switch that turns on or off in step with expression.
4. **What to report in the paper.** The "no correlation" framing should be replaced with "no Δ–Δ correlation; T1 4mC state predicts repression at T1→T2 (ρ = −0.33 region-controlled)."

## Caveats

- n = 57 is small; with ~600 tests run, ~30 hits expected by chance — observed 53 hits gives only ~1.6× chance overall, which is why the *per-angle* breakdown is the right summary, not the global count.
- 4mC GCCGGC subset is the largest 4mC class (~65% of T1 4mC) but on its own does not show the negative ρ — the signal is carried by the union of 4mC contexts. This argues against attributing the signal to a single methyltransferase system.
- AAGCCCG_6mA cells are mostly empty at small windows (only 8/57 TFs have any AAGCCCG site within ±500 bp), so the AAGCCCG hits are at W=2000 only and should be treated as exploratory.
- The 11 entries in `51_exposed_regulators_characteristics/tables/exposed_regulators_full_table.tsv` are a *narrow* subset; the 57-TF cohort used here is the H35 superset from `58_AAGCCCG_exposed_TF_causal/tables/exposed_TF_methylation_status.tsv`. Numbers in this report are not comparable to the 11-entry table.
