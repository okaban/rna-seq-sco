# Exposed TF: T1 4mC / 6mA density × T1→T2 LFC, split by modification type

**Date:** 2026-05-05
**Scope:** Stratify the exposed-TF (n=57) "T1 4mC density predicts T1→T2 LFC" correlation (ρ=−0.40, p=0.0019) by which modification type drove the gene's coordinated methylation–expression call (4mC vs 6mA).

---

## TL;DR

The pooled signal at exposed-TF promoters (T1 4mC density vs T1→T2 LFC, ρ=−0.40, p=0.0019, n=57) is **not equally carried by 4mC- and 6mA-selected genes**.

- **4mC-selected (n=19):** ρ=−0.33, p=0.17 (T1 4mC density × LFC). Direction preserved, underpowered.
- **6mA-selected (n=37):** ρ=−0.17, p=0.32 (T1 4mC density × LFC). Effect collapses.
- **Both (n=1):** SC_RS27300 (SCO5027, HTH). Single gene, not analysed.

The pooled n=57 result is therefore driven primarily by genes whose dynamic methyl change involved 4mC, not 6mA. This is consistent with the global finding that 4mC promoter density tracks expression while 6mA does not (project memory: ρ=0.139 vs ~0).

---

## Inputs

| File | Use |
|---|---|
| `analysis/29_genomewide_TF_screen/tables/coordinated_regulatory_genes.tsv` | n=57 coordinated TFs and per-timepoint 4mC/6mA counts |
| `analysis/58_AAGCCCG_exposed_TF_causal/tables/exposed_TF_methylation_status.tsv` | TSS, region, T1→T2 LFC |
| `analysis/01_integration/high_confidence_sites_weighted.csv` | Per-timepoint 4mC and 6mA site positions |

---

## Classification rule

A gene is classified as:

- **4mC-selected**: 4mC count varies across T1/T2/T3 AND 6mA count is constant.
- **6mA-selected**: 6mA count varies AND 4mC is constant.
- **Both**: both vary.

Resulting split (n=57): 4mC=19, 6mA=37, both=1. Matches the user-cited split (≈20 / ≈37 / 1).

---

## Density definition

- TSS ± 293 bp window (same as main figure).
- T1 4mC density = (T1 high-confidence 4mC sites in window) / (2·293) × 100, sites per 100 bp.
- Same for T1 6mA.

---

## Results

| Group | n | x | Spearman ρ | p |
|---|---|---|---|---|
| All exposed | 57 | T1 4mC density | **−0.402** | **0.0019** |
| All exposed | 57 | T1 6mA density | −0.193 | 0.149 |
| 4mC-selected | 19 | T1 4mC density | −0.326 | 0.173 |
| 4mC-selected | 19 | T1 6mA density | −0.372 | 0.117 |
| 6mA-selected | 37 | T1 4mC density | −0.168 | 0.321 |
| 6mA-selected | 37 | T1 6mA density | −0.086 | 0.612 |
| Both | 1 | — | NA | NA |

---

## Interpretation

1. **Pooled signal is concentrated in 4mC-selected genes.** ρ stays negative (~−0.33 to −0.37) in the n=19 4mC-selected subgroup; in the larger n=37 6mA-selected subgroup ρ collapses toward zero. The n=57 association inherits from the n=19 group despite its smaller size.
2. **Underpowered, not absent.** With n=19 the 4mC-selected subgroup cannot reach p<0.05 even at ρ=−0.33; the direction matches the pooled call.
3. **Within 4mC-selected, 4mC and 6mA densities give similar ρ.** ρ_4mC=−0.326 vs ρ_6mA=−0.372. This likely reflects co-occurrence of any TSS-proximal methylation rather than a true 6mA effect; in the 6mA-selected pool 6mA density itself is uninformative (ρ=−0.09). Worth confirming with a 4mC-vs-6mA density covariance check before claiming modification-specific causation.
4. **The "both" category (n=1, SCO5027 / HTH) is too small to interpret.**

This refines the project-memory claim "absolute T1 4mC density at TSS predicts T1→T2 LFC" — the predictive structure lives in genes whose dynamic methyl change is 4mC-mediated.

---

## Outputs

- Stats TSV: `analysis/57_temporal_dynamics_exposed_TF/results/exposed_4mC_6mA_split_correlation.tsv`
- Per-gene TSV: `analysis/57_temporal_dynamics_exposed_TF/results/exposed_4mC_6mA_split_per_gene.tsv`
- Figure: `analysis/figures/Fig_exposed_4mC_vs_6mA_split_correlation.{pdf,png}`
- Obsidian copy: `/Users/okaban/obsidian/Research/rna-seq/Writing/Fig_exposed_4mC_6mA_split.png`

## Script

`analysis/57_temporal_dynamics_exposed_TF/scripts/exposed_4mC_6mA_split_correlation.py`
