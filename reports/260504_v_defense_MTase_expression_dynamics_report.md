# V-Defense MTase expression dynamics vs whole-genome GCCGGC-m4C density (T1→T2→T3)

Date: 2026-05-04
Strain: *S. coelicolor* A3(2) M145
Hypothesis under test: down-regulation of the 12 V-Defense 4mC-hypomethylated genes (10/12 SAM-dependent methyltransferases) drives the genome-wide collapse of GCCGGC-4mC across T1→T2→T3 — a self-regulatory MTase loop.

Script: `11_epigenome_integration/analysis/57_temporal_dynamics_exposed_TF/scripts/v_defense_MTase_expression_dynamics.py`
Figures: `11_epigenome_integration/analysis/figures/v_defense_MTase_expression_dynamics.{pdf,png}`, `…/v_defense_MTase_vs_global_4mC.{pdf,png}`
Result tables: `analysis/57_temporal_dynamics_exposed_TF/results/v_defense_MTase_expression_long.tsv`, `…_per_gene_timepoint.tsv`, `…_mean_expression_summary.tsv`, `v_defense_global_4mC_density.tsv`, `v_defense_MTase_vs_global_4mC_correlation.tsv`

---

## TL;DR

The hypothesis as stated is **not supported by the mean of the 12-gene set**.

- Whole-genome GCCGGC-m4C density collapses monotonically: **T1 4.10% → T2 1.29% → T3 0.067%** (≈60× drop, T1→T3).
- The mean log₂ expression of the 12 V-Defense MTase genes does **not** mirror this collapse. It actually **peaks at T2** (7.17 → 7.89 → 7.21) because 7/12 genes are *induced* T1→T2.
- Naïve cross-timepoint correlation (n=3): Pearson r = −0.270 (p=0.83), Spearman ρ = −0.50 (p=0.67). Descriptive only at n=3.
- The MTase→methylation coupling is real but **gene-specific**, not set-level. Only **4/12 genes (SC_RS12890, SC_RS12500, SC_RS13305, SC_RS38280)** decline monotonically T1>T2>T3, and within that subset the per-timepoint mean tracks density tightly (Pearson r = +0.975, Spearman ρ = +1.0, n=3). Of these four, **SC_RS12890 (SCO2170, class I SAM-MTase)** and **SC_RS38280 (SCO7213, SAM-MTase)** are the only DNA-MTase-plausible candidates; SC_RS12500/rsmH is a 16S rRNA MTase and SC_RS13305/panB is a pantothenate biosynthesis enzyme.

Reframed conclusion: a 60× drop in GCCGGC-m4C is **not** explained by a coordinated MTase-set shutdown. It is consistent with **focal silencing of one or two specific GCCGGC MTases** (top candidate SCO2170 / SC_RS12890; secondary SCO7213 / SC_RS38280), against a backdrop of broad MTase induction at T2.

---

## Inputs and methods

| Resource | Path |
|---|---|
| 12-gene V-Defense × 4mC-hypo set | `analysis/57_temporal_dynamics_exposed_TF/results/v_defense_hypo_genes.tsv` |
| DESeq2 size-factor normalized counts | `04_deseq2/results/results/normalized_counts_M145.tsv` |
| Per-site 4mC SMRT calls (T1/T2/T3, freq in %) | `analysis/23_expanded_motif_search/4mC_final_census.csv` |
| Genome reference | `11_epigenome_integration/data/NC_003888.3.fna` |

Sample mapping: `M145_1_*` → T1, `M145_2_*` → T2, `M145_3_*` → T3, n=3 each.
Expression metric: log₂(normalized count + 1), per replicate then averaged.
GCCGGC-m4C density: count of called 4mC sites at GCCGGC motifs (TGGCCGGC / GGCCGG context, `final_motif` in census) with frequency ≥ 50%, divided by 31,452 strand-resolved palindromic GCCGGC positions in the reference (15,726 motif occurrences × 2 strands; GCCGGC is palindromic).
Correlation: Pearson and Spearman over n=3 timepoint means; reported as descriptive only.

---

## Analysis 1 — 12 V-Defense 4mC-hypo genes: T1/T2/T3 expression trajectory

Figure: `analysis/figures/v_defense_MTase_expression_dynamics.{pdf,png}`

Per-gene mean log₂(norm count + 1):

| gene_id | SCO | product | T1 | T2 | T3 | T1→T2 |
|---|---|---|---|---|---|---|
| SC_RS05450 | SCO0705 | SAM-dependent methyltransferase | 5.84 | 9.70 | 7.81 | ↑ |
| SC_RS06105 | SCO0835 | class I SAM-dependent methyltransferase | 9.32 | 9.97 | 11.07 | ↑ |
| SC_RS06595 | SCO0929 | SAM-dependent methyltransferase | 4.42 | 7.27 | 5.74 | ↑ |
| SC_RS11875 | SCO1969 | methylated-DNA–[protein]-cysteine S-MTase | 5.53 | 11.34 | 8.06 | ↑ |
| SC_RS12500 | SCO2092 | rsmH (16S rRNA m4C MTase) | 9.45 | 6.93 | 6.58 | ↓ |
| SC_RS12890 | SCO2170 | class I SAM-dependent methyltransferase | 10.27 | 8.36 | 7.19 | ↓ |
| SC_RS13205 | SCO2235 | type II TA antitoxin (Phd/YefM) | 5.78 | 7.93 | 6.69 | ↑ |
| SC_RS13305 | SCO2256 | panB (3-methyl-2-oxobutanoate hydroxymethyltransferase) | 8.59 | 6.81 | 8.32 | ↓ |
| SC_RS26345 | SCO4837 | glycine hydroxymethyltransferase | 8.38 | 9.31 | 8.45 | ↑ |
| SC_RS38275 | – | methyltransferase type 11 | 7.97 | 7.05 | 6.81 | ↓ |
| SC_RS38280 | SCO7213 | SAM-dependent methyltransferase | 6.49 | 6.27 | 5.20 | ↓ |
| SC_RS39465 | SCO7452 | methyltransferase | 4.04 | 6.34 | 6.20 | ↑ |

T1→T2 direction: 5/12 down, 7/12 up. T1→T3 direction: 5/12 down. The set is **bimodal**, not coordinately repressed.

Monotonic decliners (T1 > T2 > T3), n=4:
- **SC_RS12890 / SCO2170** — class I SAM-MTase (T1 LFC −1.93, T3 LFC −3.22). DNA-MTase-plausible.
- **SC_RS12500 / SCO2092 / rsmH** — 16S rRNA (cytosine(1402)-N4)-methyltransferase. Acts on rRNA, **not DNA** — exclude as causal driver of GCCGGC.
- **SC_RS13305 / SCO2256 / panB** — pantothenate biosynthesis, **not a DNA MTase** — exclude.
- **SC_RS38280 / SCO7213** — SAM-MTase. DNA-MTase-plausible (paralog/ adjacent to SC_RS38275, also down).

So among the 4 monotonic decliners, **only SC_RS12890 and SC_RS38280** are credible candidates for being a GCCGGC DNA MTase whose decline could shrink the 4mC pool. Combined with neighbouring SC_RS38275 (also declining T2→T3, COG V), the SC_RS38275/SC_RS38280 locus is the second strongest candidate region.

---

## Analysis 2 — Mean MTase expression vs whole-genome GCCGGC-m4C density

Figure: `analysis/figures/v_defense_MTase_vs_global_4mC.{pdf,png}`

Whole-genome GCCGGC-m4C density (sites with freq ≥ 50% / 31,452 strand-resolved GCCGGC palindromic positions):

| timepoint | sites called | mean freq (%) | density ≥50% |
|---|---|---|---|
| T1 | 1,289 | 83.7 | **4.10 %** |
| T2 | 407   | 78.4 | **1.29 %** |
| T3 | 21    | 68.3 | **0.07 %** |

Mean log₂(norm count + 1) of the 12 genes (per-replicate mean across 3 reps, then mean):

| timepoint | mean log₂ MTase expression | SD |
|---|---|---|
| T1 | 7.169 | 0.086 |
| T2 | **7.892** | 0.066 |
| T3 | 7.211 | 0.049 |

Cross-timepoint correlation (n=3, descriptive):

| metric | n | Pearson r | p | Spearman ρ | p |
|---|---|---|---|---|---|
| mean log₂ (12 genes) vs density ≥50% | 3 | **−0.270** | 0.83 | **−0.50** | 0.67 |
| mean log₂ (4 monotonic decliners) vs density ≥50% | 3 | **+0.975** | 0.14 | **+1.000** | – |

The "all-12 mean" view does not produce a coherent positive trend because the set is dominated by induced genes at T2 (e.g., SC_RS11875 LFC +5.5, SC_RS06595 LFC +2.7, SC_RS39465 LFC +2.1, SC_RS06105 LFC +0.6). Only when the analysis is restricted to the 4 monotonic decliners does the per-timepoint mean track GCCGGC-m4C density tightly (and even then n=3 means the correlation must be treated as descriptive).

---

## Interpretation

1. **The "self-regulatory MTase set" framing is a Simpson-style aggregation artefact.** The 12-gene mean hides two opposing trajectories. Treating "COG V Defense + 4mC-hypo at promoter" as a single coordinated module is misleading at the expression level.

2. **The MTase that plausibly *catalyses* GCCGGC m4C is a single gene (or a small focal cluster), not a 12-gene module.** The 60× collapse of GCCGGC-m4C is consistent with focal silencing of one or two specific DNA MTases. Within this 12-gene set the credible candidates are:
   - **SC_RS12890 / SCO2170** — class I SAM-MTase, monotonic decline (LFC −1.93 / −3.22), 4mC promoter loss −95%.
   - **SC_RS38275 / SC_RS38280 (SCO7213)** locus — both COG V SAM-MTases, both declining at T2→T3 (T3vsT2 LFC −0.20 and −0.63), 4mC promoter loss ≈ −87%.

   These are the genes to cross-check against the **REBASE GCCGGC / GGCCGG MTase candidate list** (`analysis/33_GCCGGC_RM_identification/tables/REBASE_GCCGGC_GGCCGG_entries.tsv`) and the **MTase reverse-ID hits** (`analysis/39_GCCGGC_MTase_reverse_ID/`).

3. **The 7/12 induced-at-T2 genes are not GCCGGC catalysts.** Three of them have non-DNA-MTase products (rsmH = rRNA, panB = pantothenate, glycine hydroxymethyltransferase), and the others (SCO0705, SCO0835, SCO0929, SCO1969, SCO7452) lack any current evidence linking them to GCCGGC. Their COG-V "Defense" annotation comes from generic SAM-MTase membership, not from a confirmed DNA-MTase function. Their induction at T2 is more parsimoniously explained by stress/secondary-metabolism activation than by an anti-phage program.

4. **What this analysis does *not* prove.** Per-site quantitative coupling between an individual MTase's expression trajectory and the 4mC frequency of *its* target sites would require (a) confirmed substrate assignment per MTase and (b) per-site frequency time-series, not just a binary called/not-called per timepoint. The current SMRT census drops below the detection floor for low-modification sites at T2/T3, so part of the 60× drop is also a sensitivity-floor effect, not just biological loss.

---

## Files written

| File | Description |
|---|---|
| `analysis/figures/v_defense_MTase_expression_dynamics.{pdf,png}` | per-gene T1/T2/T3 expression trajectories (12 lines) |
| `analysis/figures/v_defense_MTase_vs_global_4mC.{pdf,png}` | dual-axis: mean MTase expression (left) vs GCCGGC-m4C density (right) |
| `…/57_temporal_dynamics_exposed_TF/results/v_defense_MTase_expression_long.tsv` | per-replicate normalized counts (108 rows) |
| `…/v_defense_MTase_expression_per_gene_timepoint.tsv` | per-(gene, timepoint) mean log₂ expression |
| `…/v_defense_MTase_mean_expression_summary.tsv` | per-timepoint mean log₂ across 12 genes |
| `…/v_defense_global_4mC_density.tsv` | per-timepoint called-sites and density |
| `…/v_defense_MTase_vs_global_4mC_correlation.tsv` | Pearson/Spearman summary |

---

## Recommended follow-up

1. **Cross-reference the 4 monotonic decliners against `33_GCCGGC_RM_identification` and `39_GCCGGC_MTase_reverse_ID`** — does SCO2170 and/or SCO7213 correspond to a homolog of a known GCCGGC C5/N4 MTase (NaeI / NgoMIV / SfoI families)?
2. **Per-site frequency time-series, not just called/not-called.** Re-run the analysis using continuous methylation frequencies (kineticsTools per-position output) at every genomic GCCGGC, joining SCO2170 and SCO7213 trajectories per site rather than collapsing to genome-wide density.
3. **Disambiguate the T2 detection-floor effect.** Quantify how much of the T1→T2 site-count drop (1,289 → 407) is biological vs SMRT sensitivity loss by comparing per-site IPD ratios at T2 even where they fall below the official call threshold.
