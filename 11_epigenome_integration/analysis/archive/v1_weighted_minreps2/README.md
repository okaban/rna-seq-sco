# Archive: v1 Weighted MIN_REPS=2 Analysis Outputs

**Archived**: 2026-02-24
**Reason**: Replaced by v2 analysis (unweighted, MIN_REPS=3) after T3 resequencing

## v1 Analysis Conditions

| Parameter | Value |
|-----------|-------|
| Consensus method | Coverage-weighted average |
| MIN_REPS | 2 (2/3 replicates required) |
| Sample 3-2 coverage | 14.3x (low) |
| High-confidence sites | 11,778 |

## v2 Analysis Conditions (current)

| Parameter | Value |
|-----------|-------|
| Consensus method | **Unweighted (simple average)** |
| MIN_REPS | **3 (all 3 replicates required)** |
| Sample 3-2 coverage | **57.6x (+303%)** |
| High-confidence sites | **12,429** |

## Key Differences v1 → v2

- 4mC T2vsT1: r=0.172 → r=0.125 (maintained, robust)
- **6mA T2vsT1: r=0.169 → r=-0.038 (LOST, not robust)**
- FDR correction (67 tests): only 1/67 significant in v2

## Archived Contents

| Directory | Contents | Why archived |
|-----------|----------|-------------|
| `01_integration/` | 10 files from initial integration (Jan 29) | Superseded by v2 integration |
| `02_publication_figures/` | Seqlogo figures (Feb 4) | Regenerated with v2 data |
| `04_multiomics_tracks/` | Overview figure (Feb 2) | Regenerated with v2 data |
| `07_motif_analysis/` | MEME results + old SVG (Feb 3-4) | MEME re-run with v2 sites |
| `08_coordinated_enrichment/` | COG enrichment of coordinated genes (Feb 3) | Superseded by 14_DMG_selectivity |
| `18_tss_analyses/` | 13 files from initial TSS analysis (Feb 4) | v2 versions generated |

## Note on Reference Analyses

The following directories contain version-independent reference/comparative analyses
and remain in the main analysis/ directory (not archived):

- `10_atcc_comparative_methylome/` — ATCC strain comparison (REBASE data)
- `11_rm_system_identification/` — R-M system gene identification
- `13_sc_rs17645_analysis/` — SC_RS17645 sequence/structure homology
- `19_bgc_regulator_overview/` — BGC regulator landscape (literature-based)
- `20_atcc_real_methylome/` — REBASE real methylome comparison
- `21_genuswide_motif_conservation/` — Genus-wide motif O/E analysis
- `22_4mC_5mC_competition/` — 4mC/5mC dual modification model
