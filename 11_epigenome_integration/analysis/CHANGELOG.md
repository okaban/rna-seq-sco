# Analysis CHANGELOG

## 2026-05-06: Exposed TF cohort revision (n=62 → n=57)

### Summary
The "Exposed TF" cohort definition was revised from n=62 to n=57. The corresponding
"Shielded" count was revised from n=993 to n=998 (Total regulators 1055 − Exposed 57 = 998).
This CHANGELOG documents the text-only propagation of the revised counts through scripts,
the manuscript, and historical reports.

### Authoritative source of truth
- `29_genomewide_TF_screen/tables/coordinated_regulatory_genes.tsv` — 57 rows (header excluded)
- `29_genomewide_TF_screen/tables/all_regulatory_genes.tsv` — 1055 rows (header excluded)

### Numerical conversions
| Old value           | New value           | Note                                  |
|---------------------|---------------------|---------------------------------------|
| n=62 (Exposed)      | n=57                | Cohort size                           |
| n=993 (Shielded)    | n=998               | 1055 − 57 = 998                       |
| 5.876% (62/1055)    | 5.403% (57/1055)    | Expected exposed rate                 |
| 94.1% (993/1055)    | 94.6% (998/1055)    | Shielded fraction                     |
| 5.9% (62/1055)      | 5.4% (57/1055)      | Exposed fraction                      |

### Scope of edits — text replacement only
This update is a *text-only* propagation. Python scripts were **not re-executed**.
The `coordinated_regulatory_genes.tsv` is already the n=57 version (re-run by the
preceding session). The figures/tables under analysis directories 29, 40, 51, 52,
54-59 were also re-rendered by the preceding session.

### Files modified — analysis scripts (Python)
The following scripts had hardcoded `62`, `993`, or related label/comment text replaced:

- `31_coordinated_regulators_characterization/scripts/H8_analysis.py`
- `34_gatekeeper_model_v2/scripts/gatekeeper_model_v2.py`
- `50_coordinated_regulators_protection/scripts/H27_protection_zone_analysis.py`
- `51_exposed_regulators_characteristics/scripts/H28_exposed_regulators_analysis.py`
- `52_shielded_exposed_boundary/scripts/H29_shielded_exposed_boundary.py`
- `53_TSS_sequence_determinants/scripts/H30_TSS_sequence_determinants.py`
- `54_exposed_TF_downstream_network/scripts/H31_downstream_target_network.py`
- `55_exposed_regulatory_module/scripts/H32_exposed_regulatory_module.py`
- `56_exposed_TF_neighborhood/scripts/H33_neighborhood_impact.py`
- `57_temporal_dynamics_exposed_TF/scripts/H34_temporal_dynamics.py`
- `58_AAGCCCG_exposed_TF_causal/scripts/H35_AAGCCCG_causal_pathway.py`
- `58b_exposed_TF_functional_prediction/scripts/H35_functional_prediction.py`
- `59_exposed_TF_conservation/scripts/H36_conservation_analysis.py`
- `60_coexpression_regulon_prediction/scripts/H37_coexpression_regulon.py`

Patterns replaced inside scripts:
- `n=62`, `n=993` (plot labels, docstrings, print statements)
- `n_exp = 62`, `n_exposed = 62`, `n_shi = 993`, `n_shielded = 993`
- `62 / 1055`, `62/1055` → `57 / 1055`, `57/1055`
- `n_shi_other = 993 - n_shi_fam`, `n_shi_poly = 993 - n_shi_mono`,
  `n_shi_dynamic = 993 - n_shi_const - n_shi_nodata`, `denom_shi = 993 - n_shi_nodata`
- `5.876`, `5.876777251184834` → `5.403`, `5.4028436018957345`
- Plot labels like `'Exposed\n(n=62)'`, `'Shielded\n(n=993)'`,
  `'Coordinated\n(n=62)'`, `'Coordinated 62'`, `'62\nCoordinated'`
- Docstring/comment text such as
  "62 exposed TFs", "62 coordinated regulators", "62 regulatory genes",
  "62 TF families", "62 'exposed'", "62 vs genome", "62 vs 1,055", "62 non-literature",
  "Random sample of 62", "size=62", "from 62 genes", "Both in 62", "Coordination Types in 62",
  "Layer 3: 62 non-literature", "62 Non-Lit.", "62/62 outside literature 37"

### Files modified — manuscript (`/15_paper_figures/manuscript/`)
- `00_outline.md`
- `01_results.md`
- `02_methods.md`
- `03_discussion.md`
- `04_introduction.md`
- `05_abstract.md`
- `07_style_guide.md`
- `full_manuscript.md`

Patterns replaced (in addition to those above):
- `993/1,055`, `993/1055` → `998/1,055`, `998/1055`
- `94.1%` → `94.6%`, `(5.9%)` → `(5.4%)`
- `62 × 62` → `57 × 57`
- `62 (5.4%)` → `57 (5.4%)`
- Subset-fraction denominators: `0/62`, `2/62`, `9 of 62`, `39/62`, `61/62`,
  `38 of 62`, `~16/62`, `62/62` → `0/57`, `2/57`, `9 of 57`, `39/57`, `61/57`,
  `38 of 57`, `~16/57`, `57/57`
- Phrase replacements: "None of the 62", "the 62 are", "the 62 have",
  "of the 62", `the 62 "`, `size 62`, `993 shielded`, `993 (94`, `993 maintain`,
  `993 regulators`

### Files modified — reports (`/reports/`)
The same set of patterns was propagated to historical hypothesis reports H8, H11,
H27, H28, H29, H30, H31, H32, H33, H34, H35, H36, H37, the gatekeeper synthesis
reports v3 and v4, the TSS-based integration report, and the report `README.md`.

### Caveats / next steps
1. **Subset numerators were not recomputed.** Where a fraction has `62` only as
   the denominator (e.g., "9 of 62" → "9 of 57"), the numerator was left
   unchanged. Some subset counts may shift after re-running the analysis on
   the n=57 cohort. Concrete subset facts to re-verify:
   - "9 of 57" cognate TCS members among exposed
   - "39/57" early responders (phase ratio > 0.6)
   - "61/57" — this fraction is now invalid (>1); module coverage needs
     recomputation on the 57-gene set
   - "0/57" curated FIMO motifs — likely still 0 but should be confirmed
   - "2/57" methylated AAGCCCG at exposed TSS — needs confirmation
   - "~16/57" AAGCCCG sequence at exposed TSS — needs recomputation
   - "38 of 57" themselves targets of known TFs — needs recomputation
2. **Hardcoded constants were updated but scripts were not re-executed.**
   Specifically, in `H28_exposed_regulators_analysis.py` the previous run had
   used hardcoded `62` and `993` even though the input table contained 57
   rows (e.g., `expected_pct=5.876` baked into `TF_family_distribution.tsv`).
   The script source is now consistent with n=57 / n=998, but the existing
   output tables in `51_exposed_regulators_characteristics/tables/` and the
   corresponding figures still reflect the previous mixed-state run. A
   re-execution of H28, H29, H29c, and any downstream plot scripts is
   recommended for fully consistent numerical outputs and labels.
3. **`H29c_pfam_regulatory_gene_selection.py` was not executed.** This is a
   heavy ML script intentionally skipped per user direction. Manual
   re-execution is required to recompute the AUC metric on the new n=57
   cohort.
4. **`coordinated_62_pct` column name** in `H8_analysis.py` was not renamed
   to avoid breaking downstream column-name references; this is a cosmetic
   inconsistency only.
5. **Documenting reports `260504_review_groupA_fig_regeneration_report.md`**
   intentionally retains references to the old `n=62/955` numbers because
   those passages explicitly describe the change ("旧版 ... と異なる").
