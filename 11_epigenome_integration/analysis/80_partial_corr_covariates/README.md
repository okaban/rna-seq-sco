# 80 — partial Spearman, position covariates (ledger B04, 2026-09-13)

`partial_corr_covariates.py` -> `tables/partial_corr_covariates.tsv`, `tables/partial_corr_oriC_definition_sensitivity.tsv`

Promoter GCCGGC 4mC occupancy at T1 (sum of site frequencies within TSS ± 2 kb) vs LFC_T2vsT1,
regulatory genes with a defined LFC (n = 1,019; canonical identity table
`52_shielded_exposed_boundary/tables/SuppTable_regulatory_gene_classification_n1051.tsv`).
Partial Spearman = Pearson of rank residuals after OLS on ranked covariate(s); bootstrap 95% CI
(1,000 gene resamples, seed 1).

| covariate | r | 95% CI | p |
|---|---|---|---|
| none (raw Spearman) | −0.127 | — | 4.7e-05 |
| core/arm binary (CANON) | **−0.090** | [−0.152, −0.025] | 0.0041 |
| \|TSS − oriC\| continuous (oriC = 4,271,763 bp, dnaA midpoint; same constant as 83_oriC_confound_FIRE) | −0.070 | [−0.127, −0.010] | 0.026 |
| both jointly | −0.073 | [−0.132, −0.015] | 0.021 |

The canon value is reproduced exactly. The manuscript's oriC-covariate values (r = −0.09, p = 0.0028;
joint r = −0.10, p = 0.0018; EN L85) are NOT reproduced by this definition. Of the alternatives tried,
only log10(|TSS − oriC| + 1) entered linearly comes close for the single covariate (r = −0.094,
p = 0.0026); none reproduces the joint value. Decision on wording is the author's/editor's.
