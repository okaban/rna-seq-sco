# V-Defense (COG V) × 4mC-hypomethylated genes: composition and GCCGGC promoter scan

Date: 2026-05-04
Author: confirmation analyses on the COG-V Defense × 4mC-hypo enrichment
Script: `11_epigenome_integration/analysis/57_temporal_dynamics_exposed_TF/scripts/v_defense_4mC_hypo_GCCGGC.py`
Outputs: `11_epigenome_integration/analysis/57_temporal_dynamics_exposed_TF/results/`

## Inputs
- `11_epigenome_integration/analysis/01_integration/integrated_methyl_expression_weighted.csv` (n=8,083 genes, 4mC/6mA per timepoint + DESeq2 LFC/padj)
- `12_supplementary_figures/.../tables/gene_COG_classification.tsv` (142 genes in `V - Defense`)
- `05_annotation/.../tables/gene_annotation_basic.tsv` (gene names, products, coords, strand)
- `11_epigenome_integration/analysis/18_tss_analyses/comprehensive_tss_table.csv` (Jeong2016 dRNA-seq TSS with GFF fallback)
- `11_epigenome_integration/data/NC_003888.3.fna` (genome reference)

DMG threshold: |Δ4mC| > 10% in at least one of the three pairwise comparisons (T2vsT1, T3vsT1, T3vsT2). Same threshold as the upstream `14_DMG_functional_enrichment` pipeline.

## Result 1 — composition of the V-Defense × 4mC-hypo set

`v_defense_hypo_genes.tsv` (n = 12)

| gene_id | SCO | product | min Δ4mC (%) | dominant LFC sign |
|---|---|---|---|---|
| SC_RS12890 | SCO2170 | class I SAM-dependent methyltransferase | −95.14 | repressed (T3vsT1 LFC −3.22) |
| SC_RS12500 | SCO2092 | rsmH (16S rRNA m4C MTase) | −87.92 | repressed (T3vsT1 LFC −2.78) |
| SC_RS38275 | – | methyltransferase type 11 | −86.85 | repressed (T3vsT1 LFC −1.27) |
| SC_RS38280 | SCO7213 | SAM-dependent methyltransferase | −86.85 | repressed (T3vsT1 LFC −0.92) |
| SC_RS05450 | SCO0705 | SAM-dependent methyltransferase | −85.53 | induced T1→T2, then declines |
| SC_RS13205 | SCO2235 | type II TA system Phd/YefM antitoxin | −77.56 | induced T1→T2, then repressed |
| SC_RS13305 | SCO2256 | panB | −68.55 | repressed T1→T2, recovers T2→T3 |
| SC_RS26345 | SCO4837 | glycine hydroxymethyltransferase | −65.08 | mild |
| SC_RS39465 | SCO7452 | methyltransferase | −57.64 | induced (T2vsT1 LFC +2.12) |
| SC_RS11875 | SCO1969 | methylated-DNA–[protein]-cysteine S-MTase | −54.61 | strongly induced (T2vsT1 LFC +5.53) |
| SC_RS06105 | SCO0835 | class I SAM-dependent methyltransferase | −52.91 | induced (T3vsT1 LFC +1.72) |
| SC_RS06595 | SCO0929 | SAM-dependent methyltransferase | −51.82 | induced T1→T2 |

Composition of the 12-gene set: **10 SAM/methyltransferase-family genes**, 1 toxin–antitoxin antitoxin (SCO2235), 1 hydroxymethyltransferase (SCO2256/panB). The COG-V signal is dominated by methyltransferases, not by classical anti-phage modules.

### Pgl/BREX and methyl-specific restriction nucleases — **not present**

`v_defense_pgl_restriction_check.tsv` confirms that none of the canonical anti-phage candidates appears in the V-Defense × 4mC-hypo set:

| gene_id | label | in COG-V | 4mC-hypo (any TP) | T3vsT1 LFC |
|---|---|---|---|---|
| SC_RS35330 | pglW (SCO6626) | False | False | −0.90 |
| SC_RS35335 | pglX (SCO6627) | True  | False | −0.75 |
| SC_RS35370 | pglY (SCO6635) | False | False | −0.29 |
| SC_RS35375 | pglZ (SCO6636) | False | False | −0.36 |
| SC_RS28835 | pglX (SCO5331, BREX-2 MTase) | True  | False | −0.65 |
| SC_RS23250 | SCO4213 (restriction endonuclease) | True  | False | +2.15 |
| SC_RS25315 | SCO4631 (mcrA, type IV restriction endonuclease) | True | False | −1.46 |
| SC_RS16410 | SCO2863 (DUF3427) | False | False | −1.67 |

All BREX/Pgl loci and the methyl-directed restriction nucleases show 0 % methylation change in the integrated table — they are not detected as DMGs at all and therefore cannot drive the COG-V enrichment. Note: in the M145 NCBI annotation the canonical BREX/Pgl operon is at SCO6626–6636 (not SCO1444–1447); the SCO1444–1447 region is a chitinase/oxidoreductase/TF cluster and is unrelated to BREX.

**Interpretation.** The COG-V Defense × 4mC-hypo enrichment is driven by SAM-dependent methyltransferases of the V category, not by anti-phage restriction systems. This weakens the "phage-defense gene set is being unmasked" reading of the enrichment.

## Result 2 — GCCGGC presence in promoter (−200 / +50 around TSS)

`v_defense_promoter_GCCGGC.tsv` (n = 12; TSS source: 9 Jeong2016 dRNA-seq, 3 GFF fallback)

| gene_id | SCO | strand | TSS source | GCCGGC in promoter | positions rel. TSS |
|---|---|---|---|---|---|
| SC_RS06105 | SCO0835 | − | GFF | **True** | +17 |
| SC_RS12890 | SCO2170 | + | GFF | **True** | −133 |
| SC_RS39465 | SCO7452 | − | GFF | **True** | −64 |
| SC_RS05450 | SCO0705 | − | Jeong2016 | False | – |
| SC_RS06595 | SCO0929 | − | Jeong2016 | False | – |
| SC_RS11875 | SCO1969 | + | Jeong2016 | False | – |
| SC_RS12500 | SCO2092 | − | Jeong2016 | False | – |
| SC_RS13205 | SCO2235 | + | Jeong2016 | False | – |
| SC_RS13305 | SCO2256 | + | Jeong2016 | False | – |
| SC_RS26345 | SCO4837 | + | Jeong2016 | False | – |
| SC_RS38275 | – | + | GFF | False | – |
| SC_RS38280 | SCO7213 | + | Jeong2016 | False | – |

**3 / 12 (25 %)** V-Defense × 4mC-hypo genes carry a GCCGGC site in the −200 / +50 promoter window. The GFF-fallback hits are biased: the three motif-positive genes are exactly the three with no experimental TSS, so the TSS coordinate has higher uncertainty there; only one of the 9 experimental-TSS genes (none, in fact) has a promoter-region GCCGGC.

**Interpretation.** The GCCGGC site whose 4mC loss drives the methylation change for these V-Defense genes is rarely located in the proximal promoter window. For ≥ 75 % of these genes the responsible GCCGGC must lie in the gene body or further upstream, so a clean "promoter GCCGGC demethylation → derepression" model does not fit this set.

## Files written

```
11_epigenome_integration/analysis/57_temporal_dynamics_exposed_TF/results/
├── v_defense_hypo_genes.tsv               (12 rows; gene composition + Δ4mC + LFC)
├── v_defense_promoter_GCCGGC.tsv          (12 rows; TSS, strand, motif, positions)
└── v_defense_pgl_restriction_check.tsv    (12 canonical anti-phage candidates checked)
```

Symlinked from the analysis directory `results/260504_v_defense_4mC_hypo_GCCGGC_report.md → reports/...`.

## Caveats

- The "OR ≈ 2.38, padj ≈ 0.028" enrichment in the upstream brief most closely matches the COG-V row in `COG_DMG_T3vsT2_4mC_hypo.tsv` (fold-enrichment 2.33, p = 0.0076, padj = 0.129 in the existing per-comparison file). The set used here is the union across all three comparisons (n = 12), which is broader and so includes genes with hypomethylation in any direction. Restricting to T3vsT2 alone would give the 11-gene subset listed in that COG row.
- The promoter window (−200 / +50) is a single defensible default. Wider windows (e.g. −500 / +100) or gene-body inclusion would change the GCCGGC hit rate; rerun with adjusted bounds if needed.
- Three TSS values come from GFF annotation rather than Jeong 2016 dRNA-seq, which biases the proximal-promoter motif scan upward for those three genes.
