# Annotation Report -- M145 RNA-seq

## 1. Purpose

DESeq2 結果とゲノムアノテーション（GFF）、BGC/regulator 情報を統合したマスターテーブルを作成し、後続の enrichment 解析・BGC 発現プロファイリング・転写制御ネットワーク解析の共通土台を整える。

---

## 2. Input Data

- **GFF**: NCBI RefSeq GCF_000203835.1 `genomic.gff`（8,275 gene/pseudogene, 8,188 CDS records）
- **DESeq2 results**: 3 contrasts (2_vs_1, 3_vs_1, 3_vs_2), 7,646 genes each (after low-count filtering)
- **Normalized counts**: 7,646 genes x 9 samples

---

## 3. Output Tables

| File | Description | Rows | Cols |
|------|-------------|------|------|
| `gene_annotation_basic.tsv` | GFF-derived gene info (gene_id, old_locus_tag, gene_name, product, coordinates, gene_biotype, protein_id, ontology_term) | 8,275 | 11 |
| `gene_master_DESeq2.tsv` | Basic annotation + DESeq2 (LFC, padj x 3 contrasts) + normalized counts | 8,275 | 33 |
| `gene_master_with_BGC.tsv` | Above + bgc_name/bgc_role columns | 8,275 | 35 |
| `gene_master_with_BGC_regulators.tsv` | Above + is_regulator, regulator_name, regulator_type | 8,275 | 38 |
| `BGC_definition_manual.tsv` | Major 4 BGC gene lists (act/red/cda/cpk) | 100 | 5 |
| `BGC_major4_DE_summary.tsv` | DE summary for major 4 BGC genes | 100 | 22 |

---

## 4. Major 4 BGC Definition

| BGC | Product | SCO Range | SC_RS Range | Genes |
|-----|---------|-----------|-------------|-------|
| **act** | actinorhodin | SCO5071--5092 | SC_RS27515--SC_RS27620 | 22 |
| **red** | undecylprodigiosin | SCO5877--5898 | SC_RS31630--SC_RS31735 | 22 |
| **cda** | calcium-dependent antibiotic | SCO3210--3249 | SC_RS18165--SC_RS18360 | 40 |
| **cpk** | coelimycin P1 | SCO6273--6288 | SC_RS33615--SC_RS33690 | 16 |
| | | | **Total** | **100** |

---

## 5. BGC Differential Expression Summary

### 5.1 Overview (M145_3 vs M145_1, padj < 0.05)

| BGC | Total | Sig (padj<0.05) | Up (LFC>1) | Down (LFC<-1) |
|-----|-------|-----------------|------------|----------------|
| act | 22 | **22 (100%)** | **22** | 0 |
| red | 22 | **22 (100%)** | **21** | 0 |
| cda | 40 | **40 (100%)** | **39** | 0 |
| cpk | 16 | **16 (100%)** | **16** | 0 |

All 100 genes in the major 4 BGCs are significantly upregulated in the late growth phase (M145_3 vs M145_1).

### 5.2 Temporal Pattern (all 3 contrasts)

| BGC | 2_vs_1 sig | 2_vs_1 up(LFC>1) | 3_vs_1 sig | 3_vs_1 up(LFC>1) | 3_vs_2 sig | 3_vs_2 up(LFC>1) |
|-----|-----------|-----------------|-----------|-----------------|-----------|-----------------|
| act | 16/22 | 13 | 22/22 | 22 | 22/22 | 22 |
| red | 22/22 | 21 | 22/22 | 21 | 20/22 | 17 |
| cda | 40/40 | 39 | 40/40 | 39 | 9/40 | 0 |
| cpk | 16/16 | 16 | 16/16 | 16 | 10/16 | 0 |

**Key observations**:

- **act**: Progressive activation across timepoints. Most genes already upregulated at M145_2 (13/22 with LFC>1), and **all 22 genes further upregulated** from M145_2 to M145_3 (3_vs_2: 22/22 up). This indicates continuous escalation of actinorhodin biosynthesis from mid to late phase.
- **red**: Rapidly activated by M145_2 (21/22 up), with only **moderate further increase** from M145_2 to M145_3 (17/22 up in 3_vs_2). The red cluster is induced early and plateaus.
- **cda**: Strongly activated by M145_2 (39/40 up), but **no further increase** from M145_2 to M145_3 (0 up in 3_vs_2). The cda cluster reaches maximum expression at M145_2 and stabilizes.
- **cpk**: Strongly activated by M145_2 (16/16 up), but **no further increase** from M145_2 to M145_3 (0 up, 2 down in 3_vs_2). Similar to cda, cpk shows early-mid activation with a plateau.

### 5.3 Top Genes by |LFC| (3_vs_1)

**act** (range LFC +2.6 to +11.3):

| gene_id | SCO | gene_name | LFC | padj |
|---------|-----|-----------|-----|------|
| SC_RS27515 | SCO5071 | -- | +11.34 | 8.5e-27 |
| SC_RS27590 | SCO5086 | -- | +11.11 | 4.5e-162 |
| SC_RS27520 | SCO5072 | -- | +10.89 | 2.5e-158 |
| SC_RS27530 | SCO5074 | -- | +10.87 | 5.4e-153 |
| SC_RS27605 | SCO5089 | -- | +10.17 | 6.0e-31 |

**red** (range LFC +1.0 to +9.8):

| gene_id | SCO | gene_name | LFC | padj |
|---------|-----|-----------|-----|------|
| SC_RS31690 | SCO5889 | -- | +9.83 | 5.7e-13 |
| SC_RS31680 | SCO5887 | -- | +5.81 | 1.5e-55 |
| SC_RS31675 | SCO5886 | -- | +5.80 | 5.8e-87 |
| SC_RS31710 | SCO5893 | -- | +5.42 | 2.6e-48 |
| SC_RS31660 | SCO5883 | -- | +5.42 | 1.1e-36 |

**cda** (range LFC +0.5 to +10.8):

| gene_id | SCO | gene_name | LFC | padj |
|---------|-----|-----------|-----|------|
| SC_RS18350 | SCO3247 | -- | +10.80 | 2.5e-66 |
| SC_RS18315 | SCO3240 | -- | +9.49 | 3.4e-20 |
| SC_RS18225 | SCO3222 | -- | +9.28 | 1.6e-70 |
| SC_RS18325 | SCO3242 | -- | +9.16 | 1.0e-41 |
| SC_RS18295 | SCO3236 | asnO | +9.05 | 1.2e-136 |

**cpk** (range LFC +3.9 to +13.2):

| gene_id | SCO | gene_name | LFC | padj |
|---------|-----|-----------|-----|------|
| SC_RS33660 | SCO6282 | cpkO/kasO | **+13.19** | 4.1e-182 |
| SC_RS33645 | SCO6279 | -- | +9.30 | 1.8e-145 |
| SC_RS33630 | SCO6276 | -- | +8.99 | 2.0e-140 |
| SC_RS33640 | SCO6278 | -- | +8.52 | 3.8e-151 |
| SC_RS33665 | SCO6283 | -- | +7.81 | 4.8e-151 |

SC_RS33660 (cpkO/kasO, SCO6282) has the **highest absolute LFC (+13.19)** among all 100 BGC genes.

---

## 6. Regulators within BGCs

| BGC | gene_id | SCO | Regulator | LFC (3v1) | padj (3v1) | Note |
|-----|---------|-----|-----------|-----------|------------|------|
| act | SC_RS27570 | SCO5082 | actII-orf4 (SARP) | +2.60 | 1.2e-24 | Cluster-situated activator |
| act | SC_RS27585 | SCO5085 | SARP family (auto-detected) | +6.44 | 6.3e-136 | |
| red | SC_RS31630 | SCO5877 | redD (SARP) | +5.00 | 3.5e-65 | Cluster-situated activator |
| red | SC_RS31650 | SCO5881 | redZ (response regulator) | +0.97 | 6.2e-06 | Modest LFC |
| cda | SC_RS18200 | SCO3217 | SARP family (auto-detected) | +5.76 | 1.3e-114 | |
| cda | SC_RS18240 | SCO3225 | absA1 (sensor kinase) | +2.73 | 5.1e-26 | Pleiotropic regulator |
| cda | SC_RS18245 | SCO3226 | absA2 (response reg.) | +0.54 | 1.8e-02 | Modest change |
| cpk | SC_RS33650 | SCO6280 | SARP family (auto-detected) | +5.21 | 2.4e-91 | |
| cpk | SC_RS33660 | SCO6282 | cpkO/kasO (SARP-like) | +13.19 | 4.1e-182 | Strongest BGC regulator |
| cpk | SC_RS33680 | SCO6286 | scbR2 (butyrolactone rec.) | +5.49 | 1.4e-111 | |
| cpk | SC_RS33690 | SCO6288 | SARP family (auto-detected) | +3.90 | 2.1e-32 | |

### Regulator Detection Summary

- **Auto-detected** (product pattern matching): 873 genes genome-wide
- **Manually curated**: 9 known regulators (actII-orf4, redD, redZ, atrA, scbR, scbR2, cpkO/kasO, absA1, absA2)
- **Total**: 877 genes flagged as `is_regulator = TRUE`

---

## 7. Biological Interpretation

### 7.1 All 4 BGCs are induced during growth-phase transition

The most striking finding is that **100% of genes** in all 4 major BGCs are significantly upregulated in M145_3 vs M145_1. This is consistent with the well-known developmental switch in *Streptomyces* from vegetative growth to secondary metabolism.

### 7.2 Distinct temporal activation patterns

The 3-timepoint design reveals **two distinct temporal profiles**:

- **act**: Progressive escalation (M145_1 < M145_2 < M145_3). Actinorhodin biosynthesis continues to increase throughout the observed time window.
- **red / cda / cpk**: Early-mid activation with plateau. These three clusters are strongly induced by M145_2 and show minimal further change from M145_2 to M145_3.

This suggests that act biosynthesis may be regulated by a different or additional set of signals compared to the other three BGCs.

### 7.3 cpkO/kasO as the most dramatically regulated gene

SC_RS33660 (cpkO/kasO, SCO6282) shows LFC = +13.19, the highest of any BGC gene. This SARP-like regulator is known to be a master activator of the cpk cluster. Its extreme dynamic range (near-zero at M145_1 to very high at M145_3) makes it a key marker of the metabolic transition.

### 7.4 absA2 shows a modest response

While absA1 shows clear induction (LFC +2.73), absA2 (the response regulator component) shows only +0.54. Since AbsA1/AbsA2 is a two-component system known to negatively regulate antibiotic production in some contexts, this differential response merits further investigation.

---

## 8. Key Points

1. **All 4 major BGCs (act, red, cda, cpk) are uniformly upregulated** in the late growth phase, with 100/100 genes significant at padj < 0.05 (M145_3 vs M145_1).
2. **act shows a distinct progressive activation pattern** (continuing to increase from M145_2 to M145_3), while red/cda/cpk plateau by M145_2.
3. **cpkO/kasO (SC_RS33660)** is the most dramatically induced regulator (LFC +13.19), suggesting it transitions from silenced to fully active.
4. **Cluster-situated SARP regulators** (actII-orf4, redD, cpkO) are consistently co-induced with their respective clusters, while **pleiotropic regulators** (absA1/absA2) show variable responses.
5. **877 putative regulators** have been flagged genome-wide, providing the foundation for transcription factor-BGC correlation analysis in subsequent steps.

---

## 9. Next Steps

- **Functional enrichment**: Use `gene_master_with_BGC_regulators.tsv` for GO/KEGG/COG enrichment analysis of DEGs
- **BGC expression profiling**: Use `BGC_major4_DE_summary.tsv` to build per-cluster expression heatmaps across timepoints
- **Regulator-BGC network**: Correlate `is_regulator = TRUE` genes with BGC expression patterns to identify candidate transcriptional activators/repressors

---

## 10. Materials & Methods

Gene-level annotations were extracted from the NCBI RefSeq GFF (GCF_000203835.1, 8,275 genes) and merged with DESeq2 differential expression results (3 pairwise contrasts, 7,646 tested genes) and normalized counts. The four major biosynthetic gene clusters of *S. coelicolor* A3(2) -- actinorhodin (act, 22 genes), undecylprodigiosin (red, 22 genes), calcium-dependent antibiotic (cda, 40 genes), and coelimycin P1 (cpk, 16 genes) -- were manually defined based on literature SCO ranges and mapped to RefSeq locus_tags via the GFF `old_locus_tag` attribute. Putative transcriptional regulators were identified by keyword matching against CDS product annotations (873 genes) and supplemented with 9 manually curated regulators of known function. All annotations were integrated into a single master table for downstream enrichment and network analyses.

---

*Generated: 2026-01-28*
*Run directory: `05_annotation_260128_v1`*
