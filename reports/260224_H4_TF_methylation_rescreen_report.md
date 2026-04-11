# H4: Full Re-screen of 37 TFs for Methylation-Expression Coordination

**Date:** 2026-02-24  
**Analysis:** `11_epigenome_integration/analysis/27_TF_methylation_rescreen/`  
**Script:** `run_tf_rescreen.py`

## Objective

Re-screen all 37 literature-curated TFs for methylation-expression coordination using CORRECTED locus_tags, with special focus on SARP upstream regulators controlling the four major BGCs (Red, Act, CDA, CPK).

## Methods

- **TF set:** 37 TFs from `literature_tf_master.csv` (13 global Tier-1, 8 pleiotropic Tier-2, 6 CSR Tier-3, 10 sigma Tier-2)
- **Methylation data:** integrated_methyl_expression_weighted.csv (8,083 genes), high_confidence_sites_weighted.csv (12,429 sites)
- **Expression data:** DESeq2 results for T2vsT1, T3vsT1, T3vsT2
- **Cross-validation:** All 37 locus_tags confirmed present in GFF and integrated dataset
- **Thresholds:** |log2FC| >= 1.0, padj < 0.05 for DEG; methylation change = gain/loss of sites
- **Promoter analysis:** AAGCCCG motif scan in -300 to +50 bp from TSS; HC site overlap

## Key Findings

### 1. Methylation is Rare Among the 37 TFs

Only **3/37 TFs (8%)** carry any high-confidence methylation at any timepoint:

| TF | Locus Tag | 6mA (T1/T2/T3) | 4mC (T1/T2/T3) | DEG? |
|----|-----------|-----------------|-----------------|------|
| afsR | SC_RS24295 | 1/1/1 | 0/0/0 | No (padj>0.05) |
| bldB | SC_RS30830 | 0/0/0 | 1/1/1 | Yes (up T2vsT1) |
| sigF | SC_RS22295 | 0/0/0 | 1/1/1 | Yes (up T3vsT1) |

All three methylated TFs show **stable** methylation (same site count across timepoints), meaning there is no temporal methylation gain/loss that could coordinate with expression changes.

### 2. Coordination Is Absent

Because the 3 methylated TFs show no methylation change across timepoints, **zero TFs** show concordant or discordant methylation-expression coordination:

| Comparison | Concordant | Discordant | Meth-only | Expr-only | No change |
|------------|-----------|------------|-----------|-----------|-----------|
| T2 vs T1   | 0         | 0          | 0         | 28        | 9         |
| T3 vs T1   | 0         | 0          | 0         | 34        | 3         |
| T3 vs T2   | 0         | 0          | 0         | 29        | 8         |

The dominant pattern is **expression-only** change (28-34/37 TFs are DEGs but lack methylation changes).

### 3. SARP Regulators: Strongly Activated but Unmethylated

All four SARP CSRs show dramatic transcriptional activation but carry **zero methylation sites**:

| BGC | SARP | LFC T2v1 | LFC T3v1 | Methylation |
|-----|------|----------|----------|-------------|
| Red | redD (SC_RS31630) | **+5.79** | **+5.00** | None |
| Act | actII-ORF4 (SC_RS27585) | **+1.99** | **+6.44** | None |
| CDA | cdaR (SC_RS18200) | **+6.85** | **+5.76** | None |
| CPK | cpkO/kasO (SC_RS33650) | **+7.18** | **+5.21** | None |

Among upstream regulators, only **afsR** (SC_RS24295) carries methylation (1 stable 6mA site across all timepoints). afsR itself is not a DEG (|LFC| < 1.0 in all comparisons), suggesting its constitutive 6mA mark does not modulate its expression during secondary metabolism onset.

### 4. AAGCCCG Motif: Single Hit at bldB

Only **1/37 TFs** (bldB, SC_RS30830) has the AAGCCCG motif in its promoter region at position 6,244,200. bldB also has a HC 4mC site in its promoter. bldB is upregulated at T2 (LFC +3.07), consistent with its role in morphological differentiation, but its 4mC is stable across timepoints.

### 5. Promoter Methylation

Only **3/37 TFs** have any HC methylation site within their promoter (-300 to +50 bp):

| TF | Promoter HC Sites | AAGCCCG? | Expression |
|----|-------------------|----------|------------|
| bldB | 1 | Yes | Up at T2 |
| afsR | 1 | No | Stable |
| sigF | 1 | No | Up at T3 |

### 6. TFs with Both Methylation and Differential Expression

Only **2 TFs** have both methylation marks and significant expression changes:
- **bldB** (SC_RS30830): 4mC stable 1/1/1, up at T2vsT1 -- expr_only coordination (methylation does not change)
- **sigF** (SC_RS22295): 4mC stable 1/1/1, up at T3vsT1 -- expr_only coordination (methylation does not change)

## Interpretation

The re-screen with corrected locus_tags confirms that **DNA methylation does not play a direct regulatory role at the 37 key TF loci** in *S. coelicolor* M145 during the developmental transition:

1. **92% of TFs lack any methylation** -- the regulatory TF gene set is largely methylation-free
2. The 3 methylated TFs (afsR, bldB, sigF) have **constitutive marks** that do not change temporally
3. All 4 SARP pathway activators are **strongly upregulated (LFC +2 to +7) without any methylation** involvement
4. The massive transcriptional reprogramming of the TF network (28-34 DEGs per comparison) is driven by **transcription factor cascades, not epigenetic modification**

## Output Files

### Tables
- `tables/TF_methylation_expression_corrected.tsv` -- Full 37-TF methylation + expression data
- `tables/TF_coordination_summary.tsv` -- Coordination status per comparison
- `tables/SARP_upstream_methylation.tsv` -- SARP cascade detailed profiles (16 rows)

### Figures
- `figures/fig1_TF_methylation_expression_heatmap.png` -- 37 TF heatmap (methylation + LFC)
- `figures/fig2_SARP_cascade_methylation.png` -- SARP cascade bar charts per BGC
- `figures/fig3_coordination_summary.png` -- Coordination category distribution
- `figures/fig4_promoter_AAGCCCG_motif.png` -- Promoter HC sites + motif presence
