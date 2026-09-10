Directory: results (canonical, current — strand-specific -s 2)
Created: 2026-04-19 13:02:18

Contents:
  results/
    DESeq2_M145_2_vs_1.tsv       - T2 vs T1 DEGs (apeglm LFC shrinkage)
    DESeq2_M145_3_vs_1.tsv       - T3 vs T1 DEGs
    DESeq2_M145_3_vs_2.tsv       - T3 vs T2 DEGs
    DE_summary_M145.tsv          - Summary table (up/down counts)
    normalized_counts_M145.tsv   - DESeq2 size-factor normalized counts
    comparison_s0_vs_s2.tsv      - s0 vs s2 DEG comparison report
  figures/                        - PCA, volcano, heatmaps (PDF + SVG)
  logs/deseq2_pipeline_s2.log    - Full run log
  rds/dds_raw.rds, rld.rds       - DESeqDataSet objects

Change history:
  2026-04-19 13:02:18 - Replaced unstranded (-s 0) with strand-specific counts (-s 2).
  Reason: featureCounts strandedness correction (RF/fr-firststrand confirmed).
  Script: 04_deseq2/scripts/run_deseq2_M145_s2.R
  Counts: 03_quantification/data/featureCounts/featureCounts_M145_s2.txt

Key statistics (padj<0.05, |LFC|>1):
  T2 vs T1: 3840 DEGs (up: 2038, down: 1802)  [s0: 3848, agreement: 91.7%]
  T3 vs T1: 4674 DEGs (up: 2727, down: 1947)  [s0: 4841, agreement: 88.7%]
  T3 vs T2: 3067 DEGs (up: 1850, down: 1217)  [s0: 3507, agreement: 80.7%]

Archive: ../results_s0_archive/  (old -s 0 results)

NOTE: Downstream scripts in 11_epigenome_integration/ still reference:
  04_deseq2/analysis/04_deseq2_260128_v1/results/
  Those files are preserved; update scripts to use this directory for s2 results.
