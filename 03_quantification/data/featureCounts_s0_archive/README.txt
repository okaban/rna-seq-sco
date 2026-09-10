Directory: featureCounts_s0_archive (archived, DO NOT USE for new analyses)
Created: 2026-04-19 13:02:04

Contents:
  featureCounts_M145.txt         - Count matrix, unstranded (-s 0)
  featureCounts_M145.txt.summary - Assignment summary

Archive reason:
  2026-04-19 13:02:04 - Archived because strandedness setting was incorrect.
  The library (NEBNext Ultra II Directional, dUTP) is RF/fr-firststrand.
  Unstranded (-s 0) overcounts by including antisense reads; replaced by -s 2.
  Original run: 03_quantification/analysis/03_quant_260128_v1/ (still intact)
  Correct results: ../featureCounts/  (uses -s 2)
