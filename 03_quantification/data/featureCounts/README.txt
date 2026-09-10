Directory: featureCounts (canonical, current)
Created: 2026-04-19 13:02:04

Contents:
  featureCounts_M145_s2.txt         - Count matrix, 8275 genes x 9 samples
  featureCounts_M145_s2.txt.summary - Assignment summary

Change history:
  2026-04-19 13:02:04 - Replaced unstranded (-s 0) with strand-specific (-s 2, RF/fr-firststrand)
  Reason: Library is NEBNext Ultra II Directional (dUTP method), confirmed RF by
          infer_experiment.py (90%). Previous -s 0 caused ~2x count overestimation
          in 17-25% of DEGs due to antisense reads being counted.
  Parameters: featureCounts -T 6 -p --countReadPairs -B -C -s 2 -t gene -g gene_id
  Average assignment rate: 78.6% (vs 89.3% for -s 0 unstranded)

Archive: ../featureCounts_s0_archive/  (old -s 0 results)
