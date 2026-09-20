# SUPERSEDED — 2026-09-20

このディレクトリの出力は原稿には**使わない**。

**代替**: `11_epigenome_integration/analysis/79_comod_full_denominator`

**理由**: いずれも `07_motif_analysis/methylation_site_sequences.csv` の `sequence` 窓経由でオフセットを求めており、ラベルが 1 bp ずれる（C₅/C₃ は実際には C₄）。244（23_ 直書き）と 254（65_）はこの窓に由来する。正: 同一 instance 406 ペア。

`make_analysis_index.py` が生成した `CURRENT_ANALYSIS_INDEX.md` を参照のこと。
このファイルは生成物なので手で書き換えない。
