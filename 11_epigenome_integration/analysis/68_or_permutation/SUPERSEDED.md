# SUPERSEDED — 2026-09-21

このディレクトリの出力は原稿には**使わない**。

**代替**: `11_epigenome_integration/analysis/87_site_coloc_stats_C4`

**理由**: 68_ は 1 bp ずれた `sequence` 窓で共修飾を判定し、universe 1,238,215 を直書きしていた。OR = 138,440 / 244 in / 1,090 not はいずれも再現不能。87_ は universe を明示列挙（1,696,280 候補ペア）。

`make_analysis_index.py` が生成した `CURRENT_ANALYSIS_INDEX.md` を参照のこと。
このファイルは生成物なので手で書き換えない。
