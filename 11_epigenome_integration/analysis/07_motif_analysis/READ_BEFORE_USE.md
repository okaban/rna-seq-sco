# READ BEFORE USE — 2026-09-21

**この表は first-appearance 表である（2026-09-21 判明、BLOCKER-0）**: position で重複除去され、各部位は最初に現れた時点の行しか持たない（4mC 1,987 / 664 / 42）。**T1 の集計にしか使えない。**T2/T3 の部位数・占有率・ペア数・Jaccard をここから出すと「T2 で新たに現れた部位」を数えることになる（本文の 1,289→407→21、Jaccard 0、Exposed 62→0、27.4→6.8→3.7% の原因）。時点別の値は `01_integration/high_confidence_sites_weighted.csv`（全時点・全部位）から出すこと（`90_per_timepoint_census_audit/`）。さらに `sequence` 列の窓は真の部位より 1 塩基上流に中心がある（C₅/C₃ ラベルずれの原因）。オフセットは参照配列に position を当てて求めること。

`11_epigenome_integration/analysis/CURRENT_ANALYSIS_INDEX.md` を参照のこと。
生成物なので手で書き換えない。
