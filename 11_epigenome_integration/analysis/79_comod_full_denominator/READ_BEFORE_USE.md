# READ BEFORE USE — 2026-09-21

site census の **T2 (26) / T3 (0) ペア数は first-appearance 表由来で誤り**（正準: T2 248 / T3 171 instances、`90_`）。**per-read OR 4.84 も誤り**: `comod_full_denominator_C4.py` は `modified_bases_forward` を SEQ 座標と混用し minus 鎖を全て未修飾扱いにしていた（EPI-03）。正しい表は `comod_full_denominator_T1_C4_orientation_fixed.tsv`（OR 1.04 → 0.96；`91_` で生成）。T1 の site 値と pooled 406（定義明記）は有効。

`11_epigenome_integration/analysis/CURRENT_ANALYSIS_INDEX.md` を参照のこと。
生成物なので手で書き換えない。
