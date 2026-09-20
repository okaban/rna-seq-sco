# SUPERSEDED — 2026-09-20

このディレクトリの出力は原稿には**使わない**。

**代替**: `11_epigenome_integration/analysis/79_comod_full_denominator`

**理由**: 73_ は 2×2 を「≥1 コール行がある read」に条件付けし未修飾 read を構造的に除外していた（Berkson 選択）ため OR が全閾値で 1 未満になった。79_ は全 read を分母にし、さらに 4mC の真の位置 C₄ で走査する（C₃/C₅ 版も同ディレクトリに残すが使わない）。

`make_analysis_index.py` が生成した `CURRENT_ANALYSIS_INDEX.md` を参照のこと。
このファイルは生成物なので手で書き換えない。
