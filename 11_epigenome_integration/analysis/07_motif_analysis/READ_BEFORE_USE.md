# READ BEFORE USE — 2026-09-20

**`position` と `strand` は正しい。`sequence` 列の窓は真の部位より 1 塩基上流に中心が置かれている。**モチーフ内オフセットをこの列から求めてはいけない（C₅/C₃ ラベルずれの原因）。オフセットが要る場合は参照配列に position を当てること（79_/87_/88_ の方式）。

`11_epigenome_integration/analysis/CURRENT_ANALYSIS_INDEX.md` を参照のこと。
生成物なので手で書き換えない。
