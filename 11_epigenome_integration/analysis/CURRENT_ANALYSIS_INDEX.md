# 解析ディレクトリの現行版インデックス

**このファイルは生成物。手で書き換えない。** 判断（どれがどれを置き換えたか）は `make_analysis_index.py` の `SUPERSESSION` 定数にある。そこを直して再実行すること。

原稿に入る数値は、下表の「現行」列のディレクトリからのみ取る。「旧」列のディレクトリには `SUPERSEDED.md` が置いてある。

| 量 | 現行 | 旧（使わない） | 旧が誤っている理由 |
|---|---|---|---|
| AAGCCCG 4mC/6mA 共局在の site レベル 2×2・OR・permutation | `87_site_coloc_stats_C4` | `68_or_permutation` | 68_ は 1 bp ずれた `sequence` 窓で共修飾を判定し、universe 1,238,215 を直書きしていた。OR = 138,440 / 244 in / 1,090 not はいずれも再現不能。87_ は universe を明示列挙（1,696,280 候補ペア）。 |
| AAGCCCG per-read 共修飾（Supplementary Table 9） | `79_comod_full_denominator` | `73_comod_threshold_ROC` | 73_ は 2×2 を「≥1 コール行がある read」に条件付けし未修飾 read を構造的に除外していた（Berkson 選択）ため OR が全閾値で 1 未満になった。79_ は全 read を分母にし、さらに 4mC の真の位置 C₄ で走査する（C₃/C₅ 版も同ディレクトリに残すが使わない）。 |
| AAGCCCG 共修飾の site census（ペア数・spacing ラベル） | `79_comod_full_denominator` | `65_per_read_comod`<br>`23_expanded_motif_search` | いずれも `07_motif_analysis/methylation_site_sequences.csv` の `sequence` 窓経由でオフセットを求めており、ラベルが 1 bp ずれる（C₅/C₃ は実際には C₄）。244（23_ 直書き）と 254（65_）はこの窓に由来する。正: 同一 instance 406 ペア。 |
| AAGCCCG 占有率の時系列・Clark–Evans | `88_occupancy_series_and_CE` | — | 新規（2026-09-20）。旧系列 31.6/7.2/4.3% は分母 1,415 が追跡不能、Clark–Evans R = 1.761 / n_eff 88.8 は生成コードがリポジトリに存在しなかった。 |
| メチル化サイトの developmental dynamism（Supplementary Figure 5 / Table 4） | `70_methylation_dynamics/tables` | — | 閾値は **10 パーセントポイント**（`A5b_*_10pp.tsv`）。0.3 の旧閾値は頻度が % 格納であることを見落としたもので、ほぼ全サイトが dynamic になり主張が空虚だった。同ディレクトリの 0.3 版ファイルは使わない。 |
| promoter メチル化 × LFC の偏相関（共変量の定義） | `80_partial_corr_covariates` | — | canon の r = −0.090 は **core/arm 二値**を共変量とする。oriC 連続距離統制はr = −0.070 (p 0.026)、両方同時は −0.073 (p 0.021)。本文に以前あった p = 0.0028 / 0.0018 は生成コードが存在しなかった。 |
| 遺伝子セットのサイズと入れ子関係 | `86_geneset_counts_B10` | — | 本文の 3,197/459/448 は `62_GO_KEGG_enrichment/tables/*_s0_archive.tsv` 由来の旧集計。現行は 3,106/438/443（n = 7,374）。`*_s0_archive.tsv` は参照しない。 |

## 混同しやすいファイル

- `methylation_site_sequences.csv` — **`position` と `strand` は正しい。`sequence` 列の窓は真の部位より 1 塩基上流に中心が置かれている。**モチーフ内オフセットをこの列から求めてはいけない（C₅/C₃ ラベルずれの原因）。オフセットが要る場合は参照配列に position を当てること（79_/87_/88_ の方式）。
- `80_rebase_dualmod_search` — `80_` が 2 つある（`80_partial_corr_covariates` と `80_rebase_dualmod_search`）。番号は重複しているが別物。改番はリンク切れを招くのでしない。

## 撤回済みトークンの漏れ検査

現行ディレクトリの `.tsv` / `.csv` 出力に旧解析固有の値は無い（`138,440`、`138440`、`1,238,215`、`0.0067`、`A₁↔C₅`、`A₀↔C₃`、`31.6%`、`1.761`、`88.8`、`853/855`、`99.8%`、`0.0028`、`0.0018`、`3,197`、`459/448` を検索）。

---

_生成: 2026-09-20 17:25 · `make_analysis_index.py` · HEAD `03a7bd1 analysis: add the supersession index and the supplem` · SUPERSEDED.md 更新 0 件_
