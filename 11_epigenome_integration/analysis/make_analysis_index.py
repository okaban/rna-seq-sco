#!/usr/bin/env python3
"""Generate CURRENT_ANALYSIS_INDEX.md — which directory is canonical for each
manuscript-facing quantity, and which older ones it superseded.

Same contract as Writing/make_bootstrap.py and Writing/make_data_map.py: the
OUTPUT is generated, never hand-edited. The judgements that cannot be derived
from the filesystem (which analysis replaced which, and why) live in the
SUPERSESSION constant below — edit that, then re-run.

The script is also a tripwire: it verifies every path in the map exists, writes
a SUPERSEDED.md stub into every superseded directory, and greps the CANONICAL
directories for values that belong only to the retired analyses. A retired value
appearing in a canonical directory means the two have been mixed again.

Usage:  python make_analysis_index.py [--check]   (--check = verify only, no writes)
"""
import argparse, os, re, sys, subprocess
from datetime import datetime

H = os.path.dirname(os.path.abspath(__file__))
REPO = "/Users/okaban/bioinfo/rna-seq"

# quantity -> (canonical path, [superseded paths], why the old one is wrong, retired tokens)
SUPERSESSION = [
 dict(quantity="AAGCCCG 4mC/6mA 共局在の site レベル 2×2・OR・permutation",
      canonical="11_epigenome_integration/analysis/87_site_coloc_stats_C4",
      superseded=["11_epigenome_integration/analysis/68_or_permutation"],
      why="68_ は 1 bp ずれた `sequence` 窓で共修飾を判定し、universe 1,238,215 を直書きしていた。"
          "OR = 138,440 / 244 in / 1,090 not はいずれも再現不能。87_ は universe を明示列挙（1,696,280 候補ペア）。",
      retired=["138,440", "138440", "1,238,215"]),
 dict(quantity="AAGCCCG per-read 共修飾（Supplementary Table 9）",
      canonical="11_epigenome_integration/analysis/79_comod_full_denominator",
      superseded=["11_epigenome_integration/analysis/73_comod_threshold_ROC"],
      why="73_ は 2×2 を「≥1 コール行がある read」に条件付けし未修飾 read を構造的に除外していた"
          "（Berkson 選択）ため OR が全閾値で 1 未満になった。79_ は全 read を分母にし、"
          "さらに 4mC の真の位置 C₄ で走査する（C₃/C₅ 版も同ディレクトリに残すが使わない）。",
      retired=["0.0067"]),
 dict(quantity="AAGCCCG 共修飾の site census（ペア数・spacing ラベル）",
      canonical="11_epigenome_integration/analysis/79_comod_full_denominator",
      superseded=["11_epigenome_integration/analysis/65_per_read_comod",
                  "11_epigenome_integration/analysis/23_expanded_motif_search"],
      why="いずれも `07_motif_analysis/methylation_site_sequences.csv` の `sequence` 窓経由で"
          "オフセットを求めており、ラベルが 1 bp ずれる（C₅/C₃ は実際には C₄）。"
          "244（23_ 直書き）と 254（65_）はこの窓に由来する。正: 同一 instance 406 ペア。",
      retired=["A₁↔C₅", "A₀↔C₃"]),
 dict(quantity="AAGCCCG 占有率の時系列・Clark–Evans",
      canonical="11_epigenome_integration/analysis/88_occupancy_series_and_CE",
      superseded=[],
      why="新規（2026-09-20）。旧系列 31.6/7.2/4.3% は分母 1,415 が追跡不能、"
          "Clark–Evans R = 1.761 / n_eff 88.8 は生成コードがリポジトリに存在しなかった。",
      retired=["31.6%", "1.761", "88.8"]),
 dict(quantity="メチル化サイトの developmental dynamism（Supplementary Figure 5 / Table 4）",
      canonical="11_epigenome_integration/analysis/70_methylation_dynamics/tables",
      superseded=[],
      why="閾値は **10 パーセントポイント**（`A5b_*_10pp.tsv`）。0.3 の旧閾値は"
          "頻度が % 格納であることを見落としたもので、ほぼ全サイトが dynamic になり主張が空虚だった。"
          "同ディレクトリの 0.3 版ファイルは使わない。",
      retired=["853/855", "99.8%"]),
 dict(quantity="promoter メチル化 × LFC の偏相関（共変量の定義）",
      canonical="11_epigenome_integration/analysis/80_partial_corr_covariates",
      superseded=[],
      why="canon の r = −0.090 は **core/arm 二値**を共変量とする。oriC 連続距離統制は"
          "r = −0.070 (p 0.026)、両方同時は −0.073 (p 0.021)。"
          "本文に以前あった p = 0.0028 / 0.0018 は生成コードが存在しなかった。",
      retired=["0.0028", "0.0018"]),
 dict(quantity="遺伝子セットのサイズと入れ子関係",
      canonical="11_epigenome_integration/analysis/86_geneset_counts_B10",
      superseded=[],
      why="本文の 3,197/459/448 は `62_GO_KEGG_enrichment/tables/*_s0_archive.tsv` 由来の旧集計。"
          "現行は 3,106/438/443（n = 7,374）。`*_s0_archive.tsv` は参照しない。",
      retired=["3,197", "459/448"]),
]

# Files that are canonical for one purpose but misleading for another.
CAVEATS = [
 dict(path="11_epigenome_integration/analysis/07_motif_analysis/methylation_site_sequences.csv",
      note="**`position` と `strand` は正しい。`sequence` 列の窓は真の部位より 1 塩基上流に中心が置かれている。**"
           "モチーフ内オフセットをこの列から求めてはいけない（C₅/C₃ ラベルずれの原因）。"
           "オフセットが要る場合は参照配列に position を当てること（79_/87_/88_ の方式）。"),
 dict(path="11_epigenome_integration/analysis/80_rebase_dualmod_search",
      note="`80_` が 2 つある（`80_partial_corr_covariates` と `80_rebase_dualmod_search`）。"
           "番号は重複しているが別物。改番はリンク切れを招くのでしない。"),
]

STUB = """# SUPERSEDED — {date}

このディレクトリの出力は原稿には**使わない**。

**代替**: `{canon}`

**理由**: {why}

`make_analysis_index.py` が生成した `CURRENT_ANALYSIS_INDEX.md` を参照のこと。
このファイルは生成物なので手で書き換えない。
"""

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    problems, L = [], []
    W = L.append

    # 1. every path in the map must exist
    for e in SUPERSESSION:
        for p in [e["canonical"], *e["superseded"]]:
            if not os.path.exists(os.path.join(REPO, p)):
                problems.append(f"missing path in map: {p}")
    for c in CAVEATS:
        if not os.path.exists(os.path.join(REPO, c["path"])):
            problems.append(f"missing caveat path: {c['path']}")

    # 2. tripwire: retired tokens inside canonical directories
    leaks = []
    for e in SUPERSESSION:
        cdir = os.path.join(REPO, e["canonical"])
        for tok in e["retired"]:
            hits = []
            for root, _, files in os.walk(cdir):
                for f in files:
                    if not f.endswith((".tsv", ".csv", ".md", ".py")): continue
                    fp = os.path.join(root, f)
                    try: txt = open(fp, encoding="utf-8", errors="ignore").read()
                    except OSError: continue
                    # a retired value quoted inside an explanation is fine; flag only
                    # files that are outputs (tsv/csv), not prose or scripts
                    if tok in txt and f.endswith((".tsv", ".csv")):
                        hits.append(os.path.relpath(fp, REPO))
            if hits: leaks.append((e["canonical"], tok, hits))

    # 3. write SUPERSEDED.md stubs
    written = []
    if not a.check:
        for e in SUPERSESSION:
            for p in e["superseded"]:
                fp = os.path.join(REPO, p, "SUPERSEDED.md")
                body = STUB.format(date=datetime.now().strftime("%Y-%m-%d"),
                                   canon=e["canonical"], why=e["why"])
                old = open(fp, encoding="utf-8").read() if os.path.exists(fp) else None
                if old != body:
                    open(fp, "w", encoding="utf-8").write(body); written.append(p)

    # 3b. caveat markers sit next to the file they warn about
    CAVEAT_STUB = """# READ BEFORE USE — {date}

{note}

`11_epigenome_integration/analysis/CURRENT_ANALYSIS_INDEX.md` を参照のこと。
生成物なので手で書き換えない。
"""
    if not a.check:
        for c in CAVEATS:
            tgt = os.path.join(REPO, c["path"])
            d = tgt if os.path.isdir(tgt) else os.path.dirname(tgt)
            fp = os.path.join(d, "READ_BEFORE_USE.md")
            body = CAVEAT_STUB.format(date=datetime.now().strftime("%Y-%m-%d"), note=c["note"])
            if (open(fp, encoding="utf-8").read() if os.path.exists(fp) else None) != body:
                open(fp, "w", encoding="utf-8").write(body); written.append(c["path"])

    # 4. the index itself
    W(f"# 解析ディレクトリの現行版インデックス\n")
    W("**このファイルは生成物。手で書き換えない。** 判断（どれがどれを置き換えたか）は "
      "`make_analysis_index.py` の `SUPERSESSION` 定数にある。そこを直して再実行すること。\n")
    W("原稿に入る数値は、下表の「現行」列のディレクトリからのみ取る。"
      "「旧」列のディレクトリには `SUPERSEDED.md` が置いてある。\n")
    W("| 量 | 現行 | 旧（使わない） | 旧が誤っている理由 |")
    W("|---|---|---|---|")
    for e in SUPERSESSION:
        rel = lambda q: q.split("analysis/", 1)[-1]
        old = "<br>".join(f"`{rel(p)}`" for p in e["superseded"]) or "—"
        W(f"| {e['quantity']} | `{rel(e['canonical'])}` | {old} | {e['why']} |")
    W("\n## 混同しやすいファイル\n")
    for c in CAVEATS:
        W(f"- `{c['path'].split('analysis/')[-1]}` — {c['note']}")
    W("\n## 撤回済みトークンの漏れ検査\n")
    if leaks:
        W("**現行ディレクトリの出力ファイルに、旧解析固有の値が見つかった。混入の疑い。**\n")
        for canon, tok, hits in leaks:
            W(f"- `{os.path.basename(canon)}` に `{tok}`: " + ", ".join(f"`{h}`" for h in hits[:5]))
    else:
        W("現行ディレクトリの `.tsv` / `.csv` 出力に旧解析固有の値は無い（"
          + "、".join(f"`{t}`" for e in SUPERSESSION for t in e["retired"]) + " を検索）。")
    if problems:
        W("\n## 地図の不整合\n"); [W(f"- {p}") for p in problems]
    try:
        head = subprocess.run(["git","log","--oneline","-1"], cwd=REPO, capture_output=True, text=True).stdout.strip()[:60]
    except Exception:
        head = "?"
    W(f"\n---\n\n_生成: {datetime.now():%Y-%m-%d %H:%M} · `make_analysis_index.py` · "
      f"HEAD `{head}` · SUPERSEDED.md 更新 {len(written)} 件_")

    out = os.path.join(H, "CURRENT_ANALYSIS_INDEX.md")
    txt = "\n".join(L) + "\n"
    if not a.check: open(out, "w", encoding="utf-8").write(txt)
    print(f"{'checked' if a.check else 'written'} {out} / leaks {len(leaks)} / problems {len(problems)} / stubs {len(written)}")
    return 1 if (leaks or problems) else 0

if __name__ == "__main__":
    sys.exit(main())
