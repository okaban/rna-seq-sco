#!/usr/bin/env python3
"""Regenerate the two supplementary deliverable tables that live in Writing/ from
their canonical analysis outputs, so they cannot drift from the manuscript again.

Targets:
  Writing/Supplementary_statistical_tests.tsv      (Supplementary Table 6)
  Writing/Supplementary_threshold_sensitivity.tsv  (Supplementary Table 9)

Why this script exists: on 2026-09-20 both files still held retired analyses —
the statistical-tests table carried the co-modification OR 138,440 with the
C5/C3 offset labels, the Clark-Evans R = 1.761, and, in rows T04-T06, the
RETRACTED protection-zone classifier (AUC 0.917, n = 57 exposed) that the
manuscript no longer contains; the threshold sweep was the pre-C4 file with an
AUC column. Both are deliverables, so a stale copy would have shipped.

Superseded copies are moved to Writing/archive_superseded/<name>.<date> rather
than overwritten in place. Run after any change to 79_/87_/88_.
"""
import csv, os, shutil, sys
from datetime import datetime
import pandas as pd

REPO = "/Users/okaban/bioinfo/rna-seq"
WRITING = "/Users/okaban/obsidian/Research/rna-seq/Writing"
A = os.path.join(REPO, "11_epigenome_integration/analysis")
ARCH = os.path.join(WRITING, "archive_superseded")
DATE = datetime.now().strftime("%Y%m%d")

coloc = pd.read_csv(f"{A}/87_site_coloc_stats_C4/tables/site_coloc_2x2_C4.tsv", sep="\t").set_index("scope")
ce    = pd.read_csv(f"{A}/88_occupancy_series_and_CE/clark_evans_C4_census.tsv", sep="\t").set_index("scope")
cemc  = pd.read_csv(f"{A}/88_occupancy_series_and_CE/clark_evans_motif_conditioned_null.tsv", sep="\t").set_index("scope")
sweep = pd.read_csv(f"{A}/79_comod_full_denominator/comod_full_denominator_T1_C4.tsv", sep="\t")
p = coloc.loc["pooled"]; c = ce.loc["pooled_any_timepoint"]; m = cemc.loc["pooled_any_timepoint"]

# ---- Supplementary Table 6: statistical tests -------------------------------
UPDATES = {
 "T01": dict(description="AAGCCCG intra-motif co-modification of same-strand A/C position pairs at 3-4 bp spacing vs statistical independence (corrected C4 offsets; Fisher's exact)",
   test_type="Fisher's exact (one-sided)",
   n=f"contingency: a={int(p.in_comod)} b={int(p.in_notcomod)} c={int(p.out_comod)} d={int(p.out_notcomod)}",
   raw_p="< 1e-300", adj_p="< 1e-300",
   conclusion=(f"Significant; no co-modified pair occurs outside the motif, so the odds ratio is saturated "
               f"(Haldane-corrected OR={p.OR:,.0f} [95% CI {p.OR_CI95_lower:,.0f}-{p.OR_CI95_upper:,.0f}]) and is not an effect size; "
               f"A0-C4 (4 bp, {int(p.in_comod_spacing4_A0C4)} pairs) and A1-C4 (3 bp, {int(p.in_comod_spacing3_A1C4)} pairs) pairing is non-random")),
 "T02": dict(description="AAGCCCG co-modification vs a hypergeometric label-permutation null over all 1,696,280 candidate pairs (10,000 draws, seed 42; corrected C4 offsets)",
   test_type="Permutation (n=10000)",
   n=f"n=10000 permutations; expected in-motif co-modified pairs under H0 = {p.expected_in_comod_H0:.2f}; null max = {int(p.perm_null_max_in_comod)}",
   raw_p="< 0.0001", adj_p="< 0.0001",
   conclusion=(f"Observed {int(p.in_comod)} in-motif co-modified pairs vs a null whose maximum over 10,000 draws is "
               f"{int(p.perm_null_max_in_comod)}; {p.fold_over_expectation:.0f}-fold over expectation")),
 "T03": dict(description="Per-read co-modification of 6mA and 4mC on the same AAGCCCG instance, full read denominator (T1; corrected C4 offset)",
   test_type="Fisher's exact (two-sided), per-read",
   n=f"{int(sweep.iloc[0].total):,} instance x read pairs at probability cutoff >= 0.5",
   raw_p="< 1e-300", adj_p="< 1e-300",
   conclusion=(f"OR={sweep.iloc[0].OR:.2f}; P(6mA | 4mC on the same read)={sweep.iloc[0].P_A_given_C:.2f} vs P(6mA)={sweep.iloc[0].P_A:.2f}; "
               f"{sweep.iloc[0].fold_enrichment:.2f}-fold over independence. Replaces the earlier read-identity permutation "
               f"whose 2x2 was conditioned on reads carrying at least one modification call")),
 "T14": dict(description=f"Clark-Evans nearest-neighbour test for the spatial distribution of the {int(c.n)} AAGCCCG dual-modification events of the corrected C4 census",
   test_type="Clark-Evans (circular chromosome) + motif-conditioned permutation (10,000 draws, seed 42)",
   n=f"n={int(c.n)} events; mean nearest-neighbour {c.mean_nn_bp:,.0f} bp vs {c.expected_nn_bp:,.0f} bp expected",
   raw_p=f"{c.p_two_sided}", adj_p=f"{c.p_two_sided}",
   conclusion=(f"R={c.R} (z={c.z}), i.e. weakly CLUSTERED, not dispersed; against a null drawing the same number of instances "
               f"from the 1,334 AAGCCCG instances present, observed/null = {m.ratio_obs_null} (p = {m.p_two_sided}). "
               f"Supersedes the earlier dispersed result, which was computed on the 1-bp-misaligned census")),
}
WITHDRAW = {
 "T04": "Retracted: protection-zone classifier AUC 0.917 on n=57 exposed genes. The label was derived from the methylation data the classifier used as its predictor (circular); the manuscript reports no AUC and uses 62 Exposed / 989 Shielded / 1,051.",
 "T05": "Retracted with T04 (same classifier, T2 hold-out).",
 "T06": "Retracted with T04 (same classifier, T3 hold-out).",
}

src = os.path.join(WRITING, "Supplementary_statistical_tests.tsv")
rows = list(csv.DictReader(open(src, encoding="utf-8"), delimiter="\t"))
fields = list(rows[0].keys())
kept, dropped, updated = [], [], []
for r in rows:
    tid = r["test_id"]
    if tid in WITHDRAW:
        dropped.append((tid, WITHDRAW[tid])); continue
    if tid in UPDATES:
        r.update({k: v for k, v in UPDATES[tid].items() if k in fields}); updated.append(tid)
    kept.append(r)
missing = [t for t in UPDATES if t not in updated]

# ---- Supplementary Table 9: per-read threshold sweep ------------------------
sw = sweep.rename(columns={"OR": "odds_ratio", "P_A_given_C": "P_6mA_given_4mC", "P_A": "P_6mA",
                           "P_C": "P_4mC", "fold_enrichment": "fold_over_independence"})
sw = sw[["threshold", "n_both", "n_A_only", "n_C_only", "n_neither", "total",
         "odds_ratio", "expected_both", "fold_over_independence", "P_6mA_given_4mC", "P_6mA", "P_4mC", "p_value"]]

def main():
    check = "--check" in sys.argv
    if check:
        print(f"would update {updated}, drop {[d[0] for d in dropped]}, missing {missing}"); return 0
    os.makedirs(ARCH, exist_ok=True)
    for name in ("Supplementary_statistical_tests.tsv", "Supplementary_threshold_sensitivity.tsv"):
        f = os.path.join(WRITING, name)
        dest = os.path.join(ARCH, f"{name}.superseded_{DATE}")
        if os.path.exists(f):
            # never overwrite an archived original with a later regeneration
            shutil.move(f, dest) if not os.path.exists(dest) else os.remove(f)
    with open(os.path.join(WRITING, "Supplementary_statistical_tests.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t"); w.writeheader(); w.writerows(kept)
    sw.to_csv(os.path.join(WRITING, "Supplementary_threshold_sensitivity.tsv"), sep="\t", index=False, float_format="%.6g")
    log = os.path.join(os.path.dirname(os.path.abspath(__file__)), "withdrawn_tests.md")
    with open(log, "w", encoding="utf-8") as fh:
        fh.write(f"# Supplementary Table 6 から恒久的に除外している検定（{DATE} 更新）\n\n")
        fh.write("この一覧は WITHDRAW 定数そのもの。再実行しても消えない。\n"
                 "（このスクリプトは冪等なので、2 回目以降の実行では実際の削除は 0 件になる。）\n\n")
        for tid, why in WITHDRAW.items(): fh.write(f"- **{tid}** — {why}\n")
        fh.write(f"\n今回の実行で実際に削除した行: {', '.join(t for t, _ in dropped) or 'なし（既に除外済み）'}\n")
        fh.write(f"\n更新した行: {', '.join(updated)}\n\n本文からの参照は無い（T-id は本文・Supplementary-Materials-List のどちらにも現れない）。\n")
    print(f"statistical tests: {len(rows)} -> {len(kept)} rows (updated {updated}, dropped {[d[0] for d in dropped]}, missing {missing})")
    print(f"threshold sweep:   {len(sw)} rows, columns {list(sw.columns)}")
    print(f"superseded copies -> {ARCH}")
    return 1 if missing else 0

if __name__ == "__main__":
    sys.exit(main())
