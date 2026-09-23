#!/usr/bin/env python3
"""104 — does the fixed per-site methylated strand track the replication fork?

The Discussion leaves two mechanisms open for the per-site strand preference
(inter-strand r = -0.92, cross-library concordance 98.95%): an asymmetric
recognition/modification event, or something replication-coupled. The two make
different predictions. Replication coupling ties the modified strand to which
fork copies the locus, so the dominant strand should FLIP across oriC. An
asymmetric recognition event does not care about the fork.

Second test, scoped to the 1,051-gene regulatory set (the only gene table in the
repo carrying coordinates and strand): is the dominant strand the gene's sense
or antisense strand? A transcription-linked mechanism predicts a bias.
"""
import csv, gzip, sys
from pathlib import Path
from scipy.stats import binomtest, fisher_exact

H = Path(__file__).resolve().parent
A = H.parent
ORIC = int(sys.argv[1]) if len(sys.argv) > 1 else 4_270_000

dom = {}
with gzip.open(A / "91_comod_stratified_and_hemi/tables/REFA06_GCCGGC_per_instance_strand_pct.tsv.gz", "rt") as fh:
    for d in csv.DictReader(fh, delimiter="\t"):
        p, m = float(d["plus_pct"]), float(d["minus_pct"])
        if max(p, m) < 50:
            continue
        dom.setdefault(int(d["start"]), []).append("+" if p >= m else "-")
# one call per instance: the strand it carries in a majority of libraries
inst = {}
for pos, v in dom.items():
    if len(v) >= 3 and len(set(v)) == 1:
        inst[pos] = v[0]

rows = []
# --- test 1: does the dominant strand flip across oriC? --------------------
left = [(p, s) for p, s in inst.items() if p < ORIC]
right = [(p, s) for p, s in inst.items() if p >= ORIC]
for lab, sub in (("left of oriC", left), ("right of oriC", right)):
    k = sum(1 for _, s in sub if s == "+")
    b = binomtest(k, len(sub), 0.5)
    rows.append(dict(test="plus-strand share by replichore", stratum=lab, n=len(sub),
                     plus=k, frac_plus=round(k / len(sub), 4),
                     p_two_sided=float(f"{b.pvalue:.3g}")))
tab = [[sum(1 for _, s in left if s == "+"), sum(1 for _, s in left if s == "-")],
       [sum(1 for _, s in right if s == "+"), sum(1 for _, s in right if s == "-")]]
orr, pf = fisher_exact(tab)
rows.append(dict(test="replichore x dominant strand", stratum="left vs right of oriC",
                 n=len(inst), plus=tab[0][0] + tab[1][0],
                 frac_plus=round((tab[0][0] + tab[1][0]) / len(inst), 4),
                 p_two_sided=float(f"{pf:.3g}"), odds_ratio=round(float(orr), 3)))

# --- test 2: sense vs antisense, regulatory gene set ------------------------
genes = []
for g in csv.DictReader(open(A / "52_shielded_exposed_boundary/tables/all_genes_features_unified_n57.tsv"),
                        delimiter="\t"):
    try:
        genes.append((int(float(g["start"])), g["strand"]))
    except (ValueError, KeyError):
        continue
genes.sort()
import bisect
starts = [g[0] for g in genes]
sense = anti = 0
for pos, s in inst.items():
    j = min(bisect.bisect_left(starts, pos), len(genes) - 1)
    for cand in {max(0, j - 1), j}:
        pass
    j = min(range(max(0, j - 1), min(len(genes), j + 2)), key=lambda k: abs(starts[k] - pos))
    if abs(starts[j] - pos) > 5000:
        continue
    sense += (s == genes[j][1]); anti += (s != genes[j][1])
b = binomtest(sense, sense + anti, 0.5)
rows.append(dict(test="sense vs antisense (regulatory genes, |d|<=5 kb)",
                 stratum="nearest regulatory gene", n=sense + anti, plus=sense,
                 frac_plus=round(sense / (sense + anti), 4) if sense + anti else None,
                 p_two_sided=float(f"{b.pvalue:.3g}")))

with open(H / "tables/hemi_mechanism.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t",
                       fieldnames=["test", "stratum", "n", "plus", "frac_plus", "p_two_sided", "odds_ratio"])
    w.writeheader()
    for r in rows:
        w.writerow({k: r.get(k, "") for k in w.fieldnames})
print(f"instances with a consistent dominant strand in >=3 libraries: {len(inst)}  (oriC={ORIC:,})")
for r in rows:
    print(" ", r)
