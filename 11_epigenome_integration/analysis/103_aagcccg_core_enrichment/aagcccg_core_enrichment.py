#!/usr/bin/env python3
"""103 — is the AAGCCCG methylome core-concentrated, against a motif-conditioned null?

The manuscript reports AAGCCCG core fractions but had no null to call them
enriched: the only null in the repo was GCCGGC's (0.631, from its 32,118 genomic
C positions). A mark can only occur where its motif is, so the honest comparison
is the core fraction of ALL genomic AAGCCCG instances, not of the genome.

Matching convention is taken verbatim from 90_/aagcccg_per_timepoint.py: 0-based
positions, strand-aware, 6mA at A0/A1 and 4mC at C4 of the same instance.
"""
import collections, csv, re
from pathlib import Path
import numpy as np

H = Path(__file__).resolve().parent
A = H.parent
CORE = (1_500_000, 7_170_000)
ref = "".join(l.strip() for l in open("/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa")
              if not l.startswith(">")).upper()
inst = ([(m.start(), "+") for m in re.finditer("AAGCCCG", ref)]
        + [(m.start(), "-") for m in re.finditer("CGGGCTT", ref)])


def coords(st, s):
    return (st, st + 1, st + 4) if s == "+" else (st + 6, st + 5, st + 2)


S = collections.defaultdict(set)
for r in csv.DictReader(open(A / "01_integration/high_confidence_sites_weighted.csv")):
    mod = "4mC" if "4mC" in r["mod_type"] else "6mA"
    S[(r["timepoint"], mod)].add((int(float(r["position"])), r["strand"]))

in_core = np.array([CORE[0] <= st <= CORE[1] for st, _ in inst])
null = in_core.mean()
rng = np.random.default_rng(42)
rows = [dict(scope="motif-conditioned null", timepoint="-", mark="AAGCCCG instances",
             n=len(inst), n_core=int(in_core.sum()), core_fraction=round(float(null), 4),
             fold_vs_null="-", p_two_sided="-")]
for tp in ("T1", "T2", "T3"):
    for mark in ("6mA", "4mC", "both"):
        idx = []
        for k, (st, s) in enumerate(inst):
            a0, a1, c4 = coords(st, s)
            h6 = ((a0, s) in S[(tp, "6mA")]) or ((a1, s) in S[(tp, "6mA")])
            h4 = (c4, s) in S[(tp, "4mC")]
            hit = h6 and h4 if mark == "both" else (h6 if mark == "6mA" else h4)
            if hit:
                idx.append(k)
        n = len(idx)
        if n == 0:
            continue
        obs = in_core[idx].mean()
        draws = np.array([in_core[rng.choice(len(inst), n, replace=False)].mean()
                          for _ in range(10_000)])
        p = 2 * min((draws >= obs).mean(), (draws <= obs).mean())
        rows.append(dict(scope="observed", timepoint=tp, mark=mark, n=n,
                         n_core=int(in_core[idx].sum()), core_fraction=round(float(obs), 4),
                         fold_vs_null=round(float(obs / null), 3),
                         p_two_sided=float(f"{max(p, 1e-4):.4g}")))
with open(H / "tables/aagcccg_core_enrichment.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0]))
    w.writeheader(); w.writerows(rows)
for r in rows:
    print(r)
