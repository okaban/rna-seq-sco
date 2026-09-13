#!/usr/bin/env python3
"""AAGCCCG per-read co-modification with the FULL read denominator — CORRECTED C offset.

Derived from comod_full_denominator.py (2026-09-13). The original scanned C
positions at motif offsets 3 and 5; the modkit T1 pileup (1-1_pileup.bed) shows
the 4mC call at AAGCCCG lies at offset 4 (AAGCC*CG; mean 69%, >=50% at 88%/81%
of +/- sites) and is ~0% at offsets 3 and 5. This copy scans offset 4 only and
uses a bisect lookup of motif instances overlapping each read (same logic,
faster). Output: comod_full_denominator_T1_C4.tsv.

Unit: (motif instance x read) pair in which the read covers both A positions
(A0, A1) region and the C4 position of the same instance (at least one A and
the C aligned). Per pair: max ML prob of 6mA over covered A positions, of 4mC
at C4. A pair is A-positive / C-positive if that max >= threshold.
"""
import re, sys, csv, bisect
import pysam
from collections import defaultdict
from scipy.stats import fisher_exact

BAMS = [f"/Users/okaban/bioinfo/methyl/260102_M145/analysis/{s}_mapped.bam" for s in ("1-1","1-2","1-3")]
REF  = sys.argv[1]; CHROM = "NC_003888.3"
MOTIF = "AAGCCCG"; A_OFF = (0,1); C_OFF = (4,)
THRESH = [round(0.5+0.05*i,2) for i in range(9)]

seq = "".join(l.strip() for l in open(REF) if not l.startswith(">")).upper()
rc  = MOTIF[::-1].translate(str.maketrans("ACGT","TGCA")); L = len(MOTIF)
inst = {"+": [], "-": []}   # per strand: list of (start, A positions, C positions) sorted by start
for m in re.finditer(MOTIF, seq):
    s = m.start(); inst["+"].append((s, [s+a for a in A_OFF], [s+c for c in C_OFF]))
for m in re.finditer(rc, seq):
    s = m.start(); inst["-"].append((s, [s+L-1-a for a in A_OFF], [s+L-1-c for c in C_OFF]))
starts = {k: [t[0] for t in v] for k, v in inst.items()}
print(f"motif instances: + {len(inst['+'])}  - {len(inst['-'])}", file=sys.stderr, flush=True)

per_pair = {}   # (strand, start, read) -> [maxA, maxC]
n_reads = 0
for bam in BAMS:
    with pysam.AlignmentFile(bam) as fh:
        for r in fh.fetch(CHROM):
            if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
            strand = "-" if r.is_reverse else "+"; n_reads += 1
            if n_reads % 50000 == 0: print(f"reads {n_reads}", file=sys.stderr, flush=True)
            lo, hi = r.reference_start, r.reference_end
            i0 = bisect.bisect_left(starts[strand], lo - L); i1 = bisect.bisect_right(starts[strand], hi)
            if i0 >= i1: continue
            q2r = dict(r.get_aligned_pairs(matches_only=True)); covered = set(q2r.values())
            probA, probC = {}, {}
            for (canon, st, code), calls in (r.modified_bases_forward or {}).items():
                if (canon, code) == ("A", "a"): tgt = probA
                elif (canon, code) == ("C", 21839): tgt = probC
                else: continue
                for qpos, ml in calls:
                    rp = q2r.get(qpos)
                    if rp is not None: tgt[rp] = max(tgt.get(rp, 0), ml/255.0)
            for s, Apos, Cpos in inst[strand][i0:i1]:
                Acov = [p for p in Apos if p in covered]; Ccov = [p for p in Cpos if p in covered]
                if not Acov or not Ccov: continue
                pa = max((probA.get(p, 0.0) for p in Acov), default=0.0); pc = max((probC.get(p, 0.0) for p in Ccov), default=0.0)
                per_pair[(strand, s, r.query_name)] = (pa, pc)
print(f"reads scanned {n_reads}; (instance,read) pairs covering A and C4: {len(per_pair)}", file=sys.stderr, flush=True)

rows = []
for t in THRESH:
    nb = na = nc = nn = 0
    for pa, pc in per_pair.values():
        a = pa >= t; c = pc >= t
        if a and c: nb += 1
        elif a: na += 1
        elif c: nc += 1
        else: nn += 1
    OR, p = fisher_exact([[nb, na], [nc, nn]])
    tot = nb+na+nc+nn; exp = (nb+na)*(nb+nc)/tot
    rows.append(dict(threshold=t, n_both=nb, n_A_only=na, n_C_only=nc, n_neither=nn, total=tot, OR=round(OR,4), p_value=p,
                     expected_both=round(exp,1), fold_enrichment=round(nb/exp,3) if exp else None,
                     P_A=round((nb+na)/tot,4), P_A_given_C=round(nb/(nb+nc),4) if nb+nc else None, P_C=round((nb+nc)/tot,4)))
with open("comod_full_denominator_T1_C4.tsv","w",newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
for r in rows: print(r["threshold"], r["n_both"], r["n_A_only"], r["n_C_only"], r["n_neither"], "OR", r["OR"], "fold", r["fold_enrichment"], "p", f"{r['p_value']:.2e}")
