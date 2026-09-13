#!/usr/bin/env python3
"""AAGCCCG per-read co-modification with the FULL read denominator.

73_comod_threshold_ROC/threshold_sensitivity.tsv built its 2x2 only over reads
that had >=1 modification-call row (modonly files), so 'neither' reads were
structurally excluded and the OR is biased below 1 (Berkson selection). Here
every T1 read that covers both an A position (A0/A1) and a C position (C3/C5)
of the same motif instance on the + strand enters the table, using the MM/ML
tags directly. Same thresholds 0.50-0.90.
"""
import re, sys, csv, math
import pysam
from scipy.stats import fisher_exact

BAMS = [f"/Users/okaban/bioinfo/methyl/260102_M145/analysis/{s}_mapped.bam" for s in ("1-1","1-2","1-3")]
REF  = sys.argv[1]
CHROM = "NC_003888.3"
MOTIF = "AAGCCCG"; A_OFF = (0,1); C_OFF = (3,5)
THRESH = [round(0.5+0.05*i,2) for i in range(9)]
MODS = {("A",0,"a"):"A", ("C",0,21839):"C"}   # (canonical, strand, code) -> arm; strand 0 = same as read seq

seq = "".join(l.strip() for l in open(REF) if not l.startswith(">")).upper()
# motif instances on + strand only (same-strand claim). Reverse strand handled symmetrically below.
fwd = [m.start() for m in re.finditer(MOTIF, seq)]
rc  = MOTIF[::-1].translate(str.maketrans("ACGT","TGCA"))
rev = [m.start() for m in re.finditer(rc, seq)]
inst = {}
for s in fwd:
    for a in A_OFF:
        for c in C_OFF: inst[(s+a, s+c, "+")] = s
for s in rev:  # reverse-strand motif: A positions at end, C positions mirrored
    L=len(MOTIF)
    for a in A_OFF:
        for c in C_OFF: inst[(s+L-1-a, s+L-1-c, "-")] = s
print(f"motif instances: + {len(fwd)}  - {len(rev)}", file=sys.stderr)
A_pos = {k[0] for k in inst}; C_pos = {k[1] for k in inst}

# per (instance, read): max prob of 6mA over its A positions, 4mC over its C positions; None if position not covered
from collections import defaultdict
per_pair = defaultdict(lambda: {"A":None,"C":None,"covA":False,"covC":False})
n_reads=0
for bam in BAMS:
    with pysam.AlignmentFile(bam) as fh:
        for r in fh.fetch(CHROM):
            if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
            strand = "-" if r.is_reverse else "+"
            n_reads+=1
            mb = r.modified_bases_forward or {}
            # reference positions covered by this read
            ap = r.get_aligned_pairs(matches_only=True)
            q2r = dict(ap)
            covered = set(q2r.values())
            # probs at reference positions
            probA = {}; probC = {}
            for key, calls in mb.items():
                canon, st, code = key
                if (canon,code) == ("A","a"): tgt = probA
                elif (canon,code) == ("C",21839): tgt = probC
                else: continue
                for qpos, ml in calls:
                    rp = q2r.get(qpos)
                    if rp is not None: tgt[rp] = max(tgt.get(rp,0), ml/255.0)
            for (ap_, cp_, ms), start in inst.items():
                if ms != strand: continue
                if ap_ in covered and cp_ in covered:
                    d = per_pair[(start, r.query_name)]
                    d["covA"]=d["covC"]=True
                    pa = probA.get(ap_); pc = probC.get(cp_)
                    if pa is not None: d["A"] = max(d["A"] or 0, pa)
                    if pc is not None: d["C"] = max(d["C"] or 0, pc)
print(f"reads scanned {n_reads}; (instance,read) pairs covering both arms: {len(per_pair)}", file=sys.stderr)

rows=[]
for t in THRESH:
    nb=na=nc=nn=0
    for d in per_pair.values():
        a = (d["A"] or 0) >= t; c = (d["C"] or 0) >= t
        if a and c: nb+=1
        elif a: na+=1
        elif c: nc+=1
        else: nn+=1
    OR, p = fisher_exact([[nb,na],[nc,nn]])
    # expected both under independence
    tot=nb+na+nc+nn; exp = (nb+na)*(nb+nc)/tot if tot else float("nan")
    pA_givenC = nb/(nb+nc) if nb+nc else float("nan"); pA = (nb+na)/tot
    rows.append(dict(threshold=t,n_both=nb,n_A_only=na,n_C_only=nc,n_neither=nn,total=tot,OR=round(OR,4),p_value=p,
                     expected_both=round(exp,1),fold_enrichment=round(nb/exp,3) if exp else None,
                     P_A=round(pA,4),P_A_given_C=round(pA_givenC,4)))
with open("comod_full_denominator_T1.tsv","w",newline="") as fh:
    w=csv.DictWriter(fh,delimiter="\t",fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
for r in rows: print(r["threshold"], r["n_both"], r["n_A_only"], r["n_C_only"], r["n_neither"], "OR", r["OR"], "fold", r["fold_enrichment"], "p", f"{r['p_value']:.2e}")
