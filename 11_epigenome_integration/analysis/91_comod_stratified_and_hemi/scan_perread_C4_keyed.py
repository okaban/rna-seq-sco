#!/usr/bin/env python3
"""AAGCCCG per-read 4mC(C4)/6mA(A0/A1) scan, T1, keeping instance and replicate keys.

Same unit and same call extraction as 79_comod_full_denominator/comod_full_denominator_C4.py
(full read denominator; a (instance x read) pair is kept if the read covers >=1 A of A0/A1
and C4 of the same instance; per pair the max ML probability of 6mA over covered A's and of
4mC at C4). Difference: every pair is written out with (replicate, strand, instance_start,
read_name, pA, pC) so that downstream stats can stratify by instance and replicate and
build a within-read cross-instance null. Output: tables/perread_pairs_T1_C4.tsv.gz
"""
import re, sys, gzip, bisect
import pysam

BAMS = {s: f"/Users/okaban/bioinfo/methyl/260102_M145/analysis/{s}_mapped.bam" for s in ("1-1","1-2","1-3")}
REF  = "/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa"; CHROM = "NC_003888.3"
MOTIF = "AAGCCCG"; A_OFF = (0,1); C_OFF = (4,)
OUT = sys.argv[1] if len(sys.argv) > 1 else "tables/perread_pairs_T1_C4.tsv.gz"

seq = "".join(l.strip() for l in open(REF) if not l.startswith(">")).upper()
rc  = MOTIF[::-1].translate(str.maketrans("ACGT","TGCA")); L = len(MOTIF)
inst = {"+": [], "-": []}
for m in re.finditer(MOTIF, seq):
    s = m.start(); inst["+"].append((s, [s+a for a in A_OFF], [s+c for c in C_OFF]))
for m in re.finditer(rc, seq):
    s = m.start(); inst["-"].append((s, [s+L-1-a for a in A_OFF], [s+L-1-c for c in C_OFF]))
starts = {k: [t[0] for t in v] for k, v in inst.items()}
print(f"motif instances: + {len(inst['+'])}  - {len(inst['-'])}", file=sys.stderr, flush=True)

n_reads = 0; n_pairs = 0
with gzip.open(OUT, "wt") as out:
    out.write("rep\tstrand\tinst_start\tread\tpA\tpC\tnA_cov\n")
    for rep, bam in BAMS.items():
        with pysam.AlignmentFile(bam) as fh:
            for r in fh.fetch(CHROM):
                if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
                strand = "-" if r.is_reverse else "+"; n_reads += 1
                if n_reads % 50000 == 0: print(f"reads {n_reads} pairs {n_pairs}", file=sys.stderr, flush=True)
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
                    out.write(f"{rep}\t{strand}\t{s}\t{r.query_name}\t{pa:.4f}\t{pc:.4f}\t{len(Acov)}\n"); n_pairs += 1
print(f"reads scanned {n_reads}; (instance,read) pairs covering A and C4: {n_pairs}", file=sys.stderr, flush=True)
