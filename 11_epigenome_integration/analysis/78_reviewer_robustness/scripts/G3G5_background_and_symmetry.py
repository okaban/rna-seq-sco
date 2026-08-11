#!/usr/bin/env python3
"""GAP-3 (non-motif m4C background = empirical FDR proxy) and GAP-5 (GCCGGC
palindrome strand symmetry) computed directly from modkit pileup bedMethyl.

bedMethyl (18 col): chrom start end mod score strand ... Nvalid(10) frac%(11)
Nmod(12) Ncanon(13) ...   mod codes: m=5mC, 21839=4mC, a=6mA.

GCCGGC is palindromic; modified C at +strand offset s+2, -strand partner s+3.
"""
import sys
import subprocess
from pathlib import Path

GENOME = "/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa"
PILEUP = sys.argv[1]            # e.g. .../pileup/1-3_pileup.bed
DEPTH = 10
FREQ = 50.0                     # canonical: depth>=10 & freq>=50%

def load_genome(path: str) -> str:
    seq = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    return "".join(seq)

g = load_genome(GENOME)
N = len(g)

# --- motif footprints (exclude from background) + GCCGGC modified-C sets ---
motif_mask = bytearray(N)                 # 1 if position inside any motif footprint
gccggc_plus = set()                       # modified C on + strand (offset +2)
gccggc_minus = set()                      # modified C on - strand (offset +3)
pairs = []                                # (plus_pos, minus_pos) per occurrence

i = 0
while i <= N - 6:
    if g[i:i+6] == "GCCGGC":
        for j in range(i, i+6):
            motif_mask[j] = 1
        gccggc_plus.add(i+2)
        gccggc_minus.add(i+3)
        pairs.append((i+2, i+3))
    i += 1

# AAGCCCG (asymmetric) + reverse complement CGGGCTT footprints -> exclude too
for mtf in ("AAGCCCG", "CGGGCTT"):
    L = len(mtf)
    start = 0
    while True:
        k = g.find(mtf, start)
        if k < 0:
            break
        for j in range(k, k+L):
            motif_mask[j] = 1
        start = k + 1

print(f"[motif] GCCGGC occurrences={len(pairs)}  "
      f"AAGCCCG+rc footprints masked; total masked bp={sum(motif_mask)}",
      file=sys.stderr)

# --- stream 4mC rows (pre-filtered by awk to mod==21839 & depth>=DEPTH) ---
# GAP-3 counters (multi-threshold to confirm the background is genuinely near-zero)
bg_tot = 0
bg_meth50 = bg_meth25 = bg_meth10 = 0
bg_maxfrac = 0.0
mt_tot = mt_meth = 0          # GCCGGC modified-C cytosines
# GAP-5: frac at modified-C positions, keyed by (pos,strand)
modC_frac = {}               # (pos,strand) -> frac

awk = subprocess.Popen(
    ["awk", "-F\t", f'$4=="21839" && $10>={DEPTH}', PILEUP],
    stdout=subprocess.PIPE, text=True, bufsize=1 << 20)

for line in awk.stdout:
    f = line.split("\t")
    pos = int(f[1]); strand = f[5]; frac = float(f[10])
    is_meth = frac >= FREQ
    # GAP-5 capture: only modified-C positions of GCCGGC
    if strand == "+" and pos in gccggc_plus:
        modC_frac[(pos, "+")] = frac
        mt_tot += 1; mt_meth += is_meth
    elif strand == "-" and pos in gccggc_minus:
        modC_frac[(pos, "-")] = frac
        mt_tot += 1; mt_meth += is_meth
    else:
        # background: cytosine not inside any motif footprint
        if not motif_mask[pos]:
            bg_tot += 1
            bg_meth50 += frac >= 50.0
            bg_meth25 += frac >= 25.0
            bg_meth10 += frac >= 10.0
            if frac > bg_maxfrac:
                bg_maxfrac = frac
awk.wait()

# --- GAP-5 symmetry on paired occurrences (both strands measured) ---
import statistics as st
dp = []; dm = []
for pp, mm in pairs:
    a = modC_frac.get((pp, "+")); b = modC_frac.get((mm, "-"))
    if a is not None and b is not None:
        dp.append(a); dm.append(b)

def pearson(x, y):
    n = len(x); mx = sum(x)/n; my = sum(y)/n
    sxy = sum((xi-mx)*(yi-my) for xi, yi in zip(x, y))
    sxx = sum((xi-mx)**2 for xi in x); syy = sum((yi-my)**2 for yi in y)
    return sxy/((sxx*syy)**0.5) if sxx and syy else float("nan")

print("\n=== GAP-3  non-motif m4C background (empirical FDR proxy) ===")
print(f"pileup: {Path(PILEUP).name}  (depth>={DEPTH})")
print(f"non-motif cytosines:    n={bg_tot:>10,}")
print(f"  called m4C freq>=50%: {bg_meth50:>8,}  ({100*bg_meth50/bg_tot:.5f}%)")
print(f"  called m4C freq>=25%: {bg_meth25:>8,}  ({100*bg_meth25/bg_tot:.5f}%)")
print(f"  called m4C freq>=10%: {bg_meth10:>8,}  ({100*bg_meth10/bg_tot:.5f}%)")
print(f"  max background frac observed: {bg_maxfrac:.1f}%")
print(f"GCCGGC modified-C cyt.: n={mt_tot:>10,}  methylated(>=50%)={mt_meth:>8,}  "
      f"rate={100*mt_meth/mt_tot:.2f}%" if mt_tot else "no motif Cs")

print("\n=== GAP-5  GCCGGC palindrome strand symmetry ===")
print(f"occurrences with both strands at depth>={DEPTH}: n={len(dp):,}")
# strand-bias QC across ALL paired sites (artifact check: should be balanced)
if dp:
    print(f"strand-bias QC (all paired): mean + frac={sum(dp)/len(dp):.2f}%  "
          f"mean - frac={sum(dm)/len(dm):.2f}%  (balanced => no calling artifact)")
# biologically meaningful: among ACTIVE sites (>=50% on at least one strand)
act = [(a, b) for a, b in zip(dp, dm) if max(a, b) >= 50.0]
if act:
    ap = [a for a, b in act]; bp = [b for a, b in act]
    full = sum(1 for a, b in act if a >= 50 and b >= 50)
    hemi = sum(1 for a, b in act if (a >= 50) ^ (b >= 50))
    print(f"\nactive occurrences (>=50% on >=1 strand): n={len(act):,}")
    print(f"  full (both strands >=50%): {full:,}  ({100*full/len(act):.1f}%)")
    print(f"  hemi (one strand >=50%):   {hemi:,}  ({100*hemi/len(act):.1f}%)")
    print(f"  mean + frac={sum(ap)/len(ap):.1f}%  mean - frac={sum(bp)/len(bp):.1f}%")
    print(f"  Pearson r(+ vs -) among active = {pearson(ap, bp):.3f}")

# --- dump small data files for figure (only when OUTPREFIX env given) ---
import os
outpref = os.environ.get("OUTPREFIX")
if outpref:
    with open(outpref + "_active_pairs.tsv", "w") as fh:
        fh.write("plus_frac\tminus_frac\n")
        for a, b in zip(dp, dm):
            fh.write(f"{a}\t{b}\n")
    print(f"[wrote] {outpref}_active_pairs.tsv ({len(dp)} paired GCCGGC occurrences)",
          file=sys.stderr)
