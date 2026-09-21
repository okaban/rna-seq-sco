#!/usr/bin/env python3
"""Recompute the nearest_GCCGGC_T1/T2/T3 columns of SuppTable1_Exposed62_identity.tsv from the
canonical per-timepoint file (BLOCKER-0). The old T2/T3 columns came from the first-appearance
table (T2 median 57,884 bp; T3 196,230 bp). The old file is kept as *.superseded_260921."""
import csv, bisect, shutil, collections
from pathlib import Path
A=Path(__file__).resolve().parent.parent
ref="".join(l.strip() for l in open("/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa") if not l.startswith(">")).upper()
def gccggc(p,s): return (ref[p-2:p+4]=="GCCGGC") if s=="+" else (ref[p-3:p+3]=="GCCGGC")
S=collections.defaultdict(list)
for r in csv.DictReader(open(A/"01_integration/high_confidence_sites_weighted.csv")):
    if "4mC" in r["mod_type"]:
        p=int(float(r["position"]))
        if gccggc(p,r["strand"]): S[r["timepoint"]].append(p)
for k in S: S[k].sort()
def near(arr,x):
    i=bisect.bisect_left(arr,x); return min(abs(arr[k]-x) for k in (i-1,i) if 0<=k<len(arr))
f=A/"52_shielded_exposed_boundary/tables/SuppTable1_Exposed62_identity.tsv"
shutil.copy(f, str(f)+".superseded_260921")
rows=list(csv.DictReader(open(f),delimiter="\t")); cols=rows[0].keys()
chg=collections.Counter()
for r in rows:
    tss=int(float(r["tss"]))
    for tp in ("T1","T2","T3"):
        new=near(S[tp],tss); old=float(r[f"nearest_GCCGGC_{tp}"]); r[f"nearest_GCCGGC_{tp}"]=str(new)
        chg[tp]+= (abs(new-old)>0.5)
with open(f,"w",newline="") as fh:
    w=csv.DictWriter(fh,delimiter="\t",fieldnames=list(cols)); w.writeheader(); w.writerows(rows)
import statistics
for tp in ("T1","T2","T3"):
    d=[float(r[f"nearest_GCCGGC_{tp}"]) for r in rows]
    print(tp,"changed",chg[tp],"| within 293:",sum(x<=293 for x in d),"| median",statistics.median(d))
