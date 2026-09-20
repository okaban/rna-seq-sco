#!/usr/bin/env python3
"""90 — per-timepoint site census from the CANONICAL file, and the first-appearance
decomposition that explains the manuscript's 1,289 -> 407 -> 21 series.

Written by the main agent 2026-09-21 (BLOCKER-0). The manuscript's temporal GCCGGC
numbers came from position-deduplicated tables (07_/23_) that keep only the FIRST
timepoint a site appears in; 'sites at T2' was therefore 'sites NEW at T2'.
"""
import csv, collections, bisect, statistics, itertools
from pathlib import Path
H=Path(__file__).resolve().parent; T=H/"tables"; T.mkdir(exist_ok=True)
A=H.parent
ref="".join(l.strip() for l in open("/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa") if not l.startswith(">")).upper()
L=len(ref); CORE=(1_500_000,7_170_000)
rows=list(csv.DictReader(open(A/"01_integration/high_confidence_sites_weighted.csv")))
def motif(p,s,mod):
    if mod=="4mC":
        if (ref[p-2:p+4]=="GCCGGC") if s=="+" else (ref[p-3:p+3]=="GCCGGC"): return "GCCGGC"
        if (ref[p-4:p+3]=="AAGCCCG") if s=="+" else (ref[p-2:p+5]=="CGGGCTT"): return "AAGCCCG"
    else:
        if s=="+" and (ref[p:p+7]=="AAGCCCG" or ref[p-1:p+6]=="AAGCCCG"): return "AAGCCCG"
        if s=="-" and (ref[p-6:p+1]=="CGGGCTT" or ref[p-5:p+2]=="CGGGCTT"): return "AAGCCCG"
    return "other"
S=collections.defaultdict(set)   # (tp, mod, motif) -> set of (pos,strand)
for r in rows:
    mod="4mC" if "4mC" in r["mod_type"] else ("6mA" if "6mA" in r["mod_type"] else r["mod_type"])
    p=int(float(r["position"])); s=r["strand"]; tp=r["timepoint"]
    S[(tp,mod,"all")].add((p,s)); S[(tp,mod,motif(p,s,mod))].add((p,s))
core=lambda st: round(sum(CORE[0]<=p<=CORE[1] for p,_ in st)/len(st),3) if st else float("nan")
out=[]
for mod in ("4mC","6mA"):
    for mo in ("all","GCCGGC","AAGCCCG","other"):
        for tp in ("T1","T2","T3"):
            st=S[(tp,mod,mo)]
            if st: out.append(dict(mod=mod,motif=mo,timepoint=tp,n_sites=len(st),core_fraction=core(st)))
with open(T/"per_timepoint_census.tsv","w",newline="") as fh:
    w=csv.DictWriter(fh,delimiter="\t",fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
# overlaps + first-appearance decomposition
dec=[]
for mod,mo in (("4mC","GCCGGC"),("4mC","all"),("6mA","AAGCCCG"),("6mA","all")):
    t1,t2,t3=(S[(tp,mod,mo)] for tp in ("T1","T2","T3"))
    j=lambda a,b: round(len(a&b)/len(a|b),3)
    dec.append(dict(mod=mod,motif=mo,n_T1=len(t1),n_T2=len(t2),n_T3=len(t3),
        shared_T1_T2=len(t1&t2),jaccard_T1_T2=j(t1,t2),shared_T2_T3=len(t2&t3),jaccard_T2_T3=j(t2,t3),shared_T1_T3=len(t1&t3),jaccard_T1_T3=j(t1,t3),
        all_three=len(t1&t2&t3),new_at_T2=len(t2-t1),core_new_at_T2=core(t2-t1),new_at_T3=len(t3-t1-t2),core_new_at_T3=core(t3-t1-t2),
        lost_after_T1=len(t1-t2),core_lost_after_T1=core(t1-t2),first_appearance_series=f"{len(t1)}/{len(t2-t1)}/{len(t3-t1-t2)}"))
with open(T/"first_appearance_decomposition.tsv","w",newline="") as fh:
    w=csv.DictWriter(fh,delimiter="\t",fieldnames=list(dec[0])); w.writeheader(); w.writerows(dec)
# Exposed-62 promoter status per timepoint
ex=list(csv.DictReader(open(A/"52_shielded_exposed_boundary/tables/SuppTable1_Exposed62_identity.tsv"),delimiter="\t"))
exo=[]
for tp in ("T1","T2","T3"):
    pos=sorted(p for p,_ in S[(tp,"4mC","GCCGGC")])
    def near(x):
        i=bisect.bisect_left(pos,x); c=[abs(pos[k]-x) for k in (i-1,i) if 0<=k<len(pos)]; return min(c)
    d=[near(int(float(e["tss"]))) for e in ex]
    exo.append(dict(timepoint=tp,n_exposed=len(ex),within_293bp=sum(x<=293 for x in d),median_nearest_bp=statistics.median(d),
                    manuscript_table_within_293=sum(float(e[f"nearest_GCCGGC_{tp}"])<=293 for e in ex)))
with open(T/"exposed62_per_timepoint.tsv","w",newline="") as fh:
    w=csv.DictWriter(fh,delimiter="\t",fieldnames=list(exo[0])); w.writeheader(); w.writerows(exo)
for d_ in dec: print({k:d_[k] for k in ("mod","motif","n_T1","n_T2","n_T3","jaccard_T1_T2","first_appearance_series","core_new_at_T2")})
for e in exo: print(e)
