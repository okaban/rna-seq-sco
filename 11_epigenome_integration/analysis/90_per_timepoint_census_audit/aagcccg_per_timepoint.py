#!/usr/bin/env python3
"""AAGCCCG occupancy and same-instance 4mC(C4)+6mA(A0/A1) pairs per timepoint, from the
CANONICAL per-timepoint file. Replaces the T2/T3 values of 79_/87_/88_, which read the
first-appearance table 07_motif_analysis/methylation_site_sequences.csv (T1 unaffected)."""
import csv, re, collections
from pathlib import Path
H=Path(__file__).resolve().parent; A=H.parent
ref="".join(l.strip() for l in open("/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa") if not l.startswith(">")).upper()
inst=[(m.start(),"+") for m in re.finditer("AAGCCCG",ref)]+[(m.start(),"-") for m in re.finditer("CGGGCTT",ref)]
# instance key -> genomic positions of A0, A1, C4 (0-based, same strand)
def coords(st,s):
    return ((st,st+1,st+4) if s=="+" else (st+6,st+5,st+2))   # '-' : CGGGCTT ; A0 is last base, A1 second last, C4 at offset 2
rows=list(csv.DictReader(open(A/"01_integration/high_confidence_sites_weighted.csv")))
S=collections.defaultdict(set)
for r in rows:
    mod="4mC" if "4mC" in r["mod_type"] else "6mA"
    S[(r["timepoint"],mod)].add((int(float(r["position"])),r["strand"]))
out=[]; pairs_pooled=set()
for tp in ("T1","T2","T3"):
    n4=n6=nb=0; core=0
    for st,s in inst:
        a0,a1,c4=coords(st,s)
        has6=((a0,s) in S[(tp,"6mA")]) or ((a1,s) in S[(tp,"6mA")]); has4=(c4,s) in S[(tp,"4mC")]
        n6+=has6; n4+=has4
        if has6 and has4: nb+=1; pairs_pooled.add((st,s)); core+= 1_500_000<=st<=7_170_000
    N=len(inst); e=n4*n6/N
    out.append(dict(timepoint=tp,n_instances=N,inst_6mA=n6,occupancy_6mA_pct=round(100*n6/N,1),inst_4mC_C4=n4,pct_4mC=round(100*n4/N,1),
                    inst_both=nb,expected_both_indep=round(e,1),ratio_obs_exp=round(nb/e,2) if e else None,core_frac_both=round(core/nb,3) if nb else None))
with open(H/"tables/aagcccg_per_timepoint.tsv","w",newline="") as fh:
    w=csv.DictWriter(fh,delimiter="\t",fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
for o in out: print(o)
print("pooled instances with both marks at the SAME timepoint (union over timepoints):",len(pairs_pooled))
