#!/usr/bin/env python3
"""Per-offset mean modification frequency inside GCCGGC and AAGCCCG, from the
modkit pileups — the empirical evidence for WHICH base carries the mark.

Supplementary Figure 10 hard-coded 6mA ~23% at A0/A1 and 4mC ~35% at C3/C5.
The C3/C5 offsets came from the 1-bp-misaligned sequence window; the pileup puts
the AAGCCCG 4mC at C4. This recomputes every offset from the T1 pileups so the
figure and its legend can state measured values at the right positions.

Mean is over every genomic instance of the motif (both strands, offsets mapped
strand-aware), pooling the three T1 replicates by read-weighted average; offsets
with no pileup row contribute 0 (the base is unmodified there, not missing).
"""
import re, glob, numpy as np, pandas as pd
from pathlib import Path
REF="/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa"
PILE="/Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup"
H=Path(__file__).resolve().parent
seq="".join(l.strip() for l in open(REF) if not l.startswith(">")).upper()
MOT={"GCCGGC":"GCCGGC","AAGCCCG":"AAGCCCG"}
RC=str.maketrans("ACGT","TGCA")
inst={}
for name,m in MOT.items():
    rc=m.translate(RC)[::-1]
    plus=[(s.start(),"+") for s in re.finditer(m,seq)]
    minus=[(s.start()+len(m)-1,"-") for s in re.finditer(rc,seq)]
    inst[name]=(plus+minus,len(m))

cols=["chrom","start","end","mod","score","strand","t1","t2","t3","coverage","pct","nmod","ncanon","nother","ndel","nfail","ndiff","nnocall"]
rows=[]
for f in sorted(glob.glob(f"{PILE}/1-*_pileup.bed")):
    d=pd.read_csv(f,sep="\t",header=None,names=cols,usecols=["start","mod","strand","coverage","pct"])
    d=d[d["coverage"]>=10]
    d["rep"]=Path(f).stem.split("_")[0]
    rows.append(d)
pile=pd.concat(rows,ignore_index=True)
MODMAP={"a":"6mA","21839":"4mC","m":"5mC"}
pile["mod_type"]=pile["mod"].astype(str).map(MODMAP)
pile=pile.dropna(subset=["mod_type"])
key={(int(p),s,mt):v for p,s,mt,v in zip(pile.start,pile.strand,pile.mod_type,pile.pct)}
# read-weighted mean across replicates at each (pos,strand,mod)
agg=pile.groupby(["start","strand","mod_type"]).apply(lambda g: np.average(g.pct,weights=g["coverage"]),include_groups=False)

out=[]
for name,(positions,L) in inst.items():
    for off in range(L):
        for mt in ("4mC","6mA"):
            vals=[]
            for p0,st in positions:
                p = p0+off if st=="+" else p0-off
                vals.append(agg.get((p,st,mt),0.0))
            base = MOT[name][off]
            v=np.array(vals); called=v[v>=50]
            out.append(dict(motif=name,offset=off,base=base,mod_type=mt,
                            n_instances=len(positions),mean_pct_all=round(float(v.mean()),2),
                            n_called=int(called.size),
                            mean_pct_at_called_sites=round(float(called.mean()),1) if called.size else 0.0,
                            pct_instances_called=round(100*float((v>=50).mean()),1)))
df=pd.DataFrame(out)
df=df[~((df.mod_type=="4mC")&(df.base!="C")) & ~((df.mod_type=="6mA")&(df.base!="A"))]
df.to_csv(H/"per_position_motif_frequency_T1.tsv",sep="\t",index=False)
print(df.to_string(index=False))
