#!/usr/bin/env python3
"""REFA-06: GCCGGC strand symmetry from modkit pileups (4mC code 21839, valid coverage >=10).
For every GCCGGC instance the 4mC-carrying C on each strand is determined empirically (offset with
highest mean % modified). Definitions reported: (A) among instances 'methylated' (>=50% on >=1 strand):
fraction hemi (one strand >=50%) vs full (both >=50%); (B) per-strand frequency statistics at methylated
strands; (C) Pearson r between + and - strand % among methylated instances and across all instances;
(D) mean + vs - strand % across all instances.
"""
import re, sys, subprocess, io, pandas as pd, numpy as np
from scipy.stats import pearsonr
REF="/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa"; PD="/Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup"
seq="".join(l.strip() for l in open(REF) if not l.startswith(">")).upper()
starts=[m.start() for m in re.finditer("GCCGGC", seq)]; print("GCCGGC instances", len(starts), file=sys.stderr)
rows=[]; per_inst=[]
for samp in sys.argv[1:]:
    txt=subprocess.run(f"awk '$4==21839 && $10>=10' {PD}/{samp}_pileup.bed | cut -f2,6,10,11", shell=True, capture_output=True, text=True).stdout
    p=pd.read_csv(io.StringIO(txt), sep="\t", header=None, names=["pos","strand","cov","pct"])
    plus={ (r.pos):r.pct for r in p[p.strand=="+"].itertuples()}; minus={ (r.pos):r.pct for r in p[p.strand=="-"].itertuples()}
    # find 4mC offset: mean pct per offset per strand
    off={}
    for k in range(6):
        vp=[plus[s+k] for s in starts if s+k in plus]; vm=[minus[s+k] for s in starts if s+k in minus]
        off[k]=(np.mean(vp) if vp else np.nan, len(vp), np.mean(vm) if vm else np.nan, len(vm))
    kp=max((k for k in off if not np.isnan(off[k][0])), key=lambda k: off[k][0]); km=max((k for k in off if not np.isnan(off[k][2])), key=lambda k: off[k][2])
    d=pd.DataFrame([(s, plus.get(s+kp,np.nan), minus.get(s+km,np.nan)) for s in starts], columns=["start","plus_pct","minus_pct"]).dropna()
    meth=d[(d.plus_pct>=50)|(d.minus_pct>=50)]; full=meth[(meth.plus_pct>=50)&(meth.minus_pct>=50)]
    r_meth=pearsonr(meth.plus_pct, meth.minus_pct)[0] if len(meth)>2 else np.nan; r_all=pearsonr(d.plus_pct,d.minus_pct)[0]
    # per-strand frequency at methylated strands
    ms=np.concatenate([meth.plus_pct[meth.plus_pct>=50].values, meth.minus_pct[meth.minus_pct>=50].values])
    # opposite strand at hemi sites
    opp=np.concatenate([meth.minus_pct[(meth.plus_pct>=50)].values, meth.plus_pct[(meth.minus_pct>=50)].values])
    rows.append(dict(sample=samp, plus_4mC_offset=kp, minus_4mC_offset=km, n_instances_both_strands_cov10=len(d), n_methylated_ge50_any_strand=len(meth), n_full_both_ge50=len(full), pct_hemi=round(100*(1-len(full)/len(meth)),2), pct_full=round(100*len(full)/len(meth),2),
        median_pct_at_methylated_strands=np.median(ms), mean_pct_at_methylated_strands=round(ms.mean(),1), frac_methylated_strands_ge90=round((ms>=90).mean(),3), frac_methylated_strands_ge99=round((ms>=99).mean(),3),
        mean_pct_opposite_strand_at_methylated=round(opp.mean(),1), frac_opposite_strand_ge50=round((opp>=50).mean(),4), pearson_r_meth=round(r_meth,3), pearson_r_all=round(r_all,3), mean_plus_all=round(d.plus_pct.mean(),2), mean_minus_all=round(d.minus_pct.mean(),2)))
    d["sample"]=samp; per_inst.append(d)
R=pd.DataFrame(rows); R.to_csv("tables/REFA06_GCCGGC_strand_symmetry_pileup.tsv", sep="\t", index=False); pd.concat(per_inst).to_csv("tables/REFA06_GCCGGC_per_instance_strand_pct.tsv.gz", sep="\t", index=False)
pd.set_option("display.width",300); print(R.T.to_string())
