#!/usr/bin/env python3
"""88 — AAGCCCG 6mA occupancy series (T1/T2/T3) and Clark-Evans on the corrected census.

Two manuscript numbers still rest on the 1-bp-misaligned `sequence` window (see
79_/87_ READMEs):
  (i)  EN L42/L201 series "31.6% -> 7.2% -> 4.3%" of AAGCCCG 6mA occupancy. Its
       denominator is untraceable (447/0.316 = 1,415, which is neither the 2,668
       A0/A1 positions nor the 1,334 instances). Recomputed here on both explicit
       bases, all three timepoints, same canonical-site definition as 87_.
  (ii) EN L170 Clark-Evans "244 dual-modification events, R = 1.761, n_eff = 88.8".
       Recomputed on the corrected same-instance census (406 pooled / 214 T1).
       1D circular chromosome: under CSR the mean nearest-neighbour distance is
       L/(2n); R = observed/expected; z from the normal approximation with
       sd = L/(2n*sqrt(n)) (Clark & Evans 1954, 1D circular form).
       n_eff = 88.8 is not reproducible from any script in the repo and is not
       reproduced here.
"""
import re, numpy as np, pandas as pd
from pathlib import Path
REF="/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa"
SITES="/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/07_motif_analysis/methylation_site_sequences.csv"
H=Path(__file__).resolve().parent
seq="".join(l.strip() for l in open(REF) if not l.startswith(">")).upper(); L=len(seq)
df=pd.read_csv(SITES)
inst=[("+",m.start()) for m in re.finditer("AAGCCCG",seq)]+[("-",m.start()+6) for m in re.finditer("CGGGCTT",seq)]
assert len(inst)==1334
def A_positions(s,p): return [(p,s),(p-1,s)] if s=="-" else [(p,s),(p+1,s)]
def C4(s,p): return (p-4,s) if s=="-" else (p+4,s)
Apos={x for s,p in inst for x in A_positions(s,p)}; C4pos={C4(s,p) for s,p in inst}
assert len(Apos)==2668 and len(C4pos)==1334
rows=[]
for tp in ("T1","T2","T3"):
    sub=df[df.timepoint==tp]
    m6={(int(p),s) for p,s in zip(sub.position[sub.mod_type=="6mA"],sub.strand[sub.mod_type=="6mA"])}
    m4={(int(p),s) for p,s in zip(sub.position[sub.mod_type=="4mC"],sub.strand[sub.mod_type=="4mC"])}
    n_pos=len(Apos&m6); n_inst=sum(any(a in m6 for a in A_positions(s,p)) for s,p in inst)
    n_c4=len(C4pos&m4)
    rows.append(dict(timepoint=tp,n_6mA_positions=n_pos,pct_of_2668=round(100*n_pos/2668,1),
                     n_instances_with_6mA=n_inst,pct_of_1334=round(100*n_inst/1334,1),
                     n_4mC_at_C4=n_c4,pct_C4_of_1334=round(100*n_c4/1334,1)))
occ=pd.DataFrame(rows); occ.to_csv(H/"AAGCCCG_occupancy_series.tsv",sep="\t",index=False)

def clark_evans(pos,L):
    p=np.sort(np.asarray(pos)); n=len(p)
    if n<2: return dict(n=n,R=np.nan,z=np.nan)
    d=np.diff(np.concatenate([p,[p[0]+L]]))          # circular gaps
    nn=np.minimum(d,np.roll(d,1))                     # nearest neighbour = min(left,right) gap
    obs=nn.mean(); exp=L/(2*n); sd=L/(2*n*np.sqrt(n))
    z=(obs-exp)/sd
    from scipy.stats import norm
    return dict(n=n,mean_nn_bp=round(obs,1),expected_nn_bp=round(exp,1),R=round(obs/exp,3),
                z=round(z,2),p_two_sided=float(f"{2*norm.sf(abs(z)):.3g}"),
                pattern="dispersed" if obs>exp else "clustered")
cen=pd.read_csv(H.parent/"79_comod_full_denominator/site_census_AAGCCCG_pairs.tsv",sep="\t")
cen=cen[cen.rule=="same_instance"] if "rule" in cen.columns else cen
ce=[]
for scope,sel in [("pooled_any_timepoint",cen[cen.scope=="pooled_any_timepoint"] if "scope" in cen.columns else cen),
                  ("T1",cen[cen.scope=="T1"] if "scope" in cen.columns else cen.iloc[:0])]:
    if len(sel)==0: continue
    ce.append(dict(scope=scope,**clark_evans(sel["pos_4mC"].astype(int).values,L)))
pd.DataFrame(ce).to_csv(H/"clark_evans_C4_census.tsv",sep="\t",index=False)
print(occ.to_string(index=False)); print(); print(pd.DataFrame(ce).to_string(index=False))

# --- motif-conditioned null: the right comparison ---------------------------
# CSR treats every genomic position as an equally likely site, but a dual mark can
# only occur where an AAGCCCG instance is. The informative question is whether the
# co-modified instances are more clustered than a random draw of the same number of
# instances from the 1,334 that exist. 10,000 draws, seed 42.
rng=np.random.default_rng(42)
inst_c4=np.sort(np.array([p for p,_ in C4pos]))
def mean_nn(p,L):
    p=np.sort(p); d=np.diff(np.concatenate([p,[p[0]+L]])); return np.minimum(d,np.roll(d,1)).mean()
rows2=[]
for scope,sel in [("pooled_any_timepoint",cen[cen.scope=="pooled_any_timepoint"]),("T1",cen[cen.scope=="T1"])]:
    obs=mean_nn(sel["pos_4mC"].astype(int).values,L); n=len(sel)
    null=np.array([mean_nn(rng.choice(inst_c4,n,replace=False),L) for _ in range(10_000)])
    rows2.append(dict(scope=scope,n=n,obs_mean_nn_bp=round(obs,1),null_mean_nn_bp=round(null.mean(),1),
        ratio_obs_null=round(obs/null.mean(),3),p_two_sided=round(2*min((null<=obs).mean(),(null>=obs).mean()),4),
        pattern="more clustered than motif set" if obs<null.mean() else "more dispersed than motif set"))
mc=pd.DataFrame(rows2); mc.to_csv(H/"clark_evans_motif_conditioned_null.tsv",sep="\t",index=False)
print(); print(mc.to_string(index=False))
