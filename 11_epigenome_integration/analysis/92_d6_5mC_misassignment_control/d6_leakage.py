#!/usr/bin/env python3
"""D6 — in-silico 5mC->4mC misassignment control.

Ground truth: ONT modbase-validation_2024.10 synthetic all-5-mer constructs; every listed
site in all_5mers_5mC_sites.bed is a true 5mC (5mC_rep1), every C in control_rep1 is canonical.
Both were basecalled here with the manuscript's exact model chain
(dna_r10.4.1_e8.2_400bps_sup@v5.2.0 + 4mC_5mC@v1, dorado 1.1.1 bundled with MinKNOW 7.11.2).
Question: at a TRUE 5mC, how often does the model emit a 4mC call (code 21839) at or above the
thresholds the manuscript uses (0.5 .. 0.9)? Report overall and per sequence context.
"""
import pysam, sys, collections, csv
CODE_4mC=21839; CODE_5mC="m"
ref={}; name=None
for l in open("all_5mers.fa"):
    l=l.strip()
    if l.startswith(">"): name=l[1:].split()[0]; ref[name]=""
    else: ref[name]+=l
truth=collections.defaultdict(set)
for l in open("all_5mers_5mC_sites.bed"):
    c,s,e,*_=l.split("\t"); truth[c].add(int(s))
def scan(bam, want_truth):
    rows=[]  # (ctx5, p4, p5)
    for r in pysam.AlignmentFile(bam):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality<20: continue
        mb=r.modified_bases_forward if r.is_reverse else r.modified_bases  # keys: (base,strand,code)
        if not mb: continue
        # pysam gives query positions with probabilities; map query->ref
        q2r=dict((q,rp) for q,rp in r.get_aligned_pairs(matches_only=True))
        p4={q:p for q,p in mb.get(("C",0,CODE_4mC),[])}
        p5={q:p for q,p in mb.get(("C",0,CODE_5mC),[])}
        for q in set(p4)|set(p5):
            rp=q2r.get(q)
            if rp is None: continue
            # ref base must be C on the read's strand
            rb=ref[r.reference_name][rp]
            if r.is_reverse:
                if rb!="G": continue
            else:
                if rb!="C": continue
            is_truth = (rp in truth[r.reference_name]) if not r.is_reverse else False  # bed is + strand sites
            if want_truth and not is_truth: continue
            if not want_truth and r.is_reverse: continue
            s=ref[r.reference_name]; ctx=s[max(0,rp-2):rp+3]
            rows.append((ctx, p4.get(q,0)/255, p5.get(q,0)/255))
    return rows
out=[]
for label,bam,want in (("true_5mC","5mC_rep1.4mC_5mC.bam",True),("canonical_C","control_rep1.4mC_5mC.bam",False)):
    rows=scan(bam,want); n=len(rows)
    for t in (0.5,0.6,0.7,0.8,0.9):
        n4=sum(1 for _,a,b in rows if a>=t); n5=sum(1 for _,a,b in rows if b>=t)
        out.append(dict(sample=label,threshold=t,n_obs=n,called_4mC=n4,frac_4mC=round(n4/n,5),called_5mC=n5,frac_5mC=round(n5/n,5)))
    # per-context at 0.5, plus the CCGG stratum
    ctx=collections.defaultdict(lambda:[0,0,0])
    for c,a,b in rows: ctx[c][0]+=1; ctx[c][1]+= a>=0.5; ctx[c][2]+= b>=0.5
    with open(f"d6_{label}_per_context.tsv","w",newline="") as fh:
        w=csv.writer(fh,delimiter="\t"); w.writerow(["context5","n","frac_4mC_t0.5","frac_5mC_t0.5"])
        for c,(n_,a,b) in sorted(ctx.items(),key=lambda x:-x[1][1]/max(1,x[1][0])): w.writerow([c,n_,round(a/n_,4),round(b/n_,4)])
    cc=[(a,b) for c,a,b in rows if "CCGG" in c]
    if cc: out.append(dict(sample=label+"_CCGG_context",threshold=0.5,n_obs=len(cc),called_4mC=sum(a>=0.5 for a,b in cc),frac_4mC=round(sum(a>=0.5 for a,b in cc)/len(cc),5),called_5mC=sum(b>=0.5 for a,b in cc),frac_5mC=round(sum(b>=0.5 for a,b in cc)/len(cc),5)))
with open("d6_leakage_summary.tsv","w",newline="") as fh:
    w=csv.DictWriter(fh,delimiter="\t",fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
for o in out: print(o)
