#!/usr/bin/env python3
"""Site-level AAGCCCG 4mC/6mA co-localisation census (current), classified by
motif-offset pair (A0/A1 x C3/C4/C5), 2026-09-13.

Input: 07_motif_analysis/methylation_site_sequences.csv (canonical sites,
depth>=10, freq>=50%, per timepoint). COORDINATE NOTE: the `position` column is
0-based (modkit bedMethyl start): ref[position] is the modified base for all
four mod/strand classes; the 31-bp `sequence` column is centred one base
upstream (it matches a 1-based reading) and therefore mis-assigns motif
offsets by 1. This script maps sites to AAGCCCG offsets from the reference
directly (0-based) and does NOT use the `sequence` column.

Pairing: 4mC site and 6mA site on the same strand within the SAME AAGCCCG
instance (also reported: 65_per_read_comod rule = same strand, |d|<=10 bp,
which admits a few cross-instance pairs). Pooled = sites deduplicated by
(position, strand, mod_type) over T1-T3, as in 65_per_read_comod.
"""
import re, pandas as pd
from pathlib import Path
REF = "/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa"
SITES = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/07_motif_analysis/methylation_site_sequences.csv"
H = Path(__file__).resolve().parent
seq = "".join(l.strip() for l in open(REF) if not l.startswith(">")).upper()
inst = {}   # (0-based pos, strand) -> (instance_id, offset)
for m in re.finditer("AAGCCCG", seq):
    for k in range(7): inst[(m.start()+k, "+")] = (f"+{m.start()}", k)
for m in re.finditer("CGGGCTT", seq):
    for k in range(7): inst[(m.start()+6-k, "-")] = (f"-{m.start()}", k)
df = pd.read_csv(SITES)
df["inst"] = [inst.get((int(p), s), (None, None))[0] for p, s in zip(df.position, df.strand)]
df["offset"] = [inst.get((int(p), s), (None, None))[1] for p, s in zip(df.position, df.strand)]
aag = df[df.inst.notna()].copy()
print("offset distribution of canonical sites inside AAGCCCG (unique sites):")
print(aag.drop_duplicates(["position","strand","mod_type"]).groupby(["mod_type","strand","offset"]).size().to_string())
def pairs(sub, rule):
    d4 = sub[sub.mod_type=="4mC"]; d6 = sub[sub.mod_type=="6mA"]; rows = []
    for p4, s4, i4, o4 in d4[["position","strand","inst","offset"]].to_numpy():
        for p6, s6, i6, o6 in d6[["position","strand","inst","offset"]].to_numpy():
            if s4 != s6: continue
            if rule == "same_instance" and i4 != i6: continue
            if rule == "within10bp" and abs(int(p4)-int(p6)) > 10: continue
            rows.append(dict(pos_4mC=int(p4), pos_6mA=int(p6), strand=s4, distance=abs(int(p4)-int(p6)),
                             pair_class=f"A{int(o6)}-C{int(o4)}" if i4 == i6 else "cross-instance", inst=i4))
    return pd.DataFrame(rows)
out = []
for rule in ("same_instance", "within10bp"):
    for tp in ("T1","T2","T3"):
        p = pairs(aag[aag.timepoint==tp], rule); p.insert(0,"scope",tp); p.insert(0,"rule",rule); out.append(p)
    u = aag.drop_duplicates(["position","strand","mod_type"])
    p = pairs(u, rule); p.insert(0,"scope","pooled_any_timepoint"); p.insert(0,"rule",rule); out.append(p)
allp = pd.concat(out, ignore_index=True)
allp.to_csv(H/"site_census_AAGCCCG_pairs.tsv", sep="\t", index=False)
summ = allp.groupby(["rule","scope","pair_class"]).size().unstack(fill_value=0); summ["total"] = summ.sum(axis=1)
summ.to_csv(H/"site_census_AAGCCCG_pairs_summary.tsv", sep="\t")
print(summ.to_string())
# instance-level: AAGCCCG instances (of 1,334) with >=1 same-instance pair; and site counts per timepoint
si = allp[allp.rule=="same_instance"].groupby("scope")["inst"].nunique().rename("n_instances_with_pair")
sites = aag.groupby(["timepoint","mod_type"]).size().unstack(fill_value=0)
sites.loc["pooled_unique"] = aag.drop_duplicates(["position","strand","mod_type"]).groupby("mod_type").size()
inst_tot = pd.DataFrame({"n_instances_with_pair": si}).join(sites, how="outer")
inst_tot.to_csv(H/"site_census_AAGCCCG_instances_summary.tsv", sep="\t")
print(inst_tot.to_string()); print("AAGCCCG instances total:", len({v[0] for v in inst.values()}))
