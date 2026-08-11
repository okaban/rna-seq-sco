"""C6: Comprehensive Exposed(62) vs Shielded(989) metric scan.
Tests which gene-level metrics actually differ, to substantiate (or refute) the
permissive thesis. Canonical 62/989 labels from the classification table.
Seed 42. 2026-06-29."""
import csv, math, statistics as st
from collections import Counter
import numpy as np
from scipy import stats

np.random.seed(42)
SRC = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/52_shielded_exposed_boundary/tables/SuppTable_regulatory_gene_classification_n1051.tsv"
OUT = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/82_exposed_shielded_metric_scan/tables/exposed_vs_shielded_metric_scan.tsv"

rows = list(csv.DictReader(open(SRC), delimiter="\t"))
exp = [r for r in rows if r["class"] == "Exposed"]
shi = [r for r in rows if r["class"] == "Shielded"]
CENTER = (min(float(r["tss"]) for r in rows) + max(float(r["tss"]) for r in rows)) / 2  # chromosome-center proxy for oriC

def fnum(r, k):
    try: return float(r[k])
    except: return None

def cliffs_delta(a, b):
    a, b = np.asarray(a), np.asarray(b)
    gt = sum((a[:,None] > b[None,:]).sum(axis=1))
    lt = sum((a[:,None] < b[None,:]).sum(axis=1))
    return (gt - lt) / (len(a)*len(b))

def mwu(metric, getter):
    a = [getter(r) for r in exp if getter(r) is not None]
    b = [getter(r) for r in shi if getter(r) is not None]
    U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return dict(metric=metric, n_exp=len(a), n_shi=len(b), test="MWU",
               effect=f"Cliff_d={cliffs_delta(a,b):+.3f}",
               med_exp=round(np.median(a),3), med_shi=round(np.median(b),3), p=p)

def fisher(metric, pred):
    a11 = sum(1 for r in exp if pred(r)); a12 = len(exp)-a11
    a21 = sum(1 for r in shi if pred(r)); a22 = len(shi)-a21
    OR, p = stats.fisher_exact([[a11,a12],[a21,a22]])
    return dict(metric=metric, n_exp=len(exp), n_shi=len(shi), test="Fisher",
               effect=f"OR={OR:.2f}", med_exp=f"{a11}/{len(exp)}", med_shi=f"{a21}/{len(shi)}", p=p)

res = []
# expression metrics (the key permissive test)
res.append(mwu("baseline_expr_log10(baseMean)", lambda r: math.log10(fnum(r,"baseMean")+1) if fnum(r,"baseMean") is not None else None))
res.append(mwu("signed_LFC_T2vsT1", lambda r: fnum(r,"LFC_T2vsT1")))
res.append(mwu("signed_LFC_T3vsT1", lambda r: fnum(r,"LFC_T3vsT1")))
res.append(mwu("abs_LFC_T2vsT1(variability)", lambda r: abs(fnum(r,"LFC_T2vsT1")) if fnum(r,"LFC_T2vsT1") is not None else None))
res.append(mwu("abs_LFC_T3vsT1(variability)", lambda r: abs(fnum(r,"LFC_T3vsT1")) if fnum(r,"LFC_T3vsT1") is not None else None))
res.append(fisher("DEG_frac_T2(|LFC|>=1)", lambda r: fnum(r,"LFC_T2vsT1") is not None and abs(fnum(r,"LFC_T2vsT1"))>=1))
res.append(fisher("DEG_frac_T3(|LFC|>=1)", lambda r: fnum(r,"LFC_T3vsT1") is not None and abs(fnum(r,"LFC_T3vsT1"))>=1))
# position metrics (expected to differ)
res.append(fisher("position_core", lambda r: r["region_label"].lower().startswith("core")))
res.append(mwu("dist_to_chrom_center(bp)", lambda r: abs(fnum(r,"tss")-CENTER) if fnum(r,"tss") is not None else None))
# TF family enrichment (expected: MerR/LysR/LacI)
fams = [f for f,_ in Counter(r["tf_family"] for r in rows).most_common() if f]
for fam in ["MerR","LysR","LacI","TetR","Sigma","sensor_kinase"]:
    res.append(fisher(f"family_{fam}", lambda r,fam=fam: fam.lower() in (r["tf_family"] or "").lower()))

# BH-FDR
ps = [r["p"] for r in res]
order = np.argsort(ps); m = len(ps); bh = [None]*m
prev = 1.0
for rank, idx in enumerate(reversed(order)):
    k = m - rank
    val = min(prev, ps[idx]*m/k); bh[idx] = val; prev = val

with open(OUT,"w") as f:
    f.write("metric\ttest\tn_exp\tn_shi\teffect\tmed/frac_exp\tmed/frac_shi\tp_raw\tp_BH\tsig_FDR0.05\n")
    for r,q in zip(res,bh):
        f.write(f"{r['metric']}\t{r['test']}\t{r['n_exp']}\t{r['n_shi']}\t{r['effect']}\t{r['med_exp']}\t{r['med_shi']}\t{r['p']:.2e}\t{q:.2e}\t{'YES' if q<0.05 else 'no'}\n")

print(f"wrote {OUT}\n")
print(f"{'metric':38s} {'test':7s} {'effect':16s} {'p_raw':>9s} {'p_BH':>9s}  sig")
for r,q in zip(res,bh):
    print(f"{r['metric']:38s} {r['test']:7s} {r['effect']:16s} {r['p']:9.2e} {q:9.2e}  {'YES' if q<0.05 else 'no'}")
