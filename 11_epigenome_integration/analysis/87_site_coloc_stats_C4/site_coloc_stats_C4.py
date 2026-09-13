#!/usr/bin/env python3
"""87 — site-level AAGCCCG 4mC/6mA co-localisation statistics with CORRECT motif offsets
(ledger B01b, 2026-09-13). Replaces the numbers of 68_or_permutation/A3 (OR 138,440;
244 in / 1,090 not; universe 1,238,215) which were derived through the 1-bp-misaligned
`sequence` window (see 79_comod_full_denominator/README.md).

Design (mirrors 68_: 2x2 of candidate units inside vs outside AAGCCCG x co-modified or not;
Fisher exact one-sided 'greater'; Haldane-corrected OR; log-OR normal 95% CI; hypergeometric
label-permutation null, 10,000 draws, seed 42).

Unit = a same-strand candidate (A, C) position pair in the AAGCCCG dual-mark geometry:
  an A located 3 bp or 4 bp upstream (5') of a C on the same strand (= A1-C4 and A0-C4 spacing).
  Enumerated genome-wide from the reference on both strands (68_'s universe of 1,238,215
  candidate pairs was hard-coded and is not reproducible; this universe is explicit).
Co-modified = the C is a canonical 4mC site AND the A a canonical 6mA site (depth >= 10,
  freq >= 50%; 07_motif_analysis/methylation_site_sequences.csv, 0-based positions), on the
  same strand. 'Pooled' = sites deduplicated over T1-T3 (as 65_/68_); 'T1' = T1 sites only.
Inside AAGCCCG = the pair's C is offset 4 of an AAGCCCG instance (so the A is offset 0 or 1).
"""
import re, numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
REF = "/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa"
SITES = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/07_motif_analysis/methylation_site_sequences.csv"
H = Path(__file__).resolve().parent
seq = "".join(l.strip() for l in open(REF) if not l.startswith(">")).upper(); G = len(seq)
# canonical sites
df = pd.read_csv(SITES)
def siteset(sub, mod): return set(zip(sub.position[sub.mod_type==mod].astype(int), sub.strand[sub.mod_type==mod]))
pooled = df.drop_duplicates(["position","strand","mod_type"])
S = {"pooled": (siteset(pooled,"4mC"), siteset(pooled,"6mA")), "T1": (siteset(df[df.timepoint=="T1"],"4mC"), siteset(df[df.timepoint=="T1"],"6mA"))}
# AAGCCCG C4 positions (0-based) per strand
c4 = {("+", m.start()+4) for m in re.finditer("AAGCCCG", seq)} | {("-", m.start()+2) for m in re.finditer("CGGGCTT", seq)}
n_inst = len(c4); assert n_inst == 1334, n_inst
# enumerate candidate pairs: + strand: A at c-3 / c-4 ; - strand: A (=T on +) at c+3 / c+4
cands = []   # (strand, c_pos, a_pos, spacing, in_motif)
plus_C = np.frombuffer(seq.encode(), dtype="S1") == b"C"; plus_A = np.frombuffer(seq.encode(), dtype="S1") == b"A"
plus_G = np.frombuffer(seq.encode(), dtype="S1") == b"G"; plus_T = np.frombuffer(seq.encode(), dtype="S1") == b"T"
rows = []
for sp in (3, 4):
    idx = np.where(plus_C[sp:] & plus_A[:-sp])[0] + sp          # + strand: C at idx, A at idx-sp
    rows.append(pd.DataFrame(dict(strand="+", c_pos=idx, a_pos=idx-sp, spacing=sp)))
    idx = np.where(plus_G[:-sp] & plus_T[sp:])[0]                # - strand: C(-)=G(+) at idx, A(-)=T(+) at idx+sp
    rows.append(pd.DataFrame(dict(strand="-", c_pos=idx, a_pos=idx+sp, spacing=sp)))
cand = pd.concat(rows, ignore_index=True)
cand["in_motif"] = [ (s, c) in c4 for s, c in zip(cand.strand, cand.c_pos) ]
assert cand.in_motif.sum() == 2*n_inst, cand.in_motif.sum()
out_rows = []
for scope, (m4, m6) in S.items():
    cand["comod"] = [ ((c, s) in m4) and ((a, s) in m6) for c, a, s in zip(cand.c_pos, cand.a_pos, cand.strand) ]
    a = int((cand.in_motif & cand.comod).sum()); b = int((cand.in_motif & ~cand.comod).sum())
    c = int((~cand.in_motif & cand.comod).sum()); d = int((~cand.in_motif & ~cand.comod).sum())
    by_sp = cand[cand.in_motif & cand.comod].spacing.value_counts().to_dict()
    OR_raw, p_fisher = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
    aa, bb, cc, dd = (np.array([a, b, c, d]) + (0.5 if min(a, b, c, d) == 0 else 0)).astype(float)
    OR = aa*dd/(bb*cc); se = np.sqrt(1/aa+1/bb+1/cc+1/dd); lo, hi = np.exp(np.log(OR)-1.96*se), np.exp(np.log(OR)+1.96*se)
    n_tot, n_in, n_comod = a+b+c+d, a+b, a+c
    exp_in = n_in * n_comod / n_tot; fold = a / exp_in
    rng = np.random.default_rng(42); N = 10_000
    a_perm = rng.hypergeometric(n_comod, n_tot-n_comod, n_in, N)
    p_perm = (np.sum(a_perm >= a) + 1) / (N + 1)
    rate_in = a / n_in; rate_out = c / (c+d); rr = rate_in / rate_out if c else np.inf
    n_inst_comod = cand[cand.in_motif & cand.comod].drop_duplicates(["strand","c_pos"]).shape[0]   # distinct AAGCCCG instances (one C4 per instance) with >=1 co-mod pair
    out_rows.append(dict(scope=scope, n_instances_with_comod_pair=n_inst_comod, in_comod=a, in_notcomod=b, out_comod=c, out_notcomod=d, n_candidate_pairs=n_tot, n_in_motif_pairs=n_in,
                         in_comod_spacing4_A0C4=by_sp.get(4,0), in_comod_spacing3_A1C4=by_sp.get(3,0),
                         OR=OR, OR_CI95_lower=lo, OR_CI95_upper=hi, haldane_corrected=(min(a,b,c,d)==0), fisher_p_greater=p_fisher,
                         expected_in_comod_H0=exp_in, fold_over_expectation=fold, perm_p=f"< {1/(N+1):.1e}" if np.sum(a_perm>=a)==0 else f"{p_perm:.4f}",
                         perm_null_max_in_comod=int(a_perm.max()), rate_in=rate_in, rate_out=rate_out, rate_ratio=rr))
res = pd.DataFrame(out_rows); res.to_csv(H/"tables/site_coloc_2x2_C4.tsv", sep="\t", index=False, float_format="%.6g")
# cross-check instance counts against the independent census (79_/site_census_AAGCCCG_instances_summary.tsv)
chk = pd.read_csv(H.parent/"79_comod_full_denominator/site_census_AAGCCCG_instances_summary.tsv", sep="\t", index_col=0)["n_instances_with_pair"]
assert int(res.set_index("scope").loc["T1","n_instances_with_comod_pair"]) == int(chk["T1"]), (res.n_instances_with_comod_pair.tolist(), chk.to_dict())
assert int(res.set_index("scope").loc["pooled","n_instances_with_comod_pair"]) == int(chk["pooled_any_timepoint"])
print("instance-count cross-check vs 79_ census: OK")
print(res.T.to_string())
# 6mA occupancy at AAGCCCG, T1 (corrected denominators)
a_pos = {("+", m.start()+k) for m in re.finditer("AAGCCCG", seq) for k in (0,1)} | {("-", m.start()+6-k) for m in re.finditer("CGGGCTT", seq) for k in (0,1)}
m6T1 = S["T1"][1]; n6 = sum((p, s) in m6T1 for s, p in a_pos)
inst_of = {}
for m in re.finditer("AAGCCCG", seq): inst_of[("+", m.start())] = m.start(); inst_of[("+", m.start()+1)] = m.start()
for m in re.finditer("CGGGCTT", seq): inst_of[("-", m.start()+6)] = m.start(); inst_of[("-", m.start()+5)] = m.start()
inst6 = {inst_of[(s, p)] for s, p in a_pos if (p, s) in m6T1}
m4T1 = S["T1"][0]; n4 = sum((c, s) in m4T1 for s, c in c4)
occ = pd.DataFrame([dict(metric="6mA canonical sites at AAGCCCG A0/A1, T1", n=n6, denominator="A0+A1 positions = 2 x 1,334 = 2,668", pct=round(n6/2668*100,1)),
                    dict(metric="AAGCCCG instances with >=1 canonical 6mA (A0 or A1), T1", n=len(inst6), denominator="1,334 instances", pct=round(len(inst6)/1334*100,1)),
                    dict(metric="4mC canonical sites at AAGCCCG C4, T1", n=n4, denominator="1,334 C4 positions", pct=round(n4/1334*100,1))])
occ.to_csv(H/"tables/AAGCCCG_occupancy_T1_C4.tsv", sep="\t", index=False); print(occ.to_string(index=False))
