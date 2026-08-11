"""C1/C2: is the GCCGGC 'core enrichment' a copy-number (oriC coverage) artefact?

The linear chromosome has oriC near the centre (~4.27 Mb), inside the 'core'
(1.5-7.17 Mb). During active growth (T1) replication raises copy number near
oriC -> higher coverage in core -> more GCCGGC positions pass the depth>=10
site-calling threshold -> inflated 'site count' in core. Per-site methylation
FREQUENCY (a ratio) is copy-number independent.

Test: among ALL genomic GCCGGC cytosine positions (denominator = sequence, not
detection), compare the methylated fraction core vs arm, STRATIFIED BY COVERAGE.
If core>arm only at low coverage / collapses when coverage-matched -> artefact.
If core>arm persists within matched coverage bands -> real enrichment.

Writes a coverage-position list then streams the T1 pileup.
"""
import os, re, subprocess
import numpy as np, pandas as pd

REF = "/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa"
PILE = "/Users/okaban/bioinfo/methyl/260102_M145/analysis/cursor_results/pileup/1-1_pileup.bed"
OUT = os.path.expanduser("~/bioinfo/rna-seq/11_epigenome_integration/analysis/78_reviewer_robustness/tables")
CORE_LO, CORE_HI = 1_500_000, 7_170_000

# --- genome -> GCCGGC methylated-C positions (+ at s+2, - at s+3) ---
seq = []
with open(REF) as f:
    for ln in f:
        if not ln.startswith(">"):
            seq.append(ln.strip().upper())
seq = "".join(seq)
sites = {}  # "pos:strand" -> region
for m in re.finditer("(?=GCCGGC)", seq):
    s = m.start()
    for pos, strand in [(s+2, "+"), (s+3, "-")]:
        region = "core" if CORE_LO <= pos <= CORE_HI else "arm"
        sites[f"{pos}:{strand}"] = region
print(f"GCCGGC genomic methylated-C positions: {len(sites)} "
      f"(core {sum(v=='core' for v in sites.values())}, arm {sum(v=='arm' for v in sites.values())})")
keyfile = "/tmp/gccggc_keys.tsv"
with open(keyfile, "w") as f:
    for k, v in sites.items():
        f.write(f"{k}\t{v}\n")

# --- stream pileup, collect (coverage, fraction, region) for GCCGGC C positions ---
# pileup cols: 1 chrom 2 start 3 end 4 code 5 score 6 strand ... 10 Nvalid 11 fraction% 12 Nmod
awk = r'''
BEGIN{while((getline line < "%s")>0){split(line,a,"\t"); reg[a[1]]=a[2]}}
$4=="21839"{key=$2":"$6; if(key in reg){print $10"\t"$11"\t"reg[key]}}
''' % keyfile
res = subprocess.run(["awk", "-F", "\t", awk, PILE], capture_output=True, text=True)
rows = [l.split("\t") for l in res.stdout.strip().split("\n") if l]
df = pd.DataFrame(rows, columns=["coverage", "frac", "region"])
df["coverage"] = df["coverage"].astype(float); df["frac"] = df["frac"].astype(float)
df = df[df["coverage"] >= 10].copy()
df["meth"] = df.frac >= 50.0  # canonical site threshold
print(f"\nGCCGGC C positions with coverage>=10: {len(df)} "
      f"(core {sum(df.region=='core')}, arm {sum(df.region=='arm')})")

def summ(d):
    return dict(n=len(d), mean_cov=d["coverage"].mean(), pct_meth=100*d.meth.mean(),
                n_meth=int(d.meth.sum()))

print("\n=== OVERALL (cov>=10): methylated fraction of GCCGGC motifs, core vs arm ===")
for reg in ["core", "arm"]:
    s = summ(df[df.region == reg])
    print(f"  {reg}: n_motif={s['n']}  mean_cov={s['mean_cov']:.1f}  "
          f"%methylated={s['pct_meth']:.1f}%  (n_meth={s['n_meth']})")

print("\n=== COVERAGE-STRATIFIED: %methylated of GCCGGC motifs, core vs arm ===")
bins = [(10, 20), (20, 30), (30, 50), (50, 1e9)]
rows_out = []
for lo, hi in bins:
    sub = df[(df["coverage"] >= lo) & (df["coverage"] < hi)]
    c = sub[sub.region == "core"]; a = sub[sub.region == "arm"]
    cc, aa = summ(c), summ(a)
    print(f"  cov[{lo},{hi if hi<1e9 else '+'}): "
          f"core n={cc['n']} %meth={cc['pct_meth']:.1f}  |  "
          f"arm n={aa['n']} %meth={aa['pct_meth']:.1f}")
    rows_out.append(dict(cov_bin=f"{lo}-{hi if hi<1e9 else 'inf'}",
                         core_n=cc['n'], core_pct_meth=cc['pct_meth'],
                         arm_n=aa['n'], arm_pct_meth=aa['pct_meth']))

# --- coverage-MATCHED site density: methylated sites per Mb at a capped coverage ---
# Recount methylated sites core vs arm using only motifs in a common coverage band
band = df[(df["coverage"] >= 10) & (df["coverage"] <= 25)]
core_kb = (CORE_HI - CORE_LO) / 1000.0
arm_kb = (8_670_000 - (CORE_HI - CORE_LO)) / 1000.0
core_meth = int(band[(band.region == "core") & band.meth].shape[0])
arm_meth = int(band[(band.region == "arm") & band.meth].shape[0])
print("\n=== coverage-matched band [10,25x]: methylated GCCGGC sites per Mb ===")
print(f"  core: {core_meth} sites / {core_kb/1000:.2f} Mb = {core_meth/(core_kb/1000):.1f}/Mb")
print(f"  arm : {arm_meth} sites / {arm_kb/1000:.2f} Mb = {arm_meth/(arm_kb/1000):.1f}/Mb")

pd.DataFrame(rows_out).to_csv(os.path.join(OUT, "C1_coverage_stratified.tsv"), sep="\t", index=False)
print("\nsaved -> C1_coverage_stratified.tsv")
print("\nINTERPRETATION:")
print(" - if %methylated(core) ~ %methylated(arm) within each coverage bin,")
print("   the 'core enrichment' of SITE COUNT is a coverage/copy-number detection effect.")
print(" - if core still > arm within matched coverage bins, the enrichment is real.")
