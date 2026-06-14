"""Genome-wide cross-method overlap: Pisciotta 2023 (bisulfite, 321 upstream-methylated genes)
vs our single-molecule GCCGGC m4C proximity, and vs our Exposed/Shielded regulator classes.

Tests whether Pisciotta's promoter-methylated genes are independently closer to GCCGGC sites
in our Nanopore data than the genome background (cross-method validation), and how the
regulators among them split into our Exposed/Shielded classes.
"""
import numpy as np, pandas as pd
from pathlib import Path
from scipy.stats import mannwhitneyu, fisher_exact

BASE = Path(__file__).resolve().parents[4]
GM = BASE/'05_annotation/analysis/05_annotation_260128_v1/tables/gene_master_with_BGC.tsv'
SITES = BASE/'11_epigenome_integration/analysis/37_defense_island_GCCGGC/tables/GCCGGC_sites_by_timepoint.tsv'
REG = BASE/'11_epigenome_integration/analysis/52_shielded_exposed_boundary/tables/all_genes_features_unified_n57.tsv'
PISC = open('/tmp/pisc_sco.txt').read().split()

pisc = set(PISC)
gm = pd.read_csv(GM, sep='\t')
gm = gm[gm['old_locus_tag'].astype(str).str.match('SCO')].copy()
t1 = pd.read_csv(SITES, sep='\t'); t1 = t1[t1['timepoint']=='T1']
spos = np.sort(t1['position'].values)
print(f"genes(SCO)={len(gm)}  Pisciotta SCO={len(pisc)}  GCCGGC-T1 sites={len(spos)}")

# TSS proxy = start (+ strand) / end (- strand); nearest T1 GCCGGC distance
def tss(r): return r['start'] if r['strand']=='+' else r['end']
gm['tss'] = gm.apply(tss, axis=1)
idx = np.searchsorted(spos, gm['tss'].values)
idx = np.clip(idx, 1, len(spos)-1)
left = spos[idx-1]; right = spos[idx]
gm['nearest_gccggc'] = np.minimum(np.abs(gm['tss'].values-left), np.abs(gm['tss'].values-right))
gm['in_pisc'] = gm['old_locus_tag'].isin(pisc)

a = gm[gm['in_pisc']]['nearest_gccggc']; b = gm[~gm['in_pisc']]['nearest_gccggc']
U,p = mannwhitneyu(a,b,alternative='less')
print(f"\n=== Cross-method: Pisciotta-321 vs background, nearest GCCGGC-T1 distance ===")
print(f"  Pisciotta median={a.median():.0f} bp (n={len(a)})  vs  background median={b.median():.0f} bp (n={len(b)})")
print(f"  Mann-Whitney (Pisciotta closer): p={p:.3e}")
for w in (293,500,1000):
    fa=(a<=w).mean(); fb=(b<=w).mean()
    od,pf=fisher_exact([[(a<=w).sum(),(a>w).sum()],[(b<=w).sum(),(b>w).sum()]])
    print(f"  within {w}bp: Pisciotta {fa:.0%} vs background {fb:.0%}  (OR={od:.2f}, Fisher p={pf:.2e})")

# overlap with our regulator classes
reg = pd.read_csv(REG, sep='\t')
reg_sco = reg.dropna(subset=['old_locus_tag'])
reg_sco = reg_sco[reg_sco['old_locus_tag'].astype(str).str.match('SCO')]
inreg = reg_sco[reg_sco['old_locus_tag'].isin(pisc)]
print(f"\n=== Pisciotta ∩ our 1,055 regulators ===")
print(f"  {len(inreg)} of {len(pisc)} Pisciotta genes are in our regulator set")
print(f"  Exposed: {(inreg['is_exposed']==1).sum()}  /  Shielded: {(inreg['is_exposed']==0).sum()}")
print(f"  median nearest_methyl_distance (regulators, Pisciotta): {inreg['nearest_methyl_distance'].median():.0f} bp")
print(f"  Exposed prevalence: Pisciotta-regulators {(inreg['is_exposed']==1).mean():.1%} vs all regulators {(reg['is_exposed']==1).mean():.1%}")
