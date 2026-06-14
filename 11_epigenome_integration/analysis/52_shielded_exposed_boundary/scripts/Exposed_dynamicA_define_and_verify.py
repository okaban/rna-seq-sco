"""Exposed redefinition (LOCKED 2026-06-15, definition A) + full multi-axis verification.

Exposed = regulatory gene with TSS-nearest GCCGGC m4C <= 293 bp at T1 (vegetative).
Defined on the METHYLATION MAP ONLY (no expression) -> non-circular.
Verifies: synchronized demethylation, two-system division of labour, permissiveness,
and the weak geography-controlled modulatory bias. See PAPER-DIRECTION-LOCKED-260615.
"""
import numpy as np, pandas as pd
from scipy.stats import spearmanr, pearsonr, mannwhitneyu, fisher_exact
from pathlib import Path
B = Path(__file__).resolve().parents[4]
reg = pd.read_csv(B/'11_epigenome_integration/analysis/52_shielded_exposed_boundary/tables/all_genes_features_unified_n57.tsv', sep='\t')
g   = pd.read_csv(B/'11_epigenome_integration/analysis/37_defense_island_GCCGGC/tables/GCCGGC_sites_by_timepoint.tsv', sep='\t')
aag = pd.read_csv(B/'11_epigenome_integration/analysis/36_AAGCCCG_distribution/tables/AAGCCCG_site_gene_mapping.tsv', sep='\t')
W = 293
r = reg.dropna(subset=['tss']).copy(); r['tss'] = r['tss'].astype(int)
def nd(tss, pos):
    if len(pos) == 0: return np.nan
    i = np.clip(np.searchsorted(pos, tss), 1, len(pos)-1); return min(abs(tss-pos[i-1]), abs(tss-pos[i]))
for tp in ['T1','T2','T3']:
    pos = np.sort(g[g.timepoint==tp].position.values); r['d'+tp] = r['tss'].apply(lambda t: nd(t, pos))
r['Exposed'] = r['dT1'] <= W
E, S = r[r.Exposed], r[~r.Exposed]
print(f"Exposed={len(E)} / Shielded={len(S)} (TSS'd {len(r)}; total {len(reg)})")
print(f"synchronized demethylation: near@T1={len(E)}, @T2={(E.dT2<=W).sum()}, @T3={(E.dT3<=W).sum()}")
# location + family enrichment
print("core fraction:", round((E.region==0).mean(),2), "vs Shielded", round((S.region==0).mean(),2))
for fam in ['MerR','LysR','LacI','TetR','Sensor kinase','Sigma factor']:
    a=((r.Exposed)&(r.tf_family==fam)).sum(); c=((~r.Exposed)&(r.tf_family==fam)).sum()
    od,p=fisher_exact([[a,len(E)-a],[c,len(S)-c]]); print(f"  {fam}: {a}/{len(E)} OR={od:.2f} p={p:.4f}")
# permissiveness
u,p=mannwhitneyu(E['LFC_T2vsT1'].abs().dropna(), S['LFC_T2vsT1'].abs().dropna())
print(f"permissive: |LFC_T2| Exposed {E['LFC_T2vsT1'].abs().median():.2f} vs Shielded {S['LFC_T2vsT1'].abs().median():.2f}, MWU p={p:.2f}")
# two-system division of labour
aagprom = set(aag[(aag.timepoint=='T1')&(aag['distance']<=W)]['locus_tag'])
r['AAG'] = r['locus_tag'].isin(aagprom)
print(f"two-system: GCCGGC-Exposed={r.Exposed.sum()}, AAGCCCG-promoter={r.AAG.sum()}, overlap={(r.Exposed&r.AAG).sum()}")
# weak modulatory bias (geography-controlled)
gT1=g[g.timepoint=='T1']; P=np.sort(gT1.position.values); Fq=gT1.sort_values('position').frequency.values
def occ(tss,w):
    m=(P>=tss-w)&(P<=tss+w); return Fq[m].sum() if m.any() else 0.0
r['occ2k']=r['tss'].apply(lambda t:occ(t,2000))
d=r[['occ2k','LFC_T2vsT1','region']].dropna()
rx=d.occ2k.rank()-np.polyval(np.polyfit(d.region.rank(),d.occ2k.rank(),1),d.region.rank())
ry=d.LFC_T2vsT1.rank()-np.polyval(np.polyfit(d.region.rank(),d.LFC_T2vsT1.rank(),1),d.region.rank())
rr,pp=pearsonr(rx,ry); print(f"weak modulatory bias (occ2k vs LFC_T2, region-controlled): r={rr:.3f} p={pp:.4f}")
E[['locus_tag','old_locus_tag','gene_name','product','tf_family','region','tss','dT1','dT2','dT3','LFC_T2vsT1','LFC_T3vsT1']].to_csv(
    B/'11_epigenome_integration/analysis/52_shielded_exposed_boundary/tables/Exposed_dynamicA_T1promoter.tsv', sep='\t', index=False)
