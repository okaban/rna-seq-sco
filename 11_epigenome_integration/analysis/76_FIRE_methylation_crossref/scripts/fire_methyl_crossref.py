"""FIRE (Deng 2023 PNAS, sd04) x GCCGGC m4C methylation cross-reference.

Tests whether the vegetative methylome tracks 3D local interaction (FIRE), i.e.
whether methylation can serve as a Hi-C-free proxy for high-expression
integration neighbourhoods (Deng's FIRE-guided strategy).

Test A: genome-wide 5kb-bin Spearman(FIRE_M, GCCGGC-T1 density), raw + region-controlled.
        secondary: FIRE_L vs T2 methylation.
Cross : reproduce compartment-A enrichment vs Deng's actual PC eigenvector.
Test B: do Exposed-62 promoters sit in higher-FIRE bins than other regulators?

Deng M-phase ~14h ~ our T1 (12h); L-phase ~22h ~ our T2 (24h).
"""
import numpy as np, pandas as pd
from pathlib import Path
from scipy.stats import spearmanr, mannwhitneyu, fisher_exact, pearsonr

BASE = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis")
D = BASE / "76_FIRE_methylation_crossref"
BIN = 5000
CORE_LO, CORE_HI = 2_300_000, 6_200_000   # Deng compartment-A central core

def load_fire(sheet):
    f = pd.read_excel(D/"data/pnas.2222045120.sd04.xlsx", sheet_name=sheet, skiprows=1)
    f.columns = ['start', 'end', 'PC', 'FIRE', 'note']
    f = f[pd.to_numeric(f['start'], errors='coerce').notna()].copy()
    f = f.astype({'start': float, 'end': float, 'PC': float, 'FIRE': float})
    f['bin'] = (f['start'] // BIN).astype(int)
    return f[['bin', 'start', 'PC', 'FIRE']]

fire_m = load_fire("M-phase").rename(columns={'PC': 'PC_M', 'FIRE': 'FIRE_M'})
fire_l = load_fire("L-phase").rename(columns={'PC': 'PC_L', 'FIRE': 'FIRE_L'})
fire = fire_m.merge(fire_l[['bin', 'PC_L', 'FIRE_L']], on='bin', how='outer')
fire['center'] = fire['bin'] * BIN + BIN/2
fire['core'] = ((fire['center'] >= CORE_LO) & (fire['center'] <= CORE_HI)).astype(int)

# --- methylation density per 5kb bin ---
# 2026-09-21 (BLOCKER-0): the 37_ table is position-deduplicated across timepoints
# (its T2 rows = sites NEW at T2, n = 407), which fed m_n_T2 and hence the retracted
# panel-c rho = 0.43. Read the canonical per-timepoint GCCGGC 4mC sets
# (1,289/1,595/1,073) through 90_/canonical_sites.py. T1 is unchanged.
import importlib.util as _ilu
_spec = _ilu.spec_from_file_location("canonical_sites", BASE/"90_per_timepoint_census_audit/canonical_sites.py")
_cs = _ilu.module_from_spec(_spec); _spec.loader.exec_module(_cs)
g = _cs.gccggc_by_timepoint()
def methyl_density(tp):
    s = g[g.timepoint == tp].copy()
    s['bin'] = (s['position'] // BIN).astype(int)
    d = s.groupby('bin').agg(n=('position', 'size'), freqsum=('frequency', 'sum')).reset_index()
    return d
mt1 = methyl_density("T1").rename(columns={'n': 'm_n_T1', 'freqsum': 'm_fq_T1'})
mt2 = methyl_density("T2").rename(columns={'n': 'm_n_T2', 'freqsum': 'm_fq_T2'})
df = fire.merge(mt1, on='bin', how='left').merge(mt2, on='bin', how='left')
for c in ['m_n_T1', 'm_fq_T1', 'm_n_T2', 'm_fq_T2']:
    df[c] = df[c].fillna(0.0)

def partial_spearman(x, y, z):
    """Spearman of x,y residualised on z (all rank-transformed)."""
    rx, ry, rz = x.rank(), y.rank(), z.rank()
    rx_r = rx - np.polyval(np.polyfit(rz, rx, 1), rz)
    ry_r = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
    return pearsonr(rx_r, ry_r)

print("="*70)
print(f"TEST A — genome-wide FIRE x methylation ({len(df)} 5kb bins)")
print("="*70)
for fcol, mcol, lab in [('FIRE_M','m_n_T1','FIRE_M vs GCCGGC-T1 count'),
                        ('FIRE_M','m_fq_T1','FIRE_M vs GCCGGC-T1 freq-sum'),
                        ('FIRE_L','m_n_T2','FIRE_L vs GCCGGC-T2 count (L~T2)')]:
    d = df.dropna(subset=[fcol, mcol])
    rho, p = spearmanr(d[fcol], d[mcol])
    rr, pp = partial_spearman(d[fcol], d[mcol], d['core'])
    print(f"  {lab}:  raw rho={rho:.3f} p={p:.2e} (n={len(d)})  | core/arm-controlled r={rr:.3f} p={pp:.2e}")

print("\n  PC (compartment eigenvector) vs methylation:")
d = df.dropna(subset=['PC_M'])
rho, p = spearmanr(d['PC_M'], d['m_n_T1']); print(f"    PC_M vs GCCGGC-T1 count: rho={rho:.3f} p={p:.2e}")

print("\n" + "="*70)
print("CROSS — compartment-A enrichment vs Deng PC & core definition")
print("="*70)
tot_sites = df['m_n_T1'].sum()
core_sites = df.loc[df.core == 1, 'm_n_T1'].sum()
core_bins_frac = df.core.mean()
print(f"  GCCGGC-T1 sites in core(2.3-6.2Mb): {core_sites:.0f}/{tot_sites:.0f} = {core_sites/tot_sites:.3f}"
      f"  (core spans {core_bins_frac:.3f} of bins)")
# Fisher: methylated-heavy bins vs core
df['methyl_bin'] = (df['m_n_T1'] > 0).astype(int)
ct = pd.crosstab(df['core'], df['methyl_bin'])
orr, pf = fisher_exact(ct.values)
print(f"  core vs has-methylation bin: OR={orr:.2f} p={pf:.2e}\n  crosstab:\n{ct}")
# FIRE in methylated vs unmethylated bins
mm = df[df.m_n_T1 > 0]['FIRE_M']; uu = df[df.m_n_T1 == 0]['FIRE_M']
u, pu = mannwhitneyu(mm, uu)
print(f"  FIRE_M methylated bins median={mm.median():.3f} (n={len(mm)}) vs unmethylated {uu.median():.3f} (n={len(uu)}), MWU p={pu:.2e}")

print("\n" + "="*70)
print("TEST B — Exposed-62 promoters vs FIRE")
print("="*70)
exp = pd.read_csv(BASE/"52_shielded_exposed_boundary/tables/Exposed_dynamicA_T1promoter.tsv", sep="\t")
reg = pd.read_csv(BASE/"52_shielded_exposed_boundary/tables/all_genes_features_unified_n57.tsv", sep="\t").dropna(subset=['tss'])
reg['tss'] = reg['tss'].astype(int)
exp_lt = set(exp['locus_tag'])
fire_by_bin = df.set_index('bin')['FIRE_M'].to_dict()
def tss_fire(tss):
    return fire_by_bin.get(int(tss)//BIN, np.nan)
reg['FIRE_tss'] = reg['tss'].apply(tss_fire)
reg['Exposed'] = reg['locus_tag'].isin(exp_lt)
E = reg[reg.Exposed].dropna(subset=['FIRE_tss']); S = reg[~reg.Exposed].dropna(subset=['FIRE_tss'])
u, p = mannwhitneyu(E['FIRE_tss'], S['FIRE_tss'])
print(f"  Exposed n={len(E)} FIRE median={E['FIRE_tss'].median():.3f}  vs  Shielded n={len(S)} median={S['FIRE_tss'].median():.3f}")
print(f"  Mann-Whitney p={p:.3e}")
# vs genome-wide
gw = df['FIRE_M'].median()
print(f"  genome-wide bin FIRE median={gw:.3f}")

df.to_csv(D/"tables/fire_methyl_bins.tsv", sep="\t", index=False)
reg[['locus_tag','tf_family','region','tss','FIRE_tss','Exposed','LFC_T2vsT1']].to_csv(
    D/"tables/regulatory_FIRE_at_TSS.tsv", sep="\t", index=False)
print(f"\n[saved] tables/fire_methyl_bins.tsv ({len(df)} bins), regulatory_FIRE_at_TSS.tsv")
