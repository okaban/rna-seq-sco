"""(X) Functional test: does per-BGC T1 GCCGGC methylation density predict cluster activation?

For each BGC (grouped by bgc_name), compute T1 GCCGGC site density (sites/kb over the
gene span) and the mean activation log2FC (T2vs1, T3vs1) of its member genes, then
correlate density vs activation across clusters (Spearman). A positive/negative
association would be a (weak) functional link; a null reinforces the sequence-artifact reading.
"""
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr

BASE = Path(__file__).resolve().parents[4]
GM = BASE / '05_annotation/analysis/05_annotation_260128_v1/tables/gene_master_with_BGC.tsv'
SITES = BASE / '11_epigenome_integration/analysis/37_defense_island_GCCGGC/tables/GCCGGC_sites_by_timepoint.tsv'

BGC = BASE / '11_epigenome_integration/analysis/72_antismash_genomewide/bgc_summary.tsv'
gm = pd.read_csv(GM, sep='\t')
sites = pd.read_csv(SITES, sep='\t')
bgc = pd.read_csv(BGC, sep='\t')
t1 = sites[sites['timepoint'] == 'T1']
print(f"antiSMASH regions: {len(bgc)} | genes: {len(gm)}")

# assign genes to each antiSMASH region by coordinate overlap; compute density + mean LFC
rows = []
for _, r in bgc.iterrows():
    start, end = r['start'], r['end']
    span_kb = (end - start) / 1000.0
    g = gm[(gm['start'] <= end) & (gm['end'] >= start)]   # overlap
    g = g.dropna(subset=['log2FoldChange_2_vs_1'])
    if span_kb < 1 or len(g) < 3:
        continue
    n_sites = ((t1['position'] >= start) & (t1['position'] <= end)).sum()
    rows.append({
        'region': r['region'], 'product': str(r.get('products', ''))[:20],
        'known': str(r.get('top_known_name', ''))[:18],
        'n_genes': len(g), 'span_kb': round(span_kb, 1),
        'GCCGGC_T1_density': round(n_sites / span_kb, 3),
        'mean_LFC_T2v1': round(g['log2FoldChange_2_vs_1'].mean(), 3),
        'mean_LFC_T3v1': round(g['log2FoldChange_3_vs_1'].mean(), 3),
        'max_LFC_T2v1': round(g['log2FoldChange_2_vs_1'].max(), 3),
    })
df = pd.DataFrame(rows).sort_values('GCCGGC_T1_density', ascending=False)
print(f"\nClusters analysed: {len(df)}")
print(df.to_string(index=False))

print("\n=== Spearman: T1 GCCGGC density vs activation ===")
for col in ['mean_LFC_T2v1', 'mean_LFC_T3v1', 'max_LFC_T2v1']:
    d = df.dropna(subset=['GCCGGC_T1_density', col])
    if len(d) >= 4:
        rho, p = spearmanr(d['GCCGGC_T1_density'], d[col])
        print(f"  density vs {col}: rho={rho:.3f}, p={p:.3f}  (n={len(d)})")

out = BASE / '11_epigenome_integration/analysis/47_BGC_methylation_geographic_test/tables/H24c_BGC_density_vs_activation.tsv'
df.to_csv(out, sep='\t', index=False)
print(f"\nSaved: {out}")
