"""R6 canonical recompute (n=57, corrected blocs).

Recomputes every number used in the R6 / Figure 5 manuscript section from the
canonical n=57 data after the H34 bloc-assignment fix (data-driven by trajectory).
Blocs are the module-based labels in temporal_classification.tsv (matches Fig5).

Outputs: bloc sizes, per-sample eigengene correlation (+bootstrap CI, permutation p),
family composition per bloc, TetR tests (within-Exposed and vs genome background),
early-responder fraction.
"""
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import fisher_exact, binomtest, pearsonr

BASE = Path(__file__).resolve().parents[4]
TC = BASE / '11_epigenome_integration/analysis/57_temporal_dynamics_exposed_TF/tables/temporal_classification.tsv'
NORM = BASE / '04_deseq2/analysis/04_deseq2_260128_v1/results/normalized_counts_M145.tsv'
SAMPLES = ['M145_1_1', 'M145_1_2', 'M145_1_3',
           'M145_2_1', 'M145_2_3', 'M145_2_4',
           'M145_3_2', 'M145_3_3', 'M145_3_4']
SEED = 42
rng = np.random.default_rng(SEED)

t = pd.read_csv(TC, sep='\t')
nc = pd.read_csv(NORM, sep='\t').set_index('gene_id')

print("=== Bloc sizes (module-based, matches Fig5) ===")
print(t['bloc'].value_counts(dropna=False).to_dict())

# per-gene z-scored expression across the 9 samples
genes = [g for g in t['locus_tag'] if g in nc.index]
expr = np.log2(nc.loc[genes, SAMPLES].astype(float) + 1)
z = expr.sub(expr.mean(axis=1), axis=0).div(expr.std(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
bloc = t.set_index('locus_tag')['bloc']

def eigengene(gene_list):
    """PC1 (first right singular vector projection) of the bloc's z-matrix over 9 samples."""
    M = z.loc[[g for g in gene_list if g in z.index]].values  # genes x samples
    if M.shape[0] < 2:
        return None
    U, S, Vt = np.linalg.svd(M - M.mean(axis=0, keepdims=True), full_matrices=False)
    eg = Vt[0]                      # sample-space eigengene
    # orient so it tracks the mean trajectory
    if np.corrcoef(eg, M.mean(axis=0))[0, 1] < 0:
        eg = -eg
    return eg

act_genes = bloc[bloc == 'activation'].index.tolist()
rep_genes = bloc[bloc == 'repression'].index.tolist()
eg_a, eg_r = eigengene(act_genes), eigengene(rep_genes)
r_obs, p_obs = pearsonr(eg_a, eg_r)
print(f"\n=== Eigengene correlation (activation vs repression, 9 samples) ===")
print(f"  n_act={len(act_genes)}  n_rep={len(rep_genes)}")
print(f"  Pearson r = {r_obs:.4f}  (parametric p = {p_obs:.2e})")

# bootstrap CI over samples
boot = []
for _ in range(10000):
    idx = rng.integers(0, 9, 9)
    a, b = eg_a[idx], eg_r[idx]
    if np.std(a) > 0 and np.std(b) > 0:
        boot.append(np.corrcoef(a, b)[0, 1])
lo, hi = np.percentile(boot, [2.5, 97.5])
print(f"  bootstrap 95% CI = [{lo:.3f}, {hi:.3f}]  (n={len(boot)})")

# permutation null: shuffle bloc labels
labels = np.array(['activation'] * len(act_genes) + ['repression'] * len(rep_genes))
pool = act_genes + rep_genes
perm = []
for _ in range(10000):
    sh = rng.permutation(labels)
    a = eigengene([g for g, l in zip(pool, sh) if l == 'activation'])
    b = eigengene([g for g, l in zip(pool, sh) if l == 'repression'])
    if a is not None and b is not None:
        perm.append(np.corrcoef(a, b)[0, 1])
perm = np.array(perm)
p_perm = (np.sum(perm <= r_obs) + 1) / (len(perm) + 1)
print(f"  permutation p (r <= observed) = {p_perm:.4f}")

print("\n=== Family composition per bloc ===")
comp = pd.crosstab(t['tf_family'], t['bloc'])
print(comp.to_string())

print("\n=== TetR tests ===")
tetr = t[t['tf_family'].str.contains('TetR', case=False, na=False)]
ta = ((tetr['bloc'] == 'activation')).sum()
tr = ((tetr['bloc'] == 'repression')).sum()
na = (t['bloc'] == 'activation').sum()
nr = (t['bloc'] == 'repression').sum()
odd, p_within = fisher_exact([[tr, nr - tr], [ta, na - ta]])
print(f"  TetR: {ta} activation / {tr} repression")
print(f"  repression bloc TetR fraction: {tr}/{nr} = {tr/nr:.1%}")
print(f"  Test A (repression vs activation, Fisher): OR={odd:.3f}, p={p_within:.4f}")
for bg in (0.04,):
    bt = binomtest(tr, nr, bg, alternative='greater')
    print(f"  Test B (repression {tr}/{nr} vs genome {bg:.0%}, binomial): p={bt.pvalue:.5f}")

print("\n=== Early responder fraction (phase_ratio > 0.6) ===")
if 'phase_ratio' in t.columns:
    er = (t['phase_ratio'] > 0.6).sum()
    print(f"  {er}/{len(t)} = {er/len(t):.0%}")
