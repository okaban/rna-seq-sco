#!/usr/bin/env python3
"""Generate Figure4_shielded_exposed.png at 300 dpi (n=57 update)."""
import warnings; warnings.filterwarnings('ignore')
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.metrics import roc_curve, auc

BASE      = Path(__file__).resolve().parents[2]
EPIGENOME = BASE / '11_epigenome_integration/analysis'
FIG_DIR   = BASE / '15_paper_figures/figures/main'
TBLS      = EPIGENOME / '52_shielded_exposed_boundary/tables'

COL_EXPOSED  = '#7E57C2'
COL_SHIELDED = '#B0BEC5'

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
    'axes.titlesize':10,'axes.labelsize':10,'xtick.labelsize':9,
    'ytick.labelsize':9,'legend.fontsize':8,'figure.dpi':300,
    'savefig.dpi':300,'axes.linewidth':1.0,
    'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})

df   = pd.read_csv(TBLS/'all_genes_features_unified_n57.tsv', sep='\t'
                   ).dropna(subset=['nearest_methyl_distance'])
df_q = pd.read_csv(TBLS/'expression_quintile.tsv', sep='\t')

print(f"Data: {len(df)} genes, Exposed={int(df.is_exposed.sum())}, Shielded={(df.is_exposed==0).sum()}")

shi = df.loc[df.is_exposed==0,'nearest_methyl_distance'].values
exp = df.loc[df.is_exposed==1,'nearest_methyl_distance'].values

fig, axes = plt.subplots(1, 3, figsize=(7.09, 3.07))
fig.subplots_adjust(left=0.09, right=0.97, top=0.88, bottom=0.20, wspace=0.50)

# ── Panel A: violin ──
ax = axes[0]
vp = ax.violinplot([shi, exp], positions=[1,2], showmedians=False, showextrema=False)
for body, col in zip(vp['bodies'], [COL_SHIELDED, COL_EXPOSED]):
    body.set_facecolor(col); body.set_alpha(0.4); body.set_edgecolor('none')
rng = np.random.default_rng(42)
for i, (d, pos, col) in enumerate([(shi,1,COL_SHIELDED),(exp,2,COL_EXPOSED)]):
    q25, med, q75 = np.percentile(d, [25, 50, 75])
    ax.plot([pos-.08,pos+.08], [med,med],  color=col, lw=2.5, zorder=5)
    ax.plot([pos-.05,pos+.05], [q25,q25],  color=col, lw=1.2, zorder=5)
    ax.plot([pos-.05,pos+.05], [q75,q75],  color=col, lw=1.2, zorder=5)
    ax.plot([pos,pos],         [q25,q75],  color=col, lw=1.0, zorder=4)
n_strip = min(150, len(shi))
idx_s = rng.choice(len(shi), n_strip, replace=False)
ax.scatter(rng.normal(1,.05,n_strip), shi[idx_s], color=COL_SHIELDED, s=3, alpha=.35, zorder=3, linewidths=0)
ax.scatter(rng.normal(2,.05,len(exp)), exp,        color=COL_EXPOSED,  s=6, alpha=.70, zorder=3, linewidths=0)
ax.axhline(293, color='#E53935', lw=1.2, ls='--', alpha=.9, zorder=6)
ax.text(2.35, 293, '293 bp', color='#E53935', fontsize=7.5, va='center', fontweight='bold')
ax.set_ylim(-500, max(shi.max(), exp.max())*1.05)
ax.set_xticks([1,2])
ax.set_xticklabels([f'Shielded\n(n={len(shi)})', f'Exposed\n(n={len(exp)})'], fontsize=9)
ax.set_ylabel('Nearest methylation site\ndistance from TSS (bp)', fontsize=9)
ax.set_title('Promoter methylation\nexposure', fontsize=9, fontweight='bold')
ax.text(0.05, 1.10, 'a', transform=ax.transAxes, fontsize=13, fontweight='bold', va='top')

# ── Panel B: ROC ──
ax = axes[1]
valid = ~np.isnan(df.nearest_methyl_distance.values) & ~np.isnan(df.baseMean.values)
y  = df.is_exposed.values[valid]
d  = df.nearest_methyl_distance.values[valid]
b  = df.baseMean.values[valid]
fpr_d, tpr_d, _ = roc_curve(y, -d); auc_d = auc(fpr_d, tpr_d)
fpr_b, tpr_b, _ = roc_curve(y, -b); auc_b = auc(fpr_b, tpr_b)
ax.plot(fpr_d, tpr_d, color=COL_EXPOSED, lw=2.0, label=f'Nearest distance\nAUC = {auc_d:.3f}')
ax.plot(fpr_b, tpr_b, color=COL_SHIELDED, lw=1.5, ls='--', label=f'Expression level\nAUC = {auc_b:.3f}')
ax.plot([0,1],[0,1],'k:',lw=.6,alpha=.5)
op_fpr = 1.0-0.802
ax.scatter([op_fpr],[1.0], color='#E53935', s=50, zorder=5, marker='D', linewidths=0)
ax.annotate('293 bp\n(sens=1.00\nspec=0.80)',
    xy=(op_fpr,1.0), xytext=(op_fpr+.10,1.0-.18),
    fontsize=6.5, color='#E53935',
    arrowprops=dict(arrowstyle='->', color='#E53935', lw=.8))
ax.set_xlabel('False positive rate', fontsize=9)
ax.set_ylabel('True positive rate',  fontsize=9)
ax.set_title('Exposed/Shielded\nclassification (ROC)', fontsize=9, fontweight='bold')
ax.legend(fontsize=7, loc='lower right', frameon=True, fancybox=False,
          edgecolor='#ccc', labelspacing=.5)
ax.set_xlim(-0.02,1.02); ax.set_ylim(-0.02,1.02); ax.set_aspect('equal')
ax.text(-0.22, 1.10, 'b', transform=ax.transAxes, fontsize=13, fontweight='bold', va='top')

# ── Panel C: quintile ──
ax = axes[2]
q    = df_q.quintile.values
frac = df_q.frac_exposed.values * 100
ax.bar(q, frac, color=COL_EXPOSED, alpha=.7, edgecolor='white', lw=.5, width=.6)
for qi, fi in zip(q, frac):
    ax.text(qi, fi+.3, f'{fi:.1f}%', ha='center', va='bottom', fontsize=7)
sl, ic, *_ = stats.linregress(q, frac)
xf = np.array([.5, 5.5])
ax.plot(xf, sl*xf+ic, 'k--', lw=.8, alpha=.5)
ax.axhline(frac.mean(), color=COL_EXPOSED, lw=.8, ls=':', alpha=.6)
ax.text(0.97, 0.95,
    'Jonckheere–Terpstra\np = 0.730 (NS)\n\nbaseMean AUC = 0.543',
    transform=ax.transAxes, ha='right', va='top', fontsize=7,
    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#ccc', alpha=.9))
ax.set_xlabel('Expression quintile\n(1 = lowest, 5 = highest)', fontsize=9)
ax.set_ylabel('Exposed genes (%)', fontsize=9)
ax.set_ylim(0, max(frac)*1.55)
ax.set_xticks(q)
ax.set_title('Expression-independent\nclassification', fontsize=9, fontweight='bold')
ax.text(-0.22, 1.10, 'c', transform=ax.transAxes, fontsize=13, fontweight='bold', va='top')

for fmt in ('pdf', 'svg', 'png'):
    out = FIG_DIR / f'Figure4_shielded_exposed.{fmt}'
    fig.savefig(out, format=fmt, dpi=300, bbox_inches='tight')
    print(f'  Saved: {out}')
plt.close(fig)
print('Done.')
