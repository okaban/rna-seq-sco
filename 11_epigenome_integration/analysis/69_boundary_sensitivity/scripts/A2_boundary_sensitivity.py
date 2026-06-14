#!/usr/bin/env python3
"""
A-2: Protection Zone Boundary Sensitivity + T2/T3 Independent Classifier Models
=================================================================================

Key facts established from H29 analysis:
  - "Exposed" genes (n=11) have methylation CLOSE to TSS (median 122 bp)
    → methylation in promoter = gene expression responds to methylation changes
  - "Shielded" genes (n=353) have methylation FAR from TSS (median 819 bp)
    → no methylation in promoter = expression is insensitive to methylation
  - Original H29 used ALL-timepoint combined methylation data + 5-fold stratified CV
  - AUC=0.9233 (raw), 5-fold CV AUC=0.9166  (feature: -nearest_methyl_distance)
  - 293 bp = Youden-optimal threshold on nearest_methyl_distance (all-TP combined)

Analysis 1: Boundary sensitivity (100-500 bp sweep)
  Feature: count of (all-TP combined) methylation sites within window upstream of TSS
  Higher count → more likely exposed
  Method: 5-fold stratified CV (consistent with H29)
  Goal: verify 293 bp is stable / near-optimal

Analysis 2: T2/T3 independent classifier models
  Feature: nearest_methyl_distance computed from per-timepoint methylation sites
  Method: 5-fold stratified CV
  Goal: compare T1/T2/T3 independent AUCs

Output:
  tables/A2_boundary_sensitivity_AUC.tsv
  tables/A2_timepoint_models_AUC.tsv
  figures/A2_boundary_sensitivity.png
  figures/A2_timepoint_ROC.png
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, roc_curve
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Paths
# ============================================================
BASE = Path('/Users/okaban/bioinfo/rna-seq')
OUT_DIR = BASE / '11_epigenome_integration/analysis/69_boundary_sensitivity'
TBL_DIR = OUT_DIR / 'tables'
FIG_DIR = OUT_DIR / 'figures'

FEATURES_TSV = BASE / '11_epigenome_integration/analysis/52_shielded_exposed_boundary/tables/all_genes_features.tsv'
METHYL_CSV   = BASE / '11_epigenome_integration/analysis/01_integration/high_confidence_sites_weighted.csv'

SEED = 42
N_SPLITS = 5
np.random.seed(SEED)

BOUNDARIES = [100, 150, 200, 250, 293, 350, 400, 450, 500]
TIMEPOINTS  = ['T1', 'T2', 'T3']

# ============================================================
# Load data
# ============================================================
print("=" * 70)
print("A-2: Boundary Sensitivity + Timepoint Independent Models")
print("=" * 70)

genes  = pd.read_csv(FEATURES_TSV, sep='\t')
methyl = pd.read_csv(METHYL_CSV)

n_exp  = int(genes['is_exposed'].sum())
n_shi  = int((genes['is_exposed'] == 0).sum())
print(f"Genes: {len(genes)} total  (exposed={n_exp}, shielded={n_shi})")
print(f"Methylation sites: {len(methyl)}  (T1={len(methyl[methyl.timepoint=='T1'])}, "
      f"T2={len(methyl[methyl.timepoint=='T2'])}, T3={len(methyl[methyl.timepoint=='T3'])})")

y = genes['is_exposed'].values

# Note on feature direction:
# - Exposed genes have LOWER nearest_methyl_distance (closer methylation to TSS)
# - For logistic regression: use -nearest_methyl_distance OR window-count
# - Window-count: higher count within window → more likely exposed

# ============================================================
# Helper functions
# ============================================================

def count_sites_in_window(tss: int, strand: str, positions: np.ndarray, window: int) -> int:
    """Count methylation positions in promoter window upstream of TSS."""
    if strand == '+':
        lo, hi = tss - window, tss + 50
    else:
        lo, hi = tss - 50,     tss + window
    return int(((positions >= lo) & (positions <= hi)).sum())


def nearest_methyl_dist(tss: int, positions: np.ndarray) -> float:
    if len(positions) == 0:
        return np.nan
    return float(np.abs(positions - tss).min())


def cv5_auc(X: np.ndarray, y: np.ndarray) -> tuple[float, np.ndarray]:
    """5-fold stratified CV AUC; returns (mean_AUC, predicted_probs)."""
    skf = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=SEED)
    probs = np.zeros(len(y))
    clf = LogisticRegression(max_iter=1000, class_weight='balanced', random_state=SEED)
    for tr, te in skf.split(X, y):
        clf.fit(X[tr], y[tr])
        probs[te] = clf.predict_proba(X[te])[:, 1]
    auc_val = roc_auc_score(y, probs)
    return auc_val, probs

# ============================================================
# Analysis 1: Boundary Sensitivity (ALL-timepoints combined sites)
# ============================================================
print("\n" + "=" * 70)
print("Analysis 1: Boundary Sensitivity Sweep")
print("  (feature: site count in window, ALL timepoints combined, 5-fold CV)")
print("=" * 70)

all_positions = np.sort(methyl['position'].unique())
print(f"All-timepoint unique positions: {len(all_positions)}")

rows_sens = []
probs_per_boundary = {}

for bp in BOUNDARIES:
    counts = genes.apply(
        lambda r: count_sites_in_window(int(r['tss']), r['strand'], all_positions, bp),
        axis=1
    ).values.astype(float).reshape(-1, 1)

    n_nonzero  = int((counts.ravel() > 0).sum())
    n_exp_nonzero = int((counts.ravel() > 0)[y == 1].sum())
    n_shi_nonzero = int((counts.ravel() > 0)[y == 0].sum())

    # Raw AUC (no CV)
    raw_auc = roc_auc_score(y, counts.ravel())

    # 5-fold CV AUC
    auc5, probs = cv5_auc(counts, y)
    probs_per_boundary[bp] = probs

    rows_sens.append({
        'boundary_bp':        bp,
        'raw_AUC':            round(raw_auc, 4),
        'cv5_AUC':            round(auc5, 4),
        'n_genes_with_site':  n_nonzero,
        'n_exposed_with_site':n_exp_nonzero,
        'n_shielded_with_site':n_shi_nonzero,
        'mean_count_exposed': round(float(counts[y == 1].mean()), 3),
        'mean_count_shielded':round(float(counts[y == 0].mean()), 3),
    })
    print(f"  {bp:3d} bp | raw_AUC={raw_auc:.4f}  CV5_AUC={auc5:.4f} | "
          f"sites: exp={n_exp_nonzero}/{n_exp}  shi={n_shi_nonzero}/{n_shi}")

df_sens = pd.DataFrame(rows_sens)
df_sens.to_csv(TBL_DIR / 'A2_boundary_sensitivity_AUC.tsv', sep='\t', index=False)

best_row  = df_sens.loc[df_sens['cv5_AUC'].idxmax()]
row_293   = df_sens[df_sens['boundary_bp'] == 293].iloc[0]
auc_range = df_sens['cv5_AUC'].max() - df_sens['cv5_AUC'].min()

print(f"\nBest window:  {best_row['boundary_bp']:.0f} bp  (CV5_AUC={best_row['cv5_AUC']:.4f})")
print(f"293 bp (orig): CV5_AUC={row_293['cv5_AUC']:.4f}")
print(f"AUC range (100-500 bp): {auc_range:.4f}")

# ============================================================
# Analysis 2: Timepoint-independent Classifier Models
# ============================================================
print("\n" + "=" * 70)
print("Analysis 2: Timepoint-independent Models")
print("  (feature: -nearest_methyl_distance per timepoint, 5-fold CV)")
print("=" * 70)

rows_tp = []
roc_data = {}

# Reference: all-timepoints combined (= H29 original)
X_all = genes['nearest_methyl_distance'].values.reshape(-1, 1)
X_all_neg = -X_all.copy()
auc_all, probs_all = cv5_auc(X_all_neg, y)
print(f"  All-TP (reference, H29-consistent):  CV5_AUC={auc_all:.4f}")
fpr, tpr, _ = roc_curve(y, probs_all)
roc_data['All-TP'] = (fpr, tpr, auc_all)

for tp in TIMEPOINTS:
    pos_tp = np.sort(methyl[methyl['timepoint'] == tp]['position'].values)
    dists = genes['tss'].apply(
        lambda t: nearest_methyl_dist(int(t), pos_tp)
    ).values

    n_nan = int(np.isnan(dists).sum())
    if n_nan > 0:
        med = np.nanmedian(dists)
        dists[np.isnan(dists)] = med

    # Flip sign: lower distance = more exposed → higher -distance = more exposed
    X = (-dists).reshape(-1, 1)
    raw_auc = roc_auc_score(y, X.ravel())
    auc5, probs = cv5_auc(X, y)

    dist_exp = dists[y == 1]
    dist_shi = dists[y == 0]
    mwu_stat, mwu_p = stats.mannwhitneyu(dist_exp, dist_shi, alternative='two-sided')

    fpr, tpr, _ = roc_curve(y, probs)
    roc_data[tp] = (fpr, tpr, auc5)

    rows_tp.append({
        'timepoint':              tp,
        'raw_AUC':                round(raw_auc, 4),
        'cv5_AUC':                round(auc5, 4),
        'exposed_median_dist_bp': round(float(np.median(dist_exp)), 1),
        'shielded_median_dist_bp':round(float(np.median(dist_shi)), 1),
        'MWU_p':                  f"{mwu_p:.3e}",
        'n_methyl_sites':         int(len(pos_tp)),
        'n_nan_imputed':          n_nan,
    })
    print(f"  {tp}: raw_AUC={raw_auc:.4f}  CV5_AUC={auc5:.4f} | "
          f"exp_median={np.median(dist_exp):.0f}  shi_median={np.median(dist_shi):.0f} bp  "
          f"MWU_p={mwu_p:.2e}")

df_tp = pd.DataFrame(rows_tp)
df_tp.to_csv(TBL_DIR / 'A2_timepoint_models_AUC.tsv', sep='\t', index=False)

t1_auc = df_tp[df_tp['timepoint'] == 'T1']['cv5_AUC'].iloc[0]
t2_auc = df_tp[df_tp['timepoint'] == 'T2']['cv5_AUC'].iloc[0]
t3_auc = df_tp[df_tp['timepoint'] == 'T3']['cv5_AUC'].iloc[0]

# ============================================================
# Figures
# ============================================================
print("\n--- Generating Figures ---")

BLUE  = '#2C7BB6'
RED   = '#D7191C'
GREEN = '#1A9641'
GOLD  = '#FDAE61'
GRAY  = '#636363'

# ---- Figure 1: Boundary Sensitivity ----
fig1, axes1 = plt.subplots(1, 2, figsize=(13, 5))
fig1.suptitle('A-2: Protection Zone Boundary Sensitivity\n'
              '(All-timepoint methylation, 5-fold CV, window-count feature)',
              fontsize=12, fontweight='bold')

bps  = df_sens['boundary_bp'].values
aucs = df_sens['cv5_AUC'].values
raw_aucs = df_sens['raw_AUC'].values

# Panel A: AUC vs window
ax = axes1[0]
ax.plot(bps, raw_aucs, 's--', color=BLUE, lw=1.5, ms=7, alpha=0.5, label='Raw AUC')
ax.plot(bps, aucs, 'o-',  color=BLUE, lw=2.5, ms=8, zorder=3,    label='5-fold CV AUC')

# Highlight 293 bp
idx_293 = list(bps).index(293)
ax.axvline(x=293, color=RED, ls='--', lw=1.5, alpha=0.8, label=f'293 bp (original)')
ax.scatter([293], [aucs[idx_293]], color=RED, s=120, zorder=5)
ax.annotate(f'293 bp\nCV={aucs[idx_293]:.3f}',
            xy=(293, aucs[idx_293]),
            xytext=(310, aucs[idx_293] - 0.03),
            fontsize=9, color=RED,
            arrowprops=dict(arrowstyle='->', color=RED, lw=1.2))

best_idx = np.argmax(aucs)
if best_idx != idx_293:
    ax.scatter([bps[best_idx]], [aucs[best_idx]], color=GREEN, s=150, zorder=6,
               marker='*', label=f'Best: {bps[best_idx]} bp (CV={aucs[best_idx]:.3f})')

# All-timepoint distance AUC reference line (live cv5 value; was stale n=62 0.9166 hardcode)
ax.axhline(y=auc_all, color=GREEN, ls=':', lw=1.5, alpha=0.7, label=f'all-TP dist (cv5): {auc_all:.3f}')

ax.set_xlabel('Promoter window (bp)', fontsize=12)
ax.set_ylabel('AUC', fontsize=12)
ax.set_title('A. AUC vs Window Size', fontsize=12, fontweight='bold')
ax.set_xticks(bps)
ax.set_xticklabels([str(b) for b in bps], rotation=45, fontsize=9)
ax.set_ylim([0.5, 1.02])
ax.legend(fontsize=8, loc='lower right')
ax.grid(axis='y', alpha=0.3)

# Panel B: site coverage
ax2 = axes1[1]
ax2b = ax2.twinx()

width = 35
exp_counts = df_sens['n_exposed_with_site'].values
shi_counts = df_sens['n_shielded_with_site'].values
x_pos = np.arange(len(bps))

bar_exp = ax2.bar(x_pos - width/200, exp_counts, width=0.35,
                  color=RED, alpha=0.7, label='Exposed', edgecolor='gray', lw=0.5)
bar_shi = ax2.bar(x_pos + width/200, shi_counts, width=0.35,
                  color=BLUE, alpha=0.5, label='Shielded', edgecolor='gray', lw=0.5)
ax2b.plot(x_pos, aucs, 'o-', color='black', lw=2, ms=7, zorder=3, label='CV5 AUC')

ax2.set_xlabel('Promoter window (bp)', fontsize=12)
ax2.set_ylabel('# Genes with ≥1 site in window', fontsize=11)
ax2b.set_ylabel('5-fold CV AUC', fontsize=11, color='black')
ax2.set_xticks(x_pos)
ax2.set_xticklabels([str(b) for b in bps], rotation=45, fontsize=9)
ax2.set_title('B. Site Coverage vs Window', fontsize=12, fontweight='bold')
ax2.legend(fontsize=9, loc='upper left')
ax2b.set_ylim([0.5, 1.02])
ax2.grid(axis='y', alpha=0.3)

fig1.tight_layout()
fig1.savefig(FIG_DIR / 'A2_boundary_sensitivity.png', dpi=300, bbox_inches='tight')
plt.close(fig1)
print("  Saved: A2_boundary_sensitivity.png")

# ---- Figure 2: Timepoint ROC ----
COLORS_TP = {'T1': '#4C72B0', 'T2': '#DD8452', 'T3': '#2CA02C', 'All-TP': '#7F7F7F'}

fig2, axes2 = plt.subplots(1, 2, figsize=(13, 5))
fig2.suptitle('A-2: Timepoint-Independent Classifier Models\n'
              '(feature: nearest methylation distance per timepoint, 5-fold CV)',
              fontsize=12, fontweight='bold')

# Panel A: ROC curves
ax3 = axes2[0]
for label, (fpr, tpr, auc_v) in roc_data.items():
    ls = ':' if label == 'All-TP' else '-'
    lw = 1.5 if label == 'All-TP' else 2.5
    ax3.plot(fpr, tpr, ls=ls, lw=lw, color=COLORS_TP.get(label, GRAY),
             label=f'{label} (AUC={auc_v:.3f})')
ax3.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.4)
ax3.set_xlabel('False Positive Rate', fontsize=12)
ax3.set_ylabel('True Positive Rate', fontsize=12)
ax3.set_title('A. ROC Curves by Timepoint', fontsize=12, fontweight='bold')
ax3.legend(fontsize=10, loc='lower right')
ax3.set_xlim([-0.02, 1.02])
ax3.set_ylim([-0.02, 1.02])
ax3.set_aspect('equal')

# Panel B: AUC bar
ax4 = axes2[1]
labels_bar = list(df_tp['timepoint'].values) + ['All-TP\n(H29 orig)']
aucs_bar   = list(df_tp['cv5_AUC'].values) + [auc_all]
colors_bar = [COLORS_TP.get(t, GRAY) for t in df_tp['timepoint']] + [GRAY]
hatches    = ['', '', '', '///']

bars = ax4.bar(range(len(labels_bar)), aucs_bar, color=colors_bar,
               alpha=0.8, edgecolor='gray', lw=1, width=0.5)
for bar, h in zip(bars, hatches):
    bar.set_hatch(h)
for i, (bar, auc_v) in enumerate(zip(bars, aucs_bar)):
    ax4.text(bar.get_x() + bar.get_width() / 2,
             bar.get_height() + 0.008,
             f'{auc_v:.3f}', ha='center', va='bottom', fontsize=11, fontweight='bold')

# Add cross-timepoint reference (T1 model applied to T2/T3)
ax4.axhline(y=0.644, ls='--', lw=1.5, color='orange', alpha=0.8,
            label='T1→T2 cross-TP: 0.644')
ax4.axhline(y=0.533, ls='--', lw=1.5, color='tomato', alpha=0.8,
            label='T1→T3 cross-TP: 0.533')
ax4.axhline(y=0.7,   ls=':',  lw=1,   color='gray', alpha=0.5)

ax4.set_xlabel('Timepoint (independent model)', fontsize=12)
ax4.set_ylabel('5-fold CV AUC', fontsize=12)
ax4.set_title('B. CV AUC per Timepoint', fontsize=12, fontweight='bold')
ax4.set_xticks(range(len(labels_bar)))
ax4.set_xticklabels(labels_bar, fontsize=10)
ax4.set_ylim([0, 1.05])
ax4.legend(fontsize=8, loc='lower right')
ax4.grid(axis='y', alpha=0.3)

fig2.tight_layout()
fig2.savefig(FIG_DIR / 'A2_timepoint_ROC.png', dpi=300, bbox_inches='tight')
plt.close(fig2)
print("  Saved: A2_timepoint_ROC.png")

# ============================================================
# Final Summary
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print("\nAnalysis 1 — Boundary Sensitivity (all-TP sites, window-count, 5-fold CV):")
print(df_sens[['boundary_bp', 'raw_AUC', 'cv5_AUC',
               'n_exposed_with_site', 'n_shielded_with_site']].to_string(index=False))
print(f"\n  293 bp (original): CV5_AUC={row_293['cv5_AUC']:.4f}")
print(f"  Best window:       CV5_AUC={best_row['cv5_AUC']:.4f} @ {best_row['boundary_bp']:.0f} bp")
print(f"  AUC range (100-500): {auc_range:.4f}")

if auc_range < 0.03:
    stability = "STABLE (plateau ≤ 0.03) — 293 bp boundary is robust"
elif auc_range < 0.08:
    stability = "MODERATELY SENSITIVE (range 0.03-0.08)"
else:
    stability = "SENSITIVE (range > 0.08) — 293 bp is not uniquely optimal"
print(f"  Stability: {stability}")

print("\nAnalysis 2 — Timepoint-Independent Models (nearest_methyl_distance, 5-fold CV):")
print(df_tp[['timepoint', 'raw_AUC', 'cv5_AUC', 'exposed_median_dist_bp',
             'shielded_median_dist_bp', 'MWU_p']].to_string(index=False))
print(f"\n  All-TP combined (H29 reference):  CV5_AUC={auc_all:.4f}")
print(f"  T1 independent:                   CV5_AUC={t1_auc:.4f}  (Δ={t1_auc - auc_all:+.4f} vs All-TP)")
print(f"  T2 independent:                   CV5_AUC={t2_auc:.4f}  (Δ={t2_auc - auc_all:+.4f} vs All-TP)")
print(f"  T3 independent:                   CV5_AUC={t3_auc:.4f}  (Δ={t3_auc - auc_all:+.4f} vs All-TP)")
print(f"\n  Cross-timepoint context (T1 model applied to T2/T3 methylation):")
print(f"    T1→T2 cross-TP: AUC=0.644  |  T2 independent: {t2_auc:.4f}  (Δ={t2_auc - 0.644:+.4f})")
print(f"    T1→T3 cross-TP: AUC=0.533  |  T3 independent: {t3_auc:.4f}  (Δ={t3_auc - 0.533:+.4f})")

print("\nOutput files:")
for f in sorted(TBL_DIR.glob('A2_*.tsv')):
    print(f"  {f}")
for f in sorted(FIG_DIR.glob('A2_*.png')):
    print(f"  {f}")
print("=" * 70)
