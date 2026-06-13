#!/usr/bin/env python3
"""
Figure 4: Shielded/Exposed Binary Classification and 293 bp Boundary

(A) Violin + strip plot — nearest_methyl_distance (Shielded vs Exposed)
(B) ROC curves — distance AUC=0.917 vs expression AUC=0.547 (+ 293bp threshold)
(C) Expression-independence — exposed fraction by expression quintile
"""

import importlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from scipy import stats
from sklearn.metrics import roc_curve, auc

_utils = importlib.import_module('00_shared_utils')
for _attr in dir(_utils):
    if not _attr.startswith('_'):
        globals()[_attr] = getattr(_utils, _attr)


def panel_a(ax, df_genes):
    """Panel A: Violin + strip plot of nearest_methyl_distance (Shielded vs Exposed)."""
    shielded = df_genes.loc[df_genes['is_exposed'] == 0, 'nearest_methyl_distance'].dropna()
    exposed = df_genes.loc[df_genes['is_exposed'] == 1, 'nearest_methyl_distance'].dropna()

    data = [shielded.values, exposed.values]
    labels = [f'Shielded\n(n={len(shielded)})', f'Exposed\n(n={len(exposed)})']
    colors = [COL_SHIELDED, COL_EXPOSED]

    # Violin
    vp = ax.violinplot(data, positions=[1, 2], showmedians=False, showextrema=False)
    for i, body in enumerate(vp['bodies']):
        body.set_facecolor(colors[i])
        body.set_alpha(0.4)
        body.set_edgecolor('none')

    # Box overlay
    for i, (d, pos) in enumerate(zip(data, [1, 2])):
        q25, med, q75 = np.percentile(d, [25, 50, 75])
        ax.plot([pos - 0.08, pos + 0.08], [med, med], color=colors[i], linewidth=2.5, zorder=5)
        ax.plot([pos - 0.05, pos + 0.05], [q25, q25], color=colors[i], linewidth=1.2, zorder=5)
        ax.plot([pos - 0.05, pos + 0.05], [q75, q75], color=colors[i], linewidth=1.2, zorder=5)
        ax.plot([pos, pos], [q25, q75], color=colors[i], linewidth=1.0, zorder=4)

    # Strip plot (subsample shielded to avoid overplotting)
    rng = np.random.default_rng(42)
    n_strip = min(150, len(shielded))
    idx_s = rng.choice(len(shielded), n_strip, replace=False)
    xs = rng.normal(1, 0.05, n_strip)
    ax.scatter(xs, shielded.values[idx_s], color=COL_SHIELDED, s=3, alpha=0.35, zorder=3, linewidths=0)

    xe = rng.normal(2, 0.05, len(exposed))
    ax.scatter(xe, exposed.values, color=COL_EXPOSED, s=6, alpha=0.7, zorder=3, linewidths=0)

    # 293 bp threshold line
    ax.axhline(293, color='#E53935', linewidth=1.2, linestyle='--', alpha=0.9, zorder=6)
    ax.text(2.35, 293, '293 bp', color='#E53935', fontsize=7.5, va='center', fontweight='bold')

    # Median annotations
    for i, (d, pos) in enumerate(zip(data, [1, 2])):
        med = np.median(d)
        ax.text(pos, -420, f'median\n{med:.0f} bp', ha='center', va='top',
                fontsize=7, color=colors[i])

    # Mann-Whitney p annotation
    stat, p = stats.mannwhitneyu(shielded, exposed, alternative='greater')
    p_str = f'p = {p:.1e}' if p < 0.001 else f'p = {p:.4f}'
    ax.text(1.5, ax.get_ylim()[1] * 0.95 if ax.get_ylim()[1] > 1 else 3000,
            p_str, ha='center', va='top', fontsize=7.5,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#ccc', alpha=0.9))

    ax.set_xticks([1, 2])
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel('Nearest methylation site\ndistance from TSS (bp)', fontsize=9)
    ax.set_title('Promoter methylation\nexposure', fontsize=9, fontweight='bold')
    ax.set_xlim(0.5, 2.5)
    ax.set_ylim(-500, None)

    # Set y-axis to show actual range
    y_max = max(shielded.max(), exposed.max())
    ax.set_ylim(-500, y_max * 1.05)


def panel_b(ax, df_genes):
    """Panel B: ROC curves — nearest distance (AUC=0.908) vs expression (AUC=0.547)."""
    y_true = df_genes['is_exposed'].values

    # Distance ROC on the full distance-valid set (canonical AUC = 0.908,
    # 95% CI [0.886, 0.929]; matches M1 table and manuscript text). The
    # distance predictor does not depend on baseMean, so it is NOT gated on it.
    d_valid = ~np.isnan(df_genes['nearest_methyl_distance'].values)
    dist = df_genes.loc[d_valid, 'nearest_methyl_distance'].values
    y_d = y_true[d_valid]
    fpr_d, tpr_d, thresh_d = roc_curve(y_d, -dist)
    auc_d = auc(fpr_d, tpr_d)

    # Expression ROC on its own valid subset (lower = more exposed)
    b_valid = ~np.isnan(df_genes['baseMean'].values)
    bm = df_genes.loc[b_valid, 'baseMean'].values
    y_val = y_true[b_valid]
    fpr_b, tpr_b, _ = roc_curve(y_val, -bm)
    auc_b = auc(fpr_b, tpr_b)

    ax.plot(fpr_d, tpr_d, color=COL_EXPOSED, linewidth=2.0,
            label=f'Nearest distance\nAUC = {auc_d:.3f}')
    ax.plot(fpr_b, tpr_b, color=COL_GRAY, linewidth=1.5, linestyle='--',
            label=f'Expression level\nAUC = {auc_b:.3f}')
    ax.plot([0, 1], [0, 1], 'k:', linewidth=0.6, alpha=0.5)

    # Mark 293 bp operating point on distance ROC
    # sensitivity=1.0 (all 57 exposed detected), specificity=0.802 → FPR=0.198
    op_fpr = 1.0 - 0.802
    op_tpr = 1.0
    ax.scatter([op_fpr], [op_tpr], color='#E53935', s=50, zorder=5,
               marker='D', linewidths=0)
    ax.annotate('293 bp\n(sens=1.00\nspec=0.80)',
                xy=(op_fpr, op_tpr), xytext=(op_fpr + 0.10, op_tpr - 0.18),
                fontsize=6.5, color='#E53935',
                arrowprops=dict(arrowstyle='->', color='#E53935', lw=0.8))

    # M1 (reviewer self-review): geography-stratified AUC to show the
    # core/arm distinction is NOT a Simpson's-paradox confound. region is
    # unified 0=core / 1=arm (normalized 2026-06-13).
    reg_all = df_genes['region'].astype(str).values
    dist_all = df_genes['nearest_methyl_distance'].values
    y_all = df_genes['is_exposed'].values
    strat_lines = []
    for lab, name in [('0', 'core'), ('1', 'arm')]:
        m = reg_all == lab
        if m.sum() > 0 and len(np.unique(y_all[m])) == 2:
            fpr_s, tpr_s, _ = roc_curve(y_all[m], -dist_all[m])
            a_s = auc(fpr_s, tpr_s)
            n_e = int((y_all[m] == 1).sum())
            n_s = int((y_all[m] == 0).sum())
            strat_lines.append(f'{name}: {a_s:.3f} (n={n_e}/{n_s})')
    if strat_lines:
        txt = 'Geography-stratified\n' + '\n'.join(strat_lines)
        ax.text(0.04, 0.96, txt, transform=ax.transAxes, ha='left', va='top',
                fontsize=6.0,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='#FFF8E1',
                          edgecolor='#E0C200', alpha=0.95))

    ax.set_xlabel('False positive rate', fontsize=9)
    ax.set_ylabel('True positive rate', fontsize=9)
    ax.set_title('Exposed/Shielded\nclassification (ROC)', fontsize=9, fontweight='bold')
    ax.legend(fontsize=7, loc='lower right', frameon=True, fancybox=False,
              edgecolor='#ccc', labelspacing=0.5)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_aspect('equal')


def panel_c(ax, df_quintile):
    """Panel C: Expression-independence — exposed fraction by expression quintile."""
    q = df_quintile['quintile'].values
    frac = df_quintile['frac_exposed'].values * 100

    bars = ax.bar(q, frac, color=COL_EXPOSED, alpha=0.7,
                  edgecolor='white', linewidth=0.5, width=0.6)

    # Value labels
    for qi, fi in zip(q, frac):
        ax.text(qi, fi + 0.3, f'{fi:.1f}%', ha='center', va='bottom', fontsize=7)

    # Flat trend line
    slope, intercept, r, p, _ = stats.linregress(q, frac)
    x_fit = np.array([0.5, 5.5])
    ax.plot(x_fit, slope * x_fit + intercept, 'k--', linewidth=0.8, alpha=0.5)

    # Global fraction reference
    overall = frac.mean()
    ax.axhline(overall, color=COL_EXPOSED, linewidth=0.8, linestyle=':', alpha=0.6)

    ax.text(0.97, 0.95,
            f'Jonckheere–Terpstra\np = 0.730 (NS)\n\nbaseMean AUC = 0.543',
            transform=ax.transAxes, ha='right', va='top', fontsize=7,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='#ccc', alpha=0.9))

    ax.set_xlabel('Expression quintile\n(1 = lowest, 5 = highest)', fontsize=9)
    ax.set_ylabel('Exposed genes (%)', fontsize=9)
    ax.set_ylim(0, max(frac) * 1.55)
    ax.set_xticks(q)
    ax.set_title('Expression-independent\nclassification', fontsize=9, fontweight='bold')


def main():
    apply_style()
    print('=== Figure 4: Shielded/Exposed Binary Classification ===')
    print()

    print('Loading data...')
    # Reviewer A5 fix: use latest unified 1,055-gene set with refreshed
    # is_exposed flag (57 exposed / 998 shielded) computed from
    # exposed_regulators_full_table.tsv (51_exposed_regulators_characteristics).
    # Built by extending 260312 archive (1,017) with 38 additional genes
    # whose nearest_methyl_distance was computed from HC sites.
    import pandas as pd
    df_genes = pd.read_csv(
        EPIGENOME / '52_shielded_exposed_boundary' / 'tables' /
        'all_genes_features_unified_n57.tsv',
        sep='\t'
    )
    df_genes = df_genes.dropna(subset=['nearest_methyl_distance']).copy()
    n_exp = int(df_genes['is_exposed'].sum())
    n_shi = int((df_genes['is_exposed'] == 0).sum())
    print(f'  All regulatory genes: {len(df_genes)} '
          f'({n_exp} exposed / {n_shi} shielded)')
    df_quintile = load_expression_quintile()

    # Figure layout: 1 row × 3 panels, 180 × 75 mm
    fig, axes = plt.subplots(1, 3, figsize=(mm_to_inch(180), mm_to_inch(78)))
    fig.subplots_adjust(left=0.09, right=0.97, top=0.88, bottom=0.20, wspace=0.50)

    print('Drawing Panel A: Violin plot (Shielded vs Exposed)...')
    panel_a(axes[0], df_genes)
    add_panel_label(axes[0], 'a', x=-0.22, y=1.15)

    print('Drawing Panel B: ROC curves...')
    panel_b(axes[1], df_genes)
    add_panel_label(axes[1], 'b', x=-0.22, y=1.15)

    print('Drawing Panel C: Expression quintile...')
    panel_c(axes[2], df_quintile)
    add_panel_label(axes[2], 'c', x=-0.22, y=1.15)

    # Add p-value annotation for panel A after axes limits are set
    # (re-annotate with actual y-max)

    out_path = FIG_DIR / 'Figure4_shielded_exposed'
    save_figure(fig, out_path)
    print()
    print('=== Done ===')


if __name__ == '__main__':
    main()
