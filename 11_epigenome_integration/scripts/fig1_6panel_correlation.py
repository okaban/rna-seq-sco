#!/usr/bin/env python3
"""
Fig. 3B: 6-panel scatter plot of methylation change vs expression change
6mA/4mC × T2vsT1/T3vsT1/T3vsT2

Color-coded by methylation site dynamics (Gained/Stable/Lost) with
categorical statistical tests (Kruskal-Wallis, Mann-Whitney U).

Streptomyces coelicolor A3(2) M145
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats
from statsmodels.stats.multitest import multipletests
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# === Paths ===
BASE = Path('/Users/okaban/bioinfo/rna-seq/11_epigenome_integration')
DATA_FILE = BASE / 'analysis/01_integration/integrated_methyl_expression_weighted.csv'
OUTPUT_DIR = BASE / 'analysis/02_publication_figures'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# === Style ===
COLORS = {
    'gained': '#E53935',   # Red — site gained
    'lost': '#1565C0',     # Blue — site lost
    'stable': '#B0BEC5',   # Gray — site present in both
}

plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 1.0,
    'xtick.major.width': 1.0,
    'ytick.major.width': 1.0,
})

# === Configuration ===
MOD_TYPES = ['6mA', '4mC']
# (comp_key, label, tp_before, tp_after)
COMPARISONS = [
    ('T2_vs_T1', 'T2 vs T1', 'T1', 'T2'),
    ('T3_vs_T1', 'T3 vs T1', 'T1', 'T3'),
    ('T3_vs_T2', 'T3 vs T2', 'T2', 'T3'),
]


def load_data():
    """Load integrated methylation-expression data."""
    df = pd.read_csv(DATA_FILE)
    print(f"Loaded {len(df)} genes from {DATA_FILE.name}")
    return df


def classify_methyl_dynamics(row, mod, tp_before, tp_after):
    """Classify methylation change as Gained/Lost/Stable based on site presence."""
    cnt_before = row[f'{mod}_{tp_before}_count']
    cnt_after = row[f'{mod}_{tp_after}_count']
    if cnt_before == 0 and cnt_after > 0:
        return 'Gained'
    elif cnt_before > 0 and cnt_after == 0:
        return 'Lost'
    elif cnt_before > 0 and cnt_after > 0:
        return 'Stable'
    return None  # both zero — should not appear after filtering


def compute_categorical_stats(df):
    """Compute Kruskal-Wallis and Mann-Whitney U for all 6 conditions, apply BH FDR."""
    results = []
    for mod in MOD_TYPES:
        for comp_key, comp_label, tp_before, tp_after in COMPARISONS:
            change_col = f'{mod}_change_{comp_key}'
            fc_col = f'log2FC_{comp_key}'

            valid = df.dropna(subset=[fc_col])
            valid = valid[valid[change_col] != 0].copy()

            if len(valid) < 10:
                results.append({
                    'mod': mod, 'comp_key': comp_key, 'comp_label': comp_label,
                    'kw_p': 1.0, 'mw_p': 1.0, 'mw_r': 0.0, 'n': len(valid),
                    'n_gained': 0, 'n_lost': 0, 'n_stable': 0,
                    'med_gained': np.nan, 'med_lost': np.nan, 'med_stable': np.nan,
                    'spearman_stable_r': np.nan, 'spearman_stable_p': 1.0,
                    'spearman_stable_n': 0,
                })
                continue

            # Classify
            valid['dyn'] = valid.apply(
                lambda x: classify_methyl_dynamics(x, mod, tp_before, tp_after), axis=1
            )
            valid = valid.dropna(subset=['dyn'])

            gained = valid[valid['dyn'] == 'Gained'][fc_col]
            lost = valid[valid['dyn'] == 'Lost'][fc_col]
            stable = valid[valid['dyn'] == 'Stable'][fc_col]

            # Kruskal-Wallis (3-group)
            groups = [g for g in [gained, lost, stable] if len(g) >= 3]
            if len(groups) >= 2:
                kw_stat, kw_p = stats.kruskal(*groups)
            else:
                kw_stat, kw_p = 0.0, 1.0

            # Mann-Whitney U (Gained vs Lost)
            if len(gained) >= 3 and len(lost) >= 3:
                u_stat, mw_p = stats.mannwhitneyu(gained, lost, alternative='two-sided')
                # Rank-biserial correlation as effect size: r = 1 - 2U/(n1*n2)
                mw_r = 1 - 2 * u_stat / (len(gained) * len(lost))
            else:
                u_stat, mw_p, mw_r = 0.0, 1.0, 0.0

            # Concordance rate: proportion of genes where methylation and expression
            # move in the same direction (Gained+Up or Lost+Down)
            gained_up = ((valid['dyn'] == 'Gained') & (valid[fc_col] > 0)).sum()
            lost_down = ((valid['dyn'] == 'Lost') & (valid[fc_col] < 0)).sum()
            n_directional = len(gained) + len(lost)
            concordance_rate = (gained_up + lost_down) / n_directional if n_directional > 0 else 0.0

            # Spearman within Stable group (dose-response test)
            if len(stable) >= 10:
                stable_data = valid[valid['dyn'] == 'Stable']
                sp_r, sp_p = stats.spearmanr(stable_data[change_col], stable_data[fc_col])
            else:
                sp_r, sp_p = np.nan, 1.0

            results.append({
                'mod': mod, 'comp_key': comp_key, 'comp_label': comp_label,
                'kw_stat': kw_stat, 'kw_p': kw_p,
                'mw_u': u_stat, 'mw_p': mw_p, 'mw_r': mw_r,
                'concordance_rate': concordance_rate,
                'n_gained_up': gained_up, 'n_lost_down': lost_down,
                'n': len(valid),
                'n_gained': len(gained), 'n_lost': len(lost), 'n_stable': len(stable),
                'med_gained': gained.median() if len(gained) > 0 else np.nan,
                'med_lost': lost.median() if len(lost) > 0 else np.nan,
                'med_stable': stable.median() if len(stable) > 0 else np.nan,
                'spearman_stable_r': sp_r, 'spearman_stable_p': sp_p,
                'spearman_stable_n': len(stable),
            })

    # BH FDR correction on Mann-Whitney p-values
    p_vals = [r['mw_p'] for r in results]
    _, fdr_vals, _, _ = multipletests(p_vals, method='fdr_bh')
    for i, r in enumerate(results):
        r['mw_fdr'] = fdr_vals[i]
        r['sig'] = fdr_vals[i] < 0.05

    # BH FDR on Kruskal-Wallis
    kw_pvals = [r['kw_p'] for r in results]
    _, kw_fdr_vals, _, _ = multipletests(kw_pvals, method='fdr_bh')
    for i, r in enumerate(results):
        r['kw_fdr'] = kw_fdr_vals[i]

    return results


def plot_6panel(df, stat_results):
    """Generate 2×3 scatter plot grid with categorical coloring."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))

    # Build lookup
    stat_lookup = {}
    for r in stat_results:
        stat_lookup[(r['mod'], r['comp_key'])] = r

    for row_idx, mod in enumerate(MOD_TYPES):
        for col_idx, (comp_key, comp_label, tp_before, tp_after) in enumerate(COMPARISONS):
            ax = axes[row_idx, col_idx]
            sr = stat_lookup[(mod, comp_key)]

            change_col = f'{mod}_change_{comp_key}'
            fc_col = f'log2FC_{comp_key}'

            # Filter
            plot_df = df.dropna(subset=[fc_col]).copy()
            plot_df = plot_df[plot_df[change_col] != 0].copy()

            # Classify dynamics
            plot_df['dyn'] = plot_df.apply(
                lambda x: classify_methyl_dynamics(x, mod, tp_before, tp_after), axis=1
            )
            plot_df = plot_df.dropna(subset=['dyn'])

            # Color assignment by methylation dynamics
            color_map = {
                'Gained': COLORS['gained'],
                'Lost': COLORS['lost'],
                'Stable': COLORS['stable'],
            }
            colors = plot_df['dyn'].map(color_map)

            # Highlight significant panels
            if sr['sig']:
                ax.set_facecolor('#FFFDE7')
                for spine in ax.spines.values():
                    spine.set_linewidth(2.5)
                    spine.set_edgecolor('#D32F2F')
            else:
                for spine in ['top', 'right']:
                    ax.spines[spine].set_visible(False)

            # Scatter plot — draw Stable first (background), then Lost & Gained on top
            for dyn_cat, zorder in [('Stable', 2), ('Lost', 3), ('Gained', 3)]:
                mask = plot_df['dyn'] == dyn_cat
                if mask.sum() == 0:
                    continue
                ax.scatter(
                    plot_df.loc[mask, change_col], plot_df.loc[mask, fc_col],
                    c=color_map[dyn_cat], alpha=0.6, s=30,
                    edgecolors='white', linewidth=0.3, zorder=zorder,
                )

            # Add regression line for overall correlation visualization
            x_data = plot_df[change_col].values
            y_data = plot_df[fc_col].values
            if len(x_data) > 10:
                # Linear regression
                slope, intercept, r_value, p_value, std_err = stats.linregress(x_data, y_data)
                x_line = np.linspace(x_data.min(), x_data.max(), 100)
                y_line = slope * x_line + intercept

                # Line color based on slope direction
                if slope > 0:
                    line_color = '#2E7D32'  # Green for positive
                else:
                    line_color = '#C62828'  # Red for negative

                # Line style based on significance
                line_alpha = 0.9 if sr['sig'] else 0.5
                line_width = 2.5 if sr['sig'] else 1.5

                ax.plot(x_line, y_line, color=line_color, linewidth=line_width,
                       linestyle='-', alpha=line_alpha, zorder=4,
                       label=f'r={r_value:.2f}')

                # Add Spearman correlation coefficient annotation near the line
                # Position at right side of plot
                x_annot = x_data.max() * 0.7
                y_annot = slope * x_annot + intercept
                sp_r, sp_p = stats.spearmanr(x_data, y_data)
                ax.annotate(f'Spearman r = {sp_r:.2f}', xy=(x_annot, y_annot),
                           fontsize=8, fontweight='bold', color=line_color,
                           ha='left', va='bottom' if slope > 0 else 'top',
                           bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                                    alpha=0.9, edgecolor=line_color, linewidth=1))

            # Reference lines
            ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.5, alpha=0.4)
            ax.axvline(x=0, color='gray', linestyle='-', linewidth=0.5, alpha=0.4)

            # Gap shading to make trimodal structure explicit
            y_lo, y_hi = ax.get_ylim()
            ax.axvspan(20, 45, alpha=0.04, color='gray', zorder=0)
            ax.axvspan(-45, -20, alpha=0.04, color='gray', zorder=0)

            # Statistics annotation
            # Concordance rate: proportion where methylation & expression move same direction
            conc_pct = sr['concordance_rate'] * 100
            if sr['sig']:
                stats_text = (
                    f"MW (G vs L): p = {sr['mw_p']:.1e}\n"
                    f"FDR = {sr['mw_fdr']:.3f} *\n"
                    f"Concordance = {conc_pct:.1f}%\n"
                    f"KW (3-grp): p = {sr['kw_p']:.1e}"
                )
                bbox_color = '#FFCDD2'
            else:
                stats_text = (
                    f"MW (G vs L): p = {sr['mw_p']:.2e}\n"
                    f"FDR = {sr['mw_fdr']:.2f} (n.s.)\n"
                    f"Concordance = {conc_pct:.1f}%\n"
                    f"KW (3-grp): p = {sr['kw_p']:.2e}"
                )
                bbox_color = 'white'

            ax.text(0.05, 0.95, stats_text, transform=ax.transAxes,
                    fontsize=8, verticalalignment='top', family='monospace',
                    bbox=dict(boxstyle='round', facecolor=bbox_color,
                              alpha=0.85, edgecolor='gray'))

            # Group n and median annotations (bottom of panel)
            grp_text = (
                f"Gained: n={sr['n_gained']}, med={sr['med_gained']:+.2f}\n"
                f"Stable: n={sr['n_stable']}, med={sr['med_stable']:+.2f}\n"
                f"Lost: n={sr['n_lost']}, med={sr['med_lost']:+.2f}"
            )
            ax.text(0.95, 0.05, grp_text, transform=ax.transAxes,
                    fontsize=7, verticalalignment='bottom', ha='right',
                    family='monospace',
                    bbox=dict(boxstyle='round', facecolor='white',
                              alpha=0.8, edgecolor='#BDBDBD'))

            # Axis labels
            if row_idx == 1:
                ax.set_xlabel('Δ Methylation frequency (%)', fontweight='bold')
            if col_idx == 0:
                ax.set_ylabel('log₂ Fold Change', fontweight='bold')

            # Panel title
            ax.set_title(f'{mod}  —  {comp_label}', fontweight='bold', fontsize=11)

            # Symmetric x-axis
            x_max = max(abs(plot_df[change_col].min()), abs(plot_df[change_col].max())) * 1.1
            if x_max > 0:
                ax.set_xlim(-x_max, x_max)

    # Shared legend at bottom
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS['gained'],
               markersize=8, label='Gained (site absent → present)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS['stable'],
               markersize=8, label='Stable (site present in both)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS['lost'],
               markersize=8, label='Lost (site present → absent)'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=3,
               bbox_to_anchor=(0.5, -0.01), frameon=False, fontsize=10)

    plt.tight_layout(rect=[0, 0.03, 1, 1])

    # Save
    out_base = OUTPUT_DIR / 'Fig1_methylation_expression_correlation_6panel'
    fig.savefig(f'{out_base}.pdf')
    fig.savefig(f'{out_base}.svg')
    fig.savefig(f'{out_base}.png', dpi=300)
    plt.close()
    print(f"\nSaved: {out_base}.pdf / .svg / .png")


def main():
    df = load_data()
    stat_results = compute_categorical_stats(df)

    print("\n=== Categorical Statistics (BH FDR on Mann-Whitney) ===")
    print(f"{'Condition':<20} {'MW p':>12} {'MW FDR':>10} {'MW r':>8} "
          f"{'KW p':>12} {'nG':>5} {'nS':>5} {'nL':>5} {'Sig':>5}")
    print("-" * 90)
    for r in stat_results:
        sig_mark = '***' if r['sig'] else ''
        print(f"{r['mod']} {r['comp_label']:<14} {r['mw_p']:>12.2e} {r['mw_fdr']:>10.4f} "
              f"{r['mw_r']:>8.3f} {r['kw_p']:>12.2e} "
              f"{r['n_gained']:>5} {r['n_stable']:>5} {r['n_lost']:>5} {sig_mark:>5}")

    print("\n=== Within-Stable Spearman (dose-response test) ===")
    print(f"{'Condition':<20} {'r':>8} {'p':>12} {'n':>6}")
    print("-" * 50)
    for r in stat_results:
        print(f"{r['mod']} {r['comp_label']:<14} {r['spearman_stable_r']:>8.4f} "
              f"{r['spearman_stable_p']:>12.2e} {r['spearman_stable_n']:>6}")

    print("\n=== Group Median log2FC ===")
    print(f"{'Condition':<20} {'Gained':>10} {'Stable':>10} {'Lost':>10}")
    print("-" * 55)
    for r in stat_results:
        print(f"{r['mod']} {r['comp_label']:<14} {r['med_gained']:>+10.3f} "
              f"{r['med_stable']:>+10.3f} {r['med_lost']:>+10.3f}")

    plot_6panel(df, stat_results)

    # Save statistics CSV
    stats_df = pd.DataFrame(stat_results)
    stats_csv = OUTPUT_DIR / 'methylation_expression_categorical_stats.csv'
    stats_df.to_csv(stats_csv, index=False)
    print(f"\nSaved: {stats_csv}")


if __name__ == '__main__':
    main()
