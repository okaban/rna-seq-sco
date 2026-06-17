#!/usr/bin/env python3
"""
[RETRACTED / DO NOT USE — 2026-06-17]
This switch-model figure rests on the expression-selected "two antagonistic blocs"
claim and panel D reports the tautological distance-classifier AUC (0.917) on the
old 57/998 set — all circular under the non-circular reframe (62/989, no AUC).
The manuscript Figure 6 is now produced by `30_figure6_spatial_organizer.py`
(permissive spatial-organiser schematic). Do NOT regenerate this figure or copy its
output into Writing/fig_images/Figure6.png. Kept only for git history.

--- original docstring ---
New Figure 5: Vegetative-to-Developmental Switch Model

(A) Repression bloc — vegetative programs OFF (key genes + LFC waterfall)
(B) Activation bloc — developmental programs ON (key genes + LFC waterfall)
(C) Simultaneous switch timing (aggregate expression by bloc)
(D) Model schematic (Gatekeeper architecture summary)
"""

import importlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch
from matplotlib.lines import Line2D

_utils = importlib.import_module('00_shared_utils')
for _attr in dir(_utils):
    if not _attr.startswith('_'):
        globals()[_attr] = getattr(_utils, _attr)


# TF family → display color for bloc panels
FAMILY_COLORS = {
    'TetR':                '#E65100',   # dark orange
    'Sigma factor':        '#1565C0',   # dark blue
    'Sensor kinase':       '#6A1B9A',   # deep purple
    'Response regulator':  '#AD1457',   # pink-purple
    'WhiB':                '#2E7D32',   # dark green
    'MerR':                '#00838F',   # teal
    'LysR':                '#F9A825',   # amber
    'MarR':                '#4E342E',   # brown
    'GntR':                '#546E7A',   # blue-gray
    'ArsR':                '#795548',   # medium brown
    'IclR':                '#558B2F',   # olive
    'LacI':                '#0277BD',   # steel blue
    'DeoR':                '#6D4C41',   # warm brown
    'ROK':                 '#37474F',   # dark slate
    'HTH (other)':         '#90A4AE',   # light blue-gray
    'Other regulatory':    '#B0BEC5',   # light gray
}

# Short functional annotation for TF families
FAMILY_SHORT = {
    'TetR':                'TetR',
    'Sigma factor':        'Sigma',
    'Sensor kinase':       'SK',
    'Response regulator':  'RR',
    'WhiB':                'WhiB',
    'MerR':                'MerR',
    'LysR':                'LysR',
    'MarR':                'MarR',
    'GntR':                'GntR',
    'ArsR':                'ArsR',
    'IclR':                'IclR',
    'LacI':                'LacI',
    'DeoR':                'DeoR',
    'ROK':                 'ROK',
    'HTH (other)':         'HTH',
    'Other regulatory':    'Other',
}

# Genes with known names worth highlighting
HIGHLIGHT_GENES = {'ramR', 'nsdB', 'fasR', 'sigJ', 'fxsT', 'tcrA', 'mfd'}


def _get_label(row):
    gene = row.get('gene_name', '')
    if pd.isna(gene) or gene == '':
        return row['old_locus_tag'] if pd.notna(row.get('old_locus_tag')) else row['locus_tag']
    return gene


def _family_color(fam):
    return FAMILY_COLORS.get(str(fam), '#B0BEC5')


def panel_a(ax, df_temporal):
    """Panel A: Repression bloc — vegetative programs OFF.
    Shows grouped bars: T2vsT1 (lighter, narrow) + T3vsT1 (darker, wide).
    Bars colored by TF family.
    """
    rep = df_temporal[df_temporal['bloc'] == 'repression'].copy()
    rep = rep.sort_values('LFC_T3vsT1', ascending=True).reset_index(drop=True)

    n = len(rep)
    y_pos = np.arange(n)

    # Draw grouped bars: T3vsT1 (back, wide) and T2vsT1 (front, narrow)
    for i, row in rep.iterrows():
        fc = _family_color(row['tf_family'])
        # T3vsT1 — wider, semi-transparent
        ax.barh(y_pos[i] + 0.18, row['LFC_T3vsT1'], color=fc, alpha=0.55,
                edgecolor='none', height=0.38)
        # T2vsT1 — narrower, solid
        ax.barh(y_pos[i] - 0.18, row['LFC_T2vsT1'], color=fc, alpha=0.90,
                edgecolor='none', height=0.28)

    # Gene labels — bold if named gene
    labels = [_get_label(row) for _, row in rep.iterrows()]
    ax.set_yticks(y_pos)
    yticklabels = ax.set_yticklabels(labels, fontsize=5.5)
    for i, (tick, (_, row)) in enumerate(zip(ax.get_yticklabels(), rep.iterrows())):
        if _get_label(row) in HIGHLIGHT_GENES:
            tick.set_fontweight('bold')

    # Family annotation on right side
    xmax = ax.get_xlim()[1]
    for i, row in rep.iterrows():
        fam_short = FAMILY_SHORT.get(str(row['tf_family']), '')
        if fam_short not in ('HTH', 'Other', ''):
            ax.text(xmax * 1.02, y_pos[i], fam_short, fontsize=4.5,
                    va='center', ha='left', color=_family_color(row['tf_family']))

    ax.set_xlabel('$\\log_2$FC', fontsize=9)
    ax.axvline(0, color='gray', linewidth=0.5, linestyle=':')
    ax.set_title('Repression bloc\n(vegetative programs OFF)', fontsize=9,
                 fontweight='bold', color=COL_REPRESSION)

    # Mini legend for bars
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#555', alpha=0.9, label='T2 vs T1'),
        Patch(facecolor='#555', alpha=0.5, label='T3 vs T1'),
    ]
    ax.legend(handles=legend_elements, fontsize=6, loc='lower left',
              frameon=True, framealpha=0.9, edgecolor='#ccc')


def panel_b(ax, df_temporal):
    """Panel B: Activation bloc — developmental programs ON.
    Shows grouped bars: T2vsT1 (lighter) + T3vsT1 (darker), colored by TF family.
    """
    act = df_temporal[df_temporal['bloc'] == 'activation'].copy()
    act = act.sort_values('LFC_T3vsT1', ascending=True).reset_index(drop=True)

    n = len(act)
    y_pos = np.arange(n)

    for i, row in act.iterrows():
        fc = _family_color(row['tf_family'])
        ax.barh(y_pos[i] + 0.18, row['LFC_T3vsT1'], color=fc, alpha=0.55,
                edgecolor='none', height=0.38)
        ax.barh(y_pos[i] - 0.18, row['LFC_T2vsT1'], color=fc, alpha=0.90,
                edgecolor='none', height=0.28)

    labels = [_get_label(row) for _, row in act.iterrows()]
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=5.5)
    for tick, (_, row) in zip(ax.get_yticklabels(), act.iterrows()):
        if _get_label(row) in HIGHLIGHT_GENES:
            tick.set_fontweight('bold')

    xmax = ax.get_xlim()[1]
    for i, row in act.iterrows():
        fam_short = FAMILY_SHORT.get(str(row['tf_family']), '')
        if fam_short not in ('HTH', 'Other', ''):
            ax.text(xmax * 1.02, y_pos[i], fam_short, fontsize=4.5,
                    va='center', ha='left', color=_family_color(row['tf_family']))

    ax.set_xlabel('$\\log_2$FC', fontsize=9)
    ax.axvline(0, color='gray', linewidth=0.5, linestyle=':')
    ax.set_title('Activation bloc\n(developmental programs ON)', fontsize=9,
                 fontweight='bold', color=COL_ACTIVATION)

    # TF family legend
    from matplotlib.patches import Patch
    # Show only families present in this bloc
    families_present = act['tf_family'].unique()
    legend_els = [Patch(facecolor=_family_color(f), alpha=0.85,
                        label=FAMILY_SHORT.get(str(f), str(f)))
                  for f in sorted(families_present)]
    ax.legend(handles=legend_els, fontsize=5.5, loc='lower right',
              frameon=True, framealpha=0.9, edgecolor='#ccc',
              ncol=2, title='TF family', title_fontsize=6)


def panel_c(ax, df_temporal, df_bloc):
    """Panel C: Simultaneous switch timing — aggregate by bloc."""
    timepoints = ['T1', 'T2', 'T3']
    x = [1, 2, 3]

    for bloc, color, marker, ls in [
        ('activation', COL_ACTIVATION, 'o', '-'),
        ('repression', COL_REPRESSION, 's', '-'),
    ]:
        sub = df_temporal[df_temporal['bloc'] == bloc]
        z_cols = ['T1_z', 'T2_z', 'T3_z']

        means = [sub[c].mean() for c in z_cols]
        sems = [sub[c].std() / np.sqrt(len(sub)) for c in z_cols]

        ax.errorbar(x, means, yerr=sems, color=color, linewidth=2.5,
                    marker=marker, markersize=8, capsize=3, capthick=1.5,
                    label=f'{bloc.capitalize()} (n={len(sub)})')
        ax.fill_between(x,
                        [m - s for m, s in zip(means, sems)],
                        [m + s for m, s in zip(means, sems)],
                        color=color, alpha=0.15)

    ax.axhline(0, color='gray', linewidth=0.5, linestyle=':')

    # Mark the switch point
    ax.axvspan(1.5, 2.5, color='#F5F5F5', alpha=0.5, zorder=0)
    ax.text(2, ax.get_ylim()[0] + 0.1, 'Switch\npoint', ha='center',
            fontsize=7, fontstyle='italic', color=COL_DARK)

    # Phase ratio annotation
    ax.text(0.97, 0.50, 'Both blocs engage\nat T1→T2 transition\n(p = 0.459, simultaneous)',
            transform=ax.transAxes, ha='right', va='center', fontsize=7,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='#ccc', alpha=0.9))

    ax.set_xticks(x)
    ax.set_xticklabels(timepoints)
    ax.set_xlabel('Timepoint', fontsize=9)
    ax.set_ylabel('Mean expression (z-score)', fontsize=9)
    ax.set_title('Simultaneous switch timing', fontsize=9, fontweight='bold')
    ax.legend(fontsize=7, loc='upper left', frameon=False)


def panel_d(ax):
    """Panel D: Observation summary — what this study established and what remains open."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    # ── Observation 1: Methylation redistribution ──────────────────────────
    box1 = FancyBboxPatch((0.3, 8.0), 9.4, 1.65, boxstyle='round,pad=0.15',
                           facecolor=COL_4mC, alpha=0.10,
                           edgecolor=COL_4mC, linewidth=1.5)
    ax.add_patch(box1)
    ax.text(0.7, 9.28, 'Observation 1', ha='left', va='center',
            fontsize=7, fontweight='bold', color=COL_4mC)
    ax.text(0.7, 8.82, 'GCCGGC 4mC undergoes geographic redistribution during development',
            ha='left', va='center', fontsize=6.5, color=COL_DARK)
    ax.text(0.7, 8.35, 'T1: 83% core  →  T2: 82% arm  →  T3: residual (n=21)',
            ha='left', va='center', fontsize=6, color=COL_DARK)

    # ── Observation 2: 57 exposed regulators (structural, minimal Layer 2) ──
    box2 = FancyBboxPatch((0.3, 5.95), 9.4, 1.75, boxstyle='round,pad=0.15',
                           facecolor=COL_EXPOSED, alpha=0.07,
                           edgecolor=COL_EXPOSED, linewidth=1.2)
    ax.add_patch(box2)
    ax.text(0.7, 7.33, 'Observation 2', ha='left', va='center',
            fontsize=7, fontweight='bold', color=COL_EXPOSED)
    ax.text(0.7, 6.88, '57 of 1,055 regulatory genes lack methylation protection at TSS',
            ha='left', va='center', fontsize=6.5, color=COL_DARK)
    ax.text(0.7, 6.38, '(293 bp threshold, AUC = 0.917; 998 shielded, 57 exposed)',
            ha='left', va='center', fontsize=6, color=COL_DARK)

    # ── Observation 3: Co-regulation (with causation caveat) ───────────────
    box3 = FancyBboxPatch((0.3, 3.35), 9.4, 2.35, boxstyle='round,pad=0.15',
                           facecolor='#F5F5F5', alpha=0.6,
                           edgecolor=COL_DARK, linewidth=1.2)
    ax.add_patch(box3)
    ax.text(0.7, 5.33, 'Observation 3', ha='left', va='center',
            fontsize=7, fontweight='bold', color=COL_DARK)
    ax.text(0.7, 4.87, 'The 57 exposed regulators form two antagonistic blocs',
            ha='left', va='center', fontsize=6.5, color=COL_DARK)

    # Activation/Repression sub-labels inline
    ax.text(0.7, 4.40,
            'Activation (n=35): TCS, Sigma, WhiB — development ON',
            ha='left', va='center', fontsize=6, color=COL_ACTIVATION)
    ax.text(0.7, 3.95,
            'Repression (n=26): TetR, metabolic regulators — vegetative OFF',
            ha='left', va='center', fontsize=6, color=COL_REPRESSION)
    ax.text(0.7, 3.52,
            'Both blocs switch simultaneously at T1→T2 (MW p = 0.459)',
            ha='left', va='center', fontsize=6, color=COL_DARK)

    # ── Causation caveat box ───────────────────────────────────────────────
    box_cav = FancyBboxPatch((0.3, 0.5), 9.4, 2.55, boxstyle='round,pad=0.15',
                              facecolor='#FFF9C4', alpha=0.85,
                              edgecolor='#F9A825', linewidth=1.5)
    ax.add_patch(box_cav)
    ax.text(5.0, 2.72, 'Causation unresolved', ha='center', va='center',
            fontsize=7.5, fontweight='bold', color='#E65100')
    ax.text(5.0, 2.25,
            'Methylation timing is decoupled from expression timing',
            ha='center', va='center', fontsize=6.5, color=COL_DARK)
    ax.text(5.0, 1.78,
            r'($\rho$ = 0.136, p = 0.30 at T1$\to$T2; $\rho$ = 0.012, p = 0.93 at T2$\to$T3)',
            ha='center', va='center', fontsize=6, color=COL_DARK)
    ax.text(5.0, 1.28,
            'This study characterizes methylation distribution and motif patterns;',
            ha='center', va='center', fontsize=6, color=COL_DARK)
    ax.text(5.0, 0.82,
            'a causal role for methylation in gene expression regulation remains to be established.',
            ha='center', va='center', fontsize=6, color=COL_DARK, fontstyle='italic')

    ax.set_title('Summary of Key Observations', fontsize=10, fontweight='bold')


def main():
    apply_style()
    print('=== New Figure 5: Vegetative-to-Developmental Switch Model ===')
    print()

    # Load data
    print('Loading data...')
    df_temporal = load_temporal_classification()
    df_bloc = load_bloc_comparison()

    # Create figure: 180mm × 200mm
    fig = plt.figure(figsize=(mm_to_inch(180), mm_to_inch(210)))

    gs = fig.add_gridspec(2, 2, hspace=0.45, wspace=0.35,
                          left=0.10, right=0.97, top=0.95, bottom=0.04,
                          height_ratios=[1.3, 1])

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    print('Drawing Panel A: Repression bloc...')
    panel_a(ax_a, df_temporal)
    add_panel_label(ax_a, 'a', x=-0.15, y=1.05)

    print('Drawing Panel B: Activation bloc...')
    panel_b(ax_b, df_temporal)
    add_panel_label(ax_b, 'b', x=-0.15, y=1.05)

    print('Drawing Panel C: Simultaneous timing...')
    panel_c(ax_c, df_temporal, df_bloc)
    add_panel_label(ax_c, 'c', x=-0.15, y=1.10)

    print('Drawing Panel D: Model schematic...')
    panel_d(ax_d)
    add_panel_label(ax_d, 'd', x=-0.08, y=1.10)

    # Save
    out_path = FIG_DIR / 'new_Figure4_switch_model'
    save_figure(fig, out_path)
    print()
    print('=== Done ===')


if __name__ == '__main__':
    main()
