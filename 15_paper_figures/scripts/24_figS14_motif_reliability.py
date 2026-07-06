"""
FigS14: Motif reliability and occupancy analysis (A-1, A-2, A-3)

Panels:
  a. AAGCCCG methylated position in motif (bar chart)
  b. GCCGGC methylated positions (palindromic)
  c. Motif occupancy across timepoints (AAGCCCG / GCCGGC / GATC)
  d. 6mA motif assignment breakdown (stacked bar: AAGCCCG / GATC / Unassigned)
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from shared_utils_local import apply_style, add_panel_label, save_figure, FIG_SUP_DIR
from importlib import import_module as _im
mm_to_inch = _im('00_shared_utils').mm_to_inch

# === Data paths ===
RELIABILITY_DIR = Path('/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/61_motif_reliability')
TABLES = RELIABILITY_DIR / 'tables'

apply_style()

# ── Load data ────────────────────────────────────────────────────────────────
pos_aag = pd.read_csv(TABLES / 'A1_AAGCCCG_methylated_position.tsv', sep='\t')
pos_gcc = pd.read_csv(TABLES / 'A1_GCCGGC_methylated_position.tsv', sep='\t')
occupancy = pd.read_csv(TABLES / 'A2_motif_occupancy.tsv', sep='\t')
assignment = pd.read_csv(TABLES / 'A3_6mA_motif_assignment_by_timepoint.tsv', sep='\t')

# ── Panel a: AAGCCCG position ────────────────────────────────────────────────
aag_counts = pos_aag.groupby('pos_in_motif').size().reset_index(name='count')
aag_total = aag_counts['count'].sum()
aag_counts['pct'] = aag_counts['count'] / aag_total * 100
# Position labels: AAGCCCG = A(0)A(1)G(2)C(3)C(4)C(5)G(6)
# Only A positions (0,1) expected as 6mA
pos_labels_aag = {0: 'A (pos 0)', 1: 'A (pos 1)', 2: 'G (pos 2)',
                  3: 'C (pos 3)', 4: 'C (pos 4)', 5: 'C (pos 5)', 6: 'G (pos 6)'}

# ── Panel b: GCCGGC position ─────────────────────────────────────────────────
gcc_counts = pos_gcc.groupby('pos_in_motif').size().reset_index(name='count')
gcc_total = gcc_counts['count'].sum()
gcc_counts['pct'] = gcc_counts['count'] / gcc_total * 100
# GCCGGC = G(0)C(1)C(2)G(3)G(4)C(5) — 4mC expected at C positions
pos_labels_gcc = {0:'G(0)', 1:'C(1)', 2:'C(2)', 3:'G(3)', 4:'G(4)', 5:'C(5)'}

# ── Figure layout ─────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(mm_to_inch(174), mm_to_inch(150)))  # NAR full-width
axes = axes.flatten()

# --- Panel a ---
ax = axes[0]
x = aag_counts['pos_in_motif'].values
y = aag_counts['pct'].values
bar_colors = ['#7E57C2' if p in [0, 1] else '#B0BEC5' for p in x]
bars = ax.bar(x, y, color=bar_colors, edgecolor='white', linewidth=0.5)
for bar, pct in zip(bars, y):
    if pct > 2:
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                f'{pct:.0f}%', ha='center', va='bottom', fontsize=6)
ax.set_xticks(x)
pos_aag_label = {0: 'A (pos 0)\n1st A', 1: 'A (pos 1)\n2nd A', 2: 'G (pos 2)'}
ax.set_xticklabels([pos_aag_label.get(p, str(p)) for p in x], fontsize=6)
ax.set_ylabel('% of methylated sites', fontsize=8)
ax.set_title('AAGCCCG: methylated position in motif', fontsize=8)
ax.set_ylim(0, max(y) * 1.2)
ax.set_xlabel('Position in AAGCCCG', fontsize=7)
# Motif label
ax.text(0.5, 0.95, 'AAGCCCG', transform=ax.transAxes,
        ha='center', va='top', fontsize=8, fontfamily='monospace',
        color='#7E57C2', fontweight='bold')
add_panel_label(ax, 'a')

# --- Panel b ---
ax = axes[1]
x = gcc_counts['pos_in_motif'].values
y = gcc_counts['pct'].values
# C positions (1, 2, 5) are the 4mC sites in GCCGGC
bar_colors = ['#E53935' if p in [1, 2, 5] else '#B0BEC5' for p in x]
bars = ax.bar(x, y, color=bar_colors, edgecolor='white', linewidth=0.5)
for bar, pct in zip(bars, y):
    if pct > 2:
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                f'{pct:.0f}%', ha='center', va='bottom', fontsize=6)
ax.set_xticks(x)
pos_gcc_label = {0: 'G(0)', 1: 'C(1)', 2: 'C(2)', 3: 'G(3)\n[C RC]', 4: 'G(4)', 5: 'C(5)'}
ax.set_xticklabels([pos_gcc_label.get(p, str(p)) for p in x], fontsize=6)
ax.set_ylabel('% of methylated sites', fontsize=8)
ax.set_title('GCCGGC: methylated position in motif', fontsize=8)
ax.set_ylim(0, max(y) * 1.2)
ax.set_xlabel('Position in GCCGGC', fontsize=7)
ax.text(0.5, 0.95, 'GCCGGC (palindrome)', transform=ax.transAxes,
        ha='center', va='top', fontsize=8, fontfamily='monospace',
        color='#E53935', fontweight='bold')
add_panel_label(ax, 'b')

# --- Panel c: Motif occupancy ---
ax = axes[2]
tp_order = ['T1', 'T2', 'T3']
motifs = ['AAGCCCG', 'GCCGGC', 'GATC']
motif_colors = {'AAGCCCG': '#7E57C2', 'GCCGGC': '#E53935', 'GATC': '#1565C0'}
tp_labels = ['T1 (12 h)', 'T2 (24 h)', 'T3 (50 h)']
x = np.arange(len(tp_order))
width = 0.25

for i, motif in enumerate(motifs):
    sub = occupancy[occupancy['motif'] == motif].set_index('timepoint')
    vals = [sub.loc[tp, 'occupancy_pct'] if tp in sub.index else 0 for tp in tp_order]
    bars = ax.bar(x + i * width, vals, width, label=motif,
                  color=motif_colors[motif], alpha=0.85, edgecolor='white')
    for bar, v in zip(bars, vals):
        if v > 0.5:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.15,
                    f'{v:.1f}%', ha='center', va='bottom', fontsize=7.5)

ax.set_xticks(x + width)
ax.set_xticklabels(tp_labels, fontsize=6)
ax.set_ylabel('Motif occupancy (%)', fontsize=8)
ax.set_title('R-M recognition site occupancy\n(methylated / total genomic occurrences)', fontsize=7)
ax.legend(fontsize=6)
add_panel_label(ax, 'c')

# --- Panel d: 6mA motif assignment ---
ax = axes[3]
tp_order2 = ['T1', 'T2', 'T3']
tp_labels2 = ['T1 (12 h)', 'T2 (24 h)', 'T3 (50 h)']
asg = assignment.set_index('timepoint')

cats = ['AAGCCCG_only', 'GATC_only', 'both', 'unassigned']
cat_labels = ['AAGCCCG', 'GATC', 'Both', 'Unassigned']
cat_colors = ['#7E57C2', '#1565C0', '#43A047', '#B0BEC5']

bottoms = np.zeros(len(tp_order2))
for cat, label, color in zip(cats, cat_labels, cat_colors):
    vals = []
    for tp in tp_order2:
        row = asg.loc[tp]
        pct = row[cat] / row['total'] * 100
        vals.append(pct)
    bars = ax.bar(range(len(tp_order2)), vals, bottom=bottoms,
                  label=label, color=color, alpha=0.85, edgecolor='white')
    # Add label in bar if large enough
    for j, (bar, v, b) in enumerate(zip(bars, vals, bottoms)):
        if v > 5:
            ax.text(j, b + v / 2, f'{v:.0f}%', ha='center', va='center',
                    fontsize=6, color='white' if color not in ['#B0BEC5'] else '#37474F',
                    fontweight='bold')
    bottoms += np.array(vals)

ax.set_xticks(range(len(tp_order2)))
ax.set_xticklabels(tp_labels2, fontsize=6)
ax.set_ylabel('% of 6mA sites', fontsize=8)
ax.set_ylim(0, 105)
ax.set_title('6mA site motif assignment by timepoint', fontsize=8)
ax.legend(loc='upper left', fontsize=6, framealpha=0.8)
add_panel_label(ax, 'd')

plt.tight_layout(rect=[0, 0, 1, 0.97])
fig.suptitle('Figure S14: R-M recognition motif characterization', fontsize=9, fontweight='bold')

save_figure(fig, FIG_SUP_DIR / 'FigS14_motif_reliability.pdf')
