"""
FigS15: Genome-wide window methylation dynamics (B-1)

Panels:
  a. 4mC site density across genome (50kb windows, T1/T2/T3 overlay + T3-T1 delta)
  b. 6mA site density across genome (same)
  c. Core vs Arm density across timepoints (summary bar chart)
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

WINDOW_DIR = Path('/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/63_window_methylation')
TABLES = WINDOW_DIR / 'tables'

apply_style()

# ── Load data ────────────────────────────────────────────────────────────────
df_4mc = pd.read_csv(TABLES / 'B1_window_methylation_4mC.tsv', sep='\t')
df_6ma = pd.read_csv(TABLES / 'B1_window_methylation_6mA.tsv', sep='\t')
core_arm = pd.read_csv(TABLES / 'B1_core_arm_density_summary.tsv', sep='\t')

GENOME_MB = 8.67507
ARM_LEFT = 2.0
ARM_RIGHT = 6.0
WINDOW_MB = 0.05

TP_COLORS = {'T1': '#4C72B0', 'T2': '#DD8452', 'T3': '#55A868'}
TP_LABELS = {'T1': 'T1 (12 h)', 'T2': 'T2 (24 h)', 'T3': 'T3 (50 h)'}
MOD_COLORS = {'4mC': '#E53935', '6mA': '#1565C0'}

# ── Figure layout ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(mm_to_inch(174), mm_to_inch(214)))  # NAR full-width
gs = fig.add_gridspec(3, 2, hspace=0.55, wspace=0.30,
                      height_ratios=[2.5, 2.5, 2.0])

# Helper: shade arm regions
def shade_arms(ax):
    ax.axvspan(0, ARM_LEFT, color='#CFD8DC', alpha=0.3, zorder=0)
    ax.axvspan(ARM_RIGHT, GENOME_MB, color='#CFD8DC', alpha=0.3, zorder=0)
    ax.text(ARM_LEFT / 2, ax.get_ylim()[1] * 0.92, 'arm', ha='center', fontsize=6, color='gray')
    ax.text((ARM_RIGHT + GENOME_MB) / 2, ax.get_ylim()[1] * 0.92, 'arm', ha='center', fontsize=6, color='gray')
    ax.text((ARM_LEFT + ARM_RIGHT) / 2, ax.get_ylim()[1] * 0.92, 'core', ha='center', fontsize=6, color='gray')

def plot_density_panel(ax, df, mod, title, panel_label):
    x = df['center'].values / 1e6
    for tp in ['T1', 'T2', 'T3']:
        y = df[f'density_{tp}'].values
        ax.plot(x, y, color=TP_COLORS[tp], linewidth=1.2, label=TP_LABELS[tp])
        ax.fill_between(x, y, alpha=0.12, color=TP_COLORS[tp])
    ax.set_ylabel('Sites per kb', fontsize=7)
    ax.set_title(f'{mod} methylation density', fontsize=8)
    ax.legend(fontsize=6, loc='upper right')
    ax.set_xlim(0, GENOME_MB)
    ax.set_ylim(bottom=0)
    shade_arms(ax)
    add_panel_label(ax, panel_label)

def plot_delta_panel(ax, df, mod, panel_label):
    x = df['center'].values / 1e6
    delta = df['density_T3'].values - df['density_T1'].values
    colors = ['#C62828' if v > 0 else '#1565C0' for v in delta]
    ax.bar(x, delta, width=WINDOW_MB * 0.85, color=colors, alpha=0.8)
    ax.axhline(0, color='black', linewidth=0.7)
    ax.set_ylabel('Δ density (T3−T1)', fontsize=7)
    ax.set_title(f'{mod} density change (T3−T1)', fontsize=8)
    ax.set_xlim(0, GENOME_MB)
    shade_arms(ax)
    add_panel_label(ax, panel_label)

# --- Panels a/b: density overlay ---
ax_4mc_dens = fig.add_subplot(gs[0, 0])
ax_6ma_dens = fig.add_subplot(gs[0, 1])
ax_4mc_del = fig.add_subplot(gs[1, 0])
ax_6ma_del = fig.add_subplot(gs[1, 1])

plot_density_panel(ax_4mc_dens, df_4mc, '4mC', '4mC methylation density', 'a')
plot_density_panel(ax_6ma_dens, df_6ma, '6mA', '6mA methylation density', 'b')
plot_delta_panel(ax_4mc_del, df_4mc, '4mC', 'c')
plot_delta_panel(ax_6ma_del, df_6ma, '6mA', 'd')

for ax in [ax_4mc_del, ax_6ma_del]:
    ax.set_xlabel('Genomic position (Mb)', fontsize=7)
# Refresh arm shading after ylim is set
for ax, df, mod in [(ax_4mc_dens, df_4mc, '4mC'), (ax_6ma_dens, df_6ma, '6mA')]:
    shade_arms(ax)

# --- Panel e: core vs arm summary ---
ax_ca = fig.add_subplot(gs[2, :])
tp_order = ['T1', 'T2', 'T3']
tp_x = np.arange(len(tp_order))
width = 0.18
mods = ['4mC', '6mA']
regions = ['core', 'arm']
region_hatches = {'core': '', 'arm': '//'}
region_alphas = {'core': 0.9, 'arm': 0.6}
mod_colors_bar = {'4mC': '#E53935', '6mA': '#1565C0'}

offsets = [-1.5 * width, -0.5 * width, 0.5 * width, 1.5 * width]
combos = [('4mC', 'core'), ('4mC', 'arm'), ('6mA', 'core'), ('6mA', 'arm')]
labels = ['4mC core', '4mC arm', '6mA core', '6mA arm']
colors = [mod_colors_bar['4mC'], mod_colors_bar['4mC'], mod_colors_bar['6mA'], mod_colors_bar['6mA']]
alphas = [0.9, 0.5, 0.9, 0.5]
hatches = ['', '////', '', '////']

for (mod, region), offset, label, color, alpha, hatch in zip(combos, offsets, labels, colors, alphas, hatches):
    sub = core_arm[(core_arm['mod_type'] == mod) & (core_arm['region'] == region)]
    sub = sub.set_index('timepoint')
    vals = [sub.loc[tp, 'mean_density'] if tp in sub.index else 0 for tp in tp_order]
    ax_ca.bar(tp_x + offset, vals, width, label=label,
              color=color, alpha=alpha, hatch=hatch, edgecolor='white')

ax_ca.set_xticks(tp_x)
ax_ca.set_xticklabels(['T1 (12 h)', 'T2 (24 h)', 'T3 (50 h)'], fontsize=7)
ax_ca.set_ylabel('Mean density (sites/kb)', fontsize=8)
ax_ca.set_title('Methylation density: core vs arm regions', fontsize=8)
ax_ca.legend(fontsize=6, ncol=4, loc='upper right', framealpha=0.9)
add_panel_label(ax_ca, 'e')

fig.suptitle('Figure S15: Genome-wide methylation density dynamics (50 kb windows)',
             fontsize=9, fontweight='bold', y=1.01)

save_figure(fig, FIG_SUP_DIR / 'FigS15_window_methylation.pdf')
