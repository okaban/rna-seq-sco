"""
Figure 1 panel — GCCGGC core/arm distribution across T1/T2/T3.

Uses GCCGGC site calls from `37_defense_island_GCCGGC` (region pre-labelled
with the canonical arm cutoffs 1.5 Mb and 7.167 Mb). Reports both raw site
counts and Mb-normalised density so the geographic shift across timepoints is
read directly from the bar heights.
"""

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def mm_to_inch(mm):  # NAR width helper
    return mm / 25.4

BASE = Path('/Users/okaban/bioinfo/rna-seq')
SITES = BASE / '11_epigenome_integration/analysis/37_defense_island_GCCGGC/tables/GCCGGC_sites_by_timepoint.tsv'
OUT_MAIN = BASE / '15_paper_figures/figures/main/Figure1_GCCGGC_arm_distribution.png'
OUT_PDF = OUT_MAIN.with_suffix('.pdf')
OUT_SVG = OUT_MAIN.with_suffix('.svg')
OUT_TABLE = BASE / '15_paper_figures/tables/Figure1_GCCGGC_arm_distribution.tsv'
OUT_WRITING = BASE / 'Writing/fig_images/Figure1_arm_panel.png'

# Region sizes for density normalisation (matches defense_island_GCCGGC pipeline)
ARM_LEFT_BP = 1_500_000
ARM_RIGHT_BP = 7_167_507
GENOME_BP = 8_667_507
CORE_MB = (ARM_RIGHT_BP - ARM_LEFT_BP) / 1e6                 # 5.668 Mb
ARM_MB = (ARM_LEFT_BP + (GENOME_BP - ARM_RIGHT_BP)) / 1e6    # 3.000 Mb

TIMEPOINTS = ['T1', 'T2', 'T3']
TP_LABELS = ['T1 (12 h)', 'T2 (24 h)', 'T3 (50 h)']
COL_CORE = '#1565C0'
COL_ARM = '#FB8C00'


def main() -> None:
    sites = pd.read_csv(SITES, sep='\t')
    sites = sites[sites['region'].isin(['core', 'arm'])]

    counts = (sites.groupby(['timepoint', 'region']).size()
              .unstack(fill_value=0)
              .reindex(index=TIMEPOINTS, columns=['core', 'arm'], fill_value=0))

    density = counts.copy().astype(float)
    density['core'] = counts['core'] / CORE_MB
    density['arm'] = counts['arm'] / ARM_MB

    summary = pd.DataFrame({
        'timepoint': TIMEPOINTS,
        'core_sites': counts['core'].values,
        'arm_sites': counts['arm'].values,
        'core_density_per_Mb': density['core'].round(2).values,
        'arm_density_per_Mb': density['arm'].round(2).values,
        'core_arm_ratio': (density['core'] / density['arm'].replace(0, np.nan)).round(2).values,
    })
    OUT_TABLE.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(OUT_TABLE, sep='\t', index=False)

    fig, axes = plt.subplots(1, 2, figsize=(mm_to_inch(174), mm_to_inch(75)))  # NAR full-width
    x = np.arange(len(TIMEPOINTS))
    width = 0.38

    for ax, frame, ylabel, title in (
        (axes[0], counts, 'GCCGGC sites (count)', 'Site count'),
        (axes[1], density, 'GCCGGC sites per Mb', 'Density (per Mb)'),
    ):
        b1 = ax.bar(x - width / 2, frame['core'].values, width,
                    color=COL_CORE, edgecolor='black', linewidth=0.5,
                    label=f'Core ({CORE_MB:.2f} Mb)')
        b2 = ax.bar(x + width / 2, frame['arm'].values, width,
                    color=COL_ARM, edgecolor='black', linewidth=0.5,
                    label=f'Arm ({ARM_MB:.2f} Mb)')
        ax.set_xticks(x)
        ax.set_xticklabels(TP_LABELS)
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=7)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ymax = max(frame.values.max(), 1)
        ax.set_ylim(0, ymax * 1.18)
        for bars in (b1, b2):
            for rect in bars:
                h = rect.get_height()
                if h <= 0:
                    continue
                label = f'{int(h)}' if h >= 5 else f'{h:.1f}'
                ax.text(rect.get_x() + rect.get_width() / 2, h, label,
                        ha='center', va='bottom', fontsize=6)

    axes[0].legend(loc='upper right', frameon=False, fontsize=6)

    fig.suptitle('GCCGGC 4mC sites — core vs arm across timepoints',
                 fontsize=8, y=1.02)
    fig.tight_layout()

    OUT_MAIN.parent.mkdir(parents=True, exist_ok=True)
    OUT_WRITING.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_MAIN, dpi=300, bbox_inches='tight')
    fig.savefig(OUT_PDF, bbox_inches='tight')
    fig.savefig(OUT_SVG, bbox_inches='tight')
    fig.savefig(OUT_WRITING, dpi=300, bbox_inches='tight')
    plt.close(fig)

    print('Saved:')
    for p in (OUT_MAIN, OUT_PDF, OUT_SVG, OUT_WRITING, OUT_TABLE):
        print(' ', p)
    print('\nSummary:')
    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()
