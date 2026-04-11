#!/usr/bin/env python3
"""
Methylation-Expression Concordance Patterns across Three Timepoint Comparisons.

Horizontal grouped bar charts showing concordant (Gained_Up + Lost_Down),
discordant (Lost_Up + Gained_Down), and other (Stable/Mild/Increased/Decreased)
gene counts, split by modification type (4mC vs 6mA).

Streptomyces coelicolor A3(2) M145
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

# === Paths ===
BASE = Path('/Users/okaban/bioinfo/rna-seq/11_epigenome_integration')
OUTPUT_DIR = BASE / 'analysis/02_publication_figures'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

DATA_FILES = {
    'T2 vs T1': BASE / 'analysis/01_integration/T2vsT1_coordinated_genes.csv',
    'T3 vs T1': BASE / 'analysis/16_t3_coordinated/T3vsT1_coordinated_genes.csv',
    'T3 vs T2': BASE / 'analysis/16_t3_coordinated/T3vsT2_coordinated_genes.csv',
}

# === Style ===
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

# === Color scheme ===
# Concordant: warm tones (methylation and expression change in same functional direction)
# Discordant: cool tones
# Other: grays
COLORS = {
    'concordant_4mC': '#D32F2F',   # dark red
    'concordant_6mA': '#FF8A65',   # light orange-red
    'discordant_4mC': '#1565C0',   # dark blue
    'discordant_6mA': '#64B5F6',   # light blue
    'other_4mC':      '#616161',   # dark gray
    'other_6mA':      '#BDBDBD',   # light gray
}

# Category labels
CATEGORY_LABELS = ['Concordant', 'Discordant', 'Other']


def classify_coordination(coordination_value):
    """
    Classify a coordination string into Concordant, Discordant, or Other.

    Concordant (positive correlation - same direction):
        Gained_Up, Lost_Down, Increased_Up, Decreased_Down
    Discordant (negative correlation - opposite direction):
        Lost_Up, Gained_Down, Decreased_Up, Increased_Down
    Other (partial/mild changes):
        *_Mild, Stable_*, or anything else
    """
    c = coordination_value
    # Concordant: methylation gain + expression up, or methylation loss + expression down
    concordant = [
        'Gained_Up', 'Lost_Down',
        'Increased_Up', 'Decreased_Down',
    ]
    # Discordant: methylation loss + expression up, or methylation gain + expression down
    discordant = [
        'Lost_Up', 'Gained_Down',
        'Decreased_Up', 'Increased_Down',
    ]
    if c in concordant:
        return 'Concordant'
    elif c in discordant:
        return 'Discordant'
    else:
        return 'Other'


def load_and_classify(filepath, label):
    """Load a CSV and classify each gene into concordance category."""
    df = pd.read_csv(filepath)
    df['concordance'] = df['coordination'].apply(classify_coordination)
    df['comparison'] = label
    return df


def count_by_mod_and_concordance(df):
    """
    Count genes per (concordance category, mod_type) combination.
    Returns a dict: {(concordance, mod_type): count}
    """
    counts = {}
    for cat in CATEGORY_LABELS:
        for mod in ['4mC', '6mA']:
            n = len(df[(df['concordance'] == cat) & (df['mod_type'] == mod)])
            counts[(cat, mod)] = n
    return counts


def main():
    # --- Load all data ---
    all_data = {}
    for label, fpath in DATA_FILES.items():
        df = load_and_classify(fpath, label)
        all_data[label] = df
        print(f"\n{label}: {len(df)} total genes")
        for cat in CATEGORY_LABELS:
            sub = df[df['concordance'] == cat]
            n4 = len(sub[sub['mod_type'] == '4mC'])
            n6 = len(sub[sub['mod_type'] == '6mA'])
            print(f"  {cat:12s}: 4mC={n4:3d}  6mA={n6:3d}  total={n4+n6:3d}")

    # --- Prepare counts for plotting ---
    comparisons = ['T2 vs T1', 'T3 vs T1', 'T3 vs T2']
    counts_all = {}
    for comp in comparisons:
        counts_all[comp] = count_by_mod_and_concordance(all_data[comp])

    # --- Create figure: three panels side by side ---
    fig, axes = plt.subplots(1, 3, figsize=(11, 4.2), sharey=False)
    fig.subplots_adjust(wspace=0.45, bottom=0.22)

    bar_height = 0.35
    y_positions = np.arange(len(CATEGORY_LABELS))

    for idx, (comp, ax) in enumerate(zip(comparisons, axes)):
        counts = counts_all[comp]

        # Values for each category
        vals_4mC = [counts[(cat, '4mC')] for cat in CATEGORY_LABELS]
        vals_6mA = [counts[(cat, '6mA')] for cat in CATEGORY_LABELS]

        # Compute x-axis limit first (needed for label positioning)
        all_vals = vals_4mC + vals_6mA
        data_max = max(all_vals) if max(all_vals) > 0 else 10
        xmax = data_max * 1.55
        label_offset = data_max * 0.04   # small gap after bar end
        total_offset = data_max * 0.15   # additional gap for n= label

        # Plot horizontal bars: 4mC on top, 6mA on bottom within each category
        colors_4mC = [COLORS[f'{cat.lower()}_4mC'] for cat in CATEGORY_LABELS]
        colors_6mA = [COLORS[f'{cat.lower()}_6mA'] for cat in CATEGORY_LABELS]

        bars_4mC = ax.barh(
            y_positions + bar_height / 2, vals_4mC, bar_height,
            color=colors_4mC, edgecolor='white', linewidth=0.5,
            zorder=3,
        )
        bars_6mA = ax.barh(
            y_positions - bar_height / 2, vals_6mA, bar_height,
            color=colors_6mA, edgecolor='white', linewidth=0.5,
            zorder=3,
        )

        # Add count labels next to bars
        for bar, val in zip(bars_4mC, vals_4mC):
            if val > 0:
                ax.text(
                    bar.get_width() + label_offset,
                    bar.get_y() + bar.get_height() / 2,
                    str(val), va='center', ha='left', fontsize=8,
                    fontweight='bold', color='#333333',
                )
        for bar, val in zip(bars_6mA, vals_6mA):
            if val > 0:
                ax.text(
                    bar.get_width() + label_offset,
                    bar.get_y() + bar.get_height() / 2,
                    str(val), va='center', ha='left', fontsize=8,
                    fontweight='bold', color='#333333',
                )

        # Total counts per category (n= annotation at right)
        for i, cat in enumerate(CATEGORY_LABELS):
            total = vals_4mC[i] + vals_6mA[i]
            max_bar = max(vals_4mC[i], vals_6mA[i])
            ax.text(
                max_bar + total_offset, y_positions[i],
                f'n={total}', va='center', ha='left', fontsize=8,
                color='#555555', fontstyle='italic',
            )

        # Axes formatting
        ax.set_yticks(y_positions)
        ax.set_yticklabels(CATEGORY_LABELS, fontsize=10)
        ax.set_xlabel('Number of genes', fontsize=10)
        ax.set_title(comp, fontsize=12, fontweight='bold', pad=10)
        ax.invert_yaxis()  # Concordant at top

        # Clean spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis='y', length=0)
        ax.grid(axis='x', linestyle='--', alpha=0.3, zorder=0)
        ax.set_xlim(0, xmax)

    # --- Legend ---
    legend_elements = [
        mpatches.Patch(facecolor=COLORS['concordant_4mC'], edgecolor='white',
                       label='Concordant - 4mC'),
        mpatches.Patch(facecolor=COLORS['concordant_6mA'], edgecolor='white',
                       label='Concordant - 6mA'),
        mpatches.Patch(facecolor=COLORS['discordant_4mC'], edgecolor='white',
                       label='Discordant - 4mC'),
        mpatches.Patch(facecolor=COLORS['discordant_6mA'], edgecolor='white',
                       label='Discordant - 6mA'),
        mpatches.Patch(facecolor=COLORS['other_4mC'], edgecolor='white',
                       label='Other - 4mC'),
        mpatches.Patch(facecolor=COLORS['other_6mA'], edgecolor='white',
                       label='Other - 6mA'),
    ]
    fig.legend(
        handles=legend_elements,
        loc='lower center',
        ncol=3,
        bbox_to_anchor=(0.5, 0.0),
        frameon=False,
        fontsize=8.5,
        handlelength=1.5,
        handleheight=1.0,
        columnspacing=1.5,
    )

    fig.suptitle(
        'Methylation-Expression Concordance Patterns',
        fontsize=13, fontweight='bold', y=0.98,
    )

    # --- Save ---
    out_base = OUTPUT_DIR / 'concordance_patterns'
    for fmt in ['pdf', 'svg', 'png']:
        outpath = out_base.with_suffix(f'.{fmt}')
        fig.savefig(outpath, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"\nSaved: {outpath}")

    plt.close(fig)

    # --- Print summary table ---
    print("\n" + "=" * 70)
    print("SUMMARY TABLE")
    print("=" * 70)
    header = f"{'Comparison':<12} {'Category':<14} {'4mC':>5} {'6mA':>5} {'Total':>6} {'%':>7}"
    print(header)
    print("-" * 70)
    for comp in comparisons:
        total_comp = sum(counts_all[comp].values())
        for cat in CATEGORY_LABELS:
            n4 = counts_all[comp][(cat, '4mC')]
            n6 = counts_all[comp][(cat, '6mA')]
            nt = n4 + n6
            pct = 100 * nt / total_comp if total_comp > 0 else 0
            print(f"{comp:<12} {cat:<14} {n4:>5} {n6:>5} {nt:>6} {pct:>6.1f}%")
        print(f"{'':>12} {'TOTAL':<14} {'':>5} {'':>5} {total_comp:>6}")
        print("-" * 70)


if __name__ == '__main__':
    main()
