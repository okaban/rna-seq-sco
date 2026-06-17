#!/usr/bin/env python3
"""
Figure 6 (LOCKED reframe, 2026-06-17): Methylation as a permissive spatial organiser.

Replaces the retracted four-layer "Gatekeeper switch" figure (old panels rested on
the expression-selected two-antagonistic-bloc claim and reported a tautological
distance-classifier AUC). The non-circular reframe keeps a single capstone
schematic of the paper's actual thesis:

  T1 "core-methylated" state  --(developmental switch)-->  T2 "arm-redistributed" state
  - GCCGGC m4C concentrated in the active central core / Compartment A (2.3-6.2 Mb) at T1
  - synchronously erased and relocated to the arms at T2, tracking the 3D compartment
    refolding reported by Deng et al. 2023
  - the relocation is PERMISSIVE: it biases expression only weakly
    (geography-controlled r = -0.09, p = 0.004). No AUC, no expression-direction claim.

Outputs:
  15_paper_figures/figures/main/Figure6_spatial_organizer.{pdf,svg,png}
  Writing/fig_images/Figure6.png   (manuscript image slot; no renumbering)
"""
from pathlib import Path
import importlib
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch

sys.path.insert(0, str(Path(__file__).parent))
_utils = importlib.import_module('00_shared_utils')
for _attr in dir(_utils):
    if not _attr.startswith('_'):
        globals()[_attr] = getattr(_utils, _attr)

GENOME = 8_667_507          # NC_003888.3 length (bp)
CORE_LO, CORE_HI = 2.30e6, 6.20e6   # Deng 2023 Compartment A
RNG = np.random.default_rng(7)


def _methyl_ticks(n, lo, hi, concentrate_core, spread):
    """Deterministic tick x-positions (Mb), weighted toward core or arms."""
    if concentrate_core:
        xs = RNG.normal((CORE_LO + CORE_HI) / 2, (CORE_HI - CORE_LO) / 4, n)
    else:
        left = RNG.normal(CORE_LO * 0.5, spread, n // 2)
        right = RNG.normal(CORE_HI + (GENOME - CORE_HI) * 0.5, spread, n - n // 2)
        xs = np.concatenate([left, right])
    return np.clip(xs, 0.05e6, GENOME - 0.05e6) / 1e6


def draw_state(ax, title, sub, concentrate_core, n_marks, exposed_marked,
               show_legend):
    """One chromosome-state panel."""
    ax.set_xlim(-0.3, GENOME / 1e6 + 0.3)
    ax.set_ylim(0, 10)
    ax.axis('off')
    gmb = GENOME / 1e6
    ybar = 4.7
    h = 1.0
    top = ybar + h
    # title + subtitle (kept well clear of the lollipops)
    ax.text(gmb / 2, 9.6, title, ha='center', va='top', fontsize=9,
            fontweight='bold', color=COL_DARK)
    ax.text(gmb / 2, 8.75, sub, ha='center', va='top', fontsize=6.6, color=COL_DARK)
    # active central core / Compartment A label (above lollipops)
    ax.text((CORE_LO + CORE_HI) / 2e6, 8.05,
            'active core / Compartment A (2.3–6.2 Mb)', ha='center', va='center',
            fontsize=6.2, color='#2E6B4F')
    # GCCGGC m4C marks (vertical lollipops above the bar), capped below the label
    xs = _methyl_ticks(n_marks, 0, gmb, concentrate_core, spread=0.9e6)
    for x in xs:
        ax.vlines(x, top, top + 0.85, color=COL_4mC, lw=0.7, alpha=0.85, zorder=3)
        ax.plot(x, top + 0.85, 'o', ms=2.0, color=COL_4mC, alpha=0.9, zorder=3)
    # chromosome backbone
    ax.add_patch(Rectangle((0, ybar), gmb, h, facecolor='#ECEFF1',
                           edgecolor=COL_DARK, linewidth=1.1, zorder=2))
    # active central core / Compartment A shade
    ax.add_patch(Rectangle((CORE_LO / 1e6, ybar), (CORE_HI - CORE_LO) / 1e6, h,
                           facecolor='#CFE3D6', edgecolor='none', zorder=1.5))
    # oriC tick (~ chromosome centre); label below the bar
    ax.plot([gmb / 2, gmb / 2], [ybar - 0.12, top + 0.12], color=COL_DARK,
            lw=0.8, ls=':', zorder=3)
    ax.text(gmb / 2, ybar - 0.45, 'oriC', ha='center', va='top',
            fontsize=6, color=COL_DARK)
    # arm labels (offset from oriC)
    ax.text(CORE_LO / 2e6, ybar - 0.45, 'arm', ha='center', va='top',
            fontsize=6.2, color=COL_DARK)
    ax.text((CORE_HI / 1e6 + gmb) / 2, ybar - 0.45, 'arm', ha='center', va='top',
            fontsize=6.2, color=COL_DARK)
    # Exposed regulators (squares on the backbone)
    ereg = np.clip(RNG.normal((CORE_LO + CORE_HI) / 2, (CORE_HI - CORE_LO) / 4, 7) / 1e6,
                   0.2, gmb - 0.2)
    if exposed_marked:
        for x in ereg:
            ax.plot(x, ybar + h / 2, 's', ms=4.2, color=COL_EXPOSED,
                    markeredgecolor='w', markeredgewidth=0.4, zorder=4)
    else:
        for x in ereg:
            ax.plot(x, ybar + h / 2, 's', ms=4.2, color='none',
                    markeredgecolor=COL_EXPOSED, markeredgewidth=0.9, zorder=4)
    # per-state one-line caption under the bar
    cap = ('■ Exposed regulators carry promoter m4C' if exposed_marked
           else '□ same regulators, promoters demethylated')
    ax.text(gmb / 2, 2.7, cap, ha='center', fontsize=6, color=COL_EXPOSED)
    # shared mark legend (left panel only)
    if show_legend:
        ax.text(0.0, 1.5, '│ GCCGGC m4C site     ■ Exposed regulator',
                fontsize=6, color=COL_DARK, ha='left')


def main():
    apply_style()
    print('=== Figure 6: permissive spatial organiser (reframe) ===')
    fig = plt.figure(figsize=(mm_to_inch(180), mm_to_inch(96)))

    axL = fig.add_axes([0.015, 0.30, 0.45, 0.66])
    axR = fig.add_axes([0.535, 0.30, 0.45, 0.66])

    draw_state(axL, 'T1 — core-methylated state', 'vegetative growth (12 h)',
               concentrate_core=True, n_marks=46, exposed_marked=True, show_legend=True)
    draw_state(axR, 'T2 — arm-redistributed state', 'developmental switch (24 h)',
               concentrate_core=False, n_marks=40, exposed_marked=False, show_legend=False)

    # central transition arrow spanning the gap
    arr = FancyArrowPatch((0.467, 0.60), (0.533, 0.60), transform=fig.transFigure,
                          arrowstyle='-|>', mutation_scale=20, lw=2.2,
                          color=COL_DARK, zorder=5)
    fig.add_artist(arr)
    fig.text(0.50, 0.685, 'developmental\nswitch', ha='center', va='bottom',
             fontsize=6.6, fontweight='bold', color=COL_DARK)
    fig.text(0.50, 0.55, 'synchronous\ndemethylation\n+ core→arm\nrelocation', ha='center',
             va='top', fontsize=5.8, color=COL_DARK)

    # bottom permissive banner (full width)
    ax_b = fig.add_axes([0.015, 0.015, 0.97, 0.235]); ax_b.axis('off')
    ax_b.set_xlim(0, 10); ax_b.set_ylim(0, 10)
    banner = FancyBboxPatch((0.1, 0.6), 9.8, 8.8, boxstyle='round,pad=0.12',
                            facecolor='#FFF8E1', edgecolor='#E0A82E', linewidth=1.3)
    ax_b.add_patch(banner)
    ax_b.text(5.0, 7.2, 'Permissive, not instructive', ha='center', va='center',
              fontsize=8.5, fontweight='bold', color='#B26A00')
    ax_b.text(5.0, 4.6,
              'The core→arm relocation mirrors the active 3D chromosomal compartment '
              '(Deng et al. 2023) but biases expression only weakly\n'
              '(geography-controlled $r = -0.09$, $p = 0.004$). Methylation marks '
              'which regulatory loci are spatially organised —\n'
              'it sets transcriptional competence, not the direction of transcription. '
              'No predictive classifier (AUC) is claimed.',
              ha='center', va='center', fontsize=6.6, color=COL_DARK)

    out = FIG_DIR / 'Figure6_spatial_organizer'
    save_figure(fig, out, formats=('pdf', 'svg', 'png'))

    # also write the manuscript image slot (keep Figure6.png; no renumbering).
    # The manuscript writing layer lives in the Obsidian repo, NOT under bioinfo BASE;
    # prefer it if present, else fall back to the bioinfo-side mirror.
    import shutil
    png = out.with_suffix('.png')
    obsidian_slot = Path.home() / 'obsidian' / 'Research' / 'rna-seq' / 'Writing' / 'fig_images' / 'Figure6.png'
    targets = [obsidian_slot] if obsidian_slot.parent.is_dir() else []
    targets.append(BASE / 'Writing' / 'fig_images' / 'Figure6.png')   # bioinfo-side mirror
    for wr in targets:
        wr.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(png, wr)
        print(f'  Copied → {wr}')
    print('=== Done ===')


if __name__ == '__main__':
    main()
