#!/usr/bin/env python3
"""
Figure 5: Integrated Gatekeeper Architecture Model

4-panel concept figure (2×2):
  A: Chromosome architecture – Core vs Arms, T1→T2 methylation redistribution
  B: Gatekeeper protection zone – TSS-centred 293 bp methylation-free window
  C: 57 Exposed TF 2-block expression switch (heatmap of z-scores)
  D: Integrated model flow (4-layer cascade with key statistics)

Output:
  15_paper_figures/figures/main/Figure5_gatekeeper_model.png  (300 dpi)
  Writing/fig_images/Figure5_gatekeeper_model.png
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec

# ── Paths ────────────────────────────────────────────────────────────────────
import os, sys
# Support running from both host and sandbox mount paths
_host_base = Path('/Users/okaban/bioinfo/rna-seq')
_mnt_base  = Path('/sessions/vigilant-laughing-maxwell/mnt/rna-seq')
BASE       = _mnt_base if _mnt_base.exists() else _host_base
EPIBASE    = BASE / '11_epigenome_integration' / 'analysis'
FIG_MAIN = BASE      / '15_paper_figures' / 'figures' / 'main'
FIG_WR   = BASE      / 'Writing' / 'fig_images'
# Output paths always use host paths for final files
HOST_FIG_MAIN = _host_base / '15_paper_figures' / 'figures' / 'main'
HOST_FIG_WR   = _host_base / 'Writing' / 'fig_images'

# ── Colour palette ────────────────────────────────────────────────────────────
COL_CORE    = '#2196F3'   # blue  – core (methylation-rich T1)
COL_ARM     = '#FF5722'   # deep orange – arms (methylation-rich T2)
COL_PROTECT = '#4CAF50'   # green – protection zone / shielded
# COL_EXPOSED inherited from shared_utils as '#7E57C2' (purple)
COL_SOLO    = '#9C27B0'   # purple – solo MTase / AAGCCCG
COL_BG      = '#F5F5F5'   # light grey background for boxes
COL_ARROW   = '#455A64'   # dark slate for arrows

FONT = 'Arial'
plt.rcParams.update({
    'font.family':       FONT,
    'font.size':         8,
    'axes.labelsize':    8,
    'axes.titlesize':    9,
    'xtick.labelsize':   7,
    'ytick.labelsize':   7,
    'axes.linewidth':    0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'pdf.fonttype':      42,
    'svg.fonttype':      'none',
})

# ── Load data ─────────────────────────────────────────────────────────────────
def load_temporal_data():
    """Load 57 exposed TF z-score expression data."""
    path = EPIBASE / '57_temporal_dynamics_exposed_TF' / 'tables' / 'temporal_classification.tsv'
    df = pd.read_csv(path, sep='\t')
    # Sort: bloc (activation first), then module, then LFC_T3vsT1
    bloc_order = {'activation': 0, 'repression': 1}
    df['_bloc_rank'] = df['bloc'].map(bloc_order).fillna(2)
    df = df.sort_values(['_bloc_rank', 'module', 'LFC_T3vsT1'], ascending=[True, True, False])
    return df


# ─────────────────────────────────────────────────────────────────────────────
# Panel A: Chromosome architecture diagram
# ─────────────────────────────────────────────────────────────────────────────
def draw_panel_a(ax):
    ax.set_xlim(0, 8.7)
    ax.set_ylim(-1.8, 2.5)
    ax.axis('off')

    chrom_y   = 1.0
    chrom_h   = 0.38
    core_s    = 1.5
    core_e    = 6.5

    # ── Chromosome bar ──────────────────────────────────────────────────────
    # Arms (left + right)
    ax.add_patch(FancyBboxPatch((0, chrom_y), core_s, chrom_h,
                                boxstyle='round,pad=0.04', fc=COL_ARM,
                                ec='white', lw=1.0, zorder=3))
    ax.add_patch(FancyBboxPatch((core_e, chrom_y), 8.7 - core_e, chrom_h,
                                boxstyle='round,pad=0.04', fc=COL_ARM,
                                ec='white', lw=1.0, zorder=3))
    # Core
    ax.add_patch(FancyBboxPatch((core_s, chrom_y), core_e - core_s, chrom_h,
                                boxstyle='round,pad=0.04', fc=COL_CORE,
                                ec='white', lw=1.0, zorder=3))

    # Boundary dashed lines
    for x in [core_s, core_e]:
        ax.axvline(x=x, ymin=0.32, ymax=0.76, color='white', lw=1.2,
                   ls='--', zorder=4, alpha=0.85)

    # Region labels inside bar
    ax.text(core_s / 2, chrom_y + chrom_h / 2, 'Left arm\n<1.5 Mb',
            ha='center', va='center', fontsize=6.5, color='white',
            fontweight='bold', zorder=5)
    ax.text((core_s + core_e) / 2, chrom_y + chrom_h / 2,
            'Core  (1.5 – 6.5 Mb)',
            ha='center', va='center', fontsize=7, color='white',
            fontweight='bold', zorder=5)
    ax.text((core_e + 8.7) / 2, chrom_y + chrom_h / 2, 'Right arm\n>6.5 Mb',
            ha='center', va='center', fontsize=6.5, color='white',
            fontweight='bold', zorder=5)

    # ── Scale bar ───────────────────────────────────────────────────────────
    ax.annotate('', xy=(8.7, chrom_y - 0.25), xytext=(0, chrom_y - 0.25),
                arrowprops=dict(arrowstyle='<->', color='#546E7A', lw=1.0))
    ax.text(8.7 / 2, chrom_y - 0.45, '8.7 Mb', ha='center', va='top',
            fontsize=7, color='#546E7A')

    # ── Methylation density profiles ─────────────────────────────────────────
    x_pos = np.linspace(0, 8.7, 300)

    # T1 methylation: Gaussian centred on core
    def t1_profile(x):
        c = np.exp(-0.5 * ((x - 4.1) / 1.8) ** 2)
        return 0.15 + 0.65 * c

    # T2 methylation: bimodal (arm peaks)
    def t2_profile(x):
        left  = np.exp(-0.5 * ((x - 0.65) / 0.55) ** 2)
        right = np.exp(-0.5 * ((x - 8.05) / 0.55) ** 2)
        core_base = 0.10 * np.exp(-0.5 * ((x - 4.1) / 1.2) ** 2)
        return 0.10 + 0.58 * (left + right) + core_base

    y_t1 = t1_profile(x_pos)
    y_t2 = t2_profile(x_pos)

    profile_base = chrom_y + chrom_h + 0.05
    scale = 0.65

    ax.plot(x_pos, profile_base + y_t1 * scale, color=COL_CORE, lw=1.8,
            label='T1 (vegetative)', zorder=3)
    ax.fill_between(x_pos, profile_base, profile_base + y_t1 * scale,
                    color=COL_CORE, alpha=0.15, zorder=2)

    ax.plot(x_pos, profile_base + y_t2 * scale, color=COL_ARM, lw=1.8,
            ls='--', label='T2 (transition)', zorder=3)
    ax.fill_between(x_pos, profile_base, profile_base + y_t2 * scale,
                    color=COL_ARM, alpha=0.12, zorder=2)

    ax.text(8.7 + 0.08, profile_base + y_t1[-1] * scale, 'T1',
            va='center', fontsize=7, color=COL_CORE, fontweight='bold')
    ax.text(8.7 + 0.08, profile_base + y_t2[-1] * scale, 'T2',
            va='center', fontsize=7, color=COL_ARM, fontweight='bold')

    # ── Redistribution arrows ──────────────────────────────────────────────
    arrow_y = chrom_y + chrom_h * 0.5
    for (xs, xe, label) in [
        (0.75, 0.02, None),
        (7.95, 8.68, None),
    ]:
        ax.annotate('', xy=(xe, arrow_y), xytext=(xs, arrow_y),
                    arrowprops=dict(arrowstyle='->', color=COL_ARM, lw=1.5),
                    zorder=6)

    # ── Statistics boxes ───────────────────────────────────────────────────
    stats_y = -0.35
    for (x, text, col) in [
        (2.15, 'T1: Core-dominant\n76% GCCGGC sites', COL_CORE),
        (0.75, 'T2: Arm redistribution\nJaccard = 0.000', COL_ARM),
        (6.1, 'Solo MTase SC_RS17645\nAAGCCCG / 6mA  (stable)', COL_SOLO),
    ]:
        ax.text(x, stats_y, text, ha='center', va='top', fontsize=6.2,
                color=col, fontweight='bold',
                bbox=dict(fc='white', ec=col, lw=0.8, boxstyle='round,pad=0.25',
                          alpha=0.92))

    # ── Panel label ───────────────────────────────────────────────────────
    ax.text(-0.35, 2.45, 'A', fontsize=12, fontweight='bold', va='top')
    ax.set_title('Genome Architecture & RM System Dynamics', fontsize=8.5,
                 fontweight='bold', pad=4)


# ─────────────────────────────────────────────────────────────────────────────
# Panel B: Gatekeeper protection zone
# ─────────────────────────────────────────────────────────────────────────────
def draw_panel_b(ax):
    ax.set_xlim(-1600, 1600)
    ax.set_ylim(-1.4, 3.0)
    ax.axis('off')

    # ── Gene body ──────────────────────────────────────────────────────────
    gene_y   = 1.35
    gene_h   = 0.30
    tss_x    = 0
    gene_end = 1400

    # Promoter / upstream intergenic
    ax.add_patch(FancyBboxPatch((-1500, gene_y), 1500, gene_h,
                                boxstyle='round,pad=0.0', fc='#B0BEC5',
                                ec='#78909C', lw=0.8, zorder=2))
    # Gene body
    ax.add_patch(FancyBboxPatch((0, gene_y), gene_end, gene_h,
                                boxstyle='round,pad=0.0', fc='#90CAF9',
                                ec='#1565C0', lw=0.8, zorder=2))
    # TSS marker
    ax.annotate('TSS', xy=(0, gene_y + gene_h), xytext=(0, gene_y + gene_h + 0.28),
                arrowprops=dict(arrowstyle='->', color='#1565C0', lw=1.2),
                ha='center', va='bottom', fontsize=7, color='#1565C0',
                fontweight='bold')
    # Gene label
    ax.text(700, gene_y + gene_h / 2, 'Control gene', ha='center', va='center',
            fontsize=6.5, color='#1565C0', fontweight='bold')

    # ── GCCGGC site density (schematic) ────────────────────────────────────
    x_bp = np.linspace(-1500, 1500, 400)

    def methyl_density(x):
        # Depletion around TSS, elevated upstream and in gene body
        depletion = np.exp(-0.5 * (x / 380) ** 2)
        upstream  = 0.55 * np.exp(-0.5 * ((x + 900) / 350) ** 2)
        body      = 0.45 * (1 / (1 + np.exp(-(x - 500) / 200)))
        return 0.12 + 0.70 * (1 - depletion) * 0.6 + upstream + body * 0.25

    density = methyl_density(x_bp)
    # Normalise 0-1
    density = (density - density.min()) / (density.max() - density.min())

    base_y = 0.05
    scale  = 0.82
    ax.plot(x_bp, base_y + density * scale, color='#546E7A', lw=1.5, zorder=4)
    ax.fill_between(x_bp, base_y, base_y + density * scale,
                    color='#90A4AE', alpha=0.35, zorder=3)

    # ── Protection zone highlight ──────────────────────────────────────────
    zone_half = 146   # ≈ 293 bp / 2
    pz = mpatches.FancyBboxPatch((-zone_half, base_y - 0.05),
                                  zone_half * 2, scale + 0.12,
                                  boxstyle='round,pad=0.0',
                                  fc=COL_PROTECT, ec=COL_PROTECT,
                                  lw=1.5, alpha=0.18, zorder=1)
    ax.add_patch(pz)
    # Zone boundary lines
    for xz in [-zone_half, zone_half]:
        ax.axvline(x=xz, ymin=0.02, ymax=0.57, color=COL_PROTECT,
                   lw=1.2, ls='-', zorder=5, alpha=0.75)

    # Bracket + label
    bracket_y = base_y + scale + 0.18
    ax.annotate('', xy=(zone_half, bracket_y), xytext=(-zone_half, bracket_y),
                arrowprops=dict(arrowstyle='<->', color=COL_PROTECT, lw=1.3))
    ax.text(0, bracket_y + 0.16, 'Methylation-free zone\n293 bp  (AUC = 0.917)',
            ha='center', va='bottom', fontsize=7, color=COL_PROTECT,
            fontweight='bold')

    # ── GCCGGC site ticks (schematic) ──────────────────────────────────────
    rng = np.random.default_rng(42)
    # Outside zone – sites present
    for region_x in np.concatenate([
        rng.choice(np.arange(-1490, -zone_half, 40), 14, replace=False),
        rng.choice(np.arange(zone_half + 40, 1490, 40), 12, replace=False),
    ]):
        tick_h = 0.04 + 0.08 * rng.random()
        ax.plot([region_x, region_x],
                [base_y + methyl_density(np.array([region_x]))[0] * scale * 0.7,
                 base_y + methyl_density(np.array([region_x]))[0] * scale * 0.7 + tick_h],
                color='#37474F', lw=0.8, alpha=0.5, zorder=6)

    # ── y-axis label ──────────────────────────────────────────────────────
    ax.text(-1520, base_y + scale / 2, 'GCCGGC\ndensity', ha='right',
            va='center', fontsize=6.5, color='#546E7A')

    # ── Stats box ─────────────────────────────────────────────────────────
    ax.text(0, -0.85,
            '998 shielded genes  ·  57 exposed TFs  (OR = 921×)',
            ha='center', va='center', fontsize=7, color='#455A64',
            bbox=dict(fc='white', ec='#90A4AE', lw=0.8,
                      boxstyle='round,pad=0.3', alpha=0.95))

    ax.text(-1580, 2.90, 'B', fontsize=12, fontweight='bold', va='top')
    ax.set_title('Gatekeeper Protection Zone (TSS ± 293 bp)', fontsize=8.5,
                 fontweight='bold', pad=4)


# ─────────────────────────────────────────────────────────────────────────────
# Panel C: 57 Exposed TF expression heatmap
# ─────────────────────────────────────────────────────────────────────────────
def draw_panel_c(ax):
    df = load_temporal_data()

    # z-scores across T1, T2, T3
    z = df[['T1_z', 'T2_z', 'T3_z']].values   # (n_TF, 3)

    # Custom blue-white-red diverging colourmap
    cmap = LinearSegmentedColormap.from_list(
        'bwr_custom',
        ['#1565C0', '#5C9BD6', '#FAFAFA', '#EF9A9A', '#C62828'],
        N=256)

    vmax = 2.2
    im = ax.imshow(z, aspect='auto', cmap=cmap, vmin=-vmax, vmax=vmax,
                   interpolation='nearest')

    # Bloc divider
    n_act = (df['bloc'] == 'activation').sum()
    ax.axhline(y=n_act - 0.5, color='white', lw=1.5, zorder=5)

    # x-axis
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(['T1\n(veg.)', 'T2\n(trans.)', 'T3\n(dev.)'],
                       fontsize=7)
    ax.tick_params(axis='x', length=0)

    # y-axis – bloc labels
    ax.set_yticks([n_act / 2 - 0.5, n_act + (len(df) - n_act) / 2 - 0.5])
    ax.set_yticklabels(['Activation\nbloc', 'Repression\nbloc'],
                       fontsize=7, fontweight='bold')
    ax.tick_params(axis='y', length=0)

    # Bloc brackets on right side
    for (y0, y1, label, col) in [
        (-0.5, n_act - 0.5, f'{n_act} TFs\n↑ T3', COL_EXPOSED),
        (n_act - 0.5, len(df) - 0.5, f'{len(df)-n_act} TFs\n↓ T3', COL_CORE),
    ]:
        ax.annotate('', xy=(3.38, y0), xytext=(3.38, y1),
                    xycoords='data', textcoords='data',
                    arrowprops=dict(arrowstyle='<->', color=col, lw=1.2),
                    annotation_clip=False)
        ax.text(3.55, (y0 + y1) / 2, label, va='center', ha='left',
                fontsize=6.2, color=col, fontweight='bold',
                clip_on=False)

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.01, shrink=0.65,
                        orientation='vertical')
    cbar.set_label('z-score', fontsize=6.5)
    cbar.set_ticks([-2, 0, 2])
    cbar.ax.tick_params(labelsize=6)

    # Correlation annotation
    ax.text(1, -3.5, r'ρ = −0.995  (2-bloc eigengene switch)',
            ha='center', va='top', fontsize=7, color='#455A64',
            fontweight='bold', clip_on=False)

    ax.text(-0.55, -0.8, 'C', fontsize=12, fontweight='bold', va='top',
            transform=ax.transData, clip_on=False)
    ax.set_title('57 Exposed TF Expression Switch', fontsize=8.5,
                 fontweight='bold', pad=4)


# ─────────────────────────────────────────────────────────────────────────────
# Panel D: Integrated 4-layer model flow
# ─────────────────────────────────────────────────────────────────────────────
def draw_panel_d(ax):
    ax.set_xlim(0, 10)
    ax.set_ylim(-0.3, 10.5)
    ax.axis('off')

    # ── Layer definitions ─────────────────────────────────────────────────
    layers = [
        {
            'y': 8.8, 'label': 'Layer 1\nGenome Architecture',
            'col': COL_CORE, 'text_col': 'white',
            'detail': '8.7 Mb linear chromosome\nCore (1.5–6.5 Mb) vs Arms (<1.5, >6.5 Mb)',
        },
        {
            'y': 6.5, 'label': 'Layer 2\nRM System (GCCGGC / m4C)',
            'col': COL_ARM, 'text_col': 'white',
            'detail': 'T1: Core-dominant  (76%)\nT2: Arm redistribution  (Jaccard = 0.000)\nSolo MTase SC_RS17645: AAGCCCG/6mA (stable)',
        },
        {
            'y': 4.2, 'label': 'Layer 3\nGatekeeper Protection Zone',
            'col': COL_PROTECT, 'text_col': 'white',
            'detail': 'TSS methylation-free zone  (293 bp, AUC = 0.917)\n998 shielded genes — immune to redistribution\n57 Exposed TFs — lose protection at T2 → switch',
        },
        {
            'y': 1.9, 'label': 'Layer 4\nTranscriptional Output',
            'col': COL_EXPOSED, 'text_col': 'white',
            'detail': '57 TF 2-bloc eigengene switch  (ρ = −0.995)\nTetR enrichment  (p = 0.028)\nSecondary metabolism activation',
        },
    ]

    box_w = 8.4
    box_h = 1.65
    box_x = 0.8

    for lay in layers:
        y0 = lay['y']
        # Shadow
        ax.add_patch(FancyBboxPatch((box_x + 0.07, y0 - 0.07),
                                    box_w, box_h,
                                    boxstyle='round,pad=0.18',
                                    fc='#90A4AE', ec='none', alpha=0.25,
                                    zorder=1))
        # Main box
        ax.add_patch(FancyBboxPatch((box_x, y0), box_w, box_h,
                                    boxstyle='round,pad=0.18',
                                    fc=lay['col'], ec='white', lw=1.2,
                                    zorder=2))
        # Layer label (left side)
        ax.text(box_x + 0.28, y0 + box_h / 2, lay['label'],
                ha='left', va='center', fontsize=7.5, color=lay['text_col'],
                fontweight='bold', zorder=3, linespacing=1.35)
        # Detail text (right side)
        ax.text(box_x + 3.1, y0 + box_h / 2, lay['detail'],
                ha='left', va='center', fontsize=6.5, color=lay['text_col'],
                zorder=3, linespacing=1.4)
        # Separator line
        ax.plot([box_x + 2.85, box_x + 2.85], [y0 + 0.15, y0 + box_h - 0.15],
                color='white', lw=0.8, alpha=0.5, zorder=4)

    # ── Arrows with key statistics ────────────────────────────────────────
    arrow_specs = [
        # (from_y, to_y, stat_label, stat_col)
        (layers[0]['y'], layers[1]['y'] + box_h,
         'methylation landscape\nT1 → T2 shift', COL_ARM),
        (layers[1]['y'], layers[2]['y'] + box_h,
         'OR = 921×\nprotection enrichment', COL_PROTECT),
        (layers[2]['y'], layers[3]['y'] + box_h,
         'AUC = 0.917\nexposure selectivity', COL_EXPOSED),
    ]

    for (y_from, y_to, stat, col) in arrow_specs:
        mid_y = (y_from + y_to) / 2
        # Arrow
        ax.annotate('', xy=(5.0, y_to + 0.05), xytext=(5.0, y_from - 0.05),
                    arrowprops=dict(
                        arrowstyle='->', color=col, lw=1.8,
                        connectionstyle='arc3,rad=0.0',
                        mutation_scale=14))
        # Stat label bubble
        ax.text(6.35, mid_y, stat, ha='left', va='center', fontsize=6.2,
                color=col, fontweight='bold',
                bbox=dict(fc='white', ec=col, lw=0.8,
                          boxstyle='round,pad=0.22', alpha=0.95))

    # ── Output arrow at bottom ────────────────────────────────────────────
    output_y = layers[-1]['y'] - 0.12
    ax.annotate('', xy=(5.0, output_y - 0.65), xytext=(5.0, output_y),
                arrowprops=dict(arrowstyle='->', color='#455A64', lw=1.8,
                                mutation_scale=14))
    ax.text(5.0, output_y - 0.82,
            'Secondary metabolite activation  &  developmental reprogramming',
            ha='center', va='top', fontsize=7, color='#455A64',
            fontweight='bold',
            bbox=dict(fc='#FFF8E1', ec='#FFC107', lw=0.9,
                      boxstyle='round,pad=0.28', alpha=0.95))

    ax.text(-0.1, 10.45, 'D', fontsize=12, fontweight='bold', va='top')
    ax.set_title('Integrated Gatekeeper Architecture', fontsize=8.5,
                 fontweight='bold', pad=4)


# ─────────────────────────────────────────────────────────────────────────────
# Main: compose figure
# ─────────────────────────────────────────────────────────────────────────────
def main():
    fig = plt.figure(figsize=(14, 11))
    fig.patch.set_facecolor('white')

    gs = gridspec.GridSpec(
        2, 2,
        figure=fig,
        left=0.06, right=0.97,
        top=0.95, bottom=0.04,
        hspace=0.42, wspace=0.35,
    )

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    draw_panel_a(ax_a)
    draw_panel_b(ax_b)
    draw_panel_c(ax_c)
    draw_panel_d(ax_d)

    # ── Figure title & footnote ───────────────────────────────────────────
    fig.suptitle(
        'Figure 5  |  Integrated Gatekeeper Architecture Model\n'
        r'$\it{S. coelicolor}$ A3(2) M145 Epigenome–Transcriptome Integration',
        fontsize=10, fontweight='bold', y=0.99, va='top',
    )

    # ── Save ─────────────────────────────────────────────────────────────
    FIG_MAIN.mkdir(parents=True, exist_ok=True)
    FIG_WR.mkdir(parents=True, exist_ok=True)

    out_main = FIG_MAIN / 'Figure5_gatekeeper_model.png'
    out_wr   = FIG_WR   / 'Figure5_gatekeeper_model.png'

    fig.savefig(out_main, dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    fig.savefig(out_wr,   dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')

    print(f'Saved → {out_main}')
    print(f'Saved → {out_wr}')
    plt.close(fig)


if __name__ == '__main__':
    main()
