#!/usr/bin/env python3
"""Figure 2 — GCCGGC 4mC redistribution and promoter protection architecture.

Panel A: Genome-wide GCCGGC 4mC density at T1/T2/T3 (50 kb windows).
Panel B: TSS metagene methylation profile, T1, ±4 kb, all genes.
Panel C: Continuous TSS-to-nearest-GCCGGC-4mC distance across all regulatory
         genes (n=1,051; 62 Exposed / 989 Shielded), shown as a single histogram
         + KDE on a log10-bp axis with the 293 bp operating point marked. The
         distribution is unimodal (Hartigan's dip test does not reject), so it is
         a gradient, not two natural groups.
Panel D: |log2 fold-change| variability of Shielded vs Exposed TFs at T2 vs T1
         and T3 vs T1.

Outputs PDF/SVG/PNG to FIG_DIR. Adds uppercase panel labels A/B/C/D and
positions panel D's legend so it does not overlap the significance markers.
"""

import sys
import importlib
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
_utils = importlib.import_module('00_shared_utils')
for _attr in dir(_utils):
    if not _attr.startswith('_'):
        globals()[_attr] = getattr(_utils, _attr)


GENOME_LEN = 8_667_507
ARM_LEFT = 1_500_000
ARM_RIGHT = 7_167_507
PROTECTION_ZONE_BP = 293

# ── Unified palette (Rich & Calm, 2026-07) ───────────────────────────────────
# 4mC = rich muted red, 6mA = rich muted blue (identical across all figures).
# Timepoint line series use the grayscale ramp (T1 light → T3 dark), matching
# the SuppFig15 bar convention, since panel A is all GCCGGC 4mC across time.
COL_4mC = '#A64B44'      # unified rich muted red — 4mC
COL_6mA = '#3A6B8C'      # unified rich muted blue — 6mA
COL_GREY = '#9AA7B0'     # unified neutral grey baseline / core shading
# Fig1a core/arm shading + oriC (imported constants, not re-picked by eye) so
# Fig2 panel A matches Figure 1's colour tone exactly (figure-legibility-qc §6).
COL_CORE = '#3E7256'     # deep muted green — chromosomal core (translucent band)
COL_ARM = '#C0803A'      # calm amber — chromosomal arms (translucent band)
ORIC_POS = 4_270_777     # dnaA; drawn grey per canon (never green)
COL_ORIC = '#333333'

T_LABELS = {'T1': 'T1 (12 h)', 'T2': 'T2 (24 h)', 'T3': 'T3 (50 h)'}


def add_uppercase_label(ax, label, x=None, y=None, fontsize=13):
    """Delegate to the shared uniform panel-label helper (constant offset-points
    placement, lowercase). Name/signature kept for call-site compatibility."""
    _utils.add_panel_label(ax, label, fontsize=fontsize)


def panel_a(ax, df_sites, bin_size_bp=50_000):
    """GCCGGC 4mC core-to-arm redistribution — stacked per-timepoint landscape.

    Formerly a 3-line overlay in an Okabe-Ito trio (blue/orange/green), which
    introduced three colours that clash with the paper's 4mC-red / 6mA-blue /
    core-green / arm-amber scheme (figure-legibility-qc §6). Redrawn as three
    stacked density tracks (T1 top → T3 bottom), ALL in the canonical 4mC red,
    with Figure 1a's translucent core/arm bands and grey oriC line, so the
    developmental redistribution reads by vertical track rather than by an
    extra colour dimension. This is the both-motif SuppFig11 panel-a design,
    restricted to GCCGGC to stay coherent with panels B/C.
    """
    edges = np.arange(0, GENOME_LEN + bin_size_bp, bin_size_bp)
    centers_mb = (edges[:-1] + edges[1:]) / 2 / 1e6
    tps = ['T1', 'T2', 'T3']
    band_h, ntp = 0.80, 3
    top = ntp

    # Per-timepoint max for a shared vertical scale across tracks.
    hmax = max((np.histogram(df_sites[df_sites['timepoint'] == tp]['position'].values,
                             bins=edges)[0].max() or 1) for tp in tps)

    # Fig1a translucent genome-region bands spanning all three tracks.
    ax.add_patch(plt.Rectangle((0, 0), ARM_LEFT / 1e6, top,
                               fc=COL_ARM, alpha=0.25, ec='none', zorder=0))
    ax.add_patch(plt.Rectangle((ARM_LEFT / 1e6, 0), (ARM_RIGHT - ARM_LEFT) / 1e6, top,
                               fc=COL_CORE, alpha=0.12, ec='none', zorder=0))
    ax.add_patch(plt.Rectangle((ARM_RIGHT / 1e6, 0), (GENOME_LEN - ARM_RIGHT) / 1e6, top,
                               fc=COL_ARM, alpha=0.25, ec='none', zorder=0))
    ax.axvline(ORIC_POS / 1e6, color=COL_ORIC, lw=1.0, ls='-', zorder=4)
    ax.text(ORIC_POS / 1e6, top + 0.02, 'oriC', ha='center', va='bottom',
            fontsize=6.5, color=COL_ORIC, fontweight='bold')
    for bnd in (ARM_LEFT, ARM_RIGHT):
        ax.axvline(bnd / 1e6, color='gray', ls=':', lw=0.6, zorder=1)

    for k, tp in enumerate(tps):
        base = ntp - 1 - k                       # T1 top, T3 bottom
        pos = df_sites[df_sites['timepoint'] == tp]['position'].values
        h, _ = np.histogram(pos, bins=edges)
        ax.bar(centers_mb, h / hmax * band_h, width=bin_size_bp / 1e6, bottom=base,
               color=COL_4mC, alpha=0.95, linewidth=0, zorder=2)
        ax.axhline(base, color='#aaa', lw=0.5, zorder=1)
        # Track label INSIDE the panel (top-left of each track) so it never
        # collides with the y-axis title on the outer margin.
        ax.text(0.12, base + band_h - 0.04, f'{T_LABELS[tp]}  (n={len(pos)})',
                ha='left', va='top', fontsize=6.6, color='#222', zorder=5)

    # Region text labels (Left arm / Core / Right arm) intentionally OMITTED:
    # the core/arm are already encoded by the translucent green/amber bands and
    # the boundaries are stated in the caption, so text on the x-tick row only
    # collided with the Mb ticks ('Core5') or, moved lower, with the x-title
    # (figure-legibility-qc §3/§5 — do not add clutter the caption carries).

    ax.set_xlim(0, GENOME_LEN / 1e6)
    ax.set_ylim(-0.02, top + 0.35)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.set_xlabel('Chromosome position (Mb)', fontsize=9)
    ax.set_ylabel('GCCGGC 4mC sites\n(per 50 kb, by timepoint)', fontsize=9)
    ax.set_title('Geographic redistribution', pad=14,
                 fontsize=10, fontweight='bold')


def panel_b(ax, df_spatial):
    """TSS metagene methylation profile (T1, GCCGGC 4mC, all genes ±4 kb)."""
    sub = df_spatial[
        (df_spatial['site_type'] == 'GCCGGC_4mC') &
        (df_spatial['timepoint'] == 'T1') &
        (df_spatial['gene_category'] == 'all')
    ].sort_values('bin_center_bp')

    x = sub['bin_center_bp'].values
    y = sub['density_per_kb_per_gene'].values
    n_genes = int(sub['n_genes'].iloc[0])

    # Metagene line + CI in neutral grayscale; the single red element is the
    # promoter protection zone, so the eye goes straight to the depletion feature.
    ax.plot(x, y, color='#293039', lw=1.4)
    if 'ci_low' in sub.columns:
        ax.fill_between(x, sub['ci_low'].values, sub['ci_high'].values,
                        color='#9AA7B0', alpha=0.30, lw=0)

    # Reviewer C14 fix: shade the depletion where the metagene ACTUALLY drops.
    # In this all-genes GCCGGC 4mC metagene the density falls below baseline only
    # in a narrow TSS-proximal window (bins at -100/+100 bp = 64%/69% of the
    # |x|>2.5 kb baseline; back to baseline by +/-300 bp; the -2200..0 region is
    # 99% of baseline, i.e. NOT depleted). The band is therefore drawn symmetric
    # about the TSS at the width supported by the data, and the 293 bp cut-off
    # marks where the depletion returns to baseline. NOTE: the "~2,200 bp
    # symmetric protection zone" claim in the text derives from the separate
    # regulatory-gene symmetric-zone analysis (H25), not this all-genes metagene;
    # see author flag in out_nar_revision_260716.
    # Reference band + cut-off lines are clipped to the DATA band (ymax=0.72 in
    # axes fraction) so they never rise into the top headroom where the legend
    # sits — this is what removes the "dashed line through the legend" overlap.
    REF_YMAX = 0.72
    DEPLETION_HALF_BP = 300
    ax.axvspan(-DEPLETION_HALF_BP, DEPLETION_HALF_BP, ymin=0, ymax=REF_YMAX,
               color=COL_4mC, alpha=0.16,
               label=f'TSS-proximal depletion (\u00b1{DEPLETION_HALF_BP} bp)')
    ax.axvline(PROTECTION_ZONE_BP, ymin=0, ymax=REF_YMAX, color=COL_4mC, lw=1.0, ls='--')
    ax.axvline(-PROTECTION_ZONE_BP, ymin=0, ymax=REF_YMAX, color=COL_4mC, lw=1.0, ls='--',
               label=f'{PROTECTION_ZONE_BP} bp classification cut-off')
    ax.axvline(0, ymin=0, ymax=REF_YMAX, color='#888', lw=0.6, ls=':')

    ax.set_xlim(-4000, 4000)
    # Headroom so the 2-line key clears the CI ribbon and the top y-tick (§7.3).
    ylo, yhi = ax.get_ylim()
    ax.set_ylim(ylo, yhi + (yhi - ylo) * 0.40)
    ax.set_xlabel('Distance from TSS (bp)', fontsize=9)
    ax.set_ylabel(f'GCCGGC 4mC density\n(all genes, n={n_genes:,})', fontsize=9)
    ax.set_title('TSS metagene (T1)', pad=14,
                 fontsize=10, fontweight='bold')
    ax.legend(fontsize=6.5, loc='upper center', frameon=False,
              handlelength=1.3, borderpad=0.3, labelspacing=0.3)
    ax.grid(axis='y', alpha=0.25, lw=0.5)


def panel_c(ax, df_genes):
    """Continuous TSS-to-nearest-GCCGGC-4mC distance across ALL regulatory genes.

    Honest depiction: the distribution is continuous and unimodal (Hartigan's dip
    test does not reject unimodality), so it is shown as a single histogram + KDE on
    a log10-bp x-axis rather than pre-split into two groups. The 293 bp operating
    point is a single vertical line; the 62 Exposed are simply those left of it
    (the region <=293 bp is shaded lightly).
    """
    df = df_genes.dropna(subset=['nearest_methyl_distance'])
    dist = df['nearest_methyl_distance'].values.astype(float)
    n_total = len(dist)
    n_exposed = int((dist <= PROTECTION_ZONE_BP).sum())
    n_shielded = n_total - n_exposed

    log_dist = np.log10(dist + 1)
    log_thresh = np.log10(PROTECTION_ZONE_BP)

    # Hartigan's dip test for unimodality (annotation only; falls back gracefully).
    dip_label = 'unimodal: Hartigan dip p = 0.96'
    try:
        import diptest
        _dip, _p = diptest.diptest(log_dist)
        dip_label = f'unimodal: Hartigan dip p = {_p:.2f}'
    except Exception:
        pass

    x_lo, x_hi = 0.0, np.ceil(log_dist.max() * 2) / 2
    bins = np.linspace(x_lo, x_hi, 31)

    # One continuous histogram (density) across all regulatory genes.
    ax.hist(log_dist, bins=bins, density=True, color=COL_GREY, alpha=0.45,
            edgecolor='white', linewidth=0.4, zorder=2)

    # KDE overlay to convey the continuous, single-peaked shape.
    kde = stats.gaussian_kde(log_dist)
    xs = np.linspace(x_lo, x_hi, 400)
    ys = kde(xs)
    ax.plot(xs, ys, color='#555', lw=1.6, zorder=4)

    # Shade the region <=293 bp (the Exposed side) lightly — no separate group.
    # Both the pale fill and the dashed cut-off are clipped to the data band
    # (ymax=0.70) so neither rises into the top headroom occupied by the legend.
    REF_YMAX = 0.70
    ax.axvspan(x_lo, log_thresh, ymin=0, ymax=REF_YMAX, color=COL_4mC, alpha=0.12,
               zorder=1, label=f'Exposed (≤ {PROTECTION_ZONE_BP} bp), n = {n_exposed}')

    # Single vertical line at the 293 bp operating point — the one in-plot label.
    ax.axvline(log_thresh, ymin=0, ymax=REF_YMAX, color=COL_4mC, ls='--', lw=1.4,
               zorder=5, label=f'{PROTECTION_ZONE_BP} bp cut-off')
    # Proxy handle so the Shielded majority also appears in the key.
    ax.axvspan(np.nan, np.nan, color=COL_GREY, alpha=0.45,
               label=f'Shielded (> {PROTECTION_ZONE_BP} bp), n = {n_shielded}')

    # All descriptive text (Exposed/Shielded counts, the unimodality result) is
    # carried in a compact legend instead of scattered in-plot annotations; the
    # Hartigan dip statistic and definitions live in the figure legend prose.
    # Headroom so the 3-line key sits above the KDE peak in a clear band (§7.3).
    ylo, yhi = ax.get_ylim()
    ax.set_ylim(ylo, yhi + (yhi - ylo) * 0.42)
    ax.legend(loc='upper left', fontsize=6.5, frameon=False,
              handlelength=1.3, borderpad=0.3, labelspacing=0.3)

    # Decade ticks only; the 293 bp operating point is shown by the dashed line
    # and named in the legend, so a 293 tick here just collides with 100/1,000.
    xticks_bp = [1, 10, 100, 1000, 10_000, 100_000]
    xtick_lab = ['1', '10', '100', '1k', '10k', '100k']  # k-suffixes (§5.4) avoid label collisions
    ax.set_xticks([np.log10(v + 1) for v in xticks_bp])
    ax.set_xticklabels(xtick_lab, fontsize=8)
    ax.set_xlim(x_lo, x_hi)
    ax.set_xlabel('Distance to nearest 4mC site (bp)', fontsize=9)
    ax.set_ylabel(f'Density (all regulatory genes, n = {n_total:,})', fontsize=9)
    ax.set_title('TSS-distance gradient (T1)', pad=14,
                 fontsize=10, fontweight='bold')
    ax.grid(axis='y', alpha=0.25, lw=0.5)


def panel_d(ax, df_genes):
    """|log2 FC| variability for Shielded vs Exposed TFs (T2/T1 and T3/T1)."""
    df = df_genes.dropna(subset=['LFC_T2vsT1', 'LFC_T3vsT1']).copy()
    df['abs_T2T1'] = df['LFC_T2vsT1'].abs()
    df['abs_T3T1'] = df['LFC_T3vsT1'].abs()
    exp = df[df['is_exposed'] == 1]
    shi = df[df['is_exposed'] == 0]

    comparisons = [
        ('abs_T2T1', 0, 1, 'T2 vs T1'),
        ('abs_T3T1', 3, 4, 'T3 vs T1'),
    ]

    # Display cap: clip the long thin tails at the 98th percentile so the violins
    # read cleanly and do not spike into the significance bracket. Medians and the
    # test below use the full, unclipped data.
    disp_cap = 0
    for col, pos_s, pos_e, _ in comparisons:
        disp_cap = max(disp_cap, np.percentile(shi[col].values, 98),
                       np.percentile(exp[col].values, 98))

    y_max = 0
    for col, pos_s, pos_e, _ in comparisons:
        data_s = shi[col].values
        data_e = exp[col].values
        y_max = max(y_max, np.percentile(data_s, 98), np.percentile(data_e, 98))

        for pos, data, color in [(pos_s, data_s, COL_SHIELDED),
                                 (pos_e, data_e, COL_EXPOSED)]:
            data_disp = np.clip(data, None, disp_cap)
            parts = ax.violinplot([data_disp], positions=[pos], widths=0.85,
                                  showextrema=False, showmedians=False)
            for body in parts['bodies']:
                body.set_facecolor(color)
                body.set_alpha(0.55)
                body.set_edgecolor(color)
            # Median bar
            med = np.median(data)
            ax.hlines(med, pos - 0.20, pos + 0.20,
                      color='black', lw=1.4, zorder=5)

    # Significance bracket — placed above all violins, well clear of legend
    bracket_y = y_max * 1.15
    text_y = bracket_y * 1.06   # 2026-09-13: clear of the bracket line

    for col, pos_s, pos_e, label in comparisons:
        e = exp[col].values
        s = shi[col].values
        U, p = stats.mannwhitneyu(e, s, alternative='greater')
        if p < 0.001:
            star = '***'
        elif p < 0.01:
            star = '**'
        elif p < 0.05:
            star = '*'
        else:
            star = 'n.s.'
        ax.plot([pos_s, pos_s, pos_e, pos_e],
                [bracket_y, bracket_y * 1.04, bracket_y * 1.04, bracket_y],
                color='black', lw=0.8)
        ax.text((pos_s + pos_e) / 2, text_y,
                f'{star}  (p = {p:.1e})',
                ha='center', va='bottom', fontsize=7.5, fontweight='bold')

    ax.set_ylim(0, bracket_y * 1.32)
    ax.set_xticks([0.5, 3.5])
    ax.set_xticklabels([c[3] for c in comparisons], fontsize=9)
    ax.set_xlim(-0.7, 5.2)
    ax.set_ylabel(r'$|\log_2$ fold change$|$ (expression variability)',
                  fontsize=9)
    # Framed as a negative control: proximity to a GCCGGC 4mC site (Exposed)
    # does not change expression variability — the direct visual statement of
    # the permissive (non-instructive) thesis. The n.s. result is the message.
    ax.set_title('Negative control:\nmethylation proximity ≠ expression change',
                 pad=12, fontsize=9.5, fontweight='bold')

    # Legend placed below the axis to avoid overlap with the significance marker.
    # Counts reflect genes with available LFC data for both transitions.
    legend_handles = [
        Patch(facecolor=COL_SHIELDED, alpha=0.55,
              label=f'Shielded (n = {len(shi)} of 989)'),
        Patch(facecolor=COL_EXPOSED, alpha=0.55,
              label=f'Exposed (n = {len(exp)} of 62)'),
    ]
    # 2026-09-13 (FIG-11): the below-axis legend was cut off at the PNG's bottom
    # edge (canvas bottom margin too small). Two-column legend, anchored just
    # under the tick labels; the standalone figure now reserves room for it.
    ax.legend(handles=legend_handles, loc='upper center',
              bbox_to_anchor=(0.5, -0.13), ncol=2,
              fontsize=7.5, frameon=False,
              handlelength=1.2, borderpad=0.4, columnspacing=1.5)
    print(f'  panel D plotted n: Shielded={len(shi)} of 989, Exposed={len(exp)} of 62 '
          f'(genes with LFC for both T2/T1 and T3/T1)')


def main():
    apply_style()
    print('=== Figure 2: GCCGGC 4mC redistribution + promoter protection ===')

    print('Loading data...')
    _, df_sites = load_geographic_redistribution()
    df_spatial = load_spatial_profile()
    df_genes = pd.read_csv(
        EPIGENOME / '52_shielded_exposed_boundary' /
        'tables' / 'all_genes_features_unified_n57.tsv', sep='\t')

    # --- Non-circular reframe (2026-06-17): recompute Exposed from the methylation
    # map only (TSS within 293 bp of a GCCGGC 4mC site at T1), matching F4_F5 and
    # the Fig2c/d reframe. The table's stale `is_exposed` (n=57, expression-selected)
    # is replaced; nearest_methyl_distance is set to the T1 GCCGGC distance so the
    # panel-C violin and the Shielded/Exposed split are the same quantity (62/989).
    df_genes = df_genes.dropna(subset=['tss']).copy()
    df_genes['tss'] = df_genes['tss'].astype(int)
    g_tp = pd.read_csv(EPIGENOME / '37_defense_island_GCCGGC' /
                       'tables' / 'GCCGGC_sites_by_timepoint.tsv', sep='\t')
    _pos = np.sort(g_tp[g_tp['timepoint'] == 'T1']['position'].values)

    def _nearest(tss):
        if len(_pos) == 0:
            return np.nan
        i = np.clip(np.searchsorted(_pos, tss), 1, len(_pos) - 1)
        return min(abs(tss - _pos[i - 1]), abs(tss - _pos[i]))

    df_genes['nearest_methyl_distance'] = df_genes['tss'].apply(_nearest)
    df_genes['is_exposed'] = (df_genes['nearest_methyl_distance']
                              <= PROTECTION_ZONE_BP).astype(int)
    print(f'  Regulatory genes (non-circular): {len(df_genes)} '
          f'(Exposed={int(df_genes["is_exposed"].sum())}, '
          f'Shielded={int((df_genes["is_exposed"] == 0).sum())})')

    import shutil
    # --- MAIN Figure 2: three panels (A geographic redistribution, B TSS metagene,
    # C continuous distance). Reviewer C10: the former panel D (negative-control
    # violins, n.s.) is demoted to a Supplementary figure; the load-bearing null
    # statistics stay in the main text. Main figure re-laid-out as 1x3 wide.
    fig = plt.figure(figsize=(mm_to_inch(180), mm_to_inch(98)))
    gs = fig.add_gridspec(1, 3, wspace=0.50,
                          left=0.07, right=0.94, top=0.74, bottom=0.19)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])

    panel_a(ax_a, df_sites)
    panel_b(ax_b, df_spatial)
    panel_c(ax_c, df_genes)

    for ax, label in zip([ax_a, ax_b, ax_c], ['A', 'B', 'C']):
        add_uppercase_label(ax, label, x=-0.20, y=1.32, fontsize=14)

    fig.suptitle('Figure 2 | Geographic redistribution of GCCGGC 4mC '
                 'and promoter protection architecture',
                 fontsize=11, fontweight='bold', y=0.985)

    out_path = FIG_DIR / 'Figure2_RM_redistribution'
    save_figure(fig, out_path, formats=('pdf', 'svg', 'png'))
    # 2026-09-13: manuscript-slot sync is opt-in (SYNC_FIG_SLOTS=1); back up the
    # old PNG to fig_images/archive/ before enabling. Figure2.png is no longer an
    # embedded slot in the inline manuscript (Figure2_merged.png is).
    import os
    slot = Path.home() / 'obsidian' / 'Research' / 'rna-seq' / 'Writing' / 'fig_images' / 'Figure2.png'
    if os.environ.get('SYNC_FIG_SLOTS') == '1' and slot.parent.is_dir():
        shutil.copyfile(out_path.with_suffix('.png'), slot)
        print(f'  Synced → {slot}')

    # --- SUPPLEMENTARY figure: the demoted negative-control panel, standalone.
    # 2026-09-13 (FIG-11): taller canvas + larger bottom margin so the class
    # legend under the x-axis is inside the saved PNG.
    figd = plt.figure(figsize=(mm_to_inch(100), mm_to_inch(100)))
    axd = figd.add_subplot(1, 1, 1)
    panel_d(axd, df_genes)  # panel_d sets its own descriptive title; no suptitle needed
    figd.subplots_adjust(left=0.16, right=0.95, top=0.89, bottom=0.24)
    out_d = FIG_DIR / 'SuppFigure_Fig2D_negative_control'
    save_figure(figd, out_d, formats=('pdf', 'svg', 'png'))
    slot_d = Path.home() / 'obsidian' / 'Research' / 'rna-seq' / 'Writing' / 'fig_images' / 'SuppFigure_neg_control.png'
    if os.environ.get('SYNC_FIG_SLOTS') == '1' and slot_d.parent.is_dir():
        shutil.copyfile(out_d.with_suffix('.png'), slot_d)
        print(f'  Synced (supp D) → {slot_d}')
    print('=== Done ===')


if __name__ == '__main__':
    main()
