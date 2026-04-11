#!/usr/bin/env python3
"""
H25 Improved Figures: Smoothed profiles and refined protection zone analysis.
Applies Gaussian smoothing to reveal underlying trends while preserving the
raw data in CI ribbons.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.ndimage import gaussian_filter1d
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

BASE = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/48_TSS_methylation_gradient")
FIG_DIR = BASE / "figures"
TBL_DIR = BASE / "tables"

# Load pre-computed spatial profile data
profile_df = pd.read_csv(TBL_DIR / "spatial_profile_data.tsv", sep='\t')
region_df = pd.read_csv(TBL_DIR / "region_stratified_profiles.tsv", sep='\t')
sub_df = pd.read_csv(TBL_DIR / "subcategory_profiles.tsv", sep='\t')
perm_df = pd.read_csv(TBL_DIR / "permutation_test_results.tsv", sep='\t')
prot_df = pd.read_csv(TBL_DIR / "protection_zone_metrics.tsv", sep='\t')

BIN_SIZE = 200
SIGMA = 2.0  # Gaussian smoothing sigma in bins (= 400 bp effective)

# Color scheme
COL_REG = '#D32F2F'
COL_NONREG = '#1565C0'
COL_REG_LIGHT = '#FFCDD2'
COL_NONREG_LIGHT = '#BBDEFB'

SUBTYPE_COLORS = {
    'sigma_factor': '#E65100',
    'TCS_response_reg': '#AD1457',
    'SARP': '#6A1B9A',
    'sensor_kinase': '#00695C',
    'other_TF': '#F57F17',
    'non_regulatory': '#1565C0',
}
SUBTYPE_LABELS = {
    'sigma_factor': 'Sigma factors',
    'TCS_response_reg': 'TCS response regulators',
    'SARP': 'SARP family',
    'sensor_kinase': 'Sensor kinases',
    'other_TF': 'Other TFs',
    'non_regulatory': 'Non-regulatory',
}


def get_profile(df, site_type, timepoint, category):
    """Extract profile for given parameters."""
    mask = (df['site_type'] == site_type) & (df['timepoint'] == timepoint) & (df['gene_category'] == category)
    sub = df[mask].sort_values('bin_center_bp')
    return sub['bin_center_bp'].values, sub['density_per_kb_per_gene'].values, sub['ci_low'].values, sub['ci_high'].values, sub['n_genes'].values[0]


def smooth(y, sigma=SIGMA):
    """Apply Gaussian smoothing."""
    return gaussian_filter1d(y, sigma=sigma)


def plot_smoothed_profile(ax, bin_centers, reg_density, reg_ci_low, reg_ci_high, n_reg,
                          nonreg_density, nonreg_ci_low, nonreg_ci_high, n_nonreg,
                          title, show_legend=True, show_annotations=True):
    """Plot smoothed regulatory vs non-regulatory profile with raw CI."""
    reg_sm = smooth(reg_density)
    nonreg_sm = smooth(nonreg_density)

    # CI ribbons (raw)
    ax.fill_between(bin_centers, reg_ci_low, reg_ci_high,
                    color=COL_REG, alpha=0.1, linewidth=0)
    ax.fill_between(bin_centers, nonreg_ci_low, nonreg_ci_high,
                    color=COL_NONREG, alpha=0.1, linewidth=0)

    # Raw data as thin transparent lines
    ax.plot(bin_centers, nonreg_density, color=COL_NONREG, linewidth=0.5, alpha=0.3)
    ax.plot(bin_centers, reg_density, color=COL_REG, linewidth=0.5, alpha=0.3)

    # Smoothed lines
    ax.plot(bin_centers, nonreg_sm, color=COL_NONREG, linewidth=2.0,
            label=f"Non-regulatory (n={n_nonreg})")
    ax.plot(bin_centers, reg_sm, color=COL_REG, linewidth=2.0,
            label=f"Regulatory (n={n_reg})")

    ax.axvline(0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.set_xlabel('Distance from TSS (bp)', fontsize=10)
    ax.set_ylabel('Methylation density\n(sites/kb/gene)', fontsize=10)
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.tick_params(labelsize=9)

    if show_legend:
        ax.legend(fontsize=8, loc='upper right', framealpha=0.9)

    if show_annotations:
        ylims = ax.get_ylim()
        ax.text(-4200, ylims[1] - 0.02 * (ylims[1] - ylims[0]), 'Upstream',
                fontsize=8, ha='center', va='top', color='gray', style='italic')
        ax.text(4200, ylims[1] - 0.02 * (ylims[1] - ylims[0]), 'Gene body',
                fontsize=8, ha='center', va='top', color='gray', style='italic')


# ============================================================
# Figure 1: Main profile (improved)
# ============================================================
print("Creating improved main profile...")
bc, reg_d, reg_lo, reg_hi, n_reg = get_profile(profile_df, 'All_methylation', 'T1', 'regulatory')
_, nonreg_d, nonreg_lo, nonreg_hi, n_nonreg = get_profile(profile_df, 'All_methylation', 'T1', 'non_regulatory')

fig1, ax1 = plt.subplots(figsize=(8, 5))
plot_smoothed_profile(ax1, bc, reg_d, reg_lo, reg_hi, n_reg,
                      nonreg_d, nonreg_lo, nonreg_hi, n_nonreg,
                      'Methylation density around TSS (All sites, T1)')

# Add shaded protection zone
ax1.axvspan(-1300, 700, alpha=0.05, color='red', zorder=0)
ax1.annotate('Protection zone\n(-1300 to +700 bp)', xy=(-300, ax1.get_ylim()[0]),
             fontsize=7, ha='center', va='bottom', color='darkred', alpha=0.7)

fig1.tight_layout()
fig1.savefig(FIG_DIR / "TSS_methylation_profile.pdf", dpi=300, bbox_inches='tight')
fig1.savefig(FIG_DIR / "TSS_methylation_profile.svg", dpi=300, bbox_inches='tight')
print("  Saved TSS_methylation_profile.pdf/svg")

# ============================================================
# Figure 2: Motif-specific profiles (improved)
# ============================================================
print("Creating improved motif-specific profiles...")
fig2, axes2 = plt.subplots(2, 2, figsize=(12, 9))
motifs = [
    ('All_4mC', 'All 4mC sites (T1)'),
    ('All_6mA', 'All 6mA sites (T1)'),
    ('GCCGGC_4mC', 'GCCGGC 4mC sites (T1)'),
    ('AAGCCCG_6mA', 'AAGCCCG 6mA sites (T1)'),
]
for ax, (sn, title) in zip(axes2.flat, motifs):
    bc, rd, rlo, rhi, nr = get_profile(profile_df, sn, 'T1', 'regulatory')
    _, nd, nlo, nhi, nn = get_profile(profile_df, sn, 'T1', 'non_regulatory')
    plot_smoothed_profile(ax, bc, rd, rlo, rhi, nr, nd, nlo, nhi, nn, title,
                          show_annotations=False)

fig2.suptitle('Motif-specific methylation profiles around TSS', fontsize=13, fontweight='bold', y=1.01)
fig2.tight_layout()
fig2.savefig(FIG_DIR / "motif_specific_profiles.pdf", dpi=300, bbox_inches='tight')
fig2.savefig(FIG_DIR / "motif_specific_profiles.svg", dpi=300, bbox_inches='tight')
print("  Saved motif_specific_profiles.pdf/svg")

# ============================================================
# Figure 3: Heatmap (improved with smoothing)
# ============================================================
print("Creating improved heatmap...")
site_names_ordered = ['All_methylation', 'All_4mC', 'All_6mA', 'GCCGGC_4mC', 'AAGCCCG_6mA']
site_labels = ['All methylation', 'All 4mC', 'All 6mA', 'GCCGGC 4mC', 'AAGCCCG 6mA']
n_bins = len(bc)

ratio_matrix = np.zeros((len(site_names_ordered), n_bins))
for i, sn in enumerate(site_names_ordered):
    _, rd, _, _, _ = get_profile(profile_df, sn, 'T1', 'regulatory')
    _, nd, _, _, _ = get_profile(profile_df, sn, 'T1', 'non_regulatory')
    rd_sm = smooth(rd)
    nd_sm = smooth(nd)
    ratio_matrix[i] = np.where(nd_sm > 0, rd_sm / nd_sm, np.nan)

fig3, ax3 = plt.subplots(figsize=(10, 4))
ratio_clipped = np.clip(ratio_matrix, 0.6, 1.4)
im = ax3.imshow(ratio_clipped, aspect='auto', cmap='RdBu_r',
                vmin=0.6, vmax=1.4,
                extent=[bc[0] - BIN_SIZE/2, bc[-1] + BIN_SIZE/2,
                        len(site_names_ordered) - 0.5, -0.5])
ax3.set_yticks(range(len(site_labels)))
ax3.set_yticklabels(site_labels, fontsize=10)
ax3.set_xlabel('Distance from TSS (bp)', fontsize=10)
ax3.set_title('Regulatory / Non-regulatory density ratio (smoothed)', fontsize=11, fontweight='bold')
ax3.axvline(0, color='black', linestyle='--', linewidth=1.0, alpha=0.7)
cbar = plt.colorbar(im, ax=ax3, label='Density ratio (reg/nonreg)', shrink=0.8)
cbar.ax.axhline(y=1.0, color='black', linewidth=0.8, linestyle='-')
# Add significance dots from permutation tests
for i, sn in enumerate(site_names_ordered):
    sub = perm_df[(perm_df['site_type'] == sn) & (perm_df['timepoint'] == 'T1')]
    sig_bins = sub[sub['significant']]['bin_center'].values
    for sb in sig_bins:
        ax3.plot(sb, i, '*', color='black', markersize=6)

fig3.tight_layout()
fig3.savefig(FIG_DIR / "protection_zone_heatmap.pdf", dpi=300, bbox_inches='tight')
fig3.savefig(FIG_DIR / "protection_zone_heatmap.svg", dpi=300, bbox_inches='tight')
print("  Saved protection_zone_heatmap.pdf/svg")

# ============================================================
# Figure 4: Temporal comparison (improved)
# ============================================================
print("Creating improved temporal comparison...")
fig4, axes4 = plt.subplots(1, 2, figsize=(12, 5))

for ax, (cat, color) in zip(axes4, [('regulatory', COL_REG), ('non_regulatory', COL_NONREG)]):
    bc, t1_d, t1_lo, t1_hi, n1 = get_profile(profile_df, 'All_methylation', 'T1', cat)
    _, t2_d, t2_lo, t2_hi, n2 = get_profile(profile_df, 'All_methylation', 'T2', cat)

    t1_sm = smooth(t1_d)
    t2_sm = smooth(t2_d)

    ax.fill_between(bc, t1_lo, t1_hi, alpha=0.1, color=color)
    ax.plot(bc, t1_d, color=color, linewidth=0.5, alpha=0.3)
    ax.plot(bc, t1_sm, color=color, linewidth=2.0, label=f'T1 (24h, n={n1})')

    ax.plot(bc, t2_d, color=color, linewidth=0.5, alpha=0.2, linestyle='--')
    ax.plot(bc, t2_sm, color=color, linewidth=2.0, linestyle='--', alpha=0.8,
            label=f'T2 (36h, n={n2})')

    ax.axvline(0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.set_xlabel('Distance from TSS (bp)', fontsize=10)
    ax.set_ylabel('Methylation density\n(sites/kb/gene)', fontsize=10)
    ax.set_title(f'{"Regulatory" if cat == "regulatory" else "Non-regulatory"} genes',
                 fontsize=11, fontweight='bold')
    ax.legend(fontsize=9)

fig4.suptitle('Temporal comparison: T1 vs T2 (smoothed)', fontsize=13, fontweight='bold', y=1.01)
fig4.tight_layout()
fig4.savefig(FIG_DIR / "temporal_comparison.pdf", dpi=300, bbox_inches='tight')
fig4.savefig(FIG_DIR / "temporal_comparison.svg", dpi=300, bbox_inches='tight')
print("  Saved temporal_comparison.pdf/svg")

# ============================================================
# Figure 5: Core/arm (improved)
# ============================================================
print("Creating improved core/arm profiles...")
fig5, axes5 = plt.subplots(1, 2, figsize=(12, 5))

for ax, region in zip(axes5, ['core', 'arm']):
    for reg_status, color, label_prefix in [
        ('non_regulatory', COL_NONREG, 'Non-reg'),
        ('regulatory', COL_REG, 'Regulatory')
    ]:
        cat_name = f"{region}_{reg_status}"
        sub = region_df[region_df['region_category'] == cat_name].sort_values('bin_center_bp')
        if len(sub) == 0:
            continue
        x = sub['bin_center_bp'].values
        y = sub['density_per_kb_per_gene'].values
        lo = sub['ci_low'].values
        hi = sub['ci_high'].values
        ng = sub['n_genes'].values[0]

        y_sm = smooth(y)
        ax.fill_between(x, lo, hi, color=color, alpha=0.1)
        ax.plot(x, y, color=color, linewidth=0.5, alpha=0.3)
        ax.plot(x, y_sm, color=color, linewidth=2.0, label=f"{label_prefix} (n={ng})")

    ax.axvline(0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.set_xlabel('Distance from TSS (bp)', fontsize=10)
    ax.set_ylabel('Methylation density\n(sites/kb/gene)', fontsize=10)
    ax.set_title(f'{region.capitalize()} region', fontsize=11, fontweight='bold')
    ax.legend(fontsize=8)

fig5.suptitle('Core vs Arm: TSS methylation profiles (All sites, T1)', fontsize=13, fontweight='bold', y=1.01)
fig5.tight_layout()
fig5.savefig(FIG_DIR / "core_arm_profiles.pdf", dpi=300, bbox_inches='tight')
fig5.savefig(FIG_DIR / "core_arm_profiles.svg", dpi=300, bbox_inches='tight')
print("  Saved core_arm_profiles.pdf/svg")

# ============================================================
# Figure 6: Sub-category profiles (improved, exclude SARP n=7)
# ============================================================
print("Creating improved subcategory profiles...")
fig6, ax6 = plt.subplots(figsize=(9, 6))

# Plot non-regulatory first as thick reference line
for cat_name in ['non_regulatory', 'sigma_factor', 'TCS_response_reg', 'sensor_kinase', 'other_TF']:
    sub = sub_df[sub_df['reg_subtype'] == cat_name].sort_values('bin_center_bp')
    if len(sub) == 0:
        continue
    x = sub['bin_center_bp'].values
    y = sub['density_per_kb_per_gene'].values
    ng = sub['n_genes'].values[0]
    y_sm = smooth(y)

    color = SUBTYPE_COLORS[cat_name]
    label = SUBTYPE_LABELS[cat_name]
    lw = 2.5 if cat_name == 'non_regulatory' else 1.5
    alpha = 1.0 if cat_name == 'non_regulatory' else 0.85

    ax6.plot(x, y, color=color, linewidth=0.3, alpha=0.2)
    ax6.plot(x, y_sm, color=color, linewidth=lw, alpha=alpha,
             label=f"{label} (n={ng})")

ax6.axvline(0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
ax6.set_xlabel('Distance from TSS (bp)', fontsize=10)
ax6.set_ylabel('Methylation density (sites/kb/gene)', fontsize=10)
ax6.set_title('Regulatory sub-type methylation profiles (All sites, T1)\n(SARP excluded, n=7)',
              fontsize=11, fontweight='bold')
ax6.legend(fontsize=8, loc='upper right')

fig6.tight_layout()
fig6.savefig(FIG_DIR / "subcategory_profiles.pdf", dpi=300, bbox_inches='tight')
fig6.savefig(FIG_DIR / "subcategory_profiles.svg", dpi=300, bbox_inches='tight')
print("  Saved subcategory_profiles.pdf/svg")

# ============================================================
# Comprehensive multi-panel figure (publication quality)
# ============================================================
print("Creating comprehensive multi-panel figure...")

fig = plt.figure(figsize=(16, 20))
gs = gridspec.GridSpec(5, 2, figure=fig, hspace=0.40, wspace=0.30,
                       height_ratios=[1.2, 1, 1, 0.8, 1])

# --- Panel A: Main profile (full width) ---
ax_A = fig.add_subplot(gs[0, :])
bc_r, rd, rlo, rhi, nr = get_profile(profile_df, 'All_methylation', 'T1', 'regulatory')
bc_n, nd, nlo, nhi, nn = get_profile(profile_df, 'All_methylation', 'T1', 'non_regulatory')

rd_sm = smooth(rd)
nd_sm = smooth(nd)

ax_A.fill_between(bc_r, rlo, rhi, color=COL_REG, alpha=0.10, linewidth=0)
ax_A.fill_between(bc_n, nlo, nhi, color=COL_NONREG, alpha=0.10, linewidth=0)
ax_A.plot(bc_n, nd, color=COL_NONREG, linewidth=0.5, alpha=0.3)
ax_A.plot(bc_r, rd, color=COL_REG, linewidth=0.5, alpha=0.3)
ax_A.plot(bc_n, nd_sm, color=COL_NONREG, linewidth=2.5,
          label=f"Non-regulatory (n={nn})")
ax_A.plot(bc_r, rd_sm, color=COL_REG, linewidth=2.5,
          label=f"Regulatory (n={nr})")
ax_A.axvline(0, color='black', linestyle='--', linewidth=1.0, alpha=0.5)
ax_A.axvspan(-1300, 700, alpha=0.04, color='red', zorder=0)
ax_A.set_xlabel('Distance from TSS (bp)', fontsize=11)
ax_A.set_ylabel('Methylation density\n(sites/kb/gene)', fontsize=11)
ax_A.set_title('A. TSS-centered methylation density (All sites, T1)',
               fontsize=12, fontweight='bold', loc='left')
ax_A.legend(fontsize=10, loc='upper right', framealpha=0.9)
ylims_A = ax_A.get_ylim()
ax_A.text(-4200, ylims_A[1] - 0.03 * (ylims_A[1] - ylims_A[0]), 'Upstream',
          fontsize=9, ha='center', va='top', color='gray', style='italic')
ax_A.text(4200, ylims_A[1] - 0.03 * (ylims_A[1] - ylims_A[0]), 'Gene body',
          fontsize=9, ha='center', va='top', color='gray', style='italic')

# --- Panels B-E: Motif-specific ---
motifs_panel = [
    ('All_4mC', 'B. All 4mC sites (T1)'),
    ('All_6mA', 'C. All 6mA sites (T1)'),
    ('GCCGGC_4mC', 'D. GCCGGC 4mC sites (T1)'),
    ('AAGCCCG_6mA', 'E. AAGCCCG 6mA sites (T1)'),
]
for idx, (sn, title) in enumerate(motifs_panel):
    row = 1 + idx // 2
    col = idx % 2
    ax = fig.add_subplot(gs[row, col])
    bc_, rd_, rlo_, rhi_, nr_ = get_profile(profile_df, sn, 'T1', 'regulatory')
    _, nd_, nlo_, nhi_, nn_ = get_profile(profile_df, sn, 'T1', 'non_regulatory')

    rd_sm_ = smooth(rd_)
    nd_sm_ = smooth(nd_)

    ax.fill_between(bc_, rlo_, rhi_, color=COL_REG, alpha=0.10, linewidth=0)
    ax.fill_between(bc_, nlo_, nhi_, color=COL_NONREG, alpha=0.10, linewidth=0)
    ax.plot(bc_, nd_, color=COL_NONREG, linewidth=0.5, alpha=0.3)
    ax.plot(bc_, rd_, color=COL_REG, linewidth=0.5, alpha=0.3)
    ax.plot(bc_, nd_sm_, color=COL_NONREG, linewidth=1.8,
            label=f"Non-reg (n={nn_})")
    ax.plot(bc_, rd_sm_, color=COL_REG, linewidth=1.8,
            label=f"Regulatory (n={nr_})")
    ax.axvline(0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.set_xlabel('Distance from TSS (bp)', fontsize=10)
    ax.set_ylabel('Density (sites/kb/gene)', fontsize=10)
    ax.set_title(title, fontsize=11, fontweight='bold', loc='left')
    ax.legend(fontsize=8, loc='upper right')
    ax.tick_params(labelsize=9)

# --- Panel F: Heatmap (smoothed ratio) ---
ax_F = fig.add_subplot(gs[3, 0])
ratio_sm = np.zeros_like(ratio_matrix)
for i in range(ratio_matrix.shape[0]):
    ratio_sm[i] = smooth(ratio_matrix[i])
ratio_clipped2 = np.clip(ratio_sm, 0.6, 1.4)
im2 = ax_F.imshow(ratio_clipped2, aspect='auto', cmap='RdBu_r',
                   vmin=0.6, vmax=1.4,
                   extent=[bc[0] - BIN_SIZE/2, bc[-1] + BIN_SIZE/2,
                           len(site_names_ordered) - 0.5, -0.5])
ax_F.set_yticks(range(len(site_labels)))
ax_F.set_yticklabels(site_labels, fontsize=9)
ax_F.set_xlabel('Distance from TSS (bp)', fontsize=10)
ax_F.set_title('F. Density ratio (reg/non-reg, smoothed)', fontsize=11, fontweight='bold', loc='left')
ax_F.axvline(0, color='black', linestyle='--', linewidth=1.0, alpha=0.7)
cbar2 = plt.colorbar(im2, ax=ax_F, label='Ratio', shrink=0.8)
cbar2.ax.axhline(y=1.0, color='black', linewidth=0.8)

# --- Panel G: Protection zone bar chart ---
ax_G = fig.add_subplot(gs[3, 1])
pz_widths = prot_df['protection_zone_width_bp'].values
pz_labels_list = site_labels
colors_bar = ['#424242', '#1565C0', '#D32F2F', '#00695C', '#F57F17']
bars = ax_G.barh(range(len(pz_labels_list)), pz_widths, color=colors_bar,
                 edgecolor='black', linewidth=0.5)
ax_G.set_yticks(range(len(pz_labels_list)))
ax_G.set_yticklabels(pz_labels_list, fontsize=9)
ax_G.set_xlabel('Protection zone width (bp)', fontsize=10)
ax_G.set_title('G. Protection zone width by motif', fontsize=11, fontweight='bold', loc='left')
for i, w in enumerate(pz_widths):
    if w > 0:
        ax_G.text(w + 30, i, f'{w:.0f} bp', va='center', fontsize=9)
    else:
        ax_G.text(30, i, 'n.d.', va='center', fontsize=9, color='gray')

# --- Panel H: Core vs Arm comparison ---
ax_H1 = fig.add_subplot(gs[4, 0])
ax_H2 = fig.add_subplot(gs[4, 1])

for ax, region, panel_lbl in [(ax_H1, 'core', 'H'), (ax_H2, 'arm', 'I')]:
    for reg_status, color, label_prefix in [
        ('non_regulatory', COL_NONREG, 'Non-reg'),
        ('regulatory', COL_REG, 'Regulatory')
    ]:
        cat_name = f"{region}_{reg_status}"
        sub = region_df[region_df['region_category'] == cat_name].sort_values('bin_center_bp')
        if len(sub) == 0:
            continue
        x = sub['bin_center_bp'].values
        y = sub['density_per_kb_per_gene'].values
        lo = sub['ci_low'].values
        hi = sub['ci_high'].values
        ng = sub['n_genes'].values[0]
        y_sm = smooth(y)

        ax.fill_between(x, lo, hi, color=color, alpha=0.1)
        ax.plot(x, y, color=color, linewidth=0.5, alpha=0.3)
        ax.plot(x, y_sm, color=color, linewidth=2.0, label=f"{label_prefix} (n={ng})")

    ax.axvline(0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.set_xlabel('Distance from TSS (bp)', fontsize=10)
    ax.set_ylabel('Density (sites/kb/gene)', fontsize=10)
    ax.set_title(f'{panel_lbl}. {region.capitalize()} region (All sites, T1)',
                 fontsize=11, fontweight='bold', loc='left')
    ax.legend(fontsize=8)

fig.savefig(FIG_DIR / "H25_comprehensive_summary.pdf", dpi=300, bbox_inches='tight')
fig.savefig(FIG_DIR / "H25_comprehensive_summary.svg", dpi=300, bbox_inches='tight')
print("  Saved H25_comprehensive_summary.pdf/svg")

plt.close('all')

# ============================================================
# Quantitative summary for report
# ============================================================
print("\n" + "=" * 60)
print("KEY METRICS FOR REPORT")
print("=" * 60)

# Overall depletion at TSS for regulatory vs non-regulatory
# Use ±500bp window
tss_mask = np.abs(bc) <= 500
flank_mask = np.abs(bc) > 3000

for sn in site_names_ordered:
    _, rd, _, _, _ = get_profile(profile_df, sn, 'T1', 'regulatory')
    _, nd, _, _, _ = get_profile(profile_df, sn, 'T1', 'non_regulatory')
    _, ad, _, _, _ = get_profile(profile_df, sn, 'T1', 'all')

    reg_tss = rd[tss_mask].mean()
    nonreg_tss = nd[tss_mask].mean()
    all_tss = ad[tss_mask].mean()
    reg_flank = rd[flank_mask].mean()
    nonreg_flank = nd[flank_mask].mean()
    all_flank = ad[flank_mask].mean()

    print(f"\n{sn}:")
    print(f"  TSS (±500bp) density: reg={reg_tss:.4f}, nonreg={nonreg_tss:.4f}, all={all_tss:.4f}")
    print(f"  Flanking (>3kb) density: reg={reg_flank:.4f}, nonreg={nonreg_flank:.4f}, all={all_flank:.4f}")
    if nonreg_tss > 0:
        print(f"  Reg/nonreg ratio at TSS: {reg_tss/nonreg_tss:.3f}")
    if all_flank > 0:
        print(f"  Universal TSS depletion: {(1 - all_tss/all_flank)*100:.1f}%")
    if reg_flank > 0:
        print(f"  Regulatory TSS depletion: {(1 - reg_tss/reg_flank)*100:.1f}%")
    if nonreg_flank > 0:
        print(f"  Non-regulatory TSS depletion: {(1 - nonreg_tss/nonreg_flank)*100:.1f}%")

# Core vs arm summary
print("\n--- Core vs Arm ---")
for region in ['core', 'arm']:
    for reg_status in ['regulatory', 'non_regulatory']:
        cat_name = f"{region}_{reg_status}"
        sub = region_df[region_df['region_category'] == cat_name].sort_values('bin_center_bp')
        if len(sub) == 0:
            continue
        y = sub['density_per_kb_per_gene'].values
        tss_val = y[tss_mask].mean()
        flank_val = y[flank_mask].mean()
        depl = (1 - tss_val/flank_val)*100 if flank_val > 0 else float('nan')
        print(f"  {cat_name}: TSS={tss_val:.4f}, flank={flank_val:.4f}, depletion={depl:.1f}%")

print("\nDone!")
