#!/usr/bin/env python3
"""
SARP Binding Site – Methylation Spatial Integration Analysis

SARPs bind direct repeats of 7-mer (TCGAGC(G/C)) at/near -35 box,
replacing σHrdB R4 direct DNA contact. This script tests whether
methylation at SARP binding zones modulates BGC gene expression.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import re
from scipy import stats
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
import warnings
warnings.filterwarnings('ignore')

# ── Paths ──
TSS_TABLE = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/18_tss_analyses/comprehensive_tss_table.csv")
METHYL_SITES = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/01_integration/high_confidence_sites_weighted.csv")
INTEGRATED_CSV = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/01_integration/integrated_methyl_expression_weighted.csv")
BGC_DEF = Path("/Users/okaban/bioinfo/rna-seq/05_annotation/analysis/05_annotation_260128_v1/tables/BGC_definition_manual.tsv")
FIMO_RESULTS = Path("/Users/okaban/bioinfo/rna-seq/10_SARP_motif_scan/analysis/10_SARP_motif_scan_260128_v1/tables/fimo_BGC_promoters.tsv")
GENOME_FASTA = "/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/GCF_000203835.1_ASM20383v1_genomic.fna"
OUTPUT_DIR = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/18_tss_analyses")

plt.rcParams.update({
    'font.size': 10, 'axes.labelsize': 12, 'axes.titlesize': 13,
    'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
    'font.family': 'sans-serif',
})


def load_genome(fasta_path):
    seq_parts = []
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith('>'):
                seq_parts.append(line.strip())
    return ''.join(seq_parts).upper()


def revcomp(seq):
    return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def scan_sarp_motifs(genome_seq, genes_df, window_up=100, window_down=10):
    """
    Scan for SARP binding motifs in the -35 box vicinity.
    Motifs: TCGAGC[GC] (ActII-ORF4 consensus) and its variants.
    Also scan for TCGA-containing heptamers and direct repeat patterns.
    """
    # SARP core motif patterns
    patterns = {
        'TCGAGC_strict': 'TCGAGC[GC]',       # ActII-ORF4 strict
        'TCGA_heptamer':  'TCGA...',           # TCGA-containing 7-mer (any)
        'GCTCGAA':        'GCTCGAA',           # OTC consensus (Arias et al.)
    }

    records = []
    for _, gene in genes_df.iterrows():
        tss = gene['tss']
        strand = gene['strand']
        gid = gene['gene_id']

        # Focus on -35 box vicinity: -100 to +10 from TSS
        if strand == '+':
            start = max(0, tss - window_up - 1)
            end = min(len(genome_seq), tss + window_down)
            seq = genome_seq[start:end]
            for pname, pattern in patterns.items():
                for m in re.finditer(pattern, seq):
                    rel = (start + m.start()) - (tss - 1)
                    records.append({
                        'gene_id': gid, 'strand': strand,
                        'pattern': pname, 'match_seq': m.group(),
                        'abs_pos': start + m.start() + 1,
                        'rel_pos': rel, 'match_strand': '+'
                    })
                rc_pattern = pattern  # regex won't simple revcomp; scan both
                for m in re.finditer(revcomp(pattern.replace('[GC]', 'G').replace('...', 'NNN')).replace('N', '.'), seq):
                    rel = (start + m.start()) - (tss - 1)
                    records.append({
                        'gene_id': gid, 'strand': strand,
                        'pattern': pname, 'match_seq': m.group(),
                        'abs_pos': start + m.start() + 1,
                        'rel_pos': rel, 'match_strand': '-'
                    })
        else:
            start = max(0, tss - window_down - 1)
            end = min(len(genome_seq), tss + window_up)
            seq = genome_seq[start:end]
            for pname, pattern in patterns.items():
                for m in re.finditer(pattern, seq):
                    rel = (tss - 1) - (start + m.start())
                    records.append({
                        'gene_id': gid, 'strand': strand,
                        'pattern': pname, 'match_seq': m.group(),
                        'abs_pos': start + m.start() + 1,
                        'rel_pos': rel, 'match_strand': '+'
                    })
                for m in re.finditer(revcomp(pattern.replace('[GC]', 'G').replace('...', 'NNN')).replace('N', '.'), seq):
                    rel = (tss - 1) - (start + m.start())
                    records.append({
                        'gene_id': gid, 'strand': strand,
                        'pattern': pname, 'match_seq': m.group(),
                        'abs_pos': start + m.start() + 1,
                        'rel_pos': rel, 'match_strand': '-'
                    })

    return pd.DataFrame(records)


def main():
    print("=" * 60)
    print("SARP Binding Site – Methylation Spatial Integration")
    print("=" * 60)

    genes_df = pd.read_csv(TSS_TABLE)
    methyl_df = pd.read_csv(METHYL_SITES)
    integrated_df = pd.read_csv(INTEGRATED_CSV)
    bgc_df = pd.read_csv(BGC_DEF, sep='\t')
    fimo_df = pd.read_csv(FIMO_RESULTS, sep='\t')
    genome_seq = load_genome(GENOME_FASTA)

    bgc_genes = set(bgc_df['gene_id'].values)
    print(f"  Genes: {len(genes_df)}, BGC genes: {len(bgc_genes)}")
    print(f"  FIMO SARP hits in BGCs: {len(fimo_df)}")

    # ── 1. SARP binding zone analysis ──
    # Define SARP binding zone: -60 to -20 bp from TSS (encompasses -35 box ± margin)
    SARP_ZONE = (-60, -20)
    print(f"\n[1] SARP binding zone defined as {SARP_ZONE[0]} to {SARP_ZONE[1]} bp from TSS")

    # Map methylation sites to genes in SARP zone
    gene_tss = genes_df.set_index('gene_id')[['tss', 'strand', 'start', 'end']].to_dict('index')

    sarp_zone_methyl = []
    for _, site in methyl_df.iterrows():
        pos = site['position']
        tp = site['timepoint']
        mod = site['mod_type']
        freq = site['weighted_mod_freq']

        for gid, ginfo in gene_tss.items():
            tss = ginfo['tss']
            strand = ginfo['strand']
            rel = (pos - tss) if strand == '+' else (tss - pos)
            if SARP_ZONE[0] <= rel <= SARP_ZONE[1]:
                sarp_zone_methyl.append({
                    'gene_id': gid, 'position': pos, 'rel_pos': rel,
                    'mod_type': mod, 'timepoint': tp, 'freq': freq,
                    'is_bgc': gid in bgc_genes,
                })

    sz_df = pd.DataFrame(sarp_zone_methyl)
    print(f"  Methylation sites in SARP zone: {len(sz_df)}")
    print(f"    BGC genes: {sz_df[sz_df['is_bgc']]['gene_id'].nunique()}")
    print(f"    Non-BGC genes: {sz_df[~sz_df['is_bgc']]['gene_id'].nunique()}")

    # ── 2. BGC-specific SARP zone methylation ──
    print(f"\n[2] BGC gene SARP zone methylation detail:")

    bgc_methyl_detail = sz_df[sz_df['is_bgc']].merge(
        bgc_df[['gene_id', 'bgc_name', 'gene_name', 'role']], on='gene_id', how='left'
    )

    if len(bgc_methyl_detail) > 0:
        for bgc in ['act', 'red', 'cpk', 'cda']:
            sub = bgc_methyl_detail[bgc_methyl_detail['bgc_name'] == bgc]
            if len(sub) > 0:
                print(f"\n  {bgc.upper()} BGC ({len(sub)} methylation sites in SARP zone):")
                for _, r in sub.drop_duplicates(['gene_id', 'mod_type', 'timepoint']).iterrows():
                    print(f"    {r['gene_name']} ({r['gene_id']}): {r['mod_type']} at {r['rel_pos']:+d} bp, "
                          f"{r['timepoint']}, freq={r['freq']:.1f}%")
    else:
        print("  No methylation sites found in SARP zone of BGC genes")

    # ── 3. AAGCCCG and CCGG in SARP zone ──
    print(f"\n[3] Methylation target motifs in SARP zone (-60 to -20 bp):")

    motif_in_sarp = []
    for _, gene in genes_df.iterrows():
        tss = gene['tss']
        strand = gene['strand']
        gid = gene['gene_id']

        if strand == '+':
            zone_start = max(0, tss + SARP_ZONE[0] - 1)
            zone_end = tss + SARP_ZONE[1]
            zone_seq = genome_seq[zone_start:zone_end]
        else:
            zone_start = max(0, tss - SARP_ZONE[1] - 1)
            zone_end = tss - SARP_ZONE[0]
            zone_seq = genome_seq[zone_start:zone_end]
            zone_seq = zone_seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]

        for motif_name, motif_seq in [('AAGCCCG', 'AAGCCCG'), ('CCGG', 'CCGG')]:
            # Count on both strands
            fwd = len(re.findall(motif_seq, zone_seq))
            rev = len(re.findall(revcomp(motif_seq), zone_seq))
            if fwd + rev > 0:
                motif_in_sarp.append({
                    'gene_id': gid, 'motif': motif_name,
                    'count': fwd + rev, 'is_bgc': gid in bgc_genes,
                })

    motif_sarp_df = pd.DataFrame(motif_in_sarp)
    for motif_name in ['AAGCCCG', 'CCGG']:
        sub = motif_sarp_df[motif_sarp_df['motif'] == motif_name]
        total_genes = sub['gene_id'].nunique()
        bgc_count = sub[sub['is_bgc']]['gene_id'].nunique()
        print(f"  {motif_name} in SARP zone: {total_genes} genes ({bgc_count} BGC)")

    # ── 4. SARP binding motif scan ──
    print(f"\n[4] Scanning SARP consensus motifs (TCGAGC[GC], GCTCGAA) in -100 to +10...")
    sarp_motif_df = scan_sarp_motifs(genome_seq, genes_df)
    print(f"  SARP motif occurrences found: {len(sarp_motif_df)}")
    for p in sarp_motif_df['pattern'].unique():
        n = (sarp_motif_df['pattern'] == p).sum()
        print(f"    {p}: {n}")

    # ── 5. Triple overlap: SARP motif ∩ methylation motif ∩ -35 region ──
    print(f"\n[5] Checking SARP motif – methylation motif co-occurrence in SARP zone...")

    # Find genes where both SARP motif AND methylation motif occur in -60 to -20
    sarp_zone_sarp = sarp_motif_df[(sarp_motif_df['rel_pos'] >= SARP_ZONE[0]) & (sarp_motif_df['rel_pos'] <= SARP_ZONE[1])]
    genes_with_sarp_motif = set(sarp_zone_sarp['gene_id'])
    genes_with_methyl_motif = set(motif_sarp_df['gene_id'])
    genes_with_methyl_site = set(sz_df['gene_id'])

    overlap_sarp_methyl_motif = genes_with_sarp_motif & genes_with_methyl_motif
    overlap_sarp_methyl_site = genes_with_sarp_motif & genes_with_methyl_site
    triple_overlap = genes_with_sarp_motif & genes_with_methyl_motif & genes_with_methyl_site

    print(f"  Genes with SARP motif in zone: {len(genes_with_sarp_motif)}")
    print(f"  Genes with methyl target motif in zone: {len(genes_with_methyl_motif)}")
    print(f"  Genes with actual methylation in zone: {len(genes_with_methyl_site)}")
    print(f"  SARP motif ∩ methyl motif: {len(overlap_sarp_methyl_motif)}")
    print(f"  SARP motif ∩ actual methylation: {len(overlap_sarp_methyl_site)}")
    print(f"  Triple overlap (all three): {len(triple_overlap)}")

    # ── 6. Expression analysis: SARP zone methylation vs BGC expression ──
    print(f"\n[6] SARP zone methylation – BGC expression correlation:")

    # For BGC genes, compute methylation in SARP zone per timepoint
    bgc_expression = integrated_df[integrated_df['gene_id'].isin(bgc_genes)].copy()
    bgc_zone = sz_df[sz_df['is_bgc']].copy()

    if len(bgc_zone) > 0:
        # Aggregate SARP zone methylation per gene per timepoint
        zone_agg = bgc_zone.groupby(['gene_id', 'timepoint', 'mod_type']).agg(
            mean_freq=('freq', 'mean'), n_sites=('freq', 'count')
        ).reset_index()

        zone_pivot = zone_agg.pivot_table(
            index=['gene_id', 'mod_type'],
            columns='timepoint', values='mean_freq', aggfunc='mean'
        ).reset_index()

        for tp_a, tp_b, label in [('T2', 'T1', 'T2_vs_T1'), ('T3', 'T1', 'T3_vs_T1')]:
            if tp_a in zone_pivot.columns and tp_b in zone_pivot.columns:
                zone_pivot[f'sarp_zone_methyl_change_{label}'] = zone_pivot[tp_a] - zone_pivot[tp_b]

        zone_merged = zone_pivot.merge(
            bgc_expression[['gene_id', 'log2FC_T2_vs_T1', 'padj_T2_vs_T1',
                            'log2FC_T3_vs_T1', 'padj_T3_vs_T1']],
            on='gene_id', how='inner'
        )
        zone_merged = zone_merged.merge(
            bgc_df[['gene_id', 'bgc_name', 'gene_name']], on='gene_id', how='left'
        )

        print(f"  BGC genes with SARP zone methylation + expression data: {zone_merged['gene_id'].nunique()}")

        # Print detailed table
        if len(zone_merged) > 0:
            zone_merged.to_csv(OUTPUT_DIR / 'sarp_zone_bgc_detail.csv', index=False)
            print(f"\n  Detailed BGC SARP zone methylation-expression:")
            for _, r in zone_merged.iterrows():
                mc = r.get('sarp_zone_methyl_change_T2_vs_T1', np.nan)
                ec = r.get('log2FC_T2_vs_T1', np.nan)
                print(f"    {r.get('gene_name','?'):>15s} ({r['gene_id']}) {r.get('bgc_name','?'):>4s} {r['mod_type']}: "
                      f"Δmethyl={mc:+.1f}% log2FC={ec:+.2f}" if not np.isnan(mc) and not np.isnan(ec)
                      else f"    {r.get('gene_name','?'):>15s} ({r['gene_id']}) incomplete data")

    # ── 7. Genome-wide: SARP zone methylation vs expression ──
    print(f"\n[7] Genome-wide SARP zone methylation-expression correlation:")

    # Aggregate for all genes
    all_zone_agg = sz_df.groupby(['gene_id', 'timepoint', 'mod_type']).agg(
        mean_freq=('freq', 'mean'), n_sites=('freq', 'count')
    ).reset_index()

    all_zone_pivot = all_zone_agg.pivot_table(
        index=['gene_id', 'mod_type'],
        columns='timepoint', values='mean_freq', aggfunc='mean'
    ).reset_index()

    for tp_a, tp_b, label in [('T2', 'T1', 'T2_vs_T1'), ('T3', 'T1', 'T3_vs_T1'), ('T3', 'T2', 'T3_vs_T2')]:
        if tp_a in all_zone_pivot.columns and tp_b in all_zone_pivot.columns:
            all_zone_pivot[f'sarp_zone_change_{label}'] = all_zone_pivot[tp_a] - all_zone_pivot[tp_b]

    all_merged = all_zone_pivot.merge(
        integrated_df[['gene_id', 'log2FC_T2_vs_T1', 'padj_T2_vs_T1',
                        'log2FC_T3_vs_T1', 'padj_T3_vs_T1',
                        'log2FC_T3_vs_T2', 'padj_T3_vs_T2']],
        on='gene_id', how='inner'
    )

    corr_results = []
    for mod in ['6mA', '4mC']:
        for label in ['T2_vs_T1', 'T3_vs_T1', 'T3_vs_T2']:
            mc = f'sarp_zone_change_{label}'
            ec = f'log2FC_{label}'
            sub = all_merged[(all_merged['mod_type'] == mod)].dropna(subset=[mc, ec])
            sub = sub[sub[mc].abs() > 0]
            if len(sub) >= 10:
                r, p = stats.spearmanr(sub[mc], sub[ec])
                corr_results.append({
                    'mod_type': mod, 'comparison': label,
                    'region': 'SARP_zone_(-60_to_-20)',
                    'n_genes': len(sub), 'spearman_r': r, 'p_value': p,
                })
                sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
                print(f"  {mod} {label}: r={r:+.3f} p={p:.2e} n={len(sub)} {sig}")

    corr_df = pd.DataFrame(corr_results)
    corr_df.to_csv(OUTPUT_DIR / 'sarp_zone_correlation.csv', index=False)

    # ══════════════════════════════════════════════
    # FIGURES
    # ══════════════════════════════════════════════

    # ── Figure 1: Comprehensive SARP-methylation spatial map ──
    fig, ax = plt.subplots(figsize=(14, 8))

    y_base = 5.0
    # DNA backbone
    ax.add_patch(Rectangle((-120, y_base - 0.15), 240, 0.3, fc='#f0f0f0', ec='grey', lw=0.5))
    ax.annotate('', xy=(120, y_base), xytext=(-120, y_base),
                arrowprops=dict(arrowstyle='->', color='grey', lw=1))

    # Promoter elements
    ax.add_patch(Rectangle((-45, y_base - 0.35), 20, 0.7, fc='purple', alpha=0.3, ec='purple', lw=1.5))
    ax.text(-35, y_base + 0.55, '-35 box', ha='center', fontsize=10, color='purple', fontweight='bold')
    ax.add_patch(Rectangle((-15, y_base - 0.35), 10, 0.7, fc='orange', alpha=0.3, ec='orange', lw=1.5))
    ax.text(-10, y_base + 0.55, '-10 box', ha='center', fontsize=10, color='orange', fontweight='bold')
    ax.axvline(0, color='green', ls='-', lw=2.5)
    ax.text(3, y_base + 0.55, 'TSS', ha='left', fontsize=10, color='green', fontweight='bold')

    # SARP binding zone
    ax.add_patch(Rectangle((-60, y_base - 0.6), 40, 1.2, fc='none', ec='#e31a1c', lw=2, ls='--'))
    ax.text(-40, y_base + 0.8, 'SARP binding zone', ha='center', fontsize=10, color='#e31a1c', fontweight='bold')

    # SARP protein representation
    ax.add_patch(Rectangle((-55, y_base + 1.2), 15, 0.5, fc='#e31a1c', alpha=0.5, ec='#e31a1c', lw=1))
    ax.text(-47.5, y_base + 1.45, 'SARP₁', ha='center', fontsize=8, color='white', fontweight='bold')
    ax.add_patch(Rectangle((-38, y_base + 1.2), 15, 0.5, fc='#e31a1c', alpha=0.5, ec='#e31a1c', lw=1))
    ax.text(-30.5, y_base + 1.45, 'SARP₂', ha='center', fontsize=8, color='white', fontweight='bold')
    ax.text(-40, y_base + 1.9, '11 bp repeat\n(TCGAGC[G/C])', ha='center', fontsize=8, color='#e31a1c')

    # σHrdB sitting on SARP
    ax.add_patch(Rectangle((-42, y_base + 2.5), 20, 0.4, fc='#4daf4a', alpha=0.5, ec='#4daf4a', lw=1))
    ax.text(-32, y_base + 2.7, 'σHrdB R4', ha='center', fontsize=8, color='darkgreen', fontweight='bold')
    ax.annotate('', xy=(-35, y_base + 2.5), xytext=(-35, y_base + 1.7),
                arrowprops=dict(arrowstyle='->', color='#4daf4a', lw=1.5))

    # Methylation motif densities in SARP zone
    # AAGCCCG
    aag_in_zone = motif_sarp_df[motif_sarp_df['motif'] == 'AAGCCCG']
    ccgg_in_zone = motif_sarp_df[motif_sarp_df['motif'] == 'CCGG']

    ax.text(-115, y_base - 1.5, 'Methylation target motifs\nin SARP zone:', fontsize=10, fontweight='bold')
    ax.text(-115, y_base - 2.2, f'AAGCCCG: {len(aag_in_zone)} genes\nCCGG: {len(ccgg_in_zone)} genes',
            fontsize=9, color='#555')

    # Methylation sites in SARP zone
    if len(sz_df) > 0:
        zone_6mA = len(sz_df[sz_df['mod_type'] == '6mA'])
        zone_4mC = len(sz_df[sz_df['mod_type'] == '4mC'])
        ax.text(20, y_base - 1.5, 'Actual methylation sites\nin SARP zone:', fontsize=10, fontweight='bold')
        ax.text(20, y_base - 2.2, f'6mA: {zone_6mA} sites\n4mC: {zone_4mC} sites\n'
                f'({sz_df["gene_id"].nunique()} genes total)',
                fontsize=9, color='#555')

    # Methylation-free zone annotation
    ax.add_patch(Rectangle((-20, y_base - 0.8), 25, 0.1, fc='#2166ac', alpha=0.3))
    ax.text(-7, y_base - 1.1, '"Methyl-free zone"\n(D1: p=1.73e-12)', ha='center', fontsize=8, color='#2166ac')

    # Key question
    ax.text(-40, y_base - 3.2,
            'Key hypothesis: Methylation at SARP zone (-60 to -20)\n'
            'could block SARP binding → prevent σHrdB recruitment\n'
            '→ silence BGC gene transcription',
            fontsize=10, ha='center', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.4', fc='lightyellow', ec='grey', alpha=0.8))

    ax.set_xlim(-125, 125)
    ax.set_ylim(-4.5, 7.5)
    ax.set_xlabel('Distance from TSS (bp)')
    ax.set_yticks([])
    ax.set_title('SARP-mediated transcription activation: spatial relationship with methylation', fontsize=13)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / 'sarp_methylation_mechanism_schematic.png')
    fig.savefig(OUTPUT_DIR / 'sarp_methylation_mechanism_schematic.pdf')
    plt.close(fig)

    # ── Figure 2: SARP zone methylation density: BGC vs non-BGC ──
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for ax, mod in zip(axes, ['6mA', '4mC']):
        bgc_zone_sub = sz_df[(sz_df['is_bgc']) & (sz_df['mod_type'] == mod)]
        non_bgc_sub = sz_df[(~sz_df['is_bgc']) & (sz_df['mod_type'] == mod)]

        bins = np.arange(SARP_ZONE[0], SARP_ZONE[1] + 3, 3)

        if len(non_bgc_sub) > 0:
            ax.hist(non_bgc_sub['rel_pos'], bins=bins, color='grey', alpha=0.5,
                    label=f'Non-BGC (n={len(non_bgc_sub)})', density=True)
        if len(bgc_zone_sub) > 0:
            ax.hist(bgc_zone_sub['rel_pos'], bins=bins, color='#e31a1c', alpha=0.7,
                    label=f'BGC (n={len(bgc_zone_sub)})', density=True)

        ax.axvspan(-45, -25, alpha=0.15, color='purple', label='-35 box')
        ax.set_xlabel('Distance from TSS (bp)')
        ax.set_ylabel('Density')
        ax.set_title(f'{mod} in SARP zone')
        ax.legend(fontsize=8)

    fig.suptitle('Methylation in SARP binding zone: BGC vs non-BGC genes', fontsize=13, y=1.02)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / 'sarp_zone_bgc_vs_nonbgc.png')
    fig.savefig(OUTPUT_DIR / 'sarp_zone_bgc_vs_nonbgc.pdf')
    plt.close(fig)

    # ── Figure 3: BGC gene-level SARP zone methylation heatmap ──
    if len(bgc_methyl_detail) > 0:
        fig, ax = plt.subplots(figsize=(12, max(4, len(bgc_methyl_detail['gene_id'].unique()) * 0.4)))

        pivot_data = bgc_methyl_detail.pivot_table(
            index=['bgc_name', 'gene_name', 'gene_id'],
            columns=['timepoint', 'mod_type'],
            values='freq', aggfunc='mean'
        )

        if len(pivot_data) > 0:
            # Sort by BGC
            pivot_data = pivot_data.sort_index(level='bgc_name')
            labels = [f"{idx[0]}:{idx[1]} ({idx[2]})" for idx in pivot_data.index]

            im = ax.imshow(pivot_data.values, cmap='YlOrRd', aspect='auto', vmin=0, vmax=100)
            ax.set_yticks(range(len(labels)))
            ax.set_yticklabels(labels, fontsize=8)
            col_labels = [f"{c[0]}_{c[1]}" for c in pivot_data.columns]
            ax.set_xticks(range(len(col_labels)))
            ax.set_xticklabels(col_labels, rotation=45, ha='right', fontsize=8)

            for i in range(pivot_data.shape[0]):
                for j in range(pivot_data.shape[1]):
                    v = pivot_data.values[i, j]
                    if not np.isnan(v):
                        ax.text(j, i, f'{v:.0f}', ha='center', va='center', fontsize=7,
                                color='white' if v > 60 else 'black')

            plt.colorbar(im, ax=ax, label='Methylation freq (%)', shrink=0.6)
            ax.set_title('BGC gene SARP zone methylation (freq %)')

        fig.tight_layout()
        fig.savefig(OUTPUT_DIR / 'sarp_zone_bgc_heatmap.png')
        fig.savefig(OUTPUT_DIR / 'sarp_zone_bgc_heatmap.pdf')
        plt.close(fig)

    # ── Figure 4: Comparison of correlation strength across promoter regions ──
    # Load previous C1 results for comparison
    c1_path = OUTPUT_DIR / 'C1_distance_stratified_correlation.csv'
    if c1_path.exists():
        c1_df = pd.read_csv(c1_path)
        # Add SARP zone results
        combined = pd.concat([c1_df, corr_df], ignore_index=True)
        combined.to_csv(OUTPUT_DIR / 'sarp_zone_vs_other_regions_correlation.csv', index=False)

        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        for ax, mod in zip(axes, ['6mA', '4mC']):
            sub = combined[(combined['mod_type'] == mod) & (combined['comparison'] == 'T2_vs_T1')]
            if len(sub) > 0:
                # Custom order
                order = ['distal_upstream', 'mid_upstream', 'proximal_upstream',
                         'SARP_zone_(-60_to_-20)', 'core_promoter', 'TSS_proximal',
                         'early_gene_body', 'mid_gene_body', 'distal_gene_body']
                sub_ordered = sub.set_index('region' if 'region' in sub.columns else 'distance_bin')
                existing = [o for o in order if o in sub_ordered.index]
                sub_ordered = sub_ordered.loc[existing]

                colors = ['#2166ac' if 'SARP' not in r else '#e31a1c' for r in existing]
                bars = ax.barh(range(len(existing)), sub_ordered['spearman_r'], color=colors, edgecolor='white')
                ax.set_yticks(range(len(existing)))
                ax.set_yticklabels(existing, fontsize=9)
                ax.axvline(0, color='grey', lw=0.5)
                ax.set_xlabel('Spearman r')
                ax.set_title(f'{mod} (T2 vs T1)')

                # Add significance annotations
                for i, (idx, r) in enumerate(sub_ordered.iterrows()):
                    p = r['p_value']
                    sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
                    ax.text(r['spearman_r'] + 0.01 * np.sign(r['spearman_r']), i,
                            f"{r['spearman_r']:.2f}{sig}", va='center', fontsize=8)

        fig.suptitle('Methylation-Expression Correlation by Region (including SARP zone)', fontsize=13, y=1.02)
        fig.tight_layout()
        fig.savefig(OUTPUT_DIR / 'sarp_zone_correlation_comparison.png')
        fig.savefig(OUTPUT_DIR / 'sarp_zone_correlation_comparison.pdf')
        plt.close(fig)

    # Save summary
    summary = {
        'sarp_zone': SARP_ZONE,
        'methyl_sites_in_zone': len(sz_df),
        'genes_with_zone_methyl': sz_df['gene_id'].nunique(),
        'bgc_genes_with_zone_methyl': sz_df[sz_df['is_bgc']]['gene_id'].nunique(),
        'aagcccg_in_zone': len(aag_in_zone),
        'ccgg_in_zone': len(ccgg_in_zone),
        'genes_with_sarp_motif': len(genes_with_sarp_motif),
        'triple_overlap': len(triple_overlap),
    }

    print(f"\n  All figures saved to {OUTPUT_DIR}")
    print("=" * 60)


if __name__ == '__main__':
    main()
