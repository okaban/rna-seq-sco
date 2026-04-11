#!/usr/bin/env python3
"""
Temporal Methylation Motif Analysis with Expression Association
- T1/T2/T3 timepoint-specific motif enrichment
- Gained/Lost site motif characterization
- Expression-associated motif identification per comparison
"""

import pandas as pd
import numpy as np
from pathlib import Path
import re
from collections import Counter, defaultdict
from scipy import stats
from scipy.stats import fisher_exact
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings
warnings.filterwarnings('ignore')

METHYL_SITES = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/01_integration/high_confidence_sites_weighted.csv")
INTEGRATED_CSV = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/01_integration/integrated_methyl_expression_weighted.csv")
TSS_TABLE = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/18_tss_analyses/comprehensive_tss_table.csv")
GENOME_FASTA = "/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/GCF_000203835.1_ASM20383v1_genomic.fna"
OUTPUT_DIR = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/18_tss_analyses")

plt.rcParams.update({
    'font.size': 10, 'axes.labelsize': 11, 'axes.titlesize': 12,
    'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
    'font.family': 'sans-serif',
})

def load_genome(path):
    parts = []
    with open(path) as f:
        for line in f:
            if not line.startswith('>'):
                parts.append(line.strip())
    return ''.join(parts).upper()

def revcomp(seq):
    return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]

def extract_flanking(genome, positions, flank=10):
    """Extract flanking sequences around positions."""
    seqs = []
    for pos, strand in positions:
        start = pos - flank - 1  # 0-indexed
        end = pos + flank
        if start < 0 or end > len(genome):
            continue
        seq = genome[start:end]
        if strand == '-':
            seq = revcomp(seq)
        seqs.append(seq)
    return seqs

def count_kmers(seqs, k=4):
    """Count k-mer frequencies in sequences."""
    counts = Counter()
    total = 0
    for s in seqs:
        for i in range(len(s) - k + 1):
            kmer = s[i:i+k]
            if 'N' not in kmer:
                counts[kmer] += 1
                total += 1
    # Normalize
    freqs = {kmer: count / total if total > 0 else 0 for kmer, count in counts.items()}
    return counts, freqs, total

def known_motif_enrichment(seqs, motifs_dict):
    """Count known motif occurrences in sequences."""
    results = {}
    for name, pattern in motifs_dict.items():
        count = 0
        for s in seqs:
            if re.search(pattern, s):
                count += 1
        results[name] = count
    return results

def main():
    print("=" * 60)
    print("Temporal Methylation Motif – Expression Analysis")
    print("=" * 60)

    genome = load_genome(GENOME_FASTA)
    methyl_df = pd.read_csv(METHYL_SITES)
    integrated_df = pd.read_csv(INTEGRATED_CSV)
    genes_df = pd.read_csv(TSS_TABLE)

    # Known methylation motifs to test
    KNOWN_MOTIFS = {
        'CCGG': 'CCGG',
        'CCGG_rc': 'CCGG',
        'AAGCCCG': 'AAGCCCG',
        'CCCGCTT': 'CCCGCTT',  # revcomp of AAGCCCG
        'GCGC': 'GCGC',
        'GATC': 'GATC',
        'CTCGAG': 'CTCGAG',  # XhoI-like
        'TCGA': 'TCGA',  # TaqI-like / SARP core
    }

    # ══════════════════════════════════════════════
    # 1. Timepoint-specific motif enrichment
    # ══════════════════════════════════════════════
    print("\n[1] Timepoint-specific motif analysis")

    tp_motif_results = {}
    tp_kmer_freqs = {}

    for mod in ['6mA', '4mC']:
        for tp in ['T1', 'T2', 'T3']:
            sub = methyl_df[(methyl_df['mod_type'] == mod) & (methyl_df['timepoint'] == tp)]
            positions = list(zip(sub['position'].values, sub['strand'].values))
            seqs = extract_flanking(genome, positions, flank=10)

            # Known motifs
            motif_counts = known_motif_enrichment(seqs, KNOWN_MOTIFS)
            n_seqs = len(seqs)
            tp_motif_results[(mod, tp)] = {k: v / n_seqs * 100 if n_seqs > 0 else 0
                                            for k, v in motif_counts.items()}
            tp_motif_results[(mod, tp)]['n_sites'] = n_seqs

            # k-mer analysis
            _, freqs, _ = count_kmers(seqs, k=4)
            tp_kmer_freqs[(mod, tp)] = freqs

            print(f"  {mod} {tp}: {n_seqs} sites")
            top3 = sorted(motif_counts.items(), key=lambda x: -x[1])[:3]
            for name, cnt in top3:
                print(f"    {name}: {cnt} ({cnt/n_seqs*100:.1f}%)" if n_seqs > 0 else f"    {name}: 0")

    # ── Table: motif % per timepoint ──
    motif_table_rows = []
    for mod in ['6mA', '4mC']:
        for motif_name in ['CCGG', 'AAGCCCG', 'GCGC', 'GATC', 'TCGA']:
            row = {'mod_type': mod, 'motif': motif_name}
            for tp in ['T1', 'T2', 'T3']:
                row[f'{tp}_pct'] = tp_motif_results[(mod, tp)].get(motif_name, 0)
                row[f'{tp}_n'] = tp_motif_results[(mod, tp)].get('n_sites', 0)
            motif_table_rows.append(row)

    motif_table = pd.DataFrame(motif_table_rows)
    motif_table.to_csv(OUTPUT_DIR / 'temporal_motif_enrichment.csv', index=False)

    # ══════════════════════════════════════════════
    # 2. Gained/Lost site motif analysis
    # ══════════════════════════════════════════════
    print("\n[2] Gained/Lost site motif characterization")

    sites_by_tp = {}
    for tp in ['T1', 'T2', 'T3']:
        sub = methyl_df[methyl_df['timepoint'] == tp]
        sites_by_tp[tp] = {(r['position'], r['strand'], r['mod_type']): r
                           for _, r in sub.iterrows()}

    comparisons = [('T1', 'T2', 'T2_vs_T1'), ('T1', 'T3', 'T3_vs_T1'), ('T2', 'T3', 'T3_vs_T2')]
    gained_lost_results = []

    for tp_a, tp_b, label in comparisons:
        keys_a = {(p, s, m) for p, s, m in sites_by_tp[tp_a]}
        keys_b = {(p, s, m) for p, s, m in sites_by_tp[tp_b]}

        gained = keys_b - keys_a
        lost = keys_a - keys_b
        shared = keys_a & keys_b

        for mod in ['6mA', '4mC']:
            gained_mod = [(p, s) for p, s, m in gained if m == mod]
            lost_mod = [(p, s) for p, s, m in lost if m == mod]
            shared_mod = [(p, s) for p, s, m in shared if m == mod]

            seqs_gained = extract_flanking(genome, gained_mod, flank=10)
            seqs_lost = extract_flanking(genome, lost_mod, flank=10)
            seqs_shared = extract_flanking(genome, shared_mod, flank=10)

            for category, seqs in [('gained', seqs_gained), ('lost', seqs_lost), ('shared', seqs_shared)]:
                motif_counts = known_motif_enrichment(seqs, KNOWN_MOTIFS)
                n = len(seqs)
                for motif_name, cnt in motif_counts.items():
                    gained_lost_results.append({
                        'comparison': label, 'mod_type': mod, 'category': category,
                        'motif': motif_name, 'count': cnt, 'n_sites': n,
                        'pct': cnt / n * 100 if n > 0 else 0,
                    })

            print(f"  {label} {mod}: gained={len(gained_mod)}, lost={len(lost_mod)}, shared={len(shared_mod)}")

            # Fisher's test: is CCGG enriched in gained vs lost?
            for test_motif in ['CCGG', 'AAGCCCG', 'TCGA']:
                g_count = sum(1 for s in seqs_gained if test_motif in s)
                l_count = sum(1 for s in seqs_lost if test_motif in s)
                g_no = len(seqs_gained) - g_count
                l_no = len(seqs_lost) - l_count
                if min(len(seqs_gained), len(seqs_lost)) > 5:
                    _, p = fisher_exact([[g_count, g_no], [l_count, l_no]])
                    g_pct = g_count / len(seqs_gained) * 100 if len(seqs_gained) > 0 else 0
                    l_pct = l_count / len(seqs_lost) * 100 if len(seqs_lost) > 0 else 0
                    if p < 0.05:
                        print(f"    ** {test_motif} gained={g_pct:.1f}% vs lost={l_pct:.1f}%, Fisher p={p:.2e}")

    gl_df = pd.DataFrame(gained_lost_results)
    gl_df.to_csv(OUTPUT_DIR / 'temporal_gained_lost_motifs.csv', index=False)

    # ══════════════════════════════════════════════
    # 3. Expression-associated motif analysis
    # ══════════════════════════════════════════════
    print("\n[3] Expression-associated methylation motif analysis")

    gene_tss = genes_df.set_index('gene_id')[['tss', 'strand']].to_dict('index')

    # Assign each methylation site to nearest gene promoter
    site_gene_map = []
    for _, site in methyl_df.iterrows():
        pos = site['position']
        best_gene = None
        best_dist = 999999
        for gid, ginfo in gene_tss.items():
            tss = ginfo['tss']
            strand = ginfo['strand']
            rel = (pos - tss) if strand == '+' else (tss - pos)
            if -300 <= rel <= 50:
                dist = abs(rel)
                if dist < best_dist:
                    best_dist = dist
                    best_gene = gid
        if best_gene:
            site_gene_map.append({
                'position': pos, 'strand': site['strand'],
                'mod_type': site['mod_type'], 'timepoint': site['timepoint'],
                'freq': site['weighted_mod_freq'], 'gene_id': best_gene,
            })

    sgm_df = pd.DataFrame(site_gene_map)
    print(f"  Promoter-assigned methylation sites: {len(sgm_df)}")

    # Merge with expression
    expr_assoc_results = []

    for label, ec, pc in [('T2_vs_T1', 'log2FC_T2_vs_T1', 'padj_T2_vs_T1'),
                          ('T3_vs_T1', 'log2FC_T3_vs_T1', 'padj_T3_vs_T1'),
                          ('T3_vs_T2', 'log2FC_T3_vs_T2', 'padj_T3_vs_T2')]:
        if ec not in integrated_df.columns:
            continue

        expr_sub = integrated_df[['gene_id', ec, pc]].dropna()
        sig_up = set(expr_sub[(expr_sub[pc] < 0.05) & (expr_sub[ec] > 1)]['gene_id'])
        sig_down = set(expr_sub[(expr_sub[pc] < 0.05) & (expr_sub[ec] < -1)]['gene_id'])
        non_deg = set(expr_sub[expr_sub[pc] >= 0.05]['gene_id'])

        print(f"\n  {label}: UP={len(sig_up)}, DOWN={len(sig_down)}, non-DEG={len(non_deg)}")

        # Determine gained/lost for this comparison
        tp_map = {'T2_vs_T1': ('T1', 'T2'), 'T3_vs_T1': ('T1', 'T3'), 'T3_vs_T2': ('T2', 'T3')}
        tp_a, tp_b = tp_map[label]

        for mod in ['6mA', '4mC']:
            # Sites in tp_b but not tp_a (gained) that are in promoters
            sites_a = sgm_df[(sgm_df['timepoint'] == tp_a) & (sgm_df['mod_type'] == mod)]
            sites_b = sgm_df[(sgm_df['timepoint'] == tp_b) & (sgm_df['mod_type'] == mod)]

            pos_a = set(sites_a['position'])
            pos_b = set(sites_b['position'])

            gained_sites = sites_b[~sites_b['position'].isin(pos_a)]
            lost_sites = sites_a[~sites_a['position'].isin(pos_b)]

            # For gained sites: which genes are they in?
            for category, cat_sites, gene_set, set_name in [
                ('gained_in_UP', gained_sites, sig_up, 'UP'),
                ('gained_in_DOWN', gained_sites, sig_down, 'DOWN'),
                ('gained_in_nonDEG', gained_sites, non_deg, 'nonDEG'),
                ('lost_in_UP', lost_sites, sig_up, 'UP'),
                ('lost_in_DOWN', lost_sites, sig_down, 'DOWN'),
                ('lost_in_nonDEG', lost_sites, non_deg, 'nonDEG'),
            ]:
                cat_sub = cat_sites[cat_sites['gene_id'].isin(gene_set)]
                positions = list(zip(cat_sub['position'].values,
                                     cat_sub['strand'].values if 'strand' in cat_sub else ['+'] * len(cat_sub)))
                seqs = extract_flanking(genome, positions, flank=10)

                motif_counts = known_motif_enrichment(seqs, KNOWN_MOTIFS)
                n = len(seqs)
                for motif_name, cnt in motif_counts.items():
                    expr_assoc_results.append({
                        'comparison': label, 'mod_type': mod,
                        'site_category': category.split('_in_')[0],
                        'expr_category': set_name,
                        'motif': motif_name, 'count': cnt, 'n_sites': n,
                        'pct': cnt / n * 100 if n > 0 else 0,
                    })

            # Fisher test: motif enrichment in gained-UP vs gained-nonDEG
            for test_motif in ['CCGG', 'AAGCCCG', 'TCGA']:
                gained_up = gained_sites[gained_sites['gene_id'].isin(sig_up)]
                gained_nondeg = gained_sites[gained_sites['gene_id'].isin(non_deg)]

                seqs_up = extract_flanking(genome,
                    list(zip(gained_up['position'].values, gained_up['strand'].values)), flank=10)
                seqs_nd = extract_flanking(genome,
                    list(zip(gained_nondeg['position'].values, gained_nondeg['strand'].values)), flank=10)

                if len(seqs_up) >= 5 and len(seqs_nd) >= 5:
                    up_cnt = sum(1 for s in seqs_up if test_motif in s)
                    nd_cnt = sum(1 for s in seqs_nd if test_motif in s)
                    _, p = fisher_exact([[up_cnt, len(seqs_up) - up_cnt],
                                         [nd_cnt, len(seqs_nd) - nd_cnt]])
                    if p < 0.1:
                        print(f"    {mod} {label} gained {test_motif}: UP={up_cnt}/{len(seqs_up)} "
                              f"({up_cnt/len(seqs_up)*100:.1f}%) vs nonDEG={nd_cnt}/{len(seqs_nd)} "
                              f"({nd_cnt/len(seqs_nd)*100:.1f}%), p={p:.3f}")

    ea_df = pd.DataFrame(expr_assoc_results)
    ea_df.to_csv(OUTPUT_DIR / 'temporal_expression_associated_motifs.csv', index=False)

    # ══════════════════════════════════════════════
    # FIGURES
    # ══════════════════════════════════════════════

    # ── Fig 1: Motif prevalence by timepoint (heatmap) ──
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    key_motifs = ['CCGG', 'AAGCCCG', 'GCGC', 'GATC', 'TCGA']

    for ax, mod in zip(axes, ['6mA', '4mC']):
        data = np.zeros((len(key_motifs), 3))
        for j, tp in enumerate(['T1', 'T2', 'T3']):
            for i, motif in enumerate(key_motifs):
                data[i, j] = tp_motif_results[(mod, tp)].get(motif, 0)

        im = ax.imshow(data, cmap='YlOrRd', aspect='auto', vmin=0)
        ax.set_xticks([0, 1, 2])
        n_sites = [tp_motif_results[(mod, tp)]['n_sites'] for tp in ['T1', 'T2', 'T3']]
        ax.set_xticklabels([f'T1\n(n={n_sites[0]})', f'T2\n(n={n_sites[1]})', f'T3\n(n={n_sites[2]})'])
        ax.set_yticks(range(len(key_motifs)))
        ax.set_yticklabels(key_motifs)
        for i in range(len(key_motifs)):
            for j in range(3):
                ax.text(j, i, f'{data[i,j]:.1f}%', ha='center', va='center', fontsize=9,
                        color='white' if data[i,j] > 40 else 'black')
        ax.set_title(f'{mod} motif prevalence (% of sites)')
        plt.colorbar(im, ax=ax, shrink=0.7, label='%')

    fig.suptitle('Methylation motif prevalence by timepoint', fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / 'temporal_motif_heatmap.png')
    fig.savefig(OUTPUT_DIR / 'temporal_motif_heatmap.pdf')
    plt.close(fig)

    # ── Fig 2: Gained vs Lost motif comparison ──
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    for row_idx, mod in enumerate(['6mA', '4mC']):
        for col_idx, (_, _, label) in enumerate(comparisons):
            ax = axes[row_idx, col_idx]
            sub = gl_df[(gl_df['mod_type'] == mod) & (gl_df['comparison'] == label)]

            for motif in key_motifs:
                msub = sub[sub['motif'] == motif]
                gained_pct = msub[msub['category'] == 'gained']['pct'].values
                lost_pct = msub[msub['category'] == 'lost']['pct'].values
                gained_pct = gained_pct[0] if len(gained_pct) > 0 else 0
                lost_pct = lost_pct[0] if len(lost_pct) > 0 else 0

            # Grouped bar chart
            x = np.arange(len(key_motifs))
            width = 0.25
            gained_vals = []
            lost_vals = []
            shared_vals = []
            for motif in key_motifs:
                msub = sub[sub['motif'] == motif]
                for cat, lst in [('gained', gained_vals), ('lost', lost_vals), ('shared', shared_vals)]:
                    v = msub[msub['category'] == cat]['pct'].values
                    lst.append(v[0] if len(v) > 0 else 0)

            ax.bar(x - width, gained_vals, width, label='Gained', color='#e31a1c', alpha=0.8)
            ax.bar(x, lost_vals, width, label='Lost', color='#2166ac', alpha=0.8)
            ax.bar(x + width, shared_vals, width, label='Shared', color='grey', alpha=0.6)

            ax.set_xticks(x)
            ax.set_xticklabels(key_motifs, rotation=45, ha='right')
            ax.set_ylabel('% of sites')
            ax.set_title(f'{mod} — {label}')
            n_g = sub[sub['category'] == 'gained']['n_sites'].max() if len(sub[sub['category'] == 'gained']) > 0 else 0
            n_l = sub[sub['category'] == 'lost']['n_sites'].max() if len(sub[sub['category'] == 'lost']) > 0 else 0
            ax.text(0.95, 0.95, f'G={int(n_g)}\nL={int(n_l)}', transform=ax.transAxes,
                    va='top', ha='right', fontsize=8, bbox=dict(fc='white', alpha=0.7))
            if row_idx == 0 and col_idx == 0:
                ax.legend(fontsize=8)

    fig.suptitle('Motif prevalence in Gained vs Lost methylation sites', fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / 'temporal_gained_lost_motif_barplot.png')
    fig.savefig(OUTPUT_DIR / 'temporal_gained_lost_motif_barplot.pdf')
    plt.close(fig)

    # ── Fig 3: Expression-associated motif heatmap ──
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    for row_idx, mod in enumerate(['6mA', '4mC']):
        for col_idx, (_, _, label) in enumerate(comparisons):
            ax = axes[row_idx, col_idx]
            sub = ea_df[(ea_df['mod_type'] == mod) & (ea_df['comparison'] == label)]

            # Build matrix: rows=motifs, cols=gained_UP, gained_DOWN, lost_UP, lost_DOWN
            categories = ['gained_UP', 'gained_DOWN', 'gained_nonDEG', 'lost_UP', 'lost_DOWN', 'lost_nonDEG']
            cat_labels = ['G→UP', 'G→DOWN', 'G→ns', 'L→UP', 'L→DOWN', 'L→ns']
            data = np.zeros((len(key_motifs), len(categories)))

            for i, motif in enumerate(key_motifs):
                for j, (sc, ec) in enumerate([(c.split('_')[0], c.split('_')[1]) for c in
                        ['gained_UP','gained_DOWN','gained_nonDEG','lost_UP','lost_DOWN','lost_nonDEG']]):
                    match = sub[(sub['motif'] == motif) & (sub['site_category'] == sc) & (sub['expr_category'] == ec)]
                    if len(match) > 0:
                        data[i, j] = match.iloc[0]['pct']

            im = ax.imshow(data, cmap='YlOrRd', aspect='auto', vmin=0)
            ax.set_xticks(range(len(cat_labels)))
            ax.set_xticklabels(cat_labels, rotation=45, ha='right', fontsize=8)
            ax.set_yticks(range(len(key_motifs)))
            ax.set_yticklabels(key_motifs, fontsize=9)
            for i in range(data.shape[0]):
                for j in range(data.shape[1]):
                    if data[i, j] > 0:
                        ax.text(j, i, f'{data[i,j]:.0f}', ha='center', va='center', fontsize=7,
                                color='white' if data[i, j] > 40 else 'black')
            ax.set_title(f'{mod} — {label}', fontsize=11)

    fig.suptitle('Motif prevalence (%) by methylation change × expression change', fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / 'temporal_expression_motif_heatmap.png')
    fig.savefig(OUTPUT_DIR / 'temporal_expression_motif_heatmap.pdf')
    plt.close(fig)

    # ── Fig 4: Key motif temporal dynamics summary ──
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    for ax_idx, (mod, motif) in enumerate([('4mC', 'CCGG'), ('6mA', 'AAGCCCG'), ('6mA', 'TCGA'), ('4mC', 'GCGC')]):
        ax = axes[ax_idx // 2, ax_idx % 2]
        tps = ['T1', 'T2', 'T3']
        pcts = [tp_motif_results[(mod, tp)].get(motif, 0) for tp in tps]
        ns = [tp_motif_results[(mod, tp)]['n_sites'] for tp in tps]

        ax.bar(tps, pcts, color=['#2166ac', '#d6604d', '#4daf4a'], edgecolor='white', width=0.5)
        for i, (p, n) in enumerate(zip(pcts, ns)):
            ax.text(i, p + 0.5, f'{p:.1f}%\n(n={n})', ha='center', fontsize=9)
        ax.set_ylabel('% of sites containing motif')
        ax.set_title(f'{mod} — {motif}')
        ax.set_ylim(0, max(pcts) * 1.3 if max(pcts) > 0 else 10)

    fig.suptitle('Key motif prevalence across timepoints', fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / 'temporal_key_motifs_barplot.png')
    fig.savefig(OUTPUT_DIR / 'temporal_key_motifs_barplot.pdf')
    plt.close(fig)

    print(f"\n  All figures saved to {OUTPUT_DIR}")
    print("=" * 60)


if __name__ == '__main__':
    main()
