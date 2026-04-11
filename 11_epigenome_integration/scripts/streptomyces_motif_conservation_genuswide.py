#!/usr/bin/env python3
"""
Genus-wide Streptomyces Motif Conservation Analysis
=====================================================
Analyze 7 methylation motifs (this study + prior studies) across all
RefSeq Streptomyces species using REBASE R-M data and genome sequences.

Motifs:
  This study (Nanopore, S. coelicolor M145):
    CCGG (4mC), AAGCCCG (6mA), GATC (6mA)
  Pisciotta et al. 2023 (BS-seq, S. coelicolor M145):
    GGCCGG (5mC), GCCCG (5mC)
  Fang et al. 2022 (SMRT, S. roseosporus L30):
    GCGG (4mC), CGACNNNCTCC (6mA)

Outputs:
  - Panel A: R-M system conservation bar chart
  - Panel B: Motif site density violin plots
  - Panel C: Species × motif heatmap (REBASE)
  - Panel D: Observed/Expected density scatter

Streptomyces coelicolor A3(2) M145 Project
"""

import re
import csv
import gzip
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool, cpu_count

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from matplotlib.colors import ListedColormap
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
import warnings
warnings.filterwarnings('ignore')

# ── Paths ──
BASE = Path('/Users/okaban/bioinfo/rna-seq/11_epigenome_integration')
REBASE_PATH = BASE / 'data' / 'rebase' / 'bairoch.txt'
GENOME_DIR = BASE / 'data' / 'ncbi_genomes_all_species'
SPECIES_TSV = BASE / 'analysis' / '21_genuswide_motif_conservation' / 'species_representative_accessions.tsv'
OUTPUT_DIR = BASE / 'analysis' / '21_genuswide_motif_conservation'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

M145_GENOME = Path(
    '/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/'
    'M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/'
    'GCF_000203835.1/GCF_000203835.1_ASM20383v1_genomic.fna'
)

# ── Publication style ──
COL_THIS = '#D32F2F'      # Red – this study
COL_THIS_ALT = '#E57373'  # Light red
COL_PRIOR = '#1565C0'     # Blue – prior studies
COL_PRIOR_ALT = '#42A5F5' # Light blue
COL_HIGHLIGHT = '#D32F2F'  # Red star for M145
COL_HIGHLIGHT2 = '#1565C0' # Blue diamond for L30

plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 9,
    'axes.titlesize': 10,
    'axes.labelsize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 0.7,
    'xtick.major.width': 0.7,
    'ytick.major.width': 0.7,
    'axes.spines.top': False,
    'axes.spines.right': False,
})

# ── Motifs ──
MOTIFS = [
    # (motif, mod_type, source, species, color)
    ('CCGG',         '4mC', 'This study',      'S. coelicolor M145', COL_THIS),
    ('AAGCCCG',      '6mA', 'This study',      'S. coelicolor M145', COL_THIS),
    ('GATC',         '6mA', 'This study',      'S. coelicolor M145', COL_THIS),
    ('GGCCGG',       '5mC', 'Pisciotta 2023',  'S. coelicolor M145', COL_PRIOR),
    ('GCCCG',        '5mC', 'Pisciotta 2023',  'S. coelicolor M145', COL_PRIOR),
    ('GCGG',         '4mC', 'Fang 2022',       'S. roseosporus L30', COL_PRIOR),
    ('CGACNNNCTCC',  '6mA', 'Fang 2022',       'S. roseosporus L30', COL_PRIOR),
]

MOTIF_NAMES = [m[0] for m in MOTIFS]

# ── IUPAC ──
IUPAC = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'},
    'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'G', 'C'},
    'W': {'A', 'T'}, 'K': {'G', 'T'}, 'M': {'A', 'C'},
    'B': {'C', 'G', 'T'}, 'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'},
    'N': {'A', 'C', 'G', 'T'},
}

COMPLEMENT = {
    'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
    'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
    'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
    'D': 'H', 'H': 'D', 'N': 'N',
}


def reverse_complement(seq):
    return ''.join(COMPLEMENT.get(b, 'N') for b in reversed(seq.upper()))


def iupac_match(pattern, target):
    if len(pattern) != len(target):
        return False
    for p, t in zip(pattern.upper(), target.upper()):
        if not (IUPAC.get(p, {p}) & IUPAC.get(t, {t})):
            return False
    return True


def specificity_score(pattern):
    score = sum(1.0 / len(IUPAC.get(b, {b})) for b in pattern.upper())
    return score / len(pattern) if pattern else 0


def motif_matches_recognition(motif, recognition):
    rec_clean = re.sub(r'[^A-Za-z]', '', str(recognition)).upper()
    motif_up = motif.upper()
    rc_motif = reverse_complement(motif_up)

    if len(rec_clean) == len(motif_up):
        return iupac_match(rec_clean, motif_up) or iupac_match(rec_clean, rc_motif)

    for i in range(len(rec_clean) - len(motif_up) + 1):
        sub = rec_clean[i:i + len(motif_up)]
        if (iupac_match(sub, motif_up) or iupac_match(sub, rc_motif)):
            if specificity_score(sub) >= 0.6:
                return True
    return False


def iupac_to_regex(motif):
    """Convert IUPAC motif to regex pattern."""
    base_map = {
        'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
        'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
        'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
        'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]',
    }
    return ''.join(base_map.get(b, b) for b in motif.upper())


# ══════════════════════════════════════════
# REBASE parsing
# ══════════════════════════════════════════

def parse_bairoch(filepath):
    records = []
    current = {}
    with open(filepath) as fh:
        for line in fh:
            if line.startswith('//'):
                if current:
                    records.append(current)
                current = {}
                continue
            if line.startswith('CC') or not line.strip():
                continue
            tag = line[:5].strip()
            value = line[5:].strip()
            if tag == 'ID':
                current['enzyme_name'] = value
            elif tag == 'ET':
                current['enzyme_type'] = value
            elif tag == 'AC':
                current['accession'] = value.rstrip(';')
            elif tag == 'OS':
                current['organism'] = value
            elif tag == 'RS':
                current.setdefault('recognition_raw', '')
                current['recognition_raw'] += (' ' + value if current['recognition_raw'] else value)
            elif tag == 'MS':
                current.setdefault('methylation_raw', '')
                current['methylation_raw'] += (' ' + value if current['methylation_raw'] else value)
    return records


def filter_streptomyces(records):
    strepto = [r for r in records if r.get('organism', '').startswith('Streptomyces')]
    df = pd.DataFrame(strepto)
    if df.empty:
        return df
    df['recognition_seq'] = df['recognition_raw'].apply(
        lambda x: re.sub(r'[^A-Za-z,; ]', '', str(x)).strip().split(',')[0].strip()
        if pd.notna(x) else ''
    )
    return df


def rebase_motif_conservation(rebase_df):
    """Compute R-M conservation rate for each motif across REBASE species."""
    # Get unique species
    species_organisms = defaultdict(set)
    for _, row in rebase_df.iterrows():
        org = row.get('organism', '')
        parts = org.split()
        if len(parts) >= 2:
            binomial = f'{parts[0]} {parts[1]}'
        else:
            binomial = org
        species_organisms[binomial].add(org)

    n_species = len(species_organisms)
    print(f'  REBASE: {n_species} unique Streptomyces species')

    results = []
    for motif, mod_type, source, det_species, color in MOTIFS:
        match_species = set()
        match_enzymes = []
        for sp, organisms in species_organisms.items():
            for org in organisms:
                org_entries = rebase_df[rebase_df['organism'] == org]
                for _, entry in org_entries.iterrows():
                    rec = entry.get('recognition_seq', '')
                    if rec and motif_matches_recognition(motif, rec):
                        match_species.add(sp)
                        match_enzymes.append(entry.get('enzyme_name', ''))

        n_match = len(match_species)
        rate = n_match / n_species * 100 if n_species > 0 else 0

        # Wilson confidence interval
        p = n_match / n_species if n_species > 0 else 0
        z = 1.96  # 95% CI
        denom = 1 + z**2 / n_species
        center = (p + z**2 / (2 * n_species)) / denom
        margin = z * np.sqrt(p * (1 - p) / n_species + z**2 / (4 * n_species**2)) / denom
        ci_low = max(0, center - margin) * 100
        ci_high = min(1, center + margin) * 100

        results.append({
            'motif': motif,
            'mod_type': mod_type,
            'source': source,
            'det_species': det_species,
            'n_match': n_match,
            'n_total': n_species,
            'rate_pct': rate,
            'ci_low': ci_low,
            'ci_high': ci_high,
            'enzymes': '; '.join(sorted(set(match_enzymes))),
            'color': color,
        })

    return pd.DataFrame(results)


def rebase_species_motif_matrix(rebase_df):
    """Build species × motif binary matrix from REBASE data."""
    # Get species → recognition seqs mapping
    species_recs = defaultdict(set)
    species_orgs = defaultdict(set)
    for _, row in rebase_df.iterrows():
        org = row.get('organism', '')
        parts = org.split()
        binomial = f'{parts[0]} {parts[1]}' if len(parts) >= 2 else org
        rec = row.get('recognition_seq', '')
        if rec:
            species_recs[binomial].add(rec)
        species_orgs[binomial].add(org)

    # Build matrix
    rows = []
    for sp in sorted(species_recs.keys()):
        row = {'species': sp}
        for motif_name in MOTIF_NAMES:
            has_match = False
            for rec in species_recs[sp]:
                if motif_matches_recognition(motif_name, rec):
                    has_match = True
                    break
            row[motif_name] = 1 if has_match else 0
        row['has_any_rm'] = 1  # species is in REBASE
        rows.append(row)

    return pd.DataFrame(rows)


# ══════════════════════════════════════════
# Genome motif density
# ══════════════════════════════════════════

def load_genome(path):
    """Load FASTA -> concatenated sequence and total length."""
    sequences = []
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as f:
        seq_parts = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq_parts:
                    sequences.append(''.join(seq_parts))
                seq_parts = []
            else:
                seq_parts.append(line.upper())
        if seq_parts:
            sequences.append(''.join(seq_parts))
    full_seq = ''.join(sequences)
    return full_seq


def count_motif_both_strands(seq, motif):
    """Count motif occurrences on both strands."""
    motif_up = motif.upper()
    has_degenerate = any(b not in 'ACGT' for b in motif_up)

    if has_degenerate:
        pattern = iupac_to_regex(motif_up)
        rc_pattern = iupac_to_regex(reverse_complement(motif_up))
        count = len(re.findall(f'(?={pattern})', seq))
        if rc_pattern != pattern:
            count += len(re.findall(f'(?={rc_pattern})', seq))
    else:
        rc = reverse_complement(motif_up)
        count = seq.count(motif_up)
        if rc != motif_up:
            count += seq.count(rc)
    return count


def compute_gc_content(seq):
    gc = seq.count('G') + seq.count('C')
    total = len(seq)
    return gc / total if total > 0 else 0


def expected_motif_count(seq_len, gc, motif):
    """Compute expected motif count based on GC content."""
    freq = {'A': (1 - gc) / 2, 'T': (1 - gc) / 2,
            'G': gc / 2, 'C': gc / 2}
    motif_up = motif.upper()
    prob_fwd = 1.0
    for b in motif_up:
        bases = IUPAC.get(b, {b})
        prob_fwd *= sum(freq.get(x, 0.25) for x in bases)
    rc = reverse_complement(motif_up)
    prob_rc = 1.0
    for b in rc:
        bases = IUPAC.get(b, {b})
        prob_rc *= sum(freq.get(x, 0.25) for x in bases)
    # Both strands
    prob_total = prob_fwd + (prob_rc if rc != motif_up else 0)
    return prob_total * seq_len


def process_genome(args):
    """Process a single genome for all motif counts."""
    accession, fasta_path, species = args
    try:
        seq = load_genome(fasta_path)
        size_bp = len(seq)
        size_mb = size_bp / 1e6
        gc = compute_gc_content(seq)

        result = {
            'accession': accession,
            'species': species,
            'size_mb': round(size_mb, 2),
            'gc_pct': round(gc * 100, 1),
        }

        for motif_name in MOTIF_NAMES:
            n = count_motif_both_strands(seq, motif_name)
            exp = expected_motif_count(size_bp, gc, motif_name)
            result[f'{motif_name}_count'] = n
            result[f'{motif_name}_density'] = round(n / size_mb, 1) if size_mb > 0 else 0
            result[f'{motif_name}_expected'] = round(exp, 1)
            result[f'{motif_name}_oe'] = round(n / exp, 3) if exp > 0 else 0

        return result
    except Exception as e:
        print(f'  ERROR processing {accession}: {e}')
        return None


def find_genome_fastas(genome_dir):
    """Find downloaded genome FASTA files from datasets download."""
    fastas = {}
    # datasets download extracts to ncbi_dataset/data/<accession>/
    for fna in genome_dir.rglob('*.fna'):
        # Get accession from path
        parts = fna.parts
        for p in parts:
            if p.startswith('GCF_') or p.startswith('GCA_'):
                fastas[p] = fna
                break
    return fastas


# ══════════════════════════════════════════
# Figure generation
# ══════════════════════════════════════════

def plot_panel_a(ax, conservation_df):
    """Panel A: R-M conservation bar chart."""
    n = len(conservation_df)
    x = np.arange(n)

    # Color bars by source
    colors = conservation_df['color'].tolist()
    bars = ax.bar(x, conservation_df['rate_pct'], color=colors, alpha=0.85,
                  edgecolor='white', linewidth=0.5, width=0.7, zorder=3)

    # Error bars (Wilson CI)
    ci_low = conservation_df['rate_pct'] - conservation_df['ci_low']
    ci_high = conservation_df['ci_high'] - conservation_df['rate_pct']
    ax.errorbar(x, conservation_df['rate_pct'],
                yerr=[ci_low, ci_high],
                fmt='none', ecolor='#424242', capsize=3, capthick=0.8,
                linewidth=0.8, zorder=4)

    # Value labels above bars
    for i, row in conservation_df.iterrows():
        n_match = row['n_match']
        n_total = row['n_total']
        rate = row['rate_pct']
        y_pos = rate + ci_high.iloc[i] + 1.0
        ax.text(i, y_pos, f'{n_match}/{n_total}\n({rate:.1f}%)',
                ha='center', va='bottom', fontsize=6, fontweight='bold')

    # X-axis labels: motif (mod_type) — rotated diagonally
    labels = [f"{row['motif']} ({row['mod_type']})" for _, row in conservation_df.iterrows()]
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7, rotation=35, ha='right')

    # Separator between groups
    n_this = sum(1 for m in MOTIFS if m[2] == 'This study')
    ax.axvline(n_this - 0.5, color='#BDBDBD', linestyle='--', linewidth=0.7, zorder=1)

    # Group labels — placed at very top above value labels, using figure coordinates
    y_top = ax.get_ylim()[1]
    ax.set_ylim(0, max(conservation_df['ci_high'].max() + 22, 55))
    y_label = ax.get_ylim()[1] - 2
    ax.text(n_this / 2 - 0.5, y_label, 'This study',
            ha='center', va='top', fontsize=8, fontstyle='italic', color=COL_THIS)
    mid_prior = n_this + (n - n_this) / 2 - 0.5
    ax.text(mid_prior, y_label, 'Prior studies',
            ha='center', va='top', fontsize=8, fontstyle='italic', color=COL_PRIOR)

    ax.set_ylabel('R-M system conservation (%)')
    ax.set_title('R-M conservation in REBASE (n=82 spp.)',
                 fontweight='bold', loc='left', fontsize=9)


def plot_panel_b(ax, density_df):
    """Panel B: Motif site density — box + strip with M145 percentile.

    Uses native log scale on y-axis so tick labels show actual values
    (e.g., 100, 1000, 10000) instead of log10 exponents.
    """
    if density_df.empty:
        ax.text(0.5, 0.5, 'No density data', transform=ax.transAxes,
                ha='center', va='center', fontsize=10)
        return

    # Prepare data — use raw values (not log-transformed)
    data_list = []
    positions = []
    valid_motifs = []
    for i, motif_name in enumerate(MOTIF_NAMES):
        col = f'{motif_name}_density'
        if col in density_df.columns:
            vals = density_df[col].dropna()
            vals = vals[vals > 0]
            if len(vals) > 0:
                data_list.append(vals.values)
                positions.append(i)
                valid_motifs.append(motif_name)

    colors = [MOTIFS[i][4] for i in positions]

    # Set log scale BEFORE plotting so boxplot uses log-spaced layout
    ax.set_yscale('log')

    # Box plot (thin, no outliers)
    bp = ax.boxplot(data_list, positions=positions, widths=0.5,
                    patch_artist=True, showfliers=False,
                    medianprops=dict(color='#424242', linewidth=1.2),
                    whiskerprops=dict(color='#757575', linewidth=0.7),
                    capprops=dict(color='#757575', linewidth=0.7),
                    boxprops=dict(linewidth=0.7))
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.35)

    # Strip plot (jittered points) — subsample for readability
    rng = np.random.default_rng(42)
    n_show = min(80, len(density_df))
    idx_show = rng.choice(len(density_df), n_show, replace=False)
    for pi, (pos, motif_name, color) in enumerate(zip(positions, valid_motifs, colors)):
        col = f'{motif_name}_density'
        vals = density_df.iloc[idx_show][col].dropna()
        vals = vals[vals > 0]
        jitter = rng.uniform(-0.15, 0.15, size=len(vals))
        ax.scatter(pos + jitter, vals.values, s=6, c=color,
                   alpha=0.35, edgecolors='none', zorder=3)

    # Highlight M145 with percentile annotation
    m145 = density_df[density_df['species'].str.contains('coelicolor', case=False, na=False)]
    if len(m145) > 0:
        for pi, (pos, motif_name) in enumerate(zip(positions, valid_motifs)):
            col = f'{motif_name}_density'
            if col in m145.columns:
                val = m145.iloc[0][col]
                if pd.notna(val) and val > 0:
                    all_vals = density_df[col].dropna()
                    pctl = (all_vals < val).sum() / len(all_vals) * 100
                    ax.scatter([pos], [val], marker='*', s=120,
                               c=COL_HIGHLIGHT, zorder=12, edgecolors='black',
                               linewidth=0.5)
                    # Percentile label — offset in axes fraction to avoid overlap
                    ax.annotate(f'{pctl:.0f}%',
                                xy=(pos, val),
                                xytext=(8, 6), textcoords='offset points',
                                fontsize=6, color=COL_HIGHLIGHT, fontweight='bold',
                                ha='left', va='bottom')

    ax.set_xticks(range(len(MOTIF_NAMES)))
    ax.set_xticklabels([m[0] for m in MOTIFS], fontsize=7, rotation=35, ha='right')
    ax.set_ylabel('Motif density (sites / Mb)')
    ax.set_title('Motif density across 833 species',
                 fontweight='bold', loc='left', fontsize=9)

    # Format y-axis ticks as plain integers (100, 1000, 10000, …)
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(
        lambda x, _: f'{x:,.0f}' if x >= 1 else f'{x:.1f}'))
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())

    # Custom legend — place at upper-left (high-density region has room)
    legend_items = [
        plt.Line2D([0], [0], marker='*', color='w', markerfacecolor=COL_HIGHLIGHT,
                    markersize=10, markeredgecolor='black', markeredgewidth=0.5,
                    label='M145 (percentile)'),
        mpatches.Patch(facecolor='#999999', alpha=0.35, label='833 spp. distribution'),
    ]
    leg = ax.legend(handles=legend_items, loc='upper left', frameon=True,
                    fancybox=False, edgecolor='#BDBDBD', fontsize=6.5,
                    bbox_to_anchor=(0.0, 1.0))


def plot_panel_c(ax, species_matrix):
    """Panel C: Species × motif heatmap."""
    if species_matrix.empty:
        ax.text(0.5, 0.5, 'No REBASE data', transform=ax.transAxes,
                ha='center', va='center')
        return

    # Only species with at least one motif match
    mat = species_matrix.set_index('species')[MOTIF_NAMES]
    row_sums = mat.sum(axis=1)
    mat_filtered = mat[row_sums > 0]

    if len(mat_filtered) < 2:
        # Just show all species if too few with matches
        mat_filtered = mat

    # Hierarchical clustering on rows
    if len(mat_filtered) > 2:
        try:
            dist = pdist(mat_filtered.values, metric='hamming')
            Z = linkage(dist, method='average')
            dendro = dendrogram(Z, no_plot=True)
            order = dendro['leaves']
            mat_filtered = mat_filtered.iloc[order]
        except Exception:
            pass

    # Custom colormap: 0=gray, 1=red
    cmap = ListedColormap(['#E0E0E0', '#D32F2F'])

    im = ax.imshow(mat_filtered.values, aspect='auto', cmap=cmap,
                   interpolation='nearest')

    # Truncate species names for readability
    species_labels = []
    for sp in mat_filtered.index:
        short = sp.replace('Streptomyces ', 'S. ')
        if len(short) > 25:
            short = short[:23] + '..'
        species_labels.append(short)

    ax.set_yticks(range(len(species_labels)))
    ax.set_yticklabels(species_labels, fontsize=5.5)
    ax.set_xticks(range(len(MOTIF_NAMES)))
    ax.set_xticklabels(MOTIF_NAMES, fontsize=6.5, rotation=45, ha='right')
    ax.tick_params(axis='x', pad=2)

    # Separator
    n_this = sum(1 for m in MOTIFS if m[2] == 'This study')
    ax.axvline(n_this - 0.5, color='white', linewidth=2, zorder=5)

    ax.set_title('R-M recognition match by species\n'
                 '(REBASE recognition seq. vs motif, IUPAC match)',
                 fontweight='bold', loc='left', fontsize=9)

    # Legend
    legend_items = [
        mpatches.Patch(facecolor='#D32F2F', label='Recognition match'),
        mpatches.Patch(facecolor='#E0E0E0', label='No match / no data'),
    ]
    ax.legend(handles=legend_items, loc='lower right', frameon=True,
              fancybox=False, edgecolor='#BDBDBD', fontsize=6.5)


def plot_panel_d(ax, density_df, conservation_df):
    """Panel D: O/E ratio dot chart — genomic over/under-representation."""
    if density_df.empty:
        ax.text(0.5, 0.5, 'No data', transform=ax.transAxes,
                ha='center', va='center')
        return

    # Compute O/E statistics per motif
    oe_data = []
    for i, motif_name in enumerate(MOTIF_NAMES):
        oe_col = f'{motif_name}_oe'
        if oe_col in density_df.columns:
            vals = density_df[oe_col].dropna()
            vals = vals[vals > 0]
            if len(vals) > 0:
                oe_data.append({
                    'motif': motif_name,
                    'median': np.median(vals),
                    'q25': np.percentile(vals, 25),
                    'q75': np.percentile(vals, 75),
                    'color': MOTIFS[i][4],
                    'pos': i,
                })

    if not oe_data:
        return

    # Horizontal dot chart (motifs on Y, O/E on X)
    for d in oe_data:
        # IQR whisker
        ax.plot([d['q25'], d['q75']], [d['pos'], d['pos']],
                color=d['color'], linewidth=2.5, alpha=0.5, zorder=2,
                solid_capstyle='round')
        # Median dot
        ax.scatter([d['median']], [d['pos']], s=60, c=d['color'],
                   edgecolors='white', linewidth=0.8, zorder=5)
        # Value annotation
        ax.text(d['q75'] + 0.03, d['pos'], f"{d['median']:.2f}",
                va='center', ha='left', fontsize=7, fontweight='bold',
                color=d['color'])

    # Reference line at O/E = 1
    ax.axvline(1.0, color='#9E9E9E', linestyle='--', linewidth=0.8, zorder=1)
    ax.text(1.0, len(oe_data) - 0.3, 'O/E = 1\n(random)', ha='center',
            va='bottom', fontsize=6.5, color='#757575', fontstyle='italic')

    # Shade regions
    ax.axvspan(0, 1.0, alpha=0.04, color='#1565C0', zorder=0)  # under
    ax.axvspan(1.0, ax.get_xlim()[1] if ax.get_xlim()[1] > 1 else 3.5,
               alpha=0.04, color='#D32F2F', zorder=0)  # over
    ax.text(0.55, -0.7, 'Under-\nrepresented', fontsize=6, color='#1565C0',
            ha='center', va='center', fontstyle='italic')
    ax.text(2.2, -0.7, 'Over-\nrepresented', fontsize=6, color='#D32F2F',
            ha='center', va='center', fontstyle='italic')

    ax.set_yticks([d['pos'] for d in oe_data])
    ax.set_yticklabels([d['motif'] for d in oe_data], fontsize=7.5)
    ax.set_xlabel('Observed / Expected ratio (GC-based)')
    ax.set_xlim(0.4, 3.2)
    ax.set_ylim(-1.2, len(oe_data) - 0.3)
    ax.set_title('Genomic O/E ratio (833 spp.)',
                 fontweight='bold', loc='left', fontsize=9)
    ax.invert_yaxis()


def generate_composite_figure(conservation_df, density_df, species_matrix):
    """Generate 4-panel composite figure."""
    fig = plt.figure(figsize=(220/25.4, 270/25.4))

    # Use gridspec for flexible layout
    gs = fig.add_gridspec(2, 2, hspace=0.55, wspace=0.40,
                          height_ratios=[1, 1.3])

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    # Panel labels
    for ax, label in zip([ax_a, ax_b, ax_c, ax_d], ['A', 'B', 'C', 'D']):
        ax.text(-0.12, 1.08, label, transform=ax.transAxes,
                fontsize=14, fontweight='bold', va='top', ha='left')

    plot_panel_a(ax_a, conservation_df)
    plot_panel_b(ax_b, density_df)
    plot_panel_c(ax_c, species_matrix)
    plot_panel_d(ax_d, density_df, conservation_df)

    out_base = OUTPUT_DIR / 'fig_motif_conservation_composite'
    fig.savefig(f'{out_base}.pdf')
    fig.savefig(f'{out_base}.svg')
    fig.savefig(f'{out_base}.png', dpi=300)
    plt.close()
    print(f'Saved: {out_base}.pdf/.svg/.png')

    # Also save individual panels
    for panel_name, plot_func, args in [
        ('panel_a_conservation', plot_panel_a, (conservation_df,)),
        ('panel_b_density', plot_panel_b, (density_df,)),
        ('panel_c_heatmap', plot_panel_c, (species_matrix,)),
    ]:
        fig_s, ax_s = plt.subplots(figsize=(90/25.4, 80/25.4))
        plot_func(ax_s, *args)
        fig_s.savefig(OUTPUT_DIR / f'{panel_name}.pdf')
        plt.close()


# ══════════════════════════════════════════
# Main
# ══════════════════════════════════════════

def main():
    print('=' * 60)
    print('Genus-wide Streptomyces Motif Conservation Analysis')
    print('=' * 60)

    # --- Phase 1: REBASE analysis ---
    print('\n[Phase 1] REBASE R-M system analysis')
    print(f'  Parsing {REBASE_PATH}...')
    records = parse_bairoch(REBASE_PATH)
    print(f'  Total REBASE records: {len(records)}')

    rebase_df = filter_streptomyces(records)
    print(f'  Streptomyces entries: {len(rebase_df)}')

    conservation_df = rebase_motif_conservation(rebase_df)
    conservation_df.to_csv(OUTPUT_DIR / 'rebase_motif_conservation_matrix.csv', index=False)
    print('\n  Conservation rates:')
    for _, row in conservation_df.iterrows():
        print(f"    {row['motif']:<15s} {row['n_match']:>3d}/{row['n_total']:>3d} "
              f"({row['rate_pct']:.1f}%) [{row['ci_low']:.1f}-{row['ci_high']:.1f}%]"
              f"  {row['source']}")

    species_matrix = rebase_species_motif_matrix(rebase_df)
    species_matrix.to_csv(OUTPUT_DIR / 'rebase_species_motif_matrix.csv', index=False)

    # --- Phase 2: Genome density analysis ---
    print('\n[Phase 2] Genome motif density analysis')

    # Load species list
    species_df = pd.read_csv(SPECIES_TSV, sep='\t')
    print(f'  {len(species_df)} representative species')

    # Find downloaded genomes
    genome_fastas = find_genome_fastas(GENOME_DIR)
    print(f'  Found {len(genome_fastas)} downloaded genomes')

    # Also check for M145 reference
    if 'GCF_000203835.1' not in genome_fastas and M145_GENOME.exists():
        genome_fastas['GCF_000203835.1'] = M145_GENOME

    # Build task list
    tasks = []
    for _, row in species_df.iterrows():
        acc = row['accession']
        if acc in genome_fastas:
            tasks.append((acc, genome_fastas[acc], row['binomial']))

    print(f'  Processing {len(tasks)} genomes...')

    if len(tasks) > 0:
        n_workers = min(cpu_count(), 8)
        with Pool(n_workers) as pool:
            results = pool.map(process_genome, tasks)
        results = [r for r in results if r is not None]
        density_df = pd.DataFrame(results)
        density_df.to_csv(OUTPUT_DIR / 'motif_site_density_genuswide.csv', index=False)
        print(f'  Computed density for {len(density_df)} genomes')
    else:
        print('  No genomes available yet - generating REBASE-only figure')
        density_df = pd.DataFrame()

    # --- Phase 3: Figure generation ---
    print('\n[Phase 3] Figure generation')
    generate_composite_figure(conservation_df, density_df, species_matrix)

    # --- Summary log ---
    with open(OUTPUT_DIR / 'analysis_log.txt', 'w') as f:
        f.write('Genus-wide Streptomyces Motif Conservation Analysis\n')
        f.write(f'REBASE entries: {len(rebase_df)}\n')
        f.write(f'Representative species: {len(species_df)}\n')
        f.write(f'Genomes processed: {len(density_df) if not density_df.empty else 0}\n')
        f.write(f'Motifs analyzed: {len(MOTIFS)}\n\n')
        f.write('Conservation rates:\n')
        for _, row in conservation_df.iterrows():
            f.write(f"  {row['motif']}: {row['n_match']}/{row['n_total']} "
                    f"({row['rate_pct']:.1f}%)\n")

    print('\n' + '=' * 60)
    print('Analysis complete')
    print(f'Output: {OUTPUT_DIR}')
    print('=' * 60)


if __name__ == '__main__':
    main()
