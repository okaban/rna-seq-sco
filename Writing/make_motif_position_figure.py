#!/usr/bin/env python3
"""
Per-position methylation probability figure for S. coelicolor A3(2) M145.
mBio publication-quality version — double-column (180 mm x 80 mm), side-by-side panels.

Panel A: GCCGGC (m4C) — position 2 (+strand) and position 3 (-strand) are methylated
Panel B: AAGCCCG (6mA) — position 2 is the primary methylated A, position 1 secondary

For each motif position i, the methylation probability is:
    prob[i] = mean(frequency[i] across all N motif occurrences)
where frequency = 0 if the position is not in the census (below detection).

Census frequency column: stored as 0-100 (percent), converted to 0-1 at load time.
Timepoint: T1 only (exponential growth).

Verified lookup formula: (p + i + 1, strand), where p is 0-based genome start of
the motif occurrence and i is 0-indexed position within the motif.
"""

import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import re
import os

# ── paths ──────────────────────────────────────────────────────────────────────
BASE = '/sessions/fervent-tender-euler/mnt/rna-seq'
CENSUS_4MC = f'{BASE}/11_epigenome_integration/analysis/23_expanded_motif_search/4mC_final_census.csv'
CENSUS_6MA = f'{BASE}/11_epigenome_integration/analysis/23_expanded_motif_search/6mA_final_census.csv'
GENOME_FA  = f'{BASE}/11_epigenome_integration/data/NC_003888.3.fna'
OUT_DIR    = f'{BASE}/Writing'
os.makedirs(OUT_DIR, exist_ok=True)

TIMEPOINT = 'T1'

# ── load genome ────────────────────────────────────────────────────────────────
print("Loading genome...")
genome = {}
with open(GENOME_FA) as fh:
    chrom = None
    seqs = []
    for line in fh:
        line = line.rstrip()
        if line.startswith('>'):
            if chrom:
                genome[chrom] = ''.join(seqs).upper()
            chrom = line[1:].split()[0]
            seqs = []
        else:
            seqs.append(line)
    if chrom:
        genome[chrom] = ''.join(seqs).upper()

CHROM = 'NC_003888.3'
seq = genome[CHROM]
print(f"  Genome {CHROM}: {len(seq):,} bp")

# ── load census ────────────────────────────────────────────────────────────────
# Frequency values are 0-100; divide by 100 -> 0-1 fraction
print(f"Loading census files (timepoint={TIMEPOINT})...")

look4 = {}   # (1-based pos, strand) -> frequency 0-1
with open(CENSUS_4MC) as fh:
    for row in csv.DictReader(fh):
        if row['timepoint'] != TIMEPOINT:
            continue
        fm = row.get('final_motif', '')
        if fm not in ('TGGCCGGC', 'GCCGGC'):
            continue
        pos    = int(row['position'])
        strand = row['strand']
        freq   = float(row['frequency']) / 100.0   # convert % -> fraction
        look4[(pos, strand)] = freq

look6 = {}   # (1-based pos, strand) -> frequency 0-1
with open(CENSUS_6MA) as fh:
    for row in csv.DictReader(fh):
        if row['timepoint'] != TIMEPOINT:
            continue
        if row.get('final_motif', '') != 'AAGCCCG':
            continue
        if row.get('center_base', '') != 'A':
            continue   # skip non-adenine entries (artifact detections)
        pos    = int(row['position'])
        strand = row['strand']
        freq   = float(row['frequency']) / 100.0
        look6[(pos, strand)] = freq

print(f"  4mC GCCGGC (T1): {len(look4):,} entries")
print(f"  6mA AAGCCCG (T1): {len(look6):,} entries")

# ── reverse complement ─────────────────────────────────────────────────────────
def rc(s):
    comp = str.maketrans('ACGTacgt', 'TGCAtgca')
    return s.translate(comp)[::-1]

# ── scan genome for motif occurrences ─────────────────────────────────────────
MOTIF_4MC = 'GCCGGC'   # palindrome (RC = GCCGGC)
MOTIF_6MA = 'AAGCCCG'  # RC = CGGGCTT

def find_all(sequence, motif):
    positions = []
    pat = re.compile(motif, re.IGNORECASE)
    for m in pat.finditer(sequence):
        positions.append(m.start())
    return positions

print("Scanning genome...")
plus_gccggc     = find_all(seq, MOTIF_4MC)
plus_aagcccg    = find_all(seq, MOTIF_6MA)
plus_aagcccg_rc = find_all(seq, rc(MOTIF_6MA))   # CGGGCTT -> minus-strand AAGCCCG

print(f"  GCCGGC: {len(plus_gccggc):,} occurrences (palindrome, same on both strands)")
print(f"  AAGCCCG: {len(plus_aagcccg):,} + strand, {len(plus_aagcccg_rc):,} - strand")

# ── per-position methylation probability ──────────────────────────────────────
n4 = len(MOTIF_4MC)   # 6
n6 = len(MOTIF_6MA)   # 7
N_gccggc = len(plus_gccggc)

# GCCGGC
sum_plus_g  = np.zeros(n4)
sum_minus_g = np.zeros(n4)
for p in plus_gccggc:
    for i in range(n4):
        gp = p + i + 1   # 1-based
        sum_plus_g[i]  += look4.get((gp, '+'), 0.0)
        sum_minus_g[i] += look4.get((gp, '-'), 0.0)

prob_plus_g  = sum_plus_g  / N_gccggc   # 0-1
prob_minus_g = sum_minus_g / N_gccggc   # 0-1

print("\nGCCGGC per-position methylation probability:")
for i, base in enumerate(MOTIF_4MC):
    n_plus  = sum(1 for p in plus_gccggc if look4.get((p+i+1,'+'), 0) > 0)
    n_minus = sum(1 for p in plus_gccggc if look4.get((p+i+1,'-'), 0) > 0)
    print(f"  pos{i+1} ({base}): +strand={prob_plus_g[i]*100:.2f}% ({n_plus} sites)  "
          f"-strand={prob_minus_g[i]*100:.2f}% ({n_minus} sites)")

# AAGCCCG
N_aagcccg = len(plus_aagcccg) + len(plus_aagcccg_rc)
sum_aagcccg = np.zeros(n6)

for p in plus_aagcccg:
    for j in range(n6):
        gp = p + j + 1
        sum_aagcccg[j] += look6.get((gp, '+'), 0.0)

for q in plus_aagcccg_rc:
    for j in range(n6):
        gp = q + (n6 - 1 - j) + 1   # 1-based fwd coord of AAGCCCG pos j on - strand
        sum_aagcccg[j] += look6.get((gp, '-'), 0.0)

prob_aagcccg = sum_aagcccg / N_aagcccg   # 0-1

print("\nAAGCCCG per-position methylation probability:")
for i, base in enumerate(MOTIF_6MA):
    print(f"  pos{i+1} ({base}): {prob_aagcccg[i]*100:.2f}%")

# ── figure (mBio publication style) ───────────────────────────────────────────
print("\nGenerating mBio publication-quality figure...")

import matplotlib.font_manager as fm
_avail  = {f.name for f in fm.fontManager.ttflist}
_pref   = ['Arial', 'Helvetica', 'Liberation Sans', 'FreeSans', 'DejaVu Sans']
_font   = next((f for f in _pref if f in _avail), 'DejaVu Sans')
print(f"  Font selected: {_font}")

plt.rcParams.update({
    'font.family':       'sans-serif',
    'font.sans-serif':   _pref,
    'font.size':         9,
    'axes.linewidth':    0.8,
    'axes.spines.top':   False,
    'axes.spines.right': False,
    'xtick.major.size':  3.5,
    'ytick.major.size':  3.5,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'xtick.direction':   'out',
    'ytick.direction':   'out',
    'pdf.fonttype':      42,   # TrueType in PDF -- editable in Illustrator
    'ps.fonttype':       42,
})

# ── Colour palette (print-safe, colour-blind friendly) ────────────────────────
C_PLUS    = '#D9541E'   # orange-red  -- + strand m4C  (pos 2)
C_MINUS   = '#8B2020'   # dark crimson -- - strand m4C  (pos 3)
C_6MA_PRI = '#1A6FAD'   # steel blue  -- primary 6mA   (pos 2)
C_6MA_SEC = '#80BBD9'   # sky blue    -- secondary 6mA  (pos 1)
C_ZERO    = '#E0E0E0'   # light gray  -- unmethylated

# ── Canvas: 180 mm x 80 mm  (mBio double-column standard) ────────────────────
fig_w  = 180 / 25.4     # 7.087 in
fig_h  =  80 / 25.4     # 3.150 in
YMAX   = 12.0
BWIDTH = 0.60

fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(fig_w, fig_h))
fig.subplots_adjust(left=0.11, right=0.97, bottom=0.22, top=0.87, wspace=0.52)

# ─────────────────────────────────────────────────────────────────────────────
# Panel A  --  GCCGGC (m4C)
# ─────────────────────────────────────────────────────────────────────────────
x_a   = np.arange(n4)
vals_a = (prob_plus_g + prob_minus_g) * 100   # total per position

col_a = []
for i in range(n4):
    if i == 1:      col_a.append(C_PLUS)    # position 2 -- + strand m4C
    elif i == 2:    col_a.append(C_MINUS)   # position 3 -- - strand m4C
    else:           col_a.append(C_ZERO)

ax_a.bar(x_a, vals_a, BWIDTH, color=col_a, edgecolor='none', zorder=3)

# Annotations above the two significant bars (staggered to avoid overlap)
ANN_FS = 6.5
ax_a.annotate('+strand m⁴C',
              xy=(1, vals_a[1]), xytext=(0.5, vals_a[1] + 3.2),
              ha='center', va='bottom', fontsize=ANN_FS,
              color=C_PLUS, fontweight='bold',
              arrowprops=dict(arrowstyle='->', color=C_PLUS, lw=0.7,
                              shrinkA=0, shrinkB=2))
ax_a.annotate('−strand m⁴C\n(palindrome)',
              xy=(2, vals_a[2]), xytext=(3.5, vals_a[2] + 3.2),
              ha='center', va='bottom', fontsize=ANN_FS,
              color=C_MINUS, fontweight='bold',
              arrowprops=dict(arrowstyle='->', color=C_MINUS, lw=0.7,
                              shrinkA=0, shrinkB=2))

# Axis labels: "1\nG", "2\nC", ...
tick_labels_a = [f'{i+1}\n{b}' for i, b in enumerate(MOTIF_4MC)]
ax_a.set_xticks(x_a)
ax_a.set_xticklabels(tick_labels_a, fontsize=8.5)
ax_a.set_ylim(0, YMAX)
ax_a.set_xlim(-0.55, n4 - 0.45)
ax_a.set_xlabel('Motif position', fontsize=9, labelpad=2)
ax_a.set_ylabel('Mean methylation\nprobability (%)', fontsize=9, labelpad=3)
ax_a.set_title('GCCGGC (m⁴C)', fontsize=10, fontweight='bold', pad=5)
ax_a.tick_params(axis='both', labelsize=8.5)
ax_a.yaxis.set_major_locator(plt.MultipleLocator(4))
ax_a.yaxis.set_minor_locator(plt.MultipleLocator(2))

# Panel label A -- 12 pt bold, upper-left outside axes
ax_a.text(-0.24, 1.10, 'A', transform=ax_a.transAxes,
          fontsize=12, fontweight='bold', va='top', ha='left')

# ─────────────────────────────────────────────────────────────────────────────
# Panel B  --  AAGCCCG (6mA)
# ─────────────────────────────────────────────────────────────────────────────
x_b   = np.arange(n6)
vals_b = prob_aagcccg * 100

col_b = []
for i in range(n6):
    if i == 1:      col_b.append(C_6MA_PRI)   # position 2 -- primary 6mA
    elif i == 0:    col_b.append(C_6MA_SEC)   # position 1 -- secondary 6mA
    else:           col_b.append(C_ZERO)

ax_b.bar(x_b, vals_b, BWIDTH, color=col_b, edgecolor='none', zorder=3)

tick_labels_b = [f'{i+1}\n{b}' for i, b in enumerate(MOTIF_6MA)]
ax_b.set_xticks(x_b)
ax_b.set_xticklabels(tick_labels_b, fontsize=8.5)
ax_b.set_ylim(0, YMAX)
ax_b.set_xlim(-0.55, n6 - 0.45)
ax_b.set_xlabel('Motif position', fontsize=9, labelpad=2)
ax_b.set_ylabel('Mean methylation\nprobability (%)', fontsize=9, labelpad=3)
ax_b.set_title('AAGCCCG (⁶mA)', fontsize=10, fontweight='bold', pad=5)
ax_b.tick_params(axis='both', labelsize=8.5)
ax_b.yaxis.set_major_locator(plt.MultipleLocator(4))
ax_b.yaxis.set_minor_locator(plt.MultipleLocator(2))

# Legend
legend_b = [
    mpatches.Patch(color=C_6MA_PRI, label='Primary ⁶mA (pos 2)'),
    mpatches.Patch(color=C_6MA_SEC, label='Secondary ⁶mA (pos 1)'),
]
ax_b.legend(handles=legend_b, fontsize=7.0, frameon=False,
            loc='upper right', handlelength=1.0, handletextpad=0.4)

# Panel label B
ax_b.text(-0.24, 1.10, 'B', transform=ax_b.transAxes,
          fontsize=12, fontweight='bold', va='top', ha='left')

# ── save ───────────────────────────────────────────────────────────────────────
out_pdf = f'{OUT_DIR}/Figure_motif_methylation_position.pdf'
out_png = f'{OUT_DIR}/Figure_motif_methylation_position.png'

fig.savefig(out_pdf, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
print(f"\nSaved PDF : {out_pdf}")
print(f"Saved PNG : {out_png}")

print("\n=== Final per-position statistics ===")
for i, base in enumerate(MOTIF_4MC):
    if prob_plus_g[i] + prob_minus_g[i] > 0:
        print(f"  GCCGGC  pos{i+1} ({base}):  "
              f"+strand={prob_plus_g[i]*100:.2f}%  "
              f"-strand={prob_minus_g[i]*100:.2f}%")
for i, base in enumerate(MOTIF_6MA):
    if prob_aagcccg[i] > 0:
        print(f"  AAGCCCG pos{i+1} ({base}):  {prob_aagcccg[i]*100:.2f}%")
print("Done.")
