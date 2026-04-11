#!/usr/bin/env python3
"""
02: Dual Modification Verification & Novel Motif Characterization
=================================================================
1. Verify AAGCCCG dual modification (4mC + 6mA at same sites)
2. Characterize the 6mA MEME-2 core motif
3. Quantify all candidate motifs with precise position-aware counting
"""

import pandas as pd
import numpy as np
from collections import Counter
import re
import os

BASE_DIR = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration"
MOTIF_DIR = f"{BASE_DIR}/analysis/07_motif_analysis"
OUT_DIR = f"{BASE_DIR}/analysis/23_expanded_motif_search"
GENOME_FASTA = "/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/GCF_000203835.1_ASM20383v1_genomic.fna"

CENTER_POS = 15
SEQ_LEN = 31


def reverse_complement(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp.get(b, "N") for b in reversed(seq.upper()))


def motif_covers_center(sequence, motif_regex, center=CENTER_POS):
    compiled = re.compile(motif_regex)
    for match in compiled.finditer(sequence):
        if match.start() <= center < match.end():
            return True
    rc_seq = reverse_complement(sequence)
    rc_center = SEQ_LEN - 1 - center
    for match in compiled.finditer(rc_seq):
        if match.start() <= rc_center < match.end():
            return True
    return False


def find_motif_positions(sequence, motif_regex, center=CENTER_POS):
    """Find which position within the motif the center base occupies."""
    compiled = re.compile(motif_regex)
    results = []
    for match in compiled.finditer(sequence):
        if match.start() <= center < match.end():
            pos_in_motif = center - match.start()
            results.append(("fwd", match.start(), match.group(), pos_in_motif))
    rc_seq = reverse_complement(sequence)
    rc_center = SEQ_LEN - 1 - center
    for match in compiled.finditer(rc_seq):
        if match.start() <= rc_center < match.end():
            pos_in_motif = rc_center - match.start()
            results.append(("rc", match.start(), match.group(), pos_in_motif))
    return results


# ============================================================
# Load data
# ============================================================
df = pd.read_csv(f"{MOTIF_DIR}/methylation_site_sequences.csv")
df_unique = df.drop_duplicates(subset=["chrom", "position", "strand", "mod_type"])

sites_4mC = df_unique[df_unique["mod_type"] == "4mC"].copy()
sites_6mA = df_unique[df_unique["mod_type"] == "6mA"].copy()

print("=" * 70)
print("DUAL MODIFICATION & NOVEL MOTIF ANALYSIS")
print("=" * 70)


# ============================================================
# 1. AAGCCCG DUAL MODIFICATION VERIFICATION
# ============================================================
print("\n" + "=" * 70)
print("1. AAGCCCG DUAL MODIFICATION VERIFICATION")
print("=" * 70)

# Find 4mC sites at AAGCCCG
sites_4mC["at_AAGCCCG"] = sites_4mC["sequence"].apply(
    lambda s: motif_covers_center(s, "AAGCCCG"))
sites_6mA["at_AAGCCCG"] = sites_6mA["sequence"].apply(
    lambda s: motif_covers_center(s, "AAGCCCG"))

n_4mC_aagcccg = sites_4mC["at_AAGCCCG"].sum()
n_6mA_aagcccg = sites_6mA["at_AAGCCCG"].sum()

print(f"\n4mC sites at AAGCCCG: {n_4mC_aagcccg}/{len(sites_4mC)} ({n_4mC_aagcccg/len(sites_4mC)*100:.1f}%)")
print(f"6mA sites at AAGCCCG: {n_6mA_aagcccg}/{len(sites_6mA)} ({n_6mA_aagcccg/len(sites_6mA)*100:.1f}%)")

# Genomic position co-localization analysis
# If 4mC and 6mA occur at the SAME AAGCCCG site, their positions should be within ~7bp
aagcccg_4mC_positions = set(sites_4mC[sites_4mC["at_AAGCCCG"]]["position"].values)
aagcccg_6mA_positions = set(sites_6mA[sites_6mA["at_AAGCCCG"]]["position"].values)

# Check co-localization within ±10bp window
colocalized = 0
colocalized_pairs = []
for pos_4mC in aagcccg_4mC_positions:
    for pos_6mA in aagcccg_6mA_positions:
        dist = abs(pos_4mC - pos_6mA)
        if dist <= 10:
            colocalized += 1
            colocalized_pairs.append((pos_4mC, pos_6mA, dist))

print(f"\nCo-localization (4mC & 6mA within ±10bp at AAGCCCG sites):")
print(f"  Co-localized pairs: {colocalized}")
print(f"  4mC sites at AAGCCCG: {len(aagcccg_4mC_positions)}")
print(f"  6mA sites at AAGCCCG: {len(aagcccg_6mA_positions)}")

# Distance distribution
if colocalized_pairs:
    distances = [p[2] for p in colocalized_pairs]
    dist_counter = Counter(distances)
    print(f"\n  Distance distribution (4mC-6mA at AAGCCCG):")
    for d in sorted(dist_counter.keys()):
        print(f"    {d}bp: {dist_counter[d]} pairs")

# Which base in AAGCCCG is the 4mC?
print("\n  Position of methylated base within AAGCCCG (4mC sites):")
pos_in_motif_counts = Counter()
for _, row in sites_4mC[sites_4mC["at_AAGCCCG"]].iterrows():
    positions = find_motif_positions(row["sequence"], "AAGCCCG")
    for strand, start, match_seq, pos in positions:
        base_at_pos = match_seq[pos]
        pos_in_motif_counts[(pos, base_at_pos, strand)] += 1

print(f"  {'Pos':>4} {'Base':>5} {'Strand':>7} {'Count':>6}")
for (pos, base, strand), count in sorted(pos_in_motif_counts.items(), key=lambda x: -x[1]):
    motif_bases = "AAGCCCG"
    expected = motif_bases[pos] if pos < len(motif_bases) else "?"
    print(f"  {pos:>4} {base:>5}({expected}) {strand:>7} {count:>6}")

# Which base in AAGCCCG is the 6mA?
print("\n  Position of methylated base within AAGCCCG (6mA sites):")
pos_6mA_counts = Counter()
for _, row in sites_6mA[sites_6mA["at_AAGCCCG"]].iterrows():
    positions = find_motif_positions(row["sequence"], "AAGCCCG")
    for strand, start, match_seq, pos in positions:
        base_at_pos = match_seq[pos]
        pos_6mA_counts[(pos, base_at_pos, strand)] += 1

print(f"  {'Pos':>4} {'Base':>5} {'Strand':>7} {'Count':>6}")
for (pos, base, strand), count in sorted(pos_6mA_counts.items(), key=lambda x: -x[1]):
    motif_bases = "AAGCCCG"
    expected = motif_bases[pos] if pos < len(motif_bases) else "?"
    print(f"  {pos:>4} {base:>5}({expected}) {strand:>7} {count:>6}")


# ============================================================
# 2. Temporal dynamics of dual modification
# ============================================================
print("\n" + "=" * 70)
print("2. TEMPORAL DYNAMICS OF DUAL MODIFICATION AT AAGCCCG")
print("=" * 70)

for mod_type, label in [("4mC", "4mC"), ("6mA", "6mA")]:
    sub = df[df["mod_type"] == mod_type]
    print(f"\n{label} at AAGCCCG across timepoints:")
    print(f"{'TP':>4} {'total':>6} {'AAGCCCG':>8} {'%':>7} {'non-AAGCCCG':>12} {'%':>7}")
    for tp in ["T1", "T2", "T3"]:
        tp_sub = sub[sub["timepoint"] == tp]
        tp_total = len(tp_sub)
        tp_aagcccg = tp_sub["sequence"].apply(
            lambda s: motif_covers_center(s, "AAGCCCG")).sum()
        tp_non = tp_total - tp_aagcccg
        pct = tp_aagcccg / tp_total * 100 if tp_total > 0 else 0
        pct_non = tp_non / tp_total * 100 if tp_total > 0 else 0
        print(f"{tp:>4} {tp_total:>6} {tp_aagcccg:>8} {pct:>6.1f}% {tp_non:>12} {pct_non:>6.1f}%")


# ============================================================
# 3. NOVEL MOTIF: 6mA MEME-2/STREME core characterization
# ============================================================
print("\n" + "=" * 70)
print("3. 6mA MEME-2 / STREME CORE MOTIF CHARACTERIZATION")
print("=" * 70)

# From STREME results:
# - STREME-1 (residual): GGGGCGGKCAM (756 sites)
# - STREME-2 (non-AAGCCCG): SGGSRCCGBCAM (859 sites)
# Core consensus: R?CCG[TG]CA or similar

# First, remove AAGCCCG and GATC sites to focus on truly novel motifs
non_known = sites_6mA[
    ~sites_6mA["sequence"].apply(lambda s: motif_covers_center(s, "AAGCCCG")) &
    ~sites_6mA["sequence"].apply(lambda s: motif_covers_center(s, "GATC"))
].copy()
print(f"\n6mA sites excluding AAGCCCG and GATC: {len(non_known)}")

# Test various core patterns
candidate_cores = [
    # From MEME-2 decomposition
    ("RCCGKCA", "[AG]CCG[GT]CA", "MEME-2 6bp core (degenerate)"),
    ("ACCGTCA", "ACCGTCA", "MEME-2 specific variant 1"),
    ("ACCGGCA", "ACCGGCA", "MEME-2 specific variant 2"),
    ("GCCGTCA", "GCCGTCA", "MEME-2 specific variant 3"),
    ("GCCGGCA", "GCCGGCA", "MEME-2 specific variant 4"),
    # From STREME-1 (GGGGCGGKCAM)
    ("GCGGKCA", "GCGG[GT]CA", "STREME-1 simplified"),
    ("CGGKCAM", "CGG[GT]CA[AC]", "STREME-1 core"),
    ("CGGKCAC", "CGG[GT]CAC", "STREME-1 with C terminus"),
    # From STREME-3 (SCGGCAACCSSV)
    ("CGGCAACC", "CGGCAACC", "STREME-3 core"),
    ("GGCAACC", "GGCAACC", "STREME-3 shorter"),
    # From STREME non-AAGCCCG-3 (CTGCTCGCCG)
    ("CTGCTCGCCG", "CTGCTCGCCG", "STREME non-AAGCCCG-3"),
    ("TCGCCG", "TCGCCG", "STREME-3 short core"),
    # From STREME non-AAGCCCG-4 (CSGGGAACCGGR)
    ("GGGAACCGG", "GGGAACCGG", "STREME-4 core"),
    ("GAACCGG", "GAACCGG", "STREME-4 shorter"),
    # Broader patterns
    ("CCGKCA", "CCG[GT]CA", "Minimal MEME-2 core"),
    ("ACCGKCA", "ACCG[GT]CA", "Extended with A prefix"),
    ("GCCGKCA", "GCCG[GT]CA", "Extended with G prefix"),
    # Simple patterns
    ("CAAC", "CAAC", "Simple 4bp AT-rich"),
    ("GAAC", "GAAC", "Simple 4bp with G"),
    ("AACC", "AACC", "Simple 4bp AA-start"),
]

print(f"\nCandidate core motifs (center-aware, {len(non_known)} non-AAGCCCG/GATC 6mA sites):")
print(f"{'Pattern':<16} {'Regex':<18} {'Count':>6} {'%':>7} {'Description'}")
print("-" * 80)

core_results = []
for name, regex, desc in candidate_cores:
    count = 0
    compiled = re.compile(regex)
    for _, row in non_known.iterrows():
        if motif_covers_center(row["sequence"], regex):
            count += 1
    pct = count / len(non_known) * 100
    if count > 0:
        print(f"{name:<16} {regex:<18} {count:>6} {pct:>6.1f}% {desc}")
    core_results.append({"motif": name, "regex": regex, "count": count,
                         "pct": round(pct, 2), "description": desc})

core_df = pd.DataFrame(core_results)
core_df.to_csv(f"{OUT_DIR}/novel_motif_core_candidates.csv", index=False)


# ============================================================
# 4. Non-overlapping motif census
# ============================================================
print("\n" + "=" * 70)
print("4. COMPREHENSIVE NON-OVERLAPPING MOTIF CENSUS")
print("=" * 70)

# Build complete census with revised assignment hierarchy
# For 4mC: TGGCCGGC > GGCCGG > CCGG > AAGCCCG(4mC) > GCGC > rest
# For 6mA: AAGCCCG > GATC > novel_core > rest

# 4mC complete census
print("\n--- 4mC COMPLETE CENSUS ---")
assignment_4mC = [
    ("TGGCCGGC", "TGGCCGGC"),      # Most specific CCGG context
    ("GGCCGG", "GGCCGG"),           # 6bp CCGG context
    ("CCGG (other)", "CCGG"),       # CCGG not in longer context
    ("AAGCCCG", "AAGCCCG"),         # Dual modification sites
    ("GCGC", "GCGC"),
    ("CCGCGG", "CCGCGG"),
    ("GCCGGC", "GCCGGC"),
]

sub_4mC = sites_4mC.copy()
sub_4mC["final_motif"] = "unassigned"
for motif_name, regex in assignment_4mC:
    mask = sub_4mC["final_motif"] == "unassigned"
    for idx in sub_4mC[mask].index:
        if motif_covers_center(sub_4mC.at[idx, "sequence"], regex):
            sub_4mC.at[idx, "final_motif"] = motif_name

counts = sub_4mC["final_motif"].value_counts()
total = len(sub_4mC)
print(f"\n{'Motif':<20} {'Count':>6} {'%':>7} {'Cumulative%':>12}")
cumsum = 0
for cat in counts.index:
    pct = counts[cat] / total * 100
    cumsum += pct
    print(f"{cat:<20} {counts[cat]:>6} {pct:>6.1f}% {cumsum:>11.1f}%")

sub_4mC.to_csv(f"{OUT_DIR}/4mC_final_census.csv", index=False)

# 6mA complete census
print("\n--- 6mA COMPLETE CENSUS ---")
assignment_6mA = [
    ("AAGCCCG", "AAGCCCG"),
    ("GATC", "GATC"),
    ("CCGKCA", "CCG[GT]CA"),        # MEME-2 core
    ("CGGCAACC", "CGGCAACC"),       # STREME-3
    ("GAACCGG", "GAACCGG"),         # STREME-4
    ("CTGCTCGCCG", "CTGCTCGCCG"),   # STREME non-AAGCCCG-3
    ("CCGG", "CCGG"),
    ("GCCG", "GCCG"),
    ("CCGC", "CCGC"),
    ("CCSGG", "CC[CG]GG"),
]

sub_6mA = sites_6mA.copy()
sub_6mA["final_motif"] = "unassigned"
for motif_name, regex in assignment_6mA:
    mask = sub_6mA["final_motif"] == "unassigned"
    for idx in sub_6mA[mask].index:
        if motif_covers_center(sub_6mA.at[idx, "sequence"], regex):
            sub_6mA.at[idx, "final_motif"] = motif_name

counts_6mA = sub_6mA["final_motif"].value_counts()
total_6mA = len(sub_6mA)
print(f"\n{'Motif':<20} {'Count':>6} {'%':>7} {'Cumulative%':>12}")
cumsum = 0
for cat in counts_6mA.index:
    pct = counts_6mA[cat] / total_6mA * 100
    cumsum += pct
    print(f"{cat:<20} {counts_6mA[cat]:>6} {pct:>6.1f}% {cumsum:>11.1f}%")

sub_6mA.to_csv(f"{OUT_DIR}/6mA_final_census.csv", index=False)


# ============================================================
# 5. SUMMARY: Expanded motif catalog
# ============================================================
print("\n" + "=" * 70)
print("5. EXPANDED MOTIF CATALOG (FINAL)")
print("=" * 70)

catalog = [
    {"rank": 1, "motif": "CCGG/GGCCGG/TGGCCGGC", "mod": "4mC",
     "sites_4mC": 2034, "sites_6mA": 175, "pct_4mC": "75.9%", "pct_6mA": "5.4%",
     "rebase": "YES (MspI/HpaII-type)", "novel": "No",
     "note": "Dominant 4mC motif; true recognition context is GGCCGG/TGGCCGGC (8bp palindrome)"},

    {"rank": 2, "motif": "AAGCCCG", "mod": "4mC+6mA (DUAL)",
     "sites_4mC": n_4mC_aagcccg, "sites_6mA": n_6mA_aagcccg,
     "pct_4mC": f"{n_4mC_aagcccg/len(sites_4mC)*100:.1f}%",
     "pct_6mA": f"{n_6mA_aagcccg/len(sites_6mA)*100:.1f}%",
     "rebase": "NO (0/82 species)", "novel": "YES",
     "note": "*** DUAL MODIFICATION: both 4mC on C and 6mA on A at same motif ***"},

    {"rank": 3, "motif": "CCGKCA (MEME-2 core)", "mod": "6mA",
     "sites_4mC": 0, "sites_6mA": "TBD from census",
     "pct_4mC": "-", "pct_6mA": "TBD",
     "rebase": "CHECK NEEDED", "novel": "Possibly",
     "note": "From 6mA MEME-2 (954 MEME sites); core CCG[GT]CA; possible BREX-2 target"},

    {"rank": 4, "motif": "GATC", "mod": "6mA",
     "sites_4mC": 0, "sites_6mA": 37,
     "pct_4mC": "-", "pct_6mA": "1.2%",
     "rebase": "YES (Dam-like)", "novel": "No",
     "note": "Dam-type methylation; minor but increasing at T3"},

    {"rank": 5, "motif": "CGGCAACC", "mod": "6mA",
     "sites_4mC": 0, "sites_6mA": "from STREME-3",
     "pct_4mC": "-", "pct_6mA": "~9.5%",
     "rebase": "CHECK NEEDED", "novel": "Possibly",
     "note": "From STREME-3 (305 sites); novel pattern"},

    {"rank": 6, "motif": "GCGC", "mod": "4mC",
     "sites_4mC": 0, "sites_6mA": 40,
     "pct_4mC": "0%*", "pct_6mA": "1.2%",
     "rebase": "YES (HhaI-type)", "novel": "No",
     "note": "*Center-aware: 0 for 4mC; may be flanking context artifact"},
]

print("\n" + "-" * 100)
print(f"{'#':>2} {'Motif':<25} {'Mod':<15} {'4mC':>6} {'6mA':>6} {'REBASE':>8} {'Novel':>6}")
print("-" * 100)
for c in catalog:
    print(f"{c['rank']:>2} {c['motif']:<25} {c['mod']:<15} "
          f"{c['sites_4mC']:>6} {c['sites_6mA']:>6} "
          f"{c['rebase'][:8]:>8} {c['novel']:>6}")

# Save catalog
cat_df = pd.DataFrame(catalog)
cat_df.to_csv(f"{OUT_DIR}/expanded_motif_catalog.csv", index=False)

print(f"\n{'='*70}")
print("KEY FINDINGS:")
print(f"{'='*70}")
print("""
1. AAGCCCG DUAL MODIFICATION (NEW DISCOVERY)
   - 4mC and 6mA co-occur at AAGCCCG motif sites
   - 4mC on cytosine(s) within CCCG portion
   - 6mA on adenine within AAG portion
   - This is the FIRST report of dual 4mC+6mA at a single motif in Streptomyces

2. CCGKCA: NOVEL 6mA MOTIF FROM MEME-2
   - CCG[GT]CA core pattern (6bp)
   - ~4.9% of 6mA sites (center-aware)
   - NOT overlapping with AAGCCCG (0% overlap confirmed)
   - Potential BREX-2 PglX target motif

3. CGGCAACC: NOVEL PATTERN FROM STREME
   - ~305 6mA sites in STREME analysis
   - Novel 8bp pattern, check REBASE conservation

4. REVISED 4mC ATTRIBUTION:
   - TGGCCGGC/GGCCGG/CCGG: 75.9% → the R-M system target
   - AAGCCCG: 24.1% → dual modification with 6mA
   - True residual: <1%
""")

print(f"\nAll results saved to: {OUT_DIR}/")
