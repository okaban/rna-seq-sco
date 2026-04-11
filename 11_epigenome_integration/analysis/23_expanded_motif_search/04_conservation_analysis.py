#!/usr/bin/env python3
"""
04: Conservation Analysis for Novel Motifs
==========================================
Phase 4: Calculate O/E ratios for novel motifs in M145 genome,
estimate conservation across 833 RefSeq genomes using GC-based models,
and check REBASE coverage.
"""

import pandas as pd
import numpy as np
from collections import Counter
import re
import os

BASE_DIR = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration"
OUT_DIR = f"{BASE_DIR}/analysis/23_expanded_motif_search"
GENOME_FASTA = "/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/GCF_000203835.1_ASM20383v1_genomic.fna"
GENUSWIDE_CSV = f"{BASE_DIR}/analysis/21_genuswide_motif_conservation/motif_site_density_genuswide.csv"
REBASE_FILE = f"{BASE_DIR}/data/rebase/streptomyces_rm_systems.csv"
REBASE_CONSERVATION = f"{BASE_DIR}/analysis/21_genuswide_motif_conservation/rebase_motif_conservation_matrix.csv"

# ============================================================
# Utility functions
# ============================================================
def iupac_to_regex(seq):
    iupac = {"A": "A", "C": "C", "G": "G", "T": "T",
             "R": "[AG]", "Y": "[CT]", "S": "[CG]", "W": "[AT]",
             "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
             "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"}
    return "".join(iupac.get(ch, ch) for ch in seq.upper())


def reverse_complement(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp.get(b, "N") for b in reversed(seq.upper()))


def count_motif_both_strands(genome, motif_regex):
    """Count motif occurrences on both strands."""
    compiled = re.compile(motif_regex)
    fwd = len(compiled.findall(genome))
    rc_regex = iupac_to_regex(reverse_complement(
        motif_regex.replace("[", "").replace("]", "").replace("AG", "R").replace("CT", "Y")
        .replace("CG", "S").replace("AT", "W").replace("GT", "K").replace("AC", "M")))
    # Simplified: just count forward matches (palindromic motifs will be counted once)
    # For non-palindromic, we need both strands
    rc_compiled = re.compile(iupac_to_regex(reverse_complement(
        re.sub(r'\[([A-Z]+)\]', lambda m: m.group(1)[0], motif_regex))))
    rev = len(rc_compiled.findall(genome))
    return fwd + rev


def expected_count_gc(motif_seq, genome_size, gc_content):
    """Calculate expected count of a motif based on GC content model."""
    gc = gc_content / 100.0
    at = 1.0 - gc
    base_prob = {"A": at/2, "T": at/2, "C": gc/2, "G": gc/2}
    iupac_bases = {
        "A": ["A"], "C": ["C"], "G": ["G"], "T": ["T"],
        "R": ["A", "G"], "Y": ["C", "T"], "S": ["C", "G"], "W": ["A", "T"],
        "K": ["G", "T"], "M": ["A", "C"], "B": ["C", "G", "T"],
        "D": ["A", "G", "T"], "H": ["A", "C", "T"], "V": ["A", "C", "G"],
        "N": ["A", "C", "G", "T"],
    }

    prob = 1.0
    for ch in motif_seq.upper():
        bases = iupac_bases.get(ch, [ch])
        pos_prob = sum(base_prob[b] for b in bases)
        prob *= pos_prob

    # Both strands (approximate for non-palindromic)
    rc_motif = reverse_complement(motif_seq)
    if rc_motif == motif_seq:
        # Palindromic: only one strand
        return prob * genome_size
    else:
        return prob * genome_size * 2


def load_genome_sequence(fasta_path):
    seq_parts = []
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith(">"):
                seq_parts.append(line.strip().upper())
    return "".join(seq_parts)


# ============================================================
# Load data
# ============================================================
print("=" * 70)
print("PHASE 4: CONSERVATION ANALYSIS FOR NOVEL MOTIFS")
print("=" * 70)

genome_seq = load_genome_sequence(GENOME_FASTA)
genome_len = len(genome_seq)
gc_count = genome_seq.count("G") + genome_seq.count("C")
gc_pct = gc_count / genome_len * 100
print(f"M145 genome: {genome_len:,} bp, GC: {gc_pct:.1f}%")

# Load existing genuswide data
genuswide = pd.read_csv(GENUSWIDE_CSV)
print(f"Genuswide data: {len(genuswide)} species")
print(f"GC range: {genuswide['gc_pct'].min():.1f} - {genuswide['gc_pct'].max():.1f}%")
print(f"GC median: {genuswide['gc_pct'].median():.1f}%")


# ============================================================
# 1. Novel motif counts in M145 genome
# ============================================================
print("\n" + "=" * 70)
print("1. NOVEL MOTIF COUNTS IN M145 GENOME")
print("=" * 70)

novel_motifs = {
    # From Phase 1-2 discoveries
    "CCGKCA": "CCG[GT]CA",          # MEME-2 core
    "CGGCAACC": "CGGCAACC",          # STREME-3
    "GAACCGG": "GAACCGG",            # STREME-4
    "TGGCCGGC": "TGGCCGGC",          # 8bp palindrome
    # Original TOP 3 for comparison
    "CCGG": "CCGG",
    "AAGCCCG": "AAGCCCG",
    "GATC": "GATC",
    # Additional candidates
    "GCCG": "GCCG",
    "CCGC": "CCGC",
    "CCSGG": "CC[CG]GG",
    "GCGC": "GCGC",
}

m145_results = []
print(f"\n{'Motif':<14} {'Fwd':>8} {'Rev':>8} {'Total':>8} {'Density':>10} {'Expected':>10} {'O/E':>7}")
print("-" * 80)

for motif_name, regex in novel_motifs.items():
    compiled = re.compile(regex)
    fwd_count = len(compiled.findall(genome_seq))

    # Count on reverse complement
    rc_seq = reverse_complement(genome_seq)
    rev_count = len(compiled.findall(rc_seq))

    # For palindromic motifs, fwd == rev, so total = fwd (not fwd + rev)
    motif_clean = re.sub(r'\[([A-Z]+)\]', lambda m: m.group(1)[0], motif_name)
    is_palindrome = motif_clean == reverse_complement(motif_clean)

    if is_palindrome:
        total = fwd_count  # Same on both strands
    else:
        total = fwd_count + rev_count

    density = total / (genome_len / 1e6)
    expected = expected_count_gc(motif_name, genome_len, gc_pct)
    oe = total / expected if expected > 0 else float("inf")

    print(f"{motif_name:<14} {fwd_count:>8} {rev_count:>8} {total:>8} "
          f"{density:>9.1f} {expected:>9.1f} {oe:>7.3f}")

    m145_results.append({
        "motif": motif_name,
        "fwd_count": fwd_count,
        "rev_count": rev_count,
        "total_count": total,
        "density_per_Mb": round(density, 1),
        "expected_gc_model": round(expected, 1),
        "O_E_ratio": round(oe, 3),
        "is_palindrome": is_palindrome,
    })

m145_df = pd.DataFrame(m145_results)
m145_df.to_csv(f"{OUT_DIR}/m145_novel_motif_oe.csv", index=False)


# ============================================================
# 2. Expected counts across 833 genomes (GC-based model)
# ============================================================
print("\n" + "=" * 70)
print("2. PREDICTED O/E RANGE ACROSS 833 Streptomyces GENOMES")
print("=" * 70)
print("(Using GC-content model; actual genome sequences needed for precise O/E)")

genus_oe_results = []

for motif_name in novel_motifs:
    expected_counts = []
    for _, row in genuswide.iterrows():
        exp = expected_count_gc(motif_name, row["size_mb"] * 1e6, row["gc_pct"])
        expected_counts.append(exp)

    expected_arr = np.array(expected_counts)

    genus_oe_results.append({
        "motif": motif_name,
        "expected_median": round(np.median(expected_arr), 0),
        "expected_min": round(np.min(expected_arr), 0),
        "expected_max": round(np.max(expected_arr), 0),
        "note": "Based on GC model only; true O/E requires actual genome sequences",
    })

genus_df = pd.DataFrame(genus_oe_results)
genus_df.to_csv(f"{OUT_DIR}/genus_expected_counts_gc_model.csv", index=False)

print(f"\n{'Motif':<14} {'Med_expected':>12} {'Min':>10} {'Max':>10} {'M145_O/E':>10} {'Interpretation'}")
print("-" * 80)
for _, row in genus_df.iterrows():
    m145_oe = m145_df[m145_df["motif"] == row["motif"]]["O_E_ratio"].values
    oe_str = f"{m145_oe[0]:.3f}" if len(m145_oe) > 0 else "N/A"

    if len(m145_oe) > 0:
        oe = m145_oe[0]
        if oe < 0.7:
            interp = "AVOIDANCE (selection pressure)"
        elif oe > 1.5:
            interp = "MAINTAINED (positive selection)"
        elif oe > 1.3:
            interp = "Slight enrichment"
        elif oe < 0.85:
            interp = "Slight underrepresentation"
        else:
            interp = "Neutral"
    else:
        interp = "N/A"

    print(f"{row['motif']:<14} {row['expected_median']:>12.0f} {row['expected_min']:>10.0f} "
          f"{row['expected_max']:>10.0f} {oe_str:>10} {interp}")


# ============================================================
# 3. REBASE check for novel motifs
# ============================================================
print("\n" + "=" * 70)
print("3. REBASE CHECK FOR NOVEL MOTIFS")
print("=" * 70)

import csv
rebase_motifs = set()
rebase_data = {}
with open(REBASE_FILE) as f:
    reader = csv.DictReader(f)
    for row in reader:
        seq = row.get("recognition_seq", "").strip()
        if seq and seq != "?" and len(seq) >= 3:
            rebase_motifs.add(seq)
            if seq not in rebase_data:
                rebase_data[seq] = []
            rebase_data[seq].append(row.get("enzyme_name", ""))

print(f"\nTotal REBASE Streptomyces motifs: {len(rebase_motifs)}")

novel_check = ["CCGKCA", "CGGCAACC", "GAACCGG", "TGGCCGGC", "GCCG", "CCGC", "CCSGG"]
print(f"\nNovel motif REBASE presence check:")
print(f"{'Motif':<14} {'In REBASE?':>12} {'Matching REBASE motifs'}")
print("-" * 60)

for motif in novel_check:
    # Check exact match
    exact = motif in rebase_motifs
    # Check if any REBASE motif contains or is contained in this motif
    containing = [m for m in rebase_motifs if motif in m or m in motif]
    # Check if motif is a subsequence of any REBASE motif
    overlapping = [m for m in rebase_motifs
                   if any(motif[i:i+4] in m for i in range(len(motif)-3))]

    if exact:
        enzymes = ", ".join(rebase_data.get(motif, [])[:3])
        print(f"{motif:<14} {'YES':>12} {enzymes}")
    elif containing:
        print(f"{motif:<14} {'PARTIAL':>12} Contains/in: {', '.join(containing[:3])}")
    elif overlapping:
        print(f"{motif:<14} {'OVERLAP':>12} 4bp overlap with: {', '.join(overlapping[:3])}")
    else:
        print(f"{motif:<14} {'NO':>12} *** NOVEL (not in REBASE) ***")


# ============================================================
# 4. Existing genuswide REBASE conservation for comparison
# ============================================================
print("\n" + "=" * 70)
print("4. EXISTING REBASE CONSERVATION RATES (for reference)")
print("=" * 70)

if os.path.exists(REBASE_CONSERVATION):
    rebase_cons = pd.read_csv(REBASE_CONSERVATION)
    print(f"\n{'Motif':<14} {'Conservation rate'}")
    print("-" * 40)
    for _, row in rebase_cons.iterrows():
        print(f"{row.iloc[0]:<14} {row.iloc[1] if len(row) > 1 else 'N/A'}")
else:
    print("REBASE conservation matrix not found in expected format")


# ============================================================
# 5. COMPREHENSIVE SUMMARY TABLE
# ============================================================
print("\n" + "=" * 70)
print("5. COMPREHENSIVE MOTIF CHARACTERIZATION TABLE")
print("=" * 70)

summary = [
    {"motif": "CCGG/GGCCGG/TGGCCGGC", "mod": "4mC", "sites": "2034 (75.9%)",
     "m145_oe": "1.03", "rebase": "18/82 (22%)", "rebase_novel": "No",
     "candidate_mtase": "SC_RS19770/SC_RS36410", "status": "ESTABLISHED"},

    {"motif": "AAGCCCG", "mod": "4mC+6mA", "sites": "973+360 (DUAL)",
     "m145_oe": "0.65 (avoidance)", "rebase": "0/82 (0%)", "rebase_novel": "YES",
     "candidate_mtase": "SC_RS17645", "status": "*** NOVEL DUAL ***"},

    {"motif": "CCGKCA", "mod": "6mA", "sites": "153 (4.8%)",
     "m145_oe": m145_df[m145_df["motif"]=="CCGKCA"]["O_E_ratio"].values[0] if len(m145_df[m145_df["motif"]=="CCGKCA"]) > 0 else "N/A",
     "rebase": "Not in REBASE", "rebase_novel": "YES",
     "candidate_mtase": "SC_RS28835/SC_RS35335 (BREX-2)", "status": "*** NEW CANDIDATE ***"},

    {"motif": "GATC", "mod": "6mA", "sites": "37 (1.2%)",
     "m145_oe": "1.99 (maintained)", "rebase": "11/82 (13.4%)", "rebase_novel": "No",
     "candidate_mtase": "Unknown", "status": "ESTABLISHED"},

    {"motif": "CGGCAACC", "mod": "6mA", "sites": "57 (1.8%)",
     "m145_oe": m145_df[m145_df["motif"]=="CGGCAACC"]["O_E_ratio"].values[0] if len(m145_df[m145_df["motif"]=="CGGCAACC"]) > 0 else "N/A",
     "rebase": "Not in REBASE", "rebase_novel": "YES",
     "candidate_mtase": "Unknown orphan", "status": "NEW CANDIDATE"},

    {"motif": "GAACCGG", "mod": "6mA", "sites": "61 (1.9%)",
     "m145_oe": m145_df[m145_df["motif"]=="GAACCGG"]["O_E_ratio"].values[0] if len(m145_df[m145_df["motif"]=="GAACCGG"]) > 0 else "N/A",
     "rebase": "Not in REBASE", "rebase_novel": "YES",
     "candidate_mtase": "Unknown orphan", "status": "NEW CANDIDATE"},
]

print(f"\n{'Motif':<22} {'Mod':<8} {'Sites':<18} {'M145 O/E':>10} {'REBASE':>15} {'Novel?':>7} {'MTase':>25}")
print("-" * 120)
for s in summary:
    oe_str = f"{s['m145_oe']}" if isinstance(s['m145_oe'], str) else f"{s['m145_oe']:.3f}"
    print(f"{s['motif']:<22} {s['mod']:<8} {s['sites']:<18} {oe_str:>10} "
          f"{s['rebase']:>15} {s['rebase_novel']:>7} {s['candidate_mtase']:>25}")

# Save comprehensive summary
summary_df = pd.DataFrame(summary)
summary_df.to_csv(f"{OUT_DIR}/comprehensive_motif_summary.csv", index=False)

print(f"\nAll results saved to: {OUT_DIR}/")
print("Done.")
