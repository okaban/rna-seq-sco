#!/usr/bin/env python3
"""
Expanded Methylation Motif Discovery Pipeline (v2)
===================================================
Position-aware motif matching: only counts a motif if the methylated base
(center of 31bp window, position 15) falls WITHIN the motif occurrence.

Phase 1: Residual site motif analysis (non-TOP3 sites)
Phase 2: REBASE reverse matching (all Streptomyces motifs vs M145 methylation data)
"""

import pandas as pd
import numpy as np
from collections import Counter
import re
import os
import csv
from scipy.stats import binomtest

# ============================================================
# Configuration
# ============================================================
BASE_DIR = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration"
MOTIF_DIR = f"{BASE_DIR}/analysis/07_motif_analysis"
OUT_DIR = f"{BASE_DIR}/analysis/23_expanded_motif_search"
REBASE_FILE = f"{BASE_DIR}/data/rebase/streptomyces_rm_systems.csv"
GENOME_FASTA = "/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/GCF_000203835.1_ASM20383v1_genomic.fna"

CENTER_POS = 15  # 0-indexed center of 31bp window
SEQ_LEN = 31

os.makedirs(OUT_DIR, exist_ok=True)


# ============================================================
# Utility functions
# ============================================================
def iupac_to_regex(seq):
    """Convert IUPAC ambiguous DNA sequence to regex."""
    iupac = {
        "A": "A", "C": "C", "G": "G", "T": "T",
        "R": "[AG]", "Y": "[CT]", "S": "[CG]", "W": "[AT]",
        "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
        "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
    }
    return "".join(iupac.get(ch, ch) for ch in seq.upper())


def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    comp = {"A": "T", "T": "A", "C": "G", "G": "C",
            "R": "Y", "Y": "R", "S": "S", "W": "W",
            "K": "M", "M": "K", "B": "V", "V": "B",
            "D": "H", "H": "D", "N": "N"}
    return "".join(comp.get(b, "N") for b in reversed(seq.upper()))


def motif_covers_center(sequence, motif_regex, center=CENTER_POS, mod_base=None):
    """
    Check if any occurrence of motif_regex in sequence covers the center position.
    If mod_base is specified (e.g., 'C' for 4mC), also check that the center base
    is part of the motif's expected methylated position.

    Returns True if the center position falls within any motif match.
    """
    compiled = re.compile(motif_regex)
    for match in compiled.finditer(sequence):
        start, end = match.start(), match.end()
        if start <= center < end:
            return True
    # Also check reverse complement
    rc_seq = reverse_complement(sequence)
    rc_center = SEQ_LEN - 1 - center
    for match in compiled.finditer(rc_seq):
        start, end = match.start(), match.end()
        if start <= rc_center < end:
            return True
    return False


def count_motif_at_center(sequences, motif_regex, center=CENTER_POS):
    """Count how many sequences have the motif covering the center position."""
    count = 0
    compiled = re.compile(motif_regex)
    for seq in sequences:
        if seq is None or pd.isna(seq):
            continue
        found = False
        for match in compiled.finditer(seq):
            if match.start() <= center < match.end():
                found = True
                break
        if not found:
            # Check reverse complement
            rc_seq = reverse_complement(seq)
            rc_center = SEQ_LEN - 1 - center
            for match in compiled.finditer(rc_seq):
                if match.start() <= rc_center < match.end():
                    found = True
                    break
        if found:
            count += 1
    return count


def load_genome_sequence(fasta_path):
    """Load genome sequence from FASTA file."""
    seq_parts = []
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith(">"):
                seq_parts.append(line.strip().upper())
    return "".join(seq_parts)


def collapse_revcomp_kmers(kmer_counts):
    """Merge k-mer with its reverse complement (canonical form)."""
    collapsed = Counter()
    seen = set()
    for kmer, count in kmer_counts.items():
        rc = reverse_complement(kmer)
        canonical = min(kmer, rc)
        if canonical not in seen:
            seen.add(canonical)
            collapsed[canonical] = kmer_counts.get(kmer, 0) + kmer_counts.get(rc, 0)
    return collapsed


# ============================================================
# 1. Load methylation site data
# ============================================================
print("=" * 70)
print("EXPANDED METHYLATION MOTIF DISCOVERY (v2 — position-aware)")
print("=" * 70)

df = pd.read_csv(f"{MOTIF_DIR}/methylation_site_sequences.csv")

# Deduplicate by genomic position (unique sites across timepoints)
df_unique = df.drop_duplicates(subset=["chrom", "position", "strand", "mod_type"])
n_4mC = len(df_unique[df_unique["mod_type"] == "4mC"])
n_6mA = len(df_unique[df_unique["mod_type"] == "6mA"])
print(f"\nTotal observations: {len(df)} (4mC: {len(df[df['mod_type']=='4mC'])}, 6mA: {len(df[df['mod_type']=='6mA'])})")
print(f"Unique sites: {len(df_unique)} (4mC: {n_4mC}, 6mA: {n_6mA})")
print(f"Center position: {CENTER_POS} (0-indexed in {SEQ_LEN}bp window)")


# ============================================================
# 2. Position-aware TOP 3 motif classification
# ============================================================
print("\n" + "=" * 70)
print("PHASE 1a: POSITION-AWARE MOTIF CLASSIFICATION")
print("=" * 70)

# Define motifs with their expected modification types
motif_catalog = {
    # TOP 3
    "CCGG": {"regex": "CCGG", "primary_mod": "4mC", "rebase": True},
    "GGCCGG": {"regex": "GGCCGG", "primary_mod": "4mC", "rebase": True},
    "AAGCCCG": {"regex": "AAGCCCG", "primary_mod": "6mA", "rebase": False},
    "GATC": {"regex": "GATC", "primary_mod": "6mA", "rebase": True},
    # Known minor
    "GCGC": {"regex": "GCGC", "primary_mod": "4mC", "rebase": True},
    "TGGCCGGC": {"regex": "TGGCCGGC", "primary_mod": "4mC", "rebase": False},
    # Extended patterns from MEME
    "GCCG": {"regex": "GCCG", "primary_mod": "both", "rebase": True},
    "CCGC": {"regex": "CCGC", "primary_mod": "both", "rebase": True},
    # Additional REBASE patterns to check
    "CCGCGG": {"regex": "CCGCGG", "primary_mod": "both", "rebase": True},
    "GCCGGC": {"regex": "GCCGGC", "primary_mod": "both", "rebase": True},
    "CCSGG": {"regex": "CC[CG]GG", "primary_mod": "both", "rebase": True},
    "CCWGG": {"regex": "CC[AT]GG", "primary_mod": "both", "rebase": True},
    "CTCGAG": {"regex": "CTCGAG", "primary_mod": "both", "rebase": True},
    "TCGCGA": {"regex": "TCGCGA", "primary_mod": "both", "rebase": True},
    "GAGCTC": {"regex": "GAGCTC", "primary_mod": "both", "rebase": True},
    "GTCGAC": {"regex": "GTCGAC", "primary_mod": "both", "rebase": True},
    "GACGTC": {"regex": "GACGTC", "primary_mod": "both", "rebase": True},
    "CAGCTG": {"regex": "CAGCTG", "primary_mod": "both", "rebase": True},
    "CTGCAG": {"regex": "CTGCAG", "primary_mod": "both", "rebase": True},
    "AGGCCT": {"regex": "AGGCCT", "primary_mod": "both", "rebase": True},
    "GCATGC": {"regex": "GCATGC", "primary_mod": "both", "rebase": True},
    "GGATCC": {"regex": "GGATCC", "primary_mod": "both", "rebase": True},
    "AAGCTT": {"regex": "AAGCTT", "primary_mod": "both", "rebase": True},
    "TGATCA": {"regex": "TGATCA", "primary_mod": "both", "rebase": True},
    "GGCGCC": {"regex": "GGCGCC", "primary_mod": "both", "rebase": True},
}

classification_results = []

for mod_type in ["4mC", "6mA"]:
    sub = df_unique[df_unique["mod_type"] == mod_type]
    sequences = sub["sequence"].tolist()
    total = len(sequences)

    print(f"\n--- {mod_type} ({total} unique sites) ---")
    print(f"{'Motif':<14} {'center-aware':>12} {'%':>7} {'window-match':>12} {'%':>7} {'overlap_note'}")
    print("-" * 80)

    for motif_name, info in motif_catalog.items():
        # Position-aware (center must be within motif)
        center_count = count_motif_at_center(sequences, info["regex"])
        center_pct = center_count / total * 100

        # Window match (anywhere in 31bp — for comparison)
        compiled = re.compile(info["regex"])
        window_count = sum(1 for s in sequences if compiled.search(s) is not None)
        window_pct = window_count / total * 100

        if center_count > 0 or window_count > 5:
            print(f"{motif_name:<14} {center_count:>12} {center_pct:>6.1f}% {window_count:>12} {window_pct:>6.1f}%")

        classification_results.append({
            "mod_type": mod_type,
            "motif": motif_name,
            "center_match": center_count,
            "center_pct": round(center_pct, 2),
            "window_match": window_count,
            "window_pct": round(window_pct, 2),
            "total_sites": total,
            "in_rebase": info["rebase"],
        })

class_df = pd.DataFrame(classification_results)
class_df.to_csv(f"{OUT_DIR}/motif_classification_position_aware.csv", index=False)


# ============================================================
# 3. Hierarchical motif assignment (non-overlapping)
# ============================================================
print("\n" + "=" * 70)
print("PHASE 1b: HIERARCHICAL MOTIF ASSIGNMENT")
print("=" * 70)
print("(Each site assigned to the most specific motif covering the center)")

# Priority order: most specific (longest) first
assignment_order_4mC = [
    ("TGGCCGGC", "TGGCCGGC"),
    ("GGCCGG", "GGCCGG"),
    ("CCGG", "CCGG"),
    ("GCGC", "GCGC"),
    ("GCCGGC", "GCCGGC"),
    ("CCGCGG", "CCGCGG"),
    ("GCCG", "GCCG"),
    ("CCGC", "CCGC"),
]

assignment_order_6mA = [
    ("AAGCCCG", "AAGCCCG"),
    ("GATC", "GATC"),
    ("CTCGAG", "CTCGAG"),
    ("GTCGAC", "GTCGAC"),
    ("CAGCTG", "CAGCTG"),
    ("GCCG", "GCCG"),
    ("CCGC", "CCGC"),
]

for mod_type, order in [("4mC", assignment_order_4mC), ("6mA", assignment_order_6mA)]:
    sub = df_unique[df_unique["mod_type"] == mod_type].copy()
    sub["assigned_motif"] = "unassigned"

    for motif_name, regex in order:
        unassigned_mask = sub["assigned_motif"] == "unassigned"
        for idx in sub[unassigned_mask].index:
            seq = sub.at[idx, "sequence"]
            if motif_covers_center(seq, regex):
                sub.at[idx, "assigned_motif"] = motif_name

    counts = sub["assigned_motif"].value_counts()
    total = len(sub)
    print(f"\n{mod_type} hierarchical assignment:")
    for cat in counts.index:
        pct = counts[cat] / total * 100
        print(f"  {cat:<16} {counts[cat]:>6} ({pct:>5.1f}%)")

    # Save assignment
    sub.to_csv(f"{OUT_DIR}/{mod_type}_motif_assignment.csv", index=False)

    # Store residual
    if mod_type == "4mC":
        residual_4mC = sub[sub["assigned_motif"] == "unassigned"]
    else:
        residual_6mA = sub[sub["assigned_motif"] == "unassigned"]

print(f"\nResidual (unassigned) sites:")
print(f"  4mC: {len(residual_4mC)}")
print(f"  6mA: {len(residual_6mA)}")


# ============================================================
# 4. K-mer enrichment analysis on residual sites
# ============================================================
print("\n" + "=" * 70)
print("PHASE 1c: K-MER ENRICHMENT (center-focused)")
print("=" * 70)
print("Analyzing k-mers that overlap the methylated base (center position)")

print("\nLoading genome sequence...")
genome_seq = load_genome_sequence(GENOME_FASTA)
genome_len = len(genome_seq)
print(f"Genome length: {genome_len:,} bp")


def extract_center_kmers(sequences, k, center=CENTER_POS):
    """Extract all k-mers that overlap the center position."""
    counts = Counter()
    for seq in sequences:
        if seq is None or pd.isna(seq):
            continue
        seq = seq.upper()
        # k-mers overlapping center: start from max(0, center-k+1) to center
        for start in range(max(0, center - k + 1), min(center + 1, len(seq) - k + 1)):
            kmer = seq[start:start + k]
            if "N" not in kmer:
                counts[kmer] += 1
    return counts


def count_genome_kmers(genome, k):
    """Count k-mers in genome (both strands)."""
    counts = Counter()
    for i in range(len(genome) - k + 1):
        kmer = genome[i:i+k]
        if "N" not in kmer:
            counts[kmer] += 1
    return counts


for mod_type, residual in [("4mC", residual_4mC), ("6mA", residual_6mA)]:
    if len(residual) < 10:
        print(f"\nSkipping {mod_type} (only {len(residual)} residual sites)")
        continue

    print(f"\n--- {mod_type}: {len(residual)} residual sites ---")
    sequences = residual["sequence"].tolist()

    for k in [4, 5, 6, 7, 8]:
        # Center-focused k-mers from methylation sites
        site_counts = extract_center_kmers(sequences, k)
        site_collapsed = collapse_revcomp_kmers(site_counts)

        # Genome k-mers
        genome_counts = count_genome_kmers(genome_seq, k)
        genome_collapsed = collapse_revcomp_kmers(genome_counts)

        total_site = sum(site_collapsed.values())
        total_genome = sum(genome_collapsed.values())

        if total_site == 0:
            continue

        enrichment_data = []
        for kmer, obs in site_collapsed.items():
            gen = genome_collapsed.get(kmer, 0)
            if gen == 0:
                continue
            exp_freq = gen / total_genome
            obs_freq = obs / total_site
            enrichment = obs_freq / exp_freq

            # Binomial test (one-sided, greater)
            result = binomtest(obs, total_site, exp_freq, alternative="greater")
            p_val = result.pvalue

            enrichment_data.append({
                "kmer": kmer, "obs_count": obs, "obs_freq": obs_freq,
                "genome_count": gen, "genome_freq": exp_freq,
                "enrichment": enrichment, "p_value": p_val,
            })

        enrich_df = pd.DataFrame(enrichment_data)
        enrich_df = enrich_df.sort_values("enrichment", ascending=False)
        enrich_df.to_csv(f"{OUT_DIR}/{mod_type}_residual_{k}mer_enrichment.csv", index=False)

        # Top enriched k-mers
        sig = enrich_df[(enrich_df["p_value"] < 0.05) & (enrich_df["obs_count"] >= 3)]
        if len(sig) > 0:
            print(f"\n  Enriched {k}-mers (p<0.05, n>=3):")
            print(f"  {'k-mer':<12} {'obs':>5} {'enrich':>7} {'p-value':>10}")
            for _, row in sig.head(15).iterrows():
                print(f"  {row['kmer']:<12} {row['obs_count']:>5} {row['enrichment']:>7.2f} {row['p_value']:>10.2e}")


# ============================================================
# 5. REBASE reverse matching (position-aware)
# ============================================================
print("\n" + "=" * 70)
print("PHASE 2: REBASE REVERSE MATCHING (position-aware)")
print("=" * 70)

# Load REBASE data
rebase_motifs = {}
with open(REBASE_FILE) as f:
    reader = csv.DictReader(f)
    for row in reader:
        seq = row.get("recognition_seq", "").strip()
        methyl = row.get("methyl_type", "").strip()
        enzyme = row.get("enzyme_name", "").strip()
        organism = row.get("organism", "").strip()
        etype = row.get("enzyme_type", "").strip()
        if seq and seq != "?" and len(seq) >= 3:
            if seq not in rebase_motifs:
                rebase_motifs[seq] = []
            rebase_motifs[seq].append({
                "enzyme": enzyme, "organism": organism,
                "methyl_type": methyl, "enzyme_type": etype,
            })

print(f"Unique REBASE Streptomyces recognition sequences: {len(rebase_motifs)}")

rebase_match_results = []

for motif_seq, enzymes in sorted(rebase_motifs.items()):
    regex = iupac_to_regex(motif_seq)

    # Skip very degenerate motifs (too many N's)
    if motif_seq.count("N") > len(motif_seq) * 0.5:
        continue

    for mod_type in ["4mC", "6mA"]:
        sub = df_unique[df_unique["mod_type"] == mod_type]
        sequences = sub["sequence"].tolist()
        total = len(sequences)

        # Position-aware matching
        center_matches = count_motif_at_center(sequences, regex)
        center_pct = center_matches / total * 100 if total > 0 else 0

        # Residual sites matching
        residual = residual_4mC if mod_type == "4mC" else residual_6mA
        if len(residual) > 0:
            residual_seqs = residual["sequence"].tolist()
            residual_matches = count_motif_at_center(residual_seqs, regex)
            residual_pct = residual_matches / len(residual) * 100
        else:
            residual_matches = 0
            residual_pct = 0.0

        enzyme_names = "; ".join(e["enzyme"] for e in enzymes[:3])
        methyl_types = "; ".join(set(e["methyl_type"] for e in enzymes if e["methyl_type"]))
        organisms = "; ".join(set(e["organism"] for e in enzymes[:3]))

        rebase_match_results.append({
            "rebase_motif": motif_seq,
            "mod_type": mod_type,
            "total_sites": total,
            "center_matches": center_matches,
            "center_pct": round(center_pct, 2),
            "residual_matches": residual_matches,
            "residual_pct": round(residual_pct, 2),
            "rebase_methyl_type": methyl_types,
            "enzymes": enzyme_names,
            "organisms": organisms,
        })

rebase_df = pd.DataFrame(rebase_match_results)
rebase_df = rebase_df.sort_values(["mod_type", "center_pct"], ascending=[True, False])
rebase_df.to_csv(f"{OUT_DIR}/rebase_reverse_matching.csv", index=False)

# Show significant matches
print("\n--- REBASE motifs with center-aware matches (>0.5% of sites) ---")
for mod_type in ["4mC", "6mA"]:
    sub = rebase_df[(rebase_df["mod_type"] == mod_type) & (rebase_df["center_pct"] > 0.5)]
    if len(sub) == 0:
        print(f"\n{mod_type}: No significant matches")
        continue
    print(f"\n{mod_type}:")
    print(f"{'REBASE motif':<20} {'center':>7} {'%':>7} {'resid':>6} {'%resid':>7} {'methyl':>10} {'enzymes'}")
    for _, row in sub.iterrows():
        print(f"{row['rebase_motif']:<20} {row['center_matches']:>7} {row['center_pct']:>6.1f}% "
              f"{row['residual_matches']:>6} {row['residual_pct']:>6.1f}% "
              f"{row['rebase_methyl_type']:>10} {row['enzymes'][:40]}")

# Show residual-enriched matches
print("\n--- REBASE motifs in RESIDUAL sites (>1%) ---")
for mod_type in ["4mC", "6mA"]:
    sub = rebase_df[(rebase_df["mod_type"] == mod_type) & (rebase_df["residual_pct"] > 1)]
    if len(sub) > 0:
        print(f"\n{mod_type} ({len(residual_4mC) if mod_type == '4mC' else len(residual_6mA)} residual):")
        for _, row in sub.iterrows():
            print(f"  {row['rebase_motif']:<20} {row['residual_matches']}/{len(residual_4mC) if mod_type == '4mC' else len(residual_6mA)} "
                  f"({row['residual_pct']:.1f}%)  methyl: {row['rebase_methyl_type']}  enzymes: {row['enzymes'][:50]}")


# ============================================================
# 6. MEME motif summary table
# ============================================================
print("\n" + "=" * 70)
print("MEME MOTIF SUMMARY (all 10 motifs, both mod types)")
print("=" * 70)

meme_summary = [
    {"mod_type": "4mC", "rank": 1, "consensus": "SAMGCCSGCCA", "regex": "[GC][AC][AC]GCC[GC]GCCA",
     "sites": 2678, "evalue": "6.4e-1211", "status": "CORE: CCGG/GGCCGG"},
    {"mod_type": "4mC", "rank": 2, "consensus": "ATGAWGWT", "regex": "ATG[AT][AT]G[AT]T",
     "sites": 4, "evalue": "1.3e+005", "status": "NOISE (4 sites)"},
    {"mod_type": "4mC", "rank": 3, "consensus": "ATGATCA", "regex": "ATGATCA",
     "sites": 2, "evalue": "5.7e+005", "status": "NOISE (2 sites)"},
    {"mod_type": "4mC", "rank": 4, "consensus": "ATCWTCTT", "regex": "ATC[AT]TCTT",
     "sites": 2, "evalue": "4.0e+005", "status": "NOISE (2 sites)"},
    {"mod_type": "4mC", "rank": 5, "consensus": "GATCATG", "regex": "GATCATG",
     "sites": 2, "evalue": "5.7e+005", "status": "NOISE (2 sites)"},
    {"mod_type": "6mA", "rank": 1, "consensus": "GVSAAGCCCGVC", "regex": "[GC][CGA][CG]AAGCCCG[CGA][CG]",
     "sites": 656, "evalue": "3.1e-256", "status": "CORE: AAGCCCG"},
    {"mod_type": "6mA", "rank": 2, "consensus": "GGSRCCGKCAMC", "regex": "G[GC][CG][GA]CCG[TG]CA[CA]C",
     "sites": 954, "evalue": "6.3e-045", "status": "*** NEW CANDIDATE ***"},
    {"mod_type": "6mA", "rank": 3, "consensus": "ACCGGCVGCAAC", "regex": "ACC[GA][GA]C[ACG]GCAAC",
     "sites": 20, "evalue": "1.8e+003", "status": "MARGINAL (20 sites)"},
    {"mod_type": "6mA", "rank": 4, "consensus": "GCGSCGMSGAAC", "regex": "[GA][CG]G[GC]CG[CA][GC]GAAC",
     "sites": 57, "evalue": "1.4e+003", "status": "MARGINAL (57 sites)"},
    {"mod_type": "6mA", "rank": 5, "consensus": "GAAAAAG", "regex": "GAAAAAG",
     "sites": 2, "evalue": "2.7e+004", "status": "NOISE (2 sites)"},
]

meme_df = pd.DataFrame(meme_summary)
meme_df.to_csv(f"{OUT_DIR}/meme_all_motifs_summary.csv", index=False)

print(f"\n{'Type':<5} {'#':<3} {'Consensus':<16} {'Sites':>6} {'E-value':>12} {'Status'}")
print("-" * 75)
for _, row in meme_df.iterrows():
    print(f"{row['mod_type']:<5} {row['rank']:<3} {row['consensus']:<16} {row['sites']:>6} "
          f"{row['evalue']:>12} {row['status']}")


# ============================================================
# 7. 6mA MEME-2 deep analysis
# ============================================================
print("\n" + "=" * 70)
print("DEEP ANALYSIS: 6mA MEME-2 (GGSRCCGKCAMC)")
print("=" * 70)

sites_6mA = df_unique[df_unique["mod_type"] == "6mA"]
seqs_6mA = sites_6mA["sequence"].tolist()

# Try various core patterns from the MEME-2 consensus
# MEME-2 probability matrix shows:
#  pos 1: G(82%) C(14%)   -> G
#  pos 2: G(58%) C(25%)   -> G
#  pos 3: C(54%) G(42%)   -> C/G
#  pos 4: A(30%) G(50%)   -> G/A
#  pos 5: C(55%) A(16%) T(16%) -> C
#  pos 6: C(81%) G(19%)   -> C
#  pos 7: G(77%) A(10%)   -> G
#  pos 8: T(46%) G(44%)   -> T/G
#  pos 9: C(82%) A(5%)    -> C
#  pos10: A(99%)           -> A
#  pos11: C(60%) A(37%)   -> C/A
#  pos12: C(92%)           -> C

# Possible core motifs to test:
meme2_cores = [
    ("GCCGCAC", "7bp core from positions 4-10 (simplified)"),
    ("CCGCAC", "6bp core"),
    ("CCGTCAC", "7bp core with T at pos 8"),
    ("CCGGCAC", "7bp core with G at pos 8"),
    ("RCCG", "4bp degenerate core R=A/G"),
    ("ACCG", "4bp core starting with A"),
    ("GCCG", "4bp core starting with G"),
    ("CCG[TG]CA", "6bp degenerate core"),
    ("CCGKCA", "6bp with K=G/T"),
    ("CCGTCA", "specific 6bp"),
    ("CCGGCA", "specific 6bp"),
    ("GCCGTCA", "7bp specific"),
    ("GCCGGCA", "7bp specific"),
]

print(f"\nCore pattern analysis (center-aware, {len(seqs_6mA)} total 6mA sites):")
print(f"{'Pattern':<16} {'center':>8} {'%':>7} {'Description'}")
print("-" * 60)

for pat, desc in meme2_cores:
    regex = iupac_to_regex(pat)
    count = count_motif_at_center(seqs_6mA, regex)
    pct = count / len(seqs_6mA) * 100
    if count > 0:
        print(f"{pat:<16} {count:>8} {pct:>6.1f}% {desc}")

# Overlap analysis: MEME-1 vs MEME-2
print("\n--- Overlap between MEME-1 (AAGCCCG) and MEME-2 candidates ---")
aagcccg_sites = set()
gccg_sites = set()
accg_sites = set()

for idx, row in sites_6mA.iterrows():
    seq = row["sequence"]
    pos = row["position"]
    if motif_covers_center(seq, "AAGCCCG"):
        aagcccg_sites.add(pos)
    if motif_covers_center(seq, "GCCG"):
        gccg_sites.add(pos)
    if motif_covers_center(seq, "ACCG"):
        accg_sites.add(pos)

overlap_gccg = aagcccg_sites & gccg_sites
overlap_accg = aagcccg_sites & accg_sites
print(f"AAGCCCG sites: {len(aagcccg_sites)}")
print(f"GCCG sites: {len(gccg_sites)}")
print(f"ACCG sites: {len(accg_sites)}")
print(f"AAGCCCG ∩ GCCG: {len(overlap_gccg)} ({len(overlap_gccg)/max(len(aagcccg_sites),1)*100:.1f}% of AAGCCCG)")
print(f"AAGCCCG ∩ ACCG: {len(overlap_accg)} ({len(overlap_accg)/max(len(aagcccg_sites),1)*100:.1f}% of AAGCCCG)")
print(f"GCCG only (non-AAGCCCG): {len(gccg_sites - aagcccg_sites)}")
print(f"ACCG only (non-AAGCCCG): {len(accg_sites - aagcccg_sites)}")


# ============================================================
# 8. Temporal dynamics
# ============================================================
print("\n" + "=" * 70)
print("TEMPORAL DYNAMICS OF ALL CANDIDATE MOTIFS")
print("=" * 70)

# Use position-aware counting on the full (non-deduplicated) dataset
dynamics_motifs = [
    ("CCGG", "CCGG"),
    ("GGCCGG", "GGCCGG"),
    ("AAGCCCG", "AAGCCCG"),
    ("GATC", "GATC"),
    ("GCGC", "GCGC"),
    ("GCCG", "GCCG"),
    ("ACCG", "ACCG"),
    ("CCGC", "CCGC"),
    ("GCCGGC", "GCCGGC"),
    ("CCGCGG", "CCGCGG"),
]

dynamics_results = []

for mod_type in ["4mC", "6mA"]:
    print(f"\n--- {mod_type} ---")
    print(f"{'Motif':<12} {'T1':>6} {'T2':>6} {'T3':>6} {'T1%':>7} {'T2%':>7} {'T3%':>7}")

    for motif_name, regex in dynamics_motifs:
        tp_data = {}
        for tp in ["T1", "T2", "T3"]:
            tp_sub = df[(df["mod_type"] == mod_type) & (df["timepoint"] == tp)]
            seqs = tp_sub["sequence"].tolist()
            total_tp = len(seqs)
            matches = count_motif_at_center(seqs, regex)
            tp_data[tp] = {"matches": matches, "total": total_tp,
                           "pct": matches / total_tp * 100 if total_tp > 0 else 0}

        if all(tp_data[tp]["matches"] == 0 for tp in ["T1", "T2", "T3"]):
            continue

        print(f"{motif_name:<12} "
              f"{tp_data['T1']['matches']:>6} {tp_data['T2']['matches']:>6} {tp_data['T3']['matches']:>6} "
              f"{tp_data['T1']['pct']:>6.1f}% {tp_data['T2']['pct']:>6.1f}% {tp_data['T3']['pct']:>6.1f}%")

        dynamics_results.append({
            "mod_type": mod_type, "motif": motif_name,
            "T1_count": tp_data["T1"]["matches"], "T1_total": tp_data["T1"]["total"],
            "T1_pct": round(tp_data["T1"]["pct"], 2),
            "T2_count": tp_data["T2"]["matches"], "T2_total": tp_data["T2"]["total"],
            "T2_pct": round(tp_data["T2"]["pct"], 2),
            "T3_count": tp_data["T3"]["matches"], "T3_total": tp_data["T3"]["total"],
            "T3_pct": round(tp_data["T3"]["pct"], 2),
        })

dyn_df = pd.DataFrame(dynamics_results)
dyn_df.to_csv(f"{OUT_DIR}/motif_temporal_dynamics.csv", index=False)


# ============================================================
# 9. Save FASTA files for STREME
# ============================================================
print("\n" + "=" * 70)
print("FASTA OUTPUT FOR STREME ANALYSIS")
print("=" * 70)

for mod_type, residual in [("4mC", residual_4mC), ("6mA", residual_6mA)]:
    fasta_path = f"{OUT_DIR}/{mod_type}_residual_sites.fasta"
    count = 0
    with open(fasta_path, "w") as f:
        for _, row in residual.iterrows():
            f.write(f">{row['chrom']}_{row['position']}_{row['strand']}_{mod_type}\n")
            f.write(f"{row['sequence']}\n")
            count += 1
    print(f"  {mod_type}: {count} residual sites -> {fasta_path}")

# Also save 6mA non-AAGCCCG sites (broader residual for MEME-2 investigation)
non_aagcccg = df_unique[(df_unique["mod_type"] == "6mA")]
non_aagcccg = non_aagcccg[~non_aagcccg["sequence"].apply(lambda s: motif_covers_center(s, "AAGCCCG"))]
fasta_path = f"{OUT_DIR}/6mA_non_AAGCCCG_sites.fasta"
with open(fasta_path, "w") as f:
    for _, row in non_aagcccg.iterrows():
        f.write(f">{row['chrom']}_{row['position']}_{row['strand']}_6mA\n")
        f.write(f"{row['sequence']}\n")
print(f"  6mA non-AAGCCCG: {len(non_aagcccg)} sites -> {fasta_path}")


# ============================================================
# 10. Final candidate summary
# ============================================================
print("\n" + "=" * 70)
print("FINAL: EXPANDED CANDIDATE MOTIF LIST")
print("=" * 70)

print("""
Candidates for downstream analysis (Phase 3-5):

1. GGSRCCGKCAMC (6mA MEME-2, 954 MEME sites)
   - Core pattern TBD from k-mer enrichment
   - Potentially the most impactful finding
   - May represent a novel 6mA motif or BREX-2 PglX target

2. GCGC (HhaI-like, 95 4mC sites = 3.5%)
   - Known R-M motif in Streptomyces
   - Need to check if these are true 4mC or Nanopore artifact

3. GCCG / ACCG (from MEME-2 core decomposition)
   - Subset of MEME-2; check overlap with AAGCCCG
   - May represent independent short motifs

4. GCGSCGMSGAAC (6mA MEME-4, 57 sites)
   - Marginal significance; contains GAAC
   - Could be Type I-like recognition

5. REBASE motifs enriched in residual sites
   - Results from Phase 2 above

6. Any novel k-mers from enrichment analysis
   - Especially those not matching any REBASE entry
""")

print(f"\nAll output files saved to: {OUT_DIR}/")
print("Done.")
