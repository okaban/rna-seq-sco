#!/usr/bin/env python3
"""
Per-read co-modification analysis for AAGCCCG motif.
Determines what fraction of AAGCCCG-containing reads show simultaneous
4mC (at C3/C5) and 6mA (at A0/A1) modification on the same read.
"""

import sys
import re
from collections import defaultdict
from pathlib import Path

# ── Configuration ──────────────────────────────────────────────────────────────
REF_FA = "/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa"
EXTRACT_FILES = {
    "T1_1-1": "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/65_per_read_comod/1-1_extract.tsv",
    # Add more as generated:
    # "T1_1-2": "...",
    # "T1_1-3": "...",
}

MOTIF = "AAGCCCG"
# 0-indexed positions within AAGCCCG: A=0,A=1,G=2,C=3,C=4,C=5,G=6
MOD_6mA_POSITIONS = {0, 1}   # A0, A1
MOD_4mC_POSITIONS = {3, 5}   # C3, C5 (C4 not modified per analysis)

CALL_6mA = "a"
CALL_4mC = "21839"
MIN_PROB = 0.5   # minimum call probability to count as modified

def load_reference(fa_path: str) -> str:
    """Load first sequence from FASTA as uppercase string."""
    seq_parts = []
    with open(fa_path) as f:
        for line in f:
            if line.startswith(">"):
                if seq_parts:
                    break  # only need first contig (NC_003888.3)
            else:
                seq_parts.append(line.strip().upper())
    return "".join(seq_parts)


def find_motif_positions(genome: str, motif: str):
    """
    Find all occurrences of motif on + and - strands.
    Returns dict: ref_position -> (strand, pos_in_motif)
    for positions that correspond to 6mA or 4mC sites.
    """
    pos_map = {}   # ref_pos (0-based) -> set of (motif_start, pos_in_motif, strand, mod_type)
    motif_len = len(motif)
    rc_motif = motif[::-1].translate(str.maketrans("ACGT", "TGCA"))

    def register(motif_start, pos_in_motif, strand):
        ref_pos = motif_start + pos_in_motif
        if ref_pos not in pos_map:
            pos_map[ref_pos] = []
        if pos_in_motif in MOD_6mA_POSITIONS:
            pos_map[ref_pos].append((motif_start, pos_in_motif, strand, "6mA"))
        elif pos_in_motif in MOD_4mC_POSITIONS:
            pos_map[ref_pos].append((motif_start, pos_in_motif, strand, "4mC"))

    # + strand
    for m in re.finditer(f"(?={motif})", genome):
        start = m.start()
        for p in list(MOD_6mA_POSITIONS) + list(MOD_4mC_POSITIONS):
            register(start, p, "+")

    # - strand: the reverse complement appears on - strand
    # For - strand occurrences: ref_pos of a base = motif_start + (motif_len - 1 - pos_in_motif)
    for m in re.finditer(f"(?={rc_motif})", genome):
        rc_start = m.start()
        for p in list(MOD_6mA_POSITIONS) + list(MOD_4mC_POSITIONS):
            # actual ref position on - strand occurrence
            ref_pos = rc_start + (motif_len - 1 - p)
            if ref_pos not in pos_map:
                pos_map[ref_pos] = []
            if p in MOD_6mA_POSITIONS:
                pos_map[ref_pos].append((rc_start, p, "-", "6mA"))
            elif p in MOD_4mC_POSITIONS:
                pos_map[ref_pos].append((rc_start, p, "-", "4mC"))

    return pos_map


def analyze_extract(extract_path: str, pos_map: dict, sample_name: str):
    """
    Stream through extract file. For each read, collect 4mC/6mA calls
    at AAGCCCG positions. Report co-modification statistics.
    """
    # Per read: track which motif_starts have 6mA calls, which have 4mC calls
    # key: read_id -> motif_start -> set of mod_types confirmed
    read_motif_mods = defaultdict(lambda: defaultdict(set))

    n_lines = 0
    with open(extract_path) as f:
        header = f.readline()  # skip header
        for line in f:
            n_lines += 1
            if n_lines % 10_000_000 == 0:
                print(f"  [{sample_name}] Processed {n_lines/1e6:.0f}M lines...", flush=True)

            parts = line.split("\t")
            # columns (0-indexed): 0=read_id, 2=ref_position, 12=call_prob, 13=call_code, 19=fail
            call_code = parts[13]
            if call_code not in (CALL_6mA, CALL_4mC):
                continue

            try:
                ref_pos = int(parts[2])
            except ValueError:
                continue

            if ref_pos not in pos_map:
                continue

            try:
                call_prob = float(parts[12])
            except ValueError:
                continue

            if call_prob < MIN_PROB:
                continue

            fail = parts[19].strip()
            if fail.lower() == "true":
                continue

            read_id = parts[0]
            entries = pos_map[ref_pos]

            for (motif_start, pos_in_motif, strand, expected_mod_type) in entries:
                # Verify call type matches expected modification
                if expected_mod_type == "6mA" and call_code == CALL_6mA:
                    read_motif_mods[read_id][motif_start].add("6mA")
                elif expected_mod_type == "4mC" and call_code == CALL_4mC:
                    read_motif_mods[read_id][motif_start].add("4mC")

    # Aggregate
    total_reads_with_motif_any = len(read_motif_mods)
    reads_6mA_only = 0
    reads_4mC_only = 0
    reads_both = 0
    reads_neither = 0  # reads touching motif positions but no modification called above threshold

    motif_comod_counts = defaultdict(lambda: {"6mA": 0, "4mC": 0, "both": 0})

    for read_id, motif_dict in read_motif_mods.items():
        for motif_start, mod_set in motif_dict.items():
            has_6mA = "6mA" in mod_set
            has_4mC = "4mC" in mod_set
            motif_comod_counts[motif_start]["6mA"] += has_6mA
            motif_comod_counts[motif_start]["4mC"] += has_4mC
            if has_6mA and has_4mC:
                motif_comod_counts[motif_start]["both"] += 1

    # Summarize across all motif instances
    total_motif_read_pairs = sum(
        sum(v.values()) // 3 * 3  # rough count
        for v in motif_comod_counts.values()
    )

    # Better: count per (read, motif) pairs
    pairs_with_6mA = 0
    pairs_with_4mC = 0
    pairs_with_both = 0
    pairs_total = 0

    for read_id, motif_dict in read_motif_mods.items():
        for motif_start, mod_set in motif_dict.items():
            pairs_total += 1
            has_6mA = "6mA" in mod_set
            has_4mC = "4mC" in mod_set
            if has_6mA:
                pairs_with_6mA += 1
            if has_4mC:
                pairs_with_4mC += 1
            if has_6mA and has_4mC:
                pairs_with_both += 1

    print(f"\n{'='*60}")
    print(f"SAMPLE: {sample_name}")
    print(f"{'='*60}")
    print(f"Total reads with >=1 AAGCCCG position call: {total_reads_with_motif_any:,}")
    print(f"\n(read × motif_occurrence) pairs:")
    print(f"  Total pairs:              {pairs_total:,}")
    print(f"  Pairs with 6mA (A0/A1):  {pairs_with_6mA:,} ({100*pairs_with_6mA/max(pairs_total,1):.1f}%)")
    print(f"  Pairs with 4mC (C3/C5):  {pairs_with_4mC:,} ({100*pairs_with_4mC/max(pairs_total,1):.1f}%)")
    print(f"  Pairs with BOTH:          {pairs_with_both:,} ({100*pairs_with_both/max(pairs_total,1):.1f}%)")
    if pairs_with_6mA > 0:
        print(f"\n  Of pairs with 6mA: {100*pairs_with_both/pairs_with_6mA:.1f}% also have 4mC")
    if pairs_with_4mC > 0:
        print(f"  Of pairs with 4mC: {100*pairs_with_both/pairs_with_4mC:.1f}% also have 6mA")

    # Top co-modified motif instances
    top_comods = sorted(motif_comod_counts.items(), key=lambda x: x[1]["both"], reverse=True)[:10]
    print(f"\nTop AAGCCCG instances by co-modification count (motif_start: 6mA / 4mC / both):")
    for ms, counts in top_comods:
        if counts["both"] > 0:
            print(f"  pos {ms}: 6mA={counts['6mA']}, 4mC={counts['4mC']}, both={counts['both']}")

    return {
        "sample": sample_name,
        "total_reads": total_reads_with_motif_any,
        "pairs_total": pairs_total,
        "pairs_6mA": pairs_with_6mA,
        "pairs_4mC": pairs_with_4mC,
        "pairs_both": pairs_with_both,
    }


def main():
    print("Loading reference genome...", flush=True)
    genome = load_reference(REF_FA)
    print(f"  Genome length: {len(genome):,} bp", flush=True)

    print(f"Finding {MOTIF} motif positions...", flush=True)
    pos_map = find_motif_positions(genome, MOTIF)
    # Count unique motif instances
    motif_starts_plus = set()
    motif_starts_minus = set()
    for ref_pos, entries in pos_map.items():
        for (ms, pip, strand, mod) in entries:
            if strand == "+":
                motif_starts_plus.add(ms)
            else:
                motif_starts_minus.add(ms)
    print(f"  AAGCCCG occurrences: + strand={len(motif_starts_plus)}, - strand={len(motif_starts_minus)}", flush=True)
    print(f"  Total ref positions to track: {len(pos_map):,}", flush=True)

    results = []
    for sample_name, extract_path in EXTRACT_FILES.items():
        if not Path(extract_path).exists():
            print(f"SKIP {sample_name}: {extract_path} not found")
            continue
        print(f"\nAnalyzing {sample_name} ({extract_path})...", flush=True)
        r = analyze_extract(extract_path, pos_map, sample_name)
        results.append(r)

    print("\n\n=== SUMMARY TABLE ===")
    print(f"{'Sample':<15} {'Reads':<10} {'Pairs':<10} {'6mA%':<8} {'4mC%':<8} {'Both%':<8}")
    for r in results:
        n = r["pairs_total"]
        print(f"{r['sample']:<15} {r['total_reads']:<10,} {n:<10,} "
              f"{100*r['pairs_6mA']/max(n,1):<8.1f} "
              f"{100*r['pairs_4mC']/max(n,1):<8.1f} "
              f"{100*r['pairs_both']/max(n,1):<8.1f}")


if __name__ == "__main__":
    main()
