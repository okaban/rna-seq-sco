#!/usr/bin/env python3
"""
Fast per-read co-modification analysis using pre-filtered modonly TSV.
Input: read_id, ref_position, call_prob, call_code, fail
(from awk-filtered extract calls output)
"""
import sys
import re
from collections import defaultdict

REF_FA = "/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa"
MOTIF = "AAGCCCG"
MOD_6mA_POSITIONS = {0, 1}   # A0, A1 in motif
MOD_4mC_POSITIONS = {3, 5}   # C3, C5 in motif
CALL_6mA = "a"
CALL_4mC = "21839"
MIN_PROB = 0.5

def load_reference(fa_path):
    parts = []
    with open(fa_path) as f:
        for line in f:
            if line.startswith(">"):
                if parts:
                    break
            else:
                parts.append(line.strip().upper())
    return "".join(parts)

def find_motif_positions(genome, motif):
    pos_map = {}
    motif_len = len(motif)
    rc_motif = motif[::-1].translate(str.maketrans("ACGT", "TGCA"))
    target_positions = MOD_6mA_POSITIONS | MOD_4mC_POSITIONS

    def register(motif_start, pos_in_motif, strand):
        ref_pos = motif_start + pos_in_motif if strand == "+" else motif_start + (motif_len - 1 - pos_in_motif)
        if ref_pos not in pos_map:
            pos_map[ref_pos] = []
        mod_type = "6mA" if pos_in_motif in MOD_6mA_POSITIONS else "4mC"
        pos_map[ref_pos].append((motif_start, pos_in_motif, strand, mod_type))

    for m in re.finditer(f"(?={motif})", genome):
        start = m.start()
        for p in target_positions:
            register(start, p, "+")

    for m in re.finditer(f"(?={rc_motif})", genome):
        start = m.start()
        for p in target_positions:
            register(start, p, "-")

    return pos_map

def analyze(modonly_path, pos_map, sample_name):
    # read_id -> motif_start -> set of mod_types
    read_motif = defaultdict(lambda: defaultdict(set))
    n = 0

    with open(modonly_path) as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            read_id, ref_pos_str, prob_str, call_code, fail = parts
            if call_code not in (CALL_6mA, CALL_4mC):
                continue
            try:
                ref_pos = int(ref_pos_str)
                prob = float(prob_str)
            except ValueError:
                continue
            if prob < MIN_PROB:
                continue
            if fail.strip().lower() == "true":
                continue
            if ref_pos not in pos_map:
                continue
            for (motif_start, pos_in_motif, strand, expected_mod) in pos_map[ref_pos]:
                if expected_mod == "6mA" and call_code == CALL_6mA:
                    read_motif[read_id][motif_start].add("6mA")
                elif expected_mod == "4mC" and call_code == CALL_4mC:
                    read_motif[read_id][motif_start].add("4mC")
            n += 1

    # Aggregate (read × motif) pairs
    pairs_total = pairs_6mA = pairs_4mC = pairs_both = 0
    motif_stats = defaultdict(lambda: {"6mA": 0, "4mC": 0, "both": 0})

    for rid, mdict in read_motif.items():
        for ms, mods in mdict.items():
            has6 = "6mA" in mods
            has4 = "4mC" in mods
            pairs_total += 1
            if has6:
                pairs_6mA += 1
                motif_stats[ms]["6mA"] += 1
            if has4:
                pairs_4mC += 1
                motif_stats[ms]["4mC"] += 1
            if has6 and has4:
                pairs_both += 1
                motif_stats[ms]["both"] += 1

    print(f"\n{'='*60}")
    print(f"SAMPLE: {sample_name}")
    print(f"{'='*60}")
    print(f"Total filtered modification calls processed: {n:,}")
    print(f"Unique reads with AAGCCCG position calls: {len(read_motif):,}")
    print(f"\n(read × motif_occurrence) pairs:")
    pct = lambda x: f"{100*x/max(pairs_total,1):.1f}%"
    print(f"  Total pairs:              {pairs_total:,}")
    print(f"  With 6mA (A0/A1):         {pairs_6mA:,} ({pct(pairs_6mA)})")
    print(f"  With 4mC (C3/C5):         {pairs_4mC:,} ({pct(pairs_4mC)})")
    print(f"  With BOTH (co-modified):  {pairs_both:,} ({pct(pairs_both)})")
    if pairs_6mA:
        print(f"\n  Of 6mA-modified pairs: {100*pairs_both/pairs_6mA:.1f}% also have 4mC")
    if pairs_4mC:
        print(f"  Of 4mC-modified pairs: {100*pairs_both/pairs_4mC:.1f}% also have 6mA")

    # Top co-modified loci
    top = sorted(motif_stats.items(), key=lambda x: x[1]["both"], reverse=True)[:5]
    print(f"\nTop AAGCCCG loci by co-modification:")
    for ms, s in top:
        if s["both"] > 0:
            print(f"  ref pos {ms}: reads with 6mA={s['6mA']}, 4mC={s['4mC']}, both={s['both']}")

    return dict(sample=sample_name, pairs_total=pairs_total,
                pairs_6mA=pairs_6mA, pairs_4mC=pairs_4mC, pairs_both=pairs_both,
                n_reads=len(read_motif))

if __name__ == "__main__":
    import sys
    modonly = sys.argv[1] if len(sys.argv) > 1 else None
    sample = sys.argv[2] if len(sys.argv) > 2 else "unknown"
    if not modonly:
        print("Usage: fast_comod.py <modonly.tsv> <sample_name>")
        sys.exit(1)

    print("Loading reference...", flush=True)
    genome = load_reference(REF_FA)
    print(f"Genome: {len(genome):,} bp", flush=True)
    print("Building motif position map...", flush=True)
    pos_map = find_motif_positions(genome, MOTIF)
    plus = len({ms for entries in pos_map.values() for ms, _, s, _ in entries if s == "+"})
    minus = len({ms for entries in pos_map.values() for ms, _, s, _ in entries if s == "-"})
    print(f"AAGCCCG: +strand={plus}, -strand={minus}, tracked_positions={len(pos_map):,}", flush=True)
    print(f"\nAnalyzing {modonly}...", flush=True)
    analyze(modonly, pos_map, sample)
