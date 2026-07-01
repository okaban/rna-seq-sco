#!/usr/bin/env python3
"""Same as fast_comod.py but ignores the fail flag - uses only prob threshold."""
import sys, re
from collections import defaultdict

REF_FA = "/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa"
MOTIF = "AAGCCCG"
MOD_6mA_POSITIONS = {0, 1}
MOD_4mC_POSITIONS = {3, 5}
CALL_6mA = "a"
CALL_4mC = "21839"
MIN_PROB = 0.5

def load_reference(fa_path):
    parts = []
    with open(fa_path) as f:
        for line in f:
            if line.startswith(">"):
                if parts: break
            else:
                parts.append(line.strip().upper())
    return "".join(parts)

def find_motif_positions(genome, motif):
    pos_map = {}
    motif_len = len(motif)
    rc_motif = motif[::-1].translate(str.maketrans("ACGT", "TGCA"))
    target = MOD_6mA_POSITIONS | MOD_4mC_POSITIONS

    def reg(ms, pip, strand):
        rp = ms + pip if strand == "+" else ms + (motif_len - 1 - pip)
        pos_map.setdefault(rp, [])
        mt = "6mA" if pip in MOD_6mA_POSITIONS else "4mC"
        pos_map[rp].append((ms, pip, strand, mt))

    for m in re.finditer(f"(?={motif})", genome):
        for p in target: reg(m.start(), p, "+")
    for m in re.finditer(f"(?={rc_motif})", genome):
        for p in target: reg(m.start(), p, "-")
    return pos_map

def analyze(modonly_path, pos_map, sample_name, min_prob=MIN_PROB):
    read_motif = defaultdict(lambda: defaultdict(set))
    n = 0
    with open(modonly_path) as f:
        f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4: continue
            read_id, ref_pos_str, prob_str, call_code = parts[:4]
            if call_code not in (CALL_6mA, CALL_4mC): continue
            try:
                ref_pos = int(ref_pos_str)
                prob = float(prob_str)
            except ValueError:
                continue
            if prob < min_prob: continue
            if ref_pos not in pos_map: continue
            for (ms, pip, strand, em) in pos_map[ref_pos]:
                if em == "6mA" and call_code == CALL_6mA:
                    read_motif[read_id][ms].add("6mA")
                elif em == "4mC" and call_code == CALL_4mC:
                    read_motif[read_id][ms].add("4mC")
            n += 1

    pairs_total = pairs_6mA = pairs_4mC = pairs_both = 0
    for rid, md in read_motif.items():
        for ms, mods in md.items():
            h6 = "6mA" in mods; h4 = "4mC" in mods
            pairs_total += 1
            pairs_6mA += h6; pairs_4mC += h4; pairs_both += (h6 and h4)

    pct = lambda x: f"{100*x/max(pairs_total,1):.1f}%"
    print(f"\n=== {sample_name} (fail-ignored, prob>={min_prob}) ===")
    print(f"Reads with AAGCCCG calls: {len(read_motif):,}")
    print(f"Total (read×motif) pairs: {pairs_total:,}")
    print(f"  6mA: {pairs_6mA:,} ({pct(pairs_6mA)})")
    print(f"  4mC: {pairs_4mC:,} ({pct(pairs_4mC)})")
    print(f"  BOTH: {pairs_both:,} ({pct(pairs_both)})")
    if pairs_6mA: print(f"  Of 6mA pairs: {100*pairs_both/pairs_6mA:.1f}% have 4mC")
    if pairs_4mC: print(f"  Of 4mC pairs: {100*pairs_both/pairs_4mC:.1f}% have 6mA")
    return pairs_total, pairs_6mA, pairs_4mC, pairs_both

if __name__ == "__main__":
    modonly = sys.argv[1]
    sample = sys.argv[2] if len(sys.argv) > 2 else "unknown"
    genome = load_reference(REF_FA)
    pos_map = find_motif_positions(genome, MOTIF)
    analyze(modonly, pos_map, sample)
