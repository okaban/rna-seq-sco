#!/usr/bin/env python3
"""
Sensitivity analysis for per-read co-modification at AAGCCCG.
Three analyses:
  1. fail=false filter effect (AAGCCCG, T1 pooled)
  2. call_prob threshold sensitivity: 0.5 / 0.7 / 0.9 (AAGCCCG, T1 pooled)
  3. GCCGGC negative control (T1 pooled)
      - 4mC at GCCGGC internal C positions
      - 6mA at nearest A positions within ±4bp of the motif (flanking or adjacent)
      - Compare co-modification rate to AAGCCCG ~50%
"""
import re
from collections import defaultdict
from scipy.stats import fisher_exact
from statsmodels.stats.proportion import proportions_ztest

REF_FA = "/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa"

CALL_6mA = "a"
CALL_4mC = "21839"

# ── Reference loading ──────────────────────────────────────────────────────────
def load_reference(fa_path: str) -> str:
    parts = []
    with open(fa_path) as f:
        for line in f:
            if line.startswith(">"):
                if parts:
                    break
            else:
                parts.append(line.strip().upper())
    return "".join(parts)

def rc(seq: str) -> str:
    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))

# ── Motif position map builder ─────────────────────────────────────────────────
def build_pos_map(genome: str, motif: str,
                  mod4mC_positions: set, mod6mA_positions: set) -> dict:
    """
    Returns: ref_pos (0-based) -> list of (motif_start, pos_in_motif, strand, mod_type)
    mod_type: '4mC' or '6mA'
    """
    pos_map = {}
    motif_len = len(motif)
    rc_motif = rc(motif)
    target = mod4mC_positions | mod6mA_positions

    def reg(ms, pip, strand):
        rp = ms + pip if strand == "+" else ms + (motif_len - 1 - pip)
        pos_map.setdefault(rp, [])
        mt = "6mA" if pip in mod6mA_positions else "4mC"
        pos_map[rp].append((ms, pip, strand, mt))

    for m in re.finditer(f"(?={motif})", genome):
        for p in target:
            reg(m.start(), p, "+")
    for m in re.finditer(f"(?={rc_motif})", genome):
        for p in target:
            reg(m.start(), p, "-")

    return pos_map

# ── Core analysis ──────────────────────────────────────────────────────────────
def analyze_modonly(path: str, pos_map: dict,
                    min_prob: float = 0.5,
                    require_fail_false: bool = False) -> dict:
    """
    Returns per-(read × motif_start) co-modification stats.
    """
    read_motif = defaultdict(lambda: defaultdict(set))
    with open(path) as f:
        f.readline()  # header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            read_id = parts[0]
            ref_pos_str = parts[1]
            prob_str = parts[2]
            call_code = parts[3]
            fail = parts[4].strip() if len(parts) > 4 else "false"

            if call_code not in (CALL_6mA, CALL_4mC):
                continue
            if require_fail_false and fail.lower() != "false":
                continue
            try:
                ref_pos = int(ref_pos_str)
                prob = float(prob_str)
            except ValueError:
                continue
            if prob < min_prob:
                continue
            if ref_pos not in pos_map:
                continue
            for (ms, pip, strand, em) in pos_map[ref_pos]:
                if em == "6mA" and call_code == CALL_6mA:
                    read_motif[read_id][ms].add("6mA")
                elif em == "4mC" and call_code == CALL_4mC:
                    read_motif[read_id][ms].add("4mC")

    pairs_total = pairs_6mA = pairs_4mC = pairs_both = 0
    for rid, md in read_motif.items():
        for ms, mods in md.items():
            h6 = "6mA" in mods
            h4 = "4mC" in mods
            pairs_total += 1
            pairs_6mA += h6
            pairs_4mC += h4
            pairs_both += h6 and h4

    return dict(pairs_total=pairs_total,
                pairs_6mA=pairs_6mA,
                pairs_4mC=pairs_4mC,
                pairs_both=pairs_both)

def pool_results(results: list) -> dict:
    keys = ("pairs_total", "pairs_6mA", "pairs_4mC", "pairs_both")
    return {k: sum(r[k] for r in results) for k in keys}

def comod_rate(r: dict) -> float:
    return r["pairs_both"] / r["pairs_4mC"] if r["pairs_4mC"] > 0 else 0.0

def ztest_vs_null(r: dict, p_null: float):
    nb = r["pairs_both"]
    n4 = r["pairs_4mC"]
    if n4 == 0:
        return float("nan"), float("nan")
    z, p = proportions_ztest(nb, n4, p_null)
    return z, p

def fmt_p(p):
    if p < 1e-100:
        return "<1e-100"
    return f"{p:.2e}"

# ── GCCGGC flanking-6mA position map ──────────────────────────────────────────
def build_gccggc_pos_map(genome: str) -> dict:
    """
    GCCGGC: pos0=G, pos1=C, pos2=C, pos3=G, pos4=G, pos5=C
    4mC positions: pos1 (C), pos2 (C), pos5 (C)
    6mA "partner" positions: A residues in flanking ±4bp window around each motif
    (motif itself has no A; flanking sequence does)
    We register all A positions at genome[motif_start-4 : motif_start+10] (= ±4bp of 7bp motif)
    """
    MOTIF = "GCCGGC"
    motif_len = len(MOTIF)
    rc_motif = rc(MOTIF)
    mod4mC = {1, 2, 5}   # C positions in GCCGGC
    FLANK = 4

    pos_map = {}

    def reg_4mC(ms, pip, strand):
        rp = ms + pip if strand == "+" else ms + (motif_len - 1 - pip)
        pos_map.setdefault(rp, [])
        pos_map[rp].append((ms, pip, strand, "4mC"))

    def reg_flanking_6mA(ms, strand):
        # scan genome for A positions in [ms-FLANK, ms+motif_len+FLANK)
        start = max(0, ms - FLANK)
        end = min(len(genome), ms + motif_len + FLANK)
        for i in range(start, end):
            # on + strand: check genome[i] == A
            # on - strand: check genome[i] == T (complement of A)
            base = genome[i]
            if strand == "+" and base == "A":
                pos_map.setdefault(i, [])
                pos_map[i].append((ms, i - ms, strand, "6mA"))
            elif strand == "-" and base == "T":
                # the - strand has an A at this ref pos
                pos_map.setdefault(i, [])
                pos_map[i].append((ms, i - ms, strand, "6mA"))

    for m in re.finditer(f"(?={MOTIF})", genome):
        ms = m.start()
        for p in mod4mC:
            reg_4mC(ms, p, "+")
        reg_flanking_6mA(ms, "+")

    for m in re.finditer(f"(?={rc_motif})", genome):
        ms = m.start()
        for p in mod4mC:
            reg_4mC(ms, p, "-")
        reg_flanking_6mA(ms, "-")

    return pos_map

# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    OUTDIR = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/65_per_read_comod"
    T1_FILES = [f"{OUTDIR}/{s}_modonly.tsv" for s in ("1-1", "1-2", "1-3")]

    print("Loading reference...", flush=True)
    genome = load_reference(REF_FA)
    print(f"  Genome: {len(genome):,} bp", flush=True)

    # AAGCCCG position map
    print("Building AAGCCCG position map...", flush=True)
    aag_pos_map = build_pos_map(
        genome, "AAGCCCG",
        mod4mC_positions={3, 5},
        mod6mA_positions={0, 1}
    )

    # GCCGGC flanking position map
    print("Building GCCGGC flanking position map...", flush=True)
    gcc_pos_map = build_gccggc_pos_map(genome)

    # ── Analysis 1: fail filter effect ────────────────────────────────────────
    print("\n[1/3] fail filter sensitivity (AAGCCCG, T1 pooled)...", flush=True)
    res_nofail = pool_results([
        analyze_modonly(f, aag_pos_map, min_prob=0.5, require_fail_false=False)
        for f in T1_FILES
    ])
    res_fail_strict = pool_results([
        analyze_modonly(f, aag_pos_map, min_prob=0.5, require_fail_false=True)
        for f in T1_FILES
    ])

    # ── Analysis 2: probability threshold sensitivity ─────────────────────────
    print("[2/3] Probability threshold sensitivity (AAGCCCG, T1 pooled)...", flush=True)
    thresh_results = {}
    for thr in (0.5, 0.7, 0.9):
        thresh_results[thr] = pool_results([
            analyze_modonly(f, aag_pos_map, min_prob=thr, require_fail_false=False)
            for f in T1_FILES
        ])

    # ── Analysis 3: GCCGGC negative control ───────────────────────────────────
    print("[3/3] GCCGGC negative control (T1 pooled)...", flush=True)
    gcc_results = pool_results([
        analyze_modonly(f, gcc_pos_map, min_prob=0.5, require_fail_false=False)
        for f in T1_FILES
    ])

    # AAGCCCG T1 pooled (prob=0.5, no fail filter) — reference
    aag_ref = res_nofail
    p_null_T1 = 0.3161  # pileup 6mA occupancy T1

    # ── Print results ──────────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("SENSITIVITY ANALYSIS RESULTS")
    print("=" * 70)

    # ── 1. fail filter ────────────────────────────────────────────────────────
    print("\n━━━ Analysis 1: fail=false フィルタ効果 (AAGCCCG, T1 pooled) ━━━")
    print(f"{'条件':<30} {'4mCペア':<10} {'同時修飾数':<12} {'同時修飾率':<12} z値 / p値")
    print("-" * 75)
    for label, r in [("全コール (fail無視, prob≥0.5)", res_nofail),
                      ("fail=false のみ (prob≥0.5)", res_fail_strict)]:
        rate = comod_rate(r)
        z, p = ztest_vs_null(r, p_null_T1)
        print(f"{label:<30} {r['pairs_4mC']:<10,} {r['pairs_both']:<12,} "
              f"{rate*100:<12.1f} z={z:.2f}, p={fmt_p(p)}")

    print(f"\n  fail=false フィルタ適用後の 4mC コール保持率: "
          f"{res_fail_strict['pairs_4mC'] / max(res_nofail['pairs_4mC'], 1) * 100:.1f}%")
    diff = comod_rate(res_fail_strict) - comod_rate(res_nofail)
    print(f"  同時修飾率の変化: {comod_rate(res_nofail)*100:.1f}% → "
          f"{comod_rate(res_fail_strict)*100:.1f}% (Δ={diff*100:+.1f}%)")

    # ── 2. threshold sensitivity ──────────────────────────────────────────────
    print("\n━━━ Analysis 2: call_prob 閾値感度解析 (AAGCCCG, T1 pooled) ━━━")
    print(f"{'閾値':<8} {'4mCペア':<10} {'同時修飾数':<12} {'同時修飾率':<12} z値 / p値")
    print("-" * 60)
    for thr in (0.5, 0.7, 0.9):
        r = thresh_results[thr]
        rate = comod_rate(r)
        z, p = ztest_vs_null(r, p_null_T1)
        print(f"≥{thr:<7} {r['pairs_4mC']:<10,} {r['pairs_both']:<12,} "
              f"{rate*100:<12.1f} z={z:.2f}, p={fmt_p(p)}")

    # retention rates
    base_4mC = thresh_results[0.5]["pairs_4mC"]
    print(f"\n  4mC ペア保持率: prob≥0.7 = "
          f"{thresh_results[0.7]['pairs_4mC']/base_4mC*100:.1f}%,  "
          f"prob≥0.9 = {thresh_results[0.9]['pairs_4mC']/base_4mC*100:.1f}%")

    # ── 3. GCCGGC control ─────────────────────────────────────────────────────
    print("\n━━━ Analysis 3: GCCGGCネガティブコントロール vs AAGCCCG (T1 pooled) ━━━")
    print(f"{'モチーフ':<12} {'4mCペア':<10} {'6mA共局在':<12} {'共局在率':<12} "
          f"{'比較':<25}")
    print("-" * 72)

    aag_rate = comod_rate(aag_ref)
    gcc_rate = comod_rate(gcc_results)

    print(f"{'AAGCCCG':<12} {aag_ref['pairs_4mC']:<10,} {aag_ref['pairs_both']:<12,} "
          f"{aag_rate*100:<12.1f} {'(モチーフ内A, 3-4bp)'}")
    print(f"{'GCCGGC':<12} {gcc_results['pairs_4mC']:<10,} {gcc_results['pairs_both']:<12,} "
          f"{gcc_rate*100:<12.1f} {'(フランキング±4bp A)'}")

    # Fisher's exact test: 2×2 table
    # AAGCCCG: both=X, 4mC_only=Y
    # GCCGGC:  both=A, 4mC_only=B
    a = aag_ref["pairs_both"]
    b = aag_ref["pairs_4mC"] - a
    c = gcc_results["pairs_both"]
    d = gcc_results["pairs_4mC"] - c
    if d >= 0 and b >= 0:
        odds, p_fish = fisher_exact([[a, b], [c, d]])
        print(f"\n  Fisher's exact test (AAGCCCG vs GCCGGC 6mA共局在率):")
        print(f"    OR = {odds:.2f},  p = {fmt_p(p_fish)}")
        print(f"  倍率差: AAGCCCG/GCCGGC = {aag_rate/max(gcc_rate,1e-9):.1f}×")
    else:
        print("\n  (Fisher test skipped: negative cell)")

    # ── Summary Table ──────────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("SUMMARY (論文記述用)")
    print("=" * 70)
    print(f"[A1] fail=false限定: 同時修飾率 {comod_rate(res_fail_strict)*100:.1f}% "
          f"(全コール: {comod_rate(res_nofail)*100:.1f}%) → "
          f"{'頑健' if abs(comod_rate(res_fail_strict) - comod_rate(res_nofail)) < 0.10 else '要注意'}")
    rates = {thr: comod_rate(thresh_results[thr]) for thr in (0.5, 0.7, 0.9)}
    print(f"[A2] 閾値感度: prob≥0.5={rates[0.5]*100:.1f}%, "
          f"≥0.7={rates[0.7]*100:.1f}%, ≥0.9={rates[0.9]*100:.1f}% → "
          f"{'頑健' if max(rates.values()) - min(rates.values()) < 0.10 else '感度あり'}")
    print(f"[A3] GCCGGC対照: {gcc_rate*100:.1f}% vs AAGCCCG {aag_rate*100:.1f}% → "
          f"{aag_rate/max(gcc_rate,1e-9):.1f}× 差 (p={fmt_p(p_fish)})")


if __name__ == "__main__":
    main()
