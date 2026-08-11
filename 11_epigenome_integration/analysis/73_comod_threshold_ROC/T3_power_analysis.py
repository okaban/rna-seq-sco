#!/usr/bin/env python3
"""
Post-hoc power analysis for AAGCCCG dual-modification detection at T3
(Reviewer Major-3-1).

At T3 only n=59 AAGCCCG motif instances pass a depth >= 10 filter.
The T1 observed co-modification rate is ~31.6% (paper text).  We ask:
given a sample of n=59 sites, what is the power to reject the null
hypothesis of marginal independence (H0: p_co_mod = p_A * p_C) when the
true rate is the T1 observed (H1: p_co_mod = 0.316)?

Test framework
--------------
Each AAGCCCG site is a Bernoulli trial: co-mod (1) or not (0).
   H0:  p_co_mod = p_0 = p_A × p_C   (independence)
   H1:  p_co_mod = p_1 = 0.316       (T1 observed)
   alpha = 0.05, two-sided exact binomial.

Power = P(reject H0 | H1).
Critical region: |X − n·p_0| ≥ k* such that two-sided exact tail ≤ alpha.

Also reports the minimum n needed for power ≥ 0.80.

T3 marginals (p_A, p_C) are computed from the T3 modonly TSVs:
  - For each AAGCCCG motif instance × strand, keep sites with both
    A_total ≥ depth_min AND C_total ≥ depth_min.
  - p_A = fraction of those sites with > 50% of A reads at prob ≥ 0.5.
  - p_C = same for C.
  - Falls back to T1 marginals if T3 has < 5 retained sites.
"""
from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binom

# ── Config ─────────────────────────────────────────────────────────────────────
REF_FA      = "/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa"
MODONLY_DIR = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/65_per_read_comod")
PILEUP_DIR  = Path("/Users/okaban/bioinfo/methyl/260102_M145/analysis/cursor_results/pileup")
OUT_DIR     = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/73_comod_threshold_ROC")

T1_SAMPLES  = ["1-1", "1-2", "1-3"]
T3_SAMPLES  = ["3-2", "3-3", "3-4"]

DEPTH_MIN     = 10        # per-site coverage filter (n_valid_cov + n_fail, summed over reps)
PROB_THRESH   = 0.5       # per-read probability cutoff
FREQ_CUT      = 0.5       # site is mod if > FREQ_CUT of reads pass
ALPHA         = 0.05      # two-sided
T1_COMOD_RATE = 0.316     # H1 (paper-stated T1 observed)

CALL_6mA = "a"
CALL_4mC = "21839"
MOTIF      = "AAGCCCG"
A_POS_SET  = {0, 1}
C_POS_SET  = {3, 5}


# ── Reference + motif map (same as threshold script) ──────────────────────────
def rc(seq: str) -> str:
    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))


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


def build_pos_map(genome: str) -> dict:
    pos_map: dict[int, list] = {}
    motif_len = len(MOTIF)
    rc_motif = rc(MOTIF)
    targets = A_POS_SET | C_POS_SET

    def reg(motif_start: int, p: int, strand: str) -> None:
        ref_pos = motif_start + p if strand == "+" else motif_start + (motif_len - 1 - p)
        cls = "A" if p in A_POS_SET else "C"
        pos_map.setdefault(ref_pos, []).append((motif_start, p, strand, cls))

    for m in re.finditer(f"(?={MOTIF})", genome):
        for p in targets:
            reg(m.start(), p, "+")
    for m in re.finditer(f"(?={rc_motif})", genome):
        for p in targets:
            reg(m.start(), p, "-")
    return pos_map


# ── Coverage from pileup BED ──────────────────────────────────────────────────
def load_pileup_coverage(pileup_path: Path, target_positions: set,
                         call_code: str) -> dict:
    """{pos: total_cov} where total_cov = N_valid_cov (col10) + N_fail (col16)."""
    cov: dict[int, int] = {}
    with open(pileup_path) as f:
        for line in f:
            parts = line.split("\t")
            if len(parts) < 16:
                continue
            if parts[3] != call_code:
                continue
            pos = int(parts[1])
            if pos not in target_positions:
                continue
            n_valid = int(parts[9])
            n_fail  = int(parts[15])
            cov[pos] = cov.get(pos, 0) + n_valid + n_fail   # accumulate over reps
    return cov


def collect_pileup_coverage(samples: list[str], pos_map: dict) -> dict:
    """
    Returns {(motif_start, strand): {'A_cov': int, 'C_cov': int}}
    summed across replicates, where A_cov is the minimum coverage across
    the A-class positions (most restrictive).
    """
    a_pos_set = set()
    c_pos_set = set()
    for ref_pos, entries in pos_map.items():
        for (_, _, _, cls) in entries:
            if cls == "A":
                a_pos_set.add(ref_pos)
            else:
                c_pos_set.add(ref_pos)

    a_cov_total = defaultdict(int)
    c_cov_total = defaultdict(int)
    for s in samples:
        bed = PILEUP_DIR / f"{s}_pileup.bed"
        if not bed.exists():
            print(f"  WARNING: missing {bed.name}")
            continue
        print(f"  pileup {bed.name} ...", flush=True)
        a_cov = load_pileup_coverage(bed, a_pos_set, CALL_6mA)
        c_cov = load_pileup_coverage(bed, c_pos_set, CALL_4mC)
        for pos, c in a_cov.items():
            a_cov_total[pos] += c
        for pos, c in c_cov.items():
            c_cov_total[pos] += c

    site: dict[tuple, dict] = {}
    for ref_pos, entries in pos_map.items():
        for (motif_start, _, strand, cls) in entries:
            skey = (motif_start, strand)
            if skey not in site:
                site[skey] = {"A_cov_min": float("inf"), "C_cov_min": float("inf")}
            if cls == "A":
                cov = a_cov_total.get(ref_pos, 0)
                site[skey]["A_cov_min"] = min(site[skey]["A_cov_min"], cov)
            else:
                cov = c_cov_total.get(ref_pos, 0)
                site[skey]["C_cov_min"] = min(site[skey]["C_cov_min"], cov)
    # Replace inf with 0 (no positions of that class — shouldn't happen for AAGCCCG)
    for s in site.values():
        if s["A_cov_min"] == float("inf"):
            s["A_cov_min"] = 0
        if s["C_cov_min"] == float("inf"):
            s["C_cov_min"] = 0
    return site


# ── Site call collection per timepoint ────────────────────────────────────────
def collect_site_calls(modonly_paths: list[Path], pos_map: dict) -> dict:
    site: dict[tuple, dict] = defaultdict(
        lambda: {"A_pass": 0, "A_total": 0, "C_pass": 0, "C_total": 0})
    for p in modonly_paths:
        if not p.exists():
            print(f"  WARNING: missing {p.name}")
            continue
        print(f"  reading {p.name} ...", flush=True)
        with open(p) as f:
            f.readline()
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                try:
                    ref_pos = int(parts[1])
                    prob    = float(parts[2])
                except ValueError:
                    continue
                if ref_pos not in pos_map:
                    continue
                call_code = parts[3]
                for (motif_start, _, strand, cls) in pos_map[ref_pos]:
                    skey = (motif_start, strand)
                    s = site[skey]
                    if cls == "A" and call_code == CALL_6mA:
                        s["A_total"] += 1
                        if prob >= PROB_THRESH:
                            s["A_pass"] += 1
                    elif cls == "C" and call_code == CALL_4mC:
                        s["C_total"] += 1
                        if prob >= PROB_THRESH:
                            s["C_pass"] += 1
    return site


def filtered_marginals(site_calls: dict, site_cov: dict,
                       depth_min: int, freq_cut: float):
    """
    Returns (p_A, p_C, p_co_mod, n_sites) where a site is included if
    BOTH A_cov_min >= depth_min AND C_cov_min >= depth_min (from pileup).
    A site is A_mod if A_pass / max(A_total, 1) > freq_cut.
    """
    n_A = n_C = n_both = n = 0
    for skey, c in site_cov.items():
        if c["A_cov_min"] < depth_min or c["C_cov_min"] < depth_min:
            continue
        s = site_calls.get(skey, {"A_pass": 0, "A_total": 0,
                                  "C_pass": 0, "C_total": 0})
        a_freq = s["A_pass"] / s["A_total"] if s["A_total"] else 0.0
        c_freq = s["C_pass"] / s["C_total"] if s["C_total"] else 0.0
        n += 1
        a_mod = a_freq > freq_cut
        c_mod = c_freq > freq_cut
        n_A   += a_mod
        n_C   += c_mod
        n_both += a_mod and c_mod
    if n == 0:
        return float("nan"), float("nan"), float("nan"), 0
    return n_A / n, n_C / n, n_both / n, n


# ── Power calculation ─────────────────────────────────────────────────────────
def two_sided_critical(n: int, p0: float, alpha: float):
    """
    Returns the acceptance interval [k_low, k_high] for an exact two-sided
    binomial test of H0: p == p0 with target two-sided alpha.
    Uses equal-tailed acceptance (alpha/2 per tail).
    """
    a2 = alpha / 2
    # k_low = largest k with P(X < k | p0) <= a2  ⇒  X >= k_low + 1 in lower tail
    # k_high = smallest k with P(X > k | p0) <= a2
    k_low  = int(binom.ppf(a2, n, p0))
    while k_low > 0 and binom.cdf(k_low - 1, n, p0) >= a2:
        k_low -= 1
    k_high = int(binom.ppf(1 - a2, n, p0))
    while k_high < n and binom.sf(k_high, n, p0) >= a2:
        k_high += 1
    return k_low, k_high


def exact_two_sided_power(n: int, p0: float, p1: float, alpha: float) -> float:
    """
    Power = P(reject H0 | X ~ Binomial(n, p1)) for an exact two-sided
    binomial test of H0: p == p0.  Equal-tailed alpha/2 per tail.
    """
    k_low, k_high = two_sided_critical(n, p0, alpha)
    # Reject if X < k_low or X > k_high
    p_lo = binom.cdf(k_low - 1, n, p1) if k_low - 1 >= 0 else 0.0
    p_hi = binom.sf(k_high, n, p1)
    return float(p_lo + p_hi)


def min_n_for_power(p0: float, p1: float, alpha: float, target: float = 0.80,
                    n_max: int = 5000) -> int:
    for n in range(2, n_max + 1):
        if exact_two_sided_power(n, p0, p1, alpha) >= target:
            return n
    return -1  # not reached


# ── Main ───────────────────────────────────────────────────────────────────────
def main() -> None:
    print("=" * 72)
    print("T3 post-hoc power analysis for AAGCCCG co-modification detection")
    print("=" * 72)
    print(f"  alpha = {ALPHA}  two-sided exact binomial")
    print(f"  H1 (T1 observed) co-mod rate = {T1_COMOD_RATE}")
    print(f"  depth filter = {DEPTH_MIN}, prob threshold = {PROB_THRESH}, "
          f"site freq cut = {FREQ_CUT}")

    print("\n[1/3] loading reference + motif positions ...")
    genome = load_reference(REF_FA)
    pos_map = build_pos_map(genome)
    print(f"  reference positions tracked: {len(pos_map):,}")

    print("\n[2a/3] computing T3 pileup-based site coverage ...")
    t3_cov = collect_pileup_coverage(T3_SAMPLES, pos_map)
    n_t3_pass = sum(1 for c in t3_cov.values()
                    if c["A_cov_min"] >= DEPTH_MIN and c["C_cov_min"] >= DEPTH_MIN)
    print(f"  T3 motif instances with min A&C coverage >= {DEPTH_MIN}: n = {n_t3_pass}")

    print("\n[2b/3] computing T3 marginals on those sites ...")
    t3_paths = [MODONLY_DIR / f"{s}_modonly.tsv" for s in T3_SAMPLES]
    t3_site_calls = collect_site_calls(t3_paths, pos_map)
    p_A_t3, p_C_t3, p_comod_t3, n_t3 = filtered_marginals(
        t3_site_calls, t3_cov, DEPTH_MIN, FREQ_CUT)
    print(f"  T3 sites used: n = {n_t3}")
    print(f"  p_A_t3 = {p_A_t3:.4f}")
    print(f"  p_C_t3 = {p_C_t3:.4f}")
    print(f"  p_co_mod (observed at T3, depth-filtered) = {p_comod_t3:.4f}")

    if n_t3 < 5:
        print("  WARNING: too few T3 sites — falling back to T1 marginals")
        t1_cov = collect_pileup_coverage(T1_SAMPLES, pos_map)
        t1_paths = [MODONLY_DIR / f"{s}_modonly.tsv" for s in T1_SAMPLES]
        t1_site_calls = collect_site_calls(t1_paths, pos_map)
        p_A_t1, p_C_t1, _, _ = filtered_marginals(
            t1_site_calls, t1_cov, DEPTH_MIN, FREQ_CUT)
        p_A_used = p_A_t1
        p_C_used = p_C_t1
        marg_src = "T1"
    else:
        p_A_used = p_A_t3
        p_C_used = p_C_t3
        marg_src = "T3"

    p_0 = p_A_used * p_C_used
    p_1 = T1_COMOD_RATE
    print(f"\n  using marginals from: {marg_src}")
    print(f"  H0 expected co-mod rate (independence) = {p_0:.4f}")
    print(f"  H1 expected co-mod rate (observed T1)  = {p_1:.4f}")

    print("\n[3/3] computing power ...")
    n_used = 59
    power_n59 = exact_two_sided_power(n_used, p_0, p_1, ALPHA)
    print(f"  Power at n = {n_used} : {power_n59:.4f}  ({power_n59*100:.2f}%)")

    min_n = min_n_for_power(p_0, p_1, ALPHA, target=0.80)
    if min_n > 0:
        power_at_min = exact_two_sided_power(min_n, p_0, p_1, ALPHA)
        print(f"  Minimum n for power ≥ 0.80 : n = {min_n}  "
              f"(power = {power_at_min:.4f})")
    else:
        print(f"  Minimum n for power ≥ 0.80 : not reached within search range")

    # ── Save TSV ──────────────────────────────────────────────────────────────
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_tsv = OUT_DIR / "T3_power_analysis.tsv"
    rows = [
        ("alpha",                           ALPHA),
        ("test",                            "two-sided exact binomial"),
        ("H1_p_co_mod (T1 observed)",       T1_COMOD_RATE),
        ("depth_min_filter",                DEPTH_MIN),
        ("prob_threshold",                  PROB_THRESH),
        ("freq_cut",                        FREQ_CUT),
        ("marginals_source",                marg_src),
        ("p_A (T3 if available)",           p_A_t3),
        ("p_C (T3 if available)",           p_C_t3),
        ("p_co_mod_observed_T3_filtered",   p_comod_t3),
        ("n_T3_sites_depth_filtered",       n_t3),
        ("p_A_used_for_H0",                 p_A_used),
        ("p_C_used_for_H0",                 p_C_used),
        ("H0_p_co_mod (independence)",      p_0),
        ("n_for_power_calc",                n_used),
        ("power_at_n59",                    power_n59),
        ("min_n_for_power_0_80",            min_n),
    ]
    with open(out_tsv, "w") as fh:
        fh.write("metric\tvalue\n")
        for k, v in rows:
            if isinstance(v, float):
                fh.write(f"{k}\t{v:.6f}\n")
            else:
                fh.write(f"{k}\t{v}\n")
    print(f"\nTSV saved: {out_tsv}")

    # Sweep table for context
    print("\nPower curve over a range of n:")
    print("    n     power")
    sweep_rows = []
    for n in [10, 20, 30, 40, 50, 59, 75, 100, 150, 200, 300, 500, 1000]:
        pw = exact_two_sided_power(n, p_0, p_1, ALPHA)
        print(f"  {n:>4}    {pw:.4f}")
        sweep_rows.append({"n": n, "power": pw})
    pd.DataFrame(sweep_rows).to_csv(OUT_DIR / "T3_power_curve.tsv", sep="\t",
                                    index=False, float_format="%.6f")
    print(f"\nPower curve saved: {OUT_DIR / 'T3_power_curve.tsv'}")


if __name__ == "__main__":
    main()
