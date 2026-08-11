#!/usr/bin/env python3
"""
AAGCCCG dual-modification threshold sensitivity (Reviewer Major-2-1).

Two analyses, both at T1, both at probability thresholds 0.50–0.90 by 0.05.

(a) PER-READ pair analysis  (unit = read × motif_instance × strand):
       C3/C5 4mC pass       C3/C5 4mC fail
  A pass    n_both         n_A_only
  A fail    n_C_only       n_neither
    Tests molecular co-occurrence: do the two modifications appear on the
    same DNA molecule more often than expected from marginals?

(b) SITE-LEVEL per-motif-instance analysis  (unit = AAGCCCG motif instance):
    A motif instance is "A_pass" if ANY read covering it has 6mA at A0 or A1
    with prob >= t.  Similarly C_pass for 4mC at C3 or C5.
    Builds the same 4-cell table over the 1,334 motif instances.
    Closer in spirit to the paper's site-level OR=138,440.

Pair/site denominator: include any pair/site that has at least one modonly
row at any of {A0,A1,C3,C5} regardless of probability.  Failing rows
(fail=true) are kept.  Set FILTER_FAIL=True to drop them.

For each (analysis, threshold):
  - Fisher's exact OR + 95% CI (log-OR normal approx, Haldane 0.5 if any
    cell is 0)
  - Fisher's exact two-sided p
  - co-mod rate = n_both / total

Outputs:
  - threshold_sensitivity.tsv         (per-read pairs, requested format)
  - threshold_sensitivity_site.tsv    (per-motif-instance, site-level)
  - threshold_vs_logOR.pdf            (both analyses, side by side)
"""
from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact

# ── Config ─────────────────────────────────────────────────────────────────────
REF_FA      = "/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa"
MODONLY_DIR = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/65_per_read_comod")
OUT_DIR     = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/73_comod_threshold_ROC")
T1_SAMPLES  = ["1-1", "1-2", "1-3"]

THRESHOLDS  = [0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90]

CALL_6mA = "a"
CALL_4mC = "21839"

MOTIF      = "AAGCCCG"
A_POS_SET  = {0, 1}        # 6mA at A0, A1
C_POS_SET  = {3, 5}        # 4mC at C3, C5

FILTER_FAIL = False         # set True to drop fail=true rows

# ── Reference + motif map ──────────────────────────────────────────────────────
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
    """ref_pos -> list of (motif_start, position_in_motif, strand, mod_class)
       where mod_class is 'A' (6mA target) or 'C' (4mC target).
    """
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


# ── Per-read classification ────────────────────────────────────────────────────
def collect_pair_and_site_calls(modonly_paths: list[Path], pos_map: dict):
    """
    Returns:
      pair_max: {(read_id, motif_start, strand): {'A': maxp_A, 'C': maxp_C}}
          per-read max prob at any A and C class positions of that motif.
      site_calls: {(motif_start, strand): {'A_pass':[probs], 'A_total':int,
                                           'C_pass':[probs], 'C_total':int}}
          per-motif aggregated read counts.  At threshold t, the site is
          A_majority if (#A probs >= t) / A_total > 0.5; C similarly.
          Total counts include all rows (passing or failing).
    """
    pair: dict[tuple, dict] = defaultdict(lambda: {"A": -np.inf, "C": -np.inf})
    site: dict[tuple, dict] = defaultdict(
        lambda: {"A_probs": [], "A_total": 0, "C_probs": [], "C_total": 0})

    for p in modonly_paths:
        if not p.exists():
            print(f"  WARNING: missing {p.name}")
            continue
        print(f"  reading {p.name} ...", flush=True)
        with open(p) as f:
            f.readline()  # header
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
                if FILTER_FAIL and parts[4].strip().lower() == "true":
                    continue

                call_code = parts[3]
                read_id   = parts[0]
                for (motif_start, _, strand, cls) in pos_map[ref_pos]:
                    pkey = (read_id, motif_start, strand)
                    skey = (motif_start, strand)
                    pdict = pair[pkey]
                    sdict = site[skey]
                    if cls == "A" and call_code == CALL_6mA:
                        if prob > pdict["A"]:
                            pdict["A"] = prob
                        sdict["A_probs"].append(prob)
                        sdict["A_total"] += 1
                    elif cls == "C" and call_code == CALL_4mC:
                        if prob > pdict["C"]:
                            pdict["C"] = prob
                        sdict["C_probs"].append(prob)
                        sdict["C_total"] += 1
    return pair, site


def classify_site_at_threshold(site_calls: dict, threshold: float,
                               freq_cut: float = 0.5):
    """
    A motif instance is A_pass if (#A probs >= t) / A_total > freq_cut,
    similarly C_pass.  Sites with A_total == 0 or C_total == 0 are excluded
    from the relevant arm (treated as not_pass).
    """
    n_both = n_A_only = n_C_only = n_neither = 0
    for s in site_calls.values():
        a_pass = (s["A_total"] > 0 and
                  sum(p >= threshold for p in s["A_probs"]) / s["A_total"] > freq_cut)
        c_pass = (s["C_total"] > 0 and
                  sum(p >= threshold for p in s["C_probs"]) / s["C_total"] > freq_cut)
        if a_pass and c_pass:
            n_both += 1
        elif a_pass:
            n_A_only += 1
        elif c_pass:
            n_C_only += 1
        else:
            n_neither += 1
    return n_both, n_A_only, n_C_only, n_neither


def classify_at_threshold(pair_max: dict, threshold: float) -> tuple[int, int, int, int]:
    n_both = n_A_only = n_C_only = n_neither = 0
    for d in pair_max.values():
        a_pass = d["A"] >= threshold
        c_pass = d["C"] >= threshold
        if a_pass and c_pass:
            n_both += 1
        elif a_pass:
            n_A_only += 1
        elif c_pass:
            n_C_only += 1
        else:
            n_neither += 1
    return n_both, n_A_only, n_C_only, n_neither


# ── Statistics ─────────────────────────────────────────────────────────────────
def or_ci_p(n_both: int, n_A_only: int, n_C_only: int, n_neither: int):
    """
    Returns (OR, ci_lo, ci_hi, fisher_p_two_sided).
    2x2 layout for OR = (n_both * n_neither) / (n_A_only * n_C_only).
    Uses Haldane 0.5 correction for log-OR SE if any cell is 0.
    """
    table = np.array([[n_both, n_A_only],
                      [n_C_only, n_neither]], dtype=float)

    odds, p_two = fisher_exact(table.astype(int), alternative="two-sided")

    cells = table.flatten()
    if (cells == 0).any():
        a, b, c, d = cells + 0.5
    else:
        a, b, c, d = cells
    if (n_A_only == 0 and n_C_only == 0):
        # OR -> infinity; CI undefined
        return float(odds), float("nan"), float("nan"), float(p_two)

    log_or = np.log((a * d) / (b * c))
    se     = float(np.sqrt(1/a + 1/b + 1/c + 1/d))
    lo     = float(np.exp(log_or - 1.96 * se))
    hi     = float(np.exp(log_or + 1.96 * se))
    return float(odds), lo, hi, float(p_two)


# ── Main ───────────────────────────────────────────────────────────────────────
def main() -> None:
    print("=" * 72)
    print("AAGCCCG per-read co-modification — threshold sensitivity (T1)")
    print("=" * 72)
    print(f"  thresholds: {THRESHOLDS}")
    print(f"  filter fail rows: {FILTER_FAIL}")

    print("\n[1/4] loading reference...")
    genome = load_reference(REF_FA)
    print(f"  genome: {len(genome):,} bp")

    print("\n[2/4] building motif position map...")
    pos_map = build_pos_map(genome)
    n_motif_starts = len({(ms, st)
                          for entries in pos_map.values()
                          for (ms, _, st, _) in entries})
    print(f"  AAGCCCG motif instances: {n_motif_starts:,}")
    print(f"  reference positions tracked: {len(pos_map):,}")

    print("\n[3/4] collecting per-read & per-site calls from T1 modonly files...")
    modonly_paths = [MODONLY_DIR / f"{s}_modonly.tsv" for s in T1_SAMPLES]
    pair_max, site_calls = collect_pair_and_site_calls(modonly_paths, pos_map)
    print(f"  total (read × motif × strand) pairs: {len(pair_max):,}")
    print(f"  total motif-instance × strand sites: {len(site_calls):,}")

    print("\n[4/4] classifying at each threshold + computing OR/CI/p ...")

    def sweep(label: str, classify_fn) -> list[dict]:
        rows: list[dict] = []
        print(f"\n  --- {label} ---")
        for t in THRESHOLDS:
            n_both, n_A_only, n_C_only, n_neither = classify_fn(t)
            total = n_both + n_A_only + n_C_only + n_neither
            OR, ci_lo, ci_hi, p_two = or_ci_p(n_both, n_A_only, n_C_only, n_neither)
            comod_rate = n_both / total if total else float("nan")
            rows.append(dict(
                threshold=round(t, 2),
                n_both=n_both,
                n_A_only=n_A_only,
                n_C_only=n_C_only,
                n_neither=n_neither,
                OR=OR,
                CI_lower=ci_lo,
                CI_upper=ci_hi,
                p_value=p_two,
                co_mod_rate=comod_rate,
            ))
            print(f"  t={t:.2f}: both={n_both:>7,d} A_only={n_A_only:>7,d} "
                  f"C_only={n_C_only:>7,d} neither={n_neither:>7,d} | "
                  f"OR={OR:9.3g}  CI=[{ci_lo:9.3g}, {ci_hi:9.3g}]  "
                  f"p={p_two:.2e}  comod%={comod_rate*100:.2f}")
        return rows

    rows_pair = sweep(
        "PER-READ pair (read × motif_instance)",
        lambda t: classify_at_threshold(pair_max, t),
    )
    rows_site = sweep(
        "SITE-LEVEL per-motif-instance (frequency >= 0.5)",
        lambda t: classify_site_at_threshold(site_calls, t, freq_cut=0.5),
    )

    # ── Save TSVs ─────────────────────────────────────────────────────────────
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_pair = OUT_DIR / "threshold_sensitivity.tsv"
    out_site = OUT_DIR / "threshold_sensitivity_site.tsv"
    pd.DataFrame(rows_pair).to_csv(out_pair, sep="\t", index=False, float_format="%.6g")
    pd.DataFrame(rows_site).to_csv(out_site, sep="\t", index=False, float_format="%.6g")
    print(f"\nTSV saved (per-read): {out_pair}")
    print(f"TSV saved (site-level): {out_site}")

    # ── PDF plot: threshold vs log10(OR) for both analyses ────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.7), sharex=True)

    for ax, rows, title, color in [
        (axes[0], rows_pair, "Per-read pair OR", "#C0392B"),
        (axes[1], rows_site, "Site-level OR (freq ≥ 0.5)", "#2C5AA0"),
    ]:
        df = pd.DataFrame(rows)
        OR  = df["OR"].astype(float).values
        lo  = df["CI_lower"].astype(float).values
        hi  = df["CI_upper"].astype(float).values

        with np.errstate(invalid="ignore", divide="ignore"):
            log_or = np.log10(np.where(np.isfinite(OR) & (OR > 0), OR, np.nan))
            log_lo = np.log10(np.where(np.isfinite(lo) & (lo > 0), lo, np.nan))
            log_hi = np.log10(np.where(np.isfinite(hi) & (hi > 0), hi, np.nan))

        # Cap CIs at finite values for plotting; mark unbounded with arrows.
        finite = np.isfinite(log_or)
        thr_finite    = df["threshold"].values[finite]
        log_or_finite = log_or[finite]
        # Make sure errors are non-negative; replace nan with 0
        err_lo = np.where(np.isfinite(log_lo[finite]),
                          log_or_finite - log_lo[finite], 0.0)
        err_hi = np.where(np.isfinite(log_hi[finite]),
                          log_hi[finite] - log_or_finite, 0.0)
        err_lo = np.clip(err_lo, 0, None)
        err_hi = np.clip(err_hi, 0, None)
        err = np.vstack([err_lo, err_hi])

        ax.errorbar(thr_finite, log_or_finite, yerr=err, fmt="o-",
                    color=color, ecolor="#888", elinewidth=1.0,
                    capsize=3, markersize=5, linewidth=1.4,
                    label="log₁₀(OR) ± 95% CI")

        # Annotate non-finite ORs (inf or undefined) at top of plot.
        for tval, ORval in zip(df["threshold"].values, OR):
            if not np.isfinite(ORval):
                ax.annotate("∞ / NA", xy=(tval, ax.get_ylim()[1] * 0.9 if ax.get_ylim()[1] else 0),
                            ha="center", fontsize=7, color="#666")

        ax.axhline(0, color="grey", linewidth=0.7, linestyle="--")
        ax.set_xlabel("Per-read modification prob. threshold")
        ax.set_ylabel("log₁₀ Odds Ratio")
        ax.set_title(title)
        ax.set_xticks(THRESHOLDS)
        ax.tick_params(direction="in")
        ax.grid(axis="y", linewidth=0.4, alpha=0.4)
        ax.legend(fontsize=8, loc="best")

    fig.suptitle("AAGCCCG 4mC×6mA OR vs probability threshold (T1)", y=1.02)
    plt.tight_layout()
    pdf_path = OUT_DIR / "threshold_vs_logOR.pdf"
    fig.savefig(pdf_path, dpi=300, bbox_inches="tight")
    print(f"PDF saved: {pdf_path}")


if __name__ == "__main__":
    main()
