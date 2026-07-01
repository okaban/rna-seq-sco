#!/usr/bin/env python3
"""
Threshold sensitivity analysis for Nanopore methylation calling.

Reviewer concern (Major-2-1): "How does the per-read probability threshold
choice affect results?"

Computes at thresholds 0.5, 0.6, 0.75, 0.8, 0.9:
  1. GCCGGC mean site occupancy (T1/T2/T3)
  2. AAGCCCG mean site occupancy (T1/T2/T3)
  3. Protection zone classifier AUC (T1, 5-fold CV)
  4. OR for AAGCCCG dual 4mC+6mA modification (T1)

Method:
  - Occupancy: N_mod(T) / total_coverage per canonical motif site,
    where N_mod(T) counts per-read calls with call_prob >= T (from modonly.tsv)
    and total_coverage = N_valid_cov + N_fail from modkit pileup BED.
  - AUC: Recompute high-confidence methylation sites at each T (n_reps>=2,
    mean_freq>=10%), then recompute nearest_methyl_distance per gene TSS,
    run logistic regression (is_exposed) with 5-fold stratified CV.
  - OR: Per-read co-modification analysis (existing sensitivity_analysis.py
    approach, extended to new thresholds).
"""

import re
import sys
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict
from scipy.stats import fisher_exact
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

# ── Paths ──────────────────────────────────────────────────────────────────
REF_FA      = "/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa"
MODONLY_DIR = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/65_per_read_comod")
PILEUP_DIR  = Path("/Users/okaban/bioinfo/methyl/260102_M145/analysis/cursor_results/pileup")
FEATURES_TSV = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/52_shielded_exposed_boundary/tables/all_genes_features.tsv")
OUT_TSV     = Path("/Users/okaban/obsidian/Research/rna-seq/Writing/Supplementary_threshold_sensitivity.tsv")

# ── Parameters ─────────────────────────────────────────────────────────────
THRESHOLDS  = [0.5, 0.6, 0.75, 0.8, 0.9]
SAMPLES = {
    "T1": ["1-1", "1-2", "1-3"],
    "T2": ["2-1", "2-3", "2-4"],
    "T3": ["3-2", "3-3", "3-4"],
}

CALL_6mA = "a"
CALL_4mC = "21839"

# GCCGGC: 4mC at positions 1,2,5 within 6-mer (0-indexed: G=0,C=1,C=2,G=3,G=4,C=5)
GCCGGC_4mC_POSITIONS = {1, 2, 5}
# AAGCCCG: 6mA at positions 0,1 within 7-mer (A=0,A=1,G=2,C=3,C=4,C=5,G=6)
AAGCCCG_6mA_POSITIONS = {0, 1}

AUC_MIN_FREQ   = 10.0   # site must have mean_freq >= this % to count as methylated
AUC_MIN_REPS   = 2      # site must appear in at least this many T1 replicates
AUC_WINDOW_BP  = 2000   # search window upstream of TSS for nearest methylation
SEED           = 42
N_SPLITS       = 5


# ─────────────────────────────────────────────────────────────────────────────
# 1. Reference genome
# ─────────────────────────────────────────────────────────────────────────────

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


# ─────────────────────────────────────────────────────────────────────────────
# 2. Canonical motif positions in reference genome
# ─────────────────────────────────────────────────────────────────────────────

def build_canonical_positions(genome: str, motif: str, mod_positions: set,
                              label: str) -> dict:
    """
    Returns set of reference positions (0-based) that correspond to
    mod_positions within all canonical occurrences of motif (both strands).
    """
    positions = set()
    motif_len = len(motif)
    rc_motif  = rc(motif)

    for m in re.finditer(f"(?={motif})", genome):
        for p in mod_positions:
            positions.add(m.start() + p)

    for m in re.finditer(f"(?={rc_motif})", genome):
        for p in mod_positions:
            positions.add(m.start() + (motif_len - 1 - p))

    print(f"  Canonical {label} positions: {len(positions):,}", flush=True)
    return positions


# ─────────────────────────────────────────────────────────────────────────────
# 3. Coverage from pileup BED
# ─────────────────────────────────────────────────────────────────────────────

def load_pileup_coverage(pileup_path: Path, target_positions: set,
                         call_code: str) -> dict:
    """
    Returns: {pos: total_coverage} for positions in target_positions.
    total_coverage = N_valid_cov + N_fail (physical read depth).
    """
    cov = {}
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
            cov[pos] = n_valid + n_fail
    return cov


# ─────────────────────────────────────────────────────────────────────────────
# 4. Per-read modification counts from modonly TSV
# ─────────────────────────────────────────────────────────────────────────────

def count_mod_per_threshold(modonly_path: Path, target_positions: set,
                             call_code: str, thresholds: list) -> dict:
    """
    Returns: {threshold: {pos: count}} for positions in target_positions.
    Reads modonly file once for efficiency.
    """
    counts = {t: defaultdict(int) for t in thresholds}
    min_t = min(thresholds)

    with open(modonly_path) as f:
        f.readline()  # header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            if parts[3] != call_code:
                continue
            try:
                pos  = int(parts[1])
                prob = float(parts[2])
            except ValueError:
                continue
            if prob < min_t:
                continue
            if pos not in target_positions:
                continue
            for t in thresholds:
                if prob >= t:
                    counts[t][pos] += 1

    return counts


# ─────────────────────────────────────────────────────────────────────────────
# 5. Per-site per-sample occupancy matrix
# ─────────────────────────────────────────────────────────────────────────────

def compute_occupancy_matrix(samples_list: list,
                              target_positions: set,
                              call_code: str,
                              thresholds: list,
                              label: str) -> dict:
    """
    For each threshold, returns list of site-level methylation frequencies
    (one value per (site, sample) pair that has coverage).

    occupancy(T) = N_mod(T) / (N_valid_cov + N_fail)
    """
    # {threshold: [per-site frequencies]}
    all_freqs = {t: [] for t in thresholds}
    # {threshold: {pos: [per-sample freqs]}} for AUC site matrix
    site_sample = {t: defaultdict(list) for t in thresholds}

    for sample in samples_list:
        pileup_path  = PILEUP_DIR / f"{sample}_pileup.bed"
        modonly_path = MODONLY_DIR / f"{sample}_modonly.tsv"

        if not pileup_path.exists():
            print(f"    WARNING: missing pileup {pileup_path.name}", flush=True)
            continue
        if not modonly_path.exists():
            print(f"    WARNING: missing modonly {modonly_path.name}", flush=True)
            continue

        print(f"    [{label}] coverage: {sample}...", end=" ", flush=True)
        cov = load_pileup_coverage(pileup_path, target_positions, call_code)
        print(f"{len(cov):,} covered sites", flush=True)

        print(f"    [{label}] mod counts: {sample}...", end=" ", flush=True)
        counts = count_mod_per_threshold(modonly_path, target_positions,
                                         call_code, thresholds)
        print("done", flush=True)

        for pos, total_cov in cov.items():
            if total_cov == 0:
                continue
            for t in thresholds:
                n_mod = counts[t].get(pos, 0)
                freq  = n_mod / total_cov * 100.0
                all_freqs[t].append(freq)
                site_sample[t][pos].append(freq)

    mean_occ = {}
    for t in thresholds:
        vals = all_freqs[t]
        mean_occ[t] = np.mean(vals) if vals else float("nan")

    return mean_occ, site_sample


# ─────────────────────────────────────────────────────────────────────────────
# 6. Protection zone AUC
# ─────────────────────────────────────────────────────────────────────────────

def nearest_dist(tss: int, positions: np.ndarray) -> float:
    if len(positions) == 0:
        return float("inf")
    return float(np.min(np.abs(positions - tss)))


def compute_protection_auc(t1_site_sample: dict, threshold: float,
                            genes_df: pd.DataFrame) -> float:
    """
    Recompute protection zone AUC using methylation sites detected at a
    given probability threshold for T1.

    High-confidence sites: positions where mean_freq >= AUC_MIN_FREQ
    across >=AUC_MIN_REPS replicates.
    Returns 5-fold CV AUC.
    """
    # Build high-confidence site positions at this threshold
    site_freqs = t1_site_sample[threshold]
    active_positions = []
    for pos, freqs in site_freqs.items():
        if len(freqs) >= AUC_MIN_REPS and np.mean(freqs) >= AUC_MIN_FREQ:
            active_positions.append(pos)
    active_positions = np.array(sorted(active_positions))

    if len(active_positions) < 10:
        return float("nan")

    # Compute nearest methylation distance for each gene
    genes_copy = genes_df.copy()
    genes_copy["nearest_dist_T"] = genes_copy["tss"].apply(
        lambda tss_val: nearest_dist(int(tss_val), active_positions)
    )

    # Filter to genes with finite distance and valid is_exposed label
    valid = genes_copy[
        genes_copy["nearest_dist_T"].apply(np.isfinite) &
        genes_copy["is_exposed"].notna()
    ].copy()

    if valid["is_exposed"].sum() < 3:
        return float("nan")

    X = (-valid["nearest_dist_T"].values).reshape(-1, 1)  # lower dist → more exposed
    y = valid["is_exposed"].astype(int).values

    skf = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=SEED)
    fold_aucs = []
    for train_idx, test_idx in skf.split(X, y):
        if len(np.unique(y[test_idx])) < 2:
            continue
        lr = LogisticRegression(random_state=SEED, max_iter=1000)
        lr.fit(X[train_idx], y[train_idx])
        prob = lr.predict_proba(X[test_idx])[:, 1]
        fold_aucs.append(roc_auc_score(y[test_idx], prob))

    return float(np.mean(fold_aucs)) if fold_aucs else float("nan")


# ─────────────────────────────────────────────────────────────────────────────
# 7. OR for AAGCCCG dual-modification
#    (per-read co-modification, same logic as sensitivity_analysis.py)
# ─────────────────────────────────────────────────────────────────────────────

def build_aagcccg_pos_map(genome: str) -> dict:
    """
    Map ref_pos -> list of (motif_start, pos_in_motif, strand, mod_type)
    for 4mC positions {3,5} and 6mA positions {0,1} within AAGCCCG.
    """
    MOTIF    = "AAGCCCG"
    motif_len = len(MOTIF)
    rc_motif  = rc(MOTIF)
    mod4mC    = {3, 5}
    mod6mA    = {0, 1}
    target    = mod4mC | mod6mA
    pos_map   = {}

    def reg(ms, pip, strand):
        rp = ms + pip if strand == "+" else ms + (motif_len - 1 - pip)
        pos_map.setdefault(rp, [])
        mt = "6mA" if pip in mod6mA else "4mC"
        pos_map[rp].append((ms, pip, strand, mt))

    for m in re.finditer(f"(?={MOTIF})", genome):
        for p in target:
            reg(m.start(), p, "+")
    for m in re.finditer(f"(?={rc_motif})", genome):
        for p in target:
            reg(m.start(), p, "-")

    return pos_map


def analyze_comod(modonly_path: Path, pos_map: dict, min_prob: float) -> dict:
    read_motif = defaultdict(lambda: defaultdict(set))
    with open(modonly_path) as f:
        f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            if parts[3] not in (CALL_6mA, CALL_4mC):
                continue
            try:
                ref_pos = int(parts[1])
                prob    = float(parts[2])
            except ValueError:
                continue
            if prob < min_prob:
                continue
            if ref_pos not in pos_map:
                continue
            for (ms, pip, strand, em) in pos_map[ref_pos]:
                call_code = parts[3]
                if em == "6mA" and call_code == CALL_6mA:
                    read_motif[parts[0]][ms].add("6mA")
                elif em == "4mC" and call_code == CALL_4mC:
                    read_motif[parts[0]][ms].add("4mC")

    pairs_total = pairs_6mA = pairs_4mC = pairs_both = 0
    for rid, md in read_motif.items():
        for ms, mods in md.items():
            h6 = "6mA" in mods
            h4 = "4mC" in mods
            pairs_total += 1
            pairs_6mA   += h6
            pairs_4mC   += h4
            pairs_both  += h6 and h4

    return dict(pairs_total=pairs_total, pairs_6mA=pairs_6mA,
                pairs_4mC=pairs_4mC, pairs_both=pairs_both)


def pool(results: list) -> dict:
    keys = ("pairs_total", "pairs_6mA", "pairs_4mC", "pairs_both")
    return {k: sum(r[k] for r in results) for k in keys}


def compute_or(r: dict) -> float:
    """Compute odds ratio from contingency table."""
    a = r["pairs_both"]
    b = r["pairs_4mC"] - a           # 4mC only
    c = r["pairs_6mA"] - a           # 6mA only (approx)
    d = r["pairs_total"] - a - b - c  # neither

    if b <= 0 or c <= 0:
        return float("nan")
    try:
        odds, _ = fisher_exact([[a, b], [c, d]])
        return float(odds)
    except Exception:
        return float("nan")


def comod_rate(r: dict) -> float:
    return r["pairs_both"] / r["pairs_4mC"] if r["pairs_4mC"] > 0 else 0.0


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    print("=" * 70, flush=True)
    print("Threshold Sensitivity Analysis", flush=True)
    print("=" * 70, flush=True)

    # ── Load reference ─────────────────────────────────────────────────────
    print("\n[0/4] Loading reference genome...", flush=True)
    genome = load_reference(REF_FA)
    print(f"  Genome: {len(genome):,} bp", flush=True)

    # ── Canonical motif positions ──────────────────────────────────────────
    print("\n[1/4] Building canonical motif positions...", flush=True)
    gccggc_pos  = build_canonical_positions(genome, "GCCGGC",  GCCGGC_4mC_POSITIONS,  "GCCGGC-4mC")
    aagcccg_pos = build_canonical_positions(genome, "AAGCCCG", AAGCCCG_6mA_POSITIONS, "AAGCCCG-6mA")
    aagcccg_map = build_aagcccg_pos_map(genome)  # for OR analysis

    # ── Per-timepoint occupancy ────────────────────────────────────────────
    print("\n[2/4] Computing per-timepoint occupancy...", flush=True)

    gcc_occ  = {}   # {timepoint: {threshold: mean_freq}}
    aag_occ  = {}
    t1_gcc_site_sample  = None
    t1_aag_site_sample  = None

    for tp, samples in SAMPLES.items():
        print(f"\n  Timepoint {tp}:", flush=True)

        gcc_mean, gcc_ss = compute_occupancy_matrix(
            samples, gccggc_pos, CALL_4mC, THRESHOLDS, "GCCGGC")
        aag_mean, aag_ss = compute_occupancy_matrix(
            samples, aagcccg_pos, CALL_6mA, THRESHOLDS, "AAGCCCG")

        gcc_occ[tp]  = gcc_mean
        aag_occ[tp]  = aag_mean

        if tp == "T1":
            t1_gcc_site_sample = gcc_ss
            t1_aag_site_sample = aag_ss

    # ── Protection zone AUC (T1 only) ─────────────────────────────────────
    print("\n[3/4] Computing protection zone AUC (T1)...", flush=True)

    if not FEATURES_TSV.exists():
        print(f"  WARNING: features file not found: {FEATURES_TSV}", flush=True)
        print("  Skipping AUC computation.", flush=True)
        auc_by_threshold = {t: float("nan") for t in THRESHOLDS}
    else:
        genes_df = pd.read_csv(FEATURES_TSV, sep="\t")
        # Verify required columns
        required = {"tss", "is_exposed"}
        if not required.issubset(genes_df.columns):
            print(f"  WARNING: missing columns {required - set(genes_df.columns)}", flush=True)
            auc_by_threshold = {t: float("nan") for t in THRESHOLDS}
        else:
            # Combine GCCGGC (4mC) and AAGCCCG (6mA) site_sample for T1
            combined_site_sample = {}
            for t in THRESHOLDS:
                merged = defaultdict(list)
                for pos, freqs in t1_gcc_site_sample[t].items():
                    merged[pos].extend(freqs)
                for pos, freqs in t1_aag_site_sample[t].items():
                    merged[pos].extend(freqs)
                combined_site_sample[t] = dict(merged)

            auc_by_threshold = {}
            for t in THRESHOLDS:
                print(f"  AUC at threshold {t}...", end=" ", flush=True)
                auc_val = compute_protection_auc(
                    combined_site_sample, t, genes_df)
                auc_by_threshold[t] = auc_val
                print(f"{auc_val:.4f}", flush=True)

    # ── OR for AAGCCCG dual-modification (T1) ─────────────────────────────
    print("\n[4/4] Computing OR (AAGCCCG dual-modification, T1)...", flush=True)
    T1_FILES = [MODONLY_DIR / f"{s}_modonly.tsv" for s in SAMPLES["T1"]]

    or_by_threshold = {}
    comod_by_threshold = {}
    for t in THRESHOLDS:
        print(f"  OR at threshold {t}...", end=" ", flush=True)
        results = []
        for fp in T1_FILES:
            if fp.exists():
                results.append(analyze_comod(fp, aagcccg_map, min_prob=t))
        if results:
            pooled = pool(results)
            or_val    = compute_or(pooled)
            comod_val = comod_rate(pooled)
            or_by_threshold[t]    = or_val
            comod_by_threshold[t] = comod_val
            print(f"OR={or_val:.1f}, comod_rate={comod_val*100:.1f}%", flush=True)
        else:
            or_by_threshold[t]    = float("nan")
            comod_by_threshold[t] = float("nan")
            print("no data", flush=True)

    # ── Assemble output table ──────────────────────────────────────────────
    print("\n" + "=" * 70, flush=True)
    print("RESULTS TABLE", flush=True)
    print("=" * 70, flush=True)

    rows = []
    header = [
        "threshold",
        "GCCGGC_T1", "GCCGGC_T2", "GCCGGC_T3",
        "AAGCCCG_T1", "AAGCCCG_T2", "AAGCCCG_T3",
        "AUC_T1",
        "comod_rate_T1",
        "OR_value",
    ]
    print("\t".join(header))
    print("-" * 90)

    for t in THRESHOLDS:
        row = {
            "threshold":      t,
            "GCCGGC_T1":      round(gcc_occ["T1"][t], 2),
            "GCCGGC_T2":      round(gcc_occ["T2"][t], 2),
            "GCCGGC_T3":      round(gcc_occ["T3"][t], 2),
            "AAGCCCG_T1":     round(aag_occ["T1"][t], 2),
            "AAGCCCG_T2":     round(aag_occ["T2"][t], 2),
            "AAGCCCG_T3":     round(aag_occ["T3"][t], 2),
            "AUC_T1":         round(auc_by_threshold[t], 4)
                               if not np.isnan(auc_by_threshold[t]) else "NA",
            "comod_rate_T1":  round(comod_by_threshold[t] * 100, 2)
                               if not np.isnan(comod_by_threshold[t]) else "NA",
            "OR_value":       round(or_by_threshold[t], 1)
                               if not np.isnan(or_by_threshold[t]) else "NA",
        }
        rows.append(row)
        line = "\t".join(str(row[h]) for h in header)
        print(line)

    # ── Save to file ───────────────────────────────────────────────────────
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows, columns=header)
    df.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"\nSaved: {OUT_TSV}", flush=True)

    # ── Summary ────────────────────────────────────────────────────────────
    print("\n" + "=" * 70, flush=True)
    print("SUMMARY (for reviewer response)", flush=True)
    print("=" * 70, flush=True)

    for metric, vals in [
        ("GCCGGC_T1 occupancy (%)",  [gcc_occ["T1"][t] for t in THRESHOLDS]),
        ("AAGCCCG_T1 occupancy (%)", [aag_occ["T1"][t] for t in THRESHOLDS]),
        ("AUC (T1)",                 [auc_by_threshold[t] for t in THRESHOLDS]),
        ("OR (AAGCCCG comod)",       [or_by_threshold[t] for t in THRESHOLDS]),
    ]:
        vals_clean = [v for v in vals if not (isinstance(v, float) and np.isnan(v))]
        if vals_clean:
            range_str = f"[{min(vals_clean):.2f} – {max(vals_clean):.2f}]"
        else:
            range_str = "[N/A]"
        print(f"  {metric:<35} range {range_str}")


if __name__ == "__main__":
    main()
