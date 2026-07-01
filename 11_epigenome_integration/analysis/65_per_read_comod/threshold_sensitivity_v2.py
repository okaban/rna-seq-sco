#!/usr/bin/env python3
"""
Threshold sensitivity analysis v2 (corrected).

Fixes from v1:
  - Occupancy = mean freq of DETECTED sites (freq>=MIN_FREQ), not all canonical positions
  - OR = uses per-read counts with fixed background (c=2, d=1,236,879 from paper)
  - AUC = uses ALL high_confidence_sites positions (not only GCCGGC+AAGCCCG)

Output: Supplementary_threshold_sensitivity.tsv
"""

import re, sys, warnings
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
REF_FA       = "/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa"
MODONLY_DIR  = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/65_per_read_comod")
PILEUP_DIR   = Path("/Users/okaban/bioinfo/methyl/260102_M145/analysis/cursor_results/pileup")
FEATURES_TSV = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/52_shielded_exposed_boundary/tables/all_genes_features.tsv")
HIGHCONF_CSV = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/01_integration/high_confidence_sites_weighted.csv")
OUT_TSV      = Path("/Users/okaban/obsidian/Research/rna-seq/Writing/Supplementary_threshold_sensitivity.tsv")

# ── Parameters ─────────────────────────────────────────────────────────────
THRESHOLDS = [0.5, 0.6, 0.75, 0.8, 0.9]

SAMPLES = {
    "T1": ["1-1", "1-2", "1-3"],
    "T2": ["2-1", "2-3", "2-4"],
    "T3": ["3-2", "3-3", "3-4"],
}

CALL_6mA = "a"
CALL_4mC = "21839"

# GCCGGC: 4mC at positions 1,2,5 (0-indexed within GCCGGC)
GCCGGC_4mC_POSITIONS = {1, 2, 5}
# AAGCCCG: 6mA at positions 0,1 (0-indexed within AAGCCCG)
AAGCCCG_6mA_POSITIONS = {0, 1}

OCC_MIN_FREQ   = 20.0   # % – site must have mean freq >= this to count as "detected"
OCC_MIN_REPS   = 1      # min replicates with coverage

AUC_MIN_FREQ   = 20.0   # % – same filter for AUC site inclusion
AUC_MIN_REPS   = 2      # min replicates for AUC site
AUC_WINDOW_BP  = 2000
SEED           = 42
N_SPLITS       = 5

# OR background (from paper S-4 pileup analysis at threshold 0.75)
OR_BACKGROUND_C = 2          # co-mod reads outside AAGCCCG
OR_BACKGROUND_D = 1_236_879  # non-co-mod reads outside AAGCCCG


# ─────────────────────────────────────────────────────────────────────────────
# Reference & motif positions
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

def build_canonical_positions(genome: str, motif: str,
                              mod_positions: set) -> dict:
    """
    Returns {ref_pos: strand} where strand is '+' or '-'.
    For palindromic motifs (e.g. GCCGGC), + and - strand scan the SAME
    reference sequence but produce different ref positions.
    """
    positions = {}  # {pos: strand}
    motif_len = len(motif)
    rc_motif  = rc(motif)
    for m in re.finditer(f"(?={motif})", genome):
        for p in mod_positions:
            pos = m.start() + p
            positions[pos] = "+"
    for m in re.finditer(f"(?={rc_motif})", genome):
        for p in mod_positions:
            pos = m.start() + (motif_len - 1 - p)
            if pos not in positions:  # + strand takes priority for palindromic overlaps
                positions[pos] = "-"
    return positions

def build_aagcccg_comod_map(genome: str) -> dict:
    """Position map for per-read co-modification analysis."""
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
        for p in target: reg(m.start(), p, "+")
    for m in re.finditer(f"(?={rc_motif})", genome):
        for p in target: reg(m.start(), p, "-")
    return pos_map


# ─────────────────────────────────────────────────────────────────────────────
# Pileup coverage
# ─────────────────────────────────────────────────────────────────────────────

def load_pileup_coverage(pileup_path: Path,
                         target_pos_strand: dict,
                         call_code: str) -> dict:
    """
    Return {pos: total_cov} (strand-specific), where total_cov = N_valid_cov + N_fail.
    N_valid_cov = N_mod + N_canonical (reads that received a definitive call at the
    pileup's auto-threshold).  N_fail = ambiguous reads (probability between low and
    high thresholds).  Together they equal all reads that received any methylation
    assessment, which is the correct denominator when counting reads at an arbitrary
    threshold T from the modonly file.  Using N_valid_cov alone would be too small
    when T < auto_threshold (~0.9375) because many reads in N_fail would then qualify
    as modified, making counts > denominator (occupancy > 100%).
    target_pos_strand: {pos: strand} where strand is '+' or '-'.
    """
    cov = {}
    with open(pileup_path) as f:
        for line in f:
            parts = line.split("\t")
            if len(parts) < 16:
                continue
            if parts[3] != call_code:
                continue
            pos    = int(parts[1])
            strand = parts[5].strip()
            if pos not in target_pos_strand:
                continue
            if strand != target_pos_strand[pos]:
                continue
            n_valid = int(parts[9])
            n_fail  = int(parts[15])
            total   = n_valid + n_fail
            if total > 0:
                cov[pos] = total
    return cov


# ─────────────────────────────────────────────────────────────────────────────
# Per-read modification counts
# ─────────────────────────────────────────────────────────────────────────────

def count_mod_per_threshold(modonly_path: Path,
                             target_positions,   # set or dict keys
                             call_code: str,
                             thresholds: list) -> dict:
    """
    Return {threshold: {pos: count}}. Single-pass over file.
    No strand filtering needed: modonly calls at position X for a given
    call_code are intrinsically strand-specific (only the strand with the
    canonical base makes modification calls there).
    """
    counts  = {t: defaultdict(int) for t in thresholds}
    min_t   = min(thresholds)
    with open(modonly_path) as f:
        f.readline()
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
# Occupancy: mean freq of DETECTED sites
# ─────────────────────────────────────────────────────────────────────────────

def compute_occupancy_timepoint(samples_list: list,
                                 target_pos_strand: dict,  # {pos: strand}
                                 call_code: str,
                                 thresholds: list,
                                 label: str) -> tuple:
    """
    Returns:
      occ_mean:   {threshold: mean_freq}   (mean of detected sites)
      occ_nsites: {threshold: n_sites}     (number of detected sites)
      site_sample: {threshold: {pos: [per-sample freqs]}}   (for AUC reuse)
    """
    # site_sample[threshold][pos] = [freq_sample1, freq_sample2, ...]
    site_sample = {t: defaultdict(list) for t in thresholds}

    for sample in samples_list:
        pileup_path  = PILEUP_DIR / f"{sample}_pileup.bed"
        modonly_path = MODONLY_DIR / f"{sample}_modonly.tsv"
        if not pileup_path.exists() or not modonly_path.exists():
            print(f"    WARN: missing files for {sample}", flush=True)
            continue

        cov    = load_pileup_coverage(pileup_path, target_pos_strand, call_code)
        counts = count_mod_per_threshold(modonly_path, target_pos_strand,
                                          call_code, thresholds)

        for pos, total_cov in cov.items():
            if total_cov == 0:
                continue
            for t in thresholds:
                freq = counts[t].get(pos, 0) / total_cov * 100.0
                site_sample[t][pos].append(freq)

    occ_mean   = {}
    occ_nsites = {}
    for t in thresholds:
        detected_freqs = []
        for pos, freqs in site_sample[t].items():
            if len(freqs) >= OCC_MIN_REPS:
                mean_f = np.mean(freqs)
                if mean_f >= OCC_MIN_FREQ:
                    detected_freqs.append(mean_f)
        occ_mean[t]   = np.mean(detected_freqs) if detected_freqs else float("nan")
        occ_nsites[t] = len(detected_freqs)

    print(f"    [{label}] detected sites at each threshold: "
          + ", ".join(f"{t}→{occ_nsites[t]}" for t in thresholds), flush=True)
    return occ_mean, occ_nsites, site_sample


# ─────────────────────────────────────────────────────────────────────────────
# AUC: protection zone classifier using all high_confidence sites
# ─────────────────────────────────────────────────────────────────────────────

def load_highconf_positions(highconf_csv: Path, timepoint: str = "T1") -> dict:
    """
    Returns {pos: (call_code, strand)} for all sites in the given timepoint.
    call_code: 'a' for 6mA, '21839' for 4mC.
    strand: '+' or '-'.
    """
    df = pd.read_csv(highconf_csv)
    tp_df = df[df["timepoint"] == timepoint]
    pos_info = {}
    for _, row in tp_df.iterrows():
        code   = CALL_4mC if row["mod_type"] == "4mC" else CALL_6mA
        strand = str(row["strand"]).strip() if "strand" in row else "+"
        pos_info[int(row["position"])] = (code, strand)
    return pos_info


def compute_auc_from_sites(active_positions: np.ndarray,
                            genes_df: pd.DataFrame) -> float:
    """Run 5-fold stratified CV logistic regression using -nearest_methyl_dist."""
    def nearest_dist(tss, positions):
        if len(positions) == 0:
            return float("inf")
        return float(np.min(np.abs(positions - tss)))

    valid = genes_df[genes_df["is_exposed"].notna()].copy()
    valid["nd"] = valid["tss"].apply(
        lambda t: nearest_dist(int(t), active_positions))

    # Keep finite distances
    valid = valid[valid["nd"].apply(np.isfinite)]

    if valid["is_exposed"].sum() < 3 or len(valid) < 10:
        return float("nan")

    X = (-valid["nd"].values).reshape(-1, 1)
    y = valid["is_exposed"].astype(int).values

    skf = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=SEED)
    fold_aucs = []
    for tr, te in skf.split(X, y):
        if len(np.unique(y[te])) < 2:
            continue
        lr = LogisticRegression(max_iter=1000, random_state=SEED)
        lr.fit(X[tr], y[tr])
        fold_aucs.append(roc_auc_score(y[te], lr.predict_proba(X[te])[:, 1]))

    return float(np.mean(fold_aucs)) if fold_aucs else float("nan")


def compute_all_auc(all_pos_code: dict,   # {pos: code}
                    pileup_cov_all: dict,  # {(pos, code): total_cov} for T1
                    t1_counts_all: dict,   # {threshold: {(pos,code): n_mod}}
                    genes_df: pd.DataFrame,
                    thresholds: list) -> dict:
    """
    For each threshold T, build the set of 'active' T1 methylation positions:
    - freq(T) = n_mod(T) / total_cov >= AUC_MIN_FREQ  in >= AUC_MIN_REPS samples
    Then compute AUC.
    """
    auc_vals = {}
    for t in thresholds:
        print(f"  AUC at T={t}...", end=" ", flush=True)
        active = []
        for pos, code in all_pos_code.items():
            key = (pos, code)
            if key not in pileup_cov_all:
                continue
            total_cov = pileup_cov_all[key]
            if total_cov == 0:
                continue
            n_mod = t1_counts_all[t].get(key, 0)
            freq  = n_mod / total_cov * 100.0
            if freq >= AUC_MIN_FREQ:
                active.append(pos)
        active_arr = np.array(sorted(set(active)))
        n_act = len(active_arr)
        auc_v = compute_auc_from_sites(active_arr, genes_df)
        auc_vals[t] = auc_v
        print(f"{n_act} sites, AUC={auc_v:.4f}", flush=True)
    return auc_vals


# ─────────────────────────────────────────────────────────────────────────────
# OR: per-read co-modification with fixed background
# ─────────────────────────────────────────────────────────────────────────────

def analyze_comod(modonly_path: Path, pos_map: dict,
                  min_prob: float) -> dict:
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
                pos  = int(parts[1])
                prob = float(parts[2])
            except ValueError:
                continue
            if prob < min_prob:
                continue
            if pos not in pos_map:
                continue
            for (ms, pip, strand, em) in pos_map[pos]:
                if em == "6mA" and parts[3] == CALL_6mA:
                    read_motif[parts[0]][ms].add("6mA")
                elif em == "4mC" and parts[3] == CALL_4mC:
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


def compute_or_from_perread(r: dict) -> float:
    """
    OR = (a * d) / (b * c)
    a = dual-modified (AAGCCCG) = pairs_both
    b = not dual-modified (AAGCCCG) = pairs_4mC - pairs_both
    c = background co-mod (genome) = OR_BACKGROUND_C = 2
    d = background no co-mod = OR_BACKGROUND_D = 1,236,879
    """
    a = r["pairs_both"]
    b = r["pairs_4mC"] - a
    c = OR_BACKGROUND_C
    d = OR_BACKGROUND_D
    if b <= 0 or a <= 0:
        return float("nan")
    return float(a * d) / float(b * c)


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    print("=" * 70, flush=True)
    print("Threshold Sensitivity Analysis v2", flush=True)
    print("=" * 70, flush=True)

    # ── Reference genome ───────────────────────────────────────────────────
    print("\n[0] Loading reference genome...", flush=True)
    genome = load_reference(REF_FA)
    print(f"  {len(genome):,} bp", flush=True)

    gccggc_pos  = build_canonical_positions(genome, "GCCGGC",  GCCGGC_4mC_POSITIONS)
    aagcccg_pos = build_canonical_positions(genome, "AAGCCCG", AAGCCCG_6mA_POSITIONS)
    aagcccg_map = build_aagcccg_comod_map(genome)

    print(f"  GCCGGC canonical 4mC positions: {len(gccggc_pos):,}", flush=True)
    print(f"  AAGCCCG canonical 6mA positions: {len(aagcccg_pos):,}", flush=True)

    # ── Load high_confidence_sites for AUC ────────────────────────────────
    print("\n[1] Loading high-confidence sites for AUC...", flush=True)
    # hc_pos_info: {pos: (call_code, strand)}
    hc_pos_info = load_highconf_positions(HIGHCONF_CSV, timepoint="T1")
    all_hc_positions = set(hc_pos_info.keys())
    print(f"  T1 high-confidence positions: {len(all_hc_positions):,}", flush=True)

    # ── Occupancy per timepoint ────────────────────────────────────────────
    print("\n[2] Computing occupancy per timepoint...", flush=True)

    gcc_occ  = {}
    aag_occ  = {}
    gcc_nsites = {}
    aag_nsites = {}

    for tp, samples in SAMPLES.items():
        print(f"\n  Timepoint {tp}:", flush=True)
        g_mean, g_n, _ = compute_occupancy_timepoint(
            samples, gccggc_pos, CALL_4mC, THRESHOLDS, "GCCGGC-4mC")
        a_mean, a_n, _ = compute_occupancy_timepoint(
            samples, aagcccg_pos, CALL_6mA, THRESHOLDS, "AAGCCCG-6mA")
        gcc_occ[tp]    = g_mean
        aag_occ[tp]    = a_mean
        gcc_nsites[tp] = g_n
        aag_nsites[tp] = a_n

    # ── AUC: process T1 using ALL high-confidence positions ───────────────
    print("\n[3] Computing AUC (T1, all high-confidence sites)...", flush=True)

    # Build pileup coverage for all high-confidence positions
    # {(pos, code): total_cov}
    t1_samples = SAMPLES["T1"]
    pileup_cov_all = {}
    t1_counts_all  = {t: defaultdict(int) for t in THRESHOLDS}

    for sample in t1_samples:
        pileup_path  = PILEUP_DIR / f"{sample}_pileup.bed"
        modonly_path = MODONLY_DIR / f"{sample}_modonly.tsv"
        if not pileup_path.exists() or not modonly_path.exists():
            continue

        print(f"  Loading T1 sample {sample}...", flush=True)

        # Pileup coverage for ALL high-confidence positions (strand-aware)
        with open(pileup_path) as f:
            for line in f:
                parts = line.split("\t")
                if len(parts) < 16:
                    continue
                code   = parts[3]
                pstrand = parts[5].strip()
                if code not in (CALL_4mC, CALL_6mA):
                    continue
                pos = int(parts[1])
                if pos not in all_hc_positions:
                    continue
                exp_code, exp_strand = hc_pos_info[pos]
                if code != exp_code:
                    continue
                if pstrand != exp_strand:
                    continue
                key = (pos, code)
                total_cov = int(parts[9]) + int(parts[15])
                if key not in pileup_cov_all:
                    pileup_cov_all[key] = []
                pileup_cov_all[key].append(total_cov)

        # Per-read counts (no strand filter needed for modonly)
        pos_4mC = {p for p, (c, s) in hc_pos_info.items() if c == CALL_4mC}
        pos_6mA = {p for p, (c, s) in hc_pos_info.items() if c == CALL_6mA}

        print(f"    Counting 4mC...", end=" ", flush=True)
        cnt4 = count_mod_per_threshold(modonly_path, pos_4mC, CALL_4mC, THRESHOLDS)
        print(f"    Counting 6mA...", end=" ", flush=True)
        cnt6 = count_mod_per_threshold(modonly_path, pos_6mA, CALL_6mA, THRESHOLDS)
        print("done", flush=True)

        for t in THRESHOLDS:
            for pos, n in cnt4[t].items():
                t1_counts_all[t][(pos, CALL_4mC)] += n
            for pos, n in cnt6[t].items():
                t1_counts_all[t][(pos, CALL_6mA)] += n

    # Average coverage across T1 samples
    pileup_cov_avg = {}
    for key, covs in pileup_cov_all.items():
        pileup_cov_avg[key] = np.mean(covs) if covs else 0

    # For compute_all_auc: convert hc_pos_info to {pos: code} for AUC
    hc_pos_code = {pos: code for pos, (code, strand) in hc_pos_info.items()}

    # Genes dataframe for AUC
    genes_df = pd.read_csv(FEATURES_TSV, sep="\t") if FEATURES_TSV.exists() else None
    if genes_df is None or "is_exposed" not in genes_df.columns:
        print("  WARNING: features file missing, skipping AUC", flush=True)
        auc_vals = {t: float("nan") for t in THRESHOLDS}
    else:
        auc_vals = compute_all_auc(
            hc_pos_code, pileup_cov_avg, t1_counts_all, genes_df, THRESHOLDS)

    # ── OR: per-read co-modification ──────────────────────────────────────
    print("\n[4] Computing OR (AAGCCCG dual-modification, T1)...", flush=True)
    T1_FILES = [MODONLY_DIR / f"{s}_modonly.tsv" for s in SAMPLES["T1"]]

    comod_rate_vals = {}
    or_vals         = {}

    for t in THRESHOLDS:
        print(f"  T={t}...", end=" ", flush=True)
        results = [analyze_comod(fp, aagcccg_map, min_prob=t)
                   for fp in T1_FILES if fp.exists()]
        if results:
            pooled  = pool(results)
            cr      = pooled["pairs_both"] / pooled["pairs_4mC"] if pooled["pairs_4mC"] else float("nan")
            or_val  = compute_or_from_perread(pooled)
            comod_rate_vals[t] = cr
            or_vals[t]         = or_val
            print(f"pairs_4mC={pooled['pairs_4mC']:,}, pairs_both={pooled['pairs_both']:,}, "
                  f"comod_rate={cr*100:.1f}%, OR={or_val:,.0f}", flush=True)
        else:
            comod_rate_vals[t] = float("nan")
            or_vals[t]         = float("nan")

    # ── Output table ───────────────────────────────────────────────────────
    print("\n" + "=" * 80, flush=True)
    print("RESULTS", flush=True)
    print("=" * 80, flush=True)

    header = [
        "threshold",
        "GCCGGC_T1", "GCCGGC_T2", "GCCGGC_T3",
        "GCCGGC_nsites_T1",
        "AAGCCCG_T1", "AAGCCCG_T2", "AAGCCCG_T3",
        "AAGCCCG_nsites_T1",
        "AUC_T1",
        "comod_rate_T1_pct",
        "OR_value",
    ]

    def fmt(x, dec=2):
        if isinstance(x, float) and np.isnan(x):
            return "NA"
        if isinstance(x, float):
            return f"{x:.{dec}f}"
        return str(x)

    print("\t".join(header))
    print("-" * 100)

    rows = []
    for t in THRESHOLDS:
        row = {
            "threshold":          t,
            "GCCGGC_T1":          fmt(gcc_occ["T1"][t]),
            "GCCGGC_T2":          fmt(gcc_occ["T2"][t]),
            "GCCGGC_T3":          fmt(gcc_occ["T3"][t]),
            "GCCGGC_nsites_T1":   gcc_nsites["T1"][t],
            "AAGCCCG_T1":         fmt(aag_occ["T1"][t]),
            "AAGCCCG_T2":         fmt(aag_occ["T2"][t]),
            "AAGCCCG_T3":         fmt(aag_occ["T3"][t]),
            "AAGCCCG_nsites_T1":  aag_nsites["T1"][t],
            "AUC_T1":             fmt(auc_vals[t], 4),
            "comod_rate_T1_pct":  fmt(comod_rate_vals[t] * 100 if not np.isnan(comod_rate_vals[t]) else float("nan")),
            "OR_value":           fmt(or_vals[t], 0),
        }
        rows.append(row)
        print("\t".join(str(row[h]) for h in header))

    # Save
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows, columns=header)
    df.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"\nSaved: {OUT_TSV}", flush=True)

    # ── Summary ────────────────────────────────────────────────────────────
    print("\n" + "=" * 70, flush=True)
    print("SUMMARY (for reviewer response, all 5 thresholds)", flush=True)
    print("=" * 70, flush=True)

    for label, vals in [
        ("GCCGGC_T1 occupancy (%)",  [gcc_occ["T1"][t] for t in THRESHOLDS]),
        ("AAGCCCG_T1 occupancy (%)", [aag_occ["T1"][t] for t in THRESHOLDS]),
        ("AUC T1 (5-fold CV)",       [auc_vals[t] for t in THRESHOLDS]),
        ("Co-mod rate T1 (%)",       [comod_rate_vals[t]*100 if not np.isnan(comod_rate_vals[t]) else float("nan") for t in THRESHOLDS]),
        ("OR (vs background c=2)",   [or_vals[t] for t in THRESHOLDS]),
    ]:
        clean = [v for v in vals if isinstance(v, float) and not np.isnan(v)]
        rng = f"[{min(clean):.2f} – {max(clean):.2f}]" if clean else "[N/A]"
        by_t = ", ".join(f"{t}→{fmt(v)}" for t, v in zip(THRESHOLDS, vals))
        print(f"  {label:<35} range {rng}")
        print(f"    {by_t}")


if __name__ == "__main__":
    main()
