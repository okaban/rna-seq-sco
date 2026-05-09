#!/usr/bin/env python3
"""
Figure 2 — AAGCCCG 4mC × 6mA co-modification (Story B, main text Fig 2A/B).

Reproduces the dual-modification narrative used in
`23_expanded_motif_search/02_dual_modification_and_novel_motif_analysis.py`
and the contingency framework of
`68_or_permutation/A3_OR_CI_permutation.py`, using the validated call
set in `07_motif_analysis/methylation_site_sequences.csv` (modkit pileup
≥ 50 % modified, per-timepoint).

Panel A — Pileup co-modification at AAGCCCG vs hypergeometric null
    Histogram of permutation null OR distribution (n=10,000) with the
    observed OR + 95 % CI.  Pool of 4mC×6mA AAGCCCG sites de-duplicated
    by (chrom, position, strand, mod_type) across all timepoints to
    match the 68_or_permutation A-3 frame.

Panel B — Chromosomal distribution of co-modified AAGCCCG positions
    Each timepoint (T1, T2, T3) is a horizontal track plotting the
    genomic positions of (4mC, 6mA) co-localized pairs (≤ 10 bp,
    same strand).  Each pair is plotted at the midpoint of the two
    positions.

Co-modification call (per pair, per timepoint):
    A "co-mod pair at Tk" = a 4mC call at an AAGCCCG site (sequence
    window center inside AAGCCCG) AND a 6mA call at an AAGCCCG site,
    same strand, |Δposition| ≤ 10 bp, both called in timepoint Tk.

Outputs:
    figures/Figure2_perread_comod.{png,pdf}
    Writing/fig_images/Figure2_perread_comod.{png,pdf}
    tables/comod_pairs_by_timepoint.tsv
    tables/figure2_summary.tsv
"""

from __future__ import annotations

import re
from collections import Counter
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd
from scipy import stats

# ── Paths ────────────────────────────────────────────────────────────────────────
REF_FA = "/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa"
SITES_CSV = Path(
    "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/"
    "07_motif_analysis/methylation_site_sequences.csv"
)
OUT_DIR = Path(
    "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/65_per_read_comod"
)
FIG_DIR = OUT_DIR / "figures"
TAB_DIR = OUT_DIR / "tables"
FIG_DIR.mkdir(parents=True, exist_ok=True)
TAB_DIR.mkdir(parents=True, exist_ok=True)

WRITING_FIG_DIR = Path("/Users/okaban/bioinfo/rna-seq/Writing/fig_images")
WRITING_FIG_DIR.mkdir(parents=True, exist_ok=True)

TIMEPOINTS = ("T1", "T2", "T3")
TP_COLOR = {"T1": "#2C5AA0", "T2": "#E67E22", "T3": "#C0392B"}

# ── Calling parameters (must match 23_expanded_motif_search & 68_or_permutation) ─
MOTIF = "AAGCCCG"
CENTER_POS = 15
SEQ_LEN = 31
COLOC_DISTANCE = 10
N_UNIVERSE = 1_238_215  # genome-wide 4mC×6mA candidate pair universe
                         # (from 68_or_permutation/A3 — empirical pileup pair sweep)
N_AAG = 1_334            # AAGCCCG motif instances + and - strand combined
N_PERM = 10_000


# ── Reference & motif placement ──────────────────────────────────────────────────
def load_reference(fa_path: str) -> str:
    parts: list[str] = []
    with open(fa_path) as f:
        for line in f:
            if line.startswith(">"):
                if parts:
                    break
            else:
                parts.append(line.strip().upper())
    return "".join(parts)


def reverse_complement(seq: str) -> str:
    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))


def motif_covers_center(seq: str, motif: str = MOTIF) -> bool:
    """Return True if the methylated base (sequence center) lies inside `motif`
    on either strand of the 31 bp window (matches 23_expanded_motif_search)."""
    for m in re.finditer(motif, seq):
        if m.start() <= CENTER_POS < m.end():
            return True
    rc = reverse_complement(seq)
    rc_center = SEQ_LEN - 1 - CENTER_POS
    for m in re.finditer(motif, rc):
        if m.start() <= rc_center < m.end():
            return True
    return False


# ── Statistics ──────────────────────────────────────────────────────────────────
def or_with_ci(a, b, c, d):
    """Haldane-corrected OR + 95 % CI (log-OR normal) + Fisher's exact (greater)."""
    cells = np.array([a, b, c, d], dtype=float)
    if (cells == 0).any():
        cells = cells + 0.5
    aa, bb, cc, dd = cells
    or_val = (aa * dd) / (bb * cc)
    log_or = np.log(or_val)
    se = float(np.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd))
    ci_lo = float(np.exp(log_or - 1.96 * se))
    ci_hi = float(np.exp(log_or + 1.96 * se))
    _, p_fisher = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
    return float(or_val), ci_lo, ci_hi, float(p_fisher)


def hypergeometric_null(n_total: int, n_aag: int, n_comod: int,
                         n_perm: int = N_PERM, seed: int = 42):
    """OR null under H0: AAGCCCG label is independent of co-mod status."""
    rng = np.random.default_rng(seed)
    a_perm = rng.hypergeometric(
        ngood=n_comod,
        nbad=max(n_total - n_comod, 0),
        nsample=n_aag,
        size=n_perm,
    )
    b_perm = n_aag - a_perm
    c_perm = n_comod - a_perm
    d_perm = (n_total - n_aag) - c_perm

    cells = np.stack([a_perm, b_perm, c_perm, d_perm], axis=1).astype(float)
    zero_mask = (cells == 0).any(axis=1)
    cells[zero_mask] += 0.5
    or_vals = (cells[:, 0] * cells[:, 3]) / (cells[:, 1] * cells[:, 2])
    return or_vals, a_perm


# ── Co-localization ─────────────────────────────────────────────────────────────
def find_pairs(df_4mC: pd.DataFrame, df_6mA: pd.DataFrame,
                strand_aware: bool = True,
                distance: int = COLOC_DISTANCE) -> pd.DataFrame:
    """
    Return DataFrame of (pos_4mC, pos_6mA, distance, strand_4mC, strand_6mA, midpoint).
    Quadratic but fast enough for ≤ 1k × 1k.
    """
    rows: list = []
    arr4 = df_4mC[["position", "strand"]].to_numpy()
    arr6 = df_6mA[["position", "strand"]].to_numpy()
    for p4, s4 in arr4:
        for p6, s6 in arr6:
            if strand_aware and s4 != s6:
                continue
            d = abs(int(p4) - int(p6))
            if d <= distance:
                rows.append({
                    "pos_4mC": int(p4),
                    "pos_6mA": int(p6),
                    "strand": s4,
                    "distance": d,
                    "midpoint": (int(p4) + int(p6)) // 2,
                })
    return pd.DataFrame(rows, columns=[
        "pos_4mC", "pos_6mA", "strand", "distance", "midpoint",
    ])


# ── Main ────────────────────────────────────────────────────────────────────────
def main() -> None:
    print("=" * 72)
    print("Figure 2 — AAGCCCG 4mC × 6mA co-modification (Story B)")
    print("=" * 72)

    print("\n[1/5] Loading reference & methylation calls...")
    genome = load_reference(REF_FA)
    df = pd.read_csv(SITES_CSV)
    print(f"  rows: {len(df):,}")
    print(df.groupby(["mod_type", "timepoint"]).size().to_string())

    print("\n[2/5] Annotating each call: at AAGCCCG?")
    df["at_aag"] = df["sequence"].apply(motif_covers_center)
    aag_count = df.groupby(["mod_type", "timepoint"])["at_aag"].sum()
    print(aag_count.to_string())

    print("\n[3/5] Finding co-localized (4mC, 6mA) pairs at AAGCCCG, by timepoint...")
    pair_dfs: dict = {}
    for tp in TIMEPOINTS:
        d4 = df[(df["mod_type"] == "4mC") & (df["timepoint"] == tp) & (df["at_aag"])]
        d6 = df[(df["mod_type"] == "6mA") & (df["timepoint"] == tp) & (df["at_aag"])]
        pdf = find_pairs(d4, d6, strand_aware=True)
        pair_dfs[tp] = pdf
        dist_break = (
            sorted(Counter(pdf["distance"]).items()) if len(pdf) else []
        )
        print(f"  {tp}: 4mC={len(d4):>4d}  6mA={len(d6):>4d}  "
              f"co-mod pairs={len(pdf):>4d}  distances={dist_break}")

    # Pool: dedup by (chrom, position, strand, mod_type) across all timepoints,
    # matching the 23_expanded_motif_search/02 script exactly so the pool count
    # reproduces the 244-pair (≈) number reported in the figure plan.
    df_unique = df.drop_duplicates(
        subset=["chrom", "position", "strand", "mod_type"]
    )
    d4_pool = df_unique[(df_unique["mod_type"] == "4mC") & (df_unique["at_aag"])]
    d6_pool = df_unique[(df_unique["mod_type"] == "6mA") & (df_unique["at_aag"])]
    pdf_pool = find_pairs(d4_pool, d6_pool, strand_aware=True)
    n_comod_pool = len(pdf_pool)
    dist_break_pool = sorted(Counter(pdf_pool["distance"]).items())
    print(f"\n  POOL union (dedup by pos+strand+mod_type, ANY timepoint): "
          f"4mC={len(d4_pool)}  6mA={len(d6_pool)}  "
          f"co-mod pairs={n_comod_pool}  dist={dist_break_pool}")

    # Save the complete pair list
    rows_out: list = []
    for tp in TIMEPOINTS:
        for _, r in pair_dfs[tp].iterrows():
            rows_out.append({
                "timepoint": tp,
                "pos_4mC": r["pos_4mC"],
                "pos_6mA": r["pos_6mA"],
                "strand": r["strand"],
                "distance": r["distance"],
                "midpoint": r["midpoint"],
            })
    if rows_out:
        pd.DataFrame(rows_out).sort_values(
            ["timepoint", "midpoint"]
        ).to_csv(TAB_DIR / "comod_pairs_by_timepoint.tsv", sep="\t", index=False)
    else:
        pd.DataFrame(columns=["timepoint", "pos_4mC", "pos_6mA",
                              "strand", "distance", "midpoint"]).to_csv(
            TAB_DIR / "comod_pairs_by_timepoint.tsv", sep="\t", index=False)
    print(f"  wrote {TAB_DIR/'comod_pairs_by_timepoint.tsv'} "
          f"({len(rows_out)} rows)")

    # ── [4/5] Statistics: 2×2 contingency on POOL ──────────────────────────────
    print("\n[4/5] Statistics — pooled contingency (matching 68_or_permutation/A3)")
    a_obs = n_comod_pool
    b_obs = N_AAG - a_obs
    c_obs = 2  # genome-wide out-AAGCCCG co-mod pair count from 68_or_permutation/A3
    d_obs = N_UNIVERSE - N_AAG - c_obs
    or_val, ci_lo, ci_hi, p_fisher = or_with_ci(a_obs, b_obs, c_obs, d_obs)
    print(f"  Contingency: in={a_obs}  in_not={b_obs}  "
          f"out={c_obs}  out_not={d_obs:,}")
    print(f"  OR = {or_val:,.0f}    95% CI = [{ci_lo:,.0f}, {ci_hi:,.3e}]")
    print(f"  Fisher's exact p (one-sided) = {p_fisher:.3e}")

    expected_in_aag = (a_obs + c_obs) * N_AAG / N_UNIVERSE
    fold_enrich = a_obs / max(expected_in_aag, 1e-9)
    print(f"  Expected in-AAGCCCG co-mod under H0 = {expected_in_aag:.3f}")
    print(f"  Observed/expected fold-enrichment = {fold_enrich:.0f}×")

    print(f"\n  Hypergeometric permutation null (n={N_PERM:,})...")
    perm_ors, perm_a = hypergeometric_null(
        N_UNIVERSE, N_AAG, a_obs + c_obs, n_perm=N_PERM, seed=42,
    )
    perm_p = float(np.mean(perm_ors >= or_val))
    perm_p_str = f"< {1/N_PERM:.4f}" if perm_p == 0.0 else f"{perm_p:.4f}"
    print(f"  permutation p (OR ≥ observed) = {perm_p_str}")
    print(f"  null OR median={np.median(perm_ors):.3f}  max={perm_ors.max():.3f}")

    # Save summary
    summary = {
        "n_aagcccg_motifs": N_AAG,
        "n_universe_pairs": N_UNIVERSE,
        "co_mod_pairs_T1": len(pair_dfs["T1"]),
        "co_mod_pairs_T2": len(pair_dfs["T2"]),
        "co_mod_pairs_T3": len(pair_dfs["T3"]),
        "co_mod_pairs_pool": n_comod_pool,
        "expected_in_aag_under_H0": float(f"{expected_in_aag:.6f}"),
        "fold_enrichment": float(f"{fold_enrich:.4f}"),
        "OR": float(f"{or_val:.4f}"),
        "OR_CI_lower": float(f"{ci_lo:.4f}"),
        "OR_CI_upper": float(f"{ci_hi:.4f}"),
        "fisher_p_one_sided": p_fisher,
        "permutation_p": perm_p_str,
        "permutation_n": N_PERM,
        "perm_null_median_OR": float(f"{np.median(perm_ors):.4f}"),
        "perm_null_max_OR": float(f"{perm_ors.max():.4f}"),
    }
    pd.Series(summary).to_csv(TAB_DIR / "figure2_summary.tsv",
                              sep="\t", header=False)
    print(f"  wrote {TAB_DIR/'figure2_summary.tsv'}")

    # ── [5/5] Figure ────────────────────────────────────────────────────────────
    print("\n[5/5] Building Figure 2...")
    fig = plt.figure(figsize=(11.6, 5.0))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.55], wspace=0.30,
                          left=0.06, right=0.96, top=0.84, bottom=0.16)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])

    # ── Panel A: permutation null vs observed ────────────────────────────────────
    log10_perm = np.log10(np.maximum(perm_ors, 1e-9))
    log10_obs = np.log10(or_val)
    log10_ci_lo = np.log10(max(ci_lo, 1e-9))
    log10_ci_hi = np.log10(max(ci_hi, 1e-9))

    axA.hist(log10_perm, bins=60, color="#8FB7E5", edgecolor="white",
             linewidth=0.4, alpha=0.95,
             label=f"Null distribution\n(n = {N_PERM:,} permutations)")
    axA.axvline(log10_obs, color="#C0392B", linewidth=2.4,
                label=f"Observed OR = {or_val:,.0f}")
    axA.axvspan(log10_ci_lo, log10_ci_hi, alpha=0.15, color="#C0392B",
                label=f"95% CI [{ci_lo:,.0f}, {ci_hi:,.0f}]")

    fisher_p_text = (
        "Fisher p < 1e-300" if p_fisher < 1e-300
        else f"Fisher p = {p_fisher:.1e}"
    )
    axA.text(
        0.97, 0.96,
        f"Observed = {a_obs} pairs\n"
        f"Expected = {expected_in_aag:.2f}\n"
        f"Fold-enrichment = {fold_enrich:.0f}×\n\n"
        f"OR = {or_val:,.0f}\n"
        f"95% CI [{ci_lo:,.0f}, {ci_hi:,.1e}]\n"
        f"perm p {perm_p_str}\n"
        f"{fisher_p_text}",
        transform=axA.transAxes, fontsize=8.5,
        va="top", ha="right",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
                  edgecolor="#cccccc", alpha=0.92),
    )
    axA.set_xlabel("log₁₀(Odds Ratio)", fontsize=10)
    axA.set_ylabel(f"Permutations (of {N_PERM:,})", fontsize=10)
    axA.set_title("A. AAGCCCG 4mC × 6mA co-modification\n"
                  "vs hypergeometric null (T1∪T2∪T3 union)",
                  fontsize=10.5, loc="left")
    axA.legend(fontsize=8, loc="upper left")
    axA.grid(axis="y", linewidth=0.4, alpha=0.4)

    # ── Panel B: chromosomal distribution by timepoint ───────────────────────────
    chrom_len = len(genome)
    tracks_y = {"T1": 2.0, "T2": 1.0, "T3": 0.0}
    track_h = 0.55

    for tp in TIMEPOINTS:
        y = tracks_y[tp]
        axB.add_patch(Rectangle((0, y - track_h / 2), chrom_len, track_h,
                                facecolor="#f3f3f3", edgecolor="#bbbbbb",
                                linewidth=0.5, zorder=1))
        positions = pair_dfs[tp]["midpoint"].sort_values().values
        n_pairs = len(pair_dfs[tp])
        if n_pairs:
            axB.scatter(positions, [y] * len(positions),
                        marker="|", s=110,
                        color=TP_COLOR[tp], alpha=0.65,
                        linewidth=1.1, zorder=3)
        axB.text(chrom_len * 1.005, y, f"n = {n_pairs} pairs",
                 va="center", ha="left", fontsize=9, color=TP_COLOR[tp])

    ori_pos = 4_270_000
    axB.axvline(ori_pos, color="#444", linestyle=":", linewidth=0.9, zorder=2,
                alpha=0.6)
    axB.text(ori_pos, 2.85, "ori (~4.27 Mb)", fontsize=8,
             color="#444", ha="center")

    axB.set_yticks([tracks_y[k] for k in TIMEPOINTS])
    axB.set_yticklabels(list(TIMEPOINTS), fontsize=10)
    axB.set_ylim(-0.7, 2.95)
    axB.set_xlim(-50_000, chrom_len * 1.10)
    xticks = np.arange(0, chrom_len + 1, 1_000_000)
    axB.set_xticks(xticks)
    axB.set_xticklabels([f"{int(x / 1e6)}" for x in xticks], fontsize=9)
    axB.set_xlabel("Chromosomal position (Mb, NC_003888.3)", fontsize=10)
    axB.set_title(
        f"B. Chromosomal distribution of co-modified AAGCCCG pairs by timepoint "
        f"(pool = {n_comod_pool})",
        fontsize=10.5, loc="left",
    )
    for sp in ("top", "right"):
        axB.spines[sp].set_visible(False)
    axB.tick_params(left=False)
    axB.grid(axis="x", linewidth=0.3, alpha=0.4)

    fig.suptitle(
        "Figure 2 — Per-read co-modification at AAGCCCG (SC_RS17645 dual targeting)",
        fontsize=11.8, y=0.985, x=0.06, ha="left", fontweight="bold",
    )

    out_png = FIG_DIR / "Figure2_perread_comod.png"
    out_pdf = FIG_DIR / "Figure2_perread_comod.pdf"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    print(f"  saved {out_png}")
    print(f"  saved {out_pdf}")

    publish_png = WRITING_FIG_DIR / "Figure2_perread_comod.png"
    publish_pdf = WRITING_FIG_DIR / "Figure2_perread_comod.pdf"
    fig.savefig(publish_png, dpi=300, bbox_inches="tight")
    fig.savefig(publish_pdf, bbox_inches="tight")
    print(f"  saved {publish_png}")
    print(f"  saved {publish_pdf}")

    plt.close(fig)
    print("\nDone.")


if __name__ == "__main__":
    main()
