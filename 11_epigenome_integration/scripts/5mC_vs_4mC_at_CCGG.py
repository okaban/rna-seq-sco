#!/usr/bin/env python3
"""
5mC vs 4mC comparison at CCGG/GGCCGG sites
Using existing modkit pileup BED files.

Dorado model: dna_r10.4.1_e8.2_400bps_sup@v5.2.0_4mC_5mC@v1
  → Both 4mC and 5mC are called. Code "m" = 5mC, code "21839" = 4mC, code "a" = 6mA

Purpose: Demonstrate that CCGG sites carry 4mC (not 5mC), contradicting
bisulfite-based reports (Pisciotta et al. 2023) that classified them as 5mC.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from Bio import SeqIO
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# === Paths ===
PILEUP_DIR = Path("/Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup")
REANALYSIS_DIR = Path("/Users/okaban/bioinfo/methyl/260102_M145/analysis/reanalysis_260224")
REF_FASTA = Path("/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa")
OUT_DIR = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/24_5mC_vs_4mC_CCGG")
OUT_DIR.mkdir(parents=True, exist_ok=True)
(OUT_DIR / "figures").mkdir(exist_ok=True)
(OUT_DIR / "tables").mkdir(exist_ok=True)

SAMPLES = {
    "1-1": ("T1", PILEUP_DIR / "1-1_pileup.bed"),
    "1-2": ("T1", PILEUP_DIR / "1-2_pileup.bed"),
    "1-3": ("T1", PILEUP_DIR / "1-3_pileup.bed"),
    "2-1": ("T2", PILEUP_DIR / "2-1_pileup.bed"),
    "2-3": ("T2", PILEUP_DIR / "2-3_pileup.bed"),
    "2-4": ("T2", PILEUP_DIR / "2-4_pileup.bed"),
    "3-2": ("T3", PILEUP_DIR / "3-2_pileup.bed"),
    "3-3": ("T3", PILEUP_DIR / "3-3_pileup.bed"),
    "3-4": ("T3", PILEUP_DIR / "3-4_pileup.bed"),
}

MIN_COVERAGE = 10

# === Load reference genome ===
print("Loading reference genome...")
ref_record = next(SeqIO.parse(REF_FASTA, "fasta"))
genome_seq = str(ref_record.seq).upper()
genome_len = len(genome_seq)
print(f"  Genome: {ref_record.id}, {genome_len:,} bp")

# === Find all CCGG positions in genome ===
print("\nFinding CCGG motif positions in genome...")
ccgg_positions_plus = set()  # positions of the internal C (2nd C in CCGG) on + strand
ccgg_positions_minus = set()
ggccgg_positions_plus = set()
ggccgg_positions_minus = set()

for i in range(genome_len - 3):
    if genome_seq[i:i+4] == "CCGG":
        # Internal C is at position i+1 (0-based)
        ccgg_positions_plus.add(i + 1)
        # Also check GGCCGG context
        if i >= 2 and genome_seq[i-2:i+4] == "GGCCGG":
            ggccgg_positions_plus.add(i + 1)
    # Reverse complement: CCGG rc = CCGG (palindrome)
    # On minus strand, the internal C of CCGG is at the complement position
    if genome_seq[i:i+4] == "CCGG":
        # Complement: position i+2 on minus strand
        ccgg_positions_minus.add(i + 2)
        if i + 4 < genome_len - 1 and genome_seq[i:i+6] == "CCGGCC":
            ggccgg_positions_minus.add(i + 2)

all_ccgg = ccgg_positions_plus | ccgg_positions_minus
all_ggccgg = ggccgg_positions_plus | ggccgg_positions_minus
print(f"  CCGG internal C positions: {len(ccgg_positions_plus)} (+), {len(ccgg_positions_minus)} (-)")
print(f"  GGCCGG subset: {len(ggccgg_positions_plus)} (+), {len(ggccgg_positions_minus)} (-)")

# === Parse pileup BED files ===
BED_COLS = [
    "chrom", "start", "end", "mod_code", "coverage", "strand",
    "thick_start", "thick_end", "rgb", "n_valid", "pct_mod",
    "n_mod", "n_canonical", "n_other_mod", "n_delete", "n_fail", "n_diff", "n_nocall"
]

def parse_pileup(bed_path):
    """Parse modkit pileup BED, return 4mC and 5mC rows for C positions."""
    df = pd.read_csv(bed_path, sep="\t", header=None, names=BED_COLS,
                     dtype={"mod_code": str})
    # 5mC = "m", 4mC = "21839"
    df_5mC = df[df["mod_code"] == "m"].copy()
    df_4mC = df[df["mod_code"] == "21839"].copy()
    return df_4mC, df_5mC

# === Main analysis ===
print("\n=== Processing all samples ===\n")

all_results = []
genome_wide_summary = []

for sample_id, (timepoint, bed_path) in SAMPLES.items():
    print(f"Processing {sample_id} ({timepoint})...")

    if not bed_path.exists():
        print(f"  WARNING: {bed_path} not found, skipping")
        continue

    df_4mC, df_5mC = parse_pileup(bed_path)

    # Genome-wide summary
    n_4mC_all = len(df_4mC[df_4mC["n_valid"] >= MIN_COVERAGE])
    n_5mC_all = len(df_5mC[df_5mC["n_valid"] >= MIN_COVERAGE])
    mean_4mC_all = df_4mC[df_4mC["n_valid"] >= MIN_COVERAGE]["pct_mod"].mean()
    mean_5mC_all = df_5mC[df_5mC["n_valid"] >= MIN_COVERAGE]["pct_mod"].mean()

    # High-confidence sites (≥50% frequency, ≥10x coverage)
    n_4mC_hc = len(df_4mC[(df_4mC["n_valid"] >= MIN_COVERAGE) & (df_4mC["pct_mod"] >= 50)])
    n_5mC_hc = len(df_5mC[(df_5mC["n_valid"] >= MIN_COVERAGE) & (df_5mC["pct_mod"] >= 50)])

    genome_wide_summary.append({
        "sample": sample_id,
        "timepoint": timepoint,
        "4mC_sites_cov10": n_4mC_all,
        "5mC_sites_cov10": n_5mC_all,
        "4mC_mean_pct": round(mean_4mC_all, 2),
        "5mC_mean_pct": round(mean_5mC_all, 2),
        "4mC_HC_sites": n_4mC_hc,
        "5mC_HC_sites": n_5mC_hc,
    })

    # Focus on CCGG sites
    for context_name, pos_plus, pos_minus in [
        ("CCGG", ccgg_positions_plus, ccgg_positions_minus),
        ("GGCCGG", ggccgg_positions_plus, ggccgg_positions_minus),
    ]:
        for mod_name, df_mod, code in [("4mC", df_4mC, "21839"), ("5mC", df_5mC, "m")]:
            # Filter by coverage
            df_cov = df_mod[df_mod["n_valid"] >= MIN_COVERAGE].copy()

            # Match to motif positions
            plus_hits = df_cov[(df_cov["strand"] == "+") & (df_cov["start"].isin(pos_plus))]
            minus_hits = df_cov[(df_cov["strand"] == "-") & (df_cov["start"].isin(pos_minus))]
            hits = pd.concat([plus_hits, minus_hits])

            n_sites = len(hits)
            mean_pct = hits["pct_mod"].mean() if n_sites > 0 else 0
            median_pct = hits["pct_mod"].median() if n_sites > 0 else 0
            n_hc = len(hits[hits["pct_mod"] >= 50])

            all_results.append({
                "sample": sample_id,
                "timepoint": timepoint,
                "context": context_name,
                "mod_type": mod_name,
                "n_sites_cov10": n_sites,
                "mean_pct_mod": round(mean_pct, 2),
                "median_pct_mod": round(median_pct, 2),
                "n_HC_sites_50pct": n_hc,
            })

    print(f"  Done.")

# === Create summary tables ===
df_results = pd.DataFrame(all_results)
df_genome = pd.DataFrame(genome_wide_summary)

# Save genome-wide summary
df_genome.to_csv(OUT_DIR / "tables" / "genome_wide_4mC_vs_5mC_summary.csv", index=False)
print(f"\nGenome-wide summary saved.")

# Pivot for CCGG comparison
print("\n" + "=" * 80)
print("GENOME-WIDE: 4mC vs 5mC (coverage ≥ 10x)")
print("=" * 80)
print(df_genome.to_string(index=False))

# CCGG-specific summary: aggregate across samples by timepoint
print("\n" + "=" * 80)
print("CCGG SITES: 4mC vs 5mC by timepoint (coverage ≥ 10x)")
print("=" * 80)

for context in ["CCGG", "GGCCGG"]:
    print(f"\n--- {context} context ---")
    ctx = df_results[df_results["context"] == context]
    summary = ctx.groupby(["timepoint", "mod_type"]).agg(
        mean_sites=("n_sites_cov10", "mean"),
        mean_pct=("mean_pct_mod", "mean"),
        mean_HC=("n_HC_sites_50pct", "mean"),
    ).round(1)
    print(summary.to_string())

# Save detailed results
df_results.to_csv(OUT_DIR / "tables" / "CCGG_4mC_vs_5mC_detailed.csv", index=False)

# Timepoint-aggregated summary
tp_summary = df_results.groupby(["timepoint", "context", "mod_type"]).agg(
    mean_sites=("n_sites_cov10", "mean"),
    std_sites=("n_sites_cov10", "std"),
    mean_pct=("mean_pct_mod", "mean"),
    mean_HC=("n_HC_sites_50pct", "mean"),
).round(2).reset_index()
tp_summary.to_csv(OUT_DIR / "tables" / "CCGG_4mC_vs_5mC_by_timepoint.csv", index=False)

print(f"\nTables saved to {OUT_DIR / 'tables'}")

# === Figures ===

# Figure 1: 4mC vs 5mC at CCGG — bar chart by timepoint
fig, axes = plt.subplots(1, 3, figsize=(14, 5))

# Panel A: Number of high-confidence sites
for i, context in enumerate(["CCGG", "GGCCGG"]):
    ax = axes[i]
    ctx = tp_summary[tp_summary["context"] == context]
    timepoints = ["T1", "T2", "T3"]
    x = np.arange(len(timepoints))
    width = 0.35

    for j, mod in enumerate(["4mC", "5mC"]):
        mod_data = ctx[ctx["mod_type"] == mod].set_index("timepoint")
        vals = [mod_data.loc[t, "mean_HC"] if t in mod_data.index else 0 for t in timepoints]
        color = "#d62728" if mod == "4mC" else "#1f77b4"
        ax.bar(x + j * width, vals, width, label=mod, color=color, alpha=0.8)

    ax.set_xlabel("Timepoint")
    ax.set_ylabel("High-confidence sites (≥50% freq)")
    ax.set_title(f"{context} context")
    ax.set_xticks(x + width / 2)
    ax.set_xticklabels(timepoints)
    ax.legend()
    ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))

# Panel C: Mean modification frequency at CCGG
ax = axes[2]
for mod, color, marker in [("4mC", "#d62728", "o"), ("5mC", "#1f77b4", "s")]:
    for context, ls in [("CCGG", "-"), ("GGCCGG", "--")]:
        ctx = tp_summary[(tp_summary["context"] == context) & (tp_summary["mod_type"] == mod)]
        ctx = ctx.set_index("timepoint")
        vals = [ctx.loc[t, "mean_pct"] if t in ctx.index else 0 for t in timepoints]
        label = f"{mod} ({context})"
        ax.plot(timepoints, vals, marker=marker, color=color, ls=ls, label=label, linewidth=2)

ax.set_xlabel("Timepoint")
ax.set_ylabel("Mean modification frequency (%)")
ax.set_title("Mean frequency at CCGG/GGCCGG")
ax.legend(fontsize=8)
ax.set_ylim(bottom=0)

plt.suptitle("4mC vs 5mC at CCGG sites\n(Dorado 4mC_5mC@v1 model, both modifications called)",
             fontsize=12, fontweight="bold")
plt.tight_layout()

for ext in ["pdf", "png", "svg"]:
    fig.savefig(OUT_DIR / "figures" / f"4mC_vs_5mC_CCGG_comparison.{ext}", dpi=200, bbox_inches="tight")
plt.close()

# Figure 2: Histogram of modification frequencies at CCGG positions (all samples pooled)
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for idx, context in enumerate(["CCGG", "GGCCGG"]):
    ax = axes[idx]

    # Collect all per-site frequencies from all samples
    all_4mC_freqs = []
    all_5mC_freqs = []

    for sample_id, (timepoint, bed_path) in SAMPLES.items():
        if not bed_path.exists():
            continue
        df_4mC, df_5mC = parse_pileup(bed_path)

        pos_plus = ccgg_positions_plus if context == "CCGG" else ggccgg_positions_plus
        pos_minus = ccgg_positions_minus if context == "CCGG" else ggccgg_positions_minus

        for mod_name, df_mod, freqs_list in [("4mC", df_4mC, all_4mC_freqs), ("5mC", df_5mC, all_5mC_freqs)]:
            df_cov = df_mod[df_mod["n_valid"] >= MIN_COVERAGE]
            plus_hits = df_cov[(df_cov["strand"] == "+") & (df_cov["start"].isin(pos_plus))]
            minus_hits = df_cov[(df_cov["strand"] == "-") & (df_cov["start"].isin(pos_minus))]
            hits = pd.concat([plus_hits, minus_hits])
            freqs_list.extend(hits["pct_mod"].tolist())

    bins = np.arange(0, 105, 5)
    ax.hist(all_4mC_freqs, bins=bins, alpha=0.7, color="#d62728", label=f"4mC (n={len(all_4mC_freqs):,})")
    ax.hist(all_5mC_freqs, bins=bins, alpha=0.7, color="#1f77b4", label=f"5mC (n={len(all_5mC_freqs):,})")
    ax.set_xlabel("Modification frequency (%)")
    ax.set_ylabel("Count (site × sample)")
    ax.set_title(f"{context} sites — frequency distribution")
    ax.legend()

plt.suptitle("Distribution of 4mC vs 5mC frequencies at CCGG sites (all samples pooled)",
             fontsize=11, fontweight="bold")
plt.tight_layout()

for ext in ["pdf", "png", "svg"]:
    fig.savefig(OUT_DIR / "figures" / f"4mC_vs_5mC_frequency_histogram.{ext}", dpi=200, bbox_inches="tight")
plt.close()

print(f"\nFigures saved to {OUT_DIR / 'figures'}")
print("\n=== Analysis complete ===")
