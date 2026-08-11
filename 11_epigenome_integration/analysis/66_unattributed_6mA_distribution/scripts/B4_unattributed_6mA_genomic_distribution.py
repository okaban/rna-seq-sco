"""
B-4: Genomic Distribution of Unattributed 6mA Sites
=====================================================
Purpose: Quantify where non-AAGCCCG 6mA sites reside (BGC, promoter, GCCGGC co-localization)
         to reinforce the paper's claim about biological context of unattributed 6mA.

Input:
  - 23_expanded_motif_search/6mA_final_census.csv  (all T1 6mA with motif assignment)
  - 23_expanded_motif_search/4mC_motif_assignment.csv  (4mC sites for GCCGGC check)
  - 05_annotation/.../gene_master_with_BGC.tsv  (genes + BGC labels)
  - 39_GCCGGC_MTase_reverse_ID/GCF_000203835.1_ASM20383v1_genomic.gff  (for TSS)

Output:
  - tables/B4_summary.tsv
  - tables/B4_bgc_detail.tsv
  - tables/B4_tss_distribution.tsv
  - tables/B4_coloc_detail.tsv
  - figures/B4_tss_distance_histogram.png
  - figures/B4_bgc_pie.png
"""

import sys
import re
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = Path("/Users/okaban/bioinfo/rna-seq")
OUT  = BASE / "11_epigenome_integration/analysis/66_unattributed_6mA_distribution"
OUT_TABLES  = OUT / "tables"
OUT_FIGURES = OUT / "figures"

CENSUS_6mA     = BASE / "11_epigenome_integration/analysis/23_expanded_motif_search/6mA_final_census.csv"
ASSIGN_4mC     = BASE / "11_epigenome_integration/analysis/23_expanded_motif_search/4mC_motif_assignment.csv"
GENE_MASTER    = BASE / "05_annotation/analysis/05_annotation_260128_v1/tables/gene_master_with_BGC.tsv"
GFF            = BASE / "11_epigenome_integration/analysis/39_GCCGGC_MTase_reverse_ID/data/GCF_000203835.1_ASM20383v1_genomic.gff"

GENOME_SIZE    = 8_667_507   # NC_003888.3 chromosome (main contig)
MAIN_CHROM     = "NC_003888.3"

# Analysis parameters
TSS_WINDOW     = 500   # ±500 bp for promoter-proximal definition
COLOC_WINDOW   = 50    # ±50 bp for 4mC/6mA co-localization

# ---------------------------------------------------------------------------
# Step 1: Load 6mA sites, define unattributed set
# ---------------------------------------------------------------------------
logger.info("=== Step 1: Loading 6mA sites ===")

census = pd.read_csv(CENSUS_6mA)
# Focus on T1 (high-confidence, main timepoint), main chromosome
t1_6mA = census[(census["timepoint"] == "T1") & (census["chrom"] == MAIN_CHROM)].copy()

total_6mA   = len(t1_6mA)
aagcccg_6mA = t1_6mA[t1_6mA["at_AAGCCCG"] == True]
unattr_6mA  = t1_6mA[t1_6mA["at_AAGCCCG"] == False]

n_total     = len(t1_6mA)
n_aagcccg   = len(aagcccg_6mA)
n_unattr    = len(unattr_6mA)

logger.info(f"T1 6mA total:      {n_total}")
logger.info(f"AAGCCCG-assigned:  {n_aagcccg}")
logger.info(f"Unattributed:      {n_unattr}")

# Motif breakdown of unattributed
motif_breakdown = unattr_6mA["final_motif"].value_counts()
logger.info(f"Unattributed motif breakdown:\n{motif_breakdown.to_string()}")

# ---------------------------------------------------------------------------
# Step 2: BGC region analysis
# ---------------------------------------------------------------------------
logger.info("\n=== Step 2: BGC region analysis ===")

genes = pd.read_csv(GENE_MASTER, sep="\t")
genes_chr = genes[genes["contig"] == MAIN_CHROM].copy()

# Build BGC interval set from gene-level BGC annotations
bgc_genes = genes_chr.dropna(subset=["bgc_name"])
bgc_intervals = []
for bgc, grp in bgc_genes.groupby("bgc_name"):
    bgc_start = grp["start"].min()
    bgc_end   = grp["end"].max()
    bgc_intervals.append({"bgc_name": bgc, "start": bgc_start, "end": bgc_end,
                           "length": bgc_end - bgc_start})

bgc_df = pd.DataFrame(bgc_intervals).sort_values("start").reset_index(drop=True)
logger.info(f"BGC regions ({len(bgc_df)}):")
for _, row in bgc_df.iterrows():
    logger.info(f"  {row['bgc_name']}: {int(row['start'])}-{int(row['end'])} ({int(row['length'])} bp)")

total_bgc_bp = bgc_df["length"].sum()
bgc_fraction = total_bgc_bp / GENOME_SIZE
logger.info(f"Total BGC bp: {total_bgc_bp:,}  ({bgc_fraction:.4f} = {bgc_fraction*100:.2f}% of genome)")


def in_bgc(pos: int, bgc_df: pd.DataFrame) -> str:
    """Return bgc_name if position falls in any BGC, else None."""
    for _, row in bgc_df.iterrows():
        if row["start"] <= pos <= row["end"]:
            return row["bgc_name"]
    return None


unattr_6mA = unattr_6mA.copy()
unattr_6mA["bgc_name"] = unattr_6mA["position"].apply(lambda p: in_bgc(p, bgc_df))
unattr_6mA["in_bgc"]   = unattr_6mA["bgc_name"].notna()

n_unattr_bgc     = unattr_6mA["in_bgc"].sum()
n_unattr_non_bgc = (~unattr_6mA["in_bgc"]).sum()
obs_frac_bgc     = n_unattr_bgc / n_unattr

# Fisher's exact test: BGC enrichment among unattributed 6mA
# Contingency table:
#             in_BGC   not_in_BGC
# unattr_6mA  a        b
# rest_genome  c       d
genome_bgc_bp     = total_bgc_bp
genome_non_bgc_bp = GENOME_SIZE - total_bgc_bp
# Use site counts (unattr) vs random expectation from genome fraction
a = n_unattr_bgc
b = n_unattr_non_bgc
c = int(genome_bgc_bp)       # "bases in BGC"
d = int(genome_non_bgc_bp)   # "bases outside BGC"

odds_bgc, p_bgc = fisher_exact([[a, b], [c, d]])
expected_bgc = n_unattr * bgc_fraction
enrichment_bgc = obs_frac_bgc / bgc_fraction if bgc_fraction > 0 else float("nan")

logger.info(f"\nUnattr 6mA in BGC:     {n_unattr_bgc} / {n_unattr} ({obs_frac_bgc*100:.1f}%)")
logger.info(f"Expected (genome frac): {expected_bgc:.1f} ({bgc_fraction*100:.2f}%)")
logger.info(f"Enrichment:             {enrichment_bgc:.3f}x")
logger.info(f"Fisher p-value:         {p_bgc:.4e}")

bgc_detail_unattr = unattr_6mA[unattr_6mA["in_bgc"]]["bgc_name"].value_counts().reset_index()
bgc_detail_unattr.columns = ["bgc_name", "n_unattr_6mA"]
bgc_detail_unattr = bgc_detail_unattr.merge(bgc_df[["bgc_name","length"]], on="bgc_name")
bgc_detail_unattr["bgc_frac_genome"] = bgc_detail_unattr["length"] / GENOME_SIZE
bgc_detail_unattr["obs_frac_unattr"] = bgc_detail_unattr["n_unattr_6mA"] / n_unattr
bgc_detail_unattr["enrichment"] = bgc_detail_unattr["obs_frac_unattr"] / bgc_detail_unattr["bgc_frac_genome"]

# ---------------------------------------------------------------------------
# Step 3: TSS proximity analysis
# ---------------------------------------------------------------------------
logger.info("\n=== Step 3: TSS proximity analysis ===")

# Parse GFF to extract TSS positions from CDS/gene features on main chrom
tss_records = []
with open(GFF) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        cols = line.rstrip().split("\t")
        if len(cols) < 9:
            continue
        seqid, source, ftype, start_s, end_s, score, strand, phase, attrs = cols
        if seqid != MAIN_CHROM:
            continue
        if ftype not in ("gene", "mRNA", "CDS"):
            continue
        start_pos = int(start_s)
        end_pos   = int(end_s)
        tss = start_pos if strand == "+" else end_pos
        tss_records.append({"chrom": seqid, "tss": tss, "strand": strand, "ftype": ftype})

tss_df = pd.DataFrame(tss_records)
# Use gene-level TSS (deduplicate)
tss_gene = tss_df[tss_df["ftype"] == "gene"][["tss"]].drop_duplicates().reset_index(drop=True)
tss_positions = tss_gene["tss"].values
logger.info(f"TSS positions (gene-level): {len(tss_positions)}")

# For each unattr 6mA, find minimum TSS distance
def min_tss_dist(pos: int, tss_arr: np.ndarray) -> int:
    return int(np.min(np.abs(tss_arr - pos)))

logger.info("Computing TSS distances (may take ~10s)...")
unattr_6mA = unattr_6mA.copy()
unattr_6mA["min_tss_dist"] = unattr_6mA["position"].apply(
    lambda p: min_tss_dist(p, tss_positions)
)

# TSS-proximal = within ±500 bp
n_proximal = (unattr_6mA["min_tss_dist"] <= TSS_WINDOW).sum()
frac_proximal = n_proximal / n_unattr

# Expected fraction: 2*500 bp windows per TSS, but overlapping possible
# Use simple genome-fraction approach: (n_tss * 2 * TSS_WINDOW) / GENOME_SIZE, capped at 1
expected_proximal_frac = min(1.0, len(tss_positions) * 2 * TSS_WINDOW / GENOME_SIZE)
expected_proximal_n    = n_unattr * expected_proximal_frac

# Fisher's exact test for TSS proximity enrichment
# Rows: proximal/distal; Cols: unattr_6mA / "genome background"
a_tss = n_proximal
b_tss = n_unattr - n_proximal
c_tss = int(len(tss_positions) * 2 * TSS_WINDOW)           # approx proximal genome bp
d_tss = max(1, GENOME_SIZE - c_tss)                          # distal genome bp
odds_tss, p_tss = fisher_exact([[a_tss, b_tss], [c_tss, d_tss]])
enrichment_tss = frac_proximal / expected_proximal_frac if expected_proximal_frac > 0 else float("nan")

logger.info(f"Unattr 6mA proximal (±{TSS_WINDOW}bp): {n_proximal} / {n_unattr} ({frac_proximal*100:.1f}%)")
logger.info(f"Expected:                              {expected_proximal_n:.1f} ({expected_proximal_frac*100:.2f}%)")
logger.info(f"Enrichment:                            {enrichment_tss:.3f}x")
logger.info(f"Fisher p-value:                        {p_tss:.4e}")

# TSS distance distribution bins
dist_bins = [0, 100, 250, 500, 1000, 2000, 5000, 20000, np.inf]
bin_labels = ["0-100", "100-250", "250-500", "500-1k", "1k-2k", "2k-5k", "5k-20k", ">20k"]
unattr_6mA["dist_bin"] = pd.cut(
    unattr_6mA["min_tss_dist"], bins=dist_bins, labels=bin_labels, right=True
)
dist_counts = unattr_6mA["dist_bin"].value_counts().sort_index()

# ---------------------------------------------------------------------------
# Step 4: GCCGGC 4mC co-localization
# ---------------------------------------------------------------------------
logger.info("\n=== Step 4: GCCGGC 4mC co-localization ===")

m4c = pd.read_csv(ASSIGN_4mC)
# GCCGGC motif → look for TGGCCGGC (the observed 4mC motif containing GCCGGC)
gccggc_4mC = m4c[
    (m4c["timepoint"] == "T1") &
    (m4c["chrom"] == MAIN_CHROM) &
    (m4c["assigned_motif"].str.contains("GCCGGC", na=False))
].copy()

n_gccggc = len(gccggc_4mC)
gccggc_pos = gccggc_4mC["position"].values

logger.info(f"GCCGGC-motif 4mC sites (T1): {n_gccggc}")

# For each unattr 6mA, find nearest GCCGGC 4mC site
def min_coloc_dist(pos: int, other_arr: np.ndarray) -> int:
    if len(other_arr) == 0:
        return 999999
    return int(np.min(np.abs(other_arr - pos)))

unattr_6mA["gccggc_dist"] = unattr_6mA["position"].apply(
    lambda p: min_coloc_dist(p, gccggc_pos)
)
unattr_6mA["coloc_gccggc"] = unattr_6mA["gccggc_dist"] <= COLOC_WINDOW

n_coloc      = unattr_6mA["coloc_gccggc"].sum()
frac_coloc   = n_coloc / n_unattr

# Expected co-localization under independence:
# P(6mA within ±50bp of any 4mC) ≈ (2*50*n_gccggc) / GENOME_SIZE
expected_coloc_frac = min(1.0, 2 * COLOC_WINDOW * n_gccggc / GENOME_SIZE)
expected_coloc_n    = n_unattr * expected_coloc_frac

# Fisher's exact
a_co = n_coloc
b_co = n_unattr - n_coloc
c_co = int(2 * COLOC_WINDOW * n_gccggc)
d_co = max(1, GENOME_SIZE - c_co)
odds_co, p_co = fisher_exact([[a_co, b_co], [c_co, d_co]])
enrichment_co = frac_coloc / expected_coloc_frac if expected_coloc_frac > 0 else float("nan")

logger.info(f"Unattr 6mA co-localizing with GCCGGC 4mC (±{COLOC_WINDOW}bp): "
            f"{n_coloc} / {n_unattr} ({frac_coloc*100:.2f}%)")
logger.info(f"Expected (independence):  {expected_coloc_n:.2f} ({expected_coloc_frac*100:.4f}%)")
logger.info(f"Enrichment:              {enrichment_co:.3f}x")
logger.info(f"Fisher p-value:          {p_co:.4e}")

# ---------------------------------------------------------------------------
# Save tables
# ---------------------------------------------------------------------------
logger.info("\n=== Saving tables ===")

# Summary table
summary_rows = [
    # --- Step 1 ---
    {"analysis": "6mA_total_T1",       "category": "all",       "n": n_total,    "pct_unattr": float("nan"), "expected_pct": float("nan"), "enrichment": float("nan"), "fisher_p": float("nan")},
    {"analysis": "6mA_AAGCCCG_T1",    "category": "AAGCCCG",   "n": n_aagcccg,  "pct_unattr": float("nan"), "expected_pct": float("nan"), "enrichment": float("nan"), "fisher_p": float("nan")},
    {"analysis": "6mA_unattr_T1",      "category": "unattr",    "n": n_unattr,   "pct_unattr": float("nan"), "expected_pct": float("nan"), "enrichment": float("nan"), "fisher_p": float("nan")},
    # --- Step 2 ---
    {"analysis": "BGC_enrichment",     "category": "in_BGC",    "n": n_unattr_bgc,   "pct_unattr": round(obs_frac_bgc*100, 2),     "expected_pct": round(bgc_fraction*100, 2),           "enrichment": round(enrichment_bgc, 3),  "fisher_p": float(p_bgc)},
    {"analysis": "BGC_enrichment",     "category": "non_BGC",   "n": n_unattr_non_bgc, "pct_unattr": round((1-obs_frac_bgc)*100, 2), "expected_pct": round((1-bgc_fraction)*100, 2),      "enrichment": float("nan"),              "fisher_p": float("nan")},
    # --- Step 3 ---
    {"analysis": "TSS_proximity",      "category": "proximal",  "n": n_proximal,      "pct_unattr": round(frac_proximal*100, 2),    "expected_pct": round(expected_proximal_frac*100, 2), "enrichment": round(enrichment_tss, 3),  "fisher_p": float(p_tss)},
    {"analysis": "TSS_proximity",      "category": "distal",    "n": n_unattr-n_proximal, "pct_unattr": round((1-frac_proximal)*100, 2), "expected_pct": round((1-expected_proximal_frac)*100, 2), "enrichment": float("nan"), "fisher_p": float("nan")},
    # --- Step 4 ---
    {"analysis": "GCCGGC_coloc",       "category": "colocalized","n": n_coloc,         "pct_unattr": round(frac_coloc*100, 2),      "expected_pct": round(expected_coloc_frac*100, 4),    "enrichment": round(enrichment_co, 3),   "fisher_p": float(p_co)},
    {"analysis": "GCCGGC_coloc",       "category": "non_coloc",  "n": n_unattr-n_coloc,"pct_unattr": round((1-frac_coloc)*100, 2),  "expected_pct": round((1-expected_coloc_frac)*100, 4),"enrichment": float("nan"),              "fisher_p": float("nan")},
]
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(OUT_TABLES / "B4_summary.tsv", sep="\t", index=False)
logger.info(f"Saved: {OUT_TABLES/'B4_summary.tsv'}")

# BGC detail
bgc_detail_unattr.to_csv(OUT_TABLES / "B4_bgc_detail.tsv", sep="\t", index=False)
logger.info(f"Saved: {OUT_TABLES/'B4_bgc_detail.tsv'}")

# TSS distribution
tss_dist_df = dist_counts.reset_index()
tss_dist_df.columns = ["dist_bin", "n_unattr_6mA"]
tss_dist_df["pct"] = (tss_dist_df["n_unattr_6mA"] / n_unattr * 100).round(2)
tss_dist_df.to_csv(OUT_TABLES / "B4_tss_distribution.tsv", sep="\t", index=False)
logger.info(f"Saved: {OUT_TABLES/'B4_tss_distribution.tsv'}")

# Co-localization detail
coloc_sites = unattr_6mA[unattr_6mA["coloc_gccggc"]][
    ["chrom","position","strand","final_motif","gccggc_dist","bgc_name","min_tss_dist"]
].sort_values("gccggc_dist")
coloc_sites.to_csv(OUT_TABLES / "B4_coloc_detail.tsv", sep="\t", index=False)
logger.info(f"Saved: {OUT_TABLES/'B4_coloc_detail.tsv'}")

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
logger.info("\n=== Generating figures ===")

# Fig 1: TSS distance histogram
fig, ax = plt.subplots(figsize=(8, 4))
bins = np.logspace(np.log10(1), np.log10(GENOME_SIZE/2), 50)
ax.hist(unattr_6mA["min_tss_dist"], bins=bins, color="#4C72B0", edgecolor="white", linewidth=0.3)
ax.axvline(TSS_WINDOW, color="#DD4444", linestyle="--", linewidth=1.5, label=f"±{TSS_WINDOW} bp window")
ax.set_xscale("log")
ax.set_xlabel("Distance to nearest TSS (bp, log scale)", fontsize=11)
ax.set_ylabel("Number of unattributed 6mA sites", fontsize=11)
ax.set_title("TSS Distance Distribution of Unattributed 6mA Sites (T1)", fontsize=12)
ax.legend(fontsize=10)
fig.tight_layout()
fig.savefig(OUT_FIGURES / "B4_tss_distance_histogram.png", dpi=200)
plt.close()
logger.info(f"Saved: {OUT_FIGURES/'B4_tss_distance_histogram.png'}")

# Fig 2: BGC pie / bar
fig, axes = plt.subplots(1, 2, figsize=(10, 4))

# Pie chart
ax = axes[0]
labels_pie  = ["In BGC", "Non-BGC"]
sizes_pie   = [n_unattr_bgc, n_unattr_non_bgc]
colors_pie  = ["#E07B54", "#A8C5DA"]
ax.pie(sizes_pie, labels=labels_pie, colors=colors_pie, autopct="%1.1f%%",
       startangle=90, textprops={"fontsize": 11})
ax.set_title(f"Unattributed 6mA by BGC (n={n_unattr})", fontsize=12)

# Per-BGC bar chart
ax2 = axes[1]
bgc_names  = bgc_detail_unattr["bgc_name"].tolist() + (["(non-BGC)"] if n_unattr_non_bgc > 0 else [])
bgc_counts = bgc_detail_unattr["n_unattr_6mA"].tolist() + ([n_unattr_non_bgc] if n_unattr_non_bgc > 0 else [])
colors_bar = ["#E07B54"] * len(bgc_detail_unattr) + (["#A8C5DA"] if n_unattr_non_bgc > 0 else [])
bars = ax2.bar(bgc_names, bgc_counts, color=colors_bar, edgecolor="white")
ax2.set_xlabel("BGC", fontsize=11)
ax2.set_ylabel("Unattributed 6mA count", fontsize=11)
ax2.set_title("Per-BGC Breakdown", fontsize=12)
for bar, cnt in zip(bars, bgc_counts):
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3, str(cnt),
             ha="center", va="bottom", fontsize=9)

fig.tight_layout()
fig.savefig(OUT_FIGURES / "B4_bgc_distribution.png", dpi=200)
plt.close()
logger.info(f"Saved: {OUT_FIGURES/'B4_bgc_distribution.png'}")

# ---------------------------------------------------------------------------
# Print final summary to stdout
# ---------------------------------------------------------------------------
print("\n" + "="*65)
print("B-4: Genomic Distribution of Unattributed 6mA Sites")
print("="*65)
print(f"\n[Step 1] 6mA Site Counts (T1, {MAIN_CHROM})")
print(f"  Total 6mA:          {n_total:>6}")
print(f"  AAGCCCG-assigned:   {n_aagcccg:>6}  ({n_aagcccg/n_total*100:.1f}%)")
print(f"  Unattributed:       {n_unattr:>6}  ({n_unattr/n_total*100:.1f}%)")
print(f"\n  Unattributed motif sub-breakdown:")
for motif, cnt in motif_breakdown.items():
    print(f"    {motif:<20}: {cnt:>5}  ({cnt/n_unattr*100:.1f}%)")

print(f"\n[Step 2] BGC Region Analysis")
print(f"  Genome size:                    {GENOME_SIZE:>12,} bp")
print(f"  Total BGC coverage:             {total_bgc_bp:>12,} bp  ({bgc_fraction*100:.2f}%)")
print(f"  Unattr 6mA in any BGC:          {n_unattr_bgc:>6} / {n_unattr} ({obs_frac_bgc*100:.1f}%)")
print(f"  Expected (genome fraction):     {expected_bgc:>8.1f}  ({bgc_fraction*100:.2f}%)")
print(f"  Enrichment:                     {enrichment_bgc:.3f}x")
print(f"  Fisher's exact p-value:         {p_bgc:.3e}")
print(f"\n  Per-BGC detail:")
for _, row in bgc_detail_unattr.iterrows():
    print(f"    {row['bgc_name']:<6}: {int(row['n_unattr_6mA'])} sites  "
          f"(enrichment {row['enrichment']:.2f}x)")

print(f"\n[Step 3] TSS Proximity Analysis (±{TSS_WINDOW} bp)")
print(f"  TSS count (genes):              {len(tss_positions):>6}")
print(f"  Expected proximal fraction:     {expected_proximal_frac*100:.2f}%")
print(f"  Unattr 6mA proximal:            {n_proximal:>6} / {n_unattr} ({frac_proximal*100:.1f}%)")
print(f"  Enrichment:                     {enrichment_tss:.3f}x")
print(f"  Fisher's exact p-value:         {p_tss:.3e}")
print(f"\n  Distance distribution:")
for _, row in tss_dist_df.iterrows():
    bar = "#" * int(row["pct"] / 2)
    print(f"    {str(row['dist_bin']):<12}: {int(row['n_unattr_6mA']):>5}  ({row['pct']:>5.1f}%)  {bar}")

print(f"\n[Step 4] GCCGGC 4mC Co-localization (±{COLOC_WINDOW} bp)")
print(f"  GCCGGC 4mC sites (T1):         {n_gccggc:>6}")
print(f"  Unattr 6mA co-localizing:      {n_coloc:>6} / {n_unattr} ({frac_coloc*100:.2f}%)")
print(f"  Expected (independence):       {expected_coloc_n:.2f} ({expected_coloc_frac*100:.4f}%)")
print(f"  Enrichment:                    {enrichment_co:.3f}x")
print(f"  Fisher's exact p-value:        {p_co:.3e}")

print("\n" + "="*65)
print("Output files:")
for f in sorted(OUT_TABLES.iterdir()):
    print(f"  {f}")
for f in sorted(OUT_FIGURES.iterdir()):
    print(f"  {f}")
print("="*65)
