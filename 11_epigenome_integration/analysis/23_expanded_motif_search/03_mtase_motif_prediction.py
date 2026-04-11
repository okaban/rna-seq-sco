#!/usr/bin/env python3
"""
03: Orphan MTase → Motif Prediction
====================================
Phase 3: Correlate MTase expression dynamics with motif temporal patterns
to predict which MTase produces which motif.
"""

import pandas as pd
import numpy as np
from scipy import stats
import os

BASE_DIR = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration"
OUT_DIR = f"{BASE_DIR}/analysis/23_expanded_motif_search"
GENE_MASTER = "/Users/okaban/bioinfo/rna-seq/05_annotation/analysis/05_annotation_260128_v1/tables/gene_master_DESeq2.tsv"

# ============================================================
# Load data
# ============================================================
mtase_expr = pd.read_csv(
    f"{BASE_DIR}/analysis/11_rm_system_identification/mtase_genes_with_expression.csv")
rm_systems = pd.read_csv(
    f"{BASE_DIR}/analysis/11_rm_system_identification/rm_systems.csv")
motif_dynamics = pd.read_csv(f"{OUT_DIR}/motif_temporal_dynamics.csv")

# Load gene master for normalized counts
gene_master = pd.read_csv(GENE_MASTER, sep="\t")

print("=" * 70)
print("PHASE 3: ORPHAN MTase → MOTIF PREDICTION")
print("=" * 70)


# ============================================================
# 1. MTase expression profiles (normalized counts at T1, T2, T3)
# ============================================================
print("\n" + "=" * 70)
print("1. MTase EXPRESSION PROFILES")
print("=" * 70)

# Get normalized count columns
norm_cols_T1 = [c for c in gene_master.columns if c.startswith("norm_M145_1")]
norm_cols_T2 = [c for c in gene_master.columns if c.startswith("norm_M145_2")]
norm_cols_T3 = [c for c in gene_master.columns if c.startswith("norm_M145_3")]

# Merge MTase genes with gene master for normalized counts
mtase_with_counts = mtase_expr.merge(
    gene_master[["gene_id"] + norm_cols_T1 + norm_cols_T2 + norm_cols_T3],
    left_on="locus_tag", right_on="gene_id", how="left")

# Calculate mean expression per timepoint
mtase_with_counts["mean_T1"] = mtase_with_counts[norm_cols_T1].mean(axis=1)
mtase_with_counts["mean_T2"] = mtase_with_counts[norm_cols_T2].mean(axis=1)
mtase_with_counts["mean_T3"] = mtase_with_counts[norm_cols_T3].mean(axis=1)

# Expression trend classification
def classify_trend(row):
    t1, t2, t3 = row["mean_T1"], row["mean_T2"], row["mean_T3"]
    if pd.isna(t1) or pd.isna(t2) or pd.isna(t3):
        return "NO_DATA"
    if t1 == 0 and t2 == 0 and t3 == 0:
        return "SILENT"
    if t2 > t1 * 1.5 and t3 > t1 * 1.5:
        return "UP_sustained"
    elif t2 < t1 * 0.67 and t3 < t1 * 0.67:
        return "DOWN_sustained"
    elif t2 < t1 * 0.67 and t3 > t2 * 1.5:
        return "DOWN_then_UP"
    elif t2 > t1 * 1.5 and t3 < t2 * 0.67:
        return "UP_then_DOWN"
    elif t3 > t1 * 2:
        return "LATE_UP"
    elif t3 < t1 * 0.5:
        return "LATE_DOWN"
    else:
        return "STABLE"

mtase_with_counts["expression_trend"] = mtase_with_counts.apply(classify_trend, axis=1)

# DNA methyltransferases only (exclude repair enzymes and SAM-dependent misc)
dna_mtases = mtase_with_counts[
    mtase_with_counts["product"].str.contains(
        "methyltransferase|methylase|PglX", case=False, na=False) &
    ~mtase_with_counts["product"].str.contains(
        "glycosylase|kinase|polymerase|RadA|repair|phosphatase|coupling", case=False, na=False)
].copy()

print(f"\nDNA methyltransferases with expression data: {len(dna_mtases)}")
print(f"\n{'Locus':<14} {'Product':<45} {'T1':>8} {'T2':>8} {'T3':>8} {'Trend':<18} {'Methyl'}")
print("-" * 120)
for _, row in dna_mtases.sort_values("mean_T1", ascending=False).iterrows():
    t1 = row["mean_T1"]
    t2 = row["mean_T2"]
    t3 = row["mean_T3"]
    if pd.isna(t1):
        continue
    print(f"{row['locus_tag']:<14} {row['product'][:44]:<45} "
          f"{t1:>8.1f} {t2:>8.1f} {t3:>8.1f} {row['expression_trend']:<18} {row['predicted_methyl_type']}")


# ============================================================
# 2. Motif temporal dynamics (from Phase 1)
# ============================================================
print("\n" + "=" * 70)
print("2. MOTIF TEMPORAL DYNAMICS (SITE COUNTS)")
print("=" * 70)

# Already have motif_dynamics from Phase 1
# Also add the key new motifs from Phase 2
# Recalculate for expanded motif list using the original methylation data

print(f"\n{'Mod':<5} {'Motif':<14} {'T1':>6} {'T2':>6} {'T3':>6} {'T1%':>7} {'T2%':>7} {'T3%':>7} {'Pattern'}")
print("-" * 80)
for _, row in motif_dynamics.iterrows():
    pattern = ""
    t1p = row["T1_pct"]
    t3p = row["T3_pct"]
    if t3p > t1p * 1.3:
        pattern = "↑ increasing"
    elif t3p < t1p * 0.7:
        pattern = "↓ decreasing"
    else:
        pattern = "→ stable"
    print(f"{row['mod_type']:<5} {row['motif']:<14} "
          f"{row['T1_count']:>6} {row['T2_count']:>6} {row['T3_count']:>6} "
          f"{row['T1_pct']:>6.1f}% {row['T2_pct']:>6.1f}% {row['T3_pct']:>6.1f}% {pattern}")


# ============================================================
# 3. Expression-motif correlation matrix
# ============================================================
print("\n" + "=" * 70)
print("3. MTase EXPRESSION ↔ MOTIF DYNAMICS CORRELATION")
print("=" * 70)

# For each MTase, create a 3-point expression vector (T1, T2, T3)
# For each motif, create a 3-point site count vector (T1, T2, T3)
# Compute Pearson correlation

# MTase expression vectors
mtase_vectors = {}
for _, row in dna_mtases.iterrows():
    if pd.isna(row["mean_T1"]):
        continue
    mtase_vectors[row["locus_tag"]] = np.array([
        row["mean_T1"], row["mean_T2"], row["mean_T3"]
    ])

# Motif count vectors (aggregate across mod types for same motif)
motif_vectors = {}
for _, row in motif_dynamics.iterrows():
    key = f"{row['mod_type']}_{row['motif']}"
    motif_vectors[key] = np.array([row["T1_count"], row["T2_count"], row["T3_count"]])

# Compute correlations
corr_results = []
for mtase_id, mtase_vec in mtase_vectors.items():
    # Normalize
    if mtase_vec.std() == 0:
        continue
    for motif_key, motif_vec in motif_vectors.items():
        if motif_vec.std() == 0:
            continue
        # Pearson correlation on 3 points (note: low statistical power)
        r, p = stats.pearsonr(mtase_vec, motif_vec.astype(float))
        corr_results.append({
            "mtase": mtase_id,
            "motif": motif_key,
            "pearson_r": round(r, 3),
            "p_value": round(p, 4),
            "mtase_trend": dna_mtases[dna_mtases["locus_tag"] == mtase_id]["expression_trend"].values[0],
            "mtase_product": dna_mtases[dna_mtases["locus_tag"] == mtase_id]["product"].values[0][:40],
        })

corr_df = pd.DataFrame(corr_results)
corr_df.to_csv(f"{OUT_DIR}/mtase_motif_correlation.csv", index=False)

# Show strong correlations (|r| > 0.8)
strong = corr_df[abs(corr_df["pearson_r"]) > 0.8].sort_values("pearson_r", ascending=False)
print(f"\nStrong correlations (|r| > 0.8) between MTase expression and motif counts:")
print(f"{'MTase':<14} {'Motif':<18} {'r':>6} {'p':>8} {'Trend':<18} {'Product'}")
print("-" * 100)
for _, row in strong.iterrows():
    print(f"{row['mtase']:<14} {row['motif']:<18} {row['pearson_r']:>6.3f} {row['p_value']:>8.4f} "
          f"{row['mtase_trend']:<18} {row['mtase_product']}")


# ============================================================
# 4. Orphan MTase identification
# ============================================================
print("\n" + "=" * 70)
print("4. ORPHAN MTase IDENTIFICATION")
print("=" * 70)
print("(MTases lacking cognate restriction enzyme partners)")

# Known R-M system assignments
assigned_mtases = set(rm_systems["locus_tag"].values)

print(f"\nR-M system assigned MTases: {assigned_mtases}")
print(f"\nAll DNA MTases:")
print(f"{'Locus':<14} {'Product':<45} {'In R-M?':>8} {'Methyl':>10} {'Status'}")
print("-" * 100)

orphan_list = []
for _, row in dna_mtases.iterrows():
    in_rm = "YES" if row["locus_tag"] in assigned_mtases else "NO"
    status = "Assigned" if in_rm == "YES" else "ORPHAN"

    # Check genomic neighbors for restriction enzyme genes
    locus = row["locus_tag"]
    gene_row = gene_master[gene_master["gene_id"] == locus]
    if not gene_row.empty:
        gene_start = gene_row["start"].values[0]
        gene_end = gene_row["end"].values[0]
        # Look for restriction-enzyme related genes within ±5kb
        nearby = gene_master[
            (abs(gene_master["start"] - gene_start) < 5000) |
            (abs(gene_master["end"] - gene_end) < 5000)
        ]
        restriction_nearby = nearby[nearby["product"].str.contains(
            "restrict|endonuclease|HsdR|HsdS|specificity", case=False, na=False)]
        if len(restriction_nearby) > 0:
            neighbor_info = "; ".join(restriction_nearby["gene_id"].values[:3])
            status += f" (restriction nearby: {neighbor_info})"
    else:
        status += " (not in gene_master)"

    print(f"{row['locus_tag']:<14} {row['product'][:44]:<45} {in_rm:>8} "
          f"{row['predicted_methyl_type']:>10} {status}")

    orphan_list.append({
        "locus_tag": row["locus_tag"],
        "product": row["product"],
        "in_rm_system": in_rm == "YES",
        "predicted_methyl_type": row["predicted_methyl_type"],
        "status": status,
        "mean_T1": row["mean_T1"] if not pd.isna(row["mean_T1"]) else None,
        "mean_T2": row["mean_T2"] if not pd.isna(row["mean_T2"]) else None,
        "mean_T3": row["mean_T3"] if not pd.isna(row["mean_T3"]) else None,
        "expression_trend": row["expression_trend"],
    })

orphan_df = pd.DataFrame(orphan_list)
orphan_df.to_csv(f"{OUT_DIR}/orphan_mtase_inventory.csv", index=False)


# ============================================================
# 5. SC_RS17645 neighborhood analysis (Type I subunit search)
# ============================================================
print("\n" + "=" * 70)
print("5. SC_RS17645 GENOMIC NEIGHBORHOOD (Type I R-M subunit search)")
print("=" * 70)

sc17645 = gene_master[gene_master["gene_id"] == "SC_RS17645"]
if not sc17645.empty:
    sc_start = sc17645["start"].values[0]
    sc_end = sc17645["end"].values[0]

    # Look within ±10kb
    neighbors = gene_master[
        (gene_master["start"] > sc_start - 10000) &
        (gene_master["end"] < sc_end + 10000)
    ].sort_values("start")

    print(f"\nGenes within ±10kb of SC_RS17645 (pos {sc_start}-{sc_end}):")
    print(f"{'Gene':<14} {'Start':>8} {'End':>8} {'Strand':>7} {'Product'}")
    print("-" * 80)
    for _, row in neighbors.iterrows():
        marker = " <<<" if row["gene_id"] == "SC_RS17645" else ""
        print(f"{row['gene_id']:<14} {row['start']:>8} {row['end']:>8} {row['strand']:>7} "
              f"{str(row['product'])[:50]}{marker}")


# ============================================================
# 6. BREX-2 PglX neighborhood analysis
# ============================================================
print("\n" + "=" * 70)
print("6. BREX-2 PglX GENOMIC NEIGHBORHOODS")
print("=" * 70)

for pglx_id in ["SC_RS28835", "SC_RS35335"]:
    pglx = gene_master[gene_master["gene_id"] == pglx_id]
    if pglx.empty:
        print(f"\n{pglx_id}: Not found in gene_master")
        continue

    p_start = pglx["start"].values[0]
    p_end = pglx["end"].values[0]

    neighbors = gene_master[
        (gene_master["start"] > p_start - 15000) &
        (gene_master["end"] < p_end + 15000)
    ].sort_values("start")

    print(f"\n{pglx_id} neighborhood (±15kb, pos {p_start}-{p_end}):")
    print(f"{'Gene':<14} {'Start':>8} {'End':>8} {'Strand':>7} {'Product'}")
    print("-" * 80)
    for _, row in neighbors.iterrows():
        marker = " <<<" if row["gene_id"] == pglx_id else ""
        brex_marker = " [BREX?]" if any(kw in str(row["product"]).lower()
                                         for kw in ["brex", "pgl", "phosphatase"]) else ""
        print(f"{row['gene_id']:<14} {row['start']:>8} {row['end']:>8} {row['strand']:>7} "
              f"{str(row['product'])[:45]}{marker}{brex_marker}")


# ============================================================
# 7. MTase → Motif prediction mapping
# ============================================================
print("\n" + "=" * 70)
print("7. MTase → MOTIF PREDICTION MAPPING")
print("=" * 70)

predictions = [
    {
        "mtase": "SC_RS17645",
        "type": "Type I HsdM",
        "predicted_motif": "AAGCCCG",
        "modification": "6mA (+4mC dual)",
        "evidence": [
            "MEME-1 (6mA): E=3.1e-256, 656 sites",
            "Expression correlation: DOWN at T2, correlates with AAGCCCG site reduction",
            "Domain: N6_Mtase + TRD (Type I architecture)",
            "Structural homology to PacII (Type I M subunit)",
            "BLAST: Streptomyces-specific (100% query coverage)",
            "NEW: Dual modification (4mC+6mA) at AAGCCCG confirmed",
        ],
        "confidence": "HIGH",
    },
    {
        "mtase": "SC_RS19770 / SC_RS36410",
        "type": "Dcm-like (cytosine MTase)",
        "predicted_motif": "CCGG (as GGCCGG/TGGCCGGC)",
        "modification": "4mC (or 5mC?)",
        "evidence": [
            "Both upregulated at T3 (SC_RS19770: log2FC=+2.32; SC_RS36410: log2FC=+5.34)",
            "Known Dcm-like cytosine MTase domain",
            "CCGG sites show strong 4mC (75.6% of all 4mC)",
            "Paradox: Dcm normally produces 5mC, not 4mC",
        ],
        "confidence": "MEDIUM (4mC vs 5mC ambiguity)",
    },
    {
        "mtase": "SC_RS28835 / SC_RS35335",
        "type": "BREX-2 PglX",
        "predicted_motif": "CCGKCA or unknown 6mA motif",
        "modification": "6mA",
        "evidence": [
            "Annotated as adenine-specific DNA-methyltransferase PglX",
            "SC_RS28835: slightly downregulated (stable)",
            "SC_RS35335: downregulated at T3 (log2FC=-1.0)",
            "BREX systems typically methylate 6mA at 5-6bp sites",
            "CCGKCA (MEME-2 core) is a novel 6mA motif candidate (153 sites)",
            "Expression pattern partially consistent with CCGKCA temporal dynamics",
        ],
        "confidence": "LOW-MEDIUM (needs validation)",
    },
    {
        "mtase": "SC_RS19670 / SC_RS36625",
        "type": "DNA-methyltransferase (uncharacterized)",
        "predicted_motif": "Unknown (LATE_UP pattern)",
        "modification": "Unknown",
        "evidence": [
            "SC_RS19670: LATE_UP (log2FC=+3.51 T3vsT1)",
            "SC_RS36625: LATE_UP (log2FC=+3.87 T3vsT1)",
            "Both strongly induced at T3",
            "May contribute to T3-specific methylation patterns",
            "Near SC_RS19770 (Dcm-like) — possible operonic relationship",
        ],
        "confidence": "LOW (expression-based only)",
    },
    {
        "mtase": "SC_RS24685",
        "type": "Class I SAM-dependent",
        "predicted_motif": "Unknown",
        "modification": "Unknown",
        "evidence": [
            "Strongly downregulated: log2FC=-3.22 (T2vsT1)",
            "Pattern matches SC_RS17645 (both DOWN at T2)",
            "Could be auxiliary to AAGCCCG modification or independent",
        ],
        "confidence": "LOW",
    },
    {
        "mtase": "SC_RS03950",
        "type": "Class I SAM-dependent",
        "predicted_motif": "Unknown (stable expression)",
        "modification": "Unknown",
        "evidence": [
            "Stable expression across all timepoints",
            "May contribute to constitutive methylation patterns",
        ],
        "confidence": "LOW",
    },
]

print(f"\n{'MTase':<28} {'Type':<22} {'Predicted Motif':<20} {'Mod':<12} {'Confidence'}")
print("-" * 100)
for p in predictions:
    print(f"{p['mtase']:<28} {p['type']:<22} {p['predicted_motif']:<20} "
          f"{p['modification']:<12} {p['confidence']}")
    for ev in p["evidence"][:3]:
        print(f"  - {ev}")
    print()

# Save predictions
pred_df = pd.DataFrame([{k: v if k != "evidence" else "; ".join(v)
                          for k, v in p.items()} for p in predictions])
pred_df.to_csv(f"{OUT_DIR}/mtase_motif_predictions.csv", index=False)

print(f"\nAll results saved to: {OUT_DIR}/")
print("Done.")
