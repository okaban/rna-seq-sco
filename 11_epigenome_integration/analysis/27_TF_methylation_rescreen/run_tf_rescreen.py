#!/usr/bin/env python3
"""
H4: Full re-screen of 37 TFs for methylation-expression coordination
using CORRECTED locus_tags, with SARP upstream regulator focus.

Author: Claude Opus 4.6
Date: 2026-02-24
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import re
import os
from pathlib import Path

# ============================================================
# PATHS
# ============================================================
BASE = Path("/Users/okaban/bioinfo/rna-seq")
ANALYSIS = BASE / "11_epigenome_integration/analysis/27_TF_methylation_rescreen"
FIGURES = ANALYSIS / "figures"
TABLES = ANALYSIS / "tables"

TF_MASTER = BASE / "11_epigenome_integration/analysis/12_grn_tf_methylation/literature_tf_master.csv"
INTEGRATED = BASE / "11_epigenome_integration/analysis/01_integration/integrated_methyl_expression_weighted.csv"
HC_SITES = BASE / "11_epigenome_integration/analysis/01_integration/high_confidence_sites_weighted.csv"
GFF = Path("/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/genomic.gff")

DESEQ_DIR = BASE / "04_deseq2/analysis/04_deseq2_260128_v1/results"
DESEQ_FILES = {
    "T2vsT1": DESEQ_DIR / "DESeq2_M145_2_vs_1.tsv",
    "T3vsT1": DESEQ_DIR / "DESeq2_M145_3_vs_1.tsv",
    "T3vsT2": DESEQ_DIR / "DESeq2_M145_3_vs_2.tsv",
}

# ============================================================
# A. LOAD DATA
# ============================================================
print("=" * 70)
print("LOADING DATA")
print("=" * 70)

# 1. TF master (37 TFs with corrected locus_tags)
tf_master = pd.read_csv(TF_MASTER)
print(f"TF master: {len(tf_master)} TFs loaded")

# 2. Integrated methylation-expression
integrated = pd.read_csv(INTEGRATED)
integrated.set_index("gene_id", inplace=True)
print(f"Integrated data: {len(integrated)} genes")

# 3. High-confidence methylation sites
hc_sites = pd.read_csv(HC_SITES)
print(f"High-confidence sites: {len(hc_sites)} sites")

# 4. DESeq2 results
deseq = {}
for comp, fpath in DESEQ_FILES.items():
    df = pd.read_csv(fpath, sep="\t")
    df.set_index("gene_id", inplace=True)
    deseq[comp] = df
    print(f"DESeq2 {comp}: {len(df)} genes")

# 5. Parse GFF for gene coordinates
print("\nParsing GFF for gene coordinates...")
gene_coords = {}
with open(GFF, "r") as f:
    for line in f:
        if line.startswith("#"):
            continue
        cols = line.strip().split("\t")
        if len(cols) < 9:
            continue
        if cols[2] != "gene":
            continue
        attrs = cols[8]
        # Extract locus_tag
        lt_match = re.search(r'locus_tag=([^;]+)', attrs)
        if lt_match:
            lt = lt_match.group(1)
            gene_coords[lt] = {
                "chrom": cols[0],
                "start": int(cols[3]),
                "end": int(cols[4]),
                "strand": cols[6],
            }
            # Also extract gene name if present
            name_match = re.search(r'Name=([^;]+)', attrs)
            if name_match:
                gene_coords[lt]["name"] = name_match.group(1)
            # Extract old_locus_tag if present
            old_lt_match = re.search(r'old_locus_tag=([^;]+)', attrs)
            if old_lt_match:
                gene_coords[lt]["old_locus_tag"] = old_lt_match.group(1)

print(f"GFF parsed: {len(gene_coords)} genes with coordinates")

# ============================================================
# CROSS-VALIDATE TF LOCUS_TAGS
# ============================================================
print("\n" + "=" * 70)
print("CROSS-VALIDATING TF LOCUS_TAGS")
print("=" * 70)

for _, row in tf_master.iterrows():
    lt = row["locus_tag"]
    name = row["name"]
    if lt in gene_coords:
        old_lt = gene_coords[lt].get("old_locus_tag", "N/A")
        gff_name = gene_coords[lt].get("name", "N/A")
        in_integrated = "YES" if lt in integrated.index else "NO"
        print(f"  {name:20s} {lt:15s} -> GFF_name={gff_name:15s} old_lt={old_lt:10s} in_integrated={in_integrated}")
    else:
        print(f"  {name:20s} {lt:15s} -> NOT FOUND IN GFF")

# ============================================================
# B. EXTRACT METHYLATION + EXPRESSION FOR EACH TF
# ============================================================
print("\n" + "=" * 70)
print("EXTRACTING METHYLATION + EXPRESSION FOR 37 TFs")
print("=" * 70)

LFC_THRESHOLD = 1.0   # |log2FC| >= 1 for significant expression change
PADJ_THRESHOLD = 0.05

results = []

for _, row in tf_master.iterrows():
    lt = row["locus_tag"]
    name = row["name"]
    tier = row["tier"]
    category = row["category"]
    bgc = row.get("bgc", "")

    rec = {
        "name": name,
        "locus_tag": lt,
        "tier": tier,
        "category": category,
        "bgc": bgc if pd.notna(bgc) else "",
    }

    # Get coordinates
    if lt in gene_coords:
        rec["strand"] = gene_coords[lt]["strand"]
        rec["start"] = gene_coords[lt]["start"]
        rec["end"] = gene_coords[lt]["end"]
    else:
        rec["strand"] = "?"
        rec["start"] = 0
        rec["end"] = 0

    # Get methylation data from integrated
    if lt in integrated.index:
        irow = integrated.loc[lt]
        for tp in ["T1", "T2", "T3"]:
            for mod in ["6mA", "4mC"]:
                rec[f"{mod}_{tp}_count"] = int(irow[f"{mod}_{tp}_count"])
                rec[f"{mod}_{tp}_freq"] = float(irow[f"{mod}_{tp}_mean_freq"])
        rec["has_methylation"] = True
    else:
        for tp in ["T1", "T2", "T3"]:
            for mod in ["6mA", "4mC"]:
                rec[f"{mod}_{tp}_count"] = 0
                rec[f"{mod}_{tp}_freq"] = 0.0
        rec["has_methylation"] = False

    # Get expression data from DESeq2
    for comp in ["T2vsT1", "T3vsT1", "T3vsT2"]:
        if lt in deseq[comp].index:
            rec[f"log2FC_{comp}"] = deseq[comp].loc[lt, "log2FoldChange"]
            rec[f"padj_{comp}"] = deseq[comp].loc[lt, "padj"]
        else:
            rec[f"log2FC_{comp}"] = np.nan
            rec[f"padj_{comp}"] = np.nan

    # Calculate total methylation per timepoint
    for tp in ["T1", "T2", "T3"]:
        rec[f"total_meth_{tp}"] = rec[f"6mA_{tp}_count"] + rec[f"4mC_{tp}_count"]

    results.append(rec)

df_tf = pd.DataFrame(results)

# ============================================================
# METHYLATION CHANGE CATEGORIES
# ============================================================
for comp, tp_from, tp_to in [("T2vsT1", "T1", "T2"), ("T3vsT1", "T1", "T3"), ("T3vsT2", "T2", "T3")]:
    for mod in ["6mA", "4mC"]:
        col_from_cnt = f"{mod}_{tp_from}_count"
        col_to_cnt = f"{mod}_{tp_to}_count"
        col_from_freq = f"{mod}_{tp_from}_freq"
        col_to_freq = f"{mod}_{tp_to}_freq"
        change_col = f"{mod}_change_{comp}"

        changes = []
        for _, r in df_tf.iterrows():
            cnt_from = r[col_from_cnt]
            cnt_to = r[col_to_cnt]
            freq_from = r[col_from_freq]
            freq_to = r[col_to_freq]

            if cnt_from == 0 and cnt_to == 0:
                changes.append("absent")
            elif cnt_from == 0 and cnt_to > 0:
                changes.append("gained")
            elif cnt_from > 0 and cnt_to == 0:
                changes.append("lost")
            elif cnt_to > cnt_from:
                changes.append("gained")
            elif cnt_to < cnt_from:
                changes.append("lost")
            else:
                # Same count - check frequency change
                if abs(freq_to - freq_from) > 10:  # >10% frequency change
                    changes.append("freq_change")
                else:
                    changes.append("stable")
        df_tf[change_col] = changes

# ============================================================
# EXPRESSION CHANGE CATEGORIES
# ============================================================
for comp in ["T2vsT1", "T3vsT1", "T3vsT2"]:
    expr_changes = []
    for _, r in df_tf.iterrows():
        lfc = r[f"log2FC_{comp}"]
        padj = r[f"padj_{comp}"]
        if pd.isna(lfc) or pd.isna(padj):
            expr_changes.append("no_data")
        elif padj < PADJ_THRESHOLD and lfc >= LFC_THRESHOLD:
            expr_changes.append("up")
        elif padj < PADJ_THRESHOLD and lfc <= -LFC_THRESHOLD:
            expr_changes.append("down")
        elif padj < PADJ_THRESHOLD:
            expr_changes.append("mild_change")
        else:
            expr_changes.append("stable")
    df_tf[f"expr_change_{comp}"] = expr_changes

# ============================================================
# COORDINATION ASSESSMENT
# ============================================================
def assess_coordination(meth_change, expr_change):
    """Determine if methylation and expression changes are coordinated."""
    if meth_change in ["absent", "stable"] and expr_change in ["stable", "no_data"]:
        return "no_change"
    if meth_change in ["absent", "stable"]:
        return "expr_only"
    if expr_change in ["stable", "no_data"]:
        return "meth_only"

    # Both change
    if meth_change == "gained" and expr_change == "up":
        return "concordant_up"
    elif meth_change == "gained" and expr_change == "down":
        return "discordant_gain_down"
    elif meth_change == "lost" and expr_change == "up":
        return "discordant_loss_up"
    elif meth_change == "lost" and expr_change == "down":
        return "concordant_down"
    elif meth_change == "freq_change":
        return "freq_change_with_expr"
    elif meth_change == "gained" and expr_change == "mild_change":
        return "meth_gain_mild_expr"
    elif meth_change == "lost" and expr_change == "mild_change":
        return "meth_loss_mild_expr"
    else:
        return "other"

for comp in ["T2vsT1", "T3vsT1", "T3vsT2"]:
    for mod in ["6mA", "4mC"]:
        coord_col = f"{mod}_coordination_{comp}"
        df_tf[coord_col] = df_tf.apply(
            lambda r: assess_coordination(r[f"{mod}_change_{comp}"], r[f"expr_change_{comp}"]),
            axis=1
        )
    # Combined coordination (either modification)
    combined_col = f"combined_coordination_{comp}"
    combined = []
    for _, r in df_tf.iterrows():
        c6 = r[f"6mA_coordination_{comp}"]
        c4 = r[f"4mC_coordination_{comp}"]
        if "concordant" in c6 or "concordant" in c4:
            combined.append("concordant")
        elif "discordant" in c6 or "discordant" in c4:
            combined.append("discordant")
        elif c6 == "no_change" and c4 == "no_change":
            combined.append("no_change")
        elif "meth_only" in c6 or "meth_only" in c4:
            combined.append("meth_only")
        elif "expr_only" in c6 or "expr_only" in c4:
            combined.append("expr_only")
        else:
            combined.append("partial")
    df_tf[combined_col] = combined

# ============================================================
# C. SARP UPSTREAM FOCUS
# ============================================================
print("\n" + "=" * 70)
print("SARP UPSTREAM REGULATORY HIERARCHY")
print("=" * 70)

sarp_hierarchy = {
    "Red": {
        "sarp": ("redD", "SC_RS31630"),
        "upstream": [
            ("redZ", "SC_RS31650"),
            ("afsR", "SC_RS24295"),
            ("hrdD", "SC_RS18125"),
            ("absA2", "SC_RS18245"),
        ],
    },
    "Act": {
        "sarp": ("actII-ORF4", "SC_RS27585"),
        "upstream": [
            ("afsR", "SC_RS24295"),
            ("hrdD", "SC_RS18125"),
            ("absA2", "SC_RS18245"),
            ("afsQ1", "SC_RS26695"),
        ],
    },
    "CDA": {
        "sarp": ("cdaR", "SC_RS18200"),
        "upstream": [
            ("absA2", "SC_RS18245"),
            ("afsR", "SC_RS24295"),
        ],
    },
    "CPK": {
        "sarp": ("cpkO/kasO", "SC_RS33650"),
        "upstream": [
            ("afsR", "SC_RS24295"),
            ("absA2", "SC_RS18245"),
        ],
    },
}

sarp_rows = []
for bgc_name, info in sarp_hierarchy.items():
    sarp_name, sarp_lt = info["sarp"]
    # SARP itself
    tf_row = df_tf[df_tf["locus_tag"] == sarp_lt]
    if len(tf_row) > 0:
        r = tf_row.iloc[0]
        sarp_rows.append({
            "BGC": bgc_name,
            "role": "SARP",
            "TF_name": sarp_name,
            "locus_tag": sarp_lt,
            "6mA_T1": r["6mA_T1_count"], "6mA_T2": r["6mA_T2_count"], "6mA_T3": r["6mA_T3_count"],
            "4mC_T1": r["4mC_T1_count"], "4mC_T2": r["4mC_T2_count"], "4mC_T3": r["4mC_T3_count"],
            "6mA_T1_freq": r["6mA_T1_freq"], "6mA_T2_freq": r["6mA_T2_freq"], "6mA_T3_freq": r["6mA_T3_freq"],
            "4mC_T1_freq": r["4mC_T1_freq"], "4mC_T2_freq": r["4mC_T2_freq"], "4mC_T3_freq": r["4mC_T3_freq"],
            "LFC_T2v1": r["log2FC_T2vsT1"], "padj_T2v1": r["padj_T2vsT1"],
            "LFC_T3v1": r["log2FC_T3vsT1"], "padj_T3v1": r["padj_T3vsT1"],
            "LFC_T3v2": r["log2FC_T3vsT2"], "padj_T3v2": r["padj_T3vsT2"],
            "6mA_change_T2v1": r["6mA_change_T2vsT1"],
            "6mA_change_T3v1": r["6mA_change_T3vsT1"],
            "4mC_change_T2v1": r["4mC_change_T2vsT1"],
            "4mC_change_T3v1": r["4mC_change_T3vsT1"],
            "expr_change_T2v1": r["expr_change_T2vsT1"],
            "expr_change_T3v1": r["expr_change_T3vsT1"],
        })

    # Upstream TFs
    for up_name, up_lt in info["upstream"]:
        tf_row = df_tf[df_tf["locus_tag"] == up_lt]
        if len(tf_row) > 0:
            r = tf_row.iloc[0]
            sarp_rows.append({
                "BGC": bgc_name,
                "role": "upstream",
                "TF_name": up_name,
                "locus_tag": up_lt,
                "6mA_T1": r["6mA_T1_count"], "6mA_T2": r["6mA_T2_count"], "6mA_T3": r["6mA_T3_count"],
                "4mC_T1": r["4mC_T1_count"], "4mC_T2": r["4mC_T2_count"], "4mC_T3": r["4mC_T3_count"],
                "6mA_T1_freq": r["6mA_T1_freq"], "6mA_T2_freq": r["6mA_T2_freq"], "6mA_T3_freq": r["6mA_T3_freq"],
                "4mC_T1_freq": r["4mC_T1_freq"], "4mC_T2_freq": r["4mC_T2_freq"], "4mC_T3_freq": r["4mC_T3_freq"],
                "LFC_T2v1": r["log2FC_T2vsT1"], "padj_T2v1": r["padj_T2vsT1"],
                "LFC_T3v1": r["log2FC_T3vsT1"], "padj_T3v1": r["padj_T3vsT1"],
                "LFC_T3v2": r["log2FC_T3vsT2"], "padj_T3v2": r["padj_T3vsT2"],
                "6mA_change_T2v1": r["6mA_change_T2vsT1"],
                "6mA_change_T3v1": r["6mA_change_T3vsT1"],
                "4mC_change_T2v1": r["4mC_change_T2vsT1"],
                "4mC_change_T3v1": r["4mC_change_T3vsT1"],
                "expr_change_T2v1": r["expr_change_T2vsT1"],
                "expr_change_T3v1": r["expr_change_T3vsT1"],
            })

df_sarp = pd.DataFrame(sarp_rows)

# Print SARP cascade summary
for bgc_name in ["Red", "Act", "CDA", "CPK"]:
    sub = df_sarp[df_sarp["BGC"] == bgc_name]
    print(f"\n--- {bgc_name} BGC ---")
    for _, r in sub.iterrows():
        role_str = f"[{r['role']:8s}]"
        meth_str = (f"6mA:{int(r['6mA_T1'])}/{int(r['6mA_T2'])}/{int(r['6mA_T3'])}  "
                    f"4mC:{int(r['4mC_T1'])}/{int(r['4mC_T2'])}/{int(r['4mC_T3'])}")
        lfc_str = f"LFC(T2v1)={r['LFC_T2v1']:+.2f}  LFC(T3v1)={r['LFC_T3v1']:+.2f}" if pd.notna(r["LFC_T2v1"]) else "LFC=N/A"
        print(f"  {role_str} {r['TF_name']:15s} ({r['locus_tag']})  {meth_str}  {lfc_str}")

# ============================================================
# D. CHECK AAGCCCG MOTIF IN PROMOTER
# ============================================================
print("\n" + "=" * 70)
print("AAGCCCG MOTIF CHECK IN PROMOTER REGIONS")
print("=" * 70)

MOTIF = "AAGCCCG"
MOTIF_RC = "CGGGCTT"  # reverse complement
PROMOTER_UP = 300
PROMOTER_DOWN = 50

# Load genome sequence for motif scanning
genome_fasta = Path("/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/GCF_000203835.1_ASM20383v1_genomic.fna")

genome_seq = None
if genome_fasta.exists():
    print(f"Loading genome from {genome_fasta.name}...")
    current_chr = None
    seqs = {}
    with open(genome_fasta) as f:
        for line in f:
            if line.startswith(">"):
                current_chr = line.strip().split()[0][1:]
                seqs[current_chr] = []
            else:
                seqs[current_chr].append(line.strip())
    for k in seqs:
        seqs[k] = "".join(seqs[k])
    genome_seq = seqs
    print(f"  Loaded {len(genome_seq)} chromosome(s)")
    for k, v in genome_seq.items():
        print(f"    {k}: {len(v):,} bp")
else:
    print(f"WARNING: Genome FASTA not found at {genome_fasta}")

motif_results = []
for _, row in df_tf.iterrows():
    lt = row["locus_tag"]
    name = row["name"]

    if lt not in gene_coords:
        motif_results.append({"name": name, "locus_tag": lt, "promoter_sites": 0,
                              "has_AAGCCCG": False, "motif_positions": ""})
        continue

    gc = gene_coords[lt]
    chrom = gc["chrom"]
    strand = gc["strand"]

    # Define promoter region based on strand
    if strand == "+":
        prom_start = max(1, gc["start"] - PROMOTER_UP)
        prom_end = gc["start"] + PROMOTER_DOWN
    else:
        prom_start = max(1, gc["end"] - PROMOTER_DOWN)
        prom_end = gc["end"] + PROMOTER_UP

    # Count HC methylation sites in promoter
    prom_sites = hc_sites[
        (hc_sites["chrom"] == chrom) &
        (hc_sites["position"] >= prom_start) &
        (hc_sites["position"] <= prom_end)
    ]

    # Check for AAGCCCG motif in promoter sequence
    has_motif = False
    motif_positions = []
    if genome_seq and chrom in genome_seq:
        seq = genome_seq[chrom]
        prom_seq = seq[prom_start - 1:prom_end]  # 0-based
        # Search both strands
        for i in range(len(prom_seq) - len(MOTIF) + 1):
            subseq = prom_seq[i:i + len(MOTIF)].upper()
            if subseq == MOTIF or subseq == MOTIF_RC:
                has_motif = True
                motif_positions.append(prom_start + i)

    n_unique_sites = len(prom_sites[["position", "mod_type"]].drop_duplicates()) if len(prom_sites) > 0 else 0

    motif_results.append({
        "name": name,
        "locus_tag": lt,
        "promoter_sites": n_unique_sites,
        "has_AAGCCCG": has_motif,
        "motif_positions": ";".join(map(str, motif_positions)),
        "prom_start": prom_start,
        "prom_end": prom_end,
        "strand": strand,
    })

df_motif = pd.DataFrame(motif_results)

# Merge motif info into main TF table
df_tf = df_tf.merge(df_motif[["locus_tag", "promoter_sites", "has_AAGCCCG", "motif_positions"]],
                     on="locus_tag", how="left")

print(f"\nTFs with AAGCCCG in promoter: {df_motif['has_AAGCCCG'].sum()}/{len(df_motif)}")
for _, r in df_motif[df_motif["has_AAGCCCG"]].iterrows():
    print(f"  {r['name']:20s} ({r['locus_tag']})  promoter_sites={r['promoter_sites']}  motif_pos={r['motif_positions']}")

print(f"\nTFs with HC methylation sites in promoter:")
for _, r in df_motif[df_motif["promoter_sites"] > 0].iterrows():
    print(f"  {r['name']:20s} ({r['locus_tag']})  sites={r['promoter_sites']}  AAGCCCG={r['has_AAGCCCG']}")

# ============================================================
# E. VISUALIZATIONS
# ============================================================
print("\n" + "=" * 70)
print("GENERATING VISUALIZATIONS")
print("=" * 70)

plt.rcParams.update({
    "font.size": 9,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.dpi": 150,
})

# --------------------------------------------------------
# FIGURE 1: Heatmap of 37 TFs x methylation + expression
# --------------------------------------------------------
print("  Figure 1: TF methylation-expression heatmap...")

# Prepare heatmap data
hm_cols = [
    "6mA_T1_count", "6mA_T2_count", "6mA_T3_count",
    "4mC_T1_count", "4mC_T2_count", "4mC_T3_count",
    "log2FC_T2vsT1", "log2FC_T3vsT1",
]
hm_labels = [
    "6mA\nT1", "6mA\nT2", "6mA\nT3",
    "4mC\nT1", "4mC\nT2", "4mC\nT3",
    "LFC\nT2v1", "LFC\nT3v1",
]

hm_data = df_tf[hm_cols].copy()
hm_data.index = df_tf["name"]

# Separate methylation counts and LFC for different normalization
meth_data = hm_data.iloc[:, :6].fillna(0)
lfc_data = hm_data.iloc[:, 6:].fillna(0)

# Sort by tier then category
sort_order = df_tf.sort_values(["tier", "category", "name"]).index
hm_data = hm_data.iloc[sort_order]
meth_data = meth_data.iloc[sort_order]
lfc_data = lfc_data.iloc[sort_order]

fig, axes = plt.subplots(1, 2, figsize=(12, 12), gridspec_kw={"width_ratios": [3, 1.2]})

# Left: methylation site counts
ax1 = axes[0]
sns.heatmap(
    meth_data, ax=ax1, cmap="YlOrRd", linewidths=0.5, linecolor="white",
    xticklabels=hm_labels[:6], yticklabels=True,
    cbar_kws={"label": "# sites", "shrink": 0.4},
    annot=True, fmt=".0f",
)
ax1.set_title("Methylation Site Counts per Timepoint", fontweight="bold")
ax1.set_ylabel("")

# Right: log2FC
ax2 = axes[1]
lfc_max = max(abs(lfc_data.min().min()), abs(lfc_data.max().max()))
if lfc_max == 0:
    lfc_max = 1
sns.heatmap(
    lfc_data, ax=ax2, cmap="RdBu_r", center=0,
    vmin=-lfc_max, vmax=lfc_max,
    linewidths=0.5, linecolor="white",
    xticklabels=hm_labels[6:], yticklabels=False,
    cbar_kws={"label": "log2FC", "shrink": 0.4},
    annot=True, fmt=".1f",
)
ax2.set_title("Expression Change", fontweight="bold")

plt.suptitle("H4: 37 TF Methylation-Expression Re-screen (Corrected Locus Tags)",
             fontsize=13, fontweight="bold", y=1.01)
plt.tight_layout()
fig.savefig(FIGURES / "fig1_TF_methylation_expression_heatmap.png", bbox_inches="tight", dpi=150)
plt.close()
print("    Saved fig1_TF_methylation_expression_heatmap.png")

# --------------------------------------------------------
# FIGURE 2: SARP cascade bar charts
# --------------------------------------------------------
print("  Figure 2: SARP cascade diagrams...")

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
bgc_order = ["Red", "Act", "CDA", "CPK"]
bgc_colors = {"Red": "#e41a1c", "Act": "#377eb8", "CDA": "#4daf4a", "CPK": "#984ea3"}

for idx, bgc_name in enumerate(bgc_order):
    ax = axes[idx // 2, idx % 2]
    sub = df_sarp[df_sarp["BGC"] == bgc_name].copy()
    sub = sub.reset_index(drop=True)

    n_tf = len(sub)
    x = np.arange(n_tf)
    width = 0.25

    # Plot LFC for T2v1 and T3v1
    lfc_t2 = sub["LFC_T2v1"].fillna(0).values
    lfc_t3 = sub["LFC_T3v1"].fillna(0).values

    bars1 = ax.bar(x - width / 2, lfc_t2, width, label="LFC T2vsT1",
                   color=bgc_colors[bgc_name], alpha=0.6, edgecolor="black", linewidth=0.5)
    bars2 = ax.bar(x + width / 2, lfc_t3, width, label="LFC T3vsT1",
                   color=bgc_colors[bgc_name], alpha=1.0, edgecolor="black", linewidth=0.5)

    # Annotate methylation above bars
    for i, (_, r) in enumerate(sub.iterrows()):
        total_6mA = r["6mA_T1"] + r["6mA_T2"] + r["6mA_T3"]
        total_4mC = r["4mC_T1"] + r["4mC_T2"] + r["4mC_T3"]
        y_max = max(abs(lfc_t2[i]), abs(lfc_t3[i]), 0.5)
        marker = ""
        if total_6mA > 0:
            marker += f"6mA:{int(r['6mA_T1'])}/{int(r['6mA_T2'])}/{int(r['6mA_T3'])}"
        if total_4mC > 0:
            if marker:
                marker += "\n"
            marker += f"4mC:{int(r['4mC_T1'])}/{int(r['4mC_T2'])}/{int(r['4mC_T3'])}"
        if not marker:
            marker = "no meth"
        ax.text(i, y_max + 0.3, marker, ha="center", va="bottom", fontsize=6.5,
                fontstyle="italic", color="darkred")

    ax.set_xticks(x)
    labels = []
    for _, r in sub.iterrows():
        role_marker = "*" if r["role"] == "SARP" else ""
        labels.append(f"{r['TF_name']}{role_marker}")
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("log2 Fold Change")
    ax.set_title(f"{bgc_name} BGC Regulatory Cascade", fontweight="bold")
    ax.axhline(y=0, color="gray", linestyle="--", linewidth=0.5)
    ax.legend(fontsize=7)

plt.suptitle("SARP Upstream Regulators: Expression + Methylation Status",
             fontsize=13, fontweight="bold")
plt.tight_layout()
fig.savefig(FIGURES / "fig2_SARP_cascade_methylation.png", bbox_inches="tight", dpi=150)
plt.close()
print("    Saved fig2_SARP_cascade_methylation.png")

# --------------------------------------------------------
# FIGURE 3: Coordination summary
# --------------------------------------------------------
print("  Figure 3: Coordination summary...")

fig, axes = plt.subplots(1, 3, figsize=(14, 5))
comps = ["T2vsT1", "T3vsT1", "T3vsT2"]
comp_labels = ["T2 vs T1", "T3 vs T1", "T3 vs T2"]

coord_categories = ["concordant", "discordant", "meth_only", "expr_only", "no_change", "partial"]
coord_colors = {
    "concordant": "#2ca02c",
    "discordant": "#d62728",
    "meth_only": "#ff7f0e",
    "expr_only": "#1f77b4",
    "no_change": "#7f7f7f",
    "partial": "#9467bd",
}

for i, (comp, label) in enumerate(zip(comps, comp_labels)):
    ax = axes[i]
    col = f"combined_coordination_{comp}"
    counts = df_tf[col].value_counts()

    cats = [c for c in coord_categories if c in counts.index]
    vals = [counts[c] for c in cats]
    colors = [coord_colors[c] for c in cats]

    bars = ax.barh(cats, vals, color=colors, edgecolor="black", linewidth=0.5)
    ax.set_xlabel("Number of TFs")
    ax.set_title(label, fontweight="bold")

    for bar, val in zip(bars, vals):
        ax.text(bar.get_width() + 0.3, bar.get_y() + bar.get_height() / 2,
                str(val), va="center", fontsize=9)

    ax.set_xlim(0, max(vals) + 3 if vals else 5)

plt.suptitle("Methylation-Expression Coordination Summary (37 TFs)",
             fontsize=13, fontweight="bold")
plt.tight_layout()
fig.savefig(FIGURES / "fig3_coordination_summary.png", bbox_inches="tight", dpi=150)
plt.close()
print("    Saved fig3_coordination_summary.png")

# --------------------------------------------------------
# FIGURE 4: AAGCCCG motif + promoter methylation summary
# --------------------------------------------------------
print("  Figure 4: Promoter methylation + AAGCCCG motif...")

fig, ax = plt.subplots(figsize=(10, 8))
df_plot = df_tf[["name", "promoter_sites", "has_AAGCCCG", "tier", "category"]].copy()
df_plot = df_plot.sort_values(["tier", "category", "name"])
df_plot = df_plot.reset_index(drop=True)

colors = []
for _, r in df_plot.iterrows():
    if r["has_AAGCCCG"]:
        colors.append("#d62728")  # red for motif present
    elif r["promoter_sites"] > 0:
        colors.append("#1f77b4")  # blue for sites but no motif
    else:
        colors.append("#cccccc")  # gray for no sites

ax.barh(range(len(df_plot)), df_plot["promoter_sites"], color=colors,
        edgecolor="black", linewidth=0.3)
ax.set_yticks(range(len(df_plot)))
ax.set_yticklabels(df_plot["name"], fontsize=8)
ax.set_xlabel("HC Methylation Sites in Promoter (-300 to +50 bp)")
ax.set_title("TF Promoter Methylation + AAGCCCG Motif Presence", fontweight="bold")
ax.invert_yaxis()

legend_elements = [
    mpatches.Patch(facecolor="#d62728", edgecolor="black", label="AAGCCCG present"),
    mpatches.Patch(facecolor="#1f77b4", edgecolor="black", label="Sites present (no motif)"),
    mpatches.Patch(facecolor="#cccccc", edgecolor="black", label="No HC sites"),
]
ax.legend(handles=legend_elements, loc="lower right")

plt.tight_layout()
fig.savefig(FIGURES / "fig4_promoter_AAGCCCG_motif.png", bbox_inches="tight", dpi=150)
plt.close()
print("    Saved fig4_promoter_AAGCCCG_motif.png")

# ============================================================
# F. SAVE OUTPUT TABLES
# ============================================================
print("\n" + "=" * 70)
print("SAVING OUTPUT TABLES")
print("=" * 70)

# Table 1: Full 37-TF table
out_cols = [
    "name", "locus_tag", "tier", "category", "bgc", "strand", "start", "end",
    "6mA_T1_count", "6mA_T1_freq", "6mA_T2_count", "6mA_T2_freq", "6mA_T3_count", "6mA_T3_freq",
    "4mC_T1_count", "4mC_T1_freq", "4mC_T2_count", "4mC_T2_freq", "4mC_T3_count", "4mC_T3_freq",
    "total_meth_T1", "total_meth_T2", "total_meth_T3",
    "log2FC_T2vsT1", "padj_T2vsT1", "log2FC_T3vsT1", "padj_T3vsT1", "log2FC_T3vsT2", "padj_T3vsT2",
    "6mA_change_T2vsT1", "6mA_change_T3vsT1", "6mA_change_T3vsT2",
    "4mC_change_T2vsT1", "4mC_change_T3vsT1", "4mC_change_T3vsT2",
    "expr_change_T2vsT1", "expr_change_T3vsT1", "expr_change_T3vsT2",
    "combined_coordination_T2vsT1", "combined_coordination_T3vsT1", "combined_coordination_T3vsT2",
    "promoter_sites", "has_AAGCCCG", "motif_positions",
]
df_tf[out_cols].to_csv(TABLES / "TF_methylation_expression_corrected.tsv", sep="\t", index=False)
print(f"  Saved TF_methylation_expression_corrected.tsv ({len(df_tf)} rows)")

# Table 2: Coordination summary
coord_cols = [
    "name", "locus_tag", "tier", "category", "bgc",
    "6mA_coordination_T2vsT1", "4mC_coordination_T2vsT1", "combined_coordination_T2vsT1",
    "6mA_coordination_T3vsT1", "4mC_coordination_T3vsT1", "combined_coordination_T3vsT1",
    "6mA_coordination_T3vsT2", "4mC_coordination_T3vsT2", "combined_coordination_T3vsT2",
]
df_tf[coord_cols].to_csv(TABLES / "TF_coordination_summary.tsv", sep="\t", index=False)
print(f"  Saved TF_coordination_summary.tsv ({len(df_tf)} rows)")

# Table 3: SARP upstream methylation
df_sarp.to_csv(TABLES / "SARP_upstream_methylation.tsv", sep="\t", index=False)
print(f"  Saved SARP_upstream_methylation.tsv ({len(df_sarp)} rows)")

# ============================================================
# KEY FINDINGS SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("KEY FINDINGS SUMMARY")
print("=" * 70)

# How many TFs have any methylation?
n_any_meth = df_tf[(df_tf["total_meth_T1"] > 0) | (df_tf["total_meth_T2"] > 0) | (df_tf["total_meth_T3"] > 0)].shape[0]
n_6mA = df_tf[(df_tf["6mA_T1_count"] > 0) | (df_tf["6mA_T2_count"] > 0) | (df_tf["6mA_T3_count"] > 0)].shape[0]
n_4mC = df_tf[(df_tf["4mC_T1_count"] > 0) | (df_tf["4mC_T2_count"] > 0) | (df_tf["4mC_T3_count"] > 0)].shape[0]
n_deg = df_tf[(df_tf["expr_change_T2vsT1"].isin(["up", "down"])) |
              (df_tf["expr_change_T3vsT1"].isin(["up", "down"]))].shape[0]

print(f"\n1. METHYLATION OVERVIEW:")
print(f"   TFs with any methylation:  {n_any_meth}/37 ({100*n_any_meth/37:.0f}%)")
print(f"   TFs with 6mA:              {n_6mA}/37")
print(f"   TFs with 4mC:              {n_4mC}/37")
print(f"   TFs that are DEGs:         {n_deg}/37")

print(f"\n2. COORDINATION (combined, any comparison):")
for comp in ["T2vsT1", "T3vsT1", "T3vsT2"]:
    col = f"combined_coordination_{comp}"
    conc = (df_tf[col] == "concordant").sum()
    disc = (df_tf[col] == "discordant").sum()
    monly = (df_tf[col] == "meth_only").sum()
    eonly = (df_tf[col] == "expr_only").sum()
    nochg = (df_tf[col] == "no_change").sum()
    print(f"   {comp}: concordant={conc}, discordant={disc}, meth_only={monly}, expr_only={eonly}, no_change={nochg}")

print(f"\n3. SARP REGULATORS:")
for bgc_name in ["Red", "Act", "CDA", "CPK"]:
    sub = df_sarp[df_sarp["BGC"] == bgc_name]
    sarp_row = sub[sub["role"] == "SARP"]
    if len(sarp_row) > 0:
        sr = sarp_row.iloc[0]
        meth_total = sr["6mA_T1"] + sr["6mA_T2"] + sr["6mA_T3"] + sr["4mC_T1"] + sr["4mC_T2"] + sr["4mC_T3"]
        lfc = f"T2v1={sr['LFC_T2v1']:+.2f}" if pd.notna(sr['LFC_T2v1']) else "T2v1=N/A"
        lfc3 = f"T3v1={sr['LFC_T3v1']:+.2f}" if pd.notna(sr['LFC_T3v1']) else "T3v1=N/A"
        print(f"   {bgc_name}: {sr['TF_name']} - total_meth_sites={int(meth_total)}  {lfc}  {lfc3}")
        # Upstream summary
        ups = sub[sub["role"] == "upstream"]
        for _, ur in ups.iterrows():
            um = ur["6mA_T1"] + ur["6mA_T2"] + ur["6mA_T3"] + ur["4mC_T1"] + ur["4mC_T2"] + ur["4mC_T3"]
            ulfc = f"T2v1={ur['LFC_T2v1']:+.2f}" if pd.notna(ur['LFC_T2v1']) else "N/A"
            print(f"       upstream: {ur['TF_name']:12s}  meth={int(um)}  {ulfc}")

print(f"\n4. AAGCCCG MOTIF:")
n_motif = df_tf["has_AAGCCCG"].sum()
print(f"   TFs with AAGCCCG in promoter: {n_motif}/37")
for _, r in df_tf[df_tf["has_AAGCCCG"] == True].iterrows():
    print(f"     {r['name']} ({r['locus_tag']}) - {r['motif_positions']}")

print(f"\n5. PROMOTER METHYLATION:")
n_prom_meth = (df_tf["promoter_sites"] > 0).sum()
print(f"   TFs with HC methylation in promoter: {n_prom_meth}/37")
top_prom = df_tf[df_tf["promoter_sites"] > 0].sort_values("promoter_sites", ascending=False).head(10)
for _, r in top_prom.iterrows():
    print(f"     {r['name']:20s}  {r['promoter_sites']} sites  AAGCCCG={r['has_AAGCCCG']}")

# Highlight methylated + differentially expressed TFs
print(f"\n6. TFs WITH BOTH METHYLATION AND DIFFERENTIAL EXPRESSION:")
for _, r in df_tf.iterrows():
    has_meth = (r["total_meth_T1"] + r["total_meth_T2"] + r["total_meth_T3"]) > 0
    is_deg = r["expr_change_T2vsT1"] in ["up", "down"] or r["expr_change_T3vsT1"] in ["up", "down"]
    if has_meth and is_deg:
        meth_str = f"6mA:{r['6mA_T1_count']}/{r['6mA_T2_count']}/{r['6mA_T3_count']}  4mC:{r['4mC_T1_count']}/{r['4mC_T2_count']}/{r['4mC_T3_count']}"
        expr_str = f"T2v1={r['expr_change_T2vsT1']}  T3v1={r['expr_change_T3vsT1']}"
        coord = r.get("combined_coordination_T3vsT1", "?")
        print(f"   {r['name']:20s} ({r['locus_tag']})  {meth_str}  {expr_str}  coord={coord}")

print("\n" + "=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)
