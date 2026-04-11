#!/usr/bin/env python3
"""
run_annotation_M145.py
Build master annotation table for S. coelicolor A3(2) M145 RNA-seq.
Integrates GFF annotation, DESeq2 results, BGC definitions, and regulator tags.
"""

import sys
import re
import datetime
from pathlib import Path

import pandas as pd

# ── Paths ──────────────────────────────────────────────────────────
ANNOT_RUN_DIR = Path("/Users/okaban/bioinfo/rna-seq/05_annotation/analysis/05_annotation_260128_v1")
GFF_FILE = Path("/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/genomic.gff")
DESEQ_DIR = Path("/Users/okaban/bioinfo/rna-seq/04_deseq2/analysis/04_deseq2_260128_v1/results")
NORM_COUNTS = DESEQ_DIR / "normalized_counts_M145.tsv"
DE_FILES = {
    "2_vs_1": DESEQ_DIR / "DESeq2_M145_2_vs_1.tsv",
    "3_vs_1": DESEQ_DIR / "DESeq2_M145_3_vs_1.tsv",
    "3_vs_2": DESEQ_DIR / "DESeq2_M145_3_vs_2.tsv",
}

TABLES = ANNOT_RUN_DIR / "tables"
LOG_FILE = ANNOT_RUN_DIR / "logs" / "annotation_pipeline.log"

# ── Logging ────────────────────────────────────────────────────────
log_fh = open(LOG_FILE, "w")

def log(msg):
    line = msg
    print(line)
    log_fh.write(line + "\n")
    log_fh.flush()

log(f"=== 05_annotation started at {datetime.datetime.now():%Y-%m-%d %H:%M:%S} ===")
log(f"ANNOT_RUN_DIR: {ANNOT_RUN_DIR}")
log(f"GFF: {GFF_FILE}")

# ══════════════════════════════════════════════════════════════════
# 1. Parse GFF → gene_annotation_basic.tsv
# ══════════════════════════════════════════════════════════════════
log("\n=== Step 1: Parsing GFF ===")

def parse_attr(attr_str, key):
    """Extract value for a key from GFF attributes string."""
    m = re.search(rf'{key}=([^;]+)', attr_str)
    return m.group(1) if m else None

gene_rows = []
cds_info = {}  # locus_tag -> {product, protein_id, ontology_term}

with open(GFF_FILE) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue
        contig, source, ftype, start, end, score, strand, phase, attrs = parts

        if ftype in ("gene", "pseudogene"):
            lt = parse_attr(attrs, "locus_tag")
            if lt is None:
                continue
            gene_rows.append({
                "gene_id": lt,
                "old_locus_tag": parse_attr(attrs, "old_locus_tag"),
                "gene_name": parse_attr(attrs, "gene"),
                "gene_biotype": parse_attr(attrs, "gene_biotype"),
                "contig": contig,
                "start": int(start),
                "end": int(end),
                "strand": strand,
            })
        elif ftype == "CDS":
            lt = parse_attr(attrs, "locus_tag")
            if lt is None:
                continue
            product = parse_attr(attrs, "product")
            if product:
                product = product.replace("%2C", ",").replace("%3B", ";").replace("%3A", ":").replace("%25", "%")
            pid = parse_attr(attrs, "protein_id")
            ont = parse_attr(attrs, "Ontology_term")
            # Keep first CDS per locus_tag (avoid duplicates)
            if lt not in cds_info:
                cds_info[lt] = {"product": product, "protein_id": pid, "ontology_term": ont}

gene_df = pd.DataFrame(gene_rows)
log(f"Parsed {len(gene_df)} gene/pseudogene records from GFF")
log(f"Parsed {len(cds_info)} CDS records with product info")

# Merge product / protein_id / ontology_term from CDS
cds_df = pd.DataFrame.from_dict(cds_info, orient="index")
cds_df.index.name = "gene_id"
cds_df = cds_df.reset_index()
gene_df = gene_df.merge(cds_df, on="gene_id", how="left")

# Save
out_basic = TABLES / "gene_annotation_basic.tsv"
gene_df.to_csv(out_basic, sep="\t", index=False)
log(f"Saved gene_annotation_basic.tsv ({gene_df.shape[0]} rows, {gene_df.shape[1]} cols)")

# ══════════════════════════════════════════════════════════════════
# 2. Merge DESeq2 + normalized counts → gene_master_DESeq2.tsv
# ══════════════════════════════════════════════════════════════════
log("\n=== Step 2: Merging DESeq2 results ===")

master = gene_df.copy()

for contrast, fpath in DE_FILES.items():
    de = pd.read_csv(fpath, sep="\t")
    de = de.rename(columns={"gene_id": "gene_id"})
    suffix = f"_{contrast}"
    cols_rename = {}
    for c in de.columns:
        if c == "gene_id":
            continue
        cols_rename[c] = f"{c}{suffix}"
    de = de.rename(columns=cols_rename)
    # Keep baseMean only from first contrast to avoid duplicates
    if contrast != "2_vs_1":
        bm_col = f"baseMean{suffix}"
        if bm_col in de.columns:
            de = de.drop(columns=[bm_col])
    master = master.merge(de, on="gene_id", how="left")
    log(f"  Merged {contrast}: {de.shape[0]} genes")

# Rename baseMean from first contrast
if "baseMean_2_vs_1" in master.columns:
    master = master.rename(columns={"baseMean_2_vs_1": "baseMean"})

# Merge normalized counts
norm = pd.read_csv(NORM_COUNTS, sep="\t")
# Rename columns to norm_xxx
norm_cols = {c: f"norm_{c}" for c in norm.columns if c != "gene_id"}
norm = norm.rename(columns=norm_cols)
master = master.merge(norm, on="gene_id", how="left")
log(f"  Merged normalized counts: {norm.shape[0]} genes, {len(norm_cols)} sample columns")

out_master = TABLES / "gene_master_DESeq2.tsv"
master.to_csv(out_master, sep="\t", index=False)
log(f"Saved gene_master_DESeq2.tsv ({master.shape[0]} rows, {master.shape[1]} cols)")

# ══════════════════════════════════════════════════════════════════
# 3. Build BGC_definition_manual.tsv (4 major BGCs)
# ══════════════════════════════════════════════════════════════════
log("\n=== Step 3: Building BGC_definition_manual.tsv ===")

# Define BGC ranges by old_locus_tag (SCO number)
bgc_ranges = {
    "act": {"sco_start": 5071, "sco_end": 5092},
    "red": {"sco_start": 5877, "sco_end": 5898},
    "cda": {"sco_start": 3210, "sco_end": 3249},
    "cpk": {"sco_start": 6273, "sco_end": 6288},
}

# Known gene roles from literature
known_roles = {
    # act
    "SCO5082": ("actII-orf4", "regulator"),   # actII, SARP
    "SCO5070": (None, "biosynthesis"),         # extended region
    "SCO5080": ("actVA", "biosynthesis"),      # actVA oxygenase
    "SCO5092": ("actVB", "biosynthesis"),      # actVB
    "SCO5085": (None, "transport"),            # actA, export
    "SCO5086": (None, "transport"),            # actB, export
    "SCO5087": (None, "resistance"),           # actR / actIII
    # red
    "SCO5877": ("redD", "regulator"),          # SARP
    "SCO5881": ("redZ", "regulator"),          # response regulator
    "SCO5879": ("redW", "biosynthesis"),
    "SCO5891": ("redM", "biosynthesis"),
    "SCO5898": (None, "transport"),            # redH-like
    # cda
    "SCO3225": ("absA1", "regulator"),         # two-component sensor kinase
    "SCO3226": ("absA2", "regulator"),         # response regulator
    "SCO3229": ("hppD", "biosynthesis"),
    "SCO3234": ("hasP", "biosynthesis"),
    "SCO3236": ("asnO", "biosynthesis"),
    "SCO3238": ("eboE", "biosynthesis"),
    # cpk
    "SCO6282": ("cpkO", "regulator"),          # SARP-like / kasO
    "SCO6286": ("scbR2", "regulator"),         # gamma-butyrolactone receptor
}

# Build from gene_annotation_basic
bgc_rows = []
for bgc_name, rng in bgc_ranges.items():
    sco_start = rng["sco_start"]
    sco_end = rng["sco_end"]
    for sco_num in range(sco_start, sco_end + 1):
        sco_tag = f"SCO{sco_num}"
        # Find matching gene
        match = gene_df[gene_df["old_locus_tag"] == sco_tag]
        if match.empty:
            # Some SCO numbers may not exist in GFF (gaps)
            log(f"  WARNING: {sco_tag} not found in GFF (may not exist in this annotation)")
            continue
        row = match.iloc[0]
        gene_id = row["gene_id"]
        gene_name_gff = row["gene_name"] if pd.notna(row["gene_name"]) else None

        # Check known roles
        if sco_tag in known_roles:
            known_name, role = known_roles[sco_tag]
            gene_name = known_name if known_name else gene_name_gff
        else:
            gene_name = gene_name_gff
            role = "biosynthesis"

        bgc_rows.append({
            "bgc_name": bgc_name,
            "gene_id": gene_id,
            "old_locus_tag": sco_tag,
            "gene_name": gene_name if gene_name else "NA",
            "role": role,
        })

bgc_df = pd.DataFrame(bgc_rows)
out_bgc = TABLES / "BGC_definition_manual.tsv"
bgc_df.to_csv(out_bgc, sep="\t", index=False)

for bgc_name in bgc_ranges:
    n = len(bgc_df[bgc_df["bgc_name"] == bgc_name])
    log(f"  {bgc_name}: {n} genes")
log(f"Saved BGC_definition_manual.tsv ({len(bgc_df)} total genes)")

# ══════════════════════════════════════════════════════════════════
# 4. Merge BGC info → gene_master_with_BGC.tsv
# ══════════════════════════════════════════════════════════════════
log("\n=== Step 4: Merging BGC info into master table ===")

bgc_merge = bgc_df[["gene_id", "bgc_name", "role"]].rename(columns={"role": "bgc_role"})
master_bgc = master.merge(bgc_merge, on="gene_id", how="left")

n_bgc = master_bgc["bgc_name"].notna().sum()
log(f"Genes with BGC annotation: {n_bgc} / {len(master_bgc)}")

out_master_bgc = TABLES / "gene_master_with_BGC.tsv"
master_bgc.to_csv(out_master_bgc, sep="\t", index=False)
log(f"Saved gene_master_with_BGC.tsv ({master_bgc.shape[0]} rows, {master_bgc.shape[1]} cols)")

# ══════════════════════════════════════════════════════════════════
# 5. BGC_major4_DE_summary.tsv
# ══════════════════════════════════════════════════════════════════
log("\n=== Step 5: Building BGC_major4_DE_summary.tsv ===")

bgc_genes = bgc_df["gene_id"].tolist()
bgc_summary = master_bgc[master_bgc["gene_id"].isin(bgc_genes)].copy()

# Select key columns
summary_cols = [
    "bgc_name", "gene_id", "old_locus_tag", "gene_name", "product", "bgc_role",
    "baseMean",
    "log2FoldChange_2_vs_1", "padj_2_vs_1",
    "log2FoldChange_3_vs_1", "padj_3_vs_1",
    "log2FoldChange_3_vs_2", "padj_3_vs_2",
]
# Add norm count columns
norm_count_cols = [c for c in master_bgc.columns if c.startswith("norm_")]
summary_cols_exist = [c for c in summary_cols if c in bgc_summary.columns] + norm_count_cols

bgc_summary = bgc_summary[summary_cols_exist].sort_values(["bgc_name", "gene_id"])

out_bgc_summary = TABLES / "BGC_major4_DE_summary.tsv"
bgc_summary.to_csv(out_bgc_summary, sep="\t", index=False)

# Print per-BGC DE counts
log("\nDE summary per BGC (padj < 0.05 & |log2FC| > 1, contrast: 3_vs_1):")
for bgc_name in ["act", "red", "cda", "cpk"]:
    sub = bgc_summary[bgc_summary["bgc_name"] == bgc_name]
    total = len(sub)
    has_padj = sub["padj_3_vs_1"].notna()
    sig_up = ((sub["padj_3_vs_1"] < 0.05) & (sub["log2FoldChange_3_vs_1"] > 1)).sum()
    sig_down = ((sub["padj_3_vs_1"] < 0.05) & (sub["log2FoldChange_3_vs_1"] < -1)).sum()
    sig_any = ((sub["padj_3_vs_1"] < 0.05)).sum()
    log(f"  {bgc_name}: {total} genes, {sig_up} up, {sig_down} down, {sig_any} sig total (padj<0.05)")

log(f"\nSaved BGC_major4_DE_summary.tsv ({len(bgc_summary)} rows)")

# ══════════════════════════════════════════════════════════════════
# 6. Regulator tagging → gene_master_with_BGC_regulators.tsv
# ══════════════════════════════════════════════════════════════════
log("\n=== Step 6: Tagging regulators ===")

# 6a. Pattern-based regulator detection from product annotation
regulator_patterns = [
    r"transcriptional regulator",
    r"SARP family",
    r"LuxR",
    r"AraC",
    r"TetR",
    r"LysR",
    r"MarR",
    r"GntR",
    r"IclR",
    r"MerR",
    r"two-component",
    r"sensor histidine kinase",
    r"response regulator",
    r"sigma factor",
    r"sigma-70",
    r"anti-sigma",
    r"butyrolactone",
    r"regulatory protein",
    r"DNA-binding.*regulator",
    r"repressor",
    r"activator",
    r"WhiB",
]
pattern = "|".join(regulator_patterns)

master_reg = master_bgc.copy()
master_reg["is_regulator"] = False
has_product = master_reg["product"].notna()
master_reg.loc[has_product, "is_regulator"] = master_reg.loc[has_product, "product"].str.contains(
    pattern, case=False, regex=True, na=False
)

n_auto = master_reg["is_regulator"].sum()
log(f"Auto-detected regulators from product annotation: {n_auto}")

# 6b. Manual regulator tags
manual_regulators = {
    "SC_RS27570": {"regulator_name": "actII-orf4", "regulator_type": "SARP (act cluster-situated)"},
    "SC_RS31630": {"regulator_name": "redD", "regulator_type": "SARP (red cluster-situated)"},
    "SC_RS31650": {"regulator_name": "redZ", "regulator_type": "response regulator"},
    "SC_RS25560": {"regulator_name": "atrA", "regulator_type": "TetR family"},
    "SC_RS33575": {"regulator_name": "scbR", "regulator_type": "gamma-butyrolactone receptor"},
    "SC_RS33680": {"regulator_name": "scbR2", "regulator_type": "gamma-butyrolactone receptor"},
    "SC_RS33660": {"regulator_name": "cpkO/kasO", "regulator_type": "SARP-like (cpk cluster-situated)"},
    "SC_RS18240": {"regulator_name": "absA1", "regulator_type": "two-component sensor kinase (cda cluster-situated)"},
    "SC_RS18245": {"regulator_name": "absA2", "regulator_type": "response regulator (cda cluster-situated)"},
}

master_reg["regulator_name"] = pd.NA
master_reg["regulator_type"] = pd.NA

for gid, info in manual_regulators.items():
    mask = master_reg["gene_id"] == gid
    if mask.any():
        master_reg.loc[mask, "is_regulator"] = True
        master_reg.loc[mask, "regulator_name"] = info["regulator_name"]
        master_reg.loc[mask, "regulator_type"] = info["regulator_type"]

n_manual = sum(1 for gid in manual_regulators if (master_reg["gene_id"] == gid).any())
n_total_reg = master_reg["is_regulator"].sum()
log(f"Manually tagged regulators: {n_manual}")
log(f"Total regulators (auto + manual): {n_total_reg}")

out_reg = TABLES / "gene_master_with_BGC_regulators.tsv"
master_reg.to_csv(out_reg, sep="\t", index=False)
log(f"Saved gene_master_with_BGC_regulators.tsv ({master_reg.shape[0]} rows, {master_reg.shape[1]} cols)")

# ══════════════════════════════════════════════════════════════════
# 7. Summary statistics for report
# ══════════════════════════════════════════════════════════════════
log("\n=== Step 7: Summary statistics ===")

# Per-BGC detailed summary for each contrast
for bgc_name in ["act", "red", "cda", "cpk"]:
    sub = bgc_summary[bgc_summary["bgc_name"] == bgc_name]
    total = len(sub)
    log(f"\n--- {bgc_name} cluster ({total} genes) ---")

    for contrast in ["2_vs_1", "3_vs_1", "3_vs_2"]:
        lfc_col = f"log2FoldChange_{contrast}"
        padj_col = f"padj_{contrast}"
        if lfc_col not in sub.columns:
            continue
        sig = sub[sub[padj_col] < 0.05]
        up = sig[sig[lfc_col] > 1]
        down = sig[sig[lfc_col] < -1]
        log(f"  {contrast}: {len(sig)}/{total} sig (padj<0.05), {len(up)} up (LFC>1), {len(down)} down (LFC<-1)")

    # Top genes by |LFC| in 3_vs_1
    lfc_col = "log2FoldChange_3_vs_1"
    padj_col = "padj_3_vs_1"
    top = sub[sub[padj_col].notna()].copy()
    top["abs_lfc"] = top[lfc_col].abs()
    top = top.sort_values("abs_lfc", ascending=False).head(5)
    log(f"  Top 5 by |LFC| (3_vs_1):")
    for _, r in top.iterrows():
        name = r["gene_name"] if pd.notna(r.get("gene_name")) and r["gene_name"] != "NA" else r["gene_id"]
        log(f"    {r['gene_id']} ({r['old_locus_tag']}, {name}): LFC={r[lfc_col]:.2f}, padj={r[padj_col]:.2e}")

# Regulator summary for BGC genes
log("\n--- Regulators within major 4 BGCs ---")
bgc_regs = master_reg[(master_reg["bgc_name"].notna()) & (master_reg["is_regulator"] == True)]
for _, r in bgc_regs.iterrows():
    rname = r["regulator_name"] if pd.notna(r["regulator_name"]) else r["product"]
    lfc = r.get("log2FoldChange_3_vs_1")
    padj = r.get("padj_3_vs_1")
    lfc_str = f"LFC={lfc:.2f}" if pd.notna(lfc) else "LFC=NA"
    padj_str = f"padj={padj:.2e}" if pd.notna(padj) else "padj=NA"
    log(f"  {r['bgc_name']}: {r['gene_id']} ({r['old_locus_tag']}, {rname}) {lfc_str}, {padj_str}")

log(f"\n=== 05_annotation COMPLETED at {datetime.datetime.now():%Y-%m-%d %H:%M:%S} ===")
log_fh.close()
