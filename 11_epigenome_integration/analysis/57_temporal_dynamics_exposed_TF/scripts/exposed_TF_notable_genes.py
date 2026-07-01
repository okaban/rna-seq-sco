#!/usr/bin/env python3
"""Flag notable Exposed TF genes by methylation/expression change magnitude.

Spearman rho across the n=57 cohort is small and non-significant
(see motif_stratified_stats.txt). Instead of a correlation claim, we
list individual genes that show large changes and tag them by category
for downstream literature cross-checking.

Inputs
------
- 51_exposed_regulators_characteristics/tables/exposed_regulators_full_table.tsv
    region, product (annotation), tf_family, log2FC_T2, padj_T2
- 57_temporal_dynamics_exposed_TF/tables/exposed_TF_methyl_motif_stratified.tsv
    per-motif deltas (% mod freq): GCCGGC_delta, AAGCCCG_delta, unassigned_6mA_delta
- 57_temporal_dynamics_exposed_TF/tables/exposed_TF_methyl_change_vs_LFC.tsv
    methyl_change (0-1): combined cross-motif aggregate per gene

Outputs
-------
- tables/exposed_TF_notable_list.tsv          (all 57 genes, with category flags)
- ../figures/exposed_TF_notable_highlight.{pdf,png}
- results/notable_genes_for_literature_check.txt
"""
from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
ANA = HERE.parent.parent  # 11_epigenome_integration/analysis
DIR_57 = ANA / "57_temporal_dynamics_exposed_TF"
DIR_51 = ANA / "51_exposed_regulators_characteristics"
FIG_DIR = ANA / "figures"

PATH_FULL = DIR_51 / "tables" / "exposed_regulators_full_table.tsv"
PATH_MOTIF = DIR_57 / "tables" / "exposed_TF_methyl_motif_stratified.tsv"
PATH_COMBINED = DIR_57 / "tables" / "exposed_TF_methyl_change_vs_LFC.tsv"

OUT_TABLE = DIR_57 / "tables" / "exposed_TF_notable_list.tsv"
OUT_FIG_PDF = FIG_DIR / "exposed_TF_notable_highlight.pdf"
OUT_FIG_PNG = FIG_DIR / "exposed_TF_notable_highlight.png"
OUT_TXT = DIR_57 / "results" / "notable_genes_for_literature_check.txt"

# ---------------------------------------------------------------------
# Thresholds
# ---------------------------------------------------------------------
THR_METHYL_A = 0.10   # |methyl_change| (0-1 scale, combined across motifs)
THR_LFC_A = 1.0       # |log2FC_T2|
THR_METHYL_B = 0.15   # large methylation change regardless of expression
THR_LFC_C = 2.0       # large expression change regardless of methylation

# ---------------------------------------------------------------------
# Keyword extraction for literature cross-check
# ---------------------------------------------------------------------
FUNCTION_KEYWORDS = [
    "sigma factor", "transcription factor", "regulator", "kinase",
    "phosphatase", "polymerase", "helicase", "nuclease", "methyltransferase",
    "transferase", "synthase", "synthetase", "dehydrogenase", "reductase",
    "oxidase", "ATPase", "permease", "transporter", "receptor",
    "sporulation", "secondary metabolite", "antibiotic", "biosynthesis",
    "stress response", "heat shock", "oxidative stress", "starvation",
    "cell wall", "membrane", "ribosome", "tRNA", "rRNA",
    "two-component", "response regulator", "histidine kinase",
    "morphogenesis", "development", "aerial", "mycelium",
    "actinorhodin", "undecylprodigiosin", "CDA", "yCPK",
    "iron", "zinc", "copper", "metal",
    "DNA-binding", "DNA repair", "restriction", "modification",
]

FAMILY_KEYWORDS = {
    "lysr": "LysR-type transcription factor (LTTR)",
    "tetr": "TetR-family repressor",
    "marr": "MarR-family regulator",
    "arac": "AraC-family regulator",
    "luxr": "LuxR-family regulator",
    "gntr": "GntR-family regulator",
    "iclr": "IclR-family regulator",
    "padr": "PadR-family regulator",
    "asnc": "AsnC/Lrp-family regulator",
    "merr": "MerR-family regulator",
    "deor": "DeoR-family regulator",
    "rok": "ROK-family regulator",
    "wblc": "WblC-like WhiB-family regulator",
    "whib": "WhiB-family regulator",
    "sigma factor": "sigma factor",
    "hth (other)": "helix-turn-helix DNA-binding",
    "two-component": "two-component response regulator",
    "tcs": "two-component response regulator",
    "tcs (rr)": "two-component response regulator",
    "sarp": "SARP (Streptomyces antibiotic regulatory protein)",
    "anti-sigma": "anti-sigma factor",
}


def extract_keywords(annotation: str | float, tf_family: str | float) -> list[str]:
    """Pull functional keywords for literature cross-check."""
    kws: list[str] = []
    fam_norm = "" if not isinstance(tf_family, str) else tf_family.strip().lower()
    if fam_norm:
        mapped = FAMILY_KEYWORDS.get(fam_norm)
        if mapped is None:
            for key, label in FAMILY_KEYWORDS.items():
                if key in fam_norm:
                    mapped = label
                    break
        if mapped:
            kws.append(mapped)
        else:
            kws.append(f"{tf_family} family")
    if isinstance(annotation, str) and annotation:
        text_lower = annotation.lower()
        for kw in FUNCTION_KEYWORDS:
            if kw.lower() in text_lower:
                kws.append(kw)
    kws.append("methylation")
    seen: set[str] = set()
    out: list[str] = []
    for k in kws:
        kl = k.lower()
        if kl in seen:
            continue
        seen.add(kl)
        out.append(k)
    return out


# ---------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------
df_full = pd.read_csv(PATH_FULL, sep="\t")
df_motif = pd.read_csv(PATH_MOTIF, sep="\t")
df_comb = pd.read_csv(PATH_COMBINED, sep="\t")

df_full_keep = df_full[
    [
        "locus_tag", "gene_name", "old_locus_tag", "product", "tf_family",
        "region", "log2FC_T2", "padj_T2", "temporal_pattern",
        "coordination_T2",
    ]
].rename(columns={
    "log2FC_T2": "log2FC",
    "padj_T2": "padj",
    "product": "annotation",
})

df_motif_keep = df_motif[
    ["locus_tag", "GCCGGC_delta", "AAGCCCG_delta", "unassigned_6mA_delta"]
].rename(columns={
    "GCCGGC_delta": "delta_GCCGGC_pct",
    "AAGCCCG_delta": "delta_AAGCCCG_pct",
    "unassigned_6mA_delta": "delta_6mA_pct",
})

df_comb_keep = df_comb[["locus_tag", "methyl_change"]].rename(
    columns={"methyl_change": "methyl_change_combined"}
)

df = (
    df_full_keep.merge(df_motif_keep, on="locus_tag", how="left")
    .merge(df_comb_keep, on="locus_tag", how="left")
)

# Sanity: 57 rows
assert len(df) == 57, f"Expected 57 exposed TFs, got {len(df)}"

# ---------------------------------------------------------------------
# Categorize
# ---------------------------------------------------------------------
methyl_abs = df["methyl_change_combined"].abs()
lfc_abs = df["log2FC"].abs()
sign_match = (
    ((df["methyl_change_combined"] > 0) & (df["log2FC"] > 0))
    | ((df["methyl_change_combined"] < 0) & (df["log2FC"] < 0))
)

df["category_A_codirected"] = (
    (methyl_abs > THR_METHYL_A) & (lfc_abs > THR_LFC_A) & sign_match
)
df["category_B_largeMethyl"] = methyl_abs > THR_METHYL_B
df["category_C_largeLFC"] = lfc_abs > THR_LFC_C


def label_categories(row: pd.Series) -> str:
    tags: list[str] = []
    if row["category_A_codirected"]:
        tags.append("A")
    if row["category_B_largeMethyl"]:
        tags.append("B")
    if row["category_C_largeLFC"]:
        tags.append("C")
    return ",".join(tags) if tags else "-"


df["category"] = df.apply(label_categories, axis=1)

# ---------------------------------------------------------------------
# Output table
# ---------------------------------------------------------------------
out_cols = [
    "locus_tag", "gene_name", "old_locus_tag", "tf_family", "annotation",
    "region", "log2FC", "padj", "methyl_change_combined",
    "delta_GCCGGC_pct", "delta_AAGCCCG_pct", "delta_6mA_pct",
    "temporal_pattern", "coordination_T2",
    "category_A_codirected", "category_B_largeMethyl", "category_C_largeLFC",
    "category",
]
df_out = df[out_cols].sort_values(
    by=["category_A_codirected", "category_B_largeMethyl", "category_C_largeLFC",
        "methyl_change_combined"],
    ascending=[False, False, False, False],
)

OUT_TABLE.parent.mkdir(parents=True, exist_ok=True)
df_out.to_csv(OUT_TABLE, sep="\t", index=False)
print(f"[ok] wrote {OUT_TABLE}")

n_a = int(df_out["category_A_codirected"].sum())
n_b = int(df_out["category_B_largeMethyl"].sum())
n_c = int(df_out["category_C_largeLFC"].sum())
print(f"     n_A (codirected) = {n_a}")
print(f"     n_B (|Δmethyl|>{THR_METHYL_B}) = {n_b}")
print(f"     n_C (|log2FC|>{THR_LFC_C}) = {n_c}")

# ---------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------
FIG_DIR.mkdir(parents=True, exist_ok=True)

fig, ax = plt.subplots(figsize=(10, 8))

ax.scatter(
    df["methyl_change_combined"], df["log2FC"],
    s=40, color="lightgrey", edgecolor="grey", linewidth=0.5,
    alpha=0.7, label="other Exposed TFs", zorder=2,
)

mask_C = df["category_C_largeLFC"] & ~df["category_A_codirected"]
mask_B = df["category_B_largeMethyl"] & ~df["category_A_codirected"]
mask_A = df["category_A_codirected"]

if mask_C.any():
    ax.scatter(
        df.loc[mask_C, "methyl_change_combined"], df.loc[mask_C, "log2FC"],
        s=110, marker="D", facecolor="mediumseagreen", edgecolor="black",
        linewidth=0.7, label=f"C: |log2FC|>{THR_LFC_C} (n={int(mask_C.sum())})",
        zorder=3,
    )
if mask_B.any():
    ax.scatter(
        df.loc[mask_B, "methyl_change_combined"], df.loc[mask_B, "log2FC"],
        s=120, marker="^", facecolor="royalblue", edgecolor="black",
        linewidth=0.7, label=f"B: |Δmethyl|>{THR_METHYL_B} (n={int(mask_B.sum())})",
        zorder=4,
    )
if mask_A.any():
    ax.scatter(
        df.loc[mask_A, "methyl_change_combined"], df.loc[mask_A, "log2FC"],
        s=180, marker="*", facecolor="crimson", edgecolor="black",
        linewidth=0.7,
        label=f"A: codirected, |Δmethyl|>{THR_METHYL_A} & |log2FC|>{THR_LFC_A} (n={int(mask_A.sum())})",
        zorder=5,
    )

# Label notable genes (categories A and B)
labelled = df[mask_A | mask_B].copy()
for _, row in labelled.iterrows():
    label = row["gene_name"] if isinstance(row["gene_name"], str) and row["gene_name"] else row["locus_tag"]
    ax.annotate(
        label,
        xy=(row["methyl_change_combined"], row["log2FC"]),
        xytext=(5, 5), textcoords="offset points",
        fontsize=8, color="black",
        path_effects=None,
    )

ax.axhline(0, color="black", linewidth=0.6, linestyle="--", alpha=0.6)
ax.axvline(0, color="black", linewidth=0.6, linestyle="--", alpha=0.6)
for thr, c in [(THR_METHYL_B, "royalblue"), (-THR_METHYL_B, "royalblue")]:
    ax.axvline(thr, color=c, linewidth=0.5, linestyle=":", alpha=0.5)
for thr, c in [(THR_LFC_C, "mediumseagreen"), (-THR_LFC_C, "mediumseagreen")]:
    ax.axhline(thr, color=c, linewidth=0.5, linestyle=":", alpha=0.5)

ax.set_xlabel("Δ methylation (combined across motifs, T2-T1, fraction units)")
ax.set_ylabel("log2 fold change (T2 vs T1)")
ax.set_title(
    f"Exposed TFs (n=57): notable genes by Δmethylation × Δexpression\n"
    f"individual gene listing (correlation across cohort is NS)",
    fontsize=11,
)
ax.legend(loc="upper left", fontsize=9, frameon=True)
ax.grid(True, linewidth=0.3, alpha=0.4)

fig.tight_layout()
fig.savefig(OUT_FIG_PDF)
fig.savefig(OUT_FIG_PNG, dpi=200)
plt.close(fig)
print(f"[ok] wrote {OUT_FIG_PDF.name} and .png")

# ---------------------------------------------------------------------
# Literature cross-check text file
# ---------------------------------------------------------------------
OUT_TXT.parent.mkdir(parents=True, exist_ok=True)

cat_a = df_out[df_out["category_A_codirected"]].copy()
cat_b_only = df_out[
    df_out["category_B_largeMethyl"] & ~df_out["category_A_codirected"]
].copy()
cat_c_only = df_out[
    df_out["category_C_largeLFC"] & ~df_out["category_A_codirected"]
].copy()


def fmt_block(idx: int, row: pd.Series) -> str:
    has_name = isinstance(row["gene_name"], str) and row["gene_name"].strip()
    has_old = isinstance(row["old_locus_tag"], str) and row["old_locus_tag"].strip() and row["old_locus_tag"] != "nan"
    parts: list[str] = [str(row["locus_tag"])]
    if has_name:
        parts.append(f"({row['gene_name']})")
    if has_old:
        parts.append(f"[{row['old_locus_tag']}]")
    header = " ".join(parts)
    annot = row["annotation"] if isinstance(row["annotation"], str) else "(no annotation)"
    keys = extract_keywords(annot, row["tf_family"])
    deltas: list[str] = []
    for col, label in [
        ("delta_GCCGGC_pct", "GCCGGC"),
        ("delta_AAGCCCG_pct", "AAGCCCG"),
        ("delta_6mA_pct", "6mA"),
    ]:
        val = row[col]
        if pd.notna(val) and val != 0:
            deltas.append(f"{label}: {val:+.1f}%")
    delta_str = ", ".join(deltas) if deltas else "no per-motif delta"
    region = row["region"] if isinstance(row["region"], str) else "?"
    family = row["tf_family"] if isinstance(row["tf_family"], str) else "?"
    return (
        f"{idx}. {header} | region={region} | family={family}\n"
        f"   log2FC = {row['log2FC']:+.2f} | padj = {row['padj']:.2e} | "
        f"Δmethyl(combined) = {row['methyl_change_combined']:+.3f}\n"
        f"   per-motif: {delta_str}\n"
        f"   annotation: \"{annot}\"\n"
        f"   → 照合キーワード: [{', '.join(keys)}]\n"
    )


lines: list[str] = []
lines.append("Exposed TF (n=57): notable gene list for cross-organism literature check")
lines.append("=" * 78)
lines.append("Generated by exposed_TF_notable_genes.py")
lines.append("")
lines.append("Thresholds:")
lines.append(f"  A (codirected) : |Δmethyl_combined| > {THR_METHYL_A} AND |log2FC| > {THR_LFC_A}")
lines.append(f"                   AND sign(Δmethyl) == sign(log2FC)")
lines.append(f"  B (large Δmethyl): |Δmethyl_combined| > {THR_METHYL_B}")
lines.append(f"  C (large LFC)    : |log2FC| > {THR_LFC_C}")
lines.append("")
lines.append("Note: Spearman rho across the cohort is small and non-significant")
lines.append("      (ρ=0.11-0.26 by motif, all p>0.24). This file is a per-gene")
lines.append("      listing for downstream literature cross-check, not a")
lines.append("      correlation claim.")
lines.append("")

lines.append("=== カテゴリA: メチル化変化 × 発現変化 同方向 ===")
lines.append(f"({len(cat_a)} genes)")
lines.append("")
if cat_a.empty:
    lines.append("(none)")
    lines.append("")
else:
    for i, (_, row) in enumerate(cat_a.iterrows(), start=1):
        lines.append(fmt_block(i, row))

lines.append("=== カテゴリB only: 大きなメチル化変化（Aに含まれないもの） ===")
lines.append(f"({len(cat_b_only)} genes)")
lines.append("")
if cat_b_only.empty:
    lines.append("(none)")
    lines.append("")
else:
    for i, (_, row) in enumerate(cat_b_only.iterrows(), start=1):
        lines.append(fmt_block(i, row))

lines.append("=== カテゴリC only: 大きな発現変化（Aに含まれないもの） ===")
lines.append(f"({len(cat_c_only)} genes)")
lines.append("")
if cat_c_only.empty:
    lines.append("(none)")
    lines.append("")
else:
    for i, (_, row) in enumerate(cat_c_only.iterrows(), start=1):
        lines.append(fmt_block(i, row))

OUT_TXT.write_text("\n".join(lines), encoding="utf-8")
print(f"[ok] wrote {OUT_TXT}")

print("")
print("=== summary ===")
print(f"category A (codirected, large both)     : {len(cat_a)}")
print(f"category B only (large Δmethyl, not A)  : {len(cat_b_only)}")
print(f"category C only (large LFC, not A)      : {len(cat_c_only)}")
n_any = int((df_out["category"] != "-").sum())
print(f"any category                            : {n_any} / 57")
