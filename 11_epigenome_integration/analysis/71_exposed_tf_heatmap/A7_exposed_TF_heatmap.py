"""
A-7: 57 Exposed Regulatory Genes - Expression Heatmap (Figure 5b)
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
from pathlib import Path

# ──────────────────────────────────────────────
# Paths
# ──────────────────────────────────────────────
BASE = Path.home() / "bioinfo/rna-seq"
GENE_LIST = BASE / "11_epigenome_integration/analysis/51_exposed_regulators_characteristics/tables/exposed_regulators_full_table.tsv"
NORM_COUNTS = BASE / "04_deseq2/analysis/04_deseq2_260128_v1/results/normalized_counts_M145.tsv"
DEG_T2 = BASE / "04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_2_vs_1.tsv"
OUT_DIR = Path(__file__).parent

# ──────────────────────────────────────────────
# Load data
# ──────────────────────────────────────────────
genes = pd.read_csv(GENE_LIST, sep="\t")
print(f"Exposed regulatory genes: {len(genes)}")

counts_raw = pd.read_csv(NORM_COUNTS, sep="\t", index_col=0)
print(f"Normalized counts: {counts_raw.shape}")

deg_t2 = pd.read_csv(DEG_T2, sep="\t", index_col=0)

# ──────────────────────────────────────────────
# Filter counts to exposed genes
# ──────────────────────────────────────────────
gene_ids = genes["locus_tag"].tolist()
counts = counts_raw.loc[counts_raw.index.isin(gene_ids)].copy()
print(f"Genes found in counts: {len(counts)}")

missing = [g for g in gene_ids if g not in counts_raw.index]
if missing:
    print(f"Missing from counts: {missing}")

# ──────────────────────────────────────────────
# Column ordering: T1 x3, T2 x3, T3 x3
# ──────────────────────────────────────────────
sample_order = [
    "M145_1_1", "M145_1_2", "M145_1_3",   # T1
    "M145_2_1", "M145_2_3", "M145_2_4",   # T2
    "M145_3_2", "M145_3_3", "M145_3_4",   # T3
]
counts = counts[sample_order]

# ──────────────────────────────────────────────
# Z-score per gene (row-wise)
# ──────────────────────────────────────────────
# Add pseudocount before log to avoid log(0)
log_counts = np.log2(counts + 1)
z_scores = log_counts.subtract(log_counts.mean(axis=1), axis=0).divide(log_counts.std(axis=1) + 1e-10, axis=0)

# ──────────────────────────────────────────────
# Hierarchical clustering of genes (rows)
# ──────────────────────────────────────────────
dist_mat = pdist(z_scores.values, metric="euclidean")
link = linkage(dist_mat, method="ward")
row_order = leaves_list(link)

z_sorted = z_scores.iloc[row_order]
gene_ids_sorted = z_sorted.index.tolist()

# Build gene labels with product name
id2meta = genes.set_index("locus_tag")
def make_label(locus):
    row = id2meta.loc[locus] if locus in id2meta.index else None
    if row is None:
        return locus
    gene = row.get("gene_name", "")
    product = row.get("product", "")
    old = row.get("old_locus_tag", "")
    name = gene if (isinstance(gene, str) and gene.strip()) else old
    if not (isinstance(name, str) and name.strip()):
        name = locus
    if isinstance(product, str) and product.strip():
        # shorten product name
        prod_short = product.replace(" transcriptional regulator", " reg.")
        prod_short = prod_short.replace(" family", "")
        prod_short = prod_short.replace("helix-turn-helix domain-containing protein", "HTH protein")
        prod_short = prod_short.replace("sigma-70 family RNA polymerase sigma factor", "σ70 sigma factor")
        prod_short = prod_short.replace("two-component system", "TCS")
        prod_short = prod_short[:50]
        return f"{name} | {prod_short}"
    return name

row_labels = [make_label(g) for g in gene_ids_sorted]

# ──────────────────────────────────────────────
# DEG T2 vs T1 annotation
# ──────────────────────────────────────────────
t2_sig = set(
    deg_t2[(deg_t2["padj"] < 0.05) & (deg_t2["log2FoldChange"] > 1)].index
)
is_t2_up = [g in t2_sig for g in gene_ids_sorted]
print(f"T2 up-regulated (padj<0.05, LFC>1): {sum(is_t2_up)} / {len(gene_ids_sorted)}")

# ──────────────────────────────────────────────
# Figure
# ──────────────────────────────────────────────
n_genes = len(z_sorted)
fig_h = max(12, n_genes * 0.22)
fig, axes = plt.subplots(
    1, 2,
    figsize=(11, fig_h),
    gridspec_kw={"width_ratios": [0.04, 1]},
)

# Left: T2 DEG annotation bar
ax_bar, ax_hm = axes

bar_colors = ["#e74c3c" if up else "#ecf0f1" for up in is_t2_up]
for i, c in enumerate(bar_colors):
    ax_bar.add_patch(mpatches.Rectangle((0, i), 1, 1, color=c))
ax_bar.set_xlim(0, 1)
ax_bar.set_ylim(0, n_genes)
ax_bar.set_xticks([])
ax_bar.set_yticks([])
ax_bar.set_ylabel("")
ax_bar.set_title("T2↑", fontsize=8, pad=2)

# Right: Heatmap
vmax = min(2.5, np.percentile(np.abs(z_sorted.values), 98))
sns.heatmap(
    z_sorted,
    ax=ax_hm,
    cmap="RdBu_r",
    center=0,
    vmin=-vmax,
    vmax=vmax,
    xticklabels=["T1-r1", "T1-r2", "T1-r3",
                 "T2-r1", "T2-r2", "T2-r3",
                 "T3-r1", "T3-r2", "T3-r3"],
    yticklabels=row_labels,
    linewidths=0.0,
    cbar_kws={"label": "Z-score (log₂ normalized counts)", "shrink": 0.4},
    rasterized=True,
)

# Timepoint separator lines
ax_hm.axvline(x=3, color="black", linewidth=1.2)
ax_hm.axvline(x=6, color="black", linewidth=1.2)

ax_hm.set_yticklabels(row_labels, fontsize=6.5)
ax_hm.set_xticklabels(ax_hm.get_xticklabels(), fontsize=8, rotation=45, ha="right")
ax_hm.set_ylabel("")

# Timepoint labels above
ax_hm.text(1.5, n_genes + 0.5, "T1", ha="center", va="bottom", fontsize=9, fontweight="bold",
           transform=ax_hm.transData)
ax_hm.text(4.5, n_genes + 0.5, "T2", ha="center", va="bottom", fontsize=9, fontweight="bold",
           transform=ax_hm.transData)
ax_hm.text(7.5, n_genes + 0.5, "T3", ha="center", va="bottom", fontsize=9, fontweight="bold",
           transform=ax_hm.transData)

fig.suptitle(
    f"Exposed Regulatory Genes (n={n_genes}): Expression Heatmap\n"
    "Z-score normalized log₂ counts, Ward clustering",
    y=1.005, fontsize=10, fontweight="bold"
)

# Legend
red_patch = mpatches.Patch(color="#e74c3c", label="T2 up (padj<0.05, LFC>1)")
ax_hm.legend(
    handles=[red_patch],
    loc="lower right",
    bbox_to_anchor=(1.3, -0.06),
    fontsize=7,
    frameon=True,
)

plt.tight_layout()

out_fig = OUT_DIR / "figures/A7_exposed_TF_heatmap.png"
plt.savefig(out_fig, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print(f"Saved: {out_fig}")

# ──────────────────────────────────────────────
# Save expression table
# ──────────────────────────────────────────────
# Merge gene metadata, normalized counts, z-scores, DEG status
counts_out = counts.copy()
counts_out.columns = [f"norm_{c}" for c in counts_out.columns]

zscore_out = z_sorted.copy()
zscore_out.columns = [f"z_{c}" for c in zscore_out.columns]

meta_cols = ["locus_tag", "gene_name", "old_locus_tag", "product", "tf_family"]
meta_out = genes[meta_cols].set_index("locus_tag")

out_df = meta_out.join(counts_out, how="right").join(zscore_out)

# Add DEG info
out_df["T2_up_DEG"] = out_df.index.isin(t2_sig)
out_df["T2_LFC"] = deg_t2["log2FoldChange"].reindex(out_df.index)
out_df["T2_padj"] = deg_t2["padj"].reindex(out_df.index)
t3_deg = pd.read_csv(
    BASE / "04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_3_vs_1.tsv",
    sep="\t", index_col=0
)
out_df["T3_LFC"] = t3_deg["log2FoldChange"].reindex(out_df.index)
out_df["T3_padj"] = t3_deg["padj"].reindex(out_df.index)

# Add clustering order
out_df["ward_cluster_order"] = [row_order.tolist().index(i) for i in range(len(row_order))]

out_tsv = OUT_DIR / "tables/A7_exposed_TF_expression.tsv"
out_df.to_csv(out_tsv, sep="\t")
print(f"Saved: {out_tsv}")

# ──────────────────────────────────────────────
# Summary statistics
# ──────────────────────────────────────────────
print("\n=== Summary ===")
print(f"Total exposed regulatory genes analyzed: {n_genes}")
print(f"T2 up-regulated (padj<0.05, LFC>1): {sum(is_t2_up)}")

# Identify rough clusters from dendrogram
from scipy.cluster.hierarchy import fcluster
cluster_labels = fcluster(link, t=3, criterion="maxclust")
cluster_labels_sorted = cluster_labels[row_order]

print("\nWard cluster summary (k=3):")
for cl in sorted(set(cluster_labels_sorted)):
    idx = [i for i, c in enumerate(cluster_labels_sorted) if c == cl]
    genes_in = [gene_ids_sorted[i] for i in idx]
    t2_up_n = sum(g in t2_sig for g in genes_in)
    # mean z-score per timepoint
    mz = z_sorted.iloc[idx].mean()
    print(f"  Cluster {cl} (n={len(idx)}): T2_up={t2_up_n}, "
          f"mean_z T1={mz.iloc[:3].mean():.2f}, T2={mz.iloc[3:6].mean():.2f}, T3={mz.iloc[6:].mean():.2f}")

# Add cluster annotation to output
out_df["ward_k3_cluster"] = pd.Series(
    {gene_ids_sorted[i]: cluster_labels_sorted[i] for i in range(n_genes)}
)
out_df.to_csv(out_tsv, sep="\t")
print(f"Updated table with cluster labels: {out_tsv}")
