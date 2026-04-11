#!/usr/bin/env python3
"""
BGC Regulator Overview Figure
Shows expression changes and methylation status for all known BGC regulators,
organized by regulatory tier and mapped to target BGCs.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import TwoSlopeNorm
import matplotlib.gridspec as gridspec

# --- Data definition ---
# All regulators with their expression (LFC) and methylation data
# Source: DESeq2 results + integrated_methyl_expression.csv

regulators = [
    # name, tier, category, locus_tag,
    # LFC_2v1, padj_2v1, LFC_3v1, padj_3v1, LFC_3v2, padj_3v2,
    # 6mA_T1, 6mA_T2, 6mA_T3, 4mC_T1, 4mC_T2, 4mC_T3,
    # known_target_BGCs (act, red, cda, cpk)
    ("bldA",      1, "global",      "SC_RS26640", -1.17, 2.9e-12, -2.11, 2.5e-36, -0.92, 1.6e-7,  0,0,0, 0,0,0, [1,1,0,0]),
    ("bldB",      1, "global",      "SC_RS25655", -0.94, 6.7e-6,   0.02, 0.91,     0.96, 3.8e-6,  0,0,0, 0,0,0, [1,1,1,0]),
    ("bldD",      1, "global",      "SC_RS25165", -0.82, 3.4e-6,  -1.40, 4.1e-16, -0.57, 1.8e-3,  0,0,0, 0,0,0, [1,1,1,1]),
    ("bldG",      1, "global",      "SC_RS25070", -2.75, 1.1e-57, -2.68, 5.5e-57,  0.06, 0.75,    0,0,0, 0,0,0, [1,1,0,0]),
    ("adpA",      1, "global",      "SC_RS14000", -1.27, 6.8e-10, -1.51, 5.1e-14, -0.23, 0.31,    0,0,0, 0,0,0, [1,1,1,0]),
    ("bldN",      1, "global",      "SC_RS26120",  0.50, 0.089,    1.46, 2.2e-7,   0.93, 1.1e-3,  2,1,1, 0,0,0, [1,1,0,0]),
    ("afsR",      1, "global",      "SC_RS24295", -0.26, 0.10,    -0.92, 1.5e-9,  -0.65, 3.0e-5,  1,1,1, 0,0,0, [1,1,1,0]),
    ("afsK",      1, "global",      "SC_RS24290",  3.04, 4.0e-26,  3.07, 1.1e-26,  0.02, 0.94,    0,0,0, 0,0,0, [1,1,1,0]),
    ("afsS",      1, "global",      "SC_RS22980", -0.03, 0.91,    -1.39, 3.4e-12, -1.36, 2.3e-11, 0,0,0, 1,0,0, [1,1,1,0]),
    ("crp",       1, "global",      "SC_RS10725", -0.51, 0.099,    0.49, 0.069,    1.02, 5.5e-4,  0,0,0, 0,0,0, [1,1,0,0]),
    ("dasR",      1, "global",      "SC_RS19305", -0.10, 0.80,     0.91, 3.6e-3,   1.01, 1.8e-3,  0,0,0, 0,0,0, [1,1,1,0]),
    ("absA1",     2, "pleiotropic", "SC_RS02965",  2.90, 1.8e-32,  3.21, 1.2e-39,  0.29, 0.23,    0,0,0, 0,0,0, [1,1,1,0]),
    ("absA2",     2, "pleiotropic", "SC_RS02970",  6.29, 1.8e-81,  7.53, 6.4e-116, 1.18, 7.4e-5,  0,0,0, 0,0,0, [1,1,1,0]),
    ("absB",      2, "pleiotropic", "SC_RS27825", -1.45, 1.8e-4,   0.24, 0.47,     1.72, 8.7e-6,  0,0,0, 0,0,0, [1,1,0,0]),
    ("nsdA",      2, "pleiotropic", "SC_RS35450",  0.28, 0.41,     2.57, 1.7e-18,  2.25, 3.5e-15, 0,0,0, 0,0,0, [1,1,0,0]),
    ("nsdB",      2, "pleiotropic", "SC_RS37265",  1.45, 9.0e-9,   3.38, 3.8e-44,  1.88, 2.5e-16, 0,0,0, 0,0,0, [1,1,0,0]),
    ("wblA",      2, "pleiotropic", "SC_RS04915",  0.34, 0.12,    -0.06, 0.79,    -0.39, 0.059,   0,0,0, 0,0,0, [1,1,1,0]),
    ("hrdD",      2, "sigma",       "SC_RS25840", -2.32, 3.2e-47, -2.84, 1.4e-70, -0.51, 2.5e-3,  0,0,0, 0,0,0, [1,1,0,0]),
    ("sigE",      2, "sigma",       "SC_RS24925", -0.58, 0.21,     2.92, 1.2e-12,  3.61, 8.9e-15, 0,0,0, 0,0,0, [1,0,0,0]),
    ("sigF",      2, "sigma",       "SC_RS04925", -1.26, 6.0e-7,  -1.26, 3.6e-7,   0.00, 1.00,    0,0,0, 0,0,0, [0,1,1,1]),
    ("sigU",      2, "sigma",       "SC_RS04965",  1.49, 3.0e-8,   0.95, 3.2e-4,  -0.51, 0.058,   0,0,0, 0,0,0, [0,1,1,1]),
    ("actII-ORF4",3, "csr",         "SC_RS27570", -0.83, 1.3e-3,   2.60, 1.2e-24,  3.47, 4.7e-41, 0,0,0, 0,0,0, [1,0,0,0]),
    ("redD",      3, "csr",         "SC_RS27225",  4.77, 4.0e-51,  3.57, 2.5e-29, -1.15, 3.8e-7,  0,0,0, 0,0,0, [0,1,0,0]),
    ("redZ",      3, "csr",         "SC_RS27300", -2.25, 1.3e-17, -1.10, 1.7e-6,   1.11, 3.5e-5,  1,0,0, 0,0,0, [0,1,0,0]),
    ("cdaR",      3, "csr",         "SC_RS17785", -0.13, 0.66,     2.42, 8.8e-24,  2.56, 1.8e-23, 0,0,0, 0,0,0, [0,0,1,0]),
    ("cpkO/kasO", 3, "csr",         "SC_RS31605",  0.00, 0.98,    -0.11, 0.65,    -0.11, 0.66,    0,0,0, 0,0,0, [0,0,0,1]),
    ("papR2",     3, "csr",         "SC_RS18200",  6.85, 9.1e-162, 5.76, 1.3e-114,-1.05, 1.4e-5,  0,0,0, 0,0,0, [0,0,1,1]),
]

# Build DataFrame
cols = ["name","tier","category","locus_tag",
        "LFC_2v1","padj_2v1","LFC_3v1","padj_3v1","LFC_3v2","padj_3v2",
        "6mA_T1","6mA_T2","6mA_T3","4mC_T1","4mC_T2","4mC_T3",
        "bgc_targets"]
df = pd.DataFrame(regulators, columns=cols)

# Derived columns
df["has_methyl"] = (df["6mA_T1"] + df["6mA_T2"] + df["6mA_T3"] +
                    df["4mC_T1"] + df["4mC_T2"] + df["4mC_T3"]) > 0

df["methyl_changed"] = False
for i, row in df.iterrows():
    t1 = row["6mA_T1"] + row["4mC_T1"]
    t2 = row["6mA_T2"] + row["4mC_T2"]
    t3 = row["6mA_T3"] + row["4mC_T3"]
    if t1 != t2 or t1 != t3:
        df.at[i, "methyl_changed"] = True

df["expr_sig_3v1"] = df["padj_3v1"] < 0.05
df["expr_sig_2v1"] = df["padj_2v1"] < 0.05

# Methylation status label
def methyl_label(row):
    t1 = row["6mA_T1"] + row["4mC_T1"]
    t2 = row["6mA_T2"] + row["4mC_T2"]
    t3 = row["6mA_T3"] + row["4mC_T3"]
    if t1 == 0 and t2 == 0 and t3 == 0:
        return "Unmethylated"
    elif t1 == t2 == t3:
        return "Stable"
    elif t1 > t2 or t1 > t3:
        return "Lost"
    else:
        return "Gained"

df["methyl_status"] = df.apply(methyl_label, axis=1)

# BGC names
bgc_names = ["act", "red", "cda", "cpk"]
bgc_colors = {"act": "#2166AC", "red": "#B2182B", "cda": "#4DAF4A", "cpk": "#FF7F00"}
bgc_full = {"act": "Actinorhodin", "red": "Undecylprodigiosin", "cda": "CDA", "cpk": "Coelimycin"}

# Category display order and colors
tier_labels = {1: "Tier 1: Global", 2: "Tier 2: Pleiotropic/Sigma", 3: "Tier 3: CSR"}
tier_bg = {1: "#F0F4FF", 2: "#FFF8F0", 3: "#F0FFF0"}

# Sort: by tier, then by category, then by name
cat_order = {"global": 0, "pleiotropic": 1, "sigma": 2, "csr": 3}
df["cat_sort"] = df["category"].map(cat_order)
df = df.sort_values(["tier", "cat_sort", "name"]).reset_index(drop=True)

# ====== FIGURE ======
fig = plt.figure(figsize=(20, 15))
gs = gridspec.GridSpec(1, 3, width_ratios=[3.2, 1.0, 3.8], wspace=0.08,
                       left=0.10, right=0.93, top=0.90, bottom=0.08)

n = len(df)
y_positions = np.arange(n)

# ---- Panel A: Expression heatmap (3 comparisons) ----
ax_expr = fig.add_subplot(gs[0])

lfc_data = df[["LFC_2v1", "LFC_3v1", "LFC_3v2"]].values
padj_data = df[["padj_2v1", "padj_3v1", "padj_3v2"]].values

norm = TwoSlopeNorm(vmin=-4, vcenter=0, vmax=8)
comp_labels = ["T2 vs T1", "T3 vs T1", "T3 vs T2"]

for j in range(3):
    for i in range(n):
        val = lfc_data[i, j]
        sig = padj_data[i, j] < 0.05 if not np.isnan(padj_data[i, j]) else False
        color = plt.cm.RdBu_r(norm(np.clip(val, -4, 8)))
        rect = mpatches.FancyBboxPatch(
            (j - 0.45, i - 0.4), 0.9, 0.8,
            boxstyle="round,pad=0.02",
            facecolor=color,
            edgecolor="gray" if sig else "lightgray",
            linewidth=1.5 if sig else 0.5,
            alpha=1.0 if sig else 0.4
        )
        ax_expr.add_patch(rect)
        if sig and abs(val) >= 0.5:
            txt_color = "white" if abs(val) > 2.5 else "black"
            ax_expr.text(j, i, f"{val:+.1f}", ha="center", va="center",
                        fontsize=7, fontweight="bold", color=txt_color)

ax_expr.set_xlim(-0.6, 2.6)
ax_expr.set_ylim(-0.6, n - 0.4)
ax_expr.set_xticks([0, 1, 2])
ax_expr.set_xticklabels(comp_labels, fontsize=10, fontweight="bold")
# Bold methylated regulator names
name_labels = []
for _, row in df.iterrows():
    name_labels.append(row["name"])

ax_expr.set_yticks(y_positions)
ax_expr.set_yticklabels(name_labels, fontsize=9.5)
# Bold the names of methylated regulators
for i, row in df.iterrows():
    if row["has_methyl"]:
        ax_expr.get_yticklabels()[i].set_fontweight("bold")
        ax_expr.get_yticklabels()[i].set_color("#B22222")
ax_expr.invert_yaxis()
ax_expr.set_title("A. Expression Changes (log$_2$FC)", fontsize=14, fontweight="bold", pad=15)

# Add tier grouping with colored background bands
prev_tier = None
tier_start = 0
for i, row in df.iterrows():
    if prev_tier is not None and row["tier"] != prev_tier:
        ax_expr.axhline(y=i - 0.5, color="black", linewidth=1.5, linestyle="-")
        # Add background
        ax_expr.axhspan(tier_start - 0.5, i - 0.5, color=tier_bg[prev_tier],
                       alpha=0.3, zorder=0)
        mid = (tier_start + i - 1) / 2
        ax_expr.text(-0.62, mid, tier_labels[prev_tier], fontsize=8,
                    rotation=90, va="center", ha="right", fontstyle="italic",
                    color="#555555", fontweight="bold")
        tier_start = i
    prev_tier = row["tier"]
ax_expr.axhspan(tier_start - 0.5, n - 0.5, color=tier_bg[prev_tier],
               alpha=0.3, zorder=0)
mid = (tier_start + n - 1) / 2
ax_expr.text(-0.62, mid, tier_labels[prev_tier], fontsize=8,
            rotation=90, va="center", ha="right", fontstyle="italic",
            color="#555555", fontweight="bold")

# Add significance legend note
ax_expr.text(0, n + 0.5, "Bold border = padj < 0.05; faded = n.s.  |  Red bold name = has methylation site",
            fontsize=7.5, color="gray", ha="left")

# ---- Panel B: Methylation status ----
ax_met = fig.add_subplot(gs[1])

met_labels_tp = ["T1", "T2", "T3"]
for j, tp in enumerate(["T1", "T2", "T3"]):
    for i in range(n):
        total = df.iloc[i][f"6mA_{tp}"] + df.iloc[i][f"4mC_{tp}"]
        if total > 0:
            # Filled circle - size proportional to count
            size = 80 + total * 120
            mod_type = "6mA" if df.iloc[i][f"6mA_{tp}"] > 0 else "4mC"
            color = "#E07B39" if mod_type == "6mA" else "#7B39E0"
            ax_met.scatter(j, i, s=size, c=color, edgecolors="black",
                          linewidth=1.2, zorder=5)
            ax_met.text(j, i, str(int(total)), ha="center", va="center",
                       fontsize=7, fontweight="bold", color="white", zorder=6)
        else:
            ax_met.scatter(j, i, s=30, c="white", edgecolors="lightgray",
                          linewidth=0.5, marker="o", zorder=3)

# Draw arrows for changes
for i in range(n):
    t1 = df.iloc[i]["6mA_T1"] + df.iloc[i]["4mC_T1"]
    t3 = df.iloc[i]["6mA_T3"] + df.iloc[i]["4mC_T3"]
    if t1 > 0 and t3 < t1:
        # Lost methylation arrow
        ax_met.annotate("", xy=(2.35, i), xytext=(2.55, i),
                        arrowprops=dict(arrowstyle="->", color="red", lw=1.5))
        ax_met.text(2.7, i, "Lost", fontsize=6, color="red", va="center")
    elif t1 > 0 and t3 == t1:
        ax_met.text(2.55, i, "Stable", fontsize=6, color="teal", va="center")

ax_met.set_xlim(-0.5, 3.5)
ax_met.set_ylim(-0.6, n - 0.4)
ax_met.set_xticks([0, 1, 2])
ax_met.set_xticklabels(met_labels_tp, fontsize=10, fontweight="bold")
ax_met.set_yticks([])
ax_met.invert_yaxis()
ax_met.set_title("B. Methylation\n    Sites", fontsize=14, fontweight="bold", pad=15)

# Add tier separators + background
prev_tier = None
tier_start = 0
for i, row in df.iterrows():
    if prev_tier is not None and row["tier"] != prev_tier:
        ax_met.axhline(y=i - 0.5, color="black", linewidth=1.5, linestyle="-")
        ax_met.axhspan(tier_start - 0.5, i - 0.5, color=tier_bg[prev_tier],
                       alpha=0.3, zorder=0)
        tier_start = i
    prev_tier = row["tier"]
ax_met.axhspan(tier_start - 0.5, n - 0.5, color=tier_bg[prev_tier],
               alpha=0.3, zorder=0)

# Methylation legend
leg_6mA = mpatches.Patch(facecolor="#E07B39", edgecolor="black", label="6mA")
leg_4mC = mpatches.Patch(facecolor="#7B39E0", edgecolor="black", label="4mC")
ax_met.legend(handles=[leg_6mA, leg_4mC], loc="lower center",
             fontsize=7, frameon=True, ncol=2, bbox_to_anchor=(0.3, -0.06))

# ---- Panel C: BGC target matrix with dual annotation ----
ax_bgc = fig.add_subplot(gs[2])

for j, bgc in enumerate(bgc_names):
    for i in range(n):
        targets = df.iloc[i]["bgc_targets"]
        if targets[j] == 0:
            ax_bgc.scatter(j, i, s=20, c="white", edgecolors="#DDDDDD",
                          linewidth=0.3, marker="s", zorder=2)
            continue

        # This regulator targets this BGC
        lfc = df.iloc[i]["LFC_3v1"]
        sig = df.iloc[i]["padj_3v1"] < 0.05
        has_met_change = df.iloc[i]["methyl_changed"]
        has_met = df.iloc[i]["has_methyl"]

        # Determine category: both, expression only, methylation only, neither
        if sig and has_met_change:
            # BOTH: star marker, large
            marker = "*"
            size = 350
            edge_color = "black"
            edge_width = 2.0
        elif sig and has_met and not has_met_change:
            # Expression + stable methylation
            marker = "D"
            size = 120
            edge_color = "black"
            edge_width = 1.5
        elif sig:
            # Expression only
            marker = "o"
            size = 150
            edge_color = "black"
            edge_width = 1.0
        elif has_met_change:
            # Methylation only
            marker = "^"
            size = 120
            edge_color = "black"
            edge_width = 1.5
        else:
            # No significant change
            marker = "o"
            size = 60
            edge_color = "gray"
            edge_width = 0.5

        # Color by LFC direction/magnitude
        if sig:
            color = plt.cm.RdBu_r(norm(np.clip(lfc, -4, 8)))
        else:
            color = "#CCCCCC"

        ax_bgc.scatter(j, i, s=size, c=[color], edgecolors=edge_color,
                      linewidth=edge_width, marker=marker, zorder=5)

        # Annotate BOTH cases with highlight box
        if sig and has_met_change:
            rect_bg = mpatches.FancyBboxPatch(
                (j - 0.42, i - 0.42), 0.84, 0.84,
                boxstyle="round,pad=0.05",
                facecolor="yellow", alpha=0.2,
                edgecolor="orange", linewidth=1.5, linestyle="--",
                zorder=1
            )
            ax_bgc.add_patch(rect_bg)

ax_bgc.set_xlim(-0.6, 3.6)
ax_bgc.set_ylim(-0.6, n - 0.4)
ax_bgc.set_xticks(range(4))
ax_bgc.set_xticklabels([f"{bgc}\n({bgc_full[bgc]})" for bgc in bgc_names],
                        fontsize=9, fontweight="bold")
for j, bgc in enumerate(bgc_names):
    ax_bgc.get_xticklabels()[j].set_color(bgc_colors[bgc])
ax_bgc.set_yticks(y_positions)
ax_bgc.set_yticklabels(df["name"], fontsize=9.5)
for i, row in df.iterrows():
    if row["has_methyl"]:
        ax_bgc.get_yticklabels()[i].set_fontweight("bold")
        ax_bgc.get_yticklabels()[i].set_color("#B22222")
ax_bgc.yaxis.set_label_position("right")
ax_bgc.yaxis.tick_right()
ax_bgc.invert_yaxis()
ax_bgc.set_title("C. Regulator\u2013BGC Mapping (T3 vs T1, color = log$_2$FC)",
                  fontsize=14, fontweight="bold", pad=15)

# Tier separators + background
prev_tier = None
tier_start = 0
for i, row in df.iterrows():
    if prev_tier is not None and row["tier"] != prev_tier:
        ax_bgc.axhline(y=i - 0.5, color="black", linewidth=1.5, linestyle="-")
        ax_bgc.axhspan(tier_start - 0.5, i - 0.5, color=tier_bg[prev_tier],
                       alpha=0.3, zorder=0)
        tier_start = i
    prev_tier = row["tier"]
ax_bgc.axhspan(tier_start - 0.5, n - 0.5, color=tier_bg[prev_tier],
               alpha=0.3, zorder=0)

# Legend for Panel C
leg_items = [
    plt.scatter([], [], s=350, c="gold", edgecolors="black", linewidth=2.0,
               marker="*", label="Expression + Methylation BOTH changed"),
    plt.scatter([], [], s=120, c="steelblue", edgecolors="black", linewidth=1.5,
               marker="D", label="Expression changed + Stable methylation"),
    plt.scatter([], [], s=150, c="steelblue", edgecolors="black", linewidth=1.0,
               marker="o", label="Expression changed only"),
    plt.scatter([], [], s=60, c="#CCCCCC", edgecolors="gray", linewidth=0.5,
               marker="o", label="No significant change"),
]
dashed_patch = mpatches.FancyBboxPatch(
    (0, 0), 1, 1, boxstyle="round,pad=0.05",
    facecolor="yellow", alpha=0.3, edgecolor="orange",
    linewidth=1.5, linestyle="--"
)
ax_bgc.legend(handles=leg_items, loc="lower center",
             fontsize=7.5, frameon=True, ncol=1, bbox_to_anchor=(0.5, -0.14),
             title="Marker type", title_fontsize=8)

# Colorbar for LFC
sm = plt.cm.ScalarMappable(cmap=plt.cm.RdBu_r, norm=norm)
sm.set_array([])
cax = fig.add_axes([0.35, 0.03, 0.30, 0.015])
cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
cbar.set_label("log$_2$ Fold Change (blue = down, red = up)", fontsize=10)
cbar.ax.tick_params(labelsize=9)

# Main title
fig.suptitle("BGC Regulatory Network: Expression Dynamics and Methylation Status\n"
             "Streptomyces coelicolor M145 — Three Growth Phases",
             fontsize=16, fontweight="bold", y=0.97)

# Save
outdir = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/19_bgc_regulator_overview"
fig.savefig(f"{outdir}/bgc_regulator_overview.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{outdir}/bgc_regulator_overview.pdf", dpi=300, bbox_inches="tight")
fig.savefig(f"{outdir}/bgc_regulator_overview.svg", dpi=300, bbox_inches="tight")
print(f"Saved to {outdir}/bgc_regulator_overview.[png|pdf|svg]")

# --- Summary table ---
summary_rows = []
for _, row in df.iterrows():
    for j, bgc in enumerate(bgc_names):
        if row["bgc_targets"][j] == 1:
            sig = row["padj_3v1"] < 0.05
            mc = row["methyl_changed"]
            status = "BOTH" if sig and mc else ("Expr only" if sig else ("Methyl only" if mc else "Neither"))
            summary_rows.append({
                "Regulator": row["name"],
                "Tier": row["tier"],
                "BGC": bgc,
                "LFC_T3vT1": round(row["LFC_3v1"], 2),
                "padj_T3vT1": f"{row['padj_3v1']:.2e}",
                "Methylation_Status": row["methyl_status"],
                "Dual_Change": status
            })

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(f"{outdir}/bgc_regulator_summary.tsv", sep="\t", index=False)
print(f"Saved summary table to {outdir}/bgc_regulator_summary.tsv")

# Print dual-change summary
both = summary_df[summary_df["Dual_Change"] == "BOTH"]
if len(both) > 0:
    print("\n=== Regulators with BOTH expression AND methylation changes ===")
    print(both[["Regulator","BGC","LFC_T3vT1","Methylation_Status"]].to_string(index=False))
else:
    print("\nNo regulators with BOTH expression AND methylation changes found.")
