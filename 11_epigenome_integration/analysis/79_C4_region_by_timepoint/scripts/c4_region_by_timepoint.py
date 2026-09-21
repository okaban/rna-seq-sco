#!/usr/bin/env python3
"""C4: genomic-region composition of methylation sites by modification system
and developmental timepoint (T1/T2/T3).

Driven by the CANONICAL per-timepoint file 01_integration/high_confidence_sites_weighted.csv
(depth>=10, 3 reps, weighted freq>=50%) via 90_/canonical_sites.py, so per-system
per-TP counts match every other figure: GCCGGC 4mC 1,289/1,595/1,073; AAGCCCG-4mC (C4)
698/851/574; AAGCCCG-6mA (A0/A1) 418/451/443.
2026-09-21 (BLOCKER-0): the previous source, 23_expanded_motif_search/*_final_census.csv,
is position-deduplicated across timepoints (first-appearance), which produced the
retracted 1,289/407/21 series and the 'T3 strata small (n = 15-38)' caveat. Region categories use the SAME first-match logic as the
SuppFig 9 classifier (62_GO_KEGG_enrichment/stratified_position_enrichment.py):
promoter = TSS-500..-1; 5'UTR = in-gene 0..+100; CDS_internal = rest of gene;
else intergenic.
"""
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis")
TSS = pd.read_csv(BASE / "18_tss_analyses/comprehensive_tss_table.csv")[
    ["gene_id", "chrom", "start", "end", "strand", "tss"]].rename(
    columns={"start": "gstart", "end": "gend"})
import importlib.util as _ilu
_spec = _ilu.spec_from_file_location("canonical_sites", BASE / "90_per_timepoint_census_audit/canonical_sites.py")
_cs = _ilu.module_from_spec(_spec); _spec.loader.exec_module(_cs)
_HC = _cs.load_canonical()          # columns: position, strand, mod, timepoint, frequency, motif
C4MC = _HC[_HC["mod"] == "4mC"].rename(columns={"motif": "final_motif"})
C6MA = _HC[_HC["mod"] == "6mA"].rename(columns={"motif": "final_motif"})

CATS = ["promoter", "5UTR_approx", "CDS_internal", "intergenic"]


def assign_cat(pos, gstart, gend, strand, tss):
    dist = (pos - tss) if strand == "+" else (tss - pos)
    in_gene = gstart <= pos <= gend
    if not in_gene:
        return "promoter" if -500 <= dist < -1 else "intergenic"
    if 0 <= dist <= 100:
        return "5UTR_approx"
    return "intergenic" if dist < 0 else "CDS_internal"


_G = TSS.to_dict("records")


def classify(df):
    """In-gene first (CDS/5'UTR), then promoter (-500..-1 of a non-overlapped
    TSS), else intergenic. Order-independent; reproduces the SuppFig 9 classifier
    at 99.9% agreement (an in-gene site is CDS-internal, never a neighbour's
    promoter)."""
    out = []
    for pos in df["position"].astype(int):
        cat = None
        for gene in _G:                                   # 1) in-gene priority
            if gene["gstart"] <= pos <= gene["gend"]:
                cat = assign_cat(pos, gene["gstart"], gene["gend"], gene["strand"], gene["tss"])
                if cat == "intergenic":                   # in-gene but upstream of TSS
                    cat = "CDS_internal"
                break
        if cat is None:                                   # 2) promoter
            for gene in _G:
                dist = (pos - gene["tss"]) if gene["strand"] == "+" else (gene["tss"] - pos)
                if -500 <= dist < -1 and not (gene["gstart"] <= pos <= gene["gend"]):
                    cat = "promoter"
                    break
        out.append(cat if cat is not None else "intergenic")
    return out


# 2026-09-13 (FIG-09): locked notation is 4mC (not m4C); panel titles + table
# 'system' column follow it.
SYSTEMS = [
    ("GCCGGC 4mC", C4MC, "GCCGGC"),
    ("AAGCCCG 4mC", C4MC, "AAGCCCG"),
    ("AAGCCCG 6mA", C6MA, "AAGCCCG"),
]
TPS = ["T1", "T2", "T3"]

records = []
for sysname, census, motif in SYSTEMS:
    sub = census[(census["final_motif"] == motif) & (census["frequency"] >= 50)].copy()
    sub["category"] = classify(sub)
    for tp in TPS:
        tpdf = sub[sub["timepoint"] == tp]
        n = len(tpdf)
        counts = tpdf["category"].value_counts().to_dict()
        rec = {"system": sysname, "timepoint": tp, "n": n}
        for c in CATS:
            rec[c] = counts.get(c, 0)
        records.append(rec)

tab = pd.DataFrame(records)
tab.to_csv(BASE / "79_C4_region_by_timepoint/tables/C4_region_composition.tsv",
           sep="\t", index=False)
print(tab.to_string(index=False))

# ---------- figure: stacked % bars, 3 systems x 3 timepoints ----------
COLORS = {"promoter": "#0072B2", "5UTR_approx": "#56B4E9",
          "CDS_internal": "#E69F00", "intergenic": "#999999"}
LBL = {"promoter": "promoter (TSS −500..−1)", "5UTR_approx": "5′UTR (0..+100)",
       "CDS_internal": "CDS-internal", "intergenic": "intergenic"}

fig, axes = plt.subplots(1, 3, figsize=(12, 4.4), sharey=True)
for ax, (sysname, _, _) in zip(axes, SYSTEMS):
    sysrows = tab[tab["system"] == sysname].set_index("timepoint")
    bottom = np.zeros(len(TPS))
    for c in CATS:
        vals = np.array([sysrows.loc[tp, c] for tp in TPS], float)
        ns = np.array([sysrows.loc[tp, "n"] for tp in TPS], float)
        pct = np.where(ns > 0, 100 * vals / ns, 0)
        ax.bar(TPS, pct, bottom=bottom, color=COLORS[c],
               label=LBL[c] if ax is axes[0] else None, width=0.62, edgecolor="white")
        bottom += pct
    for i, tp in enumerate(TPS):
        ax.text(i, 102, f"n={int(sysrows.loc[tp, 'n'])}", ha="center", fontsize=8.5)
    ax.set_title(sysname, fontsize=11)
    ax.set_ylim(0, 112)
    ax.spines[["top", "right"]].set_visible(False)
axes[0].set_ylabel("% of canonical sites (freq ≥ 50%)")
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=9,
           frameon=False, bbox_to_anchor=(0.5, -0.02))
fig.suptitle("Genomic-region composition of methylation by system and developmental timepoint",
             fontsize=12)
fig.tight_layout(rect=[0, 0.05, 1, 0.95])
out = BASE / "79_C4_region_by_timepoint/figures/C4_region_by_timepoint.png"
fig.savefig(out, dpi=200)
print("[saved]", out)
