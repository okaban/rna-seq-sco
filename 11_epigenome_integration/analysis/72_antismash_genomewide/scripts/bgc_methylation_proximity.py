"""BGC × DNA methylation proximity analysis for S. coelicolor M145.

Computes per-BGC methylation site counts (4mC / 6mA × T1/T2/T3) inside the
BGC body and within ±5 kb flanking, then density (sites/kb).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path("/Users/okaban/bioinfo/rna-seq")
ANALYSIS = ROOT / "11_epigenome_integration/analysis"
BGC_TSV = ANALYSIS / "72_antismash_genomewide/bgc_summary.tsv"
SITES_CSV = ANALYSIS / "high_confidence_sites_weighted.csv"
OUT_TSV = ANALYSIS / "72_antismash_genomewide/bgc_methylation_proximity.tsv"
OUT_FIG = ANALYSIS / "figures/Fig_BGC_methylation_density.png"
OUT_FIG_OBS = Path(
    "/Users/okaban/obsidian/Research/rna-seq/Writing/Fig_BGC_methylation_density.png"
)

FLANK_BP = 5000
MOD_TYPES = ("4mC", "6mA")
TIMEPOINTS = ("T1", "T2", "T3")
MAIN_BGCS = {
    "NC_003888.3_region11": "CDA",
    "NC_003888.3_region12": "ACT",
    "NC_003888.3_region17": "RED",
    "NC_003888.3_region21": "CPK",
}


def short_label(row: pd.Series) -> str:
    """Concise label for plotting: short name if main BGC else top_known_name."""
    if row["region"] in MAIN_BGCS:
        return MAIN_BGCS[row["region"]]
    name = row.get("top_known_name")
    if isinstance(name, str) and name:
        return name.split("/")[0][:18]
    return row["products"].split(",")[0][:18] if isinstance(row["products"], str) else row["region"]


def count_in_range(sites: pd.DataFrame, lo: int, hi: int) -> int:
    return int(((sites["position"] >= lo) & (sites["position"] <= hi)).sum())


def main() -> None:
    bgc = pd.read_csv(BGC_TSV, sep="\t")
    sites = pd.read_csv(SITES_CSV)

    sites_idx = {
        (mod, tp): sites[(sites["mod_type"] == mod) & (sites["timepoint"] == tp)][
            ["chrom", "position"]
        ].reset_index(drop=True)
        for mod in MOD_TYPES
        for tp in TIMEPOINTS
    }

    rows = []
    for _, b in bgc.iterrows():
        chrom = b["record"]
        start = int(b["start"])
        end = int(b["end"])
        length_kb = (end - start + 1) / 1000.0
        flank_left_lo, flank_left_hi = max(1, start - FLANK_BP), start - 1
        flank_right_lo, flank_right_hi = end + 1, end + FLANK_BP
        flank_kb = ((flank_left_hi - flank_left_lo + 1) + (flank_right_hi - flank_right_lo + 1)) / 1000.0

        out = {
            "bgc_id": b["region"],
            "bgc_label": short_label(b),
            "is_main4": b["region"] in MAIN_BGCS,
            "bgc_type": b["products"],
            "known_cluster": b.get("top_known_name"),
            "start": start,
            "end": end,
            "length_kb": round(length_kb, 3),
        }

        for mod in MOD_TYPES:
            for tp in TIMEPOINTS:
                df = sites_idx[(mod, tp)]
                df_chr = df[df["chrom"] == chrom]
                n_int = count_in_range(df_chr, start, end)
                n_fl = count_in_range(df_chr, flank_left_lo, flank_left_hi) + count_in_range(
                    df_chr, flank_right_lo, flank_right_hi
                )
                out[f"{mod}_{tp}_internal"] = n_int
                out[f"{mod}_{tp}_flank"] = n_fl
                out[f"density_{mod}_{tp}"] = round(n_int / length_kb, 4) if length_kb > 0 else 0.0
                out[f"density_flank_{mod}_{tp}"] = round(n_fl / flank_kb, 4) if flank_kb > 0 else 0.0
        rows.append(out)

    df = pd.DataFrame(rows)
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"[write] {OUT_TSV} ({len(df)} BGCs)")

    print("\n=== Genome-wide background density (sites/kb across full chromosome) ===")
    chrom_len_kb = (sites["position"].max() - sites["position"].min()) / 1000.0
    bg = {}
    for mod in MOD_TYPES:
        for tp in TIMEPOINTS:
            n = len(sites_idx[(mod, tp)])
            bg[(mod, tp)] = n / chrom_len_kb
            print(f"  {mod} {tp}: {n} sites, density={bg[(mod, tp)]:.4f} sites/kb")

    print("\n=== Mean internal density: main4 vs others ===")
    summary = []
    for mod in MOD_TYPES:
        for tp in TIMEPOINTS:
            col = f"density_{mod}_{tp}"
            m4 = df[df["is_main4"]][col].mean()
            ot = df[~df["is_main4"]][col].mean()
            print(f"  {mod} {tp}:  main4={m4:.4f}  others={ot:.4f}  bg={bg[(mod, tp)]:.4f}")
            summary.append({
                "mod": mod, "tp": tp,
                "main4_mean": m4, "others_mean": ot, "background": bg[(mod, tp)],
            })

    print("\n=== T1→T3 4mC density change (per BGC) ===")
    df["delta_4mC_T1T3"] = df["density_4mC_T3"] - df["density_4mC_T1"]
    print(df[["bgc_id", "bgc_label", "is_main4",
              "density_4mC_T1", "density_4mC_T2", "density_4mC_T3",
              "delta_4mC_T1T3"]].sort_values("delta_4mC_T1T3").to_string(index=False))

    print("\n=== Top BGCs by combined methylation density (T2 sum 4mC+6mA) ===")
    df["combined_T2_density"] = df["density_4mC_T2"] + df["density_6mA_T2"]
    print(df.sort_values("combined_T2_density", ascending=False)[
        ["bgc_id", "bgc_label", "is_main4", "combined_T2_density",
         "density_4mC_T2", "density_6mA_T2"]
    ].head(10).to_string(index=False))

    make_figure(df, bg)
    print(f"\n[write] {OUT_FIG}")


def make_figure(df: pd.DataFrame, bg: dict) -> None:
    """Two-panel figure: (a) per-BGC density bars, (b) main4 vs others summary."""
    df_sorted = df.sort_values("start").reset_index(drop=True)
    n = len(df_sorted)
    x = np.arange(n)
    width = 0.4

    fig, axes = plt.subplots(2, 1, figsize=(13, 8), gridspec_kw={"height_ratios": [3, 2]})

    ax = axes[0]
    ax.bar(x - width/2, df_sorted["density_4mC_T2"], width, label="4mC T2", color="#1f77b4")
    ax.bar(x + width/2, df_sorted["density_6mA_T2"], width, label="6mA T2", color="#d62728")
    ax.axhline(bg[("4mC", "T2")], color="#1f77b4", linestyle="--", alpha=0.5, linewidth=1,
               label=f"4mC genome bg ({bg[('4mC','T2')]:.3f})")
    ax.axhline(bg[("6mA", "T2")], color="#d62728", linestyle="--", alpha=0.5, linewidth=1,
               label=f"6mA genome bg ({bg[('6mA','T2')]:.3f})")
    ax.set_xticks(x)
    labels = [f"{r['bgc_label']}{'*' if r['is_main4'] else ''}" for _, r in df_sorted.iterrows()]
    ax.set_xticklabels(labels, rotation=70, ha="right", fontsize=8)
    ax.set_ylabel("Density (sites/kb), T2")
    ax.set_title("Per-BGC methylation density at T2  (asterisk = main 4 BGCs: ACT, CDA, RED, CPK)")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(axis="y", alpha=0.3)

    ax2 = axes[1]
    cats = []
    main_vals = {"4mC": [], "6mA": []}
    others_vals = {"4mC": [], "6mA": []}
    bg_vals = {"4mC": [], "6mA": []}
    for tp in TIMEPOINTS:
        cats.append(tp)
        for mod in MOD_TYPES:
            main_vals[mod].append(df[df["is_main4"]][f"density_{mod}_{tp}"].mean())
            others_vals[mod].append(df[~df["is_main4"]][f"density_{mod}_{tp}"].mean())
            bg_vals[mod].append(bg[(mod, tp)])

    xc = np.arange(len(cats))
    w = 0.18
    ax2.bar(xc - 2*w, main_vals["4mC"], w, label="4mC main4", color="#1f77b4")
    ax2.bar(xc - w, others_vals["4mC"], w, label="4mC others", color="#7fb1d3")
    ax2.bar(xc, bg_vals["4mC"], w, label="4mC bg", color="#7fb1d3", alpha=0.4, hatch="//")
    ax2.bar(xc + w, main_vals["6mA"], w, label="6mA main4", color="#d62728")
    ax2.bar(xc + 2*w, others_vals["6mA"], w, label="6mA others", color="#f4a3a3")
    ax2.set_xticks(xc)
    ax2.set_xticklabels(cats)
    ax2.set_ylabel("Mean density (sites/kb)")
    ax2.set_title("Mean per-BGC density: main4 vs other 25 BGCs (with genome background)")
    ax2.legend(ncol=3, fontsize=8, loc="upper right")
    ax2.grid(axis="y", alpha=0.3)

    fig.tight_layout()
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG, dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
