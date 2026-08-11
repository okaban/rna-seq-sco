#!/usr/bin/env python3
"""Layer 2 (sequence counter-selection) vs Layer 3 (NAP occupancy) overlap.

Reviewer Minor-3-5 asks whether the ~33% counter-selection contribution and
~67% NAP occupancy contribution to the regulatory promoter protection zone
describe the same sites or independent phenomena.

Operationalisation (gene-level, in the 2-kb extended window matching H22/H30):

  expected_motifs   = genome_baseline_density * window_kb
  observed_motifs   = TGGCCGGC count in window (sequence)
  methylated_T1     = T1 4mC GCCGGC sites in window (methylation calls)

  Layer-2 deficit   = max(0, expected - observed)   # sites missing from sequence
  Layer-3 deficit   = max(0, observed - methylated) # sites present but unmethylated

A gene is a Layer-2 contributor if its Layer-2 deficit > 0.
A gene is a Layer-3 contributor if its Layer-3 deficit > 0.

These two conditions are NOT mutually exclusive at the gene level: a single
gene can simultaneously have fewer motifs than expected AND have those that
remain be unmethylated. The overlap is thus an empirical question.

A site-level note: the underlying Layer-2 and Layer-3 *positions* are
mutually exclusive by construction (a position either has the GCCGGC
sequence or does not), so site-level Jaccard = 0 trivially. The interesting
question is whether the gene populations contributing to each layer are
the same set of genes (high gene-level Jaccard) or different sets (low
gene-level Jaccard).
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

BASE = Path("/Users/okaban/bioinfo/rna-seq")
ANALYSIS = BASE / "11_epigenome_integration/analysis"
OUT_DIR = ANALYSIS / "73_layer_overlap"
TBL_DIR = OUT_DIR / "tables"
TBL_DIR.mkdir(parents=True, exist_ok=True)

EXTENDED_UPSTREAM = 2000
EXTENDED_DOWNSTREAM = 2000
GENOME_LEN = 8_667_507

# Boundary defining shielded vs exposed regulators (H29, AUC=0.917)
SHIELD_BOUNDARY = 293

# TSS-proximal protection-zone window for the secondary analysis.
# The paper reports a ~1.2 kb methylation-free promoter protection zone
# centred near the TSS; we apply a symmetric 600-bp half-window.
PROTZONE_HALF = 600


def main() -> None:
    # 1. Per-gene motif counts (H22) — gives extended_count_TGGCCGGC, the
    #    number of TGGCCGGC sites in [gene_start-2000, gene_end+2000].
    motif = pd.read_csv(
        ANALYSIS / "45_sequence_level_motif_depletion/tables/gene_motif_counts.tsv",
        sep="\t",
    )
    print(f"[load] gene_motif_counts: {len(motif):>6} rows")

    # Genome-wide TGGCCGGC density (sites/kb): use gene-body counts ÷ total
    # gene-body length (independent of any flanking-region effects).
    body_sites_total = motif["body_count_TGGCCGGC"].sum()
    body_len_total = (motif["end"] - motif["start"]).sum()
    genome_baseline_density = body_sites_total / (body_len_total / 1000.0)
    print(
        f"[baseline] gene-body TGGCCGGC density = "
        f"{genome_baseline_density:.4f} sites/kb (n={body_sites_total} sites, "
        f"len={body_len_total/1000:.1f} kb)"
    )

    # 2. Regulatory genes with TSS info (H27)
    reg = pd.read_csv(
        ANALYSIS / "50_coordinated_regulators_protection/tables/gene_level_metrics.tsv",
        sep="\t",
    )
    print(f"[load] regulatory genes with TSS+methyl_distance: {len(reg):>6}")

    reg_valid = reg.dropna(subset=["nearest_methyl_dist", "tss"]).copy()
    reg_valid["is_shielded"] = reg_valid["nearest_methyl_dist"] > SHIELD_BOUNDARY
    n_shield = int(reg_valid["is_shielded"].sum())
    n_exposed = int((~reg_valid["is_shielded"]).sum())
    print(f"[partition] shielded={n_shield}, exposed={n_exposed} (boundary={SHIELD_BOUNDARY} bp)")

    shield = reg_valid[reg_valid["is_shielded"]].copy()

    # 3. Merge extended motif count + start/end onto shielded genes
    motif_sub = motif[
        ["gene_id", "start", "end", "extended_count_TGGCCGGC"]
    ].rename(
        columns={
            "gene_id": "locus_tag",
            "start": "motif_start",
            "end": "motif_end",
            "extended_count_TGGCCGGC": "observed_motifs",
        }
    )
    shield = shield.merge(motif_sub, on="locus_tag", how="left")
    n_with_motif = shield["observed_motifs"].notna().sum()
    print(f"[merge] shielded with motif count: {n_with_motif}/{len(shield)}")
    shield = shield.dropna(subset=["observed_motifs"]).copy()
    shield["observed_motifs"] = shield["observed_motifs"].astype(int)

    # Extended window length for each gene (clamped to genome bounds)
    ext_start = (shield["motif_start"] - EXTENDED_UPSTREAM).clip(lower=0)
    ext_end = (shield["motif_end"] + EXTENDED_DOWNSTREAM).clip(upper=GENOME_LEN)
    shield["extended_kb"] = (ext_end - ext_start) / 1000.0

    # 4. Expected motifs under uniform genome model
    shield["expected_motifs"] = (
        genome_baseline_density * shield["extended_kb"]
    )

    # 5. T1 GCCGGC methylated sites in the same extended window
    sites = pd.read_csv(
        ANALYSIS / "37_defense_island_GCCGGC/tables/GCCGGC_sites_by_timepoint.tsv",
        sep="\t",
    )
    t1 = sites[sites["timepoint"] == "T1"].copy()
    t1_pos_sorted = np.sort(t1["position"].to_numpy())
    print(f"[load] T1 GCCGGC methylated sites: {len(t1_pos_sorted)}")

    def n_t1_in_extended(row: pd.Series) -> int:
        lo = max(0, int(row["motif_start"]) - EXTENDED_UPSTREAM)
        hi = min(GENOME_LEN, int(row["motif_end"]) + EXTENDED_DOWNSTREAM)
        i = np.searchsorted(t1_pos_sorted, lo)
        j = np.searchsorted(t1_pos_sorted, hi, side="right")
        return int(j - i)

    shield["methylated_T1"] = shield.apply(n_t1_in_extended, axis=1)

    # 6. Per-gene deficit decomposition
    shield["layer2_deficit"] = (shield["expected_motifs"] - shield["observed_motifs"]).clip(lower=0)
    shield["layer3_deficit"] = (shield["observed_motifs"] - shield["methylated_T1"]).clip(lower=0)
    shield["total_deficit"] = (shield["expected_motifs"] - shield["methylated_T1"]).clip(lower=0)

    # Population-level decomposition (sanity-check the 33/67 split)
    sum_l2 = shield["layer2_deficit"].sum()
    sum_l3 = shield["layer3_deficit"].sum()
    sum_total = sum_l2 + sum_l3
    pct_l2 = sum_l2 / sum_total * 100 if sum_total else float("nan")
    pct_l3 = sum_l3 / sum_total * 100 if sum_total else float("nan")
    print()
    print(f"[population deficit decomposition over n={len(shield)} shielded genes]")
    print(f"   sum Layer-2 deficit (sequence-absent) = {sum_l2:.1f}")
    print(f"   sum Layer-3 deficit (present-unmethylated) = {sum_l3:.1f}")
    print(f"   Layer-2 share = {pct_l2:.1f}% (paper target ~33%)")
    print(f"   Layer-3 share = {pct_l3:.1f}% (paper target ~67%)")

    # 7. Gene-level set membership
    shield["is_layer2"] = shield["layer2_deficit"] > 0
    shield["is_layer3"] = shield["layer3_deficit"] > 0

    L2 = set(shield.loc[shield["is_layer2"], "locus_tag"].tolist())
    L3 = set(shield.loc[shield["is_layer3"], "locus_tag"].tolist())
    universe = set(shield["locus_tag"].tolist())
    only_L2 = L2 - L3
    only_L3 = L3 - L2
    both = L2 & L3
    union = L2 | L3
    neither = universe - union

    jaccard = len(both) / len(union) if union else 0.0

    contingency = np.array(
        [
            [len(both), len(only_L2)],
            [len(only_L3), len(neither)],
        ]
    )
    odds, p_fisher = stats.fisher_exact(contingency, alternative="two-sided")

    # Expected overlap if independent
    p_l2 = len(L2) / len(universe)
    p_l3 = len(L3) / len(universe)
    expected_both = p_l2 * p_l3 * len(universe)

    print()
    print("=" * 64)
    print("LAYER 2 / LAYER 3 OVERLAP — gene-level (shielded regulators)")
    print("=" * 64)
    print(f"Universe (shielded with motif count): N = {len(universe)}")
    print(f"|Layer 2|       = {len(L2)}  ({p_l2:.1%})")
    print(f"|Layer 3|       = {len(L3)}  ({p_l3:.1%})")
    print(f"|Layer 2 only|  = {len(only_L2)}")
    print(f"|Layer 3 only|  = {len(only_L3)}")
    print(f"|Layer 2 ∩ Layer 3| = {len(both)}  (expected if indep ≈ {expected_both:.1f})")
    print(f"|Neither|        = {len(neither)}")
    print(f"|Layer 2 ∪ Layer 3| = {len(union)}")
    print(f"Jaccard         = {jaccard:.4f}")
    print(f"Fisher OR       = {odds:.3f}, two-sided p = {p_fisher:.3g}")
    print()
    if jaccard >= 0.5:
        interp = "HIGH (≥0.5): Layers describe LARGELY THE SAME GENES"
    elif jaccard >= 0.2:
        interp = "MODERATE (0.2–0.5): Partially shared gene populations"
    else:
        interp = "LOW (<0.2): Gene populations are LARGELY INDEPENDENT"
    print(f"Interpretation : {interp}")

    # 8. Save tables
    out_summary = pd.DataFrame(
        {
            "metric": [
                "n_shielded_universe",
                "n_layer2",
                "n_layer3",
                "n_layer2_only",
                "n_layer3_only",
                "n_intersection",
                "n_neither",
                "n_union",
                "expected_intersection_if_independent",
                "jaccard",
                "fisher_odds_ratio",
                "fisher_p_value",
                "layer2_share_of_universe",
                "layer3_share_of_universe",
                "population_layer2_pct_of_total_deficit",
                "population_layer3_pct_of_total_deficit",
                "genome_baseline_density_per_kb",
                "shield_boundary_bp",
                "window_kb_per_gene_each_side",
            ],
            "value": [
                len(universe),
                len(L2),
                len(L3),
                len(only_L2),
                len(only_L3),
                len(both),
                len(neither),
                len(union),
                round(expected_both, 2),
                round(jaccard, 4),
                round(odds, 4),
                p_fisher,
                round(p_l2, 4),
                round(p_l3, 4),
                round(pct_l2, 2),
                round(pct_l3, 2),
                round(genome_baseline_density, 4),
                SHIELD_BOUNDARY,
                EXTENDED_UPSTREAM / 1000,
            ],
        }
    )
    summary_path = TBL_DIR / "layer2_layer3_overlap.tsv"
    out_summary.to_csv(summary_path, sep="\t", index=False)
    print(f"\n[write] {summary_path}")

    # Per-gene assignments
    shield_out = shield[
        [
            "locus_tag",
            "gene_name",
            "tss",
            "strand",
            "region",
            "nearest_methyl_dist",
            "extended_kb",
            "expected_motifs",
            "observed_motifs",
            "methylated_T1",
            "layer2_deficit",
            "layer3_deficit",
            "is_layer2",
            "is_layer3",
        ]
    ].copy()
    shield_out["category"] = np.select(
        [
            shield_out["is_layer2"] & shield_out["is_layer3"],
            shield_out["is_layer2"] & ~shield_out["is_layer3"],
            ~shield_out["is_layer2"] & shield_out["is_layer3"],
        ],
        ["both", "layer2_only", "layer3_only"],
        default="neither",
    )
    per_gene_path = TBL_DIR / "layer2_layer3_per_gene.tsv"
    shield_out.to_csv(per_gene_path, sep="\t", index=False)
    print(f"[write] {per_gene_path}")

    # ------------------------------------------------------------------
    # 9. Secondary analysis: TSS-proximal 1.2-kb protection zone
    # ------------------------------------------------------------------
    print()
    print("=" * 64)
    print(f"SECONDARY: TSS-proximal protection zone (±{PROTZONE_HALF} bp)")
    print("=" * 64)

    # Load full genome FASTA-derived motif sites: re-scan via gene-body density
    # We use the same expected-density baseline. We must count motifs in the
    # symmetric TSS window from the TGGCCGGC site list itself, since the
    # gene_motif_counts.tsv promoter window is asymmetric (-500/+100).
    #
    # The full set of TGGCCGGC genome positions is the union of T1/T2/T3
    # methylated sites plus any unmethylated ones; we approximate by using
    # all GCCGGC sites that appear in any timepoint as the "sequence" set.
    all_motif_pos = np.sort(sites["position"].drop_duplicates().to_numpy())
    print(f"[load] union of GCCGGC sites across timepoints: {len(all_motif_pos)}")

    def n_in_window(arr: np.ndarray, lo: int, hi: int) -> int:
        i = np.searchsorted(arr, lo)
        j = np.searchsorted(arr, hi, side="right")
        return int(j - i)

    rows = []
    for _, row in shield.iterrows():
        tss = int(row["tss"])
        lo = max(0, tss - PROTZONE_HALF)
        hi = min(GENOME_LEN, tss + PROTZONE_HALF)
        win_kb = (hi - lo) / 1000.0
        observed = n_in_window(all_motif_pos, lo, hi)
        methylated = n_in_window(t1_pos_sorted, lo, hi)
        expected = genome_baseline_density * win_kb
        rows.append(
            {
                "locus_tag": row["locus_tag"],
                "expected": expected,
                "observed": observed,
                "methylated_T1": methylated,
                "layer2_deficit": max(0.0, expected - observed),
                "layer3_deficit": max(0.0, observed - methylated),
            }
        )
    pz = pd.DataFrame(rows)

    pz_l2_sum = pz["layer2_deficit"].sum()
    pz_l3_sum = pz["layer3_deficit"].sum()
    pz_total = pz_l2_sum + pz_l3_sum
    pz_pct_l2 = pz_l2_sum / pz_total * 100 if pz_total else float("nan")
    pz_pct_l3 = pz_l3_sum / pz_total * 100 if pz_total else float("nan")
    print(f"[promoter zone deficit decomposition over n={len(pz)} genes]")
    print(f"   sum Layer-2 deficit = {pz_l2_sum:.1f}")
    print(f"   sum Layer-3 deficit = {pz_l3_sum:.1f}")
    print(f"   Layer-2 share = {pz_pct_l2:.1f}%  (paper target ~33%)")
    print(f"   Layer-3 share = {pz_pct_l3:.1f}%  (paper target ~67%)")

    pz["is_layer2"] = pz["layer2_deficit"] > 0
    pz["is_layer3"] = pz["layer3_deficit"] > 0
    pzL2 = set(pz.loc[pz["is_layer2"], "locus_tag"].tolist())
    pzL3 = set(pz.loc[pz["is_layer3"], "locus_tag"].tolist())
    pzU = set(pz["locus_tag"].tolist())
    pz_only_l2 = pzL2 - pzL3
    pz_only_l3 = pzL3 - pzL2
    pz_both = pzL2 & pzL3
    pz_neither = pzU - (pzL2 | pzL3)
    pz_jacc = len(pz_both) / len(pzL2 | pzL3) if (pzL2 | pzL3) else 0.0
    pz_cont = np.array(
        [
            [len(pz_both), len(pz_only_l2)],
            [len(pz_only_l3), len(pz_neither)],
        ]
    )
    pz_or, pz_p = stats.fisher_exact(pz_cont, alternative="two-sided")
    pz_exp_both = (len(pzL2) / len(pzU)) * (len(pzL3) / len(pzU)) * len(pzU)
    print()
    print(f"|Layer 2| = {len(pzL2)}, |Layer 3| = {len(pzL3)}")
    print(f"|L2 only| = {len(pz_only_l2)}, |L3 only| = {len(pz_only_l3)}")
    print(f"|L2 ∩ L3| = {len(pz_both)}  (expected if indep ≈ {pz_exp_both:.1f})")
    print(f"Jaccard = {pz_jacc:.4f}  Fisher OR = {pz_or:.3f}, p = {pz_p:.3g}")

    # Append to summary
    pz_summary = pd.DataFrame(
        {
            "metric": [
                f"PZ_n_universe",
                f"PZ_n_layer2",
                f"PZ_n_layer3",
                f"PZ_n_intersection",
                f"PZ_jaccard",
                f"PZ_fisher_OR",
                f"PZ_fisher_p",
                f"PZ_expected_intersection_if_independent",
                f"PZ_layer2_pct_of_deficit",
                f"PZ_layer3_pct_of_deficit",
                f"PZ_window_half_bp",
            ],
            "value": [
                len(pzU),
                len(pzL2),
                len(pzL3),
                len(pz_both),
                round(pz_jacc, 4),
                round(pz_or, 4),
                pz_p,
                round(pz_exp_both, 2),
                round(pz_pct_l2, 2),
                round(pz_pct_l3, 2),
                PROTZONE_HALF,
            ],
        }
    )
    out_summary_full = pd.concat([out_summary, pz_summary], ignore_index=True)
    out_summary_full.to_csv(summary_path, sep="\t", index=False)
    print(f"\n[write] {summary_path} (with TSS-proximal block appended)")


if __name__ == "__main__":
    main()
