#!/usr/bin/env python3
"""V-Defense (COG V) x 4mC-hypomethylated gene composition + GCCGGC promoter scan.

Two confirmation analyses tied to the COG-V Defense enrichment in 4mC
hypomethylated genes (OR ~ 2.3, padj ~ 0.03 reported in upstream COG enrichment).

Outputs (TSV):
  - v_defense_hypo_genes.tsv
  - v_defense_promoter_GCCGGC.tsv
  - v_defense_pgl_restriction_check.tsv (helper: explicit Pgl/BREX + Mcr* check)
"""

from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd

PROJECT_ROOT = Path("/Users/okaban/bioinfo/rna-seq")
RESULTS_DIR = (
    PROJECT_ROOT
    / "11_epigenome_integration"
    / "analysis"
    / "57_temporal_dynamics_exposed_TF"
    / "results"
)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

INTEGRATED_METHYL = (
    PROJECT_ROOT
    / "11_epigenome_integration"
    / "analysis"
    / "01_integration"
    / "integrated_methyl_expression_weighted.csv"
)
COG_CLASS = (
    PROJECT_ROOT
    / "12_supplementary_figures"
    / "analysis"
    / "12_supplementary_260202_v1"
    / "tables"
    / "gene_COG_classification.tsv"
)
GENE_ANNOT = (
    PROJECT_ROOT
    / "05_annotation"
    / "analysis"
    / "05_annotation_260128_v1"
    / "tables"
    / "gene_annotation_basic.tsv"
)
TSS_TABLE = (
    PROJECT_ROOT
    / "11_epigenome_integration"
    / "analysis"
    / "18_tss_analyses"
    / "comprehensive_tss_table.csv"
)
GENOME_FASTA = (
    PROJECT_ROOT
    / "11_epigenome_integration"
    / "data"
    / "NC_003888.3.fna"
)

DMG_THRESHOLD = 10.0  # percent methylation change threshold (matches existing pipeline)
PROMOTER_UP = 200     # bp upstream of TSS
PROMOTER_DOWN = 50    # bp downstream of TSS
MOTIF = "GCCGGC"      # palindromic 6-mer; reverse complement is itself


def load_genome_fasta(path: Path) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    current: str | None = None
    parts: List[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current is not None:
                    seqs[current] = "".join(parts).upper()
                current = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if current is not None:
        seqs[current] = "".join(parts).upper()
    return seqs


def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(comp)[::-1]


def find_motif_positions(seq: str, motif: str) -> List[int]:
    positions: List[int] = []
    start = 0
    while True:
        i = seq.find(motif, start)
        if i == -1:
            break
        positions.append(i)
        start = i + 1
    return positions


def main() -> None:
    integrated = pd.read_csv(INTEGRATED_METHYL)
    cog = pd.read_csv(COG_CLASS, sep="\t")
    annot = pd.read_csv(GENE_ANNOT, sep="\t")
    tss = pd.read_csv(TSS_TABLE)
    genome = load_genome_fasta(GENOME_FASTA)
    chrom = "NC_003888.3"
    chrom_seq = genome[chrom]

    # ---- Step 1: hypomethylated 4mC gene set across timepoints --------------
    change_cols = [
        "4mC_change_T2_vs_T1",
        "4mC_change_T3_vs_T1",
        "4mC_change_T3_vs_T2",
    ]
    expr_cols = ["log2FC_T2_vs_T1", "log2FC_T3_vs_T1", "log2FC_T3_vs_T2"]

    methyl = integrated[["gene_id"] + change_cols + expr_cols].copy()
    methyl[change_cols] = methyl[change_cols].fillna(0.0)
    methyl["min_4mC_change"] = methyl[change_cols].min(axis=1)
    methyl["is_4mC_hypo_any"] = methyl[change_cols].lt(-DMG_THRESHOLD).any(axis=1)
    hypo_set = set(methyl.loc[methyl["is_4mC_hypo_any"], "gene_id"])
    print(f"4mC-hypomethylated genes (any pairwise comparison): {len(hypo_set)}")

    # ---- Step 2: V-Defense gene set ----------------------------------------
    cog["is_V_defense"] = cog["COG_category"].str.startswith("V")
    v_set = set(cog.loc[cog["is_V_defense"], "gene_id"])
    print(f"COG V-Defense genes: {len(v_set)}")

    # ---- Step 3: intersection ---------------------------------------------
    v_hypo = sorted(v_set & hypo_set)
    print(f"V-Defense ∩ 4mC-hypomethylated: {len(v_hypo)}")

    annot_lookup = annot.set_index("gene_id")[
        ["old_locus_tag", "gene_name", "product", "contig", "start", "end", "strand"]
    ]
    cog_lookup = cog.set_index("gene_id")["COG_category"]

    rows = []
    for gid in v_hypo:
        ann = annot_lookup.loc[gid] if gid in annot_lookup.index else None
        m = methyl.loc[methyl["gene_id"] == gid].iloc[0]
        rows.append(
            {
                "gene_id": gid,
                "old_locus_tag": ann["old_locus_tag"] if ann is not None else "",
                "gene_name": ann["gene_name"] if ann is not None else "",
                "product": ann["product"] if ann is not None else "",
                "COG_category": cog_lookup.get(gid, ""),
                "4mC_change_T2_vs_T1": m["4mC_change_T2_vs_T1"],
                "4mC_change_T3_vs_T1": m["4mC_change_T3_vs_T1"],
                "4mC_change_T3_vs_T2": m["4mC_change_T3_vs_T2"],
                "min_4mC_change": m["min_4mC_change"],
                "log2FC_T2_vs_T1": m["log2FC_T2_vs_T1"],
                "log2FC_T3_vs_T1": m["log2FC_T3_vs_T1"],
                "log2FC_T3_vs_T2": m["log2FC_T3_vs_T2"],
            }
        )
    out1 = pd.DataFrame(rows)
    # Stable ordering by largest hypomethylation magnitude (most negative first)
    out1 = out1.sort_values("min_4mC_change").reset_index(drop=True)
    out1_path = RESULTS_DIR / "v_defense_hypo_genes.tsv"
    out1.to_csv(out1_path, sep="\t", index=False)
    print(f"Saved: {out1_path}  (rows={len(out1)})")

    # ---- Step 3b: explicit Pgl/BREX + restriction nuclease check -----------
    canonical_targets = {
        # BREX / Pgl operon (from gene_annotation_basic.tsv)
        "SC_RS35330": "pglW (SCO6626)",
        "SC_RS35335": "pglX (SCO6627)",
        "SC_RS35370": "pglY (SCO6635)",
        "SC_RS35375": "pglZ (SCO6636)",
        "SC_RS28835": "pglX (SCO5331, BREX-2 MTase)",
        # User-mentioned "SCO1444-SCO1447" region (note: these are NOT pgl in M145 NCBI annotation)
        "SC_RS09195": "SCO1444 (chitinase domain)",
        "SC_RS09200": "SCO1445 (SDR oxidoreductase)",
        "SC_RS09205": "SCO1446 (hypothetical)",
        "SC_RS09210": "SCO1447 (ROK family TF)",
        # Methyl-specific restriction nucleases
        "SC_RS23250": "SCO4213 (restriction endonuclease)",
        "SC_RS25315": "SCO4631 (mcrA, type IV restriction endonuclease)",
        "SC_RS16410": "SCO2863 (DUF3427)",
    }
    chk_rows = []
    for gid, label in canonical_targets.items():
        in_v_defense = gid in v_set
        in_hypo = gid in hypo_set
        in_v_hypo = in_v_defense and in_hypo
        m_row = methyl.loc[methyl["gene_id"] == gid]
        if len(m_row):
            mr = m_row.iloc[0]
            chk_rows.append(
                {
                    "gene_id": gid,
                    "label": label,
                    "in_COG_V_Defense": in_v_defense,
                    "is_4mC_hypo_any": in_hypo,
                    "in_V_AND_hypo": in_v_hypo,
                    "4mC_change_T2_vs_T1": mr["4mC_change_T2_vs_T1"],
                    "4mC_change_T3_vs_T1": mr["4mC_change_T3_vs_T1"],
                    "4mC_change_T3_vs_T2": mr["4mC_change_T3_vs_T2"],
                    "min_4mC_change": mr["min_4mC_change"],
                    "log2FC_T3_vs_T1": mr["log2FC_T3_vs_T1"],
                }
            )
        else:
            chk_rows.append(
                {
                    "gene_id": gid,
                    "label": label,
                    "in_COG_V_Defense": in_v_defense,
                    "is_4mC_hypo_any": False,
                    "in_V_AND_hypo": False,
                    "4mC_change_T2_vs_T1": None,
                    "4mC_change_T3_vs_T1": None,
                    "4mC_change_T3_vs_T2": None,
                    "min_4mC_change": None,
                    "log2FC_T3_vs_T1": None,
                }
            )
    chk_df = pd.DataFrame(chk_rows)
    chk_path = RESULTS_DIR / "v_defense_pgl_restriction_check.tsv"
    chk_df.to_csv(chk_path, sep="\t", index=False)
    print(f"Saved: {chk_path}")

    # ---- Step 4: GCCGGC promoter scan --------------------------------------
    tss_lookup = tss.set_index("gene_id")[["tss", "strand", "tss_source"]]
    motif_rc = reverse_complement(MOTIF)  # equals MOTIF (palindrome)
    promoter_rows: List[Dict] = []
    for gid in v_hypo:
        ann = annot_lookup.loc[gid] if gid in annot_lookup.index else None
        gene_name_display = ""
        if ann is not None:
            gene_name_display = (
                ann["old_locus_tag"] or ann["gene_name"] or ""
            )
        if gid not in tss_lookup.index:
            promoter_rows.append(
                {
                    "gene_id": gid,
                    "gene_name": gene_name_display,
                    "TSS": None,
                    "strand": ann["strand"] if ann is not None else "",
                    "tss_source": "MISSING",
                    "promoter_start_genome": None,
                    "promoter_end_genome": None,
                    "GCCGGC_in_promoter": False,
                    "GCCGGC_count": 0,
                    "GCCGGC_positions_relative_to_TSS": "",
                }
            )
            continue
        t_row = tss_lookup.loc[gid]
        tss_pos = int(t_row["tss"])  # 1-based per GFF convention
        strand = t_row["strand"]
        tss_source = t_row["tss_source"]

        if strand == "+":
            prom_start = tss_pos - PROMOTER_UP   # genomic 1-based inclusive
            prom_end = tss_pos + PROMOTER_DOWN
        else:
            prom_start = tss_pos - PROMOTER_DOWN
            prom_end = tss_pos + PROMOTER_UP

        # clamp
        prom_start = max(1, prom_start)
        prom_end = min(len(chrom_seq), prom_end)

        # 0-based slicing
        seq_plus = chrom_seq[prom_start - 1 : prom_end]
        if strand == "+":
            promoter_seq = seq_plus
        else:
            promoter_seq = reverse_complement(seq_plus)

        # search both strands (palindromic, but generic)
        positions_fwd = find_motif_positions(promoter_seq, MOTIF)
        positions_rev = (
            []
            if MOTIF == motif_rc
            else find_motif_positions(promoter_seq, motif_rc)
        )
        # Promoter coordinate frame: index 0 corresponds to TSS - PROMOTER_UP
        # so position relative to TSS = i - PROMOTER_UP
        rel = sorted(set(p - PROMOTER_UP for p in positions_fwd + positions_rev))
        rel_str = ";".join(str(r) for r in rel)
        promoter_rows.append(
            {
                "gene_id": gid,
                "gene_name": gene_name_display,
                "TSS": tss_pos,
                "strand": strand,
                "tss_source": tss_source,
                "promoter_start_genome": prom_start,
                "promoter_end_genome": prom_end,
                "GCCGGC_in_promoter": bool(rel),
                "GCCGGC_count": len(rel),
                "GCCGGC_positions_relative_to_TSS": rel_str,
            }
        )
    out2 = pd.DataFrame(promoter_rows)
    out2_path = RESULTS_DIR / "v_defense_promoter_GCCGGC.tsv"
    out2.to_csv(out2_path, sep="\t", index=False)
    print(f"Saved: {out2_path}  (rows={len(out2)})")

    # ---- Console summary ----------------------------------------------------
    n_with_motif = int(out2["GCCGGC_in_promoter"].sum())
    n_total = len(out2)
    print()
    print("=== Summary ===")
    print(f"  V-Defense ∩ 4mC-hypo: {n_total}")
    print(
        f"  GCCGGC in promoter (-{PROMOTER_UP}/+{PROMOTER_DOWN}): "
        f"{n_with_motif} / {n_total} "
        f"({100.0 * n_with_motif / n_total:.1f}%)"
    )
    pgl_hits = chk_df[chk_df["in_V_AND_hypo"]]
    print(f"  Pgl/BREX/restriction-nuclease canonical targets in V∩hypo set: "
          f"{len(pgl_hits)}")
    if len(pgl_hits):
        for _, r in pgl_hits.iterrows():
            print(f"    - {r['gene_id']}  ({r['label']})")


if __name__ == "__main__":
    main()
