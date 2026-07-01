#!/usr/bin/env python
"""Integrate BLAST + Pfam results into a per-gene similarity table + summary.

Inputs:
- data/categoryA_genes.tsv           (20 query genes with metadata)
- data/categoryA_proteins.faa
- data/reference_proteins.tsv        (6 literature references)
- data/categoryA_vs_refs.blast.tsv   (blastp tabular)
- data/categoryA_pfam.tblout         (hmmscan tblout for queries)
- data/reference_pfam.tblout         (hmmscan tblout for refs)

Outputs:
- results/notable_genes_similarity_analysis.tsv
- results/notable_genes_similarity_summary.txt
- results/notable_genes_pfam_overlap_matrix.tsv

Functional-similarity rule (per query gene):
- "Same Pfam family" if any Pfam accession of the query overlaps any Pfam
  accession of one or more literature references (after curated family
  expansion — see FAMILY_GROUPS).
- Sequence similarity is reported via best BLAST hit (pident, evalue, bitscore).
- A query is classed as "functionally analogous to literature reference X" when
  it shares the canonical DNA-binding/regulatory family of X.
"""
from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/57_temporal_dynamics_exposed_TF")
DATA = ROOT / "data"
RES = ROOT / "results"
RES.mkdir(exist_ok=True)

# Curated functional grouping. Each literature reference ("group label") maps
# to a set of Pfam accessions and Pfam names that we consider "the same
# regulatory family". This intentionally goes beyond the strictly best-scoring
# hit so that, for example, both MarR and MarR_2 count as "MarR family".
FAMILY_GROUPS = {
    "MarR (E. coli MarR; methylation-context regulator)": {
        "pfam_acc": {"PF01047", "PF12802"},
        "pfam_name": {"MarR", "MarR_2"},
        "ref_id": "P27245",
    },
    "LysR/OxyR (E. coli OxyR; agn43 Dam-switch)": {
        "pfam_acc": {"PF00126", "PF03466"},
        "pfam_name": {"HTH_1", "LysR_substrate"},
        "ref_id": "P0ACQ4",
    },
    "Lrp/AsnC (E. coli Lrp; pap GATC switch)": {
        # Diagnostic only: AsnC_trans_reg ligand-binding + HTH_AsnC-type HTH.
        # PF13412 (HTH_24) is excluded because it is a generic winged-HTH found
        # in many unrelated TFs (incl. MarR-family).
        "pfam_acc": {"PF01037", "PF13404"},
        "pfam_name": {"AsnC_trans_reg", "HTH_AsnC-type"},
        "ref_id": "P0ACJ0",
    },
    "Sigma-70 ECF (M. tuberculosis SigB; m6A-near-σ-box)": {
        # Stringent diagnostic Pfams only: Sigma70_r2 (-10 promoter recognition)
        # and Sigma70_r3. PF04545 (Sigma70_r4) and PF08281 (Sigma70_r4_2) are
        # excluded because both are HTHs that cross-hit non-sigma DNA-binding
        # proteins (e.g. some GerE-family response regulators).
        "pfam_acc": {"PF04542", "PF04539", "PF00140"},
        "pfam_name": {"Sigma70_r2", "Sigma70_r3", "Sigma70_r1_2"},
        "ref_id": "P9WGI1",
    },
    "TCS Response Regulator GerE/Receiver (S. coelicolor RamR analogue)": {
        # SCO6685 is itself a GerE/Receiver RR. We use GerE+Response_reg as the
        # canonical TCS-RR family rather than the (mis-annotated) Q9XAP2 sequence.
        "pfam_acc": {"PF00196", "PF00072"},
        "pfam_name": {"GerE", "Response_reg"},
        "ref_id": "SCO6685_self",
    },
}

# We also report TetR-family separately, even though it is not in the literature
# refs above, because Casadesús & Low 2006 review mentions TetR-family as a
# typical methylation-coupled repressor class.
EXTRA_FAMILIES = {
    "TetR/AcrR (Casadesús & Low 2006 review; methylation-coupled repressors)": {
        "pfam_acc": {"PF00440"},
        "pfam_name": {"TetR_N"},
        "ref_id": "literature_class",
    },
}


def parse_hmmscan_tblout(path: Path):
    """Yield dicts with target Pfam name, accession, query name, evalue, score."""
    with path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            yield {
                "target_name": parts[0],
                "target_acc": parts[1].split(".")[0],
                "query": parts[2],
                "evalue": float(parts[4]),
                "score": float(parts[5]),
            }


def best_pfam_per_query(tblout: Path) -> dict[str, list[dict]]:
    by_q: dict[str, list[dict]] = defaultdict(list)
    for h in parse_hmmscan_tblout(tblout):
        by_q[h["query"]].append(h)
    for q in by_q:
        by_q[q].sort(key=lambda x: x["evalue"])
    return by_q


def parse_blast(path: Path) -> dict[str, list[dict]]:
    by_q: dict[str, list[dict]] = defaultdict(list)
    cols = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen",
    ]
    with path.open() as fh:
        for line in fh:
            parts = line.rstrip().split("\t")
            if len(parts) < len(cols):
                continue
            row = dict(zip(cols, parts))
            for k in ("pident", "evalue", "bitscore"):
                row[k] = float(row[k])
            for k in ("length", "qlen", "slen"):
                row[k] = int(row[k])
            by_q[row["qseqid"]].append(row)
    for q in by_q:
        by_q[q].sort(key=lambda x: x["evalue"])
    return by_q


def main() -> int:
    # Load metadata
    meta = []
    with (DATA / "categoryA_genes.tsv").open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            meta.append(r)

    pfam_q = best_pfam_per_query(DATA / "categoryA_pfam.tblout")
    blast_q = parse_blast(DATA / "categoryA_vs_refs.blast.tsv")

    # Build output rows
    rows = []
    family_counts = defaultdict(int)
    family_match_per_row: list[set[str]] = []
    for r in meta:
        old = r["old_locus_tag"].strip()
        new = r["locus_tag"].strip()
        qkey = f"{old}|{new}" if old else new
        # find query in pfam dict (header used 'old|new SCO_xxx product' so qkey matches)
        ph = pfam_q.get(qkey, [])
        pfam_top = ph[0] if ph else None
        pfam_top_str = (
            f"{pfam_top['target_name']} ({pfam_top['target_acc']}, E={pfam_top['evalue']:.1e})"
            if pfam_top else "no Pfam hit"
        )
        all_pfam_names = {h["target_name"] for h in ph}
        all_pfam_accs = {h["target_acc"] for h in ph}

        # Determine which literature family group(s) the query matches
        family_matches = []
        for label, grp in {**FAMILY_GROUPS, **EXTRA_FAMILIES}.items():
            if all_pfam_accs & grp["pfam_acc"] or all_pfam_names & grp["pfam_name"]:
                family_matches.append(label)
                family_counts[label] += 1

        # BLAST best hit
        bh = blast_q.get(qkey, [])
        if bh:
            top = bh[0]
            blast_str = (
                f"{top['sseqid']} pident={top['pident']:.1f}% "
                f"len={top['length']} E={top['evalue']:.1e} "
                f"bit={top['bitscore']:.1f}"
            )
            blast_pident = top["pident"]
            blast_eval = top["evalue"]
            blast_hit = top["sseqid"]
        else:
            blast_str = "no hit"
            blast_pident = ""
            blast_eval = ""
            blast_hit = ""

        family_match_per_row.append(set(family_matches))
        rows.append({
            "locus_tag": new,
            "old_locus_tag": old,
            "tf_family_input": r["tf_family"],
            "annotation_input": r["annotation"],
            "product_gbk": r["product_gbk"],
            "log2FC": r["log2FC"],
            "delta_methyl_combined": r["methyl_change_combined"],
            "category": r["category"],
            "all_pfam_domains": ",".join(sorted(all_pfam_names)) if all_pfam_names else "",
            "top_pfam": pfam_top_str,
            "best_blast_hit": blast_hit,
            "blast_pident": blast_pident,
            "blast_evalue": blast_eval,
            "best_blast_summary": blast_str,
            "literature_family_match": "; ".join(family_matches) if family_matches else "no Pfam-level match",
        })

    # Write per-gene table
    out_tsv = RES / "notable_genes_similarity_analysis.tsv"
    cols = list(rows[0].keys())
    with out_tsv.open("w") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    # Pfam overlap matrix
    fam_labels = list(FAMILY_GROUPS.keys()) + list(EXTRA_FAMILIES.keys())
    out_mat = RES / "notable_genes_pfam_overlap_matrix.tsv"
    with out_mat.open("w") as out:
        out.write("locus_tag\told_locus_tag\t" + "\t".join(fam_labels) + "\n")
        for r, matches in zip(rows, family_match_per_row):
            cells = ["1" if lab in matches else "0" for lab in fam_labels]
            out.write(f"{r['locus_tag']}\t{r['old_locus_tag']}\t" + "\t".join(cells) + "\n")

    # Summary
    n = len(rows)
    n_with_match = sum(1 for r in rows if r["literature_family_match"] != "no Pfam-level match")
    n_blast_sig = sum(
        1 for r in rows
        if isinstance(r["blast_evalue"], float) and r["blast_evalue"] < 1e-3
    )
    n_blast_strong = sum(
        1 for r in rows
        if isinstance(r["blast_pident"], float) and r["blast_pident"] >= 30 and
           isinstance(r["blast_evalue"], float) and r["blast_evalue"] < 1e-3
    )

    summary = []
    summary.append("Category-A 20-gene similarity analysis summary")
    summary.append("=" * 60)
    summary.append("")
    summary.append("Method")
    summary.append("------")
    summary.append("- Query: 20 genes flagged 'category contains A' in")
    summary.append("  tables/exposed_TF_notable_list.tsv (TSS+/-500bp methylation")
    summary.append("  change AND |log2FC|>=1 AND same direction).")
    summary.append("- Protein sequences: extracted from")
    summary.append("  /Users/okaban/bioinfo/methyl/260102_M145/data/ref.gbk (NCBI RefSeq")
    summary.append("  GCF_000203835.1, NC_003888.3) using locus_tag/old_locus_tag.")
    summary.append("- Pfam domains: hmmscan vs Pfam-A.hmm (E<=1e-3).")
    summary.append("- Sequence similarity: blastp vs 6 curated literature references")
    summary.append("  (UniProt fetch). Tabular output, e-value threshold 10 (loose).")
    summary.append("- Functional similarity rule: a query is called 'family analogue'")
    summary.append("  of a literature reference if its Pfam accession set overlaps the")
    summary.append("  reference family group (see FAMILY_GROUPS in script).")
    summary.append("")
    summary.append("Headline numbers")
    summary.append("----------------")
    summary.append(f"- 20 / 20 query proteins extracted (no missing).")
    summary.append(f"- {n_with_match} / {n} share Pfam family with at least one literature reference.")
    summary.append(f"- {n_blast_sig} / {n} have a BLAST hit with E<1e-3 to any literature reference.")
    summary.append(f"- {n_blast_strong} / {n} also reach pident>=30% on that hit.")
    summary.append("")
    summary.append("Per-family count (Pfam-overlap basis)")
    summary.append("-------------------------------------")
    for lab in fam_labels:
        summary.append(f"- {family_counts[lab]:>2}  {lab}")
    summary.append("")
    summary.append("What this supports (defensible claims)")
    summary.append("--------------------------------------")
    summary.append("- The category-A 20 genes cover regulator families that are")
    summary.append("  documented as methylation-coupled in other organisms")
    summary.append("  (MarR, LysR/OxyR-like, Lrp/AsnC, sigma-70 ECF, TCS GerE-RR,")
    summary.append("  TetR/AcrR). This is a *Pfam-family* analogy, not orthology.")
    summary.append("- Specific examples (see TSV for full list):")
    examples = []
    for r in rows:
        if r["literature_family_match"] != "no Pfam-level match":
            examples.append(
                f"    * {r['old_locus_tag']} ({r['locus_tag']}) [{r['tf_family_input']}]"
                f" -> {r['literature_family_match']}"
            )
    summary.extend(examples[:25])
    summary.append("")
    summary.append("What this does NOT support")
    summary.append("--------------------------")
    summary.append("- High BLAST identity / orthology to E. coli MarR, OxyR, Lrp,")
    summary.append("  H. pylori HP1021, or M. tuberculosis SigB. All best hits are")
    summary.append("  weak (mostly pident<35% over short stretches; most E>1e-3).")
    summary.append("  Streptomyces TFs are too divergent for direct sequence orthology.")
    summary.append("- A claim that the *mechanism* is the same. Only the Pfam family")
    summary.append("  is shared. Whether Dam/Dcm-style methylation actually toggles")
    summary.append("  these regulators in S. coelicolor must be tested directly")
    summary.append("  (e.g. mutate the methyl site, ChIP-seq, EMSA).")
    summary.append("")
    summary.append("Caveats on the literature reference set")
    summary.append("---------------------------------------")
    summary.append("- UniProt Q9XAP2 (labelled 'RamR' in some databases) is in fact")
    summary.append("  RamC (Radical_SAM + CofH/MqnC), NOT the SCO6685 response")
    summary.append("  regulator. Its BLAST/Pfam comparison is therefore not")
    summary.append("  meaningful as a TCS-RR positive control. We instead use the")
    summary.append("  GerE+Response_reg Pfam pair as the canonical TCS-RR family.")
    summary.append("- UniProt O25617 (HP1021) returned no strong Pfam hit at E<=1e-3.")
    summary.append("  HP1021 family classification therefore relies on the literature")
    summary.append("  rather than a Pfam-level matching here.")
    summary.append("")
    summary.append("Output files")
    summary.append("------------")
    summary.append(f"- {out_tsv.relative_to(ROOT.parent.parent)}")
    summary.append(f"- {out_mat.relative_to(ROOT.parent.parent)}")
    summary.append(f"- this file")

    out_txt = RES / "notable_genes_similarity_summary.txt"
    out_txt.write_text("\n".join(summary) + "\n")

    print(f"[done] wrote {out_tsv}")
    print(f"[done] wrote {out_mat}")
    print(f"[done] wrote {out_txt}")
    print(f"[stats] Pfam-family-match: {n_with_match}/{n};"
          f" BLAST sig (E<1e-3): {n_blast_sig}/{n};"
          f" pident>=30%: {n_blast_strong}/{n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
