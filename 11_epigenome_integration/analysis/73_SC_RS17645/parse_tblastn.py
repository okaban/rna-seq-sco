#!/usr/bin/env python3
"""
Parse tblastn output (outfmt 6 with extended columns) for SC_RS17645.
Aggregate per-genome best hit, classify conservation, and emit:
  - blastp_streptomyces_hits.tsv  (top 50 best per genome, accession-level)
  - tblastn_conservation_summary.txt  (counts at thresholds)
"""
from __future__ import annotations
import csv
from collections import defaultdict
from pathlib import Path

OUT = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/73_SC_RS17645")
TBLASTN = OUT / "tblastn_raw.tsv"
ACC_MAP = OUT / "acc_to_organism.tsv"
TOP_TSV = OUT / "blastp_streptomyces_hits.tsv"
SUMMARY = OUT / "tblastn_conservation_summary.txt"

QUERY_LEN = 679  # SC_RS17645 protein length (aa)

# Outfmt: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs
def main() -> None:
    acc2org: dict[str, tuple[str, str]] = {}
    with ACC_MAP.open() as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                acc2org[parts[0]] = (parts[1], parts[2])
            elif len(parts) == 2:
                acc2org[parts[0]] = (parts[1], "")

    # Best hit per genome (max bitscore)
    best: dict[str, dict] = {}
    with TBLASTN.open() as f:
        for line in f:
            f_ = line.rstrip("\n").split("\t")
            if len(f_) < 12:
                continue
            sseqid = f_[1]
            # sseqid format: GCF_xxx.x|<contig>
            acc = sseqid.split("|", 1)[0] if "|" in sseqid else sseqid
            pident = float(f_[2])
            qlen_aln = int(f_[3])  # alignment length (aa, includes gaps)
            qstart = int(f_[6]); qend = int(f_[7])
            evalue = float(f_[10])
            bitscore = float(f_[11])
            # True query coverage = span of query covered by HSP / query length
            qcovhsp = (qend - qstart + 1) / QUERY_LEN * 100.0
            cur = best.get(acc)
            if cur is None or bitscore > cur["bitscore"]:
                best[acc] = {
                    "accession": acc,
                    "subject_contig": sseqid.split("|", 1)[1] if "|" in sseqid else sseqid,
                    "pident": pident,
                    "aln_len_aa": qlen_aln,
                    "qcov_pct": round(qcovhsp, 1),
                    "qstart": qstart,
                    "qend": qend,
                    "evalue": evalue,
                    "bitscore": bitscore,
                }

    # Aggregate HSP coverage per genome (sum of unique query positions, capped at QUERY_LEN)
    # Simpler: take max single-HSP coverage as primary; also report best-hit HSP
    # For final TSV, attach organism info
    rows = []
    for acc, h in best.items():
        org, strain = acc2org.get(acc, ("", ""))
        rows.append({
            "accession": acc,
            "organism": org,
            "strain": strain,
            "pct_identity": round(h["pident"], 2),
            "qcov_pct": h["qcov_pct"],
            "aln_len_aa": h["aln_len_aa"],
            "evalue": h["evalue"],
            "bitscore": round(h["bitscore"], 1),
        })

    # Sort by bitscore desc
    rows.sort(key=lambda r: r["bitscore"], reverse=True)

    # Write top 50 (accession-level best)
    cols = ["accession", "organism", "strain", "pct_identity", "qcov_pct",
            "aln_len_aa", "evalue", "bitscore"]
    with TOP_TSV.open("w") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows[:50]:
            w.writerow(r)

    # Conservation classification at multiple thresholds
    n_total = len(rows)
    def count(min_id: float, min_cov: float) -> int:
        return sum(1 for r in rows if r["pct_identity"] >= min_id and r["qcov_pct"] >= min_cov)

    # Also count distinct species (first 2 tokens of organism name)
    def species(name: str) -> str:
        toks = name.split()
        return " ".join(toks[:2]) if len(toks) >= 2 else name

    species_with_hit = defaultdict(list)
    for r in rows:
        species_with_hit[species(r["organism"])].append(r)

    with SUMMARY.open("w") as f:
        f.write("SC_RS17645 conservation across 833 local Streptomyces genomes\n")
        f.write("=" * 64 + "\n\n")
        f.write(f"Total genomes with any tblastn hit:                  {n_total} / 833\n")
        f.write(f"Genomes with hit (identity >=30%, qcov >=60%):       {count(30, 60)} / 833\n")
        f.write(f"Genomes with hit (identity >=50%, qcov >=80%):       {count(50, 80)} / 833\n")
        f.write(f"Genomes with hit (identity >=80%, qcov >=80%):       {count(80, 80)} / 833\n")
        f.write(f"Genomes with hit (identity >=95%, qcov >=95%):       {count(95, 95)} / 833\n")
        f.write(f"\nDistinct species with any hit:                       {len(species_with_hit)}\n")

        # Highest-identity per distinct species (top 30)
        top_species = []
        for sp, hits in species_with_hit.items():
            best_hit = max(hits, key=lambda r: r["bitscore"])
            top_species.append((sp, best_hit))
        top_species.sort(key=lambda x: x[1]["bitscore"], reverse=True)
        f.write(f"\nTop 30 species by best-hit bitscore:\n")
        for sp, h in top_species[:30]:
            f.write(f"  {sp:<35s} ident={h['pct_identity']:5.1f}%  qcov={h['qcov_pct']:5.1f}%  bits={h['bitscore']:6.1f}  ({h['accession']})\n")

    print(f"[parse_tblastn] {n_total} genomes had a hit; wrote {TOP_TSV} and {SUMMARY}")
    # Echo summary to stdout
    print(SUMMARY.read_text())


if __name__ == "__main__":
    main()
