#!/usr/bin/env python
"""Extract category-A 20-gene protein sequences and Pfam-like product info from GenBank.

Inputs:
- tables/exposed_TF_notable_list.tsv  (rows where 'category' contains 'A')
- /Users/okaban/bioinfo/methyl/260102_M145/data/ref.gbk  (NCBI GenBank for GCF_000203835.1)

Outputs:
- data/categoryA_genes.tsv        (locus_tag, old_locus_tag, ..., protein_id, product, length)
- data/categoryA_proteins.faa     (FASTA of protein sequences, header = old_locus_tag|locus_tag)
- data/missing_categoryA.tsv      (rows we could not match)
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

from Bio import SeqIO

ROOT = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/57_temporal_dynamics_exposed_TF")
NOTABLE = ROOT / "tables" / "exposed_TF_notable_list.tsv"
GBK = Path("/Users/okaban/bioinfo/methyl/260102_M145/data/ref.gbk")
OUT_TSV = ROOT / "data" / "categoryA_genes.tsv"
OUT_FAA = ROOT / "data" / "categoryA_proteins.faa"
OUT_MISS = ROOT / "data" / "missing_categoryA.tsv"


def load_categoryA() -> list[dict]:
    rows = []
    with NOTABLE.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            cat = r.get("category", "")
            if "A" in cat:
                rows.append(r)
    return rows


def index_genbank(gbk_path: Path) -> dict:
    """Index CDS features by both new locus_tag (SC_RS...) and old_locus_tag (SCO...)."""
    idx = {}
    for rec in SeqIO.parse(str(gbk_path), "genbank"):
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            q = feat.qualifiers
            keys = []
            for k in ("locus_tag", "old_locus_tag"):
                for v in q.get(k, []):
                    keys.append(v)
            translation = q.get("translation", [None])[0]
            entry = {
                "locus_tag": q.get("locus_tag", [""])[0],
                "old_locus_tag": q.get("old_locus_tag", [""])[0],
                "protein_id": q.get("protein_id", [""])[0],
                "product": q.get("product", [""])[0],
                "translation": translation,
                "length": len(translation) if translation else 0,
            }
            for k in keys:
                idx[k] = entry
    return idx


def main() -> int:
    cat = load_categoryA()
    print(f"[info] category-A rows: {len(cat)}")
    idx = index_genbank(GBK)
    print(f"[info] indexed {len(idx)} locus_tag keys from GenBank")

    found, missing = [], []
    with OUT_FAA.open("w") as faa:
        for r in cat:
            new_lt = r["locus_tag"].strip()
            old_lt = r["old_locus_tag"].strip()
            entry = idx.get(new_lt) or idx.get(old_lt)
            if entry is None or not entry["translation"]:
                missing.append(r)
                continue
            header = f">{old_lt or new_lt}|{new_lt} {entry['product']}"
            faa.write(header + "\n")
            seq = entry["translation"]
            for i in range(0, len(seq), 60):
                faa.write(seq[i:i + 60] + "\n")
            r2 = dict(r)
            r2["protein_id"] = entry["protein_id"]
            r2["product_gbk"] = entry["product"]
            r2["protein_length"] = entry["length"]
            found.append(r2)

    if found:
        with OUT_TSV.open("w") as out:
            cols = list(found[0].keys())
            out.write("\t".join(cols) + "\n")
            for r in found:
                out.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    if missing:
        with OUT_MISS.open("w") as out:
            cols = list(missing[0].keys())
            out.write("\t".join(cols) + "\n")
            for r in missing:
                out.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    print(f"[done] proteins extracted: {len(found)} / missing: {len(missing)}")
    print(f"[out] {OUT_FAA}")
    print(f"[out] {OUT_TSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
