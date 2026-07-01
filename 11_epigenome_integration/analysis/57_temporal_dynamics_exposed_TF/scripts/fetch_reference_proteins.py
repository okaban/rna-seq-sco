#!/usr/bin/env python
"""Fetch literature reference proteins from UniProt for similarity comparison.

References (curated for methylation-coupled regulation literature):
- E. coli K-12 MarR        UniProt P27245   (Casadesús/Low; classic MarR family)
- E. coli K-12 OxyR        UniProt P0ACQ4   (Wallecha/Hale 2002 — agn43 Dam-methylation switch)
- E. coli K-12 Lrp         UniProt P0ACJ0   (Braaten/Low — pap GATC methylation memory)
- H. pylori 26695 HP1021   UniProt O25617   (Kumar 2018; CGCG m4C-coupled regulator)
- M. tuberculosis SigB     UniProt P9WGI1   (Shell 2013 — m6A near σ-binding boxes; canonical SigB)
- S. coelicolor RamR       UniProt Q9XAP2   (TCS response regulator; positive control: same gene
                                            appears in our category-A as SCO6685, expected 100% self-hit)

Output:
- data/reference_proteins.faa  (FASTA)
- data/reference_proteins.tsv  (id, name, organism, role, source_ref)
"""
from __future__ import annotations

import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/57_temporal_dynamics_exposed_TF")
OUT_FAA = ROOT / "data" / "reference_proteins.faa"
OUT_TSV = ROOT / "data" / "reference_proteins.tsv"

REFS = [
    {
        "id": "P27245",
        "name": "MarR",
        "organism": "E. coli K-12",
        "role": "Multiple antibiotic resistance repressor (MarR family); classic Dam/methylation context",
        "ref": "Alekshun & Levy 1997; Casadesús & Low 2006 (Microbiol Mol Biol Rev)",
    },
    {
        "id": "P0ACQ4",
        "name": "OxyR",
        "organism": "E. coli K-12",
        "role": "agn43 Dam-methylation-sensitive transcription factor (LysR-family); methylation memory switch",
        "ref": "Wallecha et al. 2002 (J Bacteriol); Casadesús & Low 2006",
    },
    {
        "id": "P0ACJ0",
        "name": "Lrp",
        "organism": "E. coli K-12",
        "role": "Leucine-responsive regulatory protein; pap pyelonephritis-associated pilus GATC switch",
        "ref": "Braaten et al. 1994 (Cell); Casadesús & Low 2006",
    },
    {
        "id": "O25617",
        "name": "HP1021",
        "organism": "H. pylori 26695",
        "role": "m4C/m5C-context-coupled iron/redox regulator (RR-like)",
        "ref": "Kumar et al. 2018 (Nucleic Acids Res); Estibariz et al. 2019 (NAR)",
    },
    {
        "id": "P9WGI1",
        "name": "SigB",
        "organism": "M. tuberculosis H37Rv",
        "role": "Stress-response sigma factor; m6A enriched near sigma-binding boxes",
        "ref": "Shell et al. 2013 (PLoS Genet)",
    },
    {
        "id": "Q9XAP2",
        "name": "RamR",
        "organism": "S. coelicolor",
        "role": "Positive control / same gene as SCO6685 (response regulator of ramCSAB cluster)",
        "ref": "Keijser et al. 2002; Nguyen et al. 2002 (Mol Microbiol)",
    },
]


def fetch(uid: str) -> tuple[str, str] | None:
    url = f"https://rest.uniprot.org/uniprotkb/{uid}.fasta"
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            data = resp.read().decode()
        if not data.startswith(">"):
            return None
        header, *seqlines = data.strip().splitlines()
        seq = "".join(seqlines)
        return header, seq
    except Exception as e:
        print(f"[warn] fetch {uid} failed: {e}", file=sys.stderr)
        return None


def main() -> int:
    rows = []
    with OUT_FAA.open("w") as faa:
        for r in REFS:
            res = fetch(r["id"])
            if res is None:
                print(f"[warn] skipping {r['id']}", file=sys.stderr)
                continue
            header, seq = res
            short_header = f">{r['id']}|{r['name']}_{r['organism'].replace(' ', '_')}"
            faa.write(short_header + "\n")
            for i in range(0, len(seq), 60):
                faa.write(seq[i:i + 60] + "\n")
            r2 = dict(r)
            r2["length"] = len(seq)
            rows.append(r2)
            time.sleep(0.6)
    with OUT_TSV.open("w") as out:
        cols = ["id", "name", "organism", "role", "ref", "length"]
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    print(f"[done] reference proteins: {len(rows)}")
    print(f"[out] {OUT_FAA}")
    print(f"[out] {OUT_TSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
