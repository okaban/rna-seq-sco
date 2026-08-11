#!/usr/bin/env python3
"""Build accession -> (organismName, strain) map from assembly_data_report.jsonl."""
import json, csv, sys
src = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/data/ncbi_genomes_all_species/extracted/ncbi_dataset/data/assembly_data_report.jsonl"
out = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/73_SC_RS17645/acc_to_organism.tsv"

with open(src) as fin, open(out, "w") as fout:
    w = csv.writer(fout, delimiter="\t")
    w.writerow(["accession", "organism_name", "strain"])
    n = 0
    for line in fin:
        d = json.loads(line)
        acc = d.get("accession", "")
        org = d.get("organism", {}).get("organismName", "")
        strain = d.get("organism", {}).get("infraspecificNames", {}).get("strain", "")
        w.writerow([acc, org, strain])
        n += 1
print(f"Wrote {n} rows to {out}")
