#!/usr/bin/env python3
"""
01_extract_gene_info.py
GFFファイルから全gene（CDS）の基本情報を抽出し M145_gene_basic_info.tsv を作成する。
gene_master_DESeq2.tsv のSCO_IDマッピングも活用。
"""
import pandas as pd
import re
import sys

GFF_PATH = "/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/genomic.gff"
GENE_MASTER = "/Users/okaban/bioinfo/rna-seq/05_annotation/analysis/05_annotation_260128_v1/tables/gene_master_DESeq2.tsv"
OUT_DIR = "/Users/okaban/bioinfo/rna-seq/13_TF_binding-site/analysis/01_master_TF_list_260206_v1/intermediate"

CONTIG_MAP = {
    "NC_003888.3": "chromosome",
    "NC_003903.1": "SCP1",
    "NC_003904.1": "SCP2",
}

def parse_gff_attributes(attr_str):
    attrs = {}
    for item in attr_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            attrs[k] = v
    return attrs

def main():
    # Parse GFF for gene features
    genes = []
    with open(GFF_PATH) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            if fields[2] != "gene":
                continue
            contig = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attrs = parse_gff_attributes(fields[8])
            locus_tag = attrs.get("locus_tag", "")
            gene_name = attrs.get("gene", "")
            old_locus_tag = attrs.get("old_locus_tag", "")
            gene_biotype = attrs.get("gene_biotype", "")
            genes.append({
                "gene_id": locus_tag,
                "old_locus_tag": old_locus_tag,
                "gene_name": gene_name,
                "gene_biotype": gene_biotype,
                "contig_raw": contig,
                "contig": CONTIG_MAP.get(contig, contig),
                "start": start,
                "end": end,
                "strand": strand,
            })

    df_gff = pd.DataFrame(genes)
    print(f"GFF gene entries: {len(df_gff)}")
    print(f"  protein_coding: {(df_gff['gene_biotype']=='protein_coding').sum()}")
    print(f"  pseudogene: {(df_gff['gene_biotype']=='pseudogene').sum()}")
    print(f"  other: {(~df_gff['gene_biotype'].isin(['protein_coding','pseudogene'])).sum()}")

    # Also load gene_master for SCO_ID mapping enrichment
    gm = pd.read_csv(GENE_MASTER, sep="\t", usecols=["gene_id", "old_locus_tag", "gene_name", "product", "protein_id"])
    gm = gm.rename(columns={"old_locus_tag": "gm_old_locus_tag", "gene_name": "gm_gene_name", "product": "gm_product", "protein_id": "gm_protein_id"})

    df = df_gff.merge(gm, on="gene_id", how="left")

    # Fill in old_locus_tag and gene_name from gene_master if missing
    df["old_locus_tag"] = df["old_locus_tag"].where(df["old_locus_tag"] != "", df["gm_old_locus_tag"])
    df["gene_name"] = df["gene_name"].where(df["gene_name"] != "", df["gm_gene_name"])
    df["product"] = df["gm_product"]
    df["protein_id"] = df["gm_protein_id"]

    # Create SCO_ID column (use old_locus_tag as SCO_ID)
    df["SCO_ID"] = df["old_locus_tag"].fillna("")

    # Select and output
    out_cols = ["gene_id", "SCO_ID", "gene_name", "gene_biotype", "contig", "start", "end", "strand", "product", "protein_id"]
    df_out = df[out_cols].copy()

    out_path = f"{OUT_DIR}/M145_gene_basic_info.tsv"
    df_out.to_csv(out_path, sep="\t", index=False)
    print(f"\nOutput: {out_path}")
    print(f"Total genes: {len(df_out)}")
    print(f"With SCO_ID: {(df_out['SCO_ID']!='').sum()}")
    print(f"With gene_name: {df_out['gene_name'].notna().sum() & (df_out['gene_name']!='').sum()}")
    print(f"\nContig distribution:")
    print(df_out["contig"].value_counts().to_string())

if __name__ == "__main__":
    main()
