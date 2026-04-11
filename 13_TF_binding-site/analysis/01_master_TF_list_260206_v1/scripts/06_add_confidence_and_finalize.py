#!/usr/bin/env python3
"""
06_add_confidence_and_finalize.py
master_TF_list_M145.tsv に confidence カラムを追加し、最終版を出力。

Confidence levels:
- high: Pfam TF domain + 少なくとも1つの文献ソース、またはσ因子
- medium: Pfam TF domain のみ、または複数の文献ソース
- low: Castro-Melchor のcistronメンバーのみ（Pfam TF domainなし）
"""
import pandas as pd

BASE_DIR = "/Users/okaban/bioinfo/rna-seq/13_TF_binding-site/analysis/01_master_TF_list_260206_v1"

def main():
    df = pd.read_csv(f"{BASE_DIR}/master_TF_list_M145.tsv", sep="\t")
    print(f"Input: {len(df)} TF candidates")

    def assign_confidence(row):
        sources = str(row["source_flags"]).split(";") if row["source_flags"] != "NA" else []
        has_pfam = "pfam" in sources
        has_za = "ZorroAranda" in sources
        has_cm = "CastroMelchor" in sources
        has_sigma = "SigmaReview" in sources
        n_sources = len(sources)

        if has_sigma:
            return "high"
        if has_pfam and (has_za or has_sigma):
            return "high"
        if has_za and has_cm:
            return "high"
        if has_pfam and has_cm:
            return "medium"
        if has_pfam:
            return "medium"
        if has_za:
            return "medium"
        if has_cm:
            # Castro-Melchor only - check if product suggests regulator
            product = str(row.get("product", "")).lower()
            reg_keywords = ["regulator", "repressor", "activator", "transcription",
                            "sigma", "dna-binding", "response receiver",
                            "two-component", "kinase sensor"]
            if any(kw in product for kw in reg_keywords):
                return "medium"
            return "low"
        return "low"

    df["confidence"] = df.apply(assign_confidence, axis=1)

    # Reorder columns
    out_cols = [
        "SCO_ID", "gene_id", "gene_name", "contig", "start", "end", "strand",
        "pfam_TF_flag", "tf_related_domains",
        "ZorroAranda_regulator_type", "CastroMelchor_regulator_flag",
        "sigma_flag", "sigma_group", "TF_family", "source_flags",
        "confidence", "product",
    ]
    df = df[out_cols]

    # Save final version
    out_path = f"{BASE_DIR}/master_TF_list_M145.tsv"
    df.to_csv(out_path, sep="\t", index=False)
    print(f"\nOutput: {out_path}")

    # Stats
    print(f"\n{'='*60}")
    print("CONFIDENCE LEVEL DISTRIBUTION")
    print(f"{'='*60}")
    print(df["confidence"].value_counts().to_string())

    print(f"\n--- High confidence TF family distribution ---")
    high = df[df["confidence"] == "high"]
    print(high["TF_family"].value_counts().to_string())

    print(f"\n--- Medium confidence TF family distribution ---")
    med = df[df["confidence"] == "medium"]
    print(med["TF_family"].value_counts().to_string())

    print(f"\n--- Low confidence summary ---")
    low = df[df["confidence"] == "low"]
    print(f"Low confidence count: {len(low)}")
    print(f"  With Pfam TF domain: {(low['pfam_TF_flag']=='yes').sum()}")
    print(f"  Without Pfam TF domain: {(low['pfam_TF_flag']=='no').sum()}")

    # Overall summary
    print(f"\n{'='*60}")
    print("FINAL SUMMARY")
    print(f"{'='*60}")
    print(f"Total TF/regulator candidates: {len(df)}")
    print(f"  High confidence: {len(high)}")
    print(f"  Medium confidence: {len(med)}")
    print(f"  Low confidence: {len(low)}")
    print(f"\nHigh + Medium (recommended working set): {len(high) + len(med)}")

if __name__ == "__main__":
    main()
