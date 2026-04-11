#!/usr/bin/env python3
"""
04_standardize_zorro_aranda.py
Zorro-Aranda 2022 のレギュレーターデータを gene_id (SC_RS#####) に変換し標準化。
"""
import pandas as pd

OUT_DIR = "/Users/okaban/bioinfo/rna-seq/13_TF_binding-site/analysis/01_master_TF_list_260206_v1/intermediate"
LIT_DIR = "/Users/okaban/bioinfo/rna-seq/13_TF_binding-site/analysis/01_master_TF_list_260206_v1/literature"
GENE_INFO = f"{OUT_DIR}/M145_gene_basic_info.tsv"

def main():
    # Load gene info for SCO_ID -> gene_id mapping
    df_genes = pd.read_csv(GENE_INFO, sep="\t")
    sco_to_geneid = dict(zip(df_genes["SCO_ID"], df_genes["gene_id"]))
    print(f"Gene mapping entries: {len(sco_to_geneid)}")

    # Load Zorro-Aranda regulators
    df_za = pd.read_csv(f"{LIT_DIR}/Zorro-Aranda_2022_regulators.tsv", sep="\t")
    print(f"Zorro-Aranda regulators: {len(df_za)}")
    print(f"Columns: {list(df_za.columns)}")
    print(df_za.head())

    # Map SCO_ID to gene_id
    df_za["gene_id"] = df_za["SCO_ID"].map(sco_to_geneid)

    # Check mapping
    mapped = df_za["gene_id"].notna().sum()
    unmapped = df_za["gene_id"].isna().sum()
    print(f"\nMapped: {mapped}, Unmapped: {unmapped}")
    if unmapped > 0:
        print("Unmapped SCO_IDs:")
        for _, row in df_za[df_za["gene_id"].isna()].iterrows():
            print(f"  {row['SCO_ID']} ({row.get('gene_name', 'NA')})")

    # Standardize regulator_type for the master table
    def classify_regulator_type(reg_type):
        """Map regulator_type to standardized categories."""
        if pd.isna(reg_type):
            return "TF"
        reg_type_lower = str(reg_type).lower()
        if "global" in reg_type_lower:
            return "global"
        elif "sigma" in reg_type_lower and "anti" not in reg_type_lower:
            return "sigma"
        elif "tcs" in reg_type_lower or "two-component" in reg_type_lower or "response_regulator" in reg_type_lower:
            return "TCS"
        else:
            return "TF"

    df_za["ZorroAranda_regulator_type"] = df_za["regulator_type"].apply(classify_regulator_type)

    # Output standardized version
    out_cols = ["SCO_ID", "gene_id", "gene_name", "regulator_type", "ZorroAranda_regulator_type", "is_global_regulator", "source"]
    df_out = df_za[out_cols].copy()
    out_path = f"{OUT_DIR}/ZorroAranda2022_regulators_standardized.tsv"
    df_out.to_csv(out_path, sep="\t", index=False)
    print(f"\nOutput: {out_path}")
    print(f"\nStandardized regulator type distribution:")
    print(df_out["ZorroAranda_regulator_type"].value_counts().to_string())

if __name__ == "__main__":
    main()
