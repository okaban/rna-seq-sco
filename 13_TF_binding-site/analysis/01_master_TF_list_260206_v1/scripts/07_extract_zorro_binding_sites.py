#!/usr/bin/env python3
"""
07_extract_zorro_binding_sites.py
Zorro-Aranda 2022 MOESM2 + MOESM3 から binding site 情報を抽出。
- MOESM3 MEME_BS: 推定binding site（配列、位置、p-value）
- MOESM2 Table 1: キュレーション済みTF→target相互作用（実験エビデンス分類付き）
"""
import pandas as pd
import openpyxl

BASE_DIR = "/Users/okaban/bioinfo/rna-seq/13_TF_binding-site/analysis/01_master_TF_list_260206_v1"
LIT_DIR = f"{BASE_DIR}/literature"
OUT_DIR = f"{BASE_DIR}/intermediate"
GENE_INFO = f"{OUT_DIR}/M145_gene_basic_info.tsv"

def main():
    # Load gene info for SCO_ID -> gene_id mapping
    df_genes = pd.read_csv(GENE_INFO, sep="\t")
    sco_to_geneid = {k: v for k, v in zip(df_genes["SCO_ID"], df_genes["gene_id"]) if pd.notna(k) and k != ""}
    geneid_to_sco = {v: k for k, v in sco_to_geneid.items()}

    # =============================================
    # 1. MOESM3 MEME_BS: Predicted binding sites
    # =============================================
    print("=" * 60)
    print("1. Parsing MOESM3 MEME_BS (predicted binding sites)")
    print("=" * 60)

    df_meme = pd.read_excel(
        f"{LIT_DIR}/Zorro-Aranda_2022_MOESM3_SupplementaryFile2.xlsx",
        sheet_name="MEME_BS"
    )
    print(f"MEME_BS rows: {len(df_meme)}")
    print(f"Columns: {list(df_meme.columns)}")
    print(df_meme.head())

    # Standardize column names (10 columns: 8 data + 2 extra)
    df_meme = df_meme.iloc[:, :8]  # Keep only first 8 columns
    df_meme.columns = ["TF_SCO", "TG_SCO", "site_sequence", "strand",
                       "left_pos", "right_pos", "p_value", "q_value"]

    # Map SCO_IDs to gene_ids
    df_meme["TF_gene_id"] = df_meme["TF_SCO"].map(sco_to_geneid)
    df_meme["TG_gene_id"] = df_meme["TG_SCO"].map(sco_to_geneid)

    # Add source
    df_meme["source"] = "ZorroAranda2022_MEME"
    df_meme["evidence_type"] = "computational_MEME"

    print(f"\nUnique TFs with binding sites: {df_meme['TF_SCO'].nunique()}")
    print(f"Unique TF-TG pairs: {len(df_meme.drop_duplicates(['TF_SCO','TG_SCO']))}")
    print(f"TF mapping success: {df_meme['TF_gene_id'].notna().sum()}/{len(df_meme)}")

    # Save
    meme_out = f"{OUT_DIR}/ZorroAranda2022_MEME_binding_sites.tsv"
    df_meme.to_csv(meme_out, sep="\t", index=False)
    print(f"Output: {meme_out}")

    # =============================================
    # 2. MOESM3 Inferred_BS: TF-TG pairs (BS-based inference)
    # =============================================
    print("\n" + "=" * 60)
    print("2. Parsing MOESM3 Inferred_BS")
    print("=" * 60)

    df_inferred_bs = pd.read_excel(
        f"{LIT_DIR}/Zorro-Aranda_2022_MOESM3_SupplementaryFile2.xlsx",
        sheet_name="Inferred_BS"
    )
    print(f"Inferred_BS rows: {len(df_inferred_bs)}")
    df_inferred_bs.columns = ["TF_SCO", "TG_SCO"]
    df_inferred_bs["TF_gene_id"] = df_inferred_bs["TF_SCO"].map(sco_to_geneid)
    df_inferred_bs["TG_gene_id"] = df_inferred_bs["TG_SCO"].map(sco_to_geneid)
    df_inferred_bs["source"] = "ZorroAranda2022_Inferred_BS"

    print(f"Unique TFs: {df_inferred_bs['TF_SCO'].nunique()}")
    inferred_bs_out = f"{OUT_DIR}/ZorroAranda2022_Inferred_BS_pairs.tsv"
    df_inferred_bs.to_csv(inferred_bs_out, sep="\t", index=False)
    print(f"Output: {inferred_bs_out}")

    # =============================================
    # 3. MOESM2 Table 1: Curated TF-TG interactions
    # =============================================
    print("\n" + "=" * 60)
    print("3. Parsing MOESM2 Table 1 (curated interactions)")
    print("=" * 60)

    df_curated = pd.read_excel(
        f"{LIT_DIR}/Zorro-Aranda_2022_MOESM2_SupplementaryTables.xlsx",
        sheet_name="Table 1",
        header=0
    )
    print(f"Table 1 rows: {len(df_curated)}")
    print(f"Columns: {list(df_curated.columns)}")

    # Clean up column names (first column is an index)
    cols = list(df_curated.columns)
    # The actual data columns based on our earlier inspection
    df_curated = df_curated.rename(columns={
        cols[0]: "idx",
        cols[1]: "TF_name",
        cols[2]: "TF_SCO",
        cols[3]: "TF_description",
        cols[4]: "TG_name",
        cols[5]: "TG_SCO",
        cols[6]: "experiment",
        cols[7]: "evidence",
        cols[8]: "regulatory_function",
        cols[9]: "evidence_classification",
    })

    # Drop rows where TF_SCO is NaN
    df_curated = df_curated.dropna(subset=["TF_SCO"])

    # Map to gene_ids
    df_curated["TF_gene_id"] = df_curated["TF_SCO"].map(sco_to_geneid)
    df_curated["TG_gene_id"] = df_curated["TG_SCO"].map(sco_to_geneid)
    df_curated["source"] = "ZorroAranda2022_Curated"

    print(f"\nCurated interactions: {len(df_curated)}")
    print(f"Unique TFs: {df_curated['TF_SCO'].nunique()}")
    print(f"Unique TGs: {df_curated['TG_SCO'].nunique()}")
    print(f"\nEvidence classification distribution:")
    print(df_curated["evidence_classification"].value_counts().to_string())
    print(f"\nRegulatory function distribution:")
    print(df_curated["regulatory_function"].value_counts().to_string())

    curated_out = f"{OUT_DIR}/ZorroAranda2022_curated_interactions.tsv"
    df_curated.to_csv(curated_out, sep="\t", index=False)
    print(f"Output: {curated_out}")

    # =============================================
    # 4. Summary: TFs with binding site information
    # =============================================
    print("\n" + "=" * 60)
    print("4. Summary")
    print("=" * 60)

    # Combine all TF SCO_IDs that have binding site info
    tfs_with_meme = set(df_meme["TF_SCO"].unique())
    tfs_with_inferred = set(df_inferred_bs["TF_SCO"].unique())
    tfs_curated = set(df_curated["TF_SCO"].unique())

    print(f"TFs with MEME binding sites: {len(tfs_with_meme)}")
    print(f"TFs with inferred BS pairs: {len(tfs_with_inferred)}")
    print(f"TFs in curated interactions: {len(tfs_curated)}")
    print(f"TFs with MEME AND curated: {len(tfs_with_meme & tfs_curated)}")
    print(f"All unique TFs with any BS info: {len(tfs_with_meme | tfs_with_inferred | tfs_curated)}")

if __name__ == "__main__":
    main()
