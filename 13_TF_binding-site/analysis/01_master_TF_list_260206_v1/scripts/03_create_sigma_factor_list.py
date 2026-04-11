#!/usr/bin/env python3
"""
03_create_sigma_factor_list.py
S. coelicolor A3(2) M145 のσ因子リストを作成。
情報源:
1. Sun et al. 2017 (Front Microbiol 8:2546) レビュー + 文献知識
2. Zorro-Aranda 2022 のsigma_factorフラグ
3. GFF product列からの"sigma factor"キーワードマッチ
4. Bentley et al. 2002 (Nature) ゲノム論文の情報

S. coelicolor A3(2) は Table 1 より 64個のσ因子を持つ:
- Group 1 (primary): 1 (hrdB/sigA = SCO5820)
- Group 2: 3 (hrdA = SCO3465, hrdC = SCO5683, hrdD = SCO2954 actually this needs checking)
- Group 3 (alternative): 10 (sigB, sigF, sigG, sigH, sigI, sigK, sigL, sigM, sigN, whiG)
- Group 4 (ECF): 50

主要な named sigma factors の SCO_ID は well-documented.
ECF σ因子はGFFからsigma関連product名で抽出する。
"""
import pandas as pd

OUT_DIR = "/Users/okaban/bioinfo/rna-seq/13_TF_binding-site/analysis/01_master_TF_list_260206_v1/intermediate"
GENE_INFO = f"{OUT_DIR}/M145_gene_basic_info.tsv"

# Well-characterized sigma factors from literature
# Sources: Sun et al. 2017, Bentley et al. 2002, Paget 2015
KNOWN_SIGMAS = [
    # Group 1 - Primary sigma factor
    {"SCO_ID": "SCO5820", "gene_name": "hrdB", "sigma_group": "sigma70_group1", "description": "Principal sigma factor (sigA/hrdB)"},
    # Group 2
    {"SCO_ID": "SCO3465", "gene_name": "hrdA", "sigma_group": "sigma70_group2", "description": "Group 2 sigma factor HrdA"},
    {"SCO_ID": "SCO5683", "gene_name": "hrdC", "sigma_group": "sigma70_group2", "description": "Group 2 sigma factor HrdC"},
    {"SCO_ID": "SCO0895", "gene_name": "hrdD", "sigma_group": "sigma70_group2", "description": "Group 2 sigma factor HrdD"},
    # Group 3 - Alternative sigma factors (sigB-like)
    {"SCO_ID": "SCO0600", "gene_name": "sigB", "sigma_group": "sigma70_group3", "description": "Alternative sigma factor SigB (osmotic/oxidative stress, differentiation)"},
    {"SCO_ID": "SCO1541", "gene_name": "sigF", "sigma_group": "sigma70_group3", "description": "Alternative sigma factor SigF (sporulation)"},
    {"SCO_ID": "SCO7341", "gene_name": "sigG", "sigma_group": "sigma70_group3", "description": "Alternative sigma factor SigG"},
    {"SCO_ID": "SCO5243", "gene_name": "sigH", "sigma_group": "sigma70_group3", "description": "Alternative sigma factor SigH (stress response)"},
    {"SCO_ID": "SCO3068", "gene_name": "sigI", "sigma_group": "sigma70_group3", "description": "Alternative sigma factor SigI"},
    {"SCO_ID": "SCO6520", "gene_name": "sigK", "sigma_group": "sigma70_group3", "description": "Alternative sigma factor SigK (differentiation)"},
    {"SCO_ID": "SCO7278", "gene_name": "sigL", "sigma_group": "sigma70_group3", "description": "Alternative sigma factor SigL (osmotic stress)"},
    {"SCO_ID": "SCO7314", "gene_name": "sigM", "sigma_group": "sigma70_group3", "description": "Alternative sigma factor SigM (osmotic stress)"},
    {"SCO_ID": "SCO4034", "gene_name": "sigN", "sigma_group": "sigma70_group3", "description": "Alternative sigma factor SigN (differentiation, stress)"},
    {"SCO_ID": "SCO5621", "gene_name": "whiG", "sigma_group": "sigma70_group3", "description": "Alternative sigma factor WhiG (sporulation)"},
    # Group 4 - ECF sigma factors (named)
    {"SCO_ID": "SCO5216", "gene_name": "sigR", "sigma_group": "ECF", "description": "ECF sigma factor SigR (thiol-oxidative stress)"},
    {"SCO_ID": "SCO3356", "gene_name": "sigE", "sigma_group": "ECF", "description": "ECF sigma factor SigE (cell wall stress)"},
    {"SCO_ID": "SCO2954", "gene_name": "sigU", "sigma_group": "ECF", "description": "ECF sigma factor SigU"},
    {"SCO_ID": "SCO4908", "gene_name": "bldN", "sigma_group": "ECF", "description": "ECF sigma factor BldN (aerial mycelium formation)"},
    {"SCO_ID": "SCO4005", "gene_name": "sigT", "sigma_group": "ECF", "description": "ECF sigma factor SigT (nitrogen stress, actinorhodin)"},
    {"SCO_ID": "SCO7328", "gene_name": "sigQ", "sigma_group": "ECF", "description": "ECF sigma factor SigQ"},
    {"SCO_ID": "SCO4895", "gene_name": "litS", "sigma_group": "ECF", "description": "ECF sigma factor LitS (light-induced carotenogenesis)"},
]

def main():
    # Load gene basic info
    df_genes = pd.read_csv(GENE_INFO, sep="\t")
    print(f"Total genes: {len(df_genes)}")

    # 1. Create known sigma factors dataframe
    df_known = pd.DataFrame(KNOWN_SIGMAS)
    df_known["source"] = "literature_curated"
    print(f"\nLiterature-curated named sigma factors: {len(df_known)}")

    # 2. Find additional sigma factors from GFF product annotations
    sigma_keywords = [
        "sigma factor",
        "RNA polymerase sigma",
        "sigma-70",
        "ECF sigma",
        "sigma70",
    ]

    mask = df_genes["product"].fillna("").str.lower().apply(
        lambda x: any(kw.lower() in x for kw in sigma_keywords)
    )
    df_gff_sigmas = df_genes[mask][["gene_id", "SCO_ID", "gene_name", "product"]].copy()
    print(f"\nGFF product-based sigma factor candidates: {len(df_gff_sigmas)}")

    # 3. Merge: combine known with GFF-based, preferring known annotations
    known_scos = set(df_known["SCO_ID"].values)

    # Additional sigmas from GFF not in known list
    additional_sigmas = []
    for _, row in df_gff_sigmas.iterrows():
        sco = row["SCO_ID"]
        if sco and sco not in known_scos:
            product = row["product"]
            # Determine sigma group from product
            if "ECF" in str(product) or "extracytoplasmic" in str(product).lower():
                sigma_group = "ECF"
            elif "sigma-70" in str(product).lower() or "sigma70" in str(product).lower():
                sigma_group = "sigma70_unclassified"
            else:
                sigma_group = "ECF"  # Most unnamed sigma factors in Streptomyces are ECF
            additional_sigmas.append({
                "SCO_ID": sco,
                "gene_name": row["gene_name"] if pd.notna(row["gene_name"]) else "",
                "sigma_group": sigma_group,
                "description": str(product),
                "source": "GFF_product_annotation",
            })
            known_scos.add(sco)

    df_additional = pd.DataFrame(additional_sigmas) if additional_sigmas else pd.DataFrame(
        columns=["SCO_ID", "gene_name", "sigma_group", "description", "source"]
    )
    print(f"Additional sigma factors from GFF: {len(df_additional)}")

    # 4. Combine all
    df_all_sigmas = pd.concat([df_known, df_additional], ignore_index=True)

    # 5. Map SCO_ID to gene_id
    sco_to_geneid = dict(zip(df_genes["SCO_ID"], df_genes["gene_id"]))
    df_all_sigmas["gene_id"] = df_all_sigmas["SCO_ID"].map(sco_to_geneid)

    # Check for unmapped
    unmapped = df_all_sigmas[df_all_sigmas["gene_id"].isna()]
    if len(unmapped) > 0:
        print(f"\nWARNING: {len(unmapped)} sigma factors could not be mapped to gene_id:")
        for _, row in unmapped.iterrows():
            print(f"  {row['SCO_ID']} ({row['gene_name']})")

    # Output
    out_path = f"{OUT_DIR}/SigmaFactors_M145_from_lit.tsv"
    df_all_sigmas.to_csv(out_path, sep="\t", index=False)
    print(f"\nOutput: {out_path}")
    print(f"Total sigma factors: {len(df_all_sigmas)}")
    print(f"\nSigma group distribution:")
    print(df_all_sigmas["sigma_group"].value_counts().to_string())
    print(f"\nSource distribution:")
    print(df_all_sigmas["source"].value_counts().to_string())

if __name__ == "__main__":
    main()
