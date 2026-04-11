#!/usr/bin/env python3
"""
05_parse_pfam_and_integrate.py
1. hmmscan domtblout結果をパースしてTF関連ドメインをフラグ付け
2. 全データソースを統合して master_TF_list_M145.tsv を作成

TF関連ドメインの包括的リスト（Pfam accession + domain name pattern）
"""
import pandas as pd
import re
import sys

BASE_DIR = "/Users/okaban/bioinfo/rna-seq/13_TF_binding-site/analysis/01_master_TF_list_260206_v1"
OUT_DIR = f"{BASE_DIR}/intermediate"
FINAL_DIR = BASE_DIR

GENE_INFO = f"{OUT_DIR}/M145_gene_basic_info.tsv"
DOMTBLOUT = f"{OUT_DIR}/M145_pfam_domtblout.txt"
ZORRO = f"{OUT_DIR}/ZorroAranda2022_regulators_standardized.tsv"
CASTRO = f"{OUT_DIR}/CastroMelchor2010_regulators_dedup.tsv"
SIGMA = f"{OUT_DIR}/SigmaFactors_M145_from_lit.tsv"

# --- TF-related Pfam domain patterns ---
# Comprehensive list based on literature and Pfam clan annotations
TF_DOMAIN_PATTERNS = {
    # HTH (Helix-Turn-Helix) superfamily
    "HTH": [
        "HTH_1", "HTH_3", "HTH_4", "HTH_5", "HTH_6", "HTH_7", "HTH_8",
        "HTH_9", "HTH_10", "HTH_11", "HTH_12", "HTH_13", "HTH_14", "HTH_15",
        "HTH_16", "HTH_17", "HTH_18", "HTH_19", "HTH_20", "HTH_21", "HTH_22",
        "HTH_23", "HTH_24", "HTH_25", "HTH_26", "HTH_27", "HTH_28", "HTH_29",
        "HTH_30", "HTH_31", "HTH_32", "HTH_33", "HTH_34", "HTH_35", "HTH_36",
        "HTH_37", "HTH_38", "HTH_39", "HTH_40", "HTH_41", "HTH_42", "HTH_43",
        "HTH_44", "HTH_45", "HTH_46", "HTH_47", "HTH_48", "HTH_49", "HTH_50",
        "HTH_AraC", "HTH_CodY", "HTH_Crp_2", "HTH_DeoR", "HTH_IclR",
        "HTH_Mga", "HTH_OrfB_IS605", "HTH_Tnp_1", "HTH_Tnp_4",
        "HTH_WhiA_N",
    ],
    # TF family-specific domains
    "TetR": ["TetR_N", "TetR_C", "TetR_C_2", "TetR_C_3", "TetR_C_4", "TetR_C_5",
             "TetR_C_6", "TetR_C_7", "TetR_C_8", "TetR_C_9", "TetR_C_10",
             "TetR_C_11", "TetR_C_12", "TetR_C_13"],
    "GntR": ["GntR", "FCD"],
    "MarR": ["MarR", "MarR_2"],
    "LacI": ["LacI", "Peripla_BP_1", "Peripla_BP_3", "Peripla_BP_4", "Peripla_BP_5", "Peripla_BP_6"],
    "LysR": ["LysR_substrate", "LysR"],
    "AraC": ["AraC_binding", "AraC_binding_2", "AraC_E_bind"],
    "MerR": ["MerR", "MerR_1", "MerR-HTH"],
    "IclR": ["IclR"],
    "DeoR": ["DeoR", "DeoRC"],
    "AsnC": ["AsnC_trans_reg"],
    "Crp_Fnr": ["Crp", "CAP_ED"],
    "LuxR": ["LuxR_C", "GerE", "Trans_reg_C"],
    "PadR": ["PadR"],
    "XRE": ["XRE_family", "Cro", "Phage_CI_repr"],
    # SARP and related
    "SARP": ["SARP", "HTH_SARP", "BTAD", "LAL"],
    # WhiB
    "WhiB": ["WhiB", "WhiB_2"],
    # Sigma factors
    "Sigma": [
        "Sigma70_r1_1", "Sigma70_r1_2", "Sigma70_r2", "Sigma70_r3",
        "Sigma70_r4", "Sigma70_r4_2", "Sigma70_ner", "Sigma70_ECF",
        "ECF_sigma", "Sigma54_activat", "Sigma54_CBD", "Sigma54_AID",
        "Sigma_E", "SigmaS",
    ],
    # Two-component systems
    "TCS": ["REC", "Response_reg", "Trans_reg_C", "GerE", "OmpR",
            "NarL", "LytTR", "CitB"],
    # Other DNA-binding / transcriptional regulation domains
    "Other_TF": [
        "Fur", "Fe_dep_repr_C", "CopG", "Arc",
        "MazE_antitoxin", "BolA", "FlhD", "FlhC",
        "OmpR", "NarL", "CitB", "LytTR",
        "SpoOA_C", "AlgZ_FimS", "PRD", "PRD_2",
        "Bac_DNA_binding", "DNA_binding_1", "DNA_binding_2", "DNA_binding_3",
        "DNA_binding_4", "Arg_repressor", "ArgR_DNA_bind",
        "GreA_GreB_N", "NusA_S1",
        "Mga", "RpiR", "BirA_HTH",
        "GerE", "ComK",
        "Rap_MphR",
    ],
}

# Flatten to a set for quick lookup
TF_DOMAINS_SET = set()
TF_DOMAIN_TO_FAMILY = {}
for family, domains in TF_DOMAIN_PATTERNS.items():
    for d in domains:
        TF_DOMAINS_SET.add(d)
        if d not in TF_DOMAIN_TO_FAMILY:
            TF_DOMAIN_TO_FAMILY[d] = family

def parse_domtblout(filepath):
    """Parse HMMER domtblout format."""
    records = []
    with open(filepath) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 23:
                continue
            records.append({
                "target_name": fields[0],       # Pfam domain name
                "target_acc": fields[1],         # Pfam accession
                "query_name": fields[3],         # protein_id
                "query_acc": fields[4],          # protein accession
                "full_evalue": float(fields[6]),
                "full_score": float(fields[7]),
                "dom_evalue": float(fields[12]),
                "dom_score": float(fields[13]),
                "ali_from": int(fields[17]),
                "ali_to": int(fields[18]),
                "env_from": int(fields[19]),
                "env_to": int(fields[20]),
                "description": " ".join(fields[22:]),
            })
    return pd.DataFrame(records)


def assign_tf_family(domains_list):
    """Assign TF family based on domain composition."""
    families = set()
    for d in domains_list:
        if d in TF_DOMAIN_TO_FAMILY:
            families.add(TF_DOMAIN_TO_FAMILY[d])

    # Priority-based assignment
    # Specific families take precedence over generic HTH
    specific_families = families - {"HTH", "Other_TF", "TCS"}

    if "SARP" in specific_families:
        return "SARP"
    elif "TetR" in specific_families:
        return "TetR"
    elif "GntR" in specific_families:
        return "GntR"
    elif "MarR" in specific_families:
        return "MarR"
    elif "LysR" in specific_families:
        return "LysR"
    elif "AraC" in specific_families:
        return "AraC"
    elif "LuxR" in specific_families:
        return "LuxR"
    elif "XRE" in specific_families:
        return "XRE"
    elif "WhiB" in specific_families:
        return "WhiB"
    elif "Sigma" in specific_families:
        return "Sigma"
    elif "MerR" in specific_families:
        return "MerR"
    elif "LacI" in specific_families:
        return "LacI"
    elif "IclR" in specific_families:
        return "IclR"
    elif "DeoR" in specific_families:
        return "DeoR"
    elif "AsnC" in specific_families:
        return "AsnC"
    elif "PadR" in specific_families:
        return "PadR"
    elif "Crp_Fnr" in specific_families:
        return "Crp/Fnr"
    elif "TCS" in families:
        return "TCS_response_regulator"
    elif "HTH" in families:
        return "HTH_other"
    elif "Other_TF" in families:
        return "Other_TF"
    else:
        return "Unknown_TF"


def main():
    # =========================================================
    # 1. Load gene basic info
    # =========================================================
    df_genes = pd.read_csv(GENE_INFO, sep="\t")
    print(f"Total genes: {len(df_genes)}")

    # protein_id -> gene_id mapping
    protid_to_geneid = dict(zip(df_genes["protein_id"].dropna(), df_genes["gene_id"]))
    sco_to_geneid = {k: v for k, v in zip(df_genes["SCO_ID"], df_genes["gene_id"]) if pd.notna(k) and k != ""}

    # =========================================================
    # 2. Parse Pfam domtblout
    # =========================================================
    print("\n--- Parsing Pfam domtblout ---")
    df_pfam = parse_domtblout(DOMTBLOUT)
    print(f"Total domain hits: {len(df_pfam)}")
    print(f"Unique proteins with hits: {df_pfam['query_name'].nunique()}")

    # Map protein_id -> gene_id
    df_pfam["gene_id"] = df_pfam["query_name"].map(protid_to_geneid)

    # Identify TF-related domains
    df_pfam["is_tf_domain"] = df_pfam["target_name"].isin(TF_DOMAINS_SET)

    # Also check by pattern matching for any HTH variants we might have missed
    def is_tf_by_pattern(name):
        name_upper = name.upper()
        if name_upper.startswith("HTH_"):
            return True
        if "SIGMA70" in name_upper or "ECF" in name_upper:
            return True
        if name in TF_DOMAINS_SET:
            return True
        return False

    df_pfam["is_tf_domain"] = df_pfam["target_name"].apply(is_tf_by_pattern)

    tf_hits = df_pfam[df_pfam["is_tf_domain"]]
    print(f"TF-related domain hits: {len(tf_hits)}")
    print(f"Unique proteins with TF domains: {tf_hits['query_name'].nunique()}")

    # =========================================================
    # 3. Aggregate per gene
    # =========================================================
    # All domains per gene
    gene_all_domains = df_pfam.groupby("gene_id").agg(
        all_domains=("target_name", lambda x: ";".join(sorted(set(x)))),
    ).reset_index()

    # TF domains per gene
    gene_tf_domains = tf_hits.groupby("gene_id").agg(
        tf_domains=("target_name", lambda x: ";".join(sorted(set(x)))),
    ).reset_index()

    # =========================================================
    # 4. Load literature data
    # =========================================================
    print("\n--- Loading literature data ---")

    # Zorro-Aranda
    df_za = pd.read_csv(ZORRO, sep="\t")
    za_geneid_map = dict(zip(df_za["gene_id"].dropna(), df_za["ZorroAranda_regulator_type"]))
    print(f"Zorro-Aranda regulators: {len(df_za)} ({len(za_geneid_map)} mapped)")

    # Castro-Melchor
    df_cm = pd.read_csv(CASTRO, sep="\t")
    # Need to map SCO_ID -> gene_id
    df_cm["gene_id"] = df_cm["SCO_ID"].map(sco_to_geneid)
    cm_geneids = set(df_cm["gene_id"].dropna())
    print(f"Castro-Melchor regulators: {len(df_cm)} (mapped: {len(cm_geneids)})")

    # Sigma factors
    df_sig = pd.read_csv(SIGMA, sep="\t")
    sig_geneid_to_group = dict(zip(df_sig["gene_id"].dropna(), df_sig["sigma_group"]))
    sig_geneids = set(df_sig["gene_id"].dropna())
    print(f"Sigma factors: {len(df_sig)} (mapped: {len(sig_geneids)})")

    # =========================================================
    # 5. Build master table
    # =========================================================
    print("\n--- Building master table ---")

    # Start from gene basic info (protein_coding only)
    df_master = df_genes[df_genes["gene_biotype"] == "protein_coding"].copy()
    print(f"Protein-coding genes: {len(df_master)}")

    # Add all Pfam domains
    df_master = df_master.merge(gene_all_domains, on="gene_id", how="left")
    df_master = df_master.merge(gene_tf_domains, on="gene_id", how="left")

    # pfam_TF_flag
    df_master["pfam_TF_flag"] = df_master["tf_domains"].notna().map({True: "yes", False: "no"})

    # tf_related_domains
    df_master["tf_related_domains"] = df_master["tf_domains"].fillna("")

    # ZorroAranda_regulator_type
    df_master["ZorroAranda_regulator_type"] = df_master["gene_id"].map(za_geneid_map).fillna("NA")

    # CastroMelchor_regulator_flag
    df_master["CastroMelchor_regulator_flag"] = df_master["gene_id"].isin(cm_geneids).astype(int)

    # sigma_flag
    df_master["sigma_flag"] = df_master["gene_id"].isin(sig_geneids).map({True: "yes", False: "no"})

    # sigma_group
    df_master["sigma_group"] = df_master["gene_id"].map(sig_geneid_to_group).fillna("NA")

    # TF_family (from Pfam domains)
    def get_tf_family_for_gene(row):
        if row["tf_related_domains"] == "":
            # Try from sigma or literature
            if row["sigma_flag"] == "yes":
                return "Sigma"
            elif row["ZorroAranda_regulator_type"] != "NA":
                # Check Zorro-Aranda detailed type
                za_row = df_za[df_za["gene_id"] == row["gene_id"]]
                if len(za_row) > 0:
                    reg_type = za_row.iloc[0]["regulator_type"]
                    if pd.notna(reg_type):
                        reg_lower = str(reg_type).lower()
                        if "tetr" in reg_lower:
                            return "TetR"
                        elif "gntr" in reg_lower:
                            return "GntR"
                        elif "arac" in reg_lower:
                            return "AraC"
                        elif "merr" in reg_lower:
                            return "MerR"
                        elif "lysr" in reg_lower:
                            return "LysR"
                        elif "luxr" in reg_lower:
                            return "LuxR"
                        elif "marr" in reg_lower:
                            return "MarR"
                        elif "deor" in reg_lower:
                            return "DeoR"
                        elif "asnc" in reg_lower:
                            return "AsnC"
                        elif "lsr2" in reg_lower:
                            return "Lsr2"
                        elif "sigma" in reg_lower:
                            return "Sigma"
                        elif "tcs" in reg_lower or "response_regulator" in reg_lower:
                            return "TCS_response_regulator"
                return "Unclassified_TF"
            return "NA"
        domains = row["tf_related_domains"].split(";")
        return assign_tf_family(domains)

    df_master["TF_family"] = df_master.apply(get_tf_family_for_gene, axis=1)

    # source_flags
    def build_source_flags(row):
        sources = []
        if row["pfam_TF_flag"] == "yes":
            sources.append("pfam")
        if row["ZorroAranda_regulator_type"] != "NA":
            sources.append("ZorroAranda")
        if row["CastroMelchor_regulator_flag"] == 1:
            sources.append("CastroMelchor")
        if row["sigma_flag"] == "yes":
            sources.append("SigmaReview")
        return ";".join(sources) if sources else "NA"

    df_master["source_flags"] = df_master.apply(build_source_flags, axis=1)

    # =========================================================
    # 6. Filter to TF candidates only for the master TF list
    # =========================================================
    # A gene is a TF candidate if it has at least one source flag
    tf_mask = df_master["source_flags"] != "NA"
    df_tf = df_master[tf_mask].copy()

    print(f"\nTotal TF/regulator candidates: {len(df_tf)}")

    # =========================================================
    # 7. Select output columns and save
    # =========================================================
    out_cols = [
        "SCO_ID", "gene_id", "gene_name", "contig", "start", "end", "strand",
        "pfam_TF_flag", "tf_related_domains",
        "ZorroAranda_regulator_type", "CastroMelchor_regulator_flag",
        "sigma_flag", "sigma_group", "TF_family", "source_flags",
        "product",
    ]

    df_out = df_tf[out_cols].copy()

    # Sort by contig and start position
    contig_order = {"chromosome": 0, "SCP1": 1, "SCP2": 2}
    df_out["contig_sort"] = df_out["contig"].map(contig_order).fillna(3)
    df_out = df_out.sort_values(["contig_sort", "start"]).drop(columns=["contig_sort"])

    out_path = f"{FINAL_DIR}/master_TF_list_M145.tsv"
    df_out.to_csv(out_path, sep="\t", index=False)
    print(f"\nOutput: {out_path}")

    # Also save the full table (all genes with annotations) for reference
    full_out_cols = [
        "SCO_ID", "gene_id", "gene_name", "contig", "start", "end", "strand",
        "pfam_TF_flag", "tf_related_domains", "all_domains",
        "ZorroAranda_regulator_type", "CastroMelchor_regulator_flag",
        "sigma_flag", "sigma_group", "TF_family", "source_flags", "product",
    ]
    full_path = f"{OUT_DIR}/M145_all_genes_with_TF_annotations.tsv"
    df_master[full_out_cols].to_csv(full_path, sep="\t", index=False)
    print(f"Full annotation table: {full_path}")

    # =========================================================
    # 8. Summary statistics
    # =========================================================
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)
    print(f"\nTotal protein-coding genes: {len(df_master)}")
    print(f"TF/regulator candidates: {len(df_tf)}")
    print(f"\n--- Source overlap ---")
    print(f"pfam_TF_flag == yes: {(df_tf['pfam_TF_flag']=='yes').sum()}")
    print(f"ZorroAranda (non-NA): {(df_tf['ZorroAranda_regulator_type']!='NA').sum()}")
    print(f"CastroMelchor (flag=1): {(df_tf['CastroMelchor_regulator_flag']==1).sum()}")
    print(f"sigma_flag == yes: {(df_tf['sigma_flag']=='yes').sum()}")
    print(f"\n--- TF Family distribution ---")
    print(df_tf["TF_family"].value_counts().to_string())
    print(f"\n--- Sigma group distribution (sigma_flag=yes only) ---")
    sig_only = df_tf[df_tf["sigma_flag"]=="yes"]
    print(sig_only["sigma_group"].value_counts().to_string())
    print(f"\n--- Contig distribution ---")
    print(df_tf["contig"].value_counts().to_string())

    # Source combination analysis
    print(f"\n--- Source combination analysis ---")
    print(df_tf["source_flags"].value_counts().head(20).to_string())

if __name__ == "__main__":
    main()
