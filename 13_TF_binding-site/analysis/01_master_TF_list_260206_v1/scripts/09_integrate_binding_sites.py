#!/usr/bin/env python3
"""
09_integrate_binding_sites.py
全ソースの binding site 情報を統合し、2つのテーブルを作成:

1. master_TF_binding_sites_M145.tsv
   TF × binding site の詳細テーブル（1行 = 1 binding site）

2. master_TF_list_M145_with_BS.tsv
   master_TF_list に binding site 有無とサマリを紐付けたテーブル
"""
import pandas as pd
import os

BASE_DIR = "/Users/okaban/bioinfo/rna-seq/13_TF_binding-site/analysis/01_master_TF_list_260206_v1"
OUT_DIR = f"{BASE_DIR}/intermediate"
GENE_INFO = f"{OUT_DIR}/M145_gene_basic_info.tsv"

CONTIG_MAP = {
    "NC_003888.3": "chromosome",
    "NC_003903.1": "SCP1",
    "NC_003904.1": "SCP2",
}


def main():
    # Load gene info
    df_genes = pd.read_csv(GENE_INFO, sep="\t")
    sco_to_geneid = {k: v for k, v in zip(df_genes["SCO_ID"], df_genes["gene_id"]) if pd.notna(k) and k != ""}
    geneid_to_info = {}
    for _, row in df_genes.iterrows():
        geneid_to_info[row["gene_id"]] = {
            "SCO_ID": row.get("SCO_ID", ""),
            "gene_name": row.get("gene_name", ""),
            "contig": row.get("contig", ""),
            "start": row.get("start", ""),
            "end": row.get("end", ""),
            "strand": row.get("strand", ""),
        }

    all_binding_sites = []

    # ====================================================
    # Source 1: Curated motifs (user-curated + literature)
    # ====================================================
    print("--- Source 1: Curated motifs ---")
    df_curated = pd.read_csv(f"{OUT_DIR}/curated_TF_motifs_all_sources.tsv", sep="\t")
    print(f"  {len(df_curated)} curated motif entries")

    # ====================================================
    # Source 2: FIMO genome-wide scan results
    # ====================================================
    print("--- Source 2: FIMO results ---")
    df_fimo = pd.read_csv(f"{OUT_DIR}/fimo_binding_sites_processed.tsv", sep="\t")

    # Map TF motif_id to SCO_ID
    motif_to_sco = dict(zip(df_curated["TF_name"], df_curated["SCO_ID"]))

    # Apply stricter threshold: q-value < 0.01
    df_fimo_strict = df_fimo[df_fimo["q-value"] < 0.01].copy()
    print(f"  FIMO hits (q < 0.01): {len(df_fimo_strict)}")
    print(f"  Per motif:")
    print(df_fimo_strict["motif_id"].value_counts().to_string())

    # Map to genes - find nearest gene for each FIMO hit
    # For each hit, find genes whose promoter region (-500 to +50 relative to start) overlaps
    for _, row in df_fimo_strict.iterrows():
        tf_name = row["motif_id"]
        sco_id = motif_to_sco.get(tf_name, "")
        gene_id = sco_to_geneid.get(sco_id, "")
        contig = CONTIG_MAP.get(row["sequence_name"], row["sequence_name"])

        all_binding_sites.append({
            "TF_name": tf_name,
            "TF_SCO_ID": sco_id,
            "TF_gene_id": gene_id,
            "BS_contig": contig,
            "BS_start": int(row["start"]),
            "BS_end": int(row["stop"]),
            "BS_strand": row["strand"],
            "BS_sequence": row["matched_sequence"],
            "BS_score": row["score"],
            "BS_pvalue": row["p-value"],
            "BS_qvalue": row["q-value"],
            "BS_source": "FIMO_curated_motif",
            "BS_evidence": "computational_FIMO",
        })

    # ====================================================
    # Source 3: Zorro-Aranda MEME binding sites
    # ====================================================
    print("\n--- Source 3: Zorro-Aranda MEME binding sites ---")
    df_za_meme = pd.read_csv(f"{OUT_DIR}/ZorroAranda2022_MEME_binding_sites.tsv", sep="\t")
    # Keep best hit per TF-TG pair (lowest p-value)
    df_za_best = df_za_meme.sort_values("p_value").drop_duplicates(subset=["TF_SCO", "TG_SCO"], keep="first")
    print(f"  Best MEME hits (1 per TF-TG pair): {len(df_za_best)}")

    for _, row in df_za_best.iterrows():
        tf_sco = row["TF_SCO"]
        tg_sco = row["TG_SCO"]
        tf_geneid = sco_to_geneid.get(tf_sco, "")
        tg_geneid = sco_to_geneid.get(tg_sco, "")
        tg_info = geneid_to_info.get(tg_geneid, {})

        all_binding_sites.append({
            "TF_name": tf_sco,
            "TF_SCO_ID": tf_sco,
            "TF_gene_id": tf_geneid,
            "TG_SCO_ID": tg_sco,
            "TG_gene_id": tg_geneid,
            "BS_contig": tg_info.get("contig", ""),
            "BS_start": "",  # relative positions in MEME output
            "BS_end": "",
            "BS_strand": row.get("strand", ""),
            "BS_sequence": row.get("site_sequence", ""),
            "BS_score": "",
            "BS_pvalue": row.get("p_value", ""),
            "BS_qvalue": row.get("q_value", ""),
            "BS_relative_left": row.get("left_pos", ""),
            "BS_relative_right": row.get("right_pos", ""),
            "BS_source": "ZorroAranda2022_MEME",
            "BS_evidence": "computational_MEME",
        })

    # ====================================================
    # Source 4: RegPrecise binding sites
    # ====================================================
    print("\n--- Source 4: RegPrecise binding sites ---")
    rp_path = f"{OUT_DIR}/RegPrecise_binding_sites.tsv"
    if os.path.exists(rp_path):
        df_rp = pd.read_csv(rp_path, sep="\t")
        print(f"  RegPrecise records: {len(df_rp)}")

        for _, row in df_rp.iterrows():
            tf_name = row.get("TF_name", "")
            tf_locus = row.get("TF_locus_tag", "")
            tf_geneid = sco_to_geneid.get(tf_locus, "")
            tg_locus = row.get("target_locus_tag", "")
            tg_geneid = sco_to_geneid.get(tg_locus, "")

            all_binding_sites.append({
                "TF_name": tf_name,
                "TF_SCO_ID": tf_locus,
                "TF_gene_id": tf_geneid,
                "TG_SCO_ID": tg_locus,
                "TG_gene_id": tg_geneid,
                "BS_contig": "",
                "BS_start": row.get("binding_site_position", ""),
                "BS_end": "",
                "BS_strand": "",
                "BS_sequence": row.get("binding_site_sequence", ""),
                "BS_score": row.get("score", ""),
                "BS_pvalue": "",
                "BS_qvalue": "",
                "BS_source": "RegPrecise",
                "BS_evidence": "computational_comparative",
            })
    else:
        print(f"  RegPrecise file not found")

    # ====================================================
    # Source 5: Zorro-Aranda curated interactions (no sequence, but evidence)
    # ====================================================
    print("\n--- Source 5: Zorro-Aranda curated interactions ---")
    df_za_curated = pd.read_csv(f"{OUT_DIR}/ZorroAranda2022_curated_interactions.tsv", sep="\t")
    # Keep only strong evidence interactions
    df_za_strong = df_za_curated[df_za_curated["evidence_classification"] == "Strong"]
    print(f"  Strong evidence curated interactions: {len(df_za_strong)}")

    for _, row in df_za_strong.iterrows():
        tf_sco = row.get("TF_SCO", "")
        tg_sco = row.get("TG_SCO", "")
        tf_geneid = sco_to_geneid.get(tf_sco, "")
        tg_geneid = sco_to_geneid.get(tg_sco, "")

        all_binding_sites.append({
            "TF_name": row.get("TF_name", tf_sco),
            "TF_SCO_ID": tf_sco,
            "TF_gene_id": tf_geneid,
            "TG_SCO_ID": tg_sco,
            "TG_gene_id": tg_geneid,
            "BS_contig": "",
            "BS_start": "",
            "BS_end": "",
            "BS_strand": "",
            "BS_sequence": "",
            "BS_score": "",
            "BS_pvalue": "",
            "BS_qvalue": "",
            "BS_source": "ZorroAranda2022_Curated_Strong",
            "BS_evidence": f"experimental_{row.get('experiment','')}",
            "experiment_type": row.get("experiment", ""),
            "regulatory_function": row.get("regulatory_function", ""),
        })

    # ====================================================
    # Build binding sites master table
    # ====================================================
    print("\n" + "=" * 60)
    print("Building master binding site table")
    print("=" * 60)

    df_bs = pd.DataFrame(all_binding_sites)
    # Standardize columns
    out_cols = [
        "TF_name", "TF_SCO_ID", "TF_gene_id",
        "TG_SCO_ID", "TG_gene_id",
        "BS_contig", "BS_start", "BS_end", "BS_strand",
        "BS_sequence", "BS_score", "BS_pvalue", "BS_qvalue",
        "BS_source", "BS_evidence",
    ]
    for col in out_cols:
        if col not in df_bs.columns:
            df_bs[col] = ""

    df_bs_out = df_bs[out_cols].copy()

    bs_path = f"{BASE_DIR}/master_TF_binding_sites_M145.tsv"
    df_bs_out.to_csv(bs_path, sep="\t", index=False)
    print(f"\nOutput: {bs_path}")
    print(f"Total binding site records: {len(df_bs_out)}")
    print(f"\nSource distribution:")
    print(df_bs_out["BS_source"].value_counts().to_string())

    # ====================================================
    # Create TF-level summary and merge with master_TF_list
    # ====================================================
    print("\n" + "=" * 60)
    print("Creating TF-level summary")
    print("=" * 60)

    # Load master TF list
    df_tf = pd.read_csv(f"{BASE_DIR}/master_TF_list_M145.tsv", sep="\t")

    # For each TF gene_id, summarize binding site info
    # Group by TF_gene_id
    def summarize_bs_for_tf(gene_id):
        mask = df_bs_out["TF_gene_id"] == gene_id
        sub = df_bs_out[mask]
        if len(sub) == 0:
            return pd.Series({
                "has_binding_site": "no",
                "BS_count": 0,
                "BS_sources": "",
                "BS_evidence_types": "",
                "motif_consensus": "",
            })

        sources = ";".join(sorted(sub["BS_source"].dropna().unique()))
        evidence = ";".join(sorted(sub["BS_evidence"].dropna().unique()))
        count = len(sub)

        # Get motif consensus from curated table if available
        curated_match = df_curated[df_curated["SCO_ID"].isin(
            sub["TF_SCO_ID"].dropna().unique()
        )]
        consensus = ""
        if len(curated_match) > 0:
            consensus = curated_match.iloc[0]["motif_consensus"]

        return pd.Series({
            "has_binding_site": "yes",
            "BS_count": count,
            "BS_sources": sources,
            "BS_evidence_types": evidence,
            "motif_consensus": consensus,
        })

    # Also check by SCO_ID for TFs identified by SCO name
    def summarize_bs_for_tf_sco(sco_id):
        mask = df_bs_out["TF_SCO_ID"] == sco_id
        sub = df_bs_out[mask]
        if len(sub) == 0:
            return None
        sources = ";".join(sorted(sub["BS_source"].dropna().unique()))
        evidence = ";".join(sorted(sub["BS_evidence"].dropna().unique()))
        count = len(sub)
        return {"BS_count": count, "BS_sources": sources, "BS_evidence_types": evidence}

    bs_summary = df_tf["gene_id"].apply(summarize_bs_for_tf)
    df_tf_bs = pd.concat([df_tf, bs_summary], axis=1)

    # For TFs that weren't found by gene_id, try SCO_ID
    no_bs_mask = df_tf_bs["has_binding_site"] == "no"
    for idx in df_tf_bs[no_bs_mask].index:
        sco = df_tf_bs.loc[idx, "SCO_ID"]
        if pd.notna(sco) and sco != "":
            result = summarize_bs_for_tf_sco(sco)
            if result:
                df_tf_bs.loc[idx, "has_binding_site"] = "yes"
                df_tf_bs.loc[idx, "BS_count"] = result["BS_count"]
                df_tf_bs.loc[idx, "BS_sources"] = result["BS_sources"]
                df_tf_bs.loc[idx, "BS_evidence_types"] = result["BS_evidence_types"]
                # Check curated
                curated_match = df_curated[df_curated["SCO_ID"] == sco]
                if len(curated_match) > 0:
                    df_tf_bs.loc[idx, "motif_consensus"] = curated_match.iloc[0]["motif_consensus"]

    # Save
    tf_bs_path = f"{BASE_DIR}/master_TF_list_M145_with_BS.tsv"
    df_tf_bs.to_csv(tf_bs_path, sep="\t", index=False)
    print(f"\nOutput: {tf_bs_path}")

    # Summary
    has_bs = (df_tf_bs["has_binding_site"] == "yes").sum()
    no_bs = (df_tf_bs["has_binding_site"] == "no").sum()
    has_motif = (df_tf_bs["motif_consensus"].fillna("") != "").sum()

    print(f"\n{'='*60}")
    print("FINAL SUMMARY")
    print(f"{'='*60}")
    print(f"Total TF candidates in master list: {len(df_tf_bs)}")
    print(f"TFs with any binding site info: {has_bs} ({has_bs/len(df_tf_bs)*100:.1f}%)")
    print(f"TFs without binding site info: {no_bs}")
    print(f"TFs with defined motif consensus: {has_motif}")

    print(f"\n--- TFs with BS by confidence level ---")
    for conf in ["high", "medium", "low"]:
        sub = df_tf_bs[df_tf_bs["confidence"] == conf]
        sub_bs = (sub["has_binding_site"] == "yes").sum()
        print(f"  {conf}: {sub_bs}/{len(sub)} ({sub_bs/len(sub)*100:.1f}%)")

    print(f"\n--- TFs with BS by TF family (top 15) ---")
    bs_by_family = df_tf_bs[df_tf_bs["has_binding_site"] == "yes"].groupby("TF_family").size()
    total_by_family = df_tf_bs.groupby("TF_family").size()
    family_summary = pd.DataFrame({"with_BS": bs_by_family, "total": total_by_family}).fillna(0)
    family_summary["pct"] = (family_summary["with_BS"] / family_summary["total"] * 100).round(1)
    family_summary = family_summary.sort_values("with_BS", ascending=False)
    print(family_summary.head(15).to_string())

if __name__ == "__main__":
    main()
