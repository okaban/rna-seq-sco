#!/usr/bin/env python3
"""
02_parse_castro_melchor.py
Castro-Melchor et al. 2010 の692 regulator cistrons を SCO_ID にマッピング。
AddFile3.csv (gene->operon mapping) と AddFile4 (692 regulator cistron columns) を使用。
"""
import pandas as pd
import re

LIT_DIR = "/Users/okaban/bioinfo/rna-seq/13_TF_binding-site/analysis/01_master_TF_list_260206_v1/literature"
OUT_DIR = "/Users/okaban/bioinfo/rna-seq/13_TF_binding-site/analysis/01_master_TF_list_260206_v1/intermediate"

def main():
    # 1. Load gene->operon mapping
    gene_operon = pd.read_csv(f"{LIT_DIR}/Castro-Melchor_2010_AddFile3.csv")
    print(f"Gene-operon mapping: {len(gene_operon)} entries")
    print(gene_operon.head())

    # 2. Extract 692 regulator cistron IDs from AddFile4 column headers
    # Read only the first two rows to get the headers
    with open(f"{LIT_DIR}/Castro-Melchor_2010_AddFile4_network.csv") as f:
        header_line1 = f.readline().strip()  # numeric 1..692
        header_line2 = f.readline().strip()  # cistron names

    # Parse cistron names from line 2
    parts = header_line2.split(",")
    regulator_cistrons = []
    for p in parts[1:]:  # skip first empty column
        cistron_name = p.strip().strip("'\"")
        if cistron_name:
            # Extract the numeric part
            m = re.match(r'cistron(\d+)', cistron_name)
            if m:
                regulator_cistrons.append(int(m.group(1)))

    print(f"\nRegulator cistrons: {len(regulator_cistrons)}")

    # 3. Map cistron -> operonID
    # In Castro-Melchor, cistron = operon unit. The AddFile3 maps gene->operonID.
    # However, the cistron numbering may differ from operonID.
    # Let's check: the paper uses "cistron" as the unit. Let's see if operonID matches cistron number.

    # Actually, looking at the data: genes map to operonIDs (1, 2, 3, ...).
    # The regulator cistrons are listed as cistron26, cistron31, etc.
    # These cistron numbers should correspond to the operonID in AddFile3.

    # Create operonID -> list of genes mapping
    operon_to_genes = gene_operon.groupby('operonID')['gene'].apply(list).to_dict()

    # 4. For each regulator cistron, find the genes
    regulator_genes = []
    missing_cistrons = []
    for cid in regulator_cistrons:
        if cid in operon_to_genes:
            for gene in operon_to_genes[cid]:
                # Standardize SCO_ID format (SCO0001 -> SCO0001)
                sco_id = gene.strip()
                regulator_genes.append({
                    "SCO_ID": sco_id,
                    "cistron_id": f"cistron{cid}",
                    "regulator_flag": 1,
                    "source": "CastroMelchor_2010"
                })
        else:
            missing_cistrons.append(cid)

    if missing_cistrons:
        print(f"WARNING: {len(missing_cistrons)} cistrons not found in operon mapping")
        print(f"  First 10: {missing_cistrons[:10]}")

    df_reg = pd.DataFrame(regulator_genes)
    print(f"\nRegulator genes (all genes in regulator cistrons): {len(df_reg)}")
    print(f"Unique SCO_IDs: {df_reg['SCO_ID'].nunique()}")

    # Save
    out_path = f"{OUT_DIR}/CastroMelchor2010_regulators.tsv"
    df_reg.to_csv(out_path, sep="\t", index=False)
    print(f"Output: {out_path}")

    # Also create a deduplicated version at the SCO_ID level
    df_dedup = df_reg.drop_duplicates(subset='SCO_ID')[['SCO_ID', 'regulator_flag', 'source']]
    dedup_path = f"{OUT_DIR}/CastroMelchor2010_regulators_dedup.tsv"
    df_dedup.to_csv(dedup_path, sep="\t", index=False)
    print(f"Deduplicated output: {dedup_path} ({len(df_dedup)} unique genes)")

if __name__ == "__main__":
    main()
