#!/usr/bin/env python3
"""
08_create_motif_table_and_fimo.py
1. ユーザーキュレーション + サブエージェント収集のモチーフ情報を統合
2. FIMO用MEMEフォーマットのモチーフファイルを生成
3. FIMOでゲノムワイドスキャンを実行
"""
import pandas as pd
import subprocess
import os
import re

BASE_DIR = "/Users/okaban/bioinfo/rna-seq/13_TF_binding-site/analysis/01_master_TF_list_260206_v1"
OUT_DIR = f"{BASE_DIR}/intermediate"
GENOME = "/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/GCF_000203835.1_ASM20383v1_genomic.fna"

# =====================================================================
# Expert-curated TF binding motifs from user + literature
# IUPAC ambiguity codes will be expanded for FIMO
# =====================================================================
CURATED_MOTIFS = [
    # --- 1. Nutrient response / Global regulators ---
    {
        "TF_name": "GlnR", "SCO_ID": "SCO4159",
        "motif_consensus": "GTNACNNNNNNGANAC",
        "motif_length": 16, "motif_type": "single_box",
        "evidence_type": "experimental",
        "technique": "DNase_I_footprint;EMSA",
        "reference": "Tiffert_2008;Wang_2013",
        "PMID": "18179599;23543735",
        "notes": "22bp a-site+b-site; here single box unit gTnAc-n6-GaAAc"
    },
    {
        "TF_name": "PhoP", "SCO_ID": "SCO7637",
        "motif_consensus": "GTTCACCNNNNGTTCACC",
        "motif_length": 18, "motif_type": "direct_repeat",
        "evidence_type": "experimental",
        "technique": "ChIP-chip;DNase_I_footprint;EMSA",
        "reference": "Allenby_2012;Sola-Landa_2005",
        "PMID": "22923519;16042609",
        "notes": "PHO box: GTTCACC-N4-GTTCACC direct repeat"
    },
    {
        "TF_name": "DasR", "SCO_ID": "SCO5231",
        "motif_consensus": "ANTGGTCTAGACCANT",
        "motif_length": 16, "motif_type": "palindrome",
        "evidence_type": "experimental",
        "technique": "ChIP-seq;EMSA",
        "reference": "Colson_2007;Swiatek-Polatynska_2015",
        "PMID": "17183212;25826289",
        "notes": "dre-site: A(G/C)TGGTCTAGACCA(G/C)T canonical palindrome"
    },
    {
        "TF_name": "AfsQ1", "SCO_ID": "SCO4907",
        "motif_consensus": "GTNACNNNNNNGTNAC",
        "motif_length": 16, "motif_type": "direct_repeat",
        "evidence_type": "experimental",
        "technique": "DNase_I_footprint",
        "reference": "Wang_2013_MolMicrobiol",
        "PMID": "23116219",
        "notes": "GTnAC-n6-GTnAC direct repeat"
    },
    {
        "TF_name": "DraR", "SCO_ID": "SCO3063",
        "motif_consensus": "AMAAWYMAKCA",
        "motif_length": 11, "motif_type": "consensus",
        "evidence_type": "experimental",
        "technique": "DNase_I_footprint;EMSA",
        "reference": "Yu_2012",
        "PMID": "22676800",
        "notes": "M=A/C,W=A/T,Y=C/T,K=G/T"
    },
    {
        "TF_name": "Crp", "SCO_ID": "SCO3571",
        "motif_consensus": "TGTGANNNNNNTCACA",
        "motif_length": 16, "motif_type": "palindrome",
        "evidence_type": "computational",
        "technique": "ChIP-chip;homology",
        "reference": "Chan_2013",
        "PMID": "23292766",
        "notes": "CRP-like motif; exact S. coelicolor consensus is approximate"
    },
    # --- 2. GBL / Quorum sensing ---
    {
        "TF_name": "ScbR", "SCO_ID": "SCO6265",
        "motif_consensus": "TTTGGNNNNNNNNNNNNNNNNNNCCAAA",
        "motif_length": 28, "motif_type": "inverted_repeat",
        "evidence_type": "experimental",
        "technique": "DNase_I_footprint;MEME",
        "reference": "Li_2015;Wang_2011",
        "PMID": "26452094;21771106",
        "notes": "Two inverted TTTGG half-sites ~30bp; core TTTSTT"
    },
    # --- 3. Antibiotic pathway / CSR ---
    {
        "TF_name": "AbrC3", "SCO_ID": "SCO4596",
        "motif_consensus": "GAASGSGRMS",
        "motif_length": 10, "motif_type": "consensus",
        "evidence_type": "experimental",
        "technique": "EMSA;ChIP-qPCR",
        "reference": "Rico_2014",
        "PMID": "24488312",
        "notes": "S=G/C,R=A/G,M=A/C"
    },
    # --- 4. Development / Morphogenesis ---
    {
        "TF_name": "BldD", "SCO_ID": "SCO1489",
        "motif_consensus": "AGTGANNNNNNTCACT",
        "motif_length": 16, "motif_type": "palindrome",
        "evidence_type": "experimental",
        "technique": "ChIP-chip;DNase_I_footprint;EMSA",
        "reference": "Elliot_1999;Schumacher_2017",
        "PMID": "10542188;28541410",
        "notes": "BldD box: AGTgA-(N)m-TCACt palindromic"
    },
    {
        "TF_name": "AdpA", "SCO_ID": "SCO2792",
        "motif_consensus": "TGGCSNGWWY",
        "motif_length": 10, "motif_type": "consensus",
        "evidence_type": "experimental",
        "technique": "EMSA;DNase_I_footprint;lacZ_reporter",
        "reference": "Ohnishi_2005;Wolanski_2011",
        "PMID": "15857081;21784946",
        "notes": "S=G/C,W=A/T,Y=C/T,N=any"
    },
    # --- 5. Sigma factor promoter consensus ---
    {
        "TF_name": "HrdB", "SCO_ID": "SCO5820",
        "motif_consensus": "TTGACANNNNNNNNNNNNNNNNNTAGAAT",
        "motif_length": 28, "motif_type": "promoter_-35_-10",
        "evidence_type": "experimental",
        "technique": "promoter_compilation;in_vitro_transcription",
        "reference": "Strohl_1992",
        "PMID": "1549509",
        "notes": "sigma70-like: TTGACa-N16-18-TAgaaT"
    },
    {
        "TF_name": "SigR", "SCO_ID": "SCO5216",
        "motif_consensus": "GGAACNNNNNNNNNNNNNNNNCGTT",
        "motif_length": 24, "motif_type": "promoter_-35_-10",
        "evidence_type": "experimental",
        "technique": "in_vitro_transcription;ChIP-chip",
        "reference": "Paget_2001;Kim_2012",
        "PMID": "11737643;22651816",
        "notes": "ECF sigma: GGAAC-N16/17-CGTT"
    },
    {
        "TF_name": "SigE", "SCO_ID": "SCO5147",
        "motif_consensus": "GGAACNNNNNNNNNNNNNNNNNGTT",
        "motif_length": 24, "motif_type": "promoter_-35_-10",
        "evidence_type": "experimental",
        "technique": "ChIP-seq;in_vitro_transcription",
        "reference": "Tran_2019",
        "PMID": "30900386",
        "notes": "ECF sigma: GGAAC-N16/17-GTT"
    },
    {
        "TF_name": "SigB", "SCO_ID": "SCO0600",
        "motif_consensus": "GNNTNNNNNNNNNNNNNNNNNGGGTAC",
        "motif_length": 26, "motif_type": "promoter_-35_-10",
        "evidence_type": "experimental",
        "technique": "consensus-directed_search",
        "reference": "Lee_2004",
        "PMID": "15357310",
        "notes": "SigB-like promoter: GNNTN14-16GGGTAC/T"
    },
    # --- Additional from RegPrecise ---
    {
        "TF_name": "ArgR", "SCO_ID": "SCO1576",
        "motif_consensus": "TGAATAANNNNNNTTATTCA",
        "motif_length": 20, "motif_type": "palindrome",
        "evidence_type": "experimental",
        "technique": "EMSA;DNase_I_footprint",
        "reference": "Perez-Redondo_2012",
        "PMID": "22403700",
        "notes": "ARG box: TGAATAA-n6-TTATTCA"
    },
    {
        "TF_name": "NdgR", "SCO_ID": "SCO5552",
        "motif_consensus": "BTYCANYWSYBGGAC",
        "motif_length": 15, "motif_type": "palindrome",
        "evidence_type": "experimental",
        "technique": "ChIP-seq;EMSA",
        "reference": "Kim_2015;Yang_2009",
        "PMID": "25766138;19083232",
        "notes": "IclR-family; [GT]T[CT]CAC[CA][CTA][TC][GC][TC]GGAC"
    },
]


# IUPAC to regex mapping
IUPAC = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
    'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
    'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]',
}

# IUPAC to probability for MEME minimal format
IUPAC_PROB = {
    'A': [1, 0, 0, 0], 'C': [0, 1, 0, 0], 'G': [0, 0, 1, 0], 'T': [0, 0, 0, 1],
    'R': [0.5, 0, 0.5, 0], 'Y': [0, 0.5, 0, 0.5], 'S': [0, 0.5, 0.5, 0],
    'W': [0.5, 0, 0, 0.5], 'K': [0, 0, 0.5, 0.5], 'M': [0.5, 0.5, 0, 0],
    'B': [0, 0.333, 0.333, 0.333], 'D': [0.333, 0, 0.333, 0.333],
    'H': [0.333, 0.333, 0, 0.333], 'V': [0.333, 0.333, 0.333, 0],
    'N': [0.25, 0.25, 0.25, 0.25],
}


def consensus_to_meme(motifs, output_path):
    """Convert IUPAC consensus motifs to MEME minimal format."""
    with open(output_path, 'w') as f:
        f.write("MEME version 5\n\n")
        f.write("ALPHABET= ACGT\n\n")
        f.write("strands: + -\n\n")
        # S. coelicolor has ~72% GC content
        f.write("Background letter frequencies\n")
        f.write("A 0.14 C 0.36 G 0.36 T 0.14\n\n")

        for m in motifs:
            name = m["TF_name"]
            consensus = m["motif_consensus"].upper()
            w = len(consensus)
            f.write(f"MOTIF {name} {m['SCO_ID']}\n")
            f.write(f"letter-probability matrix: alength= 4 w= {w}\n")
            for base in consensus:
                probs = IUPAC_PROB.get(base, [0.25, 0.25, 0.25, 0.25])
                f.write(f" {probs[0]:.4f} {probs[1]:.4f} {probs[2]:.4f} {probs[3]:.4f}\n")
            f.write("\n")
    print(f"MEME motif file: {output_path} ({len(motifs)} motifs)")


def main():
    # 1. Save curated motif table
    df_motifs = pd.DataFrame(CURATED_MOTIFS)
    motif_tsv = f"{OUT_DIR}/curated_TF_motifs_all_sources.tsv"
    df_motifs.to_csv(motif_tsv, sep="\t", index=False)
    print(f"Curated motifs: {len(df_motifs)} entries")
    print(f"Output: {motif_tsv}")

    # 2. Generate MEME format motif file
    meme_file = f"{OUT_DIR}/curated_TF_motifs.meme"
    consensus_to_meme(CURATED_MOTIFS, meme_file)

    # 3. Run FIMO
    fimo_outdir = f"{OUT_DIR}/fimo_results"
    os.makedirs(fimo_outdir, exist_ok=True)

    cmd = [
        "fimo",
        "--thresh", "1e-4",
        "--oc", fimo_outdir,
        "--verbosity", "1",
        meme_file,
        GENOME,
    ]
    print(f"\nRunning FIMO: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"FIMO stderr: {result.stderr}")
        print(f"FIMO returned {result.returncode}")
    else:
        print("FIMO completed successfully")

    # 4. Parse FIMO results
    fimo_tsv = f"{fimo_outdir}/fimo.tsv"
    if os.path.exists(fimo_tsv):
        df_fimo = pd.read_csv(fimo_tsv, sep="\t", comment="#")
        print(f"\nFIMO results: {len(df_fimo)} hits")
        print(f"Hits per motif:")
        print(df_fimo["motif_id"].value_counts().to_string())

        # Save processed version
        fimo_processed = f"{OUT_DIR}/fimo_binding_sites_processed.tsv"
        df_fimo.to_csv(fimo_processed, sep="\t", index=False)
        print(f"\nProcessed FIMO output: {fimo_processed}")
    else:
        print(f"WARNING: FIMO output not found at {fimo_tsv}")

if __name__ == "__main__":
    main()
