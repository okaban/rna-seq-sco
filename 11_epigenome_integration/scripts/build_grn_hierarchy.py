#!/usr/bin/env python3
"""
Build Hierarchical Gene Regulatory Network (GRN) for M145
==========================================================
Construct a comprehensive GRN based on literature and existing data.

Hierarchy:
- Tier 1: Global regulators (BldA, BldD, AdpA, AfsR, etc.)
- Tier 2: Pleiotropic regulators (AbsA, AbsB, AfsK, etc.)
- Tier 3: Cluster-situated regulators (CSRs) - SARP family, etc.
- Tier 4: BGC structural genes

References:
- Zorro-Aranda et al. (2022) Sci Rep - Global GRN reconstruction
- Lee et al. (2024) Adv Sci - ML analysis of transcriptomes
- Literature compilation (1990-2024)

Author: Claude
Date: 2026-02-03
"""

import pandas as pd
import numpy as np
from pathlib import Path
import json

# Paths
OUTPUT_DIR = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/12_grn_tf_methylation")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

SARP_FILE = Path("/Users/okaban/bioinfo/rna-seq/09_SARP_integration/analysis/09_SARP_integration_260128_v2/tables/SARP_list_M145_v2.tsv")
REGULATOR_FILE = Path("/Users/okaban/bioinfo/rna-seq/09_SARP_integration/analysis/09_SARP_integration_260128_v2/tables/regulator_master_with_SARP_v2.tsv")

# =============================================================================
# Literature-based GRN knowledge base
# =============================================================================

# Tier 1: Global Regulators
# These affect morphological and physiological differentiation globally
GLOBAL_REGULATORS = {
    # bld genes (developmental regulators)
    "SC_RS26640": {"name": "bldA", "product": "tRNA-Leu", "function": "UUA codon translation, master regulator", "tier": 1},
    "SC_RS25655": {"name": "bldB", "product": "BldB", "function": "Morphological differentiation", "tier": 1},
    "SC_RS36355": {"name": "bldC", "product": "BldC", "function": "MerR family, aerial mycelium", "tier": 1},
    "SC_RS25165": {"name": "bldD", "product": "BldD", "function": "Master regulator of development, c-di-GMP binding", "tier": 1},
    "SC_RS25070": {"name": "bldG", "product": "Anti-sigma factor antagonist", "function": "σF regulation", "tier": 1},
    "SC_RS14000": {"name": "bldH/adpA", "product": "AdpA", "function": "AraC family, pleiotropic regulator", "tier": 1},
    "SC_RS18130": {"name": "bldM", "product": "Response regulator", "function": "Orphan response regulator", "tier": 1},
    "SC_RS26120": {"name": "bldN", "product": "σBldN", "function": "ECF sigma factor, aerial development", "tier": 1},

    # AfsR system
    "SC_RS24295": {"name": "afsR", "product": "AfsR", "function": "SARP-like global regulator, phosphorylated by AfsK", "tier": 1},
    "SC_RS24290": {"name": "afsK", "product": "Ser/Thr kinase", "function": "Phosphorylates AfsR", "tier": 1},
    "SC_RS22980": {"name": "afsS", "product": "AfsS", "function": "σ-like factor, secondary metabolism", "tier": 1},

    # Other global regulators
    "SC_RS10725": {"name": "crp", "product": "CRP/FNR family", "function": "cAMP receptor protein", "tier": 1},
    "SC_RS19305": {"name": "dasR", "product": "DasR", "function": "GntR family, N-acetylglucosamine signaling", "tier": 1},
}

# Tier 2: Pleiotropic Regulators
# Affect multiple BGCs
PLEIOTROPIC_REGULATORS = {
    "SC_RS02970": {"name": "absA2", "product": "Response regulator", "function": "Two-component, represses multiple BGCs", "tier": 2},
    "SC_RS02965": {"name": "absA1", "product": "Sensor kinase", "function": "Two-component with AbsA2", "tier": 2},
    "SC_RS27825": {"name": "absB", "product": "RNase III", "function": "Post-transcriptional regulation of regulators", "tier": 2},
    "SC_RS33490": {"name": "afsQ1", "product": "Sensor kinase", "function": "Two-component, secondary metabolism", "tier": 2},
    "SC_RS33495": {"name": "afsQ2", "product": "Response regulator", "function": "Two-component with AfsQ1", "tier": 2},
    "SC_RS35450": {"name": "nsdA", "product": "NsdA", "function": "Negative regulator of antibiotic synthesis", "tier": 2},
    "SC_RS37265": {"name": "nsdB", "product": "NsdB", "function": "Negative regulator, MarR family", "tier": 2},
    "SC_RS04915": {"name": "wblA", "product": "WhiB-like", "function": "Wbl family, negative regulator", "tier": 2},
}

# Tier 3: Cluster-Situated Regulators (CSRs)
# Pathway-specific regulators within BGCs
CLUSTER_SITUATED_REGULATORS = {
    # Actinorhodin (Act) cluster
    "SC_RS27570": {"name": "actII-ORF4", "product": "SARP regulator", "function": "Act pathway activator", "tier": 3, "bgc": "act"},

    # Undecylprodigiosin (Red) cluster
    "SC_RS27225": {"name": "redD", "product": "SARP regulator", "function": "Red pathway activator", "tier": 3, "bgc": "red"},
    "SC_RS27300": {"name": "redZ", "product": "Response regulator", "function": "Activates redD transcription", "tier": 3, "bgc": "red"},

    # CDA cluster
    "SC_RS17785": {"name": "cdaR", "product": "SARP regulator", "function": "CDA pathway activator", "tier": 3, "bgc": "cda"},

    # Cpk/coelimycin cluster
    "SC_RS31605": {"name": "cpkO/kasO", "product": "SARP regulator", "function": "Cpk pathway activator", "tier": 3, "bgc": "cpk"},

    # Other SARP family
    "SC_RS18200": {"name": "papR2", "product": "AfsR/SARP", "function": "Activates multiple BGCs", "tier": 3, "bgc": "multiple"},
}

# Sigma factors
SIGMA_FACTORS = {
    "SC_RS25850": {"name": "hrdB", "product": "σHrdB", "function": "Principal sigma factor", "tier": 1},
    "SC_RS23200": {"name": "hrdA", "product": "σHrdA", "function": "Sigma-70 family", "tier": 2},
    "SC_RS25845": {"name": "hrdC", "product": "σHrdC", "function": "Sigma-70 family", "tier": 2},
    "SC_RS25840": {"name": "hrdD", "product": "σHrdD", "function": "Transcribes actII-ORF4, redD", "tier": 2},
    "SC_RS24925": {"name": "sigE", "product": "σE", "function": "ECF sigma, cell envelope stress", "tier": 2},
    "SC_RS25290": {"name": "sigH", "product": "σH", "function": "Heat shock, osmotic stress", "tier": 2},
    "SC_RS26330": {"name": "sigR", "product": "σR", "function": "Thiol-oxidative stress", "tier": 2},
    "SC_RS04925": {"name": "sigF", "product": "σF", "function": "Sporulation sigma factor", "tier": 2},
    "SC_RS26785": {"name": "sigB", "product": "σB", "function": "Osmotic stress", "tier": 2},
    "SC_RS04965": {"name": "sigU", "product": "σU", "function": "ECF sigma factor", "tier": 2},
}

# Known regulatory interactions (TF -> target)
REGULATORY_INTERACTIONS = [
    # BldD regulon
    {"regulator": "bldD", "target": "adpA", "effect": "repression", "evidence": "ChIP"},
    {"regulator": "bldD", "target": "bldN", "effect": "repression", "evidence": "ChIP"},
    {"regulator": "bldD", "target": "bldM", "effect": "repression", "evidence": "ChIP"},
    {"regulator": "bldD", "target": "whiG", "effect": "repression", "evidence": "ChIP"},

    # AdpA regulon
    {"regulator": "adpA", "target": "bldA", "effect": "activation", "evidence": "ChIP"},
    {"regulator": "adpA", "target": "ramR", "effect": "activation", "evidence": "ChIP"},
    {"regulator": "adpA", "target": "sti1", "effect": "activation", "evidence": "ChIP"},

    # BldA translational control
    {"regulator": "bldA", "target": "adpA", "effect": "translation", "evidence": "TTA codon"},
    {"regulator": "bldA", "target": "actII-ORF4", "effect": "translation", "evidence": "TTA codon"},

    # AbsA1/A2 regulon
    {"regulator": "absA2", "target": "actII-ORF4", "effect": "repression", "evidence": "genetic"},
    {"regulator": "absA2", "target": "redD", "effect": "repression", "evidence": "genetic"},
    {"regulator": "absA2", "target": "cdaR", "effect": "repression", "evidence": "genetic"},

    # CSR -> BGC genes
    {"regulator": "actII-ORF4", "target": "actI-ORF1", "effect": "activation", "evidence": "DNase footprint"},
    {"regulator": "actII-ORF4", "target": "actIII", "effect": "activation", "evidence": "genetic"},
    {"regulator": "redD", "target": "redL", "effect": "activation", "evidence": "genetic"},
    {"regulator": "redZ", "target": "redD", "effect": "activation", "evidence": "genetic"},

    # AfsR cascade
    {"regulator": "afsR", "target": "actII-ORF4", "effect": "activation", "evidence": "genetic"},
    {"regulator": "afsR", "target": "redD", "effect": "activation", "evidence": "genetic"},
    {"regulator": "afsK", "target": "afsR", "effect": "phosphorylation", "evidence": "biochemical"},
]


def build_tf_master_table():
    """Combine all TF sources into master table."""
    print("=" * 60)
    print("Building TF Master Table")
    print("=" * 60)

    # Combine literature knowledge
    all_tfs = {}

    for tf_dict, category in [(GLOBAL_REGULATORS, "global"),
                               (PLEIOTROPIC_REGULATORS, "pleiotropic"),
                               (CLUSTER_SITUATED_REGULATORS, "csr"),
                               (SIGMA_FACTORS, "sigma")]:
        for locus, info in tf_dict.items():
            info["category"] = category
            info["locus_tag"] = locus
            all_tfs[locus] = info

    # Convert to DataFrame
    lit_tf_df = pd.DataFrame.from_dict(all_tfs, orient='index').reset_index(drop=True)

    print(f"\nLiterature-based TFs: {len(lit_tf_df)}")
    print(f"  Global regulators: {len(GLOBAL_REGULATORS)}")
    print(f"  Pleiotropic: {len(PLEIOTROPIC_REGULATORS)}")
    print(f"  CSRs: {len(CLUSTER_SITUATED_REGULATORS)}")
    print(f"  Sigma factors: {len(SIGMA_FACTORS)}")

    # Load existing SARP data
    if SARP_FILE.exists():
        sarp_df = pd.read_csv(SARP_FILE, sep='\t')
        print(f"\nExisting SARP data: {len(sarp_df)} genes")
    else:
        sarp_df = pd.DataFrame()

    # Load existing regulator data
    if REGULATOR_FILE.exists():
        reg_df = pd.read_csv(REGULATOR_FILE, sep='\t')
        print(f"Existing regulator data: {len(reg_df)} genes")
    else:
        reg_df = pd.DataFrame()

    return lit_tf_df, sarp_df, reg_df


def create_hierarchical_network(lit_tf_df, reg_df):
    """Create hierarchical regulatory network structure."""
    print("\n" + "=" * 60)
    print("Creating Hierarchical Network")
    print("=" * 60)

    # Merge literature TFs with expression data
    if len(reg_df) > 0:
        # Find matching genes
        merged = lit_tf_df.merge(
            reg_df[['gene_id', 'log2FC_3_vs_1', 'padj_3_vs_1', 'corr_act', 'corr_red', 'corr_cda', 'corr_cpk']],
            left_on='locus_tag',
            right_on='gene_id',
            how='left'
        )
    else:
        merged = lit_tf_df

    # Add tier information
    tier_counts = merged.groupby('tier').size()
    print("\nNetwork Hierarchy:")
    for tier, count in tier_counts.items():
        tier_name = {1: "Global", 2: "Pleiotropic/Sigma", 3: "CSR"}
        print(f"  Tier {tier} ({tier_name.get(tier, 'Other')}): {count} genes")

    return merged


def create_regulatory_interaction_table():
    """Create table of known regulatory interactions."""
    print("\n" + "=" * 60)
    print("Creating Regulatory Interaction Table")
    print("=" * 60)

    interactions_df = pd.DataFrame(REGULATORY_INTERACTIONS)
    print(f"\nKnown interactions: {len(interactions_df)}")
    print(f"  Activation: {len(interactions_df[interactions_df['effect'] == 'activation'])}")
    print(f"  Repression: {len(interactions_df[interactions_df['effect'] == 'repression'])}")
    print(f"  Translation: {len(interactions_df[interactions_df['effect'] == 'translation'])}")
    print(f"  Phosphorylation: {len(interactions_df[interactions_df['effect'] == 'phosphorylation'])}")

    return interactions_df


def analyze_bgc_specific_regulators(merged_df):
    """Analyze regulators specific to Act and Red BGCs."""
    print("\n" + "=" * 60)
    print("BGC-Specific Regulator Analysis (Act & Red focus)")
    print("=" * 60)

    # Act cluster regulators
    act_regulators = [
        "actII-ORF4",  # Direct activator
        "afsR",        # Global activator
        "absA2",       # Repressor
        "bldA",        # Translation (TTA codon)
        "hrdD",        # Sigma factor for actII-ORF4
    ]

    # Red cluster regulators
    red_regulators = [
        "redD",        # Direct activator
        "redZ",        # Activates redD
        "afsR",        # Global activator
        "absA2",       # Repressor
        "bldA",        # Translation (redZ TTA codon)
        "hrdD",        # Sigma factor for redD
    ]

    print("\nActinorhodin (Act) Cluster Regulatory Cascade:")
    print("-" * 40)
    print("  bldA (tRNA) → [translates] → actII-ORF4 mRNA")
    print("  σHrdD → [transcribes] → actII-ORF4")
    print("  AfsR → [activates] → actII-ORF4")
    print("  AbsA2 → [represses] → actII-ORF4")
    print("  ActII-ORF4 → [activates] → act biosynthetic genes")

    print("\nUndecylprodigiosin (Red) Cluster Regulatory Cascade:")
    print("-" * 40)
    print("  bldA (tRNA) → [translates] → redZ mRNA")
    print("  RedZ → [activates] → redD")
    print("  σHrdD → [transcribes] → redD")
    print("  AfsR → [activates] → redD")
    print("  AbsA2 → [represses] → redD")
    print("  RedD → [activates] → red biosynthetic genes")

    # Find these genes in our data
    bgc_info = {
        "Act Cluster": {
            "CSR": "actII-ORF4 (SC_RS27570)",
            "global_activators": ["afsR", "bldA"],
            "repressors": ["absA2", "nsdA"],
            "sigma": "hrdD"
        },
        "Red Cluster": {
            "CSR": "redD (SC_RS27225), redZ (SC_RS27300)",
            "global_activators": ["afsR", "bldA"],
            "repressors": ["absA2", "nsdA"],
            "sigma": "hrdD"
        }
    }

    return bgc_info


def generate_grn_summary():
    """Generate summary of GRN for report."""

    summary = {
        "total_literature_tfs": len(GLOBAL_REGULATORS) + len(PLEIOTROPIC_REGULATORS) + len(CLUSTER_SITUATED_REGULATORS) + len(SIGMA_FACTORS),
        "global_regulators": len(GLOBAL_REGULATORS),
        "pleiotropic_regulators": len(PLEIOTROPIC_REGULATORS),
        "csrs": len(CLUSTER_SITUATED_REGULATORS),
        "sigma_factors": len(SIGMA_FACTORS),
        "known_interactions": len(REGULATORY_INTERACTIONS),
        "key_regulators_act": ["actII-ORF4", "afsR", "absA2", "bldA", "hrdD"],
        "key_regulators_red": ["redD", "redZ", "afsR", "absA2", "bldA", "hrdD"],
    }

    return summary


def main():
    """Main analysis pipeline."""
    print("\n" + "=" * 70)
    print("M145 Gene Regulatory Network (GRN) Construction")
    print("=" * 70)

    # 1. Build TF master table
    lit_tf_df, sarp_df, reg_df = build_tf_master_table()

    # 2. Create hierarchical network
    network_df = create_hierarchical_network(lit_tf_df, reg_df)

    # 3. Create interaction table
    interactions_df = create_regulatory_interaction_table()

    # 4. BGC-specific analysis
    bgc_info = analyze_bgc_specific_regulators(network_df)

    # 5. Generate summary
    summary = generate_grn_summary()

    # 6. Save outputs
    lit_tf_df.to_csv(OUTPUT_DIR / "literature_tf_master.csv", index=False)
    network_df.to_csv(OUTPUT_DIR / "hierarchical_network_tfs.csv", index=False)
    interactions_df.to_csv(OUTPUT_DIR / "regulatory_interactions.csv", index=False)

    with open(OUTPUT_DIR / "grn_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)

    print("\n" + "=" * 70)
    print("GRN Construction Complete")
    print("=" * 70)
    print(f"\nOutput directory: {OUTPUT_DIR}")
    print(f"\nFiles created:")
    print(f"  - literature_tf_master.csv")
    print(f"  - hierarchical_network_tfs.csv")
    print(f"  - regulatory_interactions.csv")
    print(f"  - grn_summary.json")

    print(f"\nGRN Summary:")
    print(f"  Total TFs from literature: {summary['total_literature_tfs']}")
    print(f"  Known regulatory interactions: {summary['known_interactions']}")

    return network_df, interactions_df, summary


if __name__ == "__main__":
    network_df, interactions_df, summary = main()
