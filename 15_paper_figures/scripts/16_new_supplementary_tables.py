#!/usr/bin/env python3
"""
Generate Supplementary Tables ST4–ST8 for the Gatekeeper Model paper.

ST1–ST3 are carried over from existing outputs.
ST4: 1,055 regulatory genes — shielded/exposed classification + features
ST5: 62 exposed TFs — full annotation
ST6: TCS pair analysis (7 pairs)
ST7: Simpson's paradox statistics
ST8: Hypothesis ledger (H1–H36)
"""

import sys, importlib
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent))
_utils = importlib.import_module('00_shared_utils')
for _attr in dir(_utils):
    if not _attr.startswith('_'):
        globals()[_attr] = getattr(_utils, _attr)

import numpy as np
import pandas as pd

TABLE_SUP_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# ST4: 1,055 regulatory genes — shielded/exposed classification
# ============================================================
print("Generating ST4...")
df4 = load_all_genes_features()
# Select and reorder columns for publication
cols_st4 = [
    'locus_tag', 'old_locus_tag', 'gene_name', 'product', 'tf_family',
    'region', 'is_exposed', 'nearest_methyl_distance', 'n_methyl_sites_2kb',
    'baseMean', 'LFC_T2vsT1', 'LFC_T3vsT1', 'n_FIMO_hits',
]
# Only keep columns that exist
cols_st4 = [c for c in cols_st4 if c in df4.columns]
st4 = df4[cols_st4].copy()
st4['classification'] = st4['is_exposed'].map({1: 'exposed', 0: 'shielded'})
st4 = st4.drop(columns=['is_exposed'])
# Round numeric columns
for c in ['nearest_methyl_distance']:
    if c in st4.columns:
        st4[c] = st4[c].round(0).astype(int)
for c in ['baseMean', 'LFC_T2vsT1', 'LFC_T3vsT1']:
    if c in st4.columns:
        st4[c] = st4[c].round(3)
st4.to_csv(TABLE_SUP_DIR / 'ST4_regulatory_genes_classification.tsv',
           sep='\t', index=False)
print(f"  ST4: {len(st4)} genes written")

# ============================================================
# ST5: 62 exposed TFs — full annotation
# ============================================================
print("Generating ST5...")
exp = load_exposed_regulators()
tc = load_temporal_classification()

# Merge temporal data
st5 = exp.copy()
if 'locus_tag' in tc.columns and 'locus_tag' in st5.columns:
    tc_cols = ['locus_tag', 'bloc', 'module', 'phase_ratio', 'temporal_class',
               'LFC_T2vsT1', 'LFC_T3vsT1', 'T1_z', 'T2_z', 'T3_z']
    tc_cols = [c for c in tc_cols if c in tc.columns]
    # Avoid column name conflicts
    merge_cols = [c for c in tc_cols if c not in st5.columns or c == 'locus_tag']
    st5 = st5.merge(tc[merge_cols], on='locus_tag', how='left')

# Select key columns for publication
desired_cols = [
    'locus_tag', 'old_locus_tag', 'gene_name', 'product', 'tf_family',
    'region', 'nearest_methyl_dist', 'nearest_methyl_distance',
    'bloc', 'module', 'phase_ratio', 'temporal_class',
    'baseMean', 'LFC_T2vsT1', 'LFC_T3vsT1',
    'T1_z', 'T2_z', 'T3_z',
]
cols_st5 = [c for c in desired_cols if c in st5.columns]
st5 = st5[cols_st5].copy()
# Round
for c in st5.select_dtypes(include=[np.floating]).columns:
    st5[c] = st5[c].round(3)
st5.to_csv(TABLE_SUP_DIR / 'ST5_exposed_TFs_full_annotation.tsv',
           sep='\t', index=False)
print(f"  ST5: {len(st5)} exposed TFs written")

# ============================================================
# ST6: TCS pair analysis
# ============================================================
print("Generating ST6...")
tcs = load_tcs_pairs()
tcs.to_csv(TABLE_SUP_DIR / 'ST6_TCS_pair_analysis.tsv',
           sep='\t', index=False)
print(f"  ST6: {len(tcs)} TCS pairs written")

# ============================================================
# ST7: Simpson's paradox statistics
# ============================================================
print("Generating ST7...")
simpsons_data = []

# GCCGGC suppression (H17)
simpsons_data.append({
    'instance': 'GCCGGC_suppression',
    'hypothesis': 'H17',
    'description': 'GCCGGC proximity suppresses gene expression',
    'unstratified_p': 4.1e-10,
    'unstratified_effect': 'negative correlation',
    'core_only_p': 0.87,
    'core_only_effect': 'no effect',
    'arm_only_p': np.nan,
    'arm_only_effect': '',
    'confound': 'T1 GCCGGC 83% core; core genes lower LFC',
    'verdict': 'ARTIFACT',
})
# GCCGGC de-repression (H19)
simpsons_data.append({
    'instance': 'GCCGGC_derepression',
    'hypothesis': 'H19',
    'description': 'GCCGGC loss causes temporal de-repression',
    'unstratified_p': 8.3e-8,
    'unstratified_effect': 'Lost genes DOWN (median LFC = -0.259)',
    'core_only_p': 0.87,
    'core_only_effect': 'median LFC diff = +0.022 (no effect)',
    'arm_only_p': np.nan,
    'arm_only_effect': '',
    'confound': 'Lost genes enriched in core region',
    'verdict': 'ARTIFACT',
})
# BGC enrichment (H24)
simpsons_data.append({
    'instance': 'BGC_methylation_enrichment',
    'hypothesis': 'H24',
    'description': 'Methylation enriched in biosynthetic gene clusters',
    'unstratified_p': 4.4e-8,
    'unstratified_effect': 'fold = 1.66',
    'core_only_p': np.nan,
    'core_only_effect': 'fold = 1.07 (NS)',
    'arm_only_p': np.nan,
    'arm_only_effect': '',
    'confound': '100% co-localization of BGCs and T1 sites in core',
    'verdict': 'ARTIFACT',
})
# AAGCCCG (clean null, no confound)
simpsons_data.append({
    'instance': 'AAGCCCG_derepression',
    'hypothesis': 'H21',
    'description': 'AAGCCCG loss causes temporal de-repression',
    'unstratified_p': 0.91,
    'unstratified_effect': 'r = 0.003 (null)',
    'core_only_p': np.nan,
    'core_only_effect': 'No stratification needed (minimal geographic bias)',
    'arm_only_p': np.nan,
    'arm_only_effect': '',
    'confound': 'None (4.7 pp geographic shift only)',
    'verdict': 'GENUINE_NULL',
})

st7 = pd.DataFrame(simpsons_data)
st7.to_csv(TABLE_SUP_DIR / 'ST7_simpsons_paradox_statistics.tsv',
           sep='\t', index=False)
print(f"  ST7: {len(st7)} instances written")

# ============================================================
# ST8: Hypothesis ledger (H1–H36)
# ============================================================
print("Generating ST8...")
hypotheses = [
    ('1', 'H1', 'Partially supported', '4mC promoter-biased, 6mA gene-body-dispersed', 'Distribution differs (p=1.56e-14) but 6mA more promoter-proximal'),
    ('1', 'H2', 'Not tested', 'MTase expression stability causes 4mC/6mA difference', 'Deferred to H5'),
    ('1', 'H3', 'Premise rejected', 'absA de-repression compensates for redZ', '36/37 locus_tag mapping errors; true redZ upregulated (LFC=+0.91)'),
    ('2', 'H4', 'Rejected', 'Literature TFs show methylation-expression coordination', '3/37 methylated, 0/37 coordinated'),
    ('2', 'H5', 'Partially supported', 'MTase stability correlates with methylation-expression', 'Site-level confirmed, gene-level weak (|rho|<0.16)'),
    ('3', 'H6', 'Supported', 'Genome-wide regulator screen finds coordinated genes', '62/1055 coordinated (all novel)'),
    ('3', 'H7', 'Rejected', 'GCCGGC 4mC decrease is passive dilution', 'Jaccard=0.000; active remodeling'),
    ('4', 'H8', 'Partially supported', 'Coordinated regulators arm-enriched, secondary metabolism', 'T3 discordant_gain_up 87% arm (OR=8.05)'),
    ('4', 'H9', 'Rejected', 'GCCGGC "4mC" is misclassified 5mC', '5mC 0.01%; 4mC 82.7%; 91.4% in GCCGGC'),
    ('5', 'H10', 'Partially supported', 'REBASE identifies GCCGGC N4-C MTase', '2/87 NaeI produce m4C; SC_RS24685 top by annotation'),
    ('6', 'H12', 'Rejected', 'BLAST identifies GCCGGC MTase as SC_RS24685', 'SC_RS19770 top hit (E=0.007)'),
    ('6', 'H13', 'Partially supported', 'AAGCCCG non-random, secondary metabolism-associated', 'Non-random confirmed; regulatory depleted (fold=0.43)'),
    ('7', 'H14', 'Partially supported', 'Defense island produces T3 arm GCCGGC methylation', 'Geographic shift (chi2=597); 1 T3 site in island'),
    ('7', 'H15', 'Supported', 'Cross-motif regulatory gene methylation avoidance', 'All motifs: fold=0.43-0.62, p<0.02'),
    ('8', 'H16', 'Partially supported', 'Expression dynamics reverse-ID GCCGGC MTase', 'SC_RS13615 top composite; no candidate satisfies all criteria'),
    ('8', 'H17', 'Rejected', 'GCCGGC and AAGCCCG show transcriptional division of labor', "GCCGGC suppression = geographic confound (Simpson's paradox)"),
    ('9', 'H18', 'Partially supported', 'GCCGGC dose-response on expression', 'Threshold switch (0 vs >=1), not graded; geographic confound'),
    ('9', 'H19', 'Rejected', 'GCCGGC loss causes temporal de-repression', 'core-only p=0.87; Lost genes DOWN (median LFC=-0.253)'),
    ('10', 'H20', 'Supported', 'H15 regulatory avoidance survives geographic stratification', 'CMH-adjusted p<0.007, core+arm both significant'),
    ('10', 'H21', 'Rejected', 'AAGCCCG loss causes temporal de-repression', 'p=0.91, r=0.003; clean null'),
    ('11', 'H22', 'Partially supported', 'R-M motif sequences depleted from regulatory gene DNA', 'Gene body fold=0.74-0.76; promoter NOT depleted'),
    ('11', 'H23', 'Rejected', 'Methylation avoidance extends to all essential categories', 'Only regulatory genes; ribosomal, DNA repair NS'),
    ('12', 'H24', 'Partially supported', 'BGC methylation enrichment is geographic confound', 'core-only fold=1.07 (NS)'),
    ('12', 'H25', 'Supported', 'TSS methylation spatial gradient at regulatory genes', '2,200 bp zone, +300 bp deepest (ratio=0.541)'),
    ('13', 'H26', 'Rejected', 'Individual TF binding sites show methylation depletion', 'Flat profiles; GCCGGC enriched at TFBS (fold=1.16)'),
    ('13', 'H27', 'Supported', '62 coordinated regulators have absent protection zones', '8.4x TSS enrichment, p=5.3e-28'),
    ('14', 'H28', 'Supported', 'Exposed regulators: low T1, 100% dynamic, high GC', 'constitutive=0% (OR=inf, p=5.3e-13)'),
    ('14', 'H29', 'Rejected', 'baseMean threshold separates shielded/exposed', 'baseMean AUC=0.547; distance AUC=0.917'),
    ('15', 'H30', 'Partially supported', 'TSS sequence features determine protection zone', 'AAGCCCG 5.0x enriched; combined CV AUC=0.712'),
    ('15', 'H31', 'Partially supported', 'Exposed TFs regulate downstream cascade', '0/62 have FIMO motifs; coordination is intrinsic'),
    ('16', 'H32', 'Partially supported', 'Exposed TFs form self-regulatory module', 'Not physical module; within-type rho=0.500; 2 blocs'),
    ('16', 'H33', 'Rejected', 'Exposed TF neighborhoods show correlated expression', 'Permutation p=0.857; no distance decay'),
    ('16', 'H35-TF', 'Supported', 'TF family predicts bloc membership', 'TetR OR=0.16 (p=0.028)'),
    ('17', 'H34', 'Partially supported', 'Temporal phase separation between blocs', 'No separation (p=0.459); 63% early; simultaneous'),
    ('17', 'H35-causal', 'Rejected', 'AAGCCCG methylation mediates exposed TF regulation', '2/62 have methylated AAGCCCG; sequence != methylation'),
    ('18', 'H36', 'Partially supported', 'Exposed TF evolutionary conservation', 'Equally conserved (composite p=0.276, AUC=0.459)'),
]
st8 = pd.DataFrame(hypotheses, columns=['loop', 'hypothesis_id', 'verdict', 'description', 'key_statistic'])
st8.to_csv(TABLE_SUP_DIR / 'ST8_hypothesis_ledger.tsv',
           sep='\t', index=False)
print(f"  ST8: {len(st8)} hypotheses written")

print("\nAll supplementary tables generated successfully.")
print(f"Output directory: {TABLE_SUP_DIR}")
