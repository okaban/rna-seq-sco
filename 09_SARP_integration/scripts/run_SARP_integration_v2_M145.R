#!/usr/bin/env Rscript
# 09_SARP_integration v2 — M145
# Integrates functional SARP whitelist (actII-orf4) with v1 structural SARPs

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

cat("=== 09_SARP_integration v2 started at", format(Sys.time()), "===\n")

# -------------------------------------------------------------------
# 0. Paths
# -------------------------------------------------------------------
SARP_ROOT    <- "/Users/okaban/bioinfo/rna-seq/09_SARP_integration"
ANNOT_ROOT   <- "/Users/okaban/bioinfo/rna-seq/05_annotation"
REG_ROOT     <- "/Users/okaban/bioinfo/rna-seq/07_regulator_network"
TF_ROOT      <- "/Users/okaban/bioinfo/rna-seq/08_candidate_TF_prioritization"

ANNOT_RUN  <- file.path(ANNOT_ROOT, "analysis/05_annotation_260128_v1")
REG_RUN    <- file.path(REG_ROOT,   "analysis/07_regulator_network_260128_v1")
TF_RUN     <- file.path(TF_ROOT,    "analysis/08_candidate_TF_260128_v1")
V1_RUN     <- file.path(SARP_ROOT,  "analysis/09_SARP_integration_260128_v1")
V2_RUN     <- file.path(SARP_ROOT,  "analysis/09_SARP_integration_260128_v2")

GENE_MASTER_FILE   <- file.path(ANNOT_RUN, "tables/gene_master_with_BGC_regulators.tsv")
REG_CORR_FILE      <- file.path(REG_RUN,   "tables/regulator_BGC_correlation.tsv")
TF_REG_MASTER_FILE <- file.path(TF_RUN,    "tables/regulator_master_for_prioritization.tsv")

SARP_V1_LIST     <- file.path(V1_RUN, "tables/SARP_list_M145.tsv")
SARP_V1_BGC_MAP  <- file.path(V1_RUN, "tables/SARP_BGC_mapping.tsv")
SARP_V1_ACTIVITY <- file.path(V1_RUN, "tables/SARP_BGC_activity_summary.tsv")

OUT_TABLES <- file.path(V2_RUN, "tables")

# -------------------------------------------------------------------
# 1. Read inputs
# -------------------------------------------------------------------
cat("Reading input files...\n")

sarp_v1    <- read_tsv(SARP_V1_LIST, show_col_types = FALSE)
gene_master <- read_tsv(GENE_MASTER_FILE, show_col_types = FALSE)
reg_corr   <- read_tsv(REG_CORR_FILE, show_col_types = FALSE)
reg_master <- read_tsv(TF_REG_MASTER_FILE, show_col_types = FALSE)
bgc_map_v1 <- read_tsv(SARP_V1_BGC_MAP, show_col_types = FALSE)
activity_v1 <- read_tsv(SARP_V1_ACTIVITY, show_col_types = FALSE)

cat("  v1 SARP list:", nrow(sarp_v1), "rows\n")
cat("  gene_master:", nrow(gene_master), "rows\n")
cat("  regulator_master:", nrow(reg_master), "rows\n")
cat("  regulator_BGC_correlation:", nrow(reg_corr), "rows\n")

# -------------------------------------------------------------------
# 2. Define functional SARP whitelist
# -------------------------------------------------------------------
cat("\n--- Step 2: Define functional SARP whitelist ---\n")

functional_SARP_whitelist <- tibble::tribble(
  ~gene_id,      ~old_locus_tag, ~gene_name,     ~notes,
  "SC_RS27570",  "SCO5082",      "actII-orf4",   "Founder small SARP CSR for ACT; domain detection failed in v1 but functionally defined as SARP in literature",
  "SC_RS31630",  "SCO5877",      "redD",         "small SARP CSR for RED",
  "SC_RS18200",  "SCO3217",      NA_character_,  "medium SARP CSR for CDA (literature name: cdaR)",
  "SC_RS33690",  "SCO6288",      NA_character_,  "small SARP CSR for CPK (literature name: cpkN)"
)

cat("  Functional SARP whitelist:", nrow(functional_SARP_whitelist), "genes\n")
cat("  Whitelist gene_ids:", paste(functional_SARP_whitelist$gene_id, collapse = ", "), "\n")

# -------------------------------------------------------------------
# 3. Merge v1 structural SARPs + functional whitelist → SARP_list_v2
# -------------------------------------------------------------------
cat("\n--- Step 3: Merge structural + functional → SARP_list_v2 ---\n")

# Mark which genes are in v1 (structural)
structural_ids <- sarp_v1$gene_id
functional_ids <- functional_SARP_whitelist$gene_id

# actII-orf4 is NOT in v1 structural list; the other 3 CSRs ARE
cat("  v1 structural IDs:", paste(structural_ids, collapse = ", "), "\n")
cat("  functional whitelist IDs:", paste(functional_ids, collapse = ", "), "\n")
cat("  Overlap:", paste(intersect(structural_ids, functional_ids), collapse = ", "), "\n")
cat("  New from whitelist:", paste(setdiff(functional_ids, structural_ids), collapse = ", "), "\n")

# For actII-orf4, we need to add its info from gene_master
actii_gm <- gene_master %>%
  filter(gene_id == "SC_RS27570") %>%
  select(gene_id, old_locus_tag, gene_name, product, protein_id)

# Build actII-orf4 row in the same format as v1 SARP list
actii_row <- tibble(
  gene_id       = "SC_RS27570",
  old_locus_tag = "SCO5082",
  gene_name     = actii_gm$gene_name,
  product       = actii_gm$product,
  protein_id    = actii_gm$protein_id,
  protein_len   = 259L,         # verified from protein.faa
  has_DBD       = FALSE,
  has_BTAD      = FALSE,
  has_NB_ARC    = FALSE,
  has_TPR       = FALSE,
  SARP_subtype  = "small",      # literature-based: ActII-ORF4 is a small SARP
  domains_list  = NA_character_,
  best_DBD_evalue  = NA_real_,
  best_BTAD_evalue = NA_real_
)

# Combine v1 structural SARPs + actII-orf4
sarp_v2_base <- bind_rows(sarp_v1, actii_row)

# Add classification flags
sarp_v2 <- sarp_v2_base %>%
  mutate(
    is_SARP_structural = gene_id %in% structural_ids,
    is_SARP_functional = gene_id %in% functional_ids,
    SARP_category = case_when(
      is_SARP_structural & is_SARP_functional  ~ "both",
      is_SARP_structural & !is_SARP_functional ~ "structural_only",
      !is_SARP_structural & is_SARP_functional ~ "functional_only",
      TRUE ~ NA_character_
    )
  )

cat("\nSARP_list_v2 summary:\n")
cat("  Total SARPs:", nrow(sarp_v2), "\n")
cat("  Categories:\n")
print(table(sarp_v2$SARP_category))

# Reorder: gene_name for actII-orf4 is "actII" in GFF; override to "actII-orf4" for clarity
sarp_v2 <- sarp_v2 %>%
  mutate(gene_name = ifelse(gene_id == "SC_RS27570", "actII-orf4", gene_name))

# Save
write_tsv(sarp_v2, file.path(OUT_TABLES, "SARP_list_M145_v2.tsv"))
cat("  Saved: SARP_list_M145_v2.tsv\n")

# Print full table
cat("\nSARP_list_M145_v2:\n")
sarp_v2 %>%
  select(gene_id, old_locus_tag, gene_name, SARP_subtype, is_SARP_structural, is_SARP_functional, SARP_category) %>%
  print(n = Inf, width = 200)

# -------------------------------------------------------------------
# 4. Integrate into gene_master → gene_master_with_SARP_v2
# -------------------------------------------------------------------
cat("\n--- Step 4: Integrate into gene_master ---\n")

sarp_flags <- sarp_v2 %>%
  select(gene_id, is_SARP_structural, is_SARP_functional, SARP_category, SARP_subtype_v2 = SARP_subtype)

gene_master_v2 <- gene_master %>%
  left_join(sarp_flags, by = "gene_id") %>%
  mutate(
    is_SARP_structural = replace_na(is_SARP_structural, FALSE),
    is_SARP_functional = replace_na(is_SARP_functional, FALSE)
  )

cat("  gene_master_v2:", nrow(gene_master_v2), "rows\n")
cat("  is_SARP_structural TRUE:", sum(gene_master_v2$is_SARP_structural), "\n")
cat("  is_SARP_functional TRUE:", sum(gene_master_v2$is_SARP_functional), "\n")

write_tsv(gene_master_v2, file.path(OUT_TABLES, "gene_master_with_SARP_v2.tsv"))
cat("  Saved: gene_master_with_SARP_v2.tsv\n")

# -------------------------------------------------------------------
# 5. Integrate into regulator_master → regulator_master_with_SARP_v2
# -------------------------------------------------------------------
cat("\n--- Step 5: Integrate into regulator_master ---\n")

reg_master_v2 <- reg_master %>%
  left_join(sarp_flags, by = "gene_id") %>%
  mutate(
    is_SARP_structural = replace_na(is_SARP_structural, FALSE),
    is_SARP_functional = replace_na(is_SARP_functional, FALSE)
  )

cat("  regulator_master_v2:", nrow(reg_master_v2), "rows\n")
cat("  is_SARP_structural TRUE:", sum(reg_master_v2$is_SARP_structural), "\n")
cat("  is_SARP_functional TRUE:", sum(reg_master_v2$is_SARP_functional), "\n")

# Show all SARPs in regulator_master
cat("\nSARP entries in regulator_master_v2:\n")
reg_master_v2 %>%
  filter(is_SARP_structural | is_SARP_functional) %>%
  select(gene_id, old_locus_tag, gene_name, SARP_subtype_v2, SARP_category, corr_act, corr_red, corr_cda, corr_cpk, phase_preference) %>%
  print(n = Inf, width = 200)

write_tsv(reg_master_v2, file.path(OUT_TABLES, "regulator_master_with_SARP_v2.tsv"))
cat("  Saved: regulator_master_with_SARP_v2.tsv\n")

# -------------------------------------------------------------------
# 6. SARP_BGC_mapping_v2
# -------------------------------------------------------------------
cat("\n--- Step 6: SARP_BGC_mapping_v2 ---\n")

# Start from v1 BGC mapping columns structure, add actII-orf4
bgc_map_v1_cols <- names(bgc_map_v1)

# Build v2 mapping: for structural SARPs, reuse v1; for actII-orf4, add manually
bgc_map_v2 <- bgc_map_v1 %>%
  mutate(
    is_SARP_structural = TRUE,
    is_SARP_functional = gene_id %in% functional_ids,
    SARP_category = ifelse(is_SARP_functional, "both", "structural_only")
  )

# Add actII-orf4 row
actii_bgc <- tibble(
  gene_id        = "SC_RS27570",
  old_locus_tag  = "SCO5082",
  gene_name      = "actII-orf4",
  product        = "TetR family transcriptional regulator ActII",
  protein_id     = "WP_011030045.1",
  protein_len    = 259L,
  has_DBD        = FALSE,
  has_BTAD       = FALSE,
  has_NB_ARC     = FALSE,
  has_TPR        = FALSE,
  SARP_subtype   = "small",
  domains_list   = NA_character_,
  best_DBD_evalue  = NA_real_,
  best_BTAD_evalue = NA_real_,
  bgc_name       = "act",
  role_in_BGC    = "CSR",
  is_SARP_structural = FALSE,
  is_SARP_functional = TRUE,
  SARP_category  = "functional_only"
)

bgc_map_v2 <- bind_rows(bgc_map_v2, actii_bgc)

cat("  SARP_BGC_mapping_v2:", nrow(bgc_map_v2), "rows\n")
cat("\nMapping summary:\n")
bgc_map_v2 %>%
  select(gene_id, old_locus_tag, gene_name, SARP_subtype, bgc_name, role_in_BGC, SARP_category) %>%
  print(n = Inf, width = 200)

write_tsv(bgc_map_v2, file.path(OUT_TABLES, "SARP_BGC_mapping_v2.tsv"))
cat("  Saved: SARP_BGC_mapping_v2.tsv\n")

# -------------------------------------------------------------------
# 7. SARP_BGC_activity_summary_v2
# -------------------------------------------------------------------
cat("\n--- Step 7: SARP_BGC_activity_summary_v2 ---\n")

# Start from v1 activity, add v2 flags
activity_v2 <- activity_v1 %>%
  mutate(
    is_SARP_structural = TRUE,
    is_SARP_functional = gene_id %in% functional_ids,
    SARP_category = ifelse(is_SARP_functional, "both", "structural_only")
  )

# Build actII-orf4 activity row from regulator_master + regulator_BGC_correlation
actii_reg <- reg_master %>% filter(gene_id == "SC_RS27570")
actii_corr <- reg_corr %>% filter(gene_id == "SC_RS27570")

actii_activity <- tibble(
  gene_id        = "SC_RS27570",
  old_locus_tag  = "SCO5082",
  gene_name      = "actII-orf4",
  product        = "TetR family transcriptional regulator ActII",
  protein_len    = 259L,
  SARP_subtype   = "small",
  bgc_name       = "act",
  role_in_BGC    = "CSR",
  corr_act       = actii_corr$corr_act,
  corr_red       = actii_corr$corr_red,
  corr_cda       = actii_corr$corr_cda,
  corr_cpk       = actii_corr$corr_cpk,
  abs_corr_act   = abs(actii_corr$corr_act),
  abs_corr_red   = abs(actii_corr$corr_red),
  abs_corr_cda   = abs(actii_corr$corr_cda),
  abs_corr_cpk   = abs(actii_corr$corr_cpk),
  n_high_corr    = sum(c(abs(actii_corr$corr_act), abs(actii_corr$corr_red),
                          abs(actii_corr$corr_cda), abs(actii_corr$corr_cpk)) > 0.7),
  mean_abs_corr  = mean(c(abs(actii_corr$corr_act), abs(actii_corr$corr_red),
                           abs(actii_corr$corr_cda), abs(actii_corr$corr_cpk))),
  delta_2_vs_1   = actii_reg$delta_2_vs_1,
  delta_3_vs_2   = actii_reg$delta_3_vs_2,
  delta_3_vs_1   = actii_reg$delta_3_vs_1,
  log2FC_3_vs_1  = actii_reg$log2FC_3_vs_1,
  padj_3_vs_1    = actii_reg$padj_3_vs_1,
  phase_preference = actii_reg$phase_preference,
  is_SARP_structural = FALSE,
  is_SARP_functional = TRUE,
  SARP_category  = "functional_only"
)

activity_v2 <- bind_rows(activity_v2, actii_activity)

cat("  SARP_BGC_activity_summary_v2:", nrow(activity_v2), "rows\n")
cat("\nActivity summary:\n")
activity_v2 %>%
  select(gene_id, old_locus_tag, gene_name, SARP_subtype, bgc_name, role_in_BGC,
         corr_act, corr_red, corr_cda, corr_cpk, mean_abs_corr,
         log2FC_3_vs_1, phase_preference, SARP_category) %>%
  print(n = Inf, width = 200)

write_tsv(activity_v2, file.path(OUT_TABLES, "SARP_BGC_activity_summary_v2.tsv"))
cat("  Saved: SARP_BGC_activity_summary_v2.tsv\n")

# -------------------------------------------------------------------
# Summary
# -------------------------------------------------------------------
cat("\n=== 09_SARP_integration v2 completed at", format(Sys.time()), "===\n")
cat("\nOutput files:\n")
for (f in list.files(OUT_TABLES, full.names = TRUE)) {
  cat("  ", f, " (", file.size(f), " bytes)\n")
}
