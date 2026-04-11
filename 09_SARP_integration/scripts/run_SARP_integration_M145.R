#!/usr/bin/env Rscript
# 09_SARP_integration: SARP family identification, classification, and network integration
# M145 RNA-seq (3 timepoints x 3 replicates = 9 samples)

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(pheatmap)
  library(svglite)
})

# ── Paths ──────────────────────────────────────────────────────────────────
BASE       <- "/Users/okaban/bioinfo/rna-seq"
ANNOT_RUN  <- file.path(BASE, "05_annotation/analysis/05_annotation_260128_v1")
BGC_RUN    <- file.path(BASE, "06_BGC_dynamics/analysis/06_BGC_dynamics_260128_v1")
REG_RUN    <- file.path(BASE, "07_regulator_network/analysis/07_regulator_network_260128_v1")
TF_RUN     <- file.path(BASE, "08_candidate_TF_prioritization/analysis/08_candidate_TF_260128_v1")
DESEQ_RUN  <- file.path(BASE, "04_deseq2/analysis/04_deseq2_260128_v1")
OUT_DIR    <- file.path(BASE, "09_SARP_integration/analysis/09_SARP_integration_260128_v1")
REF_DIR    <- "/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1"

LOG_FILE   <- file.path(OUT_DIR, "logs/SARP_integration_pipeline.log")
sink(LOG_FILE, split = TRUE)

cat("=== 09_SARP_integration started at", format(Sys.time()), "===\n")

# ── 1. Parse hmmsearch domain table output ────────────────────────────────
cat("\n--- Parsing hmmsearch domain hits ---\n")

domtbl_file <- file.path(OUT_DIR, "tables/hmmsearch_domtblout.txt")
domtbl_raw <- readLines(domtbl_file)
domtbl_data <- domtbl_raw[!grepl("^#", domtbl_raw)]
domtbl_data <- domtbl_data[nchar(trimws(domtbl_data)) > 0]

# Parse the domain table (space-delimited, last column is description)
parse_domtbl <- function(lines) {
  lapply(lines, function(line) {
    fields <- strsplit(trimws(line), "\\s+")[[1]]
    data.frame(
      protein_id = fields[1],
      protein_len = as.integer(fields[3]),
      domain_name = fields[4],
      domain_acc = fields[5],
      domain_len = as.integer(fields[6]),
      seq_evalue = as.numeric(fields[7]),
      seq_score = as.numeric(fields[8]),
      dom_evalue = as.numeric(fields[13]),
      dom_score = as.numeric(fields[14]),
      ali_from = as.integer(fields[18]),
      ali_to = as.integer(fields[19]),
      env_from = as.integer(fields[20]),
      env_to = as.integer(fields[21]),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
}

dom_hits <- parse_domtbl(domtbl_data)
cat("Total domain hits:", nrow(dom_hits), "\n")
cat("Domain types found:\n")
print(table(dom_hits$domain_name))

# Save raw domain hits
write_tsv(dom_hits, file.path(OUT_DIR, "tables/SARP_domain_hits_raw.tsv"))
cat("Saved: SARP_domain_hits_raw.tsv\n")

# ── 2. Map protein_id → gene_id using GFF/master table ───────────────────
cat("\n--- Mapping protein_id to gene_id ---\n")

gene_master <- read_tsv(file.path(ANNOT_RUN, "tables/gene_master_with_BGC_regulators.tsv"),
                        show_col_types = FALSE)
cat("Gene master:", nrow(gene_master), "genes\n")

# protein_id column
prot_map <- gene_master %>%
  filter(!is.na(protein_id), protein_id != "") %>%
  select(protein_id, gene_id, old_locus_tag, gene_name, product)
cat("Proteins with IDs in master:", nrow(prot_map), "\n")

# ── 3. Identify SARP candidates (DBD + BTAD) ─────────────────────────────
cat("\n--- Identifying SARP candidates ---\n")

# Summarize domains per protein
protein_domains <- dom_hits %>%
  group_by(protein_id) %>%
  summarise(
    protein_len = dplyr::first(protein_len),
    has_DBD    = any(domain_name == "Trans_reg_C"),
    has_BTAD   = any(domain_name == "BTAD"),
    has_NB_ARC = any(domain_name == "NB-ARC"),
    has_TPR    = any(domain_name == "TPR_12"),
    n_domains  = n_distinct(domain_name),
    domains_list = paste(sort(unique(domain_name)), collapse = ";"),
    best_DBD_evalue  = ifelse(any(domain_name == "Trans_reg_C"),
                               min(dom_evalue[domain_name == "Trans_reg_C"]), NA),
    best_BTAD_evalue = ifelse(any(domain_name == "BTAD"),
                               min(dom_evalue[domain_name == "BTAD"]), NA),
    .groups = "drop"
  )

cat("Proteins with any SARP-related domain:", nrow(protein_domains), "\n")
cat("  with DBD (Trans_reg_C):", sum(protein_domains$has_DBD), "\n")
cat("  with BTAD:", sum(protein_domains$has_BTAD), "\n")
cat("  with NB-ARC:", sum(protein_domains$has_NB_ARC), "\n")
cat("  with TPR_12:", sum(protein_domains$has_TPR), "\n")

# SARP candidates: at least DBD + BTAD
sarp_candidates <- protein_domains %>%
  filter(has_DBD & has_BTAD)
cat("SARP candidates (DBD + BTAD):", nrow(sarp_candidates), "\n")

# Also flag proteins with BTAD only (might be SARP-related)
btad_only <- protein_domains %>%
  filter(has_BTAD & !has_DBD)
cat("BTAD-only proteins (no DBD):", nrow(btad_only), "\n")

# Map to gene IDs
sarp_with_genes <- sarp_candidates %>%
  left_join(prot_map, by = "protein_id")
cat("\nSARP candidates mapped to genes:\n")
print(sarp_with_genes %>% select(protein_id, gene_id, old_locus_tag, gene_name, product,
                                  protein_len, has_NB_ARC, has_TPR, domains_list))

# Save candidates domain-based
write_tsv(sarp_with_genes %>%
            select(protein_id, gene_id, old_locus_tag, gene_name, product, protein_len,
                   has_DBD, has_BTAD, has_NB_ARC, has_TPR, domains_list,
                   best_DBD_evalue, best_BTAD_evalue),
          file.path(OUT_DIR, "tables/SARP_candidates_domain_based.tsv"))
cat("Saved: SARP_candidates_domain_based.tsv\n")

# Also save BTAD-only for reference
btad_with_genes <- btad_only %>%
  left_join(prot_map, by = "protein_id")
cat("\nBTAD-only proteins:\n")
print(btad_with_genes %>% select(protein_id, gene_id, old_locus_tag, gene_name, product,
                                  protein_len, domains_list))

# ── 4. SARP subtype classification ───────────────────────────────────────
cat("\n--- SARP subtype classification ---\n")

# Classification rules:
# small:    DBD + BTAD, no NB-ARC, no TPR, length ~250-350 (use <500)
# medium:   DBD + BTAD + NB-ARC, no TPR, length ~550-650 (use 400-800)
# large:    DBD + BTAD + NB-ARC + TPR, length ~900-1100 (use >800)
# SARP-LAL: DBD + BTAD, large (>800), potentially LAL domain (check later)

sarp_classified <- sarp_with_genes %>%
  mutate(
    SARP_subtype = case_when(
      has_NB_ARC & has_TPR                    ~ "large",
      has_NB_ARC & !has_TPR                   ~ "medium",
      !has_NB_ARC & !has_TPR & protein_len <= 500 ~ "small",
      !has_NB_ARC & !has_TPR & protein_len > 500  ~ "small_long",
      TRUE ~ "unclassified"
    )
  )

cat("Subtype distribution:\n")
print(table(sarp_classified$SARP_subtype))
cat("\nSARP list:\n")
print(sarp_classified %>%
        select(gene_id, old_locus_tag, gene_name, product, protein_len,
               SARP_subtype, domains_list))

# Save SARP_list_M145.tsv
sarp_list <- sarp_classified %>%
  select(gene_id, old_locus_tag, gene_name, product, protein_id, protein_len,
         has_DBD, has_BTAD, has_NB_ARC, has_TPR, SARP_subtype, domains_list,
         best_DBD_evalue, best_BTAD_evalue)

write_tsv(sarp_list, file.path(OUT_DIR, "tables/SARP_list_M145.tsv"))
cat("Saved: SARP_list_M145.tsv\n")

# ── 5. Check actII-orf4 and cpkO/kasO domain status ─────────────────────
cat("\n--- Checking actII-orf4 and cpkO/kasO ---\n")

# actII-orf4: SC_RS27570 / SCO5082 — GFF says TetR, literature says SARP
actII_prot <- prot_map %>% filter(gene_id == "SC_RS27570")
cat("actII-orf4 protein_id:", actII_prot$protein_id, "\n")
actII_dom <- protein_domains %>% filter(protein_id == actII_prot$protein_id)
if (nrow(actII_dom) > 0) {
  cat("actII-orf4 domains found:", actII_dom$domains_list, "\n")
  cat("  has_DBD:", actII_dom$has_DBD, " has_BTAD:", actII_dom$has_BTAD, "\n")
} else {
  cat("actII-orf4: NO SARP-related domains detected by hmmsearch\n")
}

# cpkO/kasO: SC_RS33660 / SCO6282 — GFF says SDR
cpkO_prot <- prot_map %>% filter(gene_id == "SC_RS33660")
cat("cpkO/kasO protein_id:", cpkO_prot$protein_id, "\n")
cpkO_dom <- protein_domains %>% filter(protein_id == cpkO_prot$protein_id)
if (nrow(cpkO_dom) > 0) {
  cat("cpkO/kasO domains found:", cpkO_dom$domains_list, "\n")
} else {
  cat("cpkO/kasO: NO SARP-related domains detected by hmmsearch\n")
}

# ── 6. Integrate into gene_master ─────────────────────────────────────────
cat("\n--- Integrating SARP info into gene_master ---\n")

sarp_ids <- sarp_list %>%
  select(gene_id, is_SARP = has_DBD, SARP_subtype)
sarp_ids$is_SARP <- TRUE

gene_master_sarp <- gene_master %>%
  left_join(sarp_ids %>% select(gene_id, is_SARP, SARP_subtype), by = "gene_id") %>%
  mutate(is_SARP = ifelse(is.na(is_SARP), FALSE, is_SARP),
         SARP_subtype = ifelse(is.na(SARP_subtype), NA_character_, SARP_subtype))

cat("Gene master with SARP:", nrow(gene_master_sarp), "genes,",
    sum(gene_master_sarp$is_SARP), "SARPs\n")

write_tsv(gene_master_sarp, file.path(OUT_DIR, "tables/gene_master_with_SARP.tsv"))
cat("Saved: gene_master_with_SARP.tsv\n")

# Integrate into regulator_master_for_prioritization
tf_master <- read_tsv(file.path(TF_RUN, "tables/regulator_master_for_prioritization.tsv"),
                      show_col_types = FALSE)

reg_master_sarp <- tf_master %>%
  left_join(sarp_ids %>% select(gene_id, is_SARP, SARP_subtype), by = "gene_id") %>%
  mutate(is_SARP = ifelse(is.na(is_SARP), FALSE, is_SARP),
         SARP_subtype = ifelse(is.na(SARP_subtype), NA_character_, SARP_subtype))

cat("Regulator master with SARP:", nrow(reg_master_sarp), "regulators,",
    sum(reg_master_sarp$is_SARP), "SARPs\n")

write_tsv(reg_master_sarp, file.path(OUT_DIR, "tables/regulator_master_with_SARP.tsv"))
cat("Saved: regulator_master_with_SARP.tsv\n")

# ── 7. SARP-BGC mapping ──────────────────────────────────────────────────
cat("\n--- SARP-BGC mapping ---\n")

bgc_def <- read_tsv(file.path(ANNOT_RUN, "tables/BGC_definition_manual.tsv"),
                    show_col_types = FALSE)
bgc_genes <- bgc_def %>% select(bgc_name, gene_id)

sarp_bgc <- sarp_list %>%
  left_join(bgc_genes, by = "gene_id") %>%
  mutate(
    bgc_name = ifelse(is.na(bgc_name), "none", bgc_name),
    role_in_BGC = ifelse(bgc_name == "none", "global", "CSR")
  )

cat("SARP-BGC mapping:\n")
print(sarp_bgc %>% select(gene_id, old_locus_tag, gene_name, product,
                            SARP_subtype, bgc_name, role_in_BGC))

write_tsv(sarp_bgc, file.path(OUT_DIR, "tables/SARP_BGC_mapping.tsv"))
cat("Saved: SARP_BGC_mapping.tsv\n")

cat("\nSummary: CSR vs Global:\n")
print(table(sarp_bgc$role_in_BGC))
cat("\nBy BGC:\n")
print(table(sarp_bgc$bgc_name))

# ── 8. SARP expression & BGC correlation summary ─────────────────────────
cat("\n--- SARP expression & BGC correlation summary ---\n")

reg_corr <- read_tsv(file.path(REG_RUN, "tables/regulator_BGC_correlation.tsv"),
                     show_col_types = FALSE)

# Get phase deltas from regulator_master_for_prioritization (step 08)
phase_data <- tf_master %>%
  select(gene_id, delta_2_vs_1, delta_3_vs_2, delta_3_vs_1,
         log2FC_3_vs_1, padj_3_vs_1, phase_preference) %>%
  distinct()

# Build SARP activity summary
sarp_activity <- sarp_bgc %>%
  select(gene_id, old_locus_tag, gene_name, product, protein_len,
         SARP_subtype, bgc_name, role_in_BGC) %>%
  left_join(reg_corr %>% select(gene_id, corr_act, corr_red, corr_cda, corr_cpk,
                                 abs_corr_act, abs_corr_red, abs_corr_cda, abs_corr_cpk,
                                 n_high_corr, mean_abs_corr),
            by = "gene_id") %>%
  left_join(phase_data, by = "gene_id")

cat("SARP activity summary:\n")
print(sarp_activity %>%
        select(gene_id, old_locus_tag, SARP_subtype, bgc_name, role_in_BGC,
               corr_act, corr_red, corr_cda, corr_cpk,
               log2FC_3_vs_1, delta_2_vs_1, delta_3_vs_2, phase_preference))

write_tsv(sarp_activity, file.path(OUT_DIR, "tables/SARP_BGC_activity_summary.tsv"))
cat("Saved: SARP_BGC_activity_summary.tsv\n")

# ── 9. Integration with step 08 priority list ────────────────────────────
cat("\n--- Integration with step 08 experimental priority ---\n")

exp_priority <- read_tsv(file.path(TF_RUN, "tables/TF_candidates_experimental_priority.tsv"),
                         show_col_types = FALSE)

# Check which priority candidates are SARPs
priority_sarp <- exp_priority %>%
  left_join(sarp_ids, by = "gene_id") %>%
  mutate(is_SARP = ifelse(is.na(is_SARP), FALSE, is_SARP))

cat("Step 08 priority candidates with SARP status:\n")
print(priority_sarp %>% select(gene_id, old_locus_tag, gene_name, target_BGC, role,
                                is_known, is_SARP, SARP_subtype))

# Identify interesting SARPs NOT in priority list
novel_sarp_candidates <- sarp_activity %>%
  filter(!gene_id %in% exp_priority$gene_id) %>%
  filter(role_in_BGC == "global") %>%
  filter(!is.na(mean_abs_corr)) %>%
  arrange(desc(mean_abs_corr))

cat("\nNovel global SARPs not in priority list:\n")
print(novel_sarp_candidates %>%
        select(gene_id, old_locus_tag, gene_name, product, SARP_subtype,
               corr_act, corr_red, corr_cda, corr_cpk, mean_abs_corr,
               log2FC_3_vs_1, phase_preference))

# Build extended priority list with SARP additions
sarp_additions <- novel_sarp_candidates %>%
  filter(mean_abs_corr >= 0.5) %>%
  head(3) %>%
  mutate(
    target_BGC = "SARP_global",
    role = case_when(
      !is.na(corr_red) & corr_red > 0 & !is.na(corr_cda) & corr_cda > 0 ~ "candidate_SARP_activator",
      !is.na(corr_red) & corr_red < 0 & !is.na(corr_cda) & corr_cda < 0 ~ "candidate_SARP_repressor",
      TRUE ~ "candidate_SARP_global"
    ),
    is_known = FALSE,
    known_name = NA_character_,
    comment = paste0("SARP (", SARP_subtype, "); mean_abs_corr=", round(mean_abs_corr, 3),
                     "; corr_act=", round(corr_act, 3),
                     "; corr_red=", round(corr_red, 3),
                     "; phase=", phase_preference)
  ) %>%
  select(gene_id, old_locus_tag, gene_name, product, target_BGC, role,
         is_known, known_name, comment)

extended_priority <- bind_rows(
  exp_priority %>%
    left_join(sarp_ids, by = "gene_id") %>%
    mutate(comment = ifelse(!is.na(SARP_subtype),
                            paste0(comment, " [SARP:", SARP_subtype, "]"),
                            comment)) %>%
    select(names(exp_priority)),
  sarp_additions
)

cat("\nExtended priority list with SARP additions:\n")
print(extended_priority %>% select(gene_id, old_locus_tag, gene_name, role, comment))

write_tsv(extended_priority,
          file.path(OUT_DIR, "tables/TF_candidates_experimental_priority_with_SARP.tsv"))
cat("Saved: TF_candidates_experimental_priority_with_SARP.tsv\n")

# ── 10. Figures ───────────────────────────────────────────────────────────
cat("\n--- Creating figures ---\n")
fig_dir <- file.path(OUT_DIR, "figures")

# 10a. SARP BGC correlation heatmap
sarp_for_heatmap <- sarp_activity %>%
  filter(!is.na(corr_act)) %>%
  mutate(label = case_when(
    !is.na(gene_name) & gene_name != "" & gene_name != "NA" ~ paste0(gene_name, " (", old_locus_tag, ")"),
    TRUE ~ paste0(old_locus_tag, " [", SARP_subtype, "]")
  ))

if (nrow(sarp_for_heatmap) >= 2) {
  mat_corr <- sarp_for_heatmap %>%
    select(label, corr_act, corr_red, corr_cda, corr_cpk) %>%
    column_to_rownames("label") %>%
    as.matrix()
  colnames(mat_corr) <- c("act", "red", "cda", "cpk")

  # Annotation row
  ann_row <- sarp_for_heatmap %>%
    select(label, SARP_subtype, role_in_BGC) %>%
    column_to_rownames("label")

  ann_colors <- list(
    SARP_subtype = c("small" = "#E41A1C", "medium" = "#377EB8",
                      "large" = "#4DAF4A", "small_long" = "#FF7F00"),
    role_in_BGC = c("CSR" = "#984EA3", "global" = "#FF7F00")
  )

  pdf(file.path(fig_dir, "SARP_BGC_correlation_heatmap_M145.pdf"), width = 7, height = 6)
  pheatmap(mat_corr,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           breaks = seq(-1, 1, length.out = 101),
           cluster_cols = FALSE,
           annotation_row = ann_row,
           annotation_colors = ann_colors,
           main = "SARP - BGC correlation (Pearson r)",
           fontsize_row = 8)
  dev.off()

  svglite(file.path(fig_dir, "SARP_BGC_correlation_heatmap_M145.svg"), width = 7, height = 6)
  pheatmap(mat_corr,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           breaks = seq(-1, 1, length.out = 101),
           cluster_cols = FALSE,
           annotation_row = ann_row,
           annotation_colors = ann_colors,
           main = "SARP - BGC correlation (Pearson r)",
           fontsize_row = 8)
  dev.off()
  cat("Saved: SARP_BGC_correlation_heatmap_M145.pdf/.svg\n")
}

# 10b. SARP expression heatmap (delta values)
sarp_expr <- sarp_for_heatmap %>%
  filter(!is.na(delta_2_vs_1)) %>%
  select(label, delta_2_vs_1, delta_3_vs_2, log2FC_3_vs_1)

if (nrow(sarp_expr) >= 2) {
  mat_expr <- sarp_expr %>%
    column_to_rownames("label") %>%
    as.matrix()
  colnames(mat_expr) <- c("Δ(2-1)", "Δ(3-2)", "log2FC(3/1)")

  ann_row2 <- sarp_for_heatmap %>%
    filter(!is.na(delta_2_vs_1)) %>%
    select(label, SARP_subtype, role_in_BGC) %>%
    column_to_rownames("label")

  pdf(file.path(fig_dir, "SARP_expression_heatmap_M145.pdf"), width = 6, height = 6)
  pheatmap(mat_expr,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           cluster_cols = FALSE,
           annotation_row = ann_row2,
           annotation_colors = ann_colors,
           main = "SARP expression dynamics (rlog deltas)",
           fontsize_row = 8)
  dev.off()

  svglite(file.path(fig_dir, "SARP_expression_heatmap_M145.svg"), width = 6, height = 6)
  pheatmap(mat_expr,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           cluster_cols = FALSE,
           annotation_row = ann_row2,
           annotation_colors = ann_colors,
           main = "SARP expression dynamics (rlog deltas)",
           fontsize_row = 8)
  dev.off()
  cat("Saved: SARP_expression_heatmap_M145.pdf/.svg\n")
}

# ── 11. Summary statistics ────────────────────────────────────────────────
cat("\n--- Summary statistics ---\n")
cat("Total SARP candidates (DBD + BTAD):", nrow(sarp_list), "\n")
cat("Subtype distribution:\n")
print(table(sarp_list$SARP_subtype))
cat("\nCSR vs Global:\n")
print(table(sarp_bgc$role_in_BGC))
cat("\nCSRs by BGC:\n")
sarp_csr <- sarp_bgc %>% filter(role_in_BGC == "CSR")
print(sarp_csr %>% select(gene_id, old_locus_tag, gene_name, SARP_subtype, bgc_name))
cat("\nGlobal SARPs:\n")
sarp_global <- sarp_bgc %>% filter(role_in_BGC == "global")
print(sarp_global %>% select(gene_id, old_locus_tag, gene_name, SARP_subtype, product))

cat("\n=== 09_SARP_integration completed at", format(Sys.time()), "===\n")
cat("Output directory:", OUT_DIR, "\n")

sink()
