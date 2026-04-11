#!/usr/bin/env Rscript
# =============================================================================
# 07_regulator_network: Regulator–BGC covariation analysis for M145 RNA-seq
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
  library(svglite)
  library(ggrepel)
})

cat("=== 07_regulator_network started at", format(Sys.time()), "===\n")

# --- Paths -------------------------------------------------------------------
REG_RUN_DIR   <- "/Users/okaban/bioinfo/rna-seq/07_regulator_network/analysis/07_regulator_network_260128_v1"
DESEQ_RUN_DIR <- "/Users/okaban/bioinfo/rna-seq/04_deseq2/analysis/04_deseq2_260128_v1"
ANNOT_RUN_DIR <- "/Users/okaban/bioinfo/rna-seq/05_annotation/analysis/05_annotation_260128_v1"
BGC_DYN_DIR   <- "/Users/okaban/bioinfo/rna-seq/06_BGC_dynamics/analysis/06_BGC_dynamics_260128_v1"

RLD_FILE           <- file.path(DESEQ_RUN_DIR, "rds", "rld.rds")
GENE_MASTER_FILE   <- file.path(ANNOT_RUN_DIR, "tables", "gene_master_with_BGC_regulators.tsv")
BGC_SCORES_FILE    <- file.path(BGC_DYN_DIR, "tables", "BGC_sample_scores.tsv")

TABLE_DIR <- file.path(REG_RUN_DIR, "tables")
FIG_DIR   <- file.path(REG_RUN_DIR, "figures")

cat("RLD_FILE:", RLD_FILE, "\n")
cat("GENE_MASTER_FILE:", GENE_MASTER_FILE, "\n")
cat("BGC_SCORES_FILE:", BGC_SCORES_FILE, "\n\n")

# =============================================================================
# 1. Load data
# =============================================================================
cat("--- Loading data ---\n")

# 1a. rlog object
rld <- readRDS(RLD_FILE)
rlog_mat <- assay(rld)
cat("rlog matrix:", nrow(rlog_mat), "genes x", ncol(rlog_mat), "samples\n")
cat("Sample names:", paste(colnames(rlog_mat), collapse = ", "), "\n")

# 1b. Master table
master <- read_tsv(GENE_MASTER_FILE, show_col_types = FALSE)
cat("Master table:", nrow(master), "genes\n")

# 1c. BGC sample scores
bgc_scores <- read_tsv(BGC_SCORES_FILE, show_col_types = FALSE)
cat("BGC scores:", nrow(bgc_scores), "samples\n")
print(as.data.frame(bgc_scores))

# Sample info
sample_ids <- colnames(rlog_mat)
sample_info <- tibble(
  sample_id = sample_ids,
  condition = sub("^(M145_[123])_.*", "\\1", sample_ids)
)

# Align BGC scores to rlog sample order
bgc_scores_aligned <- bgc_scores %>%
  filter(sample_id %in% sample_ids) %>%
  arrange(match(sample_id, sample_ids))
stopifnot(all(bgc_scores_aligned$sample_id == sample_ids))

cat("\nBGC scores aligned to rlog sample order.\n")

# =============================================================================
# 2. Extract regulator genes
# =============================================================================
cat("\n--- Extracting regulator genes ---\n")

regulators <- master %>%
  filter(is_regulator == TRUE | is_regulator == "True" | is_regulator == "TRUE") %>%
  select(gene_id, old_locus_tag, gene_name, product, regulator_name, regulator_type)

cat("Regulators in master table:", nrow(regulators), "\n")

# Subset to those present in rlog matrix
reg_in_rlog <- regulators %>% filter(gene_id %in% rownames(rlog_mat))
cat("Regulators present in rlog:", nrow(reg_in_rlog), "\n")

# Extract rlog submatrix
rlog_reg <- rlog_mat[reg_in_rlog$gene_id, , drop = FALSE]
cat("Regulator rlog matrix:", nrow(rlog_reg), "x", ncol(rlog_reg), "\n")

# =============================================================================
# 3. Known regulator list
# =============================================================================
known_regs <- tribble(
  ~gene_id,      ~known_name,      ~known_target,
  "SC_RS27570",  "actII-orf4",     "act",
  "SC_RS27585",  "SCO5085_SARP",   "act",
  "SC_RS31630",  "redD",           "red",
  "SC_RS31650",  "redZ",           "red",
  "SC_RS18200",  "SCO3217_SARP",   "cda",
  "SC_RS18240",  "absA1",          "cda",
  "SC_RS18245",  "absA2",          "cda",
  "SC_RS33660",  "cpkO/kasO",      "cpk",
  "SC_RS33650",  "SCO6280_SARP",   "cpk",
  "SC_RS33680",  "scbR2",          "cpk",
  "SC_RS33690",  "SCO6288_SARP",   "cpk",
  "SC_RS33575",  "scbR",           "cpk/pleiotropic",
  "SC_RS25560",  "atrA",           "act/pleiotropic"
)

# =============================================================================
# 4. Pearson correlation: regulator rlog vs BGC scores
# =============================================================================
cat("\n--- Computing Pearson correlations ---\n")

act_scores <- bgc_scores_aligned$act_score
red_scores <- bgc_scores_aligned$red_score
cda_scores <- bgc_scores_aligned$cda_score
cpk_scores <- bgc_scores_aligned$cpk_score

corr_results <- tibble(
  gene_id = rownames(rlog_reg),
  corr_act = apply(rlog_reg, 1, function(x) cor(x, act_scores, use = "complete.obs")),
  corr_red = apply(rlog_reg, 1, function(x) cor(x, red_scores, use = "complete.obs")),
  corr_cda = apply(rlog_reg, 1, function(x) cor(x, cda_scores, use = "complete.obs")),
  corr_cpk = apply(rlog_reg, 1, function(x) cor(x, cpk_scores, use = "complete.obs"))
) %>%
  mutate(
    abs_corr_act = abs(corr_act),
    abs_corr_red = abs(corr_red),
    abs_corr_cda = abs(corr_cda),
    abs_corr_cpk = abs(corr_cpk)
  )

# Merge annotation
corr_full <- corr_results %>%
  left_join(reg_in_rlog, by = "gene_id") %>%
  left_join(known_regs, by = "gene_id") %>%
  mutate(is_known = !is.na(known_name)) %>%
  select(gene_id, old_locus_tag, gene_name, product,
         regulator_name, regulator_type, known_name, is_known,
         corr_act, corr_red, corr_cda, corr_cpk,
         abs_corr_act, abs_corr_red, abs_corr_cda, abs_corr_cpk)

write_tsv(corr_full, file.path(TABLE_DIR, "regulator_BGC_correlation.tsv"))
cat("Saved: regulator_BGC_correlation.tsv (", nrow(corr_full), " regulators)\n")

# =============================================================================
# 5. Top hits ranking per BGC
# =============================================================================
cat("\n--- Ranking top hits per BGC ---\n")

top_n <- 30

make_ranking <- function(df, bgc_name, corr_col, abs_col) {
  df %>%
    arrange(desc(.data[[abs_col]])) %>%
    head(top_n) %>%
    mutate(
      bgc_name = bgc_name,
      corr = .data[[corr_col]],
      sign = ifelse(corr >= 0, "+", "-"),
      corr_rank = row_number()
    ) %>%
    select(bgc_name, gene_id, old_locus_tag, gene_name, product,
           known_name, is_known, corr, sign, corr_rank)
}

top_act <- make_ranking(corr_full, "act", "corr_act", "abs_corr_act")
top_red <- make_ranking(corr_full, "red", "corr_red", "abs_corr_red")
top_cda <- make_ranking(corr_full, "cda", "corr_cda", "abs_corr_cda")
top_cpk <- make_ranking(corr_full, "cpk", "corr_cpk", "abs_corr_cpk")

top_all <- bind_rows(top_act, top_red, top_cda, top_cpk)
write_tsv(top_all, file.path(TABLE_DIR, "regulator_top_hits_by_BGC.tsv"))
cat("Saved: regulator_top_hits_by_BGC.tsv (", nrow(top_all), " entries)\n")

# Print known regulators in top hits
cat("\n--- Known regulators in top 30 ---\n")
top_all %>%
  filter(is_known) %>%
  select(bgc_name, corr_rank, gene_id, known_name, corr) %>%
  as.data.frame() %>%
  print()

# =============================================================================
# 6. Known vs novel table
# =============================================================================
cat("\n--- Known vs novel summary ---\n")

known_corrs <- corr_full %>%
  filter(is_known) %>%
  select(gene_id, old_locus_tag, known_name,
         corr_act, corr_red, corr_cda, corr_cpk) %>%
  arrange(known_name)

cat("Known regulator correlations:\n")
print(as.data.frame(known_corrs))

write_tsv(known_corrs, file.path(TABLE_DIR, "regulator_BGC_known_vs_novel.tsv"))
cat("Saved: regulator_BGC_known_vs_novel.tsv\n")

# =============================================================================
# 7. Correlation heatmap of known regulators
# =============================================================================
cat("\n--- Creating known regulator heatmap ---\n")

hm_data <- known_corrs %>%
  column_to_rownames("known_name") %>%
  select(corr_act, corr_red, corr_cda, corr_cpk) %>%
  as.matrix()
colnames(hm_data) <- c("act", "red", "cda", "cpk")

# Remove gene_id and old_locus_tag rows if they exist
hm_data <- hm_data[!rownames(hm_data) %in% c("gene_id", "old_locus_tag"), , drop = FALSE]

# Color: blue (negative) -> white -> red (positive)
hm_colors <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)

pdf(file.path(FIG_DIR, "heatmap_known_regulators_BGC_corr_M145.pdf"),
    width = 7, height = 6)
pheatmap(hm_data,
         cluster_rows = TRUE, cluster_cols = FALSE,
         display_numbers = TRUE, number_format = "%.2f",
         main = "Known Regulators: Pearson Correlation with BGC Scores",
         color = hm_colors,
         breaks = seq(-1, 1, length.out = 101),
         fontsize = 11, fontsize_number = 9,
         angle_col = 0)
dev.off()

svglite(file.path(FIG_DIR, "heatmap_known_regulators_BGC_corr_M145.svg"),
        width = 7, height = 6)
pheatmap(hm_data,
         cluster_rows = TRUE, cluster_cols = FALSE,
         display_numbers = TRUE, number_format = "%.2f",
         main = "Known Regulators: Pearson Correlation with BGC Scores",
         color = hm_colors,
         breaks = seq(-1, 1, length.out = 101),
         fontsize = 11, fontsize_number = 9,
         angle_col = 0)
dev.off()
cat("Saved: heatmap_known_regulators_BGC_corr_M145.pdf/.svg\n")

# =============================================================================
# 8. Top-30 correlation heatmap per BGC
# =============================================================================
cat("\n--- Creating top-30 heatmap ---\n")

# Gather unique gene_ids from all top-30 lists
unique_top_ids <- unique(top_all$gene_id)

top_hm <- corr_full %>%
  filter(gene_id %in% unique_top_ids) %>%
  mutate(label = ifelse(!is.na(known_name), known_name,
                        ifelse(!is.na(gene_name) & gene_name != "NA",
                               paste0(old_locus_tag, " (", gene_name, ")"),
                               old_locus_tag))) %>%
  select(label, corr_act, corr_red, corr_cda, corr_cpk)

# Deduplicate labels
dup_labels <- duplicated(top_hm$label)
if (any(dup_labels)) {
  top_hm$label[dup_labels] <- paste0(top_hm$label[dup_labels], "_2")
}

top_mat <- top_hm %>%
  column_to_rownames("label") %>%
  as.matrix()
colnames(top_mat) <- c("act", "red", "cda", "cpk")

fig_h <- max(8, nrow(top_mat) * 0.22 + 2)

pdf(file.path(FIG_DIR, "heatmap_top_regulators_all_BGC_M145.pdf"),
    width = 8, height = fig_h)
pheatmap(top_mat,
         cluster_rows = TRUE, cluster_cols = FALSE,
         display_numbers = FALSE,
         main = "Top Regulators: Pearson Correlation with 4 BGC Scores",
         color = hm_colors,
         breaks = seq(-1, 1, length.out = 101),
         fontsize = 8, fontsize_row = 6,
         angle_col = 0)
dev.off()

svglite(file.path(FIG_DIR, "heatmap_top_regulators_all_BGC_M145.svg"),
        width = 8, height = fig_h)
pheatmap(top_mat,
         cluster_rows = TRUE, cluster_cols = FALSE,
         display_numbers = FALSE,
         main = "Top Regulators: Pearson Correlation with 4 BGC Scores",
         color = hm_colors,
         breaks = seq(-1, 1, length.out = 101),
         fontsize = 8, fontsize_row = 6,
         angle_col = 0)
dev.off()
cat("Saved: heatmap_top_regulators_all_BGC_M145.pdf/.svg\n")

# =============================================================================
# 9. Scatter plots: BGC score vs representative regulator rlog
# =============================================================================
cat("\n--- Creating scatter plots ---\n")

# Build sample-level data frame for plotting
scatter_df <- tibble(
  sample_id = sample_ids,
  condition = sample_info$condition,
  act_score = act_scores,
  red_score = red_scores,
  cda_score = cda_scores,
  cpk_score = cpk_scores
)

# Add rlog values for key regulators
key_reg_ids <- c(
  "SC_RS27570",  # actII-orf4
  "SC_RS31630",  # redD
  "SC_RS31650",  # redZ
  "SC_RS18240",  # absA1
  "SC_RS18245",  # absA2
  "SC_RS33660",  # cpkO/kasO
  "SC_RS33680",  # scbR2
  "SC_RS33575",  # scbR
  "SC_RS25560"   # atrA
)

for (gid in key_reg_ids) {
  if (gid %in% rownames(rlog_mat)) {
    scatter_df[[gid]] <- rlog_mat[gid, sample_ids]
  }
}

cond_colors <- c(M145_1 = "#4393C3", M145_2 = "#FDDBC7", M145_3 = "#D6604D")

make_scatter <- function(df, x_col, y_col, x_lab, y_lab, title, fname) {
  p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], color = condition)) +
    geom_point(size = 3.5, alpha = 0.9) +
    scale_color_manual(values = cond_colors) +
    labs(title = title, x = x_lab, y = y_lab) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold", size = 12))

  # Add correlation annotation
  r_val <- cor(df[[x_col]], df[[y_col]], use = "complete.obs")
  p <- p + annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
                     label = sprintf("r = %.3f", r_val), size = 4, fontface = "italic")

  ggsave(file.path(FIG_DIR, paste0(fname, ".pdf")), p, width = 6, height = 5)
  ggsave(file.path(FIG_DIR, paste0(fname, ".svg")), p, width = 6, height = 5)
  cat("Saved:", fname, ".pdf/.svg\n")
}

# act
make_scatter(scatter_df, "act_score", "SC_RS27570",
             "act BGC score", "actII-orf4 rlog expression",
             "act score vs actII-orf4 (SC_RS27570)",
             "scatter_act_score_vs_actII_orf4_M145")

# red
make_scatter(scatter_df, "red_score", "SC_RS31630",
             "red BGC score", "redD rlog expression",
             "red score vs redD (SC_RS31630)",
             "scatter_red_score_vs_redD_M145")

make_scatter(scatter_df, "red_score", "SC_RS31650",
             "red BGC score", "redZ rlog expression",
             "red score vs redZ (SC_RS31650)",
             "scatter_red_score_vs_redZ_M145")

# cda
make_scatter(scatter_df, "cda_score", "SC_RS18240",
             "cda BGC score", "absA1 rlog expression",
             "cda score vs absA1 (SC_RS18240)",
             "scatter_cda_score_vs_absA1_M145")

# cpk
make_scatter(scatter_df, "cpk_score", "SC_RS33660",
             "cpk BGC score", "cpkO/kasO rlog expression",
             "cpk score vs cpkO/kasO (SC_RS33660)",
             "scatter_cpk_score_vs_cpkO_M145")

make_scatter(scatter_df, "cpk_score", "SC_RS33680",
             "cpk BGC score", "scbR2 rlog expression",
             "cpk score vs scbR2 (SC_RS33680)",
             "scatter_cpk_score_vs_scbR2_M145")

# atrA vs act
make_scatter(scatter_df, "act_score", "SC_RS25560",
             "act BGC score", "atrA rlog expression",
             "act score vs atrA (SC_RS25560)",
             "scatter_act_score_vs_atrA_M145")

# scbR vs cpk
make_scatter(scatter_df, "cpk_score", "SC_RS33575",
             "cpk BGC score", "scbR rlog expression",
             "cpk score vs scbR (SC_RS33575)",
             "scatter_cpk_score_vs_scbR_M145")

# =============================================================================
# 10. Timecourse plots: BGC score + regulator rlog
# =============================================================================
cat("\n--- Creating timecourse plots ---\n")

# Compute condition means for rlog
rlog_cond_means <- scatter_df %>%
  group_by(condition) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(condition = factor(condition, levels = c("M145_1", "M145_2", "M145_3")))

make_timecourse <- function(cond_df, bgc_col, reg_cols, reg_labels,
                            bgc_label, title, fname) {
  # Normalize: z-score each variable for dual-axis display
  bgc_vals <- cond_df[[bgc_col]]
  bgc_z <- (bgc_vals - mean(bgc_vals)) / sd(bgc_vals)

  plot_df <- tibble(
    condition = cond_df$condition,
    value = bgc_z,
    variable = bgc_label
  )

  for (i in seq_along(reg_cols)) {
    rv <- cond_df[[reg_cols[i]]]
    rz <- (rv - mean(rv)) / sd(rv)
    plot_df <- bind_rows(plot_df, tibble(
      condition = cond_df$condition,
      value = rz,
      variable = reg_labels[i]
    ))
  }

  plot_df$variable <- factor(plot_df$variable,
                              levels = c(bgc_label, reg_labels))

  n_vars <- length(unique(plot_df$variable))
  line_colors <- c("#D6604D", "#2166AC", "#4DAF4A", "#FF7F00")[1:n_vars]
  names(line_colors) <- levels(plot_df$variable)

  p <- ggplot(plot_df, aes(x = condition, y = value,
                            color = variable, group = variable)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    scale_color_manual(values = line_colors, name = "") +
    labs(title = title,
         subtitle = "Z-scored condition means (rlog for regulators, normalized count for BGC score)",
         x = "Condition (timepoint)", y = "Z-score") +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold", size = 12),
          legend.position = "bottom")

  ggsave(file.path(FIG_DIR, paste0(fname, ".pdf")), p, width = 7, height = 5)
  ggsave(file.path(FIG_DIR, paste0(fname, ".svg")), p, width = 7, height = 5)
  cat("Saved:", fname, ".pdf/.svg\n")
}

# act + actII-orf4
make_timecourse(rlog_cond_means, "act_score",
                c("SC_RS27570"), c("actII-orf4"),
                "act BGC score", "act: BGC Score + actII-orf4 Timecourse",
                "timecourse_act_score_with_actII_orf4_M145")

# red + redD, redZ
make_timecourse(rlog_cond_means, "red_score",
                c("SC_RS31630", "SC_RS31650"), c("redD", "redZ"),
                "red BGC score", "red: BGC Score + redD/redZ Timecourse",
                "timecourse_red_score_with_redD_redZ_M145")

# cda + absA1, absA2
make_timecourse(rlog_cond_means, "cda_score",
                c("SC_RS18240", "SC_RS18245"), c("absA1", "absA2"),
                "cda BGC score", "cda: BGC Score + absA1/absA2 Timecourse",
                "timecourse_cda_score_with_absA1_absA2_M145")

# cpk + cpkO, scbR2
make_timecourse(rlog_cond_means, "cpk_score",
                c("SC_RS33660", "SC_RS33680"), c("cpkO/kasO", "scbR2"),
                "cpk BGC score", "cpk: BGC Score + cpkO/scbR2 Timecourse",
                "timecourse_cpk_score_with_cpkO_scbR2_M145")

# =============================================================================
# 11. Global regulator correlation profile
# =============================================================================
cat("\n--- Global regulator analysis ---\n")

# Find regulators highly correlated with multiple BGCs
corr_full <- corr_full %>%
  mutate(
    n_high_corr = (abs_corr_act > 0.9) + (abs_corr_red > 0.9) +
                  (abs_corr_cda > 0.9) + (abs_corr_cpk > 0.9),
    mean_abs_corr = (abs_corr_act + abs_corr_red + abs_corr_cda + abs_corr_cpk) / 4
  )

global_candidates <- corr_full %>%
  filter(n_high_corr >= 2) %>%
  arrange(desc(mean_abs_corr))

cat("Regulators with |corr| > 0.9 for >= 2 BGCs:", nrow(global_candidates), "\n")
if (nrow(global_candidates) > 0) {
  global_candidates %>%
    select(gene_id, old_locus_tag, gene_name, product, known_name,
           corr_act, corr_red, corr_cda, corr_cpk, n_high_corr) %>%
    head(20) %>%
    as.data.frame() %>%
    print()
}

# Save updated correlation table with global info
write_tsv(corr_full, file.path(TABLE_DIR, "regulator_BGC_correlation.tsv"))

# =============================================================================
# 12. Summary statistics
# =============================================================================
cat("\n--- Summary statistics ---\n")

for (bgc in c("act", "red", "cda", "cpk")) {
  col <- paste0("corr_", bgc)
  acol <- paste0("abs_corr_", bgc)
  cat(sprintf("%s: median |corr| = %.3f, n(|corr|>0.9) = %d, n(|corr|>0.8) = %d\n",
              bgc,
              median(corr_full[[acol]], na.rm = TRUE),
              sum(corr_full[[acol]] > 0.9, na.rm = TRUE),
              sum(corr_full[[acol]] > 0.8, na.rm = TRUE)))
}

# =============================================================================
# Done
# =============================================================================
cat("\n=== 07_regulator_network completed at", format(Sys.time()), "===\n")
cat("Output directory:", REG_RUN_DIR, "\n")
