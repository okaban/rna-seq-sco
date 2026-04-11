#!/usr/bin/env Rscript
# =============================================================================
# 06_BGC_dynamics: BGC expression dynamics analysis for M145 RNA-seq
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(svglite)
})

cat("=== 06_BGC_dynamics started at", format(Sys.time()), "===\n")

# --- Paths -------------------------------------------------------------------
BGC_RUN_DIR  <- "/Users/okaban/bioinfo/rna-seq/06_BGC_dynamics/analysis/06_BGC_dynamics_260128_v1"
DESEQ_RUN_DIR <- "/Users/okaban/bioinfo/rna-seq/04_deseq2/analysis/04_deseq2_260128_v1"
ANNOT_RUN_DIR <- "/Users/okaban/bioinfo/rna-seq/05_annotation/analysis/05_annotation_260128_v1"

NORM_FILE   <- file.path(DESEQ_RUN_DIR, "results", "normalized_counts_M145.tsv")
MASTER_FILE <- file.path(ANNOT_RUN_DIR, "tables", "gene_master_with_BGC_regulators.tsv")
BGC_DEF_FILE <- file.path(ANNOT_RUN_DIR, "tables", "BGC_definition_manual.tsv")

TABLE_DIR  <- file.path(BGC_RUN_DIR, "tables")
FIG_DIR    <- file.path(BGC_RUN_DIR, "figures")

cat("BGC_RUN_DIR:", BGC_RUN_DIR, "\n")
cat("Input files:\n")
cat("  NORM_FILE:", NORM_FILE, "\n")
cat("  MASTER_FILE:", MASTER_FILE, "\n")
cat("  BGC_DEF_FILE:", BGC_DEF_FILE, "\n")

# =============================================================================
# 1. Load data
# =============================================================================
cat("\n--- Loading data ---\n")

norm_counts <- read_tsv(NORM_FILE, show_col_types = FALSE)
cat("Normalized counts:", nrow(norm_counts), "genes x", ncol(norm_counts)-1, "samples\n")

master <- read_tsv(MASTER_FILE, show_col_types = FALSE)
cat("Master table:", nrow(master), "genes x", ncol(master), "columns\n")

bgc_def <- read_tsv(BGC_DEF_FILE, show_col_types = FALSE)
cat("BGC definition:", nrow(bgc_def), "genes in", length(unique(bgc_def$bgc_name)), "BGCs\n")

# =============================================================================
# 2. Sample -> condition mapping
# =============================================================================
sample_cols <- colnames(norm_counts)[-1]  # exclude gene_id
sample_info <- tibble(
  sample_id = sample_cols,
  condition = sub("^(M145_[123])_.*", "\\1", sample_cols)
)
cat("\nSample mapping:\n")
print(as.data.frame(sample_info))

# =============================================================================
# 3. Gene x condition means
# =============================================================================
cat("\n--- Computing gene x condition means ---\n")

# Pivot to long, join condition, compute means
norm_long <- norm_counts %>%
  pivot_longer(-gene_id, names_to = "sample_id", values_to = "norm_count") %>%
  left_join(sample_info, by = "sample_id")

gene_cond_means <- norm_long %>%
  group_by(gene_id, condition) %>%
  summarise(mean_count = mean(norm_count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = mean_count,
              names_glue = "{condition}_mean")

cat("Gene-condition means:", nrow(gene_cond_means), "genes\n")

write_tsv(gene_cond_means, file.path(TABLE_DIR, "gene_condition_means.tsv"))
cat("Saved: gene_condition_means.tsv\n")

# =============================================================================
# 4. BGC-level condition means
# =============================================================================
cat("\n--- Computing BGC-level condition means ---\n")

# Merge BGC definition with gene-condition means
bgc_genes <- bgc_def %>%
  select(bgc_name, gene_id) %>%
  left_join(gene_cond_means, by = "gene_id")

# Check for missing genes
n_missing <- sum(is.na(bgc_genes$M145_1_mean))
if (n_missing > 0) {
  cat("WARNING:", n_missing, "BGC genes not found in normalized counts\n")
  missing_genes <- bgc_genes %>% filter(is.na(M145_1_mean))
  print(as.data.frame(missing_genes))
}

bgc_cond_means <- bgc_genes %>%
  filter(!is.na(M145_1_mean)) %>%
  group_by(bgc_name) %>%
  summarise(
    n_genes = n(),
    M145_1_mean = mean(M145_1_mean, na.rm = TRUE),
    M145_2_mean = mean(M145_2_mean, na.rm = TRUE),
    M145_3_mean = mean(M145_3_mean, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nBGC condition means:\n")
print(as.data.frame(bgc_cond_means))

write_tsv(bgc_cond_means, file.path(TABLE_DIR, "BGC_condition_means.tsv"))
cat("Saved: BGC_condition_means.tsv\n")

# =============================================================================
# 5. Figures
# =============================================================================

# --- 5.1 BGC timecourse lineplot ---
cat("\n--- Creating BGC timecourse lineplot ---\n")

bgc_long <- bgc_cond_means %>%
  select(bgc_name, M145_1_mean, M145_2_mean, M145_3_mean) %>%
  pivot_longer(-bgc_name, names_to = "condition", values_to = "mean_expr") %>%
  mutate(
    condition = sub("_mean$", "", condition),
    condition = factor(condition, levels = c("M145_1", "M145_2", "M145_3"))
  )

# Also compute log10(mean+1) version
bgc_long <- bgc_long %>%
  mutate(log10_expr = log10(mean_expr + 1))

# Color palette for BGCs
bgc_colors <- c(act = "#2166AC", red = "#B2182B", cda = "#4DAF4A", cpk = "#FF7F00")

# Raw scale lineplot
p_line_raw <- ggplot(bgc_long, aes(x = condition, y = mean_expr,
                                    color = bgc_name, group = bgc_name)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = bgc_colors, name = "BGC") +
  labs(title = "Major 4 BGC Expression Timecourse (M145)",
       subtitle = "Mean normalized counts across cluster genes",
       x = "Condition (timepoint)",
       y = "Mean normalized count") +
  theme_bw(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"))

ggsave(file.path(FIG_DIR, "BGC_timecourse_lineplot_M145.pdf"),
       p_line_raw, width = 8, height = 5)
ggsave(file.path(FIG_DIR, "BGC_timecourse_lineplot_M145.svg"),
       p_line_raw, width = 8, height = 5)
cat("Saved: BGC_timecourse_lineplot_M145.pdf/.svg\n")

# Log10 scale lineplot
p_line_log <- ggplot(bgc_long, aes(x = condition, y = log10_expr,
                                    color = bgc_name, group = bgc_name)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = bgc_colors, name = "BGC") +
  labs(title = "Major 4 BGC Expression Timecourse (M145, log10)",
       subtitle = "log10(mean normalized count + 1)",
       x = "Condition (timepoint)",
       y = "log10(mean + 1)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"))

ggsave(file.path(FIG_DIR, "BGC_timecourse_lineplot_log10_M145.pdf"),
       p_line_log, width = 8, height = 5)
ggsave(file.path(FIG_DIR, "BGC_timecourse_lineplot_log10_M145.svg"),
       p_line_log, width = 8, height = 5)
cat("Saved: BGC_timecourse_lineplot_log10_M145.pdf/.svg\n")

# --- 5.2 BGC x condition heatmap ---
cat("\n--- Creating BGC x condition heatmap ---\n")

# Prepare matrix (log10 scale for better visibility)
hm_mat <- bgc_cond_means %>%
  select(bgc_name, M145_1_mean, M145_2_mean, M145_3_mean) %>%
  column_to_rownames("bgc_name") %>%
  as.matrix()
colnames(hm_mat) <- c("M145_1", "M145_2", "M145_3")

# Log10 transform
hm_mat_log <- log10(hm_mat + 1)

# Z-score (row-wise)
hm_mat_z <- t(scale(t(hm_mat_log)))

# Heatmap (log10 scale)
pdf(file.path(FIG_DIR, "BGC_condition_heatmap_M145.pdf"), width = 6, height = 4)
pheatmap(hm_mat_log,
         cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE, number_format = "%.2f",
         main = "BGC Mean Expression (log10) by Condition",
         color = colorRampPalette(c("#F7F7F7", "#FEE08B", "#D73027"))(100),
         fontsize = 12, fontsize_number = 10,
         angle_col = 0)
dev.off()

svglite(file.path(FIG_DIR, "BGC_condition_heatmap_M145.svg"), width = 6, height = 4)
pheatmap(hm_mat_log,
         cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE, number_format = "%.2f",
         main = "BGC Mean Expression (log10) by Condition",
         color = colorRampPalette(c("#F7F7F7", "#FEE08B", "#D73027"))(100),
         fontsize = 12, fontsize_number = 10,
         angle_col = 0)
dev.off()
cat("Saved: BGC_condition_heatmap_M145.pdf/.svg\n")

# --- 5.3 Per-gene heatmap within each BGC ---
cat("\n--- Creating per-BGC gene-level heatmaps ---\n")

# Get per-gene condition means for BGC genes
bgc_gene_means <- bgc_def %>%
  select(bgc_name, gene_id, old_locus_tag, gene_name) %>%
  left_join(gene_cond_means, by = "gene_id") %>%
  filter(!is.na(M145_1_mean))

for (bgc in c("act", "red", "cda", "cpk")) {
  sub_df <- bgc_gene_means %>% filter(bgc_name == bgc)

  # Create row labels: old_locus_tag (gene_name if available)
  row_labels <- ifelse(is.na(sub_df$gene_name) | sub_df$gene_name == "NA",
                        sub_df$old_locus_tag,
                        paste0(sub_df$old_locus_tag, " (", sub_df$gene_name, ")"))

  mat <- sub_df %>%
    select(M145_1_mean, M145_2_mean, M145_3_mean) %>%
    as.matrix()
  rownames(mat) <- row_labels
  colnames(mat) <- c("M145_1", "M145_2", "M145_3")

  # Log10 transform
  mat_log <- log10(mat + 1)

  # Figure height scales with number of genes
  fig_h <- max(4, nrow(mat) * 0.35 + 1.5)

  pdf(file.path(FIG_DIR, paste0("BGC_gene_heatmap_", bgc, "_M145.pdf")),
      width = 7, height = fig_h)
  pheatmap(mat_log,
           cluster_rows = FALSE, cluster_cols = FALSE,
           display_numbers = TRUE, number_format = "%.1f",
           main = paste0(toupper(bgc), " cluster gene expression (log10)"),
           color = colorRampPalette(c("#F7F7F7", "#FEE08B", "#D73027"))(100),
           fontsize = 10, fontsize_number = 8,
           angle_col = 0)
  dev.off()

  svglite(file.path(FIG_DIR, paste0("BGC_gene_heatmap_", bgc, "_M145.svg")),
          width = 7, height = fig_h)
  pheatmap(mat_log,
           cluster_rows = FALSE, cluster_cols = FALSE,
           display_numbers = TRUE, number_format = "%.1f",
           main = paste0(toupper(bgc), " cluster gene expression (log10)"),
           color = colorRampPalette(c("#F7F7F7", "#FEE08B", "#D73027"))(100),
           fontsize = 10, fontsize_number = 8,
           angle_col = 0)
  dev.off()

  cat("Saved: BGC_gene_heatmap_", bgc, "_M145.pdf/.svg\n", sep = "")
}

# =============================================================================
# 6. Key genes identification
# =============================================================================
cat("\n--- Identifying key genes per BGC ---\n")

# Get DESeq2 info for BGC genes from master table
master_bgc <- master %>%
  filter(!is.na(bgc_name) & bgc_name != "") %>%
  filter(bgc_name %in% c("act", "red", "cda", "cpk")) %>%
  select(gene_id, old_locus_tag, gene_name, product, bgc_name, bgc_role,
         is_regulator, regulator_name, regulator_type,
         log2FoldChange_3_vs_1, padj_3_vs_1,
         log2FoldChange_2_vs_1, padj_2_vs_1,
         log2FoldChange_3_vs_2, padj_3_vs_2)

# Join with BGC definition role
bgc_def_role <- bgc_def %>% select(gene_id, role)
master_bgc <- master_bgc %>%
  left_join(bgc_def_role, by = "gene_id") %>%
  mutate(role = coalesce(role, bgc_role))

# For each BGC: top 5 by |LFC_3_vs_1| + all regulators
key_genes_list <- list()

for (bgc in c("act", "red", "cda", "cpk")) {
  sub <- master_bgc %>% filter(bgc_name == bgc)

  # Top 5 by LFC
  top5 <- sub %>%
    filter(!is.na(log2FoldChange_3_vs_1)) %>%
    arrange(desc(abs(log2FoldChange_3_vs_1))) %>%
    head(5) %>%
    mutate(key_reason = "top_LFC")

  # All regulators in this BGC
  regs <- sub %>%
    filter(is_regulator == TRUE | role == "regulator") %>%
    mutate(key_reason = "regulator")

  # Combine, deduplicate
  combined <- bind_rows(top5, regs) %>%
    distinct(gene_id, .keep_all = TRUE)

  key_genes_list[[bgc]] <- combined
}

key_genes <- bind_rows(key_genes_list) %>%
  select(bgc_name, gene_id, old_locus_tag, gene_name, product, role,
         is_regulator, regulator_name,
         log2FoldChange_3_vs_1, padj_3_vs_1,
         log2FoldChange_2_vs_1, log2FoldChange_3_vs_2,
         key_reason) %>%
  arrange(bgc_name, desc(abs(log2FoldChange_3_vs_1)))

cat("Key genes summary:\n")
cat("  act:", sum(key_genes$bgc_name == "act"), "genes\n")
cat("  red:", sum(key_genes$bgc_name == "red"), "genes\n")
cat("  cda:", sum(key_genes$bgc_name == "cda"), "genes\n")
cat("  cpk:", sum(key_genes$bgc_name == "cpk"), "genes\n")
cat("  Total:", nrow(key_genes), "key genes\n")

write_tsv(key_genes, file.path(TABLE_DIR, "BGC_key_genes_summary.tsv"))
cat("Saved: BGC_key_genes_summary.tsv\n")

# Print key genes table
cat("\nKey genes detail:\n")
key_genes %>%
  select(bgc_name, gene_id, old_locus_tag, gene_name,
         log2FoldChange_3_vs_1, padj_3_vs_1, key_reason) %>%
  as.data.frame() %>%
  print()

# =============================================================================
# 7. Sample-level BGC scores
# =============================================================================
cat("\n--- Computing sample-level BGC scores ---\n")

bgc_sample_scores <- tibble(sample_id = sample_cols)

for (bgc in c("act", "red", "cda", "cpk")) {
  bgc_gene_ids <- bgc_def %>% filter(bgc_name == bgc) %>% pull(gene_id)

  # Filter normalized counts to BGC genes
  bgc_counts <- norm_counts %>%
    filter(gene_id %in% bgc_gene_ids) %>%
    select(-gene_id)

  # Mean across genes for each sample
  score <- colMeans(bgc_counts, na.rm = TRUE)

  bgc_sample_scores <- bgc_sample_scores %>%
    mutate(!!paste0(bgc, "_score") := score[sample_id])
}

bgc_sample_scores <- bgc_sample_scores %>%
  left_join(sample_info, by = "sample_id") %>%
  select(sample_id, condition, everything())

cat("BGC sample scores:\n")
print(as.data.frame(bgc_sample_scores))

write_tsv(bgc_sample_scores, file.path(TABLE_DIR, "BGC_sample_scores.tsv"))
cat("Saved: BGC_sample_scores.tsv\n")

# --- 7.1 BGC score boxplot ---
cat("\n--- Creating BGC score boxplot ---\n")

scores_long <- bgc_sample_scores %>%
  pivot_longer(cols = ends_with("_score"),
               names_to = "bgc", values_to = "score") %>%
  mutate(
    bgc = sub("_score$", "", bgc),
    condition = factor(condition, levels = c("M145_1", "M145_2", "M145_3"))
  )

p_box <- ggplot(scores_long, aes(x = condition, y = score, fill = condition)) +
  geom_boxplot(alpha = 0.7, outlier.size = 2) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.8) +
  facet_wrap(~ bgc, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c(M145_1 = "#4393C3", M145_2 = "#FDDBC7", M145_3 = "#D6604D")) +
  labs(title = "Per-sample BGC Expression Scores (M145)",
       subtitle = "Mean normalized count of cluster genes per sample",
       x = "Condition", y = "BGC score (mean norm count)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", size = 12))

ggsave(file.path(FIG_DIR, "BGC_sample_scores_boxplot_M145.pdf"),
       p_box, width = 8, height = 6)
ggsave(file.path(FIG_DIR, "BGC_sample_scores_boxplot_M145.svg"),
       p_box, width = 8, height = 6)
cat("Saved: BGC_sample_scores_boxplot_M145.pdf/.svg\n")

# =============================================================================
# 8. Fold-change summary statistics
# =============================================================================
cat("\n--- BGC fold-change summary ---\n")

fc_summary <- master_bgc %>%
  group_by(bgc_name) %>%
  summarise(
    n_total = n(),
    n_sig_3v1 = sum(padj_3_vs_1 < 0.05, na.rm = TRUE),
    n_up_3v1 = sum(log2FoldChange_3_vs_1 > 1 & padj_3_vs_1 < 0.05, na.rm = TRUE),
    median_LFC_3v1 = median(log2FoldChange_3_vs_1, na.rm = TRUE),
    max_LFC_3v1 = max(log2FoldChange_3_vs_1, na.rm = TRUE),
    n_sig_2v1 = sum(padj_2_vs_1 < 0.05, na.rm = TRUE),
    n_up_2v1 = sum(log2FoldChange_2_vs_1 > 1 & padj_2_vs_1 < 0.05, na.rm = TRUE),
    median_LFC_2v1 = median(log2FoldChange_2_vs_1, na.rm = TRUE),
    n_sig_3v2 = sum(padj_3_vs_2 < 0.05, na.rm = TRUE),
    n_up_3v2 = sum(log2FoldChange_3_vs_2 > 1 & padj_3_vs_2 < 0.05, na.rm = TRUE),
    median_LFC_3v2 = median(log2FoldChange_3_vs_2, na.rm = TRUE),
    .groups = "drop"
  )

cat("Fold-change summary:\n")
print(as.data.frame(fc_summary))

# =============================================================================
# 9. Compute fold-change between conditions (for reporting)
# =============================================================================
cat("\n--- BGC expression fold-changes between conditions ---\n")

bgc_fc <- bgc_cond_means %>%
  mutate(
    fc_2v1 = M145_2_mean / M145_1_mean,
    fc_3v1 = M145_3_mean / M145_1_mean,
    fc_3v2 = M145_3_mean / M145_2_mean
  )
cat("BGC-level fold changes (cluster mean):\n")
print(as.data.frame(bgc_fc))

# =============================================================================
# Done
# =============================================================================
cat("\n=== 06_BGC_dynamics completed at", format(Sys.time()), "===\n")
cat("Output directory:", BGC_RUN_DIR, "\n")
