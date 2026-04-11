#!/usr/bin/env Rscript
# =============================================================================
# 03_COG_classification.R
# COG様の機能分類解析
# Product nameとKEGGパスウェイから機能カテゴリを推定
# Project: M145 RNA-seq
# Date: 2026-02-02
# Updated: 2026-02-05 - Added T3vsT2, SVG output, combined net change plot
# =============================================================================

library(tidyverse)
library(ggplot2)

# --- パス設定 ---
deseq2_dir <- "/Users/okaban/bioinfo/rna-seq/04_deseq2/analysis/04_deseq2_260128_v1/results"
anno_file <- "/Users/okaban/bioinfo/rna-seq/05_annotation/analysis/05_annotation_260128_v1/tables/gene_annotation_basic.tsv"
output_dir <- "/Users/okaban/bioinfo/rna-seq/12_supplementary_figures/analysis/12_supplementary_260202_v1"

# --- データ読み込み ---
de_3v1 <- read_tsv(file.path(deseq2_dir, "DESeq2_M145_3_vs_1.tsv"), show_col_types = FALSE)
de_2v1 <- read_tsv(file.path(deseq2_dir, "DESeq2_M145_2_vs_1.tsv"), show_col_types = FALSE)
de_3v2 <- read_tsv(file.path(deseq2_dir, "DESeq2_M145_3_vs_2.tsv"), show_col_types = FALSE)
anno <- read_tsv(anno_file, show_col_types = FALSE)

# --- COG様カテゴリの定義（キーワードベース）---
cog_categories <- list(
  "J - Translation" = c("ribosom", "tRNA", "aminoacyl", "translation", "elongation factor"),
  "K - Transcription" = c("transcription", "RNA polymerase", "sigma factor", "regulator", "repressor", "activator", "DNA-binding"),
  "L - Replication/Repair" = c("DNA polymerase", "helicase", "recombinase", "ligase", "topoisomerase", "repair", "gyrase"),
  "D - Cell division" = c("cell division", "FtsZ", "septum", "chromosome partitioning"),
  "V - Defense" = c("restriction", "toxin-antitoxin", "defense", "CRISPR", "methyltransferase"),
  "T - Signal transduction" = c("kinase", "phosphatase", "two-component", "sensor", "response regulator", "histidine kinase"),
  "M - Cell envelope" = c("peptidoglycan", "membrane", "cell wall", "lipopolysaccharide", "penicillin-binding"),
  "N - Cell motility" = c("flagell", "pili", "motility", "chemotaxis"),
  "U - Secretion" = c("secretion", "sec", "tat", "export", "type II", "type III", "type IV"),
  "O - Protein modification" = c("protease", "chaperone", "heat shock", "GroE", "DnaK", "peptidase", "ubiquitin"),
  "C - Energy production" = c("cytochrome", "NADH", "ATP synthase", "dehydrogenase", "oxidoreductase", "respiratory", "electron transfer"),
  "G - Carbohydrate metabolism" = c("glycosyl", "sugar", "glucos", "galactos", "xylos", "arabinos", "maltose", "cellulase"),
  "E - Amino acid metabolism" = c("amino acid", "aminotransferase", "synthase", "glutam", "aspart", "alanine", "serine", "threonine", "lysine", "arginine", "histidine", "proline", "tryptophan", "tyrosine", "phenylalanine", "cysteine", "methionine", "valine", "leucine", "isoleucine"),
  "F - Nucleotide metabolism" = c("nucleotide", "purine", "pyrimidine", "nucleoside", "GMP", "AMP", "CTP", "UMP"),
  "H - Coenzyme metabolism" = c("coenzyme", "CoA", "NAD", "FAD", "FMN", "biotin", "folate", "thiamin", "riboflavin", "pyridoxal", "cobalamin"),
  "I - Lipid metabolism" = c("lipid", "fatty acid", "acyl", "phospholipid", "sterol", "lipase"),
  "P - Inorganic ion transport" = c("iron", "zinc", "copper", "manganese", "magnesium", "sulfate", "phosphate", "ABC transporter", "permease"),
  "Q - Secondary metabolism" = c("polyketide", "NRPS", "terpene", "antibiotic", "pigment", "siderophore", "actinorhodin", "undecylprodigiosin", "calcium-dependent antibiotic", "coelibactin"),
  "S - Unknown" = c("hypothetical", "uncharacterized", "unknown", "DUF")
)

# --- 機能カテゴリの割り当て ---
assign_cog <- function(product) {
  if (is.na(product) || product == "") return("S - Unknown")

  product_lower <- tolower(product)

  for (category in names(cog_categories)) {
    keywords <- cog_categories[[category]]
    for (kw in keywords) {
      if (grepl(tolower(kw), product_lower, fixed = FALSE)) {
        return(category)
      }
    }
  }

  return("R - General function")
}

# アノテーションにCOGカテゴリを追加
anno_cog <- anno %>%
  mutate(COG_category = map_chr(product, assign_cog))

cat("COG category distribution:\n")
print(table(anno_cog$COG_category))

# --- DEGs の抽出と分類 ---
get_degs_cog <- function(de_df, anno_cog, padj_thresh = 0.05, lfc_thresh = 1) {
  de_df %>%
    filter(!is.na(padj), padj < padj_thresh, abs(log2FoldChange) >= lfc_thresh) %>%
    left_join(anno_cog %>% select(gene_id, product, COG_category), by = "gene_id") %>%
    mutate(direction = ifelse(log2FoldChange > 0, "Up", "Down"))
}

degs_3v1_cog <- get_degs_cog(de_3v1, anno_cog)
degs_2v1_cog <- get_degs_cog(de_2v1, anno_cog)
degs_3v2_cog <- get_degs_cog(de_3v2, anno_cog)

# --- COGカテゴリ別の集計 ---
summarize_cog <- function(degs_df, comparison_name) {
  degs_df %>%
    group_by(COG_category, direction) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(comparison = comparison_name)
}

cog_summary_3v1 <- summarize_cog(degs_3v1_cog, "T3 vs T1")
cog_summary_2v1 <- summarize_cog(degs_2v1_cog, "T2 vs T1")
cog_summary_3v2 <- summarize_cog(degs_3v2_cog, "T3 vs T2")

cog_summary_all <- bind_rows(cog_summary_3v1, cog_summary_2v1, cog_summary_3v2)

# --- バーチャート作成 ---
create_cog_barplot <- function(cog_df, title, output_prefix) {
  # Wideフォーマットに変換
  plot_df <- cog_df %>%
    pivot_wider(names_from = direction, values_from = count, values_fill = 0) %>%
    mutate(
      Up = ifelse(is.na(Up), 0, Up),
      Down = ifelse(is.na(Down), 0, -Down)  # Downを負の値に
    ) %>%
    pivot_longer(cols = c(Up, Down), names_to = "direction", values_to = "count") %>%
    mutate(COG_category = fct_reorder(COG_category, abs(count), .fun = sum))

  p <- ggplot(plot_df, aes(x = count, y = COG_category, fill = direction)) +
    geom_bar(stat = "identity", position = "identity", width = 0.7) +
    scale_fill_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8"),
                      name = "Regulation",
                      labels = c("Up" = "Up-regulated", "Down" = "Down-regulated")) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      plot.title = element_text(size = 12, face = "bold"),
      legend.position = "top"
    ) +
    labs(
      title = title,
      x = "Number of Genes",
      y = ""
    ) +
    scale_x_continuous(labels = abs)

  ggsave(paste0(output_prefix, ".pdf"), p, width = 12, height = 10)
  ggsave(paste0(output_prefix, ".png"), p, width = 12, height = 10, dpi = 150)
  ggsave(paste0(output_prefix, ".svg"), p, width = 12, height = 10)

  return(p)
}

# Figure出力
fig_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")

# T3 vs T1
create_cog_barplot(cog_summary_3v1, "Functional Categories of DEGs (T3 vs T1)",
                   file.path(fig_dir, "COG_barplot_T3vsT1"))
cat("Saved COG barplot T3 vs T1\n")

# T2 vs T1
create_cog_barplot(cog_summary_2v1, "Functional Categories of DEGs (T2 vs T1)",
                   file.path(fig_dir, "COG_barplot_T2vsT1"))
cat("Saved COG barplot T2 vs T1\n")

# T3 vs T2
create_cog_barplot(cog_summary_3v2, "Functional Categories of DEGs (T3 vs T2)",
                   file.path(fig_dir, "COG_barplot_T3vsT2"))
cat("Saved COG barplot T3 vs T2\n")

# --- 全比較の散布図プロット ---
compare_df <- cog_summary_all %>%
  group_by(COG_category, comparison) %>%
  summarise(
    Up = sum(count[direction == "Up"]),
    Down = sum(count[direction == "Down"]),
    .groups = "drop"
  )

p_compare <- ggplot(compare_df, aes(x = Up, y = Down, color = comparison)) +
  geom_point(size = 3) +
  geom_text(aes(label = COG_category), hjust = -0.1, vjust = 0.5, size = 2.5, check_overlap = TRUE) +
  theme_bw() +
  labs(
    title = "COG Category Distribution: Up vs Down Regulated",
    x = "Up-regulated genes",
    y = "Down-regulated genes",
    color = "Comparison"
  ) +
  theme(legend.position = "top")

ggsave(file.path(fig_dir, "COG_scatter_up_vs_down.pdf"), p_compare, width = 10, height = 10)
ggsave(file.path(fig_dir, "COG_scatter_up_vs_down.png"), p_compare, width = 10, height = 10, dpi = 150)
ggsave(file.path(fig_dir, "COG_scatter_up_vs_down.svg"), p_compare, width = 10, height = 10)

# --- Combined Net Change Plot (all 3 comparisons side by side) ---
net_change_df <- cog_summary_all %>%
  pivot_wider(names_from = direction, values_from = count, values_fill = 0) %>%
  mutate(
    Up = ifelse(is.na(Up), 0, Up),
    Down = ifelse(is.na(Down), 0, Down),
    net_change = Up - Down
  ) %>%
  mutate(
    comparison = factor(comparison, levels = c("T2 vs T1", "T3 vs T1", "T3 vs T2")),
    COG_category = fct_reorder(COG_category, abs(net_change), .fun = sum)
  )

p_net <- ggplot(net_change_df, aes(x = net_change, y = COG_category, fill = comparison)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(
    values = c("T2 vs T1" = "#4DAF4A", "T3 vs T1" = "#E41A1C", "T3 vs T2" = "#377EB8"),
    name = "Comparison"
  ) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 10),
    legend.position = "top"
  ) +
  labs(
    title = "Net Change in Functional Categories Across All Comparisons",
    subtitle = "Net change = (Up-regulated) - (Down-regulated)",
    x = "Net Change (number of genes)",
    y = ""
  )

ggsave(file.path(fig_dir, "COG_net_change_all_comparisons.pdf"), p_net, width = 14, height = 10)
ggsave(file.path(fig_dir, "COG_net_change_all_comparisons.png"), p_net, width = 14, height = 10, dpi = 150)
ggsave(file.path(fig_dir, "COG_net_change_all_comparisons.svg"), p_net, width = 14, height = 10)

cat("Saved combined net change plot (all 3 comparisons)\n")

# --- テーブル出力 ---
write_tsv(cog_summary_all, file.path(table_dir, "COG_classification_summary.tsv"))
write_tsv(anno_cog %>% select(gene_id, old_locus_tag, product, COG_category),
          file.path(table_dir, "gene_COG_classification.tsv"))

cat("\nCOG classification completed.\n")
cat("Categories assigned to", nrow(anno_cog), "genes\n")
