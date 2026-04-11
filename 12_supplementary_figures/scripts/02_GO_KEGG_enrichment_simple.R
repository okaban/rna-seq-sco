#!/usr/bin/env Rscript
# =============================================================================
# 02_GO_KEGG_enrichment_simple.R
# GO/KEGGエンリッチメント解析（シンプル実装）
# Fisher's exact testを使用
# Project: M145 RNA-seq
# Date: 2026-02-02
# =============================================================================

library(tidyverse)
library(ggplot2)

# --- パス設定 ---
deseq2_dir <- "/Users/okaban/bioinfo/rna-seq/04_deseq2/analysis/04_deseq2_260128_v1/results"
anno_file <- "/Users/okaban/bioinfo/rna-seq/05_annotation/analysis/05_annotation_260128_v1/tables/gene_annotation_basic.tsv"
output_dir <- "/Users/okaban/bioinfo/rna-seq/12_supplementary_figures/analysis/12_supplementary_260202_v1"

# --- データ読み込み ---
de_2v1 <- read_tsv(file.path(deseq2_dir, "DESeq2_M145_2_vs_1.tsv"), show_col_types = FALSE)
de_3v1 <- read_tsv(file.path(deseq2_dir, "DESeq2_M145_3_vs_1.tsv"), show_col_types = FALSE)
de_3v2 <- read_tsv(file.path(deseq2_dir, "DESeq2_M145_3_vs_2.tsv"), show_col_types = FALSE)
anno <- read_tsv(anno_file, show_col_types = FALSE)

# --- GO term データの準備 ---
go_data <- anno %>%
  filter(!is.na(ontology_term), ontology_term != "") %>%
  select(gene_id, ontology_term) %>%
  separate_rows(ontology_term, sep = ",") %>%
  mutate(ontology_term = trimws(ontology_term)) %>%
  filter(grepl("^GO:", ontology_term))

cat("GO annotations found for", length(unique(go_data$gene_id)), "genes\n")
cat("Total GO term associations:", nrow(go_data), "\n")

# --- 有意なDEGsの抽出 ---
get_degs <- function(df, padj_thresh = 0.05, lfc_thresh = 1, direction = "both") {
  result <- df %>%
    filter(!is.na(padj), padj < padj_thresh, abs(log2FoldChange) >= lfc_thresh)

  if (direction == "up") {
    result <- result %>% filter(log2FoldChange > 0)
  } else if (direction == "down") {
    result <- result %>% filter(log2FoldChange < 0)
  }

  result %>% pull(gene_id)
}

# All genes as background
all_genes <- unique(anno$gene_id)
genes_with_go <- unique(go_data$gene_id)

# --- Fisher's exact test エンリッチメント ---
run_go_enrichment <- function(deg_list, all_genes, go_data, min_genes = 3) {
  # DEGsの中でGOアノテーションがあるもの
  deg_with_go <- intersect(deg_list, genes_with_go)
  bg_with_go <- intersect(all_genes, genes_with_go)

  cat("DEGs with GO annotation:", length(deg_with_go), "/", length(deg_list), "\n")

  # GOタームごとにエンリッチメントテスト
  go_terms <- unique(go_data$ontology_term)

  results <- map_df(go_terms, function(term) {
    genes_in_term <- go_data %>% filter(ontology_term == term) %>% pull(gene_id) %>% unique()

    # 2x2 contingency table
    deg_in_term <- length(intersect(deg_with_go, genes_in_term))
    deg_not_in_term <- length(deg_with_go) - deg_in_term
    bg_in_term <- length(intersect(bg_with_go, genes_in_term)) - deg_in_term
    bg_not_in_term <- length(bg_with_go) - length(deg_with_go) - bg_in_term

    if (deg_in_term < min_genes) return(NULL)

    # Fisher's exact test (one-sided, greater)
    mat <- matrix(c(deg_in_term, bg_in_term, deg_not_in_term, bg_not_in_term), nrow = 2)
    test <- fisher.test(mat, alternative = "greater")

    tibble(
      GO_term = term,
      count = deg_in_term,
      total_in_term = length(genes_in_term),
      gene_ratio = deg_in_term / length(deg_with_go),
      bg_ratio = length(genes_in_term) / length(bg_with_go),
      fold_enrichment = (deg_in_term / length(deg_with_go)) / (length(genes_in_term) / length(bg_with_go)),
      pvalue = test$p.value
    )
  })

  if (nrow(results) == 0) return(NULL)

  results %>%
    mutate(padj = p.adjust(pvalue, method = "BH")) %>%
    filter(padj < 0.05) %>%
    arrange(pvalue)
}

# --- T3 vs T1 のエンリッチメント解析 ---
cat("\n=== T3 vs T1 Up-regulated ===\n")
up_3v1 <- get_degs(de_3v1, direction = "up")
go_up_3v1 <- run_go_enrichment(up_3v1, all_genes, go_data)

cat("\n=== T3 vs T1 Down-regulated ===\n")
down_3v1 <- get_degs(de_3v1, direction = "down")
go_down_3v1 <- run_go_enrichment(down_3v1, all_genes, go_data)

cat("\n=== T2 vs T1 Up-regulated ===\n")
up_2v1 <- get_degs(de_2v1, direction = "up")
go_up_2v1 <- run_go_enrichment(up_2v1, all_genes, go_data)

cat("\n=== T2 vs T1 Down-regulated ===\n")
down_2v1 <- get_degs(de_2v1, direction = "down")
go_down_2v1 <- run_go_enrichment(down_2v1, all_genes, go_data)

# --- バブルチャートの作成 ---
create_bubble_plot <- function(enrich_df, title, output_prefix, top_n = 15) {
  if (is.null(enrich_df) || nrow(enrich_df) == 0) {
    cat("No enrichment results for:", title, "\n")
    return(NULL)
  }

  plot_df <- enrich_df %>%
    slice_head(n = top_n) %>%
    mutate(GO_term = fct_reorder(GO_term, -log10(pvalue)))

  p <- ggplot(plot_df, aes(x = gene_ratio, y = GO_term)) +
    geom_point(aes(size = count, color = -log10(padj))) +
    scale_color_gradient(low = "blue", high = "red", name = "-log10(padj)") +
    scale_size_continuous(range = c(3, 10), name = "Count") +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      plot.title = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = title,
      x = "Gene Ratio",
      y = ""
    )

  ggsave(paste0(output_prefix, ".pdf"), p, width = 10, height = 8)
  ggsave(paste0(output_prefix, ".png"), p, width = 10, height = 8, dpi = 150)

  return(p)
}

# --- バーチャートの作成 ---
create_bar_plot <- function(enrich_df, title, output_prefix, top_n = 15) {
  if (is.null(enrich_df) || nrow(enrich_df) == 0) {
    cat("No enrichment results for:", title, "\n")
    return(NULL)
  }

  plot_df <- enrich_df %>%
    slice_head(n = top_n) %>%
    mutate(GO_term = fct_reorder(GO_term, count))

  p <- ggplot(plot_df, aes(x = count, y = GO_term, fill = -log10(padj))) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "lightblue", high = "darkblue", name = "-log10(padj)") +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      plot.title = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = title,
      x = "Gene Count",
      y = ""
    )

  ggsave(paste0(output_prefix, ".pdf"), p, width = 10, height = 8)
  ggsave(paste0(output_prefix, ".png"), p, width = 10, height = 8, dpi = 150)

  return(p)
}

# --- Figure出力 ---
fig_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")

# T3 vs T1 Up
if (!is.null(go_up_3v1) && nrow(go_up_3v1) > 0) {
  create_bubble_plot(go_up_3v1, "GO Enrichment: Up-regulated (T3 vs T1)",
                     file.path(fig_dir, "GO_bubble_T3vsT1_up"))
  create_bar_plot(go_up_3v1, "GO Enrichment: Up-regulated (T3 vs T1)",
                  file.path(fig_dir, "GO_bar_T3vsT1_up"))
  write_tsv(go_up_3v1, file.path(table_dir, "GO_enrichment_T3vsT1_up.tsv"))
  cat("Saved GO enrichment T3vsT1 Up:", nrow(go_up_3v1), "terms\n")
}

# T3 vs T1 Down
if (!is.null(go_down_3v1) && nrow(go_down_3v1) > 0) {
  create_bubble_plot(go_down_3v1, "GO Enrichment: Down-regulated (T3 vs T1)",
                     file.path(fig_dir, "GO_bubble_T3vsT1_down"))
  create_bar_plot(go_down_3v1, "GO Enrichment: Down-regulated (T3 vs T1)",
                  file.path(fig_dir, "GO_bar_T3vsT1_down"))
  write_tsv(go_down_3v1, file.path(table_dir, "GO_enrichment_T3vsT1_down.tsv"))
  cat("Saved GO enrichment T3vsT1 Down:", nrow(go_down_3v1), "terms\n")
}

# T2 vs T1 Up
if (!is.null(go_up_2v1) && nrow(go_up_2v1) > 0) {
  create_bubble_plot(go_up_2v1, "GO Enrichment: Up-regulated (T2 vs T1)",
                     file.path(fig_dir, "GO_bubble_T2vsT1_up"))
  create_bar_plot(go_up_2v1, "GO Enrichment: Up-regulated (T2 vs T1)",
                  file.path(fig_dir, "GO_bar_T2vsT1_up"))
  write_tsv(go_up_2v1, file.path(table_dir, "GO_enrichment_T2vsT1_up.tsv"))
  cat("Saved GO enrichment T2vsT1 Up:", nrow(go_up_2v1), "terms\n")
}

# T2 vs T1 Down
if (!is.null(go_down_2v1) && nrow(go_down_2v1) > 0) {
  create_bubble_plot(go_down_2v1, "GO Enrichment: Down-regulated (T2 vs T1)",
                     file.path(fig_dir, "GO_bubble_T2vsT1_down"))
  create_bar_plot(go_down_2v1, "GO Enrichment: Down-regulated (T2 vs T1)",
                  file.path(fig_dir, "GO_bar_T2vsT1_down"))
  write_tsv(go_down_2v1, file.path(table_dir, "GO_enrichment_T2vsT1_down.tsv"))
  cat("Saved GO enrichment T2vsT1 Down:", nrow(go_down_2v1), "terms\n")
}

# --- 全結果統合 ---
all_go <- bind_rows(
  if(!is.null(go_up_3v1)) go_up_3v1 %>% mutate(comparison = "T3vsT1", direction = "Up") else NULL,
  if(!is.null(go_down_3v1)) go_down_3v1 %>% mutate(comparison = "T3vsT1", direction = "Down") else NULL,
  if(!is.null(go_up_2v1)) go_up_2v1 %>% mutate(comparison = "T2vsT1", direction = "Up") else NULL,
  if(!is.null(go_down_2v1)) go_down_2v1 %>% mutate(comparison = "T2vsT1", direction = "Down") else NULL
)

if (nrow(all_go) > 0) {
  write_tsv(all_go, file.path(table_dir, "GO_enrichment_all.tsv"))
}

cat("\nGO enrichment analysis completed.\n")
cat("Total enriched GO terms:", nrow(all_go), "\n")
