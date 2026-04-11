#!/usr/bin/env Rscript
# =============================================================================
# 02b_KEGG_enrichment.R
# KEGGパスウェイエンリッチメント解析
# KEGG API経由でS. coelicolor (sco)のパスウェイ情報を取得
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
de_3v1 <- read_tsv(file.path(deseq2_dir, "DESeq2_M145_3_vs_1.tsv"), show_col_types = FALSE)
de_2v1 <- read_tsv(file.path(deseq2_dir, "DESeq2_M145_2_vs_1.tsv"), show_col_types = FALSE)
de_3v2 <- read_tsv(file.path(deseq2_dir, "DESeq2_M145_3_vs_2.tsv"), show_col_types = FALSE)
anno <- read_tsv(anno_file, show_col_types = FALSE)

# --- KEGG pathway データの取得 ---
cat("Fetching KEGG pathway data for S. coelicolor (sco)...\n")

# KEGG REST API からパスウェイリストを取得
kegg_pathways <- tryCatch({
  url <- "https://rest.kegg.jp/list/pathway/sco"
  pathways <- read_tsv(url, col_names = c("pathway_id", "pathway_name"), show_col_types = FALSE)
  pathways <- pathways %>%
    mutate(pathway_id = gsub("path:", "", pathway_id))
  cat("Found", nrow(pathways), "KEGG pathways for sco\n")
  pathways
}, error = function(e) {
  cat("Error fetching KEGG pathways:", e$message, "\n")
  return(NULL)
})

if (is.null(kegg_pathways)) {
  cat("Could not fetch KEGG data. Using local approach.\n")
  quit(status = 0)
}

# --- KEGG gene-pathway mapping の取得 ---
cat("Fetching gene-pathway mappings...\n")

kegg_genes <- tryCatch({
  url <- "https://rest.kegg.jp/link/sco/pathway"
  genes <- read_tsv(url, col_names = c("pathway_id", "gene_id"), show_col_types = FALSE)
  genes <- genes %>%
    mutate(
      pathway_id = gsub("path:", "", pathway_id),
      gene_id = gsub("sco:", "", gene_id)  # SCOxxxx形式
    )
  cat("Found", nrow(genes), "gene-pathway associations\n")
  genes
}, error = function(e) {
  cat("Error fetching gene mappings:", e$message, "\n")
  return(NULL)
})

if (is.null(kegg_genes)) {
  quit(status = 0)
}

# --- old_locus_tag (SCOxxxx) マッピング ---
gene_mapping <- anno %>%
  filter(!is.na(old_locus_tag), old_locus_tag != "") %>%
  select(gene_id, old_locus_tag)

# --- DEGs の抽出 ---
get_degs <- function(df, padj_thresh = 0.05, lfc_thresh = 1, direction = "both") {
  result <- df %>%
    filter(!is.na(padj), padj < padj_thresh, abs(log2FoldChange) >= lfc_thresh) %>%
    left_join(gene_mapping, by = "gene_id") %>%
    filter(!is.na(old_locus_tag))

  if (direction == "up") {
    result <- result %>% filter(log2FoldChange > 0)
  } else if (direction == "down") {
    result <- result %>% filter(log2FoldChange < 0)
  }

  result %>% pull(old_locus_tag)
}

# Background genes (all genes with SCO ID)
background <- gene_mapping$old_locus_tag

# --- KEGG Enrichment (Fisher's exact test) ---
run_kegg_enrichment <- function(deg_list, background, kegg_genes, kegg_pathways, min_genes = 3) {
  cat("DEGs with KEGG annotation:", length(deg_list), "\n")

  pathway_ids <- unique(kegg_genes$pathway_id)

  results <- map_df(pathway_ids, function(pw_id) {
    genes_in_pw <- kegg_genes %>% filter(pathway_id == pw_id) %>% pull(gene_id)

    # 2x2 table
    deg_in_pw <- length(intersect(deg_list, genes_in_pw))
    deg_not_in_pw <- length(deg_list) - deg_in_pw
    bg_in_pw <- length(intersect(background, genes_in_pw)) - deg_in_pw
    bg_not_in_pw <- length(background) - length(deg_list) - bg_in_pw

    if (deg_in_pw < min_genes) return(NULL)

    mat <- matrix(c(deg_in_pw, bg_in_pw, deg_not_in_pw, bg_not_in_pw), nrow = 2)
    test <- fisher.test(mat, alternative = "greater")

    pw_name <- kegg_pathways %>% filter(pathway_id == pw_id) %>% pull(pathway_name)

    tibble(
      pathway_id = pw_id,
      pathway_name = pw_name,
      count = deg_in_pw,
      total_in_pathway = length(genes_in_pw),
      gene_ratio = deg_in_pw / length(deg_list),
      fold_enrichment = (deg_in_pw / length(deg_list)) / (length(genes_in_pw) / length(background)),
      pvalue = test$p.value
    )
  })

  if (nrow(results) == 0) return(NULL)

  results %>%
    mutate(padj = p.adjust(pvalue, method = "BH")) %>%
    filter(padj < 0.1) %>%  # 少し緩めの閾値
    arrange(pvalue)
}

# --- 各比較の解析実行 ---
cat("\n=== T3 vs T1 Up-regulated ===\n")
up_3v1 <- get_degs(de_3v1, direction = "up")
kegg_up_3v1 <- run_kegg_enrichment(up_3v1, background, kegg_genes, kegg_pathways)

cat("\n=== T3 vs T1 Down-regulated ===\n")
down_3v1 <- get_degs(de_3v1, direction = "down")
kegg_down_3v1 <- run_kegg_enrichment(down_3v1, background, kegg_genes, kegg_pathways)

cat("\n=== T2 vs T1 Up-regulated ===\n")
up_2v1 <- get_degs(de_2v1, direction = "up")
kegg_up_2v1 <- run_kegg_enrichment(up_2v1, background, kegg_genes, kegg_pathways)

cat("\n=== T2 vs T1 Down-regulated ===\n")
down_2v1 <- get_degs(de_2v1, direction = "down")
kegg_down_2v1 <- run_kegg_enrichment(down_2v1, background, kegg_genes, kegg_pathways)

cat("\n=== T3 vs T2 Up-regulated ===\n")
up_3v2 <- get_degs(de_3v2, direction = "up")
kegg_up_3v2 <- run_kegg_enrichment(up_3v2, background, kegg_genes, kegg_pathways)

cat("\n=== T3 vs T2 Down-regulated ===\n")
down_3v2 <- get_degs(de_3v2, direction = "down")
kegg_down_3v2 <- run_kegg_enrichment(down_3v2, background, kegg_genes, kegg_pathways)

# --- バブルチャート作成 ---
create_kegg_bubble <- function(enrich_df, title, output_prefix, top_n = 15) {
  if (is.null(enrich_df) || nrow(enrich_df) == 0) {
    cat("No enrichment for:", title, "\n")
    return(NULL)
  }

  plot_df <- enrich_df %>%
    slice_head(n = top_n) %>%
    mutate(pathway_name = str_wrap(pathway_name, width = 40)) %>%
    mutate(pathway_name = fct_reorder(pathway_name, -log10(pvalue)))

  p <- ggplot(plot_df, aes(x = gene_ratio, y = pathway_name)) +
    geom_point(aes(size = count, color = -log10(padj))) +
    scale_color_gradient(low = "blue", high = "red", name = "-log10(padj)") +
    scale_size_continuous(range = c(3, 10), name = "Count") +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 10),
      plot.title = element_text(size = 12, face = "bold")
    ) +
    labs(title = title, x = "Gene Ratio", y = "")

  ggsave(paste0(output_prefix, ".pdf"), p, width = 12, height = 8)
  ggsave(paste0(output_prefix, ".png"), p, width = 12, height = 8, dpi = 150)
  return(p)
}

# --- Figure出力 ---
fig_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")

if (!is.null(kegg_up_3v1) && nrow(kegg_up_3v1) > 0) {
  create_kegg_bubble(kegg_up_3v1, "KEGG Pathway Enrichment: Up-regulated (T3 vs T1)",
                     file.path(fig_dir, "KEGG_bubble_T3vsT1_up"))
  write_tsv(kegg_up_3v1, file.path(table_dir, "KEGG_enrichment_T3vsT1_up.tsv"))
  cat("Saved KEGG T3vsT1 Up:", nrow(kegg_up_3v1), "pathways\n")
}

if (!is.null(kegg_down_3v1) && nrow(kegg_down_3v1) > 0) {
  create_kegg_bubble(kegg_down_3v1, "KEGG Pathway Enrichment: Down-regulated (T3 vs T1)",
                     file.path(fig_dir, "KEGG_bubble_T3vsT1_down"))
  write_tsv(kegg_down_3v1, file.path(table_dir, "KEGG_enrichment_T3vsT1_down.tsv"))
  cat("Saved KEGG T3vsT1 Down:", nrow(kegg_down_3v1), "pathways\n")
}

if (!is.null(kegg_up_2v1) && nrow(kegg_up_2v1) > 0) {
  create_kegg_bubble(kegg_up_2v1, "KEGG Pathway Enrichment: Up-regulated (T2 vs T1)",
                     file.path(fig_dir, "KEGG_bubble_T2vsT1_up"))
  write_tsv(kegg_up_2v1, file.path(table_dir, "KEGG_enrichment_T2vsT1_up.tsv"))
  cat("Saved KEGG T2vsT1 Up:", nrow(kegg_up_2v1), "pathways\n")
}

if (!is.null(kegg_down_2v1) && nrow(kegg_down_2v1) > 0) {
  create_kegg_bubble(kegg_down_2v1, "KEGG Pathway Enrichment: Down-regulated (T2 vs T1)",
                     file.path(fig_dir, "KEGG_bubble_T2vsT1_down"))
  write_tsv(kegg_down_2v1, file.path(table_dir, "KEGG_enrichment_T2vsT1_down.tsv"))
  cat("Saved KEGG T2vsT1 Down:", nrow(kegg_down_2v1), "pathways\n")
}

if (!is.null(kegg_up_3v2) && nrow(kegg_up_3v2) > 0) {
  create_kegg_bubble(kegg_up_3v2, "KEGG Pathway Enrichment: Up-regulated (T3 vs T2)",
                     file.path(fig_dir, "KEGG_bubble_T3vsT2_up"))
  write_tsv(kegg_up_3v2, file.path(table_dir, "KEGG_enrichment_T3vsT2_up.tsv"))
  cat("Saved KEGG T3vsT2 Up:", nrow(kegg_up_3v2), "pathways\n")
}

if (!is.null(kegg_down_3v2) && nrow(kegg_down_3v2) > 0) {
  create_kegg_bubble(kegg_down_3v2, "KEGG Pathway Enrichment: Down-regulated (T3 vs T2)",
                     file.path(fig_dir, "KEGG_bubble_T3vsT2_down"))
  write_tsv(kegg_down_3v2, file.path(table_dir, "KEGG_enrichment_T3vsT2_down.tsv"))
  cat("Saved KEGG T3vsT2 Down:", nrow(kegg_down_3v2), "pathways\n")
}

# --- 統合テーブル ---
all_kegg <- bind_rows(
  if(!is.null(kegg_up_3v1)) kegg_up_3v1 %>% mutate(comparison = "T3vsT1", direction = "Up") else NULL,
  if(!is.null(kegg_down_3v1)) kegg_down_3v1 %>% mutate(comparison = "T3vsT1", direction = "Down") else NULL,
  if(!is.null(kegg_up_2v1)) kegg_up_2v1 %>% mutate(comparison = "T2vsT1", direction = "Up") else NULL,
  if(!is.null(kegg_down_2v1)) kegg_down_2v1 %>% mutate(comparison = "T2vsT1", direction = "Down") else NULL,
  if(!is.null(kegg_up_3v2)) kegg_up_3v2 %>% mutate(comparison = "T3vsT2", direction = "Up") else NULL,
  if(!is.null(kegg_down_3v2)) kegg_down_3v2 %>% mutate(comparison = "T3vsT2", direction = "Down") else NULL
)

if (nrow(all_kegg) > 0) {
  write_tsv(all_kegg, file.path(table_dir, "KEGG_enrichment_all.tsv"))
}

cat("\nKEGG enrichment analysis completed.\n")
cat("Total enriched KEGG pathways:", nrow(all_kegg), "\n")
