#!/usr/bin/env Rscript
# =============================================================================
# 02_GO_KEGG_enrichment.R
# GO/KEGGエンリッチメント解析
# Project: M145 RNA-seq
# Date: 2026-02-02
# Organism: Streptomyces coelicolor A3(2) M145 (KEGG: sco)
# =============================================================================

library(tidyverse)
library(clusterProfiler)
library(enrichplot)
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

# --- old_locus_tag to gene_id mapping (for KEGG) ---
# S. coelicolorのKEGGではSCOxxxxの形式を使用
gene2sco <- anno %>%
  filter(!is.na(old_locus_tag), old_locus_tag != "") %>%
  select(gene_id, old_locus_tag) %>%
  mutate(kegg_id = paste0("sco:", old_locus_tag))

cat("Gene ID mapping created:", nrow(gene2sco), "genes\n")

# --- 有意なDEGsの抽出 ---
get_degs_info <- function(df, padj_thresh = 0.05, lfc_thresh = 1) {
  df %>%
    filter(!is.na(padj), padj < padj_thresh, abs(log2FoldChange) >= lfc_thresh) %>%
    left_join(gene2sco, by = "gene_id") %>%
    filter(!is.na(kegg_id))
}

degs_2v1 <- get_degs_info(de_2v1)
degs_3v1 <- get_degs_info(de_3v1)
degs_3v2 <- get_degs_info(de_3v2)

# Up/Down分け
up_2v1 <- degs_2v1 %>% filter(log2FoldChange > 0) %>% pull(old_locus_tag)
down_2v1 <- degs_2v1 %>% filter(log2FoldChange < 0) %>% pull(old_locus_tag)
up_3v1 <- degs_3v1 %>% filter(log2FoldChange > 0) %>% pull(old_locus_tag)
down_3v1 <- degs_3v1 %>% filter(log2FoldChange < 0) %>% pull(old_locus_tag)
up_3v2 <- degs_3v2 %>% filter(log2FoldChange > 0) %>% pull(old_locus_tag)
down_3v2 <- degs_3v2 %>% filter(log2FoldChange < 0) %>% pull(old_locus_tag)

cat("\nDEGs with KEGG ID:\n")
cat("  2_vs_1 - Up:", length(up_2v1), "Down:", length(down_2v1), "\n")
cat("  3_vs_1 - Up:", length(up_3v1), "Down:", length(down_3v1), "\n")
cat("  3_vs_2 - Up:", length(up_3v2), "Down:", length(down_3v2), "\n")

# --- Background gene list ---
background_genes <- gene2sco$old_locus_tag

# --- KEGG Enrichment Analysis ---
run_kegg_enrichment <- function(gene_list, label, direction = "all") {
  cat("Running KEGG enrichment for:", label, "(", direction, ")\n")

  if (length(gene_list) < 10) {
    cat("  Skipping - too few genes\n")
    return(NULL)
  }

  tryCatch({
    kegg_res <- enrichKEGG(
      gene = gene_list,
      organism = "sco",
      keyType = "kegg",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      universe = background_genes
    )

    if (is.null(kegg_res) || nrow(kegg_res@result) == 0) {
      cat("  No significant enrichment found\n")
      return(NULL)
    }

    kegg_res@result$comparison <- label
    kegg_res@result$direction <- direction
    return(kegg_res)
  }, error = function(e) {
    cat("  Error:", e$message, "\n")
    return(NULL)
  })
}

# 各比較でKEGGエンリッチメント実行
kegg_3v1_up <- run_kegg_enrichment(up_3v1, "T3_vs_T1", "Up")
kegg_3v1_down <- run_kegg_enrichment(down_3v1, "T3_vs_T1", "Down")
kegg_2v1_up <- run_kegg_enrichment(up_2v1, "T2_vs_T1", "Up")
kegg_2v1_down <- run_kegg_enrichment(down_2v1, "T2_vs_T1", "Down")
kegg_3v2_up <- run_kegg_enrichment(up_3v2, "T3_vs_T2", "Up")
kegg_3v2_down <- run_kegg_enrichment(down_3v2, "T3_vs_T2", "Down")

# --- 結果の統合と可視化 ---
# T3 vs T1 (最も変動が大きい比較) のエンリッチメントを可視化
if (!is.null(kegg_3v1_up)) {
  # Dot plot (bubble chart)
  p1 <- dotplot(kegg_3v1_up, showCategory = 15,
                title = "KEGG Enrichment: Up-regulated (T3 vs T1)") +
    theme(axis.text.y = element_text(size = 10))
  ggsave(file.path(output_dir, "figures", "KEGG_dotplot_T3vsT1_up.pdf"), p1,
         width = 10, height = 8)
  ggsave(file.path(output_dir, "figures", "KEGG_dotplot_T3vsT1_up.png"), p1,
         width = 10, height = 8, dpi = 150)

  # Bar plot
  p2 <- barplot(kegg_3v1_up, showCategory = 15,
                title = "KEGG Enrichment: Up-regulated (T3 vs T1)") +
    theme(axis.text.y = element_text(size = 10))
  ggsave(file.path(output_dir, "figures", "KEGG_barplot_T3vsT1_up.pdf"), p2,
         width = 10, height = 8)
  ggsave(file.path(output_dir, "figures", "KEGG_barplot_T3vsT1_up.png"), p2,
         width = 10, height = 8, dpi = 150)

  # Save table
  write_tsv(kegg_3v1_up@result,
            file.path(output_dir, "tables", "KEGG_T3vsT1_up.tsv"))
}

if (!is.null(kegg_3v1_down)) {
  p3 <- dotplot(kegg_3v1_down, showCategory = 15,
                title = "KEGG Enrichment: Down-regulated (T3 vs T1)") +
    theme(axis.text.y = element_text(size = 10))
  ggsave(file.path(output_dir, "figures", "KEGG_dotplot_T3vsT1_down.pdf"), p3,
         width = 10, height = 8)
  ggsave(file.path(output_dir, "figures", "KEGG_dotplot_T3vsT1_down.png"), p3,
         width = 10, height = 8, dpi = 150)

  p4 <- barplot(kegg_3v1_down, showCategory = 15,
                title = "KEGG Enrichment: Down-regulated (T3 vs T1)") +
    theme(axis.text.y = element_text(size = 10))
  ggsave(file.path(output_dir, "figures", "KEGG_barplot_T3vsT1_down.pdf"), p4,
         width = 10, height = 8)
  ggsave(file.path(output_dir, "figures", "KEGG_barplot_T3vsT1_down.png"), p4,
         width = 10, height = 8, dpi = 150)

  write_tsv(kegg_3v1_down@result,
            file.path(output_dir, "tables", "KEGG_T3vsT1_down.tsv"))
}

# T2 vs T1
if (!is.null(kegg_2v1_up)) {
  p5 <- dotplot(kegg_2v1_up, showCategory = 15,
                title = "KEGG Enrichment: Up-regulated (T2 vs T1)") +
    theme(axis.text.y = element_text(size = 10))
  ggsave(file.path(output_dir, "figures", "KEGG_dotplot_T2vsT1_up.pdf"), p5,
         width = 10, height = 8)
  write_tsv(kegg_2v1_up@result,
            file.path(output_dir, "tables", "KEGG_T2vsT1_up.tsv"))
}

if (!is.null(kegg_2v1_down)) {
  p6 <- dotplot(kegg_2v1_down, showCategory = 15,
                title = "KEGG Enrichment: Down-regulated (T2 vs T1)") +
    theme(axis.text.y = element_text(size = 10))
  ggsave(file.path(output_dir, "figures", "KEGG_dotplot_T2vsT1_down.pdf"), p6,
         width = 10, height = 8)
  write_tsv(kegg_2v1_down@result,
            file.path(output_dir, "tables", "KEGG_T2vsT1_down.tsv"))
}

# --- 全結果の統合テーブル ---
all_kegg_results <- bind_rows(
  if(!is.null(kegg_3v1_up)) kegg_3v1_up@result else NULL,
  if(!is.null(kegg_3v1_down)) kegg_3v1_down@result else NULL,
  if(!is.null(kegg_2v1_up)) kegg_2v1_up@result else NULL,
  if(!is.null(kegg_2v1_down)) kegg_2v1_down@result else NULL,
  if(!is.null(kegg_3v2_up)) kegg_3v2_up@result else NULL,
  if(!is.null(kegg_3v2_down)) kegg_3v2_down@result else NULL
)

if (nrow(all_kegg_results) > 0) {
  write_tsv(all_kegg_results,
            file.path(output_dir, "tables", "KEGG_all_comparisons.tsv"))
  cat("\nKEGG enrichment results saved\n")
  cat("Total enriched pathways:", nrow(all_kegg_results), "\n")
} else {
  cat("\nNo significant KEGG enrichment found\n")
}

# --- GO Enrichment (using annotation file GO terms) ---
cat("\n--- GO Enrichment Analysis ---\n")

# Parse GO terms from annotation
go_data <- anno %>%
  filter(!is.na(ontology_term), ontology_term != "") %>%
  select(gene_id, ontology_term) %>%
  separate_rows(ontology_term, sep = ",") %>%
  mutate(ontology_term = trimws(ontology_term)) %>%
  filter(grepl("^GO:", ontology_term))

if (nrow(go_data) > 0) {
  cat("GO annotations found for", length(unique(go_data$gene_id)), "genes\n")
  cat("Total GO term associations:", nrow(go_data), "\n")

  # Create TERM2GENE mapping
  term2gene <- go_data %>%
    select(ontology_term, gene_id) %>%
    rename(term = ontology_term, gene = gene_id)

  # Run GO enrichment for T3 vs T1 (Up)
  up_3v1_genes <- degs_3v1 %>% filter(log2FoldChange > 0) %>% pull(gene_id)

  if (length(up_3v1_genes) >= 10) {
    go_up <- enricher(
      gene = up_3v1_genes,
      TERM2GENE = term2gene,
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )

    if (!is.null(go_up) && nrow(go_up@result) > 0) {
      p_go <- dotplot(go_up, showCategory = 15,
                      title = "GO Enrichment: Up-regulated (T3 vs T1)") +
        theme(axis.text.y = element_text(size = 9))
      ggsave(file.path(output_dir, "figures", "GO_dotplot_T3vsT1_up.pdf"), p_go,
             width = 10, height = 8)
      ggsave(file.path(output_dir, "figures", "GO_dotplot_T3vsT1_up.png"), p_go,
             width = 10, height = 8, dpi = 150)
      write_tsv(go_up@result,
                file.path(output_dir, "tables", "GO_T3vsT1_up.tsv"))
      cat("GO enrichment (Up T3vsT1): saved\n")
    }
  }

  # Down-regulated
  down_3v1_genes <- degs_3v1 %>% filter(log2FoldChange < 0) %>% pull(gene_id)

  if (length(down_3v1_genes) >= 10) {
    go_down <- enricher(
      gene = down_3v1_genes,
      TERM2GENE = term2gene,
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )

    if (!is.null(go_down) && nrow(go_down@result) > 0) {
      p_go2 <- dotplot(go_down, showCategory = 15,
                       title = "GO Enrichment: Down-regulated (T3 vs T1)") +
        theme(axis.text.y = element_text(size = 9))
      ggsave(file.path(output_dir, "figures", "GO_dotplot_T3vsT1_down.pdf"), p_go2,
             width = 10, height = 8)
      ggsave(file.path(output_dir, "figures", "GO_dotplot_T3vsT1_down.png"), p_go2,
             width = 10, height = 8, dpi = 150)
      write_tsv(go_down@result,
                file.path(output_dir, "tables", "GO_T3vsT1_down.tsv"))
      cat("GO enrichment (Down T3vsT1): saved\n")
    }
  }
}

cat("\nGO/KEGG enrichment analysis completed.\n")
cat("Output directory:", file.path(output_dir, "figures"), "\n")
