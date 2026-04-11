#!/usr/bin/env Rscript
# =============================================================================
# 01_venn_diagram_DEGs.R
# ベン図：3比較間のDEGsの重複を可視化
# Project: M145 RNA-seq
# Date: 2026-02-02
# =============================================================================

library(tidyverse)
library(VennDiagram)
library(grid)

# --- パス設定 ---
deseq2_dir <- "/Users/okaban/bioinfo/rna-seq/04_deseq2/analysis/04_deseq2_260128_v1/results"
output_dir <- "/Users/okaban/bioinfo/rna-seq/12_supplementary_figures/analysis/12_supplementary_260202_v1"

# --- DESeq2結果の読み込み ---
de_2v1 <- read_tsv(file.path(deseq2_dir, "DESeq2_M145_2_vs_1.tsv"), show_col_types = FALSE)
de_3v1 <- read_tsv(file.path(deseq2_dir, "DESeq2_M145_3_vs_1.tsv"), show_col_types = FALSE)
de_3v2 <- read_tsv(file.path(deseq2_dir, "DESeq2_M145_3_vs_2.tsv"), show_col_types = FALSE)

# --- 有意なDEGsの抽出（padj < 0.05, |log2FC| > 1）---
get_degs <- function(df, padj_thresh = 0.05, lfc_thresh = 1) {
  df %>%
    filter(!is.na(padj), padj < padj_thresh, abs(log2FoldChange) >= lfc_thresh) %>%
    pull(gene_id)
}

degs_2v1 <- get_degs(de_2v1)
degs_3v1 <- get_degs(de_3v1)
degs_3v2 <- get_degs(de_3v2)

cat("DEGs count:\n")
cat("  2_vs_1:", length(degs_2v1), "\n")
cat("  3_vs_1:", length(degs_3v1), "\n")
cat("  3_vs_2:", length(degs_3v2), "\n")

# --- Up-regulated DEGsの抽出 ---
get_up_degs <- function(df, padj_thresh = 0.05, lfc_thresh = 1) {
  df %>%
    filter(!is.na(padj), padj < padj_thresh, log2FoldChange >= lfc_thresh) %>%
    pull(gene_id)
}

get_down_degs <- function(df, padj_thresh = 0.05, lfc_thresh = 1) {
  df %>%
    filter(!is.na(padj), padj < padj_thresh, log2FoldChange <= -lfc_thresh) %>%
    pull(gene_id)
}

up_2v1 <- get_up_degs(de_2v1)
up_3v1 <- get_up_degs(de_3v1)
up_3v2 <- get_up_degs(de_3v2)

down_2v1 <- get_down_degs(de_2v1)
down_3v1 <- get_down_degs(de_3v1)
down_3v2 <- get_down_degs(de_3v2)

# --- Venn Diagram (All DEGs) ---
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

venn_all <- venn.diagram(
  x = list(
    "T2 vs T1" = degs_2v1,
    "T3 vs T1" = degs_3v1,
    "T3 vs T2" = degs_3v2
  ),
  filename = NULL,
  fill = c("#E41A1C", "#377EB8", "#4DAF4A"),
  alpha = 0.5,
  cat.cex = 1.2,
  cex = 1.5,
  cat.fontface = "bold",
  main = "Differentially Expressed Genes (|log2FC| ≥ 1, padj < 0.05)",
  main.cex = 1.3
)

pdf(file.path(output_dir, "figures", "venn_all_DEGs.pdf"), width = 8, height = 8)
grid.draw(venn_all)
dev.off()

png(file.path(output_dir, "figures", "venn_all_DEGs.png"), width = 800, height = 800, res = 150)
grid.draw(venn_all)
dev.off()

# --- Venn Diagram (Up-regulated DEGs) ---
venn_up <- venn.diagram(
  x = list(
    "T2 vs T1" = up_2v1,
    "T3 vs T1" = up_3v1,
    "T3 vs T2" = up_3v2
  ),
  filename = NULL,
  fill = c("#FC8D62", "#8DA0CB", "#66C2A5"),
  alpha = 0.5,
  cat.cex = 1.2,
  cex = 1.5,
  cat.fontface = "bold",
  main = "Up-regulated Genes (log2FC ≥ 1, padj < 0.05)",
  main.cex = 1.3
)

pdf(file.path(output_dir, "figures", "venn_up_DEGs.pdf"), width = 8, height = 8)
grid.draw(venn_up)
dev.off()

png(file.path(output_dir, "figures", "venn_up_DEGs.png"), width = 800, height = 800, res = 150)
grid.draw(venn_up)
dev.off()

# --- Venn Diagram (Down-regulated DEGs) ---
venn_down <- venn.diagram(
  x = list(
    "T2 vs T1" = down_2v1,
    "T3 vs T1" = down_3v1,
    "T3 vs T2" = down_3v2
  ),
  filename = NULL,
  fill = c("#E78AC3", "#A6D854", "#FFD92F"),
  alpha = 0.5,
  cat.cex = 1.2,
  cex = 1.5,
  cat.fontface = "bold",
  main = "Down-regulated Genes (log2FC ≤ -1, padj < 0.05)",
  main.cex = 1.3
)

pdf(file.path(output_dir, "figures", "venn_down_DEGs.pdf"), width = 8, height = 8)
grid.draw(venn_down)
dev.off()

png(file.path(output_dir, "figures", "venn_down_DEGs.png"), width = 800, height = 800, res = 150)
grid.draw(venn_down)
dev.off()

# --- 重複遺伝子のサマリーテーブル出力 ---
# Calculate overlaps
all_genes <- unique(c(degs_2v1, degs_3v1, degs_3v2))
overlap_df <- data.frame(
  gene_id = all_genes,
  in_2v1 = all_genes %in% degs_2v1,
  in_3v1 = all_genes %in% degs_3v1,
  in_3v2 = all_genes %in% degs_3v2
) %>%
  mutate(
    overlap_count = in_2v1 + in_3v1 + in_3v2,
    category = case_when(
      in_2v1 & in_3v1 & in_3v2 ~ "All three",
      in_2v1 & in_3v1 ~ "2v1 & 3v1 only",
      in_2v1 & in_3v2 ~ "2v1 & 3v2 only",
      in_3v1 & in_3v2 ~ "3v1 & 3v2 only",
      in_2v1 ~ "2v1 only",
      in_3v1 ~ "3v1 only",
      in_3v2 ~ "3v2 only"
    )
  )

write_tsv(overlap_df, file.path(output_dir, "tables", "venn_DEGs_overlap.tsv"))

# Summary statistics
summary_stats <- overlap_df %>%
  count(category) %>%
  arrange(desc(n))

write_tsv(summary_stats, file.path(output_dir, "tables", "venn_DEGs_summary.tsv"))

cat("\nOverlap summary:\n")
print(summary_stats)

cat("\nVenn diagrams saved to:", file.path(output_dir, "figures"), "\n")
