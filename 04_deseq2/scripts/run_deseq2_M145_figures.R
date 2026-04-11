#!/usr/bin/env Rscript
# run_deseq2_M145_figures.R
# Regenerate figures from saved RDS objects

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
  library(matrixStats)
  library(svglite)
})

deseq_run_dir <- "/Users/okaban/bioinfo/rna-seq/04_deseq2/analysis/04_deseq2_260128_v1"

# Reload saved objects
cat("Loading saved RDS objects...\n")
dds <- readRDS(file.path(deseq_run_dir, "rds", "dds_raw.rds"))
rld <- readRDS(file.path(deseq_run_dir, "rds", "rld.rds"))

coldata <- as.data.frame(colData(dds))

# Reload DE results from TSV
res_2_vs_1_tbl <- read_tsv(file.path(deseq_run_dir, "results", "DESeq2_M145_2_vs_1.tsv"), show_col_types = FALSE)
res_3_vs_1_tbl <- read_tsv(file.path(deseq_run_dir, "results", "DESeq2_M145_3_vs_1.tsv"), show_col_types = FALSE)
res_3_vs_2_tbl <- read_tsv(file.path(deseq_run_dir, "results", "DESeq2_M145_3_vs_2.tsv"), show_col_types = FALSE)

contrasts_tbl <- list(
  "2_vs_1" = res_2_vs_1_tbl,
  "3_vs_1" = res_3_vs_1_tbl,
  "3_vs_2" = res_3_vs_2_tbl
)

# -------------------------------------------------------------------
# PCA plot
# -------------------------------------------------------------------
cat("Generating PCA plot...\n")
pca_data <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 3) +
  geom_text_repel(size = 2.5, max.overlaps = 20) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA - M145 RNA-seq (rlog)")

ggsave(file.path(deseq_run_dir, "figures", "PCA_M145.pdf"), p_pca, width = 7, height = 5)
ggsave(file.path(deseq_run_dir, "figures", "PCA_M145.svg"), p_pca, width = 7, height = 5)
cat("Saved PCA_M145.pdf/.svg\n")
cat("PCA variance explained: PC1 =", percentVar[1], "%, PC2 =", percentVar[2], "%\n")

# -------------------------------------------------------------------
# Sample distance heatmap
# -------------------------------------------------------------------
cat("Generating sample distance heatmap...\n")
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- coldata$sample_id
colnames(sampleDistMatrix) <- coldata$sample_id
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pdf(file.path(deseq_run_dir, "figures", "sample_distance_heatmap_M145.pdf"), width = 7, height = 6)
pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  col = colors
)
dev.off()

svglite(file.path(deseq_run_dir, "figures", "sample_distance_heatmap_M145.svg"), width = 7, height = 6)
pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  col = colors
)
dev.off()
cat("Saved sample_distance_heatmap_M145.pdf/.svg\n")

# -------------------------------------------------------------------
# Volcano plots
# -------------------------------------------------------------------
cat("Generating volcano plots...\n")
contrast_labels <- c(
  "2_vs_1" = "M145_2 vs M145_1",
  "3_vs_1" = "M145_3 vs M145_1",
  "3_vs_2" = "M145_3 vs M145_2"
)

for (name in names(contrasts_tbl)) {
  res_tbl <- contrasts_tbl[[name]] %>%
    filter(!is.na(padj)) %>%
    mutate(
      neg_log10_padj = -log10(padj),
      sig = case_when(
        padj < 0.05 & log2FoldChange > 1  ~ "Up (|LFC|>1)",
        padj < 0.05 & log2FoldChange < -1 ~ "Down (|LFC|>1)",
        padj < 0.05 ~ "Sig (|LFC|<=1)",
        TRUE ~ "NS"
      )
    )

  # Cap -log10(padj) for visualization
  y_cap <- min(max(res_tbl$neg_log10_padj, na.rm = TRUE), 300)
  res_tbl$neg_log10_padj <- pmin(res_tbl$neg_log10_padj, y_cap)

  # Label top genes
  top_genes <- res_tbl %>%
    filter(sig %in% c("Up (|LFC|>1)", "Down (|LFC|>1)")) %>%
    arrange(padj) %>%
    head(10)

  p_vol <- ggplot(res_tbl, aes(x = log2FoldChange, y = neg_log10_padj, color = sig)) +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_color_manual(values = c(
      "Up (|LFC|>1)" = "red",
      "Down (|LFC|>1)" = "blue",
      "Sig (|LFC|<=1)" = "orange",
      "NS" = "grey70"
    )) +
    geom_text_repel(data = top_genes, aes(label = gene_id),
                    size = 2, max.overlaps = 20, color = "black") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    theme_bw() +
    ggtitle(contrast_labels[name]) +
    xlab("log2 Fold Change (shrunk)") +
    ylab("-log10(padj)")

  ggsave(file.path(deseq_run_dir, "figures", paste0("volcano_M145_", name, ".pdf")),
         p_vol, width = 7, height = 5)
  ggsave(file.path(deseq_run_dir, "figures", paste0("volcano_M145_", name, ".svg")),
         p_vol, width = 7, height = 5)
  cat("Saved volcano_M145_", name, ".pdf/.svg\n", sep = "")
}

# -------------------------------------------------------------------
# Top variable genes heatmap
# -------------------------------------------------------------------
cat("Generating top variable genes heatmap...\n")
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 100)
mat <- assay(rld)[topVarGenes, ]
mat <- mat - rowMeans(mat)

anno_col <- data.frame(
  condition = coldata$condition,
  row.names = coldata$sample_id
)

pdf(file.path(deseq_run_dir, "figures", "topVarGenes_heatmap_M145.pdf"), width = 8, height = 12)
pheatmap(
  mat,
  annotation_col = anno_col,
  show_rownames = TRUE,
  fontsize_row = 5,
  main = "Top 100 variable genes (rlog, row-centered)"
)
dev.off()

svglite(file.path(deseq_run_dir, "figures", "topVarGenes_heatmap_M145.svg"), width = 8, height = 12)
pheatmap(
  mat,
  annotation_col = anno_col,
  show_rownames = TRUE,
  fontsize_row = 5,
  main = "Top 100 variable genes (rlog, row-centered)"
)
dev.off()
cat("Saved topVarGenes_heatmap_M145.pdf/.svg\n")

cat("\n=== All figures generated successfully ===\n")
