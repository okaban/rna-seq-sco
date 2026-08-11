#!/usr/bin/env Rscript
# run_deseq2_M145_s2.R
# DESeq2 analysis using strand-specific counts (-s 2, RF/fr-firststrand)
# Library: NEBNext Ultra II Directional RNA Library Prep Kit (dUTP method)

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
  library(matrixStats)
})

# -------------------------------------------------------------------
# 0. Paths
# -------------------------------------------------------------------
deseq_run_dir <- "/Users/okaban/bioinfo/rna-seq/04_deseq2/results_s2"
counts_file   <- "/Users/okaban/bioinfo/rna-seq/03_quantification/data/featureCounts_s2/featureCounts_M145_s2.txt"
sample_table  <- "/Users/okaban/bioinfo/rna-seq/03_quantification/analysis/03_quant_260128_v1/summary/quant_sample_table.tsv"

log_file <- file.path(deseq_run_dir, "logs", "deseq2_pipeline_s2.log")
sink(log_file, split = TRUE)
cat("=== 04_deseq2 (strand-specific -s 2) started at", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "===\n")
cat("R version:", R.version.string, "\n")
cat("DESeq2 version:", as.character(packageVersion("DESeq2")), "\n")
cat("apeglm version:", as.character(packageVersion("apeglm")), "\n")
cat("Counts file:", counts_file, "\n")
cat("Strandedness: -s 2 (RF / fr-firststrand / reversely stranded)\n")

# -------------------------------------------------------------------
# 1. Read count matrix
# -------------------------------------------------------------------
cat("\n=== Reading count matrix ===\n")
raw <- read.delim(counts_file, comment.char = "#", header = TRUE, check.names = FALSE)
cat("Raw dimensions:", nrow(raw), "genes x", ncol(raw) - 6, "samples\n")

# Extract count matrix (columns 7+)
count_mat <- as.matrix(raw[, 7:ncol(raw)])
rownames(count_mat) <- raw$Geneid

# Clean column names: extract sample_id from BAM paths
colnames(count_mat) <- sub(".*/", "", colnames(count_mat))
colnames(count_mat) <- sub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", colnames(count_mat))
cat("Count matrix:", nrow(count_mat), "genes x", ncol(count_mat), "samples\n")
cat("Samples:", paste(colnames(count_mat), collapse = ", "), "\n")

# -------------------------------------------------------------------
# 2. Read sample metadata
# -------------------------------------------------------------------
cat("\n=== Reading sample table ===\n")
coldata <- read.delim(sample_table, header = TRUE, sep = "\t")
coldata$condition <- factor(coldata$condition, levels = c("M145_1", "M145_2", "M145_3"))
coldata$replicate <- factor(coldata$replicate)
rownames(coldata) <- coldata$sample_id
cat("colData:\n")
print(coldata[, c("sample_id", "condition", "replicate")])

# Ensure column order matches
stopifnot(all(colnames(count_mat) == rownames(coldata)))
cat("Column order check: PASS\n")

# -------------------------------------------------------------------
# 3. Create DESeqDataSet
# -------------------------------------------------------------------
cat("\n=== Creating DESeqDataSet ===\n")
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = coldata,
  design    = ~ condition
)
cat("Before filtering:", nrow(dds), "genes\n")

# Filter: keep genes with >=10 counts in at least 3 samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
cat("After filtering (>=10 counts in >=3 samples):", nrow(dds), "genes\n")

# -------------------------------------------------------------------
# 4. Run DESeq2
# -------------------------------------------------------------------
cat("\n=== Running DESeq2 ===\n")
dds <- DESeq(dds)
cat("Dispersion estimation complete\n")
cat("resultsNames:", paste(resultsNames(dds), collapse = ", "), "\n")

# Save DESeqDataSet
saveRDS(dds, file = file.path(deseq_run_dir, "rds", "dds_raw.rds"))
cat("Saved dds_raw.rds\n")

# -------------------------------------------------------------------
# 5. Normalized counts
# -------------------------------------------------------------------
cat("\n=== Exporting normalized counts ===\n")
norm_counts <- counts(dds, normalized = TRUE)
norm_df <- as.data.frame(norm_counts) %>%
  tibble::rownames_to_column("gene_id")
write_tsv(norm_df, file.path(deseq_run_dir, "results", "normalized_counts_M145.tsv"))
cat("Saved normalized_counts_M145.tsv\n")

cat("Size factors:\n")
print(sizeFactors(dds))

# -------------------------------------------------------------------
# 6. LFC shrinkage — pairwise comparisons
# -------------------------------------------------------------------
cat("\n=== Pairwise comparisons with LFC shrinkage ===\n")

cat("\n--- M145_2 vs M145_1 ---\n")
res_2_vs_1 <- lfcShrink(dds, coef = "condition_M145_2_vs_M145_1", type = "apeglm")
cat("Summary:\n")
summary(res_2_vs_1, alpha = 0.05)

cat("\n--- M145_3 vs M145_1 ---\n")
res_3_vs_1 <- lfcShrink(dds, coef = "condition_M145_3_vs_M145_1", type = "apeglm")
cat("Summary:\n")
summary(res_3_vs_1, alpha = 0.05)

cat("\n--- M145_3 vs M145_2 ---\n")
dds2 <- dds
dds2$condition <- relevel(dds2$condition, ref = "M145_2")
dds2 <- nbinomWaldTest(dds2)
cat("resultsNames (releveled):", paste(resultsNames(dds2), collapse = ", "), "\n")
res_3_vs_2 <- lfcShrink(dds2, coef = "condition_M145_3_vs_M145_2", type = "apeglm")
cat("Summary:\n")
summary(res_3_vs_2, alpha = 0.05)

# Save results
contrasts <- list(
  "2_vs_1" = res_2_vs_1,
  "3_vs_1" = res_3_vs_1,
  "3_vs_2" = res_3_vs_2
)

de_summary <- data.frame(
  contrast = character(),
  total_tested = integer(),
  up_padj005 = integer(),
  down_padj005 = integer(),
  up_padj005_lfc1 = integer(),
  down_padj005_lfc1 = integer(),
  stringsAsFactors = FALSE
)

for (name in names(contrasts)) {
  res <- contrasts[[name]]
  res_tbl <- as_tibble(as.data.frame(res), rownames = "gene_id")
  write_tsv(res_tbl, file.path(deseq_run_dir, "results", paste0("DESeq2_M145_", name, ".tsv")))
  cat("Saved DESeq2_M145_", name, ".tsv\n", sep = "")

  sig <- res_tbl %>% filter(!is.na(padj))
  up005   <- sum(sig$padj < 0.05 & sig$log2FoldChange > 0, na.rm = TRUE)
  down005 <- sum(sig$padj < 0.05 & sig$log2FoldChange < 0, na.rm = TRUE)
  up005_lfc1   <- sum(sig$padj < 0.05 & sig$log2FoldChange > 1, na.rm = TRUE)
  down005_lfc1 <- sum(sig$padj < 0.05 & sig$log2FoldChange < -1, na.rm = TRUE)

  de_summary <- rbind(de_summary, data.frame(
    contrast = name,
    total_tested = nrow(sig),
    up_padj005 = up005,
    down_padj005 = down005,
    up_padj005_lfc1 = up005_lfc1,
    down_padj005_lfc1 = down005_lfc1,
    stringsAsFactors = FALSE
  ))
}

cat("\n=== DE summary ===\n")
print(de_summary)
write_tsv(de_summary, file.path(deseq_run_dir, "results", "DE_summary_M145.tsv"))

# -------------------------------------------------------------------
# 7. QC Visualizations
# -------------------------------------------------------------------
cat("\n=== Generating QC figures ===\n")

cat("Running rlog transformation...\n")
rld <- rlog(dds, blind = FALSE)
saveRDS(rld, file = file.path(deseq_run_dir, "rds", "rld.rds"))
cat("Saved rld.rds\n")

# PCA plot
cat("Generating PCA plot...\n")
pca_data <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 3) +
  geom_text_repel(size = 2.5, max.overlaps = 20) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA — M145 RNA-seq (rlog, strand-specific -s 2)")

ggsave(file.path(deseq_run_dir, "figures", "PCA_M145.pdf"), p_pca, width = 7, height = 5)
ggsave(file.path(deseq_run_dir, "figures", "PCA_M145.svg"), p_pca, width = 7, height = 5)
cat("Saved PCA_M145.pdf/.svg\n")

# Sample distance heatmap
cat("Generating sample distance heatmap...\n")
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- coldata$sample_id
colnames(sampleDistMatrix) <- coldata$sample_id
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pdf(file.path(deseq_run_dir, "figures", "sample_distance_heatmap_M145.pdf"), width = 7, height = 6)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.off()

svg(file.path(deseq_run_dir, "figures", "sample_distance_heatmap_M145.svg"), width = 7, height = 6)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.off()
cat("Saved sample_distance_heatmap_M145.pdf/.svg\n")

# Volcano plots
cat("Generating volcano plots...\n")
contrast_labels <- c(
  "2_vs_1" = "M145_2 vs M145_1",
  "3_vs_1" = "M145_3 vs M145_1",
  "3_vs_2" = "M145_3 vs M145_2"
)

for (name in names(contrasts)) {
  res <- contrasts[[name]]
  res_tbl <- as_tibble(as.data.frame(res), rownames = "gene_id") %>%
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

  y_cap <- min(max(res_tbl$neg_log10_padj, na.rm = TRUE), 300)
  res_tbl$neg_log10_padj <- pmin(res_tbl$neg_log10_padj, y_cap)

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

# Top variable genes heatmap
cat("Generating top variable genes heatmap...\n")
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 100)
mat <- assay(rld)[topVarGenes, ]
mat <- mat - rowMeans(mat)

anno_col <- data.frame(
  condition = coldata$condition,
  row.names = coldata$sample_id
)

pdf(file.path(deseq_run_dir, "figures", "topVarGenes_heatmap_M145.pdf"), width = 8, height = 12)
pheatmap(mat, annotation_col = anno_col, show_rownames = TRUE,
         fontsize_row = 5, main = "Top 100 variable genes (rlog, row-centered, -s 2)")
dev.off()

svg(file.path(deseq_run_dir, "figures", "topVarGenes_heatmap_M145.svg"), width = 8, height = 12)
pheatmap(mat, annotation_col = anno_col, show_rownames = TRUE,
         fontsize_row = 5, main = "Top 100 variable genes (rlog, row-centered, -s 2)")
dev.off()
cat("Saved topVarGenes_heatmap_M145.pdf/.svg\n")

# -------------------------------------------------------------------
# 8. Session info
# -------------------------------------------------------------------
cat("\n=== Session info ===\n")
sessionInfo()

cat("\n=== 04_deseq2 (-s 2) COMPLETED at", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "===\n")
sink()
