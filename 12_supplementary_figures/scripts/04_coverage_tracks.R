#!/usr/bin/env Rscript
# =============================================================================
# 04_coverage_tracks.R
# ゲノムカバレッジトラックビュー
# 主要BGC領域（act, red, cpk, cda）のリードカバレッジを可視化
# Project: M145 RNA-seq
# Date: 2026-02-02
# =============================================================================

library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggrepel)

# --- パス設定 ---
bam_dir <- "/Users/okaban/bioinfo/rna-seq/02_alignment/analysis/02_alignment_260127_v1/bam"
anno_file <- "/Users/okaban/bioinfo/rna-seq/05_annotation/analysis/05_annotation_260128_v1/tables/gene_annotation_basic.tsv"
bgc_file <- "/Users/okaban/bioinfo/rna-seq/05_annotation/analysis/05_annotation_260128_v1/tables/BGC_definition_manual.tsv"
output_dir <- "/Users/okaban/bioinfo/rna-seq/12_supplementary_figures/analysis/12_supplementary_260202_v1"

# --- BGC座標の定義 ---
bgc_regions <- tibble(
  bgc = c("act", "red", "cda", "cpk"),
  chrom = "NC_003888.3",
  start = c(5513809, 6432812, 3519449, 6900898),
  end = c(5535091, 6464206, 3602320, 6948414),
  description = c("Actinorhodin", "Undecylprodigiosin", "CDA", "Coelimycin P1")
)

# サンプル情報
samples <- tibble(
  sample_id = c("M145_1_1", "M145_1_2", "M145_1_3",
                "M145_2_1", "M145_2_3", "M145_2_4",
                "M145_3_2", "M145_3_3", "M145_3_4"),
  condition = c(rep("T1", 3), rep("T2", 3), rep("T3", 3))
)

# --- アノテーション読み込み ---
anno <- read_tsv(anno_file, show_col_types = FALSE)
bgc_genes <- read_tsv(bgc_file, show_col_types = FALSE)

# --- samtoolsでカバレッジ計算 ---
calculate_coverage <- function(bam_file, region, sample_name, bin_size = 100) {
  temp_file <- tempfile(fileext = ".depth")

  cmd <- sprintf("samtools depth -r %s %s > %s 2>/dev/null",
                 region, bam_file, temp_file)
  system(cmd)

  if (file.exists(temp_file) && file.size(temp_file) > 0) {
    depth <- read_tsv(temp_file, col_names = c("chrom", "pos", "depth"),
                      col_types = "cii", progress = FALSE)
    unlink(temp_file)

    if (nrow(depth) > 0) {
      # ビン化して平滑化
      depth %>%
        mutate(bin = floor(pos / bin_size) * bin_size) %>%
        group_by(bin) %>%
        summarise(depth = mean(depth), .groups = "drop") %>%
        rename(pos = bin) %>%
        mutate(sample = sample_name)
    } else {
      return(NULL)
    }
  } else {
    unlink(temp_file)
    return(NULL)
  }
}

# --- 各BGC領域のカバレッジを計算 ---
process_bgc <- function(bgc_name, bgc_info, samples_df, bam_dir) {
  region <- sprintf("%s:%d-%d", bgc_info$chrom, bgc_info$start, bgc_info$end)
  cat("Processing", bgc_name, "region:", region, "\n")

  all_coverage <- map_df(1:nrow(samples_df), function(i) {
    sample <- samples_df$sample_id[i]
    condition <- samples_df$condition[i]
    bam_file <- file.path(bam_dir, paste0(sample, ".Aligned.sortedByCoord.out.bam"))

    if (!file.exists(bam_file)) {
      cat("  Warning: BAM file not found for", sample, "\n")
      return(NULL)
    }

    cov <- calculate_coverage(bam_file, region, sample, bin_size = 50)
    if (!is.null(cov)) {
      cov$condition <- condition
    }
    cov
  })

  return(all_coverage)
}

# --- 各BGCのカバレッジを計算 ---
coverage_data <- list()
for (i in 1:nrow(bgc_regions)) {
  bgc_name <- bgc_regions$bgc[i]
  coverage_data[[bgc_name]] <- process_bgc(
    bgc_name,
    bgc_regions[i, ],
    samples,
    bam_dir
  )
}

# --- 矢印形状の遺伝子ポリゴンを作成する関数 ---
create_arrow_polygon <- function(start, end, strand, y_pos, height = 0.35, arrow_prop = 0.15) {
  gene_length <- end - start
  arrow_len <- min(gene_length * arrow_prop, gene_length * 0.3)

  if (strand == "+") {
    # 右向き矢印
    tibble(
      x = c(start, end - arrow_len, end, end - arrow_len, start, start),
      y = c(y_pos - height, y_pos - height, y_pos, y_pos + height, y_pos + height, y_pos - height)
    )
  } else {
    # 左向き矢印
    tibble(
      x = c(start + arrow_len, end, end, start + arrow_len, start, start + arrow_len),
      y = c(y_pos - height, y_pos - height, y_pos + height, y_pos + height, y_pos, y_pos - height)
    )
  }
}

# --- カバレッジプロット作成 ---
create_coverage_plot <- function(cov_df, bgc_name, bgc_info, anno_df, bgc_genes_df) {
  if (is.null(cov_df) || nrow(cov_df) == 0) {
    cat("No coverage data for", bgc_name, "\n")
    return(NULL)
  }

  # 条件ごとの平均カバレッジを計算
  cov_mean <- cov_df %>%
    group_by(pos, condition) %>%
    summarise(
      mean_depth = mean(depth),
      sd_depth = sd(depth),
      .groups = "drop"
    )

  # 遺伝子アノテーションの取得
  genes_in_region <- anno_df %>%
    filter(contig == bgc_info$chrom,
           start >= bgc_info$start - 1000,
           end <= bgc_info$end + 1000) %>%
    left_join(bgc_genes_df %>% filter(bgc_name == !!bgc_name) %>% select(gene_id, role, mibig_annotation),
              by = "gene_id")

  # カバレッジプロット
  p1 <- ggplot(cov_mean, aes(x = pos, y = mean_depth, color = condition, fill = condition)) +
    geom_ribbon(aes(ymin = pmax(0, mean_depth - sd_depth),
                    ymax = mean_depth + sd_depth), alpha = 0.2, color = NA) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = c("T1" = "#66C2A5", "T2" = "#FC8D62", "T3" = "#8DA0CB")) +
    scale_fill_manual(values = c("T1" = "#66C2A5", "T2" = "#FC8D62", "T3" = "#8DA0CB")) +
    theme_bw() +
    labs(
      title = sprintf("%s (%s) Coverage", toupper(bgc_name), bgc_info$description),
      x = "Genomic Position",
      y = "Mean Coverage (reads)",
      color = "Timepoint",
      fill = "Timepoint"
    ) +
    theme(
      legend.position = "top",
      plot.title = element_text(size = 14, face = "bold")
    )

  # 遺伝子トラック（矢印形状）
  if (nrow(genes_in_region) > 0) {
    gene_track <- genes_in_region %>%
      mutate(
        y_pos = ifelse(strand == "+", 1, -1),
        # MiBIGアノテーションがあればそれを使用、なければgene_name、それもなければSCO番号
        gene_label = case_when(
          !is.na(mibig_annotation) & mibig_annotation != "NA" ~ mibig_annotation,
          !is.na(gene_name) & gene_name != "NA" & gene_name != "" ~ gene_name,
          !is.na(old_locus_tag) & old_locus_tag != "" ~ old_locus_tag,
          TRUE ~ ""
        ),
        gene_id_row = row_number()
      )

    # 各遺伝子の矢印ポリゴンデータを作成
    arrow_data <- gene_track %>%
      rowwise() %>%
      do({
        poly <- create_arrow_polygon(.$start, .$end, .$strand, .$y_pos)
        poly$gene_id_row <- .$gene_id_row
        poly$role <- .$role
        poly
      }) %>%
      ungroup() %>%
      left_join(gene_track %>% select(gene_id_row, role), by = c("gene_id_row", "role"))

    # ラベル用のデータ（重なり防止）
    label_data <- gene_track %>%
      filter(gene_label != "") %>%
      mutate(
        label_x = (start + end) / 2,
        label_y = ifelse(strand == "+", y_pos + 0.7, y_pos - 0.7)
      )

    p2 <- ggplot() +
      # 矢印形状の遺伝子
      geom_polygon(data = arrow_data,
                   aes(x = x, y = y, group = gene_id_row, fill = role),
                   color = "black", linewidth = 0.4) +
      # 中央線（strandを示す補助線）
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
      # 遺伝子ラベル（ggrepelで重なり防止）
      geom_label_repel(data = label_data,
                       aes(x = label_x, y = label_y, label = gene_label),
                       size = 2.8, fontface = "bold",
                       fill = "white", alpha = 0.85,
                       label.padding = unit(0.15, "lines"),
                       label.size = 0.2,
                       segment.size = 0.3, segment.color = "gray40",
                       max.overlaps = 30,
                       direction = "both",
                       nudge_y = ifelse(label_data$y_pos > 0, 0.3, -0.3),
                       min.segment.length = 0) +
      scale_fill_manual(values = c("biosynthesis" = "#3498db",
                                   "regulator" = "#e74c3c",
                                   "transport" = "#27ae60",
                                   "resistance" = "#f39c12"),
                        na.value = "gray70",
                        name = "Gene Function") +
      scale_y_continuous(limits = c(-2.5, 2.5),
                         breaks = c(-1, 1),
                         labels = c("(-) strand", "(+) strand")) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 9, face = "italic"),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 10, face = "bold")
      ) +
      labs(x = "", y = "") +
      coord_cartesian(xlim = c(bgc_info$start, bgc_info$end))

    # 結合
    combined <- p1 / p2 + plot_layout(heights = c(3, 1.5))
    return(combined)
  } else {
    return(p1)
  }
}

# --- Figure出力 ---
fig_dir <- file.path(output_dir, "figures")

for (i in 1:nrow(bgc_regions)) {
  bgc_name <- bgc_regions$bgc[i]
  bgc_info <- bgc_regions[i, ]

  p <- create_coverage_plot(
    coverage_data[[bgc_name]],
    bgc_name,
    bgc_info,
    anno,
    bgc_genes
  )

  if (!is.null(p)) {
    ggsave(file.path(fig_dir, sprintf("coverage_track_%s.pdf", bgc_name)),
           p, width = 14, height = 8)
    ggsave(file.path(fig_dir, sprintf("coverage_track_%s.png", bgc_name)),
           p, width = 14, height = 8, dpi = 150)
    cat("Saved coverage track for", bgc_name, "\n")
  }
}

# --- 4つのBGCを1つのFigureにまとめる ---
cat("\nCreating combined coverage figure...\n")

create_simple_coverage <- function(cov_df, bgc_name, bgc_info) {
  if (is.null(cov_df) || nrow(cov_df) == 0) return(NULL)

  cov_mean <- cov_df %>%
    group_by(pos, condition) %>%
    summarise(mean_depth = mean(depth), .groups = "drop")

  ggplot(cov_mean, aes(x = pos / 1000, y = mean_depth, color = condition)) +
    geom_line(linewidth = 0.6) +
    scale_color_manual(values = c("T1" = "#66C2A5", "T2" = "#FC8D62", "T3" = "#8DA0CB")) +
    theme_bw() +
    labs(
      title = sprintf("%s (%s)", bgc_name, bgc_info$description),
      x = "Position (kb)",
      y = "Coverage"
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 11, face = "bold")
    )
}

# 4パネル結合
plots <- list()
for (i in 1:nrow(bgc_regions)) {
  bgc_name <- bgc_regions$bgc[i]
  plots[[bgc_name]] <- create_simple_coverage(
    coverage_data[[bgc_name]],
    bgc_name,
    bgc_regions[i, ]
  )
}

# NULL を除去
plots <- plots[!sapply(plots, is.null)]

if (length(plots) > 0) {
  combined_plot <- wrap_plots(plots, ncol = 2) +
    plot_annotation(
      title = "RNA-seq Coverage across Major BGC Clusters",
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )

  # 凡例を追加
  legend_plot <- ggplot(data.frame(x = 1:3, condition = c("T1", "T2", "T3")),
                        aes(x = x, y = 1, color = condition)) +
    geom_point(size = 4) +
    scale_color_manual(values = c("T1" = "#66C2A5", "T2" = "#FC8D62", "T3" = "#8DA0CB"),
                       name = "Timepoint") +
    theme_void() +
    theme(legend.position = "bottom")

  final_plot <- combined_plot / guide_area() + plot_layout(guides = "collect", heights = c(10, 1))

  ggsave(file.path(fig_dir, "coverage_tracks_all_BGC.pdf"), combined_plot, width = 14, height = 10)
  ggsave(file.path(fig_dir, "coverage_tracks_all_BGC.png"), combined_plot, width = 14, height = 10, dpi = 150)

  cat("Saved combined coverage figure\n")
}

cat("\nCoverage track analysis completed.\n")
