#!/usr/bin/env Rscript
# =============================================================================
# 05_report_tables_visualization.R
# Publication-quality figures for report tables
# Project: M145 RNA-seq
# Date: 2026-02-02
# =============================================================================

library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(svglite)

# --- Path settings ---
table_dir <- "/Users/okaban/bioinfo/rna-seq/12_supplementary_figures/analysis/12_supplementary_260202_v1/tables"
fig_dir <- "/Users/okaban/bioinfo/rna-seq/12_supplementary_figures/analysis/12_supplementary_260202_v1/figures"

# --- Publication theme ---
theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
  theme(
    text = element_text(family = "Helvetica", color = "black"),
    axis.text = element_text(color = "black", size = base_size - 1),
    axis.title = element_text(size = base_size, face = "bold"),
    plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = base_size - 1, hjust = 0, color = "gray30"),
    legend.title = element_text(size = base_size - 1, face = "bold"),
    legend.text = element_text(size = base_size - 2),
    legend.key.size = unit(0.5, "cm"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.8),
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(face = "bold", size = base_size - 1)
  )
}

# --- Color palette (colorblind-friendly) ---
col_up <- "#D73027"      # Red
col_down <- "#4575B4"    # Blue
col_neutral <- "#808080" # Gray
col_t2 <- "#FDB462"      # Orange
col_t3 <- "#80B1D3"      # Light blue

# =============================================================================
# 1. Pathway Change Pattern Summary (Stacked Bar Chart)
# =============================================================================
cat("=== 1. Pathway Pattern Summary ===\n")

pattern_df <- tibble(
  Pattern = c("Both significant\n(same direction)", "Direction\nreversal", "T2 only\nsignificant", "T3 only\nsignificant"),
  Up = c(3, 2, 7, 9),
  Down = c(14, 0, 0, 28)
) %>%
  pivot_longer(cols = c(Up, Down), names_to = "Direction", values_to = "Count") %>%
  mutate(
    Pattern = factor(Pattern, levels = c("Both significant\n(same direction)", "Direction\nreversal", "T2 only\nsignificant", "T3 only\nsignificant")),
    Direction = factor(Direction, levels = c("Up", "Down"))
  )

p1 <- ggplot(pattern_df, aes(x = Pattern, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = ifelse(Count > 0, Count, "")),
            position = position_stack(vjust = 0.5), size = 4.5, color = "white", fontface = "bold") +
  scale_fill_manual(values = c("Up" = col_up, "Down" = col_down),
                    labels = c("Up-regulated", "Down-regulated")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_pub(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 10)
  ) +
  labs(
    title = "KEGG Pathway Change Patterns",
    subtitle = "Comparison of T2 vs T1 and T3 vs T1 (n = 63 pathways)",
    x = NULL, y = "Number of pathways", fill = NULL
  )

ggsave(file.path(fig_dir, "report_pathway_pattern_summary.pdf"), p1, width = 8, height = 5, device = cairo_pdf)
ggsave(file.path(fig_dir, "report_pathway_pattern_summary.svg"), p1, width = 8, height = 5)
cat("Saved: report_pathway_pattern_summary.pdf/svg\n")

# =============================================================================
# 2. Fold Enrichment Change Patterns (Arrow Plot)
# =============================================================================
cat("\n=== 2. FE Change Patterns ===\n")

fe_change_df <- tibble(
  Pathway = c(
    "Metabolic pathways", "Biosynthesis of secondary metabolites",
    "Prodigiosin biosynthesis", "Siderophore biosynthesis",
    "Ribosome",
    "Biosynthesis of cofactors", "Pyrimidine metabolism",
    "Nucleotide metabolism", "DNA replication"
  ),
  T2_FE = c(1.16, 1.35, 2.83, 3.06, 3.34, 1.64, 2.39, 2.16, 2.43),
  T3_FE = c(-1.50, -1.53, 1.83, 2.60, 2.76, 2.17, 3.06, 2.79, 2.93),
  Category = c(
    "Reversal (Up to Down)", "Reversal (Up to Down)",
    "Both Up (T2 peak)", "Both Up (T2 peak)",
    "Both Down (T2 peak)",
    "Both Down (T3 stronger)", "Both Down (T3 stronger)",
    "Both Down (T3 stronger)", "Both Down (T3 stronger)"
  )
) %>%
  mutate(
    T2_FE_plot = ifelse(grepl("Reversal|Both Up", Category), T2_FE, -T2_FE),
    T3_FE_plot = T3_FE,
    Pathway = factor(Pathway, levels = rev(Pathway))
  )

p2 <- ggplot(fe_change_df) +
  geom_segment(aes(x = T2_FE_plot, xend = T3_FE_plot, y = Pathway, yend = Pathway, color = Category),
               linewidth = 1.2, arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
  geom_point(aes(x = T2_FE_plot, y = Pathway), size = 3.5, color = col_t2, shape = 16) +
  geom_point(aes(x = T3_FE_plot, y = Pathway, color = Category), size = 3.5, shape = 16) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  scale_color_manual(values = c(
    "Reversal (Up to Down)" = "#7570B3",
    "Both Up (T2 peak)" = col_up,
    "Both Down (T2 peak)" = "#1B9E77",
    "Both Down (T3 stronger)" = col_down
  )) +
  scale_x_continuous(breaks = seq(-4, 4, 1)) +
  theme_pub(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.grid.major.y = element_blank()
  ) +
  guides(color = guide_legend(nrow = 2)) +
  labs(
    title = "Fold Enrichment Changes from T2 to T3",
    subtitle = "Orange circle: T2 value; Arrow end: T3 value",
    x = "Fold enrichment (positive: up-regulated, negative: down-regulated)",
    y = NULL, color = "Pattern"
  ) +
  annotate("text", x = 2.5, y = 0.3, label = "Up-regulated", size = 3.5, fontface = "italic", color = "gray40") +
  annotate("text", x = -2.5, y = 0.3, label = "Down-regulated", size = 3.5, fontface = "italic", color = "gray40")

ggsave(file.path(fig_dir, "report_FE_change_patterns.pdf"), p2, width = 10, height = 6, device = cairo_pdf)
ggsave(file.path(fig_dir, "report_FE_change_patterns.svg"), p2, width = 10, height = 6)
cat("Saved: report_FE_change_patterns.pdf/svg\n")

# =============================================================================
# 3. BGC Expression Dynamics Heatmap
# =============================================================================
cat("\n=== 3. BGC Expression Heatmap ===\n")

bgc_df <- tibble(
  BGC = c("act (Actinorhodin)", "red (Prodigiosin)", "cda (CDA)", "cpk (Coelimycin P1)"),
  T1 = c(1, 1, 1, 1),
  T2 = c(2, 3, 2, 2),
  T3 = c(3, 3, 2, 3)
) %>%
  pivot_longer(cols = c(T1, T2, T3), names_to = "Timepoint", values_to = "Expression") %>%
  mutate(
    BGC = factor(BGC, levels = c("act (Actinorhodin)", "red (Prodigiosin)", "cda (CDA)", "cpk (Coelimycin P1)")),
    Timepoint = factor(Timepoint, levels = c("T1", "T2", "T3")),
    Expression_label = case_when(
      Expression == 1 ~ "Low",
      Expression == 2 ~ "Mid",
      Expression == 3 ~ "High"
    )
  )

p3 <- ggplot(bgc_df, aes(x = Timepoint, y = BGC, fill = Expression)) +
  geom_tile(color = "white", linewidth = 1.5) +
  geom_text(aes(label = Expression_label), size = 4, fontface = "bold", color = "black") +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 2, limits = c(1, 3),
                       breaks = c(1, 2, 3), labels = c("Low", "Mid", "High")) +
  scale_x_discrete(labels = c("T1\n(Early)", "T2\n(Mid)", "T3\n(Late)")) +
  theme_pub(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Biosynthetic Gene Cluster Expression Dynamics",
    subtitle = "Expression levels across growth phases in S. coelicolor M145",
    x = "Growth phase", y = NULL, fill = "Expression\nlevel"
  )

ggsave(file.path(fig_dir, "report_BGC_expression_heatmap.pdf"), p3, width = 7, height = 4.5, device = cairo_pdf)
ggsave(file.path(fig_dir, "report_BGC_expression_heatmap.svg"), p3, width = 7, height = 4.5)
cat("Saved: report_BGC_expression_heatmap.pdf/svg\n")

# =============================================================================
# 4. COG Functional Category Net Changes (Horizontal Bar Chart)
# =============================================================================
cat("\n=== 4. COG Classification ===\n")

# --- 4a. T3 vs T1 (19 categories) ---
cog_t3_df <- read_tsv(file.path(table_dir, "report_COG_T3vsT1.tsv"), show_col_types = FALSE) %>%
  mutate(
    Category_short = str_extract(Category, "^[A-Z]"),
    Category_name = str_extract(Category, "(?<= - ).*"),
    Net = as.numeric(gsub("[+]", "", Net))
  ) %>%
  arrange(Net)

cog_t3_df$Category_label <- paste0(cog_t3_df$Category_short, ": ", cog_t3_df$Category_name)
cog_t3_df$Category_label <- factor(cog_t3_df$Category_label, levels = cog_t3_df$Category_label)

p4a <- ggplot(cog_t3_df, aes(x = Net, y = Category_label, fill = Net > 0)) +
  geom_bar(stat = "identity", width = 0.75, color = "black", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
  geom_text(aes(label = ifelse(Net >= 0, paste0("+", Net), Net),
                hjust = ifelse(Net >= 0, -0.15, 1.15)), size = 3.2) +
  scale_fill_manual(values = c("TRUE" = col_up, "FALSE" = col_down),
                    labels = c("Net decrease", "Net increase")) +
  scale_x_continuous(limits = c(min(cog_t3_df$Net) - 30, max(cog_t3_df$Net) + 80),
                     breaks = seq(-100, 450, 100)) +
  theme_pub(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank()
  ) +
  labs(
    title = "COG Functional Category Changes (T3 vs T1)",
    subtitle = paste0("Net change = (up-regulated) - (down-regulated); n = ", nrow(cog_t3_df), " categories"),
    x = "Net change in gene count", y = NULL
  )

ggsave(file.path(fig_dir, "report_COG_net_change.pdf"), p4a, width = 10, height = 8, device = cairo_pdf)
ggsave(file.path(fig_dir, "report_COG_net_change.svg"), p4a, width = 10, height = 8)
cat("Saved: report_COG_net_change.pdf/svg (T3 vs T1, 19 categories)\n")

# --- 4b. T2 vs T1 (20 categories) ---
cog_t2_df <- read_tsv(file.path(table_dir, "report_COG_T2vsT1.tsv"), show_col_types = FALSE) %>%
  mutate(
    Category_short = str_extract(Category, "^[A-Z]"),
    Category_name = str_extract(Category, "(?<= - ).*"),
    Net = as.numeric(gsub("[+]", "", Net))
  ) %>%
  arrange(Net)

cog_t2_df$Category_label <- paste0(cog_t2_df$Category_short, ": ", cog_t2_df$Category_name)
cog_t2_df$Category_label <- factor(cog_t2_df$Category_label, levels = cog_t2_df$Category_label)

p4b <- ggplot(cog_t2_df, aes(x = Net, y = Category_label, fill = Net > 0)) +
  geom_bar(stat = "identity", width = 0.75, color = "black", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
  geom_text(aes(label = ifelse(Net >= 0, paste0("+", Net), Net),
                hjust = ifelse(Net >= 0, -0.15, 1.15)), size = 3.2) +
  scale_fill_manual(values = c("TRUE" = col_up, "FALSE" = col_down),
                    labels = c("Net decrease", "Net increase")) +
  scale_x_continuous(limits = c(min(cog_t2_df$Net) - 20, max(cog_t2_df$Net) + 50),
                     breaks = seq(-100, 250, 50)) +
  theme_pub(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank()
  ) +
  labs(
    title = "COG Functional Category Changes (T2 vs T1)",
    subtitle = paste0("Net change = (up-regulated) - (down-regulated); n = ", nrow(cog_t2_df), " categories"),
    x = "Net change in gene count", y = NULL
  )

ggsave(file.path(fig_dir, "report_COG_net_change_T2vsT1.pdf"), p4b, width = 10, height = 8, device = cairo_pdf)
ggsave(file.path(fig_dir, "report_COG_net_change_T2vsT1.svg"), p4b, width = 10, height = 8)
cat("Saved: report_COG_net_change_T2vsT1.pdf/svg (T2 vs T1, 20 categories)\n")

# =============================================================================
# 5. Persistent Down-regulated Pathways FE Comparison (Dot Plot)
# =============================================================================
cat("\n=== 5. Persistent Down-regulated Pathways ===\n")

down_df <- read_tsv(file.path(table_dir, "report_both_down_pathways.tsv"), show_col_types = FALSE) %>%
  select(Pathway, T2vsT1_FE, T3vsT1_FE, FE_Change) %>%
  mutate(
    Pathway_short = case_when(
      Pathway == "Amino sugar/nucleotide sugar" ~ "Amino sugar/nucleotide sugar metab.",
      Pathway == "Alanine, aspartate, glutamate" ~ "Ala, Asp, Glu metabolism",
      Pathway == "Biosynthesis of nucleotide sugars" ~ "Nucleotide sugar biosynthesis",
      Pathway == "Aminoacyl-tRNA biosynthesis" ~ "Aminoacyl-tRNA biosynthesis",
      TRUE ~ Pathway
    ),
    Pathway_short = str_trunc(Pathway_short, 35),
    Change_type = case_when(
      FE_Change == "\\u2191" | FE_Change == "↑" ~ "T3 stronger",
      FE_Change == "\\u2193" | FE_Change == "↓" ~ "T2 stronger",
      TRUE ~ "Sustained"
    )
  ) %>%
  arrange(T3vsT1_FE)

down_df$Pathway_short <- factor(down_df$Pathway_short, levels = down_df$Pathway_short)

down_long <- down_df %>%
  pivot_longer(cols = c(T2vsT1_FE, T3vsT1_FE), names_to = "Comparison", values_to = "FE") %>%
  mutate(Comparison = ifelse(Comparison == "T2vsT1_FE", "T2 vs T1", "T3 vs T1"))

p5 <- ggplot(down_long, aes(x = FE, y = Pathway_short)) +
  geom_line(aes(group = Pathway_short), color = "gray60", linewidth = 0.6) +
  geom_point(aes(color = Comparison, shape = Comparison), size = 3.5) +
  scale_color_manual(values = c("T2 vs T1" = col_t2, "T3 vs T1" = col_down)) +
  scale_shape_manual(values = c("T2 vs T1" = 16, "T3 vs T1" = 17)) +
  scale_x_continuous(limits = c(0, 4), breaks = seq(0, 4, 0.5)) +
  theme_pub(base_size = 11) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_blank()
  ) +
  labs(
    title = "Persistent Down-regulated Pathways",
    subtitle = "Fold enrichment comparison between T2 vs T1 and T3 vs T1 (n = 14 pathways)",
    x = "Fold enrichment", y = NULL, color = NULL, shape = NULL
  )

ggsave(file.path(fig_dir, "report_both_down_FE_comparison.pdf"), p5, width = 9, height = 7, device = cairo_pdf)
ggsave(file.path(fig_dir, "report_both_down_FE_comparison.svg"), p5, width = 9, height = 7)
cat("Saved: report_both_down_FE_comparison.pdf/svg\n")

# =============================================================================
# 6. Key Biological Insights Summary (Horizontal Bar)
# =============================================================================
cat("\n=== 6. Key Insights Summary ===\n")

insights_df <- tibble(
  ID = 1:6,
  Insight = c(
    "T2 is the peak of transcriptional reprogramming",
    "Prodigiosin shows early induction (FE: 2.83 at T2)",
    "~1,000 core response genes across all timepoints",
    "Clear growth-secondary metabolism tradeoff",
    "Distinct BGC activation timing (red > cda > act/cpk)",
    "Nutrient starvation adaptation response at T3"
  ),
  Category = c("Dynamics", "BGC", "Core genes", "Tradeoff", "BGC", "Adaptation"),
  Evidence = c("Ribosome FE: 3.34 (T2) > 2.76 (T3)",
               "KEGG enrichment",
               "Venn diagram overlap",
               "COG: Translation J = -98, Secondary Q = +12",
               "Coverage tracks",
               "ABC transporters: 108 genes")
)

p6 <- ggplot(insights_df, aes(x = reorder(Insight, ID), y = ID)) +
  geom_segment(aes(xend = Insight, y = 0, yend = ID, color = Category), linewidth = 1.5) +
  geom_point(aes(color = Category), size = 4) +
  geom_text(aes(label = paste0("#", ID)), y = -0.3, size = 3.5, fontface = "bold") +
  coord_flip() +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(limits = c(-1, 7), breaks = NULL) +
  theme_pub(base_size = 11) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.text.y = element_text(size = 10)
  ) +
  labs(
    title = "Key Biological Insights from Functional Enrichment Analysis",
    subtitle = "S. coelicolor M145 transcriptome dynamics across growth phases",
    x = NULL, y = NULL, color = "Category"
  )

ggsave(file.path(fig_dir, "report_key_insights_summary.pdf"), p6, width = 11, height = 5, device = cairo_pdf)
ggsave(file.path(fig_dir, "report_key_insights_summary.svg"), p6, width = 11, height = 5)
cat("Saved: report_key_insights_summary.pdf/svg\n")

cat("\n=== All publication-quality figures completed ===\n")
