#!/usr/bin/env Rscript
# 08_candidate_TF_prioritization: Candidate TF scoring and prioritization
# M145 RNA-seq (3 timepoints x 3 replicates = 9 samples)

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
})

# ── Paths ──────────────────────────────────────────────────────────────────
BASE       <- "/Users/okaban/bioinfo/rna-seq"
DESEQ_RUN  <- file.path(BASE, "04_deseq2/analysis/04_deseq2_260128_v1")
ANNOT_RUN  <- file.path(BASE, "05_annotation/analysis/05_annotation_260128_v1")
BGC_RUN    <- file.path(BASE, "06_BGC_dynamics/analysis/06_BGC_dynamics_260128_v1")
REG_RUN    <- file.path(BASE, "07_regulator_network/analysis/07_regulator_network_260128_v1")
OUT_DIR    <- file.path(BASE, "08_candidate_TF_prioritization/analysis/08_candidate_TF_260128_v1")

LOG_FILE   <- file.path(OUT_DIR, "logs/TF_candidate_pipeline.log")
sink(LOG_FILE, split = TRUE)

cat("=== 08_candidate_TF_prioritization started at", format(Sys.time()), "===\n")

# ── 1. Load input data ────────────────────────────────────────────────────
cat("\n--- Loading input data ---\n")

# 1a. Regulator-BGC correlation table (856 regulators)
reg_corr <- read_tsv(file.path(REG_RUN, "tables/regulator_BGC_correlation.tsv"),
                     show_col_types = FALSE)
cat("Regulator correlation table:", nrow(reg_corr), "regulators\n")

# 1b. Known vs novel table (13 known regulators)
reg_known <- read_tsv(file.path(REG_RUN, "tables/regulator_BGC_known_vs_novel.tsv"),
                      show_col_types = FALSE)
cat("Known regulators:", nrow(reg_known), "\n")

# 1c. Top hits by BGC
reg_top <- read_tsv(file.path(REG_RUN, "tables/regulator_top_hits_by_BGC.tsv"),
                    show_col_types = FALSE)
cat("Top hits:", nrow(reg_top), "entries\n")

# 1d. DESeq2 M145_3_vs_1 results
de_3v1 <- read_tsv(file.path(DESEQ_RUN, "results/DESeq2_M145_3_vs_1.tsv"),
                   show_col_types = FALSE) %>%
  select(gene_id, log2FC_3_vs_1 = log2FoldChange, padj_3_vs_1 = padj)
cat("DESeq2 3_vs_1:", nrow(de_3v1), "genes\n")

# 1e. rlog data for phase-specific delta computation
rld <- readRDS(file.path(DESEQ_RUN, "rds/rld.rds"))
rlog_mat <- assay(rld)
cat("rlog matrix:", nrow(rlog_mat), "genes x", ncol(rlog_mat), "samples\n")

# 1f. BGC sample scores
bgc_scores <- read_tsv(file.path(BGC_RUN, "tables/BGC_sample_scores.tsv"),
                       show_col_types = FALSE)
cat("BGC scores:", nrow(bgc_scores), "samples\n")

# ── 2. Compute phase-specific deltas for regulators ───────────────────────
cat("\n--- Computing phase-specific deltas ---\n")

# Map samples to conditions
sample_cond <- data.frame(
  sample_id = colnames(rlog_mat),
  condition = sub("_[0-9]$", "", colnames(rlog_mat)),
  stringsAsFactors = FALSE
)
cat("Sample-condition mapping:\n")
print(sample_cond)

# Compute per-condition means for each regulator
reg_genes <- reg_corr$gene_id
rlog_reg  <- rlog_mat[rownames(rlog_mat) %in% reg_genes, ]
cat("Regulator rlog matrix:", nrow(rlog_reg), "x", ncol(rlog_reg), "\n")

# Mean rlog per condition for each regulator
cond_means <- data.frame(gene_id = rownames(rlog_reg))
for (cond in c("M145_1", "M145_2", "M145_3")) {
  samps <- sample_cond$sample_id[sample_cond$condition == cond]
  cond_means[[paste0("mean_rlog_", sub("M145_", "", cond))]] <-
    rowMeans(rlog_reg[, samps, drop = FALSE])
}
cond_means <- cond_means %>%
  mutate(
    delta_2_vs_1 = mean_rlog_2 - mean_rlog_1,
    delta_3_vs_2 = mean_rlog_3 - mean_rlog_2,
    delta_3_vs_1 = mean_rlog_3 - mean_rlog_1
  )
cat("Phase deltas computed for", nrow(cond_means), "regulators\n")

# Phase preference assignment
cond_means <- cond_means %>%
  mutate(
    abs_d21 = abs(delta_2_vs_1),
    abs_d32 = abs(delta_3_vs_2),
    phase_preference = case_when(
      abs_d21 < 0.5 & abs_d32 < 0.5 ~ "stable",
      abs_d21 > abs_d32 * 1.5       ~ "early_mid",
      abs_d32 > abs_d21 * 1.5       ~ "late",
      TRUE                           ~ "both"
    )
  ) %>%
  select(-abs_d21, -abs_d32)

cat("Phase preference distribution:\n")
print(table(cond_means$phase_preference))

# ── 3. BGC score phase deltas ─────────────────────────────────────────────
cat("\n--- BGC score phase deltas ---\n")
bgc_means <- bgc_scores %>%
  mutate(condition = sub("_[0-9]$", "", sample_id)) %>%
  group_by(condition) %>%
  summarise(across(ends_with("_score"), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  arrange(condition)

bgc_phase <- data.frame(
  bgc = c("act", "red", "cda", "cpk"),
  score_1 = as.numeric(bgc_means[bgc_means$condition == "M145_1",
                                  c("act_score", "red_score", "cda_score", "cpk_score")]),
  score_2 = as.numeric(bgc_means[bgc_means$condition == "M145_2",
                                  c("act_score", "red_score", "cda_score", "cpk_score")]),
  score_3 = as.numeric(bgc_means[bgc_means$condition == "M145_3",
                                  c("act_score", "red_score", "cda_score", "cpk_score")])
)
bgc_phase <- bgc_phase %>%
  mutate(
    delta_score_2_vs_1 = score_2 - score_1,
    delta_score_3_vs_2 = score_3 - score_2,
    phase_pattern = case_when(
      delta_score_3_vs_2 > delta_score_2_vs_1 * 2 ~ "late_dominant",
      delta_score_2_vs_1 > delta_score_3_vs_2 * 2 ~ "early_mid_dominant",
      TRUE ~ "gradual"
    )
  )
cat("BGC phase patterns:\n")
print(bgc_phase %>% select(bgc, delta_score_2_vs_1, delta_score_3_vs_2, phase_pattern))

# ── 4. Build master prioritization table ──────────────────────────────────
cat("\n--- Building master prioritization table ---\n")

# Start from reg_corr, add known info, DE results, phase deltas
master <- reg_corr %>%
  left_join(de_3v1, by = "gene_id") %>%
  left_join(cond_means %>% select(gene_id, mean_rlog_1, mean_rlog_2, mean_rlog_3,
                                   delta_2_vs_1, delta_3_vs_2, delta_3_vs_1,
                                   phase_preference),
            by = "gene_id")

# Ensure abs_corr columns exist (they should from reg_corr)
master <- master %>%
  mutate(
    abs_log2FC_3_vs_1 = abs(log2FC_3_vs_1),
    max_abs_corr = pmax(abs_corr_act, abs_corr_red, abs_corr_cda, abs_corr_cpk, na.rm = TRUE),
    # Act-specificity: how much stronger is act correlation vs the best of the others
    act_specificity = abs_corr_act - pmax(abs_corr_red, abs_corr_cda, abs_corr_cpk, na.rm = TRUE),
    # Mid-BGC score: average of red/cda/cpk correlation
    mid_mean_corr = (abs_corr_red + abs_corr_cda + abs_corr_cpk) / 3
  )

cat("Master table:", nrow(master), "regulators x", ncol(master), "columns\n")

# Save master table
write_tsv(master, file.path(OUT_DIR, "tables/regulator_master_for_prioritization.tsv"))
cat("Saved: regulator_master_for_prioritization.tsv\n")

# ── 5. Scoring: TF_score_act (act late-phase candidates) ─────────────────
cat("\n--- Scoring: act late-phase candidates ---\n")

# Act late-phase scoring logic:
# - High |corr_act| → regulator covaries with act BGC score
# - High act_specificity → uniquely associated with act (not red/cda/cpk)
# - Large delta_3_vs_2 → regulator upregulated in late phase (matching act)
# - abs_log2FC_3_vs_1 → overall dynamic range
# Phase weight: late or both preferred

master_act <- master %>%
  mutate(
    phase_weight_act = case_when(
      phase_preference == "late" ~ 0.3,
      phase_preference == "both" ~ 0.15,
      TRUE ~ 0
    ),
    TF_score_act = abs_corr_act +
      0.5 * pmax(act_specificity, 0) +
      0.3 * pmin(abs_log2FC_3_vs_1 / 5, 1) +
      phase_weight_act
  ) %>%
  arrange(desc(TF_score_act))

# Top act candidates (top 30)
act_candidates <- master_act %>%
  filter(abs_corr_act >= 0.7 | TF_score_act >= 1.0) %>%
  head(30) %>%
  mutate(bgc_class = "act_late") %>%
  select(bgc_class, gene_id, old_locus_tag, gene_name, product, is_known, known_name,
         TF_score_act, corr_act, abs_corr_act, act_specificity,
         abs_corr_red, abs_corr_cda, abs_corr_cpk,
         log2FC_3_vs_1, padj_3_vs_1, delta_2_vs_1, delta_3_vs_2,
         phase_preference, mean_rlog_1, mean_rlog_2, mean_rlog_3)

cat("Act late-phase candidates:", nrow(act_candidates), "\n")
cat("Top 5 act candidates:\n")
print(act_candidates %>% select(gene_id, old_locus_tag, gene_name, known_name,
                                 TF_score_act, corr_act, act_specificity) %>% head(5))

write_tsv(act_candidates, file.path(OUT_DIR, "tables/TF_candidates_act_late.tsv"))
cat("Saved: TF_candidates_act_late.tsv\n")

# ── 6. Scoring: TF_score_mid (red/cda/cpk mid-phase candidates) ──────────
cat("\n--- Scoring: red/cda/cpk mid-phase candidates ---\n")

# Mid-phase scoring logic:
# - High mean of |corr_red|, |corr_cda|, |corr_cpk|
# - Large |delta_2_vs_1| → regulator changes early in culture
# - Phase weight: early_mid or both preferred
# - abs_log2FC_3_vs_1 → overall dynamic range

master_mid <- master %>%
  mutate(
    phase_weight_mid = case_when(
      phase_preference == "early_mid" ~ 0.3,
      phase_preference == "both" ~ 0.15,
      TRUE ~ 0
    ),
    TF_score_mid = mid_mean_corr +
      0.3 * pmin(abs_log2FC_3_vs_1 / 5, 1) +
      phase_weight_mid
  ) %>%
  arrange(desc(TF_score_mid))

# Top mid candidates (top 30)
mid_candidates <- master_mid %>%
  filter(mid_mean_corr >= 0.7 | TF_score_mid >= 1.0) %>%
  head(30) %>%
  mutate(bgc_class = "red_cda_cpk_mid") %>%
  select(bgc_class, gene_id, old_locus_tag, gene_name, product, is_known, known_name,
         TF_score_mid, corr_red, corr_cda, corr_cpk,
         abs_corr_red, abs_corr_cda, abs_corr_cpk, mid_mean_corr,
         corr_act, abs_corr_act,
         log2FC_3_vs_1, padj_3_vs_1, delta_2_vs_1, delta_3_vs_2,
         phase_preference, mean_rlog_1, mean_rlog_2, mean_rlog_3)

cat("Mid-phase candidates:", nrow(mid_candidates), "\n")
cat("Top 5 mid candidates:\n")
print(mid_candidates %>% select(gene_id, old_locus_tag, gene_name, known_name,
                                 TF_score_mid, mid_mean_corr, corr_red, corr_cda, corr_cpk) %>% head(5))

write_tsv(mid_candidates, file.path(OUT_DIR, "tables/TF_candidates_red_cda_cpk_mid.tsv"))
cat("Saved: TF_candidates_red_cda_cpk_mid.tsv\n")

# ── 7. Global regulator candidates ───────────────────────────────────────
cat("\n--- Global regulator candidates ---\n")

# Global = high |corr| with >= 3 BGCs (or >= 2 with very high values)
global_candidates <- master %>%
  mutate(
    n_high_090 = (abs_corr_act > 0.9) + (abs_corr_red > 0.9) +
      (abs_corr_cda > 0.9) + (abs_corr_cpk > 0.9),
    n_high_080 = (abs_corr_act > 0.8) + (abs_corr_red > 0.8) +
      (abs_corr_cda > 0.8) + (abs_corr_cpk > 0.8),
    # Determine sign consistency: are all high correlations same sign?
    sign_act = sign(corr_act),
    sign_red = sign(corr_red),
    sign_cda = sign(corr_cda),
    sign_cpk = sign(corr_cpk),
    all_positive = (corr_act > 0) & (corr_red > 0) & (corr_cda > 0) & (corr_cpk > 0),
    all_negative = (corr_act < 0) & (corr_red < 0) & (corr_cda < 0) & (corr_cpk < 0),
    global_role = case_when(
      all_positive & n_high_080 >= 3 ~ "global_activator",
      all_negative & n_high_080 >= 3 ~ "global_repressor",
      n_high_080 >= 3               ~ "global_mixed",
      TRUE                          ~ "not_global"
    )
  ) %>%
  filter(n_high_090 >= 2 | (n_high_080 >= 3 & mean_abs_corr >= 0.8)) %>%
  arrange(desc(mean_abs_corr)) %>%
  select(global_role, gene_id, old_locus_tag, gene_name, product, is_known, known_name,
         corr_act, corr_red, corr_cda, corr_cpk,
         abs_corr_act, abs_corr_red, abs_corr_cda, abs_corr_cpk,
         n_high_corr, mean_abs_corr, n_high_090, n_high_080,
         log2FC_3_vs_1, padj_3_vs_1, delta_2_vs_1, delta_3_vs_2,
         phase_preference)

cat("Global regulator candidates:", nrow(global_candidates), "\n")
cat("By role:\n")
print(table(global_candidates$global_role))
cat("\nTop 10 global candidates:\n")
print(global_candidates %>%
        select(global_role, gene_id, old_locus_tag, gene_name, known_name,
               mean_abs_corr, n_high_090) %>% head(10))

write_tsv(global_candidates, file.path(OUT_DIR, "tables/TF_candidates_global.tsv"))
cat("Saved: TF_candidates_global.tsv\n")

# ── 8. Experimental priority list (5-10 candidates) ──────────────────────
cat("\n--- Experimental priority candidates ---\n")

# Strategy: Select diverse high-impact candidates
# - At least 1 act-specific (actII-orf4 as positive control + novel)
# - Top novel mid-phase candidates
# - Top global activator + repressor
# - Ensure mix of known (validation) and novel (discovery)

# 8a. Act late-specific: top novel + actII-orf4 as reference
act_top <- master_act %>%
  filter(abs_corr_act >= 0.9, act_specificity > 0) %>%
  head(10)

# 8b. Mid-phase: top novel
mid_top <- master_mid %>%
  filter(mid_mean_corr >= 0.9, is_known == FALSE) %>%
  head(10)

# 8c. Global activators
global_act_top <- global_candidates %>%
  filter(global_role == "global_activator") %>%
  head(5)

# 8d. Global repressors
global_rep_top <- global_candidates %>%
  filter(global_role == "global_repressor") %>%
  head(5)

# Build priority list
build_priority <- function(df, target_bgc, role, comment_fn) {
  if (nrow(df) == 0) return(tibble())
  df %>%
    head(3) %>%
    mutate(
      target_BGC = target_bgc,
      role = role,
      comment = comment_fn(.)
    ) %>%
    select(gene_id, old_locus_tag, gene_name, product, target_BGC, role,
           is_known, known_name, comment)
}

# Manually curate the priority list
priority <- bind_rows(
  # 1. actII-orf4 (known, positive control)
  master %>%
    filter(gene_id == "SC_RS27570") %>%
    mutate(target_BGC = "act",
           role = "CSR_activator",
           comment = paste0("Known act CSR; corr_act=", round(corr_act, 3),
                            ", highly act-specific (act_specificity=",
                            round(abs_corr_act - pmax(abs_corr_red, abs_corr_cda, abs_corr_cpk), 3), ")")),

  # 2. Top novel act-late candidate
  act_top %>%
    filter(is_known == FALSE) %>%
    head(2) %>%
    mutate(target_BGC = "act",
           role = "candidate_act_activator",
           comment = paste0("Novel; corr_act=", round(corr_act, 3),
                            ", act_specificity=", round(act_specificity, 3),
                            ", phase=", phase_preference)),

  # 3. redZ (known global, positive control)
  master %>%
    filter(gene_id == "SC_RS31650") %>%
    mutate(target_BGC = "red/cda/cpk",
           role = "global_activator",
           comment = paste0("Known red CSR; corr_red=", round(corr_red, 3),
                            ", corr_cda=", round(corr_cda, 3),
                            ", corr_cpk=", round(corr_cpk, 3))),

  # 4. SC_RS21215 (SCO3818, novel global positive regulator)
  master %>%
    filter(gene_id == "SC_RS21215") %>%
    mutate(target_BGC = "red/cda/cpk",
           role = "candidate_global_activator",
           comment = paste0("Novel response regulator; 3 BGC positive corr >=0.93; ",
                            "phase=", phase_preference)),

  # 5. SC_RS29775 (SCO5518, novel global negative regulator)
  master %>%
    filter(gene_id == "SC_RS29775") %>%
    mutate(target_BGC = "red/cda/cpk",
           role = "candidate_global_repressor",
           comment = paste0("Novel PucR family; 3 BGC negative corr; KO may enhance secondary metabolism; ",
                            "phase=", phase_preference)),

  # 6. SC_RS28860 (SCO5337, XRE family, 3 BGC negative)
  master %>%
    filter(gene_id == "SC_RS28860") %>%
    mutate(target_BGC = "red/cda/cpk",
           role = "candidate_global_repressor",
           comment = paste0("Novel XRE family; 3 BGC negative corr; ",
                            "phase=", phase_preference)),

  # 7. Top novel mid-phase activator (not already selected)
  mid_top %>%
    filter(!gene_id %in% c("SC_RS21215", "SC_RS29775", "SC_RS28860")) %>%
    head(1) %>%
    mutate(target_BGC = "red/cda/cpk",
           role = "candidate_mid_activator",
           comment = paste0("Novel mid-phase; mid_mean_corr=",
                            round((abs_corr_red + abs_corr_cda + abs_corr_cpk)/3, 3),
                            "; phase=", phase_preference)),

  # 8. Top novel global activator (not already selected)
  global_act_top %>%
    filter(!gene_id %in% c("SC_RS21215", "SC_RS31650")) %>%
    head(1) %>%
    mutate(target_BGC = "act/red/cda/cpk",
           role = "candidate_global_activator",
           comment = paste0("Novel; mean_abs_corr=", round(mean_abs_corr, 3),
                            "; n_high_090=", n_high_090,
                            "; phase=", phase_preference))
) %>%
  select(gene_id, old_locus_tag, gene_name, product, target_BGC, role,
         is_known, known_name, comment) %>%
  distinct(gene_id, .keep_all = TRUE)

cat("Experimental priority candidates:", nrow(priority), "\n")
cat("\n")
print(priority %>% select(gene_id, old_locus_tag, gene_name, target_BGC, role, is_known))

write_tsv(priority, file.path(OUT_DIR, "tables/TF_candidates_experimental_priority.tsv"))
cat("\nSaved: TF_candidates_experimental_priority.tsv\n")

# ── 9. Summary statistics ────────────────────────────────────────────────
cat("\n--- Summary statistics ---\n")
cat("Total regulators analyzed:", nrow(master), "\n")
cat("Act late-phase candidates (|corr_act|>=0.7 or score>=1.0):", nrow(act_candidates), "\n")
cat("Mid-phase candidates (mid_mean_corr>=0.7 or score>=1.0):", nrow(mid_candidates), "\n")
cat("Global candidates:", nrow(global_candidates), "\n")
cat("  Global activators:", sum(global_candidates$global_role == "global_activator"), "\n")
cat("  Global repressors:", sum(global_candidates$global_role == "global_repressor"), "\n")
cat("  Global mixed:", sum(global_candidates$global_role == "global_mixed"), "\n")
cat("Experimental priority:", nrow(priority), "\n")

# Known regulators in each list
cat("\nKnown regulators in act candidates:\n")
print(act_candidates %>% filter(is_known == TRUE) %>%
        select(gene_id, known_name, TF_score_act, corr_act))
cat("\nKnown regulators in mid candidates:\n")
print(mid_candidates %>% filter(is_known == TRUE) %>%
        select(gene_id, known_name, TF_score_mid, mid_mean_corr))
cat("\nKnown regulators in global candidates:\n")
print(global_candidates %>% filter(is_known == TRUE) %>%
        select(gene_id, known_name, global_role, mean_abs_corr))

# ── 10. Figures (optional bonus) ─────────────────────────────────────────
cat("\n--- Creating figures ---\n")

library(svglite)

# 10a. Score distribution: act vs mid
fig_dir <- file.path(OUT_DIR, "figures")

# Volcano-style plot: act specificity vs act correlation
p_act <- master %>%
  mutate(highlight = case_when(
    is_known == TRUE & abs_corr_act > 0.5 ~ known_name,
    abs_corr_act > 0.95 ~ old_locus_tag,
    TRUE ~ NA_character_
  )) %>%
  ggplot(aes(x = act_specificity, y = abs_corr_act)) +
  geom_point(aes(color = phase_preference), alpha = 0.4, size = 1.5) +
  geom_point(data = . %>% filter(!is.na(highlight)),
             aes(color = phase_preference), size = 2.5) +
  ggrepel::geom_text_repel(aes(label = highlight), size = 2.5,
                            max.overlaps = 20, na.rm = TRUE) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("early_mid" = "#4DAF4A", "late" = "#E41A1C",
                                 "both" = "#377EB8", "stable" = "grey70")) +
  labs(title = "Act-specific TF candidates",
       x = "Act specificity (|corr_act| - max other |corr|)",
       y = "|corr_act|",
       color = "Phase preference") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "scatter_act_specificity_M145.pdf"), p_act,
       width = 7, height = 6)
ggsave(file.path(fig_dir, "scatter_act_specificity_M145.svg"), p_act,
       width = 7, height = 6)
cat("Saved: scatter_act_specificity_M145.pdf/.svg\n")

# 10b. Mid-phase candidates: mean mid corr vs delta_2_vs_1
p_mid <- master %>%
  mutate(highlight = case_when(
    is_known == TRUE & mid_mean_corr > 0.8 ~ known_name,
    mid_mean_corr > 0.95 ~ old_locus_tag,
    TRUE ~ NA_character_
  )) %>%
  ggplot(aes(x = delta_2_vs_1, y = mid_mean_corr)) +
  geom_point(aes(color = phase_preference), alpha = 0.4, size = 1.5) +
  geom_point(data = . %>% filter(!is.na(highlight)),
             aes(color = phase_preference), size = 2.5) +
  ggrepel::geom_text_repel(aes(label = highlight), size = 2.5,
                            max.overlaps = 20, na.rm = TRUE) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("early_mid" = "#4DAF4A", "late" = "#E41A1C",
                                 "both" = "#377EB8", "stable" = "grey70")) +
  labs(title = "Mid-phase TF candidates (red/cda/cpk)",
       x = "Delta rlog (M145_2 - M145_1)",
       y = "Mean |corr| with red/cda/cpk",
       color = "Phase preference") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "scatter_mid_candidates_M145.pdf"), p_mid,
       width = 7, height = 6)
ggsave(file.path(fig_dir, "scatter_mid_candidates_M145.svg"), p_mid,
       width = 7, height = 6)
cat("Saved: scatter_mid_candidates_M145.pdf/.svg\n")

# 10c. Global candidates: heatmap-style dot plot
global_for_plot <- global_candidates %>%
  head(30) %>%
  mutate(label = ifelse(!is.na(known_name), known_name,
                        ifelse(!is.na(gene_name), gene_name, old_locus_tag))) %>%
  select(label, global_role, corr_act, corr_red, corr_cda, corr_cpk) %>%
  pivot_longer(cols = starts_with("corr_"), names_to = "bgc", values_to = "correlation") %>%
  mutate(bgc = sub("corr_", "", bgc),
         bgc = factor(bgc, levels = c("act", "red", "cda", "cpk")))

p_global <- global_for_plot %>%
  ggplot(aes(x = bgc, y = fct_rev(fct_inorder(label)),
             fill = correlation, size = abs(correlation))) +
  geom_point(shape = 21, stroke = 0.3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                        midpoint = 0, limits = c(-1, 1)) +
  scale_size_continuous(range = c(1, 5), guide = "none") +
  labs(title = "Top 30 global regulator candidates",
       x = "BGC", y = "", fill = "Pearson r") +
  theme_bw(base_size = 9) +
  theme(axis.text.y = element_text(size = 7))

ggsave(file.path(fig_dir, "dotplot_global_candidates_M145.pdf"), p_global,
       width = 6, height = 8)
ggsave(file.path(fig_dir, "dotplot_global_candidates_M145.svg"), p_global,
       width = 6, height = 8)
cat("Saved: dotplot_global_candidates_M145.pdf/.svg\n")

# 10d. Priority candidates: summary barplot
priority_for_plot <- priority %>%
  left_join(master %>% select(gene_id, corr_act, corr_red, corr_cda, corr_cpk),
            by = "gene_id") %>%
  mutate(label = ifelse(!is.na(known_name), known_name,
                        ifelse(!is.na(gene_name), gene_name, old_locus_tag))) %>%
  select(label, role, corr_act, corr_red, corr_cda, corr_cpk) %>%
  pivot_longer(cols = starts_with("corr_"), names_to = "bgc", values_to = "correlation") %>%
  mutate(bgc = sub("corr_", "", bgc),
         bgc = factor(bgc, levels = c("act", "red", "cda", "cpk")))

p_priority <- priority_for_plot %>%
  ggplot(aes(x = bgc, y = correlation, fill = bgc)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, color = "grey50") +
  facet_wrap(~label, nrow = 2, scales = "free_x") +
  scale_fill_manual(values = c("act" = "#E41A1C", "red" = "#377EB8",
                                "cda" = "#4DAF4A", "cpk" = "#984EA3")) +
  labs(title = "Experimental priority candidates: BGC correlations",
       x = "", y = "Pearson r", fill = "BGC") +
  theme_bw(base_size = 9) +
  coord_cartesian(ylim = c(-1, 1)) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "barplot_priority_candidates_M145.pdf"), p_priority,
       width = 10, height = 5)
ggsave(file.path(fig_dir, "barplot_priority_candidates_M145.svg"), p_priority,
       width = 10, height = 5)
cat("Saved: barplot_priority_candidates_M145.pdf/.svg\n")

cat("\n=== 08_candidate_TF_prioritization completed at", format(Sys.time()), "===\n")
cat("Output directory:", OUT_DIR, "\n")

sink()
