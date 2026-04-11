#!/usr/bin/env Rscript
# 10_SARP_motif_scan pipeline — M145
# Steps 5-7: FIMO scan, aggregation, and reporting

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

cat("=== 10_SARP_motif_scan Steps 5-7 started at", format(Sys.time()), "===\n")

# -------------------------------------------------------------------
# 0. Paths
# -------------------------------------------------------------------
MOTIF_ROOT   <- "/Users/okaban/bioinfo/rna-seq/10_SARP_motif_scan"
ANNOT_ROOT   <- "/Users/okaban/bioinfo/rna-seq/05_annotation"
TF_ROOT      <- "/Users/okaban/bioinfo/rna-seq/08_candidate_TF_prioritization"
SARP_ROOT    <- "/Users/okaban/bioinfo/rna-seq/09_SARP_integration"

ANNOT_RUN    <- file.path(ANNOT_ROOT, "analysis/05_annotation_260128_v1")
TF_RUN       <- file.path(TF_ROOT,    "analysis/08_candidate_TF_260128_v1")
SARP_V2_RUN  <- file.path(SARP_ROOT,  "analysis/09_SARP_integration_260128_v2")
RUN_DIR      <- file.path(MOTIF_ROOT, "analysis/10_SARP_motif_scan_260128_v1")

OUT_TABLES    <- file.path(RUN_DIR, "tables")
OUT_FIGURES   <- file.path(RUN_DIR, "figures")
OUT_PROMOTERS <- file.path(RUN_DIR, "promoters")
OUT_PWM       <- file.path(RUN_DIR, "pwm")

# Input files
COORD_FILE   <- file.path(OUT_TABLES, "BGC_and_TF_coordinates.tsv")
SARP_V2_LIST <- file.path(SARP_V2_RUN, "tables/SARP_list_M145_v2.tsv")
SARP_V2_ACT  <- file.path(SARP_V2_RUN, "tables/SARP_BGC_activity_summary_v2.tsv")
TF_CANDIDATES <- file.path(TF_RUN, "tables/TF_candidates_experimental_priority.tsv")
GENE_MASTER  <- file.path(ANNOT_RUN, "tables/gene_master_with_BGC_regulators.tsv")

# Promoter FASTA files
BGC_FASTA    <- file.path(OUT_PROMOTERS, "promoters_BGCs_500bp.fasta")
SARP_FASTA   <- file.path(OUT_PROMOTERS, "promoters_SARPs_500bp.fasta")
TF_FASTA     <- file.path(OUT_PROMOTERS, "promoters_TFs_500bp.fasta")

# PWM files
PWM_ACTII    <- file.path(OUT_PWM, "PWM_ActII_ORF4.meme")
PWM_STRICT   <- file.path(OUT_PWM, "PWM_SARP_heptamer_strict.meme")
PWM_RELAXED  <- file.path(OUT_PWM, "PWM_SARP_heptamer_relaxed.meme")

# -------------------------------------------------------------------
# 1. Read reference data
# -------------------------------------------------------------------
cat("\n--- Reading reference data ---\n")

coords <- read_tsv(COORD_FILE, show_col_types = FALSE)
sarp_v2 <- read_tsv(SARP_V2_LIST, show_col_types = FALSE)
gm <- read_tsv(GENE_MASTER, show_col_types = FALSE)
tf_cand <- read_tsv(TF_CANDIDATES, show_col_types = FALSE)

# Try to read activity summary
sarp_act <- tryCatch(
  read_tsv(SARP_V2_ACT, show_col_types = FALSE),
  error = function(e) NULL
)

cat("  Coordinates:", nrow(coords), "entries\n")
cat("  SARPs:", nrow(sarp_v2), "\n")
cat("  TF candidates:", nrow(tf_cand), "\n")

# -------------------------------------------------------------------
# 2. Run FIMO for all PWM × FASTA combinations
# -------------------------------------------------------------------
cat("\n--- Step 5: Running FIMO motif scans ---\n")

pwm_files <- list(
  ActII_ORF4        = PWM_ACTII,
  SARP_heptamer_strict  = PWM_STRICT,
  SARP_heptamer_relaxed = PWM_RELAXED
)

fasta_files <- list(
  BGC  = BGC_FASTA,
  SARP = SARP_FASTA,
  TF   = TF_FASTA
)

# Run FIMO for each combination
fimo_results <- list()

for (pwm_name in names(pwm_files)) {
  for (fasta_name in names(fasta_files)) {
    label <- paste0(pwm_name, "_vs_", fasta_name)
    out_dir <- file.path(RUN_DIR, "fimo_raw", label)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    cmd <- sprintf(
      "fimo --oc %s --thresh 1e-3 --verbosity 1 %s %s",
      shQuote(out_dir),
      shQuote(pwm_files[[pwm_name]]),
      shQuote(fasta_files[[fasta_name]])
    )

    cat("  Running:", label, "...\n")
    ret <- system(cmd, intern = FALSE)

    fimo_tsv <- file.path(out_dir, "fimo.tsv")
    if (file.exists(fimo_tsv)) {
      df <- tryCatch(
        read_tsv(fimo_tsv, comment = "#", show_col_types = FALSE),
        error = function(e) {
          cat("    WARNING: Could not parse", fimo_tsv, "\n")
          NULL
        }
      )
      if (!is.null(df) && nrow(df) > 0) {
        df <- df %>%
          mutate(pwm_name = pwm_name, target_set = fasta_name)
        fimo_results[[label]] <- df
        cat("    Hits:", nrow(df), "\n")
      } else {
        cat("    No hits above threshold\n")
      }
    } else {
      cat("    WARNING: FIMO output not found\n")
    }
  }
}

# Combine all FIMO results
if (length(fimo_results) > 0) {
  all_fimo <- bind_rows(fimo_results)
  cat("\n  Total FIMO hits across all scans:", nrow(all_fimo), "\n")
} else {
  cat("\n  WARNING: No FIMO hits found. Trying with relaxed threshold...\n")
  # Retry with more relaxed threshold
  for (pwm_name in names(pwm_files)) {
    for (fasta_name in names(fasta_files)) {
      label <- paste0(pwm_name, "_vs_", fasta_name)
      out_dir <- file.path(RUN_DIR, "fimo_raw", label)
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      cmd <- sprintf(
        "fimo --oc %s --thresh 0.01 --verbosity 1 %s %s",
        shQuote(out_dir),
        shQuote(pwm_files[[pwm_name]]),
        shQuote(fasta_files[[fasta_name]])
      )

      cat("  Retrying with --thresh 0.01:", label, "...\n")
      ret <- system(cmd, intern = FALSE)

      fimo_tsv <- file.path(out_dir, "fimo.tsv")
      if (file.exists(fimo_tsv)) {
        df <- tryCatch(
          read_tsv(fimo_tsv, comment = "#", show_col_types = FALSE),
          error = function(e) NULL
        )
        if (!is.null(df) && nrow(df) > 0) {
          df <- df %>%
            mutate(pwm_name = pwm_name, target_set = fasta_name)
          fimo_results[[label]] <- df
          cat("    Hits:", nrow(df), "\n")
        } else {
          cat("    Still no hits\n")
        }
      }
    }
  }
  if (length(fimo_results) > 0) {
    all_fimo <- bind_rows(fimo_results)
    cat("\n  Total FIMO hits (relaxed):", nrow(all_fimo), "\n")
  } else {
    all_fimo <- tibble()
    cat("\n  No FIMO hits found even with relaxed threshold\n")
  }
}

# -------------------------------------------------------------------
# 3. Parse sequence names → gene_id
# -------------------------------------------------------------------
cat("\n--- Parsing FIMO results ---\n")

if (nrow(all_fimo) > 0) {
  # bedtools getfasta with -name produces headers like "gene_id::chr:start-end(strand)"
  # Extract gene_id from sequence_name
  all_fimo <- all_fimo %>%
    mutate(
      gene_id = str_replace(sequence_name, "::.*$", ""),
      gene_id = str_replace(gene_id, "\\(.*\\)$", "")
    )

  # Join with coordinates to get annotation
  all_fimo <- all_fimo %>%
    left_join(
      coords %>% select(gene_id, old_locus_tag, gene_name, bgc_name, gene_type, product) %>%
        distinct(gene_id, .keep_all = TRUE),
      by = "gene_id"
    )

  cat("  Annotated FIMO hits:", nrow(all_fimo), "\n")
  cat("  Unique target genes:", n_distinct(all_fimo$gene_id), "\n")
  cat("  By target set:\n")
  print(table(all_fimo$target_set))
  cat("  By PWM:\n")
  print(table(all_fimo$pwm_name))
}

# -------------------------------------------------------------------
# 4. Write per-target-set FIMO tables
# -------------------------------------------------------------------
cat("\n--- Writing FIMO result tables ---\n")

write_fimo_table <- function(df, target, filename) {
  subset <- df %>% filter(target_set == target)
  if (nrow(subset) > 0) {
    write_tsv(subset, file.path(OUT_TABLES, filename))
    cat("  Saved:", filename, "(", nrow(subset), "hits )\n")
  } else {
    # Write empty table with headers
    write_tsv(subset, file.path(OUT_TABLES, filename))
    cat("  Saved:", filename, "( 0 hits )\n")
  }
  return(subset)
}

if (nrow(all_fimo) > 0) {
  fimo_bgc  <- write_fimo_table(all_fimo, "BGC",  "fimo_BGC_promoters.tsv")
  fimo_sarp <- write_fimo_table(all_fimo, "SARP", "fimo_SARP_promoters.tsv")
  fimo_tf   <- write_fimo_table(all_fimo, "TF",   "fimo_TF_promoters.tsv")
} else {
  # Write empty files
  empty_df <- tibble(
    motif_id = character(), motif_alt_id = character(),
    sequence_name = character(), start = integer(), stop = integer(),
    strand = character(), score = numeric(), `p-value` = numeric(),
    `q-value` = numeric(), matched_sequence = character(),
    pwm_name = character(), target_set = character(),
    gene_id = character(), old_locus_tag = character(),
    gene_name = character(), bgc_name = character(),
    gene_type = character(), product = character()
  )
  write_tsv(empty_df, file.path(OUT_TABLES, "fimo_BGC_promoters.tsv"))
  write_tsv(empty_df, file.path(OUT_TABLES, "fimo_SARP_promoters.tsv"))
  write_tsv(empty_df, file.path(OUT_TABLES, "fimo_TF_promoters.tsv"))
  fimo_bgc <- empty_df
  fimo_sarp <- empty_df
  fimo_tf <- empty_df
  cat("  Written empty FIMO tables (no hits found)\n")
}

cat("\n=== Step 5 (FIMO scan) completed ===\n")

# -------------------------------------------------------------------
# 5. Aggregate: BGC SARP motif summary
# -------------------------------------------------------------------
cat("\n--- Step 6: Aggregating results ---\n")

if (nrow(fimo_bgc) > 0) {
  bgc_summary <- fimo_bgc %>%
    group_by(bgc_name, gene_id, old_locus_tag, gene_name, product, pwm_name) %>%
    summarise(
      motif_count   = n(),
      top_score     = max(score, na.rm = TRUE),
      best_pvalue   = min(`p-value`, na.rm = TRUE),
      best_qvalue   = min(`q-value`, na.rm = TRUE),
      mean_start_pos = mean(start, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(bgc_name, best_pvalue)

  write_tsv(bgc_summary, file.path(OUT_TABLES, "BGC_SARP_motif_summary.tsv"))
  cat("  BGC_SARP_motif_summary:", nrow(bgc_summary), "rows\n")

  # Summary by BGC
  bgc_overview <- bgc_summary %>%
    group_by(bgc_name, pwm_name) %>%
    summarise(
      genes_with_hits = n_distinct(gene_id),
      total_hits      = sum(motif_count),
      best_pvalue     = min(best_pvalue),
      .groups = "drop"
    )
  cat("\n  BGC overview:\n")
  print(as.data.frame(bgc_overview))
} else {
  bgc_summary <- tibble()
  write_tsv(bgc_summary, file.path(OUT_TABLES, "BGC_SARP_motif_summary.tsv"))
  cat("  BGC_SARP_motif_summary: 0 rows (no hits)\n")
}

# -------------------------------------------------------------------
# 6. TF promoter SARP hits summary
# -------------------------------------------------------------------
if (nrow(fimo_tf) > 0) {
  tf_summary <- fimo_tf %>%
    group_by(gene_id, old_locus_tag, product, pwm_name) %>%
    summarise(
      hit_count   = n(),
      best_score  = max(score, na.rm = TRUE),
      best_pvalue = min(`p-value`, na.rm = TRUE),
      best_qvalue = min(`q-value`, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(best_pvalue)

  write_tsv(tf_summary, file.path(OUT_TABLES, "TF_promoters_SARP_hits.tsv"))
  cat("  TF_promoters_SARP_hits:", nrow(tf_summary), "rows\n")
} else {
  tf_summary <- tibble()
  write_tsv(tf_summary, file.path(OUT_TABLES, "TF_promoters_SARP_hits.tsv"))
  cat("  TF_promoters_SARP_hits: 0 rows (no hits)\n")
}

# -------------------------------------------------------------------
# 7. SARP promoter hits summary (auto-regulation / cross-regulation)
# -------------------------------------------------------------------
if (nrow(fimo_sarp) > 0) {
  sarp_promo_summary <- fimo_sarp %>%
    group_by(gene_id, old_locus_tag, product, pwm_name) %>%
    summarise(
      hit_count   = n(),
      best_score  = max(score, na.rm = TRUE),
      best_pvalue = min(`p-value`, na.rm = TRUE),
      best_qvalue = min(`q-value`, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(best_pvalue)

  cat("  SARP promoter hits:", nrow(sarp_promo_summary), "rows\n")
} else {
  sarp_promo_summary <- tibble()
  cat("  SARP promoter hits: 0 rows\n")
}

# -------------------------------------------------------------------
# 8. Build SARP target network
# -------------------------------------------------------------------
cat("\n--- Building SARP target network ---\n")

# For each SARP, identify which genes have SARP motif hits in their promoters
# We assign motifs to the most relevant SARP based on the PWM used:
# - ActII_ORF4 PWM hits → attributed to actII-orf4 (SC_RS27570)
# - SARP_heptamer PWMs → attributed to the resident CSR SARP of the BGC,
#   or "any SARP" for TF targets

build_network <- function(fimo_all) {
  if (nrow(fimo_all) == 0) return(tibble())

  network_rows <- list()

  # For BGC targets
  bgc_hits <- fimo_all %>% filter(target_set == "BGC")
  if (nrow(bgc_hits) > 0) {
    # ActII-ORF4 PWM hits in act cluster → actII-orf4
    act_bgc <- bgc_hits %>%
      filter(pwm_name == "ActII_ORF4", bgc_name == "act") %>%
      group_by(gene_id, old_locus_tag, gene_name, product, bgc_name) %>%
      summarise(
        hit_count = n(), best_score = max(score), best_pvalue = min(`p-value`),
        best_qvalue = min(`q-value`),
        .groups = "drop"
      ) %>%
      mutate(SARP_gene_id = "SC_RS27570", SARP_name = "actII-orf4",
             target_type = "BGC_gene", pwm_used = "ActII_ORF4")

    network_rows <- c(network_rows, list(act_bgc))

    # Heptamer hits in each BGC → assigned to resident SARP CSR
    bgc_sarp_map <- tribble(
      ~bgc_name, ~SARP_gene_id, ~SARP_name,
      "act", "SC_RS27585", "SCO5085",
      "red", "SC_RS31630", "redD",
      "cda", "SC_RS18200", "SCO3217/cdaR",
      "cpk", "SC_RS33690", "SCO6288/cpkN"
    )

    heptamer_bgc <- bgc_hits %>%
      filter(str_detect(pwm_name, "heptamer")) %>%
      group_by(gene_id, old_locus_tag, gene_name, product, bgc_name, pwm_name) %>%
      summarise(
        hit_count = n(), best_score = max(score), best_pvalue = min(`p-value`),
        best_qvalue = min(`q-value`),
        .groups = "drop"
      ) %>%
      left_join(bgc_sarp_map, by = "bgc_name") %>%
      mutate(target_type = "BGC_gene", pwm_used = pwm_name)

    network_rows <- c(network_rows, list(heptamer_bgc))

    # ActII-ORF4 PWM hits in non-act BGCs (cross-cluster regulation possibility)
    nonact_actii <- bgc_hits %>%
      filter(pwm_name == "ActII_ORF4", bgc_name != "act") %>%
      group_by(gene_id, old_locus_tag, gene_name, product, bgc_name) %>%
      summarise(
        hit_count = n(), best_score = max(score), best_pvalue = min(`p-value`),
        best_qvalue = min(`q-value`),
        .groups = "drop"
      )
    if (nrow(nonact_actii) > 0) {
      nonact_actii <- nonact_actii %>%
        mutate(SARP_gene_id = "SC_RS27570", SARP_name = "actII-orf4",
               target_type = "BGC_gene_cross", pwm_used = "ActII_ORF4")
      network_rows <- c(network_rows, list(nonact_actii))
    }
  }

  # For TF targets
  tf_hits <- fimo_all %>% filter(target_set == "TF")
  if (nrow(tf_hits) > 0) {
    # All SARP PWM hits on TF promoters → "SARP family"
    tf_network <- tf_hits %>%
      group_by(gene_id, old_locus_tag, product, pwm_name) %>%
      summarise(
        hit_count = n(), best_score = max(score), best_pvalue = min(`p-value`),
        best_qvalue = min(`q-value`),
        .groups = "drop"
      ) %>%
      mutate(
        bgc_name = NA_character_,
        gene_name = NA_character_,
        SARP_gene_id = case_when(
          pwm_name == "ActII_ORF4" ~ "SC_RS27570",
          TRUE ~ "SARP_family"
        ),
        SARP_name = case_when(
          pwm_name == "ActII_ORF4" ~ "actII-orf4",
          TRUE ~ "SARP_heptamer"
        ),
        target_type = "TF_candidate",
        pwm_used = pwm_name
      )
    network_rows <- c(network_rows, list(tf_network))
  }

  # For SARP auto/cross-regulation targets
  sarp_hits <- fimo_all %>% filter(target_set == "SARP")
  if (nrow(sarp_hits) > 0) {
    sarp_network <- sarp_hits %>%
      group_by(gene_id, old_locus_tag, product, pwm_name) %>%
      summarise(
        hit_count = n(), best_score = max(score), best_pvalue = min(`p-value`),
        best_qvalue = min(`q-value`),
        .groups = "drop"
      ) %>%
      mutate(
        bgc_name = NA_character_,
        gene_name = NA_character_,
        SARP_gene_id = case_when(
          pwm_name == "ActII_ORF4" ~ "SC_RS27570",
          TRUE ~ "SARP_family"
        ),
        SARP_name = case_when(
          pwm_name == "ActII_ORF4" ~ "actII-orf4",
          TRUE ~ "SARP_heptamer"
        ),
        target_type = "SARP_autoregulation",
        pwm_used = pwm_name
      )
    network_rows <- c(network_rows, list(sarp_network))
  }

  if (length(network_rows) > 0) {
    bind_rows(network_rows) %>%
      select(SARP_gene_id, SARP_name, gene_id, old_locus_tag, gene_name,
             product, bgc_name, target_type, pwm_used, hit_count,
             best_score, best_pvalue, best_qvalue) %>%
      arrange(best_pvalue)
  } else {
    tibble()
  }
}

network <- build_network(all_fimo)
write_tsv(network, file.path(OUT_TABLES, "SARP_target_network.tsv"))
cat("  SARP_target_network:", nrow(network), "rows\n")

if (nrow(network) > 0) {
  cat("  By target_type:\n")
  print(table(network$target_type))
  cat("  By SARP:\n")
  print(table(network$SARP_name))
}

cat("\n=== Step 6 (aggregation) completed ===\n")

# -------------------------------------------------------------------
# 9. Heatmap visualization
# -------------------------------------------------------------------
cat("\n--- Creating visualizations ---\n")

if (nrow(all_fimo) > 0) {
  library(ggplot2)

  # Heatmap: BGC genes × motif type
  if (nrow(fimo_bgc) > 0) {
    heatmap_data <- fimo_bgc %>%
      group_by(bgc_name, gene_id, old_locus_tag, pwm_name) %>%
      summarise(
        hit_count = n(),
        neg_log10_pval = -log10(min(`p-value`)),
        .groups = "drop"
      )

    # Order genes by genomic position within each BGC
    gene_order <- coords %>%
      filter(gene_type == "BGC_gene") %>%
      arrange(bgc_name, start) %>%
      pull(gene_id)

    heatmap_data <- heatmap_data %>%
      mutate(gene_id = factor(gene_id, levels = rev(gene_order)))

    p1 <- ggplot(heatmap_data, aes(x = pwm_name, y = gene_id, fill = neg_log10_pval)) +
      geom_tile(color = "white", linewidth = 0.2) +
      facet_grid(bgc_name ~ ., scales = "free_y", space = "free_y") +
      scale_fill_gradient(low = "lightyellow", high = "red3",
                          name = expression(-log[10](p))) +
      theme_minimal(base_size = 8) +
      theme(
        axis.text.y = element_text(size = 5),
        strip.text.y = element_text(angle = 0, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      labs(x = "Motif PWM", y = "BGC gene (promoter)",
           title = "SARP motif hits in BGC promoters")

    ggsave(file.path(OUT_FIGURES, "heatmap_BGC_SARP_motif_hits.pdf"),
           p1, width = 6, height = 12, limitsize = FALSE)
    ggsave(file.path(OUT_FIGURES, "heatmap_BGC_SARP_motif_hits.png"),
           p1, width = 6, height = 12, dpi = 150, limitsize = FALSE)
    cat("  Saved: heatmap_BGC_SARP_motif_hits.pdf/png\n")
  }

  # Heatmap: TF candidates × motif type
  if (nrow(fimo_tf) > 0) {
    tf_heat <- fimo_tf %>%
      group_by(gene_id, old_locus_tag, product, pwm_name) %>%
      summarise(
        hit_count = n(),
        neg_log10_pval = -log10(min(`p-value`)),
        .groups = "drop"
      ) %>%
      mutate(label = paste0(old_locus_tag, " (", str_trunc(product, 30), ")"))

    p2 <- ggplot(tf_heat, aes(x = pwm_name, y = label, fill = neg_log10_pval)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = hit_count), size = 3) +
      scale_fill_gradient(low = "lightyellow", high = "steelblue",
                          name = expression(-log[10](p))) +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Motif PWM", y = "TF candidate",
           title = "SARP motif hits in TF candidate promoters")

    ggsave(file.path(OUT_FIGURES, "heatmap_TF_SARP_motif_hits.pdf"),
           p2, width = 7, height = 5)
    ggsave(file.path(OUT_FIGURES, "heatmap_TF_SARP_motif_hits.png"),
           p2, width = 7, height = 5, dpi = 150)
    cat("  Saved: heatmap_TF_SARP_motif_hits.pdf/png\n")
  }

  # Overview barplot: hits per BGC per motif
  if (nrow(fimo_bgc) > 0) {
    bgc_bar <- fimo_bgc %>%
      group_by(bgc_name, pwm_name) %>%
      summarise(n_genes = n_distinct(gene_id), total_hits = n(), .groups = "drop")

    p3 <- ggplot(bgc_bar, aes(x = bgc_name, y = total_hits, fill = pwm_name)) +
      geom_col(position = "dodge") +
      theme_minimal(base_size = 12) +
      scale_fill_brewer(palette = "Set2", name = "PWM") +
      labs(x = "BGC", y = "Total FIMO hits",
           title = "SARP motif scan: hits per BGC")

    ggsave(file.path(OUT_FIGURES, "barplot_BGC_hits_by_motif.pdf"),
           p3, width = 7, height = 4)
    ggsave(file.path(OUT_FIGURES, "barplot_BGC_hits_by_motif.png"),
           p3, width = 7, height = 4, dpi = 150)
    cat("  Saved: barplot_BGC_hits_by_motif.pdf/png\n")
  }
} else {
  cat("  No FIMO hits to visualize\n")
}

cat("\n=== Visualization completed ===\n")

# -------------------------------------------------------------------
# 10. Save combined session data for report generation
# -------------------------------------------------------------------
save(all_fimo, bgc_summary, tf_summary, sarp_promo_summary, network,
     coords, sarp_v2, tf_cand, gm,
     file = file.path(RUN_DIR, "motif_scan_session.RData"))
cat("  Session data saved: motif_scan_session.RData\n")

cat("\n=== Steps 5-7 pipeline completed at", format(Sys.time()), "===\n")
