#!/usr/bin/env Rscript
# 10_SARP_motif_scan pipeline — M145
# Steps 5b-7: Re-aggregate FIMO results with proper coordinate mapping

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

cat("=== 10_SARP_motif_scan Aggregation (corrected) started at", format(Sys.time()), "===\n")

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

# Input files
COORD_FILE   <- file.path(OUT_TABLES, "BGC_and_TF_coordinates.tsv")
SARP_V2_LIST <- file.path(SARP_V2_RUN, "tables/SARP_list_M145_v2.tsv")
SARP_V2_ACT  <- file.path(SARP_V2_RUN, "tables/SARP_BGC_activity_summary_v2.tsv")
TF_CANDIDATES <- file.path(TF_RUN, "tables/TF_candidates_experimental_priority.tsv")
GENE_MASTER  <- file.path(ANNOT_RUN, "tables/gene_master_with_BGC_regulators.tsv")

# -------------------------------------------------------------------
# 1. Read reference data
# -------------------------------------------------------------------
cat("\n--- Reading reference data ---\n")

coords <- read_tsv(COORD_FILE, show_col_types = FALSE)
sarp_v2 <- read_tsv(SARP_V2_LIST, show_col_types = FALSE)
gm <- read_tsv(GENE_MASTER, show_col_types = FALSE)
tf_cand <- read_tsv(TF_CANDIDATES, show_col_types = FALSE)
sarp_act <- tryCatch(read_tsv(SARP_V2_ACT, show_col_types = FALSE), error = function(e) NULL)

cat("  Coordinates:", nrow(coords), "entries\n")

# -------------------------------------------------------------------
# 2. Read BED files to get promoter regions
# -------------------------------------------------------------------
cat("\n--- Reading promoter BED regions ---\n")

read_bed <- function(path, label) {
  bed <- read_tsv(path, col_names = c("chrom", "chromStart", "chromEnd", "name", "score", "strand"),
                  show_col_types = FALSE)
  bed$target_set <- label
  bed
}

bgc_bed  <- read_bed(file.path(OUT_PROMOTERS, "bgc_promoters.bed"), "BGC")
sarp_bed <- read_bed(file.path(OUT_PROMOTERS, "sarp_promoters.bed"), "SARP")
tf_bed   <- read_bed(file.path(OUT_PROMOTERS, "tf_promoters.bed"), "TF")

all_bed <- bind_rows(bgc_bed, sarp_bed, tf_bed)
cat("  Total promoter regions:", nrow(all_bed), "\n")

# -------------------------------------------------------------------
# 3. Read all raw FIMO outputs
# -------------------------------------------------------------------
cat("\n--- Reading raw FIMO outputs ---\n")

fimo_raw_dir <- file.path(RUN_DIR, "fimo_raw")
fimo_subdirs <- list.dirs(fimo_raw_dir, recursive = FALSE)

all_fimo_raw <- list()
for (subdir in fimo_subdirs) {
  fimo_tsv <- file.path(subdir, "fimo.tsv")
  if (file.exists(fimo_tsv)) {
    label <- basename(subdir)
    parts <- str_match(label, "^(.+)_vs_(.+)$")
    pwm_name <- parts[1, 2]
    target_set <- parts[1, 3]

    df <- tryCatch(
      read_tsv(fimo_tsv, comment = "#", show_col_types = FALSE),
      error = function(e) NULL
    )
    if (!is.null(df) && nrow(df) > 0) {
      df$pwm_name <- pwm_name
      df$target_set <- target_set
      all_fimo_raw[[label]] <- df
      cat("  ", label, ":", nrow(df), "hits\n")
    }
  }
}

all_fimo <- bind_rows(all_fimo_raw)
cat("  Total raw FIMO hits:", nrow(all_fimo), "\n")

# -------------------------------------------------------------------
# 4. Map FIMO hits to promoter regions (gene_id assignment)
# -------------------------------------------------------------------
cat("\n--- Mapping FIMO hits to genes via promoter coordinates ---\n")

# FIMO reports genomic coordinates. Each hit's start/stop should fall within
# a promoter region (chromStart..chromEnd in BED). Map by overlap.

map_hits_to_genes <- function(fimo_df, bed_df) {
  # For each FIMO hit, find which promoter region it falls in
  results <- list()

  for (ts in unique(fimo_df$target_set)) {
    fimo_sub <- fimo_df %>% filter(target_set == ts)
    bed_sub  <- bed_df %>% filter(target_set == ts)

    if (nrow(fimo_sub) == 0 || nrow(bed_sub) == 0) next

    # For each hit, check overlap with any promoter
    mapped <- fimo_sub %>%
      rowwise() %>%
      mutate(
        # Find the promoter that contains this hit
        matched_idx = {
          idx <- which(
            bed_sub$chrom == sequence_name &
            bed_sub$chromStart <= start &
            bed_sub$chromEnd >= stop
          )
          if (length(idx) > 0) idx[1] else NA_integer_
        }
      ) %>%
      ungroup()

    # Add gene info from matched BED entries
    mapped <- mapped %>%
      mutate(
        gene_id_mapped = ifelse(!is.na(matched_idx), bed_sub$name[matched_idx], NA_character_),
        prom_start = ifelse(!is.na(matched_idx), bed_sub$chromStart[matched_idx], NA_integer_),
        prom_end   = ifelse(!is.na(matched_idx), bed_sub$chromEnd[matched_idx], NA_integer_)
      ) %>%
      select(-matched_idx)

    results[[ts]] <- mapped
  }

  bind_rows(results)
}

all_fimo_mapped <- map_hits_to_genes(all_fimo, all_bed)

# Check mapping success
n_mapped <- sum(!is.na(all_fimo_mapped$gene_id_mapped))
n_total  <- nrow(all_fimo_mapped)
cat("  Successfully mapped:", n_mapped, "/", n_total, "hits\n")

if (n_mapped == 0) {
  cat("  WARNING: No hits mapped. Checking coordinate systems...\n")
  cat("  FIMO hit range example:", min(all_fimo$start), "-", max(all_fimo$stop), "\n")
  cat("  BED region range example:", min(all_bed$chromStart), "-", max(all_bed$chromEnd), "\n")
}

# Drop unmapped hits
all_fimo_mapped <- all_fimo_mapped %>%
  filter(!is.na(gene_id_mapped)) %>%
  rename(gene_id = gene_id_mapped)

# Add gene annotation
all_fimo_mapped <- all_fimo_mapped %>%
  left_join(
    coords %>%
      select(gene_id, old_locus_tag, gene_name, bgc_name, gene_type, product, start, end, strand) %>%
      rename(gene_start = start, gene_end = end, gene_strand = strand) %>%
      distinct(gene_id, .keep_all = TRUE),
    by = "gene_id"
  )

# Calculate distance from motif to gene start (TSS proxy)
all_fimo_mapped <- all_fimo_mapped %>%
  mutate(
    # Distance: positive = upstream of TSS, as expected for promoter
    dist_to_tss = case_when(
      gene_strand == "+" ~ gene_start - start,
      gene_strand == "-" ~ stop - gene_end,
      TRUE ~ NA_real_
    )
  )

cat("  Annotated hits:", nrow(all_fimo_mapped), "\n")
cat("  Unique target genes with hits:", n_distinct(all_fimo_mapped$gene_id), "\n")
cat("\n  By target set:\n")
print(table(all_fimo_mapped$target_set))
cat("\n  By PWM:\n")
print(table(all_fimo_mapped$pwm_name))
cat("\n  By BGC (BGC targets only):\n")
bgc_only <- all_fimo_mapped %>% filter(target_set == "BGC")
if (nrow(bgc_only) > 0) print(table(bgc_only$bgc_name))

# -------------------------------------------------------------------
# 5. Write per-target-set FIMO tables
# -------------------------------------------------------------------
cat("\n--- Writing FIMO result tables ---\n")

output_cols <- c("motif_id", "motif_alt_id", "gene_id", "old_locus_tag", "gene_name",
                 "bgc_name", "gene_type", "product", "start", "stop", "strand", "score",
                 "p-value", "q-value", "matched_sequence", "pwm_name", "target_set",
                 "dist_to_tss")

# Ensure all columns exist
for (col in output_cols) {
  if (!col %in% names(all_fimo_mapped)) {
    all_fimo_mapped[[col]] <- NA
  }
}

write_fimo <- function(target, filename) {
  subset <- all_fimo_mapped %>%
    filter(target_set == target) %>%
    select(any_of(output_cols)) %>%
    arrange(`p-value`)
  write_tsv(subset, file.path(OUT_TABLES, filename))
  cat("  Saved:", filename, "(", nrow(subset), "hits,",
      n_distinct(subset$gene_id), "genes )\n")
  subset
}

fimo_bgc  <- write_fimo("BGC",  "fimo_BGC_promoters.tsv")
fimo_sarp <- write_fimo("SARP", "fimo_SARP_promoters.tsv")
fimo_tf   <- write_fimo("TF",   "fimo_TF_promoters.tsv")

# -------------------------------------------------------------------
# 6. BGC SARP motif summary
# -------------------------------------------------------------------
cat("\n--- BGC_SARP_motif_summary ---\n")

if (nrow(fimo_bgc) > 0) {
  bgc_summary <- fimo_bgc %>%
    group_by(bgc_name, gene_id, old_locus_tag, gene_name, product, pwm_name) %>%
    summarise(
      motif_count     = n(),
      top_score       = max(score, na.rm = TRUE),
      best_pvalue     = min(`p-value`, na.rm = TRUE),
      best_qvalue     = min(`q-value`, na.rm = TRUE),
      mean_dist_to_tss = mean(dist_to_tss, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(bgc_name, best_pvalue)

  write_tsv(bgc_summary, file.path(OUT_TABLES, "BGC_SARP_motif_summary.tsv"))
  cat("  Rows:", nrow(bgc_summary), "\n")

  # Overview by BGC × PWM
  bgc_overview <- bgc_summary %>%
    group_by(bgc_name, pwm_name) %>%
    summarise(
      genes_with_hits = n_distinct(gene_id),
      total_hits      = sum(motif_count),
      best_pvalue     = min(best_pvalue),
      .groups = "drop"
    )
  cat("\n  BGC × PWM overview:\n")
  print(as.data.frame(bgc_overview))

  # Total genes per BGC
  bgc_total_genes <- coords %>% filter(gene_type == "BGC_gene") %>% count(bgc_name)
  bgc_hit_genes   <- fimo_bgc %>% distinct(bgc_name, gene_id) %>% count(bgc_name, name = "n_hit")
  bgc_coverage    <- left_join(bgc_total_genes, bgc_hit_genes, by = "bgc_name") %>%
    mutate(n_hit = replace_na(n_hit, 0),
           pct_hit = round(n_hit / n * 100, 1))
  cat("\n  BGC promoter coverage (any motif):\n")
  print(as.data.frame(bgc_coverage))
} else {
  bgc_summary <- tibble()
  write_tsv(bgc_summary, file.path(OUT_TABLES, "BGC_SARP_motif_summary.tsv"))
  cat("  No BGC hits\n")
}

# -------------------------------------------------------------------
# 7. TF promoters SARP hits summary
# -------------------------------------------------------------------
cat("\n--- TF_promoters_SARP_hits ---\n")

if (nrow(fimo_tf) > 0) {
  tf_summary <- fimo_tf %>%
    group_by(gene_id, old_locus_tag, product, pwm_name) %>%
    summarise(
      hit_count   = n(),
      best_score  = max(score, na.rm = TRUE),
      best_pvalue = min(`p-value`, na.rm = TRUE),
      best_qvalue = min(`q-value`, na.rm = TRUE),
      mean_dist_to_tss = mean(dist_to_tss, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(best_pvalue)

  write_tsv(tf_summary, file.path(OUT_TABLES, "TF_promoters_SARP_hits.tsv"))
  cat("  Rows:", nrow(tf_summary), "\n")
  cat("\n  TF × PWM:\n")
  print(as.data.frame(tf_summary %>% select(old_locus_tag, pwm_name, hit_count, best_pvalue)))
} else {
  tf_summary <- tibble()
  write_tsv(tf_summary, file.path(OUT_TABLES, "TF_promoters_SARP_hits.tsv"))
  cat("  No TF hits\n")
}

# -------------------------------------------------------------------
# 8. SARP promoter hits (auto/cross-regulation)
# -------------------------------------------------------------------
cat("\n--- SARP promoter hits ---\n")

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

  cat("  Rows:", nrow(sarp_promo_summary), "\n")
  print(as.data.frame(sarp_promo_summary %>% select(old_locus_tag, pwm_name, hit_count, best_pvalue)))
} else {
  sarp_promo_summary <- tibble()
  cat("  No SARP promoter hits\n")
}

# -------------------------------------------------------------------
# 9. Build SARP target network
# -------------------------------------------------------------------
cat("\n--- Building SARP target network ---\n")

# Map: which SARP → which target promoter
# ActII_ORF4 PWM → actII-orf4 (SC_RS27570) as the SARP
# heptamer PWMs in BGC targets → resident CSR SARP for that BGC
# heptamer PWMs in TF/SARP targets → "SARP_family" (could be any SARP)

bgc_sarp_map <- tribble(
  ~bgc_name, ~SARP_gene_id, ~SARP_old_locus_tag, ~SARP_label,
  "act", "SC_RS27585", "SCO5085", "SCO5085 (act structural SARP)",
  "red", "SC_RS31630", "SCO5877", "redD (red CSR)",
  "cda", "SC_RS18200", "SCO3217", "SCO3217/cdaR (cda CSR)",
  "cpk", "SC_RS33690", "SCO6288", "SCO6288/cpkN (cpk CSR)"
)

network_rows <- list()

# BGC targets: ActII-ORF4 PWM
if (nrow(fimo_bgc) > 0) {
  # ActII-ORF4 PWM in act cluster
  act_actii <- fimo_bgc %>%
    filter(pwm_name == "ActII_ORF4", bgc_name == "act") %>%
    group_by(gene_id, old_locus_tag, gene_name, product, bgc_name) %>%
    summarise(hit_count = n(), best_score = max(score, na.rm = TRUE),
              best_pvalue = min(`p-value`, na.rm = TRUE),
              best_qvalue = min(`q-value`, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(SARP_gene_id = "SC_RS27570", SARP_old_locus_tag = "SCO5082",
           SARP_label = "actII-orf4 (act CSR)",
           target_type = "BGC_gene", pwm_used = "ActII_ORF4",
           comment = "ActII-ORF4 PWM match in act cluster")

  network_rows <- c(network_rows, list(act_actii))

  # ActII-ORF4 PWM in other clusters (cross-cluster)
  cross_actii <- fimo_bgc %>%
    filter(pwm_name == "ActII_ORF4", bgc_name != "act") %>%
    group_by(gene_id, old_locus_tag, gene_name, product, bgc_name) %>%
    summarise(hit_count = n(), best_score = max(score, na.rm = TRUE),
              best_pvalue = min(`p-value`, na.rm = TRUE),
              best_qvalue = min(`q-value`, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(SARP_gene_id = "SC_RS27570", SARP_old_locus_tag = "SCO5082",
           SARP_label = "actII-orf4 (cross-cluster?)",
           target_type = "BGC_gene_cross", pwm_used = "ActII_ORF4",
           comment = "ActII-ORF4-like motif in non-act BGC")

  if (nrow(cross_actii) > 0) network_rows <- c(network_rows, list(cross_actii))

  # Heptamer PWMs in BGCs → resident CSR SARP
  for (pwm in c("SARP_heptamer_strict", "SARP_heptamer_relaxed")) {
    hept_bgc <- fimo_bgc %>%
      filter(pwm_name == pwm) %>%
      group_by(gene_id, old_locus_tag, gene_name, product, bgc_name) %>%
      summarise(hit_count = n(), best_score = max(score, na.rm = TRUE),
                best_pvalue = min(`p-value`, na.rm = TRUE),
                best_qvalue = min(`q-value`, na.rm = TRUE),
                .groups = "drop") %>%
      left_join(bgc_sarp_map, by = "bgc_name") %>%
      mutate(target_type = "BGC_gene", pwm_used = pwm,
             comment = paste("Heptamer match; assigned to", SARP_label))

    if (nrow(hept_bgc) > 0) network_rows <- c(network_rows, list(hept_bgc))
  }
}

# TF targets
if (nrow(fimo_tf) > 0) {
  for (pwm in unique(fimo_tf$pwm_name)) {
    tf_net <- fimo_tf %>%
      filter(pwm_name == pwm) %>%
      group_by(gene_id, old_locus_tag, product) %>%
      summarise(hit_count = n(), best_score = max(score, na.rm = TRUE),
                best_pvalue = min(`p-value`, na.rm = TRUE),
                best_qvalue = min(`q-value`, na.rm = TRUE),
                .groups = "drop") %>%
      mutate(
        bgc_name = NA_character_,
        gene_name = NA_character_,
        SARP_gene_id = ifelse(pwm == "ActII_ORF4", "SC_RS27570", "SARP_family"),
        SARP_old_locus_tag = ifelse(pwm == "ActII_ORF4", "SCO5082", NA_character_),
        SARP_label = ifelse(pwm == "ActII_ORF4", "actII-orf4", "SARP heptamer (any)"),
        target_type = "TF_candidate", pwm_used = pwm,
        comment = paste(pwm, "motif in TF candidate promoter")
      )
    network_rows <- c(network_rows, list(tf_net))
  }
}

# SARP auto/cross-regulation targets
if (nrow(fimo_sarp) > 0) {
  for (pwm in unique(fimo_sarp$pwm_name)) {
    sarp_net <- fimo_sarp %>%
      filter(pwm_name == pwm) %>%
      group_by(gene_id, old_locus_tag, product) %>%
      summarise(hit_count = n(), best_score = max(score, na.rm = TRUE),
                best_pvalue = min(`p-value`, na.rm = TRUE),
                best_qvalue = min(`q-value`, na.rm = TRUE),
                .groups = "drop") %>%
      mutate(
        bgc_name = NA_character_,
        gene_name = NA_character_,
        SARP_gene_id = ifelse(pwm == "ActII_ORF4", "SC_RS27570", "SARP_family"),
        SARP_old_locus_tag = ifelse(pwm == "ActII_ORF4", "SCO5082", NA_character_),
        SARP_label = ifelse(pwm == "ActII_ORF4", "actII-orf4", "SARP heptamer (any)"),
        target_type = "SARP_autoregulation", pwm_used = pwm,
        comment = paste(pwm, "motif in SARP gene promoter (auto/cross-reg)")
      )
    network_rows <- c(network_rows, list(sarp_net))
  }
}

if (length(network_rows) > 0) {
  network <- bind_rows(network_rows) %>%
    select(SARP_gene_id, SARP_old_locus_tag, SARP_label,
           gene_id, old_locus_tag, gene_name, product, bgc_name,
           target_type, pwm_used, hit_count, best_score, best_pvalue, best_qvalue,
           comment) %>%
    arrange(best_pvalue)
} else {
  network <- tibble()
}

write_tsv(network, file.path(OUT_TABLES, "SARP_target_network.tsv"))
cat("  SARP_target_network:", nrow(network), "rows\n")

if (nrow(network) > 0) {
  cat("  By target_type:\n")
  print(table(network$target_type))
  cat("\n  Top 20 by p-value:\n")
  top20 <- network %>%
    head(20) %>%
    select(SARP_label, old_locus_tag, gene_name, bgc_name, target_type, hit_count, best_pvalue)
  print(as.data.frame(top20))
}

# -------------------------------------------------------------------
# 10. Heatmap visualizations
# -------------------------------------------------------------------
cat("\n--- Creating visualizations ---\n")
library(ggplot2)

if (nrow(fimo_bgc) > 0) {
  # Prepare data for heatmap
  heatmap_data <- fimo_bgc %>%
    group_by(bgc_name, gene_id, old_locus_tag, gene_name, pwm_name) %>%
    summarise(
      hit_count = n(),
      neg_log10_pval = -log10(min(`p-value`)),
      .groups = "drop"
    )

  # Gene order by genomic position
  gene_order <- coords %>%
    filter(gene_type == "BGC_gene") %>%
    arrange(bgc_name, start) %>%
    mutate(label = ifelse(is.na(gene_name) | gene_name == "",
                          old_locus_tag,
                          paste0(gene_name, " (", old_locus_tag, ")")))

  heatmap_data <- heatmap_data %>%
    left_join(gene_order %>% select(gene_id, label), by = "gene_id") %>%
    mutate(label = factor(label, levels = rev(gene_order$label)))

  # Only include genes with hits
  genes_with_hits <- unique(heatmap_data$gene_id)
  gene_order_hits <- gene_order %>% filter(gene_id %in% genes_with_hits)

  heatmap_data_hits <- heatmap_data %>%
    filter(gene_id %in% genes_with_hits) %>%
    mutate(label = factor(label, levels = rev(gene_order_hits$label)))

  p1 <- ggplot(heatmap_data_hits, aes(x = pwm_name, y = label, fill = neg_log10_pval)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = hit_count), size = 2.5, color = "black") +
    facet_grid(bgc_name ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient(low = "lightyellow", high = "red3",
                        name = expression(-log[10](p))) +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.y = element_text(size = 6),
      strip.text.y = element_text(angle = 0, face = "bold", size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid = element_blank()
    ) +
    labs(x = "Motif PWM", y = "BGC gene (promoter)",
         title = "SARP motif hits in BGC promoters (-log10 p-value, number = hit count)")

  h <- max(4, nrow(gene_order_hits) * 0.25 + 2)
  ggsave(file.path(OUT_FIGURES, "heatmap_BGC_SARP_motif_hits.pdf"),
         p1, width = 7, height = h, limitsize = FALSE)
  ggsave(file.path(OUT_FIGURES, "heatmap_BGC_SARP_motif_hits.png"),
         p1, width = 7, height = h, dpi = 150, limitsize = FALSE)
  cat("  Saved: heatmap_BGC_SARP_motif_hits.pdf/png\n")

  # Bar plot: total hits per BGC
  bgc_bar <- fimo_bgc %>%
    group_by(bgc_name, pwm_name) %>%
    summarise(n_genes = n_distinct(gene_id), total_hits = n(), .groups = "drop")

  p3 <- ggplot(bgc_bar, aes(x = bgc_name, y = total_hits, fill = pwm_name)) +
    geom_col(position = "dodge") +
    geom_text(aes(label = total_hits), position = position_dodge(0.9), vjust = -0.3, size = 3) +
    theme_minimal(base_size = 12) +
    scale_fill_brewer(palette = "Set2", name = "PWM") +
    labs(x = "BGC", y = "Total FIMO hits",
         title = "SARP motif scan: hits per BGC")
  ggsave(file.path(OUT_FIGURES, "barplot_BGC_hits_by_motif.pdf"), p3, width = 7, height = 4)
  ggsave(file.path(OUT_FIGURES, "barplot_BGC_hits_by_motif.png"), p3, width = 7, height = 4, dpi = 150)
  cat("  Saved: barplot_BGC_hits_by_motif.pdf/png\n")
}

# TF heatmap
if (nrow(fimo_tf) > 0) {
  tf_heat <- fimo_tf %>%
    group_by(gene_id, old_locus_tag, product, pwm_name) %>%
    summarise(
      hit_count = n(),
      neg_log10_pval = -log10(min(`p-value`)),
      .groups = "drop"
    ) %>%
    mutate(label = paste0(old_locus_tag, " (", str_trunc(product, 30), ")"))

  p2 <- ggplot(tf_heat, aes(x = pwm_name, y = reorder(label, neg_log10_pval), fill = neg_log10_pval)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = hit_count), size = 3) +
    scale_fill_gradient(low = "lightyellow", high = "steelblue",
                        name = expression(-log[10](p))) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank()) +
    labs(x = "Motif PWM", y = "TF candidate",
         title = "SARP motif hits in TF candidate promoters")
  ggsave(file.path(OUT_FIGURES, "heatmap_TF_SARP_motif_hits.pdf"), p2, width = 8, height = 5)
  ggsave(file.path(OUT_FIGURES, "heatmap_TF_SARP_motif_hits.png"), p2, width = 8, height = 5, dpi = 150)
  cat("  Saved: heatmap_TF_SARP_motif_hits.pdf/png\n")
}

# -------------------------------------------------------------------
# 11. Save session data
# -------------------------------------------------------------------
save(all_fimo_mapped, bgc_summary, tf_summary, sarp_promo_summary, network,
     fimo_bgc, fimo_sarp, fimo_tf,
     coords, sarp_v2, tf_cand, gm, all_bed,
     file = file.path(RUN_DIR, "motif_scan_session.RData"))
cat("  Session data saved\n")

cat("\n=== Aggregation completed at", format(Sys.time()), "===\n")
