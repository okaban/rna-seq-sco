#!/usr/bin/env Rscript
# 10_SARP_motif_scan pipeline — M145
# Steps 1-2: Coordinate extraction & promoter region preparation

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

cat("=== 10_SARP_motif_scan started at", format(Sys.time()), "===\n")

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

REF_DIR      <- "/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1"
GENOME_FNA   <- file.path(REF_DIR, "GCF_000203835.1_ASM20383v1_genomic.fna")
GENOME_GFF   <- file.path(REF_DIR, "genomic.gff")

GENE_MASTER  <- file.path(ANNOT_RUN, "tables/gene_master_with_BGC_regulators.tsv")
BGC_DEF      <- file.path(ANNOT_RUN, "tables/BGC_definition_manual.tsv")
SARP_V2_LIST <- file.path(SARP_V2_RUN, "tables/SARP_list_M145_v2.tsv")
TF_CANDIDATES <- file.path(TF_RUN, "tables/TF_candidates_experimental_priority.tsv")

OUT_TABLES    <- file.path(RUN_DIR, "tables")
OUT_PROMOTERS <- file.path(RUN_DIR, "promoters")
OUT_PWM       <- file.path(RUN_DIR, "pwm")

# -------------------------------------------------------------------
# 1. Read inputs
# -------------------------------------------------------------------
cat("\n--- Step 1: Reading inputs ---\n")

gm      <- read_tsv(GENE_MASTER, show_col_types = FALSE)
bgc_def <- read_tsv(BGC_DEF, show_col_types = FALSE)
sarp_v2 <- read_tsv(SARP_V2_LIST, show_col_types = FALSE)
tf_cand <- read_tsv(TF_CANDIDATES, show_col_types = FALSE)

cat("  gene_master:", nrow(gm), "rows\n")
cat("  BGC definition:", nrow(bgc_def), "genes across",
    length(unique(bgc_def$bgc_name)), "BGCs\n")
cat("  SARP v2 list:", nrow(sarp_v2), "SARPs\n")
cat("  TF candidates:", nrow(tf_cand), "TFs\n")

# -------------------------------------------------------------------
# 2. Build coordinates table
# -------------------------------------------------------------------
cat("\n--- Step 2: Building coordinates ---\n")

# BGC genes with coordinates
bgc_coords <- bgc_def %>%
  left_join(gm %>% select(gene_id, contig, start, end, strand, product),
            by = "gene_id") %>%
  mutate(gene_type = "BGC_gene")

# SARP genes with coordinates
sarp_coords <- sarp_v2 %>%
  select(gene_id, old_locus_tag) %>%
  left_join(gm %>% select(gene_id, contig, start, end, strand, product, bgc_name),
            by = "gene_id") %>%
  mutate(gene_type = "SARP",
         role = "SARP")

# TF candidates with coordinates
tf_coords <- tf_cand %>%
  select(gene_id, old_locus_tag) %>%
  left_join(gm %>% select(gene_id, contig, start, end, strand, product, bgc_name),
            by = c("gene_id")) %>%
  mutate(gene_type = "TF_candidate",
         role = "TF_candidate")

# Combine all into one coordinates table
all_coords <- bind_rows(
  bgc_coords %>% select(bgc_name, gene_id, old_locus_tag, gene_name,
                         role, contig, start, end, strand, product, gene_type),
  sarp_coords %>% select(bgc_name, gene_id, old_locus_tag,
                          role, contig, start, end, strand, product, gene_type) %>%
    mutate(gene_name = NA_character_),
  tf_coords %>% select(bgc_name, gene_id, old_locus_tag,
                        role, contig, start, end, strand, product, gene_type) %>%
    mutate(gene_name = NA_character_)
) %>%
  distinct(gene_id, gene_type, .keep_all = TRUE)

cat("  Total coordinate entries:", nrow(all_coords), "\n")
cat("  By gene_type:\n")
print(table(all_coords$gene_type))

write_tsv(all_coords, file.path(OUT_TABLES, "BGC_and_TF_coordinates.tsv"))
cat("  Saved: BGC_and_TF_coordinates.tsv\n")

# -------------------------------------------------------------------
# 3. Define promoter regions (BED format for bedtools)
# -------------------------------------------------------------------
cat("\n--- Step 3: Defining promoter regions ---\n")

# Get genome contig lengths from genome FASTA index
fai_file <- paste0(GENOME_FNA, ".fai")
if (!file.exists(fai_file)) {
  cat("  Creating genome index...\n")
  system2("samtools", args = c("faidx", shQuote(GENOME_FNA)))
}
fai <- read_tsv(fai_file, col_names = c("contig", "length", "offset", "bases_per_line", "bytes_per_line"),
                show_col_types = FALSE)
contig_lengths <- setNames(fai$length, fai$contig)

# Helper function to create promoter BED entry
make_promoter_bed <- function(contig, start, end, strand, gene_id, upstream = 500) {
  clen <- contig_lengths[contig]
  if (strand == "+") {
    prom_start <- max(0, start - upstream)
    prom_end   <- max(0, start - 1)
  } else {
    prom_start <- min(clen, end + 1)
    prom_end   <- min(clen, end + upstream)
  }
  if (prom_start >= prom_end) return(NULL)
  tibble(chrom = contig, chromStart = prom_start, chromEnd = prom_end,
         name = gene_id, score = 0, strand = strand)
}

# BGC promoters
bgc_genes_unique <- all_coords %>%
  filter(gene_type == "BGC_gene") %>%
  distinct(gene_id, .keep_all = TRUE)

bgc_bed <- pmap_dfr(
  list(bgc_genes_unique$contig, bgc_genes_unique$start,
       bgc_genes_unique$end, bgc_genes_unique$strand,
       bgc_genes_unique$gene_id),
  function(c, s, e, st, g) make_promoter_bed(c, s, e, st, g)
)

# SARP promoters
sarp_genes <- all_coords %>%
  filter(gene_type == "SARP") %>%
  distinct(gene_id, .keep_all = TRUE)

sarp_bed <- pmap_dfr(
  list(sarp_genes$contig, sarp_genes$start,
       sarp_genes$end, sarp_genes$strand,
       sarp_genes$gene_id),
  function(c, s, e, st, g) make_promoter_bed(c, s, e, st, g)
)

# TF promoters
tf_genes <- all_coords %>%
  filter(gene_type == "TF_candidate") %>%
  distinct(gene_id, .keep_all = TRUE)

tf_bed <- pmap_dfr(
  list(tf_genes$contig, tf_genes$start,
       tf_genes$end, tf_genes$strand,
       tf_genes$gene_id),
  function(c, s, e, st, g) make_promoter_bed(c, s, e, st, g)
)

cat("  BGC promoter regions:", nrow(bgc_bed), "\n")
cat("  SARP promoter regions:", nrow(sarp_bed), "\n")
cat("  TF promoter regions:", nrow(tf_bed), "\n")

# Write BED files
bgc_bed_file  <- file.path(OUT_PROMOTERS, "bgc_promoters.bed")
sarp_bed_file <- file.path(OUT_PROMOTERS, "sarp_promoters.bed")
tf_bed_file   <- file.path(OUT_PROMOTERS, "tf_promoters.bed")

write_tsv(bgc_bed,  bgc_bed_file,  col_names = FALSE)
write_tsv(sarp_bed, sarp_bed_file, col_names = FALSE)
write_tsv(tf_bed,   tf_bed_file,   col_names = FALSE)

# -------------------------------------------------------------------
# 4. Extract promoter sequences with bedtools
# -------------------------------------------------------------------
cat("\n--- Step 4: Extracting promoter sequences ---\n")

bgc_fasta  <- file.path(OUT_PROMOTERS, "promoters_BGCs_500bp.fasta")
sarp_fasta <- file.path(OUT_PROMOTERS, "promoters_SARPs_500bp.fasta")
tf_fasta   <- file.path(OUT_PROMOTERS, "promoters_TFs_500bp.fasta")

for (info in list(
  list(bed = bgc_bed_file,  fa = bgc_fasta,  label = "BGC"),
  list(bed = sarp_bed_file, fa = sarp_fasta, label = "SARP"),
  list(bed = tf_bed_file,   fa = tf_fasta,   label = "TF")
)) {
  cmd <- sprintf("bedtools getfasta -fi %s -bed %s -name -s -fo %s",
                 shQuote(GENOME_FNA), shQuote(info$bed), shQuote(info$fa))
  cat("  Extracting", info$label, "promoters...\n")
  ret <- system(cmd, intern = FALSE)
  if (ret != 0) {
    cat("  WARNING: bedtools failed for", info$label, "\n")
  } else {
    n_seqs <- as.integer(system(paste("grep -c '^>'", shQuote(info$fa)), intern = TRUE))
    cat("  ", info$label, ":", n_seqs, "sequences extracted\n")
  }
}

cat("\n=== Steps 1-4 (coordinates & promoters) completed ===\n")
