#!/usr/bin/env bash
# Build a combined nucleotide BLAST DB from 833 Streptomyces genomes.
# Each FASTA header is rewritten to "<accession>|<original_id>" so that the
# genome of origin can be recovered from blast hit subject IDs.
set -euo pipefail

OUT=/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/73_SC_RS17645
DATA=/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/data/ncbi_genomes_all_species/extracted/ncbi_dataset/data
COMBINED="$OUT/strep833_combined.fna"
DB="$OUT/strep833_db"

: > "$COMBINED"
n=0
for d in "$DATA"/GCF_*; do
  acc=$(basename "$d")
  fna=$(ls "$d"/*_genomic.fna 2>/dev/null | head -1)
  [ -z "$fna" ] && continue
  # Prefix every contig header with the accession
  awk -v acc="$acc" '/^>/ {sub(/^>/, ">" acc "|"); print; next} {print}' "$fna" >> "$COMBINED"
  n=$((n+1))
done
echo "Combined $n genomes into $COMBINED"
ls -lh "$COMBINED"

makeblastdb -in "$COMBINED" -dbtype nucl -out "$DB" -title "Strep833"
echo "DB built at $DB"
