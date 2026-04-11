# Alignment Report — M145 RNA-seq

## 1. Project Overview

- **Organism**: *Streptomyces coelicolor* A3(2) M145
- **Experiment**: 3 conditions (M145_1, M145_2, M145_3) × 3 biological replicates = 9 samples
- **Sequencing**: Paired-end 150 bp (Illumina)
- **Alignment tool**: HISAT2 2.2.1 (see note below)
- **Reference genome**: GCF_000203835.1 (*S. coelicolor* A3(2), ~8.67 Mb)

### Note on Aligner Selection

The original prompt specified STAR 2.7.11b. However, STAR 2.7.11b is **broken on macOS arm64** — all available builds (conda bioconda builds haf7d672_6/7/8 and Homebrew) report 0 input reads for any FASTQ file. This is a confirmed platform-specific I/O bug where `ifstream` fails to read file contents despite successfully opening the file.

HISAT2 2.2.1 was used as a replacement. For prokaryotic RNA-seq (no introns/splicing in *S. coelicolor*), HISAT2 with `--no-spliced-alignment` is scientifically equivalent and commonly used in bacterial transcriptomics.

---

## 2. Reference Genome

- **FASTA**: `GCF_000203835.1_ASM20383v1_genomic.fna`
- **GTF**: `genomic.gtf`
- **Source**: NCBI RefSeq GCF_000203835.1
- **Chromosomes**: NC_003888.3 (chromosome, 8,667,507 bp), NC_003903.1 (SCP1 plasmid, 356,023 bp), NC_003904.1 (SCP2 plasmid, 31,317 bp)

---

## 3. Input FASTQ Files

All samples used **trimmed** FASTQ files (post-fastp QC).

| Sample | Condition | Replicate | Source |
|--------|-----------|-----------|--------|
| M145_1_1 | M145_1 | 1 | TRIMMED |
| M145_1_2 | M145_1 | 2 | TRIMMED |
| M145_1_3 | M145_1 | 3 | TRIMMED |
| M145_2_1 | M145_2 | 1 | TRIMMED |
| M145_2_3 | M145_2 | 3 | TRIMMED |
| M145_2_4 | M145_2 | 4 | TRIMMED |
| M145_3_2 | M145_3 | 2 | TRIMMED |
| M145_3_3 | M145_3 | 3 | TRIMMED |
| M145_3_4 | M145_3 | 4 | TRIMMED |

Trimmed FASTQ source directory: `01_qc/analysis/01_qc_260127_v1/fastp/`

---

## 4. Alignment Parameters

```bash
hisat2 \
  -x hisat2_index/M145 \
  -1 <R1_trimmed.fastq.gz> \
  -2 <R2_trimmed.fastq.gz> \
  --threads 8 \
  --no-spliced-alignment \
  --new-summary \
  --rg-id <sample_id> \
  --rg "SM:<sample_id>" \
  --rg "PL:ILLUMINA" \
| samtools sort -@ 4 -o <sample>.Aligned.sortedByCoord.out.bam -
samtools index <sample>.Aligned.sortedByCoord.out.bam
```

Key parameters:
- `--no-spliced-alignment`: Disables splice-aware alignment, appropriate for prokaryotic genomes without introns
- Piped to `samtools sort` for coordinate-sorted BAM output
- BAM indexed with `samtools index` for downstream visualization and analysis

---

## 5. Alignment QC Results

| Sample | Total Reads | Uniquely Mapped (%) | Multi-mapped (%) | Overall Rate (%) |
|--------|------------|--------------------|-----------------|-----------------:|
| M145_1_1 | 7,099,866 | 96.66% | 0.7% | 98.47% |
| M145_1_2 | 6,747,804 | 96.86% | 0.59% | 98.54% |
| M145_1_3 | 8,786,477 | 96.8% | 0.58% | 98.49% |
| M145_2_1 | 7,423,361 | 96.65% | 1.13% | 98.7% |
| M145_2_3 | 6,649,141 | 95.93% | 1.75% | 98.63% |
| M145_2_4 | 7,404,965 | 96.21% | 1.51% | 98.67% |
| M145_3_2 | 7,734,535 | 96.42% | 0.93% | 98.41% |
| M145_3_3 | 6,976,378 | 96.78% | 0.71% | 98.51% |
| M145_3_4 | 6,828,711 | 96.65% | 0.74% | 98.46% |
| **Average** | **7,294,582** | **96.55%** | **0.96%** | **98.54%** |

### Summary

- **Average overall alignment rate**: 98.54%
- **Average uniquely mapped**: 96.55%
- **Average multi-mapped**: 0.96%
- All samples show consistent, high-quality alignment (≥98.4%)
- Low multi-mapping rate (<2%) indicates good specificity
- The high alignment rate is expected for this well-characterized genome

---

## 6. Output BAM Files

- `bam/M145_1_1.Aligned.sortedByCoord.out.bam` (508 MB)
- `bam/M145_1_2.Aligned.sortedByCoord.out.bam` (484 MB)
- `bam/M145_1_3.Aligned.sortedByCoord.out.bam` (617 MB)
- `bam/M145_2_1.Aligned.sortedByCoord.out.bam` (528 MB)
- `bam/M145_2_3.Aligned.sortedByCoord.out.bam` (483 MB)
- `bam/M145_2_4.Aligned.sortedByCoord.out.bam` (526 MB)
- `bam/M145_3_2.Aligned.sortedByCoord.out.bam` (543 MB)
- `bam/M145_3_3.Aligned.sortedByCoord.out.bam` (476 MB)
- `bam/M145_3_4.Aligned.sortedByCoord.out.bam` (475 MB)

All BAMs are coordinate-sorted and indexed (.bai).

---

## 7. M145-Specific Notes

- **High GC content (~72%)**: The high GC content of *S. coelicolor* can affect read mapping. The consistently high alignment rates (>98%) indicate no significant GC-related mapping issues.
- **Operon structure**: Prokaryotic operons produce polycistronic mRNAs. This does not affect alignment but is relevant for downstream gene-level quantification (featureCounts handles this via `--fracOverlap` or overlap assignment modes).
- **Plasmid mapping**: Reads mapping to SCP1 and SCP2 plasmids are included in the BAM files. These should be considered in downstream analysis depending on the biological question.

---

## 8. Next Steps

1. **Gene quantification**: Run featureCounts (or FADU for prokaryotic-optimized counting) on the sorted BAMs using the GTF annotation.
2. **Differential expression**: Use DESeq2 with the count matrix from featureCounts.
3. **Quality assessment**: Consider running `samtools flagstat` or `samtools stats` for additional BAM-level QC if needed.

---

## 9. Materials & Methods

Trimmed paired-end reads (150 bp) were aligned to the *Streptomyces coelicolor* A3(2) reference genome (GCF_000203835.1) using HISAT2 v2.2.1 with `--no-spliced-alignment` for prokaryotic alignment. Output was piped to SAMtools v1.21 for coordinate sorting and BAM indexing. On average, 98.54% of read pairs aligned to the reference, with 96.55% mapping uniquely. Alignment statistics were extracted from HISAT2 summary reports for all nine samples (three biological replicates per condition).

---

*Generated: 2026-01-27*
*Run directory: `02_alignment_260127_v1`*
