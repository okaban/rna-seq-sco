# QC Report: M145 RNA-seq

## Overview

- **Project**: M145_RNA-seq
- **Organism**: *Streptomyces coelicolor* A3(2) M145
- **Run ID**: 01_qc_260127_v1
- **Pipeline**: FastQC (raw) → fastp (trim/filter) → FastQC (trimmed) → MultiQC
- **Samples**: 9 (3 conditions × 3 replicates, paired-end 150 bp)

| Condition | Replicates |
|-----------|------------|
| M145_1 (Timepoint 1) | M145_1_1, M145_1_2, M145_1_3 |
| M145_2 (Timepoint 2) | M145_2_1, M145_2_3, M145_2_4 |
| M145_3 (Timepoint 3) | M145_3_2, M145_3_3, M145_3_4 |

---

## Raw Data Summary

| Sample | Raw Read Pairs | Read Length | %GC |
|--------|---------------|-------------|-----|
| M145_1_1 | 7,212,652 | 150 bp | 68% |
| M145_1_2 | 6,851,465 | 150 bp | 68% |
| M145_1_3 | 8,933,580 | 150 bp | 68% |
| M145_2_1 | 7,560,906 | 150 bp | 68% |
| M145_2_3 | 6,993,324 | 150 bp | 68% |
| M145_2_4 | 7,513,310 | 150 bp | 68% |
| M145_3_2 | 7,892,330 | 150 bp | 68% |
| M145_3_3 | 7,422,280 | 150 bp | 67% |
| M145_3_4 | 7,294,192 | 150 bp | 67% |

**Total raw read pairs**: 67,674,039

---

## fastp Trimming & Filtering Results

Parameters: `--qualified_quality_phred 20 --length_required 50 --detect_adapter_for_pe --correction --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20`

| Sample | Raw Reads | Clean Reads | Pass Rate | Q20 (before) | Q20 (after) | Q30 (before) | Q30 (after) |
|--------|-----------|-------------|-----------|-------------|-------------|-------------|-------------|
| M145_1_1 | 14,425,304 | 14,199,732 | 98.4% | 98.1% | 98.7% | 93.7% | 94.8% |
| M145_1_2 | 13,702,930 | 13,495,608 | 98.5% | 98.1% | 98.7% | 93.7% | 94.8% |
| M145_1_3 | 17,867,160 | 17,572,954 | 98.4% | 98.1% | 98.7% | 93.7% | 94.8% |
| M145_2_1 | 15,121,812 | 14,846,722 | 98.2% | 98.2% | 98.8% | 94.0% | 95.0% |
| M145_2_3 | 13,986,648 | 13,298,282 | 95.1% | 98.2% | 98.8% | 94.0% | 94.9% |
| M145_2_4 | 15,026,620 | 14,809,930 | 98.6% | 98.2% | 98.8% | 94.0% | 95.0% |
| M145_3_2 | 15,784,660 | 15,469,070 | 98.0% | 98.2% | 98.8% | 94.0% | 95.0% |
| M145_3_3 | 14,844,560 | 13,952,756 | 94.0% | 98.3% | 98.9% | 94.2% | 95.1% |
| M145_3_4 | 14,588,384 | 13,657,422 | 93.6% | 98.3% | 98.8% | 94.2% | 95.1% |

**Overall**: 135,348,078 → 131,302,476 reads (97.0% retained)

### Filtering Breakdown (per sample)

| Sample | Low Quality | Too Short | Too Many N | Adapter |
|--------|------------|-----------|-----------|---------|
| M145_1_1 | 190,456 | 34,258 | 396 | 156,716 |
| M145_1_2 | 178,702 | 27,966 | 372 | 162,732 |
| M145_1_3 | 251,486 | 41,934 | 378 | 207,927 |
| M145_2_1 | 184,024 | 90,454 | 380 | 198,005 |
| M145_2_3 | 174,742 | 512,994 | 294 | 386,316 |
| M145_2_4 | 187,964 | 28,122 | 368 | 202,539 |
| M145_3_2 | 199,858 | 115,188 | 322 | 209,582 |
| M145_3_3 | 179,856 | 711,322 | 332 | 491,438 |
| M145_3_4 | 177,176 | 752,988 | 332 | 515,452 |

---

## FastQC Module Summary (Raw vs Trimmed)

### Key Observations

1. **Per base sequence quality**: All samples PASS (both raw and trimmed). Base quality is excellent.
2. **Per tile sequence quality**: Several raw samples showed FAIL; fastp trimming improved most to PASS/WARN.
3. **Per base sequence content**: WARN/FAIL — typical RNA-seq bias from random hexamer priming at read start positions. Not a data quality concern.
4. **Per sequence GC content**: WARN/FAIL — expected for *S. coelicolor* (genome GC ~72%), which deviates from FastQC's theoretical normal distribution.
5. **Adapter content**: All PASS in both raw and trimmed. No adapter contamination detected.
6. **Overrepresented sequences**: WARN/FAIL — likely rRNA fragments. Consider rRNA depletion assessment at alignment step.
7. **Sequence duplication levels**: All FAIL — expected for RNA-seq where highly expressed genes produce many identical reads.

### Conclusion

The data quality is **suitable for downstream analysis**. The FAIL/WARN flags are explained by:
- High GC content inherent to *S. coelicolor* (~72% GC)
- Normal RNA-seq duplication from highly expressed transcripts
- Standard hexamer priming bias at read starts

fastp filtering retained **93.6–98.6%** of reads with improved Q20/Q30 scores. The trimmed data is recommended for alignment.

### M145 特有の注意点

- *S. coelicolor* A3(2) のゲノム GC 含量は約 72% と極めて高く、FastQC の per sequence GC content モジュールで WARN/FAIL が出るのは正常な挙動である。
- 過剰出現配列（overrepresented sequences）の多くは rRNA 由来と考えられる。アラインメント後に rRNA マッピング率を確認し、必要に応じて rRNA 除去率を評価すること。

---

## Trimmed FASTQ Files for Alignment

以降のアラインメントステップ（02_alignment）では、以下の **fastp トリミング済みファイル** を入力として使用してください。

| Sample | R1 (trimmed) | R2 (trimmed) |
|--------|-------------|-------------|
| M145_1_1 | `fastp/M145_1_1_trimmed_R1.fastq.gz` | `fastp/M145_1_1_trimmed_R2.fastq.gz` |
| M145_1_2 | `fastp/M145_1_2_trimmed_R1.fastq.gz` | `fastp/M145_1_2_trimmed_R2.fastq.gz` |
| M145_1_3 | `fastp/M145_1_3_trimmed_R1.fastq.gz` | `fastp/M145_1_3_trimmed_R2.fastq.gz` |
| M145_2_1 | `fastp/M145_2_1_trimmed_R1.fastq.gz` | `fastp/M145_2_1_trimmed_R2.fastq.gz` |
| M145_2_3 | `fastp/M145_2_3_trimmed_R1.fastq.gz` | `fastp/M145_2_3_trimmed_R2.fastq.gz` |
| M145_2_4 | `fastp/M145_2_4_trimmed_R1.fastq.gz` | `fastp/M145_2_4_trimmed_R2.fastq.gz` |
| M145_3_2 | `fastp/M145_3_2_trimmed_R1.fastq.gz` | `fastp/M145_3_2_trimmed_R2.fastq.gz` |
| M145_3_3 | `fastp/M145_3_3_trimmed_R1.fastq.gz` | `fastp/M145_3_3_trimmed_R2.fastq.gz` |
| M145_3_4 | `fastp/M145_3_4_trimmed_R1.fastq.gz` | `fastp/M145_3_4_trimmed_R2.fastq.gz` |

---

## Materials & Methods (QC Summary)

> RNA-seq raw reads (paired-end, 150 bp) were assessed for quality using FastQC v0.12.1. Based on the initial QC results, reads were trimmed and filtered using fastp v1.1.0 with the following parameters: minimum Phred quality score of 20, minimum read length of 50 bp, automatic adapter detection for paired-end reads, overlap-based error correction, and sliding window quality trimming (window size 4, mean quality 20). After filtering, 93.6–98.6% of reads were retained across all samples with Q30 scores of 94.8–95.1%. Aggregated quality reports were generated using MultiQC v1.33. All tools were managed within a Conda environment (Python 3.11).

---

## Output Files

| Directory | Contents |
|-----------|----------|
| `raw_fastqc/` | FastQC reports on raw data (18 HTML + ZIP) |
| `multiqc/` | MultiQC aggregated reports (`multiqc_raw_fastqc.html`, `multiqc_final_qc.html`) |
| `fastp/` | Trimmed FASTQ files + fastp per-sample HTML/JSON reports + `fastp_summary.tsv` |
| `trimmed_fastqc/` | FastQC reports on trimmed data (18 HTML + ZIP) |
| `figures/` | Publication-quality QC figures (PDF/SVG) |
| `pipeline.log` | Execution log with timestamps |
| `sample_table.tsv` | Sample metadata table |

---

## Tools & Versions

| Tool | Version |
|------|---------|
| FastQC | 0.12.1 |
| fastp | 1.1.0 |
| MultiQC | 1.33 |
| Python | 3.11 |
| Conda env | rna-seq |
