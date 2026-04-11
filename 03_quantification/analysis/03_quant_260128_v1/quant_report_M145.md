# Quantification Report — M145 RNA-seq

## 1. Project Overview

- **Organism**: *Streptomyces coelicolor* A3(2) M145
- **Experiment**: 3 conditions (M145_1, M145_2, M145_3) x 3 biological replicates = 9 samples
- **Sequencing**: Paired-end 150 bp (Illumina)
- **Quantification tool**: featureCounts v2.1.1 (subread)
- **Reference genome**: GCF_000203835.1 (*S. coelicolor* A3(2), ~8.67 Mb)

---

## 2. Input BAM Files

9 coordinate-sorted BAM files from HISAT2 v2.2.1 alignment (`--no-spliced-alignment`).

Source directory: `02_alignment/analysis/02_alignment_260127_v1/bam/`

All BAMs were generated with HISAT2 (STAR was unavailable on macOS arm64) and indexed with samtools.

---

## 3. Annotation

- **GTF**: `genomic.gtf` from NCBI RefSeq GCF_000203835.1
- **Feature composition**: gene (8,276), CDS (8,200), start_codon (8,122), stop_codon (8,098), exon (87, rRNA/ncRNA only), transcript (87)
- **Feature type used**: `-t gene` (covers all annotated genes including CDS, rRNA, tRNA, ncRNA)
- **Grouping attribute**: `-g gene_id` (e.g., `SC_RS02065`)

---

## 4. Strandedness Determination

Prior to the full run, strandedness was tested on M145_1_1 with `-s 0`, `-s 1`, `-s 2`:

| `-s` value | Setting | Assigned (M145_1_1) |
|-----------|---------|-------------------|
| `-s 0` | Unstranded | 6,516,052 (89.7%) |
| `-s 1` | Forward | 984,941 (13.6%) |
| `-s 2` | Reverse | 5,674,291 (78.1%) |

**Conclusion**: `-s 0` (unstranded) yielded the highest assignment rate and was adopted for all samples. The library shows a reverse-strand bias (`-s 2` = 78.1%) but unstranded counting captures the maximum information.

---

## 5. featureCounts Command

```bash
featureCounts \
  -T 6 \
  -p --countReadPairs \
  -B \
  -C \
  -s 0 \
  -t gene \
  -g gene_id \
  -a "${REF_GTF}" \
  -o counts/featureCounts_M145.txt \
  ${BAM_DIR}/*.bam
```

Key parameters:
- `-p --countReadPairs`: Count fragments (read pairs) rather than individual reads
- `-B`: Require both ends of a pair to be mapped
- `-C`: Exclude chimeric fragments
- `-s 0`: Unstranded counting
- `-t gene`: Count at the gene feature level (8,276 genes in GTF)

---

## 6. Quantification Results

| Sample | Total Counts | Nonzero Genes | Top-10 Fraction | Assign Rate |
|--------|-------------|--------------|----------------|------------|
| M145_1_1 | 6,516,052 | 7700 / 8275 | 20.90% | 89.7% |
| M145_1_2 | 6,208,949 | 7684 / 8275 | 20.73% | 90.3% |
| M145_1_3 | 8,079,327 | 7730 / 8275 | 21.08% | 90.3% |
| M145_2_1 | 6,861,814 | 7747 / 8275 | 30.42% | 89.2% |
| M145_2_3 | 6,087,842 | 7734 / 8275 | 29.73% | 86.4% |
| M145_2_4 | 6,816,908 | 7735 / 8275 | 29.44% | 87.7% |
| M145_3_2 | 7,081,629 | 7856 / 8275 | 33.22% | 89.2% |
| M145_3_3 | 6,416,864 | 7848 / 8275 | 37.53% | 90.2% |
| M145_3_4 | 6,284,055 | 7860 / 8275 | 36.38% | 90.3% |
| **Average** | **6,705,938** | **7,766** | **28.83%** | **89.3%** |

### Key Observations

- **Average library size**: 6,705,938 assigned fragments per sample
- **Average assignment rate**: 89.3% (range: 86.4%--90.3%)
- **Average detected genes**: 7,766 / 8275 (93.8%)
- **Top-10 gene concentration**: M145_1 samples (~21%) show more even distribution than M145_3 samples (~36%), suggesting differential expression of highly expressed genes across conditions
- Library sizes are consistent across samples (6.1M--8.1M), suitable for DESeq2 normalization

---

## 7. M145-Specific Notes

- **High GC content (~72%)**: No apparent GC-related counting bias observed. Assignment rates are uniformly high.
- **Operon structure**: Prokaryotic polycistronic mRNAs can cause reads to overlap multiple gene boundaries. The Unassigned_Ambiguity category (1.8--3.7%) reflects this. featureCounts defaults to not assigning ambiguous reads, which is conservative but appropriate for DESeq2 input.
- **FADU alternative**: For follow-up analysis, FADU (Fragment Assignment Diagnosis Utility) could provide complementary counts optimized for prokaryotic gene structures, particularly for short overlapping genes in operons. This was not run in this analysis but can be added.

---

## 8. Output Files

- `counts/featureCounts_M145.txt` — Gene-level count matrix (8,275 genes x 9 samples)
- `counts/featureCounts_M145.txt.summary` — Assignment summary per sample
- `summary/quant_sample_table.tsv` — Sample metadata with BAM paths
- `summary/counts_summary.tsv` — Per-sample library size and QC metrics
- `figures/library_size_barplot.pdf/.svg` — Library size visualization
- `figures/count_distribution_boxplot.pdf/.svg` — Gene count distribution
- `figures/assignment_summary.pdf/.svg` — featureCounts assignment breakdown

---

## 9. Next Steps (DESeq2)

1. **Count matrix**: Use `featureCounts_M145.txt` columns 7+ (skip Geneid/Chr/Start/End/Strand/Length)
2. **Sample metadata**: Use `quant_sample_table.tsv` for condition/replicate information
3. **Design formula**: `~ condition` (M145_1 vs M145_2 vs M145_3) for pairwise comparisons
4. **Pre-filtering**: Consider removing genes with very low counts (e.g., < 10 total across all samples)
5. **Normalization**: DESeq2 will handle library size normalization internally (median-of-ratios)

---

## 10. Materials & Methods

Gene-level read counts were generated from coordinate-sorted BAM files using featureCounts v2.1.1 (subread package) with the NCBI RefSeq GTF annotation (GCF_000203835.1). Paired-end fragments were counted at the gene level (`-t gene -g gene_id`) in unstranded mode (`-s 0`), requiring both mates to be aligned (`-B`) and excluding chimeric fragments (`-C`). On average, 89.3% of aligned fragments were successfully assigned to 7,766 genes per sample (total genes in annotation: 8275). The resulting count matrix was used as input for differential expression analysis with DESeq2.

---

*Generated: 2026-01-28*
*Run directory: `03_quant_260128_v1`*
