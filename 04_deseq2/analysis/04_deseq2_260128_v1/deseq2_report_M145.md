# DESeq2 Report -- M145 RNA-seq

## 1. Project Overview

- **Organism**: *Streptomyces coelicolor* A3(2) M145
- **Experiment**: 3 conditions (M145_1, M145_2, M145_3 = timepoints) x 3 biological replicates = 9 samples
- **Sequencing**: Paired-end 150 bp (Illumina)
- **Analysis tool**: DESeq2 v1.46.0 (R v4.4.2)
- **LFC shrinkage**: apeglm v1.28.0
- **Reference genome**: GCF_000203835.1 (*S. coelicolor* A3(2), ~8.67 Mb)

---

## 2. Input Data

- **Count matrix**: `featureCounts_M145.txt` (featureCounts v2.1.1, `-t gene -g gene_id -s 0`)
- **Raw gene count**: 8,275 genes x 9 samples
- **After filtering** (>= 10 counts in >= 3 samples): **7,646 genes** retained for DE analysis

---

## 3. DESeq2 Analysis Settings

- **Design formula**: `~ condition` (M145_1 as reference level)
- **Normalization**: DESeq2 median-of-ratios (size factors)
- **Dispersion estimation**: DESeq2 default (gene-wise, trend, final MAP)
- **LFC shrinkage**: apeglm (`coef`-based) for M145_2 vs M145_1 and M145_3 vs M145_1; for M145_3 vs M145_2, reference level was releveled to M145_2 and Wald test was re-run before apeglm shrinkage

### Size Factors

| Sample | Size Factor |
|--------|------------|
| M145_1_1 | 0.918 |
| M145_1_2 | 0.859 |
| M145_1_3 | 1.120 |
| M145_2_1 | 0.815 |
| M145_2_3 | 0.709 |
| M145_2_4 | 0.804 |
| M145_3_2 | 1.718 |
| M145_3_3 | 0.916 |
| M145_3_4 | 1.444 |

M145_3 samples have larger size factors (0.92--1.72), reflecting the higher absolute counts in those libraries due to upregulation of highly expressed genes.

---

## 4. QC Results

### 4.1 PCA

- **PC1**: 83% variance -- clearly separates conditions along the time axis (M145_1 -> M145_2 -> M145_3)
- **PC2**: 16% variance -- captures within-condition variability
- All three biological replicates within each condition cluster tightly together
- The three conditions are well-separated, indicating a strong transcriptional shift across timepoints

### 4.2 Sample Distance Heatmap

- Replicates within the same condition show the shortest Euclidean distances
- M145_1 and M145_2 are more similar to each other than either is to M145_3
- Hierarchical clustering groups the three conditions into distinct clusters

---

## 5. Differential Expression Results

### 5.1 Summary (padj < 0.05)

| Contrast | Total Tested | Up | Down | Up (|LFC|>1) | Down (|LFC|>1) |
|----------|-------------|-----|------|--------------|----------------|
| M145_2 vs M145_1 | 7,646 | 2,563 (34%) | 2,697 (35%) | 2,003 | 1,845 |
| M145_3 vs M145_1 | 7,646 | 3,290 (43%) | 2,898 (38%) | 2,766 | 2,075 |
| M145_3 vs M145_2 | 7,646 | 2,639 (35%) | 2,904 (38%) | 2,014 | 1,493 |

### Key Observations

- The majority of genes (69--81%) show significant differential expression at padj < 0.05, consistent with large-scale transcriptional reprogramming across *S. coelicolor* growth phases
- The M145_3 vs M145_1 comparison shows the most DE genes (6,188; 81%), reflecting the cumulative transcriptional changes from early to late timepoints
- Among significant genes, approximately half show |log2FC| > 1, indicating substantial fold-change magnitudes

### 5.2 Top DE Genes

**M145_2 vs M145_1** (top by padj):

| gene_id | baseMean | log2FC | padj |
|---------|----------|--------|------|
| SC_RS13950 | 4,924 | +4.74 | 4.4e-286 |
| SC_RS32975 | 2,015 | +6.59 | 3.5e-250 |
| SC_RS11870 | 24,155 | +12.50 | 7.6e-248 |
| SC_RS26900 | 9,708 | +8.67 | 9.5e-241 |
| SC_RS10140 | 8,028 | +8.75 | 1.5e-235 |

**M145_3 vs M145_1** (top by padj):

| gene_id | baseMean | log2FC | padj |
|---------|----------|--------|------|
| SC_RS13760 | 1,314 | -5.64 | 6.4e-247 |
| SC_RS27535 | 4,664 | +8.95 | 1.6e-215 |
| SC_RS14645 | 6,832 | +7.66 | 1.8e-197 |
| SC_RS27575 | 7,441 | +7.65 | 2.0e-189 |
| SC_RS18270 | 5,357 | +5.58 | 6.4e-187 |

**M145_3 vs M145_2** (top by padj):

| gene_id | baseMean | log2FC | padj |
|---------|----------|--------|------|
| SC_RS13950 | 4,924 | -4.25 | 2.1e-239 |
| SC_RS27575 | 7,441 | +9.18 | 2.9e-217 |
| SC_RS27580 | 8,008 | +9.49 | 2.2e-212 |
| SC_RS27535 | 4,664 | +7.85 | 3.7e-192 |
| SC_RS27520 | 4,495 | +9.14 | 1.3e-185 |

### 5.3 Notable Patterns

- **SC_RS13950**: Strongly upregulated from M145_1 to M145_2 (+4.74 LFC) then downregulated from M145_2 to M145_3 (-4.25 LFC), showing peak expression at M145_2
- **SC_RS275xx cluster** (SC_RS27520, SC_RS27530, SC_RS27535, SC_RS27560, SC_RS27575, SC_RS27580, SC_RS27590): A cluster of adjacent genes consistently and strongly upregulated in M145_3 vs both M145_1 and M145_2 (LFC +7 to +9), suggesting activation of a biosynthetic gene cluster in the late growth phase
- **SC_RS11870**: Extremely high fold change (+12.5 LFC in M145_2 vs M145_1, baseMean 24,155), suggesting a gene that transitions from near-zero to very high expression

---

## 6. Biological Context

The large-scale differential expression observed across timepoints is consistent with known *S. coelicolor* biology:

- *Streptomyces* undergoes dramatic transcriptional reprogramming during its developmental life cycle, transitioning from vegetative growth to secondary metabolism and sporulation
- The progressive increase in Top-10 gene fraction observed in quantification (21% at M145_1, 30% at M145_2, 36% at M145_3) is consistent with the activation of highly expressed secondary metabolite biosynthetic gene clusters at later timepoints
- The SC_RS275xx gene cluster showing strong late-phase activation is located in a genomic region characteristic of *S. coelicolor* secondary metabolite biosynthetic clusters
- The 81% DE rate (M145_3 vs M145_1) reflects a near-global transcriptional shift, consistent with the known complexity of *Streptomyces* developmental programs

---

## 7. Output Files

### Results

- `results/normalized_counts_M145.tsv` -- DESeq2-normalized counts (7,646 genes x 9 samples)
- `results/DESeq2_M145_2_vs_1.tsv` -- M145_2 vs M145_1 DE results (apeglm-shrunk LFC)
- `results/DESeq2_M145_3_vs_1.tsv` -- M145_3 vs M145_1 DE results (apeglm-shrunk LFC)
- `results/DESeq2_M145_3_vs_2.tsv` -- M145_3 vs M145_2 DE results (apeglm-shrunk LFC)
- `results/DE_summary_M145.tsv` -- Summary counts of DE genes per contrast

### Figures

- `figures/PCA_M145.pdf/.svg` -- PCA plot (PC1 vs PC2, colored by condition)
- `figures/sample_distance_heatmap_M145.pdf/.svg` -- Euclidean distance heatmap
- `figures/volcano_M145_2_vs_1.pdf/.svg` -- Volcano plot (M145_2 vs M145_1)
- `figures/volcano_M145_3_vs_1.pdf/.svg` -- Volcano plot (M145_3 vs M145_1)
- `figures/volcano_M145_3_vs_2.pdf/.svg` -- Volcano plot (M145_3 vs M145_2)
- `figures/topVarGenes_heatmap_M145.pdf/.svg` -- Top 100 variable genes heatmap (rlog, row-centered)

### R Objects

- `rds/dds_raw.rds` -- DESeqDataSet after `DESeq()` (for downstream re-analysis)
- `rds/rld.rds` -- rlog-transformed data (for visualization)

### Logs

- `logs/deseq2_pipeline.log` -- Full pipeline execution log

---

## 8. Materials & Methods

Gene-level read counts from featureCounts v2.1.1 (8,275 genes, unstranded, `-t gene -g gene_id`) were analyzed for differential expression using DESeq2 v1.46.0 (R v4.4.2) with a `~ condition` design (M145_1, M145_2, M145_3 representing growth-phase timepoints). Genes with fewer than 10 counts in at least 3 samples were removed, retaining 7,646 genes. Library size normalization used DESeq2's median-of-ratios method. Log2 fold-change estimates were shrunk using the apeglm method (v1.28.0) via the `coef` parameter for direct coefficients, and via reference-level releveling for the M145_3 vs M145_2 contrast. Differentially expressed genes were identified at an adjusted p-value threshold of 0.05 (Benjamini-Hochberg). Sample-level QC was performed using PCA and Euclidean distance heatmaps on rlog-transformed counts.

---

*Generated: 2026-01-28*
*Run directory: `04_deseq2_260128_v1`*
