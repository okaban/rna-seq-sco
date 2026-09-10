# Nanopore methylome and strand-resolved transcriptome of *Streptomyces coelicolor* A3(2) M145

Analysis code for the manuscript

> *A dual N4-methylcytosine/N6-methyladenine methylome of Streptomyces coelicolor:
> base-resolution reassignment of a cytosine mark and a same-strand AAGCCCG duplication*

(under revision at *Nucleic Acids Research*).

Nine biological samples across three time points of the developmental cycle
(T1 = 12 h, T2 = 24 h, T3 = 50 h; sample IDs `1-1 1-2 1-3 / 2-1 2-3 2-4 /
3-2 3-3 3-4`, the first digit being the time point). Genomic DNA and total RNA
were extracted from the **same harvested culture** for each sample, so the
methylome and the transcriptome derive from one biological specimen.

## Data availability

Sequencing data are deposited under BioProject `TODO_ACCESSION` (SRA):

| assay | platform | n | files |
|---|---|---|---|
| methylome | Oxford Nanopore MinION, FLO-MIN114 (R10.4.1), SQK-RBK114-96 | 9 | modification-tagged BAM (MM/ML: 6mA, 4mC, 5mC) aligned to `NC_003888.3`, plus raw pod5 signal |
| transcriptome | Illumina NovaSeq X Plus, PE150 | 9 | raw FASTQ (Ribo-Zero Plus depletion, NEBNext Ultra II Directional) |

Sequencing data are **not** in this repository. Reference genome:
`GCF_000203835.1` / `NC_003888.3`.

## Repository layout

Numbered directories follow the order of the analysis. Each holds
`scripts/` (code), `analysis/` (dated output) and sometimes `data/`.

| | |
|---|---|
| `01_qc` – `04_deseq2` | RNA-seq: fastp QC, HISAT2 alignment, featureCounts, DESeq2 |
| `05_annotation` – `10_SARP_motif_scan` | annotation, BGC dynamics, regulator networks, motif scans |
| `11_epigenome_integration` | methylome × transcriptome integration — the core of the paper |
| `12_supplementary_figures`, `15_paper_figures` | figure generation |
| `13_TF_binding-site`, `14_DMG_functional_enrichment` | binding sites, enrichment |
| `reports/` | per-analysis narrative reports |

## Reproducing

Conda environments are exported per stage, e.g.

```
conda env create -f 01_qc/scripts/environment_rnaseq_260127.yml
```

**Scripts hard-code absolute paths under the original author's home
directory.** Roughly 200 files contain such paths; they must be adapted before
the code will run elsewhere. Absolute paths appearing inside *tool output*
(MultiQC, featureCounts, sample tables) are left as-is deliberately: they
record what was actually executed, and rewriting them would falsify the
analysis record.

Large regenerable artifacts are excluded by `.gitignore` — sequencing data,
BLAST/aligner indexes, downloaded NCBI genomes, per-read co-modification
tables. All are reproducible from the deposited data with the scripts here.

## Retracted analysis

`15_paper_figures/scripts/27_figure4_shielded_exposed.py` is **disabled by a
guard and must not be run.** It computed a circular ROC — scoring
`y_true = is_exposed` against the predictor `nearest_methyl_distance`, where
`is_exposed` is itself defined by thresholding that same distance at 293 bp, so
the classifier recovered its own labelling rule. The resulting AUC values
(0.908, 0.917) were withdrawn from the manuscript. The file is retained with
its history and the reason stated in its docstring.

## License

MIT — see [LICENSE](LICENSE). The manuscript itself is not covered by this
license.
