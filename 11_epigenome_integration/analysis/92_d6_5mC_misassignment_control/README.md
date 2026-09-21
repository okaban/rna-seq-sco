# 92 — In-silico 5mC→4mC misassignment control (D6), 2026-09-21

**Question.** Novelty (a) rests on the Dorado `4mC_5mC@v1` model calling the GCCGGC cytosine as 4mC
(5mC signal < 0.01%). Would a genuine 5mC at that position have been reported as 4mC?

**Ground truth.** ONT public modified-base validation set `s3://ont-open-data/modbase-validation_2024.10/`
(synthetic all-5-mer constructs, 32 references, 4,960 bp): `5mC_rep1.pod5` — all 256 annotated C sites are
5mC; `control_rep1.pod5` — all C canonical. Each file has 159,999 reads (32 × 5,000). Refs and BEDs copied here;
pod5 md5 in `pod5_md5.txt` (files not committed, 361/375 MB).

**Run.** MinKNOW-bundled dorado 1.1.1 (`/Applications/MinKNOW.app/Contents/Resources/dorado/bin/dorado`),
`dna_r10.4.1_e8.2_400bps_sup@v5.2.0` + `..._4mC_5mC@v1`, `--reference all_5mers.fa`, `--device cpu`
(Metal refused inside the sandbox), 3,000 reads per library sampled with seed 42 (`*.sub3000.txt`;
full files would take ~57 h on CPU at 0.64 s/read). Logs: `*.dorado.log`. `run_metal.sh` = GPU variant for Terminal.

**Analysis.** `d6_leakage.py`: forward-strand observations at annotated 5mC sites (5mC library) or any C (control),
MAPQ ≥ 20, MM/ML via pysam; 4mC code 21839, 5mC code m. Output `tables/SuppTable11_5mC_misassignment_control.tsv`.

| library | thr | n_obs | 4mC | frac_4mC | 5mC | frac_5mC |
|---|---|---|---|---|---|---|
| genuine 5mC | 0.5 | 18,010 | 162 | 0.90% | 16,724 | 92.9% |
| genuine 5mC | 0.9 | 18,010 | 12 | 0.07% | 14,313 | 79.5% |
| genuine 5mC, CCGG context | 0.5 | 334 | 0 | 0.0% | 328 | 98.2% |
| unmodified C | 0.5 | 8,925 | 420 | 4.7% | 38 | 0.43% |
| unmodified C | 0.9 | 8,925 | 41 | 0.46% | 12 | 0.13% |
| unmodified C, CCGG context | 0.5 | 51 | 4 | 7.8% | 1 | 2.0% |

**Reading.** A GCCGGC cytosine that was really 5mC would have been called 5mC (~93–98%), not 4mC (≤ 0.9%, 0/334 in
CCGG). The observed GCCGGC pattern (4mC high, 5mC < 0.01%) is therefore not a 5mC misassignment. Caveats: synthetic
templates, not S. coelicolor DNA; the per-read false-4mC rate on unmodified C is 4.7% at 0.5 (0.46% at 0.9) — relevant
to per-read analyses, not to the replicate-concordant site calls. Orthogonal chemical confirmation remains pending.
Worst leakage contexts (t=0.5): GACCA 14%, CGCTA 11%, GACTT 9% (n 57–89 each) — none resembles CCGGC.
