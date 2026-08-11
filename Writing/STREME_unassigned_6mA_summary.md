# STREME on T1 unassigned 6mA sites

- Input: 1,497 high-confidence T1 6mA sites with `assigned_motif=unassigned` (from `analysis/23_expanded_motif_search/6mA_motif_assignment.csv`), ±50 bp windows extracted from `NC_003888.3` (101-bp sequences, strand-corrected). AAGCCCG-attributed (n=260) and other-known-motif (n=174) sites excluded.
- De novo top motif: **CGGCGACGRCGA** (STREME-1, E=3.0e-03, 353 sites) — high-GC novel candidate, not matching any annotated M145 MTase target.
- Motif 2 (CAAGCCCGCCV, E=1.4e-02, 118 sites) is residual AAGCCCG context — the ±3 bp in-motif test leaks AAGCCCG occurrences within ±15 bp; tighten the assignment window in future revisions.
- Motifs 3–5 (CTGCTCGCC, CGCCCYCRCC, CCSGGTTSCYG) reflect generic Streptomyces GC-rich background (E=2.2e-02 to 9.1e-02) and are not site-specific candidates.
- Discriminative run (unassigned vs n=260 AAGCCCG positives, `streme_unassigned_T1_discriminative/`) yields no significant motif (best E=0.12 for CCGGTGCGCGG); small negative set limits power. Outputs at `analysis/74_streme_unassigned_6mA/`.
