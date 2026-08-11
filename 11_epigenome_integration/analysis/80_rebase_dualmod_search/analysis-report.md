# D2 — REBASE-wide scan for single-motif 4mC+6mA dual modification

**Question.** How unusual, across all catalogued bacterial/archaeal R-M systems,
is the AAGCCCG configuration: **both 4mC and 6mA on the same strand of one short
contiguous motif**?

**Data.** REBASE v602 flat file (`data/rebase/bairoch.txt`), 17,500 enzyme
records. Per record: organism (OS), enzyme type (ET), recognition motif (RS),
methylation type + strand sign (MS). 3,384 distinct recognition motifs carry ≥1
methylation annotation.

**Method.** `scripts/rebase_dualmod_scan.py`. Mod types parsed as substrings of
the MS label (REBASE writes N4-methyl-C as `Nm4C`; an earlier anchored regex
silently dropped all m4C — fixed). Strand from MS position sign (negative =
bottom strand). Two views: (1) single-enzyme dual specificity; (2) per-motif
dual (a motif with both an m4C and an m6A enzyme).

**Result.**
- Single-enzyme combined 4mC+6mA specificity: **51 enzymes**. **50/51 place the
  two marks on opposite strands** of the site (the Morgan 2016 / Furmanek-Błaszk
  2009 configuration). The **one** "same-strand" case (M2.Pae95I) is a gapped
  bipartite Type I motif (CCCNNNNNGTAG), not a short contiguous motif.
- Same-strand AND contiguous (no-N) AND ≤8 bp: **0**.
- Per-motif dual (across different enzymes): 55/3,384 motifs (1.63%); these are
  mostly palindromes acted on by separate m4C and m6A systems on opposite
  strands.
- **AAGCCCG is not in REBASE at all.**

**Interpretation.** Same-strand dual modification of a short contiguous motif —
the AAGCCCG configuration — is unrepresented in REBASE across all genera. This
quantitatively backs the novelty claim (previously supported only qualitatively
by Blow 2016 / Morgan 2016 / Harris & Goldman 2020).

**Caveat.** REBASE MS annotations reflect experimentally-characterised systems;
absence from REBASE means "not catalogued", not "cannot exist". The claim is
framed accordingly ("to our knowledge / not catalogued").

Tables: `tables/single_enzyme_dual.tsv`, `tables/per_motif_dual.tsv`.
Manuscript: Discussion §dual modification; Supplementary Table 10.
