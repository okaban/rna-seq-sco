# D1 — Public-data integration inventory (stock-take + feasible-analysis proposal)

Date: 2026-06-21. Purpose: confirm every feasible cross-dataset integration of
the epi-trans methylome/transcriptome has been exploited, and scope the rest.

## A. Integrated and in the manuscript

| Dataset | Type | Integration | Manuscript |
|---|---|---|---|
| **Deng et al. 2023 (PNAS) Hi-C / FIRE** | 3D genome | methylome×FIRE genome-wide (ρ=0.32, partial r=0.18, spatial-perm p=0.008); Exposed vs Shielded promoter FIRE; HCR-M1–10 integration loci; co-high hotspots + nomination resolution | Discussion §FIRE; Supp Fig 8; **Supp Fig 12** |
| **Deng 2023 compartments A/B** | 3D genome | methylome compartment occupancy (62% A at T1 → 4% at T2; OR=1.96/3.24) | Discussion §FIRE |
| **REBASE v602** | R-M catalogue | AAGCCCG absence (0/82 Streptomyces); genus O/E; cross-genus dual-mod survey | Discussion; Supp Tables 8/10; Methods |
| **833 Streptomyces RefSeq genomes** | comparative genomics | AAGCCCG/GCCGGC O/E avoidance, SC_RS17645 homolog conservation (590/833) | Discussion; Fig 9; Supp Table 5 |
| **Own RNA-seq (DESeq2)** | transcriptome | methylation×expression (permissive r=−0.09); Exposed-TF synchronous demethylation; BGC activation | Results/Discussion (core) |
| **antiSMASH BGC predictions** | annotation | GCCGGC enrichment = composition, not regulation; per-cluster timing | Discussion §BGC; Supp Fig 7 |

## B. Feasible from local data — already exhausted

All locally-available orthogonal datasets are integrated. No further
local-only integration analysis remains: the FIRE/Hi-C cross-reference
(previously blocked in the 2026-05-06 report for lack of FIRE files) was
completed once Deng Dataset S4/Table S1 were obtained.

## C. Future directions — require external data not in this environment
(network-restricted; PNAS/NCBI/GEO not allowlisted)

1. **Second, independent Hi-C** (e.g. Lioy/Szafran *Streptomyces* 3C) — would
   test whether the methylome–FIRE coupling is robust across Hi-C datasets and
   growth conditions. *Strengthens, not required.*
2. **Regulator ChIP-seq** (BldD, HrdB, GlnR, WhiA, NAPs HupA/HupS) — would test
   whether protection-zone promoters coincide with TF/NAP occupancy (mechanism
   for the protection-zone architecture, currently a stated open question).
3. **Same-sample methylome + Hi-C time course** — the only way to establish
   temporal co-dynamics (vs the current static-compartment alignment); already
   named as the decisive future experiment in the Discussion.
4. **Cross-genus single-molecule methylomes** (PacBio/ONT for other
   actinomycetes) — direct test of AAGCCCG same-strand dual-mod conservation
   beyond the REBASE-catalogue argument.

## Conclusion
The manuscript exploits every orthogonal dataset available locally. Items in
(C) are genuine future directions requiring external data acquisition and (for
3) new wet-lab experiments; they are appropriately framed as future work in the
Discussion rather than gaps in the present analysis.
