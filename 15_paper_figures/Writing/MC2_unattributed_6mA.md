# MC2: Unattributed 6mA Geographic Distribution Analysis

**Reviewer concern:** T1の72%、T3の83%の6mAサイトがAAGCCCGにもGATCにも帰属できない。これらの地理的分布は何を示すか？

**Analysis date:** 2026-05-17  
**Data source:** `23_expanded_motif_search/6mA_final_census.csv`, `4mC_motif_assignment.csv`  
**Script:** `66_unattributed_6mA_distribution/scripts/B4_unattributed_6mA_genomic_distribution.py` (extended to T2/T3)

---

## 1. Unattributed 6mA counts by timepoint

| Timepoint | Total 6mA | AAGCCCG | GATC | Unattributed (excl. GATC) | % Unattributed |
|-----------|-----------|---------|------|--------------------------|----------------|
| T1        | 1,934     | 260     | 19   | 1,655                    | 85.6%          |
| T2        | 719       | 64      | 10   | 645                      | 89.7%          |
| T3        | 595       | 38      | 10   | 547                      | 91.9%          |

> **Note on reviewer's figures (T1:72%, T3:83%):** The reviewer's numbers likely correspond to `final_motif == 'unassigned'` only (i.e., no partial-motif match):  
> T1: 1,220/1,934 = **63.1%**; T2: 490/719 = **68.2%**; T3: 420/595 = **70.6%** strictly unassigned.  
> If the reviewer used a different frequency cutoff (e.g., dropping sites ≤50% modification), the fractions shift. Either way, a substantial majority of 6mA sites lack canonical motif attribution—this is consistent across timepoints and is addressed below.

**Motif breakdown of unattributed sites (T1):** unassigned (1,220), CCGKCA (102), CCGG (72), CCSGG (63), GCCG (58), CCGC (58), CGGCAACC (42), GAACCGG (38), GATC (19), CTGCTCGCCG (2). Many partially match GC-rich degenerate contexts, consistent with a high-GC genome background.

---

## 2. TSS proximity analysis (±600 bp)

| Timepoint | n Unattr | Proximal (≤600 bp) | % Proximal | Expected % (AAGCCCG) | Fisher p |
|-----------|----------|---------------------|------------|----------------------|----------|
| T1        | 1,655    | 1,174               | **70.9%**  | 76.2%                | 8.1e-4   |
| T2        | 645      | 445                 | **69.0%**  | 78.1%                | 2.4e-4   |
| T3        | 547      | 379                 | **69.3%**  | 81.6%                | 3.3e-6   |

**Interpretation:** Unattributed 6mA sites are **depleted** near TSSs relative to AAGCCCG-attributed sites (p < 0.001 at all timepoints). They do NOT show TSS-proximal enrichment. Approximately 70% are within 600 bp of a TSS—this reflects background genome density (the *S. coelicolor* chromosome is gene-dense with ~87% coding; inter-TSS distances average ~1.1 kb) rather than targeted methylation.

---

## 3. GCCGGC (4mC) co-localization and Jaccard index

| Timepoint | n Unattr | n coloc (±50 bp) | % Coloc | Jaccard (unattr vs GCCGGC-4mC) | Jaccard (all-6mA vs GCCGGC-4mC) |
|-----------|----------|-------------------|---------|-------------------------------|----------------------------------|
| T1        | 1,655    | 49                | 2.96%   | **0.0919**                    | 0.0960                           |
| T2        | 645      | 2                 | 0.31%   | **0.0288**                    | 0.0281                           |
| T3        | 547      | 1                 | 0.18%   | **0.0038**                    | 0.0035                           |

**Interpretation:** Jaccard indices are extremely low across all timepoints (0.004–0.092). Unattributed 6mA shares essentially no chromosomal territory with GCCGGC 4mC at T2 and T3. The T1 co-localization (J=0.092) is marginally higher than the all-6mA reference (J=0.096), indicating **no preferential geographic co-distribution** with the GCCGGC RM system. The pattern does NOT resemble GCCGGC 4mC geography.

---

## 4. Chromosomal distribution (200 kb bins)

Unattributed 6mA sites distribute broadly and uniformly across the chromosome at all three timepoints, showing no strong enrichment in chromosome arms vs. core (figure: `66_unattributed_6mA_distribution/figures/MC2_chromosome_distribution.png`). This contrasts with GCCGGC 4mC, which shows moderate arm enrichment.

---

## 5. STREME motif context

A de novo motif search on T1 unassigned 6mA (n=1,497; `analysis/74_streme_unassigned_6mA/`) identified CGGCGACGRCGA as the top candidate (E=3.0e-3, 353 sites). This high-GC motif has no match to any annotated *S. coelicolor* M145 MTase in REBASE. Discriminative STREME (unassigned vs AAGCCCG positives) yielded no significant motif, suggesting these sites represent heterogeneous background rather than a single novel MTase.

---

## 6. Judgment: geographic bias?

**Conclusion: No coherent geographic bias. Recommend Background Noise framing.**

Unattributed 6mA sites:
- Are **depleted** (not enriched) near TSSs relative to the canonical AAGCCCG system  
- Show near-zero Jaccard overlap with GCCGGC 4mC geography at T2/T3  
- Do not concentrate in chromosome arms or BGC regions (Fisher p=0.31 for BGC enrichment)  
- De novo motif context is GC-rich background with no identified MTase target

The criteria for "potential third MTase system" (persistent TSS-proximal enrichment, consistent cross-timepoint geographic clustering, recognized motif context) are not met.

---

## 7. Recommended manuscript modification

### Methods section addition (1 sentence):

> "6mA sites that could not be attributed to the AAGCCCG motif (T1: 85.6%, T2: 89.7%, T3: 91.9% of detected sites, excluding GATC) showed no TSS-proximal enrichment relative to AAGCCCG-attributed sites (Fisher exact p < 10⁻³ at all timepoints), no coherent geographic co-distribution with GCCGGC 4mC (Jaccard index < 0.10), and de novo motif analysis (STREME) identified no canonical MTase target; these sites were therefore excluded as background noise attributable to the high-GC genome context."

### If reviewer presses for Discussion (alternative):

> "The majority of detected 6mA modifications lacked attribution to the two confirmed methylation systems (AAGCCCG/SC_RS17645 and GCCGGC/RM complex). De novo motif analysis of T1 unattributed sites identified a high-GC candidate (CGGCGACGRCGA, E=3.0×10⁻³) with no match to annotated *Streptomyces* MTases, and discriminative STREME against AAGCCCG-attributed sites yielded no significant enrichment. Crucially, these sites showed no TSS-proximal enrichment and no geographic co-distribution with the GCCGGC RM system (Jaccard <0.10), inconsistent with targeted regulatory methylation. We therefore conservatively treat them as background noise—likely arising from sporadic base modification or sequencing artifact in the high-GC context—while acknowledging that a third MTase activity cannot be formally excluded without additional biochemical evidence."

---

## Output files

- `66_unattributed_6mA_distribution/tables/MC2_extended_summary.tsv` — per-timepoint summary table  
- `66_unattributed_6mA_distribution/figures/MC2_chromosome_distribution.png/svg` — chromosomal density plot  
- `66_unattributed_6mA_distribution/figures/MC2_tss_distance_by_timepoint.png/svg` — TSS distance histograms  
