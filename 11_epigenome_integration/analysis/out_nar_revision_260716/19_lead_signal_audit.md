# 19 — Lead-Signal Audit: is a main-role analysis un/under-used?

**Date:** 2026-07-16
**Scope:** Adversarial scan of the ~90 numbered analysis directories under
`11_epigenome_integration/analysis/` for signals that are (i) large-effect ×
small-p × adequate-n × novel, AND (ii) currently unused / under-used in the
manuscript, AND that clear **both** hard gates:

- **Gate A — no retracted-object dependency** (no 0.908/0.917 distance AUC, no
  "57 Exposed", no two-antagonistic-blocs / 2-program model, no eigengene
  r=−0.94/ρ=−0.995, no four-layer "Gatekeeper switch", no dissolution /
  cross-timepoint transfer).
- **Gate B — current-canon consistency** (62 Exposed / 989 Shielded / 1,051
  TSS-assignable; permissive geography-controlled partial Spearman r=−0.09).
  Values computed under old canon (57 / 955 / 998 / 1,055; AUC 0.908·0.917) are
  **invalid** and never enter the candidate table regardless of effect size.

**Every number below is quoted from a named result file** (path given). No value
is from memory.

---

## 0. Bottom line (read this first)

**No new *independent* lead-level finding survives both gates.** The five known
lead candidates already cover every current-canon signal that clears Gate A+B.
What the current-canon batch (dirs 62–85, 2026-04→06) actually produced is:

1. **Strengthening evidence for existing lead #3** (same-strand dual 4mC/6mA at
   AAGCCCG) — dirs 78, 80. These *also retire* the inflated `OR=138,440`.
2. **A confound-kill of the FIRE/methylome integration proxy** — dirs 76→83.
   This is eye-catching (ρ up to 0.66) but collapses to an oriC-distance
   artefact. It must be **demoted, not promoted**.
3. **A position clarification that caps the dual-KEGG claim** — dir 62. The
   siderophore enrichment is real but **CDS-internal**, so it cannot be read as
   transcriptional regulation.
4. **Supporting characterization of the Exposed set** — dir 82 (MerR/LysR family,
   core position). Current-canon, FDR-significant, but modest and descriptive.
5. **One genuinely UNEXECUTED high-value analysis** — dir 85 G1/G2 (Deng CID
   boundary alignment / dynamic co-change). Data is in hand locally; not yet run.
   **This is the only place a new lead could still emerge**, and it is the honest
   answer to the concern that "a main-role analysis may not yet have been done."

**Conclusion on the 5 known candidates: no addition. They are exhaustive for the
current-canon evidence base.** The residual risk is not an *un-analyzed number
sitting on disk*; it is *one un-run analysis* (dir 85 CID boundaries).

---

## 1. Summary table — signals evaluated against both gates

| # | Signal (dir) | Effect | p / padj | n | Source file | Novel vs known-5? | Claimable range | Retraction risk (Gate A) | Canon (Gate B) | Recommendation |
|---|---|---|---|---|---|---|---|---|---|---|
| S1 | Dual-mod motif specificity (78, P1-1) | near/far 365× within AAGCCCG; **430×** vs generic; 19% per-read co-occ | far-rate 0.00013 vs near 0.0474 | 18,117 AAGCCCG 6mA vs 2.10 M generic | `78_reviewer_robustness/tables/P1_1_crosstalk_rates.tsv` | **Supports #3, not new** | Descriptive (basecaller QC) | ✅ none | ✅ current | **Fold into lead #3**; use to replace OR=138,440 |
| S2 | REBASE-wide same-strand dual scan (80, D2) | same-strand+contiguous+≤8bp = **0 / 17,500 enzymes**; AAGCCCG absent | — (census) | 17,500 enzymes / 3,384 motifs | `80_rebase_dualmod_search/tables/single_enzyme_dual.tsv` (report) | **Supports #3, not new** | "not catalogued" (correctly framed) | ✅ none | ✅ current | **Fold into lead #3** (novelty backing) |
| S3 | Dual-targeted siderophore KEGG (62) | sco00975 **OR=39.5**, enrichment 11.9× | padj=**9.9e-4** | 5/7 pathway genes | `62_GO_KEGG_enrichment/tables/E1_KEGG_enrichment_Dual-targeted.tsv` | = known #2 | **Descriptive only — CDS-internal** (see §5) | ✅ none | ✅ current | **Supp only; strip any transcription-regulation implication** |
| S4 | Exposed TF-family / position (82, C6) | MerR **OR=4.71**, LysR **OR=3.47**, core OR=3.87 | p_BH 0.022 / 0.035 / 3e-4 | 62 vs 989 | `82_exposed_shielded_metric_scan/tables/exposed_vs_shielded_metric_scan.tsv` | Characterizes canon set; not a new lead | Descriptive | ✅ none | ✅ current (62/989) | **One supp sentence** (permissive-consistent) |
| S5 | Methylome × FIRE integration (76) | ρ=0.319 genome-wide; **ρ=0.661** at Deng integration loci | 5.6e-42 / 0.038 | 1,724 bins / 10 loci | `76_FIRE_methylation_crossref/tables/*` | Would be new IF real | **KILLED** — see S6 | ✅ none | ✅ current | **Demote existing FIRE sub-claim; do NOT promote** |
| S6 | oriC-confound test of FIRE (83) | partial ρ **0.026** (from 0.319); Exposed–FIRE gap n.s. in all 4 oriC quartiles | Q1 0.42 / Q2 0.13 / Q3 0.055 / Q4 0.46 | 1,724 / 62 vs 985 | `83_oriC_confound_FIRE/tables/oric_confound_results.tsv`, `within_compartment.tsv` | **Negative result** (kills S5 as independent coupling) | Falsification of FIRE-proxy sub-claim | ✅ none | ✅ current | **Use to demote FIRE engineering claim; permissive-consistent** |
| S7 | Deng CID-boundary alignment (85 G1/G2) | **not computed** | — | — (Deng sd02 local) | `85_deng_integration_gap_analysis/gap_analysis.md` | **Potential new lead — unexecuted** | Unknown until run | ✅ none (proposed clean) | n/a | **RUN IT** — only remaining lead-level gap |

Effect sizes for S1–S6 read verbatim from the files above on 2026-07-16.

---

## 2. Which known lead candidates are confirmed vs. refined

- **#1 Developmental spatial redistribution** (core→arm; GCCGGC 4mC
  T1:1289→T2:407→T3:21; Jaccard=0.000). Independently reproduced in current-canon
  dir 79 (`C4_region_by_timepoint`: counts "match the canonical numbers …
  GCCGGC 1289/407/21") and dir 63. **Confirmed, current-canon.** Note the raw
  temporal counts (1289/407/21) originate in old-canon dir 37 (H14) but are
  **re-quoted unchanged** by current-canon dirs 79/57 — the site-count vector is
  canon-invariant (it is a per-motif census, not an Exposed/Shielded partition),
  so it is safe. Foreground unchanged.
- **#2 dual-mod KEGG enrichment** (siderophore sco00975 OR≈39.5 padj=0.001).
  **Confirmed numerically** but see §5 — the hits are CDS-internal → keep
  descriptive, no transcription-regulation reading.
- **#3 Two novelties** (4mC reassignment; same-strand dual 4mC/6mA at AAGCCCG,
  REBASE-absent). **Strengthened** by S1+S2. **Action:** replace the retracted-
  adjacent `OR=138,440` with motif-specificity 430× + per-read 19% + REBASE
  0/17,500 (all current-canon).
- **#4 Protection zone (±300 bp promoter avoidance).** Consistent with dir 45
  (sequence-level, current-canon) and dir 64. **Confirmed** — see §5.
- **#5 genus O/E contrast** (AAGCCCG≈0.65–0.70 vs GCCGGC≈1.0). Sits in the
  genus-conservation analyses; not re-audited here for new signal (no new
  current-canon table supersedes it). **Retained.**

---

## 3. Detailed evidence — the strengthening items (S1, S2 → lead #3)

**S1 — dual modification is motif-specific, not k-mer cross-talk**
(`78_reviewer_robustness/tables/P1_1_crosstalk_rates.tsv`):

| 6mA class | n | near 4mC rate | far rate | near/far |
|---|---|---|---|---|
| AAGCCCG | 18,117 | 0.047442 | 0.000130 | **364.99×** |
| non-AAGCCCG generic | 2,101,211 | 0.000110 | 0.000101 | 1.09× |

Motif specificity = 0.047442 / 0.000110 = **430×**; per-read co-occurrence
`frac_6mA_with_near_4mC` = **0.1897** (19%). Interpretation (from the report):
generic cross-talk would inflate near-4mC around *all* 6mA (it does not: near/far
1.09 across 2.1 M controls), so the short-range 4mC excess is confined to
AAGCCCG → genuine, motif-specific co-deposition. **The report itself states
`OR=138,440` is a saturated-2×2 artefact and must be dropped.** Claimable without
KO/biochem: yes (basecaller-level statistical argument, framed as
co-deposition/co-occurrence, not mechanism).

**S2 — REBASE-wide rarity** (`80_rebase_dualmod_search/analysis-report.md`,
tables `single_enzyme_dual.tsv` / `per_motif_dual.tsv`): across 17,500 enzyme
records / 3,384 methylated motifs, single-enzyme same-strand + contiguous + ≤8 bp
dual = **0**; the one same-strand case (M2.Pae95I) is a gapped bipartite Type I
motif, not contiguous; AAGCCCG is absent from REBASE entirely. Correctly framed
as "not catalogued," not "cannot exist." Claimable: yes.

Both are **support for existing novelty #3(b)**, not independent leads.

---

## 4. The FIRE proxy — eye-catching but a confound (S5 killed by S6)

`76_FIRE_methylation_crossref` reports methylome density tracking Deng FIRE:
genome-wide Spearman ρ=0.319 (p=5.6e-42, n=1,724) and, at Deng's 10 validated
HCR integration loci, ρ=0.661 (p=0.038). Taken alone this reads as a
main-figure "Hi-C-free integration-site proxy."

`83_oriC_confound_FIRE/tables/oric_confound_results.tsv` and
`within_compartment.tsv` dismantle it under current canon (seed 42):

- FIRE is almost entirely an oriC-distance gradient: ρ(FIRE, dist_oriC) = **−0.934**.
- Partial Spearman(meth, FIRE | dist_oriC) = **0.026** (collapses from 0.319).
- Within-core partial (| oriC-dist) = **0.024**.
- Exposed>Shielded FIRE overall p=8.6e-6, but **not significant in any
  oriC-distance quartile** (Q1 p=0.42, Q2 0.13, Q3 0.055, Q4 0.46).

**Verdict:** the standalone methylome→FIRE predictive/engineering sub-claim is
**not supported as an independent coupling** — methylation, FIRE, and compartment
all track one oriC-centred core architecture. This is **consistent with the
permissive/spatial thesis** but undercuts the engineering proxy. **Do NOT promote
S5; demote the existing FIRE sub-claim** (Discussion 3D para / Supp Fig 8) to a
shared-gradient statement. This is a null/negative result that *supports* the
permissive framing — log it as such, do not dress it as a finding.

---

## 5. Position stratification — why the dual-KEGG stays descriptive (S3)

Per the standing rule that most marks are CDS-internal and gene-body enrichment
must not be read as transcriptional control:

- **Site-position census** (`62_GO_KEGG_enrichment/tables/F1_classified_meth_sites.tsv`,
  n=2,737 classified): CDS_internal 2,293 / 5UTR_approx 183 / promoter 172 /
  intergenic 89. i.e. **~84% CDS-internal** — matches dir 79's 86/84/86% (GCCGGC).
- **The 5 siderophore sco00975 genes are CDS-internal.** The two sites mapping to
  the pathway genes are `SC_RS15995 GCCGGC CDS_internal dist_from_tss +1457` and
  `SC_RS16000 AAGCCCG CDS_internal +1037/+1038`. **No promoter-proximal dual site
  on the siderophore cluster.** → the sco00975 OR=39.5 is a **gene-body**
  enrichment; report as descriptive co-localization, **not** epigenetic
  transcriptional regulation.
- **Position-stratified KEGG** (`F1_stratified_KEGG_by_position_v2.tsv`): the only
  FDR-passing pathway terms are also CDS-internal (GCCGGC glyoxylate padj=0.0038,
  TCA padj=0.0099; AAGCCCG aminoacyl-tRNA padj=0.0202). **Promoter-stratified
  terms do not survive FDR** (best GCCGGC promoter Two-component padj=0.094; best
  AAGCCCG 5UTR purine padj=0.0002 is the sole exception and sits in 5′UTR, not
  promoter). → No robust *promoter* functional enrichment exists to promote.
- **Protection layer (#4) confirmed at sequence level, region-stratified**
  (`45_sequence_level_motif_depletion/tables/stratified_by_region.tsv`): for
  regulatory genes, both motifs are depleted in **gene_body** (TGGCCGGC core
  Fisher OR=0.586 p=3e-6; AAGCCCG core OR=0.548 p=1.4e-5) but **neutral in the
  promoter** (TGGCCGGC core OR=0.889 p=0.43; AAGCCCG core OR=1.02 p=0.88). Correct
  reading: sequence-level counter-selection acts in the body, not the promoter;
  promoter protection is occupancy-based, not sequence-encoded.

**Any enrichment-type lead must therefore stay CDS-internal-labelled and
descriptive.** None qualifies for main-figure promotion on a regulatory reading.

---

## 6. Supporting characterization (S4) — current-canon, modest

`82_exposed_shielded_metric_scan/tables/exposed_vs_shielded_metric_scan.tsv`
(62 Exposed vs 989 Shielded, seed 42, BH-FDR). FDR-significant metrics:

- `position_core` OR=**3.87** (p_BH=3e-4) — 53/62 Exposed core.
- `dist_to_chrom_center` Cliff δ=**−0.316** (p_BH=3e-4) — Exposed closer to centre.
- `family_MerR` OR=**4.71** (p_BH=0.022, 6/62); `family_LysR` OR=**3.47**
  (p_BH=0.035, 7/62).
- `signed_LFC_T2vsT1` Cliff δ=**−0.193** (p_BH=0.038) — the report notes this
  **is** the same weak repression-direction bias as the canonical r=−0.09.
- **Not significant:** baseMean, |LFC| magnitude/variability (T2 & T3), DEG
  fraction (T2 & T3). i.e. the partition is positional/family-distinct, **not**
  transcriptionally distinct by magnitude → permissive, as reported.

Also cross-checked in current-canon dir 78: `P1_3b_permissive_spatial.tsv` gives
partial ρ=**0.128**, block-perm p=**0.0021**, n=**1,019** ("SURVIVES (weak bias is
real)"), and `P1_4_equivalence.tsv` shows |LFC| Exposed≈Shielded (Cliff δ −0.010 /
+0.091; TOST-equivalent) — consistent with the locked r=−0.09 permissive result.
**Recommendation:** at most one supplementary sentence on MerR/LysR/core position;
frame strictly as characterization of the (canon) Exposed set, permissive-consistent.

---

## 7. The one un-run analysis — dir 85 (the real answer to the concern)

`85_deng_integration_gap_analysis/gap_analysis.md` (2026-06-29) inventories every
Deng-2023 integration angle and flags two **feasible-now, HIGH-value,
never-executed** analyses using local data (Deng sd02):

- **G1 — CID (interaction-domain) boundary alignment**: do GCCGGC/AAGCCCG sites,
  the protection zone, or Exposed promoters coincide-with / avoid Deng's CID
  boundaries? *"Paper itself names this as future work but data is in hand … a
  structural finding that can REPLACE the dead FIRE-density proxy."*
- **G2 — dynamic co-change**: Δmethylation(T1→T2) vs Δcompartment(PC_M→PC_L) and
  ΔFIRE(M→L). Quantifies the central "methylation relocation mirrors 3D refolding"
  claim, which is currently only a static overlap.

**This is the substantive gap.** It is not a number sitting unused on disk; it is
an analysis not yet performed, on data already local, that could (a) supply a
structural epigenome×3D finding to replace the collapsed FIRE proxy, and (b)
convert the "mirrors refolding" statement from static to quantitative. If any new
main-role result exists to be found, it is here. **Recommend running G1 (and G2)
before submission**; both clear Gate A by construction (no retracted objects) and
must be reported at current canon.

---

## 8. Excluded — eye-catching but failed a gate (completeness ledger)

Listed so the audit is auditable. Each has a one-line reason.

| Signal | Dir | Number | Excluded because |
|---|---|---|---|
| Coordinated-regulator TSS proximity | 50 (H27) | median **114 bp vs 762 bp, p=5.3e-28** ("strongest signal in the project") | **Old canon** (n=57 / 998). Partition concept survives but is recomputed at 62/989; this specific value is stale — do not quote. |
| Exposed/Shielded classifier | 52 (H29) | 293 bp threshold **AUC=0.917**, 100% sens | **Gate A fail** — distance-classifier AUC is retracted (circular). Entire dir 52 + `ROC_analysis*.tsv` / `ROC_curves.*` treated as poisoned. |
| Exposed regulator counts | 51,54–60 (H28,H31–37) | n=**57** Exposed / **955**–**998** Shielded / **1,055** | **Old canon** throughout. Superseded by 62/989/1,051. |
| Two antagonistic blocs / phase | 57 (H34) | activation vs repression bloc | **Gate A fail** — 2-program/antagonistic-blocs model is retracted. |
| Four-layer Gatekeeper model | 34 (H11) | "3 active layers + 1 excluded" | **Gate A fail** — four-layer "Gatekeeper switch" retracted. |
| redZ paradox resolution | 26 | redZ +0.91/+0.97 | Author-flagged as role-ambiguous / on hold; not a methylation lead. |
| comod threshold ROC | 73_comod_threshold_ROC | ROC/AUC | **Gate A fail** — ROC/AUC family; excluded on sight. |
| Distance-AUC eigengene etc. | 34,52,57 | AUC 0.908/0.917, ρ=−0.995 | **Gate A fail** — retracted objects. |
| STREME unassigned-6mA motif | 74 | top motif CCGGTGCGCGG **P=0.023** | Weak; not a robust novel motif. Not a lead. |
| Unattributed 6mA (1,674 sites) | 66 | TSS-proximity **0.725 depleted** (p=1.1e-149); BGC n.s. | Real but permissive-consistent (avoids TSS); a limitation (most 6mA unassigned), not a headline. |

---

## 9. Coverage / completeness

- **Directories present:** ~90 numbered analysis dirs (101 top-level entries incl.
  non-numbered helpers). **Markdown summaries read for all that carry one.**
- **Result tables read verbatim for every current-canon decision-bearing dir:**
  62, 63, 64, 66, 74, 76, 78, 79, 80, 82, 83, 85.
- **Old-canon dirs (26–60, H1–H37)** were read at summary level and **excluded as
  lead sources by rule** (Gate B), not by table-by-table re-audit — the canon
  mismatch (57/955/998/1,055) disqualifies their effect sizes regardless.
- **Read-verbatim but not a lead (added on review):** `25_4mC_6mA_differential`
  — region-resolved correlation `tables/region_specific_correlation.tsv` gives 4mC
  **Promoter** ρ=0.1393 (p=8.3e-4, n=573), 4mC Gene Body ρ=0.0109 (p=0.66); 6mA
  null in both regions. This is the raw positive promoter correlation that the
  **geography-controlled permissive r=−0.09 supersedes** — it is the pre-reframe
  simple correlation, not an independent lead, and quoting it as positive would
  contradict the locked central thesis. Region composition table matches dir 79
  (promoter fraction small). **Excluded as a lead** (superseded framing); no Gate-A
  object involved.
  - `18_tss_analyses` (20 tables) — TSS geometry; feeds #4, no standalone lead expected.
  - `23_expanded_motif_search` (37 tables) — census/motif source, not a results dir.
  - `65_per_read_comod` (12 tables) — per-read co-mod; overlaps S1 (dual-mod);
    confirm it holds nothing beyond `78/P1_1` before finalizing #3 wording.
  - `67_spatial_integration`, `68_or_permutation`, `69_boundary_sensitivity`,
    `70_methylation_dynamics` — robustness/sensitivity supports; no new lead
    surfaced in headers, tables not fully opened.
  - `77_reviewer_figures`, `81_D1_integration_inventory` — figure/inventory, no data.
- **Judgment held / to confirm:** `65_per_read_comod` is the one current-canon dir
  I did not fully table-audit that could plausibly touch lead #3; recommend a quick
  confirm it adds nothing beyond the 430× / 19% already in `78/P1_1`.

---

## 10. Recommendation to the author

1. **Add no new main-role candidate.** The known 5 are exhaustive for
   current-canon evidence. Say so explicitly in the response letter if asked
   whether a lead analysis is missing — with one exception (point 4).
2. **Lead #3:** replace `OR=138,440` everywhere with motif-specificity **430×** +
   per-read **19%** + REBASE **0/17,500** (S1, S2). Keeps the novelty, drops the
   artefact.
3. **Lead #2 / dual-KEGG:** keep in supplement, **CDS-internal-labelled**, with no
   transcriptional-regulation implication (S3, §5).
4. **Run dir 85 G1 (CID boundaries) and G2 (dynamic co-change)** before
   submission — the only place a genuine new lead can still appear, data already
   local. If G1 yields a boundary-alignment signal, it can **replace** the
   demoted FIRE proxy as the epigenome×3D structural result.
5. **FIRE proxy:** demote the standalone engineering sub-claim to a shared
   oriC-gradient statement (S6). Report the collapse (0.319→0.026) as a
   permissive-consistent negative control, not a finding.
