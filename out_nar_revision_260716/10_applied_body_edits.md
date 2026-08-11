# 10 · Applied body-text edits (EN manuscript)
File: Manuscript_EN_NAR_numbered_inline.md — edits applied in place (git-reversible; author commits).
Applied after author authorized the recommended approach. All numbers locked; canon respected.

## Investigations that unblocked judgments 5 & 6
- **J5 (transcriptome DEG enrichment, C16): already used correctly.** The manuscript (Results
  §R6, L234) uses the DEG enrichment exactly as the permissive-supporting contrast: differentially
  METHYLATED genes show 0 significant KEGG/GO/COG terms, "versus up to 18 for differentially
  EXPRESSED genes" (DEG enrichment report: dozens of significant pathways — ribosome/metabolism
  down, ABC-transport/siderophore up). No new figure needed; the contrast is the point. No change.
- **J6 (Le reference): manuscript is correct; the SURVEY doc had the error.** Manuscript ref 10 =
  "Le, T.B.K. & Laub, M.T. (2016) EMBO J. 35, 1582–1595", matching the survey §6 prose and its
  in-text use at L232 (CID-boundary sculpting by transcription). The survey reference-table row 32
  mis-tagged "Le (2013) PMID 24158908" as ref 10 — that is the survey's error, to fix in
  03_prokaryote_epigenetics_survey.md, NOT in the manuscript. No manuscript change.

## Edits applied (5)
1. **[J3/C8] Limitations header rename** — "External validation of Nanopore 4mC calling accuracy in
   high-GC contexts." → "Internal specificity control for Nanopore 4mC calling in high-GC contexts."
   (the section provides an internal control, not external validation; L283.)
2. **[J2/C14] Intro protection-zone width correction (factual)** — "(extending ~2,200 bp around gene
   transcription start sites, within which a 293 bp promoter-proximal window…)" → "(the all-gene TSS
   metagene is depleted within ~±300 bp of the transcription start site; a 293 bp promoter-proximal
   window…)". Matches the Fig 2B data (05 log); the ~2,200 bp figure was the regulatory-gene
   symmetric-zone quantity, mis-applied to the all-genes metagene.
3. **[J1/C1] Intro "three principal findings" → hierarchy** — one central finding (spatial organiser,
   relational verb) + two enabling observations (4mC reassignment enabling; dual mark independent).
   Follows 07 Proposal 2 + R2-EDIT-3 (dual mark not over-sold as load-bearing).
4. **[J1/C1] Abstract opening reframed** — foregrounds the one claim with a relational verb
   ("spatially organised relative to … rather than directing transcription"; R2-EDIT-2, guard F1),
   4mC + dual mark demoted to "two enabling observations". Follows 07 Proposal 1 (R2-revised).
5. **[J1/R2-EDIT-1] Abstract co-dynamics demoted** — the ρ=0.43 / n≈3 ΔΔ clause removed from the
   Abstract; replaced with the robust static co-localisation (Compartment-A OR=1.96, p=3.5×10⁻¹⁷).
   Temporal co-dynamics reduced to "suggested but not established … Discussion". The full ρ=0.43
   treatment remains in Discussion L210 unchanged.

## NOT changed (deliberately)
- Closing Abstract thesis "acts as a spatial organizer, not a direct regulator" — kept (title
  framing; "organizer" as noun contrasting "regulator" is the established, defensible phrasing).
  Author may relational-ise if desired.
- Discussion L210 ρ=0.43 paragraph — unchanged (that is where the hedge already lives).
- two-layer "on/off gate" (L233) — left for author (canon-adjacent; confirm intended vs retracted
  four-layer model).

## Still open (need author / reviewer material)
- Mutilka dataset identity (C16) — needs the reviewer's original (non-transcribed) comment.
- Fig2D → supplement: author-approved; figure-panel move is a layout task at figure-assembly time.
- Red fold wording (SuppFig15): quote T2→T3 ≈17× / floor-aware, not ~30× T1→T3 (05 log).

## Round-2 edits (after author sign-off on remaining judgments)
6. **[Red fold, SuppFig15] floor-aware wording applied (L30)** — "increased ~30-fold from T1 to T3,
   beginning at T2" → "rose from near the detection limit at T1 to a marked peak at T3, increasing
   ~17-fold between the detectable T2 and T3 stages (… a T1-referenced fold is not quoted because
   T1 signal is at the detection floor)." Uses the only stable figure-supported ratio (T3/T2≈17×);
   T1 sits at the detection floor (0.002–0.009 mg/L) so T1-denominated folds are unstable.
   Welch ANOVA statistic unchanged.

## Agent judgment (delegated by author): two-layer "on/off gate" vocabulary (L233) — KEEP, no change
Decision: retain as written. Rationale:
- Distinct from the retracted content: retracted model was a *four-layer* "Gatekeeper switch" brand;
  L233 is a *two-layer* "on/off gate" (different structure; the retracted brand name is not used).
  "Gate" is generic accessibility-state vocabulary, not the retracted term.
- Properly hedged and Discussion-only ("consistent with / associated with / appear to translate");
  the binary "on/off gate" maps onto the paper's defensible binary Shielded/Exposed classification.
- Confirmed the reframe did NOT promote "gate"/"two-layer" into Abstract or Intro (grep: only L233).
Guardrail enforced: vocabulary stays contained at L233 where it is carefully qualified.

## Mutilka dataset (C16) — UNRESOLVED, closed as such
Author also does not recognise the term; no external dataset matches and Fig 7 is built from the
study's own gene sets. Recorded as "[UNRESOLVED — reviewer's original comment needed]"; the
point-by-point response will state that Fig 7 uses the study's own methylation-proximal gene sets
(n = 3,197 / 459 / 448) vs the M145 KEGG annotation, i.e. not an external ("Mutilka") dataset, and
ask the reviewer to clarify if a specific external dataset was intended.

## JP manuscript sync (Manuscript_JP_NAR_numbered_inline.md) — all 6 edits mirrored
The same six edits were applied to the Japanese manuscript, matching the EN wording and register:
1. Limitations header → "Internal specificity control for Nanopore 4mC calling in high-GC contexts."
2. Intro protection-zone width → "全遺伝子のTSSメタジーンではTSSの約±300 bp以内で枯渇" (was 約2,200 bp).
3. Intro "主要な知見は3点である。第一/第二/第三" → "本研究の中心的知見は…" 1中心+2補助的観察.
4. Abstract opening reframed (防御関連methylationが制御ゲノムに対して何をなすか→空間的に組織化;
   relational framing) with 4mC + dual mark demoted to 2 enabling observations.
5. Abstract co-dynamics demoted: ρ=0.43 clause removed; robust static co-localisation
   (Compartment-A odds ratio = 1.96, p = 3.5 × 10⁻¹⁷); temporal co-variation "示唆されるにとどまり
   確立されていない…Discussion参照".
6. Red fold floor-aware: "T1では検出限界付近…検出可能なT2からT3にかけて約17倍" (was T1→T3 約30倍).
Verified: all 5 old JP phrasings removed (count 0); 5 of 6 new-phrasing probe strings matched
directly. The 6th probe (the Abstract central-clause) returned 0 only because the probe string used
"制御ゲノム" whereas the applied text reads "当該ゲノムに対して転写を指令するのではなく空間的に組織化されている"
(printed and read verbatim) — the edit itself is correct; the probe wording, not the manuscript, was
the mismatch. Abstract & Intro read cleanly.
Note: JP ρ=0.43 full treatment remains in Discussion (L212) unchanged, mirroring EN L210.

## Figure-legend + supplementary-list edits (C7 notation, C10/C12/C13) — APPLIED
Two prose files edited in place (git-reversible):
- **Figure-Legends.md**: Fig 2 legend rewritten for the 3-panel main figure (A/B/C); added a
  panel (b)/(c) complementarity note (C13: metagene density vs per-gene nearest-distance, not a
  contradiction); made the C12 point explicit (293 bp = operating point on a continuous unimodal
  distribution, not two discrete classes); unified notation m4C→4mC (9 occurrences, Fig 2–3+);
  corrected Fig 2(b) n to 2,646; recorded the former panel (d) demotion to the negative-control
  Supplementary figure with the main-text null retained.
- **Supplementary-Materials-List.md**: unified m4C→4mC (9 occurrences); added
  **Supplementary Figure 19** legend for the demoted negative control (former Fig 2d), citing the
  n.s. result (p=0.55 / p=0.12) and the retained main-text falsification statistics.

## Pre-submission consistency audit (nar-revision s2 + direct scan) — CLEAN
- nar-revision s2 on build sources: 0 CRITICAL violations; canonical 62/989/1,051/293 bp all
  present; 17 informational hits only in non-build working drafts (Key-Claims, Reviewer-Response-
  Draft, Story-Draft — archive candidates).
- Direct retracted-term scan of the edited files (_inline EN+JP, Figure-Legends,
  Supplementary-Materials-List): every hit for AUC / "Gatekeeper switch" / "four-layer" /
  "57 Exposed" is a negation, withdrawal, or "retracted…replaced" reference — none is a live
  claim. The retracted 0.908 AUC value appears nowhere. See 13_consistency_audit_card.md.
