#!/usr/bin/env bash
# 12 · Commit preparation for the NAR revision (260716)
# ⚠️ REVIEW BEFORE RUNNING. Nothing here is executed by the agent — git auth is author-side.
# Two separate repositories are involved. Run each block from its own repo root.
# Both repos have MANY unrelated pre-existing modified files (regenerated PDFs/SVGs/TSVs);
# these scripts stage ONLY the revision-relevant paths so the commit stays reviewable.
set -euo pipefail

########################################################################
# REPO 1 — figure-generator scripts   /Users/okaban/bioinfo/rna-seq  (HEAD was 49a70dd)
########################################################################
cd /Users/okaban/bioinfo/rna-seq

git add -p 15_paper_figures/scripts/01_figure1_landscape.py \
          15_paper_figures/scripts/02b_figure2_RM_redistribution.py \
          11_epigenome_integration/analysis/77_reviewer_figures/fig_modpos_and_redistribution.py
# (use -p to review each hunk; or drop -p to stage whole files)

git commit -m "Revise figures for NAR reviewer comments (C3-C14)

- Fig1 (01_figure1_landscape.py): state core/arm boundary on panel A;
  recolour oriC green->grey (#333333); re-anchor Fig1C significance
  markers to a per-region ceiling to remove asterisk/number overlap;
  set Fig1D motif-logo titles to black.
- Fig2 (02b_figure2_RM_redistribution.py): correct the promoter
  protection-zone band from ~2,200 bp to the data-supported +/-300 bp;
  demote the negative-control panel (former 2D) to a standalone
  Supplementary figure; main Fig2 re-laid-out as 3 panels (A/B/C).
- SuppFig10 (fig_modpos_and_redistribution.py): unify 'm4C'->'4mC';
  replace detached red text box with an integrated summary panel;
  correct AAGCCCG 6mA site count 447->441 (canonical, 22.80% occupancy);
  move in-panel grey prose to the caption; add GCCGGC breakdown bar."

########################################################################
# REPO 2 — manuscript + figure slots   /Users/okaban/obsidian/Research/rna-seq  (HEAD was a35d712)
########################################################################
cd /Users/okaban/obsidian/Research/rna-seq

# NOTE: the edited body file Writing/Manuscript_EN_NAR_numbered_inline.md is UNTRACKED
#       (git status '??'). The tracked file is Manuscript_EN_NAR_numbered.md.
#       Decide whether the _inline file should be committed as the working copy,
#       or whether _numbered.md should be regenerated from it (build_nar_numbered.py).
#       If _inline is the working manuscript, add it explicitly:
git add Writing/Manuscript_EN_NAR_numbered_inline.md
git add Writing/fig_images/Figure1.png \
        Writing/fig_images/Figure2.png \
        Writing/fig_images/SuppFigure10.png \
        Writing/fig_images/SuppFigure15.png \
        Writing/fig_images/SuppFigure_neg_control.png

git commit -m "Apply NAR revision body-text edits + revised figure slots

Body text (Manuscript_EN_NAR_numbered_inline.md):
- Abstract/Intro reframed to one central finding + two enabling
  observations; relational verb ('spatially organised relative to')
  in place of agentive 'organises' (reverse-causation hedge).
- Demote rho=0.43 temporal co-dynamics (n~3) out of the Abstract;
  keep robust static co-localisation (Compartment-A OR=1.96, p=3.5e-17).
- Correct protection-zone width ~2,200 bp -> +/-300 bp (Intro).
- Rename Limitations header 'External validation...' ->
  'Internal specificity control for Nanopore 4mC calling...'.
- Undecylprodigiosin (Red) fold: floor-aware wording (T2->T3 ~17x)
  in place of '~30-fold T1 to T3' (T1 at detection floor).

Figures: revised Figure1/2, SuppFigure10, SuppFigure15 (merged growth),
and the new standalone SuppFigure_neg_control (demoted former Fig2D)."

echo "Both commits staged. Push when ready (author credentials):"
echo "  (repo1) cd /Users/okaban/bioinfo/rna-seq && git push"
echo "  (repo2) cd /Users/okaban/obsidian/Research/rna-seq && git push"
