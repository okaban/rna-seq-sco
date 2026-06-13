<!-- Inherits global rules from ~/.claude/CLAUDE.md -->

# Project: epi-trans — S. coelicolor A3(2) M145 Epigenome × Transcriptome Integration

<!-- 表示名 epi-trans（旧称 rna-seq）。ディレクトリ名・project_id・conda env は配管として rna-seq のまま。 -->


## ⚡ CANONICAL LOCATIONS — read this first every session

- **Single source of truth = the dashboard:** `/Users/okaban/obsidian/Research/rna-seq/00-Hub.md`, section "🟢 LATEST". At the start of any session, read it to find the current outline / manuscript / figures. Do NOT guess by browsing folders.
- This repo `/Users/okaban/bioinfo/rna-seq/` holds **analysis, figures, manuscript PDFs** (absolute paths are hardcoded → never move/rename top-level dirs).
- **Outline, manuscript source, notes live in Obsidian:** `/Users/okaban/obsidian/Research/rna-seq/Writing/`. Current outline: `Writing/Paper-Outline-v3.md`.
- Methylation raw data: `/Users/okaban/bioinfo/methyl/260102_M145/`.
- Full cross-repo map: `/Users/okaban/obsidian/Research/rna-seq/Repository-Map.md`.
- **When you create a new version** of the outline/manuscript: (1) update the "🟢 LATEST" table in 00-Hub.md, (2) add a deprecation banner to the old file. Keep dated/versioned filenames for history.

## Operational Guidelines

- Execute analysis scripts (Python, R, bash) and file operations (mv, cp, rm, mkdir, ln) without asking for confirmation
- Read/write files freely within the project directory `/Users/okaban/bioinfo/rna-seq/`
- Web searches and fetches are pre-authorized for literature and database lookups
- When running analysis pipelines, proceed through all steps automatically unless errors occur

## Report Management Rules

- All analysis reports (`.md`) MUST be saved in `/Users/okaban/bioinfo/rna-seq/reports/`
- Do NOT save reports directly in analysis subdirectories
- File naming convention: `YYMMDD_<description>_report.md` (snake_case, common abbreviations OK: DEGs, GO, KEGG, BGC, TF, DMG, TSS, etc.)
- When a report is relevant to a specific analysis directory, create a symlink from the analysis directory to `/reports/`
- When updating an existing report, update the file in `/reports/` (symlinks will automatically reflect the change)
- Update `reports/README.md` when adding new reports (add entry to the appropriate section in the table)

# Project Overview

Research assistant for a peer-reviewed epigenomics paper on *Streptomyces coelicolor* M145. Paper targets Nature Microbiology, EMBO Journal, or Nucleic Acids Research.

Current paper status: story confirmed (Story Option B), manuscript ~12,900 words across 5 sections. Key figures: Figure 1–5 (main), Supplementary S1–S16. Remaining work: Fig 2 co-modification figure generation, Abstract language revision, Fig 5 model description.

Key claims of the paper:
- Genome-wide methylation–expression correlation is spurious (geographic confounding by core/arm chromosome structure)
- Two independent methylation systems: m4C at GCCGGC (RM) and 6mA at AAGCCCG (solo MTase SC_RS17645/SCO3104)
- TSS protection zone (293 bp threshold, AUC=0.908, 95% CI 0.886–0.929; old 0.917 was the stale n=62 value) structurally partitions 998 Shielded vs 57 Exposed regulatory genes
- 57 Exposed TFs show synchronized developmental transition with two-block z-score structure (ρ=−0.995)
- Per-read co-modification: OR=921×, n=244/254 pairs

# Project-Specific Standards

- Paper structure: Abstract / Introduction / Results / Discussion / Methods / References
- Gene names: italicized (*SCO3104*, *scbA*); protein names: roman (SCO3104)
- Critique prior art alignment: claims must be consistent with published *Streptomyces* epigenomics literature specifically
- Journal fit: framing must satisfy a Nature Microbiology reviewer skeptical of epigenomics–transcription linkages
