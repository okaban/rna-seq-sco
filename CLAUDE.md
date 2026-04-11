# Project: S. coelicolor A3(2) M145 RNA-seq + Epigenome Integration

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
