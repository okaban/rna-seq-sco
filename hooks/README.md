# Git hooks (versioned)

The sandbox the analysis agent runs in cannot write `.git/config`, so the hooks
are installed by the author once per clone:

    git config core.hooksPath hooks

`pre-commit` runs the supersession tripwire (make_analysis_index.py --check) whenever
anything under 11_epigenome_integration/analysis is staged, refuses stale/backup files, and blocks
the commit on FAIL. Bypass consciously with `git commit --no-verify` and record
why in `open_items.md`.
