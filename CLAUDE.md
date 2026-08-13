# CLAUDE.md, project rules for SpectralBiodiversity

<!--
  TEMPLATE USAGE
  This is the lab's CLAUDE.md template: a lean core that every governed repo
  keeps, followed by opt-in modules that a project switches on only if it uses
  them.

  1. Fill every <FILL: ...> placeholder. Where a value is genuinely unknown at
     init, leave the placeholder and note it under "Known pending items" rather
     than inventing one.
  2. Delete any opt-in module you are not using. Do not leave an empty module in
     place; an unused module is noise.
  3. Do not restate rules that live in SCIENCE_PRINCIPLES.md or its domain files.
     This file adds only project-specific operational content and points to the
     principles for the rest.
-->

## Lab principles source

This project inherits, and does not override, the lab principles at the versions
recorded here. These files are the governing rules. This CLAUDE.md adds only
project-specific operational detail and never restates them.

- SCIENCE_PRINCIPLES.md, version <FILL: e.g. 2.0>, at commit 1caa354
- <FILL, if a pipeline project: SCIENCE_PRINCIPLES_PIPELINES.md, same version and commit>
- <FILL, if a text-analysis project: SCIENCE_PRINCIPLES_TEXT_ANALYSIS.md, same version and commit>

Record the commit hash as it was at adoption. A published result must be
traceable to the exact version of the standards that governed it.

## Project context and stakes

<FILL: one short paragraph. What this project produces, who uses the output, and
what the cost of an error is. State the asymmetry if there is one, for example
whether a false negative or a false positive is the more harmful failure, since
that sets the conservative default.>

- PI: <FILL>
- Collaborators: <FILL>
- Repository: <FILL: URL>

## Language, packages, and environment

<FILL or keep the lab defaults below and edit as needed.>

- R only in pipeline scripts. No Python or shell in the pipeline itself.
- CRAN packages only, declared in renv.lock via renv. No GitHub-only or
  Bioconductor packages without explicit discussion.
- Preferred packages: <FILL: e.g. dplyr, purrr, tidyr, stringr, ...>
- Tidyverse style: native pipe (|>), snake_case, no library() calls inside
  functions.

Environment variables (credentials come from the environment only, never from
files in the repo):

| Variable | Purpose | Set where |
|----------|---------|-----------|
| <FILL>   | <FILL>  | <FILL>    |

## Paths and execution

- Relative paths only. Every script is independently runnable from the project
  root (for example `Rscript R/0X_scriptname.R`).
- Scripts communicate only through files in the project's data and output
  directories, never through a shared R session or global state.
- Every script checks for its required inputs at the top and fails with a clear
  stop() message if any are missing.

## Hard rules for this project

- Permitted data sources: <FILL>
- Forbidden data sources: <FILL>
- Tracked directories: <FILL: e.g. data snapshots/manifests, human override files>
- Gitignored directories: <FILL: e.g. raw data, processed data, outputs, figures, run logs>
- <FILL: any project-specific hard constraint not covered above>

(The universal constraints, no hardcoded paths or credentials, never commit data
or credentials or large binaries, live in SCIENCE_PRINCIPLES.md and are not
restated here.)

## Session journal (SESSION_LOG.md)

This project keeps a human-readable session journal at SESSION_LOG.md in the repo
root. It is the committed, durable record of what each session accomplished, and
it is distinct from the on-disk run log (see SCIENCE_PRINCIPLES.md for the run
log).

Convention:

- Every substantive task appends a dated entry, including read-and-report tasks.
- Entries are in reverse chronological order, newest at the top.
- Each entry header is `## YYYY-MM-DD HH:MM UTC — <short title>`. Use UTC, and
  include the time, not just the date, because there are often several sessions
  in a day.
- Write the entry as work completes, not only at the very end, so a dropped
  connection or a usage cap loses at most the last increment.
- Commit and push the entry before the session ends or hands back. An entry that
  exists only in the working tree, or only locally, is treated as lost. The
  committed, pushed SESSION_LOG is the durable channel out of the tool; terminal
  scrollback is not recoverable.
- Log Claude's structured outputs (reports, audits, investigation summaries),
  not the prompts or the back-and-forth.

## Output metadata and reproducibility

<FILL for a pipeline project: the output-provenance format is defined in
SCIENCE_PRINCIPLES_PIPELINES.md (companion .meta.json / manifest with run
datetime UTC, code commit hash, input sources, session info). Record here only
what is project-specific, for example the exact index URL scraped or the
manifest filename.>

## Data use and citation

<FILL: the data-use terms and citation obligations for this project's sources.>

## Known pending items

<FILL: placeholders left unresolved at init, and anything a future session should
pick up. Delete once resolved.>

---

<!--
  OPT-IN MODULES
  Keep only the modules this project uses; delete the rest, including this
  comment. Each module carries only this project's parameters and points to the
  principles file for the standing rule behind it.
-->

## MODULE, manual overrides (delete if unused)

Human overrides for this project are stored in <FILL: e.g. data/manual_overrides.csv>.

- Required columns: <FILL: e.g. party_name, flag, reviewer_name, review_date, notes>
- Overrides are marked visibly in the review output with <FILL: e.g. an orange
  badge showing reviewer name and date>.

(The rule that this file is human-only and is never modified or silently
overwritten by the pipeline lives in SCIENCE_PRINCIPLES.md.)

## MODULE, two-layer reproducibility (delete if unused)

This project maintains two analysis layers with deliberately different
reproducibility guarantees. Every output is labeled with the layer that produced
it, and the two are never given equal evidential weight.

- Layer 1, deterministic (policy-grade): <FILL: which scripts or steps. renv-locked
  R, no model calls, reproducible bit-for-bit from renv.lock.>
- Layer 2, LLM extraction (research-grade): <FILL: which steps. Reproducibility
  rests on a pinned model version string (never a floating alias), temperature 0,
  the prompt versioned as code, and responses cached by hash of (model + prompt +
  input). The response cache is a preserved artifact, committed or archived, never
  throwaway.>

A Layer 2 finding is never given the evidential weight of a Layer 1 result.

## MODULE, table extraction and QA (delete if unused)

Standing rules for every table extraction, re-extraction, or correction pass;
they hold without being restated in a task prompt.

- Provenance: model-read table values are Layer 2; record the model version
  string and the prompt sha256 per row. Deterministic parsing and assembly are
  Layer 1.
- Verbatim capture: transcribe each cell exactly as printed. Never split a
  combined cell, fill a missing value from a bound, normalize a printed
  inconsistency, or tidy footnote markers and garbled tokens away. Parsing is a
  separate deterministic pass against the captured strings, never against the
  image.
- Reading structure: locate a table by its caption confirmed verbatim in the
  source, not by a banked page number. On a shared page, read only the table
  matching the target caption. Capture the full stratification key, each label
  level in its own named column.
- Stop-clean and no-guess: a table that cannot be resolved is recorded failed
  with a reason and the pass continues; one bad table never aborts the batch and
  is never silently dropped. Flag an ambiguous span for a human rather than
  guessing.
- Write safety: passes are additive or surgical (new files, or field-level edits
  to named rows), never a read/write round-trip over a file with historical mixed
  quoting. Assert untouched records byte-identical to a pre-write backup. A
  surrogate row_id is assigned once and never renumbered.
- Metadata: derive table metadata (page, span, key) from the Layer 1 manifest,
  never hand-typed into a dispatch.
- Human sign-off: correcting a banked policy-grade value is two-phase, a
  read-only diagnosis and proposal that a human verifies against the source, then
  the write. Never single-pass correct-and-bank.

## MODULE, standards drift audit (delete if unused)

Standing rule: about every week of active development, run a standards check and
record the result as a dated SESSION_LOG.md entry. Check that renv.lock matches
the installed library, that every model version string is explicitly pinned,
that each cache key's prompt matches the current versioned prompt, that the
gold-set eval has been rerun with precision and recall compared to the prior run,
and that every script added since the last audit is labeled Layer 1 or Layer 2.

(This is the standing rule only. A dedicated drift-audit skill may own the
procedure later.)
