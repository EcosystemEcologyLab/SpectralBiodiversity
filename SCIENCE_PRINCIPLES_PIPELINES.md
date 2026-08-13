# SCIENCE_PRINCIPLES_PIPELINES.md — Data Pipeline Projects

**Version:** 1.0
**Repository:** EcosystemEcologyLab/lab-principles
**Applies to:** EcosystemEcologyLab data pipeline projects (R or Python)
**Extends:** SCIENCE_PRINCIPLES.md — read that file first

When copying this file into a project repository, record the source commit hash
in the project's CLAUDE.md alongside the hash for SCIENCE_PRINCIPLES.md.

---

## Scope

These rules apply to any project whose primary purpose is acquiring, processing,
quality-controlling, or analysing scientific datasets. They apply regardless of
whether the pipeline is written in R, Python, or a combination.

---

## Dependency management

- All R dependencies must be declared and version-locked using `renv`
- All Python dependencies must be declared and version-locked using
  `pyproject.toml` or `requirements.txt` with pinned versions
- Never install packages interactively without immediately updating the
  lockfile and committing it
- Never introduce a new dependency without flagging it for discussion first
- The lockfile is a scientific record — treat it as such

---

## Output metadata

Every output file (CSV, RDS, NetCDF, figure, HTML report) must carry provenance
metadata. The metadata may be embedded in the file (e.g. as a header block or
file attributes) or stored in a companion file with the same base name and a
`.meta.json` extension.

Required fields:

| Field | Content | Example |
|-------|---------|---------|
| `run_datetime_utc` | ISO 8601 timestamp of pipeline run | `2026-03-27T18:25:56Z` |
| `pipeline_version` | Git commit hash of the repository at run time | `34befe2` |
| `input_sources` | URLs, DOIs, or file paths of all primary inputs | snapshot CSV path + per-site DOIs |
| `r_session_info` | Saved output of `sessionInfo()` | saved to `outputs/session_info.txt` |
| `notes` | Any manual decisions, overrides, or deviations from defaults | free text |

Implementation rules:
- `pipeline_version` must be captured programmatically at run time:
  `system("git rev-parse --short HEAD", intern = TRUE)` in R
- `r_session_info` must be written to `outputs/session_info.txt` at the end
  of every pipeline run — this file is gitignored but must be present with
  every output set
- `notes` is required even if empty — an empty string is acceptable; a
  missing field is not

---

## Exclusion logging

Every record excluded from analysis must be logged. Exclusions and unknowns
are distinct categories and must be logged separately.

### Exclusion log
- File: `outputs/exclusion_log.csv` (gitignored — regenerated each run)
- Required columns:

| Column | Content |
|--------|---------|
| `site_id` | FLUXNET site ID or equivalent record identifier |
| `variable` | Variable name or `ALL` if the whole record is excluded |
| `timestamp` | Record timestamp or `ALL` if the whole site-year is excluded |
| `reason` | Human-readable reason for exclusion |
| `threshold` | The threshold or rule applied (e.g. `QC_THRESHOLD_YY=0.75`) |
| `excluded_by` | Script name that performed the exclusion |

### Unknown log
- File: `outputs/unknown_log.csv` (gitignored — regenerated each run)
- Required columns: `record_id`, `reason`, `logged_by`
- A record is UNKNOWN when it cannot be assessed — not when it fails QC.
  Failed QC → exclusion log. Missing data → unknown log.

### Rules
- Both logs must be written even if empty (zero-row CSV with headers)
- Summary counts from both logs must be printed to the console at the end
  of each QC script run
- No record may be silently dropped — every exclusion must appear in one
  of the two logs

---

## Human override files

When a scientist makes a manual decision that overrides a pipeline default
(e.g. manually including or excluding a site, overriding a QC threshold for
a specific site-year), that decision must be stored in a human override file.

- File location: `data/overrides/` (git-tracked — these are scientific decisions)
- File format: CSV with columns `record_id`, `decision`, `reason`, `date`, `author`
- Claude must never modify override files — read them, apply them, flag them
  in output metadata, but never write to them
- Override files must survive pipeline reruns — they are inputs, not outputs

---

## Script conventions

- Every script must begin with `source("R/pipeline_config.R")` and
  `check_pipeline_config()` (or the Python equivalent)
- Scripts are numbered and must be run in order — document dependencies
  between scripts explicitly
- Scripts communicate only via files in `data/` or `outputs/` — never via
  R global environment variables or Python module-level state
- Every script must write a completion message to the console on success:
  `message("Script XX complete: N records processed, M excluded, K unknown")`

---

## Configuration and thresholds

- All scientifically meaningful thresholds (QC cutoffs, minimum data coverage
  requirements, aggregation rules) must be declared as named constants in a
  configuration file — never as magic numbers inline in scripts
- Named constants must be documented with their units, valid range, default
  value, and the scientific rationale for the default
- Changing a threshold is a scientific decision — it must be committed with
  a clear commit message explaining why

---

## Data directory conventions

| Directory | Git-tracked | Content |
|-----------|-------------|---------|
| `data/snapshots/` | Yes | Timestamped manifests of input data |
| `data/overrides/` | Yes | Human override files |
| `data/raw/` | No | Downloaded source data |
| `data/extracted/` | No | Unzipped/parsed source data |
| `data/processed/` | No | Pipeline outputs |
| `outputs/` | No | Final analysis outputs and logs |
| `figures/` | No | Generated figures |
