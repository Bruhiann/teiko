# Loblaw Bio — Immune Cell Count Analysis

Tools to load, analyze, and explore immune cell population counts from a
clinical trial of the drug candidate **miraclib**, built for Bob Loblaw.

The project loads `cell-count.csv` into a normalized SQLite database (Part 1),
computes per-sample population relative frequencies (Part 2), compares
responders vs non-responders with statistics and a boxplot (Part 3), and
summarizes a baseline subset (Part 4) — all surfaced in an interactive Dash
dashboard.

## Quick start (GitHub Codespaces)

The project is graded via `make` targets and works out of the box in a
Codespace (or any Linux/macOS environment with Python 3.9+ and `make`).

```bash
make setup        # install dependencies from requirements.txt
make pipeline     # build the DB and generate all tables + plots (Parts 1-4)
make dashboard    # launch the interactive dashboard
```

| Target | What it does |
|--------|--------------|
| `make setup` | Installs all dependencies (`pip install -r requirements.txt`). |
| `make pipeline` | Runs the whole pipeline with no manual steps: `load_data.py` initializes/loads `cell-count.db` (Part 1), then `analysis.py` prints the summary tables and writes `cell_frequencies.csv` and `responder_boxplot.png` (Parts 2-4). |
| `make dashboard` | Starts the Dash server (binds `0.0.0.0:8050`). |

`make` invokes `python3` by default; override with `make pipeline PYTHON=python`
if your environment uses a different interpreter.

### Viewing the dashboard

> **Note:** The dashboard runs locally at **http://localhost:8050** via
> `make dashboard`. In GitHub Codespaces it's accessible via the forwarded
> port in the **Ports** tab.

Run `make dashboard`, then open **http://localhost:8050**.

In a Codespace, the server binds to `0.0.0.0:8050`; Codespaces auto-forwards
the port — open it from the **Ports** tab (or the toast notification) to load
the dashboard in your browser. The dashboard has three tabs:

- **Overview (Part 2):** server-side paginated frequency table with sample-id
  and population filters.
- **Responders (Part 3):** configurable filters (condition / treatment /
  sample type / timepoint, including single-timepoint views), boxplot, and the
  per-population statistics table.
- **Baseline subset (Part 4):** configurable filters with project / response /
  sex breakdown charts.

## Reproducing the outputs manually

`make pipeline` is equivalent to:

```bash
python load_data.py     # -> cell-count.db
python analysis.py      # -> prints Part 2-4 tables; writes outputs below
```

Generated outputs (the CSV and plot are committed for convenience; the `.db`
is rebuilt by the pipeline):

- `cell-count.db` — the SQLite database (binary, git-ignored).
- `cell_frequencies.csv` — the full Part 2 frequency table (52,500 rows).
- `responder_boxplot.png` — the Part 3 responders-vs-non-responders boxplot.

## Database schema

The flat CSV is normalized into **three tables** (see `load_data.py`):

| Table | Grain | Columns |
|-------|-------|---------|
| `subjects` | one row per patient | `subject` (PK), `project`, `condition`, `age`, `sex`, `treatment`, `response` |
| `samples` | one row per sample | `sample` (PK), `subject` (FK → subjects), `sample_type`, `time_from_treatment_start` |
| `cell_counts` | one row per (sample, population) | `sample` (FK → samples), `population`, `count`, PK `(sample, population)` |

### Rationale

- **Normalization removes redundancy.** Attributes that are constant across a
  subject's samples (`project`, `condition`, `age`, `sex`, `treatment`,
  `response`) live once in `subjects` instead of being repeated on every CSV
  row. This was verified against the data: none of these columns vary within a
  subject. It prevents update anomalies (e.g. a subject's response recorded
  inconsistently across rows).
- **`cell_counts` is stored long/tidy**, one row per population rather than
  five fixed columns. This is the key design choice: relative-frequency math
  (`SUM(count) … OVER (PARTITION BY sample)`) and per-population grouping/tests
  become generic SQL that never hardcodes the five population names. Adding a
  sixth population is a data change, not a schema migration.
- **Foreign keys + indexes** (`samples.subject`, `subjects.condition`,
  `subjects.treatment`, `cell_counts.population`) keep the typical
  filter/group-by queries fast and the relationships explicit.

### Scaling to hundreds of projects / thousands of samples

The schema is designed to grow along the natural dimensions of the data:

- **Volume.** Hundreds of projects and thousands of samples is still small for
  SQLite, but the structure ports directly to a server database (Postgres) with
  no model changes. The long `cell_counts` table grows as
  `samples × populations`; B-tree indexes on the filter columns plus a covering
  index on `cell_counts(population, sample)` keep aggregate queries scalable. If
  the per-sample total is needed constantly, it can be materialized (a
  `samples.total_count` column or a small `sample_totals` table) so percentages
  don't recompute the window each time.
- **New analytics / new metadata.** Because subject- and sample-level metadata
  are separated, new attributes (e.g. dose, batch, assay panel, study arm)
  attach to the right table without touching cell data. The long-format counts
  mean a new analysis is just a new query/filter, not a reshape.
- **Many populations / multi-omics.** The tidy `(entity, feature, value)` shape
  generalizes: more populations, or entirely new measurement types, are new
  rows (or sibling tables) rather than ever-widening columns — avoiding sparse,
  hard-to-evolve wide tables.
- **Partitioning by project.** A `project` dimension is already first-class, so
  data can be partitioned/sharded or row-level access-controlled per project as
  the trial portfolio grows.

## Code structure

```
load_data.py   Part 1: schema definition + CSV loader. No arguments; run
               directly to (re)build cell-count.db. Idempotent.
analysis.py    Analytics layer. Pure functions (one per question) that take a
               DB connection and return pandas DataFrames. Runnable as a CLI
               that prints all Part 2-4 tables and writes the CSV + plot.
app.py         Presentation layer. A Dash dashboard that imports analysis.py
               and renders its DataFrames as interactive tables and figures.
Makefile       setup / pipeline / dashboard orchestration.
requirements.txt  Pinned-by-minimum dependencies.
cell-count.csv    Source data.
```

### Why this design

- **Separation of concerns / single source of truth.** All analytical logic
  lives in `analysis.py` as small, importable, side-effect-free functions
  (`cell_frequencies`, `cell_frequencies_page`, `responder_frequencies`,
  `responder_stats`, `responder_boxplot`, `baseline_subset`,
  `subset_breakdowns`). The CLI (`make pipeline`) and the dashboard call the
  **same** functions, so the static deliverables and the interactive app can
  never drift out of sync.
- **SQL does the heavy lifting.** Filtering, per-sample totals (window
  functions), and pagination (`LIMIT/OFFSET`) run in the database, so the app
  never loads the full 52,500-row table into memory — the Overview tab is
  server-side paginated.
- **Reproducibility.** Each step is a plain, no-argument command, wired
  together by the Makefile so the whole pipeline runs start to finish with no
  manual intervention.

## Analyses summary

- **Part 2 — Overview.** `cell_frequencies()` returns one row per
  (sample, population) with `total_count`, `count`, and `percentage` (the
  population's share of the sample's total cells). Verified: every sample's
  percentages sum to 100.
- **Part 3 — Responders vs non-responders.** On the melanoma / miraclib / PBMC
  subset, compares relative frequencies per population with a two-sided
  Mann-Whitney U test (non-parametric) and Benjamini-Hochberg FDR correction.
  *Result: no population is significant after correction; `cd4_t_cell` is
  closest (raw p ≈ 0.013, adj p ≈ 0.067).*
- **Part 4 — Baseline subset.** `baseline_subset()` filters to baseline (day 0)
  melanoma PBMC miraclib samples (656 samples) and breaks them down by project
  (prj1=384, prj3=272), response (yes=331, no=325), and sex (M=344, F=312).

## Requirements

Python 3.9+ and the packages in `requirements.txt` (pandas, scipy, matplotlib,
seaborn, dash, plotly). Install with `make setup`.
