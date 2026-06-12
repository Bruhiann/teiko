# Loblaw Bio — Immune Cell Count Analysis

Tools to load, analyze, and explore immune cell population counts from a
clinical trial of the drug candidate **miraclib**, built for Bob Loblaw.

## Contents

| File | Purpose |
|------|---------|
| `load_data.py` | Initialize the SQLite database and load `cell-count.csv` (Part 1). |
| `analysis.py`  | Reusable analytics functions for Parts 2–4; runnable as a CLI report. |
| `app.py`       | Interactive Dash dashboard surfacing all analyses. |
| `cell-count.csv` | Source data: one row per biological sample. |

## Setup

```bash
pip install -r requirements.txt
```

## Usage

1. **Build the database** (creates `cell-count.db` in the repo root):

   ```bash
   python load_data.py
   ```

2. **Run the analysis report** (prints tables and writes
   `cell_frequencies.csv` + `responder_boxplot.png`):

   ```bash
   python analysis.py
   ```

3. **Launch the dashboard** (then open http://127.0.0.1:8050):

   ```bash
   python app.py
   ```

## Database schema

Normalized into three tables:

- **`subjects`** — one row per patient: `subject` (PK), `project`,
  `condition`, `age`, `sex`, `treatment`, `response`. These attributes are
  constant across a subject's samples.
- **`samples`** — one row per sample: `sample` (PK), `subject` (FK),
  `sample_type`, `time_from_treatment_start`.
- **`cell_counts`** — long/tidy format: `(sample, population)` PK with raw
  `count`. Storing populations as rows (not five fixed columns) lets the
  analysis compute relative frequencies generically.

## Analyses

- **Part 2 — Overview.** `cell_frequencies()` returns one row per
  (sample, population) with `total_count`, `count`, and `percentage`
  (the population's share of the sample's total cells).
- **Part 3 — Responders vs non-responders.** For a chosen subset (default
  melanoma / miraclib / PBMC), compares relative frequencies between
  responders and non-responders per population using a two-sided
  Mann-Whitney U test with Benjamini-Hochberg FDR correction, plus a
  boxplot. *Result on the default subset: no population is significant after
  correction (cd4_t_cell is closest, adj p ≈ 0.067).*
- **Part 4 — Baseline subset.** `baseline_subset()` filters to baseline
  (day 0) melanoma PBMC miraclib samples and breaks them down by project,
  response, and sex.

## Dashboard tabs

- **Overview (Part 2):** searchable/sortable frequency summary table.
- **Responders (Part 3):** configurable filters (condition, treatment,
  sample type, timepoint — including a single-timepoint view), boxplot, and
  the per-population statistics table.
- **Baseline subset (Part 4):** configurable filters with project /
  response / sex breakdown charts.
